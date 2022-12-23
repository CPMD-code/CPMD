MODULE epr_hyp_utils
  USE adat,                            ONLY: covrad,&
                                             elem
  USE atwf,                            ONLY: atrg_epr,&
                                             atwf_mod,&
                                             atwfr_epr,&
                                             atwr
  USE cnst,                            ONLY: fbohr,&
                                             fpi,&
                                             pi
  USE cppt,                            ONLY: gk,&
                                             hg,&
                                             inyh,&
                                             nzh
  USE error_handling,                  ONLY: stopgm
  USE fftmain_utils,                   ONLY: fwfftn
  USE fitpack_utils,                   ONLY: curv1,&
                                             curvi
  USE geq0mod,                         ONLY: geq0
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum,&
                                             mp_sync
  USE parac,                           ONLY: parai,&
                                             paral
  USE pslo,                            ONLY: pslo_com
  USE sfac,                            ONLY: ei1,&
                                             ei2,&
                                             ei3
  USE system,                          ONLY: fpar,&
                                             ncpw,&
                                             parm
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: epr_hyp

CONTAINS

  ! ==================================================================
  SUBROUTINE epr_hyp(c0,rhoes,psi,nstate)
    ! ==--------------------------------------------------------------==
    ! == Computes hyperfine parameters (VERY rough method!!!)         ==
    ! == with the method of C.G. Van de Walle and P.E. Blochl         ==
    ! ==                    Phys. Rev. B 47, 4244 (1993)              ==
    ! == how to use this routine                                      ==
    ! == 1) use TM PSP's                                              ==
    ! == 2) force CPMD to read in the atomic pseudo wavefunctions     ==
    ! ==    by putting the following in your input file:              ==
    ! ==         &BASIS                                               ==
    ! ==             PSEUDO AO [NUM]                                  ==
    ! ==             01 2 ... [NUM-1]                                ==
    ! ==             ... ! and this for all elements                  ==
    ! ==         &END                                                 ==
    ! == 3) for every element, place a [elementsname].ae in the       ==
    ! ==    same directory where your input file is located           ==
    ! ==    These .ae files must be in the same format as the         ==
    ! ==    USPP pseudopotential generation program of D. Vanderbilt  ==
    ! == 4) run the EPR routine with the option HYP                   ==
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    ! Arguments & common variables
    COMPLEX(real_8)                          :: c0(ncpw%ngw,*)
    REAL(real_8)                             :: rhoes(fpar%nnr1)
    COMPLEX(real_8)                          :: psi(:)
    INTEGER                                  :: nstate

    CHARACTER(*), PARAMETER                  :: procedureN = 'epr_hyp'

    COMPLEX(real_8), ALLOCATABLE             :: rhoesg(:), rhoesgtemp(:)
    INTEGER                                  :: i, ia, iat, ierr, io, is, k, &
                                                kkk
    INTEGER, ALLOCATABLE                     :: irel(:), mesh(:), ncspvs(:), &
                                                nkkk(:,:), nlm(:,:), nnlz(:,:)
    LOGICAL                                  :: dop, dos
    REAL(real_8)                             :: cte(1000), eigr_dot, &
                                                gridft(1000,6), gridt(1000), &
                                                gyrom(99,4), rtemp, &
                                                temp(1000), tot(6)
    REAL(real_8), ALLOCATABLE :: aasf(:), ahyp(:), bbsf(:), ee(:,:), &
      exfact(:), int_1_(:,:), int_2_(:,:), int_3_(:), int_4_(:), int_5_(:), &
      int_6_(:), rmax(:), ruae(:,:), snl(:,:,:), wwnl(:,:), xion(:), z(:)

    gyrom(1:99,1)=(/ 42.576_real_8,-32.433_real_8, 62.653_real_8,  5.983_real_8,  4.574_real_8,&
         10.705_real_8,  3.077_real_8, -5.772_real_8, 40.054_real_8, -3.361_real_8,&
         11.262_real_8,  2.606_real_8, 11.094_real_8, -8.458_real_8, 17.235_real_8,&
         3.268_real_8, 1._real_8, 1._real_8, 1._real_8, 1._real_8, 1._real_8, 1._real_8, 1._real_8, 1._real_8,&
         1._real_8, 1._real_8, 1._real_8, 1._real_8, 1._real_8, 1._real_8, 1._real_8, 1._real_8, 1._real_8, 1._real_8,&
         1._real_8, 1._real_8, 1._real_8, 1._real_8, 1._real_8, 1._real_8, 1._real_8, 1._real_8, 1._real_8, 1._real_8,&
         1._real_8, 1._real_8, 1._real_8, 1._real_8, 1._real_8, 1._real_8, 1._real_8, 1._real_8, 1._real_8, 1._real_8,&
         1._real_8, 1._real_8, 1._real_8, 1._real_8, 1._real_8, 1._real_8, 1._real_8, 1._real_8, 1._real_8, 1._real_8,&
         1._real_8, 1._real_8, 1._real_8, 1._real_8, 1._real_8, 1._real_8, 1._real_8, 1._real_8, 1._real_8, 1._real_8,&
         1._real_8, 1._real_8, 1._real_8, 1._real_8, 1._real_8, 1._real_8, 1._real_8, 1._real_8, 1._real_8, 1._real_8,&
         1._real_8, 1._real_8, 1._real_8, 1._real_8, 1._real_8, 1._real_8, 1._real_8, 1._real_8, 1._real_8, 1._real_8,&
         1._real_8, 1._real_8, 1._real_8, 1._real_8, 1._real_8 /)

    ALLOCATE(aasf(ions1%nsp),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(   aasf)!,ions1%nsp)
    ALLOCATE(bbsf(ions1%nsp),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(   bbsf)!,ions1%nsp)
    ALLOCATE(rmax(ions1%nsp),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(   rmax)!,ions1%nsp)
    ALLOCATE(exfact(ions1%nsp),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(   exfact)!,ions1%nsp)
    ALLOCATE(z(ions1%nsp),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(   z)!,ions1%nsp)
    ALLOCATE(xion(ions1%nsp),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(   xion)!,ions1%nsp)
    ALLOCATE(wwnl(26,ions1%nsp),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(   wwnl)!,26*ions1%nsp)
    ALLOCATE(ee(26,ions1%nsp),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(   ee)!,26*ions1%nsp)
    ALLOCATE(ruae(1000,ions1%nsp),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(   ruae)!,1000*ions1%nsp)
    ALLOCATE(snl(1000,26,ions1%nsp),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(   snl)!,1000*26*ions1%nsp)
    ALLOCATE(nnlz(26,ions1%nsp),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(   nnlz)!,SIZE(nnlz))
    ALLOCATE(nkkk(26,ions1%nsp),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(   nkkk)!,SIZE(nkkk))
    ALLOCATE(ncspvs(ions1%nsp),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(   ncspvs)!,SIZE(ncspvs))
    ALLOCATE(mesh(ions1%nsp),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(   mesh)!,SIZE(mesh))
    ALLOCATE(irel(ions1%nsp),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(   irel)!,SIZE(irel))
    ALLOCATE(nlm(26,ions1%nsp),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(   nlm)!,SIZE(nlm))
    ALLOCATE(int_1_(6,ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(   int_1_)!,6*ions1%nat)
    ALLOCATE(int_2_(6,ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(   int_2_)!,6*ions1%nat)
    ALLOCATE(int_3_(ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(   int_3_)!,ions1%nat)
    ALLOCATE(int_4_(ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(   int_4_)!,ions1%nat)
    ALLOCATE(int_5_(ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(   int_5_)!,ions1%nat)
    ALLOCATE(int_6_(ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(   int_6_)!,ions1%nat)
    ALLOCATE(ahyp(ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(   ahyp)!,ions1%nat)
    ALLOCATE(rhoesgtemp(fpar%nnr1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(   rhoesgtemp)!,SIZE(rhoesgtemp))
    ALLOCATE(rhoesg(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(   rhoesg)!,SIZE(rhoesg))

    ! Preventing division by zero
    !$omp parallel do private(i)
#ifdef __SR8000
    !poption parallel
#endif
    DO i=1,ions1%nat
       int_3_(i)=1._real_8
    ENDDO

    IF (paral%io_parent)&
         WRITE(6,'(/,A,A)') ' ********************** DEBUG INFORMATION',&
         ' ***********************'
    ! ==--------------------------------------------------------------==
    ! Conversion of spin density to reciprocal space
    ! ==--------------------------------------------------------------==
    !$omp parallel do private(i)
#ifdef __SR8000
    !poption parallel
#endif
    DO i=1,fpar%nnr1
       rhoesgtemp(i)=CMPLX(rhoes(i),0._real_8,kind=real_8)
    ENDDO
    CALL  fwfftn(rhoesgtemp,.FALSE.,parai%allgrp)

    DO i=1,ncpw%nhg
       rhoesg(i)=rhoesgtemp(nzh(i))
    ENDDO

    ! ==--------------------------------------------------------------==
    ! Mark PSEUDO AO orbitals + Check norm
    ! ==--------------------------------------------------------------==

    ! Assign a number NLM to the orbitals
    DO is=1,ions1%nsp
       IF (.NOT.pslo_com%tvan(is)) THEN
          DO i=1, atwf_mod%nshell(is)
             nlm(i,is) = atwf_mod%nqsto(i,is)*100 + atwf_mod%lshell(i,is)*10
             DO k=1,atwr%mesh_epr(is)
                gridft(k,1)=atwfr_epr(k,i,is)**2*fpi
                gridt(k)=atrg_epr(k,is)
             ENDDO
             CALL curv1(atwr%mesh_epr(is),gridt,gridft(1,1),0.0_real_8,0.0_real_8,0,cte,&
                  temp,0.0_real_8,ierr)
             IF (paral%parent) THEN
                IF (paral%io_parent)&
                     WRITE(6,'(A,11X,T54,F12.5)')&
                     ' NORM PSEUDO ATOMIC ORBITALS:',&
                     curvi(atrg_epr(1,is),atrg_epr(atwr%mesh_epr(is),is),&
                     atwr%mesh_epr(is),gridt,gridft(1,1),cte,0.0_real_8)
                IF (paral%io_parent)&
                     WRITE(6,'(A,T27,3(8X,I5))') ' IDENTIFICATION ORBTIAL: ',&
                     is,i,nlm(i,is)
             ENDIF
          ENDDO
       ENDIF
    ENDDO

    ! ==--------------------------------------------------------------==
    ! Read All Electron AO orbitals + Check norm
    ! ==--------------------------------------------------------------==

    IF (paral%parent) THEN
       DO is=1,ions1%nsp
          IF (.NOT.pslo_com%tvan(is)) THEN
             IF (paral%io_parent)&
                  OPEN( unit = 10 , file = elem%el(ions0%iatyp(is))//'.ae', status = 'old' ,&
                  form = 'unformatted' )
             IF (paral%io_parent)&
                  READ(10) z(is),xion(is),exfact(is),mesh(is),irel(is),&
                  aasf(is),bbsf(is),rmax(is)
             IF (paral%io_parent)&
                  READ(10) ncspvs(is)
             DO k=1,ncspvs(is)
                IF (paral%io_parent)&
                     READ(10) nnlz(k,is),wwnl(k,is),ee(k,is)
             ENDDO
             IF (paral%io_parent)&
                  READ(10) (ruae(i,is),i=1,mesh(is))
             DO k=1,ncspvs(is)
                IF (paral%io_parent)&
                     READ(10) kkk,(snl(i,k,is),i=1,kkk)
                nkkk(k,is)=kkk
                DO i=1,kkk
                   snl(i,k,is)=snl(i,k,is)/SQRT(fpi)
                   gridt(i) = 1/z(is)*(EXP(-aasf(is)+&
                        (i-1)/bbsf(is))-EXP(-aasf(is)))
                   gridft(i,1)=snl(i,k,is)**2*fpi
                ENDDO
                CALL curv1(kkk,gridt,gridft(1,1),0.0_real_8,0.0_real_8,0,cte,temp,&
                     0.0_real_8,ierr)
                IF (paral%io_parent)&
                     WRITE(6,'(A,11X,T54,F12.5)')&
                     ' NORM ALL ELECTRON ATOMIC ORBITALS:',&
                     curvi(gridt(1),gridt(kkk),kkk,gridt,gridft(1,1),&
                     cte,0.0_real_8)
             ENDDO
             IF (paral%io_parent)&
                  CLOSE(10)
          ENDIF
       ENDDO
    ENDIF

    ! brodacast all data
    CALL mp_bcast(aasf,ions1%nsp,parai%source,parai%allgrp)
    CALL mp_bcast(bbsf,ions1%nsp,parai%source,parai%allgrp)
    CALL mp_bcast(rmax,ions1%nsp,parai%source,parai%allgrp)
    CALL mp_bcast(exfact,ions1%nsp,parai%source,parai%allgrp)
    CALL mp_bcast(z,ions1%nsp,parai%source,parai%allgrp)
    CALL mp_bcast(xion,ions1%nsp,parai%source,parai%allgrp)
    CALL mp_bcast(wwnl,26*ions1%nsp,parai%source,parai%allgrp)
    CALL mp_bcast(ee,26*ions1%nsp,parai%source,parai%allgrp)
    CALL mp_bcast(ruae,1000*ions1%nsp,parai%source,parai%allgrp)
    CALL mp_bcast(snl,26000*ions1%nsp,parai%source,parai%allgrp)
    CALL mp_bcast(nnlz,26*ions1%nsp,parai%source,parai%allgrp)
    CALL mp_bcast(nkkk,26*ions1%nsp,parai%source,parai%allgrp)
    CALL mp_bcast(ncspvs,ions1%nsp,parai%source,parai%allgrp)
    CALL mp_bcast(mesh,ions1%nsp,parai%source,parai%allgrp)
    CALL mp_bcast(irel,ions1%nsp,parai%source,parai%allgrp)

    ! ==--------------------------------------------------------------==
    ! Compute Hyperfine Coupling Tensor
    ! ==--------------------------------------------------------------==

    CALL mp_sync(parai%allgrp)
    ! Initialize
    iat=0                     ! index atoms
    DO is = 1,ions1%nsp
       IF (.NOT.pslo_com%tvan(is)) THEN
          DO ia = 1,ions0%na(is)
             dos=.TRUE.    ! do s-part?
             dop=.TRUE.    ! do p-part?
             iat=iat+1
             IF (paral%io_parent)&
                  WRITE(6,'(A,T32,I4)') ' STARTING CALCULATION OF ATOM ',iat
             DO io = 1,atwf_mod%nshell(is)
                IF ((dos.OR.dop).AND.paral%io_parent)&
                     WRITE(6,'(A,T63,I3)') '  NOW CHECKING ORBITAL:',nlm(io,is)
                IF (MOD(nlm(io,is),100) .EQ. 0) THEN
                   DO k=1,ncspvs(is)
                      IF ((nlm(io,is).EQ.nnlz(k,is)).AND.(dos)) THEN
                         IF (paral%io_parent)&
                              WRITE(6,*) '   -- TAKING THIS ONE FOR S-PART'
                         rtemp=1/z(is)*(EXP(-aasf(is)+&
                              REAL(2-1,kind=real_8)/bbsf(is))-EXP(-aasf(is)))
                         ahyp(iat)=-104.98_real_8*gyrom(ions0%iatyp(is),1)&
                              *dens_rad_2(iat,rhoesg,0.0_real_8,0,0)&
                              *(snl(2,k,is)/rtemp)**2/&
                              (atwfr_epr(1,io,is)/atrg_epr(1,is))**2
                         dos = .FALSE.
                      ENDIF
                   ENDDO
                ELSEIF (MOD(nlm(io,is),100) .EQ. 10) THEN
                   DO k=1,ncspvs(is)
                      IF ((nlm(io,is).EQ.nnlz(k,is)).AND.(dop)) THEN
                         IF (paral%io_parent)&
                              WRITE(6,*) '   -- TAKING THIS ONE FOR P-PART'
                         ! INT1  
                         DO i=1,ncpw%nhg
                            rhoesgtemp(i)=rhoesg(i)/hg(i)*fpi*(-1._real_8)*&
                                 (gk(1,i)*gk(1,i)-hg(i)/3._real_8)
                         ENDDO
                         int_1_(1,iat)=eigr_dot(iat,rhoesgtemp(1))
                         DO i=1,ncpw%nhg
                            rhoesgtemp(i)=rhoesg(i)/hg(i)*fpi*(-1._real_8)*&
                                 (gk(1,i)*gk(2,i))
                         ENDDO
                         int_1_(2,iat)=eigr_dot(iat,rhoesgtemp(1))
                         DO i=1,ncpw%nhg
                            rhoesgtemp(i)=rhoesg(i)/hg(i)*fpi*(-1._real_8)*&
                                 (gk(1,i)*gk(3,i))
                         ENDDO
                         int_1_(3,iat)=eigr_dot(iat,rhoesgtemp(1))
                         DO i=1,ncpw%nhg
                            rhoesgtemp(i)=rhoesg(i)/hg(i)*fpi*(-1._real_8)*&
                                 (gk(2,i)*gk(2,i)-hg(i)/3._real_8)
                         ENDDO
                         int_1_(4,iat)=eigr_dot(iat,rhoesgtemp(1))
                         DO i=1,ncpw%nhg
                            rhoesgtemp(i)=rhoesg(i)/hg(i)*fpi*(-1._real_8)*&
                                 (gk(2,i)*gk(3,i))
                         ENDDO
                         int_1_(5,iat)=eigr_dot(iat,rhoesgtemp(1))
                         DO i=1,ncpw%nhg
                            rhoesgtemp(i)=rhoesg(i)/hg(i)*fpi*(-1._real_8)*&
                                 (gk(3,i)*gk(3,i)-hg(i)/3._real_8)
                         ENDDO
                         int_1_(6,iat)=eigr_dot(iat,rhoesgtemp(1))
                         ! /INT1
                         ! INT2
                         DO i=1,ncpw%nhg
                            rhoesgtemp(i)=rhoesg(i)/hg(i)*fpi*&
                                 (gk(1,i)*gk(1,i)-hg(i)/3._real_8)&
                                 *(COS(covrad(ions0%iatyp(is))*fbohr&
                                 *parm%tpiba*SQRT(hg(i)))-1._real_8)
                         ENDDO
                         int_2_(1,iat)=eigr_dot(iat,rhoesgtemp(1))
                         DO i=1,ncpw%nhg
                            rhoesgtemp(i)=rhoesg(i)/hg(i)*fpi*&
                                 (gk(1,i)*gk(2,i))&
                                 *(COS(covrad(ions0%iatyp(is))*fbohr&
                                 *parm%tpiba*SQRT(hg(i)))-1._real_8)
                         ENDDO
                         int_2_(2,iat)=eigr_dot(iat,rhoesgtemp(1))
                         DO i=1,ncpw%nhg
                            rhoesgtemp(i)=rhoesg(i)/hg(i)*fpi*&
                                 (gk(1,i)*gk(3,i))&
                                 *(COS(covrad(ions0%iatyp(is))*fbohr&
                                 *parm%tpiba*SQRT(hg(i)))-1._real_8)
                         ENDDO
                         int_2_(3,iat)=eigr_dot(iat,rhoesgtemp(1))
                         DO i=1,ncpw%nhg
                            rhoesgtemp(i)=rhoesg(i)/hg(i)*fpi*&
                                 (gk(2,i)*gk(2,i)-hg(i)/3._real_8)&
                                 *(COS(covrad(ions0%iatyp(is))*fbohr&
                                 *parm%tpiba*SQRT(hg(i)))-1._real_8)
                         ENDDO
                         int_2_(4,iat)=eigr_dot(iat,rhoesgtemp(1))
                         DO i=1,ncpw%nhg
                            rhoesgtemp(i)=rhoesg(i)/hg(i)*fpi*&
                                 (gk(2,i)*gk(3,i))&
                                 *(COS(covrad(ions0%iatyp(is))*fbohr&
                                 *parm%tpiba*SQRT(hg(i)))-1._real_8)
                         ENDDO
                         int_2_(5,iat)=eigr_dot(iat,rhoesgtemp(1))
                         DO i=1,ncpw%nhg
                            rhoesgtemp(i)=rhoesg(i)/hg(i)*fpi*&
                                 (gk(3,i)*gk(3,i)-hg(i)/3._real_8)&
                                 *(COS(covrad(ions0%iatyp(is))*fbohr&
                                 *parm%tpiba*SQRT(hg(i)))-1._real_8)
                         ENDDO
                         int_2_(6,iat)=eigr_dot(iat,rhoesgtemp(1))
                         ! /INT2
                         ! INT3,INT5
                         DO i=1,atwr%mesh_epr(is)
                            gridt(i) = atrg_epr(i,is)
                         ENDDO
                         DO i=1,atwr%mesh_epr(is)
                            gridft(i,1)=atwfr_epr(i,io,is)**2&
                                 /(gridt(i)**3)*fpi
                         ENDDO
                         CALL curv1(atwr%mesh_epr(is),gridt,gridft(1,1),0.0_real_8,&
                              0.0_real_8,2,cte,temp,0.0_real_8,ierr)
                         int_3_(iat)=curvi(0.0_real_8,covrad(ions0%iatyp(is))*fbohr,&
                              atwr%mesh_epr(is),gridt,gridft(1,1),cte,0.0_real_8)
                         int_5_(iat)=curvi(0.0_real_8,gridt(atwr%mesh_epr(is)),&
                              atwr%mesh_epr(is),gridt,gridft(1,1),cte,0.0_real_8)
                         ! /INT3,/INT5
                         ! INT4
                         DO i=1,nkkk(k,is)-1
                            gridt(i) = 1/z(is)*(EXP(-aasf(is)+&
                                 REAL(i,kind=real_8)/bbsf(is))-EXP(-aasf(is)))
                         ENDDO
                         DO i=1,nkkk(k,is)-1
                            gridft(i,1)=snl(i+1,k,is)**2/(gridt(i)**3)*fpi
                         ENDDO
                         CALL curv1(nkkk(k,is)-1,gridt,gridft(1,1),0.0_real_8,&
                              0.0_real_8,2,cte,temp,0.0_real_8,ierr)
                         int_4_(iat)=curvi(0.0_real_8,gridt(nkkk(k,is)-1),&
                              nkkk(k,is)-1,gridt,gridft(1,1),cte,0.0_real_8)
                         ! /INT4
                         dop=.FALSE.
                         ! WRITE(6,'(A,T12,6(F9.4))')'  6 INT1: ',
                         ! &                    (int_1_(i,iat),i=1,6)
                         ! WRITE(6,'(A,T12,6(F9.4))')'  6 INT2: ',
                         ! &                    (int_2_(i,iat),i=1,6)
                         ! WRITE(6,'(A,T12,1(F9.4))')'  1 INT3: ',
                         ! &                    int_3_(iat)
                         ! WRITE(6,'(A,T12,1(F9.4))')'  1 INT4: ',
                         ! &                    int_4_(iat)
                         ! WRITE(6,'(A,T12,1(F9.4))')'  1 INT5: ',
                         ! &                    int_5_(iat)
                         ! WRITE(6,*)
                      ENDIF
                   ENDDO
                ENDIF
             ENDDO         ! io
          ENDDO             ! ia
       ELSE
          DO ia = 1,ions0%na(is)
             iat=iat+1
             ! Not yet implemented (RDeclerck)
             ! 
             ! CALL EPR_HYP_VAN(c0,rhoes(1),rhoesg(1),is,iat,nstate,
             ! &        ahyp(iat),int_1_(1,iat),int_2_(1,iat),scr,lscr,psi)
             ! ahyp(iat)=-104.98*GYROM(IATYP(IS),1)*ahyp(iat)
             int_3_(iat)=1._real_8
             int_4_(iat)=1._real_8
             int_5_(iat)=0._real_8
          ENDDO
       ENDIF
    ENDDO                     ! is

    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,*) ' *** CALCULATION COMPLETED. '
       IF (paral%io_parent)&
            WRITE(6,'(1X,64("*"))')
       IF (paral%io_parent)&
            WRITE(6,*)
       IF (paral%io_parent)&
            WRITE(6,'(/,A,A)') ' ************************ FINAL RESULTS',&
            ' *************************'
    ENDIF

    CALL mp_sum(ahyp,ions1%nat,parai%allgrp)

    iat=0                     ! index atoms
    DO is = 1,ions1%nsp
       DO ia = 1,ions0%na(is)
          iat=iat+1

          DO i=1,6
             tot(i)=-12.531_real_8*gyrom(ions0%iatyp(is),1)*(int_1_(i,iat)&
                  +int_2_(i,iat)/int_3_(iat)*(int_4_(iat)-int_5_(iat)))
          ENDDO
          CALL mp_sum(tot,6,parai%allgrp)
          IF (paral%parent) THEN
             IF (paral%io_parent)&
                  WRITE(6,*)
             IF (paral%io_parent)&
                  WRITE(6,'(A,T7,I4,A,A,A)') ' ATOM ',iat,', symbol',&
                  elem%el(ions0%iatyp(is)),'         HYPERFINE TENSOR (Mhz)'
             IF (paral%io_parent)&
                  WRITE(6,'(T18,A,3(1X,F12.6),T64,A)') '[   ',&
                  ahyp(iat)+tot(1),tot(2),tot(3),' ]'
             IF (paral%io_parent)&
                  WRITE(6,'(T18,A,3(1X,F12.6),T64,A)') '[   ',&
                  tot(2),tot(4)+ahyp(iat),tot(5),' ]'
             IF (paral%io_parent)&
                  WRITE(6,'(T18,A,3(1X,F12.6),T64,A)') '[   ',&
                  tot(3),tot(5),tot(6)+ahyp(iat),' ]'
             IF (paral%io_parent)&
                  WRITE(6,'(A,T36,F9.2)')&
                  ' ISOTROPIC INTERACTION TERM (Mhz):', ahyp(iat)
          ENDIF
       ENDDO
    ENDDO
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,'(1X,64("*"))')
       IF (paral%io_parent)&
            WRITE(6,*)
    ENDIF

    DEALLOCATE(int_1_,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(int_2_,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(int_3_,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(int_4_,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(int_5_,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(int_6_,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(ahyp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    RETURN
  END SUBROUTINE epr_hyp

  ! ==================================================================
  FUNCTION dens_rad_2(isa,func,crad,cazi,cpol) RESULT(my_res)
    ! ==--------------------------------------------------------------==
    ! INPUT:   One function f(G), packed storage, dim=NHG
    ! \        Species/Atom index ISA.
    ! OUTPUT:  function value returns
    ! \        Sum_G  exp(i G R_ion[isa])  f(G)
    ! \        i.e. the Fourier transform of f(G) at the position of
    ! \        the nucleus ISA.
    ! \        ATTENTION: The G=0 component of func is set to zero!!
    ! NOT PARALLEL -- corresponds to a DOTP-call
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: isa
    COMPLEX(real_8), DIMENSION(:), &
      INTENT(in)                             :: func
    REAL(real_8)                             :: crad
    INTEGER                                  :: cazi, cpol
    REAL(real_8)                             :: my_res

    COMPLEX(real_8)                          :: expigr
    INTEGER                                  :: ig, ig1
    REAL(real_8)                             :: dot, xx1, yy1, zz1

! ==--------------------------------------------------------------==
! CALL TISET('DENS_RAD       ',isub)

    xx1=crad*COS(cazi/20._real_8*pi)*SIN(cpol/20._real_8*pi)
    yy1=crad*SIN(cazi/20._real_8*pi)*SIN(cpol/20._real_8*pi)
    zz1=crad*COS(cpol/20._real_8*pi)

    my_res = 0._real_8
    ig1=1
    IF (geq0) THEN
       ig1=2
       expigr = ei1(isa,inyh(1,1))&
            * ei2(isa,inyh(2,1))&
            * ei3(isa,inyh(3,1))&
            * CMPLX(COS(parm%tpiba*(xx1*gk(1,1)+yy1*gk(2,1)+&
            zz1*gk(3,1))),SIN(parm%tpiba*(xx1*gk(1,1)+yy1*gk(2,1)+&
            zz1*gk(3,1))),kind=real_8)
       dot    = REAL(expigr,KIND=real_8)*REAL(func(1),KIND=real_8)&
            + AIMAG(expigr)*AIMAG(func(1))
       !dot    = real(expigr)*func(1,1) + aimag(expigr)*func(2,1)
       my_res = 0.5_real_8 * dot
    ENDIF

    !$omp parallel do private(ig,expigr,dot) reduction(+:my_res)
    DO ig=ig1,ncpw%nhg
       expigr = ei1(isa,inyh(1,ig))&
            * ei2(isa,inyh(2,ig))&
            * ei3(isa,inyh(3,ig))&
            * CMPLX(COS(parm%tpiba*(xx1*gk(1,ig)+yy1*gk(2,ig)+&
            zz1*gk(3,ig))),SIN(parm%tpiba*(xx1*gk(1,ig)+yy1*gk(2,ig)+&
            zz1*gk(3,ig))),kind=real_8)
       dot    = REAL(expigr,KIND=real_8)*REAL(func(ig),KIND=real_8)&
            + AIMAG(expigr)*AIMAG(func(ig))
       !dot    = real(expigr)*func(1,ig)&
       !     + aimag(expigr)*func(2,ig)
       my_res = my_res + dot
    ENDDO
    my_res = 2._real_8 * my_res
    ! CALL TIHALT('DENS_RAD       ',isub)
    RETURN
  END FUNCTION dens_rad_2
  ! ==================================================================

END MODULE epr_hyp_utils
