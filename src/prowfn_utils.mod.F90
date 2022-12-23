MODULE prowfn_utils
  USE adat,                            ONLY: elem
  USE atwf,                            ONLY: atwf_mod,&
                                             atwp,&
                                             catom,&
                                             loadc_foc_array_size,&
                                             xsmat,&
                                             xxmat
  USE augchg_utils,                    ONLY: augchg
  USE cmaos_utils,                     ONLY: cmaos,&
                                             give_scr_cmaos,&
                                             give_scr_satch,&
                                             satch,&
                                             satch3,&
                                             satch4
  USE dotp_utils,                      ONLY: dotp
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_def,&
                                             fo_ufo
  USE forcep_utils,                    ONLY: rhoe_psi_size
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE machine,                         ONLY: m_flush
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum,&
                                             mp_sync
  USE ovlap_utils,                     ONLY: ovlap
  USE parac,                           ONLY: parai,&
                                             paral
  USE prden,                           ONLY: numpr
  USE prmem_utils,                     ONLY: prmem
  USE prop,                            ONLY: prop1,&
                                             prop2,&
                                             prop3
  USE pslo,                            ONLY: pslo_com
  USE rnlsm_utils,                     ONLY: give_scr_rnlsm,&
                                             rnlsm
  USE setbasis_utils,                  ONLY: loadc
  USE sfac,                            ONLY: fnl,&
                                             fnl2
  USE spin,                            ONLY: spin_mod
  USE summat_utils,                    ONLY: give_scr_summat,&
                                             summat
  USE system,                          ONLY: cntl,&
                                             maxsys,&
                                             ncpw
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: dspevy,&
                                             unitmx
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: prowfn
  PUBLIC :: mklabel
  !public :: prtmat
  !public :: poploc
  !public :: locfun
  PUBLIC :: give_scr_prowfn

CONTAINS

  ! ==================================================================
  SUBROUTINE prowfn(c0,tau0,nstate)
    ! ==--------------------------------------------------------------==
    ! == PROJECT WAVEFUNCTIONS ON PSEUDO ATOMIC ORBITALS              ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(3*maxsys%nax*maxsys%nsx)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(ncpw%ngw,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'prowfn'

    CHARACTER(len=15)                        :: label(10000)
    CHARACTER(len=30)                        :: tag
    COMPLEX(real_8), ALLOCATABLE             :: cscr(:,:), psi(:,:)
    EXTERNAL                                 :: ddot
    INTEGER :: i, i1, i2, i3max, i4max, ia, ia1, ia2, iao, iao1, iao2, iaorb, &
      iat, iat1, iat2, iat3, iat3m, iat4, iat4m, ib, ic, id, ierr, ii, ijk, &
      il_psi_1d, il_psi_2d, il_rhoe_1d, il_rhoe_2d, ip, is, is1, is2, isub, &
      isw, ixx, izero, j, k, ki, kl, l, lscr, msweep, n3, n4, nao, nao1, &
      nao2, natst, nomax, numin
    LOGICAL                                  :: ferror, tlsd2
    REAL(real_8)                             :: abcm, ddot, fac, &
                                                foc(loadc_foc_array_size), &
                                                ql, qm, rlen, sfc, unac
    REAL(real_8), ALLOCATABLE :: a(:), abc(:), abcd(:), aux(:), bab(:,:), &
      comp(:), qa(:,:), rhoe(:,:), scr(:), smat(:,:), w(:), xtmat(:,:), z(:,:)

! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    IF (cntl%tfdist) CALL stopgm(procedureN,'TFDIST NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    nomax = MAX(prop2%numorb,atwp%nattot)
    ALLOCATE(comp(prop2%numorb),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(xtmat(atwp%nattot,nomax),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(xsmat(atwp%nattot,nomax),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(xxmat(atwp%nattot,nomax),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(catom(ncpw%ngw,atwp%nattot),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(cscr(ncpw%ngw,atwp%nattot),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(bab(ions1%nat,ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(qa(ions1%nat,6),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(qa)!,6*ions1%nat)
    IF (prop1%dpan) THEN
       ALLOCATE(smat(atwp%nattot,atwp%nattot),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    n3=ions1%nat*ions1%nat*ions1%nat/6 + 1
    IF (prop1%dpan.AND.prop2%ncen.GE.3) THEN
       ALLOCATE(abc(n3),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    n4=ions1%nat*ions1%nat*ions1%nat*ions1%nat/24 + 1
    IF (prop1%dpan.AND.prop2%ncen.GE.4) THEN
       ALLOCATE(abcd(n4),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! Allocation of RHOE and PSI
    CALL rhoe_psi_size(il_psi_1d=il_psi_1d,il_psi_2d=il_psi_2d,&
         il_rhoe_1d=il_rhoe_1d,il_rhoe_2d=il_rhoe_2d)
    IF (prop1%locl.AND.(numpr.NE.0)) THEN
       ALLOCATE(psi(il_psi_1d,il_psi_2d),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(rhoe(il_rhoe_1d,il_rhoe_2d),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ELSEIF (pslo_com%tivan) THEN
       ALLOCATE(psi(il_psi_1d,il_psi_2d),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! SCR ALLOCATION AND PARTITION (SCRATCH ARRAY).
    CALL give_scr_prowfn(lscr,tag)
    ALLOCATE(scr(lscr),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    IF (paral%parent) CALL prmem(procedureN)
    ! ==--------------------------------------------------------------==
    ! Load AO in PW basis
    iaorb=1
    iat=0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          iat=iat+1
          CALL loadc(catom(1,iaorb),foc,ncpw%ngw,ncpw%ngw,atwp%nattot,SIZE(foc),&
               is,iat,natst)
          DO ixx=iaorb,iaorb+natst-1
             sfc=dotp(ncpw%ngw,catom(:,ixx),catom(:,ixx))
             CALL mp_sum(sfc,parai%allgrp)
             sfc=1._real_8/SQRT(sfc)
             CALL zdscal(ncpw%ngw,sfc,catom(1,ixx),1)
          ENDDO
          iaorb=iaorb+natst
       ENDDO
    ENDDO
    ! Overlap matrix
    tlsd2=cntl%tlsd
    cntl%tlsd=.FALSE.
    CALL ovlap(atwp%nattot,xsmat,catom,catom)
    cntl%tlsd=tlsd2
    CALL summat(xsmat,atwp%nattot)
    IF (paral%io_parent) THEN
       CALL fileopen(41,'OVERLAP',fo_def+fo_ufo,ferror)
       WRITE(41)((xsmat(kl,ki),kl=1,atwp%nattot),ki=1,atwp%nattot)
       CALL fileclose(41)
    ENDIF
    CALL dcopy(atwp%nattot*atwp%nattot,xsmat(1,1),1,xtmat(1,1),1)
    ! O**(-1/2)
    ALLOCATE(a(atwp%nattot*atwp%nattot),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(w(atwp%nattot),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(z(atwp%nattot, atwp%nattot),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(aux(3*atwp%nattot),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    k=0
    DO j=1,atwp%nattot
       DO i=j,atwp%nattot
          k=k+1
          a(k)=xsmat(i,j)
       ENDDO
    ENDDO
    CALL dspevy(1,a,w,z,atwp%nattot,atwp%nattot,aux,3*atwp%nattot)
    izero=0
    DO i=1,atwp%nattot
       IF (ABS(w(i)).LT.1.e-5_real_8) THEN
          w(i)=1.0e30_real_8
          izero=izero+1
          DO j=1,atwp%nattot
             z(j,i)=0.0_real_8
          ENDDO
       ENDIF
    ENDDO
    IF ((izero.NE.0).AND.paral%io_parent)&
         WRITE(6,'(A,I3)')&
         ' NUMBER OF LINEARLY DEPENDENT ATOMIC BASIS FUNCTIONS :',izero
    DO i=1,atwp%nattot
       w(i)=1._real_8/SQRT(w(i))
    ENDDO
    CALL zeroing(xsmat)!,atwp%nattot*atwp%nattot)
    DO k=1,atwp%nattot
       DO j=1,atwp%nattot
          fac=w(k)*z(j,k)
          CALL daxpy(atwp%nattot,fac,z(1,k),1,xsmat(1,j),1)
       ENDDO
    ENDDO
    ! Orthonormalize atomic orbitals
    CALL dgemm('N','N',2*ncpw%ngw,atwp%nattot,atwp%nattot,1.0_real_8,catom(1,1),2*ncpw%ngw,&
         xsmat(1,1),atwp%nattot,0.0_real_8,cscr(1,1),2*ncpw%ngw)
    CALL dcopy(2*atwp%nattot*ncpw%ngw,cscr(1,1),1,catom(1,1),1)
    ! Overlap AO with wavefunctions
    CALL ovlap2(ncpw%ngw,atwp%nattot,prop2%numorb,xxmat,catom,c0,.TRUE.)
    CALL mp_sum(xxmat,atwp%nattot*prop2%numorb,parai%allgrp)
    IF (pslo_com%tivan) THEN
       ! Calculate augmentation charges
       CALL augchg(fnl,crge%f,qa(1,4),prop2%numorb)
       ALLOCATE(fnl2(1,ions1%nat,maxsys%nhxs,prop2%numorb,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(fnl,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(fnl(1,ions1%nat,maxsys%nhxs,atwp%nattot,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL rnlsm(catom,atwp%nattot,1,1,.FALSE.)
       CALL zeroing(fnl2)!,ions1%nat*maxsys%nhxs*prop2%numorb)
       DO i=1,prop2%numorb
          DO j=1,atwp%nattot
             CALL daxpy(ions1%nat*maxsys%nhxs,xxmat(j,i),&
                  fnl(1,1,1,j,1),1,fnl2(1,1,1,i,1),1)
          ENDDO
       ENDDO
       CALL augchg(fnl2,crge%f,qa(1,5),prop2%numorb)
    ENDIF
    DO i=1,prop2%numorb
       rlen=dotp(ncpw%ngw,c0(:,i),c0(:,i))
       CALL mp_sum(rlen,parai%allgrp)
       comp(i)=ddot(atwp%nattot,xxmat(1,i),1,xxmat(1,i),1)/rlen
    ENDDO
    IF (paral%parent) THEN
       ! Lowdin Population Analysis
       iat=0
       iao=0
       DO is=1,ions1%nsp
          nao=atwf_mod%numaor(is)
          DO ia=1,ions0%na(is)
             iat=iat+1
             ! QA(.,4) is the true augmentation charge!
             qa(iat,2)=ions0%zv(is)-qa(iat,4)
             DO i=iao+1,iao+nao
                DO j=1,nstate
                   qa(iat,2)=qa(iat,2)-crge%f(j,1)*xxmat(i,j)**2
                ENDDO
             ENDDO
             iao=iao+nao
          ENDDO
       ENDDO
    ENDIF
    ! Rotate back to AO basis
    CALL dgemm('N','N',atwp%nattot,prop2%numorb,atwp%nattot,1.0_real_8,xsmat(1,1),atwp%nattot,&
         xxmat(1,1),atwp%nattot,0.0_real_8,scr(1),atwp%nattot)
    CALL dcopy(atwp%nattot*prop2%numorb,scr(1),1,xxmat(1,1),1)
    IF (paral%io_parent) THEN
       ! Print the wavefunctions
       CALL mklabel(label)
       IF (prop1%prto) THEN
          WRITE(6,'(/,A,/)') ' WAVEFUNCTIONS IN ATOMIC ORBITAL BASIS'
          IF (cntl%tlsd) THEN
             WRITE(6,'(21X,A)') ' ****** ALPHA SPIN ******'
             CALL prtmat(xxmat(1,1),atwp%nattot,spin_mod%nsup,label,comp,crge%f)
             WRITE(6,'(/,21X,A)') ' ****** BETA  SPIN ******'
             CALL prtmat(xxmat(1,spin_mod%nsup+1),atwp%nattot,spin_mod%nsdown,label,&
                  comp(spin_mod%nsup+1),crge%f(spin_mod%nsup+1,1))
          ELSE
             CALL prtmat(xxmat,atwp%nattot,prop2%numorb,label,comp,crge%f)
          ENDIF
       ENDIF
       CALL fileopen(41,'WFNCOEF',fo_def+fo_ufo,ferror)
       WRITE(41) atwp%nattot,ions1%nsp,(ions0%zv(is),ions0%na(is),atwf_mod%numaor(is),is=1,ions1%nsp)
       WRITE(41)((xxmat(kl,ki),kl=1,atwp%nattot),ki=1,prop2%numorb)
       CALL fileclose(41)
       IF (pslo_com%tivan) THEN
          ! Print Augmentation Charges
          WRITE(6,'(/,A,/)')&
               ' VANDERBILT AUGMENTATION CHARGES '
          WRITE(6,'(A,A)') '       ATOM        FULL BASIS    ATOM BASIS'
          iat=0
          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                iat=iat+1
                WRITE(6,'(I4,4X,A,9X,F10.3,4X,F10.3)')&
                     iat,elem%el(ions0%iatyp(is)),qa(iat,4),qa(iat,5)
             ENDDO
          ENDDO
       ENDIF
       IF (prop1%mpan) THEN
          ! Mulliken Population Analysis
          DO i=1,atwp%nattot
             DO j=1,atwp%nattot
                xsmat(i,j)=0.0_real_8
                DO k=1,nstate
                   xsmat(i,j)=xsmat(i,j)+crge%f(k,1)*xxmat(i,k)*xxmat(j,k)
                ENDDO
             ENDDO
          ENDDO
          CALL zeroing(scr(1:atwp%nattot))!,atwp%nattot)
          DO j=1,atwp%nattot
             DO i=1,atwp%nattot
                scr(i)=scr(i)+xsmat(i,j)*xtmat(i,j)
             ENDDO
          ENDDO
          iat=0
          iao=0
          DO is=1,ions1%nsp
             nao=atwf_mod%numaor(is)
             DO ia=1,ions0%na(is)
                iat=iat+1
                ! QA(.,4) is the true augmentation charge!
                qa(iat,1)=ions0%zv(is)-qa(iat,4)
                DO i=iao+1,iao+nao
                   qa(iat,1)=qa(iat,1)-scr(i)
                ENDDO
                iao=iao+nao
             ENDDO
          ENDDO
          ! Bond Orders
          CALL dcopy(atwp%nattot*prop2%numorb,xxmat(1,1),1,scr(1),1)
          CALL dgemm('N','N',atwp%nattot,atwp%nattot,atwp%nattot,1.0_real_8,xsmat(1,1),&
               atwp%nattot,xtmat(1,1),atwp%nattot,0.0_real_8,xxmat(1,1),atwp%nattot)
          CALL zeroing(bab)!,ions1%nat*ions1%nat)
          iat1=0
          iao1=0
          DO is1=1,ions1%nsp
             nao1=atwf_mod%numaor(is1)
             DO ia1=1,ions0%na(is1)
                iat1=iat1+1
                DO i1=iao1+1,iao1+nao1
                   iat2=0
                   iao2=0
                   DO is2=1,ions1%nsp
                      nao2=atwf_mod%numaor(is2)
                      DO ia2=1,ions0%na(is2)
                         iat2=iat2+1
                         IF (iat2.LT.iat1) THEN
                            DO i2=iao2+1,iao2+nao2
                               bab(iat1,iat2)=bab(iat1,iat2)+&
                                    xxmat(i1,i2)*xxmat(i2,i1)
                            ENDDO
                         ENDIF
                         iao2=iao2+nao2
                      ENDDO
                   ENDDO
                ENDDO
                iao1=iao1+nao1
             ENDDO
          ENDDO
          CALL dcopy(atwp%nattot*prop2%numorb,scr(1),1,xxmat(1,1),1)
          DO i=1,ions1%nat
             DO j=i+1,ions1%nat
                bab(i,j)=bab(j,i)
             ENDDO
          ENDDO
          DO i=1,ions1%nat
             qa(i,3)=0._real_8
             DO j=1,ions1%nat
                qa(i,3)=qa(i,3)+bab(j,i)
             ENDDO
          ENDDO
          ! SPIN POPULATION
          IF (cntl%tlsd) THEN
             DO i=1,atwp%nattot
                DO j=1,atwp%nattot
                   xsmat(i,j)=0.0_real_8
                   DO k=1,nstate
                      sfc = crge%f(k,1)
                      IF (k.GT.spin_mod%nsup) sfc=-sfc
                      xsmat(i,j)=xsmat(i,j)+sfc*xxmat(i,k)*xxmat(j,k)
                   ENDDO
                ENDDO
             ENDDO
             CALL zeroing(scr(1:atwp%nattot))!,atwp%nattot)
             DO j=1,atwp%nattot
                DO i=1,atwp%nattot
                   scr(i)=scr(i)+xsmat(i,j)*xtmat(i,j)
                ENDDO
             ENDDO
             iat=0
             iao=0
             DO is=1,ions1%nsp
                nao=atwf_mod%numaor(is)
                DO ia=1,ions0%na(is)
                   iat=iat+1
                   DO i=iao+1,iao+nao
                      qa(iat,6)=qa(iat,6)-scr(i)
                   ENDDO
                   iao=iao+nao
                ENDDO
             ENDDO
          ENDIF
          ! Print Population Analysis
          WRITE(6,'(/,A,/)')&
               ' POPULATION ANALYSIS FROM PROJECTED WAVEFUNCTIONS'
          IF (cntl%tlsd) THEN
             WRITE(6,'(A)')&
                  '       ATOM      MULLIKEN(SPIN)      LOWDIN       VALENCE'
             iat=0
             DO is=1,ions1%nsp
                DO ia=1,ions0%na(is)
                   iat=iat+1
                   WRITE(6,'(I4,4X,A,4X,F7.3,F9.3,3X,F10.3,3X,F10.3)')&
                        iat,elem%el(ions0%iatyp(is)),qa(iat,1),qa(iat,6),&
                        qa(iat,2),qa(iat,3)
                ENDDO
             ENDDO
             WRITE(6,'(A,E16.8)') ' ChkSum(POP_MUL) = ',SUM(ABS(qa(:,1:3)))
          ELSE
             WRITE(6,'(A)')&
                  '       ATOM          MULLIKEN        LOWDIN       VALENCE'
             iat=0
             DO is=1,ions1%nsp
                DO ia=1,ions0%na(is)
                   iat=iat+1
                   WRITE(6,'(I4,4X,A,9X,F10.3,4X,F10.3,4X,F10.3)')&
                        iat,elem%el(ions0%iatyp(is)),qa(iat,1),qa(iat,2),qa(iat,3)
                ENDDO
             ENDDO
             WRITE(6,'(A,E16.8)') ' ChkSum(POP_MUL) = ',SUM(ABS(qa(:,1:3)))
          ENDIF
          qm=-crge%charge
          ql=-crge%charge
          DO i=1,ions1%nat
             qm=qm+qa(i,1)
             ql=ql+qa(i,2)
          ENDDO
          WRITE(6,'(A,F10.3,4X,F10.3)') '  UNASSIGNED CHARGE',qm,qL
          ! Print Bond Orders
          WRITE(6,'(/,A)')&
               ' MAYER BOND ORDERS FROM PROJECTED WAVEFUNCTIONS'
          iat=0
          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                iat=iat+1
                WRITE(label(iat),'(I4,1X,A2)') iat,elem%el(ions0%iatyp(is))
             ENDDO
          ENDDO
          msweep=ions1%nat/8
          IF (MOD(ions1%nat,8).NE.0) msweep=msweep+1
          DO isw=1,msweep
             i1=(isw-1)*8+1
             i2=MIN(8*isw,ions1%nat)
             WRITE(6,'(/,9(A8))') '        ',(label(ii),ii=i1,i2)
             DO ia=1,ions1%nat
                WRITE(6,'(A8,8(F7.3,1X))') label(ia),(bab(ia,ii),ii=i1,i2)
             ENDDO
          ENDDO
          WRITE(6,*)
       ENDIF
       IF (prop1%dpan) THEN
          ! Davidson Population Analysis
          numin=0
          DO is=1,ions1%nsp
             numin=numin+ions0%na(is)*prop2%maos(is)
             IF (atwf_mod%numaor(is).LT.prop2%maos(is)) THEN
                WRITE(6,'(A,I4)')&
                     ' INCONSISTENT ATOMIC BASIS FOR SPECIES ',iS
                CALL stopgm(procedureN,' ',& 
                     __LINE__,__FILE__)
             ENDIF
          ENDDO
          ! Density Matrix
          DO i=1,atwp%nattot
             DO j=1,atwp%nattot
                xsmat(i,j)=0.0_real_8
                DO k=1,nstate
                   xsmat(i,j)=xsmat(i,j)+crge%f(k,1)*xxmat(i,k)*xxmat(j,k)
                ENDDO
             ENDDO
          ENDDO
          unac=0._real_8
          DO iat=1,ions1%nat
             unac=unac+qa(iat,4)
          ENDDO
          ! Calculate modified atomic orbitals
          CALL dcopy(atwp%nattot*atwp%nattot,xtmat(1,1),1,smat(1,1),1)
          CALL cmaos(xsmat,smat,numin,unac)
          unac=unac/REAL(ions1%nat,kind=real_8)
          ! Calculate shared atomic charges
          CALL satch(xsmat,xtmat,smat,numin,bab)
          ! Calculate 3-center shared atomic charges
          IF (prop2%ncen.GE.3) CALL satch3(xsmat,xtmat,smat,numin,&
               bab,abc)
          ! Calculate 4-center shared atomic charges
          IF (prop2%ncen.GE.4) CALL satch4(xsmat,xtmat,smat,numin,&
               bab,abc,abcd)
          IF (prop2%ncen.GE.5)&
               WRITE(6,*) prop2%ncen,' CENTER TERMS NOT IMPLEMENTED '
          ! Charges
          iat=0
          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                iat=iat+1
                ! QA(.,4) is the true augmentation charge!
                qa(iat,1)=ions0%zv(is)-qa(iat,4)-unac
             ENDDO
          ENDDO
          iat3=0
          iat4=0
          DO ia=1,ions1%nat
             qa(ia,1)=qa(ia,1)-bab(ia,ia)
             DO ib=ia+1,ions1%nat
                qa(ia,1)=qa(ia,1)+bab(ia,ib)/2._real_8
                qa(ib,1)=qa(ib,1)+bab(ia,ib)/2._real_8
                IF (prop2%ncen.GE.3) THEN
                   DO ic=ib+1,ions1%nat
                      iat3=iat3+1
                      qa(ia,1)=qa(ia,1)-abc(iat3)/3._real_8
                      qa(ib,1)=qa(ib,1)-abc(iat3)/3._real_8
                      qa(ic,1)=qa(ic,1)-abc(iat3)/3._real_8
                      IF (prop2%ncen.GE.4) THEN
                         DO id=ic+1,ions1%nat
                            iat4=iat4+1
                            qa(ia,1)=qa(ia,1)+abcd(iat4)/4._real_8
                            qa(ib,1)=qa(ib,1)+abcd(iat4)/4._real_8
                            qa(ic,1)=qa(ic,1)+abcd(iat4)/4._real_8
                            qa(id,1)=qa(id,1)+abcd(iat4)/4._real_8
                         ENDDO
                      ENDIF
                   ENDDO
                ENDIF
             ENDDO
          ENDDO
          ! Print Population Analysis
          WRITE(6,'(/,A,/)')&
               ' POPULATION ANALYSIS FROM PROJECTED WAVEFUNCTIONS'
          WRITE(6,'(A)')&
               '       ATOM          DAVIDSON        LOWDIN       '
          iat=0
          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                iat=iat+1
                WRITE(6,'(I4,4X,A,9X,F10.3,4X,F10.3)')&
                     iat,elem%el(ions0%iatyp(is)),qa(iat,1),qa(iat,2)
             ENDDO
          ENDDO
          qm=0.0_real_8
          ql=0.0_real_8
          DO i=1,ions1%nat
             qm=qm+qa(i,1)
             ql=ql+qa(i,2)
          ENDDO
          WRITE(6,'(A,F10.3,4X,F10.3)') '  UNASSIGNED CHARGE',qm,qL
          WRITE(6,'(A,E16.8)') ' ChkSum(POP_DAV) = ',SUM(ABS(qa(:,1:2)))
          ! Print Shared Electron Numbers
          WRITE(6,'(/,A)')&
               ' SHARED ELECTRON NUMBERS FROM PROJECTED WAVEFUNCTIONS'
          iat=0
          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                iat=iat+1
                WRITE(label(iat),'(I4,1X,A2)') iat,elem%el(ions0%iatyp(is))
             ENDDO
          ENDDO
          msweep=ions1%nat/8
          IF (MOD(ions1%nat,8).NE.0) msweep=msweep+1
          DO isw=1,msweep
             i1=(isw-1)*8+1
             i2=MIN(8*isw,ions1%nat)
             WRITE(6,'(/,9(A8))') '        ',(label(ii),ii=i1,i2)
             DO ia=1,ions1%nat
                WRITE(6,'(A8,8(F7.3,1X))') label(ia),(bab(ia,ii),ii=i1,i2)
             ENDDO
          ENDDO
          WRITE(6,*)
          ! Print 3-Center Shared Electron Numbers
          IF (prop2%ncen.GE.3) THEN
             WRITE(6,'(/,A)')&
                  ' 3-CENTER SHARED ELECTRON NUMBERS FROM PROJECTED WAVEFUNCTIONS'
             WRITE(6,'(A,F10.4)') ' CUTOFF FOR PRINTING IS :',prop3%cut3o
             ip=0
             abcm=0.0_real_8
             iat3=0
             DO i=1,ions1%nat
                DO j=i+1,ions1%nat
                   DO k=j+1,ions1%nat
                      iat3=iat3+1
                      IF (ABS(abc(iat3)).GT.abcm) THEN
                         i3max=i+(j-1)*ions1%nat+(k-1)*ions1%nat*ions1%nat
                         iat3m=iat3
                         abcm=ABS(abc(iat3))
                      ENDIF
                      IF (abc(iat3).GT.prop3%cut3o) THEN
                         ip=ip+1
                         WRITE(6,'(3A8,F10.4)')&
                              label(i),label(j),label(k),abc(iat3)
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
             IF (ip.EQ.0.AND.abcm.GT.0._real_8) THEN
                ijk=i3max
                i = MOD(ijk-1,ions1%nat)+1
                k = ijk/(ions1%nat*ions1%nat)+1
                j = (ijk-(k-1)*ions1%nat*ions1%nat)/ions1%nat+1
                WRITE(6,'(3A8,F10.4)')&
                     label(i),label(j),label(k),abc(iat3m)
             ENDIF
          ENDIF
          ! Print 4-Center Shared Electron Numbers
          IF (prop2%ncen.GE.4) THEN
             WRITE(6,'(/,A)')&
                  ' 4-CENTER SHARED ELECTRON NUMBERS FROM PROJECTED WAVEFUNCTIONS'
             WRITE(6,'(A,F10.4)') ' CUTOFF FOR PRINTING IS :',prop3%cut4o
             ip=0
             abcm=0.0_real_8
             iat4=0
             DO i=1,ions1%nat
                DO j=i+1,ions1%nat
                   DO k=j+1,ions1%nat
                      DO l=k+1,ions1%nat
                         iat4=iat4+1
                         IF (ABS(abcd(iat4)).GT.abcm) THEN
                            i4max=i+(j-1)*ions1%nat+(k-1)*ions1%nat*ions1%nat+(l-1)*ions1%nat*ions1%nat*ions1%nat
                            iat4m=iat4
                            abcm=ABS(abcd(iat4))
                         ENDIF
                         IF (abcd(iat4).GT.prop3%cut4o) THEN
                            ip=ip+1
                            WRITE(6,'(4A8,F10.4)')&
                                 label(i),label(j),label(k),label(l),abcd(iat4)
                         ENDIF
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO
             IF (ip.EQ.0.AND.abcm.GT.0.0_real_8) THEN
                ijk=i4max
                i = MOD(ijk-1,ions1%nat)+1
                k = ijk/(ions1%nat*ions1%nat)+1
                j = (ijk-(k-1)*ions1%nat*ions1%nat)/ions1%nat+1
                l = (ijk-(k-1)*ions1%nat*ions1%nat-(j-1)*ions1%nat*ions1%nat*ions1%nat)/ions1%nat+1
                WRITE(6,'(4A8,F10.4)')&
                     label(i),label(j),label(k),label(l),abcd(iat4m)
             ENDIF
          ENDIF
       ENDIF
    ENDIF
    ! Clean up the memory
    DEALLOCATE(comp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(catom,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(cscr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(xsmat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(xxmat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(xtmat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(bab,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(qa,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(scr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (pslo_com%tivan) THEN
       DEALLOCATE(fnl2,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    CALL tihalt(procedureN,isub)
    CALL mp_sync(parai%allgrp)
    CALL m_flush(6)
    ! ==--------------------------------------------------------------==
    DEALLOCATE(a,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(w,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(z,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(aux,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE prowfn
  ! ==================================================================
  SUBROUTINE mklabel(label)
    ! ==--------------------------------------------------------------==
    CHARACTER(len=15)                        :: label(:)

    CHARACTER(len=1), DIMENSION(0:4), &
      PARAMETER                              :: ilab = (/'S','P','D','F','G'/)

    CHARACTER(len=3)                         :: jlab(0:4,-4:4)
    INTEGER                                  :: ia, iat, im, is, ish, k, l

! ==--------------------------------------------------------------==

    jlab(0,0)='   '
    jlab(1,-1)='x  '
    jlab(1,0)='z  '
    jlab(1,1)='y  '
    jlab(2,-2)='xy '
    jlab(2,-1)='xz '
    jlab(2,0)='z^2'
    jlab(2,1)='yz '
    jlab(2,2)='x-y'
    jlab(3,-3)='-3 '
    jlab(3,-2)='-2 '
    jlab(3,-1)='-1 '
    jlab(3,0)=' 0 '
    jlab(3,1)='+1 '
    jlab(3,2)='+2 '
    jlab(3,3)='+3 '
    jlab(4,-3)='-4 '
    jlab(4,-3)='-3 '
    jlab(4,-2)='-2 '
    jlab(4,-1)='-1 '
    jlab(4,0)=' 0 '
    jlab(4,1)='+1 '
    jlab(4,2)='+2 '
    jlab(4,3)='+3 '
    jlab(4,4)='+4 '
    iat=0
    k=0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          iat=iat+1
          DO ish=1,atwf_mod%nshell(is)
             l=atwf_mod%lshell(ish,is)
             DO im=-l,l
                k=k+1
                IF (ish.EQ.1._real_8 .AND. im.EQ.-l) THEN
                   IF (paral%io_parent)&
                        WRITE(label(k),'(I3,2X,A2,2X,A1,A3)')&
                        iat,elem%el(ions0%iatyp(is)),ilab(l),jlab(l,im)
                ELSE
                   IF (paral%io_parent)&
                        WRITE(label(k),'(9X,A1,A3)') ilab(l),jlab(l,im)
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
  END SUBROUTINE mklabel
  ! ==================================================================
  SUBROUTINE prtmat(xxmat,nattot,numorb,label,comp,f)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nattot
    REAL(real_8)                             :: xxmat(nattot,*)
    INTEGER                                  :: numorb
    CHARACTER(len=15)                        :: label(nattot)
    REAL(real_8)                             :: comp(*), f(*)

    INTEGER                                  :: i1, i2, ia, ii, isw, msweep

    msweep=numorb/8
    IF (MOD(numorb,8).NE.0) msweep=msweep+1
    DO isw=1,msweep
       i1=(isw-1)*8+1
       i2=MIN(8*isw,numorb)
       IF (paral%io_parent)&
            WRITE(6,'(/,A,8(I5,3X))') '      ORBITAL  ',(ii,ii=i1,i2)
       IF (paral%io_parent)&
            WRITE(6,'(A,8(F7.3,1X))')&
            '  COMPLETNESS  ',(comp(ii),ii=i1,i2)
       IF (paral%io_parent)&
            WRITE(6,'(A,8(F7.3,1X),/)')&
            '  OCCUPATION   ',(f(ii),ii=i1,i2)
       DO ia=1,nattot
          IF (paral%io_parent)&
               WRITE(6,'(A15,8(F7.3,1X))') label(ia),&
               (xxmat(ia,ii),ii=i1,i2)
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
  END SUBROUTINE prtmat
  ! ==================================================================
  SUBROUTINE poploc(c,s,r,norb,nao,nloc)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: norb
    REAL(real_8)                             :: r(norb,norb)
    INTEGER                                  :: nao
    REAL(real_8)                             :: c(nao,norb), s(nao,nao)
    INTEGER                                  :: nloc

    INTEGER, PARAMETER                       :: imax = 100
    REAL(real_8), PARAMETER                  :: tol = 1.e-8_real_8 

    INTEGER                                  :: i, iter, j
    REAL(real_8)                             :: a, ab, alfa, b, cg, dfun, fi, &
                                                fun, funa, sg, sin4a

! ==--------------------------------------------------------------==

    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,'(/,A,/)') '  ORBITAL LOCALIZATION '
       CALL unitmx(r,norb)
       IF (paral%io_parent)&
            WRITE(6,'(A)') '  ITER       FUNCTIONAL         DIFFERENCE '
       funa=0.0_real_8
       DO iter=1,imax
          DO i=1,nloc
             DO j=i+1,nloc
                CALL locfun(c(1,i),c(1,j),s,nao,a,b,fi)
                ab=SQRT((a*a+b*b))
                IF (ab.GT.1.e-4_real_8) THEN
                   sin4a=b/ab
                   alfa=0.25_real_8*ASIN(sin4a)
                   cg=COS(alfa)
                   sg=SIN(alfa)
                   CALL drot(nao,c(1,i),1,c(1,j),1,cg,sg)
                   CALL drot(norb,r(1,i),1,r(1,j),1,cg,sg)
                ENDIF
             ENDDO
          ENDDO
          fun=0.0_real_8
          DO i=1,nloc
             CALL locfun(c(1,i),c(1,i),s,nao,a,b,fi)
             fun=fun+fi
          ENDDO
          dfun=funa-fun
          IF (paral%io_parent)&
               WRITE(6,'(I5,5X,F12.6,5X,F16.10)') iter,fun,dfun
          funa=fun
          IF (ABS(dfun).LT.tol) GOTO 100
       ENDDO
       IF (paral%io_parent)&
            WRITE(6,'(A)') ' NO CONVERGENCE IN LOCALIZATION PROCEDURE'
100    CONTINUE
    ENDIF
    CALL mp_sync(parai%allgrp)
    CALL mp_bcast(c,SIZE(c),parai%source,parai%allgrp)
    CALL mp_bcast(s,SIZE(s),parai%source,parai%allgrp)
    CALL mp_bcast(r,SIZE(r),parai%source,parai%allgrp)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE poploc
  ! ==================================================================
  SUBROUTINE locfun(cs,ct,s,nao,a,b,f)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: cs(*), ct(*)
    INTEGER                                  :: nao
    REAL(real_8)                             :: s(nao,nao), a, b, f

    INTEGER                                  :: i, ia, iao, is, j, nsao
    REAL(real_8)                             :: pss, pst, ptt

    a = 0.0_real_8
    b = 0.0_real_8
    f = 0.0_real_8
    iao=0
    DO is=1,ions1%nsp
       nsao=atwf_mod%numaor(is)
       DO ia=1,ions0%na(is)
          pss=0.0_real_8
          pst=0.0_real_8
          ptt=0.0_real_8
          DO j=1,nao
             DO i=iao+1,iao+nsao
                pss=pss+0.5_real_8*(cs(j)*s(j,i)*cs(i)+cs(i)*s(i,j)*cs(j))
                pst=pst+0.5_real_8*(cs(j)*s(j,i)*ct(i)+cs(i)*s(i,j)*ct(j))
                ptt=ptt+0.5_real_8*(ct(j)*s(j,i)*ct(i)+ct(i)*s(i,j)*ct(j))
             ENDDO
          ENDDO
          iao=iao+nsao
          a=a+pst*pst-0.25_real_8*(pss-ptt)**2
          b=b+pst*(pss-ptt)
          f=f+pss*pss
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
  END SUBROUTINE locfun
  ! ==================================================================
  SUBROUTINE give_scr_prowfn(lprowfn,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lprowfn
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: is, lcmaos, lrnlsm, lsatch, &
                                                lsummat, lwfnrho, nstate, &
                                                numin

    nstate=crge%n
    lrnlsm=0
    lcmaos=0
    lsatch=0
    lwfnrho=0
    CALL give_scr_summat(lsummat,tag,atwp%nattot)
    IF (pslo_com%tivan) CALL give_scr_rnlsm(lrnlsm,tag,atwp%nattot,.FALSE.)
    IF (prop1%dpan) THEN
       ! Davidson Population Analysis
       numin=0
       DO is=1,ions1%nsp
          numin=numin+ions0%na(is)*prop2%maos(is)
       ENDDO
       CALL give_scr_cmaos(lcmaos,tag,numin)
       CALL give_scr_satch(lsatch,tag,numin)
    ENDIF
    ! PROWFN
    lprowfn=MAX(2*atwp%nattot*atwp%nattot+4*atwp%nattot,&
         atwp%nattot*prop2%numorb,&
         lrnlsm,lsummat,lcmaos,lsatch,lwfnrho)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE give_scr_prowfn
  ! ==================================================================

END MODULE prowfn_utils
