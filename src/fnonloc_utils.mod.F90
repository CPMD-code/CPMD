#include "cpmd_global.h"

MODULE fnonloc_utils
  USE cp_grp_utils,                    ONLY: cp_grp_get_sizes,&
                                             cp_grp_redist
  USE cppt,                            ONLY: twnl
  USE cvan,                            ONLY: deeq,&
                                             dvan
  USE ener,                            ONLY: ener_d
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE kpnt,                            ONLY: eigkr
  USE kpts,                            ONLY: tkpts
  USE mp_interface,                    ONLY: mp_sum
  USE nlps,                            ONLY: imagp,&
                                             nghtol,&
                                             nlps_com,&
                                             wsg
  USE nvtx_utils
  USE parac,                           ONLY: parai
  USE pslo,                            ONLY: pslo_com
  USE reshaper,                        ONLY: reshape_inplace
  USE sfac,                            ONLY: eigr,&
                                             fnl,&
                                             fnl2
  USE sgpp,                            ONLY: sgpp1,&
                                             sgpp2
  USE spin,                            ONLY: clsd,&
                                             lspin2,&
                                             spin_mod
  USE system,                          ONLY: cntl,&
                                             iatpe,&
                                             ipept,&
                                             maxsys,&
                                             ncpw,&
                                             nkpt,&
                                             parap
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: fnonloc
  PUBLIC :: give_scr_fnonloc
  !public :: fcasnl

CONTAINS

  ! ==================================================================
  SUBROUTINE fnonloc(c2,f,nstate,ikind,ispin,redist_c2)
    ! ==--------------------------------------------------------------==
    ! == CALCULATES THE NON-LOCAL PP CONTRIBUTION OF THE HAMILTONIAN  ==
    ! ==--------------------------------------------------------------==
    ! == INPUT/OUPUT:                                                 ==
    ! ==   C2(NGWK,NSTATE) Wavefunctions, output C2 = C2 + |N.L.PART> ==
    ! == INPUT:                                                       ==
    ! ==   NSTATE      Number of states                               ==
    ! ==   F(1:NSTATE) Occupation numbers                             ==
    ! ==   IKIND  Index of k point                                    ==
    ! ==   ISPIN  Need with LSD option for diagonalization scheme     ==
    ! ==          Does not work with cntl%tdiag (ISPIN=1) and TIVAN        ==
    ! == SCRATCH                                                      ==
    ! ==   AUXC(NGWK)                                                 ==
    ! ==   DDIA(IMAGP*maxsys%nax)                                            ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    REAL(real_8)                             :: f(nstate)
    COMPLEX(real_8)                          :: c2(nkpt%ngwk,nstate)
    INTEGER                                  :: ikind, ispin
    LOGICAL                                  :: redist_c2

    CHARACTER(*), PARAMETER                  :: procedureN = 'fnonloc'
    COMPLEX(real_8), PARAMETER               :: zone = (1._real_8,0._real_8), &
                                                zzero = (0._real_8,0._real_8)

    COMPLEX(real_8)                          :: ctm
    COMPLEX(real_8), ALLOCATABLE             :: auxc(:,:), C2_local(:,:)
    INTEGER :: i, i_g0_l, ia, iaa, ibeg_c0, iend_c0, ierr, ig, is, isa, isa0, &
      isub, isub2, isub3, iv, jv, ki, kj, l, l2, li, lj, lspin, NGW_local, &
      NGWK_local
    LOGICAL                                  :: GEQ0_local
    REAL(real_8)                             :: ffi, t, t1, ti, tr
    REAL(real_8), ALLOCATABLE, TARGET        :: ddia(:,:)
    REAL(real_8), POINTER                    :: ddki(:,:,:)

! Variables
! ==--------------------------------------------------------------==
! ==--------------------------------------------------------------==
! == Compute the force on the electronic degrees of freedom due   ==
! == to the non-local part of the potential, and add it to the    ==
! == other piece, coming from the local contribution.             ==
! ==--------------------------------------------------------------==

    IF (nkpt%ngwk.EQ.0 .AND. .NOT.cntl%tfdist) RETURN
    CALL tiset(procedureN,isub)
    __NVTX_TIMER_START ( procedureN )

    ALLOCATE(auxc(nkpt%ngwk, nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(ddia(imagp*maxsys%nax, nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL reshape_inplace(ddia, (/2, maxsys%nax, nstate/), ddki)
    ! >>>>>>> cp_grp trick
    CALL tiset(procedureN//'_grps_a',isub2)
    CALL cp_grp_get_sizes(NGWK_l=NGWK_local,NGW_l=NGW_local,&
         GEQ0_l=GEQ0_local,FIRSTK_g=ibeg_c0,LASTK_g=iend_c0,&
         i_g0_l=i_g0_l)
    ALLOCATE(C2_local(NGWK_local,nstate),stat=ierr)
    IF (ierr.NE.0) CALL stopgm(procedureN,'Allocation problem',& 
         __LINE__,__FILE__)
    CALL zeroing( C2_local)!, NGWK_local * nstate )
    ! vw      CALL cp_grp_copy_wfn_to_local(C2,NGWK,C2_local,NGWK_local,
    ! vw     &     ibeg_c0,NGWK_local,NSTATE)
    CALL tihalt(procedureN//'_grps_a',isub2)
    ! <<<<<<<

    isa0=0
    DO is=1,ions1%nsp
       IF (pslo_com%tvan(is)) THEN
          ! ==--------------------------------------------------------------==
          ! ==  VANDERBILT PP                                               ==
          ! ==--------------------------------------------------------------==
          lspin=1
          DO iv=1,nlps_com%ngh(is)
             CALL zeroing(ddia)!,imagp*maxsys%nax*nstate)
#if defined(__ES)
             !$omp parallel do private(I,IA,JV,LSPIN,ISA,IAA)
#endif
#ifdef __SR8000
             !poption parallel, tlocal(I,IA,JV,LSPIN,ISA,IAA)
#endif 
             DO i=1,nstate
                IF (cntl%tlsd.AND.ispin.EQ.2) THEN
                   IF (i.LE.spin_mod%nsup) THEN
                      lspin=1
                   ELSE
                      lspin=2
                   ENDIF
                ENDIF
                DO ia=1,ions0%na(is)
                   isa=isa0+ia
                   IF (cntl%tfdist) THEN
                      IF (iatpe(isa).EQ.parai%mepos) THEN
                         iaa=isa-ipept(1,parai%mepos)+1
                         DO jv=1,nlps_com%ngh(is)
                            ddia(ia,i)=ddia(ia,i) +&
                                 fnl(1,iaa,jv,i,ikind)*(deeq(isa,iv,jv,lspin)+&
                                 dvan(iv,jv,is))
                         ENDDO
                      ENDIF
                   ELSE
                      DO jv=1,nlps_com%ngh(is)
                         ddia(ia,i)=ddia(ia,i)+&
                              fnl(1,isa,jv,i,ikind)*(deeq(isa,iv,jv,lspin)+&
                              dvan(iv,jv,is))
                      ENDDO
                   ENDIF
                ENDDO
             ENDDO
             IF (cntl%tfdist) CALL mp_sum(ddia,imagp*maxsys%nax*nstate,parai%allgrp)
             IF (NGWK_local.GT.0) THEN
                IF (tkpts%tkpnt) THEN
                   CALL zgemm("N","N",NGWK_local,nstate,ions0%na(is),zone,&
                        eigkr(ibeg_c0,isa0+1,ikind),nkpt%ngwk,ddki(1,1,1),&
                        maxsys%nax,zzero,auxc,nkpt%ngwk)
                ELSE
                   CALL dgemm("N","N",2*NGW_local,nstate,ions0%na(is),1._real_8,&
                        eigr(ibeg_c0,isa0+1,1),2*ncpw%ngw,ddia(1,1),maxsys%nax,0._real_8,&
                        auxc,2*ncpw%ngw)
                ENDIF
                !$omp parallel do private(I,IG,FFI,CTM,T1,TR,TI)
#ifdef __SR8000
                !poption parallel, tlocal(I,IG,FFI,CTM,T1,TR,TI)
#endif 
                DO i=1,nstate
                   ffi=f(i)
                   IF (ffi.LT.1.e-5_real_8) ffi=1.0_real_8
                   ctm=-ffi*(0.0_real_8,-1.0_real_8)**nghtol(iv,is)
                   ! MAKE USE OF THE SPECIAL STRUCTURE OF CTM
                   IF (ABS(REAL(ctm)).GT.0.01_real_8) THEN
                      ! CTM IS REAL
                      DO ig=1,NGWK_local
                         t1=twnl(ig+ibeg_c0-1,iv,is,ikind)*REAL(ctm)
                         tr=REAL(auxc(ig,i))
                         ti=AIMAG(auxc(ig,i))
                         C2_local(ig,i)=c2_local(ig,i)+CMPLX(t1*tr,t1*ti,kind=real_8)
                      ENDDO
                   ELSE
                      ! CTM IS IMAGINARY
                      DO ig=1,NGWK_local
                         t1=twnl(ig+ibeg_c0-1,iv,is,ikind)*AIMAG(ctm)
                         tr=REAL(auxc(ig,i))
                         ti=AIMAG(auxc(ig,i))
                         C2_local(ig,i)=c2_local(ig,i)+CMPLX(-t1*ti,t1*tr,kind=real_8)
                      ENDDO
                   ENDIF
                   ! IF(TKPNT.AND.GEQ0_local) C2_local(1+NGW_local,I)=ZZERO
                   IF (tkpts%tkpnt.AND.GEQ0_local) C2_local(i_g0_l,i)=zzero
                ENDDO
             ENDIF
          ENDDO
       ELSE IF (sgpp1%tsgp(is)) THEN
          ! ==--------------------------------------------------------------==
          ! ==  Stefan Goedeckers PP                                        ==
          ! ==--------------------------------------------------------------==
          DO iv=1,nlps_com%ngh(is)
             CALL zeroing(ddia)!,imagp*maxsys%nax*nstate)
             l=nghtol(iv,is)+1
             ki=sgpp2%lfval(iv,is)
             li=sgpp2%lpval(iv,is)
             DO jv=1,nlps_com%ngh(is)
                l2=nghtol(jv,is)+1
                lj=sgpp2%lpval(jv,is)
                IF (l.EQ.l2.AND.li.EQ.lj) THEN
                   kj=sgpp2%lfval(jv,is)
                   !$omp parallel do private(I,IA,ISA,IAA)
#ifdef __SR8000
                   !poption parallel, tlocal(I,IA,ISA,IAA)
#endif 
                   DO i=1,nstate
#ifdef _vpp_
                      !OCL NOALIAS 
#endif 
                      DO ia=1,ions0%na(is)
                         isa=isa0+ia
                         IF (cntl%tfdist) THEN
                            IF (iatpe(isa).EQ.parai%mepos) THEN
                               iaa=isa-ipept(1,parai%mepos)+1
                               IF (tkpts%tkpnt) THEN
                                  ddki(1,ia,i)=ddki(1,ia,i)+&
                                       fnl(1,iaa,jv,i,ikind)*sgpp2%hlsg(ki,kj,l,is)
                                  ddki(2,ia,i)=ddki(2,ia,i)+&
                                       fnl(2,iaa,jv,i,ikind)*sgpp2%hlsg(ki,kj,l,is)
                               ELSE
                                  ddia(ia,i)=ddia(ia,i)+&
                                       fnl(1,iaa,jv,i,ikind)*sgpp2%hlsg(ki,kj,l,is)
                               ENDIF
                            ENDIF
                         ELSE
                            IF (tkpts%tkpnt) THEN
                               ddki(1,ia,i)=ddki(1,ia,i)+&
                                    fnl(1,isa,jv,i,ikind)*sgpp2%hlsg(ki,kj,l,is)
                               ddki(2,ia,i)=ddki(2,ia,i)+&
                                    fnl(2,isa,jv,i,ikind)*sgpp2%hlsg(ki,kj,l,is)
                            ELSE
                               ddia(ia,i)=ddia(ia,i)+&
                                    fnl(1,isa,jv,i,ikind)*sgpp2%hlsg(ki,kj,l,is)
                            ENDIF
                         ENDIF
                      ENDDO
                   ENDDO
                ENDIF
             ENDDO
             IF (cntl%tfdist) CALL mp_sum(ddia,imagp*maxsys%nax*nstate,parai%allgrp)

             IF (NGWK_local.GT.0) THEN
                ! vw here we can build only the submatrix
                IF (tkpts%tkpnt) THEN
                   CALL zgemm("N","N",NGWK_local,nstate,ions0%na(is),zone,&
                        eigkr(ibeg_c0,isa0+1,ikind),nkpt%ngwk,ddki(1,1,1),&
                        maxsys%nax,zzero,auxc,nkpt%ngwk)
                ELSE
                   CALL dgemm("N","N",2*NGW_local,nstate,ions0%na(is),1._real_8,&
                        eigr(ibeg_c0,isa0+1,1),2*ncpw%ngw,ddia(1,1),maxsys%nax,0._real_8,&
                        auxc,2*ncpw%ngw)
                ENDIF
                !$omp parallel do private(I,IG,FFI,CTM,T1,TR,TI)
#ifdef __SR8000
                !poption parallel, tlocal(I,IG,FFI,CTM,T1,TR,TI)
#endif 
                DO i=1,nstate
                   ffi=f(i)
                   IF (ffi.LT.1.e-5_real_8) ffi=1.0_real_8
                   ctm=-ffi*(0.0_real_8,-1.0_real_8)**nghtol(iv,is)
                   ! MAKE USE OF THE SPECIAL STRUCTURE OF CTM
                   IF (ABS(REAL(ctm)).GT.0.01_real_8) THEN
                      ! CTM IS REAL
                      DO ig=1,NGWK_local
                         t1=twnl(ig+ibeg_c0-1,iv,is,ikind)*REAL(ctm)
                         tr=REAL(auxc(ig,i))
                         ti=AIMAG(auxc(ig,i))
                         C2_local(ig,i)=c2_local(ig,i)+CMPLX(t1*tr,t1*ti,kind=real_8)
                      ENDDO
                   ELSE
                      ! CTM IS IMAGINARY
                      DO ig=1,NGWK_local
                         t1=twnl(ig+ibeg_c0-1,iv,is,ikind)*AIMAG(ctm)
                         tr=REAL(auxc(ig,i))
                         ti=AIMAG(auxc(ig,i))
                         C2_local(ig,i)=c2_local(ig,i)+CMPLX(-t1*ti,t1*tr,kind=real_8)
                      ENDDO
                   ENDIF
                   ! IF(TKPNT.AND.GEQ0_local) 
                   ! &              C2_local(1+NGW_local,I)=ZZERO
                   IF (tkpts%tkpnt.AND.GEQ0_local) C2_local(i_g0_l,i)=zzero
                ENDDO        ! NSTATE
             ENDIF          ! NGWK.GT.0
          ENDDO              ! IV
       ELSE
          ! ==--------------------------------------------------------------==
          ! ==  BACHELET HAMANN SCHLUTER AND TROULLIER MARTINS              ==
          ! ==--------------------------------------------------------------==
          DO iv=1,nlps_com%ngh(is)
             IF (cntl%tfdist) THEN
                CALL zeroing(ddia)!,imagp*maxsys%nax*nstate)
                DO i=1,nstate
                   IF (i.GE.parap%nst12(parai%mepos,1).AND.i.LE.parap%nst12(parai%mepos,2)) THEN
                      iaa=i-parap%nst12(parai%mepos,1)+1
                      CALL dcopy(imagp*ions0%na(is),fnl2(1,isa0+1,iv,iaa,ikind),1,&
                           ddia(1,i),1)
                   ENDIF
                ENDDO
                CALL mp_sum(ddia,imagp*maxsys%nax*nstate,parai%allgrp)
             ELSE
                DO i=1,nstate
                   CALL dcopy(imagp*ions0%na(is),fnl(1,isa0+1,iv,i,ikind),1,&
                        ddia(1,i),1)
                ENDDO
             ENDIF
             IF (NGWK_local.GT.0) THEN
                IF (tkpts%tkpnt) THEN
                   CALL zgemm('N','N',NGWK_local,nstate,ions0%na(is),zone,&
                        eigkr(ibeg_c0,isa0+1,ikind),nkpt%ngwk,ddia,maxsys%nax,&
                        zzero,auxc,nkpt%ngwk)
                ELSE
                   CALL dgemm('N','N',2*NGW_local,nstate,ions0%na(is),1._real_8,&
                        eigr(ibeg_c0,isa0+1,1),2*ncpw%ngw,ddia,maxsys%nax,0.0_real_8,&
                        auxc,2*ncpw%ngw)
                ENDIF
                !$omp parallel do private(I,IG,CTM,FFI,T)
#ifdef __SR8000
                !poption parallel, tlocal(I,IG,CTM,FFI,T)
#endif 
                DO i=1,nstate
                   ctm=(0.0_real_8,-1.0_real_8)**nghtol(iv,is)
                   ffi=f(i)
                   IF (ffi.LT.1.e-5_real_8) ffi=1.0_real_8
                   t=-ffi*wsg(is,iv)
                   ctm=ctm*t
                   DO ig=1,NGWK_local
                      auxc(ig,i)=ctm*auxc(ig,i)
                      C2_local(ig,i)=c2_local(ig,i)+&
                           twnl(ig+ibeg_c0-1,iv,is,ikind)*auxc(ig,i)
                   ENDDO
                ENDDO
             ENDIF
          ENDDO
          ! ==--------------------------------------------------------------==
       ENDIF
       isa0=isa0+ions0%na(is)
    ENDDO

    ! NL PP FORCE CONTRIBUTION TO CAS22 METHOD
    IF (lspin2%tlse .AND. lspin2%tcas22) THEN
       IF (parai%cp_nogrp.GT.1) CALL stopgm(procedureN,'more work here',& 
            __LINE__,__FILE__)
       CALL fcasnl(C2_local,auxc,ddia)
    ENDIF

    ! >>>>>>> cp_grp trick
    CALL tiset(procedureN//'_grps_b',isub3)
    ! we need to zero the C0 that we can do the reduce
    ! vw      CALL cp_grp_copy_local_to_wfn(C2_local,NGWK_local,C2,NGWK,
    ! vw     &     ibeg_c0,NGWK_local,NSTATE)
    CALL cp_grp_add_local_to_wfn(C2_local,NGWK_local,c2,nkpt%ngwk,&
         ibeg_c0,NGWK_local,nstate)
    ! vw    IF(TKPNT.AND.GEQ0) write(6,*) SUM(ABS(C2(1+NGW,1:NSTATE)))
    IF (redist_c2) THEN
       CALL cp_grp_zero_g(c2,nkpt%ngwk,nstate,ibeg_c0,iend_c0)
       ! we reduce so that we get back the cp_grp distribution of C0
       CALL cp_grp_redist(c2,nkpt%ngwk,nstate)
    ENDIF
    DEALLOCATE(C2_local,stat=ierr)
    IF (ierr.NE.0) CALL stopgm(procedureN,'Deallocation problem',& 
         __LINE__,__FILE__)
    CALL tihalt(procedureN//'_grps_b',isub3)
    ! <<<<<<<

    ! 
    DEALLOCATE(auxc,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(ddia,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)

    __NVTX_TIMER_STOP
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE fnonloc
  ! ==================================================================
  SUBROUTINE give_scr_fnonloc(il_auxc,il_ddia,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: il_auxc, il_ddia, nstate

! ==--------------------------------------------------------------==

    il_auxc=2*nkpt%ngwk*nstate+2
    il_ddia=2*maxsys%nax*nstate+2
    ! ==--------------------------------------------------------------==
  END SUBROUTINE give_scr_fnonloc
  ! ==================================================================
  SUBROUTINE fcasnl(c2,auxc,ddia)
    ! ==--------------------------------------------------------------==
    ! == CALCULATES THE NON-LOCAL PP CONTRIBUTION OF THE HAMILTONIAN  ==
    ! == SPECIAL TERMS FOR THE CAS22 METHOD                           ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c2(nkpt%ngwk,*), &
                                                auxc(nkpt%ngwk)
    REAL(real_8)                             :: ddia(maxsys%nax)

    COMPLEX(real_8)                          :: ctm
    INTEGER                                  :: ia, iaa, ig, is, isa, isa0, &
                                                iv, jv, ki, kj, l, l2, li, lj
    REAL(real_8)                             :: ca, ctest, dd, sa, t, t1, ti, &
                                                tr

    IF (nkpt%ngwk.EQ.0) RETURN
    IF (tkpts%tkpnt) CALL stopgm("FCASNL","NOT IMPLEMENTED",& 
         __LINE__,__FILE__)
    IF (pslo_com%tivan) CALL stopgm("FCASNL","NOT IMPLEMENTED",& 
         __LINE__,__FILE__)
    ca=COS(ener_d%casang)
    sa=SIN(ener_d%casang)
    isa0=0
    DO is=1,ions1%nsp
       IF (sgpp1%tsgp(is)) THEN
          ! Stefan Goedeckers PP
          DO iv=1,nlps_com%ngh(is)
             ! STATE=ALPHA
             CALL zeroing(auxc)!,nkpt%ngwk)
             DO ia=1,ions0%na(is)
                isa=isa0+ia
                l=nghtol(iv,is)+1
                ki=sgpp2%lfval(iv,is)
                li=sgpp2%lpval(iv,is)
                dd=0.0_real_8
                DO jv=1,nlps_com%ngh(is)
                   l2=nghtol(jv,is)+1
                   lj=sgpp2%lpval(jv,is)
                   IF (l.EQ.l2.AND.li.EQ.lj) THEN
                      kj=sgpp2%lfval(jv,is)
                      IF (cntl%tfdist) THEN
                         IF (iatpe(isa).EQ.parai%mepos) THEN
                            iaa=isa-ipept(1,parai%mepos)+1
                            dd=dd+(-0.5_real_8*ca*sa*fnl(1,iaa,jv,clsd%ibeta,1)+sa*sa*&
                                 fnl(1,iaa,jv,clsd%ialpha,1))*sgpp2%hlsg(ki,kj,l,is)
                         ENDIF
                      ELSE
                         dd=dd+(-0.5_real_8*ca*sa*fnl(1,isa,jv,clsd%ibeta,1)+&
                              sa*sa*fnl(1,isa,jv,clsd%ialpha,1))*sgpp2%hlsg(ki,kj,l,is)
                      ENDIF
                   ENDIF
                ENDDO
                IF (cntl%tfdist) CALL mp_sum(dd,parai%allgrp)
                CALL daxpy(2*ncpw%ngw,dd,eigr(1,isa,1),1,auxc(1),1)
             ENDDO
             ctm=-(0.0_real_8,-1.0_real_8)**nghtol(iv,is)
             ! MAKE USE OF THE SPECIAL STRUCTURE OF CTM
             IF (ABS(REAL(ctm)).GT.0.01_real_8) THEN
                ! CTM IS REAL
                !$omp  parallel do private(IG,T1,TR,TI) 
                DO ig=1,nkpt%ngwk
                   t1=twnl(ig,iv,is,1)*REAL(ctm)
                   tr=REAL(auxc(ig))
                   ti=AIMAG(auxc(ig))
                   c2(ig,clsd%ialpha)=c2(ig,clsd%ialpha)+CMPLX(t1*tr,t1*ti,kind=real_8)
                ENDDO
             ELSE
                ! CTM IS IMAGINARY
                !$omp  parallel do private(IG,T1,TR,TI) 
                DO ig=1,nkpt%ngwk
                   t1=twnl(ig,iv,is,1)*AIMAG(ctm)
                   tr=REAL(auxc(ig))
                   ti=AIMAG(auxc(ig))
                   c2(ig,clsd%ialpha)=c2(ig,clsd%ialpha)+CMPLX(-t1*ti,t1*tr,kind=real_8)
                ENDDO
             ENDIF
             ! STATE=BETA
             CALL zeroing(auxc)!,nkpt%ngwk)
             DO ia=1,ions0%na(is)
                isa=isa0+ia
                l=nghtol(iv,is)+1
                ki=sgpp2%lfval(iv,is)
                li=sgpp2%lpval(iv,is)
                dd=0.0_real_8
                DO jv=1,nlps_com%ngh(is)
                   l2=nghtol(jv,is)+1
                   lj=sgpp2%lpval(jv,is)
                   IF (l.EQ.l2.AND.li.EQ.lj) THEN
                      kj=sgpp2%lfval(jv,is)
                      IF (cntl%tfdist) THEN
                         IF (iatpe(isa).EQ.parai%mepos) THEN
                            iaa=isa-ipept(1,parai%mepos)+1
                            dd=dd+(-0.5_real_8*ca*sa*fnl(1,iaa,jv,clsd%ialpha,1)-&
                                 sa*sa*fnl(1,iaa,jv,clsd%ibeta,1))*sgpp2%hlsg(ki,kj,l,is)
                         ENDIF
                      ELSE
                         dd=dd+(-0.5_real_8*sa*ca*fnl(1,isa,jv,clsd%ialpha,1)-&
                              sa*sa*fnl(1,isa,jv,clsd%ibeta,1))*sgpp2%hlsg(ki,kj,l,is)
                      ENDIF
                   ENDIF
                ENDDO
                IF (cntl%tfdist) CALL mp_sum(dd,parai%allgrp)
                CALL daxpy(2*ncpw%ngw,dd,eigr(1,isa,1),1,auxc(1),1)
             ENDDO
             ctm=-(0.0_real_8,-1.0_real_8)**nghtol(iv,is)
             ! MAKE USE OF THE SPECIAL STRUCTURE OF CTM
             IF (ABS(REAL(ctm)).GT.0.01_real_8) THEN
                ! CTM IS REAL
                !$omp  parallel do private(IG,T1,TR,TI) 
                DO ig=1,nkpt%ngwk
                   t1=twnl(ig,iv,is,1)*REAL(ctm)
                   tr=REAL(auxc(ig))
                   ti=AIMAG(auxc(ig))
                   c2(ig,clsd%ibeta)=c2(ig,clsd%ibeta)+CMPLX(t1*tr,t1*ti,kind=real_8)
                ENDDO
             ELSE
                ! CTM IS IMAGINARY
                !$omp  parallel do private(IG,T1,TR,TI) 
                DO ig=1,nkpt%ngwk
                   t1=twnl(ig,iv,is,1)*AIMAG(ctm)
                   tr=REAL(auxc(ig))
                   ti=AIMAG(auxc(ig))
                   c2(ig,clsd%ibeta)=c2(ig,clsd%ibeta)+CMPLX(-t1*ti,t1*tr,kind=real_8)
                ENDDO
             ENDIF
          ENDDO
       ELSE
          ! BACHELET HAMANN SCHLUTER
          DO iv=1,nlps_com%ngh(is)
             ! STATE=ALPHA
             IF (cntl%tfdist) THEN
                CALL zeroing(ddia(1:ions0%na(is)))!,ions0%na(is))
                IF (clsd%ibeta.GE.parap%nst12(parai%mepos,1).AND.&
                     clsd%ibeta.LE.parap%nst12(parai%mepos,2)) THEN
                   iaa=clsd%ibeta-parap%nst12(parai%mepos,1)+1
                   CALL dcopy(ions0%na(is),fnl2(1,isa0+1,iv,iaa,1),1,ddia,1)
                   CALL dscal(ions0%na(is),-0.5_real_8*ca*sa,ddia,1)
                ENDIF
                IF (clsd%ialpha.GE.parap%nst12(parai%mepos,1).AND.&
                     clsd%ialpha.LE.parap%nst12(parai%mepos,2)) THEN
                   iaa=clsd%ialpha-parap%nst12(parai%mepos,1)+1
                   CALL daxpy(ions0%na(is),sa*sa,fnl2(1,isa0+1,iv,iaa,1),&
                        1,ddia,1)
                ENDIF
                CALL mp_sum(ddia,ions0%na(is),parai%allgrp)
             ELSE
                CALL dcopy(ions0%na(is),fnl(1,isa0+1,iv,clsd%ibeta,1),1,ddia,1)
                CALL dscal(ions0%na(is),-0.5_real_8*ca*sa,ddia,1)
                CALL daxpy(ions0%na(is),sa*sa,fnl(1,isa0+1,iv,clsd%ialpha,1),&
                     1,ddia,1)
             ENDIF
             CALL dgemv('N',2*ncpw%ngw,ions0%na(is),1.0_real_8,eigr(1,isa0+1,1),2*ncpw%ngw,&
                  ddia,1,0.0_real_8,auxc(1),1)
             t=-wsg(is,iv)
             ctm=(0.0_real_8,-1.0_real_8)**nghtol(iv,is)
             ctest=REAL(ctm)
             ! MAKE USE OF THE SPECIAL STRUCTURE OF CTM
             IF (ABS(ctest).GT.1.e-2_real_8) THEN
                ! CTM IS REAL
                ctm=t*ctm
                !$omp  parallel do private(IG,T1,TR,TI) 
                DO ig=1,nkpt%ngwk
                   t1=twnl(ig,iv,is,1)*REAL(ctm)
                   tr=REAL(auxc(ig))
                   ti=AIMAG(auxc(ig))
                   c2(ig,clsd%ialpha)=c2(ig,clsd%ialpha)+CMPLX(t1*tr,t1*ti,kind=real_8)
                ENDDO
             ELSE
                ! CTM IS IMAGINARY
                ctm=t*ctm
                !$omp  parallel do private(IG,T1,TR,TI) 
                DO ig=1,nkpt%ngwk
                   t1=twnl(ig,iv,is,1)*AIMAG(ctm)
                   tr=REAL(auxc(ig))
                   ti=AIMAG(auxc(ig))
                   c2(ig,clsd%ialpha)=c2(ig,clsd%ialpha)+CMPLX(-t1*ti,t1*tr,kind=real_8)
                ENDDO
             ENDIF
             ! STATE=BETA
             IF (cntl%tfdist) THEN
                CALL zeroing(ddia(1:ions0%na(is)))!,ions0%na(is))
                IF (clsd%ialpha.GE.parap%nst12(parai%mepos,1).AND.&
                     clsd%ialpha.LE.parap%nst12(parai%mepos,2)) THEN
                   iaa=clsd%ialpha-parap%nst12(parai%mepos,1)+1
                   CALL dcopy(ions0%na(is),fnl2(1,isa0+1,iv,iaa,1),1,ddia,1)
                   CALL dscal(ions0%na(is),-0.5_real_8*ca*sa,ddia,1)
                ENDIF
                IF (clsd%ibeta.GE.parap%nst12(parai%mepos,1).AND.&
                     clsd%ibeta.LE.parap%nst12(parai%mepos,2)) THEN
                   iaa=clsd%ibeta-parap%nst12(parai%mepos,1)+1
                   CALL daxpy(ions0%na(is),-sa*sa,&
                        fnl2(1,isa0+1,iv,iaa,1),1,ddia,1)
                ENDIF
                CALL mp_sum(ddia,ions0%na(is),parai%allgrp)
             ELSE
                CALL dcopy(ions0%na(is),fnl(1,isa0+1,iv,clsd%ialpha,1),1,ddia,1)
                CALL dscal(ions0%na(is),-0.5_real_8*sa*ca,ddia,1)
                CALL daxpy(ions0%na(is),-sa*sa,fnl(1,isa0+1,iv,clsd%ibeta,1),&
                     1,ddia,1)
             ENDIF
             CALL dgemv('N',2*ncpw%ngw,ions0%na(is),1.0_real_8,eigr(1,isa0+1,1),2*ncpw%ngw,&
                  ddia,1,0.0_real_8,auxc(1),1)
             t=-wsg(is,iv)
             ctm=(0.0_real_8,-1.0_real_8)**nghtol(iv,is)
             ctest=REAL(ctm)
             ! MAKE USE OF THE SPECIAL STRUCTURE OF CTM
             IF (ABS(ctest).GT.1.e-2_real_8) THEN
                ! CTM IS REAL
                ctm=t*ctm
                !$omp  parallel do private(IG,T1,TR,TI) 
                DO ig=1,nkpt%ngwk
                   t1=twnl(ig,iv,is,1)*REAL(ctm)
                   tr=REAL(auxc(ig))
                   ti=AIMAG(auxc(ig))
                   c2(ig,clsd%ibeta)=c2(ig,clsd%ibeta)+CMPLX(t1*tr,t1*ti,kind=real_8)
                ENDDO
             ELSE
                ! CTM IS IMAGINARY
                ctm=t*ctm
                !$omp  parallel do private(IG,T1,TR,TI) 
                DO ig=1,nkpt%ngwk
                   t1=twnl(ig,iv,is,1)*AIMAG(ctm)
                   tr=REAL(auxc(ig))
                   ti=AIMAG(auxc(ig))
                   c2(ig,clsd%ibeta)=c2(ig,clsd%ibeta)+CMPLX(-t1*ti,t1*tr,kind=real_8)
                ENDDO
             ENDIF
          ENDDO
          ! ==--------------------------------------------------------------==
       ENDIF
       isa0=isa0+ions0%na(is)
    ENDDO
    ! ==--------------------------------------------------------------==
  END SUBROUTINE fcasnl
  ! ==================================================================

END MODULE fnonloc_utils
