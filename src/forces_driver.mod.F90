#include "cpmd_global.h"

MODULE forces_driver
  USE cp_grp_utils,                    ONLY: cp_grp_get_cp_rank,&
                                             cp_grp_get_sizes,&
                                             cp_grp_redist
  USE csize_utils,                     ONLY: csize
  USE dotp_utils,                      ONLY: dotp
  USE efield_utils,                    ONLY: extfield
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE fft_maxfft,                      ONLY: maxfft
  USE fnonloc_utils,                   ONLY: fnonloc,&
                                             give_scr_fnonloc
  USE func,                            ONLY: func1
  USE geq0mod,                         ONLY: geq0
  USE gsize_utils,                     ONLY: gsize
  USE hfx_drivers,                     ONLY: hfx
  USE hfxmod,                          ONLY: hfxc3
  USE hnlmat_utils,                    ONLY: hnlmat
  USE hubbardu,                        ONLY: c2u0,hubbu
  USE hubbardu_utils,                  ONLY: add_hubbardu 
  USE jrotation_utils,                 ONLY: set_orbdist
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE mp_interface,                    ONLY: mp_sum
  USE nlforce_utils,                   ONLY: give_scr_nlforce,&
                                             nlforce
  USE nlps,                            ONLY: imagp
  USE norm,                            ONLY: cnorm,&
                                             gemax,&
                                             gnmax,&
                                             gnorm
  USE nvtx_utils
  USE opeigr_utils,                    ONLY: give_scr_opeigr
  USE ovlap_utils,                     ONLY: ovlap,&
                                             ovlap2_dist
  USE parac,                           ONLY: parai,&
                                             paral
  USE pslo,                            ONLY: pslo_com
  USE puttau_utils,                    ONLY: taucl
  USE reshaper,                        ONLY: reshape_inplace
  USE rnlfl_utils,                     ONLY: rnlfl
  USE rnlsm_utils,                     ONLY: rnlsm
  USE ropt,                            ONLY: ropt_mod
  USE rotate_utils,                    ONLY: rotate
  USE rscpot_utils,                    ONLY: rscpot
  USE rswfmod,                         ONLY: rsactive
  USE spin,                            ONLY: clsd,&
                                             lspin2,&
                                             lspin3,&
                                             spin_mod
  USE summat_utils,                    ONLY: give_scr_summat,&
                                             summat
  USE symtrz_utils,                    ONLY: symvec
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             maxsys,&
                                             ncpw,&
                                             nkpt,&
                                             parap,&
                                             paraw
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: zclean,&
                                             zclean_k
  USE vpsi_utils,                      ONLY: vpsi
!!use rotate_utils, only : rotate_c
!!use ovlap_utils, only : ovlap_c
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: forces
  PUBLIC :: gscal
  PUBLIC :: gscal_c

CONTAINS

  ! ==================================================================
  SUBROUTINE forces(c0,c2,tau0,fion,rhoe,psi,nstate,nkpoint,lproj,tfor)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==                     THE TOTAL ENERGY                         ==
    ! ==                  THE ELECTRONIC FORCES                       ==
    ! ==                  THE FORCES ON THE IONS                      ==
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8)                          :: c0(:,:,:), c2(:,:)
    REAL(real_8)                             :: tau0(:,:,:), fion(:,:,:), &
                                                rhoe(:,:)
    COMPLEX(real_8)                          :: psi(:,:)
    INTEGER                                  :: nstate, nkpoint
    LOGICAL                                  :: lproj, tfor

    CHARACTER(*), PARAMETER                  :: procedureN = 'forces'

    CHARACTER(len=30)                        :: tag
    COMPLEX(real_8)                          :: zee
    COMPLEX(real_8), ALLOCATABLE             :: psiab(:)
    COMPLEX(real_8), ALLOCATABLE, TARGET     :: auxc(:)
    COMPLEX(real_8), EXTERNAL                :: zdotc
    COMPLEX(real_8), POINTER                 :: cgam(:)
    INTEGER :: first_g, i, ierr, ik, ikind, il_auxc, il_ddia, il_fsc, il_gam, &
      il_psiab, il_scrdip, ipp, isub, isub2, isub3, j, jj, last_g, lsummat, &
      naa, ngw_l, nstx, NSTX_grp
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: NWA12_grp
    LOGICAL                                  :: debug, redist_c2
    REAL(real_8)                             :: ee, ehfx, fi, fj, vhfx
    REAL(real_8), ALLOCATABLE                :: a1mat(:,:), a2mat(:,:), &
                                                a3mat(:,:), ddia(:), fsc(:), &
                                                scrdip(:)
    REAL(real_8), ALLOCATABLE, TARGET        :: gam(:,:)
    REAL(real_8), EXTERNAL                   :: dasum
    REAL(real_8), POINTER                    :: aux(:)

    CALL tiset(procedureN,isub)
    CALL tiset(procedureN//'_a',isub2)
    gnmax=0.0_real_8
    gnorm=0.0_real_8
    ehfx =0.0_real_8
    vhfx =0.0_real_8
    ! ==--------------------------------------------------------------==
    redist_c2=.FALSE.
    IF (func1%mhfx.EQ.1.AND..NOT.hfxc3%use_new_hfx) redist_c2=.TRUE.
    IF (pslo_com%tivan.OR.lspin2%tlse) redist_c2=.TRUE.
    IF (cntl%ttau) redist_c2=.TRUE.
    ! ==--------------------------------------------------------------==
    debug=.FALSE.
    IF (pslo_com%tivan .AND. lproj .AND. cnti%iproj.NE.0) THEN
       CALL give_scr_nlforce(il_gam,il_auxc,il_ddia,nstate)
    ELSE
       CALL give_scr_fnonloc(il_auxc,il_ddia,nstate)
       il_gam = imagp*nstate*nstate
    ENDIF
    IF (cntl%tfield) CALL give_scr_opeigr(il_scrdip,tag,nstate)
    IF (cntl%tdmal) THEN
       ALLOCATE(NWA12_grp(0:parai%cp_nproc-1,2),stat=ierr)
       IF (ierr.NE.0) CALL stopgm(procedureN,'Allocation problem',& 
            __LINE__,__FILE__)
       CALL set_orbdist(nstate,cnti%nstblk,parai%cp_nproc,NSTX_grp)
       NWA12_grp(0:parai%cp_nproc-1,1:2)=paraw%nwa12(0:parai%cp_nproc-1,1:2)
       CALL set_orbdist(nstate,cnti%nstblk,parai%nproc,nstx)
       il_gam=1
    ELSE
       CALL give_scr_summat(lsummat,tag,nstate)
       il_auxc=MAX(il_auxc,lsummat)
       il_auxc=MAX(il_auxc,nstate*nstate)! AUXC space for OVLAP (Laio A.)
    ENDIF
    il_fsc=nstate
    il_psiab=2
    IF (lspin2%tlse) il_psiab=2*maxfft
    ! ==--------------------------------------------------------------==
    ! ==   CALCULATE THE POTENTIAL AND THE FORCE ON THE IONS          ==
    ! ==--------------------------------------------------------------==
    DO ikind=1,nkpoint
       IF (nstate >0 ) THEN
          CALL rnlsm(c0(:,:,ikind),nstate,1,ikind,tfor)
       ENDIF
    ENDDO
    rsactive = cntl%krwfn
    CALL rscpot(c0(:,:,1),tau0,fion,rhoe,psi,tfor,ropt_mod%calste,nstate,nkpoint)
    CALL tihalt(procedureN//'_a',isub2)
    ! ==--------------------------------------------------------------==
    ! ==   CALCULATE THE ELECTRONIC FORCE                             ==
    ! ==--------------------------------------------------------------==
    IF (nstate.EQ.0) THEN
       gemax=0._real_8
       cnorm=0._real_8
       GOTO 2000
    ENDIF
    CALL zeroing(c2)!,nkpt%ngwk*nstate)
    ! EHR[
    ! cntl%tmdeh only needs forces acting on the nuclei
    ! skip calculation of C2 (forces on KS orbitals)
    IF (cntl%tmdeh) GOTO 2000
    ! EHR]
    ! ==--------------------------------------------------------------==
    ! == COMPUTE THE FORCE ON THE ELECTRONIC DEGREES OF FREEDOM DUE   ==
    ! == TO THE LOCAL POTENTIAL (STORED IN RHOE)                      ==
    ! ==--------------------------------------------------------------==
    IF (cntl%tdmal) THEN
       ALLOCATE(a1mat(nstate, nstx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(a2mat(nstate, nstx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(a3mat(nstate,nstx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
    ENDIF

    ALLOCATE(auxc(il_auxc),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL reshape_inplace(auxc, (/2*il_auxc/), aux)

    ALLOCATE(gam(nstate, nstate*2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL reshape_inplace(gam, (/nstate*nstate/), cgam)

    ALLOCATE(ddia(il_ddia),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(fsc(il_fsc),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(psiab(il_psiab),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    IF (cntl%tfield) THEN
       ALLOCATE(scrdip(il_scrdip),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
    ENDIF

    DO ik=1,nkpoint
       CALL vpsi(c0(:,:,ik),c2,crge%f(:,1),rhoe,psi(:,1),nstate,ik,clsd%nlsd,redist_c2)
       ! c2u0 is calculated in uprho or rscpot
       IF (hubbu%debug) THEN
            IF (paral%io_parent) write(6,*) procedureN,"| starting add_hubbardu"
       ENDIF
       IF (cntl%thubb) CALL add_hubbardu(c2,c2u0,nstate)
       CALL hfx(c0(:,:,ik),c2,crge%f(:,1),psi(:,1),nstate,ehfx,vhfx,redist_c2)
       ! McB    IF (cntl%tfield) CALL EFIELD(TAU0,FION,C0,C2,SCRDIP,IL_SCRDIP,NSTATE)
       IF (cntl%tfield) THEN
          CALL extfield(tau0,fion,c0,c2,nstate)
       ENDIF
       ener_com%exc=ener_com%exc+ehfx
       ener_com%etot=ener_com%etot+ehfx
       CALL tiset(procedureN//'_b',isub3)
       __NVTX_TIMER_START ( procedureN//'_b' )

       IF (pslo_com%tivan) THEN
          IF (tkpts%tkpnt) CALL stopgm(procedureN,'K-POINTS NOT IMPLEMENTED',& 
               __LINE__,__FILE__)
          IF (lproj.AND.cnti%iproj.NE.0) THEN
             CALL ovlap(nstate,gam,c2,c0(:,:,ik))
             CALL hnlmat(gam,crge%f,nstate)
             CALL mp_sum(gam,nstate*nstate,parai%allgrp)
             IF (tfor) CALL rnlfl(fion,gam,nstate,nkpoint)
             IF (cnti%iproj.EQ.1) THEN
                DO i=1,nstate
                   CALL daxpy(2*nkpt%ngwk,-gam(i,i),c0(1,i,ik),1,c2(1,1),1)
                ENDDO
             ELSEIF (cnti%iproj.EQ.2) THEN
                CALL rotate(-1.0_real_8,c0(:,:,ik),1.0_real_8,c2,gam,&
                     nstate,2*nkpt%ngwk,cntl%tlsd,spin_mod%nsup,spin_mod%nsdown)
             ENDIF
             CALL nlforce(c2,crge%f,gam,auxc,ddia,nstate)
          ELSE
             CALL fnonloc(c2,crge%f,nstate,ik,clsd%nlsd,.TRUE.)
          ENDIF
       ELSEIF (lspin2%tlse) THEN
          CALL fnonloc(c2,crge%f,nstate,ik,clsd%nlsd,.TRUE.)
          CALL ovlap(nstate,aux,c0(:,:,ik),c2)
          IF (lspin2%tros) THEN
             CALL addup2(nstate,gam,aux)
          ELSE
             CALL addup(nstate,gam,aux)
          ENDIF
          IF (lspin2%tmgoe.OR.lspin2%tros) THEN
             CALL mp_sum(gam,nstate*nstate,parai%allgrp)
          ELSE
             CALL summat(gam,nstate)
          ENDIF
          lspin3%hablse=gam(clsd%ialpha,clsd%ibeta)
          lspin3%haolse=dasum(clsd%ialpha-1,gam(1,clsd%ialpha),1)
          lspin3%hbolse=dasum(clsd%ialpha-1,gam(1,clsd%ibeta),1)
          lspin3%rotab=0.5_real_8*ATAN(2._real_8*lspin3%hablse/&
               (gam(clsd%ialpha,clsd%ialpha)-gam(clsd%ibeta,clsd%ibeta)))
          IF (lproj)&
               CALL rotate(-0.5_real_8,c0(:,:,ik),1._real_8,c2,gam,&
               nstate,2*nkpt%ngwk,.FALSE.,spin_mod%nsup,spin_mod%nsdown)
       ELSE
          CALL fnonloc(c2,crge%f,nstate,ik,clsd%nlsd,redist_c2)
          IF (.NOT.redist_c2) CALL cp_grp_redist(c2,nkpt%ngwk,nstate)
          IF (lproj) THEN
             DO i=1,nstate
                IF (crge%f(i,ik).EQ.0._real_8) THEN
                   fsc(i)=0._real_8
                ELSE
                   fsc(i)=1._real_8
                ENDIF
             ENDDO
             IF (cnti%iproj.EQ.0) THEN
             ELSEIF (cnti%iproj.EQ.1) THEN
                IF (tkpts%tkpnt) THEN
                   IF (nkpt%ngwk.GT.0) THEN
                      DO i=1,nstate
                         zee=fsc(i)*zdotc(nkpt%ngwk,c0(1,i,ik),1,c2(1,1),1)
                         CALL mp_sum(zee,parai%allgrp)
                         CALL zaxpy(nkpt%ngwk,-zee,c0(1,i,ik),1,c2(1,1),1)
                      ENDDO
                   ENDIF
                ELSE
                   IF (ncpw%ngw.GT.0) THEN
                      DO i=1,nstate
                         ee=fsc(i)*dotp(ncpw%ngw,c0(:,i,ik),c2(:,1))
                         CALL mp_sum(ee,parai%allgrp)
                         CALL daxpy(2*ncpw%ngw,-ee,c0(1,i,ik),1,c2(1,1),1)
                      ENDDO
                   ENDIF
                ENDIF
             ELSEIF (cnti%iproj.EQ.2) THEN
                IF (tkpts%tkpnt) THEN
                   IF (cntl%tlsd) THEN
                      naa=spin_mod%nsup+1
                      CALL zeroing(gam)!,2*nstate*nstate)
                      CALL ovlap_c(spin_mod%nsup,gam,c2,c0(:,:,ik))
                      CALL ovlap_c(spin_mod%nsdown,gam(2*naa,naa),&
                           c2(1,naa),c0(1,naa,ik))
                   ELSE
                      CALL ovlap_c(nstate,gam,c2,c0(:,:,ik))
                   ENDIF
                   ! AK FIXME: can this work? GAM is complex.
                   CALL mp_sum(gam,2*nstate*nstate,parai%allgrp)
                   CALL gscal_c(nstate,cgam,fsc)
                   CALL rotate_c(CMPLX(-1.0_real_8,0._real_8,kind=real_8),c0(:,:,ik),&
                        CMPLX( 1.0_real_8,0._real_8,kind=real_8),c2,gam,nstate)
                ELSE
                   IF (cntl%tdmal) THEN

                      CALL ovlap2_dist( nstate, a1mat, c2, c0(:,:,ik) )

                      ipp = cp_grp_get_cp_rank(parai%me,parai%cp_inter_me)
                      IF (parai%cp_me.NE.ipp) CALL stopgm(procedureN,&
                           ' something wrong here',& 
                           __LINE__,__FILE__)
                      DO i=1,nstate
                         fi= fsc(i)
                         DO j=NWA12_grp(parai%cp_me,1),nwa12_grp(parai%cp_me,2)
                            jj=j-NWA12_grp(parai%cp_me,1)+1
                            fj= fsc(j)
                            a1mat(i,jj)=(fi*fj)*a1mat(i,jj)
                            ! >>>>
                            ! SDF the overlap above is not taking care of the 
                            ! SDF orthogonality of the spins. So we do it here...
                            ! SDF some speedup could be gained in ovlap2 and rotate_da.
                            IF (cntl%tlsd) THEN
                               IF (j.GT.spin_mod%nsup.AND.i.LE.spin_mod%nsup) a1mat(i,jj)=0.0_real_8
                               IF (j.LE.spin_mod%nsup.AND.i.GT.spin_mod%nsup) a1mat(i,jj)=0.0_real_8
                            ENDIF
                            ! <<<
                         ENDDO
                      ENDDO

                      CALL cp_grp_get_sizes(ngw_l=ngw_l,first_g=first_g,&
                           last_g=last_g)

                      CALL rotate_da(-1._real_8,c0(first_g,1,ik),1._real_8,&
                           c2(first_g,1),a1mat,2*ncpw%ngw,2*ngw_l,&
                           nstate,NWA12_grp(0,1),nwa12_grp(0,2),&
                           NSTX_grp,parai%cp_me,parap%pgroup,parai%cp_nproc,parai%cp_grp,&
                           cntl%tlsd,spin_mod%nsup,spin_mod%nsdown)

                      CALL cp_grp_zero_g(c2,nkpt%ngwk,nstate,first_g,last_g)
                      CALL cp_grp_redist(c2,nkpt%ngwk,nstate)

                   ELSE
                      CALL ovlap(nstate,gam,c2,c0(:,:,ik))
                      CALL mp_sum(gam,nstate*nstate,parai%allgrp)
                      CALL gscal(nstate,gam,fsc)
                      CALL rotate(-1.0_real_8,c0(:,:,ik),1.0_real_8,c2,gam,&
                           nstate,2*nkpt%ngwk,cntl%tlsd,spin_mod%nsup,spin_mod%nsdown)
                   ENDIF
                ENDIF
             ENDIF
          ENDIF
       ENDIF
       IF (geq0) THEN
          IF (tkpts%tkpnt) THEN
             CALL zclean_k(c2,nstate,ncpw%ngw)
          ELSE
             CALL zclean(c2,nstate,ncpw%ngw)
          ENDIF
       ENDIF
       CALL csize(c2,nstate,gemax,cnorm)

       __NVTX_TIMER_STOP
       CALL tihalt(procedureN//'_b',isub3)
    ENDDO

    IF (cntl%tdmal) THEN
       DEALLOCATE(a1mat,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(a2mat,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(a3mat,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    ENDIF
    DEALLOCATE(auxc,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(gam,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(ddia,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(fsc,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(psiab,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    IF (cntl%tfield) THEN
       DEALLOCATE(scrdip,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    ENDIF

2000 CONTINUE
    rsactive = .FALSE.
    IF (tfor) THEN
       CALL mp_sum(fion,3*maxsys%nax*maxsys%nsx,parai%allgrp)
       IF (paral%parent) THEN
          CALL symvec(fion)
          CALL taucl(fion)
          CALL gsize(fion,gnmax,gnorm)
       ENDIF
    ENDIF
    ! 
    IF (ALLOCATED(NWA12_grp)) THEN
       DEALLOCATE(NWA12_grp,stat=ierr)
       IF (ierr.NE.0) CALL stopgm(procedureN,'Deallocation problem',& 
            __LINE__,__FILE__)
    ENDIF
    ! 
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==

  END SUBROUTINE forces
  ! ==================================================================
  SUBROUTINE addup(nstate,a1,a2)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: a1(nstate,*), a2(nstate,*)

    INTEGER                                  :: i, j

! ==--------------------------------------------------------------==

    DO i=1,nstate
       DO j=i,nstate
          a1(i,j)=a2(i,j)+a2(j,i)
          a1(j,i)=a1(i,j)
       ENDDO
    ENDDO
    IF (lspin2%tmgoe) THEN
       a1(clsd%ibeta,clsd%ialpha) = 2._real_8*((1._real_8-lspin3%mgmab)*a2(clsd%ibeta,clsd%ialpha)&
            +lspin3%mgmab *A2(clsd%ialpha,clsd%ibeta))
       a1(clsd%ialpha,clsd%ibeta) = 2._real_8*((1._real_8-lspin3%mgmba)*a2(clsd%ialpha,clsd%ibeta)&
            +lspin3%mgmba *A2(clsd%ibeta,clsd%ialpha))
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE addup
  ! ==================================================================
  SUBROUTINE addup2(nstate,a1,a2)
    ! NEW VERSION: MODIFIED GOEDECKER
    INTEGER                                  :: nstate
    REAL(real_8)                             :: a1(nstate,*), a2(nstate,*)

    INTEGER                                  :: i, j

! ==--------------------------------------------------------------==

    DO i=1,nstate-2
       a1(i,i) = 2.0_real_8 * a2(i,i)
       DO j=i+1,nstate-2
          a1(i,j) = a2(i,j) + a2(j,i)
          a1(j,i) = a2(j,i) + a2(i,j)
       ENDDO
       a1(nstate-1,i) = 2.0_real_8*(lspin3%mgab(1)*a2(nstate-1,i)&
            +lspin3%mgab(2)*A2(I,NSTATE-1))
       a1(nstate,i) = 2.0_real_8*(lspin3%mgab(3)*a2(nstate,i)&
            +lspin3%mgab(4)*A2(I,NSTATE))
       a1(i,nstate-1) = 2.0_real_8*(lspin3%mgab(5)*a2(i,nstate-1)&
            +lspin3%mgab(6)*A2(NSTATE-1,I))
       a1(i,nstate) = 2.0_real_8*(lspin3%mgab(7)*a2(i,nstate)&
            +lspin3%mgab(8)*A2(NSTATE,I))
    ENDDO
    a1(nstate,nstate-1) = 2.0_real_8*(lspin3%mgab(9)*a2(nstate,nstate-1)&
         +lspin3%mgab(10)*A2(NSTATE-1,NSTATE))
    a1(nstate-1,nstate) = 2.0_real_8*(lspin3%mgab(11)*a2(nstate-1,nstate)&
         +lspin3%mgab(12)*A2(NSTATE,NSTATE-1))
    a1(nstate-1,nstate-1) = 2.0_real_8 * a2(nstate-1,nstate-1)
    a1(nstate,nstate) = 2.0_real_8 * a2(nstate,nstate)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE addup2
  ! ==================================================================
  SUBROUTINE gscal(nstate,gam,f)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: gam(nstate,nstate), f(nstate)

    INTEGER                                  :: i, j

! ==--------------------------------------------------------------==

    !$omp parallel do private(I,J)
    DO j=1,nstate
       DO i=1,nstate
          gam(i,j)=f(i)*f(j)*gam(i,j)
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gscal
  ! ==================================================================
  SUBROUTINE gscal_c(nstate,gam,f)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: gam(nstate,nstate)
    REAL(real_8)                             :: f(nstate)

    INTEGER                                  :: i, j

! ==--------------------------------------------------------------==

    !$omp parallel do private(I,J)
    DO j=1,nstate
       DO i=1,nstate
          gam(i,j)=f(i)*f(j)*gam(i,j)
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gscal_c
  ! ==================================================================

END MODULE forces_driver
