MODULE ainitwf_utils
  USE atomwf_utils,                    ONLY: dsygvx
  USE atrho_utils,                     ONLY: atrho,&
                                             give_scr_atrho
  USE atwf,                            ONLY: atwp,&
                                             catom,&
                                             xsmat,&
                                             xxmat
  USE error_handling,                  ONLY: stopgm
  USE gsortho_utils,                   ONLY: gs_ortho
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE ksmat_dist_utils,                ONLY: dist_ksmat
  USE ksmat_utils,                     ONLY: give_scr_ksmat,&
                                             ksmat
  USE mp_interface,                    ONLY: mp_bcast
  USE ortho_utils,                     ONLY: ortho
  USE ovlap_utils,                     ONLY: ovlap
  USE parac,                           ONLY: parai,&
                                             paral
  USE phfac_utils,                     ONLY: phfac
  USE prng_utils,                      ONLY: repprngu_vec_cmplx
  USE pslo,                            ONLY: pslo_com
  USE rscpot_utils,                    ONLY: give_scr_rscpot,&
                                             rscpot
  USE spin,                            ONLY: lspin2,&
                                             spin_mod
  USE summat_utils,                    ONLY: give_scr_summat,&
                                             summat
  USE system,                          ONLY: cntl,&
                                             nkpt
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE wfnio_utils,                     ONLY: hasrands,&
                                             queryrands
  USE zeroing_utils,                   ONLY: zeroing
#include "sizeof.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ainitwf
  PUBLIC :: give_scr_ainitwf

CONTAINS

  ! ==================================================================
  SUBROUTINE ainitwf(c0,nstate,tau0,fion,rhoe,psi)
    ! ==--------------------------------------------------------------==
    ! == ATOMIC WAVEFUNCTIONS AS INITIAL GUESS FOR ADDITIONAL STATES  ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(nkpt%ngwk,nstate,*)
    REAL(real_8)                             :: tau0(:,:,:), fion(:,:,:), &
                                                rhoe(:,:)
    COMPLEX(real_8)                          :: psi(:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'ainitwf'

    COMPLEX(real_8), ALLOCATABLE             :: cscr(:), cx(:,:)
    INTEGER                                  :: i, ierr, ii, iopt, isub
    LOGICAL                                  :: debug, tlsd2
    REAL(real_8), ALLOCATABLE                :: eivat(:), workat(:), &
                                                xmatat(:,:)

! Variables
! ==--------------------------------------------------------------==

    tlsd2 = cntl%tlsd
    IF (nstate.EQ.0) RETURN
    IF (tkpts%tkpnt) RETURN
    IF (pslo_com%tivan) RETURN
    IF (cntl%ttau) RETURN
    IF (lspin2%tlse) RETURN
    IF (.NOT.hasrands()) RETURN
    ! ==--------------------------------------------------------------==
    CALL tiset('   AINITWF',isub)
    ! ==--------------------------------------------------------------==
    debug=.FALSE.
    ! Phase factors
    CALL phfac(tau0)
    ! Set up catom

    ALLOCATE(catom(nkpt%ngwk,MAX(nstate,atwp%nattot)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL atrho(rhoe,psi(:,1),nstate)
    ! Orthogonalize atomic wavefunctions
    IF (cntl%tdmal) THEN
       IF (tkpts%tkpnt) CALL stopgm('AINITWF','DISTRIBUTED INITIALIZATION '&
            // 'NOT AVAILABLE WITH K-POINTS',& 
            __LINE__,__FILE__)
       ALLOCATE(cscr(nkpt%ngwk*MAX(nstate,atwp%nattot)),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(cscr)!,MAX(atwp%nattot,nstate)*nkpt%ngwk)
       CALL ortho(atwp%nattot, catom(:,:), cscr(1))
       DEALLOCATE(cscr,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! local potential
    CALL rscpot(c0(:,:,1),tau0,fion,rhoe,psi,.FALSE.,.FALSE.,&
         nstate,1)
    ! ==--------------------------------------------------------------==
    ALLOCATE(cx(nkpt%ngwk,atwp%numaormax),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(cx)
    IF (.NOT.cntl%tdmal) THEN
       ALLOCATE(xxmat(atwp%nattot,atwp%nattot),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(xxmat)
       ALLOCATE(xsmat(atwp%nattot,atwp%nattot),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (cntl%tdmal) THEN
       CALL dist_ksmat(cx,rhoe,psi(:,1),nstate,1,c0(:,:,1:1),MIN(nstate,&
            atwp%nattot),tlsd2,.TRUE.)
       GOTO 100
    ENDIF
    IF (cntl%tlsd) THEN
       cntl%tlsd=.FALSE.
       ! alpha spin
       CALL ksmat(cx,rhoe(:,1:1),psi(:,1),spin_mod%nsup,1)
       CALL ovlap(atwp%nattot,xsmat,catom,catom)
       CALL summat(xxmat,atwp%nattot)
       CALL summat(xsmat,atwp%nattot)
       ALLOCATE(xmatat(atwp%nattot, atwp%nattot),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(workat(10*atwp%nattot),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(eivat(atwp%nattot),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)

       IF (paral%parent) THEN
          iopt=1
          CALL dsygvx(iopt,xxmat,atwp%nattot,xsmat,atwp%nattot,eivat,xmatat,&
               atwp%nattot,atwp%nattot,workat,SIZE(workat))
       ENDIF
       CALL mp_bcast(xmatat,SIZE(xmatat),parai%source,parai%allgrp)
       DO i=1,spin_mod%nsup
          IF (queryrands(i)) THEN
             IF(i<=atwp%nattot) THEN !vw need to protect vs out of bound
                IF (nkpt%ngwk.GT.0) THEN
                   CALL dgemv('N',2*nkpt%ngwk,atwp%nattot,1.0_real_8,catom(1,1),2*nkpt%ngwk,&
                        xmatat(1,i),1,0.0_real_8,c0(1,i,1),1)
                ENDIF
             ELSE
                CALL repprngu_vec_cmplx(nkpt%ngwk,c0(:,i,1))
             ENDIF
             CALL gs_ortho(c0(1,1,1),i-1,c0(1,i,1),1)
          ENDIF
       ENDDO
       DEALLOCATE(xmatat,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(workat,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(eivat,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)

       ! beta spin
       CALL ksmat(cx,rhoe(:,2:2),psi(:,1),nstate,1)
       CALL ovlap(atwp%nattot,xsmat,catom,catom)
       CALL summat(xxmat,atwp%nattot)
       CALL summat(xsmat,atwp%nattot)
       ALLOCATE(xmatat(atwp%nattot, atwp%nattot),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(workat(10*atwp%nattot),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(eivat(atwp%nattot),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)

       IF (paral%parent) THEN
          iopt=1
          CALL dsygvx(iopt,xxmat,atwp%nattot,xsmat,atwp%nattot,eivat,xmatat,&
               atwp%nattot,atwp%nattot,workat,SIZE(workat))
       ENDIF
       CALL mp_bcast(xmatat,SIZE(xmatat),parai%source,parai%allgrp)
       DO i=spin_mod%nsup+1,nstate
          IF (queryrands(i)) THEN
             ii=i-spin_mod%nsup
             IF(ii<=atwp%nattot) THEN !vw need to protect vs out of bound
                IF (nkpt%ngwk.GT.0) THEN
                   CALL dgemv('N',2*nkpt%ngwk,atwp%nattot,1.0_real_8,catom(1,1),2*nkpt%ngwk,&
                        xmatat(1,ii),1,0.0_real_8,c0(1,i,1),1)
                ENDIF
             ELSE
                CALL repprngu_vec_cmplx(nkpt%ngwk,c0(:,i,1))
             ENDIF
             CALL gs_ortho(c0(1,spin_mod%nsup+1,1),ii-1,c0(1,i,1),1)
          ENDIF
       ENDDO
       DEALLOCATE(xmatat,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(workat,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(eivat,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       cntl%tlsd=.TRUE.
    ELSE
       CALL ksmat(cx,rhoe,psi(:,1),nstate,1)
       CALL ovlap(atwp%nattot,xsmat,catom,catom)
       CALL summat(xxmat,atwp%nattot)
       CALL summat(xsmat,atwp%nattot)
       ALLOCATE(xmatat(atwp%nattot, atwp%nattot),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(workat(10*atwp%nattot),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(eivat(atwp%nattot),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)

       IF (paral%parent) THEN
          iopt=1
          CALL dsygvx(iopt,xxmat,atwp%nattot,xsmat,atwp%nattot,eivat,xmatat,&
               atwp%nattot,atwp%nattot,workat,SIZE(workat))
       ENDIF
       CALL mp_bcast(xmatat,SIZE(xmatat),parai%source,parai%allgrp)
       DO i=1,nstate
          IF (queryrands(i)) THEN
             IF(i<=atwp%nattot) THEN !vw need to protect vs out of bound
                CALL dgemv('N',2*nkpt%ngwk,atwp%nattot,1.0_real_8,catom(1,1),2*nkpt%ngwk, &
                     xmatat(1,i),1,0.0_real_8,c0(1,i,1),1)
             ELSE
                CALL repprngu_vec_cmplx(nkpt%ngwk,c0(:,i,1))
             ENDIF
             CALL gs_ortho(c0(1,1,1),i-1,c0(1,i,1),1)
          ENDIF
       ENDDO
       DEALLOCATE(xmatat,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(workat,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(eivat,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    ENDIF
100 CONTINUE
    ! ==--------------------------------------------------------------==
    IF (.NOT.cntl%tdmal) THEN
       DEALLOCATE(xsmat,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(xxmat,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    DEALLOCATE(catom,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(cx,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL tihalt('   AINITWF',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE ainitwf
  ! ==================================================================
  SUBROUTINE give_scr_ainitwf(lainitwf,tag,nstate)
    ! TODO: NSTATE is unused var, but is used in call from initrun.F
    ! should we keep it for uniformity?..
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: lainitwf
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate

    INTEGER                                  :: latrho, lksmat, lrscpot, &
                                                lsummat

    lainitwf=atwp%nattot*atwp%nattot+3*atwp%nattot+atwp%nattot
    CALL give_scr_atrho(latrho,tag)
    CALL give_scr_rscpot(lrscpot,tag,.FALSE.)
    CALL give_scr_ksmat(lksmat,tag)
    CALL give_scr_summat(lsummat,tag,atwp%nattot)
    lainitwf=MAX(lainitwf,lrscpot,lksmat,lsummat,latrho)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_ainitwf
  ! ==================================================================

END MODULE ainitwf_utils
