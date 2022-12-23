MODULE pcgrad_utils
  USE forcedr_utils,                   ONLY: give_scr_forcedr
  USE mm_dim_utils,                    ONLY: mm_dim
  USE mm_dimmod,                       ONLY: mm_go_mm,&
                                             mm_revert
  USE mm_input,                        ONLY: lqmmm
  USE ortho_utils,                     ONLY: give_scr_ortho
  USE rhoofr_utils,                    ONLY: give_scr_rhoofr
  USE rnlsm_utils,                     ONLY: give_scr_rnlsm
  USE rscpot_utils,                    ONLY: give_scr_rscpot
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             maxsys

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: give_scr_pcgrad
  !public :: give_scr_linesr
  !public :: give_scr_xetot

CONTAINS

  ! ==================================================================
  SUBROUTINE give_scr_pcgrad(lpcgrad,tag,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lpcgrad
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate

    LOGICAL                                  :: oldstatus

! ==--------------------------------------------------------------==

    IF (lqmmm%qmmm) CALL mm_dim(mm_go_mm,oldstatus)
    CALL give_scr_linesr(lpcgrad,tag,nstate)
    lpcgrad=lpcgrad+3*maxsys%nax*maxsys%nsx
    IF (lqmmm%qmmm) CALL mm_dim(mm_revert,oldstatus)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_pcgrad
  ! ==================================================================
  SUBROUTINE give_scr_linesr(llinesr,tag,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: llinesr
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate

    INTEGER                                  :: lforcedr, lrhoofr

! ==--------------------------------------------------------------==

    CALL give_scr_xetot(llinesr,tag,nstate)
    CALL give_scr_forcedr(lforcedr,tag,nstate,.TRUE.,.FALSE.)
    CALL give_scr_rhoofr(lrhoofr,tag)
    llinesr=MAX(llinesr,lforcedr,lrhoofr)
    IF (lqmmm%qmmm)THEN
       llinesr=MAX(llinesr,fpar%kr1*fpar%kr2s*fpar%kr3s)
       llinesr=MAX(llinesr,maxsys%nax*maxsys%nsx*3)
    ENDIF

    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_linesr
  ! ==================================================================
  SUBROUTINE give_scr_xetot(lxetot,tag,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lxetot
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate

    INTEGER                                  :: lortho, lrnlsm, lrscpot

    IF (cntl%nonort.AND.(.NOT.cntl%quenchb)) THEN
       lortho=0
       lrnlsm=0
       lrscpot=0
    ELSE
       CALL give_scr_ortho(lortho,tag,nstate)
       CALL give_scr_rnlsm(lrnlsm,tag,nstate,.FALSE.)
       CALL give_scr_rscpot(lrscpot,tag,.FALSE.)
    ENDIF
    lxetot=MAX(lortho,lrnlsm,lrscpot)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_xetot
  ! ==================================================================

END MODULE pcgrad_utils
