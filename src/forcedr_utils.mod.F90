MODULE forcedr_utils
  USE forces_utils,                    ONLY: give_scr_forces
  USE k_forces_utils,                  ONLY: give_scr_kforces
  USE kpts,                            ONLY: tkpts
  USE mm_dim_utils,                    ONLY: mm_dim
  USE mm_dimmod,                       ONLY: mm_go_mm,&
                                             mm_revert
  USE noforce_utils,                   ONLY: give_scr_noforce
  USE symtrz_utils,                    ONLY: give_scr_symvec
  USE system,                          ONLY: cntl,&
                                             maxsys
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: give_scr_forcedr

CONTAINS

  ! ==================================================================
  SUBROUTINE give_scr_forcedr(lforcedr,tag,nstate,lproj,tfor)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lforcedr
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate
    LOGICAL                                  :: lproj, tfor

    CHARACTER(*), PARAMETER :: procedureN = 'give_scr_forcedr'

    INTEGER                                  :: il_auxc, il_ddia, il_gam, &
                                                il_smat, isub, lsymvec, ltscr
    LOGICAL                                  :: oldstatus

! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    CALL mm_dim(mm_go_mm,oldstatus)
    ! ivano 4/3/2011
    lforcedr=0
    ! ivano 4/3/2011
    IF (cntl%nonort) THEN
       CALL give_scr_noforce(lforcedr,il_gam,il_auxc,il_smat,il_ddia,&
            tag,nstate,tfor)
    ELSE
       ! EHR[
       IF (tkpts%tkpnt.AND.(.NOT.(cntl%tmdeh.OR.cntl%tpspec.OR.cntl%tpdist))) THEN
          CALL give_scr_kforces(lforcedr,tag,nstate,lproj,tfor)
       ELSE
          CALL give_scr_forces(lforcedr,tag,nstate,lproj,tfor)
       ENDIF
    ENDIF
    ! EHR]
    ltscr=3*maxsys%nax*maxsys%nsx
    lsymvec=0
    CALL give_scr_symvec(lsymvec,tag)
    lforcedr=MAX(lforcedr,ltscr,lsymvec)
    tag='MAX(LFORCEDR,LTSCR,LSYMVEC)'
    CALL mm_dim(mm_revert,oldstatus)
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_forcedr
  ! ==================================================================


END MODULE forcedr_utils
