MODULE rinitwf_utils
  USE atomwf_utils,                    ONLY: give_scr_atomwf
  USE copot_utils,                     ONLY: give_scr_copot
  USE forcedr_utils,                   ONLY: give_scr_forcedr
  USE newcell_utils,                   ONLY: give_scr_newcell
  USE ortho_utils,                     ONLY: give_scr_ortho
  USE system,                          ONLY: cnti,&
                                             cntl
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: give_scr_rinitwf

CONTAINS

  ! ==================================================================
  SUBROUTINE give_scr_randwf(lrandwf,tag,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lrandwf
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate

    CHARACTER(*), PARAMETER                  :: procedureN = 'give_scr_randwf'

    INTEGER                                  :: isub, lcopot, lforces, lortho

    CALL tiset(procedureN,isub)
    CALL give_scr_ortho(lortho,tag,nstate)
    CALL give_scr_copot(lcopot,tag)
    CALL give_scr_forcedr(lforces,tag,nstate,.TRUE.,.FALSE.)
    lrandwf=MAX(lortho,lcopot,lforces)
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_randwf
  ! ==================================================================
  SUBROUTINE give_scr_rinitwf(lrinitwf,tag,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lrinitwf
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate

    CHARACTER(*), PARAMETER :: procedureN = 'give_scr_rinitwf'

    INTEGER                                  :: isub, lnewcell, lother

    CALL tiset(procedureN,isub)
    lrinitwf=0
    lnewcell=0
    lother=0
    IF (cntl%tprcp.OR.cntl%tpres) THEN
       CALL give_scr_newcell(lnewcell,tag)
    ELSE
       lnewcell=0
    ENDIF
    IF(cnti%inwfun.EQ.1.OR.cnti%inwfun.EQ.3) THEN
       CALL give_scr_randwf(lother,tag,nstate)
    ELSEIF (cnti%inwfun.EQ.2) THEN
       CALL give_scr_atomwf(lother,tag,nstate)
    ELSE
       lother=0
    ENDIF
    lrinitwf=MAX(lnewcell,lother)
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_rinitwf
  ! ==================================================================
END MODULE rinitwf_utils
