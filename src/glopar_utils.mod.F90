MODULE glopar_utils
  USE loadpa_utils,                    ONLY: leadim
  USE system,                          ONLY: fpar,&
                                             ncpw,&
                                             nkpt,&
                                             parm,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: glopar

CONTAINS

  ! ==================================================================
  SUBROUTINE glopar
    ! ==--------------------------------------------------------------==
    CHARACTER(*), PARAMETER                  :: procedureN = 'glopar'

    INTEGER                                  :: isub

! ==--------------------------------------------------------------==
! ==  SAVE THE GLOBAL SYSTEM PARAMETERS, SO THEY CAN BE CHANGED   ==
! ==  IN THE PARALLEL VERSION OF THE CODE                         ==
! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    spar%nhgs    = ncpw%nhg
    spar%nhgks   = nkpt%nhgk
    spar%nhgls   = -1
    spar%ngws    = ncpw%ngw
    spar%ngwks   = nkpt%ngwk
    spar%ngwls   = -1
    spar%nr1s    = parm%nr1
    spar%nr2s    = parm%nr2
    spar%nr3s    = parm%nr3
    CALL leadim(spar%nr1s,spar%nr2s,spar%nr3s,fpar%kr1s,fpar%kr2s,fpar%kr3s)
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE glopar
  ! ==================================================================

END MODULE glopar_utils
