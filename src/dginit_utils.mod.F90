MODULE dginit_utils
  USE dg,                              ONLY: edgcomm,&
                                             ipooldg,&
                                             real_8,&
                                             tdgcomm
  USE fftnew_utils,                    ONLY: addfftnset
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: dg_init

CONTAINS

  ! CASPUR 02/05/2004
  ! ================================================================
  SUBROUTINE dg_init
    ! ================================================================
    ! Initialize the FFT data structure for the DOUBLEGRID.
    ! The new FFT grid is identified by the index IPOOLDG in the 
    ! array FPOOLV.
    ! ================================================================
    ! Variables
    CHARACTER(*), PARAMETER                  :: procedureN = 'dg_init'

    INTEGER                                  :: isub

! ================================================================

    CALL tiset(procedureN,isub)
    IF (tdgcomm%tdg) THEN
       CALL addfftnset(edgcomm%ecutdg,edgcomm%ecutwdg,ipooldg)
    ELSE
       CALL addfftnset(-1._real_8,-1._real_8,ipooldg)
       ipooldg=0
    ENDIF
    CALL tihalt(procedureN,isub)
    ! ================================================================
    RETURN
  END SUBROUTINE dg_init
  ! ================================================================

END MODULE dginit_utils
