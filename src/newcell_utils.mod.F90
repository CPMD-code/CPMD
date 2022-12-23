MODULE newcell_utils
  USE initclust_utils,                 ONLY: gf_periodic
  USE nlccset_utils,                   ONLY: nlccset
  USE rggen_utils,                     ONLY: gvector
  USE rinforce_utils,                  ONLY: give_scr_putwnl,&
                                             putps,&
                                             putwnl,&
                                             testspline

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: newcell
  PUBLIC :: give_scr_newcell

CONTAINS

  ! ==================================================================
  SUBROUTINE newcell
    ! ==--------------------------------------------------------------==
    ! ==           INITIALIZE ALL ARRAYS THAT DEPEND ON HT            ==
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    CALL gvector
    CALL testspline
    CALL putps
    CALL putwnl
    CALL nlccset
    CALL gf_periodic
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE newcell
  ! ==================================================================
  SUBROUTINE give_scr_newcell(lnewcell,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lnewcell
    CHARACTER(len=30)                        :: tag

! ==--------------------------------------------------------------==

    CALL give_scr_putwnl(lnewcell,tag)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_newcell
  ! ==================================================================


END MODULE newcell_utils
