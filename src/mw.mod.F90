! ==================================================================
MODULE mw
  IMPLICIT NONE

  ! == MULTIPLE WALKER METADYNAMCIS                                 ==
  ! ==================================================================
  LOGICAL :: tmw
  TYPE :: mwi_t
     INTEGER :: nwalk
     INTEGER :: walker_id
  END TYPE mwi_t
  TYPE(mwi_t), SAVE :: mwi
  ! ==--------------------------------------------------------------==

END MODULE mw
