MODULE soft
  IMPLICIT NONE

  ! ==================================================================
  ! == INCLUDE FILE FOR SOFT EXIT                                   ==
  ! ==================================================================

  ! ==================================================================
  TYPE :: soft_com_t
     LOGICAL :: exsoft
     LOGICAL :: exnomore
  END TYPE soft_com_t
  TYPE(soft_com_t) :: soft_com

END MODULE soft
