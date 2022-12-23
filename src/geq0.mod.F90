MODULE geq0mod
  IMPLICIT NONE

  ! ==================================================================
  ! == Used in parallel version otherwise GEQ0 is always .TRUE.     ==
  ! == GEQ0 = .TRUE. if the processors are the 0 component          ==
  ! ==================================================================
  LOGICAL :: geq0
  ! ==================================================================

END MODULE geq0mod
