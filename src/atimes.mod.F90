MODULE atimesmod
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==================================================================
  INTEGER :: ikind_atimes
  REAL(real_8) :: atimes_eval
  TYPE :: atic_t
     REAL(real_8) :: atimes_eval
     INTEGER :: ikind_atimes
  END TYPE atic_t
  TYPE(atic_t) :: atic
  ! ==================================================================

END MODULE atimesmod
