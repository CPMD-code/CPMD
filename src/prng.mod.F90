MODULE prng
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  TYPE :: prng_com_t
     REAL(real_8) :: repseed(2,3) = 0.0_real_8
     REAL(real_8) :: paraseed(2,3)  = 0.0_real_8
  END TYPE prng_com_t
  TYPE(prng_com_t), SAVE :: prng_com
  TYPE(prng_com_t), SAVE, ALLOCATABLE :: pi_prng_com(:)

END MODULE prng
