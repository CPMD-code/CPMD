MODULE dipomod
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  TYPE :: moment_t
     REAL(real_8) :: dmom(3)
     REAL(real_8) :: dnuc(3)
     REAL(real_8) :: qmom(3,3)
  END TYPE moment_t
  TYPE(moment_t) :: moment

END MODULE dipomod
