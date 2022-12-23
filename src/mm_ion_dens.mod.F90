MODULE mm_ion_dens
  USE kinds,                           ONLY: real_8
  USE system,                          ONLY: maxsp

  IMPLICIT NONE

  REAL(real_8), ALLOCATABLE :: mm_RHOPS(:,:)

  REAL(real_8) :: mm_RAGGIO(maxsp)

END MODULE mm_ion_dens
