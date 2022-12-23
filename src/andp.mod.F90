MODULE andp
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==================================================================
  ! == DENSITY  AND ALEXANDER MIXING POINTERS INCLUDE FILE          ==
  ! ==================================================================
  ! == DENSITY ARRAYS                                               ==
  ! ==--------------------------------------------------------------==
  REAL(real_8), ALLOCATABLE, TARGET :: rin0(:,:)
  REAL(real_8), ALLOCATABLE :: rout0(:,:)
  REAL(real_8), ALLOCATABLE :: rmix(:,:)
  ! ==================================================================

END MODULE andp
