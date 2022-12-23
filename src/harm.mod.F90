MODULE harm
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==================================================================
  INTEGER :: ngharm
  REAL(real_8), ALLOCATABLE :: freq(:)
  REAL(real_8), ALLOCATABLE :: dcgt(:)
  REAL(real_8), ALLOCATABLE :: dsgtw(:)
  REAL(real_8), ALLOCATABLE :: wdsgt(:)
  REAL(real_8), ALLOCATABLE :: dtan2w(:)
  REAL(real_8), ALLOCATABLE :: dtan2c(:)
  REAL(real_8), ALLOCATABLE :: xmu(:)


  ! ==================================================================

END MODULE harm
