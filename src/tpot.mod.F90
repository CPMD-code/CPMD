MODULE tpot
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  COMPLEX(real_8), POINTER :: c0y(:,:,:)

  REAL(real_8), ALLOCATABLE :: eigval(:)
  REAL(real_8), ALLOCATABLE :: foccp(:)
  REAL(real_8), ALLOCATABLE :: eigref(:)


  INTEGER :: mstate

END MODULE tpot
