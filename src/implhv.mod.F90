MODULE implhv
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==-----------------------------------------------------------------==
  ! quantities needed for implicit Hessian.vector product
  ! sd0(nodim,nodim) = part of the hessian independent on c1's
  ! rs_v(nodim,*)    = vectors generated with the product (1st random)
  ! xma(nodim)       = vector of the inverses of the sqrt of masses
  ! rhoe0(maxfft)    = array for unperturbed electronic density
  ! ==-----------------------------------------------------------------==
  REAL(real_8), ALLOCATABLE :: sd0(:,:)
  REAL(real_8), ALLOCATABLE :: rs_v(:,:)
  REAL(real_8), ALLOCATABLE :: xma(:)

END MODULE implhv
