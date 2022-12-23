MODULE poin
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE
  ! ==================================================================
  REAL(real_8), ALLOCATABLE :: potr(:,:)

  REAL(real_8), POINTER :: rhoo(:,:)
  ! use pointer, because RHOO is used as pointer in many places, refactoring postponed

  COMPLEX(real_8), ALLOCATABLE :: cstate(:)
  REAL(real_8), ALLOCATABLE :: ptau(:,:)
  REAL(real_8), ALLOCATABLE :: rtau(:,:)

  ! ==================================================================
END MODULE poin
