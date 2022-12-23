MODULE prden
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==================================================================
  INTEGER, ALLOCATABLE :: mwfn(:)
  INTEGER :: numpr


  ! ==================================================================
  TYPE :: elfcb_t
     REAL(real_8) :: elfcut
     REAL(real_8) :: elfeps
     LOGICAL :: telf
  END TYPE elfcb_t
  TYPE(elfcb_t) :: elfcb

END MODULE prden
