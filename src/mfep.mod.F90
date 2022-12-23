MODULE mfep
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==================================================================
  REAL(real_8) :: factor,alpha_pmin,tolds

  TYPE :: mfepl_t
     LOGICAL :: projout
     LOGICAL :: rlate_string
  END TYPE mfepl_t
  TYPE(mfepl_t) :: mfepl

  INTEGER, ALLOCATABLE :: irest(:)
  LOGICAL, ALLOCATABLE :: lcvsp(:)
  ! ==================================================================
  TYPE :: mfepi_t
     INTEGER :: nloop
     INTEGER :: nequi
     INTEGER :: nloopold
     INTEGER :: ncvsp
     INTEGER :: istring
  END TYPE mfepi_t
  TYPE(mfepi_t) :: mfepi
  ! ==================================================================
END MODULE mfep
