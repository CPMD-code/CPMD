MODULE tbxc
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==================================================================
  ! == EXCHANGE-CORRELATION INCLUDE FILE                            ==
  ! ==================================================================
  LOGICAL :: toldcode
  ! ==================================================================


  ! ==================================================================
  REAL(real_8), ALLOCATABLE :: eexc(:) ! TODO must be allocated as (0:)
  REAL(real_8), ALLOCATABLE :: vvxc(:) ! TODO must be allocated as (0:)
  ! ==================================================================
  TYPE :: tabx_t
     REAL(real_8) :: rmaxxc
     REAL(real_8) :: ddro
     REAL(real_8) :: ddgc
     REAL(real_8) :: rmaxbx
     REAL(real_8) :: bero
     INTEGER :: narray
  END TYPE tabx_t
  TYPE(tabx_t) :: tabx

END MODULE tbxc
