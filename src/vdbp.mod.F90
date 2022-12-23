MODULE vdbp
  USE kinds,                           ONLY: real_8
  USE system,                          ONLY: maxsp

  IMPLICIT NONE

  ! ==================================================================
  ! == VANDERBILT PSEUDOPOTENTALS                                   ==
  ! ==================================================================
  REAL(real_8), ALLOCATABLE :: rscore(:,:)
  REAL(real_8), ALLOCATABLE :: dion(:,:,:)
  REAL(real_8), ALLOCATABLE :: betar(:,:,:)
  REAL(real_8), ALLOCATABLE :: qqq(:,:,:)
  REAL(real_8), ALLOCATABLE :: qfunc(:,:,:,:)
  REAL(real_8), ALLOCATABLE :: qrl(:,:,:,:,:)
  REAL(real_8), ALLOCATABLE :: r(:,:)
  REAL(real_8), ALLOCATABLE :: rucore(:,:,:)
  REAL(real_8), ALLOCATABLE :: ru(:,:)
  REAL(real_8), ALLOCATABLE :: rab(:,:)
  REAL(real_8), ALLOCATABLE :: rsatom(:,:)
  REAL(real_8), ALLOCATABLE :: vdb_pawf(:,:,:)
  REAL(real_8), ALLOCATABLE :: vdb_r(:,:)

  ! ==================================================================


  ! ==================================================================
  TYPE :: ncpr1_t
     REAL(real_8) :: cmesh(maxsp)
     INTEGER :: nbeta(maxsp)
     INTEGER :: nvales(maxsp)
     INTEGER :: ifpcor(maxsp)
     INTEGER :: kkbeta(maxsp)
     INTEGER :: meshva(maxsp)
  END TYPE ncpr1_t
  TYPE(ncpr1_t) :: ncpr1

END MODULE vdbp
