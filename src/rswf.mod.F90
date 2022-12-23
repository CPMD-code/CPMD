MODULE rswfmod
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==================================================================
  INTEGER :: lwdim,maxstates
  COMPLEX(real_8), ALLOCATABLE, TARGET :: rswf(:,:)

  LOGICAL :: rsactive
  ! ==================================================================
  INTEGER :: lwdimx,maxstatesx
  COMPLEX(real_8), POINTER :: rswfx(:,:)

  ! ==================================================================

END MODULE rswfmod
