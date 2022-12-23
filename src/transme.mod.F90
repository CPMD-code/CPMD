! ==================================================================
MODULE transmemod
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! Include file of the CDFT transition matrix element calculation code
  ! H. Oberhofer (ho246@cam.ac.uk) 2009
  ! ==================================================================

  ! ==================================================================
  REAL(real_8), PARAMETER :: transme_tol=1.0e-3_real_8 
  ! ==================================================================
  INTEGER, ALLOCATABLE :: gindex(:,:,:)
  INTEGER :: mxhg,biga
  INTEGER, ALLOCATABLE :: minyh(:,:)

  COMPLEX(real_8), ALLOCATABLE :: fullw(:)

  COMPLEX(real_8), ALLOCATABLE :: fullp(:)

  COMPLEX(real_8), ALLOCATABLE :: w(:)

END MODULE transmemod
