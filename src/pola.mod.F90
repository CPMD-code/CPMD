MODULE pola
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==================================================================
  ! == Calculation of polarisability (Ali Alavi - Oct 1997)         ==
  ! ==--------------------------------------------------------------==
  INTEGER, PARAMETER :: iupola=64 
  INTEGER, PARAMETER :: itol=1,itmax=10000 
  REAL(real_8), PARAMETER :: tol=1.e-7_real_8,tollocc=1.e-1_real_8 
  ! ==--------------------------------------------------------------==
  REAL(real_8), ALLOCATABLE :: zeff(:,:,:,:)
  REAL(real_8), ALLOCATABLE :: zeffav(:,:,:,:)
  REAL(real_8) :: alphap(3,3),alphapav(3,3)

  INTEGER :: ipolarise,nfpolar
  LOGICAL :: tpolarb,tzeff
  ! ==================================================================


END MODULE pola
