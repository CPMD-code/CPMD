MODULE mm_extrap
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! storage for wavefunction extrapolation in BOMD.
  COMPLEX(real_8), ALLOCATABLE, SAVE :: cold(:,:,:,:)
  INTEGER, SAVE :: numcold, nnow

  ! MTS high level functional quantities
  COMPLEX(real_8), ALLOCATABLE, SAVE :: cold_high(:,:,:,:)
  INTEGER, SAVE :: numcold_high, nnow_high

END MODULE mm_extrap
