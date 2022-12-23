MODULE extrap_utils
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: extrap

CONTAINS

  ! ==================================================================
  SUBROUTINE extrap(nnr1,alex,rhom1,rho0,rinp)
    ! ==--------------------------------------------------------------==
    ! == ALEXANDER MIXING FOR MOLECULAR DYNAMICS                      ==
    ! == OR GEOMETRY OPTIMIZATION                                     ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nnr1
    REAL(real_8)                             :: alex
    REAL(real_8), DIMENSION(nnr1)            :: rhom1, rho0, rinp

! ==--------------------------------------------------------------==

    CALL daxpy(nnr1,-1._real_8,rho0(1),1,rhom1(1),1)
    CALL dcopy(nnr1,rho0(1),1,rinp(1),1)
    CALL daxpy(nnr1,-alex,rhom1(1),1,rinp(1),1)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE extrap
  ! ==================================================================

END MODULE extrap_utils
