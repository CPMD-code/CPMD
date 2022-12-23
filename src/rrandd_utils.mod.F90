MODULE rrandd_utils
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: paral
  USE prng_utils,                      ONLY: repprngu
  USE spin,                            ONLY: clsd
  USE system,                          ONLY: fpar

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rrandd

CONTAINS

  ! ==================================================================
  SUBROUTINE rrandd(rho,amp)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rho(fpar%nnr1,clsd%nlsd), amp

    REAL(real_8), PARAMETER                  :: tiny = 1.e-8_real_8 

    INTEGER                                  :: i, ilsd
    REAL(real_8)                             :: dt, dt2
    REAL(real_8), EXTERNAL                   :: dasum

! ==--------------------------------------------------------------==
! ==           Randomization of the electronic density.           ==
! ==--------------------------------------------------------------==

    IF (paral%io_parent)&
         WRITE(6,'(A)') ' RANDOMIZATION - DENSITY'
    DO ilsd = 1 , clsd%nlsd
       dt = dasum(fpar%nnr1,rho(1,ilsd),1)
       DO i = 1 , fpar%nnr1
          rho(i,ilsd)=rho(i,ilsd)+amp*(repprngu()-0.5_real_8)
          IF (rho(i,ilsd).LE.0._real_8) rho(i,ilsd) = tiny
       ENDDO
       dt2=dasum(fpar%nnr1,rho(1,ilsd),1)
       CALL dscal(fpar%nnr1,dt/dt2,rho(1,ilsd),1)
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rrandd
  ! ==================================================================

END MODULE rrandd_utils
