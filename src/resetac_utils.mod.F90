MODULE resetac_utils
  USE kinds,                           ONLY: real_8
  USE system,                          ONLY: acc,&
                                             maxsys
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: resetac

CONTAINS

  ! ==================================================================
  SUBROUTINE resetac(tau0,taui,nfi)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:), taui(:,:,:)
    INTEGER                                  :: nfi

! ==--------------------------------------------------------------==
! STATISTICAL ACCUMULATORS

    CALL zeroing(acc)!,nacc)
    ! SAVE INITIAL POSITIONS
    CALL dcopy(3*maxsys%nax*maxsys%nsx,tau0(1,1,1),1,taui(1,1,1),1)
    ! RESET THE TIME COUNTER
    nfi=0
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE resetac
  ! ==================================================================

END MODULE resetac_utils
