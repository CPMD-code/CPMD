MODULE ranp_utils
  USE cotr,                            ONLY: lskcor
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: paral
  USE prng_utils,                      ONLY: repprngu
  USE system,                          ONLY: cntr

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ranp

CONTAINS

  ! ==================================================================
  SUBROUTINE ranp(tau0)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:)

    INTEGER                                  :: i, ia, iat, is

! Variables
! ==--------------------------------------------------------------==
! ==           Randomization of the atomic coordinates.           ==
! ==--------------------------------------------------------------==

    IF (paral%io_parent)&
         WRITE(6,'(A)') ' RANDOMIZATION - ATOMIC POSITIONS'
    DO i=1,3
       iat=0
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             iat=iat+1
             IF (lskcor(i,iat).NE.0)&
                  tau0(i,ia,is)=tau0(i,ia,is)+cntr%amprp*(repprngu()-0.5_real_8)
          ENDDO
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE ranp
  ! ==================================================================

END MODULE ranp_utils
