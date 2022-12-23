MODULE noseinit_utils
  USE kinds,                           ONLY: real_8
  USE nose,                            ONLY: eta,&
                                             etadot,&
                                             nche,&
                                             nedof,&
                                             qnosee
  USE parac,                           ONLY: paral
  USE system,                          ONLY: cntr
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: noseinit

CONTAINS

  ! ==================================================================
  SUBROUTINE noseinit(ip)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ip

    INTEGER                                  :: i

! Variables
! ==--------------------------------------------------------------==
! ==  ELECTRON NOSE PARAMETERS                                    ==
! ==--------------------------------------------------------------==

    IF (paral%parent) THEN
       WRITE(6,*) 'NOSEINIT| INITIALIZATION OF NOSE VELOCITIES'
    ENDIF
    CALL zeroing(eta(:,ip))!,nche)
    etadot(1,ip) = SQRT(2._real_8*cntr%ekinw/qnosee(1))
    DO i=2,nche
       etadot(i,ip) = SQRT(2._real_8*cntr%ekinw/qnosee(i)/REAL(nedof,kind=real_8))
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE noseinit
  ! ==================================================================



END MODULE noseinit_utils
