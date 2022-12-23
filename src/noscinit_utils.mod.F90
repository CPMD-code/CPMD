MODULE noscinit_utils
  USE cnst,                            ONLY: factem
  USE nose,                            ONLY: etc,&
                                             etcdot,&
                                             nchc,&
                                             qnoscc
  USE system,                          ONLY: cntr
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: noscinit

CONTAINS

  ! ==================================================================
  SUBROUTINE noscinit(ip)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ip

    INTEGER                                  :: i

! Variables
! ==--------------------------------------------------------------==
! ==  CELL NOSE PARAMETERS                                        ==
! ==--------------------------------------------------------------==

    CALL zeroing(etc(:,ip))!,nchc)
    DO i=1,nchc
       etcdot(i,ip) = SQRT(cntr%tempc/(factem*qnoscc(i)))
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE noscinit
  ! ==================================================================



END MODULE noscinit_utils
