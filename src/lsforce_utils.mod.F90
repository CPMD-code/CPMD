MODULE lsforce_utils
  USE bsym,                            ONLY: cnstwgt
  USE kinds,                           ONLY: real_8
  USE system,                          ONLY: maxsys

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: lsforce

CONTAINS

  ! ==================================================================
  SUBROUTINE lsforce(fnbs,fion)
    ! CB: Calculates projected LS forces in BS runs      
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: fnbs(3*maxsys%nax*maxsys%nsx),&
                                                fion(3*maxsys%nax*maxsys%nsx)

    INTEGER                                  :: i

    DO i=1,3*maxsys%nax*maxsys%nsx
       fion(i) = (1+cnstwgt) * fnbs(i) - cnstwgt * fion(i)
    ENDDO
    ! ==--------------------------------------------------------------==
  END SUBROUTINE lsforce


END MODULE lsforce_utils
