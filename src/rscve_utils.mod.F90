MODULE rscve_utils
  USE geq0mod,                         ONLY: geq0
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: paral

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rscve

CONTAINS

  ! ==================================================================
  SUBROUTINE rscve(ekin1,ekin2,ekincp,ekinw,cm,nstate,ngw)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: ekin1, ekin2, ekincp, ekinw
    INTEGER                                  :: nstate, ngw
    COMPLEX(real_8)                          :: cm(ngw*nstate)

    REAL(real_8)                             :: alfae

! Variables
! ==--------------------------------------------------------------==
! ==  This is to perform an annealing on the electronic degrees   ==
! ==  of freedom. It is simply a control of the electronic        ==
! ==  "temperature", by rescaling the velocities.                 ==
! ==  WFS ARE ORTHOGONALIZED AFTER RESCALING                      ==
! ==--------------------------------------------------------------==

    IF (ekincp.GT.ekin1.OR.ekincp.LT.ekin2.AND.ekincp.NE.0._real_8) THEN
       alfae=SQRT(ekinw/ekincp)
       IF ((geq0).AND.paral%io_parent)&
            WRITE(6,'(A)') ' RESCALING ELEC. VELOCITIES '
       CALL dscal(2*nstate*ngw,alfae,cm(1),1)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rscve
  ! ==================================================================

END MODULE rscve_utils
