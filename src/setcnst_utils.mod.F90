MODULE setcnst_utils
  USE atoms_utils,                     ONLY: atoms
  USE cnst,                            ONLY: unit_txt
  USE parac,                           ONLY: paral
  USE soft,                            ONLY: soft_com

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: setcnst

CONTAINS

  ! ==================================================================
  SUBROUTINE setcnst
    ! ==--------------------------------------------------------------==
    ! == DEFINE CONSTANTS                                             ==
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent)&
         WRITE(6,*) 'SETCNST| USING: ',unit_txt,' UNITS'
    ! ==--------------------------------------------------------------==
    soft_com%exsoft=.FALSE.
    soft_com%exnomore=.FALSE.
    ! ==--------------------------------------------------------------==
    CALL atoms
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE setcnst
  ! ==================================================================

END MODULE setcnst_utils
