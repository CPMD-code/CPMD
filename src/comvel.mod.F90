MODULE comvelmod
  IMPLICIT NONE

  ! ==================================================================
  ! common block to store parameters for handling of the removal of
  ! the center of mass and rotational velocities.
  ! ==================================================================
  TYPE :: comvl_t
     LOGICAL :: tsubcom
     LOGICAL :: tsubrot
     LOGICAL :: subcom
     LOGICAL :: subrot
     INTEGER :: ncomv
     INTEGER :: nrotv
  END TYPE comvl_t
  TYPE(comvl_t) :: comvl

END MODULE comvelmod
