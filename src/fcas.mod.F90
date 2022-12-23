MODULE fcas
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==================================================================
  ! ==   DYNAMIC ALLOCATION OF ARRAYS FOR CAS22 FORCES              ==
  ! ==================================================================
  ! == FION_A(3,NAX,NSX) Atom Forces Ground State                   ==
  ! == FION_AB(3,NAX,NSX) Atom Forces Interaction Energy            ==
  ! ==--------------------------------------------------------------==
  INTEGER :: init_flag
  REAL(real_8), ALLOCATABLE :: fion_a(:,:,:)
  REAL(real_8), ALLOCATABLE :: fion_ab(:,:,:)


  ! ==================================================================

END MODULE fcas
