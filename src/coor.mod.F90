MODULE coor
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==================================================================
  ! ==   DYNAMIC ALLOCATION OF ARRAYS RELATED TO IONIC COORDINATES  ==
  ! ==================================================================
  ! == TAU0(3,NAX,NSX) Atom coordinates                             ==
  ! == TAUP(3,NAX,NSX) Old atom coordinates                         ==
  ! == FION(3,NAX,NSX) Atom forces                                  ==
  ! == VELP(3,NAX,NSX) Atom velocities                              ==
  ! == LVELINI(0:NAX,NSX) Use to initialize velocities from input   ==
  ! ==                    file (velocitinp)                         ==
  ! ==--------------------------------------------------------------==
  REAL(real_8), ALLOCATABLE, SAVE :: tau0(:,:,:)
  REAL(real_8), ALLOCATABLE, TARGET, SAVE :: taup(:,:,:)
  REAL(real_8), ALLOCATABLE, TARGET, SAVE :: fion(:,:,:)
  REAL(real_8), ALLOCATABLE, TARGET, SAVE :: velp(:,:,:)

  ! 1D to be included into
  ! detsp.F
  ! mm_init.F
  ! 2D to be included into
  ! setsys.F
  ! velocitinp.F 
  LOGICAL, ALLOCATABLE, SAVE :: lvelini(:,:)

  ! ==================================================================

END MODULE coor
