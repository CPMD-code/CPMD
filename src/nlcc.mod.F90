MODULE nlcc
  USE kinds,                           ONLY: real_8
  USE system,                          ONLY: maxsp

  IMPLICIT NONE

  ! ==================================================================
  ! == NON LINEAR CORE CORRECTION (calculated in COPOT)             ==
  ! ==================================================================
  ! == TINLCC True if Non Linear Core Correction is used            ==
  ! ==              for a species                                   ==
  ! == tnlcc(MAXSP) True if for the species use NLCC                ==
  ! ==--------------------------------------------------------------==
  ! ==================================================================
  ! ==================================================================
  ! == RHOC: core charges per species                               ==
  ! ==--------------------------------------------------------------==
  REAL(real_8), ALLOCATABLE, SAVE :: rcgrid(:,:)
  REAL(real_8), ALLOCATABLE, SAVE :: corecg(:,:)
  REAL(real_8), ALLOCATABLE, SAVE :: rhoc(:,:)
  REAL(real_8), ALLOCATABLE, SAVE :: drhoc(:,:,:)
  REAL(real_8), ALLOCATABLE, SAVE :: rhocspl(:,:,:)

  COMPLEX(real_8), ALLOCATABLE, SAVE :: vnlcc(:,:)
  COMPLEX(real_8), ALLOCATABLE, SAVE :: vnlt(:)


  ! ==================================================================
  REAL(real_8), ALLOCATABLE, SAVE :: roct(:)

  ! ==================================================================
  TYPE :: corei_t
     INTEGER :: nlcct(maxsp)
     INTEGER :: meshcc(maxsp)
  END TYPE corei_t
  TYPE(corei_t), SAVE :: corei
  TYPE :: corel_t
     LOGICAL :: tinlc
     LOGICAL :: tnlcc(maxsp)
  END TYPE corel_t
  TYPE(corel_t), SAVE :: corel
  TYPE :: corer_t
     REAL(real_8) :: anlcc(maxsp)=0.0_real_8
     REAL(real_8) :: bnlcc(maxsp)=0.0_real_8
     REAL(real_8) :: enlcc(maxsp)=0.0_real_8
     REAL(real_8) :: clogcc(maxsp)=0.0_real_8
     REAL(real_8) :: excco=0.0_real_8
     REAL(real_8) :: egcxco=0.0_real_8
     REAL(real_8) :: egccco=0.0_real_8
  END TYPE corer_t
  TYPE(corer_t), SAVE :: corer

END MODULE nlcc
