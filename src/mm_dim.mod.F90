MODULE mm_dimmod
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! cpmd dimensions of the QM and the MM systems
  INTEGER, ALLOCATABLE :: NAq(:)

  INTEGER, ALLOCATABLE :: NAm(:)

  INTEGER, PARAMETER :: mm_go_qm=1, mm_go_mm=2, mm_revert=3 
  LOGICAL :: mm_stat         ! if true we are qm

  ! from classical order to cpmd order
  INTEGER, ALLOCATABLE :: cpat(:)

  INTEGER, ALLOCATABLE :: cpsp(:)

  INTEGER, ALLOCATABLE :: NAT_cpmd(:) ! must be allocated as (0:)
  INTEGER, ALLOCATABLE :: NAT_grm(:) ! must be allocated as (0:)

  INTEGER, ALLOCATABLE :: gratom(:)


  ! classical charges of all the atoms
  REAL(real_8), ALLOCATABLE :: mm_charge(:,:)

  INTEGER :: NCAG_l  ! number of chege groups
  INTEGER, ALLOCATABLE :: INC_l(:)

  ! classical box
  TYPE :: clsaabox_t
     REAL(real_8) :: BOX_au(3)
     REAL(real_8) :: mm_c_trans(3)
  END TYPE clsaabox_t
  TYPE(clsaabox_t) :: clsaabox

  TYPE :: mmdim_t
     INTEGER :: NSPq            ! number of quantum species
     INTEGER :: NATq            ! total number of quantum atoms
     INTEGER :: NAXq            ! maximum number of quantum atoms per species.
     INTEGER :: NSPm            ! total number of species
     INTEGER :: NATm            ! total number of atoms
  END TYPE mmdim_t
  TYPE(mmdim_t) :: mmdim

  TYPE :: solsolv_t
     INTEGER :: nrpt
     INTEGER :: nsolv
  END TYPE solsolv_t
  TYPE(solsolv_t) :: solsolv

  ! constraints
  TYPE :: solvvv_t
     INTEGER :: NCONS_gr
     INTEGER :: NRAM_gr
  END TYPE solvvv_t
  TYPE(solvvv_t) :: solvvv

END MODULE mm_dimmod
