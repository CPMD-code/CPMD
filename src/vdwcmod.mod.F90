MODULE vdwcmod
  USE dftd3_api,                       ONLY: dftd3_input,&
                                             dftd3_calc
  USE kinds,                           ONLY: real_8
  USE readsr_utils,                    ONLY: input_string_len
  USE system,                          ONLY: maxsp

  IMPLICIT NONE

  ! ==================================================================
  ! ==  DYNAMIC ALLOCATION OF ARRAYS RELATED TO VAN DER WAALS STUFF ==
  ! ==================================================================
  TYPE :: vdwl_t
     LOGICAL :: dcacp
     LOGICAL :: vdwc
     LOGICAL :: vdwd
  END TYPE vdwl_t
  TYPE(vdwl_t), SAVE :: vdwl
  TYPE :: vdwi_t
     INTEGER :: nxvdw
     INTEGER :: nyvdw
     INTEGER :: nzvdw
  END TYPE vdwi_t
  TYPE(vdwi_t), SAVE :: vdwi
  TYPE :: vdwr_t
     REAL(real_8) :: evdw=0.0_real_8
     REAL(real_8) :: devdw(6)=0.0_real_8
  END TYPE vdwr_t
  TYPE(vdwr_t), SAVE :: vdwr
  ! ==================================================================
  ! == AUXILIARY VARIABLES FOR EMPIRICAL vdW CORRECTION             ==
  ! ==================================================================
  ! == NVDW : Maximum number of ion pairs                           ==
  ! ==================================================================

  INTEGER, ALLOCATABLE :: idvdw(:)
  INTEGER, ALLOCATABLE :: ivdw(:)
  INTEGER, ALLOCATABLE :: jvdw(:)

  REAL(real_8), ALLOCATABLE :: vdwst(:)
  REAL(real_8), ALLOCATABLE :: vdwrm(:)
  REAL(real_8), ALLOCATABLE :: vdwbe(:)
  ! ==================================================================
  ! ==--------------------------------------------------------------==
  TYPE :: empvdwi_t
     INTEGER :: nvdw
  END TYPE empvdwi_t
  TYPE(empvdwi_t) :: empvdwi
  TYPE :: empvdwr_t
     REAL(real_8) :: vdweps
     REAL(real_8) :: s6grim
  END TYPE empvdwr_t
  TYPE(empvdwr_t) :: empvdwr
  TYPE :: empvdwc_t
     CHARACTER(len=8) :: s6_str
     CHARACTER(len=8) :: dft_func
  END TYPE empvdwc_t
  TYPE(empvdwc_t) :: empvdwc
  ! dft-d3
  LOGICAL :: tdftd3
  LOGICAL, DIMENSION(maxsp) :: tdftd3noc = .FALSE.
  TYPE :: dftd3f_t 
     CHARACTER(len=input_string_len) :: func
     INTEGER :: version
     LOGICAL :: tz
  END TYPE dftd3f_t
  TYPE(dftd3f_t) :: dftd3f
  TYPE(dftd3_input) :: dftd3i
  TYPE(dftd3_calc) :: dftd3c
  ! ==================================================================
  ! = AUXILIARY VARIABLES FOR vdW CORRECTION USING WANNIER FUNCTIONS =
  ! ==================================================================
  ! == NATWFCX : Maximum number of ions within the spread of one WF ==
  ! == NWFCX   : Maximum number of WFs belonging to one fragment    ==
  ! == NFRAGX  : Maximum number of fragments included in the system ==
  ! ==================================================================
  ! INTEGER, PARAMETER :: natwfcx=101,nwfcx=1000,nfragx=500 
  INTEGER, PARAMETER :: natwfcx=1001
  INTEGER :: nwfcx = HUGE(0),nfragx = HUGE(0)
  INTEGER, ALLOCATABLE :: ifragdata(:),boadwf(:,:),icontfragw(:,:),&
       icontfragwi(:,:,:),iwfcref(:,:,:,:),ifragw(:,:,:),nfrags(:)

  REAL(real_8), ALLOCATABLE :: rwfc(:,:,:,:),spr(:,:,:),radfrag(:),&
       tauref(:,:,:,:),taufrag(:,:,:,:),rwann(:,:,:),swann(:,:),&
       wwfcref(:,:,:,:)

  LOGICAL, ALLOCATABLE :: twannupx(:),trwanncx(:)
  ! ==--------------------------------------------------------------==
  TYPE :: vdwwfl_t
     LOGICAL :: twannup
     LOGICAL :: trwannc
     LOGICAL :: tpinfo
     LOGICAL :: tpfrag
     LOGICAL :: tpc6
     LOGICAL :: tpforce
     LOGICAL :: treffit
     LOGICAL :: tdampda
  END TYPE vdwwfl_t
  TYPE(vdwwfl_t) :: vdwwfl
  TYPE :: vdwwfi_t
     INTEGER :: iswitchvdw
     INTEGER :: icriteri
     INTEGER :: nelpwf
     INTEGER :: multifrag
     INTEGER :: nboadwf
  END TYPE vdwwfi_t
  TYPE(vdwwfi_t) :: vdwwfi
  TYPE :: vdwwfr_t
     REAL(real_8) :: a6
     REAL(real_8) :: zlevel
     REAL(real_8) :: enmonomer
     REAL(real_8) :: xmfacwf
     REAL(real_8) :: tolwann
     REAL(real_8) :: tolref
  END TYPE vdwwfr_t
  TYPE(vdwwfr_t) :: vdwwfr

  INTEGER, DIMENSION(:,:), ALLOCATABLE :: npt12!(0:maxcpu,2)
  ! ==================================================================


END MODULE vdwcmod
