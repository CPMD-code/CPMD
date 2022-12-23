MODULE cdftmod
  USE kinds,                           ONLY: real_8
  USE system,                          ONLY: maxsp

  IMPLICIT NONE
  ! ==================================================================
  ! CDFT main include file containing important CDFT variables
  ! H. Oberhofer (ho246@cam.ac.uk) 2009
  ! ==================================================================

  ! ==================================================================
  INTEGER, PARAMETER :: cdft_mat=250 ! Maximum CDFT Donor/Acceptor group size
  INTEGER, PARAMETER :: maxsavev=50 ! Maximum order of V predictor in MD
  ! ------------------------------------------------------------------
  ! Parameter Integers
  TYPE :: cdftpi_t
     INTEGER :: cdft_a(cdft_mat)
     INTEGER :: cdft_d(cdft_mat)
     INTEGER :: spa(cdft_mat)
     INTEGER :: spd(cdft_mat)
     INTEGER :: sna(cdft_mat)
     INTEGER :: snd(cdft_mat)
     INTEGER :: mmx_a
     INTEGER :: mmx_d
     INTEGER :: ndon
     INTEGER :: naccr
  END TYPE cdftpi_t
  TYPE(cdftpi_t) :: cdftpi
  ! ------------------------------------------------------------------
  TYPE :: cdftci_t
     INTEGER :: cdft_end
     INTEGER :: wstep!Control Integers
  END TYPE cdftci_t
  TYPE(cdftci_t) :: cdftci
  ! ------------------------------------------------------------------
  TYPE :: cdftlog_t
     LOGICAL :: recw       ! Recalculate the weights?
     LOGICAL :: newgrad    ! Optimisation options
     LOGICAL :: dekk       ! Optimisation options
     LOGICAL :: tczones    ! Convergence zones
     LOGICAL :: reswf      ! Reset Wavefunction after each CDFT step
     LOGICAL :: ccor       ! Cutoff Correction in MD
     LOGICAL :: rcdft      ! try to load initial V from CDFT_RESTART file
     LOGICAL :: thda       ! Calculate the transition matrix
     LOGICAL :: thdaproj   ! Calculate the transition matrix
     LOGICAL :: tauto      ! Automatically choose NC from initial guess
     LOGICAL :: tpred      ! Try to predict next V in an MD simulation
     LOGICAL :: tphio      ! Write out the overlap matrices
     LOGICAL :: tspinc     ! Constrain the spin density
     LOGICAL :: tcall      ! Constrain the density and the spin density
     LOGICAL :: tvmirr     ! if true always apply CDFT_V = - CDFT_V for second state
  END TYPE cdftlog_t
  TYPE(cdftlog_t) :: cdftlog
  ! ------------------------------------------------------------------
  TYPE :: cdftrun_t
     LOGICAL :: reconvv
  END TYPE cdftrun_t
  TYPE(cdftrun_t) :: cdftrun
  ! ------------------------------------------------------------------
  TYPE :: cdftcom_t
     REAL(real_8) :: cdft_v(2) = 0.0_real_8
     REAL(real_8) :: oldv(2)
     REAL(real_8) :: vgrad(2) = 0.0_real_8
     REAL(real_8) :: vgrado(2)
     REAL(real_8) :: vgrada
     REAL(real_8) :: oldva
     REAL(real_8) :: wslice
     REAL(real_8) :: cdft_nc = 0.0_real_8
     REAL(real_8) :: nother
     REAL(real_8) :: cdft_ns = 0.0_real_8
     REAL(real_8) :: nsother
     REAL(real_8) :: vccon
     REAL(real_8) :: vcconu
  END TYPE cdftcom_t
  TYPE(cdftcom_t), SAVE :: cdftcom
  TYPE :: cdfthess_t
     REAL(real_8) :: chess(2,2)
  END TYPE cdfthess_t
  TYPE(cdfthess_t) :: cdfthess
  TYPE :: cdftvgrs_t
     REAL(real_8) :: cdft_v2(2)
     REAL(real_8) :: vgstep
     REAL(real_8) :: maxvmov
  END TYPE cdftvgrs_t
  TYPE(cdftvgrs_t) :: cdftvgrs
  TYPE :: wcfix_t
     REAL(real_8) :: wcut
  END TYPE wcfix_t
  TYPE(wcfix_t) :: wcfix
  ! ------------------------------------------------------------------
  REAL(real_8), ALLOCATABLE :: wa(:)
  REAL(real_8), ALLOCATABLE :: wd(:)
  REAL(real_8), ALLOCATABLE :: wdiff(:)
  REAL(real_8), ALLOCATABLE :: rhol(:)
  REAL(real_8), ALLOCATABLE :: wder(:,:,:)

  REAL(real_8) :: czones(3,2)
  REAL(real_8) :: finalchrg(2)=0.0_real_8
  ! ------------------------------------------------------------------
  TYPE :: cdftmd_t
     REAL(real_8) :: cdft_rc(maxsp)
     REAL(real_8) :: rhocut(maxsp)
     INTEGER :: cdft_mmax(maxsp)
     REAL(real_8) :: cdft_shell
  END TYPE cdftmd_t
  TYPE(cdftmd_t) :: cdftmd
  ! ------------------------------------------------------------------
  TYPE :: cdftpred_t
     REAL(real_8) :: vpredbuf(maxsavev)
     INTEGER :: predord          ! Order of the extrapolation scheme
  END TYPE cdftpred_t
  TYPE(cdftpred_t) :: cdftpred
  ! ------------------------------------------------------------------
  TYPE :: cdfthda_t
     LOGICAL :: hdafirst
     LOGICAL :: hdaresb
     REAL(real_8) :: vbuff(2)
  END TYPE cdfthda_t
  TYPE(cdfthda_t) :: cdfthda
  ! ------------------------------------------------------------------
  TYPE :: sccomm_t
     INTEGER :: n_s0
     INTEGER :: n_s1
     INTEGER :: n_s0up
     INTEGER :: n_s1up
     INTEGER :: tsysk
     INTEGER :: tsysl
  END TYPE sccomm_t
  TYPE(sccomm_t) :: sccomm
  TYPE :: projlog_t
     LOGICAL :: projnow
  END TYPE projlog_t
  TYPE(projlog_t) :: projlog
  ! ------------------------------------------------------------------
  TYPE :: wgaussl_t
     LOGICAL :: twgauss! Construct the weight from Gaussians
     LOGICAL :: thdas ! HDA calculation with Single Acceptor and no Donor
     LOGICAL :: thdawm ! HDA calculation with two distinct difference weights
  END TYPE wgaussl_t
  TYPE(wgaussl_t) :: wgaussl
  INTEGER :: wg_n
  REAL(real_8) :: wg_sigma(maxsp) = 0.0_real_8
  ! ------------------------------------------------------------------
  REAL(real_8), ALLOCATABLE :: atchg2(:)
  LOGICAL :: chgset
  ! ------------------------------------------------------------------
  TYPE :: dummycom_t
     REAL(real_8) :: dondum(3)
     REAL(real_8) :: accdum(3)
  END TYPE dummycom_t
  TYPE(dummycom_t) :: dummycom
  ! ------------------------------------------------------------------
  REAL(real_8) :: cm_dr
  INTEGER :: cm_dir
  ! ==================================================================
END MODULE cdftmod
