MODULE linres
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==================================================================
  ! == INCLUDE FILE FOR LINEAR RESPONSE                             ==
  ! ==================================================================
  INTEGER :: lrsym
  ! ==================================================================
  COMPLEX(real_8), POINTER :: clrwf(:,:,:,:) ! TODO refactor CLRWF

  COMPLEX(real_8), ALLOCATABLE :: clrv(:,:,:,:)

  REAL(real_8), ALLOCATABLE :: urot(:,:,:)
  REAL(real_8), ALLOCATABLE :: wcent(:,:)

  INTEGER, ALLOCATABLE :: ireor(:)

  INTEGER :: nlinw,nlinr,nolr,nua,nub
  LOGICAL :: lractive
  ! ==================================================================
  ! ==================================================================
  ! ..TSH[
  ! ==================================================================
  CHARACTER (len=30) :: shfilen1,shfilen2
  CHARACTER (len=12) :: shfilen3

  COMPLEX(real_8), ALLOCATABLE :: ck(:,:)
  COMPLEX(real_8), ALLOCATABLE :: c1k(:,:,:)
  COMPLEX(real_8), ALLOCATABLE :: c_sh(:,:)
  COMPLEX(real_8), ALLOCATABLE :: cdotlct(:,:)

  REAL(real_8), ALLOCATABLE :: shsigma(:,:,:)
  REAL(real_8), ALLOCATABLE :: shsigmaf(:,:,:)
  REAL(real_8), ALLOCATABLE :: v_sh(:,:)

  REAL(real_8) :: dt_sh,tempp_sh
  REAL(real_8), PARAMETER :: step_diff=1.e-5_real_8 
  ! ..TSH]
  ! ==================================================================
  TYPE :: lr01_t
     INTEGER :: int1
     INTEGER :: nmaxlr
     INTEGER :: lopti
     INTEGER :: ndpoint
     INTEGER :: mldiis
  END TYPE lr01_t
  TYPE(lr01_t) :: lr01
  TYPE :: lr02_t
     REAL(real_8) :: rea1
     REAL(real_8) :: xc_eps
     REAL(real_8) :: lr_hthrs
     REAL(real_8) :: dlrs
     REAL(real_8) :: tol_lr
     REAL(real_8) :: tol_qs
     REAL(real_8) :: thauto(2)
  END TYPE lr02_t
  TYPE(lr02_t) :: lr02
  TYPE :: lr03_t
     LOGICAL :: log1
     LOGICAL :: txc_analytic
     LOGICAL :: txc_dd_ana
     LOGICAL :: tpara_gauge
     LOGICAL :: tgauge_all
  END TYPE lr03_t
  TYPE(lr03_t) :: lr03
  TYPE :: lrf1_t
     INTEGER :: td_x
     INTEGER :: td_c
     INTEGER :: td_gx
     INTEGER :: td_gc
     INTEGER :: td_hf
     INTEGER :: td_mtau
  END TYPE lrf1_t
  TYPE(lrf1_t) :: lrf1
  TYPE :: lrf2_t
     LOGICAL :: td_tgc
     LOGICAL :: td_tgcx
     LOGICAL :: td_tgcc
     LOGICAL :: td_code
     LOGICAL :: td_hybrid
     LOGICAL :: td_ttau
  END TYPE lrf2_t
  TYPE(lrf2_t) :: lrf2
  TYPE :: lrf3_t
     REAL(real_8) :: tdpxlda
     REAL(real_8) :: tdpclda
     REAL(real_8) :: tdpxgc
     REAL(real_8) :: tdpcgc
     REAL(real_8) :: tdphfx
  END TYPE lrf3_t
  TYPE(lrf3_t) :: lrf3
  TYPE :: lrf4_t
     INTEGER :: td_method
     INTEGER :: ks_functional
     INTEGER :: td_functional
     INTEGER :: kxc_functional
  END TYPE lrf4_t
  TYPE(lrf4_t) :: lrf4
  TYPE :: lrhd_t
     LOGICAL :: local_orb
     LOGICAL :: diagonal
  END TYPE lrhd_t
  TYPE(lrhd_t) :: lrhd
  TYPE :: lrhi_t
     INTEGER :: numo
     INTEGER :: refat
  END TYPE lrhi_t
  TYPE(lrhi_t) :: lrhi
  TYPE :: td01_t
     INTEGER :: int2
     INTEGER :: ns_sin
     INTEGER :: ns_tri
     INTEGER :: ns_mix
     INTEGER :: ldiag
     INTEGER :: ndavmax
     INTEGER :: ndavspace
     INTEGER :: msubta
     INTEGER :: tdpcgiter
     INTEGER :: rotit
     INTEGER :: fstate
     INTEGER :: store_fstate
     INTEGER :: nreor
     INTEGER :: zatom
     INTEGER :: npstate
     INTEGER :: ntdiismax
     INTEGER :: nrestdmax
     INTEGER :: ioutput
  END TYPE td01_t
  TYPE(td01_t) :: td01
  TYPE :: td02_t
     REAL(real_8) :: rea2
     REAL(real_8) :: epstdav
     REAL(real_8) :: tdrandom
     REAL(real_8) :: tdpcgstep
     REAL(real_8) :: epsrot
     REAL(real_8) :: rotstep
     REAL(real_8) :: rdiistin
  END TYPE td02_t
  TYPE(td02_t) :: td02
  TYPE :: td03_t
     LOGICAL :: log2
     LOGICAL :: tda
     LOGICAL :: tdacanon
     LOGICAL :: tdpcgmin
     LOGICAL :: treorder
     LOGICAL :: tdlocal
     LOGICAL :: los
     LOGICAL :: lqp
     LOGICAL :: lbf
     LOGICAL :: lbfh
     LOGICAL :: lberry
     LOGICAL :: tprop
     LOGICAL :: molstat
     LOGICAL :: tdlz
     LOGICAL :: sh_diag
  END TYPE td03_t
  TYPE(td03_t) :: td03
  TYPE :: tshf_t
     REAL(real_8) :: apot
     REAL(real_8) :: etpot
     REAL(real_8) :: aampl
     REAL(real_8) :: afreq
     REAL(real_8) :: apara1
  END TYPE tshf_t
  TYPE(tshf_t) :: tshf
  TYPE :: tshi_t
     INTEGER :: nroot_sh
     INTEGER :: shstep
     INTEGER :: nshs
     INTEGER :: adir
  END TYPE tshi_t
  TYPE(tshi_t) :: tshi
  TYPE :: tshl_t
     LOGICAL :: tdtully
     LOGICAL :: tully_sh
     LOGICAL :: txfmqc
     LOGICAL :: nacv_direct_only
     LOGICAL :: s0_sh
     LOGICAL :: tdextpot
     LOGICAL :: sh_phex
     LOGICAL :: isc
     LOGICAL :: isc_2ndlr
  END TYPE tshl_t
  TYPE(tshl_t) :: tshl
  TYPE :: shlct_t
     LOGICAL :: sh_lcontrol
     INTEGER :: adir
     INTEGER :: lct_state
     REAL(real_8) :: sh_lct_seed_t
     REAL(real_8) :: sh_lct_seed_m
     REAL(real_8) :: sh_lct_l
     REAL(real_8) :: trdipl
  END TYPE shlct_t
  TYPE(shlct_t) :: shlct
  TYPE :: xfmqc_t
     INTEGER      :: n_xftraj
     REAL(real_8),ALLOCATABLE  :: sigma(:,:,:,:) !FEDE sigma is an array
     REAL(real_8) :: threshold ! MIN: threshold for determining gaussian groups
     REAL(real_8), ALLOCATABLE :: eigv(:)
     REAL(real_8), ALLOCATABLE :: nacv(:,:,:,:)  
     !     REAL(real_8), ALLOCATABLE :: qmoment(:,:,:,:,:,:)
     REAL(real_8), ALLOCATABLE :: qmoment(:,:,:,:) ! MIN: reducing the array dimension 
     REAL(real_8), ALLOCATABLE :: tauall(:,:,:,:)
     REAL(real_8), ALLOCATABLE :: w_ij(:,:,:,:,:)
     REAL(real_8), ALLOCATABLE :: k_li(:,:)
     REAL(real_8), ALLOCATABLE :: f_nli(:,:,:,:)
     REAL(real_8), ALLOCATABLE :: bf_nli(:,:,:,:,:)
     REAL(real_8), ALLOCATABLE :: r0_ni(:,:,:,:)
     REAL(real_8), ALLOCATABLE :: d0_ni(:,:,:,:)
     REAL(real_8), ALLOCATABLE :: bw0_nkli(:,:,:,:,:) ! MIN: buffer for weights
     REAL(real_8), ALLOCATABLE :: fion_state(:,:,:,:)
     REAL(real_8), ALLOCATABLE :: bfion_state(:,:,:,:,:)
     REAL(real_8), ALLOCATABLE :: fa_ni(:,:,:)
     REAL(real_8), ALLOCATABLE :: bfa_ni(:,:,:,:)
     REAL(real_8), ALLOCATABLE :: brho_l(:,:) ! MIN: buffer for |C_l|^2
     COMPLEX(real_8), ALLOCATABLE :: cf(:)
     COMPLEX(real_8), ALLOCATABLE :: bcf(:,:)
  END TYPE xfmqc_t
  TYPE(xfmqc_t) :: xfmqc

END MODULE linres
