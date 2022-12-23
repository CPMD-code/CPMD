MODULE cp_xc_utils
  USE cp_gga_correlation_utils,        ONLY: &
       cp_gga_c_check, cp_gga_c_param, cp_gga_c_hcth, cp_gga_c_lyp, cp_gga_c_p86, &
       cp_gga_c_pbe, cp_gga_c_pbe_flex, cp_gga_c_pbe_sol, &
       cp_spin_gga_c_hcth, cp_spin_gga_c_lyp, cp_spin_gga_c_p86, cp_spin_gga_c_pbe, &
       cp_spin_gga_c_pbe_flex, cp_spin_gga_c_pbe_sol
  USE cp_gga_exchange_utils,           ONLY: &
       cp_gga_x_check, cp_gga_x_param, cp_gga_screening_cam, cp_gga_x_b88, cp_gga_x_optx, &
       cp_gga_x_hcth, cp_gga_x_pbe, cp_gga_x_pbe_flex, cp_gga_x_pbe_sol, cp_gga_x_revpbe, &
       cp_spin_gga_x_b88, cp_spin_gga_x_hcth, cp_spin_gga_x_optx, cp_spin_gga_x_pbe, &
       cp_spin_gga_x_pbe_flex, cp_spin_gga_x_pbe_sol, cp_spin_gga_x_revpbe
  USE cp_dgga_exchange_utils,          ONLY: cp_dgga_x_b88,&
                                             cp_dgga_x_pbe,&
                                             cp_spin_dgga_x_pbe,&
                                             cp_spin_dgga_x_b88
  USE cp_dgga_correlation_utils,       ONLY: cp_dgga_c_lyp,&
                                             cp_dgga_c_p86,&
                                             cp_dgga_c_pbe,&
                                             cp_spin_dgga_c_lyp,&
                                             cp_spin_dgga_c_p86,&
                                             cp_spin_dgga_c_pbe
  USE cp_lda_correlation_utils,        ONLY: cp_lda_c_check,&
                                             cp_lda_c_ob_pw,&
                                             cp_lda_c_pw,&
                                             cp_lda_c_pz,&
                                             cp_lda_c_vwn,&
                                             cp_spin_lda_c_ob_pw,&
                                             cp_spin_lda_c_pw,&
                                             cp_spin_lda_c_pz,&
                                             cp_spin_lda_c_vwn
  USE cp_lda_exchange_utils,           ONLY: cp_lda_x_check,&
                                             cp_lda_x,&
                                             cp_spin_lda_x
  USE cp_mgga_exchange_utils,          ONLY: cp_mgga_x_check,&
                                             cp_mgga_x_param,&
                                             cp_mgga_x_tpss,&
                                             cp_mgga_x_m05_m06,&
                                             cp_mgga_x_m08_m11,&
                                             cp_mgga_x_mn12_mn15,&
                                             ! cp_mgga_x_vs98, &
                                             ! cp_spin_mgga_x_vs98, &
                                             cp_spin_mgga_x_tpss,&
                                             cp_spin_mgga_x_m05_m06,&
                                             cp_spin_mgga_x_m08_m11,&
                                             cp_spin_mgga_x_mn12_mn15
  USE cp_mgga_correlation_utils,       ONLY: cp_mgga_c_check,&
                                             cp_mgga_c_param,&
                                             cp_mgga_c_tpss,&
                                             cp_mgga_c_m05_m06,&
                                             cp_mgga_c_m08_m11,&
                                             ! cp_mgga_c_vs98, &
                                             ! cp_spin_mgga_c_vs98, &
                                             cp_spin_mgga_c_tpss,&
                                             cp_spin_mgga_c_m05_m06, &
                                             cp_spin_mgga_c_m08_m11
  USE cpfunc_types,                    ONLY: cp_xc_functional_t,&
                                             cp_xc_scratch_t,&
                                             cp_dxc_functional_t,&
                                             cp_dxc_scratch_t,&
                                             cp_xc_spin_components,&
                                             cp_xc_spin_pairs,&
                                             cp_xc_abs_xyz,&
                                             cp_xc_xyz
  USE error_handling,                  ONLY: stopgm
  USE func,                            ONLY: mgcx_is_skipped, &
                                             mgcx_is_becke88, &
                                             mgcx_is_pbex, &
                                             mgcx_is_revpbex, &
                                             mgcx_is_optx, &
                                             mgcc_is_skipped, &
                                             mgcc_is_lyp, &
                                             mgcc_is_perdew86, &
                                             mgcc_is_pbec, &
                                             msrx_is_CAM,&
                                             msrx_is_ashcroft,&
                                             msrx_is_erfc,&
                                             msrx_is_exp
  USE kinds,                           ONLY: default_string_length,&
                                             real_8
  USE lxc_utils,                       ONLY: XC_POLARIZED,&
                                             XC_UNPOLARIZED,&
                                             lxc_func_get_info,&
                                             lxc_func_init,&
                                             lxc_functional_get_number,&
                                             xc_f03_func_info_t,&
                                             xc_f03_func_t


  USE, INTRINSIC :: ISO_C_BINDING,  ONLY: C_INT

  IMPLICIT NONE
  PRIVATE

  INTEGER, PARAMETER, PUBLIC       :: cp_xc_max_nbr_funcs = 32
  INTEGER, PARAMETER, PRIVATE      :: a                   =  1
  INTEGER, PARAMETER, PRIVATE      :: b                   =  2
  INTEGER, PARAMETER, PRIVATE      :: ab                  =  3

  TYPE, PUBLIC :: cp_xc_env_t
     LOGICAL :: init        = .FALSE.
     LOGICAL :: xc_driver   = .FALSE.
     LOGICAL :: libxc       = .FALSE.
     LOGICAL :: polarized   = .FALSE.
     LOGICAL :: get_CAM_GGA = .FALSE.
     !
     CHARACTER(default_string_length) :: set_name       = 'EXCHANGE-CORRELATION FUNCTIONALS'
     !
     LOGICAL                          :: get_hfx        = .FALSE.
     LOGICAL                          :: overwrite_hfx  = .FALSE.
     CHARACTER(default_string_length) :: hfx_operator   = ' '
     REAL(real_8)                     :: hfx_screening  = 0.0_real_8
     REAL(real_8)                     :: hfx_constant   = 0.0_real_8
     REAL(real_8)                     :: hfx_attenuated = 0.0_real_8
     !
     LOGICAL :: use_compatibility_mode = .FALSE.
     !
     INTEGER :: n_funcs   = 0
     CHARACTER(default_string_length), DIMENSION(cp_xc_max_nbr_funcs) :: funcs=''
     REAL(real_8), DIMENSION(cp_xc_max_nbr_funcs)                     :: scales = 1.0_real_8
     LOGICAL, DIMENSION(cp_xc_max_nbr_funcs)                          :: via_libxc = .FALSE.
  CONTAINS
     PROCEDURE, PASS, PUBLIC :: report              => cp_xc_env_report
     PROCEDURE, PASS, PUBLIC :: get_pseudo_code     => cp_xc_env_get_pseudo_code
     PROCEDURE, PASS, PUBLIC :: set_hybrid_defaults => cp_xc_env_set_hybrid_defaults
  END TYPE cp_xc_env_t

  !
  ! Libxc
  !
  TYPE, PRIVATE :: cp_lxc_t
     INTEGER, PUBLIC :: n_funcs = 0
     INTEGER( C_INT ), DIMENSION(cp_xc_max_nbr_funcs), PUBLIC         :: funcs_id = 0_C_INT
     TYPE(xc_f03_func_t), DIMENSION(cp_xc_max_nbr_funcs), PUBLIC      :: funcs
     TYPE(xc_f03_func_info_t), DIMENSION(cp_xc_max_nbr_funcs), PUBLIC :: infos
     REAL(real_8), DIMENSION(cp_xc_max_nbr_funcs), PUBLIC             :: scales  = 1.0_real_8
  END TYPE cp_lxc_t

  !
  ! HFX parameters
  !
  TYPE, PRIVATE :: cp_hfx_t
     LOGICAL, PUBLIC      :: is_coulomb_attenuated ! CAM
     LOGICAL, PUBLIC      :: is_screened_erfc      ! HSE etc.
     LOGICAL, PUBLIC      :: is_screened_exp       ! work by Hutter and co.
     LOGICAL, PUBLIC      :: is_screened_ashcroft  ! (obvious)
     REAL(real_8), PUBLIC :: scale = 0.0_real_8
     REAL(real_8), PUBLIC :: alpha = 0.0_real_8 
     REAL(real_8), PUBLIC :: beta  = 0.0_real_8
     REAL(real_8), PUBLIC :: gamma = 0.0_real_8
  CONTAINS
     PROCEDURE, PASS, PUBLIC :: init => hfx_func_init
  END TYPE cp_hfx_t

  !
  ! Internal functionals
  !
  TYPE, PRIVATE :: cp_func_t
     INTEGER, PUBLIC :: n_funcs = 0
     LOGICAL, PUBLIC :: is_coulomb_attenuated = .FALSE.
     LOGICAL, PUBLIC :: use_compatibility_mode = .FALSE.
     REAL(real_8), DIMENSION(cp_xc_max_nbr_funcs), PUBLIC :: scales = 1.0_real_8
     CHARACTER(default_string_length), DIMENSION(cp_xc_max_nbr_funcs), PUBLIC :: funcs=''
  CONTAINS
     PROCEDURE, PASS, PUBLIC :: init => cpfunc_func_init
  END TYPE cp_func_t

  TYPE, PUBLIC :: cp_xc_t
     LOGICAL, PRIVATE :: init = .FALSE.

     LOGICAL, PUBLIC  :: use_cpfunc = .FALSE.
     LOGICAL, PUBLIC  :: use_libxc  = .FALSE.
     LOGICAL, PUBLIC  :: polarized  = .FALSE.

     LOGICAL, PRIVATE :: is_hybrid = .FALSE.
     LOGICAL, PRIVATE :: is_mgga   = .FALSE.
     LOGICAL, PRIVATE :: is_gga    = .FALSE.

     TYPE(cp_hfx_t), PUBLIC  :: hfx
     TYPE(cp_lxc_t), PUBLIC  :: libxc
     TYPE(cp_func_t), PUBLIC :: cpfunc
  CONTAINS
     PROCEDURE, PASS, PUBLIC   :: create                       => cp_xc_create
     PROCEDURE, PASS, PUBLIC   :: get                          => cp_xc_get
     PROCEDURE, PASS, PUBLIC   :: set_functional               => cp_xc_set_functional
     PROCEDURE, PASS, PUBLIC   :: set_functional_derivatives   => cp_xc_set_functional_derivatives
     PROCEDURE, NOPASS, PUBLIC :: purge_functional             => cp_xc_purge_functional
     PROCEDURE, NOPASS, PUBLIC :: purge_functional_derivatives => cp_xc_purge_functional_derivatives
     PROCEDURE, NOPASS, PUBLIC :: set_scratch                  => cp_xc_set_scratch
     PROCEDURE, NOPASS, PUBLIC :: set_scratch_derivatives      => cp_xc_set_scratch_derivatives
     PROCEDURE, NOPASS, PUBLIC :: purge_scratch                => cp_xc_purge_scratch
     PROCEDURE, NOPASS, PUBLIC :: purge_scratch_derivatives    => cp_xc_purge_scratch_derivatives
  END TYPE cp_xc_t

  ! ==================================================================
  !
  ! XC functional pointers for internal functionals and enhancment factors
  ! K:   Spin-restricted
  ! Ks:  Spin-unrestricted
  ! dK:  Analytical linres derivatives, spin-restricted
  ! dKs: Analytical linres derivatives, spin-unrestricted
  !
  TYPE, PRIVATE :: cp_xc_K_p_t
     PROCEDURE(cp_xc_functional_p), POINTER, NOPASS       :: compute   => NULL()
  END TYPE cp_xc_K_p_t
  !
  TYPE, PRIVATE :: cp_dxc_K_p_t
     PROCEDURE(cp_dxc_functional_p), POINTER, NOPASS      :: compute   => NULL()
  END TYPE cp_dxc_K_p_t
  !
  TYPE, PRIVATE :: cp_xc_Ks_p_t
     PROCEDURE(cp_xc_spin_functional_p), POINTER, NOPASS  :: compute   => NULL()
  END TYPE cp_xc_Ks_p_t
  !
  TYPE, PRIVATE :: cp_dxc_Ks_p_t
     PROCEDURE(cp_dxc_spin_functional_p), POINTER, NOPASS :: compute   => NULL()
  END TYPE cp_dxc_Ks_p_t
  !
  TYPE, PUBLIC :: cp_xc_functional_p_t
     TYPE(cp_xc_K_p_t),  DIMENSION(cp_xc_max_nbr_funcs)   :: K
     TYPE(cp_xc_Ks_p_t), DIMENSION(cp_xc_max_nbr_funcs)   :: Ks
     PROCEDURE(cp_xc_functional_p), POINTER, NOPASS       :: attenuate => NULL()
  END TYPE cp_xc_functional_p_t
  !
  TYPE, PUBLIC :: cp_dxc_functional_p_t
     TYPE(cp_dxc_K_p_t),  DIMENSION(cp_xc_max_nbr_funcs)  :: dK
     TYPE(cp_dxc_Ks_p_t), DIMENSION(cp_xc_max_nbr_funcs)  :: dKs
  END TYPE cp_dxc_functional_p_t
  !
  INTERFACE
     PURE SUBROUTINE cp_xc_functional_p( scratch,functional )
       IMPORT                                     cp_xc_scratch_t
       IMPORT                                     cp_xc_functional_t
       TYPE(cp_xc_scratch_t), INTENT(in)       :: scratch
       TYPE(cp_xc_functional_t), INTENT(inout) :: functional
     END SUBROUTINE cp_xc_functional_p
     PURE SUBROUTINE cp_xc_spin_functional_p( scratch,functional )
       IMPORT                                     cp_xc_scratch_t
       IMPORT                                     cp_xc_functional_t
       IMPORT                                     cp_xc_spin_components
       TYPE(cp_xc_scratch_t), DIMENSION(cp_xc_spin_components), &
            INTENT(in)    :: scratch
       TYPE(cp_xc_functional_t), DIMENSION(cp_xc_spin_components), &
            INTENT(inout) :: functional
     END SUBROUTINE cp_xc_spin_functional_p
     !
     PURE SUBROUTINE cp_dxc_functional_p( scratch,functional )
       IMPORT                                      cp_dxc_scratch_t
       IMPORT                                      cp_dxc_functional_t
       TYPE(cp_dxc_scratch_t), INTENT(in)       :: scratch
       TYPE(cp_dxc_functional_t), INTENT(inout) :: functional
     END SUBROUTINE cp_dxc_functional_p
     PURE SUBROUTINE cp_dxc_spin_functional_p( scratch,functional )
       IMPORT                                      cp_dxc_scratch_t
       IMPORT                                      cp_dxc_functional_t
       IMPORT                                      cp_xc_spin_pairs
       TYPE(cp_dxc_scratch_t), DIMENSION(cp_xc_spin_pairs), &
            INTENT(in)    :: scratch
       TYPE(cp_dxc_functional_t), DIMENSION(cp_xc_spin_pairs), &
            INTENT(inout) :: functional
     END SUBROUTINE cp_dxc_spin_functional_p
  END INTERFACE

  ! ==================================================================
  !
  ! Scratch space for internal functionals
  !
  TYPE, PUBLIC :: cp_xc_scratch_p_t
     PROCEDURE(cp_xc_scratch_p), POINTER, NOPASS       :: restricted   => NULL()
     PROCEDURE(cp_xc_spin_scratch_p), POINTER, NOPASS  :: unrestricted => NULL()
  END TYPE cp_xc_scratch_p_t
  !
  TYPE, PUBLIC :: cp_dxc_scratch_p_t
     PROCEDURE(cp_dxc_scratch_p), POINTER, NOPASS      :: restricted   => NULL()
     PROCEDURE(cp_dxc_spin_scratch_p), POINTER, NOPASS :: unrestricted => NULL()
  END TYPE cp_dxc_scratch_p_t
  !
  INTERFACE
     PURE SUBROUTINE cp_xc_scratch_p(scratch,n_in,grad_squared_in,tau_in)
       IMPORT                                real_8
       IMPORT                                cp_xc_scratch_t
       REAL(real_8), INTENT(in)           :: n_in, grad_squared_in, tau_in
       TYPE(cp_xc_scratch_t), INTENT(out) :: scratch
     END SUBROUTINE cp_xc_scratch_p
     PURE SUBROUTINE cp_xc_spin_scratch_p(scratch,n_in,grad_squared_in,tau_in)
       IMPORT                                real_8
       IMPORT                                cp_xc_scratch_t
       IMPORT                                cp_xc_spin_components
       IMPORT                                cp_xc_spin_pairs
       REAL(real_8), DIMENSION(cp_xc_spin_components), &
            INTENT(in)           :: n_in, grad_squared_in
       REAL(real_8), DIMENSION(cp_xc_spin_pairs), &
            INTENT(in)           :: tau_in
       TYPE(cp_xc_scratch_t), DIMENSION(cp_xc_spin_components), &
            INTENT(out) :: scratch
     END SUBROUTINE cp_xc_spin_scratch_p
     !
     PURE SUBROUTINE cp_dxc_scratch_p(scratch,n_in,dn_in,grad_squared_in,dgrad_in)
       IMPORT                                 real_8
       IMPORT                                 cp_dxc_scratch_t
       IMPORT                                 cp_xc_xyz
       IMPORT                                 cp_xc_abs_xyz
       REAL(real_8), INTENT(in)            :: n_in, dn_in
       REAL(real_8), DIMENSION(cp_xc_xyz), &
                     INTENT(in)            :: dgrad_in
       REAL(real_8), DIMENSION(cp_xc_abs_xyz), &
                     INTENT(in)            :: grad_squared_in
       TYPE(cp_dxc_scratch_t), INTENT(out) :: scratch
     END SUBROUTINE cp_dxc_scratch_p
     PURE SUBROUTINE cp_dxc_spin_scratch_p(scratch,n_in,dn_in,grad_squared_in,dgrad_in)
       IMPORT                                real_8
       IMPORT                                cp_dxc_scratch_t
       IMPORT                                cp_xc_spin_pairs
       IMPORT                                cp_xc_xyz
       IMPORT                                cp_xc_abs_xyz
       REAL(real_8), DIMENSION(cp_xc_spin_pairs), &
                     INTENT(in)  :: n_in, dn_in
       REAL(real_8), DIMENSION(cp_xc_spin_pairs*cp_xc_xyz), &
                     INTENT(in)  :: dgrad_in
       REAL(real_8), DIMENSION(cp_xc_spin_pairs*cp_xc_abs_xyz), &
                     INTENT(in)  :: grad_squared_in
       TYPE(cp_dxc_scratch_t), DIMENSION(cp_xc_spin_pairs), &
                     INTENT(out) :: scratch
     END SUBROUTINE cp_dxc_spin_scratch_p
  END INTERFACE


  TYPE(cp_xc_env_t), SAVE, PUBLIC :: cp_xc_functional_env
  TYPE(cp_xc_env_t), SAVE, PUBLIC :: cp_xc_kernel_env
  TYPE(cp_xc_env_t), SAVE, PUBLIC :: cp_xc_mts_low_func_env

  TYPE(cp_xc_t), SAVE, TARGET,  PUBLIC :: cp_xc_functional
  TYPE(cp_xc_t), SAVE, TARGET,  PUBLIC :: cp_xc_mts_low_func
  TYPE(cp_xc_t), SAVE, TARGET,  PUBLIC :: cp_xc_kernel
  TYPE(cp_xc_t), SAVE, POINTER, PUBLIC :: cp_xc => NULL()

  PRIVATE :: cp_xc_env_report
  PRIVATE :: cp_xc_env_get_pseudo_code
  PRIVATE :: cp_xc_create
  PRIVATE :: cp_xc_get
  PRIVATE :: cp_xc_set_functional
  PRIVATE :: cp_xc_purge_functional
  PRIVATE :: cp_xc_set_functional_derivatives
  PRIVATE :: cp_xc_purge_functional_derivatives
  PRIVATE :: cp_xc_set_scratch
  PRIVATE :: cp_xc_purge_scratch
  PRIVATE :: cp_xc_set_scratch_derivatives
  PRIVATE :: cp_xc_purge_scratch_derivatives

  PRIVATE :: get_lda_scratch
  PRIVATE :: get_gga_scratch
  PRIVATE :: get_meta_scratch
  PRIVATE :: get_dgga_scratch

  PRIVATE :: cpfunc_func_init

CONTAINS

  ! ==================================================================
  ! Initialisation of the cp_xc type
  ! Splitting of input into cpfunc and libxc types
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE cp_xc_create ( this, cp_xc_env )
    CLASS(cp_xc_t), INTENT(INOUT)            :: this
    TYPE(cp_xc_env_t), INTENT(IN)            :: cp_xc_env

    CHARACTER(*), PARAMETER                  :: procedureN = 'cp_xc_create'

    CHARACTER(len=default_string_length)     :: functional_tag
    INTEGER                                  :: i_func, i_func_libxc
    LOGICAL                                  :: compatibility_mode, &
                                                gga_coulomb_attenuated, is_hybrid, &
                                                is_polarized

    IF( this%init ) CALL stopgm(procedureN,'already initialized',&
         __LINE__,__FILE__)

    is_polarized   = cp_xc_env%polarized
    this%polarized = is_polarized

    this%use_cpfunc = .FALSE.
    this%use_libxc  = .FALSE.

    this%libxc%n_funcs  = 0
    this%cpfunc%n_funcs = 0

    !
    ! Since cp_xc_driver is called in gcener, the GGA flag has to be .true.,
    ! even if one uses only LDA...
    !
    this%is_gga = .TRUE.
    IF (ANY( INDEX(cp_xc_env%funcs(:),"MGGA_") /= 0 )) this%is_mgga = .TRUE.

    !
    ! Init generic info
    !
    is_hybrid                          = cp_xc_env%get_hfx
    this%is_hybrid                     = is_hybrid

    gga_coulomb_attenuated             = cp_xc_env%get_CAM_GGA
    this%cpfunc%is_coulomb_attenuated  = gga_coulomb_attenuated
    IF (gga_coulomb_attenuated .AND. this%is_mgga) CALL stopgm(procedureN, &
       'CAM and MGGA are incompatible', __LINE__, __FILE__)

    compatibility_mode                 = cp_xc_env%use_compatibility_mode
    this%cpfunc%use_compatibility_mode = compatibility_mode

    !
    ! Init xc functional
    !
    i_func_libxc = 0
    DO i_func = 1, cp_xc_env%n_funcs
       IF (cp_xc_env%via_libxc( i_func )) THEN
          !
          i_func_libxc   = i_func_libxc + 1
          !
          functional_tag = cp_xc_env%funcs( i_func )
          this%libxc%scales( i_func_libxc ) = cp_xc_env%scales( i_func )
          this%libxc%funcs_id( i_func_libxc ) = lxc_functional_get_number( functional_tag )
          IF( is_polarized ) THEN
             CALL lxc_func_init( this%libxc%funcs( i_func_libxc ), this%libxc%funcs_id( i_func_libxc ), XC_POLARIZED)
          ELSE
             CALL lxc_func_init( this%libxc%funcs( i_func_libxc ), this%libxc%funcs_id( i_func_libxc ), XC_UNPOLARIZED)
          ENDIF
          this%libxc%infos( i_func_libxc ) = lxc_func_get_info( this%libxc%funcs( i_func_libxc ) )
          !BLOCK
          ! real(real_8) omega, alpha, beta
          ! call xc_f03_hyb_cam_coef(this%funcs( i_func_libxc ),omega,alpha, beta)
          !END BLOCK
          !
       ELSE
          ! Internal xc functionals
          CALL this%cpfunc%init ( i_func, cp_xc_env )
       ENDIF
    ENDDO
    CALL this%hfx%init ( cp_xc_env )

    this%libxc%n_funcs = i_func_libxc

    IF (this%libxc%n_funcs > 0)  this%use_libxc  = .TRUE.
    IF (this%cpfunc%n_funcs > 0) this%use_cpfunc = .TRUE.

    !>

!!$    select case(xc_f03_func_info_get_kind(this%x_info))
!!$    case(XC_EXCHANGE)
!!$       write(*, '(a)') 'Exchange'
!!$    case(XC_CORRELATION)
!!$       write(*, '(a)') 'Correlation'
!!$    case(XC_EXCHANGE_CORRELATION)
!!$       write(*, '(a)') 'Exchange-correlation'
!!$    case(XC_KINETIC)
!!$       write(*, '(a)') 'Kinetic'
!!$    end select
!!$    write(*,*) trim(xc_f03_func_info_get_name(this%x_info))
!!$    select case(xc_f03_func_info_get_family(this%x_info))
!!$    case (XC_FAMILY_LDA);       write(*,'(a)') "LDA"
!!$    case (XC_FAMILY_GGA);       write(*,'(a)') "GGA"
!!$    case (XC_FAMILY_HYB_GGA);   write(*,'(a)') "Hybrid GGA"
!!$    case (XC_FAMILY_MGGA);      write(*,'(a)') "MGGA"
!!$    case (XC_FAMILY_HYB_MGGA);  write(*,'(a)') "Hybrid MGGA"
!!$    case (XC_FAMILY_LCA);       write(*,'(a)') "LCA"
!!$    end select
!!$
!!$
!!$    select case(xc_f03_func_info_get_kind(this%c_info))
!!$    case(XC_EXCHANGE)
!!$       write(*, '(a)') 'Exchange'
!!$    case(XC_CORRELATION)
!!$       write(*, '(a)') 'Correlation'
!!$    case(XC_EXCHANGE_CORRELATION)
!!$       write(*, '(a)') 'Exchange-correlation'
!!$    case(XC_KINETIC)
!!$       write(*, '(a)') 'Kinetic'
!!$    end select
!!$    write(*,*) trim(xc_f03_func_info_get_name(this%c_info))
!!$    select case(xc_f03_func_info_get_family(this%c_info))
!!$    case (XC_FAMILY_LDA);       write(*,'(a)') "LDA"
!!$    case (XC_FAMILY_GGA);       write(*,'(a)') "GGA"
!!$    case (XC_FAMILY_HYB_GGA);   write(*,'(a)') "Hybrid GGA"
!!$    case (XC_FAMILY_MGGA);      write(*,'(a)') "MGGA"
!!$    case (XC_FAMILY_HYB_MGGA);  write(*,'(a)') "Hybrid MGGA"
!!$    case (XC_FAMILY_LCA);       write(*,'(a)') "LCA"
!!$    end select
!!$
    !<

    this%init = .TRUE.

  END SUBROUTINE cp_xc_create
  ! ==================================================================
  SUBROUTINE cp_xc_get( cp_xc, tgc, tgcx, tgcc, ttau, thybrid, mhfx, msrx, &
                        phfx, srxa, cam_alpha, cam_beta )
    ! ==--------------------------------------------------------------==
    ! Gets value of a chosen (sub)set of private members and translates
    ! them into variables compatible with func1, func2 etc.
    !                               14.08.2017 M. P. Bircher @ LCBC/EPFL
 
    CLASS(cp_xc_t), INTENT(in)               :: cp_xc
    LOGICAL, OPTIONAL, INTENT(out)           :: tgc, tgcx, tgcc, ttau, thybrid
    INTEGER, OPTIONAL, INTENT(out)           :: mhfx, msrx
    REAL(real_8), OPTIONAL, INTENT(out)      :: phfx, srxa, cam_alpha, cam_beta

    IF (present(tgc)) THEN
       tgc = (cp_xc%is_gga .OR. cp_xc%is_mgga)
    ENDIF
    IF (present(tgcx)) THEN
       tgcx = (cp_xc%is_gga .OR. cp_xc%is_mgga)
    ENDIF
    IF (present(tgcc)) THEN
       tgcc = (cp_xc%is_gga .OR. cp_xc%is_mgga)
    ENDIF
    IF (present(ttau)) THEN
       ttau = cp_xc%is_mgga
    ENDIF

    !
    ! HFX
    !
    IF (present(thybrid)) THEN
       thybrid = cp_xc%is_hybrid
    ENDIF
    IF (present(mhfx)) THEN
       IF (cp_xc%is_hybrid) THEN
          mhfx = 1
       ELSE
          mhfx = 0
       ENDIF
    ENDIF
    IF (present(phfx)) THEN
       phfx = cp_xc%hfx%scale
    ENDIF

    !
    ! Screened exchange
    !
    IF (present(msrx)) THEN
       IF (cp_xc%is_hybrid) THEN
          IF (cp_xc%hfx%is_coulomb_attenuated) THEN
             msrx = msrx_is_CAM
          ELSEIF (cp_xc%hfx%is_screened_ashcroft) THEN
             msrx = msrx_is_ashcroft
          ELSEIF (cp_xc%hfx%is_screened_erfc) THEN
             msrx = msrx_is_erfc
          ELSEIF (cp_xc%hfx%is_screened_exp) THEN
             msrx = msrx_is_exp
          ELSE
             msrx = 0
          ENDIF
       ELSE
          msrx = 0
       ENDIF
    ENDIF
    IF (present(srxa)) THEN
       srxa      = cp_xc%hfx%gamma
    ENDIF
    IF (present(cam_alpha)) THEN
       cam_alpha = cp_xc%hfx%alpha
    ENDIF
    IF (present(cam_beta)) THEN
       cam_beta  = cp_xc%hfx%beta
    ENDIF
 
  END SUBROUTINE cp_xc_get
  ! ==================================================================

  ! ==================================================================
  ! Subtype initialisation
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE cpfunc_func_init( cpfunc, env_index, cp_xc_env )
    ! ==--------------------------------------------------------------==
    ! Initialisation of internal xc functionals
    !                               09.10.2016 M. P. Bircher @ LCBC/EPFL

    CLASS(cp_func_t), INTENT(INOUT)          :: cpfunc
    TYPE(cp_xc_env_t), INTENT(IN)            :: cp_xc_env
    INTEGER, INTENT(IN)                      :: env_index

    REAL(real_8), PARAMETER                  :: axlsda = -0.9305257363491_real_8

    CHARACTER(len=*), PARAMETER :: procedureN = 'cpfunc_func_init'

    SELECT CASE(TRIM(ADJUSTL(cp_xc_env%funcs( env_index ))))
       !
       ! Special cases
       ! GGA Exchange
    CASE ("GGA_X_PBE_R")
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_GGA_X_REVPBE"
       cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
       !
       ! GGA Correlation
    CASE ("GGA_C_PBE")
       IF (cpfunc%use_compatibility_mode) THEN
          ! Follows CPMD OLDCODE
          cpfunc%n_funcs                  = cpfunc%n_funcs + 1
          cpfunc%funcs( cpfunc%n_funcs )  = "CP_GGA_C_PBE_FLEX"
          cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
          cp_gga_c_param%pbe%lda_c        = "CP_LDA_C_PZ"
          cp_gga_c_param%pbe%pbe_beta     = 0.06672455060314922_real_8
          cp_gga_c_param%pbe%pbe_gamma    = 0.031090690869654895_real_8 
          cp_gga_c_param%init             = .true.
       ELSE
          cpfunc%n_funcs                  = cpfunc%n_funcs + 1
          cpfunc%funcs( cpfunc%n_funcs )  = "CP_GGA_C_PBE"
          cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
       END IF
    CASE ("GGA_C_PBE_SOL")
       IF (cpfunc%use_compatibility_mode) THEN
          ! Follows CPMD OLDCODE
          cpfunc%n_funcs                  = cpfunc%n_funcs + 1
          cpfunc%funcs( cpfunc%n_funcs )  = "CP_GGA_C_PBE_FLEX"
          cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
          cp_gga_c_param%pbe%lda_c        = "CP_LDA_C_PZ"
          cp_gga_c_param%pbe%pbe_beta     = 0.0460_real_8
          cp_gga_c_param%pbe%pbe_gamma    = 0.031090690869654895_real_8 
          cp_gga_c_param%init             = .true.
       ELSE
          cpfunc%n_funcs                  = cpfunc%n_funcs + 1
          cpfunc%funcs( cpfunc%n_funcs )  = "CP_GGA_C_PBE_SOL"
          cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
       END IF
       !
       ! GGA
       ! Combinations
    CASE ("GGA_XC_HCTH_93") 
       cpfunc%n_funcs                     = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )     = "CP_GGA_X_HCTH"
       cpfunc%scales( cpfunc%n_funcs )    = cp_xc_env%scales( env_index )
       cp_gga_x_param%hcth%gamma          =  0.004_real_8
       cp_gga_x_param%hcth%c_0            =  0.109320e+01_real_8
       cp_gga_x_param%hcth%c              = (/ -0.744056e+00_real_8, &
                                                0.559920e+01_real_8, &
                                               -0.678549e+01_real_8, &
                                                0.449357e+01_real_8 /)
       cp_gga_x_param%init                = .true.
       !
       cpfunc%n_funcs                     = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )     = "CP_GGA_C_HCTH"
       cpfunc%scales( cpfunc%n_funcs )    = cp_xc_env%scales( env_index )
       cp_gga_c_param%hcth%gamma_parallel =  0.200000_real_8
       cp_gga_c_param%hcth%gamma_opposite =  0.006000_real_8
       cp_gga_c_param%hcth%c_0_parallel   =  0.222601_real_8
       cp_gga_c_param%hcth%c_0_opposite   =  0.729974_real_8
       cp_gga_c_param%hcth%c_parallel     =  (/ -0.338622e-01_real_8, &
                                                -0.125170e-01_real_8, &
                                                -0.802496e+00_real_8, &
                                                 0.155396e+01_real_8 /)
       cp_gga_c_param%hcth%c_opposite     =  (/  0.335287e+01_real_8, &
                                                -0.115430e+02_real_8, &
                                                 0.808564e+01_real_8, &
                                                -0.447857e+01_real_8 /)
       cp_gga_c_param%init                = .true.
    CASE ("GGA_XC_HCTH_120") 
       cpfunc%n_funcs                     = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )     = "CP_GGA_X_HCTH"
       cpfunc%scales( cpfunc%n_funcs )    = cp_xc_env%scales( env_index )
       cp_gga_x_param%hcth%gamma          =  0.004_real_8
       cp_gga_x_param%hcth%c_0            =  0.109163e+01_real_8
       cp_gga_x_param%hcth%c              = (/ -0.747215e+00_real_8, &
                                                0.507833e+01_real_8, &
                                               -0.410746e+01_real_8, &
                                                0.117173e+01_real_8 /)
       cp_gga_x_param%init                = .true.
       !
       cpfunc%n_funcs                     = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )     = "CP_GGA_C_HCTH"
       cpfunc%scales( cpfunc%n_funcs )    = cp_xc_env%scales( env_index )
       cp_gga_c_param%hcth%gamma_parallel =  0.200000_real_8
       cp_gga_c_param%hcth%gamma_opposite =  0.006000_real_8
       cp_gga_c_param%hcth%c_0_parallel   =  0.489508_real_8
       cp_gga_c_param%hcth%c_0_opposite   =  0.514730_real_8
       cp_gga_c_param%hcth%c_parallel     =  (/ -0.260699e+00_real_8, &
                                                 0.432917e+00_real_8, &
                                                -0.199274e+01_real_8, &
                                                 0.248531e+01_real_8 /)
       cp_gga_c_param%hcth%c_opposite     =  (/  0.692982e+01_real_8, &
                                                -0.247073e+02_real_8, &
                                                 0.231098e+02_real_8, &
                                                -0.113234e+02_real_8 /)
       cp_gga_c_param%init                = .true.
    CASE ("GGA_XC_HCTH_147") 
       cpfunc%n_funcs                     = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )     = "CP_GGA_X_HCTH"
       cpfunc%scales( cpfunc%n_funcs )    = cp_xc_env%scales( env_index )
       cp_gga_x_param%hcth%gamma          =  0.004_real_8
       cp_gga_x_param%hcth%c_0            =  0.109025e+01_real_8
       cp_gga_x_param%hcth%c              = (/ -0.799194e+00_real_8, &
                                                0.557212e+01_real_8, &
                                               -0.586760e+01_real_8, &
                                                0.304544e+01_real_8 /)
       cp_gga_x_param%init                = .true.
       !
       cpfunc%n_funcs                     = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )     = "CP_GGA_C_HCTH"
       cpfunc%scales( cpfunc%n_funcs )    = cp_xc_env%scales( env_index )
       cp_gga_c_param%hcth%gamma_parallel =  0.200000_real_8
       cp_gga_c_param%hcth%gamma_opposite =  0.006000_real_8
       cp_gga_c_param%hcth%c_0_parallel   =  0.562576_real_8
       cp_gga_c_param%hcth%c_0_opposite   =  0.542352_real_8
       cp_gga_c_param%hcth%c_parallel     =  (/  0.171436e-01_real_8, &
                                                -0.130636e+01_real_8, &
                                                 0.105747e+01_real_8, &
                                                 0.885429e+00_real_8 /)
       cp_gga_c_param%hcth%c_opposite     =  (/  0.701464e+01_real_8, &
                                                -0.283822e+02_real_8, &
                                                 0.350329e+02_real_8, &
                                                -0.204284e+02_real_8 /)
       cp_gga_c_param%init                = .true.
    CASE ("GGA_XC_HCTH_407") 
       cpfunc%n_funcs                     = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )     = "CP_GGA_X_HCTH"
       cpfunc%scales( cpfunc%n_funcs )    = cp_xc_env%scales( env_index )
       cp_gga_x_param%hcth%gamma          =  0.004_real_8
       cp_gga_x_param%hcth%c_0            =  0.108184e+01_real_8
       cp_gga_x_param%hcth%c              = (/ -0.518339e+00_real_8, &
                                                0.342562e+01_real_8, &
                                               -0.262901e+01_real_8, &
                                                0.228855e+01_real_8 /)
       cp_gga_x_param%init                = .true.
       !
       cpfunc%n_funcs                     = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )     = "CP_GGA_C_HCTH"
       cpfunc%scales( cpfunc%n_funcs )    = cp_xc_env%scales( env_index )
       cp_gga_c_param%hcth%gamma_parallel =  0.200000_real_8
       cp_gga_c_param%hcth%gamma_opposite =  0.006000_real_8
       cp_gga_c_param%hcth%c_0_parallel   =  1.18777_real_8
       cp_gga_c_param%hcth%c_0_opposite   =  0.589076_real_8
       cp_gga_c_param%hcth%c_parallel     =  (/ -0.240292e+01_real_8, &
                                                 0.561741e+01_real_8, &
                                                -0.917923e+01_real_8, &
                                                 0.624798e+01_real_8 /)
       cp_gga_c_param%hcth%c_opposite     =  (/  0.442374e+01_real_8, &
                                                -0.192218e+02_real_8, &
                                                 0.425721e+02_real_8, &
                                                -0.420052e+02_real_8 /)
       cp_gga_c_param%init                = .true.
    CASE ("GGA_XC_PBE") 
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_GGA_X_PBE"
       cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
       !
       IF (cpfunc%use_compatibility_mode) THEN
          ! Follows CPMD OLDCODE
          cpfunc%n_funcs                  = cpfunc%n_funcs + 1
          cpfunc%funcs( cpfunc%n_funcs )  = "CP_GGA_C_PBE_FLEX"
          cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
          cp_gga_c_param%pbe%lda_c        = "CP_LDA_C_PZ"
          cp_gga_c_param%pbe%pbe_beta     = 0.06672455060314922_real_8
          cp_gga_c_param%pbe%pbe_gamma    = 0.031090690869654895_real_8 
          cp_gga_c_param%init             = .true.
       ELSE
          cpfunc%n_funcs                  = cpfunc%n_funcs + 1
          cpfunc%funcs( cpfunc%n_funcs )  = "CP_GGA_C_PBE"
          cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
       END IF
    CASE ("GGA_XC_REVPBE","GGA_XC_PBE_R") 
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_GGA_X_REVPBE"
       cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
       !
       IF (cpfunc%use_compatibility_mode) THEN
          ! Follows CPMD OLDCODE
          cpfunc%n_funcs                  = cpfunc%n_funcs + 1
          cpfunc%funcs( cpfunc%n_funcs )  = "CP_GGA_C_PBE_FLEX"
          cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
          cp_gga_c_param%pbe%lda_c        = "CP_LDA_C_PZ"
          cp_gga_c_param%pbe%pbe_beta     = 0.06672455060314922_real_8
          cp_gga_c_param%pbe%pbe_gamma    = 0.031090690869654895_real_8 
          cp_gga_c_param%init             = .true.
       ELSE
          cpfunc%n_funcs                  = cpfunc%n_funcs + 1
          cpfunc%funcs( cpfunc%n_funcs )  = "CP_GGA_C_PBE"
          cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
       END IF
    CASE ("GGA_XC_PBE_SOL") 
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_GGA_X_PBE_SOL"
       cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
       !
       IF (cpfunc%use_compatibility_mode) THEN
          ! Follows CPMD OLDCODE
          cpfunc%n_funcs                  = cpfunc%n_funcs + 1
          cpfunc%funcs( cpfunc%n_funcs )  = "CP_GGA_C_PBE_FLEX"
          cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
          cp_gga_c_param%pbe%lda_c        = "CP_LDA_C_PZ"
          cp_gga_c_param%pbe%pbe_beta     = 0.0460_real_8
          cp_gga_c_param%pbe%pbe_gamma    = 0.031090690869654895_real_8 
          cp_gga_c_param%init             = .true.
       ELSE
          cpfunc%n_funcs                  = cpfunc%n_funcs + 1
          cpfunc%funcs( cpfunc%n_funcs )  = "CP_GGA_C_PBE_SOL"
          cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
       END IF
    CASE ("GGA_XC_BLYP") 
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_GGA_X_B88"
       cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
       !
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_GGA_C_LYP"
       cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
    CASE ("GGA_XC_BP","GGA_XC_BP86") 
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_GGA_X_B88"
       cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
       !
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_GGA_C_P86"
       cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
    CASE ("GGA_XC_OLYP") 
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_GGA_X_OPTX"
       cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
       !
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_GGA_C_LYP"
       cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
    CASE ("GGA_XC_OPBE") 
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_GGA_X_OPTX"
       cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
       !
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_GGA_C_PBE"
       cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
       !
       ! MGGA
       ! Fallback options
       ! CASE ("MGGA_X_VS98")
       !    cpfunc%n_funcs                    = cpfunc%n_funcs + 1
       !    cpfunc%funcs( cpfunc%n_funcs )    = "CP_MGGA_X_VS98"
       !    cpfunc%scales( cpfunc%n_funcs )   = cp_xc_env%scales( env_index )
       !    cp_mgga_x_param%VS98%r1           = -9.800683e-01_real_8
       !    cp_mgga_x_param%VS98%r2           = -3.556788e-03_real_8
       !    cp_mgga_x_param%VS98%r3           =  6.250326e-03_real_8
       !    cp_mgga_x_param%VS98%r4           = -2.354518e-05_real_8
       !    cp_mgga_x_param%VS98%r5           = -1.282732e-04_real_8
       !    cp_mgga_x_param%VS98%r6           =  3.574822e-04_real_8
       ! CASE ("MGGA_C_VS98")
       !    cpfunc%n_funcs                    = cpfunc%n_funcs + 1
       !    cpfunc%funcs( cpfunc%n_funcs )    = "CP_MGGA_C_VS98"
       !    cpfunc%scales( cpfunc%n_funcs )   = cp_xc_env%scales( env_index )
       !    cp_mgga_c_param%VS98%r7           =  7.035010e-01_real_8
       !    cp_mgga_c_param%VS98%r8           =  7.694574e-03_real_8
       !    cp_mgga_c_param%VS98%r9           =  5.152765e-02_real_8
       !    cp_mgga_c_param%VS98%r10          =  3.394308e-05_real_8
       !    cp_mgga_c_param%VS98%r11          = -1.269420e-03_real_8
       !    cp_mgga_c_param%VS98%r12          =  1.296118e-03_real_8
       !    cp_mgga_c_param%VS98%r13          =  3.270912e-01_real_8
       !    cp_mgga_c_param%VS98%r14          = -3.228915e-02_real_8
       !    cp_mgga_c_param%VS98%r15          = -2.942406e-02_real_8
       !    cp_mgga_c_param%VS98%r16          =  2.134222e-03_real_8
       !    cp_mgga_c_param%VS98%r17          = -5.451559e-03_real_8
       !    cp_mgga_c_param%VS98%r18          =  1.577575e-02_real_8
    CASE ("MGGA_C_TPSS")
       IF (cpfunc%use_compatibility_mode) THEN
          ! Follows CPMD OLDCODE
          cpfunc%n_funcs                  = cpfunc%n_funcs + 1
          cpfunc%funcs( cpfunc%n_funcs )  = "CP_MGGA_C_TPSS"
          cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
          cp_gga_c_param%pbe%lda_c        = "CP_LDA_C_PZ"
          cp_gga_c_param%pbe%pbe_beta     = 0.06672455060314922_real_8
          cp_gga_c_param%pbe%pbe_gamma    = 0.031090690869654895_real_8 
          cp_gga_c_param%init             = .true.
       ELSE
          cpfunc%n_funcs                  = cpfunc%n_funcs + 1
          cpfunc%funcs( cpfunc%n_funcs )  = "CP_MGGA_C_TPSS"
          cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
          cp_gga_c_param%pbe%lda_c        = "CP_LDA_C_PW"
          cp_gga_c_param%pbe%pbe_beta     = 0.06672455060314922_real_8
          cp_gga_c_param%pbe%pbe_gamma    = 0.031090690869654895_real_8 
          cp_gga_c_param%init             = .true.
       END IF
       !
       ! MGGA
       ! Combinations
    CASE ("MGGA_XC_TPSS") 
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_MGGA_X_TPSS"
       cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
       !
       IF (cpfunc%use_compatibility_mode) THEN
          ! Follows CPMD OLDCODE
          cpfunc%n_funcs                  = cpfunc%n_funcs + 1
          cpfunc%funcs( cpfunc%n_funcs )  = "CP_MGGA_C_TPSS"
          cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
          cp_gga_c_param%pbe%lda_c        = "CP_LDA_C_PZ"
          cp_gga_c_param%pbe%pbe_beta     = 0.06672455060314922_real_8
          cp_gga_c_param%pbe%pbe_gamma    = 0.031090690869654895_real_8 
          cp_gga_c_param%init             = .true.
       ELSE
          cpfunc%n_funcs                  = cpfunc%n_funcs + 1
          cpfunc%funcs( cpfunc%n_funcs )  = "CP_MGGA_C_TPSS"
          cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
          cp_gga_c_param%pbe%lda_c        = "CP_LDA_C_PW"
          cp_gga_c_param%pbe%pbe_beta     = 0.06672455060314922_real_8
          cp_gga_c_param%pbe%pbe_gamma    = 0.031090690869654895_real_8 
          cp_gga_c_param%init             = .true.
       END IF
    CASE ("MGGA_XC_M06_L")
       !
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_MGGA_X_M05_M06"
       cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
       cp_mgga_x_param%M05_M06%add_vs98=  .true.
       cp_mgga_x_param%M05_M06%at00    =  3.987756e-01_real_8
       cp_mgga_x_param%M05_M06%at01    =  2.548219e-01_real_8
       cp_mgga_x_param%M05_M06%at02    =  3.923994e-01_real_8
       cp_mgga_x_param%M05_M06%at03    = -2.103655e+00_real_8
       cp_mgga_x_param%M05_M06%at04    = -6.302147e+00_real_8
       cp_mgga_x_param%M05_M06%at05    =  1.097615e+01_real_8
       cp_mgga_x_param%M05_M06%at06    =  3.097273e+01_real_8
       cp_mgga_x_param%M05_M06%at07    = -2.318489e+01_real_8
       cp_mgga_x_param%M05_M06%at08    = -5.673480e+01_real_8
       cp_mgga_x_param%M05_M06%at09    =  2.160364e+01_real_8
       cp_mgga_x_param%M05_M06%at10    =  3.421814e+01_real_8
       cp_mgga_x_param%M05_M06%at11    = -9.049762e+00_real_8
       cp_mgga_x_param%VS98%r1         =  6.012244e-01_real_8*axlsda
       cp_mgga_x_param%VS98%r2         =  4.748822e-03_real_8*axlsda
       cp_mgga_x_param%VS98%r3         = -8.635108e-03_real_8*axlsda
       cp_mgga_x_param%VS98%r4         = -9.308062e-06_real_8*axlsda
       cp_mgga_x_param%VS98%r5         =  4.482811e-05_real_8*axlsda
       cp_mgga_x_param%VS98%r6         =  0.000000e+00_real_8
       cp_mgga_x_param%init            =  .true.
       !
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_MGGA_C_M05_M06"
       cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
       cp_mgga_c_param%M05_M06%add_vs98=  .true.
       cp_mgga_c_param%M05_M06%sopp0   =  6.042374e-01_real_8
       cp_mgga_c_param%M05_M06%sopp1   =  1.776783e+02_real_8
       cp_mgga_c_param%M05_M06%sopp2   = -2.513252e+02_real_8
       cp_mgga_c_param%M05_M06%sopp3   =  7.635173e+01_real_8
       cp_mgga_c_param%M05_M06%sopp4   = -1.255699e+01_real_8
       cp_mgga_c_param%M05_M06%sss0    =  5.349466e-01_real_8
       cp_mgga_c_param%M05_M06%sss1    =  5.396620e-01_real_8
       cp_mgga_c_param%M05_M06%sss2    = -3.161217e+01_real_8
       cp_mgga_c_param%M05_M06%sss3    =  5.149592e+01_real_8
       cp_mgga_c_param%M05_M06%sss4    = -2.919613e+01_real_8
       cp_mgga_c_param%VS98%r7         =  3.957626e-01_real_8
       cp_mgga_c_param%VS98%r8         = -5.614546e-01_real_8
       cp_mgga_c_param%VS98%r9         =  1.403963e-02_real_8
       cp_mgga_c_param%VS98%r10        =  9.831442e-04_real_8
       cp_mgga_c_param%VS98%r11        = -3.577176e-03_real_8
       cp_mgga_c_param%VS98%r12        =  0.000000e+00_real_8
       cp_mgga_c_param%VS98%r13        =  4.650534e-01_real_8
       cp_mgga_c_param%VS98%r14        =  1.617589e-01_real_8
       cp_mgga_c_param%VS98%r15        =  1.833657e-01_real_8
       cp_mgga_c_param%VS98%r16        =  4.692100e-04_real_8
       cp_mgga_c_param%VS98%r17        = -4.990573e-03_real_8
       cp_mgga_c_param%VS98%r18        =  0.000000e+00_real_8
       cp_mgga_c_param%init            =  .true.
    CASE ("MGGA_XC_REVM06_L")
       !
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_MGGA_X_M05_M06"
       cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
       cp_mgga_x_param%M05_M06%add_vs98=  .true.
       cp_mgga_x_param%M05_M06%at00    =  1.423227252e+00_real_8 
       cp_mgga_x_param%M05_M06%at01    =  4.718204380e-01_real_8
       cp_mgga_x_param%M05_M06%at02    = -1.675557010e-01_real_8
       cp_mgga_x_param%M05_M06%at03    = -2.501542620e-01_real_8
       cp_mgga_x_param%M05_M06%at04    =  6.248758800e-02_real_8
       cp_mgga_x_param%M05_M06%at05    =  7.335012400e-01_real_8
       cp_mgga_x_param%M05_M06%at06    = -2.359736776e+00_real_8
       cp_mgga_x_param%M05_M06%at07    = -1.436594372e+00_real_8
       cp_mgga_x_param%M05_M06%at08    =  4.446437930e-01_real_8
       cp_mgga_x_param%M05_M06%at09    =  1.529925054e+00_real_8
       cp_mgga_x_param%M05_M06%at10    =  2.053941717e+00_real_8
       cp_mgga_x_param%M05_M06%at11    = -3.653603100e-02_real_8
       cp_mgga_x_param%VS98%r1         = -4.232272520e-01_real_8*axlsda
       cp_mgga_x_param%VS98%r2         =  0.000000000e+00_real_8*axlsda
       cp_mgga_x_param%VS98%r3         =  3.724234000e-03_real_8*axlsda
       cp_mgga_x_param%VS98%r4         =  0.000000000e+00_real_8*axlsda
       cp_mgga_x_param%VS98%r5         =  0.000000000e+00_real_8*axlsda
       cp_mgga_x_param%init            =  .true.
       !
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_MGGA_C_M05_M06"
       cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
       cp_mgga_c_param%M05_M06%add_vs98=  .true.
       cp_mgga_c_param%M05_M06%sopp0   =  3.443606960e-01_real_8 
       cp_mgga_c_param%M05_M06%sopp1   = -5.570802420e-01_real_8
       cp_mgga_c_param%M05_M06%sopp2   = -2.009821162e+00_real_8
       cp_mgga_c_param%M05_M06%sopp3   = -1.857641887e+00_real_8
       cp_mgga_c_param%M05_M06%sopp4   = -1.076639864e+00_real_8
       cp_mgga_c_param%M05_M06%sss0    =  1.227659748e+00_real_8
       cp_mgga_c_param%M05_M06%sss1    =  8.552012830e-01_real_8
       cp_mgga_c_param%M05_M06%sss2    = -3.113346677e+00_real_8
       cp_mgga_c_param%M05_M06%sss3    = -2.239678026e+00_real_8
       cp_mgga_c_param%M05_M06%sss4    =  3.546389620e-01_real_8
       cp_mgga_c_param%VS98%r7         =  4.00714600e-01_real_8 
       cp_mgga_c_param%VS98%r8         =  1.57965690e-02_real_8
       cp_mgga_c_param%VS98%r9         = -3.26809840e-02_real_8
       cp_mgga_c_param%VS98%r10        =  0.00000000e+00_real_8
       cp_mgga_c_param%VS98%r11        =  0.00000000e+00_real_8
       cp_mgga_c_param%VS98%r12        =  1.26013200e-03_real_8
       cp_mgga_c_param%VS98%r13        = -5.38821292e-01_real_8
       cp_mgga_c_param%VS98%r14        = -2.82960300e-02_real_8
       cp_mgga_c_param%VS98%r15        =  2.38896960e-02_real_8
       cp_mgga_c_param%VS98%r16        =  0.00000000e+00_real_8
       cp_mgga_c_param%VS98%r17        =  0.00000000e+00_real_8
       cp_mgga_c_param%VS98%r18        = -2.43790200e-03_real_8
       cp_mgga_c_param%init            =  .true.
    CASE ("MGGA_XC_M11_L")
       !
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_MGGA_X_M08_M11"
       cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
       cp_mgga_x_param%M08_M11%lda_has_lc = .true.
       cp_mgga_x_param%M08_M11%at00    =  8.121131e-01_real_8
       cp_mgga_x_param%M08_M11%at01    =  1.738124e+01_real_8
       cp_mgga_x_param%M08_M11%at02    =  1.154007e+00_real_8
       cp_mgga_x_param%M08_M11%at03    =  6.869556e+01_real_8
       cp_mgga_x_param%M08_M11%at04    =  1.016864e+02_real_8
       cp_mgga_x_param%M08_M11%at05    = -5.887467e+00_real_8
       cp_mgga_x_param%M08_M11%at06    =  4.517409e+01_real_8
       cp_mgga_x_param%M08_M11%at07    = -2.773149e+00_real_8
       cp_mgga_x_param%M08_M11%at08    = -2.617211e+01_real_8
       cp_mgga_x_param%M08_M11%at09    =  0.000000e+00_real_8
       cp_mgga_x_param%M08_M11%at10    =  0.000000e+00_real_8
       cp_mgga_x_param%M08_M11%at11    =  0.000000e+00_real_8
       cp_mgga_x_param%M08_M11%bt00    =  1.878869e-01_real_8
       cp_mgga_x_param%M08_M11%bt01    = -1.653877e+01_real_8
       cp_mgga_x_param%M08_M11%bt02    =  6.755753e-01_real_8
       cp_mgga_x_param%M08_M11%bt03    = -7.567572e+01_real_8
       cp_mgga_x_param%M08_M11%bt04    = -1.040272e+02_real_8
       cp_mgga_x_param%M08_M11%bt05    =  1.831853e+01_real_8
       cp_mgga_x_param%M08_M11%bt06    = -5.573352e+01_real_8
       cp_mgga_x_param%M08_M11%bt07    = -3.520210e+00_real_8
       cp_mgga_x_param%M08_M11%bt08    =  3.724276e+01_real_8
       cp_mgga_x_param%M08_M11%bt09    =  0.000000e+00_real_8
       cp_mgga_x_param%M08_M11%bt10    =  0.000000e+00_real_8
       cp_mgga_x_param%M08_M11%bt11    =  0.000000e+00_real_8
       cp_mgga_x_param%M08_M11%ct00    = -4.386615e-01_real_8
       cp_mgga_x_param%M08_M11%ct01    = -1.214016e+02_real_8
       cp_mgga_x_param%M08_M11%ct02    = -1.393573e+02_real_8
       cp_mgga_x_param%M08_M11%ct03    = -2.046649e+00_real_8
       cp_mgga_x_param%M08_M11%ct04    =  2.804098e+01_real_8
       cp_mgga_x_param%M08_M11%ct05    = -1.312258e+01_real_8
       cp_mgga_x_param%M08_M11%ct06    = -6.361819e+00_real_8
       cp_mgga_x_param%M08_M11%ct07    = -8.055758e-01_real_8
       cp_mgga_x_param%M08_M11%ct08    =  3.736551e+00_real_8
       cp_mgga_x_param%M08_M11%ct09    =  0.000000e+00_real_8
       cp_mgga_x_param%M08_M11%ct10    =  0.000000e+00_real_8
       cp_mgga_x_param%M08_M11%ct11    =  0.000000e+00_real_8
       cp_mgga_x_param%M08_M11%dt00    =  1.438662e+00_real_8
       cp_mgga_x_param%M08_M11%dt01    =  1.209465e+02_real_8
       cp_mgga_x_param%M08_M11%dt02    =  1.328252e+02_real_8
       cp_mgga_x_param%M08_M11%dt03    =  1.296355e+01_real_8
       cp_mgga_x_param%M08_M11%dt04    =  5.854866e+00_real_8
       cp_mgga_x_param%M08_M11%dt05    = -3.378162e+00_real_8
       cp_mgga_x_param%M08_M11%dt06    = -4.423393e+01_real_8
       cp_mgga_x_param%M08_M11%dt07    =  6.844475e+00_real_8
       cp_mgga_x_param%M08_M11%dt08    =  1.949541e+01_real_8
       cp_mgga_x_param%M08_M11%dt09    =  0.000000e+00_real_8
       cp_mgga_x_param%M08_M11%dt10    =  0.000000e+00_real_8
       cp_mgga_x_param%M08_M11%dt11    =  0.000000e+00_real_8
       cp_mgga_x_param%init            =  .true.
       !
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_MGGA_C_M08_M11"
       cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
       cp_gga_c_param%pbe%lda_c        = "CP_LDA_C_PW"
       cp_gga_c_param%pbe%pbe_beta     = 0.06672455060314922_real_8
       cp_gga_c_param%pbe%pbe_gamma    = 0.031090690869654895_real_8
       cp_gga_c_param%init             = .true. 
       cp_mgga_c_param%M08_M11%at00    =  1.000000e+00_real_8
       cp_mgga_c_param%M08_M11%at01    =  0.000000e+00_real_8
       cp_mgga_c_param%M08_M11%at02    =  2.750880e+00_real_8
       cp_mgga_c_param%M08_M11%at03    = -1.562287e+01_real_8
       cp_mgga_c_param%M08_M11%at04    =  9.363381e+00_real_8
       cp_mgga_c_param%M08_M11%at05    =  2.141024e+01_real_8
       cp_mgga_c_param%M08_M11%at06    = -1.424975e+01_real_8
       cp_mgga_c_param%M08_M11%at07    = -1.134712e+01_real_8
       cp_mgga_c_param%M08_M11%at08    =  1.022365e+01_real_8
       cp_mgga_c_param%M08_M11%at09    =  0.000000e+00_real_8
       cp_mgga_c_param%M08_M11%at10    =  0.000000e+00_real_8
       cp_mgga_c_param%M08_M11%at11    =  0.000000e+00_real_8
       cp_mgga_c_param%M08_M11%bt00    =  1.000000e+00_real_8
       cp_mgga_c_param%M08_M11%bt01    = -9.082060e+00_real_8
       cp_mgga_c_param%M08_M11%bt02    =  6.134682e+00_real_8
       cp_mgga_c_param%M08_M11%bt03    = -1.333216e+01_real_8
       cp_mgga_c_param%M08_M11%bt04    = -1.464115e+01_real_8
       cp_mgga_c_param%M08_M11%bt05    =  1.713143e+01_real_8
       cp_mgga_c_param%M08_M11%bt06    =  2.480738e+00_real_8
       cp_mgga_c_param%M08_M11%bt07    = -1.007036e+01_real_8
       cp_mgga_c_param%M08_M11%bt08    = -1.117521e-01_real_8
       cp_mgga_c_param%M08_M11%bt09    =  0.000000e+00_real_8
       cp_mgga_c_param%M08_M11%bt10    =  0.000000e+00_real_8
       cp_mgga_c_param%M08_M11%bt11    =  0.000000e+00_real_8
       cp_mgga_c_param%init            =  .true.
    CASE ("MGGA_XC_MN12_L")
       !
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_MGGA_X_MN12_MN15"
       cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
       cp_mgga_x_param%MN12_MN15%cc000 =  6.735981e-01_real_8
       cp_mgga_x_param%MN12_MN15%cc001 = -2.270598e+00_real_8
       cp_mgga_x_param%MN12_MN15%cc002 = -2.613712e+00_real_8
       cp_mgga_x_param%MN12_MN15%cc003 =  3.993609e+00_real_8
       cp_mgga_x_param%MN12_MN15%cc004 =  4.635575e+00_real_8
       cp_mgga_x_param%MN12_MN15%cc005 =  1.250676e+00_real_8
       cp_mgga_x_param%MN12_MN15%cc010 =  8.444920e-01_real_8
       cp_mgga_x_param%MN12_MN15%cc011 = -1.301173e+01_real_8
       cp_mgga_x_param%MN12_MN15%cc012 = -1.777730e+01_real_8
       cp_mgga_x_param%MN12_MN15%cc013 = -4.627211e+00_real_8
       cp_mgga_x_param%MN12_MN15%cc014 =  5.976605e+00_real_8
       cp_mgga_x_param%MN12_MN15%cc020 =  1.142897e+00_real_8
       cp_mgga_x_param%MN12_MN15%cc021 = -2.040226e+01_real_8
       cp_mgga_x_param%MN12_MN15%cc022 = -2.382843e+01_real_8
       cp_mgga_x_param%MN12_MN15%cc023 =  7.119109e+00_real_8
       cp_mgga_x_param%MN12_MN15%cc030 = -2.335726e+01_real_8
       cp_mgga_x_param%MN12_MN15%cc031 = -1.622633e+01_real_8
       cp_mgga_x_param%MN12_MN15%cc032 =  1.482732e+01_real_8
       cp_mgga_x_param%MN12_MN15%cc100 =  1.449285e+00_real_8
       cp_mgga_x_param%MN12_MN15%cc101 =  1.020598e+01_real_8
       cp_mgga_x_param%MN12_MN15%cc102 =  4.407450e+00_real_8
       cp_mgga_x_param%MN12_MN15%cc103 = -2.008193e+01_real_8
       cp_mgga_x_param%MN12_MN15%cc104 = -1.253561e+01_real_8
       cp_mgga_x_param%MN12_MN15%cc110 = -5.435031e+00_real_8
       cp_mgga_x_param%MN12_MN15%cc111 =  1.656736e+01_real_8
       cp_mgga_x_param%MN12_MN15%cc112 =  2.000229e+01_real_8
       cp_mgga_x_param%MN12_MN15%cc113 = -2.513105e+00_real_8
       cp_mgga_x_param%MN12_MN15%cc120 =  9.658436e+00_real_8
       cp_mgga_x_param%MN12_MN15%cc121 = -3.825281e+00_real_8
       cp_mgga_x_param%MN12_MN15%cc122 = -2.500000e+01_real_8
       cp_mgga_x_param%MN12_MN15%cc200 = -2.070080e+00_real_8
       cp_mgga_x_param%MN12_MN15%cc201 = -9.951913e+00_real_8
       cp_mgga_x_param%MN12_MN15%cc202 =  8.731211e-01_real_8
       cp_mgga_x_param%MN12_MN15%cc203 =  2.210891e+01_real_8
       cp_mgga_x_param%MN12_MN15%cc210 =  8.822633e+00_real_8
       cp_mgga_x_param%MN12_MN15%cc211 =  2.499949e+01_real_8
       cp_mgga_x_param%MN12_MN15%cc212 =  2.500000e+01_real_8
       cp_mgga_x_param%MN12_MN15%cc300 =  6.851693e-01_real_8
       cp_mgga_x_param%MN12_MN15%cc301 = -7.406948e-02_real_8
       cp_mgga_x_param%MN12_MN15%cc302 = -6.788000e-01_real_8
       cp_mgga_x_param%init            =  .true.
       !
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_MGGA_C_M08_M11"
       cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
       cp_gga_c_param%pbe%lda_c        = "CP_LDA_C_PW"
       cp_gga_c_param%pbe%pbe_beta     = 0.06672455060314922_real_8
       cp_gga_c_param%pbe%pbe_gamma    = 0.031090690869654895_real_8
       cp_gga_c_param%init             = .true. 
       cp_mgga_c_param%M08_M11%at00    =  8.844610e-01_real_8
       cp_mgga_c_param%M08_M11%at01    = -2.202279e-01_real_8
       cp_mgga_c_param%M08_M11%at02    =  5.701372e+00_real_8
       cp_mgga_c_param%M08_M11%at03    = -2.562378e+00_real_8
       cp_mgga_c_param%M08_M11%at04    = -9.646827e-01_real_8
       cp_mgga_c_param%M08_M11%at05    =  1.982183e-01_real_8
       cp_mgga_c_param%M08_M11%at06    =  1.019976e+01_real_8
       cp_mgga_c_param%M08_M11%at07    =  9.789352e-01_real_8
       cp_mgga_c_param%M08_M11%at08    = -1.512722e+00_real_8
       cp_mgga_c_param%M08_M11%at09    =  0.000000e+00_real_8
       cp_mgga_c_param%M08_M11%at10    =  0.000000e+00_real_8
       cp_mgga_c_param%M08_M11%at11    =  0.000000e+00_real_8
       cp_mgga_c_param%M08_M11%bt00    =  5.323948e-01_real_8
       cp_mgga_c_param%M08_M11%bt01    = -5.831909e+00_real_8
       cp_mgga_c_param%M08_M11%bt02    =  3.882386e+00_real_8
       cp_mgga_c_param%M08_M11%bt03    =  5.878488e+00_real_8
       cp_mgga_c_param%M08_M11%bt04    =  1.493228e+01_real_8
       cp_mgga_c_param%M08_M11%bt05    = -1.374636e+01_real_8
       cp_mgga_c_param%M08_M11%bt06    = -8.492327e+00_real_8
       cp_mgga_c_param%M08_M11%bt07    = -2.486548e+00_real_8
       cp_mgga_c_param%M08_M11%bt08    = -1.822346e+01_real_8
       cp_mgga_c_param%M08_M11%bt09    =  0.000000e+00_real_8
       cp_mgga_c_param%M08_M11%bt10    =  0.000000e+00_real_8
       cp_mgga_c_param%M08_M11%bt11    =  0.000000e+00_real_8
       cp_mgga_c_param%init            =  .true.
    CASE ("MGGA_XC_MN15_L")
       !
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_MGGA_X_MN12_MN15"
       cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
       cp_mgga_x_param%MN12_MN15%cc000 =  0.670864162e0_real_8
       cp_mgga_x_param%MN12_MN15%cc001 = -0.822003903e0_real_8
       cp_mgga_x_param%MN12_MN15%cc002 = -1.022407046e0_real_8
       cp_mgga_x_param%MN12_MN15%cc003 =  1.689460986e0_real_8
       cp_mgga_x_param%MN12_MN15%cc004 = -0.005620320e0_real_8
       cp_mgga_x_param%MN12_MN15%cc005 = -0.110293849e0_real_8
       cp_mgga_x_param%MN12_MN15%cc010 =  0.972245178e0_real_8
       cp_mgga_x_param%MN12_MN15%cc011 = -6.697641991e0_real_8
       cp_mgga_x_param%MN12_MN15%cc012 = -4.322814495e0_real_8
       cp_mgga_x_param%MN12_MN15%cc013 = -6.786641376e0_real_8
       cp_mgga_x_param%MN12_MN15%cc014 = -5.687461462e0_real_8
       cp_mgga_x_param%MN12_MN15%cc020 =  9.419643818e0_real_8
       cp_mgga_x_param%MN12_MN15%cc021 =  11.83939406e0_real_8
       cp_mgga_x_param%MN12_MN15%cc022 =  5.086951311e0_real_8
       cp_mgga_x_param%MN12_MN15%cc023 =  4.302369948e0_real_8
       cp_mgga_x_param%MN12_MN15%cc030 = -8.073440650e0_real_8
       cp_mgga_x_param%MN12_MN15%cc031 =  2.429988978e0_real_8
       cp_mgga_x_param%MN12_MN15%cc032 =  11.09485698e0_real_8
       cp_mgga_x_param%MN12_MN15%cc100 =  1.247333909e0_real_8
       cp_mgga_x_param%MN12_MN15%cc101 =  3.700485291e0_real_8
       cp_mgga_x_param%MN12_MN15%cc102 =  0.867791614e0_real_8
       cp_mgga_x_param%MN12_MN15%cc103 = -0.591190518e0_real_8
       cp_mgga_x_param%MN12_MN15%cc104 = -0.295305435e0_real_8
       cp_mgga_x_param%MN12_MN15%cc110 = -5.825759145e0_real_8
       cp_mgga_x_param%MN12_MN15%cc111 =  2.537532196e0_real_8
       cp_mgga_x_param%MN12_MN15%cc112 =  3.143390933e0_real_8
       cp_mgga_x_param%MN12_MN15%cc113 =  2.939126332e0_real_8
       cp_mgga_x_param%MN12_MN15%cc120 =  0.599342114e0_real_8
       cp_mgga_x_param%MN12_MN15%cc121 =  2.241702738e0_real_8
       cp_mgga_x_param%MN12_MN15%cc122 =  2.035713838e0_real_8
       cp_mgga_x_param%MN12_MN15%cc200 = -1.525344043e0_real_8
       cp_mgga_x_param%MN12_MN15%cc201 = -2.325875691e0_real_8
       cp_mgga_x_param%MN12_MN15%cc202 =  1.141940663e0_real_8
       cp_mgga_x_param%MN12_MN15%cc203 = -1.563165026e0_real_8
       cp_mgga_x_param%MN12_MN15%cc210 =  7.882032871e0_real_8
       cp_mgga_x_param%MN12_MN15%cc211 =  11.93400684e0_real_8
       cp_mgga_x_param%MN12_MN15%cc212 =  9.852928303e0_real_8
       cp_mgga_x_param%MN12_MN15%cc300 =  0.584030245e0_real_8
       cp_mgga_x_param%MN12_MN15%cc301 = -0.720941131e0_real_8
       cp_mgga_x_param%MN12_MN15%cc302 = -2.836037078e0_real_8
       cp_mgga_x_param%init            =  .true.
       !
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_MGGA_C_M08_M11"
       cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
       cp_gga_c_param%pbe%lda_c        = "CP_LDA_C_PW"
       cp_gga_c_param%pbe%pbe_beta     = 0.06672455060314922_real_8
       cp_gga_c_param%pbe%pbe_gamma    = 0.031090690869654895_real_8
       cp_gga_c_param%init             = .true. 
       cp_mgga_c_param%M08_M11%at00    =  0.952058087e+00_real_8
       cp_mgga_c_param%M08_M11%at01    = -0.756954364e+00_real_8
       cp_mgga_c_param%M08_M11%at02    =  5.677396094e+00_real_8
       cp_mgga_c_param%M08_M11%at03    = -5.017104782e+00_real_8
       cp_mgga_c_param%M08_M11%at04    = -5.106540710e+00_real_8
       cp_mgga_c_param%M08_M11%at05    = -4.812053335e+00_real_8
       cp_mgga_c_param%M08_M11%at06    =  3.397640087e+00_real_8
       cp_mgga_c_param%M08_M11%at07    =  1.980041517e+00_real_8
       cp_mgga_c_param%M08_M11%at08    =  1.012310460e+01_real_8
       cp_mgga_c_param%M08_M11%at09    =  0.000000000e+00_real_8
       cp_mgga_c_param%M08_M11%at10    =  0.000000000e+00_real_8
       cp_mgga_c_param%M08_M11%at11    =  0.000000000e+00_real_8
       cp_mgga_c_param%M08_M11%bt00    =  0.819504932e+00_real_8
       cp_mgga_c_param%M08_M11%bt01    = -7.689358913e+00_real_8
       cp_mgga_c_param%M08_M11%bt02    = -0.705326630e+00_real_8
       cp_mgga_c_param%M08_M11%bt03    = -0.600096421e+00_real_8
       cp_mgga_c_param%M08_M11%bt04    =  1.103332527e+01_real_8
       cp_mgga_c_param%M08_M11%bt05    =  5.861969337e+00_real_8
       cp_mgga_c_param%M08_M11%bt06    =  8.913865465e+00_real_8
       cp_mgga_c_param%M08_M11%bt07    =  5.745298760e+00_real_8
       cp_mgga_c_param%M08_M11%bt08    =  4.254880837e+00_real_8
       cp_mgga_c_param%M08_M11%bt09    =  0.000000000e+00_real_8
       cp_mgga_c_param%M08_M11%bt10    =  0.000000000e+00_real_8
       cp_mgga_c_param%M08_M11%bt11    =  0.000000000e+00_real_8
       cp_mgga_c_param%init            =  .true.
       !
       ! Hybrid-GGA
       ! Combinations
    CASE ("HYB_GGA_XC_B3LYP")
       !
       ! This is the exact definition of B3LYP and backwards-compatible with prior
       ! versions of CPMD (OLDCODE). Funnily enough, the energies from libxc are
       ! substantially different.
       ! 
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_LDA_X"
       cpfunc%scales( cpfunc%n_funcs ) = 0.08_real_8*cp_xc_env%scales( env_index )
       !
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_GGA_X_B88"
       cpfunc%scales( cpfunc%n_funcs ) = 0.72_real_8*cp_xc_env%scales( env_index )
       !
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_LDA_C_VWN"
       cpfunc%scales( cpfunc%n_funcs ) = 0.19_real_8*cp_xc_env%scales( env_index )
       !
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_GGA_C_LYP"
       cpfunc%scales( cpfunc%n_funcs ) = 0.81_real_8*cp_xc_env%scales( env_index )
    CASE ("HYB_GGA_XC_O3LYP")
       !
       ! Multiple definitions of \Delta OPTX possible; we stick to the one used in
       ! libxc, which can reproduce Handy's reference energies:
       !
       ! %E(LDA) = [0.9262-1.05151*0.81330]*100
       !
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_LDA_X"
       cpfunc%scales( cpfunc%n_funcs ) = 0.0710069170_real_8*cp_xc_env%scales( env_index )
       !
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_GGA_X_OPTX"
       cpfunc%scales( cpfunc%n_funcs ) = 0.81330_real_8*cp_xc_env%scales( env_index )
       !
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_LDA_C_VWN"
       cpfunc%scales( cpfunc%n_funcs ) = 0.19_real_8*cp_xc_env%scales( env_index )
       !
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_GGA_C_LYP"
       cpfunc%scales( cpfunc%n_funcs ) = 0.81_real_8*cp_xc_env%scales( env_index )
    CASE ("HYB_GGA_XC_CAM_B3LYP") 
       !
       ! Force attenuation here in case it is overwritten in the input
       !
       cpfunc%is_coulomb_attenuated    = .TRUE.
       !
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_GGA_X_B88"
       cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
       !
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_LDA_C_VWN"
       cpfunc%scales( cpfunc%n_funcs ) = 0.19_real_8*cp_xc_env%scales( env_index )
       !
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_GGA_C_LYP"
       cpfunc%scales( cpfunc%n_funcs ) = 0.81_real_8*cp_xc_env%scales( env_index )
    CASE ("HYB_GGA_XC_REVPBE0","HYB_GGA_XC_PBE0_R") 
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_GGA_X_REVPBE"
       cpfunc%scales( cpfunc%n_funcs ) = 0.75_real_8*cp_xc_env%scales( env_index )
       !
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_GGA_C_PBE"
       cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
       ! 
       ! Cases that were different in OLDCODE
    CASE ("HYB_GGA_XC_PBE0") 
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_GGA_X_PBE"
       cpfunc%scales( cpfunc%n_funcs ) = 0.75_real_8*cp_xc_env%scales( env_index )
       !
       IF (cpfunc%use_compatibility_mode) THEN
          ! Follows CPMD OLDCODE
          cpfunc%n_funcs                  = cpfunc%n_funcs + 1
          cpfunc%funcs( cpfunc%n_funcs )  = "CP_GGA_C_PBE_FLEX"
          cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
          cp_gga_c_param%pbe%lda_c        = "CP_LDA_C_VWN"
          cp_gga_c_param%pbe%pbe_beta     = 0.06672455060314922_real_8
          cp_gga_c_param%pbe%pbe_gamma    = 0.031090690869654895_real_8 
          cp_gga_c_param%init             = .true.
       ELSE
          cpfunc%n_funcs                  = cpfunc%n_funcs + 1
          cpfunc%funcs( cpfunc%n_funcs )  = "CP_GGA_C_PBE"
          cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
       END IF
    CASE ("HYB_GGA_XC_PBE0_SOL") 
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_GGA_X_PBE_SOL"
       cpfunc%scales( cpfunc%n_funcs ) = 0.75_real_8*cp_xc_env%scales( env_index )
       !
       IF (cpfunc%use_compatibility_mode) THEN
          CALL stopgm(procedureN,'Old definition n/a; Pade xc not coded',&
               __LINE__,__FILE__)
       ELSE
          cpfunc%n_funcs                  = cpfunc%n_funcs + 1
          cpfunc%funcs( cpfunc%n_funcs )  = "CP_GGA_C_PBE_SOL"
          cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
       END IF
       !
       ! Meta hybrids
    CASE ("HYB_MGGA_XC_M05")
       !
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_MGGA_X_M05_M06"
       cpfunc%scales( cpfunc%n_funcs ) = 0.72*cp_xc_env%scales( env_index )
       cp_mgga_x_param%M05_M06%add_vs98= .false.
       cp_mgga_x_param%M05_M06%at00    =    1.00000_real_8
       cp_mgga_x_param%M05_M06%at01    =    0.08151_real_8
       cp_mgga_x_param%M05_M06%at02    =   -0.43956_real_8
       cp_mgga_x_param%M05_M06%at03    =   -3.22422_real_8
       cp_mgga_x_param%M05_M06%at04    =    2.01819_real_8
       cp_mgga_x_param%M05_M06%at05    =    8.79431_real_8
       cp_mgga_x_param%M05_M06%at06    =   -0.00295_real_8
       cp_mgga_x_param%M05_M06%at07    =    9.82029_real_8
       cp_mgga_x_param%M05_M06%at08    =   -4.82351_real_8
       cp_mgga_x_param%M05_M06%at09    =  -48.17574_real_8
       cp_mgga_x_param%M05_M06%at10    =    3.64802_real_8
       cp_mgga_x_param%M05_M06%at11    =   34.02248_real_8
       cp_mgga_x_param%init            =  .true.
       !
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_MGGA_C_M05_M06"
       cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
       cp_mgga_c_param%M05_M06%add_vs98= .false.
       cp_mgga_c_param%M05_M06%sopp0   =   1.00000_real_8
       cp_mgga_c_param%M05_M06%sopp1   =   3.78569_real_8
       cp_mgga_c_param%M05_M06%sopp2   = -14.15261_real_8
       cp_mgga_c_param%M05_M06%sopp3   =  -7.46589_real_8
       cp_mgga_c_param%M05_M06%sopp4   =  17.94491_real_8
       cp_mgga_c_param%M05_M06%sss0    =   1.00000_real_8
       cp_mgga_c_param%M05_M06%sss1    =   3.77344_real_8
       cp_mgga_c_param%M05_M06%sss2    = -26.04463_real_8
       cp_mgga_c_param%M05_M06%sss3    =  30.69913_real_8
       cp_mgga_c_param%M05_M06%sss4    =  -9.22695_real_8
       cp_mgga_c_param%init            =  .true.
    CASE ("HYB_MGGA_XC_M05_2X")
       !
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_MGGA_X_M05_M06"
       cpfunc%scales( cpfunc%n_funcs ) = 0.44*cp_xc_env%scales( env_index )
       cp_mgga_x_param%M05_M06%add_vs98= .false.
       cp_mgga_x_param%M05_M06%at00    =   1.00000_real_8
       cp_mgga_x_param%M05_M06%at01    =  -0.56833_real_8
       cp_mgga_x_param%M05_M06%at02    =  -1.30057_real_8
       cp_mgga_x_param%M05_M06%at03    =   5.50070_real_8
       cp_mgga_x_param%M05_M06%at04    =   9.06402_real_8
       cp_mgga_x_param%M05_M06%at05    = -32.21075_real_8
       cp_mgga_x_param%M05_M06%at06    = -23.73298_real_8
       cp_mgga_x_param%M05_M06%at07    =  70.22996_real_8
       cp_mgga_x_param%M05_M06%at08    =  29.88614_real_8
       cp_mgga_x_param%M05_M06%at09    = -60.25778_real_8
       cp_mgga_x_param%M05_M06%at10    = -13.22205_real_8
       cp_mgga_x_param%M05_M06%at11    =  15.23694_real_8
       cp_mgga_x_param%init            =  .true.
       !
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_MGGA_C_M05_M06"
       cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
       cp_mgga_c_param%M05_M06%add_vs98= .false.
       cp_mgga_c_param%M05_M06%sopp0   =   1.00000_real_8
       cp_mgga_c_param%M05_M06%sopp1   =   1.09297_real_8
       cp_mgga_c_param%M05_M06%sopp2   =  -3.79171_real_8
       cp_mgga_c_param%M05_M06%sopp3   =   2.82810_real_8
       cp_mgga_c_param%M05_M06%sopp4   = -10.58909_real_8
       cp_mgga_c_param%M05_M06%sss0    =   1.00000_real_8
       cp_mgga_c_param%M05_M06%sss1    =  -3.05430_real_8
       cp_mgga_c_param%M05_M06%sss2    =   7.61854_real_8
       cp_mgga_c_param%M05_M06%sss3    =   1.47665_real_8
       cp_mgga_c_param%M05_M06%sss4    = -11.92365_real_8
       cp_mgga_c_param%init            =  .true.
    CASE ("HYB_MGGA_XC_M06")
       !
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_MGGA_X_M05_M06"
       ! Definition differs from M05; HFX scaling is included in
       ! at00
       cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
       cp_mgga_x_param%M05_M06%add_vs98=  .true.
       cp_mgga_x_param%M05_M06%at00    =  5.877943e-01_real_8
       cp_mgga_x_param%M05_M06%at01    = -1.371776e-01_real_8
       cp_mgga_x_param%M05_M06%at02    =  2.682367e-01_real_8
       cp_mgga_x_param%M05_M06%at03    = -2.515898e+00_real_8
       cp_mgga_x_param%M05_M06%at04    = -2.978892e+00_real_8
       cp_mgga_x_param%M05_M06%at05    =  8.710679e+00_real_8
       cp_mgga_x_param%M05_M06%at06    =  1.688195e+01_real_8
       cp_mgga_x_param%M05_M06%at07    = -4.489724e+00_real_8
       cp_mgga_x_param%M05_M06%at08    = -3.299983e+01_real_8
       cp_mgga_x_param%M05_M06%at09    = -1.449050e+01_real_8
       cp_mgga_x_param%M05_M06%at10    =  2.043747e+01_real_8
       cp_mgga_x_param%M05_M06%at11    =  1.256504e+01_real_8
       cp_mgga_x_param%VS98%r1         =  1.422057e-01_real_8*axlsda
       cp_mgga_x_param%VS98%r2         =  7.370319e-04_real_8*axlsda
       cp_mgga_x_param%VS98%r3         = -1.601373e-02_real_8*axlsda
       cp_mgga_x_param%VS98%r4         =  0.000000e+00_real_8
       cp_mgga_x_param%VS98%r5         =  0.000000e+00_real_8
       cp_mgga_x_param%VS98%r6         =  0.000000e+00_real_8
       cp_mgga_x_param%init            =  .true.
       !
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_MGGA_C_M05_M06"
       cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
       cp_mgga_c_param%M05_M06%add_vs98=  .true.
       cp_mgga_c_param%M05_M06%sopp0   =  3.741539e+00_real_8
       cp_mgga_c_param%M05_M06%sopp1   =  2.187098e+02_real_8
       cp_mgga_c_param%M05_M06%sopp2   = -4.531252e+02_real_8
       cp_mgga_c_param%M05_M06%sopp3   =  2.936479e+02_real_8
       cp_mgga_c_param%M05_M06%sopp4   = -6.287470e+01_real_8
       cp_mgga_c_param%M05_M06%sss0    =  5.094055e-01_real_8
       cp_mgga_c_param%M05_M06%sss1    = -1.491085e+00_real_8
       cp_mgga_c_param%M05_M06%sss2    =  1.723922e+01_real_8
       cp_mgga_c_param%M05_M06%sss3    = -3.859018e+01_real_8
       cp_mgga_c_param%M05_M06%sss4    =  2.845044e+01_real_8
       cp_mgga_c_param%VS98%r7         = -2.741539e+00_real_8
       cp_mgga_c_param%VS98%r8         = -6.720113e-01_real_8
       cp_mgga_c_param%VS98%r9         = -7.932688e-02_real_8
       cp_mgga_c_param%VS98%r10        =  1.918681e-03_real_8
       cp_mgga_c_param%VS98%r11        = -2.032902e-03_real_8
       cp_mgga_c_param%VS98%r12        =  0.000000e+00_real_8
       cp_mgga_c_param%VS98%r13        =  4.905945e-01_real_8
       cp_mgga_c_param%VS98%r14        = -1.437348e-01_real_8
       cp_mgga_c_param%VS98%r15        =  2.357824e-01_real_8
       cp_mgga_c_param%VS98%r16        =  1.871015e-03_real_8
       cp_mgga_c_param%VS98%r17        = -3.788963e-03_real_8
       cp_mgga_c_param%VS98%r18        =  0.000000e+00_real_8
       cp_mgga_c_param%init            =  .true.
    CASE ("HYB_MGGA_XC_M06_2X")
       !
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_MGGA_X_M05_M06"
       ! Definition differs from M05-2X; HFX scaling is included in
       ! at00
       cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
       cp_mgga_x_param%M05_M06%add_vs98=  .false.
       cp_mgga_x_param%M05_M06%at00    =  4.600000e-01_real_8
       cp_mgga_x_param%M05_M06%at01    = -2.206052e-01_real_8
       cp_mgga_x_param%M05_M06%at02    = -9.431788e-02_real_8
       cp_mgga_x_param%M05_M06%at03    =  2.164494e+00_real_8
       cp_mgga_x_param%M05_M06%at04    = -2.556466e+00_real_8
       cp_mgga_x_param%M05_M06%at05    = -1.422133e+01_real_8
       cp_mgga_x_param%M05_M06%at06    =  1.555044e+01_real_8
       cp_mgga_x_param%M05_M06%at07    =  3.598078e+01_real_8
       cp_mgga_x_param%M05_M06%at08    = -2.722754e+01_real_8
       cp_mgga_x_param%M05_M06%at09    = -3.924093e+01_real_8
       cp_mgga_x_param%M05_M06%at10    =  1.522808e+01_real_8
       cp_mgga_x_param%M05_M06%at11    =  1.522227e+01_real_8
       cp_mgga_x_param%init            =  .true.
       !
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_MGGA_C_M05_M06"
       cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
       cp_mgga_c_param%M05_M06%add_vs98=  .true.
       cp_mgga_c_param%M05_M06%sopp0   =  8.833596e-01_real_8
       cp_mgga_c_param%M05_M06%sopp1   =  3.357972e+01_real_8
       cp_mgga_c_param%M05_M06%sopp2   = -7.043548e+01_real_8
       cp_mgga_c_param%M05_M06%sopp3   =  4.978271e+01_real_8
       cp_mgga_c_param%M05_M06%sopp4   = -1.852891e+01_real_8
       cp_mgga_c_param%M05_M06%sss0    =  3.097855e-01_real_8
       cp_mgga_c_param%M05_M06%sss1    = -5.528642e+00_real_8
       cp_mgga_c_param%M05_M06%sss2    =  1.347420e+01_real_8
       cp_mgga_c_param%M05_M06%sss3    = -3.213623e+01_real_8
       cp_mgga_c_param%M05_M06%sss4    =  2.846742e+01_real_8
       cp_mgga_c_param%VS98%r7         =  1.166404e-01_real_8
       cp_mgga_c_param%VS98%r8         = -9.120847e-02_real_8
       cp_mgga_c_param%VS98%r9         = -6.726189e-02_real_8
       cp_mgga_c_param%VS98%r10        =  6.720580e-05_real_8
       cp_mgga_c_param%VS98%r11        =  8.448011e-04_real_8
       cp_mgga_c_param%VS98%r12        =  0.000000e+00_real_8
       cp_mgga_c_param%VS98%r13        =  6.902145e-01_real_8
       cp_mgga_c_param%VS98%r14        =  9.847204e-02_real_8
       cp_mgga_c_param%VS98%r15        =  2.214797e-01_real_8
       cp_mgga_c_param%VS98%r16        = -1.968264e-03_real_8
       cp_mgga_c_param%VS98%r17        = -6.775479e-03_real_8
       cp_mgga_c_param%VS98%r18        =  0.000000e+00_real_8
       cp_mgga_c_param%init            =  .true.
    CASE ("HYB_MGGA_XC_M06_HF")
       !
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_MGGA_X_M05_M06"
       cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
       cp_mgga_x_param%M05_M06%add_vs98=  .true.
       cp_mgga_x_param%M05_M06%at00    =  1.179732e-01_real_8
       cp_mgga_x_param%M05_M06%at01    = -1.066708e+00_real_8
       cp_mgga_x_param%M05_M06%at02    = -1.462405e-01_real_8
       cp_mgga_x_param%M05_M06%at03    =  7.481848e+00_real_8
       cp_mgga_x_param%M05_M06%at04    =  3.776679e+00_real_8
       cp_mgga_x_param%M05_M06%at05    = -4.436118e+01_real_8
       cp_mgga_x_param%M05_M06%at06    = -1.830962e+01_real_8
       cp_mgga_x_param%M05_M06%at07    =  1.003903e+02_real_8
       cp_mgga_x_param%M05_M06%at08    =  3.864360e+01_real_8
       cp_mgga_x_param%M05_M06%at09    = -9.806018e+01_real_8
       cp_mgga_x_param%M05_M06%at10    = -2.557716e+01_real_8
       cp_mgga_x_param%M05_M06%at11    =  3.590404e+01_real_8
       cp_mgga_x_param%VS98%r1         = -1.179732e-01_real_8*axlsda
       cp_mgga_x_param%VS98%r2         = -2.500000e-03_real_8*axlsda
       cp_mgga_x_param%VS98%r3         = -1.180065e-02_real_8*axlsda
       cp_mgga_x_param%VS98%r4         =  0.000000e+00_real_8
       cp_mgga_x_param%VS98%r5         =  0.000000e+00_real_8
       cp_mgga_x_param%VS98%r6         =  0.000000e+00_real_8
       cp_mgga_x_param%init            =  .true.
       !
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_MGGA_C_M05_M06"
       cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
       cp_mgga_c_param%M05_M06%add_vs98=  .true.
       cp_mgga_c_param%M05_M06%sopp0   =  1.674634e+00_real_8
       cp_mgga_c_param%M05_M06%sopp1   =  5.732017e+01_real_8
       cp_mgga_c_param%M05_M06%sopp2   =  5.955416e+01_real_8
       cp_mgga_c_param%M05_M06%sopp3   = -2.311007e+02_real_8
       cp_mgga_c_param%M05_M06%sopp4   =  1.255199e+02_real_8
       cp_mgga_c_param%M05_M06%sss0    =  1.023254e-01_real_8
       cp_mgga_c_param%M05_M06%sss1    = -2.453783e+00_real_8
       cp_mgga_c_param%M05_M06%sss2    =  2.913180e+01_real_8
       cp_mgga_c_param%M05_M06%sss3    = -3.494358e+01_real_8
       cp_mgga_c_param%M05_M06%sss4    =  2.315955e+01_real_8
       cp_mgga_c_param%VS98%r7         = -6.746338e-01_real_8
       cp_mgga_c_param%VS98%r8         = -1.534002e-01_real_8
       cp_mgga_c_param%VS98%r9         = -9.021521e-02_real_8
       cp_mgga_c_param%VS98%r10        = -1.292037e-03_real_8
       cp_mgga_c_param%VS98%r11        = -2.352983e-04_real_8
       cp_mgga_c_param%VS98%r12        =  0.000000e+00_real_8
       cp_mgga_c_param%VS98%r13        =  8.976746e-01_real_8
       cp_mgga_c_param%VS98%r14        = -2.345830e-01_real_8
       cp_mgga_c_param%VS98%r15        =  2.368173e-01_real_8
       cp_mgga_c_param%VS98%r16        = -9.913890e-04_real_8
       cp_mgga_c_param%VS98%r17        = -1.146165e-02_real_8
       cp_mgga_c_param%VS98%r18        =  0.000000e+00_real_8
       cp_mgga_c_param%init            =  .true.
    CASE ("HYB_MGGA_XC_M08_SO")
       !
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_MGGA_X_M08_M11"
       cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
       cp_mgga_x_param%M08_M11%lda_has_lc = .false.
       cp_mgga_x_param%M08_M11%at00    = -3.4888428e-01_real_8
       cp_mgga_x_param%M08_M11%at01    = -5.8157416e+00_real_8
       cp_mgga_x_param%M08_M11%at02    =  3.7550810e+01_real_8
       cp_mgga_x_param%M08_M11%at03    =  6.3727406e+01_real_8
       cp_mgga_x_param%M08_M11%at04    = -5.3742313e+01_real_8
       cp_mgga_x_param%M08_M11%at05    = -9.8595529e+01_real_8
       cp_mgga_x_param%M08_M11%at06    =  1.6282216e+01_real_8
       cp_mgga_x_param%M08_M11%at07    =  1.7513468e+01_real_8
       cp_mgga_x_param%M08_M11%at08    = -6.7627553e+00_real_8
       cp_mgga_x_param%M08_M11%at09    =  1.1106658e+01_real_8
       cp_mgga_x_param%M08_M11%at10    =  1.5663545e+00_real_8
       cp_mgga_x_param%M08_M11%at11    =  8.7603470e+00_real_8
       cp_mgga_x_param%M08_M11%bt00    =  7.8098428e-01_real_8
       cp_mgga_x_param%M08_M11%bt01    =  5.4538178e+00_real_8
       cp_mgga_x_param%M08_M11%bt02    = -3.7853348e+01_real_8
       cp_mgga_x_param%M08_M11%bt03    = -6.2295080e+01_real_8
       cp_mgga_x_param%M08_M11%bt04    =  4.6713254e+01_real_8
       cp_mgga_x_param%M08_M11%bt05    =  8.7321376e+01_real_8
       cp_mgga_x_param%M08_M11%bt06    =  1.6053446e+01_real_8
       cp_mgga_x_param%M08_M11%bt07    =  2.0126920e+01_real_8
       cp_mgga_x_param%M08_M11%bt08    = -4.0343695e+01_real_8
       cp_mgga_x_param%M08_M11%bt09    = -5.8577565e+01_real_8
       cp_mgga_x_param%M08_M11%bt10    =  2.0890272e+01_real_8
       cp_mgga_x_param%M08_M11%bt11    =  1.0946903e+01_real_8
       cp_mgga_x_param%init            =  .true.
       !
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_MGGA_C_M08_M11"
       cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
       cp_gga_c_param%pbe%lda_c        = "CP_LDA_C_PW"
       cp_gga_c_param%pbe%pbe_beta     = 0.06672455060314922_real_8
       cp_gga_c_param%pbe%pbe_gamma    = 0.031090690869654895_real_8
       cp_gga_c_param%init             = .true.
       cp_mgga_c_param%M08_M11%at00    =   1.0000000e+00_real_8
       cp_mgga_c_param%M08_M11%at01    =   0.0000000e+00_real_8
       cp_mgga_c_param%M08_M11%at02    =  -3.9980886e+00_real_8
       cp_mgga_c_param%M08_M11%at03    =   1.2982340e+01_real_8
       cp_mgga_c_param%M08_M11%at04    =   1.0117507e+02_real_8
       cp_mgga_c_param%M08_M11%at05    =  -8.9541984e+01_real_8
       cp_mgga_c_param%M08_M11%at06    =  -3.5640242e+02_real_8
       cp_mgga_c_param%M08_M11%at07    =   2.0698803e+02_real_8
       cp_mgga_c_param%M08_M11%at08    =   4.6037780e+02_real_8
       cp_mgga_c_param%M08_M11%at09    =  -2.4510559e+02_real_8
       cp_mgga_c_param%M08_M11%at10    =  -1.9638425e+02_real_8
       cp_mgga_c_param%M08_M11%at11    =   1.1881459e+02_real_8
       cp_mgga_c_param%M08_M11%bt00    =   1.0000000e+00_real_8
       cp_mgga_c_param%M08_M11%bt01    =  -4.4117403e+00_real_8
       cp_mgga_c_param%M08_M11%bt02    =  -6.4128622e+00_real_8
       cp_mgga_c_param%M08_M11%bt03    =   4.7583635e+01_real_8
       cp_mgga_c_param%M08_M11%bt04    =   1.8630053e+02_real_8
       cp_mgga_c_param%M08_M11%bt05    =  -1.2800784e+02_real_8
       cp_mgga_c_param%M08_M11%bt06    =  -5.5385258e+02_real_8
       cp_mgga_c_param%M08_M11%bt07    =   1.3873727e+02_real_8
       cp_mgga_c_param%M08_M11%bt08    =   4.1646537e+02_real_8
       cp_mgga_c_param%M08_M11%bt09    =  -2.6626577e+02_real_8
       cp_mgga_c_param%M08_M11%bt10    =   5.6676300e+01_real_8
       cp_mgga_c_param%M08_M11%bt11    =   3.1673746e+02_real_8
       cp_mgga_c_param%init            =  .true.
    CASE ("HYB_MGGA_XC_M08_HX")
       !
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_MGGA_X_M08_M11"
       cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
       cp_mgga_x_param%M08_M11%lda_has_lc = .false.
       cp_mgga_x_param%M08_M11%at00    =  1.3340172e+00_real_8
       cp_mgga_x_param%M08_M11%at01    = -9.4751087e+00_real_8
       cp_mgga_x_param%M08_M11%at02    = -1.2541893e+01_real_8
       cp_mgga_x_param%M08_M11%at03    =  9.1369974e+00_real_8
       cp_mgga_x_param%M08_M11%at04    =  3.4717204e+01_real_8
       cp_mgga_x_param%M08_M11%at05    =  5.8831807e+01_real_8
       cp_mgga_x_param%M08_M11%at06    =  7.1369574e+01_real_8
       cp_mgga_x_param%M08_M11%at07    =  2.3312961e+01_real_8
       cp_mgga_x_param%M08_M11%at08    =  4.8314679e+00_real_8
       cp_mgga_x_param%M08_M11%at09    = -6.5044167e+00_real_8
       cp_mgga_x_param%M08_M11%at10    = -1.4058265e+01_real_8
       cp_mgga_x_param%M08_M11%at11    =  1.2880570e+01_real_8
       cp_mgga_x_param%M08_M11%bt00    = -8.5631823e-01_real_8
       cp_mgga_x_param%M08_M11%bt01    =  9.2810354e+00_real_8
       cp_mgga_x_param%M08_M11%bt02    =  1.2260749e+01_real_8
       cp_mgga_x_param%M08_M11%bt03    = -5.5189665e+00_real_8
       cp_mgga_x_param%M08_M11%bt04    = -3.5534989e+01_real_8
       cp_mgga_x_param%M08_M11%bt05    = -8.2049996e+01_real_8
       cp_mgga_x_param%M08_M11%bt06    = -6.8586558e+01_real_8
       cp_mgga_x_param%M08_M11%bt07    =  3.6085694e+01_real_8
       cp_mgga_x_param%M08_M11%bt08    = -9.3740983e+00_real_8
       cp_mgga_x_param%M08_M11%bt09    = -5.9731688e+01_real_8
       cp_mgga_x_param%M08_M11%bt10    =  1.6587868e+01_real_8
       cp_mgga_x_param%M08_M11%bt11    =  1.3993203e+01_real_8
       cp_mgga_x_param%init            =  .true.
       !
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_MGGA_C_M08_M11"
       cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
       cp_gga_c_param%pbe%lda_c        = "CP_LDA_C_PW"
       cp_gga_c_param%pbe%pbe_beta     = 0.06672455060314922_real_8
       cp_gga_c_param%pbe%pbe_gamma    = 0.031090690869654895_real_8
       cp_gga_c_param%init             = .true.
       cp_mgga_c_param%M08_M11%at00    =   1.0000000e+00_real_8
       cp_mgga_c_param%M08_M11%at01    =  -4.0661387e-01_real_8
       cp_mgga_c_param%M08_M11%at02    =  -3.3232530e+00_real_8
       cp_mgga_c_param%M08_M11%at03    =   1.5540980e+00_real_8
       cp_mgga_c_param%M08_M11%at04    =   4.4248033e+01_real_8
       cp_mgga_c_param%M08_M11%at05    =  -8.4351930e+01_real_8
       cp_mgga_c_param%M08_M11%at06    =  -1.1955581e+02_real_8
       cp_mgga_c_param%M08_M11%at07    =   3.9147081e+02_real_8
       cp_mgga_c_param%M08_M11%at08    =   1.8363851e+02_real_8
       cp_mgga_c_param%M08_M11%at09    =  -6.3268223e+02_real_8
       cp_mgga_c_param%M08_M11%at10    =  -1.1297403e+02_real_8
       cp_mgga_c_param%M08_M11%at11    =   3.3629312e+02_real_8
       cp_mgga_c_param%M08_M11%bt00    =   1.3812334e+00_real_8
       cp_mgga_c_param%M08_M11%bt01    =  -2.4683806e+00_real_8
       cp_mgga_c_param%M08_M11%bt02    =  -1.1901501e+01_real_8
       cp_mgga_c_param%M08_M11%bt03    =  -5.4112667e+01_real_8
       cp_mgga_c_param%M08_M11%bt04    =   1.0055846e+01_real_8
       cp_mgga_c_param%M08_M11%bt05    =   1.4800687e+02_real_8
       cp_mgga_c_param%M08_M11%bt06    =   1.1561420e+02_real_8
       cp_mgga_c_param%M08_M11%bt07    =   2.5591815e+02_real_8
       cp_mgga_c_param%M08_M11%bt08    =   2.1320772e+02_real_8
       cp_mgga_c_param%M08_M11%bt09    =  -4.8412067e+02_real_8
       cp_mgga_c_param%M08_M11%bt10    =  -4.3430813e+02_real_8
       cp_mgga_c_param%M08_M11%bt11    =   5.6627964e+01_real_8
       cp_mgga_c_param%init            =  .true.
    CASE ("HYB_MGGA_XC_M11")
       !
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_MGGA_X_M08_M11"
       cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
       cp_mgga_x_param%M08_M11%lda_has_lc = .true.
       cp_mgga_x_param%M08_M11%at00    = -0.18399900e+00_real_8
       cp_mgga_x_param%M08_M11%at01    = -1.39046703e+01_real_8
       cp_mgga_x_param%M08_M11%at02    =  1.18206837e+01_real_8
       cp_mgga_x_param%M08_M11%at03    =  3.10098465e+01_real_8
       cp_mgga_x_param%M08_M11%at04    = -5.19625696e+01_real_8
       cp_mgga_x_param%M08_M11%at05    =  1.55750312e+01_real_8
       cp_mgga_x_param%M08_M11%at06    = -6.94775730e+00_real_8
       cp_mgga_x_param%M08_M11%at07    = -1.58465014e+02_real_8
       cp_mgga_x_param%M08_M11%at08    = -1.48447565e+00_real_8
       cp_mgga_x_param%M08_M11%at09    =  5.51042124e+01_real_8
       cp_mgga_x_param%M08_M11%at10    = -1.34714184e+01_real_8
       cp_mgga_x_param%M08_M11%at11    =  0.00000000e+00_real_8
       cp_mgga_x_param%M08_M11%bt00    =  0.75599900e+00_real_8
       cp_mgga_x_param%M08_M11%bt01    =  1.37137944e+01_real_8
       cp_mgga_x_param%M08_M11%bt02    = -1.27998304e+01_real_8
       cp_mgga_x_param%M08_M11%bt03    = -2.93428814e+01_real_8
       cp_mgga_x_param%M08_M11%bt04    =  5.91075674e+01_real_8
       cp_mgga_x_param%M08_M11%bt05    = -2.27604866e+01_real_8
       cp_mgga_x_param%M08_M11%bt06    = -1.02769340e+01_real_8
       cp_mgga_x_param%M08_M11%bt07    =  1.64752731e+02_real_8
       cp_mgga_x_param%M08_M11%bt08    =  1.85349258e+01_real_8
       cp_mgga_x_param%M08_M11%bt09    = -5.56825639e+01_real_8
       cp_mgga_x_param%M08_M11%bt10    =  7.47980859e+00_real_8
       cp_mgga_x_param%M08_M11%bt11    =  0.00000000e+00_real_8
       cp_mgga_x_param%init            =  .true.
       !
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_MGGA_C_M08_M11"
       cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
       cp_gga_c_param%pbe%lda_c        = "CP_LDA_C_PW"
       cp_gga_c_param%pbe%pbe_beta     = 0.06672455060314922_real_8
       cp_gga_c_param%pbe%pbe_gamma    = 0.031090690869654895_real_8
       cp_gga_c_param%init             = .true. 
       cp_mgga_c_param%M08_M11%at00    =   1.0000000e+00_real_8
       cp_mgga_c_param%M08_M11%at01    =   0.0000000e+00_real_8
       cp_mgga_c_param%M08_M11%at02    =  -3.8933250e+00_real_8
       cp_mgga_c_param%M08_M11%at03    =  -2.1688455e+00_real_8
       cp_mgga_c_param%M08_M11%at04    =   9.3497200e+00_real_8
       cp_mgga_c_param%M08_M11%at05    =  -1.9845140e+01_real_8
       cp_mgga_c_param%M08_M11%at06    =   2.3455253e+00_real_8
       cp_mgga_c_param%M08_M11%at07    =   7.9246513e+01_real_8
       cp_mgga_c_param%M08_M11%at08    =   9.6042757e+00_real_8
       cp_mgga_c_param%M08_M11%at09    =  -6.7856719e+01_real_8
       cp_mgga_c_param%M08_M11%at10    =  -9.1841067e+00_real_8
       cp_mgga_c_param%M08_M11%at11    =   0.0000000e+00_real_8
       cp_mgga_c_param%M08_M11%bt00    =   7.2239798e-01_real_8
       cp_mgga_c_param%M08_M11%bt01    =   4.3730564e-01_real_8
       cp_mgga_c_param%M08_M11%bt02    =  -1.6088809e+01_real_8
       cp_mgga_c_param%M08_M11%bt03    =  -6.5542437e+01_real_8
       cp_mgga_c_param%M08_M11%bt04    =   3.2057230e+01_real_8
       cp_mgga_c_param%M08_M11%bt05    =   1.8617888e+02_real_8
       cp_mgga_c_param%M08_M11%bt06    =   2.0483468e+01_real_8
       cp_mgga_c_param%M08_M11%bt07    =  -7.0853739e+01_real_8
       cp_mgga_c_param%M08_M11%bt08    =   4.4483915e+01_real_8
       cp_mgga_c_param%M08_M11%bt09    =  -9.4484747e+01_real_8
       cp_mgga_c_param%M08_M11%bt10    =  -1.1459868e+02_real_8
       cp_mgga_c_param%M08_M11%bt11    =   0.0000000e+00_real_8
       cp_mgga_c_param%init            =  .true.
    CASE ("HYB_MGGA_XC_MN12_SX")
       !
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_MGGA_X_MN12_MN15"
       cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
       cp_mgga_x_param%MN12_MN15%cc000 =  5.226556e-01_real_8
       cp_mgga_x_param%MN12_MN15%cc001 = -2.681208e-01_real_8
       cp_mgga_x_param%MN12_MN15%cc002 = -4.670705e+00_real_8
       cp_mgga_x_param%MN12_MN15%cc003 =  3.067320e+00_real_8
       cp_mgga_x_param%MN12_MN15%cc004 =  4.095370e+00_real_8
       cp_mgga_x_param%MN12_MN15%cc005 =  2.653023e+00_real_8
       cp_mgga_x_param%MN12_MN15%cc010 =  5.165969e-01_real_8
       cp_mgga_x_param%MN12_MN15%cc011 = -2.035442e+01_real_8
       cp_mgga_x_param%MN12_MN15%cc012 = -9.946472e+00_real_8
       cp_mgga_x_param%MN12_MN15%cc013 =  2.938637e+00_real_8
       cp_mgga_x_param%MN12_MN15%cc014 =  1.131100e+01_real_8
       cp_mgga_x_param%MN12_MN15%cc020 =  4.752452e+00_real_8
       cp_mgga_x_param%MN12_MN15%cc021 = -3.061331e+00_real_8
       cp_mgga_x_param%MN12_MN15%cc022 = -2.523173e+01_real_8
       cp_mgga_x_param%MN12_MN15%cc023 =  1.710903e+01_real_8
       cp_mgga_x_param%MN12_MN15%cc030 = -2.357480e+01_real_8
       cp_mgga_x_param%MN12_MN15%cc031 = -2.727754e+01_real_8
       cp_mgga_x_param%MN12_MN15%cc032 =  1.603291e+01_real_8
       cp_mgga_x_param%MN12_MN15%cc100 =  1.842503e+00_real_8
       cp_mgga_x_param%MN12_MN15%cc101 =  1.927120e+00_real_8
       cp_mgga_x_param%MN12_MN15%cc102 =  1.107987e+01_real_8
       cp_mgga_x_param%MN12_MN15%cc103 = -1.182087e+01_real_8
       cp_mgga_x_param%MN12_MN15%cc104 = -1.117768e+01_real_8
       cp_mgga_x_param%MN12_MN15%cc110 = -5.821000e+00_real_8
       cp_mgga_x_param%MN12_MN15%cc111 =  2.266545e+01_real_8
       cp_mgga_x_param%MN12_MN15%cc112 =  8.246708e+00_real_8
       cp_mgga_x_param%MN12_MN15%cc113 = -4.778364e+00_real_8
       cp_mgga_x_param%MN12_MN15%cc120 =  5.329122e-01_real_8
       cp_mgga_x_param%MN12_MN15%cc121 = -6.666755e+00_real_8
       cp_mgga_x_param%MN12_MN15%cc122 =  1.671429e+00_real_8
       cp_mgga_x_param%MN12_MN15%cc200 = -3.311409e+00_real_8
       cp_mgga_x_param%MN12_MN15%cc201 =  3.415913e-01_real_8
       cp_mgga_x_param%MN12_MN15%cc202 = -6.413076e+00_real_8
       cp_mgga_x_param%MN12_MN15%cc203 =  1.038584e+01_real_8
       cp_mgga_x_param%MN12_MN15%cc210 =  9.026277e+00_real_8
       cp_mgga_x_param%MN12_MN15%cc211 =  1.929689e+01_real_8
       cp_mgga_x_param%MN12_MN15%cc212 =  2.669232e+01_real_8
       cp_mgga_x_param%MN12_MN15%cc300 =  1.517278e+00_real_8
       cp_mgga_x_param%MN12_MN15%cc301 = -3.442503e+00_real_8
       cp_mgga_x_param%MN12_MN15%cc302 =  1.100161e+00_real_8
       cp_mgga_x_param%init            =  .true.
       !
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_MGGA_C_M08_M11"
       cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
       cp_gga_c_param%pbe%lda_c        = "CP_LDA_C_PW"
       cp_gga_c_param%pbe%pbe_beta     = 0.06672455060314922_real_8
       cp_gga_c_param%pbe%pbe_gamma    = 0.031090690869654895_real_8
       cp_gga_c_param%init             = .true. 
       cp_mgga_c_param%M08_M11%at00    =  7.171161e-01_real_8
       cp_mgga_c_param%M08_M11%at01    = -2.380914e+00_real_8
       cp_mgga_c_param%M08_M11%at02    =  5.793565e+00_real_8
       cp_mgga_c_param%M08_M11%at03    = -1.243624e+00_real_8
       cp_mgga_c_param%M08_M11%at04    =  1.364920e+01_real_8
       cp_mgga_c_param%M08_M11%at05    = -2.110812e+01_real_8
       cp_mgga_c_param%M08_M11%at06    = -1.598767e+01_real_8
       cp_mgga_c_param%M08_M11%at07    =  1.429208e+01_real_8
       cp_mgga_c_param%M08_M11%at08    =  6.149191e+00_real_8
       cp_mgga_c_param%M08_M11%at09    =  0.000000e+00_real_8
       cp_mgga_c_param%M08_M11%at10    =  0.000000e+00_real_8
       cp_mgga_c_param%M08_M11%at11    =  0.000000e+00_real_8
       cp_mgga_c_param%M08_M11%bt00    =  4.663699e-01_real_8
       cp_mgga_c_param%M08_M11%bt01    = -9.110685e+00_real_8
       cp_mgga_c_param%M08_M11%bt02    =  8.705051e+00_real_8
       cp_mgga_c_param%M08_M11%bt03    = -1.813949e+00_real_8
       cp_mgga_c_param%M08_M11%bt04    = -4.147211e-01_real_8
       cp_mgga_c_param%M08_M11%bt05    = -1.021527e+01_real_8
       cp_mgga_c_param%M08_M11%bt06    =  8.240270e-01_real_8
       cp_mgga_c_param%M08_M11%bt07    =  4.993815e+00_real_8
       cp_mgga_c_param%M08_M11%bt08    = -2.563930e+01_real_8
       cp_mgga_c_param%M08_M11%bt09    =  0.000000e+00_real_8
       cp_mgga_c_param%M08_M11%bt10    =  0.000000e+00_real_8
       cp_mgga_c_param%M08_M11%bt11    =  0.000000e+00_real_8
       cp_mgga_c_param%init            =  .true.
       !
       ! In case there is nothing special to be done...
    CASE default
       cpfunc%n_funcs                  = cpfunc%n_funcs + 1
       cpfunc%funcs( cpfunc%n_funcs )  = "CP_"//TRIM(ADJUSTL(cp_xc_env%funcs( env_index )))
       cpfunc%scales( cpfunc%n_funcs ) = cp_xc_env%scales( env_index )
    END SELECT

  END SUBROUTINE cpfunc_func_init
  ! ==================================================================
  SUBROUTINE hfx_func_init( hfx, cp_xc_env )
    ! ==--------------------------------------------------------------==
    ! Initialisation of hfx parameters
    !                               17.08.2017 M. P. Bircher @ LCBC/EPFL

    CLASS(cp_hfx_t), INTENT(INOUT)           :: hfx
    TYPE(cp_xc_env_t), INTENT(IN)            :: cp_xc_env

    CHARACTER(len=*), PARAMETER :: procedureN = 'hfx_func_init'

    SELECT CASE(cp_xc_env%hfx_operator)
    CASE("CAM")
       hfx%is_coulomb_attenuated = .TRUE.
       hfx%scale = 1.0_real_8
       hfx%alpha = cp_xc_env%hfx_constant
       hfx%beta  = cp_xc_env%hfx_attenuated
       hfx%gamma = cp_xc_env%hfx_screening
       IF ( hfx%alpha + hfx%beta > 1.0_real_8 ) CALL stopgm(procedureN, &
          'CAM alpha and beta sum up to more than 100%', &
          __LINE__, __FILE__ )
    CASE("erfc","ERFC")
       hfx%is_screened_erfc = .TRUE.
       hfx%scale = cp_xc_env%hfx_constant
       hfx%gamma = cp_xc_env%hfx_screening
    CASE("exp","EXP")
       hfx%is_screened_exp = .TRUE.
       hfx%scale = cp_xc_env%hfx_constant
       hfx%gamma = cp_xc_env%hfx_screening
    CASE("Ashcroft","ashcroft","ASHCROFT")
       hfx%is_screened_ashcroft = .TRUE.
       hfx%scale = cp_xc_env%hfx_constant
       hfx%gamma = cp_xc_env%hfx_screening
    CASE DEFAULT
       hfx%scale = cp_xc_env%hfx_constant
       hfx%gamma = cp_xc_env%hfx_screening
    END SELECT

    ! ==--------------------------------------------------------------==
  END SUBROUTINE hfx_func_init
  ! ==================================================================

  ! ==================================================================
  ! Environment (input) parameters
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE cp_xc_env_report(this, salpha, gc_epsilon, write_to_unit)
    ! ==--------------------------------------------------------------==
    ! Writes chosen functionals to stdout BEFORE internal conversion
    !                               07.10.2016 M. P. Bircher @ LCBC/EPFL


    CLASS(cp_xc_env_t), INTENT(in)           :: this
    REAL(real_8), INTENT(in)                 :: gc_epsilon, salpha
    INTEGER, INTENT(in)                      :: write_to_unit

    CHARACTER(len=*), PARAMETER :: procedureN = 'cp_xc_env_report'

    INTEGER                                  :: i_func

    WRITE(write_to_unit,*) 
    WRITE(write_to_unit,'(1x,A)') TRIM(ADJUSTL(this%set_name))//':'
    WRITE(write_to_unit,*) 
    !
    ! Minnesota functional? Print a warning
    IF ( is_a_minnesota_functional( this%funcs(:) )) THEN
          WRITE(write_to_unit,'(1X,A)') 'MINNESOTA FUNCTIONALS: TO AVOID OSCILLATORY BEHAVIOUR, IT IS IM-'
          WRITE(write_to_unit,'(1X,A)') 'PERATIVE THAT THE DUAL BE SET TO A VALUE LARGER THAN THE DEFAULT'
          WRITE(write_to_unit,'(1X,A)') 'PLEASE REFER TO AND CITE:'
          WRITE(write_to_unit,'(1X,A)') '              M.P. Bircher, P. Lopez-Tarifa and U. Rothlisberger'
          WRITE(write_to_unit,'(1X,A)') '               J. Chem. Theory Comput., 2019, 15 (1), pp 557-571'
          WRITE(write_to_unit,'(1X,A)') '                                   DOI: 10.1021/acs.jctc.8b00897'
          WRITE(write_to_unit,'(A)')
    ENDIF
    !
    ! (Semi-)local functional
    DO i_func=1, this%n_funcs
       IF (this%via_libxc( i_func )) THEN
          WRITE(write_to_unit,'(4x,F5.1,A1,1x,A5,1x,A)') 100_real_8*this%scales( i_func ) , '%', &
               'LIBXC', TRIM(ADJUSTL(this%funcs( i_func )))
       ELSE
          WRITE(write_to_unit,'(4x,F5.1,A1,1x,A2,4x,A)') 100_real_8*this%scales( i_func ) , '%', &
               'CP', TRIM(ADJUSTL(this%funcs( i_func )))
       END IF
    END DO
    !
    ! Compatibility mode
    IF (this%use_compatibility_mode) &
         WRITE(write_to_unit,'(11x,A,1x,F5.3)') 'CP FUNCTIONALS FOLLOW DEPRECATED DEFINITIONS (OLDCODE)'
    !
    ! Non-standard S_alpha
    IF (salpha /= 2.0_real_8/3.0_real_8) &
         WRITE(write_to_unit,'(11x,A,1x,F5.3)') 'CP USES CUSTOMISED SLATER ALPHA:', salpha
    !
    ! Screening of (semi-)local part
    IF (this%get_CAM_GGA) THEN
       WRITE(write_to_unit,'(11x,A)')&
            'CP EXCHANGE FUNCTIONALS ARE COULOMB-ATTENUATED (CAM)'
       WRITE(write_to_unit,'(11x,A10,F8.5,1x,A8,F8.5,1x,A7,F8.5)')&
            'MU [BOHR]:', this%hfx_screening, &
            '; ALPHA:', this%hfx_constant, &
            '; BETA:', this%hfx_attenuated
    END IF
    !
    ! GC-Cutoff
    IF ( ANY(INDEX(this%funcs(:),'GGA') /= 0) ) THEN
       WRITE(write_to_unit,'(11x,A,3x,ES12.4)') 'GC-CUTOFF VALUE FOR GRADIENT SMOOTHING:', gc_epsilon
    END IF
    !
    ! HFX
    IF ( this%get_hfx ) THEN
       WRITE(write_to_unit,*) 

       SELECT CASE( this%hfx_operator)
       CASE("exp","EXP")
          WRITE(write_to_unit,'(4x,F5.1,A1,1x,A,3x,A)') 100_real_8*this%hfx_constant , '%', 'HFX', 'EXACT EXCHANGE'
          WRITE(write_to_unit,'(11x,A)') 'HFX IS SCREENED USING EXP(-MU*R)/R'
          WRITE(write_to_unit,'(11x,A13,F8.5)') 'MU [BOHR^-1]:', this%hfx_screening

       CASE("erfc","ERFC")
          WRITE(write_to_unit,'(4x,F5.1,A1,1x,A,3x,A)') 100_real_8*this%hfx_constant , '%', 'HFX', 'EXACT EXCHANGE'
          WRITE(write_to_unit,'(11x,A)') 'HFX IS SCREENED USING ERFC(MU*R)/R'
          WRITE(write_to_unit,'(11x,A13,F8.5)') 'MU [BOHR^-1]:', this%hfx_screening

       CASE("Ashcroft","ashcroft","ASHCROFT") 
          WRITE(write_to_unit,'(4x,F5.1,A1,1x,A,3x,A)') 100_real_8*this%hfx_constant , '%', 'HFX', 'EXACT EXCHANGE'
          WRITE(write_to_unit,'(11x,A)') 'HFX IS SCREENED USING ASHCROFT'

       CASE("CAM")
          WRITE(write_to_unit,'(4x,F5.1,A1,1x,A,4x,A)') 100_real_8*this%hfx_constant, '%', &
                                         'SR', 'EXACT EXCHANGE'
          WRITE(write_to_unit,'(4x,F5.1,A1,1x,A,4x,A)') 100_real_8*(this%hfx_constant + &
                                         this%hfx_attenuated), '%', 'LR', 'EXACT EXCHANGE'
          WRITE(write_to_unit,'(11x,A)')&
               'HFX IS COULOMB-ATTENUATED (CAM)'
          WRITE(write_to_unit,'(11x,A10,F8.5,1x,A8,F8.5,1x,A7,F8.5)')&
               'MU [BOHR]:', this%hfx_screening, &
               '; ALPHA:', this%hfx_constant, &
               '; BETA:', this%hfx_attenuated
          WRITE(write_to_unit,'(A)')
          WRITE(write_to_unit,'(1X,A)') 'FOR THE RECIPROCAL SPACE IMPLEMENTATION OF CAM, PLEASE REFER TO:'
          WRITE(write_to_unit,'(1X,A)') '                               M.P. Bircher and U. Rothlisberger'
          WRITE(write_to_unit,'(1X,A)') '             J. Chem. Theory Comput., 2018, 14 (6), pp 3184–3195'
          WRITE(write_to_unit,'(1X,A)') '                                   DOI: 10.1021/acs.jctc.8b00069'
          WRITE(write_to_unit,'(A)')

       CASE DEFAULT
          WRITE(write_to_unit,'(4x,F5.1,A1,1x,A,3x,A)') 100_real_8*this%hfx_constant , '%', 'HFX', 'EXACT EXCHANGE'
       END SELECT

       IF (this%overwrite_hfx) WRITE(write_to_unit,'(11x,A)') 'HFX PARAMETERS ARE NON-DEFAULT'

    END IF

    CONTAINS

    ! ==--------------------------------------------------------------==
    LOGICAL FUNCTION is_a_minnesota_functional( funcs )

        CHARACTER(len=*), DIMENSION(:), INTENT(in) :: funcs
        INTEGER                                    :: i

        is_a_minnesota_functional = .FALSE.

        DO i=1,size(funcs)
           IF (  trim(adjustl( funcs(i) )) == 'HYB_MGGA_XC_M05'    &
            .OR. trim(adjustl( funcs(i) )) == 'HYB_MGGA_XC_M05_2X' &
            .OR. trim(adjustl( funcs(i) )) == 'HYB_MGGA_XC_M06'    &
            .OR. trim(adjustl( funcs(i) )) == 'HYB_MGGA_XC_M06_2X' &
            .OR. trim(adjustl( funcs(i) )) == 'HYB_MGGA_XC_M06_HF' &
            .OR. trim(adjustl( funcs(i) )) == 'MGGA_XC_M06_L'      &
            .OR. trim(adjustl( funcs(i) )) == 'MGGA_XC_REVM06_L'   &
            .OR. trim(adjustl( funcs(i) )) == 'HYB_MGGA_XC_M08_HX' &
            .OR. trim(adjustl( funcs(i) )) == 'HYB_MGGA_XC_M08_SO' &
            .OR. trim(adjustl( funcs(i) )) == 'HYB_MGGA_XC_M11'    &
            .OR. trim(adjustl( funcs(i) )) == 'MGGA_XC_M11_L'      &
            .OR. trim(adjustl( funcs(i) )) == 'HYB_MGGA_XC_MN12'   &
            .OR. trim(adjustl( funcs(i) )) == 'MGGA_XC_MN12_L'     &
            .OR. trim(adjustl( funcs(i) )) == 'MGGA_XC_MN15_L' ) THEN
              is_a_minnesota_functional = .TRUE.
              EXIT
           ENDIF
        ENDDO

    END FUNCTION
    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_xc_env_report
  ! ==================================================================
  SUBROUTINE cp_xc_env_set_hybrid_defaults(this,default_divergence)
    ! ==--------------------------------------------------------------==
    ! Library with hybrid parameters. This is separated from
    ! cpfunc_func_init for outputting/compatiblity with oldcode.
    ! Can be changed/refactored once the old drivers are deprecated.
    !                               31.08.2017 M. P. Bircher @ LCBC/EPFL


    CLASS(cp_xc_env_t), INTENT(inout)        :: this
    INTEGER                                  :: j
    LOGICAL, OPTIONAL, INTENT(inout)         :: default_divergence

    IF (.not. this%overwrite_hfx) THEN
       DO j=1,this%n_funcs
          SELECT CASE(this%funcs(j))
          CASE('HYB_GGA_XC_B3LYP','HYB_GGA_XC_B3LYP5')
             this%get_hfx      = .TRUE. 
             this%hfx_constant = 0.2_real_8
          CASE('HYB_GGA_XC_PBEH')
             this%get_hfx      = .TRUE. 
             this%hfx_constant = 0.25_real_8
          CASE('HYB_GGA_XC_CAM_B3LYP')
             this%get_hfx        = .TRUE.
             this%get_CAM_GGA    = .TRUE.
             this%hfx_operator   = 'CAM'
             this%hfx_constant   = 0.190_real_8
             this%hfx_attenuated = 0.460_real_8
             this%hfx_screening  = 0.330_real_8
             IF (present(default_divergence)) default_divergence  = .TRUE.
             IF (this%via_libxc(j) .and. present(default_divergence)) default_divergence = .FALSE.
          CASE('HYB_GGA_XC_O3LYP')
             this%get_hfx      = .TRUE.
             this%hfx_constant = 0.11610_real_8
          CASE('HYB_GGA_XC_PBE0', &
               'HYB_GGA_XC_REVPBE0', 'HYB_GGA_XC_PBE0_R', &
               'HYB_GGA_XC_PBE0_SOL')
             this%get_hfx      = .TRUE.
             this%hfx_constant = 0.25_real_8
          CASE('HYB_MGGA_XC_M05')
             this%get_hfx      = .TRUE.
             this%hfx_constant = 0.28_real_8
          CASE('HYB_MGGA_XC_M05_2X')
             this%get_hfx      = .TRUE.
             this%hfx_constant = 0.56_real_8
          CASE('HYB_MGGA_XC_M06')
             this%get_hfx      = .TRUE.
             this%hfx_constant = 0.27_real_8
          CASE('HYB_MGGA_XC_M06_2X')
             this%get_hfx      = .TRUE.
             this%hfx_constant = 0.54_real_8
          CASE('HYB_MGGA_XC_M06_HF')
             this%get_hfx      = .TRUE.
             this%hfx_constant = 1.00_real_8
          CASE('HYB_MGGA_XC_M08_SO')
             this%get_hfx      = .TRUE.
             this%hfx_constant = 0.5679_real_8
          CASE('HYB_MGGA_XC_M08_HX')
             this%get_hfx      = .TRUE.
             this%hfx_constant = 0.5223_real_8
          CASE('MGGA_XC_M11_L')
             ! not a true hybrid, but with range separation parameter
             this%get_CAM_GGA    = .false. ! unfortunately hard-coded
             this%hfx_screening  = 0.250_real_8 ! hard-coded in MGGA part, too!
          CASE('HYB_MGGA_XC_M11')
             IF (present(default_divergence)) default_divergence = .TRUE.
             this%get_hfx        = .TRUE. ! unfortunately hard-coded
             this%get_CAM_GGA    = .FALSE. ! unfortunately hard-coded
             this%hfx_operator   = 'CAM'
             this%hfx_screening  = 0.250_real_8 ! hard-coded in MGGA part, too!
             this%hfx_constant   = 0.428_real_8
             this%hfx_attenuated = 0.572_real_8
          CASE('HYB_MGGA_XC_MN12_SX')
             this%get_hfx       = .TRUE.
             this%hfx_operator  = 'erfc'
             this%hfx_screening = 0.110_real_8 ! hard-coded in MGGA part, too!
             this%hfx_constant  = 0.25_real_8 
          END SELECT
       ENDDO
    ENDIF

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_xc_env_set_hybrid_defaults
  ! ==================================================================
  PURE FUNCTION cp_xc_env_get_pseudo_code(this) &
       RESULT   (pseudo_code)
    ! ==--------------------------------------------------------------==
    ! Quick hack to make the pseudo-potential check work with the new 
    ! driver
    !                               20.04.2017 M. P. Bircher @ LCBC/EPFL

    CLASS(cp_xc_env_t), INTENT(in)           :: this

    INTEGER                                  :: i_func
    INTEGER                                  :: old_mgcx, old_mgcc

    INTEGER                                  :: pseudo_code


    old_mgcx = 0.0_real_8 ; old_mgcc = 0.0_real_8

    DO i_func=1, this%n_funcs
       !
       ! Exchange
       !
       IF (TRIM(ADJUSTL(this%funcs( i_func ))) == 'GGA_X_B88') THEN
          IF (old_mgcx == 0.0_real_8) THEN
              old_mgcx = mgcx_is_becke88
          ELSE
              old_mgcx = -1.0_real_8
              EXIT
          END IF
       ELSE IF (TRIM(ADJUSTL(this%funcs( i_func ))) == 'GGA_X_OPTX') THEN
          IF (old_mgcx == 0.0_real_8) THEN
              old_mgcx = mgcx_is_optx
          ELSE
              old_mgcx = -1.0_real_8
              EXIT
          END IF
       ELSE IF (TRIM(ADJUSTL(this%funcs( i_func ))) == 'GGA_X_PBE') THEN
          IF (old_mgcx == 0.0_real_8) THEN
              old_mgcx = mgcx_is_pbex
          ELSE
              old_mgcx = -1.0_real_8
              EXIT
          END IF
       ELSE IF (TRIM(ADJUSTL(this%funcs( i_func ))) == 'GGA_X_PBE_R' .OR. &
                TRIM(ADJUSTL(this%funcs( i_func ))) == 'GGA_X_REVPBE') THEN !bugfix
          IF (old_mgcx == 0.0_real_8) THEN
              old_mgcx = mgcx_is_revpbex
          ELSE
              old_mgcx = -1.0_real_8
              EXIT
          END IF
       !
       ! Correlation
       !
       ELSE IF (TRIM(ADJUSTL(this%funcs( i_func ))) == 'GGA_C_LYP') THEN
          IF (old_mgcc == 0.0_real_8) THEN
              old_mgcc = mgcc_is_lyp
          ELSE
              old_mgcc = 0.0_real_8
              EXIT
          END IF
       ELSE IF (TRIM(ADJUSTL(this%funcs( i_func ))) == 'GGA_C_P86') THEN
          IF (old_mgcc == 0.0_real_8) THEN
              old_mgcc = mgcc_is_perdew86
          ELSE
              old_mgcc = 0.0_real_8
              EXIT
          END IF
       ELSE IF (TRIM(ADJUSTL(this%funcs( i_func ))) == 'GGA_C_PBE') THEN
          IF (old_mgcc == 0.0_real_8) THEN
              old_mgcc = mgcc_is_pbec
          ELSE
              old_mgcc = 0.0_real_8
              EXIT
          END IF
       !
       ! Combinations
       !
       ELSE IF (TRIM(ADJUSTL(this%funcs( i_func ))) == 'GGA_XC_BLYP') THEN
          IF (old_mgcx == 0.0_real_8) THEN
              old_mgcx = mgcx_is_becke88
          ELSE
              old_mgcx = 0.0_real_8
              EXIT
          END IF
          IF (old_mgcc == 0.0_real_8) THEN
              old_mgcc = mgcc_is_lyp
          ELSE
              old_mgcc = 0.0_real_8
              EXIT
          END IF
       ELSE IF (TRIM(ADJUSTL(this%funcs( i_func ))) == 'GGA_XC_BP86' .OR. &
                TRIM(ADJUSTL(this%funcs( i_func ))) == 'GGA_XC_BP') THEN
          IF (old_mgcx == 0.0_real_8) THEN
              old_mgcx = mgcx_is_becke88
          ELSE
              old_mgcx = 0.0_real_8
              EXIT
          END IF
          IF (old_mgcc == 0.0_real_8) THEN
              old_mgcc = mgcc_is_perdew86
          ELSE
              old_mgcc = 0.0_real_8
              EXIT
          END IF
       ELSE IF (TRIM(ADJUSTL(this%funcs( i_func ))) == 'GGA_XC_OLYP') THEN
          IF (old_mgcx == 0.0_real_8) THEN
              old_mgcx = mgcx_is_optx
          ELSE
              old_mgcx = 0.0_real_8
              EXIT
          END IF
          IF (old_mgcc == 0.0_real_8) THEN
              old_mgcc = mgcc_is_lyp
          ELSE
              old_mgcc = 0.0_real_8
              EXIT
          END IF
       ELSE IF (TRIM(ADJUSTL(this%funcs( i_func ))) == 'GGA_XC_OPBE') THEN
          IF (old_mgcx == 0.0_real_8) THEN
              old_mgcx = mgcx_is_optx
          ELSE
              old_mgcx = 0.0_real_8
              EXIT
          END IF
          IF (old_mgcc == 0.0_real_8) THEN
              old_mgcc = mgcc_is_pbec
          ELSE
              old_mgcc = 0.0_real_8
              EXIT
          END IF
       ELSE IF (TRIM(ADJUSTL(this%funcs( i_func ))) == 'GGA_XC_PBE') THEN
          IF (old_mgcx == 0.0_real_8) THEN
              old_mgcx = mgcx_is_pbex
          ELSE
              old_mgcx = 0.0_real_8
              EXIT
          END IF
          IF (old_mgcc == 0.0_real_8) THEN
              old_mgcc = mgcc_is_pbec
          ELSE
              old_mgcc = 0.0_real_8
              EXIT
          END IF
       ELSE IF (TRIM(ADJUSTL(this%funcs( i_func ))) == 'GGA_XC_PBE_R' .OR. &
                TRIM(ADJUSTL(this%funcs( i_func ))) == 'GGA_XC_REVPBE') THEN !bugfix
          IF (old_mgcx == 0.0_real_8) THEN
              old_mgcx = mgcx_is_revpbex
          ELSE
              old_mgcx = 0.0_real_8
              EXIT
          END IF
          IF (old_mgcc == 0.0_real_8) THEN
              old_mgcc = mgcc_is_pbec
          ELSE
              old_mgcc = 0.0_real_8
              EXIT
          END IF
       END IF
    END DO 

    pseudo_code = 10._real_8*old_mgcx + old_mgcc

    ! ==--------------------------------------------------------------==
  END FUNCTION cp_xc_env_get_pseudo_code
  ! ==================================================================

  ! ==================================================================
  ! Utilities for setting and nullifying procedure pointers to
  ! internal functionals and precomputation of scratch space
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE cp_xc_set_functional( cp_xc, functional, dual_is_set)
    ! ==--------------------------------------------------------------==
    ! Sets a functional pointer to the procedure specified in
    ! cp_xc%cpfunc for n_funcs functionals; sets Coulomb attenuation
    !                               07.10.2016 M. P. Bircher @ LCBC/EPFL


    CLASS(cp_xc_t), INTENT(in)               :: cp_xc
    TYPE(cp_xc_functional_p_t), INTENT(out)  :: functional
    LOGICAL, INTENT(in), OPTIONAL            :: dual_is_set

    CHARACTER(len=*), PARAMETER              :: procedureN = 'set_functional'

    INTEGER                                  :: i_func

    LOGICAL                                  :: initialisation_OK
    LOGICAL                                  :: dual_set_manually
    LOGICAL                                  :: dual_OK

    dual_set_manually = .false.
    IF (present(dual_is_set)) dual_set_manually = dual_is_set

    DO i_func=1,cp_xc%cpfunc%n_funcs

       initialisation_OK = .false.
       dual_OK           = .true.

       SELECT CASE(TRIM(ADJUSTL(cp_xc%cpfunc%funcs( i_func ))))
       !
       ! Exchange functionals
       !
       CASE("CP_LDA_X")
          functional%K( i_func )%compute  => CP_LDA_X
          functional%Ks( i_func )%compute => CP_SPIN_LDA_X
          initialisation_OK               =  cp_lda_x_check(cp_xc%cpfunc%funcs( i_func ))
       CASE("CP_GGA_X_B88")
          functional%K( i_func )%compute  => CP_GGA_X_B88
          functional%Ks( i_func )%compute => CP_SPIN_GGA_X_B88
          initialisation_OK               =  cp_gga_x_check(cp_xc%cpfunc%funcs( i_func ))
       CASE("CP_GGA_X_HCTH")
          functional%K( i_func )%compute  => CP_GGA_X_HCTH
          functional%Ks( i_func )%compute => CP_SPIN_GGA_X_HCTH
          initialisation_OK               =  cp_gga_x_check(cp_xc%cpfunc%funcs( i_func ))
       CASE("CP_GGA_X_OPTX")
          functional%K( i_func )%compute  => CP_GGA_X_OPTX
          functional%Ks( i_func )%compute => CP_SPIN_GGA_X_OPTX
          initialisation_OK               =  cp_gga_x_check(cp_xc%cpfunc%funcs( i_func ))
       CASE("CP_GGA_X_PBE")
          functional%K( i_func )%compute  => CP_GGA_X_PBE
          functional%Ks( i_func )%compute => CP_SPIN_GGA_X_PBE
          initialisation_OK               =  cp_gga_x_check(cp_xc%cpfunc%funcs( i_func ))
       CASE("CP_GGA_X_PBE_SOL")
          functional%K( i_func )%compute  => CP_GGA_X_PBE_SOL
          functional%Ks( i_func )%compute => CP_SPIN_GGA_X_PBE_SOL
          initialisation_OK               =  .true.
          initialisation_OK               =  cp_gga_x_check(cp_xc%cpfunc%funcs( i_func ))
       CASE("CP_GGA_X_REVPBE")
          functional%K( i_func )%compute  => CP_GGA_X_REVPBE
          functional%Ks( i_func )%compute => CP_SPIN_GGA_X_REVPBE
          initialisation_OK               =  .true.
          initialisation_OK               =  cp_gga_x_check(cp_xc%cpfunc%funcs( i_func ))
       CASE("CP_GGA_X_PBE_FLEX")
          functional%K( i_func )%compute  => CP_GGA_X_PBE_FLEX
          functional%Ks( i_func )%compute => CP_SPIN_GGA_X_PBE_FLEX
          initialisation_OK               =  cp_gga_x_param%init
          initialisation_OK               =  cp_gga_x_check(cp_xc%cpfunc%funcs( i_func ))
       CASE("CP_MGGA_X_TPSS")
          functional%K( i_func )%compute  => CP_MGGA_X_TPSS
          functional%Ks( i_func )%compute => CP_SPIN_MGGA_X_TPSS
          initialisation_OK               =  cp_mgga_x_check(cp_xc%cpfunc%funcs( i_func ))
       !  CASE("CP_MGGA_X_VS98")
       !     functional%K( i_func )%compute  => CP_MGGA_X_VS98
       !     functional%Ks( i_func )%compute => CP_SPIN_MGGA_X_VS98
       CASE("CP_MGGA_X_M05_M06")
          functional%K( i_func )%compute  => CP_MGGA_X_M05_M06
          functional%Ks( i_func )%compute => CP_SPIN_MGGA_X_M05_M06
          initialisation_OK               =  cp_mgga_x_check(cp_xc%cpfunc%funcs( i_func ))
          dual_OK                         =  dual_set_manually
       CASE("CP_MGGA_X_M08_M11")
          functional%K( i_func )%compute  => CP_MGGA_X_M08_M11
          functional%Ks( i_func )%compute => CP_SPIN_MGGA_X_M08_M11
          initialisation_OK               =  cp_mgga_x_check(cp_xc%cpfunc%funcs( i_func ))
          dual_OK                         =  dual_set_manually
       CASE("CP_MGGA_X_MN12_MN15")
          functional%K( i_func )%compute  => CP_MGGA_X_MN12_MN15
          functional%Ks( i_func )%compute => CP_SPIN_MGGA_X_MN12_MN15
          initialisation_OK               =  cp_mgga_x_check(cp_xc%cpfunc%funcs( i_func ))
          dual_OK                         =  dual_set_manually
       !
       ! Correlation functionals
       !
       CASE("CP_LDA_C_VWN")
          functional%K( i_func )%compute  => CP_LDA_C_VWN
          functional%Ks( i_func )%compute => CP_SPIN_LDA_C_VWN
          initialisation_OK               =  cp_lda_c_check(cp_xc%cpfunc%funcs( i_func ))
       CASE("CP_LDA_C_PZ")
          functional%K( i_func )%compute  => CP_LDA_C_PZ
          functional%Ks( i_func )%compute => CP_SPIN_LDA_C_PZ
          initialisation_OK               =  cp_lda_c_check(cp_xc%cpfunc%funcs( i_func ))
       CASE("CP_LDA_C_PW")
          functional%K( i_func )%compute  => CP_LDA_C_PW
          functional%Ks( i_func )%compute => CP_SPIN_LDA_C_PW
          initialisation_OK               =  cp_lda_c_check(cp_xc%cpfunc%funcs( i_func ))
       CASE("CP_LDA_C_OB_PW")
          functional%K( i_func )%compute  => CP_LDA_C_OB_PW
          functional%Ks( i_func )%compute => CP_SPIN_LDA_C_OB_PW
          initialisation_OK               =  cp_lda_c_check(cp_xc%cpfunc%funcs( i_func ))
       CASE("CP_GGA_C_HCTH")
          functional%K( i_func )%compute  => CP_GGA_C_HCTH
          functional%Ks( i_func )%compute => CP_SPIN_GGA_C_HCTH
          initialisation_OK               =  cp_gga_c_check(cp_xc%cpfunc%funcs( i_func ))
       CASE("CP_GGA_C_LYP")
          functional%K( i_func )%compute  => CP_GGA_C_LYP
          functional%Ks( i_func )%compute => CP_SPIN_GGA_C_LYP
          initialisation_OK               =  cp_gga_c_check(cp_xc%cpfunc%funcs( i_func ))
       CASE("CP_GGA_C_P86")
          functional%K( i_func )%compute  => CP_GGA_C_P86
          functional%Ks( i_func )%compute => CP_SPIN_GGA_C_P86
          initialisation_OK               =  cp_gga_c_check(cp_xc%cpfunc%funcs( i_func ))
       CASE("CP_GGA_C_PBE")
          functional%K( i_func )%compute  => CP_GGA_C_PBE
          functional%Ks( i_func )%compute => CP_SPIN_GGA_C_PBE
          initialisation_OK               =  cp_gga_c_check(cp_xc%cpfunc%funcs( i_func ))
       CASE("CP_GGA_C_PBE_SOL")
          functional%K( i_func )%compute  => CP_GGA_C_PBE_SOL
          functional%Ks( i_func )%compute => CP_SPIN_GGA_C_PBE_SOL
          initialisation_OK               =  cp_gga_c_check(cp_xc%cpfunc%funcs( i_func ))
       CASE("CP_GGA_C_PBE_FLEX")
          functional%K( i_func )%compute  => CP_GGA_C_PBE_FLEX
          functional%Ks( i_func )%compute => CP_SPIN_GGA_C_PBE_FLEX
          initialisation_OK               =  cp_gga_c_check(cp_xc%cpfunc%funcs( i_func ))
       CASE("CP_MGGA_C_TPSS")
          functional%K( i_func )%compute  => CP_MGGA_C_TPSS
          functional%Ks( i_func )%compute => CP_SPIN_MGGA_C_TPSS
          initialisation_OK               =  cp_mgga_c_check(cp_xc%cpfunc%funcs( i_func ))
       !  CASE("CP_MGGA_C_VS98")
       !     functional%K( i_func )%compute  => CP_MGGA_C_VS98
       !     functional%Ks( i_func )%compute => CP_SPIN_MGGA_C_VS98
       !     initialisation_OK               =  cp_mgga_c_check(cp_xc%cpfunc%funcs( i_func ))
       CASE("CP_MGGA_C_M05_M06")
          functional%K( i_func )%compute  => CP_MGGA_C_M05_M06
          functional%Ks( i_func )%compute => CP_SPIN_MGGA_C_M05_M06
          initialisation_OK               =  cp_mgga_c_check(cp_xc%cpfunc%funcs( i_func ))
          dual_OK                         =  dual_set_manually
       CASE("CP_MGGA_C_M08_M11")
          functional%K( i_func )%compute  => CP_MGGA_C_M08_M11
          functional%Ks( i_func )%compute => CP_SPIN_MGGA_C_M08_M11
          initialisation_OK               =  cp_mgga_c_check(cp_xc%cpfunc%funcs( i_func ))
          dual_OK                         =  dual_set_manually
          !
          ! This should not happen
          !
       CASE default
          CALL stopgm(procedureN,'Could not set functional pointer; &
               &functional '//trim(adjustl(cp_xc%cpfunc%funcs( i_func )))//' n/a',&
               __LINE__,__FILE__)
       END SELECT
       
       IF (.not. dual_OK)  CALL stopgm(procedureN,'Invalid setup for Minnesota functional. Please set DUAL to 8-12 '// &
                                 'and refer to the manual for more information. See also: '// &
                                 'dx.doi.org/10.1021/acs.jctc.8b00897', &
                                 __LINE__,__FILE__)

       IF (.not. initialisation_OK) CALL stopgm(procedureN,'Functional parameters have not been properly set; &
               &concerns '//trim(adjustl(cp_xc%cpfunc%funcs( i_func ))), &
               __LINE__,__FILE__)

    END DO
    !
    ! Coulomb attenuation (only option for the time being)
    !      
    SELECT CASE(cp_xc%cpfunc%is_coulomb_attenuated)
    CASE(.TRUE.)
       functional%attenuate => CP_GGA_SCREENING_CAM
    CASE(.FALSE.)
       functional%attenuate => NULL()
    END SELECT

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_xc_set_functional
  ! ==================================================================
  ELEMENTAL SUBROUTINE cp_xc_purge_functional( functional ) 
    ! ==--------------------------------------------------------------==

    TYPE(cp_xc_functional_p_t), INTENT(out)  :: functional

    INTEGER                                  :: i_func

    DO i_func=1,SIZE(functional%K)
       functional%K( i_func )%compute  => NULL()
       functional%Ks( i_func )%compute => NULL()
    END DO
    functional%attenuate => NULL()

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_xc_purge_functional
  ! ==================================================================
  SUBROUTINE cp_xc_set_functional_derivatives( cp_xc, functional )
    ! ==--------------------------------------------------------------==
    ! Sets a functional pointer to the procedure specified in
    ! cp_xc%cpfunc for n_funcs functionals; sets Coulomb attenuation
    !                               07.10.2016 M. P. Bircher @ LCBC/EPFL


    CLASS(cp_xc_t), INTENT(in)               :: cp_xc
    TYPE(cp_dxc_functional_p_t), INTENT(out) :: functional

    CHARACTER(len=*), PARAMETER              :: procedureN = 'set_functional_derivatives'

    INTEGER                                  :: i_func

    LOGICAL                                  :: initialisation_OK

    DO i_func=1,cp_xc%cpfunc%n_funcs

       initialisation_OK = .false.

       SELECT CASE(TRIM(ADJUSTL(cp_xc%cpfunc%funcs( i_func ))))
          !
          ! Exchange
          !
       CASE("CP_GGA_X_B88")
          functional%dK( i_func )%compute  => CP_dGGA_X_B88
          functional%dKs( i_func )%compute => CP_SPIN_dGGA_X_B88
          initialisation_OK                =  cp_gga_x_check(cp_xc%cpfunc%funcs( i_func ))
          !
       CASE("CP_GGA_X_PBE")
          functional%dK( i_func )%compute  => CP_dGGA_X_PBE
          functional%dKs( i_func )%compute => CP_SPIN_dGGA_X_PBE
          initialisation_OK                =  cp_gga_x_check(cp_xc%cpfunc%funcs( i_func ))
          !
          ! Correlation
          !
       CASE("CP_GGA_C_LYP")
          functional%dK( i_func )%compute  => CP_dGGA_C_LYP
          functional%dKs( i_func )%compute => CP_SPIN_dGGA_C_LYP
          initialisation_OK                =  cp_gga_c_check(cp_xc%cpfunc%funcs( i_func ))
          !
       CASE("CP_GGA_C_P86")
          functional%dK( i_func )%compute  => CP_dGGA_C_P86
          functional%dKs( i_func )%compute => CP_SPIN_dGGA_C_P86
          initialisation_OK                =  cp_gga_c_check(cp_xc%cpfunc%funcs( i_func ))
          !
       CASE("CP_GGA_C_PBE")
          functional%dK( i_func )%compute  => CP_dGGA_C_PBE
          functional%dKs( i_func )%compute => CP_SPIN_dGGA_C_PBE
          initialisation_OK                =  cp_gga_c_check(cp_xc%cpfunc%funcs( i_func ))
          !
          ! This should not happen
          !
       CASE default
          CALL stopgm(procedureN,'Could not set functional pointer; &
               &analytical derivatives for functional '//trim(adjustl(cp_xc%cpfunc%funcs( i_func )))//' n/a',&
               __LINE__,__FILE__)
       END SELECT

       IF (.not. initialisation_OK) CALL stopgm(procedureN,'Functional parameters have not been properly set; &
               &concerns '//trim(adjustl(cp_xc%cpfunc%funcs( i_func ))), &
               __LINE__,__FILE__)

    END DO
    !
    ! Coulomb attenuation is NYI
    !      
    IF (cp_xc%cpfunc%is_coulomb_attenuated) CALL stopgm(procedureN,'Analytical derivatives for linear response are &
               &unavailable for Coulomb-attenuated functionals',&
                __LINE__,__FILE__)

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_xc_set_functional_derivatives
  ! ==================================================================
  ELEMENTAL SUBROUTINE cp_xc_purge_functional_derivatives( functional ) 
    ! ==--------------------------------------------------------------==


    TYPE(cp_dxc_functional_p_t), INTENT(out) :: functional

    INTEGER                                  :: i_func

    DO i_func=1,SIZE(functional%dK)
       functional%dK( i_func )%compute  => NULL()
       functional%dKs( i_func )%compute => NULL()
    END DO

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_xc_purge_functional_derivatives
  ! ==================================================================
  
  ! ==================================================================
  ! Generic procedure pointer targets for scratch generation
  ! ==================================================================

  ! ==================================================================
  ELEMENTAL SUBROUTINE cp_xc_set_scratch(get_scratch, has_gradient, has_tau)
    ! ==--------------------------------------------------------------==
    ! Sets appropriate scratch procedure to ensure unneeded variables are
    ! zero


    TYPE(cp_xc_scratch_p_t), INTENT(out)     :: get_scratch
    LOGICAL, INTENT(in)                      :: has_gradient, has_tau

    IF (has_tau) THEN
       get_scratch%restricted   => get_meta_scratch
       get_scratch%unrestricted => get_meta_spin_scratch
    ELSE IF (has_gradient .AND. .NOT. has_tau) THEN
       get_scratch%restricted   => get_GGA_scratch
       get_scratch%unrestricted => get_GGA_spin_scratch
    ELSE
       get_scratch%restricted   => get_LDA_scratch
       get_scratch%unrestricted => get_LDA_spin_scratch
    END IF

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_xc_set_scratch
  ! ==================================================================
  ELEMENTAL SUBROUTINE cp_xc_purge_scratch(get_scratch)
    ! ==--------------------------------------------------------------==


    TYPE(cp_xc_scratch_p_t), INTENT(inout)   :: get_scratch

    get_scratch%restricted   => NULL()
    get_scratch%unrestricted => NULL()

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_xc_purge_scratch
  ! ==================================================================
  ELEMENTAL SUBROUTINE cp_xc_set_scratch_derivatives(get_scratch, has_gradient, has_tau)
    ! ==--------------------------------------------------------------==
    ! Sets appropriate scratch procedure to ensure unneeded variables are
    ! zero


    TYPE(cp_dxc_scratch_p_t), INTENT(out)    :: get_scratch
    LOGICAL, INTENT(in)                      :: has_gradient, has_tau

    get_scratch%restricted   => get_dGGA_scratch
    get_scratch%unrestricted => get_dGGA_spin_scratch

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_xc_set_scratch_derivatives
  ! ==================================================================
  ELEMENTAL SUBROUTINE cp_xc_purge_scratch_derivatives(get_scratch)
    ! ==--------------------------------------------------------------==


    TYPE(cp_dxc_scratch_p_t), INTENT(inout)  :: get_scratch

    get_scratch%restricted   => NULL()
    get_scratch%unrestricted => NULL()

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_xc_purge_scratch_derivatives
  ! ==================================================================
  PURE SUBROUTINE get_lda_scratch(scratch,n_in,abs_grad_in,tau_in)
    ! ==--------------------------------------------------------------==
    ! Converts LDA to LSDA, empty gradients
    ! ==--------------------------------------------------------------==


    TYPE(cp_xc_scratch_t), INTENT(out)       :: scratch
    REAL(real_8), INTENT(in)                 :: n_in, abs_grad_in, tau_in

    REAL(real_8)                             :: n_2_3


    scratch%n     = 0.5_real_8*n_in
    scratch%n_1_3 = scratch%n**(1.0_real_8/3.0_real_8)
    n_2_3         = scratch%n_1_3 * scratch%n_1_3
    scratch%n_4_3 = n_2_3 * n_2_3
    !
    scratch%abs_g = 0.0_real_8 
    scratch%g_2   = 0.0_real_8 
    !
    scratch%tau   = 0.0_real_8

    ! ==--------------------------------------------------------------==
  END SUBROUTINE get_lda_scratch
  ! ==================================================================
  PURE SUBROUTINE get_lda_spin_scratch(scratch,n_in,g_2_in,tau_in)
    ! ==--------------------------------------------------------------==
    ! LSDA utilities, no gradients.
    ! ==--------------------------------------------------------------==


    TYPE(cp_xc_scratch_t), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(out)                            :: scratch
    REAL(real_8), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(in)                             :: n_in, g_2_in
    REAL(real_8), &
      DIMENSION(cp_xc_spin_pairs), &
      INTENT(in)                             :: tau_in

    REAL(real_8), &
      DIMENSION(cp_xc_spin_components)       :: n_2_3


    scratch(:)%n     = n_in(:)
    scratch(:)%n_1_3 = scratch(:)%n**(1.0_real_8/3.0_real_8)
    n_2_3(:)         = scratch(:)%n_1_3 * scratch(:)%n_1_3
    scratch(:)%n_4_3 = n_2_3(:) * n_2_3(:)
    !
    scratch(:)%abs_g = 0.0_real_8 
    scratch(:)%g_2   = 0.0_real_8 
    !
    scratch(:)%tau   = 0.0_real_8

    ! ==--------------------------------------------------------------==
  END SUBROUTINE get_lda_spin_scratch
  ! ==================================================================
  PURE SUBROUTINE get_gga_scratch(scratch,n_in,g_2_in,tau_in)
    ! ==--------------------------------------------------------------==
    ! Converts LDA to LSDA densities. Including gradients.
    ! ==--------------------------------------------------------------==


    TYPE(cp_xc_scratch_t), INTENT(out)       :: scratch
    REAL(real_8), INTENT(in)                 :: n_in, g_2_in, tau_in

    REAL(real_8)                             :: n_2_3


    scratch%n     = 0.5_real_8*n_in
    scratch%n_1_3 = scratch%n**(1.0_real_8/3.0_real_8)
    n_2_3         = scratch%n_1_3 * scratch%n_1_3
    scratch%n_4_3 = n_2_3 * n_2_3
    !
    scratch%g_2   = 0.25_real_8*g_2_in
    scratch%abs_g = SQRT(scratch%g_2)
    !
    scratch%tau   = 0.0_real_8

    ! ==--------------------------------------------------------------==
  END SUBROUTINE get_gga_scratch
  ! ==================================================================
  PURE SUBROUTINE get_gga_spin_scratch(scratch,n_in,g_2_in,tau_in)
    ! ==--------------------------------------------------------------==
    ! LSDA utilities. Including gradients.
    ! ==--------------------------------------------------------------==


    TYPE(cp_xc_scratch_t), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(out)                            :: scratch
    REAL(real_8), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(in)                             :: n_in, g_2_in
    REAL(real_8), &
      DIMENSION(cp_xc_spin_pairs), &
      INTENT(in)                             :: tau_in

    REAL(real_8), &
      DIMENSION(cp_xc_spin_components)       :: n_2_3


    scratch(:)%n       = n_in(:)
    scratch(:)%n_1_3   = scratch(:)%n**(1.0_real_8/3.0_real_8)
    n_2_3(:)           = scratch(:)%n_1_3 * scratch(:)%n_1_3
    scratch(:)%n_4_3   = n_2_3(:) * n_2_3(:)
    !
    scratch(:)%g_2     = g_2_in(:)
    scratch(a:b)%abs_g = SQRT(g_2_in(a:b))
    !
    ! This term is never needed and must not be used
    ! (proper sqrt(g_2(ab)) may be NaN, so we set it to 0 instead)
    !
    scratch(ab)%abs_g  = 0.0_real_8
    !
    scratch(:)%tau     = 0.0_real_8

    ! ==--------------------------------------------------------------==
  END SUBROUTINE get_gga_spin_scratch
  ! ==================================================================
  PURE SUBROUTINE get_meta_scratch(scratch,n_in,g_2_in,tau_in)
    ! ==--------------------------------------------------------------==
    ! Converts LDA to LSDA densities. Including gradients and tau.
    ! ==--------------------------------------------------------------==


    TYPE(cp_xc_scratch_t), INTENT(out)       :: scratch
    REAL(real_8), INTENT(in)                 :: n_in, g_2_in, tau_in

    REAL(real_8)                             :: n_2_3


    scratch%n     = 0.5_real_8*n_in
    scratch%n_1_3 = scratch%n**(1.0_real_8/3.0_real_8)
    n_2_3         = scratch%n_1_3 * scratch%n_1_3
    scratch%n_4_3 = n_2_3 * n_2_3
    !
    scratch%g_2   = 0.25_real_8*g_2_in
    scratch%abs_g = SQRT(scratch%g_2)
    !
    scratch%tau   = 0.5_real_8*tau_in

    ! ==--------------------------------------------------------------==
  END SUBROUTINE get_meta_scratch
  ! ==================================================================
  PURE SUBROUTINE get_meta_spin_scratch(scratch,n_in,g_2_in,tau_in)
    ! ==--------------------------------------------------------------==
    ! Meta-GGA utilities. Including gradients and tau.
    ! ==--------------------------------------------------------------==


    TYPE(cp_xc_scratch_t), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(out)                            :: scratch
    REAL(real_8), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(in)                             :: n_in, g_2_in
    REAL(real_8), &
      DIMENSION(cp_xc_spin_pairs), &
      INTENT(in)                             :: tau_in

    REAL(real_8), &
      DIMENSION(cp_xc_spin_components)       :: n_2_3


    scratch(:)%n       = n_in(:)
    scratch(:)%n_1_3   = scratch(:)%n**(1.0_real_8/3.0_real_8)
    n_2_3(:)           = scratch(:)%n_1_3 * scratch(:)%n_1_3
    scratch(:)%n_4_3   = n_2_3(:) * n_2_3(:)
    !
    scratch(:)%g_2     = g_2_in(:)
    scratch(a:b)%abs_g = SQRT(g_2_in(a:b))
    !
    scratch(a:b)%tau   = tau_in(a:b)
    !
    ! This term is never needed and must not be used
    ! (proper sqrt(g_2(ab)) may be NaN, so we set it to 0 instead)
    !
    scratch(ab)%abs_g  = 0.0_real_8
    scratch(ab)%tau    = 0.0_real_8

    ! ==--------------------------------------------------------------==
  END SUBROUTINE get_meta_spin_scratch
  ! ==================================================================
  PURE SUBROUTINE get_dgga_scratch(scratch,n_in,dn_in,g_2_in,dgrad_in)
    ! ==--------------------------------------------------------------==
    ! Converts LDA to LSDA densities. Including gradients.
    ! ==--------------------------------------------------------------==


    TYPE(cp_dxc_scratch_t), INTENT(out)      :: scratch
    REAL(real_8), INTENT(in)                 :: n_in, dn_in
    REAL(real_8), &
      DIMENSION(cp_xc_xyz), INTENT(in)       :: dgrad_in
    REAL(real_8), &
      DIMENSION(cp_xc_abs_xyz), INTENT(in)   :: g_2_in

    scratch%n       = n_in        ! These have all been multiplied by
    scratch%dn      = dn_in       ! 0.5 long before being passed
    !
    scratch%g_2     = g_2_in(1)
    scratch%gx_2    = g_2_in(2)
    scratch%gy_2    = g_2_in(3)
    scratch%gz_2    = g_2_in(4)
    !
    scratch%dgx     = dgrad_in(1)
    scratch%dgy     = dgrad_in(2)
    scratch%dgz     = dgrad_in(3)

    ! ==--------------------------------------------------------------==
  END SUBROUTINE get_dgga_scratch
  ! ==================================================================
  PURE SUBROUTINE get_dgga_spin_scratch(scratch,n_in,dn_in,g_2_in,dgrad_in)
    ! ==--------------------------------------------------------------==
    ! LSDA utilities. Including gradients.
    ! ==--------------------------------------------------------------==


    TYPE(cp_dxc_scratch_t), &
      DIMENSION(cp_xc_spin_pairs), &
      INTENT(out)                            :: scratch
    REAL(real_8), &
      DIMENSION(cp_xc_spin_pairs), &
      INTENT(in)                             :: n_in, dn_in
    REAL(real_8), &
      DIMENSION(cp_xc_spin_pairs*cp_xc_xyz), & 
      INTENT(in)                             :: dgrad_in
    REAL(real_8), &
      DIMENSION(cp_xc_spin_pairs*cp_xc_abs_xyz), &
      INTENT(in)                             :: g_2_in

    scratch(:)%n    = n_in(:)
    scratch(:)%dn   = dn_in(:)
    !
    scratch(a)%g_2  = g_2_in(1)
    scratch(a)%gx_2 = g_2_in(2)
    scratch(a)%gy_2 = g_2_in(3)
    scratch(a)%gz_2 = g_2_in(4)
    scratch(b)%g_2  = g_2_in(5)
    scratch(b)%gx_2 = g_2_in(6)
    scratch(b)%gy_2 = g_2_in(7)
    scratch(b)%gz_2 = g_2_in(8)
    !
    scratch(a)%dgx  = dgrad_in(1)
    scratch(a)%dgy  = dgrad_in(2)
    scratch(a)%dgz  = dgrad_in(3)
    scratch(b)%dgx  = dgrad_in(4)
    scratch(b)%dgy  = dgrad_in(5)
    scratch(b)%dgz  = dgrad_in(6)

    ! ==--------------------------------------------------------------==
  END SUBROUTINE get_dgga_spin_scratch
  ! ==================================================================
END MODULE cp_xc_utils
