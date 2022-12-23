! ==================================================================
! Provides: - GGA exchange, rewritten or adapted from OLDCODE
!
! Input:  Spin-polarised density (1/2 n for closed shell systems)
! Output: For spin-component only
! Convention: E_x = -0.5 \sum_s rho_s**(4/3)*K_sx
! Most functionals: E_x = 0.5 \sum_s (2rho_s)**(4/3) K_x
! Therefore   K_sx = 2**(4/3)K_x
!
!                             02.10.2016 - M. P. Bircher @ LCBC/EPFL
! ==================================================================
MODULE cp_gga_exchange_utils
  USE cnst,                            ONLY: pi
  USE cpfunc_types,                    ONLY: cp_xc_functional_t,&
                                             cp_xc_scratch_t,&
                                             cp_xc_spin_components
  USE func,                            ONLY: func2
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER :: a                          = 1
  INTEGER, PARAMETER :: b                          = 2
  INTEGER, PARAMETER :: ab                         = 3

  INTEGER, PARAMETER :: hcth_gm                    = 4

  REAL(real_8), PARAMETER :: slater_prefactor      = 1.10783814957303361_real_8
  REAL(real_8), PARAMETER :: sqrt_of_pi            = SQRT(pi)
  REAL(real_8), PARAMETER :: four_thirds           = 4._real_8/3._real_8
  REAL(real_8), PARAMETER :: one_third             = 1.0_real_8/3.0_real_8
  REAL(real_8), PARAMETER :: two_to_one_third      = 2._real_8**(1._real_8/3._real_8)
  REAL(real_8), PARAMETER :: two_to_four_thirds    = 2._real_8**(4._real_8/3._real_8)
  REAL(real_8), PARAMETER :: oo_two_to_four_thirds = 1._real_8/two_to_four_thirds

  !
  ! Customisable functionals (FLEX type)
  !
  TYPE, PRIVATE :: hcth_x_t
    REAL(real_8)                      :: gamma = 0.004_real_8
    REAL(real_8)                      :: c_0   =    0.109320e+01_real_8
    REAL(real_8), DIMENSION(hcth_gm)  :: c     = (/-0.744056e+00_real_8, &
                                                               0.559920e+01_real_8, &
                                                              -0.678549e+01_real_8, &
                                                               0.449357e+01_real_8 /)
  END TYPE hcth_x_t
  ! 
  TYPE, PRIVATE :: pbe_x_t
     REAL(real_8)                     :: kappa = 0.8040_real_8
     REAL(real_8)                     :: mu    = 0.2195149727645171_real_8
  END TYPE pbe_x_t
  !
  TYPE, PRIVATE :: cp_gga_x_param_t
     LOGICAL                          :: init = .false.
     TYPE(hcth_x_t)                   :: hcth
     TYPE(pbe_x_t)                    :: pbe
  END TYPE cp_gga_x_param_t
  !
  TYPE(cp_gga_x_param_t), SAVE, PUBLIC :: cp_gga_x_param

  PUBLIC :: CP_GGA_X_B88
  PUBLIC :: CP_GGA_X_HCTH
  PUBLIC :: CP_GGA_X_PBE
  PUBLIC :: CP_GGA_X_REVPBE
  PUBLIC :: CP_GGA_X_PBE_SOL
  PUBLIC :: CP_GGA_X_PBE_FLEX
  PUBLIC :: CP_GGA_X_OPTX
  PUBLIC :: CP_SPIN_GGA_X_B88
  PUBLIC :: CP_SPIN_GGA_X_HCTH
  PUBLIC :: CP_SPIN_GGA_X_PBE
  PUBLIC :: CP_SPIN_GGA_X_REVPBE
  PUBLIC :: CP_SPIN_GGA_X_PBE_SOL
  PUBLIC :: CP_SPIN_GGA_X_PBE_FLEX
  PUBLIC :: CP_SPIN_GGA_X_OPTX
  PUBLIC :: CP_GGA_SCREENING_CAM

  PUBLIC :: CP_GGA_X_CHECK

CONTAINS

  ! ==================================================================
  PURE SUBROUTINE cp_spin_gga_x_b88( scratch,functional )
    ! ==--------------------------------------------------------------==


    TYPE(cp_xc_scratch_t), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(in)                             :: scratch
    TYPE(cp_xc_functional_t), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(inout)                          :: functional

    IF (scratch(a)%n >= 1.0e-10_real_8) &
         CALL cp_gga_x_b88( scratch(a), functional(a) )
    IF (scratch(b)%n >= 1.0e-10_real_8) &
         CALL cp_gga_x_b88( scratch(b), functional(b) )

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_spin_gga_x_b88
  ! ==================================================================
  PURE SUBROUTINE cp_gga_x_b88( scratch,functional )
    ! ==--------------------------------------------------------------==
    ! Rewritten from scratch (not ported from oldcode)
    ! Output term by term identical to libxc
    !                                     09.09.2016 M.P. Bircher @ EPFL
    ! ==--------------------------------------------------------------==


    TYPE(cp_xc_scratch_t), INTENT(in)        :: scratch
    TYPE(cp_xc_functional_t), INTENT(inout)  :: functional

    REAL(real_8) :: b88_K, ddenominator_dx, denominator, dK_dx, dsinhfrac_dx, &
      dsinhln_dx, dx_dg, dx_dn, lda_K, sinh_frac, sinh_ln, x, x_2

!
! UEG
!

    lda_K             = two_to_four_thirds*slater_prefactor*func2%salpha
    !
    ! Gradient correction
    !
    x                 = scratch%abs_g / scratch%n_4_3
    x_2               = x*x
    !
    sinh_frac         = SQRT(1._real_8+x_2)
    sinh_ln           = LOG(x+sinh_frac)
    denominator       = 1._real_8 + 6._real_8*func2%bbeta*x*sinh_ln
    !
    dx_dn             = -four_thirds*scratch%abs_g / (scratch%n_4_3*scratch%n)
    dx_dg             = 1.0_real_8/scratch%n_4_3
    !
    dsinhfrac_dx      = x / sinh_frac
    dsinhln_dx        = (1.0_real_8+dsinhfrac_dx) / (x+sinh_frac)
    ddenominator_dx   = 6.0_real_8*func2%bbeta*sinh_ln &
         + 6.0_real_8*func2%bbeta*x*dsinhln_dx
    !
    b88_K             = 2.0_real_8*func2%bbeta*x_2 / denominator
    !
    ! Derivatives
    !
    ! Formally:
    !
    ! denominator_squared  = denominator*denominator
    ! functional%dK_dx = 2.0_real_8*func2%bbeta*(2.0_real_8*x/denominator &
    !                         - x_2*ddenominator_dx/denominator_squared)
    !
    ! ... we can use less terms for this:
    !
    dK_dx             = b88_K*(2.0_real_8/x - ddenominator_dx/denominator)
    !
    functional%dK_dn  = dK_dx*dx_dn
    functional%dK_dg  = dK_dx*dx_dg
    !
    functional%K      = ( b88_K + lda_K )
    !
    ! Proper conventions for readability
    !
    functional%sx_sigma = -0.5_real_8*( scratch%n_4_3*functional%K )
    functional%dsx_dn   = -0.5_real_8*( four_thirds*scratch%n_1_3*functional%K + scratch%n_4_3*functional%dK_dn )
    functional%dsx_dg   = -0.5_real_8*( scratch%n_4_3*functional%dK_dg )
    !

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_gga_x_b88
  ! ==================================================================
  PURE SUBROUTINE cp_spin_gga_x_pbe( scratch,functional )
    ! ==--------------------------------------------------------------==


    TYPE(cp_xc_scratch_t), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(in)                             :: scratch
    TYPE(cp_xc_functional_t), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(inout)                          :: functional

    IF (scratch(a)%n >= 1.0e-10_real_8) &
         CALL cp_gga_x_pbe( scratch(a), functional(a) )
    IF (scratch(b)%n >= 1.0e-10_real_8) &
         CALL cp_gga_x_pbe( scratch(b), functional(b) )

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_spin_gga_x_pbe
  ! ==================================================================
  PURE SUBROUTINE cp_gga_x_pbe( scratch,functional )
    ! ==--------------------------------------------------------------==
    ! Rewritten from scratch (not ported from oldcode)
    ! Output term by term identical to libxc
    !                                     09.09.2016 M.P. Bircher @ EPFL
    ! ==--------------------------------------------------------------==


    TYPE(cp_xc_scratch_t), INTENT(in)        :: scratch
    TYPE(cp_xc_functional_t), INTENT(inout)  :: functional

    REAL(real_8), PARAMETER :: kappa = 0.8040_real_8, &
      mu = 0.2195149727645171_real_8, mu_by_kappa = mu/kappa , &
      us = 0.161620459673995492_real_8

    REAL(real_8)                             :: ddenominator_ds, denominator, &
                                                denominator_2, dK_ds, ds_dg, &
                                                ds_dn, Fx, lda_K, s, s_2

!
! UEG
!

    lda_K            = two_to_four_thirds*slater_prefactor*func2%salpha
    !
    ! Gradient correction
    !
    s                = us/two_to_one_third*scratch%abs_g / scratch%n_4_3
    s_2              = s*s
    !
    ! Work with the spin-polarized enhancement factor (cf. PBE paper, p. 3867)
    ! Prefactor incorporated into UEG
    !
    denominator      = 1.0_real_8 + mu_by_kappa*s_2
    denominator_2    = denominator*denominator
    ddenominator_ds  = 2.0_real_8*mu_by_kappa*s
    !
    Fx               = ( 1.0_real_8 + kappa - kappa/denominator )
    functional%K     = lda_K*Fx
    !
    dK_ds            = lda_K*( ddenominator_ds*kappa/denominator_2 )
    !
    ds_dn            = -four_thirds*us/two_to_one_third*scratch%abs_g / (scratch%n_4_3*scratch%n)
    ds_dg            = us /( two_to_one_third * scratch%n_4_3 )
    !
    functional%dK_dn = dK_ds*ds_dn
    functional%dK_dg = dK_ds*ds_dg
    !
    functional%sx_sigma = -0.5_real_8*( scratch%n_4_3*functional%K )
    functional%dsx_dn   = -0.5_real_8*( four_thirds*scratch%n_1_3*functional%K + scratch%n_4_3*functional%dK_dn )
    functional%dsx_dg   = -0.5_real_8*( scratch%n_4_3*functional%dK_dg )
    !

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_gga_x_pbe
  ! ==================================================================
  PURE SUBROUTINE cp_spin_gga_x_revpbe( scratch,functional )
    ! ==--------------------------------------------------------------==


    TYPE(cp_xc_scratch_t), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(in)                             :: scratch
    TYPE(cp_xc_functional_t), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(inout)                          :: functional

    IF (scratch(a)%n >= 1.0e-10_real_8) &
         CALL cp_gga_x_revpbe( scratch(a), functional(a) )
    IF (scratch(b)%n >= 1.0e-10_real_8) &
         CALL cp_gga_x_revpbe( scratch(b), functional(b) )

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_spin_gga_x_revpbe
  ! ==================================================================
  PURE SUBROUTINE cp_gga_x_revpbe( scratch,functional )
    ! ==--------------------------------------------------------------==
    ! Rewritten from scratch (not ported from oldcode)
    ! Output term by term identical to libxc
    !                                     09.09.2016 M.P. Bircher @ EPFL
    ! ==--------------------------------------------------------------==


    TYPE(cp_xc_scratch_t), INTENT(in)        :: scratch
    TYPE(cp_xc_functional_t), INTENT(inout)  :: functional

    REAL(real_8), PARAMETER :: kappa = 1.2450_real_8, &
      mu = 0.2195149727645171_real_8, mu_by_kappa = mu/kappa , &
      us = 0.161620459673995492_real_8

    REAL(real_8)                             :: ddenominator_ds, denominator, &
                                                denominator_2, dK_ds, ds_dg, &
                                                ds_dn, Fx, lda_K, s, s_2

!
! UEG
!

    lda_K            = two_to_four_thirds*slater_prefactor*func2%salpha
    !
    ! Gradient correction
    !
    s                = us/two_to_one_third*scratch%abs_g / scratch%n_4_3
    s_2              = s*s
    !
    ! Work with the spin-polarized enhancement factor (cf. PBE paper, p. 3867)
    ! Prefactor incorporated into UEG
    !
    denominator      = 1.0_real_8 + mu_by_kappa*s_2
    denominator_2    = denominator*denominator
    ddenominator_ds  = 2.0_real_8*mu_by_kappa*s
    !
    Fx               = ( 1.0_real_8 + kappa - kappa/denominator )
    functional%K     = lda_K*Fx
    !
    dK_ds            = lda_K*( ddenominator_ds*kappa/denominator_2 )
    !
    ds_dn            = -four_thirds*us/two_to_one_third*scratch%abs_g / (scratch%n_4_3*scratch%n)
    ds_dg            = us /( two_to_one_third * scratch%n_4_3 )
    !
    functional%dK_dn = dK_ds*ds_dn
    functional%dK_dg = dK_ds*ds_dg
    !
    functional%sx_sigma = -0.5_real_8*( scratch%n_4_3*functional%K )
    functional%dsx_dn   = -0.5_real_8*( four_thirds*scratch%n_1_3*functional%K + scratch%n_4_3*functional%dK_dn )
    functional%dsx_dg   = -0.5_real_8*( scratch%n_4_3*functional%dK_dg )
    !

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_gga_x_revpbe
  ! ==================================================================
  PURE SUBROUTINE cp_spin_gga_x_pbe_sol( scratch,functional )
    ! ==--------------------------------------------------------------==


    TYPE(cp_xc_scratch_t), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(in)                             :: scratch
    TYPE(cp_xc_functional_t), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(inout)                          :: functional

    IF (scratch(a)%n >= 1.0e-10_real_8) &
         CALL cp_gga_x_pbe_sol( scratch(a), functional(a) )
    IF (scratch(b)%n >= 1.0e-10_real_8) &
         CALL cp_gga_x_pbe_sol( scratch(b), functional(b) )

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_spin_gga_x_pbe_sol
  ! ==================================================================
  PURE SUBROUTINE cp_gga_x_pbe_sol( scratch,functional )
    ! ==--------------------------------------------------------------==
    ! Rewritten from scratch (not ported from oldcode)
    ! Output term by term identical to libxc
    !                                     09.09.2016 M.P. Bircher @ EPFL
    ! ==--------------------------------------------------------------==


    TYPE(cp_xc_scratch_t), INTENT(in)        :: scratch
    TYPE(cp_xc_functional_t), INTENT(inout)  :: functional

    REAL(real_8), PARAMETER :: kappa = 0.8040_real_8, &
      mu = 0.123456790123456789_real_8, mu_by_kappa = mu/kappa , &
      us = 0.161620459673995492_real_8

    REAL(real_8)                             :: ddenominator_ds, denominator, &
                                                denominator_2, dK_ds, ds_dg, &
                                                ds_dn, Fx, lda_K, s, s_2

!
! UEG
!

    lda_K            = two_to_four_thirds*slater_prefactor*func2%salpha
    !
    ! Gradient correction
    !
    s                = us/two_to_one_third*scratch%abs_g / scratch%n_4_3
    s_2              = s*s
    !
    ! Work with the spin-polarized enhancement factor (cf. PBE paper, p. 3867)
    ! Prefactor incorporated into UEG
    !
    denominator      = 1.0_real_8 + mu_by_kappa*s_2
    denominator_2    = denominator*denominator
    ddenominator_ds  = 2.0_real_8*mu_by_kappa*s
    !
    Fx               = ( 1.0_real_8 + kappa - kappa/denominator )
    functional%K     = lda_K*Fx
    !
    dK_ds            = lda_K*( ddenominator_ds*kappa/denominator_2 )
    !
    ds_dn            = -four_thirds*us/two_to_one_third*scratch%abs_g / (scratch%n_4_3*scratch%n)
    ds_dg            = us /( two_to_one_third * scratch%n_4_3 )
    !
    functional%dK_dn = dK_ds*ds_dn
    functional%dK_dg = dK_ds*ds_dg
    !
    functional%sx_sigma = -0.5_real_8*( scratch%n_4_3*functional%K )
    functional%dsx_dn   = -0.5_real_8*( four_thirds*scratch%n_1_3*functional%K + scratch%n_4_3*functional%dK_dn )
    functional%dsx_dg   = -0.5_real_8*( scratch%n_4_3*functional%dK_dg )
    !

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_gga_x_pbe_sol
  ! ==================================================================
  PURE SUBROUTINE cp_spin_gga_x_pbe_flex( scratch,functional )
    ! ==--------------------------------------------------------------==


    TYPE(cp_xc_scratch_t), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(in)                             :: scratch
    TYPE(cp_xc_functional_t), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(inout)                          :: functional

    IF (scratch(a)%n >= 1.0e-10_real_8) &
         CALL cp_gga_x_pbe_flex( scratch(a), functional(a) )
    IF (scratch(b)%n >= 1.0e-10_real_8) &
         CALL cp_gga_x_pbe_flex( scratch(b), functional(b) )

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_spin_gga_x_pbe_flex
  ! ==================================================================
  PURE SUBROUTINE cp_gga_x_pbe_flex( scratch,functional )
    ! ==--------------------------------------------------------------==
    ! Rewritten from scratch (not ported from oldcode)
    ! Output term by term identical to libxc
    !                                     09.09.2016 M.P. Bircher @ EPFL
    ! ==--------------------------------------------------------------==


    TYPE(cp_xc_scratch_t), INTENT(in)        :: scratch
    TYPE(cp_xc_functional_t), INTENT(inout)  :: functional

    REAL(real_8), PARAMETER :: us = 0.161620459673995492_real_8

    REAL(real_8) :: ddenominator_ds, denominator, denominator_2, dK_ds, &
      ds_dg, ds_dn, Fx, kappa, lda_K, mu, mu_by_kappa, s, s_2

    kappa            = cp_gga_x_param%pbe%kappa
    mu               = cp_gga_x_param%pbe%mu
    mu_by_kappa      = mu/kappa
    !
    ! UEG
    !
    lda_K            = two_to_four_thirds*slater_prefactor*func2%salpha
    !
    ! Gradient correction
    !
    s                = us/two_to_one_third*scratch%abs_g / scratch%n_4_3
    s_2              = s*s
    !
    ! Work with the spin-polarized enhancement factor (cf. PBE paper, p. 3867)
    ! Prefactor incorporated into UEG
    !
    denominator      = 1.0_real_8 + mu_by_kappa*s_2
    denominator_2    = denominator*denominator
    ddenominator_ds  = 2.0_real_8*mu_by_kappa*s
    !
    Fx               = ( 1.0_real_8 + kappa - kappa/denominator )
    functional%K     = lda_K*Fx
    !
    dK_ds            = lda_K*( ddenominator_ds*kappa/denominator_2 )
    !
    ds_dn            = -four_thirds*us/two_to_one_third*scratch%abs_g / (scratch%n_4_3*scratch%n)
    ds_dg            = us /( two_to_one_third * scratch%n_4_3 )
    !
    functional%dK_dn = dK_ds*ds_dn
    functional%dK_dg = dK_ds*ds_dg
    !
    functional%sx_sigma = -0.5_real_8*( scratch%n_4_3*functional%K )
    functional%dsx_dn   = -0.5_real_8*( four_thirds*scratch%n_1_3*functional%K + scratch%n_4_3*functional%dK_dn )
    functional%dsx_dg   = -0.5_real_8*( scratch%n_4_3*functional%dK_dg )
    !

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_gga_x_pbe_flex
  ! ==================================================================
  PURE SUBROUTINE cp_spin_gga_x_optx( scratch,functional )
    ! ==--------------------------------------------------------------==


    TYPE(cp_xc_scratch_t), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(in)                             :: scratch
    TYPE(cp_xc_functional_t), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(inout)                          :: functional

    IF (scratch(a)%n >= 1.0e-10_real_8) &
         CALL cp_gga_x_optx( scratch(a), functional(a) )
    IF (scratch(b)%n >= 1.0e-10_real_8) &
         CALL cp_gga_x_optx( scratch(b), functional(b) )

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_spin_gga_x_optx
  ! ==================================================================
  PURE SUBROUTINE cp_gga_x_optx( scratch,functional )
    ! ==--------------------------------------------------------------==
    ! Rewritten from scratch
    !                                     09.09.2016 M.P. Bircher @ EPFL
    ! ==--------------------------------------------------------------==

    TYPE(cp_xc_scratch_t), INTENT(in)        :: scratch
    TYPE(cp_xc_functional_t), INTENT(inout)  :: functional

    REAL(real_8), PARAMETER :: optx_c_u = 1.43169_real_8, &
      optx_c_x = 1.05151_real_8, optx_gamma = 0.00600_real_8

    REAL(real_8) :: denominator, doptx_dx, doptx_K_dx, dx_dg, dx_dn, gamma_x, &
      gamma_x_2, lda_K, optx_K, u_optx, u_optx_2, x

!
! UEG
!

    lda_K              = optx_c_x*two_to_four_thirds*slater_prefactor*func2%salpha

    !
    ! OPTX
    !
    x                  = scratch%abs_g / scratch%n_4_3
    dx_dn              = -four_thirds*scratch%abs_g / (scratch%n_4_3*scratch%n)
    dx_dg              = 1.0_real_8 / scratch%n_4_3
    !
    gamma_x            = optx_gamma*x
    gamma_x_2          = gamma_x*x
    denominator        = 1.0_real_8 + gamma_x_2
    u_optx             = gamma_x_2 / denominator
    u_optx_2           = u_optx*u_optx
    optx_K             = 2.0_real_8*optx_c_u*u_optx_2
    !
    ! Formally:
    ! dgamma_x_2_dx      = 2.0_real_8*gamma_x
    ! ddenominator_dx    = -dgamma_x_2_dx/(denominator*denominator)
    ! doptx_dx           = dgamma_x_2_dx / denominator +  gamma_x_2*ddenominator_dx
    ! Simplifies to:
    doptx_dx           = 2.0_real_8*gamma_x / denominator * (1.0_real_8 - u_optx)

    ! Formally:                         
    ! doptx_2_dx         = 2.0_real_8*doptx_dx*u_optx
    ! doptx_K_dx         = 2.0_real_8*optx_c_u*doptx_2_dx
    ! Simplifies to:
    doptx_K_dx         = 4.0_real_8*optx_c_u*doptx_dx*u_optx
    !
    functional%dK_dn   = doptx_K_dx*dx_dn
    functional%dK_dg   = doptx_K_dx*dx_dg
    !
    functional%K       = ( lda_K + optx_K )
    !
    ! Proper conventions for readability
    !
    functional%sx_sigma = -0.5_real_8*( scratch%n_4_3*functional%K )
    functional%dsx_dn   = -0.5_real_8*( four_thirds*scratch%n_1_3*functional%K + scratch%n_4_3*functional%dK_dn )
    functional%dsx_dg   = -0.5_real_8*( scratch%n_4_3*functional%dK_dg )

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_gga_x_optx
  ! ==================================================================
  PURE SUBROUTINE cp_spin_gga_x_hcth( scratch,functional )
    ! ==--------------------------------------------------------------==


    TYPE(cp_xc_scratch_t), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(in)                             :: scratch
    TYPE(cp_xc_functional_t), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(inout)                          :: functional

    IF (scratch(a)%n >= 1.0e-10_real_8) &
         CALL cp_gga_x_hcth( scratch(a), functional(a) )
    IF (scratch(b)%n >= 1.0e-10_real_8) &
         CALL cp_gga_x_hcth( scratch(b), functional(b) )

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_spin_gga_x_hcth
  ! ==================================================================
  PURE SUBROUTINE cp_gga_x_hcth( scratch,functional )
    ! ==--------------------------------------------------------------==
    ! Rewritten from scratch.
    !                                     19.07.2017 M.P. Bircher @ EPFL
    ! ==--------------------------------------------------------------==

    TYPE(cp_xc_scratch_t), INTENT(IN)       :: scratch
    TYPE(cp_xc_functional_t), INTENT(INOUT) :: functional

    INTEGER, PARAMETER                      :: gm = hcth_gm

    REAL(real_8), PARAMETER                 :: sqrt_small = 1.0e-12_real_8

    REAL(real_8)                            :: lda_K, hcth_K
    REAL(real_8), DIMENSION(gm)             :: hcth_g, u_hcth, dhcth_u_ds
    REAL(real_8)                            :: abs_g, s, ds_dn, ds_dg
    REAL(real_8)                            :: gamma_s, gamma_s_2, denominator
    REAL(real_8)                            :: dhcth_1_ds, dhcth_K_ds


    abs_g              = max(scratch%abs_g,sqrt_small)

    !
    ! UEG
    !
    lda_K              = two_to_four_thirds*slater_prefactor*func2%salpha

    !
    ! HCTH / B97
    !
    s                  = abs_g / scratch%n_4_3
    ds_dn              = -four_thirds* s / scratch%n
    ds_dg              = 1.0_real_8 / scratch%n_4_3
    !
    gamma_s            = cp_gga_x_param%hcth%gamma*s
    gamma_s_2          = gamma_s*s
    denominator        = 1.0_real_8 + gamma_s_2
    !
    u_hcth(1)          = gamma_s_2 / denominator
    u_hcth(2)          = u_hcth(1)*u_hcth(1)
    u_hcth(3)          = u_hcth(1)*u_hcth(2)
    u_hcth(4)          = u_hcth(2)*u_hcth(2)
    !
    hcth_g(:)          = cp_gga_x_param%hcth%c(:)*u_hcth(:)
    hcth_K             = lda_K*( cp_gga_x_param%hcth%c_0 + sum(hcth_g(:)) )
    !
    ! dhcth_0_ds       = 0
    dhcth_1_ds         = 2.0_real_8 * gamma_s / denominator * (1.0_real_8 - u_hcth(1))
    dhcth_u_ds(1)      =            dhcth_1_ds
    dhcth_u_ds(2)      = 2.0_real_8*dhcth_1_ds*u_hcth(1)
    dhcth_u_ds(3)      = 3.0_real_8*dhcth_1_ds*u_hcth(2)
    dhcth_u_ds(4)      = 4.0_real_8*dhcth_1_ds*u_hcth(3)
    !
    dhcth_K_ds         = lda_K*sum(cp_gga_x_param%hcth%c(:)*dhcth_u_ds(:))
    !
    functional%dK_dn   = dhcth_K_ds*ds_dn
    functional%dK_dg   = dhcth_K_ds*ds_dg
    !
    functional%K       = ( hcth_K )
    !
    ! Proper conventions for readability
    !
    functional%sx_sigma = -0.5_real_8*( scratch%n_4_3*functional%K )
    functional%dsx_dn   = -0.5_real_8*( four_thirds*scratch%n_1_3*functional%K + scratch%n_4_3*functional%dK_dn )
    functional%dsx_dg   = -0.5_real_8*( scratch%n_4_3*functional%dK_dg )

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_gga_x_hcth
  ! ==================================================================

  ! ==================================================================
  ! Coulomb attenuation method for LDA & GGA exchange functionals
  ! ==================================================================

  ! ==================================================================
  PURE SUBROUTINE cp_gga_screening_cam( scratch,functional )
    ! ==--------------------------------------------------------------==
    ! Screening routine (spin-dependent)
    !                                     09.09.2016 M.P. Bircher @ EPFL
    ! ==--------------------------------------------------------------==


    TYPE(cp_xc_scratch_t), INTENT(in)        :: scratch
    TYPE(cp_xc_functional_t), INTENT(inout)  :: functional

    REAL(real_8) :: a_times_b_minus_c, cam, cam_dsx_dg, cam_dsx_dn, cam_sx, &
      da_dg, da_dn, db_da, dc_da, dcam_da, dcam_dg, dcam_dn, derf_da, dk_dg, &
      dk_dn, erf_of_a, grad_dependent_attenuation, sqrt_K, term_a, &
      term_a_squared, term_b, term_bexp, term_c

!

    sqrt_K         = SQRT(functional%K)
    !
    term_a         = func2%srxa*sqrt_K/(6._real_8*sqrt_of_pi*scratch%n_1_3)
    term_a_squared = term_a*term_a
    !
    IF (term_a < 1.0e6_real_8) THEN
       term_bexp = EXP(-1._real_8/(4._real_8*term_a_squared))
       term_b    = term_bexp - 1._real_8
    ELSE
       term_b    = -1.0_real_8/(4.0_real_8*term_a_squared) 
       term_bexp = EXP(-1._real_8/(4._real_8*term_a_squared))
       ! term_bexp = term_b + 1.0_real_8
    END IF
    !
    term_c = 2._real_8*term_a_squared*term_b + 0.5_real_8
    !
    erf_of_a = sqrt_of_pi*ERF(0.5_real_8/term_a)
    !
    a_times_b_minus_c = 2.0_real_8*term_a*(term_b-term_c)
    !
    grad_dependent_attenuation  = 8._real_8/3._real_8*term_a*(erf_of_a + a_times_b_minus_c)
    cam = 1._real_8 - func2%cam_alpha - func2%cam_beta*grad_dependent_attenuation
    !
    dk_dn   = functional%dK_dn
    dk_dg   = functional%dK_dg
    !
    derf_da = -1.0_real_8*term_bexp/term_a_squared
    db_da   = 0.5_real_8/(term_a_squared*term_a)*term_bexp
    dc_da   = 4._real_8*term_a*term_b + 2._real_8*term_a_squared*db_da
    !
    dcam_da = -8.0_real_8/3.0_real_8*func2%cam_beta &
         *( (erf_of_a + a_times_b_minus_c) + &
         term_a*(derf_da + (2.0_real_8*(term_b - term_c) + 2.0_real_8*term_a*(db_da -dc_da) ) ) )
    !
    !
    da_dn = func2%srxa/(6._real_8*sqrt_of_pi)*( 0.5_real_8*dk_dn/sqrt_K / scratch%n_1_3 &
         - one_third*sqrt_K/scratch%n_4_3 ) 
    da_dg = func2%srxa/(6._real_8*sqrt_of_pi*scratch%n_1_3) * ( 0.5_real_8*dk_dg/sqrt_K )
    !
    !
    dcam_dn = dcam_da*da_dn 
    dcam_dg = dcam_da*da_dg
    !
    !
    cam_sx     = cam*(functional%sx_sigma)
    cam_dsx_dn = cam*functional%dsx_dn + dcam_dn*functional%sx_sigma
    cam_dsx_dg = cam*functional%dsx_dg + functional%sx_sigma*dcam_dg
    !
    ! overwrite...
    !
    functional%sx_sigma = cam_sx
    functional%dsx_dn   = cam_dsx_dn
    functional%dsx_dg   = cam_dsx_dg
    !

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_gga_screening_cam
  ! ==================================================================

  ! ==================================================================
  ! Check parameter compatibility
  ! ==================================================================

  ! ==================================================================
  ELEMENTAL FUNCTION cp_gga_x_check(tag) &
  RESULT (OK)
    ! ==--------------------------------------------------------------==
    ! Checks whether parameters are properly initialised and reasonable
    !
    !                                 25.07.2017 M.P.Bircher @ LCBC/EPFL
    ! ==--------------------------------------------------------------==

    CHARACTER(len=*), INTENT(in)             :: tag
    LOGICAL                                  :: OK 

    OK = (func2%salpha /= 0.0_real_8)

    SELECT CASE( tag )
    !
    ! Functionals with flexible parameters
    !
    CASE( "CP_GGA_X_HCTH","CP_SPIN_GGA_X_HCTH" )
      OK = ( cp_gga_x_param%init .and. OK )
    CASE( "CP_GGA_X_PBE_FLEX","CP_SPIN_GGA_X_PBE_FLEX" )
      OK = ( (cp_gga_x_param%init) .and. OK )
    !
    ! Functionals with hard-coded parameters
    !
    CASE( "CP_GGA_X_B88",  "CP_SPIN_GGA_X_B88" )
      OK = ( (func2%bbeta  /= 0.0_real_8) .and. OK )
    CASE( "CP_GGA_X_OPTX", "CP_SPIN_GGA_X_OPTX",&
          "CP_GGA_X_PBE",  "CP_SPIN_GGA_X_PBE",&
          "CP_GGA_X_REVPBE",  "CP_SPIN_GGA_C_REVPBE",&
          "CP_GGA_X_PBE_SOL", "CP_SPIN_GGA_X_PBE_SOL" )
      OK = OK
    !
    ! Any other tag
    !
    CASE DEFAULT
      OK = .false.
    END SELECT

    ! ==--------------------------------------------------------------==
  END FUNCTION cp_gga_x_check
  ! ==================================================================
END MODULE cp_gga_exchange_utils
