! ==================================================================
! Provides: - Generic types for xc functionals used in xc_functionals
!             and xc_driver, storing enhancment factors and
!             their derivatives (used for Coulomb attenuation)
!           - Scratch type for storage of frequently used powers of
!             rho and the reduced gradient
!           
! analytical linres deriv.    19.06.2017 - M. P. Bircher @ LCBC/EPFL
! included meta functionals   18.04.2017 - M. P. Bircher @ LCBC/EPFL
!                             02.10.2016 - M. P. Bircher @ LCBC/EPFL
! ==================================================================
MODULE cpfunc_types
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  PRIVATE

  !
  ! Parameters (important for spin-dependent case)
  !
  INTEGER, PARAMETER, PUBLIC :: cp_xc_spin_pairs      = 2 ! a, b
  INTEGER, PARAMETER, PUBLIC :: cp_xc_spin_components = 3 ! a, b, ab
  INTEGER, PARAMETER, PUBLIC :: cp_xc_abs_xyz         = 4 ! cartesian components + norm
  INTEGER, PARAMETER, PUBLIC :: cp_xc_xyz             = 3 ! cartesian components
  !
  ! Work type
  !
  TYPE, PUBLIC ::  cp_xc_scratch_t
     !
     ! n: rho
     REAL(real_8) :: n 
     REAL(real_8) :: n_1_3
     REAL(real_8) :: n_4_3
     !
     ! g: gradient
     REAL(real_8) :: g_2
     REAL(real_8) :: abs_g
     !
     ! tau (t): kinetic energy density
     REAL(real_8) :: tau
  END TYPE cp_xc_scratch_t
  !
  ! Work type for analytical derivatives
  !
  TYPE, PUBLIC ::  cp_dxc_scratch_t
     !
     ! n: rho
     REAL(real_8) :: n 
     REAL(real_8) :: dn 
     REAL(real_8) :: g_2
     REAL(real_8) :: gx_2
     REAL(real_8) :: gy_2
     REAL(real_8) :: gz_2
     REAL(real_8) :: dgx
     REAL(real_8) :: dgy
     REAL(real_8) :: dgz
  END TYPE cp_dxc_scratch_t

  !
  ! Energy density and derivatives
  !
  TYPE, PUBLIC :: cp_xc_functional_t
     !
     ! Exchange
     REAL(real_8) :: sx_sigma
     REAL(real_8) :: K
     REAL(real_8) :: dK_dn
     REAL(real_8) :: dK_dg
     REAL(real_8) :: dsx_dn
     REAL(real_8) :: dsx_dg
     REAL(real_8) :: dsx_dt
     !
     ! Correlation
     REAL(real_8) :: sc
     REAL(real_8) :: v1c
     REAL(real_8) :: v2c
     REAL(real_8) :: vtc
  END TYPE cp_xc_functional_t
  !
  ! Linres derivatives
  ! (lazy rewrite, kept old names from dxc_ana)
  !
  TYPE, PUBLIC :: cp_dxc_functional_t
     !
     ! Exchange-correlation
     REAL(real_8) :: tmp1
     REAL(real_8) :: tmp2
     REAL(real_8) :: tmp3
     REAL(real_8) :: tmp4
     REAL(real_8) :: tmp5
  END TYPE cp_dxc_functional_t

END MODULE cpfunc_types
