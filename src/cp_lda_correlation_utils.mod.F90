! ==================================================================
! Provides: - LDA correlation, adapted from OLDCODE
!
! Input:  Spin-polarised density (1/2 n for closed shell systems)
! Output: For CLOSED-SHELL system for spin-unpolarised functionals
!         For OPEN-SHELL system for the SPIN functionals
!
!                             24.03.2017 - M. P. Bircher @ LCBC/EPFL
! ==================================================================
MODULE cp_lda_correlation_utils
  USE cnst,                            ONLY: pi
  USE cpfunc_types,                    ONLY: cp_xc_functional_t,&
                                             cp_xc_scratch_t,&
                                             cp_xc_spin_components
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER      :: alpha                  = 1
  INTEGER, PARAMETER      :: beta                   = 2
  INTEGER, PARAMETER      :: alpha_beta             = 3

  REAL(real_8), PARAMETER :: three_quarters_by_pi   = 0.75_real_8/pi
  REAL(real_8), PARAMETER :: tqbp_to_one_third      = three_quarters_by_pi**(1.0_real_8/3.0_real_8)
  REAL(real_8), PARAMETER :: two_to_one_third       = 2.0_real_8**(1.0_real_8/3.0_real_8)
  REAL(real_8), PARAMETER :: one_third              = 1._real_8/3._real_8
  REAL(real_8), PARAMETER :: two_thirds             = 2._real_8/3._real_8
  REAL(real_8), PARAMETER :: four_thirds            = 4._real_8/3._real_8
  REAL(real_8), PARAMETER :: two_to_four_thirds     = 2._real_8**(4._real_8/3._real_8)

  !
  ! Functional wrapping routines
  !
  PUBLIC :: CP_LDA_C_VWN
  PUBLIC :: CP_LDA_C_PZ
  PUBLIC :: CP_LDA_C_PW
  PUBLIC :: CP_LDA_C_OB_PW
  PUBLIC :: CP_SPIN_LDA_C_VWN
  PUBLIC :: CP_SPIN_LDA_C_PZ
  PUBLIC :: CP_SPIN_LDA_C_PW
  PUBLIC :: CP_SPIN_LDA_C_OB_PW

  !
  ! Check
  !
  PUBLIC :: CP_LDA_C_CHECK

  !
  ! These routines have to be available to cp_gga_correlation_utils for PBE
  !
  PUBLIC :: CP_LDA_C_PW_P
  PUBLIC :: CP_LDA_C_PW_F
  PUBLIC :: CP_LDA_C_PW_alpha_fit
  !
  ! These are for TPSS
  !
  PUBLIC :: CP_LDA_C_OB_PW_P
  PUBLIC :: CP_LDA_C_PZ_P
  PUBLIC :: CP_LDA_C_VWN_P

  PUBLIC :: zeta_of
  PUBLIC :: f_of
  PUBLIC :: df_dz_of

  !
  ! Routines called from the wrappers
  !

  PRIVATE :: CP_LDA_C_VWN_F
  PRIVATE :: CP_LDA_C_VWN_alpha_fit
  PRIVATE :: CP_LDA_C_PZ_F
  PRIVATE :: CP_LDA_C_OB_PW_F

CONTAINS

  ! ==================================================================
  ! Wrapping routines
  ! ==================================================================

  ! ==================================================================
  PURE SUBROUTINE cp_spin_lda_c_vwn(scratch,functional)
    ! ==--------------------------------------------------------------==
    ! Ported from OLDCODE
    !                                            M.P.Bircher @ LCBC/EPFL
    ! ==--------------------------------------------------------------==

    TYPE(cp_xc_scratch_t), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(in)                             :: scratch
    TYPE(cp_xc_functional_t), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(inout)                          :: functional

    REAL(real_8), PARAMETER :: &
      fpp = 4.0_real_8 / (9.0_real_8 * (two_to_one_third - 1.0_real_8) )

    REAL(real_8) :: dDFP_dnab, df_of_zeta_dna, df_of_zeta_dnb, df_of_zeta_dz, &
      DFP, dipo_alpha_dna, dipo_alpha_dnb, dipo_delta_dna, dipo_delta_dnb, &
      dzeta_4_dna, dzeta_4_dnb, dzeta_dna, dzeta_dnb, f_of_zeta, ipo_alpha, &
      ipo_delta, n_ab, nab_2, rs_ab, scvwn, zeta, zeta_2, zeta_3, zeta_4
    TYPE(cp_xc_functional_t)                 :: vwn_alpha, vwn_F, vwn_P

! \Delta (e_F e_P)

    n_ab  = scratch(alpha_beta)%n
    rs_ab = tqbp_to_one_third/(scratch(alpha_beta)%n_1_3)

    CALL cp_lda_c_vwn_P( n_ab, rs_ab, vwn_P )
    CALL cp_lda_c_vwn_F( n_ab, rs_ab, vwn_F )
    CALL cp_lda_c_vwn_alpha_fit( n_ab, rs_ab, vwn_alpha )

    zeta            = zeta_of( scratch(:) )
    zeta_2          = zeta   * zeta
    zeta_3          = zeta_2 * zeta
    zeta_4          = zeta_2 * zeta_2

    nab_2           = scratch(alpha_beta)%n * scratch(alpha_beta)%n 

    dzeta_dna       =  2.0_real_8 * scratch(beta)%n  / nab_2
    dzeta_dnb       = -2.0_real_8 * scratch(alpha)%n / nab_2

    dzeta_4_dna     = 4.0_real_8*zeta_3*dzeta_dna
    dzeta_4_dnb     = 4.0_real_8*zeta_3*dzeta_dnb

    f_of_zeta       = f_of( zeta )
    df_of_zeta_dz   = df_dz_of( zeta )
    df_of_zeta_dna  = df_of_zeta_dz * dzeta_dna
    df_of_zeta_dnb  = df_of_zeta_dz * dzeta_dnb

    DFP             = vwn_F%sc  - vwn_P%sc
    !
    ! No explicit a, b dependency
    !
    dDFP_dnab       = vwn_F%v1c - vwn_P%v1c

    ipo_alpha       = f_of_zeta / fpp * (1.0_real_8 - zeta_4)
    dipo_alpha_dna  = df_of_zeta_dna / fpp * (1.0_real_8 - zeta_4) &
         - f_of_zeta / fpp * (dzeta_4_dna)
    dipo_alpha_dnb  = df_of_zeta_dnb / fpp * (1.0_real_8 - zeta_4) &
         - f_of_zeta / fpp * (dzeta_4_dnb)

    ipo_delta       = f_of_zeta * zeta_4
    dipo_delta_dna  = df_of_zeta_dna * zeta_4 + f_of_zeta * dzeta_4_dna
    dipo_delta_dnb  = df_of_zeta_dnb * zeta_4 + f_of_zeta * dzeta_4_dnb

    scvwn           = vwn_P%sc + ipo_delta*DFP + ipo_alpha*vwn_alpha%sc

    functional(alpha)%sc       = 0.0_real_8      
    functional(beta)%sc        = 0.0_real_8      
    functional(alpha_beta)%sc  = scvwn

    functional(alpha)%v1c      = vwn_P%v1c + dipo_delta_dna*DFP + ipo_delta*dDFP_dnab &
         + dipo_alpha_dna*vwn_alpha%sc  + ipo_alpha*vwn_alpha%v1c
    functional(beta)%v1c       = vwn_P%v1c + dipo_delta_dnb*DFP + ipo_delta*dDFP_dnab &
         + dipo_alpha_dnb*vwn_alpha%sc  + ipo_alpha*vwn_alpha%v1c
    functional(alpha_beta)%v1c = 0.0_real_8

    functional(alpha)%v2c      = 0.0_real_8
    functional(beta)%v2c       = 0.0_real_8 
    functional(alpha_beta)%v2c = 0.0_real_8

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_spin_lda_c_vwn
  ! ==================================================================
  PURE SUBROUTINE cp_lda_c_vwn(scratch,functional)
    ! ==--------------------------------------------------------------==
    ! Ported from OLDCODE
    !                                            M.P.Bircher @ LCBC/EPFL
    ! ==--------------------------------------------------------------==

    TYPE(cp_xc_scratch_t), INTENT(in)        :: scratch
    TYPE(cp_xc_functional_t), INTENT(inout)  :: functional

    REAL(real_8)                             :: n_ab, rs_ab

    n_ab  = 2.0_real_8*scratch%n
    rs_ab = tqbp_to_one_third/(two_to_one_third*scratch%n_1_3)

    CALL cp_lda_c_vwn_P( n_ab, rs_ab, functional )

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_lda_c_vwn
  ! ==================================================================
  PURE SUBROUTINE cp_spin_lda_c_pw(scratch,functional)
    ! ==--------------------------------------------------------------==
    ! Rewritten from scratch
    !                                            M.P.Bircher @ LCBC/EPFL
    ! ==--------------------------------------------------------------==

    TYPE(cp_xc_scratch_t), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(in)                             :: scratch
    TYPE(cp_xc_functional_t), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(inout)                          :: functional

    REAL(real_8), PARAMETER                  :: fpp = 1.709921_real_8

    REAL(real_8) :: dDFP_dnab, df_of_zeta_dna, df_of_zeta_dnb, df_of_zeta_dz, &
      DFP, dipo_alpha_dna, dipo_alpha_dnb, dipo_delta_dna, dipo_delta_dnb, &
      dzeta_4_dna, dzeta_4_dnb, dzeta_dna, dzeta_dnb, f_of_zeta, ipo_alpha, &
      ipo_delta, n_ab, nab_2, rs_ab, scpw, zeta, zeta_2, zeta_3, zeta_4
    TYPE(cp_xc_functional_t)                 :: pw_alpha, pw_F, pw_P

! \Delta (e_F - e_P)

    n_ab  = scratch(alpha_beta)%n
    rs_ab = tqbp_to_one_third/(scratch(alpha_beta)%n_1_3)

    CALL cp_lda_c_pw_P( n_ab, rs_ab, pw_P )
    CALL cp_lda_c_pw_F( n_ab, rs_ab, pw_F )
    CALL cp_lda_c_pw_alpha_fit( n_ab, rs_ab, pw_alpha )

    zeta            = zeta_of( scratch(:) )
    zeta_2          = zeta   * zeta
    zeta_3          = zeta_2 * zeta
    zeta_4          = zeta_2 * zeta_2

    nab_2           = n_ab*n_ab 

    dzeta_dna       =  2.0_real_8 * scratch(beta)%n  / nab_2
    dzeta_dnb       = -2.0_real_8 * scratch(alpha)%n / nab_2

    dzeta_4_dna     = 4.0_real_8*zeta_3*dzeta_dna
    dzeta_4_dnb     = 4.0_real_8*zeta_3*dzeta_dnb

    f_of_zeta       = f_of( zeta )
    df_of_zeta_dz   = df_dz_of( zeta )
    df_of_zeta_dna  = df_of_zeta_dz * dzeta_dna
    df_of_zeta_dnb  = df_of_zeta_dz * dzeta_dnb

    DFP             = pw_F%sc  - pw_P%sc
    !
    ! No explicit a, b dependency
    !
    dDFP_dnab       = pw_F%v1c - pw_P%v1c

    !
    ! Alpha takes a negative sign; fix this by changing the sign of ipo_alpha
    !
    ipo_alpha       = - ( f_of_zeta / fpp * (1.0_real_8 - zeta_4) )
    dipo_alpha_dna  = - ( df_of_zeta_dna / fpp * (1.0_real_8 - zeta_4) &
         - f_of_zeta / fpp * (dzeta_4_dna) )
    dipo_alpha_dnb  = - ( df_of_zeta_dnb / fpp * (1.0_real_8 - zeta_4) &
         - f_of_zeta / fpp * (dzeta_4_dnb) )

    ipo_delta       = f_of_zeta * zeta_4
    dipo_delta_dna  = df_of_zeta_dna * zeta_4 + f_of_zeta * dzeta_4_dna
    dipo_delta_dnb  = df_of_zeta_dnb * zeta_4 + f_of_zeta * dzeta_4_dnb

    scpw            = pw_P%sc + ipo_delta*DFP + ipo_alpha*pw_alpha%sc

    functional(alpha)%sc       = 0.0_real_8      
    functional(beta)%sc        = 0.0_real_8      
    functional(alpha_beta)%sc  = scpw

    functional(alpha)%v1c      = pw_P%v1c + dipo_delta_dna*DFP + ipo_delta*dDFP_dnab &
         + dipo_alpha_dna*pw_alpha%sc  + ipo_alpha*pw_alpha%v1c
    functional(beta)%v1c       = pw_P%v1c + dipo_delta_dnb*DFP + ipo_delta*dDFP_dnab &
         + dipo_alpha_dnb*pw_alpha%sc  + ipo_alpha*pw_alpha%v1c
    functional(alpha_beta)%v1c = 0.0_real_8

    functional(alpha)%v2c      = 0.0_real_8
    functional(beta)%v2c       = 0.0_real_8 
    functional(alpha_beta)%v2c = 0.0_real_8

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_spin_lda_c_pw
  ! ==================================================================
  PURE SUBROUTINE cp_lda_c_pw(scratch,functional)
    ! ==--------------------------------------------------------------==
    ! Rewritten from scratch
    !                                            M.P.Bircher @ LCBC/EPFL
    ! ==--------------------------------------------------------------==


    TYPE(cp_xc_scratch_t), INTENT(in)        :: scratch
    TYPE(cp_xc_functional_t), INTENT(inout)  :: functional

    REAL(real_8)                             :: n_ab, rs_ab

    n_ab  = 2.0_real_8*scratch%n
    rs_ab = tqbp_to_one_third/(two_to_one_third*scratch%n_1_3)

    CALL cp_lda_c_pw_P( n_ab, rs_ab, functional )

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_lda_c_pw
  ! ==================================================================
  PURE SUBROUTINE cp_spin_lda_c_ob_pw(scratch,functional)
    ! ==--------------------------------------------------------------==
    ! Rewritten from scratch
    !                                            M.P.Bircher @ LCBC/EPFL
    ! ==--------------------------------------------------------------==

    TYPE(cp_xc_scratch_t), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(in)                             :: scratch
    TYPE(cp_xc_functional_t), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(inout)                          :: functional

    REAL(real_8), PARAMETER                  :: fpp = 1.709921_real_8

    REAL(real_8) :: dDFP_dnab, df_of_zeta_dna, df_of_zeta_dnb, df_of_zeta_dz, &
      DFP, dipo_alpha_dna, dipo_alpha_dnb, dipo_delta_dna, dipo_delta_dnb, &
      dzeta_4_dna, dzeta_4_dnb, dzeta_dna, dzeta_dnb, f_of_zeta, ipo_alpha, &
      ipo_delta, n_ab, nab_2, rs_ab, scobpw, zeta, zeta_2, zeta_3, zeta_4
    TYPE(cp_xc_functional_t)                 :: obpw_alpha, obpw_F, obpw_P

! \Delta (e_F - e_P)

    n_ab  = scratch(alpha_beta)%n
    rs_ab = tqbp_to_one_third/(scratch(alpha_beta)%n_1_3)

    CALL cp_lda_c_ob_pw_P( n_ab, rs_ab, obpw_P )
    CALL cp_lda_c_ob_pw_F( n_ab, rs_ab, obpw_F )
    CALL cp_lda_c_pw_alpha_fit( n_ab, rs_ab, obpw_alpha )

    zeta            = zeta_of( scratch(:) )
    zeta_2          = zeta   * zeta
    zeta_3          = zeta_2 * zeta
    zeta_4          = zeta_2 * zeta_2

    nab_2           = n_ab*n_ab 

    dzeta_dna       =  2.0_real_8 * scratch(beta)%n  / nab_2
    dzeta_dnb       = -2.0_real_8 * scratch(alpha)%n / nab_2

    dzeta_4_dna     = 4.0_real_8*zeta_3*dzeta_dna
    dzeta_4_dnb     = 4.0_real_8*zeta_3*dzeta_dnb

    f_of_zeta       = f_of( zeta )
    df_of_zeta_dz   = df_dz_of( zeta )
    df_of_zeta_dna  = df_of_zeta_dz * dzeta_dna
    df_of_zeta_dnb  = df_of_zeta_dz * dzeta_dnb

    DFP             = obpw_F%sc  - obpw_P%sc
    !
    ! No explicit a, b dependency
    !
    dDFP_dnab       = obpw_F%v1c - obpw_P%v1c

    !
    ! Alpha takes a negative sign; fix this by changing the sign of ipo_alpha
    !
    ipo_alpha       = - ( f_of_zeta / fpp * (1.0_real_8 - zeta_4) )
    dipo_alpha_dna  = - ( df_of_zeta_dna / fpp * (1.0_real_8 - zeta_4) &
         - f_of_zeta / fpp * (dzeta_4_dna) )
    dipo_alpha_dnb  = - ( df_of_zeta_dnb / fpp * (1.0_real_8 - zeta_4) &
         - f_of_zeta / fpp * (dzeta_4_dnb) )

    ipo_delta       = f_of_zeta * zeta_4
    dipo_delta_dna  = df_of_zeta_dna * zeta_4 + f_of_zeta * dzeta_4_dna
    dipo_delta_dnb  = df_of_zeta_dnb * zeta_4 + f_of_zeta * dzeta_4_dnb

    scobpw          = obpw_P%sc + ipo_delta*DFP + ipo_alpha*obpw_alpha%sc

    functional(alpha)%sc       = 0.0_real_8      
    functional(beta)%sc        = 0.0_real_8      
    functional(alpha_beta)%sc  = scobpw

    functional(alpha)%v1c      = obpw_P%v1c + dipo_delta_dna*DFP + ipo_delta*dDFP_dnab &
         + dipo_alpha_dna*obpw_alpha%sc  + ipo_alpha*obpw_alpha%v1c
    functional(beta)%v1c       = obpw_P%v1c + dipo_delta_dnb*DFP + ipo_delta*dDFP_dnab &
         + dipo_alpha_dnb*obpw_alpha%sc  + ipo_alpha*obpw_alpha%v1c
    functional(alpha_beta)%v1c = 0.0_real_8

    functional(alpha)%v2c      = 0.0_real_8
    functional(beta)%v2c       = 0.0_real_8 
    functional(alpha_beta)%v2c = 0.0_real_8

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_spin_lda_c_ob_pw
  ! ==================================================================
  PURE SUBROUTINE cp_lda_c_ob_pw(scratch,functional)
    ! ==--------------------------------------------------------------==
    ! Rewritten from scratch
    !                                            M.P.Bircher @ LCBC/EPFL
    ! ==--------------------------------------------------------------==


    TYPE(cp_xc_scratch_t), INTENT(in)        :: scratch
    TYPE(cp_xc_functional_t), INTENT(inout)  :: functional

    REAL(real_8)                             :: n_ab, rs_ab

    n_ab  = 2.0_real_8*scratch%n
    rs_ab = tqbp_to_one_third/(two_to_one_third*scratch%n_1_3)

    CALL cp_lda_c_ob_pw_P( n_ab, rs_ab, functional )

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_lda_c_ob_pw
  ! ==================================================================
  PURE SUBROUTINE cp_spin_lda_c_pz(scratch,functional)
    ! ==--------------------------------------------------------------==
    ! Ported from OLDCODE
    !                                            M.P.Bircher @ LCBC/EPFL
    ! ==--------------------------------------------------------------==


    TYPE(cp_xc_scratch_t), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(in)                             :: scratch
    TYPE(cp_xc_functional_t), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(inout)                          :: functional

    REAL(real_8)                             :: df_of_zeta_dz, f_of_zeta, &
                                                n_ab, rs_ab, scpz, vpza, &
                                                vpzb, zeta
    TYPE(cp_xc_functional_t)                 :: pz_F, pz_P

    n_ab  = scratch(alpha_beta)%n
    rs_ab = tqbp_to_one_third/(scratch(alpha_beta)%n_1_3)

    CALL cp_lda_c_pz_P( n_ab, rs_ab, pz_P )
    CALL cp_lda_c_pz_F( n_ab, rs_ab, pz_F )

    zeta          = zeta_of( scratch(:) )
    f_of_zeta     = f_of( zeta )
    df_of_zeta_dz = df_dz_of( zeta )

    scpz          = pz_P%sc  + f_of_zeta*(pz_F%sc  - pz_P%sc)
    vpza          = pz_P%v1c + f_of_zeta*(pz_F%v1c - pz_P%v1c) &
         + (pz_F%sc - pz_P%sc)/n_ab*(1._real_8-zeta)*df_of_zeta_dz
    vpzb          = pz_P%v1c + f_of_zeta*(pz_F%v1c - pz_P%v1c) &
         - (pz_F%sc - pz_P%sc)/n_ab*(1._real_8+zeta)*df_of_zeta_dz

    functional(alpha)%sc       = 0.0_real_8      
    functional(beta)%sc        = 0.0_real_8      
    functional(alpha_beta)%sc  = scpz

    functional(alpha)%v1c      = vpza 
    functional(beta)%v1c       = vpzb
    functional(alpha_beta)%v1c = 0.0_real_8

    functional(alpha)%v2c      = 0.0_real_8
    functional(beta)%v2c       = 0.0_real_8 
    functional(alpha_beta)%v2c = 0.0_real_8

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_spin_lda_c_pz
  ! ==================================================================
  PURE SUBROUTINE cp_lda_c_pz(scratch,functional)
    ! ==--------------------------------------------------------------==
    ! Ported from OLDCODE
    !                                            M.P.Bircher @ LCBC/EPFL
    ! ==--------------------------------------------------------------==


    TYPE(cp_xc_scratch_t), INTENT(in)        :: scratch
    TYPE(cp_xc_functional_t), INTENT(inout)  :: functional

    REAL(real_8)                             :: n_ab, rs_ab

    n_ab  = 2.0_real_8*scratch%n
    rs_ab = tqbp_to_one_third/(two_to_one_third*scratch%n_1_3)

    CALL cp_lda_c_pz_P( n_ab, rs_ab, functional )

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_lda_c_pz
  ! ==================================================================

  ! ==================================================================
  ! Routines for ferromagnetic, paramagnetic state and interpolation
  ! ==================================================================

  ! ==================================================================
  ELEMENTAL SUBROUTINE cp_lda_c_pz_P(n,rs,functional)
    ! ==--------------------------------------------------------------==
    ! Eqns from OLDCODE; split into ferromagnetic and paramagnetic part
    !                                            M.P.Bircher @ LCBC/EPFL
    ! ==--------------------------------------------------------------==


    REAL(real_8), INTENT(in)                 :: n, rs
    TYPE(cp_xc_functional_t), INTENT(inout)  :: functional

    REAL(real_8), PARAMETER :: a = 0.0311_real_8, b = -0.048_real_8, &
      b1 = 1.0529_real_8, b2 = 0.3334_real_8, c = 0.0020_real_8, &
      d = -0.0116_real_8, gc = -0.1423_real_8

    REAL(real_8)                             :: dox, epz, ox, vpz, x, xln

    IF (rs < 1.0_real_8) THEN
       !
       ! High density eq.
       !
       xln = LOG(rs)
       epz = a*xln + b + c*rs*xln + d*rs
       vpz = a*xln + (b-a/3._real_8) + 2._real_8/3._real_8*c*rs*xln + (2._real_8*d-c)/3._real_8*rs
    ELSE
       !
       ! Low density interpolation
       !
       x   = SQRT(rs)
       ox  = 1._real_8 + b1*x + b2*rs
       dox = 1._real_8 + 7._real_8/6._real_8*b1*x + 4._real_8/3._real_8*b2*rs
       epz = gc/ox
       vpz = epz*dox/ox
    ENDIF

    functional%sc  = n*epz
    functional%v1c = vpz
    functional%v2c = 0.0_real_8

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_lda_c_pz_P
  ! ==================================================================
  ELEMENTAL SUBROUTINE cp_lda_c_pz_F(n,rs,functional)
    ! ==--------------------------------------------------------------==
    ! Eqns from OLDCODE; split into ferromagnetic and paramagnetic part
    !                                            M.P.Bircher @ LCBC/EPFL
    ! ==--------------------------------------------------------------==


    REAL(real_8), INTENT(in)                 :: n, rs
    TYPE(cp_xc_functional_t), INTENT(inout)  :: functional

    REAL(real_8), PARAMETER :: a = 0.01555_real_8, b = -0.02690_real_8, &
      b1 = 1.39810_real_8, b2 = 0.26110_real_8, c = 0.00070_real_8, &
      d = -0.00480_real_8, gc = -0.08430_real_8

    REAL(real_8)                             :: dox, epz, ox, vpz, x, xln

    IF (rs < 1.0_real_8) THEN
       !
       ! High density eq.
       !
       xln = LOG(rs)
       epz = a*xln + b + c*rs*xln + d*rs
       vpz = a*xln + (b-a/3._real_8) + 2._real_8/3._real_8*c*rs*xln + (2._real_8*d-c)/3._real_8*rs
    ELSE
       !
       ! Low density interpolation
       !
       x   = SQRT(rs)
       ox  = 1._real_8 + b1*x + b2*rs
       dox = 1._real_8 + 7._real_8/6._real_8*b1*x + 4._real_8/3._real_8*b2*rs
       epz = gc/ox
       vpz = epz*dox/ox
    ENDIF

    functional%sc  = n*epz
    functional%v1c = vpz
    functional%v2c = 0.0_real_8

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_lda_c_pz_F
  ! ==================================================================
  ELEMENTAL SUBROUTINE cp_lda_c_vwn_P(n,rs,functional)
    ! ==--------------------------------------------------------------==
    ! Bosonic limit (p. 1207), eqs based on OLDCODE
    !                                            M.P.Bircher @ LCBC/EPFL
    ! ==--------------------------------------------------------------==


    REAL(real_8), INTENT(in)                 :: n, rs
    TYPE(cp_xc_functional_t), INTENT(inout)  :: functional

    REAL(real_8), PARAMETER :: a = 0.0310907_real_8, b = 3.72744_real_8, &
      c = 12.9352_real_8, q = SQRT(4._real_8*c-b*b), f1 = 2.0_real_8*b/q, &
      x0 = -0.10498_real_8 , f3 = 2.0_real_8*(2.0_real_8*x0 + b)/q, &
      f2 = b*x0/(x0*x0 + b*x0 + c)

    REAL(real_8)                             :: evwn, fx, qx, ttqq, txpb, &
                                                vvwn, x, xx0

    x    = SQRT(rs)

    fx   = rs + b*x + c
    qx   = ATAN(q/(2.0_real_8*x + b))
    xx0  = x-x0
    evwn = a*(LOG(rs/fx) + f1*qx-f2*(LOG(xx0*xx0/fx) + f3*qx))

    txpb = 2.0_real_8*x + b
    ttqq = txpb*txpb + q*q
    vvwn = evwn - x*a/6._real_8*(2.0_real_8/x-txpb/fx-4._real_8*b/ttqq-f2*(2.0_real_8/xx0 &
         -txpb/fx-4._real_8*(2.0_real_8*x0 + b)/ttqq))

    functional%sc  = evwn*n
    functional%v1c = vvwn
    functional%v2c = 0.0_real_8

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_lda_c_vwn_P
  ! ==================================================================
  ELEMENTAL SUBROUTINE cp_lda_c_vwn_F(n,rs,functional)
    ! ==--------------------------------------------------------------==
    ! Fermionic limit (p. 1207), eqs based on OLDCODE
    !                                            M.P.Bircher @ LCBC/EPFL
    ! ==--------------------------------------------------------------==


    REAL(real_8), INTENT(in)                 :: n, rs
    TYPE(cp_xc_functional_t), INTENT(inout)  :: functional

    REAL(real_8), PARAMETER :: a = 0.01554535_real_8, b = 7.06042_real_8, &
      c = 18.0578_real_8, q = SQRT(4._real_8*c-b*b), f1 = 2.0_real_8*b/q, &
      x0 = -0.32500_real_8 , f3 = 2.0_real_8*(2.0_real_8*x0 + b)/q, &
      f2 = b*x0/(x0*x0 + b*x0 + c)

    REAL(real_8)                             :: evwn, fx, qx, ttqq, txpb, &
                                                vvwn, x, xx0

    x    = SQRT(rs)

    fx   = rs + b*x + c
    qx   = ATAN(q/(2.0_real_8*x + b))
    xx0  = x-x0
    evwn = a*(LOG(rs/fx) + f1*qx-f2*(LOG(xx0*xx0/fx) + f3*qx))

    txpb = 2.0_real_8*x + b
    ttqq = txpb*txpb + q*q
    vvwn = evwn - x*a/6._real_8*(2.0_real_8/x-txpb/fx-4._real_8*b/ttqq-f2*(2.0_real_8/xx0 &
         -txpb/fx-4._real_8*(2.0_real_8*x0 + b)/ttqq))

    functional%sc  = evwn*n
    functional%v1c = vvwn
    functional%v2c = 0.0_real_8

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_lda_c_vwn_F
  ! ==================================================================
  ELEMENTAL SUBROUTINE cp_lda_c_vwn_alpha_fit(n,rs,functional)
    ! ==--------------------------------------------------------------==
    ! Fit to \alpha_c (p. 1209) eqs based on OLDCODE
    !                                            M.P.Bircher @ LCBC/EPFL
    ! ==--------------------------------------------------------------==


    REAL(real_8), INTENT(in)                 :: n, rs
    TYPE(cp_xc_functional_t), INTENT(inout)  :: functional

    REAL(real_8), PARAMETER :: a = -1.0_real_8/(6.0_real_8 * pi * pi), &
      b = 1.13107_real_8, c = 13.0045_real_8, q = SQRT(4._real_8*c-b*b), &
      f1 = 2.0_real_8*b/q, x0 = -0.0047584_real_8 , &
      f3 = 2.0_real_8*(2.0_real_8*x0 + b)/q, f2 = b*x0/(x0*x0 + b*x0 + c)

    REAL(real_8)                             :: evwn, fx, qx, ttqq, txpb, &
                                                vvwn, x, xx0

    x    = SQRT(rs)

    fx   = rs + b*x + c
    qx   = ATAN(q/(2.0_real_8*x + b))
    xx0  = x-x0
    evwn = a*(LOG(rs/fx) + f1*qx-f2*(LOG(xx0*xx0/fx) + f3*qx))

    txpb = 2.0_real_8*x + b
    ttqq = txpb*txpb + q*q
    vvwn = evwn - x*a/6._real_8*(2.0_real_8/x-txpb/fx-4._real_8*b/ttqq-f2*(2.0_real_8/xx0 &
         -txpb/fx-4._real_8*(2.0_real_8*x0 + b)/ttqq))

    functional%sc  = evwn*n
    functional%v1c = vvwn
    functional%v2c = 0.0_real_8

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_lda_c_vwn_alpha_fit
  ! ==================================================================
  ELEMENTAL SUBROUTINE cp_lda_c_pw_P(n,rs,functional)
    ! ==--------------------------------------------------------------==
    ! Paramagnetic limit, eqs from OLDCODE
    !                                            M.P.Bircher @ LCBC/EPFL
    ! ==--------------------------------------------------------------==


    REAL(real_8), INTENT(in)                 :: n, rs
    TYPE(cp_xc_functional_t), INTENT(inout)  :: functional

    REAL(real_8), PARAMETER :: a = 0.0310907_real_8, a1 = 0.21370_real_8, &
      b1 = 7.5957_real_8, b2 = 3.5876_real_8, b3 = 1.6382_real_8, &
      b4 = 0.49294_real_8, c0 = a, c1 = 0.046644_real_8, c2 = 0.00664_real_8, &
      c3 = 0.01043_real_8, d0 = 0.4335_real_8, d1 = 1.4408_real_8

    REAL(real_8)                             :: a1rs, b1rs1, b2rs2, b3rs3, &
                                                b4rs4, c0xln, c2rs, d0_rs, &
                                                dom, epwc, olog, om, rs1, &
                                                rs3, rs4, vpwc, xln

    IF (rs < 0.5_real_8) THEN
       !
       ! High density formula
       !
       xln   = LOG(rs)
       c0xln = c0*xln
       c2rs  = c2*rs
       epwc  = c0xln-c1 + c2rs*xln-c3*rs
       vpwc  = c0xln-(c1 + c0/3._real_8) + 2._real_8/3._real_8*c2rs*xln &
            -(2._real_8*c3 + c2)/3._real_8*rs
    ELSE IF (rs > 100._real_8) THEN
       !
       ! Low density formula
       !
       d0_rs = d0/rs
       epwc  = -d0_rs + d1/rs**1.5_real_8
       vpwc  = -4._real_8/3._real_8*d0_rs + 1.5_real_8*d1/(rs*SQRT(rs))
    ELSE
       !
       ! Interpolation formula
       !
       rs1   = SQRT(rs)
       rs3   = rs*rs1
       rs4   = rs*rs
       b1rs1 = b1*rs1
       b2rs2 = b2*rs
       b3rs3 = b3*rs3
       b4rs4 = b4*rs4
       a1rs  = a1*rs
       om    = 2._real_8*a*(b1rs1 + b2rs2 + b3rs3 + b4rs4)
       dom   = 2._real_8*a*(0.5_real_8*b1rs1 + b2rs2 + 1.5_real_8*b3rs3 + 2._real_8*b4rs4)
       olog  = LOG(1._real_8 + 1.0_real_8/om)
       epwc  = -2._real_8*a*(1.0_real_8 + a1rs)*olog
       vpwc  = -2._real_8*a*(1._real_8 + 2._real_8/3._real_8*a1rs)*olog &
            -2._real_8/3._real_8*a*(1._real_8 + a1rs)*dom/(om*(om + 1._real_8))
    END IF

    functional%sc  = epwc*n
    functional%v1c = vpwc
    functional%v2c = 0.0_real_8

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_lda_c_pw_P
  ! ==================================================================
  ELEMENTAL SUBROUTINE cp_lda_c_pw_F(n,rs,functional)
    ! ==--------------------------------------------------------------==
    ! Ferromagnetic limit, eqs from OLDCODE
    !                                            M.P.Bircher @ LCBC/EPFL
    ! ==--------------------------------------------------------------==


    REAL(real_8), INTENT(in)                 :: n, rs
    TYPE(cp_xc_functional_t), INTENT(inout)  :: functional

    REAL(real_8), PARAMETER :: a = 0.015545_real_8, a1 = 0.20548_real_8, &
      b1 = 14.1189_real_8, b2 = 6.1977_real_8, b3 = 3.3662_real_8, &
      b4 = 0.62517_real_8, c0 = a, c1 = 0.025599_real_8, c2 = 0.00319_real_8, &
      c3 = 0.00384_real_8, d0 = 0.3287_real_8, d1 = 1.7697_real_8

    REAL(real_8)                             :: a1rs, b1rs1, b2rs2, b3rs3, &
                                                b4rs4, c0xln, c2rs, d0_rs, &
                                                dom, epwc, olog, om, rs1, &
                                                rs3, rs4, vpwc, xln

    IF (rs < 0.5_real_8) THEN
       !
       ! High density formula
       !
       xln   = LOG(rs)
       c0xln = c0*xln
       c2rs  = c2*rs
       epwc  = c0xln-c1 + c2rs*xln-c3*rs
       vpwc  = c0xln-(c1 + c0/3._real_8) + 2._real_8/3._real_8*c2rs*xln &
            -(2._real_8*c3 + c2)/3._real_8*rs
    ELSE IF (rs > 100._real_8) THEN
       !
       ! Low density formula
       !
       d0_rs = d0/rs
       epwc  = -d0_rs + d1/rs**1.5_real_8
       vpwc  = -4._real_8/3._real_8*d0_rs + 1.5_real_8*d1/(rs*SQRT(rs))
    ELSE
       !
       ! Interpolation formula
       !
       rs1   = SQRT(rs)
       rs3   = rs*rs1
       rs4   = rs*rs
       b1rs1 = b1*rs1
       b2rs2 = b2*rs
       b3rs3 = b3*rs3
       b4rs4 = b4*rs4
       a1rs  = a1*rs
       om    = 2._real_8*a*(b1rs1 + b2rs2 + b3rs3 + b4rs4)
       dom   = 2._real_8*a*(0.5_real_8*b1rs1 + b2rs2 + 1.5_real_8*b3rs3 + 2._real_8*b4rs4)
       olog  = LOG(1._real_8 + 1.0_real_8/om)
       epwc  = -2._real_8*a*(1.0_real_8 + a1rs)*olog
       vpwc  = -2._real_8*a*(1._real_8 + 2._real_8/3._real_8*a1rs)*olog &
            -2._real_8/3._real_8*a*(1._real_8 + a1rs)*dom/(om*(om + 1._real_8))
    END IF

    functional%sc  = epwc*n
    functional%v1c = vpwc
    functional%v2c = 0.0_real_8

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_lda_c_pw_F
  ! ==================================================================
  ELEMENTAL SUBROUTINE cp_lda_c_pw_alpha_fit(n,rs,functional)
    ! ==--------------------------------------------------------------==
    ! Fit to \alpha_c, eqs from OLDCODE
    !                                            M.P.Bircher @ LCBC/EPFL
    ! ==--------------------------------------------------------------==


    REAL(real_8), INTENT(in)                 :: n, rs
    TYPE(cp_xc_functional_t), INTENT(inout)  :: functional

    REAL(real_8), PARAMETER :: a = 0.016887_real_8, a1 = 0.11125_real_8, &
      b1 = 10.357_real_8, b2 = 3.6231_real_8, b3 = 0.88026_real_8, &
      b4 = 0.49671_real_8, c0 = a, c1 = 0.035475_real_8, c2 = 0.00188_real_8, &
      c3 = 0.00521_real_8, d0 = 0.2240_real_8, d1 = 0.3969_real_8

    REAL(real_8)                             :: a1rs, b1rs1, b2rs2, b3rs3, &
                                                b4rs4, c0xln, c2rs, d0_rs, &
                                                dom, epwc, olog, om, rs1, &
                                                rs3, rs4, vpwc, xln

    IF (rs < 0.5_real_8) THEN
       !
       ! High density formula
       !
       xln   = LOG(rs)
       c0xln = c0*xln
       c2rs  = c2*rs
       epwc  = c0xln-c1 + c2rs*xln-c3*rs
       vpwc  = c0xln-(c1 + c0/3._real_8) + 2._real_8/3._real_8*c2rs*xln &
            -(2._real_8*c3 + c2)/3._real_8*rs
    ELSE IF (rs > 100._real_8) THEN
       !
       ! Low density formula
       !
       d0_rs = d0/rs
       epwc  = -d0_rs + d1/rs**1.5_real_8
       vpwc  = -4._real_8/3._real_8*d0_rs + 1.5_real_8*d1/(rs*SQRT(rs))
    ELSE
       !
       ! Interpolation formula
       !
       rs1   = SQRT(rs)
       rs3   = rs*rs1
       rs4   = rs*rs
       b1rs1 = b1*rs1
       b2rs2 = b2*rs
       b3rs3 = b3*rs3
       b4rs4 = b4*rs4
       a1rs  = a1*rs
       om    = 2._real_8*a*(b1rs1 + b2rs2 + b3rs3 + b4rs4)
       dom   = 2._real_8*a*(0.5_real_8*b1rs1 + b2rs2 + 1.5_real_8*b3rs3 + 2._real_8*b4rs4)
       olog  = LOG(1._real_8 + 1.0_real_8/om)
       epwc  = -2._real_8*a*(1.0_real_8 + a1rs)*olog
       vpwc  = -2._real_8*a*(1._real_8 + 2._real_8/3._real_8*a1rs)*olog &
            -2._real_8/3._real_8*a*(1._real_8 + a1rs)*dom/(om*(om + 1._real_8))
    END IF

    functional%sc  = epwc*n
    functional%v1c = vpwc
    functional%v2c = 0.0_real_8

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_lda_c_pw_alpha_fit
  ! ==================================================================
  ELEMENTAL SUBROUTINE cp_lda_c_ob_pw_P(n,rs,functional)
    ! ==--------------------------------------------------------------==
    ! Paramagnetic limit, eqs from OLDCODE
    !                                            M.P.Bircher @ LCBC/EPFL
    ! ==--------------------------------------------------------------==


    REAL(real_8), INTENT(in)                 :: n, rs
    TYPE(cp_xc_functional_t), INTENT(inout)  :: functional

    REAL(real_8), PARAMETER :: a = 0.031091_real_8, a1 = 0.026481_real_8, &
      b1 = 7.5957_real_8    , b2 = 3.5876_real_8, b3 = -0.46647_real_8, &
      b4 = 0.13354_real_8, c0 = a, c1 = 0.046644_real_8, c2 = 0.00664_real_8, &
      c3 = 0.01043_real_8, d0 = 0.4335_real_8, d1 = 1.4408_real_8 

    REAL(real_8)                             :: a1rs, b1rs1, b2rs2, b3rs3, &
                                                b4rs4, c0xln, c2rs, d0_rs, &
                                                dom, epwc, olog, om, rs1, &
                                                rs3, rs4, vpwc, xln

    IF (rs < 0.5_real_8) THEN
       !
       ! High density formula
       !
       xln   = LOG(rs)
       c0xln = c0*xln
       c2rs  = c2*rs
       epwc  = c0xln-c1 + c2rs*xln-c3*rs
       vpwc  = c0xln-(c1 + c0/3._real_8) + 2._real_8/3._real_8*c2rs*xln &
            -(2._real_8*c3 + c2)/3._real_8*rs
    ELSE IF (rs > 100._real_8) THEN
       !
       ! Low density formula
       !
       d0_rs = d0/rs
       epwc  = -d0_rs + d1/rs**1.5_real_8
       vpwc  = -4._real_8/3._real_8*d0_rs + 1.5_real_8*d1/(rs*SQRT(rs))
    ELSE
       !
       ! Interpolation formula
       !
       rs1   = SQRT(rs)
       rs3   = rs*rs1
       rs4   = rs*rs
       b1rs1 = b1*rs1
       b2rs2 = b2*rs
       b3rs3 = b3*rs3
       b4rs4 = b4*rs4
       a1rs  = a1*rs
       om    = 2._real_8*a*(b1rs1 + b2rs2 + b3rs3 + b4rs4)
       dom   = 2._real_8*a*(0.5_real_8*b1rs1 + b2rs2 + 1.5_real_8*b3rs3 + 2._real_8*b4rs4)
       olog  = LOG(1._real_8 + 1.0_real_8/om)
       epwc  = -2._real_8*a*(1.0_real_8 + a1rs)*olog
       vpwc  = -2._real_8*a*(1._real_8 + 2._real_8/3._real_8*a1rs)*olog &
            -2._real_8/3._real_8*a*(1._real_8 + a1rs)*dom/(om*(om + 1._real_8))
    END IF

    functional%sc  = epwc*n
    functional%v1c = vpwc
    functional%v2c = 0.0_real_8

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_lda_c_ob_pw_P
  ! ==================================================================
  ELEMENTAL SUBROUTINE cp_lda_c_ob_pw_F(n,rs,functional)
    ! ==--------------------------------------------------------------==
    ! Ferromagnetic limit, eqs from OLDCODE
    !                                            M.P.Bircher @ LCBC/EPFL
    ! ==--------------------------------------------------------------==


    REAL(real_8), INTENT(in)                 :: n, rs
    TYPE(cp_xc_functional_t), INTENT(inout)  :: functional

    REAL(real_8), PARAMETER :: a = 0.015545_real_8, a1 = 0.022465_real_8, &
      b1 = 14.1189_real_8, b2 = 6.1977_real_8, b3 = -0.56043_real_8, &
      b4 = 0.11313_real_8, c0 = a, c1 = 0.025599_real_8, c2 = 0.00319_real_8, &
      c3 = 0.00384_real_8, d0 = 0.3287_real_8, d1 = 1.7697_real_8

    REAL(real_8)                             :: a1rs, b1rs1, b2rs2, b3rs3, &
                                                b4rs4, c0xln, c2rs, d0_rs, &
                                                dom, epwc, olog, om, rs1, &
                                                rs3, rs4, vpwc, xln

    IF (rs < 0.5_real_8) THEN
       !
       ! High density formula
       !
       xln   = LOG(rs)
       c0xln = c0*xln
       c2rs  = c2*rs
       epwc  = c0xln-c1 + c2rs*xln-c3*rs
       vpwc  = c0xln-(c1 + c0/3._real_8) + 2._real_8/3._real_8*c2rs*xln &
            -(2._real_8*c3 + c2)/3._real_8*rs
    ELSE IF (rs > 100._real_8) THEN
       !
       ! Low density formula
       !
       d0_rs = d0/rs
       epwc  = -d0_rs + d1/rs**1.5_real_8
       vpwc  = -4._real_8/3._real_8*d0_rs + 1.5_real_8*d1/(rs*SQRT(rs))
    ELSE
       !
       ! Interpolation formula
       !
       rs1   = SQRT(rs)
       rs3   = rs*rs1
       rs4   = rs*rs
       b1rs1 = b1*rs1
       b2rs2 = b2*rs
       b3rs3 = b3*rs3
       b4rs4 = b4*rs4
       a1rs  = a1*rs
       om    = 2._real_8*a*(b1rs1 + b2rs2 + b3rs3 + b4rs4)
       dom   = 2._real_8*a*(0.5_real_8*b1rs1 + b2rs2 + 1.5_real_8*b3rs3 + 2._real_8*b4rs4)
       olog  = LOG(1._real_8 + 1.0_real_8/om)
       epwc  = -2._real_8*a*(1.0_real_8 + a1rs)*olog
       vpwc  = -2._real_8*a*(1._real_8 + 2._real_8/3._real_8*a1rs)*olog &
            -2._real_8/3._real_8*a*(1._real_8 + a1rs)*dom/(om*(om + 1._real_8))
    END IF

    functional%sc  = epwc*n
    functional%v1c = vpwc
    functional%v2c = 0.0_real_8

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_lda_c_ob_pw_F
  ! ==================================================================

  ! ==================================================================
  ! Recurring functions
  ! ==================================================================

  ! ==================================================================
  PURE FUNCTION zeta_of(scratch) &
       RESULT  (zeta)
    ! ==--------------------------------------------------------------==
    ! Calculates zeta(n_a,n_b)
    !                                            M.P.Bircher @ LCBC/EPFL
    ! ==--------------------------------------------------------------==

    TYPE(cp_xc_scratch_t), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(in)                             :: scratch
    REAL(real_8)                             :: zeta

    zeta = (scratch(alpha)%n - scratch(beta)%n) / scratch(alpha_beta)%n

    ! ==--------------------------------------------------------------==
  END FUNCTION zeta_of
  ! ==================================================================
  ELEMENTAL FUNCTION f_of(zeta) &
       RESULT  (f_of_zeta)
    ! ==--------------------------------------------------------------==
    ! Calculates f(zeta)
    !                                            M.P.Bircher @ LCBC/EPFL
    ! ==--------------------------------------------------------------==

    REAL(real_8), INTENT(in)                 :: zeta
    REAL(real_8)                             :: f_of_zeta

    REAL(real_8), PARAMETER :: &
      f_denominator = 2._real_8*(two_to_one_third - 1._real_8)

    f_of_zeta = ( (1.0_real_8 + zeta)**four_thirds &
         + (1.0_real_8 - zeta)**four_thirds &
         - 2.0_real_8 ) / f_denominator 

    ! ==--------------------------------------------------------------==
  END FUNCTION f_of
  ! ==================================================================
  ELEMENTAL FUNCTION df_dz_of(zeta) &
       RESULT  (df_dz_of_zeta)
    ! ==--------------------------------------------------------------==
    ! Derivatives of f(zeta)
    !                                            M.P.Bircher @ LCBC/EPFL
    ! ==--------------------------------------------------------------==

    REAL(real_8), INTENT(in)                 :: zeta
    REAL(real_8)                             :: df_dz_of_zeta

    REAL(real_8), PARAMETER :: &
      f_denominator = 2._real_8*(two_to_one_third - 1._real_8)

    df_dz_of_zeta = ( four_thirds*(1.0_real_8 + zeta)**one_third   &
         - four_thirds*(1.0_real_8 - zeta)**one_third ) &
         / f_denominator

    ! ==--------------------------------------------------------------==
  END FUNCTION df_dz_of
  ! ==================================================================

  ! ==================================================================
  ! Check parameter compatibility
  ! ==================================================================
 
  ! ==================================================================
  ELEMENTAL FUNCTION cp_lda_c_check(tag) &
  RESULT (OK)
    ! ==--------------------------------------------------------------==
    ! Does not really do anything for a pure LDA functional, but is
    ! important for flexible GGA/MGGA routines
    !
    !                                 25.07.2017 M.P.Bircher @ LCBC/EPFL
    ! ==--------------------------------------------------------------==
 
    CHARACTER(len=*), INTENT(in)             :: tag
    LOGICAL                                  :: OK 
 
    OK = .false.
 
    SELECT CASE( tag )
    !
    ! Functionals with flexible parameters
    !
    CASE( "CP_LDA_C_VWN","CP_SPIN_LDA_C_VWN",&
          "CP_LDA_C_PZ", "CP_SPIN_LDA_C_PZ",&
          "CP_LDA_C_PW", "CP_SPIN_LDA_C_PW",&
          "CP_LDA_C_OB_PW", "CP_SPIN_LDA_C_OB_PW" )
      OK = .true.
    !
    ! Any other tag
    !
    CASE DEFAULT
      OK = .false.
    END SELECT
 
    ! ==--------------------------------------------------------------==
  END FUNCTION cp_lda_c_check
  ! ==================================================================
END MODULE cp_lda_correlation_utils
