! ==================================================================
! Provides: - GGA correlation, adapted from OLDCODE
!
! Input:  Spin-polarised density (1/2 n for closed shell systems)
! Output: Converted to CLOSED-SHELL system for the non-SPIN routines
!         Alpha, beta and alpha+beta for the SPIN routines
!
!                             02.10.2016 - M. P. Bircher @ LCBC/EPFL
! ==================================================================
MODULE cp_gga_correlation_utils
  USE cnst,                            ONLY: pi
  USE cp_lda_correlation_utils,        ONLY: &
       cp_lda_c_check, cp_lda_c_ob_pw, cp_lda_c_pw, cp_lda_c_pw_F, cp_lda_c_pw_P, &
       cp_lda_c_pw_alpha_fit, cp_lda_c_pz, cp_lda_c_vwn, cp_spin_lda_c_ob_pw, &
       cp_spin_lda_c_pw, cp_spin_lda_c_pz, cp_spin_lda_c_vwn, df_dz_of, f_of, &
       zeta_of
  USE cpfunc_types,                    ONLY: cp_xc_functional_t,&
                                             cp_xc_scratch_t,&
                                             cp_xc_spin_components,&
                                             cp_xc_spin_pairs
  USE kinds,                           ONLY: default_string_length,&
                                             real_8

  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER      :: alpha                  = 1
  INTEGER, PARAMETER      :: beta                   = 2
  INTEGER, PARAMETER      :: alpha_beta             = 3

  INTEGER, PARAMETER      :: hcth_gm                = 4

  REAL(real_8), PARAMETER :: three_quarters_by_pi   =  0.75_real_8/pi
  REAL(real_8), PARAMETER :: one_third              =  1.0_real_8/3.0_real_8
  REAL(real_8), PARAMETER :: two_thirds             =  2.0_real_8/3.0_real_8
  REAL(real_8), PARAMETER :: four_thirds            =  4.0_real_8/3.0_real_8
  REAL(real_8), PARAMETER :: seven_thirds           =  7.0_real_8/3.0_real_8
  REAL(real_8), PARAMETER :: eight_thirds           =  8.0_real_8/3.0_real_8
  REAL(real_8), PARAMETER :: eleven_thirds          = 11.0_real_8/3.0_real_8
  REAL(real_8), PARAMETER :: tqbp_to_one_third      = three_quarters_by_pi**(one_third)
  REAL(real_8), PARAMETER :: nine_quarters_times_pi =  2.25_real_8*pi
  REAL(real_8), PARAMETER :: nqtp_to_one_third      = nine_quarters_times_pi**(one_third)
  REAL(real_8), PARAMETER :: three_pi_to_one_third  = (3.0_real_8*pi)**(one_third)
  REAL(real_8), PARAMETER :: two_to_one_third       =  2.0_real_8**(one_third)
  REAL(real_8), PARAMETER :: oo_two_to_one_third    =  1.0_real_8/two_to_one_third
  REAL(real_8), PARAMETER :: two_to_four_thirds     =  2.0_real_8**(4.0_real_8/3.0_real_8)


  !
  ! Customisable functionals (FLEX type)
  !
  TYPE, PRIVATE :: hcth_c_t 
    REAL(real_8)                       :: gamma_parallel   =  0.200000_real_8
    REAL(real_8)                       :: gamma_opposite   =  0.006000_real_8
    REAL(real_8)                       :: c_0_parallel =  0.222601_real_8
    REAL(real_8)                       :: c_0_opposite =  0.729974_real_8
    REAL(real_8), DIMENSION(hcth_gm)   :: c_parallel   =  (/ -0.338622e-01_real_8, &
                                                             -0.125170e-01_real_8, &
                                                             -0.802496e+00_real_8, &
                                                              0.155396e+01_real_8 /)
    REAL(real_8), DIMENSION(hcth_gm)   :: c_opposite   =  (/  0.335287e+01_real_8, &
                                                             -0.115430e+02_real_8, &
                                                              0.808564e+01_real_8, &
                                                             -0.447857e+01_real_8 /)
  END TYPE hcth_c_t
  !
  TYPE, PRIVATE :: pbe_c_t
     REAL(real_8)                     :: pbe_beta  = 0.06672455060314922_real_8
     REAL(real_8)                     :: pbe_gamma = 0.031090690869654895_real_8 
     CHARACTER(default_string_length) :: lda_c = 'CP_LDA_C_PW'   
  END TYPE pbe_c_t
  !
  TYPE, PRIVATE :: cp_gga_c_param_t
     LOGICAL                          :: init = .false.
     TYPE(hcth_c_t)                   :: hcth
     TYPE(pbe_c_t)                    :: pbe
  END TYPE cp_gga_c_param_t
  !
  TYPE(cp_gga_c_param_t), PUBLIC, SAVE :: cp_gga_c_param

  PUBLIC :: CP_GGA_C_HCTH
  PUBLIC :: CP_GGA_C_LYP
  PUBLIC :: CP_GGA_C_P86
  PUBLIC :: CP_GGA_C_PBE
  PUBLIC :: CP_GGA_C_PBE_SOL
  PUBLIC :: CP_GGA_C_PBE_FLEX
  PUBLIC :: CP_SPIN_GGA_C_HCTH
  PUBLIC :: CP_SPIN_GGA_C_LYP
  PUBLIC :: CP_SPIN_GGA_C_P86
  PUBLIC :: CP_SPIN_GGA_C_PBE
  PUBLIC :: CP_SPIN_GGA_C_PBE_SOL
  PUBLIC :: CP_SPIN_GGA_C_PBE_FLEX

  PUBLIC :: CP_GGA_C_CHECK

  PRIVATE :: phi_of

CONTAINS

  ! ==================================================================
  PURE SUBROUTINE cp_spin_gga_c_lyp(scratch,functional)
    ! ==--------------------------------------------------------------==
    ! Ported and adapted from OLDCODE (combined lsd_lyp and lsd_glyp)
    !                                 07.04.2017 M.P.Bircher @ LCBC/EPFL
    ! ==--------------------------------------------------------------==


    TYPE(cp_xc_scratch_t), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(in)                             :: scratch
    TYPE(cp_xc_functional_t), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(inout)                          :: functional

    REAL(real_8), PARAMETER :: a = 0.04918_real_8, b = 0.132_real_8, &
      ab = a*b, c = 0.2533_real_8, cf = 2.87123400018819108_real_8, &
      d = 0.349_real_8 , small = 1.0e-24_real_8

    REAL(real_8) :: ab_or, cf_a_b, d_n_ab_4_3_dr, dder, de1a, de1b, de2a, &
      de2b, der, dlaa, dlaaa, dlaab, dlab, dlaba, dlabb, dlbb, dlbba, dlbbb, &
      dor, dor_or, dr, e1, e2, n_a, n_a_b, n_a_by_ab, n_ab, n_ab_11_3, &
      n_ab_1_3, n_ab_4_3, n_ab_8_3, n_ab_ab, n_b, n_b_by_ab, or, ra, ra_8_3, &
      ra_rb, rb, rb_8_3, v1ca, v1cb, v2ca, v2cab, v2cb, valyp, vblyp, zeta

    n_a       = scratch(alpha)%n
    n_b       = scratch(beta)%n
    n_ab      = scratch(alpha_beta)%n
    n_ab      = MAX(n_ab,small)

    n_a_b     = n_a*n_b
    n_ab_ab   = n_ab*n_ab

    n_a_by_ab = n_a/n_ab
    n_b_by_ab = n_b/n_ab

    n_ab_1_3  = 1.0_real_8 / scratch(alpha_beta)%n_1_3
    n_ab_4_3  = 1.0_real_8 / scratch(alpha_beta)%n_4_3
    n_ab_8_3  = n_ab_4_3*n_ab_4_3
    n_ab_11_3 = n_ab_8_3/n_ab

    zeta          = zeta_of(scratch)
    ra            = n_ab*0.5_real_8*(1._real_8 + zeta)
    ra            = MAX(ra,small)
    rb            = n_ab*0.5_real_8*(1._real_8 - zeta)
    rb            = MAX(rb,small)
    ra_rb         = ra * rb
    ra_8_3        = ra**(eight_thirds)
    rb_8_3        = rb**(eight_thirds)
    dr            = (1._real_8 + d*n_ab_1_3)
    d_n_ab_4_3_dr = d*n_ab_4_3/dr
    cf_a_b        = cf*a*b

    e1       = 4._real_8*a*ra_rb/(n_ab*dr)
    or       = EXP( - c*n_ab_1_3)/dr*n_ab_11_3
    dor      = -one_third*n_ab_4_3*or*(11._real_8/n_ab_1_3 - c - d/dr)
    e2       = 2._real_8**(eleven_thirds)*cf_a_b*or*ra_rb*(ra_8_3 + rb_8_3)
    de1a     = -e1*(one_third*d_n_ab_4_3_dr + 1._real_8/ra - 1._real_8/n_ab)
    de1b     = -e1*(one_third*d_n_ab_4_3_dr + 1._real_8/rb - 1._real_8/n_ab)
    de2a     = -2._real_8**(eleven_thirds)*cf_a_b*( dor*ra_rb*(ra_8_3 + &
         rb_8_3) + or*rb*(eleven_thirds*ra_8_3 + &
         rb_8_3) )
    de2b     = -2._real_8**(eleven_thirds)*cf_a_b*( dor*ra_rb*(ra_8_3  +  &
         rb_8_3)  +  or*ra*(eleven_thirds*rb_8_3  +  &
         ra_8_3) )
    valyp    = de1a + de2a
    vblyp    = de1b + de2b

    dor_or   = dor/or
    ab_or    = ab*or
    der      = c*n_ab_1_3 + d*n_ab_1_3/dr
    dder     = 1._real_8/3._real_8*(d*d_n_ab_4_3_dr*n_ab_1_3/(dr) - der/n_ab)
    dlaa     = -ab_or*(n_a_b/9._real_8*(1._real_8 - 3.0_real_8*der - (der - 11._real_8)*n_a/n_ab) - n_b*n_b)
    dlab     = -ab_or*(n_a_b/9._real_8*(47._real_8 - 7._real_8*der) - 4._real_8/3._real_8*n_ab_ab)
    dlbb     = -ab_or*(n_a_b/9._real_8*(1._real_8 - 3.0_real_8*der - (der - 11._real_8)*n_b/n_ab) - n_a*n_a)
    dlaaa    = dor_or*dlaa - ab_or*(n_b/9._real_8*(1._real_8 - 3.0_real_8*der - (der - 11._real_8)*n_a_by_ab) - &
         n_a_b/9._real_8*((3._real_8 + n_a_by_ab)*dder + (der - 11._real_8)*n_b/(n_ab_ab)))
    dlaab    = dor_or*dlaa - ab_or*(n_a/9._real_8*(1._real_8 - 3._real_8*der&
         - (der - 11._real_8)*n_a_by_ab) - n_a_b/9._real_8*((3._real_8 + n_a_by_ab)*dder&
         - (der - 11._real_8)*n_a/(n_ab_ab)) - 2._real_8*n_b)
    dlaba    = dor_or*dlab - ab_or*(n_b/9._real_8*(47._real_8 - 7._real_8*der)&
         - 7._real_8/9._real_8*n_a_b*dder - 8._real_8/3._real_8*n_ab)
    dlabb    = dor_or*dlab - ab_or*(n_a/9._real_8*(47._real_8 - 7._real_8*der)&
         - 7._real_8/9._real_8*n_a_b*dder - 8._real_8/3._real_8*n_ab)
    dlbba    = dor_or*dlbb - ab_or*(n_b/9._real_8*(1._real_8 - 3._real_8*der&
         - (der - 11._real_8)*n_b_by_ab) - n_a_b/9._real_8*((3._real_8 + n_b_by_ab)*dder&
         - (der - 11._real_8)*n_b/(n_ab_ab)) - 2._real_8*n_a)
    dlbbb    = dor_or*dlbb - ab_or*(n_a/9._real_8*(1._real_8 - 3.0_real_8*der&
         - (der - 11._real_8)*n_b_by_ab) - n_a_b/9._real_8*((3._real_8 + n_b_by_ab)*dder&
         + (der - 11._real_8)*n_a/(n_ab_ab)))
    v1ca     = dlaaa*scratch(alpha)%g_2+ dlaba*scratch(alpha_beta)%g_2 + dlbba*scratch(beta)%g_2 
    v1cb     = dlaab*scratch(alpha)%g_2+ dlabb*scratch(alpha_beta)%g_2 + dlbbb*scratch(beta)%g_2 
    v2ca     = 2.0_real_8*dlaa
    v2cb     = 2.0_real_8*dlbb
    v2cab    = dlab

    functional(alpha)%sc       = -e1 + dlaa*scratch(alpha)%g_2
    functional(beta)%sc        = -e2 + dlbb*scratch(beta)%g_2 
    functional(alpha_beta)%sc  = dlab*scratch(alpha_beta)%g_2

    functional(alpha)%v1c      = valyp + v1ca
    functional(beta)%v1c       = vblyp + v1cb
    functional(alpha_beta)%v1c = 0.0_real_8

    functional(alpha)%v2c      = v2ca
    functional(beta)%v2c       = v2cb
    functional(alpha_beta)%v2c = v2cab

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_spin_gga_c_lyp
  ! ==================================================================
  PURE SUBROUTINE cp_gga_c_lyp(scratch,functional)
    ! ==--------------------------------------------------------------==
    ! Ported and adapted from OLDCODE (combined lyp and glyp)
    !                                 14.10.2016 M.P.Bircher @ LCBC/EPFL
    ! ==--------------------------------------------------------------==

    TYPE(cp_xc_scratch_t), INTENT(in)        :: scratch
    TYPE(cp_xc_functional_t), INTENT(inout)  :: functional

    REAL(real_8), PARAMETER :: a = 0.04918_real_8, b = 0.132_real_8, &
      ab = a*b, c = 0.2533_real_8, cf = 2.87123400018819108_real_8, &
      d = 0.349_real_8 

    REAL(real_8) :: cdoo_n_1_3, coo_n_1_3, dom, doo_n_1_3, doo_n_1_3_p_1, &
      dr5, dxl, ecrs, elyp, expcn, ff, grho, n, n_4_3, om, omxl, oo_n_1_3, &
      oo_n_5_3, ox, sc, slyp, v1c, v2c, vlyp, xl

    oo_n_1_3      = oo_two_to_one_third/scratch%n_1_3
    n             = 2.0_real_8*scratch%n
    n_4_3         = two_to_four_thirds*scratch%n_4_3

    doo_n_1_3     = d*oo_n_1_3
    coo_n_1_3     = c*oo_n_1_3
    cdoo_n_1_3    = c*doo_n_1_3
    doo_n_1_3_p_1 = doo_n_1_3 + 1.0_real_8
    expcn         = EXP(-coo_n_1_3)

    ecrs     = b*cf*expcn
    ox       = 1._real_8/(doo_n_1_3_p_1)
    elyp     = -a*ox*(1._real_8 + ecrs)
    slyp     = n*elyp
    vlyp     = elyp - oo_n_1_3/3._real_8*a*ox*(d*ox + ecrs*(d*ox + c))

    grho     = 4.0_real_8*scratch%g_2
    om       = expcn*ox
    oo_n_5_3 = oo_n_1_3/n_4_3
    !
    xl       = 1._real_8 + seven_thirds*(coo_n_1_3  +  doo_n_1_3*ox)
    omxl     = om*xl       
    ff       = ab*grho/24._real_8
    sc       = ff*oo_n_5_3*omxl
    !
    dr5      = 5._real_8/n_4_3
    dom      = -om*(c + d + cdoo_n_1_3)/(doo_n_1_3_p_1)
    dxl      = seven_thirds*(c + d + 2._real_8*cdoo_n_1_3 + cdoo_n_1_3*doo_n_1_3)/(doo_n_1_3_p_1)**2
    v1c      = -ff/(3._real_8*n_4_3)*( dr5*omxl  +  oo_n_5_3*(dom*xl + om*dxl))
    v2c      = ab*oo_n_5_3*omxl/12._real_8

    functional%sc  = slyp + sc
    functional%v1c = vlyp + v1c
    functional%v2c = v2c

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_gga_c_lyp
  ! ==================================================================
  PURE SUBROUTINE cp_spin_gga_c_p86(scratch,functional)
    ! ==--------------------------------------------------------------==
    ! Ported and adapted from OLDCODE (combined lsd_p86 and lsd_pz)
    !                                 24.05.2017 M.P.Bircher @ LCBC/EPFL
    ! ==--------------------------------------------------------------==


    TYPE(cp_xc_scratch_t), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(in)                             :: scratch
    TYPE(cp_xc_functional_t), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(inout)                          :: functional

    REAL(real_8), PARAMETER :: p1 = 0.023266_real_8, p2 = 7.389e-6_real_8, &
      p3 = 8.723_real_8, p4 = 0.472_real_8 , pc1 = 0.001667_real_8, &
      pc2 = 0.002568_real_8, pci = pc1 + pc2, cdrs = -one_third*(three_quarters_by_pi)**one_third, &
      seven_sixths = 7.0_real_8 / 6.0_real_8

    REAL(real_8)                             :: abs_g, g_2, na_5_3, nb_5_3, nab_5_6, &
                                                br1, br2, br4, cn, &
                                                cna, cnb, d, dcn, dcna, dcnb, &
                                                dda, ddb, drs, ephi, phi, &
                                                rho, rs, rs2, rs3, s1, s2, &
                                                sc_p86, v1c_p86, v1ca_p86, v1cb_p86, v2c_p86


    CALL cp_spin_lda_c_pz( scratch(:), functional(:) )
    
    rs    = tqbp_to_one_third/(scratch(alpha_beta)%n_1_3)
    g_2   = scratch(alpha)%g_2 + 2.0_real_8*scratch(alpha_beta)%g_2 + scratch(beta)%g_2
    abs_g = sqrt(g_2)

    na_5_3  = scratch(alpha)%n_1_3 * scratch(alpha)%n_4_3
    nb_5_3  = scratch(beta)%n_1_3  * scratch(beta)%n_4_3
    nab_5_6 = sqrt(scratch(alpha_beta)%n_1_3 * scratch(alpha_beta)%n_4_3)

    br1   = scratch(alpha_beta)%n_1_3
    br2   = br1*br1
    br4   = br2*br2
    rs2   = rs*rs
    rs3   = rs*rs2
    cna   = pc2 + p1*rs + p2*rs2
    cnb   = 1.0_real_8 + p3*rs + p4*rs2 + 1.e4_real_8*p2*rs3
    cn    = pc1  +  cna/cnb
    drs   = cdrs / br4
    dcna  = (p1 + 2.0_real_8*p2*rs)*drs
    dcnb  = (p3 + 2.0_real_8*p4*rs + 3.e4_real_8*p2*rs2)*drs
    dcn   = dcna/cnb - cna/(cnb*cnb)*dcnb
    s1    = SQRT(na_5_3 + nb_5_3)
    s2    = two_to_one_third/nab_5_6
    d     = s1*s2
    dda   = s2*5.0_real_8/6.0_real_8*(scratch(alpha)%n_1_3 * scratch(alpha)%n_1_3/s1 - s1/scratch(alpha_beta)%n)
    ddb   = s2*5.0_real_8/6.0_real_8*(scratch(beta)%n_1_3  * scratch(beta)%n_1_3/s1  - s1/scratch(alpha_beta)%n)
    phi   = 0.192_real_8*pci/cn*abs_g/(scratch(alpha_beta)%n_1_3*nab_5_6)
    ephi  = EXP(-phi)

    sc_p86    = g_2/br4*cn*ephi/d
    v1c_p86   = sc_p86*((1.0_real_8 + phi)*dcn/cn  - (four_thirds - seven_sixths*phi)/scratch(alpha_beta)%n)
    v1ca_p86  =  - sc_p86/d*dda + v1c_p86
    v1cb_p86  =  - sc_p86/d*ddb + v1c_p86
    v2c_p86   = cn*ephi/br4*(2.0_real_8 - phi)/d

    functional(alpha)%sc       = 0.0_real_8
    functional(beta)%sc        = 0.0_real_8
    functional(alpha_beta)%sc  = functional(alpha_beta)%sc + sc_p86

    functional(alpha)%v1c      = functional(alpha)%v1c + v1ca_p86
    functional(beta)%v1c       = functional(beta)%v1c  + v1cb_p86
    functional(alpha_beta)%v1c = 0.0_real_8

    functional(alpha)%v2c      = v2c_p86 
    functional(beta)%v2c       = v2c_p86
    functional(alpha_beta)%v2c = v2c_p86

    ! ==--------------------------------------------------------------==
    END SUBROUTINE cp_spin_gga_c_p86
  ! ==================================================================
  PURE SUBROUTINE cp_gga_c_p86(scratch,functional)
    ! ==--------------------------------------------------------------==
    ! Ported and adapted from OLDCODE (combined p86 and cp_lda_c_pz)
    !                                 24.05.2017 M.P.Bircher @ LCBC/EPFL
    ! ==--------------------------------------------------------------==

    TYPE(cp_xc_scratch_t), INTENT(in)        :: scratch
    TYPE(cp_xc_functional_t), INTENT(inout)  :: functional

    !
    ! P86
    REAL(real_8), PARAMETER :: p1 = 0.023266_real_8, p2 = 7.389e-6_real_8, &
      p3 = 8.723_real_8, p4 = 0.472_real_8 , pc1 = 0.001667_real_8, &
      pc2 = 0.002568_real_8, cdrs = -one_third*(three_quarters_by_pi)**one_third, &
      pci = pc1 + pc2, seven_sixths = 7.0_real_8 / 6.0_real_8
    !
    ! PZ
    REAL(real_8), PARAMETER :: a = 0.0311_real_8, b = -0.048_real_8, &
      b1 = 1.0529_real_8, b2 = 0.3334_real_8, c = 0.0020_real_8, &
      d = -0.0116_real_8, gc = -0.1423_real_8

    REAL(real_8)                             :: g_2, abs_g, rho, br1, br2, br4, cn, &
                                                cna, cnb, dcn, dcna, dcnb, &
                                                drs, ephi, phi, rs, rs2, rs3
    REAL(real_8)                             :: dox, epz, ox, vpz, x, xln, sc_p86, v1_p86

    g_2   = 4.0_real_8*scratch%g_2
    abs_g = 2.0_real_8*scratch%abs_g
    rho   = 2.0_real_8*scratch%n
    rs    = tqbp_to_one_third/(two_to_one_third*scratch%n_1_3)

    IF (rs < 1.0_real_8) THEN
       !
       ! High density eq.
       !
       xln = LOG(rs)
       epz = a*xln + b + c*rs*xln + d*rs
       vpz = a*xln + (b-a/3._real_8) + two_thirds*c*rs*xln + (2._real_8*d-c)/3._real_8*rs
    ELSE
       !
       ! Low density interpolation
       !
       x   = SQRT(rs)
       ox  = 1._real_8 + b1*x + b2*rs
       dox = 1._real_8 + seven_sixths*b1*x + four_thirds*b2*rs
       epz = gc/ox
       vpz = epz*dox/ox
    ENDIF

    br1   = two_to_one_third * scratch%n_1_3
    br2   = br1*br1
    br4   = br2*br2
    rs2   = rs*rs
    rs3   = rs*rs2
    cna   = pc2 + p1*rs + p2*rs2
    cnb   = 1.0_real_8 + p3*rs + p4*rs2 + 1.e4_real_8*p2*rs3
    cn    = pc1  +  cna/cnb
    drs   = cdrs / br4
    dcna  = (p1 + 2.0_real_8*p2*rs)*drs
    dcnb  = (p3 + 2.0_real_8*p4*rs + 3.e4_real_8*p2*rs2)*drs
    dcn   = dcna/cnb  -  cna/(cnb*cnb)*dcnb
    phi   = 0.192_real_8*pci/cn*abs_g*rho**(-seven_sixths)
    ephi  = EXP( -phi)
    sc_p86 = g_2/br4*cn*ephi 
    v1_p86 = sc_p86*((1.0_real_8 + phi)*dcn/cn  - (four_thirds - (seven_sixths)*phi)/rho)

    functional%sc  = rho*epz + sc_p86
    functional%v1c = vpz + v1_p86
    functional%v2c = cn*ephi/br4*(2.0_real_8 - phi)

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_gga_c_p86
  ! ==================================================================
  PURE SUBROUTINE cp_spin_gga_c_pbe(scratch,functional)
    ! ==--------------------------------------------------------------==
    ! PBE correlation using PW as in the original PBE routine.
    !
    !                                 11.04.2017 M.P.Bircher @ LCBC/EPFL
    ! ==--------------------------------------------------------------==


    TYPE(cp_xc_scratch_t), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(in)                             :: scratch
    TYPE(cp_xc_functional_t), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(inout)                          :: functional

    REAL(real_8), PARAMETER :: fpp = 1.709921_real_8, &
      pbe_beta = 0.06672455060314922_real_8, &
      pbe_gamma = 0.031090690869654895_real_8 , &
      beta_by_gamma = pbe_beta/pbe_gamma , small = 1.0e-24_real_8, &
      xk = (9._real_8*pi/4._real_8)**one_third

    REAL(real_8) :: a, aa, af, da_dra, da_drb, dDFP_dnab, de_dra, de_drb, &
      df_of_zeta_dna, df_of_zeta_dnb, df_of_zeta_dz, DFP, dh_dra, dh_drb, &
      dh_dt, dipo_alpha_dna, dipo_alpha_dnb, dipo_delta_dna, dipo_delta_dnb, &
      dphi_da, dphi_db, dphi_de, ds_da, ds_dra, ds_drb, ds_dt, dt_dphi, &
      dt_dra, dt_drb, dzeta_4_dna, dzeta_4_dnb, dzeta_dna, dzeta_dnb, ecpw, &
      expe, f_of_zeta, g_ab, h0, ipo_alpha, ipo_delta, n_a, n_ab, n_ab_2, &
      n_b, phi, phi_3, rs, s1, sc, scpw, t, t_2, v1apw, v1bpw, v1ca, v1cb, &
      v2ca, v2cab, v2cb, xkf, xks, xy, y, y_2, zeta, zeta_2, zeta_3, zeta_4
    TYPE(cp_xc_functional_t)                 :: pw_alpha, pw_F, pw_P

!
! LSDA/PW correlation contribution
!
! \Delta (e_F - e_P)
!
! GGA/PBE correlation contribution
!

    n_a  = scratch(alpha)%n
    n_b  = scratch(beta)%n
    n_ab = scratch(alpha_beta)%n
    rs   = tqbp_to_one_third/(scratch(alpha_beta)%n_1_3)

    !
    ! PW (semi-inlined, keeps us from recomputing certain terms)
    !
    CALL cp_lda_c_pw_P( n_ab, rs, pw_P )
    CALL cp_lda_c_pw_F( n_ab, rs, pw_F )
    CALL cp_lda_c_pw_alpha_fit( n_ab, rs, pw_alpha )

    zeta            = zeta_of( scratch(:) )
    zeta_2          = zeta   * zeta
    zeta_3          = zeta_2 * zeta
    zeta_4          = zeta_2 * zeta_2

    n_ab_2          = n_ab * n_ab

    dzeta_dna       =  2.0_real_8 * scratch(beta)%n  / n_ab_2
    dzeta_dnb       = -2.0_real_8 * scratch(alpha)%n / n_ab_2

    dzeta_4_dna     = 4.0_real_8*zeta_3*dzeta_dna
    dzeta_4_dnb     = 4.0_real_8*zeta_3*dzeta_dnb

    f_of_zeta       = f_of( zeta )
    df_of_zeta_dz   = df_dz_of( zeta )
    df_of_zeta_dna  = df_of_zeta_dz * dzeta_dna
    df_of_zeta_dnb  = df_of_zeta_dz * dzeta_dnb

    DFP             = pw_F%sc  - pw_P%sc
    dDFP_dnab       = pw_F%v1c - pw_P%v1c

    ipo_alpha       = - ( f_of_zeta / fpp * (1.0_real_8 - zeta_4) )
    dipo_alpha_dna  = - ( df_of_zeta_dna / fpp * (1.0_real_8 - zeta_4) &
         - f_of_zeta / fpp * (dzeta_4_dna) )
    dipo_alpha_dnb  = - ( df_of_zeta_dnb / fpp * (1.0_real_8 - zeta_4) &
         - f_of_zeta / fpp * (dzeta_4_dnb) )

    ipo_delta       = f_of_zeta * zeta_4
    dipo_delta_dna  = df_of_zeta_dna * zeta_4 + f_of_zeta * dzeta_4_dna
    dipo_delta_dnb  = df_of_zeta_dnb * zeta_4 + f_of_zeta * dzeta_4_dnb

    scpw            = pw_P%sc + ipo_delta*DFP + ipo_alpha*pw_alpha%sc
    ecpw            = scpw/n_ab

    v1apw           = pw_P%v1c + dipo_delta_dna*DFP + ipo_delta*dDFP_dnab &
         + dipo_alpha_dna*pw_alpha%sc  + ipo_alpha*pw_alpha%v1c
    v1bpw           = pw_P%v1c + dipo_delta_dnb*DFP + ipo_delta*dDFP_dnab &
         + dipo_alpha_dnb*pw_alpha%sc  + ipo_alpha*pw_alpha%v1c

    !
    ! PBEc
    !
    phi   = phi_of( zeta )
    phi_3 = phi*phi*phi
    g_ab  = scratch(alpha)%g_2 + 2._real_8*scratch(alpha_beta)%g_2 + scratch(beta)%g_2

    aa    = MAX(g_ab,small)
    a     = SQRT(aa)
    xkf   = xk/rs
    xks   = SQRT(4._real_8*xkf/pi)
    t     = a/(2._real_8*xks*n_ab*phi)
    t_2   = t*t
    expe  = EXP( -ecpw/(phi_3*pbe_gamma))
    af    = pbe_beta/pbe_gamma * (1._real_8/(expe - 1._real_8))
    y     = af*t_2
    y_2   = y*y
    xy    = (1._real_8 + y)/(1._real_8 + y + y_2)
    s1    = 1._real_8 + pbe_beta/pbe_gamma*t_2*xy
    h0    = pbe_gamma*phi_3 * LOG(s1)
    sc    = n_ab*h0

    IF (zeta > -1.0_real_8 .AND. zeta < 1.0_real_8) THEN
       dt_dphi = -t/phi
       dphi_de =  1._real_8/(3._real_8*(1._real_8 + zeta)**one_third) - 1._real_8/(3._real_8*(1._real_8 - zeta)**one_third)
       de_dra  =  2._real_8*n_b/(n_ab_2)
       de_drb  = -2._real_8*n_a/(n_ab_2)
       dphi_da =  dphi_de*de_dra
       dphi_db =  dphi_de*de_drb
       dt_dra  = -t*(dphi_da/phi + 7._real_8/(6._real_8*n_ab))
       dt_drb  = -t*(dphi_db/phi + 7._real_8/(6._real_8*n_ab))
       da_dra  =  af*af*expe/( - pbe_beta*phi_3)*(3._real_8*ecpw/phi*dphi_da - (v1apw - ecpw)/n_ab)
       da_drb  =  af*af*expe/( - pbe_beta*phi_3)*(3._real_8*ecpw/phi*dphi_db - (v1bpw - ecpw)/n_ab)
       ds_da   = -pbe_beta/pbe_gamma * af * t_2**3 * (2._real_8 + y) / (1._real_8 + y + y_2)**2
       ds_dt   =  2._real_8*pbe_beta/pbe_gamma * t * (1._real_8 + 2._real_8*y) / (1._real_8 + y + y_2)**2
       ds_dra  =  ds_da*da_dra  +  ds_dt*dt_dra
       ds_drb  =  ds_da*da_drb  +  ds_dt*dt_drb
       dh_dt   =  pbe_gamma*phi_3/s1*ds_dt
       dh_dra  =  3._real_8*h0/phi*dphi_da  +  pbe_gamma*phi_3/s1*ds_dra
       dh_drb  =  3._real_8*h0/phi*dphi_db  +  pbe_gamma*phi_3/s1*ds_drb
    ELSE
       dt_dra  = -t*(7._real_8/(6._real_8*n_ab))
       dt_drb  = -t*(7._real_8/(6._real_8*n_ab))
       da_dra  =  af*af*expe/( -pbe_beta*phi_3)*( -(v1apw - ecpw)/n_ab)
       da_drb  =  af*af*expe/( -pbe_beta*phi_3)*( -(v1bpw - ecpw)/n_ab)
       ds_da   = -pbe_beta/pbe_gamma * af * t_2**3 * (2._real_8 + y) / (1._real_8 + y + y_2)**2
       ds_dt   =  2._real_8*pbe_beta/pbe_gamma * t * (1._real_8 + 2._real_8*y) / (1._real_8 + y + y_2)**2
       ds_dra  =  ds_da*da_dra  +  ds_dt*dt_dra
       ds_drb  =  ds_da*da_drb  +  ds_dt*dt_drb
       dh_dt   =  pbe_gamma*phi_3/s1*ds_dt
       dh_dra  =  pbe_gamma*phi_3/s1*ds_dra
       dh_drb  =  pbe_gamma*phi_3/s1*ds_drb
    ENDIF

    v1ca   =  h0  +  n_ab*dh_dra
    v1cb   =  h0  +  n_ab*dh_drb
    v2ca   =  n_ab*dh_dt*t/aa
    v2cb   =  v2ca 
    v2cab  =  v2ca

    functional(alpha)%sc       = 0.0_real_8
    functional(beta)%sc        = 0.0_real_8
    functional(alpha_beta)%sc  = sc + scpw

    functional(alpha)%v1c      = v1ca + v1apw
    functional(beta)%v1c       = v1cb + v1bpw
    functional(alpha_beta)%v1c = 0.0_real_8

    functional(alpha)%v2c      = v2ca 
    functional(beta)%v2c       = v2cb
    functional(alpha_beta)%v2c = v2cab

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_spin_gga_c_pbe
  ! ==================================================================
  PURE SUBROUTINE cp_gga_c_pbe(scratch,functional)
    ! ==--------------------------------------------------------------==
    ! PBE correlation using PW as in the original PBE routine.
    !                                 14.10.2016 M.P.Bircher @ LCBC/EPFL
    ! ==--------------------------------------------------------------==


    TYPE(cp_xc_scratch_t), INTENT(in)        :: scratch
    TYPE(cp_xc_functional_t), INTENT(inout)  :: functional

    REAL(real_8), PARAMETER :: pbe_beta = 0.06672455060314922_real_8, &
      pbe_gamma = 0.031090690869654895_real_8 , &
      beta_by_gamma = pbe_beta/pbe_gamma , pw_a = 0.0310907_real_8, &
      pw_a1 = 0.21370_real_8, pw_b1 = 7.5957_real_8, pw_b2 = 3.5876_real_8, &
      pw_b3 = 1.6382_real_8, pw_b4 = 0.49294_real_8, pw_c0 = pw_a, &
      pw_c1 = 0.046644_real_8, pw_c2 = 0.00664_real_8, &
      pw_c3 = 0.01043_real_8, pw_d0 = 0.4335_real_8, pw_d1 = 1.4408_real_8

    REAL(real_8) :: a, a1rs, aa, af, b1rs1, b2rs2, b3rs3, b4rs4, c0xln, c2rs, &
      d0_rs, dadr, dhdr, dhdt, dom, dsda, dsdr, dsdt, dtdr, epwc, expe, h0, &
      olog, om, rho, rs, rs1, rs3, rs4, s1, sc, spwc, t, t_2, v1c, v2c, vpwc, &
      xkf, xks, xln, xy, y, y_2, yden, yden_2

    rs   = tqbp_to_one_third/(two_to_one_third*scratch%n_1_3)
    aa   = 4.0_real_8*scratch%g_2
    a    = 2.0_real_8*scratch%abs_g
    rho  = 2.0_real_8*scratch%n

    !
    ! PW (inlined - just because)
    !
    IF (rs < 0.5_real_8) THEN
       !
       ! High density formula
       !
       xln   = LOG(rs)
       c0xln = pw_c0*xln
       c2rs  = pw_c2*rs
       epwc  = c0xln-pw_c1 + c2rs*xln-pw_c3*rs
       spwc  = 2.0_real_8*scratch%n*epwc
       vpwc  = c0xln-(pw_c1 + pw_c0/3._real_8) + 2._real_8/3._real_8*c2rs*xln &
            -(2._real_8*pw_c3 + pw_c2)/3._real_8*rs
    ELSE IF (rs > 100._real_8) THEN
       !
       ! Low density formula
       !
       d0_rs = pw_d0/rs
       epwc  = -d0_rs + pw_d1/rs**1.5_real_8
       spwc  = 2.0_real_8*scratch%n*epwc
       vpwc  = -4._real_8/3._real_8*d0_rs + 1.5_real_8*pw_d1/(rs*SQRT(rs))
    ELSE
       !
       ! Interpolation formula
       !
       rs1   = SQRT(rs)
       rs3   = rs*rs1
       rs4   = rs*rs
       b1rs1 = pw_b1*rs1
       b2rs2 = pw_b2*rs
       b3rs3 = pw_b3*rs3
       b4rs4 = pw_b4*rs4
       a1rs  = pw_a1*rs
       om    = 2._real_8*pw_a*(b1rs1 + b2rs2 + b3rs3 + b4rs4)
       dom   = 2._real_8*pw_a*(0.5_real_8*b1rs1 + b2rs2 + 1.5_real_8*b3rs3 + 2._real_8*b4rs4)
       olog  = LOG(1._real_8 + 1.0_real_8/om)
       epwc  = -2._real_8*pw_a*(1.0_real_8 + a1rs)*olog
       spwc  = 2.0_real_8*scratch%n*epwc
       vpwc  = -2._real_8*pw_a*(1._real_8 + 2._real_8/3._real_8*a1rs)*olog &
            -2._real_8/3._real_8*pw_a*(1._real_8 + a1rs)*dom/(om*(om + 1._real_8))
    END IF

    !
    ! PBE-GGAc
    !
    xkf    = nqtp_to_one_third/rs
    xks    = SQRT(4._real_8*xkf/pi)
    t      = a/(2._real_8*xks*rho)
    t_2    = t*t
    expe   = EXP(-epwc/pbe_gamma)
    af     = beta_by_gamma * (1._real_8/(expe-1._real_8))
    y      = af*t_2
    y_2    = y*y
    yden   = (1._real_8 + y + y_2)
    yden_2 = yden*yden
    xy     = (1._real_8 + y)/yden
    s1     = 1._real_8 + beta_by_gamma*t_2*xy
    h0     = pbe_gamma * LOG(s1)
    dtdr   = -7._real_8*t/(6._real_8*rho)
    dadr   = af*af*expe/pbe_beta * (vpwc-epwc)/rho
    dsda   = -beta_by_gamma * af * t_2 * t_2 * t_2 * (2._real_8 + y) / yden_2
    dsdt   = 2._real_8*beta_by_gamma * t * (1._real_8 + 2._real_8*y) / yden_2
    dsdr   = dsda*dadr  +  dsdt*dtdr
    dhdt   = pbe_gamma/s1*dsdt
    dhdr   = pbe_gamma/s1*dsdr
    sc     = rho*h0
    v1c    = h0 + dhdr*rho
    v2c    = rho*dhdt*t/aa
    !
    functional%sc  = spwc + sc
    functional%v1c = vpwc + v1c
    functional%v2c = v2c
    !

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_gga_c_pbe
  ! ==================================================================
  PURE SUBROUTINE cp_spin_gga_c_pbe_sol(scratch,functional)
    ! ==--------------------------------------------------------------==
    ! PBEsol correlation using PW as in the original PBE routine.
    ! Copy from cp_spin_gga_c_pbe with changed parameters.
    !
    !                                 12.04.2017 M.P.Bircher @ LCBC/EPFL
    ! ==--------------------------------------------------------------==


    TYPE(cp_xc_scratch_t), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(in)                             :: scratch
    TYPE(cp_xc_functional_t), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(inout)                          :: functional

    REAL(real_8), PARAMETER :: fpp = 1.709921_real_8, &
      pbe_beta = 0.0460_real_8, pbe_gamma = 0.031090690869654895_real_8 , &
      beta_by_gamma = pbe_beta/pbe_gamma , small = 1.0e-24_real_8, &
      xk = (9._real_8*pi/4._real_8)**one_third

    REAL(real_8) :: a, aa, af, da_dra, da_drb, dDFP_dnab, de_dra, de_drb, &
      df_of_zeta_dna, df_of_zeta_dnb, df_of_zeta_dz, DFP, dh_dra, dh_drb, &
      dh_dt, dipo_alpha_dna, dipo_alpha_dnb, dipo_delta_dna, dipo_delta_dnb, &
      dphi_da, dphi_db, dphi_de, ds_da, ds_dra, ds_drb, ds_dt, dt_dphi, &
      dt_dra, dt_drb, dzeta_4_dna, dzeta_4_dnb, dzeta_dna, dzeta_dnb, ecpw, &
      expe, f_of_zeta, g_ab, h0, ipo_alpha, ipo_delta, n_a, n_ab, n_ab_2, &
      n_b, phi, phi_3, rs, s1, sc, scpw, t, t_2, v1apw, v1bpw, v1ca, v1cb, &
      v2ca, v2cab, v2cb, xkf, xks, xy, y, y_2, zeta, zeta_2, zeta_3, zeta_4
    TYPE(cp_xc_functional_t)                 :: pw_alpha, pw_F, pw_P

!
! LSDA/PW correlation contribution
!
! \Delta (e_F - e_P)
!
! GGA/PBE correlation contribution
!

    n_a  = scratch(alpha)%n
    n_b  = scratch(beta)%n
    n_ab = scratch(alpha_beta)%n
    rs   = tqbp_to_one_third/(scratch(alpha_beta)%n_1_3)

    !
    ! PW (semi-inlined, keeps us from recomputing certain terms)
    !
    CALL cp_lda_c_pw_P( n_ab, rs, pw_P )
    CALL cp_lda_c_pw_F( n_ab, rs, pw_F )
    CALL cp_lda_c_pw_alpha_fit( n_ab, rs, pw_alpha )

    zeta            = zeta_of( scratch(:) )
    zeta_2          = zeta   * zeta
    zeta_3          = zeta_2 * zeta
    zeta_4          = zeta_2 * zeta_2

    n_ab_2          = n_ab * n_ab

    dzeta_dna       =  2.0_real_8 * scratch(beta)%n  / n_ab_2
    dzeta_dnb       = -2.0_real_8 * scratch(alpha)%n / n_ab_2

    dzeta_4_dna     = 4.0_real_8*zeta_3*dzeta_dna
    dzeta_4_dnb     = 4.0_real_8*zeta_3*dzeta_dnb

    f_of_zeta       = f_of( zeta )
    df_of_zeta_dz   = df_dz_of( zeta )
    df_of_zeta_dna  = df_of_zeta_dz * dzeta_dna
    df_of_zeta_dnb  = df_of_zeta_dz * dzeta_dnb

    DFP             = pw_F%sc  - pw_P%sc
    dDFP_dnab       = pw_F%v1c - pw_P%v1c

    ipo_alpha       = - ( f_of_zeta / fpp * (1.0_real_8 - zeta_4) )
    dipo_alpha_dna  = - ( df_of_zeta_dna / fpp * (1.0_real_8 - zeta_4) &
         - f_of_zeta / fpp * (dzeta_4_dna) )
    dipo_alpha_dnb  = - ( df_of_zeta_dnb / fpp * (1.0_real_8 - zeta_4) &
         - f_of_zeta / fpp * (dzeta_4_dnb) )

    ipo_delta       = f_of_zeta * zeta_4
    dipo_delta_dna  = df_of_zeta_dna * zeta_4 + f_of_zeta * dzeta_4_dna
    dipo_delta_dnb  = df_of_zeta_dnb * zeta_4 + f_of_zeta * dzeta_4_dnb

    scpw            = pw_P%sc + ipo_delta*DFP + ipo_alpha*pw_alpha%sc
    ecpw            = scpw/n_ab

    v1apw           = pw_P%v1c + dipo_delta_dna*DFP + ipo_delta*dDFP_dnab &
         + dipo_alpha_dna*pw_alpha%sc  + ipo_alpha*pw_alpha%v1c
    v1bpw           = pw_P%v1c + dipo_delta_dnb*DFP + ipo_delta*dDFP_dnab &
         + dipo_alpha_dnb*pw_alpha%sc  + ipo_alpha*pw_alpha%v1c

    !
    ! PBEc
    !
    phi   = phi_of( zeta )
    phi_3 = phi*phi*phi
    g_ab  = scratch(alpha)%g_2 + 2._real_8*scratch(alpha_beta)%g_2 + scratch(beta)%g_2

    aa    = MAX(g_ab,small)
    a     = SQRT(aa)
    xkf   = xk/rs
    xks   = SQRT(4._real_8*xkf/pi)
    t     = a/(2._real_8*xks*n_ab*phi)
    t_2   = t*t
    expe  = EXP( -ecpw/(phi_3*pbe_gamma))
    af    = pbe_beta/pbe_gamma * (1._real_8/(expe - 1._real_8))
    y     = af*t_2
    y_2   = y*y
    xy    = (1._real_8 + y)/(1._real_8 + y + y_2)
    s1    = 1._real_8 + pbe_beta/pbe_gamma*t_2*xy
    h0    = pbe_gamma*phi_3 * LOG(s1)
    sc    = n_ab*h0

    IF (zeta > -1.0_real_8 .AND. zeta < 1.0_real_8) THEN
       dt_dphi = -t/phi
       dphi_de =  1._real_8/(3._real_8*(1._real_8 + zeta)**one_third) - 1._real_8/(3._real_8*(1._real_8 - zeta)**one_third)
       de_dra  =  2._real_8*n_b/(n_ab_2)
       de_drb  = -2._real_8*n_a/(n_ab_2)
       dphi_da =  dphi_de*de_dra
       dphi_db =  dphi_de*de_drb
       dt_dra  = -t*(dphi_da/phi + 7._real_8/(6._real_8*n_ab))
       dt_drb  = -t*(dphi_db/phi + 7._real_8/(6._real_8*n_ab))
       da_dra  =  af*af*expe/( - pbe_beta*phi_3)*(3._real_8*ecpw/phi*dphi_da - (v1apw - ecpw)/n_ab)
       da_drb  =  af*af*expe/( - pbe_beta*phi_3)*(3._real_8*ecpw/phi*dphi_db - (v1bpw - ecpw)/n_ab)
       ds_da   = -pbe_beta/pbe_gamma * af * t_2**3 * (2._real_8 + y) / (1._real_8 + y + y_2)**2
       ds_dt   =  2._real_8*pbe_beta/pbe_gamma * t * (1._real_8 + 2._real_8*y) / (1._real_8 + y + y_2)**2
       ds_dra  =  ds_da*da_dra  +  ds_dt*dt_dra
       ds_drb  =  ds_da*da_drb  +  ds_dt*dt_drb
       dh_dt   =  pbe_gamma*phi_3/s1*ds_dt
       dh_dra  =  3._real_8*h0/phi*dphi_da  +  pbe_gamma*phi_3/s1*ds_dra
       dh_drb  =  3._real_8*h0/phi*dphi_db  +  pbe_gamma*phi_3/s1*ds_drb
    ELSE
       dt_dra  = -t*(7._real_8/(6._real_8*n_ab))
       dt_drb  = -t*(7._real_8/(6._real_8*n_ab))
       da_dra  =  af*af*expe/( -pbe_beta*phi_3)*( -(v1apw - ecpw)/n_ab)
       da_drb  =  af*af*expe/( -pbe_beta*phi_3)*( -(v1bpw - ecpw)/n_ab)
       ds_da   = -pbe_beta/pbe_gamma * af * t_2**3 * (2._real_8 + y) / (1._real_8 + y + y_2)**2
       ds_dt   =  2._real_8*pbe_beta/pbe_gamma * t * (1._real_8 + 2._real_8*y) / (1._real_8 + y + y_2)**2
       ds_dra  =  ds_da*da_dra  +  ds_dt*dt_dra
       ds_drb  =  ds_da*da_drb  +  ds_dt*dt_drb
       dh_dt   =  pbe_gamma*phi_3/s1*ds_dt
       dh_dra  =  pbe_gamma*phi_3/s1*ds_dra
       dh_drb  =  pbe_gamma*phi_3/s1*ds_drb
    ENDIF

    v1ca   =  h0  +  n_ab*dh_dra
    v1cb   =  h0  +  n_ab*dh_drb
    v2ca   =  n_ab*dh_dt*t/aa
    v2cb   =  v2ca 
    v2cab  =  v2ca

    functional(alpha)%sc       = 0.0_real_8
    functional(beta)%sc        = 0.0_real_8
    functional(alpha_beta)%sc  = sc + scpw

    functional(alpha)%v1c      = v1ca + v1apw
    functional(beta)%v1c       = v1cb + v1bpw
    functional(alpha_beta)%v1c = 0.0_real_8

    functional(alpha)%v2c      = v2ca 
    functional(beta)%v2c       = v2cb
    functional(alpha_beta)%v2c = v2cab

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_spin_gga_c_pbe_sol
  ! ==================================================================
  PURE SUBROUTINE cp_gga_c_pbe_sol(scratch,functional)
    ! ==--------------------------------------------------------------==
    ! PBEsol correlation using PW as in the original PBE routine.
    !                                 14.10.2016 M.P.Bircher @ LCBC/EPFL
    ! ==--------------------------------------------------------------==


    TYPE(cp_xc_scratch_t), INTENT(in)        :: scratch
    TYPE(cp_xc_functional_t), INTENT(inout)  :: functional

    REAL(real_8), PARAMETER :: pbe_beta = 0.0460_real_8, &
      pbe_gamma = 0.031090690869654895_real_8 , &
      beta_by_gamma = pbe_beta/pbe_gamma , pw_a = 0.0310907_real_8, &
      pw_a1 = 0.21370_real_8, pw_b1 = 7.5957_real_8, pw_b2 = 3.5876_real_8, &
      pw_b3 = 1.6382_real_8, pw_b4 = 0.49294_real_8, pw_c0 = pw_a, &
      pw_c1 = 0.046644_real_8, pw_c2 = 0.00664_real_8, &
      pw_c3 = 0.01043_real_8, pw_d0 = 0.4335_real_8, pw_d1 = 1.4408_real_8

    REAL(real_8) :: a, a1rs, aa, af, b1rs1, b2rs2, b3rs3, b4rs4, c0xln, c2rs, &
      d0_rs, dadr, dhdr, dhdt, dom, dsda, dsdr, dsdt, dtdr, epwc, expe, h0, &
      olog, om, rho, rs, rs1, rs3, rs4, s1, sc, spwc, t, t_2, v1c, v2c, vpwc, &
      xkf, xks, xln, xy, y, y_2, yden, yden_2

    rs   = tqbp_to_one_third/(two_to_one_third*scratch%n_1_3)
    aa   = 4.0_real_8*scratch%g_2
    a    = 2.0_real_8*scratch%abs_g
    rho  = 2.0_real_8*scratch%n

    !
    ! PW
    !
    IF (rs < 0.5_real_8) THEN
       !
       ! High density formula
       !
       xln   = LOG(rs)
       c0xln = pw_c0*xln
       c2rs  = pw_c2*rs
       epwc  = c0xln-pw_c1 + c2rs*xln-pw_c3*rs
       spwc  = 2.0_real_8*scratch%n*epwc
       vpwc  = c0xln-(pw_c1 + pw_c0/3._real_8) + 2._real_8/3._real_8*c2rs*xln &
            -(2._real_8*pw_c3 + pw_c2)/3._real_8*rs
    ELSE IF (rs > 100._real_8) THEN
       !
       ! Low density formula
       !
       d0_rs = pw_d0/rs
       epwc  = -d0_rs + pw_d1/rs**1.5_real_8
       spwc  = 2.0_real_8*scratch%n*epwc
       vpwc  = -4._real_8/3._real_8*d0_rs + 1.5_real_8*pw_d1/(rs*SQRT(rs))
    ELSE
       !
       ! Interpolation formula
       !
       rs1   = SQRT(rs)
       rs3   = rs*rs1
       rs4   = rs*rs
       b1rs1 = pw_b1*rs1
       b2rs2 = pw_b2*rs
       b3rs3 = pw_b3*rs3
       b4rs4 = pw_b4*rs4
       a1rs  = pw_a1*rs
       om    = 2._real_8*pw_a*(b1rs1 + b2rs2 + b3rs3 + b4rs4)
       dom   = 2._real_8*pw_a*(0.5_real_8*b1rs1 + b2rs2 + 1.5_real_8*b3rs3 + 2._real_8*b4rs4)
       olog  = LOG(1._real_8 + 1.0_real_8/om)
       epwc  = -2._real_8*pw_a*(1.0_real_8 + a1rs)*olog
       spwc  = 2.0_real_8*scratch%n*epwc
       vpwc  = -2._real_8*pw_a*(1._real_8 + 2._real_8/3._real_8*a1rs)*olog &
            -2._real_8/3._real_8*pw_a*(1._real_8 + a1rs)*dom/(om*(om + 1._real_8))
    END IF

    !
    ! PBE-GGAc
    !
    xkf    = nqtp_to_one_third/rs
    xks    = SQRT(4._real_8*xkf/pi)
    t      = a/(2._real_8*xks*rho)
    t_2    = t*t
    expe   = EXP(-epwc/pbe_gamma)
    af     = beta_by_gamma * (1._real_8/(expe-1._real_8))
    y      = af*t_2
    y_2    = y*y
    yden   = (1._real_8 + y + y_2)
    yden_2 = yden*yden
    xy     = (1._real_8 + y)/yden
    s1     = 1._real_8 + beta_by_gamma*t_2*xy
    h0     = pbe_gamma * LOG(s1)
    dtdr   = -7._real_8*t/(6._real_8*rho)
    dadr   = af*af*expe/pbe_beta * (vpwc-epwc)/rho
    dsda   = -beta_by_gamma * af * t_2 * t_2 * t_2 * (2._real_8 + y) / yden_2
    dsdt   = 2._real_8*beta_by_gamma * t * (1._real_8 + 2._real_8*y) / yden_2
    dsdr   = dsda*dadr  +  dsdt*dtdr
    dhdt   = pbe_gamma/s1*dsdt
    dhdr   = pbe_gamma/s1*dsdr
    sc     = rho*h0
    v1c    = h0 + dhdr*rho
    v2c    = rho*dhdt*t/aa
    !
    functional%sc  = spwc + sc
    functional%v1c = vpwc + v1c
    functional%v2c = v2c
    !

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_gga_c_pbe_sol
  ! ==================================================================
  PURE SUBROUTINE cp_spin_gga_c_pbe_flex(scratch,functional)
    ! ==--------------------------------------------------------------==
    ! Mixing routine for PBE correlation.
    ! Backwards-compatible with old CPMD. Energies differ from libxc;
    !
    !           PBEC in CPMD uses PZ, PBEC in libxc uses PW
    !
    ! energies are identical to OLDCODE, but final energies may differ in
    ! the last digit due to the GC-cutoff being applied to the LDA terms
    !
    !                                 13.04.2017 M.P.Bircher @ LCBC/EPFL
    ! ==--------------------------------------------------------------==


    TYPE(cp_xc_scratch_t), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(in)                             :: scratch
    TYPE(cp_xc_functional_t), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(inout)                          :: functional

    REAL(real_8), PARAMETER :: fpp = 1.709921_real_8, small = 1.0e-24_real_8, &
      xk = (9._real_8*pi/4._real_8)**one_third

    REAL(real_8) :: a, aa, af, beta_by_gamma, da_dra, da_drb, de_dra, de_drb, &
      dh_dra, dh_drb, dh_dt, dphi_da, dphi_db, dphi_de, ds_da, ds_dra, &
      ds_drb, ds_dt, dt_dphi, dt_dra, dt_drb, eclsd, expe, g_ab, h0, n_a, &
      n_ab, n_ab_2, n_b, pbe_beta, pbe_gamma, phi, phi_3, rs, s1, sc, sclsd, &
      t, t_2, v1alsd, v1blsd, v1ca, v1cb, v2ca, v2cab, v2cb, xkf, xks, xy, y, &
      y_2, zeta
    TYPE(cp_xc_functional_t), &
      DIMENSION(cp_xc_spin_components)       :: lsd

    pbe_beta      = cp_gga_c_param%pbe%pbe_beta
    pbe_gamma     = cp_gga_c_param%pbe%pbe_gamma
    beta_by_gamma = pbe_beta / pbe_gamma

    n_a    = scratch(alpha)%n
    n_b    = scratch(beta)%n
    n_ab   = scratch(alpha_beta)%n
    n_ab_2 = n_ab * n_ab
    rs     = tqbp_to_one_third/(scratch(alpha_beta)%n_1_3)
    zeta   = zeta_of( scratch(:) )

    ! Call any correlation functional for the LDA part
    !
    SELECT CASE(TRIM(ADJUSTL(cp_gga_c_param%pbe%lda_c)))
    CASE("CP_LDA_C_VWN")
       CALL cp_spin_lda_c_vwn( scratch(:),lsd(:) )
    CASE("CP_LDA_C_PZ")
       CALL cp_spin_lda_c_pz( scratch(:),lsd(:) )
    CASE("CP_LDA_C_OB_PW")
       CALL cp_spin_lda_c_ob_pw( scratch(:),lsd(:) )
    CASE default
       CALL cp_spin_lda_c_pw( scratch(:), lsd(:) )
    END SELECT
    sclsd  = lsd(alpha_beta)%sc
    eclsd  = lsd(alpha_beta)%sc/n_ab
    v1alsd = lsd(alpha)%v1c
    v1blsd = lsd(beta)%v1c

    phi   = phi_of( zeta )
    phi_3 = phi*phi*phi
    g_ab  = scratch(alpha)%g_2 + 2._real_8*scratch(alpha_beta)%g_2 + scratch(beta)%g_2

    aa    = MAX(g_ab,small)
    a     = SQRT(aa)
    xkf   = xk/rs
    xks   = SQRT(4._real_8*xkf/pi)
    t     = a/(2._real_8*xks*n_ab*phi)
    t_2   = t*t
    expe  = EXP( -eclsd/(phi_3*pbe_gamma))
    af    = pbe_beta/pbe_gamma * (1._real_8/(expe - 1._real_8))
    y     = af*t_2
    y_2   = y*y
    xy    = (1._real_8 + y)/(1._real_8 + y + y_2)
    s1    = 1._real_8 + pbe_beta/pbe_gamma*t_2*xy
    h0    = pbe_gamma*phi_3 * LOG(s1)
    sc    = n_ab*h0

    IF (zeta > -1.0_real_8 .AND. zeta < 1.0_real_8) THEN
       dt_dphi = -t/phi
       dphi_de =  1._real_8/(3._real_8*(1._real_8 + zeta)**one_third) - 1._real_8/(3._real_8*(1._real_8 - zeta)**one_third)
       de_dra  =  2._real_8*n_b/(n_ab_2)
       de_drb  = -2._real_8*n_a/(n_ab_2)
       dphi_da =  dphi_de*de_dra
       dphi_db =  dphi_de*de_drb
       dt_dra  = -t*(dphi_da/phi + 7._real_8/(6._real_8*n_ab))
       dt_drb  = -t*(dphi_db/phi + 7._real_8/(6._real_8*n_ab))
       da_dra  =  af*af*expe/( - pbe_beta*phi_3)*(3._real_8*eclsd/phi*dphi_da - (v1alsd - eclsd)/n_ab)
       da_drb  =  af*af*expe/( - pbe_beta*phi_3)*(3._real_8*eclsd/phi*dphi_db - (v1blsd - eclsd)/n_ab)
       ds_da   = -pbe_beta/pbe_gamma * af * t_2**3 * (2._real_8 + y) / (1._real_8 + y + y_2)**2
       ds_dt   =  2._real_8*pbe_beta/pbe_gamma * t * (1._real_8 + 2._real_8*y) / (1._real_8 + y + y_2)**2
       ds_dra  =  ds_da*da_dra  +  ds_dt*dt_dra
       ds_drb  =  ds_da*da_drb  +  ds_dt*dt_drb
       dh_dt   =  pbe_gamma*phi_3/s1*ds_dt
       dh_dra  =  3._real_8*h0/phi*dphi_da  +  pbe_gamma*phi_3/s1*ds_dra
       dh_drb  =  3._real_8*h0/phi*dphi_db  +  pbe_gamma*phi_3/s1*ds_drb
    ELSE
       dt_dra  = -t*(7._real_8/(6._real_8*n_ab))
       dt_drb  = -t*(7._real_8/(6._real_8*n_ab))
       da_dra  =  af*af*expe/( -pbe_beta*phi_3)*( -(v1alsd - eclsd)/n_ab)
       da_drb  =  af*af*expe/( -pbe_beta*phi_3)*( -(v1blsd - eclsd)/n_ab)
       ds_da   = -pbe_beta/pbe_gamma * af * t_2**3 * (2._real_8 + y) / (1._real_8 + y + y_2)**2
       ds_dt   =  2._real_8*pbe_beta/pbe_gamma * t * (1._real_8 + 2._real_8*y) / (1._real_8 + y + y_2)**2
       ds_dra  =  ds_da*da_dra  +  ds_dt*dt_dra
       ds_drb  =  ds_da*da_drb  +  ds_dt*dt_drb
       dh_dt   =  pbe_gamma*phi_3/s1*ds_dt
       dh_dra  =  pbe_gamma*phi_3/s1*ds_dra
       dh_drb  =  pbe_gamma*phi_3/s1*ds_drb
    ENDIF

    v1ca   =  h0  +  n_ab*dh_dra
    v1cb   =  h0  +  n_ab*dh_drb
    v2ca   =  n_ab*dh_dt*t/aa
    v2cb   =  v2ca 
    v2cab  =  v2ca

    functional(alpha)%sc       = 0.0_real_8
    functional(beta)%sc        = 0.0_real_8
    functional(alpha_beta)%sc  = sc + sclsd

    functional(alpha)%v1c      = v1ca + v1alsd
    functional(beta)%v1c       = v1cb + v1blsd
    functional(alpha_beta)%v1c = 0.0_real_8

    functional(alpha)%v2c      = v2ca 
    functional(beta)%v2c       = v2cb
    functional(alpha_beta)%v2c = v2cab

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_spin_gga_c_pbe_flex
  ! ==================================================================
  PURE SUBROUTINE cp_gga_c_pbe_flex(scratch,functional)
    ! ==--------------------------------------------------------------==
    ! Mixing routine for PBE correlation.
    ! Backwards-compatible with old CPMD. Energies differ from libxc;
    !
    !           PBEC in CPMD uses PZ, PBEC in libxc uses PW
    !
    ! energies are identical to OLDCODE, but final energies may differ in
    ! the last digit due to the GC-cutoff being applied to the LDA terms
    !                                 14.10.2016 M.P.Bircher @ LCBC/EPFL
    ! ==--------------------------------------------------------------==


    TYPE(cp_xc_scratch_t), INTENT(in)        :: scratch
    TYPE(cp_xc_functional_t), INTENT(inout)  :: functional

    REAL(real_8) :: a, aa, af, beta_by_gamma, dadr, dhdr, dhdt, dsda, dsdr, &
      dsdt, dtdr, elda, expe, h0, pbe_beta, pbe_gamma, rho, rs, s1, sc, t, &
      t_2, v1c, v2c, xkf, xks, xy, y, y_2, yden, yden_2
    TYPE(cp_xc_functional_t)                 :: lda

    pbe_beta      = cp_gga_c_param%pbe%pbe_beta
    pbe_gamma     = cp_gga_c_param%pbe%pbe_gamma
    beta_by_gamma = pbe_beta / pbe_gamma

    rs   = tqbp_to_one_third/(two_to_one_third*scratch%n_1_3)
    aa   = 4.0_real_8*scratch%g_2
    a    = 2.0_real_8*scratch%abs_g
    rho  = 2.0_real_8*scratch%n

    ! Call any correlation functional for the LDA part
    !
    lda%sc  = 0.0_real_8
    lda%v1c = 0.0_real_8
    SELECT CASE(TRIM(ADJUSTL(cp_gga_c_param%pbe%lda_c)))
    CASE("CP_LDA_C_VWN")
       CALL cp_lda_c_vwn(scratch,lda)
    CASE("CP_LDA_C_PZ")
       CALL cp_lda_c_pz(scratch,lda)
    CASE("CP_LDA_C_OB_PW")
       CALL cp_lda_c_ob_pw(scratch,lda)
    CASE default
       CALL cp_lda_c_pw(scratch,lda)
    END SELECT
    elda   = lda%sc/rho
    !
    ! PBE-GGAc
    !
    xkf    = nqtp_to_one_third/rs
    xks    = SQRT(4._real_8*xkf/pi)
    t      = a/(2._real_8*xks*rho)
    t_2    = t*t
    expe   = EXP(-elda/pbe_gamma)
    af     = beta_by_gamma * (1._real_8/(expe-1._real_8))
    y      = af*t_2
    y_2    = y*y
    yden   = (1._real_8 + y + y_2)
    yden_2 = yden*yden
    xy     = (1._real_8 + y)/yden
    s1     = 1._real_8 + beta_by_gamma*t_2*xy
    h0     = pbe_gamma * LOG(s1)
    dtdr   = -7._real_8*t/(6._real_8*rho)
    dadr   = af*af*expe/pbe_beta * (lda%v1c-elda)/rho
    dsda   = -beta_by_gamma * af * t_2 * t_2 * t_2 * (2._real_8 + y) / yden_2
    dsdt   = 2._real_8*beta_by_gamma * t * (1._real_8 + 2._real_8*y) / yden_2
    dsdr   = dsda*dadr  +  dsdt*dtdr
    dhdt   = pbe_gamma/s1*dsdt
    dhdr   = pbe_gamma/s1*dsdr
    sc     = rho*h0
    v1c    = h0 + dhdr*rho
    v2c    = rho*dhdt*t/aa
    !
    functional%sc  = lda%sc + sc
    functional%v1c = lda%v1c + v1c
    functional%v2c = v2c
    !

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_gga_c_pbe_flex
  ! ==================================================================
  PURE SUBROUTINE cp_spin_gga_c_hcth(scratch,functional)
    ! ==--------------------------------------------------------------==
    ! HCTH correlation rewritten from scratch
    !
    !                                 19.07.2017 M.P.Bircher @ LCBC/EPFL
    ! ==--------------------------------------------------------------==


    TYPE(cp_xc_scratch_t), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(in)                             :: scratch
    TYPE(cp_xc_functional_t), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(inout)                          :: functional

    INTEGER, PARAMETER                       :: aa = alpha
    INTEGER, PARAMETER                       :: bb = beta
    INTEGER, PARAMETER                       :: ab = alpha_beta
    INTEGER, PARAMETER                       :: gm = hcth_gm
    REAL(real_8), PARAMETER                  :: sqrt_small        =  1.0e-12_real_8

    REAL(real_8), &
      DIMENSION(cp_xc_spin_components)       :: lda_sc, lda_v1a, lda_v1b, &
                                                sc, dsc_dna, dsc_dnb, dsc_dga, dsc_dgb, &
                                                abs_g
    REAL(real_8), &
      DIMENSION(cp_xc_spin_components)       :: hcth_g, s, s_2, denominator, ds_dna, ds_dnb, &
                                                ds_dga, ds_dgb, gamma_s, gamma_s_2, dhcth_1_ds, &
                                                dhcth_g_ds
    REAL(real_8), &
      DIMENSION(gm,cp_xc_spin_components)    :: u_hcth, dhcth_u_ds

    TYPE(cp_xc_scratch_t), &
      DIMENSION(cp_xc_spin_components)       :: parallel
    
    !
    ! Parallel correlation
    !
    IF ( scratch(aa)%n >= sqrt_small ) THEN
       parallel(aa)%n     = scratch(aa)%n
       parallel(aa)%n_1_3 = scratch(aa)%n_1_3
       parallel(aa)%n_4_3 = scratch(aa)%n_4_3
       parallel(bb)%n     = 0.0_real_8 
       parallel(bb)%n_1_3 = 0.0_real_8 
       parallel(bb)%n_4_3 = 0.0_real_8 
       parallel(ab)%n     = scratch(aa)%n
       parallel(ab)%n_1_3 = scratch(aa)%n_1_3
       parallel(ab)%n_4_3 = scratch(aa)%n_4_3
       !
       CALL cp_spin_lda_c_pw( parallel, functional )
       lda_sc(aa)  = functional(ab)%sc
       lda_v1a(aa) = functional(aa)%v1c
       lda_v1b(aa) = 0.0_real_8
    ELSE
       lda_sc(aa)  = 0.0_real_8 
       lda_v1a(aa) = 0.0_real_8 
       lda_v1b(aa) = 0.0_real_8 
    ENDIF
    !
    ! Idem
    !
    IF ( scratch(bb)%n >= sqrt_small ) THEN
       parallel(aa)%n     = 0.0_real_8
       parallel(aa)%n_1_3 = 0.0_real_8
       parallel(aa)%n_4_3 = 0.0_real_8
       parallel(bb)%n     = scratch(bb)%n
       parallel(bb)%n_1_3 = scratch(bb)%n_1_3
       parallel(bb)%n_4_3 = scratch(bb)%n_4_3
       parallel(ab)%n     = scratch(bb)%n
       parallel(ab)%n_1_3 = scratch(bb)%n_1_3
       parallel(ab)%n_4_3 = scratch(bb)%n_4_3
       !
       CALL cp_spin_lda_c_pw( parallel, functional )
       lda_sc(bb)  = functional(ab)%sc
       lda_v1a(bb) = 0.0_real_8
       lda_v1b(bb) = functional(bb)%v1c
    ELSE
       lda_sc(bb)  = 0.0_real_8 
       lda_v1a(bb) = 0.0_real_8 
       lda_v1b(bb) = 0.0_real_8 
    ENDIF
    !
    ! Antiparallel correlation
    ! (no if, since n_ab will always be > 1.0e-10, cf. cp_xc_driver)
    !
    ! Eq. according to Stoll
    !
    CALL cp_spin_lda_c_pw( scratch, functional )
    lda_sc(ab)  = functional(ab)%sc  - lda_sc(aa)  - lda_sc(bb) 
    lda_v1a(ab) = functional(aa)%v1c - lda_v1a(aa) - lda_v1a(bb)
    lda_v1b(ab) = functional(bb)%v1c - lda_v1b(aa) - lda_v1b(bb)

    !
    ! HCTH preliminaries
    !
    abs_g(:)    = max( scratch(:)%abs_g, sqrt_small )
    s(aa:bb)    = abs_g(aa:bb) / scratch(aa:bb)%n_4_3
    s_2(aa:bb)  = s(aa:bb) * s(aa:bb) 
    s_2(ab)     = 0.5_real_8 * (s_2(aa) + s_2(bb))
    s(ab)       = sqrt( s_2(ab) )
    !
    ds_dna(aa)  = -four_thirds*s(aa) / scratch(aa)%n
    ds_dna(bb)  = 0.0_real_8 
    ds_dna(ab)  = 0.5_real_8 * s(aa) * ds_dna(aa) / s(ab)
    !
    ds_dnb(aa)  = 0.0_real_8 
    ds_dnb(bb)  = -four_thirds*s(bb) / scratch(bb)%n
    ds_dnb(ab)  = 0.5_real_8 * s(bb) * ds_dnb(bb) / s(ab)
    !
    ds_dga(aa)  = 1.0_real_8 / scratch(aa)%n_4_3 
    ds_dga(bb)  = 0.0_real_8 
    ds_dga(ab)  = 0.5_real_8 * s(aa)* ds_dga(aa) / s(ab)
    !
    ds_dgb(aa)  = 0.0_real_8
    ds_dgb(bb)  = 1.0_real_8 / scratch(bb)%n_4_3 
    ds_dgb(ab)  = 0.5_real_8 * s(bb)* ds_dgb(bb) / s(ab)
    !
    ! Enhancement factors
    !
    gamma_s(aa:bb)   = cp_gga_c_param%hcth%gamma_parallel*s(aa:bb)
    gamma_s_2(aa:bb) = cp_gga_c_param%hcth%gamma_parallel*s_2(aa:bb)
    gamma_s(ab)      = cp_gga_c_param%hcth%gamma_opposite*s(ab)
    gamma_s_2(ab)    = cp_gga_c_param%hcth%gamma_opposite*s_2(ab)
    denominator(:)   = 1.0_real_8 + gamma_s_2(:)
    !
    u_hcth(1,:) = gamma_s_2(:) / denominator(:)
    u_hcth(2,:) = u_hcth(1,:) * u_hcth(1,:) 
    u_hcth(3,:) = u_hcth(1,:) * u_hcth(2,:) 
    u_hcth(4,:) = u_hcth(2,:) * u_hcth(2,:) 
    hcth_g(aa)  = cp_gga_c_param%hcth%c_0_parallel + sum(cp_gga_c_param%hcth%c_parallel(:) * u_hcth(:,aa))
    hcth_g(bb)  = cp_gga_c_param%hcth%c_0_parallel + sum(cp_gga_c_param%hcth%c_parallel(:) * u_hcth(:,bb))
    hcth_g(ab)  = cp_gga_c_param%hcth%c_0_opposite + sum(cp_gga_c_param%hcth%c_opposite(:) * u_hcth(:,ab))
    !
    ! Derivatives
    !
    dhcth_1_ds(:)   = 2.0_real_8*gamma_s(:) / denominator(:) * (1.0_real_8 - u_hcth(1,:)) 
    dhcth_u_ds(1,:) =            dhcth_1_ds(:)
    dhcth_u_ds(2,:) = 2.0_real_8*dhcth_1_ds(:)*u_hcth(1,:)
    dhcth_u_ds(3,:) = 3.0_real_8*dhcth_1_ds(:)*u_hcth(2,:)
    dhcth_u_ds(4,:) = 4.0_real_8*dhcth_1_ds(:)*u_hcth(3,:)
    !
    dhcth_g_ds(aa)  = sum(cp_gga_c_param%hcth%c_parallel(:)*dhcth_u_ds(:,aa))
    dhcth_g_ds(bb)  = sum(cp_gga_c_param%hcth%c_parallel(:)*dhcth_u_ds(:,bb))
    dhcth_g_ds(ab)  = sum(cp_gga_c_param%hcth%c_opposite(:)*dhcth_u_ds(:,ab))
    !
    sc(:)           = lda_sc(:)*hcth_g(:)
    !
    dsc_dna(:)      = lda_v1a(:)*hcth_g(:) + lda_sc(:)*dhcth_g_ds(:)*ds_dna(:)
    dsc_dnb(:)      = lda_v1b(:)*hcth_g(:) + lda_sc(:)*dhcth_g_ds(:)*ds_dnb(:)
    !
    dsc_dga(:)      = lda_sc(:)*dhcth_g_ds(:)*ds_dga(:)
    dsc_dgb(:)      = lda_sc(:)*dhcth_g_ds(:)*ds_dgb(:)
    !
    functional(aa)%sc  = sc(aa)
    functional(bb)%sc  = sc(bb)
    functional(ab)%sc  = sc(ab)
    !
    functional(aa)%v1c = sum(dsc_dna(:))
    functional(bb)%v1c = sum(dsc_dnb(:))
    !
    functional(aa)%v2c = sum(dsc_dga(:)) / scratch(aa)%abs_g
    functional(bb)%v2c = sum(dsc_dgb(:)) / scratch(bb)%abs_g
    functional(ab)%v2c = 0.0_real_8

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_spin_gga_c_hcth
  ! ==================================================================
  PURE SUBROUTINE cp_gga_c_hcth(scratch,functional)
    ! ==--------------------------------------------------------------==
    ! HCTH correlation rewritten from scratch, simplification of the
    ! SPIN routine
    ! o is used to represent lowercase sigma (o = a = b)
    !                                 24.06.2017 M.P.Bircher @ LCBC/EPFL
    ! ==--------------------------------------------------------------==


    TYPE(cp_xc_scratch_t), INTENT(in)        :: scratch
    TYPE(cp_xc_functional_t), INTENT(inout)  :: functional

    INTEGER, PARAMETER                       :: oo = 1
    INTEGER, PARAMETER                       :: ab = 2
    INTEGER, PARAMETER                       :: gm = 4
    REAL(real_8), PARAMETER                  :: lsda_to_lda       =  2.0_real_8
    REAL(real_8), PARAMETER                  :: sqrt_small        =  1.0e-12_real_8

    REAL(real_8)                             :: abs_g, s, s_2
    REAL(real_8), &
      DIMENSION(cp_xc_spin_pairs)            :: lda_sc, lda_v1, &
                                                sc, dsc_dno, dsc_dgo
    REAL(real_8), &
      DIMENSION(cp_xc_spin_pairs)            :: hcth_g, denominator, ds_dno, ds_dnb, &
                                                ds_dgo, ds_dgb, gamma_s, gamma_s_2, dhcth_1_ds, &
                                                dhcth_g_ds
    REAL(real_8), &
      DIMENSION(gm,cp_xc_spin_pairs)         :: u_hcth, dhcth_u_ds

    TYPE(cp_xc_scratch_t), &
      DIMENSION(cp_xc_spin_components)       :: parallel
    TYPE(cp_xc_functional_t), &
      DIMENSION(cp_xc_spin_components)       :: spin_functional

    !
    ! Parallel correlation
    !
    parallel(alpha)%n          = scratch%n
    parallel(alpha)%n_1_3      = scratch%n_1_3
    parallel(alpha)%n_4_3      = scratch%n_4_3
    parallel(beta )%n          = 0.0_real_8 
    parallel(beta )%n_1_3      = 0.0_real_8 
    parallel(beta )%n_4_3      = 0.0_real_8 
    parallel(alpha_beta)%n     = scratch%n
    parallel(alpha_beta)%n_1_3 = scratch%n_1_3
    parallel(alpha_beta)%n_4_3 = scratch%n_4_3
    !
    CALL cp_spin_lda_c_pw( parallel, spin_functional )
    lda_sc(oo) = spin_functional(alpha_beta)%sc
    lda_v1(oo) = spin_functional(alpha)%v1c
    !
    ! Antiparallel correlation
    !
    ! Eq. according to Stoll
    !
    CALL cp_lda_c_pw( scratch, functional )
    lda_sc(ab) = functional%sc  - lsda_to_lda*lda_sc(oo)
    lda_v1(ab) = functional%v1c - lda_v1(oo)

    !
    ! HCTH preliminaries (s(ab) = 0.5( s(aa) + s(bb) ))
    !
    abs_g      = max( scratch%abs_g, sqrt_small )
    s          = abs_g / scratch%n_4_3
    s_2        = s*s 
    !
    ds_dno(oo) = -four_thirds*s / scratch%n
    ds_dno(ab) = 0.5_real_8 * ds_dno(oo)  ! cf. above for the 0.5
    !
    ds_dgo(oo) = 1.0_real_8 / scratch%n_4_3 
    ds_dgo(ab) = 0.5_real_8 * ds_dgo(oo)  ! cf. above for the 0.5
    !
    ! Enhancement factors
    !
    gamma_s(oo)      = cp_gga_c_param%hcth%gamma_parallel*s
    gamma_s(ab)      = cp_gga_c_param%hcth%gamma_opposite*s
    gamma_s_2(oo)    = cp_gga_c_param%hcth%gamma_parallel*s_2
    gamma_s_2(ab)    = cp_gga_c_param%hcth%gamma_opposite*s_2
    denominator(:)   = 1.0_real_8 + gamma_s_2(:)
    !
    u_hcth(1,:) = gamma_s_2(:) / denominator(:)
    u_hcth(2,:) = u_hcth(1,:) * u_hcth(1,:) 
    u_hcth(3,:) = u_hcth(1,:) * u_hcth(2,:) 
    u_hcth(4,:) = u_hcth(2,:) * u_hcth(2,:) 
    hcth_g(oo)  = cp_gga_c_param%hcth%c_0_parallel + sum(cp_gga_c_param%hcth%c_parallel(:) * u_hcth(:,oo))
    hcth_g(ab)  = cp_gga_c_param%hcth%c_0_opposite + sum(cp_gga_c_param%hcth%c_opposite(:) * u_hcth(:,ab))
    !
    ! Derivatives
    !
    dhcth_1_ds(:)   = 2.0_real_8*gamma_s(:) / denominator(:) * (1.0_real_8 - u_hcth(1,:)) 
    dhcth_u_ds(1,:) =            dhcth_1_ds(:)
    dhcth_u_ds(2,:) = 2.0_real_8*dhcth_1_ds(:)*u_hcth(1,:)
    dhcth_u_ds(3,:) = 3.0_real_8*dhcth_1_ds(:)*u_hcth(2,:)
    dhcth_u_ds(4,:) = 4.0_real_8*dhcth_1_ds(:)*u_hcth(3,:)
    !
    dhcth_g_ds(oo) = sum(cp_gga_c_param%hcth%c_parallel(:)*dhcth_u_ds(:,oo))
    dhcth_g_ds(ab) = sum(cp_gga_c_param%hcth%c_opposite(:)*dhcth_u_ds(:,ab))
    !
    sc(oo)         = lsda_to_lda*lda_sc(oo)*hcth_g(oo)
    sc(ab)         =             lda_sc(ab)*hcth_g(ab)
    dsc_dno(:)     = lda_v1(:)*hcth_g(:) + lda_sc(:)*dhcth_g_ds(:)*ds_dno(:)
    dsc_dgo(:)     = lda_sc(:)*dhcth_g_ds(:)*ds_dgo(:)
    !
    functional%sc  = sum(sc(:))
    functional%v1c = sum(dsc_dno(:))
    functional%v2c = sum(dsc_dgo(:)) / (lsda_to_lda * scratch%abs_g)

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_gga_c_hcth
  ! ==================================================================

  ! ==================================================================
  ! Recurring functions
  ! ==================================================================

  ! ==================================================================
  ELEMENTAL FUNCTION phi_of(zeta) &
  RESULT  (phi)
    ! ==--------------------------------------------------------------==
    ! Phi function for PBE
    !                                            M.P.Bircher @ LCBC/EPFL
    ! ==--------------------------------------------------------------==

    REAL(real_8), INTENT(in)                 :: zeta
    REAL(real_8)                             :: phi

    phi = 0.5_real_8*((1._real_8 + zeta)**two_thirds + (1._real_8 - zeta)**two_thirds)

    ! ==--------------------------------------------------------------==
  END FUNCTION phi_of
  ! ==================================================================

  ! ==================================================================
  ! Check parameter compatibility
  ! ==================================================================

  ! ==================================================================
  ELEMENTAL FUNCTION cp_gga_c_check(tag) &
  RESULT  (OK)
    ! ==--------------------------------------------------------------==
    ! Important for pbe_flex and HCTH to prevent meaningless results,
    ! purely cosmetic for the other functionals
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
    CASE( "CP_GGA_C_HCTH","CP_SPIN_GGA_C_HCTH" )
      OK = cp_gga_c_param%init
    !
    ! Functionals with flexible parameters and dependencies on other
    ! routines
    !
    CASE( "CP_GGA_C_PBE_FLEX","CP_SPIN_GGA_C_PBE_FLEX" )
      OK = ( cp_gga_c_param%init .and. &
             cp_lda_c_check(cp_gga_c_param%pbe%lda_c) )
    !
    ! Functionals with hard-coded parameters
    !
    CASE( "CP_GGA_C_LYP", "CP_SPIN_GGA_C_LYP",&
          "CP_GGA_C_P86", "CP_SPIN_GGA_C_P86",&
          "CP_GGA_C_PBE", "CP_SPIN_GGA_C_PBE",&
          "CP_GGA_C_PBE_SOL", "CP_SPIN_GGA_C_PBE_SOL" )
      OK = .true.
    !
    ! Any other tag
    !
    CASE DEFAULT
      OK = .false.
    END SELECT

    ! ==--------------------------------------------------------------==
  END FUNCTION cp_gga_c_check
  ! ==================================================================
END MODULE cp_gga_correlation_utils
