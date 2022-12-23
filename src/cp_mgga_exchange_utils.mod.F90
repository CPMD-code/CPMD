! ==================================================================
! Provides: - Meta-GGA exchange, rewritten or adapted from OLDCODE
!
! Input:  Spin-polarised density (1/2 n for closed shell systems)
! Output: For spin-component only
! Convention: GGA conventions do not apply, since CAM is not
!             implemented for meta-functionals 
!
!                             18.04.2017 - M. P. Bircher @ LCBC/EPFL
! ==================================================================
MODULE cp_mgga_exchange_utils
  USE cnst,                            ONLY: pi
  USE cpfunc_types,                    ONLY: cp_xc_functional_t,&
                                             cp_xc_scratch_t,&
                                             cp_xc_spin_components
  USE func,                            ONLY: func2
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER :: a                           = 1
  INTEGER, PARAMETER :: b                           = 2
  INTEGER, PARAMETER :: ab                          = 3

  REAL(real_8), PARAMETER :: one_third              = 1.0_real_8/3.0_real_8
  REAL(real_8), PARAMETER :: two_thirds             = 2.0_real_8/3.0_real_8
  REAL(real_8), PARAMETER :: four_thirds            = 4.0_real_8/3.0_real_8
  REAL(real_8), PARAMETER :: five_thirds            = 5.0_real_8/3.0_real_8
  REAL(real_8), PARAMETER :: eight_thirds           = 8.0_real_8/3.0_real_8
  REAL(real_8), PARAMETER :: three_halfs            = 3.0_real_8/2.0_real_8
  REAL(real_8), PARAMETER :: four_ninths            = 4.0_real_8/9.0_real_8 
  REAL(real_8), PARAMETER :: three_fifths           = 3.0_real_8/5.0_real_8
  REAL(real_8), PARAMETER :: two_to_one_third       = 2.0_real_8**(one_third)
  REAL(real_8), PARAMETER :: two_to_four_thirds     = 2.0_real_8**(four_thirds)

  !
  ! Customisable functionals (FLEX type)
  ! 
  TYPE, PRIVATE :: vs98_x_t
     REAL(real_8)                :: r1 = 0.0_real_8 
     REAL(real_8)                :: r2 = 0.0_real_8
     REAL(real_8)                :: r3 = 0.0_real_8
     REAL(real_8)                :: r4 = 0.0_real_8
     REAL(real_8)                :: r5 = 0.0_real_8
     REAL(real_8)                :: r6 = 0.0_real_8
     REAL(real_8)                :: r7 = 0.0_real_8
  END TYPE vs98_x_t
  !
  TYPE, PRIVATE :: m05_m06_x_t
     REAL(real_8)                :: at00 = 0.0_real_8 
     REAL(real_8)                :: at01 = 0.0_real_8
     REAL(real_8)                :: at02 = 0.0_real_8
     REAL(real_8)                :: at03 = 0.0_real_8
     REAL(real_8)                :: at04 = 0.0_real_8
     REAL(real_8)                :: at05 = 0.0_real_8
     REAL(real_8)                :: at06 = 0.0_real_8
     REAL(real_8)                :: at07 = 0.0_real_8
     REAL(real_8)                :: at08 = 0.0_real_8
     REAL(real_8)                :: at09 = 0.0_real_8
     REAL(real_8)                :: at10 = 0.0_real_8
     REAL(real_8)                :: at11 = 0.0_real_8
     LOGICAL                     :: add_vs98 = .false.
  END TYPE m05_m06_x_t
  !
  TYPE, PRIVATE :: m08_m11_x_t
     REAL(real_8)                :: at00 = 0.0_real_8 
     REAL(real_8)                :: at01 = 0.0_real_8 
     REAL(real_8)                :: at02 = 0.0_real_8 
     REAL(real_8)                :: at03 = 0.0_real_8 
     REAL(real_8)                :: at04 = 0.0_real_8 
     REAL(real_8)                :: at05 = 0.0_real_8 
     REAL(real_8)                :: at06 = 0.0_real_8 
     REAL(real_8)                :: at07 = 0.0_real_8 
     REAL(real_8)                :: at08 = 0.0_real_8 
     REAL(real_8)                :: at09 = 0.0_real_8 
     REAL(real_8)                :: at10 = 0.0_real_8 
     REAL(real_8)                :: at11 = 0.0_real_8 
     REAL(real_8)                :: bt00 = 0.0_real_8 
     REAL(real_8)                :: bt01 = 0.0_real_8 
     REAL(real_8)                :: bt02 = 0.0_real_8 
     REAL(real_8)                :: bt03 = 0.0_real_8 
     REAL(real_8)                :: bt04 = 0.0_real_8 
     REAL(real_8)                :: bt05 = 0.0_real_8 
     REAL(real_8)                :: bt06 = 0.0_real_8 
     REAL(real_8)                :: bt07 = 0.0_real_8 
     REAL(real_8)                :: bt08 = 0.0_real_8 
     REAL(real_8)                :: bt09 = 0.0_real_8 
     REAL(real_8)                :: bt10 = 0.0_real_8 
     REAL(real_8)                :: bt11 = 0.0_real_8 
     REAL(real_8)                :: ct00 = 0.0_real_8 
     REAL(real_8)                :: ct01 = 0.0_real_8 
     REAL(real_8)                :: ct02 = 0.0_real_8 
     REAL(real_8)                :: ct03 = 0.0_real_8 
     REAL(real_8)                :: ct04 = 0.0_real_8 
     REAL(real_8)                :: ct05 = 0.0_real_8 
     REAL(real_8)                :: ct06 = 0.0_real_8 
     REAL(real_8)                :: ct07 = 0.0_real_8 
     REAL(real_8)                :: ct08 = 0.0_real_8 
     REAL(real_8)                :: ct09 = 0.0_real_8 
     REAL(real_8)                :: ct10 = 0.0_real_8 
     REAL(real_8)                :: ct11 = 0.0_real_8 
     REAL(real_8)                :: dt00 = 0.0_real_8 
     REAL(real_8)                :: dt01 = 0.0_real_8 
     REAL(real_8)                :: dt02 = 0.0_real_8 
     REAL(real_8)                :: dt03 = 0.0_real_8 
     REAL(real_8)                :: dt04 = 0.0_real_8 
     REAL(real_8)                :: dt05 = 0.0_real_8 
     REAL(real_8)                :: dt06 = 0.0_real_8 
     REAL(real_8)                :: dt07 = 0.0_real_8 
     REAL(real_8)                :: dt08 = 0.0_real_8 
     REAL(real_8)                :: dt09 = 0.0_real_8 
     REAL(real_8)                :: dt10 = 0.0_real_8 
     REAL(real_8)                :: dt11 = 0.0_real_8 
     LOGICAL                     :: lda_has_lc = .false.
  END TYPE m08_m11_x_t
  !
  TYPE, PRIVATE :: mn12_mn15_x_t
     REAL(real_8)                :: cc000 = 0.0_real_8 
     REAL(real_8)                :: cc001 = 0.0_real_8 
     REAL(real_8)                :: cc002 = 0.0_real_8 
     REAL(real_8)                :: cc003 = 0.0_real_8 
     REAL(real_8)                :: cc004 = 0.0_real_8 
     REAL(real_8)                :: cc005 = 0.0_real_8 
     REAL(real_8)                :: cc010 = 0.0_real_8 
     REAL(real_8)                :: cc011 = 0.0_real_8 
     REAL(real_8)                :: cc012 = 0.0_real_8 
     REAL(real_8)                :: cc013 = 0.0_real_8 
     REAL(real_8)                :: cc014 = 0.0_real_8 
     REAL(real_8)                :: cc020 = 0.0_real_8 
     REAL(real_8)                :: cc021 = 0.0_real_8 
     REAL(real_8)                :: cc022 = 0.0_real_8 
     REAL(real_8)                :: cc023 = 0.0_real_8 
     REAL(real_8)                :: cc030 = 0.0_real_8 
     REAL(real_8)                :: cc031 = 0.0_real_8 
     REAL(real_8)                :: cc032 = 0.0_real_8 
     REAL(real_8)                :: cc100 = 0.0_real_8 
     REAL(real_8)                :: cc101 = 0.0_real_8 
     REAL(real_8)                :: cc102 = 0.0_real_8 
     REAL(real_8)                :: cc103 = 0.0_real_8 
     REAL(real_8)                :: cc104 = 0.0_real_8 
     REAL(real_8)                :: cc110 = 0.0_real_8 
     REAL(real_8)                :: cc111 = 0.0_real_8 
     REAL(real_8)                :: cc112 = 0.0_real_8 
     REAL(real_8)                :: cc113 = 0.0_real_8 
     REAL(real_8)                :: cc120 = 0.0_real_8 
     REAL(real_8)                :: cc121 = 0.0_real_8 
     REAL(real_8)                :: cc122 = 0.0_real_8 
     REAL(real_8)                :: cc200 = 0.0_real_8 
     REAL(real_8)                :: cc201 = 0.0_real_8 
     REAL(real_8)                :: cc202 = 0.0_real_8 
     REAL(real_8)                :: cc203 = 0.0_real_8 
     REAL(real_8)                :: cc210 = 0.0_real_8 
     REAL(real_8)                :: cc211 = 0.0_real_8 
     REAL(real_8)                :: cc212 = 0.0_real_8 
     REAL(real_8)                :: cc300 = 0.0_real_8 
     REAL(real_8)                :: cc301 = 0.0_real_8 
     REAL(real_8)                :: cc302 = 0.0_real_8 
  END TYPE mn12_mn15_x_t
  !
  TYPE, PRIVATE :: cp_mgga_x_param_t
     LOGICAL                     :: init = .false.
     TYPE(vs98_x_t)              :: VS98
     TYPE(m05_m06_x_t)           :: M05_M06
     TYPE(m08_m11_x_t)           :: M08_M11
     TYPE(mn12_mn15_x_t)         :: MN12_MN15
  END TYPE cp_mgga_x_param_t
  !
  TYPE(cp_mgga_x_param_t), PUBLIC, SAVE :: cp_mgga_x_param

  PUBLIC :: CP_MGGA_X_TPSS
  PUBLIC :: CP_MGGA_X_M05_M06
  PUBLIC :: CP_MGGA_X_M08_M11
  PUBLIC :: CP_MGGA_X_MN12_MN15
  PUBLIC :: CP_SPIN_MGGA_X_TPSS
  PUBLIC :: CP_SPIN_MGGA_X_M05_M06
  PUBLIC :: CP_SPIN_MGGA_X_M08_M11
  PUBLIC :: CP_SPIN_MGGA_X_MN12_MN15

  PUBLIC :: CP_MGGA_X_CHECK

  PUBLIC :: CP_MGGA_X_VS98_GVT4

  PRIVATE :: CP_MGGA_X_VS98
  PRIVATE :: CP_SPIN_MGGA_X_VS98

CONTAINS

  ! ==================================================================
  PURE SUBROUTINE cp_spin_mgga_x_tpss( scratch,functional )
    ! ==--------------------------------------------------------------==


    TYPE(cp_xc_scratch_t), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(in)                             :: scratch
    TYPE(cp_xc_functional_t), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(inout)                          :: functional

    IF (scratch(a)%n >= 1.0e-10_real_8) &
         CALL cp_mgga_x_tpss( scratch(a), functional(a) )
    IF (scratch(b)%n >= 1.0e-10_real_8) &
         CALL cp_mgga_x_tpss( scratch(b), functional(b) )

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_spin_mgga_x_tpss
  ! ==================================================================
  PURE SUBROUTINE cp_mgga_x_tpss( scratch, functional )
    ! ==--------------------------------------------------------------==
    !
    ! Adapted from OLDCODE
    !                                     18.04.2017 M.P. Bircher @ EPFL
    ! ==--------------------------------------------------------------==

    TYPE(cp_xc_scratch_t), INTENT(in)        :: scratch
    TYPE(cp_xc_functional_t), INTENT(inout)  :: functional

    REAL(real_8), PARAMETER :: b = 0.40_real_8, c = 1.59096_real_8, &
      cexunif = -3._real_8/(4._real_8*pi) * (3._real_8*pi*pi)**one_third, &
      cp = 4._real_8*(3._real_8*pi*pi)**two_thirds, e = 1.537_real_8, &
      f1081 = 10._real_8/81._real_8, kappa = 0.804_real_8, &
      mu = 0.21951_real_8, se = SQRT(e), small = 1.0e-14_real_8 

    REAL(real_8) :: alpha, dfxdg, dfxdr, dfxdt, dfxdx, dpdg, dpdr, dqbtda, &
      dqbtdp, dqbtdz, dxdp, dxdqbt, dxdz, dzdg, dzdr, dzdt, emup_2, exunif, &
      fx, grho, opsep_2, p, p_2, qbt, rho, rho_1_3, rho_8_3, sep, sx, tau, &
      v1x, v2x, vtt, x, z, z_2

    rho     = 2.0_real_8 * scratch%n
    rho_1_3 = two_to_one_third * scratch%n_1_3
    grho    = 4.0_real_8 * scratch%g_2
    tau     = 2.0_real_8 * scratch%tau

    exunif  = cexunif * rho_1_3

    IF ( ABS(tau) > small .AND. grho > small ) THEN
       rho_8_3 = two_to_four_thirds * scratch%n_4_3
       rho_8_3 = rho_8_3 * rho_8_3
       p       = grho/(cp*rho_8_3)
       p_2     = p*p
       emup_2  = e*mu*p_2
       sep     = se*p
       opsep_2 = (1.0_real_8 + sep)**2
       z       = 0.125_real_8*grho/rho/tau
       z_2     = z*z
       alpha   = 5._real_8*one_third*p * (1._real_8/z - 1._real_8)
       qbt     = 0.45_real_8*(alpha - 1._real_8)/SQRT(1._real_8 + b*alpha*(alpha - 1._real_8)) &
            + two_thirds*p
       x       = ( f1081 + c*z_2/(1._real_8 + z_2)**2 ) * p + 146._real_8/2025._real_8*qbt**2&
            - 73._real_8/405._real_8 * qbt * SQRT ( 0.5_real_8*(0.6_real_8*z)**2 + 0.5_real_8*p_2)&
            + f1081**2*p_2/kappa + 2._real_8*se*f1081*(0.6_real_8*z)**2 + emup_2*p
       x       = x/opsep_2
       fx      = 1._real_8 + kappa - kappa/(1._real_8 + x/kappa)
       sx      = rho*exunif*fx

       dfxdx  = 1._real_8/(1._real_8 + x/kappa)**2
       dqbtda = 0.225_real_8*(2._real_8 + b*alpha - b) &
            /(1._real_8 + b*alpha*alpha - b*alpha)**1.5_real_8
       dqbtdp = dqbtda*5._real_8*one_third*(1._real_8/z-1._real_8) + two_thirds
       dqbtdz = dqbtda*5._real_8*one_third*p*(-1._real_8/(z_2))
       dxdqbt = (146._real_8/2025._real_8 * 2._real_8*qbt &
            - 73._real_8/405._real_8 * SQRT ( 0.18_real_8*z_2 + 0.5_real_8*p_2 )) &
            /opsep_2
       dxdz   = dxdqbt*dqbtdz + ( c*p*(2._real_8*z/(1._real_8+z_2)**2 &
            - 4._real_8*z**3/(1._real_8+z_2)**3 ) - 73._real_8/405._real_8 * qbt*0.5_real_8 &
            * 0.36_real_8*z / SQRT ( 0.18_real_8*z_2 + 0.5_real_8*p_2 ) &
            + 4._real_8*se*f1081*0.36_real_8*z ) / opsep_2
       dxdp   = -2._real_8*se*x/(1._real_8+sep) + dxdqbt*dqbtdp &
            + ( ( f1081 + c*z_2/(1._real_8+z_2)**2 ) - 73._real_8/405._real_8*qbt &
            * 0.5_real_8*p/SQRT ( 0.5_real_8*(0.6_real_8*z)**2 + 0.5_real_8*p_2 ) &
            + 2._real_8*f1081**2*p/kappa + 3._real_8*emup_2 )/opsep_2
       dpdr   = -eight_thirds*p/rho
       dpdg   = 2._real_8*p/SQRT(grho)
       dzdr   = -z/rho
       dzdg   = 2._real_8*z/SQRT(grho)
       dzdt   = -z/tau
       dfxdr  = dfxdx*(dxdp*dpdr + dxdz*dzdr)
       dfxdg  = dfxdx*(dxdp*dpdg + dxdz*dzdg)
       dfxdt  = dfxdx*dxdz*dzdt
       ! 
       v1x    = exunif*(fx + one_third*fx + rho*dfxdr)
       v2x    = rho*exunif*dfxdg
       vtt    = rho*exunif*dfxdt
    ELSE
       fx     = 1.0_real_8
       sx     = rho*exunif*fx
       v1x    = exunif + one_third*sx/rho
       v2x    = 0.0_real_8
       vtt    = 0.0_real_8
    ENDIF

    !
    ! K, dK_dn, dK_dg are not needed, since the CAM routine cannot
    ! deal with an MGGA, anyway - so we zero them
    ! (Maybe the above part can be rewritten in a later release, such that the terms are there)
    !
    functional%K        = 0.0_real_8 ! if ever needed, one could convert from F[x] -> K[x] here, too 
    functional%dK_dn    = 0.0_real_8 ! if ever needed, one could convert from F[x] -> K[x] here, too
    functional%dK_dg    = 0.0_real_8 ! if ever needed, one could convert from F[x] -> K[x] here, too
    !
    functional%sx_sigma = 0.5_real_8*sx
    functional%dsx_dn   = v1x
    functional%dsx_dg   = v2x
    functional%dsx_dt   = vtt

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_mgga_x_tpss
  ! ==================================================================
  PURE SUBROUTINE cp_spin_mgga_x_m05_m06( scratch,functional )
    ! ==--------------------------------------------------------------==
    ! Generic spin-unrestricted interface to Minnesota 05/06
    !                                     24.04.2017 M.P. Bircher @ EPFL

    TYPE(cp_xc_scratch_t), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(in)                             :: scratch
    TYPE(cp_xc_functional_t), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(inout)                          :: functional

    IF (scratch(a)%n >= 1.0e-10_real_8) &
         CALL cp_mgga_x_m05_m06( scratch(a), functional(a) )
    IF (scratch(b)%n >= 1.0e-10_real_8) &
         CALL cp_mgga_x_m05_m06( scratch(b), functional(b) )

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_spin_mgga_x_m05_m06
  ! ==================================================================
  PURE SUBROUTINE cp_mgga_x_m05_m06( scratch,functional )
    ! ==--------------------------------------------------------------==
    ! Combined M05 and M06 into one single routine
    ! Cleaned up, prettified
    ! Adapted from OLDCODE to xc_driver format
    !                                     24.04.2017 M.P. Bircher @ EPFL
    ! Porting based on MFM 1.9 at http://comp.chem.umn.edu/mfm/
    !                                        2015 P. Lopez-Tarifa @ EPFL
    ! ==--------------------------------------------------------------==

    TYPE(cp_xc_scratch_t), INTENT(in)        :: scratch
    TYPE(cp_xc_functional_t), INTENT(inout)  :: functional

    REAL(real_8), PARAMETER                  :: c1     =  3.36116E-03_real_8
    REAL(real_8), PARAMETER                  :: c2     =  4.49267e-03_real_8
    REAL(real_8), PARAMETER                  :: ax     = -1.5_real_8*(four_thirds*pi)**(-one_third)
    REAL(real_8), PARAMETER                  :: ctau   =  ((6.0_real_8*pi*pi)**two_thirds) 
    REAL(real_8), PARAMETER                  :: ctaueg =  three_fifths*ctau 
    REAL(real_8), PARAMETER                  :: small  =  1.0e-08_real_8

    REAL(real_8)                             :: tau, n_5_3, taueg, tsig, wsig, fsig
    REAL(real_8)                             :: sx, dsx_dn, dsx_dg, dsx_dt
    REAL(real_8)                             :: x, x2, en, ed, e
    REAL(real_8)                             :: den, ded, de, dfdw, dwdt, dtdr, dtdtau
    REAL(real_8)                             :: dggadr, dfdr, dfdtau, dggadg

    IF (cp_mgga_x_param%M05_M06%add_vs98) CALL cp_mgga_x_vs98( scratch, functional )

    sx     = functional%sx_sigma
    dsx_dn = functional%dsx_dn
    dsx_dg = functional%dsx_dg
    dsx_dt = functional%dsx_dt

    tau = 2.0_real_8*scratch%tau

    if ( scratch%n > small ) then 
       n_5_3  = scratch%n_1_3*scratch%n_4_3
       taueg  = ctaueg*n_5_3
       tsig   = taueg/tau
       wsig   = (tsig-1.0_real_8)/(tsig+1.0_real_8)
       fsig   = (cp_mgga_x_param%M05_M06%at00+wsig*(cp_mgga_x_param%M05_M06%at01+wsig &
                *(cp_mgga_x_param%M05_M06%at02+wsig*(cp_mgga_x_param%M05_M06%at03+wsig &
                 *(cp_mgga_x_param%M05_M06%at04+wsig*(cp_mgga_x_param%M05_M06%at05+wsig &
                  *(cp_mgga_x_param%M05_M06%at06+wsig*(cp_mgga_x_param%M05_M06%at07+wsig &
                   *(cp_mgga_x_param%M05_M06%at08+wsig*(cp_mgga_x_param%M05_M06%at09+wsig &
                    *(cp_mgga_x_param%M05_M06%at10+wsig*cp_mgga_x_param%M05_M06%at11)))))))))))
       x      = scratch%abs_g/scratch%n_4_3
       x2     = x*x
       en     = c1*x2
       ed     = 1.0_real_8+c2*x2
       e      = -en/ed
       sx     = sx + (ax+e)*fsig*scratch%n_4_3

       den    = 2.0_real_8*c1*x
       ded    = 2.0_real_8*c2*x
       de     = -(den*ed-en*ded)/(ed*ed)
       dfdw   = (cp_mgga_x_param%M05_M06%at01 + wsig*(2.0_real_8*cp_mgga_x_param%M05_M06%at02 + wsig &
                *(3.0_real_8*cp_mgga_x_param%M05_M06%at03 + wsig &
                 *(4.0_real_8*cp_mgga_x_param%M05_M06%at04 + wsig*(5.0_real_8*cp_mgga_x_param%M05_M06%at05 + wsig &
                  *(6.0_real_8*cp_mgga_x_param%M05_M06%at06 + wsig &
                   *(7.0_real_8*cp_mgga_x_param%M05_M06%at07 + wsig*(8.0_real_8*cp_mgga_x_param%M05_M06%at08 + wsig &
                    *(9.0_real_8*cp_mgga_x_param%M05_M06%at09 + wsig &
                     *(10.0_real_8*cp_mgga_x_param%M05_M06%at10 + wsig*11.0_real_8*cp_mgga_x_param%M05_M06%at11))))))))))
       dwdt   = 2.0_real_8/((1.0_real_8 + tsig)**2)
       dtdr   = ctau*(scratch%n_1_3 * scratch%n_1_3)/tau
       dtdtau = -taueg/(tau*tau)
       dggadr = four_thirds*scratch%n_1_3*(ax + (e-x*de))
       dfdr   = dfdw*dwdt*dtdr
       dfdtau = dfdw*dwdt*dtdtau
       dggadg = de
       dsx_dn = dsx_dn + dggadr*fsig + (ax + e)*scratch%n_4_3*dfdr
       dsx_dg = dsx_dg + dggadg*fsig
       dsx_dt = dsx_dt + 2.0_real_8*scratch%n_4_3*(ax + e)*dfdtau
    end if

    functional%sx_sigma = sx
    functional%dsx_dn   = dsx_dn 
    functional%dsx_dg   = dsx_dg
    functional%dsx_dt   = dsx_dt
    !
    functional%K        = 0.0_real_8 ! if ever needed, one could convert from F[x] -> K[x] here
    functional%dK_dn    = 0.0_real_8 ! if ever needed, one could convert from F[x] -> K[x] here
    functional%dK_dg    = 0.0_real_8 ! if ever needed, one could convert from F[x] -> K[x] here

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_mgga_x_m05_m06
  ! ==================================================================
  PURE SUBROUTINE cp_spin_mgga_x_m08_m11( scratch,functional )
    ! ==--------------------------------------------------------------==
    ! Generic spin-unrestricted interface to Minnesota 08/11
    !                                     08.05.2017 M.P. Bircher @ EPFL

    TYPE(cp_xc_scratch_t), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(in)                             :: scratch
    TYPE(cp_xc_functional_t), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(inout)                          :: functional

    IF (scratch(a)%n >= 1.0e-10_real_8) &
         CALL cp_mgga_x_m08_m11( scratch(a), functional(a) )
    IF (scratch(b)%n >= 1.0e-10_real_8) &
         CALL cp_mgga_x_m08_m11( scratch(b), functional(b) )

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_spin_mgga_x_m08_m11
  ! ==================================================================
  PURE SUBROUTINE cp_mgga_x_m08_m11( scratch, functional )
    ! ==--------------------------------------------------------------==
    ! Cleaned up, prettified, simplified
    ! Adapted from OLDCODE to xc_driver format
    !                                     12.05.2017 M.P. Bircher @ EPFL
    ! Porting based on MFM 1.9 at http://comp.chem.umn.edu/mfm/
    !                                        2015 P. Lopez-Tarifa @ EPFL
    ! ==--------------------------------------------------------------==

    TYPE(cp_xc_scratch_t), INTENT(in)        :: scratch
    TYPE(cp_xc_functional_t), INTENT(inout)  :: functional

    REAL(real_8), PARAMETER                  :: ax         = -1.5_real_8*(four_thirds*pi)**(-one_third)
    REAL(real_8), PARAMETER                  :: mu         =  0.21951_real_8
    REAL(real_8), PARAMETER                  :: mus        = 10.0_real_8 / 81.0_real_8
    REAL(real_8), PARAMETER                  :: kappa      =  0.804_real_8
    REAL(real_8), PARAMETER                  :: kappas     =  0.552_real_8
    REAL(real_8), PARAMETER                  :: ctau       =  ((6.0_real_8*pi*pi)**two_thirds) 
    REAL(real_8), PARAMETER                  :: ctaueg     =  three_fifths*ctau 
    REAL(real_8), PARAMETER                  :: cs         =  (48.0_real_8*pi*pi)**one_third
    REAL(real_8), PARAMETER                  :: small_tau  =  1.00e-10_real_8
    REAL(real_8), PARAMETER                  :: small_g_2  =  0.25e-10_real_8

    REAL(real_8)                             :: tau, sx, dsx_dn, dsx_dg, dsx_dt
    REAL(real_8)                             :: tauueg, tsig, wsig, fsig1, fsig2, fsig3, fsig4
    REAL(real_8)                             :: w1, w2, w3, w4, w5, w6, w7, w8, w9, w10, w11
    REAL(real_8)                             :: x, y, s, deno, fx1, fx2
    REAL(real_8)                             :: ellr, elsr, delsrdn, dellrdn
    REAL(real_8)                             :: gga1, gga2, gga3, gga4, dydn, dydg
    REAL(real_8)                             :: dfx1dy, dfx1dn, dfx1dg
    REAL(real_8)                             :: dfx2dy, dfx2dn, dfx2dg
    REAL(real_8)                             :: df1dw, df2dw, df3dw, df4dw, dwdt, dtdn, dtdtau
    REAL(real_8)                             :: df1dn, df1dtau, df2dn, df2dtau, df3dn, df3dtau, df4dn, df4dtau
    REAL(real_8)                             :: dgga1dn, dgga2dn, dgga3dn, dgga4dn
    REAL(real_8)                             :: dgga1dg, dgga2dg, dgga3dg, dgga4dg


    tau  = 2.0_real_8 * scratch%tau

    if( abs(tau) > small_tau .and. scratch%g_2 > small_g_2 ) then
       tauueg = ctaueg * scratch%n_4_3 * scratch%n_1_3
       tsig   = tauueg / tau
       wsig   = (tsig-1.0_real_8)/(tsig+1.0_real_8)

       w1     = wsig
       w2     = wsig*w1
       w3     = wsig*w2
       w4     = wsig*w3
       w5     = wsig*w4
       w6     = wsig*w5
       w7     = wsig*w6
       w8     = wsig*w7
       w9     = wsig*w8
       w10    = wsig*w9
       w11    = wsig*w10
       fsig1  =    cp_mgga_x_param%M08_M11%at00 &
                 + cp_mgga_x_param%M08_M11%at01*w1 &
                 + cp_mgga_x_param%M08_M11%at02*w2 &
                 + cp_mgga_x_param%M08_M11%at03*w3 & 
                 + cp_mgga_x_param%M08_M11%at04*w4 &
                 + cp_mgga_x_param%M08_M11%at05*w5 &
                 + cp_mgga_x_param%M08_M11%at06*w6 &
                 + cp_mgga_x_param%M08_M11%at07*w7 & 
                 + cp_mgga_x_param%M08_M11%at08*w8 &
                 + cp_mgga_x_param%M08_M11%at09*w9 &
                 + cp_mgga_x_param%M08_M11%at10*w10 &
                 + cp_mgga_x_param%M08_M11%at11*w11
       fsig2  =    cp_mgga_x_param%M08_M11%bt00 &
                 + cp_mgga_x_param%M08_M11%bt01*w1 &
                 + cp_mgga_x_param%M08_M11%bt02*w2 &
                 + cp_mgga_x_param%M08_M11%bt03*w3 & 
                 + cp_mgga_x_param%M08_M11%bt04*w4 &
                 + cp_mgga_x_param%M08_M11%bt05*w5 &
                 + cp_mgga_x_param%M08_M11%bt06*w6 &
                 + cp_mgga_x_param%M08_M11%bt07*w7 & 
                 + cp_mgga_x_param%M08_M11%bt08*w8 &
                 + cp_mgga_x_param%M08_M11%bt09*w9 &
                 + cp_mgga_x_param%M08_M11%bt10*w10 &
                 + cp_mgga_x_param%M08_M11%bt11*w11
       fsig3  =    cp_mgga_x_param%M08_M11%ct00 &
                 + cp_mgga_x_param%M08_M11%ct01*w1 &
                 + cp_mgga_x_param%M08_M11%ct02*w2 &
                 + cp_mgga_x_param%M08_M11%ct03*w3 & 
                 + cp_mgga_x_param%M08_M11%ct04*w4 &
                 + cp_mgga_x_param%M08_M11%ct05*w5 &
                 + cp_mgga_x_param%M08_M11%ct06*w6 &
                 + cp_mgga_x_param%M08_M11%ct07*w7 & 
                 + cp_mgga_x_param%M08_M11%ct08*w8 &
                 + cp_mgga_x_param%M08_M11%ct09*w9 &
                 + cp_mgga_x_param%M08_M11%ct10*w10 &
                 + cp_mgga_x_param%M08_M11%ct11*w11
       fsig4  =    cp_mgga_x_param%M08_M11%dt00 &
                 + cp_mgga_x_param%M08_M11%dt01*w1 &
                 + cp_mgga_x_param%M08_M11%dt02*w2 &
                 + cp_mgga_x_param%M08_M11%dt03*w3 & 
                 + cp_mgga_x_param%M08_M11%dt04*w4 &
                 + cp_mgga_x_param%M08_M11%dt05*w5 &
                 + cp_mgga_x_param%M08_M11%dt06*w6 &
                 + cp_mgga_x_param%M08_M11%dt07*w7 & 
                 + cp_mgga_x_param%M08_M11%dt08*w8 &
                 + cp_mgga_x_param%M08_M11%dt09*w9 &
                 + cp_mgga_x_param%M08_M11%dt10*w10 &
                 + cp_mgga_x_param%M08_M11%dt11*w11

       x      = scratch%abs_g/scratch%n_4_3
       s      = x/cs
       y      = s*s
       deno   = (1.0_real_8 + mu*y/kappa)
       fx1    = 1.0_real_8 + kappa*(1.0_real_8 - 1.0_real_8/deno)
       fx2    = 1.0_real_8 + kappas*(1.0_real_8 - exp( - mus*y/kappas))

       if (cp_mgga_x_param%M08_M11%lda_has_lc) then 
          CALL cp_mgga_x_m08_m11_lc(scratch,elsr,delsrdn)
          ellr    = ax*scratch%n_4_3 - elsr
          dellrdn = ax*four_thirds*scratch%n_1_3 - delsrdn
       else 
          elsr    = ax*scratch%n_4_3
          ellr    = 0.0_real_8
          delsrdn = ax*four_thirds*scratch%n_1_3
          dellrdn = 0.0_real_8
       endif

       gga1 = elsr*fx1
       gga2 = elsr*fx2
       gga3 = ellr*fx1
       gga4 = ellr*fx2
       sx   = (gga1*fsig1 + gga2*fsig2 + gga3*fsig3 + gga4*fsig4)

       dydn   =  -(8.0_real_8/3.0_real_8)*y/scratch%n
       dydg   = 2.0_real_8*y/scratch%abs_g
       dfx1dy = mu/(deno*deno)
       dfx1dn = dfx1dy*dydn
       dfx1dg = dfx1dy*dydg
       dfx2dy = mus*exp(-mus*y/kappas)
       dfx2dn = dfx2dy*dydn
       dfx2dg = dfx2dy*dydg

       df1dw  =               cp_mgga_x_param%M08_M11%at01 &
                +  2.0_real_8*cp_mgga_x_param%M08_M11%at02*w1 &
                +  3.0_real_8*cp_mgga_x_param%M08_M11%at03*w2 &
                +  4.0_real_8*cp_mgga_x_param%M08_M11%at04*w3 & 
                +  5.0_real_8*cp_mgga_x_param%M08_M11%at05*w4 &
                +  6.0_real_8*cp_mgga_x_param%M08_M11%at06*w5 &
                +  7.0_real_8*cp_mgga_x_param%M08_M11%at07*w6 &
                +  8.0_real_8*cp_mgga_x_param%M08_M11%at08*w7 &
                +  9.0_real_8*cp_mgga_x_param%M08_M11%at09*w8 &
                + 10.0_real_8*cp_mgga_x_param%M08_M11%at10*w9 &
                + 11.0_real_8*cp_mgga_x_param%M08_M11%at11*w10
       df2dw  =               cp_mgga_x_param%M08_M11%bt01 &
                +  2.0_real_8*cp_mgga_x_param%M08_M11%bt02*w1 &
                +  3.0_real_8*cp_mgga_x_param%M08_M11%bt03*w2 &
                +  4.0_real_8*cp_mgga_x_param%M08_M11%bt04*w3 & 
                +  5.0_real_8*cp_mgga_x_param%M08_M11%bt05*w4 &
                +  6.0_real_8*cp_mgga_x_param%M08_M11%bt06*w5 &
                +  7.0_real_8*cp_mgga_x_param%M08_M11%bt07*w6 &
                +  8.0_real_8*cp_mgga_x_param%M08_M11%bt08*w7 &
                +  9.0_real_8*cp_mgga_x_param%M08_M11%bt09*w8 &
                + 10.0_real_8*cp_mgga_x_param%M08_M11%bt10*w9 &
                + 11.0_real_8*cp_mgga_x_param%M08_M11%bt11*w10
       df3dw  =               cp_mgga_x_param%M08_M11%ct01 &
                +  2.0_real_8*cp_mgga_x_param%M08_M11%ct02*w1 &
                +  3.0_real_8*cp_mgga_x_param%M08_M11%ct03*w2 &
                +  4.0_real_8*cp_mgga_x_param%M08_M11%ct04*w3 & 
                +  5.0_real_8*cp_mgga_x_param%M08_M11%ct05*w4 &
                +  6.0_real_8*cp_mgga_x_param%M08_M11%ct06*w5 &
                +  7.0_real_8*cp_mgga_x_param%M08_M11%ct07*w6 &
                +  8.0_real_8*cp_mgga_x_param%M08_M11%ct08*w7 &
                +  9.0_real_8*cp_mgga_x_param%M08_M11%ct09*w8 &
                + 10.0_real_8*cp_mgga_x_param%M08_M11%ct10*w9 &
                + 11.0_real_8*cp_mgga_x_param%M08_M11%ct11*w10
       df4dw  =               cp_mgga_x_param%M08_M11%dt01 &
                +  2.0_real_8*cp_mgga_x_param%M08_M11%dt02*w1 &
                +  3.0_real_8*cp_mgga_x_param%M08_M11%dt03*w2 &
                +  4.0_real_8*cp_mgga_x_param%M08_M11%dt04*w3 & 
                +  5.0_real_8*cp_mgga_x_param%M08_M11%dt05*w4 &
                +  6.0_real_8*cp_mgga_x_param%M08_M11%dt06*w5 &
                +  7.0_real_8*cp_mgga_x_param%M08_M11%dt07*w6 &
                +  8.0_real_8*cp_mgga_x_param%M08_M11%dt08*w7 &
                +  9.0_real_8*cp_mgga_x_param%M08_M11%dt09*w8 &
                + 10.0_real_8*cp_mgga_x_param%M08_M11%dt10*w9 &
                + 11.0_real_8*cp_mgga_x_param%M08_M11%dt11*w10
 
       dwdt   = 2.0_real_8/((1.0_real_8+tsig)*(1.0_real_8+tsig))                                          
       dtdn   = ctau*(scratch%n_1_3 * scratch%n_1_3)/tau                  
       dtdtau = -tauueg / (tau*tau)

       dgga1dn = delsrdn*fx1 + elsr*dfx1dn
       dgga2dn = delsrdn*fx2 + elsr*dfx2dn
       dgga3dn = dellrdn*fx1 + ellr*dfx1dn
       dgga4dn = dellrdn*fx2 + ellr*dfx2dn

       df1dn   = df1dw*dwdt*dtdn
       df2dn   = df2dw*dwdt*dtdn
       df3dn   = df3dw*dwdt*dtdn
       df4dn   = df4dw*dwdt*dtdn
       df1dtau = df1dw*dwdt*dtdtau
       df2dtau = df2dw*dwdt*dtdtau
       df3dtau = df3dw*dwdt*dtdtau
       df4dtau = df4dw*dwdt*dtdtau

       dgga1dg = elsr*dfx1dg
       dgga2dg = elsr*dfx2dg
       dgga3dg = ellr*dfx1dg
       dgga4dg = ellr*dfx2dg

       dsx_dn  =   dgga1dn*fsig1 + gga1*df1dn + dgga2dn*fsig2 + gga2*df2dn & 
                 + dgga3dn*fsig3 + gga3*df3dn + dgga4dn*fsig4 + gga4*df4dn
       dsx_dg  = dgga1dg*fsig1 + dgga2dg*fsig2 + dgga3dg*fsig3 + dgga4dg*fsig4
       dsx_dt  = 2.0_real_8*(gga1*df1dtau + gga2*df2dtau + gga3*df3dtau + gga4*df4dtau)
    else 
       sx      = ax*scratch%n_4_3
       dsx_dn  = four_thirds*ax*scratch%n_1_3
       dsx_dg  = 0.0_real_8 
       dsx_dt  = 0.0_real_8 
    end if 

    functional%sx_sigma = sx
    functional%dsx_dn   = dsx_dn 
    functional%dsx_dg   = dsx_dg
    functional%dsx_dt   = dsx_dt
    !
    functional%K        = 0.0_real_8 ! if ever needed, one could convert from F[x] -> K[x] here
    functional%dK_dn    = 0.0_real_8 ! if ever needed, one could convert from F[x] -> K[x] here
    functional%dK_dg    = 0.0_real_8 ! if ever needed, one could convert from F[x] -> K[x] here

    CONTAINS

    ! =--------------------------------------------------------------==
    ELEMENTAL SUBROUTINE cp_mgga_x_m08_m11_lc(scratch,f,d1f)

      TYPE(cp_xc_scratch_t), INTENT(in)        :: scratch
      REAL(real_8), INTENT(out)                :: f, d1f

      REAL(real_8), PARAMETER                  :: sqrt_pi = sqrt(pi)
      REAL(real_8), PARAMETER                  :: cmu = (6.0_real_8*pi*pi)**one_third
      REAL(real_8), PARAMETER                  :: emu = 0.25_real_8

      REAL(real_8)                             :: fsr, tmu, tmu2, tmu3, tmu4
      REAL(real_8)                             :: w, erfv, dtmudr, dfsrdtmu
 
      tmu      = func2%srxa/(2.0_real_8*cmu*scratch%n_1_3)
      tmu2     = tmu*tmu
      tmu3     = tmu*tmu2
      tmu4     = tmu*tmu3
      w        = exp(-0.25_real_8/tmu2)
      erfv     = erf( 0.50_real_8/tmu )
      dtmudr   = -one_third*tmu/scratch%n
      fsr      = 1.0_real_8 - four_thirds*tmu*(-6.0_real_8*tmu + 8.0_real_8*tmu3 &
                 + w*(4.0_real_8*tmu - 8.0_real_8*tmu3) + 2.0_real_8*sqrt_pi*erfv)
      dfsrdtmu = eight_thirds*(2.0_real_8*tmu*(3.0_real_8 - 8.0_real_8*tmu2 &
                 + w*(-1.0_real_8 + 8.0_real_8*tmu2)) - sqrt_pi*erfv)
      f        = ax*scratch%n_4_3*fsr
      d1f      = four_thirds*ax*scratch%n_1_3*fsr &
                 + ax*scratch%n_4_3*(dfsrdtmu*dtmudr)

    END SUBROUTINE cp_mgga_x_m08_m11_lc
    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_mgga_x_m08_m11
  ! ==================================================================
  PURE SUBROUTINE cp_spin_mgga_x_mn12_mn15( scratch,functional )
    ! ==--------------------------------------------------------------==
    ! Generic spin-unrestricted interface to Minnesota 12
    !                                     12.05.2017 M.P. Bircher @ EPFL

    TYPE(cp_xc_scratch_t), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(in)                             :: scratch
    TYPE(cp_xc_functional_t), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(inout)                          :: functional

    IF (scratch(a)%n >= 1.0e-10_real_8) &
         CALL cp_mgga_x_mn12_mn15( scratch(a), functional(a) )
    IF (scratch(b)%n >= 1.0e-10_real_8) &
         CALL cp_mgga_x_mn12_mn15( scratch(b), functional(b) )

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_spin_mgga_x_mn12_mn15
  ! ==================================================================
  PURE SUBROUTINE cp_mgga_x_mn12_mn15( scratch, functional )
    ! ==--------------------------------------------------------------==
    ! Cleaned up, prettified, reintroduced thresholds
    ! Adapted from OLDCODE to xc_driver format
    !                                     12.05.2017 M.P. Bircher @ EPFL
    ! Porting based on MFM 1.9 at http://comp.chem.umn.edu/mfm/
    !                                        2015 P. Lopez-Tarifa @ EPFL
    ! ==--------------------------------------------------------------==

    TYPE(cp_xc_scratch_t), INTENT(in)        :: scratch
    TYPE(cp_xc_functional_t), INTENT(inout)  :: functional

    REAL(real_8), PARAMETER                  :: ax        = -1.5_real_8*(four_thirds*pi)**(-one_third)
    REAL(real_8), PARAMETER                  :: small_tau =  0.50e-14_real_8
    REAL(real_8), PARAMETER                  :: small_g_2 =  0.25e-14_real_8
    REAL(real_8), PARAMETER                  :: g         =  0.004_real_8
    REAL(real_8), PARAMETER                  :: ome       =  2.5_real_8
    REAL(real_8), PARAMETER                  :: ctau      =  ((6.0_real_8*pi*pi)**two_thirds) 
    REAL(real_8), PARAMETER                  :: ctaueg    =  three_fifths*ctau 

    REAL(real_8)                             :: tau, s2, e, tauueg, tsig
    REAL(real_8)                             :: fu, fu_2, fu_3, fu_4, fu_5
    REAL(real_8)                             :: ft, ft_2, ft_3, ft_4, ft_5
    REAL(real_8)                             :: fv, fv_2, fv_3, fv_4, fv_5
    REAL(real_8)                             :: fmn12_mn15, er, s, sr, sg, us, deno
    REAL(real_8)                             :: dfvdn, dwdt, dtdn, dtdtau, dftdn, dftdtau
    REAL(real_8)                             :: dfmn12_mn15dfv, dfmn12_mn15dfy, dfmn12_mn15dft, dfmn12_mn15dfu
    REAL(real_8)                             :: dfmn12_mn15dn, dfmn12_mn15dg, dfmn12_mn15dt

    tau = 2.0_real_8*scratch%tau

    if ( abs(tau) > small_tau .and. scratch%g_2 > small_g_2 ) then
       s2     = scratch%g_2/(scratch%n_4_3 * scratch%n_4_3)
       e      = ax*scratch%n_4_3
       tauueg = ctaueg * (scratch%n_4_3 * scratch%n_1_3)
       tsig   = tauueg/tau

       fu     = g*s2/(1.0_real_8 + g*s2)
       ft     = (tsig - 1.0_real_8)/(tsig + 1.0_real_8)
       fv     = ome*scratch%n_1_3/(1.0_real_8 + ome*scratch%n_1_3)

       fu_2   = fu*fu
       fu_3   = fu*fu_2
       fu_4   = fu*fu_3
       fu_5   = fu*fu_4
       ft_2   = ft*ft
       ft_3   = ft*ft_2
       ft_4   = ft*ft_3
       ft_5   = ft*ft_4
       fv_2   = fv*fv
       fv_3   = fv*fv_2
       fv_4   = fv*fv_3
       fv_5   = fv*fv_4

       fmn12_mn15  =                cp_mgga_x_param%MN12_MN15%cc000 &
                +           ft*cp_mgga_x_param%MN12_MN15%cc001 &
                +         ft_2*cp_mgga_x_param%MN12_MN15%cc002 &
                +         ft_3*cp_mgga_x_param%MN12_MN15%cc003 &
                +         ft_4*cp_mgga_x_param%MN12_MN15%cc004 &
                +         ft_5*cp_mgga_x_param%MN12_MN15%cc005 &
                +           fu*cp_mgga_x_param%MN12_MN15%cc010 &
                +        fu*ft*cp_mgga_x_param%MN12_MN15%cc011 &
                +      ft_2*fu*cp_mgga_x_param%MN12_MN15%cc012 &
                +      ft_3*fu*cp_mgga_x_param%MN12_MN15%cc013 &
                +      ft_4*fu*cp_mgga_x_param%MN12_MN15%cc014 &
                +         fu_2*cp_mgga_x_param%MN12_MN15%cc020 &
                +      ft*fu_2*cp_mgga_x_param%MN12_MN15%cc021 &
                +    ft_2*fu_2*cp_mgga_x_param%MN12_MN15%cc022 &
                +    ft_3*fu_2*cp_mgga_x_param%MN12_MN15%cc023 &
                +         fu_3*cp_mgga_x_param%MN12_MN15%cc030 &
                +      ft*fu_3*cp_mgga_x_param%MN12_MN15%cc031 &
                +    ft_2*fu_3*cp_mgga_x_param%MN12_MN15%cc032 &
                +           fv*cp_mgga_x_param%MN12_MN15%cc100 &
                +        fv*ft*cp_mgga_x_param%MN12_MN15%cc101 &
                +      fv*ft_2*cp_mgga_x_param%MN12_MN15%cc102 &
                +      fv*ft_3*cp_mgga_x_param%MN12_MN15%cc103 &
                +      fv*ft_4*cp_mgga_x_param%MN12_MN15%cc104 &
                +        fv*fu*cp_mgga_x_param%MN12_MN15%cc110 &
                +     fv*ft*fu*cp_mgga_x_param%MN12_MN15%cc111 &
                +   fv*ft_2*fu*cp_mgga_x_param%MN12_MN15%cc112 &
                +   fv*ft_3*fu*cp_mgga_x_param%MN12_MN15%cc113 &
                +      fv*fu_2*cp_mgga_x_param%MN12_MN15%cc120 &
                +   fv*ft*fu_2*cp_mgga_x_param%MN12_MN15%cc121 &
                + fv*ft_2*fu_2*cp_mgga_x_param%MN12_MN15%cc122 &
                +         fv_2*cp_mgga_x_param%MN12_MN15%cc200 &
                +      fv_2*ft*cp_mgga_x_param%MN12_MN15%cc201 &
                +    fv_2*ft_2*cp_mgga_x_param%MN12_MN15%cc202 &
                +    fv_2*ft_3*cp_mgga_x_param%MN12_MN15%cc203 &
                +      fv_2*fu*cp_mgga_x_param%MN12_MN15%cc210 &
                +   fv_2*ft*fu*cp_mgga_x_param%MN12_MN15%cc211 &
                + fv_2*ft_2*fu*cp_mgga_x_param%MN12_MN15%cc212 &
                +         fv_3*cp_mgga_x_param%MN12_MN15%cc300 &
                +      fv_3*ft*cp_mgga_x_param%MN12_MN15%cc301 &
                +    fv_3*ft_2*cp_mgga_x_param%MN12_MN15%cc302

       er      = four_thirds*e/scratch%n
       s       = sqrt(s2)
       sr      = -four_thirds*s/scratch%n
       sg      = s/scratch%abs_g
       us      = 2.0_real_8*g*s/((1.0_real_8 + g*s2)*(1.0_real_8 + g*s2))
       deno    = (1.0_real_8 + ome*scratch%n_1_3)
       deno    = deno*deno
       dfvdn   = (ome/(3.0_real_8*scratch%n_1_3*scratch%n_1_3))/deno
       dwdt    = 2.0_real_8/((1.0_real_8 + tsig)*(1.0_real_8 + tsig))
       dtdn    = ctau*(scratch%n_1_3*scratch%n_1_3)/tau
       dtdtau  = -tauueg/(tau*tau)
       dftdn   = dwdt*dtdn
       dftdtau = dwdt*dtdtau

       dfmn12_mn15dfv =                         cp_mgga_x_param%MN12_MN15%cc100 &
                   +                    ft*cp_mgga_x_param%MN12_MN15%cc101 &
                   +                  ft_2*cp_mgga_x_param%MN12_MN15%cc102 &
                   +                  ft_3*cp_mgga_x_param%MN12_MN15%cc103 &
                   +                  ft_4*cp_mgga_x_param%MN12_MN15%cc104 &
                   +                    fu*cp_mgga_x_param%MN12_MN15%cc110 &
                   +                 ft*fu*cp_mgga_x_param%MN12_MN15%cc111 &
                   +               ft_2*fu*cp_mgga_x_param%MN12_MN15%cc112 &
                   +               ft_3*fu*cp_mgga_x_param%MN12_MN15%cc113 &
                   +                  fu_2*cp_mgga_x_param%MN12_MN15%cc120 &
                   +               ft*fu_2*cp_mgga_x_param%MN12_MN15%cc121 &
                   +             ft_2*fu_2*cp_mgga_x_param%MN12_MN15%cc122 &
                   +         2.0_real_8*fv*cp_mgga_x_param%MN12_MN15%cc200 &
                   +      2.0_real_8*fv*ft*cp_mgga_x_param%MN12_MN15%cc201 &
                   +    2.0_real_8*fv*ft_2*cp_mgga_x_param%MN12_MN15%cc202 &
                   +    2.0_real_8*fv*ft_3*cp_mgga_x_param%MN12_MN15%cc203 &
                   +      2.0_real_8*fv*fu*cp_mgga_x_param%MN12_MN15%cc210 &
                   +   2.0_real_8*fv*ft*fu*cp_mgga_x_param%MN12_MN15%cc211 &
                   + 2.0_real_8*fv*ft_2*fu*cp_mgga_x_param%MN12_MN15%cc212 &
                   +       3.0_real_8*fv_2*cp_mgga_x_param%MN12_MN15%cc300 &
                   +    3.0_real_8*fv_2*ft*cp_mgga_x_param%MN12_MN15%cc301 &
                   +  3.0_real_8*fv_2*ft_2*cp_mgga_x_param%MN12_MN15%cc302

       dfmn12_mn15dfu =                         cp_mgga_x_param%MN12_MN15%cc010 &
                   +                    ft*cp_mgga_x_param%MN12_MN15%cc011 &
                   +                  ft_2*cp_mgga_x_param%MN12_MN15%cc012 &
                   +                  ft_3*cp_mgga_x_param%MN12_MN15%cc013 &
                   +                  ft_4*cp_mgga_x_param%MN12_MN15%cc014 &
                   +         2.0_real_8*fu*cp_mgga_x_param%MN12_MN15%cc020 &
                   +      2.0_real_8*ft*fu*cp_mgga_x_param%MN12_MN15%cc021 &
                   +    2.0_real_8*ft_2*fu*cp_mgga_x_param%MN12_MN15%cc022 &
                   +    2.0_real_8*ft_3*fu*cp_mgga_x_param%MN12_MN15%cc023 &
                   +       3.0_real_8*fu_2*cp_mgga_x_param%MN12_MN15%cc030 & 
                   +    3.0_real_8*ft*fu_2*cp_mgga_x_param%MN12_MN15%cc031 &
                   +  3.0_real_8*ft_2*fu_2*cp_mgga_x_param%MN12_MN15%cc032 &
                   +                    fv*cp_mgga_x_param%MN12_MN15%cc110 &
                   +                 fv*ft*cp_mgga_x_param%MN12_MN15%cc111 & 
                   +               fv*ft_2*cp_mgga_x_param%MN12_MN15%cc112 &
                   +               fv*ft_3*cp_mgga_x_param%MN12_MN15%cc113 &  
                   +      2.0_real_8*fv*fu*cp_mgga_x_param%MN12_MN15%cc120 &
                   +   2.0_real_8*fv*ft*fu*cp_mgga_x_param%MN12_MN15%cc121 &
                   + 2.0_real_8*fv*ft_2*fu*cp_mgga_x_param%MN12_MN15%cc122 & 
                   +                  fv_2*cp_mgga_x_param%MN12_MN15%cc210 &
                   +               fv_2*ft*cp_mgga_x_param%MN12_MN15%cc211 &
                   +             fv_2*ft_2*cp_mgga_x_param%MN12_MN15%cc212

       dfmn12_mn15dft =                         cp_mgga_x_param%MN12_MN15%cc001 &
                   +         2.0_real_8*ft*cp_mgga_x_param%MN12_MN15%cc002 &
                   +       3.0_real_8*ft_2*cp_mgga_x_param%MN12_MN15%cc003 &
                   +       4.0_real_8*ft_3*cp_mgga_x_param%MN12_MN15%cc004 &
                   +       5.0_real_8*ft_4*cp_mgga_x_param%MN12_MN15%cc005 &
                   +                    fu*cp_mgga_x_param%MN12_MN15%cc011 &
                   +      2.0_real_8*ft*fu*cp_mgga_x_param%MN12_MN15%cc012 &
                   +    3.0_real_8*ft_2*fu*cp_mgga_x_param%MN12_MN15%cc013 &
                   +    4.0_real_8*ft_3*fu*cp_mgga_x_param%MN12_MN15%cc014 &
                   +                  fu_2*cp_mgga_x_param%MN12_MN15%cc021 &
                   +    2.0_real_8*ft*fu_2*cp_mgga_x_param%MN12_MN15%cc022 &
                   +  3.0_real_8*ft_2*fu_2*cp_mgga_x_param%MN12_MN15%cc023 &
                   +                  fu_3*cp_mgga_x_param%MN12_MN15%cc031 &
                   +    2.0_real_8*ft*fu_3*cp_mgga_x_param%MN12_MN15%cc032 &
                   +                    fv*cp_mgga_x_param%MN12_MN15%cc101 &
                   +      2.0_real_8*fv*ft*cp_mgga_x_param%MN12_MN15%cc102 &
                   +    3.0_real_8*fv*ft_2*cp_mgga_x_param%MN12_MN15%cc103 &
                   +    4.0_real_8*fv*ft_3*cp_mgga_x_param%MN12_MN15%cc104 &
                   +                 fv*fu*cp_mgga_x_param%MN12_MN15%cc111 &
                   +   2.0_real_8*fv*ft*fu*cp_mgga_x_param%MN12_MN15%cc112 &
                   + 3.0_real_8*fv*ft_2*fu*cp_mgga_x_param%MN12_MN15%cc113 &
                   +               fv*fu_2*cp_mgga_x_param%MN12_MN15%cc121 &
                   + 2.0_real_8*fv*ft*fu_2*cp_mgga_x_param%MN12_MN15%cc122 &
                   +                  fv_2*cp_mgga_x_param%MN12_MN15%cc201 &
                   +    2.0_real_8*fv_2*ft*cp_mgga_x_param%MN12_MN15%cc202 &
                   +  3.0_real_8*fv_2*ft_2*cp_mgga_x_param%MN12_MN15%cc203 &
                   +               fv_2*fu*cp_mgga_x_param%MN12_MN15%cc211 &
                   + 2.0_real_8*fv_2*ft*fu*cp_mgga_x_param%MN12_MN15%cc212 &
                   +                  fv_3*cp_mgga_x_param%MN12_MN15%cc301 &
                   +    2.0_real_8*fv_3*ft*cp_mgga_x_param%MN12_MN15%cc302

       dfmn12_mn15dn = dfmn12_mn15dfv*dfvdn + dfmn12_mn15dfu*us*sr + dfmn12_mn15dft*dftdn
       dfmn12_mn15dg = dfmn12_mn15dfu*us*sg
       dfmn12_mn15dt = dfmn12_mn15dft*dftdtau

       functional%sx_sigma =  e*fmn12_mn15
       functional%dsx_dn   = er*fmn12_mn15 + e*dfmn12_mn15dn
       functional%dsx_dg   = e*dfmn12_mn15dg
       functional%dsx_dt   = 2.0_real_8*e*dfmn12_mn15dt
    else
       functional%sx_sigma = 0.0_real_8 
       functional%dsx_dn   = 0.0_real_8 
       functional%dsx_dg   = 0.0_real_8 
       functional%dsx_dt   = 0.0_real_8 
    end if
    !
    functional%K        = 0.0_real_8 ! if ever needed, one could convert from F[x] -> K[x] here
    functional%dK_dn    = 0.0_real_8 ! if ever needed, one could convert from F[x] -> K[x] here
    functional%dK_dg    = 0.0_real_8 ! if ever needed, one could convert from F[x] -> K[x] here

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_mgga_x_mn12_mn15
  ! ==================================================================
  PURE SUBROUTINE cp_spin_mgga_x_vs98( scratch,functional )
    ! ==--------------------------------------------------------------==

    TYPE(cp_xc_scratch_t), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(in)                             :: scratch
    TYPE(cp_xc_functional_t), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(inout)                          :: functional

    IF (scratch(a)%n >= 1.0e-10_real_8) &
         CALL cp_mgga_x_vs98( scratch(a), functional(a) )
    IF (scratch(b)%n >= 1.0e-10_real_8) &
         CALL cp_mgga_x_vs98( scratch(b), functional(b) )

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_spin_mgga_x_vs98
  ! ==================================================================
  PURE SUBROUTINE cp_mgga_x_vs98( scratch, functional )
    ! ==--------------------------------------------------------------==
    ! Cleaned up, prettified
    ! Adapted from OLDCODE to xc_driver format
    !                                     24.04.2017 M.P. Bircher @ EPFL
    ! Porting based on MFM 1.9 at http://comp.chem.umn.edu/mfm/
    !                                        2015 P. Lopez-Tarifa @ EPFL
    ! ==--------------------------------------------------------------==

    TYPE(cp_xc_scratch_t), INTENT(in)        :: scratch
    TYPE(cp_xc_functional_t), INTENT(inout)  :: functional

    REAL(real_8), PARAMETER                  :: cf     =  9.115599720_real_8 
    REAL(real_8), PARAMETER                  :: gg     =  0.00186726_real_8
    REAL(real_8), PARAMETER                  :: small  =  1.0e-08_real_8

    REAL(real_8)                             :: tau, n_5_3, n_8_3
    REAL(real_8)                             :: x, dxdn, dxdg, z, dzdn, dzdt
    REAL(real_8)                             :: kx
    REAL(real_8)                             :: gx, dgdx, dgdz

    tau = 2.0_real_8 * scratch%tau 

    if ( scratch%n > small ) then 
       n_5_3 = scratch%n_1_3 * scratch%n_4_3
       n_8_3 = scratch%n_4_3 * scratch%n_4_3

       x     = scratch%g_2/n_8_3

       dxdn  = -eight_thirds*x / scratch%n
       dxdg  = 1.0_real_8/n_8_3
       z     = tau/n_5_3 - cf
       dzdn  = -five_thirds*tau/n_8_3
       dzdt  = 1.0_real_8/n_5_3
       kx    = 1.0_real_8 + gg*x + gg*z

       CALL cp_mgga_x_vs98_gvt4(x,z,kx,gg,&
                                cp_mgga_x_param%VS98%r1,&
                                cp_mgga_x_param%VS98%r2,&
                                cp_mgga_x_param%VS98%r3,&
                                cp_mgga_x_param%VS98%r4,&
                                cp_mgga_x_param%VS98%r5,&
                                cp_mgga_x_param%VS98%r6, &
                                gx,dgdx,dgdz)

       functional%sx_sigma = scratch%n_4_3*gx
       functional%dsx_dn   = four_thirds*scratch%n_1_3*gx+scratch%n_4_3*(dgdx*dxdn+dgdz*dzdn)
       functional%dsx_dg   = 2.0_real_8*scratch%abs_g*scratch%n_4_3*dgdx*dxdg
       functional%dsx_dt   = 2.0_real_8*scratch%n_4_3*(dgdz*dzdt)
       !
       functional%K        = 0.0_real_8 ! if ever needed, one could convert from F[x] -> K[x] here
       functional%dK_dn    = 0.0_real_8 ! if ever needed, one could convert from F[x] -> K[x] here
       functional%dK_dg    = 0.0_real_8 ! if ever needed, one could convert from F[x] -> K[x] here
    end if

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_mgga_x_vs98
  ! ==================================================================

  ! ==================================================================
  ! Basic utilities
  ! ==================================================================

  ! ==================================================================
  ELEMENTAL SUBROUTINE cp_mgga_x_vs98_gvt4(x,z,g,ct,a,b,c,d,e,f,gvt,dgvt_dx,dgvt_dz)
    ! ==--------------------------------------------------------------==
    ! Coded from scratch based on J. Chem. Phys 109, 400 (1998)
    ! GVT4 enhancment factor for Van Voorhis-Scuseria exchange functional
    ! WARNING: This routine is highly prone to be badly affected by noise
    !          terms. The present form should be better at balancing this
    !          out than the official routine from Truhlar's group.
    !                                     24.04.2017 M.P. Bircher @ EPFL
    ! ==--------------------------------------------------------------==

      real(real_8)              :: g_2, g_3, dg_dx, dg_dz
    
      real(real_8)              :: defxz, bcxz

      real(real_8), intent(in)  :: x, z, g, ct, a, b, c, d, e, f
      real(real_8), intent(out) :: gvt, dgvt_dx, dgvt_dz

      g_2   = g*g
      g_3   = g_2*g

      bcxz  = (b*x + c*z)
      defxz = (d*x*x + e*x*z + f*z*z)

      gvt   = a/g + bcxz/g_2 + defxz/g_3

      dg_dx = ct 
      dg_dz = ct 

      dgvt_dx = -a*dg_dx &
                + (b) + bcxz*(-2.0_real_8/g)*dg_dx &
                + (2.0_real_8*d*x + e*z)/g + defxz*(-3.0_real_8/g_2)*dg_dx
      dgvt_dx = dgvt_dx / g_2
                                                                                              
      dgvt_dz = -a*dg_dz &                                                                
                + (c) + bcxz*(-2.0_real_8/g)*dg_dz &      
                + (e*x + 2.0_real_8*f*z)/g + defxz*(-3.0_real_8/g_2)*dg_dz 
      dgvt_dz = dgvt_dz / g_2

    ! ==--------------------------------------------------------------== 
  END SUBROUTINE cp_mgga_x_vs98_gvt4 
  ! ==================================================================

  ! ==================================================================
  ! Check parameter compatibility
  ! ==================================================================

  ! ==================================================================
  ELEMENTAL FUNCTION cp_mgga_x_check(tag) &
  RESULT (OK)
    ! ==--------------------------------------------------------------==
    ! Checks whether parameters are properly initialised and reasonable
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
    CASE( "CP_MGGA_X_M05_M06",  "CP_SPIN_MGGA_X_M05_M06", &
          "CP_MGGA_X_M08_M11",  "CP_SPIN_MGGA_X_M08_M11", &
          "CP_MGGA_X_MN12_MN15","CP_SPIN_MGGA_X_MN12_MN15", & 
          "CP_MGGA_X_VS98",     "CP_SPIN_MGGA_X_VS98"  )
      OK = cp_mgga_x_param%init
    !
    ! Functionals with hard-coded parameters
    !
    CASE( "CP_MGGA_X_TPSS","CP_SPIN_GGA_X_TPSS" )
      OK = .true.
    !
    ! Any other tag
    !
    CASE DEFAULT
      OK = .false.
    END SELECT

    ! ==--------------------------------------------------------------==
  END FUNCTION cp_mgga_x_check
  ! ==================================================================
END MODULE cp_mgga_exchange_utils
