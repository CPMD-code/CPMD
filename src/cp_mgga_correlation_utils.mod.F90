! ==================================================================
! Provides: - Meta-GGA correlation, rewritten or adapted from OLDCODE
!
! Input:  Spin-polarised density (1/2 n for closed shell systems)
! Output: For spin-component only
! Convention: GGA conventions do not apply, since CAM is not
!             implemented for meta-functionals 
!
!                             18.04.2017 - M. P. Bircher @ LCBC/EPFL
! ==================================================================
MODULE cp_mgga_correlation_utils
  USE kinds,                           ONLY: real_8
  USE cnst,                            ONLY: pi
  USE func,                            ONLY: func1, &
                                             func2
  USE cpfunc_types,                    ONLY: cp_xc_scratch_t, &
                                             cp_xc_functional_t, &
                                             cp_xc_spin_components
  USE cp_lda_correlation_utils,        ONLY: cp_lda_c_pz_P, &
                                             cp_lda_c_pw_P, &
                                             cp_lda_c_ob_pw_P, &
                                             cp_lda_c_vwn_P, &
                                             cp_spin_lda_c_pz, &
                                             cp_spin_lda_c_pw, &
                                             cp_spin_lda_c_ob_pw, &
                                             cp_spin_lda_c_vwn, &
                                             zeta_of
  USE cp_gga_correlation_utils,        ONLY: cp_gga_c_pbe_flex, &
                                             cp_spin_gga_c_pbe_flex, &
                                             cp_gga_c_param, &
                                             cp_gga_c_check
  USE cp_mgga_exchange_utils,          ONLY: cp_mgga_x_vs98_gvt4
 
  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER :: a                           = 1
  INTEGER, PARAMETER :: b                           = 2
  INTEGER, PARAMETER :: ab                          = 3

  REAL(real_8), PARAMETER :: three_quarters_by_pi   = 0.75_real_8/pi
  REAL(real_8), PARAMETER :: one_third              = 1.0_real_8/3.0_real_8
  REAL(real_8), PARAMETER :: two_thirds             = 2.0_real_8/3.0_real_8
  REAL(real_8), PARAMETER :: four_thirds            = 4.0_real_8/3.0_real_8
  REAL(real_8), PARAMETER :: five_thirds            = 5.0_real_8/3.0_real_8
  REAL(real_8), PARAMETER :: seven_thirds           = 7.0_real_8/3.0_real_8
  REAL(real_8), PARAMETER :: eight_thirds           = 8.0_real_8/3.0_real_8
  REAL(real_8), PARAMETER :: eleven_thirds          =11.0_real_8/3.0_real_8
  REAL(real_8), PARAMETER :: three_halfs            = 3.0_real_8/2.0_real_8
  REAL(real_8), PARAMETER :: four_ninths            = 4.0_real_8/9.0_real_8 
  REAL(real_8), PARAMETER :: three_fifths           = 3.0_real_8/5.0_real_8
  REAL(real_8), PARAMETER :: nine_quarters_times_pi = 2.25_real_8*pi
  REAL(real_8), PARAMETER :: nqtp_to_one_third      = nine_quarters_times_pi**(one_third)
  REAL(real_8), PARAMETER :: two_to_one_third       = 2.0_real_8**(one_third)
  REAL(real_8), PARAMETER :: two_to_four_thirds     = 2.0_real_8**(four_thirds)
  REAL(real_8), PARAMETER :: tqbp_to_one_third      = three_quarters_by_pi**(one_third)

  !
  ! Customisable functionals (FLEX type)
  ! 
  TYPE, PRIVATE :: vs98_c_t
     REAL(real_8)                :: r7  = 0.0_real_8 
     REAL(real_8)                :: r8  = 0.0_real_8 
     REAL(real_8)                :: r9  = 0.0_real_8 
     REAL(real_8)                :: r10 = 0.0_real_8
     REAL(real_8)                :: r11 = 0.0_real_8
     REAL(real_8)                :: r12 = 0.0_real_8
     REAL(real_8)                :: r13 = 0.0_real_8
     REAL(real_8)                :: r14 = 0.0_real_8
     REAL(real_8)                :: r15 = 0.0_real_8
     REAL(real_8)                :: r16 = 0.0_real_8
     REAL(real_8)                :: r17 = 0.0_real_8
     REAL(real_8)                :: r18 = 0.0_real_8
  END TYPE vs98_c_t
  !
  TYPE, PRIVATE :: m05_m06_c_t
     REAL(real_8)                :: sopp0 = 0.0_real_8 
     REAL(real_8)                :: sopp1 = 0.0_real_8 
     REAL(real_8)                :: sopp2 = 0.0_real_8 
     REAL(real_8)                :: sopp3 = 0.0_real_8 
     REAL(real_8)                :: sopp4 = 0.0_real_8 
     REAL(real_8)                :: sss0  = 0.0_real_8
     REAL(real_8)                :: sss1  = 0.0_real_8
     REAL(real_8)                :: sss2  = 0.0_real_8
     REAL(real_8)                :: sss3  = 0.0_real_8
     REAL(real_8)                :: sss4  = 0.0_real_8
     LOGICAL                     :: add_vs98 = .false.
  END TYPE m05_m06_c_t
  !
  TYPE, PRIVATE :: m08_m11_c_t
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
  END TYPE m08_m11_c_t
  !
  TYPE, PRIVATE :: cp_mgga_c_param_t
     LOGICAL                     :: init = .false.
     TYPE(vs98_c_t)              :: VS98
     TYPE(m05_m06_c_t)           :: M05_M06
     TYPE(m08_m11_c_t)           :: M08_M11
  END TYPE cp_mgga_c_param_t
  !
  TYPE(cp_mgga_c_param_t), PUBLIC, SAVE :: cp_mgga_c_param

  PUBLIC :: CP_MGGA_C_TPSS
  PUBLIC :: CP_MGGA_C_M05_M06
  PUBLIC :: CP_MGGA_C_M08_M11
  PUBLIC :: CP_SPIN_MGGA_C_TPSS
  PUBLIC :: CP_SPIN_MGGA_C_M05_M06
  PUBLIC :: CP_SPIN_MGGA_C_M08_M11

  PUBLIC :: CP_MGGA_C_CHECK

  PRIVATE :: CP_MGGA_C_VS98
  PRIVATE :: CP_SPIN_MGGA_C_VS98

  PRIVATE :: CP_MGGA_C_REV_PKZB_DERIVATIVES
  PRIVATE :: CP_SPIN_MGGA_C_REV_PKZB_DERIVATIVES
  PRIVATE :: CP_SPIN_MGGA_C_M05_M06_PARALLEL
  PRIVATE :: CP_SPIN_MGGA_C_VS98_PARALLEL

CONTAINS

  ! ==================================================================
  PURE SUBROUTINE cp_spin_mgga_c_tpss( scratch, functional )
    ! ==--------------------------------------------------------------==
    !
    ! Adapted from OLDCODE

      real(real_8), parameter                 :: d     = 2.8_real_8
      real(real_8), parameter                 :: small = 1.e-14_real_8 

      real(real_8)                            :: dedga, dedgab, dedgb
      real(real_8)                            :: dedra, dedrb, dedz, dzdg, dzdr
      real(real_8)                            :: e, ec, edaz, edbz, eta
      real(real_8)                            :: v2sum, vtsum, nedz, eop
      real(real_8)                            :: sc, v1ca, v1cb, v2ca, v2cb, v2cab, vtta, vttb
      real(real_8)                            :: n_ab, grho, tau, rs
      real(real_8)                            :: op, rdzdt, ro, to
      real(real_8)                            :: vca, vcb, z, z_2, z_3

      type(cp_xc_scratch_t), dimension(cp_xc_spin_components), &
                             intent(in)       :: scratch
      type(cp_xc_functional_t), dimension(cp_xc_spin_components), &
                                intent(inout) :: functional


      grho = scratch(a)%g_2 + 2.0_real_8*scratch(ab)%g_2 + scratch(b)%g_2
      tau  = scratch(a)%tau + scratch(b)%tau

      if ( abs(tau) > small .and. grho > small) then
         z     = 0.125_real_8*grho/scratch(ab)%n/tau
         z_2   = z*z
         z_3   = z*z_2
         dzdr  = -z/scratch(ab)%n
         rdzdt = -z*scratch(ab)%n/tau
         dzdg  = z/grho
         !
         CALL cp_spin_mgga_c_rev_pkzb_derivatives(scratch,z,e,dedra,dedrb, &
                                                  dedga,dedgb,dedgab,dedz)
         !
         op    = 1.0_real_8 + d * e * z_3
         eop   = e*op
         sc    = scratch(ab)%n * eop
         edaz  = dedra + dedz*dzdr
         nedz  = scratch(ab)%n*e*d*z_2
         v1ca  = eop + scratch(ab)%n*edaz*op + nedz*(edaz*z + 3.0_real_8*e*dzdr)
         edbz  = dedrb + dedz*dzdr
         v1cb  = eop + scratch(ab)%n*edbz*op + nedz*(edbz*z + 3.0_real_8*e*dzdr)
         ro    = e*d*e*3.0_real_8*z_2
         to    = scratch(ab)%n * ( 1.0_real_8 + 2.0_real_8*d*e*z_3 )
         !
         v2sum = 2.0_real_8*scratch(ab)%n*ro*dzdg + 2.0_real_8*to*dedz*dzdg
         v2ca  = v2sum + to*dedga
         v2cb  = v2sum + to*dedgb
         v2cab = v2sum + to*dedgab
         vtsum = dedz*rdzdt*op + e*d*dedz*rdzdt*z_3
         vtta  = vtsum + ro*rdzdt
         vttb  = vtsum + ro*rdzdt
         !
         functional(a)%sc   = 0.0_real_8
         functional(b)%sc   = 0.0_real_8
         functional(ab)%sc  = sc
         !
         functional(a)%v1c  = v1ca
         functional(b)%v1c  = v1cb
         functional(ab)%v1c = 0.0_real_8
         !
         functional(a)%v2c  = v2ca
         functional(b)%v2c  = v2cb
         functional(ab)%v2c = v2cab
         !
         functional(a)%vtc  = vtta
         functional(b)%vtc  = vttb
         functional(ab)%vtc = 0.0_real_8
      else
         select case(trim(adjustl(cp_gga_c_param%pbe%lda_c)))
         case("CP_LDA_C_VWN")
            CALL cp_spin_lda_c_vwn(scratch,functional)
         case("CP_LDA_C_PZ")
            CALL cp_spin_lda_c_pz(scratch,functional)
         case("CP_LDA_C_OB_PW")
            CALL cp_spin_lda_c_ob_pw(scratch,functional)
         case default
            CALL cp_spin_lda_c_pw(scratch,functional)
         end select
      endif

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_spin_mgga_c_tpss
  ! ==================================================================
  PURE SUBROUTINE cp_mgga_c_tpss( scratch, functional )
    ! ==--------------------------------------------------------------==
    !
    ! Adapted from OLDCODE

      real(real_8), parameter                 :: d     = 2.8_real_8 
      real(real_8), parameter                 :: small = 1.e-14_real_8 

      real(real_8)                            :: dedg, dedr, dedt, dedz, dzdg, &
                                                 dzdr, dzdt, e, edgt, edrt, &
                                                 op, ro, rs_ab, z, z_2
      real(real_8)                            :: n_ab, grho, tau, abs_g

      type(cp_xc_scratch_t), intent(in)       :: scratch
      type(cp_xc_functional_t), intent(inout) :: functional

      n_ab  = 2.0_real_8 * scratch%n
      grho  = 4.0_real_8 * scratch%g_2
      abs_g = 2.0_real_8 * scratch%abs_g
      tau   = 2.0_real_8 * scratch%tau

      if ( abs(tau) > small .and. grho > small ) then
         z    = 0.125_real_8*grho/n_ab/tau
         z_2  = z*z
         dzdr = -z/n_ab
         dzdg = 2.0_real_8*z/abs_g
         dzdt = -z/tau

         CALL cp_mgga_c_rev_pkzb_derivatives(scratch,z,e,dedr,dedg,dedz)

         dedt = dedz*dzdt
         edrt = dedr + dedz*dzdr
         edgt = dedg + dedz*dzdg
         op   = 1.0_real_8 + d * e * z_2*z
         ro   = n_ab * e * d * z_2

         functional%sc  = n_ab * e * op
         functional%v1c = (e + n_ab*edrt)*op + ro * (edrt*z + 3.0_real_8*e*dzdr)
         functional%v2c = ( n_ab*edgt*op + ro * (edgt*z + 3.0_real_8*e*dzdg) ) / abs_g
         functional%vtc = n_ab*dedt*op + ro * (dedt*z + 3.0_real_8*e*dzdt)

      else
         rs_ab = tqbp_to_one_third/(two_to_one_third*scratch%n_1_3)
         select case(trim(adjustl(cp_gga_c_param%pbe%lda_c)))
         case("CP_SPIN_LDA_C_VWN")
            CALL cp_lda_c_vwn_P( n_ab, rs_ab, functional )
         case("CP_SPIN_LDA_C_PZ")
            CALL cp_lda_c_pz_P( n_ab, rs_ab, functional )
         case("CP_SPIN_LDA_C_OB_PW")
            CALL cp_lda_c_ob_pw_P( n_ab, rs_ab, functional )
         case default
            CALL cp_lda_c_pw_P( n_ab, rs_ab, functional )
         end select
      endif

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_mgga_c_tpss
  ! ==================================================================
  PURE SUBROUTINE cp_spin_mgga_c_m05_m06( scratch, functional )
    ! ==--------------------------------------------------------------==
    ! Generalisation to the spin-dependent case
    !                                     24.04.2017 M.P. Bircher @ EPFL

    real(real_8), parameter                 :: copp  = 0.0031_real_8 
    real(real_8), parameter                 :: small = 1.0e-14_real_8 

    real(real_8)                            :: sc, scta, sctb, sctab
    real(real_8)                            :: v1ct, v2ct, v3ct
    real(real_8)                            :: v1cta, v2cta, v3cta, v1ctb, v2ctb, v3ctb
    real(real_8)                            :: taua,eua,euega,chia,eupa,chiap,chiag
    real(real_8)                            :: taub,eub,euegb,chib,eupb,chibp,chibg
    real(real_8)                            :: eueg,u,w,dudchia,dudchib,dudpa,dudpb,dudga,dudgb
    real(real_8)                            :: dwdu, dwdpa, dwdpb, dwdga, dwdgb, euegpa, euegpb

    type(cp_xc_scratch_t), dimension(cp_xc_spin_components), &
                           intent(in)       :: scratch
    type(cp_xc_functional_t), dimension(cp_xc_spin_components), &
                              intent(inout) :: functional


    if (cp_mgga_c_param%M05_M06%add_vs98) CALL cp_spin_mgga_c_vs98( scratch, functional )

    taua = 2.0_real_8 * scratch(a)%tau
    taub = 2.0_real_8 * scratch(b)%tau

    !
    ! Note: functional(:) is always properly zeroed before this routine is called; if no contribution
    ! from VS98 is added, this is simply naught.
    !
    sctab = functional(ab)%sc

    if (scratch(a)%n > small) then
       CALL cp_spin_mgga_c_m05_m06_parallel(scratch(a),scta,v1cta,v2cta,v3cta,eua,chia,eupa,chiap,chiag)
       scta  = functional(a)%sc  + scta
       v1cta = functional(a)%v1c + v1cta
       v2cta = 0.5_real_8*functional(a)%v2c + v2cta
       v3cta = functional(a)%vtc + 2.0_real_8*v3cta
    else
       scta  = functional(a)%sc
       v1cta = functional(a)%v1c
       v2cta = 0.5_real_8*functional(a)%v2c ! Cf. convention of output in cp_spin_mgga_c_VS98
       v3cta = functional(a)%vtc
       eua   = 0.0_real_8
       chia  = 0.0_real_8
       eupa  = 0.0_real_8
       chiap = 0.0_real_8
       chiag = 0.0_real_8
    end if

    if (scratch(b)%n > small) then
       CALL cp_spin_mgga_c_m05_m06_parallel(scratch(b),sctb,v1ctb,v2ctb,v3ctb,eub,chib,eupb,chibp,chibg)
       sctb  = functional(b)%sc  + sctb
       v1ctb = functional(b)%v1c + v1ctb
       v2ctb = 0.5_real_8*functional(b)%v2c + v2ctb
       v3ctb = functional(b)%vtc + 2.0_real_8*v3ctb
    else
       sctb  = functional(b)%sc
       v1ctb = functional(b)%v1c
       v2ctb = 0.5_real_8*functional(b)%v2c ! Cf. convention of output in cp_spin_mgga_c_VS98
       v3ctb = functional(b)%vtc
       eub   = 0.0_real_8
       chib  = 0.0_real_8
       eupb  = 0.0_real_8
       chibp = 0.0_real_8
       chibg = 0.0_real_8
    end if

    if (scratch(a)%n > small .and. scratch(b)%n > small) then 
       !
       ! Note: functional is overwritten by a call to cp_spin_lda_c_pw
       !
       CALL cp_spin_lda_c_pw( scratch, functional )

       eueg = functional(ab)%sc - eua - eub
       u    = copp*(chia + chib)/(1.0_real_8 + copp*(chia + chib))
       w    = cp_mgga_c_param%M05_M06%sopp0 + u*(cp_mgga_c_param%M05_M06%sopp1 &
              + u*(cp_mgga_c_param%M05_M06%sopp2 + u*(cp_mgga_c_param%M05_M06%sopp3 + u*cp_mgga_c_param%M05_M06%sopp4)))
       sc   = eueg*w

       dudchia = copp/(1.0_real_8 + copp*(chia + chib))**2.0_real_8
       dudchib = dudchia 
       dudpa   = dudchia*chiap
       dudpb   = dudchib*chibp
       dudga   = dudchia*chiag
       dudgb   = dudchib*chibg
       dwdu    = cp_mgga_c_param%M05_M06%sopp1 + u*(2.0_real_8*cp_mgga_c_param%M05_M06%sopp2 &
                 + u*(3.0_real_8*cp_mgga_c_param%M05_M06%sopp3 + u*4.0_real_8*cp_mgga_c_param%M05_M06%sopp4))
       dwdpa   = dwdu*dudpa
       dwdpb   = dwdu*dudpb
       dwdga   = dwdu*dudga
       dwdgb   = dwdu*dudgb

       euegpa  = functional(a)%v1c - eupa
       euegpb  = functional(b)%v1c - eupb
       v1cta   = v1cta + (euegpa*w + eueg*dwdpa)
       v2cta   = v2cta + (eueg*dwdga)
       v1ctb   = v1ctb + (euegpb*w + eueg*dwdpb)
       v2ctb   = v2ctb + (eueg*dwdgb)
    else
       sc      = 0.0_real_8
    end if

    functional(a)%sc   = scta 
    functional(b)%sc   = sctb
    functional(ab)%sc  = sc + sctab
    !
    functional(a)%v1c  = v1cta
    functional(b)%v1c  = v1ctb
    functional(ab)%v1c = 0.0_real_8
    !
    ! For open-shell systems:
    ! d/d_abs_g_sigma F = d/d_abs_g_sigma_2 d_abs_g_sigma_2/d_abs_g_sigma F = 2 * d/d_abs_g_sigma_2 F
    !
    functional(a)%v2c  = 2.0_real_8 * v2cta
    functional(b)%v2c  = 2.0_real_8 * v2ctb
    functional(ab)%v2c = 0.0_real_8
    !
    functional(a)%vtc  = v3cta
    functional(b)%vtc  = v3ctb
    functional(ab)%vtc = 0.0_real_8

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_spin_mgga_c_m05_m06
  ! ==================================================================
  PURE SUBROUTINE cp_mgga_c_m05_m06( scratch, functional )
    ! ==--------------------------------------------------------------==
    ! Integration into xc driver, cleanup, removed redundancies
    !                                     24.04.2017 M.P. Bircher @ EPFL
    ! Porting based on MFM 1.9 at http://comp.chem.umn.edu/mfm/
    !                                        2015 P. Lopez-Tarifa @ EPFL

    real(real_8), parameter                 :: copp  =  0.0031_real_8 
    real(real_8), parameter                 :: small =  1.0e-8_real_8 

    real(real_8)                            :: sct, v1ct, v2ct, v3ct, cocia
    real(real_8)                            :: taua,eua,euega,chia,eupa,chiap,chiag
    real(real_8)                            :: n_ab, rs,eueg,u,w,dudchia,dudchib,dudpa,dudpb,dudga,dudgb
    real(real_8)                            :: dwdu, dwdpa, dwdpb, dwdga, dwdgb, euegpa, euegpb
    type(cp_xc_scratch_t), intent(in)       :: scratch
    type(cp_xc_functional_t), intent(inout) :: functional

    if (cp_mgga_c_param%M05_M06%add_vs98) CALL cp_mgga_c_vs98( scratch, functional )

    taua = 2.0_real_8 * scratch%tau

    if (scratch%n > small) then
       CALL cp_spin_mgga_c_m05_m06_parallel(scratch,sct,v1ct,v2ct,v3ct,eua,chia,eupa,chiap,chiag)
       sct   =  2.0_real_8*sct + functional%sc
       v1ct  =            v1ct + functional%v1c
       v2ct  =            v2ct + functional%v2c
       v3ct  = 2.0_real_8*v3ct + functional%vtc
    else
       sct   = functional%sc 
       v1ct  = functional%v1c
       v2ct  = functional%v2c
       v3ct  = functional%vtc
       eua   = 0.0_real_8
       chia  = 0.0_real_8
       eupa  = 0.0_real_8
       chiap = 0.0_real_8
       chiag = 0.0_real_8
    end if

    if (scratch%n > small) then
       n_ab = 2.0_real_8*scratch%n 
       rs   = tqbp_to_one_third/(two_to_one_third * scratch%n_1_3)

       CALL cp_lda_c_pw_P( n_ab, rs, functional )

       eueg  = functional%sc - 2.0_real_8*eua
       cocia = (1.0_real_8 + copp*(2.0_real_8*chia))
       u     = copp*(2.0_real_8*chia)/cocia
       w     = cp_mgga_c_param%M05_M06%sopp0 + u*(cp_mgga_c_param%M05_M06%sopp1 &
               + u*(cp_mgga_c_param%M05_M06%sopp2 + u*(cp_mgga_c_param%M05_M06%sopp3 + u*cp_mgga_c_param%M05_M06%sopp4)))
       sct   = sct + eueg*w

       dudchia = copp/(cocia*cocia)
       dudpa   = dudchia*chiap
       dudga   = dudchia*chiag
       dwdu    = cp_mgga_c_param%M05_M06%sopp1 + u*(2.0_real_8*cp_mgga_c_param%M05_M06%sopp2 &
               + u*(3.0_real_8*cp_mgga_c_param%M05_M06%sopp3 + u*4.0_real_8*cp_mgga_c_param%M05_M06%sopp4))
       dwdpa   = dwdu*dudpa
       dwdga   = dwdu*dudga

       euegpa  = functional%v1c - eupa
       v1ct    = v1ct + (euegpa*w + eueg*dwdpa)
       v2ct    = v2ct + (eueg*dwdga)
       v3ct    = v3ct
    end if 

    functional%sc  = sct
    !
    ! For closed shell systems:
    ! d/d_n_ab F = d/d(2*n_a) F
    !
    functional%v1c = v1ct
    !
    ! For closed-shell systems:
    ! d/d_abs_g F = d/d_abs_g_sigma_2 d_abs_g_sigma_2/d_abs_g F = d/d_abs_g_sigma_2 F
    !
    functional%v2c = v2ct
    functional%vtc = v3ct

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_mgga_c_m05_m06
  ! ==================================================================
  PURE SUBROUTINE cp_spin_mgga_c_m08_m11( scratch, functional )
    ! ==--------------------------------------------------------------==
    ! Generalisation to open-shell systems
    !                                     11.05.2017 M.P. Bircher @ EPFL

    real(real_8), parameter                 :: small_tau  = 1.00e-12_real_8 
    real(real_8), parameter                 :: small_g_2  = 0.25e-12_real_8 
    real(real_8), parameter                 :: ctau       = ((3.0_real_8*pi*pi)**two_thirds) 
    real(real_8), parameter                 :: ctaueg     = 0.5_real_8*three_fifths*ctau 

    real(real_8)                            :: tau
    real(real_8)                            :: tauueg, tsig, wsig
    real(real_8)                            :: w1, w2, w3, w4, w5, w6, w7, w8, w9, w10, w11
    real(real_8)                            :: fsig1, fsig2
    real(real_8)                            :: sc_lda_ab, v1_lda_a, v1_lda_b
    real(real_8)                            :: sc_pbe_ab, v1_pbe_a, v1_pbe_b, v2_pbe_a, v2_pbe_b, v2_pbe_ab
    real(real_8)                            :: dwdt, dtdn, dtdt, df1dn, df2dn, df1dt, df2dt, df1dw, df2dw
    real(real_8)                            :: e1, e2, de1dna, de1dnb, de2dna, de2dnb 

    type(cp_xc_scratch_t), dimension(cp_xc_spin_components), &
                           intent(in)       :: scratch
    type(cp_xc_functional_t), dimension(cp_xc_spin_components), &
                              intent(inout) :: functional


    tau   = scratch(a)%tau + scratch(b)%tau

    if( abs(tau) > small_tau .and. scratch(ab)%g_2 > small_g_2 ) then
       tauueg = ctaueg * scratch(ab)%n_4_3 * scratch(ab)%n_1_3
       tsig   = tauueg/tau
       wsig   = (tsig - 1.0_real_8)/(tsig + 1.0_real_8)
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
       fsig1  =    cp_mgga_c_param%M08_M11%at00 &
                 + cp_mgga_c_param%M08_M11%at01*w1 &
                 + cp_mgga_c_param%M08_M11%at02*w2 &
                 + cp_mgga_c_param%M08_M11%at03*w3 & 
                 + cp_mgga_c_param%M08_M11%at04*w4 &
                 + cp_mgga_c_param%M08_M11%at05*w5 &
                 + cp_mgga_c_param%M08_M11%at06*w6 &
                 + cp_mgga_c_param%M08_M11%at07*w7 & 
                 + cp_mgga_c_param%M08_M11%at08*w8 &
                 + cp_mgga_c_param%M08_M11%at09*w9 &
                 + cp_mgga_c_param%M08_M11%at10*w10 &
                 + cp_mgga_c_param%M08_M11%at11*w11
       fsig2  =    cp_mgga_c_param%M08_M11%bt00 &
                 + cp_mgga_c_param%M08_M11%bt01*w1 &
                 + cp_mgga_c_param%M08_M11%bt02*w2 &
                 + cp_mgga_c_param%M08_M11%bt03*w3 & 
                 + cp_mgga_c_param%M08_M11%bt04*w4 &
                 + cp_mgga_c_param%M08_M11%bt05*w5 &
                 + cp_mgga_c_param%M08_M11%bt06*w6 &
                 + cp_mgga_c_param%M08_M11%bt07*w7 & 
                 + cp_mgga_c_param%M08_M11%bt08*w8 &
                 + cp_mgga_c_param%M08_M11%bt09*w9 &
                 + cp_mgga_c_param%M08_M11%bt10*w10 &
                 + cp_mgga_c_param%M08_M11%bt11*w11

       CALL cp_spin_lda_c_pw( scratch, functional )
       sc_lda_ab = functional(ab)%sc
       v1_lda_a  = functional(a)%v1c
       v1_lda_b  = functional(b)%v1c

       CALL cp_spin_gga_c_pbe_flex( scratch, functional )
       sc_pbe_ab = functional(ab)%sc - sc_lda_ab
       v1_pbe_a  = functional(a)%v1c - v1_lda_a
       v1_pbe_b  = functional(b)%v1c - v1_lda_b
       v2_pbe_a  = functional(a)%v2c
       v2_pbe_b  = functional(b)%v2c
       v2_pbe_ab = functional(ab)%v2c

       e1     = fsig1*sc_lda_ab
       e2     = fsig2*sc_pbe_ab

       df1dw  =               cp_mgga_c_param%M08_M11%at01 &
                +  2.0_real_8*cp_mgga_c_param%M08_M11%at02*w1 &
                +  3.0_real_8*cp_mgga_c_param%M08_M11%at03*w2 &
                +  4.0_real_8*cp_mgga_c_param%M08_M11%at04*w3 & 
                +  5.0_real_8*cp_mgga_c_param%M08_M11%at05*w4 &
                +  6.0_real_8*cp_mgga_c_param%M08_M11%at06*w5 &
                +  7.0_real_8*cp_mgga_c_param%M08_M11%at07*w6 &
                +  8.0_real_8*cp_mgga_c_param%M08_M11%at08*w7 &
                +  9.0_real_8*cp_mgga_c_param%M08_M11%at09*w8 &
                + 10.0_real_8*cp_mgga_c_param%M08_M11%at10*w9 &
                + 11.0_real_8*cp_mgga_c_param%M08_M11%at11*w10
       df2dw  =               cp_mgga_c_param%M08_M11%bt01 &
                +  2.0_real_8*cp_mgga_c_param%M08_M11%bt02*w1 &
                +  3.0_real_8*cp_mgga_c_param%M08_M11%bt03*w2 &
                +  4.0_real_8*cp_mgga_c_param%M08_M11%bt04*w3 & 
                +  5.0_real_8*cp_mgga_c_param%M08_M11%bt05*w4 &
                +  6.0_real_8*cp_mgga_c_param%M08_M11%bt06*w5 &
                +  7.0_real_8*cp_mgga_c_param%M08_M11%bt07*w6 &
                +  8.0_real_8*cp_mgga_c_param%M08_M11%bt08*w7 &
                +  9.0_real_8*cp_mgga_c_param%M08_M11%bt09*w8 &
                + 10.0_real_8*cp_mgga_c_param%M08_M11%bt10*w9 &
                + 11.0_real_8*cp_mgga_c_param%M08_M11%bt11*w10
       dwdt    = 2.0_real_8/((1.0_real_8 + tsig)*(1.0_real_8 + tsig))
       dtdn    = 5.0_real_8*tsig/(3.0_real_8*scratch(ab)%n)
       dtdt    = -tsig/tau
       df1dn   = df1dw*dwdt*dtdn
       df2dn   = df2dw*dwdt*dtdn
       df1dt   = df1dw*dwdt*dtdt
       df2dt   = df2dw*dwdt*dtdt
       de1dna  = v1_lda_a*fsig1 + sc_lda_ab*df1dn
       de1dnb  = v1_lda_b*fsig1 + sc_lda_ab*df1dn
       de2dna  = v1_pbe_a*fsig2 + sc_pbe_ab*df2dn
       de2dnb  = v1_pbe_b*fsig2 + sc_pbe_ab*df2dn

       functional(a)%sc   = 0.0_real_8
       functional(b)%sc   = 0.0_real_8
       functional(ab)%sc  = e1 + e2 
       functional(a)%v1c  = de1dna + de2dna
       functional(b)%v1c  = de1dnb + de2dnb
       functional(a)%v2c  = v2_pbe_a*fsig2
       functional(b)%v2c  = v2_pbe_b*fsig2
       functional(ab)%v2c = v2_pbe_ab*fsig2
       functional(a)%vtc  = sc_lda_ab*df1dt + sc_pbe_ab*df2dt
       functional(b)%vtc  = functional(a)%vtc
    else
       CALL cp_spin_lda_c_pw( scratch, functional )
    end if

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_spin_mgga_c_m08_m11
  ! ==================================================================
  PURE SUBROUTINE cp_mgga_c_m08_m11( scratch, functional )
    ! ==--------------------------------------------------------------==
    ! Integration into xc driver, cleanup, removed redundancies
    !                                     10.05.2017 M.P. Bircher @ EPFL
    ! Porting based on MFM 1.9 at http://comp.chem.umn.edu/mfm/
    !                                        2015 P. Lopez-Tarifa @ EPFL

    real(real_8), parameter                 :: small_tau  = 1.00e-12_real_8 
    real(real_8), parameter                 :: small_g_2  = 0.25e-12_real_8 
    real(real_8), parameter                 :: ctau       = ((3.0_real_8*pi*pi)**two_thirds) 
    real(real_8), parameter                 :: ctaueg     = 0.5_real_8*three_fifths*ctau & ! convert n_a > n_ab
                                                             *two_to_four_thirds*two_to_one_third

    real(real_8)                            :: n_ab, abs_g, tau, rs_ab
    real(real_8)                            :: tauueg, tsig, wsig
    real(real_8)                            :: w1, w2, w3, w4, w5, w6, w7, w8, w9, w10, w11
    real(real_8)                            :: fsig1, fsig2
    real(real_8)                            :: sc_lda, v1_lda, sc_pbe, v1_pbe, v2_pbe
    real(real_8)                            :: dwdt, dtdn, dtdt, df1dn, df2dn, df1dt, df2dt, df1dw, df2dw
    real(real_8)                            :: e1, e2, de1dn, de2dn 

    type(cp_xc_scratch_t), intent(in)       :: scratch
    type(cp_xc_functional_t), intent(inout) :: functional


    n_ab  = 2.0_real_8*scratch%n
    abs_g = 2.0_real_8*scratch%abs_g
    tau   = 2.0_real_8*scratch%tau

    rs_ab = tqbp_to_one_third/(two_to_one_third*scratch%n_1_3)

    if( abs(tau) > small_tau .and. scratch%g_2 > small_g_2 ) then
       tauueg = ctaueg * scratch%n_4_3 * scratch%n_1_3
       tsig   = tauueg/tau
       wsig   = (tsig - 1.0_real_8)/(tsig + 1.0_real_8)
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
       fsig1  =    cp_mgga_c_param%M08_M11%at00 &
                 + cp_mgga_c_param%M08_M11%at01*w1 &
                 + cp_mgga_c_param%M08_M11%at02*w2 &
                 + cp_mgga_c_param%M08_M11%at03*w3 & 
                 + cp_mgga_c_param%M08_M11%at04*w4 &
                 + cp_mgga_c_param%M08_M11%at05*w5 &
                 + cp_mgga_c_param%M08_M11%at06*w6 &
                 + cp_mgga_c_param%M08_M11%at07*w7 & 
                 + cp_mgga_c_param%M08_M11%at08*w8 &
                 + cp_mgga_c_param%M08_M11%at09*w9 &
                 + cp_mgga_c_param%M08_M11%at10*w10 &
                 + cp_mgga_c_param%M08_M11%at11*w11
       fsig2  =    cp_mgga_c_param%M08_M11%bt00 &
                 + cp_mgga_c_param%M08_M11%bt01*w1 &
                 + cp_mgga_c_param%M08_M11%bt02*w2 &
                 + cp_mgga_c_param%M08_M11%bt03*w3 & 
                 + cp_mgga_c_param%M08_M11%bt04*w4 &
                 + cp_mgga_c_param%M08_M11%bt05*w5 &
                 + cp_mgga_c_param%M08_M11%bt06*w6 &
                 + cp_mgga_c_param%M08_M11%bt07*w7 & 
                 + cp_mgga_c_param%M08_M11%bt08*w8 &
                 + cp_mgga_c_param%M08_M11%bt09*w9 &
                 + cp_mgga_c_param%M08_M11%bt10*w10 &
                 + cp_mgga_c_param%M08_M11%bt11*w11

       CALL cp_lda_c_pw_P( n_ab, rs_ab, functional )
       sc_lda = functional%sc
       v1_lda = functional%v1c

       CALL cp_gga_c_pbe_flex( scratch, functional )
       sc_pbe = functional%sc  - sc_lda 
       v1_pbe = functional%v1c - v1_lda
       v2_pbe = functional%v2c
       e1     = fsig1*sc_lda
       e2     = fsig2*sc_pbe

       df1dw  =               cp_mgga_c_param%M08_M11%at01 &
                +  2.0_real_8*cp_mgga_c_param%M08_M11%at02*w1 &
                +  3.0_real_8*cp_mgga_c_param%M08_M11%at03*w2 &
                +  4.0_real_8*cp_mgga_c_param%M08_M11%at04*w3 & 
                +  5.0_real_8*cp_mgga_c_param%M08_M11%at05*w4 &
                +  6.0_real_8*cp_mgga_c_param%M08_M11%at06*w5 &
                +  7.0_real_8*cp_mgga_c_param%M08_M11%at07*w6 &
                +  8.0_real_8*cp_mgga_c_param%M08_M11%at08*w7 &
                +  9.0_real_8*cp_mgga_c_param%M08_M11%at09*w8 &
                + 10.0_real_8*cp_mgga_c_param%M08_M11%at10*w9 &
                + 11.0_real_8*cp_mgga_c_param%M08_M11%at11*w10
       df2dw  =               cp_mgga_c_param%M08_M11%bt01 &
                +  2.0_real_8*cp_mgga_c_param%M08_M11%bt02*w1 &
                +  3.0_real_8*cp_mgga_c_param%M08_M11%bt03*w2 &
                +  4.0_real_8*cp_mgga_c_param%M08_M11%bt04*w3 & 
                +  5.0_real_8*cp_mgga_c_param%M08_M11%bt05*w4 &
                +  6.0_real_8*cp_mgga_c_param%M08_M11%bt06*w5 &
                +  7.0_real_8*cp_mgga_c_param%M08_M11%bt07*w6 &
                +  8.0_real_8*cp_mgga_c_param%M08_M11%bt08*w7 &
                +  9.0_real_8*cp_mgga_c_param%M08_M11%bt09*w8 &
                + 10.0_real_8*cp_mgga_c_param%M08_M11%bt10*w9 &
                + 11.0_real_8*cp_mgga_c_param%M08_M11%bt11*w10
       dwdt    = 2.0_real_8/((1.0_real_8 + tsig)*(1.0_real_8 + tsig))
       dtdn    = 5.0_real_8*tsig/(3.0_real_8*n_ab)
       dtdt    = -tsig/tau
       df1dn   = df1dw*dwdt*dtdn
       df2dn   = df2dw*dwdt*dtdn
       df1dt   = df1dw*dwdt*dtdt
       df2dt   = df2dw*dwdt*dtdt
       de1dn   = v1_lda*fsig1 + sc_lda*df1dn
       de2dn   = v1_pbe*fsig2 + sc_pbe*df2dn

       functional%sc  = e1 + e2 
       functional%v1c = de1dn + de2dn
       functional%v2c = v2_pbe*fsig2
       functional%vtc = sc_lda*df1dt + sc_pbe*df2dt
    else
       CALL cp_lda_c_pw_P( n_ab, rs_ab, functional )
    end if

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_mgga_c_m08_m11
  ! ==================================================================
  PURE SUBROUTINE cp_spin_mgga_c_vs98( scratch, functional )
    ! ==--------------------------------------------------------------==
    ! Generalisation of cp_mgga_c_vs98
    !                                     24.04.2017 M.P. Bircher @ EPFL

    real(real_8), parameter                 :: small =  1.0e-8_real_8 
    real(real_8), parameter                 :: cf    =  9.115599720_real_8
    real(real_8), parameter                 :: gab   =  0.003049660_real_8

    real(real_8)                            :: taua, taub, sc
    real(real_8)                            :: scta, v1cta, v2cta, v3cta
    real(real_8)                            :: sctb, v1ctb, v2ctb, v3ctb
    real(real_8)                            :: eua, za, chia, eupa, chiap, chiag, zap, zat
    real(real_8)                            :: eub, zb, chib, eupb, chibp, chibg, zbp, zbt
    real(real_8)                            :: eueg, zab, xab, kab, xk, zk
    real(real_8)                            :: gc, dgdx, dgdz
    real(real_8)                            :: dgdrhoa, dgdgrada, dgdtaua, euegpa
    real(real_8)                            :: dgdrhob, dgdgradb, dgdtaub, euegpb

    type(cp_xc_scratch_t), dimension(cp_xc_spin_components), &
                           intent(in)       :: scratch
    type(cp_xc_functional_t), dimension(cp_xc_spin_components), &
                              intent(inout) :: functional


    taua = 2.0_real_8 * scratch(a)%tau
    taub = 2.0_real_8 * scratch(b)%tau

    sc   = 0.0_real_8

    if ( scratch(a)%n > small ) then
       CALL cp_spin_mgga_c_vs98_parallel(scratch(a), scta, v1cta, v2cta, v3cta, &
                                         eua, za, chia, eupa, chiap, chiag, zap, zat)
       v3cta = 2.0_real_8*v3cta
    else
       scta  = 0.0_real_8
       v1cta = 0.0_real_8
       v2cta = 0.0_real_8
       v3cta = 0.0_real_8
       eua   = 0.0_real_8
       za    = 0.0_real_8
       chia  = 0.0_real_8
       eupa  = 0.0_real_8
       chiap = 0.0_real_8
       chiag = 0.0_real_8
       zap   = 0.0_real_8
       zat   = 0.0_real_8
    end if

    if ( scratch(b)%n > small ) then
       CALL cp_spin_mgga_c_vs98_parallel(scratch(b), sctb, v1ctb, v2ctb, v3ctb, &
                                         eub, zb, chib, eupb, chibp, chibg, zbp, zbt)
       v3ctb = 2.0_real_8*v3ctb
    else
       sctb  = 0.0_real_8
       v1ctb = 0.0_real_8
       v2ctb = 0.0_real_8
       v3ctb = 0.0_real_8
       eub   = 0.0_real_8
       zb    = 0.0_real_8
       chib  = 0.0_real_8
       eupb  = 0.0_real_8
       chibp = 0.0_real_8
       chibg = 0.0_real_8
       zbp   = 0.0_real_8
       zbt   = 0.0_real_8
    end if

    if (scratch(a)%n > small .and. scratch(b)%n > small) then 
       CALL cp_spin_lda_c_pw( scratch, functional )

       eueg = functional(ab)%sc - eua - eub

       zab   = za   + zb
       xab   = chia + chib
       kab   = 1.0_real_8 + gab*(xab + zab)
       xk    = xab/kab
       zk    = zab/kab

       CALL cp_mgga_x_vs98_gvt4(xab,zab,kab,gab, &
                                cp_mgga_c_param%VS98%r7, &
                                cp_mgga_c_param%VS98%r8, &
                                cp_mgga_c_param%VS98%r9, &
                                cp_mgga_c_param%VS98%r10,&
                                cp_mgga_c_param%VS98%r11,&
                                cp_mgga_c_param%VS98%r12,&
                                gc,dgdx,dgdz)

       sc       = gc*eueg

       dgdrhoa  = dgdx*chiap + dgdz*zap
       dgdrhob  = dgdx*chibp + dgdz*zbp
       dgdgrada = dgdx*chiag
       dgdgradb = dgdx*chibg
       dgdtaua  = dgdz*zat
       dgdtaub  = dgdz*zbt
       euegpa   = functional(a)%v1c - eupa
       euegpb   = functional(b)%v1c - eupb
       v1cta    = v1cta + (euegpa*gc + eueg*dgdrhoa)
       v1ctb    = v1ctb + (euegpb*gc + eueg*dgdrhob)
       v2cta    = v2cta + (eueg*dgdgrada)
       v2ctb    = v2ctb + (eueg*dgdgradb)
       v3cta    = v3cta + 2.0_real_8*eueg*dgdtaua
       v3ctb    = v3ctb + 2.0_real_8*eueg*dgdtaub
    end if

    functional(a)%sc   = scta 
    functional(b)%sc   = sctb
    functional(ab)%sc  = sc
    !
    functional(a)%v1c  = v1cta
    functional(b)%v1c  = v1ctb
    functional(ab)%v1c = 0.0_real_8
    !
    ! For open-shell systems:
    ! d/d_abs_g_sigma F = d/d_abs_g_sigma_2 d_abs_g_sigma_2/d_abs_g_sigma F = 2 * d/d_abs_g_sigma_2 F
    !
    functional(a)%v2c  = 2.0_real_8 * v2cta
    functional(b)%v2c  = 2.0_real_8 * v2ctb
    functional(ab)%v2c = 0.0_real_8
    !
    functional(a)%vtc  = v3cta
    functional(b)%vtc  = v3ctb
    functional(ab)%vtc = 0.0_real_8

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_spin_mgga_c_vs98
  ! ==================================================================
  PURE SUBROUTINE cp_mgga_c_vs98( scratch, functional )
    ! ==--------------------------------------------------------------==
    ! Integration into xc driver, cleanup, removed redundancies
    !                                     24.04.2017 M.P. Bircher @ EPFL
    ! Porting based on MFM 1.9 at http://comp.chem.umn.edu/mfm/
    !                                        2015 P. Lopez-Tarifa @ EPFL

    real(real_8), parameter                 :: small =  1.0e-8_real_8 
    real(real_8), parameter                 :: cf    =  9.115599720_real_8
    real(real_8), parameter                 :: gab   =  0.003049660_real_8

    real(real_8)                            :: tau, n_ab, rs
    real(real_8)                            :: sc, v1c, v2c, v3c, sct, v1ct, v2ct, v3ct
    real(real_8)                            :: eua, za, chia, eupa, chiap, chiag, zap, zat
    real(real_8)                            :: eueg, zab, xab, kab, xk, zk
    real(real_8)                            :: gc, dgdx, dgdz
    real(real_8)                            :: dgdrhoa, dgdgrada, dgdtaua, euegpa

    type(cp_xc_scratch_t), intent(in)       :: scratch
    type(cp_xc_functional_t), intent(inout) :: functional

    tau = 2.0_real_8 * scratch%tau

    if ( scratch%n > small ) then
       CALL cp_spin_mgga_c_vs98_parallel(scratch, sct, v1ct, v2ct, v3ct, &
                                         eua, za, chia, eupa, chiap, chiag, zap, zat)
       sct  = 2.0_real_8*sct
       v3ct = 2.0_real_8*v3ct

       n_ab = 2.0_real_8*scratch%n 
       rs   = tqbp_to_one_third/(two_to_one_third * scratch%n_1_3)

       CALL cp_lda_c_pw_P( n_ab, rs, functional )

       eueg  = functional%sc - 2.0_real_8*eua
       zab   = 2.0_real_8*za
       xab   = 2.0_real_8*chia
       kab   = 1.0_real_8 + gab*(xab + zab)
       xk    = xab/kab
       zk    = zab/kab

       CALL cp_mgga_x_vs98_gvt4(xab,zab,kab,gab, &
                                cp_mgga_c_param%VS98%r7, &
                                cp_mgga_c_param%VS98%r8, &
                                cp_mgga_c_param%VS98%r9, &
                                cp_mgga_c_param%VS98%r10,&
                                cp_mgga_c_param%VS98%r11,&
                                cp_mgga_c_param%VS98%r12,&
                                gc,dgdx,dgdz)

       sc    = sct + gc*eueg

       dgdrhoa  = dgdx*chiap + dgdz*zap
       dgdgrada = dgdx*chiag
       dgdtaua  = dgdz*zat
       euegpa   = functional%v1c - eupa
       v1c      = v1ct + (euegpa*gc + eueg*dgdrhoa)
       v2c      = v2ct + (eueg*dgdgrada)
       v3c      = v3ct + 2.0_real_8*eueg*dgdtaua

       functional%sc  = sc
       !
       ! For closed shell systems:
       ! d/d_n_ab F = d/d(2*n_a) F
       !
       functional%v1c = v1c 
       !
       ! For closed-shell systems:
       ! d/d_abs_g F = d/d_abs_g_sigma_2 d_abs_g_sigma_2/d_abs_g F = d/d_abs_g_sigma_2 F
       !
       functional%v2c = v2c
       functional%vtc = v3c
    end if

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_mgga_c_vs98
  ! ==================================================================

  ! ==================================================================
  ! Basic routines (revPKZB, m0X parallel correlation)
  ! ==================================================================

  ! ==================================================================
  PURE SUBROUTINE cp_spin_mgga_c_rev_pkzb_derivatives(scratch,z,e, &
                                   dedra,dedrb,dedga,dedgb,dedgab,dedz)
    ! ==--------------------------------------------------------------==
    !
    ! Adapted from OLDCODE

      real(real_8), parameter                  :: c00   = 0.54_real_8
      real(real_8), parameter                  :: small = 1.e-14_real_8 

      real(real_8)                             :: cxe, dcdgaa, dcdgab, dcdgbb, dcdra, dcdrb
      real(real_8)                             :: e2, eca, ecb, ef1, ef2
      real(real_8)                             :: ecpbe, ecpbe_spin_a, ecpbe_spin_b
      real(real_8)                             :: t1a, t1b, t2a, t2aba, t2abb, t2b
      real(real_8)                             :: tv1ca, tv1cb, tv2ca, tv2cb
      real(real_8)                             :: txaba, txabb, tyaba, tyabb
      real(real_8)                             :: tec, teca, tecb, v1ca, v1cb, v2ca, v2cab, v2cb
      real(real_8)                             :: abs_g

      real(real_8), intent(in)                 :: z
      real(real_8), intent(out)                :: e, dedra, dedrb, dedga, dedgb, dedgab, dedz

      type(cp_xc_scratch_t), dimension(cp_xc_spin_components) &
                                               :: spin_scratch_a, spin_scratch_b
      type(cp_xc_scratch_t), dimension(cp_xc_spin_components), &
                             intent(in)        :: scratch
      type(cp_xc_functional_t), dimension(cp_xc_spin_components) &
                                               :: pbe, pbe_a, pbe_b


      CALL get_pbe_spin_scratch_spaces(scratch,spin_scratch_a,spin_scratch_b)
      CALL cp_spin_gga_c_pbe_flex(scratch,pbe)
      CALL cp_spin_gga_c_pbe_flex(spin_scratch_a,pbe_a)
      CALL cp_spin_gga_c_pbe_flex(spin_scratch_b,pbe_b)

      !
      ! Conventions differ from cp_mgga_c_rev_pkzb_derivatives due to differences in
      ! OLDCODE (never touch a running system...)
      !
      tec     = pbe(ab)%sc   / scratch(ab)%n
      teca    = pbe_a(ab)%sc / spin_scratch_a(ab)%n
      tecb    = pbe_b(ab)%sc / spin_scratch_b(ab)%n

      abs_g   = sqrt( scratch(a)%g_2 + 2.0_real_8*scratch(ab)%g_2 + scratch(b)%g_2 )
      v1ca    = ( pbe(a)%v1c - tec ) / scratch(ab)%n
      v1cb    = ( pbe(b)%v1c - tec ) / scratch(ab)%n
      v2ca    = ( pbe(a)%v2c )  / scratch(ab)%n
      v2cb    = ( pbe(b)%v2c )  / scratch(ab)%n
      v2cab   = ( pbe(ab)%v2c ) / scratch(ab)%n

      if (teca > tec) then
         eca = teca
         t1a = ( pbe_a(a)%v1c - teca ) / spin_scratch_a(a)%n
         if (spin_scratch_a(a)%g_2 > small) then
            t2a = ( pbe_a(a)%v2c ) / spin_scratch_a(a)%n  
         else
            t2a = 0.0_real_8
         end if
         t2aba = 0.0_real_8
         txaba = 0.0_real_8
         tyaba = 0.0_real_8
      else
         eca   = tec
         t1a   = v1ca
         t2a   = v2ca
         t2aba = v2cab
         txaba = v1cb
         tyaba = v2cb
      end if

      if (tecb > tec) then
         ecb = tecb
         t1b = ( pbe_b(a)%v1c - tecb ) / spin_scratch_b(a)%n            
         if (spin_scratch_b(a)%g_2 > small) then
            t2b = ( pbe_b(a)%v2c ) / spin_scratch_b(a)%n  
         else
            t2b = 0.0_real_8
         end if
         t2abb = 0.0_real_8
         txabb = 0.0_real_8
         tyabb = 0.0_real_8
      else
         ecb   = tec
         t1b   = v1cb
         t2b   = v2cb
         t2abb = v2cab
         txabb = v1ca
         tyabb = v2ca
      end if

      CALL ccfun(scratch,cxe,dcdra,dcdrb,dcdgaa,dcdgbb,dcdgab)

      ef1    =   1.0_real_8 + cxe * z * z
      ef2    = ( 1.0_real_8 + cxe ) * z * z
      e2     = scratch(a)%n/scratch(ab)%n * eca + scratch(b)%n/scratch(ab)%n * ecb
      e      = tec * ef1 - ef2 * e2
      dedra  = v1ca * ef1 - (-e2 + eca + scratch(a)%n*t1a + scratch(b)%n*txabb)/scratch(ab)%n * ef2 &
               + dcdra*z*z * ( tec - e2 )
      dedrb  = v1cb * ef1 - (-e2 + ecb + scratch(b)%n*t1b + scratch(a)%n*txaba)/scratch(ab)%n * ef2 &
               + dcdrb*z*z * ( tec - e2 )
      dedga  = v2ca * ef1 - ef2 * ( scratch(a)%n/scratch(ab)%n*t2a + scratch(b)%n/scratch(ab)%n*tyabb ) &
               + 2.0_real_8*dcdgaa*z*z * ( tec - e2 )
      dedgb  = v2cb * ef1 - ef2 * ( scratch(b)%n/scratch(ab)%n*t2b + scratch(a)%n/scratch(ab)%n*tyaba ) &
               + 2.0_real_8*dcdgbb*z*z * ( tec - e2 )
      dedgab = v2cab * ef1 - ef2 * (scratch(a)%n/scratch(ab)%n*t2aba + scratch(b)%n/scratch(ab)%n*t2abb) &
               + 2.0_real_8*dcdgab*z*z * ( tec - e2 )
      dedz   = 2.0_real_8 * tec * cxe * z - 2.0_real_8 * ( 1.0_real_8 + cxe ) * z * e2

      CONTAINS

      ! ==--------------------------------------------------------------==
      PURE SUBROUTINE get_pbe_spin_scratch_spaces(scratch,spin_scratch_a, &
                                                          spin_scratch_b)
        !
        ! Does as it says...

        implicit none

        type(cp_xc_scratch_t), dimension(cp_xc_spin_components), &
                               intent(out)       :: spin_scratch_a, spin_scratch_b
        type(cp_xc_scratch_t), dimension(cp_xc_spin_components), &
                               intent(in)        :: scratch

        spin_scratch_a(a)%n      = scratch(a)%n      
        spin_scratch_a(a)%n_1_3  = scratch(a)%n_1_3  
        spin_scratch_a(a)%n_4_3  = scratch(a)%n_4_3  
        spin_scratch_a(a)%g_2    = scratch(a)%g_2   
        spin_scratch_a(a)%abs_g  = scratch(a)%abs_g 
        spin_scratch_a(a)%tau    = scratch(a)%tau   

        spin_scratch_a(b)%n      = 0.0_real_8 
        spin_scratch_a(b)%n_1_3  = 0.0_real_8 
        spin_scratch_a(b)%n_4_3  = 0.0_real_8 
        spin_scratch_a(b)%g_2    = 0.0_real_8 
        spin_scratch_a(b)%abs_g  = 0.0_real_8 
        spin_scratch_a(b)%tau    = 0.0_real_8 

        spin_scratch_a(ab)%n     = scratch(a)%n    
        spin_scratch_a(ab)%n_1_3 = scratch(a)%n_1_3
        spin_scratch_a(ab)%n_4_3 = scratch(a)%n_4_3
        spin_scratch_a(ab)%g_2   = 0.0_real_8 
        spin_scratch_a(ab)%abs_g = 0.0_real_8 
        spin_scratch_a(ab)%tau   = 0.0_real_8 

        spin_scratch_b(a)%n      = scratch(b)%n      
        spin_scratch_b(a)%n_1_3  = scratch(b)%n_1_3  
        spin_scratch_b(a)%n_4_3  = scratch(b)%n_4_3  
        spin_scratch_b(a)%g_2    = scratch(b)%g_2   
        spin_scratch_b(a)%abs_g  = scratch(b)%abs_g 
        spin_scratch_b(a)%tau    = scratch(b)%tau   

        spin_scratch_b(b)%n      = 0.0_real_8 
        spin_scratch_b(b)%n_1_3  = 0.0_real_8 
        spin_scratch_b(b)%n_4_3  = 0.0_real_8 
        spin_scratch_b(b)%g_2    = 0.0_real_8 
        spin_scratch_b(b)%abs_g  = 0.0_real_8 
        spin_scratch_b(b)%tau    = 0.0_real_8 

        spin_scratch_b(ab)%n     = scratch(b)%n    
        spin_scratch_b(ab)%n_1_3 = scratch(b)%n_1_3
        spin_scratch_b(ab)%n_4_3 = scratch(b)%n_4_3
        spin_scratch_b(ab)%g_2   = 0.0_real_8 
        spin_scratch_b(ab)%abs_g = 0.0_real_8 
        spin_scratch_b(ab)%tau   = 0.0_real_8 


      END SUBROUTINE get_pbe_spin_scratch_spaces
      ! ==--------------------------------------------------------------==
      PURE SUBROUTINE ccfun(scratch,cxe,dcdra,dcdrb,dcdgaa,dcdgbb,dcdgab)
        !
        ! Taken from OLDCODE, slightly adapted (left old comments in...)

        real(real_8), parameter                  :: small_2 = small*small
        real(real_8), parameter                  :: pi_2    = pi*pi

        real(real_8)                             :: cnn, dcde, dcdx
        real(real_8)                             :: dedgaa, dedgab, dedgbb
        real(real_8)                             :: dedra, dedrb, dedx, dxdra, dxdrb
        real(real_8)                             :: dzeta0, eta2, oo, zeta, n_ab_2

        real(real_8), intent(out)                :: cxe, dcdra, dcdrb, dcdgaa, dcdgbb, dcdgab
  
        type(cp_xc_scratch_t), dimension(cp_xc_spin_components), &
                               intent(in)        :: scratch
  
   
        n_ab_2 = scratch(ab)%n*scratch(ab)%n
 
        zeta  = zeta_of( scratch )
        if (abs(zeta) > 1.0_real_8) zeta = sign( 1.0_real_8,zeta )

        ! deb
        ! dzeta0  = (zeta-1.0_real_8)**2 * scratch(a)%g_2 + 2.0_real_8*(zeta*zeta-1.0_real_8)*scratch(ab)%g_2
        ! *           + (zeta + 1.0_real_8)**2 * scratch(b)%g_2

        dzeta0 = (zeta-1.0_real_8)**2 * scratch(a)%g_2 + (zeta + 1.0_real_8)**2 * scratch(b)%g_2
        oo     = 0.25_real_8/(n_ab_2)/(3.0_real_8*pi_2*scratch(ab)%n)**two_thirds
        eta2   = dzeta0*oo

        CALL cgefun(zeta,eta2,cnn,dcdx,dcde)

        dxdra  =  2.0_real_8*scratch(b)%n/(n_ab_2)
        dxdrb  = -2.0_real_8*scratch(a)%n/(n_ab_2)

        if (dzeta0 < small_2) then
           dedx   = 0.0_real_8
           dedgaa = 0.0_real_8
           dedgbb = 0.0_real_8
           dedgab = 0.0_real_8
        else
           ! deb    dedx   = (2.0_real_8*(zeta-1.0_real_8) * scratch(a)%g_2 + 4.0_real_8*zeta*scratch(ab)%g_2
           ! deb *            + 2.0_real_8*(zeta + 1.0_real_8) * scratch(b)%g_2)*oo

           dedx   = ( 2.0_real_8*(zeta-1.0_real_8)*scratch(a)%g_2 &
                    + 2.0_real_8*(zeta + 1.0_real_8)*scratch(b)%g_2 ) * oo
           dedgaa = (zeta-1.0_real_8)**2 * oo
           dedgbb = (zeta + 1.0_real_8)**2 * oo

           ! deb    dedgab = 2.0_real_8*oo*(zeta*zeta-1.0_real_8)

           dedgab = 0.0_real_8
        end if

        dedra  = dedx*dxdra - eight_thirds*eta2/scratch(ab)%n
        dedrb  = dedx*dxdrb - eight_thirds*eta2/scratch(ab)%n
        dcdra  = dcdx*dxdra + dcde*dedra
        dcdrb  = dcdx*dxdrb + dcde*dedrb
        dcdgaa = dcde*dedgaa
        dcdgbb = dcde*dedgbb
        dcdgab = dcde*dedgab
        cxe    = cnn

      END SUBROUTINE ccfun
      ! ==--------------------------------------------------------------==
      ELEMENTAL SUBROUTINE cgefun(zeta,eta2,cxe,dcdx,dcde)
        !
        ! Taken from OLDCODE, slightly adapted

        real(real_8), parameter                  :: small_2  =  small*small
        real(real_8), parameter                  :: et_to_st = -eight_thirds**(-seven_thirds) 

        real(real_8)                             :: cx0, dno, dx0, dx, xx, dxx
        real(real_8)                             :: zeta_2, dno_4
 
        real(real_8), intent(in)                 :: zeta, eta2
        real(real_8), intent(out)                :: cxe, dcdx, dcde
  
        if (abs(zeta) < small .and. abs(eta2) < small_2) then
           cxe  = 0.53_real_8
           dcdx = 0.0_real_8
           dcde = 0.0_real_8
        else if( abs(abs(zeta)-1.0_real_8) < small .and. abs(eta2) < small_2 ) then
           zeta_2 = zeta*zeta
           cx0    = 0.53_real_8 + 0.87_real_8*zeta_2 + 0.50_real_8 * zeta_2*zeta_2 + 2.26_real_8 * zeta_2*zeta_2*zeta_2
           dx0    = 2.0_real_8*zeta * ( 0.87_real_8 + zeta_2 + 3.0_real_8**2.26_real_8 * zeta_2*zeta_2 )
           cxe    = cx0
           dcdx   = dx0
           dcde   = 0.0_real_8
        else
           zeta_2  = zeta*zeta
           cx0 = 0.53_real_8 + 0.87_real_8*zeta_2 + 0.50_real_8 * zeta_2*zeta_2 + 2.26_real_8 * zeta_2*zeta_2*zeta_2
           dx0 = 2.0_real_8*zeta * ( 0.87_real_8 + zeta_2 + 3.0_real_8**2.26_real_8 * zeta_2*zeta_2 )
           if (abs(abs(zeta)-1.0_real_8) < small) then
              xx  = 1.0_real_8 / two_to_four_thirds
              dxx = et_to_st*sign(1.0_real_8,zeta)
           else
              xx  = (1.0_real_8 + zeta)**(-four_thirds) + (1.0_real_8 - zeta)**(-four_thirds)
              dxx = (-four_thirds)*(1.0_real_8 + zeta)**(-seven_thirds) - (-four_thirds)*(1.0_real_8 - zeta)**(-seven_thirds)
           end if
           dno   = (1.0_real_8 + 0.5_real_8*eta2*xx)
           dno_4 = dno * dno
           dno_4 = dno_4 * dno_4
           cxe   = cx0/dno_4
           dcdx  = dx0/dno_4 - 4.0_real_8*cxe*0.5_real_8*eta2*dxx/dno
           dcde  = -2.0_real_8*cxe*xx/dno
        end if

      END SUBROUTINE cgefun
      ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_spin_mgga_c_rev_pkzb_derivatives
  ! ==================================================================
  PURE SUBROUTINE cp_mgga_c_rev_pkzb_derivatives(scratch,z,e,dedr,dedg,dedz)
    ! ==--------------------------------------------------------------==
    !
    ! revPKZB taken from OLDCODE and updated.

      real(real_8), parameter                  :: c00 = 0.53_real_8
      real(real_8), parameter                  :: c0p = 1.53_real_8

      real(real_8)                             :: ecpbe, ecpbe_spin, tecpbe, pbe_spin_aa
      real(real_8)                             :: t1, tv1c, tv2c, pv1c, pv2c, z_2

      real(real_8), intent(in)                 :: z
      real(real_8), intent(out)                :: e, dedr, dedg, dedz

      type(cp_xc_scratch_t), dimension(cp_xc_spin_components) &
                                               :: spin_scratch
      type(cp_xc_scratch_t), intent(in)        :: scratch

      type(cp_xc_functional_t), dimension(cp_xc_spin_components) &
                                               :: pbe_spin
      type(cp_xc_functional_t)                 :: pbe

      spin_scratch(a)%n      = scratch%n    
      spin_scratch(a)%n_1_3  = scratch%n_1_3
      spin_scratch(a)%n_4_3  = scratch%n_4_3
      spin_scratch(a)%g_2    = scratch%g_2  
      spin_scratch(a)%abs_g  = scratch%abs_g
      spin_scratch(a)%tau    = scratch%tau  

      spin_scratch(b)%n      = 0.0_real_8 
      spin_scratch(b)%n_1_3  = 0.0_real_8 
      spin_scratch(b)%n_4_3  = 0.0_real_8 
      spin_scratch(b)%g_2    = 0.0_real_8 
      spin_scratch(b)%abs_g  = 0.0_real_8 
      spin_scratch(b)%tau    = 0.0_real_8 

      spin_scratch(ab)%n     = scratch%n    
      spin_scratch(ab)%n_1_3 = scratch%n_1_3
      spin_scratch(ab)%n_4_3 = scratch%n_4_3
      spin_scratch(ab)%g_2   = 0.0_real_8 
      spin_scratch(ab)%abs_g = 0.0_real_8 
      spin_scratch(ab)%tau   = 0.0_real_8 

      CALL cp_gga_c_pbe_flex(scratch,pbe)
      CALL cp_spin_gga_c_pbe_flex(spin_scratch,pbe_spin)

      ecpbe      = pbe%sc  / (2.0_real_8*scratch%n)
      ecpbe_spin = pbe_spin(ab)%sc / spin_scratch(ab)%n

      pv1c       = ( pbe%v1c - ecpbe ) / (2.0_real_8*scratch%n)
      pv2c       = scratch%abs_g*( pbe%v2c ) / scratch%n

      !
      ! tv1c and tv2c are not CPMD-v1c and -v2c, but de_dn and de_dg, respectively (!)
      !
      if ( ecpbe_spin > ecpbe ) then
         tecpbe = ecpbe_spin
         tv1c   = 0.5_real_8 * ( pbe_spin(a)%v1c - ecpbe_spin ) / spin_scratch(a)%n            
         ! We use abs_g_a directly, since in c_pbe, aa = abs_g_a + abs_g_b + 2*abs_g_ab = abs_g_a
         tv2c   = 0.5_real_8 * spin_scratch(a)%abs_g * ( pbe_spin(a)%v2c ) / spin_scratch(a)%n  
      else
         tecpbe = ecpbe
         tv1c   = pv1c
         tv2c   = pv2c
      end if

      z_2  = z*z
      t1   = 1.0_real_8 + c00*z_2
      e    = ( ecpbe*t1 - c0p*z_2*tecpbe )
      dedr = ( pv1c*t1 - c0p*z_2*tv1c )
      dedg = ( pv2c*t1 - c0p*z_2*tv2c )
      dedz = ( 2.0_real_8*z * (ecpbe*c00 - c0p*tecpbe) )

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_mgga_c_rev_pkzb_derivatives
  ! ==================================================================
  PURE SUBROUTINE cp_spin_mgga_c_m05_m06_parallel(scratch,f,fp,fg,ft,eueg,chi,euegp,chip,chig)
    ! ==--------------------------------------------------------------==
    ! Integration into xc driver, cleanup, removed redundancies
    !                                     24.04.2017 M.P. Bircher @ EPFL
    ! Porting based on MFM 1.9 at http://comp.chem.umn.edu/mfm/
    !                                        2015 P. Lopez-Tarifa @ EPFL

      real(real_8), parameter                :: small = 1.0e-8_real_8
      real(real_8), parameter                :: css   = 0.06_real_8

      real(real_8)                           :: tx, d, fscc, dFsccP, dFsccG
      real(real_8)                           :: w, u, dfscct, dudchi, dwdu, dwdp, dWdG
      real(real_8)                           :: potlc, dlds, dldz

      real(real_8), intent(out)              :: f, fp, fg, ft
      real(real_8), intent(out)              :: eueg, chi, euegp, chip, chig

      type(cp_xc_scratch_t), dimension(cp_xc_spin_components) &
                                             :: copy
      type(cp_xc_functional_t), dimension(cp_xc_spin_components) &
                                             :: lda

      type(cp_xc_scratch_t), intent(in)      :: scratch

      copy(a)%n      = scratch%n
      copy(a)%n_1_3  = scratch%n_1_3
      copy(a)%n_4_3  = scratch%n_4_3

      copy(b)%n      = 0.0_real_8 
      copy(b)%n_1_3  = 0.0_real_8 
      copy(b)%n_4_3  = 0.0_real_8 

      copy(ab)%n     = scratch%n
      copy(ab)%n_1_3 = scratch%n_1_3
      copy(ab)%n_4_3 = scratch%n_4_3

      tx = 2.0_real_8 * scratch%tau

      if (scratch%n < small) then 
         eueg   = 0.0_real_8
         chi    = 0.0_real_8
         euegp  = 0.0_real_8
         chip   = 0.0_real_8
         chig   = 0.0_real_8
         f      = 0.0_real_8
         fp     = 0.0_real_8
         fg     = 0.0_real_8
         ft     = 0.0_real_8
      else 
         CALL cp_spin_lda_c_pw(copy,lda)
         eueg   = lda(ab)%sc
         d      = tx-0.25_real_8*scratch%g_2/scratch%n
         chi    = scratch%g_2/(scratch%n_4_3*scratch%n_4_3)
         u      = css*chi/(1.0_real_8 + css*chi)
         w      = cp_mgga_c_param%M05_M06%sss0 + u*(cp_mgga_c_param%M05_M06%sss1 &
                  + u*(cp_mgga_c_param%M05_M06%sss2 + u*(cp_mgga_c_param%M05_M06%sss3 + u*cp_mgga_c_param%M05_M06%sss4)))
         fscc   = d/tx
         f      = fscc*w*eueg
         chig   = 1.0_real_8/(scratch%n_4_3*scratch%n_4_3)
         chip   = -eight_thirds*chi/scratch%n
         dfsccp = 0.25_real_8*scratch%g_2/(tx*scratch%n**2)
         dfsccg = -0.25_real_8/(tx*scratch%n)
         dfscct = 0.25_real_8*scratch%g_2/(scratch%n*tx**2)
         dudchi = css/((1.0_real_8 + css*chi)**2)
         dwdu   = cp_mgga_c_param%M05_M06%sss1 + u*(2.0_real_8*cp_mgga_c_param%M05_M06%sss2 &
                  + u*(3.0_real_8*cp_mgga_c_param%M05_M06%sss3 + u*4.0_real_8*cp_mgga_c_param%M05_M06%sss4))
         dwdp   = dwdu*dudchi*chip
         dwdg   = dwdu*dudchi*chig
         euegp  = lda(a)%v1c 
         fp     = dfsccp*w*eueg + fscc*dwdp*eueg + fscc*w*euegp
         fg     = dfsccg*w*eueg + fscc*dwdg*eueg
         ft     = dfscct*w*eueg
      end if

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_spin_mgga_c_m05_m06_parallel
  ! ==================================================================
  PURE SUBROUTINE cp_spin_mgga_c_vs98_parallel(scratch,f,fp,fg,ft,eueg, &
                                               z,chi,euegp,chip,chig,zp,zt)
    ! ==--------------------------------------------------------------==
    ! Integration into xc driver, cleanup, removed redundancies
    !                                     24.04.2017 M.P. Bircher @ EPFL
    ! Porting based on MFM 1.9 at http://comp.chem.umn.edu/mfm/
    !                                        2015 P. Lopez-Tarifa @ EPFL

      real(real_8), parameter                :: small =  1.0e-8_real_8
      real(real_8), parameter                :: cf    =  9.115599720_real_8
      real(real_8), parameter                :: gcc   =  0.005150880_real_8

      real(real_8)                           :: tx, n_5_3, n_8_3
      real(real_8)                           :: kc, xk, zk, d
      real(real_8)                           :: dz, dx, dp, dg, dt
      real(real_8)                           :: dgdp, dgdg, dgdt
      real(real_8)                           :: gc, dgdx, dgdz

      real(real_8), intent(out)              :: f, fp, fg, ft, eueg, z, chi
      real(real_8), intent(out)              :: euegp, chip, chig, zp, zt

      type(cp_xc_scratch_t), dimension(cp_xc_spin_components) &
                                             :: copy
      type(cp_xc_functional_t), dimension(cp_xc_spin_components) &
                                             :: lda

      type(cp_xc_scratch_t), intent(in)      :: scratch

      copy(a)%n      = scratch%n
      copy(a)%n_1_3  = scratch%n_1_3
      copy(a)%n_4_3  = scratch%n_4_3

      copy(b)%n      = 0.0_real_8 
      copy(b)%n_1_3  = 0.0_real_8 
      copy(b)%n_4_3  = 0.0_real_8 

      copy(ab)%n     = scratch%n
      copy(ab)%n_1_3 = scratch%n_1_3
      copy(ab)%n_4_3 = scratch%n_4_3

      tx = 2.0_real_8 * scratch%tau

      if (scratch%n < small) then 
         eueg   = 0.0_real_8
         chi    = 0.0_real_8
         z      = 0.0_real_8
         zt     = 0.0_real_8
         euegp  = 0.0_real_8
         chip   = 0.0_real_8
         chig   = 0.0_real_8
         f      = 0.0_real_8
         fp     = 0.0_real_8
         fg     = 0.0_real_8
         ft     = 0.0_real_8
      else 
         n_5_3 = scratch%n_1_3 * scratch%n_4_3
         n_8_3 = scratch%n_4_3 * scratch%n_4_3
         CALL cp_spin_lda_c_pw(copy,lda)
         eueg  = lda(ab)%sc
         chi   = scratch%g_2/(n_8_3)
         z     = (tx/n_5_3)-cf
         kc    = 1.0_real_8 + gcc * (chi + z)
         xk    = chi/kc
         zk    = z/kc
         d     = 1.0_real_8-chi/(4.0_real_8*(z + cf))
         CALL cp_mgga_x_vs98_gvt4(chi,z,kc,gcc, &
                                  cp_mgga_c_param%VS98%r13,&
                                  cp_mgga_c_param%VS98%r14,&
                                  cp_mgga_c_param%VS98%r15,&
                                  cp_mgga_c_param%VS98%r16,&
                                  cp_mgga_c_param%VS98%r17,&
                                  cp_mgga_c_param%VS98%r18,&
                                  gc,dgdx,dgdz)
         f     = d*eueg*gc
         chig  =  1.0_real_8/n_8_3
         chip  = -eight_thirds*chi/scratch%n
         zp    = -five_thirds*tx/n_8_3
         zt    = 1.0_real_8/n_5_3
         dz    = chi/(4.0_real_8*(z + cf)*(z + cf))
         dx    = -0.25_real_8/(z + cf)
         dp    = dz*zp + dx*chip
         dg    = dx*chig
         dt    = dz*zt
         dgdp  = dgdx*chip + dgdz*zp
         dgdg  = dgdx*chig
         dgdt  = dgdz*zt
         euegp = lda(a)%v1c
         fp    = dp*eueg*gc + d*euegp*gc + d*eueg*dgdp
         fg    = dg*eueg*gc + d*eueg*dgdg
         ft    = dt*eueg*gc + d*eueg*dgdt
      end if

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_spin_mgga_c_vs98_parallel
  ! ==================================================================

  ! ==================================================================
  ! Check parameter compatibility
  ! ==================================================================

  ! ==================================================================
  ELEMENTAL FUNCTION cp_mgga_c_check(tag) &
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
    CASE( "CP_MGGA_C_M05_M06", "CP_SPIN_MGGA_C_M05_M06",&
          "CP_MGGA_C_VS98",    "CP_SPIN_MGGA_C_VS98" )
      OK = cp_mgga_c_param%init
    CASE( "CP_MGGA_C_M08_M11", "CP_SPIN_MGGA_C_M08_M11" )
      OK = ( cp_mgga_c_param%init .and. &
             cp_gga_c_check("CP_GGA_C_PBE_FLEX") )
    CASE( "CP_MGGA_C_TPSS","CP_SPIN_GGA_C_TPSS" )
      OK = cp_gga_c_check("CP_GGA_C_PBE_FLEX")
    !
    ! Any other tag
    !
    CASE DEFAULT
      OK = .false.
    END SELECT

    ! ==--------------------------------------------------------------==
  END FUNCTION cp_mgga_c_check
  ! ==================================================================
END MODULE cp_mgga_correlation_utils
