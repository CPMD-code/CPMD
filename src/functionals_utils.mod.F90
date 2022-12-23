#ifdef __SR8000
!option MP(P(0)), LANGLVL(SAVE(0))
#endif

MODULE functionals_utils
  USE cnst,                            ONLY: fpi,&
                                             pi
  USE error_handling,                  ONLY: stopgm
  USE func,                            ONLY: &
       func1, func2, func3, mfxcc_is_hedin, mfxcc_is_lyp, mfxcc_is_obpw, &
       mfxcc_is_obpz, mfxcc_is_pade, mfxcc_is_pw, mfxcc_is_pz, &
       mfxcc_is_skipped, mfxcc_is_vwn, mfxcc_is_wigner, mfxcx_is_slaterx, &
       mgcc_is_dfr_zpbec, mgcc_is_ggac, mgcc_is_hse, mgcc_is_lyp, &
       mgcc_is_optc, mgcc_is_pbec, mgcc_is_pbesc, mgcc_is_perdew86, &
       mgcc_is_skipped, mgcx_is_becke88, mgcx_is_dfr_xpbex, &
       mgcx_is_dfr_xpbex_hybrid, mgcx_is_dfr_zpbex, mgcx_is_ggax, &
       mgcx_is_hcth, mgcx_is_optx, mgcx_is_ox, mgcx_is_ox_hybrid, &
       mgcx_is_pbesx, mgcx_is_pbex, mgcx_is_revpbex, mgcx_is_skipped
  USE kinds,                           ONLY: real_8
  USE special_functions,               ONLY: cp_erf,&
                                             cp_erfc
  USE system,                          ONLY: cntl

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: xc
  PUBLIC :: gcxc
  PUBLIC :: slaterx
  PUBLIC :: pz
  PUBLIC :: vwn
  PUBLIC :: lyp
  PUBLIC :: pw
  PUBLIC :: wigner
  PUBLIC :: hedin
  PUBLIC :: obpz
  PUBLIC :: obpw
  PUBLIC :: pade
  PUBLIC :: becke88
  PUBLIC :: ggax
  PUBLIC :: pbesx
  PUBLIC :: pbex
  PUBLIC :: rks_x_zpbe
  PUBLIC :: uks_x_zpbe
  PUBLIC :: revpbex
  PUBLIC :: perdew86
  PUBLIC :: glyp
  PUBLIC :: ggac
  PUBLIC :: pbesc
  PUBLIC :: pbec
  PUBLIC :: hcth
  PUBLIC :: pwcorr
  PUBLIC :: optx
  PUBLIC :: pbexsr_lsd
  PUBLIC :: pbexsr
  PUBLIC :: wpbe_analytical_erfc_approx_grad
  PUBLIC :: expint

CONTAINS

  ! ==================================================================
  SUBROUTINE xc(rho,ex,ec,vx,vc)
    ! ==--------------------------------------------------------------==
    ! ==  LDA EXCHANGE AND CORRELATION FUNCTIONALS                    ==
    ! ==                                                              ==
    ! ==  EXCHANGE  :  SLATER alpha                                   ==
    ! ==  CORRELATION : CEPERLEY & ALDER (PERDEW-ZUNGER PARAMETERS)   ==
    ! ==                VOSKO, WILK & NUSSAIR                         ==
    ! ==                LEE, YANG & PARR                              ==
    ! ==                PERDEW & WANG                                 ==
    ! ==                WIGNER                                        ==
    ! ==                HEDIN & LUNDQVIST                             ==
    ! ==                ORTIZ & BALLONE (PERDEW-ZUNGER FORMULA)       ==
    ! ==                ORTIZ & BALLONE (PERDEW-WANG FORMULA)         ==
    ! ==                HCTH/120                                      ==
    ! ==--------------------------------------------------------------==

    IMPLICIT NONE

    REAL(real_8)                             :: rho, ex, ec, vx, vc

    CHARACTER(len=*), PARAMETER              :: procedureN = 'xc'
    REAL(real_8), PARAMETER                  :: pi34 = 0.75_real_8/pi, &
         small = 1.e-10_real_8 , &
         third = 1._real_8/3._real_8 

    INTEGER                                  :: iflg
    REAL(real_8)                             :: rs

    ! ==--------------------------------------------------------------==
    ! ..Exchange

    IF (func1%mfxcx == mfxcx_is_slaterx) THEN
       CALL slaterx(rho,ex,vx,func2%salpha)
       ex=ex*func3%pxlda
       vx=vx*func3%pxlda
    ELSE
       ex=0.0_real_8
       vx=0.0_real_8
    ENDIF
    ! shouldn't this be checked before, such that slaterx is not calculated in vein?
    IF (rho.LE.small) THEN
       ec = 0.0_real_8
       vc = 0.0_real_8
       ex = 0.0_real_8
       vx = 0.0_real_8
    ELSE
       ! take care of the fact that the LDA part of LYP has to be calculated
       ! separately if LYP is the GGA correlation functional (mgcc=2)
       ! thybrid is in the if statement for prophylactic reasons
       IF (cntl%thybrid .AND. func1%mgcc==2 .AND. .NOT. func1%mfxcc==3) THEN
          CALL mixed_lda_correlation()
          ! This is the old loop.
       ELSE
          IF (func1%mfxcc == mfxcc_is_skipped) THEN
             ec = 0.0_real_8
             vc =  0.0_real_8
          ELSEIF (func1%mfxcc == mfxcc_is_pz) THEN
             rs=(pi34/rho)**third
             iflg=2
             IF (rs.LT.1.0_real_8) iflg=1
             CALL pz(rs,ec,vc,iflg)
          ELSEIF (func1%mfxcc == mfxcc_is_vwn) THEN
             rs = (pi34/rho)**third
             CALL vwn(rs,ec,vc)
          ELSEIF (func1%mfxcc == mfxcc_is_lyp) THEN
             CALL lyp(rho,ec,vc)
          ELSEIF (func1%mfxcc == mfxcc_is_pw) THEN
             rs=(pi34/rho)**third
             iflg=2
             IF (rs.LT.0.5_real_8) iflg=1
             IF (rs.GT.100._real_8) iflg=3
             CALL pw(rs,ec,vc,iflg)
          ELSEIF (func1%mfxcc == mfxcc_is_wigner) THEN
             CALL wigner(rho,ec,vc)
          ELSEIF (func1%mfxcc == mfxcc_is_hedin) THEN
             CALL hedin(rho,ec,vc)
          ELSEIF (func1%mfxcc == mfxcc_is_obpz) THEN
             rs=(pi34/rho)**third
             iflg=2
             IF (rs.LT.1.0_real_8) iflg=1
             CALL obpz(rs,ec,vc,iflg)
          ELSEIF (func1%mfxcc == mfxcc_is_obpw) THEN
             rs=(pi34/rho)**third
             iflg=2
             IF (rs.LT.0.5_real_8) iflg=1
             IF (rs.GT.100._real_8) iflg=3
             CALL obpw(rs,ec,vc,iflg)
          ELSEIF (func1%mfxcc == mfxcc_is_pade) THEN
             rs=(pi34/rho)**third
             CALL pade(rs,ec,vc)
          ELSE
             ! vw MFXCC should be valid at this point!
             CALL stopgm(procedureN,"Correlation functional "//&
                  "not implemented",& 
                  __LINE__,__FILE__)
          ENDIF
          ec=ec*func3%pclda
          vc=vc*func3%pclda
       END IF
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN

  CONTAINS

    !   --------------------------------------------------------------------
    SUBROUTINE mixed_lda_correlation()
      !   --------------------------------------------------------------------
    CHARACTER(len=*), PARAMETER :: procedureN = 'mixed_lda_correlation'

    REAL(real_8)                             :: ec_aux, ec_lyp, vc_aux, vc_lyp

!   --------------------------------------------------------------------

      ec_aux=0.0_real_8 ; vc_aux=0.0_real_8
      ec_lyp=0.0_real_8 ; vc_lyp=0.0_real_8
      IF (func1%mfxcc == mfxcc_is_skipped) THEN
         ec_aux = 0.0_real_8
         vc_aux =  0.0_real_8
      ELSE IF (func1%mfxcc == mfxcc_is_pz) THEN
         rs=(pi34/rho)**third
         iflg=2
         IF (rs.LT.1.0_real_8) iflg=1
         CALL pz(rs,ec_aux,vc_aux,iflg)
      ELSE IF (func1%mfxcc == mfxcc_is_vwn) THEN
         rs = (pi34/rho)**third
         CALL vwn(rs,ec_aux,vc_aux)
      ELSE IF (func1%mfxcc == mfxcc_is_pw) THEN
         rs=(pi34/rho)**third
         iflg=2
         IF (rs.LT.0.5_real_8) iflg=1
         IF (rs.GT.100._real_8) iflg=3
         CALL pw(rs,ec_aux,vc_aux,iflg)
      ELSE IF (func1%mfxcc == mfxcc_is_wigner) THEN
         CALL wigner(rho,ec_aux,vc_aux)
      ELSE IF (func1%mfxcc == mfxcc_is_hedin) THEN
         CALL hedin(rho,ec_aux,vc_aux)
      ELSE IF (func1%mfxcc == mfxcc_is_obpz) THEN
         rs=(pi34/rho)**third
         iflg=2
         IF (rs.LT.1.0_real_8) iflg=1
         CALL obpz(rs,ec_aux,vc_aux,iflg)
      ELSE IF (func1%mfxcc == mfxcc_is_obpw) THEN
         rs=(pi34/rho)**third
         iflg=2
         IF (rs.LT.0.5_real_8) iflg=1
         IF (rs.GT.100._real_8) iflg=3
         CALL obpw(rs,ec_aux,vc_aux,iflg)
      ELSE IF (func1%mfxcc == mfxcc_is_pade) THEN
         rs=(pi34/rho)**third
         CALL pade(rs,ec_aux,vc_aux)
      ELSE
         CALL stopgm(procedureN,"Correlation functional "//&
              "not implemented",& 
              __LINE__,__FILE__)
      END IF
      ! Calculate the LDA-contribution of LYP and add to the real LDA correlation using the appropriate weight
      CALL lyp(rho,ec_lyp,vc_lyp)
      ec=ec_aux*func3%pclda
      vc=vc_aux*func3%pclda
      ec=ec+ec_lyp*func3%pcgc
      vc=vc+vc_lyp*func3%pcgc
      RETURN
      !   --------------------------------------------------------------------
    END SUBROUTINE mixed_lda_correlation
    !   --------------------------------------------------------------------

  END SUBROUTINE xc
  ! ==================================================================
  SUBROUTINE gcxc(rho,grho,sx,sc,v1x,v2x,v1c,v2c)
    ! ==--------------------------------------------------------------==
    ! ==  GRADIENT CORRECTIONS FOR EXCHANGE AND CORRELATION           ==
    ! ==                                                              ==
    ! ==  EXCHANGE  :  BECKE88                                        ==
    ! ==               GGAX                                           ==
    ! ==               PBEX                                           ==
    ! ==               PBESX                                           ==
    ! ==               revPBEX                                        ==
    ! ==               HCTH/120                                       ==
    ! ==               OPTX                                           ==
    ! ==  CORRELATION : PERDEW86                                      ==
    ! ==                LEE, YANG & PARR                              ==
    ! ==                GGAC                                          ==
    ! ==                PBEC                                          ==
    ! ==                PBESC                                          ==
    ! ==                HCTH/120                                      ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rho, grho, sx, sc, v1x, v2x, &
                                                v1c, v2c

    CHARACTER(len=*), PARAMETER              :: procedureN = 'gcxc'
    REAL(real_8), PARAMETER                  :: small = 1.e-10_real_8 

    REAL(real_8)                             :: sxa, sxb, v1xa, v1xb, v2xa, &
                                                v2xb, w1

! ==--------------------------------------------------------------==
! ..Exchange

    IF (rho.LE.small) THEN
       sx  = 0.0_real_8
       v1x = 0.0_real_8
       v2x = 0.0_real_8
    ELSE
       IF (func1%mgcx == mgcx_is_skipped) THEN
          sx  = 0.0_real_8
          v1x = 0.0_real_8
          v2x = 0.0_real_8
       ELSEIF (func1%mgcx == mgcx_is_becke88) THEN
          CALL becke88(func2%bbeta,rho,grho,sx,v1x,v2x)
       ELSEIF (func1%mgcx == mgcx_is_ggax) THEN
          CALL ggax(rho,grho,sx,v1x,v2x)
       ELSEIF (func1%mgcx == mgcx_is_pbex) THEN
          CALL pbex(rho,grho,sx,v1x,v2x)
       ELSEIF (func1%mgcx == mgcx_is_revpbex) THEN
          CALL revpbex(rho,grho,sx,v1x,v2x)
       ELSEIF (func1%mgcx == mgcx_is_hcth.AND.func1%mgcc == mgcc_is_hse) THEN
          CALL hcth(rho,grho,sx,v1x,v2x)! x&c
          sc=0.0_real_8
          v1c=0.0_real_8
          v2c=0.0_real_8
       ELSEIF (func1%mgcx == mgcx_is_optx) THEN
          CALL optx(rho,grho,sx,v1x,v2x)
       ELSEIF (func1%mgcx == mgcx_is_ox) THEN
          CALL becke88(func2%bbeta,rho,grho,sxa,v1xa,v2xa)
          CALL ggax(rho,grho,sxb,v1xb,v2xb)
          sx=0.722_real_8*sxa+0.347_real_8*sxb
          v1x=0.722_real_8*v1xa+0.347_real_8*v1xb
          v2x=0.722_real_8*v2xa+0.347_real_8*v2xb
       ELSEIF (func1%mgcx == mgcx_is_ox_hybrid) THEN
          CALL becke88(func2%bbeta,rho,grho,sxa,v1xa,v2xa)
          CALL ggax(rho,grho,sxb,v1xb,v2xb)
          sx=0.542_real_8*sxa+0.167_real_8*sxb
          v1x=0.542_real_8*v1xa+0.167_real_8*v1xb
          v2x=0.542_real_8*v2xa+0.167_real_8*v2xb
       ELSEIF (func1%mgcx == mgcx_is_pbesx) THEN
          CALL pbesx(rho,grho,sx,v1x,v2x)
       ELSEIF (func1%mgcx == mgcx_is_dfr_zpbex) THEN
          CALL stopgm(procedureN,'you need the DFRepository.F',&
               __LINE__,__FILE__)
       ELSEIF (func1%mgcx == mgcx_is_dfr_xpbex_hybrid) THEN
          CALL stopgm(procedureN,'you need the DFRepository.F',&
               __LINE__,__FILE__)
       ELSEIF (func1%mgcx == mgcx_is_dfr_xpbex) THEN
          CALL stopgm(procedureN,'you need the DFRepository.F',&
               __LINE__,__FILE__)
       ELSE
          ! vw MFXCX should be valide at this point!
          CALL stopgm(procedureN,"exchange functional "//&
               "not implemented",& 
               __LINE__,__FILE__)
       ENDIF
    ENDIF
    ! ..Correlation
    IF (rho.LE.small) THEN
       sc  = 0.0_real_8
       v1c = 0.0_real_8
       v2c = 0.0_real_8
    ELSE
       IF (func1%mgcc == mgcc_is_skipped) THEN
          sc  = 0.0_real_8
          v1c = 0.0_real_8
          v2c = 0.0_real_8
       ELSEIF (func1%mgcx == mgcx_is_hcth.AND.func1%mgcc == mgcc_is_hse) THEN
          ! HCTH
          sc  = 0.0_real_8
          v1c = 0.0_real_8
          v2c = 0.0_real_8
       ELSEIF (func1%mgcc == mgcc_is_perdew86) THEN
          CALL perdew86(rho,grho,sc,v1c,v2c)
       ELSEIF (func1%mgcc == mgcc_is_lyp) THEN
          CALL glyp(rho,grho,sc,v1c,v2c)
       ELSEIF (func1%mgcc == mgcc_is_ggac) THEN
          CALL ggac(rho,grho,sc,v1c,v2c)
       ELSEIF (func1%mgcc == mgcc_is_pbec) THEN
          w1=1._real_8
          CALL pbec(rho,grho,w1,sc,v1c,v2c)
       ELSEIF (func1%mgcc == mgcc_is_optc) THEN
          w1=0.74_real_8
          CALL pbec(rho,grho,w1,sc,v1c,v2c)
       ELSEIF (func1%mgcc == mgcc_is_pbesc) THEN
          w1=1.0_real_8
          CALL pbesc(rho,grho,w1,sc,v1c,v2c)
       ELSEIF (func1%mgcc == mgcc_is_dfr_zpbec) THEN
          CALL stopgm(procedureN,'you need the DFRepository.F',&
               __LINE__,__FILE__)
       ELSE
          ! vw MFCC should be valide at this point!
          CALL stopgm(procedureN,"correlation functional "//&
               "not implemented",& 
               __LINE__,__FILE__)
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gcxc
  ! ==================================================================
  SUBROUTINE slaterx(rho,ex,vx,alpha)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rho, ex, vx, alpha

    REAL(real_8), PARAMETER :: f1 = -1.10783814957303361_real_8 , &
      f43 = 4._real_8/3._real_8 , small = 1.e-10_real_8 , &
      third = 1._real_8/3._real_8

    REAL(real_8)                             :: rs

! ==--------------------------------------------------------------==

    IF (rho.LE.small) THEN
       ex = 0.0_real_8
       vx = 0.0_real_8
    ELSE
       rs = rho**third
       ex = f1*alpha*rs
       vx = f43*f1*alpha*rs
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE slaterx
  ! ==================================================================
  SUBROUTINE pz(rs,epz,vpz,iflg)
    ! ==--------------------------------------------------------------==
    ! ==  J.P. PERDEW AND ALEX ZUNGER PRB 23, 5048 (1981)             ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rs, epz, vpz
    INTEGER                                  :: iflg

    REAL(real_8), PARAMETER :: a = 0.0311_real_8, b = -0.048_real_8, &
      b1 = 1.0529_real_8, b2 = 0.3334_real_8 , c = 0.0020_real_8, &
      d = -0.0116_real_8, gc = -0.1423_real_8

    REAL(real_8)                             :: dox, ox, rs1, rs2, xln

! ==--------------------------------------------------------------==

    IF (iflg.EQ.1) THEN
       ! ..High density formula
       xln=LOG(rs)
       epz=a*xln+b+c*rs*xln+d*rs
       vpz=a*xln+(b-a/3._real_8)+2._real_8/3._real_8*c*rs*xln+&
            (2._real_8*d-c)/3._real_8*rs
    ELSEIF (iflg.EQ.2) THEN
       ! ..Interpolation formula
       rs1=SQRT(rs)
       rs2=rs
       ox=1._real_8+b1*rs1+b2*rs2
       dox=1._real_8+7._real_8/6._real_8*b1*rs1+4._real_8/3._real_8*b2*rs2
       epz=gc/ox
       vpz=epz*dox/ox
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE pz
  ! ==================================================================
  SUBROUTINE vwn(rs,evwn,vvwn)
    ! ==--------------------------------------------------------------==
    ! ==  S.H VOSKO, L.WILK, AND M. NUSAIR,                           ==
    ! ==                 CAN. J. PHYS. 581200  (1980)                ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rs, evwn, vvwn

    REAL(real_8), PARAMETER :: a = 0.0310907_real_8, b = 3.72744_real_8, &
      c = 12.9352_real_8, two = 2.0_real_8 , x0 = -0.10498_real_8 

    REAL(real_8)                             :: f1, f2, f3, fx, q, qx, ttqq, &
                                                txpb, x

! ==--------------------------------------------------------------==

    q  = SQRT(4._real_8*c-b*b)
    f1 = two*b/q
    f2 = b*x0/(x0*x0+b*x0+c)
    f3 = two*(two*x0+b)/q
    x  = SQRT(rs)
    fx = x*x+b*x+c
    qx = ATAN(q/(two*x+b))
    evwn=a*(LOG(rs/fx)+f1*qx-f2*(LOG((x-x0)**2/fx)+f3*qx))
    txpb=two*x+b
    ttqq=txpb*txpb+q*q
    vvwn=evwn - x*a/6._real_8*(two/x-txpb/fx-4._real_8*b/ttqq-f2*(two/(x-x0)&
         -txpb/fx-4._real_8*(two*x0+b)/ttqq))
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE vwn
  ! ==================================================================
  SUBROUTINE lyp(rho,elyp,vlyp)
    ! ==--------------------------------------------------------------==
    ! ==  C. LEE, W. YANG, AND R.G. PARR, PRB 37, 785 (1988)          ==
    ! ==  THIS IS ONLY THE LDA PART                                   ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rho, elyp, vlyp

    REAL(real_8), PARAMETER :: a = 0.04918_real_8, b = 0.132_real_8, &
      c = 0.2533_real_8, cf = 2.87123400018819108_real_8 , d = 0.349_real_8 

    REAL(real_8)                             :: ecrs, ox, rs

! ==--------------------------------------------------------------==

    rs=rho**(-1._real_8/3._real_8)
    ecrs=b*cf*EXP(-c*rs)
    ox=1._real_8/(1._real_8+d*rs)
    elyp=-a*ox*(1._real_8+ecrs)
    vlyp=elyp-rs/3._real_8*a*ox*(d*ox+ecrs*(d*ox+c))
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE lyp
  ! ==================================================================
  SUBROUTINE pw(rs,epwc,vpwc,iflg)
    ! ==--------------------------------------------------------------==
    ! ==  J.P. PERDEW AND YUE WANG PRB 45, 13244 (1992)               ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rs, epwc, vpwc
    INTEGER                                  :: iflg

    REAL(real_8), PARAMETER :: a = 0.031091_real_8, a1 = 0.21370_real_8, &
      b1 = 7.5957_real_8, b2 = 3.5876_real_8, b3 = 1.6382_real_8, &
      b4 = 0.49294_real_8, c0 = a, c1 = 0.046644_real_8, c2 = 0.00664_real_8, &
      c3 = 0.01043_real_8, d0 = 0.4335_real_8, d1 = 1.4408_real_8 

    REAL(real_8)                             :: dom, olog, om, rs1, rs2, rs3, &
                                                rs4, xln

! ==--------------------------------------------------------------==

    epwc=0.0_real_8
    vpwc=0.0_real_8
    IF (iflg.EQ.1) THEN
       ! ..High density formula
       xln=LOG(rs)
       epwc=c0*xln-c1+c2*rs*xln-c3*rs
       vpwc=c0*xln-(c1+c0/3._real_8)+2._real_8/3._real_8*c2*rs*xln-&
            (2._real_8*c3+c2)/3._real_8*rs
    ELSEIF (iflg.EQ.2) THEN
       ! ..Interpolation formula
       rs1=SQRT(rs)
       rs2=rs
       rs3=rs2*rs1
       rs4=rs2*rs2
       om=2._real_8*a*(b1*rs1+b2*rs2+b3*rs3+b4*rs4)
       dom=2._real_8*a*(0.5_real_8*b1*rs1+b2*rs2+1.5_real_8*b3*rs3+2._real_8*b4*rs4)
       olog=LOG(1._real_8+1.0_real_8/om)
       epwc=-2._real_8*a*(1._real_8+a1*rs)*olog
       vpwc=-2._real_8*a*(1._real_8+2._real_8/3._real_8*a1*rs)*olog&
            -2._real_8/3._real_8*a*(1._real_8+a1*rs)*dom/(om*(om+1._real_8))
    ELSEIF (iflg.EQ.3) THEN
       ! ..Low density formula
       epwc=-d0/rs+d1/rs**1.5_real_8
       vpwc=-4._real_8/3._real_8*d0/rs+1.5_real_8*d1/rs**1.5_real_8
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE pw
  ! ==================================================================
  SUBROUTINE wigner(rho,exc,fxc)
    REAL(real_8)                             :: rho, exc, fxc

    REAL(real_8)                             :: rh, x

    rh=rho
    x=rh**0.33333333333333333_real_8
    fxc=-x*((0.943656_real_8+8.8963_real_8*x)/(1.0_real_8+12.57_real_8*x)**2)
    exc=-0.738_real_8*x*(0.959_real_8/(1.0_real_8+12.57_real_8*x))
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE wigner
  ! ==================================================================
  SUBROUTINE hedin(rho,ecp,fcp)
    REAL(real_8)                             :: rho, ecp, fcp

    REAL(real_8)                             :: aln, rh, rsm1, x

! mb-ike   VARIABLES with SAVE attribute hinder the vectorization
! mb-ike      SAVE RH
! mb-ike      IF(RH .EQ. 0.0_real_8) RETURN

    rh=rho
    rsm1=0.62035049_real_8*rh**(0.3333333333333333_real_8)
    aln=LOG(1.0_real_8 + 21.0_real_8*rsm1)
    x=21.0_real_8/rsm1
    ecp = aln+(x**3*aln-x*x)+x/2.0_real_8-1.0_real_8/3.0_real_8
    ecp = -0.0225_real_8*ecp
    fcp = -0.0225_real_8*aln
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE hedin
  ! ==================================================================
  SUBROUTINE obpz(rs,epz,vpz,iflg)
    ! ==--------------------------------------------------------------==
    ! ==  G.ORTIZ AND P. BALLONE PRB 50, 1391 (1994)                  ==
    ! ==  PERDEW-ZUNGER FORMULA                                       ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rs, epz, vpz
    INTEGER                                  :: iflg

    REAL(real_8), PARAMETER :: a = 0.031091_real_8, b = -0.046644_real_8, &
      b1 = 0.56371_real_8, b2 = 0.27358_real_8 , c = 0.00419_real_8, &
      d = -0.00983_real_8, gc = -0.103756_real_8

    REAL(real_8)                             :: dox, ox, rs1, rs2, xln

! ==--------------------------------------------------------------==

    IF (iflg.EQ.1) THEN
       ! ..High density formula
       xln=LOG(rs)
       epz=a*xln+b+c*rs*xln+d*rs
       vpz=a*xln+(b-a/3._real_8)+2._real_8/3._real_8*c*rs*xln+&
            (2._real_8*d-c)/3._real_8*rs
    ELSEIF (iflg.EQ.2) THEN
       ! ..Interpolation formula
       rs1=SQRT(rs)
       rs2=rs
       ox=1._real_8+b1*rs1+b2*rs2
       dox=1._real_8+7._real_8/6._real_8*b1*rs1+4._real_8/3._real_8*b2*rs2
       epz=gc/ox
       vpz=epz*dox/ox
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE obpz
  ! ==================================================================
  SUBROUTINE obpw(rs,epwc,vpwc,iflg)
    ! ==--------------------------------------------------------------==
    ! ==  G.ORTIZ AND P. BALLONE PRB 50, 1391 (1994)                  ==
    ! ==  PERDEW-WANG FORMULA                                         ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rs, epwc, vpwc
    INTEGER                                  :: iflg

    REAL(real_8), PARAMETER :: a = 0.031091_real_8, a1 = 0.026481_real_8, &
      b1 = 7.5957_real_8, b2 = 3.5876_real_8, b3 = -0.46647_real_8, &
      b4 = 0.13354_real_8, c0 = a, c1 = 0.046644_real_8, c2 = 0.00664_real_8, &
      c3 = 0.01043_real_8, d0 = 0.4335_real_8, d1 = 1.4408_real_8 

    REAL(real_8)                             :: dom, olog, om, rs1, rs2, rs3, &
                                                rs4, xln

! ==--------------------------------------------------------------==

    epwc=0.0_real_8
    vpwc=0.0_real_8
    IF (iflg.EQ.1) THEN
       ! ..High density formula
       xln=LOG(rs)
       epwc=c0*xln-c1+c2*rs*xln-c3*rs
       vpwc=c0*xln-(c1+c0/3._real_8)+2._real_8/3._real_8*c2*rs*xln-&
            (2._real_8*c3+c2)/3._real_8*rs
    ELSEIF (iflg.EQ.2) THEN
       ! ..Interpolation formula
       rs1=SQRT(rs)
       rs2=rs
       rs3=rs2*rs1
       rs4=rs2*rs2
       om=2._real_8*a*(b1*rs1+b2*rs2+b3*rs3+b4*rs4)
       dom=2._real_8*a*(0.5_real_8*b1*rs1+b2*rs2+1.5_real_8*b3*rs3+2._real_8*b4*rs4)
       olog=LOG(1._real_8+1.0_real_8/om)
       epwc=-2._real_8*a*(1.0_real_8+a1*rs)*olog
       vpwc=-2._real_8*a*(1._real_8+2._real_8/3._real_8*a1*rs)*olog&
            -2._real_8/3._real_8*a*(1._real_8+a1*rs)*dom/(om*(om+1._real_8))
    ELSEIF (iflg.EQ.3) THEN
       ! ..Low density formula
       epwc=-d0/rs+d1/rs**1.5_real_8
       vpwc=-4._real_8/3._real_8*d0/rs+1.5_real_8*d1/rs**1.5_real_8
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE obpw
  ! ==================================================================
  SUBROUTINE pade(rs,ec,vc)
    ! ==--------------------------------------------------------------==
    ! ==  PADE APPROXIMATION                                          ==
    ! ==  S. GOEDECKER, M. TETER, J. HUTTER, PRB in press             ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rs, ec, vc

    REAL(real_8), PARAMETER :: a0 = 0.4581652932831429_real_8, &
      a1 = 2.217058676663745_real_8, a2 = 0.7405551735357053_real_8, &
      a3 = 0.01968227878617998_real_8 , b1 = 1.0000000000000000_real_8, &
      b2 = 4.504130959426697_real_8, b3 = 1.110667363742916_real_8, &
      b4 = 0.02359291751427506_real_8 , o3 = 1._real_8/3._real_8 

    REAL(real_8)                             :: bot, dbot, dtop, top

! ==--------------------------------------------------------------==

    top=a0+rs*(a1+rs*(a2+rs*a3))
    dtop=a1+rs*(2._real_8*a2+3._real_8*a3*rs)
    bot=rs*(b1+rs*(b2+rs*(b3+rs*b4)))
    dbot=b1+rs*(2._real_8*b2+rs*(3._real_8*b3+rs*4._real_8*b4))
    ec=-top/bot
    vc=ec+rs*o3*(dtop/bot-top*dbot/(bot*bot))
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE pade
  ! ==================================================================
  SUBROUTINE becke88(b1,rho,grho,sx,v1x,v2x)
    ! ==--------------------------------------------------------------==
    ! BECKE EXCHANGE: PRA 38, 3098 (1988)
    REAL(real_8)                             :: b1, rho, grho, sx, v1x, v2x

    REAL(real_8), PARAMETER                  :: ob3 = 1._real_8/3._real_8 

    REAL(real_8)                             :: a, aa, br1, br2, br4, dd, &
                                                dd2, ee, sa2b8, shm1, two13, &
                                                xs, xs2

! ==--------------------------------------------------------------==

    two13 = 2.0_real_8**(1._real_8/3._real_8)
    aa    = grho
    a     = SQRT(aa)
    br1   = rho**ob3
    br2   = br1*br1
    br4   = br2*br2
    xs    = two13*a/br4
    xs2   = xs*xs
    sa2b8 = SQRT(1.0_real_8+xs2)
    shm1  = LOG(xs+sa2b8)
    dd    = 1.0_real_8 + 6.0_real_8*b1*xs*shm1
    dd2   = dd*dd
    ee    = 6.0_real_8*b1*xs2/sa2b8 - 1._real_8
    sx    = two13*aa/br4*(-b1/dd)
    v1x   = -(4._real_8/3._real_8)/two13*xs2*b1*br1*ee/dd2
    v2x   = two13*b1*(ee-dd)/(br4*dd2)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE becke88
  ! ==================================================================
  SUBROUTINE ggax(rho,grho,sx,v1x,v2x)
    ! ==--------------------------------------------------------------==
    ! J.P.PERDEW ET AL. PRB 466671 (1992)
    REAL(real_8)                             :: rho, grho, sx, v1x, v2x

    REAL(real_8), PARAMETER :: f1 = 0.19645_real_8, f2 = 7.7956_real_8, &
      f3 = 0.2743_real_8, f4 = 0.1508_real_8, f5 = 0.004e9_real_8 

    REAL(real_8)                             :: a, aa, as, bs, das, dbs, dls, &
                                                exps, fp1, fp2, rr, s, s2, &
                                                s3, s4, sa2b8, shm1

! ==--------------------------------------------------------------==

    fp1   = -3._real_8/(16._real_8*pi)*(3._real_8*pi*pi)**(-1._real_8/3._real_8)
    fp2   = 0.5_real_8*(3._real_8*pi*pi)**(-1._real_8/3._real_8)
    aa    = grho
    a     = SQRT(aa)
    rr    = rho**(-4._real_8/3._real_8)
    s     = fp2*a*rr
    s2    = s*s
    s3    = s2*s
    s4    = s2*s2
    exps  = f4*EXP(-100._real_8*s2)
    as    = f3-exps-f5*s2
    sa2b8 = SQRT(1.0_real_8+f2*f2*s2)
    shm1  = LOG(f2*s+sa2b8)
    bs    = 1._real_8+f1*s*shm1+f5*s4
    das   = 200._real_8*s*exps-2._real_8*s*f5
    dbs   = f1*(shm1+f2*s/sa2b8)+4._real_8*f5*s3
    dls   = (das/as-dbs/bs)
    sx    = fp1*aa*rr*as/bs
    v1x   = -4._real_8/3._real_8*sx/rho*(1._real_8+s*dls)
    v2x   = fp1*rr*as/bs*(2._real_8+s*dls)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE ggax
  ! ==================================================================
  SUBROUTINE pbesx(rho,grho,sx,v1x,v2x)
    ! ==--------------------------------------------------------------==
    ! PBESol functional
    REAL(real_8)                             :: rho, grho, sx, v1x, v2x

    REAL(real_8), PARAMETER :: ax = -0.738558766382022406_real_8, &
      uk = 0.8040_real_8, um = 0.123456790123456789_real_8, ul = um/uk , &
      us = 0.161620459673995492_real_8

    REAL(real_8)                             :: aa, dfx, ex, fx, po, rr, s2

! ==--------------------------------------------------------------==

    aa    = grho
    rr    = rho**(-4._real_8/3._real_8)
    ex    = ax/rr
    s2    = aa*rr*rr*us*us
    po    = 1._real_8/(1._real_8 + ul*s2)
    fx    = uk-uk*po
    sx    = ex*fx
    dfx   = 2._real_8*uk*ul*po*po
    v1x   = 1.33333333333333_real_8*ax*rho**0.333333333333_real_8*(fx-s2*dfx)
    v2x   = ex*dfx*(us*rr)**2
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE pbesx
  ! ==================================================================
  SUBROUTINE pbex(rho,grho,sx,v1x,v2x)
    ! ==--------------------------------------------------------------==
    ! J.P.PERDEW ET AL. PRL 773865 (1996)
    REAL(real_8)                             :: rho, grho, sx, v1x, v2x

    REAL(real_8), PARAMETER :: ax = -0.738558766382022406_real_8, &
      uk = 0.8040_real_8, um = 0.2195149727645171_real_8, ul = um/uk , &
      us = 0.161620459673995492_real_8

    REAL(real_8)                             :: aa, dfx, ex, fx, po, rr, s2

! ==--------------------------------------------------------------==

    aa    = grho
    rr    = rho**(-4._real_8/3._real_8)
    ex    = ax/rr
    s2    = aa*rr*rr*us*us
    po    = 1._real_8/(1._real_8 + ul*s2)
    fx    = uk-uk*po
    sx    = ex*fx
    dfx   = 2._real_8*uk*ul*po*po
    v1x   = 1.33333333333333_real_8*ax*rho**0.333333333333_real_8*(fx-s2*dfx)
    v2x   = ex*dfx*(us*rr)**2
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE pbex
  ! ==================================================================
  SUBROUTINE rks_x_zpbe(r0,rho,grho,sx,v1x,v2x)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: r0, rho, grho, sx, v1x, v2x

    CHARACTER(len=*), PARAMETER              :: procedureN = 'rks_x_zpbe'
    REAL(real_8), PARAMETER                  :: kappa = 0.8040_real_8 , &
                                                mu = 0.2195149727645171_real_8

! ==--------------------------------------------------------------==

    IF (r0.LT.1.0e-5_real_8) CALL stopgm(procedureN,'RO too small',& 
         __LINE__,__FILE__)
    sx = -0.00806288360829988_real_8*(kappa - kappa/(0.0261211729852336_real_8*&
         GRHO*mu/(RHO**(8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8)**3*(&
         9.5707800006273_real_8*R0**2*RHO**(2.0_real_8/3.0_real_8)/(kappa - kappa/(&
         0.0261211729852336_real_8*GRHO*mu/(RHO**(8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8)&
         + 1.0_real_8) - 3.09366772628014_real_8*R0*RHO**(1.0_real_8/3.0_real_8)*SIN(&
         6.18733545256027_real_8*R0*RHO**(1.0_real_8/3.0_real_8)/SQRT(kappa - kappa/(&
         0.0261211729852336_real_8*GRHO*mu/(RHO**(8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8)&
         + 1.0_real_8))/SQRT(kappa - kappa/(0.0261211729852336_real_8*GRHO*mu/(RHO**&
         (8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8) - 0.5_real_8*COS(&
         6.18733545256027_real_8*R0*RHO**(1.0_real_8/3.0_real_8)/SQRT(kappa - kappa/(&
         0.0261211729852336_real_8*GRHO*mu/(RHO**(8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8)&
         + 1.0_real_8)) + 0.5_real_8)/R0**4
    v1x = 0.00168489581993764_real_8*grho*mu*(kappa - kappa/(&
         0.0261211729852336_real_8*GRHO*mu/(RHO**(8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8)&
         + 1.0_real_8)**2*(9.5707800006273_real_8*R0**2*RHO**(2.0_real_8/3.0_real_8)/(kappa -&
         kappa/(0.0261211729852336_real_8*GRHO*mu/(RHO**(8.0_real_8/3.0_real_8)*kappa) +&
         1.0_real_8) + 1.0_real_8) - 3.09366772628014_real_8*R0*RHO**(1.0_real_8/3.0_real_8)*SIN(&
         6.18733545256027_real_8*R0*RHO**(1.0_real_8/3.0_real_8)/SQRT(kappa - kappa/(&
         0.0261211729852336_real_8*GRHO*mu/(RHO**(8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8)&
         + 1.0_real_8))/SQRT(kappa - kappa/(0.0261211729852336_real_8*GRHO*mu/(RHO**&
         (8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8) - 0.5_real_8*COS(&
         6.18733545256027_real_8*R0*RHO**(1.0_real_8/3.0_real_8)/SQRT(kappa - kappa/(&
         0.0261211729852336_real_8*GRHO*mu/(RHO**(8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8)&
         + 1.0_real_8)) + 0.5_real_8)/(R0**4*RHO**(11.0_real_8/3.0_real_8)*(&
         0.0261211729852336_real_8*GRHO*mu/(RHO**(8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8)&
         **2) - 0.00806288360829988_real_8*(kappa - kappa/(0.0261211729852336_real_8&
         *GRHO*mu/(RHO**(8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8)**3*(&
         0.666666666666667_real_8*GRHO*R0**2*mu/(RHO**3*(0.0261211729852336_real_8*&
         GRHO*mu/(RHO**(8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8)**2*(kappa - kappa/(&
         0.0261211729852336_real_8*GRHO*mu/(RHO**(8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8)&
         + 1.0_real_8)**2) - 0.107746973115997_real_8*GRHO*R0*mu*SIN(&
         6.18733545256027_real_8*R0*RHO**(1.0_real_8/3.0_real_8)/SQRT(kappa - kappa/(&
         0.0261211729852336_real_8*GRHO*mu/(RHO**(8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8)&
         + 1.0_real_8))/(RHO**(10.0_real_8/3.0_real_8)*(0.0261211729852336_real_8*GRHO*mu/(RHO&
         **(8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8)**2*(kappa - kappa/(&
         0.0261211729852336_real_8*GRHO*mu/(RHO**(8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8)&
         + 1.0_real_8)**(3.0_real_8/2.0_real_8)) + 6.3805200004182_real_8*R0**2/(RHO**(1.0_real_8/&
         3.0_real_8)*(kappa - kappa/(0.0261211729852336_real_8*GRHO*mu/(RHO**(8.0_real_8/&
         3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8)) - 3.09366772628014_real_8*R0*RHO**(&
         1.0_real_8/3.0_real_8)*(0.215493946231994_real_8*GRHO*R0*mu/(RHO**(10.0_real_8/3.0_real_8)&
         *(0.0261211729852336_real_8*GRHO*mu/(RHO**(8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8&
         )**2*(kappa - kappa/(0.0261211729852336_real_8*GRHO*mu/(RHO**(8.0_real_8/&
         3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8)**(3.0_real_8/2.0_real_8)) +&
         2.06244515085342_real_8*R0/(RHO**(2.0_real_8/3.0_real_8)*SQRT(kappa - kappa/(&
         0.0261211729852336_real_8*GRHO*mu/(RHO**(8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8)&
         + 1.0_real_8)))*COS(6.18733545256027_real_8*R0*RHO**(1.0_real_8/3.0_real_8)/SQRT(&
         kappa - kappa/(0.0261211729852336_real_8*GRHO*mu/(RHO**(8.0_real_8/3.0_real_8)*&
         kappa) + 1.0_real_8) + 1.0_real_8))/SQRT(kappa - kappa/(&
         0.0261211729852336_real_8*GRHO*mu/(RHO**(8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8)&
         + 1.0_real_8) - 1.03122257542671_real_8*R0*SIN(6.18733545256027_real_8*R0*RHO**(&
         1.0_real_8/3.0_real_8)/SQRT(kappa - kappa/(0.0261211729852336_real_8*GRHO*mu/(&
         RHO**(8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8))/(RHO**(2.0_real_8/3.0_real_8)*&
         SQRT(kappa - kappa/(0.0261211729852336_real_8*GRHO*mu/(RHO**(8.0_real_8/&
         3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8)) + 0.5_real_8*(0.215493946231994_real_8*&
         GRHO*R0*mu/(RHO**(10.0_real_8/3.0_real_8)*(0.0261211729852336_real_8*GRHO*mu/(&
         RHO**(8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8)**2*(kappa - kappa/(&
         0.0261211729852336_real_8*GRHO*mu/(RHO**(8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8)&
         + 1.0_real_8)**(3.0_real_8/2.0_real_8)) + 2.06244515085342_real_8*R0/(RHO**(2.0_real_8/&
         3.0_real_8)*SQRT(kappa - kappa/(0.0261211729852336_real_8*GRHO*mu/(RHO**(&
         8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8)))*SIN(6.18733545256027_real_8*R0&
         *RHO**(1.0_real_8/3.0_real_8)/SQRT(kappa - kappa/(0.0261211729852336_real_8*GRHO&
         *mu/(RHO**(8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8)))/R0**4
    v2x = -0.0161257672165998_real_8*(0.25_real_8*r0**2*mu*COS(&
         6.18733545256027_real_8*R0*RHO**(1.0_real_8/3.0_real_8)/SQRT(kappa - kappa/(&
         0.0261211729852336_real_8*GRHO*mu/(RHO**(8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8)&
         + 1.0_real_8))/(RHO**2*(0.0261211729852336_real_8*GRHO*mu/(RHO**(8.0_real_8/&
         3.0_real_8)*kappa) + 1.0_real_8)**2*(kappa - kappa/(0.0261211729852336_real_8*&
         GRHO*mu/(RHO**(8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8)**2) - 0.25_real_8&
         *R0**2*mu/(RHO**2*(0.0261211729852336_real_8*GRHO*mu/(RHO**(8.0_real_8/&
         3.0_real_8)*kappa) + 1.0_real_8)**2*(kappa - kappa/(0.0261211729852336_real_8*&
         GRHO*mu/(RHO**(8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8)**2))*(kappa&
         - kappa/(0.0261211729852336_real_8*GRHO*mu/(RHO**(8.0_real_8/3.0_real_8)*kappa)&
         + 1.0_real_8) + 1.0_real_8)**3/R0**4 - 0.00126367186495323_real_8*mu*(kappa -&
         kappa/(0.0261211729852336_real_8*GRHO*mu/(RHO**(8.0_real_8/3.0_real_8)*kappa) +&
         1.0_real_8) + 1.0_real_8)**2*(9.5707800006273_real_8*R0**2*RHO**(2.0_real_8/3.0_real_8)/(&
         kappa - kappa/(0.0261211729852336_real_8*GRHO*mu/(RHO**(8.0_real_8/3.0_real_8)*&
         kappa) + 1.0_real_8) + 1.0_real_8) - 3.09366772628014_real_8*R0*RHO**(1.0_real_8/&
         3.0_real_8)*SIN(6.18733545256027_real_8*R0*RHO**(1.0_real_8/3.0_real_8)/SQRT(kappa -&
         kappa/(0.0261211729852336_real_8*GRHO*mu/(RHO**(8.0_real_8/3.0_real_8)*kappa) +&
         1.0_real_8) + 1.0_real_8))/SQRT(kappa - kappa/(0.0261211729852336_real_8*GRHO*mu&
         /(RHO**(8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8) - 0.5_real_8*COS(&
         6.18733545256027_real_8*R0*RHO**(1.0_real_8/3.0_real_8)/SQRT(kappa - kappa/(&
         0.0261211729852336_real_8*GRHO*mu/(RHO**(8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8)&
         + 1.0_real_8)) + 0.5_real_8)/(R0**4*RHO**(8.0_real_8/3.0_real_8)*(&
         0.0261211729852336_real_8*GRHO*mu/(RHO**(8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8)&
         **2)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rks_x_zpbe

  ! ==================================================================
  SUBROUTINE uks_x_zpbe(r0,rhoa,rhob,grhoaa,grhobb,&
       sx,v1xa,v1xb,v2xa,v2xb,v2xab)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: r0, rhoa, rhob, grhoaa, &
                                                grhobb, sx, v1xa, v1xb, v2xa, &
                                                v2xb, v2xab

    CHARACTER(len=*), PARAMETER              :: procedureN = 'uks_x_zpbe'
    REAL(real_8), PARAMETER                  :: kappa = 0.8040_real_8, mu = &
                                                0.2195149727645171_real_8, &
                                                tol = 1.0e-20_real_8 

! ==--------------------------------------------------------------==

    IF (r0.LT.1.0e-5_real_8) CALL stopgm(procedureN,'RO too small',& 
         __LINE__,__FILE__)
    IF (rhoa.GT.tol) THEN
       IF (rhob.GT.tol) THEN
          sx = -0.00403144180414994_real_8*(kappa - kappa/(0.0164553078460206_real_8*&
               GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8)**3*(&
               15.192666241152_real_8*R0**2*RHOA**(2.0_real_8/3.0_real_8)/(kappa - kappa/(&
               0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8) - 3.89777708972075_real_8*R0*RHOA**(1.0_real_8/3.0_real_8)*SIN(&
               7.79555417944151_real_8*R0*RHOA**(1.0_real_8/3.0_real_8)/SQRT(kappa - kappa/(&
               0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8))/SQRT(kappa - kappa/(0.0164553078460206_real_8*GRHOAA*&
               mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8) - 0.5_real_8*COS(&
               7.79555417944151_real_8*R0*RHOA**(1.0_real_8/3.0_real_8)/SQRT(kappa - kappa/(&
               0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8)) + 0.5_real_8)/R0**4 - 0.00403144180414994_real_8*(kappa -&
               kappa/(0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8)*kappa&
               ) + 1.0_real_8) + 1.0_real_8)**3*(15.192666241152_real_8*R0**2*RHOB**(2.0_real_8/&
               3.0_real_8)/(kappa - kappa/(0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(&
               8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8) - 3.89777708972075_real_8*R0*&
               RHOB**(1.0_real_8/3.0_real_8)*SIN(7.79555417944151_real_8*R0*RHOB**(1.0_real_8/3.0_real_8)&
               /SQRT(kappa - kappa/(0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8&
               /3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8))/SQRT(kappa - kappa/(&
               0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8) - 0.5_real_8*COS(7.79555417944151_real_8*R0*RHOB**(1.0_real_8/&
               3.0_real_8)/SQRT(kappa - kappa/(0.0164553078460206_real_8*GRHOBB*mu/(RHOB**&
               (8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8)) + 0.5_real_8)/R0**4
          v1xa = 0.00053070892760483_real_8*grhoaa*mu*(kappa - kappa/(&
               0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8)**2*(15.192666241152_real_8*R0**2*RHOA**(2.0_real_8/3.0_real_8)/(&
               kappa - kappa/(0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8&
               )*kappa) + 1.0_real_8) + 1.0_real_8) - 3.89777708972075_real_8*R0*RHOA**(1.0_real_8/&
               3.0_real_8)*SIN(7.79555417944151_real_8*R0*RHOA**(1.0_real_8/3.0_real_8)/SQRT(kappa -&
               kappa/(0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa&
               ) + 1.0_real_8) + 1.0_real_8))/SQRT(kappa - kappa/(0.0164553078460206_real_8*&
               GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8) - 0.5_real_8*&
               COS(7.79555417944151_real_8*R0*RHOA**(1.0_real_8/3.0_real_8)/SQRT(kappa - kappa/&
               (0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8)) + 0.5_real_8)/(R0**4*RHOA**(11.0_real_8/3.0_real_8)*(&
               0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8)**2) - 0.00403144180414994_real_8*(kappa - kappa/(&
               0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8)**3*(0.666666666666667_real_8*GRHOAA*R0**2*mu/(RHOA**3*&
               (0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8)**2*(kappa - kappa/(0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(&
               8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8)**2) - 0.0855188292353615_real_8*&
               GRHOAA*R0*mu*SIN(7.79555417944151_real_8*R0*RHOA**(1.0_real_8/3.0_real_8)/SQRT(&
               kappa - kappa/(0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8&
               )*kappa) + 1.0_real_8) + 1.0_real_8))/(RHOA**(10.0_real_8/3.0_real_8)*(&
               0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8)**2*(kappa - kappa/(0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(&
               8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8)**(3.0_real_8/2.0_real_8)) +&
               10.128444160768_real_8*R0**2/(RHOA**(1.0_real_8/3.0_real_8)*(kappa - kappa/(&
               0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8)) - 3.89777708972075_real_8*R0*RHOA**(1.0_real_8/3.0_real_8)*(&
               0.171037658470723_real_8*GRHOAA*R0*mu/(RHOA**(10.0_real_8/3.0_real_8)*(&
               0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8)**2*(kappa - kappa/(0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(&
               8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8)**(3.0_real_8/2.0_real_8)) +&
               2.59851805981384_real_8*R0/(RHOA**(2.0_real_8/3.0_real_8)*SQRT(kappa - kappa/(&
               0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8)))*COS(7.79555417944151_real_8*R0*RHOA**(1.0_real_8/3.0_real_8)/&
               SQRT(kappa - kappa/(0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/&
               3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8))/SQRT(kappa - kappa/(&
               0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8) - 1.29925902990692_real_8*R0*SIN(7.79555417944151_real_8*R0&
               *RHOA**(1.0_real_8/3.0_real_8)/SQRT(kappa - kappa/(0.0164553078460206_real_8*&
               GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8))/(RHOA**(&
               2.0_real_8/3.0_real_8)*SQRT(kappa - kappa/(0.0164553078460206_real_8*GRHOAA*mu/(&
               RHOA**(8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8)) + 0.5_real_8*(&
               0.171037658470723_real_8*GRHOAA*R0*mu/(RHOA**(10.0_real_8/3.0_real_8)*(&
               0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8)**2*(kappa - kappa/(0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(&
               8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8)**(3.0_real_8/2.0_real_8)) +&
               2.59851805981384_real_8*R0/(RHOA**(2.0_real_8/3.0_real_8)*SQRT(kappa - kappa/(&
               0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8)))*SIN(7.79555417944151_real_8*R0*RHOA**(1.0_real_8/3.0_real_8)/&
               SQRT(kappa - kappa/(0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/&
               3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8)))/R0**4
          v1xb = 0.00053070892760483_real_8*grhobb*mu*(kappa - kappa/(&
               0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8)**2*(15.192666241152_real_8*R0**2*RHOB**(2.0_real_8/3.0_real_8)/(&
               kappa - kappa/(0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8&
               )*kappa) + 1.0_real_8) + 1.0_real_8) - 3.89777708972075_real_8*R0*RHOB**(1.0_real_8/&
               3.0_real_8)*SIN(7.79555417944151_real_8*R0*RHOB**(1.0_real_8/3.0_real_8)/SQRT(kappa -&
               kappa/(0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8)*kappa&
               ) + 1.0_real_8) + 1.0_real_8))/SQRT(kappa - kappa/(0.0164553078460206_real_8*&
               GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8) - 0.5_real_8*&
               COS(7.79555417944151_real_8*R0*RHOB**(1.0_real_8/3.0_real_8)/SQRT(kappa - kappa/&
               (0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8)) + 0.5_real_8)/(R0**4*RHOB**(11.0_real_8/3.0_real_8)*(&
               0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8)**2) - 0.00403144180414994_real_8*(kappa - kappa/(&
               0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8)**3*(0.666666666666667_real_8*GRHOBB*R0**2*mu/(RHOB**3*&
               (0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8)**2*(kappa - kappa/(0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(&
               8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8)**2) - 0.0855188292353615_real_8*&
               GRHOBB*R0*mu*SIN(7.79555417944151_real_8*R0*RHOB**(1.0_real_8/3.0_real_8)/SQRT(&
               kappa - kappa/(0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8&
               )*kappa) + 1.0_real_8) + 1.0_real_8))/(RHOB**(10.0_real_8/3.0_real_8)*(&
               0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8)**2*(kappa - kappa/(0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(&
               8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8)**(3.0_real_8/2.0_real_8)) +&
               10.128444160768_real_8*R0**2/(RHOB**(1.0_real_8/3.0_real_8)*(kappa - kappa/(&
               0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8)) - 3.89777708972075_real_8*R0*RHOB**(1.0_real_8/3.0_real_8)*(&
               0.171037658470723_real_8*GRHOBB*R0*mu/(RHOB**(10.0_real_8/3.0_real_8)*(&
               0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8)**2*(kappa - kappa/(0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(&
               8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8)**(3.0_real_8/2.0_real_8)) +&
               2.59851805981384_real_8*R0/(RHOB**(2.0_real_8/3.0_real_8)*SQRT(kappa - kappa/(&
               0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8)))*COS(7.79555417944151_real_8*R0*RHOB**(1.0_real_8/3.0_real_8)/&
               SQRT(kappa - kappa/(0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/&
               3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8))/SQRT(kappa - kappa/(&
               0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8) - 1.29925902990692_real_8*R0*SIN(7.79555417944151_real_8*R0&
               *RHOB**(1.0_real_8/3.0_real_8)/SQRT(kappa - kappa/(0.0164553078460206_real_8*&
               GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8))/(RHOB**(&
               2.0_real_8/3.0_real_8)*SQRT(kappa - kappa/(0.0164553078460206_real_8*GRHOBB*mu/(&
               RHOB**(8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8)) + 0.5_real_8*(&
               0.171037658470723_real_8*GRHOBB*R0*mu/(RHOB**(10.0_real_8/3.0_real_8)*(&
               0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8)**2*(kappa - kappa/(0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(&
               8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8)**(3.0_real_8/2.0_real_8)) +&
               2.59851805981384_real_8*R0/(RHOB**(2.0_real_8/3.0_real_8)*SQRT(kappa - kappa/(&
               0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8)))*SIN(7.79555417944151_real_8*R0*RHOB**(1.0_real_8/3.0_real_8)/&
               SQRT(kappa - kappa/(0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/&
               3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8)))/R0**4
          v2xa = -0.00806288360829988_real_8*(0.25_real_8*r0**2*mu*COS(&
               7.79555417944151_real_8*R0*RHOA**(1.0_real_8/3.0_real_8)/SQRT(kappa - kappa/(&
               0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8))/(RHOA**2*(0.0164553078460206_real_8*GRHOAA*mu/(RHOA**&
               (8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8)**2*(kappa - kappa/(&
               0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8)**2) - 0.25_real_8*R0**2*mu/(RHOA**2*(&
               0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8)**2*(kappa - kappa/(0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(&
               8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8)**2))*(kappa - kappa/(&
               0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8)**3/R0**4 - 0.000398031695703623_real_8*mu*(kappa -&
               kappa/(0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa&
               ) + 1.0_real_8) + 1.0_real_8)**2*(15.192666241152_real_8*R0**2*RHOA**(2.0_real_8/&
               3.0_real_8)/(kappa - kappa/(0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(&
               8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8) - 3.89777708972075_real_8*R0*&
               RHOA**(1.0_real_8/3.0_real_8)*SIN(7.79555417944151_real_8*R0*RHOA**(1.0_real_8/3.0_real_8)&
               /SQRT(kappa - kappa/(0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8&
               /3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8))/SQRT(kappa - kappa/(&
               0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8) - 0.5_real_8*COS(7.79555417944151_real_8*R0*RHOA**(1.0_real_8/&
               3.0_real_8)/SQRT(kappa - kappa/(0.0164553078460206_real_8*GRHOAA*mu/(RHOA**&
               (8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8)) + 0.5_real_8)/(R0**4*RHOA**(&
               8.0_real_8/3.0_real_8)*(0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)&
               *kappa) + 1.0_real_8)**2)
          v2xb = -0.00806288360829988_real_8*(0.25_real_8*r0**2*mu*COS(&
               7.79555417944151_real_8*R0*RHOB**(1.0_real_8/3.0_real_8)/SQRT(kappa - kappa/(&
               0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8))/(RHOB**2*(0.0164553078460206_real_8*GRHOBB*mu/(RHOB**&
               (8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8)**2*(kappa - kappa/(&
               0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8)**2) - 0.25_real_8*R0**2*mu/(RHOB**2*(&
               0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8)**2*(kappa - kappa/(0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(&
               8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8)**2))*(kappa - kappa/(&
               0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8)**3/R0**4 - 0.000398031695703623_real_8*mu*(kappa -&
               kappa/(0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8)*kappa&
               ) + 1.0_real_8) + 1.0_real_8)**2*(15.192666241152_real_8*R0**2*RHOB**(2.0_real_8/&
               3.0_real_8)/(kappa - kappa/(0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(&
               8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8) - 3.89777708972075_real_8*R0*&
               RHOB**(1.0_real_8/3.0_real_8)*SIN(7.79555417944151_real_8*R0*RHOB**(1.0_real_8/3.0_real_8)&
               /SQRT(kappa - kappa/(0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8&
               /3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8))/SQRT(kappa - kappa/(&
               0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8) - 0.5_real_8*COS(7.79555417944151_real_8*R0*RHOB**(1.0_real_8/&
               3.0_real_8)/SQRT(kappa - kappa/(0.0164553078460206_real_8*GRHOBB*mu/(RHOB**&
               (8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8)) + 0.5_real_8)/(R0**4*RHOB**(&
               8.0_real_8/3.0_real_8)*(0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8)&
               *kappa) + 1.0_real_8)**2)
          v2xab = 0.0_real_8
       ELSE
          sx = -0.00403144180414994_real_8*(kappa - kappa/(0.0164553078460206_real_8*&
               GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8)**3*(&
               15.192666241152_real_8*R0**2*RHOA**(2.0_real_8/3.0_real_8)/(kappa - kappa/(&
               0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8) - 3.89777708972075_real_8*R0*RHOA**(1.0_real_8/3.0_real_8)*SIN(&
               7.79555417944151_real_8*R0*RHOA**(1.0_real_8/3.0_real_8)/SQRT(kappa - kappa/(&
               0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8))/SQRT(kappa - kappa/(0.0164553078460206_real_8*GRHOAA*&
               mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8) - 0.5_real_8*COS(&
               7.79555417944151_real_8*R0*RHOA**(1.0_real_8/3.0_real_8)/SQRT(kappa - kappa/(&
               0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8)) + 0.5_real_8)/R0**4
          v1xa = 0.00053070892760483_real_8*grhoaa*mu*(kappa - kappa/(&
               0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8)**2*(15.192666241152_real_8*R0**2*RHOA**(2.0_real_8/3.0_real_8)/(&
               kappa - kappa/(0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8&
               )*kappa) + 1.0_real_8) + 1.0_real_8) - 3.89777708972075_real_8*R0*RHOA**(1.0_real_8/&
               3.0_real_8)*SIN(7.79555417944151_real_8*R0*RHOA**(1.0_real_8/3.0_real_8)/SQRT(kappa -&
               kappa/(0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa&
               ) + 1.0_real_8) + 1.0_real_8))/SQRT(kappa - kappa/(0.0164553078460206_real_8*&
               GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8) - 0.5_real_8*&
               COS(7.79555417944151_real_8*R0*RHOA**(1.0_real_8/3.0_real_8)/SQRT(kappa - kappa/&
               (0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8)) + 0.5_real_8)/(R0**4*RHOA**(11.0_real_8/3.0_real_8)*(&
               0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8)**2) - 0.00403144180414994_real_8*(kappa - kappa/(&
               0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8)**3*(0.666666666666667_real_8*GRHOAA*R0**2*mu/(RHOA**3*&
               (0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8)**2*(kappa - kappa/(0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(&
               8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8)**2) - 0.0855188292353615_real_8*&
               GRHOAA*R0*mu*SIN(7.79555417944151_real_8*R0*RHOA**(1.0_real_8/3.0_real_8)/SQRT(&
               kappa - kappa/(0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8&
               )*kappa) + 1.0_real_8) + 1.0_real_8))/(RHOA**(10.0_real_8/3.0_real_8)*(&
               0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8)**2*(kappa - kappa/(0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(&
               8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8)**(3.0_real_8/2.0_real_8)) +&
               10.128444160768_real_8*R0**2/(RHOA**(1.0_real_8/3.0_real_8)*(kappa - kappa/(&
               0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8)) - 3.89777708972075_real_8*R0*RHOA**(1.0_real_8/3.0_real_8)*(&
               0.171037658470723_real_8*GRHOAA*R0*mu/(RHOA**(10.0_real_8/3.0_real_8)*(&
               0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8)**2*(kappa - kappa/(0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(&
               8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8)**(3.0_real_8/2.0_real_8)) +&
               2.59851805981384_real_8*R0/(RHOA**(2.0_real_8/3.0_real_8)*SQRT(kappa - kappa/(&
               0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8)))*COS(7.79555417944151_real_8*R0*RHOA**(1.0_real_8/3.0_real_8)/&
               SQRT(kappa - kappa/(0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/&
               3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8))/SQRT(kappa - kappa/(&
               0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8) - 1.29925902990692_real_8*R0*SIN(7.79555417944151_real_8*R0&
               *RHOA**(1.0_real_8/3.0_real_8)/SQRT(kappa - kappa/(0.0164553078460206_real_8*&
               GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8))/(RHOA**(&
               2.0_real_8/3.0_real_8)*SQRT(kappa - kappa/(0.0164553078460206_real_8*GRHOAA*mu/(&
               RHOA**(8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8)) + 0.5_real_8*(&
               0.171037658470723_real_8*GRHOAA*R0*mu/(RHOA**(10.0_real_8/3.0_real_8)*(&
               0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8)**2*(kappa - kappa/(0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(&
               8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8)**(3.0_real_8/2.0_real_8)) +&
               2.59851805981384_real_8*R0/(RHOA**(2.0_real_8/3.0_real_8)*SQRT(kappa - kappa/(&
               0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8)))*SIN(7.79555417944151_real_8*R0*RHOA**(1.0_real_8/3.0_real_8)/&
               SQRT(kappa - kappa/(0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/&
               3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8)))/R0**4
          v1xb = 0.0_real_8
          v2xa = -0.00806288360829988_real_8*(0.25_real_8*r0**2*mu*COS(&
               7.79555417944151_real_8*R0*RHOA**(1.0_real_8/3.0_real_8)/SQRT(kappa - kappa/(&
               0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8))/(RHOA**2*(0.0164553078460206_real_8*GRHOAA*mu/(RHOA**&
               (8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8)**2*(kappa - kappa/(&
               0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8)**2) - 0.25_real_8*R0**2*mu/(RHOA**2*(&
               0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8)**2*(kappa - kappa/(0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(&
               8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8)**2))*(kappa - kappa/(&
               0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8)**3/R0**4 - 0.000398031695703623_real_8*mu*(kappa -&
               kappa/(0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa&
               ) + 1.0_real_8) + 1.0_real_8)**2*(15.192666241152_real_8*R0**2*RHOA**(2.0_real_8/&
               3.0_real_8)/(kappa - kappa/(0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(&
               8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8) - 3.89777708972075_real_8*R0*&
               RHOA**(1.0_real_8/3.0_real_8)*SIN(7.79555417944151_real_8*R0*RHOA**(1.0_real_8/3.0_real_8)&
               /SQRT(kappa - kappa/(0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8&
               /3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8))/SQRT(kappa - kappa/(&
               0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8) - 0.5_real_8*COS(7.79555417944151_real_8*R0*RHOA**(1.0_real_8/&
               3.0_real_8)/SQRT(kappa - kappa/(0.0164553078460206_real_8*GRHOAA*mu/(RHOA**&
               (8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8)) + 0.5_real_8)/(R0**4*RHOA**(&
               8.0_real_8/3.0_real_8)*(0.0164553078460206_real_8*GRHOAA*mu/(RHOA**(8.0_real_8/3.0_real_8)&
               *kappa) + 1.0_real_8)**2)
          v2xb = 0.0_real_8
          v2xab = 0.0_real_8
       ENDIF
    ELSE
       IF (rhob.GT.tol) THEN
          sx = -0.00403144180414994_real_8*(kappa - kappa/(0.0164553078460206_real_8*&
               GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8)**3*(&
               15.192666241152_real_8*R0**2*RHOB**(2.0_real_8/3.0_real_8)/(kappa - kappa/(&
               0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8) - 3.89777708972075_real_8*R0*RHOB**(1.0_real_8/3.0_real_8)*SIN(&
               7.79555417944151_real_8*R0*RHOB**(1.0_real_8/3.0_real_8)/SQRT(kappa - kappa/(&
               0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8))/SQRT(kappa - kappa/(0.0164553078460206_real_8*GRHOBB*&
               mu/(RHOB**(8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8) - 0.5_real_8*COS(&
               7.79555417944151_real_8*R0*RHOB**(1.0_real_8/3.0_real_8)/SQRT(kappa - kappa/(&
               0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8)) + 0.5_real_8)/R0**4
          v1xa = 0.0_real_8
          v1xb = 0.00053070892760483_real_8*grhobb*mu*(kappa - kappa/(&
               0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8)**2*(15.192666241152_real_8*R0**2*RHOB**(2.0_real_8/3.0_real_8)/(&
               kappa - kappa/(0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8&
               )*kappa) + 1.0_real_8) + 1.0_real_8) - 3.89777708972075_real_8*R0*RHOB**(1.0_real_8/&
               3.0_real_8)*SIN(7.79555417944151_real_8*R0*RHOB**(1.0_real_8/3.0_real_8)/SQRT(kappa -&
               kappa/(0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8)*kappa&
               ) + 1.0_real_8) + 1.0_real_8))/SQRT(kappa - kappa/(0.0164553078460206_real_8*&
               GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8) - 0.5_real_8*&
               COS(7.79555417944151_real_8*R0*RHOB**(1.0_real_8/3.0_real_8)/SQRT(kappa - kappa/&
               (0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8)) + 0.5_real_8)/(R0**4*RHOB**(11.0_real_8/3.0_real_8)*(&
               0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8)**2) - 0.00403144180414994_real_8*(kappa - kappa/(&
               0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8)**3*(0.666666666666667_real_8*GRHOBB*R0**2*mu/(RHOB**3*&
               (0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8)**2*(kappa - kappa/(0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(&
               8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8)**2) - 0.0855188292353615_real_8*&
               GRHOBB*R0*mu*SIN(7.79555417944151_real_8*R0*RHOB**(1.0_real_8/3.0_real_8)/SQRT(&
               kappa - kappa/(0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8&
               )*kappa) + 1.0_real_8) + 1.0_real_8))/(RHOB**(10.0_real_8/3.0_real_8)*(&
               0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8)**2*(kappa - kappa/(0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(&
               8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8)**(3.0_real_8/2.0_real_8)) +&
               10.128444160768_real_8*R0**2/(RHOB**(1.0_real_8/3.0_real_8)*(kappa - kappa/(&
               0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8)) - 3.89777708972075_real_8*R0*RHOB**(1.0_real_8/3.0_real_8)*(&
               0.171037658470723_real_8*GRHOBB*R0*mu/(RHOB**(10.0_real_8/3.0_real_8)*(&
               0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8)**2*(kappa - kappa/(0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(&
               8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8)**(3.0_real_8/2.0_real_8)) +&
               2.59851805981384_real_8*R0/(RHOB**(2.0_real_8/3.0_real_8)*SQRT(kappa - kappa/(&
               0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8)))*COS(7.79555417944151_real_8*R0*RHOB**(1.0_real_8/3.0_real_8)/&
               SQRT(kappa - kappa/(0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/&
               3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8))/SQRT(kappa - kappa/(&
               0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8) - 1.29925902990692_real_8*R0*SIN(7.79555417944151_real_8*R0&
               *RHOB**(1.0_real_8/3.0_real_8)/SQRT(kappa - kappa/(0.0164553078460206_real_8*&
               GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8))/(RHOB**(&
               2.0_real_8/3.0_real_8)*SQRT(kappa - kappa/(0.0164553078460206_real_8*GRHOBB*mu/(&
               RHOB**(8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8)) + 0.5_real_8*(&
               0.171037658470723_real_8*GRHOBB*R0*mu/(RHOB**(10.0_real_8/3.0_real_8)*(&
               0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8)**2*(kappa - kappa/(0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(&
               8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8)**(3.0_real_8/2.0_real_8)) +&
               2.59851805981384_real_8*R0/(RHOB**(2.0_real_8/3.0_real_8)*SQRT(kappa - kappa/(&
               0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8)))*SIN(7.79555417944151_real_8*R0*RHOB**(1.0_real_8/3.0_real_8)/&
               SQRT(kappa - kappa/(0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/&
               3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8)))/R0**4
          v2xa = 0.0_real_8
          v2xb = -0.00806288360829988_real_8*(0.25_real_8*r0**2*mu*COS(&
               7.79555417944151_real_8*R0*RHOB**(1.0_real_8/3.0_real_8)/SQRT(kappa - kappa/(&
               0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8))/(RHOB**2*(0.0164553078460206_real_8*GRHOBB*mu/(RHOB**&
               (8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8)**2*(kappa - kappa/(&
               0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8)**2) - 0.25_real_8*R0**2*mu/(RHOB**2*(&
               0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8)**2*(kappa - kappa/(0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(&
               8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8)**2))*(kappa - kappa/(&
               0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8)**3/R0**4 - 0.000398031695703623_real_8*mu*(kappa -&
               kappa/(0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8)*kappa&
               ) + 1.0_real_8) + 1.0_real_8)**2*(15.192666241152_real_8*R0**2*RHOB**(2.0_real_8/&
               3.0_real_8)/(kappa - kappa/(0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(&
               8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8) - 3.89777708972075_real_8*R0*&
               RHOB**(1.0_real_8/3.0_real_8)*SIN(7.79555417944151_real_8*R0*RHOB**(1.0_real_8/3.0_real_8)&
               /SQRT(kappa - kappa/(0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8&
               /3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8))/SQRT(kappa - kappa/(&
               0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8)*kappa) +&
               1.0_real_8) + 1.0_real_8) - 0.5_real_8*COS(7.79555417944151_real_8*R0*RHOB**(1.0_real_8/&
               3.0_real_8)/SQRT(kappa - kappa/(0.0164553078460206_real_8*GRHOBB*mu/(RHOB**&
               (8.0_real_8/3.0_real_8)*kappa) + 1.0_real_8) + 1.0_real_8)) + 0.5_real_8)/(R0**4*RHOB**(&
               8.0_real_8/3.0_real_8)*(0.0164553078460206_real_8*GRHOBB*mu/(RHOB**(8.0_real_8/3.0_real_8)&
               *kappa) + 1.0_real_8)**2)
          v2xab = 0.0_real_8
       ELSE
          sx=0.0_real_8
          v1xa=0.0_real_8
          v1xb=0.0_real_8
          v2xa=0.0_real_8
          v2xb=0.0_real_8
          v2xab=0.0_real_8
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE uks_x_zpbe
  ! ==================================================================
  SUBROUTINE revpbex(rho,grho,sx,v1x,v2x)
    ! ==--------------------------------------------------------------==
    ! Y. ZHANG ET AL. PRL 80890 (1998)
    REAL(real_8)                             :: rho, grho, sx, v1x, v2x

    REAL(real_8), PARAMETER :: ax = -0.738558766382022406_real_8, &
      uk = 1.2450_real_8, um = 0.2195149727645171_real_8, ul = um/uk , &
      us = 0.161620459673995492_real_8

    REAL(real_8)                             :: aa, dfx, ex, fx, po, rr, s2

! ==--------------------------------------------------------------==

    aa    = grho
    rr    = rho**(-4._real_8/3._real_8)
    ex    = ax/rr
    s2    = aa*rr*rr*us*us
    po    = 1._real_8/(1._real_8 + ul*s2)
    fx    = uk-uk*po
    sx    = ex*fx
    dfx   = 2._real_8*uk*ul*po*po
    v1x   = 1.33333333333333_real_8*ax*rho**0.333333333333_real_8*(fx-s2*dfx)
    v2x   = ex*dfx*(us*rr)**2
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE revpbex
  ! ==================================================================
  SUBROUTINE perdew86(rho,grho,sc,v1c,v2c)
    ! ==--------------------------------------------------------------==
    ! PERDEW CORRELATION: PRB 33, 8822 (1986)
    REAL(real_8)                             :: rho, grho, sc, v1c, v2c

    REAL(real_8), PARAMETER :: ob3 = 1._real_8/3._real_8 , &
      p1 = 0.023266_real_8, p2 = 7.389e-6_real_8, p3 = 8.723_real_8, &
      p4 = 0.472_real_8 , pc1 = 0.001667_real_8, pc2 = 0.002568_real_8, &
      pci = pc1+pc2 

    REAL(real_8)                             :: a, aa, br1, br2, br4, cn, &
                                                cna, cnb, dcn, dcna, dcnb, &
                                                drs, ephi, phi, rs, rs2, rs3

! ==--------------------------------------------------------------==

    aa    = grho
    a     = SQRT(aa)
    br1   = rho**ob3
    br2   = br1*br1
    br4   = br2*br2
    rs    = (3._real_8/(fpi*rho))**ob3
    rs2   = rs*rs
    rs3   = rs*rs2
    cna   = pc2+p1*rs+p2*rs2
    cnb   = 1._real_8+p3*rs+p4*rs2+1.e4_real_8*p2*rs3
    cn    = pc1 + cna/cnb
    drs   = -ob3*(3._real_8/fpi)**ob3 / br4
    dcna  = (p1+2._real_8*p2*rs)*drs
    dcnb  = (p3+2._real_8*p4*rs+3.e4_real_8*p2*rs2)*drs
    dcn   = dcna/cnb - cna/(cnb*cnb)*dcnb
    phi   = 0.192_real_8*pci/cn*a*rho**(-7._real_8/6._real_8)
    ephi  = EXP(-phi)
    sc    = aa/br4*cn*ephi
    v1c   = sc*((1._real_8+phi)*dcn/cn -((4._real_8/3._real_8)-(7._real_8/6._real_8)*phi)/rho)
    v2c   = cn*ephi/br4*(2._real_8-phi)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE perdew86
  ! ==================================================================
  SUBROUTINE glyp(rho,grho,sc,v1c,v2c)
    ! ==--------------------------------------------------------------==
    ! LEE, YANG PARR: GRADIENT CORRECTION PART
    REAL(real_8)                             :: rho, grho, sc, v1c, v2c

    REAL(real_8), PARAMETER                  :: a = 0.04918_real_8, &
                                                b = 0.132_real_8, &
                                                c = 0.2533_real_8, &
                                                d = 0.349_real_8 

    REAL(real_8)                             :: aa, dom, dr5, dxl, ff, om, r, &
                                                r5, xl

! ==--------------------------------------------------------------==

    aa    = grho
    r     = rho**(-1._real_8/3._real_8)
    om    = EXP(-c*r)/(1._real_8+d*r)
    r5    = r**5
    xl    = 1._real_8+(7._real_8/3._real_8)*(c*r + d*r/(1._real_8+d*r))
    ff    = a*b*aa/24._real_8
    sc    = ff*r5*om*xl
    dr5   = 5._real_8*r*r*r*r
    dom   = -om*(c+d+c*d*r)/(1._real_8+d*r)
    dxl   = (7._real_8/3._real_8)*(c+d+2._real_8*c*d*r+c*d*d*r*r)/(1._real_8+d*r)**2
    v1c   = -ff*(r*r*r*r)/3._real_8*( dr5*om*xl + r5*dom*xl + r5*om*dxl)
    v2c   = a*b*r5*om*xl/12._real_8
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE glyp
  ! ==================================================================
  SUBROUTINE ggac(rho,grho,sc,v1c,v2c)
    ! ==--------------------------------------------------------------==
    ! PERDEW & WANG GGA CORRELATION PART      
    REAL(real_8)                             :: rho, grho, sc, v1c, v2c

    REAL(real_8), PARAMETER :: al = 0.09_real_8, cx = -0.001667_real_8, &
      cxc0 = 0.002568_real_8, cc0 = -cx+cxc0 , ob3 = 1._real_8/3._real_8 , &
      pa = 0.023266_real_8, pb = 7.389e-6_real_8, pc = 8.723_real_8, &
      pd = 0.472_real_8

    REAL(real_8) :: a, aa, af, be, bf, cn, cna, cnb, dcn, dcna, dcnb, ddh0, &
      ddh1, dh0, dh1, ec, ee, ex, expe, h0, h1, qy, rs, rs2, rs3, s1, t, vc, &
      vx, xkf, xks, xnu, xy, y

! ==--------------------------------------------------------------==

    xnu   = 16._real_8/pi*(3._real_8*pi*pi)**ob3
    be    = xnu*cc0
    CALL xc(rho,ex,ec,vx,vc)
    aa    = grho
    a     = SQRT(aa)
    rs    = (3._real_8/(4._real_8*pi*rho))**ob3
    rs2   = rs*rs
    rs3   = rs*rs2
    xkf   = (9._real_8*pi/4._real_8)**ob3/rs
    xks   = SQRT(4._real_8*xkf/pi)
    t     = a/(2._real_8*xks*rho)
    expe  = EXP(-2._real_8*al*ec/(be*be))
    af    = 2._real_8*al/be * (1._real_8/(expe-1._real_8))
    bf    = expe*(vc-ec)
    y     = af*t*t
    xy    = (1._real_8+y)/(1._real_8+y+y*y)
    qy    = y*y*(2._real_8+y)/(1._real_8+y+y*y)**2
    s1    = 1._real_8+2._real_8*al/be*t*t*xy
    h0    = be*be/(2._real_8*al) * LOG(s1)
    dh0   = be*t*t/s1*(-7._real_8/3._real_8*xy-qy*(af*bf/be-7._real_8/3._real_8))
    ddh0  = be/(2._real_8*xks*xks*rho)*(xy-qy)/s1
    ee    = -100._real_8*(xks/xkf*t)**2
    cna   = cxc0+pa*rs+pb*rs2
    dcna  = -(pa*rs+2._real_8*pb*rs2)/3._real_8
    cnb   = 1._real_8+pc*rs+pd*rs2+1.e4_real_8*pb*rs3
    dcnb  = -(pc*rs+2._real_8*pd*rs2+3.e4_real_8*pb*rs3)/3._real_8
    cn    = cna/cnb - cx
    dcn   = dcna/cnb - cna*dcnb/(cnb*cnb)
    h1    = xnu*(cn-cc0-3._real_8/7._real_8*cx)*t*t*EXP(ee)
    dh1   = -ob3*(h1*(7._real_8+8._real_8*ee)+xnu*t*t*EXP(ee)*dcn)
    ddh1  = 2._real_8*h1*(1._real_8+ee)*rho/aa
    sc    = rho*(h0+h1)
    v1c   = h0+h1+dh0+dh1
    v2c   = ddh0+ddh1
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE ggac
  ! ==================================================================
  SUBROUTINE pbesc(rho,grho,w1,sc,v1c,v2c)
    ! ==--------------------------------------------------------------==
    ! PBESol Correlation functional
    REAL(real_8)                             :: rho, grho, w1, sc, v1c, v2c

    REAL(real_8), PARAMETER :: be = 0.046_real_8, &
      ga = 0.031090690869654895_real_8 , ob3 = 1._real_8/3._real_8 

    REAL(real_8)                             :: a, aa, af, dadr, dhdr, dhdt, &
                                                dsda, dsdr, dsdt, dtdr, ec, &
                                                ex, expe, h0, rs, s1, t, vc, &
                                                vx, xkf, xks, xy, y

! ==--------------------------------------------------------------==

    CALL xc(rho,ex,ec,vx,vc)
    aa    = grho
    a     = SQRT(aa)
    rs    = (3._real_8/(4._real_8*pi*rho))**ob3
    xkf   = (9._real_8*pi/4._real_8)**ob3/rs
    xks   = SQRT(4._real_8*xkf/pi)
    t     = a/(2._real_8*xks*rho)
    expe  = EXP(-ec/ga)
    af    = be/ga * (1._real_8/(expe-1._real_8))
    y     = af*t*t
    xy    = (1._real_8+y)/(1._real_8+y+y*y)
    s1    = 1._real_8+be/ga*t*t*xy
    h0    = ga * LOG(s1)
    dtdr  = -t*7._real_8/(6._real_8*rho)
    dadr  = af*af*expe/be*(vc-ec)/rho
    dsda  = -be/ga * af * t**6 * (2._real_8+y) / (1._real_8+y+y*y)**2
    dsdt  = 2._real_8*be/ga * t * (1._real_8+2._real_8*y) / (1._real_8+y+y*y)**2
    dsdr  = dsda*dadr + dsdt*dtdr
    dhdt  = ga/s1*dsdt
    dhdr  = ga/s1*dsdr
    sc    = w1*rho*h0
    v1c   = w1*h0+w1*dhdr*rho
    v2c   = w1*rho*dhdt*t/aa
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE pbesc
  ! ==================================================================
  SUBROUTINE pbec(rho,grho,w1,sc,v1c,v2c)
    ! ==--------------------------------------------------------------==
    ! PBE Correlation functional
    REAL(real_8)                             :: rho, grho, w1, sc, v1c, v2c

    REAL(real_8), PARAMETER :: be = 0.06672455060314922_real_8, &
      ga = 0.031090690869654895_real_8 , ob3 = 1._real_8/3._real_8 

    REAL(real_8)                             :: a, aa, af, dadr, dhdr, dhdt, &
                                                dsda, dsdr, dsdt, dtdr, ec, &
                                                ex, expe, h0, rs, s1, t, vc, &
                                                vx, xkf, xks, xy, y

! ==--------------------------------------------------------------==

    CALL xc(rho,ex,ec,vx,vc)
    aa    = grho
    a     = SQRT(aa)
    rs    = (3._real_8/(4._real_8*pi*rho))**ob3
    xkf   = (9._real_8*pi/4._real_8)**ob3/rs
    xks   = SQRT(4._real_8*xkf/pi)
    t     = a/(2._real_8*xks*rho)
    expe  = EXP(-ec/ga)
    af    = be/ga * (1._real_8/(expe-1._real_8))
    y     = af*t*t
    xy    = (1._real_8+y)/(1._real_8+y+y*y)
    s1    = 1._real_8+be/ga*t*t*xy
    h0    = ga * LOG(s1)
    dtdr  = -t*7._real_8/(6._real_8*rho)
    dadr  = af*af*expe/be*(vc-ec)/rho
    dsda  = -be/ga * af * t**6 * (2._real_8+y) / (1._real_8+y+y*y)**2
    dsdt  = 2._real_8*be/ga * t * (1._real_8+2._real_8*y) / (1._real_8+y+y*y)**2
    dsdr  = dsda*dadr + dsdt*dtdr
    dhdt  = ga/s1*dsdt
    dhdr  = ga/s1*dsdr
    sc    = w1*rho*h0
    v1c   = w1*h0+w1*dhdr*rho
    v2c   = w1*rho*dhdt*t/aa
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE pbec
  ! ==================================================================
  SUBROUTINE hcth(rho,grho,sx,v1x,v2x)
    ! HCTH, JCP 109, 6264 (1998)
    ! Parameters set-up after N.L. Doltsisnis & M. Sprik (1999)
    ! Present release: Tsukuba, 09/02/2005
    ! --------------------------------------------------------------------------
    ! rhoa = rhob = 0.5 * rho
    ! grho is the SQUARE of the gradient of rho! --> gr=sqrt(grho)
    ! sx  : total exchange correlation energy at point r 
    ! v1x : d(sx)/drho  (eq. dfdra = dfdrb in original)
    ! v2x : 1/gr*d(sx)/d(gr) (eq. 0.5 * dfdza = 0.5 * dfdzb in original)
    ! --------------------------------------------------------------------------
    REAL(real_8)                             :: rho, grho, sx, v1x, v2x

    REAL(real_8), PARAMETER                  :: fr83 = 8._real_8/3._real_8 , &
                                                o3 = 1.0_real_8/3.0_real_8

    REAL(real_8) :: bygr, denaa, denab, denx, dera1_dra, derab0_drab, &
      dex_drho, dffaa_drho, dffab_drho, dg, dgaa_dgr, dgaa_drho, dgab_dgr, &
      dgab_drho, dgx_dgr, dgx_drho, dra_drho, drab_drho, era1, erab0, ex, &
      f83rho, ffaa, ffab, g, gaa, gab, gr, gx, r3pi, r3q2, ra, rab, rho_o3, &
      rho_o34, taa, tab, txx, uaa, uab, ux, xa, xa2
    REAL(real_8), DIMENSION(6)               :: caa, cab, cg0, cg1, cx

    r3q2=EXP(-o3*0.69314718055994531_real_8)
    r3pi=EXP(-o3*0.04611759718129048_real_8)
    ! .....coefficients for PW correlation......................................
    cg0(1)= 0.031091_real_8
    cg0(2)= 0.213700_real_8
    cg0(3)= 7.595700_real_8
    cg0(4)= 3.587600_real_8
    cg0(5)= 1.638200_real_8
    cg0(6)= 0.492940_real_8
    cg1(1)= 0.015545_real_8
    cg1(2)= 0.205480_real_8
    cg1(3)=14.118900_real_8
    cg1(4)= 6.197700_real_8
    cg1(5)= 3.366200_real_8
    cg1(6)= 0.625170_real_8
    ! ......HCTH-19-4.....................................
    caa(1)=  0.489508e+00_real_8
    caa(2)= -0.260699e+00_real_8
    caa(3)=  0.432917e+00_real_8
    caa(4)= -0.199247e+01_real_8
    caa(5)=  0.248531e+01_real_8
    caa(6)=  0.200000e+00_real_8
    cab(1)=  0.514730e+00_real_8
    cab(2)=  0.692982e+01_real_8
    cab(3)= -0.247073e+02_real_8
    cab(4)=  0.231098e+02_real_8
    cab(5)= -0.113234e+02_real_8
    cab(6)=  0.006000e+00_real_8
    cx(1) =  0.109163e+01_real_8
    cx(2) = -0.747215e+00_real_8
    cx(3) =  0.507833e+01_real_8
    cx(4) = -0.410746e+01_real_8
    cx(5) =  0.117173e+01_real_8
    cx(6) =  0.004000e+00_real_8
    ! ...........................................................................
    gr=SQRT(grho)
    rho_o3=rho**(o3)
    rho_o34=rho*rho_o3
    xa=1.25992105_real_8*gr/rho_o34
    xa2=xa*xa
    ra=0.781592642_real_8/rho_o3
    rab=r3q2*ra
    dra_drho=-0.260530881_real_8/rho_o34
    drab_drho=r3q2*dra_drho
    CALL pwcorr(ra,cg1,g,dg)
    era1=g
    dera1_dra=dg
    CALL pwcorr(rab,cg0,g,dg)
    erab0=g
    derab0_drab=dg
    ex=-0.75_real_8*r3pi*rho_o34
    dex_drho=-r3pi*rho_o3
    uaa=caa(6)*xa2
    uaa=uaa/(1.0_real_8+uaa)
    uab=cab(6)*xa2
    uab=uab/(1.0_real_8+uab)
    ux=cx(6)*xa2
    ux=ux/(1.0_real_8+ux)
    ffaa=rho*era1
    ffab=rho*erab0-ffaa
    dffaa_drho=era1+rho*dera1_dra*dra_drho
    dffab_drho=erab0+rho*derab0_drab*drab_drho-dffaa_drho
    ! mb-> i-loop removed
    denaa=1._real_8/(1.0_real_8+caa(6)*xa2)
    denab=1._real_8/(1.0_real_8+cab(6)*xa2)
    denx =1._real_8/(1.0_real_8+cx(6)*xa2)
    f83rho=fr83/rho
    bygr=2.0_real_8/gr
    gaa=caa(1)+uaa*(caa(2)+uaa*(caa(3)+uaa*(caa(4)+uaa*caa(5))))
    gab=cab(1)+uab*(cab(2)+uab*(cab(3)+uab*(cab(4)+uab*cab(5))))
    gx=cx(1)+ux*(cx(2)+ux*(cx(3)+ux*(cx(4)+ux*cx(5))))
    taa=denaa*uaa*(caa(2)+uaa*(2._real_8*caa(3)+uaa&
         *(3._real_8*caa(4)+uaa*4._real_8*caa(5))))
    tab=denab*uab*(cab(2)+uab*(2._real_8*cab(3)+uab&
         *(3._real_8*cab(4)+uab*4._real_8*cab(5))))
    txx=denx*ux*(cx(2)+ux*(2._real_8*cx(3)+ux&
         *(3._real_8*cx(4)+ux*4._real_8*cx(5))))
    dgaa_drho=-f83rho*taa
    dgab_drho=-f83rho*tab
    dgx_drho=-f83rho*txx
    dgaa_dgr=bygr*taa
    dgab_dgr=bygr*tab
    dgx_dgr=bygr*txx
    ! mb
    sx=ex*gx+ffaa*gaa+ffab*gab
    v1x=dex_drho*gx+ex*dgx_drho&
         +dffaa_drho*gaa+ffaa*dgaa_drho&
         +dffab_drho*gab+ffab*dgab_drho
    v2x=(ex*dgx_dgr+ffaa*dgaa_dgr+ffab*dgab_dgr)/gr
    RETURN
  END SUBROUTINE hcth
  ! =-------------------------------------------------------------------=
  SUBROUTINE pwcorr(r,c,g,dg)
    REAL(real_8)                             :: r
    REAL(real_8), DIMENSION(6)               :: c
    REAL(real_8)                             :: g, dg

    REAL(real_8)                             :: drb, r12, r2, r32, rb, sb

    r12=SQRT(r)
    r32=r*r12
    r2=r*r
    rb=c(3)*r12+c(4)*r+c(5)*r32+c(6)*r2
    sb=1.0_real_8+1.0_real_8/(2.0_real_8*c(1)*rb)
    g=-2.0_real_8*c(1)*(1.0_real_8+c(2)*r)*LOG(sb)
    drb=c(3)/(2.0_real_8*r12)+c(4)+1.5_real_8*c(5)*r12+2.0_real_8*c(6)*r
    dg=(1.0_real_8+c(2)*r)*drb/(rb*rb*sb)-2.0_real_8*c(1)*c(2)*LOG(sb)
    RETURN
  END SUBROUTINE pwcorr
  ! ==================================================================
  SUBROUTINE optx(rho,grho,sx,v1x,v2x)
    ! OPTX, Handy et al. JCP 116, p. 5411 (2002) and refs. therein
    ! Present release: Tsukuba, 20/6/2002
    ! --------------------------------------------------------------------------
    ! rhoa = rhob = 0.5 * rho in LDA implementation
    ! grho is the SQUARE of the gradient of rho! --> gr=sqrt(grho)
    ! sx  : total exchange correlation energy at point r
    ! v1x : d(sx)/drho
    ! v2x : 1/gr*d(sx)/d(gr)
    ! --------------------------------------------------------------------------
    REAL(real_8)                             :: rho, grho, sx, v1x, v2x

    REAL(real_8), PARAMETER :: a1cx = 0.9784571170284421_real_8, &
      a2 = 1.43169_real_8 , gam = 0.006_real_8 , o43 = 4.0_real_8/3.0_real_8, &
      smal2 = 1.e-08_real_8 , small = 1.e-20_real_8, &
      two13 = 1.259921049894873_real_8 , two53 = 3.174802103936399_real_8

    REAL(real_8)                             :: gamx2, gr, rho43, uden, uu, xa

! .......coefficients and exponents....................
! .......OPTX in compact form..........................

    IF (rho.LE.small) THEN
       sx=0.0_real_8
       v1x=0.0_real_8
       v2x=0.0_real_8
    ELSE
       gr=MAX(grho,smal2)
       rho43=rho**o43
       xa=two13*SQRT(gr)/rho43
       gamx2=gam*xa*xa
       uden=1.e+00_real_8/(1.e+00_real_8+gamx2)
       uu=a2*gamx2*gamx2*uden*uden
       uden=rho43*uu*uden
       sx=-rho43*(a1cx+uu)/two13
       v1x=o43*(sx+two53*uden)/rho
       v2x=-two53*uden/gr
    ENDIF
    ! 
    RETURN
  END SUBROUTINE optx
  ! =-------------------------------------------------------------------=



  ! These are wrappers to the reference implementation
  ! below by Heyd & Scuseria.
  ! For more information about implementation, see also
  ! Komsa, Broqvist, Pasquarello, Phys. Rev. B 81, 205118 (2010)
  ! 
  ! Maybe this should be moved to lsd_func.F

  ! ==================================================================
  SUBROUTINE pbexsr_lsd(rhoa,rhob,grhoaa,grhobb,sx,&
       v1xa,v2xa,v1xb,v2xb,v2xab,omega)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rhoa, rhob, grhoaa, grhobb, &
                                                sx, v1xa, v2xa, v1xb, v2xb, &
                                                v2xab, omega

    REAL(real_8), PARAMETER                  :: small = 1.e-20_real_8 

    REAL(real_8)                             :: GRHOAA_tmp, GRHOBB_tmp, &
                                                RHOA_tmp, RHOB_tmp, sxa, sxb

! ==--------------------------------------------------------------==

    sxa=0.0_real_8
    sxb=0.0_real_8
    v1xa=0.0_real_8
    v2xa=0.0_real_8
    v1xb=0.0_real_8
    v2xb=0.0_real_8
    v2xab=0.0_real_8

    GRHOAA_tmp=4._real_8*grhoaa
    GRHOBB_tmp=4._real_8*grhobb
    RHOA_tmp=2._real_8*rhoa
    RHOB_tmp=2._real_8*rhob

    IF (ABS(rhoa).GT.small) THEN
       CALL pbexsr(RHOA_tmp, GRHOAA_tmp, sxa, v1xa, v2xa, omega)
    ENDIF
    IF (ABS(rhob).GT.small) THEN
       CALL pbexsr(RHOB_tmp, GRHOBB_tmp, sxb, v1xb, v2xb, omega)
    ENDIF
    sx = 0.5_real_8*(sxa+sxb)
    v2xa = 2._real_8*v2xa
    v2xb = 2._real_8*v2xb          ! I REALLY HOPE THIS WORKS JUST LIKE THIS

    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE pbexsr_lsd

  ! ==================================================================
  SUBROUTINE pbexsr(rho,grho,sx,v1x,v2x,omega)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rho, grho, sx, v1x, v2x, omega

    REAL(real_8), PARAMETER :: ax = -0.738558766382022406_real_8, &
      smal2 = 1.e-08_real_8 , small = 1.e-20_real_8, uk = 0.8040_real_8, &
      um = 0.2195149727645171_real_8, ul = um/uk , &
      us = 0.161620459673995492_real_8

    REAL(real_8)                             :: aa, d1x, d2x, dsdg, dsdn, ec, &
                                                ex, fx, rr, s, s2, vc, vx

! C     ==--------------------------------------------------------------==

    CALL xc(rho,ex,ec,vx,vc)

    ! AA    = MAX(GRHO,SMAL2)
    aa    = grho

    rr    = rho**(-4._real_8/3._real_8)
    ex    = ax/rr
    s2    = aa*rr*rr*us*us

    s = SQRT(s2)
    ! Lieb-Oxford bound
    IF (s.GT.8.3_real_8) THEN
       s = 8.572844_real_8 - 18.796223_real_8/s2
    ENDIF
    CALL wpbe_analytical_erfc_approx_grad(rho,s,omega,fx,d1x,d2x)
    sx = ex*fx        ! - EX
    dsdn = -4._real_8/3._real_8*s/rho
    v1x = vx*fx + (dsdn*d2x+d1x)*ex   ! - VX
    dsdg = us*rr
    v2x = ex*1._real_8/SQRT(aa)*dsdg*d2x

    ! NOTE, in PBEX() SX is only the gradient correction energy density(?),
    ! but here it is the total energy density
    ! And same for the potential V1X

    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE pbexsr


  ! --------------------------------------------------------------------
  ! 
  ! wPBE Enhancement Factor (erfc approx.,analytical, gradients)
  ! 
  ! --------------------------------------------------------------------

  SUBROUTINE wpbe_analytical_erfc_approx_grad(rho,s,omega,Fx_wpbe,&
       d1rfx,d1sfx)


    REAL(real_8)                             :: rho, s, omega, Fx_wpbe, &
                                                d1rfx, d1sfx

    REAL(real_8), PARAMETER :: Eight = 8._real_8, Fifteen = 1.5e1_real_8, &
      Five = 5._real_8, Four = 4._real_8, Nine = 9._real_8, One = 1._real_8, &
      r105 = 1.05e2_real_8, r12 = 1.2e1_real_8, r120 = 1.2e2_real_8, &
      r1215 = 1.215e3_real_8, r128 = 1.28e2_real_8, r135 = 1.35e2_real_8, &
      r144 = 1.44e2_real_8, r15309 = 1.5309e4_real_8, r1944 = 1.944e3_real_8, &
      r20 = 2.0e1_real_8, r2187 = 2.187e3_real_8, r243 = 2.43e2_real_8, &
      r25 = 2.5e1_real_8, r256 = 2.56e2_real_8 , r27 = 2.7e1_real_8, &
      r288 = 2.88e2_real_8, r30 = 3.0e1_real_8, r32 = 3.2e1_real_8, &
      r324 = 3.24e2_real_8, r36 = 3.6e1_real_8, r384 = 3.84e2_real_8
    REAL(real_8), PARAMETER :: r40 = 4.0e1_real_8 , r4374 = 4.374e3_real_8 , &
      r48 = 4.8e1_real_8, r512 = 5.12e2_real_8, r54 = 5.4e1_real_8, &
      r64 = 6.4e1_real_8, r6561 = 6.561e3_real_8, r729 = 7.29e2_real_8, &
      r75 = 7.5e1_real_8 , r81 = 8.1e1_real_8, r864 = 8.64e2_real_8, &
      Seven = 7._real_8, Six = 6._real_8, Sixteen = 1.6e1_real_8, &
      Ten = 10._real_8, Three = 3._real_8, Two = 2._real_8, Zero = 0._real_8

    REAL(real_8) :: a, a12, a2, a3, a32, a4, a52, a72, b, c, d, d1rDHsw, &
      d1rexpei, d1rf2, d1rf3, d1rf4, d1rf5, d1rf6, d1rf7, d1rf8, d1rf9, &
      d1rHsbw, d1rnp1, d1rnp2, d1rpiexperf, d1rt1, d1rt10, d1rt2t9, d1rterm1, &
      d1rterm3, d1rterm4, d1rterm5, d1rw, d1sDHs, d1sEG, d1sexpei, d1sF, &
      d1sf2, d1sf3, d1sf4, d1sf5, d1sf6, d1sf7, d1sf8, d1sf9, d1sG_a, d1sG_b, &
      d1sH, d1sHden, d1sHnum, d1sHsbw, d1spiexperf, d1st1, d1st10, d1st2t9, &
      d1sterm1, d1sterm2, d1sterm3, d1sterm4, d1sterm5, DHs, DHs2, DHs3, &
      DHs4, DHs72, DHs92, DHsbw, DHsbw12, DHsbw2, DHsbw3, DHsbw32, DHsbw4, &
      DHsbw5, DHsbw52, DHsbw72, DHsbw92, DHsw, DHsw2
    REAL(real_8) :: DHsw52, DHsw72, e, ea1, ea2, ea3, ea4, ea5, ea6, ea7, &
      ea8, eb1, eg, EGa1, EGa2, EGa3, EGscut, expei, expei1, expei2, expei3, &
      expei4, expeid1, expfcutoff, f, f12, f13, f14, f1516, f18, f2, f23, &
      f2d1, f3, f32, f34, f3d1, f4, f43, f4d1, f5, f5d1, f6, f6d1, f7, f72, &
      f8, f8d1, f9, f94, f98, f9d1, Fc1, Fc2, G_a, G_b, h, Ha1, Ha2, Ha3, &
      Ha4, Ha5, Hden, Hnum, Hsbw, Hsbw12, Hsbw2, Hsbw3, Hsbw32, Hsbw4, &
      Hsbw52, Hsbw72, HsbwA94, HsbwA9412, HsbwA942, HsbwA943, HsbwA945, np1, &
      np2, pi, pi2, pi_23, piexperf, piexperfd1, s2, s3, s4, s5, s6, srpi, &
      t1, t10, t10d1, t2t9, term1, term1d1, term2
    REAL(real_8) :: term3, term4, term5, Three_13, w, w2, w3, w4, w5, w6, w7, &
      w8, wcutoff, x, xkf, xkfrho

    f12    = 0.5_real_8
    f13    = One/Three
    f14    = 0.25_real_8
    f18    = 0.125_real_8

    f23    = two * f13
    f43    = two * f23

    f32    = 1.5_real_8
    f72    = 3.5_real_8
    f34    = 0.75_real_8
    f94    = 2.25_real_8
    f98    = 1.125_real_8
    f1516  = Fifteen / Sixteen

    pi     = ACOS(-One)
    pi2    = pi*pi
    pi_23  = pi2**f13
    srpi   = SQRT(pi)

    Three_13 = Three**f13

    ! Constants from fit

    ea1 = -1.128223946706117_real_8
    ea2 = 1.452736265762971_real_8
    ea3 = -1.243162299390327_real_8
    ea4 = 0.971824836115601_real_8
    ea5 = -0.568861079687373_real_8
    ea6 = 0.246880514820192_real_8
    ea7 = -0.065032363850763_real_8
    ea8 = 0.008401793031216_real_8

    eb1 = 1.455915450052607_real_8

    ! Constants for PBE hole

    a      =  1.0161144_real_8
    b      = -3.7170836e-1_real_8
    c      = -7.7215461e-2_real_8
    d      =  5.7786348e-1_real_8
    e      = -5.1955731e-2_real_8
    x      = - Eight/Nine

    ! Constants for fit of H(s) (PBE)

    Ha1    = 9.79681e-3_real_8
    Ha2    = 4.10834e-2_real_8
    Ha3    = 1.87440e-1_real_8
    Ha4    = 1.20824e-3_real_8
    Ha5    = 3.47188e-2_real_8

    ! Constants for F(H) (PBE)

    Fc1    = 6.4753871_real_8
    Fc2    = 4.7965830e-1_real_8

    ! Constants for polynomial expansion for EG for small s

    EGa1   = -2.628417880e-2_real_8
    EGa2   = -7.117647788e-2_real_8
    EGa3   =  8.534541323e-2_real_8

    ! Constants for large x expansion of exp(x)*ei(-x)

    expei1 = 4.03640_real_8
    expei2 = 1.15198_real_8
    expei3 = 5.03627_real_8
    expei4 = 4.19160_real_8

    ! Cutoff criterion below which to use polynomial expansion

    EGscut     = 8.0e-2_real_8
    wcutoff    = 1.4e1_real_8
    expfcutoff = 7.0e2_real_8

    ! Calculate prelim variables

    xkf    = (Three*pi2*rho) ** f13
    xkfrho = xkf * rho

    a2 = a*a
    a3 = a2*a
    a4 = a3*a
    a12 = SQRT(a)
    a32 = a12*a
    a52 = a32*a
    a72 = a52*a

    w      = omega / xkf
    w2    = w * w
    w3    = w2 * w
    w4    = w2 * w2
    w5    = w3 * w2
    w6    = w5 * w
    w7    = w6 * w
    w8    = w7 * w

    d1rw  = -(One/(Three*rho))*w

    x      = - Eight/Nine

    s2     = s*s
    s3     = s2*s
    s4     = s2*s2
    s5     = s4*s
    s6     = s5*s

    ! Calculate wPBE enhancement factor

    Hnum    = Ha1*s2 + Ha2*s4
    Hden    = One + Ha3*s4 + Ha4*s5 + Ha5*s6

    h       = Hnum/Hden

    d1sHnum = two*Ha1*s + Four*Ha2*s3
    d1sHden = Four*Ha3*s3 + Five*Ha4*s4 + Six*Ha5*s5

    d1sH    = (Hden*d1sHnum - Hnum*d1sHden) / (hden*hden)

    f      = Fc1*h + Fc2
    d1sF   = Fc1*d1sH

    ! Change exponent of Gaussian if we're using the simple approx.

    IF (w .GT. wcutoff) THEN

       eb1 = 2.0_real_8

    ENDIF

    ! Calculate helper variables (should be moved later on...)

    Hsbw = s2*h + eb1*w2
    Hsbw2 = Hsbw*hsbw
    Hsbw3 = Hsbw2*Hsbw
    Hsbw4 = Hsbw3*Hsbw
    Hsbw12 = SQRT(Hsbw)
    Hsbw32 = Hsbw12*Hsbw
    Hsbw52 = Hsbw32*Hsbw
    Hsbw72 = Hsbw52*Hsbw

    d1sHsbw  = d1sH*s2 + two*s*h
    d1rHsbw  = two*eb1*d1rw*w

    DHsbw = d + s2*h + eb1*w2
    DHsbw2 = DHsbw*dhsbw
    DHsbw3 = DHsbw2*DHsbw
    DHsbw4 = DHsbw3*DHsbw
    DHsbw5 = DHsbw4*DHsbw
    DHsbw12 = SQRT(DHsbw)
    DHsbw32 = DHsbw12*DHsbw
    DHsbw52 = DHsbw32*DHsbw
    DHsbw72 = DHsbw52*DHsbw
    DHsbw92 = DHsbw72*DHsbw

    HsbwA94   = f94 * Hsbw / a
    HsbwA942  = HsbwA94*hsbwa94
    HsbwA943  = HsbwA942*HsbwA94
    HsbwA945  = HsbwA943*HsbwA942
    HsbwA9412 = SQRT(HsbwA94)

    DHs    = d + s2*h
    DHs2   = DHs*dhs
    DHs3   = DHs2*DHs
    DHs4   = DHs3*DHs
    DHs72  = DHs3*SQRT(DHs)
    DHs92  = DHs72*DHs

    d1sDHs = two*s*h + s2*d1sH

    DHsw   = DHs + w2
    DHsw2  = DHsw*dhsw
    DHsw52 = SQRT(DHsw)*DHsw2
    DHsw72 = DHsw52*DHsw

    d1rDHsw = two*d1rw*w

    IF (s .GT. EGscut) THEN

       G_a    = srpi * (Fifteen*e + Six*c*(One+f*s2)*DHs +&
            Four*b*(DHs2) + Eight*a*(DHs3))&
            * (One / (Sixteen * DHs72))&
            - f34*pi*SQRT(a) * EXP(f94*h*s2/a) *&
            (One - cp_erf(f32*s*SQRT(h/a)))

       d1sG_a = (One/r32)*srpi *&
            ((r36*(two*h + d1sH*s) / (a12*SQRT(h/a)))&
            + (One/DHs92) *&
            (-Eight*a*d1sDHs*DHs3 - r105*d1sdhs*e&
            -r30*c*d1sDHs*DHs*(One+s2*f)&
            +r12*DHs2*(-b*d1sDHs + c*s*(d1sF*s + two*f)))&
            - ((r54*EXP(f94*h*s2/a)*srpi*s*(two*h+d1sH*s)*&
            cp_erfc(f32*SQRT(h/a)*s))&
            / a12))

       G_b    = (f1516 * srpi * s2) / DHs72

       d1sG_b = (Fifteen*srpi*s*(Four*DHs - Seven*d1sDHs*s))&
            / (r32*DHs92)

       eg     = - (f34*pi + G_a) / G_b

       d1sEG  = (-Four*d1sG_a*G_b + d1sG_b*(four*G_a + Three*pi))&
            / (Four*G_b*g_b)

    ELSE

       eg    = EGa1 + EGa2*s2 + EGa3*s4
       d1sEG = two*EGa2*s + Four*EGa3*s3

    ENDIF

    ! Calculate the terms needed in any case

    term2 =       (DHs2*b + DHs*c + two*e + dhs*s2*c*f + two*s2*eg) /&
         (two*DHs3)

    d1sterm2 = (-Six*d1sDHs*(eg*s2 + e)&
         + DHs2 * (-d1sDHs*b + s*c*(d1sF*s + two*f))&
         + two*DHs * (two*eg*s - d1sDHs*c&
         + s2 * (d1sEG - d1sDHs*c*f)))&
         / (two*DHs4)

    term3 = - w  * (Four*DHsw2*b + Six*DHsw*c + Fifteen*e&
         + Six*DHsw*s2*c*f + Fifteen*s2*eg) /&
         (Eight*DHs*DHsw52)

    d1sterm3 = w * (two*d1sDHs*DHsw * (Four*DHsw2*b&
         + Six*DHsw*c + Fifteen*e&
         + Three*s2*(Five*eg + two*DHsw*c*f))&
         + DHs * (r75*d1sDHs*(eg*s2 + e)&
         + Four*DHsw2*(d1sDHs*b&
         - Three*s*c*(d1sF*s + two*f))&
         - Six*DHsw*(-Three*d1sDHs*c&
         + s*(Ten*eg + Five*d1sEG*s&
         - Three*d1sDHs*s*c*f))))&
         / (Sixteen*DHs2*DHsw72)

    d1rterm3 = (-two*d1rw*DHsw * (Four*DHsw2*b&
         + Six*DHsw*c + Fifteen*e&
         + Three*s2*(Five*eg + two*DHsw*c*f))&
         + w * d1rDHsw * (r75*(eg*s2 + e)&
         + two*DHsw*(two*dhsw*b + Nine*c&
         + Nine*s2*c*f)))&
         / (Sixteen*DHs*DHsw72)

    term4 = - w3 * (DHsw*c + Five*e + dhsw*s2*c*f + five*s2*eg) /&
         (two*DHs2*DHsw52)

    d1sterm4 = (w3 * (Four*d1sDHs*DHsw * (dhsw*c + Five*e&
         + s2 * (Five*eg + DHsw*c*f))&
         + DHs * (r25*d1sDHs*(eg*s2 + e)&
         - two*DHsw2*s*c*(d1sF*s + two*f)&
         + DHsw * (Three*d1sDHs*c + s*(-r20*eg&
         - Ten*d1sEG*s&
         + Three*d1sDHs*s*c*f)))))&
         / (Four*DHs3*DHsw72)

    d1rterm4 = (w2 * (-Six*d1rw*DHsw * (dhsw*c + Five*e&
         + s2 * (Five*eg + DHsw*c*f))&
         + w * d1rDHsw * (r25*(eg*s2 + e) +&
         Three*DHsw*c*(One + s2*f))))&
         / (Four*DHs2*DHsw72)

    term5 = - w5 * (e + s2*eg) /&
         (DHs3*DHsw52)

    d1sterm5 = (w5 * (Six*d1sDHs*DHsw*(eg*s2 + e)&
         + DHs * (-two*DHsw*s * (two*eg + d1sEG*s)&
         + Five*d1sDHs * (eg*s2 + e))))&
         / (two*DHs4*DHsw72)

    d1rterm5 = (w4 * Five*(eg*s2 + e) * (-two*d1rw*DHsw&
         + d1rDHsw * w))&
         / (two*DHs3*DHsw72)


    IF ((s.GT.0.0_real_8).OR.(w.GT.0.0_real_8)) THEN

       t10    = (f12)*a*LOG(Hsbw / DHsbw)
       t10d1  = f12*a*(One/Hsbw - one/DHsbw)
       d1st10 = d1sHsbw*t10d1
       d1rt10 = d1rHsbw*t10d1

    ENDIF

    ! Calculate exp(x)*f(x) depending on size of x

    IF (HsbwA94 .LT. expfcutoff) THEN

       piexperf = pi*EXP(HsbwA94)*cp_erfc(HsbwA9412)
       ! expei    = Exp(HsbwA94)*Ei(-HsbwA94)
       expei    = EXP(HsbwA94)*(-expint(1,hsbwa94))

    ELSE

       piexperf = pi*(One/(srpi*HsbwA9412)&
            - One/(two*SQRT(pi*HsbwA943))&
            + Three/(Four*SQRT(pi*HsbwA945)))

       expei  = - (One/HsbwA94) *&
            (HsbwA942 + expei1*HsbwA94 + expei2) /&
            (HsbwA942 + expei3*HsbwA94 + expei4)

    ENDIF

    ! Calculate the derivatives (based on the orig. expression)
    ! --> Is this ok? ==> seems to be ok...

    piexperfd1  = - (Three*srpi*SQRT(Hsbw/a))/(two*hsbw)&
         + (Nine*piexperf)/(Four*a)
    d1spiexperf = d1sHsbw*piexperfd1
    d1rpiexperf = d1rHsbw*piexperfd1

    expeid1  = f14*(Four/Hsbw + (Nine*expei)/a)
    d1sexpei = d1sHsbw*expeid1
    d1rexpei = d1rHsbw*expeid1

    IF (w .EQ. Zero) THEN

       ! Fall back to original expression for the PBE hole

       t1 = -f12*a*expei
       d1st1 = -f12*a*d1sexpei
       d1rt1 = -f12*a*d1rexpei

       IF (s .GT. 0.0_real_8) THEN

          term1    = t1 + t10
          d1sterm1 = d1st1 + d1st10
          d1rterm1 = d1rt1 + d1rt10

          Fx_wpbe = x * (term1 + term2)

          d1sfx = x * (d1sterm1 + d1sterm2)
          d1rfx = x * d1rterm1

       ELSE

          Fx_wpbe = 1.0_real_8

          ! TODO    This is checked to be true for term1
          ! How about the other terms???

          d1sfx   = 0.0_real_8
          d1rfx   = 0.0_real_8

       ENDIF


    ELSEIF (w .GT. wcutoff) THEN

       ! Use simple Gaussian approximation for large w

       ! WRITE(6,*)rho,s," LARGE w"

       term1 = -f12*a*(expei+LOG(DHsbw)-LOG(Hsbw))

       term1d1  = - a/(two*DHsbw) - f98*expei
       d1sterm1 = d1sHsbw*term1d1
       d1rterm1 = d1rHsbw*term1d1

       Fx_wpbe = x * (term1 + term2 + term3 + term4 + term5)

       d1sfx = x * (d1sterm1 + d1sterm2 + d1sterm3&
            + d1sterm4 + d1sterm5)

       d1rfx = x * (d1rterm1 + d1rterm3 + d1rterm4 + d1rterm5)

    ELSE

       ! For everything else, use the full blown expression

       ! First, we calculate the polynomials for the first term

       np1    = -f32*ea1*a12*w + r27*ea3*w3/(Eight*a12)&
            - r243*ea5*w5/(r32*A32) + r2187*ea7*w7/(r128*A52)

       d1rnp1 = - f32*ea1*d1rw*a12 + (r81*ea3*d1rw*w2)/(Eight*a12)&
            - (r1215*ea5*d1rw*w4)/(r32*a32)&
            + (r15309*ea7*d1rw*w6)/(r128*a52)

       np2 = -a + f94*ea2*w2 - r81*ea4*w4/(Sixteen*a)&
            + r729*ea6*w6/(r64*A2) - r6561*ea8*w8/(r256*A3)

       d1rnp2 =   f12*(Nine*ea2*d1rw*w)&
            - (r81*ea4*d1rw*w3)/(Four*a)&
            + (r2187*ea6*d1rw*w5)/(r32*a2)&
            - (r6561*ea8*d1rw*w7)/(r32*a3)

       ! The first term is

       t1    = f12*(np1*piexperf + np2*expei)
       d1st1 = f12*(d1spiexperf*np1 + d1sexpei*np2)
       d1rt1 = f12*(d1rnp2*expei + d1rpiexperf*np1 +&
            d1rexpei*np2 + d1rnp1*piexperf)

       ! The factors for the main polynomoal in w and their derivatives

       f2    = (f12)*ea1*srpi*a / DHsbw12
       f2d1  = - ea1*srpi*a / (Four*DHsbw32)
       d1sf2 = d1sHsbw*f2d1
       d1rf2 = d1rHsbw*f2d1

       f3    = (f12)*ea2*a / DHsbw
       f3d1  = - ea2*a / (two*DHsbw2)
       d1sf3 = d1sHsbw*f3d1
       d1rf3 = d1rHsbw*f3d1

       f4    =  ea3*srpi*(-f98 / Hsbw12&
            + f14*A / DHsbw32)
       f4d1  = ea3*srpi*((Nine/(Sixteen*Hsbw32))-&
            (Three*A/(Eight*DHsbw52)))
       d1sf4 = d1sHsbw*f4d1
       d1rf4 = d1rHsbw*f4d1

       f5    = ea4*(One/r128) * (-r144*(one/Hsbw)&
            + r64*(One/DHsbw2)*A)
       f5d1  = ea4*((f98/Hsbw2)-(a/DHsbw3))
       d1sf5 = d1sHsbw*f5d1
       d1rf5 = d1rHsbw*f5d1

       f6    = ea5*(Three*srpi*(three*DHsbw52*(Nine*Hsbw-two*a)&
            + Four*Hsbw32*A2))&
            / (r32*DHsbw52*Hsbw32*A)
       f6d1  = ea5*srpi*((r27/(r32*Hsbw52))-&
            (r81/(r64*Hsbw32*A))-&
            ((Fifteen*A)/(Sixteen*DHsbw72)))
       d1sf6 = d1sHsbw*f6d1
       d1rf6 = d1rHsbw*f6d1

       f7    = ea6*(((r32*a)/DHsbw3&
            + (-r36 + (r81*s2*H)/A)/Hsbw2)) / r32
       d1sf7 = ea6*(Three*(r27*d1sH*DHsbw4*Hsbw*s2 +&
            Eight*d1sHsbw*A*(Three*DHsbw4 - Four*Hsbw3*A) + &
            r54*DHsbw4*s*(Hsbw - d1sHsbw*s)*H))/&
            (r32*DHsbw4*Hsbw3*A)
       d1rf7 = ea6*d1rHsbw*((f94/Hsbw3)-((Three*a)/DHsbw4)&
            -((r81*s2*H)/(Sixteen*Hsbw3*A)))

       f8    = ea7*(-Three*srpi*(-r40*Hsbw52*a3&
            +Nine*DHsbw72*(r27*Hsbw2-Six*Hsbw*A+Four*A2)))&
            / (r128 * DHsbw72*Hsbw52*A2)
       f8d1  = ea7*srpi*((r135/(r64*Hsbw72)) + (r729/(r256*Hsbw32*a2))&
            -(r243/(r128*Hsbw52*A)) &
            -((r105*A)/(r32*DHsbw92)))
       d1sf8 = d1sHsbw*f8d1
       d1rf8 = d1rHsbw*f8d1

       f9    = (r324*ea6*eb1*DHsbw4*Hsbw*a&
            + ea8*(r384*Hsbw3*A3 + DHsbw4*(-r729*Hsbw2&
            + r324*Hsbw*A - r288*A2))) / (r128*DHsbw4*Hsbw3*A2)
       f9d1  = -((r81*ea6*eb1)/(Sixteen*Hsbw3*a))&
            + ea8*((r27/(Four*Hsbw4))+(r729/(r128*Hsbw2*A2))&
            -(r81/(Sixteen*Hsbw3*A))&
            -((r12*A/DHsbw5)))
       d1sf9 = d1sHsbw*f9d1
       d1rf9 = d1rHsbw*f9d1

       t2t9    = f2*w  + f3*w2 + f4*w3 + f5*w4 + f6*w5&
            + f7*w6 + f8*w7 + f9*w8
       d1st2t9 = d1sf2*w + d1sf3*w2 + d1sf4*w3 + d1sf5*w4&
            + d1sf6*w5 + d1sf7*w6 + d1sf8*w7&
            + d1sf9*w8
       d1rt2t9 = d1rw*f2 + d1rf2*w + two*d1rw*f3*w&
            + d1rf3*w2 + Three*d1rw*f4*w2&
            + d1rf4*w3 + Four*d1rw*f5*w3&
            + d1rf5*w4 + Five*d1rw*f6*w4&
            + d1rf6*w5 + Six*d1rw*f7*w5&
            + d1rf7*w6 + Seven*d1rw*f8*w6&
            + d1rf8*w7 + Eight*d1rw*f9*w7 + d1rf9*w8

       ! The final value of term1 for 0 < omega < wcutoff is:

       term1 = t1 + t2t9 + t10

       d1sterm1 = d1st1 + d1st2t9 + d1st10
       d1rterm1 = d1rt1 + d1rt2t9 + d1rt10

       ! The final value for the enhancement factor and its
       ! derivatives is:

       Fx_wpbe = x * (term1 + term2 + term3 + term4 + term5)

       d1sfx = x * (d1sterm1 + d1sterm2 + d1sterm3&
            + d1sterm4 + d1sterm5)

       d1rfx = x * (d1rterm1 + d1rterm3 + d1rterm4 + d1rterm5)

    ENDIF

  END SUBROUTINE wpbe_analytical_erfc_approx_grad


  ! ==================================================================
  FUNCTION expint(n,x)
    ! ==================================================================
    INTEGER                                  :: n
    REAL(real_8)                             :: x, expint

    INTEGER, PARAMETER                       :: maxit = 100 
    REAL(real_8), PARAMETER :: eps = 1.e-9_real_8, &
      euler = 0.5772156649015328_real_8 , fpmin = 1.e-30_real_8 

    INTEGER                                  :: i, j, n1
    REAL(real_8)                             :: a, b, c, d, del, h, vact, vsi

! =-------------------------------------------------------------------=

    IF ((n<0).OR.(x<0.0_real_8).OR.((x==0._real_8).AND.(n==0.OR.n==1))) THEN
       CALL stopgm('EXPINT',"Passing bad arguments in EXPINT! ",& 
            __LINE__,__FILE__)
    ENDIF

    n1=n-1
    IF (n==0) THEN
       expint=EXP(-x)/x
    ELSE
       del=HUGE(0.0_real_8)
       IF (x==0.0_real_8) THEN
          expint=1._real_8/n1
       ELSE IF (x>1.0_real_8) THEN
          b=x+n
          c=1._real_8/fpmin
          d=1.0_real_8/b
          h=d
          i=0
          DO WHILE ((i.LE.maxit).AND.(ABS(del-1._real_8)>=eps))
             i=i+1
             a=-i*(n1+i)
             b=b+2._real_8
             d=1.0_real_8/(a*d+b)
             c=b+a/c
             del=c*d
             h=h*del
          ENDDO
          IF (i.GT.maxit) CALL stopgm('EXPINT',&
               "Error in EXPINT! (X>1.0) : Not converged!",& 
               __LINE__,__FILE__)
          expint=h*EXP(-x)
       ELSE
          IF (n1==0) THEN
             expint=-LOG(x)-euler
          ELSE
             expint=1.0_real_8/n1
          ENDIF
          vact=1._real_8
          i=0
          DO WHILE ((i.LE.maxit).AND.(ABS(del)>=ABS(expint)*eps))
             i=i+1
             vact=-vact*x/i
             IF (i==n1) THEN
                vsi=-euler
                DO j=1,n1
                   vsi=vsi+1._real_8/j
                ENDDO
                del=vact*(-LOG(x)+vsi)
             ELSE
                del=-vact/(i-n1)
             ENDIF
             expint=expint+del
          ENDDO
          IF (i.GT.maxit) CALL stopgm('EXPINT',&
               "Error in EXPINT! (X<1.0) : Not converged!",& 
               __LINE__,__FILE__)
       ENDIF
    ENDIF
  END FUNCTION expint
  ! =-------------------------------------------------------------------=

END MODULE functionals_utils
