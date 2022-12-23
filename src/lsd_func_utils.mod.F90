#ifdef __SR8000
!option MP(P(0)), LANGLVL(SAVE(0))
#endif

MODULE lsd_func_utils
  USE cnst,                            ONLY: fpi,&
                                             pi
  USE error_handling,                  ONLY: stopgm
  USE func,                            ONLY: &
       func1, func2, func3, mfxcc_is_lyp, mfxcc_is_pade, mfxcc_is_pz, &
       mfxcc_is_skipped, mfxcx_is_slaterx, mgcc_is_dfr_zpbec, mgcc_is_hse, &
       mgcc_is_lyp, mgcc_is_pbec, mgcc_is_pbesc, mgcc_is_perdew86, &
       mgcc_is_skipped, mgcx_is_becke88, mgcx_is_dfr_xpbex, &
       mgcx_is_dfr_xpbex_hybrid, mgcx_is_dfr_zpbex, mgcx_is_ggax, &
       mgcx_is_hcth, mgcx_is_optx, mgcx_is_ox, mgcx_is_ox_hybrid, &
       mgcx_is_pbesx, mgcx_is_pbex, mgcx_is_revpbex, mgcx_is_skipped
  USE functionals_utils,               ONLY: ggax
  USE kinds,                           ONLY: real_8
  USE system,                          ONLY: cntl

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: xc_lsd
  PUBLIC :: gc_lsd
  PUBLIC :: lsd_sx
  PUBLIC :: lsd_pz
  PUBLIC :: lsd_lyp
  PUBLIC :: lsd_pade
  PUBLIC :: lsd_b88
  PUBLIC :: lsd_ggax
  PUBLIC :: lsd_pbesx
  PUBLIC :: lsd_pbex
  PUBLIC :: lsd_revpbex
  PUBLIC :: lsd_p86
  PUBLIC :: lsd_glyp
  PUBLIC :: lsd_pbesc
  PUBLIC :: lsd_pbec
  PUBLIC :: lsd_hcth
  PUBLIC :: pwcorr2
  PUBLIC :: lsd_optx

CONTAINS

  ! ==================================================================
  SUBROUTINE xc_lsd(rho,eta,ex,ec,vxa,vca,vxb,vcb)
    ! ==--------------------------------------------------------------==
    ! ==  LSD EXCHANGE AND CORRELATION FUNCTIONALS                    ==
    ! ==                                                              ==
    ! ==  EXCHANGE  :  SLATER alpha                                   ==
    ! ==               OPTX                                           ==
    ! ==  CORRELATION : PERDEW & ZUNGER                               ==
    ! ==                VOSKO, WILK & NUSSAIR                         ==
    ! ==                LEE, YANG & PARR                              ==
    ! ==                PERDEW & WANG                                 ==
    ! ==                HCTH/120                                      ==
    ! ==--------------------------------------------------------------==

    IMPLICIT NONE

    REAL(real_8)                             :: rho, eta, ex, ec, vxa, vca, &
         vxb, vcb

    CHARACTER(len=*), PARAMETER              :: procedureN = 'xc_lsd'
    REAL(real_8), PARAMETER                  :: pi34 = 0.75_real_8/pi, &
         small = 1.e-20_real_8 , &
         third = 1._real_8/3._real_8 

    INTEGER                                  :: iflg
    REAL(real_8)                             :: rs

    ! ==--------------------------------------------------------------==
    ! ..Exchange

    IF (func1%mfxcx == mfxcx_is_slaterx .AND. rho.GT.small) THEN
       CALL lsd_sx(rho,eta,ex,vxa,vxb,func2%salpha)
       ex=ex*func3%pxlda
       vxa=vxa*func3%pxlda
       vxb=vxb*func3%pxlda
    ELSE
       ex=0.0_real_8
       vxa=0.0_real_8
       vxb=0.0_real_8
    ENDIF
    ! ..Correlation
    IF (rho.LE.small) THEN
       ec  = 0.0_real_8
       vca = 0.0_real_8
       vcb = 0.0_real_8
       ex  = 0.0_real_8
       vxa = 0.0_real_8
       vxb = 0.0_real_8
    ELSE
       ! Quick LYP correlation fix. thybrid is included for safety reasons (don't want to break anything)
       IF (cntl%thybrid .AND. func1%mgcc==2 .AND. .NOT. func1%mfxcc==3) THEN
          CALL mixed_lsd_correlation()
          ! This is the old loop. Unchanged except for the weights.
       ELSE
          IF (func1%mfxcc == mfxcc_is_skipped) THEN
             ec  = 0.0_real_8
             vca = 0.0_real_8
             vcb = 0.0_real_8
          ELSE IF (func1%mfxcc == mfxcc_is_pz) THEN
             rs=(pi34/rho)**third
             iflg=2
             IF (rs.LT.1.0_real_8) iflg=1
             CALL lsd_pz(rs,eta,ec,vca,vcb,iflg)
          ELSEIF (func1%mfxcc == mfxcc_is_lyp) THEN
             CALL lsd_lyp(rho,eta,ec,vca,vcb)
          ELSEIF (func1%mfxcc == mfxcc_is_pade) THEN
             CALL lsd_pade(rho,eta,ec,vca,vcb)
          ELSE
             ! vw MFXCC should be valide at this point!
             CALL stopgm(procedureN,"correlation functional "//&
                  "not implemented",& 
                  __LINE__,__FILE__)
          ENDIF
          ec = ec*func3%pclda
          vca = vca*func3%pclda ; vcb = vcb*func3%pclda
       END IF
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN

  CONTAINS

    !-------------------------------------------------------------------
    SUBROUTINE mixed_lsd_correlation()
      !-------------------------------------------------------------------
    CHARACTER(len=*), PARAMETER :: procedureN = 'mixed_lsd_correlation'

    REAL(real_8)                             :: ec_aux, ec_lyp, vca_aux, &
                                                vca_lyp, vcb_aux, vcb_lyp

!------------------------------------------------------------------

      ec_aux=0.0_real_8 ; ec_lyp=0.0_real_8
      vca_aux=0.0_real_8 ; vcb_aux=0.0_real_8
      vca_lyp=0.0_real_8 ; vcb_lyp=0.0_real_8
      IF (func1%mfxcc == mfxcc_is_skipped) THEN
         ec  = 0.0_real_8
         vca = 0.0_real_8
         vcb = 0.0_real_8
      ELSE IF (func1%mfxcc == mfxcc_is_pz) THEN
         rs=(pi34/rho)**third
         iflg=2
         IF (rs.LT.1.0_real_8) iflg=1
         CALL lsd_pz(rs,eta,ec_aux,vca_aux,vcb_aux,iflg)
      ELSEIF (func1%mfxcc == mfxcc_is_pade) THEN
         CALL lsd_pade(rho,eta,ec_aux,vca_aux,vcb_aux)
      ELSE
         ! vw MFXCC should be valide at this point!
         CALL stopgm(procedureN,"correlation functional "//&
              "not implemented",& 
              __LINE__,__FILE__)
      ENDIF
      ! Calculate LYP contribution and add it to the real LDA part.
      CALL lsd_lyp(rho,eta,ec_lyp,vca_lyp,vcb_lyp)
      ec=ec_aux*func3%pclda
      vca=vca_aux*func3%pclda ; vcb=vcb_aux*func3%pclda
      ec=ec+ec_lyp*func3%pcgc
      vca=vca+vca_lyp*func3%pcgc ; vcb=vcb+vcb_lyp*func3%pcgc
      RETURN
      !-------------------------------------------------------------------
    END SUBROUTINE mixed_lsd_correlation
    !-------------------------------------------------------------------

  END SUBROUTINE xc_lsd
  ! ==================================================================
  SUBROUTINE gc_lsd(rhoa,rhob,grhoaa,grhobb,grhoab,sx,sc,v1xa,v2xa,&
       v1xb,v2xb,v1ca,v2ca,v1cb,v2cb,v2xab,v2cab)
    ! ==--------------------------------------------------------------==
    ! ==  GRADIENT CORRECTIONS FOR EXCHANGE AND CORRELATION           ==
    ! ==                                                              ==
    ! ==  EXCHANGE  :  BECKE88                                        ==
    ! ==               GGAX                                           ==
    ! ==               PBEX                                           ==
    ! ==               revPBEX                                        ==
    ! ==               HCTH/120                                       ==
    ! ==               OPTX                                           ==
    ! ==  CORRELATION : PERDEW86                                      ==
    ! ==                LEE, YANG & PARR                              ==
    ! ==                GGAC                                          ==
    ! ==                PBEC                                          ==
    ! ==                HCTH/120                                      ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rhoa, rhob, grhoaa, grhobb, &
                                                grhoab, sx, sc, v1xa, v2xa, &
                                                v1xb, v2xb, v1ca, v2ca, v1cb, &
                                                v2cb, v2xab, v2cab

    CHARACTER(len=*), PARAMETER              :: procedureN = 'gc_lsd'
    REAL(real_8), PARAMETER                  :: small = 1.e-20_real_8 

    REAL(real_8)                             :: grho, rho, sxa, sxb, v1xaa, &
                                                v1xab, v1xba, v1xbb, v2xaa, &
                                                v2xba, v2xbb

! ==--------------------------------------------------------------==
! ..Exchange

    IF (func1%mgcx == mgcx_is_skipped) THEN
       sx=0.0_real_8
       v1xa=0.0_real_8
       v2xa=0.0_real_8
       v1xb=0.0_real_8
       v2xb=0.0_real_8
       v2xab=0.0_real_8
    ELSEIF (func1%mgcx == mgcx_is_becke88) THEN
       CALL lsd_b88(func2%bbeta,rhoa,rhob,grhoaa,grhobb,&
            sx,v1xa,v2xa,v1xb,v2xb)
       v2xab=0.0_real_8
    ELSEIF (func1%mgcx == mgcx_is_ggax) THEN
       CALL lsd_ggax(rhoa,rhob,grhoaa,grhobb,sx,v1xa,v2xa,v1xb,v2xb)
       v2xab=0.0_real_8
    ELSEIF (func1%mgcx == mgcx_is_pbex) THEN
       CALL lsd_pbex(rhoa,rhob,grhoaa,grhobb,sx,v1xa,v2xa,v1xb,v2xb)
       v2xab=0.0_real_8
    ELSEIF (func1%mgcx == mgcx_is_revpbex) THEN
       CALL lsd_revpbex(rhoa,rhob,grhoaa,grhobb,sx,v1xa,v2xa,v1xb,v2xb)
       v2xab=0.0_real_8
    ELSEIF (func1%mgcx == mgcx_is_hcth.AND.func1%mgcc == mgcc_is_hse) THEN
       CALL lsd_hcth(rhoa,rhob,grhoaa,grhobb,sx,v1xa,v2xa,v1xb,v2xb)
       sc=0
       v1ca=0.0_real_8
       v2ca=0.0_real_8
       v1cb=0.0_real_8
       v2cb=0.0_real_8
       v2xab=0.0_real_8
       v2cab=0.0_real_8
    ELSEIF (func1%mgcx == mgcx_is_optx) THEN
       CALL lsd_optx(rhoa,rhob,grhoaa,grhobb,sx,v1xa,v2xa,v1xb,v2xb)
       v2xab=0.0_real_8
    ELSEIF (func1%mgcx == mgcx_is_ox) THEN
       CALL lsd_b88(func2%bbeta,rhoa,rhob,grhoaa,grhobb,&
            sxa,v1xaa,v2xaa,v1xba,v2xba)
       CALL lsd_ggax(rhoa,rhob,grhoaa,grhobb,&
            sxb,v1xab,v2xab,v1xbb,v2xbb)
       sx=0.722_real_8*sxa+0.347_real_8*sxb
       v1xa=0.722_real_8*v1xaa+0.347_real_8*v1xab
       v2xa=0.722_real_8*v2xaa+0.347_real_8*v2xab
       v1xb=0.722_real_8*v1xba+0.347_real_8*v1xbb
       v2xb=0.722_real_8*v2xba+0.347_real_8*v2xbb
    ELSEIF (func1%mgcx == mgcx_is_ox_hybrid) THEN
       CALL lsd_b88(func2%bbeta,rhoa,rhob,grhoaa,grhobb,&
            sxa,v1xaa,v2xaa,v1xba,v2xba)
       CALL lsd_ggax(rhoa,rhob,grhoaa,grhobb,&
            sxb,v1xab,v2xab,v1xbb,v2xbb)
       sx=0.542_real_8*sxa+0.167_real_8*sxb
       v1xa=0.542_real_8*v1xaa+0.167_real_8*v1xab
       v2xa=0.542_real_8*v2xaa+0.167_real_8*v2xab
       v1xb=0.542_real_8*v1xba+0.167_real_8*v1xbb
       v2xb=0.542_real_8*v2xba+0.167_real_8*v2xbb
    ELSEIF (func1%mgcx == mgcx_is_pbesx) THEN
       CALL lsd_pbesx(rhoa,rhob,grhoaa,grhobb,sx,v1xa,v2xa,v1xb,v2xb)
       v2xab=0.0_real_8
    ELSEIF (func1%mgcx == mgcx_is_dfr_zpbex) THEN
       CALL stopgm(procedureN,'you need the DFRepository.F',&
            __LINE__,__FILE__)
    ELSEIF (func1%mgcx == mgcx_is_dfr_xpbex) THEN
       CALL stopgm(procedureN,'you need the DFRepository.F',&
            __LINE__,__FILE__)
    ELSEIF (func1%mgcx == mgcx_is_dfr_xpbex_hybrid) THEN
       CALL stopgm(procedureN,'you need the DFRepository.F',&
            __LINE__,__FILE__)
    ELSE
       ! vw MFXCX should be valide at this point!
       CALL stopgm(procedureN,"exchange functional "//&
            "not implemented",& 
            __LINE__,__FILE__)
    ENDIF
    ! ..Correlation
    sc=0.0_real_8
    v1ca=0.0_real_8
    v2ca=0.0_real_8
    v1cb=0.0_real_8
    v2cb=0.0_real_8
    v2cab=0.0_real_8
    IF (func1%mgcc == mgcc_is_skipped) THEN
       ! 
    ELSEIF (func1%mgcc == mgcc_is_perdew86) THEN
       rho=rhoa+rhob
       grho=grhoaa+2._real_8*grhoab+grhobb
       IF (rho.GT.small)&
            CALL lsd_p86(rhoa,rhob,grho,sc,v1ca,v2ca,v1cb,v2cb,v2cab)
    ELSEIF (func1%mgcc == mgcc_is_lyp) THEN
       rho=rhoa+rhob
       IF (rho.GT.small)&
            CALL lsd_glyp(rhoa,rhob,grhoaa,grhoab,grhobb,sc,&
            v1ca,v2ca,v1cb,v2cb,v2cab)
       ! mb   ELSEIF(MGCC.EQ.3) THEN
       ! mb     STOP "LSD_GGAC"
    ELSEIF (func1%mgcc == mgcc_is_hse) THEN
       ! Just do nothing - already taken into account in the exchange part above
    ELSEIF (func1%mgcc == mgcc_is_pbec) THEN
       rho=rhoa+rhob
       IF (rho.GT.small)&
            CALL lsd_pbec(rhoa,rhob,grhoaa,grhoab,grhobb,sc,&
            v1ca,v2ca,v1cb,v2cb,v2cab)
    ELSEIF (func1%mgcc == mgcc_is_pbesc) THEN
       rho=rhoa+rhob
       IF (rho.GT.small)&
            CALL lsd_pbesc(rhoa,rhob,grhoaa,grhoab,grhobb,sc,&
            v1ca,v2ca,v1cb,v2cb,v2cab)
    ELSEIF (func1%mgcc == mgcc_is_dfr_zpbec) THEN
       CALL stopgm(procedureN,'you need the DFRepository.F',&
            __LINE__,__FILE__)
    ELSE
       ! vw MFCC should be valide at this point!
       CALL stopgm(procedureN,"correlation functional "//&
            "not implemented",& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gc_lsd
  ! ==================================================================
  SUBROUTINE lsd_sx(rho,eta,ex,vxa,vxb,alpha)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rho, eta, ex, vxa, vxb, alpha

    REAL(real_8), PARAMETER :: f = -1.39578858466194911_real_8 , &
      f43 = 4._real_8/3._real_8 , third = 1._real_8/3._real_8

    REAL(real_8)                             :: rhoa, rhob, rsa, rsb

! ==--------------------------------------------------------------==

    rhoa=0.5_real_8*rho*(1._real_8+eta)
    rhob=0.5_real_8*rho*(1._real_8-eta)
    rsa = rhoa**third
    rsb = rhob**third
    ex = f*alpha*(rhoa*rsa+rhob*rsb)/rho
    vxa = f43*f*alpha*rsa
    vxb = f43*f*alpha*rsb
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE lsd_sx
  ! ==================================================================
  SUBROUTINE lsd_pz(rs,eta,epz,vpza,vpzb,iflg)
    ! ==--------------------------------------------------------------==
    ! ==  J.P. PERDEW AND ALEX ZUNGER PRB 23, 5048 (1981)             ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rs, eta, epz, vpza, vpzb
    INTEGER                                  :: iflg

    REAL(real_8), PARAMETER :: ap = 0.01555_real_8, au = 0.0311_real_8, &
      b1p = 1.3981_real_8, b1u = 1.0529_real_8, b2p = 0.2611_real_8 , &
      b2u = 0.3334_real_8 , bp = -0.0269_real_8, bu = -0.048_real_8, &
      cp = 0.0007_real_8, cu = 0.0020_real_8, dp = -0.0048_real_8, &
      du = -0.0116_real_8, f13 = 1._real_8/3._real_8, &
      f43 = 4._real_8/3._real_8 , gcp = -0.0843_real_8, gcu = -0.1423_real_8

    REAL(real_8)                             :: dfeta, doxp, doxu, ep, eu, &
                                                feta, oxp, oxu, rs1, rs2, vp, &
                                                vu, xln

! ==--------------------------------------------------------------==

    IF (iflg.EQ.1) THEN
       ! ..High density formula
       xln=LOG(rs)
       eu=au*xln+bu+cu*rs*xln+du*rs
       ep=ap*xln+bp+cp*rs*xln+dp*rs
       vu=au*xln+(bu-au/3._real_8)+2._real_8/3._real_8*cu*rs*xln+&
            (2._real_8*du-cu)/3._real_8*rs
       vp=ap*xln+(bp-ap/3._real_8)+2._real_8/3._real_8*cp*rs*xln+&
            (2._real_8*dp-cp)/3._real_8*rs
    ELSEIF (iflg.EQ.2) THEN
       ! ..Interpolation formula
       rs1=SQRT(rs)
       rs2=rs
       oxu=1._real_8+b1u*rs1+b2u*rs2
       doxu=1._real_8+7._real_8/6._real_8*b1u*rs1+4._real_8/3._real_8*b2u*rs2
       eu=gcu/oxu
       vu=eu*doxu/oxu
       oxp=1._real_8+b1p*rs1+b2p*rs2
       doxp=1._real_8+7._real_8/6._real_8*b1p*rs1+4._real_8/3._real_8*b2p*rs2
       ep=gcp/oxp
       vp=ep*doxp/oxp
    ENDIF
    feta=((1._real_8+eta)**f43+(1._real_8-eta)**f43-2._real_8)/(2._real_8**f43-2._real_8)
    dfeta=f43*((1._real_8+eta)**f13-(1._real_8-eta)**f13)/(2._real_8**f43-2._real_8)
    epz=eu+feta*(ep-eu)
    vpza=vu+feta*(vp-vu)+(ep-eu)*(1._real_8-eta)*dfeta
    vpzb=vu+feta*(vp-vu)+(ep-eu)*(-1._real_8-eta)*dfeta
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE lsd_pz
  ! ==================================================================
  SUBROUTINE lsd_lyp(rho,eta,elyp,valyp,vblyp)
    ! ==--------------------------------------------------------------==
    ! ==  C. LEE, W. YANG, AND R.G. PARR, PRB 37, 785 (1988)          ==
    ! ==  THIS IS ONLY THE LDA PART                                   ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rho, eta, elyp, valyp, vblyp

    REAL(real_8), PARAMETER :: a = 0.04918_real_8, b = 0.132_real_8, &
      c = 0.2533_real_8, cf = 2.87123400018819108_real_8 , d = 0.349_real_8 , &
      small = 1.e-24_real_8 

    REAL(real_8)                             :: de1a, de1b, de2a, de2b, dor, &
                                                dr, e1, e2, or, ra, rb, rm3

! ==--------------------------------------------------------------==

    ra=rho*0.5_real_8*(1._real_8+eta)
    ra=MAX(ra,small)
    rb=rho*0.5_real_8*(1._real_8-eta)
    rb=MAX(rb,small)
    rm3=rho**(-1._real_8/3._real_8)
    dr=(1._real_8+d*rm3)
    e1=4._real_8*a*ra*rb/rho/dr
    or=EXP(-c*rm3)/dr*rm3**11._real_8
    dor=-1._real_8/3._real_8*rm3**4*or*(11._real_8/rm3-c-d/dr)
    e2=2._real_8**(11._real_8/3._real_8)*cf*a*b*or*ra*rb*(ra**(8._real_8/3._real_8)+&
         rb**(8._real_8/3._real_8))
    elyp=(-e1-e2)/rho
    de1a=-e1*(1._real_8/3._real_8*d*rm3**4/dr+1._real_8/ra-1._real_8/rho)
    de1b=-e1*(1._real_8/3._real_8*d*rm3**4/dr+1._real_8/rb-1._real_8/rho)
    de2a=-2._real_8**(11._real_8/3._real_8)*cf*a*b*(dor*ra*rb*(ra**(8._real_8/3._real_8)+&
         rb**(8._real_8/3._real_8))+or*rb*(11._real_8/3._real_8*ra**(8._real_8/3._real_8)+&
         rb**(8._real_8/3._real_8)))
    de2b=-2._real_8**(11._real_8/3._real_8)*cf*a*b*(dor*ra*rb*(ra**(8._real_8/3._real_8)+&
         rb**(8._real_8/3._real_8))+or*ra*(11._real_8/3._real_8*rb**(8._real_8/3._real_8)+&
         ra**(8._real_8/3._real_8)))
    valyp=de1a+de2a
    vblyp=de1b+de2b
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE lsd_lyp
  ! ==================================================================
  SUBROUTINE lsd_pade(rho,eta,ec,vca,vcb)
    ! ==--------------------------------------------------------------==
    ! ==  PADE APPROXIMATION                                          ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rho, eta, ec, vca, vcb

    REAL(real_8), PARAMETER :: a0 = 0.4581652932831429_real_8, &
      a1 = 2.217058676663745_real_8, a2 = 0.7405551735357053_real_8, &
      a3 = 0.01968227878617998_real_8 , b1 = 1.0_real_8, &
      b2 = 4.504130959426697_real_8, b3 = 1.110667363742916_real_8, &
      b4 = 0.02359291751427506_real_8 , da0 = 0.119086804055547_real_8, &
      da1 = 0.6157402568883345_real_8, da2 = 0.1574201515892867_real_8, &
      da3 = 0.003532336663397157_real_8 , db1 = 0.0_real_8, &
      db2 = 0.2673612973836267_real_8, db3 = 0.2052004607777787_real_8, &
      db4 = 0.004200005045691381_real_8 , &
      fsfac = 1.92366105093153617_real_8 , rsfac = 0.6203504908994000_real_8 

    REAL(real_8)                             :: a0p, a1p, a2p, a3p, b1p, b2p, &
                                                b3p, b4p, bot, botx, dbot, &
                                                dfs, dfsa, dfsb, dtop, dx, &
                                                fs, rs, top, topx, vc

! ==--------------------------------------------------------------==

    rs=rsfac*rho**(-1._real_8/3._real_8)
    fs=fsfac*((1._real_8+eta)**(4._real_8/3._real_8)+(1._real_8-eta)**(4._real_8/3._real_8)-2._real_8)
    dfs=fsfac*4._real_8/3._real_8*&
         ((1._real_8+eta)**(1._real_8/3._real_8)-(1._real_8-eta)**(1._real_8/3._real_8))
    dfsa=dfs*(1._real_8-eta)
    dfsb=dfs*(-1._real_8-eta)
    a0p=a0+fs*da0
    a1p=a1+fs*da1
    a2p=a2+fs*da2
    a3p=a3+fs*da3
    b1p=b1+fs*db1
    b2p=b2+fs*db2
    b3p=b3+fs*db3
    b4p=b4+fs*db4
    top=a0p+rs*(a1p+rs*(a2p+rs*a3p))
    dtop=a1p+rs*(2._real_8*a2p+rs*3._real_8*a3p)
    topx=da0+rs*(da1+rs*(da2+rs*da3))
    bot=rs*(b1p+rs*(b2p+rs*(b3p+rs*b4p)))
    dbot=b1p+rs*(2._real_8*b2p+rs*(3._real_8*b3p+rs*4._real_8*b4p))
    botx=rs*(db1+rs*(db2+rs*(db3+rs*db4)))
    ec=-top/bot
    vc=ec+rs*(dtop/bot-top*dbot/(bot*bot))/3._real_8
    dx=-(topx/bot-top*botx/(bot*bot))
    vca=vc+dx*dfsa
    vcb=vc+dx*dfsb
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE lsd_pade
  ! ==================================================================
  SUBROUTINE lsd_b88(b1,rhoa,rhob,grhoa,grhob,&
       sx,v1xa,v2xa,v1xb,v2xb)
    ! ==--------------------------------------------------------------==
    ! BECKE EXCHANGE: PRA 38, 3098 (1988)
    REAL(real_8)                             :: b1, rhoa, rhob, grhoa, grhob, &
                                                sx, v1xa, v2xa, v1xb, v2xb

    REAL(real_8), PARAMETER                  :: ob3 = 1._real_8/3._real_8, &
                                                small = 1.e-20_real_8 

    REAL(real_8)                             :: a, aa, br1, br2, br4, dd, &
                                                dd2, ddd, dgf, gf, sa2b8, &
                                                shm1, xs, xs2

! ==--------------------------------------------------------------==

    sx=0.0_real_8
    v1xa=0.0_real_8
    v2xa=0.0_real_8
    v1xb=0.0_real_8
    v2xb=0.0_real_8
    IF (ABS(rhoa).GT.small) THEN
       aa    = grhoa
       a     = SQRT(aa)
       a     = MAX(a,small)
       br1   = rhoa**ob3
       br2   = br1*br1
       br4   = br2*br2
       xs    = a/br4
       xs2   = xs*xs
       sa2b8 = SQRT(1.0_real_8+xs2)
       shm1  = LOG(xs+sa2b8)
       dd    = 1.0_real_8 + 6.0_real_8*b1*xs*shm1
       dd2   = dd*dd
       ddd   = 6.0_real_8*b1*(shm1+xs/sa2b8)
       gf    = -b1*xs2/dd
       dgf   = (-2.0_real_8*b1*xs*dd + b1*xs2*ddd)/dd2
       sx    = gf*br4
       v1xa  = 4._real_8/3._real_8*br1*(gf-xs*dgf)
       v2xa  = dgf/a
    ENDIF
    IF (ABS(rhob).GT.small) THEN
       aa    = grhob
       a     = SQRT(aa)
       a     = MAX(a,small)
       br1   = rhob**ob3
       br2   = br1*br1
       br4   = br2*br2
       xs    = a/br4
       xs2   = xs*xs
       sa2b8 = SQRT(1.0_real_8+xs2)
       shm1  = LOG(xs+sa2b8)
       dd    = 1.0_real_8 + 6.0_real_8*b1*xs*shm1
       dd2   = dd*dd
       ddd   = 6.0_real_8*b1*(shm1+xs/sa2b8)
       gf    = -b1*xs2/dd
       dgf   = (-2.0_real_8*b1*xs*dd + b1*xs2*ddd)/dd2
       sx    = sx+gf*br4
       v1xb  = 4._real_8/3._real_8*br1*(gf-xs*dgf)
       v2xb  = dgf/a
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE lsd_b88
  ! ==================================================================
  SUBROUTINE lsd_ggax(rhoa,rhob,grhoaa,grhobb,&
       sx,v1xa,v2xa,v1xb,v2xb)
    ! J.P.PERDEW ET AL. PRB 466671 (1992)
    REAL(real_8)                             :: rhoa, rhob, grhoaa, grhobb, &
                                                sx, v1xa, v2xa, v1xb, v2xb

    REAL(real_8)                             :: grho, rho, sx1, sx2

! ==--------------------------------------------------------------==

    rho=2._real_8*rhoa
    grho=4._real_8*grhoaa
    CALL ggax(rho,grho,sx1,v1xa,v2xa)
    rho=2._real_8*rhob
    grho=4._real_8*grhobb
    CALL ggax(rho,grho,sx2,v1xb,v2xb)
    sx = 0.5_real_8*sx1 + 0.5_real_8*sx2
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE lsd_ggax
  ! ==================================================================
  SUBROUTINE lsd_pbesx(rhoa,rhob,grhoa,grhob,sx,v1xa,v2xa,v1xb,v2xb)
    ! J.P.PERDEW ET AL.  (2007)
    REAL(real_8)                             :: rhoa, rhob, grhoa, grhob, sx, &
                                                v1xa, v2xa, v1xb, v2xb

    REAL(real_8), PARAMETER :: ax = -0.738558766382022406_real_8, &
      small = 1.e-20_real_8 , uk = 0.8040_real_8, &
      um = 0.1234567901234568_real_8, ul = um/uk , &
      us = 0.161620459673995492_real_8

    REAL(real_8)                             :: aa, dfx, ex, fx, po, rho, rr, &
                                                s2, sxa, sxb

! ==--------------------------------------------------------------==

    sxa=0.0_real_8
    sxb=0.0_real_8
    v1xa=0.0_real_8
    v2xa=0.0_real_8
    v1xb=0.0_real_8
    v2xb=0.0_real_8
    IF (ABS(rhoa).GT.small) THEN
       rho   = 2._real_8*rhoa
       aa    = 4._real_8*grhoa
       rr    = rho**(-4._real_8/3._real_8)
       ex    = ax/rr
       s2    = aa*rr*rr*us*us
       po    = 1._real_8/(1._real_8 + ul*s2)
       fx    = uk-uk*po
       sxa   = ex*fx
       dfx   = 2._real_8*uk*ul*po*po
       v1xa  = 1.33333333333333_real_8*ax*rho**0.333333333333_real_8*(fx-s2*dfx)
       v2xa  = ex*dfx*(us*rr)**2
    ENDIF
    IF (ABS(rhob).GT.small) THEN
       rho   = 2._real_8*rhob
       aa    = 4._real_8*grhob
       rr    = rho**(-4._real_8/3._real_8)
       ex    = ax/rr
       s2    = aa*rr*rr*us*us
       po    = 1._real_8/(1._real_8 + ul*s2)
       fx    = uk-uk*po
       sxb   = ex*fx
       dfx   = 2._real_8*uk*ul*po*po
       v1xb  = 1.33333333333333_real_8*ax*rho**0.333333333333_real_8*(fx-s2*dfx)
       v2xb  = ex*dfx*(us*rr)**2
    ENDIF
    sx    = 0.5_real_8*(sxa+sxb)
    v2xa  = 2._real_8*v2xa
    v2xb  = 2._real_8*v2xb
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE lsd_pbesx
  ! ==================================================================
  SUBROUTINE lsd_pbex(rhoa,rhob,grhoa,grhob,sx,v1xa,v2xa,v1xb,v2xb)
    ! J.P.PERDEW ET AL. PRL 773865 (1996)
    REAL(real_8)                             :: rhoa, rhob, grhoa, grhob, sx, &
                                                v1xa, v2xa, v1xb, v2xb

    REAL(real_8), PARAMETER :: ax = -0.738558766382022406_real_8, &
      small = 1.e-20_real_8 , uk = 0.8040_real_8, &
      um = 0.2195149727645171_real_8, ul = um/uk , &
      us = 0.161620459673995492_real_8

    REAL(real_8)                             :: aa, dfx, ex, fx, po, rho, rr, &
                                                s2, sxa, sxb

! ==--------------------------------------------------------------==

    sxa=0.0_real_8
    sxb=0.0_real_8
    v1xa=0.0_real_8
    v2xa=0.0_real_8
    v1xb=0.0_real_8
    v2xb=0.0_real_8
    IF (ABS(rhoa).GT.small) THEN
       rho   = 2._real_8*rhoa
       aa    = 4._real_8*grhoa
       rr    = rho**(-4._real_8/3._real_8)
       ex    = ax/rr
       s2    = aa*rr*rr*us*us
       po    = 1._real_8/(1._real_8 + ul*s2)
       fx    = uk-uk*po
       sxa   = ex*fx
       dfx   = 2._real_8*uk*ul*po*po
       v1xa  = 1.33333333333333_real_8*ax*rho**0.333333333333_real_8*(fx-s2*dfx)
       v2xa  = ex*dfx*(us*rr)**2
    ENDIF
    IF (ABS(rhob).GT.small) THEN
       rho   = 2._real_8*rhob
       aa    = 4._real_8*grhob
       rr    = rho**(-4._real_8/3._real_8)
       ex    = ax/rr
       s2    = aa*rr*rr*us*us
       po    = 1._real_8/(1._real_8 + ul*s2)
       fx    = uk-uk*po
       sxb   = ex*fx
       dfx   = 2._real_8*uk*ul*po*po
       v1xb  = 1.33333333333333_real_8*ax*rho**0.333333333333_real_8*(fx-s2*dfx)
       v2xb  = ex*dfx*(us*rr)**2
    ENDIF
    sx    = 0.5_real_8*(sxa+sxb)
    v2xa  = 2._real_8*v2xa
    v2xb  = 2._real_8*v2xb
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE lsd_pbex
  ! ==================================================================
  SUBROUTINE lsd_revpbex(RHOA,RHOB,GRHOA,GRHOB,SX,V1XA,V2XA,V1XB,V2XB)
    ! Y. ZHANG ET AL. PRL 80890 (1998)
    REAL(real_8)                             :: rhoa, rhob, grhoa, grhob, sx, &
                                                v1xa, v2xa, v1xb, v2xb

    REAL(real_8), PARAMETER :: ax = -0.738558766382022406_real_8, &
      small = 1.e-20_real_8 , uk = 1.2450_real_8, &
      um = 0.2195149727645171_real_8, ul = um/uk , &
      us = 0.161620459673995492_real_8

    REAL(real_8)                             :: aa, dfx, ex, fx, po, rho, rr, &
                                                s2, sxa, sxb

! ==--------------------------------------------------------------==

    sxa=0.0_real_8
    sxb=0.0_real_8
    v1xa=0.0_real_8
    v2xa=0.0_real_8
    v1xb=0.0_real_8
    v2xb=0.0_real_8
    IF (ABS(rhoa).GT.small) THEN
       rho   = 2._real_8*rhoa
       aa    = 4._real_8*grhoa
       rr    = rho**(-4._real_8/3._real_8)
       ex    = ax/rr
       s2    = aa*rr*rr*us*us
       po    = 1._real_8/(1._real_8 + ul*s2)
       fx    = uk-uk*po
       sxa   = ex*fx
       dfx   = 2._real_8*uk*ul*po*po
       v1xa  = 1.33333333333333_real_8*ax*rho**0.333333333333_real_8*(fx-s2*dfx)
       v2xa  = ex*dfx*(us*rr)**2
    ENDIF
    IF (ABS(rhob).GT.small) THEN
       rho   = 2._real_8*rhob
       aa    = 4._real_8*grhob
       rr    = rho**(-4._real_8/3._real_8)
       ex    = ax/rr
       s2    = aa*rr*rr*us*us
       po    = 1._real_8/(1._real_8 + ul*s2)
       fx    = uk-uk*po
       sxb   = ex*fx
       dfx   = 2._real_8*uk*ul*po*po
       v1xb  = 1.33333333333333_real_8*ax*rho**0.333333333333_real_8*(fx-s2*dfx)
       v2xb  = ex*dfx*(us*rr)**2
    ENDIF
    sx    = 0.5_real_8*(sxa+sxb)
    v2xa  = 2._real_8*v2xa
    v2xb  = 2._real_8*v2xb
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE lsd_revpbex
  ! ==================================================================
  SUBROUTINE lsd_p86(rhoa,rhob,grho,sc,v1ca,v2ca,v1cb,v2cb,v2cab)
    ! ==--------------------------------------------------------------==
    ! PERDEW CORRELATION: PRB 33, 8822 (1986)
    REAL(real_8)                             :: rhoa, rhob, grho, sc, v1ca, &
                                                v2ca, v1cb, v2cb, v2cab

    REAL(real_8), PARAMETER :: ob3 = 1._real_8/3._real_8 , &
      p1 = 0.023266_real_8, p2 = 7.389e-6_real_8, p3 = 8.723_real_8, &
      p4 = 0.472_real_8 , pc1 = 0.001667_real_8, pc2 = 0.002568_real_8, &
      pci = pc1+pc2 

    REAL(real_8)                             :: a, aa, br1, br2, br4, cn, &
                                                cna, cnb, d, dcn, dcna, dcnb, &
                                                dda, ddb, drs, ephi, phi, &
                                                rho, rs, rs2, rs3, s1, s2, &
                                                v1c, v2c

! ==--------------------------------------------------------------==

    rho=rhoa+rhob
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
    s1    = SQRT(rhoa**(5._real_8/3._real_8)+rhob**(5._real_8/3._real_8))
    s2    = 2._real_8**ob3/rho**(5._real_8/6._real_8)
    d     = s1*s2
    dda   = s2*5._real_8/6._real_8*(rhoa**(2._real_8/3._real_8)/s1-s1/rho)
    ddb   = s2*5._real_8/6._real_8*(rhob**(2._real_8/3._real_8)/s1-s1/rho)
    phi   = 0.192_real_8*pci/cn*a*rho**(-7._real_8/6._real_8)
    ephi  = EXP(-phi)
    sc    = aa/br4*cn*ephi/d
    v1c   = sc*((1._real_8+phi)*dcn/cn -((4._real_8/3._real_8)-(7._real_8/6._real_8)*phi)/rho)
    v1ca  = -sc/d*dda+v1c
    v1cb  = -sc/d*ddb+v1c
    v2c   = cn*ephi/br4*(2._real_8-phi)/d
    v2ca  = v2c
    v2cb  = v2c
    v2cab = v2c
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE lsd_p86
  ! ==================================================================
  SUBROUTINE lsd_glyp(ra,rb,grhoaa,grhoab,grhobb,sc,&
       v1ca,v2ca,v1cb,v2cb,v2cab)
    ! ==--------------------------------------------------------------==
    ! LEE, YANG PARR: GRADIENT CORRECTION PART
    REAL(real_8)                             :: ra, rb, grhoaa, grhoab, &
                                                grhobb, sc, v1ca, v2ca, v1cb, &
                                                v2cb, v2cab

    REAL(real_8), PARAMETER                  :: a = 0.04918_real_8, &
                                                b = 0.132_real_8, &
                                                c = 0.2533_real_8, &
                                                d = 0.349_real_8 

    REAL(real_8)                             :: dder, der, dlaa, dlaaa, &
                                                dlaab, dlab, dlaba, dlabb, &
                                                dlbb, dlbba, dlbbb, dor, dr, &
                                                or, rho, rm3

! ==--------------------------------------------------------------==

    rho=ra+rb
    rm3=rho**(-1._real_8/3._real_8)
    dr=(1._real_8+d*rm3)
    or=EXP(-c*rm3)/dr*rm3**11._real_8
    dor=-1._real_8/3._real_8*rm3**4*or*(11._real_8/rm3-c-d/dr)
    der=c*rm3+d*rm3/dr
    dder=1._real_8/3._real_8*(d*d*rm3**5/dr/dr-der/rho)
    dlaa=-a*b*or*(ra*rb/9._real_8*(1._real_8-3*der-(der-11._real_8)*ra/rho)-rb*rb)
    dlab=-a*b*or*(ra*rb/9._real_8*(47._real_8-7._real_8*der)-4._real_8/3._real_8*rho*rho)
    dlbb=-a*b*or*(ra*rb/9._real_8*(1._real_8-3*der-(der-11._real_8)*rb/rho)-ra*ra)
    dlaaa=dor/or*dlaa-a*b*or*(rb/9._real_8*(1._real_8-3*der-(der-11._real_8)*ra/rho)-&
         ra*rb/9._real_8*((3._real_8+ra/rho)*dder+(der-11._real_8)*rb/rho/rho))
    dlaab=dor/or*dlaa-a*b*or*(ra/9._real_8*(1._real_8-3._real_8*der&
         -(der-11._real_8)*ra/rho)-ra*rb/9._real_8*((3._real_8+ra/rho)*dder&
         -(der-11._real_8)*ra/rho/rho)-2._real_8*rb)
    dlaba=dor/or*dlab-a*b*or*(rb/9._real_8*(47._real_8-7._real_8*der)&
         -7._real_8/9._real_8*ra*rb*dder-8._real_8/3._real_8*rho)
    dlabb=dor/or*dlab-a*b*or*(ra/9._real_8*(47._real_8-7._real_8*der)&
         -7._real_8/9._real_8*ra*rb*dder-8._real_8/3._real_8*rho)
    dlbba=dor/or*dlbb-a*b*or*(rb/9._real_8*(1._real_8-3._real_8*der&
         -(der-11._real_8)*rb/rho)-ra*rb/9._real_8*((3._real_8+rb/rho)*dder&
         -(der-11._real_8)*rb/rho/rho)-2._real_8*ra)
    dlbbb=dor/or*dlbb-a*b*or*(ra/9._real_8*(1._real_8-3.0_real_8*der&
         -(der-11._real_8)*rb/rho)-ra*rb/9._real_8*((3._real_8+rb/rho)*dder&
         +(der-11._real_8)*ra/rho/rho))
    sc=dlaa*grhoaa+dlab*grhoab+dlbb*grhobb
    v1ca=dlaaa*grhoaa+dlaba*grhoab+dlbba*grhobb
    v1cb=dlaab*grhoaa+dlabb*grhoab+dlbbb*grhobb
    v2ca=2._real_8*dlaa
    v2cb=2._real_8*dlbb
    v2cab=dlab
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE lsd_glyp
  ! ==================================================================
  SUBROUTINE lsd_pbesc(rhoa,rhob,grhoaa,grhoab,grhobb,sc,&
       v1ca,v2ca,v1cb,v2cb,v2cab)
    ! PBE Correlation functional
    REAL(real_8)                             :: rhoa, rhob, grhoaa, grhoab, &
                                                grhobb, sc, v1ca, v2ca, v1cb, &
                                                v2cb, v2cab

    REAL(real_8), PARAMETER :: be = 0.046_real_8, &
      ga = 0.031090690869654895_real_8 , ob3 = 1._real_8/3._real_8 , &
      small = 1.e-24_real_8 

    REAL(real_8) :: a, aa, af, dadra, dadrb, dedra, dedrb, dhdra, dhdrb, &
      dhdt, dphida, dphidb, dphide, dsda, dsdra, dsdrb, dsdt, dtdphi, dtdra, &
      dtdrb, ec, eta, ex, expe, grho, h0, phi, phi3, rho, rs, s1, t, vca, &
      vcb, vxa, vxb, xkf, xks, xy, y

! ==--------------------------------------------------------------==

    rho=rhoa+rhob
    eta=(rhoa-rhob)/rho
    phi=0.5_real_8*((1._real_8+eta)**(2._real_8/3._real_8)+(1._real_8-eta)**(2._real_8/3._real_8))
    phi3=phi*phi*phi
    grho=grhoaa+2._real_8*grhoab+grhobb
    CALL xc_lsd(rho,eta,ex,ec,vxa,vca,vxb,vcb)
    IF (ABS(ec).LT.small) THEN
       sc = 0._real_8
       v1ca  = 0._real_8
       v2ca  = 0._real_8
       v1cb  = 0._real_8
       v2cb  = 0._real_8
       v2cab = 0._real_8
       RETURN
    ENDIF
    aa    = MAX(grho,small)
    a     = SQRT(aa)
    rs    = (3._real_8/(4._real_8*pi*rho))**ob3
    xkf   = (9._real_8*pi/4._real_8)**ob3/rs
    xks   = SQRT(4._real_8*xkf/pi)
    t     = a/(2._real_8*xks*rho*phi)
    expe  = EXP(-ec/(phi3*ga))
    af    = be/ga * (1._real_8/(expe-1._real_8))
    y     = af*t*t
    xy    = (1._real_8+y)/(1._real_8+y+y*y)
    s1    = 1._real_8+be/ga*t*t*xy
    h0    = ga*phi3 * LOG(s1)
    sc    = rho*h0
    ! ..
    IF (eta.LT.0.999999999999_real_8 .AND. eta.GT.-0.999999999999_real_8) THEN
       dtdphi= -t/phi
       dphide= 1._real_8/(3._real_8*(1._real_8+eta)**ob3)-1._real_8/(3._real_8*(1._real_8-eta)**ob3)
       dedra = 2._real_8*rhob/(rho*rho)
       dedrb = -2._real_8*rhoa/(rho*rho)
       dphida= dphide*dedra
       dphidb= dphide*dedrb
       dtdra = -t*(dphida/phi+7._real_8/(6._real_8*rho))
       dtdrb = -t*(dphidb/phi+7._real_8/(6._real_8*rho))
       dadra = af*af*expe/(-be*phi3)*(3._real_8*ec/phi*dphida-(vca-ec)/rho)
       dadrb = af*af*expe/(-be*phi3)*(3._real_8*ec/phi*dphidb-(vcb-ec)/rho)
       dsda  = -be/ga * af * t**6 * (2._real_8+y) / (1._real_8+y+y*y)**2
       dsdt  = 2._real_8*be/ga * t * (1._real_8+2._real_8*y) / (1._real_8+y+y*y)**2
       dsdra = dsda*dadra + dsdt*dtdra
       dsdrb = dsda*dadrb + dsdt*dtdrb
       dhdt  = ga*phi3/s1*dsdt
       dhdra = 3._real_8*h0/phi*dphida + ga*phi3/s1*dsdra
       dhdrb = 3._real_8*h0/phi*dphidb + ga*phi3/s1*dsdrb
    ELSE
       dtdra = -t*(7._real_8/(6._real_8*rho))
       dtdrb = -t*(7._real_8/(6._real_8*rho))
       dadra = af*af*expe/(-be*phi3)*(-(vca-ec)/rho)
       dadrb = af*af*expe/(-be*phi3)*(-(vcb-ec)/rho)
       dsda  = -be/ga * af * t**6 * (2._real_8+y) / (1._real_8+y+y*y)**2
       dsdt  = 2._real_8*be/ga * t * (1._real_8+2._real_8*y) / (1._real_8+y+y*y)**2
       dsdra = dsda*dadra + dsdt*dtdra
       dsdrb = dsda*dadrb + dsdt*dtdrb
       dhdt  = ga*phi3/s1*dsdt
       dhdra = ga*phi3/s1*dsdra
       dhdrb = ga*phi3/s1*dsdrb
    ENDIF
    ! ..
    v1ca  = h0 + rho*dhdra
    v2ca  = rho*dhdt*t/aa
    v1cb  = h0 + rho*dhdrb
    v2cb  = rho*dhdt*t/aa
    v2cab = rho*dhdt*t/aa
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE lsd_pbesc
  ! ==================================================================
  SUBROUTINE lsd_pbec(rhoa,rhob,grhoaa,grhoab,grhobb,sc,&
       v1ca,v2ca,v1cb,v2cb,v2cab)
    ! PBE Correlation functional
    REAL(real_8)                             :: rhoa, rhob, grhoaa, grhoab, &
                                                grhobb, sc, v1ca, v2ca, v1cb, &
                                                v2cb, v2cab

    REAL(real_8), PARAMETER :: be = 0.06672455060314922_real_8, &
      ga = 0.031090690869654895_real_8 , ob3 = 1._real_8/3._real_8 , &
      small = 1.e-24_real_8 

    REAL(real_8) :: a, aa, af, dadra, dadrb, dedra, dedrb, dhdra, dhdrb, &
      dhdt, dphida, dphidb, dphide, dsda, dsdra, dsdrb, dsdt, dtdphi, dtdra, &
      dtdrb, ec, eta, ex, expe, grho, h0, phi, phi3, rho, rs, s1, t, vca, &
      vcb, vxa, vxb, xkf, xks, xy, y

! ==--------------------------------------------------------------==

    rho=rhoa+rhob
    eta=(rhoa-rhob)/rho
    phi=0.5_real_8*((1._real_8+eta)**(2._real_8/3._real_8)+(1._real_8-eta)**(2._real_8/3._real_8))
    phi3=phi*phi*phi
    grho=grhoaa+2._real_8*grhoab+grhobb
    CALL xc_lsd(rho,eta,ex,ec,vxa,vca,vxb,vcb)
    IF (ABS(ec).LT.small) THEN
       sc = 0._real_8
       v1ca  = 0._real_8
       v2ca  = 0._real_8
       v1cb  = 0._real_8
       v2cb  = 0._real_8
       v2cab = 0._real_8
       RETURN
    ENDIF
    aa    = MAX(grho,small)
    a     = SQRT(aa)
    rs    = (3._real_8/(4._real_8*pi*rho))**ob3
    xkf   = (9._real_8*pi/4._real_8)**ob3/rs
    xks   = SQRT(4._real_8*xkf/pi)
    t     = a/(2._real_8*xks*rho*phi)
    expe  = EXP(-ec/(phi3*ga))
    af    = be/ga * (1._real_8/(expe-1._real_8))
    y     = af*t*t
    xy    = (1._real_8+y)/(1._real_8+y+y*y)
    s1    = 1._real_8+be/ga*t*t*xy
    h0    = ga*phi3 * LOG(s1)
    sc    = rho*h0
    ! ..
    IF (eta.LT.0.999999999999_real_8 .AND. eta.GT.-0.999999999999_real_8) THEN
       dtdphi= -t/phi
       dphide= 1._real_8/(3._real_8*(1._real_8+eta)**ob3)-1._real_8/(3._real_8*(1._real_8-eta)**ob3)
       dedra = 2._real_8*rhob/(rho*rho)
       dedrb = -2._real_8*rhoa/(rho*rho)
       dphida= dphide*dedra
       dphidb= dphide*dedrb
       dtdra = -t*(dphida/phi+7._real_8/(6._real_8*rho))
       dtdrb = -t*(dphidb/phi+7._real_8/(6._real_8*rho))
       dadra = af*af*expe/(-be*phi3)*(3._real_8*ec/phi*dphida-(vca-ec)/rho)
       dadrb = af*af*expe/(-be*phi3)*(3._real_8*ec/phi*dphidb-(vcb-ec)/rho)
       dsda  = -be/ga * af * t**6 * (2._real_8+y) / (1._real_8+y+y*y)**2
       dsdt  = 2._real_8*be/ga * t * (1._real_8+2._real_8*y) / (1._real_8+y+y*y)**2
       dsdra = dsda*dadra + dsdt*dtdra
       dsdrb = dsda*dadrb + dsdt*dtdrb
       dhdt  = ga*phi3/s1*dsdt
       dhdra = 3._real_8*h0/phi*dphida + ga*phi3/s1*dsdra
       dhdrb = 3._real_8*h0/phi*dphidb + ga*phi3/s1*dsdrb
    ELSE
       dtdra = -t*(7._real_8/(6._real_8*rho))
       dtdrb = -t*(7._real_8/(6._real_8*rho))
       dadra = af*af*expe/(-be*phi3)*(-(vca-ec)/rho)
       dadrb = af*af*expe/(-be*phi3)*(-(vcb-ec)/rho)
       dsda  = -be/ga * af * t**6 * (2._real_8+y) / (1._real_8+y+y*y)**2
       dsdt  = 2._real_8*be/ga * t * (1._real_8+2._real_8*y) / (1._real_8+y+y*y)**2
       dsdra = dsda*dadra + dsdt*dtdra
       dsdrb = dsda*dadrb + dsdt*dtdrb
       dhdt  = ga*phi3/s1*dsdt
       dhdra = ga*phi3/s1*dsdra
       dhdrb = ga*phi3/s1*dsdrb
    ENDIF
    ! ..
    v1ca  = h0 + rho*dhdra
    v2ca  = rho*dhdt*t/aa
    v1cb  = h0 + rho*dhdrb
    v2cb  = rho*dhdt*t/aa
    v2cab = rho*dhdt*t/aa
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE lsd_pbec
  ! ==================================================================
  SUBROUTINE lsd_hcth(rhoa,rhob,grhoa,grhob,sx,v1xa,v2xa,v1xb,v2xb)
    ! HCTH, JCP 109, 6264 (1998)
    ! Parameters set-up after N.L. Doltsisnis & M. Sprik (1999)
    ! Present release: Tsukuba, 09/02/2005
    ! -----------------------------------------------------------------------
    ! rho=rhoa+rhob
    ! grhoa is the SQUARE of the gradient of rhoa --> gra=sqrt(grhoa)
    ! grhob is the SQUARE of the gradient of rhob --> grb=sqrt(grhob)
    ! sx  : total exchange correlation energy at point r 
    ! v1xa: d(sx)/drhoa 
    ! v1xb: d(sx)/drhob
    ! v2xa: 1/gra*d(sx)/d(gra) 
    ! v2xb: 1/grb*d(sx)/d(grb) 
    ! Note: d(fab)/drhoa = d(fab)/drho * d(rho)/drhoa = d(fab)/drho
    ! and analogously for rhob. We make use of this in the present
    ! implementation whenever possible. V2CAB = 0
    ! --------------------------------------------------------------------------
    REAL(real_8)                             :: rhoa, rhob, grhoa, grhob, sx, &
                                                v1xa, v2xa, v1xb, v2xb

    REAL(real_8), PARAMETER :: f43 = 4._real_8/3._real_8 , &
      f83 = 8._real_8/3._real_8, o3 = 1.0_real_8/3.0_real_8, &
      o43 = 4.0_real_8*o3, pi34 = 3._real_8/(4._real_8*pi), &
      pi34_o3 = 0.620350490899400017_real_8, pi36 = 36._real_8*pi, &
      pi36_o3 = 0.206783496966466672_real_8, pi4 = 4._real_8*pi, &
      sixpi = 6._real_8/pi, sixpi_o3 = 1.240700981798800033_real_8  , &
      smalg = 1.e-20_real_8 , small = 1.e-15_real_8, &
      two_o3 = 1.259921049894873165_real_8, &
      yuu = 81.0_real_8/(32.0_real_8*pi), yuu_o3 = 0.930525736349100025_real_8

    REAL(real_8) :: alpha, caa(6), cab(6), cal(6), cg0(6), cg1(6), cx(6), &
      d2ffun0, dalpha_drab, dera1_dra, derab0_drab, derab1_drab, derab_drab, &
      derab_drhoa, derab_drhob, derab_dzeta, derb1_drb, dexa_drhoa, &
      dexb_drhob, dffaa_drhoa, dffab_drhoa, dffab_drhob, dffbb_drhob, dffun, &
      dg, dgaa_dgra, dgaa_drhoa, dgab_dgra, dgab_dgrb, dgab_drhoa, &
      dgab_drhob, dgbb_dgrb, dgbb_drhob, dgxa_dgra, dgxa_drhoa, dgxb_dgrb, &
      dgxb_drhob, dra_drhoa, drab_drho, drb_drhob, dzeta_drhoa, dzeta_drhob, &
      era1, erab, erab0, erab1, erb1, exa, exb, f2ga, f2gb, f4ra, f4rb, f8ra, &
      f8rb, f8xa, f8xb, fact, ffaa, ffab, ffbb, ffun
    REAL(real_8) :: fgaa, fgbb, fxaa, fxbb, g, gaa, gab, gbb, gra, grb, gxa, &
      gxb, qr1mz, qrzp1, ra, rab, rb, rho, rho_o3, rho_o43, rhoa_o3, &
      rhoa_o43, rhob_o3, rhob_o43, two13, uaa, uaai, uab, uabi, ubb, ubbi, &
      uxa, uxai, uxb, uxbi, xa, xa2, xab2, xb, xb2, zeta, zeta3, zeta4

! 

    two13=EXP(o3*0.69314718055994531_real_8)
    fact=0.5_real_8/(two13-1.0_real_8)
    sx=0.0e+00_real_8
    v1xa=0.0e+00_real_8
    v2xa=0.0e+00_real_8
    v1xb=0.0e+00_real_8
    v2xb=0.0e+00_real_8
    ! 
    ! .....coefficients for PW correlation...........................
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
    cal(1)= 0.016887_real_8
    cal(2)= 0.111250_real_8
    cal(3)= 10.35700_real_8
    cal(4)= 3.623100_real_8
    cal(5)= 0.880260_real_8
    cal(6)= 0.496710_real_8
    ! .....coefficients for HCTH-19-4/HCTH-120 after Nikos............
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
    ! -> spin up + down treated simultaneously
    ! mb-> put a threshold to fix numerical instability
    rhoa=MAX(rhoa,small)
    rhob=MAX(rhob,small)
    rho=rhoa+rhob
    gra=SQRT(grhoa)
    grb=SQRT(grhob)
    ! mb-> put a threshold to fix numerical instability
    gra=MAX(gra,smalg)
    grb=MAX(grb,smalg)
    rhoa_o3=rhoa**o3
    rhob_o3=rhob**o3
    rho_o3=rho**o3
    rhoa_o43=rhoa*rhoa_o3
    rhob_o43=rhob*rhob_o3
    rho_o43=rho*rho_o3
    xa=gra/rhoa_o43
    xb=grb/rhob_o43
    xa2=xa*xa
    xb2=xb*xb
    xab2=0.5_real_8*(xa2+xb2)
    ra=pi34_o3/rhoa_o3
    rb=pi34_o3/rhob_o3
    rab=pi34_o3/rho_o3
    dra_drhoa=-pi36_o3/rhoa_o43
    drb_drhob=-pi36_o3/rhob_o43
    drab_drho=-pi36_o3/rho_o43
    CALL pwcorr2(ra,cg1,g,dg)
    era1=g
    dera1_dra=dg
    CALL pwcorr2(rb,cg1,g,dg)
    erb1=g
    derb1_drb=dg
    CALL pwcorr2(rab,cg0,g,dg)
    erab0=g
    derab0_drab=dg
    CALL pwcorr2(rab,cg1,g,dg)
    erab1=g
    derab1_drab=dg
    CALL pwcorr2(rab,cal,g,dg)
    alpha=-g
    dalpha_drab=-dg
    zeta=(rhoa-rhob)/rho
    ! 
    qrzp1=two_o3*rhoa_o3/rho_o3
    qr1mz=two_o3*rhob_o3/rho_o3
    ffun=((1._real_8+zeta)*qrzp1+(1._real_8-zeta)*qr1mz-2._real_8)*fact
    dffun=o43*(qrzp1-qr1mz)*fact
    d2ffun0=8._real_8/9._real_8*fact
    zeta3=zeta*zeta*zeta
    zeta4=zeta3*zeta
    erab=erab0+alpha*ffun*(1._real_8-zeta4)/d2ffun0
    erab=erab+(erab1-erab0)*ffun*zeta4
    derab_drab=derab0_drab+dalpha_drab*ffun*(1._real_8-zeta4)/d2ffun0&
         +(derab1_drab-derab0_drab)*ffun*zeta4
    derab_dzeta=alpha/d2ffun0*(dffun*(1._real_8-zeta4)-4._real_8*zeta3*ffun)&
         +(erab1-erab0)*(dffun*zeta4+4._real_8*zeta3*ffun)
    ! 
    exa=-yuu_o3*rhoa_o43
    exb=-yuu_o3*rhob_o43
    dexa_drhoa=-sixpi_o3*rhoa_o3
    dexb_drhob=-sixpi_o3*rhob_o3
    dzeta_drhoa=2._real_8*rhob/(rho*rho)
    dzeta_drhob=-2._real_8*rhoa/(rho*rho)
    ! 
    uaa=caa(6)*xa2
    uaa=uaa/(1.0_real_8+uaa)
    ubb=caa(6)*xb2
    ubb=ubb/(1.0_real_8+ubb)
    uab=cab(6)*xab2
    uab=uab/(1.0_real_8+uab)
    uxa=cx(6)*xa2
    uxa=uxa/(1.0_real_8+uxa)
    uxb=cx(6)*xb2
    uxb=uxb/(1.0_real_8+uxb)
    ffaa=rhoa*era1
    ffbb=rhob*erb1
    ffab=rho*erab0-ffaa-ffbb
    derab_drhoa=derab_drab*drab_drho+derab_dzeta*dzeta_drhoa
    derab_drhob=derab_drab*drab_drho+derab_dzeta*dzeta_drhob
    dffaa_drhoa=era1+rhoa*dera1_dra*dra_drhoa
    dffbb_drhob=erb1+rhob*derb1_drb*drb_drhob
    dffab_drhoa=erab+rho*derab_drhoa-dffaa_drhoa
    dffab_drhob=erab+rho*derab_drhob-dffbb_drhob
    ! mb-> i-loop removed
    f8ra=(f83/rhoa)/(1.0_real_8+caa(6)*xa2)
    f8rb=(f83/rhob)/(1.0_real_8+caa(6)*xb2)
    f4ra=(f43/rhoa)*(xa2/xab2)/(1.0_real_8+cab(6)*xab2)
    f4rb=(f43/rhob)*(xb2/xab2)/(1.0_real_8+cab(6)*xab2)
    f8xa=(f83/rhoa)/(1.0_real_8+cx(6)*xa2)
    f8xb=(f83/rhob)/(1.0_real_8+cx(6)*xb2)
    f2ga=2.0_real_8/gra/(1.0_real_8+caa(6)*xa2)
    f2gb=2.0_real_8/grb/(1.0_real_8+caa(6)*xb2)
    fgaa=(xa2/xab2)/gra/(1.0_real_8+cab(6)*xab2)
    fgbb=(xb2/xab2)/grb/(1.0_real_8+cab(6)*xab2)
    fxaa=2.0_real_8/gra/(1.0_real_8+cx(6)*xa2)
    fxbb=2.0_real_8/grb/(1.0_real_8+cx(6)*xb2)
    ! 
    gaa=caa(1)+uaa*(caa(2)+uaa*(caa(3)+uaa*(caa(4)+uaa*caa(5))))
    gbb=caa(1)+ubb*(caa(2)+ubb*(caa(3)+ubb*(caa(4)+ubb*caa(5))))
    gab=cab(1)+uab*(cab(2)+uab*(cab(3)+uab*(cab(4)+uab*cab(5))))
    gxa=cx(1)+uxa*(cx(2)+uxa*(cx(3)+uxa*(cx(4)+uxa*cx(5))))
    gxb=cx(1)+uxb*(cx(2)+uxb*(cx(3)+uxb*(cx(4)+uxb*cx(5))))
    uaai=uaa*(caa(2)+uaa*(2.0_real_8*caa(3)+uaa*(3.0_real_8*caa(4)&
         +uaa*4.0_real_8*caa(5))))
    ubbi=ubb*(caa(2)+ubb*(2.0_real_8*caa(3)+ubb*(3.0_real_8*caa(4)&
         +ubb*4.0_real_8*caa(5))))
    uabi=uab*(cab(2)+uab*(2.0_real_8*cab(3)+uab*(3.0_real_8*cab(4)&
         +uab*4.0_real_8*cab(5))))
    uxai=uxa*(cx(2)+uxa*(2.0_real_8*cx(3)+uxa*(3.0_real_8*cx(4)&
         +uxa*4.0_real_8*cx(5))))
    uxbi=uxb*(cx(2)+uxb*(2.0_real_8*cx(3)+uxb*(3.0_real_8*cx(4)&
         +uxb*4.0_real_8*cx(5))))
    dgaa_drhoa=-f8ra*uaai
    dgbb_drhob=-f8rb*ubbi
    dgab_drhoa=-f4ra*uabi
    dgab_drhob=-f4rb*uabi
    dgxa_drhoa=-f8xa*uxai
    dgxb_drhob=-f8xb*uxbi
    dgaa_dgra=f2ga*uaai
    dgbb_dgrb=f2gb*ubbi
    dgab_dgra=fgaa*uabi
    dgab_dgrb=fgbb*uabi
    dgxa_dgra=fxaa*uxai
    dgxb_dgrb=fxbb*uxbi
    ! 
    sx=exa*gxa+exb*gxb+ffaa*gaa+ffbb*gbb+ffab*gab
    v1xa=dexa_drhoa*gxa+exa*dgxa_drhoa&
         +dffaa_drhoa*gaa+ffaa*dgaa_drhoa&
         +dffab_drhoa*gab+ffab*dgab_drhoa
    v2xa=(exa*dgxa_dgra+ffaa*dgaa_dgra+ffab*dgab_dgra)/gra
    v1xb=dexb_drhob*gxb+exb*dgxb_drhob&
         +dffbb_drhob*gbb+ffbb*dgbb_drhob&
         +dffab_drhob*gab+ffab*dgab_drhob
    v2xb=(exb*dgxb_dgrb+ffbb*dgbb_dgrb+ffab*dgab_dgrb)/grb
    ! 
    RETURN
  END SUBROUTINE lsd_hcth
  ! =-------------------------------------------------------------------=
  SUBROUTINE pwcorr2(r,c,g,dg)
    REAL(real_8)                             :: r, c(6), g, dg

    REAL(real_8)                             :: drb, r12, r2, r32, rb, sb

    r12=SQRT(r)
    r32=r*r12
    r2=r*r
    rb=c(3)*r12+c(4)*r+c(5)*r32+c(6)*r2
    sb=1.0_real_8+1.0_real_8/(2.0_real_8*c(1)*rb)
    g=-2.0_real_8*c(1)*(1.0_real_8+c(2)*r)*LOG(sb)
    drb=c(3)/(2.0_real_8*r12)+c(4)+1.5_real_8*c(5)*r12+2.0_real_8*c(6)*r
    dg=(1.0_real_8+c(2)*r)*drb/(rb*rb*sb)&
         -2.0_real_8*c(1)*c(2)*LOG(sb)
    RETURN
  END SUBROUTINE pwcorr2
  ! ==================================================================
  SUBROUTINE lsd_optx(rhoa,rhob,grhoa,grhob&
       ,sx,v1xa,v2xa,v1xb,v2xb)
    ! OPTX, Handy et al. JCP 116, p. 5411 (2002) and refs. therein
    ! Present release: Tsukuba, 21/6/2002
    ! --------------------------------------------------------------------------
    ! grho_a,b is the SQUARE of the gradient of rho_a,b ! 
    ! sx  : total exchange correlation energy at point r
    ! v1x_a,b : d(sx)/drho_a,b
    ! v2x_a,b : 1/gr_a,b * d(sx)/d(gr_a,b)
    ! --------------------------------------------------------------------------
    REAL(real_8)                             :: rhoa, rhob, grhoa, grhob, sx, &
                                                v1xa, v2xa, v1xb, v2xb

    REAL(real_8), PARAMETER :: a1cx = 0.9784571170284421_real_8, &
      a2 = 1.43169_real_8 , gam = 0.006_real_8, o43 = 4.0_real_8/3.0_real_8 , &
      smal2 = 1.e-08_real_8 , small = 1.e-20_real_8

    REAL(real_8)                             :: exa, exb, gamxa2, gamxb2, &
                                                gra, grb, rhoa43, rhob43, &
                                                udena, udenb, uua, uub, xa, xb

    IF (rhoa.LE.small) THEN
       exa=0.0_real_8
       v1xa=0.0_real_8
       v2xa=0.0_real_8
    ELSE
       gra=MAX(grhoa,smal2)
       rhoa43=rhoa**o43
       xa=SQRT(gra)/rhoa43
       gamxa2=gam*xa*xa
       udena=1._real_8/(1._real_8+gamxa2)
       uua=a2*gamxa2*gamxa2*udena*udena
       udena=4.00_real_8*rhoa43*uua*udena
       exa=-rhoa43*(a1cx+uua)
       v1xa=o43*(exa+udena)/rhoa
       v2xa=-udena/gra
    ENDIF
    ! --> BETA
    IF (rhob.LE.small) THEN
       exb=0.0_real_8
       v1xb=0.0_real_8
       v2xb=0.0_real_8
    ELSE
       grb=MAX(grhob,smal2)
       rhob43=rhob**o43
       xb=SQRT(grb)/rhob43
       gamxb2=gam*xb*xb
       udenb=1._real_8/(1._real_8+gamxb2)
       uub=a2*gamxb2*gamxb2*udenb*udenb
       udenb=4.00_real_8*rhob43*uub*udenb
       exb=-rhob43*(a1cx+uub)
       v1xb=o43*(exb+udenb)/rhob
       v2xb=-udenb/grb
    ENDIF
    sx=exa+exb
    ! 
    RETURN
  END SUBROUTINE lsd_optx
  ! =-------------------------------------------------------------------=

END MODULE lsd_func_utils
