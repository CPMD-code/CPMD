#ifdef __SR8000
!option MP(P(0)), LANGLVL(SAVE(0))
#endif


MODULE metafun_utils
  USE cnst,                            ONLY: pi
  USE error_handling,                  ONLY: stopgm
  USE func,                            ONLY: func1,&
                                             mfxcc_is_skipped,&
                                             mfxcx_is_skipped,&
                                             mgcc_is_skipped,&
                                             mgcx_is_skipped,&
                                             mtau_is_tpss
  USE functionals_utils,               ONLY: pz
  USE kinds,                           ONLY: real_8
  USE lsd_func_utils,                  ONLY: lsd_pz

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: taufun
  PUBLIC :: taufuns
  PUBLIC :: tpss
  PUBLIC :: tpss_lsd
  PUBLIC :: tpssx
  PUBLIC :: tpssc
  PUBLIC :: tpssc_lsd
  PUBLIC :: revpkzb
  PUBLIC :: tpbec
  PUBLIC :: tpbeca
  PUBLIC :: cgefun
  PUBLIC :: ccfun
  PUBLIC :: revpkzb_lsd
  PUBLIC :: tpbecs

CONTAINS

  ! ==================================================================
  SUBROUTINE taufun(rho,grho,tau,sx,sc,v1x,v2x,v1c,v2c,vtt)
    ! ==--------------------------------------------------------------==
    ! ==  META FUNCTIONALS                                            ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rho, grho, tau, sx, sc, v1x, &
                                                v2x, v1c, v2c, vtt

    REAL(real_8), PARAMETER                  :: small = 1.e-20_real_8 

! ==--------------------------------------------------------------==

    IF ( func1%mfxcx /= mfxcx_is_skipped .OR. func1%mfxcc /= mfxcc_is_skipped .OR.&
         func1%mgcx /= mgcx_is_skipped .OR. func1%mgcc /= mgcc_is_skipped )&
         CALL stopgm("TAUFUN","FUNCTIONALS DO NOT MATCH",& 
         __LINE__,__FILE__)
    IF (func1%mtau == mtau_is_tpss .AND. rho.GT.small) THEN
       CALL tpss(rho,grho,tau,sx,sc,v1x,v2x,v1c,v2c,vtt)
    ELSE
       sx=0._real_8
       sc=0._real_8
       v1x=0._real_8
       v2x=0._real_8
       v1c=0._real_8
       v2c=0._real_8
       vtt=0._real_8
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE taufun
  ! ==================================================================
  SUBROUTINE taufuns(rhoa,rhob,grhoaa,grhobb,grhoab,taua,taub,&
       sx,sc,v1xa,v2xa,v1xb,v2xb,v1ca,v2ca,&
       v1cb,v2cb,v2xab,v2cab,vtta,vttb)
    ! ==--------------------------------------------------------------==
    ! ==  META FUNCTIONALS                                            ==
    ! ==--------------------------------------------------------------==
    REAL(real_8) :: rhoa, rhob, grhoaa, grhobb, grhoab, taua, taub, sx, sc, &
      v1xa, v2xa, v1xb, v2xb, v1ca, v2ca, v1cb, v2cb, v2xab, v2cab, vtta, vttb

    REAL(real_8), PARAMETER                  :: small = 1.e-20_real_8 

! ==--------------------------------------------------------------==

    IF ( func1%mfxcx /= mfxcx_is_skipped .OR. func1%mfxcc /= mfxcc_is_skipped .OR.&
         func1%mgcx /= mgcx_is_skipped .OR. func1%mgcc /= mgcc_is_skipped )&
         CALL stopgm("TAUFUN","FUNCTIONALS DO NOT MATCH",& 
         __LINE__,__FILE__)
    IF (func1%mtau == mtau_is_tpss .AND. (rhoa.GT.small.OR.rhob.GT.small)) THEN
       CALL tpss_lsd(rhoa,rhob,grhoaa,grhobb,grhoab,taua,taub,&
            sx,sc,v1xa,v2xa,v1xb,v2xb,v1ca,v2ca,&
            v1cb,v2cb,v2xab,v2cab,vtta,vttb)
    ELSE
       sx=0._real_8
       sc=0._real_8
       v1xa=0._real_8
       v1xb=0._real_8
       v2xa=0._real_8
       v2xb=0._real_8
       v1ca=0._real_8
       v1cb=0._real_8
       v2ca=0._real_8
       v2cb=0._real_8
       v2xab=0._real_8
       v2cab=0._real_8
       vtta=0._real_8
       vttb=0._real_8
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE taufuns
  ! ==================================================================
  SUBROUTINE tpss(rho,grho,tau,sx,sc,v1x,v2x,v1c,v2c,vtt)
    ! ==--------------------------------------------------------------==
    ! ==  TPSS FUNCTIONAL                                             ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rho, grho, tau, sx, sc, v1x, &
                                                v2x, v1c, v2c, vtt

! ==--------------------------------------------------------------==
! Exchange part

    CALL tpssx(rho,grho,tau,sx,v1x,v2x,vtt)
    ! Correlation
    CALL tpssc(rho,grho,tau,sc,v1c,v2c,vtt)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE tpss
  ! ==================================================================
  SUBROUTINE tpss_lsd(rhoa,rhob,grhoaa,grhobb,grhoab,taua,taub,&
       sx,sc,v1xa,v2xa,v1xb,v2xb,v1ca,v2ca,&
       v1cb,v2cb,v2xab,v2cab,vtta,vttb)
    ! ==--------------------------------------------------------------==
    ! ==  TPSS FUNCTIONAL  LSD Version                                ==
    ! ==--------------------------------------------------------------==
    REAL(real_8) :: rhoa, rhob, grhoaa, grhobb, grhoab, taua, taub, sx, sc, &
      v1xa, v2xa, v1xb, v2xb, v1ca, v2ca, v1cb, v2cb, v2xab, v2cab, vtta, vttb

    REAL(real_8), PARAMETER                  :: small = 1.e-20_real_8 

    REAL(real_8)                             :: grho, rho, sxx, tau, v1x, &
                                                v2x, vtt

! ==--------------------------------------------------------------==
! Exchange part
! alpha spin

    rho=2._real_8*rhoa
    grho=4._real_8*grhoaa
    tau=2._real_8*taua
    IF (rho.GT.small) THEN
       CALL tpssx(rho,grho,tau,sxx,v1x,v2x,vtt)
       sx=0.5_real_8*sxx
       v1xa=v1x
       v2xa=2._real_8*v2x
       vtta=vtt
    ELSE
       sx=0._real_8
       v1xa=0._real_8
       v2xa=0._real_8
       vtta=0._real_8
    ENDIF
    ! beta spin
    rho=2._real_8*rhob
    grho=4._real_8*grhobb
    tau=2._real_8*taub
    IF (rho.GT.small) THEN
       CALL tpssx(rho,grho,tau,sxx,v1x,v2x,vtt)
       sx=sx + 0.5_real_8*sxx
       v1xb=v1x
       v2xb=2._real_8*v2x
       vttb=vtt
    ELSE
       v1xb=0._real_8
       v2xb=0._real_8
       vttb=0._real_8
    ENDIF
    v2xab=0._real_8
    ! Correlation
    CALL tpssc_lsd(rhoa,rhob,grhoaa,grhobb,grhoab,taua,taub,&
         sc,v1ca,v2ca,v1cb,v2cb,v2cab,vtta,vttb)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE tpss_lsd
  ! ==================================================================
  SUBROUTINE tpssx(rho,grho,tau,sx,v1x,v2x,vtt)
    ! ==--------------------------------------------------------------==
    ! ==  EXCHANGE PART OF TPSS FUNCTIONAL                            ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rho, grho, tau, sx, v1x, v2x, &
                                                vtt

    REAL(real_8), PARAMETER :: b = 0.40_real_8 , c = 1.59096_real_8, &
      e = 1.537_real_8, kappa = 0.804_real_8, mu = 0.21951_real_8 , &
      small = 1.e-14_real_8 

    REAL(real_8) :: alpha, dfxdg, dfxdr, dfxdt, dfxdx, dpdg, dpdr, dqbtda, &
      dqbtdp, dqbtdz, dxdp, dxdqbt, dxdz, dzdg, dzdr, dzdt, exunif, f1081, &
      f13, f23, f83, fx, p, qbt, se, x, z

! ==--------------------------------------------------------------==

    f13 = 1._real_8/3._real_8
    f23 = 2._real_8*f13
    f83 = 8._real_8*f13
    f1081 = 10._real_8/81._real_8
    se = SQRT(e)
    ! PI  = ACOS(-1._real_8)
    ! 
    exunif = -3._real_8/(4._real_8*pi) * (3._real_8*pi*pi*rho)**f13
    IF ( ABS(tau) .GT. small .AND. grho .GT. small ) THEN
       p=grho/(4._real_8*(3._real_8*pi*pi)**f23*rho**f83)
       z=0.125_real_8*grho/rho/tau
       alpha = 5._real_8*f13*p * (1._real_8/z - 1._real_8)
       qbt = 0.45_real_8*(alpha-1._real_8)/SQRT(1._real_8+b*alpha*(alpha-1._real_8))+f23*p
       x = ( f1081 + c*z*z/(1._real_8+z*z)**2 ) * p + 146._real_8/2025._real_8*qbt**2&
            - 73._real_8/405._real_8 * qbt * SQRT ( 0.5_real_8*(0.6_real_8*z)**2+0.5_real_8*p*p)&
            + f1081**2*p*p/kappa+2._real_8*se*f1081*(0.6_real_8*z)**2+e*mu*p**3
       x = x/(1._real_8+se*p)**2
       fx = 1._real_8 + kappa - kappa/(1._real_8 + x/kappa)
    ELSE
       fx = 1._real_8
    ENDIF
    sx = rho*exunif*fx
    ! 
    IF ( ABS(tau) .GT. small .AND. grho .GT. small ) THEN
       dfxdx  = 1._real_8/(1._real_8 + x/kappa)**2
       dqbtda = 0.225_real_8*(2._real_8 + b*alpha - b)/&
            (1._real_8 + b*alpha*alpha - b*alpha)**1.5_real_8
       dqbtdp = dqbtda*5._real_8*f13*(1._real_8/z-1._real_8) + f23
       dqbtdz = dqbtda*5._real_8*f13*p*(-1._real_8/(z*z))
       dxdqbt = (146._real_8/2025._real_8 * 2._real_8*qbt -&
            73._real_8/405._real_8 * SQRT ( 0.5_real_8*(0.6_real_8*z)**2+0.5_real_8*p*p ))/&
            (1._real_8+se*p)**2
       dxdz   = dxdqbt*dqbtdz + ( c*p*(2._real_8*z/(1._real_8+z*z)**2 -&
            4._real_8*z**3/(1._real_8+z*z)**3 ) - 73._real_8/405._real_8 * qbt*0.5_real_8 *&
            0.36_real_8*z / SQRT ( 0.5_real_8*(0.6_real_8*z)**2 + 0.5_real_8*p*p ) +&
            4._real_8*se*f1081*0.36_real_8*z ) / (1._real_8+se*p)**2
       dxdp   = -2._real_8*se*x/(1._real_8+se*p) + dxdqbt*dqbtdp +&
            ( ( f1081 + c*z*z/(1._real_8+z*z)**2 ) - 73._real_8/405._real_8*qbt *&
            0.5_real_8*p/SQRT ( 0.5_real_8*(0.6_real_8*z)**2 + 0.5_real_8*p*p ) +&
            2._real_8*f1081**2*p/kappa + 3._real_8*e*mu*p**2 )/(1._real_8+se*p)**2
       dpdr   = -f83*p/rho
       dpdg   = 2._real_8*p/SQRT(grho)
       dzdr   = -z/rho
       dzdg   = 2._real_8*z/SQRT(grho)
       dzdt   = -z/tau
       dfxdr  = dfxdx*(dxdp*dpdr + dxdz*dzdr)
       dfxdg  = dfxdx*(dxdp*dpdg + dxdz*dzdg)
       dfxdt  = dfxdx*dxdz*dzdt
       ! 
       v1x = exunif*(fx + f13*fx + rho*dfxdr)
       v2x = rho*exunif*dfxdg/SQRT(grho)
       vtt = rho*exunif*dfxdt
    ELSE
       v1x = exunif + f13*sx/rho
       v2x = 0._real_8
       vtt = 0._real_8
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE tpssx
  ! ==================================================================
  SUBROUTINE tpssc(rho,grho,tau,sc,v1c,v2c,vtt)
    ! ==--------------------------------------------------------------==
    ! ==  CORRELATION PART OF TPSS FUNCTIONAL                         ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rho, grho, tau, sc, v1c, v2c, &
                                                vtt

    REAL(real_8), PARAMETER                  :: d = 2.8_real_8 
    INTEGER                                  :: iflg
    REAL(real_8)                             :: dedg, dedr, dedt, dedz, dzdg, &
                                                dzdr, dzdt, e, edgt, edrt, &
                                                op, ro, rs, z
    REAL(real_8), PARAMETER                  :: small = 1.e-14_real_8 

! ==--------------------------------------------------------------==

    IF ( ABS(tau) .GT. small .AND. grho .GT. small ) THEN
       z    = 0.125_real_8*grho/rho/tau
       dzdr = -z/rho
       dzdg = 2._real_8*z/SQRT(grho)
       dzdt = -z/tau
       CALL revpkzb(rho,grho,z,e,dedr,dedg,dedz)
       dedt = dedz*dzdt
       edrt = dedr + dedz*dzdr
       edgt = dedg + dedz*dzdg
       op  = 1._real_8 + d * e * z**3
       ro  = rho * e * d * z**2
       sc  = rho * e * op
       v1c = (e + rho*edrt)*op + ro * (edrt*z+3._real_8*e*dzdr)
       v2c = rho*edgt*op + ro * (edgt*z+3._real_8*e*dzdg)
       v2c = v2c/SQRT(grho)
       vtt = vtt  + rho*dedt*op + ro * (dedt*z+3._real_8*e*dzdt)
    ELSE
       ! PI   = ACOS(-1._real_8)
       rs   = (3._real_8/(4._real_8*pi*rho))**(1._real_8/3._real_8)
       iflg = 2
       IF (rs.LT.1.0_real_8) iflg=1
       CALL pz(rs,e,v1c,iflg)
       sc   = rho * e
       v2c  = 0._real_8
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE tpssc
  ! ==================================================================
  SUBROUTINE tpssc_lsd(rhoa,rhob,grhoaa,grhobb,grhoab,taua,taub,&
       sc,v1ca,v2ca,v1cb,v2cb,v2cab,vtta,vttb)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rhoa, rhob, grhoaa, grhobb, &
                                                grhoab, taua, taub, sc, v1ca, &
                                                v2ca, v1cb, v2cb, v2cab, &
                                                vtta, vttb

    REAL(real_8), PARAMETER                  :: d = 2.8_real_8 , &
                                                ob3 = 1._real_8/3._real_8 
    INTEGER                                  :: iflg
    REAL(real_8) :: dedga, dedgab, dedgb, dedra, dedrb, dedz, dzdg, dzdr, e, &
      ec, edaz, edbz, eta, grho, op, rdzdt, rho, ro, rs, tau, to, vca, vcb, z
    REAL(real_8), PARAMETER                  :: small = 1.e-14_real_8 

! ==--------------------------------------------------------------==

    rho = rhoa + rhob
    grho = grhoaa + 2._real_8*grhoab + grhobb
    tau = taua + taub
    IF ( ABS(tau) .GT. small .AND. grho .GT. small) THEN
       z    = 0.125_real_8*grho/rho/tau
       dzdr = -z/rho
       rdzdt = -z*rho/tau
       dzdg = z/grho
       CALL revpkzb_lsd(rhoa,rhob,grhoaa,grhobb,grhoab,z,e,&
            dedra,dedrb,dedga,dedgb,dedgab,dedz)
       op   = 1._real_8 + d * e * z**3
       sc   = rho * e * op
       edaz = dedra + dedz*dzdr
       v1ca = e*op + rho*edaz*op + rho*e*d*z*z*&
            (edaz*z+3._real_8*e*dzdr)
       edbz = dedrb + dedz*dzdr
       v1cb = e*op + rho*edbz*op + rho*e*d*z*z*&
            (edbz*z+3._real_8*e*dzdr)
       ro   = e*d*e*3._real_8*z*z
       to   = rho * ( 1._real_8 + 2._real_8*d*e*z**3 )
       v2ca = 2._real_8*rho*ro*dzdg + 2._real_8*to*dedz*dzdg + to*dedga
       v2cb = 2._real_8*rho*ro*dzdg + 2._real_8*to*dedz*dzdg + to*dedgb
       v2cab= 2._real_8*rho*ro*dzdg + 2._real_8*to*dedz*dzdg + to*dedgab
       vtta = vtta + dedz*rdzdt*op +&
            e*d*dedz*rdzdt*z**3 + ro*rdzdt
       vttb = vttb + dedz*rdzdt*op +&
            e*d*dedz*rdzdt*z**3 + ro*rdzdt
    ELSE
       eta=(rhoa-rhob)/rho
       IF (ABS(eta).GT.1._real_8) eta=SIGN(1.0_real_8,eta)
       rs = (0.75_real_8/(pi*rho))**ob3
       iflg=2
       IF (rs.LT.1.0_real_8) iflg=1
       CALL lsd_pz(rs,eta,ec,vca,vcb,iflg)
       sc   = rho * ec
       v1ca = vca
       v2ca = 0._real_8
       v1cb = vcb
       v2cb = 0._real_8
       v2cab= 0._real_8
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE tpssc_lsd
  ! ==================================================================
  SUBROUTINE revpkzb(rho,grho,z,e,dedr,dedg,dedz)
    ! ==--------------------------------------------------------------==
    ! ==  REVISED PKZB TAU CORRELATION FUNCTIONAL                     ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rho, grho, z, e, dedr, dedg, &
                                                dedz

    REAL(real_8), PARAMETER                  :: c00 = 0.53_real_8, &
                                                c0p = 1.53_real_8

    REAL(real_8)                             :: ecpbe, ecpbea, t1, tecpbe, &
                                                tv1c, tv2c, v1c, v1ca, v2c, &
                                                v2ca, z2

! ==--------------------------------------------------------------==

    CALL tpbec(rho,grho,ecpbe,v1c,v2c)
    CALL tpbeca(0.5_real_8*rho,0.25_real_8*grho,ecpbea,v1ca,v2ca)
    IF ( ecpbea .GT. ecpbe ) THEN
       tecpbe = ecpbea
       tv1c   = 0.5_real_8*v1ca
       tv2c   = 0.5_real_8*v2ca
    ELSE
       tecpbe = ecpbe
       tv1c   = v1c
       tv2c   = v2c
    ENDIF
    z2 = z*z
    t1 = 1._real_8 + c00*z2
    e    = ecpbe*t1 - c0p*z2*tecpbe
    dedr = v1c*t1 - c0p*z2*tv1c
    dedg = v2c*t1 - c0p*z2*tv2c
    dedz = 2._real_8*z * (ecpbe*c00 - c0p*tecpbe)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE revpkzb
  ! ==================================================================
  SUBROUTINE tpbec(rho,grho,sc,v1c,v2c)
    ! ==--------------------------------------------------------------==
    ! PBE Correlation functional
    REAL(real_8)                             :: rho, grho, sc, v1c, v2c

    REAL(real_8), PARAMETER :: be = 0.06672455060314922_real_8, &
      cx = -0.001667_real_8, cxc0 = 0.002568_real_8, cc0 = -cx+cxc0 , &
      ga = 0.031090690869654895_real_8 , ob3 = 1._real_8/3._real_8 

    INTEGER                                  :: iflg
    REAL(real_8)                             :: a, aa, af, dadr, decdr, dhdr, &
                                                dhdt, dsda, dsdr, dsdt, dtdr, &
                                                ec, expe, h0, rs, s1, t, vc, &
                                                xkf, xks, xy, y

! ==--------------------------------------------------------------==

    rs    = (3._real_8/(4._real_8*pi*rho))**ob3
    iflg  = 2
    IF (rs.LT.1.0_real_8) iflg=1
    CALL pz(rs,ec,vc,iflg)
    decdr = (vc-ec)/rho
    aa    = grho
    a     = SQRT(aa)
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
    dadr  = af*af*expe/be*decdr
    dsda  = -be/ga * af * t**6 * (2._real_8+y) / (1._real_8+y+y*y)**2
    dsdt  = 2._real_8*be/ga * t * (1._real_8+2._real_8*y) / (1._real_8+y+y*y)**2
    dsdr  = dsda*dadr + dsdt*dtdr
    dhdt  = ga/s1*dsdt
    dhdr  = ga/s1*dsdr
    sc    = ec + h0
    v1c   = decdr + dhdr
    v2c   = dhdt*t/a
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE tpbec
  ! ==================================================================
  SUBROUTINE tpbeca(rhoa,grhoaa,sc,v1ca,v2ca)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rhoa, grhoaa, sc, v1ca, v2ca

    REAL(real_8), PARAMETER :: be = 0.06672455060314922_real_8, &
      ga = 0.031090690869654895_real_8 , ob3 = 1._real_8/3._real_8 , &
      small = 1.e-24_real_8 

    INTEGER                                  :: iflg
    REAL(real_8) :: a, aa, af, dadra, decdra, dhdra, dhdt, dsda, dsdra, dsdt, &
      dtdra, ec, expe, grho, h0, phi, phi3, rho, rs, s1, t, vca, vcb, xkf, &
      xks, xy, y

! ==--------------------------------------------------------------==

    rho   = rhoa
    IF ( rho .GT. small ) THEN
       phi   = 0.5_real_8*2._real_8**(2._real_8/3._real_8)
       phi3  = phi*phi*phi
       grho  = grhoaa
       rs    = (3._real_8/(4._real_8*pi*rho))**ob3
       iflg=2
       IF (rs.LT.1.0_real_8) iflg=1
       CALL lsd_pz(rs,1._real_8,ec,vca,vcb,iflg)
       ! 
       aa    = MAX(grho,small)
       a     = SQRT(aa)
       xkf   = (9._real_8*pi/4._real_8)**ob3/rs
       xks   = SQRT(4._real_8*xkf/pi)
       t     = a/(2._real_8*xks*rho*phi)
       expe  = EXP(-ec/(phi3*ga))
       af    = be/ga * (1._real_8/(expe-1._real_8))
       y     = af*t*t
       xy    = (1._real_8+y)/(1._real_8+y+y*y)
       s1    = 1._real_8+be/ga*t*t*xy
       h0    = ga*phi3 * LOG(s1)
       ! 
       decdra= (vca-ec)/rho
       dtdra = -t*(7._real_8/(6._real_8*rho))
       dadra = af*af*expe/(-be*phi3)*(-decdra)
       dsda  = -be/ga * af * t**6 * (2._real_8+y) / (1._real_8+y+y*y)**2
       dsdt  = 2._real_8*be/ga * t * (1._real_8+2._real_8*y) / (1._real_8+y+y*y)**2
       dsdra = dsda*dadra + dsdt*dtdra
       dhdt  = ga*phi3/s1*dsdt
       dhdra = ga*phi3/s1*dsdra
       ! 
       sc    = ec + h0
       v1ca  = decdra + dhdra
       v2ca  = dhdt*t/a
    ELSE
       sc    = 0._real_8
       v1ca  = 0._real_8
       v2ca  = 0._real_8
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE tpbeca
  ! ==================================================================
  SUBROUTINE cgefun(xsi,eta2,cxe,dcdx,dcde)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: xsi, eta2, cxe, dcdx, dcde

    REAL(real_8), PARAMETER                  :: small = 1.e-14_real_8 

    REAL(real_8)                             :: cx0, dno, dx0, dxx, f43, f73, &
                                                x2, xx

! ==--------------------------------------------------------------==

    IF (ABS(xsi).LT.small.AND.ABS(eta2).LT.small*small) THEN
       cxe  = 0.53_real_8
       dcdx = 0._real_8
       dcde = 0._real_8
    ELSEIF(ABS(ABS(xsi)-1._real_8).LT.small.AND.&
         ABS(eta2).LT.small*small) THEN
       x2  = xsi*xsi
       cx0 = 0.53_real_8 + 0.87_real_8*x2 + 0.50_real_8 * x2*x2 + 2.26_real_8 * x2*x2*x2
       dx0 = 2._real_8*xsi * ( 0.87_real_8 + x2 + 3._real_8**2.26_real_8 * x2*x2 )
       cxe  = cx0
       dcdx = dx0
       dcde = 0._real_8
    ELSE
       f43 = -4._real_8/3._real_8
       f73 = -7._real_8/3._real_8
       x2  = xsi*xsi
       cx0 = 0.53_real_8 + 0.87_real_8*x2 + 0.50_real_8 * x2*x2 + 2.26_real_8 * x2*x2*x2
       dx0 = 2._real_8*xsi * ( 0.87_real_8 + x2 + 3._real_8**2.26_real_8 * x2*x2 )
       IF (ABS(ABS(xsi)-1._real_8).LT.small) THEN
          xx  = 2._real_8**f43
          dxx = f43*2._real_8**f73
          dxx = dxx*SIGN(1._real_8,xsi)
       ELSE
          xx  = (1._real_8 + xsi)**f43 + (1._real_8 - xsi)**f43
          dxx = f43*(1._real_8 + xsi)**f73 - f43*(1._real_8 - xsi)**f73
       ENDIF
       dno = (1._real_8 + 0.5_real_8*eta2*xx)
       cxe = cx0/dno**4
       dcdx = dx0/dno**4 - 4._real_8*cxe*0.5_real_8*eta2*dxx/dno
       dcde = -2._real_8*cxe*xx/dno
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE cgefun
  ! ==================================================================
  SUBROUTINE ccfun(rhoa,rhob,grhoaa,grhobb,grhoab,cxe,dcdra,dcdrb,&
       dcdgaa,dcdgbb,dcdgab)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rhoa, rhob, grhoaa, grhobb, &
                                                grhoab, cxe, dcdra, dcdrb, &
                                                dcdgaa, dcdgbb, dcdgab

    REAL(real_8), PARAMETER                  :: f43 = 4._real_8/3._real_8 , &
                                                ob23 = 2._real_8/3._real_8 , &
                                                small = 1.e-14_real_8 

    REAL(real_8)                             :: cnn, dcde, dcdx, dedgaa, &
                                                dedgab, dedgbb, dedra, dedrb, &
                                                dedx, dxdra, dxdrb, dxsi0, &
                                                eta2, oo, rho, xsi

! ==--------------------------------------------------------------==

    rho    = rhoa + rhob
    xsi    = ( rhoa - rhob ) / rho
    IF (ABS(xsi).GT.1._real_8) xsi=SIGN(1.0_real_8,xsi)
    ! deb
    ! DXSI0  = (XSI-1._real_8)**2 * GRHOAA + 2._real_8*(XSI*XSI-1._real_8)*GRHOAB
    ! *           + (XSI+1._real_8)**2 * GRHOBB
    dxsi0  = (xsi-1._real_8)**2 * grhoaa + (xsi+1._real_8)**2 * grhobb
    ! deb
    oo     = 0.25_real_8/(rho*rho)/(3._real_8*pi*pi*rho)**ob23
    eta2   = dxsi0*oo
    CALL cgefun(xsi,eta2,cnn,dcdx,dcde)
    dxdra  =  2._real_8*rhob/(rho*rho)
    dxdrb  = -2._real_8*rhoa/(rho*rho)
    IF (dxsi0.LT.small*small) THEN
       dedx   = 0._real_8
       dedgaa = 0._real_8
       dedgbb = 0._real_8
       dedgab = 0._real_8
    ELSE
       ! deb    DEDX   = (2._real_8*(XSI-1._real_8) * GRHOAA + 4._real_8*XSI*GRHOAB
       ! deb *            + 2._real_8*(XSI+1._real_8) * GRHOBB)*OO
       dedx   = (2._real_8*(xsi-1._real_8)*grhoaa+2._real_8*(xsi+1._real_8)*grhobb) * oo
       dedgaa = (xsi-1._real_8)**2 * oo
       dedgbb = (xsi+1._real_8)**2 * oo
       ! deb    DEDGAB = 2._real_8*OO*(XSI*XSI-1._real_8)
       dedgab = 0._real_8
    ENDIF
    dedra  = dedx*dxdra - 2._real_8*f43*eta2/rho
    dedrb  = dedx*dxdrb - 2._real_8*f43*eta2/rho
    ! 
    dcdra  = dcdx*dxdra + dcde*dedra
    dcdrb  = dcdx*dxdrb + dcde*dedrb
    dcdgaa = dcde*dedgaa
    dcdgbb = dcde*dedgbb
    dcdgab = dcde*dedgab
    cxe    = cnn
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE ccfun
  ! ==================================================================
  SUBROUTINE revpkzb_lsd(rhoa,rhob,grhoaa,grhobb,grhoab,z,e,&
       dedra,dedrb,dedga,dedgb,dedgab,dedz)
    REAL(real_8)                             :: rhoa, rhob, grhoaa, grhobb, &
                                                grhoab, z, e, dedra, dedrb, &
                                                dedga, dedgb, dedgab, dedz

    REAL(real_8), PARAMETER                  :: c00 = 0.54_real_8 , &
                                                ob3 = 1._real_8/3._real_8 , &
                                                small = 1.e-14_real_8 

    REAL(real_8) :: cxe, dcdgaa, dcdgab, dcdgbb, dcdra, dcdrb, e2, eca, ecb, &
      ef1, ef2, rho, sc, sca, scb, t1a, t1b, t2a, t2aba, t2abb, t2b, tv1ca, &
      tv1cb, tv2ca, tv2cb, txaba, txabb, tyaba, tyabb, v1ca, v1cb, v2ca, &
      v2cab, v2cb

! ==--------------------------------------------------------------==

    CALL tpbecs(rhoa,rhob,grhoaa,grhoab,grhobb,sc,&
         v1ca,v2ca,v1cb,v2cb,v2cab)
    CALL tpbeca(rhoa,grhoaa,sca,tv1ca,tv2ca)
    CALL tpbeca(rhob,grhobb,scb,tv1cb,tv2cb)
    IF (sca.GT.sc) THEN
       eca = sca
       t1a = tv1ca
       IF (grhoaa.GT.small) THEN
          t2a = tv2ca/SQRT(grhoaa)
       ELSE
          t2a = 0._real_8
       ENDIF
       t2aba = 0._real_8
       txaba = 0._real_8
       tyaba = 0._real_8
    ELSE
       eca = sc
       t1a = v1ca
       t2a = v2ca
       t2aba = v2cab
       txaba = v1cb
       tyaba = v2cb
    ENDIF
    IF (scb.GT.sc) THEN
       ecb = scb
       t1b = tv1cb
       IF (grhobb.GT.small) THEN
          t2b = tv2cb/SQRT(grhobb)
       ELSE
          t2b = 0._real_8
       ENDIF
       t2abb = 0._real_8
       txabb = 0._real_8
       tyabb = 0._real_8
    ELSE
       ecb = sc
       t1b = v1cb
       t2b = v2cb
       t2abb = v2cab
       txabb = v1ca
       tyabb = v2ca
    ENDIF
    ! 
    CALL ccfun(rhoa,rhob,grhoaa,grhobb,grhoab,cxe,dcdra,dcdrb,&
         dcdgaa,dcdgbb,dcdgab)
    ! 
    rho    = rhoa + rhob
    ef1    = 1._real_8 + cxe * z * z
    ef2    = ( 1._real_8 + cxe ) * z * z
    e2     = rhoa/rho * eca + rhob/rho * ecb
    e      = sc * ef1 - ef2 * e2
    dedra  = v1ca * ef1 - (-e2+eca + rhoa*t1a + rhob*txabb)/rho * ef2&
         + dcdra*z*z * ( sc - e2 )
    dedrb  = v1cb * ef1 - (-e2+ecb + rhob*t1b + rhoa*txaba)/rho * ef2&
         + dcdrb*z*z * ( sc - e2 )
    dedga  = v2ca * ef1 - ef2 * ( rhoa/rho*t2a + rhob/rho*tyabb )&
         + 2._real_8*dcdgaa*z*z * ( sc - e2 )
    dedgb  = v2cb * ef1 - ef2 * ( rhob/rho*t2b + rhoa/rho*tyaba )&
         + 2._real_8*dcdgbb*z*z * ( sc - e2 )
    dedgab = v2cab * ef1 - ef2 * (rhoa/rho*t2aba + rhob/rho*t2abb)&
         + 2._real_8*dcdgab*z*z * ( sc - e2 )
    dedz   = 2._real_8 * sc * cxe * z - 2._real_8 * ( 1._real_8 + cxe ) * z * e2
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE revpkzb_lsd
  ! ==================================================================
  SUBROUTINE tpbecs(rhoa,rhob,grhoaa,grhoab,grhobb,sc,&
       v1ca,v2ca,v1cb,v2cb,v2cab)
    REAL(real_8)                             :: rhoa, rhob, grhoaa, grhoab, &
                                                grhobb, sc, v1ca, v2ca, v1cb, &
                                                v2cb, v2cab

    REAL(real_8), PARAMETER :: be = 0.06672455060314922_real_8, &
      ga = 0.031090690869654895_real_8 , ob3 = 1._real_8/3._real_8 , &
      small = 1.e-24_real_8 

    INTEGER                                  :: iflg
    REAL(real_8) :: a, aa, af, dadra, dadrb, decdra, decdrb, dedra, dedrb, &
      dhdra, dhdrb, dhdt, dphida, dphidb, dphide, dsda, dsdra, dsdrb, dsdt, &
      dtdphi, dtdra, dtdrb, ec, eta, expe, grho, h0, phi, phi3, rho, rs, s1, &
      t, vca, vcb, xkf, xks, xy, y

! ==--------------------------------------------------------------==

    rho=rhoa+rhob
    eta=(rhoa-rhob)/rho
    grho=grhoaa+2._real_8*grhoab+grhobb
    IF (ABS(eta).GT.1._real_8) eta=SIGN(1.0_real_8,eta)
    rs = (0.75_real_8/(pi*rho))**ob3
    iflg=2
    IF (rs.LT.1.0_real_8) iflg=1
    CALL lsd_pz(rs,eta,ec,vca,vcb,iflg)
    decdra= (vca-ec)/rho
    decdrb= (vcb-ec)/rho
    phi   = 0.5_real_8*((1._real_8+eta)**(2._real_8/3._real_8)+(1._real_8-eta)**(2._real_8/3._real_8))
    phi3  = phi*phi*phi
    aa    = MAX(grho,small)
    a     = SQRT(aa)
    xkf   = (9._real_8*pi/4._real_8)**ob3/rs
    xks   = SQRT(4._real_8*xkf/pi)
    t     = a/(2._real_8*xks*rho*phi)
    expe  = EXP(-ec/(phi3*ga))
    af    = be/ga * (1._real_8/(expe-1._real_8))
    y     = af*t*t
    xy    = (1._real_8+y)/(1._real_8+y+y*y)
    s1    = 1._real_8+be/ga*t*t*xy
    h0    = ga*phi3 * LOG(s1)
    sc    = ec + h0
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
       dadra = af*af*expe/(-be*phi3)*(3._real_8*ec/phi*dphida-decdra)
       dadrb = af*af*expe/(-be*phi3)*(3._real_8*ec/phi*dphidb-decdrb)
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
       dadra = af*af*expe/(-be*phi3)*(-decdra)
       dadrb = af*af*expe/(-be*phi3)*(-decdrb)
       dsda  = -be/ga * af * t**6 * (2._real_8+y) / (1._real_8+y+y*y)**2
       dsdt  = 2._real_8*be/ga * t * (1._real_8+2._real_8*y) / (1._real_8+y+y*y)**2
       dsdra = dsda*dadra + dsdt*dtdra
       dsdrb = dsda*dadrb + dsdt*dtdrb
       dhdt  = ga*phi3/s1*dsdt
       dhdra = ga*phi3/s1*dsdra
       dhdrb = ga*phi3/s1*dsdrb
    ENDIF
    ! ..
    v1ca  = decdra + dhdra
    v2ca  = dhdt*t/aa
    v1cb  = decdrb + dhdrb
    v2cb  = dhdt*t/aa
    v2cab = dhdt*t/aa
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE tpbecs
  ! ==================================================================


END MODULE metafun_utils
