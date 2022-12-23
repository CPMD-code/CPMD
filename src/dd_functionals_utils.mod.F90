MODULE dd_functionals_utils
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: b88_x
  PUBLIC :: lyp88_c
  PUBLIC :: lyp88_c_loc
  PUBLIC :: slater_x
  PUBLIC :: p86_c
  PUBLIC :: p86_c_loc
  PUBLIC :: pbe96_c
  PUBLIC :: pbe96_c_loc
  PUBLIC :: pbe96_x

CONTAINS

  ! ==================================================================
  SUBROUTINE b88_x(dens,degrad,dfdn,dfdg,d2fdn2,d2fdndg,d2fdg2,lda)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: dens, degrad, dfdn, dfdg, &
                                                d2fdn2, d2fdndg, d2fdg2, lda

    REAL(real_8), PARAMETER :: &
      a = -0.930525736349100025002010218071667_real_8 , beta = 0.0042_real_8 

    REAL(real_8) :: beta6, d1, d12, d13, D1p, D1pp, D1s, D1sp, D1ss, degrado, &
      denso, f3, n1, n2, n_103, n_13, n_23, n_43, n_53, n_73, n_83, s, sgrad, &
      sp, sp3, spp, ss, ssp, sss, t1, t2, t3, t4, t5, t53, x, x2, xp, xs, &
      xsp, xss

! ==--------------------------------------------------------------==

    denso=MAX(dens,1.e-24_real_8)
    degrado=MAX(degrad,1.e-49_real_8)

    n2=denso*denso
    n_13=denso**(1._real_8/3._real_8)
    n_23=n_13**2
    n_43=n_13*denso
    n_53=n_23*denso
    n_73=n_43*denso
    n_83=n_53*denso
    n_103=n_73*denso
    n1=1._real_8/denso
    n_73=1._real_8/n_73
    f3=4._real_8/3._real_8
    beta6=6._real_8*beta

    sgrad=SQRT(degrado)
    x=sgrad/n_43
    x2=x*x
    t1=SQRT(x2+1._real_8)
    t5=1._real_8/t1
    t53=t5**3
    s=LOG(x+t1)
    t2=1._real_8+beta6*x*s
    d1=1._real_8/t2
    d12=d1*d1
    d13=d12*d1
    xs=-f3*sgrad*n_73
    xp=1._real_8/n_43
    sp=1._real_8/SQRT(n_83+degrado)
    sp3=sp**3
    ss=t5*xs
    D1s=-d12*beta6*(xs*s+x*ss)
    D1p=-d12*beta6*(xp*s+x*sp)
    dfdn=-beta*sgrad*(d1*(f3*n1*x+2._real_8*xs)+x*D1s)
    dfdg=-beta*(d1*x*2._real_8+n_43*x2*D1p)

    xss=28._real_8/9._real_8*sgrad/n_103
    t3=beta6*(xs*s+x*ss)
    sss=-t53*x*xs**2+t5*xss
    D1ss=-2._real_8*d1*D1s*t3-d12*beta6*(xss*s+2._real_8*xs*ss+x*sss)
    d2fdn2=-beta*sgrad*(xs*(f3*n1*d1+3._real_8*D1s)+x*(-f3/n2*d1+f3*n1&
         *D1s+D1ss)+2._real_8*xss*d1)

    xsp=-f3*n_73
    ssp=-sp3*f3*n_53
    D1sp=-2._real_8*d1*D1p*t3-d12*beta6*(xsp*s+xs*sp+xp*ss+x*ssp)
    d2fdndg=dfdn/sgrad-beta*sgrad*((xp*d1+x*D1p)*f3*n1&
         +2._real_8*xsp*d1+2._real_8*xs*D1p+xp*D1s+x*D1sp)

    spp=-sp3*sgrad
    t4=beta6*(xp*s+x*sp)
    D1pp=-2._real_8*d1*D1p*t4-d12*beta6*(2._real_8*xp*sp+x*spp)
    d2fdg2=-beta*(D1p*2._real_8*x+d1*2._real_8*xp+n_43&
         *(2._real_8*x*xp*D1p+x2*D1pp))

    lda=f3/3._real_8*a/n_23
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE b88_x
  ! ==================================================================
  SUBROUTINE lyp88_c(densa,densb,grada,gradax,graday,gradaz,gradb,&
       gradbx,gradby,gradbz,dfdga,nana,nanb,naga,granb,gagb,nagb,&
       nagrb,ganb)
    ! ==--------------------------------------------------------------==
    REAL(real_8) :: densa, densb, grada, gradax, graday, gradaz, gradb, &
      gradbx, gradby, gradbz, dfdga, nana, nanb, naga, granb, gagb, nagb, &
      nagrb, ganb

    REAL(real_8), PARAMETER :: a = 0.04918_real_8 , b = 0.132_real_8 , &
      c = 0.2533_real_8 , Cf = 36.462398978764777098305310076_real_8 , &
      d = 0.349_real_8 

    REAL(real_8) :: ab, br, c13, D1, d13, ddens, dens, densao, densao2, &
      densbo, densbo2, dn2, e, e3, es, ess, f3, grad, gradao, gradbo, kl, &
      n113, n143, n2, n3, n43, n_113, n_13, n_143, n_173, n_43, n_73, na_23, &
      na_53, na_83, nb_53, nb_83, o18, o3, o9, s18, skalar, t1, t12, t1s, &
      t1ss, t2, t2s, t2ss, t2sst, t2st, t3, t3s, t3ss, t4, t4p, t4s, t4sp, &
      t4ss, t4sst, t4st, t5, t5k, t5p, t5pst, t5s, t5sk, t5sp, t5ss, t5sst, &
      t5st, t6, t6s, t6ss, t7, t7p, t7s, t7ss, t7sst, t7st, t91, t92, tw3

! ==--------------------------------------------------------------==

    densao=MAX(densa,1.e-24_real_8)
    densbo=MAX(densb,1.e-24_real_8)
    gradao=MAX(grada,1.e-48_real_8)
    gradbo=MAX(gradb,1.e-48_real_8)
    dens=densao+densbo
    skalar=gradax*gradbx+graday*gradby+gradaz*gradbz
    grad=gradao+gradbo+2*skalar
    n2=dens*dens
    n3=1._real_8/n2/dens
    densao2=densao*densao
    densbo2=densbo*densbo
    o3=1._real_8/3._real_8
    tw3=o3+o3
    f3=tw3+tw3
    e3=1._real_8+f3+f3
    o9=o3*o3
    o18=o9/2._real_8
    s18=o3+o18
    n_13=dens**o3
    n_43=n_13*dens
    n_73=1._real_8/n_43/dens
    n_113=n_43*n_43*dens
    n_143=n_113*dens
    n_173=n_143*dens
    n43=1._real_8/n_43
    n113=1._real_8/n_113
    n143=1._real_8/n_143
    na_23=densao**tw3
    na_53=na_23*densao
    na_83=na_53*densao
    nb_53=densbo**(5._real_8/3._real_8)
    nb_83=nb_53*densbo
    c13=c/n_13
    e=EXP(-c13)
    d13=d/n_13
    t1=1._real_8/(1._real_8+d13)
    t12=t1*t1
    ab=densao*densbo
    ddens=1._real_8/dens
    t2=ab*ddens
    t3=a*b*e*t1*n113
    t7=(densao*gradao+densbo*gradbo)*ddens
    t6=c13+d13*t1
    br=o9*(t6-11._real_8)
    kl=(47._real_8*o18-s18*t6)
    t5=Cf*(na_83+nb_83)+kl*grad-(2.5_real_8-o18*t6)*(gradao+gradbo)-br*t7
    t4=ab*t5+tw3*n2*(gradao+gradbo-grad)-densao2*gradbo-densbo2*gradao
    es=c*e*o3*n43
    t1s=o3*d*t12*n43
    dn2=1._real_8/n2
    t2s=densbo2*dn2
    t2st=densao2*dn2
    t3s=a*b*(t1s*e*n113+t1*es*n113-t1*e*e3*n143)
    t7s=densbo*dn2*(gradao-gradbo)
    t7st=densao*dn2*(gradbo-gradao)
    t6s=-o3*n43*(c+d*t1)+d13*t1s
    t91=-s18*grad*t6s-o9*t6s*t7+o18*t6s*(gradao+gradbo)
    t5s=Cf*8._real_8*o3*na_53-br*t7s+t91
    t5st=Cf*8._real_8*o3*nb_53-br*t7st+t91
    t4s=densbo*t5+ab*t5s+(f3*dens-2._real_8*densao)*gradbo&
         +f3*dens*(gradao-grad)
    t4st=densao*t5+ab*t5st+(f3*dens-2._real_8*densbo)&
         *gradao+f3*dens*(gradbo-grad)
    d1=-4._real_8*a*(t1s*t2+t1*t2s)-t3s*t4-t3*t4s
    t7p=densao*ddens
    t5p=o9-o3*t6-br*t7p
    t5k=o9-o3*t6-br*densbo*ddens
    t5sp=-o3*t6s*(1+o3*t7p)-br*densbo*dn2
    t5sk=-o3*t6s*(1+o3*densbo*ddens)+br*densbo*dn2
    t4p=ab*t5p-densbo2
    t4sp=densbo*t5p+ab*t5sp
    dfdga=-t3*t4p
    naga=-t3s*t4p-t3*t4sp
    ess=c*o9*e*n_73*(c13-4._real_8)
    t1ss=o3*d*(2._real_8*t1*t1s*n43-f3*n_73*t12)
    t2ss=-2._real_8*densbo2*n3
    t2sst=2._real_8*ab*n3
    t3ss=a*b*((t1ss*e+2._real_8*t1s*es+t1*ess)*n113-2._real_8*e3*(t1s*e&
         *n143+t1*es*n143-7._real_8*o3*t1*e/n_173))
    t7ss=-2._real_8*densbo*n3*(gradao-gradbo)
    t7sst=(densao-densbo)*n3*(gradao-gradbo)
    t6ss=4._real_8*o9*n_73*(c+d*t1)-tw3*d*n43*t1s+d13*t1ss
    t92=t6ss*(o18*(gradao+gradbo)-o9*t7-s18*grad)
    t5ss=Cf*o9*40._real_8*na_23-2._real_8*o9*t6s*t7s-br*t7ss+t92
    t5sst=-o9*t6s/n2*(gradao-gradbo)*(densbo-densao)-br*t7sst+t92
    t4ss=2._real_8*densbo*t5s+ab*t5ss-f3*(grad+0.5_real_8*gradbo-gradao)
    t4sst=densbo*t5st+t5+densao*t5s+ab*t5sst+f3*(-grad+gradbo+gradao)
    nana=-4._real_8*a*(t1ss*t2+2._real_8*t1s*t2s+t1*t2ss)-t3ss*t4&
         -2._real_8*t3s*t4s-t3*t4ss
    t5pst=-o3*t6s-densao*(o9*t6s*ddens-br*dn2)
    granb=-t3s*t4p-t3*(densao*t5p+ab*t5pst-2._real_8*densbo)
    nanb=-4._real_8*a*(t1ss*t2+t1s*t2st+t1s*t2s+t1*t2sst)&
         -t3ss*t4-t3s*t4st-t3s*t4s-t3*t4sst
    gagb=-t3*(ab*2._real_8*kl-f3*n2)
    ganb=2._real_8*(-t3s*(ab*kl-tw3*n2)-t3*(densao*kl-ab*s18*t6s-f3*dens))
    nagrb=-t3s*(ab*t5k-densao2)-t3*(densbo*t5k+ab*t5sk-2._real_8*densao)
    nagb=2._real_8*(-t3s*(ab*kl-tw3*n2)-t3*(-f3*dens+densbo*kl-ab*s18*t6s))

    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE lyp88_c
  ! ==================================================================
  SUBROUTINE lyp88_c_loc(densa,densb,loca,locb)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: densa, densb, loca, locb

    REAL(real_8), PARAMETER :: a = 0.04918_real_8 , b = 0.132_real_8 , &
      balpha = -0.930525736349100025002010218071667_real_8 , &
      c = 0.2533_real_8 , Cf = 36.462398978764777098305310076_real_8 , &
      d = 0.349_real_8 

    REAL(real_8) :: ab, c13, d13, ddens, dens, densao, densao2, densbo, &
      densbo2, dn2, e, e3, es, ess, f3, n113, n143, n2, n3, n43, n_113, n_13, &
      n_143, n_173, n_43, n_73, na_113, na_23, na_53, na_83, nb_113, nb_83, &
      o3, o9, t1, t12, t1s, t1ss, t2, t2s, t2ss, t2sst, t2st, t3, t3s, t3ss

! ==--------------------------------------------------------------==

    densao=MAX(densa,1.e-24_real_8)
    densbo=MAX(densb,1.e-24_real_8)
    o3=1._real_8/3._real_8
    o9=o3*o3
    f3=4._real_8*o3
    e3=11._real_8*o3
    dens=densao+densbo
    n2=dens*dens
    n3=n2*dens
    densao2=densao*densao
    densbo2=densbo*densbo
    n_13=dens**o3
    n_43=n_13*dens
    n_73=1._real_8/n_43/dens
    n_113=n_43*n_43*dens
    n_143=n_113*dens
    n_173=n_143*dens
    n43=1._real_8/n_43
    n113=1._real_8/n_113
    n143=1._real_8/n_143
    na_23=densao**(2._real_8*o3)
    na_53=na_23*densao
    na_83=na_53*densao
    na_113=na_83*densao
    nb_83=densbo**(8._real_8*o3)
    nb_113=nb_83*densbo
    c13=c/n_13
    e=EXP(-c13)
    d13=d/n_13
    t1=1._real_8/(1._real_8+d13)
    t12=t1*t1
    ab=densao*densbo
    ddens=1._real_8/dens
    t2=ab*ddens
    t3=a*b*e*t1*n113
    es=c*e*o3*n43
    t1s=o3*d*t12*n43
    dn2=1._real_8/n2
    t2s=densbo2*dn2
    t2st=densao2*dn2
    t3s=a*b*(t1s*e*n113+t1*es*n113-t1*e*e3*n143)
    ess=c*o9*e*n_73*(c13-4._real_8)
    t1ss=o3*d*(2._real_8*t1*t1s*n43-f3*n_73*t12)
    t2ss=-2._real_8*densbo2/n3
    t2sst=2._real_8*ab/n3
    t3ss=a*b*((t1ss*e+2._real_8*t1s*es+t1*ess)*n113-2._real_8*e3*(t1s*e&
         *n143+t1*es*n143-7._real_8*o3*t1*e/n_173))

    loca=-4._real_8*a*(t1ss*t2+2._real_8*t1s*t2s+t1*t2ss)-Cf*t3ss*(na_113*densbo&
         +densao*nb_113)-2._real_8*Cf*t3s*(11._real_8/3._real_8*na_83*densbo+nb_113)&
         -Cf*t3*88._real_8/9._real_8*na_53*densbo
    locb=-4._real_8*a*(t1ss*t2+t1s*t2st+t1s*t2s+t1*t2sst)&
         -Cf*t3ss*(na_113*densbo+densao*nb_113)-cf*t3s*(na_113+nb_113&
         +e3*(densao*nb_83+densbo*na_83))-Cf*t3*e3*(na_83+nb_83)

    loca=loca+4._real_8*o9*balpha/na_23
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE lyp88_c_loc
  ! ==================================================================
  SUBROUTINE slater_x(dens,lda)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: dens, lda

    REAL(real_8), PARAMETER :: &
      a = -0.930525736349100025002010218071667_real_8 

    REAL(real_8)                             :: rho

    rho=MAX(dens,1.e-24_real_8)
    rho=rho**(2._real_8/3._real_8)
    lda=4._real_8/9._real_8*a/rho
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE slater_x
  ! ==================================================================
  SUBROUTINE p86_c(densa,densb,grada,gradax,graday,gradaz,gradb,&
       gradbx,gradby,gradbz,nana,naga,gaga,nask,gask,sk,ga,sksk,gagb&
       ,nagb,skgb,sknb,ganb,nanb)
    ! ==--------------------------------------------------------------==
    REAL(real_8) :: densa, densb, grada, gradax, graday, gradaz, gradb, &
      gradbx, gradby, gradbz, nana, naga, gaga, nask, gask, sk, ga, sksk, &
      gagb, nagb, skgb, sknb, ganb, nanb

    REAL(real_8), PARAMETER :: Af = 0.01555_real_8, Ap = 0.0311_real_8, &
      B1f = 1.3981_real_8, B1p = 1.0529_real_8, B2f = 0.2611_real_8 , &
      B2p = 0.3334_real_8, Bf = -0.0269_real_8, Bp = -0.048_real_8, &
      Cf = 0.0007_real_8, Cp = 0.002_real_8, Df = -0.0048_real_8 , &
      Dp = -0.0116_real_8, f3 = 4._real_8/3._real_8 , &
      fakt = 0.6203504908994000_real_8 , fakt2 = 0.519842099789746_real_8 , &
      fakt3 = 0.629960524947436_real_8 , gf = -0.0843_real_8, &
      gp = -0.1423_real_8, k = 8.1290825e-4_real_8, o3 = 1._real_8/3._real_8, &
      p0 = 0.002568_real_8, p1 = 0.023266_real_8, p2 = 7.389e-6_real_8, &
      p3 = 8.723_real_8
    REAL(real_8), PARAMETER :: P4 = 0.472_real_8, P5 = 7.389e-2_real_8 

    REAL(real_8) :: C, C2, C3, Cs, Cs2, Css, d, d2, d3, dens, densao, densbo, &
      ds, ds2, dsb, dss, dssb, e, epsf, epsfs, epsfss, epsp, epsps, epspss, &
      f, fs, fsb, fss, fssb, grad, gradao, gradbo, n2, n3, n_103, n_13, &
      n_136, n_196, n_43, n_73, n_76, Phi, Phip, Phipb, Phipp, Phippb, Phis, &
      Phisch, Phischsch, Phisp, Phiss, Phissch, rs, rs2, rs3, rsq, rss, rsss, &
      sgrad, sgrada, sgradb, skalar, t1, t12, t2, t3, t4, t5, t5s, t5sb, t6f, &
      t6p, Z, zeta, zetas, zetas2, zetasb, zetass, zetassb, zm13, zm23, zm43, &
      zm53, zp13, zp23, zp43, zp53, Zs, Zsb, Zss, Zssb

    densao=MAX(densa,1.e-24_real_8)
    densbo=MAX(densb,1.e-24_real_8)
    gradao=MAX(grada,1.e-48_real_8)
    gradbo=MAX(gradb,1.e-48_real_8)

    dens=densao+densbo
    n2=dens*dens
    n3=n2*dens
    n_13=dens**o3
    n_43=n_13*dens
    n_73=n_43*dens
    n_103=n_73*dens
    n_76=SQRT(n_73)
    n_136=n_76*dens
    n_196=n_136*dens
    skalar=gradax*gradbx+graday*gradby+gradaz*gradbz
    grad=gradao+gradbo+2._real_8*skalar

    zeta=(densao-densbo)/dens
    zp13=(1._real_8+zeta)**o3
    zm13=(1._real_8-zeta)**o3
    zp23=zp13*zp13
    zm23=zm13*zm13
    zp43=zp23*zp23
    zm43=zm23*zm23
    zp53=zp43*zp13
    zm53=zm43*zm13
    rs=fakt/n_13
    rs2=rs*rs
    rs3=rs2*rs
    rss=-o3*fakt/n_43
    rsss=4._real_8/9._real_8*fakt/n_73
    f=(zp43+zm43-2._real_8)/fakt2

    IF (.TRUE.) THEN
       ! ==------------------------------------------------------------==
       ! Discerns between high- and low-density (as in subroutine LSD_PZ)
       ! ==------------------------------------------------------------==
       ! Low-density formula
       IF (rs.GE.1._real_8) THEN
          rsq=SQRT(rs)
          epsp=gp/(1._real_8+B1p*rsq+B2p*rs)
          epsf=gf/(1._real_8+B1f*rsq+B2f*rs)
          t6p=B1p/(2._real_8*rsq)
          t6f=B1f/(2._real_8*rsq)
          epsps=-epsp**2/gp*(t6p+B2p)*rss
          epsfs=-epsf**2/gf*(t6f+B2f)*rss
          epspss=-1._real_8/gp*epsp*(2._real_8*epsps*(t6p+B2p)*rss+epsp*(-t6p/&
               (2._real_8*rs)*rss+(t6p+B2p)*rsss))
          epsfss=-1._real_8/gf*epsf*(2._real_8*epsfs*(t6f+B2f)*rss+epsf*(-t6f/&
               (2._real_8*rs)*rss+(t6f+B2f)*rsss))
          ! High-density formula
       ELSE
          epsp=Ap*LOG(rs)+Bp+Cp*rs*LOG(rs)+Dp*rs
          epsf=Af*LOG(rs)+Bf+Cf*rs*LOG(rs)+Df*rs
          epsps=(Ap/rs+Cp*(LOG(rs)+1._real_8)+Dp)*rss
          epsfs=(Af/rs+Cf*(LOG(rs)+1._real_8)+Df)*rss
          epspss=(-Ap/rs2+Cp/rs)*rss*rss+epsps/rss*rsss
          epsfss=(-Af/rs2+Cf/rs)*rss*rss+epsfs/rss*rsss
       ENDIF
    ELSE
       ! ==------------------------------------------------------------==
       ! High-density formula only (as in subroutine GCSXP86)
       ! ==------------------------------------------------------------==
       epsp=Ap*LOG(rs)+Bp+Cp*rs*LOG(rs)+Dp*rs
       epsf=Af*LOG(rs)+Bf+Cf*rs*LOG(rs)+Df*rs
       epsps=(Ap/rs+Cp*(LOG(rs)+1._real_8)+Dp)*rss
       epsfs=(Af/rs+Cf*(LOG(rs)+1._real_8)+Df)*rss
       epspss=(-Ap/rs2+Cp/rs)*rss*rss+epsps/rss*rsss
       epsfss=(-Af/rs2+Cf/rs)*rss*rss+epsfs/rss*rsss
    ENDIF
    z=epsp+f*(epsf-epsp)
    t1=1._real_8+p3*rs+p4*rs2+p5*rs3
    t12=t1*t1
    t2=p0+p1*rs+p2*rs2
    c=0.001667_real_8+t2/t1
    c2=c*c
    c3=c2*c
    d=2._real_8*fakt3*SQRT(0.5_real_8*fakt3*(zp53+zm53))
    d2=d*d
    d3=d2*d
    sgrada=SQRT(gradao)
    sgradb=SQRT(gradbo)
    sgrad=SQRT(grad)
    Phi=k*sgrad/n_76/c
    e=EXP(-Phi)

    zetas=2._real_8*densbo/n2
    zetas2=zetas*zetas
    zetasb=-zetas*densao/densbo

    fs=f3*zetas*(zp13-zm13)/fakt2
    fsb=-fs*densao/densbo
    Zs=epsps+fs*(epsf-epsp)+f*(epsfs-epsps)
    Zsb=epsps+fsb*(epsf-epsp)+f*(epsfs-epsps)
    Cs=(p1+2._real_8*p2*rs)/t1-t2*(p3+2._real_8*p4*rs+3*p5*rs2)&
         /t12
    Cs=cs*rss
    Cs2=Cs*cs
    ds=1._real_8/(d*2._real_8)*5._real_8/6._real_8*zetas*(zp23-zm23)
    ds2=ds*ds

    dsb=-ds*densao/densbo
    Phis=-k*sgrad*(7._real_8/6._real_8/(c*n_136)+Cs/(c2*n_76))

    Phip=k/(n_76*c*sgrad)*sgrada
    Phipb=k/(n_76*c*sgrad)*sgradb
    t4=e*(-Phip*grad+2._real_8*sgrada)
    ga=c/d/n_43*t4

    zetass=-4._real_8*densbo/n3
    zetassb=2._real_8*(densao-densbo)/n3

    fss=f3/fakt2*(o3/zp23*zetas2+zp13*zetass+o3/zm23*zetas2-&
         zm13*zetass)
    fssb=f3/fakt2*(zetassb*(zp13-zm13)+o3*zetas*zetasb*&
         (1._real_8/zp23+1._real_8/zm23))

    Zss=epspss+fss*(epsf-epsp)+2._real_8*fs*(epsfs-epsps)+f*(epsfss-epspss)
    t3=p3+2._real_8*p4*rs+3._real_8*p5*rs2
    Zssb=epspss+fssb*(epsf-epsp)+fs*(epsfs-epsps)*(1._real_8-densao/densbo)&
         +f*(epsfss-epspss)
    Css=(2*p2/t1-(2._real_8*(p1+2._real_8*p2*rs)*t3+t2*(2._real_8*p4+6._real_8*p5*rs))&
         /t12+1/(t12*t1)*2._real_8*t2*t3*t3)*rss*rss+Cs/rss*rsss
    dss=-ds**2/d+0.5_real_8/(fakt3*d)*5._real_8/3._real_8*(2._real_8*fakt3*o3/2._real_8*&
         zetas**2*(1._real_8/zp13+1._real_8/zm13)+fakt3*(zp23-zm23)*zetass/2._real_8)
    dssb=5._real_8/12._real_8*(zetas/d*(zp23-zm23)*(-dsb/d+zetassb/zetas)+&
         zetas/d*2._real_8*o3*zetasb*(1._real_8/zp13+1._real_8/zm13))
    Phiss=k*sgrad*(91._real_8/36._real_8/n_196/c+7._real_8*o3/n_136/c2*Cs+&
         2._real_8/C3*Cs2/n_76-Css/C2/n_76)
    t5=1._real_8/d*(-Phis*c+Cs-c*(ds/d+f3/dens))
    t5s=-ds/d2*(-Phis*c+Cs-c/d*ds-f3*c/dens)+(-Phiss*c-phis*cs+Css-cs*&
         (ds/d+f3/dens)-C*(-ds2/d2+dss/d-f3/n2))/d
    nana=2._real_8*Zs+dens*Zss+grad*e*(-Phis/n_43*t5-f3/n_73*t5+t5s/n_43)
    t5sb=-dsb/d2*(-Phis*c+Cs-c/d*ds-f3*c/dens)+(-Phiss*c-phis*cs+Css-&
         Cs*(ds/d+f3/dens)-C*(-ds*dsb/d2+dssb/d-f3/n2))/d

    nanb=Zsb+Zs+dens*Zssb+grad*e*(-Phis/n_43*t5-f3/n_73*t5+t5sb/n_43)

    Phisp=k*(-7._real_8/6._real_8/n_136/c-Cs/c2/n_76)/sgrad*sgrada
    naga=1._real_8/(n_43*d)*(Cs*t4+c*e*(Phis*Phip*grad-Phisp*grad-&
         Phis*2._real_8*sgrada))+C*t4*(-f3/n_73/d-ds/d2/n_43)

    ganb=1._real_8/(n_43*d)*(Cs*t4+c*e*(Phis*Phip*grad-Phisp*grad-&
         Phis*2._real_8*sgrada))+C*t4*(-f3/n_73/d-dsb/d2/n_43)

    nagb=naga*sgradb/sgrada

    Phipp=k/(n_76*c)*(1._real_8/sgrad-sgrada**2/(grad*sgrad))
    gaga=c/(n_43*d)*e*((Phip**2-Phipp)*grad-phip*2._real_8*sgrada+2._real_8-&
         2._real_8*sgrada*Phip)

    Phippb=k/(n_76*c)*(-sgrada*sgradb/(grad*sgrad))
    gagb=c/(n_43*d)*e*((Phipb*Phip-Phippb)*grad-phip*2._real_8*sgradb&
         -Phipb*2._real_8*sgrada)

    Phisch=k/(c*n_76*sgrad)
    sk=ga/sgrada

    Phissch=k*(-7._real_8/6._real_8/n_136/c-Cs/c2/n_76)/sgrad
    nask=nagb/sgradb

    sknb=1._real_8/(n_43*d)*e*(-Phisch*grad+2._real_8)*(Cs-f3*c/dens-c/d*dsb&
         -C*Phis)-C/(n_43*d)*e*Phissch*grad

    Phischsch=-k/(n_76*c)/(grad*sgrad)
    sksk=c/(n_43*d)*e*(Phisch**2*grad-Phischsch*grad-4._real_8*phisch)

    gask=gagb/sgradb

    skgb=sksk*sgradb
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE p86_c
  ! ==================================================================
  SUBROUTINE p86_c_loc(densa,densb,p86a,p86b)
    ! ==================================================================
    REAL(real_8)                             :: densa, densb, p86a, p86b

    REAL(real_8), PARAMETER :: Af = 0.01555_real_8, Ap = 0.0311_real_8, &
      B1f = 1.3981_real_8, B1p = 1.0529_real_8, B2f = 0.2611_real_8 , &
      B2p = 0.3334_real_8, Bf = -0.0269_real_8, Bp = -0.048_real_8, &
      Cf = 0.0007_real_8, Cp = 0.002_real_8, Df = -0.0048_real_8 , &
      Dp = -0.0116_real_8, f3 = 4._real_8/3._real_8 , &
      fakt = 0.6203504908994000_real_8 , fakt2 = 0.519842099789746_real_8 , &
      gf = -0.0843_real_8, gp = -0.1423_real_8, o3 = 1._real_8/3._real_8, &
      p0 = 0.002568_real_8, p1 = 0.023266_real_8, p2 = 7.389e-6_real_8, &
      p3 = 8.723_real_8, P4 = 0.472_real_8, P5 = 7.389e-2_real_8 

    REAL(real_8) :: dens, densao, densbo, epsf, epsfs, epsfss, epsp, epsps, &
      epspss, f, fs, fsb, fss, fssb, n2, n3, n_13, n_43, n_73, rs, rs2, rsq, &
      rss, rsss, t6f, t6p, zeta, zetas, zetas2, zetasb, zetass, zetassb, &
      zm13, zm23, zm43, zp13, zp23, zp43, Zs, Zsb, Zss, Zssb

    densao=MAX(densa,1.e-24_real_8)
    densbo=MAX(densb,1.e-24_real_8)
    dens=densao+densbo
    n2=dens*dens
    n3=n2*dens
    n_13=dens**o3
    n_43=n_13*dens
    n_73=n_43*dens
    rs=fakt/n_13
    rs2=rs*rs
    rss=-o3*fakt/n_43
    rsss=4._real_8/9._real_8*fakt/n_73
    zeta=(densao-densbo)/dens
    zetas=2._real_8*densbo/n2
    zetas2=zetas*zetas
    zetasb=-zetas*densao/densbo
    zp13=(1._real_8+zeta)**o3
    zm13=(1._real_8-zeta)**o3
    zp23=zp13*zp13
    zm23=zm13*zm13
    zp43=zp23*zp23
    zm43=zm23*zm23
    f=(zp43+zm43-2._real_8)/fakt2
    fs=f3*zetas*(zp13-zm13)/fakt2
    fsb=-fs*densao/densbo
    rsq=SQRT(rs)
    epsp=gp/(1._real_8+B1p*rsq+B2p*rs)
    epsf=gf/(1._real_8+B1f*rsq+B2f*rs)
    t6p=B1p/(2._real_8*rsq)
    t6f=B1f/(2._real_8*rsq)
    epsps=-epsp**2/gp*(t6p+B2p)*rss
    epsfs=-epsf**2/gf*(t6f+B2f)*rss
    epspss=-1._real_8/gp*epsp*(2._real_8*epsps*(t6p+B2p)*rss+epsp*(-t6p/&
         (2._real_8*rs)*rss+(t6p+B2p)*rsss))
    epsfss=-1._real_8/gf*epsf*(2._real_8*epsfs*(t6f+B2f)*rss+epsf*(-t6f/&
         (2._real_8*rs)*rss+(t6f+B2f)*rsss))

    Zs=epsps+fs*(epsf-epsp)+f*(epsfs-epsps)
    Zsb=epsps+fsb*(epsf-epsp)+f*(epsfs-epsps)

    zetass=-4._real_8*densbo/n3
    zetassb=2._real_8*(densao-densbo)/n3

    fss=f3/fakt2*(o3/zp23*zetas2+zp13*zetass+o3/zm23*zetas2-&
         zm13*zetass)
    fssb=f3/fakt2*(zetassb*(zp13-zm13)+o3*zetas*zetasb*&
         (1._real_8/zp23+1._real_8/zm23))

    Zss=epspss+fss*(epsf-epsp)+2._real_8*fs*(epsfs-epsps)+f*(epsfss-epspss)
    Zssb=epspss+fssb*(epsf-epsp)+fs*(epsfs-epsps)*(1._real_8-densao/densbo)&
         +f*(epsfss-epspss)

    p86a=2._real_8*Zs+dens*Zss
    p86b=Zsb+Zs+dens*Zssb
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE p86_c_loc
  ! ==================================================================
  SUBROUTINE pbe96_c(densa,densb,grada,gradax,graday,gradaz,gradb,&
       gradbx,gradby,gradbz,nana,naga,gaga,nask,gask,sk,ga,sksk,gagb&
       ,nagb,skgb,sknb,ganb,nanb)
    ! ==================================================================
    REAL(real_8) :: densa, densb, grada, gradax, graday, gradaz, gradb, &
      gradbx, gradby, gradbz, nana, naga, gaga, nask, gask, sk, ga, sksk, &
      gagb, nagb, skgb, sknb, ganb, nanb

    REAL(real_8), PARAMETER :: c = 1.709921_real_8, &
      fakt = 0.6203504908994000_real_8, &
      iota = 0.071599657785951923081_real_8, &
      k1 = 1.0077158813689795508_real_8, k2 = 0.031090690869654895_real_8 , &
      k3 = -41.156312099041145348_real_8, k4 = 2.1461263399673646085_real_8 , &
      lambda = 0.06672455060314922_real_8, Lg = 8.723_real_8, &
      ny = 15.755501913376439197_real_8, &
      o = 6.46309613581743018602522_real_8, &
      o3 = 0.3333333333333333333333_real_8, Phig = 0.007389_real_8, &
      t1 = 0.031091_real_8, t2 = 0.015545_real_8, T3 = 0.016887_real_8, &
      U1 = 0.2137_real_8, U2 = 0.20548_real_8, U3 = 0.11125_real_8
    REAL(real_8), PARAMETER :: V1 = 7.5957_real_8, V2 = 14.1189_real_8, &
      V3 = 10.357_real_8, W1 = 3.5876_real_8, W2 = 6.1977_real_8, &
      W3 = 3.6231_real_8, X1 = 1.6382_real_8, X2 = 3.3662_real_8, &
      X3 = 0.88026_real_8, Xi = 23.266_real_8, Y = 0.472_real_8, &
      Y1 = 0.49294_real_8, Y2 = 0.62517_real_8, Y3 = 0.49671_real_8, &
      z = -0.001667_real_8

    REAL(real_8) :: A, A2, Aa, Aaa, Aab, Ab, b, B2, Ba, Baa, Bab, Bb, Bga, &
      Bgaga, Bganb, Bnaga, d, d2, d3, d4, da, daa, dab, db, dens, densao, &
      densbo, e1, e1a, e1aa, e2, e2a, e2aa, e3, e3a, e3aa, eps, epsa, epsaa, &
      epsab, epsb, ex, exa, exaa, exab, exb, f1, f12, f1a, f1aa, f2, f22, &
      f2a, f2aa, f3, f32, f3a, f3aa, g1, g2, g3, grad, grad2, gradao, gradbo, &
      Ja, Jaa, Jab, Jb, Jga, Jgaga, Jganb, Jnaga, kl1, kl1a, kl1b, kl2, kl2a, &
      kl2b, L, La, Laa, Lab, Lb, Lga, Lgaga, Lganb, Lnaga, n2, n3, n_103, &
      n_13, n_136, n_196, n_43, n_73, n_76, na_13, na_23, nb_13, nb_23, ne1, &
      ne12, ne1a, ne1aa, ne1ab, ne1b, ne1ga
    REAL(real_8) :: ne1gaga, ne1ganb, ne1naga, ome, omea, omeaa, omeab, omeb, &
      rs, rs12, rs2, rs3, rsa, rsaa, sgrad, skalar, t, tt, ttb, u, ua, uaa, &
      uab, ub, uh2, uh3, uh4, uh5, ze1, ze1a, ze1aa, ze1ab, ze1b, ze1ga, &
      ze1gaga, ze1ganb, ze1naga, zeta, zeta2, zeta3, zeta4, zetaa, zetaaa, &
      zetaab, zetab, zm13, zm23, zm43, zm53, zp13, zp23, zp43, zp53

! ==--------------------------------------------------------------==

    densao=MAX(densa,1.e-24_real_8)
    densbo=MAX(densb,1.e-24_real_8)
    gradao=MAX(grada,1.e-48_real_8)
    gradbo=MAX(gradb,1.e-48_real_8)

    dens=densao+densbo
    n2=dens*dens
    n3=n2*dens
    n_13=dens**o3
    n_43=n_13*dens
    n_73=n_43*dens
    n_76=SQRT(n_73)
    n_136=n_76*dens
    n_196=n_136*dens
    n_103=n_73*dens
    na_13=densao**o3
    na_23=na_13*na_13
    nb_13=densbo**o3
    nb_23=nb_13*nb_13
    skalar=gradax*gradbx+graday*gradby+gradaz*gradbz
    grad=gradao+gradbo+2._real_8*skalar
    grad2=grad*grad
    sgrad=SQRT(grad)

    rs=fakt/n_13
    rs2=rs*rs
    rs3=rs2*rs
    rs12=SQRT(rs)
    rsa=-o3*fakt/n_43
    rsaa=4._real_8/9._real_8*fakt/n_73
    zeta=(densao-densbo)/dens
    zeta2=zeta*zeta
    zeta3=zeta2*zeta
    zeta4=zeta3*zeta
    zetaa=2._real_8*densbo/n2
    zp13=(1._real_8+zeta)**o3
    zm13=(1._real_8-zeta)**o3
    zp23=zp13*zp13
    zm23=zm13*zm13
    zp43=zp23*zp23
    zm43=zm23*zm23
    zp53=zp43*zp13
    zm53=zm43*zm13
    ome=(zp43+zm43-2._real_8)/(2._real_8**(4._real_8/3._real_8)-2._real_8)
    omea=o*densbo*(na_13-nb_13)/n_73

    f1=2._real_8*t1*(v1*rs12+w1*rs+x1*rs*rs12+y1*rs2)
    f12=f1*f1
    g1=1._real_8+1._real_8/f1
    f1a=2._real_8*t1*(0.5_real_8*v1/rs12+w1+3._real_8/2._real_8*x1*rs12+2._real_8*y1*rs)*rsa
    e1=-2._real_8*t1*(1._real_8+u1*rs)*LOG(g1)
    e1a=-2._real_8*t1*(u1*rsa*LOG(g1)-(1._real_8+u1*rs)/g1/f12*f1a)

    f2=2._real_8*t2*(v2*rs12+w2*rs+x2*rs*rs12+y2*rs2)
    f22=f2*f2
    g2=(f2+1._real_8)/f2
    f2a=2._real_8*t2*(0.5_real_8*v2/rs12+w2+3._real_8/2._real_8*x2*rs12+2._real_8*y2*rs)*rsa
    e2=-2._real_8*t2*(1._real_8+u2*rs)*LOG(g2)
    e2a=-2._real_8*t2*(u2*rsa*LOG(g2)-(1._real_8+u2*rs)/g2/f22*f2a)

    f3=2._real_8*t3*(v3*rs12+w3*rs+x3*rs*rs12+y3*rs2)
    f32=f3*f3
    g3=(f3+1._real_8)/f3
    f3a=2._real_8*t3*(0.5_real_8*v3/rs12+w3+3._real_8/2._real_8*x3*rs12+2._real_8*y3*rs)*rsa
    e3=-2._real_8*t3*(1._real_8+u3*rs)*LOG(g3)
    e3a=-2._real_8*t3*(u3*rsa*LOG(g3)-(1._real_8+u3*rs)/g3/f32*f3a)

    zetaaa=-4._real_8*densbo/n3
    zetaab=2._real_8*zeta/n2
    omeaa=o*densbo*(o3/na_23/n_73-(na_13-nb_13)*7._real_8/3._real_8/n_103)
    omeab=o*((na_13-nb_13)/n_73+densbo*(-o3/nb_23/n_73-(na_13-nb_13)&
         *7._real_8/3._real_8/n_103))
    t=(1._real_8-zeta4)/c
    tt=4._real_8*zeta3*zetaa

    eps=e1-e3*ome*t+(e2-e1)*ome*zeta4

    kl1=omea*t-ome*tt/c
    kl2=omea*zeta4+ome*tt

    epsa=e1a-e3a*ome*t-e3*kl1+(e2a-e1a)*ome*zeta4+(e2-e1)*kl2

    f1aa=2._real_8*t1*((-1._real_8/4._real_8*v1/rs/rs12+3._real_8/4._real_8*x1/rs12+2._real_8*y1)&
         *rsa**2+f1a/(2._real_8*T1*rsa)*rsaa)
    e1aa=-2._real_8*t1*(u1*rsaa*LOG(g1)-2._real_8*u1*rsa/g1/f12*f1a-(1._real_8+u1*rs&
         )*((f1a/f12/g1)**2+1._real_8/g1*(-2._real_8/f1/f12*f1a**2+f1aa/f12)))

    f2aa=2._real_8*t2*((-1._real_8/4._real_8*v2/rs/rs12+3._real_8/4._real_8*x2/rs12+2._real_8*y2)&
         *rsa**2+f2a/(2._real_8*T2*rsa)*rsaa)
    e2aa=-2._real_8*t2*(u2*rsaa*LOG(g2)-2._real_8*u2*rsa/g2/f22*f2a-(1._real_8+u2*rs&
         )*((f2a/f22/g2)**2+(-2._real_8/g2/f2/f22*f2a**2+f2aa/g2/f22)))

    f3aa=2._real_8*t3*((-1._real_8/4._real_8*v3/rs/rs12+3._real_8/4._real_8*x3/rs12+2._real_8*y3)&
         *rsa**2+f3a/(2._real_8*T3*rsa)*rsaa)
    e3aa=-2._real_8*t3*(u3*rsaa*LOG(g3)-2._real_8*u3*rsa/g3/f32*f3a-(1._real_8+u3*rs&
         )*((f3a/f32/g3)**2+(-2._real_8/g3/f3/f32*f3a**2+f3aa/g3/f32)))

    kl1a=omeaa*t-2._real_8*omea*tt/c-ome/c*(12._real_8*zeta2*zetaa**2+4._real_8*zeta3&
         *zetaaa)
    kl2a=omeaa*zeta4+2._real_8*omea*tt+ome*(12._real_8*zeta2*zetaa**2+4._real_8*zeta3&
         *zetaaa)
    epsaa=e1aa-e3aa*ome*t-e3a*(omea*t-ome*tt/c)-e3a*kl1&
         -e3*kl1a+(e2aa-e1aa)*ome*zeta4+(e2a-e1a)*(omea*zeta4+ome*tt)&
         +(e2a-e1a)*kl2+(e2-e1)*kl2a

    nana=2._real_8*epsa+dens*epsaa

    zetab=-densao/densbo*zetaa
    omeb=-densao/densbo*omea
    epsb=e1a-e3a*ome*(1-zeta4)/c-e3*(omeb*(1-zeta4)/c-ome*4._real_8*zeta3*&
         zetab/c)+(e2a-e1a)*ome*zeta4+(e2-e1)*&
         (omeb*zeta4+ome*4._real_8*zeta3*zetab)

    ttb=4._real_8*zeta3*zetab
    kl1b=omeab*t-omea*ttb/c-omeb*tt/c-ome/c*(12._real_8*zeta2*zetaa*zetab+&
         4._real_8*zeta3*zetaab)
    kl2b=omeab*zeta4+omea*ttb+omeb*tt+ome*(12._real_8*zeta2*zetaa*zetab+&
         4._real_8*zeta3*zetaab)
    epsab=e1aa-e3aa*ome*t-e3a*(omeb*t-ome*ttb/c)-e3a*kl1&
         -e3*kl1b+(e2aa-e1aa)*ome*zeta4+(e2a-e1a)*(omeb*zeta4+ome*ttb)&
         +(e2a-e1a)*kl2+(e2-e1)*kl2b

    nanb=epsa+epsb+dens*epsab

    ! -----------------------------------
    ! non local part
    ! -----------------------------------

    u=0.5_real_8*(zp23+zm23)
    uh2=u*u
    uh3=uh2*u
    uh4=uh3*u
    uh5=uh3*uh2
    ua=(o3/zp13-o3/zm13)*zetaa
    ub=(o3/zp13-o3/zm13)*zetab
    d=sgrad/(4._real_8*u)*k1/n_76
    d2=d*d
    d3=d2*d
    d4=d2*d2
    da=sgrad/4._real_8*k1*(-1._real_8/uh2*ua/n_76-7._real_8/6._real_8/u/n_136)
    db=sgrad/4._real_8*k1*(-1._real_8/uh2*ub/n_76-7._real_8/6._real_8/u/n_136)

    ex=EXP(-2._real_8*iota*eps/uh3/lambda**2)
    exa=-1._real_8/k2*(epsa/uh3-eps*3._real_8*ua/uh4)
    exb=-1._real_8/k2*(epsb/uh3-eps*3._real_8*ub/uh4)
    a=k4/(ex-1._real_8)
    a2=a*a
    Aa=-k4*1._real_8/(ex-1._real_8)**2*ex*exa
    ab=-k4*1._real_8/(ex-1._real_8)**2*ex*exb

    ze1=d2+a*d4
    ze1a=2._real_8*d*da+Aa*d4+a*4._real_8*d3*da
    ze1b=2._real_8*d*db+ab*d4+a*4._real_8*d3*db
    ne1=1._real_8+a*d2+a2*d4
    ne1a=Aa*d2+a*2._real_8*d*da+2._real_8*a*aa*d4+a2*4._real_8*d3*da
    ne1b=ab*d2+a*2._real_8*d*db+2._real_8*a*ab*d4+a2*4._real_8*d3*db
    ne12=ne1*ne1
    b=1._real_8+k4*ze1/ne1
    b2=b*b
    Ba=k4*(ze1a/ne1-ze1*ne1a/ne12)
    Bb=k4*(ze1b/ne1-ze1*ne1b/ne12)

    l=uh3*k2*LOG(b)
    La=k2*(3._real_8*uh2*ua*LOG(b)+uh3/b*Ba)
    Lb=k2*(3._real_8*uh2*ub*LOG(b)+uh3/b*Bb)

    uaa=zetaaa*ua/zetaa+0.5_real_8*zetaa**2*&
         (-2._real_8/9._real_8/zp43-2._real_8/9._real_8/zm43)
    uab=zetaab*ua/zetaa+0.5_real_8*zetaa*zetab*&
         (-2._real_8/9._real_8/zp43-2._real_8/9._real_8/zm43)

    daa=sgrad/4._real_8*k1*(2._real_8/uh3*ua**2/n_76-1._real_8/uh2*(uaa/n_76-&
         7._real_8/6._real_8*ua/n_136)-7._real_8/6._real_8*(-1._real_8/uh2*ua/n_136-13._real_8/&
         6._real_8/u/n_196))
    dab=sgrad/4._real_8*k1*(2._real_8/uh3*ua*ub/n_76-1._real_8/uh2*(uab/n_76-&
         7._real_8/6._real_8*ua/n_136)-7._real_8/6._real_8*(-1._real_8/uh2*ub/n_136-13._real_8/&
         6._real_8/u/n_196))

    exaa=-1._real_8/k2*(epsaa/uh3-epsa*6._real_8*ua/uh4-&
         3._real_8*eps*(uaa/uh4-ua**2*4._real_8/uh5))
    exab=-1._real_8/k2*(epsab/uh3-epsa*3._real_8*ub/uh4-&
         epsb*3._real_8*ua/uh4-3._real_8*eps*(uab/uh4-ua*ub*4._real_8/uh5))

    Aaa=-k4*(-2._real_8/(ex-1._real_8)**3*ex**2*exa**2+&
         1._real_8/(ex-1._real_8)**2*(ex*exa**2+ex*exaa))
    Aab=-k4*(-2._real_8/(ex-1._real_8)**3*ex**2*exa*exb+&
         1._real_8/(ex-1._real_8)**2*(ex*exa*exb+ex*exab))

    ze1aa=2._real_8*(da**2+d*daa)+Aaa*d4+Aa*8._real_8*d3*da+&
         4._real_8*A*(3._real_8*d2*da**2+d3*daa)

    ze1ab=2._real_8*(db*da+d*dab)+Aab*d4+Aa*4._real_8*d3*db+4._real_8*ab*d3*da+&
         4._real_8*A*(3._real_8*d2*db*da+d3*dab)

    ne1aa=Aaa*d2+Aa*2._real_8*d*da+2._real_8*aa*d*da+2._real_8*a*(da**2+d*daa)+&
         2._real_8*Aa*Aa*d4+2._real_8*A*(Aaa*d4+Aa*4._real_8*d3*da)&
         +8._real_8*A*Aa*d3*da+4._real_8*A2*(3._real_8*d2*da**2+d3*daa)

    ne1ab=Aab*d2+Aa*2._real_8*d*db+2._real_8*ab*d*da+2._real_8*a*(db*da+d*dab)+&
         2._real_8*Ab*Aa*d4+2._real_8*A*(Aab*d4+Aa*4._real_8*d3*db)&
         +8._real_8*A*Ab*d3*da+4._real_8*A2*(3._real_8*d2*db*da+d3*dab)

    Baa=k4*(ze1aa/ne1-ze1a*ne1a/ne12-(ze1a*ne1a+ze1*&
         ne1aa)/ne12+ze1*ne1a**2*2._real_8*ne1/ne12/ne12)

    Bab=k4*(ze1ab/ne1-ze1a*ne1b/ne12-(ze1b*ne1a+ze1*&
         ne1ab)/ne12+ze1*ne1a*ne1b*2._real_8*ne1/ne12/ne12)

    Laa=k2*(6._real_8*uh2*ua/b*Ba+LOG(b)*(6._real_8*u*ua**2+&
         3._real_8*uh2*uaa)+uh3*(-Ba**2/B2+Baa/B))

    Lab=k2*(3._real_8*uh2*ua/b*Bb+LOG(b)*(6._real_8*u*ub*ua+&
         3._real_8*uh2*uab)+3._real_8*uh2*ub/B*Ba+uh3*(-Ba/B2*Bb+Bab/B))

    ze1ga=(d2+a*2._real_8*d4)/grad

    ne1ga=(a*d2+a2*2._real_8*d4)/grad
    Bga=k4*(ze1ga/ne1-ze1*ne1ga/ne12)


    Lga=k2*uh3*Bga/b

    ze1gaga=a*2._real_8*d4/grad2
    ne1gaga=a2*2._real_8*d4/grad2

    Bgaga=k4*(ze1gaga/ne1-ze1ga*ne1ga/ne12-(ze1ga*ne1ga+ze1*&
         ne1gaga)/ne12+ze1*ne1ga**2*2._real_8*ne1/ne12/ne12)

    Lgaga=k2*uh3*(-Bga**2/b2+Bgaga/b)

    ze1naga=(2._real_8*d*da+Aa*2._real_8*d4+8._real_8*a*d3*da)/grad
    ze1ganb=(2._real_8*d*db+ab*2._real_8*d4+8._real_8*a*d3*db)/grad
    ne1naga=(Aa*d2+a*2._real_8*d*da+4._real_8*a*aa*d4+8._real_8*a2*d3*da)/grad
    ne1ganb=(ab*d2+a*2._real_8*d*db+4._real_8*a*ab*d4+8._real_8*a2*d3*db)/grad

    Bnaga=k4*(ze1naga/ne1-ze1ga*ne1a/ne12-(ze1a*ne1ga+ze1*&
         ne1naga)/ne12+ze1*ne1ga*ne1a*2._real_8*ne1/ne12/ne12)
    Bganb=k4*(ze1ganb/ne1-ze1ga*ne1b/ne12-(ze1b*ne1ga+ze1*&
         ne1ganb)/ne12+ze1*ne1ga*ne1b*2._real_8*ne1/ne12/ne12)

    Lnaga=k2*(3._real_8*uh2*ua*Bga/b+uh3*(-bga/b2*Ba+Bnaga/b))
    Lganb=k2*(3._real_8*uh2*ub*Bga/b+uh3*(-bga/b2*Bb+Bganb/b))

    Ja=0._real_8
    Jaa=Ja
    Jab=Ja
    Jga=Ja
    Jnaga=Ja
    Jganb=Ja
    Jgaga=Ja
    Jb=Ja

    nana=nana+2._real_8*(La+Ja)+dens*(Laa+Jaa)
    nanb=nanb+Lb+Jb+Ja+La+dens*(Lab+Jab)
    ga=dens*(Lga+Jga)
    gaga=dens*(Lgaga+Jgaga)
    naga=Lga+Jga+dens*(Lnaga+Jnaga)
    ganb=Lga+Jga+dens*(Lganb+Jganb)

    nask=2._real_8*naga
    gask=2._real_8*gaga
    sk=2._real_8*ga
    sksk=4._real_8*gaga
    gagb=gaga
    nagb=naga
    skgb=2._real_8*gaga
    sknb=2._real_8*ganb
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE pbe96_c
  ! ==================================================================
  SUBROUTINE pbe96_c_loc(densa,densb,nana,nanb)
    ! ==================================================================
    REAL(real_8)                             :: densa, densb, nana, nanb

    REAL(real_8), PARAMETER :: c = 1.709921_real_8, &
      fakt = 0.6203504908994000_real_8, o = 6.46309613581743018602522_real_8, &
      o3 = 0.3333333333333333333333_real_8 , t1 = 0.031091_real_8, &
      t2 = 0.015545_real_8, T3 = 0.016887_real_8, U1 = 0.2137_real_8, &
      U2 = 0.20548_real_8, U3 = 0.11125_real_8, V1 = 7.5957_real_8, &
      V2 = 14.1189_real_8, V3 = 10.357_real_8, W1 = 3.5876_real_8, &
      W2 = 6.1977_real_8, W3 = 3.6231_real_8, X1 = 1.6382_real_8, &
      X2 = 3.3662_real_8, X3 = 0.88026_real_8, Y1 = 0.49294_real_8, &
      Y2 = 0.62517_real_8, Y3 = 0.49671_real_8

    REAL(real_8) :: dens, densao, densbo, e1, e1a, e1aa, e2, e2a, e2aa, e3, &
      e3a, e3aa, eps, epsa, epsaa, epsab, epsb, f1, f12, f1a, f1aa, f2, f22, &
      f2a, f2aa, f3, f32, f3a, f3aa, g1, g2, g3, kl1, kl1a, kl1b, kl2, kl2a, &
      kl2b, n2, n3, n_103, n_13, n_43, n_73, na_13, na_23, nb_13, nb_23, ome, &
      omea, omeaa, omeab, omeb, rs, rs12, rs2, rs3, rsa, rsaa, t, tt, ttb, &
      zeta, zeta2, zeta3, zeta4, zetaa, zetaaa, zetaab, zetab, zm43, zp43

! ==--------------------------------------------------------------==

    densao=MAX(densa,1.e-24_real_8)
    densbo=MAX(densb,1.e-24_real_8)

    dens=densao+densbo
    n2=dens*dens
    n3=n2*dens
    n_13=dens**o3
    n_43=n_13*dens
    n_73=n_43*dens
    n_103=n_73*dens
    na_13=densao**o3
    na_23=na_13*na_13
    nb_13=densbo**o3
    nb_23=nb_13*nb_13

    rs=fakt/n_13
    rs2=rs*rs
    rs3=rs2*rs
    rs12=SQRT(rs)
    rsa=-o3*fakt/n_43
    rsaa=4._real_8/9._real_8*fakt/n_73
    zeta=(densao-densbo)/dens
    zeta2=zeta*zeta
    zeta3=zeta2*zeta
    zeta4=zeta3*zeta
    zetaa=2._real_8*densbo/n2
    zp43=(1._real_8+zeta)**(4._real_8/3._real_8)
    zm43=(1._real_8-zeta)**(4._real_8/3._real_8)
    ome=(zp43+zm43-2._real_8)/(2._real_8**(4._real_8/3._real_8)-2._real_8)
    omea=o*densbo*(na_13-nb_13)/n_73

    f1=2._real_8*t1*(v1*rs12+w1*rs+x1*rs*rs12+y1*rs2)
    f12=f1*f1
    g1=(f1+1._real_8)/f1
    f1a=2._real_8*t1*(0.5_real_8*v1/rs12+w1+3._real_8/2._real_8*x1*rs12+2._real_8*y1*rs)*rsa
    e1=-2._real_8*t1*(1._real_8+u1*rs)*LOG(g1)
    e1a=-2._real_8*t1*(u1*rsa*LOG(g1)-(1._real_8+u1*rs)/g1/f12*f1a)

    f2=2._real_8*t2*(v2*rs12+w2*rs+x2*rs*rs12+y2*rs2)
    f22=f2*f2
    g2=(f2+1._real_8)/f2
    f2a=2._real_8*t2*(0.5_real_8*v2/rs12+w2+3._real_8/2._real_8*x2*rs12+2._real_8*y2*rs)*rsa
    e2=-2._real_8*t2*(1._real_8+u2*rs)*LOG(g2)
    e2a=-2._real_8*t2*(u2*rsa*LOG(g2)-(1._real_8+u2*rs)/g2/f22*f2a)

    f3=2._real_8*t3*(v3*rs12+w3*rs+x3*rs*rs12+y3*rs2)
    f32=f3*f3
    g3=(f3+1._real_8)/f3
    f3a=2._real_8*t3*(0.5_real_8*v3/rs12+w3+3._real_8/2._real_8*x3*rs12+2._real_8*y3*rs)*rsa
    e3=-2._real_8*t3*(1._real_8+u3*rs)*LOG(g3)
    e3a=-2._real_8*t3*(u3*rsa*LOG(g3)-(1._real_8+u3*rs)/g3/f32*f3a)

    zetaaa=-4._real_8*densbo/n3
    zetaab=2._real_8*zeta/n2
    omeaa=o*densbo*(o3/na_23/n_73-(na_13-nb_13)*7._real_8/3._real_8/n_103)
    omeab=o*((na_13-nb_13)/n_73+densbo*(-o3/nb_23/n_73-(na_13-nb_13)&
         *7._real_8/3._real_8/n_103))
    t=(1._real_8-zeta4)/c
    tt=4._real_8*zeta3*zetaa

    eps=e1-e3*ome*t+(e2-e1)*ome*zeta4

    kl1=omea*t-ome*tt/c
    kl2=omea*zeta4+ome*tt

    epsa=e1a-e3a*ome*t-e3*kl1+(e2a-e1a)*ome*zeta4+(e2-e1)*kl2

    f1aa=2._real_8*t1*((-1._real_8/4._real_8*v1/rs/rs12+3._real_8/4._real_8*x1/rs12+2._real_8*y1)&
         *rsa**2+f1a/(2._real_8*T1*rsa)*rsaa)
    e1aa=-2._real_8*t1*(u1*rsaa*LOG(g1)-2._real_8*u1*rsa/g1/f12*f1a-(1._real_8+u1*rs&
         )*((f1a/f12/g1)**2+1._real_8/g1*(-2._real_8/f1/f12*f1a**2+f1aa/f12)))

    f2aa=2._real_8*t2*((-1._real_8/4._real_8*v2/rs/rs12+3._real_8/4._real_8*x2/rs12+2._real_8*y2)&
         *rsa**2+f2a/(2._real_8*T2*rsa)*rsaa)
    e2aa=-2._real_8*t2*(u2*rsaa*LOG(g2)-2._real_8*u2*rsa/g2/f22*f2a-(1._real_8+u2*rs&
         )*((f2a/f22/g2)**2+(-2._real_8/g2/f2/f22*f2a**2+f2aa/g2/f22)))

    f3aa=2._real_8*t3*((-1._real_8/4._real_8*v3/rs/rs12+3._real_8/4._real_8*x3/rs12+2._real_8*y3)&
         *rsa**2+f3a/(2._real_8*T3*rsa)*rsaa)
    e3aa=-2._real_8*t3*(u3*rsaa*LOG(g3)-2._real_8*u3*rsa/g3/f32*f3a-(1._real_8+u3*rs&
         )*((f3a/f32/g3)**2+(-2._real_8/g3/f3/f32*f3a**2+f3aa/g3/f32)))

    kl1a=omeaa*t-2._real_8*omea*tt/c-ome/c*(12._real_8*zeta2*zetaa**2+4._real_8*zeta3&
         *zetaaa)
    kl2a=omeaa*zeta4+2._real_8*omea*tt+ome*(12._real_8*zeta2*zetaa**2+4._real_8*zeta3&
         *zetaaa)
    epsaa=e1aa-e3aa*ome*t-e3a*(omea*t-ome*tt/c)-e3a*kl1&
         -e3*kl1a+(e2aa-e1aa)*ome*zeta4+(e2a-e1a)*(omea*zeta4+ome*tt)&
         +(e2a-e1a)*kl2+(e2-e1)*kl2a

    nana=2._real_8*epsa+dens*epsaa

    zetab=-densao/densbo*zetaa
    omeb=-densao/densbo*omea
    epsb=e1a-e3a*ome*(1-zeta4)/c-e3*(omeb*(1-zeta4)/c-ome*4._real_8*zeta3*&
         zetab/c)+(e2a-e1a)*ome*zeta4+(e2-e1)*&
         (omeb*zeta4+ome*4._real_8*zeta3*zetab)

    ttb=4._real_8*zeta3*zetab
    kl1b=omeab*t-omea*ttb/c-omeb*tt/c-ome/c*(12._real_8*zeta2*zetaa*zetab+&
         4._real_8*zeta3*zetaab)
    kl2b=omeab*zeta4+omea*ttb+omeb*tt+ome*(12._real_8*zeta2*zetaa*zetab+&
         4._real_8*zeta3*zetaab)
    epsab=e1aa-e3aa*ome*t-e3a*(omeb*t-ome*ttb/c)-e3a*kl1&
         -e3*kl1b+(e2aa-e1aa)*ome*zeta4+(e2a-e1a)*(omeb*zeta4+ome*ttb)&
         +(e2a-e1a)*kl2+(e2-e1)*kl2b

    nanb=epsa+epsb+dens*epsab
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE pbe96_c_loc
  ! ==================================================================
  SUBROUTINE pbe96_x(dens,grad,dfdg,d2fdn2,d2fdndg,d2fdg2)
    ! ==================================================================
    REAL(real_8)                             :: dens, grad, dfdg, d2fdn2, &
                                                d2fdndg, d2fdg2

    REAL(real_8), PARAMETER :: a = 0.16162045967399548133_real_8, &
      b = 0.98474502184269654115_real_8 , m = 0.2195164512208958179_real_8, &
      r = 0.804_real_8

    REAL(real_8) :: a2, dens0, f, Fg, Fgg, Fng, fs, Fss, grad0, kl, kl2, kl3, &
      n_103, n_13, n_23, n_43, n_73, n_83, o3, s, S2, sgrad, ss, sss, x, xs, &
      xss

! ==--------------------------------------------------------------==

    dens0=2._real_8*MAX(dens,1.e-24_real_8) ! !!E(2n)
    grad0=4._real_8*MAX(grad,1.e-48_real_8)
    sgrad=SQRT(grad0)

    a2=a*a
    o3=1._real_8/3._real_8
    n_13=dens0**o3
    n_23=n_13*n_13
    n_43=n_13*dens0
    n_73=n_43*dens0
    n_83=n_73*n_13
    n_103=n_73*dens0

    x=sgrad/n_43
    s=x*a

    s2=s*s
    kl=1._real_8+m*s2/r
    kl2=kl*kl
    kl3=kl2*kl
    f=1._real_8+r-r/kl

    xs=-4._real_8/3._real_8*sgrad/n_73
    ss=a*xs
    fs=2._real_8/kl2*m*s*ss

    xss=28._real_8/9._real_8*sgrad/n_103
    sss=a*xss
    fss=-m*(8._real_8/r*m/kl3*s2*ss**2-2._real_8/kl2*(ss**2+s*sss))
    d2fdn2=2._real_8*b*(-o3/n_23*f-2._real_8*n_13*fs-3._real_8/4._real_8*n_43*fss)

    Fg=m/kl2*a2/n_83
    Fgg=-m**2*a2**2/n_83*2._real_8/kl3/r/n_83
    Fng=m/grad0*(2._real_8*s*ss/kl2-s2*4._real_8*m/r*s*ss/kl3)

    d2fdndg=4._real_8*(-b*n_13*Fg-3._real_8/4._real_8*b*n_43*Fng)

    dfdg=-3._real_8/2._real_8*b*n_43*Fg
    d2fdg2=-6._real_8*b*n_43*Fgg
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE pbe96_x
  ! ==================================================================





END MODULE dd_functionals_utils
