! ==================================================================
! Provides: - GGA correlation derivatives for linear response
!
! Input:  Spin-polarised densities from dd_xc_ana_utils
! Output: tmp1...5 analogous to dd_xc_ana_utils
!
!                             02.10.2016 - M. P. Bircher @ LCBC/EPFL
! ==================================================================
MODULE cp_dgga_correlation_utils
  USE cpfunc_types,                    ONLY: cp_dxc_functional_t,&
                                             cp_dxc_scratch_t,&
                                             cp_xc_spin_pairs
  USE func,                            ONLY: func2
  USE kinds,                           ONLY: real_8
  USE system,                          ONLY: cntr

  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER      :: alpha       = 1
  INTEGER, PARAMETER      :: beta        = 2

  REAL(real_8), PARAMETER :: four_thirds = 4._real_8/3._real_8

  PUBLIC :: CP_dGGA_C_LYP
  PUBLIC :: CP_dGGA_C_P86
  PUBLIC :: CP_dGGA_C_PBE
  PUBLIC :: CP_SPIN_dGGA_C_LYP
  PUBLIC :: CP_SPIN_dGGA_C_P86
  PUBLIC :: CP_SPIN_dGGA_C_PBE

CONTAINS

  ! ==================================================================
  PURE SUBROUTINE cp_spin_dgga_c_lyp( scratch,functional )
    ! ==--------------------------------------------------------------==
    ! Open-shell wrapper to actual LYP derivatives
    !                                     19.06.2017 M.P. Bircher @ EPFL
    ! ==--------------------------------------------------------------==

    TYPE(cp_dxc_scratch_t), &
      DIMENSION(cp_xc_spin_pairs), &
      INTENT(in)                             :: scratch
    TYPE(cp_dxc_functional_t), &
      DIMENSION(cp_xc_spin_pairs), &
      INTENT(inout)                          :: functional

    CALL cp_spin_dgga_c_lyp_compute( scratch(alpha), scratch(beta), functional(alpha) )
    CALL cp_spin_dgga_c_lyp_compute( scratch(beta), scratch(alpha), functional(beta ) )

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_spin_dgga_c_lyp
  ! ==================================================================
  PURE SUBROUTINE cp_dgga_c_lyp( scratch,functional )
    ! ==--------------------------------------------------------------==
    ! Closed-shell wrapper to actual LYP derivatives
    !                                     19.06.2017 M.P. Bircher @ EPFL
    ! ==--------------------------------------------------------------==

    TYPE(cp_dxc_scratch_t), INTENT(in)       :: scratch
    TYPE(cp_dxc_functional_t), INTENT(inout) :: functional

    CALL cp_spin_dgga_c_lyp_compute( scratch, scratch, functional )

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_dgga_c_lyp
  ! ==================================================================
  PURE SUBROUTINE cp_spin_dgga_c_p86( scratch,functional )
    ! ==--------------------------------------------------------------==
    ! Open-shell wrapper to actual P86 derivatives
    !                                     20.09.2016 M.P. Bircher @ EPFL
    ! ==--------------------------------------------------------------==

    TYPE(cp_dxc_scratch_t), &
      DIMENSION(cp_xc_spin_pairs), &
      INTENT(in)                             :: scratch
    TYPE(cp_dxc_functional_t), &
      DIMENSION(cp_xc_spin_pairs), &
      INTENT(inout)                          :: functional

    CALL cp_spin_dgga_c_p86_compute( scratch(alpha), scratch(beta), functional(alpha) )
    CALL cp_spin_dgga_c_p86_compute( scratch(beta), scratch(alpha), functional(beta ) )

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_spin_dgga_c_p86
  ! ==================================================================
  PURE SUBROUTINE cp_dgga_c_p86( scratch,functional )
    ! ==--------------------------------------------------------------==
    ! Closed-shell wrapper to actual P86 derivatives
    !                                     20.06.2017 M.P. Bircher @ EPFL
    ! ==--------------------------------------------------------------==

    TYPE(cp_dxc_scratch_t), INTENT(in)       :: scratch
    TYPE(cp_dxc_functional_t), INTENT(inout) :: functional

    CALL cp_spin_dgga_c_p86_compute( scratch, scratch, functional )

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_dgga_c_p86
  ! ==================================================================
  PURE SUBROUTINE cp_spin_dgga_c_pbe( scratch,functional )
    ! ==--------------------------------------------------------------==
    ! Open-shell wrapper to actual PBE derivatives
    !                                     19.06.2017 M.P. Bircher @ EPFL
    ! ==--------------------------------------------------------------==

    TYPE(cp_dxc_scratch_t), &
      DIMENSION(cp_xc_spin_pairs), &
      INTENT(in)                             :: scratch
    TYPE(cp_dxc_functional_t), &
      DIMENSION(cp_xc_spin_pairs), &
      INTENT(inout)                          :: functional

    CALL cp_spin_dgga_c_pbe_compute( scratch(alpha), scratch(beta), functional(alpha) )
    CALL cp_spin_dgga_c_pbe_compute( scratch(beta), scratch(alpha), functional(beta ) )

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_spin_dgga_c_pbe
  ! ==================================================================
  PURE SUBROUTINE cp_dgga_c_pbe( scratch,functional )
    ! ==--------------------------------------------------------------==
    ! Closed-shell wrapper to actual PBE derivatives
    !                                     19.06.2017 M.P. Bircher @ EPFL
    ! ==--------------------------------------------------------------==

    TYPE(cp_dxc_scratch_t), INTENT(in)       :: scratch
    TYPE(cp_dxc_functional_t), INTENT(inout) :: functional

    CALL cp_spin_dgga_c_pbe_compute( scratch, scratch, functional )

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_dgga_c_pbe
  ! ==================================================================

  ! ==================================================================
  ! Modular routines accessed by the wrappers
  ! ==================================================================

  ! ==================================================================
  ELEMENTAL SUBROUTINE cp_spin_dgga_c_lyp_compute( scratch_alpha, scratch_beta, functional )
    ! ==--------------------------------------------------------------==
    ! Ported from "old" xc_ana code, not optimised or prettified
    !                                     19.06.2017 M.P. Bircher @ EPFL
    ! ==--------------------------------------------------------------==

    TYPE(cp_dxc_scratch_t), INTENT(in)       :: scratch_alpha, scratch_beta
    TYPE(cp_dxc_functional_t), INTENT(inout) :: functional

    REAL(real_8), PARAMETER :: a = 0.04918_real_8 , b = 0.132_real_8 , &
      balpha = -0.930525736349100025002010218071667_real_8 , &
      c = 0.2533_real_8 , cf = 36.462398978764777098305310076_real_8 , &
      d = 0.349_real_8 

    REAL(real_8) :: nana, nanb, naga, granb, gagb, nagb, nagrb, ganb
    REAL(real_8) :: ab, br, c13, D1, d13, ddens, dens, densao, densao2, &
      densbo, densbo2, dn2, e, e3, es, ess, f3, grad, gradao, gradbo, kl, &
      n113, n143, n2, n3, n43, n_113, n_13, n_143, n_173, n_43, n_73, na_23, &
      na_53, na_83, nb_53, nb_83, o18, o3, o9, s18, skalar, t1, t12, t1s, &
      t1ss, t2, t2s, t2ss, t2sst, t2st, t3, t3s, t3ss, t4, t4p, t4s, t4sp, &
      t4ss, t4sst, t4st, t5, t5k, t5p, t5pst, t5s, t5sk, t5sp, t5ss, t5sst, &
      t5st, t6, t6s, t6ss, t7, t7p, t7s, t7ss, t7sst, t7st, t91, t92, tw3, &
      skaa, skab, skba, skbb, sgrad, dfdga, loca, locb, na_113, nb_113

    densao = MAX(scratch_alpha%n,1.e-24_real_8)
    densbo = MAX(scratch_beta%n, 1.e-24_real_8)
    gradao = MAX(scratch_alpha%g_2,1.e-48_real_8)
    gradbo = MAX(scratch_beta%g_2, 1.e-48_real_8)
    dens   = densao+densbo
    skalar = scratch_alpha%gx_2*scratch_beta%gx_2 + &
             scratch_alpha%gy_2*scratch_beta%gy_2 + &
             scratch_alpha%gz_2*scratch_beta%gz_2 
    grad = gradao + gradbo + 2.0_real_8*skalar
    n2   = dens*dens
    n3   = 1._real_8/n2/dens
    densao2 = densao*densao
    densbo2 = densbo*densbo
    o3  = 1._real_8/3._real_8
    tw3 = o3+o3
    f3  = tw3+tw3
    e3  = 1._real_8+f3+f3
    o9  = o3*o3
    o18 = o9/2._real_8
    s18 = o3+o18
    n_13 = dens**o3
    n_43 = n_13*dens
    n_73 = 1._real_8/n_43/dens
    n_113 = n_43*n_43*dens
    n_143 = n_113*dens
    n_173 = n_143*dens
    n43 = 1._real_8/n_43
    n113 = 1._real_8/n_113
    n143 = 1._real_8/n_143
    na_23 = densao**tw3
    na_53 = na_23*densao
    na_83 = na_53*densao
    na_113= na_83*densao
    nb_53 = densbo**(5._real_8/3._real_8)
    nb_83 = nb_53*densbo
    nb_113= nb_83*densbo
    c13 = c/n_13
    e = EXP(-c13)
    d13 = d/n_13
    t1 = 1._real_8/(1._real_8+d13)
    t12 = t1*t1
    ab = densao*densbo
    ddens = 1._real_8/dens
    t2 = ab*ddens
    t3 = a*b*e*t1*n113
    t7 = (densao*gradao+densbo*gradbo)*ddens
    t6 = c13+d13*t1
    br = o9*(t6-11._real_8)
    kl = (47._real_8*o18-s18*t6)
    t5 = Cf*(na_83+nb_83)+kl*grad-(2.5_real_8-o18*t6)*(gradao+gradbo)-br*t7
    t4 = ab*t5+tw3*n2*(gradao+gradbo-grad)-densao2*gradbo-densbo2*gradao
    es = c*e*o3*n43
    t1s = o3*d*t12*n43
    dn2 = 1._real_8/n2
    t2s = densbo2*dn2
    t2st = densao2*dn2
    t3s = a*b*(t1s*e*n113+t1*es*n113-t1*e*e3*n143)
    t7s = densbo*dn2*(gradao-gradbo)
    t7st = densao*dn2*(gradbo-gradao)
    t6s = -o3*n43*(c+d*t1)+d13*t1s
    t91 = -s18*grad*t6s-o9*t6s*t7+o18*t6s*(gradao+gradbo)
    t5s = Cf*8._real_8*o3*na_53-br*t7s+t91
    t5st = Cf*8._real_8*o3*nb_53-br*t7st+t91
    t4s = densbo*t5+ab*t5s+(f3*dens-2._real_8*densao)*gradbo&
         +f3*dens*(gradao-grad)
    t4st = densao*t5+ab*t5st+(f3*dens-2._real_8*densbo)&
         *gradao+f3*dens*(gradbo-grad)
    d1 = -4._real_8*a*(t1s*t2+t1*t2s)-t3s*t4-t3*t4s
    t7p = densao*ddens
    t5p = o9-o3*t6-br*t7p
    t5k = o9-o3*t6-br*densbo*ddens
    t5sp = -o3*t6s*(1+o3*t7p)-br*densbo*dn2
    t5sk = -o3*t6s*(1+o3*densbo*ddens)+br*densbo*dn2
    t4p = ab*t5p-densbo2
    t4sp = densbo*t5p+ab*t5sp
    dfdga = -t3*t4p
    naga = -t3s*t4p-t3*t4sp
    ess = c*o9*e*n_73*(c13-4._real_8)
    t1ss = o3*d*(2._real_8*t1*t1s*n43-f3*n_73*t12)
    t2ss = -2._real_8*densbo2*n3
    t2sst = 2._real_8*ab*n3
    t3ss = a*b*((t1ss*e+2._real_8*t1s*es+t1*ess)*n113-2._real_8*e3*(t1s*e&
         *n143+t1*es*n143-7._real_8*o3*t1*e/n_173))
    t7ss = -2._real_8*densbo*n3*(gradao-gradbo)
    t7sst = (densao-densbo)*n3*(gradao-gradbo)
    t6ss = 4._real_8*o9*n_73*(c+d*t1)-tw3*d*n43*t1s+d13*t1ss
    t92 = t6ss*(o18*(gradao+gradbo)-o9*t7-s18*grad)
    t5ss = Cf*o9*40._real_8*na_23-2._real_8*o9*t6s*t7s-br*t7ss+t92
    t5sst = -o9*t6s/n2*(gradao-gradbo)*(densbo-densao)-br*t7sst+t92
    t4ss = 2._real_8*densbo*t5s+ab*t5ss-f3*(grad+0.5_real_8*gradbo-gradao)
    t4sst = densbo*t5st+t5+densao*t5s+ab*t5sst+f3*(-grad+gradbo+gradao)
    nana = -4._real_8*a*(t1ss*t2+2._real_8*t1s*t2s+t1*t2ss)-t3ss*t4&
         -2._real_8*t3s*t4s-t3*t4ss
    t5pst = -o3*t6s-densao*(o9*t6s*ddens-br*dn2)
    granb = -t3s*t4p-t3*(densao*t5p+ab*t5pst-2._real_8*densbo)
    nanb = -4._real_8*a*(t1ss*t2+t1s*t2st+t1s*t2s+t1*t2sst)&
         -t3ss*t4-t3s*t4st-t3s*t4s-t3*t4sst
    gagb = -t3*(ab*2._real_8*kl-f3*n2)
    ganb = 2._real_8*(-t3s*(ab*kl-tw3*n2)-t3*(densao*kl-ab*s18*t6s-f3*dens))
    nagrb = -t3s*(ab*t5k-densao2)-t3*(densbo*t5k+ab*t5sk-2._real_8*densao)
    nagb = 2._real_8*(-t3s*(ab*kl-tw3*n2)-t3*(-f3*dens+densbo*kl-ab*s18*t6s))

    IF (scratch_alpha%n > 1.e-24_real_8 .AND. (scratch_alpha%n < cntr%gceps .OR. &
                                                scratch_alpha%g_2 <= 1.e-48_real_8)) THEN
       loca = -4._real_8*a*(t1ss*t2+2._real_8*t1s*t2s+t1*t2ss)-Cf*t3ss*(na_113*densbo&
            +densao*nb_113)-2._real_8*Cf*t3s*(11._real_8/3._real_8*na_83*densbo+nb_113)&
            -Cf*t3*88._real_8/9._real_8*na_53*densbo
       locb = -4._real_8*a*(t1ss*t2+t1s*t2st+t1s*t2s+t1*t2sst)&
            -Cf*t3ss*(na_113*densbo+densao*nb_113)-cf*t3s*(na_113+nb_113&
            +e3*(densao*nb_83+densbo*na_83))-Cf*t3*e3*(na_83+nb_83)
       loca = loca+4._real_8*o9*balpha/na_23
       functional%tmp1 = loca*scratch_alpha%dn+locb*scratch_beta%dn
    ELSEIF (scratch_alpha%n > 1.e-24_real_8) THEN
       skaa = scratch_alpha%gx_2*scratch_alpha%dgx+scratch_alpha%gy_2*scratch_alpha%dgy&
             +scratch_alpha%gz_2*scratch_alpha%dgz
       skab = scratch_alpha%gx_2*scratch_beta%dgx+scratch_alpha%gy_2*scratch_beta%dgy&
             +scratch_alpha%gz_2*scratch_beta%dgz
       skba = scratch_beta%gx_2*scratch_alpha%dgx+scratch_beta%gy_2*scratch_alpha%dgy&
             +scratch_beta%gz_2*scratch_alpha%dgz
       skbb = scratch_beta%gx_2*scratch_beta%dgx+scratch_beta%gy_2*scratch_beta%dgy&
             +scratch_beta%gz_2*scratch_beta%dgz
       sgrad  = SQRT(scratch_alpha%g_2)
       functional%tmp1 = (nana)*scratch_alpha%dn&
                         +nanb*scratch_beta%dn+naga*2._real_8*skaa+nagb*skab&
                         +nagrb*2._real_8*skbb+nagb*skba
       functional%tmp2 = 2._real_8*(naga*scratch_alpha%dn+granb*scratch_beta%dn)
       functional%tmp3 = dfdga*2._real_8
       functional%tmp4 = nagb*scratch_alpha%dn+ganb*scratch_beta%dn
       functional%tmp5 = gagb
    ENDIF

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_spin_dgga_c_lyp_compute
  ! ==================================================================
  ELEMENTAL SUBROUTINE cp_spin_dgga_c_p86_compute( scratch_alpha, scratch_beta, functional )
    ! ==--------------------------------------------------------------==
    ! Ported from "old" xc_ana code, not optimised or prettified
    !                                     20.06.2017 M.P. Bircher @ EPFL
    ! ==--------------------------------------------------------------==

    TYPE(cp_dxc_scratch_t), INTENT(in)       :: scratch_alpha, scratch_beta
    TYPE(cp_dxc_functional_t), INTENT(inout) :: functional

    REAL(real_8), PARAMETER :: af = 0.01555_real_8, ap = 0.0311_real_8, &
      b1f = 1.3981_real_8, b1p = 1.0529_real_8, b2f = 0.2611_real_8 , &
      b2p = 0.3334_real_8, bf = -0.0269_real_8, bp = -0.048_real_8, &
      cf = 0.0007_real_8, cp = 0.002_real_8, df = -0.0048_real_8 , &
      dp = -0.0116_real_8, f3 = 4._real_8/3._real_8 , &
      fakt = 0.6203504908994000_real_8 , fakt2 = 0.519842099789746_real_8 , &
      fakt3 = 0.629960524947436_real_8 , gf = -0.0843_real_8, &
      gp = -0.1423_real_8, k = 8.1290825e-4_real_8, o3 = 1._real_8/3._real_8, &
      p0 = 0.002568_real_8, p1 = 0.023266_real_8, p2 = 7.389e-6_real_8, &
      p3 = 8.723_real_8,  p4 = 0.472_real_8, p5 = 7.389e-2_real_8 

    REAL(real_8) :: nana, naga, gaga, nask, gask, sk, ga, sksk, &
      gagb, nagb, skgb, sknb, ganb, nanb, skaa, skab, skba, skbb

    REAL(real_8) :: c, c2, c3, cs, cs2, css, d, d2, d3, dens, densao, densbo, &
      ds, ds2, dsb, dss, dssb, e, epsf, epsfs, epsfss, epsp, epsps, epspss, &
      f, fs, fsb, fss, fssb, grad, gradao, gradbo, n2, n3, n_103, n_13, &
      n_136, n_196, n_43, n_73, n_76, phi, phip, phipb, phipp, phippb, phis, &
      phisch, phischsch, phisp, phiss, phissch, rs, rs2, rs3, rsq, rss, rsss, &
      sgrad, sgrada, sgradb, skalar, t1, t12, t2, t3, t4, t5, t5s, t5sb, t6f, &
      t6p, z, zeta, zetas, zetas2, zetasb, zetass, zetassb, zm13, zm23, zm43, &
      zm53, zp13, zp23, zp43, zp53, Zs, Zsb, Zss, Zssb, p86a, p86b

    densao = MAX(scratch_alpha%n,1.0e-24_real_8)
    densbo = MAX(scratch_beta%n ,1.0e-24_real_8)
    gradao = MAX(scratch_alpha%g_2,1.0e-48_real_8)
    gradbo = MAX(scratch_beta%g_2 ,1.0e-48_real_8)

    dens = densao+densbo
    n2   = dens*dens
    n3   = n2*dens
    n_13 = dens**o3
    n_43 = n_13*dens
    n_73 = n_43*dens
    n_103 = n_73*dens
    n_76  = SQRT(n_73)
    n_136 = n_76*dens
    n_196 = n_136*dens
    skalar = scratch_alpha%gx_2*scratch_beta%gx_2 + &
             scratch_alpha%gy_2*scratch_beta%gy_2 + &
             scratch_alpha%gz_2*scratch_beta%gz_2 
    grad   = gradao+gradbo+2._real_8*skalar

    zeta = (densao-densbo)/dens
    zp13 = (1._real_8+zeta)**o3
    zm13 = (1._real_8-zeta)**o3
    zp23 = zp13*zp13
    zm23 = zm13*zm13
    zp43 = zp23*zp23
    zm43 = zm23*zm23
    zp53 = zp43*zp13
    zm53 = zm43*zm13
    rs = fakt/n_13
    rs2 = rs*rs
    rs3 = rs2*rs
    rss = -o3*fakt/n_43
    rsss = 4._real_8/9._real_8*fakt/n_73
    f = (zp43+zm43-2._real_8)/fakt2

    ! 
    ! Discerns between high- and low-density (cf. CP_LDA_C_PZ)
    !
    ! Low-density formula
    !
    IF (rs >= 1.0_real_8) THEN
       rsq = SQRT(rs)
       epsp = gp/(1._real_8+b1p*rsq+b2p*rs)
       epsf = gf/(1._real_8+b1f*rsq+b2f*rs)
       t6p = b1p/(2._real_8*rsq)
       t6f = b1f/(2._real_8*rsq)
       epsps = -epsp**2/gp*(t6p+b2p)*rss
       epsfs = -epsf**2/gf*(t6f+b2f)*rss
       epspss = -1._real_8/gp*epsp*(2._real_8*epsps*(t6p+b2p)*rss+epsp*(-t6p/&
            (2._real_8*rs)*rss+(t6p+b2p)*rsss))
       epsfss = -1._real_8/gf*epsf*(2._real_8*epsfs*(t6f+b2f)*rss+epsf*(-t6f/&
            (2._real_8*rs)*rss+(t6f+b2f)*rsss))
    !
    ! High-density formula
    !
    ELSE
       epsp = ap*LOG(rs)+bp+cp*rs*LOG(rs)+dp*rs
       epsf = af*LOG(rs)+bf+cf*rs*LOG(rs)+df*rs
       epsps = (ap/rs+cp*(LOG(rs)+1._real_8)+dp)*rss
       epsfs = (af/rs+cf*(LOG(rs)+1._real_8)+df)*rss
       epspss = (-ap/rs2+cp/rs)*rss*rss+epsps/rss*rsss
       epsfss = (-af/rs2+cf/rs)*rss*rss+epsfs/rss*rsss
    ENDIF

    z = epsp+f*(epsf-epsp)
    t1 = 1._real_8+p3*rs+p4*rs2+p5*rs3
    t12 = t1*t1
    t2 = p0+p1*rs+p2*rs2
    c = 0.001667_real_8+t2/t1
    c2 = c*c
    c3 = c2*c
    d = 2._real_8*fakt3*SQRT(0.5_real_8*fakt3*(zp53+zm53))
    d2 = d*d
    d3 = d2*d
    sgrada = SQRT(gradao)
    sgradb = SQRT(gradbo)
    sgrad  = SQRT(grad)
    phi = k*sgrad/n_76/c
    e = EXP(-Phi)

    zetas = 2._real_8*densbo/n2
    zetas2 = zetas*zetas
    zetasb = -zetas*densao/densbo

    fs = f3*zetas*(zp13-zm13)/fakt2
    fsb = -fs*densao/densbo
    zs = epsps+fs*(epsf-epsp)+f*(epsfs-epsps)
    zsb = epsps+fsb*(epsf-epsp)+f*(epsfs-epsps)
    cs = (p1+2._real_8*p2*rs)/t1-t2*(p3+2._real_8*p4*rs+3*p5*rs2)&
         /t12
    cs = cs*rss
    cs2 = cs*cs
    ds = 1._real_8/(d*2._real_8)*5._real_8/6._real_8*zetas*(zp23-zm23)
    ds2 = ds*ds

    dsb = -ds*densao/densbo
    phis = -k*sgrad*(7._real_8/6._real_8/(c*n_136)+cs/(c2*n_76))

    phip = k/(n_76*c*sgrad)*sgrada
    phipb = k/(n_76*c*sgrad)*sgradb
    t4 = e*(-phip*grad+2._real_8*sgrada)
    ga = c/d/n_43*t4

    zetass = -4._real_8*densbo/n3
    zetassb = 2._real_8*(densao-densbo)/n3

    fss = f3/fakt2*(o3/zp23*zetas2+zp13*zetass+o3/zm23*zetas2-&
         zm13*zetass)
    fssb = f3/fakt2*(zetassb*(zp13-zm13)+o3*zetas*zetasb*&
         (1._real_8/zp23+1._real_8/zm23))

    zss = epspss+fss*(epsf-epsp)+2._real_8*fs*(epsfs-epsps)+f*(epsfss-epspss)
    t3 = p3+2._real_8*p4*rs+3._real_8*p5*rs2
    zssb = epspss+fssb*(epsf-epsp)+fs*(epsfs-epsps)*(1._real_8-densao/densbo)&
         +f*(epsfss-epspss)
    css = (2*p2/t1-(2._real_8*(p1+2._real_8*p2*rs)*t3+t2*(2._real_8*p4+6._real_8*p5*rs))&
         /t12+1/(t12*t1)*2._real_8*t2*t3*t3)*rss*rss+cs/rss*rsss
    dss = -ds**2/d+0.5_real_8/(fakt3*d)*5._real_8/3._real_8*(2._real_8*fakt3*o3/2._real_8*&
         zetas**2*(1._real_8/zp13+1._real_8/zm13)+fakt3*(zp23-zm23)*zetass/2._real_8)
    dssb = 5._real_8/12._real_8*(zetas/d*(zp23-zm23)*(-dsb/d+zetassb/zetas)+&
         zetas/d*2._real_8*o3*zetasb*(1._real_8/zp13+1._real_8/zm13))
    phiss = k*sgrad*(91._real_8/36._real_8/n_196/c+7._real_8*o3/n_136/c2*cs+&
         2._real_8/c3*cs2/n_76-css/c2/n_76)
    t5 = 1._real_8/d*(-phis*c+cs-c*(ds/d+f3/dens))
    t5s = -ds/d2*(-phis*c+cs-c/d*ds-f3*c/dens)+(-phiss*c-phis*cs+css-cs*&
         (ds/d+f3/dens)-c*(-ds2/d2+dss/d-f3/n2))/d
    nana = 2._real_8*zs+dens*zss+grad*e*(-phis/n_43*t5-f3/n_73*t5+t5s/n_43)
    t5sb = -dsb/d2*(-phis*c+cs-c/d*ds-f3*c/dens)+(-phiss*c-phis*cs+css-&
         cs*(ds/d+f3/dens)-c*(-ds*dsb/d2+dssb/d-f3/n2))/d

    nanb = zsb+zs+dens*zssb+grad*e*(-phis/n_43*t5-f3/n_73*t5+t5sb/n_43)

    phisp = k*(-7._real_8/6._real_8/n_136/c-cs/c2/n_76)/sgrad*sgrada
    naga = 1._real_8/(n_43*d)*(cs*t4+c*e*(phis*phip*grad-phisp*grad-&
         phis*2._real_8*sgrada))+c*t4*(-f3/n_73/d-ds/d2/n_43)

    ganb = 1._real_8/(n_43*d)*(cs*t4+c*e*(phis*phip*grad-phisp*grad-&
         phis*2._real_8*sgrada))+c*t4*(-f3/n_73/d-dsb/d2/n_43)

    nagb = naga*sgradb/sgrada

    phipp = k/(n_76*c)*(1._real_8/sgrad-sgrada**2/(grad*sgrad))
    gaga = c/(n_43*d)*e*((phip**2-phipp)*grad-phip*2._real_8*sgrada+2._real_8-&
         2._real_8*sgrada*phip)

    phippb = k/(n_76*c)*(-sgrada*sgradb/(grad*sgrad))
    gagb = c/(n_43*d)*e*((phipb*phip-phippb)*grad-phip*2._real_8*sgradb&
         -phipb*2._real_8*sgrada)

    phisch = k/(c*n_76*sgrad)
    sk = ga/sgrada

    phissch = k*(-7._real_8/6._real_8/n_136/c-cs/c2/n_76)/sgrad
    nask = nagb/sgradb

    sknb = 1._real_8/(n_43*d)*e*(-phisch*grad+2._real_8)*(cs-f3*c/dens-c/d*dsb&
         -c*phis)-c/(n_43*d)*e*phissch*grad

    phischsch = -k/(n_76*c)/(grad*sgrad)
    sksk = c/(n_43*d)*e*(phisch**2*grad-phischsch*grad-4._real_8*phisch)

    gask = gagb/sgradb
    skgb = sksk*sgradb

    IF (scratch_alpha%n > 1.e-24_real_8 .AND. (scratch_alpha%n < cntr%gceps .OR. &
                                               scratch_alpha%g_2 <= 1.e-48_real_8)) THEN
       p86a = 2.0_real_8*zs+dens*zss
       p86b = zsb+zs+dens*zssb
       functional%tmp1 = p86a*scratch_alpha%dn + p86b*scratch_beta%dn
    ELSEIF (scratch_alpha%n > 1.e-24_real_8) THEN
       skaa = scratch_alpha%gx_2*scratch_alpha%dgx+scratch_alpha%gy_2*scratch_alpha%dgy&
             +scratch_alpha%gz_2*scratch_alpha%dgz
       skab = scratch_alpha%gx_2*scratch_beta%dgx+scratch_alpha%gy_2*scratch_beta%dgy&
             +scratch_alpha%gz_2*scratch_beta%dgz
       skba = scratch_beta%gx_2*scratch_alpha%dgx+scratch_beta%gy_2*scratch_alpha%dgy&
             +scratch_beta%gz_2*scratch_alpha%dgz
       skbb = scratch_beta%gx_2*scratch_beta%dgx+scratch_beta%gy_2*scratch_beta%dgy&
             +scratch_beta%gz_2*scratch_beta%dgz
       skalar = skaa/sgrada
       functional%tmp1 = nana*scratch_alpha%dn + nanb*scratch_beta%dn&
                         + naga*skalar + nask*(skba+skab) + nagb*skbb/sgradb 
       functional%tmp2 = (naga*scratch_alpha%dn + ganb*scratch_beta%dn &
                         + gask*(skba+skab) + gaga*skalar + gagb*skbb/sgradb&
                         - ga*skalar/sgrada)/sgrada
       functional%tmp3 = (ga)/sgrada
       functional%tmp4 = nask*scratch_alpha%dn + sknb*scratch_beta%dn + gask*skalar &
                         + skgb*skbb/sgradb + sksk*(skba+skab)
       functional%tmp5 = sk
    ENDIF

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_spin_dgga_c_p86_compute
  ! ==================================================================
  ELEMENTAL SUBROUTINE cp_spin_dgga_c_pbe_compute( scratch_alpha, scratch_beta, functional )
    ! ==--------------------------------------------------------------==
    ! Ported from "old" xc_ana code, not optimised or prettified
    !                                     19.06.2017 M.P. Bircher @ EPFL
    ! ==--------------------------------------------------------------==

    TYPE(cp_dxc_scratch_t), INTENT(in)       :: scratch_alpha, scratch_beta
    TYPE(cp_dxc_functional_t), INTENT(inout) :: functional


    REAL(real_8), PARAMETER :: c = 1.709921_real_8, &
      fakt = 0.6203504908994000_real_8, &
      iota = 0.071599657785951923081_real_8, &
      k1 = 1.0077158813689795508_real_8, k2 = 0.031090690869654895_real_8 , &
      k3 = -41.156312099041145348_real_8, k4 = 2.1461263399673646085_real_8 , &
      lambda = 0.06672455060314922_real_8, lg = 8.723_real_8, &
      ny = 15.755501913376439197_real_8, &
      o = 6.46309613581743018602522_real_8, &
      o3 = 0.3333333333333333333333_real_8, phig = 0.007389_real_8, &
      t1 = 0.031091_real_8, t2 = 0.015545_real_8, t3 = 0.016887_real_8, &
      u1 = 0.2137_real_8, u2 = 0.20548_real_8, u3 = 0.11125_real_8, &
      v1 = 7.5957_real_8, v2 = 14.1189_real_8, &
      v3 = 10.357_real_8, w1 = 3.5876_real_8, w2 = 6.1977_real_8, &
      w3 = 3.6231_real_8, x1 = 1.6382_real_8, x2 = 3.3662_real_8, &
      x3 = 0.88026_real_8, xi = 23.266_real_8, y = 0.472_real_8, &
      y1 = 0.49294_real_8, y2 = 0.62517_real_8, y3 = 0.49671_real_8, &
      z = -0.001667_real_8

    REAL(real_8) :: nana, naga, gaga, nask, gask, sk, ga, sksk, &
      gagb, nagb, skgb, sknb, ganb, nanb, skaa, skab, skba, skbb
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


    densao = MAX(scratch_alpha%n,1.0e-24_real_8)
    densbo = MAX(scratch_beta%n, 1.0e-24_real_8)
    gradao = MAX(scratch_alpha%g_2,1.0e-48_real_8)
    gradbo = MAX(scratch_beta%g_2, 1.0e-48_real_8)

    dens   = densao+densbo
    n2     = dens*dens
    n3     = n2*dens
    n_13   = dens**o3
    n_43   = n_13*dens
    n_73   = n_43*dens
    n_76   = SQRT(n_73)
    n_136  = n_76*dens
    n_196  = n_136*dens
    n_103  = n_73*dens
    na_13  = densao**o3
    na_23  = na_13*na_13
    nb_13  = densbo**o3
    nb_23  = nb_13*nb_13
    skalar = scratch_alpha%gx_2*scratch_beta%gx_2 + &
             scratch_alpha%gy_2*scratch_beta%gy_2 + &
             scratch_alpha%gz_2*scratch_beta%gz_2 
    grad   = gradao+gradbo+2._real_8*skalar
    grad2  = grad*grad
    sgrad  = SQRT(grad)

    !
    ! Local part
    !
    rs = fakt/n_13
    rs2 = rs*rs
    rs3 = rs2*rs
    rs12 = SQRT(rs)
    rsa = -o3*fakt/n_43
    rsaa = 4._real_8/9._real_8*fakt/n_73
    zeta = (densao-densbo)/dens
    zeta2 = zeta*zeta
    zeta3 = zeta2*zeta
    zeta4 = zeta3*zeta
    zetaa = 2._real_8*densbo/n2
    zp13 = (1._real_8+zeta)**o3
    zm13 = (1._real_8-zeta)**o3
    zp23 = zp13*zp13
    zm23 = zm13*zm13
    zp43 = zp23*zp23
    zm43 = zm23*zm23
    zp53 = zp43*zp13
    zm53 = zm43*zm13
    ome  = (zp43+zm43-2._real_8)/(2._real_8**(4._real_8/3._real_8)-2._real_8)
    omea = o*densbo*(na_13-nb_13)/n_73

    f1 = 2._real_8*t1*(v1*rs12+w1*rs+x1*rs*rs12+y1*rs2)
    f12 = f1*f1
    g1 = 1._real_8+1._real_8/f1
    f1a = 2._real_8*t1*(0.5_real_8*v1/rs12+w1+3._real_8/2._real_8*x1*rs12+2._real_8*y1*rs)*rsa
    e1 = -2._real_8*t1*(1._real_8+u1*rs)*LOG(g1)
    e1a = -2._real_8*t1*(u1*rsa*LOG(g1)-(1._real_8+u1*rs)/g1/f12*f1a)

    f2 = 2._real_8*t2*(v2*rs12+w2*rs+x2*rs*rs12+y2*rs2)
    f22 = f2*f2
    g2 = (f2+1._real_8)/f2
    f2a = 2._real_8*t2*(0.5_real_8*v2/rs12+w2+3._real_8/2._real_8*x2*rs12+2._real_8*y2*rs)*rsa
    e2 = -2._real_8*t2*(1._real_8+u2*rs)*LOG(g2)
    e2a = -2._real_8*t2*(u2*rsa*LOG(g2)-(1._real_8+u2*rs)/g2/f22*f2a)

    f3 = 2._real_8*t3*(v3*rs12+w3*rs+x3*rs*rs12+y3*rs2)
    f32 = f3*f3
    g3 = (f3+1._real_8)/f3
    f3a = 2._real_8*t3*(0.5_real_8*v3/rs12+w3+3._real_8/2._real_8*x3*rs12+2._real_8*y3*rs)*rsa
    e3 = -2._real_8*t3*(1._real_8+u3*rs)*LOG(g3)
    e3a = -2._real_8*t3*(u3*rsa*LOG(g3)-(1._real_8+u3*rs)/g3/f32*f3a)

    zetaaa = -4._real_8*densbo/n3
    zetaab = 2._real_8*zeta/n2
    omeaa = o*densbo*(o3/na_23/n_73-(na_13-nb_13)*7._real_8/3._real_8/n_103)
    omeab = o*((na_13-nb_13)/n_73+densbo*(-o3/nb_23/n_73-(na_13-nb_13)&
         *7._real_8/3._real_8/n_103))
    t = (1._real_8-zeta4)/c
    tt = 4._real_8*zeta3*zetaa

    eps = e1-e3*ome*t+(e2-e1)*ome*zeta4

    kl1 = omea*t-ome*tt/c
    kl2 = omea*zeta4+ome*tt

    epsa = e1a-e3a*ome*t-e3*kl1+(e2a-e1a)*ome*zeta4+(e2-e1)*kl2

    f1aa = 2._real_8*t1*((-1._real_8/4._real_8*v1/rs/rs12+3._real_8/4._real_8*x1/rs12+2._real_8*y1)&
         *rsa**2+f1a/(2._real_8*T1*rsa)*rsaa)
    e1aa = -2._real_8*t1*(u1*rsaa*LOG(g1)-2._real_8*u1*rsa/g1/f12*f1a-(1._real_8+u1*rs&
         )*((f1a/f12/g1)**2+1._real_8/g1*(-2._real_8/f1/f12*f1a**2+f1aa/f12)))

    f2aa = 2._real_8*t2*((-1._real_8/4._real_8*v2/rs/rs12+3._real_8/4._real_8*x2/rs12+2._real_8*y2)&
         *rsa**2+f2a/(2._real_8*T2*rsa)*rsaa)
    e2aa = -2._real_8*t2*(u2*rsaa*LOG(g2)-2._real_8*u2*rsa/g2/f22*f2a-(1._real_8+u2*rs&
         )*((f2a/f22/g2)**2+(-2._real_8/g2/f2/f22*f2a**2+f2aa/g2/f22)))

    f3aa = 2._real_8*t3*((-1._real_8/4._real_8*v3/rs/rs12+3._real_8/4._real_8*x3/rs12+2._real_8*y3)&
         *rsa**2+f3a/(2._real_8*T3*rsa)*rsaa)
    e3aa = -2._real_8*t3*(u3*rsaa*LOG(g3)-2._real_8*u3*rsa/g3/f32*f3a-(1._real_8+u3*rs&
         )*((f3a/f32/g3)**2+(-2._real_8/g3/f3/f32*f3a**2+f3aa/g3/f32)))

    kl1a = omeaa*t-2._real_8*omea*tt/c-ome/c*(12._real_8*zeta2*zetaa**2+4._real_8*zeta3&
         *zetaaa)
    kl2a = omeaa*zeta4+2._real_8*omea*tt+ome*(12._real_8*zeta2*zetaa**2+4._real_8*zeta3&
         *zetaaa)
    epsaa = e1aa-e3aa*ome*t-e3a*(omea*t-ome*tt/c)-e3a*kl1&
         -e3*kl1a+(e2aa-e1aa)*ome*zeta4+(e2a-e1a)*(omea*zeta4+ome*tt)&
         +(e2a-e1a)*kl2+(e2-e1)*kl2a

    nana = 2._real_8*epsa+dens*epsaa

    zetab = -densao/densbo*zetaa
    omeb = -densao/densbo*omea
    epsb = e1a-e3a*ome*(1-zeta4)/c-e3*(omeb*(1-zeta4)/c-ome*4._real_8*zeta3*&
         zetab/c)+(e2a-e1a)*ome*zeta4+(e2-e1)*&
         (omeb*zeta4+ome*4._real_8*zeta3*zetab)

    ttb = 4._real_8*zeta3*zetab
    kl1b = omeab*t-omea*ttb/c-omeb*tt/c-ome/c*(12._real_8*zeta2*zetaa*zetab+&
         4._real_8*zeta3*zetaab)
    kl2b = omeab*zeta4+omea*ttb+omeb*tt+ome*(12._real_8*zeta2*zetaa*zetab+&
         4._real_8*zeta3*zetaab)
    epsab = e1aa-e3aa*ome*t-e3a*(omeb*t-ome*ttb/c)-e3a*kl1&
         -e3*kl1b+(e2aa-e1aa)*ome*zeta4+(e2a-e1a)*(omeb*zeta4+ome*ttb)&
         +(e2a-e1a)*kl2+(e2-e1)*kl2b

    nanb = epsa+epsb+dens*epsab

    IF ( scratch_alpha%n > 1.0e-24_real_8 .AND. (scratch_alpha%n < cntr%gceps .OR. &
                                                 scratch_alpha%g_2 <= 1.e-48_real_8) ) THEN  
       functional%tmp1 = nana*scratch_alpha%dn + nanb*scratch_beta%dn
    ELSEIF ( scratch_alpha%n > 1.0e-24_real_8 ) THEN
       !
       ! Non-local part
       !
       u = 0.5_real_8*(zp23+zm23)
       uh2 = u*u
       uh3 = uh2*u
       uh4 = uh3*u
       uh5 = uh3*uh2
       ua = (o3/zp13-o3/zm13)*zetaa
       ub = (o3/zp13-o3/zm13)*zetab
       d = sgrad/(4._real_8*u)*k1/n_76
       d2 = d*d
       d3 = d2*d
       d4 = d2*d2
       da = sgrad/4._real_8*k1*(-1._real_8/uh2*ua/n_76-7._real_8/6._real_8/u/n_136)
       db = sgrad/4._real_8*k1*(-1._real_8/uh2*ub/n_76-7._real_8/6._real_8/u/n_136)

       ex = EXP(-2._real_8*iota*eps/uh3/lambda**2)
       exa = -1._real_8/k2*(epsa/uh3-eps*3._real_8*ua/uh4)
       exb = -1._real_8/k2*(epsb/uh3-eps*3._real_8*ub/uh4)
       a = k4/(ex-1._real_8)
       a2 = a*a
       Aa = -k4*1._real_8/(ex-1._real_8)**2*ex*exa
       ab = -k4*1._real_8/(ex-1._real_8)**2*ex*exb

       ze1 = d2+a*d4
       ze1a = 2._real_8*d*da+Aa*d4+a*4._real_8*d3*da
       ze1b = 2._real_8*d*db+ab*d4+a*4._real_8*d3*db
       ne1 = 1._real_8+a*d2+a2*d4
       ne1a = Aa*d2+a*2._real_8*d*da+2._real_8*a*aa*d4+a2*4._real_8*d3*da
       ne1b = ab*d2+a*2._real_8*d*db+2._real_8*a*ab*d4+a2*4._real_8*d3*db
       ne12 = ne1*ne1
       b = 1._real_8+k4*ze1/ne1
       b2 = b*b
       Ba = k4*(ze1a/ne1-ze1*ne1a/ne12)
       Bb = k4*(ze1b/ne1-ze1*ne1b/ne12)

       l = uh3*k2*LOG(b)
       La = k2*(3._real_8*uh2*ua*LOG(b)+uh3/b*Ba)
       Lb = k2*(3._real_8*uh2*ub*LOG(b)+uh3/b*Bb)

       uaa = zetaaa*ua/zetaa+0.5_real_8*zetaa**2*&
            (-2._real_8/9._real_8/zp43-2._real_8/9._real_8/zm43)
       uab = zetaab*ua/zetaa+0.5_real_8*zetaa*zetab*&
            (-2._real_8/9._real_8/zp43-2._real_8/9._real_8/zm43)

       daa = sgrad/4._real_8*k1*(2._real_8/uh3*ua**2/n_76-1._real_8/uh2*(uaa/n_76-&
            7._real_8/6._real_8*ua/n_136)-7._real_8/6._real_8*(-1._real_8/uh2*ua/n_136-13._real_8/&
            6._real_8/u/n_196))
       dab = sgrad/4._real_8*k1*(2._real_8/uh3*ua*ub/n_76-1._real_8/uh2*(uab/n_76-&
            7._real_8/6._real_8*ua/n_136)-7._real_8/6._real_8*(-1._real_8/uh2*ub/n_136-13._real_8/&
            6._real_8/u/n_196))

       exaa = -1._real_8/k2*(epsaa/uh3-epsa*6._real_8*ua/uh4-&
            3._real_8*eps*(uaa/uh4-ua**2*4._real_8/uh5))
       exab = -1._real_8/k2*(epsab/uh3-epsa*3._real_8*ub/uh4-&
            epsb*3._real_8*ua/uh4-3._real_8*eps*(uab/uh4-ua*ub*4._real_8/uh5))

       Aaa = -k4*(-2._real_8/(ex-1._real_8)**3*ex**2*exa**2+&
            1._real_8/(ex-1._real_8)**2*(ex*exa**2+ex*exaa))
       Aab = -k4*(-2._real_8/(ex-1._real_8)**3*ex**2*exa*exb+&
            1._real_8/(ex-1._real_8)**2*(ex*exa*exb+ex*exab))

       ze1aa = 2._real_8*(da**2+d*daa)+Aaa*d4+Aa*8._real_8*d3*da+&
            4._real_8*A*(3._real_8*d2*da**2+d3*daa)

       ze1ab = 2._real_8*(db*da+d*dab)+Aab*d4+Aa*4._real_8*d3*db+4._real_8*ab*d3*da+&
            4._real_8*A*(3._real_8*d2*db*da+d3*dab)

       ne1aa = Aaa*d2+Aa*2._real_8*d*da+2._real_8*aa*d*da+2._real_8*a*(da**2+d*daa)+&
            2._real_8*Aa*Aa*d4+2._real_8*A*(Aaa*d4+Aa*4._real_8*d3*da)&
            +8._real_8*A*Aa*d3*da+4._real_8*A2*(3._real_8*d2*da**2+d3*daa)

       ne1ab = Aab*d2+Aa*2._real_8*d*db+2._real_8*ab*d*da+2._real_8*a*(db*da+d*dab)+&
            2._real_8*Ab*Aa*d4+2._real_8*A*(Aab*d4+Aa*4._real_8*d3*db)&
            +8._real_8*A*Ab*d3*da+4._real_8*A2*(3._real_8*d2*db*da+d3*dab)

       Baa = k4*(ze1aa/ne1-ze1a*ne1a/ne12-(ze1a*ne1a+ze1*&
            ne1aa)/ne12+ze1*ne1a**2*2._real_8*ne1/ne12/ne12)

       Bab = k4*(ze1ab/ne1-ze1a*ne1b/ne12-(ze1b*ne1a+ze1*&
            ne1ab)/ne12+ze1*ne1a*ne1b*2._real_8*ne1/ne12/ne12)

       Laa = k2*(6._real_8*uh2*ua/b*Ba+LOG(b)*(6._real_8*u*ua**2+&
            3._real_8*uh2*uaa)+uh3*(-Ba**2/B2+Baa/B))

       Lab = k2*(3._real_8*uh2*ua/b*Bb+LOG(b)*(6._real_8*u*ub*ua+&
            3._real_8*uh2*uab)+3._real_8*uh2*ub/B*Ba+uh3*(-Ba/B2*Bb+Bab/B))

       ze1ga = (d2+a*2._real_8*d4)/grad

       ne1ga = (a*d2+a2*2._real_8*d4)/grad
       Bga = k4*(ze1ga/ne1-ze1*ne1ga/ne12)


       Lga = k2*uh3*Bga/b

       ze1gaga = a*2._real_8*d4/grad2
       ne1gaga = a2*2._real_8*d4/grad2

       Bgaga = k4*(ze1gaga/ne1-ze1ga*ne1ga/ne12-(ze1ga*ne1ga+ze1*&
            ne1gaga)/ne12+ze1*ne1ga**2*2._real_8*ne1/ne12/ne12)

       Lgaga = k2*uh3*(-Bga**2/b2+Bgaga/b)

       ze1naga = (2._real_8*d*da+Aa*2._real_8*d4+8._real_8*a*d3*da)/grad
       ze1ganb = (2._real_8*d*db+ab*2._real_8*d4+8._real_8*a*d3*db)/grad
       ne1naga = (Aa*d2+a*2._real_8*d*da+4._real_8*a*aa*d4+8._real_8*a2*d3*da)/grad
       ne1ganb = (ab*d2+a*2._real_8*d*db+4._real_8*a*ab*d4+8._real_8*a2*d3*db)/grad

       Bnaga = k4*(ze1naga/ne1-ze1ga*ne1a/ne12-(ze1a*ne1ga+ze1*&
            ne1naga)/ne12+ze1*ne1ga*ne1a*2._real_8*ne1/ne12/ne12)
       Bganb = k4*(ze1ganb/ne1-ze1ga*ne1b/ne12-(ze1b*ne1ga+ze1*&
            ne1ganb)/ne12+ze1*ne1ga*ne1b*2._real_8*ne1/ne12/ne12)

       Lnaga = k2*(3._real_8*uh2*ua*Bga/b+uh3*(-bga/b2*Ba+Bnaga/b))
       Lganb = k2*(3._real_8*uh2*ub*Bga/b+uh3*(-bga/b2*Bb+Bganb/b))

       Ja = 0._real_8
       Jaa = Ja
       Jab = Ja
       Jga = Ja
       Jnaga = Ja
       Jganb = Ja
       Jgaga = Ja
       Jb = Ja

       nana = nana+2._real_8*(La+Ja)+dens*(Laa+Jaa)
       nanb = nanb+Lb+Jb+Ja+La+dens*(Lab+Jab)
       ga = dens*(Lga+Jga)
       gaga = dens*(Lgaga+Jgaga)
       naga = Lga+Jga+dens*(Lnaga+Jnaga)
       ganb = Lga+Jga+dens*(Lganb+Jganb)

       nask = 2._real_8*naga
       gask = 2._real_8*gaga
       sk = 2._real_8*ga
       sksk = 4._real_8*gaga
       gagb = gaga
       nagb = naga
       skgb = 2._real_8*gaga
       sknb = 2._real_8*ganb

       skaa = scratch_alpha%gx_2*scratch_alpha%dgx+scratch_alpha%gy_2*scratch_alpha%dgy&
             +scratch_alpha%gz_2*scratch_alpha%dgz
       skab = scratch_alpha%gx_2*scratch_beta%dgx+scratch_alpha%gy_2*scratch_beta%dgy&
             +scratch_alpha%gz_2*scratch_beta%dgz
       skba = scratch_beta%gx_2*scratch_alpha%dgx+scratch_beta%gy_2*scratch_alpha%dgy&
             +scratch_beta%gz_2*scratch_alpha%dgz
       skbb = scratch_beta%gx_2*scratch_beta%dgx+scratch_beta%gy_2*scratch_beta%dgy&
             +scratch_beta%gz_2*scratch_beta%dgz

       functional%tmp1 = (nana)*scratch_alpha%dn + nanb*scratch_beta%dn &
                         +naga*2._real_8*skaa+nask*(skba+skab)+nagb*2._real_8*skbb
       functional%tmp2 = 2._real_8*(naga*scratch_alpha%dn + ganb*scratch_beta%dn +gask*&
                         (skba+skab)+gaga*2._real_8*skaa+gagb*2._real_8*skbb)
       functional%tmp3 = 2._real_8*(ga)
       functional%tmp4 = nask*scratch_alpha%dn + sknb*scratch_beta%dn + gask*2._real_8*skaa+&
                         skgb*2._real_8*skbb+sksk*(skba+skab)
       functional%tmp5 = sk
    ENDIF

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_spin_dgga_c_pbe_compute
  ! ==================================================================
END MODULE cp_dgga_correlation_utils
