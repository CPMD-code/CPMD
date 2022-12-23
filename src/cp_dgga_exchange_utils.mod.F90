! ==================================================================
! Provides: - GGA exchange derivatives for linear response
!
! Input:  Spin-polarised densities from dd_xc_ana_utils
! Output: tmp1...5 analogous to dd_xc_ana_utils
!
!                             02.10.2016 - M. P. Bircher @ LCBC/EPFL
! ==================================================================
MODULE cp_dgga_exchange_utils
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

  PUBLIC :: CP_dGGA_X_B88
  PUBLIC :: CP_dGGA_X_PBE
  PUBLIC :: CP_SPIN_dGGA_X_B88
  PUBLIC :: CP_SPIN_dGGA_X_PBE

CONTAINS

  ! ==================================================================
  PURE SUBROUTINE cp_spin_dgga_x_b88( scratch,functional )
    ! ==--------------------------------------------------------------==

    TYPE(cp_dxc_scratch_t), &
      DIMENSION(cp_xc_spin_pairs), &
      INTENT(in)                             :: scratch
    TYPE(cp_dxc_functional_t), &
      DIMENSION(cp_xc_spin_pairs), &
      INTENT(inout)                          :: functional

    CALL cp_dgga_x_b88( scratch(alpha), functional(alpha) )
    CALL cp_dgga_x_b88( scratch(beta ), functional(beta ) )

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_spin_dgga_x_b88
  ! ==================================================================
  PURE SUBROUTINE cp_dgga_x_b88( scratch,functional )
    ! ==--------------------------------------------------------------==
    ! Ported from "old" xc_ana code, not optimised or prettified
    !                                     19.06.2017 M.P. Bircher @ EPFL
    ! ==--------------------------------------------------------------==

    TYPE(cp_dxc_scratch_t), &
      INTENT(in)                             :: scratch
    TYPE(cp_dxc_functional_t), &
      INTENT(inout)                          :: functional

    REAL(real_8), PARAMETER                  :: a = -0.930525736349100025002010218071667_real_8 

    REAL(real_8)                             :: dfdn, dfdg, &
                                                d2fdn2, d2fdndg, d2fdg2, lda

    REAL(real_8) :: beta6, d1, d12, d13, D1p, D1pp, D1s, D1sp, D1ss, degrado, &
      denso, n1, n2, n_103, n_13, n_23, n_43, n_53, n_73, n_83, s, sgrad, &
      sp, sp3, spp, ss, ssp, sss, t1, t2, t3, t4, t5, t53, x, x2, xp, xs, &
      xsp, xss, skaa, skalar


    denso = MAX(scratch%n,1.e-24_real_8)
    degrado = MAX(scratch%g_2,1.e-49_real_8)

    n2 = denso*denso
    n_13 = denso**(1._real_8/3._real_8)
    n_23 = n_13**2
    n_43 = n_13*denso
    n_53 = n_23*denso
    n_73 = n_43*denso
    n_83 = n_53*denso
    n_103 = n_73*denso
    n1 = 1._real_8/denso
    n_73 = 1._real_8/n_73

    beta6 = 6._real_8*func2%bbeta

    sgrad = SQRT(degrado)
    x = sgrad/n_43
    x2 = x*x
    t1 = SQRT(x2+1._real_8)
    t5 = 1._real_8/t1
    t53 = t5**3
    s = LOG(x+t1)
    t2 = 1._real_8+beta6*x*s
    d1 = 1._real_8/t2
    d12 = d1*d1
    d13 = d12*d1
    xs = -four_thirds*sgrad*n_73
    xp = 1._real_8/n_43
    sp = 1._real_8/SQRT(n_83+degrado)
    sp3 = sp**3
    ss = t5*xs
    D1s = -d12*beta6*(xs*s+x*ss)
    D1p = -d12*beta6*(xp*s+x*sp)
    dfdn = -func2%bbeta*sgrad*(d1*(four_thirds*n1*x+2._real_8*xs)+x*D1s)
    dfdg = -func2%bbeta*(d1*x*2._real_8+n_43*x2*D1p)

    xss = 28._real_8/9._real_8*sgrad/n_103
    t3 = beta6*(xs*s+x*ss)
    sss = -t53*x*xs**2+t5*xss
    D1ss = -2._real_8*d1*D1s*t3-d12*beta6*(xss*s+2._real_8*xs*ss+x*sss)
    d2fdn2 = -func2%bbeta*sgrad*(xs*(four_thirds*n1*d1+3._real_8*D1s)+x*(-four_thirds/n2*d1+four_thirds*n1&
         *D1s+D1ss)+2._real_8*xss*d1)

    xsp = -four_thirds*n_73
    ssp = -sp3*four_thirds*n_53
    D1sp = -2._real_8*d1*D1p*t3-d12*beta6*(xsp*s+xs*sp+xp*ss+x*ssp)
    d2fdndg = dfdn/sgrad-func2%bbeta*sgrad*((xp*d1+x*D1p)*four_thirds*n1&
         +2._real_8*xsp*d1+2._real_8*xs*D1p+xp*D1s+x*D1sp)

    spp = -sp3*sgrad
    t4 = beta6*(xp*s+x*sp)
    D1pp = -2._real_8*d1*D1p*t4-d12*beta6*(2._real_8*xp*sp+x*spp)
    d2fdg2 = -func2%bbeta*(D1p*2._real_8*x+d1*2._real_8*xp+n_43&
         *(2._real_8*x*xp*D1p+x2*D1pp))

    lda = four_thirds/3._real_8*a/n_23

    IF ( scratch%n > 0.0_real_8 .AND. (scratch%n < cntr%gceps .OR. &
                                       scratch%g_2 <=  0.0_real_8) ) THEN                           
       functional%tmp1 = lda*scratch%dn
    ELSE
       skaa = scratch%gx_2*scratch%dgx+scratch%gy_2*scratch%dgy&
             +scratch%gz_2*scratch%dgz
       skalar = skaa/sgrad
       functional%tmp1 = (d2fdn2+lda)*scratch%dn+d2fdndg*skalar
       functional%tmp2 = 1._real_8/sgrad*(d2fdndg*scratch%dn+&
                         d2fdg2*skalar-dfdg*skalar/sgrad)
       functional%tmp3 = dfdg/sgrad
    ENDIF

  END SUBROUTINE cp_dgga_x_b88
  ! ==================================================================
  PURE SUBROUTINE cp_spin_dgga_x_pbe( scratch,functional )
    ! ==--------------------------------------------------------------==

    TYPE(cp_dxc_scratch_t), &
      DIMENSION(cp_xc_spin_pairs), &
      INTENT(in)                             :: scratch
    TYPE(cp_dxc_functional_t), &
      DIMENSION(cp_xc_spin_pairs), &
      INTENT(inout)                          :: functional

    CALL cp_dgga_x_pbe( scratch(alpha), functional(alpha) )
    CALL cp_dgga_x_pbe( scratch(beta ), functional(beta ) )

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_spin_dgga_x_pbe
  ! ==================================================================
  PURE SUBROUTINE cp_dgga_x_pbe( scratch,functional )
    ! ==--------------------------------------------------------------==
    ! Ported from "old" xc_ana code, not optimised or prettified
    !                                     19.06.2017 M.P. Bircher @ EPFL
    ! ==--------------------------------------------------------------==

    TYPE(cp_dxc_scratch_t), &
      INTENT(in)                             :: scratch
    TYPE(cp_dxc_functional_t), &
      INTENT(inout)                          :: functional

    REAL(real_8), PARAMETER :: a = 0.16162045967399548133_real_8, &
      b = 0.98474502184269654115_real_8 , m = 0.2195164512208958179_real_8, &
      r = 0.804_real_8, ax = -0.930525736349100025002010218071667_real_8 

    REAL(real_8) :: dfdg, d2fdn2, d2fdndg, d2fdg2
    REAL(real_8) :: a2, dens0, f, Fg, Fgg, Fng, fs, Fss, grad0, kl, kl2, kl3, &
      n_103, n_13, n_23, n_43, n_73, n_83, o3, s, S2, sgrad, ss, sss, x, xs, &
      xss, skaa, lda


    dens0 = 2.0_real_8*MAX(scratch%n,  1.0e-24_real_8) 
    grad0 = 4.0_real_8*MAX(scratch%g_2,1.0e-48_real_8)
    sgrad = SQRT(grad0)

    a2 = a*a
    o3 = 1._real_8/3._real_8
    n_13 = dens0**o3
    n_23 = n_13*n_13
    n_43 = n_13*dens0
    n_73 = n_43*dens0
    n_83 = n_73*n_13
    n_103 = n_73*dens0

    x = sgrad/n_43
    s = x*a

    s2 = s*s
    kl = 1._real_8+m*s2/r
    kl2 = kl*kl
    kl3 = kl2*kl
    f = 1._real_8+r-r/kl

    xs = -4._real_8/3._real_8*sgrad/n_73
    ss = a*xs
    fs = 2._real_8/kl2*m*s*ss

    xss = 28._real_8/9._real_8*sgrad/n_103
    sss = a*xss
    fss = -m*(8._real_8/r*m/kl3*s2*ss**2-2._real_8/kl2*(ss**2+s*sss))
    d2fdn2 = 2._real_8*b*(-o3/n_23*f-2._real_8*n_13*fs-3._real_8/4._real_8*n_43*fss)

    Fg  = m/kl2*a2/n_83
    Fgg = -m**2*a2**2/n_83*2._real_8/kl3/r/n_83
    Fng = m/grad0*(2._real_8*s*ss/kl2-s2*4._real_8*m/r*s*ss/kl3)

    d2fdndg = 4._real_8*(-b*n_13*Fg-3._real_8/4._real_8*b*n_43*Fng)

    dfdg = -3._real_8/2._real_8*b*n_43*Fg
    d2fdg2 = -6._real_8*b*n_43*Fgg

    IF ( scratch%n > 1.0e-24_real_8 .AND. (scratch%n < cntr%gceps .OR. &
                                           scratch%g_2 <= 1.0e-48_real_8) ) THEN  
       lda             = four_thirds/3._real_8*ax/n_23
       functional%tmp1 = lda*scratch%dn
    ELSEIF ( scratch%n > 1.e-24_real_8 ) THEN
       skaa = scratch%gx_2*scratch%dgx + scratch%gy_2*scratch%dgy &
              + scratch%gz_2*scratch%dgz
       functional%tmp1 = d2fdn2*scratch%dn + d2fdndg*2._real_8*skaa
       functional%tmp2 = 2._real_8*(d2fdndg*scratch%dn + 2._real_8*skaa*d2fdg2)
       functional%tmp3 = 2._real_8*(dfdg)
    ENDIF

  END SUBROUTINE cp_dgga_x_pbe
  ! ==================================================================
END MODULE cp_dgga_exchange_utils
