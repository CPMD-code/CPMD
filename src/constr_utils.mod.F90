MODULE constr_utils
  USE kinds,                           ONLY: real_8
  USE pbc_utils,                       ONLY: dist_pbc, dist_vec

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: funcr
  PUBLIC :: funcd
  PUBLIC :: funcdd
  PUBLIC :: funct
  PUBLIC :: funco
  PUBLIC :: normalize
  PUBLIC :: funcp
  PUBLIC :: diffd
  PUBLIC :: diffdd
  PUBLIC :: diffr
  PUBLIC :: difft
  PUBLIC :: diffo
  PUBLIC :: diffp
  PUBLIC :: getscal
  PUBLIC :: getnorm
  PUBLIC :: normvec
  PUBLIC :: vecmul
  PUBLIC :: addvec
  PUBLIC :: copyvec
  PUBLIC :: zerovec
  PUBLIC :: vecprod
  PUBLIC :: vecrotz
  PUBLIC :: vecroty
  PUBLIC :: vecrotx
  PUBLIC :: vecrotvec

CONTAINS

  ! ----------------STRETCH-----------------------------------------------C
  SUBROUTINE funcr(fr,r,r0,x,y)
    ! FUNCTION: R*R - R0*R0
    ! Argument variables
    REAL(real_8)                             :: fr, r, r0, x(3), y(3)

    REAL(real_8)                             :: t11, diff

    diff = dist_pbc(x,y,.TRUE.)
    fr = diff * diff-r0*r0
    r = diff
    RETURN
  END SUBROUTINE funcr
  ! ----------------DISTANCE----------------------------------------------C
  SUBROUTINE funcd(fd,d,r0,x,y)
    ! FUNCTION: R - R0
    ! Argument variables
    REAL(real_8)                             :: fd, d, r0, x(3), y(3)

    d  = dist_pbc(x,y,.TRUE.)
    fd = d-r0
    RETURN
  END SUBROUTINE funcd
  ! ----------------DIFFERENCE OF DISTANCES-------------------------------C
  SUBROUTINE funcdd(fd,d,r0,x,y,z)
    ! FUNCTION: dR - R0
    ! Argument variables
    REAL(real_8)                             :: fd, d, r0, x(3), y(3), z(3)

    REAL(real_8)                             :: d1, d2

    d1 = dist_pbc(x,y,.TRUE.)
    d2 = dist_pbc(y,z,.TRUE.)
    d  = d1-d2
    fd = d-r0
    RETURN
  END SUBROUTINE funcdd
  ! ----------------ANGLE-------------------------------------------------C
  SUBROUTINE funct(ft,t,t0,x,y,z)
    ! FUNCTION: COS(theta) - COS(theta0)
    ! Argument variables
    REAL(real_8)                             :: ft, t, t0, x(3), y(3), z(3)

    REAL(real_8), PARAMETER                  :: epsilon = 1.e-12_real_8

    REAL(real_8)                             :: ac, caa, cab, cbb, dRxy(3), dRyz(3)

    dRxy = dist_vec(x, y, .TRUE.)
    dRyz = dist_vec(y, z, .TRUE.)

    caa = dot_product(dRxy, dRxy)
    cab = dot_product(dRxy, dRyz)
    cbb = dot_product(dRyz, dRyz)
    IF (caa.LT.epsilon.OR.cbb.LT.epsilon) THEN
       ! 2 points at the same place (T.D. 16/12/1999).
       ! We give the desired value (T0) for satisfied constraint.
       ! ac=-cos(T0)
       ft=0._real_8
       t = t0
    ELSE
       ac = cab/SQRT(caa*cbb)
       IF (ac.LT.-1._real_8) ac=-1._real_8
       IF (ac.GT.1._real_8) ac=1._real_8
       ft = ac+COS(t0)
       t = ACOS(-ac)
    ENDIF
    RETURN
  END SUBROUTINE funct
  ! ----------------DIHEDRAL----------------------------------------------C
  SUBROUTINE funco(fo,o,o0,x,y,z,w,sign0)
    ! FUNCTION: phi - phi0
    ! Argument variables
    REAL(real_8)                             :: fo, o, o0, x(3), y(3), z(3), &
                                                w(3), sign0

    INTEGER                                  :: i
    REAL(real_8)                             :: a(3), ac, avb(3), b(3), &
                                                bvc(3), c(3), pi, s1, vv(3)

    pi=ACOS(-1._real_8)
    a = dist_vec(x, y, .TRUE.)
    b = dist_vec(z, y, .TRUE.)

    CALL normalize(a)
    CALL normalize(b)
    CALL vecprod(a,b,avb)
    CALL normalize(avb)

    b = dist_vec(y, z, .TRUE.)
    c = dist_vec(w, z, .TRUE.)

    CALL normalize(b)
    CALL normalize(c)
    CALL vecprod(b,c,bvc)
    CALL normalize(bvc)

    CALL getscal(avb,bvc,ac)

    CALL vecprod(avb,bvc,vv)

    b = dist_vec(z, y, .TRUE.)

    CALL getscal(b,vv,s1)
    ! avoid NaN from acos
    IF (ac.GT.1._real_8) THEN
       ac=1._real_8
    ELSEIF (ac.LT.-1._real_8)THEN
       ac=-1._real_8
    ENDIF
    o=ACOS(ac)

    ! a = x-y
    ! b = z-y (or y-z)
    ! c = w-z
    ! avb =  a x -b / || a x -b ||
    ! bvc =  b x  c / || b x  c ||
    ! ac = avb . bvc = ( a x -b ) . ( b x c ) / ( || a x -b || * || b x c || )
    ! = - ( ( a x b ) x b ) . c / ( || a x -b || * || b x c || )
    ! = [ ( b . b ) ( a . c ) - ( a . b ) ( b . c ) ] / ( || || * ||  || )

    sign0=1._real_8
    IF (s1.LT.0._real_8)THEN
       sign0=-1._real_8
       o=-o
    ENDIF

    fo=MOD(2.0e4_real_8*pi+o-o0,2.0_real_8*pi)

    IF ( fo .GT. +pi ) fo = fo - 2._real_8*pi

    o=o0+fo

    RETURN
  END SUBROUTINE funco
  ! 
  SUBROUTINE normalize(b)
    REAL(real_8)                             :: b(3)

    REAL(real_8)                             :: norm

    norm=SQRT(b(1)**2+b(2)**2+b(3)**2)
    b(1)=b(1)/norm
    b(2)=b(2)/norm
    b(3)=b(3)/norm
    RETURN
  END SUBROUTINE normalize

  ! ----------------OUT OF PLANE------------------------------------------C
  SUBROUTINE funcp(diff,teta,teta0,x,y,z,w,d)

    ! FUNCTION: phi - phi0
    ! Argument variables
    REAL(real_8)                             :: diff, teta, teta0, x(3), &
                                                y(3), z(3), w(3), d(12)

    INTEGER                                  :: i
    REAL(real_8)                             :: const, cosine, den, nn(3), &
                                                nn2, nn_norm, num, o(3), o2, &
                                                o_norm, pi, q(3), scal, &
                                                sign0, sinus, t(3), t_norm

    sign0=1._real_8
    pi=ACOS(-1._real_8)

    t = dist_vec(x, y, .TRUE.)
    q = dist_vec(z, y, .TRUE.)

    DO i=1,3
       t(i) = (x(i)-y(i))
       q(i) = (z(i)-y(i))
    ENDDO
    t_norm = t(1)*t(1)+t(2)*t(2)+t(3)*t(3)
    t_norm = SQRT(t_norm)
    CALL vecprod(t,q,nn)
    nn2 = nn(1)*nn(1)+nn(2)*nn(2)+nn(3)*nn(3)
    nn_norm = SQRT(nn2)
    DO i=1,3
       o(i) = (w(i)-y(i))
    ENDDO
    o2 = o(1)*o(1)+o(2)*o(2)+o(3)*o(3)
    o_norm = SQRT(o2)

    CALL getscal(o,t,scal)
    CALL getscal(o,nn,cosine)
    num = cosine
    den = 1.0_real_8/(o_norm*nn_norm)
    cosine = cosine*den
    sinus = SQRT(1-cosine*cosine)
    teta = ACOS(cosine)
    IF (scal .LT. 0._real_8) THEN
       sign0=-1._real_8
       teta=-teta+2._real_8*pi
    ENDIF

    diff=teta-MOD(2._real_8*pi+teta0,2._real_8*pi)

    const = -1.0_real_8/sinus*sign0

    d(1) = const*den*(-q(3)*o(2)+q(2)*o(3)&
         -num*(-q(3)*nn(2)+q(2)*nn(3))/nn2)
    d(2) = const*den*(-q(1)*o(3)+q(3)*o(1)&
         -num*(-q(1)*nn(3)+q(3)*nn(1))/nn2)
    d(3) = const*den*(-q(2)*o(1)+q(1)*o(2)&
         -num*(-q(2)*nn(1)+q(1)*nn(2))/nn2)

    d(4) = const*den*(-nn(1)-t(3)*o(2)+q(3)*o(2)-q(2)*o(3)+t(2)*o(3)&
         -num*(-o(1)/o2+&
         ((-t(3)+q(3))*nn(2)+(-q(2)+t(2))*nn(3))/nn2))
    d(5) = const*den*(-nn(2)-q(3)*o(1)+t(3)*o(1)-t(1)*o(3)+q(1)*o(3)&
         -num*(-o(2)/o2+&
         ((-q(3)+t(3))*nn(1)+(-t(1)+q(1))*nn(3))/nn2))
    d(6) = const*den*(-nn(3)-t(2)*o(1)+q(2)*o(1)-q(1)*o(2)+t(1)*o(2)&
         -num*(-o(3)/o2+&
         ((-t(2)+q(2))*nn(1)+(-q(1)+t(1))*nn(2))/nn2))

    d(7) = const*den*(t(3)*o(2)-t(2)*o(3)&
         -num*(t(3)*nn(2)-t(2)*nn(3))/nn2)
    d(8) = const*den*(t(1)*o(3)-t(3)*o(1)&
         -num*(t(1)*nn(3)-t(3)*nn(1))/nn2)
    d(9) = const*den*(t(2)*o(1)-t(1)*o(2)&
         -num*(t(2)*nn(1)-t(1)*nn(2))/nn2)

    d(10) = const*den*(nn(1)-num*o(1)/o2)
    d(11) = const*den*(nn(2)-num*o(2)/o2)
    d(12) = const*den*(nn(3)-num*o(3)/o2)

    RETURN
  END SUBROUTINE funcp

  ! ----------------------------------------------------------------------C
  SUBROUTINE diffd(dR,x,y)
    ! Argument variables
    REAL(real_8)                             :: dR(6), x(3), y(3)

    REAL(real_8), PARAMETER                  :: epsilon = 1.e-12_real_8 

    REAL(real_8)                             :: r, diff(3)

    diff = dist_vec(x, y, .TRUE.)
    r  = SQRT(dot_product(diff, diff))
    IF (r.LT.epsilon) THEN
       dR(:) = 0._real_8
    ELSE
      dR(1:3) = diff / r
      dR(4:6) = -diff / r
    ENDIF
    RETURN
  END SUBROUTINE diffd
  ! ----------------------------------------------------------------------C
  SUBROUTINE diffdd(dR,x,y,z)
    ! Argument variables
    REAL(real_8)                             :: dR(9), x(3), y(3), z(3)

    REAL(real_8), PARAMETER                  :: epsilon = 1.e-12_real_8 

    REAL(real_8)                             :: r1, r2, dRxy(3), dRyz(3)

    dRxy = dist_vec(x, y, .TRUE.)
    dRyz = dist_vec(y, z, .TRUE.)
    r1 = SQRT(dot_product(dRxy, dRxy))
    r2 = SQRT(dot_product(dRyz, dRyz))
    IF (r1.LT.epsilon) THEN
       dR(1:3) = 0._real_8
    ELSE
       dR(1:3) = dRxy/r1
    ENDIF
    IF (r2.LT.epsilon) THEN
       dR(7:9) = 0._real_8
    ELSE
       dR(7:9) = dRyz/r2
    ENDIF
    dR(4:6) = -dR(1:3) - dR(7:9)
    RETURN
  END SUBROUTINE diffdd
  ! ----------------------------------------------------------------------C
  SUBROUTINE diffr(dR,x,y)
    ! Argument variables
    REAL(real_8)                             :: dR(6), x(3), y(3)

    REAL(real_8)                             :: diff(3)

    diff = dist_vec(x, y, .TRUE.)
    dR(1:3) = 2._real_8 * diff
    dR(4:6) = -2._real_8 * diff
    RETURN
  END SUBROUTINE diffr
  ! ----------------------------------------------------------------------C
  SUBROUTINE difft(dT,x,y,z)
    ! Argument variables
    REAL(real_8)                             :: dT(9), x(3), y(3), z(3)

    REAL(real_8), PARAMETER                  :: epsilon = 1.e-12_real_8 

    INTEGER                                  :: i
    REAL(real_8)                             :: caa, cab, cbb, ccc, dRxy(3), &
                                                dRyz(3)

    dRxy = dist_vec(x, y, .TRUE.)
    dRyz = dist_vec(y, z, .TRUE.)

    caa = dot_product(dRxy, dRxy)
    cab = dot_product(dRxy, dRyz)
    cbb = dot_product(dRyz, dRyz)

    IF (caa.LT.epsilon.OR.cbb.LT.epsilon) THEN
       DO i=1,9
          dT(i)=0._real_8
       ENDDO
    ELSE
       ccc = -1._real_8/SQRT(caa*cbb)
       dT(1) =  ccc * (cab / caa * dRxy(1) - dRyz(1))
       dT(2) =  ccc * (cab / caa * dRxy(2) - dRyz(2))
       dT(3) =  ccc * (cab / caa * dRxy(3) - dRyz(3))
       dT(4) = -ccc * (cab / caa * dRxy(1) - cab / cbb * dRyz(1) + dRxy(1) - dRyz(1))
       dT(5) = -ccc * (cab / caa * dRxy(2) - cab / cbb * dRyz(2) + dRxy(2) - dRyz(2))
       dT(6) = -ccc * (cab / caa * dRxy(3) - cab / cbb * dRyz(3) + dRxy(3) - dRyz(3))
       dT(7) = -ccc * (cab / cbb * dRyz(1) - dRxy(1))
       dT(8) = -ccc * (cab / cbb * dRyz(2) - dRxy(2))
       dT(9) = -ccc * (cab / cbb * dRyz(3) - dRxy(3))
    ENDIF
    RETURN
  END SUBROUTINE difft
  ! ----------------------------------------------------------------------C
  SUBROUTINE diffo(mdo,x,y,z,w,sign0)
    ! Argument variables
    REAL(real_8)                             :: mdo(12), x(3), y(3), z(3), &
                                                w(3), sign0

    REAL(real_8), PARAMETER                  :: epsilon = 1.e-12_real_8 

    INTEGER                                  :: i
    REAL(real_8)                             :: ac, caa, cab, cac, cbb, cbc, &
                                                ccc, cost, dab, dbc, ddd, &
                                                dR1(3), dR2(3), dR3(3)

    dR1 = dist_vec(x, y, .TRUE.)
    dR2 = dist_vec(y, z, .TRUE.)
    dR3 = dist_vec(z, w, .TRUE.)

    caa = dot_product(dR1, dR1)
    cbb = dot_product(dR2, dR2)
    ccc = dot_product(dR3, dR3)
    cab = dot_product(dR1, dR2)
    cac = dot_product(dR1, dR3)
    cbc = dot_product(dR2, dR3)
    dab = caa*cbb-cab*cab
    dbc = cbb*ccc-cbc*cbc

    ! || a x b ||^2 = ( a x b ) . ( a x b )
    ! = a . ( b x ( a x b ) )
    ! = a . ( ( b . b ) a - ( b . a ) b
    ! = ( a . a ) ( b . b ) - ( a . b ) ( a . b )
    ! dab = || a x b ||^2

    ! ac = - ( ( a . b ) ( b . c ) - ( a . c ) ( b . b ) )
    ! / ( || a x -b || * || b x c || )
    ! 
    ! D[arccos(ac)] = -1 / sqrt ( 1 - ac^2 )
    ! d[arccos(ac)]/dx1
    ! = -1 / sqrt ( 1 - ac^2 ) * d(ac)/dx1
    ! = 1 / sqrt ( 1 - ac^2 ) * [ b1 ( b . c ) - ( b . b ) c1 ] / ( || || || || )
    ! - 1 / sqrt ( 1 - ac^2 ) * [ ( a . b ) ( b . c ) - ( a . c ) ( b . b ) ]
    ! * d(|| a x b ||) / dx1 * || -b x c || / ( || || || || )^2
    ! = 1 / sqrt ( 1 - ac^2 ) * [ b1 ( b . c ) - ( b . b ) c1 ] / ( || || || || )
    ! - 1 / sqrt ( 1 - ac^2 ) * [ ( a . b ) ( b . c ) - ( a . c ) ( b . b ) ]
    ! * d(|| a x b ||^2) / dx1 / ( 2 * || a x b ||^3 * || -b x c || )
    ! = 1 / sqrt ( 1 - ac^2 ) * [ b1 ( b . c ) - ( b . b ) c1 ] / ( || || || || )
    ! - 1 / sqrt ( 1 - ac^2 ) * [ ( a . b ) ( b . c ) - ( a . c ) ( b . b ) ]
    ! * ( 2 * a1 * ( b . b ) - 2 * b1 * ( a . b ) )
    ! / ( 2 * || a x b ||^3 * || -b x c || )
    ! = 1 / sqrt ( 1 - ac^2 ) * [ b1 ( b . c ) - ( b . b ) c1 ] / ( || || || || )
    ! - 1 / sqrt ( 1 - ac^2 ) * [ ( a . b ) ( b . c ) - ( a . c ) ( b . b ) ]
    ! * ( a1 * ( b . b ) - b1 * ( a . b ) ) / ( || a x b ||^3 * || -b x c || )
    ! = ddd * [ b1 ( b . c ) - ( b . b ) c1 ]
    ! - ddd * [ ( a . b ) ( b . c ) - ( a . c ) ( b . b ) ]
    ! * ( a1 * ( b . b ) - b1 * ( a . b ) ) / || a x b ||^2
    ! = ddd * [ t4 * cbc - cbb * t7 ]
    ! - ddd * [ cab * cbc - cac * cbb ] * ( t1 * cbb - t4 * cab ) / dab

    IF (dab.LT.epsilon.OR.dbc.LT.epsilon) THEN
       DO i=1,12
          mdo(i)=0._real_8
       ENDDO
    ELSE

       ! LAIO
       ac = -(cab*cbc-cac*cbb)/SQRT(dab*dbc)
       ! if(abs(ac).gt.0.99999999_real_8) ac=0.99999999_real_8
       IF (ABS(ac).GT.0.999999999999_real_8) ac=0.999999999999_real_8
       cost= 1._real_8/SQRT(1._real_8-ac*ac)*sign0
       ddd = 1._real_8/SQRT(dab*dbc)*cost

       mdo(1) = -ddd * &
            (cbc * dR2(1) - cbb * dR3(1) - &
            (cab * cbc - cac * cbb) / dab * (cbb * dR1(1) - cab * dR2(1)))

       mdo(2) = -ddd * &
            (cbc * dR2(2) - cbb * dR3(2) - (cab * cbc - cac * cbb) / dab * &
            (cbb * dR1(2) - cab * dR2(2)))

       mdo(3) = -ddd*&
            (cbc * dR2(3) - cbb * dR3(3) - (cab * cbc - cac * cbb) / dab * &
            (cbb * dR1(3) - cab * dR2(3)))

       mdo(4) = -ddd * (cbc * dR1(1) - cbc * dR2(1) + cab * dR3(1) + &
            cbb * dR3(1) - 2 * cac * dR2(1) - &
            (cab * cbc - cac * cbb) / dbc * (ccc * dR2(1) - cbc * dR3(1)) - &
            (cab * cbc - cac * cbb) / dab * (caa * dR2(1) - cbb * dR1(1) - &
            cab * dR1(1) + cab * dR2(1)))

       mdo(5) = -ddd * (cbc * dR1(2) - cbc * dR2(2) + cab * dR3(2) + &
            cbb * dR3(2) - 2 * cac * dR2(2) - &
            (cab * cbc - cac * cbb) / dbc * (ccc * dR2(2) - cbc * dR3(2)) - &
            (cab * cbc - cac * cbb) / dab * (caa * dR2(2) - cbb * dR1(2) - &
            cab * dR1(2) + cab * dR2(2)))

       mdo(6) = -ddd * (cbc * dR1(3) - cbc * dR2(3) + cab * dR3(3) + &
            cbb * dR3(3) - 2 * cac * dR2(3) - &
            (cab * cbc - cac * cbb) / dbc * (ccc * dR2(3) - cbc * dR3(3)) - &
            (cab * cbc - cac * cbb) / dab * (caa * dR2(3) - cbb * dR1(3) - &
            cab * dR1(3) + cab * dR2(3)))

       mdo(7) = -ddd * (-cbc * dR1(1) + cab * dR2(1) - cab * dR3(1) - &
            cbb * dR1(1) + 2._real_8 * cac * dR2(1) - &
            (cab * cbc - cac * cbb) / dbc * (cbb * dR3(1) - ccc * dR2(1) - &
            cbc * dR2(1) + cbc * dR3(1)) - &
            (cab * cbc - cac * cbb) / dab * (-caa * dR2(1) + cab * dR1(1)))

       mdo(8) = -ddd * (-cbc * dR1(2) + cab * dR2(2) - cab * dR3(2) - &
            cbb * dR1(2) + 2._real_8 * cac * dR2(2) - &
            (cab * cbc - cac * cbb) / dbc * (cbb * dR3(2) - &
            ccc * dR2(2) - cbc * dR2(2) + cbc * dR3(2)) - &
            (cab * cbc - cac * cbb) / dab * (-caa * dR2(2) + cab * dR1(2)))

       mdo(9) = -ddd * (-cbc * dR1(3) + cab * dR2(3) - cab * dR3(3) - &
            cbb * dR1(3) + 2._real_8 * cac * dR2(3) - &
            (cab * cbc - cac * cbb) / dbc * (cbb * dR3(3) - &
            ccc * dR2(3) - cbc * dR2(3) + cbc * dR3(3)) - &
            (cab * cbc - cac * cbb) / dab * (-caa * dR2(3) + cab * dR1(3)))

       mdo(10)= -ddd * (-cab * dR2(1) + cbb * dR1(1) - &
                (cab * cbc - cac * cbb) / dbc * &
                (-cbb * dR3(1) + cbc * dR2(1)))

       mdo(11)= -ddd * (-cab * dR2(2) + cbb * dR1(2) - &
                (cab * cbc - cac * cbb) / dbc * &
                (-cbb * dR3(2) + cbc * dR2(2)))

       mdo(12)= -ddd * (-cab * dR2(3) + cbb * dR1(3) - &
                (cab * cbc - cac * cbb) / dbc * &
                (-cbb*dR3(3)+cbc*dR2(3)))
    ENDIF

    RETURN
  END SUBROUTINE diffo
  ! ----------------------------------------------------------------------C
  SUBROUTINE diffp(dP,x,y,z,w)
    ! Argument variables
    REAL(real_8)                             :: dP(12), x(3), y(3), z(3), w(3)

    REAL(real_8), PARAMETER                  :: epsilon = 1.e-12_real_8 

    INTEGER                                  :: i
    REAL(real_8)                             :: cab, ccc, dab, ddd, &
                                                dR1(3), dR2(3), dR3(3)

    dR1 = dist_vec(x, y, .TRUE.)
    dR2 = dist_vec(y, z, .TRUE.)
    dR3 = dist_vec(w, y, .TRUE.)

    ccc = dot_product(dR3, dR3)
    cab = dR3(1) * (dR1(2) * dR2(3) - dR1(3) * dR2(2)) + &
          dR3(2) * (dR1(3) * dR2(1) - dR1(1) * dR2(3)) + &
          dR3(3) * (dR1(1) * dR2(2) - dR1(2) * dR2(1))
    dab = (dR1(2) * dR2(3) - dR1(3) * dR2(2)) * (dR1(2) * dR2(3) - dR1(3) * dR2(2)) + &
          (dR1(3) * dR2(1) - dR1(1) * dR2(3)) * (dR1(3) * dR2(1) - dR1(1) * dR2(3)) + &
          (dR1(1) * dR2(2) - dR1(2) * dR2(2)) * (dR1(1) * dR2(2) - dR1(2) * dR2(2))
    IF (ccc.LT.epsilon.OR.dab.LT.epsilon) THEN
       DO i=1,12
          dP(i)=0._real_8
       ENDDO
    ELSE
       ddd = 1._real_8/SQRT(ccc*dab)
       dP(1) = ddd * (-dR3(2) * dR2(3) + dR3(3) * dR2(2) - &
               cab * (dR1(1) * (dR2(2) * dR2(2) + dR2(3) * dR2(3)) - &
               dR2(1) * (dR1(2) * dR2(2) + dR1(3) *dR2(3))) / dab)
       dP(2) = ddd * (dR3(1) * dR2(3) - dR3(3) * dR2(1) - &
               cab * (dR1(2) * (dR2(3) * dR2(3) + dR2(1) * dR2(1)) - &
               dR2(2) * (dR1(3) * dR2(3) + dR1(1) * dR2(1))) / dab)
       dP(3) = ddd * (-dR3(1) * dR2(2) + dR3(2) * dR2(1) - &
               cab * (dR1(3) * (dR2(1) * dR2(1) + dR2(2) * dR2(2)) - &
               dR2(3) * (dR1(1) * dR2(1) + dR1(2) * dR2(2))) / dab)
       dP(7) = -ddd * (dR3(2) * dR1(3) - dR3(3) * dR1(2) - &
               cab * (dR2(1) * (dR1(2) * dR1(2) + dR1(3) * dR1(3)) - &
               dR1(1) * (dR1(2) * dR2(2) + dR1(3) * dR2(3))) / dab)
       dP(8) = -ddd * (-dR3(1) * dR1(3) + dR3(3) * dR1(1) - &
               cab * (dR2(2) * (dR1(3) * dR1(3) + dR1(1) * dR1(1)) - &
               dR1(2) * (dR1(3) * dR2(3) + dR1(1) * dR2(1))) / dab)
       dP(9) = -ddd * (dR3(1) * dR1(2) - dR3(2) * dR1(1) - &
               cab * (dR2(3) * (dR1(1) * dR1(1) + dR1(2) * dR1(2)) - &
               dR1(3) * (dR1(1) * dR2(1) + dR1(2) * dR2(2))) / dab)
       dP(10)= ddd * (dR1(2) * dR2(3) - dR1(3) * dR2(3) - cab * dR3(1) / ccc)
       dP(11)= ddd * (dR1(3) * dR2(1) - dR1(1) * dR2(3) - cab * dR3(2) / ccc)
       dP(12)= ddd * (dR1(1) * dR2(2) - dR1(2) * dR2(1) - cab * dR3(3) / ccc)
       dP(4) = -dp(1) - dp(7) - dp(10)
       dP(5) = -dp(2) - dp(8) - dp(11)
       dP(6) = -dp(3) - dp(9) - dp(12)
    ENDIF
    RETURN
  END SUBROUTINE diffp
  ! ----------------------------------------------------------------------C

  ! ==================================================================
  SUBROUTINE getscal(vec1,vec2,scal)
    REAL(real_8), DIMENSION(3)               :: vec1, vec2
    REAL(real_8)                             :: scal

! ==--------------------------------------------------------------==

    scal = vec1(1)*vec2(1) + vec1(2)*vec2(2) + vec1(3)*vec2(3)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE getscal
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE getnorm(vec,dnorm)
    REAL(real_8), DIMENSION(3)               :: vec
    REAL(real_8)                             :: dnorm

! ==--------------------------------------------------------------==

    dnorm = SQRT ( vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3) )
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE getnorm
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE normvec(vec_old,vec_norm)
    REAL(real_8), DIMENSION(3)               :: vec_old, vec_norm

    REAL(real_8)                             :: dnorm

! ==--------------------------------------------------------------==

    CALL getnorm(vec_old, dnorm)
    vec_norm(1) = vec_old(1) / dnorm
    vec_norm(2) = vec_old(2) / dnorm
    vec_norm(3) = vec_old(3) / dnorm
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE normvec
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE vecmul(a,vec1,vec2)
    REAL(real_8)                             :: a
    REAL(real_8), DIMENSION(3)               :: vec1, vec2

! ==--------------------------------------------------------------==

    vec2(1) = a * vec1(1)
    vec2(2) = a * vec1(2)
    vec2(3) = a * vec1(3)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE vecmul
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE addvec(a,vec1,b,vec2,vec3)
    REAL(real_8)                             :: a
    REAL(real_8), DIMENSION(3)               :: vec1
    REAL(real_8)                             :: b
    REAL(real_8), DIMENSION(3)               :: vec2, vec3

! ==--------------------------------------------------------------==

    vec3(1)=a*vec1(1)+b*vec2(1)
    vec3(2)=a*vec1(2)+b*vec2(2)
    vec3(3)=a*vec1(3)+b*vec2(3)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE addvec
  ! ==================================================================
  ! 

  ! ==================================================================
  SUBROUTINE copyvec(vec1,vec2)
    REAL(real_8), DIMENSION(3)               :: vec1, vec2

! ==--------------------------------------------------------------==

    vec2(1) = vec1(1)
    vec2(2) = vec1(2)
    vec2(3) = vec1(3)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE copyvec
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE zerovec(vec1)
    REAL(real_8), DIMENSION(3)               :: vec1

! ==--------------------------------------------------------------==

    vec1(1) = 0.0_real_8
    vec1(2) = 0.0_real_8
    vec1(3) = 0.0_real_8
    ! ==--------------------------------------------------------------==

    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE zerovec
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE vecprod(vec1,vec2,vec3)
    REAL(real_8), DIMENSION(3)               :: vec1, vec2, vec3

! ==--------------------------------------------------------------==

    vec3(1) = vec1(2) * vec2(3) - vec1(3) * vec2(2)
    vec3(2) = vec1(3) * vec2(1) - vec1(1) * vec2(3)
    vec3(3) = vec1(1) * vec2(2) - vec1(2) * vec2(1)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE vecprod
  ! ==================================================================
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE vecrotz(angle,vec1,vec2)
    REAL(real_8)                             :: angle
    REAL(real_8), DIMENSION(3)               :: vec1, vec2

    INTEGER                                  :: i, j
    REAL(real_8)                             :: r(3,3)

! ==--------------------------------------------------------------==

    DO i=1,3
       DO j=1,3
          r(i,j)=0.0_real_8
       ENDDO
       vec2(i)=0.0_real_8
       r(i,i)=1.0_real_8
    ENDDO
    r(1,1)=+COS(angle)
    r(2,2)=+COS(angle)
    r(1,2)=-SIN(angle)
    r(2,1)=+SIN(angle)
    DO i=1,3
       DO j=1,3
          vec2(i)=vec2(i)+r(i,j)*vec1(j)
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE vecrotz
  ! ==================================================================
  ! 
  ! ==================================================================
  SUBROUTINE vecroty(angle,vec1,vec2)
    REAL(real_8)                             :: angle
    REAL(real_8), DIMENSION(3)               :: vec1, vec2

    INTEGER                                  :: i, j
    REAL(real_8)                             :: r(3,3)

! ==--------------------------------------------------------------==

    DO i=1,3
       DO j=1,3
          r(i,j)=0.0_real_8
       ENDDO
       vec2(i)=0.0_real_8
       r(i,i)=1.0_real_8
    ENDDO
    r(1,1)=+COS(angle)
    r(3,3)=+COS(angle)
    r(1,3)=+SIN(angle)
    r(3,1)=-SIN(angle)
    DO i=1,3
       DO j=1,3
          vec2(i)=vec2(i)+r(i,j)*vec1(j)
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE vecroty
  ! ==================================================================
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE vecrotx(angle,vec1,vec2)
    REAL(real_8)                             :: angle
    REAL(real_8), DIMENSION(3)               :: vec1, vec2

    INTEGER                                  :: i, j
    REAL(real_8)                             :: r(3,3)

! ==--------------------------------------------------------------==

    DO i=1,3
       DO j=1,3
          r(i,j)=0.0_real_8
       ENDDO
       r(i,i)=1.0_real_8
       vec2(i)=0.0_real_8
    ENDDO
    r(2,2)=+COS(angle)
    r(3,3)=+COS(angle)
    r(2,3)=-SIN(angle)
    r(3,2)=+SIN(angle)
    DO i=1,3
       DO j=1,3
          vec2(i)=vec2(i)+r(i,j)*vec1(j)
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE vecrotx
  ! ==================================================================
  ! 
  ! ==================================================================
  SUBROUTINE vecrotvec(vec1,angle,vec2,vec3)
    REAL(real_8), DIMENSION(3)               :: vec1
    REAL(real_8)                             :: angle
    REAL(real_8), DIMENSION(3)               :: vec2, vec3

    REAL(real_8)                             :: a1, a2, dummy(3,6), pi

! ==--------------------------------------------------------------==

    pi=ACOS(-1.0_real_8)
    a1=(pi/2.0_real_8)-ATAN2(vec1(2),vec1(1)) ! first y than x
    CALL vecrotz(a1,vec1,dummy(1,1))          ! rotate in the yz plane
    CALL vecrotz(a1,vec2,dummy(1,2))
    a2=(pi/2.0_real_8)-ATAN2(dummy(3,1),dummy(2,1))
    CALL vecrotx(a2,dummy(1,1),dummy(1,3))    ! put on the z axis
    CALL vecrotx(a2,dummy(1,2),dummy(1,4))
    ! rotate now around z for the angle
    CALL vecrotz(angle,dummy(1,4),dummy(1,5))
    ! and go back
    CALL vecrotx(-a2,dummy(1,5),dummy(1,6))
    CALL vecrotz(-a1,dummy(1,6),vec3)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE vecrotvec
  ! ==================================================================

END MODULE constr_utils


