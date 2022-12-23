MODULE rk4ov_utils
  USE kinds,                           ONLY: real_8
  USE shop,                            ONLY: sh02
  USE shop_rest,                       ONLY: prob1

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rk4ov_new
  PUBLIC :: rk4ov_old

CONTAINS

  SUBROUTINE rk4ov_new(h,x,d,dd,e,dei,s,ds)
    REAL(real_8)                             :: h, x(3*sh02%nsurf), &
                                                d(sh02%nsurf,sh02%nsurf), &
                                                dd(2,2), e(sh02%nsurf), &
                                                dei(2), s, ds

    COMPLEX(real_8)                          :: a(3), adot(3), b(3), bdot(3), &
                                                d1, d1dot, d2, d2dot, k1(3), &
                                                k2(3), k3(3), k4(3), p(2), &
                                                pdot(2), ui
    INTEGER                                  :: i
    REAL(real_8)                             :: ds1, ds2, ds3, ds4, oms, sdot

    ui=(0.0_real_8,1.0_real_8)
    a(1)=CMPLX(x(1),x(3))
    a(2)=CMPLX(x(2),x(4))
    a(3)=CMPLX(x(5),x(6))
    CALL derivov_new(a,adot,d,e,s,sdot)
    DO i=1,3
       k1(i)=h*adot(i)
       a(i)=a(i)+k1(i)/2.0_real_8
    ENDDO
    ! ...update coupling, overlap, energies
    ds1=h*sdot
    s=s+ds1/2.0_real_8
    ! 
    d(1,2)=d(1,2)+dd(1,2)/2.0_real_8
    d(2,1)=d(2,1)+dd(2,1)/2.0_real_8
    e(1)=e(1)+dei(1)/2.0_real_8
    e(2)=e(2)+dei(2)/2.0_real_8
    ! ...end update coupling
    CALL derivov_new(a,adot,d,e,s,sdot)
    DO i=1,3
       a(i)=a(i)-k1(i)/2.0_real_8
       k2(i)=h*adot(i)
       a(i)=a(i)+k2(i)/2.0_real_8
    ENDDO
    s=s-ds1/2.0_real_8
    ds2=h*sdot
    s=s+ds2/2.0_real_8
    CALL derivov_new(a,adot,d,e,s,sdot)
    DO i=1,3
       a(i)=a(i)-k2(i)/2.0_real_8
       k3(i)=h*adot(i)
       a(i)=a(i)+k3(i)
    ENDDO
    s=s-ds2/2.0_real_8
    ds3=h*sdot
    s=s+ds3
    ! ...update coupling, overlap, energies
    d(1,2)=d(1,2)+dd(1,2)/2.0_real_8
    d(2,1)=d(2,1)+dd(2,1)/2.0_real_8
    e(1)=e(1)+dei(1)/2.0_real_8
    e(2)=e(2)+dei(2)/2.0_real_8
    ! ...end update coupling
    CALL derivov_new(a,adot,d,e,s,sdot)
    DO i=1,3
       a(i)=a(i)-k3(i)
       k4(i)=h*adot(i)
       a(i)=a(i)+k1(i)/6.0_real_8+k2(i)/3.0_real_8+k3(i)/3.0_real_8+k4(i)/6.0_real_8
    ENDDO
    s=s-ds3
    ds4=h*sdot
    s=s+ds1/6.0_real_8+ds2/3.0_real_8+ds3/3.0_real_8+ds4/6.0_real_8
    ! CONVERT BACK TO X
    x(1)=REAL(a(1))
    x(3)=AIMAG(a(1))
    x(2)=REAL(a(2))
    x(4)=AIMAG(a(2))
    x(5)=REAL(a(3))
    x(6)=AIMAG(a(3))
    p(1)=EXP(-ui*x(5))
    p(2)=EXP(-ui*x(6))
    DO i=1,2
       pdot(i)=-p(i)*ui*e(i)
       b(i)=a(i)*p(i)
       bdot(i)=adot(i)*p(i)+a(i)*pdot(i)
    ENDDO
    oms=SQRT(1.0_real_8-s**2.0_real_8)
    d1=b(1)+s*b(2)
    d2=b(2)*oms
    ! elisa added prob1
    prob1%d1sq=CONJG(d1)*d1
    prob1%d2sq=CONJG(d2)*d2
    ! COMPUTE TRANSITION PROBABILITIES
    sdot=d(1,2)+d(2,1)
    d1dot=bdot(1)+s*bdot(2)+sdot*b(2)
    d2dot=bdot(2)*oms-s*sdot/oms*b(2)
    prob1%d11dot=CONJG(d1dot)*d1+CONJG(d1)*d1dot
    prob1%d22dot=CONJG(d2dot)*d2+CONJG(d2)*d2dot
    ! 11   FORMAT(A3,I,F12.6,2X,F12.6,F12.6)
    ! 12   format(i6,4(1x,f12.8))
    ! McB
    ! if (parent) then
    ! write (6,*) 'rk4ov_fin: ', ' b(1) = ',b(1)
    ! write (6,*) 'rk4ov_fin: ', ' b(2) = ',b(2)
    ! write (6,*) 'rk4ov_fin: ', ' s = ',s 
    ! write (6,*) 'rk4ov_fin: ', ' oms = ',oms 
    ! write (6,*) 'rk4ov_fin: ', ' d1sq = ',d1sq
    ! write (6,*) 'rk4ov_fin: ', ' d2sq = ',d2sq
    ! write (6,*) 'rk4ov_fin: ', ' d1dot = ',d1dot
    ! write (6,*) 'rk4ov_fin: ', ' d2dot = ',d2dot
    ! endif
    ! McB
    RETURN
  END SUBROUTINE rk4ov_new
  ! ----------------------------------------------------------------------
  SUBROUTINE derivov_new(a,adot,d,e,s,sdot)
    COMPLEX(real_8)                          :: a(3), adot(3)
    REAL(real_8)                             :: d(sh02%nsurf,sh02%nsurf), &
                                                e(sh02%nsurf), s, sdot

    COMPLEX(real_8)                          :: p(2), p12, p21, ui
    INTEGER                                  :: n, n3
    REAL(real_8)                             :: de, s2m1

    ui=(0.0_real_8,1.0_real_8)
    n=sh02%nsurf
    n3=3*sh02%nsurf
    ! a(1)=(x(1),x(3))
    ! a(2)=(x(2),x(4))
    ! a(3)=(x(5),x(6))
    ! x(1)=x1
    ! x(2)=x2
    ! x(3)=y1
    ! x(4)=y2
    ! x(5)=eps(1)=int(e1)dt
    ! x(6)=eps(2)=int(e2)dt
    de=e(1)-e(2)
    p(1)=EXP(-ui*REAL(a(3)))
    p(2)=EXP(-ui*AIMAG(a(3)))
    p12=p(1)/p(2)
    p21=p(2)/p(1)
    s2m1=s**2.0_real_8 - 1.0_real_8
    adot(1)=(ui*a(2)*p21*s*de+a(2)*d(1,2)*p21-a(1)*d(2,1)*s)/s2m1
    adot(2)=(a(1)*d(2,1)*p12-a(2)*d(1,2)*s-ui*a(2)*s**2*de)/s2m1
    adot(3)=CMPLX(e(1),e(2),kind=real_8)
    sdot=d(1,2)+d(2,1)
    RETURN
  END SUBROUTINE derivov_new

  ! 
  ! #####################################################################
  ! 
  ! McB
  ! McB   OLD real(8) :: -mapped version of rk4ov
  ! McB
  SUBROUTINE rk4ov_old(h,x,d,dd,e,dei,s,ds)
    ! McB
    ! McB
    REAL(real_8) :: h, x(3*sh02%nsurf), d(sh02%nsurf,sh02%nsurf), &
      dd(sh02%nsurf,sh02%nsurf), e(sh02%nsurf), dei(sh02%nsurf), s, ds

    INTEGER                                  :: i, ii, n3, nzaehl
    REAL(real_8) :: csq(sh02%nsurf), d1, d2, de, dtot, k1(3*sh02%nsurf), &
      k2(3*sh02%nsurf), k3(3*sh02%nsurf), k4(3*sh02%nsurf), oms, sdot, &
      xdot(3*sh02%nsurf)

    nzaehl=0

    n3=3*sh02%nsurf
    nzaehl=nzaehl+1
    CALL derivov_old(x,xdot,d,e,s)
    DO i=1,n3
       k1(i)=h*xdot(i)
       x(i)=x(i)+k1(i)/2.0
    ENDDO
    CALL derivov_old(x,xdot,d,e,s)
    DO i=1,n3
       x(i)=x(i)-k1(i)/2.0_real_8
       k2(i)=h*xdot(i)
       x(i)=x(i)+k2(i)/2.0_real_8
    ENDDO
    CALL derivov_old(x,xdot,d,e,s)
    DO i=1,n3
       x(i)=x(i)-k2(i)/2.0_real_8
       k3(i)=h*xdot(i)
       x(i)=x(i)+k3(i)
    ENDDO
    CALL derivov_old(x,xdot,d,e,s)
    DO i=1,n3
       x(i)=x(i)-k3(i)
       k4(i)=h*xdot(i)
       x(i)=x(i)+k1(i)/6.0_real_8+k2(i)/3.0_real_8+k3(i)/3.0_real_8+k4(i)/6.0_real_8
    ENDDO
    DO i=1,sh02%nsurf
       ii=i+sh02%nsurf
       csq(i)=x(i)**2.0_real_8+x(ii)**2.0_real_8
       ! if (mod(nzaehl,40).eq.1) then
       ! write(6,11)'rk4',i,csq(i),x(i),x(ii)
       ! endif
    ENDDO
    de=x(5)-x(6)
    dtot=2.0_real_8*s*(COS(de)*(x(1)*x(2)+x(3)*x(4))+SIN(de)&
         *(x(3)*x(2)-x(1)*x(4)))
    d1=csq(1)+csq(2)*s**2.0_real_8+dtot
    d2=csq(2)*(1.0_real_8-s**2.0_real_8)
    ! McB
    prob1%d1sq=d1
    prob1%d2sq=d2
    ! McB
    ! COMPUTE TRANSITION PROBABILITIES
    sdot=d(1,2)+d(2,1)
    oms=SQRT(1.0_real_8-s**2.0_real_8)
    prob1%d22dot=2.0_real_8*((x(2)*xdot(2)+x(4)*xdot(4))*oms**2-csq(2)*s*sdot)
    ! WRITE(6,*)'D22DOT [2]=',D22DOT
    prob1%d11dot=2.0_real_8*(x(1)*xdot(1)+x(3)*xdot(3)+s**2*(x(2)*xdot(2)+&
         X(4)*XDOT(4))+S*SDOT*CSQ(2)+SDOT*0.5_real_8*DTOT/S+S*(COS(DE)*&
         (XDOT(1)*X(2)+X(1)*XDOT(2)+XDOT(3)*X(4)+X(3)*XDOT(4)+&
         (E(1)-E(2))*(X(3)*X(2)-X(1)*X(4)))+SIN(DE)*(XDOT(3)*X(2)+&
         X(3)*XDOT(2)-XDOT(1)*X(4)-X(1)*XDOT(4)-(E(1)-E(2))*(X(1)*X(2)&
         +X(3)*X(4)))))
    ! WRITE(6,*)'D11DOT [2]=',D11DOT
    ! if (mod(nzaehl,40).eq.1) then
    ! WRITE(6,*)'tot',csq(1)+csq(2)+dtot
    ! WRITE(6,*)'dtot',dtot
    ! WRITE(6,*)'POPULATION 1',d1
    ! WRITE(6,*)'POPULATION 2',d2
    ! ---- WRITE COEFFICIENTS TO FILE --------------------------
    ! open(789,file='prob.dat',position='append')
    ! write(789,12)infi,csq(1),csq(2),d1,d2
    ! close(789)
    ! --------------------------------------------------------
    ! endif
11  FORMAT(a3,i6,f12.6,2x,f12.6,f12.6)
12  FORMAT(i6,4(1x,f12.8))
    RETURN
  END SUBROUTINE rk4ov_old
  ! ----------------------------------------------------------------------
  SUBROUTINE derivov_old(x,xdot,d,e,s)
    REAL(real_8)                             :: x(3*sh02%nsurf), &
                                                xdot(3*sh02%nsurf), &
                                                d(sh02%nsurf,sh02%nsurf), &
                                                e(sh02%nsurf), s

    INTEGER                                  :: n, n3
    REAL(real_8)                             :: de

    n=sh02%nsurf
    n3=3*sh02%nsurf
    de=x(5)-x(6)

    ! x(1)=x1
    ! x(2)=x2
    ! x(3)=y1
    ! x(4)=y2
    ! x(5)=eps(1)=int(e1)dt
    ! x(6)=eps(2)=int(e2)dt
    ! ---- i=1 --------------------------------------------------------------
    xdot(1)=(x(1)*d(2,1)*s-x(2)*d(1,2)*COS(de)+x(4)*d(1,2)*SIN(de)&
         +s*(e(1)-e(2))*(x(2)*SIN(de)+x(4)*COS(de))&
         )/(1.0_real_8-s**2.0_real_8) 

    xdot(3)=(x(3)*d(2,1)*s-x(2)*d(1,2)*SIN(de)-x(4)*d(1,2)*COS(de)&
         +s*(e(1)-e(2))*(x(4)*SIN(de)-x(2)*COS(de))&
         )/(1.0_real_8-s**2.0_real_8) 

    xdot(5)=e(1)
    ! ----i=2 ---------------------------------------------------------------
    xdot(2)=(x(2)*d(1,2)*s-x(1)*d(2,1)*COS(-de)+x(3)*d(2,1)*SIN(-de)&
         -s**2*(e(1)-e(2))*x(4))/(1.0_real_8-s**2.0_real_8) 

    xdot(4)=(x(4)*d(1,2)*s-x(1)*d(2,1)*SIN(-de)-x(3)*d(2,1)*COS(-de)&
         +s**2*(e(1)-e(2))*x(2))/(1.0_real_8-s**2.0_real_8) 

    xdot(6)=e(2)
    RETURN
  END SUBROUTINE derivov_old

END MODULE rk4ov_utils
