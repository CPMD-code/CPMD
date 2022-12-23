#ifdef __SR8000
!option MP(P(0)), LANGLVL(SAVE(0))
#endif

MODULE bessm_utils
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: paral

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: getf
  PUBLIC :: bessov
  PUBLIC :: bessl

CONTAINS

  ! 
  ! .....................................................GETF.........
  SUBROUTINE getf(nr,f,x,res)
    ! **                                                              **
    ! **  CALCULATES THE VALUE RES OF THE FUNCTION F AT X             **
    ! **  BY INTERPOLATION BY A THIRD ORDER POLYNOM                   **
    ! **  F(X=I)=F(I)                                                 **
    ! **                                                              **
    ! **          P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1991) **
    ! **                                                              **
    INTEGER                                  :: nr
    REAL(real_8), DIMENSION(nr)              :: f
    REAL(real_8)                             :: x, res

    INTEGER                                  :: incr
    REAL(real_8)                             :: p1, p2, p21, p3, p32, p321, &
                                                p4, p43, p432, p4321, xx1, &
                                                xx2, xx3, xx4

    incr=INT(x)
    incr=MIN(incr,nr-2)
    incr=MAX(incr,2)
    xx1=x-REAL(incr-1,kind=real_8)
    xx2=xx1-1._real_8
    xx3=xx2-1._real_8
    xx4=xx3-1._real_8
    p1=f(incr-1)
    p2=f(incr)
    p3=f(incr+1)
    p4=f(incr+2)
    p21=-xx2*p1+xx1*p2
    p32=-xx3*p2+xx2*p3
    p43=-xx4*p3+xx3*p4
    p321=0.5_real_8*( -xx3*p21+xx1*p32 )
    p432=0.5_real_8*( -xx4*p32+xx2*p43 )
    p4321=1._real_8/3._real_8 * ( -xx4*p321+xx1*p432 )
    res=p4321
    RETURN
  END SUBROUTINE getf
  ! 
  ! .....................................................BESSOV.......
  SUBROUTINE bessov(r1,dex,nr,f,l,g,rmax,res)
    ! **                                                              **
    ! **  CALCULATES THE OVERLAP BETWEEN A SPHERICAL BESSEL FUNCTION  **
    ! **  WITH THE FUNCTION F GIVEN ON AN LOGARITHMIC GRID            **
    ! **  (3-D OVERLAPP)                                              **
    ! **                                                              **
    ! **          P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1991) **
    ! **                                                              **
    REAL(real_8)                             :: r1, dex
    INTEGER                                  :: nr
    REAL(real_8), DIMENSION(nr)              :: f
    INTEGER                                  :: l
    REAL(real_8)                             :: g, rmax, res

    INTEGER                                  :: i, ir, nstep
    REAL(real_8)                             :: r, rstep, val, x
    REAL(real_8), DIMENSION(4)               :: fac

    nstep=INT(20._real_8*g/6._real_8)
    nstep=MAX(nstep,400)
    rstep=rmax/nstep
    fac(1)=rstep*17._real_8/48._real_8
    fac(2)=rstep*59._real_8/48._real_8
    fac(3)=rstep*43._real_8/48._real_8
    fac(4)=rstep*49._real_8/48._real_8
    x=1._real_8+LOG(rmax/r1) / dex
    CALL getf(nr,f,x,val)
    res=fac(1)*val*bessl(l,rmax*g)*rmax**2
    DO i=2,4
       r=rstep*REAL(i-1,kind=real_8)
       x=1._real_8+LOG(r/r1) / dex
       CALL getf(nr,f,x,val)
       res=res+fac(i)*val*bessl(l,r*g)*r**2
       r=rmax-r
       x=1._real_8+LOG(r/r1) / dex
       CALL getf(nr,f,x,val)
       res=res+fac(i)*val*bessl(l,r*g)*r**2
    ENDDO
    DO ir=5,nstep-3
       r=rstep*REAL(ir-1,kind=real_8)
       x=1._real_8+LOG(r/r1) / dex
       CALL getf(nr,f,x,val)
       res=res+rstep*val*bessl(l,r*g)*r**2
    ENDDO
    RETURN
  END SUBROUTINE bessov
  ! 
  ! ..............................................BESSL...............
  FUNCTION bessl(l,x) RESULT(res)
    ! **                                                              **
    ! **  CALCULATES THE SPHERICAL BESSEL FUNCTION                    **
    ! **  ACCORDING TO ABRAMOWITZ AND STEGUN                          **
    ! **    FORMULA 10.1.2              FOR   X < 8                   **
    ! **    FORMULA 10.1.8 AND  10.1.9  FOR   X > 8                   **
    ! **                                                              **
    ! **          P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1991) **
    ! **          
    INTEGER                                  :: l
    REAL(real_8)                             :: x, res

    REAL(real_8), PARAMETER                  :: tol = 1.0e-10_real_8

    INTEGER                                  :: i, ii, il, isvar, k
    REAL(real_8)                             :: arg, fac, pi, xsq
    REAL(real_8), DIMENSION(0:100)           :: facul
    REAL(real_8), DIMENSION(4)               :: trig

    IF (x.GT.REAL(l,kind=real_8)) THEN
       pi=4._real_8*ATAN(1._real_8)
       arg=x-0.5_real_8*REAL(l,kind=real_8)*pi
       trig(1)=SIN(arg)/x
       trig(2)=COS(arg)/x
       trig(3)=-trig(1)
       trig(4)=-trig(2)
       res=trig(1)
       IF (l.EQ.0) RETURN
       facul(0)=1._real_8
       DO i=1,2*l
          facul(i)=facul(i-1)*REAL(i,kind=real_8)
       ENDDO
       xsq=0.5_real_8/x
       fac=1._real_8
#ifdef _vpp_
       !OCL SCALAR 
#endif 
       DO k=1,l
          ii=MOD(k,4)+1
          fac=facul(k+l)/facul(k)/facul(l-k)*xsq**k
          ! FAC=FAC*XSQ*real(L+K,kind=real_8)/real(K*(L-K),kind=real_8)
          res=res+fac*trig(ii)
       ENDDO
       ! II=MOD(L,4)+1
       ! FAC=FAC*XSQ*real(2*L,kind=real_8)/real(L,kind=real_8)
       ! BESSL=BESSL+FAC*TRIG(II)
       RETURN
    ENDIF
    ! ==================================================================
    ! ==  TAYLOR EXPANSION FOR SMALL ARGUMENTS                        ==
    ! ==================================================================
    isvar=1
    DO il=1,l
       isvar=isvar*(2*il+1)
    ENDDO
    IF (l.NE.0._real_8) THEN
       fac=x**l/REAL(isvar,kind=real_8)
    ELSE
       fac=1._real_8/REAL(isvar,kind=real_8)
    ENDIF
    res=fac
    xsq=-0.5_real_8*x*x
    isvar=2*l+1
    DO i=1,1000
       isvar=isvar+2
       fac=fac*xsq/REAL(i*isvar,kind=real_8)
       res=res+fac
       IF (ABS(fac).LT.tol) GOTO 9999
    ENDDO
    IF (paral%io_parent)&
         WRITE(6,*)'BESSL NOT CONVERGED X=',x,' L= ',l
    CALL stopgm ('BESSL',' ',& 
         __LINE__,__FILE__)
9999 CONTINUE
    RETURN
  END FUNCTION bessl

END MODULE bessm_utils
