MODULE ghermit_utils
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ghermit

CONTAINS

  SUBROUTINE ghermit(ngh,xgh,wgh)
    INTEGER                                  :: ngh
    REAL(real_8)                             :: xgh(*), wgh(*)

    INTEGER, PARAMETER                       :: nmax = 180

    INTEGER                                  :: i, index, m, n, niter, np, nup
    REAL(real_8) :: a(0:nmax), ap(0:nmax), da(0:nmax), dy, fact, pi, r, rapp, &
      root(nmax), root0(nmax), sqpi, w, x, y, y0, yp, yp1

! input

    n=2*ngh
    ! compute rhe coefs.
    IF (n.GT.nmax)&
         CALL stopgm('GHERMIT',' NUMBER OF GH-POINTS',& 
         __LINE__,__FILE__)
    nup=n/2
    DO i=0,n
       a(i)=0.0000000000000000000000000000000_real_8
       da(i)=0.0000000000000000000000000000000_real_8
       ap(i)=0.0000000000000000000000000000000_real_8
    ENDDO
    fact=nup+1
    DO i=nup+2,n
       fact=fact*i
    ENDDO
    a(0)=(-1.0000000000000000000000_real_8)**nup*fact
    fact=1.00000000000000000000000000000000_real_8
    r=0.00000000000000000000000000_real_8
    DO i=1,n
       r=r+1.00000000000000000000000_real_8
       fact=fact*SQRT(r)
    ENDDO
    pi=ACOS(-1.0000000000000000_real_8)
    sqpi=SQRT(pi)
    sqpi=SQRT(sqpi)
    a(0)=a(0)/(2.000000000000_real_8**nup)/fact/sqpi
    DO i=nup-1,0,-1
       m=i+1
       index=n-2*i
       a(index)=-a(index-2)*4.0000000000000000000000000_real_8&
            *REAL(m,kind=real_8)/&
            REAL(n-m-m+2,kind=real_8)/REAL(n-m-m+1,kind=real_8)
    ENDDO
    ! derivative
    DO  i=0,n-1
       da(i)=(i+1.0000000000000000000000_real_8)*a(i+1)
    ENDDO
    np=n+1
    nup=np/2
    fact=nup+1
    DO i=nup+2,np
       fact=fact*i
    ENDDO
    ap(1)=2.000000000000000000000_real_8*&
         (-1.000000000000000000000_real_8)**nup*fact
    fact=1.0000000000000000000000_real_8
    r=0.00000000000000000000000000_real_8
    DO i=1,np
       r=r+1.00000000000000000000000_real_8
       fact=fact*SQRT(r)
    ENDDO
    ap(1)=ap(1)/(2.000000000000_real_8**nup)/fact/sqpi&
         /SQRT(2.000000000000_real_8)
    DO i=nup-1,0,-1
       m=i+1
       index=np-2*i
       ap(index)=-ap(index-2)*4.000000000000000000_real_8&
            *REAL(m,kind=real_8)/&
            REAL(np-m-m+2,kind=real_8)/REAL(np-m-m+1,kind=real_8)
    ENDDO
    ! first approximation to roots
    CALL eval(a,n,0._real_8,y0)
    index=0
    DO i=1,5000
       x=0.00333_real_8*i
       CALL eval(a,n,x,y)
       IF (y*y0.LT.0._real_8) THEN
          index=index+1
          IF (index.GT.nmax) CALL stopgm('GHERMIT','too many roots',& 
               __LINE__,__FILE__)
          root0(index)=x
          y0=y
       ENDIF
    ENDDO

    ! newton search for the roots
    DO i=1,index
       niter=0
18     CALL eval(a,n,root0(i),y)
       IF (ABS(y).GT.1.e-16_real_8) THEN
          niter=niter+1
          IF (niter.GT.20)go to 19
          CALL eval(da,n-1,root0(i),yp)
          root0(i)=root0(i)-y/yp
          go to 18
       ENDIF
19     root(i)=root0(i)
       xgh(i)=root(i)
    ENDDO
    ! weights
    rapp=-ap(n+1)/a(n)
    DO i=1,index
       CALL eval(ap,n+1,root(i),yp1)
       CALL eval(da,n-1,root(i),dy)
       w=rapp/yp1/dy
       wgh(i)=w
    ENDDO
    RETURN
  END SUBROUTINE ghermit
  SUBROUTINE eval(coef,n,x,y)
    INTEGER                                  :: n
    REAL(real_8)                             :: coef(0:n), x, y

    INTEGER                                  :: i

    y=coef(n)
    DO i=n-1,0,-1
       y=coef(i)+x*y
    ENDDO
    RETURN
  END SUBROUTINE eval

END MODULE ghermit_utils
