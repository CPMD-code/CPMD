MODULE legendre_p_utils
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: paral

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: d_binom
  PUBLIC :: rsphgen

CONTAINS

  ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  ! C
  ! Computes the real spherical harmonic of order lm           C
  ! Here |m|<=L, L>=0                                          C
  ! Cohen-Tannoudji convention: Y0L (0) = SQRT(2L+1/4Pi)       C
  ! C
  ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  REAL(real_8) FUNCTION rsphgen(l,m,Vx,Vy,Vz)
    IMPLICIT NONE

    INTEGER :: l,m
    REAL(real_8) :: Vx,Vy,Vz

    ! including the constants(PI) definition

    REAL(real_8) :: Prefactor
    REAL(real_8) :: RESULT
    REAL(real_8) :: r_xyz,ctheta,phi
    REAL(real_8) :: psi_of_M,plm

    CALL stopgm('LEGENDRE_P','RSPHGEN not parallel.',& 
         __LINE__,__FILE__)

    ! First calculating the norm of the vector
    r_xyz=SQRT(Vx**2+Vy**2+Vz**2)

    ! is that a zero vector?
    IF (r_xyz.EQ.0._real_8) THEN
       ! YES!
       IF (m.EQ.0) THEN
          ! same value for all L, exit
          ! RSPHGEN=1._real_8
          rsphgen=0._real_8
          RETURN
       ELSE
          ! otherwise, zero...
          rsphgen=0._real_8
          RETURN
       ENDIF
    ENDIF

    ! From here onwards, R>0
    ! calculating cos(theta) [z=R cos(theta)]
    ctheta=Vz/r_xyz

    ! special case if theta=0 (Merzbacher p251)
    IF (ABS(ctheta).EQ.1._real_8) THEN
       IF (m.EQ.0) THEN
          rsphgen=1._real_8
          RETURN
       ELSE
          rsphgen=0._real_8
          RETURN
       ENDIF
    ENDIF

    ! calculating psi of m 
    ! m=0
    IF (m.EQ.0) THEN
       psi_of_M=1._real_8
       Prefactor=1._real_8
       plm=plgndr(l,0,ctheta)
    ENDIF

    ! m>0
    IF (m.GT.0) THEN
       ! calculating phi [x=R sin(theta) cos(phi), y=R sin(theta) sin(phi)]
       phi=datan2(Vy,Vx)

       psi_of_M=SQRT(2._real_8)*COS(REAL(m,kind=real_8)*phi)
       Prefactor=(-1)**m *SQRT(d_fac_rat(l-m,l+m))
       plm=plgndr(l,m,ctheta)
    ENDIF

    ! m<0
    IF (m.LT.0) THEN
       ! calculating phi
       phi=datan2(Vy,Vx)

       psi_of_M=SQRT(2._real_8)*SIN(REAL(-m,kind=real_8)*phi)
       Prefactor=(-1)**m *SQRT(d_fac_rat(l+m,l-m))
       plm=plgndr(l,-m,ctheta)
    ENDIF

    ! Now we can calculate the value of the spherical harmonic...

    RESULT=Prefactor*plm*psi_of_M

    ! copying the result
    rsphgen=RESULT

    ! IF the resulting value is a NAN, output debug informations
    IF (RESULT.NE.RESULT) THEN
       IF (paral%io_parent)&
            WRITE(6,*) '-SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS-'
       IF (paral%io_parent)&
            WRITE(6,*) 'legendre: L=',l,'M=',m,&
            'V(',Vx,',',Vy,',',Vz,')'

       IF (m.NE.0) THEN
          IF (paral%io_parent)&
               WRITE(6,*) 'phi =',phi
       ENDIF
       IF (paral%io_parent)&
            WRITE(6,*) 'ctheta =',ctheta
       IF (paral%io_parent)&
            WRITE(6,*) 'prefac =',Prefactor
       IF (paral%io_parent)&
            WRITE(6,*) 'psiofM=',psi_of_M
       IF (paral%io_parent)&
            WRITE(6,*) 'Result = ',RESULT
       IF (paral%io_parent)&
            WRITE(6,*) 'L-M! =',d_fac(l-m)
       IF (paral%io_parent)&
            WRITE(6,*) 'L+M! =',d_fac(l+m)
       IF (paral%io_parent)&
            WRITE(6,*) 'D_FAC_RAT(L+M,L-M)=',d_fac_rat(l+m,l-m)
       IF (paral%io_parent)&
            WRITE(6,*) 'D_FAC_RAT(L-M,L+M)=',d_fac_rat(l-m,l+m)
       IF (paral%io_parent)&
            WRITE(6,*) '-EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE-'
    ENDIF

    RETURN
  END FUNCTION rsphgen


  ! calculates the factorial of a returns real
  REAL(real_8) FUNCTION d_fac(a)
    IMPLICIT NONE
    INTEGER :: a

    INTEGER :: i
    REAL(real_8) :: resl

    resl=1._real_8
    IF (a.GT.1) THEN
       DO i=1,a
          resl=resl*REAL(i,kind=real_8)
       ENDDO
    ENDIF

    d_fac=resl

    RETURN
  END FUNCTION d_fac

  ! calculates the binomial coefficient n \over k
  REAL(real_8) FUNCTION d_binom(n,k)
    IMPLICIT NONE
    INTEGER :: n,k

    d_binom=0.0_real_8
    IF ((k.GE.0).AND.(k.LE.n)) THEN
       d_binom=d_fac(n)/(d_fac(n-k)*d_fac(k))
    ENDIF
    RETURN
  END FUNCTION d_binom

  ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  ! C
  ! calculates the ratio of (a!)/(b!)     C
  ! and returns the result as a real      C
  ! C
  ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  REAL(real_8) FUNCTION d_fac_rat(a,b)
    IMPLICIT NONE
    INTEGER :: a,b

    INTEGER :: i
    REAL(real_8) :: resl

    resl=1._real_8

    IF (a.LT.0)THEN
       IF (paral%io_parent)&
            WRITE(6,*) 'error in D_FAC_RAT, negative A'
    ENDIF
    IF (b.LT.0)THEN
       IF (paral%io_parent)&
            WRITE(6,*) 'error in D_FAC_RAT, negative B'
    ENDIF

    ! A>B then we only need to calculate A*(A-1)*...*(B+1)
    IF (a.GT.b) THEN
       DO i=(b+1),a
          resl=resl*REAL(i,kind=real_8)
       ENDDO
       ! B<A then we need to calculate B*(B-1)*...*(A+1) and invert it
    ELSEIF (b.GT.a) THEN
       DO i=(a+1),b
          resl=resl*REAL(i,kind=real_8)
       ENDDO
       resl=1.0_real_8/resl
    ENDIF

    d_fac_rat=resl
    RETURN
  END FUNCTION d_fac_rat

  ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  ! C
  ! Computes the associated Legendre polynomial Pml (x).       C
  ! Here m and l are integers satisfying 0 <= m <= l,          C
  ! while x lies in the range -1 <= x <= 1                     C
  ! Note:                                                      C
  ! Numerical recipe subroutine, changed to double precision   C
  ! C
  ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  REAL(real_8) FUNCTION plgndr(l,m,x)
    IMPLICIT NONE

    INTEGER :: l,m
    REAL(real_8) :: x

    INTEGER :: i,ll
    REAL(real_8) :: fact,pll,pmm,pmmp1,somx2
    CALL stopgm('LEGENDRE_P','PLGNDR not parallel.',& 
         __LINE__,__FILE__)

    IF (m.LT.0.OR.m.GT.l.OR.ABS(x).GT.1._real_8) THEN
       CALL stopgm('LEGENDRE_P','bad arguments in plgndr',& 
            __LINE__,__FILE__)
    ENDIF

    ! Compute Pmm
    pmm=1._real_8

    IF (m.GT.0) THEN
       somx2=SQRT((1._real_8-x)*(1._real_8+x))
       fact=1._real_8
       DO i=1,m
          pmm=-pmm*fact*somx2
          fact=fact+2._real_8
       ENDDO
    ENDIF

    IF (l.EQ.m) THEN
       plgndr=pmm
    ELSE

       ! Compute P mm+1
       pmmp1=x*REAL(2*m+1,kind=real_8)*pmm

       IF (l.EQ.m+1) THEN
          plgndr=pmmp1
       ELSE

          ! Compute P ml, l >m + 1
          DO ll=m+2,l
             pll=(x*REAL(2*ll-1,kind=real_8)*pmmp1-REAL(ll+m-1,kind=real_8)*pmm)&
                  /REAL(ll-m,kind=real_8)
             pmm=pmmp1
             pmmp1=pll
          ENDDO
          plgndr=pll
       ENDIF
    ENDIF
    RETURN
  END FUNCTION plgndr

END MODULE legendre_p_utils
