MODULE mtin_utils
  USE cnst,                            ONLY: fpi,&
                                             pi
  USE cppt,                            ONLY: gk,&
                                             hg,&
                                             nzh,&
                                             scg
  USE error_handling,                  ONLY: stopgm
  USE fft_maxfft,                      ONLY: maxfft
  USE fftmain_utils,                   ONLY: fwfftn
  USE geq0mod,                         ONLY: geq0
  USE isos,                            ONLY: isos1,&
                                             isos2,&
                                             isos3
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: parai,&
                                             paral
  USE prmem_utils,                     ONLY: prmem
  USE special_functions,               ONLY: cp_erf
  USE system,                          ONLY: fpar,&
                                             ncpw,&
                                             parap,&
                                             parm,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mtin
  PUBLIC :: moin

CONTAINS

  ! ==================================================================
  SUBROUTINE mtin
    ! ==--------------------------------------------------------------==
    ! == SCREENING FUNCTION FOR TUCKERMAN POISSON SOLVER              ==
    ! ==--------------------------------------------------------------==
    CHARACTER(*), PARAMETER                  :: procedureN = 'mtin'
    REAL(real_8), PARAMETER                  :: eeps = 1.e-20_real_8 

    COMPLEX(real_8), ALLOCATABLE             :: vr(:)
    INTEGER                                  :: ierr

    ALLOCATE(vr(maxfft),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(scg(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    IF (isos1%toned.OR.isos1%ttwod) CALL stopgm('MTIN','not implemented',& 
         __LINE__,__FILE__)
    ! ==------------------------------------------------------------==
    CALL screena(eeps)
    CALL screenb(vr,eeps)
    DEALLOCATE(vr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (paral%io_parent) CALL prmem('   CLUSTER')
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE mtin
  ! ==================================================================
  SUBROUTINE moin
    ! ==--------------------------------------------------------------==
    ! == SCREENING FUNCTION FOR MORTENSEN POISSON SOLVER              ==
    ! ==--------------------------------------------------------------==
    CHARACTER(*), PARAMETER                  :: procedureN = 'moin'

    INTEGER                                  :: ierr, ig, ig1, nz
    REAL(real_8)                             :: alpha, bj0, bj1, bk0, bk1, &
                                                g2, gx, gxy, gy, gyz, gz, z

! ==--------------------------------------------------------------==

    ALLOCATE(scg(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    IF (isos1%toned) THEN
       ! ..we assume periodicity in x, square box in y and z, no total charge
       z=parm%a3(3)*0.5_real_8
       ig1=1
       IF (geq0) ig1=2
       DO ig=ig1,ncpw%nhg
          g2=parm%tpiba2*hg(ig)
          gx=gk(1,ig)*parm%tpiba
          gy=gk(2,ig)
          gz=gk(3,ig)
          gyz=parm%tpiba*SQRT(gz*gz+gy*gy)
          bj0=bessj0(z*gyz)
          bj1=bessj1(z*gyz)
          bk0=bessk0(z*gx)
          bk1=bessk1(z*gx)
          scg(ig)=fpi/g2 * (1._real_8 + z*(gyz*bj1*bk0-gx*bj0*bk1))
       ENDDO
       IF (geq0) scg(1)=0._real_8
    ELSEIF (isos1%ttwod) THEN
       ! ..we assume orthorhobmic box with surface in specified plane, no total charge
       IF (isos3%snormal.EQ.1) THEN
          z=parm%a1(1)
       ELSE IF (isos3%snormal.EQ.2) THEN
          z=parm%a2(2)
       ELSE IF (isos3%snormal.EQ.3) THEN
          z=parm%a3(3)
       ENDIF
       ig1=1
       IF (geq0) ig1=2
       DO ig=ig1,ncpw%nhg
          g2=parm%tpiba2*hg(ig)
          gx=gk(1,ig)
          gy=gk(2,ig)
          gz=gk(3,ig)
          IF (isos3%snormal.EQ.1) THEN       ! YZ PLANE
             gxy=parm%tpiba*SQRT(gy*gy+gz*gz)
             nz=NINT(z/parm%alat*gx)
          ELSE IF (isos3%snormal.EQ.2) THEN  ! ZX PLANE
             gxy=parm%tpiba*SQRT(gz*gz+gx*gx)
             nz=NINT(z/parm%alat*gy)
          ELSE IF (isos3%snormal.EQ.3) THEN  ! XY PLANE
             gxy=parm%tpiba*SQRT(gx*gx+gy*gy)
             nz=NINT(z/parm%alat*gz)
          ENDIF
          scg(ig)=fpi/g2 * (1._real_8 - (-1._real_8)**nz * EXP(-0.5_real_8*z*gxy))
          ! ..following lines would be rigorous: uncomment in case of problems
          ! GZ=TPIBA*GZ
          ! SCG(IG)=FPI/G2 * (1._real_8 - cos(0.5_real_8*Z*GZ)*exp(-0.5*Z*GXY))
       ENDDO
       IF (geq0) scg(1)=0._real_8
    ELSE
       alpha=parm%alat*0.5_real_8
       ig1=1
       IF (geq0) ig1=2
       DO ig=ig1,ncpw%nhg
          g2=parm%tpiba2*hg(ig)
          scg(ig)=fpi/g2 * (1._real_8 - COS(alpha*SQRT(g2)))
       ENDDO
       IF (geq0) scg(1)=0.5_real_8*fpi*alpha*alpha
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE moin
  ! ==================================================================
  SUBROUTINE screena(eeps)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: eeps

    INTEGER                                  :: ig, ig1
    REAL(real_8)                             :: alpha, beta, g, g2, rmin, t0, &
                                                t1, t2, t3

    rmin=MIN(parm%a1(1),parm%a2(2),parm%a3(3))
    alpha=isos2%alphal/rmin
    beta=SQRT(-LOG(eeps))/rmin
    beta=MAX(0.2_real_8,beta)
    t0=2._real_8*pi*alpha/(beta*beta)/SQRT(alpha*alpha+beta*beta)
    t1=2._real_8*pi**(1.5_real_8)/beta
    t2=-0.25_real_8/(beta*beta)
    t3=SQRT(alpha*alpha/(alpha*alpha+beta*beta))*0.5_real_8/beta
    ig1=1
    IF (geq0) ig1=2
    DO ig=ig1,ncpw%nhg
       g2=parm%tpiba2*hg(ig)
       g=SQRT(g2)
       scg(ig)=t1/g*erfi(g*t3,t2*g2)
    ENDDO
    IF (geq0) scg(1)=t0
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE screena
  ! ==================================================================
  SUBROUTINE screenb(vr_,eeps)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8), TARGET                  :: vr_(:)
    REAL(real_8)                             :: eeps

    COMPLEX(real_8), DIMENSION(:, :, :), &
      POINTER                                :: vr
    INTEGER                                  :: i, ig, ig1, ii, isub, j, k
    REAL(real_8)                             :: alpha, ax, ay, az, beta, dx, &
                                                dy, dz, fb, g2, oa, op, r, &
                                                r2, rmin, x, y, z

    CALL tiset('   SCREENB',isub)

    vr(1:fpar%kr1,1:fpar%kr2s,1:fpar%kr3s) => vr_

    CALL zeroing(vr)!,maxfft)
    rmin=MIN(parm%a1(1),parm%a2(2),parm%a3(3))
    alpha=isos2%alphal/rmin
    beta=SQRT(-LOG(eeps))/rmin
    beta=MAX(0.2_real_8,beta)
    dx=parm%a1(1)/REAL(spar%nr1s,kind=real_8)
    dy=parm%a2(2)/REAL(spar%nr2s,kind=real_8)
    dz=parm%a3(3)/REAL(spar%nr3s,kind=real_8)
    ax=parm%a1(1)*0.5_real_8
    ay=parm%a2(2)*0.5_real_8
    az=parm%a3(3)*0.5_real_8
    op=parm%omega
    DO k=1,spar%nr3s
       z=REAL(k-1,kind=real_8)*dz
       IF (z.GT.az) z=2._real_8*az-z
       DO j=1,spar%nr2s
          y=REAL(j-1,kind=real_8)*dy
          IF (y.GT.ay) y=2._real_8*ay-y
          DO ii=1,parm%nr1
             i=ii+parap%nrxpl(parai%mepos,1)-1
             x=REAL(i-1,kind=real_8)*dx
             IF (x.GT.ax) x=2._real_8*ax-x
             r2=x*x+y*y+z*z
             r=SQRT(r2)
             fb=1._real_8-EXP(-beta*beta*r2)
             IF (r.GT.1.e-10_real_8) THEN
                vr(ii,j,k)=op*cp_erf(alpha*r)/r*fb
             ELSE
                vr(ii,j,k)=2._real_8*op*alpha/SQRT(pi)*fb
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    CALL fwfftn(vr_,.FALSE.,parai%allgrp)

    CALL hack_screenb(ncpw%nhg,vr,nzh,scg)
    ! 
    op=fpi
    oa=1._real_8/(4._real_8*alpha*alpha)
    ig1=1
    IF (geq0) ig1=2
    DO ig=ig1,ncpw%nhg
       g2=parm%tpiba2*hg(ig)
       scg(ig) = scg(ig) + fpi/g2 - op*EXP(-g2*oa)/g2
    ENDDO
    IF (geq0) scg(1)=scg(1) + op*oa
    CALL tihalt('   SCREENB',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE screenb

  SUBROUTINE hack_screenb(nhg,vr,nzh,scg)
    INTEGER                                  :: nhg
    COMPLEX(real_8)                          :: vr(*)
    INTEGER                                  :: nzh(nhg)
    COMPLEX(real_8)                          :: scg(nhg)

    INTEGER                                  :: ig

    DO ig=1,nhg
       scg(ig) = scg(ig) + vr(nzh(ig))
    ENDDO
  END SUBROUTINE hack_screenb
  ! ==================================================================
  FUNCTION erfi(x,y)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: x, y, erfi

    LOGICAL                                  :: flag
    REAL(real_8)                             :: u, v

! Variables
! complex(8) :: A
! ==--------------------------------------------------------------==

    CALL wofz(x,0._real_8,u,v,flag)
    ! A=CMPLX(0._real_8,1._real_8)*(1._real_8-EXP(X*X)*CMPLX(U,V))
    ! ERFI=real(A)
    erfi=v*EXP(x*x+y)
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION erfi
  ! ==================================================================
  ! ALGORITHM 680, COLLECTED ALGORITHMS FROM ACM.
  ! THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
  ! VOL. 16, NO. 1, PP. 47.
  SUBROUTINE wofz (xi, yi, u, v, flag)
    ! 
    ! GIVEN A COMPLEX NUMBER Z = (XI,YI), THIS SUBROUTINE COMPUTES
    ! THE VALUE OF THE FADDEEVA-FUNCTION W(Z) = EXP(-Z**2)*ERFC(-I*Z),
    ! WHERE ERFC IS THE COMPLEX COMPLEMENTARY ERROR-FUNCTION AND I
    ! MEANS SQRT(-1).
    ! THE ACCURACY OF THE ALGORITHM FOR Z IN THE 1ST AND 2ND QUADRANT
    ! IS 14 SIGNIFICANT DIGITS; IN THE 3RD AND 4TH IT IS 13 SIGNIFICANT
    ! DIGITS OUTSIDE A CIRCULAR REGION WITH RADIUS 0.126 AROUND A ZERO
    ! OF THE FUNCTION.
    ! ALL REAL VARIABLES IN THE PROGRAM ARE DOUBLE PRECISION.
    ! 
    ! 
    ! THE CODE CONTAINS A FEW COMPILER-DEPENDENT PARAMETERS :
    ! RMAXREAL = THE MAXIMUM VALUE OF RMAXREAL EQUALS THE ROOT OF
    ! RMAX = THE LARGEST NUMBER WHICH CAN STILL BE
    ! IMPLEMENTED ON THE COMPUTER IN DOUBLE PRECISION
    ! FLOATING-POINT ARITHMETIC
    ! RMAXEXP  = LN(RMAX) - LN(2)
    ! RMAXGONI = THE LARGEST POSSIBLE ARGUMENT OF A DOUBLE PRECISION
    ! GONIOMETRIC FUNCTION (cos, sin, ...)
    ! THE REASON WHY THESE PARAMETERS ARE NEEDED AS THEY ARE DEFINED WILL
    ! BE EXPLAINED IN THE CODE BY MEANS OF COMMENTS
    ! 
    ! 
    ! PARAMETER LIST
    ! XI     = REAL      PART OF Z
    ! YI     = IMAGINARY PART OF Z
    ! U      = REAL      PART OF W(Z)
    ! V      = IMAGINARY PART OF W(Z)
    ! FLAG   = AN ERROR FLAG INDICATING WHETHER OVERFLOW WILL
    ! OCCUR OR NOT; TYPE LOGICAL;
    ! THE VALUES OF THIS VARIABLE HAVE THE FOLLOWING
    ! MEANING :
    ! FLAG=.FALSE. : NO ERROR CONDITION
    ! FLAG=.TRUE.  : OVERFLOW WILL OCCUR, THE ROUTINE
    ! BECOMES INACTIVE
    ! XI, YI      ARE THE INPUT-PARAMETERS
    ! U, V, FLAG  ARE THE OUTPUT-PARAMETERS
    ! 
    ! FURTHERMORE THE PARAMETER FACTOR EQUALS 2/SQRT(PI)
    ! 
    ! THE ROUTINE IS NOT UNDERFLOW-PROTECTED BUT ANY VARIABLE CAN BE
    ! PUT TO 0 UPON UNDERFLOW;
    ! 
    ! REFERENCE - GPM POPPE, CMJ WIJERS; MORE EFFICIENT COMPUTATION OF
    ! THE COMPLEX ERROR-FUNCTION, ACM TRANS. MATH. SOFTWARE.
    ! 
    ! 
    ! 
    ! 
    ! 
    ! 
    REAL(real_8)                             :: xi, yi, u, v
    LOGICAL                                  :: flag

    REAL(real_8), PARAMETER :: factor = 1.12837916709551257388_real_8, &
      rmaxexp = 708.503061461606_real_8, &
      rmaxgoni = 3.53711887601422e+15_real_8 , rmaxreal = 0.5e+154_real_8

    INTEGER                                  :: i, j, kapn, n, np1, nu
    LOGICAL                                  :: a, b
    REAL(real_8) :: c, daux, h, h2, qlambda, qrho, rx, ry, sx, sy, tx, ty, &
      u1, u2, v1, v2, w1, x, xabs, xabsq, xaux, xquad, xsum, y, yabs, yquad, &
      ysum

! 

    flag = .FALSE.
    ! 
    xabs = ABS(xi)
    yabs = ABS(yi)
    x    = xabs/6.3_real_8
    y    = yabs/4.4_real_8
    ! 
    ! 
    ! THE FOLLOWING IF-STATEMENT PROTECTS
    ! QRHO = (X**2 + Y**2) AGAINST OVERFLOW
    ! 
    IF ((xabs.GT.rmaxreal).OR.(yabs.GT.rmaxreal)) GOTO 100
    ! 
    qrho = x**2 + y**2
    ! 
    xabsq = xabs**2
    xquad = xabsq - yabs**2
    yquad = 2*xabs*yabs
    ! 
    a     = qrho.LT.0.085264_real_8
    ! 
    IF (a) THEN
       ! 
       ! IF (QRHO.LT.0.085264_real_8) THEN THE FADDEEVA-FUNCTION IS EVALUATED
       ! USING A POWER-SERIES (ABRAMOWITZ/STEGUN, EQUATION (7.1.5), P.297)
       ! N IS THE MINIMUM NUMBER OF TERMS NEEDED TO OBTAIN THE REQUIRED
       ! ACCURACY
       ! 
       qrho  = (1._real_8-0.85_real_8*y)*SQRT(qrho)
       n     = idnint(6._real_8 + 72._real_8*qrho)
       j     = 2*n+1
       xsum  = 1.0_real_8/REAL(j,kind=real_8)
       ysum  = 0.0_real_8
       DO i=n, 1, -1
          j    = j - 2
          xaux = (xsum*xquad - ysum*yquad)/REAL(i,kind=real_8)
          ysum = (xsum*yquad + ysum*xquad)/REAL(i,kind=real_8)
          xsum = xaux + 1._real_8/REAL(j,kind=real_8)
       END DO
       u1   = -factor*(xsum*yabs + ysum*xabs) + 1.0_real_8
       v1   =  factor*(xsum*xabs - ysum*yabs)
       daux =  EXP(-xquad)
       u2   =  daux*COS(yquad)
       v2   = -daux*SIN(yquad)
       ! 
       u    = u1*u2 - v1*v2
       v    = u1*v2 + v1*u2
       ! 
    ELSE
       ! 
       ! IF (QRHO.GT.1.O) THEN W(Z) IS EVALUATED USING THE LAPLACE
       ! CONTINUED FRACTION
       ! NU IS THE MINIMUM NUMBER OF TERMS NEEDED TO OBTAIN THE REQUIRED
       ! ACCURACY
       ! 
       ! IF ((QRHO.GT.0.085264_real_8).AND.(QRHO.LT.1.0)) THEN W(Z) IS EVALUATED
       ! BY A TRUNCATED TAYLOR EXPANSION, WHERE THE LAPLACE CONTINUED FRACTION
       ! IS USED TO CALCULATE THE DERIVATIVES OF W(Z)
       ! KAPN IS THE MINIMUM NUMBER OF TERMS IN THE TAYLOR EXPANSION NEEDED
       ! TO OBTAIN THE REQUIRED ACCURACY
       ! NU IS THE MINIMUM NUMBER OF TERMS OF THE CONTINUED FRACTION NEEDED
       ! TO CALCULATE THE DERIVATIVES WITH THE REQUIRED ACCURACY
       ! 
       ! 
       IF (qrho.GT.1._real_8) THEN
          h    = 0.0_real_8
          kapn = 0
          qrho = SQRT(qrho)
          nu   = idint(3._real_8 + (1442._real_8/(26._real_8*qrho+77._real_8)))
       ELSE
          qrho = (1._real_8-y)*SQRT(1._real_8-qrho)
          h    = 1.88_real_8*qrho
          h2   = 2.0_real_8*h
          kapn = idnint(7._real_8  + 34._real_8*qrho)
          nu   = idnint(16._real_8 + 26._real_8*qrho)
       ENDIF
       ! 
       b = (h.GT.0.0_real_8)
       ! 
       IF (b) qlambda = h2**kapn
       ! 
       rx = 0.0_real_8
       ry = 0.0_real_8
       sx = 0.0_real_8
       sy = 0.0_real_8
       ! 
       DO n=nu, 0, -1
          np1 = n + 1
          tx  = yabs + h + np1*rx
          ty  = xabs - np1*ry
          c   = 0.5_real_8/(tx**2 + ty**2)
          rx  = c*tx
          ry  = c*ty
          IF ((b).AND.(n.LE.kapn)) THEN
             tx = qlambda + sx
             sx = rx*tx - ry*sy
             sy = ry*tx + rx*sy
             qlambda = qlambda/h2
          ENDIF
       END DO
       ! 
       IF (h.EQ.0._real_8) THEN
          u = factor*rx
          v = factor*ry
       ELSE
          u = factor*sx
          v = factor*sy
       ENDIF
       ! 
       IF (yabs.EQ.0._real_8) u = EXP(-xabs**2)
       ! 
    ENDIF
    ! 
    ! 
    ! 
    ! EVALUATION OF W(Z) IN THE OTHER QUADRANTS
    ! 
    ! 
    IF (yi.LT.0._real_8) THEN
       IF (a) THEN
          u2    = 2._real_8*u2
          v2    = 2._real_8*v2
       ELSE
          xquad =  -xquad
          ! 
          ! 
          ! THE FOLLOWING IF-STATEMENT PROTECTS 2*EXP(-Z**2)
          ! AGAINST OVERFLOW
          ! 
          IF ((yquad.GT.rmaxgoni).OR.&
               (xquad.GT.rmaxexp)) GOTO 100
          ! 
          w1 =  2._real_8*EXP(xquad)
          u2  =  w1*COS(yquad)
          v2  = -w1*SIN(yquad)
       ENDIF
       ! 
       u = u2 - u
       v = v2 - v
       IF (xi.GT.0._real_8) v = -v
    ELSE
       IF (xi.LT.0._real_8) v = -v
    ENDIF
    ! 
    RETURN
    ! 
100 flag = .TRUE.
    RETURN
    ! 
  END SUBROUTINE wofz
  ! ==================================================================
  FUNCTION bessj0(x)
    REAL(real_8)                             :: x, bessj0

    REAL(real_8), PARAMETER :: p1 = 1._real_8, p2 = -0.1098628627e-2_real_8, &
      p3 = 0.2734510407e-4_real_8, p4 = -0.2073370639e-5_real_8, &
      p5 = 0.2093887211e-6_real_8, q1 = -0.1562499995e-1_real_8, &
      q2 = 0.1430488765e-3_real_8, q3 = -0.6911147651e-5_real_8, &
      q4 = 0.7621095161e-6_real_8, q5 = -0.934945152e-7_real_8, &
      r1 = 57568490574._real_8, r2 = -13362590354._real_8, &
      r3 = 651619640.7_real_8, r4 = -11214424.18_real_8, &
      r5 = 77392.33017_real_8, r6 = -184.9052456_real_8, &
      s1 = 57568490411._real_8, s2 = 1029532985._real_8, &
      s3 = 9494680.718_real_8, s4 = 59272.64853_real_8, &
      s5 = 267.8532712_real_8
    REAL(real_8), PARAMETER :: s6 = 1._real_8

    REAL(real_8)                             :: ax, xx, y, z

    IF (ABS(x).LT.8._real_8)THEN
       y=x**2
       bessj0=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+y*&
            (s4+y*(s5+y*s6)))))
    ELSE
       ax=ABS(x)
       z=8._real_8/ax
       y=z**2
       xx=ax-0.785398164_real_8
       bessj0=SQRT(0.636619772_real_8/ax)*(COS(xx)*(p1+y*(p2+y*(p3+y*(p4+y*&
            p5))))-z*SIN(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
    ENDIF
    RETURN
  END FUNCTION bessj0
  ! ==================================================================
  FUNCTION bessj1(x)
    REAL(real_8)                             :: x, bessj1

    REAL(real_8), PARAMETER :: p1 = 1._real_8, p2 = 0.183105e-2_real_8, &
      p3 = -0.3516396496e-4_real_8, p4 = 0.2457520174e-5_real_8, &
      p5 = -0.240337019e-6_real_8, q1 = 0.04687499995_real_8, &
      q2 = -0.2002690873e-3_real_8, q3 = 0.8449199096e-5_real_8, &
      q4 = -0.88228987e-6_real_8, q5 = 0.105787412e-6_real_8, &
      r1 = 72362614232._real_8, r2 = -7895059235._real_8, &
      r3 = 242396853.1_real_8, r4 = -2972611.439_real_8, &
      r5 = 15704.48260_real_8, r6 = -30.16036606_real_8, &
      s1 = 144725228442._real_8, s2 = 2300535178._real_8, &
      s3 = 18583304.74_real_8, s4 = 99447.43394_real_8, &
      s5 = 376.9991397_real_8, s6 = 1._real_8

    REAL(real_8)                             :: ax, xx, y, z

    IF (ABS(x).LT.8._real_8)THEN
       y=x**2
       bessj1=x*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/&
            (s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))
    ELSE
       ax=ABS(x)
       z=8._real_8/ax
       y=z**2
       xx=ax-2.356194491_real_8
       bessj1=SQRT(0.636619772_real_8/ax)*(COS(xx)*(p1+y*(p2+y*(p3+y*(p4+&
            y*p5))))-z*SIN(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))*dsign(&
            1._real_8,x)
    ENDIF
    RETURN
  END FUNCTION bessj1
  ! ==================================================================
  FUNCTION bessk0(x)
    REAL(real_8)                             :: x, bessk0

    REAL(real_8), PARAMETER :: p1 = -0.57721566_real_8, &
      p2 = 0.42278420_real_8, p3 = 0.23069756_real_8, &
      p4 = 0.3488590e-1_real_8, p5 = 0.262698e-2_real_8, &
      p6 = 0.10750e-3_real_8, p7 = 0.74e-5_real_8, q1 = 1.25331414_real_8, &
      q2 = -0.7832358e-1_real_8, q3 = 0.2189568e-1_real_8, &
      q4 = -0.1062446e-1_real_8, q5 = 0.587872e-2_real_8, &
      q6 = -0.251540e-2_real_8, q7 = 0.53208e-3_real_8

    REAL(real_8)                             :: y

    IF (x.LE.2.0_real_8) THEN
       y=x*x/4.0_real_8
       bessk0=(-LOG(x/2._real_8)*bessi0(x))+(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(&
            p6+y*p7))))))
    ELSE
       y=(2.0_real_8/x)
       bessk0=(EXP(-x)/SQRT(x))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*&
            q7))))))
    ENDIF
    RETURN
  END FUNCTION bessk0
  ! ==================================================================
  FUNCTION bessk1(x)
    REAL(real_8)                             :: x, bessk1

    REAL(real_8), PARAMETER :: p1 = 1.0_real_8, p2 = 0.15443144_real_8, &
      p3 = -0.67278579_real_8, p4 = -0.18156897_real_8, &
      p5 = -0.1919402e-1_real_8, p6 = -0.110404e-2_real_8, &
      p7 = -0.4686e-4_real_8, q1 = 1.25331414_real_8, q2 = 0.23498619_real_8, &
      q3 = -0.3655620e-1_real_8, q4 = 0.1504268e-1_real_8, &
      q5 = -0.780353e-2_real_8, q6 = 0.325614e-2_real_8, &
      q7 = -0.68245e-3_real_8

    REAL(real_8)                             :: y

    IF (x.LE.2._real_8) THEN
       y=x*x/4._real_8
       bessk1=(LOG(x/2._real_8)*bessi1(x))+(1._real_8/x)*(p1+y*(p2+y*(p3+y*(p4+&
            y*(p5+y*(p6+y*p7))))))
    ELSE
       y=2._real_8/x
       bessk1=(EXP(-x)/SQRT(x))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*&
            q7))))))
    ENDIF
    RETURN
  END FUNCTION bessk1
  ! ==================================================================
  FUNCTION bessi0(x)
    REAL(real_8)                             :: x, bessi0

    REAL(real_8), PARAMETER :: p1 = 1.0_real_8, p2 = 3.5156229_real_8, &
      p3 = 3.0899424_real_8, p4 = 1.2067492_real_8, p5 = 0.2659732_real_8, &
      p6 = 0.360768e-1_real_8, p7 = 0.45813e-2_real_8, &
      q1 = 0.39894228_real_8, q2 = 0.1328592e-1_real_8, &
      q3 = 0.225319e-2_real_8, q4 = -0.157565e-2_real_8, &
      q5 = 0.916281e-2_real_8, q6 = -0.2057706e-1_real_8, &
      q7 = 0.2635537e-1_real_8, q8 = -0.1647633e-1_real_8, &
      q9 = 0.392377e-2_real_8

    REAL(real_8)                             :: ax, y

    IF (ABS(x).LT.3.75_real_8) THEN
       y=(x/3.75_real_8)**2
       bessi0=p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))))
    ELSE
       ax=ABS(x)
       y=3.75_real_8/ax
       bessi0=(EXP(ax)/SQRT(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*&
            (q7+y*(q8+y*q9))))))))
    ENDIF
    RETURN
  END FUNCTION bessi0
  ! ==================================================================
  FUNCTION bessi1(x)
    REAL(real_8)                             :: x, bessi1

    REAL(real_8), PARAMETER :: p1 = 0.5_real_8, p2 = 0.87890594_real_8, &
      p3 = 0.51498869_real_8, p4 = 0.15084934_real_8, &
      p5 = 0.2658733e-1_real_8, p6 = 0.301532e-2_real_8, &
      p7 = 0.32411e-3_real_8, q1 = 0.39894228_real_8, &
      q2 = -0.3988024e-1_real_8, q3 = -0.362018e-2_real_8, &
      q4 = 0.163801e-2_real_8, q5 = -0.1031555e-1_real_8, &
      q6 = 0.2282967e-1_real_8, q7 = -0.2895312e-1_real_8, &
      q8 = 0.1787654e-1_real_8, q9 = -0.420059e-2_real_8

    REAL(real_8)                             :: ax, y

    IF (ABS(x).LT.3.75_real_8) THEN
       y=(x/3.75_real_8)**2
       bessi1=x*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
    ELSE
       ax=ABS(x)
       y=3.75_real_8/ax
       bessi1=(EXP(ax)/SQRT(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*&
            (q7+y*(q8+y*q9))))))))
       IF (x.LT.0._real_8)bessi1=-bessi1
    ENDIF
    RETURN
  END FUNCTION bessi1

END MODULE mtin_utils
