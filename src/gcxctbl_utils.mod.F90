#include "cpmd_global.h"

MODULE gcxctbl_utils
  USE cnst,                            ONLY: fpi,&
                                             pi
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE system,                          ONLY: cntr,&
                                             fpar,&
                                             parm

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: gcxonly
  PUBLIC :: gcxlyp
  PUBLIC :: gcxp86
  PUBLIC :: gcgga
  PUBLIC :: gcpbe
  PUBLIC :: gcrevpbe
  PUBLIC :: gcsxonly
  PUBLIC :: gcsxp86
  PUBLIC :: gcsxlyp
  PUBLIC :: gcspbe
  PUBLIC :: gcsrevpbe

CONTAINS

#ifndef __VECTOR 
  ! ==================================================================
  SUBROUTINE gcxonly(b1,sgcx,sgcc,rhoe,v,vtmp,grad,flops)
    ! ==--------------------------------------------------------------==
    REAL(real_8) :: b1, sgcx, sgcc, rhoe(fpar%kr1,fpar%kr2s,fpar%kr3s)
    COMPLEX(KIND=real_8), DIMENSION(fpar%kr1&
      , fpar%kr2s, fpar%kr3s)                :: v
    REAL(real_8) :: vtmp(fpar%kr1,fpar%kr2s,fpar%kr3s), &
      grad(fpar%kr1,fpar%kr2s,fpar%kr3s,*), flops

    REAL(real_8), PARAMETER                  :: eps2 = 1.e-20_real_8, &
                                                small = 1.e-24_real_8

    INTEGER                                  :: i1, i2, i3, ic, id, ii2, ii3, &
                                                ncount
    REAL(real_8) :: d1, d2, d3, d4, dd1, dd2, dd3, dd4, drho, dsmoo1, dsmoo2, &
      dsmoo3, dsmoo4, ee1, ee2, ee3, ee4, eexc_1, eexc_2, eexc_3, eexc_4, &
      fadd, fdiv, fmult, fspez, grho1, grho2, grho3, grho4, ogceps, rho431, &
      rho432, rho433, rho434, rhou1, rhou2, rhou3, rhou4, sa2b8_1, sa2b8_2, &
      sa2b8_3, sa2b8_4, shm1_1, shm1_2, shm1_3, shm1_4, smoo1, smoo2, smoo3, &
      smoo4, two13, vexc11, vexc12, vexc13, vexc14, vexc21, vexc22, vexc23, &
      vexc24, x1, x2, x3, x4, xx1, xx2, xx3, xx4, yy1, yy2, yy3, yy4

! ==--------------------------------------------------------------==

    two13=2._real_8**0.333333333333333_real_8
    sgcx = 0.0_real_8
    sgcc = 0.0_real_8
    ogceps=1._real_8/cntr%gceps
    x1=1._real_8
    x2=1._real_8
    x3=1._real_8
    x4=1._real_8
    ic=0
    id=0
    ii3=MOD(parm%nr3,2)
    ii2=MOD(parm%nr2,2)
    IF (ii3.NE.0.OR.ii2.NE.0) CALL stopgm('GCXONLY','ODD DIMENSION',& 
         __LINE__,__FILE__)
    !$omp parallel do private(I1,I2,I3,RHOU1,RHOU2,RHOU3,RHOU4,D1,D2,D3,D4) &
    !$omp  firstprivate(X1,X2,X3,X4) &
    !$omp  private(NCOUNT,RHO431,RHO432,RHO433,RHO434,SMOO1,DSMOO1,DRHO) &
    !$omp  private(SMOO2,DSMOO2,SMOO3,DSMOO3,SMOO4,DSMOO4,GRHO1,GRHO2) &
    !$omp  private(GRHO3,GRHO4,XX1,XX2,XX3,XX4,YY1,YY2,YY3,YY4) &
    !$omp  private(SA2B8_1,SA2B8_2,SA2B8_3,SA2B8_4) &
    !$omp  private(SHM1_1,SHM1_2,SHM1_3,SHM1_4) &
    !$omp  private(DD1,DD2,DD3,DD4,EEXC_1,EEXC_2,EEXC_3,EEXC_4) &
    !$omp  private(EE1,EE2,EE3,EE4,VEXC11,VEXC12,VEXC13,VEXC14,VEXC21) &
    !$omp  private(VEXC22,VEXC23,VEXC24) &
    !$omp  reduction(+:SGCX,ID,IC) __COLLAPSE3
    DO i3=1,parm%nr3,2
       DO i2=1,parm%nr2,2
          DO i1=1,parm%nr1
             rhou1=MAX(rhoe(i1,i2+0,i3+0),small)
             rhou2=MAX(rhoe(i1,i2+1,i3+0),small)
             rhou3=MAX(rhoe(i1,i2+0,i3+1),small)
             rhou4=MAX(rhoe(i1,i2+1,i3+1),small)
             IF ((rhou1+rhou2+rhou3+rhou4).GT.5._real_8*small) THEN
                id=id+1
                DO ncount=1,100
                   d1=x1-rhou1/(x1*x1)
                   d2=x2-rhou2/(x2*x2)
                   d3=x3-rhou3/(x3*x3)
                   d4=x4-rhou4/(x4*x4)
                   ic=ic+4
                   x1=x1-0.333333333333333_real_8*d1
                   x2=x2-0.333333333333333_real_8*d2
                   x3=x3-0.333333333333333_real_8*d3
                   x4=x4-0.333333333333333_real_8*d4
                   IF (d1**2+d2**2+d3**2+d4**2.LT.eps2) GOTO 10
                ENDDO
                x1=rhou1**(0.33333333333333_real_8)
                x2=rhou2**(0.33333333333333_real_8)
                x3=rhou3**(0.33333333333333_real_8)
                x4=rhou4**(0.33333333333333_real_8)
10              CONTINUE
                rho431=x1*rhou1
                rho432=x2*rhou2
                rho433=x3*rhou3
                rho434=x4*rhou4
                IF (rhou1.GT.2.25_real_8*cntr%gceps) THEN
                   smoo1=1._real_8
                   dsmoo1=0._real_8
                ELSEIF (rhou1.LT.0.25_real_8*cntr%gceps) THEN
                   smoo1=0._real_8
                   dsmoo1=0._real_8
                ELSE
                   drho=(rhou1-0.25e9_real_8*cntr%gceps)*ogceps
                   smoo1=drho*drho*(0.75_real_8-0.25_real_8*drho)
                   dsmoo1=1.5_real_8*drho*ogceps*(1._real_8-0.5_real_8*drho)
                ENDIF
                IF (rhou2.GT.2.25_real_8*cntr%gceps) THEN
                   smoo2=1._real_8
                   dsmoo2=0._real_8
                ELSEIF (rhou2.LT.0.25_real_8*cntr%gceps) THEN
                   smoo2=0._real_8
                   dsmoo2=0._real_8
                ELSE
                   drho=(rhou2-0.25_real_8*cntr%gceps)*ogceps
                   smoo2=drho*drho*(0.75_real_8-0.25_real_8*drho)
                   dsmoo2=1.5_real_8*drho*ogceps*(1._real_8-0.5_real_8*drho)
                ENDIF
                IF (rhou3.GT.2.25_real_8*cntr%gceps) THEN
                   smoo3=1._real_8
                   dsmoo3=0._real_8
                ELSEIF (rhou3.LT.0.25_real_8*cntr%gceps) THEN
                   smoo3=0._real_8
                   dsmoo3=0._real_8
                ELSE
                   drho=(rhou3-0.25_real_8*cntr%gceps)*ogceps
                   smoo3=drho*drho*(0.75_real_8-0.25_real_8*drho)
                   dsmoo3=1.5_real_8*drho*ogceps*(1._real_8-0.5_real_8*drho)
                ENDIF
                IF (rhou4.GT.2.25_real_8*cntr%gceps) THEN
                   smoo4=1._real_8
                   dsmoo4=0._real_8
                ELSEIF (rhou4.LT.0.25_real_8*cntr%gceps) THEN
                   smoo4=0._real_8
                   dsmoo4=0._real_8
                ELSE
                   drho=(rhou4-0.25_real_8*cntr%gceps)*ogceps
                   smoo4=drho*drho*(0.75_real_8-0.25_real_8*drho)
                   dsmoo4=1.5_real_8*drho*ogceps*(1._real_8-0.5_real_8*drho)
                ENDIF
                grho1=SQRT(grad(i1,i2+0,i3+0,1))
                grho2=SQRT(grad(i1,i2+1,i3+0,1))
                grho3=SQRT(grad(i1,i2+0,i3+1,1))
                grho4=SQRT(grad(i1,i2+1,i3+1,1))
                xx1=two13*grho1/rho431
                xx2=two13*grho2/rho432
                xx3=two13*grho3/rho433
                xx4=two13*grho4/rho434
                yy1=xx1*xx1
                yy2=xx2*xx2
                yy3=xx3*xx3
                yy4=xx4*xx4
                sa2b8_1=SQRT(1._real_8+yy1)
                sa2b8_2=SQRT(1._real_8+yy2)
                sa2b8_3=SQRT(1._real_8+yy3)
                sa2b8_4=SQRT(1._real_8+yy4)
                shm1_1=LOG(xx1+sa2b8_1)
                shm1_2=LOG(xx2+sa2b8_2)
                shm1_3=LOG(xx3+sa2b8_3)
                shm1_4=LOG(xx4+sa2b8_4)
                dd1=1._real_8+6._real_8*b1*xx1*shm1_1
                dd2=1._real_8+6._real_8*b1*xx2*shm1_2
                dd3=1._real_8+6._real_8*b1*xx3*shm1_3
                dd4=1._real_8+6._real_8*b1*xx4*shm1_4

                eexc_1=-b1*yy1/dd1/two13
                eexc_2=-b1*yy2/dd2/two13
                eexc_3=-b1*yy3/dd3/two13
                eexc_4=-b1*yy4/dd4/two13

                ee1=6._real_8*b1*yy1/sa2b8_1 - 1._real_8
                ee2=6._real_8*b1*yy2/sa2b8_2 - 1._real_8
                ee3=6._real_8*b1*yy3/sa2b8_3 - 1._real_8
                ee4=6._real_8*b1*yy4/sa2b8_4 - 1._real_8
                vexc11=-1.333333333333_real_8/two13*b1*yy1*ee1/(dd1*dd1)
                vexc12=-1.333333333333_real_8/two13*b1*yy2*ee2/(dd2*dd2)
                vexc13=-1.333333333333_real_8/two13*b1*yy3*ee3/(dd3*dd3)
                vexc14=-1.333333333333_real_8/two13*b1*yy4*ee4/(dd4*dd4)
                vexc21=two13*b1*(ee1-dd1)/(dd1*dd1)
                vexc22=two13*b1*(ee2-dd2)/(dd2*dd2)
                vexc23=two13*b1*(ee3-dd3)/(dd3*dd3)
                vexc24=two13*b1*(ee4-dd4)/(dd4*dd4)

                sgcx=sgcx+(smoo1*eexc_1*rho431+smoo2*eexc_2*rho432+&
                     smoo3*eexc_3*rho433+smoo4*eexc_4*rho434)
                v(i1,i2+0,i3+0)=dsmoo1*eexc_1*rho431+smoo1*vexc11*x1
                v(i1,i2+1,i3+0)=dsmoo2*eexc_2*rho432+smoo2*vexc12*x2
                v(i1,i2+0,i3+1)=dsmoo3*eexc_3*rho433+smoo3*vexc13*x3
                v(i1,i2+1,i3+1)=dsmoo4*eexc_4*rho434+smoo4*vexc14*x4
                vtmp(i1,i2+0,i3+0)=smoo1*vexc21/rho431
                vtmp(i1,i2+1,i3+0)=smoo2*vexc22/rho432
                vtmp(i1,i2+0,i3+1)=smoo3*vexc23/rho433
                vtmp(i1,i2+1,i3+1)=smoo4*vexc24/rho434
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    ! ..counting floating point operations
    fadd   = ic * 11._real_8 + 4._real_8*id *  7._real_8
    fmult  = ic * 12._real_8 + 4._real_8*id * 24._real_8
    fdiv   = ic *  4._real_8 + 4._real_8*id *  7._real_8
    fspez  =              4._real_8*id *  3._real_8
    flops=fadd+fmult+fdiv+fspez
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gcxonly
  ! ==================================================================
  SUBROUTINE gcxlyp(b1,sgcx,sgcc,rhoe,v,vtmp,grad,flops)
    ! ==--------------------------------------------------------------==
    REAL(real_8) :: b1, sgcx, sgcc, rhoe(fpar%kr1,fpar%kr2s,fpar%kr3s)
    COMPLEX(real_8) :: v(fpar%kr1,fpar%kr2s,fpar%kr3s)
    REAL(real_8) :: vtmp(fpar%kr1,fpar%kr2s,fpar%kr3s), &
      grad(fpar%kr1,fpar%kr2s,fpar%kr3s,*), flops

    REAL(real_8), PARAMETER :: a = 0.04918_real_8, b = 0.132_real_8, &
      c = 0.2533_real_8, d = 0.349_real_8, eps2 = 1.e-20_real_8, &
      small = 1.e-24_real_8

    INTEGER                                  :: i1, i2, i3, ic, id, ii2, ii3, &
                                                ncount
    REAL(real_8) :: d1, d2, d3, d4, dd1, dd2, dd3, dd4, dom1, dom2, dom3, &
      dom4, drho, dsmoo1, dsmoo2, dsmoo3, dsmoo4, dxl1, dxl2, dxl3, dxl4, &
      ee1, ee2, ee3, ee4, eecc_1, eecc_2, eecc_3, eecc_4, eexc_1, eexc_2, &
      eexc_3, eexc_4, fadd, fdiv, ff1, ff2, ff3, ff4, fmult, fspez, grho1, &
      grho2, grho3, grho4, grhos1, grhos2, grhos3, grhos4, ogceps, om1, om2, &
      om3, om4, r1, r2, r21, r22, r23, r24, r3, r4, r41, r42, r43, r44, r51, &
      r52, r53, r54, rho431, rho432, rho433, rho434, rhou1, rhou2, rhou3, &
      rhou4, sa2b8_1, sa2b8_2, sa2b8_3, sa2b8_4, shm1_1, shm1_2, shm1_3, &
      shm1_4, smoo1, smoo2, smoo3, smoo4, two13
    REAL(real_8) :: v1cc_1, v1cc_2, v1cc_3, v1cc_4, v2cc_1, v2cc_2, v2cc_3, &
      v2cc_4, vexc11, vexc12, vexc13, vexc14, vexc21, vexc22, vexc23, vexc24, &
      x1, x2, x3, x4, xl1, xl2, xl3, xl4, xx1, xx2, xx3, xx4, yy1, yy2, yy3, &
      yy4

! ==--------------------------------------------------------------==

    two13=2._real_8**0.333333333333333_real_8
    sgcx = 0.0_real_8
    sgcc = 0.0_real_8
    ogceps=1._real_8/cntr%gceps
    x1=1._real_8
    x2=1._real_8
    x3=1._real_8
    x4=1._real_8
    ic=0
    id=0
    ii3=MOD(parm%nr3,2)
    ii2=MOD(parm%nr2,2)
    IF (ii3.NE.0.OR.ii2.NE.0) CALL stopgm('GCXLYP','ODD DIMENSION',& 
         __LINE__,__FILE__)
    !$omp parallel do private(I1,I2,I3,RHOU1,RHOU2,RHOU3,RHOU4,D1,D2,D3,D4) &
    !$omp  firstprivate(X1,X2,X3,X4) &
    !$omp  private(NCOUNT,RHO431,RHO432,RHO433,RHO434,SMOO1,DSMOO1,DRHO) &
    !$omp  private(SMOO2,DSMOO2,SMOO3,DSMOO3,SMOO4,DSMOO4,GRHO1,GRHO2) &
    !$omp  private(GRHO3,GRHO4,GRHOS1,GRHOS2,GRHOS3,GRHOS4,XX1,XX2) &
    !$omp  private(XX3,XX4,YY1,YY2,YY3,YY4,SA2B8_1,SA2B8_2,SA2B8_3) &
    !$omp  private(SA2B8_4,SHM1_1,SHM1_2,SHM1_3,SHM1_4,DD1,DD2,DD3) &
    !$omp  private(DD4,EEXC_1,EEXC_2,EEXC_3,EEXC_4,EE1,EE2,EE3,EE4) &
    !$omp  private(VEXC11,VEXC12,VEXC13,VEXC14,VEXC21,VEXC22,VEXC23) &
    !$omp  private(VEXC24,R1,R2,R3,R4,R21,R22,R23,R24,R41,R42,R43) &
    !$omp  private(R44,R51,R52,R53,R54,OM1,OM2,OM3,OM4,XL1,XL2,XL3) &
    !$omp  private(XL4,FF1,FF2,FF3,FF4,EECC_1,EECC_2,EECC_3,EECC_4) &
    !$omp  private(DOM1,DOM2,DOM3,DOM4,DXL1,DXL2,DXL3,DXL4,V1CC_1) &
    !$omp  private(V1CC_2,V1CC_3,V1CC_4,V2CC_1,V2CC_2,V2CC_3,V2CC_4) &
    !$omp  reduction(+:SGCX,SGCC,IC,ID) __COLLAPSE3
    DO i3=1,parm%nr3,2
       DO i2=1,parm%nr2,2
          DO i1=1,parm%nr1
             rhou1=MAX(rhoe(i1,i2+0,i3+0),small)
             rhou2=MAX(rhoe(i1,i2+1,i3+0),small)
             rhou3=MAX(rhoe(i1,i2+0,i3+1),small)
             rhou4=MAX(rhoe(i1,i2+1,i3+1),small)
             IF ((rhou1+rhou2+rhou3+rhou4).GT.5._real_8*small) THEN
                id=id+1
                DO ncount=1,100
                   d1=x1-rhou1/x1**2
                   d2=x2-rhou2/x2**2
                   d3=x3-rhou3/x3**2
                   d4=x4-rhou4/x4**2
                   ic=ic+4
                   x1=x1-0.333333333333333_real_8*d1
                   x2=x2-0.333333333333333_real_8*d2
                   x3=x3-0.333333333333333_real_8*d3
                   x4=x4-0.333333333333333_real_8*d4
                   IF (d1**2+d2**2+d3**2+d4**2.LT.eps2) GOTO 10
                ENDDO
                x1=rhou1**(1._real_8/3._real_8)
                x2=rhou2**(1._real_8/3._real_8)
                x3=rhou3**(1._real_8/3._real_8)
                x4=rhou4**(1._real_8/3._real_8)
10              CONTINUE
                rho431=x1*rhou1
                rho432=x2*rhou2
                rho433=x3*rhou3
                rho434=x4*rhou4
                IF (rhou1.GT.2.25_real_8*cntr%gceps) THEN
                   smoo1=1._real_8
                   dsmoo1=0._real_8
                ELSEIF (rhou1.LT.0.25_real_8*cntr%gceps) THEN
                   smoo1=0._real_8
                   dsmoo1=0._real_8
                ELSE
                   drho=(rhou1-0.25_real_8*cntr%gceps)*ogceps
                   smoo1=drho*drho*(0.75_real_8-0.25_real_8*drho)
                   dsmoo1=1.5_real_8*drho*ogceps*(1._real_8-0.5_real_8*drho)
                ENDIF
                IF (rhou2.GT.2.25_real_8*cntr%gceps) THEN
                   smoo2=1._real_8
                   dsmoo2=0._real_8
                ELSEIF (rhou2.LT.0.25_real_8*cntr%gceps) THEN
                   smoo2=0._real_8
                   dsmoo2=0._real_8
                ELSE
                   drho=(rhou2-0.25_real_8*cntr%gceps)*ogceps
                   smoo2=drho*drho*(0.75_real_8-0.25_real_8*drho)
                   dsmoo2=1.5_real_8*drho*ogceps*(1._real_8-0.5_real_8*drho)
                ENDIF
                IF (rhou3.GT.2.25_real_8*cntr%gceps) THEN
                   smoo3=1._real_8
                   dsmoo3=0._real_8
                ELSEIF (rhou3.LT.0.25_real_8*cntr%gceps) THEN
                   smoo3=0._real_8
                   dsmoo3=0._real_8
                ELSE
                   drho=(rhou3-0.25_real_8*cntr%gceps)*ogceps
                   smoo3=drho*drho*(0.75_real_8-0.25_real_8*drho)
                   dsmoo3=1.5_real_8*drho*ogceps*(1._real_8-0.5_real_8*drho)
                ENDIF
                IF (rhou4.GT.2.25_real_8*cntr%gceps) THEN
                   smoo4=1._real_8
                   dsmoo4=0._real_8
                ELSEIF (rhou4.LT.0.25_real_8*cntr%gceps) THEN
                   smoo4=0._real_8
                   dsmoo4=0._real_8
                ELSE
                   drho=(rhou4-0.25_real_8*cntr%gceps)*ogceps
                   smoo4=drho*drho*(0.75_real_8-0.25_real_8*drho)
                   dsmoo4=1.5_real_8*drho*ogceps*(1._real_8-0.5_real_8*drho)
                ENDIF
                grho1=grad(i1,i2+0,i3+0,1)
                grho2=grad(i1,i2+1,i3+0,1)
                grho3=grad(i1,i2+0,i3+1,1)
                grho4=grad(i1,i2+1,i3+1,1)
                grhos1=SQRT(grho1)
                grhos2=SQRT(grho2)
                grhos3=SQRT(grho3)
                grhos4=SQRT(grho4)
                ! ..exchange
                xx1=two13*grhos1/rho431
                xx2=two13*grhos2/rho432
                xx3=two13*grhos3/rho433
                xx4=two13*grhos4/rho434
                yy1=xx1*xx1
                yy2=xx2*xx2
                yy3=xx3*xx3
                yy4=xx4*xx4
                sa2b8_1=SQRT(1._real_8+yy1)
                sa2b8_2=SQRT(1._real_8+yy2)
                sa2b8_3=SQRT(1._real_8+yy3)
                sa2b8_4=SQRT(1._real_8+yy4)
                shm1_1=LOG(xx1+sa2b8_1)
                shm1_2=LOG(xx2+sa2b8_2)
                shm1_3=LOG(xx3+sa2b8_3)
                shm1_4=LOG(xx4+sa2b8_4)
                dd1=1._real_8+6._real_8*b1*xx1*shm1_1
                dd2=1._real_8+6._real_8*b1*xx2*shm1_2
                dd3=1._real_8+6._real_8*b1*xx3*shm1_3
                dd4=1._real_8+6._real_8*b1*xx4*shm1_4

                eexc_1=-b1*yy1/dd1/two13
                eexc_2=-b1*yy2/dd2/two13
                eexc_3=-b1*yy3/dd3/two13
                eexc_4=-b1*yy4/dd4/two13

                ee1=6._real_8*b1*yy1/sa2b8_1 - 1._real_8
                ee2=6._real_8*b1*yy2/sa2b8_2 - 1._real_8
                ee3=6._real_8*b1*yy3/sa2b8_3 - 1._real_8
                ee4=6._real_8*b1*yy4/sa2b8_4 - 1._real_8
                vexc11=-1.333333333333_real_8/two13*b1*yy1*ee1/(dd1*dd1)
                vexc12=-1.333333333333_real_8/two13*b1*yy2*ee2/(dd2*dd2)
                vexc13=-1.333333333333_real_8/two13*b1*yy3*ee3/(dd3*dd3)
                vexc14=-1.333333333333_real_8/two13*b1*yy4*ee4/(dd4*dd4)
                vexc21=two13*b1*(ee1-dd1)/(dd1*dd1)
                vexc22=two13*b1*(ee2-dd2)/(dd2*dd2)
                vexc23=two13*b1*(ee3-dd3)/(dd3*dd3)
                vexc24=two13*b1*(ee4-dd4)/(dd4*dd4)
                ! ..correlation lyp
                r1=1._real_8/x1
                r2=1._real_8/x2
                r3=1._real_8/x3
                r4=1._real_8/x4
                r21=r1*r1
                r22=r2*r2
                r23=r3*r3
                r24=r4*r4
                r41=r21*r21
                r42=r22*r22
                r43=r23*r23
                r44=r24*r24
                r51=r41*r1
                r52=r42*r2
                r53=r43*r3
                r54=r44*r4
                om1=EXP(-c*r1)/(1._real_8+d*r1)
                om2=EXP(-c*r2)/(1._real_8+d*r2)
                om3=EXP(-c*r3)/(1._real_8+d*r3)
                om4=EXP(-c*r4)/(1._real_8+d*r4)
                xl1=1._real_8+2.3333333333333_real_8*(c*r1+d*r1/(1._real_8+d*r1))
                xl2=1._real_8+2.3333333333333_real_8*(c*r2+d*r2/(1._real_8+d*r2))
                xl3=1._real_8+2.3333333333333_real_8*(c*r3+d*r3/(1._real_8+d*r3))
                xl4=1._real_8+2.3333333333333_real_8*(c*r4+d*r4/(1._real_8+d*r4))
                ff1=0.04166666666667_real_8*a*b*grho1
                ff2=0.04166666666667_real_8*a*b*grho2
                ff3=0.04166666666667_real_8*a*b*grho3
                ff4=0.04166666666667_real_8*a*b*grho4
                eecc_1=ff1*r51*om1*xl1
                eecc_2=ff2*r52*om2*xl2
                eecc_3=ff3*r53*om3*xl3
                eecc_4=ff4*r54*om4*xl4

                dom1=-om1*(c+d+c*d*r1)/(1._real_8+d*r1)
                dom2=-om2*(c+d+c*d*r2)/(1._real_8+d*r2)
                dom3=-om3*(c+d+c*d*r3)/(1._real_8+d*r3)
                dom4=-om4*(c+d+c*d*r4)/(1._real_8+d*r4)
                dxl1=2.3333333333333_real_8*(c+d+2._real_8*c*d*r1+c*d*d*r21)/&
                     (1._real_8+d*r1)**2
                dxl2=2.3333333333333_real_8*(c+d+2._real_8*c*d*r2+c*d*d*r22)/&
                     (1._real_8+d*r2)**2
                dxl3=2.3333333333333_real_8*(c+d+2._real_8*c*d*r3+c*d*d*r23)/&
                     (1._real_8+d*r3)**2
                dxl4=2.3333333333333_real_8*(c+d+2._real_8*c*d*r4+c*d*d*r24)/&
                     (1._real_8+d*r4)**2

                v1cc_1=-0.3333333333333_real_8*ff1*r41*(5._real_8*r41*om1*xl1+&
                     r51*dom1*xl1+r51*om1*dxl1)
                v1cc_2=-0.3333333333333_real_8*ff2*r42*(5._real_8*r42*om2*xl2+&
                     r52*dom2*xl2+r52*om2*dxl2)
                v1cc_3=-0.3333333333333_real_8*ff3*r43*(5._real_8*r43*om3*xl3+&
                     r53*dom3*xl3+r53*om3*dxl3)
                v1cc_4=-0.3333333333333_real_8*ff4*r44*(5._real_8*r44*om4*xl4+&
                     r54*dom4*xl4+r54*om4*dxl4)
                v2cc_1=0.08333333333333_real_8*a*b*r51*om1*xl1
                v2cc_2=0.08333333333333_real_8*a*b*r52*om2*xl2
                v2cc_3=0.08333333333333_real_8*a*b*r53*om3*xl3
                v2cc_4=0.08333333333333_real_8*a*b*r54*om4*xl4
                ! ..sum up results
                sgcx=sgcx+(smoo1*eexc_1*rho431+smoo2*eexc_2*rho432+&
                     smoo3*eexc_3*rho433+smoo4*eexc_4*rho434)
                sgcc=sgcc+(smoo1*eecc_1+smoo2*eecc_2+&
                     smoo3*eecc_3+smoo4*eecc_4)
                v(i1,i2+0,i3+0)=dsmoo1*eexc_1*rho431+smoo1*vexc11*x1
                v(i1,i2+1,i3+0)=dsmoo2*eexc_2*rho432+smoo2*vexc12*x2
                v(i1,i2+0,i3+1)=dsmoo3*eexc_3*rho433+smoo3*vexc13*x3
                v(i1,i2+1,i3+1)=dsmoo4*eexc_4*rho434+smoo4*vexc14*x4
                v(i1,i2+0,i3+0)=v(i1,i2+0,i3+0)+dsmoo1*eecc_1+smoo1*v1cc_1
                v(i1,i2+1,i3+0)=v(i1,i2+1,i3+0)+dsmoo2*eecc_2+smoo2*v1cc_2
                v(i1,i2+0,i3+1)=v(i1,i2+0,i3+1)+dsmoo3*eecc_3+smoo3*v1cc_3
                v(i1,i2+1,i3+1)=v(i1,i2+1,i3+1)+dsmoo4*eecc_4+smoo4*v1cc_4
                vtmp(i1,i2+0,i3+0)=smoo1*(vexc21/rho431+v2cc_1)
                vtmp(i1,i2+1,i3+0)=smoo2*(vexc22/rho432+v2cc_2)
                vtmp(i1,i2+0,i3+1)=smoo3*(vexc23/rho433+v2cc_3)
                vtmp(i1,i2+1,i3+1)=smoo4*(vexc24/rho434+v2cc_4)
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    ! ..counting floating point operations
    fadd   = ic * 11._real_8 + 4._real_8*id * 21._real_8
    fmult  = ic * 12._real_8 + 4._real_8*id * 68._real_8
    fdiv   = ic *  4._real_8 + 4._real_8*id * 12._real_8
    fspez  =              4._real_8*id *  4._real_8
    flops=fadd+fmult+fdiv+fspez
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gcxlyp
  ! ==================================================================
  SUBROUTINE gcxp86(b1,sgcx,sgcc,rhoe,v,vtmp,grad,flops)
    ! ==--------------------------------------------------------------==
    REAL(real_8) :: b1, sgcx, sgcc, rhoe(fpar%kr1,fpar%kr2s,fpar%kr3s)
    COMPLEX(real_8) :: v(fpar%kr1,fpar%kr2s,fpar%kr3s)
    REAL(real_8) :: vtmp(fpar%kr1,fpar%kr2s,fpar%kr3s), &
      grad(fpar%kr1,fpar%kr2s,fpar%kr3s,*), flops

    REAL(real_8), PARAMETER :: eps2 = 1.e-20_real_8, &
      fp86 = -0.206783495048375038_real_8, p1 = 0.023266_real_8, &
      p2 = 7.389e-6_real_8, p3 = 8.723_real_8, p4 = 0.472_real_8, &
      pc1 = 0.001667_real_8, pc2 = 0.002568_real_8, pci = pc1+pc2, &
      rsfac = 0.6203504908994000_real_8, small = 1.e-24_real_8

    INTEGER                                  :: i1, i2, i3, ic, id, ii2, ii3, &
                                                ncount
    REAL(real_8) :: cn1, cn2, cn3, cn4, cna1, cna2, cna3, cna4, cnb1, cnb2, &
      cnb3, cnb4, d1, d2, d3, d4, dcn1, dcn2, dcn3, dcn4, dcna1, dcna2, &
      dcna3, dcna4, dcnb1, dcnb2, dcnb3, dcnb4, dd1, dd2, dd3, dd4, drho, &
      drs1, drs2, drs3, drs4, dsmoo1, dsmoo2, dsmoo3, dsmoo4, ee1, ee2, ee3, &
      ee4, eexc_1, eexc_2, eexc_3, eexc_4, ephi1, ephi2, ephi3, ephi4, fadd, &
      fdiv, fmult, fspez, grho1, grho2, grho3, grho4, ogceps, orho431, &
      orho432, orho433, orho434, phi1, phi2, phi3, phi4, rho431, rho432, &
      rho433, rho434, rhou1, rhou2, rhou3, rhou4, rs1, rs2, rs3, rs4, &
      sa2b8_1, sa2b8_2, sa2b8_3, sa2b8_4, sc1, sc2, sc3, sc4
    REAL(real_8) :: shm1_1, shm1_2, shm1_3, shm1_4, smoo1, smoo2, smoo3, &
      smoo4, two13, v1c1, v1c2, v1c3, v1c4, v2c1, v2c2, v2c3, v2c4, vexc11, &
      vexc12, vexc13, vexc14, vexc21, vexc22, vexc23, vexc24, x1, x2, x3, x4, &
      xx1, xx2, xx3, xx4, yy1, yy2, yy3, yy4

    two13=2._real_8**0.333333333333333_real_8
    sgcx = 0.0_real_8
    sgcc = 0.0_real_8
    ogceps=1._real_8/cntr%gceps
    x1=1._real_8
    x2=1._real_8
    x3=1._real_8
    x4=1._real_8
    ic=0
    id=0
    ii3=MOD(parm%nr3,2)
    ii2=MOD(parm%nr2,2)
    IF (ii3.NE.0.OR.ii2.NE.0) CALL stopgm('GCXP86','ODD DIMENSION',& 
         __LINE__,__FILE__)
    !$omp parallel do private(I1,I2,I3,RHOU1,RHOU2,RHOU3,RHOU4,D1,D2,D3,D4) &
    !$omp  firstprivate(X1,X2,X3,X4) &
    !$omp  private(NCOUNT,RHO431,RHO432,RHO433,RHO434,ORHO431,ORHO432) &
    !$omp  private(ORHO433,ORHO434,SMOO1,DSMOO1,DRHO,SMOO2,DSMOO2) &
    !$omp  private(SMOO3,DSMOO3,SMOO4,DSMOO4,GRHO1,GRHO2,GRHO3,GRHO4) &
    !$omp  private(XX1,XX2,XX3,XX4,YY1,YY2,YY3,YY4,SA2B8_1,SA2B8_2) &
    !$omp  private(SA2B8_3,SA2B8_4,SHM1_1,SHM1_2,SHM1_3,SHM1_4,DD1,DD2) &
    !$omp  private(DD3,DD4,EEXC_1,EEXC_2,EEXC_3,EEXC_4,EE1,EE2,EE3,EE4) &
    !$omp  private(VEXC11,VEXC12,VEXC13,VEXC14,VEXC21,VEXC22,VEXC23,VEXC24) &
    !$omp  private(RS1,RS2,RS3,RS4,CNA1,CNB1,CNA2,CNB2,CNA3,CNB3,CNA4) &
    !$omp  private(CNB4,CN1,DRS1,DCNA1,DCNB1,CN2,DRS2,DCNA2,DCNB2,CN3) &
    !$omp  private(DRS3,DCNA3,DCNB3,CN4,DRS4,DCNA4,DCNB4,DCN1,PHI1,DCN2) &
    !$omp  private(PHI2,DCN3,PHI3,DCN4,PHI4,EPHI1,EPHI2,EPHI3,EPHI4) &
    !$omp  private(SC1,SC2,SC3,SC4,V1C1,V2C1,V1C2,V2C2,V1C3,V2C3,V1C4,V2C4) &
    !$omp  reduction(+:SGCX,SGCC,IC,ID) __COLLAPSE3
    DO i3=1,parm%nr3,2
       DO i2=1,parm%nr2,2
          DO i1=1,parm%nr1
             rhou1=MAX(rhoe(i1,i2+0,i3+0),small)
             rhou2=MAX(rhoe(i1,i2+1,i3+0),small)
             rhou3=MAX(rhoe(i1,i2+0,i3+1),small)
             rhou4=MAX(rhoe(i1,i2+1,i3+1),small)
             IF ((rhou1+rhou2+rhou3+rhou4).GT.5._real_8*small) THEN
                id=id+1
                DO ncount=1,100
                   d1=x1-rhou1/x1**2
                   d2=x2-rhou2/x2**2
                   d3=x3-rhou3/x3**2
                   d4=x4-rhou4/x4**2
                   ic=ic+4
                   x1=x1-0.333333333333333_real_8*d1
                   x2=x2-0.333333333333333_real_8*d2
                   x3=x3-0.333333333333333_real_8*d3
                   x4=x4-0.333333333333333_real_8*d4
                   IF (d1**2+d2**2+d3**2+d4**2.LT.eps2) GOTO 10
                ENDDO
                x1=rhou1**(1._real_8/3._real_8)
                x2=rhou2**(1._real_8/3._real_8)
                x3=rhou3**(1._real_8/3._real_8)
                x4=rhou4**(1._real_8/3._real_8)
10              CONTINUE
                rho431=x1*rhou1
                rho432=x2*rhou2
                rho433=x3*rhou3
                rho434=x4*rhou4
                orho431=1._real_8/rho431
                orho432=1._real_8/rho432
                orho433=1._real_8/rho433
                orho434=1._real_8/rho434
                IF (rhou1.GT.2.25_real_8*cntr%gceps) THEN
                   smoo1=1._real_8
                   dsmoo1=0._real_8
                ELSEIF (rhou1.LT.0.25_real_8*cntr%gceps) THEN
                   smoo1=0._real_8
                   dsmoo1=0._real_8
                ELSE
                   drho=(rhou1-0.25_real_8*cntr%gceps)*ogceps
                   smoo1=drho*drho*(0.75_real_8-0.25_real_8*drho)
                   dsmoo1=1.5_real_8*drho*ogceps*(1._real_8-0.5_real_8*drho)
                ENDIF
                IF (rhou2.GT.2.25_real_8*cntr%gceps) THEN
                   smoo2=1._real_8
                   dsmoo2=0._real_8
                ELSEIF (rhou2.LT.0.25_real_8*cntr%gceps) THEN
                   smoo2=0._real_8
                   dsmoo2=0._real_8
                ELSE
                   drho=(rhou2-0.25_real_8*cntr%gceps)*ogceps
                   smoo2=drho*drho*(0.75_real_8-0.25_real_8*drho)
                   dsmoo2=1.5_real_8*drho*ogceps*(1._real_8-0.5_real_8*drho)
                ENDIF
                IF (rhou3.GT.2.25_real_8*cntr%gceps) THEN
                   smoo3=1._real_8
                   dsmoo3=0._real_8
                ELSEIF (rhou3.LT.0.25_real_8*cntr%gceps) THEN
                   smoo3=0._real_8
                   dsmoo3=0._real_8
                ELSE
                   drho=(rhou3-0.25_real_8*cntr%gceps)*ogceps
                   smoo3=drho*drho*(0.75_real_8-0.25_real_8*drho)
                   dsmoo3=1.5_real_8*drho*ogceps*(1._real_8-0.5_real_8*drho)
                ENDIF
                IF (rhou4.GT.2.25_real_8*cntr%gceps) THEN
                   smoo4=1._real_8
                   dsmoo4=0._real_8
                ELSEIF (rhou4.LT.0.25_real_8*cntr%gceps) THEN
                   smoo4=0._real_8
                   dsmoo4=0._real_8
                ELSE
                   drho=(rhou4-0.25_real_8*cntr%gceps)*ogceps
                   smoo4=drho*drho*(0.75_real_8-0.25_real_8*drho)
                   dsmoo4=1.5_real_8*drho*ogceps*(1._real_8-0.5_real_8*drho)
                ENDIF
                grho1=SQRT(grad(i1,i2+0,i3+0,1))
                grho2=SQRT(grad(i1,i2+1,i3+0,1))
                grho3=SQRT(grad(i1,i2+0,i3+1,1))
                grho4=SQRT(grad(i1,i2+1,i3+1,1))
                ! ..exchange
                xx1=two13*grho1/rho431
                xx2=two13*grho2/rho432
                xx3=two13*grho3/rho433
                xx4=two13*grho4/rho434
                yy1=xx1*xx1
                yy2=xx2*xx2
                yy3=xx3*xx3
                yy4=xx4*xx4
                sa2b8_1=SQRT(1._real_8+yy1)
                sa2b8_2=SQRT(1._real_8+yy2)
                sa2b8_3=SQRT(1._real_8+yy3)
                sa2b8_4=SQRT(1._real_8+yy4)
                shm1_1=LOG(xx1+sa2b8_1)
                shm1_2=LOG(xx2+sa2b8_2)
                shm1_3=LOG(xx3+sa2b8_3)
                shm1_4=LOG(xx4+sa2b8_4)
                dd1=1._real_8+6._real_8*b1*xx1*shm1_1
                dd2=1._real_8+6._real_8*b1*xx2*shm1_2
                dd3=1._real_8+6._real_8*b1*xx3*shm1_3
                dd4=1._real_8+6._real_8*b1*xx4*shm1_4

                eexc_1=-b1*yy1/dd1/two13
                eexc_2=-b1*yy2/dd2/two13
                eexc_3=-b1*yy3/dd3/two13
                eexc_4=-b1*yy4/dd4/two13

                ee1=6._real_8*b1*yy1/sa2b8_1 - 1._real_8
                ee2=6._real_8*b1*yy2/sa2b8_2 - 1._real_8
                ee3=6._real_8*b1*yy3/sa2b8_3 - 1._real_8
                ee4=6._real_8*b1*yy4/sa2b8_4 - 1._real_8
                vexc11=-1.333333333333_real_8/two13*b1*yy1*ee1/(dd1*dd1)
                vexc12=-1.333333333333_real_8/two13*b1*yy2*ee2/(dd2*dd2)
                vexc13=-1.333333333333_real_8/two13*b1*yy3*ee3/(dd3*dd3)
                vexc14=-1.333333333333_real_8/two13*b1*yy4*ee4/(dd4*dd4)
                vexc21=two13*b1*(ee1-dd1)/(dd1*dd1)
                vexc22=two13*b1*(ee2-dd2)/(dd2*dd2)
                vexc23=two13*b1*(ee3-dd3)/(dd3*dd3)
                vexc24=two13*b1*(ee4-dd4)/(dd4*dd4)
                ! .. perdew 86 correlation
                xx1=grho1*orho431
                xx2=grho2*orho432
                xx3=grho3*orho433
                xx4=grho4*orho434
                rs1   = rsfac/x1
                rs2   = rsfac/x2
                rs3   = rsfac/x3
                rs4   = rsfac/x4
                cna1  = pc2+rs1*(p1+p2*rs1)
                cnb1  = 1._real_8+rs1*(p3+rs1*(p4+1.e4_real_8*p2*rs1))
                cna2  = pc2+rs2*(p1+p2*rs2)
                cnb2  = 1._real_8+rs2*(p3+rs2*(p4+1.e4_real_8*p2*rs2))
                cna3  = pc2+rs3*(p1+p2*rs3)
                cnb3  = 1._real_8+rs3*(p3+rs3*(p4+1.e4_real_8*p2*rs3))
                cna4  = pc2+rs4*(p1+p2*rs4)
                cnb4  = 1._real_8+rs4*(p3+rs4*(p4+1.e4_real_8*p2*rs4))

                cn1   = pc1 + cna1/cnb1
                drs1  = fp86*orho431
                dcna1 = (p1+2._real_8*p2*rs1)*drs1
                dcnb1 = (p3+rs1*(2._real_8*p4+3.e4_real_8*p2*rs1))*drs1
                cn2   = pc1 + cna2/cnb2
                drs2  = fp86*orho432
                dcna2 = (p1+2._real_8*p2*rs2)*drs2
                dcnb2 = (p3+rs2*(2._real_8*p4+3.e4_real_8*p2*rs2))*drs2
                cn3   = pc1 + cna3/cnb3
                drs3  = fp86*orho433
                dcna3 = (p1+2._real_8*p2*rs3)*drs3
                dcnb3 = (p3+rs3*(2._real_8*p4+3.e4_real_8*p2*rs3))*drs3
                cn4   = pc1 + cna4/cnb4
                drs4  = fp86*orho434
                dcna4 = (p1+2._real_8*p2*rs4)*drs4
                dcnb4 = (p3+rs4*(2._real_8*p4+3.e4_real_8*p2*rs4))*drs4

                dcn1  = dcna1/cnb1 - cna1/(cnb1*cnb1)*dcnb1
                phi1  = 0.192_real_8*pci/cn1*grho1*orho431*SQRT(x1)
                dcn2  = dcna2/cnb2 - cna2/(cnb2*cnb2)*dcnb2
                phi2  = 0.192_real_8*pci/cn2*grho2*orho432*SQRT(x2)
                dcn3  = dcna3/cnb3 - cna3/(cnb3*cnb3)*dcnb3
                phi3  = 0.192_real_8*pci/cn3*grho3*orho433*SQRT(x3)
                dcn4  = dcna4/cnb4 - cna4/(cnb4*cnb4)*dcnb4
                phi4  = 0.192_real_8*pci/cn4*grho4*orho434*SQRT(x4)

                ephi1 = EXP(-phi1)
                ephi2 = EXP(-phi2)
                ephi3 = EXP(-phi3)
                ephi4 = EXP(-phi4)
                sc1   = xx1*grho1*cn1*ephi1
                sc2   = xx2*grho2*cn2*ephi2
                sc3   = xx3*grho3*cn3*ephi3
                sc4   = xx4*grho4*cn4*ephi4

                v1c1  = sc1*((1._real_8+phi1)*dcn1/cn1 -&
                     (1.33333333333333_real_8-1.166666666666667_real_8*phi1)/rhou1)
                v2c1  = cn1*ephi1*orho431*(2._real_8-phi1)
                v1c2  = sc2*((1._real_8+phi2)*dcn2/cn2 -&
                     (1.33333333333333_real_8-1.166666666666667_real_8*phi2)/rhou2)
                v2c2  = cn2*ephi2*orho432*(2._real_8-phi2)
                v1c3  = sc3*((1._real_8+phi3)*dcn3/cn3 -&
                     (1.33333333333333_real_8-1.166666666666667_real_8*phi3)/rhou3)
                v2c3  = cn3*ephi3*orho433*(2._real_8-phi3)
                v1c4  = sc4*((1._real_8+phi4)*dcn4/cn4 -&
                     (1.33333333333333_real_8-1.166666666666667_real_8*phi4)/rhou4)
                v2c4  = cn4*ephi4*orho434*(2._real_8-phi4)

                sgcx=sgcx+(smoo1*eexc_1*rho431+smoo2*eexc_2*rho432+&
                     smoo3*eexc_3*rho433+smoo4*eexc_4*rho434)
                sgcc=sgcc+(smoo1*sc1+smoo2*sc2+smoo3*sc3+smoo4*sc4)
                v(i1,i2+0,i3+0)=dsmoo1*(sc1+eexc_1*rho431)+&
                     smoo1*(vexc11*x1+v1c1)
                v(i1,i2+1,i3+0)=dsmoo2*(sc2+eexc_2*rho432)+&
                     smoo2*(vexc12*x2+v1c2)
                v(i1,i2+0,i3+1)=dsmoo3*(sc3+eexc_3*rho433)+&
                     smoo3*(vexc13*x3+v1c3)
                v(i1,i2+1,i3+1)=dsmoo4*(sc4+eexc_4*rho434)+&
                     smoo4*(vexc14*x4+v1c4)
                vtmp(i1,i2+0,i3+0)=smoo1*(vexc21*orho431+v2c1)
                vtmp(i1,i2+1,i3+0)=smoo2*(vexc22*orho432+v2c2)
                vtmp(i1,i2+0,i3+1)=smoo3*(vexc23*orho433+v2c3)
                vtmp(i1,i2+1,i3+1)=smoo4*(vexc24*orho434+v2c4)
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    ! ..counting floating point operations
    fadd   = ic * 11._real_8 + id * 100._real_8
    fmult  = ic * 12._real_8 + id * 224._real_8
    fdiv   = ic *  4._real_8 + id *  54._real_8
    fspez  =              id *  20._real_8
    flops=fadd+fmult+fdiv+fspez
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gcxp86
  ! ==================================================================
  SUBROUTINE gcgga(sgcx,sgcc,rhoe,v,vtmp,grad,flops)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: sgcx, sgcc
    REAL(real_8), DIMENSION(fpar%kr1, fpar%&
      kr2s, fpar%kr3s)                       :: rhoe
    COMPLEX(real_8) :: v(fpar%kr1,fpar%kr2s,fpar%kr3s)
    REAL(real_8), DIMENSION(fpar%kr1, fpar%&
      kr2s, fpar%kr3s)                       :: vtmp
    REAL(real_8), DIMENSION(fpar%kr1, fpar%&
      kr2s, fpar%kr3s, *)                    :: grad
    REAL(real_8)                             :: flops

    REAL(real_8), PARAMETER :: a0u = 0.4581652932831429_real_8, &
      a1u = 2.217058676663745_real_8, a2u = 0.7405551735357053_real_8, &
      a3u = 0.01968227878617998_real_8 , al = 0.09_real_8, b1u = 1.0_real_8, &
      b2u = 4.504130959426697_real_8, b3u = 1.110667363742916_real_8, &
      b4u = 0.02359291751427506_real_8 , c1u = 4._real_8*a0u*b1u/3.0_real_8, &
      c2u = 5.0_real_8*a0u*b2u/3.0_real_8+a1u*b1u, c3u = 2.0_real_8*a0u*&
      b3u+4.0_real_8*a1u*b2u/3.0_real_8+2.0_real_8*a2u*b1u/3.0_real_8, c4u = &
      7.0_real_8*a0u*b4u/3.0_real_8+5.0_real_8*a1u*b3u/3.0_real_8+a2u*b2u+a3u*&
      b1u/3.0_real_8
    REAL(real_8), PARAMETER :: c5u = 2.0_real_8*a1u*b4u+4.0_real_8*a2u*b3u/&
      3.0_real_8+2.0_real_8*a3u*b2u/3.0_real_8, &
      c6u = 5.0_real_8*a2u*b4u/3.0_real_8+a3u*b3u, &
      c7u = 4.0_real_8*a3u*b4u/3.0_real_8 , cx = -0.001667_real_8, &
      cxc0 = 0.002568_real_8, cc0 = -cx+cxc0 , eps2 = 1.e-20_real_8, &
      f1 = 0.19645_real_8, f1x = -1.10783814957303361_real_8 , &
      f2 = 7.7956_real_8, f3 = 0.2743_real_8, f4 = 0.1508_real_8, &
      f5 = 0.004_real_8 , pa = 0.023266_real_8, pb = 7.389e-6_real_8, &
      pc = 8.723_real_8, pd = 0.472_real_8 , &
      rsfac = 0.6203504908994000_real_8 , small = 1.e-24_real_8 

    INTEGER                                  :: i1, i2, i3, ic, id, ii2, ii3, &
                                                ncount
    REAL(real_8) :: af1, af2, af3, af4, as1, as2, as3, as4, be, bf1, bf2, &
      bf3, bf4, botu1, botu2, botu3, botu4, bs1, bs2, bs3, bs4, cn1, cn2, &
      cn3, cn4, cna1, cna2, cna3, cna4, cnb1, cnb2, cnb3, cnb4, d1, d2, d3, &
      d4, das1, das2, das3, das4, dbs1, dbs2, dbs3, dbs4, dcn1, dcn2, dcn3, &
      dcn4, dcna1, dcna2, dcna3, dcna4, dcnb1, dcnb2, dcnb3, dcnb4, ddh01, &
      ddh02, ddh03, ddh04, ddh11, ddh12, ddh13, ddh14, dh01, dh02, dh03, &
      dh04, dh11, dh12, dh13, dh14, dls1, dls2, dls3, dls4, drho, dsmoo1, &
      dsmoo2, dsmoo3, dsmoo4, dtopu1, dtopu2, dtopu3, dtopu4, ee1, ee2, ee3, &
      ee4, epsxcu1, epsxcu2, epsxcu3, epsxcu4
    REAL(real_8) :: expe1, expe2, expe3, expe4, exps1, exps2, exps3, exps4, &
      fa1, fa2, fadd, fdiv, fmult, fp1, fp2, fspez, grho1, grho2, grho3, &
      grho4, grhos1, grhos2, grhos3, grhos4, h01, h02, h03, h04, h11, h12, &
      h13, h14, ogceps, qy1, qy2, qy3, qy4, rho431, rho432, rho433, rho434, &
      rhou1, rhou2, rhou3, rhou4, rr1, rr2, rr3, rr4, rs21, rs22, rs23, rs24, &
      rs31, rs32, rs33, rs34, rsu1, rsu2, rsu3, rsu4, s1, s11, s12, s13, s14, &
      s2, s21, s22, s23, s24, s3, s31, s32, s33, s34, s4, s41, s42, s43, s44, &
      sa2b8_1, sa2b8_2, sa2b8_3, sa2b8_4, sc1, sc2, sc3, sc4, shm11, shm12, &
      shm13, shm14, smoo1, smoo2, smoo3
    REAL(real_8) :: smoo4, sx1, sx2, sx3, sx4, t1, t2, t3, t4, topu1, topu2, &
      topu3, topu4, v1c1, v1c2, v1c3, v1c4, v1x1, v1x2, v1x3, v1x4, v2c1, &
      v2c2, v2c3, v2c4, v2x1, v2x2, v2x3, v2x4, vxc1, vxc2, vxc3, vxc4, x1, &
      x2, x3, x4, xkf1, xkf2, xkf3, xkf4, xkff, xks1, xks2, xks3, xks4, xnu, &
      xy1, xy2, xy3, xy4, y1, y2, y3, y4

! ==--------------------------------------------------------------==

    fa1=f1x*2._real_8/3._real_8
    fa2=fa1*4._real_8/3._real_8
    fp1   = -3._real_8/(16._real_8*pi)*(3._real_8*pi*pi)**(-1._real_8/3._real_8)
    fp2   = 0.5_real_8*(3._real_8*pi*pi)**(-1._real_8/3._real_8)
    xnu   = 16._real_8/pi*(3._real_8*pi*pi)**0.333333333333333_real_8
    be    = xnu*cc0
    xkff  = (9._real_8*pi/4._real_8)**0.333333333333333_real_8
    sgcx = 0.0_real_8
    sgcc = 0.0_real_8
    ogceps=1._real_8/cntr%gceps
    x1=1._real_8
    x2=1._real_8
    x3=1._real_8
    x4=1._real_8
    ic=0
    id=0
    ii3=MOD(parm%nr3,2)
    ii2=MOD(parm%nr2,2)
    IF (ii3.NE.0.OR.ii2.NE.0) CALL stopgm('GCGGA','ODD DIMENSION',& 
         __LINE__,__FILE__)
    !$omp parallel do private(I1,I2,I3,RHOU1,RHOU2,RHOU3,RHOU4,D1,D2,D3,D4) &
    !$omp  firstprivate(X1,X2,X3,X4) &
    !$omp  private(NCOUNT,RHO431,RHO432,RHO433,RHO434,SMOO1,DSMOO1) &
    !$omp  private(DRHO,SMOO2,DSMOO2,SMOO3,DSMOO3,SMOO4,DSMOO4) &
    !$omp  private(GRHO1,GRHO2,GRHO3,GRHO4,GRHOS1,GRHOS2,GRHOS3) &
    !$omp  private(GRHOS4,RR1,RR2,RR3,RR4,S11,S12,S13,S14,S21) &
    !$omp  private(S22,S23,S24,S31,S32,S33,S34,S41,S42,S43,S44) &
    !$omp  private(EXPS1,EXPS2,EXPS3,EXPS4,AS1,AS2,AS3,AS4,SA2B8_1) &
    !$omp  private(SA2B8_2,SA2B8_3,SA2B8_4,SHM11,SHM12,SHM13,SHM14) &
    !$omp  private(BS1,BS2,BS3,BS4,DAS1,DAS2,DAS3,DAS4,DBS1,DBS2) &
    !$omp  private(DBS3,DBS4,DLS1,DLS2,DLS3,DLS4,SX1,SX2,SX3,SX4) &
    !$omp  private(V1X1,V1X2,V1X3,V1X4,V2X1,V2X2,V2X3,V2X4,RSU1) &
    !$omp  private(RSU2,RSU3,RSU4,TOPU1,TOPU2,TOPU3,TOPU4,DTOPU1) &
    !$omp  private(DTOPU2,DTOPU3,DTOPU4,BOTU1,BOTU2,BOTU3,BOTU4) &
    !$omp  private(EPSXCU1,EPSXCU2,EPSXCU3,EPSXCU4,VXC1,VXC2,VXC3) &
    !$omp  private(VXC4,RS21,RS22,RS23,RS24,RS31,RS32,RS33,RS34) &
    !$omp  private(XKF1,XKF2,XKF3,XKF4,XKS1,XKS2,XKS3,XKS4,T1,T2) &
    !$omp  private(T3,T4,EXPE1,EXPE2,EXPE3,EXPE4,AF1,AF2,AF3,AF4) &
    !$omp  private(BF1,BF2,BF3,BF4,Y1,Y2,Y3,Y4,XY1,XY2,XY3,XY4) &
    !$omp  private(QY1,QY2,QY3,QY4,S1,S2,S3,S4,H01,H02,H03,H04,DH01) &
    !$omp  private(DH02,DH03,DH04,DDH01,DDH02,DDH03,DDH04,EE1,EE2) &
    !$omp  private(EE3,EE4,CNA1,CNA2,CNA3,CNA4,DCNA1,DCNA2,DCNA3,DCNA4) &
    !$omp  private(CNB1,CNB2,CNB3,CNB4,DCNB1,DCNB2,DCNB3,DCNB4,CN1,CN2) &
    !$omp  private(CN3,CN4,DCN1,DCN2,DCN3,DCN4,H11,H12,H13,H14,DH11) &
    !$omp  private(DH12,DH13,DH14,DDH11,DDH12,DDH13,DDH14,SC1,SC2) &
    !$omp  private(SC3,SC4,V1C1,V1C2,V1C3,V1C4,V2C1,V2C2,V2C3,V2C4) &
    !$omp  reduction(+:SGCX,SGCC,IC,ID) __COLLAPSE3
    DO i3=1,parm%nr3,2
       DO i2=1,parm%nr2,2
          DO i1=1,parm%nr1
             rhou1=MAX(rhoe(i1,i2+0,i3+0),small)
             rhou2=MAX(rhoe(i1,i2+1,i3+0),small)
             rhou3=MAX(rhoe(i1,i2+0,i3+1),small)
             rhou4=MAX(rhoe(i1,i2+1,i3+1),small)
             IF ((rhou1+rhou2+rhou3+rhou4).GT.5._real_8*small) THEN
                id=id+1
                DO ncount=1,100
                   d1=x1-rhou1/x1**2
                   d2=x2-rhou2/x2**2
                   d3=x3-rhou3/x3**2
                   d4=x4-rhou4/x4**2
                   ic=ic+4
                   x1=x1-0.333333333333333_real_8*d1
                   x2=x2-0.333333333333333_real_8*d2
                   x3=x3-0.333333333333333_real_8*d3
                   x4=x4-0.333333333333333_real_8*d4
                   IF (d1**2+d2**2+d3**2+d4**2.LT.eps2) GOTO 10
                ENDDO
                x1=rhou1**(1._real_8/3._real_8)
                x2=rhou2**(1._real_8/3._real_8)
                x3=rhou3**(1._real_8/3._real_8)
                x4=rhou4**(1._real_8/3._real_8)
10              CONTINUE
                rho431=x1*rhou1
                rho432=x2*rhou2
                rho433=x3*rhou3
                rho434=x4*rhou4
                IF (rhou1.GT.2.25_real_8*cntr%gceps) THEN
                   smoo1=1._real_8
                   dsmoo1=0._real_8
                ELSEIF (rhou1.LT.0.25_real_8*cntr%gceps) THEN
                   smoo1=0._real_8
                   dsmoo1=0._real_8
                ELSE
                   drho=(rhou1-0.25_real_8*cntr%gceps)*ogceps
                   smoo1=drho*drho*(0.75_real_8-0.25_real_8*drho)
                   dsmoo1=1.5_real_8*drho*ogceps*(1._real_8-0.5_real_8*drho)
                ENDIF
                IF (rhou2.GT.2.25_real_8*cntr%gceps) THEN
                   smoo2=1._real_8
                   dsmoo2=0._real_8
                ELSEIF (rhou2.LT.0.25_real_8*cntr%gceps) THEN
                   smoo2=0._real_8
                   dsmoo2=0._real_8
                ELSE
                   drho=(rhou2-0.25_real_8*cntr%gceps)*ogceps
                   smoo2=drho*drho*(0.75_real_8-0.25_real_8*drho)
                   dsmoo2=1.5_real_8*drho*ogceps*(1._real_8-0.5_real_8*drho)
                ENDIF
                IF (rhou3.GT.2.25_real_8*cntr%gceps) THEN
                   smoo3=1._real_8
                   dsmoo3=0._real_8
                ELSEIF (rhou3.LT.0.25_real_8*cntr%gceps) THEN
                   smoo3=0._real_8
                   dsmoo3=0._real_8
                ELSE
                   drho=(rhou3-0.25_real_8*cntr%gceps)*ogceps
                   smoo3=drho*drho*(0.75_real_8-0.25_real_8*drho)
                   dsmoo3=1.5_real_8*drho*ogceps*(1._real_8-0.5_real_8*drho)
                ENDIF
                IF (rhou4.GT.2.25_real_8*cntr%gceps) THEN
                   smoo4=1._real_8
                   dsmoo4=0._real_8
                ELSEIF (rhou4.LT.0.25_real_8*cntr%gceps) THEN
                   smoo4=0._real_8
                   dsmoo4=0._real_8
                ELSE
                   drho=(rhou4-0.25_real_8*cntr%gceps)*ogceps
                   smoo4=drho*drho*(0.75_real_8-0.25_real_8*drho)
                   dsmoo4=1.5_real_8*drho*ogceps*(1._real_8-0.5_real_8*drho)
                ENDIF
                grho1=grad(i1,i2+0,i3+0,1)
                grho2=grad(i1,i2+1,i3+0,1)
                grho3=grad(i1,i2+0,i3+1,1)
                grho4=grad(i1,i2+1,i3+1,1)
                grhos1=SQRT(grho1)
                grhos2=SQRT(grho2)
                grhos3=SQRT(grho3)
                grhos4=SQRT(grho4)
                ! ..exchange part
                rr1=1._real_8/rho431
                rr2=1._real_8/rho432
                rr3=1._real_8/rho433
                rr4=1._real_8/rho434
                s11=fp2*grhos1*rr1
                s12=fp2*grhos2*rr2
                s13=fp2*grhos3*rr3
                s14=fp2*grhos4*rr4
                s21=s11*s11
                s22=s12*s12
                s23=s13*s13
                s24=s14*s14
                s31=s11*s21
                s32=s12*s22
                s33=s13*s23
                s34=s14*s24
                s41=s21*s21
                s42=s22*s22
                s43=s23*s23
                s44=s24*s24
                exps1=f4*EXP(-100._real_8*s21)
                exps2=f4*EXP(-100._real_8*s22)
                exps3=f4*EXP(-100._real_8*s23)
                exps4=f4*EXP(-100._real_8*s24)
                as1=f3-exps1-f5*s21
                as2=f3-exps2-f5*s22
                as3=f3-exps3-f5*s23
                as4=f3-exps4-f5*s24
                sa2b8_1=SQRT(1._real_8+f2*f2*s21)
                sa2b8_2=SQRT(1._real_8+f2*f2*s22)
                sa2b8_3=SQRT(1._real_8+f2*f2*s23)
                sa2b8_4=SQRT(1._real_8+f2*f2*s24)
                shm11=LOG(f2*s11+sa2b8_1)
                shm12=LOG(f2*s12+sa2b8_2)
                shm13=LOG(f2*s13+sa2b8_3)
                shm14=LOG(f2*s14+sa2b8_4)
                bs1=1._real_8+f1*s11*shm11+f5*s41
                bs2=1._real_8+f1*s12*shm12+f5*s42
                bs3=1._real_8+f1*s13*shm13+f5*s43
                bs4=1._real_8+f1*s14*shm14+f5*s44
                das1=(200._real_8*exps1-2._real_8*f5)*s11
                das2=(200._real_8*exps2-2._real_8*f5)*s12
                das3=(200._real_8*exps3-2._real_8*f5)*s13
                das4=(200._real_8*exps4-2._real_8*f5)*s14
                dbs1=f1*(shm11+f2*s11/sa2b8_1)+4._real_8*f5*s31
                dbs2=f1*(shm12+f2*s12/sa2b8_2)+4._real_8*f5*s32
                dbs3=f1*(shm13+f2*s13/sa2b8_3)+4._real_8*f5*s33
                dbs4=f1*(shm14+f2*s14/sa2b8_4)+4._real_8*f5*s34
                dls1=(das1/as1-dbs1/bs1)
                dls2=(das2/as2-dbs2/bs2)
                dls3=(das3/as3-dbs3/bs3)
                dls4=(das4/as4-dbs4/bs4)
                sx1=fp1*grho1*rr1*as1/bs1
                sx2=fp1*grho2*rr2*as2/bs2
                sx3=fp1*grho3*rr3*as3/bs3
                sx4=fp1*grho4*rr4*as4/bs4
                v1x1=-1.333333333333333_real_8*sx1/rhou1*(1._real_8+s11*dls1)
                v1x2=-1.333333333333333_real_8*sx2/rhou2*(1._real_8+s12*dls2)
                v1x3=-1.333333333333333_real_8*sx3/rhou3*(1._real_8+s13*dls3)
                v1x4=-1.333333333333333_real_8*sx4/rhou4*(1._real_8+s14*dls4)
                v2x1=fp1*rr1*as1/bs1*(2._real_8+s11*dls1)
                v2x2=fp1*rr2*as2/bs2*(2._real_8+s12*dls2)
                v2x3=fp1*rr3*as3/bs3*(2._real_8+s13*dls3)
                v2x4=fp1*rr4*as4/bs4*(2._real_8+s14*dls4)
                ! ..lda functional
                rsu1=rsfac/x1
                rsu2=rsfac/x2
                rsu3=rsfac/x3
                rsu4=rsfac/x4
                topu1=a2u+rsu1*a3u
                topu2=a2u+rsu2*a3u
                topu3=a2u+rsu3*a3u
                topu4=a2u+rsu4*a3u
                topu1=a1u+rsu1*topu1
                topu2=a1u+rsu2*topu2
                topu3=a1u+rsu3*topu3
                topu4=a1u+rsu4*topu4
                topu1=a0u+rsu1*topu1
                topu2=a0u+rsu2*topu2
                topu3=a0u+rsu3*topu3
                topu4=a0u+rsu4*topu4
                dtopu1=c6u+rsu1*c7u
                dtopu2=c6u+rsu2*c7u
                dtopu3=c6u+rsu3*c7u
                dtopu4=c6u+rsu4*c7u
                dtopu1=c5u+rsu1*dtopu1
                dtopu2=c5u+rsu2*dtopu2
                dtopu3=c5u+rsu3*dtopu3
                dtopu4=c5u+rsu4*dtopu4
                dtopu1=c4u+rsu1*dtopu1
                dtopu2=c4u+rsu2*dtopu2
                dtopu3=c4u+rsu3*dtopu3
                dtopu4=c4u+rsu4*dtopu4
                dtopu1=c3u+rsu1*dtopu1
                dtopu2=c3u+rsu2*dtopu2
                dtopu3=c3u+rsu3*dtopu3
                dtopu4=c3u+rsu4*dtopu4
                dtopu1=c2u+rsu1*dtopu1
                dtopu2=c2u+rsu2*dtopu2
                dtopu3=c2u+rsu3*dtopu3
                dtopu4=c2u+rsu4*dtopu4
                dtopu1=c1u+rsu1*dtopu1
                dtopu2=c1u+rsu2*dtopu2
                dtopu3=c1u+rsu3*dtopu3
                dtopu4=c1u+rsu4*dtopu4
                dtopu1=-rsu1*dtopu1
                dtopu2=-rsu2*dtopu2
                dtopu3=-rsu3*dtopu3
                dtopu4=-rsu4*dtopu4
                botu1=b3u+rsu1*b4u
                botu2=b3u+rsu2*b4u
                botu3=b3u+rsu3*b4u
                botu4=b3u+rsu4*b4u
                botu1=b2u+rsu1*botu1
                botu2=b2u+rsu2*botu2
                botu3=b2u+rsu3*botu3
                botu4=b2u+rsu4*botu4
                botu1=b1u+rsu1*botu1
                botu2=b1u+rsu2*botu2
                botu3=b1u+rsu3*botu3
                botu4=b1u+rsu4*botu4
                botu1=rsu1*botu1
                botu2=rsu2*botu2
                botu3=rsu3*botu3
                botu4=rsu4*botu4
                epsxcu1=topu1*(-1._real_8/botu1) - fa1*x1
                epsxcu2=topu2*(-1._real_8/botu2) - fa1*x2
                epsxcu3=topu3*(-1._real_8/botu3) - fa1*x3
                epsxcu4=topu4*(-1._real_8/botu4) - fa1*x4
                vxc1=dtopu1*(-1._real_8/botu1)*(-1._real_8/botu1) - fa2*x1
                vxc2=dtopu2*(-1._real_8/botu2)*(-1._real_8/botu2) - fa2*x2
                vxc3=dtopu3*(-1._real_8/botu3)*(-1._real_8/botu3) - fa2*x3
                vxc4=dtopu4*(-1._real_8/botu4)*(-1._real_8/botu4) - fa2*x4
                ! ..correlation part
                rs21=rsu1*rsu1
                rs22=rsu2*rsu2
                rs23=rsu3*rsu3
                rs24=rsu4*rsu4
                rs31=rs21*rsu1
                rs32=rs22*rsu2
                rs33=rs23*rsu3
                rs34=rs24*rsu4
                xkf1=xkff/rsu1
                xkf2=xkff/rsu2
                xkf3=xkff/rsu3
                xkf4=xkff/rsu4
                xks1=SQRT(4._real_8*xkf1/pi)
                xks2=SQRT(4._real_8*xkf2/pi)
                xks3=SQRT(4._real_8*xkf3/pi)
                xks4=SQRT(4._real_8*xkf4/pi)
                t1  =grhos1/(2._real_8*xks1*rhou1)
                t2  =grhos2/(2._real_8*xks2*rhou2)
                t3  =grhos3/(2._real_8*xks3*rhou3)
                t4  =grhos4/(2._real_8*xks4*rhou4)
                expe1=EXP(-2._real_8*al*epsxcu1/(be*be))
                expe2=EXP(-2._real_8*al*epsxcu2/(be*be))
                expe3=EXP(-2._real_8*al*epsxcu3/(be*be))
                expe4=EXP(-2._real_8*al*epsxcu4/(be*be))
                af1=2._real_8*al/be*(1._real_8/(expe1-1._real_8))
                af2=2._real_8*al/be*(1._real_8/(expe2-1._real_8))
                af3=2._real_8*al/be*(1._real_8/(expe3-1._real_8))
                af4=2._real_8*al/be*(1._real_8/(expe4-1._real_8))
                bf1=expe1*(vxc1-epsxcu1)
                bf2=expe2*(vxc2-epsxcu2)
                bf3=expe3*(vxc3-epsxcu3)
                bf4=expe4*(vxc4-epsxcu4)
                y1 = af1*t1*t1
                y2 = af2*t2*t2
                y3 = af3*t3*t3
                y4 = af4*t4*t4
                xy1 = (1._real_8+y1)/(1._real_8+y1+y1*y1)
                xy2 = (1._real_8+y2)/(1._real_8+y2+y2*y2)
                xy3 = (1._real_8+y3)/(1._real_8+y3+y3*y3)
                xy4 = (1._real_8+y4)/(1._real_8+y4+y4*y4)
                qy1 = y1*y1*(2._real_8+y1)/(1._real_8+y1+y1*y1)**2
                qy2 = y2*y2*(2._real_8+y2)/(1._real_8+y2+y2*y2)**2
                qy3 = y3*y3*(2._real_8+y3)/(1._real_8+y3+y3*y3)**2
                qy4 = y4*y4*(2._real_8+y4)/(1._real_8+y4+y4*y4)**2
                s1  = 1._real_8+2._real_8*al/be*t1*t1*xy1
                s2  = 1._real_8+2._real_8*al/be*t2*t2*xy2
                s3  = 1._real_8+2._real_8*al/be*t3*t3*xy3
                s4  = 1._real_8+2._real_8*al/be*t4*t4*xy4
                h01 = be*be/(2._real_8*al) * LOG(s1)
                h02 = be*be/(2._real_8*al) * LOG(s2)
                h03 = be*be/(2._real_8*al) * LOG(s3)
                h04 = be*be/(2._real_8*al) * LOG(s4)
                dh01= be*t1*t1/s1*(-2.33333333333333_real_8*xy1-qy1*&
                     (af1*bf1/be-2.3333333333333333_real_8))
                dh02= be*t2*t2/s2*(-2.33333333333333_real_8*xy2-qy2*&
                     (af2*bf2/be-2.3333333333333333_real_8))
                dh03= be*t3*t3/s3*(-2.33333333333333_real_8*xy3-qy3*&
                     (af3*bf3/be-2.3333333333333333_real_8))
                dh04= be*t4*t4/s4*(-2.33333333333333_real_8*xy4-qy4*&
                     (af4*bf4/be-2.3333333333333333_real_8))
                ddh01=be/(2._real_8*xks1*xks1*rhou1)*(xy1-qy1)/s1
                ddh02=be/(2._real_8*xks2*xks2*rhou2)*(xy2-qy2)/s2
                ddh03=be/(2._real_8*xks3*xks3*rhou3)*(xy3-qy3)/s3
                ddh04=be/(2._real_8*xks4*xks4*rhou4)*(xy4-qy4)/s4
                ee1 = -100._real_8*(xks1/xkf1*t1)**2
                ee2 = -200._real_8*(xks2/xkf2*t2)**2
                ee3 = -300._real_8*(xks3/xkf3*t3)**2
                ee4 = -400._real_8*(xks4/xkf4*t4)**2
                cna1= cxc0+pa*rsu1+pb*rs21
                cna2= cxc0+pa*rsu2+pb*rs22
                cna3= cxc0+pa*rsu3+pb*rs23
                cna4= cxc0+pa*rsu4+pb*rs24
                dcna1= -(pa*rsu1+2._real_8*pb*rs21)*0.33333333333333_real_8
                dcna2= -(pa*rsu2+2._real_8*pb*rs22)*0.33333333333333_real_8
                dcna3= -(pa*rsu3+2._real_8*pb*rs23)*0.33333333333333_real_8
                dcna4= -(pa*rsu4+2._real_8*pb*rs24)*0.33333333333333_real_8
                cnb1  = 1._real_8+pc*rsu1+pd*rs21+1.e4_real_8*pb*rs31
                cnb2  = 1._real_8+pc*rsu2+pd*rs22+1.e4_real_8*pb*rs32
                cnb3  = 1._real_8+pc*rsu3+pd*rs23+1.e4_real_8*pb*rs33
                cnb4  = 1._real_8+pc*rsu4+pd*rs24+1.e4_real_8*pb*rs34
                dcnb1 = -(pc*rsu1+2._real_8*pd*rs21+3.e4_real_8*pb*rs31)*&
                     0.3333333333333333_real_8
                dcnb2 = -(pc*rsu2+2._real_8*pd*rs22+3.e4_real_8*pb*rs32)*&
                     0.3333333333333333_real_8
                dcnb3 = -(pc*rsu3+2._real_8*pd*rs23+3.e4_real_8*pb*rs33)*&
                     0.3333333333333333_real_8
                dcnb4 = -(pc*rsu4+2._real_8*pd*rs24+3.e4_real_8*pb*rs34)*&
                     0.3333333333333333_real_8
                cn1   = cna1/cnb1 - cx
                cn2   = cna2/cnb2 - cx
                cn3   = cna3/cnb3 - cx
                cn4   = cna4/cnb4 - cx
                dcn1  = dcna1/cnb1 - cna1*dcnb1/(cnb1*cnb1)
                dcn2  = dcna2/cnb2 - cna2*dcnb2/(cnb2*cnb2)
                dcn3  = dcna3/cnb3 - cna3*dcnb3/(cnb3*cnb3)
                dcn4  = dcna4/cnb4 - cna4*dcnb4/(cnb4*cnb4)
                h11   = xnu*(cn1-cc0-3._real_8/7._real_8*cx)*t1*t1*EXP(ee1)
                h12   = xnu*(cn2-cc0-3._real_8/7._real_8*cx)*t2*t2*EXP(ee2)
                h13   = xnu*(cn3-cc0-3._real_8/7._real_8*cx)*t3*t3*EXP(ee3)
                h14   = xnu*(cn4-cc0-3._real_8/7._real_8*cx)*t4*t4*EXP(ee4)
                dh11  = -0.3333333333333_real_8*(h11*(7._real_8+8._real_8*ee1)+&
                     xnu*t1*t1*EXP(ee1)*dcn1)
                dh12  = -0.3333333333333_real_8*(h12*(7._real_8+8._real_8*ee2)+&
                     xnu*t2*t2*EXP(ee2)*dcn2)
                dh13  = -0.3333333333333_real_8*(h13*(7._real_8+8._real_8*ee3)+&
                     xnu*t3*t3*EXP(ee3)*dcn3)
                dh14  = -0.3333333333333_real_8*(h14*(7._real_8+8._real_8*ee4)+&
                     xnu*t4*t4*EXP(ee4)*dcn4)
                ddh11 = 2._real_8*h11*(1._real_8+ee1)*rhou1/grho1
                ddh12 = 2._real_8*h12*(1._real_8+ee2)*rhou2/grho2
                ddh13 = 2._real_8*h13*(1._real_8+ee3)*rhou3/grho3
                ddh14 = 2._real_8*h14*(1._real_8+ee4)*rhou4/grho4
                sc1   = rhou1*(h01+h11)
                sc2   = rhou2*(h02+h12)
                sc3   = rhou3*(h03+h13)
                sc4   = rhou4*(h04+h14)
                v1c1  = h01+h11+dh01+dh11
                v1c2  = h02+h12+dh02+dh12
                v1c3  = h03+h13+dh03+dh13
                v1c4  = h04+h14+dh04+dh14
                v2c1  = ddh01+ddh11
                v2c2  = ddh02+ddh12
                v2c3  = ddh03+ddh13
                v2c4  = ddh04+ddh14
                ! ..sum up terms
                sgcx=sgcx+(smoo1*sx1+smoo2*sx2+smoo3*sx3+smoo4*sx4)
                sgcc=sgcc+(smoo1*sc1+smoo2*sc2+smoo3*sc3+smoo4*sc4)
                v(i1,i2+0,i3+0)=dsmoo1*(sx1+sc1)+smoo1*(v1x1+v1c1)
                v(i1,i2+1,i3+0)=dsmoo2*(sx2+sc2)+smoo2*(v1x2+v1c2)
                v(i1,i2+0,i3+1)=dsmoo3*(sx3+sc3)+smoo3*(v1x3+v1c3)
                v(i1,i2+1,i3+1)=dsmoo4*(sx4+sc4)+smoo4*(v1x4+v1c4)
                vtmp(i1,i2+0,i3+0)=smoo1*(v2x1+v2c1)
                vtmp(i1,i2+1,i3+0)=smoo2*(v2x2+v2c2)
                vtmp(i1,i2+0,i3+1)=smoo3*(v2x3+v2c3)
                vtmp(i1,i2+1,i3+1)=smoo4*(v2x4+v2c4)
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    ! ..counting floating point operations
    fadd   = ic * 11._real_8 + 4._real_8*id *  62._real_8
    fmult  = ic * 12._real_8 + 4._real_8*id * 125._real_8
    fdiv   = ic *  4._real_8 + 4._real_8*id *  26._real_8
    fspez  =              4._real_8*id *   8._real_8
    flops=fadd+fmult+fdiv+fspez
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gcgga
  ! ==================================================================
  SUBROUTINE gcpbe(sgcx,sgcc,rhoe,v,vtmp,grad,flops)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: sgcx, sgcc
    REAL(real_8), DIMENSION(fpar%kr1, fpar%&
      kr2s, fpar%kr3s)                       :: rhoe
    COMPLEX(real_8) :: v(fpar%kr1,fpar%kr2s,fpar%kr3s)
    REAL(real_8), DIMENSION(fpar%kr1, fpar%&
      kr2s, fpar%kr3s)                       :: vtmp
    REAL(real_8), DIMENSION(fpar%kr1, fpar%&
      kr2s, fpar%kr3s, *)                    :: grad
    REAL(real_8)                             :: flops

    REAL(real_8), PARAMETER :: a0u = 0.4581652932831429_real_8, &
      a1u = 2.217058676663745_real_8, a2u = 0.7405551735357053_real_8, &
      a3u = 0.01968227878617998_real_8 , &
      al = 0.715996577859519256e-01_real_8, &
      ax = -0.738558766382022406_real_8, b1u = 1.0_real_8, &
      b2u = 4.504130959426697_real_8, b3u = 1.110667363742916_real_8, &
      b4u = 0.02359291751427506_real_8 , be = 0.06672455060314922_real_8 , &
      c1u = 4._real_8*a0u*b1u/3.0_real_8, &
      c2u = 5.0_real_8*a0u*b2u/3.0_real_8+a1u*b1u, c3u = 2.0_real_8*a0u*&
      b3u+4.0_real_8*a1u*b2u/3.0_real_8+2.0_real_8*a2u*b1u/3.0_real_8
    REAL(real_8), PARAMETER :: c4u = 7.0_real_8*a0u*b4u/3.0_real_8+5.0_real_8*&
      a1u*b3u/3.0_real_8+a2u*b2u+a3u*b1u/3.0_real_8, c5u = 2.0_real_8*a1u*&
      b4u+4.0_real_8*a2u*b3u/3.0_real_8+2.0_real_8*a3u*b2u/3.0_real_8, &
      c6u = 5.0_real_8*a2u*b4u/3.0_real_8+a3u*b3u, &
      c7u = 4.0_real_8*a3u*b4u/3.0_real_8 , eps2 = 1.e-20_real_8, &
      f1x = -1.10783814957303361_real_8 , rsfac = 0.6203504908994000_real_8 , &
      small = 1.e-24_real_8 , uk = 0.8040_real_8, &
      um = 0.2195149727645171_real_8, ul = um/uk , &
      us = 0.161620459673995492_real_8

    INTEGER                                  :: i1, i2, i3, ic, id, ii2, ii3, &
                                                ncount
    REAL(real_8) :: af1, af2, af3, af4, bf1, bf2, bf3, bf4, botu1, botu2, &
      botu3, botu4, d1, d2, d3, d4, ddh01, ddh02, ddh03, ddh04, dfx1, dfx2, &
      dfx3, dfx4, dh01, dh02, dh03, dh04, drho, dsmoo1, dsmoo2, dsmoo3, &
      dsmoo4, dtopu1, dtopu2, dtopu3, dtopu4, epsxcu1, epsxcu2, epsxcu3, &
      epsxcu4, ex1, ex2, ex3, ex4, expe1, expe2, expe3, expe4, fa1, fa2, &
      fadd, fdiv, fmult, fspez, fx1, fx2, fx3, fx4, grho1, grho2, grho3, &
      grho4, grhos1, grhos2, grhos3, grhos4, h01, h02, h03, h04, ogceps, po1, &
      po2, po3, po4, qy1, qy2, qy3, qy4, rho431, rho432, rho433, rho434, &
      rhou1, rhou2, rhou3, rhou4, rr1, rr2, rr3, rr4, rsu1
    REAL(real_8) :: rsu2, rsu3, rsu4, s1, s2, s21, s22, s23, s24, s3, s4, &
      sc1, sc2, sc3, sc4, smoo1, smoo2, smoo3, smoo4, sx1, sx2, sx3, sx4, t1, &
      t2, t3, t4, topu1, topu2, topu3, topu4, v1c1, v1c2, v1c3, v1c4, v1x1, &
      v1x2, v1x3, v1x4, v2c1, v2c2, v2c3, v2c4, v2x1, v2x2, v2x3, v2x4, vxc1, &
      vxc2, vxc3, vxc4, x1, x2, x3, x4, xkf1, xkf2, xkf3, xkf4, xkff, xks1, &
      xks2, xks3, xks4, xy1, xy2, xy3, xy4, y1, y2, y3, y4

! ==--------------------------------------------------------------==

    fa1=f1x*2._real_8/3._real_8
    fa2=fa1*4._real_8/3._real_8
    xkff  = (9._real_8*pi/4._real_8)**0.333333333333333_real_8
    sgcx = 0.0_real_8
    sgcc = 0.0_real_8
    ogceps=1._real_8/cntr%gceps
    x1=1._real_8
    x2=1._real_8
    x3=1._real_8
    x4=1._real_8
    ic=0
    id=0
    ii3=MOD(parm%nr3,2)
    ii2=MOD(parm%nr2,2)
    IF (ii3.NE.0.OR.ii2.NE.0) CALL stopgm('GCPBE','ODD DIMENSION',& 
         __LINE__,__FILE__)
    !$omp parallel do private(I1,I2,I3,RHOU1,RHOU2,RHOU3,RHOU4,D1,D2,D3,D4) &
    !$omp  firstprivate(X1,X2,X3,X4) &
    !$omp  private(NCOUNT,RHO431,RHO432,RHO433,RHO434,SMOO1,DSMOO1) &
    !$omp  private(DRHO,SMOO2,DSMOO2,SMOO3,DSMOO3,SMOO4,DSMOO4) &
    !$omp  private(GRHO1,GRHO2,GRHO3,GRHO4,GRHOS1,GRHOS2,GRHOS3) &
    !$omp  private(GRHOS4,RR1,RR2,RR3,RR4,EX1,EX2,EX3,EX4,S21) &
    !$omp  private(S22,S23,S24,PO1,PO2,PO3,PO4,FX1,FX2,FX3,FX4) &
    !$omp  private(SX1,SX2,SX3,SX4,DFX1,DFX2,DFX3,DFX4,V1X1) &
    !$omp  private(V1X2,V1X3,V1X4,V2X1,V2X2,V2X3,V2X4,RSU1,RSU2) &
    !$omp  private(RSU3,RSU4,TOPU1,TOPU2,TOPU3,TOPU4,DTOPU1,DTOPU2) &
    !$omp  private(DTOPU3,DTOPU4,BOTU1,BOTU2,BOTU3,BOTU4,EPSXCU1) &
    !$omp  private(EPSXCU2,EPSXCU3,EPSXCU4,VXC1,VXC2,VXC3,VXC4,XKF1) &
    !$omp  private(XKF2,XKF3,XKF4,XKS1,XKS2,XKS3,XKS4,T1,T2,T3,T4) &
    !$omp  private(EXPE1,EXPE2,EXPE3,EXPE4,AF1,AF2,AF3,AF4,BF1,BF2) &
    !$omp  private(BF3,BF4,Y1,Y2,Y3,Y4,XY1,XY2,XY3,XY4,QY1,QY2) &
    !$omp  private(QY3,QY4,S1,S2,S3,S4,H01,H02,H03,H04,DH01,DH02) &
    !$omp  private(DH03,DH04,DDH01,DDH02,DDH03,DDH04,SC1,SC2,SC3,SC4) &
    !$omp  private(V1C1,V1C2,V1C3,V1C4) &
    !$omp  private(V2C1,V2C2,V2C3,V2C4) &
    !$omp  reduction(+:SGCX,SGCC,IC,ID) __COLLAPSE3
    DO i3=1,parm%nr3,2
       DO i2=1,parm%nr2,2
          DO i1=1,parm%nr1
             rhou1=MAX(rhoe(i1,i2+0,i3+0),small)
             rhou2=MAX(rhoe(i1,i2+1,i3+0),small)
             rhou3=MAX(rhoe(i1,i2+0,i3+1),small)
             rhou4=MAX(rhoe(i1,i2+1,i3+1),small)
             IF ((rhou1+rhou2+rhou3+rhou4).GT.5._real_8*small) THEN
                id=id+1
                DO ncount=1,100
                   d1=x1-rhou1/x1**2
                   d2=x2-rhou2/x2**2
                   d3=x3-rhou3/x3**2
                   d4=x4-rhou4/x4**2
                   ic=ic+4
                   x1=x1-0.333333333333333_real_8*d1
                   x2=x2-0.333333333333333_real_8*d2
                   x3=x3-0.333333333333333_real_8*d3
                   x4=x4-0.333333333333333_real_8*d4
                   IF (d1**2+d2**2+d3**2+d4**2.LT.eps2) GOTO 10
                ENDDO
                x1=rhou1**(1._real_8/3._real_8)
                x2=rhou2**(1._real_8/3._real_8)
                x3=rhou3**(1._real_8/3._real_8)
                x4=rhou4**(1._real_8/3._real_8)
10              CONTINUE
                rho431=x1*rhou1
                rho432=x2*rhou2
                rho433=x3*rhou3
                rho434=x4*rhou4
                IF (rhou1.GT.2.25_real_8*cntr%gceps) THEN
                   smoo1=1._real_8
                   dsmoo1=0._real_8
                ELSEIF (rhou1.LT.0.25_real_8*cntr%gceps) THEN
                   smoo1=0._real_8
                   dsmoo1=0._real_8
                ELSE
                   drho=(rhou1-0.25_real_8*cntr%gceps)*ogceps
                   smoo1=drho*drho*(0.75_real_8-0.25_real_8*drho)
                   dsmoo1=1.5_real_8*drho*ogceps*(1._real_8-0.5_real_8*drho)
                ENDIF
                IF (rhou2.GT.2.25_real_8*cntr%gceps) THEN
                   smoo2=1._real_8
                   dsmoo2=0._real_8
                ELSEIF (rhou2.LT.0.25_real_8*cntr%gceps) THEN
                   smoo2=0._real_8
                   dsmoo2=0._real_8
                ELSE
                   drho=(rhou2-0.25_real_8*cntr%gceps)*ogceps
                   smoo2=drho*drho*(0.75_real_8-0.25_real_8*drho)
                   dsmoo2=1.5_real_8*drho*ogceps*(1._real_8-0.5_real_8*drho)
                ENDIF
                IF (rhou3.GT.2.25_real_8*cntr%gceps) THEN
                   smoo3=1._real_8
                   dsmoo3=0._real_8
                ELSEIF (rhou3.LT.0.25_real_8*cntr%gceps) THEN
                   smoo3=0._real_8
                   dsmoo3=0._real_8
                ELSE
                   drho=(rhou3-0.25_real_8*cntr%gceps)*ogceps
                   smoo3=drho*drho*(0.75_real_8-0.25_real_8*drho)
                   dsmoo3=1.5_real_8*drho*ogceps*(1._real_8-0.5_real_8*drho)
                ENDIF
                IF (rhou4.GT.2.25_real_8*cntr%gceps) THEN
                   smoo4=1._real_8
                   dsmoo4=0._real_8
                ELSEIF (rhou4.LT.0.25_real_8*cntr%gceps) THEN
                   smoo4=0._real_8
                   dsmoo4=0._real_8
                ELSE
                   drho=(rhou4-0.25_real_8*cntr%gceps)*ogceps
                   smoo4=drho*drho*(0.75_real_8-0.25_real_8*drho)
                   dsmoo4=1.5_real_8*drho*ogceps*(1._real_8-0.5_real_8*drho)
                ENDIF
                grho1=grad(i1,i2+0,i3+0,1)
                grho2=grad(i1,i2+1,i3+0,1)
                grho3=grad(i1,i2+0,i3+1,1)
                grho4=grad(i1,i2+1,i3+1,1)
                grhos1=SQRT(grho1)
                grhos2=SQRT(grho2)
                grhos3=SQRT(grho3)
                grhos4=SQRT(grho4)
                ! ..exchange part
                rr1   = 1._real_8/rho431
                rr2   = 1._real_8/rho432
                rr3   = 1._real_8/rho433
                rr4   = 1._real_8/rho434
                ex1   = ax/rr1
                ex2   = ax/rr2
                ex3   = ax/rr3
                ex4   = ax/rr4
                s21   = grho1*rr1*rr1*us*us
                s22   = grho2*rr2*rr2*us*us
                s23   = grho3*rr3*rr3*us*us
                s24   = grho4*rr4*rr4*us*us
                po1   = 1._real_8/(1._real_8 + ul*s21)
                po2   = 1._real_8/(1._real_8 + ul*s22)
                po3   = 1._real_8/(1._real_8 + ul*s23)
                po4   = 1._real_8/(1._real_8 + ul*s24)
                fx1   = uk-uk*po1
                fx2   = uk-uk*po2
                fx3   = uk-uk*po3
                fx4   = uk-uk*po4
                sx1   = ex1*fx1
                sx2   = ex2*fx2
                sx3   = ex3*fx3
                sx4   = ex4*fx4
                dfx1  = 2._real_8*uk*ul*po1*po1
                dfx2  = 2._real_8*uk*ul*po2*po2
                dfx3  = 2._real_8*uk*ul*po3*po3
                dfx4  = 2._real_8*uk*ul*po4*po4
                v1x1  = 1.33333333333333_real_8*ax*x1*(fx1-s21*dfx1)
                v1x2  = 1.33333333333333_real_8*ax*x2*(fx2-s22*dfx2)
                v1x3  = 1.33333333333333_real_8*ax*x3*(fx3-s23*dfx3)
                v1x4  = 1.33333333333333_real_8*ax*x4*(fx4-s24*dfx4)
                v2x1  = ex1*dfx1*(us*rr1)**2
                v2x2  = ex2*dfx2*(us*rr2)**2
                v2x3  = ex3*dfx3*(us*rr3)**2
                v2x4  = ex4*dfx4*(us*rr4)**2
                ! ..lda functional
                rsu1=rsfac/x1
                rsu2=rsfac/x2
                rsu3=rsfac/x3
                rsu4=rsfac/x4
                topu1=a2u+rsu1*a3u
                topu2=a2u+rsu2*a3u
                topu3=a2u+rsu3*a3u
                topu4=a2u+rsu4*a3u
                topu1=a1u+rsu1*topu1
                topu2=a1u+rsu2*topu2
                topu3=a1u+rsu3*topu3
                topu4=a1u+rsu4*topu4
                topu1=a0u+rsu1*topu1
                topu2=a0u+rsu2*topu2
                topu3=a0u+rsu3*topu3
                topu4=a0u+rsu4*topu4
                dtopu1=c6u+rsu1*c7u
                dtopu2=c6u+rsu2*c7u
                dtopu3=c6u+rsu3*c7u
                dtopu4=c6u+rsu4*c7u
                dtopu1=c5u+rsu1*dtopu1
                dtopu2=c5u+rsu2*dtopu2
                dtopu3=c5u+rsu3*dtopu3
                dtopu4=c5u+rsu4*dtopu4
                dtopu1=c4u+rsu1*dtopu1
                dtopu2=c4u+rsu2*dtopu2
                dtopu3=c4u+rsu3*dtopu3
                dtopu4=c4u+rsu4*dtopu4
                dtopu1=c3u+rsu1*dtopu1
                dtopu2=c3u+rsu2*dtopu2
                dtopu3=c3u+rsu3*dtopu3
                dtopu4=c3u+rsu4*dtopu4
                dtopu1=c2u+rsu1*dtopu1
                dtopu2=c2u+rsu2*dtopu2
                dtopu3=c2u+rsu3*dtopu3
                dtopu4=c2u+rsu4*dtopu4
                dtopu1=c1u+rsu1*dtopu1
                dtopu2=c1u+rsu2*dtopu2
                dtopu3=c1u+rsu3*dtopu3
                dtopu4=c1u+rsu4*dtopu4
                dtopu1=-rsu1*dtopu1
                dtopu2=-rsu2*dtopu2
                dtopu3=-rsu3*dtopu3
                dtopu4=-rsu4*dtopu4
                botu1=b3u+rsu1*b4u
                botu2=b3u+rsu2*b4u
                botu3=b3u+rsu3*b4u
                botu4=b3u+rsu4*b4u
                botu1=b2u+rsu1*botu1
                botu2=b2u+rsu2*botu2
                botu3=b2u+rsu3*botu3
                botu4=b2u+rsu4*botu4
                botu1=b1u+rsu1*botu1
                botu2=b1u+rsu2*botu2
                botu3=b1u+rsu3*botu3
                botu4=b1u+rsu4*botu4
                botu1=rsu1*botu1
                botu2=rsu2*botu2
                botu3=rsu3*botu3
                botu4=rsu4*botu4
                epsxcu1=topu1*(-1._real_8/botu1) - fa1*x1
                epsxcu2=topu2*(-1._real_8/botu2) - fa1*x2
                epsxcu3=topu3*(-1._real_8/botu3) - fa1*x3
                epsxcu4=topu4*(-1._real_8/botu4) - fa1*x4
                vxc1=dtopu1*(-1._real_8/botu1)*(-1._real_8/botu1) - fa2*x1
                vxc2=dtopu2*(-1._real_8/botu2)*(-1._real_8/botu2) - fa2*x2
                vxc3=dtopu3*(-1._real_8/botu3)*(-1._real_8/botu3) - fa2*x3
                vxc4=dtopu4*(-1._real_8/botu4)*(-1._real_8/botu4) - fa2*x4
                ! ..correlation part
                xkf1=xkff/rsu1
                xkf2=xkff/rsu2
                xkf3=xkff/rsu3
                xkf4=xkff/rsu4
                xks1=SQRT(4._real_8*xkf1/pi)
                xks2=SQRT(4._real_8*xkf2/pi)
                xks3=SQRT(4._real_8*xkf3/pi)
                xks4=SQRT(4._real_8*xkf4/pi)
                t1  =grhos1/(2._real_8*xks1*rhou1)
                t2  =grhos2/(2._real_8*xks2*rhou2)
                t3  =grhos3/(2._real_8*xks3*rhou3)
                t4  =grhos4/(2._real_8*xks4*rhou4)
                expe1=EXP(-2._real_8*al*epsxcu1/(be*be))
                expe2=EXP(-2._real_8*al*epsxcu2/(be*be))
                expe3=EXP(-2._real_8*al*epsxcu3/(be*be))
                expe4=EXP(-2._real_8*al*epsxcu4/(be*be))
                af1=2._real_8*al/be*(1._real_8/(expe1-1._real_8))
                af2=2._real_8*al/be*(1._real_8/(expe2-1._real_8))
                af3=2._real_8*al/be*(1._real_8/(expe3-1._real_8))
                af4=2._real_8*al/be*(1._real_8/(expe4-1._real_8))
                bf1=expe1*(vxc1-epsxcu1)
                bf2=expe2*(vxc2-epsxcu2)
                bf3=expe3*(vxc3-epsxcu3)
                bf4=expe4*(vxc4-epsxcu4)
                y1 = af1*t1*t1
                y2 = af2*t2*t2
                y3 = af3*t3*t3
                y4 = af4*t4*t4
                xy1 = (1._real_8+y1)/(1._real_8+y1+y1*y1)
                xy2 = (1._real_8+y2)/(1._real_8+y2+y2*y2)
                xy3 = (1._real_8+y3)/(1._real_8+y3+y3*y3)
                xy4 = (1._real_8+y4)/(1._real_8+y4+y4*y4)
                qy1 = y1*y1*(2._real_8+y1)/(1._real_8+y1+y1*y1)**2
                qy2 = y2*y2*(2._real_8+y2)/(1._real_8+y2+y2*y2)**2
                qy3 = y3*y3*(2._real_8+y3)/(1._real_8+y3+y3*y3)**2
                qy4 = y4*y4*(2._real_8+y4)/(1._real_8+y4+y4*y4)**2
                s1  = 1._real_8+2._real_8*al/be*t1*t1*xy1
                s2  = 1._real_8+2._real_8*al/be*t2*t2*xy2
                s3  = 1._real_8+2._real_8*al/be*t3*t3*xy3
                s4  = 1._real_8+2._real_8*al/be*t4*t4*xy4
                h01 = be*be/(2._real_8*al) * LOG(s1)
                h02 = be*be/(2._real_8*al) * LOG(s2)
                h03 = be*be/(2._real_8*al) * LOG(s3)
                h04 = be*be/(2._real_8*al) * LOG(s4)
                dh01= be*t1*t1/s1*(-2.33333333333333_real_8*xy1-qy1*&
                     (af1*bf1/be-2.3333333333333333_real_8))
                dh02= be*t2*t2/s2*(-2.33333333333333_real_8*xy2-qy2*&
                     (af2*bf2/be-2.3333333333333333_real_8))
                dh03= be*t3*t3/s3*(-2.33333333333333_real_8*xy3-qy3*&
                     (af3*bf3/be-2.3333333333333333_real_8))
                dh04= be*t4*t4/s4*(-2.33333333333333_real_8*xy4-qy4*&
                     (af4*bf4/be-2.3333333333333333_real_8))
                ddh01=be/(2._real_8*xks1*xks1*rhou1)*(xy1-qy1)/s1
                ddh02=be/(2._real_8*xks2*xks2*rhou2)*(xy2-qy2)/s2
                ddh03=be/(2._real_8*xks3*xks3*rhou3)*(xy3-qy3)/s3
                ddh04=be/(2._real_8*xks4*xks4*rhou4)*(xy4-qy4)/s4
                sc1   = rhou1*h01
                sc2   = rhou2*h02
                sc3   = rhou3*h03
                sc4   = rhou4*h04
                v1c1  = h01+dh01
                v1c2  = h02+dh02
                v1c3  = h03+dh03
                v1c4  = h04+dh04
                v2c1  = ddh01
                v2c2  = ddh02
                v2c3  = ddh03
                v2c4  = ddh04
                ! ..sum up terms
                sgcx=sgcx+(smoo1*sx1+smoo2*sx2+smoo3*sx3+smoo4*sx4)
                sgcc=sgcc+(smoo1*sc1+smoo2*sc2+smoo3*sc3+smoo4*sc4)
                v(i1,i2+0,i3+0)=dsmoo1*(sx1+sc1)+smoo1*(v1x1+v1c1)
                v(i1,i2+1,i3+0)=dsmoo2*(sx2+sc2)+smoo2*(v1x2+v1c2)
                v(i1,i2+0,i3+1)=dsmoo3*(sx3+sc3)+smoo3*(v1x3+v1c3)
                v(i1,i2+1,i3+1)=dsmoo4*(sx4+sc4)+smoo4*(v1x4+v1c4)
                vtmp(i1,i2+0,i3+0)=smoo1*(v2x1+v2c1)
                vtmp(i1,i2+1,i3+0)=smoo2*(v2x2+v2c2)
                vtmp(i1,i2+0,i3+1)=smoo3*(v2x3+v2c3)
                vtmp(i1,i2+1,i3+1)=smoo4*(v2x4+v2c4)
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    ! ..counting floating point operations
    fadd   = ic * 11._real_8 + 4._real_8*id *  23._real_8
    fmult  = ic * 12._real_8 + 4._real_8*id *  80._real_8
    fdiv   = ic *  4._real_8 + 4._real_8*id *  19._real_8
    fspez  =              4._real_8*id *   4._real_8
    flops=fadd+fmult+fdiv+fspez
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gcpbe
  ! ==================================================================
  SUBROUTINE gcrevpbe(sgcx,sgcc,rhoe,v,vtmp,grad,flops)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: sgcx, sgcc
    REAL(real_8), DIMENSION(fpar%kr1, fpar%&
      kr2s, fpar%kr3s)                       :: rhoe
    COMPLEX(real_8) :: v(fpar%kr1,fpar%kr2s,fpar%kr3s)
    REAL(real_8), DIMENSION(fpar%kr1, fpar%&
      kr2s, fpar%kr3s)                       :: vtmp
    REAL(real_8), DIMENSION(fpar%kr1, fpar%&
      kr2s, fpar%kr3s, *)                    :: grad
    REAL(real_8)                             :: flops

    REAL(real_8), PARAMETER :: a0u = 0.4581652932831429_real_8, &
      a1u = 2.217058676663745_real_8, a2u = 0.7405551735357053_real_8, &
      a3u = 0.01968227878617998_real_8 , &
      al = 0.715996577859519256e-01_real_8, &
      ax = -0.738558766382022406_real_8, b1u = 1.0_real_8, &
      b2u = 4.504130959426697_real_8, b3u = 1.110667363742916_real_8, &
      b4u = 0.02359291751427506_real_8 , be = 0.06672455060314922_real_8 , &
      c1u = 4._real_8*a0u*b1u/3.0_real_8, &
      c2u = 5.0_real_8*a0u*b2u/3.0_real_8+a1u*b1u, c3u = 2.0_real_8*a0u*&
      b3u+4.0_real_8*a1u*b2u/3.0_real_8+2.0_real_8*a2u*b1u/3.0_real_8
    REAL(real_8), PARAMETER :: c4u = 7.0_real_8*a0u*b4u/3.0_real_8+5.0_real_8*&
      a1u*b3u/3.0_real_8+a2u*b2u+a3u*b1u/3.0_real_8, c5u = 2.0_real_8*a1u*&
      b4u+4.0_real_8*a2u*b3u/3.0_real_8+2.0_real_8*a3u*b2u/3.0_real_8, &
      c6u = 5.0_real_8*a2u*b4u/3.0_real_8+a3u*b3u, &
      c7u = 4.0_real_8*a3u*b4u/3.0_real_8 , eps2 = 1.e-20_real_8, &
      f1x = -1.10783814957303361_real_8 , rsfac = 0.6203504908994000_real_8 , &
      small = 1.e-24_real_8 , uk = 1.2450_real_8, &
      um = 0.2195149727645171_real_8, ul = um/uk , &
      us = 0.161620459673995492_real_8

    INTEGER                                  :: i1, i2, i3, ic, id, ii2, ii3, &
                                                ncount
    REAL(real_8) :: af1, af2, af3, af4, bf1, bf2, bf3, bf4, botu1, botu2, &
      botu3, botu4, d1, d2, d3, d4, ddh01, ddh02, ddh03, ddh04, dfx1, dfx2, &
      dfx3, dfx4, dh01, dh02, dh03, dh04, drho, dsmoo1, dsmoo2, dsmoo3, &
      dsmoo4, dtopu1, dtopu2, dtopu3, dtopu4, epsxcu1, epsxcu2, epsxcu3, &
      epsxcu4, ex1, ex2, ex3, ex4, expe1, expe2, expe3, expe4, fa1, fa2, &
      fadd, fdiv, fmult, fspez, fx1, fx2, fx3, fx4, grho1, grho2, grho3, &
      grho4, grhos1, grhos2, grhos3, grhos4, h01, h02, h03, h04, ogceps, po1, &
      po2, po3, po4, qy1, qy2, qy3, qy4, rho431, rho432, rho433, rho434, &
      rhou1, rhou2, rhou3, rhou4, rr1, rr2, rr3, rr4, rsu1
    REAL(real_8) :: rsu2, rsu3, rsu4, s1, s2, s21, s22, s23, s24, s3, s4, &
      sc1, sc2, sc3, sc4, smoo1, smoo2, smoo3, smoo4, sx1, sx2, sx3, sx4, t1, &
      t2, t3, t4, topu1, topu2, topu3, topu4, v1c1, v1c2, v1c3, v1c4, v1x1, &
      v1x2, v1x3, v1x4, v2c1, v2c2, v2c3, v2c4, v2x1, v2x2, v2x3, v2x4, vxc1, &
      vxc2, vxc3, vxc4, x1, x2, x3, x4, xkf1, xkf2, xkf3, xkf4, xkff, xks1, &
      xks2, xks3, xks4, xy1, xy2, xy3, xy4, y1, y2, y3, y4

! ==--------------------------------------------------------------==

    fa1=f1x*2._real_8/3._real_8
    fa2=fa1*4._real_8/3._real_8
    xkff  = (9._real_8*pi/4._real_8)**0.333333333333333_real_8
    sgcx = 0.0_real_8
    sgcc = 0.0_real_8
    ogceps=1._real_8/cntr%gceps
    x1=1._real_8
    x2=1._real_8
    x3=1._real_8
    x4=1._real_8
    ic=0
    id=0
    ii3=MOD(parm%nr3,2)
    ii2=MOD(parm%nr2,2)
    IF (ii3.NE.0.OR.ii2.NE.0) CALL stopgm('GCrevPBE','ODD DIMENSION',& 
         __LINE__,__FILE__)
    !$omp parallel do private(I1,I2,I3,RHOU1,RHOU2,RHOU3,RHOU4,D1,D2,D3,D4) &
    !$omp  firstprivate(X1,X2,X3,X4) &
    !$omp  private(NCOUNT,RHO431,RHO432,RHO433,RHO434,SMOO1,DSMOO1) &
    !$omp  private(DRHO,SMOO2,DSMOO2,SMOO3,DSMOO3,SMOO4,DSMOO4) &
    !$omp  private(GRHO1,GRHO2,GRHO3,GRHO4,GRHOS1,GRHOS2,GRHOS3) &
    !$omp  private(GRHOS4,RR1,RR2,RR3,RR4,EX1,EX2,EX3,EX4,S21) &
    !$omp  private(S22,S23,S24,PO1,PO2,PO3,PO4,FX1,FX2,FX3) &
    !$omp  private(FX4,SX1,SX2,SX3,SX4,DFX1,DFX2,DFX3,DFX4,V1X1) &
    !$omp  private(V1X2,V1X3,V1X4,V2X1,V2X2,V2X3,V2X4,RSU1,RSU2) &
    !$omp  private(RSU3,RSU4,TOPU1,TOPU2,TOPU3,TOPU4,DTOPU1,DTOPU2) &
    !$omp  private(DTOPU3,DTOPU4,BOTU1,BOTU2,BOTU3,BOTU4,EPSXCU1) &
    !$omp  private(EPSXCU2,EPSXCU3,EPSXCU4,VXC1,VXC2,VXC3,VXC4,XKF1) &
    !$omp  private(XKF2,XKF3,XKF4,XKS1,XKS2,XKS3,XKS4,T1,T2,T3,T4) &
    !$omp  private(EXPE1,EXPE2,EXPE3,EXPE4,AF1,AF2,AF3,AF4,BF1) &
    !$omp  private(BF2,BF3,BF4,Y1,Y2,Y3,Y4,XY1,XY2,XY3,XY4,QY1) &
    !$omp  private(QY2,QY3,QY4,S1,S2,S3,S4,H01,H02,H03,H04,DH01) &
    !$omp  private(DH02,DH03,DH04,DDH01,DDH02,DDH03,DDH04,SC1,SC2) &
    !$omp  private(SC3,SC4,V1C1,V1C2,V1C3,V1C4) &
    !$omp  private(V2C1,V2C2,V2C3,V2C4) &
    !$omp  reduction(+:SGCX,SGCC,IC,ID) __COLLAPSE3
    DO i3=1,parm%nr3,2
       DO i2=1,parm%nr2,2
          DO i1=1,parm%nr1
             rhou1=MAX(rhoe(i1,i2+0,i3+0),small)
             rhou2=MAX(rhoe(i1,i2+1,i3+0),small)
             rhou3=MAX(rhoe(i1,i2+0,i3+1),small)
             rhou4=MAX(rhoe(i1,i2+1,i3+1),small)
             IF ((rhou1+rhou2+rhou3+rhou4).GT.5._real_8*small) THEN
                id=id+1
                DO ncount=1,100
                   d1=x1-rhou1/x1**2
                   d2=x2-rhou2/x2**2
                   d3=x3-rhou3/x3**2
                   d4=x4-rhou4/x4**2
                   ic=ic+4
                   x1=x1-0.333333333333333_real_8*d1
                   x2=x2-0.333333333333333_real_8*d2
                   x3=x3-0.333333333333333_real_8*d3
                   x4=x4-0.333333333333333_real_8*d4
                   IF (d1**2+d2**2+d3**2+d4**2.LT.eps2) GOTO 10
                ENDDO
                x1=rhou1**(1._real_8/3._real_8)
                x2=rhou2**(1._real_8/3._real_8)
                x3=rhou3**(1._real_8/3._real_8)
                x4=rhou4**(1._real_8/3._real_8)
10              CONTINUE
                rho431=x1*rhou1
                rho432=x2*rhou2
                rho433=x3*rhou3
                rho434=x4*rhou4
                IF (rhou1.GT.2.25_real_8*cntr%gceps) THEN
                   smoo1=1._real_8
                   dsmoo1=0._real_8
                ELSEIF (rhou1.LT.0.25_real_8*cntr%gceps) THEN
                   smoo1=0._real_8
                   dsmoo1=0._real_8
                ELSE
                   drho=(rhou1-0.25_real_8*cntr%gceps)*ogceps
                   smoo1=drho*drho*(0.75_real_8-0.25_real_8*drho)
                   dsmoo1=1.5_real_8*drho*ogceps*(1._real_8-0.5_real_8*drho)
                ENDIF
                IF (rhou2.GT.2.25_real_8*cntr%gceps) THEN
                   smoo2=1._real_8
                   dsmoo2=0._real_8
                ELSEIF (rhou2.LT.0.25_real_8*cntr%gceps) THEN
                   smoo2=0._real_8
                   dsmoo2=0._real_8
                ELSE
                   drho=(rhou2-0.25_real_8*cntr%gceps)*ogceps
                   smoo2=drho*drho*(0.75_real_8-0.25_real_8*drho)
                   dsmoo2=1.5_real_8*drho*ogceps*(1._real_8-0.5_real_8*drho)
                ENDIF
                IF (rhou3.GT.2.25_real_8*cntr%gceps) THEN
                   smoo3=1._real_8
                   dsmoo3=0._real_8
                ELSEIF (rhou3.LT.0.25_real_8*cntr%gceps) THEN
                   smoo3=0._real_8
                   dsmoo3=0._real_8
                ELSE
                   drho=(rhou3-0.25_real_8*cntr%gceps)*ogceps
                   smoo3=drho*drho*(0.75_real_8-0.25_real_8*drho)
                   dsmoo3=1.5_real_8*drho*ogceps*(1._real_8-0.5_real_8*drho)
                ENDIF
                IF (rhou4.GT.2.25_real_8*cntr%gceps) THEN
                   smoo4=1._real_8
                   dsmoo4=0._real_8
                ELSEIF (rhou4.LT.0.25_real_8*cntr%gceps) THEN
                   smoo4=0._real_8
                   dsmoo4=0._real_8
                ELSE
                   drho=(rhou4-0.25_real_8*cntr%gceps)*ogceps
                   smoo4=drho*drho*(0.75_real_8-0.25_real_8*drho)
                   dsmoo4=1.5_real_8*drho*ogceps*(1._real_8-0.5_real_8*drho)
                ENDIF
                grho1=grad(i1,i2+0,i3+0,1)
                grho2=grad(i1,i2+1,i3+0,1)
                grho3=grad(i1,i2+0,i3+1,1)
                grho4=grad(i1,i2+1,i3+1,1)
                grhos1=SQRT(grho1)
                grhos2=SQRT(grho2)
                grhos3=SQRT(grho3)
                grhos4=SQRT(grho4)
                ! ..exchange part
                rr1   = 1._real_8/rho431
                rr2   = 1._real_8/rho432
                rr3   = 1._real_8/rho433
                rr4   = 1._real_8/rho434
                ex1   = ax/rr1
                ex2   = ax/rr2
                ex3   = ax/rr3
                ex4   = ax/rr4
                s21   = grho1*rr1*rr1*us*us
                s22   = grho2*rr2*rr2*us*us
                s23   = grho3*rr3*rr3*us*us
                s24   = grho4*rr4*rr4*us*us
                po1   = 1._real_8/(1._real_8 + ul*s21)
                po2   = 1._real_8/(1._real_8 + ul*s22)
                po3   = 1._real_8/(1._real_8 + ul*s23)
                po4   = 1._real_8/(1._real_8 + ul*s24)
                fx1   = uk-uk*po1
                fx2   = uk-uk*po2
                fx3   = uk-uk*po3
                fx4   = uk-uk*po4
                sx1   = ex1*fx1
                sx2   = ex2*fx2
                sx3   = ex3*fx3
                sx4   = ex4*fx4
                dfx1  = 2._real_8*uk*ul*po1*po1
                dfx2  = 2._real_8*uk*ul*po2*po2
                dfx3  = 2._real_8*uk*ul*po3*po3
                dfx4  = 2._real_8*uk*ul*po4*po4
                v1x1  = 1.33333333333333_real_8*ax*x1*(fx1-s21*dfx1)
                v1x2  = 1.33333333333333_real_8*ax*x2*(fx2-s22*dfx2)
                v1x3  = 1.33333333333333_real_8*ax*x3*(fx3-s23*dfx3)
                v1x4  = 1.33333333333333_real_8*ax*x4*(fx4-s24*dfx4)
                v2x1  = ex1*dfx1*(us*rr1)**2
                v2x2  = ex2*dfx2*(us*rr2)**2
                v2x3  = ex3*dfx3*(us*rr3)**2
                v2x4  = ex4*dfx4*(us*rr4)**2
                ! ..lda functional
                rsu1=rsfac/x1
                rsu2=rsfac/x2
                rsu3=rsfac/x3
                rsu4=rsfac/x4
                topu1=a2u+rsu1*a3u
                topu2=a2u+rsu2*a3u
                topu3=a2u+rsu3*a3u
                topu4=a2u+rsu4*a3u
                topu1=a1u+rsu1*topu1
                topu2=a1u+rsu2*topu2
                topu3=a1u+rsu3*topu3
                topu4=a1u+rsu4*topu4
                topu1=a0u+rsu1*topu1
                topu2=a0u+rsu2*topu2
                topu3=a0u+rsu3*topu3
                topu4=a0u+rsu4*topu4
                dtopu1=c6u+rsu1*c7u
                dtopu2=c6u+rsu2*c7u
                dtopu3=c6u+rsu3*c7u
                dtopu4=c6u+rsu4*c7u
                dtopu1=c5u+rsu1*dtopu1
                dtopu2=c5u+rsu2*dtopu2
                dtopu3=c5u+rsu3*dtopu3
                dtopu4=c5u+rsu4*dtopu4
                dtopu1=c4u+rsu1*dtopu1
                dtopu2=c4u+rsu2*dtopu2
                dtopu3=c4u+rsu3*dtopu3
                dtopu4=c4u+rsu4*dtopu4
                dtopu1=c3u+rsu1*dtopu1
                dtopu2=c3u+rsu2*dtopu2
                dtopu3=c3u+rsu3*dtopu3
                dtopu4=c3u+rsu4*dtopu4
                dtopu1=c2u+rsu1*dtopu1
                dtopu2=c2u+rsu2*dtopu2
                dtopu3=c2u+rsu3*dtopu3
                dtopu4=c2u+rsu4*dtopu4
                dtopu1=c1u+rsu1*dtopu1
                dtopu2=c1u+rsu2*dtopu2
                dtopu3=c1u+rsu3*dtopu3
                dtopu4=c1u+rsu4*dtopu4
                dtopu1=-rsu1*dtopu1
                dtopu2=-rsu2*dtopu2
                dtopu3=-rsu3*dtopu3
                dtopu4=-rsu4*dtopu4
                botu1=b3u+rsu1*b4u
                botu2=b3u+rsu2*b4u
                botu3=b3u+rsu3*b4u
                botu4=b3u+rsu4*b4u
                botu1=b2u+rsu1*botu1
                botu2=b2u+rsu2*botu2
                botu3=b2u+rsu3*botu3
                botu4=b2u+rsu4*botu4
                botu1=b1u+rsu1*botu1
                botu2=b1u+rsu2*botu2
                botu3=b1u+rsu3*botu3
                botu4=b1u+rsu4*botu4
                botu1=rsu1*botu1
                botu2=rsu2*botu2
                botu3=rsu3*botu3
                botu4=rsu4*botu4
                epsxcu1=topu1*(-1._real_8/botu1) - fa1*x1
                epsxcu2=topu2*(-1._real_8/botu2) - fa1*x2
                epsxcu3=topu3*(-1._real_8/botu3) - fa1*x3
                epsxcu4=topu4*(-1._real_8/botu4) - fa1*x4
                vxc1=dtopu1*(-1._real_8/botu1)*(-1._real_8/botu1) - fa2*x1
                vxc2=dtopu2*(-1._real_8/botu2)*(-1._real_8/botu2) - fa2*x2
                vxc3=dtopu3*(-1._real_8/botu3)*(-1._real_8/botu3) - fa2*x3
                vxc4=dtopu4*(-1._real_8/botu4)*(-1._real_8/botu4) - fa2*x4
                ! ..correlation part
                xkf1=xkff/rsu1
                xkf2=xkff/rsu2
                xkf3=xkff/rsu3
                xkf4=xkff/rsu4
                xks1=SQRT(4._real_8*xkf1/pi)
                xks2=SQRT(4._real_8*xkf2/pi)
                xks3=SQRT(4._real_8*xkf3/pi)
                xks4=SQRT(4._real_8*xkf4/pi)
                t1  =grhos1/(2._real_8*xks1*rhou1)
                t2  =grhos2/(2._real_8*xks2*rhou2)
                t3  =grhos3/(2._real_8*xks3*rhou3)
                t4  =grhos4/(2._real_8*xks4*rhou4)
                expe1=EXP(-2._real_8*al*epsxcu1/(be*be))
                expe2=EXP(-2._real_8*al*epsxcu2/(be*be))
                expe3=EXP(-2._real_8*al*epsxcu3/(be*be))
                expe4=EXP(-2._real_8*al*epsxcu4/(be*be))
                af1=2._real_8*al/be*(1._real_8/(expe1-1._real_8))
                af2=2._real_8*al/be*(1._real_8/(expe2-1._real_8))
                af3=2._real_8*al/be*(1._real_8/(expe3-1._real_8))
                af4=2._real_8*al/be*(1._real_8/(expe4-1._real_8))
                bf1=expe1*(vxc1-epsxcu1)
                bf2=expe2*(vxc2-epsxcu2)
                bf3=expe3*(vxc3-epsxcu3)
                bf4=expe4*(vxc4-epsxcu4)
                y1 = af1*t1*t1
                y2 = af2*t2*t2
                y3 = af3*t3*t3
                y4 = af4*t4*t4
                xy1 = (1._real_8+y1)/(1._real_8+y1+y1*y1)
                xy2 = (1._real_8+y2)/(1._real_8+y2+y2*y2)
                xy3 = (1._real_8+y3)/(1._real_8+y3+y3*y3)
                xy4 = (1._real_8+y4)/(1._real_8+y4+y4*y4)
                qy1 = y1*y1*(2._real_8+y1)/(1._real_8+y1+y1*y1)**2
                qy2 = y2*y2*(2._real_8+y2)/(1._real_8+y2+y2*y2)**2
                qy3 = y3*y3*(2._real_8+y3)/(1._real_8+y3+y3*y3)**2
                qy4 = y4*y4*(2._real_8+y4)/(1._real_8+y4+y4*y4)**2
                s1  = 1._real_8+2._real_8*al/be*t1*t1*xy1
                s2  = 1._real_8+2._real_8*al/be*t2*t2*xy2
                s3  = 1._real_8+2._real_8*al/be*t3*t3*xy3
                s4  = 1._real_8+2._real_8*al/be*t4*t4*xy4
                h01 = be*be/(2._real_8*al) * LOG(s1)
                h02 = be*be/(2._real_8*al) * LOG(s2)
                h03 = be*be/(2._real_8*al) * LOG(s3)
                h04 = be*be/(2._real_8*al) * LOG(s4)
                dh01= be*t1*t1/s1*(-2.33333333333333_real_8*xy1-qy1*&
                     (af1*bf1/be-2.3333333333333333_real_8))
                dh02= be*t2*t2/s2*(-2.33333333333333_real_8*xy2-qy2*&
                     (af2*bf2/be-2.3333333333333333_real_8))
                dh03= be*t3*t3/s3*(-2.33333333333333_real_8*xy3-qy3*&
                     (af3*bf3/be-2.3333333333333333_real_8))
                dh04= be*t4*t4/s4*(-2.33333333333333_real_8*xy4-qy4*&
                     (af4*bf4/be-2.3333333333333333_real_8))
                ddh01=be/(2._real_8*xks1*xks1*rhou1)*(xy1-qy1)/s1
                ddh02=be/(2._real_8*xks2*xks2*rhou2)*(xy2-qy2)/s2
                ddh03=be/(2._real_8*xks3*xks3*rhou3)*(xy3-qy3)/s3
                ddh04=be/(2._real_8*xks4*xks4*rhou4)*(xy4-qy4)/s4
                sc1   = rhou1*h01
                sc2   = rhou2*h02
                sc3   = rhou3*h03
                sc4   = rhou4*h04
                v1c1  = h01+dh01
                v1c2  = h02+dh02
                v1c3  = h03+dh03
                v1c4  = h04+dh04
                v2c1  = ddh01
                v2c2  = ddh02
                v2c3  = ddh03
                v2c4  = ddh04
                ! ..sum up terms
                sgcx=sgcx+(smoo1*sx1+smoo2*sx2+smoo3*sx3+smoo4*sx4)
                sgcc=sgcc+(smoo1*sc1+smoo2*sc2+smoo3*sc3+smoo4*sc4)
                v(i1,i2+0,i3+0)=dsmoo1*(sx1+sc1)+smoo1*(v1x1+v1c1)
                v(i1,i2+1,i3+0)=dsmoo2*(sx2+sc2)+smoo2*(v1x2+v1c2)
                v(i1,i2+0,i3+1)=dsmoo3*(sx3+sc3)+smoo3*(v1x3+v1c3)
                v(i1,i2+1,i3+1)=dsmoo4*(sx4+sc4)+smoo4*(v1x4+v1c4)
                vtmp(i1,i2+0,i3+0)=smoo1*(v2x1+v2c1)
                vtmp(i1,i2+1,i3+0)=smoo2*(v2x2+v2c2)
                vtmp(i1,i2+0,i3+1)=smoo3*(v2x3+v2c3)
                vtmp(i1,i2+1,i3+1)=smoo4*(v2x4+v2c4)
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    ! ..counting floating point operations
    fadd   = ic * 11._real_8 + 4._real_8*id *  23._real_8
    fmult  = ic * 12._real_8 + 4._real_8*id *  80._real_8
    fdiv   = ic *  4._real_8 + 4._real_8*id *  19._real_8
    fspez  =              4._real_8*id *   4._real_8
    flops=fadd+fmult+fdiv+fspez
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gcrevpbe
  ! ==================================================================
  SUBROUTINE gcsxonly(b1,sgcx,sgcc,rhoe,v1,v2,vtmp,grad,flops)
    REAL(real_8)                             :: b1, sgcx, sgcc
    REAL(real_8), DIMENSION(fpar%kr1, fpar%&
      kr2s, fpar%kr3s, 2)                    :: rhoe
    COMPLEX(real_8) :: v1(fpar%kr1,fpar%kr2s,fpar%kr3s), &
      v2(fpar%kr1,fpar%kr2s,fpar%kr3s)
    REAL(real_8), DIMENSION(fpar%kr1, fpar%&
      kr2s, fpar%kr3s, 3)                    :: vtmp
    REAL(real_8), DIMENSION(fpar%kr1, fpar%&
      kr2s, fpar%kr3s, 8)                    :: grad
    REAL(real_8)                             :: flops

    REAL(real_8), PARAMETER                  :: eps2 = 1.e-20_real_8, &
                                                small = 1.e-24_real_8 

    INTEGER                                  :: i1, i2, i3, ic, id, ii2, &
                                                ncount
    REAL(real_8) :: d1, d2, dda1, dda2, ddb1, ddb2, ddda1, ddda2, dddb1, &
      dddb2, dgfa1, dgfa2, dgfb1, dgfb2, drho, dsmoo1, dsmoo2, fadd, fdiv, &
      fmult, fspez, gfa1, gfa2, gfb1, gfb2, grhoa1, grhoa2, grhob1, grhob2, &
      ogceps, rho1, rho2, rhoa1, rhoa2, rhob1, rhob2, rsa1, rsa2, rsb1, rsb2, &
      sa2a1, sa2a2, sa2b1, sa2b2, shma1, shma2, shmb1, shmb2, smoo1, smoo2, &
      sxa1, sxa2, sxb1, sxb2, v1xa1, v1xa2, v1xb1, v1xb2, v2xa1, v2xa2, &
      v2xb1, v2xb2, x1, x2, xsa1, xsa2, xsaa1, xsaa2, xsb1, xsb2, xsbb1, &
      xsbb2, y1, y2

! ==--------------------------------------------------------------==

    ogceps=1._real_8/cntr%gceps
    sgcx=0._real_8
    sgcc=0._real_8
    x1=1._real_8
    x2=1._real_8
    y1=1._real_8
    y2=1._real_8
    ic=0
    id=0
    ii2=MOD(parm%nr2,2)
    IF (ii2.NE.0) CALL stopgm('GCSXONLY','ODD DIMENSIONS',& 
         __LINE__,__FILE__)
    !$omp parallel do private(I1,I2,I3,RHOA1,RHOA2,D1,D2) &
    !$omp  firstprivate(X1,X2,Y1,Y2) &
    !$omp  private(NCOUNT,RHOB1,RHOB2,RHO1,RHO2,RSA1,RSA2) &
    !$omp  private(RSB1,RSB2,SMOO1,DSMOO1,DRHO,SMOO2,DSMOO2) &
    !$omp  private(GRHOA1,GRHOA2,GRHOB1,GRHOB2,XSA1,XSA2,XSB1,XSB2) &
    !$omp  private(XSAA1,XSAA2,XSBB1,XSBB2,SA2A1,SA2A2,SA2B1,SA2B2) &
    !$omp  private(SHMA1,SHMA2,SHMB1,SHMB2,DDA1,DDA2,DDB1,DDB2,DDDA1) &
    !$omp  private(DDDA2,DDDB1,DDDB2,GFA1,GFA2,GFB1,GFB2,DGFA1,DGFA2) &
    !$omp  private(DGFB1,DGFB2,SXA1,SXA2,SXB1,SXB2,V1XA1,V1XA2,V1XB1) &
    !$omp  private(V1XB2,V2XA1,V2XA2,V2XB1,V2XB2) &
    !$omp  reduction(+:SGCX,IC,ID) __COLLAPSE3
    DO i3=1,parm%nr3
       DO i2=1,parm%nr2,2
          DO i1=1,parm%nr1
             rhoa1=MAX(rhoe(i1,i2+0,i3,1),small)
             rhoa2=MAX(rhoe(i1,i2+1,i3,1),small)
             rhob1=MAX(rhoe(i1,i2+0,i3,2),small)
             rhob2=MAX(rhoe(i1,i2+1,i3,2),small)
             rho1=rhoa1+rhob1
             rho2=rhoa2+rhob2
             IF ((rho1+rho2).GT.5._real_8*small) THEN
                id=id+1
                DO ncount=1,100
                   d1=x1-rhoa1/x1**2
                   d2=x2-rhoa2/x2**2
                   ic=ic+1
                   x1=x1-0.333333333333333_real_8*d1
                   x2=x2-0.333333333333333_real_8*d2
                   IF (d1**2+d2**2.LT.eps2) GOTO 10
                ENDDO
                x1=rhoa1**(1._real_8/3._real_8)
                x2=rhoa2**(1._real_8/3._real_8)
10              CONTINUE
                rsa1=1._real_8/x1
                rsa2=1._real_8/x2
                DO ncount=1,100
                   d1=y1-rhob1/y1**2
                   d2=y2-rhob2/y2**2
                   ic=ic+1
                   y1=y1-0.333333333333333_real_8*d1
                   y2=y2-0.333333333333333_real_8*d2
                   IF (d1**2+d2**2.LT.eps2) GOTO 20
                ENDDO
                y1=rhob1**(1._real_8/3._real_8)
                y2=rhob2**(1._real_8/3._real_8)
20              CONTINUE
                rsb1=1._real_8/y1
                rsb2=1._real_8/y2
                IF (rho1.GT.2.25_real_8*cntr%gceps) THEN
                   smoo1=1._real_8
                   dsmoo1=0._real_8
                ELSEIF (rho1.LT.0.25_real_8*cntr%gceps) THEN
                   smoo1=0._real_8
                   dsmoo1=0._real_8
                ELSE
                   drho=(rho1-0.25_real_8*cntr%gceps)*ogceps
                   smoo1=drho*drho*(0.75_real_8-0.25_real_8*drho)
                   dsmoo1=1.5_real_8*drho*ogceps*(1._real_8-0.5_real_8*drho)
                ENDIF
                IF (rho2.GT.2.25_real_8*cntr%gceps) THEN
                   smoo2=1._real_8
                   dsmoo2=0._real_8
                ELSEIF (rho2.LT.0.25_real_8*cntr%gceps) THEN
                   smoo2=0._real_8
                   dsmoo2=0._real_8
                ELSE
                   drho=(rho2-0.25_real_8*cntr%gceps)*ogceps
                   smoo2=drho*drho*(0.75_real_8-0.25_real_8*drho)
                   dsmoo2=1.5_real_8*drho*ogceps*(1._real_8-0.5_real_8*drho)
                ENDIF
                grhoa1=SQRT(grad(i1,i2+0,i3,1))
                grhoa2=SQRT(grad(i1,i2+1,i3,1))
                grhob1=SQRT(grad(i1,i2+0,i3,5))
                grhob2=SQRT(grad(i1,i2+1,i3,5))
                grhoa1=MAX(grhoa1,small)
                grhoa2=MAX(grhoa2,small)
                grhob1=MAX(grhob1,small)
                grhob2=MAX(grhob2,small)
                xsa1  = grhoa1*(rsa1*rsa1*rsa1*rsa1)
                xsa2  = grhoa2*(rsa2*rsa2*rsa2*rsa2)
                xsb1  = grhob1*(rsb1*rsb1*rsb1*rsb1)
                xsb2  = grhob2*(rsb2*rsb2*rsb2*rsb2)
                IF (rhoa1.LT.2._real_8*small) xsa1=0._real_8
                IF (rhoa2.LT.2._real_8*small) xsa2=0._real_8
                IF (rhob1.LT.2._real_8*small) xsb1=0._real_8
                IF (rhob2.LT.2._real_8*small) xsb2=0._real_8
                xsaa1 = xsa1*xsa1
                xsaa2 = xsa2*xsa2
                xsbb1 = xsb1*xsb1
                xsbb2 = xsb2*xsb2
                sa2a1 = SQRT(1._real_8+xsaa1)
                sa2a2 = SQRT(1._real_8+xsaa2)
                sa2b1 = SQRT(1._real_8+xsbb1)
                sa2b2 = SQRT(1._real_8+xsbb2)
                shma1 = LOG(xsa1+sa2a1)
                shma2 = LOG(xsa2+sa2a2)
                shmb1 = LOG(xsb1+sa2b1)
                shmb2 = LOG(xsb2+sa2b2)
                dda1  = 1._real_8+6._real_8*b1*xsa1*shma1
                dda2  = 1._real_8+6._real_8*b1*xsa2*shma2
                ddb1  = 1._real_8+6._real_8*b1*xsb1*shmb1
                ddb2  = 1._real_8+6._real_8*b1*xsb2*shmb2
                ddda1 = 6._real_8*b1*(shma1+xsa1/sa2a1)
                ddda2 = 6._real_8*b1*(shma2+xsa2/sa2a2)
                dddb1 = 6._real_8*b1*(shmb1+xsb1/sa2b1)
                dddb2 = 6._real_8*b1*(shmb2+xsb2/sa2b2)
                gfa1  = -b1*xsaa1/dda1
                gfa2  = -b1*xsaa2/dda2
                gfb1  = -b1*xsbb1/ddb1
                gfb2  = -b1*xsbb2/ddb2
                dgfa1 = (-2._real_8*b1*xsa1*dda1+b1*xsaa1*ddda1)/(dda1*dda1)
                dgfa2 = (-2._real_8*b1*xsa2*dda2+b1*xsaa2*ddda2)/(dda2*dda2)
                dgfb1 = (-2._real_8*b1*xsb1*ddb1+b1*xsbb1*dddb1)/(ddb1*ddb1)
                dgfb2 = (-2._real_8*b1*xsb2*ddb2+b1*xsbb2*dddb2)/(ddb2*ddb2)
                sxa1  = gfa1*x1*rhoa1
                sxa2  = gfa2*x2*rhoa2
                sxb1  = gfb1*y1*rhob1
                sxb2  = gfb2*y2*rhob2
                v1xa1 = 1.33333333333333_real_8*x1*(gfa1-xsa1*dgfa1)
                v1xa2 = 1.33333333333333_real_8*x2*(gfa2-xsa2*dgfa2)
                v1xb1 = 1.33333333333333_real_8*y1*(gfb1-xsb1*dgfb1)
                v1xb2 = 1.33333333333333_real_8*y2*(gfb2-xsb2*dgfb2)
                v2xa1 = dgfa1/grhoa1
                v2xa2 = dgfa2/grhoa2
                v2xb1 = dgfb1/grhob1
                v2xb2 = dgfb2/grhob2
                ! ..sum up terms
                sgcx = sgcx + (smoo1*(sxa1+sxb1) + smoo2*(sxa2+sxb2))
                v1(i1,i2+0,i3)=dsmoo1*sxa1+smoo1*v1xa1
                v1(i1,i2+1,i3)=dsmoo2*sxa2+smoo2*v1xa2
                v2(i1,i2+0,i3)=dsmoo1*sxb1+smoo1*v1xb1
                v2(i1,i2+1,i3)=dsmoo2*sxb2+smoo2*v1xb2
                vtmp(i1,i2+0,i3,1)=smoo1*v2xa1
                vtmp(i1,i2+1,i3,1)=smoo2*v2xa2
                vtmp(i1,i2+0,i3,2)=smoo1*v2xb1
                vtmp(i1,i2+1,i3,2)=smoo2*v2xb2
                vtmp(i1,i2+0,i3,3)=0._real_8
                vtmp(i1,i2+1,i3,3)=0._real_8
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    ! ..counting floating point operations
    fadd   = ic * 11._real_8 + id *  32._real_8
    fmult  = ic * 12._real_8 + id * 102._real_8
    fdiv   = ic *  4._real_8 + id *  16._real_8
    fspez  =              id *  12._real_8
    flops=fadd+fmult+fdiv+fspez
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gcsxonly
  ! ==================================================================
  SUBROUTINE gcsxp86(b1,sgcx,sgcc,rhoe,v1,v2,vtmp,grad,flops)
    REAL(real_8)                             :: b1, sgcx, sgcc
    REAL(real_8), DIMENSION(fpar%kr1, fpar%&
      kr2s, fpar%kr3s, 2)                    :: rhoe
    COMPLEX(real_8) :: v1(fpar%kr1,fpar%kr2s,fpar%kr3s), &
      v2(fpar%kr1,fpar%kr2s,fpar%kr3s)
    REAL(real_8), DIMENSION(fpar%kr1, fpar%&
      kr2s, fpar%kr3s, 3)                    :: vtmp
    REAL(real_8), DIMENSION(fpar%kr1, fpar%&
      kr2s, fpar%kr3s, 8)                    :: grad
    REAL(real_8)                             :: flops

    REAL(real_8), PARAMETER :: eps2 = 1.e-20_real_8, p1 = 0.023266_real_8, &
      p2 = 7.389e-6_real_8, p3 = 8.723_real_8, p4 = 0.472_real_8 , &
      pc1 = 0.001667_real_8, pc2 = 0.002568_real_8, pci = pc1+pc2 , &
      rsfac = 0.6203504908994000_real_8 , small = 1.e-24_real_8 

    INTEGER                                  :: i1, i2, i3, ic, id, ii2, &
                                                ncount
    REAL(real_8) :: cn1, cn2, cna1, cna2, cnb1, cnb2, d1, d2, dcn1, dcn2, &
      dcna1, dcna2, dcnb1, dcnb2, dda1, dda2, ddb1, ddb2, ddda1, ddda2, &
      dddb1, dddb2, dfx, dgfa1, dgfa2, dgfb1, dgfb2, drho, drs1, drs2, &
      dsmoo1, dsmoo2, ephi1, ephi2, fadd, fdiv, fmult, fspez, gfa1, gfa2, &
      gfb1, gfb2, grho1, grho2, grhoa1, grhoa2, grhob1, grhob2, ogceps, phi1, &
      phi2, rho1, rho2, rhoa1, rhoa2, rhob1, rhob2, rs1, rs2, rs21, rs22, &
      rs31, rs32, rsa1, rsa2, rsb1, rsb2, s11, s12, s21, s22, sa2a1, sa2a2, &
      sa2b1, sa2b2, sc1, sc2, shma1, shma2, shmb1, shmb2, smoo1, smoo2, sxa1, &
      sxa2, sxb1, sxb2, two13, v1c1, v1c2, v1ca1
    REAL(real_8) :: v1ca2, v1cb1, v1cb2, v1xa1, v1xa2, v1xb1, v1xb2, v2c1, &
      v2c2, v2ca1, v2ca2, v2cab1, v2cab2, v2cb1, v2cb2, v2xa1, v2xa2, v2xb1, &
      v2xb2, x1, x2, xsa1, xsa2, xsaa1, xsaa2, xsb1, xsb2, xsbb1, xsbb2, y1, &
      y2, z1, z2, zz1, zz2

! ==--------------------------------------------------------------==

    two13=2._real_8**0.333333333333333_real_8
    ogceps=1._real_8/cntr%gceps
    dfx = -0.3333333333333_real_8*(3._real_8/fpi)**0.3333333333333_real_8
    sgcx=0._real_8
    sgcc=0._real_8
    x1=1._real_8
    x2=1._real_8
    y1=1._real_8
    y2=1._real_8
    z1=1._real_8
    z2=1._real_8
    ic=0
    id=0
    ii2=MOD(parm%nr2,2)
    IF (ii2.NE.0) CALL stopgm('GCSXP86','ODD DIMENSIONS',& 
         __LINE__,__FILE__)
    !$omp parallel do private(I1,I2,I3,RHOA1,RHOA2,D1,D2) &
    !$omp  firstprivate(X1,X2,Y1,Y2,Z1,Z2) &
    !$omp  private(NCOUNT, RHOB1, RHOB2, RHO1, RHO2, RSA1, RSA2) &
    !$omp  private(RSB1, RSB2, ZZ1, ZZ2, RS1, RS2, SMOO1) &
    !$omp  private(DSMOO1, DRHO, SMOO2, DSMOO2, GRHOA1, GRHOA2, GRHOB1) &
    !$omp  private(GRHOB2, XSA1, XSA2, XSB1, XSB2, XSAA1, XSAA2, XSBB1) &
    !$omp  private(XSBB2, SA2A1, SA2A2, SA2B1, SA2B2, SHMA1, SHMA2, SHMB1) &
    !$omp  private(SHMB2, DDA1, DDA2, DDB1, DDB2, DDDA1, DDDA2, DDDB1) &
    !$omp  private(DDDB2, GFA1, GFA2, GFB1, GFB2, DGFA1, DGFA2, DGFB1) &
    !$omp  private(DGFB2, SXA1, SXA2, SXB1, SXB2, V1XA1, V1XA2, V1XB1) &
    !$omp  private(V1XB2, V2XA1, V2XA2, V2XB1, V2XB2, GRHO1, GRHO2, RS21) &
    !$omp  private(RS22, RS31, RS32, CNA1, CNA2, CNB1, CNB2, CN1, CN2) &
    !$omp  private(DRS1, DRS2, DCNA1, DCNA2, DCNB1, DCNB2, DCN1, DCN2) &
    !$omp  private(S11, S12, S21, S22, PHI1, PHI2, EPHI1, EPHI2, SC1, SC2) &
    !$omp  private(V1C1, V1C2, V1CA1, V1CA2, V1CB1, V1CB2, V2C1, V2C2) &
    !$omp  private(V2CA1,V2CA2,V2CB1,V2CB2,V2CAB1,V2CAB2) &
    !$omp  reduction(+:SGCX,SGCC,IC,ID) __COLLAPSE3
    DO i3=1,parm%nr3
       DO i2=1,parm%nr2,2
          DO i1=1,parm%nr1
             rhoa1=MAX(rhoe(i1,i2+0,i3,1),small)
             rhoa2=MAX(rhoe(i1,i2+1,i3,1),small)
             rhob1=MAX(rhoe(i1,i2+0,i3,2),small)
             rhob2=MAX(rhoe(i1,i2+1,i3,2),small)
             rho1=rhoa1+rhob1
             rho2=rhoa2+rhob2
             IF ((rho1+rho2).GT.5._real_8*small) THEN
                id=id+1
                DO ncount=1,100
                   d1=x1-rhoa1/x1**2
                   d2=x2-rhoa2/x2**2
                   ic=ic+1
                   x1=x1-0.333333333333333_real_8*d1
                   x2=x2-0.333333333333333_real_8*d2
                   IF (d1**2+d2**2.LT.eps2) GOTO 10
                ENDDO
                x1=rhoa1**(1._real_8/3._real_8)
                x2=rhoa2**(1._real_8/3._real_8)
10              CONTINUE
                rsa1=1._real_8/x1
                rsa2=1._real_8/x2
                DO ncount=1,100
                   d1=y1-rhob1/y1**2
                   d2=y2-rhob2/y2**2
                   ic=ic+1
                   y1=y1-0.333333333333333_real_8*d1
                   y2=y2-0.333333333333333_real_8*d2
                   IF (d1**2+d2**2.LT.eps2) GOTO 20
                ENDDO
                y1=rhob1**(1._real_8/3._real_8)
                y2=rhob2**(1._real_8/3._real_8)
20              CONTINUE
                rsb1=1._real_8/y1
                rsb2=1._real_8/y2
                DO ncount=1,100
                   d1=z1-rho1/z1**2
                   d2=z2-rho2/z2**2
                   ic=ic+1
                   z1=z1-0.333333333333333_real_8*d1
                   z2=z2-0.333333333333333_real_8*d2
                   IF (d1**2+d2**2.LT.eps2) GOTO 30
                ENDDO
                z1=rho1**(1._real_8/3._real_8)
                z2=rho2**(1._real_8/3._real_8)
30              CONTINUE
                zz1=SQRT(z1)
                zz2=SQRT(z2)
                rs1=rsfac/z1
                rs2=rsfac/z2
                IF (rho1.GT.2.25_real_8*cntr%gceps) THEN
                   smoo1=1._real_8
                   dsmoo1=0._real_8
                ELSEIF (rho1.LT.0.25_real_8*cntr%gceps) THEN
                   smoo1=0._real_8
                   dsmoo1=0._real_8
                ELSE
                   drho=(rho1-0.25_real_8*cntr%gceps)*ogceps
                   smoo1=drho*drho*(0.75_real_8-0.25_real_8*drho)
                   dsmoo1=1.5_real_8*drho*ogceps*(1._real_8-0.5_real_8*drho)
                ENDIF
                IF (rho2.GT.2.25_real_8*cntr%gceps) THEN
                   smoo2=1._real_8
                   dsmoo2=0._real_8
                ELSEIF (rho2.LT.0.25_real_8*cntr%gceps) THEN
                   smoo2=0._real_8
                   dsmoo2=0._real_8
                ELSE
                   drho=(rho2-0.25_real_8*cntr%gceps)*ogceps
                   smoo2=drho*drho*(0.75_real_8-0.25_real_8*drho)
                   dsmoo2=1.5_real_8*drho*ogceps*(1._real_8-0.5_real_8*drho)
                ENDIF
                grhoa1=SQRT(grad(i1,i2+0,i3,1))
                grhoa2=SQRT(grad(i1,i2+1,i3,1))
                grhob1=SQRT(grad(i1,i2+0,i3,5))
                grhob2=SQRT(grad(i1,i2+1,i3,5))
                grhoa1=MAX(grhoa1,small)
                grhoa2=MAX(grhoa2,small)
                grhob1=MAX(grhob1,small)
                grhob2=MAX(grhob2,small)
                xsa1  = grhoa1*(rsa1*rsa1*rsa1*rsa1)
                xsa2  = grhoa2*(rsa2*rsa2*rsa2*rsa2)
                xsb1  = grhob1*(rsb1*rsb1*rsb1*rsb1)
                xsb2  = grhob2*(rsb2*rsb2*rsb2*rsb2)
                IF (rhoa1.LT.2._real_8*small) xsa1=0._real_8
                IF (rhoa2.LT.2._real_8*small) xsa2=0._real_8
                IF (rhob1.LT.2._real_8*small) xsb1=0._real_8
                IF (rhob2.LT.2._real_8*small) xsb2=0._real_8
                xsaa1 = xsa1*xsa1
                xsaa2 = xsa2*xsa2
                xsbb1 = xsb1*xsb1
                xsbb2 = xsb2*xsb2
                sa2a1 = SQRT(1._real_8+xsaa1)
                sa2a2 = SQRT(1._real_8+xsaa2)
                sa2b1 = SQRT(1._real_8+xsbb1)
                sa2b2 = SQRT(1._real_8+xsbb2)
                shma1 = LOG(xsa1+sa2a1)
                shma2 = LOG(xsa2+sa2a2)
                shmb1 = LOG(xsb1+sa2b1)
                shmb2 = LOG(xsb2+sa2b2)
                dda1  = 1._real_8+6._real_8*b1*xsa1*shma1
                dda2  = 1._real_8+6._real_8*b1*xsa2*shma2
                ddb1  = 1._real_8+6._real_8*b1*xsb1*shmb1
                ddb2  = 1._real_8+6._real_8*b1*xsb2*shmb2
                ddda1 = 6._real_8*b1*(shma1+xsa1/sa2a1)
                ddda2 = 6._real_8*b1*(shma2+xsa2/sa2a2)
                dddb1 = 6._real_8*b1*(shmb1+xsb1/sa2b1)
                dddb2 = 6._real_8*b1*(shmb2+xsb2/sa2b2)
                gfa1  = -b1*xsaa1/dda1
                gfa2  = -b1*xsaa2/dda2
                gfb1  = -b1*xsbb1/ddb1
                gfb2  = -b1*xsbb2/ddb2
                dgfa1 = (-2._real_8*b1*xsa1*dda1+b1*xsaa1*ddda1)/(dda1*dda1)
                dgfa2 = (-2._real_8*b1*xsa2*dda2+b1*xsaa2*ddda2)/(dda2*dda2)
                dgfb1 = (-2._real_8*b1*xsb1*ddb1+b1*xsbb1*dddb1)/(ddb1*ddb1)
                dgfb2 = (-2._real_8*b1*xsb2*ddb2+b1*xsbb2*dddb2)/(ddb2*ddb2)
                sxa1  = gfa1*x1*rhoa1
                sxa2  = gfa2*x2*rhoa2
                sxb1  = gfb1*y1*rhob1
                sxb2  = gfb2*y2*rhob2
                v1xa1 = 1.33333333333333_real_8*x1*(gfa1-xsa1*dgfa1)
                v1xa2 = 1.33333333333333_real_8*x2*(gfa2-xsa2*dgfa2)
                v1xb1 = 1.33333333333333_real_8*y1*(gfb1-xsb1*dgfb1)
                v1xb2 = 1.33333333333333_real_8*y2*(gfb2-xsb2*dgfb2)
                v2xa1 = dgfa1/grhoa1
                v2xa2 = dgfa2/grhoa2
                v2xb1 = dgfb1/grhob1
                v2xb2 = dgfb2/grhob2
                ! ..perdew 86 correlation functional
                grho1=grad(i1,i2+0,i3,1)+grad(i1,i2+0,i3,5)+&
                     2._real_8*grad(i1,i2+0,i3,2)*grad(i1,i2+0,i3,6)+&
                     2._real_8*grad(i1,i2+0,i3,3)*grad(i1,i2+0,i3,7)+&
                     2._real_8*grad(i1,i2+0,i3,4)*grad(i1,i2+0,i3,8)
                grho2=grad(i1,i2+1,i3,1)+grad(i1,i2+1,i3,5)+&
                     2._real_8*grad(i1,i2+1,i3,2)*grad(i1,i2+1,i3,6)+&
                     2._real_8*grad(i1,i2+1,i3,3)*grad(i1,i2+1,i3,7)+&
                     2._real_8*grad(i1,i2+1,i3,4)*grad(i1,i2+1,i3,8)
                rs21 = rs1*rs1
                rs22 = rs2*rs2
                rs31 = rs21*rs1
                rs32 = rs22*rs2
                cna1 = pc2+p1*rs1+p2*rs21
                cna2 = pc2+p1*rs2+p2*rs22
                cnb1 = 1._real_8+p3*rs1+p4*rs21+1.e4_real_8*p2*rs31
                cnb2 = 1._real_8+p3*rs2+p4*rs22+1.e4_real_8*p2*rs32
                cn1  = pc1 + cna1/cnb1
                cn2  = pc1 + cna2/cnb2
                drs1 = dfx/(z1*rho1)
                drs2 = dfx/(z2*rho2)
                dcna1= (p1+2._real_8*p2*rs1)*drs1
                dcna2= (p1+2._real_8*p2*rs2)*drs2
                dcnb1= (p3+2._real_8*p4*rs1+3.e4_real_8*p2*rs21)*drs1
                dcnb2= (p3+2._real_8*p4*rs2+3.e4_real_8*p2*rs22)*drs2
                dcn1 = dcna1/cnb1 - cna1/(cnb1*cnb1)*dcnb1
                dcn2 = dcna2/cnb2 - cna2/(cnb2*cnb2)*dcnb2
                s11  = SQRT(x1*x1*rhoa1+y1*y1*rhob1)
                s12  = SQRT(x2*x2*rhoa2+y2*y2*rhob2)
                s21  = two13*zz1/rho1
                s22  = two13*zz2/rho2
                d1   = s11*s21
                d2   = s12*s22
                dda1 = s21*0.83333333333333_real_8*(x1*x1/s11-s11/rho1)
                dda2 = s22*0.83333333333333_real_8*(x2*x2/s12-s12/rho2)
                ddb1 = s21*0.83333333333333_real_8*(y1*y1/s11-s11/rho1)
                ddb2 = s22*0.83333333333333_real_8*(y2*y2/s12-s12/rho2)
                phi1 = 0.192_real_8*pci/cn1*SQRT(grho1)/(zz1*rho1)
                phi2 = 0.192_real_8*pci/cn2*SQRT(grho2)/(zz2*rho2)
                ephi1= EXP(-phi1)
                ephi2= EXP(-phi2)
                sc1  = grho1/(z1*rho1)*cn1*ephi1/d1
                sc2  = grho2/(z2*rho2)*cn2*ephi2/d2
                v1c1 = sc1*((1._real_8+phi1)*dcn1/cn1-(1.33333333333333_real_8-&
                     1.1666666666666_real_8*phi1)/rho1)
                v1c2 = sc2*((1._real_8+phi2)*dcn2/cn2-(1.33333333333333_real_8-&
                     1.1666666666666_real_8*phi2)/rho2)
                v1ca1= -sc1/d1*dda1+v1c1
                v1ca2= -sc2/d2*dda2+v1c2
                v1cb1= -sc1/d1*ddb1+v1c1
                v1cb2= -sc2/d2*ddb2+v1c2
                v2c1 = cn1*ephi1/(z1*rho1)*(2._real_8-phi1)/d1
                v2c2 = cn2*ephi2/(z2*rho2)*(2._real_8-phi2)/d2
                v2ca1= v2c1
                v2ca2= v2c2
                v2cb1= v2c1
                v2cb2= v2c2
                v2cab1= v2c1
                v2cab2= v2c2
                ! ..sum up terms
                sgcx = sgcx + (smoo1*(sxa1+sxb1) + smoo2*(sxa2+sxb2))
                sgcc = sgcc + (smoo1*sc1 + smoo2*sc2)
                v1(i1,i2+0,i3)=dsmoo1*(sxa1+sc1)+smoo1*(v1xa1+v1ca1)
                v1(i1,i2+1,i3)=dsmoo2*(sxa2+sc2)+smoo2*(v1xa2+v1ca2)
                v2(i1,i2+0,i3)=dsmoo1*(sxb1+sc1)+smoo1*(v1xb1+v1cb1)
                v2(i1,i2+1,i3)=dsmoo2*(sxb2+sc2)+smoo2*(v1xb2+v1cb2)
                vtmp(i1,i2+0,i3,1)=smoo1*(v2xa1+v2ca1)
                vtmp(i1,i2+1,i3,1)=smoo2*(v2xa2+v2ca2)
                vtmp(i1,i2+0,i3,2)=smoo1*(v2xb1+v2cb1)
                vtmp(i1,i2+1,i3,2)=smoo2*(v2xb2+v2cb2)
                vtmp(i1,i2+0,i3,3)=smoo1*v2cab1
                vtmp(i1,i2+1,i3,3)=smoo2*v2cab2
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    ! ..counting floating point operations
    fadd   = ic * 11._real_8 + id *  78._real_8
    fmult  = ic * 12._real_8 + id * 204._real_8
    fdiv   = ic *  4._real_8 + id *  52._real_8
    fspez  =              id *  16._real_8
    flops=fadd+fmult+fdiv+fspez
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gcsxp86
  ! ==================================================================
  SUBROUTINE gcsxlyp(b1,sgcx,sgcc,rhoe,v1,v2,vtmp,grad,flops)
    REAL(real_8)                             :: b1, sgcx, sgcc
    REAL(real_8), DIMENSION(fpar%kr1, fpar%&
      kr2s, fpar%kr3s, 2)                    :: rhoe
    COMPLEX(real_8) :: v1(fpar%kr1,fpar%kr2s,fpar%kr3s), &
      v2(fpar%kr1,fpar%kr2s,fpar%kr3s)
    REAL(real_8), DIMENSION(fpar%kr1, fpar%&
      kr2s, fpar%kr3s, 3)                    :: vtmp
    REAL(real_8), DIMENSION(fpar%kr1, fpar%&
      kr2s, fpar%kr3s, 8)                    :: grad
    REAL(real_8)                             :: flops

    REAL(real_8), PARAMETER :: a = 0.04918_real_8, b = 0.132_real_8, &
      c = 0.2533_real_8, d = 0.349_real_8 , eps2 = 1.e-20_real_8, &
      small = 1.e-24_real_8 

    INTEGER                                  :: i1, i2, i3, ic, id, ii2, &
                                                ncount
    REAL(real_8) :: d1, d2, dda1, dda2, ddb1, ddb2, ddda1, ddda2, dddb1, &
      dddb2, dder1, dder2, der1, der2, dgfa1, dgfa2, dgfb1, dgfb2, dlaa1, &
      dlaa2, dlaaa1, dlaaa2, dlaab1, dlaab2, dlab1, dlab2, dlaba1, dlaba2, &
      dlabb1, dlabb2, dlbb1, dlbb2, dlbba1, dlbba2, dlbbb1, dlbbb2, dor1, &
      dor2, dr1, dr2, drho, dsmoo1, dsmoo2, fadd, fdiv, fmult, fspez, gfa1, &
      gfa2, gfb1, gfb2, grhoa1, grhoa2, grhoaa1, grhoaa2, grhoab1, grhoab2, &
      grhob1, grhob2, grhobb1, grhobb2, ogceps, or1, or2, rho1, rho2, rhoa1, &
      rhoa2, rhob1, rhob2, rs1, rs2, rsa1, rsa2, rsb1, rsb2, sa2a1, sa2a2, &
      sa2b1, sa2b2, sc1, sc2, shma1, shma2, shmb1
    REAL(real_8) :: shmb2, smoo1, smoo2, sxa1, sxa2, sxb1, sxb2, v1ca1, &
      v1ca2, v1cb1, v1cb2, v1xa1, v1xa2, v1xb1, v1xb2, v2ca1, v2ca2, v2cab1, &
      v2cab2, v2cb1, v2cb2, v2xa1, v2xa2, v2xb1, v2xb2, x1, x2, xsa1, xsa2, &
      xsaa1, xsaa2, xsb1, xsb2, xsbb1, xsbb2, y1, y2, z1, z2

! ==--------------------------------------------------------------==

    ogceps=1._real_8/cntr%gceps
    sgcx=0._real_8
    sgcc=0._real_8
    x1=1._real_8
    x2=1._real_8
    y1=1._real_8
    y2=1._real_8
    z1=1._real_8
    z2=1._real_8
    ic=0
    id=0
    ii2=MOD(parm%nr2,2)
    IF (ii2.NE.0) CALL stopgm('GCSXLYP','ODD DIMENSIONS',& 
         __LINE__,__FILE__)
    !$omp parallel do private(I1,I2,I3,RHOA1,RHOA2,D1,D2) &
    !$omp  firstprivate(X1,X2,Y1,Y2,Z1,Z2) &
    !$omp  private(NCOUNT,RHOB1, RHOB2, RHO1, RHO2, RSA1, RSA2) &
    !$omp  private(RSB1, RSB2, RS1, RS2, SMOO1, DSMOO1, DRHO) &
    !$omp  private(SMOO2, DSMOO2, GRHOAA1, GRHOAA2, GRHOAB1, GRHOAB2) &
    !$omp  private(GRHOBB1, GRHOBB2, GRHOA1, GRHOA2, GRHOB1, GRHOB2) &
    !$omp  private(XSA1, XSA2, XSB1, XSB2, XSAA1, XSAA2, XSBB1, XSBB2) &
    !$omp  private(SA2A1, SA2A2, SA2B1, SA2B2, SHMA1, SHMA2, SHMB1, SHMB2) &
    !$omp  private(DDA1, DDA2, DDB1, DDB2, DDDA1, DDDA2, DDDB1, DDDB2) &
    !$omp  private(GFA1, GFA2, GFB1, GFB2, DGFA1, DGFA2, DGFB1, DGFB2, SXA1) &
    !$omp  private(SXA2, SXB1, SXB2, V1XA1, V1XA2,V1XB1,V1XB2,V2XA1,V2XA2) &
    !$omp  private(V2XB1, V2XB2, DR1, DR2, OR1, OR2, DOR1, DOR2, DER1, DER2) &
    !$omp  private(DDER1, DDER2, DLAA1, DLAA2, DLAB1, DLAB2, DLBB1, DLBB2) &
    !$omp  private(DLAAA1, DLAAA2, DLAAB1, DLAAB2, DLABA1, DLABA2, DLABB1) &
    !$omp  private(DLABB2, DLBBA1, DLBBA2, DLBBB1, DLBBB2, SC1, SC2, V1CA1) &
    !$omp  private(V1CA2, V1CB1, V1CB2, V2CA1, V2CA2, V2CB1, V2CB2) &
    !$omp  private(V2CAB1,V2CAB2) &
    !$omp  reduction(+:SGCX,SGCC,IC,ID) __COLLAPSE3
    DO i3=1,parm%nr3
       DO i2=1,parm%nr2,2
          DO i1=1,parm%nr1
             rhoa1=MAX(rhoe(i1,i2+0,i3,1),small)
             rhoa2=MAX(rhoe(i1,i2+1,i3,1),small)
             rhob1=MAX(rhoe(i1,i2+0,i3,2),small)
             rhob2=MAX(rhoe(i1,i2+1,i3,2),small)
             rho1=rhoa1+rhob1
             rho2=rhoa2+rhob2
             IF ((rho1+rho2).GT.5._real_8*small) THEN
                id=id+1
                DO ncount=1,100
                   d1=x1-rhoa1/x1**2
                   d2=x2-rhoa2/x2**2
                   ic=ic+1
                   x1=x1-0.333333333333333_real_8*d1
                   x2=x2-0.333333333333333_real_8*d2
                   IF (d1**2+d2**2.LT.eps2) GOTO 10
                ENDDO
                x1=rhoa1**(1._real_8/3._real_8)
                x2=rhoa2**(1._real_8/3._real_8)
10              CONTINUE
                rsa1=1._real_8/x1
                rsa2=1._real_8/x2
                DO ncount=1,100
                   d1=y1-rhob1/y1**2
                   d2=y2-rhob2/y2**2
                   ic=ic+1
                   y1=y1-0.333333333333333_real_8*d1
                   y2=y2-0.333333333333333_real_8*d2
                   IF (d1**2+d2**2.LT.eps2) GOTO 20
                ENDDO
                y1=rhob1**(1._real_8/3._real_8)
                y2=rhob2**(1._real_8/3._real_8)
20              CONTINUE
                rsb1=1._real_8/y1
                rsb2=1._real_8/y2
                DO ncount=1,100
                   d1=z1-rho1/z1**2
                   d2=z2-rho2/z2**2
                   ic=ic+1
                   z1=z1-0.333333333333333_real_8*d1
                   z2=z2-0.333333333333333_real_8*d2
                   IF (d1**2+d2**2.LT.eps2) GOTO 30
                ENDDO
                z1=rho1**(1._real_8/3._real_8)
                z2=rho2**(1._real_8/3._real_8)
30              CONTINUE
                rs1=1._real_8/z1
                rs2=1._real_8/z2
                IF (rho1.GT.2.25_real_8*cntr%gceps) THEN
                   smoo1=1._real_8
                   dsmoo1=0._real_8
                ELSEIF (rho1.LT.0.25_real_8*cntr%gceps) THEN
                   smoo1=0._real_8
                   dsmoo1=0._real_8
                ELSE
                   drho=(rho1-0.25_real_8*cntr%gceps)*ogceps
                   smoo1=drho*drho*(0.75_real_8-0.25_real_8*drho)
                   dsmoo1=1.5_real_8*drho*ogceps*(1._real_8-0.5_real_8*drho)
                ENDIF
                IF (rho2.GT.2.25_real_8*cntr%gceps) THEN
                   smoo2=1._real_8
                   dsmoo2=0._real_8
                ELSEIF (rho2.LT.0.25_real_8*cntr%gceps) THEN
                   smoo2=0._real_8
                   dsmoo2=0._real_8
                ELSE
                   drho=(rho2-0.25_real_8*cntr%gceps)*ogceps
                   smoo2=drho*drho*(0.75_real_8-0.25_real_8*drho)
                   dsmoo2=1.5_real_8*drho*ogceps*(1._real_8-0.5_real_8*drho)
                ENDIF
                grhoaa1=grad(i1,i2+0,i3,1)
                grhoaa2=grad(i1,i2+1,i3,1)
                grhoab1=grad(i1,i2+0,i3,2)*grad(i1,i2+0,i3,6)+&
                     grad(i1,i2+0,i3,3)*grad(i1,i2+0,i3,7)+&
                     grad(i1,i2+0,i3,4)*grad(i1,i2+0,i3,8)
                grhoab2=grad(i1,i2+1,i3,2)*grad(i1,i2+1,i3,6)+&
                     grad(i1,i2+1,i3,3)*grad(i1,i2+1,i3,7)+&
                     grad(i1,i2+1,i3,4)*grad(i1,i2+1,i3,8)
                grhobb1=grad(i1,i2+0,i3,5)
                grhobb2=grad(i1,i2+1,i3,5)
                grhoa1=SQRT(grad(i1,i2+0,i3,1))
                grhoa2=SQRT(grad(i1,i2+1,i3,1))
                grhob1=SQRT(grad(i1,i2+0,i3,5))
                grhob2=SQRT(grad(i1,i2+1,i3,5))
                grhoa1=MAX(grhoa1,small)
                grhoa2=MAX(grhoa2,small)
                grhob1=MAX(grhob1,small)
                grhob2=MAX(grhob2,small)
                xsa1  = grhoa1*(rsa1*rsa1*rsa1*rsa1)
                xsa2  = grhoa2*(rsa2*rsa2*rsa2*rsa2)
                xsb1  = grhob1*(rsb1*rsb1*rsb1*rsb1)
                xsb2  = grhob2*(rsb2*rsb2*rsb2*rsb2)
                IF (rhoa1.LT.2._real_8*small) xsa1=0._real_8
                IF (rhoa2.LT.2._real_8*small) xsa2=0._real_8
                IF (rhob1.LT.2._real_8*small) xsb1=0._real_8
                IF (rhob2.LT.2._real_8*small) xsb2=0._real_8
                xsaa1 = xsa1*xsa1
                xsaa2 = xsa2*xsa2
                xsbb1 = xsb1*xsb1
                xsbb2 = xsb2*xsb2
                sa2a1 = SQRT(1._real_8+xsaa1)
                sa2a2 = SQRT(1._real_8+xsaa2)
                sa2b1 = SQRT(1._real_8+xsbb1)
                sa2b2 = SQRT(1._real_8+xsbb2)
                shma1 = LOG(xsa1+sa2a1)
                shma2 = LOG(xsa2+sa2a2)
                shmb1 = LOG(xsb1+sa2b1)
                shmb2 = LOG(xsb2+sa2b2)
                dda1  = 1._real_8+6._real_8*b1*xsa1*shma1
                dda2  = 1._real_8+6._real_8*b1*xsa2*shma2
                ddb1  = 1._real_8+6._real_8*b1*xsb1*shmb1
                ddb2  = 1._real_8+6._real_8*b1*xsb2*shmb2
                ddda1 = 6._real_8*b1*(shma1+xsa1/sa2a1)
                ddda2 = 6._real_8*b1*(shma2+xsa2/sa2a2)
                dddb1 = 6._real_8*b1*(shmb1+xsb1/sa2b1)
                dddb2 = 6._real_8*b1*(shmb2+xsb2/sa2b2)
                gfa1  = -b1*xsaa1/dda1
                gfa2  = -b1*xsaa2/dda2
                gfb1  = -b1*xsbb1/ddb1
                gfb2  = -b1*xsbb2/ddb2
                dgfa1 = (-2._real_8*b1*xsa1*dda1+b1*xsaa1*ddda1)/(dda1*dda1)
                dgfa2 = (-2._real_8*b1*xsa2*dda2+b1*xsaa2*ddda2)/(dda2*dda2)
                dgfb1 = (-2._real_8*b1*xsb1*ddb1+b1*xsbb1*dddb1)/(ddb1*ddb1)
                dgfb2 = (-2._real_8*b1*xsb2*ddb2+b1*xsbb2*dddb2)/(ddb2*ddb2)
                sxa1  = gfa1*x1*rhoa1
                sxa2  = gfa2*x2*rhoa2
                sxb1  = gfb1*y1*rhob1
                sxb2  = gfb2*y2*rhob2
                v1xa1 = 1.33333333333333_real_8*x1*(gfa1-xsa1*dgfa1)
                v1xa2 = 1.33333333333333_real_8*x2*(gfa2-xsa2*dgfa2)
                v1xb1 = 1.33333333333333_real_8*y1*(gfb1-xsb1*dgfb1)
                v1xb2 = 1.33333333333333_real_8*y2*(gfb2-xsb2*dgfb2)
                v2xa1 = dgfa1/grhoa1
                v2xa2 = dgfa2/grhoa2
                v2xb1 = dgfb1/grhob1
                v2xb2 = dgfb2/grhob2
                ! ..LYP correlation functional
                dr1=(1._real_8+d*rs1)
                dr2=(1._real_8+d*rs2)
                or1=EXP(-c*rs1)/dr1*(rs1**11)
                or1=MAX(or1,small)
                or2=EXP(-c*rs2)/dr2*(rs2**11)
                or2=MAX(or2,small)
                dor1=-0.3333333333333_real_8*(rs1**4)*or1*(11._real_8/rs1-c-d/dr1)
                dor2=-0.3333333333333_real_8*(rs2**4)*or2*(11._real_8/rs2-c-d/dr2)
                der1=c*rs1+d*rs1/dr1
                der2=c*rs2+d*rs2/dr2
                dder1=0.3333333333333_real_8*(d*d*(rs1**5)/(dr1*dr1)-der1/rho1)
                dder2=0.3333333333333_real_8*(d*d*(rs2**5)/(dr2*dr2)-der2/rho2)
                dlaa1=-a*b*or1*(0.1111111111111_real_8*rhoa1*rhob1*(1._real_8-&
                     3._real_8*der1-(der1-11._real_8)*rhoa1/rho1)-rhob1*rhob1)
                dlaa2=-a*b*or2*(0.1111111111111_real_8*rhoa2*rhob2*(1._real_8-&
                     3._real_8*der2-(der2-11._real_8)*rhoa2/rho2)-rhob2*rhob2)
                dlab1=-a*b*or1*(0.1111111111111_real_8*rhoa1*rhob1*(47._real_8-&
                     7._real_8*der1)-1.3333333333333_real_8*rho1*rho1)
                dlab2=-a*b*or2*(0.1111111111111_real_8*rhoa2*rhob2*(47._real_8-&
                     7._real_8*der2)-1.3333333333333_real_8*rho2*rho2)
                dlbb1=-a*b*or1*(0.1111111111111_real_8*rhoa1*rhob1*(1._real_8-&
                     3._real_8*der1-(der1-11._real_8)*rhob1/rho1)-rhoa1*rhoa1)
                dlbb2=-a*b*or2*(0.1111111111111_real_8*rhoa2*rhob2*(1._real_8-&
                     3._real_8*der2-(der2-11._real_8)*rhob2/rho2)-rhoa2*rhoa2)
                dlaaa1=dor1/or1*dlaa1-a*b*or1*(0.1111111111111_real_8*rhob1*&
                     (1._real_8-3._real_8*der1-(der1-11._real_8)*rhoa1/rho1)-&
                     0.1111111111111_real_8*rhoa1*rhob1*((3._real_8+rhoa1/rho1)*dder1+&
                     (der1-11._real_8)*rhob1/(rho1*rho1)))
                dlaaa2=dor2/or2*dlaa2-a*b*or2*(0.1111111111111_real_8*rhob2*&
                     (1._real_8-3._real_8*der2-(der2-11._real_8)*rhoa2/rho2)-&
                     0.1111111111111_real_8*rhoa2*rhob2*((3._real_8+rhoa2/rho2)*dder2+&
                     (der2-11._real_8)*rhob2/(rho2*rho2)))
                dlaab1=dor1/or1*dlaa1-a*b*or1*(0.1111111111111_real_8*rhoa1*&
                     (1._real_8-3._real_8*der1-(der1-11._real_8)*rhoa1/rho1)-&
                     0.1111111111111_real_8*rhoa1*rhob1*((3._real_8+rhoa1/rho1)*dder1-&
                     (der1-11._real_8)*rhoa1/(rho1*rho1))-2._real_8*rhob1)
                dlaab2=dor2/or2*dlaa2-a*b*or2*(0.1111111111111_real_8*rhoa2*&
                     (1._real_8-3._real_8*der2-(der2-11._real_8)*rhoa2/rho2)-&
                     0.1111111111111_real_8*rhoa2*rhob2*((3._real_8+rhoa2/rho2)*dder2-&
                     (der2-11._real_8)*rhoa2/(rho2*rho2))-2._real_8*rhob2)
                dlaba1=dor1/or1*dlab1-a*b*or1*(0.1111111111111_real_8*rhob1*&
                     (47._real_8-7._real_8*der1)-0.7777777777778_real_8*rhoa1*rhob1*dder1-&
                     2.66666666666667_real_8*rho1)
                dlaba2=dor2/or2*dlab2-a*b*or2*(0.1111111111111_real_8*rhob2*&
                     (47._real_8-7._real_8*der2)-0.7777777777778_real_8*rhoa2*rhob2*dder2-&
                     2.66666666666667_real_8*rho2)
                dlabb1=dor1/or1*dlab1-a*b*or1*(0.1111111111111_real_8*rhoa1*&
                     (47._real_8-7._real_8*der1)-0.7777777777778_real_8*rhoa1*rhob1*dder1-&
                     2.66666666666667_real_8*rho1)
                dlabb2=dor2/or2*dlab2-a*b*or2*(0.1111111111111_real_8*rhoa2*&
                     (47._real_8-7._real_8*der2)-0.7777777777778_real_8*rhoa2*rhob2*dder2-&
                     2.66666666666667_real_8*rho2)
                dlbba1=dor1/or1*dlbb1-a*b*or1*(0.1111111111111_real_8*rhob1*&
                     (1._real_8-3._real_8*der1-(der1-11._real_8)*rhob1/rho1)-&
                     0.1111111111111_real_8*rhoa1*rhob1*((3._real_8+rhob1/rho1)*dder1-&
                     (der1-11._real_8)*rhob1/(rho1*rho1))-2._real_8*rhoa1)
                dlbba2=dor2/or2*dlbb2-a*b*or2*(0.1111111111111_real_8*rhob2*&
                     (1._real_8-3._real_8*der2-(der2-11._real_8)*rhob2/rho2)-&
                     0.1111111111111_real_8*rhoa2*rhob2*((3._real_8+rhob2/rho2)*dder2-&
                     (der2-11._real_8)*rhob2/(rho2*rho2))-2._real_8*rhoa2)
                dlbbb1=dor1/or1*dlbb1-a*b*or1*(0.1111111111111_real_8*rhoa1*&
                     (1._real_8-3._real_8*der1-(der1-11._real_8)*rhob1/rho1)-&
                     0.1111111111111_real_8*rhoa1*rhob1*((3._real_8+rhob1/rho1)*dder1+&
                     (der1-11._real_8)*rhoa1/(rho1*rho1)))
                dlbbb2=dor2/or2*dlbb2-a*b*or2*(0.1111111111111_real_8*rhoa2*&
                     (1._real_8-3._real_8*der2-(der2-11._real_8)*rhob2/rho2)-&
                     0.1111111111111_real_8*rhoa2*rhob2*((3._real_8+rhob2/rho2)*dder2+&
                     (der2-11._real_8)*rhoa2/(rho2*rho2)))
                sc1=dlaa1*grhoaa1+dlab1*grhoab1+dlbb1*grhobb1
                sc2=dlaa2*grhoaa2+dlab2*grhoab2+dlbb2*grhobb2
                v1ca1=dlaaa1*grhoaa1+dlaba1*grhoab1+dlbba1*grhobb1
                v1ca2=dlaaa2*grhoaa2+dlaba2*grhoab2+dlbba2*grhobb2
                v1cb1=dlaab1*grhoaa1+dlabb1*grhoab1+dlbbb1*grhobb1
                v1cb2=dlaab2*grhoaa2+dlabb2*grhoab2+dlbbb2*grhobb2
                v2ca1=2._real_8*dlaa1
                v2ca2=2._real_8*dlaa2
                v2cb1=2._real_8*dlbb1
                v2cb2=2._real_8*dlbb2
                v2cab1=dlab1
                v2cab2=dlab2
                ! ..sum up terms
                sgcx = sgcx + (smoo1*(sxa1+sxb1) + smoo2*(sxa2+sxb2))
                sgcc = sgcc + (smoo1*sc1 + smoo2*sc2)
                v1(i1,i2+0,i3)=dsmoo1*(sxa1+sc1)+smoo1*(v1xa1+v1ca1)
                v1(i1,i2+1,i3)=dsmoo2*(sxa2+sc2)+smoo2*(v1xa2+v1ca2)
                v2(i1,i2+0,i3)=dsmoo1*(sxb1+sc1)+smoo1*(v1xb1+v1cb1)
                v2(i1,i2+1,i3)=dsmoo2*(sxb2+sc2)+smoo2*(v1xb2+v1cb2)
                vtmp(i1,i2+0,i3,1)=smoo1*(v2xa1+v2ca1)
                vtmp(i1,i2+1,i3,1)=smoo2*(v2xa2+v2ca2)
                vtmp(i1,i2+0,i3,2)=smoo1*(v2xb1+v2cb1)
                vtmp(i1,i2+1,i3,2)=smoo2*(v2xb2+v2cb2)
                vtmp(i1,i2+0,i3,3)=smoo1*v2cab1
                vtmp(i1,i2+1,i3,3)=smoo2*v2cab2
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    ! ..counting floating point operations
    fadd   = ic * 11._real_8 + id * 170._real_8
    fmult  = ic * 12._real_8 + id * 372._real_8
    fdiv   = ic *  4._real_8 + id *  76._real_8
    fspez  =              id *  20._real_8
    flops=fadd+fmult+fdiv+fspez
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gcsxlyp
  ! ==================================================================
  SUBROUTINE gcspbe(sgcx,sgcc,rhoe,v1,v2,vtmp,grad,flops)
    REAL(real_8)                             :: sgcx, sgcc
    REAL(real_8), DIMENSION(fpar%kr1, fpar%&
      kr2s, fpar%kr3s, 2)                    :: rhoe
    COMPLEX(real_8) :: v1(fpar%kr1,fpar%kr2s,fpar%kr3s), &
      v2(fpar%kr1,fpar%kr2s,fpar%kr3s)
    REAL(real_8), DIMENSION(fpar%kr1, fpar%&
      kr2s, fpar%kr3s, 3)                    :: vtmp
    REAL(real_8), DIMENSION(fpar%kr1, fpar%&
      kr2s, fpar%kr3s, 8)                    :: grad
    REAL(real_8)                             :: flops

    REAL(real_8), PARAMETER :: a0u = 0.4581652932831429_real_8, &
      a1u = 2.217058676663745_real_8, a2u = 0.7405551735357053_real_8, &
      a3u = 0.01968227878617998_real_8 , ax = -0.738558766382022406_real_8, &
      b1u = 1.0_real_8, b2u = 4.504130959426697_real_8, &
      b3u = 1.110667363742916_real_8, b4u = 0.02359291751427506_real_8 , &
      be = 0.06672455060314922_real_8, da0 = 0.119086804055547_real_8, &
      da1 = 0.6157402568883345_real_8, da2 = 0.1574201515892867_real_8, &
      da3 = 0.003532336663397157_real_8 , db1 = 0.0_real_8, &
      db2 = 0.2673612973836267_real_8, db3 = 0.2052004607777787_real_8, &
      db4 = 0.004200005045691381_real_8 
    REAL(real_8), PARAMETER :: eps2 = 1.e-20_real_8, &
      f1 = -1.39578858466194911_real_8 , ga = 0.031090690869654895_real_8 , &
      rsfac = 0.6203504908994000_real_8 , small = 1.e-24_real_8, &
      ssmall = 1.e-12_real_8, t1 = 0.854960467080683_real_8, &
      t2 = 0.07916300621117439_real_8, t3 = 0.02580127609845684_real_8, &
      t4 = 0.0121839359353824_real_8, t5 = 0.006919272259599879_real_8, &
      t6 = 0.004391524649611372_real_8, t7 = 0.003002751897170168_real_8, &
      t8 = 0.00216587382212552_real_8, t9 = 0.01141189204579585_real_8 , &
      two3 = 2._real_8/3._real_8 , uk = 0.8040_real_8, &
      um = 0.2195149727645171_real_8, ul = um/uk 
    REAL(real_8), PARAMETER :: us = 0.161620459673995492_real_8

    INTEGER                                  :: i1, i2, i3, ic, id, ii2, &
                                                ncount
    REAL(real_8) :: a0p1, a0p2, a1p1, a1p2, a2p1, a2p2, a3p1, a3p2, aa1, aa2, &
      af1, af2, as1, as2, b1p1, b1p2, b2p1, b2p2, b3p1, b3p2, b4p1, b4p2, &
      bot1, bot2, d1, d2, dadra1, dadra2, dadrb1, dadrb2, dbot1, dbot2, &
      dedra1, dedra2, dedrb1, dedrb2, dfx1, dfx2, dfxs1, dfxs2, dhdra1, &
      dhdra2, dhdrb1, dhdrb2, dhdt1, dhdt2, dphida1, dphida2, dphidb1, &
      dphidb2, dphide1, dphide2, drho, dsda1, dsda2, dsdra1, dsdra2, dsdrb1, &
      dsdrb2, dsdt1, dsdt2, dsmoo1, dsmoo2, dtdra1, dtdra2, dtdrb1, dtdrb2, &
      dtop1, dtop2, dvxc1, dvxc2, dxa1, dxa2, dxb1, dxb2, epsxcu1, epsxcu2, &
      esla1, esla2, et1, et2, eta1, eta2, ex1, ex2
    REAL(real_8) :: expe1, expe2, f1s, f1u, fadd, fdiv, fmult, fspez, fx1, &
      fx2, fxs1, fxs2, grho1, grho2, grhoa1, grhoa2, grhob1, grhob2, h01, &
      h02, ob3, obot1, obot2, ogceps, phi1, phi2, phi31, phi32, po1, po2, &
      rho1, rho2, rhoa1, rhoa2, rhob1, rhob2, rr1, rr2, rs1, rs2, rsa1, rsa2, &
      rsb1, rsb2, s11, s12, s21, s22, sc1, sc2, smoo1, smoo2, sxa1, sxa2, &
      sxb1, sxb2, t13, t43, top1, top2, tt1, tt2, v1ca1, v1ca2, v1cb1, v1cb2, &
      v1xa1, v1xa2, v1xb1, v1xb2, v2ca1, v2ca2, v2cab1, v2cab2, v2cb1, v2cb2, &
      v2xa1, v2xa2, v2xb1, v2xb2, vsla1, vsla2, vslb1, vslb2, vxc1, vxc2, x1, &
      x2, xbot1, xbot2, xkf1, xkf2, xks1
    REAL(real_8) :: xks2, xtop1, xtop2, xy1, xy2, y1, y2, yt1, yt2, z1, z2

    t43=2._real_8**(-4._real_8/3._real_8)
    t13=2._real_8**(1._real_8/3._real_8)
    ob3=1._real_8/3._real_8
    f1s=f1*2._real_8/3._real_8
    f1u=f1s*1.3333333333333333_real_8
    ogceps=1._real_8/cntr%gceps
    sgcx=0._real_8
    sgcc=0._real_8
    x1=1._real_8
    x2=1._real_8
    y1=1._real_8
    y2=1._real_8
    z1=1._real_8
    z2=1._real_8
    ic=0
    id=0
    ii2=MOD(parm%nr2,2)
    IF (ii2.NE.0) CALL stopgm('GCSPBE','ODD DIMENSIONS',& 
         __LINE__,__FILE__)
    !$omp parallel do private(I1,I2,I3,RHOA1,RHOA2,D1,D2) &
    !$omp  firstprivate(X1,X2,Y1,Y2,Z1,Z2) &
    !$omp  private(NCOUNT,RHOB1, RHOB2, RHO1, RHO2, RSA1, RSA2) &
    !$omp  private(RSB1, RSB2, RS1, RS2, SMOO1, DSMOO1, DRHO) &
    !$omp  private(SMOO2, DSMOO2, GRHOA1, GRHOA2, GRHOB1, GRHOB2, AA1) &
    !$omp  private(AA2, RR1, RR2, EX1, EX2, S21, S22, PO1, PO2, FX1) &
    !$omp  private(FX2, SXA1, SXA2, DFX1, DFX2, V1XA1, V1XA2, V2XA1, V2XA2) &
    !$omp  private(SXB1, SXB2, V1XB1, V1XB2, V2XB1, V2XB2, ETA1, ETA2) &
    !$omp  private(ET1, FXS1, DFXS1, ET2, FXS2, DFXS2, A0P1, A1P1, A2P1) &
    !$omp  private(A3P1, B2P1, B3P1, B4P1, A0P2, A1P2, A2P2, A3P2, B2P2) &
    !$omp  private(B3P2, B4P2, TOP1, TOP2, DTOP1, DTOP2, XTOP1, XTOP2, BOT1) &
    !$omp  private(BOT2, DBOT1, DBOT2, XBOT1, XBOT2, OBOT1, OBOT2, EPSXCU1) &
    !$omp  private(EPSXCU2, VXC1, DVXC1, DXA1, DXB1, VXC2, DVXC2, DXA2) &
    !$omp  private(DXB2, ESLA1, ESLA2, VSLA1, VSLB1, VSLA2, VSLB2, GRHO1) &
    !$omp  private(GRHO2, AS1, AS2, PHI1, PHI31, PHI2, PHI32, XKF1, XKF2) &
    !$omp  private(XKS1, XKS2, TT1, TT2, EXPE1, EXPE2, AF1, AF2, YT1, YT2) &
    !$omp  private(XY1, XY2, S11, S12, H01, H02, SC1, SC2, DPHIDE1, DPHIDE2) &
    !$omp  private(DEDRA1, DEDRB1, DEDRA2, DEDRB2, DPHIDA1,DPHIDB1,DPHIDA2) &
    !$omp  private(DPHIDB2, DTDRA1, DTDRB1, DTDRA2, DTDRB2, DADRA1, DADRB1) &
    !$omp  private(DADRA2, DADRB2, DSDA1, DSDA2, DSDT1,DSDT2,DSDRA1,DSDRB1) &
    !$omp  private(DSDRA2, DSDRB2, DHDT1, DHDT2, DHDRA1, DHDRB1, DHDRA2) &
    !$omp  private(DHDRB2, V1CA1, V2CA1, V1CB1, V2CB1, V2CAB1, V1CA2) &
    !$omp  private(V2CA2, V1CB2, V2CB2, V2CAB2) &
    !$omp  private(B1P1,B1P2) &
    !$omp  reduction(+:SGCX,SGCC,IC,ID) __COLLAPSE3
    DO i3=1,parm%nr3
       DO i2=1,parm%nr2,2
          DO i1=1,parm%nr1
             rhoa1=MAX(rhoe(i1,i2+0,i3,1),small)
             rhoa2=MAX(rhoe(i1,i2+1,i3,1),small)
             rhob1=MAX(rhoe(i1,i2+0,i3,2),small)
             rhob2=MAX(rhoe(i1,i2+1,i3,2),small)
             rho1=rhoa1+rhob1
             rho2=rhoa2+rhob2
             IF ((rho1+rho2).GT.5._real_8*small) THEN
                id=id+1
                DO ncount=1,100
                   d1=x1-rhoa1/x1**2
                   d2=x2-rhoa2/x2**2
                   ic=ic+1
                   x1=x1-0.333333333333333_real_8*d1
                   x2=x2-0.333333333333333_real_8*d2
                   IF (d1**2+d2**2.LT.eps2) GOTO 10
                ENDDO
                x1=rhoa1**(1._real_8/3._real_8)
                x2=rhoa2**(1._real_8/3._real_8)
10              CONTINUE
                rsa1=1._real_8/x1
                rsa2=1._real_8/x2
                DO ncount=1,100
                   d1=y1-rhob1/y1**2
                   d2=y2-rhob2/y2**2
                   ic=ic+1
                   y1=y1-0.333333333333333_real_8*d1
                   y2=y2-0.333333333333333_real_8*d2
                   IF (d1**2+d2**2.LT.eps2) GOTO 20
                ENDDO
                y1=rhob1**(1._real_8/3._real_8)
                y2=rhob2**(1._real_8/3._real_8)
20              CONTINUE
                rsb1=1._real_8/y1
                rsb2=1._real_8/y2
                DO ncount=1,100
                   d1=z1-rho1/z1**2
                   d2=z2-rho2/z2**2
                   ic=ic+1
                   z1=z1-0.333333333333333_real_8*d1
                   z2=z2-0.333333333333333_real_8*d2
                   IF (d1**2+d2**2.LT.eps2) GOTO 30
                ENDDO
                z1=rho1**(1._real_8/3._real_8)
                z2=rho2**(1._real_8/3._real_8)
30              CONTINUE
                rs1=rsfac/z1
                rs2=rsfac/z2
                IF (rho1.GT.2.25_real_8*cntr%gceps) THEN
                   smoo1=1._real_8
                   dsmoo1=0._real_8
                ELSEIF (rho1.LT.0.25_real_8*cntr%gceps) THEN
                   smoo1=0._real_8
                   dsmoo1=0._real_8
                ELSE
                   drho=(rho1-0.25_real_8*cntr%gceps)*ogceps
                   smoo1=drho*drho*(0.75_real_8-0.25_real_8*drho)
                   dsmoo1=1.5_real_8*drho*ogceps*(1._real_8-0.5_real_8*drho)
                ENDIF
                IF (rho2.GT.2.25_real_8*cntr%gceps) THEN
                   smoo2=1._real_8
                   dsmoo2=0._real_8
                ELSEIF (rho2.LT.0.25_real_8*cntr%gceps) THEN
                   smoo2=0._real_8
                   dsmoo2=0._real_8
                ELSE
                   drho=(rho2-0.25_real_8*cntr%gceps)*ogceps
                   smoo2=drho*drho*(0.75_real_8-0.25_real_8*drho)
                   dsmoo2=1.5_real_8*drho*ogceps*(1._real_8-0.5_real_8*drho)
                ENDIF

                IF (rhoa1.GT.small) THEN
                   grhoa1=grad(i1,i2+0,i3,1)
                   grhoa1=MAX(grhoa1,small)
                ELSE
                   grhoa1=small
                ENDIF
                IF (rhoa2.GT.small) THEN
                   grhoa2=grad(i1,i2+1,i3,1)
                   grhoa2=MAX(grhoa2,small)
                ELSE
                   grhoa2=small
                ENDIF
                IF (rhob1.GT.small) THEN
                   grhob1=grad(i1,i2+0,i3,5)
                   grhob1=MAX(grhob1,small)
                ELSE
                   grhob1=small
                ENDIF
                IF (rhob2.GT.small) THEN
                   grhob2=grad(i1,i2+1,i3,5)
                   grhob2=MAX(grhob2,small)
                ELSE
                   grhob2=small
                ENDIF

                ! ..exchange
                aa1   = 4._real_8*grhoa1
                aa2   = 4._real_8*grhoa2
                rr1   = t43*rsa1*rsa1*rsa1*rsa1
                rr2   = t43*rsa2*rsa2*rsa2*rsa2
                ex1   = ax/rr1
                ex2   = ax/rr2
                s21   = aa1*rr1*rr1*us*us
                s22   = aa2*rr2*rr2*us*us
                po1   = 1._real_8/(1._real_8 + ul*s21)
                po2   = 1._real_8/(1._real_8 + ul*s22)
                fx1   = uk-uk*po1
                fx2   = uk-uk*po2
                sxa1  = 0.5_real_8*ex1*fx1
                sxa2  = 0.5_real_8*ex2*fx2
                dfx1  = 2._real_8*uk*ul*po1*po1
                dfx2  = 2._real_8*uk*ul*po2*po2
                v1xa1 = t13*1.33333333333333_real_8*ax*x1*(fx1-s21*dfx1)
                v1xa2 = t13*1.33333333333333_real_8*ax*x2*(fx2-s22*dfx2)
                v2xa1 = 2._real_8*ex1*dfx1*(us*rr1)**2
                v2xa2 = 2._real_8*ex2*dfx2*(us*rr2)**2
                aa1   = 4._real_8*grhob1
                aa2   = 4._real_8*grhob2
                rr1   = t43*rsb1*rsb1*rsb1*rsb1
                rr2   = t43*rsb2*rsb2*rsb2*rsb2
                ex1   = ax/rr1
                ex2   = ax/rr2
                s21   = aa1*rr1*rr1*us*us
                s22   = aa2*rr2*rr2*us*us
                po1   = 1._real_8/(1._real_8 + ul*s21)
                po2   = 1._real_8/(1._real_8 + ul*s22)
                fx1   = uk-uk*po1
                fx2   = uk-uk*po2
                sxb1  = 0.5_real_8*ex1*fx1
                sxb2  = 0.5_real_8*ex2*fx2
                dfx1  = 2._real_8*uk*ul*po1*po1
                dfx2  = 2._real_8*uk*ul*po2*po2
                v1xb1 = t13*1.33333333333333_real_8*ax*y1*(fx1-s21*dfx1)
                v1xb2 = t13*1.33333333333333_real_8*ax*y2*(fx2-s22*dfx2)
                v2xb1 = 2._real_8*ex1*dfx1*(us*rr1)**2
                v2xb2 = 2._real_8*ex2*dfx2*(us*rr2)**2
                ! ..LDA correlation
                eta1=(rho1-2._real_8*rhob1)/rho1
                eta2=(rho2-2._real_8*rhob2)/rho2
                et1=eta1*eta1
                fxs1=et1*(t1+et1*(t2+et1*(t3+et1*(t4+et1*(t5+et1*&
                     (t6+et1*(t7+et1*(t8+et1*t9))))))))
                dfxs1=(2._real_8*eta1)*(t1+et1*(2._real_8*t2+et1*(3._real_8*t3+et1*&
                     (4._real_8*t4+et1*(5._real_8*t5+et1*(6._real_8*t6+et1*(7._real_8*t7+et1*&
                     (8._real_8*t8+et1*9._real_8*t9))))))))
                et2=eta2*eta2
                fxs2=et2*(t1+et2*(t2+et2*(t3+et2*(t4+et2*(t5+et2*&
                     (t6+et2*(t7+et2*(t8+et2*t9))))))))
                dfxs2=(2._real_8*eta2)*(t1+et2*(2._real_8*t2+et2*(3._real_8*t3+et2*&
                     (4._real_8*t4+et2*(5._real_8*t5+et2*(6._real_8*t6+et2*(7._real_8*t7+et2*&
                     (8._real_8*t8+et2*9._real_8*t9))))))))
                ! 
                a0p1=a0u+fxs1*da0
                a1p1=a1u+fxs1*da1
                a2p1=a2u+fxs1*da2
                a3p1=a3u+fxs1*da3
                b1p1=b1u+fxs1*db1
                b2p1=b2u+fxs1*db2
                b3p1=b3u+fxs1*db3
                b4p1=b4u+fxs1*db4
                ! 
                a0p2=a0u+fxs2*da0
                a1p2=a1u+fxs2*da1
                a2p2=a2u+fxs2*da2
                a3p2=a3u+fxs2*da3
                b1p2=b1u+fxs2*db1
                b2p2=b2u+fxs2*db2
                b3p2=b3u+fxs2*db3
                b4p2=b4u+fxs2*db4
                ! 
                top1=a0p1+rs1*(a1p1+rs1*(a2p1+rs1*a3p1))
                top2=a0p2+rs2*(a1p2+rs2*(a2p2+rs2*a3p2))
                dtop1= a1p1+rs1*(2._real_8*a2p1+rs1*3._real_8*a3p1)
                dtop2= a1p2+rs2*(2._real_8*a2p2+rs2*3._real_8*a3p2)
                xtop1=da0+rs1*(da1+rs1*(da2+rs1*da3))
                xtop2=da0+rs2*(da1+rs2*(da2+rs2*da3))
                bot1=rs1*(b1p1+rs1*(b2p1+rs1*(b3p1+rs1*b4p1)))
                bot2=rs2*(b1p2+rs2*(b2p2+rs2*(b3p2+rs2*b4p2)))
                dbot1=b1p1+rs1*(2._real_8*b2p1+rs1*(3._real_8*b3p1+rs1*4._real_8*b4p1))
                dbot2=b1p2+rs2*(2._real_8*b2p2+rs2*(3._real_8*b3p2+rs2*4._real_8*b4p2))
                xbot1=rs1*(db1+rs1*(db2+rs1*(db3+rs1*db4)))
                xbot2=rs2*(db1+rs2*(db2+rs2*(db3+rs2*db4)))
                obot1=1._real_8/bot1
                obot2=1._real_8/bot2
                epsxcu1=-top1*obot1
                epsxcu2=-top2*obot2
                vxc1=epsxcu1+rs1*(dtop1*obot1-dbot1*(top1*obot1*obot1))*&
                     0.333333333333333_real_8
                dvxc1=-(xtop1*obot1-xbot1*(top1*obot1*obot1))
                dxa1=vxc1+dvxc1*dfxs1*(1._real_8-eta1)
                dxb1=vxc1-dvxc1*dfxs1*(1._real_8+eta1)
                vxc2=epsxcu2+rs2*(dtop2*obot2-dbot2*(top2*obot2*obot2))*&
                     0.333333333333333_real_8
                dvxc2=-(xtop2*obot2-xbot2*(top2*obot2*obot2))
                dxa2=vxc2+dvxc2*dfxs2*(1._real_8-eta2)
                dxb2=vxc2-dvxc2*dfxs2*(1._real_8+eta2)
                esla1=f1s*(x1*rhoa1+y1*rhob1)/rho1
                esla2=f1s*(x2*rhoa2+y2*rhob2)/rho2
                vsla1=f1u*x1
                vslb1=f1u*y1
                vsla2=f1u*x2
                vslb2=f1u*y2
                epsxcu1=epsxcu1-esla1
                epsxcu2=epsxcu2-esla2
                dxa1=dxa1-vsla1
                dxb1=dxb1-vslb1
                dxa2=dxa2-vsla2
                dxb2=dxb2-vslb2
                ! ..correlation
                grho1=grad(i1,i2+0,i3,1)+grad(i1,i2+0,i3,5)+&
                     2._real_8*grad(i1,i2+0,i3,2)*grad(i1,i2+0,i3,6)+&
                     2._real_8*grad(i1,i2+0,i3,3)*grad(i1,i2+0,i3,7)+&
                     2._real_8*grad(i1,i2+0,i3,4)*grad(i1,i2+0,i3,8)
                grho2=grad(i1,i2+1,i3,1)+grad(i1,i2+1,i3,5)+&
                     2._real_8*grad(i1,i2+1,i3,2)*grad(i1,i2+1,i3,6)+&
                     2._real_8*grad(i1,i2+1,i3,3)*grad(i1,i2+1,i3,7)+&
                     2._real_8*grad(i1,i2+1,i3,4)*grad(i1,i2+1,i3,8)
                as1=SQRT(grho1)
                as2=SQRT(grho2)
                grho1=MAX(grho1,small)
                grho2=MAX(grho2,small)
                phi1=0.5_real_8*((1._real_8+eta1)**(two3)+(1._real_8-eta1)**(two3))
                phi31=phi1*phi1*phi1
                phi2=0.5_real_8*((1._real_8+eta2)**(two3)+(1._real_8-eta2)**(two3))
                phi32=phi2*phi2*phi2
                rs1   = (3._real_8/(4._real_8*pi*rho1))**ob3
                rs2   = (3._real_8/(4._real_8*pi*rho2))**ob3
                xkf1  = (9._real_8*pi/4._real_8)**ob3/rs1
                xkf2  = (9._real_8*pi/4._real_8)**ob3/rs2
                xks1  = SQRT(4._real_8*xkf1/pi)
                xks2  = SQRT(4._real_8*xkf2/pi)
                tt1   = as1/(2._real_8*xks1*rho1*phi1)
                tt2   = as2/(2._real_8*xks2*rho2*phi2)
                expe1 = EXP(-epsxcu1/(phi31*ga))
                expe2 = EXP(-epsxcu2/(phi32*ga))
                af1   = be/ga * (1._real_8/(expe1-1._real_8))
                af2   = be/ga * (1._real_8/(expe2-1._real_8))
                yt1   = af1*tt1*tt1
                yt2   = af2*tt2*tt2
                xy1   = (1._real_8+yt1)/(1._real_8+yt1+yt1*yt1)
                xy2   = (1._real_8+yt2)/(1._real_8+yt2+yt2*yt2)
                s11   = 1._real_8+be/ga*tt1*tt1*xy1
                s12   = 1._real_8+be/ga*tt2*tt2*xy2
                h01   = ga*phi31 * LOG(s11)
                h02   = ga*phi32 * LOG(s12)
                sc1   = rho1*h01
                sc2   = rho2*h02
                ! 
                dphide1= 1._real_8/(3._real_8*(1._real_8+eta1)**ob3+ssmall)-&
                     1._real_8/(3._real_8*(1._real_8-eta1)**ob3+ssmall)
                dphide2= 1._real_8/(3._real_8*(1._real_8+eta2)**ob3+ssmall)-&
                     1._real_8/(3._real_8*(1._real_8-eta2)**ob3+ssmall)
                dedra1 = 2._real_8*rhob1/(rho1*rho1)
                dedrb1 = -2._real_8*rhoa1/(rho1*rho1)
                dedra2 = 2._real_8*rhob2/(rho2*rho2)
                dedrb2 = -2._real_8*rhoa2/(rho2*rho2)
                dphida1= dphide1*dedra1
                dphidb1= dphide1*dedrb1
                dphida2= dphide2*dedra2
                dphidb2= dphide2*dedrb2
                dtdra1 = -tt1*(dphida1/phi1+7._real_8/(6._real_8*rho1))
                dtdrb1 = -tt1*(dphidb1/phi1+7._real_8/(6._real_8*rho1))
                dtdra2 = -tt2*(dphida2/phi2+7._real_8/(6._real_8*rho2))
                dtdrb2 = -tt2*(dphidb2/phi2+7._real_8/(6._real_8*rho2))
                dadra1 = af1*af1*expe1/(-be*phi31)*(3._real_8*epsxcu1/phi1*&
                     dphida1-(dxa1-epsxcu1)/rho1)
                dadrb1 = af1*af1*expe1/(-be*phi31)*(3._real_8*epsxcu1/phi1*&
                     dphidb1-(dxb1-epsxcu1)/rho1)
                dadra2 = af2*af2*expe2/(-be*phi32)*(3._real_8*epsxcu2/phi2*&
                     dphida2-(dxa2-epsxcu2)/rho2)
                dadrb2 = af2*af2*expe2/(-be*phi32)*(3._real_8*epsxcu2/phi2*&
                     dphidb2-(dxb2-epsxcu2)/rho2)
                dsda1  = -be/ga*af1*tt1**6*(2._real_8+yt1)&
                     /(1._real_8+yt1+yt1*yt1)**2
                dsda2  = -be/ga*af2*tt2**6*(2._real_8+yt2)&
                     /(1._real_8+yt2+yt2*yt2)**2
                dsdt1  = 2._real_8*be/ga*tt1*(1._real_8+2._real_8*yt1)&
                     /(1._real_8+yt1+yt1*yt1)**2
                dsdt2  = 2._real_8*be/ga*tt2*(1._real_8+2._real_8*yt2)&
                     /(1._real_8+yt2+yt2*yt2)**2
                dsdra1 = dsda1*dadra1 + dsdt1*dtdra1
                dsdrb1 = dsda1*dadrb1 + dsdt1*dtdrb1
                dsdra2 = dsda2*dadra2 + dsdt2*dtdra2
                dsdrb2 = dsda2*dadrb2 + dsdt2*dtdrb2
                dhdt1  = ga*phi31/s11*dsdt1
                dhdt2  = ga*phi32/s12*dsdt2
                dhdra1 = 3._real_8*h01/phi1*dphida1 + ga*phi31/s11*dsdra1
                dhdrb1 = 3._real_8*h01/phi1*dphidb1 + ga*phi31/s11*dsdrb1
                dhdra2 = 3._real_8*h02/phi2*dphida2 + ga*phi32/s12*dsdra2
                dhdrb2 = 3._real_8*h02/phi2*dphidb2 + ga*phi32/s12*dsdrb2
                ! ..
                v1ca1 = h01 + rho1*dhdra1
                v2ca1 = rho1*dhdt1*tt1/grho1
                v1cb1 = h01 + rho1*dhdrb1
                v2cb1 = rho1*dhdt1*tt1/grho1
                v2cab1= rho1*dhdt1*tt1/grho1
                v1ca2 = h02 + rho2*dhdra2
                v2ca2 = rho2*dhdt2*tt2/grho2
                v1cb2 = h02 + rho2*dhdrb2
                v2cb2 = rho2*dhdt2*tt2/grho2
                v2cab2= rho2*dhdt2*tt2/grho2
                ! ..sum up terms
                sgcx = sgcx + (smoo1*(sxa1+sxb1) + smoo2*(sxa2+sxb2))
                sgcc = sgcc + (smoo1*sc1 + smoo2*sc2)

                IF (rhoa1.GT.small) THEN
                   v1(i1,i2+0,i3)=dsmoo1*(sxa1+sc1)+smoo1*(v1xa1+v1ca1)
                   vtmp(i1,i2+0,i3,1)=smoo1*(v2xa1+v2ca1)
                ENDIF
                IF (rhoa2.GT.small) THEN
                   v1(i1,i2+1,i3)=dsmoo2*(sxa2+sc2)+smoo2*(v1xa2+v1ca2)
                   vtmp(i1,i2+1,i3,1)=smoo2*(v2xa2+v2ca2)
                ENDIF
                IF (rhob1.GT.small) THEN
                   v2(i1,i2+0,i3)=dsmoo1*(sxb1+sc1)+smoo1*(v1xb1+v1cb1)
                   vtmp(i1,i2+0,i3,2)=smoo1*(v2xb1+v2cb1)
                ENDIF
                IF (rhob2.GT.small) THEN
                   v2(i1,i2+1,i3)=dsmoo2*(sxb2+sc2)+smoo2*(v1xb2+v1cb2)
                   vtmp(i1,i2+1,i3,2)=smoo2*(v2xb2+v2cb2)
                ENDIF
                IF ((rhoa1.GT.small).AND.(rhob1.GT.small)) THEN
                   vtmp(i1,i2+0,i3,3)=smoo1*v2cab1
                ENDIF
                IF ((rhoa2.GT.small).AND.(rhob2.GT.small)) THEN
                   vtmp(i1,i2+1,i3,3)=smoo2*v2cab2
                ENDIF
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    ! ..counting floating point operations
    fadd   = ic *  6._real_8 + id * 195._real_8
    fmult  = ic *  6._real_8 + id * 445._real_8
    fdiv   = ic *  2._real_8 + id *  76._real_8
    fspez  =              id *  18._real_8
    flops=fadd+fmult+fdiv+fspez
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gcspbe
  ! ==================================================================
  SUBROUTINE gcsrevpbe(sgcx,sgcc,rhoe,v1,v2,vtmp,grad,flops)
    REAL(real_8)                             :: sgcx, sgcc
    REAL(real_8), DIMENSION(fpar%kr1, fpar%&
      kr2s, fpar%kr3s, 2)                    :: rhoe
    COMPLEX(real_8) :: v1(fpar%kr1,fpar%kr2s,fpar%kr3s), &
      v2(fpar%kr1,fpar%kr2s,fpar%kr3s)
    REAL(real_8), DIMENSION(fpar%kr1, fpar%&
      kr2s, fpar%kr3s, 3)                    :: vtmp
    REAL(real_8), DIMENSION(fpar%kr1, fpar%&
      kr2s, fpar%kr3s, 8)                    :: grad
    REAL(real_8)                             :: flops

    REAL(real_8), PARAMETER :: a0u = 0.4581652932831429_real_8, &
      a1u = 2.217058676663745_real_8, a2u = 0.7405551735357053_real_8, &
      a3u = 0.01968227878617998_real_8 , ax = -0.738558766382022406_real_8, &
      b1u = 1.0_real_8, b2u = 4.504130959426697_real_8, &
      b3u = 1.110667363742916_real_8, b4u = 0.02359291751427506_real_8 , &
      be = 0.06672455060314922_real_8, da0 = 0.119086804055547_real_8, &
      da1 = 0.6157402568883345_real_8, da2 = 0.1574201515892867_real_8, &
      da3 = 0.003532336663397157_real_8 , db1 = 0.0_real_8, &
      db2 = 0.2673612973836267_real_8, db3 = 0.2052004607777787_real_8, &
      db4 = 0.004200005045691381_real_8 
    REAL(real_8), PARAMETER :: eps2 = 1.e-20_real_8, &
      f1 = -1.39578858466194911_real_8 , ga = 0.031090690869654895_real_8 , &
      rsfac = 0.6203504908994000_real_8 , small = 1.e-24_real_8, &
      ssmall = 1.e-12_real_8, t1 = 0.854960467080683_real_8, &
      t2 = 0.07916300621117439_real_8, t3 = 0.02580127609845684_real_8, &
      t4 = 0.0121839359353824_real_8, t5 = 0.006919272259599879_real_8, &
      t6 = 0.004391524649611372_real_8, t7 = 0.003002751897170168_real_8, &
      t8 = 0.00216587382212552_real_8, t9 = 0.01141189204579585_real_8 , &
      two3 = 2._real_8/3._real_8 , uk = 1.2450_real_8, &
      um = 0.2195149727645171_real_8, ul = um/uk 
    REAL(real_8), PARAMETER :: us = 0.161620459673995492_real_8

    INTEGER                                  :: i1, i2, i3, ic, id, ii2, &
                                                ncount
    REAL(real_8) :: a0p1, a0p2, a1p1, a1p2, a2p1, a2p2, a3p1, a3p2, aa1, aa2, &
      af1, af2, as1, as2, b1p1, b1p2, b2p1, b2p2, b3p1, b3p2, b4p1, b4p2, &
      bot1, bot2, d1, d2, dadra1, dadra2, dadrb1, dadrb2, dbot1, dbot2, &
      dedra1, dedra2, dedrb1, dedrb2, dfx1, dfx2, dfxs1, dfxs2, dhdra1, &
      dhdra2, dhdrb1, dhdrb2, dhdt1, dhdt2, dphida1, dphida2, dphidb1, &
      dphidb2, dphide1, dphide2, drho, dsda1, dsda2, dsdra1, dsdra2, dsdrb1, &
      dsdrb2, dsdt1, dsdt2, dsmoo1, dsmoo2, dtdra1, dtdra2, dtdrb1, dtdrb2, &
      dtop1, dtop2, dvxc1, dvxc2, dxa1, dxa2, dxb1, dxb2, epsxcu1, epsxcu2, &
      esla1, esla2, et1, et2, eta1, eta2, ex1, ex2
    REAL(real_8) :: expe1, expe2, f1s, f1u, fadd, fdiv, fmult, fspez, fx1, &
      fx2, fxs1, fxs2, grho1, grho2, grhoa1, grhoa2, grhob1, grhob2, h01, &
      h02, ob3, obot1, obot2, ogceps, phi1, phi2, phi31, phi32, po1, po2, &
      rho1, rho2, rhoa1, rhoa2, rhob1, rhob2, rr1, rr2, rs1, rs2, rsa1, rsa2, &
      rsb1, rsb2, s11, s12, s21, s22, sc1, sc2, smoo1, smoo2, sxa1, sxa2, &
      sxb1, sxb2, t13, t43, top1, top2, tt1, tt2, v1ca1, v1ca2, v1cb1, v1cb2, &
      v1xa1, v1xa2, v1xb1, v1xb2, v2ca1, v2ca2, v2cab1, v2cab2, v2cb1, v2cb2, &
      v2xa1, v2xa2, v2xb1, v2xb2, vsla1, vsla2, vslb1, vslb2, vxc1, vxc2, x1, &
      x2, xbot1, xbot2, xkf1, xkf2, xks1
    REAL(real_8) :: xks2, xtop1, xtop2, xy1, xy2, y1, y2, yt1, yt2, z1, z2

! ==--------------------------------------------------------------==

    t43=2._real_8**(-4._real_8/3._real_8)
    t13=2._real_8**(1._real_8/3._real_8)
    ob3=1._real_8/3._real_8
    f1s=f1*2._real_8/3._real_8
    f1u=f1s*1.3333333333333333_real_8
    ogceps=1._real_8/cntr%gceps
    sgcx=0._real_8
    sgcc=0._real_8
    x1=1._real_8
    x2=1._real_8
    y1=1._real_8
    y2=1._real_8
    z1=1._real_8
    z2=1._real_8
    ic=0
    id=0
    ii2=MOD(parm%nr2,2)
    IF (ii2.NE.0) CALL stopgm('GCSrevPBE','ODD DIMENSIONS',& 
         __LINE__,__FILE__)
    !$omp parallel do private(I1,I2,I3,RHOA1,RHOA2,D1,D2) &
    !$omp  firstprivate(X1,X2,Y1,Y2,Z1,Z2) &
    !$omp  private(NCOUNT, RHOB1, RHOB2, RHO1, RHO2, RSA1, RSA2) &
    !$omp  private(RSB1, RSB2, RS1, RS2, SMOO1, DSMOO1, DRHO) &
    !$omp  private(SMOO2, DSMOO2, GRHOA1, GRHOA2, GRHOB1, GRHOB2, AA1) &
    !$omp  private(AA2, RR1, RR2, EX1, EX2, S21, S22, PO1, PO2, FX1, FX2) &
    !$omp  private(SXA1, SXA2, DFX1, DFX2, V1XA1, V1XA2, V2XA1, V2XA2,SXB1) &
    !$omp  private(SXB2, V1XB1, V1XB2, V2XB1, V2XB2, ETA1, ETA2, ET1, FXS1) &
    !$omp  private(DFXS1, ET2, FXS2, DFXS2, A0P1, A1P1, A2P1, A3P1, B2P1) &
    !$omp  private(B3P1, B4P1, A0P2, A1P2, A2P2, A3P2, B2P2, B3P2, B4P2) &
    !$omp  private(TOP1, TOP2, DTOP1, DTOP2, XTOP1, XTOP2, BOT1, BOT2) &
    !$omp  private(DBOT1, DBOT2, XBOT1, XBOT2, OBOT1,OBOT2,EPSXCU1,EPSXCU2) &
    !$omp  private(VXC1, DVXC1, DXA1, DXB1, VXC2, DVXC2, DXA2, DXB2, ESLA1) &
    !$omp  private(ESLA2, VSLA1, VSLB1, VSLA2, VSLB2, GRHO1,GRHO2,AS1,AS2) &
    !$omp  private(PHI1, PHI31, PHI2, PHI32, XKF1, XKF2,XKS1,XKS2,TT1,TT2) &
    !$omp  private(EXPE1, EXPE2, AF1, AF2, YT1, YT2, XY1, XY2,S11,S12,H01) &
    !$omp  private(H02, SC1, SC2, DPHIDE1, DPHIDE2, DEDRA1, DEDRB1, DEDRA2) &
    !$omp  private(DEDRB2, DPHIDA1, DPHIDB1, DPHIDA2,DPHIDB2,DTDRA1,DTDRB1) &
    !$omp  private(DTDRA2, DTDRB2, DADRA1,DADRB1,DADRA2,DADRB2,DSDA1,DSDA2) &
    !$omp  private(DSDT1, DSDT2, DSDRA1, DSDRB1, DSDRA2,DSDRB2,DHDT1,DHDT2) &
    !$omp  private(DHDRA1, DHDRB1, DHDRA2, DHDRB2, V1CA1,V2CA1,V1CB1,V2CB1) &
    !$omp  private(V2CAB1, V1CA2, V2CA2, V1CB2, V2CB2, V2CAB2) &
    !$omp  private(B1P1,B1P2) &
    !$omp  reduction(+:SGCX,SGCC,IC,ID) __COLLAPSE3
    DO i3=1,parm%nr3
       DO i2=1,parm%nr2,2
          DO i1=1,parm%nr1
             rhoa1=MAX(rhoe(i1,i2+0,i3,1),small)
             rhoa2=MAX(rhoe(i1,i2+1,i3,1),small)
             rhob1=MAX(rhoe(i1,i2+0,i3,2),small)
             rhob2=MAX(rhoe(i1,i2+1,i3,2),small)
             rho1=rhoa1+rhob1
             rho2=rhoa2+rhob2
             IF ((rho1+rho2).GT.5._real_8*small) THEN
                id=id+1
                DO ncount=1,100
                   d1=x1-rhoa1/x1**2
                   d2=x2-rhoa2/x2**2
                   ic=ic+1
                   x1=x1-0.333333333333333_real_8*d1
                   x2=x2-0.333333333333333_real_8*d2
                   IF (d1**2+d2**2.LT.eps2) GOTO 10
                ENDDO
                x1=rhoa1**(1._real_8/3._real_8)
                x2=rhoa2**(1._real_8/3._real_8)
10              CONTINUE
                rsa1=1._real_8/x1
                rsa2=1._real_8/x2
                DO ncount=1,100
                   d1=y1-rhob1/y1**2
                   d2=y2-rhob2/y2**2
                   ic=ic+1
                   y1=y1-0.333333333333333_real_8*d1
                   y2=y2-0.333333333333333_real_8*d2
                   IF (d1**2+d2**2.LT.eps2) GOTO 20
                ENDDO
                y1=rhob1**(1._real_8/3._real_8)
                y2=rhob2**(1._real_8/3._real_8)
20              CONTINUE
                rsb1=1._real_8/y1
                rsb2=1._real_8/y2
                DO ncount=1,100
                   d1=z1-rho1/z1**2
                   d2=z2-rho2/z2**2
                   ic=ic+1
                   z1=z1-0.333333333333333_real_8*d1
                   z2=z2-0.333333333333333_real_8*d2
                   IF (d1**2+d2**2.LT.eps2) GOTO 30
                ENDDO
                z1=rho1**(1._real_8/3._real_8)
                z2=rho2**(1._real_8/3._real_8)
30              CONTINUE
                rs1=rsfac/z1
                rs2=rsfac/z2
                IF (rho1.GT.2.25_real_8*cntr%gceps) THEN
                   smoo1=1._real_8
                   dsmoo1=0._real_8
                ELSEIF (rho1.LT.0.25_real_8*cntr%gceps) THEN
                   smoo1=0._real_8
                   dsmoo1=0._real_8
                ELSE
                   drho=(rho1-0.25_real_8*cntr%gceps)*ogceps
                   smoo1=drho*drho*(0.75_real_8-0.25_real_8*drho)
                   dsmoo1=1.5_real_8*drho*ogceps*(1._real_8-0.5_real_8*drho)
                ENDIF
                IF (rho2.GT.2.25_real_8*cntr%gceps) THEN
                   smoo2=1._real_8
                   dsmoo2=0._real_8
                ELSEIF (rho2.LT.0.25_real_8*cntr%gceps) THEN
                   smoo2=0._real_8
                   dsmoo2=0._real_8
                ELSE
                   drho=(rho2-0.25_real_8*cntr%gceps)*ogceps
                   smoo2=drho*drho*(0.75_real_8-0.25_real_8*drho)
                   dsmoo2=1.5_real_8*drho*ogceps*(1._real_8-0.5_real_8*drho)
                ENDIF

                IF (rhoa1.GT.small) THEN
                   grhoa1=grad(i1,i2+0,i3,1)
                   grhoa1=MAX(grhoa1,small)
                ELSE
                   grhoa1=small
                ENDIF
                IF (rhoa2.GT.small) THEN
                   grhoa2=grad(i1,i2+1,i3,1)
                   grhoa2=MAX(grhoa2,small)
                ELSE
                   grhoa2=small
                ENDIF
                IF (rhob1.GT.small) THEN
                   grhob1=grad(i1,i2+0,i3,5)
                   grhob1=MAX(grhob1,small)
                ELSE
                   grhob1=small
                ENDIF
                IF (rhob2.GT.small) THEN
                   grhob2=grad(i1,i2+1,i3,5)
                   grhob2=MAX(grhob2,small)
                ELSE
                   grhob2=small
                ENDIF

                ! ..exchange
                aa1   = 4._real_8*grhoa1
                aa2   = 4._real_8*grhoa2
                rr1   = t43*rsa1*rsa1*rsa1*rsa1
                rr2   = t43*rsa2*rsa2*rsa2*rsa2
                ex1   = ax/rr1
                ex2   = ax/rr2
                s21   = aa1*rr1*rr1*us*us
                s22   = aa2*rr2*rr2*us*us
                po1   = 1._real_8/(1._real_8 + ul*s21)
                po2   = 1._real_8/(1._real_8 + ul*s22)
                fx1   = uk-uk*po1
                fx2   = uk-uk*po2
                sxa1  = 0.5_real_8*ex1*fx1
                sxa2  = 0.5_real_8*ex2*fx2
                dfx1  = 2._real_8*uk*ul*po1*po1
                dfx2  = 2._real_8*uk*ul*po2*po2
                v1xa1 = t13*1.33333333333333_real_8*ax*x1*(fx1-s21*dfx1)
                v1xa2 = t13*1.33333333333333_real_8*ax*x2*(fx2-s22*dfx2)
                v2xa1 = 2._real_8*ex1*dfx1*(us*rr1)**2
                v2xa2 = 2._real_8*ex2*dfx2*(us*rr2)**2
                aa1   = 4._real_8*grhob1
                aa2   = 4._real_8*grhob2
                rr1   = t43*rsb1*rsb1*rsb1*rsb1
                rr2   = t43*rsb2*rsb2*rsb2*rsb2
                ex1   = ax/rr1
                ex2   = ax/rr2
                s21   = aa1*rr1*rr1*us*us
                s22   = aa2*rr2*rr2*us*us
                po1   = 1._real_8/(1._real_8 + ul*s21)
                po2   = 1._real_8/(1._real_8 + ul*s22)
                fx1   = uk-uk*po1
                fx2   = uk-uk*po2
                sxb1  = 0.5_real_8*ex1*fx1
                sxb2  = 0.5_real_8*ex2*fx2
                dfx1  = 2._real_8*uk*ul*po1*po1
                dfx2  = 2._real_8*uk*ul*po2*po2
                v1xb1 = t13*1.33333333333333_real_8*ax*y1*(fx1-s21*dfx1)
                v1xb2 = t13*1.33333333333333_real_8*ax*y2*(fx2-s22*dfx2)
                v2xb1 = 2._real_8*ex1*dfx1*(us*rr1)**2
                v2xb2 = 2._real_8*ex2*dfx2*(us*rr2)**2
                ! ..LDA correlation
                eta1=(rho1-2._real_8*rhob1)/rho1
                eta2=(rho2-2._real_8*rhob2)/rho2
                et1=eta1*eta1
                fxs1=et1*(t1+et1*(t2+et1*(t3+et1*(t4+et1*(t5+et1*&
                     (t6+et1*(t7+et1*(t8+et1*t9))))))))
                dfxs1=(2._real_8*eta1)*(t1+et1*(2._real_8*t2+et1*(3._real_8*t3+et1*&
                     (4._real_8*t4+et1*(5._real_8*t5+et1*(6._real_8*t6+et1*(7._real_8*t7+et1*&
                     (8._real_8*t8+et1*9._real_8*t9))))))))
                et2=eta2*eta2
                fxs2=et2*(t1+et2*(t2+et2*(t3+et2*(t4+et2*(t5+et2*&
                     (t6+et2*(t7+et2*(t8+et2*t9))))))))
                dfxs2=(2._real_8*eta2)*(t1+et2*(2._real_8*t2+et2*(3._real_8*t3+et2*&
                     (4._real_8*t4+et2*(5._real_8*t5+et2*(6._real_8*t6+et2*(7._real_8*t7+et2*&
                     (8._real_8*t8+et2*9._real_8*t9))))))))
                ! 
                a0p1=a0u+fxs1*da0
                a1p1=a1u+fxs1*da1
                a2p1=a2u+fxs1*da2
                a3p1=a3u+fxs1*da3
                b1p1=b1u+fxs1*db1
                b2p1=b2u+fxs1*db2
                b3p1=b3u+fxs1*db3
                b4p1=b4u+fxs1*db4
                ! 
                a0p2=a0u+fxs2*da0
                a1p2=a1u+fxs2*da1
                a2p2=a2u+fxs2*da2
                a3p2=a3u+fxs2*da3
                b1p2=b1u+fxs2*db1
                b2p2=b2u+fxs2*db2
                b3p2=b3u+fxs2*db3
                b4p2=b4u+fxs2*db4
                ! 
                top1=a0p1+rs1*(a1p1+rs1*(a2p1+rs1*a3p1))
                top2=a0p2+rs2*(a1p2+rs2*(a2p2+rs2*a3p2))
                dtop1= a1p1+rs1*(2._real_8*a2p1+rs1*3._real_8*a3p1)
                dtop2= a1p2+rs2*(2._real_8*a2p2+rs2*3._real_8*a3p2)
                xtop1=da0+rs1*(da1+rs1*(da2+rs1*da3))
                xtop2=da0+rs2*(da1+rs2*(da2+rs2*da3))
                bot1=rs1*(b1p1+rs1*(b2p1+rs1*(b3p1+rs1*b4p1)))
                bot2=rs2*(b1p2+rs2*(b2p2+rs2*(b3p2+rs2*b4p2)))
                dbot1=b1p1+rs1*(2._real_8*b2p1+rs1*(3._real_8*b3p1+rs1*4._real_8*b4p1))
                dbot2=b1p2+rs2*(2._real_8*b2p2+rs2*(3._real_8*b3p2+rs2*4._real_8*b4p2))
                xbot1=rs1*(db1+rs1*(db2+rs1*(db3+rs1*db4)))
                xbot2=rs2*(db1+rs2*(db2+rs2*(db3+rs2*db4)))
                obot1=1._real_8/bot1
                obot2=1._real_8/bot2
                epsxcu1=-top1*obot1
                epsxcu2=-top2*obot2
                vxc1=epsxcu1+rs1*(dtop1*obot1-dbot1*(top1*obot1*obot1))*&
                     0.333333333333333_real_8
                dvxc1=-(xtop1*obot1-xbot1*(top1*obot1*obot1))
                dxa1=vxc1+dvxc1*dfxs1*(1._real_8-eta1)
                dxb1=vxc1-dvxc1*dfxs1*(1._real_8+eta1)
                vxc2=epsxcu2+rs2*(dtop2*obot2-dbot2*(top2*obot2*obot2))*&
                     0.333333333333333_real_8
                dvxc2=-(xtop2*obot2-xbot2*(top2*obot2*obot2))
                dxa2=vxc2+dvxc2*dfxs2*(1._real_8-eta2)
                dxb2=vxc2-dvxc2*dfxs2*(1._real_8+eta2)
                esla1=f1s*(x1*rhoa1+y1*rhob1)/rho1
                esla2=f1s*(x2*rhoa2+y2*rhob2)/rho2
                vsla1=f1u*x1
                vslb1=f1u*y1
                vsla2=f1u*x2
                vslb2=f1u*y2
                epsxcu1=epsxcu1-esla1
                epsxcu2=epsxcu2-esla2
                dxa1=dxa1-vsla1
                dxb1=dxb1-vslb1
                dxa2=dxa2-vsla2
                dxb2=dxb2-vslb2
                ! ..correlation
                grho1=grad(i1,i2+0,i3,1)+grad(i1,i2+0,i3,5)+&
                     2._real_8*grad(i1,i2+0,i3,2)*grad(i1,i2+0,i3,6)+&
                     2._real_8*grad(i1,i2+0,i3,3)*grad(i1,i2+0,i3,7)+&
                     2._real_8*grad(i1,i2+0,i3,4)*grad(i1,i2+0,i3,8)
                grho2=grad(i1,i2+1,i3,1)+grad(i1,i2+1,i3,5)+&
                     2._real_8*grad(i1,i2+1,i3,2)*grad(i1,i2+1,i3,6)+&
                     2._real_8*grad(i1,i2+1,i3,3)*grad(i1,i2+1,i3,7)+&
                     2._real_8*grad(i1,i2+1,i3,4)*grad(i1,i2+1,i3,8)
                as1=SQRT(grho1)
                as2=SQRT(grho2)
                phi1=0.5_real_8*((1._real_8+eta1)**(two3)+(1._real_8-eta1)**(two3))
                phi31=phi1*phi1*phi1
                phi2=0.5_real_8*((1._real_8+eta2)**(two3)+(1._real_8-eta2)**(two3))
                phi32=phi2*phi2*phi2
                rs1   = (3._real_8/(4._real_8*pi*rho1))**ob3
                rs2   = (3._real_8/(4._real_8*pi*rho2))**ob3
                xkf1  = (9._real_8*pi/4._real_8)**ob3/rs1
                xkf2  = (9._real_8*pi/4._real_8)**ob3/rs2
                xks1  = SQRT(4._real_8*xkf1/pi)
                xks2  = SQRT(4._real_8*xkf2/pi)
                tt1   = as1/(2._real_8*xks1*rho1*phi1)
                tt2   = as2/(2._real_8*xks2*rho2*phi2)
                expe1 = EXP(-epsxcu1/(phi31*ga))
                expe2 = EXP(-epsxcu2/(phi32*ga))
                af1   = be/ga * (1._real_8/(expe1-1._real_8))
                af2   = be/ga * (1._real_8/(expe2-1._real_8))
                yt1   = af1*tt1*tt1
                yt2   = af2*tt2*tt2
                xy1   = (1._real_8+yt1)/(1._real_8+yt1+yt1*yt1)
                xy2   = (1._real_8+yt2)/(1._real_8+yt2+yt2*yt2)
                s11   = 1._real_8+be/ga*tt1*tt1*xy1
                s12   = 1._real_8+be/ga*tt2*tt2*xy2
                h01   = ga*phi31 * LOG(s11)
                h02   = ga*phi32 * LOG(s12)
                sc1   = rho1*h01
                sc2   = rho2*h02
                ! 
                dphide1= 1._real_8/(3._real_8*(1._real_8+eta1)**ob3+ssmall)-&
                     1._real_8/(3._real_8*(1._real_8-eta1)**ob3+ssmall)
                dphide2= 1._real_8/(3._real_8*(1._real_8+eta2)**ob3+ssmall)-&
                     1._real_8/(3._real_8*(1._real_8-eta2)**ob3+ssmall)
                dedra1 = 2._real_8*rhob1/(rho1*rho1)
                dedrb1 = -2._real_8*rhoa1/(rho1*rho1)
                dedra2 = 2._real_8*rhob2/(rho2*rho2)
                dedrb2 = -2._real_8*rhoa2/(rho2*rho2)
                dphida1= dphide1*dedra1
                dphidb1= dphide1*dedrb1
                dphida2= dphide2*dedra2
                dphidb2= dphide2*dedrb2
                dtdra1 = -tt1*(dphida1/phi1+7._real_8/(6._real_8*rho1))
                dtdrb1 = -tt1*(dphidb1/phi1+7._real_8/(6._real_8*rho1))
                dtdra2 = -tt2*(dphida2/phi2+7._real_8/(6._real_8*rho2))
                dtdrb2 = -tt2*(dphidb2/phi2+7._real_8/(6._real_8*rho2))
                dadra1 = af1*af1*expe1/(-be*phi31)*(3._real_8*epsxcu1/phi1*&
                     dphida1-(dxa1-epsxcu1)/rho1)
                dadrb1 = af1*af1*expe1/(-be*phi31)*(3._real_8*epsxcu1/phi1*&
                     dphidb1-(dxb1-epsxcu1)/rho1)
                dadra2 = af2*af2*expe2/(-be*phi32)*(3._real_8*epsxcu2/phi2*&
                     dphida2-(dxa2-epsxcu2)/rho2)
                dadrb2 = af2*af2*expe2/(-be*phi32)*(3._real_8*epsxcu2/phi2*&
                     dphidb2-(dxb2-epsxcu2)/rho2)
                dsda1  = -be/ga*af1*tt1**6*(2._real_8+yt1)&
                     /(1._real_8+yt1+yt1*yt1)**2
                dsda2  = -be/ga*af2*tt2**6*(2._real_8+yt2)&
                     /(1._real_8+yt2+yt2*yt2)**2
                dsdt1  = 2._real_8*be/ga*tt1*(1._real_8+2._real_8*yt1)&
                     /(1._real_8+yt1+yt1*yt1)**2
                dsdt2  = 2._real_8*be/ga*tt2*(1._real_8+2._real_8*yt2)&
                     /(1._real_8+yt2+yt2*yt2)**2
                dsdra1 = dsda1*dadra1 + dsdt1*dtdra1
                dsdrb1 = dsda1*dadrb1 + dsdt1*dtdrb1
                dsdra2 = dsda2*dadra2 + dsdt2*dtdra2
                dsdrb2 = dsda2*dadrb2 + dsdt2*dtdrb2
                dhdt1  = ga*phi31/s11*dsdt1
                dhdt2  = ga*phi32/s12*dsdt2
                dhdra1 = 3._real_8*h01/phi1*dphida1 + ga*phi31/s11*dsdra1
                dhdrb1 = 3._real_8*h01/phi1*dphidb1 + ga*phi31/s11*dsdrb1
                dhdra2 = 3._real_8*h02/phi2*dphida2 + ga*phi32/s12*dsdra2
                dhdrb2 = 3._real_8*h02/phi2*dphidb2 + ga*phi32/s12*dsdrb2
                ! ..
                v1ca1 = h01 + rho1*dhdra1
                v2ca1 = rho1*dhdt1*tt1/grho1
                v1cb1 = h01 + rho1*dhdrb1
                v2cb1 = rho1*dhdt1*tt1/grho1
                v2cab1= rho1*dhdt1*tt1/grho1
                v1ca2 = h02 + rho2*dhdra2
                v2ca2 = rho2*dhdt2*tt2/grho2
                v1cb2 = h02 + rho2*dhdrb2
                v2cb2 = rho2*dhdt2*tt2/grho2
                v2cab2= rho2*dhdt2*tt2/grho2
                ! ..sum up terms
                sgcx = sgcx + (smoo1*(sxa1+sxb1) + smoo2*(sxa2+sxb2))
                sgcc = sgcc + (smoo1*sc1 + smoo2*sc2)

                IF (rhoa1.GT.small) THEN
                   v1(i1,i2+0,i3)=dsmoo1*(sxa1+sc1)+smoo1*(v1xa1+v1ca1)
                   vtmp(i1,i2+0,i3,1)=smoo1*(v2xa1+v2ca1)
                ENDIF
                IF (rhoa2.GT.small) THEN
                   v1(i1,i2+1,i3)=dsmoo2*(sxa2+sc2)+smoo2*(v1xa2+v1ca2)
                   vtmp(i1,i2+1,i3,1)=smoo2*(v2xa2+v2ca2)
                ENDIF
                IF (rhob1.GT.small) THEN
                   v2(i1,i2+0,i3)=dsmoo1*(sxb1+sc1)+smoo1*(v1xb1+v1cb1)
                   vtmp(i1,i2+0,i3,2)=smoo1*(v2xb1+v2cb1)
                ENDIF
                IF (rhob2.GT.small) THEN
                   v2(i1,i2+1,i3)=dsmoo2*(sxb2+sc2)+smoo2*(v1xb2+v1cb2)
                   vtmp(i1,i2+1,i3,2)=smoo2*(v2xb2+v2cb2)
                ENDIF
                IF ((rhoa1.GT.small).AND.(rhob1.GT.small)) THEN
                   vtmp(i1,i2+0,i3,3)=smoo1*v2cab1
                ENDIF
                IF ((rhoa2.GT.small).AND.(rhob2.GT.small)) THEN
                   vtmp(i1,i2+1,i3,3)=smoo2*v2cab2
                ENDIF

             ENDIF
          ENDDO
       ENDDO
    ENDDO
    ! ..counting floating point operations
    fadd   = ic *  6._real_8 + id * 195._real_8
    fmult  = ic *  6._real_8 + id * 445._real_8
    fdiv   = ic *  2._real_8 + id *  76._real_8
    fspez  =              id *  18._real_8
    flops=fadd+fmult+fdiv+fspez
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gcsrevpbe
  ! ==================================================================
#else 
  ! defined(__VECTOR)
  !!cmb - revised Tokyo/Strasbourg, 3/3/2014
  ! ==================================================================
  SUBROUTINE gcxonly(b1,sgcx,sgcc,rhoe,v,vtmp,grad,flops)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: b1, sgcx, sgcc
    REAL(real_8), DIMENSION(fpar%kr1*fpar%&
      kr2s*fpar%kr3s)                        :: rhoe
    COMPLEX(real_8) :: v(fpar%kr1*fpar%kr2s*fpar%kr3s)
    REAL(real_8), DIMENSION(fpar%kr1*fpar%&
      kr2s*fpar%kr3s)                        :: vtmp
    REAL(real_8), DIMENSION(fpar%kr1*fpar%&
      kr2s*fpar%kr3s, *)                     :: grad
    REAL(real_8)                             :: flops

    REAL(real_8), PARAMETER                  :: eps2 = 1.e-20_real_8, &
                                                small = 1.e-24_real_8 

    INTEGER                                  :: i, id, kkk
    REAL(real_8) :: dd, dd2on, ddon, drho, dsmoo, ee, eex, fadd, fdiv, fmult, &
      fspez, grho, ogceps, rho43, rho43on, rhou, sa2b8, shm1, smoo, two13, &
      two13on, vexc1, vexc2, x, xx, yy

! ==--------------------------------------------------------------==

    two13=2._real_8**0.333333333333333_real_8
    two13on=1.0_real_8/two13
    sgcx = 0.0_real_8
    sgcc = 0.0_real_8
    ogceps=1._real_8/cntr%gceps
    id=0
    kkk=fpar%kr1*fpar%kr2s*fpar%kr3s
    !$omp parallel do private(I) &
    !$omp  private(RHOU,X,RHO43,RHO43ON,SMOO,DSMOO,DRHO,GRHO,XX,YY) &
    !$omp  private(SA2B8,SHM1,DD,DDON,DD2ON,EEX,EE,VEXC1,VEXC2) &
    !$omp  reduction(+:SGCX,ID)
    DO i=1,kkk
       rhou=MAX(rhoe(i),small)
       IF (rhou.GT.small) THEN
          x=rhou**(0.33333333333333_real_8)
          rho43=x*rhou
          rho43on=1.0_real_8/rho43
          IF (rhou.GT.2.25_real_8*cntr%gceps) THEN
             smoo=1._real_8
             dsmoo=0._real_8
          ELSEIF (rhou.LT.0.25_real_8*cntr%gceps) THEN
             smoo=0._real_8
             dsmoo=0._real_8
          ELSE
             drho=(rhou-0.25_real_8*cntr%gceps)*ogceps
             smoo=drho*drho*(0.75_real_8-0.25_real_8*drho)
             dsmoo=1.5_real_8*drho*ogceps*(1._real_8-0.5_real_8*drho)
          ENDIF
          grho=SQRT(grad(i,1))
          xx=two13*grho*rho43on
          yy=xx*xx
          sa2b8=SQRT(1._real_8+yy)
          shm1=LOG(xx+sa2b8)
          dd=1._real_8+6._real_8*b1*xx*shm1
          ddon=1.0_real_8/dd
          dd2on=ddon*ddon
          eex=-b1*yy*ddon*two13on
          ee=6._real_8*b1*yy/sa2b8 - 1._real_8
          vexc1=-1.333333333333_real_8*two13on*b1*yy*ee*dd2on
          vexc2=two13*b1*(ee-dd)*dd2on
          v(i)=dsmoo*eex*rho43+smoo*vexc1*x
          vtmp(i)=smoo*vexc2/rho43
          sgcx=sgcx+(smoo*eex*rho43)
          id=id+1
       ENDIF
    ENDDO
    ! ..counting floating point operations
    fadd   = id *  7._real_8
    fmult  = id * 25._real_8
    fdiv   = id *  7._real_8
    fspez  = id *  3._real_8
    flops=fadd+fmult+fdiv+fspez
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gcxonly
  ! ==================================================================
  SUBROUTINE gcxlyp(b1,sgcx,sgcc,rhoe,v,vtmp,grad,flops)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: b1, sgcx, sgcc
    REAL(real_8), DIMENSION(fpar%kr1*fpar%&
      kr2s*fpar%kr3s)                        :: rhoe
    COMPLEX(real_8) :: v(fpar%kr1*fpar%kr2s*fpar%kr3s)
    REAL(real_8), DIMENSION(fpar%kr1*fpar%&
      kr2s*fpar%kr3s)                        :: vtmp, grad
    REAL(real_8)                             :: flops

    REAL(real_8), PARAMETER :: a = 0.04918_real_8, b = 0.132_real_8, &
      c = 0.2533_real_8, d = 0.349_real_8 , eps2 = 1.e-20_real_8, &
      small = 1.e-24_real_8 

    INTEGER                                  :: i, id
    REAL(real_8) :: dd, dom, drho, dsmoo, dxl, ee, eecc, eexc1, fadd, fdiv, &
      ff, fmult, fspez, grho, grhos, o1d, ogceps, om, r1, r2, r4, r5, rho43, &
      rhou, sa2b8, shm1, smoo, two13, v1cc, v2cc, vexc1, vexc2, x, xl, xx, yy

! ==--------------------------------------------------------------==

    two13=2._real_8**0.333333333333333_real_8
    sgcx = 0.0_real_8
    sgcc = 0.0_real_8
    ogceps=1._real_8/cntr%gceps
    id=0
    !$omp parallel do private(I) &
    !$omp  private(RHOU,X, RHO43, SMOO, DSMOO, DRHO, GRHO, GRHOS, XX, YY) &
    !$omp  private(SA2B8,SHM1, DD, EEXC1, EE, VEXC1, VEXC2, R1, R2, R4, R5) &
    !$omp  private(O1D,OM, XL, FF, EECC, DOM, DXL, V1CC, V2CC) &
    !$omp  reduction(+:SGCX,SGCC,ID)
    DO i=1,fpar%kr1*fpar%kr2s*fpar%kr3s
       rhou=MAX(rhoe(i),small)
       IF (rhou.GT.small) THEN
          x=rhou**(0.333333333333333_real_8)
          rho43=x*rhou
          IF (rhou.GT.2.25_real_8*cntr%gceps) THEN
             smoo=1._real_8
             dsmoo=0._real_8
          ELSEIF (rhou.LT.0.25_real_8*cntr%gceps) THEN
             smoo=0._real_8
             dsmoo=0._real_8
          ELSE
             drho=(rhou-0.25_real_8*cntr%gceps)*ogceps
             smoo=drho*drho*(0.75_real_8-0.25_real_8*drho)
             dsmoo=1.5_real_8*drho*ogceps*(1._real_8-0.5_real_8*drho)
          ENDIF
          grho=grad(i)
          grhos=SQRT(grho)
          ! ..exchange
          xx=two13*grhos/rho43
          yy=xx*xx
          sa2b8=SQRT(1._real_8+yy)
          shm1=LOG(xx+sa2b8)
          dd=1._real_8+6._real_8*b1*xx*shm1
          eexc1=-b1*yy/dd/two13
          ee=6._real_8*b1*yy/sa2b8 - 1._real_8
          vexc1=-1.333333333333_real_8/two13*b1*yy*ee/(dd*dd)
          vexc2=two13*b1*(ee-dd)/(dd*dd)
          ! ..correlation lyp
          r1=1._real_8/x
          r2=r1*r1
          r4=r2*r2
          r5=r4*r1
          o1d=1._real_8/(1._real_8+d*r1)
          om=EXP(-c*r1)*o1d
          xl=1._real_8+2.3333333333333_real_8*(c*r1+d*r1*o1d)
          ff=0.04166666666667_real_8*a*b*grho
          eecc=ff*r5*om*xl
          dom=-om*(c+d+c*d*r1)*o1d
          dxl=2.3333333333333_real_8*(c+d+2._real_8*c*d*r1+c*d*d*r2)*o1d*o1d
          v1cc=-0.3333333333333_real_8*ff*r4*(5._real_8*r4*om*xl+&
               r5*dom*xl+r5*om*dxl)
          v2cc=0.08333333333333_real_8*a*b*r5*om*xl
          ! ..sum up results
          id=id+1
          sgcx=sgcx+(smoo*eexc1*rho43)
          sgcc=sgcc+(smoo*eecc)
          v(i)=dsmoo*eexc1*rho43+smoo*vexc1*x
          v(i)=v(i)+dsmoo*eecc+smoo*v1cc
          vtmp(i)=smoo*(vexc2/rho43+v2cc)
       ENDIF
    ENDDO
    ! ..counting floating point operations
    fadd   = id * 23._real_8
    fmult  = id * 70._real_8
    fdiv   = id *  9._real_8
    fspez  = id *  4._real_8
    flops=fadd+fmult+fdiv+fspez
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gcxlyp
  ! ==================================================================
  SUBROUTINE gcxp86(b1,sgcx,sgcc,rhoe,v,vtmp,grad,flops)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: b1, sgcx, sgcc
    REAL(real_8), DIMENSION(fpar%kr1*fpar%&
      kr2s*fpar%kr3s)                        :: rhoe
    COMPLEX(real_8) :: v(fpar%kr1*fpar%kr2s*fpar%kr3s)
    REAL(real_8), DIMENSION(fpar%kr1*fpar%&
      kr2s*fpar%kr3s)                        :: vtmp
    REAL(real_8), DIMENSION(fpar%kr1*fpar%&
      kr2s*fpar%kr3s, *)                     :: grad
    REAL(real_8)                             :: flops

    REAL(real_8), PARAMETER :: eps2 = 1.e-20_real_8, &
      fp86 = -0.206783495048375038_real_8 , p1 = 0.023266_real_8, &
      p2 = 7.389e-6_real_8, p3 = 8.723_real_8, p4 = 0.472_real_8 , &
      pc1 = 0.001667_real_8, pc2 = 0.002568_real_8, pci = pc1+pc2 , &
      rsfac = 0.6203504908994000_real_8, small = 1.e-24_real_8 

    INTEGER                                  :: i, id, kkk
    REAL(real_8) :: cn, cna, cnb, cnbon, dcn, dcna, dcnb, dd, dd2on, ddon, &
      drho, drs, dsmoo, ee, eex, ephi, fadd, fdiv, fmult, fspez, grho, &
      ogceps, orho43, phi, rho43, rhou, rs, sa2b8, sc, shm1, smoo, two13, &
      two13on, v1c, v2c, vexc1, vexc2, x, xx, yy

! ==--------------------------------------------------------------==

    two13=2._real_8**(1.0_real_8/3.0_real_8)
    two13on=1.0_real_8/two13
    sgcx = 0.0_real_8
    sgcc = 0.0_real_8
    ogceps=1._real_8/cntr%gceps
    id=0
    kkk=fpar%kr1*fpar%kr2s*fpar%kr3s
    !$omp parallel do private(I) &
    !$omp  private(RHOU,X,RHO43,ORHO43,SMOO,DSMOO,DRHO,GRHO,XX,YY) &
    !$omp  private(SA2B8,SHM1,DD,DDON,DD2ON,EEX,EE,VEXC1,VEXC2,RS) &
    !$omp  private(CNA,CNB,CNBON,CN,DRS,DCNA,DCNB,DCN,PHI,EPHI,SC) &
    !$omp  private(V1C,V2C) &
    !$omp  reduction(+:SGCX,SGCC,ID)
    DO i=1,kkk
       rhou=MAX(rhoe(i),small)
       IF (rhou.GT.small) THEN
          x=rhou**(0.333333333333333_real_8)
          rho43=x*rhou
          orho43=1._real_8/rho43
          IF (rhou.GT.2.25_real_8*cntr%gceps) THEN
             smoo=1._real_8
             dsmoo=0._real_8
          ELSEIF (rhou.LT.0.25_real_8*cntr%gceps) THEN
             smoo=0._real_8
             dsmoo=0._real_8
          ELSE
             drho=(rhou-0.25_real_8*cntr%gceps)*ogceps
             smoo=drho*drho*(0.75_real_8-0.25_real_8*drho)
             dsmoo=1.5_real_8*drho*ogceps*(1._real_8-0.5_real_8*drho)
          ENDIF
          grho=SQRT(grad(i,1))
          ! ..exchange
          xx=two13*grho*orho43
          yy=xx*xx
          sa2b8=SQRT(1._real_8+yy)
          shm1=LOG(xx+sa2b8)
          dd=1._real_8+6._real_8*b1*xx*shm1
          ddon=1.0_real_8/dd
          dd2on=ddon*ddon
          eex=-b1*yy*ddon*two13on
          ee=6._real_8*b1*yy/sa2b8 - 1._real_8
          vexc1=-1.333333333333_real_8*two13on*b1*yy*ee*dd2on
          vexc2=two13*b1*(ee-dd)*dd2on
          ! .. perdew 86 correlation
          xx   = grho*orho43
          rs   = rsfac/x
          cna  = pc2+rs*(p1+p2*rs)
          cnb  = 1._real_8+rs*(p3+rs*(p4+1.e4_real_8*p2*rs))
          cnbon=1.0_real_8/cnb
          cn   = pc1 + cna*cnbon
          drs  = fp86*orho43
          dcna = (p1+2._real_8*p2*rs)*drs
          dcnb = (p3+rs*(2._real_8*p4+3.e4_real_8*p2*rs))*drs
          dcn  = dcna*cnbon-cna*cnbon*cnbon*dcnb
          phi  = 0.192_real_8*pci/cn*grho*orho43*SQRT(x)
          ephi = EXP(-phi)
          sc   = xx*grho*cn*ephi
          v1c  = sc*((1._real_8+phi)*dcn/cn -&
               (1.33333333333333_real_8-1.166666666666667_real_8*phi)/rhou)
          v2c  = cn*ephi*orho43*(2._real_8-phi)
          v(i) = dsmoo*(sc+eex*rho43)+smoo*(vexc1*x+v1c)
          vtmp(i)=smoo*(vexc2*orho43+v2c)
          id=id+1
          sgcx=sgcx+(smoo*eex*rho43)
          sgcc=sgcc+(smoo*sc)
       ENDIF
    ENDDO
    ! ..counting floating point operations
    fadd   = id * 24._real_8
    fmult  = id * 56._real_8
    fdiv   = id * 13._real_8
    fspez  = id *  5._real_8
    flops=fadd+fmult+fdiv+fspez
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gcxp86
  ! ==================================================================
  SUBROUTINE gcgga(sgcx,sgcc,rhoe,v,vtmp,grad,flops)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: sgcx, sgcc
    REAL(real_8), DIMENSION(fpar%kr1*fpar%&
      kr2s*fpar%kr3s)                        :: rhoe
    COMPLEX(real_8) :: v(fpar%kr1*fpar%kr2s*fpar%kr3s)
    REAL(real_8), DIMENSION(fpar%kr1*fpar%&
      kr2s*fpar%kr3s)                        :: vtmp
    REAL(real_8), DIMENSION(fpar%kr1*fpar%&
      kr2s*fpar%kr3s, *)                     :: grad
    REAL(real_8)                             :: flops

    REAL(real_8), PARAMETER :: a0u = 0.4581652932831429_real_8, &
      a1u = 2.217058676663745_real_8, a2u = 0.7405551735357053_real_8, &
      a3u = 0.01968227878617998_real_8 , al = 0.09_real_8, b1u = 1.0_real_8, &
      b2u = 4.504130959426697_real_8, b3u = 1.110667363742916_real_8, &
      b4u = 0.02359291751427506_real_8 , c1u = 4._real_8*a0u*b1u/3.0_real_8, &
      c2u = 5.0_real_8*a0u*b2u/3.0_real_8+a1u*b1u, c3u = 2.0_real_8*a0u*&
      b3u+4.0_real_8*a1u*b2u/3.0_real_8+2.0_real_8*a2u*b1u/3.0_real_8, c4u = &
      7.0_real_8*a0u*b4u/3.0_real_8+5.0_real_8*a1u*b3u/3.0_real_8+a2u*b2u+a3u*&
      b1u/3.0_real_8
    REAL(real_8), PARAMETER :: c5u = 2.0_real_8*a1u*b4u+4.0_real_8*a2u*b3u/&
      3.0_real_8+2.0_real_8*a3u*b2u/3.0_real_8, &
      c6u = 5.0_real_8*a2u*b4u/3.0_real_8+a3u*b3u, &
      c7u = 4.0_real_8*a3u*b4u/3.0_real_8 , cx = -0.001667_real_8, &
      cxc0 = 0.002568_real_8, cc0 = -cx+cxc0 , eps2 = 1.e-20_real_8, &
      f1 = 0.19645_real_8, f1x = -1.10783814957303361_real_8 , &
      f2 = 7.7956_real_8, f3 = 0.2743_real_8, f4 = 0.1508_real_8, &
      f5 = 0.004_real_8 , pa = 0.023266_real_8, pb = 7.389e-6_real_8, &
      pc = 8.723_real_8, pd = 0.472_real_8 , &
      rsfac = 0.6203504908994000_real_8 , small = 1.e-24_real_8 

    INTEGER                                  :: i, id, kkk
    REAL(real_8) :: af, as, be, beon, bf, botu, botuon, bs, bson, cn, cna, &
      cnb, cnbon, das, dbs, dcn, dcna, dcnb, ddh0, ddh1, dh0, dh1, dls, drho, &
      dsmoo, dtopu, ee, epsxcu, expe, exps, fa1, fa2, fadd, fdiv, fmult, fp1, &
      fp2, fspez, grho, grhos, h0, h1, ogceps, qy, rho43, rhou, rr, rs2, rs3, &
      rsu, s, s1, s2, s3, s4, sa2b8, sc, shm1, smoo, son, sx, t, topu, v1c, &
      v1x, v2c, v2x, vxc, x, xkf, xkff, xks, xnu, xy, y

! ==--------------------------------------------------------------==

    fa1=f1x*0.66666666666666_real_8
    fa2=fa1*1.33333333333333_real_8
    fp1   = -3._real_8/(16._real_8*pi)*(3._real_8*pi*pi)**(-0.33333333333333_real_8)
    fp2   = 0.5_real_8*(3._real_8*pi*pi)**(-0.33333333333333_real_8)
    xnu   = 16._real_8/pi*(3._real_8*pi*pi)**0.333333333333333_real_8
    be    = xnu*cc0
    beon=1.0_real_8/be
    xkff  = (2.25_real_8*pi)**0.333333333333333_real_8
    sgcx = 0.0_real_8
    sgcc = 0.0_real_8
    ogceps=1._real_8/cntr%gceps
    id=0
    kkk=fpar%kr1*fpar%kr2s*fpar%kr3s
    !$omp parallel do private(I) &
    !$omp  private(RHOU,X,RHO43,SMOO,DSMOO,DRHO,GRHO,GRHOS,RR,S1,S2) &
    !$omp  private(S3,S4,EXPS,AS,SA2B8,SHM1,BS,BSON,DAS,DBS,DLS,SX) &
    !$omp  private(V1X,V2X,RSU,TOPU,DTOPU,BOTU,BOTUON,EPSXCU,VXC,RS2) &
    !$omp  private(RS3,XKF,XKS,T,EXPE,AF,BF,Y,XY,QY,S,SON,H0,DH0) &
    !$omp  private(DDH0,EE,CNA,DCNA,CNB,CNBON,DCNB,CN,DCN,H1,DH1) &
    !$omp  private(DDH1,SC,V1C,V2C) &
    !$omp  reduction(+:SGCX,SGCC,ID)
    DO i=1,kkk
       rhou=MAX(rhoe(i),small)
       IF (rhou.GT.small) THEN
          x=rhou**(0.33333333333333_real_8)
          rho43=x*rhou
          IF (rhou.GT.2.25_real_8*cntr%gceps) THEN
             smoo=1._real_8
             dsmoo=0._real_8
          ELSEIF (rhou.LT.0.25_real_8*cntr%gceps) THEN
             smoo=0._real_8
             dsmoo=0._real_8
          ELSE
             drho=(rhou-0.25_real_8*cntr%gceps)*ogceps
             smoo=drho*drho*(0.75_real_8-0.25_real_8*drho)
             dsmoo=1.5_real_8*drho*ogceps*(1._real_8-0.5_real_8*drho)
          ENDIF
          grho=grad(i,1)
          grhos=SQRT(grho)
          ! ..exchange part
          rr=1._real_8/rho43
          s1=fp2*grhos*rr
          s2=s1*s1
          s3=s1*s2
          s4=s2*s2
          exps=f4*EXP(-100._real_8*s2)
          as=f3-exps-f5*s2
          sa2b8=SQRT(1._real_8+f2*f2*s2)
          shm1=LOG(f2*s1+sa2b8)
          bs=1._real_8+f1*s1*shm1+f5*s4
          bson=1.0_real_8/bs
          das=(200._real_8*exps-2._real_8*f5)*s1
          dbs=f1*(shm1+f2*s1/sa2b8)+4._real_8*f5*s3
          dls=(das/as-dbs*bson)
          sx=fp1*grho*rr*as*bson
          v1x=-1.333333333333333_real_8*sx/rhou*(1._real_8+s1*dls)
          v2x=fp1*rr*as*bson*(2._real_8+s1*dls)
          ! ..lda functional
          rsu=rsfac/x
          topu=a0u+rsu*(a1u+rsu*(a2u+rsu*a3u))
          dtopu=-rsu*(c1u+rsu*(c2u+rsu*(c3u+rsu*(c4u+rsu*&
               (c5u+rsu*(c6u+rsu*c7u))))))
          botu=rsu*(b1u+rsu*(b2u+rsu*(b3u+rsu*b4u)))
          botuon=1.0_real_8/botu
          epsxcu=-topu*botuon - fa1*x
          vxc=dtopu*botuon*botuon - fa2*x
          ! ..correlation part
          rs2=rsu*rsu
          rs3=rs2*rsu
          xkf=xkff/rsu
          xks=SQRT(4._real_8*xkf/pi)
          t  =grhos/(2._real_8*xks*rhou)
          expe=EXP(-2._real_8*al*epsxcu*beon*beon)
          af=2._real_8*al*beon/(expe-1._real_8)
          bf=expe*(vxc-epsxcu)
          y = af*t*t
          xy = (1._real_8+y)/(1._real_8+y+y*y)
          qy = y*y*(2._real_8+y)/(1._real_8+y+y*y)**2
          s  = 1._real_8+2._real_8*al*beon*t*t*xy
          son=1.0_real_8/s
          h0 = be*be/(2._real_8*al)*LOG(s)
          dh0= be*t*t*son*(-2.33333333333333_real_8*xy-qy*&
               (af*bf*beon-2.3333333333333333_real_8))
          ddh0=be/(2._real_8*xks*xks*rhou)*(xy-qy)*son
          ee = -100._real_8*(xks/xkf*t)**2
          cna= cxc0+pa*rsu+pb*rs2
          dcna= -(pa*rsu+2._real_8*pb*rs2)*0.33333333333333_real_8
          cnb  = 1._real_8+pc*rsu+pd*rs2+1.e4_real_8*pb*rs3
          cnbon=1.0_real_8/cnb
          dcnb = -(pc*rsu+2._real_8*pd*rs2+3.e4_real_8*pb*rs3)*&
               0.3333333333333333_real_8
          cn   = cna*cnbon - cx
          dcn  = dcna*cnbon - cna*dcnb*cnbon*cnbon
          h1   = xnu*(cn-cc0-3._real_8/7._real_8*cx)*t*t*EXP(ee)
          dh1  = -0.3333333333333_real_8*(h1*(7._real_8+8._real_8*ee)+&
               xnu*t*t*EXP(ee)*dcn)
          ddh1 = 2._real_8*h1*(1._real_8+ee)*rhou/grho
          sc   = rhou*(h0+h1)
          v1c  = h0+h1+dh0+dh1
          v2c  = ddh0+ddh1
          ! ..sum up terms
          v(i)=dsmoo*(sx+sc)+smoo*(v1x+v1c)
          vtmp(i)=smoo*(v2x+v2c)
          sgcx=sgcx+(smoo*sx)
          sgcc=sgcc+(smoo*sc)
          id=id+1
       ENDIF
    ENDDO
    ! ..counting floating point operations
    fadd   = id *  48._real_8
    fmult  = id * 108._real_8
    fdiv   = id *  27._real_8
    fspez  = id *   9._real_8
    flops=fadd+fmult+fdiv+fspez
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gcgga
  ! ==================================================================
  SUBROUTINE gcpbe(sgcx,sgcc,rhoe,v,vtmp,grad,flops)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: sgcx, sgcc
    REAL(real_8), DIMENSION(fpar%kr1*fpar%&
      kr2s*fpar%kr3s)                        :: rhoe
    COMPLEX(real_8) :: v(fpar%kr1*fpar%kr2s*fpar%kr3s)
    REAL(real_8), DIMENSION(fpar%kr1*fpar%&
      kr2s*fpar%kr3s)                        :: vtmp
    REAL(real_8), DIMENSION(fpar%kr1*fpar%&
      kr2s*fpar%kr3s, *)                     :: grad
    REAL(real_8)                             :: flops

    REAL(real_8), PARAMETER :: a0u = 0.4581652932831429_real_8, &
      a1u = 2.217058676663745_real_8, a2u = 0.7405551735357053_real_8, &
      a3u = 0.01968227878617998_real_8 , &
      al = 0.715996577859519256e-01_real_8, &
      ax = -0.738558766382022406_real_8, b1u = 1.0_real_8, &
      b2u = 4.504130959426697_real_8, b3u = 1.110667363742916_real_8, &
      b4u = 0.02359291751427506_real_8 , be = 0.06672455060314922_real_8 , &
      c1u = 4._real_8*a0u*b1u/3.0_real_8, &
      c2u = 5.0_real_8*a0u*b2u/3.0_real_8+a1u*b1u, c3u = 2.0_real_8*a0u*&
      b3u+4.0_real_8*a1u*b2u/3.0_real_8+2.0_real_8*a2u*b1u/3.0_real_8
    REAL(real_8), PARAMETER :: c4u = 7.0_real_8*a0u*b4u/3.0_real_8+5.0_real_8*&
      a1u*b3u/3.0_real_8+a2u*b2u+a3u*b1u/3.0_real_8, c5u = 2.0_real_8*a1u*&
      b4u+4.0_real_8*a2u*b3u/3.0_real_8+2.0_real_8*a3u*b2u/3.0_real_8, &
      c6u = 5.0_real_8*a2u*b4u/3.0_real_8+a3u*b3u, &
      c7u = 4.0_real_8*a3u*b4u/3.0_real_8 , eps2 = 1.e-20_real_8, &
      f1 = 0.19645_real_8, f1x = -1.10783814957303361_real_8 , &
      f2 = 7.7956_real_8, f3 = 0.2743_real_8, f4 = 0.1508_real_8, &
      f5 = 0.004_real_8 , rsfac = 0.6203504908994000_real_8 , &
      small = 1.e-24_real_8 , uk = 0.8040_real_8, &
      um = 0.2195149727645171_real_8, ul = um/uk 
    REAL(real_8), PARAMETER :: us = 0.161620459673995492_real_8

    INTEGER                                  :: i, id, kkk
    REAL(real_8) :: af, beon, bf, botu, botuon, ddh0, dfx, dh0, drho, dsmoo, &
      dtopu, epsxcu, ex, expe, fa1, fa2, fadd, fdiv, fmult, fp1, fp2, fspez, &
      fx, grho, grhos, h0, ogceps, po, qy, rho43, rhou, rr, rsu, s, s2, sc, &
      smoo, son, sx, t, topu, v1c, v1x, v2c, v2x, vxc, x, xkf, xkff, xks, xy, &
      y, yyy1

! ==--------------------------------------------------------------==

    fa1=f1x*0.66666666666666_real_8
    fa2=fa1*1.33333333333333_real_8
    fp1   = -3._real_8/(16._real_8*pi)*(3._real_8*pi*pi)**(-0.33333333333333_real_8)
    fp2   = 0.5_real_8*(3._real_8*pi*pi)**(-0.333333333333333_real_8)
    xkff  = (2.25_real_8*pi)**0.333333333333333_real_8
    beon=1.0_real_8/be
    sgcx = 0.0_real_8
    sgcc = 0.0_real_8
    ogceps=1._real_8/cntr%gceps
    id=0
    kkk=fpar%kr1*fpar%kr2s*fpar%kr3s
    !$omp parallel do private(I) &
    !$omp  private(RHOU,X,RHO43,SMOO,DSMOO,DRHO,GRHO,GRHOS,RR,EX,S2) &
    !$omp  private(PO,FX,SX,DFX,V1X,V2X,RSU,TOPU,DTOPU,BOTU,BOTUON) &
    !$omp  private(EPSXCU,VXC,XKF,XKS,T,EXPE,AF,BF,Y,YYY1,XY,QY) &
    !$omp  private(S,SON,H0,DH0,DDH0,SC,V1C) &
    !$omp  private(V2C) &
    !$omp  reduction(+:SGCX,SGCC,ID)
    DO i=1,kkk
       rhou=MAX(rhoe(i),small)
       IF (rhou.GT.small) THEN
          x=rhou**(0.33333333333333_real_8)
          rho43=x*rhou
          IF (rhou.GT.2.25_real_8*cntr%gceps) THEN
             smoo=1._real_8
             dsmoo=0._real_8
          ELSEIF (rhou.LT.0.25_real_8*cntr%gceps) THEN
             smoo=0._real_8
             dsmoo=0._real_8
          ELSE
             drho=(rhou-0.25_real_8*cntr%gceps)*ogceps
             smoo=drho*drho*(0.75_real_8-0.25_real_8*drho)
             dsmoo=1.5_real_8*drho*ogceps*(1._real_8-0.5_real_8*drho)
          ENDIF
          grho=grad(i,1)
          grhos=SQRT(grho)
          ! ..exchange part
          rr    = 1._real_8/rho43
          ex    = ax*rho43
          s2    = grho*rr*rr*us*us
          po    = 1._real_8/(1._real_8 + ul*s2)
          fx    = uk-uk*po
          sx    = ex*fx
          dfx   = 2._real_8*uk*ul*po*po
          v1x   = 1.33333333333333_real_8*ax*x*(fx-s2*dfx)
          v2x   = ex*dfx*(us*rr)**2
          ! ..lda functional
          rsu=rsfac/x
          topu=a0u+rsu*(a1u+rsu*(a2u+rsu*a3u))
          dtopu=-rsu*(c1u+rsu*(c2u+rsu*(c3u+rsu*(c4u+rsu*&
               (c5u+rsu*(c6u+rsu*c7u))))))
          botu=rsu*(b1u+rsu*(b2u+rsu*(b3u+rsu*b4u)))
          botuon=1.0_real_8/botu
          epsxcu=topu*(-botuon) - fa1*x
          vxc=dtopu*botuon*botuon - fa2*x
          ! ..correlation part
          xkf=xkff/rsu
          xks=SQRT(4._real_8*xkf/pi)
          t  =grhos/(2._real_8*xks*rhou)
          expe=EXP(-2._real_8*al*epsxcu*beon*beon)
          af=2._real_8*al*beon/(expe-1._real_8)
          bf=expe*(vxc-epsxcu)
          y = af*t*t
          yyy1=1.0_real_8/(1.0_real_8+y+y*y)
          xy = (1.0_real_8+y)*yyy1
          qy = y*y*(2.0_real_8+y)*yyy1**2
          s  = 1._real_8+2._real_8*al*beon*t*t*xy
          son=1.0_real_8/s
          h0 = be*be/(2._real_8*al)*LOG(s)
          dh0= be*t*t*son*(-2.33333333333333_real_8*xy-qy*&
               (af*bf*beon-2.3333333333333333_real_8))
          ddh0=be/(2._real_8*xks*xks*rhou)*(xy-qy)*son
          sc   = rhou*h0
          v1c  = h0+dh0
          v2c  = ddh0
          ! ..sum up terms
          v(i)=dsmoo*(sx+sc)+smoo*(v1x+v1c)
          vtmp(i)=smoo*(v2x+v2c)
          sgcx=sgcx+(smoo*sx)
          sgcc=sgcc+(smoo*sc)
          id=id+1
       ENDIF
    ENDDO
    ! ..counting floating point operations
    fadd   = id *  33._real_8
    fmult  = id *  77._real_8
    fdiv   = id *  18._real_8
    fspez  = id *   4._real_8
    flops=fadd+fmult+fdiv+fspez
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gcpbe
  ! ==================================================================
  SUBROUTINE gcrevpbe(sgcx,sgcc,rhoe,v,vtmp,grad,flops)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: sgcx, sgcc
    REAL(real_8), DIMENSION(fpar%kr1*fpar%&
      kr2s*fpar%kr3s)                        :: rhoe
    COMPLEX(real_8) :: v(fpar%kr1*fpar%kr2s*fpar%kr3s)
    REAL(real_8), DIMENSION(fpar%kr1*fpar%&
      kr2s*fpar%kr3s)                        :: vtmp
    REAL(real_8), DIMENSION(fpar%kr1*fpar%&
      kr2s*fpar%kr3s, *)                     :: grad
    REAL(real_8)                             :: flops

    REAL(real_8), PARAMETER :: a0u = 0.4581652932831429_real_8, &
      a1u = 2.217058676663745_real_8, a2u = 0.7405551735357053_real_8, &
      a3u = 0.01968227878617998_real_8 , &
      al = 0.715996577859519256e-01_real_8, &
      ax = -0.738558766382022406_real_8, b1u = 1.0_real_8, &
      b2u = 4.504130959426697_real_8, b3u = 1.110667363742916_real_8, &
      b4u = 0.02359291751427506_real_8 , be = 0.06672455060314922_real_8 , &
      c1u = 4._real_8*a0u*b1u/3.0_real_8, &
      c2u = 5.0_real_8*a0u*b2u/3.0_real_8+a1u*b1u, c3u = 2.0_real_8*a0u*&
      b3u+4.0_real_8*a1u*b2u/3.0_real_8+2.0_real_8*a2u*b1u/3.0_real_8
    REAL(real_8), PARAMETER :: c4u = 7.0_real_8*a0u*b4u/3.0_real_8+5.0_real_8*&
      a1u*b3u/3.0_real_8+a2u*b2u+a3u*b1u/3.0_real_8, c5u = 2.0_real_8*a1u*&
      b4u+4.0_real_8*a2u*b3u/3.0_real_8+2.0_real_8*a3u*b2u/3.0_real_8, &
      c6u = 5.0_real_8*a2u*b4u/3.0_real_8+a3u*b3u, &
      c7u = 4.0_real_8*a3u*b4u/3.0_real_8 , eps2 = 1.e-20_real_8, &
      f1 = 0.19645_real_8, f1x = -1.10783814957303361_real_8 , &
      f2 = 7.7956_real_8, f3 = 0.2743_real_8, f4 = 0.1508_real_8, &
      f5 = 0.004_real_8 , rsfac = 0.6203504908994000_real_8 , &
      small = 1.e-24_real_8 , uk = 1.2450_real_8, &
      um = 0.2195149727645171_real_8, ul = um/uk 
    REAL(real_8), PARAMETER :: us = 0.161620459673995492_real_8

    INTEGER                                  :: i, id, kkk
    REAL(real_8) :: af, beon, bf, botu, botuon, ddh0, dfx, dh0, drho, dsmoo, &
      dtopu, epsxcu, ex, expe, fa1, fa2, fadd, fdiv, fmult, fp1, fp2, fspez, &
      fx, grho, grhos, h0, ogceps, po, qy, rho43, rhou, rr, rsu, s, s2, sc, &
      smoo, son, sx, t, topu, v1c, v1x, v2c, v2x, vxc, x, xkf, xkff, xks, xy, &
      y, yyy1

! ==--------------------------------------------------------------==

    fa1=f1x*0.66666666666666_real_8
    fa2=fa1*1.33333333333333_real_8
    fp1   = -3._real_8/(16._real_8*pi)*(3._real_8*pi*pi)**(-0.33333333333333_real_8)
    fp2   = 0.5_real_8*(3._real_8*pi*pi)**(-0.333333333333333_real_8)
    xkff  = (2.25_real_8*pi)**0.333333333333333_real_8
    beon=1.0_real_8/be
    sgcx = 0.0_real_8
    sgcc = 0.0_real_8
    ogceps=1._real_8/cntr%gceps
    id=0
    kkk=fpar%kr1*fpar%kr2s*fpar%kr3s
    !$omp parallel do private(I) &
    !$omp  private(RHOU,X,RHO43,SMOO,DSMOO,DRHO,GRHO,GRHOS,RR,EX,S2) &
    !$omp  private(PO,FX,SX,DFX,V1X,V2X,RSU,TOPU,DTOPU,BOTU,BOTUON) &
    !$omp  private(EPSXCU,VXC,XKF,XKS,T,EXPE,AF,BF,Y,YYY1,XY,QY) &
    !$omp  private(S,SON,H0,DH0,DDH0,SC,V1C) &
    !$omp  private(V2C) &
    !$omp  reduction(+:SGCX,SGCC,ID)
    DO i=1,kkk
       rhou=MAX(rhoe(i),small)
       IF (rhou.GT.small) THEN
          x=rhou**(0.33333333333333_real_8)
          rho43=x*rhou
          IF (rhou.GT.2.25_real_8*cntr%gceps) THEN
             smoo=1._real_8
             dsmoo=0._real_8
          ELSEIF (rhou.LT.0.25_real_8*cntr%gceps) THEN
             smoo=0._real_8
             dsmoo=0._real_8
          ELSE
             drho=(rhou-0.25_real_8*cntr%gceps)*ogceps
             smoo=drho*drho*(0.75_real_8-0.25_real_8*drho)
             dsmoo=1.5_real_8*drho*ogceps*(1._real_8-0.5_real_8*drho)
          ENDIF
          grho=grad(i,1)
          grhos=SQRT(grho)
          ! ..exchange part
          rr    = 1._real_8/rho43
          ex    = ax*rho43
          s2    = grho*rr*rr*us*us
          po    = 1._real_8/(1._real_8 + ul*s2)
          fx    = uk-uk*po
          sx    = ex*fx
          dfx   = 2._real_8*uk*ul*po*po
          v1x   = 1.33333333333333_real_8*ax*x*(fx-s2*dfx)
          v2x   = ex*dfx*(us*rr)**2
          ! ..lda functional
          rsu=rsfac/x
          topu=a0u+rsu*(a1u+rsu*(a2u+rsu*a3u))
          dtopu=-rsu*(c1u+rsu*(c2u+rsu*(c3u+rsu*(c4u+rsu*&
               (c5u+rsu*(c6u+rsu*c7u))))))
          botu=rsu*(b1u+rsu*(b2u+rsu*(b3u+rsu*b4u)))
          botuon=1.0_real_8/botu
          epsxcu=topu*(-botuon) - fa1*x
          vxc=dtopu*botuon*botuon - fa2*x
          ! ..correlation part
          xkf=xkff/rsu
          xks=SQRT(4._real_8*xkf/pi)
          t  =grhos/(2._real_8*xks*rhou)
          expe=EXP(-2._real_8*al*epsxcu*beon*beon)
          af=2._real_8*al*beon/(expe-1._real_8)
          bf=expe*(vxc-epsxcu)
          y = af*t*t
          yyy1=1.0_real_8/(1.0_real_8+y+y*y)
          xy = (1.0_real_8+y)*yyy1
          qy = y*y*(2.0_real_8+y)*yyy1**2
          s  = 1._real_8+2._real_8*al*beon*t*t*xy
          son=1.0_real_8/s
          h0 = be*be/(2._real_8*al)*LOG(s)
          dh0= be*t*t*son*(-2.33333333333333_real_8*xy-qy*&
               (af*bf*beon-2.3333333333333333_real_8))
          ddh0=be/(2._real_8*xks*xks*rhou)*(xy-qy)*son
          sc   = rhou*h0
          v1c  = h0+dh0
          v2c  = ddh0
          ! ..sum up terms
          v(i)=dsmoo*(sx+sc)+smoo*(v1x+v1c)
          vtmp(i)=smoo*(v2x+v2c)
          sgcx=sgcx+(smoo*sx)
          sgcc=sgcc+(smoo*sc)
          id=id+1
       ENDIF
    ENDDO
    ! ..counting floating point operations
    fadd   = id *  33._real_8
    fmult  = id *  77._real_8
    fdiv   = id *  18._real_8
    fspez  = id *   4._real_8
    flops=fadd+fmult+fdiv+fspez
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gcrevpbe
  ! ==================================================================
  SUBROUTINE gcsxonly(b1,sgcx,sgcc,rhoe,v1,v2,vtmp,grad,flops)
    REAL(real_8)                             :: b1, sgcx, sgcc
    REAL(real_8), DIMENSION(fpar%kr1*fpar%&
      kr2s*fpar%kr3s, 2)                     :: rhoe
    COMPLEX(real_8) :: v1(fpar%kr1*fpar%kr2s*fpar%kr3s), &
      v2(fpar%kr1*fpar%kr2s*fpar%kr3s)
    REAL(real_8), DIMENSION(fpar%kr1*fpar%&
      kr2s*fpar%kr3s, 3)                     :: vtmp
    REAL(real_8), DIMENSION(fpar%kr1*fpar%&
      kr2s*fpar%kr3s, 8)                     :: grad
    REAL(real_8)                             :: flops

    REAL(real_8), PARAMETER                  :: eps2 = 1.e-20_real_8, &
                                                small = 1.e-24_real_8 

    INTEGER                                  :: i, id, kkk
    REAL(real_8) :: cn, cna, cnb, cnbon, d, dcn, dcna, dcnb, dda, ddaon, ddb, &
      ddbon, ddda, dddb, dgfa, dgfb, don, drho, drs, dsmoo, ephi, fadd, fdiv, &
      fmult, fspez, gfa, gfb, grho, grhoa, grhob, ogceps, phi, rho, rhoa, &
      rhob, rhoon, rs, rs2, rs3, rsa, rsb, s1, s1on, s2, sa2a, sa2b, sc, &
      shma, shmb, smoo, sxa, sxb, v1c, v1ca, v1cb, v1xa, v1xb, v2c, v2xa, &
      v2xb, x, xsa, xsaa, xsb, xsbb, y, z, zon, zz

! ==--------------------------------------------------------------==

    ogceps=1._real_8/cntr%gceps
    sgcx=0._real_8
    sgcc=0._real_8
    id=0
    kkk=fpar%kr1*fpar%kr2s*fpar%kr3s
    !$omp parallel do private(I) &
    !$omp  private(RHOA,RHOB,RHO,X,RSA,Y,RSB,Z,ZON,RHOON,ZZ,RS) &
    !$omp  private(SMOO,DSMOO,DRHO,GRHOA,GRHOB,XSA,XSB,XSAA,XSBB) &
    !$omp  private(SA2A,SA2B,SHMA,SHMB,DDA,DDB,DDAON,DDBON,DDDA,DDDB) &
    !$omp  private(GFA,GFB,DGFA,DGFB,SXA,SXB,V1XA,V1XB,V2XA,V2XB,GRHO) &
    !$omp  private(RS2,RS3,CNA,CNB,CNBON,CN,DRS,DCNA,DCNB,DCN,S1,S1ON) &
    !$omp  private(S2,D,DON,PHI,EPHI,SC,V1C,V1CA,V1CB,V2C) &
    !$omp  reduction(+:SGCX,ID)
    DO i=1,kkk
       rhoa=MAX(rhoe(i,1),small)
       rhob=MAX(rhoe(i,2),small)
       rho=rhoa+rhob
       IF (rho.GT.3._real_8*small) THEN
          x=rhoa**(0.33333333333333_real_8)
          rsa=1._real_8/x
          y=rhob**(0.33333333333333_real_8)
          rsb=1._real_8/y
          IF (rho.GT.2.25_real_8*cntr%gceps) THEN
             smoo=1._real_8
             dsmoo=0._real_8
          ELSEIF (rho.LT.0.25_real_8*cntr%gceps) THEN
             smoo=0._real_8
             dsmoo=0._real_8
          ELSE
             drho=(rho-0.25_real_8*cntr%gceps)*ogceps
             smoo=drho*drho*(0.75_real_8-0.25_real_8*drho)
             dsmoo=1.5_real_8*drho*ogceps*(1._real_8-0.5_real_8*drho)
          ENDIF
          grhoa=SQRT(grad(i,1))
          grhob=SQRT(grad(i,5))
          grhoa=MAX(grhoa,small)
          grhob=MAX(grhob,small)
          xsa  = grhoa*(rsa*rsa*rsa*rsa)
          xsb  = grhob*(rsb*rsb*rsb*rsb)
          IF (rhoa.LT.2._real_8*small) xsa=0._real_8
          IF (rhob.LT.2._real_8*small) xsb=0._real_8
          xsaa = xsa*xsa
          xsbb = xsb*xsb
          sa2a = SQRT(1._real_8+xsaa)
          sa2b = SQRT(1._real_8+xsbb)
          shma = LOG(xsa+sa2a)
          shmb = LOG(xsb+sa2b)
          dda  = 1._real_8+6._real_8*b1*xsa*shma
          ddb  = 1._real_8+6._real_8*b1*xsb*shmb
          ddaon=1.0_real_8/dda
          ddbon=1.0_real_8/ddb
          ddda = 6._real_8*b1*(shma+xsa/sa2a)
          dddb = 6._real_8*b1*(shmb+xsb/sa2b)
          gfa  = -b1*xsaa*ddaon
          gfb  = -b1*xsbb*ddbon
          dgfa = (-2._real_8*b1*xsa*dda+b1*xsaa*ddda)*(ddaon*ddaon)
          dgfb = (-2._real_8*b1*xsb*ddb+b1*xsbb*dddb)*(ddbon*ddbon)
          sxa  = gfa*x*rhoa
          sxb  = gfb*y*rhob
          v1xa = 1.33333333333333_real_8*x*(gfa-xsa*dgfa)
          v1xb = 1.33333333333333_real_8*y*(gfb-xsb*dgfb)
          v2xa = dgfa/grhoa
          v2xb = dgfb/grhob
          ! ..sum up terms
          v1(i)=dsmoo*sxa+smoo*v1xa
          v2(i)=dsmoo*sxb+smoo*v1xb
          vtmp(i,1)=smoo*v2xa
          vtmp(i,2)=smoo*v2xb
          vtmp(i,3)=0._real_8
          sgcx=sgcx + (smoo*(sxa+sxb))
          id=id+1
       ENDIF
    ENDDO
    ! ..counting floating point operations
    fadd   = id *  17._real_8
    fmult  = id *  52._real_8
    fdiv   = id *   7._real_8
    fspez  = id *   7._real_8
    flops=fadd+fmult+fdiv+fspez
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gcsxonly
  ! ==================================================================
  SUBROUTINE gcsxp86(b1,sgcx,sgcc,rhoe,v1,v2,vtmp,grad,flops)
    REAL(real_8)                             :: b1, sgcx, sgcc
    REAL(real_8), DIMENSION(fpar%kr1*fpar%&
      kr2s*fpar%kr3s, 2)                     :: rhoe
    COMPLEX(real_8) :: v1(fpar%kr1*fpar%kr2s*fpar%kr3s), &
      v2(fpar%kr1*fpar%kr2s*fpar%kr3s)
    REAL(real_8), DIMENSION(fpar%kr1*fpar%&
      kr2s*fpar%kr3s, 3)                     :: vtmp
    REAL(real_8), DIMENSION(fpar%kr1*fpar%&
      kr2s*fpar%kr3s, 8)                     :: grad
    REAL(real_8)                             :: flops

    REAL(real_8), PARAMETER :: eps2 = 1.e-20_real_8, p1 = 0.023266_real_8, &
      p2 = 7.389e-6_real_8, p3 = 8.723_real_8, p4 = 0.472_real_8 , &
      pc1 = 0.001667_real_8, pc2 = 0.002568_real_8, pci = pc1+pc2 , &
      rsfac = 0.6203504908994000_real_8 , small = 1.e-24_real_8 

    INTEGER                                  :: i, id, kkk
    REAL(real_8) :: cn, cna, cnb, cnbon, d, dcn, dcna, dcnb, dda, ddaon, ddb, &
      ddbon, ddda, dddb, dfx, dgfa, dgfb, don, drho, drs, dsmoo, ephi, fadd, &
      fdiv, fmult, fspez, gfa, gfb, grho, grhoa, grhob, ogceps, phi, rho, &
      rhoa, rhob, rhoon, rs, rs2, rs3, rsa, rsb, s1, s1on, s2, sa2a, sa2b, &
      sc, shma, shmb, smoo, sxa, sxb, two13, v1c, v1ca, v1cb, v1xa, v1xb, &
      v2c, v2ca, v2cab, v2cb, v2xa, v2xb, x, xsa, xsaa, xsb, xsbb, y, z, zon, &
      zz

! ==--------------------------------------------------------------==

    two13=2._real_8**0.333333333333333_real_8
    ogceps=1._real_8/cntr%gceps
    dfx = -0.3333333333333_real_8*(3._real_8/fpi)**0.3333333333333_real_8
    sgcx=0._real_8
    sgcc=0._real_8
    id=0
    kkk=fpar%kr1*fpar%kr2s*fpar%kr3s
    !$omp parallel do private(I) &
    !$omp  private(RHOA,RHOB,RHO,X,RSA,Y,RSB,Z,ZON,RHOON,ZZ,RS) &
    !$omp  private(SMOO,DSMOO,DRHO,GRHOA,GRHOB,XSA,XSB,XSAA,XSBB) &
    !$omp  private(SA2A,SA2B,SHMA,SHMB,DDA,DDB,DDAON,DDBON,DDDA) &
    !$omp  private(DDDB,GFA,GFB,DGFA,DGFB,SXA,SXB,V1XA,V1XB,V2XA) &
    !$omp  private(V2XB,GRHO,RS2,RS3,CNA,CNB,CNBON,CN,DRS,DCNA,DCNB) &
    !$omp  private(DCN,S1,S1ON,S2,D,DON,PHI,EPHI,SC,V1C,V1CA,V1CB,V2C) &
    !$omp  private(V2CA,V2CB,V2CAB) &
    !$omp  reduction(+:SGCX,SGCC,ID)
    DO i=1,kkk
       rhoa=MAX(rhoe(i,1),small)
       rhob=MAX(rhoe(i,2),small)
       rho=rhoa+rhob
       IF (rho.GT.3._real_8*small) THEN
          x=rhoa**(0.33333333333333333_real_8)
          rsa=1._real_8/x
          y=rhob**(0.33333333333333333_real_8)
          rsb=1._real_8/y
          z=rho**(0.33333333333333333_real_8)
          zon=1.0_real_8/z
          rhoon=1.0_real_8/rho
          zz=SQRT(z)
          rs=rsfac/z
          IF (rho.GT.2.25_real_8*cntr%gceps) THEN
             smoo=1._real_8
             dsmoo=0._real_8
          ELSEIF (rho.LT.0.25_real_8*cntr%gceps) THEN
             smoo=0._real_8
             dsmoo=0._real_8
          ELSE
             drho=(rho-0.25_real_8*cntr%gceps)*ogceps
             smoo=drho*drho*(0.75_real_8-0.25_real_8*drho)
             dsmoo=1.5_real_8*drho*ogceps*(1._real_8-0.5_real_8*drho)
          ENDIF
          grhoa=SQRT(grad(i,1))
          grhob=SQRT(grad(i,5))
          grhoa=MAX(grhoa,small)
          grhob=MAX(grhob,small)
          xsa  = grhoa*(rsa*rsa*rsa*rsa)
          xsb  = grhob*(rsb*rsb*rsb*rsb)
          IF (rhoa.LT.2._real_8*small) xsa=0._real_8
          IF (rhob.LT.2._real_8*small) xsb=0._real_8
          xsaa = xsa*xsa
          xsbb = xsb*xsb
          sa2a = SQRT(1._real_8+xsaa)
          sa2b = SQRT(1._real_8+xsbb)
          shma = LOG(xsa+sa2a)
          shmb = LOG(xsb+sa2b)
          dda  = 1._real_8+6._real_8*b1*xsa*shma
          ddb  = 1._real_8+6._real_8*b1*xsb*shmb
          ddaon=1.0_real_8/dda
          ddbon=1.0_real_8/ddb
          ddda = 6._real_8*b1*(shma+xsa/sa2a)
          dddb = 6._real_8*b1*(shmb+xsb/sa2b)
          gfa  = -b1*xsaa*ddaon
          gfb  = -b1*xsbb*ddbon
          dgfa = (-2._real_8*b1*xsa*dda+b1*xsaa*ddda)*ddaon*ddaon
          dgfb = (-2._real_8*b1*xsb*ddb+b1*xsbb*dddb)*ddbon*ddbon
          sxa  = gfa*x*rhoa
          sxb  = gfb*y*rhob
          v1xa = 1.33333333333333_real_8*x*(gfa-xsa*dgfa)
          v1xb = 1.33333333333333_real_8*y*(gfb-xsb*dgfb)
          v2xa = dgfa/grhoa
          v2xb = dgfb/grhob
          ! ..perdew 86 correlation functional
          grho=grad(i,1)+grad(i,5)+2._real_8*grad(i,2)*grad(i,6)+&
               2._real_8*grad(i,3)*grad(i,7)+2._real_8*grad(i,4)*grad(i,8)
          rs2 = rs*rs
          rs3 = rs2*rs
          cna = pc2+p1*rs+p2*rs2
          cnb = 1._real_8+p3*rs+p4*rs2+1.e4_real_8*p2*rs3
          cnbon=1.0_real_8/cnb
          cn  = pc1 + cna*cnbon
          drs = dfx/(z*rho)
          dcna= (p1+2._real_8*p2*rs)*drs
          dcnb= (p3+2._real_8*p4*rs+3.e4_real_8*p2*rs2)*drs
          dcn = dcna*cnbon - cna*cnbon*cnbon*dcnb
          s1  = SQRT(x*x*rhoa+y*y*rhob)
          s1on=1.0_real_8/s1
          s2  = two13*zz*rhoon
          d   = s1*s2
          don=1.0_real_8/d
          dda = s2*0.83333333333333_real_8*(x*x*s1on-s1*rhoon)
          ddb = s2*0.83333333333333_real_8*(y*y*s1on-s1*rhoon)
          phi = 0.192_real_8*pci*SQRT(grho)/(zz*rho*cn)
          ephi= EXP(-phi)
          sc  = grho*zon*rhoon*cn*ephi*don
          v1c = sc*((1._real_8+phi)*dcn/cn-(1.33333333333333_real_8-&
               1.1666666666666_real_8*phi)*rhoon)
          v1ca= -sc*don*dda+v1c
          v1cb= -sc*don*ddb+v1c
          v2c = cn*ephi*zon*rhoon*(2._real_8-phi)*don
          v2ca= v2c
          v2cb= v2c
          v2cab= v2c
          ! ..sum up terms
          v1(i)=dsmoo*(sxa+sc)+smoo*(v1xa+v1ca)
          v2(i)=dsmoo*(sxb+sc)+smoo*(v1xb+v1cb)
          vtmp(i,1)=smoo*(v2xa+v2ca)
          vtmp(i,2)=smoo*(v2xb+v2cb)
          vtmp(i,3)=smoo*v2cab
          sgcx = sgcx + (smoo*(sxa+sxb))
          sgcc = sgcc + (smoo*sc)
          id = id + 1
       ENDIF
    ENDDO
    ! ..counting floating point operations
    fadd   = id *  39._real_8
    fmult  = id * 105._real_8
    fdiv   = id *  29._real_8
    fspez  = id *  11._real_8
    flops=fadd+fmult+fdiv+fspez
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gcsxp86
  ! ==================================================================
  SUBROUTINE gcsxlyp(b1,sgcx,sgcc,rhoe,v1,v2,vtmp,grad,flops)
    REAL(real_8)                             :: b1, sgcx, sgcc
    REAL(real_8), DIMENSION(fpar%kr1*fpar%&
      kr2s*fpar%kr3s, 2)                     :: rhoe
    COMPLEX(real_8) :: v1(fpar%kr1*fpar%kr2s*fpar%kr3s), &
      v2(fpar%kr1*fpar%kr2s*fpar%kr3s)
    REAL(real_8), DIMENSION(fpar%kr1*fpar%&
      kr2s*fpar%kr3s, 3)                     :: vtmp
    REAL(real_8), DIMENSION(fpar%kr1*fpar%&
      kr2s*fpar%kr3s, 8)                     :: grad
    REAL(real_8)                             :: flops

    REAL(real_8), PARAMETER :: a = 0.04918_real_8, b = 0.132_real_8, &
      c = 0.2533_real_8, d = 0.349_real_8 , eps2 = 1.e-20_real_8, &
      small = 1.e-24_real_8 

    INTEGER                                  :: i, id, kkk
    REAL(real_8) :: dda, ddaon, ddb, ddbon, ddda, dddb, dder, der, dgfa, &
      dgfb, dlaa, dlaaa, dlaab, dlab, dlaba, dlabb, dlbb, dlbba, dlbbb, dor, &
      dr, drho, dron, dsmoo, fadd, fdiv, fmult, fspez, gfa, gfb, grhoa, &
      grhoaa, grhoab, grhob, grhobb, ogceps, or, oron, rho, rhoa, rhob, &
      rhoon, rs, rsa, rsb, sa2a, sa2b, sc, shma, shmb, smoo, sxa, sxb, two13, &
      v1ca, v1cb, v1xa, v1xb, v2ca, v2cab, v2cb, v2xa, v2xb, x, xsa, xsaa, &
      xsb, xsbb, y, z

! ==--------------------------------------------------------------==

    two13=2._real_8**0.333333333333333_real_8
    ogceps=1._real_8/cntr%gceps
    sgcx=0._real_8
    sgcc=0._real_8
    id=0
    kkk=fpar%kr1*fpar%kr2s*fpar%kr3s
    !$omp parallel do private(I) &
    !$omp  private(RHOA,RHOB,RHO,RHOON,X,RSA,Y,RSB,Z,RS) &
    !$omp  private(SMOO,DSMOO,DRHO,GRHOAA,GRHOAB,GRHOBB,GRHOA) &
    !$omp  private(GRHOB,XSA,XSB,XSAA,XSBB,SA2A,SA2B,SHMA,SHMB) &
    !$omp  private(DDA,DDB,DDAON,DDBON,DDDA,DDDB,GFA,GFB,DGFA) &
    !$omp  private(DGFB,SXA,SXB,V1XA,V1XB,V2XA,V2XB,DR,DRON,OR) &
    !$omp  private(ORON,DOR,DER,DDER,DLAA,DLAB,DLBB,DLAAA,DLAAB) &
    !$omp  private(DLABA,DLABB,DLBBA,DLBBB,SC,V1CA,V1CB,V2CA,V2CB) &
    !$omp  private(V2CAB) &
    !$omp  reduction(+:SGCX,SGCC,ID)
    DO i=1,kkk
       rhoa=MAX(rhoe(i,1),small)
       rhob=MAX(rhoe(i,2),small)
       rho=rhoa+rhob
       rhoon=1.0_real_8/rho
       IF (rho.GT.3._real_8*small) THEN
          x=rhoa**(0.33333333333333333_real_8)
          rsa=1._real_8/x
          y=rhob**(0.33333333333333333_real_8)
          rsb=1._real_8/y
          z=rho**(0.33333333333333333_real_8)
          rs=1._real_8/z
          IF (rho.GT.2.25_real_8*cntr%gceps) THEN
             smoo=1._real_8
             dsmoo=0._real_8
          ELSEIF (rho.LT.0.25_real_8*cntr%gceps) THEN
             smoo=0._real_8
             dsmoo=0._real_8
          ELSE
             drho=(rho-0.25_real_8*cntr%gceps)*ogceps
             smoo=drho*drho*(0.75_real_8-0.25_real_8*drho)
             dsmoo=1.5_real_8*drho*ogceps*(1._real_8-0.5_real_8*drho)
          ENDIF
          grhoaa=grad(i,1)
          grhoab=grad(i,2)*grad(i,6)+grad(i,3)*grad(i,7)+&
               grad(i,4)*grad(i,8)
          grhobb=grad(i,5)
          grhoa=SQRT(grad(i,1))
          grhob=SQRT(grad(i,5))
          grhoa=MAX(grhoa,small)
          grhob=MAX(grhob,small)
          xsa  = grhoa*(rsa*rsa*rsa*rsa)
          xsb  = grhob*(rsb*rsb*rsb*rsb)
          IF (rhoa.LT.2._real_8*small) xsa=0._real_8
          IF (rhob.LT.2._real_8*small) xsb=0._real_8
          xsaa = xsa*xsa
          xsbb = xsb*xsb
          sa2a = SQRT(1._real_8+xsaa)
          sa2b = SQRT(1._real_8+xsbb)
          shma = LOG(xsa+sa2a)
          shmb = LOG(xsb+sa2b)
          dda  = 1._real_8+6._real_8*b1*xsa*shma
          ddb  = 1._real_8+6._real_8*b1*xsb*shmb
          ddaon=1.0_real_8/dda
          ddbon=1.0_real_8/ddb
          ddda = 6._real_8*b1*(shma+xsa/sa2a)
          dddb = 6._real_8*b1*(shmb+xsb/sa2b)
          gfa  = -b1*xsaa*ddaon
          gfb  = -b1*xsbb*ddbon
          dgfa = (-2._real_8*b1*xsa*dda+b1*xsaa*ddda)*ddaon*ddaon
          dgfb = (-2._real_8*b1*xsb*ddb+b1*xsbb*dddb)*ddbon*ddbon
          sxa  = gfa*x*rhoa
          sxb  = gfb*y*rhob
          v1xa = 1.33333333333333_real_8*x*(gfa-xsa*dgfa)
          v1xb = 1.33333333333333_real_8*y*(gfb-xsb*dgfb)
          v2xa = dgfa/grhoa
          v2xb = dgfb/grhob
          ! ..LYP correlation functional
          dr=(1._real_8+d*rs)
          dron=1.0_real_8/dr
          or=EXP(-c*rs)*dron*(rs**11)
          or=MAX(or,small)
          oron=1.0_real_8/or
          dor=-0.3333333333333_real_8*(rs**4)*or*(11._real_8*z-c-d*dron)
          der=c*rs+d*rs*dron
          dder=0.3333333333333_real_8*(d*d*(rs**5)*dron*dron-der*rhoon)
          dlaa=-a*b*or*(0.1111111111111_real_8*rhoa*rhob*(1._real_8-&
               3._real_8*der-(der-11._real_8)*rhoa*rhoon)-rhob*rhob)
          dlab=-a*b*or*(0.1111111111111_real_8*rhoa*rhob*(47._real_8-&
               7._real_8*der)-1.3333333333333_real_8*rho*rho)
          dlbb=-a*b*or*(0.1111111111111_real_8*rhoa*rhob*(1._real_8-&
               3._real_8*der-(der-11._real_8)*rhob*rhoon)-rhoa*rhoa)
          dlaaa=dor*oron*dlaa-a*b*or*(0.1111111111111_real_8*rhob*&
               (1._real_8-3._real_8*der-(der-11._real_8)*rhoa*rhoon)-&
               0.1111111111111_real_8*rhoa*rhob*((3._real_8+rhoa*rhoon)*&
               dder+(der-11._real_8)*rhob*rhoon*rhoon))
          dlaab=dor*oron*dlaa-a*b*or*(0.1111111111111_real_8*rhoa*&
               (1._real_8-3._real_8*der-(der-11._real_8)*rhoa*rhoon)-&
               0.1111111111111_real_8*rhoa*rhob*((3._real_8+rhoa*rhoon)*&
               dder-(der-11._real_8)*rhoa*rhoon*rhoon)-2._real_8*rhob)
          dlaba=dor*oron*dlab-a*b*or*(0.1111111111111_real_8*rhob*&
               (47._real_8-7._real_8*der)-0.7777777777778_real_8*rhoa*rhob*dder-&
               2.66666666666667_real_8*rho)
          dlabb=dor*oron*dlab-a*b*or*(0.1111111111111_real_8*rhoa*&
               (47._real_8-7._real_8*der)-0.7777777777778_real_8*rhoa*rhob*dder-&
               2.66666666666667_real_8*rho)
          dlbba=dor*oron*dlbb-a*b*or*(0.1111111111111_real_8*rhob*&
               (1._real_8-3._real_8*der-(der-11._real_8)*rhob*rhoon)-&
               0.1111111111111_real_8*rhoa*rhob*((3._real_8+rhob*rhoon)*&
               dder-(der-11._real_8)*rhob*rhoon*rhoon)-2._real_8*rhoa)
          dlbbb=dor*oron*dlbb-a*b*or*(0.1111111111111_real_8*rhoa*&
               (1._real_8-3._real_8*der-(der-11._real_8)*rhob*rhoon)-&
               0.1111111111111_real_8*rhoa*rhob*((3._real_8+rhob*rhoon)*&
               dder+(der-11._real_8)*rhoa*rhoon*rhoon))
          sc=dlaa*grhoaa+dlab*grhoab+dlbb*grhobb
          v1ca=dlaaa*grhoaa+dlaba*grhoab+dlbba*grhobb
          v1cb=dlaab*grhoaa+dlabb*grhoab+dlbbb*grhobb
          v2ca=2._real_8*dlaa
          v2cb=2._real_8*dlbb
          v2cab=dlab
          ! ..sum up terms
          v1(i)=dsmoo*(sxa+sc)+smoo*(v1xa+v1ca)
          v2(i)=dsmoo*(sxb+sc)+smoo*(v1xb+v1cb)
          vtmp(i,1)=smoo*(v2xa+v2ca)
          vtmp(i,2)=smoo*(v2xb+v2cb)
          vtmp(i,3)=smoo*v2cab
          sgcx = sgcx + (smoo*(sxa+sxb))
          sgcc = sgcc + (smoo*sc)
          id = id + 1
       ENDIF
    ENDDO
    ! ..counting floating point operations
    fadd   = id *  77._real_8
    fmult  = id * 185._real_8
    fdiv   = id *  37._real_8
    fspez  = id *  13._real_8
    flops=fadd+fmult+fdiv+fspez
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gcsxlyp
  ! ==================================================================
  SUBROUTINE gcspbe(sgcx,sgcc,rhoe,v1,v2,vtmp,grad,flops)
    REAL(real_8)                             :: sgcx, sgcc
    REAL(real_8), DIMENSION(fpar%kr1*fpar%&
      kr2s*fpar%kr3s, 2)                     :: rhoe
    COMPLEX(real_8) :: v1(fpar%kr1*fpar%kr2s*fpar%kr3s), &
      v2(fpar%kr1*fpar%kr2s*fpar%kr3s)
    REAL(real_8), DIMENSION(fpar%kr1*fpar%&
      kr2s*fpar%kr3s, 3)                     :: vtmp
    REAL(real_8), DIMENSION(fpar%kr1*fpar%&
      kr2s*fpar%kr3s, 8)                     :: grad
    REAL(real_8)                             :: flops

    REAL(real_8), PARAMETER :: a0u = 0.4581652932831429_real_8, &
      a1u = 2.217058676663745_real_8, a2u = 0.7405551735357053_real_8, &
      a3u = 0.01968227878617998_real_8 , ax = -0.738558766382022406_real_8, &
      b1u = 1.0_real_8, b2u = 4.504130959426697_real_8, &
      b3u = 1.110667363742916_real_8, b4u = 0.02359291751427506_real_8 , &
      be = 0.06672455060314922_real_8, da0 = 0.119086804055547_real_8, &
      da1 = 0.6157402568883345_real_8, da2 = 0.1574201515892867_real_8, &
      da3 = 0.003532336663397157_real_8 , db1 = 0.0_real_8, &
      db2 = 0.2673612973836267_real_8, db3 = 0.2052004607777787_real_8, &
      db4 = 0.004200005045691381_real_8 
    REAL(real_8), PARAMETER :: eps2 = 1.e-20_real_8, &
      f1 = -1.39578858466194911_real_8 , ga = 0.031090690869654895_real_8 , &
      rsfac = 0.6203504908994000_real_8 , small = 1.e-24_real_8, &
      ssmall = 1.e-12_real_8 , t1 = 0.854960467080683_real_8, &
      t2 = 0.07916300621117439_real_8, t3 = 0.02580127609845684_real_8, &
      t4 = 0.0121839359353824_real_8, t5 = 0.006919272259599879_real_8, &
      t6 = 0.004391524649611372_real_8, t7 = 0.003002751897170168_real_8, &
      t8 = 0.00216587382212552_real_8, t9 = 0.01141189204579585_real_8 , &
      uk = 0.8040_real_8, um = 0.2195149727645171_real_8, ul = um/uk , &
      us = 0.161620459673995492_real_8

    INTEGER                                  :: i, id, kkk
    REAL(real_8) :: a0p, a1p, a2p, a3p, aa, af, as, b1p, b2p, b3p, b4p, beon, &
      bot, dadra, dadrb, dbot, dedra, dedrb, dfx, dfxs, dhdra, dhdrb, dhdt, &
      dphida, dphidb, dphide, drho, dsda, dsdra, dsdrb, dsdt, dsmoo, dtdphi, &
      dtdra, dtdrb, dtop, dvxc, dxa, dxb, epsxcu, esla, et, eta, eta1ob3, ex, &
      expe, f1s, f1u, fadd, fdiv, fmult, fspez, fx, fxs, gaon, grho, grhoa, &
      grhob, grhoon, h0, ob3, obot, ogceps, phi, phi3, phi3on, phion, pion, &
      po, rho, rhoa, rhob, rhoon, rr, rs, rsa, rsb, s1, s1on, s2, sc, smoo, &
      sxa, sxb, t13, t43, top, tt, v1ca, v1cb, v1xa, v1xb, v2ca, v2cab, v2cb, &
      v2xa, v2xb, vsla, vslb, vxc
    REAL(real_8) :: x, xbot, xkf, xks, xtop, xy, y, yt, yt3, z, zz

! ==--------------------------------------------------------------==

    gaon=1.0_real_8/ga
    beon=1.0_real_8/be
    t43=2.0_real_8**(-1.33333333333333_real_8)
    t13=2.0_real_8**0.33333333333333_real_8
    ob3=0.33333333333333_real_8
    f1s=f1*0.66666666666666_real_8
    f1u=f1s*1.3333333333333333_real_8
    ogceps=1._real_8/cntr%gceps
    sgcx=0._real_8
    sgcc=0._real_8
    id=0
    kkk=fpar%kr1*fpar%kr2s*fpar%kr3s
    !$omp parallel do private(I) &
    !$omp  private(RHOA,RHOB,RHO,RHOON,X,Y,Z,RS,ETA,ET,FXS) &
    !$omp  private(DFXS,A0P,A1P,A2P,A3P,B2P,B3P,B4P,TOP,DTOP) &
    !$omp  private(XTOP,BOT,DBOT,XBOT,OBOT,EPSXCU,VXC,DVXC) &
    !$omp  private(DXA,DXB,ESLA,VSLA,VSLB) &
    !$omp  private(SMOO,DSMOO,DRHO,B1P)
    DO i=1,kkk
       rhoa=MAX(rhoe(i,1),small)
       rhob=MAX(rhoe(i,2),small)
       rho=rhoa+rhob
       rhoon=1.0_real_8/rho
       IF (rho.GT.3._real_8*small) THEN
          x=rhoa**(0.33333333333333333_real_8)
          y=rhob**(0.33333333333333333_real_8)
          z=rho**(0.33333333333333333_real_8)
          rs=rsfac/z
          IF (rho.GT.2.25_real_8*cntr%gceps) THEN
             smoo=1._real_8
             dsmoo=0._real_8
          ELSEIF (rho.LT.0.25_real_8*cntr%gceps) THEN
             smoo=0._real_8
             dsmoo=0._real_8
          ELSE
             drho=(rho-0.25_real_8*cntr%gceps)*ogceps
             smoo=drho*drho*(0.75_real_8-0.25_real_8*drho)
             dsmoo=1.5_real_8*drho*ogceps*(1._real_8-0.5_real_8*drho)
          ENDIF
          ! ..LDA correlation
          eta=(rho-2._real_8*rhob)*rhoon
          et=eta*eta
          fxs=et*(t1+et*(t2+et*(t3+et*(t4+et*(t5+et*&
               (t6+et*(t7+et*(t8+et*t9))))))))
          dfxs=(2._real_8*eta)*(t1+et*(2._real_8*t2+et*(3._real_8*t3+et*&
               (4._real_8*t4+et*(5._real_8*t5+et*(6._real_8*t6+et*(7._real_8*t7+et*&
               (8._real_8*t8+et*9._real_8*t9))))))))
          ! 
          a0p=a0u+fxs*da0
          a1p=a1u+fxs*da1
          a2p=a2u+fxs*da2
          a3p=a3u+fxs*da3
          b1p=b1u+fxs*db1
          b2p=b2u+fxs*db2
          b3p=b3u+fxs*db3
          b4p=b4u+fxs*db4
          ! 
          top=a0p+rs*(a1p+rs*(a2p+rs*a3p))
          dtop=a1p+rs*(2._real_8*a2p+rs*3._real_8*a3p)
          xtop=da0+rs*(da1+rs*(da2+rs*da3))
          bot=rs*(b1p+rs*(b2p+rs*(b3p+rs*b4p)))
          dbot=b1p+rs*(2._real_8*b2p+rs*(3._real_8*b3p+rs*4._real_8*b4p))
          xbot=rs*(db1+rs*(db2+rs*(db3+rs*db4)))
          obot=1._real_8/bot
          epsxcu=-top*obot
          vxc=epsxcu+rs*(dtop*obot-dbot*(top*obot*obot))*&
               0.333333333333333_real_8
          dvxc=-(xtop*obot-xbot*(top*obot*obot))
          dxa=vxc+dvxc*dfxs*(1._real_8-eta)
          dxb=vxc-dvxc*dfxs*(1._real_8+eta)
          esla=f1s*(x*rhoa+y*rhob)*rhoon
          vsla=f1u*x
          vslb=f1u*y
          epsxcu=epsxcu-esla
          dxa=dxa-vsla
          dxb=dxb-vslb
          ! 
          vtmp(i,1)=epsxcu
          vtmp(i,2)=dxa
          vtmp(i,3)=dxb
       ENDIF
    ENDDO
    !$omp parallel do private(I) &
    !$omp  private(RSA,RSB,SMOO,DSMOO,DRHO,GRHOA,GRHOB,AA,RR) &
    !$omp  private(EX,S2,PO,FX,SXA,DFX,V1XA,V2XA,SXB,V1XB,V2XB) &
    !$omp  private(GRHO,GRHOON,AS,PHI,PHION,PHI3,PHI3ON,PION,XKF) &
    !$omp  private(XKS,TT,EXPE,AF,YT,XY,S1,S1ON,H0,SC,ETA1OB3) &
    !$omp  private(DPHIDE,DEDRA,DEDRB,DPHIDA,DPHIDB,DTDRA,DTDRB) &
    !$omp  private(DADRA,DADRB,YT3,DSDA,DSDT,DSDRA,DSDRB,DHDT) &
    !$omp  private(DHDRA,DHDRB,V1CA,V2CA,V1CB,V2CB,V2CAB) &
    !$omp  private(RHOA,RHOB,RHO,RHOON,X,Y,Z,ZZ,RS) &
    !$omp  private(ETA,ET,EPSXCU,DXA,DXB,DTDPHI) &
    !$omp  reduction(+:SGCX,SGCC,ID)
    DO i=1,kkk
       rhoa=MAX(rhoe(i,1),small)
       rhob=MAX(rhoe(i,2),small)
       rho=rhoa+rhob
       rhoon=1.0_real_8/rho
       IF (rho.GT.3._real_8*small) THEN
          x=rhoa**(0.33333333333333333_real_8)
          rsa=1._real_8/x
          y=rhob**(0.33333333333333333_real_8)
          rsb=1._real_8/y
          z=rho**(0.33333333333333333_real_8)
          zz=SQRT(z)
          rs=rsfac/z
          IF (rho.GT.2.25_real_8*cntr%gceps) THEN
             smoo=1._real_8
             dsmoo=0._real_8
          ELSEIF (rho.LT.0.25_real_8*cntr%gceps) THEN
             smoo=0._real_8
             dsmoo=0._real_8
          ELSE
             drho=(rho-0.25_real_8*cntr%gceps)*ogceps
             smoo=drho*drho*(0.75_real_8-0.25_real_8*drho)
             dsmoo=1.5_real_8*drho*ogceps*(1._real_8-0.5_real_8*drho)
          ENDIF
          eta=(rho-2._real_8*rhob)*rhoon
          et=eta*eta

          IF (rhoa.GT.small) THEN
             grhoa=MAX(grad(i,1),small)
          ELSE
             grhoa=small
          ENDIF
          IF (rhob.GT.small) THEN
             grhob=MAX(grad(i,5),small)
          ELSE
             grhob=small
          ENDIF
          ! 
          epsxcu=vtmp(i,1)
          dxa=vtmp(i,2)
          dxb=vtmp(i,3)
          ! ..exchange
          aa    = 4._real_8*grhoa
          rr    = t43*rsa*rsa*rsa*rsa
          ex    = ax/rr
          s2    = aa*rr*rr*us*us
          po    = 1._real_8/(1._real_8 + ul*s2)
          fx    = uk-uk*po
          sxa   = 0.5_real_8*ex*fx
          dfx   = 2._real_8*uk*ul*po*po
          v1xa  = t13*1.33333333333333_real_8*ax*x*(fx-s2*dfx)
          v2xa  = 2._real_8*ex*dfx*(us*rr)**2
          aa    = 4._real_8*grhob
          rr    = t43*rsb*rsb*rsb*rsb
          ex    = ax/rr
          s2    = aa*rr*rr*us*us
          po    = 1._real_8/(1._real_8 + ul*s2)
          fx    = uk-uk*po
          sxb   = 0.5_real_8*ex*fx
          dfx   = 2._real_8*uk*ul*po*po
          v1xb  = t13*1.33333333333333_real_8*ax*y*(fx-s2*dfx)
          v2xb  = 2._real_8*ex*dfx*(us*rr)**2
          ! ..correlation
          grho=grad(i,1)+grad(i,5)+2._real_8*grad(i,2)*grad(i,6)+&
               2._real_8*grad(i,3)*grad(i,7)+2._real_8*grad(i,4)*grad(i,8)
          grhoon=1.0_real_8/grho
          as=SQRT(grho)
          phi=0.5_real_8*((1._real_8+eta)**0.66666666666666_real_8+(1._real_8-eta)&
               **0.66666666666666_real_8)
          phion=1.0_real_8/phi
          phi3=phi*phi*phi
          phi3on=phion*phion*phion
          pion=1.0_real_8/pi
          rs    = (3.0_real_8*0.25_real_8*pion*rhoon)**ob3
          xkf   = (2.25_real_8*pi)**ob3/rs
          xks   = SQRT(4._real_8*xkf*pion)
          tt    = as/(2._real_8*xks*rho*phi)
          expe  = EXP(-epsxcu*phi3on*gaon)
          af    = be*gaon * (1._real_8/(expe-1._real_8))
          yt    = af*tt*tt
          xy    = (1._real_8+yt)/(1._real_8+yt+yt*yt)
          s1    = 1._real_8+be/ga*tt*tt*xy
          s1on=1.0_real_8/s1
          h0    = ga*phi3*LOG(s1)
          sc    = rho*h0
          ! 
          dtdphi = -tt*phion
          eta1ob3= 1._real_8/(3._real_8*(1._real_8-eta)**ob3+ssmall)
          dphide = 1._real_8*eta1ob3*eta1ob3
          dedra  = 2._real_8*rhob*rhoon*rhoon
          dedrb  = -2._real_8*rhoa*rhoon*rhoon
          dphida = dphide*dedra
          dphidb = dphide*dedrb
          dtdra  = -tt*(dphida*phion+1.16666666666666_real_8*rhoon)
          dtdrb  = -tt*(dphidb*phion+1.16666666666666_real_8*rhoon)
          dadra  = af*af*expe*(-beon*phi3on)*(3._real_8*epsxcu*phion*&
               dphida-(dxa-epsxcu)*rhoon)
          dadrb  = af*af*expe*(-beon*phi3on)*(3._real_8*epsxcu*phion*&
               dphidb-(dxb-epsxcu)*rhoon)
          yt3    = 1.0_real_8/(1._real_8+yt+yt*yt)**2
          dsda   = -be*gaon*af*tt**6*(2._real_8+yt)*yt3
          dsdt   = 2._real_8*be*gaon*tt*(1._real_8+2._real_8*yt)*yt3
          dsdra  = dsda*dadra + dsdt*dtdra
          dsdrb  = dsda*dadrb + dsdt*dtdrb
          dhdt   = ga*phi3*s1on*dsdt
          dhdra  = 3._real_8*h0*phion*dphida + ga*phi3*s1on*dsdra
          dhdrb  = 3._real_8*h0*phion*dphidb + ga*phi3*s1on*dsdrb
          ! ..
          v1ca  = h0 + rho*dhdra
          v2ca  = rho*dhdt*tt*grhoon
          v1cb  = h0 + rho*dhdrb
          v2cb  = rho*dhdt*tt*grhoon
          v2cab = rho*dhdt*tt*grhoon
          ! ..sum up terms
          IF (rhoa.GT.small) THEN
             v1(i)=dsmoo*(sxa+sc)+smoo*(v1xa+v1ca)
             vtmp(i,1)=smoo*(v2xa+v2ca)
          ENDIF
          IF (rhob.GT.small) THEN
             v2(i)=dsmoo*(sxb+sc)+smoo*(v1xb+v1cb)
             vtmp(i,2)=smoo*(v2xb+v2cb)
          ENDIF
          IF ((rhoa.GT.small).AND.(rhob.GT.small)) THEN
             vtmp(i,3)=smoo*v2cab
          ENDIF
          id = id + 1
          sgcx = sgcx + (smoo*(sxa+sxb))
          sgcc = sgcc + (smoo*sc)
       ENDIF
    ENDDO
    ! ..counting floating point operations
    fadd   = id * 100._real_8
    fmult  = id * 227._real_8
    fdiv   = id *  42._real_8
    fspez  = id *  12._real_8
    flops=fadd+fmult+fdiv+fspez
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gcspbe
  ! ==================================================================
  SUBROUTINE gcsrevpbe(sgcx,sgcc,rhoe,v1,v2,vtmp,grad,flops)
    REAL(real_8)                             :: sgcx, sgcc
    REAL(real_8), DIMENSION(fpar%kr1*fpar%&
      kr2s*fpar%kr3s, 2)                     :: rhoe
    COMPLEX(real_8) :: v1(fpar%kr1*fpar%kr2s*fpar%kr3s), &
      v2(fpar%kr1*fpar%kr2s*fpar%kr3s)
    REAL(real_8), DIMENSION(fpar%kr1*fpar%&
      kr2s*fpar%kr3s, 3)                     :: vtmp
    REAL(real_8), DIMENSION(fpar%kr1*fpar%&
      kr2s*fpar%kr3s, 8)                     :: grad
    REAL(real_8)                             :: flops

    REAL(real_8), PARAMETER :: a0u = 0.4581652932831429_real_8, &
      a1u = 2.217058676663745_real_8, a2u = 0.7405551735357053_real_8, &
      a3u = 0.01968227878617998_real_8 , ax = -0.738558766382022406_real_8, &
      b1u = 1.0_real_8, b2u = 4.504130959426697_real_8, &
      b3u = 1.110667363742916_real_8, b4u = 0.02359291751427506_real_8 , &
      be = 0.06672455060314922_real_8, da0 = 0.119086804055547_real_8, &
      da1 = 0.6157402568883345_real_8, da2 = 0.1574201515892867_real_8, &
      da3 = 0.003532336663397157_real_8 , db1 = 0.0_real_8, &
      db2 = 0.2673612973836267_real_8, db3 = 0.2052004607777787_real_8, &
      db4 = 0.004200005045691381_real_8 
    REAL(real_8), PARAMETER :: eps2 = 1.e-20_real_8, &
      f1 = -1.39578858466194911_real_8 , ga = 0.031090690869654895_real_8 , &
      rsfac = 0.6203504908994000_real_8 , small = 1.e-24_real_8, &
      ssmall = 1.e-12_real_8 , t1 = 0.854960467080683_real_8, &
      t2 = 0.07916300621117439_real_8, t3 = 0.02580127609845684_real_8, &
      t4 = 0.0121839359353824_real_8, t5 = 0.006919272259599879_real_8, &
      t6 = 0.004391524649611372_real_8, t7 = 0.003002751897170168_real_8, &
      t8 = 0.00216587382212552_real_8, t9 = 0.01141189204579585_real_8 , &
      uk = 1.2450_real_8, um = 0.2195149727645171_real_8, ul = um/uk , &
      us = 0.161620459673995492_real_8

    INTEGER                                  :: i, id, kkk
    REAL(real_8) :: a0p, a1p, a2p, a3p, aa, af, as, b1p, b2p, b3p, b4p, beon, &
      bot, dadra, dadrb, dbot, dedra, dedrb, dfx, dfxs, dhdra, dhdrb, dhdt, &
      dphida, dphidb, dphide, drho, dsda, dsdra, dsdrb, dsdt, dsmoo, dtdphi, &
      dtdra, dtdrb, dtop, dvxc, dxa, dxb, epsxcu, esla, et, eta, eta1ob3, ex, &
      expe, f1s, f1u, fadd, fdiv, fmult, fspez, fx, fxs, gaon, grho, grhoa, &
      grhob, grhoon, h0, ob3, obot, ogceps, phi, phi3, phi3on, phion, pion, &
      po, rho, rhoa, rhob, rhoon, rr, rs, rsa, rsb, s1, s1on, s2, sc, smoo, &
      sxa, sxb, t13, t43, top, tt, v1ca, v1cb, v1xa, v1xb, v2ca, v2cab, v2cb, &
      v2xa, v2xb, vsla, vslb, vxc
    REAL(real_8) :: x, xbot, xkf, xks, xtop, xy, y, yt, yt3, z, zz

! ==--------------------------------------------------------------==

    gaon=1.0_real_8/ga
    beon=1.0_real_8/be
    t43=2.0_real_8**(-1.33333333333333_real_8)
    t13=2.0_real_8**0.33333333333333_real_8
    ob3=0.33333333333333_real_8
    f1s=f1*0.66666666666666_real_8
    f1u=f1s*1.3333333333333333_real_8
    ogceps=1._real_8/cntr%gceps
    sgcx=0._real_8
    sgcc=0._real_8
    id=0
    kkk=fpar%kr1*fpar%kr2s*fpar%kr3s
    !$omp parallel do private(I) &
    !$omp  private(RHOA,RHOB,RHO,RHOON,X,Y,Z,RS,ETA,ET,FXS) &
    !$omp  private(DFXS,A0P,A1P,A2P,A3P,B2P,B3P,B4P,TOP,DTOP) &
    !$omp  private(XTOP,BOT,DBOT,XBOT,OBOT,EPSXCU,VXC,DVXC) &
    !$omp  private(DXA,DXB,ESLA,VSLA,VSLB) &
    !$omp  private(SMOO,DSMOO,DRHO,B1P)
    DO i=1,kkk
       rhoa=MAX(rhoe(i,1),small)
       rhob=MAX(rhoe(i,2),small)
       rho=rhoa+rhob
       rhoon=1.0_real_8/rho
       IF (rho.GT.3._real_8*small) THEN
          x=rhoa**(0.33333333333333333_real_8)
          y=rhob**(0.33333333333333333_real_8)
          z=rho**(0.33333333333333333_real_8)
          rs=rsfac/z
          IF (rho.GT.2.25_real_8*cntr%gceps) THEN
             smoo=1._real_8
             dsmoo=0._real_8
          ELSEIF (rho.LT.0.25_real_8*cntr%gceps) THEN
             smoo=0._real_8
             dsmoo=0._real_8
          ELSE
             drho=(rho-0.25_real_8*cntr%gceps)*ogceps
             smoo=drho*drho*(0.75_real_8-0.25_real_8*drho)
             dsmoo=1.5_real_8*drho*ogceps*(1._real_8-0.5_real_8*drho)
          ENDIF
          ! ..LDA correlation
          eta=(rho-2._real_8*rhob)*rhoon
          et=eta*eta
          fxs=et*(t1+et*(t2+et*(t3+et*(t4+et*(t5+et*&
               (t6+et*(t7+et*(t8+et*t9))))))))
          dfxs=(2._real_8*eta)*(t1+et*(2._real_8*t2+et*(3._real_8*t3+et*&
               (4._real_8*t4+et*(5._real_8*t5+et*(6._real_8*t6+et*(7._real_8*t7+et*&
               (8._real_8*t8+et*9._real_8*t9))))))))
          ! 
          a0p=a0u+fxs*da0
          a1p=a1u+fxs*da1
          a2p=a2u+fxs*da2
          a3p=a3u+fxs*da3
          b1p=b1u+fxs*db1
          b2p=b2u+fxs*db2
          b3p=b3u+fxs*db3
          b4p=b4u+fxs*db4
          ! 
          top=a0p+rs*(a1p+rs*(a2p+rs*a3p))
          dtop=a1p+rs*(2._real_8*a2p+rs*3._real_8*a3p)
          xtop=da0+rs*(da1+rs*(da2+rs*da3))
          bot=rs*(b1p+rs*(b2p+rs*(b3p+rs*b4p)))
          dbot=b1p+rs*(2._real_8*b2p+rs*(3._real_8*b3p+rs*4._real_8*b4p))
          xbot=rs*(db1+rs*(db2+rs*(db3+rs*db4)))
          obot=1._real_8/bot
          epsxcu=-top*obot
          vxc=epsxcu+rs*(dtop*obot-dbot*(top*obot*obot))*&
               0.333333333333333_real_8
          dvxc=-(xtop*obot-xbot*(top*obot*obot))
          dxa=vxc+dvxc*dfxs*(1._real_8-eta)
          dxb=vxc-dvxc*dfxs*(1._real_8+eta)
          esla=f1s*(x*rhoa+y*rhob)*rhoon
          vsla=f1u*x
          vslb=f1u*y
          epsxcu=epsxcu-esla
          dxa=dxa-vsla
          dxb=dxb-vslb
          ! 
          vtmp(i,1)=epsxcu
          vtmp(i,2)=dxa
          vtmp(i,3)=dxb
       ENDIF
    ENDDO
    !$omp parallel do private(I) &
    !$omp  private(RSA,RSB,SMOO,DSMOO,DRHO,GRHOA,GRHOB,AA,RR) &
    !$omp  private(EX,S2,PO,FX,SXA,DFX,V1XA,V2XA,SXB,V1XB,V2XB) &
    !$omp  private(GRHO,GRHOON,AS,PHI,PHION,PHI3,PHI3ON,PION,XKF) &
    !$omp  private(XKS,TT,EXPE,AF,YT,XY,S1,S1ON,H0,SC,ETA1OB3) &
    !$omp  private(DPHIDE,DEDRA,DEDRB,DPHIDA,DPHIDB,DTDRA,DTDRB) &
    !$omp  private(DADRA,DADRB,YT3,DSDA,DSDT,DSDRA,DSDRB,DHDT) &
    !$omp  private(DHDRA,DHDRB,V1CA,V2CA,V1CB,V2CB,V2CAB) &
    !$omp  private(RHOA,RHOB,RHO,RHOON,X,Y,Z,ZZ,RS) &
    !$omp  private(ETA,ET,EPSXCU,DXA,DXB,DTDPHI) &
    !$omp  reduction(+:SGCX,SGCC,ID)
    DO i=1,kkk
       rhoa=MAX(rhoe(i,1),small)
       rhob=MAX(rhoe(i,2),small)
       rho=rhoa+rhob
       rhoon=1.0_real_8/rho
       IF (rho.GT.3._real_8*small) THEN
          x=rhoa**(0.33333333333333333_real_8)
          rsa=1._real_8/x
          y=rhob**(0.33333333333333333_real_8)
          rsb=1._real_8/y
          z=rho**(0.33333333333333333_real_8)
          zz=SQRT(z)
          rs=rsfac/z
          IF (rho.GT.2.25_real_8*cntr%gceps) THEN
             smoo=1._real_8
             dsmoo=0._real_8
          ELSEIF (rho.LT.0.25_real_8*cntr%gceps) THEN
             smoo=0._real_8
             dsmoo=0._real_8
          ELSE
             drho=(rho-0.25_real_8*cntr%gceps)*ogceps
             smoo=drho*drho*(0.75_real_8-0.25_real_8*drho)
             dsmoo=1.5_real_8*drho*ogceps*(1._real_8-0.5_real_8*drho)
          ENDIF
          eta=(rho-2._real_8*rhob)*rhoon
          et=eta*eta
          grhoa=MAX(grad(i,1),small)
          grhob=MAX(grad(i,5),small)
          ! 
          epsxcu=vtmp(i,1)
          dxa=vtmp(i,2)
          dxb=vtmp(i,3)
          ! ..exchange
          aa    = 4._real_8*grhoa
          rr    = t43*rsa*rsa*rsa*rsa
          ex    = ax/rr
          s2    = aa*rr*rr*us*us
          po    = 1._real_8/(1._real_8 + ul*s2)
          fx    = uk-uk*po
          sxa   = 0.5_real_8*ex*fx
          dfx   = 2._real_8*uk*ul*po*po
          v1xa  = t13*1.33333333333333_real_8*ax*x*(fx-s2*dfx)
          v2xa  = 2._real_8*ex*dfx*(us*rr)**2
          aa    = 4._real_8*grhob
          rr    = t43*rsb*rsb*rsb*rsb
          ex    = ax/rr
          s2    = aa*rr*rr*us*us
          po    = 1._real_8/(1._real_8 + ul*s2)
          fx    = uk-uk*po
          sxb   = 0.5_real_8*ex*fx
          dfx   = 2._real_8*uk*ul*po*po
          v1xb  = t13*1.33333333333333_real_8*ax*y*(fx-s2*dfx)
          v2xb  = 2._real_8*ex*dfx*(us*rr)**2
          ! ..correlation
          grho=grad(i,1)+grad(i,5)+2._real_8*grad(i,2)*grad(i,6)+&
               2._real_8*grad(i,3)*grad(i,7)+2._real_8*grad(i,4)*grad(i,8)
          grhoon=1.0_real_8/grho
          as=SQRT(grho)
          phi=0.5_real_8*((1._real_8+eta)**0.66666666666666_real_8+(1._real_8-eta)&
               **0.66666666666666_real_8)
          phion=1.0_real_8/phi
          phi3=phi*phi*phi
          phi3on=phion*phion*phion
          pion=1.0_real_8/pi
          rs    = (3.0_real_8*0.25_real_8*pion*rhoon)**ob3
          xkf   = (2.25_real_8*pi)**ob3/rs
          xks   = SQRT(4._real_8*xkf*pion)
          tt    = as/(2._real_8*xks*rho*phi)
          expe  = EXP(-epsxcu*phi3on*gaon)
          af    = be*gaon * (1._real_8/(expe-1._real_8))
          yt    = af*tt*tt
          xy    = (1._real_8+yt)/(1._real_8+yt+yt*yt)
          s1    = 1._real_8+be/ga*tt*tt*xy
          s1on=1.0_real_8/s1
          h0    = ga*phi3*LOG(s1)
          sc    = rho*h0
          ! 
          dtdphi = -tt*phion
          eta1ob3= 1._real_8/(3._real_8*(1._real_8-eta)**ob3+ssmall)
          dphide = 1._real_8*eta1ob3*eta1ob3
          dedra  = 2._real_8*rhob*rhoon*rhoon
          dedrb  = -2._real_8*rhoa*rhoon*rhoon
          dphida = dphide*dedra
          dphidb = dphide*dedrb
          dtdra  = -tt*(dphida*phion+1.16666666666666_real_8*rhoon)
          dtdrb  = -tt*(dphidb*phion+1.16666666666666_real_8*rhoon)
          dadra  = af*af*expe*(-beon*phi3on)*(3._real_8*epsxcu*phion*&
               dphida-(dxa-epsxcu)*rhoon)
          dadrb  = af*af*expe*(-beon*phi3on)*(3._real_8*epsxcu*phion*&
               dphidb-(dxb-epsxcu)*rhoon)
          yt3    = 1.0_real_8/(1._real_8+yt+yt*yt)**2
          dsda   = -be*gaon*af*tt**6*(2._real_8+yt)*yt3
          dsdt   = 2._real_8*be*gaon*tt*(1._real_8+2._real_8*yt)*yt3
          dsdra  = dsda*dadra + dsdt*dtdra
          dsdrb  = dsda*dadrb + dsdt*dtdrb
          dhdt   = ga*phi3*s1on*dsdt
          dhdra  = 3._real_8*h0*phion*dphida + ga*phi3*s1on*dsdra
          dhdrb  = 3._real_8*h0*phion*dphidb + ga*phi3*s1on*dsdrb
          ! ..
          v1ca  = h0 + rho*dhdra
          v2ca  = rho*dhdt*tt*grhoon
          v1cb  = h0 + rho*dhdrb
          v2cb  = rho*dhdt*tt*grhoon
          v2cab = rho*dhdt*tt*grhoon
          ! ..sum up terms

          IF (rhoa.GT.small) THEN
             v1(i)=dsmoo*(sxa+sc)+smoo*(v1xa+v1ca)
             vtmp(i,1)=smoo*(v2xa+v2ca)
          ENDIF
          IF (rhob.GT.small) THEN
             v2(i)=dsmoo*(sxb+sc)+smoo*(v1xb+v1cb)
             vtmp(i,2)=smoo*(v2xb+v2cb)
          ENDIF
          IF ((rhoa.GT.small).AND.(rhob.GT.small)) THEN
             vtmp(i,3)=smoo*v2cab
          ENDIF

          id = id + 1
          sgcx = sgcx + (smoo*(sxa+sxb))
          sgcc = sgcc + (smoo*sc)
       ENDIF
    ENDDO
    ! ..counting floating point operations
    fadd   = id * 100._real_8
    fmult  = id * 227._real_8
    fdiv   = id *  42._real_8
    fspez  = id *  12._real_8
    flops=fadd+fmult+fdiv+fspez
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gcsrevpbe
  ! ==================================================================
#endif

END MODULE gcxctbl_utils
