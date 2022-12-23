MODULE nfunc_utils

  IMPLICIT NONE

  PRIVATE

  !public :: lypuu
  !public :: xalpu
  !public :: mikeu
  !public :: mikesp
  !public :: slypsp
  !public :: xalpsp

  !!contains

END MODULE nfunc_utils

#ifndef __VECTOR
! ==================================================================
SUBROUTINE lypuu(nr1,nr2,nr3,km1,km2,km3,rho,uion,eexcu,flop)
  ! ==--------------------------------------------------------------==
  ! adds unpolarized exc potential to uion and sums exc energy
  ! written by S. Goedecker, Stuttgart, march 1996
  ! modified by J. Hutter, Stuttgart, may/august 1996
  ! non-gradient part of Lee,Yang and Parr functional (non-spin polarised)
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE timer,                           ONLY: tihalt,&
       tiset
  IMPLICIT NONE
  INTEGER                                    :: nr1, nr2, nr3, km1, km2, km3
  REAL(real_8)                               :: rho(km1,km2,km3), &
                                                uion(2,km1,km2,km3), eexcu, &
                                                flop

  CHARACTER(*), PARAMETER                    :: procedureN = 'lypuu'
  REAL(real_8), PARAMETER :: a = 0.04918_real_8, b = 0.132_real_8, &
      c = 0.2533_real_8, cf = 2.87123400018819108_real_8 , d = 0.349_real_8 , &
      eps2 = 1.e-20_real_8, f1 = -1.10783814957303361_real_8 , &
      f43 = 4._real_8/3._real_8 , small = 1.e-24_real_8 , &
      third = 1._real_8/3._real_8

  INTEGER                                    :: i1, i2, i3, ic, ii2, ii3, &
                                                isub, ncount
  REAL(real_8) :: bcf, D1, D2, D3, D4, ECRS1, ECRS2, ECRS3, ECRS4, ELYP1, &
      ELYP2, ELYP3, ELYP4, EPSXCU1, EPSXCU2, EPSXCU3, EPSXCU4, fa1, fa2, &
      fadd, fdiv, fmult, fspez, OX1, OX2, OX3, OX4, rhou1, rhou2, rhou3, &
      rhou4, RSX1, RSX2, RSX3, RSX4, VXC1, VXC2, VXC3, VXC4, x1, x2, x3, x4

  CALL tiset(procedureN,isub)
  ! ==--------------------------------------------------------------==
  fa1=f1*2._real_8/3._real_8
  fa2=fa1*f43
  bcf=b*cf
  eexcu=0._real_8
  x1=1._real_8
  x2=1._real_8
  x3=1._real_8
  x4=1._real_8
  ic=0
  ii3=MOD(nr3,2)
  ii2=MOD(nr2,2)
  IF (ii3.NE.0.OR.ii2.NE.0) CALL stopgm('LYPUU','ODD DIMENSION',& 
       __LINE__,__FILE__)
  !$omp parallel do private(I1,I2,I3,RHOU1,RHOU2,RHOU3,RHOU4,D1,D2,D3,D4) &
  !$omp  firstprivate(X1,X2,X3,X4) &
  !$omp  private(NCOUNT,RSX1,RSX2,RSX3,RSX4,ECRS1,ECRS2,ECRS3,ECRS4) &
  !$omp  private(OX1,OX2,OX3,OX4,ELYP1,ELYP2,ELYP3,ELYP4) &
  !$omp  private(EPSXCU1,EPSXCU2,EPSXCU3,EPSXCU4) &
  !$omp  private(VXC1,VXC2,VXC3,VXC4) &
  !$omp  reduction(+:EEXCU,IC)
  ! ..loops are four times unrolled
  DO i3=1,nr3,2
     DO i2=1,nr2,2
        DO i1=1,nr1
           rhou1=MAX(rho(i1,i2+0,i3+0),small)
           rhou2=MAX(rho(i1,i2+1,i3+0),small)
           rhou3=MAX(rho(i1,i2+0,i3+1),small)
           rhou4=MAX(rho(i1,i2+1,i3+1),small)
           ! ..iterative calculation of rho**(1/3) (using Newton method)
           DO ncount=1,100
              d1=x1-rhou1/(x1*x1)
              d2=x2-rhou2/(x2*x2)
              d3=x3-rhou3/(x3*x3)
              d4=x4-rhou4/(x4*x4)
              ic=ic+1
              x1=x1-0.333333333333333_real_8*d1
              x2=x2-0.333333333333333_real_8*d2
              x3=x3-0.333333333333333_real_8*d3
              x4=x4-0.333333333333333_real_8*d4
              IF ((d1*d1+d2*d2+d3*d3+d4*d4).LT.eps2) GOTO 10
           ENDDO
           x1=rhou1**(1._real_8/3._real_8)
           x2=rhou2**(1._real_8/3._real_8)
           x3=rhou3**(1._real_8/3._real_8)
           x4=rhou4**(1._real_8/3._real_8)
10         CONTINUE
           rsx1=1._real_8/x1
           rsx2=1._real_8/x2
           rsx3=1._real_8/x3
           rsx4=1._real_8/x4
           ! 
           ecrs1=bcf*EXP(-c*rsx1)
           ecrs2=bcf*EXP(-c*rsx2)
           ecrs3=bcf*EXP(-c*rsx3)
           ecrs4=bcf*EXP(-c*rsx4)
           ox1=1._real_8/(1._real_8+d*rsx1)
           ox2=1._real_8/(1._real_8+d*rsx2)
           ox3=1._real_8/(1._real_8+d*rsx3)
           ox4=1._real_8/(1._real_8+d*rsx4)
           elyp1=-a*ox1*(1._real_8+ecrs1)
           elyp2=-a*ox2*(1._real_8+ecrs2)
           elyp3=-a*ox3*(1._real_8+ecrs3)
           elyp4=-a*ox4*(1._real_8+ecrs4)
           epsxcu1=fa1*x1+elyp1
           epsxcu2=fa1*x2+elyp2
           epsxcu3=fa1*x3+elyp3
           epsxcu4=fa1*x4+elyp4
           vxc1=fa2*x1+elyp1-0.33333333333333_real_8*rsx1*a*ox1*&
                (d*ox1+ecrs1*(d*ox1+c))
           vxc2=fa2*x2+elyp2-0.33333333333333_real_8*rsx2*a*ox2*&
                (d*ox2+ecrs2*(d*ox2+c))
           vxc3=fa2*x3+elyp3-0.33333333333333_real_8*rsx3*a*ox3*&
                (d*ox3+ecrs3*(d*ox3+c))
           vxc4=fa2*x4+elyp4-0.33333333333333_real_8*rsx4*a*ox4*&
                (d*ox4+ecrs4*(d*ox4+c))
           eexcu=eexcu+epsxcu1*rhou1
           uion(1,i1,i2+0,i3+0)=uion(1,i1,i2+0,i3+0)+vxc1
           eexcu=eexcu+epsxcu2*rhou2
           uion(1,i1,i2+1,i3+0)=uion(1,i1,i2+1,i3+0)+vxc2
           eexcu=eexcu+epsxcu3*rhou3
           uion(1,i1,i2+0,i3+1)=uion(1,i1,i2+0,i3+1)+vxc3
           eexcu=eexcu+epsxcu4*rhou4
           uion(1,i1,i2+1,i3+1)=uion(1,i1,i2+1,i3+1)+vxc4
        ENDDO
     ENDDO
  ENDDO
  ! ..counting floating point operations
  fadd   = ic * 11._real_8 + (nr1*nr2*nr3) *  9._real_8
  fmult  = ic * 12._real_8 + (nr1*nr2*nr3) * 15._real_8
  fdiv   = ic *  4._real_8 + (nr1*nr2*nr3) *  2._real_8
  fspez  =              (nr1*nr2*nr3) *  1._real_8
  flop=fadd+fmult+fdiv+fspez
  ! ==--------------------------------------------------------------==
  CALL tihalt(procedureN,isub)
  ! ==--------------------------------------------------------------==
END SUBROUTINE lypuu
! ==================================================================
SUBROUTINE xalpu(salpha,nr1,nr2,nr3,km1,km2,km3,rho,uion,&
     eexcu,flop)
  ! ==--------------------------------------------------------------==
  ! ..X-Alpha potential               
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE timer,                           ONLY: tihalt,&
       tiset
  IMPLICIT NONE
  REAL(real_8)                               :: salpha
  INTEGER                                    :: nr1, nr2, nr3, km1, km2, km3
  REAL(real_8)                               :: rho(km1,km2,km3), &
                                                uion(2,km1,km2,km3), eexcu, &
                                                flop

  CHARACTER(*), PARAMETER                    :: procedureN = 'xalpu'
  REAL(real_8), PARAMETER :: eps2 = 1.e-20_real_8, &
      f1 = -1.10783814957303361_real_8 , small = 1.e-24_real_8 

  INTEGER                                    :: i1, i2, i3, ic, ii2, ii3, &
                                                isub, ncount
  REAL(real_8) :: D1, D2, D3, D4, EPSXCU1, EPSXCU2, EPSXCU3, EPSXCU4, fa1, &
      fa2, fadd, fdiv, fmult, fspez, rhou1, rhou2, rhou3, rhou4, VXC1, VXC2, &
      VXC3, VXC4, x1, x2, x3, x4

! ==--------------------------------------------------------------==

  CALL tiset(procedureN,isub)
  ! ==--------------------------------------------------------------==
  fa1=f1*salpha
  fa2=fa1*4._real_8/3.0_real_8
  eexcu=0._real_8
  x1=1._real_8
  x2=1._real_8
  x3=1._real_8
  x4=1._real_8
  ic=0
  ii3=MOD(nr3,2)
  ii2=MOD(nr2,2)
  IF (ii3.NE.0.OR.ii2.NE.0) CALL stopgm('XALPU','ODD DIMENSION',& 
       __LINE__,__FILE__)
  !$omp parallel do private(I1,I2,I3,RHOU1,RHOU2,RHOU3,RHOU4,D1,D2,D3,D4) &
  !$omp  firstprivate(X1,X2,X3,X4) &
  !$omp  private(NCOUNT,EPSXCU1,EPSXCU2,EPSXCU3,EPSXCU4) &
  !$omp  private(VXC1,VXC2,VXC3,VXC4) &
  !$omp  reduction(+:EEXCU,IC)
  ! ..loops are four times unrolled
  DO i3=1,nr3,2
     DO i2=1,nr2,2
        DO i1=1,nr1
           rhou1=MAX(rho(i1,i2+0,i3+0),small)
           rhou2=MAX(rho(i1,i2+1,i3+0),small)
           rhou3=MAX(rho(i1,i2+0,i3+1),small)
           rhou4=MAX(rho(i1,i2+1,i3+1),small)
           ! ..iterative calculation of rho**(1/3) (using Newton method)
           DO ncount=1,100
              d1=x1-rhou1/(x1*x1)
              d2=x2-rhou2/(x2*x2)
              d3=x3-rhou3/(x3*x3)
              d4=x4-rhou4/(x4*x4)
              ic=ic+1
              x1=x1-0.333333333333333_real_8*d1
              x2=x2-0.333333333333333_real_8*d2
              x3=x3-0.333333333333333_real_8*d3
              x4=x4-0.333333333333333_real_8*d4
              IF ((d1*d1+d2*d2+d3*d3+d4*d4).LT.eps2) GOTO 10
           ENDDO
           x1=rhou1**(1._real_8/3._real_8)
           x2=rhou2**(1._real_8/3._real_8)
           x3=rhou3**(1._real_8/3._real_8)
           x4=rhou4**(1._real_8/3._real_8)
10         CONTINUE
           epsxcu1=fa1*x1
           epsxcu2=fa1*x2
           epsxcu3=fa1*x3
           epsxcu4=fa1*x4
           vxc1=fa2*x1
           vxc2=fa2*x2
           vxc3=fa2*x3
           vxc4=fa2*x4
           eexcu=eexcu+epsxcu1*rhou1
           uion(1,i1,i2+0,i3+0)=uion(1,i1,i2+0,i3+0)+vxc1
           eexcu=eexcu+epsxcu2*rhou2
           uion(1,i1,i2+1,i3+0)=uion(1,i1,i2+1,i3+0)+vxc2
           eexcu=eexcu+epsxcu3*rhou3
           uion(1,i1,i2+0,i3+1)=uion(1,i1,i2+0,i3+1)+vxc3
           eexcu=eexcu+epsxcu4*rhou4
           uion(1,i1,i2+1,i3+1)=uion(1,i1,i2+1,i3+1)+vxc4
        ENDDO
     ENDDO
  ENDDO
  ! ..counting floating point operations
  fadd   = ic * 11._real_8 + (nr1*nr2*nr3) *  2._real_8
  fmult  = ic * 12._real_8 + (nr1*nr2*nr3) *  3._real_8
  fdiv   = ic *  4._real_8 + (nr1*nr2*nr3) *  0._real_8
  fspez  =              (nr1*nr2*nr3) *  0._real_8
  flop=fadd+fmult+fdiv+fspez
  ! ==--------------------------------------------------------------==
  CALL tihalt(procedureN,isub)
  ! ==--------------------------------------------------------------==
END SUBROUTINE xalpu
! ==================================================================
SUBROUTINE mikeu(nr1,nr2,nr3,km1,km2,km3,rho,uion,eexcu,flop)
  ! ==--------------------------------------------------------------==
  ! ..Pade approximation to LDA functional
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE timer,                           ONLY: tihalt,&
       tiset
  IMPLICIT NONE
  INTEGER                                    :: nr1, nr2, nr3, km1, km2, km3
  REAL(real_8)                               :: rho(km1,km2,km3), &
                                                uion(2,km1,km2,km3), eexcu, &
                                                flop

  CHARACTER(*), PARAMETER                    :: procedureN = 'mikeu'
  REAL(real_8), PARAMETER :: a0u = 0.4581652932831429_real_8, &
      a1u = 2.217058676663745_real_8, a2u = 0.7405551735357053_real_8, &
      a3u = 0.01968227878617998_real_8 , b1u = 1.0_real_8, &
      b2u = 4.504130959426697_real_8, b3u = 1.110667363742916_real_8, &
      b4u = 0.02359291751427506_real_8 , c1u = 4._real_8*a0u*b1u/3._real_8, &
      c2u = 5._real_8*a0u*b2u/3._real_8+a1u*b1u, c3u = 2._real_8*a0u*&
      b3u+4._real_8*a1u*b2u/3._real_8+2._real_8*a2u*b1u/3._real_8, c4u = &
      7._real_8*a0u*b4u/3._real_8+5._real_8*a1u*b3u/3._real_8+a2u*b2u+a3u*b1u/&
      3._real_8
  REAL(real_8), PARAMETER :: c5u = 2._real_8*a1u*b4u+4._real_8*a2u*b3u/&
      3._real_8+2._real_8*a3u*b2u/3._real_8, &
      c6u = 5._real_8*a2u*b4u/3._real_8+a3u*b3u, &
      c7u = 4._real_8*a3u*b4u/3._real_8 , eps2 = 1.e-20_real_8, &
      rsfac = 0.6203504908994000_real_8 , small = 1.e-24_real_8 

  INTEGER                                    :: i1, i2, i3, ic, ii2, ii3, &
                                                isub, ncount
  REAL(real_8) :: BOTU1, BOTU2, BOTU3, BOTU4, d1, d2, d3, d4, DTOPU1, DTOPU2, &
      DTOPU3, DTOPU4, EPSXCU1, EPSXCU2, EPSXCU3, EPSXCU4, fadd, fdiv, fmult, &
      fspez, rhou1, rhou2, rhou3, rhou4, RSU1, RSU2, RSU3, RSU4, TOPU1, &
      TOPU2, TOPU3, TOPU4, VXC1, VXC2, VXC3, VXC4, x1, x2, x3, x4

  CALL tiset(procedureN,isub)
  ! ==--------------------------------------------------------------==
  eexcu=0._real_8
  x1=1._real_8
  x2=1._real_8
  x3=1._real_8
  x4=1._real_8
  ic=0
  ii3=MOD(nr3,2)
  ii2=MOD(nr2,2)
  IF (ii3.NE.0.OR.ii2.NE.0) CALL stopgm('MIKEU','ODD DIMENSION',& 
       __LINE__,__FILE__)
  ! ..loops are four times unrolled
  !$omp parallel do private(I1,I2,I3,RHOU1,RHOU2,RHOU3,RHOU4,D1,D2,D3,D4) &
  !$omp  firstprivate(X1,X2,X3,X4) &
  !$omp  private(NCOUNT,RSU1,RSU2,RSU3,RSU4,TOPU1,TOPU2,TOPU3,TOPU4) &
  !$omp  private(DTOPU1,DTOPU2,DTOPU3,DTOPU4,BOTU1,BOTU2,BOTU3,BOTU4) &
  !$omp  private(EPSXCU1,EPSXCU2,EPSXCU3,EPSXCU4,VXC1,VXC2,VXC3,VXC4) &
  !$omp  reduction(+:EEXCU,IC)
  DO i3=1,nr3,2
     DO i2=1,nr2,2
        DO i1=1,nr1
           rhou1=MAX(rho(i1,i2+0,i3+0),small)
           rhou2=MAX(rho(i1,i2+1,i3+0),small)
           rhou3=MAX(rho(i1,i2+0,i3+1),small)
           rhou4=MAX(rho(i1,i2+1,i3+1),small)
           ! ..iterative calculation of rho**(1/3) (using Newton method)
           DO ncount=1,100
              d1=x1-rhou1/(x1*x1)
              d2=x2-rhou2/(x2*x2)
              d3=x3-rhou3/(x3*x3)
              d4=x4-rhou4/(x4*x4)
              ic=ic+1
              x1=x1-0.333333333333333_real_8*d1
              x2=x2-0.333333333333333_real_8*d2
              x3=x3-0.333333333333333_real_8*d3
              x4=x4-0.333333333333333_real_8*d4
              IF ((d1*d1+d2*d2+d3*d3+d4*d4).LT.eps2) GOTO 10
           ENDDO
           x1=rhou1**(1._real_8/3._real_8)
           x2=rhou2**(1._real_8/3._real_8)
           x3=rhou3**(1._real_8/3._real_8)
           x4=rhou4**(1._real_8/3._real_8)
10         CONTINUE
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
           epsxcu1=topu1*(-1._real_8/botu1)
           epsxcu2=topu2*(-1._real_8/botu2)
           epsxcu3=topu3*(-1._real_8/botu3)
           epsxcu4=topu4*(-1._real_8/botu4)
           eexcu=eexcu+epsxcu1*rhou1
           vxc1=dtopu1*(-1._real_8/botu1)*(-1._real_8/botu1)
           uion(1,i1,i2+0,i3+0)=uion(1,i1,i2+0,i3+0)+vxc1
           eexcu=eexcu+epsxcu2*rhou2
           vxc2=dtopu2*(-1._real_8/botu2)*(-1._real_8/botu2)
           uion(1,i1,i2+1,i3+0)=uion(1,i1,i2+1,i3+0)+vxc2
           eexcu=eexcu+epsxcu3*rhou3
           vxc3=dtopu3*(-1._real_8/botu3)*(-1._real_8/botu3)
           uion(1,i1,i2+0,i3+1)=uion(1,i1,i2+0,i3+1)+vxc3
           eexcu=eexcu+epsxcu4*rhou4
           vxc4=dtopu4*(-1._real_8/botu4)*(-1._real_8/botu4)
           uion(1,i1,i2+1,i3+1)=uion(1,i1,i2+1,i3+1)+vxc4
        ENDDO
     ENDDO
  ENDDO
  ! ..counting floating point operations
  fadd   = ic * 11._real_8 + (nr1*nr2*nr3) * 14._real_8
  fmult  = ic * 12._real_8 + (nr1*nr2*nr3) * 18._real_8
  fdiv   = ic *  4._real_8 + (nr1*nr2*nr3) *  2._real_8
  fspez  =              (nr1*nr2*nr3) *  0._real_8
  flop=fadd+fmult+fdiv+fspez
  ! ==--------------------------------------------------------------==
  CALL tihalt(procedureN,isub)
  ! ==--------------------------------------------------------------==
END SUBROUTINE mikeu
! ==================================================================
SUBROUTINE mikesp(nr1,nr2,nr3,km1,km2,km3,rho,uion1,uion2,&
     eexcu,flop)
  ! ==--------------------------------------------------------------==
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE timer,                           ONLY: tihalt,&
       tiset
  IMPLICIT NONE
  INTEGER                                    :: nr1, nr2, nr3, km1, km2, km3
  REAL(real_8)                               :: rho(km1,km2,km3,2), &
                                                uion1(2,km1,km2,km3), &
                                                uion2(2,km1,km2,km3), eexcu, &
                                                flop

  CHARACTER(*), PARAMETER                    :: procedureN = 'mikesp'
  REAL(real_8), PARAMETER :: a0u = 0.4581652932831429_real_8, &
      a1u = 2.217058676663745_real_8, a2u = 0.7405551735357053_real_8, &
      a3u = 0.01968227878617998_real_8 , b1u = 1.0_real_8, &
      b2u = 4.504130959426697_real_8, b3u = 1.110667363742916_real_8, &
      b4u = 0.02359291751427506_real_8 , da0 = 0.119086804055547_real_8, &
      da1 = 0.6157402568883345_real_8, da2 = 0.1574201515892867_real_8, &
      da3 = 0.003532336663397157_real_8 , db1 = 0.0_real_8, &
      db2 = 0.2673612973836267_real_8, db3 = 0.2052004607777787_real_8, &
      db4 = 0.004200005045691381_real_8 , eps2 = 1.e-20_real_8, &
      rsfac = 0.6203504908994000_real_8 
  REAL(real_8), PARAMETER :: small = 1.e-24_real_8 , &
      t1 = 0.854960467080683_real_8, t2 = 0.07916300621117439_real_8, &
      t3 = 0.02580127609845684_real_8, t4 = 0.0121839359353824_real_8, &
      t5 = 0.006919272259599879_real_8, t6 = 0.004391524649611372_real_8, &
      t7 = 0.003002751897170168_real_8, t8 = 0.00216587382212552_real_8, &
      t9 = 0.01141189204579585_real_8 

  INTEGER                                    :: i1, i2, i3, ic, ii2, ii3, &
                                                isub, ncount
  REAL(real_8) :: A0P1, A0P2, A0P3, A0P4, A1P1, A1P2, A1P3, A1P4, A2P1, A2P2, &
      A2P3, A2P4, A3P1, A3P2, A3P3, A3P4, b1p1, b1p2, b1p3, b1p4, B2P1, B2P2, &
      B2P3, B2P4, B3P1, B3P2, B3P3, B3P4, B4P1, B4P2, B4P3, B4P4, BOT1, BOT2, &
      BOT3, BOT4, d1, d2, d3, d4, DBOT1, DBOT2, DBOT3, DBOT4, DFXS1, DFXS2, &
      DFXS3, DFXS4, DTOP1, DTOP2, DTOP3, DTOP4, DVXC1, DVXC2, DVXC3, DVXC4, &
      DXA1, DXA2, DXA3, DXA4, DXB1, DXB2, DXB3, DXB4, EPSXCU1, EPSXCU2, &
      EPSXCU3, EPSXCU4, ET1, ET2, ET3, ET4, ETA1, ETA2, ETA3, ETA4, fadd, &
      fdiv, fmult, fspez, FXS1, FXS2, FXS3, FXS4, OBOT1, OBOT2, OBOT3, OBOT4, &
      RHOA1, RHOA2, RHOA3, RHOA4, RHOB1
  REAL(real_8) :: RHOB2, RHOB3, RHOB4, rhou1, rhou2, rhou3, rhou4, RSU1, &
      RSU2, RSU3, RSU4, TOP1, TOP2, TOP3, TOP4, VXC1, VXC2, VXC3, VXC4, x1, &
      x2, x3, x4, XBOT1, XBOT2, XBOT3, XBOT4, XTOP1, XTOP2, XTOP3, XTOP4

! ==--------------------------------------------------------------==

  CALL tiset(procedureN,isub)
  ! ==--------------------------------------------------------------==
  eexcu=0._real_8
  x1=1._real_8
  x2=1._real_8
  x3=1._real_8
  x4=1._real_8
  ic=0
  ii3=MOD(nr3,2)
  ii2=MOD(nr2,2)
  IF (ii3.NE.0.OR.ii2.NE.0) CALL stopgm('MIKESP','ODD DIMENSION',& 
       __LINE__,__FILE__)
  !$omp parallel do private(I1,I2,I3,RHOU1,RHOU2,RHOU3,RHOU4,D1,D2,D3,D4) &
  !$omp  firstprivate(X1,X2,X3,X4) &
  !$omp  private(NCOUNT,RHOA1,RHOA2,RHOA3,RHOA4,RHOB1,RHOB2,RHOB3,RHOB4) &
  !$omp  private(ETA1,ETA2,ETA3,ETA4,ET1,FXS1,DFXS1,ET2,FXS2) &
  !$omp  private(DFXS2,ET3,FXS3,DFXS3,ET4,FXS4,DFXS4,A0P1,A1P1) &
  !$omp  private(A2P1,A3P1,B2P1,B3P1,B4P1,A0P2,A1P2,A2P2,A3P2) &
  !$omp  private(B2P2,B3P2,B4P2,A0P3,A1P3,A2P3,A3P3,B2P3,B3P3) &
  !$omp  private(B4P3,A0P4,A1P4,A2P4,A3P4,B2P4,B3P4,B4P4,RSU1) &
  !$omp  private(RSU2,RSU3,RSU4,TOP1,TOP2,TOP3,TOP4,DTOP1,DTOP2) &
  !$omp  private(DTOP3,DTOP4,XTOP1,XTOP2,XTOP3,XTOP4,BOT1,BOT2) &
  !$omp  private(BOT3,BOT4,DBOT1,DBOT2,DBOT3,DBOT4,XBOT1,XBOT2) &
  !$omp  private(XBOT3,XBOT4,OBOT1,OBOT2,OBOT3,OBOT4,EPSXCU1) &
  !$omp  private(EPSXCU2,EPSXCU3,EPSXCU4,DVXC1,DXA1,DXB1,DVXC2) &
  !$omp  private(DXA2,DXB2,DVXC3,DXA3,DXB3,DVXC4,DXA4,DXB4) &
  !$omp  private(VXC1,VXC2,VXC3,VXC4) &
  !$omp  private(B1P1,B1P2,B1P3,B1P4)  &
  !$omp  reduction(+:EEXCU,IC)
  DO i3=1,nr3,2
     DO i2=1,nr2,2
        DO i1=1,nr1
           rhoa1=MAX(rho(i1,i2+0,i3+0,1)-rho(i1,i2+0,i3+0,2),small)
           rhoa2=MAX(rho(i1,i2+1,i3+0,1)-rho(i1,i2+1,i3+0,2),small)
           rhoa3=MAX(rho(i1,i2+0,i3+1,1)-rho(i1,i2+0,i3+1,2),small)
           rhoa4=MAX(rho(i1,i2+1,i3+1,1)-rho(i1,i2+1,i3+1,2),small)
           rhob1=MAX(rho(i1,i2+0,i3+0,2),small)
           rhob2=MAX(rho(i1,i2+1,i3+0,2),small)
           rhob3=MAX(rho(i1,i2+0,i3+1,2),small)
           rhob4=MAX(rho(i1,i2+1,i3+1,2),small)
           rhou1=rhoa1+rhob1
           rhou2=rhoa2+rhob2
           rhou3=rhoa3+rhob3
           rhou4=rhoa4+rhob4
           eta1=(rhoa1-rhob1)/rhou1
           eta2=(rhoa2-rhob2)/rhou2
           eta3=(rhoa3-rhob3)/rhou3
           eta4=(rhoa4-rhob4)/rhou4
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
           et3=eta3*eta3
           fxs3=et3*(t1+et3*(t2+et3*(t3+et3*(t4+et3*(t5+et3*&
                (t6+et3*(t7+et3*(t8+et3*t9))))))))
           dfxs3=(2._real_8*eta3)*(t1+et3*(2._real_8*t2+et3*(3._real_8*t3+et3*&
                (4._real_8*t4+et3*(5._real_8*t5+et3*(6._real_8*t6+et3*(7._real_8*t7+et3*&
                (8._real_8*t8+et3*9._real_8*t9))))))))
           et4=eta4*eta4
           fxs4=et4*(t1+et4*(t2+et4*(t3+et4*(t4+et4*(t5+et4*&
                (t6+et4*(t7+et4*(t8+et4*t9))))))))
           dfxs4=(2._real_8*eta4)*(t1+et4*(2._real_8*t2+et4*(3._real_8*t3+et4*&
                (4._real_8*t4+et4*(5._real_8*t5+et4*(6._real_8*t6+et4*(7._real_8*t7+et4*&
                (8._real_8*t8+et4*9._real_8*t9))))))))
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
           a0p3=a0u+fxs3*da0
           a1p3=a1u+fxs3*da1
           a2p3=a2u+fxs3*da2
           a3p3=a3u+fxs3*da3
           b1p3=b1u+fxs3*db1
           b2p3=b2u+fxs3*db2
           b3p3=b3u+fxs3*db3
           b4p3=b4u+fxs3*db4
           ! 
           a0p4=a0u+fxs4*da0
           a1p4=a1u+fxs4*da1
           a2p4=a2u+fxs4*da2
           a3p4=a3u+fxs4*da3
           b1p4=b1u+fxs4*db1
           b2p4=b2u+fxs4*db2
           b3p4=b3u+fxs4*db3
           b4p4=b4u+fxs4*db4
           ! 
           ! ..iterative calculation of rho**(1/3) (using Newton method)
           DO ncount=1,100
              d1=x1-rhou1/(x1*x1)
              d2=x2-rhou2/(x2*x2)
              d3=x3-rhou3/(x3*x3)
              d4=x4-rhou4/(x4*x4)
              ic=ic+1
              x1=x1-0.333333333333333_real_8*d1
              x2=x2-0.333333333333333_real_8*d2
              x3=x3-0.333333333333333_real_8*d3
              x4=x4-0.333333333333333_real_8*d4
              IF ((d1*d1+d2*d2+d3*d3+d4*d4).LT.eps2) GOTO 10
           ENDDO
           x1=rhou1**(1._real_8/3._real_8)
           x2=rhou2**(1._real_8/3._real_8)
           x3=rhou3**(1._real_8/3._real_8)
           x4=rhou4**(1._real_8/3._real_8)
10         CONTINUE
           rsu1=rsfac/x1
           rsu2=rsfac/x2
           rsu3=rsfac/x3
           rsu4=rsfac/x4
           top1=a0p1+rsu1*(a1p1+rsu1*(a2p1+rsu1*a3p1))
           top2=a0p2+rsu2*(a1p2+rsu2*(a2p2+rsu2*a3p2))
           top3=a0p3+rsu3*(a1p3+rsu3*(a2p3+rsu3*a3p3))
           top4=a0p4+rsu4*(a1p4+rsu4*(a2p4+rsu4*a3p4))
           dtop1= a1p1+rsu1*(2._real_8*a2p1+rsu1*3._real_8*a3p1)
           dtop2= a1p2+rsu2*(2._real_8*a2p2+rsu2*3._real_8*a3p2)
           dtop3= a1p3+rsu3*(2._real_8*a2p3+rsu3*3._real_8*a3p3)
           dtop4= a1p4+rsu4*(2._real_8*a2p4+rsu4*3._real_8*a3p4)
           xtop1=da0+rsu1*(da1+rsu1*(da2+rsu1*da3))
           xtop2=da0+rsu2*(da1+rsu2*(da2+rsu2*da3))
           xtop3=da0+rsu3*(da1+rsu3*(da2+rsu3*da3))
           xtop4=da0+rsu4*(da1+rsu4*(da2+rsu4*da3))
           bot1=rsu1*(b1p1+rsu1*(b2p1+rsu1*(b3p1+rsu1*b4p1)))
           bot2=rsu2*(b1p2+rsu2*(b2p2+rsu2*(b3p2+rsu2*b4p2)))
           bot3=rsu3*(b1p3+rsu3*(b2p3+rsu3*(b3p3+rsu3*b4p3)))
           bot4=rsu4*(b1p4+rsu4*(b2p4+rsu4*(b3p4+rsu4*b4p4)))
           dbot1=b1p1+rsu1*(2._real_8*b2p1+rsu1*(3._real_8*b3p1+rsu1*4._real_8*b4p1))
           dbot2=b1p2+rsu2*(2._real_8*b2p2+rsu2*(3._real_8*b3p2+rsu2*4._real_8*b4p2))
           dbot3=b1p3+rsu3*(2._real_8*b2p3+rsu3*(3._real_8*b3p3+rsu3*4._real_8*b4p3))
           dbot4=b1p4+rsu4*(2._real_8*b2p4+rsu4*(3._real_8*b3p4+rsu4*4._real_8*b4p4))
           xbot1=rsu1*(db1+rsu1*(db2+rsu1*(db3+rsu1*db4)))
           xbot2=rsu2*(db1+rsu2*(db2+rsu2*(db3+rsu2*db4)))
           xbot3=rsu3*(db1+rsu3*(db2+rsu3*(db3+rsu3*db4)))
           xbot4=rsu4*(db1+rsu4*(db2+rsu4*(db3+rsu4*db4)))
           obot1=1._real_8/bot1
           obot2=1._real_8/bot2
           obot3=1._real_8/bot3
           obot4=1._real_8/bot4
           epsxcu1=-top1*obot1
           epsxcu2=-top2*obot2
           epsxcu3=-top3*obot3
           epsxcu4=-top4*obot4
           eexcu=eexcu+epsxcu1*rhou1
           vxc1=epsxcu1+rsu1*(dtop1*obot1-dbot1*(top1*obot1*obot1))*&
                0.333333333333333_real_8
           dvxc1=-(xtop1*obot1-xbot1*(top1*obot1*obot1))
           dxa1=vxc1+dvxc1*dfxs1*(1._real_8-eta1)
           dxb1=vxc1-dvxc1*dfxs1*(1._real_8+eta1)
           uion1(1,i1,i2+0,i3+0)=uion1(1,i1,i2+0,i3+0)+dxa1
           uion2(1,i1,i2+0,i3+0)=uion2(1,i1,i2+0,i3+0)+dxb1
           eexcu=eexcu+epsxcu2*rhou2
           vxc2=epsxcu2+rsu2*(dtop2*obot2-dbot2*(top2*obot2*obot2))*&
                0.333333333333333_real_8
           dvxc2=-(xtop2*obot2-xbot2*(top2*obot2*obot2))
           dxa2=vxc2+dvxc2*dfxs2*(1._real_8-eta2)
           dxb2=vxc2-dvxc2*dfxs2*(1._real_8+eta2)
           uion1(1,i1,i2+1,i3+0)=uion1(1,i1,i2+1,i3+0)+dxa2
           uion2(1,i1,i2+1,i3+0)=uion2(1,i1,i2+1,i3+0)+dxb2
           eexcu=eexcu+epsxcu3*rhou3
           vxc3=epsxcu3+rsu3*(dtop3*obot3-dbot3*(top3*obot3*obot3))*&
                0.333333333333333_real_8
           dvxc3=-(xtop3*obot3-xbot3*(top3*obot3*obot3))
           dxa3=vxc3+dvxc3*dfxs3*(1._real_8-eta3)
           dxb3=vxc3-dvxc3*dfxs3*(1._real_8+eta3)
           uion1(1,i1,i2+0,i3+1)=uion1(1,i1,i2+0,i3+1)+dxa3
           uion2(1,i1,i2+0,i3+1)=uion2(1,i1,i2+0,i3+1)+dxb3
           eexcu=eexcu+epsxcu4*rhou4
           vxc4=epsxcu4+rsu4*(dtop4*obot4-dbot4*(top4*obot4*obot4))*&
                0.333333333333333_real_8
           dvxc4=-(xtop4*obot4-xbot4*(top4*obot4*obot4))
           dxa4=vxc4+dvxc4*dfxs4*(1._real_8-eta4)
           dxb4=vxc4-dvxc4*dfxs4*(1._real_8+eta4)
           uion1(1,i1,i2+1,i3+1)=uion1(1,i1,i2+1,i3+1)+dxa4
           uion2(1,i1,i2+1,i3+1)=uion2(1,i1,i2+1,i3+1)+dxb4
        ENDDO
     ENDDO
  ENDDO
  ! ..counting floating point operations
  fadd   = ic * 11._real_8 + (nr1*nr2*nr3) * 53._real_8
  fmult  = ic * 12._real_8 + (nr1*nr2*nr3) * 76._real_8
  fdiv   = ic *  4._real_8 + (nr1*nr2*nr3) *  3._real_8
  fspez  =              (nr1*nr2*nr3) *  0._real_8
  flop=fadd+fmult+fdiv+fspez
  ! ==--------------------------------------------------------------==
  CALL tihalt(procedureN,isub)
  ! ==--------------------------------------------------------------==
END SUBROUTINE mikesp
! ==================================================================
SUBROUTINE slypsp(nr1,nr2,nr3,km1,km2,km3,rho,uion1,uion2,salpha,&
     eexcu,flop)
  ! ==--------------------------------------------------------------==
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE timer,                           ONLY: tihalt,&
       tiset
  IMPLICIT NONE
  INTEGER                                    :: nr1, nr2, nr3, km1, km2, km3
  REAL(real_8)                               :: rho(km1,km2,km3,2), &
                                                uion1(2,km1,km2,km3), &
                                                uion2(2,km1,km2,km3), salpha, &
                                                eexcu, flop

  CHARACTER(*), PARAMETER                    :: procedureN = 'slypsp'
  REAL(real_8), PARAMETER :: a = 0.04918_real_8, b = 0.132_real_8, &
      c = 0.2533_real_8, cf = 2.87123400018819108_real_8 , d = 0.349_real_8, &
      eps2 = 1.e-20_real_8, f1 = -1.39578858466194911_real_8 , &
      small = 1.e-24_real_8 

  INTEGER                                    :: i1, i2, i3, ic, ii2, ii3, &
                                                isub, ncount
  REAL(real_8) :: cfx, d1, d2, d3, d4, DD1, DD2, DD3, DD4, DOR1, DOR2, DOR3, &
      DOR4, EC1, EC2, EC3, EC4, ECR1, ECR2, ECR3, ECR4, EX1, EX2, EX3, EX4, &
      f1s, f1u, fadd, fdiv, fmult, fspez, OR1, OR2, OR3, OR4, RHO1, RHO2, &
      RHO3, RHO4, RHOA1, RHOA2, RHOA3, RHOA4, RHOAT1, RHOAT2, RHOAT3, RHOAT4, &
      RHOB1, RHOB2, RHOB3, RHOB4, RHOBT1, RHOBT2, RHOBT3, RHOBT4, RHOU1, &
      RHOU2, RHOU3, RHOU4, rs1, rs2, rs3, rs4, RSA1, RSA2, RSA3, RSA4, rsb1, &
      rsb2, rsb3, rsb4, VCA1, VCA2, VCA3, VCA4, VCB1, VCB2, VCB3, VCB4, VXA1, &
      VXA2, VXA3, VXA4, VXB1, VXB2, VXB3, VXB4, x1, x2, x3, x4, y1, y2, y3, y4

! ==--------------------------------------------------------------==

  CALL tiset(procedureN,isub)
  ! ==--------------------------------------------------------------==
  cfx=cf*2._real_8**(11._real_8/3._real_8)
  f1s=f1*salpha
  f1u=f1s*1.3333333333333333_real_8
  eexcu=0._real_8
  x1=1._real_8
  x2=1._real_8
  x3=1._real_8
  x4=1._real_8
  y1=1._real_8
  y2=1._real_8
  y3=1._real_8
  y4=1._real_8
  ic=0
  ii3=MOD(nr3,2)
  ii2=MOD(nr2,2)
  IF (ii3.NE.0.OR.ii2.NE.0) CALL stopgm('SLYPSP','ODD DIMENSION',& 
       __LINE__,__FILE__)
  !$omp parallel do private(I1,I2,I3,RHOU1,RHOU2,RHOU3,RHOU4,D1,D2,D3,D4) &
  !$omp  firstprivate(X1,X2,X3,X4,Y1,Y2,Y3,Y4) &
  !$omp  private(RHOA1,RHOA2,RHOA3,RHOA4,RHOB1,RHOB2,RHOB3,RHOB4) &
  !$omp  private(NCOUNT,RHO1,RHO2,RHO3,RHO4,RSA1,RSA2,RSA3,RSA4) &
  !$omp  private(EX1,EX2,EX3,EX4,VXA1,VXA2,VXA3,VXA4,VXB1,VXB2,VXB3,VXB4) &
  !$omp  private(DD1,DD2,DD3,DD4,ECR1,ECR2,ECR3,ECR4,OR1,OR2,OR3,OR4) &
  !$omp  private(DOR1,DOR2,DOR3,DOR4,RHOAT1,RHOAT2,RHOAT3,RHOAT4) &
  !$omp  private(RHOBT1,RHOBT2,RHOBT3,RHOBT4,EC1,EC2,EC3,EC4) &
  !$omp  private(VCA1,VCA2,VCA3,VCA4,VCB1,VCB2,VCB3,VCB4) &
  !$omp  private(RSB1,RSB2,RSB3,RSB4,RS1,RS2,RS3,RS4)  &
  !$omp  reduction(+:EEXCU,IC)
  DO i3=1,nr3,2
     DO i2=1,nr2,2
        DO i1=1,nr1
           rhoa1=rho(i1,i2+0,i3+0,1)-rho(i1,i2+0,i3+0,2)
           rhoa2=rho(i1,i2+1,i3+0,1)-rho(i1,i2+1,i3+0,2)
           rhoa3=rho(i1,i2+0,i3+1,1)-rho(i1,i2+0,i3+1,2)
           rhoa4=rho(i1,i2+1,i3+1,1)-rho(i1,i2+1,i3+1,2)
           rhoa1=MAX(rhoa1,small)
           rhoa2=MAX(rhoa2,small)
           rhoa3=MAX(rhoa3,small)
           rhoa4=MAX(rhoa4,small)
           rhob1=MAX(rho(i1,i2+0,i3+0,2),small)
           rhob2=MAX(rho(i1,i2+1,i3+0,2),small)
           rhob3=MAX(rho(i1,i2+0,i3+1,2),small)
           rhob4=MAX(rho(i1,i2+1,i3+1,2),small)
           rho1=rhoa1+rhob1
           rho2=rhoa2+rhob2
           rho3=rhoa3+rhob3
           rho4=rhoa4+rhob4
           ! 
           ! ..iterative calculation of rho**(1/3) (using Newton method)
           DO ncount=1,100
              d1=x1-rhoa1/(x1*x1)
              d2=x2-rhoa2/(x2*x2)
              d3=x3-rhoa3/(x3*x3)
              d4=x4-rhoa4/(x4*x4)
              ic=ic+1
              x1=x1-0.333333333333333_real_8*d1
              x2=x2-0.333333333333333_real_8*d2
              x3=x3-0.333333333333333_real_8*d3
              x4=x4-0.333333333333333_real_8*d4
              IF ((d1*d1+d2*d2+d3*d3+d4*d4).LT.eps2) GOTO 10
           ENDDO
           x1=rhoa1**(1._real_8/3._real_8)
           x2=rhoa2**(1._real_8/3._real_8)
           x3=rhoa3**(1._real_8/3._real_8)
           x4=rhoa4**(1._real_8/3._real_8)
10         CONTINUE
           rsa1=x1
           rsa2=x2
           rsa3=x3
           rsa4=x4
           DO ncount=1,100
              d1=x1-rhob1/(x1*x1)
              d2=x2-rhob2/(x2*x2)
              d3=x3-rhob3/(x3*x3)
              d4=x4-rhob4/(x4*x4)
              ic=ic+1
              x1=x1-0.333333333333333_real_8*d1
              x2=x2-0.333333333333333_real_8*d2
              x3=x3-0.333333333333333_real_8*d3
              x4=x4-0.333333333333333_real_8*d4
              IF ((d1*d1+d2*d2+d3*d3+d4*d4).LT.eps2) GOTO 20
           ENDDO
           x1=rhob1**(1._real_8/3._real_8)
           x2=rhob2**(1._real_8/3._real_8)
           x3=rhob3**(1._real_8/3._real_8)
           x4=rhob4**(1._real_8/3._real_8)
20         CONTINUE
           rsb1=x1
           rsb2=x2
           rsb3=x3
           rsb4=x4
           DO ncount=1,100
              d1=y1-rho1/(y1*y1)
              d2=y2-rho2/(y2*y2)
              d3=y3-rho3/(y3*y3)
              d4=y4-rho4/(y4*y4)
              ic=ic+1
              y1=y1-0.333333333333333_real_8*d1
              y2=y2-0.333333333333333_real_8*d2
              y3=y3-0.333333333333333_real_8*d3
              y4=y4-0.333333333333333_real_8*d4
              IF ((d1*d1+d2*d2+d3*d3+d4*d4).LT.eps2) GOTO 30
           ENDDO
           y1=rho1**(1._real_8/3._real_8)
           y2=rho2**(1._real_8/3._real_8)
           y3=rho3**(1._real_8/3._real_8)
           y4=rho4**(1._real_8/3._real_8)
30         CONTINUE
           rs1=y1
           rs2=y2
           rs3=y3
           rs4=y4
           ! ..Slater exchange
           ex1=f1s*(rhoa1*rsa1+rhob1*rsb1)
           ex2=f1s*(rhoa2*rsa2+rhob2*rsb2)
           ex3=f1s*(rhoa3*rsa3+rhob3*rsb3)
           ex4=f1s*(rhoa4*rsa4+rhob4*rsb4)
           vxa1=f1u*rsa1
           vxa2=f1u*rsa2
           vxa3=f1u*rsa3
           vxa4=f1u*rsa4
           vxb1=f1u*rsb1
           vxb2=f1u*rsb2
           vxb3=f1u*rsb3
           vxb4=f1u*rsb4
           ! ..Lee Yang Parr
           dd1=1._real_8/(1._real_8+d/rs1)
           dd2=1._real_8/(1._real_8+d/rs2)
           dd3=1._real_8/(1._real_8+d/rs3)
           dd4=1._real_8/(1._real_8+d/rs4)
           ecr1=EXP(-c/rs1)
           ecr2=EXP(-c/rs2)
           ecr3=EXP(-c/rs3)
           ecr4=EXP(-c/rs4)
           or1=ecr1*dd1/(rho1*rho1*rho1*rs1*rs1)
           or2=ecr2*dd2/(rho2*rho2*rho2*rs2*rs2)
           or3=ecr3*dd3/(rho3*rho3*rho3*rs3*rs3)
           or4=ecr4*dd4/(rho4*rho4*rho4*rs4*rs4)
           dor1=-or1*0.333333333333333_real_8/(rho1*rs1)*(11._real_8*rs1-c-d*dd1)
           dor2=-or2*0.333333333333333_real_8/(rho2*rs2)*(11._real_8*rs2-c-d*dd2)
           dor3=-or3*0.333333333333333_real_8/(rho3*rs3)*(11._real_8*rs3-c-d*dd3)
           dor4=-or4*0.333333333333333_real_8/(rho4*rs4)*(11._real_8*rs4-c-d*dd4)
           rhoat1=rhoa1*rhoa1*rsa1*rsa1
           rhobt1=rhob1*rhob1*rsb1*rsb1
           rhoat2=rhoa2*rhoa2*rsa2*rsa2
           rhobt2=rhob2*rhob2*rsb2*rsb2
           rhoat3=rhoa3*rhoa3*rsa3*rsa3
           rhobt3=rhob3*rhob3*rsb3*rsb3
           rhoat4=rhoa4*rhoa4*rsa4*rsa4
           rhobt4=rhob4*rhob4*rsb4*rsb4
           ec1=-4._real_8*a*rhoa1*rhob1/rho1*dd1-cfx*a*b*or1*rhoa1*rhob1*&
                (rhoat1+rhobt1)
           vca1=-4._real_8*a*dd1*rhoa1*rhob1/rho1*(d*0.33333333333333_real_8*dd1/&
                (rho1*rs1)+(1._real_8/rhoa1)-(1._real_8/rho1))-cfx*a*b*(dor1*&
                rhoa1*rhob1*(rhoat1+rhobt1)+or1*rhob1*&
                (3.6666666666666667_real_8*rhoat1+rhobt1))
           vcb1=-4._real_8*a*dd1*rhoa1*rhob1/rho1*(d*0.33333333333333_real_8*&
                dd1/(rho1*rs1)+(1._real_8/rhob1)-(1._real_8/rho1))-cfx*a*b*&
                (dor1*rhoa1*rhob1*(rhoat1+rhobt1)+or1*rhoa1*&
                (3.6666666666666667_real_8*rhobt1+rhoat1))
           ec2=-4._real_8*a*rhoa2*rhob2/rho2*dd2-cfx*a*b*or2*rhoa2*rhob2*&
                (rhoat2+rhobt2)
           vca2=-4._real_8*a*dd2*rhoa2*rhob2/rho2*(d*0.33333333333333_real_8*dd2/&
                (rho2*rs2)+(1._real_8/rhoa2)-(1._real_8/rho2))-cfx*a*b*(dor2*&
                rhoa2*rhob2*(rhoat2+rhobt2)+or2*rhob2*&
                (3.6666666666666667_real_8*rhoat2+rhobt2))
           vcb2=-4._real_8*a*dd2*rhoa2*rhob2/rho2*(d*0.33333333333333_real_8*dd2/&
                (rho2*rs2)+(1._real_8/rhob2)-(1._real_8/rho2))-cfx*a*b*(dor2*&
                rhoa2*rhob2*(rhoat2+rhobt2)+or2*rhoa2*&
                (3.6666666666666667_real_8*rhobt2+rhoat2))
           ec3=-4._real_8*a*rhoa3*rhob3/rho3*dd3-cfx*a*b*or3*rhoa3*rhob3*&
                (rhoat3+rhobt3)
           vca3=-4._real_8*a*dd3*rhoa3*rhob3/rho3*(d*0.33333333333333_real_8*dd3/&
                (rho3*rs3)+(1._real_8/rhoa3)-(1._real_8/rho3))-cfx*a*b*(dor3*&
                rhoa3*rhob3*(rhoat3+rhobt3)+or3*rhob3*&
                (3.6666666666666667_real_8*rhoat3+rhobt3))
           vcb3=-4._real_8*a*dd3*rhoa3*rhob3/rho3*(d*0.33333333333333_real_8*dd3/&
                (rho3*rs3)+(1._real_8/rhob3)-(1._real_8/rho3))-cfx*a*b*(dor3*&
                rhoa3*rhob3*(rhoat3+rhobt3)+or3*rhoa3*&
                (3.6666666666666667_real_8*rhobt3+rhoat3))
           ec4=-4._real_8*a*rhoa4*rhob4/rho4*dd4-cfx*a*b*or4*rhoa4*rhob4*&
                (rhoat4+rhobt4)
           vca4=-4._real_8*a*dd4*rhoa4*rhob4/rho4*(d*0.33333333333333_real_8*dd4/&
                (rho4*rs4)+(1._real_8/rhoa4)-(1._real_8/rho4))-cfx*a*b*(dor4*&
                rhoa4*rhob4*(rhoat4+rhobt4)+or4*rhob4*&
                (3.6666666666666667_real_8*rhoat4+rhobt4))
           vcb4=-4._real_8*a*dd4*rhoa4*rhob4/rho4*(d*0.33333333333333_real_8*dd4/&
                (rho4*rs4)+(1._real_8/rhob4)-(1._real_8/rho4))-cfx*a*b*(dor4*&
                rhoa4*rhob4*(rhoat4+rhobt4)+or4*rhoa4*&
                (3.6666666666666667_real_8*rhobt4+rhoat4))
           eexcu=eexcu+(ex1+ec1)
           eexcu=eexcu+(ex2+ec2)
           eexcu=eexcu+(ex3+ec3)
           eexcu=eexcu+(ex4+ec4)
           uion1(1,i1,i2+0,i3+0)=uion1(1,i1,i2+0,i3+0)+vxa1+vca1
           uion1(1,i1,i2+1,i3+0)=uion1(1,i1,i2+1,i3+0)+vxa2+vca2
           uion1(1,i1,i2+0,i3+1)=uion1(1,i1,i2+0,i3+1)+vxa3+vca3
           uion1(1,i1,i2+1,i3+1)=uion1(1,i1,i2+1,i3+1)+vxa4+vca4
           uion2(1,i1,i2+0,i3+0)=uion2(1,i1,i2+0,i3+0)+vxb1+vcb1
           uion2(1,i1,i2+1,i3+0)=uion2(1,i1,i2+1,i3+0)+vxb2+vcb2
           uion2(1,i1,i2+0,i3+1)=uion2(1,i1,i2+0,i3+1)+vxb3+vcb3
           uion2(1,i1,i2+1,i3+1)=uion2(1,i1,i2+1,i3+1)+vxb4+vcb4
        ENDDO
     ENDDO
  ENDDO
  ! ..counting floating point operations
  fadd   = ic * 11._real_8 + (nr1*nr2*nr3) * 23._real_8
  fmult  = ic * 12._real_8 + (nr1*nr2*nr3) * 65._real_8
  fdiv   = ic *  4._real_8 + (nr1*nr2*nr3) * 13._real_8
  fspez  =              (nr1*nr2*nr3) *  1._real_8
  flop=fadd+fmult+fdiv+fspez
  ! ==--------------------------------------------------------------==
  CALL tihalt(procedureN,isub)
  ! ==--------------------------------------------------------------==
END SUBROUTINE slypsp
! ==================================================================
SUBROUTINE xalpsp(salpha,nr1,nr2,nr3,km1,km2,km3,rho,uion1,uion2,&
     eexcu,flop)
  ! ==--------------------------------------------------------------==
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE timer,                           ONLY: tihalt,&
       tiset
  IMPLICIT NONE
  REAL(real_8)                               :: salpha
  INTEGER                                    :: nr1, nr2, nr3, km1, km2, km3
  REAL(real_8)                               :: rho(km1,km2,km3,2), &
                                                uion1(2,km1,km2,km3), &
                                                uion2(2,km1,km2,km3), eexcu, &
                                                flop

  CHARACTER(*), PARAMETER                    :: procedureN = 'xalpsp'
  REAL(real_8), PARAMETER :: eps2 = 1.e-20_real_8, &
      f1 = -1.39578858466194911_real_8 , small = 1.e-24_real_8 

  INTEGER                                    :: i1, i2, i3, ic, ii2, ii3, &
                                                isub, ncount
  REAL(real_8) :: d1, d2, d3, d4, EX1, EX2, EX3, EX4, f1s, f1u, fadd, fdiv, &
      fmult, fspez, rho1, rho2, rho3, rho4, RHOA1, RHOA2, RHOA3, RHOA4, &
      RHOB1, RHOB2, RHOB3, RHOB4, rsa1, rsa2, rsa3, rsa4, rsb1, rsb2, rsb3, &
      rsb4, VXA1, VXA2, VXA3, VXA4, VXB1, VXB2, VXB3, VXB4, x1, x2, x3, x4, &
      y1, y2, y3, y4

! ==--------------------------------------------------------------==

  CALL tiset(procedureN,isub)
  f1s=f1*salpha
  f1u=f1s*1.3333333333333333_real_8
  eexcu=0._real_8
  x1=1._real_8
  x2=1._real_8
  x3=1._real_8
  x4=1._real_8
  y1=1._real_8
  y2=1._real_8
  y3=1._real_8
  y4=1._real_8
  ic=0
  ii3=MOD(nr3,2)
  ii2=MOD(nr2,2)
  IF (ii3.NE.0.OR.ii2.NE.0) CALL stopgm('XALPSP','ODD DIMENSION',& 
       __LINE__,__FILE__)
  !$omp parallel do private(I1,I2,I3,D1,D2,D3,D4) &
  !$omp  firstprivate(X1,X2,X3,X4,Y1,Y2,Y3,Y4) &
  !$omp  private(NCOUNT,RHOA1,RHOA2,RHOA3,RHOA4,RHOB1,RHOB2,RHOB3,RHOB4) &
  !$omp  private(EX1,EX2,EX3,EX4,VXA1,VXA2,VXA3,VXA4) &
  !$omp  private(VXB1,VXB2,VXB3,VXB4) &
  !$omp  private(RHO1,RHO2,RHO3,RHO4)   &
  !$omp  private(RSA1,RSA2,RSA3,RSA4,RSB1,RSB2,RSB3,RSB4)   &
  !$omp  reduction(+:EEXCU,IC)
  DO i3=1,nr3,2
     DO i2=1,nr2,2
        DO i1=1,nr1
           rhoa1=rho(i1,i2+0,i3+0,1)-rho(i1,i2+0,i3+0,2)
           rhoa2=rho(i1,i2+1,i3+0,1)-rho(i1,i2+1,i3+0,2)
           rhoa3=rho(i1,i2+0,i3+1,1)-rho(i1,i2+0,i3+1,2)
           rhoa4=rho(i1,i2+1,i3+1,1)-rho(i1,i2+1,i3+1,2)
           rhoa1=MAX(rhoa1,small)
           rhoa2=MAX(rhoa2,small)
           rhoa3=MAX(rhoa3,small)
           rhoa4=MAX(rhoa4,small)
           rhob1=MAX(rho(i1,i2+0,i3+0,2),small)
           rhob2=MAX(rho(i1,i2+1,i3+0,2),small)
           rhob3=MAX(rho(i1,i2+0,i3+1,2),small)
           rhob4=MAX(rho(i1,i2+1,i3+1,2),small)
           rho1=rhoa1+rhob1
           rho2=rhoa2+rhob2
           rho3=rhoa3+rhob3
           rho4=rhoa4+rhob4
           ! 
           ! ..iterative calculation of rho**(1/3) (using Newton method)
           DO ncount=1,100
              d1=x1-rhoa1/(x1*x1)
              d2=x2-rhoa2/(x2*x2)
              d3=x3-rhoa3/(x3*x3)
              d4=x4-rhoa4/(x4*x4)
              ic=ic+1
              x1=x1-0.333333333333333_real_8*d1
              x2=x2-0.333333333333333_real_8*d2
              x3=x3-0.333333333333333_real_8*d3
              x4=x4-0.333333333333333_real_8*d4
              IF ((d1*d1+d2*d2+d3*d3+d4*d4).LT.eps2) GOTO 10
           ENDDO
           x1=rhoa1**(1._real_8/3._real_8)
           x2=rhoa2**(1._real_8/3._real_8)
           x3=rhoa3**(1._real_8/3._real_8)
           x4=rhoa4**(1._real_8/3._real_8)
10         CONTINUE
           rsa1=x1
           rsa2=x2
           rsa3=x3
           rsa4=x4
           DO ncount=1,100
              d1=y1-rhob1/(y1*y1)
              d2=y2-rhob2/(y2*y2)
              d3=y3-rhob3/(y3*y3)
              d4=y4-rhob4/(y4*y4)
              ic=ic+1
              y1=y1-0.333333333333333_real_8*d1
              y2=y2-0.333333333333333_real_8*d2
              y3=y3-0.333333333333333_real_8*d3
              y4=y4-0.333333333333333_real_8*d4
              IF ((d1*d1+d2*d2+d3*d3+d4*d4).LT.eps2) GOTO 20
           ENDDO
           y1=rhob1**(1._real_8/3._real_8)
           y2=rhob2**(1._real_8/3._real_8)
           y3=rhob3**(1._real_8/3._real_8)
           y4=rhob4**(1._real_8/3._real_8)
20         CONTINUE
           rsb1=y1
           rsb2=y2
           rsb3=y3
           rsb4=y4
           ! ..Slater exchange
           ex1=f1s*(rhoa1*rsa1+rhob1*rsb1)
           ex2=f1s*(rhoa2*rsa2+rhob2*rsb2)
           ex3=f1s*(rhoa3*rsa3+rhob3*rsb3)
           ex4=f1s*(rhoa4*rsa4+rhob4*rsb4)
           vxa1=f1u*rsa1
           vxa2=f1u*rsa2
           vxa3=f1u*rsa3
           vxa4=f1u*rsa4
           vxb1=f1u*rsb1
           vxb2=f1u*rsb2
           vxb3=f1u*rsb3
           vxb4=f1u*rsb4
           ! 
           eexcu=eexcu+ex1
           eexcu=eexcu+ex2
           eexcu=eexcu+ex3
           eexcu=eexcu+ex4
           uion1(1,i1,i2+0,i3+0)=uion1(1,i1,i2+0,i3+0)+vxa1
           uion1(1,i1,i2+1,i3+0)=uion1(1,i1,i2+1,i3+0)+vxa2
           uion1(1,i1,i2+0,i3+1)=uion1(1,i1,i2+0,i3+1)+vxa3
           uion1(1,i1,i2+1,i3+1)=uion1(1,i1,i2+1,i3+1)+vxa4
           uion2(1,i1,i2+0,i3+0)=uion2(1,i1,i2+0,i3+0)+vxb1
           uion2(1,i1,i2+1,i3+0)=uion2(1,i1,i2+1,i3+0)+vxb2
           uion2(1,i1,i2+0,i3+1)=uion2(1,i1,i2+0,i3+1)+vxb3
           uion2(1,i1,i2+1,i3+1)=uion2(1,i1,i2+1,i3+1)+vxb4
        ENDDO
     ENDDO
  ENDDO
  ! ..counting floating point operations
  fadd   = ic * 11._real_8 + (nr1*nr2*nr3) *  5._real_8
  fmult  = ic * 12._real_8 + (nr1*nr2*nr3) *  5._real_8
  fdiv   = ic *  4._real_8 + (nr1*nr2*nr3) *  0._real_8
  fspez  =              (nr1*nr2*nr3) *  0._real_8
  flop=fadd+fmult+fdiv+fspez
  ! ==--------------------------------------------------------------==
  CALL tihalt(procedureN,isub)
  ! ==--------------------------------------------------------------==
END SUBROUTINE xalpsp
! ==================================================================
#else
! ==================================================================
SUBROUTINE lypuu(nr1,nr2,nr3,km1,km2,km3,rho,uion,eexcu,flop)
  ! ==--------------------------------------------------------------==
  ! adds unpolarized exc potential to uion and sums exc energy
  ! written by S. Goedecker, Stuttgart, march 1996
  ! modified by J. Hutter, Stuttgart, may/august 1996
  ! non-gradient part of Lee,Yang and Parr functional (non-spin polarised)
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE timer,                           ONLY: tihalt,&
       tiset
  IMPLICIT NONE
  INTEGER                                    :: nr1, nr2, nr3, km1, km2, km3
  REAL(real_8)                               :: rho(km1*km2*km3), &
                                                uion(2,km1*km2*km3), eexcu, &
                                                flop

  CHARACTER(*), PARAMETER                    :: procedureN = 'lypuu'
  REAL(real_8), PARAMETER :: a = 0.04918_real_8, b = 0.132_real_8, &
      c = 0.2533_real_8, cf = 2.87123400018819108_real_8 , d = 0.349_real_8 , &
      eps2 = 1.e-20_real_8, f1 = -1.10783814957303361_real_8 , &
      f43 = 4._real_8/3._real_8 , small = 1.e-24_real_8 , &
      third = 1._real_8/3._real_8

  INTEGER                                    :: i, ic, isub, kkk
  REAL(real_8)                               :: bcf, ecrs, EEXC1X, elyp, &
                                                EPSXCU, fa1, fa2, fadd, fdiv, &
                                                fmult, fspez, ox, rhou, rsx, &
                                                VXC, x

! ==--------------------------------------------------------------==

  CALL tiset(procedureN,isub)
  ! ==--------------------------------------------------------------==
  fa1=f1*0.66666666666666_real_8
  fa2=fa1*f43
  bcf=b*cf
  eexcu=0._real_8
  ic=0
  kkk=km1*km2*km3
  !$omp parallel do private(I,RHOU,X,RSX,ECRS,OX,ELYP,EPSXCU,VXC,EEXC1X) &
  !$omp  reduction(+:EEXCU,IC)
  DO i=1,kkk
     rhou=ABS(rho(i))
     IF (rhou.GT.small)THEN
        x=rhou**(0.333333333333333_real_8)
        rsx=1._real_8/x
        ecrs=bcf*EXP(-c*rsx)
        ox=1._real_8/(1._real_8+d*rsx)
        elyp=-a*ox*(1._real_8+ecrs)
        epsxcu=fa1*x+elyp
        vxc=fa2*x+elyp-0.33333333333333_real_8*rsx*a*ox*&
             (d*ox+ecrs*(d*ox+c))
        uion(1,i)=uion(1,i)+vxc
        ic = ic + 1
        eexcu = eexcu + epsxcu*rhou
     ENDIF
  ENDDO
  ! ..counting floating point operations
  fadd   = ic *  9._real_8
  fmult  = ic * 15._real_8
  fdiv   = ic *  2._real_8
  fspez  = ic *  2._real_8
  flop=fadd+fmult+fdiv+fspez
  ! ==--------------------------------------------------------------==
  CALL tihalt(procedureN,isub)
  ! ==--------------------------------------------------------------==
END SUBROUTINE lypuu
! ==================================================================
SUBROUTINE xalpu(salpha,nr1,nr2,nr3,km1,km2,km3,rho,uion,&
     eexcu,flop)
  ! ==--------------------------------------------------------------==
  ! ..X-Alpha potential
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE timer,                           ONLY: tihalt,&
       tiset
  IMPLICIT NONE
  REAL(real_8)                               :: salpha
  INTEGER                                    :: nr1, nr2, nr3, km1, km2, km3
  REAL(real_8)                               :: rho(km1*km2*km3), &
                                                uion(2,km1*km2*km3), eexcu, &
                                                flop

  CHARACTER(*), PARAMETER                    :: procedureN = 'xalpu'
  REAL(real_8), PARAMETER :: eps2 = 1.e-20_real_8, &
      f1 = -1.10783814957303361_real_8 , small = 1.e-24_real_8 

  INTEGER                                    :: i, ic, isub, kkk
  REAL(real_8)                               :: eexc1x, epsxcu, fa1, fa2, &
                                                fadd, fdiv, fmult, fspez, &
                                                rhou, vxc, x

! ==--------------------------------------------------------------==

  CALL tiset(procedureN,isub)
  ! ==--------------------------------------------------------------==
  fa1=f1*salpha
  fa2=fa1*1.33333333333333_real_8
  eexcu=0._real_8
  ic=0
  kkk=km1*km2*km3
  !$omp parallel do private(I,RHOU,X,EPSXCU,VXC,EEXC1X) &
  !$omp  reduction(+:EEXCU,IC)
  DO i=1,kkk
     rhou=ABS(rho(i))
     IF (rhou.GT.small)THEN
        x=rhou**(0.33333333333333_real_8)
        epsxcu=fa1*x
        vxc=fa2*x
        uion(1,i)=uion(1,i)+vxc
        ic=ic+1
        eexcu=eexcu+epsxcu*rhou
     ENDIF
  ENDDO
  ! ..counting floating point operations
  fadd   = ic *  2._real_8
  fmult  = ic *  3._real_8
  fdiv   = ic *  0._real_8
  fspez  = ic *  1._real_8
  flop=fadd+fmult+fdiv+fspez
  ! ==--------------------------------------------------------------==
  CALL tihalt(procedureN,isub)
  ! ==--------------------------------------------------------------==
END SUBROUTINE xalpu
! ==================================================================
SUBROUTINE mikeu(nr1,nr2,nr3,km1,km2,km3,rho,uion,eexcu,flop)
  ! ==--------------------------------------------------------------==
  ! ..Pade approximation to LDA functional
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE timer,                           ONLY: tihalt,&
       tiset
  IMPLICIT NONE
  INTEGER                                    :: nr1, nr2, nr3, km1, km2, km3
  REAL(real_8)                               :: rho(km1*km2*km3), &
                                                uion(2,km1*km2*km3), eexcu, &
                                                flop

  CHARACTER(*), PARAMETER                    :: procedureN = 'mikeu'
  REAL(real_8), PARAMETER :: a0u = 0.4581652932831429_real_8, &
      a1u = 2.217058676663745_real_8, a2u = 0.7405551735357053_real_8, &
      a3u = 0.01968227878617998_real_8 , b1u = 1.0_real_8, &
      b2u = 4.504130959426697_real_8, b3u = 1.110667363742916_real_8, &
      b4u = 0.02359291751427506_real_8 , c1u = 4._real_8*a0u*b1u/3._real_8, &
      c2u = 5._real_8*a0u*b2u/3._real_8+a1u*b1u, c3u = 2._real_8*a0u*&
      b3u+4._real_8*a1u*b2u/3._real_8+2._real_8*a2u*b1u/3._real_8, c4u = &
      7._real_8*a0u*b4u/3._real_8+5._real_8*a1u*b3u/3._real_8+a2u*b2u+a3u*b1u/&
      3._real_8
  REAL(real_8), PARAMETER :: c5u = 2._real_8*a1u*b4u+4._real_8*a2u*b3u/&
      3._real_8+2._real_8*a3u*b2u/3._real_8, &
      c6u = 5._real_8*a2u*b4u/3._real_8+a3u*b3u, &
      c7u = 4._real_8*a3u*b4u/3._real_8 , eps2 = 1.e-20_real_8, &
      rsfac = 0.6203504908994000_real_8 , small = 1.e-24_real_8 

  INTEGER                                    :: i, ic, isub, kkk
  REAL(real_8)                               :: botu, dtopu, eexc1x, epsxcu, &
                                                fadd, fdiv, fmult, fspez, &
                                                rhou, rsu, topu, vxc

! ==--------------------------------------------------------------==

  CALL tiset(procedureN,isub)
  ! ==--------------------------------------------------------------==
  eexcu=0._real_8
  ic=0
  kkk=km1*km2*km3
  !$omp parallel do private(I,RHOU,RSU,TOPU,DTOPU,BOTU,EPSXCU,VXC,EEXC1X) &
  !$omp  reduction(+:EEXCU,IC)
  DO i=1,kkk
     rhou=ABS(rho(i))
     IF (rhou.GT.small)THEN
        rsu=rsfac*rhou**(-0.33333333333333_real_8)
        topu=a0u+rsu*(a1u+rsu*(a2u+rsu*a3u))
        dtopu=-rsu*(c1u+rsu*(c2u+rsu*(c3u+rsu*&
             (c4u+rsu*(c5u+rsu*(c6u+rsu*c7u))))))
        botu=1._real_8/(rsu*(b1u+rsu*(b2u+rsu*(b3u+rsu*b4u))))
        epsxcu=-topu*botu
        vxc=dtopu*botu*botu
        uion(1,i)=uion(1,i)+vxc
        ic=ic+1
        eexcu=eexcu+epsxcu*rhou
     ENDIF
  ENDDO
  ! ..counting floating point operations
  fadd   = ic * 14._real_8
  fmult  = ic * 19._real_8
  fdiv   = ic *  1._real_8
  fspez  = ic *  1._real_8
  flop=fadd+fmult+fdiv+fspez
  ! ==--------------------------------------------------------------==
  CALL tihalt(procedureN,isub)
  ! ==--------------------------------------------------------------==
END SUBROUTINE mikeu
! ==================================================================
SUBROUTINE mikesp(nr1,nr2,nr3,km1,km2,km3,rho,uion1,uion2,&
     eexcu,flop)
  ! ==--------------------------------------------------------------==
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE timer,                           ONLY: tihalt,&
       tiset
  IMPLICIT NONE
  INTEGER                                    :: nr1, nr2, nr3, km1, km2, km3
  REAL(real_8)                               :: rho(km1*km2*km3,2), &
                                                uion1(2,km1*km2*km3), &
                                                uion2(2,km1*km2*km3), eexcu, &
                                                flop

  CHARACTER(*), PARAMETER                    :: procedureN = 'mikesp'
  REAL(real_8), PARAMETER :: a0u = 0.4581652932831429_real_8, &
      a1u = 2.217058676663745_real_8, a2u = 0.7405551735357053_real_8, &
      a3u = 0.01968227878617998_real_8 , b1u = 1.0_real_8, &
      b2u = 4.504130959426697_real_8, b3u = 1.110667363742916_real_8, &
      b4u = 0.02359291751427506_real_8 , da0 = 0.119086804055547_real_8, &
      da1 = 0.6157402568883345_real_8, da2 = 0.1574201515892867_real_8, &
      da3 = 0.003532336663397157_real_8 , db1 = 0.0_real_8, &
      db2 = 0.2673612973836267_real_8, db3 = 0.2052004607777787_real_8, &
      db4 = 0.004200005045691381_real_8 , eps2 = 1.e-20_real_8, &
      rsfac = 0.6203504908994000_real_8 
  REAL(real_8), PARAMETER :: small = 1.e-24_real_8 , &
      t1 = 0.854960467080683_real_8, t2 = 0.07916300621117439_real_8, &
      t3 = 0.02580127609845684_real_8, t4 = 0.0121839359353824_real_8, &
      t5 = 0.006919272259599879_real_8, t6 = 0.004391524649611372_real_8, &
      t7 = 0.003002751897170168_real_8, t8 = 0.00216587382212552_real_8, &
      t9 = 0.01141189204579585_real_8 

  INTEGER                                    :: i, ic, isub
  REAL(real_8) :: a0p, a1p, a2p, a3p, B1P, B2P, B3P, B4P, BOT, DBOT, dfxs, &
      DTOP, DVXC, DXA, DXB, EEXC1X, EPSXCU, et, eta, fadd, fdiv, fmult, &
      fspez, fxs, OBOT, rhoa, rhob, rhou, RSU, TOP, VXC, X, XBOT, XTOP

! ==--------------------------------------------------------------==

  CALL tiset(procedureN,isub)
  ! ==--------------------------------------------------------------==
  eexcu=0._real_8
  ic=0
  !$omp parallel do private(I,RHOA,RHOB,RHOU,ETA,ET,FXS) &
  !$omp  private(DFXS,A0P,A1P,A2P,A3P,B1P,B2P,B3P,B4P) &
  !$omp  private(X,RSU,TOP,DTOP,XTOP,BOT,DBOT) &
  !$omp  private(XBOT,OBOT,EPSXCU,VXC,DVXC,DXA,DXB,EEXC1X) &
  !$omp  reduction(+:EEXCU,IC)
  DO i=1,km1*km2*km3
     rhoa=MAX(rho(i,1)-rho(i,2),small)
     rhob=MAX(rho(i,2),small)
     rhou=rhoa+rhob
     IF (rhou.GT.2.25_real_8*small) THEN
        eta=(rhoa-rhob)/rhou
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
        x=rhou**(0.333333333333_real_8)
        rsu=rsfac/x
        top=a0p+rsu*(a1p+rsu*(a2p+rsu*a3p))
        dtop= a1p+rsu*(2._real_8*a2p+rsu*3._real_8*a3p)
        xtop=da0+rsu*(da1+rsu*(da2+rsu*da3))
        bot=rsu*(b1p+rsu*(b2p+rsu*(b3p+rsu*b4p)))
        dbot=b1p+rsu*(2._real_8*b2p+rsu*(3._real_8*b3p+rsu*4._real_8*b4p))
        xbot=rsu*(db1+rsu*(db2+rsu*(db3+rsu*db4)))
        obot=1._real_8/bot
        epsxcu=-top*obot
        ic=ic+1
        eexcu=eexcu+epsxcu*rhou
        vxc=epsxcu+rsu*(dtop*obot-dbot*(top*obot*obot))*&
             0.333333333333333_real_8
        dvxc=-(xtop*obot-xbot*(top*obot*obot))
        dxa=vxc+dvxc*dfxs*(1._real_8-eta)
        dxb=vxc-dvxc*dfxs*(1._real_8+eta)
        uion1(1,i)=uion1(1,i)+dxa
        uion2(1,i)=uion2(1,i)+dxb
     ENDIF
  ENDDO
  ! ..counting floating point operations
  fadd   = ic * 55._real_8
  fmult  = ic * 76._real_8
  fdiv   = ic *  3._real_8
  fspez  = ic *  1._real_8
  flop=fadd+fmult+fdiv+fspez
  ! ==--------------------------------------------------------------==
  CALL tihalt(procedureN,isub)
  ! ==--------------------------------------------------------------==
END SUBROUTINE mikesp
! ==================================================================
SUBROUTINE slypsp(nr1,nr2,nr3,km1,km2,km3,rho,uion1,uion2,salpha,&
     eexcu,flop)
  ! ==--------------------------------------------------------------==
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE timer,                           ONLY: tihalt,&
       tiset
  IMPLICIT NONE
  INTEGER                                    :: nr1, nr2, nr3, km1, km2, km3
  REAL(real_8)                               :: rho(km1*km2*km3,2), &
                                                uion1(2,km1*km2*km3), &
                                                uion2(2,km1*km2*km3), salpha, &
                                                eexcu, flop

  CHARACTER(*), PARAMETER                    :: procedureN = 'slypsp'
  REAL(real_8), PARAMETER :: a = 0.04918_real_8, b = 0.132_real_8, &
      c = 0.2533_real_8, cf = 2.87123400018819108_real_8 , d = 0.349_real_8, &
      eps2 = 1.e-20_real_8, f1 = -1.39578858466194911_real_8 , &
      small = 1.e-24_real_8 

  INTEGER                                    :: i, ic, isub, kkk
  REAL(real_8) :: cfx, DD, DOR, EC, ECR, EEXC1X, EX, f1s, f1u, fadd, fdiv, &
      fmult, fspez, OR, orho, ors, rhoa, RHOAT, rhob, RHOBT, rhot, rs, rsa, &
      rsb, VCA, VCB, VXA, VXB

! ==--------------------------------------------------------------==

  CALL tiset(procedureN,isub)
  cfx=cf*2._real_8**3.66666666666666_real_8
  f1s=f1*salpha
  f1u=f1s*1.3333333333333333_real_8
  eexcu=0._real_8
  ic=0
  kkk=km1*km2*km3
  !$omp parallel do private(I,RHOA,RHOB,RHOT,RSA,RSB,RS,ORS,ORHO,EX) &
  !$omp  private(VXA,VXB,DD,ECR,OR,DOR,RHOAT,RHOBT,EC,VCA,VCB,EEXC1X) &
  !$omp  reduction(+:EEXCU,IC)
  DO i=1,kkk
     rhoa=rho(i,1)-rho(i,2)
     rhoa=MAX(rhoa,small)
     rhob=MAX(rho(i,2),small)
     rhot=rhoa+rhob
     IF (rhot.GT.2._real_8*small)THEN
        rhoa=MAX(rhoa,small)
        rhob=MAX(rhob,small)
        rsa=rhoa**(0.3333333333333_real_8)
        rsb=rhob**(0.3333333333333_real_8)
        rs=rhot**(0.3333333333333_real_8)
        ors=1._real_8/rs
        orho=1._real_8/rhot
        ! ..Slater exchange
        ex=f1s*(rhoa*rsa+rhob*rsb)
        vxa=f1u*rsa
        vxb=f1u*rsb
        ! ..Lee Yang Parr
        dd=1._real_8/(1._real_8+d*ors)
        ecr=EXP(-c*ors)
        or=ecr*dd*orho*orho*orho*ors*ors
        dor=-or*0.333333333333333_real_8*orho*ors*(11._real_8*rs-c-d*dd)
        rhoat=rhoa*rhoa*rsa*rsa
        rhobt=rhob*rhob*rsb*rsb
        ec=-4._real_8*a*rhoa*rhob*orho*dd-cfx*a*b*or*rhoa*rhob*&
             (rhoat+rhobt)
        vca=-4._real_8*a*dd*rhoa*rhob*orho*(d*0.333333333333333_real_8*dd*&
             (orho*ors)+(1._real_8/rhoa)-orho)-cfx*a*b*(dor*&
             rhoa*rhob*(rhoat+rhobt)+or*rhob*&
             (3.6666666666666667_real_8*rhoat+rhobt))
        vcb=-4._real_8*a*dd*rhoa*rhob*orho*(d*0.333333333333333_real_8*dd*&
             (orho*ors)+(1._real_8/rhob)-orho)-cfx*a*b*(dor*&
             rhoa*rhob*(rhoat+rhobt)+or*rhoa*&
             (3.6666666666666667_real_8*rhobt+rhoat))
        uion1(1,i)=uion1(1,i)+vxa+vca
        uion2(1,i)=uion2(1,i)+vxb+vcb
        ic=ic+1
        eexcu=eexcu+ex+ec
     ENDIF
  ENDDO
  ! ..counting floating point operations
  fadd   = ic * 24._real_8
  fmult  = ic * 74._real_8
  fdiv   = ic *  5._real_8
  fspez  = ic *  4._real_8
  flop=fadd+fmult+fdiv+fspez
  ! ==--------------------------------------------------------------==
  CALL tihalt(procedureN,isub)
  ! ==--------------------------------------------------------------==
END SUBROUTINE slypsp
! ==================================================================
SUBROUTINE xalpsp(salpha,nr1,nr2,nr3,km1,km2,km3,rho,uion1,uion2,&
     eexcu,flop)
  ! ==--------------------------------------------------------------==
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE timer,                           ONLY: tihalt,&
       tiset
  IMPLICIT NONE
  REAL(real_8)                               :: salpha
  INTEGER                                    :: nr1, nr2, nr3, km1, km2, km3
  REAL(real_8)                               :: rho(km1*km2*km3,2), &
                                                uion1(2,km1*km2*km3), &
                                                uion2(2,km1*km2*km3), eexcu, &
                                                flop

  CHARACTER(*), PARAMETER                    :: procedureN = 'xalpsp'
  REAL(real_8), PARAMETER :: eps2 = 1.e-20_real_8, &
      f1 = -1.39578858466194911_real_8 , small = 1.e-24_real_8 

  INTEGER                                    :: i, ic, isub, kkk
  REAL(real_8)                               :: eexc1x, ex, f1s, f1u, fadd, &
                                                fdiv, fmult, fspez, rhoa, &
                                                rhob, rhot, rsa, rsb, vxa, vxb

! ==--------------------------------------------------------------==

  CALL tiset(procedureN,isub)
  ! ==--------------------------------------------------------------==
  f1s=f1*salpha
  f1u=f1s*1.3333333333333333_real_8
  eexcu=0._real_8
  ic=0
  kkk=km1*km2*km3
  !$omp parallel do private(I,RHOA,RHOB,RHOT,RSA,RSB,EX,VXA,VXB,EEXC1X) &
  !$omp  reduction(+:EEXCU,IC)
  DO i=1,kkk
     rhoa=rho(i,1)-rho(i,2)
     rhoa=MAX(rhoa,small)
     rhob=MAX(rho(i,2),small)
     rhot=rhoa+rhob
     IF (rhot.GT.2._real_8*small) THEN
        rsa=rhoa**(0.33333333333333_real_8)
        rsb=rhob**(0.33333333333333_real_8)
        ! ..Slater exchange
        ex=f1s*(rhoa*rsa+rhob*rsb)
        vxa=f1u*rsa
        vxb=f1u*rsb
        uion1(1,i)=uion1(1,i)+vxa
        uion2(1,i)=uion2(1,i)+vxb
        ic=ic+1
        eexcu=eexcu+ex
     ENDIF
  ENDDO
  ! ..counting floating point operations
  fadd   = ic *  5._real_8
  fmult  = ic *  5._real_8
  fdiv   = ic *  0._real_8
  fspez  = ic *  2._real_8
  flop=fadd+fmult+fdiv+fspez
  ! ==--------------------------------------------------------------==
  CALL tihalt(procedureN,isub)
  ! ==--------------------------------------------------------------==
END SUBROUTINE xalpsp
! ==================================================================
#endif
