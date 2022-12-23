MODULE aainit_utils
  USE aavan,                           ONLY: ap,&
                                             lpl,&
                                             lpx
  USE cnst,                            ONLY: fpi
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: paral
  USE system,                          ONLY: lix,&
                                             lx,&
                                             mix,&
                                             mx,&
                                             nlx
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: aainit

CONTAINS

  ! ==================================================================
  SUBROUTINE aainit(LLi)
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    ! NLX= COMBINED ANGULAR MOMENTUM FOR LLi
    ! (FOR S,P AND D STATES NLX=9)
    ! LX = MAX 2*LLi-1
    ! 25 = COMBINED ANGULAR MOMENTUM FOR 2*LLi-1 (FOR S,P AND D STATES)
    ! LiX= MAX LLi
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: LLi

    COMPLEX(real_8)                          :: aa(lx,mx,LiX,MiX,lix,mix), &
                                                sum, u(lx,mx,mx)
    INTEGER                                  :: ii, ik, il, ilp, k, l, &
                                                length, Li, Lj, ll, lp, m, &
                                                Mi, Mj, mp, mpp, n, np, p
    REAL(real_8)                             :: cc(LiX,MiX,lix,mix,lx), f, &
                                                fs, fts, s, ss, t, ts, y

! Variables
! ==--------------------------------------------------------------==
! ARRAY AA
! 
! CC = < LM | LiMi | LjMj > = sqrt(2L+1/4pi) (-1)^Mj c^L(Li-Mi,LjMj)
! 
! COEFF. c^L FROM WEISSBLUTH 'ATOMS AND MOLEC..'  PAGE 246
! CC(Li,Mi,Lj,Mj,L)  L  IS FREE INDEX, M=Mj+Mj   (INPUT LiMi LjMj)
! ==--------------------------------------------------------------==

    ll=2*LLi-1
    IF (LLi.GT.LiX) CALL stopgm(' AAINIT ',' LLi .GT. LiX ',& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    length=LiX*MiX*lix*mix*lx
    CALL zeroing(cc)!,length)
    ! LiMi  LjMj   -   LM
    cc(1,1,1,1,1)= 1._real_8       ! s     s     -   s

    IF (LLi.GT.1) THEN
       cc(1,1,2,1,2)= 1._real_8   ! s    p_-1   -  p_-1
       cc(1,1,2,2,2)= 1._real_8   ! s    p_0    -  p_0
       cc(1,1,2,3,2)= 1._real_8   ! s    p_1    -  p_1

       cc(2,1,2,1,3)= SQRT(6._real_8/5._real_8)! p_-1  p_-1   -  d_-2
       cc(2,1,2,2,3)= SQRT(3._real_8/5._real_8)! p_-1  p_0    -  d_-1
       cc(2,1,2,3,1)=-1._real_8      ! p_-1  p_1    -   s
       cc(2,1,2,3,3)= SQRT(1._real_8/5._real_8)! p_-1  p_1    -  d_0
       cc(2,2,2,2,1)= 1._real_8      ! p_0   p_0    -   s
       cc(2,2,2,2,3)= SQRT(4._real_8/5._real_8)! p_0   p_0    -  d_0
       cc(2,2,2,3,3)= SQRT(3._real_8/5._real_8)! p_0   p_1    -  d_1
       cc(2,3,2,3,3)= SQRT(6._real_8/5._real_8)! p_1   p_1    -  d_2
    ENDIF

    IF (LLi.GT.2) THEN
       cc(1,1,3,1,3)= 1._real_8   ! s    d_-2   -  d_-2
       cc(1,1,3,2,3)= 1._real_8   ! s    d_-1   -  d_-1
       cc(1,1,3,3,3)= 1._real_8   ! s    d_0    -  d_0
       cc(1,1,3,4,3)= 1._real_8   ! s    d_1    -  d_1
       cc(1,1,3,5,3)= 1._real_8   ! s    d_2    -  d_2

       cc(2,1,3,1,4)= SQRT(45._real_8/35._real_8)! p_-1  d_-2  -  f_-3
       cc(2,1,3,2,4)= SQRT(30._real_8/35._real_8)! p_-1  d_-1  -  f_-2
       cc(2,1,3,3,2)=-SQRT( 1._real_8/ 5._real_8)! p_-1  d_0   -  p_-1
       cc(2,1,3,3,4)= SQRT(18._real_8/35._real_8)! p_-1  d_0   -  f_-1
       cc(2,1,3,4,2)=-SQRT( 3._real_8/ 5._real_8)! p_-1  d_1   -  p_0
       cc(2,1,3,4,4)= SQRT( 9._real_8/35._real_8)! p_-1  d_1   -  f_0
       cc(2,1,3,5,2)=-SQRT( 6._real_8/ 5._real_8)! p_-1  d_2   -  p_1
       cc(2,1,3,5,4)= SQRT( 3._real_8/35._real_8)! p_-1  d_2   -  f_1
       cc(2,2,3,1,4)= SQRT(15._real_8/35._real_8)! p_0   d_-2  -  f_-2
       cc(2,2,3,2,2)= SQRT( 3._real_8/ 5._real_8)! p_0   d_-1  -  p_-1
       cc(2,2,3,2,4)= SQRT(24._real_8/35._real_8)! p_0   d_-1  -  f_-1
       cc(2,2,3,3,2)= SQRT( 4._real_8/ 5._real_8)! p_0   d_0   -  p_0
       cc(2,2,3,3,4)= SQRT(27._real_8/35._real_8)! p_0   d_0   -  f_0
       cc(2,2,3,4,2)= SQRT( 3._real_8/ 5._real_8)! p_0   d_1   -  p_1
       cc(2,2,3,4,4)= SQRT(24._real_8/35._real_8)! p_0   d_1   -  f_1
       cc(2,2,3,5,4)= SQRT(15._real_8/35._real_8)! p_0   d_2   -  f_2
       cc(2,3,3,1,2)=-SQRT( 6._real_8/ 5._real_8)! p_1   d_-2  -  p_-1
       cc(2,3,3,1,4)= SQRT( 3._real_8/35._real_8)! p_1   d_-2  -  f_-1
       cc(2,3,3,2,2)=-SQRT( 3._real_8/ 5._real_8)! p_1   d_-1  -  p_0
       cc(2,3,3,2,4)= SQRT( 9._real_8/35._real_8)! p_1   d_-1  -  f_0
       cc(2,3,3,3,2)=-SQRT( 1._real_8/ 5._real_8)! p_1   d_0   -  p_1
       cc(2,3,3,3,4)= SQRT(18._real_8/35._real_8)! p_1   d_0   -  f_1
       cc(2,3,3,4,4)= SQRT(30._real_8/35._real_8)! p_1   d_1   -  f_2
       cc(2,3,3,5,4)= SQRT(45._real_8/35._real_8)! p_1   d_2   -  f_3

       cc(3,1,3,1,5)= SQRT(70._real_8/49._real_8)! d_-2  d_-2  -  g_-4
       cc(3,1,3,2,5)= SQRT(35._real_8/49._real_8)! d_-2  d_-1  -  g_-3
       cc(3,1,3,3,3)=-SQRT(20._real_8/49._real_8)! d_-2  d_0   -  d_-2
       cc(3,1,3,3,5)= SQRT(15._real_8/49._real_8)! d_-2  d_0   -  g_-2
       cc(3,1,3,4,3)=-SQRT(30._real_8/49._real_8)! d_-2  d_1   -  d_1
       cc(3,1,3,4,5)= SQRT( 5._real_8/49._real_8)! d_-2  d_1   -  g_1
       cc(3,1,3,5,1)= 1._real_8             ! d_-2  d_2   -  s
       cc(3,1,3,5,3)=-SQRT(20._real_8/49._real_8)! d_-2  d_2   -  d_0
       cc(3,1,3,5,5)= SQRT( 1._real_8/49._real_8)! d_-2  d_2   -  g_0

       cc(3,2,3,2,3)= SQRT(30._real_8/49._real_8)! d_-1  d_-1  -  d_-2
       cc(3,2,3,2,5)= SQRT(40._real_8/49._real_8)! d_-1  d_-1  -  g_-2
       cc(3,2,3,3,3)= SQRT( 5._real_8/49._real_8)! d_-1  d_0   -  d_-1
       cc(3,2,3,3,5)= SQRT(30._real_8/49._real_8)! d_-1  d_0   -  g_-1
       cc(3,2,3,4,1)=-1._real_8             ! d_-1  d_1   -  s
       cc(3,2,3,4,3)=-SQRT( 5._real_8/49._real_8)! d_-1  d_1   -  d_0
       cc(3,2,3,4,5)= SQRT(16._real_8/49._real_8)! d_-1  d_1   -  g_0
       cc(3,2,3,5,3)=-SQRT(30._real_8/49._real_8)! d_-1  d_2   -  d_1
       cc(3,2,3,5,5)= SQRT( 5._real_8/49._real_8)! d_-1  d_2   -  g_1

       cc(3,3,3,3,1)= 1._real_8             ! d_0   d_0   -  s
       cc(3,3,3,3,3)= SQRT(20._real_8/49._real_8)! d_0   d_0   -  d_0
       cc(3,3,3,3,5)= SQRT(36._real_8/49._real_8)! d_0   d_0   -  g_0
       cc(3,3,3,4,3)= SQRT( 5._real_8/49._real_8)! d_0   d_1   -  d_1
       cc(3,3,3,4,5)= SQRT(30._real_8/49._real_8)! d_0   d_1   -  g_1
       cc(3,3,3,5,3)=-SQRT(20._real_8/49._real_8)! d_0   d_2   -  d_2
       cc(3,3,3,5,5)= SQRT(15._real_8/49._real_8)! d_0   d_2   -  g_2

       cc(3,4,3,4,3)= SQRT(30._real_8/49._real_8)! d_1   d_1   -  d_2
       cc(3,4,3,4,5)= SQRT(40._real_8/49._real_8)! d_1   d_1   -  g_2
       cc(3,4,3,5,5)= SQRT(35._real_8/49._real_8)! d_1   d_2   -  g_3

       cc(3,5,3,5,5)= SQRT(70._real_8/49._real_8)! d_2   d_2   -  g_4
    ENDIF
    ! ==--------------------------------------------------------------==
    DO l=1,ll
       DO Li=1,LLi
          DO Mi=1,2*Li-1 ! SYMMETRY   LiMi  <-> LjMj
             DO Lj=1,Li-1
                DO Mj=1,2*Lj-1
                   cc(Li,Mi,Lj,Mj,l)=cc(lj,mj,li,mi,l)
                ENDDO
             ENDDO
             DO Mj=1,Mi-1
                cc(Li,Mi,li,Mj,l)=cc(li,mj,li,mi,l)! Li = Lj
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    CALL dscal(length,SQRT(1._real_8/fpi),cc,1)
    ! ==--------------------------------------------------------------==
    ! TRANSFORM BETWEEN REAL SPHERICAL HARMONICS AND THE ORIGINAL ONES
    ! (  Y^R_lm = Sum_n U(l,m,n) Y_ln )
    ! 
    ! U(l,m,n)  see  Weissbluth  pages 128 - 130
    ! errors of Weissbluth have been corrected:
    ! 1.)   i/2 has been changed to i/4 for u(4,6,*)
    ! 2.)  -i/4 has been changed to i/4 for u(5,8,*)
    ! ==--------------------------------------------------------------==
    CALL zeroing(u)!,lx*mx*mx)
    u(1,1,1)=CMPLX(1._real_8,0._real_8,kind=real_8) ! L = 0

    IF (LLi.GT.1) THEN
       y=1._real_8/SQRT(2._real_8)
       u(2,1,1)=CMPLX(y,0._real_8,kind=real_8)  ! L = 1    ORDERING   X   Z   Y
       u(2,1,3)=CMPLX(-y,0._real_8,kind=real_8) ! M     -1   0   1
       u(2,2,2)=CMPLX(1._real_8,0._real_8,kind=real_8)! M      1   2   3
       u(2,3,1)=CMPLX(0._real_8,y,kind=real_8)
       u(2,3,3)=CMPLX(0._real_8,y,kind=real_8)

       u(3,1,1)=CMPLX(0._real_8,y,kind=real_8) ! L = 2
       u(3,1,5)=CMPLX(0._real_8,-y,kind=real_8)! ORDERING  XY   XZ  Z^2  YZ  X^2-Y^2
       u(3,2,2)=CMPLX(y,0._real_8,kind=real_8) ! M      -2   -1   0    1     2
       u(3,2,4)=CMPLX(-y,0._real_8,kind=real_8)! M       1    2   3    4     5
       u(3,3,3)=CMPLX(1._real_8,0._real_8,kind=real_8)
       u(3,4,2)=CMPLX(0._real_8,y,kind=real_8)
       u(3,4,4)=CMPLX(0._real_8,y,kind=real_8)
       u(3,5,1)=CMPLX(y,0._real_8,kind=real_8)
       u(3,5,5)=CMPLX(y,0._real_8,kind=real_8)
    ENDIF

    IF (LLi.GT.2) THEN
       f=SQRT(5._real_8)/4._real_8
       t=SQRT(3._real_8)/4._real_8
       u(4,1,1)=CMPLX(f,0._real_8,kind=real_8)! L = 3
       u(4,1,3)=CMPLX(-t,0._real_8,kind=real_8)! ORDER
       u(4,1,5)=CMPLX(t,0._real_8,kind=real_8)! X  Y  XYZ  Z  Z(X^2-Y^2) Y(Z^2-X^2) X(Y^2-Z^2)
       u(4,1,7)=CMPLX(-f,0._real_8,kind=real_8)! 1  2   3   4       5          6          7
       u(4,2,1)=CMPLX(0._real_8,-f,kind=real_8)
       u(4,2,3)=CMPLX(0._real_8,-t,kind=real_8)
       u(4,2,5)=CMPLX(0._real_8,-t,kind=real_8)
       u(4,2,7)=CMPLX(0._real_8,-f,kind=real_8)
       u(4,3,2)=CMPLX(0._real_8, y,kind=real_8)
       u(4,3,6)=CMPLX(0._real_8,-y,kind=real_8)
       u(4,4,4)=CMPLX(1._real_8,0._real_8,kind=real_8)
       u(4,5,2)=CMPLX(y,0._real_8,kind=real_8)
       u(4,5,6)=CMPLX(y,0._real_8,kind=real_8)
       u(4,6,1)=CMPLX(0._real_8,-t,kind=real_8)! MISSPRINT IN WEISSBLUTH (PAGE 129)
       u(4,6,3)=CMPLX(0._real_8, f,kind=real_8)! MISSPRINT IN WEISSBLUTH (PAGE 129)
       u(4,6,5)=CMPLX(0._real_8, f,kind=real_8)! MISSPRINT IN WEISSBLUTH (PAGE 129)
       u(4,6,7)=CMPLX(0._real_8,-t,kind=real_8)! MISSPRINT IN WEISSBLUTH (PAGE 129)
       u(4,7,1)=CMPLX(-t,0._real_8,kind=real_8)
       u(4,7,3)=CMPLX(-f,0._real_8,kind=real_8)
       u(4,7,5)=CMPLX( f,0._real_8,kind=real_8)
       u(4,7,7)=CMPLX( t,0._real_8,kind=real_8)

       s  =SQRT( 7._real_8     )/4.0_real_8
       fs =SQRT( 5._real_8/6._real_8)/2.0_real_8
       ss =SQRT( 7._real_8/6._real_8)/2.0_real_8
       fts=SQRT(14._real_8/6._real_8)/2.0_real_8
       ts =SQRT(10._real_8/6._real_8)/2.0_real_8
       u(5,1,1)=CMPLX(fs,0._real_8,kind=real_8)! L = 4
       u(5,1,5)=CMPLX(fts,0._real_8,kind=real_8)
       u(5,1,9)=CMPLX(fs,0._real_8,kind=real_8)
       u(5,2,2)=CMPLX(0._real_8,-0.25_real_8,kind=real_8)
       u(5,2,4)=CMPLX(0._real_8,-s,kind=real_8)
       u(5,2,6)=CMPLX(0._real_8,-s,kind=real_8)
       u(5,2,8)=CMPLX(0._real_8,-0.25_real_8,kind=real_8)
       u(5,3,2)=CMPLX(-0.25_real_8,0._real_8,kind=real_8)
       u(5,3,4)=CMPLX( s,0._real_8,kind=real_8)
       u(5,3,6)=CMPLX(-s,0._real_8,kind=real_8)
       u(5,3,8)=CMPLX( 0.25_real_8,0._real_8,kind=real_8)
       u(5,4,3)=CMPLX(-y,0._real_8,kind=real_8)
       u(5,4,7)=CMPLX(-y,0._real_8,kind=real_8)
       u(5,5,1)=CMPLX(0._real_8, y,kind=real_8)
       u(5,5,9)=CMPLX(0._real_8,-y,kind=real_8)
       u(5,6,3)=CMPLX(0._real_8, y,kind=real_8)
       u(5,6,7)=CMPLX(0._real_8,-y,kind=real_8)
       u(5,7,2)=CMPLX(-s,0._real_8,kind=real_8)
       u(5,7,4)=CMPLX(-0.25_real_8,0._real_8,kind=real_8)
       u(5,7,6)=CMPLX( 0.25_real_8,0._real_8,kind=real_8)
       u(5,7,8)=CMPLX( s,0._real_8,kind=real_8)
       u(5,8,2)=CMPLX(0._real_8, s,kind=real_8)    ! MISSPRINT IN WEISSBLUTH (PAGE 129)
       u(5,8,4)=CMPLX(0._real_8,-0.25_real_8,kind=real_8)! MISSPRINT IN WEISSBLUTH (PAGE 129)
       u(5,8,6)=CMPLX(0._real_8,-0.25_real_8,kind=real_8)! MISSPRINT IN WEISSBLUTH (PAGE 129)
       u(5,8,8)=CMPLX(0._real_8, s,kind=real_8)    ! MISSPRINT IN WEISSBLUTH (PAGE 129)
       u(5,9,1)=CMPLX(-ss,0._real_8,kind=real_8)
       u(5,9,5)=CMPLX( ts,0._real_8,kind=real_8)
       u(5,9,9)=CMPLX(-ss,0._real_8,kind=real_8)
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL zeroing(aa)!,LiX*MiX*lix*mix*lx*mx)
    DO lp=1,ll
       DO l=1,LLi
          DO m=1,(2*l-1)
             DO k=1,LLi
                DO n=1,(2*k-1)
                   DO np=1,(2*k-1)
                      DO p=1,(2*l-1)
                         mp=(p-l)+(np-k)+lp! M' = P + N'
                         IF ((mp.GE.1).AND.(mp.LE.(2*lp-1))) THEN
                            aa(lp,mp,l,m,k,n)=aa(lp,mp,l,m,k,n)+&
                                 u(l,m,p)*u(k,n,np)*cc(l,p,k,np,lp)
                         ENDIF
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    DO ik=1,nlx
       DO il=1,nlx
          DO ii=1,mx
             lpl(il,ik,ii)=0
          ENDDO
          lpx(il,ik)=0
       ENDDO
    ENDDO
    DO lp=1,ll
       DO mp=1,(2*lp-1)
          DO l=1,LLi
             DO m=1,(2*l-1)
                DO k=1,LLi
                   DO n=1,(2*k-1)
                      sum=CMPLX(0._real_8,0._real_8,kind=real_8)
                      DO mpp=1,(2*lp-1)
                         sum = sum + CONJG(u(lp,mp,mpp))*aa(lp,mpp,l,m,k,n)
                      ENDDO
                      IF (ABS(sum).GT.0.001_real_8) THEN
                         il =(l-1)*(l-1)+m
                         ik =(k-1)*(k-1)+n
                         ilp=(lp-1)*(lp-1)+mp
                         IF (ABS(AIMAG(sum)).GT.0.001_real_8) THEN
                            IF (paral%io_parent)&
                                 WRITE(6,'(A,2F8.5,1X,3I3)')&
                                 ' ! !!  ERROR  !!!   AP(ILP,IL,IK) =  ',&
                                 sum,ilp,il,ik
                            CALL stopgm('AAINIT',' ',& 
                                 __LINE__,__FILE__)
                         ENDIF
                         ap(ilp,il,ik)=REAL(sum)
                         lpx(il,ik)=lpx(il,ik)+1
                         lpl(il,ik,lpx(il,ik))=ilp
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE aainit
  ! ==================================================================

END MODULE aainit_utils
