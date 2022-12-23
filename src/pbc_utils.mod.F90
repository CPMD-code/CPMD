MODULE pbc_utils
  USE bc,                              ONLY: bc_com
  USE clas,                            ONLY: clas8
  USE error_handling,                  ONLY: stopgm
  USE isos,                            ONLY: isos1
  USE kinds,                           ONLY: real_8
  USE metr,                            ONLY: metr_com
  USE mm_dimmod,                       ONLY: clsaabox,&
                                             mm_stat
  USE system,                          ONLY: cntl, parm

  IMPLICIT NONE


  PRIVATE

  PUBLIC :: dist_vec
  PUBLIC :: dist_pbc
  PUBLIC :: pbc
  PUBLIC :: pbc3
  PUBLIC :: pbc_clas

CONTAINS

  FUNCTION dist_vec(r1,r2,tpbc)
    ! ==--------------------------------------------------------------==
    ! == THIS ROUTINE COMPUTES THE VECTOR DIFFERENCE R1-R2            ==
    ! == WITH (TPBC=.TRUE.) OR WITHOUT (TPBC=.FALSE.) CONSIDERING PBC ==
    ! ==--------------------------------------------------------------==
    ! Input vectors
    REAL(real_8), intent(in)                 :: r1(3), r2(3)
    ! Flag to choose if PBC should be applied
    LOGICAL, intent(in)                      :: tpbc
    ! Resulting vector
    REAL(real_8)                             :: dist_vec(3),dist_vec_(3)

! ==--------------------------------------------------------------==

    dist_vec_(1) = r1(1) - r2(1)
    dist_vec_(2) = r1(2) - r2(2)
    dist_vec_(3) = r1(3) - r2(3)
    IF (tpbc) THEN
       CALL pbc(dist_vec_(1), dist_vec_(2), dist_vec_(3), &
            dist_vec(1), dist_vec(2), dist_vec(3), &
            1, parm%apbc, parm%ibrav)
    ELSE
       dist_vec(:) = dist_vec_(:)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION dist_vec

  FUNCTION dist_pbc(r1,r2,tpbc)
    ! ==--------------------------------------------------------------==
    ! == THIS ROUTINE COMPUTES THE DISTANCE |R1-R2| WITH (TPBC=.TRUE.)==
    ! == OR WITHOUT (TPBC=.FALSE.) CONSIDERING PBC                    ==
    ! ==--------------------------------------------------------------==
    ! Input vectors
    REAL(real_8), intent(in)                 :: r1(3), r2(3)
    ! Flag to choose if PBC should be applied
    LOGICAL, intent(in)                      :: tpbc
    ! Distance output
    REAL(real_8)                             :: dist_pbc

    REAL(real_8)                             :: dR(3)

! ==--------------------------------------------------------------==

    dR = dist_vec(r1, r2, tpbc)
    dist_pbc=SQRT(dot_product(dR, dR))
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION dist_pbc

  ! ==================================================================
  SUBROUTINE pbc(x1,y1,z1,x2,y2,z2,m,ap,ibrv)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: x1, y1, z1, x2, y2, z2
    INTEGER                                  :: m
    REAL(real_8)                             :: ap(4)
    INTEGER                                  :: ibrv

    REAL(real_8)                             :: a(4)

! ==--------------------------------------------------------------==
! AK: QM/MM bugfix. 2005/11/15
! if we are in MM dimensions, we need to apply PBC for the MM box.
! problem spotted by bernd ensing and marco de vivo.

    IF (mm_stat) THEN  ! true means QM dimensions. just copy
       a(1)=ap(1)
       a(2)=ap(2)
       a(3)=ap(3)
       a(4)=ap(4)
    ELSE
       IF (.NOT.cntl%tqmmm)CALL stopgm('PBC','BUG:mm_stat=.false. w/o QM/MM'&
            ,& 
            __LINE__,__FILE__)
       IF ((ibrv.EQ.1).OR.(ibrv.EQ.8)) THEN! MM supports only orthorhombic boxes.
          a(1)=clsaabox%box_au(1)/2.0_real_8
          a(2)=clsaabox%box_au(2)/2.0_real_8
          a(3)=clsaabox%box_au(3)/2.0_real_8
       ELSE
          CALL stopgm('PBC','QM/MM NEEDS ORTHORHOMBIC BOX',& 
               __LINE__,__FILE__)
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (isos1%tclust.AND.mm_stat) THEN
       IF (isos1%toned) THEN
          CALL pbc1(x1,y1,z1,x2,y2,z2,m,a,ibrv)
       ELSEIF (isos1%ttwod) THEN
          CALL pbc2(x1,y1,z1,x2,y2,z2,m,a,ibrv)
       ELSE
          x2=x1
          y2=y1
          z2=z1
       ENDIF
    ELSE
       CALL pbc3(x1,y1,z1,x2,y2,z2,m,a,ibrv)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE pbc
  ! ==================================================================
  SUBROUTINE pbc3(x1,y1,z1,x2,y2,z2,m,a,ibrav)
    ! ==--------------------------------------------------------------==
    ! == CALCULATES THE PERIODIC BOUNDARY CONDITIONS FOR SOME         ==
    ! == STRUCTURES:                                                  ==
    ! ==             SIMPLE CUBIC          (IBRAV=1)                  ==
    ! ==             FACE CENTERED CUBIC   (IBRAV=2)                  ==
    ! ==             ORTHOROMBIC           (IBRAV=8)                  ==
    ! ==             TRIGONAL              (IBRAV=12)                 ==
    ! ==                                                              ==
    ! == OTHER THAN THESE, WE USE A GENERAL ROUTINE                   ==
    ! ==                                                              ==
    ! == FOR FCC THE ROUTINE IS TAKEN FROM                            ==
    ! == CCP5 INF. QUAT. #10,SEPT(1983)                               ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: x1, y1, z1, x2, y2, z2
    INTEGER                                  :: m
    REAL(real_8)                             :: a(4)
    INTEGER                                  :: ibrav

    INTEGER                                  :: i
    REAL(real_8)                             :: aa, b(4), c(4), dx, r(3), &
                                                s(3), s2, ta, x0, x2t, x3, &
                                                xic, xt, y0, y2t, y3

#ifdef __SGI
    ! ..bug with the optimizer
    x2 = x1
    y2 = y1
    z2 = z1
#endif  
    ta=2._real_8*a(1)*m
    aa=a(1)*m
    go to (5,10,99,25,99,99,99,15,99,99,99,20,99,25) ibrav
    ! ==--------------------------------------------------------------==
    ! == Simple Cubic                                                 ==
    ! ==--------------------------------------------------------------==
5   CONTINUE
    x2=dmod(x1,ta)
    y2=dmod(y1,ta)
    z2=dmod(z1,ta)
    x2=x2-ta*dint(x2/aa)
    y2=y2-ta*dint(y2/aa)
    z2=z2-ta*dint(z2/aa)
    RETURN
    ! ==--------------------------------------------------------------==
    ! == Face Centered Cubic                                          ==
    ! ==--------------------------------------------------------------==
10  CONTINUE
    s2=SQRT(2.0_real_8)
    x2t=( x1+y1)/s2
    y2t=(-x1+y1)/s2
    x2=x2t
    y2=y2t
    x2=dmod(x2,ta)
    y2=dmod(y2,ta)
    z2=dmod(z1,ta*s2)
    x2=x2-ta*dint(x2/aa)
    y2=y2-ta*dint(y2/aa)
    z2=z2-ta*s2*dint(z2/(s2*aa))
    IF (ABS(x2)+ABS(y2)+s2*ABS(z2).LT.ta) GOTO 11
    x2=x2-dsign(aa,x2)
    y2=y2-dsign(aa,y2)
    z2=z2-dsign(s2*aa,z2)
11  CONTINUE
    xt=(x2-y2)/s2
    y2=(x2+y2)/s2
    x2=xt
    RETURN
    ! ==--------------------------------------------------------------==
    ! == Orthorhombic                                                 ==
    ! ==--------------------------------------------------------------==
15  CONTINUE
    DO i=1,3
       b(i)=2*a(i)*m
       c(i)=a(i)*m
    ENDDO
    x2=dmod(x1,b(1))
    y2=dmod(y1,b(2))
    z2=dmod(z1,b(3))
    x2=x2-b(1)*dint(x2/c(1))
    y2=y2-b(2)*dint(y2/c(2))
    z2=z2-b(3)*dint(z2/c(3))
    RETURN
    ! ==--------------------------------------------------------------==
    ! == Hexagonal and special case of triclinic: monocl., gamma<>90  ==
    ! ==--------------------------------------------------------------==
25  CONTINUE
    IF (a(2).EQ.0.0_real_8) GOTO 99
    DO i=1,3
       c(i)=a(i)*m
       b(i)=2*c(i)
    ENDDO
    y0=y1
    y2=dmod(y1,b(2))
    z2=dmod(z1,b(3))
    y2=y2-b(2)*dint(y2/c(2))
    z2=z2-b(3)*dint(z2/c(3))
    dx=(y2-y0)*a(4)           ! A(4)=SQRT(3)/3
    x2=dmod(x1-dx,b(1))
    x2=x2-b(1)*dint(x2/c(1))
    RETURN
    ! ==--------------------------------------------------------------==
    ! == Trigonal                                                     ==
    ! ==--------------------------------------------------------------==
20  CONTINUE
    DO i=1,4
       b(i)=2*a(i)*m
       c(i)=a(i)*m
    END DO
    z2=dmod(z1,b(3))
    z2=z2-b(3)*dint(z2/c(3))
    CALL rot(x1,y1,x0,y0,bc_com%v1crys,1)
    x0=dmod(x0,b(1))
    x0=x0-b(1)*dint(x0/c(1))
    CALL rot(x0,y0,x3,y3,bc_com%v1crys,-1)
    CALL rot(x3,y3,x0,y0,bc_com%v2crys,1)
    x0=dmod(x0,b(2))
    x0=x0-b(2)*dint(x0/c(2))
    CALL rot(x0,y0,x3,y3,bc_com%v2crys,-1)
    CALL rot(x3,y3,x0,y0,bc_com%v3crys,1)
    x0=dmod(x0,b(4))
    x0=x0-b(4)*dint(x0/c(4))
    CALL rot(x0,y0,x3,y3,bc_com%v3crys,-1)
    CALL rot(x3,y3,x0,y0,bc_com%v1crys,1)
    x0=dmod(x0,b(1))
    x0=x0-b(1)*dint(x0/c(1))
    CALL rot(x0,y0,x2,y2,bc_com%v1crys,-1)
    RETURN
    ! ==--------------------------------------------------------------==
    ! == The rest : use a general procedure                           ==
    ! ==--------------------------------------------------------------==
99  CONTINUE
    r(1)=x1
    r(2)=y1
    r(3)=z1
    CALL dgemv('T',3,3,1.0_real_8,metr_com%htm1,3,r,1,0.0_real_8,s,1)
    xic=REAL(m,kind=real_8)
    s(1)=s(1)-NINT(s(1)/xic)*xic
    s(2)=s(2)-NINT(s(2)/xic)*xic
    s(3)=s(3)-NINT(s(3)/xic)*xic
    CALL dgemv('T',3,3,1.0_real_8,metr_com%ht,3,s,1,0.0_real_8,r,1)
    x2=r(1)
    y2=r(2)
    z2=r(3)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE pbc3
  ! ==================================================================
  SUBROUTINE rot(x1,y1,x0,y0,vv,isign)
    ! ==--------------------------------------------------------------==
    ! == THIS ROUTINE TRANSFORMS THE POINT (X1,Y1) INTO THE POINT     ==
    ! == (X0,Y0) CORRESPONDING TO A ROTATION OF AN ANGLE WHICH COS    ==
    ! == AND SIN ARE GIVEN BY VV(3) AND VV(4) RESPECTIVELY. THE       ==
    ! == SENSE OF THE ROTATION IS GIVEN BY ISIGN (+ FOR CLOCKWISE,    ==
    ! == AND - FOR ANTICLOCKWISE)                                     ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: x1, y1, x0, y0, vv(4)
    INTEGER                                  :: isign

    x0=x1*vv(3)+y1*vv(4)*REAL(isign,kind=real_8)
    y0=-x1*vv(4)*REAL(isign,kind=real_8)+y1*vv(3)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rot
  ! ==================================================================
  SUBROUTINE pbc1(x1,y1,z1,x2,y2,z2,m,a,ibrav)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: x1, y1, z1, x2, y2, z2
    INTEGER                                  :: m
    REAL(real_8)                             :: a(3)
    INTEGER                                  :: ibrav

    REAL(real_8)                             :: aa, ta

#ifdef __SGI
    x2 = x1
    y2 = y1
    z2 = z1
#endif
    ta=2._real_8*a(1)*m
    aa=a(1)*m
    x2=dmod(x1,ta)
    x2=x2-ta*dint(x2/aa)
    y2=y1
    z2=z1
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE pbc1
  ! ==================================================================
  SUBROUTINE pbc2(x1,y1,z1,x2,y2,z2,m,a,ibrav)
    REAL(real_8)                             :: x1, y1, z1, x2, y2, z2
    INTEGER                                  :: m
    REAL(real_8)                             :: a(4)
    INTEGER                                  :: ibrav

    INTEGER                                  :: i
    REAL(real_8)                             :: aa, b(3), c(3), dx, ta, y0

#ifdef __SGI
    x2 = x1
    y2 = y1
    z2 = z1
#endif
    ta=2._real_8*a(1)*m
    aa=a(1)*m
    IF (ibrav.EQ.1) THEN
       x2=dmod(x1,ta)
       y2=dmod(y1,ta)
       x2=x2-ta*dint(x2/aa)
       y2=y2-ta*dint(y2/aa)
       z2=z1
    ELSE IF (ibrav.EQ.4) THEN
       DO i=1,2
          c(i)=a(i)*m
          b(i)=2*c(i)
       ENDDO
       y0=y1
       y2=dmod(y1,b(2))
       y2=y2-b(2)*dint(y2/c(2))
       dx=(y2-y0)*a(4)       ! A(4)=SQRT(3)/3
       x2=dmod(x1-dx,b(1))
       x2=x2-b(1)*dint(x2/c(1))
       z2=z1
    ELSE
       DO i=1,3
          b(i)=2*a(i)*m
          c(i)=a(i)*m
       ENDDO
       x2=dmod(x1,b(1))
       y2=dmod(y1,b(2))
       x2=x2-b(1)*dint(x2/c(1))
       y2=y2-b(2)*dint(y2/c(2))
       z2=z1
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE pbc2
  ! ==================================================================
  SUBROUTINE pbc_clas(x1,y1,z1,x2,y2,z2)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: x1, y1, z1, x2, y2, z2

    REAL(real_8)                             :: r(3), s(3)

    r(1)=x1
    r(2)=y1
    r(3)=z1
    CALL dgemv('T',3,3,1.0_real_8,clas8%clhtm1,3,r,1,0.0_real_8,s,1)
    s(1)=s(1)-NINT(s(1))
    s(2)=s(2)-NINT(s(2))
    s(3)=s(3)-NINT(s(3))
    CALL dgemv('T',3,3,1.0_real_8,clas8%clht,3,s,1,0.0_real_8,r,1)
    x2=r(1)
    y2=r(2)
    z2=r(3)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE pbc_clas
  ! ==================================================================

END MODULE pbc_utils
