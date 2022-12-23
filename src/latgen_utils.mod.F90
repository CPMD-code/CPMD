MODULE latgen_utils
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: paral

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: latgen
  PUBLIC :: genlat
  PUBLIC :: omegagen

CONTAINS

  ! ==================================================================
  SUBROUTINE latgen(ibrav,celldm,a1,a2,a3,omega)
    ! ==--------------------------------------------------------------==
    ! ==  SETS UP THE CRYSTALLOGRAPHIC VECTORS A1,A2, AND A3.         ==
    ! ==  AND THE VOLUME OF THE CELL                                  ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ibrav
    REAL(real_8)                             :: celldm(6), a1(3), a2(3), &
                                                a3(3), omega

    REAL(real_8), PARAMETER :: sr3 = 1.73205080756887729352_real_8 

    INTEGER                                  :: ir
    REAL(real_8)                             :: cbya, singam, sinval, term, &
                                                term1, term2

! ==--------------------------------------------------------------==

    DO ir=1,3
       a1(ir)=0._real_8
       a2(ir)=0._real_8
       a3(ir)=0._real_8
    ENDDO
    IF (celldm(1).EQ.0)&
         CALL stopgm('LATGEN','THE LATTICE CONSTANT IS NULL! ',& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! ==  IBRAV =  1   Simple Cubic                                   ==
    ! ==           2   Face Centered Cubic                            ==
    ! ==           3   Body Centered Cubic                            ==
    ! ==           4   Hexagonal                                      ==
    ! ==           5   Rhombohedral or trigonal                       ==
    ! ==           6   Tetragonal                                     ==
    ! ==           7   Body Centred Tetragonal                        ==
    ! ==           8   Orthorhombic                                   ==
    ! ==          12   Monoclinic (1 angle different from 90)         ==
    ! ==          14   Triclinic (3 angles different from 90)         ==
    ! ==--------------------------------------------------------------==
    IF (ibrav.EQ.1) THEN
       ! ==------------------------------------------------------------==
       ! == Simple Cubic                                               ==
       ! ==------------------------------------------------------------==
       a1(1)=celldm(1)
       a2(2)=celldm(1)
       a3(3)=celldm(1)
    ELSEIF (ibrav.EQ.2) THEN
       ! ==------------------------------------------------------------==
       ! == Face Centered Cubic                                        ==
       ! ==------------------------------------------------------------==
       term=celldm(1)/2._real_8
       a1(1)=-term
       a1(3)= term
       a2(2)= term
       a2(3)= term
       a3(1)=-term
       a3(2)= term
    ELSEIF (ibrav.EQ.3) THEN
       ! ==------------------------------------------------------------==
       ! == Body Centered Cubic                                        ==
       ! ==------------------------------------------------------------==
       term=celldm(1)/2._real_8
       DO ir=1,3
          a1(ir)=term
          a2(ir)=term
          a3(ir)=term
       ENDDO
       a2(1)=-term
       a3(1)=-term
       a3(2)=-term
    ELSEIF (ibrav.EQ.4) THEN
       ! ==------------------------------------------------------------==
       ! == Hexagonal                                                  ==
       ! ==------------------------------------------------------------==
       IF (celldm(3).EQ.0._real_8)&
            CALL stopgm('LATGEN', 'C/A (THIRD ARGUMENT) IS NULL',& 
            __LINE__,__FILE__)
       cbya=celldm(3)
       a1(1)=celldm(1)
       a2(1)=-celldm(1)/2._real_8
       a2(2)=celldm(1)*sr3/2._real_8
       a3(3)=celldm(1)*cbya
    ELSEIF (ibrav.EQ.5) THEN
       ! ==------------------------------------------------------------==
       ! == Trigonal or rhombohedral                                   ==
       ! ==------------------------------------------------------------==
       IF (celldm(4).EQ.1._real_8)&
            CALL stopgm('LATGEN','THE ANGLE IS NULL! ',& 
            __LINE__,__FILE__)
       term1=SQRT(1._real_8+2._real_8*celldm(4))
       term2=SQRT(1._real_8-celldm(4))
       a2(2)=SQRT(2._real_8)*celldm(1)*term2/sr3
       a2(3)=celldm(1)*term1/sr3
       a1(1)=celldm(1)*term2/SQRT(2._real_8)
       a1(2)=-a1(1)/sr3
       a1(3)=a2(3)
       a3(1)=-a1(1)
       a3(2)=a1(2)
       a3(3)=a2(3)
    ELSEIF (ibrav.EQ.6) THEN
       ! ==------------------------------------------------------------==
       ! == Tetragonal                                                 ==
       ! ==------------------------------------------------------------==
       IF (celldm(3).EQ.0._real_8)&
            CALL stopgm('LATGEN', 'C/A (THIRD ARGUMENT) IS NULL',& 
            __LINE__,__FILE__)
       cbya=celldm(3)
       a1(1)=celldm(1)
       a2(2)=celldm(1)
       a3(3)=celldm(1)*cbya
    ELSEIF (ibrav.EQ.7) THEN
       ! ==------------------------------------------------------------==
       ! == Body Centred Tetragonal                                    ==
       ! ==------------------------------------------------------------==
       cbya=celldm(3)
       a2(1)=celldm(1)/2._real_8
       a2(2)=a2(1)
       a2(3)=cbya*celldm(1)/2._real_8
       a1(1)=a2(1)
       a1(2)=-a2(1)
       a1(3)=a2(3)
       a3(1)=-a2(1)
       a3(2)=-a2(1)
       a3(3)=a2(3)
    ELSEIF (ibrav.EQ.8) THEN
       ! ==------------------------------------------------------------==
       ! == Orthorhombic                                               ==
       ! ==------------------------------------------------------------==
       IF (celldm(2).EQ.0._real_8)&
            CALL stopgm('LATGEN', 'THE SECOND ARGUMENT IS NULL',& 
            __LINE__,__FILE__)
       IF (celldm(3).EQ.0._real_8)&
            CALL stopgm('LATGEN', 'THE THIRD ARGUMENT IS NULL',& 
            __LINE__,__FILE__)
       a1(1)=celldm(1)
       a2(2)=celldm(1)*celldm(2)
       a3(3)=celldm(1)*celldm(3)
    ELSEIF (ibrav.EQ.12) THEN
       ! ==------------------------------------------------------------==
       ! == Monoclinic (1 angle different from 90)                     ==
       ! ==------------------------------------------------------------==
       IF (celldm(2).EQ.0._real_8)&
            CALL stopgm('LATGEN', 'THE SECOND ARGUMENT IS NULL',& 
            __LINE__,__FILE__)
       IF (celldm(3).EQ.0._real_8)&
            CALL stopgm('LATGEN', 'THE THIRD ARGUMENT IS NULL',& 
            __LINE__,__FILE__)
       IF (celldm(4).EQ.1._real_8)&
            CALL stopgm('LATGEN','THE ANGLE IS NULL! ',& 
            __LINE__,__FILE__)
       sinval=SQRT(1._real_8-celldm(4)**2)
       a1(1)=celldm(1)
       a2(1)=celldm(1)*celldm(2)*celldm(4)
       a2(2)=celldm(1)*celldm(2)*sinval
       a3(3)=celldm(1)*celldm(3)
    ELSEIF (ibrav.EQ.14) THEN
       ! ==------------------------------------------------------------==
       ! == Triclinic (3 angles different from 90)                     ==
       ! ==------------------------------------------------------------==
       IF (celldm(2).EQ.0._real_8)&
            CALL stopgm('LATGEN', 'THE SECOND ARGUMENT IS NULL',& 
            __LINE__,__FILE__)
       IF (celldm(3).EQ.0._real_8)&
            CALL stopgm('LATGEN', 'THE THIRD ARGUMENT IS NULL',& 
            __LINE__,__FILE__)
       IF (celldm(4).EQ.1._real_8)&
            CALL stopgm('LATGEN','THE FIRST ANGLE IS NULL! ',& 
            __LINE__,__FILE__)
       IF (celldm(5).EQ.1._real_8)&
            CALL stopgm('LATGEN','THE SECOND ANGLE IS NULL! ',& 
            __LINE__,__FILE__)
       IF (celldm(6).EQ.1._real_8)&
            CALL stopgm('LATGEN','THE THIRD ANGLE IS NULL! ',& 
            __LINE__,__FILE__)
       ! Compatiblity between angle
       IF ( ACOS(celldm(4))+ACOS(celldm(5)).LT.ACOS(celldm(6)) )&
            CALL stopgm('LATGEN','ALPHA + BETA >= GAMMA',& 
            __LINE__,__FILE__)
       IF ( ACOS(celldm(5))+ACOS(celldm(6)).LT.ACOS(celldm(4)) )&
            CALL stopgm('LATGEN','BETA + GAMMA >= ALPHA',& 
            __LINE__,__FILE__)
       IF ( ACOS(celldm(6))+ACOS(celldm(4)).LT.ACOS(celldm(5)) )&
            CALL stopgm('LATGEN','GAMMA + ALPHA >= BETA',& 
            __LINE__,__FILE__)
       singam=SQRT(1._real_8-celldm(6)**2)
       a1(1)=celldm(1)
       a2(1)=celldm(1)*celldm(2)*celldm(6)
       a2(2)=celldm(1)*celldm(2)*singam
       a3(1)=celldm(1)*celldm(3)*celldm(5)
       a3(2)=celldm(1)*celldm(3)*(celldm(4)-celldm(5)*celldm(6))/singam
       term=SQRT((1._real_8+2._real_8*celldm(4)*celldm(5)*celldm(6)&
            -celldm(4)**2-celldm(5)**2&
            -celldm(6)**2)/(1._real_8-celldm(6)**2))
       a3(3)=celldm(1)*celldm(3)*term
    ELSEIF(ibrav.EQ.9 .OR.&
         ibrav.EQ.10.OR.&
         ibrav.EQ.11.OR.&
         ibrav.EQ.13) THEN
       ! ==------------------------------------------------------------==
       ! == Not defined IBRAV numbers (9, 10, 11, 13)                  ==
       ! ==------------------------------------------------------------==
       IF (paral%io_parent)&
            WRITE(6,'(A,I3,A)') ' BRAVAIS LATTICE',ibrav,&
            ' NOT PROGRAMMED. STOPPING'
       CALL stopgm('LATGEN',' ',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ==  Compute the volume of the Supercell (OMEGA)                 ==
    ! ==--------------------------------------------------------------==
    CALL omegagen(a1,a2,a3,omega)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE latgen
  ! ==================================================================
  SUBROUTINE genlat(a1,a2,a3,celldm,omega)
    ! ==--------------------------------------------------------------==
    ! == COMPUTE CELLDM(1:6) Dx Dy Dz alpha(yz) beta(zx) gamma(xy)    ==
    ! == FROM A1(3), A2(3), A3(3)                                     ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: a1(3), a2(3), a3(3), &
                                                celldm(6), omega

! ==--------------------------------------------------------------==
! CELLDM(1) = |A1|

    celldm(1) = SQRT(a1(1)*a1(1)+a1(2)*a1(2)+a1(3)*a1(3))
    IF (celldm(1).EQ.0._real_8) THEN
       IF (paral%io_parent)&
            WRITE(6,fmt=*) 'GENLAT! THE FIRST CELL VECTOR A1=',a1
       CALL stopgm('GENLAT','|A1| == 0.0',& 
            __LINE__,__FILE__)
    ENDIF
    ! CELLDM(2) = |A2|/|A1|
    celldm(2) = SQRT(a2(1)*a2(1)+a2(2)*a2(2)+a2(3)*a2(3))
    IF (celldm(2).EQ.0._real_8) THEN
       IF (paral%io_parent)&
            WRITE(6,fmt=*) 'GENLAT! THE SECOND CELL VECTOR A2=',a2
       CALL stopgm('GENLAT','|A2| == 0.0',& 
            __LINE__,__FILE__)
    ENDIF
    ! CELLDM(3) = |A3|/|A1|
    celldm(3) = SQRT(a3(1)*a3(1)+a3(2)*a3(2)+a3(3)*a3(3))
    IF (celldm(3).EQ.0._real_8) THEN
       IF (paral%io_parent)&
            WRITE(6,fmt=*) 'GENLAT! THE THIRD CELL VECTOR A3=',a3
       CALL stopgm('GENLAT','|A3| == 0.0',& 
            __LINE__,__FILE__)
    ENDIF
    ! CELLDM(4) = A2.A3/(|A2|.|A3|)
    celldm(4) =&
         (a2(1)*a3(1)+a2(2)*a3(2)+a2(3)*a3(3))/(celldm(2)*celldm(3))
    IF (celldm(4).EQ.1._real_8) THEN
       CALL stopgm('GENLAT','A2 AND A3 ARE IDENTICAL! ',& 
            __LINE__,__FILE__)
    ENDIF
    ! CELLDM(5) = A3.A1/(|A3|.|A1|)
    celldm(5) =&
         (a3(1)*a1(1)+a3(2)*a1(2)+a3(3)*a1(3))/(celldm(3)*celldm(1))
    IF (celldm(5).EQ.1._real_8) THEN
       CALL stopgm('GENLAT','A3 AND A1 ARE IDENTICAL! ',& 
            __LINE__,__FILE__)
    ENDIF
    ! CELLDM(6) = A1.A2/(|A1|.|A2|)
    celldm(6) =&
         (a1(1)*a2(1)+a1(2)*a2(2)+a1(3)*a2(3))/(celldm(1)*celldm(2))
    IF (celldm(6).EQ.1._real_8) THEN
       CALL stopgm('GENLAT','A1 AND A2 ARE IDENTICAL! ',& 
            __LINE__,__FILE__)
    ENDIF
    celldm(2) = celldm(2)/celldm(1)
    celldm(3) = celldm(3)/celldm(1)
    ! ==--------------------------------------------------------------==
    ! ==  Compute the volume of the Supercell (OMEGA)                 ==
    ! ==--------------------------------------------------------------==
    CALL omegagen(a1,a2,a3,omega)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE genlat
  ! ==================================================================
  SUBROUTINE omegagen(a1,a2,a3,omega)
    ! ==--------------------------------------------------------------==
    ! == COMPUTE THE VOLUME OF THE SUPERCELL (OMEGA)                  ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: a1(3), a2(3), a3(3), omega

    INTEGER                                  :: i, iperm, ivec, j, k, l
    REAL(real_8)                             :: s

! ==--------------------------------------------------------------==

    omega=0._real_8
    DO ivec = 1,2
       IF (ivec.EQ.1) THEN
          s=1._real_8
          i=1
          j=2
          k=3
       ELSE
          s=-1._real_8
          i=2
          j=1
          k=3
       ENDIF
       DO iperm=1,3
          omega=omega+s*a1(i)*a2(j)*a3(k)
          l=i
          i=j
          j=k
          k=l
       ENDDO
    ENDDO
    omega=ABS(omega)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE omegagen
  ! ==================================================================

END MODULE latgen_utils
