! ==================================================================
SUBROUTINE compress(n,a,ia,scale,icomp,nx)
  ! ==--------------------------------------------------------------==
  ! == COMPRESS THE DATA FROM A TO IA                               ==
  ! ==--------------------------------------------------------------==
  ! == INPUT:                                                       ==
  ! ==   ICOMP      compression type                                ==
  ! ==              64 bits -> 64/ICOMP bits                        ==
  ! ==   N          Size of A array                                 ==
  ! ==   A(1:N)     Array to compress                               ==
  ! ==   SCALE      Scale of data                                   ==
  ! == OUTPUT:                                                      ==
  ! ==   NX         Size of given array IA (in integer*8)           ==
  ! ==   IA(1:NX)   Compressed array          (integer*8)           ==
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE reshaper, ONLY: reshape_inplace, type_cast
  ! corresponding reshaper procedures
#if defined(__NEC) || defined(__SUN) || defined( __Linux) || defined(__SR8000) || defined( __WINNT)
  IMPLICIT NONE
  ! Arguments
  INTEGER :: n,icomp,nx
  REAL(real_8) :: scale
  INTEGER, TARGET :: ia(*) ! IA(2*N)
  REAL(real_8) :: a(*) ! A(N)
  ! Variables
  REAL(real_8), POINTER :: a0(:)
  INTEGER :: nb
  REAL(real_8) :: fx
  ! ==--------------------------------------------------------------==
  IF (icomp.EQ.1) THEN
     CALL reshape_inplace(ia, (/n/), a0)
     CALL dcopy(n,a(1),1,a0(1),1)
     nx=n
  ELSE
     nb=64/icomp
     ! Overflow: Integer 2._real_8**(NB-1)-1
     fx=2._real_8**(nb-1)-1._real_8
     fx=fx/scale
     CALL spack(ia,a,n,nx,nb,fx)
  ENDIF
#else
  IMPLICIT NONE
  INTEGER :: n,icomp,nx
  REAL(real_8), TARGET :: a(*) ! A(N)
  REAL(real_8) :: scale
  INTEGER(int_8) :: ia(*) ! IA(N)

  ! Variables
  INTEGER(int_1), POINTER :: ia1(:)
  INTEGER(int_2), POINTER :: ia2(:)
  INTEGER(int_4), POINTER :: ia4(:)
  REAL(real_8), POINTER :: a0(:)
  INTEGER :: nb,i
  REAL(real_8) :: fx
  ! ==--------------------------------------------------------------==
  IF (icomp.EQ.1) THEN
     CALL dcopy(n,a(1),1,ia(1),1)
     nx=n
  ELSE
     nb=64/icomp
     ! Overflow in integer*4 max=2**31-1
     ! NN=2**(NB-1)
     fx=2._real_8**(nb-1)-1._real_8
     fx=fx/scale
     nx=n/icomp+1
     IF (icomp.EQ.2) THEN
        CALL type_cast(ia, n, ia4)
        DO i=1,n
           ia4(i)=NINT(a(i)*fx)
        ENDDO
     ELSEIF (icomp.EQ.4) THEN
        CALL type_cast(ia, n, ia2)
        DO i=1,n
           ia2(i)=NINT(a(i)*fx)
        ENDDO
     ELSEIF (icomp.EQ.8) THEN
        CALL type_cast(ia, n, ia1)
        DO i=1,n
           ia1(i)=NINT(a(i)*fx)
        ENDDO
     ENDIF
  ENDIF
#endif
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE compress
! ==================================================================
SUBROUTINE decompr(na,a,ia,scale,icomp,nx)
  ! ==--------------------------------------------------------------==
  ! == DECOMPRESS THE DATA FROM IA TO A                             ==
  ! ==--------------------------------------------------------------==
  ! == INPUT:                                                       ==
  ! ==   ICOMP      Compression type                                ==
  ! ==              64 bits -> 64/ICOMP bits                        ==
  ! ==   NA         Maximum size of A                               ==
  ! ==   NX         Size of IA array (in integer*8)                 ==
  ! ==   IA(1:NX)   Array to decompress (integer*8)                 ==
  ! ==   SCALE      Scale for the decompressed data                 ==
  ! == OUTPUT:                                                      ==
  ! ==   A(1:N)     Decompressed array                              ==
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE reshaper, ONLY: type_cast
#if defined(__NEC) || defined(__SUN) || defined(__Linux) || defined(__SR8000)
  USE string_utils, ONLY: int2str
  USE parac, ONLY : paral,parai
  IMPLICIT NONE
  ! Arguments
  INTEGER :: na,icomp,nx
  INTEGER :: ia(*) ! IA(2*NX)
  REAL(real_8) :: scale
  REAL(real_8) :: a(*) ! A(NA)
  ! Variables
  REAL(real_8), POINTER :: a0(:)
  INTEGER :: n,nb
  REAL(real_8) :: fx
  ! ==--------------------------------------------------------------==
  IF (icomp.EQ.1) THEN
     n=nx
  ELSE
     n=(nx-1)*icomp
  ENDIF
  IF (n.GT.na) THEN
     CALL stopgm('DECOMPR',&
          'DECOMPR! THE SIZE OF THE DECOMPRESSED ARRAY '//&
          TRIM(int2str(n))//&
          ' ARE BIGGER THAN THE SIZE OF THE GIVEN ARRAY '//&
          TRIM(int2str(na)),& 
          __LINE__,__FILE__)
  ENDIF
  IF (icomp.EQ.1) THEN
     n=nx
     CALL dcopy(n,ia,1,a,1)
  ELSE
     nb=64/icomp
     ! Overflowfor integer: max=2._real_8**(NB-1)-1._real_8
     fx=2._real_8**(nb-1)-1._real_8
     fx=scale/fx
     n=(nx-1)*icomp
     CALL sunpack(ia,a,n,nx,nb,fx)
  ENDIF
#else
  USE string_utils, ONLY: int2str
  USE parac, ONLY : paral,parai
  IMPLICIT NONE
  INTEGER :: na,icomp,nx
  REAL(real_8) :: scale
  REAL(real_8) :: a(*) ! A(NA)
  INTEGER(int_8) :: ia(*) ! IA(NX)
  ! Variables

  INTEGER(int_1), POINTER :: ia1(:)
  INTEGER(int_2), POINTER :: ia2(:)
  INTEGER(int_4), POINTER :: ia4(:)
  REAL(real_8), POINTER :: a0(:)

  INTEGER :: nb,n,i
  REAL(real_8) :: fx
  ! ==--------------------------------------------------------------==
  IF (icomp.EQ.1) THEN
     n=nx
  ELSE
     n=(nx-1)*icomp
  ENDIF
  IF (n.GT.na) THEN
     CALL stopgm('DECOMPR',&
          'DECOMPR! THE SIZE OF THE DECOMPRESSED ARRAY '//&
          TRIM(int2str(n))//&
          ' ARE BIGGER THAN THE SIZE OF THE GIVEN ARRAY '//&
          TRIM(int2str(na)),& 
          __LINE__,__FILE__)
  ENDIF
  IF (icomp.EQ.1) THEN
     n=nx
     CALL dcopy(n,ia,1,a,1)
  ELSE
     nb=64/icomp
     ! Overflow in integer*4 max=2**31-1
     ! NN=2**(NB-1)
     fx=2._real_8**(nb-1)-1._real_8
     fx=scale/fx
     n=(nx-1)*icomp
     IF (icomp.EQ.2) THEN
        CALL type_cast(ia, n, ia4)
        DO i=1,n
           a(i)=ia4(i)*fx
        ENDDO
     ELSEIF (icomp.EQ.4) THEN
        CALL type_cast(ia, n, ia2)
        DO i=1,n
           a(i)=ia2(i)*fx
        ENDDO
     ELSEIF (icomp.EQ.8) THEN
        CALL type_cast(ia, n, ia1)
        DO i=1,n
           a(i)=ia1(i)*fx
        ENDDO
     ENDIF
  ENDIF
#endif
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE decompr
! ==================================================================
#if defined(__NEC)
SUBROUTINE spack(ia,a,n,nx,ibyte,fx)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  IMPLICIT NONE
  ! Arguments
  INTEGER :: ia(*),n,nx,ibyte
  REAL(real_8) :: a(*),fx
  ! Variables
  INTEGER :: nn,m,i,ii,&
       ib1,ib2,ib3,ib4,ib5,ib6,ib7,ib8
  ! ==--------------------------------------------------------------==
  nn=2**(ibyte-1)
  IF (ibyte.EQ.32) THEN
     m=n+MOD(n,2_8)
     nx=m/2
     DO i=1,m,2
        ii=(i+1)/2
        ib1=NINT(a(i)*fx)+nn
        ib2=NINT(a(i+1)*fx)+nn
        ia(ii)=ISHFT(ib1,32)+ib2
     ENDDO
  ELSEIF (ibyte.EQ.16) THEN
     m=n+MOD(n,4_8)
     nx=m/4
     DO i=1,m,4
        ii=(i+3)/4
        ib1=NINT(a(i)*fx)+nn
        ib2=NINT(a(i+1)*fx)+nn
        ib3=NINT(a(i+2)*fx)+nn
        ib4=NINT(a(i+3)*fx)+nn
        ia(ii)=ISHFT(ib1,48)+ISHFT(ib2,32)+ISHFT(ib3,16)+ib4
     ENDDO
  ELSEIF (ibyte.EQ.8) THEN
     m=n+MOD(n,8_8)
     nx=m/8
     DO i=1,m,8
        ii=(i+7)/8
        ib1=NINT(a(i)*fx)+nn
        ib2=NINT(a(i+1)*fx)+nn
        ib3=NINT(a(i+2)*fx)+nn
        ib4=NINT(a(i+3)*fx)+nn
        ib5=NINT(a(i+4)*fx)+nn
        ib6=NINT(a(i+5)*fx)+nn
        ib7=NINT(a(i+6)*fx)+nn
        ib8=NINT(a(i+7)*fx)+nn
        ia(ii)=ISHFT(ib1,56)+ISHFT(ib2,48)+ISHFT(ib3,40)+&
             ISHFT(ib4,32)+ISHFT(ib5,24)+ISHFT(ib6,16)+&
             ISHFT(ib7,8)+ib8
     ENDDO
  ELSE
     CALL stopgm('SPACK','IBYTE NOT IMPLEMENTED',& 
          __LINE__,__FILE__)
  ENDIF
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE spack
! ==================================================================
SUBROUTINE sunpack(ia,a,n,nx,ibyte,fx)
  ! ==--------------------------------------------------------------==
  ! == USED BY DECOMPR TO DECOMPRESS DATA FROM IA TO A              ==
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  IMPLICIT NONE
  ! Arguments
  INTEGER(int_8) :: ia(*),n,nx,ibyte
  ! mb      INTEGER IA(*),N,NX,IBYTE
  REAL(real_8) :: a(*),fx
  ! Variables
  INTEGER :: nn,i,ii,&
       ib1,ib2,ib3,ib4,ib5,ib6,ib7,ib8
  ! ==--------------------------------------------------------------==
  nn=2**(ibyte-1)
  IF (ibyte.EQ.32) THEN
     DO i=1,n,2
        ii=(i+1)/2
        ib1=ISHFT(ia(ii),-32_8)
        ib2=ISHFT(ISHFT(ia(ii),32_8),-32_8)
        a(i)=(ib1-nn)*fx
        a(i+1)=(ib2-nn)*fx
     ENDDO
  ELSEIF (ibyte.EQ.16) THEN
     DO i=1,n,4
        ii=(i+3)/4
        ib1=ISHFT(ia(ii),-48_8)
        ib2=ISHFT(ISHFT(ia(ii),16_8),-48_8)
        ib3=ISHFT(ISHFT(ia(ii),32_8),-48_8)
        ib4=ISHFT(ISHFT(ia(ii),48_8),-48_8)
        a(i)=(ib1-nn)*fx
        a(i+1)=(ib2-nn)*fx
        a(i+2)=(ib3-nn)*fx
        a(i+3)=(ib4-nn)*fx
     ENDDO
  ELSEIF (ibyte.EQ.8) THEN
     DO i=1,n,8
        ii=(i+7)/8
        ib1=ISHFT(ia(ii),-56_8)
        ib2=ISHFT(ISHFT(ia(ii), 8_8),-56_8)
        ib3=ISHFT(ISHFT(ia(ii),16_8),-56_8)
        ib4=ISHFT(ISHFT(ia(ii),24_8),-56_8)
        ib5=ISHFT(ISHFT(ia(ii),32_8),-56_8)
        ib6=ISHFT(ISHFT(ia(ii),40_8),-56_8)
        ib7=ISHFT(ISHFT(ia(ii),48_8),-56_8)
        ib8=ISHFT(ISHFT(ia(ii),56_8),-56_8)
        a(i)=(ib1-nn)*fx
        a(i+1)=(ib2-nn)*fx
        a(i+2)=(ib3-nn)*fx
        a(i+3)=(ib4-nn)*fx
        a(i+4)=(ib5-nn)*fx
        a(i+5)=(ib6-nn)*fx
        a(i+6)=(ib7-nn)*fx
        a(i+7)=(ib8-nn)*fx
     ENDDO
  ELSE
     CALL stopgm('SUNPACK','IBYTE NOT IMPLEMENTED',& 
          __LINE__,__FILE__)
  ENDIF
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE sunpack
! ==================================================================
#elif defined(__SUN) || defined(__Linux) || defined (__SR8000) || defined(__WINNT)
SUBROUTINE spack(ia,a,n,nx,ibyte,fx)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  IMPLICIT NONE
  ! Arguments
  INTEGER :: ia(*),n,nx,ibyte
  REAL(real_8) :: a(*),fx
  ! Variables
  INTEGER :: ib1,ib2,ib3,ib4,nn,m,i,ii
  ! ==--------------------------------------------------------------==
  nn=2**(ibyte-1)
  IF (ibyte.EQ.32) THEN
     m=n+MOD(n,2)
     nx=m/2
     DO i=1,m
        ii=i
        ib1=NINT(a(i)*fx)+nn
        ia(ii)=ib1
     ENDDO
  ELSEIF (ibyte.EQ.16) THEN
     m=n+MOD(n,4)
     nx=m/4
     DO i=1,m,2
        ii=(i+1)/2
        ib1=NINT(a(i)*fx)+nn
        ib2=NINT(a(i+1)*fx)+nn
        ia(ii)=ISHFT(ib1,16)+ib2
     ENDDO
  ELSEIF (ibyte.EQ.8) THEN
     m=n+MOD(n,8)
     nx=m/8
     DO i=1,m,4
        ii=(i+3)/4
        ib1=NINT(a(i)*fx)+nn
        ib2=NINT(a(i+1)*fx)+nn
        ib3=NINT(a(i+2)*fx)+nn
        ib4=NINT(a(i+3)*fx)+nn
        ia(ii)=ISHFT(ib1,24)+ISHFT(ib2,16)+ISHFT(ib3,8)+ib4
     ENDDO
  ELSE
     CALL stopgm('SPACK','IBYTE NOT IMPLEMENTED',& 
          __LINE__,__FILE__)
  ENDIF
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE spack
! ==================================================================
SUBROUTINE sunpack(ia,a,n,nx,ibyte,fx)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  IMPLICIT NONE
  ! Arguments
  INTEGER :: ia(*),n,nx,ibyte
  ! Variables
  REAL(real_8) :: a(*),fx
  INTEGER :: ib1,ib2,ib3,ib4,nn,i,ii
  ! ==--------------------------------------------------------------==
  nn=2**(ibyte-1)
  IF (ibyte.EQ.32) THEN
     DO i=1,n
        ii=i
        ib1=ia(ii)
        a(i)=(ib1-nn)*fx
     ENDDO
  ELSEIF (ibyte.EQ.16) THEN
     DO i=1,n,2
        ii=(i+1)/2
        ib1=ISHFT(ia(ii),-16)
        ib2=ISHFT(ISHFT(ia(ii),16),-16)
        a(i)=(ib1-nn)*fx
        a(i+1)=(ib2-nn)*fx
     ENDDO
  ELSEIF (ibyte.EQ.8) THEN
     DO i=1,n,4
        ii=(i+3)/4
        ib1=ISHFT(ia(ii),-24)
        ib2=ISHFT(ISHFT(ia(ii), 8),-24)
        ib3=ISHFT(ISHFT(ia(ii),16),-24)
        ib4=ISHFT(ISHFT(ia(ii),24),-24)
        a(i)=(ib1-nn)*fx
        a(i+1)=(ib2-nn)*fx
        a(i+2)=(ib3-nn)*fx
        a(i+3)=(ib4-nn)*fx
     ENDDO
  ELSE
     CALL stopgm('SUNPACK','IBYTE NOT IMPLEMENTED',& 
          __LINE__,__FILE__)
  ENDIF
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE sunpack
! ==================================================================
#endif
