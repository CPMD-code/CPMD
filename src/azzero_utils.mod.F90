!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! THIS MODULE IS DEPRECATED, USE zeroing_utils INSTEAD !!
!! THIS MODULE IS DEPRECATED, USE zeroing_utils INSTEAD !!
!! THIS MODULE IS DEPRECATED, USE zeroing_utils INSTEAD !!
!! THIS MODULE IS DEPRECATED, USE zeroing_utils INSTEAD !!
!! THIS MODULE IS DEPRECATED, USE zeroing_utils INSTEAD !!
!! THIS MODULE IS DEPRECATED, USE zeroing_utils INSTEAD !!
!! THIS MODULE IS DEPRECATED, USE zeroing_utils INSTEAD !!
!! THIS MODULE IS DEPRECATED, USE zeroing_utils INSTEAD !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#if defined(__SR11000)
!option OPT(O(ss))
#endif


MODULE azzero_utils
  USE kinds,                           ONLY: int_8,&
                                             real_8

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: azzero_dont_use_me_anymore
  PUBLIC :: zazzero_dont_use_me_anymore
  PUBLIC :: iazzero_dont_use_me_anymore
  PUBLIC :: i8azzero_dont_use_me_anymore

CONTAINS

  ! ==================================================================
  ! AZZERO(A,N) -> A(I)=0 I=1..N for AC     
  ! ==================================================================
#if defined(__BG) && defined(__HAS_IBM_QPX_INTRINSIC)
  ! ==--------------------------------------------------------------==
  SUBROUTINE azzero_dont_use_me_anymore(x,n)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: n
    REAL(real_8), DIMENSION(n)               :: x

    vector(REAL(8)) :: v_zero
    INTEGER :: i,m,l
    INTEGER(int_8) :: addr
    INTEGER(int_8),PARAMETER :: v_size=32,v_nelem=4,d_size=8
    REAL(real_8) :: zero
    ! ==--------------------------------------------------------------==
    zero = 0.0_real_8
    v_zero = vec_lds(0,zero)

    addr = loc(x)
    l = MOD( ( v_size - MOD( addr, v_size ) ) / d_size, v_nelem )
    DO i = 1, MIN( l, n )
       x(i) = zero
    ENDDO
    m = n - MOD( n-l, INT(v_nelem, 4) )
    !$omp parallel do
    DO i = l+1, m, v_nelem
       CALL vec_st( v_zero, 0, x(i) )
    ENDDO
    !$omp end parallel do

    DO i = m+1, n
       x(i) = zero
    ENDDO
    ! ==--------------------------------------------------------------==
  END SUBROUTINE azzero_dont_use_me_anymore
  ! ==--------------------------------------------------------------==
#elif defined(__VECTOR)
  ! ==--------------------------------------------------------------==
  SUBROUTINE azzero_dont_use_me_anymore(a,n)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: n
    REAL(real_8)                             :: a(n)

    INTEGER                                  :: i

    !$omp parallel do private (I)
#ifdef __SR8000
    !poption parallel
#endif
#ifdef _vpp_
    !OCL NOALIAS
#endif
    DO i=1,n
       a(i)=0.0_real_8
    ENDDO
    ! ==--------------------------------------------------------------==
  END SUBROUTINE azzero_dont_use_me_anymore
  ! ==--------------------------------------------------------------==
#else
  ! ==--------------------------------------------------------------==
  SUBROUTINE azzero_dont_use_me_anymore(a,n)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: n
    REAL(real_8)                             :: a(n)

    INTEGER                                  :: i, ii

    ii=MOD(n,4)
    !$omp parallel private(I)
    !$omp do
    DO i=1,ii
       a(i)=0.0_real_8
    ENDDO
    !$omp end do
    !$omp do
    DO i=1+ii,n,4
       a(i)  =0.0_real_8
       a(i+1)=0.0_real_8
       a(i+2)=0.0_real_8
       a(i+3)=0.0_real_8
    ENDDO
    !$omp end do
    !$omp end parallel
    ! ==--------------------------------------------------------------==
  END SUBROUTINE azzero_dont_use_me_anymore
  ! ==--------------------------------------------------------------==
#endif
  ! ==================================================================
  ! ZAZZERO(A,N) -> A(I)=0 I=1..N for A complex*16
  ! ==================================================================
#if defined(__HP) 
  ! ==--------------------------------------------------------------==
  ! do nothing, special routine in sysdepend.c
  ! ==--------------------------------------------------------------==
#elif defined(__BG) && defined(__HAS_IBM_QPX_INTRINSIC)
  ! ==--------------------------------------------------------------==
  SUBROUTINE zazzero_dont_use_me_anymore(x,n)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: n
    COMPLEX(real_8), DIMENSION(n)            :: x

    vector(REAL(8)) :: v_zero
    INTEGER :: i,m,l
    INTEGER(int_8) :: addr
    INTEGER(int_8),PARAMETER :: v_size=32,v_nelem=2,c_size=16
    COMPLEX(real_8) :: zero
    ! ==--------------------------------------------------------------==
    zero = (0.0_real_8,0.0_real_8)
    v_zero = vec_lds(0,zero)

    addr = loc(x)
    l = MOD( ( v_size - MOD( addr, v_size ) ) / c_size, v_nelem )
    DO i = 1, MIN( l, n )
       x(i) = zero
    ENDDO
    m = n - MOD( n-l, INT(v_nelem, 4) )
    !$omp parallel do
    DO i = l+1, m, v_nelem
       CALL vec_st( v_zero, 0, x(i) )
    ENDDO
    !$omp end parallel do

    DO i = m+1, n
       x(i) = zero
    ENDDO
    ! ==--------------------------------------------------------------==
  END SUBROUTINE zazzero_dont_use_me_anymore
  ! ==--------------------------------------------------------------==
#elif defined(__VECTOR)
  SUBROUTINE zazzero_dont_use_me_anymore(a,n)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: n
    COMPLEX(real_8)                          :: a(n)

    COMPLEX(real_8), PARAMETER               :: zzero = (0._real_8,0._real_8)

    INTEGER                                  :: i

    !$omp parallel do private (I)
    DO i=1,n
       a(i)=zzero
    ENDDO
    ! ==--------------------------------------------------------------==
  END SUBROUTINE zazzero_dont_use_me_anymore
  ! ==--------------------------------------------------------------==
#else
  SUBROUTINE zazzero_dont_use_me_anymore(a,n)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: n
    COMPLEX(real_8)                          :: a(n)

    COMPLEX(real_8), PARAMETER               :: zzero = (0._real_8,0._real_8)

    INTEGER                                  :: i, ii

    ii=MOD(n,4)
    DO i=1,ii
       a(i)=zzero
    ENDDO
    !$omp parallel default(none) private(I) shared(ii,a,n)
    !$omp do
    DO i=1+ii,n,4
       a(i)  =zzero
       a(i+1)=zzero
       a(i+2)=zzero
       a(i+3)=zzero
    ENDDO
    !$omp end do
    !$omp end parallel
    ! ==--------------------------------------------------------------==
  END SUBROUTINE zazzero_dont_use_me_anymore
  ! ==--------------------------------------------------------------==
#endif
  ! ==================================================================
  ! IAZZERO(IA,N) -> IA(I)=0 I=1..N for IA default integer
  ! ==================================================================
  SUBROUTINE iazzero_dont_use_me_anymore(ia,n)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: n, ia(n)

    INTEGER                                  :: i

    !$omp parallel do private(I)
    DO i=1,n
       ia(i)=0
    ENDDO
    ! ==--------------------------------------------------------------==
  END SUBROUTINE iazzero_dont_use_me_anymore
  ! ==--------------------------------------------------------------==
  ! ==================================================================
  ! I8AZZERO(IA,N) -> IA(I)=0 I=1..N for IA integer*8
  ! ==================================================================
  SUBROUTINE i8azzero_dont_use_me_anymore(ia,n)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: n
    INTEGER(int_8)                           :: ia(n)

    INTEGER                                  :: i

    !$omp parallel do private(I)
    DO i=1,n
       ia(i)=0_int_8
    ENDDO
    ! ==--------------------------------------------------------------==
  END SUBROUTINE i8azzero_dont_use_me_anymore
  ! ==--------------------------------------------------------------==

END MODULE azzero_utils
