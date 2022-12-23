#include "cpmd_global.h"

MODULE cublas_utils

  USE cublas_interfaces,               ONLY: &
       CUBLAS_DIAG_NON_UNIT, CUBLAS_DIAG_UNIT, CUBLAS_FILL_MODE_LOWER, &
       CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N, CUBLAS_OP_T, CUBLAS_SIDE_LEFT, &
       CUBLAS_SIDE_RIGHT, CUBLAS_STATUS_SUCCESS, cublasCreate, cublasDasum, &
       cublasDcopy, cublasDestroy, cublasDgemm, cublasDger, cublasDscal, &
       cublasDsyrk, cublasDtrsm, cublasDzasum, cublasGetStream, &
       cublasGetVersion, cublasSetStream, cublasZcopy, cublasZdscal
  USE cublas_types,                    ONLY: cublas_handle_t
  USE cuda_types,                      ONLY: cuda_device_null,&
                                             cuda_memory_t,&
                                             cuda_stream_t
  USE cuda_utils,                      ONLY: cuda_check_device,&
                                             cuda_set_device
  USE error_handling,                  ONLY: stopgm

  USE, INTRINSIC :: ISO_C_BINDING,  ONLY: C_INT, C_SIZE_T, &
       C_CHAR, C_NULL_CHAR, C_NULL_PTR, C_DOUBLE
  USE kinds,                        ONLY: real_8, int_8, int_4
  USE sizeof_kinds,                 ONLY: sizeof_complex_8, sizeof_real_8, sizeof_int_4
  USE string_utils,                 ONLY: int2str


  IMPLICIT NONE

  PRIVATE

  PUBLIC :: cublas_create
  PUBLIC :: cublas_destroy
  PUBLIC :: cublas_get_version
  PUBLIC :: cublas_dcopy
  PUBLIC :: cublas_zcopy
  PUBLIC :: cublas_dscal
  PUBLIC :: cublas_zdscal
  PUBLIC :: cublas_dgemm
  PUBLIC :: cublas_dger
  PUBLIC :: cublas_dsyrk
  PUBLIC :: cublas_dasum
  PUBLIC :: cublas_dzasum
  PUBLIC :: cublas_dtrsm
  PUBLIC :: cublas_set_stream
  PUBLIC :: cublas_get_stream

CONTAINS


  SUBROUTINE cublas_create ( handle, device )
    TYPE(cublas_handle_t), INTENT(inout)     :: handle
    INTEGER, INTENT(IN)                      :: device

    CHARACTER(*), PARAMETER                  :: procedureN = 'cublas_create'

    INTEGER(KIND(CUBLAS_STATUS_SUCCESS))     :: c_status

#if defined(_HAS_CUDA)

    IF( handle%init ) CALL stopgm(procedureN,'handle already created',&
         __LINE__,__FILE__)
    CALL cuda_check_device ( device, procedureN )

    CALL cuda_set_device ( device )
    c_status = cublasCreate( handle%h )
    IF( c_status /= CUBLAS_STATUS_SUCCESS ) CALL stopgm(procedureN,"cublas error: "//TRIM(int2str( INT( c_status ) )),&
         __LINE__,__FILE__)

    handle%device = device
    handle%init = .TRUE.

#else

    CALL stopgm(procedureN,"no cublas available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE cublas_create


  SUBROUTINE cublas_destroy ( handle )
    TYPE(cublas_handle_t), INTENT(inout)     :: handle

    CHARACTER(*), PARAMETER                  :: procedureN = 'cublas_destroy'

    INTEGER(KIND(CUBLAS_STATUS_SUCCESS))     :: c_status

#if defined(_HAS_CUDA)

    IF( .NOT. handle%init ) CALL stopgm(procedureN,'handle not created',&
         __LINE__,__FILE__)

    CALL cuda_set_device ( handle%device )
    c_status = cublasDestroy( handle%h )
    IF( c_status /= CUBLAS_STATUS_SUCCESS ) CALL stopgm(procedureN,"cublas error: "//TRIM(int2str( INT( c_status ) )),&
         __LINE__,__FILE__)

    handle%device = cuda_device_null
    handle%init = .FALSE.

#else

    CALL stopgm(procedureN,"no cublas available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE cublas_destroy


  SUBROUTINE cublas_get_version ( handle, version )
    TYPE(cublas_handle_t), INTENT(in)        :: handle
    INTEGER, INTENT(out)                     :: version

    CHARACTER(*), PARAMETER :: procedureN = 'cublas_get_version'

    INTEGER(C_INT)                           :: c_version
    INTEGER(KIND(CUBLAS_STATUS_SUCCESS))     :: c_status

    version = 0

#if defined(_HAS_CUDA)

    IF( .NOT. handle%init ) CALL stopgm(procedureN,'handle not created',&
         __LINE__,__FILE__)

    c_status = cublasGetVersion( handle%h, c_version )
    IF( c_status /= CUBLAS_STATUS_SUCCESS ) CALL stopgm(procedureN,"cublas error: "//TRIM(int2str( INT( c_status ) )),&
         __LINE__,__FILE__)

    version = INT( c_version )

#else

    CALL stopgm(procedureN,"no cublas available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE cublas_get_version

  SUBROUTINE cublas_dcopy ( handle, n, x, incx, y, incy )
    TYPE(cublas_handle_t), INTENT(in)        :: handle
    INTEGER, INTENT(IN)                      :: n
    TYPE(cuda_memory_t)                      :: x
    INTEGER, INTENT(IN)                      :: incx
    TYPE(cuda_memory_t)                      :: y
    INTEGER, INTENT(IN)                      :: incy

    CHARACTER(*), PARAMETER                  :: procedureN = 'cublas_dcopy'

    INTEGER(C_INT)                           :: c_incx, c_incy, c_n
    INTEGER(KIND(CUBLAS_STATUS_SUCCESS))     :: c_status

#if defined(_HAS_CUDA)

    IF( .NOT. handle%init ) CALL stopgm(procedureN,'handle not created',&
         __LINE__,__FILE__)

    IF( .NOT. x%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT. y%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF(  handle%device /= x%device .OR. &
         handle%device /= y%device ) CALL stopgm(procedureN,'devices not consistent',&
         __LINE__,__FILE__)

    CALL cuda_check_device ( handle%device, procedureN )

    c_n = INT( n, C_INT )
    c_incx = INT( incx, C_INT )
    c_incy = INT( incy, C_INT )

    CALL cuda_set_device ( handle%device )
    c_status = cublasDcopy( handle%h, c_n, x%ptr, c_incx, y%ptr, c_incy )
    IF( c_status /= CUBLAS_STATUS_SUCCESS ) CALL stopgm(procedureN,"cublas error: "//TRIM(int2str( INT( c_status ) )),&
         __LINE__,__FILE__)

#else

    CALL stopgm(procedureN,"no cublas available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE cublas_dcopy


  SUBROUTINE cublas_zcopy ( handle, n, x, incx, y, incy )
    TYPE(cublas_handle_t), INTENT(in)        :: handle
    INTEGER, INTENT(IN)                      :: n
    TYPE(cuda_memory_t)                      :: x
    INTEGER, INTENT(IN)                      :: incx
    TYPE(cuda_memory_t)                      :: y
    INTEGER, INTENT(IN)                      :: incy

    CHARACTER(*), PARAMETER                  :: procedureN = 'cublas_zcopy'

    INTEGER(C_INT)                           :: c_incx, c_incy, c_n
    INTEGER(KIND(CUBLAS_STATUS_SUCCESS))     :: c_status

#if defined(_HAS_CUDA)

    IF( .NOT. handle%init ) CALL stopgm(procedureN,'handle not created',&
         __LINE__,__FILE__)

    IF( .NOT. x%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT. y%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF(  handle%device /= x%device .OR. &
         handle%device /= y%device ) CALL stopgm(procedureN,'devices not consistent',&
         __LINE__,__FILE__)

    CALL cuda_check_device ( handle%device, procedureN )

    c_n = INT( n, C_INT )
    c_incx = INT( incx, C_INT )
    c_incy = INT( incy, C_INT )

    CALL cuda_set_device ( handle%device )
    c_status = cublasZcopy( handle%h, c_n, x%ptr, c_incx, y%ptr, c_incy )
    IF( c_status /= CUBLAS_STATUS_SUCCESS ) CALL stopgm(procedureN,"cublas error: "//TRIM(int2str( INT( c_status ) )),&
         __LINE__,__FILE__)

#else

    CALL stopgm(procedureN,"no cublas available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE cublas_zcopy


  SUBROUTINE cublas_dscal ( handle, n, alpha, x, incx )
    TYPE(cublas_handle_t), INTENT(in)        :: handle
    INTEGER, INTENT(IN)                      :: n
    REAL(real_8), INTENT(IN)                 :: alpha
    TYPE(cuda_memory_t)                      :: x
    INTEGER, INTENT(IN)                      :: incx

    CHARACTER(*), PARAMETER                  :: procedureN = 'cublas_dscal'

    INTEGER(C_INT)                           :: c_incx, c_n
    INTEGER(KIND(CUBLAS_STATUS_SUCCESS))     :: c_status
    REAL(C_DOUBLE)                           :: c_alpha

#if defined(_HAS_CUDA)

    IF( .NOT. handle%init ) CALL stopgm(procedureN,'handle not created',&
         __LINE__,__FILE__)

    IF( .NOT. x%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF(  handle%device /= x%device ) CALL stopgm(procedureN,'devices not consistent',&
         __LINE__,__FILE__)

    CALL cuda_check_device ( handle%device, procedureN )

    c_n = INT( n, C_INT )
    c_incx = INT( incx, C_INT )
    c_alpha = REAL( alpha, C_DOUBLE )

    CALL cuda_set_device ( handle%device )
    c_status = cublasDscal( handle%h, c_n, c_alpha, x%ptr, c_incx )
    IF( c_status /= CUBLAS_STATUS_SUCCESS ) CALL stopgm(procedureN,"cublas error: "//TRIM(int2str( INT( c_status ) )),&
         __LINE__,__FILE__)

#else

    CALL stopgm(procedureN,"no cublas available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE cublas_dscal


  SUBROUTINE cublas_zdscal ( handle, n, alpha, x, incx )
    TYPE(cublas_handle_t), INTENT(in)        :: handle
    INTEGER, INTENT(IN)                      :: n
    REAL(real_8), INTENT(IN)                 :: alpha
    TYPE(cuda_memory_t)                      :: x
    INTEGER, INTENT(IN)                      :: incx

    CHARACTER(*), PARAMETER                  :: procedureN = 'cublas_zdscal'

    INTEGER(C_INT)                           :: c_incx, c_n
    INTEGER(KIND(CUBLAS_STATUS_SUCCESS))     :: c_status
    REAL(C_DOUBLE)                           :: c_alpha

#if defined(_HAS_CUDA)

    IF( .NOT. handle%init ) CALL stopgm(procedureN,'handle not created',&
         __LINE__,__FILE__)

    IF( .NOT. x%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF(  handle%device /= x%device ) CALL stopgm(procedureN,'devices not consistent',&
         __LINE__,__FILE__)

    CALL cuda_check_device ( handle%device, procedureN )

    c_n = INT( n, C_INT )
    c_incx = INT( incx, C_INT )
    c_alpha = REAL( alpha, C_DOUBLE )

    CALL cuda_set_device ( handle%device )
    c_status = cublasZdscal( handle%h, c_n, c_alpha, x%ptr, c_incx )
    IF( c_status /= CUBLAS_STATUS_SUCCESS ) CALL stopgm(procedureN,"cublas error: "//TRIM(int2str( INT( c_status ) )),&
         __LINE__,__FILE__)

#else

    CALL stopgm(procedureN,"no cublas available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE cublas_zdscal


  SUBROUTINE cublas_dgemm ( handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc )
    TYPE(cublas_handle_t), INTENT(in)        :: handle
    CHARACTER(len=1), INTENT(IN)             :: transa, transb
    INTEGER, INTENT(IN)                      :: m, n, k
    REAL(real_8), INTENT(IN)                 :: alpha
    TYPE(cuda_memory_t)                      :: A
    INTEGER, INTENT(IN)                      :: lda
    TYPE(cuda_memory_t)                      :: B
    INTEGER, INTENT(IN)                      :: ldb
    REAL(real_8), INTENT(IN)                 :: beta
    TYPE(cuda_memory_t)                      :: C
    INTEGER, INTENT(IN)                      :: ldc

    CHARACTER(*), PARAMETER                  :: procedureN = 'cublas_dgemm'

    INTEGER(C_INT)                           :: c_k, c_lda, c_ldb, c_ldc, &
                                                c_m, c_n
    INTEGER(KIND(CUBLAS_OP_N))               :: c_transa, c_transb
    INTEGER(KIND(CUBLAS_STATUS_SUCCESS))     :: c_status
    REAL(C_DOUBLE)                           :: c_alpha, c_beta

#if defined(_HAS_CUDA)

    IF( .NOT. handle%init ) CALL stopgm(procedureN,'handle not created',&
         __LINE__,__FILE__)

    IF( .NOT. A%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT. B%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT. C%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF(  handle%device /= A%device .OR. &
         handle%device /= B%device .OR. &
         handle%device /= C%device ) CALL stopgm(procedureN,'devices not consistent',&
         __LINE__,__FILE__)

    CALL cuda_check_device ( handle%device, procedureN )

    c_m = INT( m, C_INT )
    c_n = INT( n, C_INT )
    c_k = INT( k, C_INT )
    c_alpha = REAL( alpha, C_DOUBLE )
    c_beta = REAL( beta, C_DOUBLE )
    c_lda = INT( lda, C_INT )
    c_ldb = INT( ldb, C_INT )
    c_ldc = INT( ldc, C_INT )
    SELECT CASE( transa )
    CASE( 'N', 'n' ); c_transa = CUBLAS_OP_N
    CASE( 'T', 't' ); c_transa = CUBLAS_OP_T
    END SELECT
    SELECT CASE( transb )
    CASE( 'N', 'n' ); c_transb = CUBLAS_OP_N
    CASE( 'T', 't' ); c_transb = CUBLAS_OP_T
    END SELECT

    CALL cuda_set_device ( handle%device )
    c_status = cublasDgemm( handle%h, c_transa, c_transb, c_m, c_n, c_k, c_alpha, A%ptr, c_lda, B%ptr, c_ldb, c_beta, C%ptr, c_ldc )
    IF( c_status /= CUBLAS_STATUS_SUCCESS ) CALL stopgm(procedureN,"cublas error: "//TRIM(int2str( INT( c_status ) )),&
         __LINE__,__FILE__)

#else

    CALL stopgm(procedureN,"no cublas available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE cublas_dgemm


  SUBROUTINE cublas_dger ( handle, m, n, alpha, x, incx, y, incy, A, lda )
    TYPE(cublas_handle_t), INTENT(in)        :: handle
    INTEGER, INTENT(IN)                      :: m, n
    REAL(real_8), INTENT(IN)                 :: alpha
    TYPE(cuda_memory_t)                      :: x
    INTEGER, INTENT(IN)                      :: incx
    TYPE(cuda_memory_t)                      :: y
    INTEGER, INTENT(IN)                      :: incy
    TYPE(cuda_memory_t)                      :: A
    INTEGER, INTENT(IN)                      :: lda

    CHARACTER(*), PARAMETER                  :: procedureN = 'cublas_dger'

    INTEGER(C_INT)                           :: c_incx, c_incy, c_lda, c_m, &
                                                c_n
    INTEGER(KIND(CUBLAS_STATUS_SUCCESS))     :: c_status
    REAL(C_DOUBLE)                           :: c_alpha

#if defined(_HAS_CUDA)

    IF( .NOT. handle%init ) CALL stopgm(procedureN,'handle not created',&
         __LINE__,__FILE__)

    IF( .NOT. x%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT. y%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT. A%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF(  handle%device /= x%device .OR. &
         handle%device /= y%device .OR. &
         handle%device /= A%device ) CALL stopgm(procedureN,'devices not consistent',&
         __LINE__,__FILE__)

    CALL cuda_check_device ( handle%device, procedureN )

    c_m = INT( m, C_INT )
    c_n = INT( n, C_INT )
    c_alpha = REAL( alpha, C_DOUBLE )
    c_incx = INT( incx, C_INT )
    c_incy = INT( incy, C_INT )
    c_lda = INT( lda, C_INT )

    CALL cuda_set_device ( handle%device )
    c_status = cublasDger( handle%h, c_m, c_n, c_alpha, x%ptr, c_incx, y%ptr, c_incy, A%ptr, c_lda )
    IF( c_status /= CUBLAS_STATUS_SUCCESS ) CALL stopgm(procedureN,"cublas error: "//TRIM(int2str( INT( c_status ) )),&
         __LINE__,__FILE__)

#else

    CALL stopgm(procedureN,"no cublas available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE cublas_dger


  SUBROUTINE cublas_dsyrk ( handle, uplo, trans, n, k, alpha, A, lda, beta, C, ldc )
    TYPE(cublas_handle_t), INTENT(in)        :: handle
    CHARACTER(len=1), INTENT(IN)             :: uplo, trans
    INTEGER, INTENT(IN)                      :: n, k
    REAL(real_8), INTENT(IN)                 :: alpha
    TYPE(cuda_memory_t), INTENT(INOUT)       :: A
    INTEGER, INTENT(IN)                      :: lda
    REAL(real_8), INTENT(IN)                 :: beta
    TYPE(cuda_memory_t), INTENT(INOUT)       :: C
    INTEGER, INTENT(IN)                      :: ldc

    CHARACTER(*), PARAMETER                  :: procedureN = 'cublas_dsyrk'

    INTEGER(C_INT)                           :: c_k, c_lda, c_ldc, c_n
    INTEGER(KIND(CUBLAS_OP_N))               :: c_trans, c_uplo
    INTEGER(KIND(CUBLAS_STATUS_SUCCESS))     :: c_status
    REAL(C_DOUBLE)                           :: c_alpha, c_beta

#if defined(_HAS_CUDA)

    IF( .NOT. handle%init ) CALL stopgm(procedureN,'handle not created',&
         __LINE__,__FILE__)

    IF( .NOT. A%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT. C%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF(  handle%device /= A%device .OR. &
         handle%device /= C%device ) CALL stopgm(procedureN,'devices not consistent',&
         __LINE__,__FILE__)

    CALL cuda_check_device ( handle%device, procedureN )

    c_n = INT( n, C_INT )
    c_k = INT( k, C_INT )
    c_alpha = REAL( alpha, C_DOUBLE )
    c_beta = REAL( beta, C_DOUBLE )
    c_lda = INT( lda, C_INT )
    c_ldc = INT( ldc, C_INT )
    SELECT CASE( uplo )
    CASE( 'U', 'u' ); c_uplo = CUBLAS_FILL_MODE_UPPER
    CASE( 'L', 'l' ); c_uplo = CUBLAS_FILL_MODE_LOWER
    END SELECT
    SELECT CASE( trans )
    CASE( 'N', 'n' ); c_trans = CUBLAS_OP_N
    CASE( 'T', 't' ); c_trans = CUBLAS_OP_T
    END SELECT

    CALL cuda_set_device ( handle%device )
    c_status = cublasDsyrk( handle%h, c_uplo, c_trans, c_n, c_k, c_alpha, A%ptr, c_lda, c_beta, C%ptr, c_ldc )
    IF( c_status /= CUBLAS_STATUS_SUCCESS ) CALL stopgm(procedureN,"cublas error: "//TRIM(int2str( INT( c_status ) )),&
         __LINE__,__FILE__)

#else

    CALL stopgm(procedureN,"no cublas available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE cublas_dsyrk


  SUBROUTINE cublas_dasum ( handle, n, x, incx, RESULT )
    TYPE(cublas_handle_t), INTENT(in)        :: handle
    INTEGER, INTENT(IN)                      :: n
    TYPE(cuda_memory_t)                      :: x
    INTEGER, INTENT(IN)                      :: incx
    REAL(real_8), INTENT(OUT)                :: RESULT

    CHARACTER(*), PARAMETER                  :: procedureN = 'cublas_dasum'

    INTEGER(C_INT)                           :: c_incx, c_n
    INTEGER(KIND(CUBLAS_STATUS_SUCCESS))     :: c_status
    REAL(C_DOUBLE)                           :: c_result

    RESULT = 0.0_real_8

#if defined(_HAS_CUDA)

    IF( .NOT. handle%init ) CALL stopgm(procedureN,'handle not created',&
         __LINE__,__FILE__)

    IF( .NOT. x%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF(  handle%device /= x%device ) CALL stopgm(procedureN,'devices not consistent',&
         __LINE__,__FILE__)

    CALL cuda_check_device ( handle%device, procedureN )

    c_n = INT( n, C_INT )
    c_incx = INT( incx, C_INT )

    CALL cuda_set_device ( handle%device )
    c_status = cublasDasum( handle%h, c_n, x%ptr, c_incx, c_result )
    IF( c_status /= CUBLAS_STATUS_SUCCESS ) CALL stopgm(procedureN,"cublas error: "//TRIM(int2str( INT( c_status ) )),&
         __LINE__,__FILE__)

    RESULT = REAL( c_result, real_8 )

#else

    CALL stopgm(procedureN,"no cublas available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE cublas_dasum


  SUBROUTINE cublas_dzasum ( handle, n, x, incx, RESULT )
    TYPE(cublas_handle_t), INTENT(in)        :: handle
    INTEGER, INTENT(IN)                      :: n
    TYPE(cuda_memory_t)                      :: x
    INTEGER, INTENT(IN)                      :: incx
    REAL(real_8), INTENT(OUT)                :: RESULT

    CHARACTER(*), PARAMETER                  :: procedureN = 'cublas_dzasum'

    INTEGER(C_INT)                           :: c_incx, c_n
    INTEGER(KIND(CUBLAS_STATUS_SUCCESS))     :: c_status
    REAL(C_DOUBLE)                           :: c_result

    RESULT = 0.0_real_8

#if defined(_HAS_CUDA)

    IF( .NOT. handle%init ) CALL stopgm(procedureN,'handle not created',&
         __LINE__,__FILE__)

    IF( .NOT. x%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF(  handle%device /= x%device ) CALL stopgm(procedureN,'devices not consistent',&
         __LINE__,__FILE__)

    CALL cuda_check_device ( handle%device, procedureN )

    c_n = INT( n, C_INT )
    c_incx = INT( incx, C_INT )

    CALL cuda_set_device ( handle%device )
    c_status = cublasDzasum( handle%h, c_n, x%ptr, c_incx, c_result )
    IF( c_status /= CUBLAS_STATUS_SUCCESS ) CALL stopgm(procedureN,"cublas error: "//TRIM(int2str( INT( c_status ) )),&
         __LINE__,__FILE__)

    RESULT = REAL( c_result, real_8 )

#else

    CALL stopgm(procedureN,"no cublas available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE cublas_dzasum


  SUBROUTINE cublas_dtrsm ( handle, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb )
    TYPE(cublas_handle_t), INTENT(in)        :: handle
    CHARACTER(len=1), INTENT(IN)             :: side, uplo, trans, diag
    INTEGER, INTENT(IN)                      :: m, n
    REAL(real_8), INTENT(IN)                 :: alpha
    TYPE(cuda_memory_t)                      :: A
    INTEGER, INTENT(IN)                      :: lda
    TYPE(cuda_memory_t)                      :: B
    INTEGER, INTENT(IN)                      :: ldb

    CHARACTER(*), PARAMETER                  :: procedureN = 'cublas_dtrsm'

    INTEGER(C_INT)                           :: c_lda, c_ldb, c_m, c_n
    INTEGER(KIND(CUBLAS_DIAG_UNIT))          :: c_diag
    INTEGER(KIND(CUBLAS_FILL_MODE_UPPER))    :: c_uplo
    INTEGER(KIND(CUBLAS_OP_N))               :: c_trans
    INTEGER(KIND(CUBLAS_SIDE_LEFT))          :: c_side
    INTEGER(KIND(CUBLAS_STATUS_SUCCESS))     :: c_status
    REAL(C_DOUBLE)                           :: c_alpha

#if defined(_HAS_CUDA)

    IF( .NOT. handle%init ) CALL stopgm(procedureN,'handle not created',&
         __LINE__,__FILE__)

    IF( .NOT. A%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT. B%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF(  handle%device /= A%device .OR. &
         handle%device /= B%device ) CALL stopgm(procedureN,'devices not consistent',&
         __LINE__,__FILE__)

    CALL cuda_check_device ( handle%device, procedureN )

    SELECT CASE( side )
    CASE( 'L', 'l' ); c_side = CUBLAS_SIDE_LEFT
    CASE( 'R', 'r' ); c_side = CUBLAS_SIDE_RIGHT
    END SELECT
    SELECT CASE( uplo )
    CASE( 'U', 'u' ); c_uplo = CUBLAS_FILL_MODE_UPPER
    CASE( 'L', 'l' ); c_uplo = CUBLAS_FILL_MODE_LOWER
    END SELECT
    SELECT CASE( trans )
    CASE( 'N', 'n' ); c_trans = CUBLAS_OP_N
    CASE( 'T', 't' ); c_trans = CUBLAS_OP_T
    END SELECT
    SELECT CASE( diag )
    CASE( 'U', 'u' ); c_diag = CUBLAS_DIAG_UNIT
    CASE( 'N', 'n' ); c_diag = CUBLAS_DIAG_NON_UNIT
    END SELECT
    c_m = INT( m, C_INT )
    c_n = INT( n, C_INT )
    c_alpha = REAL( alpha, C_DOUBLE )
    c_lda = INT( lda, C_INT )
    c_ldb = INT( ldb, C_INT )

    CALL cuda_set_device ( handle%device )
    c_status = cublasDtrsm( handle%h, c_side, c_uplo, c_trans, c_diag, c_m, c_n, c_alpha, A%ptr, c_lda, B%ptr, c_ldb )
    IF( c_status /= CUBLAS_STATUS_SUCCESS ) CALL stopgm(procedureN,"cublas error: "//TRIM(int2str( INT( c_status ) )),&
         __LINE__,__FILE__)

#else

    CALL stopgm(procedureN,"no cublas available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE cublas_dtrsm


  SUBROUTINE cublas_set_stream ( handle, stream )
    TYPE(cublas_handle_t), INTENT(in)        :: handle
    TYPE(cuda_stream_t)                      :: stream

    CHARACTER(*), PARAMETER :: procedureN = 'cublas_set_stream'

    INTEGER(KIND(CUBLAS_STATUS_SUCCESS))     :: c_status

#if defined(_HAS_CUDA)

    IF( .NOT. handle%init ) CALL stopgm(procedureN,'handle not created',&
         __LINE__,__FILE__)

    IF( .NOT. stream%init ) CALL stopgm(procedureN,'stream not allocated',&
         __LINE__,__FILE__)

    IF( handle%device /= stream%device ) CALL stopgm(procedureN,'stream and handle dont share the same device',&
         __LINE__,__FILE__)

    CALL cuda_check_device ( handle%device, procedureN )

    c_status = cublasSetStream ( handle%h, stream%s )
    IF( c_status /= CUBLAS_STATUS_SUCCESS ) CALL stopgm(procedureN,"cublas error: "//TRIM(int2str( INT( c_status ) )),&
         __LINE__,__FILE__)

#else

    CALL stopgm(procedureN,"no cublas available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE cublas_set_stream


  SUBROUTINE cublas_get_stream ( handle, stream )
    TYPE(cublas_handle_t), INTENT(in)        :: handle
    TYPE(cuda_stream_t)                      :: stream

    CHARACTER(*), PARAMETER :: procedureN = 'cublas_get_stream'

    INTEGER(KIND(CUBLAS_STATUS_SUCCESS))     :: c_status

#if defined(_HAS_CUDA)

    IF( .NOT. handle%init ) CALL stopgm(procedureN,'handle not created',&
         __LINE__,__FILE__)

    IF( stream%init ) CALL stopgm(procedureN,'stream already allocated',&
         __LINE__,__FILE__)

    CALL cuda_check_device ( handle%device, procedureN )

    c_status = cublasGetStream ( handle%h, stream%s )
    IF( c_status /= CUBLAS_STATUS_SUCCESS ) CALL stopgm(procedureN,"cublas error: "//TRIM(int2str( INT( c_status ) )),&
         __LINE__,__FILE__)

    stream%device = handle%device
    stream%init = .TRUE.

#else

    CALL stopgm(procedureN,"no cublas available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE cublas_get_stream


END MODULE cublas_utils
