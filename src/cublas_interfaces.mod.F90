MODULE cublas_interfaces

  USE cuda_interfaces,                 ONLY: cudastream_t

  USE, INTRINSIC :: ISO_C_BINDING,  ONLY: C_INT, C_SIZE_T, &
       C_PTR,&
       C_DOUBLE,&
       C_DOUBLE_COMPLEX,&
       C_CHAR, C_NULL_CHAR, C_NULL_PTR

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: cublasCreate
  PUBLIC :: cublasDestroy
  PUBLIC :: cublasGetVersion
  PUBLIC :: cublasDcopy
  PUBLIC :: cublasZcopy
  PUBLIC :: cublasDgemm
  PUBLIC :: cublasDger
  PUBLIC :: cublasDsyrk
  PUBLIC :: cublasDscal
  PUBLIC :: cublasZdscal
  PUBLIC :: cublasDasum
  PUBLIC :: cublasDzasum
  PUBLIC :: cublasDtrsm
  PUBLIC :: cublasSetStream
  PUBLIC :: cublasGetStream

  !cublasStatus_t
  PUBLIC :: CUBLAS_STATUS_SUCCESS

  !cublasOperation_t
  PUBLIC :: CUBLAS_OP_N
  PUBLIC :: CUBLAS_OP_T

  !cublasFillMode_t
  PUBLIC :: CUBLAS_FILL_MODE_LOWER
  PUBLIC :: CUBLAS_FILL_MODE_UPPER

  !cublasDiagType_t                                                                                                                
  PUBLIC :: CUBLAS_DIAG_NON_UNIT
  PUBLIC :: CUBLAS_DIAG_UNIT

  !cublasSideMode_t                                                                                                                
  PUBLIC :: CUBLAS_SIDE_LEFT
  PUBLIC :: CUBLAS_SIDE_RIGHT

  ENUM, BIND ( C ) ! :: cublasStatus_t
     ENUMERATOR :: CUBLAS_STATUS_SUCCESS = 0
     ENUMERATOR :: CUBLAS_STATUS_NOT_INITIALIZED = 1
     ENUMERATOR :: CUBLAS_STATUS_ALLOC_FAILED = 3
     ENUMERATOR :: CUBLAS_STATUS_INVALID_VALUE = 7
     ENUMERATOR :: CUBLAS_STATUS_ARCH_MISMATCH = 8
     ENUMERATOR :: CUBLAS_STATUS_MAPPING_ERROR = 11
     ENUMERATOR :: CUBLAS_STATUS_EXECUTION_FAILED = 13
     ENUMERATOR :: CUBLAS_STATUS_INTERNAL_ERROR = 14
     ENUMERATOR :: CUBLAS_STATUS_NOT_SUPPORTED = 15
     ENUMERATOR :: CUBLAS_STATUS_LICENSE_ERROR = 16
  END ENUM

  ENUM, BIND( C ) !:: cublasOperation_t
     ENUMERATOR :: CUBLAS_OP_N = 0
     ENUMERATOR :: CUBLAS_OP_T = 1
     ENUMERATOR :: CUBLAS_OP_C = 2
  END ENUM

  ENUM, BIND( C ) !:: cublasFillMode_t
     ENUMERATOR :: CUBLAS_FILL_MODE_LOWER = 0
     ENUMERATOR :: CUBLAS_FILL_MODE_UPPER = 1
  END ENUM

  ENUM, BIND( C ) !:: cublasDiagType_t
     ENUMERATOR :: CUBLAS_DIAG_NON_UNIT = 0
     ENUMERATOR :: CUBLAS_DIAG_UNIT = 1
  END ENUM

  ENUM, BIND( C ) !:: cublasSideMode_t
     ENUMERATOR :: CUBLAS_SIDE_LEFT = 0
     ENUMERATOR :: CUBLAS_SIDE_RIGHT = 1
  END ENUM

  TYPE, BIND( C ), PUBLIC :: cublasHandle_t
     PRIVATE
     TYPE( C_PTR ) :: cublasHandle = C_NULL_PTR
  END TYPE cublasHandle_t


  INTERFACE

     !cublasStatus_t cublasCreate(cublasHandle_t *handle)
     INTEGER( KIND( CUBLAS_STATUS_SUCCESS ) ) FUNCTION cublasCreate( handle ) BIND( C, NAME='cublasCreate_v2' )
       IMPORT :: cublasHandle_t, CUBLAS_STATUS_SUCCESS
       IMPLICIT NONE
       TYPE( cublasHandle_t ) :: handle
     END FUNCTION cublasCreate

     !cublasStatus_t cublasDestroy(cublasHandle_t handle)
     INTEGER( KIND( CUBLAS_STATUS_SUCCESS ) ) FUNCTION cublasDestroy( handle ) BIND( C, NAME='cublasDestroy_v2' )
       IMPORT :: cublasHandle_t, CUBLAS_STATUS_SUCCESS
       IMPLICIT NONE
       TYPE( cublasHandle_t ), VALUE :: handle
     END FUNCTION cublasDestroy

     !cublasStatus_t cublasGetVersion(cublasHandle_t handle, int *version)
     INTEGER( KIND( CUBLAS_STATUS_SUCCESS ) ) FUNCTION cublasGetVersion( handle, version ) BIND( C, NAME='cublasGetVersion_v2' )
       IMPORT :: C_INT, cublasHandle_t, CUBLAS_STATUS_SUCCESS
       IMPLICIT NONE
       TYPE( cublasHandle_t ), VALUE :: handle
       INTEGER( C_INT ) :: version
     END FUNCTION cublasGetVersion

     !cublasStatus_t cublasDcopy(cublasHandle_t handle, int n,
     !                           const double          *x, int incx,
     !                           double                *y, int incy)
     INTEGER( KIND( CUBLAS_STATUS_SUCCESS ) ) FUNCTION cublasDcopy( handle, n, x, incx, y, incy ) BIND( C, NAME='cublasDcopy_v2' )
       IMPORT :: C_INT, C_PTR, cublasHandle_t, CUBLAS_STATUS_SUCCESS
       IMPLICIT NONE
       TYPE( cublasHandle_t ), VALUE :: handle
       INTEGER( C_INT), VALUE :: n
       TYPE( C_PTR ), VALUE :: x
       INTEGER( C_INT), VALUE :: incx
       TYPE( C_PTR ), VALUE :: y
       INTEGER( C_INT), VALUE :: incy
     END FUNCTION cublasDcopy

     INTEGER( KIND( CUBLAS_STATUS_SUCCESS ) ) FUNCTION cublasZcopy( handle, n, x, incx, y, incy ) BIND( C, NAME='cublasZcopy_v2' )
       IMPORT :: C_INT, C_PTR, cublasHandle_t, CUBLAS_STATUS_SUCCESS
       IMPLICIT NONE
       TYPE( cublasHandle_t ), VALUE :: handle
       INTEGER( C_INT), VALUE :: n
       TYPE( C_PTR ), VALUE :: x
       INTEGER( C_INT), VALUE :: incx
       TYPE( C_PTR ), VALUE :: y
       INTEGER( C_INT), VALUE :: incy
     END FUNCTION cublasZcopy

     !cublasStatus_t cublasDgemm(cublasHandle_t handle,
     !                           cublasOperation_t transa, cublasOperation_t transb,
     !                           int m, int n, int k,
     !                           const double          *alpha,
     !                           const double          *A, int lda,
     !                           const double          *B, int ldb,
     !                           const double          *beta,
     !                           double          *C, int ldc)
     INTEGER( KIND( CUBLAS_STATUS_SUCCESS ) ) FUNCTION cublasDgemm( handle, transa, transb, &
          & m, n, k, alpha, A, lda, B, ldb, beta, C, ldc ) BIND( C, NAME='cublasDgemm_v2' )
       IMPORT :: C_INT, C_DOUBLE, C_PTR, CUBLAS_OP_N, cublasHandle_t, CUBLAS_STATUS_SUCCESS
       IMPLICIT NONE
       TYPE( cublasHandle_t ), VALUE :: handle
       INTEGER( KIND( CUBLAS_OP_N ) ), VALUE :: transa, transb
       INTEGER( C_INT), VALUE :: m, n, k
       REAL( C_DOUBLE ) :: alpha
       TYPE( C_PTR ), VALUE :: A
       INTEGER( C_INT), VALUE :: lda
       TYPE( C_PTR ), VALUE :: B
       INTEGER( C_INT), VALUE :: ldb
       REAL( C_DOUBLE ) :: beta
       TYPE( C_PTR ), VALUE :: C
       INTEGER( C_INT), VALUE :: ldc
     END FUNCTION cublasDgemm

     !cublasStatus_t  cublasDger(cublasHandle_t handle, int m, int n,
     !                           const double          *alpha,
     !                           const double          *x, int incx,
     !                           const double          *y, int incy,
     !                           double          *A, int lda)
     INTEGER( KIND( CUBLAS_STATUS_SUCCESS ) ) FUNCTION cublasDger( handle, m, n, alpha, &
          & x, incx, y, incy, A, lda ) BIND( C, NAME='cublasDger_v2' )
       IMPORT :: C_INT, C_DOUBLE, C_PTR, cublasHandle_t, CUBLAS_STATUS_SUCCESS
       IMPLICIT NONE
       TYPE( cublasHandle_t ), VALUE :: handle
       INTEGER( C_INT), VALUE :: m, n
       REAL( C_DOUBLE ) :: alpha
       TYPE( C_PTR ), VALUE :: x
       INTEGER( C_INT), VALUE :: incx
       TYPE( C_PTR ), VALUE :: y
       INTEGER( C_INT), VALUE :: incy
       TYPE( C_PTR ), VALUE :: A
       INTEGER( C_INT), VALUE :: lda
     END FUNCTION cublasDger

     !cublasStatus_t cublasDsyrk(cublasHandle_t handle,
     !                           cublasFillMode_t uplo, cublasOperation_t trans,
     !                           int n, int k,
     !                           const double          *alpha,
     !                           const double          *A, int lda,
     !                           const double          *beta,
     !                           double          *C, int ldc)
     INTEGER( KIND( CUBLAS_STATUS_SUCCESS ) ) FUNCTION cublasDsyrk( handle, uplo, trans, n, k, alpha, &
          & A, lda, beta, C, ldc ) BIND( C, NAME='cublasDsyrk_v2' )
       IMPORT :: C_INT, C_DOUBLE, C_PTR, cublasHandle_t, CUBLAS_FILL_MODE_LOWER, CUBLAS_STATUS_SUCCESS, CUBLAS_OP_N
       IMPLICIT NONE
       TYPE( cublasHandle_t ), VALUE :: handle
       INTEGER( KIND( CUBLAS_FILL_MODE_LOWER ) ), VALUE :: uplo
       INTEGER( KIND( CUBLAS_OP_N ) ), VALUE :: trans
       INTEGER( C_INT ), VALUE :: n, k
       REAL( C_DOUBLE ) :: alpha
       TYPE( C_PTR ), VALUE :: A
       INTEGER( C_INT), VALUE :: lda
       REAL( C_DOUBLE ) :: beta
       TYPE( C_PTR ), VALUE :: C
       INTEGER( C_INT), VALUE :: ldc
     END FUNCTION cublasDsyrk

     !cublasStatus_t  cublasDscal(cublasHandle_t handle, int n,
     !                            const double          *alpha,
     !                            double          *x, int incx)
     INTEGER( KIND( CUBLAS_STATUS_SUCCESS ) ) FUNCTION cublasDscal( handle, n, alpha, &
          & x, incx ) BIND( C, NAME='cublasDscal_v2' )
       IMPORT :: C_INT, C_DOUBLE, C_PTR, cublasHandle_t, CUBLAS_STATUS_SUCCESS
       IMPLICIT NONE
       TYPE( cublasHandle_t ), VALUE :: handle
       INTEGER( C_INT), VALUE :: n
       REAL( C_DOUBLE ) :: alpha
       TYPE( C_PTR ), VALUE :: x
       INTEGER( C_INT), VALUE :: incx
     END FUNCTION cublasDscal

     !cublasStatus_t cublasZdscal(cublasHandle_t handle, int n,
     !                            const double          *alpha,
     !                            cuDoubleComplex *x, int incx)
     INTEGER( KIND( CUBLAS_STATUS_SUCCESS ) ) FUNCTION cublasZdscal( handle, n, alpha, &
          & x, incx ) BIND( C, NAME='cublasZdscal_v2' )
       IMPORT :: C_INT, C_DOUBLE, C_PTR, cublasHandle_t, CUBLAS_STATUS_SUCCESS
       IMPLICIT NONE
       TYPE( cublasHandle_t ), VALUE :: handle
       INTEGER( C_INT), VALUE :: n
       REAL( C_DOUBLE ) :: alpha
       TYPE( C_PTR ), VALUE :: x
       INTEGER( C_INT), VALUE :: incx
     END FUNCTION cublasZdscal

     !cublasStatus_t cublasDasum(cublasHandle_t handle, int n,
     !                            const double          *x, int incx, double *result
     INTEGER( KIND( CUBLAS_STATUS_SUCCESS ) ) FUNCTION cublasDasum( handle, n, x, incx, RESULT ) BIND( C, NAME='cublasDasum_v2' )
       IMPORT :: C_INT, C_DOUBLE, C_PTR, cublasHandle_t, CUBLAS_STATUS_SUCCESS
       IMPLICIT NONE
       TYPE( cublasHandle_t ), VALUE :: handle
       INTEGER( C_INT), VALUE :: n
       TYPE( C_PTR ), VALUE :: x
       INTEGER( C_INT), VALUE :: incx
       REAL( C_DOUBLE ) :: RESULT
     END FUNCTION cublasDasum

     !cublasStatus_t cublasDzasum(cublasHandle_t handle, int n,
     !                            const cuDoubleComplex *x, int incx, double *result)
     INTEGER( KIND( CUBLAS_STATUS_SUCCESS ) ) FUNCTION cublasDzasum( handle, n, x, incx, RESULT ) BIND( C, NAME='cublasDzasum_v2' )
       IMPORT :: C_INT, C_DOUBLE, C_PTR, cublasHandle_t, CUBLAS_STATUS_SUCCESS
       IMPLICIT NONE
       TYPE( cublasHandle_t ), VALUE :: handle
       INTEGER( C_INT), VALUE :: n
       TYPE( C_PTR ), VALUE :: x
       INTEGER( C_INT), VALUE :: incx
       REAL( C_DOUBLE ) :: RESULT
     END FUNCTION cublasDzasum

     !cublasStatus_t cublasDtrsm (cublasHandle_t handle, cublasSideMode_t side, cublasFillMode_t uplo,
     !                            cublasOperation_t trans, cublasDiagType_t diag, int m, int n,
     !                            const double *alpha, const double *A, int lda, double *B, int ldb);
     INTEGER( KIND( CUBLAS_STATUS_SUCCESS ) ) FUNCTION cublasDtrsm( handle, side, uplo, trans, diag, m, n, &
          & alpha, A, lda, B, ldb )  BIND( C, NAME='cublasDtrsm_v2' )
       IMPORT :: cublasHandle_t, CUBLAS_STATUS_SUCCESS, C_DOUBLE, C_PTR, C_INT, CUBLAS_OP_N, &
            & CUBLAS_FILL_MODE_LOWER, CUBLAS_DIAG_UNIT, CUBLAS_SIDE_LEFT
       TYPE( cublasHandle_t ), VALUE :: handle
       INTEGER( KIND( CUBLAS_SIDE_LEFT ) ), VALUE :: side
       INTEGER( KIND( CUBLAS_FILL_MODE_LOWER ) ), VALUE :: uplo
       INTEGER( KIND( CUBLAS_OP_N ) ), VALUE :: trans
       INTEGER( KIND( CUBLAS_DIAG_UNIT ) ), VALUE :: diag
       INTEGER( C_INT ), VALUE :: m, n
       REAL( C_DOUBLE ) :: alpha
       TYPE( C_PTR ), VALUE :: A
       INTEGER( C_INT ), VALUE :: lda
       TYPE( C_PTR ), VALUE :: B
       INTEGER( C_INT ), VALUE :: ldb
     END FUNCTION cublasDtrsm

     !cublasStatus_t cublasSetStream (cublasHandle_t handle, cudaStream_t streamId);
     INTEGER( KIND( CUBLAS_STATUS_SUCCESS ) ) FUNCTION cublasSetStream ( handle, streamId ) BIND( C, NAME='cublasSetStream_v2' )
       IMPORT :: cublasHandle_t, CUBLAS_STATUS_SUCCESS, cudaStream_t
       IMPLICIT NONE
       TYPE( cublasHandle_t ), VALUE :: handle
       TYPE( cudaStream_t ), VALUE :: streamId
     END FUNCTION cublasSetStream

     !cublasStatus_t cublasGetStream(cublasHandle_t handle, cudaStream_t *streamId)
     INTEGER( KIND( CUBLAS_STATUS_SUCCESS ) ) FUNCTION cublasGetStream ( handle, streamId ) BIND( C, NAME='cublasGetStream_v2' )
       IMPORT :: cublasHandle_t, CUBLAS_STATUS_SUCCESS, cudaStream_t
       TYPE( cublasHandle_t ), VALUE :: handle
       TYPE( cudaStream_t ) :: streamId
     END FUNCTION cublasGetStream

  END INTERFACE

END MODULE cublas_interfaces
