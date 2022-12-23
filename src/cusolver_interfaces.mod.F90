MODULE cusolver_interfaces

  USE cublas_interfaces,               ONLY: CUBLAS_FILL_MODE_LOWER
  USE cuda_interfaces,                 ONLY: cudaStream_t

  USE, INTRINSIC :: ISO_C_BINDING,  ONLY: C_INT, C_SIZE_T, &
       C_PTR,&
       C_DOUBLE,&
       C_DOUBLE_COMPLEX,&
       C_CHAR, C_NULL_CHAR, C_NULL_PTR

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: cusolverDnCreate
  PUBLIC :: cusolverDnDestroy
  PUBLIC :: cusolverDnDpotrf
  PUBLIC :: cusolverDnDpotrf_bufferSize
  PUBLIC :: cusolverDnSetStream
  PUBLIC :: cusolverDnGetStream

  !cusolverStatus_t
  PUBLIC :: CUSOLVER_STATUS_SUCCESS


  ENUM, BIND( C ) !:: cusolverStatus_t
     ENUMERATOR :: CUSOLVER_STATUS_SUCCESS = 0
     ENUMERATOR :: CUSOLVER_STATUS_NOT_INITIALIZED = 1
     ENUMERATOR :: CUSOLVER_STATUS_ALLOC_FAILED = 2
     ENUMERATOR :: CUSOLVER_STATUS_INVALID_VALUE = 3
     ENUMERATOR :: CUSOLVER_STATUS_ARCH_MISMATCH = 4
     ENUMERATOR :: CUSOLVER_STATUS_MAPPING_ERROR = 5
     ENUMERATOR :: CUSOLVER_STATUS_EXECUTION_FAILED = 6
     ENUMERATOR :: CUSOLVER_STATUS_INTERNAL_ERROR = 7
     ENUMERATOR :: CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED = 8
     ENUMERATOR :: CUSOLVER_STATUS_NOT_SUPPORTED = 9
     ENUMERATOR :: CUSOLVER_STATUS_ZERO_PIVOT = 10
     ENUMERATOR :: CUSOLVER_STATUS_INVALID_LICENSE = 11
  END ENUM

  TYPE, BIND( C ), PUBLIC :: cusolverDnHandle_t
     PRIVATE
     TYPE( C_PTR ) :: cusolverDnHandle = C_NULL_PTR
  END TYPE cusolverDnHandle_t


  INTERFACE

     !cusolverStatus_t cusolverDnCreate(cusolverDnHandle_t *handle);
     INTEGER( KIND( CUSOLVER_STATUS_SUCCESS ) ) FUNCTION cusolverDnCreate( handle ) BIND( C, NAME='cusolverDnCreate' )
       IMPORT :: cusolverDnHandle_t, CUSOLVER_STATUS_SUCCESS
       IMPLICIT NONE
       TYPE( cusolverDnHandle_t ) :: handle
     END FUNCTION cusolverDnCreate

     !cusolverStatus_t cusolverDnDestroy(cusolverDnHandle_t handle);
     INTEGER( KIND( CUSOLVER_STATUS_SUCCESS ) ) FUNCTION cusolverDnDestroy( handle ) BIND( C, NAME='cusolverDnDestroy' )
       IMPORT :: cusolverDnHandle_t, CUSOLVER_STATUS_SUCCESS
       IMPLICIT NONE
       TYPE( cusolverDnHandle_t ), VALUE :: handle
     END FUNCTION cusolverDnDestroy

     !cusolverStatus_t cusolverDnSetStream (cusolverDnHandle_t handle, cudaStream_t streamId);
     INTEGER( KIND( CUSOLVER_STATUS_SUCCESS ) ) FUNCTION cusolverDnSetStream( handle, streamId ) &
          & BIND( C, NAME='cusolverDnSetStream' )
       IMPORT :: cusolverDnHandle_t, CUSOLVER_STATUS_SUCCESS, cudaStream_t
       IMPLICIT NONE
       TYPE( cusolverDnHandle_t ), VALUE :: handle
       TYPE( cudaStream_t ), VALUE :: streamId
     END FUNCTION cusolverDnSetStream

     !cusolverStatus_t cusolverDnGetStream(cusolverDnHandle_t handle, cudaStream_t *streamId);
     INTEGER( KIND( CUSOLVER_STATUS_SUCCESS ) ) FUNCTION cusolverDnGetStream( handle, streamId ) &
          & BIND( C, NAME='cusolverDnGetStream' )
       IMPORT :: cusolverDnHandle_t, CUSOLVER_STATUS_SUCCESS, cudaStream_t
       IMPLICIT NONE
       TYPE( cusolverDnHandle_t ), VALUE :: handle
       TYPE( cudaStream_t ) :: streamId
     END FUNCTION cusolverDnGetStream

     !cusolverStatus_t cusolverDnDpotrf_bufferSize(cusolveDnHandle_t handle,
     !                                             cublasFillMode_t uplo,
     !                                             int n,
     !                                             double *A,
     !                                             int lda,
     !                                             int *Lwork );
     INTEGER( KIND( CUSOLVER_STATUS_SUCCESS ) ) FUNCTION cusolverDnDpotrf_bufferSize( handle, uplo, &
          & n, A, lda, Lwork )  BIND( C, NAME='cusolverDnDpotrf_bufferSize' )
       IMPORT :: cusolverDnHandle_t, CUSOLVER_STATUS_SUCCESS, CUBLAS_FILL_MODE_LOWER, C_INT, C_PTR
       IMPLICIT NONE
       TYPE( cusolverDnHandle_t ), VALUE :: handle
       INTEGER( KIND( CUBLAS_FILL_MODE_LOWER ) ), VALUE :: uplo
       INTEGER( C_INT ), VALUE :: n, lda
       TYPE( C_PTR ), VALUE :: A
       INTEGER( C_INT ) :: Lwork
     END FUNCTION cusolverDnDpotrf_bufferSize

     !cusolverStatus_t cusolverDnDpotrf(cusolverDnHandle_t handle, 
     !                                  cublasFillMode_t uplo,
     !                                  int n,
     !                                  double *A,
     !                                  int lda, 
     !                                  double *Workspace,
     !                                  int Lwork,
     !                                  int *devInfo );
     INTEGER( KIND( CUSOLVER_STATUS_SUCCESS ) ) FUNCTION cusolverDnDpotrf( handle, uplo, n, A, lda, &
          & Workspace, Lwork, devInfo ) BIND( C, NAME='cusolverDnDpotrf' )
       IMPORT :: cusolverDnHandle_t, CUSOLVER_STATUS_SUCCESS, CUBLAS_FILL_MODE_LOWER, C_INT, C_PTR
       IMPLICIT NONE
       TYPE( cusolverDnHandle_t ), VALUE :: handle
       INTEGER( KIND( CUBLAS_FILL_MODE_LOWER ) ), VALUE :: uplo
       INTEGER( C_INT ), VALUE :: n, lda, Lwork
       TYPE( C_PTR ), VALUE :: A, Workspace
       INTEGER( C_INT ) :: devInfo
     END FUNCTION cusolverDnDpotrf

  END INTERFACE

END MODULE cusolver_interfaces
