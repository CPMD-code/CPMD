#include "cpmd_global.h"

MODULE cusolver_utils

  USE cublas_interfaces,               ONLY: CUBLAS_FILL_MODE_LOWER,&
                                             CUBLAS_FILL_MODE_UPPER
  USE cuda_types,                      ONLY: cuda_device_null,&
                                             cuda_memory_t,&
                                             cuda_stream_t
  USE cuda_utils,                      ONLY: cuda_check_device,&
                                             cuda_get_device,&
                                             cuda_set_device
  USE cusolver_interfaces,             ONLY: CUSOLVER_STATUS_SUCCESS,&
                                             cusolverDnCreate,&
                                             cusolverDnDestroy,&
                                             cusolverDnDpotrf,&
                                             cusolverDnDpotrf_bufferSize,&
                                             cusolverDnGetStream,&
                                             cusolverDnSetStream
  USE cusolver_types,                  ONLY: cusolver_handle_t
  USE error_handling,                  ONLY: stopgm

  USE, INTRINSIC :: ISO_C_BINDING,  ONLY: C_INT, C_SIZE_T, &
       C_CHAR, C_NULL_CHAR, C_NULL_PTR, C_DOUBLE
  USE kinds,                        ONLY: real_8, int_8, int_4
  USE sizeof_kinds,                 ONLY: sizeof_complex_8, sizeof_real_8, sizeof_int_4
  USE string_utils,                 ONLY: int2str


  IMPLICIT NONE

  PRIVATE

  PUBLIC :: cusolver_create
  PUBLIC :: cusolver_destroy
  PUBLIC :: cusolver_dpotrf
  PUBLIC :: cusolver_dpotrf_buffersize
  PUBLIC :: cusolver_set_stream
  PUBLIC :: cusolver_get_stream


CONTAINS


  SUBROUTINE cusolver_create ( handle, device )
    TYPE(cusolver_handle_t), INTENT(inout)   :: handle
    INTEGER, INTENT(IN)                      :: device

    CHARACTER(*), PARAMETER                  :: procedureN = 'cusolver_create'

    INTEGER                                  :: current_device
    INTEGER(KIND(CUSOLVER_STATUS_SUCCESS))   :: c_status

#if defined(_HAS_CUDA)

    IF( handle%init ) CALL stopgm(procedureN,'handle already created',&
         __LINE__,__FILE__)

    CALL cuda_check_device ( device, procedureN )

    CALL cuda_get_device ( current_device )
    IF( current_device == device ) CALL cuda_set_device ( device )
    c_status = cusolverDnCreate( handle%h )
    IF( c_status /= CUSOLVER_STATUS_SUCCESS ) CALL stopgm(procedureN,"cusolver error: "//TRIM(int2str( INT( c_status ) )),&
         __LINE__,__FILE__)

    handle%device = device
    handle%init = .TRUE.

#else

    CALL stopgm(procedureN,"no cusolver available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE cusolver_create


  SUBROUTINE cusolver_destroy ( handle )
    TYPE(cusolver_handle_t), INTENT(inout)   :: handle

    CHARACTER(*), PARAMETER :: procedureN = 'cusolver_destroy'

    INTEGER(KIND(CUSOLVER_STATUS_SUCCESS))   :: c_status

#if defined(_HAS_CUDA)

    IF( .NOT. handle%init ) CALL stopgm(procedureN,'handle not created',&
         __LINE__,__FILE__)

    c_status = cusolverDnDestroy( handle%h )
    IF( c_status /= CUSOLVER_STATUS_SUCCESS ) CALL stopgm(procedureN,"cusolver error: "//TRIM(int2str( INT( c_status ) )),&
         __LINE__,__FILE__)

    handle%device = cuda_device_null
    handle%init = .FALSE.

#else

    CALL stopgm(procedureN,"no cusolver available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE cusolver_destroy


  SUBROUTINE cusolver_dpotrf ( handle, uplo, n, A, lda, Workspace, Lwork, devInfo )
    TYPE(cusolver_handle_t), INTENT(in)      :: handle
    CHARACTER(len=1), INTENT(in)             :: uplo
    INTEGER, INTENT(IN)                      :: n
    TYPE(cuda_memory_t)                      :: A
    INTEGER, INTENT(IN)                      :: lda
    TYPE(cuda_memory_t)                      :: Workspace
    INTEGER, INTENT(IN)                      :: Lwork
    INTEGER, INTENT(OUT)                     :: devInfo

    CHARACTER(*), PARAMETER                  :: procedureN = 'cusolver_dpotrf'

    INTEGER(C_INT)                           :: c_devInfo, c_lda, c_Lwork, &
                                                c_n, c_uplo
    INTEGER(KIND(CUSOLVER_STATUS_SUCCESS))   :: c_status

    devInfo = 0

#if defined(_HAS_CUDA)

    IF( .NOT. handle%init ) CALL stopgm(procedureN,'handle not created',&
         __LINE__,__FILE__)

    IF( .NOT. A%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT. Workspace%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF(  handle%device /= A%device .OR. &
         handle%device /=Workspace%device ) CALL stopgm(procedureN,'stream and handle dont share the same device',&
         __LINE__,__FILE__)

    CALL cuda_check_device ( handle%device, procedureN )

    SELECT CASE( uplo )
    CASE( 'U', 'u' ); c_uplo = CUBLAS_FILL_MODE_UPPER
    CASE( 'L', 'l' ); c_uplo = CUBLAS_FILL_MODE_LOWER
    CASE DEFAULT; c_uplo = HUGE( 0_C_INT )
    END SELECT
    c_n = INT( n, C_INT )
    c_lda = INT( lda, C_INT )
    c_Lwork = INT( Lwork, C_INT )
    !write(*,*) 'c_uplo',c_uplo,' c_n',c_n,' c_lda',c_lda,' c_Lwork',c_Lwork,' A',A%n_bytes,' Workspace',Workspace%n_bytes
    c_status = cusolverDnDpotrf ( handle%h, c_uplo, c_n, A%ptr, c_lda, Workspace%ptr, c_Lwork, c_devInfo )
    IF( c_status /= CUSOLVER_STATUS_SUCCESS ) CALL stopgm(procedureN,"cusolver error: "//TRIM(int2str( INT( c_status ) ))&
         & //' devInfo '//int2str( INT( c_devInfo ) ),&
         __LINE__,__FILE__)

    devInfo = INT( c_devInfo )

#else

    CALL stopgm(procedureN,"no cusolver available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE cusolver_dpotrf


  SUBROUTINE cusolver_dpotrf_buffersize ( handle, uplo, n, A, lda, Lwork )
    TYPE(cusolver_handle_t), INTENT(in)      :: handle
    CHARACTER(len=1), INTENT(in)             :: uplo
    INTEGER, INTENT(IN)                      :: n
    TYPE(cuda_memory_t)                      :: A
    INTEGER, INTENT(IN)                      :: lda
    INTEGER, INTENT(OUT)                     :: Lwork

    CHARACTER(*), PARAMETER :: procedureN = 'cusolver_dpotrf_buffersize'

    INTEGER(C_INT)                           :: c_lda, c_Lwork, c_n, c_uplo
    INTEGER(KIND(CUSOLVER_STATUS_SUCCESS))   :: c_status

#if defined(_HAS_CUDA)

    IF( .NOT. handle%init ) CALL stopgm(procedureN,'handle not created',&
         __LINE__,__FILE__)

    IF( .NOT. A%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( handle%device /= A%device ) CALL stopgm(procedureN,'stream and handle dont share the same device',&
         __LINE__,__FILE__)

    CALL cuda_check_device ( handle%device, procedureN )

    SELECT CASE( uplo )
    CASE( 'U', 'u' ); c_uplo = CUBLAS_FILL_MODE_UPPER
    CASE( 'L', 'l' ); c_uplo = CUBLAS_FILL_MODE_LOWER
    END SELECT
    c_n = INT( n, C_INT )
    c_lda = INT( lda, C_INT )

    c_status = cusolverDnDpotrf_bufferSize( handle%h, c_uplo, c_n, A%ptr, c_lda, c_Lwork )
    IF( c_status /= CUSOLVER_STATUS_SUCCESS ) CALL stopgm(procedureN,"cusolver error: "//TRIM(int2str( INT( c_status ) )),&
         __LINE__,__FILE__)

    Lwork = INT( c_Lwork )

#else

    CALL stopgm(procedureN,"no cusolver available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE cusolver_dpotrf_buffersize


  SUBROUTINE cusolver_set_stream ( handle, stream )
    TYPE(cusolver_handle_t), INTENT(in)      :: handle
    TYPE(cuda_stream_t)                      :: stream

    CHARACTER(*), PARAMETER :: procedureN = 'cusolver_set_stream'

    INTEGER(KIND(CUSOLVER_STATUS_SUCCESS))   :: c_status

#if defined(_HAS_CUDA)

    IF( .NOT. handle%init ) CALL stopgm(procedureN,'handle not created',&
         __LINE__,__FILE__)

    IF( .NOT. stream%init ) CALL stopgm(procedureN,'stream not allocated',&
         __LINE__,__FILE__)

    IF( handle%device /= stream%device ) CALL stopgm(procedureN,'stream and handle dont share the same device',&
         __LINE__,__FILE__)

    CALL cuda_check_device ( handle%device, procedureN )

    c_status = cusolverDnSetStream ( handle%h, stream%s )
    IF( c_status /= CUSOLVER_STATUS_SUCCESS ) CALL stopgm(procedureN,"cusolver error: "//TRIM(int2str( INT( c_status ) )),&
         __LINE__,__FILE__)

#else

    CALL stopgm(procedureN,"no cusolver available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE cusolver_set_stream


  SUBROUTINE cusolver_get_stream ( handle, stream )
    TYPE(cusolver_handle_t), INTENT(in)      :: handle
    TYPE(cuda_stream_t)                      :: stream

    CHARACTER(*), PARAMETER :: procedureN = 'cusolver_get_stream'

    INTEGER(KIND(CUSOLVER_STATUS_SUCCESS))   :: c_status

#if defined(_HAS_CUDA)

    IF( .NOT. handle%init ) CALL stopgm(procedureN,'handle not created',&
         __LINE__,__FILE__)

    IF( stream%init ) CALL stopgm(procedureN,'stream already allocated',&
         __LINE__,__FILE__)

    CALL cuda_check_device ( handle%device, procedureN )

    c_status = cusolverDnGetStream ( handle%h, stream%s )
    IF( c_status /= CUSOLVER_STATUS_SUCCESS ) CALL stopgm(procedureN,"cusolver error: "//TRIM(int2str( INT( c_status ) )),&
         __LINE__,__FILE__)

    stream%device = handle%device
    stream%init = .TRUE.

#else

    CALL stopgm(procedureN,"no cusolver available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE cusolver_get_stream


END MODULE cusolver_utils
