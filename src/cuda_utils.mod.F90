#include "cpmd_global.h"

MODULE cuda_utils

  USE cuda_interfaces,                 ONLY: &
       cGetMemAddrs, cMemPointsToDbl, cudaComputeModeDefault, &
       cudaComputeModeExclusive, cudaComputeModeExclusiveProcess, &
       cudaComputeModeProhibited, cudaCpuDeviceId, &
       cudaDeviceGetStreamPriorityRange, cudaDeviceProp_80_t, &
       cudaDeviceReset, cudaDeviceSynchronize, cudaDriverGetVersion, &
       cudaErrorHostMemoryAlreadyRegistered, cudaErrorNotReady, &
       cudaEventCreate, cudaEventCreateWithFlags, cudaEventDestroy, &
       cudaEventDisableTiming, cudaEventElapsedTime, cudaEventRecord, &
       cudaEventSynchronize, cudaFree, cudaFreeHost, cudaGetDevice, &
       cudaGetDeviceCount, cudaGetDeviceProperties_80, cudaHostAlloc, &
       cudaHostAllocPortable, cudaHostRegister, cudaHostRegisterPortable, &
       cudaHostUnregister, cudaMalloc, cudaMallocHost, cudaMallocManaged, &
       cudaMallocPitch, cudaMemAdvise, cudaMemAdviseSetReadMostly, &
       cudaMemAttachGlobal, cudaMemPrefetchAsync, cudaMemcpy, cudaMemcpy2D, &
       cudaMemcpy2DAsync, cudaMemcpyAsync, cudaMemcpyDeviceToDevice, &
       cudaMemcpyDeviceToHost, cudaMemcpyHostToDevice, cudaMemoryTypeDevice, &
       cudaMemoryTypeHost, cudaMemset, cudaPointerAttributes_t, &
       cudaPointerGetAttributes, cudaRuntimeGetVersion, cudaSetDevice, &
       cudaStreamCreate, cudaStreamCreateWithPriority, cudaStreamDestroy, &
       cudaStreamNonBlocking, cudaStreamQuery, cudaStreamSynchronize, &
       cudaStreamWaitEvent, cudaSuccess, cudaGetErrorString, cudaDeviceGetPCIBusId, &
       cudaMemGetInfo
  USE cuda_types,                      ONLY: cuda_device_null,&
                                             cuda_event_t,&
                                             cuda_memory_t,&
                                             cuda_stream_t
  USE error_handling,                  ONLY: stopgm

  USE, INTRINSIC :: ISO_C_BINDING,  ONLY: C_INT, C_SIZE_T, &
       C_CHAR, C_NULL_CHAR, C_NULL_PTR, C_LOC, C_F_POINTER, C_PTR, &
       C_FLOAT
  USE kinds,                        ONLY: real_8, int_8, int_4, int_8
  USE sizeof_kinds,                 ONLY: sizeof_complex_8, sizeof_real_8, sizeof_int_4
  USE string_utils,                 ONLY: int2str, c_strlen, remove_c_null_char

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: cuda_driver_get_version
  PUBLIC :: cuda_runtime_get_version
  PUBLIC :: cuda_alloc_bytes
  PUBLIC :: cuda_dealloc
  PUBLIC :: cuda_alloc_pitch_bytes
  PUBLIC :: cuda_alloc_host
  PUBLIC :: cuda_alloc_host_portable
  PUBLIC :: cuda_dealloc_host
  PUBLIC :: cuda_host_register
  PUBLIC :: cuda_host_unregister
  PUBLIC :: cuda_memcpy_host_to_device
  PUBLIC :: cuda_memcpy_device_to_host
  PUBLIC :: cuda_memcpy2d_host_to_device
  PUBLIC :: cuda_memcpy2d_device_to_host
  PUBLIC :: cuda_memcpy_device_to_device
  PUBLIC :: cuda_memcpy_async_host_to_device
  PUBLIC :: cuda_memcpy_async_device_to_host
  PUBLIC :: cuda_memcpy2d_async_host_to_device
  PUBLIC :: cuda_memcpy2d_async_device_to_host
  PUBLIC :: cuda_mem_zero_bytes
  PUBLIC :: cuda_d_points_to
  PUBLIC :: cuda_z_points_to
  PUBLIC :: cuda_get_address
  PUBLIC :: cuda_stream_create
  PUBLIC :: cuda_stream_create_with_priority
  PUBLIC :: cuda_stream_destroy
  PUBLIC :: cuda_stream_synchronize
  PUBLIC :: cuda_stream_query
  PUBLIC :: cuda_stream_wait_event
  PUBLIC :: cuda_event_create
  PUBLIC :: cuda_event_create_with_no_timing
  PUBLIC :: cuda_event_destroy
  PUBLIC :: cuda_event_synchronize
  PUBLIC :: cuda_event_elapsed_time
  PUBLIC :: cuda_event_record
  PUBLIC :: cuda_pointer_get_attributes
  PUBLIC :: cuda_get_device_properties
  PUBLIC :: cuda_device_reset
  PUBLIC :: cuda_device_synchronize
  PUBLIC :: cuda_get_device_count
  PUBLIC :: cuda_get_device
  PUBLIC :: cuda_set_device
  PUBLIC :: cuda_check_device
  PUBLIC :: cuda_device_get_stream_priority_range
  PUBLIC :: cuda_mem_prefetch_async
  PUBLIC :: cuda_mem_advise
  PUBLIC :: cuda_alloc_managed
  PUBLIC :: cuda_dealloc_managed
  PUBLIC :: cuda_device_get_pci_bus_id

  INTERFACE cuda_host_register
     MODULE PROCEDURE cuda_host_register_int4_r1
     MODULE PROCEDURE cuda_host_register_int4_r2
     MODULE PROCEDURE cuda_host_register_real8_r1
     MODULE PROCEDURE cuda_host_register_real8_r2
     MODULE PROCEDURE cuda_host_register_real8_r3
     MODULE PROCEDURE cuda_host_register_complex8_r1
     MODULE PROCEDURE cuda_host_register_complex8_r2
  END INTERFACE cuda_host_register

  INTERFACE cuda_host_unregister
     MODULE PROCEDURE cuda_host_unregister_int4_r1
     MODULE PROCEDURE cuda_host_unregister_int4_r2
     MODULE PROCEDURE cuda_host_unregister_real8_r1
     MODULE PROCEDURE cuda_host_unregister_real8_r2
     MODULE PROCEDURE cuda_host_unregister_real8_r3
     MODULE PROCEDURE cuda_host_unregister_complex8_r1
     MODULE PROCEDURE cuda_host_unregister_complex8_r2
  END INTERFACE cuda_host_unregister

  INTERFACE cuda_alloc_host
     MODULE PROCEDURE cuda_alloc_host_int4_r1
     MODULE PROCEDURE cuda_alloc_host_int4_r2
     MODULE PROCEDURE cuda_alloc_host_real8_r1
     MODULE PROCEDURE cuda_alloc_host_real8_r2
     MODULE PROCEDURE cuda_alloc_host_real8_r3
     MODULE PROCEDURE cuda_alloc_host_complex8_r1
     MODULE PROCEDURE cuda_alloc_host_complex8_r2
  END INTERFACE cuda_alloc_host

  INTERFACE cuda_alloc_host_portable
     MODULE PROCEDURE cuda_alloc_host_portable_int4_r1
     MODULE PROCEDURE cuda_alloc_host_portable_int4_r2
     MODULE PROCEDURE cuda_alloc_host_portable_real8_r1
     MODULE PROCEDURE cuda_alloc_host_portable_real8_r2
     MODULE PROCEDURE cuda_alloc_host_portable_real8_r3
     MODULE PROCEDURE cuda_alloc_host_portable_complex8_r1
     MODULE PROCEDURE cuda_alloc_host_portable_complex8_r2
  END INTERFACE cuda_alloc_host_portable

  INTERFACE cuda_dealloc_host
     MODULE PROCEDURE cuda_dealloc_host_int4_r1
     MODULE PROCEDURE cuda_dealloc_host_int4_r2
     MODULE PROCEDURE cuda_dealloc_host_real8_r1
     MODULE PROCEDURE cuda_dealloc_host_real8_r2
     MODULE PROCEDURE cuda_dealloc_host_real8_r3
     MODULE PROCEDURE cuda_dealloc_host_complex8_r1
     MODULE PROCEDURE cuda_dealloc_host_complex8_r2
  END INTERFACE cuda_dealloc_host


  INTERFACE cuda_memcpy_host_to_device
     MODULE PROCEDURE cuda_memcpy_host_to_device_int4_r1
     MODULE PROCEDURE cuda_memcpy_host_to_device_int4_r2
     MODULE PROCEDURE cuda_memcpy_host_to_device_real8_r1
     MODULE PROCEDURE cuda_memcpy_host_to_device_real8_r2
     MODULE PROCEDURE cuda_memcpy_host_to_device_complex8_r1
     MODULE PROCEDURE cuda_memcpy_host_to_device_complex8_r2
  END INTERFACE cuda_memcpy_host_to_device


  INTERFACE cuda_memcpy_device_to_host
     MODULE PROCEDURE cuda_memcpy_device_to_host_int4_r1
     MODULE PROCEDURE cuda_memcpy_device_to_host_int4_r2
     MODULE PROCEDURE cuda_memcpy_device_to_host_real8_r1
     MODULE PROCEDURE cuda_memcpy_device_to_host_real8_r2
     MODULE PROCEDURE cuda_memcpy_device_to_host_complex8_r1
     MODULE PROCEDURE cuda_memcpy_device_to_host_complex8_r2
  END INTERFACE cuda_memcpy_device_to_host


  INTERFACE cuda_memcpy2d_host_to_device
     MODULE PROCEDURE cuda_memcpy2d_host_to_device_int4_r1
     MODULE PROCEDURE cuda_memcpy2d_host_to_device_int4_r2
     MODULE PROCEDURE cuda_memcpy2d_host_to_device_real8_r1
     MODULE PROCEDURE cuda_memcpy2d_host_to_device_real8_r2
     MODULE PROCEDURE cuda_memcpy2d_host_to_device_complex8_r1
     MODULE PROCEDURE cuda_memcpy2d_host_to_device_complex8_r2
  END INTERFACE cuda_memcpy2d_host_to_device


  INTERFACE cuda_memcpy2d_device_to_host
     MODULE PROCEDURE cuda_memcpy2d_device_to_host_int4_r1
     MODULE PROCEDURE cuda_memcpy2d_device_to_host_int4_r2
     MODULE PROCEDURE cuda_memcpy2d_device_to_host_real8_r1
     MODULE PROCEDURE cuda_memcpy2d_device_to_host_real8_r2
     MODULE PROCEDURE cuda_memcpy2d_device_to_host_complex8_r1
     MODULE PROCEDURE cuda_memcpy2d_device_to_host_complex8_r2
  END INTERFACE cuda_memcpy2d_device_to_host


  INTERFACE cuda_memcpy_async_host_to_device
     MODULE PROCEDURE cuda_memcpy_async_host_to_device_int4_r1
     MODULE PROCEDURE cuda_memcpy_async_host_to_device_int4_r2
     MODULE PROCEDURE cuda_memcpy_async_host_to_device_real8_r1
     MODULE PROCEDURE cuda_memcpy_async_host_to_device_real8_r2
     MODULE PROCEDURE cuda_memcpy_async_host_to_device_complex8_r1
     MODULE PROCEDURE cuda_memcpy_async_host_to_device_complex8_r2
  END INTERFACE cuda_memcpy_async_host_to_device


  INTERFACE cuda_memcpy_async_device_to_host
     MODULE PROCEDURE  cuda_memcpy_async_device_to_host_int4_r1
     MODULE PROCEDURE  cuda_memcpy_async_device_to_host_int4_r2
     MODULE PROCEDURE  cuda_memcpy_async_device_to_host_real8_r1
     MODULE PROCEDURE  cuda_memcpy_async_device_to_host_real8_r2
     MODULE PROCEDURE  cuda_memcpy_async_device_to_host_complex8_r1
     MODULE PROCEDURE  cuda_memcpy_async_device_to_host_complex8_r2
  END INTERFACE cuda_memcpy_async_device_to_host


  INTERFACE cuda_memcpy2d_async_host_to_device
     MODULE PROCEDURE cuda_memcpy2d_async_host_to_device_int4_r1
     MODULE PROCEDURE cuda_memcpy2d_async_host_to_device_int4_r2
     MODULE PROCEDURE cuda_memcpy2d_async_host_to_device_real8_r1
     MODULE PROCEDURE cuda_memcpy2d_async_host_to_device_real8_r2
     MODULE PROCEDURE cuda_memcpy2d_async_host_to_device_complex8_r1
     MODULE PROCEDURE cuda_memcpy2d_async_host_to_device_complex8_r2
  END INTERFACE cuda_memcpy2d_async_host_to_device


  INTERFACE cuda_memcpy2d_async_device_to_host
     MODULE PROCEDURE cuda_memcpy2d_async_device_to_host_int4_r1
     MODULE PROCEDURE cuda_memcpy2d_async_device_to_host_int4_r2
     MODULE PROCEDURE cuda_memcpy2d_async_device_to_host_real8_r1
     MODULE PROCEDURE cuda_memcpy2d_async_device_to_host_real8_r2
     MODULE PROCEDURE cuda_memcpy2d_async_device_to_host_complex8_r1
     MODULE PROCEDURE cuda_memcpy2d_async_device_to_host_complex8_r2
  END INTERFACE cuda_memcpy2d_async_device_to_host


  INTERFACE cuda_alloc_managed
     MODULE PROCEDURE cuda_alloc_managed_int4_r1
     MODULE PROCEDURE cuda_alloc_managed_int4_r2
     MODULE PROCEDURE cuda_alloc_managed_real8_r1
     MODULE PROCEDURE cuda_alloc_managed_real8_r2
     MODULE PROCEDURE cuda_alloc_managed_real8_r3
     MODULE PROCEDURE cuda_alloc_managed_complex8_r1
     MODULE PROCEDURE cuda_alloc_managed_complex8_r2
  END INTERFACE cuda_alloc_managed

  INTERFACE cuda_dealloc_managed
     MODULE PROCEDURE cuda_dealloc_managed_int4_r1
     MODULE PROCEDURE cuda_dealloc_managed_int4_r2
     MODULE PROCEDURE cuda_dealloc_managed_real8_r1
     MODULE PROCEDURE cuda_dealloc_managed_real8_r2
     MODULE PROCEDURE cuda_dealloc_managed_real8_r3
     MODULE PROCEDURE cuda_dealloc_managed_complex8_r1
     MODULE PROCEDURE cuda_dealloc_managed_complex8_r2
  END INTERFACE cuda_dealloc_managed

  integer, private, parameter :: cuda_max_len_error_message = 1024
  logical, private, parameter :: get_cuda_mem_info = .true.
  integer(int_8), private, save :: cuda_mem_alloc_total = 0_c_size_t

CONTAINS


  function cuda_get_error_string ( c_status ) result( message )

    INTEGER(KIND(cudaSuccess)), INTENT(IN) :: c_status
    CHARACTER(len=cuda_max_len_error_message) :: message

    CHARACTER(*), PARAMETER :: procedureN = 'cuda_get_error_string'

    integer :: n, i
    CHARACTER(C_CHAR), pointer :: f_message(:)

    TYPE(C_PTR) :: c_message

    message = ""

#if defined(_HAS_CUDA)

    c_message = cudaGetErrorString ( c_status )

    n = INT( c_strlen ( c_message ) )

    CALL C_F_POINTER( c_message, f_message, [ n ] )

    do i = 1, n
       message(i:i) = f_message(i)
    enddo

#endif

  end function cuda_get_error_string


  SUBROUTINE cuda_check_device ( device, calling_procedure )
    INTEGER, INTENT(IN)                      :: device
    CHARACTER(*), INTENT(IN)                 :: calling_procedure

    CHARACTER(*), PARAMETER :: procedureN = 'cuda_check_device'

    INTEGER                                  :: n_devices

    IF( device == cuda_device_null ) CALL stopgm(procedureN,calling_procedure//' device is NULL',&
         __LINE__,__FILE__)

    IF( device < 0 ) CALL stopgm(procedureN,calling_procedure//' device is smaller than 0',&
         __LINE__,__FILE__)

    CALL cuda_get_device_count ( n_devices )
    IF( device >= n_devices ) CALL stopgm(procedureN,calling_procedure//' device is greater than max number device',&
         __LINE__,__FILE__)

  END SUBROUTINE cuda_check_device


  SUBROUTINE cuda_driver_get_version ( driver_version )
    INTEGER, INTENT(OUT)                     :: driver_version

    CHARACTER(*), PARAMETER :: procedureN = 'cuda_driver_get_version'

    INTEGER(C_INT)                           :: c_driver_version
    INTEGER(KIND(cudaSuccess))               :: c_status

    driver_version = 0

#if defined(_HAS_CUDA)

    c_status = cudaDriverGetVersion( c_driver_version )
    IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))&
         //". "//trim(cuda_get_error_string(c_status)),&
         __LINE__,__FILE__)

    driver_version = INT( c_driver_version )

#endif

  END SUBROUTINE cuda_driver_get_version


  SUBROUTINE cuda_runtime_get_version ( runtime_version )
    INTEGER, INTENT(OUT)                     :: runtime_version

    CHARACTER(*), PARAMETER :: procedureN = 'cuda_runtime_get_version'

    INTEGER(C_INT)                           :: c_runtime_version
    INTEGER(KIND(cudaSuccess))               :: c_status

    runtime_version = 0

#if defined(_HAS_CUDA)

    c_status = cudaRuntimeGetVersion( c_runtime_version )
    IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))&
         //". "//trim(cuda_get_error_string(c_status)),&
         __LINE__,__FILE__)

    runtime_version = INT( c_runtime_version )

#endif

  END SUBROUTINE cuda_runtime_get_version


  SUBROUTINE cuda_pointer_get_attributes ( mem )
    TYPE(cuda_memory_t), INTENT(IN)          :: mem

    CHARACTER(*), PARAMETER :: procedureN = 'cuda_pointer_get_attributes'

    INTEGER(KIND(cudaSuccess))               :: c_status
    TYPE(cudaPointerAttributes_t)            :: attributes

#if defined(_HAS_CUDA)

    IF( .NOT. mem%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    c_status = cudaPointerGetAttributes ( attributes, mem%ptr )
    IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))&
         //". "//trim(cuda_get_error_string(c_status)),&
         __LINE__,__FILE__)

    WRITE(*,*) 'CUDA: memory attribute'
    SELECT CASE( attributes%memoryType )
    CASE( cudaMemoryTypeHost )
       WRITE(*,*) 'CUDA:    type                   ', 'cudaMemoryTypeHost'
    CASE( cudaMemoryTypeDevice )
       WRITE(*,*) 'CUDA:    type                   ', 'cudaMemoryTypeDevice'
    CASE DEFAULT
       WRITE(*,*) 'CUDA:    type                   ', 'unknown'
    END SELECT
    WRITE(*,*)    'CUDA:    device                 ', attributes%memoryType
    WRITE(*,*)    'CUDA:    isManaged              ', attributes%isManaged == 1

#else

    CALL stopgm(procedureN,"no cuda available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE cuda_pointer_get_attributes


  SUBROUTINE cuda_mem_prefetch_async ( mem, n_bytes, stream, is_CPUDevice )
    TYPE(cuda_memory_t), INTENT(IN)          :: mem
    INTEGER(int_8), INTENT(IN)               :: n_bytes
    TYPE(cuda_stream_t), INTENT(IN)          :: stream
    LOGICAL, INTENT(IN), OPTIONAL            :: is_CPUDevice

    CHARACTER(*), PARAMETER :: procedureN = 'cuda_mem_prefetch_async'

    INTEGER(C_SIZE_T)                        :: c_n_bytes
    INTEGER(KIND(cudaSuccess))               :: c_status
    LOGICAL                                  :: my_is_CPUDevice

    my_is_CPUDevice = .FALSE.
    IF( PRESENT( is_CPUDevice ) ) my_is_CPUDevice = is_CPUDevice

#if defined(_HAS_CUDA)

    IF( .NOT. mem%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT. stream%init ) CALL stopgm(procedureN,'stream not initialized',&
         __LINE__,__FILE__)

    IF( mem%device /= stream%device ) CALL stopgm(procedureN,'stream and memory dont share the same device',&
         __LINE__,__FILE__)

    c_n_bytes = INT( n_bytes, C_SIZE_T )

    IF( my_is_CPUDevice ) THEN
       c_status = cudaMemPrefetchAsync ( mem%ptr, c_n_bytes, cudaCpuDeviceId, stream%s )
    ELSE
       c_status = cudaMemPrefetchAsync ( mem%ptr, c_n_bytes, stream%device, stream%s )
    ENDIF
    IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))&
         //". "//trim(cuda_get_error_string(c_status)),&
         __LINE__,__FILE__)

#else

    CALL stopgm(procedureN,"no cuda available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE cuda_mem_prefetch_async


  SUBROUTINE cuda_mem_advise ( mem, n_bytes, advice )
    TYPE(cuda_memory_t), INTENT(INOUT)       :: mem
    INTEGER(int_8), INTENT(IN)               :: n_bytes
    INTEGER&
      (KIND(cudaMemAdviseSetReadMostly)), &
      INTENT(IN)                             :: advice

    CHARACTER(*), PARAMETER                  :: procedureN = 'cuda_mem_advise'

    INTEGER(C_SIZE_T)                        :: c_n_bytes
    INTEGER(KIND(cudaSuccess))               :: c_status

#if defined(_HAS_CUDA)

    IF( .NOT. mem%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    c_n_bytes = INT( n_bytes, C_SIZE_T )

    c_status = cudaMemAdvise ( mem%ptr, c_n_bytes, advice, mem%device )
    IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))&
         //". "//trim(cuda_get_error_string(c_status)),&
         __LINE__,__FILE__)

#else

    CALL stopgm(procedureN,"no cuda available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE cuda_mem_advise


  SUBROUTINE cuda_get_device_properties ( device )
    INTEGER, INTENT(IN)                      :: device

    CHARACTER(*), PARAMETER :: procedureN = 'cuda_get_device_properties'

    TYPE(cudaDeviceProp_80_t)                :: prop_80
    INTEGER(C_INT)                           :: c_device, c_status
    INTEGER                                  :: runtime_version

    c_device = INT( device, C_INT )


#if defined(_HAS_CUDA)

    CALL cuda_runtime_get_version ( runtime_version )

    SELECT CASE( runtime_version / 1000 )
    CASE( :8 )
       c_status = cudaGetDeviceProperties_80( prop_80 , c_device )
       IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))&
            //". "//trim(cuda_get_error_string(c_status)),&
            __LINE__,__FILE__)

       !vw for the moment we print the info
       WRITE(*,*) 'CUDA: device properties'
       WRITE(*,*) 'CUDA:    name                          ', TRIM( remove_c_null_char ( prop_80%name ) )
       WRITE(*,*) 'CUDA:    totalGlobalMem                ', prop_80%totalGlobalMem
       WRITE(*,*) 'CUDA:    sharedMemPerBlock             ', prop_80%sharedMemPerBlock
       WRITE(*,*) 'CUDA:    regsPerBlock                  ', prop_80%regsPerBlock
       WRITE(*,*) 'CUDA:    warpSize                      ', prop_80%warpSize
       WRITE(*,*) 'CUDA:    memPitch                      ', prop_80%memPitch
       WRITE(*,*) 'CUDA:    maxThreadsPerBlock            ', prop_80%maxThreadsPerBlock
       WRITE(*,*) 'CUDA:    maxThreadsDim                 ', prop_80%maxThreadsDim
       WRITE(*,*) 'CUDA:    maxGridSize                   ', prop_80%maxGridSize
       WRITE(*,*) 'CUDA:    clockRate                     ', prop_80%clockRate
       WRITE(*,*) 'CUDA:    totalConstMem                 ', prop_80%totalConstMem
       WRITE(*,*) 'CUDA:    major                         ', prop_80%major
       WRITE(*,*) 'CUDA:    minor                         ', prop_80%minor
       WRITE(*,*) 'CUDA:    deviceOverlap                 ', prop_80%deviceOverlap
       WRITE(*,*) 'CUDA:    multiProcessorCount           ', prop_80%multiProcessorCount
       WRITE(*,*) 'CUDA:    kernelExecTimeoutEnabled      ', prop_80%kernelExecTimeoutEnabled == 1
       WRITE(*,*) 'CUDA:    integrated                    ', prop_80%integrated == 1
       WRITE(*,*) 'CUDA:    canMapHostMemory              ', prop_80%canMapHostMemory == 1
       SELECT CASE( prop_80%computeMode )
       CASE( cudaComputeModeDefault )
          WRITE(*,*) 'CUDA:    computeMode                   ','cudaComputeModeDefault'
       CASE( cudaComputeModeExclusive )
          WRITE(*,*) 'CUDA:    computeMode                   ','cudaComputeModeExclusive'
       CASE( cudaComputeModeProhibited )
          WRITE(*,*) 'CUDA:    computeMode                   ','cudaComputeModeProhibited'
       CASE( cudaComputeModeExclusiveProcess )
          WRITE(*,*) 'CUDA:    computeMode                   ','cudaComputeModeExclusiveProcess'
       CASE DEFAULT
          WRITE(*,*) 'CUDA:    computeMode                   ','unknown'
       END SELECT
       WRITE(*,*) 'CUDA:    concurrentKernels             ', prop_80%concurrentKernels
       WRITE(*,*) 'CUDA:    unifiedAddressing             ', prop_80%unifiedAddressing == 1
       WRITE(*,*) 'CUDA:    memoryClockRate               ', prop_80%memoryClockRate
       WRITE(*,*) 'CUDA:    memoryBusWidth                ', prop_80%memoryBusWidth
       WRITE(*,*) 'CUDA:    l2CacheSize                   ', prop_80%l2CacheSize
       WRITE(*,*) 'CUDA:    managedMemory                 ', prop_80%managedMemory == 1
       WRITE(*,*) 'CUDA:    isMultiGpuBoard               ', prop_80%isMultiGpuBoard == 1
       WRITE(*,*) 'CUDA:    multiGpuBoardGroupID          ', prop_80%multiGpuBoardGroupID
       WRITE(*,*) 'CUDA:    concurrentManagedAccess       ', prop_80%concurrentManagedAccess == 1
    CASE DEFAULT

       WRITE(*,*) 'CUDA: device properties (not available for runtime version '// int2str( runtime_version ) // ')'

    END SELECT

#endif

  END SUBROUTINE cuda_get_device_properties


  SUBROUTINE cuda_get_device_count ( n_device )
    INTEGER, INTENT(OUT)                     :: n_device

    CHARACTER(*), PARAMETER :: procedureN = 'cuda_get_device_count'

    INTEGER(C_INT)                           :: c_n_device, c_status

    n_device = 0

#if defined(_HAS_CUDA)

    c_status = cudaGetDeviceCount( c_n_device )
    IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))&
         //". "//trim(cuda_get_error_string(c_status)),&
         __LINE__,__FILE__)
    n_device = INT( c_n_device )

#endif

  END SUBROUTINE cuda_get_device_count


  SUBROUTINE cuda_device_reset ( )

    CHARACTER(*), PARAMETER :: procedureN = 'cuda_device_reset'

    INTEGER(C_INT)                           :: c_status

#if defined(_HAS_CUDA)

    c_status = cudaDeviceReset( )
    IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))&
         //". "//trim(cuda_get_error_string(c_status)),&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE cuda_device_reset


  SUBROUTINE cuda_device_synchronize ( )

    CHARACTER(*), PARAMETER :: procedureN = 'cuda_device_synchronize'

    INTEGER(C_INT)                           :: c_status

#if defined(_HAS_CUDA)

    c_status = cudaDeviceSynchronize( )
    IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))&
         //". "//trim(cuda_get_error_string(c_status)),&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE cuda_device_synchronize


  SUBROUTINE cuda_device_get_stream_priority_range ( leastPriority, greatestPriority )
    INTEGER, INTENT(OUT)                     :: leastPriority, &
                                                greatestPriority

    CHARACTER(*), PARAMETER :: &
      procedureN = 'cuda_device_get_stream_priority_range'

    INTEGER(C_INT)                           :: c_greatestPriority, &
                                                c_leastPriority, c_status

#if defined(_HAS_CUDA)

    c_status = cudaDeviceGetStreamPriorityRange( c_leastPriority, c_greatestPriority )
    IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))&
         //". "//trim(cuda_get_error_string(c_status)),&
         __LINE__,__FILE__)

    leastPriority = INT( c_leastPriority )
    greatestPriority = INT( c_greatestPriority )

#else

    leastPriority = 0
    greatestPriority = 1

#endif

  END SUBROUTINE cuda_device_get_stream_priority_range


  SUBROUTINE cuda_get_device ( device )
    INTEGER, INTENT(OUT)                     :: device

    CHARACTER(*), PARAMETER                  :: procedureN = 'cuda_get_device'

    INTEGER(C_INT)                           :: c_device, c_status

    device = cuda_device_null

#if defined(_HAS_CUDA)

    c_status = cudaGetDevice( c_device )
    IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))&
         //". "//trim(cuda_get_error_string(c_status)),&
         __LINE__,__FILE__)

    device = INT( c_device )

#endif

  END SUBROUTINE cuda_get_device


  SUBROUTINE cuda_set_device ( device )
    INTEGER, INTENT(IN)                      :: device

    CHARACTER(*), PARAMETER                  :: procedureN = 'cuda_set_device'

    INTEGER(C_INT)                           :: c_device, c_status

#if defined(_HAS_CUDA)

    CALL cuda_check_device ( device, procedureN )

    c_device = INT( device, C_INT )

    c_status = cudaSetDevice( c_device )
    IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))&
         //". "//trim(cuda_get_error_string(c_status)),&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE cuda_set_device


  SUBROUTINE cuda_alloc_bytes ( mem, n_bytes, device )
    TYPE(cuda_memory_t)                      :: mem
    INTEGER(int_8), INTENT(IN)               :: n_bytes
    INTEGER, INTENT(IN)                      :: device

    CHARACTER(*), PARAMETER :: procedureN = 'cuda_alloc_bytes'

    INTEGER(C_SIZE_T)                        :: c_n_bytes, c_free, c_total
    INTEGER(KIND(cudaSuccess))               :: c_status, c_status_info

#if defined(_HAS_CUDA)

    IF( mem%init ) CALL stopgm(procedureN,'memory already allocated',&
         __LINE__,__FILE__)

    c_n_bytes = INT( n_bytes, C_SIZE_T )

    CALL cuda_set_device ( device )

    c_status = cudaMalloc( mem%ptr, c_n_bytes )
    IF( c_status /= cudaSuccess ) THEN
       c_status_info = cudaMemGetInfo ( c_free, c_total )
       CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))//&
            ' trying to allocate '//TRIM(int2str( n_bytes ))//' bytes'&
            //'. '//trim(cuda_get_error_string(c_status))//'. '// &
            'Memory free=' //TRIM(int2str( INT( c_free, int_8 ) ))//', total='//TRIM(int2str( INT( c_total, int_8 ) ))//&
            ', alloc total='//TRIM(int2str(cuda_mem_alloc_total)),&
            __LINE__,__FILE__)
    ENDIF

    cuda_mem_alloc_total = cuda_mem_alloc_total + c_n_bytes

    mem%device = device
    mem%n_bytes = n_bytes
    mem%init = .TRUE.

#else

    CALL stopgm(procedureN,"no cuda available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE cuda_alloc_bytes


  SUBROUTINE cuda_dealloc ( mem )
    TYPE(cuda_memory_t)                      :: mem

    CHARACTER(*), PARAMETER                  :: procedureN = 'cuda_dealloc'

    INTEGER(KIND(cudaSuccess))               :: c_status

#if defined(_HAS_CUDA)

    IF( .NOT.mem%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    CALL cuda_set_device ( mem%device )
    c_status = cudaFree( mem%ptr )
    IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))&
         //". "//trim(cuda_get_error_string(c_status)),&
         __LINE__,__FILE__)

    cuda_mem_alloc_total = cuda_mem_alloc_total - mem%n_bytes

    mem%ptr = C_NULL_PTR
    mem%device = cuda_device_null
    mem%n_bytes = 0_int_8
    mem%init = .FALSE.

#else

    CALL stopgm(procedureN,"no cuda available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE cuda_dealloc


  SUBROUTINE cuda_alloc_pitch_bytes ( mem, n_bytes_row, n_col, n_bytes_pitch, device )
    TYPE(cuda_memory_t)                      :: mem
    INTEGER(int_8), INTENT(IN)               :: n_bytes_row
    INTEGER, INTENT(IN)                      :: n_col
    INTEGER(int_8), INTENT(OUT)              :: n_bytes_pitch
    INTEGER, INTENT(IN)                      :: device

    CHARACTER(*), PARAMETER :: procedureN = 'cuda_alloc_pitch_bytes'

    INTEGER(C_SIZE_T)                        :: c_n_bytes_pitch, &
                                                c_n_bytes_row, c_n_col, c_free, c_total
    INTEGER(KIND(cudaSuccess))               :: c_status, c_status_info

    c_n_bytes_pitch = 0_int_8

#if defined(_HAS_CUDA)

    IF( mem%init ) CALL stopgm(procedureN,'memory already allocated',&
         __LINE__,__FILE__)

    c_n_bytes_row = INT( n_bytes_row, C_SIZE_T )
    c_n_col = INT( n_col, C_SIZE_T )

    CALL cuda_set_device ( device )
    c_status = cudaMallocPitch ( mem%ptr, c_n_bytes_pitch, c_n_bytes_row, c_n_col )
    IF( c_status /= cudaSuccess ) THEN
       c_status_info = cudaMemGetInfo ( c_free, c_total )
       CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))// &
            ' trying to allocate '//TRIM(int2str( n_bytes_row * INT(n_col,int_8) ))//' bytes'&
            //'. '//trim(cuda_get_error_string(c_status))//'. '// &
            'Memory free=' //TRIM(int2str(INT(c_free,int_8)))//', total='//TRIM(int2str(INT(c_total,int_8)))//&
            ', alloc total='//TRIM(int2str(cuda_mem_alloc_total)),&
            __LINE__,__FILE__)
    ENDIF

    cuda_mem_alloc_total = cuda_mem_alloc_total + c_n_bytes_pitch * INT( n_col, int_8 )

    n_bytes_pitch = INT( c_n_bytes_pitch, int_8 )

    mem%device = device
    mem%n_bytes = INT( n_col, int_8 ) * n_bytes_pitch
    mem%init = .TRUE.

#else

    CALL stopgm(procedureN,"no cuda available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE cuda_alloc_pitch_bytes


  SUBROUTINE cuda_memcpy_device_to_device ( a_devc, b_devc, n_bytes )
    TYPE(cuda_memory_t)                      :: a_devc, b_devc
    INTEGER(int_8), INTENT(IN)               :: n_bytes

    CHARACTER(*), PARAMETER :: procedureN = 'cuda_memcpy_device_to_device'

    INTEGER(C_SIZE_T)                        :: c_n_bytes
    INTEGER(KIND(cudaSuccess))               :: c_status

#if defined(_HAS_CUDA)

    IF( .NOT.a_devc%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT.b_devc%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( n_bytes > 0_int_8 ) THEN

       c_n_bytes = INT( n_bytes, C_SIZE_T )

       c_status = cudaMemcpy( b_devc%ptr, a_devc%ptr, c_n_bytes, cudaMemcpyDeviceToDevice )
       IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))&
            //". "//trim(cuda_get_error_string(c_status)),&
            __LINE__,__FILE__)

    ENDIF
#else

    CALL stopgm(procedureN,"no cuda available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE cuda_memcpy_device_to_device


  SUBROUTINE cuda_stream_create ( stream, device )
    TYPE(cuda_stream_t)                      :: stream
    INTEGER, INTENT(IN)                      :: device

    CHARACTER(*), PARAMETER :: procedureN = 'cuda_stream_create'

    INTEGER(KIND(cudaSuccess))               :: c_status

#if defined(_HAS_CUDA)

    IF( stream%init ) CALL stopgm(procedureN,'stream already initialized',&
         __LINE__,__FILE__)

    CALL cuda_set_device ( device )
    c_status = cudaStreamCreate ( stream%s )
    IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))&
         //". "//trim(cuda_get_error_string(c_status)),&
         __LINE__,__FILE__)

    stream%device = device
    stream%init = .TRUE.

#else

    CALL stopgm(procedureN,"no cuda available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE cuda_stream_create


  SUBROUTINE cuda_stream_create_with_priority ( stream, priority, device )
    TYPE(cuda_stream_t)                      :: stream
    INTEGER, INTENT(IN)                      :: priority, device

    CHARACTER(*), PARAMETER :: procedureN = 'cuda_stream_create_with_priority'

    INTEGER(C_INT)                           :: c_priority
    INTEGER(KIND(cudaSuccess))               :: c_status

#if defined(_HAS_CUDA)

    IF( stream%init ) CALL stopgm(procedureN,'stream already initialized',&
         __LINE__,__FILE__)

    c_priority = INT( priority, C_INT )

    CALL cuda_set_device ( device )
    c_status = cudaStreamCreateWithPriority ( stream%s, cudaStreamNonBlocking, c_priority )
    IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))&
         //". "//trim(cuda_get_error_string(c_status)),&
         __LINE__,__FILE__)

    stream%device = device
    stream%init = .TRUE.

#else

    CALL stopgm(procedureN,"no cuda available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE cuda_stream_create_with_priority


  SUBROUTINE cuda_stream_destroy ( stream )
    TYPE(cuda_stream_t)                      :: stream

    CHARACTER(*), PARAMETER :: procedureN = 'cuda_stream_destroy'

    INTEGER(KIND(cudaSuccess))               :: c_status

#if defined(_HAS_CUDA)

    IF( .NOT. stream%init ) CALL stopgm(procedureN,'stream not initialized',&
         __LINE__,__FILE__)

    c_status = cudaStreamDestroy ( stream%s )
    IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))&
         //". "//trim(cuda_get_error_string(c_status)),&
         __LINE__,__FILE__)

    stream%device = cuda_device_null
    stream%init = .FALSE.

#else

    CALL stopgm(procedureN,"no cuda available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE cuda_stream_destroy


  SUBROUTINE cuda_stream_synchronize ( stream )
    TYPE(cuda_stream_t)                      :: stream

    CHARACTER(*), PARAMETER :: procedureN = 'cuda_stream_synchronize'

    INTEGER(KIND(cudaSuccess))               :: c_status

#if defined(_HAS_CUDA)

    IF( .NOT. stream%init ) CALL stopgm(procedureN,'stream not initialized',&
         __LINE__,__FILE__)

    c_status = cudaStreamSynchronize ( stream%s )
    IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))&
         //". "//trim(cuda_get_error_string(c_status)),&
         __LINE__,__FILE__)

#else

    CALL stopgm(procedureN,"no cuda available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE cuda_stream_synchronize


  SUBROUTINE cuda_stream_query ( stream, ready )
    TYPE(cuda_stream_t)                      :: stream
    LOGICAL, INTENT(OUT)                     :: ready

    CHARACTER(*), PARAMETER :: procedureN = 'cuda_stream_query'

    INTEGER(KIND(cudaSuccess))               :: c_status

    ready = .FALSE.

#if defined(_HAS_CUDA)

    IF( .NOT. stream%init ) CALL stopgm(procedureN,'stream not initialized',&
         __LINE__,__FILE__)

    c_status = cudaStreamQuery( stream%s ) 
    IF( c_status /= cudaSuccess .AND. c_status /= cudaErrorNotReady ) THEN
       CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))&
            //". "//trim(cuda_get_error_string(c_status)),&
            __LINE__,__FILE__)
    ENDIF

    ready = c_status == cudaSuccess

#else

    CALL stopgm(procedureN,"no cuda available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE cuda_stream_query


  SUBROUTINE cuda_stream_wait_event ( stream, event )
    TYPE(cuda_stream_t)                      :: stream
    TYPE(cuda_event_t)                       :: event

    CHARACTER(*), PARAMETER :: procedureN = 'cuda_stream_wait_event'

    INTEGER(KIND(cudaSuccess))               :: c_status

#if defined(_HAS_CUDA)

    IF( .NOT. stream%init ) CALL stopgm(procedureN,'stream not initialized',&
         __LINE__,__FILE__)

    IF( .NOT. event%init ) CALL stopgm(procedureN,'event not initialized',&
         __LINE__,__FILE__)

    c_status = cudaStreamWaitEvent( stream%s, event%e, 0_C_INT )    
    IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))&
         //". "//trim(cuda_get_error_string(c_status)),&
         __LINE__,__FILE__)

#else

    CALL stopgm(procedureN,"no cuda available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE cuda_stream_wait_event


  SUBROUTINE cuda_mem_zero_bytes ( mem, n_bytes )
    TYPE(cuda_memory_t), INTENT(INOUT)       :: mem
    INTEGER(int_8), INTENT(IN)               :: n_bytes

    CHARACTER(*), PARAMETER :: procedureN = 'cuda_mem_zero_bytes'

    INTEGER(C_SIZE_T)                        :: c_n_bytes
    INTEGER(KIND(cudaSuccess))               :: c_status

#if defined(_HAS_CUDA)

    IF( .NOT.mem%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    c_n_bytes = INT( n_bytes, C_SIZE_T )

    CALL cuda_set_device ( mem%device )
    c_status = cudaMemset( mem%ptr, 0_C_INT, c_n_bytes )
    IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))&
         //". "//trim(cuda_get_error_string(c_status)),&
         __LINE__,__FILE__)

#else

    CALL stopgm(procedureN,"no cuda available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE cuda_mem_zero_bytes


  SUBROUTINE cuda_event_create ( event )
    TYPE(cuda_event_t)                       :: event

    CHARACTER(*), PARAMETER :: procedureN = 'cuda_event_create'

    INTEGER(KIND(cudaSuccess))               :: c_status

#if defined(_HAS_CUDA)

    IF( event%init ) CALL stopgm(procedureN,'event already created',&
         __LINE__,__FILE__)

    c_status = cudaEventCreate( event%e )
    IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))&
         //". "//trim(cuda_get_error_string(c_status)),&
         __LINE__,__FILE__)

    event%init = .TRUE.

#else

    CALL stopgm(procedureN,"no cuda available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE cuda_event_create


  SUBROUTINE cuda_event_create_with_no_timing ( event )
    TYPE(cuda_event_t)                       :: event

    CHARACTER(*), PARAMETER :: procedureN = 'cuda_event_create_with_no_timing'

    INTEGER(KIND(cudaSuccess))               :: c_status

#if defined(_HAS_CUDA)

    IF( event%init ) CALL stopgm(procedureN,'event already created',&
         __LINE__,__FILE__)

    c_status = cudaEventCreateWithFlags ( event%e, cudaEventDisableTiming )
    IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))&
         //". "//trim(cuda_get_error_string(c_status)),&
         __LINE__,__FILE__)

    event%init = .TRUE.

#else

    CALL stopgm(procedureN,"no cuda available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE cuda_event_create_with_no_timing


  SUBROUTINE cuda_event_destroy ( event )
    TYPE(cuda_event_t)                       :: event

    CHARACTER(*), PARAMETER :: procedureN = 'cuda_event_destroy'

    INTEGER(KIND(cudaSuccess))               :: c_status

#if defined(_HAS_CUDA)

    IF( .NOT.event%init ) CALL stopgm(procedureN,'event not created',&
         __LINE__,__FILE__)

    c_status = cudaEventDestroy( event%e )
    IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))&
         //". "//trim(cuda_get_error_string(c_status)),&
         __LINE__,__FILE__)

    event%init = .FALSE.

#else

    CALL stopgm(procedureN,"no cuda available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE cuda_event_destroy


  SUBROUTINE cuda_event_synchronize ( event )
    TYPE(cuda_event_t)                       :: event

    CHARACTER(*), PARAMETER :: procedureN = 'cuda_event_synchronize'

    INTEGER(KIND(cudaSuccess))               :: c_status

#if defined(_HAS_CUDA)

    IF( .NOT.event%init ) CALL stopgm(procedureN,'event not created',&
         __LINE__,__FILE__)

    c_status = cudaEventSynchronize( event%e )
    IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))&
         //". "//trim(cuda_get_error_string(c_status)),&
         __LINE__,__FILE__)

#else

    CALL stopgm(procedureN,"no cuda available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE cuda_event_synchronize


  SUBROUTINE cuda_event_elapsed_time ( event_start, event_stop, time )
    TYPE(cuda_event_t)                       :: event_start, event_stop
    REAL(real_8), INTENT(OUT)                :: time

    CHARACTER(*), PARAMETER :: procedureN = 'cuda_event_elapsed_time'

    INTEGER(KIND(cudaSuccess))               :: c_status
    REAL(C_FLOAT)                            :: c_time

    time = 0.0_real_8

#if defined(_HAS_CUDA)

    IF( .NOT.event_start%init .OR. .NOT.event_stop%init ) CALL stopgm(procedureN,'event not created',&
         __LINE__,__FILE__)

    c_status = cudaEventElapsedTime( c_time, event_start%e, event_stop%e )
    IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))&
         //". "//trim(cuda_get_error_string(c_status)),&
         __LINE__,__FILE__)

    time = REAL( c_time, real_8 ) / 1000.0_real_8 ! conver ms to s

#else

    CALL stopgm(procedureN,"no cuda available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE cuda_event_elapsed_time


  SUBROUTINE cuda_event_record( event, stream )
    TYPE(cuda_event_t)                       :: event
    TYPE(cuda_stream_t)                      :: stream

    CHARACTER(*), PARAMETER :: procedureN = 'cuda_event_record'

    INTEGER(KIND(cudaSuccess))               :: c_status

#if defined(_HAS_CUDA)

    IF( .NOT.event%init ) CALL stopgm(procedureN,'event not created',&
         __LINE__,__FILE__)

    IF( .NOT.stream%init ) CALL stopgm(procedureN,'stream not created',&
         __LINE__,__FILE__)


    c_status = cudaEventRecord( event%e, stream%s )
    IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))&
         //". "//trim(cuda_get_error_string(c_status)),&
         __LINE__,__FILE__)

#else

    CALL stopgm(procedureN,"no cuda available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE cuda_event_record

  SUBROUTINE cuda_device_get_pci_bus_id( device, pciBusId )
    INTEGER, INTENT(IN) :: device
    CHARACTER(*), INTENT(OUT) :: pciBusId

    CHARACTER(*), PARAMETER :: procedureN = 'cuda_device_get_pci_bus_id'

    INTEGER :: i
    CHARACTER(C_CHAR), DIMENSION(:), ALLOCATABLE :: c_pciBusId
    INTEGER(KIND(cudaSuccess))                   :: c_status

    pciBusId = ""

#if defined(_HAS_CUDA)

    ALLOCATE(c_pciBusId(LEN(pciBusId)+1))

    CALL cuda_set_device ( device )

    c_status = cudaDeviceGetPCIBusId( c_pciBusId, SIZE( c_pciBusId ), device )
    IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))&
         //". "//trim(cuda_get_error_string(c_status)),&
         __LINE__,__FILE__)

    DO i = 1, SIZE( c_pciBusId ) - 1
       IF( c_pciBusId(i) == C_NULL_CHAR ) EXIT
       pciBusId(i:i) = c_pciBusId(i)
    ENDDO

    DEALLOCATE(c_pciBusId)

#else

    CALL stopgm(procedureN,"no cuda available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE cuda_device_get_pci_bus_id

  FUNCTION cuda_d_points_to ( mem, shift ) RESULT( reslt )
    TYPE(cuda_memory_t), INTENT(IN)          :: mem
    INTEGER, INTENT(IN)                      :: shift
    TYPE(cuda_memory_t)                      :: reslt

    INTEGER(C_INT)                           :: c_elem_size, c_shift

    c_elem_size = 1_C_INT ! double
    c_shift = INT( shift - 1, C_INT )! F to C indexing

    reslt%ptr = cMemPointsToDbl ( mem%ptr, c_elem_size, c_shift )
    reslt%init = .TRUE.
    reslt%n_bytes = 0_int_8
    reslt%device = mem%device

  END FUNCTION cuda_d_points_to

  FUNCTION cuda_z_points_to ( mem, shift ) RESULT( reslt )
    TYPE(cuda_memory_t), INTENT(IN)          :: mem
    INTEGER, INTENT(IN)                      :: shift
    TYPE(cuda_memory_t)                      :: reslt

    INTEGER(C_INT)                           :: c_elem_size, c_shift

    c_elem_size = 2_C_INT ! double complex
    c_shift = INT( shift - 1, C_INT )! F to C indexing

    reslt%ptr = cMemPointsToDbl ( mem%ptr, c_elem_size, c_shift )
    reslt%init = .TRUE.
    reslt%n_bytes = 0_int_8
    reslt%device = mem%device

  END FUNCTION cuda_z_points_to

  FUNCTION cuda_get_address ( mem ) RESULT( reslt )
    TYPE(cuda_memory_t), INTENT(IN)          :: mem
    INTEGER(int_8)                           :: reslt

    reslt = INT( cGetMemAddrs ( mem%ptr ), int_8 )

  END FUNCTION cuda_get_address

  !
  ! include file for the interfaces
  !

#include "cuda_memcpy.inc"

#include "cuda_mem_host.inc"

#include "cuda_memcpy2d.inc"

#include "cuda_mem_managed.inc"

END MODULE cuda_utils
