MODULE cuda_interfaces

  USE, INTRINSIC :: ISO_C_BINDING,  ONLY: C_INT, C_SIZE_T, &
       C_PTR,&
       C_DOUBLE, C_FLOAT, &
       C_DOUBLE_COMPLEX,&
       C_CHAR, C_NULL_CHAR, C_NULL_PTR

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: cudaMalloc
  PUBLIC :: cudaFree
  PUBLIC :: cudaMallocHost
  PUBLIC :: cudaHostAlloc
  PUBLIC :: cudaFreeHost
  PUBLIC :: cudaMallocPitch
  PUBLIC :: cudaHostRegister
  PUBLIC :: cudaHostUnregister
  PUBLIC :: cudaMallocManaged
  PUBLIC :: cudaMemcpy
  PUBLIC :: cudaMemcpy2D
  PUBLIC :: cudaMemcpyAsync
  PUBLIC :: cudaMemcpy2DAsync
  PUBLIC :: cudaMemset
  PUBLIC :: cudaStreamCreate
  PUBLIC :: cudaStreamCreateWithPriority
  PUBLIC :: cudaStreamDestroy
  PUBLIC :: cudaStreamSynchronize
  PUBLIC :: cudaStreamQuery
  PUBLIC :: cudaStreamWaitEvent
  PUBLIC :: cudaStreamAttachMemAsync
  PUBLIC :: cudaEventCreate
  PUBLIC :: cudaEventCreateWithFlags
  PUBLIC :: cudaEventRecord
  PUBLIC :: cudaEventQuery
  PUBLIC :: cudaEventSynchronize
  PUBLIC :: cudaEventElapsedTime
  PUBLIC :: cudaEventDestroy
  PUBLIC :: cudaDriverGetVersion
  PUBLIC :: cudaRuntimeGetVersion
  PUBLIC :: cudaGetErrorString
  PUBLIC :: cudaGetDeviceProperties_80
  PUBLIC :: cudaGetDeviceCount
  PUBLIC :: cudaGetDevice
  PUBLIC :: cudaSetDevice
  PUBLIC :: cudaDeviceReset
  PUBLIC :: cudaDeviceSynchronize
  PUBLIC :: cudaDeviceGetStreamPriorityRange
  PUBLIC :: cudaDeviceGetPCIBusId
  PUBLIC :: cudaPointerGetAttributes
  PUBLIC :: cudaMemPrefetchAsync
  PUBLIC :: cudaMemAdvise
  PUBLIC :: cudaMemGetInfo
  PUBLIC :: cMemPointsToDbl
  PUBLIC :: cGetMemAddrs

  !cudaError
  PUBLIC :: cudaSuccess
  PUBLIC :: cudaErrorNotReady
  PUBLIC :: cudaErrorHostMemoryAlreadyRegistered

  !cudaMemcpyKind
  PUBLIC :: cudaMemcpyHostToDevice
  PUBLIC :: cudaMemcpyDeviceToHost
  PUBLIC :: cudaMemcpyDeviceToDevice

  !cudaComputeMode
  PUBLIC :: cudaComputeModeDefault
  PUBLIC :: cudaComputeModeExclusive
  PUBLIC :: cudaComputeModeProhibited
  PUBLIC :: cudaComputeModeExclusiveProcess

  !cudaMemoryType                                                                                                                      
  PUBLIC :: cudaMemoryTypeHost
  PUBLIC :: cudaMemoryTypeDevice

  !cudaMemoryAdvise
  PUBLIC :: cudaMemAdviseSetReadMostly
  PUBLIC :: cudaMemAdviseUnsetReadMostly
  PUBLIC :: cudaMemAdviseSetPreferredLocation
  PUBLIC :: cudaMemAdviseUnsetPreferredLocation
  PUBLIC :: cudaMemAdviseSetAccessedBy
  PUBLIC :: cudaMemAdviseUnsetAccessedBy

  !
  ! vw the BOZ arguments to INT are valide from f2003
  !
  INTEGER( C_INT ), PUBLIC, PARAMETER :: cudaHostAllocDefault       = INT( z'00', C_INT ) !Default page-locked allocation flag
  INTEGER( C_INT ), PUBLIC, PARAMETER :: cudaHostAllocPortable      = INT( z'01', C_INT ) !Pinned memory accessible by all CUDA contexts
  INTEGER( C_INT ), PUBLIC, PARAMETER :: cudaHostAllocMapped        = INT( z'02', C_INT ) !Map allocation into device space
  INTEGER( C_INT ), PUBLIC, PARAMETER :: cudaHostAllocWriteCombined  = INT( z'04', C_INT ) !Write-combined memory

  INTEGER( C_INT ), PUBLIC, PARAMETER :: cudaHostRegisterDefault = 0
  INTEGER( C_INT ), PUBLIC, PARAMETER :: cudaHostRegisterPortable = 1
  INTEGER( C_INT ), PUBLIC, PARAMETER :: cudaHostRegisterMapped = 2
  INTEGER( C_INT ), PUBLIC, PARAMETER :: cudaHostRegisterIoMemory = 4

  INTEGER( C_INT ), PUBLIC, PARAMETER :: cudaStreamDefault = INT( z'00', C_INT )
  INTEGER( C_INT ), PUBLIC, PARAMETER :: cudaStreamNonBlocking = INT( z'01', C_INT )

  INTEGER( C_INT ), PUBLIC, PARAMETER :: cudaEventDefault       = INT( z'00', C_INT ) !Default event flag
  INTEGER( C_INT ), PUBLIC, PARAMETER :: cudaEventBlockingSync  = INT( z'01', C_INT ) !Event uses blocking synchronization
  INTEGER( C_INT ), PUBLIC, PARAMETER :: cudaEventDisableTiming = INT( z'02', C_INT ) !Event will not record timing data
  INTEGER( C_INT ), PUBLIC, PARAMETER :: cudaEventInterprocess  = INT( z'04', C_INT ) !Event is suitable for interprocess use. cudaEventDisableTiming must be set  

  INTEGER( C_INT ), PUBLIC, PARAMETER :: cudaMemAttachGlobal  = INT( z'01', C_INT ) !Memory can be accessed by any stream on any device
  INTEGER( C_INT ), PUBLIC, PARAMETER :: cudaMemAttachHost    = INT( z'02', C_INT ) !Memory cannot be accessed by any stream on any device
  INTEGER( C_INT ), PUBLIC, PARAMETER :: cudaMemAttachSingle  = INT( z'04', C_INT ) !Memory can only be accessed by a single stream on the associated device

  INTEGER( C_INT ), PUBLIC, PARAMETER :: cudaCpuDeviceId      = -1 !Device id that represents the CPU
  INTEGER( C_INT ), PUBLIC, PARAMETER :: cudaInvalidDeviceId  = -2 !Device id that represents an invalid device


  ENUM, BIND( C ) !:: cudaError
     ENUMERATOR :: cudaSuccess = 0
     ENUMERATOR :: cudaErrorMissingConfiguration = 1
     ENUMERATOR :: cudaErrorMemoryAllocation = 2
     ENUMERATOR :: cudaErrorInitializationError = 3
     ENUMERATOR :: cudaErrorLaunchFailure = 4
     ENUMERATOR :: cudaErrorPriorLaunchFailure = 5
     ENUMERATOR :: cudaErrorLaunchTimeout = 6
     ENUMERATOR :: cudaErrorLaunchOutOfResources = 7
     ENUMERATOR :: cudaErrorInvalidDeviceFunction = 8
     ENUMERATOR :: cudaErrorInvalidConfiguration = 9
     ENUMERATOR :: cudaErrorInvalidDevice = 10
     ENUMERATOR :: cudaErrorInvalidValue = 11
     ENUMERATOR :: cudaErrorInvalidPitchValue = 12
     ENUMERATOR :: cudaErrorInvalidSymbol = 13
     ENUMERATOR :: cudaErrorMapBufferObjectFailed = 14
     ENUMERATOR :: cudaErrorUnmapBufferObjectFailed = 15
     ENUMERATOR :: cudaErrorInvalidHostPointer = 16
     ENUMERATOR :: cudaErrorInvalidDevicePointer = 17
     ENUMERATOR :: cudaErrorInvalidTexture = 18
     ENUMERATOR :: cudaErrorInvalidTextureBinding = 19
     ENUMERATOR :: cudaErrorInvalidChannelDescriptor = 20
     ENUMERATOR :: cudaErrorInvalidMemcpyDirection = 21
     ENUMERATOR :: cudaErrorAddressOfConstant = 22
     ENUMERATOR :: cudaErrorTextureFetchFailed = 23
     ENUMERATOR :: cudaErrorTextureNotBound = 24
     ENUMERATOR :: cudaErrorSynchronizationError = 25
     ENUMERATOR :: cudaErrorInvalidFilterSetting = 26
     ENUMERATOR :: cudaErrorInvalidNormSetting = 27
     ENUMERATOR :: cudaErrorMixedDeviceExecution = 28
     ENUMERATOR :: cudaErrorCudartUnloading = 29
     ENUMERATOR :: cudaErrorUnknown = 30
     ENUMERATOR :: cudaErrorNotYetImplemented = 31
     ENUMERATOR :: cudaErrorMemoryValueTooLarge = 32
     ENUMERATOR :: cudaErrorInvalidResourceHandle = 33
     ENUMERATOR :: cudaErrorNotReady = 34
     ENUMERATOR :: cudaErrorInsufficientDriver = 35
     ENUMERATOR :: cudaErrorSetOnActiveProcess = 36
     ENUMERATOR :: cudaErrorInvalidSurface = 37
     ENUMERATOR :: cudaErrorNoDevice = 38
     ENUMERATOR :: cudaErrorECCUncorrectable = 39
     ENUMERATOR :: cudaErrorSharedObjectSymbolNotFound = 40
     ENUMERATOR :: cudaErrorSharedObjectInitFailed = 41
     ENUMERATOR :: cudaErrorUnsupportedLimit = 42
     ENUMERATOR :: cudaErrorDuplicateVariableName = 43
     ENUMERATOR :: cudaErrorDuplicateTextureName = 44
     ENUMERATOR :: cudaErrorDuplicateSurfaceName = 45
     ENUMERATOR :: cudaErrorDevicesUnavailable = 46
     ENUMERATOR :: cudaErrorInvalidKernelImage           =     47
     ENUMERATOR :: cudaErrorNoKernelImageForDevice       =     48
     ENUMERATOR :: cudaErrorIncompatibleDriverContext    =     49
     ENUMERATOR :: cudaErrorPeerAccessAlreadyEnabled     =     50
     ENUMERATOR :: cudaErrorPeerAccessNotEnabled         =     51
     ENUMERATOR :: cudaErrorDeviceAlreadyInUse           =     54
     ENUMERATOR :: cudaErrorProfilerDisabled             =     55
     ENUMERATOR :: cudaErrorProfilerNotInitialized       =     56
     ENUMERATOR :: cudaErrorProfilerAlreadyStarted       =     57
     ENUMERATOR :: cudaErrorProfilerAlreadyStopped       =     58
     ENUMERATOR :: cudaErrorAssert                       =     59
     ENUMERATOR :: cudaErrorTooManyPeers                 =     60  
     ENUMERATOR :: cudaErrorHostMemoryAlreadyRegistered  =     61
     ENUMERATOR :: cudaErrorHostMemoryNotRegistered      =     62
     ENUMERATOR :: cudaErrorOperatingSystem              =     63
     ENUMERATOR :: cudaErrorPeerAccessUnsupported        =     64
     ENUMERATOR :: cudaErrorLaunchMaxDepthExceeded       =     65
     ENUMERATOR :: cudaErrorLaunchFileScopedTex          =     66
     ENUMERATOR :: cudaErrorLaunchFileScopedSurf         =     67
     ENUMERATOR :: cudaErrorSyncDepthExceeded            =     68
     ENUMERATOR :: cudaErrorLaunchPendingCountExceeded   =     69
     ENUMERATOR :: cudaErrorNotPermitted                 =     70
     ENUMERATOR :: cudaErrorNotSupported                 =     71
     ENUMERATOR :: cudaErrorHardwareStackError           =     72
     ENUMERATOR :: cudaErrorIllegalInstruction           =     73
     ENUMERATOR :: cudaErrorMisalignedAddress            =     74
     ENUMERATOR :: cudaErrorInvalidAddressSpace          =     75
     ENUMERATOR :: cudaErrorInvalidPc                    =     76
     ENUMERATOR :: cudaErrorIllegalAddress               =     77
     ENUMERATOR :: cudaErrorInvalidPtx                   =     78
     ENUMERATOR :: cudaErrorInvalidGraphicsContext       =     79
     ENUMERATOR :: cudaErrorNvlinkUncorrectable          =     80
     ENUMERATOR :: cudaErrorStartupFailure = 127
     ENUMERATOR :: cudaErrorApiFailureBase = 10000
  END ENUM

  ENUM, BIND( C ) !:: cudaMemcpyKind
     ENUMERATOR :: cudaMemcpyHostToHost = 0
     ENUMERATOR :: cudaMemcpyHostToDevice = 1
     ENUMERATOR :: cudaMemcpyDeviceToHost = 2
     ENUMERATOR :: cudaMemcpyDeviceToDevice = 3
  END ENUM


  ENUM, BIND( C ) !:: cudaComputeMode
     ENUMERATOR :: cudaComputeModeDefault          = 0
     ENUMERATOR :: cudaComputeModeExclusive        = 1
     ENUMERATOR :: cudaComputeModeProhibited       = 2
     ENUMERATOR :: cudaComputeModeExclusiveProcess = 3
  END ENUM

  ENUM, BIND(c) !:: cudaMemoryType
     ENUMERATOR :: cudaMemoryTypeHost = 1
     ENUMERATOR :: cudaMemoryTypeDevice = 2
  END ENUM


  ENUM, BIND(c) !:: cudaMemoryAdvise
     ENUMERATOR :: cudaMemAdviseSetReadMostly          = 1 !Data will mostly be read and only occassionally be written to
     ENUMERATOR :: cudaMemAdviseUnsetReadMostly        = 2 !Undo the effect of cudaMemAdviseSetReadMostly
     ENUMERATOR :: cudaMemAdviseSetPreferredLocation   = 3 !Set the preferred location for the data as the specified device
     ENUMERATOR :: cudaMemAdviseUnsetPreferredLocation = 4 !Clear the preferred location for the data
     ENUMERATOR :: cudaMemAdviseSetAccessedBy          = 5 !Data will be accessed by the specified device, so prevent page faults as much as possible */
     ENUMERATOR :: cudaMemAdviseUnsetAccessedBy        = 6 !Let the Unified Memory subsystem decide on the page faulting policy for the specified device
  END ENUM



  TYPE, BIND( C ), PUBLIC :: cudaStream_t
     TYPE( C_PTR ) :: s
  END TYPE cudaStream_t

  TYPE, BIND( C ), PUBLIC :: cudaEvent_t
     TYPE( C_PTR ) :: event
  END TYPE cudaEvent_t

  TYPE, BIND(C), PUBLIC :: cudaPointerAttributes_t
     INTEGER( KIND( cudaMemoryTypeHost ) ) :: memoryType
     INTEGER(c_int) :: device
     TYPE(c_ptr) :: devicePointer
     TYPE(c_ptr) :: hostPointer
     INTEGER(c_int) :: isManaged
  END TYPE cudaPointerAttributes_t

  TYPE, BIND(C), PUBLIC :: cudaDeviceProp_80_t
     CHARACTER(c_char) :: name(256)
     INTEGER(c_size_t) :: totalGlobalMem
     INTEGER(c_size_t) :: sharedMemPerBlock
     INTEGER(c_int) :: regsPerBlock
     INTEGER(c_int) :: warpSize
     INTEGER(c_size_t) :: memPitch
     INTEGER(c_int) :: maxThreadsPerBlock
     INTEGER(c_int) :: maxThreadsDim(3)
     INTEGER(c_int) :: maxGridSize(3)
     INTEGER(c_int) :: clockRate
     INTEGER(c_size_t) :: totalConstMem
     INTEGER(c_int) :: major
     INTEGER(c_int) :: minor
     INTEGER(c_size_t) :: textureAlignment
     INTEGER(c_size_t) :: texturePitchAlignment
     INTEGER(c_int) :: deviceOverlap
     INTEGER(c_int) :: multiProcessorCount
     INTEGER(c_int) :: kernelExecTimeoutEnabled
     INTEGER(c_int) :: integrated
     INTEGER(c_int) :: canMapHostMemory
     INTEGER(c_int) :: computeMode
     INTEGER(c_int) :: maxTexture1D
     INTEGER(c_int) :: maxTexture1DMipmap
     INTEGER(c_int) :: maxTexture1DLinear
     INTEGER(c_int) :: maxTexture2D(2)
     INTEGER(c_int) :: maxTexture2DMipmap(2)
     INTEGER(c_int) :: maxTexture2DLinear(3)
     INTEGER(c_int) :: maxTexture2DGather(2)
     INTEGER(c_int) :: maxTexture3D(3)
     INTEGER(c_int) :: maxTexture3DAlt(3)
     INTEGER(c_int) :: maxTextureCubemap
     INTEGER(c_int) :: maxTexture1DLayered(2)
     INTEGER(c_int) :: maxTexture2DLayered(3)
     INTEGER(c_int) :: maxTextureCubemapLayered(2)
     INTEGER(c_int) :: maxSurface1D
     INTEGER(c_int) :: maxSurface2D(2)
     INTEGER(c_int) :: maxSurface3D(3)
     INTEGER(c_int) :: maxSurface1DLayered(2)
     INTEGER(c_int) :: maxSurface2DLayered(3)
     INTEGER(c_int) :: maxSurfaceCubemap
     INTEGER(c_int) :: maxSurfaceCubemapLayered(2)
     INTEGER(c_size_t) :: surfaceAlignment
     INTEGER(c_int) :: concurrentKernels
     INTEGER(c_int) :: ECCEnabled
     INTEGER(c_int) :: pciBusID
     INTEGER(c_int) :: pciDeviceID
     INTEGER(c_int) :: pciDomainID
     INTEGER(c_int) :: tccDriver
     INTEGER(c_int) :: asyncEngineCount
     INTEGER(c_int) :: unifiedAddressing
     INTEGER(c_int) :: memoryClockRate
     INTEGER(c_int) :: memoryBusWidth
     INTEGER(c_int) :: l2CacheSize
     INTEGER(c_int) :: maxThreadsPerMultiProcessor
     INTEGER(c_int) :: streamPrioritiesSupported
     INTEGER(c_int) :: globalL1CacheSupported
     INTEGER(c_int) :: localL1CacheSupported
     INTEGER(c_size_t) :: sharedMemPerMultiprocessor
     INTEGER(c_int) :: regsPerMultiprocessor
     INTEGER(c_int) :: managedMemory
     INTEGER(c_int) :: isMultiGpuBoard
     INTEGER(c_int) :: multiGpuBoardGroupID
     !cuda 8
     INTEGER(c_int) :: hostNativeAtomicSupported
     INTEGER(c_int) :: singleToDoublePrecisionPerfRatio
     INTEGER(c_int) :: pageableMemoryAccess
     INTEGER(c_int) :: concurrentManagedAccess
  END TYPE cudaDeviceProp_80_t


  INTERFACE

     INTEGER( KIND( cudaSuccess ) ) FUNCTION cudaMalloc( devPtr,size ) BIND( C, NAME="cudaMalloc" )
       IMPORT :: C_INT, C_SIZE_T, C_PTR, cudaSuccess
       IMPLICIT NONE
       TYPE( C_PTR ) :: devPtr
       INTEGER( C_SIZE_T ), VALUE :: size
     END FUNCTION cudaMalloc

     INTEGER( KIND( cudaSuccess ) ) FUNCTION cudaFree( devPtr ) BIND( C, NAME="cudaFree" )
       IMPORT :: C_INT, C_PTR, cudaSuccess
       IMPLICIT NONE
       TYPE( C_PTR ), VALUE :: devPtr
     END FUNCTION cudaFree

     INTEGER( KIND( cudaSuccess ) ) FUNCTION cudaMallocHost( ptr, size ) BIND( C, NAME="cudaMallocHost" )
       IMPORT :: C_SIZE_T, C_PTR, cudaSuccess
       IMPLICIT NONE
       TYPE( C_PTR ) :: ptr
       INTEGER( C_SIZE_T ), VALUE :: size
     END FUNCTION cudaMallocHost

     !cudaError_t cudaHostAlloc ( void **  pHost, size_t  size, unsigned int  flags ) 
     INTEGER( KIND( cudaSuccess ) ) FUNCTION cudaHostAlloc ( pHost, bytes, flags ) BIND( C, NAME='cudaHostAlloc' )
       IMPORT :: C_SIZE_T, C_PTR, C_INT, cudaSuccess
       IMPLICIT NONE
       TYPE( C_PTR) :: pHost
       INTEGER( C_SIZE_T ), VALUE :: bytes
       INTEGER( C_INT ), VALUE :: flags
     END FUNCTION cudaHostAlloc

     INTEGER( KIND( cudaSuccess ) ) FUNCTION cudaFreeHost( ptr ) BIND( C, NAME="cudaFreeHost" )
       IMPORT :: C_PTR, cudaSuccess
       IMPLICIT NONE
       TYPE( C_PTR ), VALUE :: ptr
     END FUNCTION cudaFreeHost

     !cudaError_t cudaMallocPitch ( void **  devPtr, size_t *  pitch, size_t  width, size_t  height )
     INTEGER( KIND( cudaSuccess ) ) FUNCTION cudaMallocPitch ( devPtr, pitch, width, &
          & height ) BIND( C, NAME='cudaMallocPitch' )
       IMPORT :: C_SIZE_T, C_PTR, cudaSuccess
       IMPLICIT NONE
       TYPE( C_PTR ) :: devPtr
       INTEGER( C_SIZE_T ) :: pitch
       INTEGER( C_SIZE_T ), VALUE :: width
       INTEGER( C_SIZE_T ), VALUE :: height
     END FUNCTION cudaMallocPitch

     !cudaError_t cudaHostRegister(void *ptr, size_t size, unsigned int flags)
     INTEGER( KIND( cudaSuccess ) ) FUNCTION cudaHostRegister ( ptr, size, flags ) BIND( C, NAME='cudaHostRegister' )
       IMPORT :: C_SIZE_T, C_INT, cudaSuccess, C_PTR
       IMPLICIT NONE
       TYPE( C_PTR ), VALUE :: ptr
       INTEGER( C_SIZE_T ), VALUE :: size
       INTEGER( C_INT ), VALUE :: flags
     END FUNCTION cudaHostRegister

     !cudaError_t cudaHostUnregister(void *ptr)
     INTEGER( KIND( cudaSuccess ) ) FUNCTION cudaHostUnregister ( ptr ) BIND( C, NAME='cudaHostUnregister' )
       IMPORT :: cudaSuccess, C_PTR
       IMPLICIT NONE
       TYPE( C_PTR ), VALUE :: ptr
     END FUNCTION cudaHostUnregister

     !     cudaError_t cudaMallocManaged(void **devPtr,
     !                                   size_t size,
     !                                   unsigned int flags=0);
     INTEGER( KIND( cudaSuccess ) ) FUNCTION cudaMallocManaged ( devPtr, size, flags ) BIND( C, NAME='cudaMallocManaged' )
       IMPORT :: C_SIZE_T, C_PTR, C_INT, cudaSuccess
       IMPLICIT NONE
       TYPE( C_PTR ) :: devPtr
       INTEGER( C_SIZE_T ), VALUE :: size
       INTEGER( C_INT ), VALUE :: flags
     END FUNCTION cudaMallocManaged

     !     cudaError_t cudaMemPrefetchAsync(const void *devPtr,
     !                                      size_t count,
     !                                      int dstDevice,
     !                                      cudaStream_t stream);
     INTEGER( KIND( cudaSuccess ) ) FUNCTION cudaMemPrefetchAsync ( devPtr, count, dstDevice, stream ) &
          BIND( C, NAME='cudaMemPrefetchAsync' )
       IMPORT :: C_SIZE_T, C_PTR, C_INT, cudaSuccess, cudaStream_t
       IMPLICIT NONE
       TYPE( C_PTR ), VALUE :: devPtr
       INTEGER( C_SIZE_T ), VALUE :: count
       INTEGER( C_INT ), VALUE :: dstDevice
       TYPE( cudaStream_t ), VALUE :: stream
     END FUNCTION cudaMemPrefetchAsync

     !     cudaError_t cudaMemcpy ( void *        dst,
     !                              const void *  src,
     !                              size_t        count,
     !                              enum cudaMemcpyKind  kind )
     INTEGER( KIND( cudaSuccess ) ) FUNCTION cudaMemcpy( dst, src, count, kind_arg ) BIND( C, name="cudaMemcpy" )
       IMPORT :: C_SIZE_T, C_PTR, cudaSuccess, cudaMemcpyHostToHost
       IMPLICIT NONE
       TYPE( C_PTR ), VALUE :: dst
       TYPE( C_PTR ), VALUE :: src
       INTEGER( C_SIZE_T ), VALUE :: count
       INTEGER( KIND( cudaMemcpyHostToHost ) ), VALUE :: kind_arg
     END FUNCTION cudaMemcpy

     !     cudaError_t cudaMemcpy2D ( void *  dst,
     !                                size_t  dpitch,
     !                                const void *  src,
     !                                size_t  spitch,
     !                                size_t  width,
     !                                size_t  height,
     !                                enum cudaMemcpyKind  kind ) 
     INTEGER( KIND( cudaSuccess ) ) FUNCTION cudaMemcpy2D( dst, dpitch, src, spitch, width, &
          & height, kind_arg ) BIND( C, NAME="cudaMemcpy2D" )
       IMPORT :: C_SIZE_T, C_PTR, cudaSuccess, cudaMemcpyHostToHost
       IMPLICIT NONE
       TYPE( C_PTR ), VALUE :: dst
       INTEGER( C_SIZE_T ), VALUE :: dpitch
       TYPE( C_PTR ), VALUE :: src
       INTEGER( C_SIZE_T ), VALUE :: spitch
       INTEGER( C_SIZE_T ), VALUE :: width
       INTEGER( C_SIZE_T ), VALUE :: height
       INTEGER( KIND( cudaMemcpyHostToHost ) ), VALUE :: kind_arg
     END FUNCTION cudaMemcpy2D

     INTEGER( KIND( cudaSuccess ) ) FUNCTION cudaMemcpyAsync( dst, src, count, kind_arg, stream ) BIND( C, name="cudaMemcpyAsync" )
       IMPORT :: C_SIZE_T, C_PTR, cudaSuccess, cudaMemcpyHostToHost, cudaStream_t
       IMPLICIT NONE
       TYPE( C_PTR ), VALUE :: dst
       TYPE( C_PTR ), VALUE :: src
       INTEGER( C_SIZE_T ), VALUE :: count
       INTEGER( KIND( cudaMemcpyHostToHost ) ), VALUE :: kind_arg
       TYPE( cudaStream_t ), VALUE :: stream
     END FUNCTION cudaMemcpyAsync

     !     cudaError_t cudaMemcpy2DAsync ( void *  dst,
     !                                     size_t  dpitch,
     !                                     const void *  src,
     !                                     size_t  spitch,
     !                                     size_t  width,
     !                                     size_t  height,
     !                                     enum cudaMemcpyKind  kind,
     !                                     cudaStream_t  stream = 0 ) 
     INTEGER( KIND( cudaSuccess ) ) FUNCTION cudaMemcpy2DAsync( dst, dpitch, src, spitch, width, &
          & height, kind_arg, stream ) BIND( C, NAME="cudaMemcpy2DAsync" )
       IMPORT :: C_SIZE_T, C_PTR, cudaSuccess, cudaMemcpyHostToHost, cudaStream_t
       IMPLICIT NONE
       TYPE( C_PTR ), VALUE :: dst
       INTEGER( C_SIZE_T ), VALUE :: dpitch
       TYPE( C_PTR ), VALUE :: src
       INTEGER( C_SIZE_T ), VALUE :: spitch
       INTEGER( C_SIZE_T ), VALUE :: width
       INTEGER( C_SIZE_T ), VALUE :: height
       INTEGER( KIND( cudaMemcpyHostToHost ) ), VALUE :: kind_arg
       TYPE( cudaStream_t ), VALUE :: stream
     END FUNCTION cudaMemcpy2DAsync

     INTEGER( KIND( cudaSuccess ) ) FUNCTION cudaMemset( devPtr, VALUE, count ) BIND( C, NAME="cudaMemset" )
       IMPORT :: C_SIZE_T, C_INT, C_PTR, cudaSuccess
       IMPLICIT NONE
       TYPE( C_PTR ), VALUE :: devPtr
       INTEGER( C_INT ), VALUE :: VALUE
       INTEGER( C_SIZE_T ), VALUE :: count
     END FUNCTION cudaMemset

     INTEGER( KIND( cudaSuccess ) ) FUNCTION cudaStreamCreate( pstream ) BIND( C, NAME="cudaStreamCreate" )
       IMPORT :: cudaSuccess, cudaStream_t
       IMPLICIT NONE
       TYPE( cudaStream_t ) :: pstream
     END FUNCTION cudaStreamCreate

     !cudaError_t cudaStreamCreateWithPriority ( cudaStream_t* pStream, unsigned int  flags, int  priority ) 
     INTEGER( KIND( cudaSuccess ) ) FUNCTION cudaStreamCreateWithPriority( pStream, flags, priority ) &
          BIND( C, NAME="cudaStreamCreateWithPriority" )
       IMPORT :: cudaSuccess, cudaStream_t, C_INT
       IMPLICIT NONE
       TYPE( cudaStream_t ) :: pstream
       INTEGER( C_INT ), VALUE :: flags
       INTEGER( C_INT ), VALUE :: priority
     END FUNCTION cudaStreamCreateWithPriority

     INTEGER( KIND( cudaSuccess ) ) FUNCTION cudaStreamDestroy( stream ) BIND( C, NAME="cudaStreamDestroy" )
       IMPORT :: cudaSuccess, cudaStream_t
       IMPLICIT NONE
       TYPE( cudaStream_t ), VALUE :: stream
     END FUNCTION cudaStreamDestroy

     INTEGER( KIND( cudaSuccess ) ) FUNCTION cudaStreamSynchronize( stream ) BIND( C, NAME="cudaStreamSynchronize" )
       IMPORT :: cudaSuccess, cudaStream_t
       IMPLICIT NONE
       TYPE( cudaStream_t ), VALUE :: stream
     END FUNCTION cudaStreamSynchronize

     INTEGER( KIND( cudaSuccess ) ) FUNCTION cudaStreamQuery( stream ) BIND( C, NAME="cudaStreamQuery" )
       IMPORT :: cudaSuccess, cudaStream_t
       IMPLICIT NONE
       TYPE( cudaStream_t ), VALUE :: stream
     END FUNCTION cudaStreamQuery

     !cudaError_t cudaStreamWaitEvent ( cudaStream_t stream, cudaEvent_t event, unsigned int  flags ) 
     INTEGER( KIND( cudaSuccess ) ) FUNCTION cudaStreamWaitEvent( stream, event, flags ) BIND( C, NAME='cudaStreamWaitEvent' )
       IMPORT :: cudaSuccess, cudaStream_t, C_INT, cudaEvent_t
       IMPLICIT NONE
       TYPE( cudaStream_t ), VALUE :: stream
       TYPE( cudaEvent_t ), VALUE :: event
       INTEGER( C_INT ), VALUE :: flags
     END FUNCTION cudaStreamWaitEvent

     !cudaError_t cudaStreamAttachMemAsync ( cudaStream_t stream, void* devPtr, size_t length = 0, unsigned int  flags = cudaMemAttachSingle ) 
     INTEGER( KIND( cudaSuccess ) ) FUNCTION cudaStreamAttachMemAsync( stream, devPtr, length, flags )
       IMPORT :: cudaSuccess, cudaStream_t, C_INT, C_PTR, C_SIZE_T
       IMPLICIT NONE
       TYPE( cudaStream_t ), VALUE :: stream
       TYPE( C_PTR ), VALUE :: devPtr
       INTEGER( C_SIZE_T ), VALUE :: length
       INTEGER( C_INT ), VALUE :: flags
     END FUNCTION cudaStreamAttachMemAsync

     INTEGER( KIND( cudaSuccess ) ) FUNCTION cudaEventCreate( event ) BIND( C, NAME="cudaEventCreate" )
       IMPORT :: cudaSuccess, cudaevent_t
       IMPLICIT NONE
       TYPE( cudaEvent_t ) :: event
     END FUNCTION cudaEventCreate

     INTEGER( KIND( cudaSuccess ) ) FUNCTION cudaEventDestroy( event ) BIND( C, NAME="cudaEventDestroy" )
       IMPORT :: cudaSuccess, cudaEvent_t
       IMPLICIT NONE
       TYPE( cudaEvent_t ), VALUE :: event
     END FUNCTION cudaEventDestroy

     INTEGER( KIND( cudaSuccess ) ) FUNCTION cudaEventCreateWithFlags( event, flags ) BIND( C, NAME="cudaEventCreateWithFlags" )
       IMPORT :: cudaSuccess, cudaEvent_t, C_INT
       IMPLICIT NONE
       TYPE( cudaEvent_t ) :: event
       INTEGER( C_INT ), VALUE :: flags
     END FUNCTION cudaEventCreateWithFlags

     INTEGER( KIND( cudaSuccess ) ) FUNCTION cudaEventRecord( event, stream ) BIND( C, NAME="cudaEventRecord" )
       IMPORT :: cudaSuccess, cudaStream_t, cudaEvent_t
       IMPLICIT NONE
       TYPE( cudaEvent_t ), VALUE :: event
       TYPE( cudaStream_t ), VALUE :: stream
     END FUNCTION cudaEventRecord

     INTEGER( KIND( cudaSuccess ) ) FUNCTION cudaEventQuery( event ) BIND( C, NAME="cudaEventQuery" )
       IMPORT :: cudaSuccess, cudaEvent_t
       IMPLICIT NONE
       TYPE( cudaEvent_t ), VALUE :: event
     END FUNCTION cudaEventQuery

     INTEGER( KIND( cudaSuccess ) ) FUNCTION cudaEventSynchronize( event ) BIND( C, NAME="cudaEventSynchronize" )
       IMPORT :: cudaSuccess, cudaEvent_t
       IMPLICIT NONE
       TYPE( cudaEvent_t ), VALUE :: event
     END FUNCTION cudaEventSynchronize

     !cudaError_t cudaEventElapsedTime ( float* ms, cudaEvent_t start, cudaEvent_t end ) 
     INTEGER( KIND( cudaSuccess ) ) FUNCTION cudaEventElapsedTime( ms, start, END ) BIND( C, NAME="cudaEventElapsedTime" )
       IMPORT :: cudaSuccess, cudaEvent_t, C_FLOAT
       IMPLICIT NONE
       REAL( C_FLOAT ) :: ms
       TYPE( cudaEvent_t ), VALUE :: start
       TYPE( cudaEvent_t ), VALUE :: END
     END FUNCTION cudaEventElapsedTime

     INTEGER( KIND( cudaSuccess ) ) FUNCTION cudaDriverGetVersion( driverVersion ) BIND( C, NAME="cudaDriverGetVersion" )
       IMPORT :: C_INT, cudaSuccess
       IMPLICIT NONE
       INTEGER( C_INT ) :: driverVersion
     END FUNCTION cudaDriverGetVersion

     INTEGER( KIND( cudaSuccess ) ) FUNCTION cudaRuntimeGetVersion( runtimeVersion ) BIND( C, NAME="cudaRuntimeGetVersion" )
       IMPORT :: C_INT, cudaSuccess
       IMPLICIT NONE
       INTEGER( C_INT ) :: runtimeVersion
     END FUNCTION cudaRuntimeGetVersion

     TYPE( C_PTR ) FUNCTION cudaGetErrorString( error ) BIND( C, NAME='cudaGetErrorString' )
       IMPORT :: C_PTR, cudaSuccess
       !IMPLICIT NONE !xl bug
       INTEGER( KIND( cudaSuccess ) ), VALUE :: error
     END FUNCTION cudaGetErrorString

     INTEGER( KIND( cudaSuccess ) ) FUNCTION cudaGetDeviceProperties_80( prop, device ) BIND( C, NAME='cudaGetDeviceProperties' )
       IMPORT :: C_INT, cudaSuccess, cudaDeviceProp_80_t
       IMPLICIT NONE
       TYPE(cudaDeviceProp_80_t) :: prop
       INTEGER(C_INT), VALUE :: device
     END FUNCTION cudaGetDeviceProperties_80

     INTEGER( KIND( cudaSuccess ) ) FUNCTION cudaGetDeviceCount( count ) BIND( C, NAME='cudaGetDeviceCount' )
       IMPORT :: C_INT, cudaSuccess
       IMPLICIT NONE
       INTEGER(C_INT) :: count
     END FUNCTION cudaGetDeviceCount

     !cudaError_t cudaDeviceReset(void)
     INTEGER( KIND( cudaSuccess ) ) FUNCTION cudaDeviceReset( ) BIND( C, NAME='cudaDeviceReset' )
       IMPORT :: cudaSuccess
       IMPLICIT NONE
     END FUNCTION cudaDeviceReset

     !cudaError_t cudaDeviceSynchronize(void);
     INTEGER( KIND( cudaSuccess ) ) FUNCTION cudaDeviceSynchronize( ) BIND( C, NAME='cudaDeviceSynchronize' )
       IMPORT :: cudaSuccess
       IMPLICIT NONE
     END FUNCTION cudaDeviceSynchronize

     !cudaError_t cudaDeviceGetStreamPriorityRange ( int* leastPriority, int* greatestPriority ) 
     INTEGER( KIND( cudaSuccess ) ) FUNCTION cudaDeviceGetStreamPriorityRange( leastPriority, greatestPriority ) &
          BIND( C, NAME='cudaDeviceGetStreamPriorityRange' )
       IMPORT :: cudaSuccess, C_INT
       IMPLICIT NONE
       INTEGER( C_INT ) :: leastPriority
       INTEGER( C_INT ) :: greatestPriority
     END FUNCTION cudaDeviceGetStreamPriorityRange

     !cudaError_t cudaSetDevice(int device)
     INTEGER( KIND( cudaSuccess ) ) FUNCTION cudaSetDevice( device ) BIND( C, NAME='cudaSetDevice' )
       IMPORT cudaSuccess, C_INT
       IMPLICIT NONE
       INTEGER(C_INT), VALUE :: device
     END FUNCTION cudaSetDevice

     !cudaError_t cudaGetDevice( int *device')
     INTEGER( KIND( cudaSuccess ) ) FUNCTION cudaGetDevice( device ) BIND( C, name='cudaGetDevice' )
       IMPORT cudaSuccess, C_INT
       IMPLICIT NONE
       INTEGER(C_INT) :: device
     END FUNCTION cudaGetDevice

     !cudaError_t cudaPointerGetAttributes ( cudaPointerAttributes* attributes, const void* ptr ) 
     INTEGER( KIND( cudaSuccess ) ) FUNCTION cudaPointerGetAttributes ( attributes, ptr ) &
          BIND( C, NAME='cudaPointerGetAttributes' )
       IMPORT :: cudaSuccess, cudaPointerAttributes_t, C_PTR
       IMPLICIT NONE
       TYPE(cudaPointerAttributes_t) :: attributes
       TYPE(C_PTR), VALUE :: ptr
     END FUNCTION cudaPointerGetAttributes

     !cudaError_t cudaDeviceGetPCIBusId ( char* pciBusId, int  len, int  device ) 
     INTEGER( KIND( cudaSuccess ) ) FUNCTION cudaDeviceGetPCIBusId( pciBusId, len, device ) &
          BIND( C, NAME='cudaDeviceGetPCIBusId' )
       IMPORT :: cudaSuccess, C_CHAR, C_INT
       IMPLICIT NONE
       CHARACTER(C_CHAR) :: pciBusId(*)
       INTEGER(C_INT), VALUE :: len
       INTEGER(C_INT), VALUE :: device
     END FUNCTION cudaDeviceGetPCIBusId

     !cudaError_t cudaMemAdvise(const void *devPtr, size_t count, enum cudaMemoryAdvise advice, int device)
     INTEGER( KIND( cudaSuccess ) ) FUNCTION cudaMemAdvise ( devPtr, count, advice, device ) &
          BIND( C, NAME='cudaMemAdvise' )
       IMPORT :: cudaSuccess, C_PTR, C_INT, C_SIZE_T, cudaMemAdviseSetReadMostly
       IMPLICIT NONE
       TYPE(C_PTR), VALUE :: devPtr
       INTEGER(C_SIZE_T), VALUE :: count
       INTEGER(KIND(cudaMemAdviseSetReadMostly)), VALUE :: advice
       INTEGER(C_INT), VALUE :: device
     END FUNCTION cudaMemAdvise

     !cudaError_t cudaMemGetInfo ( size_t* free, size_t* total )
     INTEGER( KIND( cudaSuccess ) ) FUNCTION cudaMemGetInfo ( free, total ) BIND( C, NAME='cudaMemGetInfo' )
       IMPORT :: cudaSuccess, C_SIZE_T
       IMPLICIT NONE
       INTEGER(C_SIZE_T) :: free, total
     END FUNCTION cudaMemGetInfo


     !cudaError_t cudaMemRangeGetAttribute ( void* data, size_t dataSize, enum cudaMemRangeAttribute attribute, const void* devPtr, size_t count ) 
     !....

     SUBROUTINE CuIdentityDouble( devPtr, m ) BIND( C, NAME='CuIdentityDouble' )
       IMPORT :: C_PTR, C_INT
       TYPE( C_PTR ), VALUE :: devPtr
       INTEGER( C_INT ), VALUE :: m
     END SUBROUTINE CuIdentityDouble

     TYPE( C_PTR ) FUNCTION cMemPointsToDbl( p, sizeElem, shift ) BIND( C, NAME="cMemPointsToDbl" )
       IMPORT :: C_PTR, C_INT
       !IMPLICIT NONE !xl bug
       TYPE( C_PTR ), VALUE :: p
       INTEGER( C_INT ), VALUE :: sizeElem, shift
     END FUNCTION cMemPointsToDbl

     !size_t cGetMemAddrs ( void *p )
     INTEGER( C_SIZE_T ) FUNCTION cGetMemAddrs( p ) BIND( C, NAME='cGetMemAddrs' )
       IMPORT :: C_PTR, C_SIZE_T
       TYPE( C_PTR ), VALUE :: p
     END FUNCTION cGetMemAddrs

  END INTERFACE


END MODULE cuda_interfaces
