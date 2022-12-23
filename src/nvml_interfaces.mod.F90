MODULE nvml_interfaces

  USE, INTRINSIC :: ISO_C_BINDING,  ONLY: C_INT, C_SIZE_T, &
       C_PTR,&
       C_DOUBLE, C_FLOAT, &
       C_DOUBLE_COMPLEX,&
       C_CHAR, C_NULL_CHAR, C_NULL_PTR


  IMPLICIT NONE


  PRIVATE


  PUBLIC :: nvmlInit
  PUBLIC :: nvmlShutdown
  PUBLIC :: nvmlDeviceGetUUID
  PUBLIC :: nvmlDeviceGetCount
  PUBLIC :: nvmlDeviceGetHandleByIndex
  PUBLIC :: nvmlDeviceGetPciInfo
  PUBLIC :: nvmlDeviceGetName
  PUBLIC :: nvmlDeviceGetBoardId
  !PUBLIC :: 
  !PUBLIC :: 


  ! nvmlReturn_t
  PUBLIC :: NVML_SUCCESS


  TYPE, BIND( C ), PUBLIC :: nvmlDevice_t
     TYPE( C_PTR ) :: s
  END TYPE nvmlDevice_t

  INTEGER, PARAMETER, PUBLIC :: NVML_DEVICE_PCI_BUS_ID_BUFFER_SIZE = 16

  TYPE, BIND(C), PUBLIC :: nvmlPciInfo_t
     character(c_char) :: busId(NVML_DEVICE_PCI_BUS_ID_BUFFER_SIZE) !< The tuple domain:bus:device.function PCI identifier (&amp; NULL terminator)
     integer(c_int) :: domain                                       !< The PCI domain on which the device's bus resides, 0 to 0xffff
     integer(c_int) :: bus                                          !< The bus on which the device resides, 0 to 0xff
     integer(c_int) :: device                                       !< The device's id on the bus, 0 to 31
     integer(c_int) :: pciDeviceId                                  !< The combined 16-bit device id and 16-bit vendor id
     integer(c_int) :: pciSubSystemId                               !< The 32-bit Sub System Device ID
     integer(c_int) :: reserved0, reserved1, reserved2, reserved3   !< NVIDIA reserved for internal use only
  END TYPE nvmlPciInfo_t

  ENUM, BIND(c) !:: nvmlReturn_t
     ENUMERATOR :: NVML_SUCCESS = 0                   !< The operation was successful
     ENUMERATOR :: NVML_ERROR_UNINITIALIZED = 1       !< NVML was not first initialized with nvmlInit()
     ENUMERATOR :: NVML_ERROR_INVALID_ARGUMENT = 2    !< A supplied argument is invalid
     ENUMERATOR :: NVML_ERROR_NOT_SUPPORTED = 3       !< The requested operation is not available on target device
     ENUMERATOR :: NVML_ERROR_NO_PERMISSION = 4       !< The current user does not have permission for operation
     ENUMERATOR :: NVML_ERROR_ALREADY_INITIALIZED = 5 !< Deprecated: Multiple initializations are now allowed through ref counting
     ENUMERATOR :: NVML_ERROR_NOT_FOUND = 6           !< A query to find an object was unsuccessful
     ENUMERATOR :: NVML_ERROR_INSUFFICIENT_SIZE = 7   !< An input argument is not large enough
     ENUMERATOR :: NVML_ERROR_INSUFFICIENT_POWER = 8  !< A device's external power cables are not properly attached
     ENUMERATOR :: NVML_ERROR_DRIVER_NOT_LOADED = 9   !< NVIDIA driver is not loaded
     ENUMERATOR :: NVML_ERROR_TIMEOUT = 10            !< User provided timeout passed
     ENUMERATOR :: NVML_ERROR_IRQ_ISSUE = 11          !< NVIDIA Kernel detected an interrupt issue with a GPU
     ENUMERATOR :: NVML_ERROR_LIBRARY_NOT_FOUND = 12  !< NVML Shared Library couldn't be found or loaded
     ENUMERATOR :: NVML_ERROR_FUNCTION_NOT_FOUND = 13 !< Local version of NVML doesn't implement this function
     ENUMERATOR :: NVML_ERROR_CORRUPTED_INFOROM = 14  !< infoROM is corrupted
     ENUMERATOR :: NVML_ERROR_GPU_IS_LOST = 15        !< The GPU has fallen off the bus or has otherwise become inaccessible
     ENUMERATOR :: NVML_ERROR_RESET_REQUIRED = 16     !< The GPU requires a reset before it can be used again
     ENUMERATOR :: NVML_ERROR_OPERATING_SYSTEM = 17   !< The GPU control device has been blocked by the operating system/cgroups
     ENUMERATOR :: NVML_ERROR_LIB_RM_VERSION_MISMATCH = 18   !< RM detects a driver/library version mismatch
     ENUMERATOR :: NVML_ERROR_IN_USE = 19             !< An operation cannot be performed because the GPU is currently in use
     ENUMERATOR :: NVML_ERROR_NO_DATA = 20            !< No data
     ENUMERATOR :: NVML_ERROR_UNKNOWN = 999           !< An internal driver error occurred
  END ENUM

  INTERFACE

     ! nvmlReturn_t nvmlInit ( void )
     INTEGER( KIND ( NVML_SUCCESS ) ) FUNCTION nvmlInit ( ) BIND( C, NAME="nvmlInit" )
       IMPORT :: NVML_SUCCESS
       IMPLICIT NONE
     END FUNCTION nvmlInit

     ! nvmlReturn_t nvmlShutdown ( void ) 
     INTEGER( KIND ( NVML_SUCCESS ) ) FUNCTION nvmlShutdown ( ) BIND( C, NAME="nvmlShutdown" )
       IMPORT :: NVML_SUCCESS
       IMPLICIT NONE
     END FUNCTION nvmlShutdown

     ! nvmlReturn_t nvmlDeviceGetUUID ( nvmlDevice_t device, char* uuid, unsigned int  length ) 
     INTEGER( KIND ( NVML_SUCCESS ) ) FUNCTION nvmlDeviceGetUUID( device, uuid, length ) BIND( C, NAME="nvmlDeviceGetUUID" )
       IMPORT :: C_CHAR, nvmlDevice_t, C_INT, NVML_SUCCESS
       IMPLICIT NONE
       TYPE( nvmlDevice_t ), VALUE :: device
       CHARACTER(C_CHAR) :: uuid(*)
       INTEGER( C_INT ), VALUE :: length
     END FUNCTION nvmlDeviceGetUUID

     ! nvmlReturn_t nvmlDeviceGetCount ( unsigned int* deviceCount )
     INTEGER( KIND ( NVML_SUCCESS ) ) FUNCTION nvmlDeviceGetCount( deviceCount ) BIND( C, NAME="nvmlDeviceGetCount" )
       IMPORT :: C_INT, NVML_SUCCESS
       IMPLICIT NONE
       INTEGER( C_INT ) :: deviceCount
     END FUNCTION nvmlDeviceGetCount

     ! nvmlReturn_t nvmlDeviceGetHandleByIndex ( unsigned int  index, nvmlDevice_t* device )
     INTEGER( KIND ( NVML_SUCCESS ) ) FUNCTION nvmlDeviceGetHandleByIndex( index, device ) &
          & BIND( C, NAME="nvmlDeviceGetHandleByIndex" )
       IMPORT :: C_INT, nvmlDevice_t, NVML_SUCCESS
       IMPLICIT NONE
       INTEGER( C_INT ), VALUE :: index
       TYPE( nvmlDevice_t ) :: device
     END FUNCTION nvmlDeviceGetHandleByIndex

     ! nvmlReturn_t nvmlDeviceGetPciInfo ( nvmlDevice_t device, nvmlPciInfo_t* pci )
     INTEGER( KIND ( NVML_SUCCESS ) ) FUNCTION nvmlDeviceGetPciInfo( device, pci ) BIND( C, NAME="nvmlDeviceGetPciInfo" )
       IMPORT :: nvmlDevice_t, NVML_SUCCESS, nvmlPciInfo_t
       IMPLICIT NONE
       TYPE( nvmlDevice_t ), VALUE :: device
       TYPE( nvmlPciInfo_t ) :: pci
     END FUNCTION nvmlDeviceGetPciInfo

     ! nvmlReturn_t nvmlDeviceGetName ( nvmlDevice_t device, char* name, unsigned int  length )
     INTEGER( KIND ( NVML_SUCCESS ) ) FUNCTION nvmlDeviceGetName( device, name, length ) BIND( C, NAME="nvmlDeviceGetName" )
       IMPORT :: nvmlDevice_t, NVML_SUCCESS, C_INT, C_CHAR
       IMPLICIT NONE
       TYPE( nvmlDevice_t ), VALUE :: device
       CHARACTER( C_CHAR ) :: name(*)
       INTEGER( C_INT ), VALUE :: length
     END FUNCTION nvmlDeviceGetName

     ! nvmlReturn_t nvmlDeviceGetBoardId ( nvmlDevice_t device, unsigned int* boardId )
     INTEGER( KIND ( NVML_SUCCESS ) ) FUNCTION nvmlDeviceGetBoardId( device, boardId ) BIND( C, NAME="nvmlDeviceGetBoardId" )
       IMPORT :: nvmlDevice_t, NVML_SUCCESS, C_INT
       IMPLICIT NONE
       TYPE( nvmlDevice_t ), VALUE :: device
       INTEGER( C_INT ) :: boardId
     END FUNCTION nvmlDeviceGetBoardId

     ! nvmlReturn_t nvmlDeviceOnSameBoard ( nvmlDevice_t device1, nvmlDevice_t device2, int* onSameBoard )

     ! nvmlReturn_t nvmlDeviceGetComputeMode ( nvmlDevice_t device, nvmlComputeMode_t* mode )
     !INTEGER( KIND ( NVML_SUCCESS ) ) FUNCTION nvmlDeviceGetComputeMode( device, mode ) BIND( C, NAME="nvmlDeviceGetComputeMode" )
     !END FUNCTION nvmlDeviceGetComputeMode

     ! nvmlReturn_t nvmlDeviceGetCpuAffinity ( nvmlDevice_t device, unsigned int  cpuSetSize, unsignedlong* cpuSet )


  END INTERFACE

END MODULE nvml_interfaces
