#include "cpmd_global.h"

MODULE nvml_utils

  USE nvml_types, ONLY: nvml_device_t
  USE nvml_interfaces, ONLY: nvmlDeviceGetUUID, nvmlDeviceGetCount, nvmlDeviceGetHandleByIndex, NVML_SUCCESS,&
       nvmlShutdown, nvmlInit, nvmlDeviceGetPciInfo, nvmlPciInfo_t, nvmlDeviceGetName, nvmlDeviceGetBoardId

  USE error_handling,                  ONLY: stopgm
  USE string_utils,                 ONLY: int2str, remove_c_null_char
  USE, INTRINSIC :: ISO_C_BINDING,  ONLY: C_INT, C_SIZE_T, &
       C_PTR,&
       C_DOUBLE, C_FLOAT, &
       C_DOUBLE_COMPLEX,&
       C_CHAR, C_NULL_CHAR, C_NULL_PTR
  
  IMPLICIT NONE

  PRIVATE


  PUBLIC :: nvml_init
  PUBLIC :: nvml_shutdown
  PUBLIC :: nvml_device_get_count
  PUBLIC :: nvml_device_get_handle_by_index
  PUBLIC :: nvml_device_get_uuid
  PUBLIC :: nvml_device_get_pci_info
  PUBLIC :: nvml_device_get_name
  PUBLIC :: nvml_device_get_board_id


CONTAINS


  SUBROUTINE nvml_init ( )

    CHARACTER(*), PARAMETER :: procedureN = 'nvml_init'
    INTEGER(KIND(NVML_SUCCESS)) :: c_status

#if defined(_HAS_CUDA)

    c_status = nvmlInit(  )
    IF( c_status /= NVML_SUCCESS ) CALL stopgm(procedureN,"nvml error: "//TRIM(int2str( INT( c_status ) )),&
         __LINE__,__FILE__)

#else

    CALL stopgm(procedureN,"no cuda available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE nvml_init


  SUBROUTINE nvml_shutdown ( )

    CHARACTER(*), PARAMETER :: procedureN = 'nvml_shutdown'
    INTEGER(KIND(NVML_SUCCESS)) :: c_status

#if defined(_HAS_CUDA)

    c_status = nvmlShutdown (  )
    IF( c_status /= NVML_SUCCESS ) CALL stopgm(procedureN,"nvml error: "//TRIM(int2str( INT( c_status ) )),&
         __LINE__,__FILE__)

#else

    CALL stopgm(procedureN,"no cuda available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE nvml_shutdown


  SUBROUTINE nvml_device_get_uuid ( device, uuid )
    TYPE( nvml_device_t ), INTENT(IN) :: device
    CHARACTER(*), INTENT(OUT) :: uuid

    INTEGER :: length, i
    INTEGER(C_INT) :: c_length
    INTEGER(KIND(NVML_SUCCESS)) :: c_status
    CHARACTER(C_CHAR), DIMENSION(:), ALLOCATABLE :: c_uuid
    CHARACTER(*), PARAMETER :: procedureN = 'nvml_device_get_uuid'


    uuid = ""

#if defined(_HAS_CUDA)

    IF( .NOT. device%init ) CALL stopgm(procedureN,'device not initialized',&
         __LINE__,__FILE__)

    ALLOCATE(c_uuid(LEN(uuid)+1))
    
    c_status = nvmlDeviceGetUUID( device%d, c_uuid, SIZE( c_uuid, KIND=C_INT) )
    IF( c_status /= NVML_SUCCESS ) CALL stopgm(procedureN,"nvml error: "//TRIM(int2str( INT( c_status ) )),&
         __LINE__,__FILE__)

    DO i = 1, SIZE( c_uuid ) - 1
       IF( c_uuid(i) == C_NULL_CHAR ) EXIT
       uuid(i:i) = c_uuid(i)
    ENDDO

    DEALLOCATE(c_uuid)

#else

    CALL stopgm(procedureN,"no cuda available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE nvml_device_get_uuid

  
  SUBROUTINE nvml_device_get_count ( deviceCount )
    INTEGER, INTENT(OUT) :: deviceCount

    INTEGER(C_INT) :: c_deviceCount
    INTEGER(KIND(NVML_SUCCESS)) :: c_status
    CHARACTER(*), PARAMETER :: procedureN = 'nvml_device_get_count'


    deviceCount = 0

#if defined(_HAS_CUDA)

    c_status = nvmlDeviceGetCount( c_deviceCount )
    IF( c_status /= NVML_SUCCESS ) CALL stopgm(procedureN,"nvml error: "//TRIM(int2str( INT( c_status ) )),&
         __LINE__,__FILE__)

    deviceCount = INT( c_deviceCount )

#else

    CALL stopgm(procedureN,"no cuda available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE nvml_device_get_count


  SUBROUTINE nvml_device_get_handle_by_index ( index, device )
    INTEGER, INTENT(IN) :: index
    TYPE( nvml_device_t ), INTENT(OUT) :: device

    INTEGER(KIND(NVML_SUCCESS)) :: c_status
    INTEGER(C_INT) :: c_index
    CHARACTER(*), PARAMETER :: procedureN = 'nvml_device_get_handle_by_index'

#if defined(_HAS_CUDA)

    c_index = INT( index, KIND=C_INT )

    c_status = nvmlDeviceGetHandleByIndex( c_index, device%d )
    IF( c_status /= NVML_SUCCESS ) CALL stopgm(procedureN,"nvml error: "//TRIM(int2str( INT( c_status ) )),&
         __LINE__,__FILE__)

    device%init = .TRUE.

#else

    CALL stopgm(procedureN,"no cuda available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE nvml_device_get_handle_by_index


  SUBROUTINE nvml_device_get_pci_info ( device )
    TYPE( nvml_device_t ), INTENT(IN) :: device

    INTEGER(KIND(NVML_SUCCESS)) :: c_status
    CHARACTER(*), PARAMETER :: procedureN = 'nvml_device_get_pci_info'
    TYPE( nvmlPciInfo_t ) :: pci

#if defined(_HAS_CUDA)

    IF( .NOT. device%init ) CALL stopgm(procedureN,'device not initialized',&
         __LINE__,__FILE__)

    c_status = nvmlDeviceGetPciInfo( device%d, pci )
    IF( c_status /= NVML_SUCCESS ) CALL stopgm(procedureN,"nvml error: "//TRIM(int2str( INT( c_status ) )),&
         __LINE__,__FILE__)

    !vw for the moment we print the info
    WRITE(*,*) 'NVML:    busId                         ', TRIM( remove_c_null_char ( pci%busId ) )
    !WRITE(*,*) 'NVML:    domain                       ', pci%domain
    !WRITE(*,*) 'NVML:    bus                          ', pci%bus
    !WRITE(*,*) 'NVML:    device                       ', pci%device
    !WRITE(*,*) 'NVML: pciDeviceId    ', pci%pciDeviceId
    !WRITE(*,*) 'NVML: pciSubSystemId ', pci%pciSubSystemId

#else

    CALL stopgm(procedureN,"no cuda available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE nvml_device_get_pci_info


  SUBROUTINE nvml_device_get_name ( device, name )
    TYPE( nvml_device_t ), INTENT(IN) :: device
    CHARACTER(*), INTENT(OUT) :: name

    CHARACTER(C_CHAR), DIMENSION(:), ALLOCATABLE :: c_name
    INTEGER :: i
    INTEGER(KIND(NVML_SUCCESS)) :: c_status
    CHARACTER(*), PARAMETER :: procedureN = 'nvml_device_get_name'

    name = ""

#if defined(_HAS_CUDA)

    IF( .NOT. device%init ) CALL stopgm(procedureN,'device not initialized',&
         __LINE__,__FILE__)

    ALLOCATE(c_name(LEN(name)+1))

    c_status = nvmlDeviceGetName( device%d, c_name, SIZE( c_name, KIND=C_INT ) )
    IF( c_status /= NVML_SUCCESS ) CALL stopgm(procedureN,"nvml error: "//TRIM(int2str( INT( c_status ) )),&
         __LINE__,__FILE__)

    DO i = 1, SIZE( c_name ) - 1
       IF( c_name(i) == C_NULL_CHAR ) EXIT
       name(i:i) = c_name(i)
    ENDDO

    DEALLOCATE(c_name)

#else

    CALL stopgm(procedureN,"no cuda available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE nvml_device_get_name


  SUBROUTINE nvml_device_get_board_id ( device, boardId )
    TYPE( nvml_device_t ), INTENT(IN) :: device
    INTEGER, INTENT(OUT) :: boardId

    INTEGER(KIND(NVML_SUCCESS)) :: c_status
    INTEGER(C_INT) :: c_boardId
    CHARACTER(*), PARAMETER :: procedureN = 'nvml_device_get_board_id'

    boardId = -1

#if defined(_HAS_CUDA)

    IF( .NOT. device%init ) CALL stopgm(procedureN,'device not initialized',&
         __LINE__,__FILE__)

    c_status = nvmlDeviceGetBoardId( device%d, c_boardId )
    IF( c_status /= NVML_SUCCESS ) CALL stopgm(procedureN,"nvml error: "//TRIM(int2str( INT( c_status ) )),&
         __LINE__,__FILE__)

    boardId = INT( c_boardId )

#else

    CALL stopgm(procedureN,"no cuda available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE nvml_device_get_board_id


END MODULE nvml_utils
