MODULE cp_curho_types

  USE cuda_types,                      ONLY: cuda_memory_t
  USE error_handling,                  ONLY: stopgm

  IMPLICIT NONE

  PRIVATE

  TYPE, PUBLIC :: cp_curho_stream_t
     LOGICAL :: init = .FALSE.

     TYPE( cuda_memory_t ) :: rho_d

  END TYPE cp_curho_stream_t


  TYPE, PUBLIC :: cp_curho_device_t
     LOGICAL :: init = .FALSE.

     !device id
     INTEGER :: id = HUGE(0)

     TYPE( cp_curho_stream_t ), DIMENSION(:), ALLOCATABLE :: streams

  END TYPE cp_curho_device_t


  TYPE, PUBLIC :: cp_curho_t
     LOGICAL :: init = .FALSE.

     TYPE( cp_curho_device_t ), DIMENSION(:), ALLOCATABLE :: devices

  END TYPE cp_curho_t

  TYPE( cp_curho_t ), SAVE, TARGET, PUBLIC :: cp_curho

  PUBLIC :: cp_curho_stream_get_ptrs

CONTAINS


  SUBROUTINE cp_curho_stream_get_ptrs ( this, device_idx, stream_idx, rho_d )
    TYPE(cp_curho_t), INTENT(IN), TARGET     :: this
    INTEGER, INTENT(IN)                      :: device_idx, stream_idx
    TYPE(cuda_memory_t), INTENT(OUT), &
      OPTIONAL, POINTER                      :: rho_d

    CHARACTER(*), PARAMETER :: procedureN = 'cp_curho_stream_get_ptrs'

    IF( .NOT. this%init ) CALL stopgm(procedureN,"this not initialized",&
         __LINE__,__FILE__)

    IF(  device_idx > UBOUND( this%devices, 1 ) .OR. &
         device_idx < LBOUND( this%devices, 1 ) ) THEN
       CALL stopgm(procedureN,"device_idx is wrong",&
            __LINE__,__FILE__)
    ENDIF
    IF(  stream_idx > UBOUND( this%devices( device_idx )%streams, 1 ) .OR. &
         stream_idx < LBOUND( this%devices( device_idx )%streams, 1 ) ) THEN
       CALL stopgm(procedureN,"stream_idx is wrong",&
            __LINE__,__FILE__)
    ENDIF

    IF(PRESENT( rho_d )) rho_d => this%devices(device_idx)%streams(stream_idx)%rho_d

  END SUBROUTINE cp_curho_stream_get_ptrs


END MODULE cp_curho_types
