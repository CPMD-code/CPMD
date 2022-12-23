MODULE cp_cuwfn_types

  USE cuda_types,                      ONLY: cuda_memory_t
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  PRIVATE

  TYPE, PUBLIC :: cp_cuwfn_device_t
     LOGICAL :: init = .FALSE.

     !device id
     INTEGER :: id = HUGE(0)

     TYPE( cuda_memory_t ) :: c0_d

     !TYPE( cp_curho_stream_t ), DIMENSION(:), ALLOCATABLE :: streams

  END TYPE cp_cuwfn_device_t

  TYPE, PUBLIC :: cp_cuwfn_t
     LOGICAL :: init = .FALSE.

     LOGICAL :: ready = .FALSE.
     REAL(real_8) :: chksum = 0.0_real_8

     TYPE( cp_cuwfn_device_t ), DIMENSION(:), ALLOCATABLE :: devices

  END TYPE cp_cuwfn_t

  TYPE( cp_cuwfn_t ), SAVE, TARGET, PUBLIC :: cp_cuwfn

  PUBLIC :: cp_cuwfn_get
  PUBLIC :: cp_cuwfn_put
  PUBLIC :: cp_cuwfn_device_get_ptrs
  PUBLIC :: cp_cuwfn_is_init


CONTAINS


  FUNCTION cp_cuwfn_is_init ( this ) RESULT( reslt )
    TYPE(cp_cuwfn_t), INTENT(IN)             :: this
    LOGICAL                                  :: reslt

    reslt = this%init

  END FUNCTION cp_cuwfn_is_init


  SUBROUTINE cp_cuwfn_get ( this, ready, chksum )
    TYPE(cp_cuwfn_t), INTENT(IN), TARGET     :: this
    LOGICAL, INTENT(OUT), OPTIONAL           :: ready
    REAL(real_8), INTENT(OUT), OPTIONAL      :: chksum

    CHARACTER(*), PARAMETER                  :: procedureN = 'cp_cuwfn_get'

    IF(PRESENT( ready ) ) ready = this%ready
    IF(PRESENT( chksum ) ) chksum = this%chksum
  END SUBROUTINE cp_cuwfn_get


  SUBROUTINE cp_cuwfn_put ( this, ready, chksum )
    TYPE(cp_cuwfn_t), INTENT(INOUT), TARGET  :: this
    LOGICAL, INTENT(IN), OPTIONAL            :: ready
    REAL(real_8), INTENT(IN), OPTIONAL       :: chksum

    CHARACTER(*), PARAMETER                  :: procedureN = 'cp_cuwfn_put'

    IF(PRESENT( ready ) ) this%ready = ready
    IF(PRESENT( chksum ) ) this%chksum = chksum
  END SUBROUTINE cp_cuwfn_put


  SUBROUTINE cp_cuwfn_device_get_ptrs ( this, device_idx, c0_d )
    TYPE(cp_cuwfn_t), INTENT(IN), TARGET     :: this
    INTEGER, INTENT(IN)                      :: device_idx
    TYPE(cuda_memory_t), INTENT(OUT), &
      OPTIONAL, POINTER                      :: c0_d

    CHARACTER(*), PARAMETER :: procedureN = 'cp_cuwfn_device_get_ptrs'

    IF( .NOT. this%init ) CALL stopgm(procedureN,"this not initialized",&
         __LINE__,__FILE__)

    IF(  device_idx > UBOUND( this%devices, 1 ) .OR. &
         device_idx < LBOUND( this%devices, 1 ) ) THEN
       CALL stopgm(procedureN,"device_idx is wrong",&
            __LINE__,__FILE__)
    ENDIF

    IF(PRESENT( c0_d )) c0_d => this%devices(device_idx)%c0_d

  END SUBROUTINE cp_cuwfn_device_get_ptrs

END MODULE cp_cuwfn_types
