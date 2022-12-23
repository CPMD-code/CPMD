MODULE thread_view_types

  USE error_handling,                  ONLY: stopgm

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: thread_view_get

  TYPE, PUBLIC :: thread_view_t
     LOGICAL :: init = .FALSE.

     INTEGER :: device_idx = HUGE(0)
     INTEGER :: stream_idx = HUGE(0)
     INTEGER :: host_buff_ptr = HUGE(0)     
  END TYPE thread_view_t

CONTAINS

  SUBROUTINE thread_view_get ( this, stream_idx, device_idx )
    TYPE(thread_view_t), INTENT(IN)          :: this
    INTEGER, INTENT(OUT), OPTIONAL           :: stream_idx, device_idx

    CHARACTER(*), PARAMETER                  :: procedureN = 'thread_view_get'

    IF( .NOT. this%init ) CALL stopgm(procedureN,"this not initialized",&
         __LINE__,__FILE__)

    IF( PRESENT( stream_idx ) ) stream_idx = this%stream_idx - 1 !vw FIX that -1       
    IF( PRESENT( device_idx ) ) device_idx = this%device_idx

  END SUBROUTINE thread_view_get

END MODULE thread_view_types
