MODULE thread_view_utils

  USE error_handling,                  ONLY: stopgm
  USE thread_view_types,               ONLY: thread_view_t

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: thread_views_init
  PUBLIC :: thread_views_finalize

CONTAINS

  SUBROUTINE thread_views_init ( n_devices, n_streams_per_device, thread_views )
    INTEGER, INTENT(IN)                      :: n_devices, &
                                                n_streams_per_device
    TYPE(thread_view_t), ALLOCATABLE, &
      DIMENSION(:), INTENT(INOUT)            :: thread_views

    CHARACTER(*), PARAMETER :: procedureN = 'thread_views_init'

    INTEGER                                  :: device_idx, ierr, &
                                                max_n_threads, stream_idx, &
                                                thread_id

    max_n_threads = n_devices * n_streams_per_device

    IF( ALLOCATED( thread_views ) ) CALL thread_views_finalize ( thread_views )

    ALLOCATE( thread_views( 0:max_n_threads-1 ), STAT=ierr )
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', __LINE__,__FILE__)

    thread_id = 0
    DO device_idx = 1, n_devices
       DO stream_idx = 1, n_streams_per_device
          thread_views( thread_id )%device_idx = device_idx
          thread_views( thread_id )%stream_idx = stream_idx
          thread_views( thread_id )%host_buff_ptr = ( device_idx - 1 ) * n_streams_per_device + stream_idx
          thread_views( thread_id )%init = .TRUE.
          thread_id = thread_id + 1
       ENDDO
    ENDDO

  END SUBROUTINE thread_views_init


  SUBROUTINE thread_views_finalize ( thread_views )
    TYPE(thread_view_t), ALLOCATABLE, &
      DIMENSION(:), INTENT(INOUT)            :: thread_views

    CHARACTER(*), PARAMETER :: procedureN = 'thread_views_finalize'

    INTEGER                                  :: ierr

    IF( ALLOCATED( thread_views ) ) THEN
       DEALLOCATE( thread_views, STAT=ierr )
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', __LINE__,__FILE__)
    ENDIF

  END SUBROUTINE thread_views_finalize

END MODULE thread_view_utils
