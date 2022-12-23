MODULE cp_curho_utils

  USE cp_cuda_types,                   ONLY: cp_cuda_devices_t,&
                                             cp_cuda_env_t
  USE cp_curho_types,                  ONLY: cp_curho_device_t,&
                                             cp_curho_stream_t,&
                                             cp_curho_t
  USE cuda_utils,                      ONLY: cuda_alloc_bytes,&
                                             cuda_dealloc,&
                                             cuda_mem_zero_bytes
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: int_8
  USE sizeof_kinds,                    ONLY: sizeof_real_8
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: cp_curho_init
  PUBLIC :: cp_curho_finalize
  PUBLIC :: cp_curho_alloc_buffers
  PUBLIC :: cp_curho_dealloc_buffers

CONTAINS


  SUBROUTINE cp_curho_init ( cp_cuda_env, cp_cuda_devices_fft, cp_curho )
    TYPE(cp_cuda_env_t), INTENT(IN)          :: cp_cuda_env
    TYPE(cp_cuda_devices_t), INTENT(IN)      :: cp_cuda_devices_fft
    TYPE(cp_curho_t), INTENT(INOUT)          :: cp_curho

    CHARACTER(*), PARAMETER                  :: procedureN = 'cp_curho_init'

    INTEGER                                  :: device_id, i_device, ierr, &
                                                isub

    CALL tiset(procedureN,isub)

    IF( cp_curho%init ) CALL stopgm(procedureN,"cp_curho already initialized",&
         __LINE__,__FILE__)

    ALLOCATE( cp_curho%devices( cp_cuda_env%fft_n_devices_per_task ), STAT=ierr )
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', __LINE__,__FILE__)

    DO i_device = 1, cp_cuda_env%fft_n_devices_per_task
       device_id = cp_cuda_devices_fft%ids( i_device )
       CALL cp_curho_init_device( cp_cuda_env, cp_curho%devices( i_device ), device_id )
    ENDDO

    cp_curho%init = .TRUE.

    CALL tihalt(procedureN,isub)

  END SUBROUTINE cp_curho_init


  SUBROUTINE cp_curho_alloc_buffers ( cp_curho, n )
    TYPE(cp_curho_t), INTENT(INOUT)          :: cp_curho
    INTEGER, INTENT(IN)                      :: n

    CHARACTER(*), PARAMETER :: procedureN = 'cp_curho_alloc_buffers'

    INTEGER                                  :: i_device, isub

    CALL tiset(procedureN,isub)

    IF( .NOT. cp_curho%init ) CALL stopgm(procedureN,"cp_curho not initialized",&
         __LINE__,__FILE__)

    DO i_device = 1, SIZE( cp_curho%devices, 1 )
       CALL cp_curho_device_alloc_buffers ( cp_curho%devices( i_device ), n )
    ENDDO

    CALL tihalt(procedureN,isub)

  END SUBROUTINE cp_curho_alloc_buffers


  SUBROUTINE cp_curho_finalize ( cp_curho )
    TYPE(cp_curho_t), INTENT(INOUT)          :: cp_curho

    CHARACTER(*), PARAMETER :: procedureN = 'cp_curho_finalize'

    INTEGER                                  :: i_device, ierr, isub

    CALL tiset(procedureN,isub)

    IF( .NOT.cp_curho%init ) CALL stopgm(procedureN,"cp_curho not initialized",&
         __LINE__,__FILE__)

    DO i_device = 1, SIZE( cp_curho%devices, 1 )
       CALL cp_curho_finalize_device ( cp_curho%devices( i_device ) )
    ENDDO

    DEALLOCATE( cp_curho%devices, STAT=ierr )
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', __LINE__,__FILE__)

    cp_curho%init = .FALSE.

    CALL tihalt(procedureN,isub)

  END SUBROUTINE cp_curho_finalize


  SUBROUTINE cp_curho_dealloc_buffers ( cp_curho )
    TYPE(cp_curho_t), INTENT(INOUT)          :: cp_curho

    CHARACTER(*), PARAMETER :: procedureN = 'cp_curho_dealloc_buffers'

    INTEGER                                  :: i_device, isub

    CALL tiset(procedureN,isub)

    IF( .NOT.cp_curho%init ) CALL stopgm(procedureN,"cp_curho not initialized",&
         __LINE__,__FILE__)

    DO i_device = 1, SIZE( cp_curho%devices, 1 )
       CALL cp_curho_device_dealloc_buffers ( cp_curho%devices( i_device ) )
    ENDDO

    CALL tihalt(procedureN,isub)

  END SUBROUTINE cp_curho_dealloc_buffers


  SUBROUTINE cp_curho_init_device( cp_cuda_env, device, device_id )
    TYPE(cp_cuda_env_t), INTENT(IN)          :: cp_cuda_env
    TYPE(cp_curho_device_t), INTENT(INOUT)   :: device
    INTEGER, INTENT(IN)                      :: device_id

    CHARACTER(*), PARAMETER :: procedureN = 'cp_curho_init_device'

    INTEGER                                  :: i_stream, ierr

    IF( device%init ) CALL stopgm(procedureN,"device already initialized",&
         __LINE__,__FILE__)

    device%id = device_id

    ALLOCATE( device%streams( 0:cp_cuda_env%fft_n_streams_per_device-1 ), STAT=ierr )
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', __LINE__,__FILE__)

    DO i_stream = 0, cp_cuda_env%fft_n_streams_per_device-1
       CALL cp_curho_init_stream ( device%streams( i_stream ), device_id )
    ENDDO

    device%init = .TRUE.

  END SUBROUTINE cp_curho_init_device


  SUBROUTINE cp_curho_device_alloc_buffers ( device, n )
    TYPE(cp_curho_device_t), INTENT(INOUT)   :: device
    INTEGER, INTENT(IN)                      :: n

    CHARACTER(*), PARAMETER :: procedureN = 'cp_curho_device_alloc_buffers'

    INTEGER                                  :: device_id, i_stream

    IF( .NOT. device%init ) CALL stopgm(procedureN,"device not initialized",&
         __LINE__,__FILE__)

    device_id = device%id

    DO i_stream = LBOUND( device%streams, 1 ), UBOUND( device%streams, 1 )
       CALL cp_curho_stream_alloc_buffers ( device%streams( i_stream ), device_id, n )
    ENDDO

  END SUBROUTINE cp_curho_device_alloc_buffers


  SUBROUTINE cp_curho_finalize_device ( device )
    TYPE(cp_curho_device_t), INTENT(INOUT)   :: device

    CHARACTER(*), PARAMETER :: procedureN = 'cp_curho_finalize_device'

    INTEGER                                  :: i_stream, ierr, isub

    CALL tiset(procedureN,isub)

    IF( .NOT.device%init ) CALL stopgm(procedureN,"device not initialized",&
         __LINE__,__FILE__)

    DO i_stream = LBOUND( device%streams,1 ), UBOUND( device%streams,1 )
       CALL cp_curho_finalize_stream ( device%streams( i_stream ) )
    ENDDO

    DEALLOCATE( device%streams, STAT=ierr )
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', __LINE__,__FILE__)

    device%id = HUGE(0)
    device%init = .FALSE.

    CALL tihalt(procedureN,isub)

  END SUBROUTINE cp_curho_finalize_device


  SUBROUTINE cp_curho_device_dealloc_buffers ( device )
    TYPE(cp_curho_device_t), INTENT(INOUT)   :: device

    CHARACTER(*), PARAMETER :: procedureN = 'cp_curho_device_dealloc_buffers'

    INTEGER                                  :: i_stream, isub

    CALL tiset(procedureN,isub)

    IF( .NOT.device%init ) CALL stopgm(procedureN,"device not initialized",&
         __LINE__,__FILE__)

    DO i_stream = LBOUND( device%streams, 1 ), UBOUND( device%streams, 1 )
       CALL cp_curho_stream_dealloc_buffers ( device%streams( i_stream ) )
    ENDDO

    CALL tihalt(procedureN,isub)

  END SUBROUTINE cp_curho_device_dealloc_buffers


  SUBROUTINE cp_curho_init_stream( stream, device_id )
    TYPE(cp_curho_stream_t), INTENT(INOUT)   :: stream
    INTEGER, INTENT(IN)                      :: device_id

    CHARACTER(*), PARAMETER :: procedureN = 'cp_curho_init_stream'

    INTEGER                                  :: isub

    CALL tiset(procedureN,isub)

    IF( stream%init ) CALL stopgm(procedureN,"stream already initialized",&
         __LINE__,__FILE__)

    !CALL cuda_stream_create ( stream%stream_1, device_id )

    stream%init = .TRUE.

    CALL tihalt(procedureN,isub)

  END SUBROUTINE cp_curho_init_stream


  SUBROUTINE cp_curho_stream_alloc_buffers( stream, device_id, n )
    TYPE(cp_curho_stream_t), INTENT(INOUT)   :: stream
    INTEGER, INTENT(IN)                      :: device_id, n

    CHARACTER(*), PARAMETER :: procedureN = 'cp_curho_stream_alloc_buffers'

    INTEGER                                  :: isub
    INTEGER(int_8)                           :: n_bytes

    CALL tiset(procedureN,isub)

    IF( .NOT. stream%init ) CALL stopgm(procedureN,"stream not initialized",&
         __LINE__,__FILE__)

    ! allocate CUDA memory for rho
    n_bytes = INT( n, int_8 ) * INT( sizeof_real_8, int_8 )
    CALL cuda_alloc_bytes ( stream%rho_d, n_bytes, device_id )
    CALL cuda_mem_zero_bytes ( stream%rho_d, n_bytes )

    CALL tihalt(procedureN,isub)

  END SUBROUTINE cp_curho_stream_alloc_buffers


  SUBROUTINE cp_curho_finalize_stream ( stream )
    TYPE(cp_curho_stream_t), INTENT(INOUT)   :: stream

    CHARACTER(*), PARAMETER :: procedureN = 'cp_curho_finalize_stream'

    INTEGER                                  :: isub

    CALL tiset(procedureN,isub)

    IF( .NOT. stream%init ) CALL stopgm(procedureN,"stream not initialized",&
         __LINE__,__FILE__)

    !CALL cuda_stream_destroy ( stream%stream_1 )

    IF( stream%rho_d%init ) CALL stopgm(procedureN,"buffer not deallocated",&
         __LINE__,__FILE__)

    stream%init = .FALSE.

    CALL tihalt(procedureN,isub)

  END SUBROUTINE cp_curho_finalize_stream


  SUBROUTINE cp_curho_stream_dealloc_buffers ( stream )
    TYPE(cp_curho_stream_t), INTENT(INOUT)   :: stream

    CHARACTER(*), PARAMETER :: procedureN = 'cp_curho_stream_dealloc_buffers'

    INTEGER                                  :: isub

    CALL tiset(procedureN,isub)

    IF( .NOT. stream%init ) CALL stopgm(procedureN,"stream not initialized",&
         __LINE__,__FILE__)

    CALL cuda_dealloc ( stream%rho_d )

    CALL tihalt(procedureN,isub)

  END SUBROUTINE cp_curho_stream_dealloc_buffers


END MODULE cp_curho_utils
