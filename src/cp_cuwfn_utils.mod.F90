MODULE cp_cuwfn_utils

  USE cp_cuda_types,                   ONLY: cp_cuda_devices_t,&
                                             cp_cuda_env_t
  USE cp_cuwfn_types,                  ONLY: cp_cuwfn_device_t,&
                                             cp_cuwfn_t
  USE cuda_utils,                      ONLY: cuda_alloc_bytes,&
                                             cuda_dealloc
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: int_8
  USE sizeof_kinds,                    ONLY: sizeof_complex_8
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: cp_cuwfn_init
  PUBLIC :: cp_cuwfn_finalize
  PUBLIC :: cp_cuwfn_alloc_buffers
  PUBLIC :: cp_cuwfn_dealloc_buffers


CONTAINS


  SUBROUTINE cp_cuwfn_init ( cp_cuda_env, cp_cuda_devices_blas, cp_cuwfn )
    TYPE(cp_cuda_env_t), INTENT(IN)          :: cp_cuda_env
    TYPE(cp_cuda_devices_t), INTENT(IN)      :: cp_cuda_devices_blas
    TYPE(cp_cuwfn_t), INTENT(INOUT)          :: cp_cuwfn

    CHARACTER(*), PARAMETER                  :: procedureN = 'cp_cuwfn_init'

    INTEGER                                  :: device_id, i_device, ierr, &
                                                isub

    CALL tiset(procedureN,isub)

    IF( cp_cuwfn%init ) CALL stopgm(procedureN,"cp_cuwfn already initialized",&
         __LINE__,__FILE__)

    ALLOCATE( cp_cuwfn%devices( cp_cuda_env%blas_n_devices_per_task ), STAT=ierr )
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', __LINE__,__FILE__)

    DO i_device = 1, cp_cuda_env%blas_n_devices_per_task
       device_id = cp_cuda_devices_blas%ids( i_device )
       CALL cp_cuwfn_init_device( cp_cuda_env, cp_cuwfn%devices( i_device ), device_id )
    ENDDO

    cp_cuwfn%init = .TRUE.

    CALL tihalt(procedureN,isub)

  END SUBROUTINE cp_cuwfn_init


  SUBROUTINE cp_cuwfn_alloc_buffers ( cp_cuwfn, n )
    TYPE(cp_cuwfn_t), INTENT(INOUT)          :: cp_cuwfn
    INTEGER, INTENT(IN)                      :: n

    CHARACTER(*), PARAMETER :: procedureN = 'cp_cuwfn_alloc_buffers'

    INTEGER                                  :: i_device, isub

    CALL tiset(procedureN,isub)

    IF( .NOT. cp_cuwfn%init ) CALL stopgm(procedureN,"cp_cuwfn not initialized",&
         __LINE__,__FILE__)

    DO i_device = 1, SIZE( cp_cuwfn%devices, 1 )
       CALL cp_cuwfn_device_alloc_buffers ( cp_cuwfn%devices( i_device ), n )
    ENDDO

    CALL tihalt(procedureN,isub)

  END SUBROUTINE cp_cuwfn_alloc_buffers


  SUBROUTINE cp_cuwfn_finalize ( cp_cuwfn )
    TYPE(cp_cuwfn_t), INTENT(INOUT)          :: cp_cuwfn

    CHARACTER(*), PARAMETER :: procedureN = 'cp_cuwfn_finalize'

    INTEGER                                  :: i_device, ierr, isub

    CALL tiset(procedureN,isub)

    IF( .NOT.cp_cuwfn%init ) CALL stopgm(procedureN,"cp_cuwfn not initialized",&
         __LINE__,__FILE__)

    DO i_device = 1, SIZE( cp_cuwfn%devices, 1 )
       CALL cp_cuwfn_finalize_device ( cp_cuwfn%devices( i_device ) )
    ENDDO

    DEALLOCATE( cp_cuwfn%devices, STAT=ierr )
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', __LINE__,__FILE__)

    cp_cuwfn%init = .FALSE.

    CALL tihalt(procedureN,isub)

  END SUBROUTINE cp_cuwfn_finalize


  SUBROUTINE cp_cuwfn_dealloc_buffers ( cp_cuwfn )
    TYPE(cp_cuwfn_t), INTENT(INOUT)          :: cp_cuwfn

    CHARACTER(*), PARAMETER :: procedureN = 'cp_cuwfn_dealloc_buffers'

    INTEGER                                  :: i_device, isub

    CALL tiset(procedureN,isub)

    IF( .NOT.cp_cuwfn%init ) CALL stopgm(procedureN,"cp_cuwfn not initialized",&
         __LINE__,__FILE__)

    DO i_device = 1, SIZE( cp_cuwfn%devices, 1 )
       CALL cp_cuwfn_device_dealloc_buffers ( cp_cuwfn%devices( i_device ) )
    ENDDO

    CALL tihalt(procedureN,isub)

  END SUBROUTINE cp_cuwfn_dealloc_buffers


  SUBROUTINE cp_cuwfn_init_device( cp_cuda_env, device, device_id )
    TYPE(cp_cuda_env_t), INTENT(IN)          :: cp_cuda_env
    TYPE(cp_cuwfn_device_t), INTENT(INOUT)   :: device
    INTEGER, INTENT(IN)                      :: device_id

    CHARACTER(*), PARAMETER :: procedureN = 'cp_cuwfn_init_device'

    IF( device%init ) CALL stopgm(procedureN,"device already initialized",&
         __LINE__,__FILE__)

    device%id = device_id

    !ALLOCATE( device%streams( 0:cp_cuda_env%fft_n_streams_per_device-1 ), STAT=ierr )
    !IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', __LINE__,__FILE__)

    !DO i_stream = 0, cp_cuda_env%fft_n_streams_per_device-1
    !   CALL cp_cuwfn_init_stream ( device%streams( i_stream ), device_id )
    !ENDDO

    device%init = .TRUE.

  END SUBROUTINE cp_cuwfn_init_device


  SUBROUTINE cp_cuwfn_device_alloc_buffers ( device, n )
    TYPE(cp_cuwfn_device_t), INTENT(INOUT)   :: device
    INTEGER, INTENT(IN)                      :: n

    CHARACTER(*), PARAMETER :: procedureN = 'cp_cuwfn_device_alloc_buffers'

    INTEGER                                  :: device_id
    INTEGER(int_8)                           :: n_bytes

    IF( .NOT. device%init ) CALL stopgm(procedureN,"device not initialized",&
         __LINE__,__FILE__)

    device_id = device%id

    n_bytes = INT( n, int_8 ) * INT( sizeof_complex_8, int_8 )
    CALL cuda_alloc_bytes ( device%c0_d, n_bytes, device_id )

    !DO i_stream = LBOUND( device%streams, 1 ), UBOUND( device%streams, 1 )
    !   CALL cp_cuwfn_stream_alloc_buffers ( device%streams( i_stream ), device_id, n )
    !ENDDO

  END SUBROUTINE cp_cuwfn_device_alloc_buffers


  SUBROUTINE cp_cuwfn_finalize_device ( device )
    TYPE(cp_cuwfn_device_t), INTENT(INOUT)   :: device

    CHARACTER(*), PARAMETER :: procedureN = 'cp_cuwfn_finalize_device'

    INTEGER                                  :: isub

    CALL tiset(procedureN,isub)

    IF( .NOT.device%init ) CALL stopgm(procedureN,"device not initialized",&
         __LINE__,__FILE__)

    !DO i_stream = LBOUND( device%streams,1 ), UBOUND( device%streams,1 )
    !   CALL cp_cuwfn_finalize_stream ( device%streams( i_stream ) )
    !ENDDO

    !DEALLOCATE( device%streams, STAT=ierr )
    !IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', __LINE__,__FILE__)

    device%id = HUGE(0)
    device%init = .FALSE.

    CALL tihalt(procedureN,isub)

  END SUBROUTINE cp_cuwfn_finalize_device


  SUBROUTINE cp_cuwfn_device_dealloc_buffers ( device )
    TYPE(cp_cuwfn_device_t), INTENT(INOUT)   :: device

    CHARACTER(*), PARAMETER :: procedureN = 'cp_cuwfn_device_dealloc_buffers'

    INTEGER                                  :: isub

    CALL tiset(procedureN,isub)

    IF( .NOT.device%init ) CALL stopgm(procedureN,"device not initialized",&
         __LINE__,__FILE__)

    CALL cuda_dealloc ( device%c0_d )

    !DO i_stream = LBOUND( device%streams, 1 ), UBOUND( device%streams, 1 )
    !   CALL cp_cuwfn_stream_dealloc_buffers ( device%streams( i_stream ) )
    !ENDDO

    CALL tihalt(procedureN,isub)

  END SUBROUTINE cp_cuwfn_device_dealloc_buffers


END MODULE cp_cuwfn_utils
