MODULE cp_cuortho_utils

  USE cp_cuda_types,                   ONLY: cp_cuda_devices_t,&
                                             cp_cuda_env_t
  USE cp_cuortho_types,                ONLY: cp_cuortho_device_t,&
                                             cp_cuortho_stream_t,&
                                             cp_cuortho_t
  USE cublas_utils,                    ONLY: cublas_create,&
                                             cublas_destroy,&
                                             cublas_set_stream
  USE cuda_utils,                      ONLY: cuda_alloc_bytes,&
                                             cuda_alloc_pitch_bytes,&
                                             cuda_dealloc,&
                                             cuda_event_create_with_no_timing,&
                                             cuda_event_destroy,&
                                             cuda_stream_create,&
                                             cuda_stream_destroy
  USE cusolver_utils,                  ONLY: cusolver_create,&
                                             cusolver_destroy,&
                                             cusolver_dpotrf_buffersize,&
                                             cusolver_set_stream
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: int_8
  USE sizeof_kinds,                    ONLY: sizeof_complex_8,&
                                             sizeof_real_8
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: cp_cuortho_init
  PUBLIC :: cp_cuortho_finalize
  PUBLIC :: cp_cuortho_alloc_buffers
  PUBLIC :: cp_cuortho_dealloc_buffers

CONTAINS


  SUBROUTINE cp_cuortho_init ( cp_cuda_env, cp_cuda_devices_blas, cp_cuortho )
    TYPE(cp_cuda_env_t), INTENT(IN)          :: cp_cuda_env
    TYPE(cp_cuda_devices_t), INTENT(IN)      :: cp_cuda_devices_blas
    TYPE(cp_cuortho_t), INTENT(INOUT)        :: cp_cuortho

    CHARACTER(*), PARAMETER                  :: procedureN = 'cp_cuortho_init'

    INTEGER                                  :: device_id, i_device, ierr, &
                                                isub

    CALL tiset(procedureN,isub)

    IF( cp_cuortho%init ) CALL stopgm(procedureN,"cp_cuortho already initialized",&
         __LINE__,__FILE__)

    ALLOCATE( cp_cuortho%devices( cp_cuda_env%blas_n_devices_per_task ), STAT=ierr )
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', __LINE__,__FILE__)

    DO i_device = 1, cp_cuda_env%blas_n_devices_per_task
       device_id = cp_cuda_devices_blas%ids( i_device )
       CALL cp_cuortho_init_device( cp_cuda_env, cp_cuortho%devices( i_device ), device_id )
    ENDDO

    cp_cuortho%init = .TRUE.

    CALL tihalt(procedureN,isub)

  END SUBROUTINE cp_cuortho_init


  SUBROUTINE cp_cuortho_alloc_buffers ( cp_cuortho, NGWK_local, nstate, norbx, ldc0l_d )
    TYPE(cp_cuortho_t), INTENT(INOUT)        :: cp_cuortho
    INTEGER, INTENT(IN)                      :: NGWK_local, nstate, norbx
    INTEGER, INTENT(OUT)                     :: ldc0l_d

    CHARACTER(*), PARAMETER :: procedureN = 'cp_cuortho_alloc_buffers'

    INTEGER                                  :: i_device, isub

    CALL tiset(procedureN,isub)

    IF( .NOT. cp_cuortho%init ) CALL stopgm(procedureN,"cp_cuortho not initialized",&
         __LINE__,__FILE__)

    DO i_device = 1, SIZE( cp_cuortho%devices, 1 )
       CALL cp_cuortho_device_alloc_buffers ( cp_cuortho%devices( i_device ), NGWK_local, nstate, norbx, ldc0l_d  )
    ENDDO

    cp_cuortho%ldc0l_d = ldc0l_d

    CALL tihalt(procedureN,isub)

  END SUBROUTINE cp_cuortho_alloc_buffers


  SUBROUTINE cp_cuortho_finalize ( cp_cuortho )
    TYPE(cp_cuortho_t), INTENT(INOUT)        :: cp_cuortho

    CHARACTER(*), PARAMETER :: procedureN = 'cp_cuortho_finalize'

    INTEGER                                  :: i_device, ierr, isub

    CALL tiset(procedureN,isub)

    IF( .NOT.cp_cuortho%init ) CALL stopgm(procedureN,"cp_cuortho not initialized",&
         __LINE__,__FILE__)

    DO i_device = 1, SIZE( cp_cuortho%devices, 1 )
       CALL cp_cuortho_finalize_device ( cp_cuortho%devices( i_device ) )
    ENDDO

    DEALLOCATE( cp_cuortho%devices, STAT=ierr )
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', __LINE__,__FILE__)

    cp_cuortho%init = .FALSE.

    CALL tihalt(procedureN,isub)

  END SUBROUTINE cp_cuortho_finalize


  SUBROUTINE cp_cuortho_dealloc_buffers ( cp_cuortho )
    TYPE(cp_cuortho_t), INTENT(INOUT)        :: cp_cuortho

    CHARACTER(*), PARAMETER :: procedureN = 'cp_cuortho_dealloc_buffers'

    INTEGER                                  :: i_device, isub

    CALL tiset(procedureN,isub)

    IF( .NOT.cp_cuortho%init ) CALL stopgm(procedureN,"cp_cuortho not initialized",&
         __LINE__,__FILE__)

    DO i_device = 1, SIZE( cp_cuortho%devices, 1 )
       CALL cp_cuortho_device_dealloc_buffers ( cp_cuortho%devices( i_device ) )
    ENDDO
    cp_cuortho%ldc0l_d = HUGE(0)

    CALL tihalt(procedureN,isub)

  END SUBROUTINE cp_cuortho_dealloc_buffers


  SUBROUTINE cp_cuortho_init_device( cp_cuda_env, device, device_id )
    TYPE(cp_cuda_env_t), INTENT(IN)          :: cp_cuda_env
    TYPE(cp_cuortho_device_t), INTENT(INOUT) :: device
    INTEGER, INTENT(IN)                      :: device_id

    CHARACTER(*), PARAMETER :: procedureN = 'cp_cuortho_init_device'

    INTEGER                                  :: i_stream, ierr

    IF( device%init ) CALL stopgm(procedureN,"device already initialized",&
         __LINE__,__FILE__)

    device%id = device_id

    ALLOCATE( device%streams( 0:cp_cuda_env%blas_n_streams_per_device-1 ), STAT=ierr )
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', __LINE__,__FILE__)

    DO i_stream = 0, cp_cuda_env%blas_n_streams_per_device-1
       CALL cp_cuortho_init_stream ( device%streams( i_stream ), device_id )
    ENDDO

    device%init = .TRUE.

  END SUBROUTINE cp_cuortho_init_device


  SUBROUTINE cp_cuortho_device_alloc_buffers ( device, NGWK_local, nstate, norbx, ldc0l_d )
    TYPE(cp_cuortho_device_t), INTENT(INOUT) :: device
    INTEGER, INTENT(IN)                      :: NGWK_local, nstate, norbx
    INTEGER, INTENT(OUT)                     :: ldc0l_d

    CHARACTER(*), PARAMETER :: procedureN = 'cp_cuortho_device_alloc_buffers'

    INTEGER                                  :: device_id, i_stream

    IF( .NOT. device%init ) CALL stopgm(procedureN,"device not initialized",&
         __LINE__,__FILE__)

    device_id = device%id

    DO i_stream = LBOUND( device%streams, 1 ), UBOUND( device%streams, 1 )
       CALL cp_cuortho_stream_alloc_buffers ( device%streams( i_stream ), device_id, NGWK_local, nstate, norbx, ldc0l_d )
    ENDDO

  END SUBROUTINE cp_cuortho_device_alloc_buffers


  SUBROUTINE cp_cuortho_finalize_device ( device )
    TYPE(cp_cuortho_device_t), INTENT(INOUT) :: device

    CHARACTER(*), PARAMETER :: procedureN = 'cp_cuortho_finalize_device'

    INTEGER                                  :: i_stream, ierr, isub

    CALL tiset(procedureN,isub)

    IF( .NOT.device%init ) CALL stopgm(procedureN,"device not initialized",&
         __LINE__,__FILE__)

    DO i_stream = LBOUND( device%streams,1 ), UBOUND( device%streams,1 )
       CALL cp_cuortho_finalize_stream ( device%streams( i_stream ) )
    ENDDO

    DEALLOCATE( device%streams, STAT=ierr )
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', __LINE__,__FILE__)

    device%id = HUGE(0)
    device%init = .FALSE.

    CALL tihalt(procedureN,isub)

  END SUBROUTINE cp_cuortho_finalize_device


  SUBROUTINE cp_cuortho_device_dealloc_buffers ( device )
    TYPE(cp_cuortho_device_t), INTENT(INOUT) :: device

    CHARACTER(*), PARAMETER :: &
      procedureN = 'cp_cuortho_device_dealloc_buffers'

    INTEGER                                  :: i_stream, isub

    CALL tiset(procedureN,isub)

    IF( .NOT.device%init ) CALL stopgm(procedureN,"device not initialized",&
         __LINE__,__FILE__)

    DO i_stream = LBOUND( device%streams, 1 ), UBOUND( device%streams, 1 )
       CALL cp_cuortho_stream_dealloc_buffers ( device%streams( i_stream ) )
    ENDDO

    CALL tihalt(procedureN,isub)

  END SUBROUTINE cp_cuortho_device_dealloc_buffers


  SUBROUTINE cp_cuortho_init_stream( stream, device_id )
    TYPE(cp_cuortho_stream_t), INTENT(INOUT) :: stream
    INTEGER, INTENT(IN)                      :: device_id

    CHARACTER(*), PARAMETER :: procedureN = 'cp_cuortho_init_stream'

    INTEGER                                  :: isub

    CALL tiset(procedureN,isub)

    IF( stream%init ) CALL stopgm(procedureN,"stream already initialized",&
         __LINE__,__FILE__)

    CALL cuda_stream_create ( stream%stream_1, device_id )
    CALL cuda_stream_create ( stream%stream_2, device_id )
    CALL cuda_stream_create ( stream%stream_3, device_id )

    CALL cuda_event_create_with_no_timing ( stream%event_copy_h2d )
    CALL cuda_event_create_with_no_timing ( stream%event_copy_d2h )

    CALL cublas_create ( stream%blas_handle, device_id )
    CALL cublas_set_stream ( stream%blas_handle, stream%stream_3 )

    CALL cusolver_create ( stream%solver_handle, device_id )
    CALL cusolver_set_stream ( stream%solver_handle, stream%stream_3 )

    stream%init = .TRUE.

    CALL tihalt(procedureN,isub)

  END SUBROUTINE cp_cuortho_init_stream


  SUBROUTINE cp_cuortho_stream_alloc_buffers( stream, device_id, NGWK_local, nstate, norbx, ldc0l_d )
    TYPE(cp_cuortho_stream_t), INTENT(INOUT) :: stream
    INTEGER, INTENT(IN)                      :: device_id, NGWK_local, &
                                                nstate, norbx
    INTEGER, INTENT(OUT)                     :: ldc0l_d

    CHARACTER(*), PARAMETER :: procedureN = 'cp_cuortho_stream_alloc_buffers'

    INTEGER                                  :: c0l_col, isub, lwork
    INTEGER(int_8)                           :: c0l_row_bytes, ldc0l_bytes_d, &
                                                n_bytes

    CALL tiset(procedureN,isub)

    IF( .NOT. stream%init ) CALL stopgm(procedureN,"stream not initialized",&
         __LINE__,__FILE__)

    ! allocate CUDA memory for disortho
    c0l_row_bytes = INT( NGWK_local, int_8 ) * INT( sizeof_complex_8, int_8 )
    c0l_col = nstate
    !vw dont need to allocate that much memory
    CALL cuda_alloc_pitch_bytes ( stream%C0_local_d, c0l_row_bytes, c0l_col, ldc0l_bytes_d, device_id )
    ldc0l_d = INT ( ldc0l_bytes_d / INT( sizeof_complex_8, int_8 ) )

    n_bytes = INT( (nstate-((nstate/norbx)-1)*norbx)*ldc0l_d, int_8 ) * INT( sizeof_complex_8, int_8 )
    CALL cuda_alloc_bytes ( stream%localw_d, n_bytes, device_id )

    n_bytes = INT( (nstate-((nstate/norbx)-1)*norbx)*nstate, int_8 ) * INT( sizeof_real_8, int_8 )
    CALL cuda_alloc_bytes ( stream%local_h_d, n_bytes, device_id )
    CALL cuda_alloc_bytes ( stream%global_h_d, n_bytes, device_id )

    CALL cusolver_dpotrf_buffersize ( stream%solver_handle, 'U', norbx, stream%global_h_d, norbx, lwork )
    n_bytes = INT( lwork, int_8 ) * INT( sizeof_real_8, int_8 )
    stream%lwork = lwork
    CALL cuda_alloc_bytes ( stream%work_d, n_bytes, device_id )

    CALL tihalt(procedureN,isub)

  END SUBROUTINE cp_cuortho_stream_alloc_buffers


  SUBROUTINE cp_cuortho_finalize_stream ( stream )
    TYPE(cp_cuortho_stream_t), INTENT(INOUT) :: stream

    CHARACTER(*), PARAMETER :: procedureN = 'cp_cuortho_finalize_stream'

    INTEGER                                  :: isub

    CALL tiset(procedureN,isub)

    IF( .NOT. stream%init ) CALL stopgm(procedureN,"stream not initialized",&
         __LINE__,__FILE__)

    CALL cuda_stream_destroy ( stream%stream_1 )
    CALL cuda_stream_destroy ( stream%stream_2 )
    CALL cuda_stream_destroy ( stream%stream_3 )

    CALL cuda_event_destroy ( stream%event_copy_h2d )
    CALL cuda_event_destroy ( stream%event_copy_d2h )

    CALL cublas_destroy ( stream%blas_handle )

    stream%lwork = HUGE(0)
    CALL cusolver_destroy ( stream%solver_handle )

    IF( stream%C0_local_d%init .OR. stream%localw_d%init   .OR. &
         stream%local_h_d%init  .OR. stream%global_h_d%init .OR. &
         stream%work_d%init ) CALL stopgm(procedureN,"buffer not deallocated",&
         __LINE__,__FILE__)

    stream%init = .FALSE.

    CALL tihalt(procedureN,isub)

  END SUBROUTINE cp_cuortho_finalize_stream


  SUBROUTINE cp_cuortho_stream_dealloc_buffers ( stream )
    TYPE(cp_cuortho_stream_t), INTENT(INOUT) :: stream

    CHARACTER(*), PARAMETER :: &
      procedureN = 'cp_cuortho_stream_dealloc_buffers'

    INTEGER                                  :: isub

    CALL tiset(procedureN,isub)

    IF( .NOT. stream%init ) CALL stopgm(procedureN,"stream not initialized",&
         __LINE__,__FILE__)

    CALL cuda_dealloc ( stream%C0_local_d )
    CALL cuda_dealloc ( stream%localw_d )
    CALL cuda_dealloc ( stream%local_h_d )
    CALL cuda_dealloc ( stream%global_h_d )
    CALL cuda_dealloc ( stream%work_d )

    CALL tihalt(procedureN,isub)

  END SUBROUTINE cp_cuortho_stream_dealloc_buffers


END MODULE cp_cuortho_utils
