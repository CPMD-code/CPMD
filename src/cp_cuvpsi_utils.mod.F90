MODULE cp_cuvpsi_utils

  USE cp_cuda_types,                   ONLY: cp_cuda_devices_t,&
                                             cp_cuda_env_t
  USE cp_cufft_types,                  ONLY: cp_cufft_stream_get_ptrs,&
                                             cp_cufft_t
  USE cp_cuvpsi_types,                 ONLY: cp_cuvpsi_device_t,&
                                             cp_cuvpsi_t
  USE cuda_types,                      ONLY: cuda_memory_t,&
                                             cuda_stream_t
  USE cuda_utils,                      ONLY: cuda_alloc_bytes,&
                                             cuda_dealloc,&
                                             cuda_mem_zero_bytes,&
                                             cuda_memcpy_host_to_device
  USE cuuser_utils,                    ONLY: CuUser_Pointwise_CxR
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: int_8,&
                                             real_8
  USE sizeof_kinds,                    ONLY: sizeof_real_8
  USE thread_view_types,               ONLY: thread_view_t
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE


  PRIVATE


  PUBLIC :: cp_cuapply_potential
  PUBLIC :: cp_cuvpotdg_copy_to_device
  PUBLIC :: cp_cuvpotdg_get_ptrs
  PUBLIC :: cp_cuvpsi_init_and_alloc
  PUBLIC :: cp_cuvpsi_init_and_alloc_device
  PUBLIC :: cp_cuvpsi_finalize_and_dealloc
  PUBLIC :: cp_cuvpsi_finalize_and_dealloc_device


CONTAINS


  SUBROUTINE cp_cuapply_potential ( n, cp_cufft, cp_cuvpsi, thread_view )
    INTEGER, INTENT(IN)                      :: n
    TYPE(cp_cufft_t), INTENT(INOUT)          :: cp_cufft
    TYPE(cp_cuvpsi_t), INTENT(IN)            :: cp_cuvpsi
    TYPE(thread_view_t), INTENT(IN)          :: thread_view

    INTEGER                                  :: device_idx, stream_idx
    TYPE(cuda_memory_t), POINTER             :: psi_d, vpotdg_d
    TYPE(cuda_stream_t), POINTER             :: stream_p

    device_idx = thread_view%device_idx
    stream_idx = thread_view%stream_idx - 1 !vw FIX that -1

    !acm here t1 shall contain the data !
    CALL cp_cufft_stream_get_ptrs ( cp_cufft, device_idx, stream_idx, stream=stream_p, t1_d=psi_d )
    CALL cp_cuvpotdg_get_ptrs ( cp_cuvpsi, device_idx, vpotdg_d=vpotdg_d )
    CALL CuUser_Pointwise_CxR ( psi_d, vpotdg_d, n, stream_p )

  END SUBROUTINE cp_cuapply_potential


  SUBROUTINE cp_cuvpotdg_copy_to_device ( cp_cuda_env, cp_cuvpsi, vpotdg )
    TYPE(cp_cuda_env_t), INTENT(IN)          :: cp_cuda_env
    TYPE(cp_cuvpsi_t), INTENT(IN)            :: cp_cuvpsi
    REAL(real_8), DIMENSION(:, :), &
      INTENT(IN)                             :: vpotdg

    INTEGER                                  :: device_idx
    TYPE(cuda_memory_t), POINTER             :: vpotdg_d

!acm copy FFT arrays to GPU memory

    DO device_idx = 1, cp_cuda_env%fft_n_devices_per_task
       CALL cp_cuvpotdg_get_ptrs ( cp_cuvpsi, device_idx, vpotdg_d=vpotdg_d )
       CALL cuda_memcpy_host_to_device (vpotdg, vpotdg_d )
    ENDDO

  END SUBROUTINE cp_cuvpotdg_copy_to_device


  SUBROUTINE cp_cuvpotdg_get_ptrs ( this, device_idx, vpotdg_d )
    TYPE(cp_cuvpsi_t), INTENT(IN), TARGET    :: this
    INTEGER, INTENT(IN)                      :: device_idx
    TYPE(cuda_memory_t), INTENT(OUT), &
      OPTIONAL, POINTER                      :: vpotdg_d

    CHARACTER(*), PARAMETER :: procedureN = 'cp_cuvpotdg_get_ptrs'

    IF( .NOT. this%init ) CALL stopgm(procedureN,"this not initialized",&
         __LINE__,__FILE__)

    IF(  device_idx > UBOUND( this%devices, 1 ) .OR. &
         device_idx < LBOUND( this%devices, 1 ) ) THEN
       CALL stopgm(procedureN,"device_idx is wrong",&
            __LINE__,__FILE__)
    ENDIF

    IF(PRESENT( vpotdg_d )) vpotdg_d => this%devices(device_idx)%vpotdg_d

  END SUBROUTINE cp_cuvpotdg_get_ptrs


  SUBROUTINE cp_cuvpsi_init_and_alloc ( cp_cuda_env, cp_cuda_devices_fft, cp_cuvpsi, n )
    TYPE(cp_cuda_env_t), INTENT(IN)          :: cp_cuda_env
    TYPE(cp_cuda_devices_t), INTENT(IN)      :: cp_cuda_devices_fft
    TYPE(cp_cuvpsi_t), INTENT(INOUT)         :: cp_cuvpsi
    INTEGER, INTENT(IN)                      :: n

    CHARACTER(*), PARAMETER :: procedureN = 'cp_cuvpsi_init_and_alloc'

    INTEGER                                  :: device_id, i_device, ierr, &
                                                isub

    CALL tiset(procedureN,isub)

    IF( cp_cuvpsi%init ) CALL stopgm(procedureN,"cp_cuvpsi already initialized",&
         __LINE__,__FILE__)

    ALLOCATE( cp_cuvpsi%devices( cp_cuda_env%fft_n_devices_per_task ), STAT=ierr )
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', __LINE__,__FILE__)

    DO i_device = 1, cp_cuda_env%fft_n_devices_per_task
       device_id = cp_cuda_devices_fft%ids( i_device )
       CALL cp_cuvpsi_init_and_alloc_device( cp_cuvpsi%devices( i_device ), device_id, n )
    ENDDO

    cp_cuvpsi%init = .TRUE.

    CALL tihalt(procedureN,isub)

  END SUBROUTINE cp_cuvpsi_init_and_alloc


  SUBROUTINE cp_cuvpsi_init_and_alloc_device( device, device_id, n )
    TYPE(cp_cuvpsi_device_t), INTENT(INOUT)  :: device
    INTEGER, INTENT(IN)                      :: device_id, n

    CHARACTER(*), PARAMETER :: procedureN = 'cp_cuvpsi_init_and_alloc_device'

    INTEGER(int_8)                           :: n_bytes

    IF( device%init ) CALL stopgm(procedureN,"device already initialized",&
         __LINE__,__FILE__)

    device%id = device_id
    device%init = .TRUE.

    n_bytes = INT( n, int_8 ) * INT( sizeof_real_8, int_8 )
    CALL cuda_alloc_bytes ( device%vpotdg_d, n_bytes, device_id )
    CALL cuda_mem_zero_bytes ( device%vpotdg_d, n_bytes )

  END SUBROUTINE cp_cuvpsi_init_and_alloc_device


  SUBROUTINE cp_cuvpsi_finalize_and_dealloc ( cp_cuvpsi )
    TYPE(cp_cuvpsi_t), INTENT(INOUT)         :: cp_cuvpsi

    CHARACTER(*), PARAMETER :: procedureN = 'cp_cuvpsi_finalize_and_dealloc'

    INTEGER                                  :: i_device, ierr, isub

    CALL tiset(procedureN,isub)

    IF( .NOT.cp_cuvpsi%init ) CALL stopgm(procedureN,"cp_cuvpsi not initialized",&
         __LINE__,__FILE__)

    DO i_device = 1, SIZE( cp_cuvpsi%devices, 1 )
       CALL cp_cuvpsi_finalize_and_dealloc_device ( cp_cuvpsi%devices( i_device ) )
    ENDDO

    DEALLOCATE( cp_cuvpsi%devices, STAT=ierr )
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', __LINE__,__FILE__)

    cp_cuvpsi%init = .FALSE.

    CALL tihalt(procedureN,isub)

  END SUBROUTINE cp_cuvpsi_finalize_and_dealloc


  SUBROUTINE cp_cuvpsi_finalize_and_dealloc_device ( device )
    TYPE(cp_cuvpsi_device_t), INTENT(INOUT)  :: device

    CHARACTER(*), PARAMETER :: &
      procedureN = 'cp_cuvpsi_finalize_and_dealloc_device'

    INTEGER                                  :: isub

    CALL tiset(procedureN,isub)

    IF( .NOT.device%init ) CALL stopgm(procedureN,"device not initialized",&
         __LINE__,__FILE__)

    CALL cuda_dealloc ( device%vpotdg_d )

    device%id = HUGE(0)
    device%init = .FALSE.

    CALL tihalt(procedureN,isub)

  END SUBROUTINE cp_cuvpsi_finalize_and_dealloc_device

END MODULE cp_cuvpsi_utils
