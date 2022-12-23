MODULE cp_cufft_utils

  USE cp_cuda_types,                   ONLY: cp_cuda_devices_t,&
                                             cp_cuda_env_t
  USE cp_cufft_types,                  ONLY: cp_cufft_device_t,&
                                             cp_cufft_plans_t,&
                                             cp_cufft_stream_t,&
                                             cp_cufft_t
  USE cublas_utils,                    ONLY: cublas_create,&
                                             cublas_destroy,&
                                             cublas_set_stream
  USE cuda_types,                      ONLY: cuda_stream_t
  USE cuda_utils,                      ONLY: cuda_alloc_bytes,&
                                             cuda_dealloc,&
                                             cuda_mem_zero_bytes,&
                                             cuda_stream_create,&
                                             cuda_stream_destroy
  USE cufft_interfaces,                ONLY: CUFFT_Z2Z
  USE cufft_types,                     ONLY: cufft_plan_t
  USE cufft_utils,                     ONLY: cufft_plan_create,&
                                             cufft_plan_destroy,&
                                             cufft_set_stream
  USE error_handling,                  ONLY: stopgm
  USE fft,                             ONLY: lmsqmax,&
                                             lnzs
  USE fft_maxfft,                      ONLY: maxfft
  USE isos,                            ONLY: isos1,&
                                             isos3
  USE kinds,                           ONLY: int_8
  USE parac,                           ONLY: parai,&
                                             paral
  USE sizeof_kinds,                    ONLY: sizeof_complex_8,&
                                             sizeof_int
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: cp_cufft_init_devices
  PUBLIC :: cp_cufft_finalize_devices
  PUBLIC :: cp_cufft_get_plan

  INTEGER, PARAMETER, PRIVATE :: idx_transa=1, idx_transb=2, idx_ldax=3, idx_ldbx=4, idx_m=5, idx_n=6

CONTAINS

  SUBROUTINE cp_cufft_init_devices ( cp_cuda_env, cp_cuda_devices_fft, cp_cufft )
    TYPE(cp_cuda_env_t), INTENT(IN)          :: cp_cuda_env
    TYPE(cp_cuda_devices_t), INTENT(IN)      :: cp_cuda_devices_fft
    TYPE(cp_cufft_t), INTENT(INOUT)          :: cp_cufft

    CHARACTER(*), PARAMETER :: procedureN = 'cp_cufft_init_devices'

    INTEGER                                  :: device_id, i_device, ierr, &
                                                isub

    CALL tiset(procedureN,isub)

    IF(paral%io_parent) WRITE(*,*) procedureN//': n_devices_per_task',cp_cuda_env%fft_n_devices_per_task

    IF( cp_cufft%init ) CALL stopgm(procedureN,"cp_cufft already initialized",&
         __LINE__,__FILE__)

    ALLOCATE( cp_cufft%devices( cp_cuda_env%fft_n_devices_per_task ), STAT=ierr )
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', __LINE__,__FILE__)

    DO i_device = 1, cp_cuda_env%fft_n_devices_per_task
       device_id = cp_cuda_devices_fft%ids( i_device )
       CALL cp_cufft_init_device( cp_cuda_env, cp_cufft%devices( i_device ), device_id )
    ENDDO

    cp_cufft%init = .TRUE.

    CALL tihalt(procedureN,isub)

  END SUBROUTINE cp_cufft_init_devices


  SUBROUTINE cp_cufft_finalize_devices ( cp_cufft )
    TYPE(cp_cufft_t), INTENT(INOUT)          :: cp_cufft

    CHARACTER(*), PARAMETER :: procedureN = 'cp_cufft_finalize_devices'

    INTEGER                                  :: i_device, ierr, isub

    CALL tiset(procedureN,isub)

    IF( .NOT.cp_cufft%init ) CALL stopgm(procedureN,"cp_cufft not initialized",&
         __LINE__,__FILE__)

    DO i_device = 1, SIZE( cp_cufft%devices )
       CALL cp_cufft_finalize_device ( cp_cufft%devices( i_device ) )
    ENDDO

    DEALLOCATE( cp_cufft%devices, STAT=ierr )
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    cp_cufft%init = .FALSE.

    CALL tihalt(procedureN,isub)

  END SUBROUTINE cp_cufft_finalize_devices


  SUBROUTINE cp_cufft_init_device ( cp_cuda_env, device, device_id )
    TYPE(cp_cuda_env_t), INTENT(IN)          :: cp_cuda_env
    TYPE(cp_cufft_device_t), INTENT(INOUT)   :: device
    INTEGER, INTENT(IN)                      :: device_id

    CHARACTER(*), PARAMETER :: procedureN = 'cp_cufft_init_device'

    INTEGER                                  :: i_stream, ierr
    INTEGER(int_8)                           :: n_bytes

    IF(paral%io_parent) WRITE(*,*) procedureN//': device_id',device_id

    IF( device%init ) CALL stopgm(procedureN,"device already initialized",&
         __LINE__,__FILE__)

    device%id = device_id

    ! allocate memory for global data needed to perform FFTs
    !parai%nproc or parai%cp_nproc ?
    n_bytes = INT(parai%cp_nproc,int_8) * INT(sizeof_int,int_8)
    CALL cuda_alloc_bytes ( device%sp5, n_bytes, device_id )
    CALL cuda_alloc_bytes ( device%sp8, n_bytes, device_id )
    CALL cuda_alloc_bytes ( device%sp9, n_bytes, device_id )

    CALL cuda_mem_zero_bytes( device%sp5, n_bytes )
    CALL cuda_mem_zero_bytes( device%sp8, n_bytes )
    CALL cuda_mem_zero_bytes( device%sp9, n_bytes )

    n_bytes = INT(2*parai%cp_nproc,int_8) * INT(sizeof_int,int_8)
    CALL cuda_alloc_bytes ( device%lrxpl, n_bytes, device_id )

    CALL cuda_mem_zero_bytes( device%lrxpl, n_bytes )

    !l=(lmsqmax*parai%nproc)
    n_bytes = INT(lmsqmax*parai%cp_nproc,int_8) * INT(sizeof_int,int_8)
    CALL cuda_alloc_bytes ( device%msqs, n_bytes, device_id )

    CALL cuda_mem_zero_bytes( device%msqs, n_bytes )

    !l=(lmsqmax*parai%nproc)
    n_bytes = INT(lmsqmax*parai%cp_nproc,int_8) * INT(sizeof_int,int_8)
    CALL cuda_alloc_bytes ( device%msqf, n_bytes, device_id )

    CALL cuda_mem_zero_bytes( device%msqf, n_bytes )

    n_bytes = INT(lnzs,int_8) * INT(sizeof_int,int_8)
    CALL cuda_alloc_bytes ( device%nzfs, n_bytes, device_id )
    CALL cuda_alloc_bytes ( device%inzs, n_bytes, device_id )

    CALL cuda_mem_zero_bytes( device%nzfs, n_bytes )
    CALL cuda_mem_zero_bytes( device%inzs, n_bytes )

    ALLOCATE( device%streams( 0:cp_cuda_env%fft_n_streams_per_device-1 ), STAT=ierr )
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', __LINE__,__FILE__)

    DO i_stream = 0, cp_cuda_env%fft_n_streams_per_device-1
       CALL cp_cufft_init_stream( device%streams( i_stream ), device_id )
    ENDDO

    device%init = .TRUE.

  END SUBROUTINE cp_cufft_init_device


  SUBROUTINE cp_cufft_finalize_device ( device )
    TYPE(cp_cufft_device_t), INTENT(INOUT)   :: device

    CHARACTER(*), PARAMETER :: procedureN = 'cp_cufft_finalize_device'

    INTEGER                                  :: i_stream, ierr, isub

    CALL tiset(procedureN,isub)

    IF( .NOT.device%init ) CALL stopgm(procedureN,"device not initialized",&
         __LINE__,__FILE__)


    CALL cuda_dealloc ( device%sp5 )
    CALL cuda_dealloc ( device%sp8 )
    CALL cuda_dealloc ( device%sp9 )

    CALL cuda_dealloc ( device%lrxpl )

    CALL cuda_dealloc ( device%msqs )
    CALL cuda_dealloc ( device%msqf )

    CALL cuda_dealloc ( device%nzfs )
    CALL cuda_dealloc ( device%inzs )

    DO i_stream = LBOUND( device%streams, 1 ), UBOUND( device%streams, 1 )
       CALL cp_cufft_finalize_stream( device%streams( i_stream ) )
    ENDDO

    DEALLOCATE( device%streams, STAT=ierr )
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)

    device%id = HUGE(0)
    device%init = .FALSE.

    CALL tihalt(procedureN,isub)

  END SUBROUTINE cp_cufft_finalize_device


  SUBROUTINE cp_cufft_init_stream( stream, device_id )
    TYPE(cp_cufft_stream_t), INTENT(INOUT)   :: stream
    INTEGER, INTENT(IN)                      :: device_id

    CHARACTER(*), PARAMETER :: procedureN = 'cp_cufft_init_stream'

    INTEGER                                  :: isub
    INTEGER(int_8)                           :: n_bytes

    CALL tiset(procedureN,isub)

    IF(paral%io_parent) WRITE(*,*) procedureN//': alloc stream for device_id',device_id

    IF( stream%init ) CALL stopgm(procedureN,"stream already initialized",&
         __LINE__,__FILE__)

    ! allocate CUDA memory for FFTs
    n_bytes = INT(maxfft,int_8) * INT(sizeof_complex_8,int_8)
    !vw for special usages need to have more memoery allocated on the device
    IF (isos1%tclust.AND. isos3%ps_type.EQ.1) n_bytes = n_bytes * 2_int_8

    CALL cuda_alloc_bytes ( stream%t1, n_bytes, device_id )
    CALL cuda_alloc_bytes ( stream%t2, n_bytes, device_id )

    CALL cuda_mem_zero_bytes( stream%t1, n_bytes )
    CALL cuda_mem_zero_bytes( stream%t2, n_bytes )

    CALL cuda_stream_create ( stream%stream, device_id )

    CALL cublas_create ( stream%blas_handle, device_id )
    CALL cublas_set_stream ( stream%blas_handle, stream%stream )

    CALL cp_cufft_init_plans ( stream%plans, device_id )

    stream%init = .TRUE.

    CALL tihalt(procedureN,isub)

  END SUBROUTINE cp_cufft_init_stream


  SUBROUTINE cp_cufft_init_plans ( plans, device_id )
    TYPE(cp_cufft_plans_t), INTENT(INOUT)    :: plans
    INTEGER, INTENT(IN)                      :: device_id

    CHARACTER(*), PARAMETER :: procedureN = 'cp_cufft_init_plans'

    INTEGER                                  :: i, isub

    CALL tiset(procedureN,isub)

    IF( plans%init ) CALL stopgm(procedureN,"plans already initialized",&
         __LINE__,__FILE__)

    plans%n_plans = 0
    plans%init = .TRUE.
    DO i = 1, SIZE(plans%p)
       plans%p( i )%init = .FALSE.
       plans%i( i )%info(:) = 0
       plans%i( i )%init = .FALSE.
    ENDDO

    CALL tihalt(procedureN,isub)

  END SUBROUTINE cp_cufft_init_plans


  SUBROUTINE cp_cufft_finalize_stream ( stream )
    TYPE(cp_cufft_stream_t), INTENT(INOUT)   :: stream

    CHARACTER(*), PARAMETER :: procedureN = 'cp_cufft_finalize_stream'

    INTEGER                                  :: isub

    CALL tiset(procedureN,isub)

    IF( .NOT. stream%init ) CALL stopgm(procedureN,"stream not initialized",&
         __LINE__,__FILE__)

    CALL cuda_dealloc ( stream%t1 )
    CALL cuda_dealloc ( stream%t2 )

    CALL cuda_stream_destroy ( stream%stream )
    CALL cublas_destroy ( stream%blas_handle )

    CALL cp_cufft_finalize_plans ( stream%plans )

    stream%init = .FALSE.

    CALL tihalt(procedureN,isub)

  END SUBROUTINE cp_cufft_finalize_stream

  SUBROUTINE cp_cufft_finalize_plans ( plans )
    TYPE(cp_cufft_plans_t), INTENT(INOUT)    :: plans

    CHARACTER(*), PARAMETER :: procedureN = 'cp_cufft_finalize_plans'

    INTEGER                                  :: i, isub

    CALL tiset(procedureN,isub)

    IF( .NOT. plans%init ) CALL stopgm(procedureN,"plans not initialized",&
         __LINE__,__FILE__)

    IF(paral%io_parent) WRITE(*,*) procedureN//': n_plans',plans%n_plans
    DO i = 1, SIZE( plans%p )
       IF( plans%p( i )%init ) THEN
          IF(paral%io_parent) WRITE(*,*) procedureN//': ',plans%i( i )%info(:)
          CALL cufft_plan_destroy ( plans%p( i ) )
          plans%i( i )%info(:) = 0
          plans%i( i )%init = .FALSE.
       ENDIF
    ENDDO

    plans%init = .FALSE.
    plans%n_plans = 0

    CALL tihalt(procedureN,isub)

  END SUBROUTINE cp_cufft_finalize_plans


  SUBROUTINE fft_create_cufft_plan ( plan_info, plan, device_id )
    INTEGER, DIMENSION(:), INTENT(IN)        :: plan_info
    TYPE(cufft_plan_t), INTENT(INOUT)        :: plan
    INTEGER, INTENT(IN)                      :: device_id

    CHARACTER(*), PARAMETER :: procedureN = 'fft_create_cufft_plan'

    CHARACTER(len=1)                         :: transa, transb
    INTEGER                                  :: isub, ldax, ldbx, m, n

    CALL tiset(procedureN,isub)


    transa = CHAR( plan_info( idx_transa ) )
    transb = CHAR( plan_info( idx_transb ) )
    ldax = plan_info( idx_ldax )
    ldbx = plan_info( idx_ldbx )
    m = plan_info( idx_m )
    n = plan_info( idx_n )
    !write(*,*) 'fft_create_cufft_plan:',transa,transb,ldax,ldbx,m,n
    IF(transa.EQ.'N'.OR.transa.EQ.'n') THEN
       IF(transb.EQ.'N'.OR.transb.EQ.'n') THEN
          CALL cufft_plan_create(plan, 1, n, n, 1, ldax, n, 1, ldbx, CUFFT_Z2Z, m, device_id )
       ELSE
          CALL cufft_plan_create(plan, 1, n, n, 1, ldax, n, ldbx, 1, CUFFT_Z2Z, m, device_id )
       ENDIF
    ELSE
       IF(transb.EQ.'N'.OR.transb.EQ.'n') THEN
          CALL cufft_plan_create(plan, 1, n, n, ldax, 1, n, 1, ldbx, CUFFT_Z2Z, m, device_id )
       ELSE
          CALL cufft_plan_create(plan, 1, n, n, ldax, 1, n, ldbx, 1, CUFFT_Z2Z, m, device_id )
       ENDIF
    ENDIF

    CALL tihalt(procedureN,isub)

  END SUBROUTINE fft_create_cufft_plan


  SUBROUTINE cp_cufft_get_plan ( transa, transb, ldax, ldbx, n, m, plan_p, plan_d, stream )
    USE cp_cuda_types, ONLY: cp_cuda_env
    CHARACTER(len=1), INTENT(IN)             :: transa, transb
    INTEGER, INTENT(IN)                      :: ldax, ldbx, n, m
    TYPE(cufft_plan_t), INTENT(OUT), POINTER :: plan_p
    TYPE(cp_cufft_plans_t), INTENT(INOUT), &
      TARGET                                 :: plan_d
    TYPE(cuda_stream_t), INTENT(IN)          :: stream

    INTEGER                                  :: i_plan
    INTEGER, DIMENSION(6)                    :: plan_info
    LOGICAL                                  :: found

    NULLIFY( plan_p )

    IF( cp_cuda_env%use_fft ) THEN

       ! create plan info vector
       plan_info( idx_transa ) = ICHAR('N')
       IF( transa.EQ.'T'.OR.transa.EQ.'t' ) plan_info( idx_transa ) = ICHAR('T')
       plan_info( idx_transb ) = ICHAR('N')
       IF( transb.EQ.'T'.OR.transb.EQ.'t' ) plan_info( idx_transb ) = ICHAR('T')
       plan_info( idx_ldax ) = ldax
       plan_info( idx_ldbx ) = ldbx
       plan_info( idx_m ) = m
       plan_info( idx_n ) = n

       ! search if present
       found = .FALSE.
       DO i_plan = 1, plan_d%n_plans
          IF( .NOT. plan_d%i( i_plan )%init ) EXIT
          IF( ALL( plan_d%i( i_plan )%info == plan_info ) ) THEN
             plan_p => plan_d%p( i_plan )
             found = .TRUE.
             EXIT
          ENDIF
       ENDDO

       ! if not present create it
       IF( .NOT. found ) THEN
          plan_d%n_plans = plan_d%n_plans + 1
          i_plan = plan_d%n_plans
          !write(*,*) 'fft_get_cufft_plan: plan_info',plan_info
          CALL fft_create_cufft_plan ( plan_info, plan_d%p( i_plan ), stream%device )
          CALL cufft_set_stream(plan_d%p( i_plan ), stream )
          plan_p => plan_d%p( i_plan )
          plan_d%i( i_plan )%info = plan_info
          plan_d%i( i_plan )%init = .TRUE.
       ENDIF

    ENDIF

  END SUBROUTINE cp_cufft_get_plan

END MODULE cp_cufft_utils
