MODULE cp_cufft_types

  USE cublas_types,                    ONLY: cublas_handle_t
  USE cuda_types,                      ONLY: cuda_memory_t,&
                                             cuda_stream_t
  USE cufft_types,                     ONLY: cufft_plan_t
  USE error_handling,                  ONLY: stopgm

  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER, PRIVATE :: max_cufft_plans = 48

  TYPE, PUBLIC :: plan_info_t
     LOGICAL :: init = .FALSE.
     INTEGER, DIMENSION(6) :: info = 0
  END TYPE plan_info_t

  TYPE, PUBLIC :: cp_cufft_plans_t
     LOGICAL :: init = .FALSE.
     INTEGER :: n_plans = 0
     TYPE( cufft_plan_t ), DIMENSION(max_cufft_plans) :: p
     TYPE( plan_info_t ), DIMENSION(max_cufft_plans) :: i
  END TYPE cp_cufft_plans_t

  !>
  TYPE, PUBLIC :: cp_cufft_stream_t
     LOGICAL :: init = .FALSE.

     ! cuda stream for fft
     TYPE( cuda_stream_t ) :: stream

     ! cublas handle for fft
     TYPE(cublas_handle_t) :: blas_handle

     ! plans for the different ffts
     TYPE( cp_cufft_plans_t ) :: plans

     ! memory needed on gpu to perform ffts
     ! cmplx
     TYPE( cuda_memory_t ) :: t1, t2

  END TYPE cp_cufft_stream_t

  TYPE, PUBLIC :: cp_cufft_device_t
     LOGICAL :: init = .FALSE.

     !device id
     INTEGER :: id = HUGE(0)

     ! global data need on the gpu
     ! int
     TYPE( cuda_memory_t ) :: sp5, sp8, sp9, msqs, msqf, lrxpl, nzfs, inzs

     TYPE( cp_cufft_stream_t ), DIMENSION(:), ALLOCATABLE :: streams

  END TYPE cp_cufft_device_t
  !<

  TYPE, PUBLIC :: cp_cufft_t
     LOGICAL :: init = .FALSE.

     TYPE( cp_cufft_device_t ), DIMENSION(:), ALLOCATABLE :: devices

  END TYPE cp_cufft_t

  TYPE( cp_cufft_t ), PUBLIC, TARGET, SAVE :: cp_cufft

  PUBLIC :: cp_cufft_device_get_ptrs
  PUBLIC :: cp_cufft_stream_get_ptrs

CONTAINS

  SUBROUTINE cp_cufft_device_get_ptrs ( this, device_idx, sp5_d, sp8_d, sp9_d, msqs_d, &
       msqf_d, lrxpl_d, nzfs_d, inzs_d )
    TYPE(cp_cufft_t), INTENT(IN), TARGET     :: this
    INTEGER, INTENT(IN)                      :: device_idx
    TYPE(cuda_memory_t), INTENT(OUT), &
      OPTIONAL, POINTER                      :: sp5_d, sp8_d, sp9_d, msqs_d, &
                                                msqf_d, lrxpl_d, nzfs_d, &
                                                inzs_d

    CHARACTER(*), PARAMETER :: procedureN = 'cp_cufft_device_get_ptrs'

    IF( .NOT. this%init ) CALL stopgm(procedureN,"this not initialized",&
         __LINE__,__FILE__)

    IF(  device_idx > UBOUND( this%devices, 1 ) .OR. &
         device_idx < LBOUND( this%devices, 1 ) ) THEN
       CALL stopgm(procedureN,"device_idx is wrong",&
            __LINE__,__FILE__)
    ENDIF

    IF(PRESENT( sp5_d     )) sp5_d     => this%devices(device_idx)%sp5
    IF(PRESENT( sp8_d     )) sp8_d     => this%devices(device_idx)%sp8
    IF(PRESENT( sp9_d     )) sp9_d     => this%devices(device_idx)%sp9
    IF(PRESENT( msqs_d    )) msqs_d    => this%devices(device_idx)%msqs
    IF(PRESENT( msqf_d    )) msqf_d    => this%devices(device_idx)%msqf
    IF(PRESENT( lrxpl_d   )) lrxpl_d   => this%devices(device_idx)%lrxpl
    IF(PRESENT( nzfs_d    )) nzfs_d    => this%devices(device_idx)%nzfs
    IF(PRESENT( inzs_d    )) inzs_d    => this%devices(device_idx)%inzs
    !IF(PRESENT(      ))      => this%devices(device_idx)%

  END SUBROUTINE cp_cufft_device_get_ptrs


  SUBROUTINE cp_cufft_stream_get_ptrs ( this, device_idx, stream_idx, stream, t1_d, t2_d, blas_handle, plans )
    TYPE(cp_cufft_t), INTENT(IN), TARGET     :: this
    INTEGER, INTENT(IN)                      :: device_idx, stream_idx
    TYPE(cuda_stream_t), INTENT(OUT), &
      OPTIONAL, POINTER                      :: stream
    TYPE(cuda_memory_t), INTENT(OUT), &
      OPTIONAL, POINTER                      :: t1_d, t2_d
    TYPE(cublas_handle_t), INTENT(OUT), &
      OPTIONAL, POINTER                      :: blas_handle
    TYPE(cp_cufft_plans_t), INTENT(OUT), &
      OPTIONAL, POINTER                      :: plans

    CHARACTER(*), PARAMETER :: procedureN = 'cp_cufft_stream_get_ptrs'

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

    IF(PRESENT( t1_d        )) t1_d     => this%devices(device_idx)%streams(stream_idx)%t1
    IF(PRESENT( t2_d        )) t2_d     => this%devices(device_idx)%streams(stream_idx)%t2
    IF(PRESENT( stream      )) stream   => this%devices(device_idx)%streams(stream_idx)%stream
    IF(PRESENT( blas_handle )) blas_handle => this%devices(device_idx)%streams(stream_idx)%blas_handle
    IF(PRESENT( plans       )) plans => this%devices(device_idx)%streams(stream_idx)%plans

  END SUBROUTINE cp_cufft_stream_get_ptrs

END MODULE cp_cufft_types
