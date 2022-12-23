MODULE cp_cuortho_types

  USE cublas_types,                    ONLY: cublas_handle_t
  USE cuda_types,                      ONLY: cuda_event_t,&
                                             cuda_memory_t,&
                                             cuda_stream_t
  USE cusolver_types,                  ONLY: cusolver_handle_t
  USE error_handling,                  ONLY: stopgm

  IMPLICIT NONE

  PRIVATE

  TYPE, PUBLIC :: cp_cuortho_stream_t
     LOGICAL :: init = .FALSE.

     ! cuda stream for ortho
     TYPE( cuda_stream_t ) :: stream_1
     TYPE( cuda_stream_t ) :: stream_2
     TYPE( cuda_stream_t ) :: stream_3

     ! events
     TYPE( cuda_event_t ) :: event_copy_d2h
     TYPE( cuda_event_t ) :: event_copy_h2d

     ! cublas handle for ortho
     TYPE( cublas_handle_t ) :: blas_handle

     ! cusolver handle for ortho
     TYPE( cusolver_handle_t ) :: solver_handle

     ! memory needed on gpu to perform ortho
     ! cmplx
     TYPE( cuda_memory_t ) :: C0_local_d, localw_d, local_h_d, global_h_d

     ! work for solver
     INTEGER :: lwork = HUGE(0)
     TYPE( cuda_memory_t ) :: work_d

  END TYPE cp_cuortho_stream_t

  TYPE, PUBLIC :: cp_cuortho_device_t
     LOGICAL :: init = .FALSE.

     !device id
     INTEGER :: id = HUGE(0)

     TYPE( cp_cuortho_stream_t ), DIMENSION(:), ALLOCATABLE :: streams

  END TYPE cp_cuortho_device_t

  TYPE, PUBLIC :: cp_cuortho_t
     !PRIVATE

     LOGICAL :: init = .FALSE.

     INTEGER :: ldc0l_d = HUGE(0)

     TYPE( cp_cuortho_device_t ), DIMENSION(:), ALLOCATABLE :: devices

  END TYPE cp_cuortho_t

  TYPE(cp_cuortho_t), SAVE, TARGET, PUBLIC :: cp_cuortho

  PUBLIC :: cp_cuortho_get
  PUBLIC :: cp_cuortho_get_buffer_ptrs

CONTAINS

  SUBROUTINE cp_cuortho_get ( this, ldc0l_d )
    TYPE(cp_cuortho_t), INTENT(IN), TARGET   :: this
    INTEGER, INTENT(OUT), OPTIONAL           :: ldc0l_d

    IF(PRESENT(ldc0l_d)) ldc0l_d = this%ldc0l_d

  END SUBROUTINE cp_cuortho_get


  SUBROUTINE cp_cuortho_get_buffer_ptrs ( this, device_idx, stream_idx, stream_1, stream_2, stream_3, &
       event_copy_d2h, event_copy_h2d, blas_handle, solver_handle, C0_local_d, localw_d, local_h_d, &
       global_h_d, work_d, lwork )
    TYPE(cp_cuortho_t), INTENT(IN), TARGET   :: this
    INTEGER, INTENT(IN)                      :: device_idx, stream_idx
    TYPE(cuda_stream_t), INTENT(OUT), &
      OPTIONAL, POINTER                      :: stream_1, stream_2, stream_3
    TYPE(cuda_event_t), INTENT(OUT), &
      OPTIONAL, POINTER                      :: event_copy_d2h, event_copy_h2d
    TYPE(cublas_handle_t), INTENT(OUT), &
      OPTIONAL, POINTER                      :: blas_handle
    TYPE(cusolver_handle_t), INTENT(OUT), &
      OPTIONAL, POINTER                      :: solver_handle
    TYPE(cuda_memory_t), INTENT(OUT), &
      OPTIONAL, POINTER                      :: C0_local_d, localw_d, &
                                                local_h_d, global_h_d, work_d
    INTEGER, INTENT(OUT), OPTIONAL           :: lwork

    CHARACTER(*), PARAMETER :: procedureN = 'cp_cuortho_get_buffer_ptrs'

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

    IF(PRESENT( C0_local_d     )) C0_local_d     => this%devices(device_idx)%streams(stream_idx)%C0_local_d
    IF(PRESENT( localw_d       )) localw_d       => this%devices(device_idx)%streams(stream_idx)%localw_d
    IF(PRESENT( local_h_d      )) local_h_d      => this%devices(device_idx)%streams(stream_idx)%local_h_d
    IF(PRESENT( global_h_d     )) global_h_d     => this%devices(device_idx)%streams(stream_idx)%global_h_d
    IF(PRESENT( work_d         )) work_d         => this%devices(device_idx)%streams(stream_idx)%work_d
    IF(PRESENT( lwork          )) lwork          =  this%devices(device_idx)%streams(stream_idx)%lwork
    IF(PRESENT( stream_1       )) stream_1       => this%devices(device_idx)%streams(stream_idx)%stream_1
    IF(PRESENT( stream_2       )) stream_2       => this%devices(device_idx)%streams(stream_idx)%stream_2
    IF(PRESENT( stream_3       )) stream_3       => this%devices(device_idx)%streams(stream_idx)%stream_3
    IF(PRESENT( blas_handle    )) blas_handle    => this%devices(device_idx)%streams(stream_idx)%blas_handle
    IF(PRESENT( solver_handle  )) solver_handle  => this%devices(device_idx)%streams(stream_idx)%solver_handle
    IF(PRESENT( event_copy_h2d )) event_copy_h2d => this%devices(device_idx)%streams(stream_idx)%event_copy_h2d
    IF(PRESENT( event_copy_d2h )) event_copy_d2h => this%devices(device_idx)%streams(stream_idx)%event_copy_d2h

  END SUBROUTINE cp_cuortho_get_buffer_ptrs

END MODULE cp_cuortho_types
