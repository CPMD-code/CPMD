PROGRAM GPUBandwidth

  USE cuda_types,                      ONLY: cuda_event_t,&
                                             cuda_memory_t,&
                                             cuda_stream_t
  USE cuda_utils,                      ONLY: &
       cuda_alloc_bytes, cuda_alloc_host_portable, cuda_dealloc, &
       cuda_dealloc_host, cuda_device_synchronize, cuda_event_create, &
       cuda_event_destroy, cuda_event_elapsed_time, cuda_event_record, &
       cuda_event_synchronize, cuda_memcpy_async_host_to_device, &
       cuda_stream_create, cuda_stream_destroy
  USE kinds,                           ONLY: int_8,&
                                             real_8
  USE sizeof_kinds,                    ONLY: sizeof_real_8

  IMPLICIT NONE

  INTEGER :: n_streams_per_device, i_stream, i_device, n_devices, max_mem_pow2_real8, i_pow2
  INTEGER :: max_mem_real8, min_mem_pow2_real8, message_size_real8
  INTEGER(int_8) :: max_mem_bytes

  REAL(real_8) :: dt
  REAL(real_8), DIMENSION(:,:,:), POINTER, CONTIGUOUS :: h
  TYPE( cuda_event_t ) :: event_start, event_stop
  TYPE( cuda_memory_t ), DIMENSION(:,:), ALLOCATABLE :: d
  TYPE( cuda_stream_t ), DIMENSION(:,:), ALLOCATABLE :: s


  n_devices = 4
  n_streams_per_device = 4

  min_mem_pow2_real8 = 1
  max_mem_pow2_real8 = 24


  max_mem_real8 = 2**max_mem_pow2_real8


  max_mem_bytes = INT( max_mem_real8, int_8 ) * INT( sizeof_real_8, int_8 )

  CALL cuda_event_create ( event_start )
  CALL cuda_event_create ( event_stop )

  ALLOCATE( d( n_streams_per_device, n_devices ), s( n_streams_per_device, n_devices ) )
  DO i_device = 1, n_devices
     DO i_stream = 1, n_streams_per_device
        CALL cuda_alloc_bytes( d( i_stream, i_device ), max_mem_bytes, i_device-1 )
        CALL cuda_stream_create ( s( i_stream, i_device ), i_device-1 )
     ENDDO
  ENDDO
  CALL cuda_alloc_host_portable ( h, [max_mem_real8,n_streams_per_device,n_devices] )
  h(:,:,:) = 2.0_real_8



  DO i_pow2 = min_mem_pow2_real8, max_mem_pow2_real8
     message_size_real8 = 2**i_pow2

     CALL cuda_device_synchronize ( )
     CALL cuda_event_record( event_start, s( 1, 1 ) )

     DO i_device = 1, n_devices
        DO i_stream = 1, n_streams_per_device
           CALL cuda_memcpy_async_host_to_device ( h(1:message_size_real8,i_stream,i_device), &
                & d( i_stream, i_device ), s( i_stream, i_device ) )
        ENDDO
     ENDDO

     CALL cuda_device_synchronize ( )
     CALL cuda_event_record( event_stop, s( 1, 1 ) )
     CALL cuda_event_synchronize ( event_stop )
     CALL cuda_event_elapsed_time ( event_start, event_stop, dt )

     WRITE(*,'(a,i12,a,f10.6,a,f10.6,a)') 'message_size per stream ',INT(message_size_real8,int_8)*INT(sizeof_real_8, int_8), &
          ' bytes, time = ', dt, &
          ' s, Bwidth = ',REAL( n_devices * n_streams_per_device, real_8 ) * &
          REAL(message_size_real8, real_8) * &
          REAL(sizeof_real_8,real_8 ) / dt * 1.0e-9_real_8, &
          ' GB/s'

  ENDDO





  CALL cuda_event_destroy ( event_start )
  CALL cuda_event_destroy ( event_stop )

  DO i_device = 1, n_devices
     DO i_stream = 1, n_streams_per_device
        CALL cuda_dealloc( d( i_stream, i_device ) )
        CALL cuda_stream_destroy ( s( i_stream, i_device ) )
     ENDDO
  ENDDO
  DEALLOCATE( d )
  CALL cuda_dealloc_host ( h )


END PROGRAM GPUBandwidth
