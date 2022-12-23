MODULE cp_cuda_utils

  USE nvml_utils, ONLY: nvml_device_get_count, nvml_init, nvml_shutdown, nvml_device_get_handle_by_index, &
       nvml_device_get_uuid, nvml_device_get_pci_info, nvml_device_get_name, nvml_device_get_board_id
  USE cp_cuda_types,                   ONLY: cp_cuda_devices_blas,&
                                             cp_cuda_devices_fft,&
                                             cp_cuda_env
  USE cublas_types,                    ONLY: cublas_handle_t
  USE cublas_utils,                    ONLY: cublas_create,&
                                             cublas_destroy,&
                                             cublas_get_version
  USE cuda_utils,                      ONLY: cuda_device_get_stream_priority_range,&
                                             cuda_driver_get_version,&
                                             cuda_get_device_count,&
                                             cuda_get_device_properties,&
                                             cuda_runtime_get_version, &
                                             cuda_device_get_pci_bus_id
  USE cufft_utils,                     ONLY: cufft_get_version
  USE error_handling,                  ONLY: stopgm
  USE mp_interface,                    ONLY: mp_bcast, mp_sum,&
                                             mp_get_node_env,&
                                             mp_get_processor_name,&
                                             mp_max_processor_name,&
                                             mp_sync
  USE parac,                           ONLY: parai,&
                                             paral
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE string_utils,                    ONLY: int2str
#if defined(_HAS_CUDA)
  USE nvml_types, ONLY: nvml_device_t
#endif


  !$ USE omp_lib, ONLY: omp_get_max_threads

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: cp_cuda_init
  PUBLIC :: cp_cuda_finalize

CONTAINS

  SUBROUTINE cp_cuda_init ( )

    CHARACTER(len=*), PARAMETER              :: procedureN = 'cp_cuda_init'

    CHARACTER(mp_max_processor_name)         :: node_io_name, proc_name
    INTEGER :: cublas_version, cufft_version, device_id, driver_version, &
      greatestPriority, i_device, ierr, iproc, isub, jj, leastPriority, from, &
      max_threads, n_devices, node_numtasks, node_taskid, runtime_version, nvml_device_count
    LOGICAL                                  :: node_io_parent
    TYPE(cublas_handle_t)                    :: blas_handle
    CHARACTER(128) :: pciBusId
    INTEGER, DIMENSION(:), ALLOCATABLE :: io_node_taskid_list, ids
#if defined(_HAS_CUDA)
    TYPE( nvml_device_t ) :: device
    INTEGER :: nvml_device, board_id
    CHARACTER(128) :: uuid, name
#endif
    

    CALL tiset(procedureN,isub)

    max_threads = 1
    !$ max_threads = omp_get_max_threads( )

    ! init cuda
    IF( cp_cuda_env%use_blas .OR. cp_cuda_env%use_fft ) THEN

       !vw use cp_grp here, it shouldnt matter as long as all the tasks used to lunch cpmd are present
       CALL mp_get_node_env ( parai%cp_grp, node_numtasks, node_taskid )
       !>vw get which node is the io node
       CALL mp_get_processor_name ( proc_name )
       IF( paral%io_parent ) node_io_name = proc_name
       CALL mp_bcast( node_io_name, parai%io_source, parai%cp_grp )
       node_io_parent = proc_name == node_io_name
       ALLOCATE(io_node_taskid_list(0:node_numtasks-1))
       io_node_taskid_list(:) = 0
       IF(node_io_parent) io_node_taskid_list(node_taskid) = parai%cp_me
       CALL mp_sum( io_node_taskid_list, SIZE(io_node_taskid_list), parai%cp_grp )

       !<vw
       CALL cuda_get_device_count ( n_devices )

       !vw sannity checks
       IF (paral%io_parent) THEN
          ! FFT
          !vw cyclic distribution in case we attach more tasks than devices
          IF( node_numtasks * cp_cuda_env%fft_n_devices_per_task > n_devices ) &
                                !     & CALL stopgm(procedureN,&
                                !     & 'the number of device requested exceeds the number of device available',&
                                !     & __LINE__,__FILE__)
               & WRITE(6,*) 'WARNING: the number of device requested exceeds the number of device available'

          IF( cp_cuda_env%fft_n_devices_per_task * cp_cuda_env%fft_n_streams_per_device > max_threads ) &
               & CALL stopgm(procedureN,&
               & 'the total number of streams requested exceeds the number of available threads',&
               & __LINE__,__FILE__)

          !BLAS
          !vw cyclic distribution in case we attach more tasks than devices
          IF( node_numtasks * cp_cuda_env%blas_n_devices_per_task > n_devices ) &
                                !     & CALL stopgm(procedureN,&
                                !     & 'the number of device requested exceeds the number of device available',&
                                !     & __LINE__,__FILE__)
               & WRITE(6,*) 'WARNING: the number of device requested exceeds the number of device available'

          IF( cp_cuda_env%blas_n_devices_per_task * cp_cuda_env%blas_n_streams_per_device > max_threads ) &
               & CALL stopgm(procedureN,&
               & 'the total number of streams requested exceeds the number of available threads',&
               & __LINE__,__FILE__)
       ENDIF


       !vw assign devices to tasks
       !FFT
       ALLOCATE(cp_cuda_devices_fft%ids( cp_cuda_env%fft_n_devices_per_task ),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', __LINE__,__FILE__)
       DO i_device = 1, cp_cuda_env%fft_n_devices_per_task
          device_id = i_device - 1
          !vw cyclic distribution in case we attach more tasks than devices
          cp_cuda_devices_fft%ids( i_device ) = MOD( node_taskid * cp_cuda_env%fft_n_devices_per_task + device_id, n_devices )
          !cp_cuda_devices_fft%ids( i_device ) = node_taskid * cp_cuda_env%fft_n_devices_per_task + device_id
       ENDDO
       !BLAS
       ALLOCATE(cp_cuda_devices_blas%ids( cp_cuda_env%blas_n_devices_per_task ),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', __LINE__,__FILE__)
       DO i_device = 1, cp_cuda_env%blas_n_devices_per_task
          device_id = i_device - 1
          !vw cyclic distribution in case we attach more tasks than devices
          cp_cuda_devices_blas%ids( i_device ) = MOD( node_taskid * cp_cuda_env%blas_n_devices_per_task + device_id, n_devices )
          !cp_cuda_devices_blas%ids( i_device ) = node_taskid * cp_cuda_env%blas_n_devices_per_task + device_id
       ENDDO

#if defined(_HAS_CUDA)
       CALL nvml_init ( )
       CALL nvml_device_get_count ( nvml_device_count )
       IF (paral%io_parent) THEN
          WRITE(6,'(1X,A,I0)') 'NVML: number devices per node ', nvml_device_count
       ENDIF
       DO jj=1,100
          FLUSH(6); CALL mp_sync(parai%cp_grp)
       ENDDO
       DO nvml_device = 0, nvml_device_count - 1
          IF (paral%io_parent) THEN
             WRITE(6,'(1X,A,I0)') 'NVML: device id ', nvml_device
          ENDIF
          CALL nvml_device_get_handle_by_index ( nvml_device, device )
          IF (paral%io_parent) THEN
             CALL nvml_device_get_name ( device, name )
             CALL nvml_device_get_uuid ( device, uuid )
             CALL nvml_device_get_board_id ( device, board_id )
             WRITE(6,  '(1X,2A)') 'NVML:    name                          ', TRIM( name )
             WRITE(6,  '(1X,2A)') 'NVML:    UUID                          ', TRIM( uuid )
             WRITE(6,'(1x,A,i0)') 'NVML:    board id                      ', board_id
             CALL nvml_device_get_pci_info ( device )
          ENDIF
          DO jj=1,100
             FLUSH(6); CALL mp_sync(parai%cp_grp)
          ENDDO
       ENDDO
       CALL nvml_shutdown ( )
#endif

       !vw printing device properties
       IF( paral%io_parent ) THEN
          WRITE(6,'(1X,A,I0)') 'CUDA: number devices per node ', n_devices
          WRITE(6,'(1X,A,I0)') 'CUDA: number tasks per node ', node_numtasks
          CALL cuda_device_get_stream_priority_range ( leastPriority, greatestPriority )
          WRITE(6,'(1X,2(A,I0),A)') 'CUDA: stream priority range [', greatestPriority,',', leastPriority,']'
       ENDIF

       DO iproc=0,node_numtasks-1
          from = io_node_taskid_list( iproc )
          IF( paral%io_parent ) THEN
             WRITE(6,'(1X,2(A,I0))') 'CUDA: MPI local ', iproc,' global ', from
          ENDIF
          DO device_id = 0, n_devices - 1
             CALL cuda_device_get_pci_bus_id ( device_id, pciBusId  )
             CALL mp_bcast( pciBusId, from, parai%cp_grp )
             IF( paral%io_parent ) THEN
                WRITE(6,'(1X,A,I0)') 'CUDA:    device id                      ', device_id
                WRITE(6,  '(1X,2A)') 'CUDA:    busId                          ', TRIM( pciBusId )
                CALL cuda_get_device_properties ( device_id )!vw this is wrong!  here *from* should send the info to io_node !
             ENDIF
          ENDDO
       ENDDO

       !blas
       IF (paral%io_parent) THEN
          WRITE(6,'(1X,A,L1)') 'CUDA: BLAS use cublas ', cp_cuda_env%use_blas
          WRITE(6,'(1X,A,I0)') 'CUDA: BLAS number device per task ', cp_cuda_env%blas_n_devices_per_task
          WRITE(6,'(1X,A,I0)') 'CUDA: BLAS number streams per device ', cp_cuda_env%blas_n_streams_per_device

       ENDIF
       CALL mp_sync(parai%cp_grp)


       DO iproc=0,node_numtasks-1
          from = io_node_taskid_list( iproc )
          IF( cp_cuda_env%blas_n_devices_per_task > 0 ) THEN
             ALLOCATE(ids(SIZE(cp_cuda_devices_blas%ids)))
             ids(:) = cp_cuda_devices_blas%ids(:)
             CALL mp_bcast( ids, SIZE(ids), from, parai%cp_grp )
             IF( paral%io_parent ) THEN
                WRITE(6,'(1X,A,'//TRIM(int2str(SIZE(ids)))//'(1X,I0),A,I0)') 'CUDA: BLAS device ids', ids,' for task ', from
             ENDIF
             DEALLOCATE(ids)
          ENDIF
       ENDDO

       !fft
       IF (paral%io_parent) THEN
          WRITE(6,'(1X,A,L1)') 'CUDA: FFT use cufft ', cp_cuda_env%use_fft 
          WRITE(6,'(1X,A,I0)') 'CUDA: FFT number device per task ', cp_cuda_env%fft_n_devices_per_task
          WRITE(6,'(1X,A,I0)') 'CUDA: FFT number streams per device ', cp_cuda_env%fft_n_streams_per_device
       ENDIF

       DO iproc=0,node_numtasks-1
          from = io_node_taskid_list( iproc )
          IF( cp_cuda_env%fft_n_devices_per_task > 0 ) THEN
             ALLOCATE(ids(SIZE(cp_cuda_devices_fft%ids)))
             ids(:) = cp_cuda_devices_fft%ids(:)
             CALL mp_bcast( ids, SIZE(ids), from, parai%cp_grp )
             IF( paral%io_parent ) THEN
                WRITE(6,'(1X,A,'//TRIM(int2str(SIZE(ids)))//'(1X,I0),A,I0)') 'CUDA: FFT device ids', ids,' for task ', from
             ENDIF
             DEALLOCATE(ids)
          ENDIF
       ENDDO

       !vw for the moment stop if n_devices_per_task /= 1 (more work needed)
       !IF (paral%io_parent .AND. cp_cuda_env%n_devices_per_task > 1 ) CALL stopgm(procedureN,&
       !     & 'cp_cuda_env%n_devices_per_task>1 NYI', __LINE__,__FILE__) 

       IF (paral%io_parent) WRITE(6,'(1X,A)') 'CUDA: initialization'

       CALL cuda_driver_get_version ( driver_version )
       CALL cuda_runtime_get_version ( runtime_version )

       IF (paral%io_parent) THEN
          WRITE(6,'(1X,A,I0)') 'CUDA: driver version ', driver_version
          WRITE(6,'(1X,A,I0)') 'CUDA: runtime version ', runtime_version
       ENDIF

    ENDIF

    ! init cublas
    IF( cp_cuda_env%use_blas ) THEN
       IF (paral%io_parent) WRITE(6,'(1X,A)') 'CUBLAS: initialization'

       CALL cublas_create ( blas_handle, 0 )
       CALL cublas_get_version ( blas_handle, cublas_version )
       CALL cublas_destroy ( blas_handle )

       IF (paral%io_parent) WRITE(6,'(1X,A,I0)') 'CUBLAS: version ', cublas_version

       DEALLOCATE(io_node_taskid_list)

    ENDIF

    ! init cufft
    IF( cp_cuda_env%use_fft ) THEN
       IF (paral%io_parent) WRITE(6,'(1X,A)') 'CUFFT: initialization'
       CALL cufft_get_version ( cufft_version )
       IF (paral%io_parent) WRITE(6,'(1X,A,I0)') 'CUFFT: version ', cufft_version
    END IF

    CALL tihalt(procedureN,isub)

  END SUBROUTINE cp_cuda_init

  SUBROUTINE cp_cuda_finalize ( )

    CHARACTER(len=*), PARAMETER :: procedureN = 'cp_cuda_finalize'

    INTEGER                                  :: ierr, isub

    CALL tiset(procedureN,isub)

    ! finalize cublas
    IF( cp_cuda_env%use_blas ) THEN
       IF (paral%io_parent) WRITE(6,*) 'CUBLAS: finalization'
    ENDIF

    ! finalize cufft
    IF( cp_cuda_env%use_fft ) THEN
       IF (paral%io_parent) WRITE(6,*) 'CUFFT: finalization'
    ENDIF

    ! finalize cuda
    IF( cp_cuda_env%use_blas .OR. cp_cuda_env%use_fft ) THEN
       IF (paral%io_parent) WRITE(6,*) 'CUDA: finalization'

       DEALLOCATE(cp_cuda_devices_fft%ids, cp_cuda_devices_blas%ids,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            & __LINE__,__FILE__)

    ENDIF

    CALL tihalt(procedureN,isub)

  END SUBROUTINE cp_cuda_finalize

END MODULE cp_cuda_utils
