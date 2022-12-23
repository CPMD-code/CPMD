#include "cpmd_global.h"

MODULE disortho_utils
  USE cp_cuda_types,                   ONLY: cp_cuda_env
  USE cp_cuortho_types,                ONLY: cp_cuortho,&
                                             cp_cuortho_get,&
                                             cp_cuortho_get_buffer_ptrs
  USE cp_grp_utils,                    ONLY: cp_grp_get_sizes,&
                                             cp_grp_redist
  USE cublas_types,                    ONLY: cublas_handle_t
  USE cublas_utils,                    ONLY: cublas_dgemm,&
                                             cublas_dger,&
                                             cublas_dsyrk,&
                                             cublas_dtrsm,&
                                             cublas_zcopy
  USE cuda_types,                      ONLY: cuda_event_t,&
                                             cuda_memory_t,&
                                             cuda_stream_t
  USE cuda_utils,                      ONLY: &
       cuda_alloc_host, cuda_dealloc_host, cuda_event_record, &
       cuda_mem_zero_bytes, cuda_memcpy2d_async_device_to_host, &
       cuda_memcpy2d_async_host_to_device, cuda_memcpy_async_device_to_host, &
       cuda_memcpy_async_host_to_device, cuda_stream_synchronize, &
       cuda_stream_wait_event, cuda_z_points_to
  USE cusolver_types,                  ONLY: cusolver_handle_t
  USE cuuser_utils,                    ONLY: cuuser_identity
  USE error_handling,                  ONLY: stopgm
  USE gsortho_utils,                   ONLY: gsortho
  USE jrotation_utils,                 ONLY: set_orbdist
  USE kinds,                           ONLY: real_8
  USE linalg_utils,                    ONLY: trans_da
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum,&
                                             mp_sync
  USE ovlap_utils,                     ONLY: dmatc
  USE parac,                           ONLY: parai,&
                                             paral
  USE part_1d,                         ONLY: part_1d_get_blk_bounds
  USE string_utils,                    ONLY: int2str
  USE system,                          ONLY: cnti,&
                                             ncpw,&
                                             nkpt,&
                                             parap,&
                                             paraw
  USE thread_view_types,               ONLY: thread_view_t
  USE thread_view_utils,               ONLY: thread_views_finalize,&
                                             thread_views_init
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  !$ USE omp_lib, ONLY: omp_get_max_threads, omp_get_thread_num, &
  !$                    omp_get_num_threads, omp_set_num_threads
#ifdef _HASNT_OMP_SET_NESTED
  !$ USE omp_lib, ONLY: omp_get_max_active_levels, &
  !$                    omp_set_max_active_levels
#else
  !$ USE omp_lib, ONLY: omp_get_max_active_levels, omp_get_nested, &
  !$                    omp_set_max_active_levels, omp_set_nested
#endif

  USE nvtx_utils

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: disortho_new
  PUBLIC :: disortho_old

CONTAINS

  SUBROUTINE disortho_new ( nstate,c0 )
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(:,:)

    CHARACTER(len=*), PARAMETER              :: procedureN = 'disortho_new'

    IF( cp_cuda_env%use_blas ) THEN
       CALL disortho_new_cuda_low(nstate,c0)
    ELSE
       CALL disortho_new_low(nstate,c0)
    ENDIF

  END SUBROUTINE disortho_new

  ! ==================================================================
  SUBROUTINE disortho_new_low(nstate,c0)
    ! ==--------------------------------------------------------------==
    ! ==     ORTHOGONALIZE A SET OF WAVEFUNCTIONS C0                  ==
    ! ==     USING A DISTRIBUTED OVERLAP MATRIX ALGORITHM             ==
    ! ==     Cache Optimized Blocked Grahm-Schmidt                    ==
    ! ==     C.Bekas and A.Curioni Comp. Phys. Comm. (2010)           ==
    ! ==--------------------------------------------------------------==
    ! 
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(:,:)

    CHARACTER(len=*), PARAMETER :: procedureN = 'disortho_new_low'

    COMPLEX(real_8), ALLOCATABLE             :: C0_local(:,:), localw(:)
    INTEGER                                  :: ibeg_c0, iend_c0, ierr, &
                                                IGEQ0_local, index, index1, &
                                                isub, isub2, isub3, &
                                                NGWK_local, norbx, norbx1
    LOGICAL                                  :: debug, GEQ0_local
    REAL(real_8), ALLOCATABLE                :: global_h(:), local_h(:)

! ==--------------------------------------------------------------==

    IF (nstate.LE.0) RETURN
    CALL tiset(procedureN,isub)
    debug=.FALSE.

    ! >>>>>>> cp_grp trick
    CALL tiset(procedureN//'_grps_a',isub2)
    CALL cp_grp_get_sizes(ngwk_l=NGWK_local,geq0_l=GEQ0_local,&
         igeq0_l=IGEQ0_local,firstk_g=ibeg_c0,lastk_g=iend_c0)
    ALLOCATE(C0_local(NGWK_local,nstate),stat=ierr)
    IF (ierr.NE.0) CALL stopgm(procedureN,'Allocation problem',& 
         __LINE__,__FILE__)
    CALL cp_grp_copy_wfn_to_local(c0,nkpt%ngwk,C0_local,NGWK_local,&
         ibeg_c0,NGWK_local,nstate)
    CALL tihalt(procedureN//'_grps_a',isub2)
    ! <<<<<<<

    IF(cnti%disortho_bsize /= 0 ) THEN
       norbx = MIN(cnti%disortho_bsize,nstate)
    ELSE
       CALL set_orbdist(nstate,cnti%nstblk,parai%cp_nproc,norbx)
    ENDIF
    norbx1 = norbx

    IF (norbx.GT.nstate) THEN
       CALL stopgm(procedureN,'BLOCK SIZE LARGER THAN # OF STATES',& 
            __LINE__,__FILE__)
    ENDIF

    ALLOCATE(localw((nstate-((nstate/norbx)-1)*norbx)&
         *NGWK_local),&
         local_h((nstate-((nstate/norbx)-1)*norbx)&
         *nstate),&
         global_h((nstate-((nstate/norbx)-1)*norbx)&
         *nstate),stat=ierr)
    IF (ierr/=0) CALL stopgm(procedureN,'Allocation problem',& 
         __LINE__,__FILE__)
    CALL zeroing(local_h)!, (nstate-((nstate/norbx)-1)*norbx)*nstate)
    CALL zeroing(global_h)!, (nstate-((nstate/norbx)-1)*norbx)*nstate)
    CALL zeroing(localw)!, (nstate-((nstate/norbx)-1)*norbx)*NGWK_local)

    DO index = 1, CEILING(REAL(nstate,real_8)/REAL(norbx,real_8))
       ! COPY LOCAL BLOCK TO ORTHOGONALIZE
       CALL dcopy(2*NGWK_local*norbx1,C0_local(1,(index-1)*norbx+1),1,&
            localw, 1)
       ! COMPUTE PROJECTIONS
       IF (index.GT.1) THEN
          IF (NGWK_local>0) THEN
             CALL dgemm("T", "N", (index-1)*norbx, norbx1,&
                  2*NGWK_local,&
                  2.0_real_8, C0_local(1,1), 2*NGWK_local, localw,&
                  2*NGWK_local, 0.0_real_8, local_h, norbx*(index-1))
          ENDIF
          IF (GEQ0_local) THEN
             CALL dger((index-1)*norbx, norbx1, -1.0_real_8, C0_local(1,1),&
                  2*NGWK_local, localw, 2*ngwk_local, local_h,&
                  (index-1)*norbx)
          ENDIF
          CALL mp_sum(local_h, global_h,norbx1*norbx*(index-1), parai%cp_grp)
          IF (NGWK_local>0) THEN
             CALL dgemm("N", "N", 2*NGWK_local, norbx1,&
                  norbx*(index-1),&
                  -1.0_real_8, C0_local(1,1),&
                  2*NGWK_local, global_h, norbx*(index-1),&
                  1.0_real_8, localw, 2*NGWK_local)
          ENDIF
       ENDIF
       ! ORTHOGONALIZE W AND STORE IT BACK TO C0
       IF (NGWK_local>0) THEN
          CALL dsyrk("U", "T", norbx1, 2*NGWK_local, 2._real_8,&
               localw, 2*NGWK_local, 0._real_8,&
               local_h, norbx1)
       ENDIF
       IF (GEQ0_local) THEN
          CALL dger(norbx1, norbx1, -1.0_real_8, localw, 2*NGWK_local,&
               localw, 2*NGWK_local, local_h, norbx1)
       ENDIF
       CALL dmatc("U", norbx1, local_h, norbx1)
       CALL mp_sum(local_h, global_h, norbx1*norbx1, parai%cp_grp)
       ! CALCULATE CHOLESKY FACTORIZATION
       CALL dpotrf("U", norbx1, global_h, norbx1, ierr)
       IF (ierr.LT.0) THEN
          IF (paral%io_parent)WRITE(6,*) "WRONG ENTRY IN CHOL:", -ierr
          CALL stopgm(procedureN,'DPOTRF: WRONG ENTRY',& 
               __LINE__,__FILE__)
       ENDIF
       IF (ierr.GT.0) THEN
          ! LAST BLOCK IS LINEARLY DEPENDENT. USE MOD. GRAM SCHMIDT.
          CALL gsortho(localw,NGWK_local,1,norbx1)
          CALL dcopy(2*NGWK_local*norbx1, localw, 1,&
               C0_local(1,(index-1)*norbx+1), 1)
          GOTO 999
       ENDIF
       DO index1=1, norbx1*norbx1
          local_h(index1) = 0._real_8
       ENDDO
       DO index1=1, norbx1
          local_h((index1-1)*norbx1+index1) = 1._real_8
       ENDDO
       ! CALCULATE INVERSE OF CHOLESKY FACTOR BY SOLVING
       DO index1=1, norbx1
          CALL dtrsv("U", "N", "N", norbx1, global_h, norbx1,&
               local_h(norbx1*(index1-1)+1), 1)
       ENDDO
       ! MULTIPLY BACK TO ORIGINAL MATRIX
       IF (NGWK_local>0) THEN
          CALL dgemm("N", "N", 2*NGWK_local, norbx1, norbx1, 1.0_real_8,&
               localw, 2*NGWK_local, local_h, norbx1, 0._real_8,&
               C0_local(1,(index-1)*norbx+1), 2*NGWK_local)
       ENDIF
999    IF (index.EQ.(CEILING(REAL(nstate,real_8)/REAL(norbx,real_8))-1)) THEN
          norbx1 = nstate - index*norbx
       ENDIF
    ENDDO

    DEALLOCATE(localw,local_h,global_h,stat=ierr)
    IF (ierr/=0) CALL stopgm(procedureN,'Deallocation problem',& 
         __LINE__,__FILE__)

    ! >>>>>>> cp_grp trick
    CALL tiset(procedureN//'_grps_b',isub3)
    ! we need to zero the C0 that we can do the reduce
    CALL cp_grp_copy_local_to_wfn(C0_local,NGWK_local,c0,nkpt%ngwk,&
         ibeg_c0,NGWK_local,nstate)
    CALL cp_grp_zero_g(c0,nkpt%ngwk,nstate,ibeg_c0,iend_c0)
    CALL cp_grp_redist(c0,nkpt%ngwk,nstate)
    DEALLOCATE(C0_local,stat=ierr)
    IF (ierr.NE.0) CALL stopgm(procedureN,'Deallocation problem',& 
         __LINE__,__FILE__)
    CALL tihalt(procedureN//'_grps_b',isub3)
    ! <<<<<<<

    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE disortho_new_low


  ! ==================================================================
  SUBROUTINE disortho_new_cuda_low(nstate,c0)
    ! ==--------------------------------------------------------------==
    ! ==     ORTHOGONALIZE A SET OF WAVEFUNCTIONS C0                  ==
    ! ==     USING A DISTRIBUTED OVERLAP MATRIX ALGORITHM             ==
    ! ==     Cache Optimized Blocked Grahm-Schmidt                    ==
    ! ==     C.Bekas and A.Curioni Comp. Phys. Comm. (2010)           ==
    ! ==--------------------------------------------------------------==    ! 
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(:,:)

    CHARACTER(len=*), PARAMETER :: procedureN = 'disortho_new_cuda_low'
    LOGICAL, PARAMETER                       :: USE_DSYRK = .FALSE.

    COMPLEX&
      (real_8) __ASYNCHRONOUS __CONTIGUOUS, &
      DIMENSION(:, :), POINTER               :: C0_local
    COMPLEX(real_8), DIMENSION(:), &
      POINTER __CONTIGUOUS                   :: C0_local_1d_p
    INTEGER :: beg_pos, device_idx, end_pos, i_g0_l, i_stream, ibeg_c0, &
      ibeg_c0_1d, ibeg_c0_loc_d, iend_c0, iend_c0_loc_d, ierr, IGEQ0_local, &
      index, isub, isub2, isub3, j_stream, ldc0_loc_d, lwork, m_c0_loc_d, &
      n_max_threads, n_nested_threads, n_streams_per_task, NGWK_local, norbx, &
      norbx1, stream_idx
    LOGICAL                                  :: debug, GEQ0_loc_d, GEQ0_local
    REAL(real_8), ALLOCATABLE                :: global_h(:), local_h(:,:)
    TYPE(cublas_handle_t), POINTER           :: blas_handle
    TYPE(cuda_event_t), POINTER              :: event_copy_d2h, event_copy_h2d
    TYPE(cuda_memory_t), POINTER             :: C0_local_d, global_h_d, &
                                                local_h_d, localw_d, work_d
    TYPE(cuda_stream_t), POINTER             :: stream_1, stream_2, stream_3
    TYPE(cusolver_handle_t), POINTER         :: solver_handle
    TYPE(thread_view_t)                      :: thread_view
    TYPE(thread_view_t), ALLOCATABLE, &
      DIMENSION(:)                           :: thread_views

    !$ LOGICAL :: nested_orig
    !$ INTEGER :: max_level_orig
    ! ==--------------------------------------------------------------==

    IF (nstate.LE.0) RETURN
    CALL tiset(procedureN,isub)
    __NVTX_TIMER_START( 'disortho' )
    __NVTX_TIMER_START( 'disortho_a' )

    debug=.FALSE.

    IF (paral%io_parent)WRITE(6,*) procedureN,' USE_DSYRK',USE_DSYRK

    ! >>>>>>> cp_grp trick
    CALL tiset(procedureN//'_grps_a',isub2)
    CALL cp_grp_get_sizes(ngwk_l=NGWK_local,geq0_l=GEQ0_local,&
         igeq0_l=IGEQ0_local,firstk_g=ibeg_c0,lastk_g=iend_c0,i_g0_l=i_g0_l)

    CALL cuda_alloc_host ( C0_local, [NGWK_local,nstate] )

    CALL cp_grp_copy_wfn_to_local(c0,nkpt%ngwk,C0_local,NGWK_local,&
         ibeg_c0,NGWK_local,nstate)
    CALL tihalt(procedureN//'_grps_a',isub2)
    ! <<<<<<<

    !vw n_streams_per_task is by default 1, can be changed in the input
    n_streams_per_task = MAX( cp_cuda_env%blas_n_streams_per_device * cp_cuda_env%blas_n_devices_per_task, 1 )!vw fix the min
    IF( .NOT. cp_cuda_env%use_blas ) n_streams_per_task = 1

    n_max_threads = 1
    !$ n_max_threads = omp_get_max_threads()
    IF( n_streams_per_task > n_max_threads ) CALL stopgm(procedureN,&
         'number of streams creater than max number of threads',&
         __LINE__,__FILE__)

    n_nested_threads = n_max_threads / n_streams_per_task

    CALL thread_views_init ( cp_cuda_env%blas_n_devices_per_task, cp_cuda_env%blas_n_streams_per_device, thread_views )


    IF( cp_cuda_env%use_blas ) THEN
#ifndef _HASNT_OMP_SET_NESTED
       !$ nested_orig = omp_get_nested ( )
#endif
       !$ max_level_orig = omp_get_max_active_levels( )
#ifndef _HASNT_OMP_SET_NESTED
       !$ CALL omp_set_nested( .TRUE. )
#endif
       !$ CALL omp_set_max_active_levels( 2 )
    ENDIF


    IF(cnti%disortho_bsize /= 0 ) THEN
       norbx = MIN(cnti%disortho_bsize,nstate)
    ELSE
       CALL set_orbdist(nstate,cnti%nstblk,parai%cp_nproc,norbx)
    ENDIF
    norbx1 = norbx

    IF (norbx.GT.nstate) THEN
       CALL stopgm(procedureN,'BLOCK SIZE LARGER THAN # OF STATES',& 
            __LINE__,__FILE__)
    ENDIF

    ALLOCATE(local_h((nstate-((nstate/norbx)-1)*norbx)*nstate, 0:n_streams_per_task-1),&
         global_h((nstate-((nstate/norbx)-1)*norbx)*nstate),stat=ierr)
    IF (ierr/=0) CALL stopgm(procedureN,'Allocation problem',& 
         __LINE__,__FILE__)
    CALL zeroing(local_h)
    CALL zeroing(global_h)


    !>
    !=
    !vw get max number of rows of C0 that needs to be on streams/devices
    !max_m_c0_loc_d = 0
    !DO i_stream = 0, n_streams_per_task - 1
    !   CALL part_1d_get_blk_bounds( NGWK_local, i_stream,n_streams_per_task, ibeg_c0_loc_d, iend_c0_loc_d )
    !   m_c0_loc_d = iend_c0_loc_d - ibeg_c0_loc_d + 1
    !   max_m_c0_loc_d = MAX( max_m_c0_loc_d, m_c0_loc_d )
    !ENDDO
    !
    !CALL cp_cuortho_init ( cp_cuda_env, cp_cuda_devices_blas, cp_cuortho )
    !CALL cp_cuortho_alloc_buffers ( cp_cuortho, max_m_c0_loc_d, nstate, norbx, ldc0_loc_d )
    CALL cp_cuortho_get ( cp_cuortho, ldc0l_d=ldc0_loc_d )

    __NVTX_TIMER_STOP

    !vw multi streams / devices
    !$omp parallel   if( cp_cuda_env%use_blas ) &
    !$omp            num_threads( n_streams_per_task ) &
    !$omp            default( none ) &
    !$omp            private( C0_local_d, localw_d, local_h_d, global_h_d, &
    !$omp                     work_d, lwork, stream_1, stream_2, stream_3, blas_handle, &
    !$omp                     solver_handle, &
    !$omp                     event_copy_h2d, event_copy_d2h, index, beg_pos, end_pos, &
    !$omp                     ierr, i_stream, j_stream, ibeg_c0_loc_d, ibeg_c0_1d, &
    !$omp                     iend_c0_loc_d, m_c0_loc_d, C0_local_1d_p, thread_view, &
    !$omp                     device_idx, stream_idx, GEQ0_loc_d ) &
    !$omp            firstprivate( norbx1 ) &
    !$omp            shared( cp_cuortho, parai, c0_local, local_h, global_h, thread_views, &
    !$omp                    NGWK_local, nstate, norbx, n_streams_per_task, &
    !$omp                    n_nested_threads, cp_cuda_env, ldc0_loc_d, i_g0_l )

    !vw set number of children threads
    !$ CALL omp_set_num_threads ( n_nested_threads )

    i_stream = 0
    !$ i_stream = omp_get_thread_num()
    IF(cp_cuda_env%blas_n_devices_per_task>0) THEN
       thread_view = thread_views( i_stream ) !vw need to fix that
       device_idx = thread_view%device_idx
       stream_idx = thread_view%stream_idx - 1 !vw FIX that -1
    ELSE
       device_idx = 1
       stream_idx = 0
    ENDIF

#if defined( _HASNT_F08_POINTER_REMAPPING )
    CALL stopgm(procedureN,'compiler needs to support pointer remapping!',&
         __LINE__,__FILE__)
#else
    C0_local_1d_p(1:SIZE(C0_local)) => C0_local
#endif

    CALL part_1d_get_blk_bounds( NGWK_local, i_stream, n_streams_per_task, ibeg_c0_loc_d, iend_c0_loc_d )
    m_c0_loc_d = iend_c0_loc_d - ibeg_c0_loc_d + 1
    GEQ0_loc_d = i_g0_l >= ibeg_c0_loc_d .AND. i_g0_l <= iend_c0_loc_d

    CALL cp_cuortho_get_buffer_ptrs ( cp_cuortho, device_idx, stream_idx, C0_local_d=C0_local_d, &
         localw_d=localw_d, local_h_d=local_h_d, global_h_d=global_h_d, &
         work_d=work_d, lwork=lwork, stream_1=stream_1, stream_2=stream_2, stream_3=stream_3, &
         blas_handle=blas_handle, solver_handle=solver_handle, event_copy_h2d=event_copy_h2d, &
         event_copy_d2h=event_copy_d2h )

    !vw> this isnt needed
    CALL cuda_mem_zero_bytes ( C0_local_d, C0_local_d%n_bytes )
    CALL cuda_mem_zero_bytes ( localw_d, localw_d%n_bytes )
    CALL cuda_mem_zero_bytes ( local_h_d, local_h_d%n_bytes )
    CALL cuda_mem_zero_bytes ( global_h_d, global_h_d%n_bytes )

    !
    !vw allocate cuda memory
    !vw copy to GPU C0_local (perhaps async, first chunk outside loop, then index+1 for next cycle)
    ibeg_c0_1d = ibeg_c0_loc_d
    CALL cuda_memcpy2d_async_host_to_device ( C0_local_1d_p(ibeg_c0_1d:), NGWK_local, &
         & C0_local_d, ldc0_loc_d, m_c0_loc_d, norbx1, stream_1 )
    !<

    DO index = 1, CEILING(REAL(nstate,real_8)/REAL(norbx,real_8))
       ! COPY LOCAL BLOCK TO ORTHOGONALIZE

       ! record an event on stream_1, so stream_3 waits for completion of the copy
       CALL cuda_event_record ( event_copy_h2d, stream_1 )
       ! stream_3 waits that stream_1 is done with copy
       CALL cuda_stream_wait_event ( stream_3, event_copy_h2d )
       IF( index < CEILING(REAL(nstate,real_8)/REAL(norbx,real_8)) ) THEN
          beg_pos = index*norbx+1
          end_pos = beg_pos + norbx - 1
          IF( index == CEILING(REAL(nstate,real_8)/REAL(norbx,real_8)) - 1 ) end_pos = nstate
          ibeg_c0_1d = ibeg_c0_loc_d + ( beg_pos - 1 ) * NGWK_local
          CALL cuda_memcpy2d_async_host_to_device ( C0_local_1d_p(ibeg_c0_1d:), NGWK_local, &
               & cuda_z_points_to( C0_local_d, ldc0_loc_d*index*norbx+1 ), ldc0_loc_d, &
               & m_c0_loc_d, end_pos-beg_pos+1, stream_1 )
       ENDIF

       !vw copy C0_local to localw
       CALL cublas_zcopy ( blas_handle, ldc0_loc_d*norbx1, &
            & cuda_z_points_to( C0_local_d, ldc0_loc_d*(index-1)*norbx+1 ), 1, localw_d, 1 )

       ! COMPUTE PROJECTIONS
       IF (index.GT.1) THEN
          IF (m_c0_loc_d>0) THEN
             CALL cublas_dgemm ( blas_handle, 'T', 'N', (index-1)*norbx, norbx1, 2*m_c0_loc_d, &
                  & 2.0_real_8, C0_local_d, 2*ldc0_loc_d, &
                  & localw_d, 2*ldc0_loc_d, &
                  & 0.0_real_8, local_h_d, norbx*(index-1) )
          ENDIF
          IF (GEQ0_loc_d) THEN
             CALL cublas_dger ( blas_handle, (index-1)*norbx, norbx1, -1.0_real_8, &
                  & C0_local_d, 2*ldc0_loc_d, localw_d, 2*ldc0_loc_d, local_h_d, (index-1)*norbx )
          ENDIF

          ! copy to PROC local_h
          CALL cuda_memcpy_async_device_to_host ( local_h_d, local_h(1:norbx1*norbx*(index-1),i_stream), stream_3 )

          CALL cuda_stream_synchronize ( stream_3 )

          !vw sync the threads and master accumulate the data
          __NVTX_TIMER_START( 'disortho_b' )
          !$omp barrier
          !$omp master
          DO j_stream = 1, n_streams_per_task - 1
             local_h(1:norbx1*norbx*(index-1), 0) = local_h(1:norbx1*norbx*(index-1), 0) + &
                  local_h(1:norbx1*norbx*(index-1), j_stream )
          ENDDO
          CALL mp_sum(local_h(:,0), global_h,norbx1*norbx*(index-1), parai%cp_grp)
          !$omp end master
          !$omp barrier
          __NVTX_TIMER_STOP

          !vw copy to GPU global_h
          CALL cuda_memcpy_async_host_to_device ( global_h(1:norbx1*norbx*(index-1)), global_h_d, stream_3 )
          CALL cuda_stream_synchronize ( stream_3 )

          IF (m_c0_loc_d>0) THEN
             CALL cublas_dgemm ( blas_handle, 'N', 'N', 2*m_c0_loc_d, norbx1, norbx*(index-1), &
                  & -1.0_real_8, C0_local_d, 2*ldc0_loc_d, &
                  & global_h_d, norbx*(index-1), &
                  & 1.0_real_8, localw_d, 2*ldc0_loc_d )
          ENDIF
       ENDIF

       ! ORTHOGONALIZE W AND STORE IT BACK TO C0
       IF (m_c0_loc_d>0) THEN
          !vw cublas_dsyrk performs poorly, we should rewrite it at some point
          IF( USE_DSYRK ) THEN
             CALL cublas_dsyrk ( blas_handle, 'U', 'T', norbx1, 2*m_c0_loc_d, &
                  & 2._real_8, localw_d, 2*ldc0_loc_d, &
                  & 0._real_8, local_h_d, norbx1 )
          ELSE
             CALL cublas_dgemm ( blas_handle, 'T','N', norbx1, norbx1, 2*m_c0_loc_d, &
                  & 2._real_8, localw_d, 2*ldc0_loc_d, localw_d, 2*ldc0_loc_d, &
                  & 0._real_8, local_h_d, norbx1 )
          ENDIF
       ENDIF

       IF (GEQ0_loc_d) THEN
          CALL cublas_dger ( blas_handle, norbx1, norbx1, -1.0_real_8, &
               & localw_d, 2*ldc0_loc_d, localw_d, 2*ldc0_loc_d, local_h_d, norbx1 )
       ENDIF

       ! copy to PROC local_h
       CALL cuda_memcpy_async_device_to_host ( local_h_d, local_h(1:norbx1**2,i_stream), stream_3 )
       CALL cuda_stream_synchronize ( stream_3 )

       !vw sync the threads and master accumulate the data
       __NVTX_TIMER_START( 'disortho_c' )
       !$omp barrier
       !$omp master
       DO j_stream = 1, n_streams_per_task - 1
          local_h(1:norbx1**2, 0) = local_h(1:norbx1**2, 0) + local_h(1:norbx1**2, j_stream )
       ENDDO
       CALL dmatc("U", norbx1, local_h, norbx1)
       CALL mp_sum(local_h(:,0), global_h, norbx1*norbx1, parai%cp_grp)

       ! CALCULATE CHOLESKY FACTORIZATION
       !>
       CALL dpotrf("U", norbx1, global_h, norbx1, ierr)
       IF (ierr.LT.0) CALL stopgm(procedureN,'DPOTRF: WRONG ENTRY IN CHOL:'//int2str(-ierr),& 
            __LINE__,__FILE__)
       !$omp end master
       !$omp barrier
       __NVTX_TIMER_STOP

       CALL cuda_memcpy_async_host_to_device ( global_h(1:norbx1**2), global_h_d, stream_3 )
       !=
       !vw copy to GPU global_h
       !   CALL cuda_memcpy_async_host_to_device ( global_h(1:norbx1**2), global_h_d, stream_3 )
       !vw doesnt work so far
       !CALL cusolver_dpotrf ( solver_handle, 'U', norbx1, global_h_d, norbx1, work_d, lwork, ierr )
       !<
       !IF (ierr.LT.0) CALL stopgm(procedureN,'DPOTRF: WRONG ENTRY IN CHOL:'//int2str(-ierr),& 
       !     __LINE__,__FILE__)
       !>
       !DO index1=1, norbx1*norbx1
       !   local_h(index1) = 0._real_8
       !ENDDO
       !DO index1=1, norbx1
       !   local_h((index1-1)*norbx1+index1) = 1._real_8
       !ENDDO
       !! CALCULATE INVERSE OF CHOLESKY FACTOR BY SOLVING
       !!DO index1=1, norbx1
       !!   CALL dtrsv("U", "N", "N", norbx1, global_h, norbx1, local_h(norbx1*(index1-1)+1), 1)
       !!ENDDO
       !CALL dtrsm ( 'L', 'U', 'N', 'N', norbx1, norbx1, 1.0_real_8, global_h, norbx1, local_h, norbx1 )
       !
       !   CALL cuda_memcpy_async_host_to_device ( local_h(1:norbx1**2), local_h_d, stream_3 )
       !=
       CALL cuuser_identity ( local_h_d, norbx1, stream_3 )
       CALL cublas_dtrsm ( blas_handle, 'L', 'U', 'N', 'N', norbx1, norbx1, &
            & 1.0_real_8, global_h_d, norbx1, local_h_d, norbx1  )
       !<

       ! MULTIPLY BACK TO ORIGINAL MATRIX
       IF (m_c0_loc_d>0) THEN
          CALL cublas_dgemm ( blas_handle, 'N', 'N', 2*m_c0_loc_d, norbx1, norbx1, &
               & 1.0_real_8, localw_d, 2*ldc0_loc_d, &
               & local_h_d, norbx1, &
               & 0.0_real_8, cuda_z_points_to( C0_local_d, ldc0_loc_d*(index-1)*norbx+1 ), 2*ldc0_loc_d )
       ENDIF

       beg_pos = (index-1)*norbx+1
       end_pos = beg_pos + norbx1 - 1
       CALL cuda_event_record ( event_copy_d2h, stream_3 )
       CALL cuda_stream_wait_event ( stream_2, event_copy_d2h )

       ibeg_c0_1d = ibeg_c0_loc_d + ( beg_pos - 1 ) * NGWK_local
       CALL cuda_memcpy2d_async_device_to_host ( cuda_z_points_to( C0_local_d, ldc0_loc_d*(index-1)*norbx+1 ), ldc0_loc_d, &
            C0_local_1d_p(ibeg_c0_1d:), NGWK_local, m_c0_loc_d, norbx1, stream_2 )

       IF (index.EQ.((CEILING(REAL(nstate,real_8)/REAL(norbx,real_8)))-1)) THEN
          norbx1 = nstate - index*norbx
       ENDIF

    ENDDO


    !vw copy to PROC C0_local
    CALL cuda_stream_synchronize ( stream_2 )

    !$omp end parallel

    __NVTX_TIMER_START( 'disortho_d' )

    CALL thread_views_finalize ( thread_views )

    IF( cp_cuda_env%use_blas ) THEN
#ifndef _HASNT_OMP_SET_NESTED
       !$ CALL omp_set_nested( nested_orig )
#endif
       !$ CALL omp_set_max_active_levels( max_level_orig )
    ENDIF

    DEALLOCATE(local_h,global_h,stat=ierr)
    IF (ierr/=0) CALL stopgm(procedureN,'Deallocation problem',& 
         __LINE__,__FILE__)


    ! >>>>>>> cp_grp trick
    CALL tiset(procedureN//'_grps_b',isub3)
    ! we need to zero the C0 that we can do the reduce
    CALL cp_grp_copy_local_to_wfn(C0_local,NGWK_local,c0,nkpt%ngwk,&
         ibeg_c0,NGWK_local,nstate)
    CALL cp_grp_zero_g(c0,nkpt%ngwk,nstate,ibeg_c0,iend_c0)
    CALL cp_grp_redist(c0,nkpt%ngwk,nstate)

    CALL cuda_dealloc_host ( C0_local )

    !CALL cp_cuortho_dealloc_buffers ( cp_cuortho )
    !CALL cp_cuortho_finalize ( cp_cuortho )

    CALL tihalt(procedureN//'_grps_b',isub3)
    ! <<<<<<<

    __NVTX_TIMER_STOP
    __NVTX_TIMER_STOP
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==

  END SUBROUTINE disortho_new_cuda_low


  ! ==================================================================
  SUBROUTINE disortho_old(nstate,c0,cscr)
    ! ==--------------------------------------------------------------==
    ! ==     ORTHOGONALIZE A SET OF WAVEFUNCTIONS C0                  ==
    ! ==     USING A DISTRIBUTED OVERLAP MATRIX ALGORITHM             ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(:,:), cscr(ncpw%ngw,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'disortho_old'

    INTEGER                                  :: i, i1, i2, ierr, ij, info, &
                                                ip, ipi, ipj, isub, j, j1, &
                                                mmx, n1, n2, norb, norbx, &
                                                npi, npj, nx
    LOGICAL                                  :: debug
    REAL(real_8), ALLOCATABLE                :: smat(:,:), tmat(:,:)
    REAL(real_8), ALLOCATABLE, TARGET        :: rmat(:)
    REAL(real_8), POINTER                    :: rmat_ptr(:,:)

    debug=.FALSE.
    IF (nstate.LE.0) RETURN
    CALL tiset(procedureN,isub)
    CALL set_orbdist(nstate,cnti%nstblk,parai%nproc,norbx)
    norb=paraw%nwa12(parai%mepos,2)-paraw%nwa12(parai%mepos,1)+1
    norb=MAX(0,norb)
    ! ..overlap matrix
    mmx=MAX(nstate*norbx,norbx*norbx*parai%nproc)
    ALLOCATE(smat(nstate, norbx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(tmat(nstate, norbx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(rmat(mmx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    CALL zeroing(smat)!,nstate*norb)
    CALL zeroing(rmat)!,nstate*norb)
    rmat_ptr(1:nstate,1:norbx)=>rmat
    DO ip=0,parai%nproc-1
       nx=paraw%nwa12(ip,2)-paraw%nwa12(ip,1)+1
       CALL zeroing(tmat)!,nstate*norbx)
       IF (nx.GT.0) THEN
          CALL ovlap2(ncpw%ngw,nstate,nx,tmat,c0,c0(1,paraw%nwa12(ip,1)),.TRUE.)
          CALL mp_sum(tmat,rmat_ptr,nstate*norbx,parap%pgroup(ip+1),parai%allgrp)
          IF (parai%mepos.EQ.parap%pgroup(ip+1)) THEN
             CALL dcopy(nstate*norbx,rmat,1,smat,1)
          ENDIF
       ENDIF
    ENDDO
    ! ..Cholesky decomposition
    DO ipj=0,parai%nproc-1
       npj=paraw%nwa12(ipj,2)-paraw%nwa12(ipj,1)+1
       j1=paraw%nwa12(ipj,1)
       IF (npj.GT.0) THEN
          DO ipi=ipj,parai%nproc-1
             npi=paraw%nwa12(ipi,2)-paraw%nwa12(ipi,1)+1
             i1=paraw%nwa12(ipi,1)
             IF (npi.GT.0) THEN
                IF (parai%mepos.LT.ipj .AND. norb.GT.0) THEN
                   CALL dgemm('N','T',npi,norb,npj,-1._real_8,smat(i1,1),nstate,&
                        smat(j1,1),nstate,0._real_8,rmat,npi)
                ELSE
                   CALL zeroing(rmat(1:npi*npj))!,npi*npj)
                ENDIF
                CALL mp_sum(rmat,npi*npj,parai%allgrp)
                IF (parai%mepos.EQ.ipj) THEN
                   DO j=1,npj
                      DO i=1,npi
                         ij=(j-1)*npi+i
                         rmat(ij)=rmat(ij)+smat(i1+i-1,j)
                      ENDDO
                   ENDDO
                   IF (ipj.EQ.ipi) THEN
                      CALL dpotrf('L',npi,rmat,npi,info)
                   ELSE
                      CALL dtrsm('R','L','T','N',npi,npj,1._real_8,smat(j1,1),&
                           nstate,rmat,npi)
                   ENDIF
                   DO j=1,npj
                      DO i=1,npi
                         ij=(j-1)*npi+i
                         smat(i1+i-1,j)=rmat(ij)
                      ENDDO
                   ENDDO
                ENDIF
                CALL mp_sync(parai%allgrp)
             ENDIF
          ENDDO
       ENDIF
    ENDDO
    ! ..Inversion of Cholesky-Matrix
    CALL zeroing(tmat)!,nstate*norbx)
    ! ....on site block
    DO i=paraw%nwa12(parai%mepos,1),paraw%nwa12(parai%mepos,2)
       i1=i-paraw%nwa12(parai%mepos,1)+1
       tmat(i,i1)=1._real_8
    ENDDO
    i1=paraw%nwa12(parai%mepos,1)
    IF (norb.GT.0)&
         CALL dtrtrs('L','N','N',norb,norb,smat(i1,1),nstate,&
         tmat(i1,1),nstate,info)
    IF (paraw%nwa12(parai%mepos,2)+1.LE.nstate) THEN
       i2=paraw%nwa12(parai%mepos,2)+1
       n2=nstate-paraw%nwa12(parai%mepos,2)
       IF (norb.GT.0)&
            CALL dgemm('N','N',n2,norb,norb,-1._real_8,smat(i2,1),nstate,&
            tmat(i1,1),nstate,0._real_8,tmat(i2,1),nstate)
    ENDIF
    ! ....other blocks
    CALL dcopy(nstate*norbx,smat,1,rmat,1)
    DO ip=1,parai%nproc-1
       n1=paraw%nwa12(ip,2)-paraw%nwa12(ip,1)+1
       i1=paraw%nwa12(ip,1)
       IF (n1.GT.0) THEN
          CALL dcopy(nstate*norbx,rmat,1,smat,1)
          CALL mp_bcast(smat,SIZE(smat),parap%pgroup(ip+1),parai%allgrp)
       ENDIF
       IF (ip.GT.parai%mepos) THEN
          ! ......solve for block mepos+ip
          IF (n1.GT.0) THEN
             IF (norb.GT.0)&
                  CALL dtrtrs('L','N','N',n1,norb,smat(i1,1),nstate,&
                  tmat(i1,1),nstate,info)
             IF (paraw%nwa12(ip,2)+1.LE.nstate) THEN
                i2=paraw%nwa12(ip,2)+1
                n2=nstate-paraw%nwa12(ip,2)
                IF (norb.GT.0)&
                     CALL dgemm('N','N',n2,n1,norb,-1._real_8,smat(i2,1),&
                     nstate,tmat(i1,1),nstate,1._real_8,tmat(i2,1),nstate)
             ENDIF
          ENDIF
       ENDIF
    ENDDO
    ! ..Transpose TMAT
    CALL trans_da(tmat,nstate,paraw%nwa12(0,1),paraw%nwa12(0,2),&
         norbx,parai%mepos,parai%nproc,parai%allgrp)
    ! ..Rotate orbitals
    CALL rotate_da(1._real_8,c0,0._real_8,cscr,tmat,2*ncpw%ngw,2*ncpw%ngw,nstate,&
         paraw%nwa12(0,1),paraw%nwa12(0,2),norbx,parai%mepos,parap%pgroup,parai%nproc,parai%allgrp,&
         .FALSE.,1,1)
    CALL dcopy(2*nstate*ncpw%ngw,cscr,1,c0,1)
    ! ==--------------------------------------------------------------==
    DEALLOCATE(smat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(tmat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(rmat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
  END SUBROUTINE disortho_old
  ! ==================================================================

END MODULE disortho_utils
