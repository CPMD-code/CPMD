#include "cpmd_global.h"

MODULE rhoofr_utils
  USE augchg_utils,                    ONLY: augchg
  USE cnst,                            ONLY: uimag
  USE cp_cuda_types,                   ONLY: cp_cuda_env
  USE cp_cudensity_utils,              ONLY: cp_cubuild_density_copy_to_host,&
                                             cp_cubuild_density_sum
  USE cp_cufft_types,                  ONLY: cp_cufft,&
                                             cp_cufft_device_get_ptrs,&
                                             cp_cufft_stream_get_ptrs
  USE cp_curho_types,                  ONLY: cp_curho,&
                                             cp_curho_stream_get_ptrs
  USE cp_curho_utils,                  ONLY: cp_curho_alloc_buffers,&
                                             cp_curho_dealloc_buffers
  USE cp_cuwfn_types,                  ONLY: cp_cuwfn,&
                                             cp_cuwfn_device_get_ptrs,&
                                             cp_cuwfn_is_init,&
                                             cp_cuwfn_put
  USE cp_grp_utils,                    ONLY: cp_grp_redist
  USE cppt,                            ONLY: indz,&
                                             nzh
  USE cuda_types,                      ONLY: cuda_memory_t,&
                                             cuda_stream_t
  USE cuda_utils,                      ONLY: cuda_alloc_host,&
                                             cuda_dealloc_host,&
                                             cuda_mem_zero_bytes,&
                                             cuda_memcpy_host_to_device,&
                                             cuda_stream_synchronize,&
                                             cuda_z_points_to
  USE cuuser_utils,                    ONLY: CuUser_setpsi_1_state_g,&
                                             CuUser_setpsi_2_states_g
  USE density_utils,                   ONLY: build_density_imag,&
                                             build_density_real,&
                                             build_density_sum
  USE dg,                              ONLY: ipooldg,&
                                             tdgcomm
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: chrg
  USE error_handling,                  ONLY: stopgm
  USE fft,                             ONLY: inzf,&
                                             jgw,&
                                             jhg,&
                                             llr1,&
                                             lr1s,&
                                             lr2s,&
                                             lr3s,&
                                             nzff
  USE fftmain_utils,                   ONLY: fwfftn,&
                                             invfftn
  USE fftnew_utils,                    ONLY: setfftn
  USE geq0mod,                         ONLY: geq0
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kin_energy_utils,                ONLY: kin_energy
  USE kinds,                           ONLY: real_8
  USE moverho_utils,                   ONLY: give_scr_moverho,&
                                             moverho
  USE mp_interface,                    ONLY: mp_comm_dup,&
                                             mp_comm_free,&
                                             mp_comm_null,&
                                             mp_sum
  USE parac,                           ONLY: parai,&
                                             paral
  USE part_1d,                         ONLY: part_1d_get_el_in_blk,&
                                             part_1d_nbr_el_in_blk
  USE pslo,                            ONLY: pslo_com
  USE rho1ofr_utils,                   ONLY: rhoabofr
  USE rhov_utils,                      ONLY: give_scr_rhov,&
                                             rhov
  USE ropt,                            ONLY: ropt_mod
  USE rswfmod,                         ONLY: maxstates,&
                                             rsactive,&
                                             rswf
  USE sfac,                            ONLY: fnl
  USE spin,                            ONLY: clsd,&
                                             lspin2,&
                                             spin_mod
  USE state_utils,                     ONLY: set_psi_1_state_g,&
                                             set_psi_2_states_g
  USE symm,                            ONLY: symmi
  USE symtrz_utils,                    ONLY: give_scr_symrho
  USE system,                          ONLY: cntl,&
                                             dual00,&
                                             fpar,&
                                             group,&
                                             ncpw,&
                                             parm,&
                                             spar
  USE tauofr_utils,                    ONLY: tauofr
  USE thread_view_types,               ONLY: thread_view_get,&
                                             thread_view_t
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

  PUBLIC :: rhoofr
  PUBLIC :: give_scr_rhoofr
  !public :: movepsih

CONTAINS

  ! ==================================================================
  SUBROUTINE rhoofr(c0,rhoe,psi,nstate)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==  THE NORMALIZED ELECTRON DENSITY RHOE IN REAL SPACE          ==
    ! ==  THE KINETIC ENERGY EKIN. IT IS DONE IN RECIPROCAL SPACE     ==
    ! ==  WHERE THE ASSOCIATED OPERATORS ARE DIAGONAL.                ==
    ! ==  RHOE IS OBTAINED FOURIER TRANSFORMING THE WFN TO REAL       ==
    ! ==  SPACE (PSI).                                                ==
    ! ==--------------------------------------------------------------==
    ! == WARNING: ALL WAVEFUNCTIONS C0 HAVE TO BE ORTHOGONAL          ==
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8)                          :: c0(:,:)
    REAL(real_8), TARGET __CONTIGUOUS        :: rhoe(:,:)
    COMPLEX(real_8), TARGET __CONTIGUOUS     :: psi(:)
    INTEGER                                  :: nstate

    CHARACTER(*), PARAMETER                  :: procedureN = 'rhoofr'
    COMPLEX(real_8), PARAMETER               :: zone = (1.0_real_8,0.0_real_8)
    REAL(real_8), PARAMETER :: delta = 1.e-6_real_8, &
      o3 = 0.33333333333333333_real_8

    COMPLEX(real_8), DIMENSION(:), &
      POINTER __CONTIGUOUS                   :: psi_p
    COMPLEX(real_8), DIMENSION(:, :), &
      POINTER __CONTIGUOUS                   :: psis
    INTEGER :: device_idx, fft_comm, i, i_stream, ia, iat, ierr, ir, is, is1, &
      is2, ispin1, ispin2, isub, isub2, isub3, isub4, iwf, n_max_threads, &
      n_nested_threads, n_streams_per_task, stream_idx
    LOGICAL                                  :: copy_data_to_device, &
                                                copy_data_to_host, tfcal
    REAL(real_8)                             :: chksum, coef3, coef4, ral, &
                                                rbe, rsp, rsum, rsum1, &
                                                rsum1abs, rsumv, rto
    REAL(real_8), ALLOCATABLE                :: qa(:)
    REAL(real_8), ALLOCATABLE, &
      DIMENSION(:, :, :), TARGET             :: rhoes
    REAL(real_8), DIMENSION(:, :), &
      POINTER __CONTIGUOUS                   :: rhoe_p
    TYPE(cuda_memory_t), POINTER             :: c0_d, inzs_d, nzfs_d, psi_d, &
                                                rho_d
    TYPE(cuda_stream_t), POINTER             :: stream_p
    TYPE(thread_view_t)                      :: thread_view
    TYPE(thread_view_t), ALLOCATABLE, &
      DIMENSION(:)                           :: thread_views

    !$ LOGICAL :: nested_orig
    !$ INTEGER :: max_level_orig
    LOGICAL, PARAMETER :: wfn_on_gpu = .FALSE.
    ! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)

    __NVTX_TIMER_START ( procedureN )
    __NVTX_TIMER_START ( procedureN//'_a' )
    ! ==--------------------------------------------------------------==
    CALL kin_energy(c0,nstate,rsum)

    ! ==--------------------------------------------------------------==
    ! CASPUR 2/5/04
    ! Initialize FFT datastructure
    IF (group%nogrp.GT.1)CALL stopgm(procedureN,&
         'OLD TASK GROUPS NOT SUPPORTED ANYMORE ',&
         __LINE__,__FILE__)

    IF (tdgcomm%tdg) THEN
       IF (lspin2%tlse)CALL stopgm(procedureN,&
            'TDG TOGETHER WITH LSE IS NOT YET IMPLEMENTED ',&
            __LINE__,__FILE__)
       CALL setfftn(ipooldg)
    ELSE
       CALL setfftn(0)
    ENDIF
    ! ==--------------------------------------------------------------==

    ! Initialize
    CALL zeroing(rhoe)!,clsd%nlsd*nnr1)

    ! Loop over the electronic states

    !vw gpu version:
    !vw if possible do the accumulation of the density on the gpu
    copy_data_to_host = rsactive .OR. cntl%tlsd .OR. lspin2%tlse .OR. ( .NOT. cp_cuda_env%use_fft )
    copy_data_to_device = .NOT. ( wfn_on_gpu .AND. cp_cuda_env%use_fft .AND. cp_cuwfn_is_init ( cp_cuwfn ) )
    !write(*,*) 'rhoofr_utils.mod.F90: copy_data_to_device=',copy_data_to_device
    !vw n_streams_per_task is by default 1, can be changed in the input
    n_streams_per_task = MAX( cp_cuda_env%fft_n_streams_per_device * cp_cuda_env%fft_n_devices_per_task, 1 )!vw fix the min
    IF( .NOT. cp_cuda_env%use_fft ) n_streams_per_task = 1

    n_max_threads = 1
    !$ n_max_threads = omp_get_max_threads()
    IF( n_streams_per_task > n_max_threads ) CALL stopgm(procedureN,&
         'number of streams creater than max number of threads',&
         __LINE__,__FILE__)

    n_nested_threads = n_max_threads / n_streams_per_task

    fft_comm = mp_comm_null
    CALL thread_views_init ( cp_cuda_env%fft_n_devices_per_task, cp_cuda_env%fft_n_streams_per_device, thread_views )

    IF( cp_cuda_env%use_fft ) THEN
       CALL tiset(procedureN//'_gpu_init',isub2)
#ifndef _HASNT_OMP_SET_NESTED
       !$ nested_orig = omp_get_nested ( )
#endif
       !$ max_level_orig = omp_get_max_active_levels( )
#ifndef _HASNT_OMP_SET_NESTED
       !$ CALL omp_set_nested( .TRUE. )
#endif
       !$ CALL omp_set_max_active_levels( 2 )
       ALLOCATE(rhoes(SIZE(rhoe,1),SIZE(rhoe,2),0:n_streams_per_task-1))
       CALL cuda_alloc_host(psis,[SIZE(psi),n_streams_per_task],lbs=[1,0])
       CALL zeroing(rhoes)

       CALL cp_curho_alloc_buffers ( cp_curho, SIZE(rhoe) ) !vw for the moment only for RKS. need to pass the extra dim.

       !>>
       IF( .NOT. copy_data_to_device ) THEN
          DO device_idx = 1, cp_cuda_env%fft_n_devices_per_task
             CALL cp_cuwfn_device_get_ptrs ( cp_cuwfn, device_idx, c0_d=c0_d )
             CALL cuda_memcpy_host_to_device ( c0, c0_d )
          ENDDO
          chksum = SUM( ABS( c0(:,1) ) )
          CALL cp_cuwfn_put ( cp_cuwfn, ready=.TRUE., chksum=chksum  )
       ENDIF
       !<<
       CALL tihalt(procedureN//'_gpu_init',isub2)
    ENDIF

    __NVTX_TIMER_STOP

    !$omp parallel   if( cp_cuda_env%use_fft ) &
    !$omp            num_threads( n_streams_per_task ) &
    !$omp            default( none ) &
    !$omp            private( i, is1, is2, iwf, tfcal, i_stream, psi_p, rhoe_p, &
    !$omp                     coef3, coef4, ispin1, ispin2, thread_view, device_idx, stream_idx, &
    !$omp                     rho_d, stream_p, c0_d, nzfs_d, inzs_d, psi_d ) &
    !$omp            firstprivate( fft_comm ) &
    !$omp            shared( rsactive, nstate, maxstates, parai, n_streams_per_task, cntl, &
    !$omp                    jgw, geq0, llr1, lspin2, clsd, crge, c0, rswf, fpar, parm, spin_mod, &
    !$omp                    rhoe, psi, rhoes, psis, cp_cuda_env, cp_cufft, cp_curho, cp_cuwfn, n_nested_threads, &
    !$omp                    copy_data_to_host, copy_data_to_device, thread_views )

    !vw set number of children threads
    !$ CALL omp_set_num_threads ( n_nested_threads )

    i_stream = 0
    !$ i_stream = omp_get_thread_num()
    IF(cp_cuda_env%fft_n_devices_per_task>0) THEN
       thread_view = thread_views( i_stream ) !vw need to fix that
    ENDIF

    !vw create local comm for each stream/thread
    DO i = 0, n_streams_per_task - 1
       IF( i == i_stream ) THEN
          CALL mp_comm_dup ( parai%allgrp, fft_comm )
       ENDIF
       !$omp barrier
    ENDDO

    !vw assign pointer to data/buffer
    IF( cp_cuda_env%use_fft ) THEN
       rhoe_p => rhoes( :, :, i_stream )
       psi_p => psis( :, i_stream )
    ELSE
       rhoe_p => rhoe( :, : )
       psi_p => psi( : )
    ENDIF

    !vw zero memory on device if needed
    IF( .NOT. copy_data_to_host ) THEN
       CALL thread_view_get( thread_view, device_idx=device_idx, stream_idx=stream_idx )
       CALL cp_curho_stream_get_ptrs ( cp_curho, device_idx, stream_idx, rho_d=rho_d )       
       CALL cuda_mem_zero_bytes ( rho_d, rho_d%n_bytes )
    ENDIF

    IF( .NOT. copy_data_to_device ) THEN
       CALL thread_view_get( thread_view, device_idx=device_idx, stream_idx=stream_idx )
       CALL cp_cuwfn_device_get_ptrs ( cp_cuwfn, device_idx, c0_d=c0_d )
       CALL cp_cufft_device_get_ptrs ( cp_cufft, device_idx, nzfs_d=nzfs_d, inzs_d=inzs_d )
       CALL cp_cufft_stream_get_ptrs ( cp_cufft, device_idx, stream_idx, stream=stream_p, t1_d=psi_d )
    ENDIF

    !$omp do schedule( static, 1 )
    DO i = 1,part_1d_nbr_el_in_blk(nstate,parai%cp_inter_me,parai%cp_nogrp),2
       is1 = part_1d_get_el_in_blk(i,nstate,parai%cp_inter_me,parai%cp_nogrp)
       is2 = nstate+1
       IF (i+1.LE.part_1d_nbr_el_in_blk(nstate,parai%cp_inter_me,parai%cp_nogrp))&
            is2 = part_1d_get_el_in_blk(i+1,nstate,parai%cp_inter_me,parai%cp_nogrp)

       tfcal=rsactive
       tfcal=tfcal.OR.(crge%f(is1,1).NE.0._real_8)
       IF (is2.LE.nstate) tfcal=tfcal.OR.(crge%f(is2,1).NE.0._real_8)

       IF (tfcal) THEN
          __NVTX_TIMER_START ( procedureN//'_pre' )
          IF( .NOT. copy_data_to_device ) THEN
             CALL cuda_mem_zero_bytes ( psi_d, psi_d%n_bytes )
             IF (is2.GT.nstate) THEN
                CALL CuUser_setpsi_1_state_g  ( zone, cuda_z_points_to( c0_d, (is1-1)*SIZE(c0,1)+1 ), &
                     & psi_d, jgw, nzfs_d, inzs_d, geq0, stream_p )
             ELSE
                CALL CuUser_setpsi_2_states_g ( cuda_z_points_to( c0_d, (is1-1)*SIZE(c0,1)+1 ), &
                     & cuda_z_points_to( c0_d, (is2-1)*SIZE(c0,1)+1 ), psi_d, jgw, nzfs_d, inzs_d, geq0, stream_p )
             ENDIF
          ELSE
             CALL zeroing(psi_p)!,maxfftn)
             IF (is2.GT.nstate) THEN
                CALL set_psi_1_state_g(zone,c0(:,is1),psi_p)
             ELSE
                CALL set_psi_2_states_g(c0(:,is1),c0(:,is2),psi_p)
             ENDIF
          ENDIF
          __NVTX_TIMER_STOP

          ! ==--------------------------------------------------------------==
          ! ==  Fourier transform the wave functions to real space.         ==
          ! ==  In the array PSI was used also the fact that the wave       ==
          ! ==  functions at Gamma are real, to form a complex array (PSI)  ==
          ! ==  with the wave functions corresponding to two different      ==
          ! ==  states (i and i+1) as the real and imaginary part. This     ==
          ! ==  allows to call the FFT routine 1/2 of the times and save    ==
          ! ==  time.                                                       ==
          ! ==--------------------------------------------------------------==
          CALL invfftn(psi_p,.TRUE., fft_comm, thread_view=thread_view, copy_data_to_host=copy_data_to_host, &
               & copy_data_to_device=copy_data_to_device )

          ! Store real space wavefunctions
          IF (rsactive) THEN
             ! this is a quick fix (wont work for the cp_groups)
             ! we should have a mapping i -> state that we can get
             ! the correct state while reading rswf
             !
             ! SuperDirtyFix (SDF): rswf is allocated for all the
             ! states in the case we use more that 1 group.
             ! IF(ID.LE.MAXSTATES) THEN
             ! IF(CP_NOGRP.GT.1) CALL STOPGM(procedureN,'fix me')
             IF (is1.LE.maxstates) THEN
                iwf=is1/2+1
                CALL dcopy(2*fpar%nnr1,psi_p(1),1,rswf(1,iwf),1)
             ENDIF
          ENDIF

          __NVTX_TIMER_START ( procedureN//'_acc' )

          ! Compute the charge density from the wave functions
          ! in real space
          coef3=crge%f(is1,1)/parm%omega
          IF (is2.GT.nstate) THEN
             coef4=0.0_real_8
          ELSE
             coef4=crge%f(is2,1)/parm%omega
          ENDIF
          IF (cntl%tlsd) THEN
             ispin1=1
             ispin2=1
             IF (is1.GT.spin_mod%nsup) ispin1=2
             IF (is2.GT.spin_mod%nsup) ispin2=2
             IF (ispin1.EQ.ispin2) THEN
                CALL build_density_sum(coef3,coef4,psi_p,rhoe_p(:,ispin1),llr1)
             ELSE
                CALL build_density_real(coef3,psi_p,rhoe_p(:,ispin1),llr1)
                CALL build_density_imag(coef4,psi_p,rhoe_p(:,ispin2),llr1)
             ENDIF
          ELSEIF (lspin2%tlse) THEN
             CALL build_density_sum(coef3,coef4,psi_p,rhoe_p(:,1),llr1)
             IF (is1.EQ.clsd%ialpha) THEN
                CALL build_density_real(coef3,psi_p,rhoe_p(:,3),llr1)
             ELSEIF (is1.EQ.clsd%ibeta) THEN
                CALL build_density_real(coef3,psi_p,rhoe_p(:,2),llr1)
             ENDIF
             IF (is2.EQ.clsd%ialpha) THEN
                CALL build_density_imag(coef4,psi_p,rhoe_p(:,3),llr1)
             ELSEIF (is2.EQ.clsd%ibeta) THEN
                CALL build_density_imag(coef4,psi_p,rhoe_p(:,2),llr1)
             ENDIF
          ELSE
             IF( .NOT. copy_data_to_host ) THEN
                CALL cp_cubuild_density_sum ( coef3, coef4, llr1, cp_cufft, cp_curho, thread_view )
             ELSE
                CALL build_density_sum(coef3,coef4,psi_p,rhoe_p(:,1),llr1)
             ENDIF
          ENDIF
       ENDIF                 ! endif TFCAL

       __NVTX_TIMER_STOP

    ENDDO                     ! End loop over the electronic states
    !$omp end do


    !vw sync the streams
    IF( .NOT. copy_data_to_host ) THEN
       CALL thread_view_get( thread_view, device_idx=device_idx, stream_idx=stream_idx )
       CALL cp_cufft_stream_get_ptrs ( cp_cufft, device_idx, stream_idx, stream=stream_p )
       CALL cuda_stream_synchronize ( stream_p )
    ENDIF

    !vw free local comm for each stream/thread
    CALL mp_comm_free ( fft_comm )

    !$omp end parallel

    __NVTX_TIMER_START ( procedureN//'_c' )

    !vw merge the rhoes to rhoe
    IF( cp_cuda_env%use_fft ) THEN
       CALL tiset(procedureN//'_gpu_acc',isub4)
       IF( .NOT. copy_data_to_host ) THEN
          !vw special case accumulate rhoes to cpu
          CALL cp_cubuild_density_copy_to_host ( cp_curho, rhoes, thread_views )
       ENDIF
       DO i_stream = 0, n_streams_per_task - 1
          CALL daxpy ( SIZE(rhoe), 1.0_real_8, rhoes( 1, 1, i_stream ), 1, rhoe, 1 )
       ENDDO

       DEALLOCATE(rhoes,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       ! acm TODO - FIX THIS FOR RHOES
       CALL cuda_dealloc_host( psis )
       CALL cp_curho_dealloc_buffers ( cp_curho )

#ifndef _HASNT_OMP_SET_NESTED
       !$ CALL omp_set_nested( nested_orig )
#endif
       !$ CALL omp_set_max_active_levels( max_level_orig )
       CALL tihalt(procedureN//'_gpu_acc',isub4)
    ENDIF

    CALL thread_views_finalize ( thread_views )

    ! ==--------------------------------------------------------------==
    ! redistribute RHOE over the groups if needed
    !
    IF (parai%cp_nogrp.GT.1) THEN
       CALL tiset(procedureN//'_grps_b',isub3)
       CALL cp_grp_redist(rhoe,fpar%nnr1,clsd%nlsd)
       CALL tihalt(procedureN//'_grps_b',isub3)
    ENDIF

    ! CASPUR 2/9/04
    IF (tdgcomm%tdg) THEN
       IF (lspin2%tlse) THEN
          CALL stopgm(procedureN,&
               'TDG TOGETHER WITH LSE IS NOT YET IMPLEMENTED ',&
               __LINE__,__FILE__)
       ELSE IF (cntl%tlsd) THEN
          !$omp parallel do private(IR)
          DO ir=1,llr1
             psi(ir)=CMPLX(rhoe(ir,1), rhoe(ir,2),kind=real_8)
          ENDDO
          CALL fwfftn(psi,.FALSE.,parai%allgrp)
          CALL movepsih(psi,2)
          CALL setfftn(0)
          CALL invfftn(psi,.FALSE.,parai%allgrp)
          !$omp parallel do private(IR)
          DO ir=1,fpar%nnr1
             rhoe(ir,1)=REAL(psi(ir))
             rhoe(ir,2)=AIMAG(psi(ir))
          ENDDO
       ELSE
          !$omp parallel do private(IR)
          DO ir=1,llr1
             psi(ir)=CMPLX(rhoe(ir,1), 0._real_8,kind=real_8)
          ENDDO
          CALL fwfftn(psi,.FALSE.,parai%allgrp)
          CALL movepsih(psi,1)
          CALL setfftn(0)
          CALL invfftn(psi,.FALSE.,parai%allgrp)
          CALL dcopy(fpar%nnr1, psi, 2, rhoe, 1)
       ENDIF
    ENDIF
    ! MOVE DENSITY ACCORDING TO MOVEMENT OF ATOMS
    IF (ropt_mod%modens) CALL moverho(rhoe,psi)
    ! CONTRIBUTION OF THE VANDERBILT PP TO RHOE
    IF (pslo_com%tivan) THEN
       IF (cntl%tlsd) THEN
          ! ALPHA SPIN
          CALL rhov(nstate,1,spin_mod%nsup,rsumv,psi)
          rsum=rsum+parm%omega*rsumv
          !$omp parallel do private(I)
          DO i=1,fpar%nnr1
             rhoe(i,1)=rhoe(i,1)+REAL(psi(i))
          ENDDO
          ! BETA SPIN
          CALL rhov(nstate,spin_mod%nsup+1,nstate,rsumv,psi)
          rsum=rsum+parm%omega*rsumv
          !$omp parallel do private(I)
          DO i=1,fpar%nnr1
             rhoe(i,2)=rhoe(i,2)+REAL(psi(i))
          ENDDO
       ELSE
          CALL rhov(nstate,1,nstate,rsumv,psi)
          rsum=rsum+parm%omega*rsumv
          !$omp parallel do private(I)
          DO i=1,fpar%nnr1
             rhoe(i,1)=rhoe(i,1)+REAL(psi(i))
          ENDDO
       ENDIF
       ! Vanderbilt Charges
       IF (paral%parent) THEN
          ALLOCATE(qa(ions1%nat),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
          CALL zeroing(qa)!,ions1%nat)
          CALL augchg(fnl,crge%f,qa,nstate)
          iat=0
          DO is=1,ions1%nsp
             chrg%vdbchg(is)=0._real_8
             DO ia=1,ions0%na(is)
                iat=iat+1
                chrg%vdbchg(is)=chrg%vdbchg(is)+qa(iat)
             ENDDO
             chrg%vdbchg(is)=chrg%vdbchg(is)/REAL(ions0%na(is),kind=real_8)
          ENDDO
          DEALLOCATE(qa,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
               __LINE__,__FILE__)
       ENDIF
    ENDIF

    ! ALPHA+BETA DENSITY IN RHOE(*,1), BETA DENSITY IN RHOE(*,2)
    chrg%csums=0._real_8
    chrg%csumsabs=0._real_8
    IF (cntl%tlsd) THEN
       rsum1=0._real_8
       rsum1abs=0._real_8
       !$omp parallel do private(I) shared(fpar,RHOE) &
       !$omp  reduction(+:RSUM1,RSUM1ABS)
       DO i=1,fpar%nnr1
          rsum1 = rsum1 + (rhoe(i,1) - rhoe(i,2))
          rsum1abs = rsum1abs + ABS(rhoe(i,1) - rhoe(i,2))
          rhoe(i,1) = rhoe(i,1) + rhoe(i,2)
       ENDDO
       chrg%csums=rsum1*parm%omega/REAL(lr1s*lr2s*lr3s,kind=real_8)
       chrg%csumsabs=rsum1abs*parm%omega/REAL(lr1s*lr2s*lr3s,kind=real_8)
       ! ALPHA+BETA DENSITY IN RHOE(*,1), BETA DENSITY IN RHOE(*,2) ; M STATE
       ! ALPHA+BETA DENSITY IN RHOE(*,3), BETA DENSITY IN RHOE(*,4) ; T STATE
    ELSEIF (lspin2%tlse) THEN
       IF (lspin2%tcas22) THEN
          ! Calculate "ground state" density and D-EX DENSITY
          !$omp parallel do private(I,RTO,RBE,RAL)
          DO i=1,fpar%nnr1
             rto=rhoe(i,1)
             rbe=rhoe(i,2)
             ral=rhoe(i,3)
             rhoe(i,6)=rto-rbe+ral
             rhoe(i,7)=rto+rbe-ral
          ENDDO
       ENDIF
       IF (lspin2%tlsets) THEN
          !$omp parallel do private(I,RTO,RBE,RAL,RSP)
          DO i=1,fpar%nnr1
             rto=rhoe(i,1)
             rbe=rhoe(i,2)
             ral=rhoe(i,3)
             rsp=0.5_real_8*(rto-ral-rbe)
             rhoe(i,2)=rsp+rbe+o3*ral
             rhoe(i,3)=rto
             rhoe(i,4)=rsp+rbe+2._real_8*o3*ral
          ENDDO
       ELSE
          !$omp parallel do private(I,RTO,RBE,RAL,RSP)
          DO i=1,fpar%nnr1
             rto=rhoe(i,1)
             rbe=rhoe(i,2)
             ral=rhoe(i,3)
             rsp=0.5_real_8*(rto-ral-rbe)
             rhoe(i,2)=rsp+rbe
             rhoe(i,3)=rto
             rhoe(i,4)=rsp+ral+rbe
          ENDDO
       ENDIF
       IF (lspin2%tross.OR.lspin2%tcas22.OR.lspin2%tpenal) THEN
          ! WE ALSO NEED THE A*B DENSITY FOR THE EXCHANGE CONTRIBUTION
          ! OR THE OFF DIAGONAL ELEMENTS IN THE CAS22 METHOD
          CALL zeroing(rhoe(:,5))!,nnr1)
          CALL rhoabofr(1,c0(:,clsd%ialpha:clsd%ialpha),c0(:,clsd%ibeta:clsd%ibeta),rhoe(:,5),psi)
       ENDIF
    ENDIF

    ! HERE TO CHECK THE INTEGRAL OF THE CHARGE DENSITY
    ! RSUM1=DASUM(NNR1,RHOE(1,1),1)
    ! --> with VDB PP RHOE might be negative in some points
    rsum1=0._real_8
#if defined(__SR8000)
    !poption parallel, tlocal(I), psum(RSUM1)
#else
    !$omp parallel do private(I) shared(fpar,RHOE) &
    !$omp  reduction(+:RSUM1)
#endif
    DO i=1,fpar%nnr1
       rsum1=rsum1+rhoe(i,1)
    ENDDO
    rsum1=rsum1*parm%omega/REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)
    chrg%csumg=rsum
    chrg%csumr=rsum1
    CALL mp_sum(chrg%csumg,parai%allgrp)
    CALL mp_sum(chrg%csumr,parai%allgrp)
    CALL mp_sum(chrg%csums,parai%allgrp)
    CALL mp_sum(chrg%csumsabs,parai%allgrp)

    IF (paral%parent.AND.ABS(chrg%csumr-chrg%csumg).GT.delta) THEN
       IF (paral%io_parent)&
            WRITE(6,'(A,T46,F20.12)') ' IN FOURIER SPACE:', chrg%csumg
       IF (paral%io_parent)&
            WRITE(6,'(A,T46,F20.12)') ' IN REAL SPACE:', chrg%csumr
       IF ((symmi%indpg.NE.0.AND.dual00%cdual.LT.4._real_8).AND.paral%io_parent)&
            WRITE(6,*) 'YOUR DUAL NUMBER ',dual00%cdual,&
            ' COULD BE TOO SMALL WITH DENSITY SYMMETRISATION'
       CALL stopgm(procedureN,'TOTAL DENSITY SUMS ARE NOT EQUAL',&
            __LINE__,__FILE__)
    ENDIF
    ! TAU FUNCTION
    IF (cntl%ttau) CALL tauofr(c0,psi,nstate)
    !
    __NVTX_TIMER_STOP
    __NVTX_TIMER_STOP
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rhoofr

  ! ==================================================================
  SUBROUTINE give_scr_rhoofr(lrhoofr,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lrhoofr
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: ldg, lmoverho, lrhov, lsymrho

    ldg=0
    IF (tdgcomm%tdg) ldg=2*ncpw%nhg*clsd%nlsd+100
    lmoverho=0
    IF (pslo_com%tivan) THEN ! RHOOFR and RHOOFR_C
       CALL give_scr_rhov(lrhov,tag)
       lrhoofr=MAX(ions1%nat,lrhov)! AUGCHG
    ELSE
       lrhoofr=0
    ENDIF
    CALL give_scr_symrho(lsymrho,tag)
    IF (ropt_mod%modens) CALL give_scr_moverho(lmoverho,tag)
    lrhoofr=MAX(lrhoofr,lsymrho,lmoverho)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_rhoofr
  ! ==================================================================
  SUBROUTINE movepsih(a,itask)
    COMPLEX(real_8)                          :: a(:)
    INTEGER                                  :: itask

    CHARACTER(*), PARAMETER                  :: procedureN = 'movepsih'

    COMPLEX(real_8)                          :: aq1, aq2, c1, c2
    COMPLEX(real_8), ALLOCATABLE             :: b1(:), b2(:)
    INTEGER                                  :: ierr, ig
    LOGICAL                                  :: debug

! ==--------------------------------------------------------------==

    debug=.FALSE.
    !
    IF (itask.EQ.1) THEN
       ALLOCATE(b1(2*ncpw%nhg),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       !$omp parallel do private(IG)
       DO ig=1,jhg
          b1(ig)=a(nzff(ig))
       ENDDO
       CALL zeroing(a)!,maxfft)
       DO ig=1,jhg
          a(nzh(ig))=b1(ig)
          a(indz(ig))=CONJG(b1(ig))
       ENDDO
       IF (geq0) a(nzh(1))=b1(1)

       DEALLOCATE(b1,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    ELSEIF (itask.EQ.2) THEN
       ALLOCATE(b1(2*ncpw%nhg),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(b2(2*ncpw%nhg),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)

       c1=CMPLX(0.5_real_8, 0._real_8,kind=real_8)
       c2=CMPLX(0.0_real_8,-0.5_real_8,kind=real_8)
       ! THE TWO COMPONENT OF DENSITY (REAL RUNCTION IN R SPACE)
       ! TRANSFORMED AT ONCE
       !$omp parallel do private(IG,AQ1,AQ2) shared(C1,C2)
       DO ig=1,jhg
          ! mb           AQ1=C1*(A(NZFF(IG))+CONJG(A(INZF(IG))))
          ! mb           AQ2=C2*(A(NZFF(IG))-CONJG(A(INZF(IG))))
          aq1=a(nzff(ig))    ! access only once
          aq2=CONJG(a(inzf(ig)))! memory locations
          b1(ig)=c1*(aq1+aq2)
          b2(ig)=c2*(aq1-aq2)
       ENDDO
       IF (geq0) THEN
          b1(1) = CMPLX(REAL(a(nzff(1))),0._real_8,kind=real_8)
          b2(1) = CMPLX(AIMAG(a(nzff(1))),0._real_8,kind=real_8)
       ENDIF
       CALL zeroing(a)!,maxfft)
       DO ig=1,jhg
          a(nzh(ig))=b1(ig)+uimag*b2(ig)
          a(indz(ig))=CONJG(b1(ig))+uimag*CONJG(b2(ig))
       ENDDO
       IF (geq0) a(nzh(1))=b1(1)+uimag*b2(1)

       DEALLOCATE(b1,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(b2,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    ELSE
       CALL stopgm('MOVEPSIH',' ITASK=? ',&
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE movepsih
  ! ==================================================================

END MODULE rhoofr_utils
