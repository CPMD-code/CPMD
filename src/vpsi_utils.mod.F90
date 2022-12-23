#include "cpmd_global.h"

MODULE vpsi_utils
  USE cdftmod,                         ONLY: cdftcom,&
                                             cdftlog,&
                                             wdiff
  USE cnst,                            ONLY: uimag
  USE cp_cuda_types,                   ONLY: cp_cuda_devices_fft,&
                                             cp_cuda_env
  USE cp_cufft_types,                  ONLY: cp_cufft,&
                                             cp_cufft_device_get_ptrs,&
                                             cp_cufft_stream_get_ptrs
  USE cp_cuvpsi_types,                 ONLY: cp_cuvpsi_t
  USE cp_cuvpsi_utils,                 ONLY: cp_cuapply_potential,&
                                             cp_cuvpotdg_copy_to_device,&
                                             cp_cuvpsi_finalize_and_dealloc,&
                                             cp_cuvpsi_init_and_alloc
  USE cp_cuwfn_types,                  ONLY: cp_cuwfn,&
                                             cp_cuwfn_device_get_ptrs,&
                                             cp_cuwfn_get
  USE cp_grp_utils,                    ONLY: cp_grp_redist
  USE cppt,                            ONLY: gk,&
                                             hg,&
                                             indzs,&
                                             nzh,&
                                             nzhs
  USE cuda_types,                      ONLY: cuda_memory_t,&
                                             cuda_stream_t
  USE cuda_utils,                      ONLY: cuda_alloc_host,&
                                             cuda_dealloc_host,&
                                             cuda_mem_zero_bytes,&
                                             cuda_z_points_to
  USE cuuser_utils,                    ONLY: CuUser_setpsi_1_state_g,&
                                             CuUser_setpsi_2_states_g
  USE dg,                              ONLY: ipooldg,&
                                             tdgcomm
  USE efld,                            ONLY: extf
  USE error_handling,                  ONLY: stopgm
  USE fft,                             ONLY: inzf,&
                                             jgw,&
                                             jhg,&
                                             llr1,&
                                             ngrm,&
                                             nzff
  USE fftmain_utils,                   ONLY: fwfftn,&
                                             invfftn
  USE fftnew_utils,                    ONLY: setfftn
  USE geq0mod,                         ONLY: geq0
  USE kinds,                           ONLY: real_8
  USE kpclean_utils,                   ONLY: c_clean
  USE kpnt,                            ONLY: hgkm,&
                                             hgkp,&
                                             rk
  USE kpts,                            ONLY: tkpts
  USE mp_interface,                    ONLY: mp_comm_dup,&
                                             mp_comm_free,&
                                             mp_comm_null
  USE parac,                           ONLY: parai
  USE part_1d,                         ONLY: part_1d_get_el_in_blk,&
                                             part_1d_nbr_el_in_blk
  USE prcp,                            ONLY: prcp_com
  USE reshaper,                        ONLY: reshape_inplace
  USE rswfmod,                         ONLY: maxstates,&
                                             rsactive,&
                                             rswf
  USE special_functions,               ONLY: cp_erf
  USE spin,                            ONLY: clsd,&
                                             lspin2,&
                                             spin_mod
  USE state_utils,                     ONLY: add_wfn,&
                                             set_psi_1_state_g,&
                                             set_psi_1_state_g_kpts,&
                                             set_psi_2_states_g
  USE switch_functionals_utils,        ONLY: switch_functionals
  USE system,                          ONLY: &
       cntl, fpar, group, locpot2, ncpw, nkpt, parap, parm, spar
  USE td_input,                        ONLY: gampl,&
                                             gndir,&
                                             gpot,&
                                             gpotv,&
                                             td_prop
  USE thread_view_types,               ONLY: thread_view_get,&
                                             thread_view_t
  USE thread_view_utils,               ONLY: thread_views_finalize,&
                                             thread_views_init
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE vpsi_lse_utils,                  ONLY: vpsi_lse
  USE vtaupsi_utils,                   ONLY: vtaupsi
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

  PUBLIC :: vpsi
  PUBLIC :: vpsimt
  !public :: movepsid

CONTAINS

  ! if (.not.tdg) r*8 vpotdg(nnr1,*) => r*8 vpot(nnr1,*)
  ! else vpotdg = copy of psi(:: 2)
  !
  ! r*8 vpotx(*) => r*8 vpotdg(1,1)
  ! r*8 psix(*) => c*8 psi(*) for use in PARALLEL sections
  ! ==================================================================
  SUBROUTINE vpsi(c0,c2,f,vpot,psi,nstate,ikind,ispin,redist_c2)
    ! ==================================================================
    ! == K-POINT AND NOT K-POINT VERSION OF VPSI.                     ==
    ! ==--------------------------------------------------------------==
    ! == VPOT:   IN INPUT POTENTIAL                                   ==
    ! == ISPIN:  Need with LSD option for diagonalization scheme      ==
    ! ==         dimension of VPOT(NNR1,ISPIN)                        ==
    ! ==--------------------------------------------------------------==
    ! EHR[
    ! EHR]
    COMPLEX(real_8)                          :: c0(:,:), c2(:,:)
    REAL(real_8)                             :: f(:)
    REAL(real_8), TARGET __CONTIGUOUS        :: vpot(:,:)
    COMPLEX(real_8), TARGET __CONTIGUOUS     :: psi(:)
    INTEGER                                  :: nstate, ikind, ispin
    LOGICAL                                  :: redist_c2

    CHARACTER(*), PARAMETER                  :: procedureN = 'vpsi'
    COMPLEX(real_8), PARAMETER               :: zone = (1.0_real_8,0.0_real_8)

    COMPLEX(real_8)                          :: fm, fp, psii, psin
    COMPLEX(real_8), ALLOCATABLE             :: C2_vpsi(:,:)
    COMPLEX(real_8), DIMENSION(:), &
      POINTER __CONTIGUOUS                   :: psi_p
    COMPLEX(real_8), DIMENSION(:, :), &
      POINTER __CONTIGUOUS                   :: psis
    INTEGER :: device_idx, fft_comm, i, i_stream, iclpot = 0, id, ierr, ig, &
      ir, is1, is2, isub, isub2, isub3, isub4, iwf, ixx, ixxs, iyy, izz, jj, &
      leadx, lspin, n_max_threads, n_nested_threads, n_streams_per_task, &
      njump, nnrx, nostat, nrxyz1s, nrxyz2, stream_idx
    INTEGER, SAVE                            :: dgfirst = 0
    LOGICAL :: copy_data_to_device, copy_data_to_device_inv, &
      copy_data_to_host, lg_vpotx3a, lg_vpotx3b, ready, wfcopy
    REAL(real_8)                             :: chksum, csmult, fi, fip1, g2, &
                                                xskin
    REAL(real_8), ALLOCATABLE                :: vpotx3a(:,:,:), vpotx3b(:,:,:)
    REAL(real_8), POINTER __CONTIGUOUS       :: VPOTX(:)
    REAL(real_8), POINTER, SAVE __CONTIGUOUS :: vpotdg(:,:)
    TYPE(cp_cuvpsi_t)                        :: cp_cuvpsi
    TYPE(cuda_memory_t), POINTER             :: c0_d, inzs_d, nzfs_d, psi_d
    TYPE(cuda_stream_t), POINTER             :: stream_p
    TYPE(thread_view_t)                      :: thread_view
    TYPE(thread_view_t), ALLOCATABLE, &
      DIMENSION(:)                           :: thread_views

    !$ LOGICAL :: nested_orig
    !$ INTEGER :: max_level_orig
    ! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    __NVTX_TIMER_START ( procedureN )
    __NVTX_TIMER_START ( procedureN//'_a' )

    IF (group%nogrp.GT.1)CALL stopgm(procedureN,&
         'OLD TASK GROUPS NOT SUPPORTED ANYMORE ',&
         __LINE__,__FILE__)
    IF (tdgcomm%tdg) THEN
       IF (lspin2%tlse)CALL stopgm(procedureN,&
            'TDG TOGETHER WITH  LSE IS NOT YET IMPLEMENTED ',&
            __LINE__,__FILE__)
       CALL setfftn(ipooldg)
       IF (dgfirst.EQ.0) THEN
          ALLOCATE(vpotdg(fpar%nnr1,2),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
          dgfirst=1
       ENDIF
    ELSE
       CALL setfftn(0)
    ENDIF
    IF (lspin2%tlse.AND.lspin2%tlsets.AND.(lspin2%tross.OR.lspin2%tcas22.OR.lspin2%tpenal.OR.lspin2%troot)) THEN
       CALL stopgm(procedureN,&
            'NO SLATER TS WITH ROSS, CAS22, PENALTY, ROOTHAAN',&
            __LINE__,__FILE__)
    ENDIF

    !
    ! Need to allocate a buffer to accumulate the vpsi part of C2
    !
    ALLOCATE(C2_vpsi(nkpt%ngwk,nstate),stat=ierr)
    IF (ierr.NE.0) CALL stopgm(procedureN,'Allocation problem',&
         __LINE__,__FILE__)

    CALL zeroing(C2_vpsi)!, nkpt%ngwk*nstate)

    ! For the moment we synchronize the C0 and C2 if the cp groups are used
    ! That should be done else where
    ! IF(cp_nogrp.GT.1) THEN
    ! CALL MY_BCAST(C0,2*NGW*NSTATE*8/IRAT,0,CP_INTER_GRP)
    ! CALL MY_BCAST(C2,2*NGW*NSTATE*8/IRAT,0,CP_INTER_GRP)
    ! ENDIF

    leadx = fpar%nnr1
    nnrx  = llr1
    !
    IF (tdgcomm%tdg) THEN
       CALL zeroing(psi)!, maxfft)
       CALL dcopy (fpar%nnr1, vpot, 1, psi, 2)
       CALL fwfftn(psi,.FALSE.,parai%allgrp)
       CALL movepsid(psi)
       CALL invfftn(psi, .FALSE.,parai%allgrp)
       CALL dcopy(llr1, psi, 2, vpotdg, 1)
       IF (cntl%tlsd.AND.ispin.EQ.2) THEN
          CALL zeroing(psi)!, maxfft)
          CALL dcopy (fpar%nnr1, vpot(1,2), 1, psi, 2)
          CALL fwfftn(psi,.FALSE.,parai%allgrp)
          CALL movepsid(psi)
          CALL invfftn(psi, .FALSE.,parai%allgrp)
          CALL dcopy(llr1, psi, 2, vpotdg(1,2), 1)
       ENDIF
    ELSE
       vpotdg => vpot
    ENDIF
    !
    !vpotx(1:SIZE(vpotdg)) => vpotdg
    CALL reshape_inplace(vpotdg, (/fpar%nnr1*ispin/), vpotx)

    njump=2
    IF (tkpts%tkpnt) njump=1
    IF (lspin2%tlse) THEN
       nostat = clsd%ialpha-1
    ELSE
       nostat = nstate
    ENDIF

    IF (cntl%cdft)THEN
       IF (.NOT.cdftlog%tcall)THEN
          csmult=cdftcom%cdft_v(1)
       ELSE
          csmult=cdftcom%cdft_v(1)+cdftcom%cdft_v(2)
       ENDIF
       !$omp parallel do private(IR)
       DO ir=1,fpar%nnr1
          vpotdg(ir,1)=vpotdg(ir,1)+csmult*wdiff(ir)
       ENDDO
       IF (.NOT.cdftlog%tcall)THEN
          IF (cdftlog%tspinc) THEN
             csmult=-cdftcom%cdft_v(1)
          ELSE
             csmult=cdftcom%cdft_v(1)
          ENDIF
       ELSE
          csmult=cdftcom%cdft_v(1)-cdftcom%cdft_v(2)
       ENDIF
       IF (cntl%tlsd.AND.ispin.EQ.2) THEN
          !$omp parallel do private(IR)
          DO ir=1,fpar%nnr1
             vpotdg(ir,2)=vpotdg(ir,2)+csmult*wdiff(ir)
          ENDDO
       ENDIF
    ENDIF

    ! ==--------------------------------------------------------------==

    !acm gpu version:
    !acm if possible do the accumulation of the potential on the gpu
    copy_data_to_host   = (cntl%tlsd.AND.ispin.EQ.2) .OR. ( .NOT. cp_cuda_env%use_fft )
    copy_data_to_device = copy_data_to_host
    CALL cp_cuwfn_get ( cp_cuwfn, ready=ready, chksum=chksum )
    copy_data_to_device_inv = .NOT. ( ready .AND. cp_cuda_env%use_fft .AND. &
         &                           .NOT. tkpts%tkpnt .AND. &
         &                            ABS( SUM( ABS( c0(:,1) ) ) - chksum ) < 1.0e-12_real_8 )
    !write(*,*) 'vpsi_utils.mod.F90: copy_data_to_device_inv=',copy_data_to_device_inv
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
       CALL cuda_alloc_host(psis,[SIZE(psi),n_streams_per_task], lbs=[1,0])

       IF( .NOT. copy_data_to_host ) THEN
          ! TODO acm: the init+allocation/deallocation should be done outside (once for entire code)
          CALL cp_cuvpsi_init_and_alloc( cp_cuda_env, cp_cuda_devices_fft, cp_cuvpsi, SIZE(vpotdg) )
          CALL cp_cuvpotdg_copy_to_device( cp_cuda_env, cp_cuvpsi, vpotdg )
       ENDIF

       IF( .NOT. copy_data_to_device_inv ) THEN
          DO device_idx = 1, cp_cuda_env%fft_n_devices_per_task
             !CALL cp_cuwfn_device_get_ptrs ( cp_cuwfn, device_idx, c0_d=c0_d )
             !CALL cuda_memcpy_host_to_device ( c0, c0_d )
          ENDDO
       ENDIF
       CALL tihalt(procedureN//'_gpu_init',isub2)
    ENDIF

    ! acm copy vpotdg on device
    __NVTX_TIMER_STOP

    !$omp parallel   if( cp_cuda_env%use_fft ) &
    !$omp            num_threads( n_streams_per_task ) &
    !$omp            default( none ) &
    !$omp            private( i, is1, is2, iwf, wfcopy, i_stream, psi_p, &
    !$omp                     xskin, fip1, fm, fp, fi, ixxs, izz, iyy, jj, ixx, &
    !$omp                     lg_vpotx3b, lg_vpotx3a, nrxyz2, nrxyz1s, iclpot, &
    !$omp                     ierr, lspin, vpotx3a, vpotx3b, g2, psii, psin, thread_view, &
    !$omp                     stream_idx, device_idx, stream_p, psi_d, c0_d, nzfs_d, &
    !$omp                     inzs_d, ig, id, ir ) &
    !$omp            firstprivate( fft_comm ) &
    !$omp            shared( rsactive, hg, jgw, prcp_com, geq0, gpotv, gampl, &
    !$omp                    gpot, gk, rk, hgkm, hgkp, indzs, nzhs, ncpw, gndir, &
    !$omp                    f, parap, parm, spar, locpot2, vpotdg, td_prop, spin_mod, &
    !$omp                    cntl, ispin, psi, psis, parai, nostat, njump, n_streams_per_task, &
    !$omp                    ikind, nstate, maxstates, tkpts, nnrx, rswf, c0, vpotx, &
    !$omp                    leadx, extf, c2_vpsi, cp_cufft, cp_cuvpsi, cp_cuda_env, n_nested_threads, &
    !$omp                    copy_data_to_host, copy_data_to_device, copy_data_to_device_inv, &
    !$omp                    thread_views, cp_cuwfn )

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

    IF( cp_cuda_env%use_fft ) THEN
       psi_p => psis( :, i_stream )
    ELSE
       psi_p => psi( : )
    ENDIF

    IF( .NOT. copy_data_to_device_inv ) THEN
       CALL thread_view_get( thread_view, device_idx=device_idx, stream_idx=stream_idx )
       CALL cp_cuwfn_device_get_ptrs ( cp_cuwfn, device_idx, c0_d=c0_d )
       CALL cp_cufft_device_get_ptrs ( cp_cufft, device_idx, nzfs_d=nzfs_d, inzs_d=inzs_d )
       CALL cp_cufft_stream_get_ptrs ( cp_cufft, device_idx, stream_idx, stream=stream_p, t1_d=psi_d )
    ENDIF

    !$omp do schedule( static, 1 )
    DO i = 1,part_1d_nbr_el_in_blk(nostat,parai%cp_inter_me,parai%cp_nogrp),njump
       is1 = part_1d_get_el_in_blk(i,nostat,parai%cp_inter_me,parai%cp_nogrp)
       is2 = nostat+1
       IF (njump.EQ.2) THEN
          IF (i+1.LE.part_1d_nbr_el_in_blk(nostat,parai%cp_inter_me,parai%cp_nogrp))&
               is2 = part_1d_get_el_in_blk(i+1,nostat,parai%cp_inter_me,parai%cp_nogrp)
       ENDIF

       wfcopy=.FALSE.
       IF (rsactive) THEN
          iwf = is1+(ikind-1)*nstate
          IF (iwf.LE.maxstates) THEN
             wfcopy=.TRUE.
             copy_data_to_host=.TRUE.
             copy_data_to_device=.TRUE.
             copy_data_to_device_inv=.TRUE.
          ENDIF
       ENDIF
       IF (wfcopy) THEN
          iwf=is1/2+1
          IF (tkpts%tkpnt) iwf = is1+(ikind-1)*nstate
          CALL dcopy(2*nnrx,rswf(1,iwf),1,psi_p(1),1)
       ELSE
          __NVTX_TIMER_START ( procedureN//'_pre' )
          ! ==----------------------------------------------------------==
          ! == Store the wavefunctions in array PSI                     ==
          ! == in the following way:                                    ==
          ! == Since we use the Gamma point only, the wavefunctions in  ==
          ! == real space are real. This fact is used to save time by   ==
          ! == using complex arrays where the wavefunctions of          ==
          ! == 2 different states (i and i+1) are combined, and         ==
          ! == performing complex Fourier Transforms.                   ==
          ! == Warning: If you use K-points other than Gamma, or zone   ==
          ! == boundary, this symmetry is broken and the trick is not   ==
          ! == valid anymore. First step is to build the array PSI.     ==
          ! ==----------------------------------------------------------==
          IF( .NOT. copy_data_to_device_inv ) THEN
             CALL cuda_mem_zero_bytes ( psi_d, psi_d%n_bytes )
             IF (tkpts%tkpnt) THEN
                CALL stopgm(procedureN,'k-points with wfn on GPU NYI',&
                     __LINE__,__FILE__)
             ELSE
                IF (is2.GT.nostat) THEN
                   CALL CuUser_setpsi_1_state_g ( zone, cuda_z_points_to( c0_d, (is1-1)*SIZE(c0,1)+1 ), &
                        & psi_d, jgw, nzfs_d, inzs_d, geq0, stream_p )
                ELSE
                   CALL CuUser_setpsi_2_states_g ( cuda_z_points_to( c0_d, (is1-1)*SIZE(c0,1)+1 ), &
                        & cuda_z_points_to( c0_d, (is2-1)*SIZE(c0,1)+1 ), psi_d, jgw, nzfs_d, inzs_d, geq0, stream_p )
                ENDIF
             ENDIF
          ELSE
             CALL zeroing(psi_p)!, maxfft)
             IF (tkpts%tkpnt) THEN
                CALL set_psi_1_state_g_kpts(zone,c0(:,is1),psi_p)
             ELSE
                IF (is2.GT.nostat) THEN
                   CALL set_psi_1_state_g(zone,c0(:,is1),psi_p)
                ELSE
                   CALL set_psi_2_states_g(c0(:,is1),c0(:,is2),psi_p)
                ENDIF
             ENDIF
          ENDIF
          __NVTX_TIMER_STOP
          ! ==----------------------------------------------------------==
          ! == Transform the wavefunction to real space                 ==
          ! ==----------------------------------------------------------==
          CALL invfftn(psi_p,.TRUE., fft_comm, thread_view=thread_view, copy_data_to_host=copy_data_to_host, &
               & copy_data_to_device=copy_data_to_device_inv )
       ENDIF
       ! ==------------------------------------------------------------==
       ! == Apply the potential (V), which acts in real space.         ==
       ! ==------------------------------------------------------------==
       __NVTX_TIMER_START ( procedureN//'_apply' )
       IF (cntl%tlsd.AND.ispin.EQ.2) THEN
          IF (is1.EQ.spin_mod%nsup) THEN
#ifdef __SR8000
             !poption parallel
#endif
             ! EHR[
             IF (td_prop%td_extpot) THEN
                !$omp parallel do private(IR)
                DO ir=1,nnrx
                   psi_p(ir)=(vpotx(ir)+extf(ir))*REAL(psi_p(ir))&
                        +uimag*(vpotx(leadx+ir)+extf(ir))*AIMAG(psi_p(ir))
                ENDDO
             ELSE
                !$omp parallel do private(IR)
                DO ir=1,nnrx
                   psi_p(ir)=vpotx(ir)*REAL(psi_p(ir))&
                        +uimag*vpotx(leadx+ir)*AIMAG(psi_p(ir))
                ENDDO
             ENDIF
             ! EHR]
          ELSE
             !lspin=0
             !IF (is1.GT.spin_mod%nsup) lspin=leadx
             lspin=1
             IF (is1.GT.spin_mod%nsup) lspin=2
#ifdef __SR8000
             !poption parallel
#endif
             !$omp parallel do private(IR)
             DO ir=1,nnrx
                psi_p(ir)=vpotdg(ir,lspin)*psi_p(ir)
             ENDDO
          ENDIF
       ELSE
          IF( .NOT. copy_data_to_host ) THEN
             CALL cp_cuapply_potential ( nnrx, cp_cufft, cp_cuvpsi, thread_view )
          ELSE
             !$omp parallel do private(IR) shared(NNRX,PSI_p)
#ifdef __SR8000
             !poption parallel
#endif
             DO ir=1,nnrx
                psi_p(ir)=vpotdg(ir,1)*psi_p(ir)
             ENDDO
          ENDIF
          ! kk-mb === print local potential (start) ===
          IF (locpot2%tlpot) THEN
             iclpot=iclpot+1
             nrxyz1s=spar%nr1s*spar%nr2s*spar%nr3s
             nrxyz2=0
             IF (iclpot .EQ. 1) THEN
                ALLOCATE(vpotx3a(spar%nr1s,spar%nr2s,spar%nr3s),STAT=ierr)
                IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                     __LINE__,__FILE__)
                CALL zeroing(vpotx3a)!,nrxyz1s)
                lg_vpotx3a=.TRUE.
                lg_vpotx3b=.FALSE.
             ELSE
                IF (lg_vpotx3a) THEN
                   ALLOCATE(vpotx3b(spar%nr1s,spar%nr2s,spar%nr3s),STAT=ierr)
                   IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                        __LINE__,__FILE__)
                   CALL zeroing(vpotx3b)!,nrxyz1s)
                   DEALLOCATE(vpotx3a,STAT=ierr)
                   IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
                        __LINE__,__FILE__)
                   lg_vpotx3a=.FALSE.
                   lg_vpotx3b=.TRUE.
                ELSE IF (lg_vpotx3b) THEN
                   ALLOCATE(vpotx3a(spar%nr1s,spar%nr2s,spar%nr3s),STAT=ierr)
                   IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                        __LINE__,__FILE__)
                   CALL zeroing(vpotx3a)!,nrxyz1s)
                   DEALLOCATE(vpotx3b,STAT=ierr)
                   IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
                        __LINE__,__FILE__)
                   lg_vpotx3a=.TRUE.
                   lg_vpotx3b=.FALSE.
                ENDIF
             ENDIF
             DO ir=1,nnrx
                IF (vpotx(ir) .NE. 0) THEN
                   nrxyz2=nrxyz2+1
                   ixx=MOD(nrxyz2-1,parm%nr1)+1
                   jj=(nrxyz2-1)/parm%nr1+1
                   iyy=MOD(jj-1,spar%nr2s)+1
                   izz=(nrxyz2-1)/(parm%nr1*spar%nr2s)+1
                   ixxs=ixx+parap%nrxpl(parai%mepos,1)-1
                   IF (lg_vpotx3a) THEN
                      vpotx3a(ixxs,iyy,izz)=vpotx(ir)
                   ELSE IF (lg_vpotx3b) THEN
                      vpotx3b(ixxs,iyy,izz)=vpotx(ir)
                   ENDIF
                ENDIF
             ENDDO
          ENDIF
          ! kk-mb === print local potential (end)  ===
       ENDIF
       __NVTX_TIMER_STOP
       ! ==------------------------------------------------------------==
       ! == Back transform to reciprocal space the product V.PSI       ==
       ! ==------------------------------------------------------------==
       CALL fwfftn(psi_p,.TRUE., fft_comm, thread_view=thread_view, copy_data_to_device=copy_data_to_device )
       ! ==------------------------------------------------------------==
       ! == Decode the combination of wavefunctions (now multiplied by ==
       ! == the potential), back to those of states i and i+1. Multiply==
       ! == also by the occupation numbers f(i) and other constants in ==
       ! == order to obtain the force on the electronic degrees of     ==
       ! == freedom, stored in array FC.                               ==
       ! == Here we also add the force from the kinetic energy.        ==
       ! ==------------------------------------------------------------==
       __NVTX_TIMER_START ( procedureN//'_acc' )
       IF (tkpts%tkpnt) THEN
          fi=f(is1)
          IF (fi.EQ.0._real_8) fi=2._real_8
          ! EHR[
          IF (cntl%tgaugep.OR.cntl%tgaugef) THEN
             IF (gndir.EQ.0) THEN
                DO ig=1,ncpw%ngw
                   fp=psi_p(nzhs(ig))
                   fm=psi_p(indzs(ig))
                   C2_vpsi(ig,is1)=-fi*&
                        (0.5_real_8*parm%tpiba2*hgkp(ig,ikind)*c0(ig,is1)+fp)
                   C2_vpsi(ig+ncpw%ngw,is1)=-fi*&
                        (0.5_real_8*parm%tpiba2*hgkm(ig,ikind)*c0(ig+ncpw%ngw,is1)+fm)
                   C2_vpsi(ig,is1)=c2_vpsi(ig,is1)&
                        -fi*gpot*(rk(td_prop%pertdir,1)&
                        +gk(td_prop%pertdir,ig))*parm%tpiba*c0(ig ,is1)&
                        -fi*0.5_real_8*(gpot**2)*              c0(ig,is1)
                   C2_vpsi(ig+ncpw%ngw,is1)=c2_vpsi(ig+ncpw%ngw,is1)&
                        -fi*gpot*(rk(td_prop%pertdir,1)&
                        -gk(td_prop%pertdir,ig))*parm%tpiba*c0(ig+ncpw%ngw,is1)&
                        -fi*0.5_real_8*(gpot**2)*              c0(ig+ncpw%ngw,is1)
                ENDDO
             ELSE
                DO ig=1,ncpw%ngw
                   fp=psi_p(nzhs(ig))
                   fm=psi_p(indzs(ig))
                   C2_vpsi(ig,is1)=-fi*&
                        (0.5_real_8*parm%tpiba2*hgkp(ig,ikind)*c0(ig,is1)+fp)
                   C2_vpsi(ig+ncpw%ngw,is1)=-fi*&
                        (0.5_real_8*parm%tpiba2*hgkm(ig,ikind)*c0(ig+ncpw%ngw,is1)+fm)
                   ! elisa
                   DO id=1,3
                      IF (gampl(id).NE.0._real_8) THEN
                         C2_vpsi(ig,is1)=c2_vpsi(ig,is1)&
                              -fi*gpotv(id)*(rk(id,1)&
                              +gk(id,ig))*parm%tpiba*c0(ig ,is1)&
                              -fi*0.5_real_8*(gpotv(id)**2)*   c0(ig,is1)
                         C2_vpsi(ig+ncpw%ngw,is1)=c2_vpsi(ig+ncpw%ngw,is1)&
                              -fi*gpotv(id)*(rk(id,1)&
                              -gk(id,ig))*parm%tpiba*c0(ig+ncpw%ngw,is1)&
                              -fi*0.5_real_8*(gpotv(id)**2)*   c0(ig+ncpw%ngw,is1)
                      ENDIF
                   ENDDO
                   ! elisa
                ENDDO
             ENDIF

          ELSE
             !$omp parallel do private(IG,FP,FM) &
             !$omp   shared(ncpw,IS1,FI,IKIND,PSI_p,C0,C2_vpsi)
#ifdef __SR8000
             !poption parallel
#endif
             DO ig=1,ncpw%ngw
                fp=psi_p(nzhs(ig))
                fm=psi_p(indzs(ig))
                C2_vpsi(ig,    is1)=-fi*&
                     (0.5_real_8*parm%tpiba2*hgkp(ig,ikind)*c0(ig,    is1)+fp)
                C2_vpsi(ig+ncpw%ngw,is1)=-fi*&
                     (0.5_real_8*parm%tpiba2*hgkm(ig,ikind)*c0(ig+ncpw%ngw,is1)+fm)
             ENDDO
          ENDIF
          ! EHR]
          IF (geq0)C2_vpsi(1+ncpw%ngw,is1)=CMPLX(0._real_8,0._real_8,kind=real_8)
       ELSE
          fi=f(is1)*0.5_real_8
          IF (fi.EQ.0._real_8.AND..NOT.cntl%tksham) fi=1._real_8
          IF (fi.EQ.0._real_8.AND.cntl%tksham) fi=0.5_real_8
          fip1=0._real_8
          IF (is2.LE.nostat) fip1=f(is2)*0.5_real_8
          IF (fip1.EQ.0._real_8.AND..NOT.cntl%tksham) fip1=1._real_8
          IF (fip1.EQ.0._real_8.AND.cntl%tksham) fip1=0.5_real_8
          IF (prcp_com%akin.GT.1.e-10_real_8) THEN
             xskin=1._real_8/prcp_com%gskin
             !$omp parallel do private(IG,FP,FM,G2,psin,psii) shared(JGW,PSI_p) &
             !$omp   shared(XSKIN,FI,C0,C2_vpsi,IS1) &
             !$omp   shared(NOSTAT,FIP1,IS2)
#ifdef __SR8000
             !poption parallel
#endif
             !dir$ concurrent
             DO ig=1,jgw
                psin=psi_p(nzhs(ig))! access only once
                psii=psi_p(indzs(ig))! these mem locations
                fp=psin+psii
                fm=psin-psii
                g2=parm%tpiba2*(hg(ig)+&
                     prcp_com%gakin*(1._real_8+cp_erf((hg(ig)-prcp_com%gckin)*xskin)))
                C2_vpsi(ig,is1)=-fi*(g2*c0(ig,is1)+&
                     CMPLX(REAL(fp),AIMAG(fm),kind=real_8))
                IF (is2.LE.nostat) C2_vpsi(ig,is2)=-fip1*(g2*c0(ig,is2)+&
                     CMPLX(AIMAG(fp),-REAL(fm),kind=real_8))
             ENDDO
          ELSE
             !$omp parallel do private(IG,FP,FM,psin,psii) shared(JGW,PSI_p) &
             !$omp    shared(C2_vpsi,C0,FI,IS1,IS2,NOSTAT,FIP1)
#ifdef __SR8000
             !poption parallel
#endif
             DO ig=1,jgw
                psin=psi_p(nzhs(ig))! access only once
                psii=psi_p(indzs(ig))! these mem locations
                fp=psin+psii
                fm=psin-psii
                C2_vpsi(ig,is1)=-fi*((parm%tpiba2*hg(ig))*c0(ig,is1)+&
                     CMPLX(REAL(fp),AIMAG(fm),kind=real_8))
                IF (is2.LE.nostat) C2_vpsi(ig,is2)=-fip1*&
                     ((parm%tpiba2*hg(ig))*&
                     c0(ig,is2)+CMPLX(AIMAG(fp),-REAL(fm),kind=real_8))
             ENDDO
          ENDIF
       ENDIF
       __NVTX_TIMER_STOP
    ENDDO                     ! Loop over IA
    !$omp end do

    !vw free local comm for each stream/thread
    CALL mp_comm_free ( fft_comm )

    !$omp end parallel

    __NVTX_TIMER_START ( procedureN//'_c' )

    CALL thread_views_finalize ( thread_views )

    IF( cp_cuda_env%use_fft ) THEN
       CALL tiset(procedureN//'_gpu_fin',isub4)
       CALL cuda_dealloc_host(psis)
#ifndef _HASNT_OMP_SET_NESTED
       !$ CALL omp_set_nested( nested_orig )
#endif
       !$ CALL omp_set_max_active_levels( max_level_orig )
       IF( .NOT. copy_data_to_host ) THEN
          ! TODO acm: the init+allocation/deallocation should be done outside (once for entire code)
          CALL cp_cuvpsi_finalize_and_dealloc( cp_cuvpsi )
       ENDIF
       CALL tihalt(procedureN//'_gpu_fin',isub4)
    ENDIF

    ! Special terms for MGGA. Parallelised over states, call needs
    ! to be before redist & sum to avoid useless allocation calls
    !
    IF (cntl%ttau) CALL vtaupsi(c0,c2_vpsi,f,psi,nstate,ispin)
    !
    ! redistribute C2_vpsi over the groups if needed
    !
    IF (redist_c2) THEN
       CALL tiset(procedureN//'_grps_b',isub3)
       CALL cp_grp_redist(C2_vpsi,nkpt%ngwk,nstate)
       CALL tihalt(procedureN//'_grps_b',isub3)
    ENDIF

    !
    ! add up the vpsi contribution to C2
    !
    CALL add_wfn(nkpt%ngwk,nstate,zone,C2_vpsi,nkpt%ngwk,c2,nkpt%ngwk)

    DEALLOCATE(C2_vpsi,stat=ierr)
    IF (ierr.NE.0) CALL stopgm(procedureN,'Deallocation problem',&
         __LINE__,__FILE__)

    IF (tkpts%tkpnt) CALL c_clean(c2,nstate,ikind)
    ! SPECIAL TERMS FOR LSE METHODS
    IF (lspin2%tlse) CALL vpsi_lse(c0,c2,f,vpot,psi,nstate,.TRUE.)
    !
    __NVTX_TIMER_STOP
    __NVTX_TIMER_STOP
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE vpsi

  ! ==================================================================
  SUBROUTINE vpsimt(c0,c2,f,vpot,psi,nstate,ispin,swfun)
    ! ==================================================================
    ! == No K-Points and no LSE  VERSION OF VPSI.                     ==
    ! == EXACTLY LIKE VPSI BUT NO CONTRIBUTION FROM KINETIC ENERGY    ==
    ! ==--------------------------------------------------------------==
    ! == VPOT:   IN INPUT POTENTIAL                                   ==
    ! == ISPIN:  Need with LSD option for diagonalization scheme      ==
    ! ==         dimension of VPOT(NNR1,ISPIN)                        ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(:,:), c2(:,:)
    REAL(real_8)                             :: f(:)
    REAL(real_8), TARGET                     :: vpot(:,:)
    COMPLEX(real_8)                          :: psi(:)
    INTEGER                                  :: nstate, ispin
    LOGICAL                                  :: swfun

    CHARACTER(*), PARAMETER                  :: procedureN = 'vpsimt'

    COMPLEX(real_8)                          :: c0tmp1, c0tmp2, fm, fp, psii, &
                                                psin
    INTEGER                                  :: i, ia, ib, ibb, ifft, ig, &
                                                is1, is2, isub, l, lead, &
                                                leadx, lspin, njump, nnrx, &
                                                nostat, nsta
    REAL(real_8)                             :: fi, fip1

!(nkpt%ngwk,nstate)
! Variables
! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    ! Other illegal combinations already caught in VPSI
    IF (lspin2%tlse.AND.(lspin2%tcas22.OR.lspin2%tpenal)) THEN
       CALL stopgm('VPSIMT','LSE (CAS22, PENALTY) NOT IMPLEMENTED',&
            __LINE__,__FILE__)
    ENDIF
    IF (tkpts%tkpnt) CALL stopgm('VPSIMT','K-POINTS NOT IMPLEMENTED',&
         __LINE__,__FILE__)
    IF (group%nogrp.EQ.1) THEN
       lead  = fpar%kr1s*parai%ngrays
       leadx = fpar%nnr1
       ifft  = 1
       nnrx=fpar%nnr1
    ELSE
       CALL stopgm(procedureN,'this doesnt work anymore....',&
            __LINE__,__FILE__)
       lead  = fpar%kr1s*ngrm
       leadx = fpar%krx*fpar%kr2s*fpar%kr3s
       ifft=2
       nnrx=fpar%krx*fpar%kr2s*fpar%kr3s
    ENDIF
    njump=2*group%nogrp
    IF (lspin2%tlse) THEN
       nostat = clsd%ialpha-1
    ELSE
       nostat = nstate
    ENDIF
    DO ia=1,nostat,njump
       CALL zeroing(psi)!,maxfft*group%nogrp)
       ! ==----------------------------------------------------------==
       ! == Store the wavefunctions in array PSI                     ==
       ! == in the following way:                                    ==
       ! == Since we use the Gamma point only, the wavefunctions in  ==
       ! == real space are real. This fact is used to save time by   ==
       ! == using complex arrays where the wavefunctions of          ==
       ! == 2 different states (i and i+1) are combined, and         ==
       ! == performing complex Fourier Transforms.                   ==
       ! == Warning: If you use K-points other than Gamma, or zone   ==
       ! == boundary, this symmetry is broken and the trick is not   ==
       ! == valid anymore. First step is to build the array PSI.     ==
       ! ==----------------------------------------------------------==
       nsta=MIN((nostat-ia+2)/2,group%nogrp)
       DO ib=1,nsta
          i=ia+2*(ib-1)
          ibb=(ib-1)*lead
          is1 = i
          is2 = i+1
          IF (is2.GT.nostat) THEN
             !CDIR NODEP
#ifdef __SR8000
             !poption parallel
#endif
#if defined(__VECTOR)
             !$omp parallel do private(IG)
             DO ig=1,ncpw%ngw
                psi(nzhs(ig)+ibb)=c0(ig,is1)
                psi(indzs(ig)+ibb)=CONJG(c0(ig,is1))
             ENDDO
#else
             !$omp parallel do private(IG,C0TMP1) schedule(static)
             DO ig=1,ncpw%ngw
                c0tmp1=c0(ig,is1)
                psi(nzhs(ig)+ibb)=c0tmp1
                psi(indzs(ig)+ibb)=CONJG(c0tmp1)
             ENDDO
#endif
             IF (geq0) psi(nzhs(1)+ibb)=c0(1,is1)
          ELSE
             !CDIR NODEP
#ifdef __SR8000
             !poption parallel
#endif
#if defined(__VECTOR)
             !$omp parallel do private(IG)
             DO ig=1,ncpw%ngw
                psi(nzhs(ig)+ibb)=c0(ig,is1)+uimag*c0(ig,is2)
                psi(indzs(ig)+ibb)=CONJG(c0(ig,is1))+&
                     uimag*CONJG(c0(ig,is2))
             ENDDO
#else
             !$omp parallel do private(IG,C0TMP1,C0TMP2) schedule(static)
             DO ig=1,ncpw%ngw
                c0tmp1=c0(ig,is1)
                c0tmp2=c0(ig,is2)
                psi(nzhs(ig)+ibb)=c0tmp1+uimag*c0tmp2
                psi(indzs(ig)+ibb)=CONJG(c0tmp1)+uimag*CONJG(c0tmp2)
             ENDDO
#endif
             IF (geq0) psi(nzhs(1)+ibb)=c0(1,is1)+uimag*c0(1,is2)
          ENDIF
       ENDDO
       ! ==----------------------------------------------------------==
       ! == Transform the wavefunction to real space                 ==
       ! ==----------------------------------------------------------==
       IF (ifft.EQ.1) THEN
          CALL  invfftn(psi,.TRUE.,parai%allgrp)
       ELSE
          CALL stopgm("THIS","SG_INVFFT NOT AVAILABLE ANYMORE",&
               __LINE__,__FILE__)
       ENDIF
       ! ==------------------------------------------------------------==
       ! == Apply the potential (V), which acts in real space.         ==
       ! ==------------------------------------------------------------==
       IF (cntl%tlsd.AND.ispin.EQ.2) THEN
          i=ia
          IF (i.EQ.spin_mod%nsup) THEN
#ifdef __SR8000
             !poption parallel
#endif
             !$omp parallel do private(L) schedule(static)
             DO l=1,nnrx
                psi(l)=vpot(l,1)*REAL(psi(l))&
                     +uimag*vpot(l,2)*AIMAG(psi(l))
             ENDDO
          ELSE
             lspin=1
             IF (i.GT.spin_mod%nsup) lspin=2
             !$omp parallel do private(L) schedule(static)
             DO l=1,nnrx
                psi(l)=vpot(l,lspin)*psi(l)
             ENDDO
          ENDIF
       ELSE
          !$omp parallel do private(L) schedule(static)
          DO l=1,nnrx
             psi(l)=vpot(l,1)*psi(l)
          ENDDO
       ENDIF
       ! ==------------------------------------------------------------==
       ! == Back transform to reciprocal space the product V.PSI       ==
       ! ==------------------------------------------------------------==
       IF (ifft.EQ.1) THEN
          CALL  fwfftn(psi,.TRUE.,parai%allgrp)
       ELSE
          CALL stopgm("THIS","SG_FWFFT NOT AVAILABLE ANYMORE",&
               __LINE__,__FILE__)
       ENDIF
       ! ==------------------------------------------------------------==
       ! == Decode the combination of wavefunctions (now multiplied by ==
       ! == the potential), back to those of states i and i+1. Multiply==
       ! == also by the occupation numbers f(i) and other constants in ==
       ! == order to obtain the force on the electronic degrees of     ==
       ! == freedom, stored in array FC.                               ==
       ! == Here we also add the force from the kinetic energy.        ==
       ! ==------------------------------------------------------------==
#ifdef __SR8000
       !poption parallel
       !poption tlocal(FIP1)
       !voption indep(NZHS,INDZS,C2)
#endif
       DO ib=1,nsta
          i=ia+2*(ib-1)
          ibb=(ib-1)*lead
          is1 = i
          is2 = i+1
          fi=f(is1)*0.5_real_8
          IF (fi.EQ.0._real_8) fi=1._real_8
          fip1=0._real_8
          IF (is2.LE.nostat) fip1=f(is2)*0.5_real_8
          IF (fip1.EQ.0._real_8) fip1=1._real_8
          !$omp parallel do private(IG,psin,psii,FP,FM) schedule(static)
          DO ig=1,ncpw%ngw
             psin=psi(nzhs(ig)+ibb)
             psii=psi(indzs(ig)+ibb)
             fp=psin+psii
             fm=psin-psii
             ! mb            FP=PSI(NZHS(IG)+IBB)+PSI(INDZS(IG)+IBB)
             ! mb            FM=PSI(NZHS(IG)+IBB)-PSI(INDZS(IG)+IBB)
             c2(ig,is1)=c2(ig,is1)-fi*CMPLX(REAL(fp),AIMAG(fm),kind=real_8)
             IF (is2.LE.nostat)&
                  c2(ig,is2)=c2(ig,is2)-fip1*CMPLX(AIMAG(fp),-REAL(fm),kind=real_8)
          ENDDO
       ENDDO
    ENDDO                     ! Loop over IA
    ! Reestablish initial conditions for the potential
    ! SPECIAL TERMS FOR LSE METHODS
    IF (lspin2%tlse) THEN
       CALL vpsi_lse(c0,c2,f,vpot,psi,nstate,.FALSE.)
    ENDIF
    ! META FUNCTIONALS NEED SPECIAL TERM
    IF (cntl%ttau) THEN
       IF (swfun) CALL switch_functionals
       IF (cntl%ttau) CALL vtaupsi(c0,c2,f,psi,nstate,ispin)
       IF (swfun) CALL switch_functionals
    ENDIF
    !
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE vpsimt
  ! ==================================================================
  SUBROUTINE movepsid(a)
    COMPLEX(real_8)                          :: a(:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'movepsid'

    COMPLEX(real_8), ALLOCATABLE             :: b(:)
    INTEGER                                  :: ierr, ig

! ==--------------------------------------------------------------==

    ALLOCATE(b(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    !$omp parallel do private(IG) schedule(static)
    DO ig=1,ncpw%nhg
       b(ig)=a(nzh(ig))
    ENDDO
    CALL zeroing(a)!,maxfftn)
    !$omp parallel do private(IG) schedule(static)
    DO ig=1,jhg
       a(nzff(ig))=b(ig)
       a(inzf(ig))=CONJG(b(ig))
    ENDDO
    IF (geq0) a(nzff(1))=b(1)
    DEALLOCATE(b,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE movepsid
  ! ==================================================================

END MODULE vpsi_utils
