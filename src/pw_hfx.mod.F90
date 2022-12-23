#include "cpmd_global.h"

#define PRINT_GROUP_INFOS .FALSE.

! ================================================================== 
MODULE pw_hfx
  USE cnst,                            ONLY: uimag
  USE cp_grp_utils,                    ONLY: cp_grp_redist
  USE cppt,                            ONLY: scgx
  USE dotp_utils,                      ONLY: dotp
  USE error_handling,                  ONLY: stopgm
  USE fft,                             ONLY: inzf,&
                                             inzs,&
                                             jgw,&
                                             jhg,&
                                             llr1,&
                                             nzff,&
                                             nzfs
  USE fft_maxfft,                      ONLY: maxfftn
  USE fftmain_utils,                   ONLY: fwfftn,&
                                             invfftn
  USE fftnew_utils,                    ONLY: setfftn
  USE func,                            ONLY: func1,&
                                             func3
  USE geq0mod,                         ONLY: geq0
  USE hfx_utils,                       ONLY: get_wannier_separation
  USE hfxmod,                          ONLY: hfxc3,&
                                             hfxc4,&
                                             hfxc5,&
                                             ipoolhfx
  USE kinds,                           ONLY: int_8,&
                                             real_8
  USE kpts,                            ONLY: tkpts
  USE machine,                         ONLY: m_cputime,&
                                             m_flush
  USE min_heap,                        ONLY: heap_fill,&
                                             heap_get_first,&
                                             heap_new,&
                                             heap_release,&
                                             heap_reset_first,&
                                             heap_t
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum
  USE parac,                           ONLY: parai,&
                                             paral
  USE part_1d,                         ONLY: part_1d_get_elem,&
                                             part_1d_nbr_elems
  USE pw_hfx_input_cnst,               ONLY: hfx_dist_block_cyclic,&
                                             hfx_dist_dynamic
  USE ropt,                            ONLY: infi,&
                                             infw,&
                                             iteropt
  USE scex_utils,                      ONLY: scex,&
                                             scex_lambda,&
                                             scex_ID_parent,&
                                             scex_ID_scaled
  USE spin,                            ONLY: lspin2,&
                                             spin_mod
  USE state_utils,                     ONLY: add_wfn,&
                                             copy_im_to_im,&
                                             copy_im_to_re,&
                                             copy_re_to_im,&
                                             copy_re_to_re,&
                                             set_psi_1_state_g,&
                                             set_psi_2_states_g
  USE system,                          ONLY: cntl,&
                                             group,&
                                             ncpw,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset
! ================================================================== 
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  ! ==--------------------------------------------------------------==

  PRIVATE

  PUBLIC :: hfx_new
  PUBLIC :: block_invfft_old
  PUBLIC :: build_block_info
  PUBLIC :: build_dist
  PUBLIC :: release_dist
  PUBLIC :: iterator_start
  PUBLIC :: iterator_blocks
  PUBLIC :: iterator_next_block
  PUBLIC :: ind2sub_rect
  PUBLIC :: dist_t
  PUBLIC :: iterator_t

  INTEGER,DIMENSION(:,:),ALLOCATABLE,PRIVATE :: tick_blks
  INTEGER(int_8),PRIVATE :: nbr_rwfn_precomputed = 0
  INTEGER(int_8),PRIVATE :: nbr_integrals = 0

  TYPE :: iterator_t
     INTEGER :: blk_count,max_blk_count
     INTEGER,DIMENSION(:),POINTER :: dist
  END TYPE iterator_t

  TYPE :: dist_t
     INTEGER :: TYPE,nblks_loc,nblks
     INTEGER,DIMENSION(:),POINTER :: dist
  END TYPE dist_t

  TYPE :: pp_t
     ! pair of pairs type
     !
     ! nbr of pair of pairs
     INTEGER :: n_pp=0
     !
     !
     INTEGER :: r_offset, c_offset
     !
     ! ptr to the pairs
     INTEGER, ALLOCATABLE, DIMENSION(:) :: pp_ptr
     !
     ! the pair of pairs ( ((r1,c1)), ((r2,c1),(r4,c1)), ... )
     INTEGER, ALLOCATABLE, DIMENSION(:,:) :: pp
  END TYPE pp_t

  TYPE :: list_t
     !
     ! nbr entries in the list
     INTEGER :: n=0
     !
     ! list of absolute positions
     INTEGER, ALLOCATABLE, DIMENSION(:) :: a
     !
     ! list of relative positions
     INTEGER, ALLOCATABLE, DIMENSION(:) :: r
     !
     ! where entries are stored (>0 : real, <0 : complex)
     ! the size is at most max_block_size
     INTEGER, ALLOCATABLE, DIMENSION(:) :: stored_in
  END TYPE list_t


  TYPE :: thresh_t
     INTEGER :: int_local_count=0
     ! decide if the integrals need to be recomputed
     LOGICAL :: init_ints
     ! local data, allocated with thresh_init/clean_local
     INTEGER,DIMENSION(:),ALLOCATABLE :: new_vals!, bin_vals
     REAL(real_8),DIMENSION(:),ALLOCATABLE :: acc_vals!, max_vals

     ! global data, should be allocated outside
     REAL(real_8),DIMENSION(:),ALLOCATABLE :: int_vals
  END TYPE thresh_t


  TYPE(thresh_t), PRIVATE, SAVE :: thresh

CONTAINS

  ! ================================================================== 
  SUBROUTINE hfx_new(c0,c2,f,psic,nstate,ehfx,vhfx,redist_c2)
    ! ==--------------------------------------------------------------==
    ! - Refactored: New procedures for ehfx, potentials, pair densities
    !               Scaled exact exchange (ScEX), some prettifications
    !                             21.11.2018 - M. P. Bircher @ LCBC/EPFL
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8)                          :: c0(:,:), c2(:,:)
    REAL(real_8)                             :: f(:)
    COMPLEX(real_8)                          :: psic(:)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: ehfx, vhfx
    LOGICAL, INTENT(in)                      :: redist_c2

    CHARACTER(*), PARAMETER                  :: procedureN = 'hfx_new'
    COMPLEX(real_8), PARAMETER               :: zone = (1.0_real_8,0.0_real_8)

    COMPLEX(real_8), ALLOCATABLE, &
      DIMENSION(:)                           :: psi_row_pack
    COMPLEX(real_8), ALLOCATABLE, &
      DIMENSION(:, :)                        :: C2_hfx, psi_col, psi_row, &
                                                vpotg
    INTEGER :: blk, col_blk, col_offset, col_size, hfx_dist, i, ierr, ispin, &
      isub, isub6, max_block_size, max_nblks, nblks, nbr_blk, nspins, nst, &
      nstates(2), row_blk, row_offset, row_size, st_offst, tick_1, tick_2
    INTEGER(int_8)                           :: i8_buff(2)
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: new_tick_blks, offsets, sizes
    REAL(real_8)                             :: r8_buff(2), t_outer, &
                                                t_outer_1, t_outer_2
    REAL(real_8), ALLOCATABLE, &
      DIMENSION(:, :)                        :: vpotr
    TYPE(dist_t)                             :: dist
    TYPE(iterator_t)                         :: iter
    TYPE(list_t)                             :: c_list, r_list
    TYPE(pp_t)                               :: pp

    ! KPT replace NGW -> NGWK (maybe not C2....)
    ! ==--------------------------------------------------------------==

    ehfx = 0.0_real_8
    vhfx = 0.0_real_8

    IF (func1%mhfx /= 1) RETURN

    ! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)

    ! ==--------------------------------------------------------------==

    IF (tkpts%tkpnt) CALL stopgm(procedureN,'K-POINTS NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    IF (lspin2%tlse) CALL stopgm(procedureN,'LSE NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    IF (group%nogrp > 1) CALL stopgm(procedureN,'TASK GROUPS NOT IMPLEMENTED',& 
         __LINE__,__FILE__)

    ! ==--------------------------------------------------------------==

    IF (cntl%use_scaled_hfx) THEN
       IF (.NOT. scex%init) CALL stopgm(procedureN,'ScEX requested, but not initialised',&
                                        __LINE__,__FILE__)
       CALL setfftn(scex_ID_parent)
    ELSE
       CALL setfftn(ipoolhfx)
    ENDIF

    ! ==--------------------------------------------------------------==

    nbr_rwfn_precomputed = 0
    nbr_integrals = 0

    nstates = 0
    IF (cntl%tlsd) THEN
       nspins = 2
       nstates(1) = spin_mod%nsup
       nstates(2) = spin_mod%nsdown
    ELSE
       nspins = 1
       nstates(1) = nstate
    ENDIF

    ! ==--------------------------------------------------------------==

    !
    ! initialize integral buffer
    !
    IF (hfxc3%twscr) THEN
       IF (.NOT.ALLOCATED(thresh%int_vals)) THEN
          ! this should be allocated before and deallocated
          ALLOCATE(thresh%int_vals( nstate*(nstate+1)/2 ),stat=ierr)
          IF (ierr /= 0) CALL stopgm(procedureN,'allocation problem',& 
               __LINE__,__FILE__)
          CALL zeroing(thresh%int_vals)
       ENDIF
    ENDIF

    CALL thresh_init_local(thresh,nstate)

    ! ==--------------------------------------------------------------==

    hfx_dist = hfxc5%hfx_distribution
    max_block_size = hfxc5%hfx_block_size
    IF (max_block_size < 1) CALL stopgm( procedureN,&
         'HFX_BLOCK_SIZE should be greater than zero' ,& 
         __LINE__,__FILE__)

    ! get the max number of block pairs
    nbr_blk = CEILING(REAL(MAXVAL(nstates),kind=real_8)/REAL(max_block_size,kind=real_8))
    max_nblks = nbr_blk * ( nbr_blk + 1 ) / 2

    ! ==--------------------------------------------------------------==

    IF (.NOT.ALLOCATED(tick_blks)) THEN
       ALLOCATE(tick_blks(max_nblks,nspins),stat=ierr)
       IF (ierr /= 0) CALL stopgm( procedureN, 'Allocation problem' ,& 
            __LINE__,__FILE__)
       tick_blks(:,:) = 1
    ELSE
       IF( SIZE(tick_blks,1) /= max_nblks ) THEN
          DEALLOCATE(tick_blks,stat=ierr)
          IF (ierr /= 0) CALL stopgm( procedureN, 'Deallocation problem' ,&
               __LINE__,__FILE__)
          ALLOCATE(tick_blks(max_nblks,nspins),stat=ierr)
          IF (ierr /= 0) CALL stopgm( procedureN, 'Allocation problem' ,&
               __LINE__,__FILE__)
          tick_blks(:,:) = 1
       ENDIF
    ENDIF

    ALLOCATE(vpotg(jhg,2),vpotr(llr1,2),&
         psi_row(maxfftn,(max_block_size+1)/2),&
         psi_col(maxfftn,(max_block_size+1)/2),&
         C2_hfx(ncpw%ngw,nstate),&
         sizes(nbr_blk),offsets(nbr_blk), & ! FIXME change to max_nbr_blk
         new_tick_blks(max_nblks), &
         psi_row_pack(maxfftn), &
         stat=ierr)
    IF (ierr /= 0) CALL stopgm( procedureN, 'Allocation problem' ,& 
         __LINE__,__FILE__)

    CALL pp_init( pp, max_block_size )
    CALL list_init( r_list, max_block_size)
    CALL list_init( c_list, max_block_size)

    CALL zeroing(C2_hfx)
    CALL zeroing(psi_row_pack)

    ! ==--------------------------------------------------------------==

    t_outer_1=m_cputime()

    IF (cntl%use_scaled_hfx) THEN
       CALL process_scaled_blocks()
    ELSE
       CALL process_blocks()
    ENDIF

    t_outer_2=m_cputime()
    t_outer = t_outer_2 - t_outer_1

    ! ==--------------------------------------------------------------==

    DEALLOCATE(vpotg,vpotr,psi_row,psi_col,new_tick_blks,&
         sizes,offsets,psi_row_pack,STAT=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN,'Deallocation problem',& 
         __LINE__,__FILE__)

    CALL pp_destroy( pp )
    CALL list_destroy( r_list )
    CALL list_destroy( c_list )


    !
    ! compute new radius if needed
    !
    CALL thresh_get_new_radius(thresh)

    !
    ! finalize integral buffer
    !
    CALL thresh_clean_local(thresh)

    ! ==--------------------------------------------------------------==

    ! 
    ! redistribute EHFX and C2_hfx over the groups if needed
    ! 
    IF (parai%cp_nogrp > 1) THEN
       CALL tiset(procedureN//'_b',isub6)
       CALL mp_sum(ehfx,parai%cp_inter_grp)
       IF (redist_c2) CALL cp_grp_redist(C2_hfx,ncpw%ngw,nstate)
       CALL tihalt(procedureN//'_b',isub6)
    ENDIF

    ! 
    ! add up the hfx contribution to C2
    ! 
    CALL add_wfn(jgw,nstate,zone,C2_hfx,ncpw%ngw,c2,ncpw%ngw)


    ! Energy
    ehfx = ehfx*parm%omega

    ! Potential
    DO i = 1,nstate
       vhfx = vhfx + dotp(ncpw%ngw,c0(:,i),c2(:,i))
    ENDDO

    r8_buff(1) = EHFX; r8_buff(2) = VHFX
    CALL mp_sum(r8_buff,2,parai%allgrp)
    EHFX = r8_buff(1); VHFX = r8_buff(2)

    ! ==--------------------------------------------------------------==

    ! Redistribute timings if needed
    IF (parai%cp_nogrp > 1.AND.hfx_dist /= hfx_dist_block_cyclic) THEN
       CALL mp_sum(tick_blks,max_nblks*nspins,parai%cp_inter_grp)
    ENDIF
    ! ==--------------------------------------------------------------==

    IF (PRINT_GROUP_INFOS) THEN
       IF (parai%me == 0) THEN
          WRITE(6,'(1X,3(A,I0),A,F0.2,A,F0.2)')&
               procedureN//'| group ',parai%cp_inter_me,&
               ' computed ',nbr_integrals,&
               ' integrals, precomputed ',nbr_rwfn_precomputed,&
               ', t_loc ',t_outer,&
               ', t_per_1k_ints ',1.0e3_real_8*t_outer/(REAL(nbr_integrals,kind=real_8)+1.0e-6_real_8) ! to avoid NANs
          CALL m_flush(6)
       ENDIF
       IF (parai%cp_nogrp > 1) THEN
          i8_buff(1) = nbr_integrals
          i8_buff(2) = nbr_rwfn_precomputed
          CALL mp_sum(i8_buff,2,parai%cp_inter_grp)
          nbr_integrals = i8_buff(1)
          nbr_rwfn_precomputed = i8_buff(2)
          IF (parai%cp_me == 0) THEN
             WRITE(6,'(1X,3(A,I0))')&
                  procedureN//'| all the groups computed ',&
                  nbr_integrals,' integrals, precomputed ',&
                  nbr_rwfn_precomputed
             CALL m_flush(6)
          ENDIF
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==

    DEALLOCATE(C2_hfx,stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN,'Deallocation problem',& 
         __LINE__,__FILE__)

    ! ==--------------------------------------------------------------==

    CALL tihalt(procedureN,isub)

    CONTAINS

    ! ================================================================== 
    SUBROUTINE process_blocks()
      ! ==--------------------------------------------------------------==
      ! Conventional algorithm for full HFX; no storage of real-space
      ! pair potentials.

      DO ispin = 1,nspins

         nst = nstates(ispin)
         st_offst = 0
         IF (ispin == 2) st_offst = nstates(1)

         ! get the number of block pairs
         nbr_blk = CEILING(REAL(nst,kind=real_8)/REAL(max_block_size,kind=real_8))
         nblks = nbr_blk * ( nbr_blk + 1 ) / 2

         CALL build_block_info(nst,nbr_blk,max_block_size,offsets,sizes)
         CALL build_dist(dist,nblks,hfx_dist,tick_blks(:,ispin))

         new_tick_blks(:) = 0
         CALL iterator_start(iter,dist)
         DO WHILE (iterator_blocks(iter))
            CALL iterator_next_block(iter,blk)
            CALL ind2sub(blk,row_blk,col_blk)
            row_offset = offsets(row_blk) + st_offst
            col_offset = offsets(col_blk) + st_offst
            row_size = sizes(row_blk); col_size = sizes(col_blk)

            tick_1 = get_ticks()

            CALL block_filter(thresh,row_offset,row_size,col_offset,col_size,pp,r_list,c_list)

            CALL block_invfft(c0,psi_row,r_list)
            CALL block_invfft(c0,psi_col,c_list)
            !call block_invfft_old(c0,psi_row,row_offset,row_size)
            !call block_invfft_old(c0,psi_col,col_offset,col_size)

            CALL hfx_compute_block_new( C2_hfx, ehfx, f, vpotg, vpotr, psic, pp, &
                 psi_row, psi_col, r_list, c_list, psi_row_pack, thresh )

            !call hfx_compute_block( C2_hfx, ehfx, f, vpotg, vpotr, psic,&
            !     psi_row, row_offset, row_size,&
            !     psi_col, col_offset, col_size)

            tick_2 = get_ticks()

            new_tick_blks(blk) = tick_2 - tick_1

         ENDDO
         CALL release_dist(dist)
         tick_blks(:,ispin) = new_tick_blks(:)
      ENDDO

      ! ==--------------------------------------------------------------==
    END SUBROUTINE process_blocks
    ! ================================================================== 
    SUBROUTINE process_scaled_blocks()
      ! ==--------------------------------------------------------------==
      ! Iterator for coordinate-scaled exchange; with storage of real-space
      ! pair potentials.
      !                             20.11.2018 - M. P. Bircher @ LCBC/EPFL

      COMPLEX(real_8), DIMENSION(:,:), &
                       ALLOCATABLE :: C2_in_real_space

      CHARACTER(len=*), PARAMETER  :: subprocedureN = procedureN//'_process_scaled_blocks'

      INTEGER         :: istate, ig, ir
      INTEGER         :: blk_state, blk_resid

      COMPLEX(real_8) :: fp, fm

      !
      ! Allocate auxiliary quantities for scaled exact exchange
      !
      ALLOCATE(C2_in_real_space(maxfftn,nstate),&
               stat=ierr)
      IF (ierr /= 0) CALL stopgm(subprocedureN,'Allocation problem',&
                                 __LINE__,__FILE__)
      CALL zeroing(C2_in_real_space)

      CALL scex%start_density_scaling()

      DO ispin = 1,nspins

         nst = nstates(ispin)
         st_offst = 0
         IF (ispin == 2) st_offst = nstates(1)

         ! get the number of block pairs
         nbr_blk = CEILING(REAL(nst,kind=real_8)/REAL(max_block_size,kind=real_8))
         nblks = nbr_blk * ( nbr_blk + 1 ) / 2

         CALL build_block_info(nst,nbr_blk,max_block_size,offsets,sizes)
         CALL build_dist(dist,nblks,hfx_dist,tick_blks(:,ispin))

         new_tick_blks(:) = 0
         CALL iterator_start(iter,dist)
         DO WHILE (iterator_blocks(iter))
            CALL iterator_next_block(iter,blk)
            CALL ind2sub(blk,row_blk,col_blk)
            row_offset = offsets(row_blk) + st_offst
            col_offset = offsets(col_blk) + st_offst
            row_size = sizes(row_blk); col_size = sizes(col_blk)

            tick_1 = get_ticks()

            CALL block_filter(thresh,row_offset,row_size,col_offset,col_size,pp,r_list,c_list)

            !
            ! Perform batch of operations on full grid
            !
            CALL setfftn(scex_ID_parent)
            CALL block_invfft(c0,psi_row,r_list)
            CALL block_invfft(c0,psi_col,c_list)

            !
            ! Calculate pair potentials on auxiliary mesh
            !
            CALL setfftn(scex_ID_scaled)
            CALL hfx_compute_block_new( C2_in_real_space, ehfx, f, vpotg, vpotr, psic, pp, &
                 psi_row, psi_col, r_list, c_list, psi_row_pack, thresh )

            tick_2 = get_ticks()

            new_tick_blks(blk) = tick_2 - tick_1

         ENDDO
         CALL release_dist(dist)
         tick_blks(:,ispin) = new_tick_blks(:)
      ENDDO

      !
      ! Calculate C2 on the full wavefunction cutoff
      !
      CALL scex%switch_density_scaling()
      CALL setfftn(scex_ID_parent)

      blk_state = nstate / 2
      blk_resid = nstate - 2*blk_state
      !
      ! Pack two states into the FFT if possible
      !
      DO istate=1,blk_state
         DO ir=1,scex%llr1
            C2_in_real_space(ir,istate) = REAL(C2_in_real_space(ir,istate)) &
                                          + uimag*REAL(C2_in_real_space(ir,istate+blk_state))
         ENDDO
         CALL scex%undo_density_scaling(C2_in_real_space(:,istate))

         CALL fwfftn(C2_in_real_space(:,istate),.TRUE.,parai%allgrp)
         !$omp parallel do private(IG,FP,FM)
         DO ig=1,jgw
            fp = C2_in_real_space(nzfs(ig),istate) + C2_in_real_space(inzs(ig),istate)
            fm = C2_in_real_space(nzfs(ig),istate) - C2_in_real_space(inzs(ig),istate)
            C2_hfx(ig,istate)           = -CMPLX(REAL(fp),AIMAG(fm),kind=real_8)
            C2_hfx(ig,istate+blk_state) = -CMPLX(AIMAG(fp),-REAL(fm),kind=real_8)
         ENDDO
      ENDDO
      !
      ! Treat residue (HOMO) if number of orbitals is odd
      !
      IF (blk_resid > 0) THEN
         CALL scex%undo_density_scaling(C2_in_real_space(:,nstate))

         CALL fwfftn(C2_in_real_space(:,nstate),.TRUE.,parai%allgrp)
         !$omp parallel do private(IG,FP,FM)
         DO ig=1,jgw
            fp = C2_in_real_space(nzfs(ig),nstate) + C2_in_real_space(inzs(ig),nstate)
            fm = C2_in_real_space(nzfs(ig),nstate) - C2_in_real_space(inzs(ig),nstate)
            C2_hfx(ig,nstate) = -CMPLX(REAL(fp),AIMAG(fm),kind=real_8)
         ENDDO
      ENDIF
      !
      CALL scex%annihilate_density_scaling()

      DEALLOCATE(C2_in_real_space,stat=ierr)
      IF (ierr /= 0) CALL stopgm(subprocedureN,'Deallocation problem',&
                                 __LINE__,__FILE__)

      ! ==--------------------------------------------------------------==
    END SUBROUTINE process_scaled_blocks
    ! ================================================================== 
  END SUBROUTINE hfx_new
  ! ================================================================== 

  ! ================================================================== 
  ! Thresholding and list construction
  ! ================================================================== 

  ! ================================================================== 
  SUBROUTINE thresh_init_local(thresh,nstate)
    ! ==--------------------------------------------------------------==

    TYPE(thresh_t), INTENT(inout)            :: thresh
    INTEGER, INTENT(in)                      :: nstate

    CHARACTER(*), PARAMETER :: procedureN = 'thresh_init_local'

    INTEGER                                  :: geo_iter, ierr
    INTEGER, SAVE                            :: prev_geo_id = -HUGE(0)

    IF (hfxc3%twscr) THEN


       ! decide if reinit intergrals
       ! need to use the *right* geo iter counter (in the cpmd spirit)!
       thresh%init_ints=.FALSE.
       geo_iter=0
       IF(.NOT.cntl%wfopt) geo_iter=infi
       IF (prev_geo_id/=geo_iter) THEN
          prev_geo_id=geo_iter
          thresh%init_ints=.TRUE.
       ENDIF

       thresh%init_ints=MOD(thresh%int_local_count,hfxc5%recomp_two_int_list_every)==0 &
            & .OR.thresh%init_ints


       ALLOCATE(thresh%acc_vals( nstate*(nstate+1)/2 ),&
            thresh%new_vals( nstate*(nstate+1)/2 ), STAT=ierr )
       IF (ierr /= 0) CALL stopgm(procedureN,'allocation problem',& 
            __LINE__,__FILE__)
       CALL zeroing(thresh%acc_vals)
       CALL zeroing(thresh%new_vals)

       IF (paral%io_parent.AND..FALSE.) THEN
          WRITE(6,*) '>>>>>>>>>>     WFOPT=',cntl%wfopt
          WRITE(6,*) '>>>>>>>>>>       NFI=',iteropt%nfi
          WRITE(6,*) '>>>>>>>>>>      INFI=',infi!<
          WRITE(6,*) '>>>>>>>>>>     IINFI=',iteropt%iinfi
          WRITE(6,*) '>>>>>>>>>>      INFW=',infw
          WRITE(6,*) '>>>>>>>>>> init_ints=',thresh%init_ints
       ENDIF

       IF (paral%io_parent.AND..FALSE.)THEN
          WRITE(6,*) 'TWSCR=',hfxc3%twscr
          WRITE(6,*) 'TWFT=',hfxc3%twft,&
               'DWF_INTEGRAL_THRESH=',hfxc4%dwf_integral_thresh
          WRITE(6,*) 'TWFC=',hfxc3%twfc,' DWFC=',hfxc4%dwfc
          WRITE(6,*) 'TDIAW=',hfxc3%tdiaw
       ENDIF


       thresh%int_local_count=thresh%int_local_count+1

    ENDIF


    ! ==--------------------------------------------------------------==
  END SUBROUTINE thresh_init_local
  ! ================================================================== 
  SUBROUTINE thresh_clean_local(thresh)
    ! ==--------------------------------------------------------------==
    TYPE(thresh_t), INTENT(inout)            :: thresh

    CHARACTER(*), PARAMETER :: procedureN = 'thresh_clean_local'

    INTEGER                                  :: ierr

    IF(ALLOCATED(thresh%acc_vals)) THEN
       DEALLOCATE(thresh%acc_vals, STAT=ierr )
       IF (ierr /= 0) CALL stopgm(procedureN,'deallocation problem',& 
            __LINE__,__FILE__)
    ENDIF
    IF(ALLOCATED(thresh%new_vals)) THEN
       DEALLOCATE(thresh%new_vals, STAT=ierr )
       IF (ierr /= 0) CALL stopgm(procedureN,'deallocation problem',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
  END SUBROUTINE thresh_clean_local
  ! ================================================================== 
  SUBROUTINE thresh_get_new_radius(thresh)
    ! ==--------------------------------------------------------------==
    TYPE(thresh_t), INTENT(inout)            :: thresh

    CHARACTER(*), PARAMETER :: procedureN = 'thresh_get_new_radius'
    REAL(real_8), PARAMETER                  :: bin_max = 200.0_real_8, &
                                                bin_range = 0.5_real_8

    INTEGER                                  :: i, ibin, ierr, int_i, int_j, &
                                                n_int
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: bin_vals
    REAL(real_8)                             :: dab, max_dab, old_DWFC
    REAL(real_8), ALLOCATABLE, DIMENSION(:)  :: max_vals

    IF (hfxc3%twscr) THEN

       !
       ! accumulate integral values and new integral
       !
       CALL mp_sum(thresh%acc_vals,SIZE(thresh%acc_vals),parai%cp_grp)
       CALL mp_sum(thresh%new_vals,SIZE(thresh%new_vals),parai%cp_inter_grp)

       ALLOCATE( max_vals( INT( bin_max / bin_range ) + 1 ),&
            bin_vals( INT( bin_max / bin_range ) + 1 ), STAT=ierr )
       IF (ierr /= 0) CALL stopgm(procedureN,'allocation problem',& 
            __LINE__,__FILE__)

       max_vals(:) = 0.0_real_8
       bin_vals(:) = 0
       n_int = 0
       max_dab = 0.0_real_8

       ! The following compiler has a problem with the reduction in the omp section
       ! IBM XL Fortran Advanced Edition for Blue Gene/P, V11.1
       ! Version: 11.01.0000.0011
       !$omp parallel do &
       !$omp             default(none) &
       !$omp             private(i,ibin,int_i,int_j,dab) &
       !$omp             shared(thresh) &
       !$omp             reduction(+:n_int,bin_vals) &
       !$omp             reduction(max:max_vals,max_dab)
       DO i=1,SIZE(thresh%int_vals)
          IF (thresh%new_vals(i) > 0) THEN
             thresh%int_vals(i) = ABS( thresh%acc_vals(i) )
             n_int = n_int + 1
          ENDIF
          CALL ind2sub(i,int_i,int_j)
          CALL get_wannier_separation(int_i,int_j,dab)
          max_dab=MAX(max_dab,dab)
          ibin=INT(dab/bin_range)+1
          max_vals(ibin)=MAX(max_vals(ibin),thresh%int_vals(i))
          bin_vals(ibin)=bin_vals(ibin)+1
       ENDDO
       !$omp end parallel do


       !
       ! find the new radius
       !
       old_DWFC = hfxc4%dwfc
       IF (hfxc3%twfc.AND.hfxc3%twft.AND..NOT.hfxc3%keep_radius_fixed) THEN
          ibin = INT( hfxc4%dwfc / bin_range ) + 1
          IF(ibin<3) CALL stopgm(procedureN,'shouldnt get there: ibin<3',&
               __LINE__,__FILE__)
          IF(ibin>SIZE(max_vals)) CALL stopgm(procedureN,'shouldnt get there: ibin>size(max_vals)',&
               __LINE__,__FILE__)
          IF( max_vals(ibin  ) > hfxc4%dwf_integral_thresh ) hfxc4%dwfc = hfxc4%dwfc + bin_range
          IF( max_vals(ibin-1) < hfxc4%dwf_integral_thresh.AND. &
               max_vals(ibin-2) < hfxc4%dwf_integral_thresh ) hfxc4%dwfc = hfxc4%dwfc - bin_range
       ENDIF

       IF (paral%io_parent.AND..FALSE.) THEN
          WRITE(6,*) 'we have ',n_int,' intergrals (re)set'
          WRITE(6,*) 'new DWFC=',hfxc4%dwfc,' old DWFC=',old_DWFC
          DO i=1,INT(max_dab/bin_range)+1
             WRITE(6,'(2F6.1,E9.2,I8)') (i-1)*bin_range,i*bin_range,max_vals(i),bin_vals(i)
          ENDDO
       ENDIF

!!$    if (paral%io_parent.and..false.) then
!!$       do i=1,size(int_vals,1)
!!$          int_j=ceiling((-1.0_real_8+sqrt(1.0_real_8+8.0_real_8*real(i,kind=real_8)))/2.0_real_8)
!!$          int_i=i-int_j*(int_j-1)/2
!!$          call get_wannier_separation(int_i,int_j,dab)
!!$          write(6,'(A,I0,5E12.4)') 'INT',int_count,dab,&
!!$               wcentx(4,int_i),wcentx(4,int_j),&
!!$               int_vals(i,:)
!!$       enddo
!!$    endif

       DEALLOCATE(max_vals,bin_vals,stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN,'deallocation problem',& 
            __LINE__,__FILE__)

    ENDIF

    ! ==--------------------------------------------------------------==
  END SUBROUTINE thresh_get_new_radius
  ! ================================================================== 
  SUBROUTINE build_dist(dist,nblks,TYPE,dist_costs)
    ! ==--------------------------------------------------------------==
    TYPE(dist_t)                             :: dist
    INTEGER, INTENT(in)                      :: nblks, TYPE
    INTEGER, DIMENSION(:), INTENT(in)        :: dist_costs

    CHARACTER(*), PARAMETER                  :: procedureN = 'build_dist'

    INTEGER                                  :: blk, grp, i, ierr, root_bin, &
                                                root_cost
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: nblks_per_bin
    TYPE(heap_t)                             :: heap

    dist%TYPE = TYPE
    dist%nblks_loc = 0
    dist%nblks = nblks
    NULLIFY(dist%dist)

    SELECT CASE(TYPE)
    CASE(hfx_dist_block_cyclic)

       dist%nblks_loc = part_1d_nbr_elems(dist%nblks,&
            parai%cp_inter_me,parai%cp_nogrp)
       ALLOCATE(dist%dist(dist%nblks_loc),stat=ierr)
       IF (ierr /= 0)CALL stopgm(procedureN,'Allocation problem',& 
            __LINE__,__FILE__)
       DO i = 1,dist%nblks_loc
          blk = part_1d_get_elem(i,parai%cp_inter_me,parai%cp_nogrp)
          dist%dist(i) = blk
       ENDDO

    CASE(hfx_dist_dynamic)

       ALLOCATE( nblks_per_bin( parai%cp_nogrp ), stat=ierr )
       IF (ierr /= 0)CALL stopgm(procedureN,'Allocation problem',& 
            __LINE__,__FILE__)

       CALL heap_new( heap, parai%cp_nogrp )
       ! 
       ! Get the local number of blocks
       ! 
       CALL heap_fill( heap, (/( grp, grp = 1,parai%cp_nogrp )/),&
            (/( 0, grp = 1,parai%cp_nogrp )/) )

       nblks_per_bin(:) = 0
       IF (parai%me == parai%source) THEN
          DO blk = 1,nblks
             CALL heap_get_first(heap,root_bin,root_cost)
             nblks_per_bin(root_bin) = nblks_per_bin(root_bin) + 1
             root_cost = root_cost + dist_costs(blk)
             CALL heap_reset_first(heap,root_cost)
          ENDDO
       ENDIF

       CALL mp_bcast(nblks_per_bin,SIZE(nblks_per_bin),&
            parai%source,parai%allgrp)

       IF (PRINT_GROUP_INFOS) THEN
          IF (parai%cp_me == 0)WRITE(6,'(2A,2(I0,1X),F7.2)')procedureN,&
               ' nblks_per_bin min/max/avg : ',&
               MINVAL(nblks_per_bin),MAXVAL(nblks_per_bin),&
               SUM(REAL(nblks_per_bin,kind=real_8))/REAL(SIZE(nblks_per_bin),kind=real_8)
          IF (parai%me == parai%source.AND.&
               parai%cp_inter_me == 0)WRITE(6,'(2A,2(I0,1X),F0.2)')procedureN,&
               ' costs min/max/avg : ',&
               MINVAL(heap%vals),MAXVAL(heap%vals),&
               SUM(REAL(heap%vals,kind=real_8))/REAL(SIZE(heap%vals),kind=real_8)
       ENDIF

       dist%nblks_loc = nblks_per_bin( parai%cp_inter_me + 1 )
       ALLOCATE( dist%dist( dist%nblks_loc ), stat=ierr )
       IF (ierr /= 0)CALL stopgm(procedureN,'Allocation problem',& 
            __LINE__,__FILE__)

       ! 
       ! Get the local number of block ids
       ! 
       CALL heap_fill( heap, (/( grp, grp = 1,parai%cp_nogrp )/),&
            (/( 0, grp = 1,parai%cp_nogrp )/) )

       IF (parai%me == parai%source) THEN
          i = 1
          DO blk = 1,nblks
             CALL heap_get_first(heap,root_bin,root_cost)
             IF ( parai%cp_inter_me + 1  ==  root_bin) THEN
                dist%dist( i ) = blk
                i = i + 1
             ENDIF
             root_cost = root_cost + dist_costs(blk)
             CALL heap_reset_first(heap,root_cost)
          ENDDO
       ENDIF

       CALL mp_bcast(dist%dist,SIZE(dist%dist),parai%source,parai%allgrp)

       CALL heap_release( heap )

       DEALLOCATE( nblks_per_bin, stat=ierr )
       IF (ierr /= 0)CALL stopgm(procedureN,'Deallocation problem',& 
            __LINE__,__FILE__)

    CASE default

       CALL stopgm(procedureN,'not a valide distribution',& 
            __LINE__,__FILE__)

    END SELECT

    ! ==--------------------------------------------------------------==
  END SUBROUTINE build_dist
  ! ================================================================== 
  SUBROUTINE release_dist(dist)
    ! ==--------------------------------------------------------------==
    TYPE(dist_t)                             :: dist

    CHARACTER(*), PARAMETER                  :: procedureN = 'release_dist'

    INTEGER                                  :: ierr

    IF (ASSOCIATED(dist%dist)) THEN
       DEALLOCATE(dist%dist,stat=ierr)
       IF (ierr /= 0)CALL stopgm(procedureN,'Deallocation problem',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
  END SUBROUTINE release_dist
  ! ================================================================== 
  SUBROUTINE iterator_start(iter,dist)
    ! ==--------------------------------------------------------------==
    TYPE(iterator_t)                         :: iter
    TYPE(dist_t)                             :: dist

    iter%blk_count = 0
    iter%max_blk_count = dist%nblks_loc
    iter%dist => dist%dist
    ! ==--------------------------------------------------------------==
  END SUBROUTINE iterator_start
  ! ================================================================== 
  LOGICAL FUNCTION iterator_blocks(iter)
    ! ==--------------------------------------------------------------==
    TYPE(iterator_t) :: iter
    iter%blk_count = iter%blk_count + 1
    iterator_blocks = iter%blk_count <= iter%max_blk_count
    ! ==--------------------------------------------------------------==
  END FUNCTION iterator_blocks
  ! ================================================================== 
  SUBROUTINE iterator_stop(iter)
    ! ==--------------------------------------------------------------==
    TYPE(iterator_t)                         :: iter

    iter%blk_count = 0
    iter%max_blk_count = 0
    NULLIFY(iter%dist)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE iterator_stop
  ! ================================================================== 
  SUBROUTINE iterator_next_block(iter,blk)
    ! ==--------------------------------------------------------------==
    TYPE(iterator_t)                         :: iter
    INTEGER, INTENT(out)                     :: blk

    blk = iter%dist(iter%blk_count)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE iterator_next_block
  ! ================================================================== 
  SUBROUTINE build_block_info(n,nbr_blk,blk_size,&
       offsets,sizes)
    ! ==--------------------------------------------------------------==
    INTEGER, INTENT(in)                      :: n, nbr_blk, blk_size
    INTEGER, DIMENSION(:), INTENT(out)       :: offsets, sizes

    CHARACTER(*), PARAMETER :: procedureN = 'build_block_info'

    INTEGER                                  :: blk

    IF (SIZE(offsets) < nbr_blk) CALL stopgm(procedureN,&
         'size of array not valid',& 
         __LINE__,__FILE__)
    IF (SIZE(sizes) < nbr_blk) CALL stopgm(procedureN,&
         'size of array not valid',& 
         __LINE__,__FILE__)

    offsets(1) = 1
    DO blk = 1,nbr_blk-1
       offsets(blk+1) = offsets(blk) + blk_size
       sizes(blk) = blk_size
    ENDDO
    ! need to protect vs nbr_blk=0
    IF (nbr_blk > 0) sizes(nbr_blk) = n - offsets(nbr_blk) + 1

    ! ==--------------------------------------------------------------==
  END SUBROUTINE build_block_info
  ! ================================================================== 
  SUBROUTINE pp_init( pp, max_entry )
    ! ==--------------------------------------------------------------==
    TYPE(pp_t), INTENT(inout)                :: pp
    INTEGER, INTENT(in)                      :: max_entry

    CHARACTER(*), PARAMETER                  :: procedureN = 'pp_init'

    INTEGER                                  :: ierr

    pp%n_pp=0
    ALLOCATE(pp%pp_ptr(max_entry**2+1),pp%pp(2,max_entry**2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'Allocation problem',&
         & __LINE__,__FILE__)
    pp%pp_ptr(:)=HUGE(0);pp%pp(:,:)=HUGE(0)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE pp_init
  ! ================================================================== 
  SUBROUTINE pp_destroy( pp )
    ! ==--------------------------------------------------------------==
    TYPE(pp_t), INTENT(inout)                :: pp

    CHARACTER(*), PARAMETER                  :: procedureN = 'pp_destroy'

    INTEGER                                  :: ierr

    pp%n_pp=0
    DEALLOCATE(pp%pp_ptr,pp%pp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'Allocation problem',&
         & __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE pp_destroy
  ! ================================================================== 
  SUBROUTINE list_init( list, max_list_size )
    ! ==--------------------------------------------------------------==
    TYPE(list_t), INTENT(inout)              :: list
    INTEGER, INTENT(in)                      :: max_list_size

    CHARACTER(*), PARAMETER                  :: procedureN = 'list_init'

    INTEGER                                  :: ierr

    list%n=0
    ALLOCATE(list%a(max_list_size),list%r(max_list_size),&
         & list%stored_in(max_list_size),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'Allocation problem',&
         & __LINE__,__FILE__)
    list%a(:) = HUGE(0)
    list%r(:) = HUGE(0)
    list%stored_in(:) = HUGE(0)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE list_init
  ! ================================================================== 
  SUBROUTINE list_destroy( list )
    ! ==--------------------------------------------------------------==
    TYPE(list_t), INTENT(inout)              :: list

    CHARACTER(*), PARAMETER                  :: procedureN = 'list_destroy'

    INTEGER                                  :: ierr

    list%n=0
    DEALLOCATE(list%a,list%r,&
         & list%stored_in,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'Deallocation problem',&
         & __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE list_destroy
  ! ================================================================== 
  INTEGER FUNCTION get_ticks()
    ! ==--------------------------------------------------------------==
    INTEGER :: count
    INTEGER,SAVE :: count_max,count_rate,cycles=-666,last_count
    IF (cycles == -666) THEN
       CALL SYSTEM_CLOCK(count,count_rate,count_max)
       cycles=0;last_count=0
    ENDIF
    CALL SYSTEM_CLOCK(count)
    IF (count < last_count) THEN
       IF (last_count-count < count_max/100) THEN
          count=last_count
       ELSE
          cycles=cycles+1
       ENDIF
    ENDIF
    last_count=count
    get_ticks = count + cycles * count_max
    ! ==--------------------------------------------------------------==
  END FUNCTION get_ticks
  ! ================================================================== 
  SUBROUTINE block_filter(thresh,row_offset,row_size,col_offset,col_size,pp,r_list,c_list)
    ! ==--------------------------------------------------------------==

    TYPE(thresh_t), INTENT(inout)            :: thresh
    INTEGER, INTENT(in)                      :: row_offset, row_size, &
                                                col_offset, col_size
    TYPE(pp_t), INTENT(inout)                :: pp
    TYPE(list_t), INTENT(inout)              :: r_list, c_list

    CHARACTER(*), PARAMETER                  :: procedureN = 'block_filter'
    LOGICAL, PARAMETER                       :: FORCE_DENSE = .FALSE.

    INTEGER                                  :: c, col_abs, i_batch, ij, &
                                                isub, n_batch, p, r, row_abs, &
                                                row_batch(row_size)
    LOGICAL                                  :: col_is_present(col_size), &
                                                row_is_present(row_size)
    REAL(real_8)                             :: dist_ab

    CALL tiset(procedureN,isub)

    IF(FORCE_DENSE) THEN
       !dense case

       r_list%n = row_size
       DO r = 1, r_list%n
          r_list%r(r) = r
          r_list%a(r) = row_offset + r - 1
       ENDDO

       c_list%n = col_size
       DO c = 1, c_list%n
          c_list%r(c) = c
          c_list%a(c) = col_offset + c - 1
       ENDDO


       pp%pp(:,:) = HUGE(0)
       pp%pp_ptr(:) = HUGE(0)

       pp%r_offset = row_offset
       pp%c_offset = col_offset
       pp%n_pp = 0
       p = 1
       DO c = 1, c_list%n
          DO r = 1, r_list%n, 2
             pp%n_pp = pp%n_pp + 1
             pp%pp_ptr(pp%n_pp) = p
             IF( r == r_list%n ) THEN
                pp%pp(1,p  ) = r
                pp%pp(2,p  ) = c
                p = p + 1
             ELSE
                pp%pp(1,p  ) = r
                pp%pp(1,p+1) = r + 1
                pp%pp(2,p  ) = c
                pp%pp(2,p+1) = c
                p = p + 2
             ENDIF
          ENDDO
       ENDDO
       pp%pp_ptr(pp%n_pp+1) = p

    ELSE
       !sparse case

       row_is_present(:) = .FALSE.
       col_is_present(:) = .FALSE.

       pp%pp(:,:) = HUGE(0)
       pp%pp_ptr(:) = HUGE(0)

       pp%r_offset = row_offset
       pp%c_offset = col_offset

       pp%n_pp = 0
       p = 1
       col_loop: DO c = 1, col_size
          col_abs = col_offset + c - 1

          n_batch = 0
          row_loop: DO r = 1, row_size
             row_abs = row_offset + r - 1

             !
             ! screen if needed
             !
             IF(hfxc3%twscr) THEN
                !
                ! First: check that the separation distance is smaller that the max distance allowed
                !
                IF(hfxc3%twfc) THEN
                   CALL get_wannier_separation(row_abs,col_abs,dist_ab)
                   IF(dist_ab > hfxc4%dwfc) THEN
                      CALL sub2ind(row_abs,col_abs,ij)
                      ! Set the int_vals to something that if the
                      ! pair comes within the radius any time later
                      ! we properly recompute the integral.
                      thresh%int_vals(ij) = 1.0_real_8
                      !nbr_int_skiped_1 = nbr_int_skiped_1 + 1
                      CYCLE row_loop
                   ENDIF
                ENDIF
                !
                ! Second: check that the integral is larger than the integral threshold
                !
                IF(.NOT.thresh%init_ints) THEN
                   CALL sub2ind(row_abs,col_abs,ij)

                   !if(paral%parent)write(6,'(A,3I4,2E12.3,L6)') 'row_abs,col_abs,ij,int,thresh',row_abs,col_abs,ij,&
                   !     thresh%int_vals(ij),hfxc4%dwf_integral_thresh,thresh%int_vals(ij) < hfxc4%dwf_integral_thresh

                   IF (thresh%int_vals(ij) < hfxc4%dwf_integral_thresh) THEN
                      !nbr_int_skiped_1 = nbr_int_skiped_1 + 1
                      CYCLE row_loop
                   ENDIF
                ENDIF
             ENDIF

             !
             ! mark the pair
             !
             n_batch = n_batch + 1
             row_batch(n_batch) = r
             row_is_present(r) = .TRUE.
             col_is_present(c) = .TRUE.

          ENDDO row_loop


          batch_loop: DO i_batch=1,n_batch,2

             pp%n_pp = pp%n_pp + 1
             pp%pp_ptr(pp%n_pp) = p
             IF( i_batch == n_batch ) THEN
                pp%pp(1,p  ) = row_batch(i_batch)
                pp%pp(2,p  ) = c
                p = p + 1
             ELSE
                pp%pp(1,p  ) = row_batch(i_batch)
                pp%pp(1,p+1) = row_batch(i_batch+1)
                pp%pp(2,p  ) = c
                pp%pp(2,p+1) = c
                p = p + 2
             ENDIF

          ENDDO batch_loop

       ENDDO col_loop
       pp%pp_ptr(pp%n_pp+1) = p

       !
       ! unique rows
       !
       r_list%n = 0
       DO r = 1, row_size
          IF(row_is_present(r)) THEN
             r_list%n = r_list%n + 1
             r_list%r(r_list%n) = r
             r_list%a(r_list%n) = row_offset + r - 1
          ENDIF
       ENDDO

       !
       ! unique cols
       !
       c_list%n = 0
       DO c = 1, col_size
          IF(col_is_present(c)) THEN
             c_list%n = c_list%n + 1
             c_list%r(c_list%n) = c
             c_list%a(c_list%n) = col_offset + c - 1
          ENDIF
       ENDDO

    ENDIF


!!$    write(6,*) 'r_list%n',r_list%n
!!$    write(6,*) 'r_list%a',r_list%a(1:r_list%n)
!!$    write(6,*) 'r_list%r',r_list%r(1:r_list%n)
!!$
!!$    write(6,*) 'c_list%n',c_list%n
!!$    write(6,*) 'c_list%a',c_list%a(1:c_list%n)
!!$    write(6,*) 'c_list%r',c_list%r(1:c_list%n)
!!$
!!$    write(6,*) 'pp%n_pp',pp%n_pp
!!$    write(6,*) 'pp%pp_ptr',pp%pp_ptr
!!$    write(6,*) 'pp%pp(1,:)',pp%pp(1,:)
!!$    write(6,*) 'pp%pp(2,:)',pp%pp(2,:)

    CALL tihalt(procedureN,isub)

    ! ==--------------------------------------------------------------==
  END SUBROUTINE block_filter
  ! ================================================================== 

  ! ================================================================== 
  ! Block FFT 
  ! ================================================================== 

  ! ================================================================== 
  SUBROUTINE block_invfft_old(c0,psis,offset,size)
    ! ==--------------------------------------------------------------==
    ! KPT replace NGW -> NGWK
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8)                          :: c0(:,:), psis(:,:)
    INTEGER, INTENT(in)                      :: offset, size

    CHARACTER(*), PARAMETER                  :: procedureN = 'block_invfft_old'
    COMPLEX(real_8), PARAMETER               :: zone = (1.0_real_8,0.0_real_8)

    INTEGER                                  :: i, isub, nbr_ffts, ptr

    ! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)

    ! ==--------------------------------------------------------------==

    ptr = 1
    nbr_ffts = 0
    ! KPT replace 2 -> jmup (=1 for kpt)
    DO i = offset,offset+size-1,2
       CALL zeroing(psis(:,ptr))
       ! KPT          IF(TKPNT) THEN
       ! KPT            CALL SET_PSI_1_STATE_G_KPTS(ZONE,C0(1,IS1),PSI)
       ! KPT          ELSE
       IF (i == offset+size-1) THEN
          CALL set_psi_1_state_g(zone,c0(:,i),psis(:,ptr))
       ELSE
          CALL set_psi_2_states_g(c0(:,i),c0(:,i+1),&
               psis(:,ptr))
       ENDIF
       nbr_rwfn_precomputed = nbr_rwfn_precomputed + 1
       nbr_ffts = nbr_ffts + 1
       ptr = ptr + 1
    ENDDO

    ! fft all
    DO ptr=1,nbr_ffts
       CALL invfftn(psis(:,ptr),.TRUE.,parai%allgrp)
       IF (cntl%use_scaled_hfx) CALL scex%do_density_scaling(psis(:,ptr))
    ENDDO

    ! ==--------------------------------------------------------------==

    CALL tihalt(procedureN,isub)

    ! ==--------------------------------------------------------------==
  END SUBROUTINE block_invfft_old
  ! ================================================================== 
  SUBROUTINE block_invfft(c0,psis,list)
    ! ==--------------------------------------------------------------==
    ! KPT replace NGW -> NGWK
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8)                          :: c0(:,:), psis(:,:)
    TYPE(list_t), INTENT(inout)              :: list

    CHARACTER(*), PARAMETER                  :: procedureN = 'block_invfft'
    COMPLEX(real_8), PARAMETER               :: zone = (1.0_real_8,0.0_real_8)

    INTEGER                                  :: i0_a, i0_r, i1_a, i1_r, isub, &
                                                l, nbr_ffts, ptr


    CALL tiset(procedureN,isub)

    ! ==--------------------------------------------------------------==

    ptr = 1
    nbr_ffts = 0
    ! KPT replace 2 -> jmup (=1 for kpt)
    DO l = 1,list%n,2
       CALL zeroing(psis(:,ptr))
       ! KPT          IF(TKPNT) THEN
       ! KPT            CALL SET_PSI_1_STATE_G_KPTS(ZONE,C0(1,IS1),PSI)
       ! KPT          ELSE
       IF (l == list%n) THEN
          i0_a = list%a(l)
          i0_r = list%r(l)
          list%stored_in(i0_r) = ptr ! real
          CALL set_psi_1_state_g(zone,c0(:,i0_a),psis(:,ptr))
       ELSE
          i0_a = list%a(l)
          i1_a = list%a(l+1)
          i0_r = list%r(l)
          i1_r = list%r(l+1)
          list%stored_in(i0_r) =  ptr ! real
          list%stored_in(i1_r) = -ptr ! imag
          CALL set_psi_2_states_g(c0(:,i0_a),c0(:,i1_a),&
               psis(:,ptr))
       ENDIF
       nbr_rwfn_precomputed = nbr_rwfn_precomputed + 1
       nbr_ffts = nbr_ffts + 1
       ptr = ptr + 1
    ENDDO

    ! fft all
    DO ptr=1,nbr_ffts
       CALL invfftn(psis(:,ptr),.TRUE.,parai%allgrp)
       IF (cntl%use_scaled_hfx) CALL scex%do_density_scaling(psis(:,ptr))
    ENDDO

    ! ==--------------------------------------------------------------==

    CALL tihalt(procedureN,isub)

    ! ==--------------------------------------------------------------==
  END SUBROUTINE block_invfft
  ! ================================================================== 



  ! ================================================================== 
  SUBROUTINE hfx_compute_block_new(C2_hfx,ehfx,f,vpotg,vpotr,psic,pp,&
       psi_row,psi_col,r_list,c_list,psi_row_pack,thresh)
    ! ==--------------------------------------------------------------==
    ! KPT replace NGW -> NGWK
    COMPLEX(real_8) :: C2_hfx(:,:)
    REAL(real_8), INTENT(inout) :: ehfx
    ! KPT F may need to be in F(N,NKPT)
    REAL(real_8) :: f(:),vpotr(:,:)
    COMPLEX(real_8) :: vpotg(:,:),psic(:)
    COMPLEX(real_8) :: psi_col(:,:),psi_row(:,:)
    TYPE(pp_t), INTENT(in) :: pp
    TYPE(list_t), INTENT(in) :: r_list, c_list
    COMPLEX(real_8), DIMENSION(:), INTENT(inout) :: psi_row_pack
    TYPE(thresh_t), INTENT(inout) :: thresh
    ! ==--------------------------------------------------------------==

    INTEGER :: isub,r1,r2,c,i0,i1,n,p,ierr,r1_r,r2_r,c_r
    INTEGER :: c_pt,r1_pt,r2_pt
    REAL(real_8) :: pfl,pfx1,pfx2,EHFX_loc_1,EHFX_loc_2
    REAL(real_8) :: EHFX_nothresh_1,EHFX_nothresh_2
    CHARACTER(*),PARAMETER :: procedureN='hfx_compute_block'
    LOGICAL :: c_in_re, r1_in_re, r2_in_re

    ! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)

    ! ==--------------------------------------------------------------==

    IF (cntl%tlsd) THEN
       pfl = 0.50_real_8
    ELSE
       pfl = 0.25_real_8
    ENDIF

    IF (cntl%thybrid) pfl = pfl * func3%phfx

    ! ==--------------------------------------------------------------==

    psi_row_pack(:)=HUGE(0.0_real_8)         ! needed?

    DO p = 1, pp%n_pp
       n = pp%pp_ptr(p+1) - pp%pp_ptr(p)
       r1 = 0; r2 = 0
       IF( n == 1 ) THEN
          i0 = pp%pp_ptr(p)
          !
          r1_r = pp%pp(1,i0)
          r1 = r1_r + pp%r_offset - 1
          r1_pt = ABS( r_list%stored_in(r1_r) )
          r1_in_re = r_list%stored_in(r1_r) > 0
          !
          c_r = pp%pp(2,i0)
          c = c_r + pp%c_offset - 1
          c_pt = ABS( c_list%stored_in(c_r) )
          c_in_re = c_list%stored_in(c_r) > 0
          !
          ! copy row(r1) to tmp
          IF(r1_in_re) THEN
             CALL copy_re_to_re(llr1,psi_row(:,r1_pt),psi_row_pack)
          ELSE
             CALL copy_im_to_re(llr1,psi_row(:,r1_pt),psi_row_pack)
          ENDIF
       ELSEIF( n == 2 ) THEN
          i0 = pp%pp_ptr(p)
          i1 = pp%pp_ptr(p)+1
          !
          r1_r = pp%pp(1,i0)
          r2_r = pp%pp(1,i1)
          r1 = r1_r + pp%r_offset - 1
          r2 = r2_r + pp%r_offset - 1
          r1_pt = ABS( r_list%stored_in(r1_r) )
          r2_pt = ABS( r_list%stored_in(r2_r) )
          r1_in_re = r_list%stored_in(r1_r) > 0
          r2_in_re = r_list%stored_in(r2_r) > 0
          !
          c_r = pp%pp(2,i0)
          c = c_r + pp%c_offset - 1
          c_pt = ABS( c_list%stored_in(c_r) )
          c_in_re = c_list%stored_in(c_r) > 0
          !
          ! copy row(r1) to tmp
          IF(r1_in_re) THEN
             CALL copy_re_to_re(llr1,psi_row(:,r1_pt),psi_row_pack)
          ELSE
             CALL copy_im_to_re(llr1,psi_row(:,r1_pt),psi_row_pack)
          ENDIF
          !
          ! copy row(r2) to buff
          IF(r2_in_re) THEN
             CALL copy_re_to_im(llr1,psi_row(:,r2_pt),psi_row_pack)
          ELSE
             CALL copy_im_to_im(llr1,psi_row(:,r2_pt),psi_row_pack)
          ENDIF
       ELSE
          CALL stopgm(procedureN,'wrong logic',&
               __LINE__,__FILE__)
       ENDIF

       ! need to:
       !  DONE 1) copy r1,r2 to some buff
       !  DONE 2) change the var name r_ptr, c_stored_in_real and add the row buff

       IF (r1 > c) r1 = 0
       IF (r2 > c) r2 = 0
       IF (r1 /= 0.OR.r2 /= 0) CALL evaluate_pair()

    ENDDO

    ! ==--------------------------------------------------------------==

    CALL tihalt(procedureN,isub)

    ! ==--------------------------------------------------------------==

  CONTAINS

    ! ================================================================== 
    SUBROUTINE evaluate_pair()
      ! ==--------------------------------------------------------------==
      CHARACTER(*), PARAMETER                  :: procedureN = 'evaluate_pair'

      EHFX_loc_1 = 0.0_real_8
      EHFX_loc_2 = 0.0_real_8
      EHFX_nothresh_1 = 0.0_real_8
      EHFX_nothresh_2 = 0.0_real_8

      IF (r1 /= 0.AND.r2 /= 0) THEN
         nbr_integrals = nbr_integrals+2
         pfx1 = pfl * f(c) * f(r1)
         pfx2 = pfl * f(c) * f(r2)
         ! This is needed to avoid aliasing
         IF (c == r1) THEN
            CALL hfxaa_new(EHFX_loc_1,EHFX_nothresh_1,pfx1,psi_col(:,c_pt),c_in_re,vpotg(:,1),vpotr(:,1),psic,C2_hfx(:,c))
            IF(hfxc3%twscr) CALL add_int_to_list(EHFX_nothresh_1,c,c,thresh%acc_vals,thresh%new_vals)
            ehfx = ehfx + EHFX_loc_1
            CALL hfxab_new(EHFX_loc_1,EHFX_nothresh_1,pfx2,psi_col(:,c_pt),c_in_re,psi_row_pack,.FALSE.,&
                           vpotg(:,1),vpotr(:,1),psic,C2_hfx(:,c),C2_hfx(:,r2))
            IF(hfxc3%twscr) CALL add_int_to_list(EHFX_nothresh_1,c,r2,thresh%acc_vals,thresh%new_vals)
            ehfx = ehfx + EHFX_loc_1
         ELSEIF (c == r2) THEN
            CALL hfxaa_new(EHFX_loc_1,EHFX_nothresh_1,pfx2,psi_col(:,c_pt),c_in_re,vpotg(:,1),vpotr(:,1),psic,C2_hfx(:,c))
            IF(hfxc3%twscr) CALL add_int_to_list(EHFX_nothresh_1,c,c,thresh%acc_vals,thresh%new_vals)
            ehfx = ehfx + EHFX_loc_1
            CALL hfxab_new(EHFX_loc_1,EHFX_nothresh_1,pfx1,psi_col(:,c_pt),c_in_re,psi_row_pack,.TRUE.,&
                           vpotg(:,1),vpotr(:,1),psic,C2_hfx(:,c),C2_hfx(:,r1))
            IF(hfxc3%twscr) CALL add_int_to_list(EHFX_nothresh_1,c,r1,thresh%acc_vals,thresh%new_vals)
            ehfx = ehfx + EHFX_loc_1

         ELSE
            CALL hfxab2_new(EHFX_loc_1,EHFX_loc_2,EHFX_nothresh_1,EHFX_nothresh_2,&
                 pfx1,pfx2,psi_col(:,c_pt),c_in_re,psi_row_pack,vpotg,vpotr,&
                 psic,C2_hfx(:,c),C2_hfx(:,r1),C2_hfx(:,r2))
            IF(hfxc3%twscr) CALL add_int_to_list(EHFX_nothresh_1,c,r1,thresh%acc_vals,thresh%new_vals)
            IF(hfxc3%twscr) CALL add_int_to_list(EHFX_nothresh_2,c,r2,thresh%acc_vals,thresh%new_vals)
            ehfx = ehfx + EHFX_loc_1 + EHFX_loc_2
         ENDIF
      ELSEIF (r1 /= 0.AND.r2 == 0) THEN
         nbr_integrals = nbr_integrals+1
         pfx1 = pfl * f(c) * f(r1)
         IF (c == r1) THEN
            CALL hfxaa_new(EHFX_loc_1,EHFX_nothresh_1,pfx1,psi_col(:,c_pt),c_in_re,vpotg(:,1),vpotr(:,1),psic,C2_hfx(:,c))
            IF(hfxc3%twscr) CALL add_int_to_list(EHFX_nothresh_1,c,c,thresh%acc_vals,thresh%new_vals)
            ehfx = ehfx + EHFX_loc_1
         ELSE
            CALL hfxab_new(EHFX_loc_1,EHFX_nothresh_1,pfx1,psi_col(:,c_pt),c_in_re,psi_row_pack,.TRUE.,&
                 vpotg(:,1),vpotr(:,1),psic,C2_hfx(:,c),C2_hfx(:,r1))
            IF(hfxc3%twscr) CALL add_int_to_list(EHFX_nothresh_1,c,r1,thresh%acc_vals,thresh%new_vals)
            ehfx = ehfx + EHFX_loc_1
         ENDIF
      ELSEIF (r1 == 0.AND.r2 /= 0) THEN
         nbr_integrals = nbr_integrals+1
         pfx2 = pfl * f(c) * f(r2)
         IF (c == r2) THEN
            CALL hfxaa_new(EHFX_loc_1,EHFX_nothresh_1,pfx1,psi_col(:,c_pt),c_in_re,vpotg(:,1),vpotr(:,1),psic,C2_hfx(:,c))
            IF(hfxc3%twscr) CALL add_int_to_list(EHFX_nothresh_1,c,c,thresh%acc_vals,thresh%new_vals)
            ehfx = ehfx + EHFX_loc_1
         ELSE
            CALL hfxab_new(EHFX_loc_1,EHFX_nothresh_1,pfx2,psi_col(:,c_pt),c_in_re,psi_row_pack,.FALSE.,&
                 vpotg(:,1),vpotr(:,1),psic,C2_hfx(:,c),C2_hfx(:,r2))
            IF(hfxc3%twscr) CALL add_int_to_list(EHFX_nothresh_1,c,r2,thresh%acc_vals,thresh%new_vals)
            ehfx = ehfx + EHFX_loc_1
         ENDIF
      ELSE
         CALL stopgm(procedureN,'wrong logic',&
              __LINE__,__FILE__)
      ENDIF


      ! ==--------------------------------------------------------------==
    END SUBROUTINE evaluate_pair
    ! ================================================================== 
    SUBROUTINE add_int_to_list(v,i,j,avs,nvs)
      ! ==--------------------------------------------------------------==
      REAL(real_8), INTENT(in)                 :: v
      INTEGER, INTENT(in)                      :: i, j
      REAL(real_8), DIMENSION(:), &
        INTENT(inout)                          :: avs
      INTEGER, DIMENSION(:), INTENT(inout)     :: nvs

      INTEGER                                  :: ij

        CALL sub2ind(i,j,ij)
        avs( ij ) = v
        nvs( ij ) = 1
      ! ==--------------------------------------------------------------==
    END SUBROUTINE add_int_to_list
    ! ==--------------------------------------------------------------==
  END SUBROUTINE hfx_compute_block_new
  ! ================================================================== 
  SUBROUTINE hfx_compute_block(C2_hfx,ehfx,f,vpotg,vpotr,psic,&
       psi_row,row_offset,row_size,psi_col,col_offset,col_size)
    ! ==--------------------------------------------------------------==
    ! KPT replace NGW -> NGWK
    COMPLEX(real_8) :: C2_hfx(:,:)
    REAL(real_8), INTENT(inout) :: ehfx
    ! KPT F may need to be in F(N,NKPT)
    REAL(real_8) :: f(:),vpotr(:,:)
    COMPLEX(real_8) :: vpotg(:,:),psic(:)
    COMPLEX(real_8) :: psi_col(:,:),psi_row(:,:)
    INTEGER, INTENT(in) :: row_offset, row_size, col_offset, col_size

    ! ==--------------------------------------------------------------==

    INTEGER :: isub,r,r1,r2,c,row_beg,row_end,col_beg,col_end
    INTEGER :: c_pt,r_ptr
    REAL(real_8) :: pfl,pfx1,pfx2,EHFX_loc,EHFX_loc_1,EHFX_loc_2
    CHARACTER(*),PARAMETER :: procedureN='hfx_compute_block'
    LOGICAL :: c_stored_in_real

    ! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)

    ! ==--------------------------------------------------------------==

    IF (cntl%tlsd) THEN
       pfl = 0.50_real_8
    ELSE
       pfl = 0.25_real_8
    ENDIF

    IF (cntl%thybrid) pfl = pfl * func3%phfx

    ! ==--------------------------------------------------------------==

    col_beg = col_offset; col_end = col_offset + col_size - 1
    row_beg = row_offset; row_end = row_offset + row_size - 1
    c_pt = 1
    c_stored_in_real = .TRUE.
    DO c = col_beg, col_end
       r_ptr = 1
       DO r = row_beg, row_end, 2
          r1 = r
          r2 = r + 1
          IF (r1 > c) r1 = 0
          IF (r2 > c.OR.r2 > row_end) r2 = 0
          IF (r1 /= 0.OR.r2 /= 0) CALL evaluate_pair()
          r_ptr = r_ptr + 1
       ENDDO
       c_pt = c_pt + MOD(c - col_beg, 2)
       c_stored_in_real = .NOT.c_stored_in_real
    ENDDO

    ! ==--------------------------------------------------------------==

    CALL tihalt(procedureN,isub)

    ! ==--------------------------------------------------------------==

  CONTAINS

    ! ================================================================== 
    SUBROUTINE evaluate_pair()
      ! ==--------------------------------------------------------------==
    CHARACTER(*), PARAMETER                  :: procedureN = 'evaluate_pair'

      EHFX_loc = 0.0_real_8
      IF (r1 /= 0.AND.r2 /= 0) THEN
         nbr_integrals = nbr_integrals+2
         pfx1 = pfl * f(c) * f(r1)
         pfx2 = pfl * f(c) * f(r2)
         ! This is needed to avoid aliasing
         IF (c == r1) THEN
            CALL hfxaa(EHFX_loc,pfx1,&
                 psi_col(:,c_pt),c_stored_in_real,&
                 vpotg(:,1),vpotr(:,1),psic,C2_hfx(:,c))
            ehfx = ehfx + EHFX_loc
            CALL hfxab(EHFX_loc,pfx2,&
                 psi_col(:,c_pt),c_stored_in_real,&
                 psi_row(:,r_ptr),.FALSE.,&
                 vpotg(:,1),vpotr(:,1),psic,C2_hfx(:,c),c2_hfx(:,r2))
            ehfx = ehfx + EHFX_loc
         ELSEIF (c == r2) THEN
            CALL hfxaa(EHFX_loc,pfx2,&
                 psi_col(:,c_pt),c_stored_in_real,&
                 vpotg(:,1),vpotr(:,1),psic,C2_hfx(:,c))
            ehfx = ehfx + EHFX_loc
            CALL hfxab(EHFX_loc,pfx1,&
                 psi_col(:,c_pt),c_stored_in_real,&
                 psi_row(:,r_ptr),.TRUE.,&
                 vpotg(:,1),vpotr(:,1),psic,C2_hfx(:,c),c2_hfx(:,r1))
            ehfx = ehfx + EHFX_loc
         ELSE
            CALL hfxab2(EHFX_loc_1,EHFX_loc_2,pfx1,pfx2,&
                 psi_col(:,c_pt),c_stored_in_real,&
                 psi_row(:,r_ptr),&
                 vpotg,vpotr,psic,&
                 C2_hfx(:,c),c2_hfx(:,r1),c2_hfx(:,r2))
            ehfx = ehfx + EHFX_loc_1 + EHFX_loc_2
         ENDIF
      ELSEIF (r1 /= 0.AND.r2 == 0) THEN
         nbr_integrals = nbr_integrals+1
         pfx1 = pfl * f(c) * f(r1)
         IF (c == r1) THEN
            CALL hfxaa(EHFX_loc,pfx1,&
                 psi_col(:,c_pt),c_stored_in_real,&
                 vpotg(:,1),vpotr(:,1),psic,C2_hfx(:,c))
         ELSE
            CALL hfxab(EHFX_loc,pfx1,&
                 psi_col(:,c_pt),c_stored_in_real,&
                 psi_row(:,r_ptr),.TRUE.,&
                 vpotg(:,1),vpotr(:,1),psic,C2_hfx(:,c),c2_hfx(:,r1))
         ENDIF
         ehfx = ehfx + EHFX_loc
      ELSEIF (r1 == 0.AND.r2 /= 0) THEN
         nbr_integrals = nbr_integrals+1
         pfx2 = pfl * f(c) * f(r2)
         IF (c == r2) THEN
            CALL hfxaa(EHFX_loc,pfx1,&
                 psi_col(:,c_pt),c_stored_in_real,&
                 vpotg(:,1),vpotr(:,1),psic,C2_hfx(:,c))
         ELSE
            CALL hfxab(EHFX_loc,pfx2,&
                 psi_col(:,c_pt),c_stored_in_real,&
                 psi_row(:,r_ptr),.FALSE.,&
                 vpotg(:,1),vpotr(:,1),psic,C2_hfx(:,c),c2_hfx(:,r2))
         ENDIF
         ehfx = ehfx + EHFX_loc
      ELSE
         CALL stopgm(procedureN,'wrong logic',& 
              __LINE__,__FILE__)
      ENDIF

      ! ==--------------------------------------------------------------==
    END SUBROUTINE evaluate_pair
    ! ================================================================== 

    ! ==--------------------------------------------------------------==
  END SUBROUTINE hfx_compute_block
  ! ================================================================== 

  ! ==================================================================
  ! Coulomb potential between orbital pairs: Old routines
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE hfxaa(ehfx,pf,psia,a_stored_in_real,&
       vpotg,vpotr,psic,c2a)
    ! ==================================================================
    REAL(real_8)                             :: ehfx, pf
    COMPLEX(real_8)                          :: psia(:)
    LOGICAL                                  :: a_stored_in_real
    COMPLEX(real_8)                          :: vpotg(:)
    REAL(real_8)                             :: vpotr(:)
    COMPLEX(real_8)                          :: psic(:), c2a(:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'hfxaa'

    COMPLEX(real_8)                          :: fm, fp
    INTEGER                                  :: ig, ir, isub

    CALL tiset(procedureN,isub)
    ehfx=0.0_real_8
    IF (a_stored_in_real) THEN
       !$omp parallel do private(IR)
       DO ir=1,llr1
          psic(ir)=REAL(psia(ir))*REAL(psia(ir))
       ENDDO
    ELSE
       !$omp parallel do private(IR)
       DO ir=1,llr1
          psic(ir)=AIMAG(psia(ir))*AIMAG(psia(ir))
       ENDDO
    ENDIF
    CALL dscal(2*llr1,1._real_8/parm%omega,psic,1)
    CALL fwfftn(psic,.FALSE.,parai%allgrp)
    !$omp parallel do private(IG,FP) &
    !$omp  reduction(+:ehfx)
    DO ig=1,jhg
       fp=psic(nzff(ig))
       vpotg(ig)=-pf*scgx(ig)*fp
       ehfx=ehfx+REAL(2._real_8*vpotg(ig)*CONJG(fp))
    ENDDO
    IF (geq0) ehfx=ehfx-REAL(vpotg(1)*CONJG(psic(nzff(1))))
    CALL zeroing(psic)
    !ocl novrec
    !$omp parallel do private(IG)
    DO ig=1,jhg
       psic(nzff(ig))=vpotg(ig)
       psic(inzf(ig))=CONJG(vpotg(ig))
    ENDDO
    IF (geq0) psic(nzff(1))=vpotg(1)
    CALL invfftn(psic,.FALSE.,parai%allgrp)
    !$omp parallel do private(IR)
    DO ir=1,llr1
       vpotr(ir)=REAL(psic(ir))
    ENDDO
    IF (a_stored_in_real) THEN
       !$omp parallel do private(IR)
       DO ir=1,llr1
          psic(ir)=vpotr(ir)*REAL(psia(ir))
       ENDDO
    ELSE
       !$omp parallel do private(IR)
       DO ir=1,llr1
          psic(ir)=vpotr(ir)*AIMAG(psia(ir))
       ENDDO
    ENDIF
    CALL fwfftn(psic,.TRUE.,parai%allgrp)
    !$omp parallel do private(IG,FP,FM)
    DO ig=1,jgw
       fp=psic(nzfs(ig))+psic(inzs(ig))
       fm=psic(nzfs(ig))-psic(inzs(ig))
       c2a(ig)=c2a(ig)-CMPLX(REAL(fp),AIMAG(fm),kind=real_8)
    ENDDO
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE hfxaa
  ! ==================================================================
  SUBROUTINE hfxab(ehfx,pf,psia,a_stored_in_real,&
       psib,b_stored_in_real,vpotg,vpotr,psic,c2a,c2b)
    ! ==================================================================
    REAL(real_8)                             :: ehfx, pf
    COMPLEX(real_8)                          :: psia(:)
    LOGICAL                                  :: a_stored_in_real
    COMPLEX(real_8)                          :: psib(:)
    LOGICAL                                  :: b_stored_in_real
    COMPLEX(real_8)                          :: vpotg(:)
    REAL(real_8)                             :: vpotr(:)
    COMPLEX(real_8)                          :: psic(:), c2a(:), c2b(:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'hfxab'

    COMPLEX(real_8)                          :: fm, fp
    INTEGER                                  :: ig, ir, isub

    ! KPT replace NGW -> NGWK
    ! Variables
    ! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    ehfx=0.0_real_8
    IF (a_stored_in_real) THEN
       IF (b_stored_in_real) THEN
          !$omp parallel do private(IR)
          DO ir=1,llr1
             psic(ir)=REAL(psia(ir))*REAL(psib(ir))
          ENDDO
       ELSE
          !$omp parallel do private(IR)
          DO ir=1,llr1
             psic(ir)=REAL(psia(ir))*AIMAG(psib(ir))
          ENDDO
       ENDIF
    ELSE
       IF (b_stored_in_real) THEN
          !$omp parallel do private(IR)
          DO ir=1,llr1
             psic(ir)=AIMAG(psia(ir))*REAL(psib(ir))
          ENDDO
       ELSE
          !$omp parallel do private(IR)
          DO ir=1,llr1
             psic(ir)=AIMAG(psia(ir))*AIMAG(psib(ir))
          ENDDO
       ENDIF
    ENDIF
    CALL dscal(2*llr1,1._real_8/parm%omega,psic,1)


    CALL fwfftn(psic,.FALSE.,parai%allgrp)


    !$omp parallel do private(IG) &
    !$omp  reduction(+:ehfx)
    DO ig=1,jhg
       vpotg(ig)=-pf*scgx(ig)*psic(nzff(ig))
       ehfx=ehfx+4._real_8*REAL(vpotg(ig)*CONJG(psic(nzff(ig))))
    ENDDO
    IF (geq0) ehfx=ehfx-2._real_8*REAL(vpotg(1)*CONJG(psic(nzff(1))))
    CALL zeroing(psic)
    !ocl novrec
    !$omp parallel do private(IG)
    DO ig=1,jhg
       psic(nzff(ig))=vpotg(ig)
       psic(inzf(ig))=CONJG(vpotg(ig))
    ENDDO
    IF (geq0) psic(nzff(1))=vpotg(1)


    CALL invfftn(psic,.FALSE.,parai%allgrp)


    !$omp parallel do private(IR)
    DO ir=1,llr1
       vpotr(ir)=REAL(psic(ir))
    ENDDO

    IF (a_stored_in_real) THEN
       IF (b_stored_in_real) THEN
          !$omp parallel do private(IR)
          DO ir=1,llr1
             psic(ir)=vpotr(ir)*(REAL(psia(ir))&
                  +uimag*REAL(psib(ir)))
          ENDDO
       ELSE
          !$omp parallel do private(IR)
          DO ir=1,llr1
             psic(ir)=vpotr(ir)*(REAL(psia(ir))&
                  +uimag*AIMAG(psib(ir)))
          ENDDO
       ENDIF
    ELSE
       IF (b_stored_in_real) THEN
          !$omp parallel do private(IR)
          DO ir=1,llr1
             psic(ir)=vpotr(ir)*(AIMAG(psia(ir))&
                  +uimag*REAL(psib(ir)))
          ENDDO
       ELSE
          !$omp parallel do private(IR)
          DO ir=1,llr1
             psic(ir)=vpotr(ir)*(AIMAG(psia(ir))&
                  +uimag*AIMAG(psib(ir)))
          ENDDO
       ENDIF
    ENDIF


    CALL fwfftn(psic,.TRUE.,parai%allgrp)

    !$omp parallel do private(IG,FP,FM)
    DO ig=1,jgw
       fp=psic(nzfs(ig))+psic(inzs(ig))
       fm=psic(nzfs(ig))-psic(inzs(ig))
       c2b(ig)=c2b(ig)-CMPLX(REAL(fp),AIMAG(fm),kind=real_8)
       c2a(ig)=c2a(ig)-CMPLX(AIMAG(fp),-REAL(fm),kind=real_8)
    ENDDO


    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE hfxab
  ! ==================================================================
  SUBROUTINE hfxab2(ehfx_1,ehfx_2,pf1,pf2,psia,a_stored_in_real,&
       psib,vpotg,vpotr,psic,c2a,c2b1,c2b2)
    ! ==================================================================
    REAL(real_8)                             :: ehfx_1, ehfx_2, pf1, pf2
    COMPLEX(real_8)                          :: psia(:)
    LOGICAL                                  :: a_stored_in_real
    COMPLEX(real_8)                          :: psib(:), vpotg(:,:)
    REAL(real_8)                             :: vpotr(llr1,2)
    COMPLEX(real_8)                          :: psic(:), c2a(:), c2b1(:), &
                                                c2b2(:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'hfxab2'

    COMPLEX(real_8)                          :: fm, fp
    INTEGER                                  :: ig, ir, isub

    CALL tiset(procedureN,isub)
    ehfx_1=0.0_real_8
    ehfx_2=0.0_real_8
    IF (a_stored_in_real) THEN
       !$omp parallel do private(IR)
       DO ir=1,llr1
          psic(ir)=REAL(psia(ir))*psib(ir)
       ENDDO
    ELSE
       !$omp parallel do private(IR)
       DO ir=1,llr1
          psic(ir)=AIMAG(psia(ir))*psib(ir)
       ENDDO
    ENDIF
    CALL dscal(2*llr1,1._real_8/parm%omega,psic,1)


    CALL fwfftn(psic,.FALSE.,parai%allgrp)

    !$omp parallel do private(IG,FP,FM) &
    !$omp  reduction(+:ehfx_1,ehfx_2)
    DO ig=1,jhg
       fp=psic(nzff(ig))+psic(inzf(ig))
       fm=psic(nzff(ig))-psic(inzf(ig))
       vpotg(ig,1)=-pf1*scgx(ig)*0.5_real_8*CMPLX(REAL(fp),AIMAG(fm),kind=real_8)
       vpotg(ig,2)=-pf2*scgx(ig)*0.5_real_8*CMPLX(AIMAG(fp),-REAL(fm),kind=real_8)
       ehfx_1=ehfx_1+2._real_8*REAL(vpotg(ig,1)&
            *CMPLX(REAL(fp),-AIMAG(fm),kind=real_8))
       ehfx_2=ehfx_2+2._real_8*REAL(vpotg(ig,2)&
            *CMPLX(AIMAG(fp),REAL(fm),kind=real_8))
    ENDDO
    IF (geq0) THEN
       fp=psic(nzff(1))+psic(inzf(1))
       fm=psic(nzff(1))-psic(inzf(1))
       ehfx_1=ehfx_1-REAL(vpotg(1,1)*CMPLX(REAL(fp),-AIMAG(fm),kind=real_8))
       ehfx_2=ehfx_2-REAL(vpotg(1,2)*CMPLX(AIMAG(fp),REAL(fm),kind=real_8))
    ENDIF
    CALL zeroing(psic)
    !ocl novrec
    !$omp parallel do private(IG)
    DO ig=1,jhg
       psic(nzff(ig))=vpotg(ig,1)+uimag*vpotg(ig,2)
       psic(inzf(ig))=CONJG(vpotg(ig,1))+uimag*CONJG(vpotg(ig,2))
    ENDDO
    IF (geq0) psic(nzff(1))=vpotg(1,1)+uimag*vpotg(1,2)


    CALL invfftn(psic,.FALSE.,parai%allgrp)


    !$omp parallel do private(IR)
    DO ir=1,llr1
       vpotr(ir,1)=REAL(psic(ir))
       vpotr(ir,2)=AIMAG(psic(ir))
    ENDDO



    IF (a_stored_in_real) THEN
       !$omp parallel do private(IR)
       DO ir=1,llr1
          psic(ir)=vpotr(ir,1)*(REAL(psia(ir))&
               +uimag*REAL(psib(ir)))
       ENDDO
    ELSE
       !$omp parallel do private(IR)
       DO ir=1,llr1
          psic(ir)=vpotr(ir,1)*(AIMAG(psia(ir))&
               +uimag*REAL(psib(ir)))
       ENDDO
    ENDIF


    CALL fwfftn(psic,.TRUE.,parai%allgrp)

    !$omp parallel do private(IG,FP,FM)
    DO ig=1,jgw
       fp=psic(nzfs(ig))+psic(inzs(ig))
       fm=psic(nzfs(ig))-psic(inzs(ig))
       c2b1(ig)=c2b1(ig)-CMPLX(REAL(fp),AIMAG(fm),kind=real_8)
       c2a(ig)=c2a(ig)-CMPLX(AIMAG(fp),-REAL(fm),kind=real_8)
    ENDDO


    IF (a_stored_in_real) THEN
       !$omp parallel do private(IR)
       DO ir=1,llr1
          psic(ir)=vpotr(ir,2)*(REAL(psia(ir))&
               +uimag*AIMAG(psib(ir)))
       ENDDO
    ELSE
       !$omp parallel do private(IR)
       DO ir=1,llr1
          psic(ir)=vpotr(ir,2)*(AIMAG(psia(ir))&
               +uimag*AIMAG(psib(ir)))
       ENDDO
    ENDIF


    CALL fwfftn(psic,.TRUE.,parai%allgrp)

    !$omp parallel do private(IG,FP,FM)
    DO ig=1,jgw
       fp=psic(nzfs(ig))+psic(inzs(ig))
       fm=psic(nzfs(ig))-psic(inzs(ig))
       c2b2(ig)=c2b2(ig)-CMPLX(REAL(fp),AIMAG(fm),kind=real_8)
       c2a(ig)=c2a(ig)-CMPLX(AIMAG(fp),-REAL(fm),kind=real_8)
    ENDDO
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE hfxab2
  ! ==================================================================

  ! ==================================================================
  ! New routines for Coulomb potential between orbital pairs
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE hfxaa_new(ehfx,ehfx_nothresh,pf,psia,a_stored_in_real,&
       vpotg,vpotr,psic,c2a)
    ! ==--------------------------------------------------------------==

    REAL(real_8), INTENT(out)                :: ehfx, ehfx_nothresh
    REAL(real_8), INTENT(in)                 :: pf
    COMPLEX(real_8)                          :: psia(:)
    LOGICAL, INTENT(in)                      :: a_stored_in_real
    COMPLEX(real_8)                          :: vpotg(:)
    REAL(real_8)                             :: vpotr(:)
    COMPLEX(real_8)                          :: psic(:), c2a(:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'hfxaa_new'
    REAL(real_8), PARAMETER                  :: ef = 1.0_real_8

    INTEGER                                  :: isub
    LOGICAL                                  :: too_small_int
    REAL(real_8)                             :: ehfx_sum

    CALL tiset(procedureN,isub)
    ehfx_nothresh = 0.0_real_8
    too_small_int = .FALSE.
   
    CALL hfx_get_pair_density(psia,psia,psic,a_stored_in_real,b_stored_in_real=a_stored_in_real)
    CALL fwfftn(psic,.FALSE.,parai%allgrp)
    CALL hfxaa_ab_get_ehfx(psic,vpotg,pf,ef,ehfx_nothresh)
    IF (hfxc3%twscr) too_small_int = check_int(ehfx_nothresh)

    IF(too_small_int) THEN
       ehfx = 0.0_real_8
    ELSE
       ehfx = ehfx_nothresh
       CALL hfx_set_vpotg(psic,vpotg)
       CALL invfftn(psic,.FALSE.,parai%allgrp)
       CALL hfx_set_vpotr(psic,vpotr)
       CALL hfx_get_potential(psia,psia,psic,vpotr,a_stored_in_real)
       IF (cntl%use_scaled_hfx) THEN
          CALL hfx_get_c2_real_space(psic,c2a)
          ehfx = scex_lambda*ehfx
       ELSE
          CALL fwfftn(psic,.TRUE.,parai%allgrp)
          CALL hfx_get_c2(psic,c2a)
       ENDIF
    ENDIF

    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE hfxaa_new
  ! ==================================================================
  SUBROUTINE hfxab_new(ehfx,ehfx_nothresh,pf,psia,a_stored_in_real,psib,b_stored_in_real,&
       vpotg,vpotr,psic,c2a,c2b)
    ! ==--------------------------------------------------------------==

    REAL(real_8), INTENT(out)                :: ehfx, ehfx_nothresh
    REAL(real_8), INTENT(in)                 :: pf
    COMPLEX(real_8)                          :: psia(:)
    LOGICAL                                  :: a_stored_in_real
    COMPLEX(real_8)                          :: psib(:)
    LOGICAL                                  :: b_stored_in_real
    COMPLEX(real_8)                          :: vpotg(:)
    REAL(real_8)                             :: vpotr(:)
    COMPLEX(real_8)                          :: psic(:), c2a(:), c2b(:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'hfxab_new'
    REAL(real_8), PARAMETER                  :: ef = 2.0_real_8

    INTEGER                                  :: isub
    LOGICAL                                  :: too_small_int
    REAL(real_8)                             :: ehfx_sum

    CALL tiset(procedureN,isub)
    ehfx_nothresh = 0.0_real_8
    too_small_int = .FALSE.

    CALL hfx_get_pair_density(psia,psib,psic,a_stored_in_real,b_stored_in_real=b_stored_in_real)
    CALL fwfftn(psic,.FALSE.,parai%allgrp)
    CALL hfxaa_ab_get_ehfx(psic,vpotg,pf,ef,ehfx_nothresh)
    IF (hfxc3%twscr) too_small_int = check_int(ehfx_nothresh)

    IF(too_small_int) THEN
       ehfx = 0.0_real_8
    ELSE
       ehfx = ehfx_nothresh
       CALL hfx_set_vpotg(psic,vpotg)
       CALL invfftn(psic,.FALSE.,parai%allgrp)
       CALL hfx_set_vpotr(psic,vpotr)
       CALL hfx_get_potential(psia,psib,psic,vpotr,a_stored_in_real,b_stored_in_real=b_stored_in_real)
       IF (cntl%use_scaled_hfx) THEN
          CALL hfx_get_c2_real_space(psic,c2b,c2_2=c2a)
          ehfx = scex_lambda*ehfx
       ELSE
          CALL fwfftn(psic,.TRUE.,parai%allgrp)
          CALL hfx_get_c2(psic,c2b,c2_2=c2a)
       ENDIF
    ENDIF

    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE hfxab_new
  ! ==================================================================
  SUBROUTINE hfxab2_new(ehfx_1,ehfx_2,ehfx_nothresh_1,ehfx_nothresh_2,&
       pf1,pf2,psia,a_stored_in_real,psib,vpotg,vpotr,psic,c2a,c2b1,c2b2)
    ! ==--------------------------------------------------------------==

    REAL(real_8), INTENT(out)                :: ehfx_1, ehfx_2, &
                                                ehfx_nothresh_1, &
                                                ehfx_nothresh_2
    REAL(real_8), INTENT(in)                 :: pf1, pf2
    COMPLEX(real_8)                          :: psia(:)
    LOGICAL                                  :: a_stored_in_real
    COMPLEX(real_8)                          :: psib(:), vpotg(:,:)
    REAL(real_8)                             :: vpotr(:,:)
    COMPLEX(real_8)                          :: psic(:), c2a(:), c2b1(:), &
                                                c2b2(:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'hfxab2_new'
    LOGICAL, PARAMETER                       :: b_stored_in_real_1 = .TRUE.
    LOGICAL, PARAMETER                       :: b_stored_in_real_2 = .FALSE.

    COMPLEX(real_8)                          :: fm, fp
    INTEGER                                  :: ig, ir, isub
    LOGICAL                                  :: too_small_int_1, &
                                                too_small_int_2
    REAL(real_8)                             :: ehfx_sum(2)

    CALL tiset(procedureN,isub)
    ehfx_nothresh_1 = 0.0_real_8
    ehfx_nothresh_2 = 0.0_real_8
    too_small_int_1 = .FALSE.
    too_small_int_2 = .FALSE.

    CALL hfx_get_pair_density(psia,psib,psic,a_stored_in_real)
    CALL fwfftn(psic,.FALSE.,parai%allgrp)
    CALL hfxab2_get_ehfx(psic,vpotg,pf1,pf2,ehfx_nothresh_1,ehfx_nothresh_2)

    too_small_int_1 = check_int(ehfx_nothresh_1)
    too_small_int_2 = check_int(ehfx_nothresh_2)

    IF (too_small_int_1 .AND. too_small_int_2) THEN
       ehfx_1 = 0.0_real_8
       ehfx_2 = 0.0_real_8
    ELSE
       CALL hfx_set_vpotg(psic,vpotg(:,1),vpotg_2=vpotg(:,2))
       CALL invfftn(psic,.FALSE.,parai%allgrp)
       CALL hfx_set_vpotr(psic,vpotr(:,1),vpotr_2=vpotr(:,2))

       IF(too_small_int_1) THEN
          ehfx_1 = 0.0_real_8
       ELSE
          ehfx_1 = ehfx_nothresh_1
          CALL hfx_get_potential(psia,psib,psic,vpotr(:,1),a_stored_in_real,b_stored_in_real=b_stored_in_real_1)
          IF (cntl%use_scaled_hfx) THEN
             CALL hfx_get_c2_real_space(psic,c2b1,c2_2=c2a)
             ehfx_1 = scex_lambda*ehfx_1
          ELSE
             CALL fwfftn(psic,.TRUE.,parai%allgrp)
             CALL hfx_get_c2(psic,c2b1,c2_2=c2a)
          ENDIF
       ENDIF

       IF(too_small_int_2) THEN
          ehfx_2 = 0.0_real_8
       ELSE
          ehfx_2 = ehfx_nothresh_2
          CALL hfx_get_potential(psia,psib,psic,vpotr(:,2),a_stored_in_real,b_stored_in_real=b_stored_in_real_2)
          IF (cntl%use_scaled_hfx) THEN
             CALL hfx_get_c2_real_space(psic,c2b2,c2_2=c2a)
             ehfx_2 = scex_lambda*ehfx_2
          ELSE
             CALL fwfftn(psic,.TRUE.,parai%allgrp)
             CALL hfx_get_c2(psic,c2b2,c2_2=c2a)
          ENDIF
       ENDIF
    ENDIF

    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE hfxab2_new
  ! ==================================================================

  ! ==================================================================
  ! Generic routines for hfxaa, hfxab, hfxab2
  ! Refactored out of hfx..._new
  !                             21.11.2018 - M. P. Bircher @ LCBC/EPFL
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE hfx_set_vpotg(psic,vpotg,vpotg_2)
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8), DIMENSION(:), &
                     INTENT(inout)  :: psic, vpotg
    COMPLEX(real_8), DIMENSION(:), OPTIONAL, &
                     INTENT(inout)  :: vpotg_2

    INTEGER                         :: ig

    CALL zeroing(psic)

    IF (present(vpotg_2)) THEN
       !ocl novrec
       !$omp parallel do private(IG)
       DO ig=1,jhg
          psic(nzff(ig))=vpotg(ig)+uimag*vpotg_2(ig)
          psic(inzf(ig))=CONJG(vpotg(ig))+uimag*CONJG(vpotg_2(ig))
       ENDDO
       IF (geq0) psic(nzff(1))=vpotg(1)+uimag*vpotg_2(1)
    ELSE
       !ocl novrec
       !$omp parallel do private(IG)
       DO ig=1,jhg
          psic(nzff(ig))=vpotg(ig)
          psic(inzf(ig))=CONJG(vpotg(ig))
       ENDDO
       IF (geq0) psic(nzff(1))=vpotg(1)
    ENDIF

    ! ==--------------------------------------------------------------==
  END SUBROUTINE hfx_set_vpotg
  ! ==================================================================
  SUBROUTINE hfx_get_pair_density(psia,psib,psic,a_stored_in_real,b_stored_in_real)
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8), DIMENSION(:), &
                     INTENT(in)     :: psia, psib
    COMPLEX(real_8), DIMENSION(:), &
                     INTENT(inout)  :: psic
    LOGICAL, INTENT(in)             :: a_stored_in_real
    LOGICAL, INTENT(in), OPTIONAL   :: b_stored_in_real

    INTEGER                         :: ir

    IF (a_stored_in_real) THEN
       IF (present(b_stored_in_real)) THEN
          IF (b_stored_in_real) THEN
             !$omp parallel do private(IR)
             DO ir=1,llr1
                psic(ir)=REAL(psia(ir))*REAL(psib(ir))
             ENDDO
          ELSE
             !$omp parallel do private(IR)
             DO ir=1,llr1
                psic(ir)=REAL(psia(ir))*AIMAG(psib(ir))
             ENDDO
          ENDIF
       ELSE
          !$omp parallel do private(IR)
          DO ir=1,llr1
             psic(ir)=REAL(psia(ir))*psib(ir)
          ENDDO
       ENDIF
    ELSE
       IF (present(b_stored_in_real)) THEN
          IF (b_stored_in_real) THEN
             !$omp parallel do private(IR)
             DO ir=1,llr1
                psic(ir)=AIMAG(psia(ir))*REAL(psib(ir))
             ENDDO
          ELSE
             !$omp parallel do private(IR)
             DO ir=1,llr1
                psic(ir)=AIMAG(psia(ir))*AIMAG(psib(ir))
             ENDDO
          ENDIF
       ELSE
          !$omp parallel do private(IR)
          DO ir=1,llr1
             psic(ir)=AIMAG(psia(ir))*psib(ir)
          ENDDO
       ENDIF
    ENDIF
    CALL dscal(2*llr1,1._real_8/parm%omega,psic,1)

    ! ==--------------------------------------------------------------==
  END SUBROUTINE hfx_get_pair_density
  ! ==================================================================
  SUBROUTINE hfx_set_vpotr(psic,vpotr,vpotr_2)
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8), DIMENSION(:), &
                     INTENT(in)     :: psic
    REAL(real_8), DIMENSION(:), &
                  INTENT(inout)     :: vpotr
    REAL(real_8), DIMENSION(:), OPTIONAL, &
                  INTENT(inout)     :: vpotr_2

    INTEGER                         :: ir

    IF (present(vpotr_2)) THEN
       !$omp parallel do private(IR)
       DO ir=1,llr1
          vpotr(ir)   = REAL(psic(ir))
          vpotr_2(ir) = AIMAG(psic(ir))
       ENDDO
    ELSE
       !$omp parallel do private(IR)
       DO ir=1,llr1
          vpotr(ir) = REAL(psic(ir))
       ENDDO
    ENDIF

    ! ==--------------------------------------------------------------==
  END SUBROUTINE hfx_set_vpotr
  ! ==================================================================
  SUBROUTINE hfx_get_potential(psia,psib,psic,vpotr,a_stored_in_real,b_stored_in_real)
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8), DIMENSION(:), &
                     INTENT(in)     :: psia, psib
    COMPLEX(real_8), DIMENSION(:), &
                     INTENT(inout)  :: psic
    REAL(real_8), DIMENSION(:), &
                  INTENT(in)        :: vpotr
    LOGICAL, INTENT(in)             :: a_stored_in_real
    LOGICAL, INTENT(in), OPTIONAL   :: b_stored_in_real

    INTEGER                         :: ir

    IF (a_stored_in_real) THEN
       IF (present(b_stored_in_real)) THEN
          IF (b_stored_in_real) THEN
             !$omp parallel do private(IR)
             DO ir=1,llr1
                psic(ir)=vpotr(ir)*(REAL(psia(ir))&
                     +uimag*REAL(psib(ir)))
             ENDDO
          ELSE
             !$omp parallel do private(IR)
             DO ir=1,llr1
                psic(ir)=vpotr(ir)*(REAL(psia(ir))&
                     +uimag*AIMAG(psib(ir)))
             ENDDO
          ENDIF
       ELSE
          !$omp parallel do private(IR)
          DO ir=1,llr1
             psic(ir)=vpotr(ir)*REAL(psia(ir))
          ENDDO
       ENDIF
    ELSE
       IF (present(b_stored_in_real)) THEN
          IF (b_stored_in_real) THEN
             !$omp parallel do private(IR)
             DO ir=1,llr1
                psic(ir)=vpotr(ir)*(AIMAG(psia(ir))&
                     +uimag*REAL(psib(ir)))
             ENDDO
          ELSE
             !$omp parallel do private(IR)
             DO ir=1,llr1
                psic(ir)=vpotr(ir)*(AIMAG(psia(ir))&
                     +uimag*AIMAG(psib(ir)))
             ENDDO
          ENDIF
       ELSE
          !$omp parallel do private(IR)
          DO ir=1,llr1
             psic(ir)=vpotr(ir)*AIMAG(psia(ir))
          ENDDO
       ENDIF
    ENDIF

    ! ==--------------------------------------------------------------==
  END SUBROUTINE hfx_get_potential
  ! ==================================================================
  SUBROUTINE hfx_get_c2_real_space(psic,c2,c2_2)
    ! ==--------------------------------------------------------------==
          
    COMPLEX(real_8), DIMENSION(:), &
                     INTENT(in)       :: psic
    COMPLEX(real_8), DIMENSION(:), &
                     INTENT(inout)    :: c2
    COMPLEX(real_8), DIMENSION(:), OPTIONAL, &
                     INTENT(inout)    :: c2_2

    INTEGER         :: ir

    IF (present(c2_2)) THEN
       !$omp parallel do private(ir)
       DO ir=1,llr1
          c2(ir)   = c2(ir)   + REAL(psic(ir)) 
          c2_2(ir) = c2_2(ir) + AIMAG(psic(ir)) 
       ENDDO
    ELSE
       !$omp parallel do private(ir)
       DO ir=1,llr1
          c2(ir)   = c2(ir) + REAL(psic(ir))
       ENDDO
    ENDIF

    ! ==--------------------------------------------------------------==
  END SUBROUTINE hfx_get_c2_real_space
  ! ================================================================== 
  SUBROUTINE hfx_get_c2(psic,c2,c2_2)
    ! ==--------------------------------------------------------------==
          
    COMPLEX(real_8), DIMENSION(:), &
                     INTENT(in)       :: psic
    COMPLEX(real_8), DIMENSION(:), &
                     INTENT(inout)    :: c2
    COMPLEX(real_8), DIMENSION(:), OPTIONAL, &
                     INTENT(inout)    :: c2_2

    COMPLEX(real_8) :: fp, fm
    INTEGER         :: ig

    IF (present(c2_2)) THEN
       !$omp parallel do private(IG,FP,FM)
       DO ig=1,jgw
          fp = psic(nzfs(ig)) + psic(inzs(ig))
          fm = psic(nzfs(ig)) - psic(inzs(ig))
          c2(ig)   = c2(ig)   - CMPLX(REAL(fp),AIMAG(fm),kind=real_8)
          c2_2(ig) = c2_2(ig) - CMPLX(AIMAG(fp),-REAL(fm),kind=real_8)
       ENDDO
    ELSE
       !$omp parallel do private(IG,FP,FM)
       DO ig=1,jgw
          fp = psic(nzfs(ig)) + psic(inzs(ig))
          fm = psic(nzfs(ig)) - psic(inzs(ig))
          c2(ig)   = c2(ig)   - CMPLX(REAL(fp),AIMAG(fm),kind=real_8)
       ENDDO
    ENDIF

    ! ==--------------------------------------------------------------==
  END SUBROUTINE hfx_get_c2
  ! ================================================================== 
  LOGICAL FUNCTION check_int(ehfx_nothresh) RESULT(too_small_int)
    ! ==--------------------------------------------------------------==

    REAL(real_8), INTENT(in) :: ehfx_nothresh
    REAL(real_8)             :: ehfx_sum

    ehfx_sum = ehfx_nothresh
    CALL mp_sum(ehfx_sum,parai%allgrp)
    too_small_int = ABS(ehfx_sum) < hfxc4%dwf_integral_thresh

    ! ==--------------------------------------------------------------==
  END FUNCTION
  ! ================================================================== 

  ! ==================================================================
  ! Routines that are not generalised
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE hfxaa_ab_get_ehfx(psic,vpotg,pf,ef,ehfx_nothresh)
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8), DIMENSION(:), &
                     INTENT(in)     :: psic
    COMPLEX(real_8), DIMENSION(:), &
                     INTENT(inout)  :: vpotg
    REAL(real_8), INTENT(inout)     :: ehfx_nothresh
    REAL(real_8), INTENT(in)        :: pf, ef

    INTEGER                         :: ig
    REAL(real_8)                    :: ef_2
    COMPLEX(real_8)                 :: fp

    ef_2 = ef + ef

    !$omp parallel do private(IG,FP) &
    !$omp  reduction(+:ehfx_nothresh)
    DO ig=1,jhg
       fp=psic(nzff(ig))
       vpotg(ig)=-pf*scgx(ig)*fp
       ehfx_nothresh=ehfx_nothresh+REAL(vpotg(ig)*CONJG(fp))
    ENDDO
    ehfx_nothresh = ehfx_nothresh*ef_2
    IF (geq0) ehfx_nothresh=ehfx_nothresh-ef*REAL(vpotg(1)*CONJG(psic(nzff(1))))

    ! ==--------------------------------------------------------------==
  END SUBROUTINE hfxaa_ab_get_ehfx
  ! ==================================================================
  SUBROUTINE hfxab2_get_ehfx(psic,vpotg,pf1,pf2,ehfx_nothresh_1,ehfx_nothresh_2)
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8), DIMENSION(:), &
                     INTENT(in)     :: psic
    COMPLEX(real_8), DIMENSION(:,:), &
                     INTENT(inout)  :: vpotg
    REAL(real_8), INTENT(inout)     :: ehfx_nothresh_1, ehfx_nothresh_2
    REAL(real_8), INTENT(in)        :: pf1,pf2

    INTEGER                         :: ig
    COMPLEX(real_8)                 :: fp,fm

    !$omp parallel do private(IG,FP,FM) &
    !$omp  reduction(+:ehfx_nothresh_1,ehfx_nothresh_2)
    DO ig=1,jhg
       fp=psic(nzff(ig))+psic(inzf(ig))
       fm=psic(nzff(ig))-psic(inzf(ig))
       vpotg(ig,1)=-pf1*scgx(ig)*0.5_real_8*CMPLX(REAL(fp),AIMAG(fm),kind=real_8)
       vpotg(ig,2)=-pf2*scgx(ig)*0.5_real_8*CMPLX(AIMAG(fp),-REAL(fm),kind=real_8)
       ehfx_nothresh_1=ehfx_nothresh_1+2._real_8*REAL(vpotg(ig,1)*CMPLX(REAL(fp),-AIMAG(fm),kind=real_8))
       ehfx_nothresh_2=ehfx_nothresh_2+2._real_8*REAL(vpotg(ig,2)*CMPLX(AIMAG(fp),REAL(fm),kind=real_8))
    ENDDO
    IF (geq0) THEN
       fp=psic(nzff(1))+psic(inzf(1))
       fm=psic(nzff(1))-psic(inzf(1))
       ehfx_nothresh_1=ehfx_nothresh_1-REAL(vpotg(1,1)*CMPLX(REAL(fp),-AIMAG(fm),kind=real_8))
       ehfx_nothresh_2=ehfx_nothresh_2-REAL(vpotg(1,2)*CMPLX(AIMAG(fp),REAL(fm),kind=real_8))
    ENDIF

    ! ==--------------------------------------------------------------==
  END SUBROUTINE hfxab2_get_ehfx
  ! ==================================================================

  ! ================================================================== 
  ! Routines for list
  ! ================================================================== 

  ! ================================================================== 
  PURE SUBROUTINE sub2ind(r,c,ij)
    ! ==--------------------------------------------------------------==
    INTEGER, INTENT(in)                      :: r, c
    INTEGER, INTENT(out)                     :: ij

    INTEGER                                  :: i, j

    i = MIN(c,r)
    j = MAX(c,r)
    ij = j*(j-1)/2+i
    ! ==--------------------------------------------------------------==
  END SUBROUTINE sub2ind
  ! ================================================================== 
  PURE SUBROUTINE ind2sub(i,r,c)
    ! ==--------------------------------------------------------------==
    INTEGER, INTENT(in)                      :: i
    INTEGER, INTENT(out)                     :: r, c

    c = CEILING((-1.0_real_8 + SQRT(1.0_real_8+8.0_real_8*REAL(i,kind=real_8)))/2.0_real_8)
    r = i - c * ( c - 1 ) / 2
    ! ==--------------------------------------------------------------==
  END SUBROUTINE ind2sub
  ! ================================================================== 
  SUBROUTINE ind2sub_rect(n,i,r,c)
    ! ==--------------------------------------------------------------==
    INTEGER, INTENT(in)                      :: n, i
    INTEGER, INTENT(out)                     :: r, c

    c = CEILING( REAL(i,kind=real_8) / REAL(n,kind=real_8) )
    r = i - ( c - 1 ) * n
    ! ==--------------------------------------------------------------==
  END SUBROUTINE ind2sub_rect
  ! ================================================================== 


  ! ==--------------------------------------------------------------==
END MODULE pw_hfx
! ================================================================== 

