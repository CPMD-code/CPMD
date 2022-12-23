! ==================================================================
MODULE pw_hfx_resp
  USE cnst,                            ONLY: uimag
  USE cp_grp_utils,                    ONLY: cp_grp_redist
  USE cppt,                            ONLY: scgx
  USE elct,                            ONLY: crge
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
  USE hfxmod,                          ONLY: hfxc5,&
                                             ipoolhfx
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai,&
                                             paral
  USE pslo,                            ONLY: pslo_com
  USE pw_hfx,                          ONLY: &
       block_invfft_old, build_block_info, build_dist, dist_t, ind2sub_rect, &
       iterator_blocks, iterator_next_block, iterator_start, iterator_t, &
       release_dist
  USE pw_hfx_input_cnst,               ONLY: hfx_dist_dynamic
  USE pw_hfx_resp_types,               ONLY: hfx_lin2,&
                                             hfx_resp,&
                                             hfx_resp_env
  USE pw_hfx_resp_utils,               ONLY: pw_hfx_resp_create,&
                                             pw_hfx_resp_destroy
  USE scex_utils,                      ONLY: scex,&
                                             scex_lambda,&
                                             scex_ID_parent,&
                                             scex_ID_scaled
  USE spin,                            ONLY: lspin2,&
                                             spin_mod
  USE state_utils,                     ONLY: add_wfn
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             group,&
                                             ncpw,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  ! ==--------------------------------------------------------------==

  PRIVATE

  PUBLIC :: hfxpsi_new
  PUBLIC :: hfxrpa_new
  PUBLIC :: hfxrpa_incore
  PUBLIC :: hfx_resp_init
  PUBLIC :: hfx_resp_finalize

CONTAINS

  ! ==================================================================
  SUBROUTINE hfx_resp_init
    INTEGER                                  :: nstate

    nstate=crge%n
    IF (hfx_resp_env%store_v_ab) CALL pw_hfx_resp_create(hfx_resp,fpar%nnr1,nstate)
  END SUBROUTINE hfx_resp_init

  ! ==================================================================
  SUBROUTINE hfx_resp_finalize
    IF (hfx_resp_env%store_v_ab) CALL pw_hfx_resp_destroy(hfx_resp)
  END SUBROUTINE hfx_resp_finalize
  ! ==================================================================
  SUBROUTINE hfxpsi_new(c0,cpsi,c2,f,sign,psic,nstate,norb)
    ! ==--------------------------------------------------------------==
    ! Hartree-Fock exchange
    ! Refactored to include coordinate-scaling on 12.02.2019
    !                                           M.P. Bircher @ LCBC/EPFL
    ! ==--------------------------------------------------------------==
    ! c0   : ground state orbitals (1,nstate+nvir) 
    ! cpsi : LR-orbitals (1,nstate=occ states)
    ! c2   : Potential applied to cpsi
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8)                          :: c0(:,:), cpsi(:,:), c2(:,:)
    REAL(real_8)                             :: f(:), sign
    COMPLEX(real_8)                          :: psic(:)
    INTEGER                                  :: nstate, norb

    CHARACTER(*), PARAMETER                  :: procedureN = 'hfxpsi_new'
    COMPLEX(real_8), PARAMETER               :: zone = (1.0_real_8,0.0_real_8)

    COMPLEX(real_8), ALLOCATABLE             :: C2_x(:,:)
    COMPLEX(real_8), ALLOCATABLE, &
      DIMENSION(:)                           :: vpotg
    COMPLEX(real_8), ALLOCATABLE, &
      DIMENSION(:, :)                        :: C2_hfx, psi_col, psi_row
    INTEGER :: blk, c, col_blk, col_offset, col_size, hfx_dist, i, ia, ia1, &
      ia2, ib, ib1, ib2, icount, ierr, ik, isub, isub2, j, lwork, &
      max_block_size, nblks, nbr_blk_orb, nbr_blk_state, nvir, r, row_blk, &
      row_offset, row_size
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: offsets_orb, offsets_state, &
                                                sizes_orb, sizes_state
    INTEGER, ALLOCATABLE, DIMENSION(:), SAVE :: tick_blks
    INTEGER, SAVE                            :: icall = 0
    LOGICAL                                  :: cholesky_decomp = .FALSE., &
                                                svd_decomp = .TRUE.
    LOGICAL, SAVE                            :: is_first_call = .TRUE.
    REAL(real_8)                             :: max_sigma, min_sigma, norm, &
                                                pfl, pfx
    REAL(real_8), ALLOCATABLE, DIMENSION(:)  :: sigma, vpotr, work
    REAL(real_8), ALLOCATABLE, &
      DIMENSION(:, :)                        :: local_h, u, vt
    REAL(real_8), EXTERNAL                   :: dotp
    TYPE(dist_t)                             :: dist
    TYPE(iterator_t)                         :: iter

! containes normalized c2 for testing

    svd_decomp=.NOT.cholesky_decomp
    hfx_resp_env%sigma_eps=0.00001_real_8
    hfx_resp_env%interval_restart_lin2=4

    IF (func1%mhfx.EQ.0) RETURN

    nvir=crge%n-nstate
    !
    ! Commented by Martin, keep output compact for release
    !
    ! IF (paral%io_parent) THEN
    !    WRITE(*,*) 'nstate', nstate
    !    WRITE(*,*) 'nvir',   nvir
    !    WRITE(*,*) 'norb',   norb
    !    WRITE(*,*) 'icall', icall
    ! ENDIF

    icall=icall+1
    IF (MOD(icall,hfx_resp_env%interval_restart_lin2).EQ.0)  &
         hfx_resp_env%lin2_restart=.TRUE.
    ! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)

    IF (func1%mhfx.EQ.2) CALL stopgm(procedureN,&
         'HARTREE METHOD NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    IF (tkpts%tkpnt) CALL stopgm(procedureN,'K-POINTS NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    IF (lspin2%tlse) CALL stopgm(procedureN,'NO LSE POSSIBLE',& 
         __LINE__,__FILE__)
    IF (pslo_com%tivan) CALL stopgm(procedureN,'NO VDB PP POSSIBLE',& 
         __LINE__,__FILE__)
    IF (group%nogrp.GT.1) CALL stopgm(procedureN,&
         'OLD TASK GROUPS NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    !
    IF (cntl%use_scaled_hfx) THEN
       IF (.NOT. scex%init) CALL stopgm(procedureN,'ScEX requested, but not initialised',&
                                        __LINE__,__FILE__)
       CALL setfftn(scex_ID_parent)
    ELSE
       CALL setfftn(ipoolhfx)
    ENDIF

    ! ==--------------------------------------------------------------==

    hfx_dist = hfxc5%hfx_distribution
    max_block_size = hfxc5%hfx_block_size
    IF (is_first_call) THEN
       IF (hfx_dist.EQ.hfx_dist_dynamic.AND.paral%io_parent) THEN
          WRITE(6,*)
          WRITE(6,'(A)') ' WARNING: In procedure '//procedureN&
               // ' HFX_DISTRIBUTION DYNAMIC is not '
          WRITE(6,'(A)') ' WARNING: supported the distribution is reset to'&
               // ' the default'
          WRITE(6,*)
       ENDIF
    ENDIF
    is_first_call = .FALSE.
    IF (max_block_size.LT.1) CALL stopgm( procedureN,&
         'HFX_BLOCK_SIZE should be greater than zero' ,& 
         __LINE__,__FILE__)

    ! get the max number of block pairs

    nbr_blk_state = CEILING(REAL(nstate,kind=real_8)/REAL(max_block_size,kind=real_8))
    nbr_blk_orb = CEILING(REAL(norb,kind=real_8)/REAL(max_block_size,kind=real_8))
    nblks = nbr_blk_orb * nbr_blk_state

    ! ==--------------------------------------------------------------==

    IF (.NOT.ALLOCATED(tick_blks)) THEN
       ALLOCATE(tick_blks(nblks),stat=ierr)! vw should be max(nstate_alpha,nstate_beta)
       IF (ierr.NE.0) CALL stopgm( procedureN, 'Allocation problem' ,& 
            __LINE__,__FILE__)
       tick_blks(:) = 1
    ENDIF

    ALLOCATE(vpotg(jhg),&
         sizes_orb(nbr_blk_orb),offsets_orb(nbr_blk_orb),&
         sizes_state(nbr_blk_state),offsets_state(nbr_blk_state),&
         psi_row(maxfftn,(max_block_size+1)/2),&
         psi_col(maxfftn,(max_block_size+1)/2),&
         C2_hfx(ncpw%ngw,norb),&
         vpotr(llr1),stat=ierr)
    IF (ierr.NE.0) CALL stopgm( procedureN, 'Allocation problem' ,& 
         __LINE__,__FILE__)
    ALLOCATE(local_h(norb,norb),stat=ierr)
    IF (ierr.NE.0) CALL stopgm( procedureN, 'Allocation problem' ,&
         __LINE__,__FILE__)
    ALLOCATE(sigma(norb),stat=ierr)
    IF (ierr.NE.0) CALL stopgm( procedureN, 'Allocation problem' ,&
         __LINE__,__FILE__)
    ALLOCATE(u(norb,norb),stat=ierr)
    IF (ierr.NE.0) CALL stopgm( procedureN, 'Allocation problem' ,&
         __LINE__,__FILE__)
    ALLOCATE(vt(norb,norb),stat=ierr)
    IF (ierr.NE.0) CALL stopgm( procedureN, 'Allocation problem' ,&
         __LINE__,__FILE__)
    ALLOCATE(c2_x(ncpw%ngw,norb),stat=ierr)
    IF (ierr.NE.0) CALL stopgm( procedureN, 'Allocation problem' ,&
         __LINE__,__FILE__)

    CALL zeroing( C2_hfx)!, ngw * norb )

    ! ==--------------------------------------------------------------==

    IF (hfx_resp_env%use_lin_lin.AND.hfx_resp_env%lin2_restart) THEN
       DO ia=1,norb
          norm=dotp(ncpw%ngw,cpsi(1,ia),cpsi(1,ia))
          CALL mp_sum(norm,parai%allgrp)
          IF (paral%io_parent) WRITE(222,*)' Norm of LR-orbital ',ia, ' >>>> ',norm
       ENDDO
    ENDIF

    ! HARTREE-FOCK

    CALL build_block_info(norb,nbr_blk_orb,max_block_size,&
         offsets_orb,sizes_orb)
    CALL build_block_info(nstate,nbr_blk_state,max_block_size,&
         offsets_state,sizes_state)
    CALL build_dist(dist,nblks,hfx_dist,tick_blks(:))

    CALL iterator_start(iter,dist)

    IF (hfx_resp_env%use_lin_lin .AND. &
         .NOT.hfx_resp_env%lin2_restart) THEN

       IF (paral%io_parent) &
            WRITE(6,*) 'LIN_LIN: apply the Xi vectors'

       pfl = 0.5_real_8
       IF (cntl%thybrid) pfl=pfl*func3%phfx
       pfl = pfl*sign
       CALL zeroing(C2_hfx)
       DO ia=1,norb
          DO ik=1,norb
             pfx=pfl !*f(ik) ! to check
             CALL hfxpb_lin2(pfx,cpsi(:,ia),hfx_lin2%xi(:,ik),C2_hfx(:,ia))
          ENDDO
       ENDDO

    ELSE
       !
       ! Conventional iterators over states
       !
       IF (hfx_resp_env%use_lin_lin .AND. paral%io_parent) &
            WRITE(6,*) 'LIN_LIN: restart the Xi vectors'

       IF (cntl%use_scaled_hfx) THEN
          CALL process_scaled_blocks()
       ELSE
          CALL process_blocks()
       ENDIF

    ENDIF
    ! 
    ! redistribute C2_hfx over the groups if needed
    ! 
    IF (parai%cp_nogrp.GT.1) THEN
       CALL tiset(procedureN//'_b',isub2)
       CALL cp_grp_redist(C2_hfx,ncpw%ngw,norb)
       CALL tihalt(procedureN//'_b',isub2)
    ENDIF
    !
    IF (hfx_resp_env%use_lin_lin.AND.hfx_resp_env%lin2_restart) THEN
       hfx_resp_env%lin2_restart=.FALSE.
       ! next time use a decomposition

       ! -- start test material
       DO ia=1,norb
          norm=dotp(ncpw%ngw,c2_hfx(1,ia),c2_hfx(1,ia))
          CALL mp_sum(norm,parai%allgrp)
       ENDDO
       DO ia=1,norb
          DO i=1,ncpw%ngw
             c2_x(i,ia)=c2_hfx(i,ia)/SQRT(norm)
          ENDDO
       ENDDO
       norm=0.0_real_8
       DO ia=1,norb
          norm=dotp(ncpw%ngw,cpsi(1,ia),cpsi(1,ia))
          CALL mp_sum(norm,parai%allgrp)
          IF (paral%io_parent) WRITE(6,'(A,i2,A,F12.5)')' Norm of LR ',ia, ' >>>> ',norm
       ENDDO

       CALL zeroing(hfx_lin2%m)
       DO i=1,norb
          DO j=1,norb !nstate
             hfx_lin2%m(i,j)=-dotp(ncpw%ngw,c0(1,i),c2_x(1,j))
          ENDDO
       ENDDO
       CALL mp_sum(hfx_lin2%m,SIZE(hfx_lin2%m),parai%allgrp)
       IF (paral%io_parent) WRITE(6,*) 'overlap of c2 with c0'
       IF (paral%io_parent)  THEN
          DO i=1,norb
             WRITE(6,'(100F12.5)') (hfx_lin2%m(i,j),j=1,norb)
          ENDDO
       ENDIF
       ! end test material

       CALL zeroing(hfx_lin2%m)
       DO i=1,norb !nstate
          DO j=1,norb !nstate
             hfx_lin2%m(i,j)=-dotp(ncpw%ngw,cpsi(1,i),c2_hfx(1,j)) !this is the correct one
             !hfx_lin2%m(i,j)=-dotp(ncpw%ngw,cpsi(1,i),c2_x(1,j))   !this is the
             !normalized one for testing
          ENDDO
       ENDDO
       CALL mp_sum(hfx_lin2%m,SIZE(hfx_lin2%m),parai%allgrp)

       IF (paral%io_parent) WRITE(6,*) 'overlap of c2 with cpsi'
       IF (paral%io_parent)  THEN
          DO i=1,norb
             WRITE(6,'(100F12.5)') (hfx_lin2%m(i,j),j=1,norb)
          ENDDO
       ENDIF

       IF (cholesky_decomp) THEN
          !----------------------------------------------------------------
          ! Cholesky decomposition -> L
          CALL dpotrf("U",norb,hfx_lin2%m,SIZE(hfx_lin2%m,1),ierr)
          IF (ierr.LT.0) THEN
             IF (paral%io_parent)WRITE(6,*) "WRONG ENTRY IN CHOL:", -ierr
             CALL stopgm(procedureN,'DPOTRF: WRONG ENTRY',&
                  __LINE__,__FILE__)
          ENDIF
          DO i=1, norb
             DO j=1, norb
                local_h(i,j) = 0._real_8
             ENDDO
          ENDDO
          DO i=1, norb
             local_h(i,i) = 1._real_8
          ENDDO
          ! CALCULATE INVERSE OF CHOLESKY FACTOR BY SOLVING
          !if (debug)  write(6,*) 'Test LIN2 h',local_h       
          DO i=1, norb
             CALL dtrsv("U", "N", "N", norb, hfx_lin2%m, SIZE(hfx_lin2%m,1),&
                  local_h(1,i), 1)
          ENDDO
          CALL zeroing(hfx_lin2%xi)
          CALL DGEMM('N','T',2*ncpw%ngw,norb,norb,1.0_real_8,hfx_lin2%c2,2*ncpw%ngw,&
               local_h,SIZE(local_h,1),0.0_real_8,hfx_lin2%xi,2*ncpw%ngw)
          !----------------------------------------------------------------
       ELSEIF (svd_decomp) THEN
          !----------------------------------------------------------------
          ALLOCATE(work(1),stat=ierr)
          IF (ierr.NE.0) CALL stopgm( procedureN, 'Allocation problem' ,&
               __LINE__,__FILE__)
          lwork=-1
          CALL dgesvd("A","A",norb,norb,hfx_lin2%m,SIZE(hfx_lin2%m,1),&
               sigma,u,SIZE(u,1),vt,SIZE(vt,1),work,lwork,ierr)
          lwork=work(1)
          DEALLOCATE(work,stat=ierr)
          IF (ierr.NE.0) CALL stopgm( procedureN, 'Deallocation problem' ,&
               __LINE__,__FILE__)
          ALLOCATE(work(lwork),stat=ierr)
          IF (ierr.NE.0) CALL stopgm( procedureN, 'Allocation problem' ,&
               __LINE__,__FILE__)
          ! A=U SIGMA V^T 
          ! SIGMA(1..norb)
          ! U (1..norb, 1...,norb)
          ! V the same as U
          CALL dgesvd("A","A",norb,norb,hfx_lin2%m,SIZE(hfx_lin2%m,1),&
               sigma,u,SIZE(u,1),vt,SIZE(vt,1),work,lwork,ierr)
          IF (ierr.NE.0) CALL stopgm( procedureN, 'Allocation problem' ,&
               __LINE__,__FILE__)
          DEALLOCATE(work,stat=ierr)
          IF (ierr.NE.0) CALL stopgm( procedureN, 'Deallocation problem' ,&
               __LINE__,__FILE__)
          IF (paral%io_parent) THEN
             DO i=1,norb
                WRITE(101,'(A,I2,A,E25.16,A)') 'Sigma(',i,')=', sigma(i), ';'
             ENDDO
          ENDIF
          max_sigma=MAXVAL(sigma)
          min_sigma=max_sigma*hfx_resp_env%sigma_eps
          icount=0
          DO ia=1,norb
             IF (sigma(ia).GT.min_sigma) icount=icount+1
          ENDDO
          !norb=icount
          IF (paral%io_parent) WRITE(*,*) 'use SVD with LIN2: icount, norb',icount, norb

          CALL zeroing(hfx_lin2%xi)
          CALL DGEMM('N','N',2*ncpw%ngw,icount,norb,1.0_real_8,hfx_lin2%c2,2*ncpw%ngw,&
               u,SIZE(u,1),0.0_real_8,hfx_lin2%xi,2*ncpw%ngw)
          DO ia=1,icount
             CALL dscal(2*ncpw%ngw,1.0_real_8/SQRT(sigma(ia)),hfx_lin2%xi(:,ia),1)
          ENDDO
          ! this can be optimized
          hfx_lin2%xi(:,icount+1:norb)=0.0_real_8
       ENDIF
       !----------------------------------------------------------------
    ENDIF
    ! 
    ! add up the hfx contribution to C2
    ! 
    CALL add_wfn(jgw,norb,zone,C2_hfx,ncpw%ngw,c2,ncpw%ngw)

    ! ==--------------------------------------------------------------==

    CALL release_dist(dist)

    DEALLOCATE(vpotg, vpotr, sizes_orb, offsets_orb,&
         sizes_state, offsets_state, psi_row, psi_col, C2_hfx,&
         stat=ierr)
    IF (ierr.NE.0) CALL stopgm( procedureN, 'Deallocation problem' ,& 
         __LINE__,__FILE__)

    IF (ALLOCATED(tick_blks)) THEN
       DEALLOCATE(tick_blks,stat=ierr)
       IF (ierr.NE.0) CALL stopgm( procedureN, 'Deallocation problem',& 
            __LINE__,__FILE__)
    ENDIF

    DEALLOCATE(local_h,stat=ierr)
    IF (ierr.NE.0) CALL stopgm( procedureN, 'Deallocation problem' ,&
         __LINE__,__FILE__)
    DEALLOCATE(sigma,stat=ierr)
    IF (ierr.NE.0) CALL stopgm( procedureN, 'Deallocation problem' ,&
         __LINE__,__FILE__)
    DEALLOCATE(u,stat=ierr)
    IF (ierr.NE.0) CALL stopgm( procedureN, 'Deallocation problem' ,&
         __LINE__,__FILE__)
    DEALLOCATE(vt,stat=ierr)
    IF (ierr.NE.0) CALL stopgm( procedureN, 'Deallocation problem' ,&
         __LINE__,__FILE__)
    DEALLOCATE(c2_x,stat=ierr)
    IF (ierr.NE.0) CALL stopgm( procedureN, 'Deallocation problem' ,&
         __LINE__,__FILE__)

    ! ==--------------------------------------------------------------==

    CALL tihalt(procedureN,isub)


    CONTAINS

    ! ==================================================================
    SUBROUTINE process_blocks()
       ! ==--------------------------------------------------------------==
       ! Conventional loop for full HFX, no storage of real-space quantities

       DO WHILE (iterator_blocks(iter))
          CALL iterator_next_block(iter,blk)
          CALL ind2sub_rect(nbr_blk_state,blk,row_blk,col_blk)
          row_offset = offsets_state(row_blk)
          col_offset = offsets_orb(col_blk)
          row_size = sizes_state(row_blk)
          col_size = sizes_orb(col_blk)

          CALL block_invfft_old(c0  ,psi_row,row_offset,row_size)
          CALL block_invfft_old(cpsi,psi_col,col_offset,col_size)


          c = 0
          DO ia = col_offset,col_offset+col_size-1,2
             c = c + 1
             ia1=ia
             ia2=ia+1
             IF (ia2.GT.norb) ia2=0

             r = 0
             DO ib = row_offset,row_offset+row_size-1,2
                r = r + 1
                ib1=ib
                ib2=ib+1
                IF (f(ib1).LT.1.e-6_real_8) ib1=0
                IF (ib2.GT.nstate) THEN
                   ib2=0
                ELSEIF (f(ib2).LT.1.e-6_real_8) THEN
                   ib2=0
                ENDIF

                pfl = 0.5_real_8
                IF (cntl%thybrid) pfl=pfl*func3%phfx
                pfl = pfl*sign

                IF (ib1.NE.0 .OR. ib2.NE.0) THEN
                   IF (ib1.NE.0) THEN
                      pfx=pfl*f(ib1)
                      CALL hfxpb(pfx,psi_col(:,c),.TRUE.,&
                           psi_row(:,r),.TRUE.,vpotg,vpotr,psic,&
                           C2_hfx(:,ia1))
                   ENDIF
                   IF (ib2.NE.0) THEN
                      pfx=pfl*f(ib2)
                      CALL hfxpb(pfx,psi_col(:,c),.TRUE.,&
                           psi_row(:,r),.FALSE.,vpotg,vpotr,psic,&
                           C2_hfx(:,ia1))
                   ENDIF
                   IF (ia2.GT.0) THEN
                      IF (ib1.NE.0) THEN
                         pfx=pfl*f(ib1)
                         CALL hfxpb(pfx,psi_col(:,c),.FALSE.,&
                              psi_row(:,r),.TRUE.,vpotg,vpotr,psic,&
                              C2_hfx(:,ia2))
                      ENDIF
                      IF (ib2.NE.0) THEN
                         pfx=pfl*f(ib2)
                         CALL hfxpb(pfx,psi_col(:,c),.FALSE.,&
                              psi_row(:,r),.FALSE.,vpotg,vpotr,psic,&
                              C2_hfx(:,ia2))
                      ENDIF
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDDO ! iter

       ! ==--------------------------------------------------------------==
    END SUBROUTINE process_blocks
    ! ==================================================================
    SUBROUTINE process_scaled_blocks()
       ! ==--------------------------------------------------------------==
       ! Loop for scaled exchange, including storage.

       COMPLEX(real_8), DIMENSION(:,:), &
                        ALLOCATABLE :: C2_in_real_space
 
       CHARACTER(len=*), PARAMETER  :: subprocedureN = procedureN//'_process_scaled_blocks'
 
       INTEGER         :: iorb, ig, ir
       INTEGER         :: blk_orb, blk_resid
       
       COMPLEX(real_8) :: fp, fm

       ALLOCATE(C2_in_real_space(maxfftn,norb),&
                stat=ierr)
       IF (ierr /= 0) CALL stopgm(subprocedureN,'Allocation problem',&
                                 __LINE__,__FILE__)
       CALL zeroing(C2_in_real_space)

       CALL scex%start_density_scaling()

       DO WHILE (iterator_blocks(iter))
          CALL iterator_next_block(iter,blk)
          CALL ind2sub_rect(nbr_blk_state,blk,row_blk,col_blk)
          row_offset = offsets_state(row_blk)
          col_offset = offsets_orb(col_blk)
          row_size = sizes_state(row_blk)
          col_size = sizes_orb(col_blk)

          !
          ! See comments in pw_hfx.mod.F90
          !
          CALL block_invfft_old(c0  ,psi_row,row_offset,row_size)
          CALL block_invfft_old(cpsi,psi_col,col_offset,col_size)
          
          !
          ! Scaled orbitals
          !
          CALL setfftn(scex_ID_scaled)
          c = 0
          DO ia = col_offset,col_offset+col_size-1,2
             c = c + 1
             ia1=ia
             ia2=ia+1
             IF (ia2.GT.norb) ia2=0

             r = 0
             DO ib = row_offset,row_offset+row_size-1,2
                r = r + 1
                ib1=ib
                ib2=ib+1
                IF (f(ib1).LT.1.e-6_real_8) ib1=0
                IF (ib2.GT.nstate) THEN
                   ib2=0
                ELSEIF (f(ib2).LT.1.e-6_real_8) THEN
                   ib2=0
                ENDIF

                pfl = 0.5_real_8
                IF (cntl%thybrid) pfl=pfl*func3%phfx
                pfl = pfl*sign

                IF (ib1.NE.0 .OR. ib2.NE.0) THEN
                   IF (ib1.NE.0) THEN
                      pfx=pfl*f(ib1)
                      CALL hfxpb(pfx,psi_col(:,c),.TRUE.,&
                           psi_row(:,r),.TRUE.,vpotg,vpotr,psic,&
                           C2_in_real_space(:,ia1))
                   ENDIF
                   IF (ib2.NE.0) THEN
                      pfx=pfl*f(ib2)
                      CALL hfxpb(pfx,psi_col(:,c),.TRUE.,&
                           psi_row(:,r),.FALSE.,vpotg,vpotr,psic,&
                           C2_in_real_space(:,ia1))
                   ENDIF
                   IF (ia2.GT.0) THEN
                      IF (ib1.NE.0) THEN
                         pfx=pfl*f(ib1)
                         CALL hfxpb(pfx,psi_col(:,c),.FALSE.,&
                              psi_row(:,r),.TRUE.,vpotg,vpotr,psic,&
                              C2_in_real_space(:,ia2))
                      ENDIF
                      IF (ib2.NE.0) THEN
                         pfx=pfl*f(ib2)
                         CALL hfxpb(pfx,psi_col(:,c),.FALSE.,&
                              psi_row(:,r),.FALSE.,vpotg,vpotr,psic,&
                              C2_in_real_space(:,ia2))
                      ENDIF
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
          !
          ! Conventional orbitals
          !
          CALL setfftn(scex_ID_parent)
          !
       ENDDO ! iter

       CALL scex%switch_density_scaling()

       blk_orb   = norb / 2
       blk_resid = norb - blk_orb*2

       DO iorb=1,blk_orb
          DO ir=1,scex%llr1
            C2_in_real_space(ir,iorb) = REAL(C2_in_real_space(ir,iorb))&
                                        + uimag*REAL(C2_in_real_space(ir,iorb+blk_orb))
          ENDDO
          CALL scex%undo_density_scaling(C2_in_real_space(:,iorb))
          CALL fwfftn(C2_in_real_space(:,iorb),.TRUE.,parai%allgrp)
          !$omp parallel do private(IG,FP,FM)
          DO ig=1,jgw
             fp = C2_in_real_space(nzfs(ig),iorb) + C2_in_real_space(inzs(ig),iorb)
             fm = C2_in_real_space(nzfs(ig),iorb) - C2_in_real_space(inzs(ig),iorb)
             C2_hfx(ig,iorb)         = 0.5_real_8*CMPLX(REAL(fp),AIMAG(fm),kind=real_8)
             C2_hfx(ig,iorb+blk_orb) = 0.5_real_8*CMPLX(AIMAG(fp),-REAL(fm),kind=real_8)
          ENDDO
 
       ENDDO
       IF (blk_resid > 0) THEN
          CALL scex%undo_density_scaling(C2_in_real_space(:,norb))
          CALL fwfftn(C2_in_real_space(:,norb),.TRUE.,parai%allgrp)
          !$omp parallel do private(IG,FP,FM)
          DO ig = 1,jgw
             fp = C2_in_real_space(nzfs(ig),norb) + C2_in_real_space(inzs(ig),norb)
             fm = C2_in_real_space(nzfs(ig),norb) - C2_in_real_space(inzs(ig),norb)
             C2_hfx(ig,norb) = 0.5_real_8*CMPLX(REAL(fp),AIMAG(fm),kind=real_8)
          ENDDO
       ENDIF

       CALL scex%annihilate_density_scaling()

       DEALLOCATE(C2_in_real_space,stat=ierr)
       IF (ierr /= 0) CALL stopgm(subprocedureN,'Deallocation problem',&
                                  __LINE__,__FILE__)

       ! ==--------------------------------------------------------------==
    END SUBROUTINE process_scaled_blocks
    ! ==================================================================
  END SUBROUTINE hfxpsi_new
  ! ==================================================================
  SUBROUTINE hfxpb(pf,psia,a_stored_in_real,psib,b_stored_in_real,&
       vpotg,vpotr,psic,c2a)
    ! ==================================================================
    ! Simplifications as in pw_hfx; density scaling
    !                                 12.02.2019 M.P.Bircher @ LCBC/EPFL
    ! ==--------------------------------------------------------------==

    REAL(real_8)                             :: pf
    COMPLEX(real_8)                          :: psia(:)
    LOGICAL                                  :: a_stored_in_real
    COMPLEX(real_8)                          :: psib(:)
    LOGICAL                                  :: b_stored_in_real
    COMPLEX(real_8)                          :: vpotg(:)
    REAL(real_8)                             :: vpotr(:)
    COMPLEX(real_8)                          :: psic(:), c2a(:)

    INTEGER                                  :: ig, ir

    CALL hfxresp_get_pair_density(psia,psib,psic,a_stored_in_real,b_stored_in_real)
    CALL fwfftn(psic,.FALSE.,parai%allgrp)
    CALL hfxresp_get_coulomb(psic,vpotg,pf)
    CALL hfxresp_set_vpotg(psic,vpotg)
    CALL invfftn(psic,.FALSE.,parai%allgrp)
    CALL hfxresp_set_vpotr(psic,vpotr)
    CALL hfxresp_get_potential(psib,psic,vpotr,b_stored_in_real)

    IF (cntl%use_scaled_hfx) THEN
       CALL hfxresp_get_c2_real_space(psic,c2a)
    ELSE
       CALL fwfftn(psic,.TRUE.,parai%allgrp)
       CALL hfxpb_get_c2(psic,c2a)
    ENDIF

    ! ==--------------------------------------------------------------==
  END SUBROUTINE hfxpb
  ! ==================================================================
  SUBROUTINE hfxpb_lin2(pf,psia,psib,c2a)
    ! ==================================================================

    ! ==--------------------------------------------------------------==

    REAL(real_8)                             :: pf
    COMPLEX(real_8)                          :: psia(:), psib(:), c2a(:)

    INTEGER                                  :: ig
    REAL(real_8)                             :: dotp, z_ka

    z_ka=dotp(ncpw%ngw,psia(1),psib(1))
    CALL mp_sum(z_ka,parai%allgrp)
    DO ig=1,ncpw%ngw
       !c2a(ig)=c2a(ig)+(-1._real_8)*pf*z_ka*psib(ig)
       c2a(ig)=c2a(ig)+pf*z_ka*psib(ig)
    ENDDO

    ! ==--------------------------------------------------------------==
  END SUBROUTINE hfxpb_lin2
  ! ==================================================================

  ! ==================================================================
  ! RPA
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE hfxrpa_new(c0,c1,c2,psia,nstate,tcis)
    ! ==================================================================
    ! Hartree-Fock exchange contribution to RPA
    ! Updated with coordinate-scaling and refactored on 12.02.2019
    !                                           M.P. Bircher @ LCBC/EPFL
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(:,:), c1(:,:), c2(:,:), &
                                                psia(:)
    INTEGER                                  :: nstate
    LOGICAL                                  :: tcis

    CHARACTER(*), PARAMETER                  :: procedureN = 'hfxrpa_new'
    COMPLEX(real_8), PARAMETER               :: zone = (1.0_real_8,0.0_real_8)

    COMPLEX(real_8), ALLOCATABLE, &
      DIMENSION(:)                           :: psic, vpotg
    COMPLEX(real_8), ALLOCATABLE, &
      DIMENSION(:, :)                        :: C2_hfx, psi_row
    INTEGER :: blk, col_blk, col_offset, col_size, hfx_dist, ia, ib, ib1, &
      ib2, ierr, ig, ispin, isub, isub2, max_block_size, max_nblks, nblks, &
      nbr_blk, nspins, nst, nstates(2), r, row_blk, row_offset, row_size, &
      st_offst
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: offsets, sizes
    INTEGER, ALLOCATABLE, DIMENSION(:), SAVE :: tick_blks
    LOGICAL, SAVE                            :: is_first_call = .TRUE.
    REAL(real_8)                             :: pfl
    REAL(real_8), ALLOCATABLE, DIMENSION(:)  :: vpotr
    TYPE(dist_t)                             :: dist
    TYPE(iterator_t)                         :: iter

    IF (func1%mhfx.EQ.0) RETURN

    ! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)

    IF (func1%mhfx.EQ.2) CALL stopgm(procedureN,&
         'HARTREE METHOD NOT POSSIBLE',& 
         __LINE__,__FILE__)
    IF (tkpts%tkpnt) CALL stopgm(procedureN,'K-POINTS NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    IF (lspin2%tlse) CALL stopgm(procedureN,'NO LSE POSSIBLE',& 
         __LINE__,__FILE__)
    IF (pslo_com%tivan) CALL stopgm(procedureN,'NO VDB PP POSSIBLE',& 
         __LINE__,__FILE__)
    IF (group%nogrp.GT.1) CALL stopgm(procedureN,&
         'OLD TASK GROUPS NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    ! 
    IF (cntl%use_scaled_hfx) THEN
       IF (.NOT. scex%init) CALL stopgm(procedureN,'ScEX requested, but not initialised',&
                                        __LINE__,__FILE__)
       CALL setfftn(scex_ID_parent)
    ELSE
       CALL setfftn(ipoolhfx)
    ENDIF

    ! ==--------------------------------------------------------------==

    nstates = 0
    IF (cntl%tlsd) THEN
       nspins = 2
       nstates(1) = spin_mod%nsup
       nstates(2) = spin_mod%nsdown
       pfl = 2.0_real_8
    ELSE
       nspins = 1
       nstates(1) = nstate
       pfl = 1.0_real_8
    ENDIF

    IF (.NOT.tcis) pfl=2._real_8*pfl
    IF (cntl%thybrid) pfl=pfl*func3%phfx

    ! ==--------------------------------------------------------------==

    hfx_dist = hfxc5%hfx_distribution
    max_block_size = hfxc5%hfx_block_size
    IF (is_first_call) THEN
       IF (hfx_dist.EQ.hfx_dist_dynamic.AND.paral%io_parent) THEN
          WRITE(6,*)
          WRITE(6,'(A)') ' WARNING: In procedure '//procedureN&
               // ' HFX_DISTRIBUTION DYNAMIC is not supported '
          WRITE(6,'(A)') ' WARNING: the distribution is reset to'&
               // 'the default'
          WRITE(6,*)
       ENDIF
    ENDIF
    is_first_call = .FALSE.
    IF (max_block_size.LT.1) CALL stopgm( procedureN,&
         'HFX_BLOCK_SIZE should be greater than zero' ,& 
         __LINE__,__FILE__)

    ! get the max number of block pairs

    nbr_blk = CEILING(REAL(MAXVAL(nstates),kind=real_8)/REAL(max_block_size,kind=real_8))
    max_nblks = nbr_blk**2

    ! ==--------------------------------------------------------------==

    IF (.NOT.ALLOCATED(tick_blks)) THEN
       ALLOCATE(tick_blks(max_nblks),stat=ierr)
       IF (ierr.NE.0) CALL stopgm( procedureN, 'Allocation problem' ,& 
            __LINE__,__FILE__)
       tick_blks(:) = 1
    ENDIF

    ! ==--------------------------------------------------------------==

    ALLOCATE(vpotg(jhg),&
         psic(maxfftn),&
         sizes(nbr_blk),offsets(nbr_blk),&
         psi_row(maxfftn,(max_block_size+1)/2),&
         C2_hfx(ncpw%ngw,nstate),&
         vpotr(llr1),stat=ierr)
    IF (ierr.NE.0) CALL stopgm( procedureN, 'Allocation problem' ,& 
         __LINE__,__FILE__)

    CALL zeroing( C2_hfx)!, ngw * nstate )

    ! ==--------------------------------------------------------------==

    IF (hfx_resp_env%write_v_ab) THEN
       OPEN(2000+parai%mepos,form="unformatted") !,access="sequential")
       REWIND(2000+parai%mepos)
    ENDIF

    !
    ! Main loop
    !
    IF (cntl%use_scaled_hfx) THEN
       CALL process_scaled_blocks() 
    ELSE
       CALL process_blocks()
    ENDIF


    IF (hfx_resp_env%write_v_ab)  CLOSE(2000+parai%mepos)
    ! 
    ! redistribute C2_hfx over the groups if needed
    ! 
    IF (parai%cp_nogrp.GT.1) THEN
       CALL tiset(procedureN//'_b',isub2)
       CALL cp_grp_redist(C2_hfx,ncpw%ngw,nstate)
       CALL tihalt(procedureN//'_b',isub2)
    ENDIF

    ! 
    ! add up the hfx contribution to C2
    ! 
    CALL add_wfn(jgw,nstate,zone,C2_hfx,ncpw%ngw,c2,ncpw%ngw)

    ! ==--------------------------------------------------------------==

    CALL release_dist(dist)

    DEALLOCATE(vpotg,psic,sizes,offsets,psi_row,C2_hfx,vpotr,&
         stat=ierr)
    IF (ierr.NE.0) CALL stopgm( procedureN, 'Deallocation problem' ,& 
         __LINE__,__FILE__)

    IF (ALLOCATED(tick_blks)) THEN
       DEALLOCATE(tick_blks,stat=ierr)
       IF (ierr.NE.0) CALL stopgm( procedureN, 'Deallocation problem',& 
            __LINE__,__FILE__)
    ENDIF

    ! ==--------------------------------------------------------------==

    CALL tihalt(procedureN,isub)


    CONTAINS

    ! ==================================================================
    SUBROUTINE process_blocks()
       ! ==--------------------------------------------------------------==
       ! Conventional loop for full HFX, no storage of real-space quantities

       DO ispin = 1,nspins

          nst = nstates(ispin)
          st_offst = 0
          IF (ispin.EQ.2) st_offst = nstates(1)

          ! get the number of block pairs

          nbr_blk = CEILING(REAL(nst,kind=real_8)/REAL(max_block_size,kind=real_8))
          nblks = nbr_blk**2

          CALL build_block_info(nst,nbr_blk,max_block_size,&
               offsets,sizes)
          CALL build_dist(dist,nblks,hfx_dist,tick_blks(:))

          CALL iterator_start(iter,dist)
          DO WHILE (iterator_blocks(iter))
             CALL iterator_next_block(iter,blk)
             CALL ind2sub_rect(nbr_blk,blk,row_blk,col_blk)
             row_offset = offsets(row_blk) + st_offst
             col_offset = offsets(col_blk) + st_offst
             row_size = sizes(row_blk); col_size = sizes(col_blk)

             CALL block_invfft_old(c0,psi_row,row_offset,row_size)

             DO ia = col_offset,col_offset+col_size-1
                ! >>>>
                CALL zeroing(psia)!,maxfftn)
                !ocl novrec
                DO ig=1,jgw
                   psia(nzfs(ig))=c1(ig,ia)+uimag*c0(ig,ia)
                   psia(inzs(ig))=CONJG(c1(ig,ia))+&
                        uimag*CONJG(c0(ig,ia))
                ENDDO
                IF (geq0) psia(nzfs(1))=c1(1,ia)+uimag*c0(1,ia)

                ! Transform the wavefunctions to real space
                CALL invfftn(psia,.TRUE.,parai%allgrp)
                ! <<<< need a block_invfft for this operation

                r = 0
                DO ib = row_offset,row_offset+row_size-1,2
                   r = r + 1

                   ib1=ib
                   ib2=ib+1
                   IF (ib2.GT.nst+st_offst) ib2=0

                   CALL hfxrpav(pfl,ia,ib1,psia,psi_row(:,r),1,vpotg,vpotr,&
                        psic,C2_hfx(:,ib1))

                   IF (ib2.NE.0) THEN
                      CALL hfxrpav(pfl,ia,ib2,psia,psi_row(:,r),2,vpotg,vpotr,&
                           psic,C2_hfx(:,ib2))
                   ENDIF

                ENDDO! IB
             ENDDO! IA
          ENDDO! iter
       ENDDO ! ispin


       ! ==--------------------------------------------------------------==
    END SUBROUTINE process_blocks
    ! ==================================================================
    SUBROUTINE process_scaled_blocks()
       ! ==--------------------------------------------------------------==
       ! Loop for scaled exchange, including storage.

       COMPLEX(real_8), DIMENSION(:,:), &
                        ALLOCATABLE :: C2_in_real_space
 
       CHARACTER(len=*), PARAMETER  :: subprocedureN = procedureN//'_process_scaled_blocks'
 
       INTEGER         :: istate, ig, ir
       INTEGER         :: blk_state, blk_resid
       
       COMPLEX(real_8) :: fp, fm

       ALLOCATE(C2_in_real_space(maxfftn,nstate),&
                stat=ierr)
       IF (ierr /= 0) CALL stopgm(subprocedureN,'Allocation problem',&
                                 __LINE__,__FILE__)
       CALL zeroing(C2_in_real_space)

       CALL scex%start_density_scaling()

       DO ispin = 1,nspins

          nst = nstates(ispin)
          st_offst = 0
          IF (ispin.EQ.2) st_offst = nstates(1)

          ! get the number of block pairs

          nbr_blk = CEILING(REAL(nst,kind=real_8)/REAL(max_block_size,kind=real_8))
          nblks = nbr_blk**2

          CALL build_block_info(nst,nbr_blk,max_block_size,&
               offsets,sizes)
          CALL build_dist(dist,nblks,hfx_dist,tick_blks(:))

          CALL iterator_start(iter,dist)
          DO WHILE (iterator_blocks(iter))
             CALL iterator_next_block(iter,blk)
             CALL ind2sub_rect(nbr_blk,blk,row_blk,col_blk)
             row_offset = offsets(row_blk) + st_offst
             col_offset = offsets(col_blk) + st_offst
             row_size = sizes(row_blk); col_size = sizes(col_blk)
        
             CALL block_invfft_old(c0,psi_row,row_offset,row_size)

             DO ia = col_offset,col_offset+col_size-1
                ! >>>>
                CALL zeroing(psia)!,maxfftn)
                !ocl novrec
                DO ig=1,jgw
                   psia(nzfs(ig))=c1(ig,ia)+uimag*c0(ig,ia)
                   psia(inzs(ig))=CONJG(c1(ig,ia))+&
                        uimag*CONJG(c0(ig,ia))
                ENDDO
                IF (geq0) psia(nzfs(1))=c1(1,ia)+uimag*c0(1,ia)

                ! Transform the wavefunctions to real space
                CALL invfftn(psia,.TRUE.,parai%allgrp)
                ! <<<< need a block_invfft for this operation

                CALL scex%do_density_scaling(psia)

                !
                ! Scaled orbitals
                !
                CALL setfftn(scex_ID_scaled)
                r = 0
                DO ib = row_offset,row_offset+row_size-1,2
                   r = r + 1

                   ib1=ib
                   ib2=ib+1
                   IF (ib2.GT.nst+st_offst) ib2=0

                   CALL hfxrpav(pfl,ia,ib1,psia,psi_row(:,r),1,vpotg,vpotr,&
                        psic,C2_in_real_space(:,ib1))

                   IF (ib2.NE.0) THEN
                      CALL hfxrpav(pfl,ia,ib2,psia,psi_row(:,r),2,vpotg,vpotr,&
                           psic,C2_in_real_space(:,ib2))
                   ENDIF

                ENDDO! IB
                !
                ! Conventional orbitals
                !
                CALL setfftn(scex_ID_parent)
                !
             ENDDO! IA
          ENDDO! iter
       ENDDO ! ispin

       CALL scex%switch_density_scaling()

       blk_state = nstate / 2
       blk_resid = nstate - blk_state*2

       DO istate=1,blk_state
          DO ir=1,scex%llr1
            C2_in_real_space(ir,istate) = REAL(C2_in_real_space(ir,istate))&
                                        + uimag*REAL(C2_in_real_space(ir,istate+blk_state))
          ENDDO
          CALL scex%undo_density_scaling(C2_in_real_space(:,istate))
          CALL fwfftn(C2_in_real_space(:,istate),.TRUE.,parai%allgrp)
          !$omp parallel do private(IG,FP,FM)
          DO ig=1,jgw
             fp = C2_in_real_space(nzfs(ig),istate) + C2_in_real_space(inzs(ig),istate)
             fm = C2_in_real_space(nzfs(ig),istate) - C2_in_real_space(inzs(ig),istate)
             C2_hfx(ig,istate)           = -0.5_real_8*CMPLX(REAL(fp),AIMAG(fm),kind=real_8)
             C2_hfx(ig,istate+blk_state) = -0.5_real_8*CMPLX(AIMAG(fp),-REAL(fm),kind=real_8)
          ENDDO
 
       ENDDO
       IF (blk_resid > 0) THEN
          CALL scex%undo_density_scaling(C2_in_real_space(:,nstate))
          CALL fwfftn(C2_in_real_space(:,nstate),.TRUE.,parai%allgrp)
          !$omp parallel do private(IG,FP,FM)
          DO ig = 1,jgw
             fp = C2_in_real_space(nzfs(ig),nstate) + C2_in_real_space(inzs(ig),nstate)
             fm = C2_in_real_space(nzfs(ig),nstate) - C2_in_real_space(inzs(ig),nstate)
             C2_hfx(ig,nstate) = -0.5_real_8*CMPLX(REAL(fp),AIMAG(fm),kind=real_8)
          ENDDO
       ENDIF

       CALL scex%annihilate_density_scaling()

       DEALLOCATE(C2_in_real_space,stat=ierr)
       IF (ierr /= 0) CALL stopgm(subprocedureN,'Deallocation problem',&
                                  __LINE__,__FILE__)

       ! ==--------------------------------------------------------------==
    END SUBROUTINE process_scaled_blocks
    ! ==================================================================
  END SUBROUTINE hfxrpa_new
  ! ==================================================================
  SUBROUTINE hfxrpa_incore(c0,c1,c2,psia,nstate,tcis)
    ! ==================================================================
    ! == Hartree-Fock exchange contribution to RPA                    ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(:,:), c1(:,:), c2(:,:), &
                                                psia(:)
    INTEGER                                  :: nstate
    LOGICAL                                  :: tcis

    CHARACTER(*), PARAMETER                  :: procedureN = 'hfxrpa_incore'
    COMPLEX(real_8), PARAMETER               :: zone = (1.0_real_8,0.0_real_8)

    COMPLEX(real_8), ALLOCATABLE, &
      DIMENSION(:)                           :: psic, vpotg
    COMPLEX(real_8), ALLOCATABLE, &
      DIMENSION(:, :)                        :: C2_hfx, psi_row
    INTEGER :: blk, col_blk, col_offset, col_size, hfx_dist, ia, ib, ib1, &
      ib2, ierr, ig, ispin, isub, isub2, max_block_size, max_nblks, nblks, &
      nbr_blk, nspins, nst, nstates(2), r, row_blk, row_offset, row_size, &
      st_offst
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: offsets, sizes
    INTEGER, ALLOCATABLE, DIMENSION(:), SAVE :: tick_blks
    LOGICAL, SAVE                            :: is_first_call = .TRUE.
    REAL(real_8)                             :: pfl
    REAL(real_8), ALLOCATABLE, DIMENSION(:)  :: vpotr
    TYPE(dist_t)                             :: dist
    TYPE(iterator_t)                         :: iter

    IF (func1%mhfx.EQ.0) RETURN

    ! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)

    IF (func1%mhfx.EQ.2) CALL stopgm(procedureN,&
         'HARTREE METHOD NOT POSSIBLE',&
         __LINE__,__FILE__)
    IF (tkpts%tkpnt) CALL stopgm(procedureN,'K-POINTS NOT IMPLEMENTED',&
         __LINE__,__FILE__)
    IF (lspin2%tlse) CALL stopgm(procedureN,'NO LSE POSSIBLE',&
         __LINE__,__FILE__)
    IF (pslo_com%tivan) CALL stopgm(procedureN,'NO VDB PP POSSIBLE',&
         __LINE__,__FILE__)
    IF (group%nogrp.GT.1) CALL stopgm(procedureN,&
         'OLD TASK GROUPS NOT IMPLEMENTED',&
         __LINE__,__FILE__)
    ! 
    CALL setfftn(ipoolhfx)

    ! ==--------------------------------------------------------------==

    nstates = 0
    IF (cntl%tlsd) THEN
       nspins = 2
       nstates(1) = spin_mod%nsup
       nstates(2) = spin_mod%nsdown
       pfl = 2.0_real_8
    ELSE
       nspins = 1
       nstates(1) = nstate
       pfl = 1.0_real_8
    ENDIF

    IF (.NOT.tcis) pfl=2._real_8*pfl
    IF (cntl%thybrid) pfl=pfl*func3%phfx

    ! ==--------------------------------------------------------------==

    hfx_dist = hfxc5%hfx_distribution
    max_block_size = hfxc5%hfx_block_size
    IF (is_first_call) THEN
       IF (hfx_dist.EQ.hfx_dist_dynamic.AND.paral%io_parent) THEN
          WRITE(6,*)
          WRITE(6,'(A)') ' WARNING: In procedure '//procedureN&
               // ' HFX_DISTRIBUTION DYNAMIC is not supported '
          WRITE(6,'(A)') ' WARNING: the distribution is reset to'&
               // 'the default'
          WRITE(6,*)
       ENDIF
    ENDIF
    is_first_call = .FALSE.
    IF (max_block_size.LT.1) CALL stopgm( procedureN,&
         'HFX_BLOCK_SIZE should be greater than zero' ,&
         __LINE__,__FILE__)

    ! get the max number of block pairs

    nbr_blk = CEILING(REAL(MAXVAL(nstates),kind=real_8)/REAL(max_block_size,kind=real_8))
    max_nblks = nbr_blk**2

    ! ==--------------------------------------------------------------==

    IF (.NOT.ALLOCATED(tick_blks)) THEN
       ALLOCATE(tick_blks(max_nblks),stat=ierr)
       IF (ierr.NE.0) CALL stopgm( procedureN, 'Allocation problem' ,&
            __LINE__,__FILE__)
       tick_blks(:) = 1
    ENDIF

    ! ==--------------------------------------------------------------==

    ALLOCATE(vpotg(jhg),&
         psic(maxfftn),&
         sizes(nbr_blk),offsets(nbr_blk),&
         psi_row(maxfftn,(max_block_size+1)/2),&
         C2_hfx(ncpw%ngw,nstate),&
         vpotr(llr1),stat=ierr)
    IF (ierr.NE.0) CALL stopgm( procedureN, 'Allocation problem' ,&
         __LINE__,__FILE__)

    CALL zeroing( C2_hfx)!, ngw * nstate )

    ! ==--------------------------------------------------------------==

    IF (hfx_resp_env%write_v_ab) THEN
       OPEN(2000+parai%mepos,form="unformatted") !,access="sequential")
       REWIND(2000+parai%mepos)
    ENDIF

    DO ispin = 1,nspins

       nst = nstates(ispin)
       st_offst = 0
       IF (ispin.EQ.2) st_offst = nstates(1)

       ! get the number of block pairs

       nbr_blk = CEILING(REAL(nst,kind=real_8)/REAL(max_block_size,kind=real_8))
       nblks = nbr_blk**2

       CALL build_block_info(nst,nbr_blk,max_block_size,&
            offsets,sizes)
       CALL build_dist(dist,nblks,hfx_dist,tick_blks(:))

       CALL iterator_start(iter,dist)
       DO WHILE (iterator_blocks(iter))
          CALL iterator_next_block(iter,blk)
          CALL ind2sub_rect(nbr_blk,blk,row_blk,col_blk)
          row_offset = offsets(row_blk) + st_offst
          col_offset = offsets(col_blk) + st_offst
          row_size = sizes(row_blk); col_size = sizes(col_blk)

          DO ia = col_offset,col_offset+col_size-1
             ! >>>>
             CALL zeroing(psia)!,maxfftn)
             !ocl novrec
             DO ig=1,jgw
                psia(nzfs(ig))=c1(ig,ia)+uimag*c0(ig,ia)
                psia(inzs(ig))=CONJG(c1(ig,ia))+&
                     uimag*CONJG(c0(ig,ia))
             ENDDO
             IF (geq0) psia(nzfs(1))=c1(1,ia)+uimag*c0(1,ia)

             ! Transform the wavefunctions to real space
             CALL invfftn(psia,.TRUE.,parai%allgrp)
             ! <<<< need a block_invfft for this operation

             r = 0
             DO ib = row_offset,row_offset+row_size-1,2
                r = r + 1

                ib1=ib
                ib2=ib+1
                IF (ib2.GT.nst+st_offst) ib2=0

                CALL hfxrpav_incore(pfl,ia,ib1,psia,psi_row(:,r),1,vpotg,vpotr,&
                     psic,C2_hfx(:,ib1))

                IF (ib2.NE.0) THEN
                   CALL hfxrpav_incore(pfl,ia,ib2,psia,psi_row(:,r),2,vpotg,vpotr,&
                        psic,C2_hfx(:,ib2))
                ENDIF

             ENDDO! IB
          ENDDO! IA
       ENDDO! iter
    ENDDO ! ispin

    IF (hfx_resp_env%write_v_ab) CLOSE(2000+parai%mepos)
    ! 
    ! redistribute C2_hfx over the groups if needed
    ! 
    IF (parai%cp_nogrp.GT.1) THEN
       CALL tiset(procedureN//'_b',isub2)
       CALL cp_grp_redist(C2_hfx,ncpw%ngw,nstate)
       CALL tihalt(procedureN//'_b',isub2)
    ENDIF

    ! 
    ! add up the hfx contribution to C2
    ! 
    CALL add_wfn(jgw,nstate,zone,C2_hfx,ncpw%ngw,c2,ncpw%ngw)

    ! ==--------------------------------------------------------------==

    CALL release_dist(dist)

    DEALLOCATE(vpotg,psic,sizes,offsets,psi_row,C2_hfx,vpotr,&
         stat=ierr)
    IF (ierr.NE.0) CALL stopgm( procedureN, 'Deallocation problem' ,&
         __LINE__,__FILE__)

    IF (ALLOCATED(tick_blks)) THEN
       DEALLOCATE(tick_blks,stat=ierr)
       IF (ierr.NE.0) CALL stopgm( procedureN, 'Deallocation problem',&
            __LINE__,__FILE__)
    ENDIF

    ! ==--------------------------------------------------------------==

    CALL tihalt(procedureN,isub)

    ! ==--------------------------------------------------------------==

    ! ==================================================================
  END SUBROUTINE hfxrpa_incore
  ! ==================================================================
  SUBROUTINE hfxrpav(pf,ia,ib,psia,psib,iran,vpotg,vpotr,psic,c2a)
    ! ==================================================================
    REAL(real_8)                             :: pf
    INTEGER                                  :: ia, ib
    COMPLEX(real_8)                          :: psia(:), psib(:)
    INTEGER                                  :: iran
    COMPLEX(real_8)                          :: vpotg(:)
    REAL(real_8)                             :: vpotr(:)
    COMPLEX(real_8)                          :: psic(:), c2a(:)

    LOGICAL, PARAMETER                       :: a_stored_in_real = .FALSE.

    COMPLEX(real_8)                          :: fm, fp
    INTEGER                                  :: ig, ir
    LOGICAL                                  :: b_stored_in_real

    b_stored_in_real = (iran == 1)

    CALL hfxresp_get_pair_density(psia,psib,psic,a_stored_in_real,b_stored_in_real)

    CALL fwfftn(psic,.FALSE.,parai%allgrp)
    CALL hfxresp_get_coulomb(psic,vpotg,pf)
    CALL hfxresp_set_vpotg(psic,vpotg)
    CALL invfftn(psic,.FALSE.,parai%allgrp)
    CALL hfxresp_set_vpotr(psic,vpotr)

    IF (hfx_resp_env%store_v_ab) THEN
       IF (hfx_resp_env%write_v_ab) THEN
          DO ir=1,llr1
             WRITE(2000+parai%mepos) vpotr(ir)
          ENDDO
       ELSE
          DO ir=1,llr1
             hfx_resp%v_ab(ir,ia,ib)=vpotr(ir)
          ENDDO
       ENDIF
       hfx_resp%is_set=.TRUE.
    ENDIF

    CALL hfxresp_get_potential(psia,psic,vpotr,(.NOT. a_stored_in_real))
    IF (cntl%use_scaled_hfx) THEN
       CALL hfxresp_get_c2_real_space(psic,c2a)
    ELSE
       CALL fwfftn(psic,.TRUE.,parai%allgrp)
       CALL hfxrpav_get_c2(psic,c2a) 
    ENDIF

    ! ==--------------------------------------------------------------==
  END SUBROUTINE hfxrpav
  ! ==================================================================
  SUBROUTINE hfxrpav_incore(pf,ia,ib,psia,psib,iran,vpotg,vpotr,psic,c2a)
    ! ==================================================================
    REAL(real_8)                             :: pf
    INTEGER                                  :: ia, ib
    COMPLEX(real_8)                          :: psia(:), psib(:)
    INTEGER                                  :: iran
    COMPLEX(real_8)                          :: vpotg(:)
    REAL(real_8)                             :: vpotr(:)
    COMPLEX(real_8)                          :: psic(:), c2a(:)

    COMPLEX(real_8)                          :: fm, fp
    INTEGER                                  :: ig, ir
    REAL(real_8)                             :: rhelp

    DO ir=1,llr1
       ! goto position (ia-1)*NSTATE*LLR1 + (ib-1)* LLR1 + 1
       IF (hfx_resp_env%write_v_ab) THEN
          READ(2000+parai%mepos) rhelp
          psic(ir)=rhelp*REAL(psia(ir))
       ELSE
          psic(ir)=hfx_resp%v_ab(ir,ia,ib)*REAL(psia(ir))
       ENDIF
    ENDDO
    CALL fwfftn(psic,.TRUE.,parai%allgrp)
    DO ig=1,jgw
       fp=psic(nzfs(ig))+psic(inzs(ig))
       fm=psic(nzfs(ig))-psic(inzs(ig))
       c2a(ig)=c2a(ig)-0.5_real_8*CMPLX(REAL(fp),AIMAG(fm),kind=real_8)
    ENDDO
    ! ==--------------------------------------------------------------==

  END SUBROUTINE hfxrpav_incore
  ! ==================================================================

  ! ==================================================================
  ! Routines for hfxrpa and hfxpsi
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE hfxresp_get_pair_density(psia,psib,psic,a_stored_in_real,b_stored_in_real)
    ! ==--------------------------------------------------------------==
 
    COMPLEX(real_8), DIMENSION(:), &
                     INTENT(in)     :: psia, psib
    COMPLEX(real_8), DIMENSION(:), &
                     INTENT(inout)  :: psic
    LOGICAL, INTENT(in)             :: a_stored_in_real
    LOGICAL, INTENT(in)             :: b_stored_in_real
 
    INTEGER                         :: ir

    CALL zeroing(psic)
    IF (a_stored_in_real) THEN
       IF (b_stored_in_real) THEN
          !$omp parallel do private(ir)
          DO ir=1,llr1
             psic(ir)=REAL(psia(ir))*REAL(psib(ir))
          ENDDO
       ELSE
          !$omp parallel do private(ir)
          DO ir=1,llr1
             psic(ir)=REAL(psia(ir))*AIMAG(psib(ir))
          ENDDO
       ENDIF
    ELSE
       IF (b_stored_in_real) THEN
          !$omp parallel do private(ir)
          DO ir=1,llr1
             psic(ir)=AIMAG(psia(ir))*REAL(psib(ir))
          ENDDO
       ELSE
          !$omp parallel do private(ir)
          DO ir=1,llr1
             psic(ir)=AIMAG(psia(ir))*AIMAG(psib(ir))
          ENDDO
       ENDIF
    ENDIF
    CALL dscal(2*llr1,1._real_8/parm%omega,psic,1)
 
    ! ==--------------------------------------------------------------==
  END SUBROUTINE hfxresp_get_pair_density
  ! ==================================================================
  SUBROUTINE hfxresp_get_coulomb(psic,vpotg,pf)
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8), DIMENSION(:), &
                     INTENT(inout)  :: psic, vpotg

    REAL(real_8), INTENT(in)        :: pf

    INTEGER                         :: ig

    !$omp parallel do private(IG)
    DO ig=1,jhg
       vpotg(ig)=-pf*scgx(ig)*psic(nzff(ig))
    ENDDO

    ! ==--------------------------------------------------------------==
  END SUBROUTINE hfxresp_get_coulomb
  ! ==================================================================
  SUBROUTINE hfxresp_set_vpotg(psic,vpotg)
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8), DIMENSION(:), &
                     INTENT(inout)  :: psic, vpotg

    INTEGER                         :: ig

    CALL zeroing(psic)

    !ocl novrec
    !$omp parallel do private(IG)
    DO ig=1,jhg
       psic(nzff(ig))=vpotg(ig)
       psic(inzf(ig))=CONJG(vpotg(ig))
    ENDDO
    IF (geq0) psic(nzff(1))=vpotg(1)

    ! ==--------------------------------------------------------------==
  END SUBROUTINE hfxresp_set_vpotg
  ! ==================================================================
  SUBROUTINE hfxresp_set_vpotr(psic,vpotr)
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8), DIMENSION(:), &
                     INTENT(in)     :: psic
    REAL(real_8), DIMENSION(:), &
                  INTENT(inout)     :: vpotr

    INTEGER                         :: ir

    !$omp parallel do private(IR)
    DO ir=1,llr1
       vpotr(ir) = REAL(psic(ir))
    ENDDO

    ! ==--------------------------------------------------------------==
  END SUBROUTINE hfxresp_set_vpotr
  ! ==================================================================
  SUBROUTINE hfxresp_get_potential(psi,psic,vpotr,stored_in_real)
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8), DIMENSION(:), &
                     INTENT(in)     :: psi
    COMPLEX(real_8), DIMENSION(:), &
                     INTENT(inout)  :: psic
    REAL(real_8), DIMENSION(:), &
                  INTENT(in)        :: vpotr
    LOGICAL, INTENT(in)             :: stored_in_real

    INTEGER                         :: ir

    CALL zeroing(psic)

    IF (stored_in_real) THEN
       !$omp parallel do private(ir)
       DO ir=1,llr1
          psic(ir) = vpotr(ir)*REAL(psi(ir))
       ENDDO
    ELSE
       !$omp parallel do private(ir)
       DO ir=1,llr1
          psic(ir) = vpotr(ir)*AIMAG(psi(ir))
       ENDDO
    ENDIF

    ! ==--------------------------------------------------------------==
  END SUBROUTINE hfxresp_get_potential
  ! ==================================================================
  SUBROUTINE hfxresp_get_c2_real_space(psic,c2)
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8), DIMENSION(:), &
                     INTENT(in)       :: psic
    COMPLEX(real_8), DIMENSION(:), &
                     INTENT(inout)    :: c2

    INTEGER         :: ir

    !$omp parallel do private(ir)
    DO ir=1,llr1
       c2(ir) = c2(ir) + REAL(psic(ir))
    ENDDO

    ! ==--------------------------------------------------------------==
  END SUBROUTINE hfxresp_get_c2_real_space
  ! ================================================================== 
  SUBROUTINE hfxpb_get_c2(psic,c2)
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8), DIMENSION(:), &
                     INTENT(in)       :: psic
    COMPLEX(real_8), DIMENSION(:), &
                     INTENT(inout)    :: c2

    COMPLEX(real_8) :: fp, fm
    INTEGER         :: ig

    !$omp parallel do private(ig,fp,fm)
    DO ig=1,jgw
       fp     = psic(nzfs(ig))+psic(inzs(ig))
       fm     = psic(nzfs(ig))-psic(inzs(ig))
       c2(ig) = c2(ig)+0.5_real_8*CMPLX(REAL(fp),AIMAG(fm),kind=real_8)
    ENDDO

    ! ==--------------------------------------------------------------==
  END SUBROUTINE hfxpb_get_c2
  ! ==================================================================
  SUBROUTINE hfxrpav_get_c2(psic,c2)
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8), DIMENSION(:), &
                     INTENT(in)       :: psic
    COMPLEX(real_8), DIMENSION(:), &
                     INTENT(inout)    :: c2

    COMPLEX(real_8) :: fp, fm
    INTEGER         :: ig

    !$omp parallel do private(ig,fp,fm)
    DO ig=1,jgw
       fp     = psic(nzfs(ig))+psic(inzs(ig))
       fm     = psic(nzfs(ig))-psic(inzs(ig))
       c2(ig) = c2(ig)-0.5_real_8*CMPLX(REAL(fp),AIMAG(fm),kind=real_8)
    ENDDO

    ! ==--------------------------------------------------------------==
  END SUBROUTINE hfxrpav_get_c2
  ! ==================================================================
END MODULE pw_hfx_resp
! ==================================================================
