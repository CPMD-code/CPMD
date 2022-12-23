MODULE ksmat_dist_utils
  USE atwf,                            ONLY: atwf_mod,&
                                             atwp,&
                                             catom,&
                                             zxmat
  USE error_handling,                  ONLY: stopgm
  USE fnlalloc_utils,                  ONLY: fnl_set,&
                                             fnlalloc,&
                                             fnldealloc
  USE fnonloc_utils,                   ONLY: fnonloc,&
                                             give_scr_fnonloc
  USE gsortho_utils,                   ONLY: gs_ortho
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE jrotation_utils,                 ONLY: set_orbdist
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE ksmatmod,                        ONLY: eig_vec,&
                                             z
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_group,&
                                             mp_sum,&
                                             mp_sync
  USE parac,                           ONLY: parai
  USE prng_utils,                      ONLY: repprngu_vec
  USE rnlsm_utils,                     ONLY: give_scr_rnlsm,&
                                             rnlsm
  USE spin,                            ONLY: lspin2,&
                                             spin_mod
  USE system,                          ONLY: cnti,&
                                             nkpt,&
                                             paraw
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE vpsi_utils,                      ONLY: vpsi
  USE wfnio_utils,                     ONLY: queryrands
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: dist_ksmat
  PUBLIC :: give_scr_dist_ksmat

CONTAINS

  ! ==================================================================
  SUBROUTINE dist_ksmat(c2,vpot,psi,nstate,ikind,c0,nn,&
       tlsd2, ainitwfflag)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==  THE KOHN-SHAM MATRIX IN THE ATOMIC WAVEFUNCTION BASIS       ==
    ! ==  IT'S EIGENDECOMPOSITION AND THE APPROXIMATE INITIAL         ==
    ! ==  WAVEFUNCTIONS                                               ==
    ! ==  USES: DISTRIBUTED LANCZOS WITH GS REORTHOGONALIZATION       ==
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8)                          :: c2(:,:)
    REAL(real_8)                             :: vpot(:,:)
    COMPLEX(real_8)                          :: psi(:)
    INTEGER                                  :: nstate, ikind
    COMPLEX(real_8)                          :: c0(:,:,:)
    INTEGER                                  :: nn
    LOGICAL                                  :: tlsd2, ainitwfflag

    CHARACTER(*), PARAMETER                  :: procedureN = 'dist_ksmat'

    COMPLEX(real_8)                          :: pab(1)
    INTEGER :: chunk_begin, chunk_begin_e, chunk_end, chunk_end_e, chunk_new, &
      chunk_new_e, from_beg, from_end, i, ia, ibeg, ierr, ifail(atwp%nattot), &
      index1, index2, index4, is, ist, isub, iwork(5*atwp%nattot), lan_max, &
      langrp, lanlist(parai%nproc), lc_index(atwp%nattot), max_n, msglen, &
      n_max, n_seed, natst, nfound, nolan, norb, norbx, num_eig, &
      old_send_cnt_e(parai%nproc), rank, send_cnt(parai%nproc), &
      send_cnt_e(parai%nproc), send_displ(parai%nproc), &
      send_displ_e(parai%nproc), start, start_col, start_spd, start_spu, sz, &
      sz_old, tm1, tm2
    INTEGER, ALLOCATABLE                     :: seed(:)
    LOGICAL                                  :: conv_flag, first_check, &
                                                flag_spdown, flag_spu
    REAL(real_8) :: abtol, alpha, beta, conv_tol, d(atwp%nattot), dnrm2, &
      e(atwp%nattot), eig(atwp%nattot), ld(atwp%nattot), le(atwp%nattot), &
      local_alpha, local_beta, nrm, sm_old, vl, vu, work(5*atwp%nattot)
    REAL(real_8), ALLOCATABLE :: e_temp(:), lanv(:), local_pw(:), &
      local_xxmat(:,:), pw(:), rnd(:), temp1_xxmat(:,:), temp_xxmat(:,:), &
      v0(:), v_prev(:), vcur(:), w(:), z_temp(:)
    REAL(real_8), DIMENSION(20)              :: foc = 1.0_real_8
    REAL(real_8), EXTERNAL                   :: ddot

    CALL tiset(procedureN,isub)
    tm1 = 0
    tm2 = 0
    d(:) = 0.0_real_8
    ! FIND OUT SIZE AND RANK
    rank = parai%me
    sz = parai%nproc
    conv_flag = .FALSE.
    first_check = .TRUE.
    sz_old = sz
    ! =================================================================
    ! DETERMINE DISTRIBUTION OF ATOMIC STATES AND COLUMNS OF
    ! RESTRICTED HAMILTONIAN
    ! =================================================================
    CALL set_orbdist(atwp%nattot,cnti%nstblk,parai%nproc,norbx)
    norb = 0
    DO index1=0, sz-1
       IF ((paraw%nwa12(index1,2)-paraw%nwa12(index1,1)+1 ).GT.0) norb=norb+1
    ENDDO
    chunk_new = INT(REAL(atwp%nattot,kind=real_8)/REAL(norb,kind=real_8))
    IF (rank.LT.norb-1) THEN
       chunk_begin = rank*chunk_new+1
       chunk_end = (rank+1)*chunk_new
    ELSE IF (rank.EQ.norb-1) THEN
       chunk_end=atwp%nattot
       chunk_begin = rank*chunk_new+1
    ENDIF
    IF (rank.GE.norb) THEN
       chunk_begin=atwp%nattot+1
       chunk_end=atwp%nattot
    ENDIF
    chunk_new = chunk_end - chunk_begin + 1
    msglen=8/2
    CALL my_concat(chunk_begin, paraw%nwa12(0,1), msglen, parai%allgrp)
    CALL my_concat(chunk_end, paraw%nwa12(0,2), msglen, parai%allgrp)
    ! ================================================================= 
    ! CREATE LANCZOS GROUP
    ! ================================================================= 
    IF (parai%me.EQ.parai%source) THEN ! vw here should be PARENT
       DO index4=1,sz
          lanlist(index4) = 0
       ENDDO
       nolan = 0
       DO index4=1, sz
          IF ((paraw%nwa12(index4-1,2)-paraw%nwa12(index4-1,1)+1).NE.0) THEN
             lanlist(index4) = index4-1
             nolan = nolan + 1
          ENDIF
       ENDDO
    ENDIF
    CALL mp_bcast(lanlist, sz, parai%source, parai%allgrp)
    CALL mp_bcast(nolan, parai%source, parai%allgrp)
    CALL mp_group(nolan,lanlist,langrp,parai%allgrp)
    ! ================================================================= 
    ! ALLOCATE LOCAL_XXMAT AND TEMP_XXMAT, TEMP1_XXMAT
    ALLOCATE(local_xxmat(atwp%nattot,chunk_new),&
         temp_xxmat(atwp%nattot,atwp%numaormax),&
         temp1_xxmat(atwp%nattot,atwp%numaormax),&
         stat=ierr)
    IF (ierr.NE.0) CALL stopgm( procedureN, 'Allocation problem 1' ,& 
         __LINE__,__FILE__)
    CALL zeroing(local_xxmat)!, atwp%nattot*chunk_new)
    CALL zeroing(temp1_xxmat)
    CALL zeroing(temp_xxmat)
    ! ==--------------------------------------------------------------==
    IF (lspin2%tlse) CALL stopgm('KSMAT','NO LSE ALLOWED HERE',& 
         __LINE__,__FILE__)
    ! AK      CALL Zeroing(PAB)!,NNR1)
    pab(1)=0.0_real_8
    ! ==--------------------------------------------------------------==
    CALL fnl_set('SAVE')
    CALL fnlalloc(atwp%numaormax,.FALSE.,.FALSE.)
    ist=1
    index4 = 1
    DO is=1,ions1%nsp
       natst=atwf_mod%numaor(is)
       DO ia=1,ions0%na(is)
          CALL rnlsm(catom(:,ist:ist+natst-1),natst,1,ikind,.FALSE.)
          CALL zeroing(c2(:,1:natst))!,nkpt%ngwk*natst)
          CALL vpsi(catom(:,ist:ist+natst-1),c2,foc,vpot,psi,natst,ikind,1,.TRUE.)
          CALL fnonloc(c2,foc,natst,ikind,1,.TRUE.)
          IF (tkpts%tkpnt) THEN
             CALL ovlap2_c(nkpt%ngwk,atwp%nattot,natst,zxmat(1,ist),catom,c2)
          ELSE
             CALL ovlap2(nkpt%ngwk,atwp%nattot,natst,temp1_xxmat(1,1),catom,c2,&
                  .TRUE.)
             ! SUMMATION ACROSS PROCESSORS
             CALL mp_sum(temp1_xxmat,temp_xxmat,atwp%nattot*natst,parai%allgrp)
             ! INTIALIZE TO ZERO
             CALL zeroing(temp1_xxmat(:,1:atwp%numaormax))!, atwp%nattot*atwp%numaormax)

             ! ==--------------------------------------------------------------==
             from_beg = MAX(chunk_begin,ist)
             from_end = MIN(chunk_end,ist+natst-1)
             n_max  = MAX(0,from_end-from_beg+1)
             ibeg = from_beg - ist + 1
             IF (n_max.GT.0) THEN
                IF (index4.GT.SIZE(local_xxmat,2))&
                     CALL stopgm( procedureN, 'PROBLEM HERE 1' ,& 
                     __LINE__,__FILE__)
                IF (index4+n_max-1.GT.SIZE(local_xxmat,2))&
                     CALL stopgm( procedureN, 'PROBLEM HERE 2' ,& 
                     __LINE__,__FILE__)
                IF (ibeg.GT.SIZE(temp_xxmat,2))&
                     CALL stopgm( procedureN, 'PROBLEM HERE 3' ,& 
                     __LINE__,__FILE__)
                IF (ibeg+n_max-1.GT.SIZE(temp_xxmat,2))&
                     CALL stopgm( procedureN, 'PROBLEM HERE 4' ,& 
                     __LINE__,__FILE__)
                CALL dcopy(atwp%nattot*n_max, temp_xxmat(1,ibeg), 1,&
                     local_xxmat(1,index4), 1)
             ENDIF
             index4 = index4 + n_max
          ENDIF
          CALL zeroing(temp_xxmat(:,atwp%numaormax))!, atwp%nattot*atwp%numaormax)
          CALL mp_sync(parai%allgrp)
          ist=ist+natst
       ENDDO
    ENDDO
    CALL fnldealloc(.FALSE.,.FALSE.)
    CALL fnl_set('RECV')
    IF (tkpts%tkpnt) THEN
       CALL dscal(2*atwp%nattot*atwp%nattot,-1._real_8,zxmat(1,1),1)
    ELSE
       CALL dscal(atwp%nattot*chunk_new,-1._real_8,local_xxmat,1)
    ENDIF
    ! ==================================================================
    ! ONLY PROCESSORS THAT HAVE A NON-ZERO CHUNK OF THE RESTRICTED
    ! HAMILTONIAN PARTICIPATE IN LANCZOS
    ! ==================================================================   
    IF (chunk_new.NE.0) THEN
       ! ==================================================================
       ! START THE LANCZOS ITERATION WITH FULL GRAM-SCHIDT REORTH. 
       ! ==================================================================
       ! --------------------------------------
       ! MEMORY ALLOCATIONS AND INITIALIZATIONS
       ! --------------------------------------
       ! EACH PROCESSOR ALLOCATES MEMORY FOR ROW-WISE DISTRUBUTED LANCZOS
       ! BASIS
       sz = nolan
       lan_max = atwp%nattot+1
       ALLOCATE(vcur(atwp%nattot),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(v0(atwp%nattot),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(lanv(chunk_new*lan_max),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(w(chunk_new),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(pw(lan_max),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(local_pw(lan_max),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(v_prev(chunk_new),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(pw)!,lan_max)
       CALL zeroing(lanv)!,chunk_new*lan_max)
       CALL zeroing(vcur)!,atwp%nattot)
       CALL zeroing(v0)!,atwp%nattot)
       CALL zeroing(w)!,chunk_new)
       CALL zeroing(local_pw)!,lan_max)
       CALL zeroing(v_prev)!,chunk_new)
       ! DETERMINE SEND AND RECEIVE BUFFERS
       msglen=8/2
       CALL my_concat(chunk_new, send_cnt(1), msglen, langrp)
       send_displ(1) = 0
       DO index4=2, sz
          send_displ(index4) = send_displ(index4-1) + send_cnt(index4-1)
       ENDDO
       ! DETERMINE A UNIT NORM INITIAL VECTOR V0 
       IF (parai%me.EQ.parai%source) THEN ! vw here should be PARENT
          CALL repprngu_vec(atwp%nattot,v0) 
          nrm = dnrm2(atwp%nattot, v0(1), 1)
          nrm = 1/nrm
          CALL dscal(atwp%nattot, nrm, v0(1), 1)
       ENDIF
       ! BROADCAST INITIAL VECTOR TO ALL PROCESSORS
       CALL mp_bcast(v0, atwp%nattot, 0, langrp)
       CALL dcopy(atwp%nattot, v0(1), 1, vcur(1), 1)
       ! INITIALIZE V_PREV and BETA
       DO is=1, chunk_new
          v_prev(is) = 0._real_8
          w(is)=0._real_8
       ENDDO
       beta = 0._real_8
       ! INITIALIZE FIRST LANCZOS VECTOR
       CALL dcopy(chunk_new, vcur(chunk_begin), 1, lanv(1), 1)
       ! INITIALIZE INITIAL SUM OF WANTED EIGENVALUES
       SM_OLD = 0.0_real_8
       ! INITIALIZE CONVERGENCE TOLERANCE
       conv_tol = 1.e-8_real_8
       ! ==================================================================
       ! LANCZOS LOOP
       ! ==================================================================   
       num_eig = nn
       DO index4 = 1, lan_max-1
          ! DISTRIBUTED MATRIX-VECTOR PRODUCT
          ! W: Holds the result
          CALL dgemv("T", atwp%nattot,chunk_new,1._real_8,local_xxmat(1,1),atwp%nattot,&
               vcur(1), 1, 0._real_8, w(1), 1)
          ! DAXPY: W = W - BETA*V_(J-1)
          CALL daxpy(chunk_new, -beta, v_prev(1), 1, w(1), 1)
          ! CALCULATE ALPHA: INNER PRODUCT WITH CURRENT LANCZOS VECTOR
          local_alpha = ddot(chunk_new,w(1),1,&
               lanv((index4-1)*chunk_new+1), 1)
          CALL mp_sum(local_alpha, alpha, langrp)
          ! DAXPY: W = W - ALPHA*V_J
          CALL daxpy(chunk_new,-alpha,&
               lanv((index4-1)*chunk_new+1),1,w(1),1)
          ! ================================================================== 
          ! REORTHOGONALIZE AGAINST ALL PREVIOUS BASIS VECTORS
          ! ==================================================================
          ! CALL DGEMV("T",CHUNK_NEW,INDEX4,1._real_8,LANV(1),CHUNK_NEW, W(1), 1,
          ! &              0._real_8, LOCAL_PW(1),1)
          DO index1=1,index4
             local_pw(index1)=ddot(chunk_new,&
                  lanv((index1-1)*chunk_new+1), 1, w(1), 1)
          ENDDO
          CALL mp_sum(local_pw, pw, index4, langrp)
          DO index1=1, index4
             CALL daxpy(chunk_new, -pw(index1), lanv((index1-1)*&
                  chunk_new+1), 1, w(1), 1)
          ENDDO
          ! ==================================================================  
          ! CALCULATE BETA: NORM OF VECTOR W
          local_beta = ddot(chunk_new, w(1), 1, w(1), 1)
          CALL mp_sum(local_beta, beta, langrp)
          beta = SQRT(beta)
          ! vw beta can be zero in some special cases, so we protect it.
          IF ( beta == 0.0_real_8 ) THEN
             conv_flag=.TRUE.
             GOTO 666
          ENDIF
          ! NORMALIZE NEXT LANCZOS VECTOR
          CALL dscal(chunk_new, 1.0_real_8/beta, w(1), 1)
          CALL dcopy(chunk_new, w(1), 1, lanv(index4*chunk_new+1), 1)
          ! UPDATE TRIDIAGONAL MATRIX T (DIAGONAL AND SUP_DIAGONAL ENTRIES) 
          d(index4) = alpha
          e(index4) = beta
          ! UPDATE AUXILIARY VECTORS
          CALL my_concatv(lanv((index4)*chunk_new+1),vcur(1),chunk_new,&
               send_cnt(1),send_displ(1),langrp)
          CALL dcopy(chunk_new,lanv((index4-1)*chunk_new+1),&
               1,v_prev(1),1)
          ! CHECK CONVERGENCE FOR EIGENVALUES
          CALL dcopy(index4, d(1), 1, ld(1), 1)
          CALL dcopy(index4-1, e(1), 1, le(1), 1)
          abtol = -1
          ! IF (PARENT) WRITE(6,*) ALPHA, BETA
          ! ================================================================
          ! THIS TEST IS PRACTICALLY NEVER DONE: FOR FUTURE USE
          !>
          !          IF (index4.GT.num_eig*2) THEN
          !             ! to avoid 'not allocated' runtime error
          !             ALLOCATE(z(1),STAT=ierr)
          !             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
          !                  __LINE__,__FILE__)
          !             CALL  dstevx('N', 'I', index4, ld, le, vl, vu, 1,&
          !                  num_eig, abtol,nfound, eig, z, 1, work,&
          !                  iwork, ifail, ierr)
          !             DEALLOCATE(z,STAT=ierr)
          !             IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
          !                  __LINE__,__FILE__)
          !             sm = dasum(index4, eig(1), 1)
          !             CALL mp_bcast(sm, 0, langrp )
          !             IF( first_check ) THEN
          !                first_check = .FALSE.
          !             ELSE
          !                IF ( (ABS(sm-sm_old)/sm_old).LE.conv_tol ) THEN
          !                   conv_flag=.TRUE.
          !                   GOTO 666
          !                ENDIF
          !             ENDIF
          !             sm_old = sm
          !          ENDIF
          !<
          ! ================================================================
       ENDDO
666    CONTINUE
       IF( conv_flag .AND. index4 /= atwp%nattot ) CALL stopgm(procedureN,'conv_flag shoudnt be TRUE', &
            __LINE__,__FILE__)
       ! ===================
       ! END OF LANCZOS LOOP
       ! ===================
       IF (.NOT.conv_flag) index4 = index4-1
       ! ==================================================================
       ! DETERMINE THE EIGENVECTORS THAT EACH PROCESSOR IS GOING TO COMPUTE
       ! ==================================================================
       norb=paraw%nwa12(parai%mepos,2)-paraw%nwa12(parai%mepos,1)+1
       norb=MAX(0,norb)
       chunk_begin_e = paraw%nwa12(parai%mepos,1)
       chunk_end_e   = paraw%nwa12(parai%mepos,2)
       chunk_new_e   = chunk_end_e - chunk_begin_e +1
       ! ==================================================================
       ! DETERMINE SEND  BUFFERS
       ! ==================================================================
       msglen=8/2
       CALL my_concat(chunk_new_e, send_cnt(1), msglen, langrp)
       send_displ(1) = 0
       DO index1=2, sz
          send_displ(index1) = send_displ(index1-1) + send_cnt(index1-1)
       ENDDO
       ! SDF  Add random noise to the diangonal element of the matrix to break
       ! SDF  degenerate eigenvalues.
       ALLOCATE(rnd(index4),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       CALL RANDOM_SEED(size=n_seed)
       ALLOCATE(seed(n_seed),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       seed = (/ (i-1, i=1,n_seed) /)
       CALL RANDOM_SEED(put=seed)
       DEALLOCATE(seed,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       CALL RANDOM_NUMBER(rnd)
       d=d*(1.0_real_8+1.0e-12_real_8*(2.0_real_8*rnd-1.0_real_8))
       DEALLOCATE(rnd,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       ! ALLOCATE MEMORY FOR EIGENVECTORS
       ALLOCATE(z(index4*(chunk_end_e-chunk_begin_e+1)),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL  dstevx('V','I',index4,d,e,vl,vu,chunk_begin_e,&
            chunk_end_e,&
            abtol, nfound, eig, z, index4, work,&
            iwork, ifail, ierr)
       ! ==================================================================
       ! CALCULATION OF APPROXIMATE EIGENVECTORS
       ! EACH PROCESSOR WILL STORE COMPLETE EIGENVECTORS 
       ! EIGENVECTORS ARE DISTRIBUTED ACROSS PROCESSORS
       ! ==================================================================
       ! ALLOCATION OF MEMORY FOR APPROXIMATE EIGENVECTORS
       ALLOCATE(eig_vec(atwp%nattot*chunk_new_e),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ! ALLOCATE ADDITIONAL MEMORY FOR BOOK KEEPING
       msglen=8/2
       CALL my_concat(chunk_new, old_send_cnt_e(1), msglen, langrp)
       ALLOCATE(z_temp(MAXVAL(send_cnt(1:sz)*atwp%nattot)),&
            e_temp(chunk_new*MAXVAL(send_cnt(1:sz))),&
            stat=ierr)
       IF (ierr.NE.0) CALL stopgm( procedureN, 'Allocation problem 2' ,& 
            __LINE__,__FILE__)
       DO index1=1, nolan
          IF (rank.EQ.(index1-1) ) THEN
             CALL dcopy(index4*chunk_new_e, z(1), 1, z_temp(1), 1)
             DO index2 = 1, sz
                send_cnt_e(index2) = old_send_cnt_e(index2) *&
                     send_cnt(rank+1)
                ! IF (SEND_CNT(INDEX2).EQ.0) SEND_CNT_E(INDEX2) = 0
             ENDDO
             send_displ_e(1) = 0
             DO index2=2, sz
                send_displ_e(index2) = send_displ_e(index2-1) +&
                     send_cnt_e(index2-1)
             ENDDO
          ELSE
             send_cnt_e=-HUGE(0)
             send_displ_e=-HUGE(0)
          ENDIF
          ! BROADCAST COUNTS AND DISPL. DATA FOR ALLGATHER
          CALL mp_bcast(send_cnt_e, sz, index1-1, langrp)
          CALL mp_bcast(send_displ_e, sz, index1-1, langrp)
          ! BROADCAST MY EIGENVECTORS OF THE TRIDIAGONAL MATRIX
          CALL mp_bcast(z_temp, send_cnt(index1)*index4, index1-1,&
               langrp)
          ! MULTIPLY RECEIVED EIGENVECTORS WITH MY PART OF THE LANCZOS
          ! BASIS
          CALL dgemm("N","N",chunk_new,send_cnt(index1),index4,1._real_8,&
               lanv(1),chunk_new,z_temp(1),index4,&
               0._real_8, e_temp, chunk_new)

          ! BRING RESULTS BACK TO THE CORESSPONDING PROCESSOR 
          CALL my_source_concatv(e_temp,eig_vec(1),&
               chunk_new*send_cnt(index1),&
               send_cnt_e, send_displ_e,&
               index1-1, langrp)

       ENDDO
       ! ==================================================================
       DO index1 = 1, atwp%nattot
          lc_index(index1) = 0
       ENDDO
       start     = 1
       start_col = 1
       DO index1 = 1, sz
          DO index2 = 1, chunk_new_e
             CALL dcopy(old_send_cnt_e(index1), eig_vec(start_col), 1,&
                  z_temp(start+lc_index(1)), 1)
             start_col = start_col + old_send_cnt_e(index1)
             start     = start + atwp%nattot
          ENDDO
          lc_index(1) = lc_index(1) + old_send_cnt_e(index1)
          start = 1
       ENDDO
       CALL dcopy(atwp%nattot*chunk_new_e, z_temp(1), 1, eig_vec(1), 1)
       ! ==================================================================
       ! END LANCZOS ITERATION ON LANGRP
       ! ==================================================================
    ELSE
       ! PROCESSORS THAT DID NOT PARTICIPATE IN LANCZOS
       ! WILL WAIT IN SUBSEQUENT COLLECTIVE 
    ENDIF
    ! SEND CHUNK INFORMATION TO ALL PROCESSORS
    msglen=8/2
    CALL my_concat(chunk_new_e, send_cnt(1), msglen, parai%allgrp)
    ! FIND LARGEST NUMBER OF APPROX. EIGENVECTORS HELP BY PROCEESORS
    IF (chunk_new.EQ.0) THEN
       chunk_new_e = send_cnt(1)
       DO index1=2,nolan
          IF (send_cnt(index1).GT.chunk_new_e)&
               chunk_new_e = send_cnt(index1)
       ENDDO
       ALLOCATE(z_temp(chunk_new_e*atwp%nattot), stat=ierr)
       IF (ierr.NE.0) CALL stopgm( procedureN, 'Allocation problem 3' &
            ,& 
            __LINE__,__FILE__)
    ENDIF
    ! ==================================================================
    ! COMPUTE APPROXIMATE "ATOMIC" WAVEFUNCTIONS
    ! ==================================================================
    start_col = 0
    start = 1
    ! ===========================
    ! NEEDED FOR LSD CALCULATIONS
    start_spu = 1
    start_spd = 1
    flag_spu = .TRUE.
    flag_spdown = .TRUE.
    ! ===========================
    DO index1=1, nolan
       IF (rank.EQ.(index1-1) ) THEN
          ! INITIALIZE TEMP. SPACE TO ZERO
          CALL zeroing(z_temp)!, atwp%nattot*chunk_new_e)
          CALL dcopy(atwp%nattot*chunk_new_e, eig_vec(1), 1, z_temp(1), 1)
       ENDIF
       ! CURRENT ROOT BROADCASTS ITS APPROX. EIGENVECTORS 
       CALL mp_bcast(z_temp, send_cnt(index1)*atwp%nattot, index1-1,&
            parai%allgrp)
       ! EACH PROCESSOR MULTIPLIES THE APPROXIMATE EIGVECS
       ! RECEIVED AND STORES THEM AT THE APPROPRIATE POSITION
       ! IN MATRIX C0
       ! ==============================================================
       ! HANDLE CASE FOR LOCAL SPIN DENSITY CALCULATION (LSD)
       ! TLSD2 IS READ AS AN ARGUMENT
       ! RESULT FOR SPIN-UP FIRST THEN FOR SPIN-DOWN
       ! ==============================================================
       IF (tlsd2) THEN
          ! SPIN UP
          IF (start_spu.LE.spin_mod%nsup) THEN
             max_n = send_cnt(index1)
             IF (start_spu+send_cnt(index1).GT.spin_mod%nsup)&
                  max_n = spin_mod%nsup-start_spu+1
             CALL dgemm("N","N",2*nkpt%ngwk,max_n,&
                  atwp%nattot,1._real_8,&
                  catom(1,1),2*nkpt%ngwk,z_temp(1),atwp%nattot,&
                  0._real_8, c0(1,start_spu,ikind), 2*nkpt%ngwk)
             start_spu = start_spu + send_cnt(index1)
          ENDIF
          ! SPIN DOWN 
          IF (start_spd.LE.spin_mod%nsdown) THEN
             max_n = send_cnt(index1)
             IF (start_spd+send_cnt(index1).GT.spin_mod%nsdown)&
                  max_n = spin_mod%nsdown-start_spd+1
             CALL dgemm("N","N",2*nkpt%ngwk,max_n,&
                  atwp%nattot,1._real_8,&
                  catom(1,1),2*nkpt%ngwk,z_temp(1),atwp%nattot,&
                  0._real_8, c0(1,spin_mod%nsup+start_spd,ikind), 2*nkpt%ngwk)
             start_spd = start_spd + send_cnt(index1)
          ENDIF
       ELSE
          ! ===================
          ! HANDLE NON LSD CASE
          ! ===================
          IF (start.LE.nstate) THEN
             max_n = send_cnt(index1)
             IF (start+send_cnt(index1).GT.nstate)&
                  max_n = nstate-start+1
             IF (.NOT.ainitwfflag) THEN
                CALL dgemm("N","N",2*nkpt%ngwk,max_n,atwp%nattot,&
                     1._real_8,&
                     catom(1,1),2*nkpt%ngwk,z_temp(1),atwp%nattot,&
                     0._real_8, c0(1,start,ikind), 2*nkpt%ngwk)
                start = start + send_cnt(index1)
             ELSE
                DO index2=1, max_n
                   IF (queryrands(start + index2 -1)) THEN
                      CALL dgemv("N",2*nkpt%ngwk,atwp%nattot,1._real_8,&
                           catom(1,1),2*nkpt%ngwk,&
                           z_temp((index2-1)*atwp%nattot+1),1,&
                           0._real_8, c0(1,start+index2-1,ikind), 1)
                      CALL gs_ortho(c0(:,:,ikind),start+index2-2,&
                           c0(:,start+index2-1:,ikind),1)

                   ENDIF
                ENDDO
                start = start + send_cnt(index1)
             ENDIF
          ENDIF
       ENDIF
    ENDDO
    ! ==================================================================
    ! DEALLOCATIONS OF MEMORY
    ! ==================================================================
    DEALLOCATE(z_temp,stat=ierr)
    IF (ierr.NE.0) CALL stopgm(procedureN,'Deallocation problem 3',& 
         __LINE__,__FILE__)
    IF (chunk_new.NE.0) THEN
       DEALLOCATE(lanv,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(w,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(pw,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(local_pw,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(v_prev,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(z,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(e_temp,stat=ierr)
       IF (ierr.NE.0) CALL stopgm(procedureN,'Deallocation problem 2',& 
            __LINE__,__FILE__)
       DEALLOCATE(eig_vec,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(vcur,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(v0,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! GLOBAL DEALLOCATIONS
    DEALLOCATE(local_xxmat,temp_xxmat,temp1_xxmat,stat=ierr)
    IF (ierr.NE.0) CALL stopgm( procedureN, 'Deallocation problem 1',& 
         __LINE__,__FILE__)
    ! =================================================================
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE dist_ksmat
  ! ==================================================================
  SUBROUTINE give_scr_dist_ksmat(lksmat,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lksmat
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: il_auxc, il_ddia, lfnonloc, &
                                                lrnlsm

    CALL give_scr_rnlsm(lrnlsm,tag,atwp%nattot,.FALSE.)
    CALL give_scr_fnonloc(il_auxc,il_ddia,atwp%numaormax)
    lfnonloc = il_auxc + il_ddia
    lksmat=MAX(lfnonloc,lrnlsm)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_dist_ksmat
  ! ==================================================================

END MODULE ksmat_dist_utils
