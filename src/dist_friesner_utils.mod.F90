MODULE dist_friesner_utils
  USE adjmu_utils,                     ONLY: convfrie
  USE dotp_utils,                      ONLY: dotp
  USE error_handling,                  ONLY: stopgm
  USE fint,                            ONLY: fint1
  USE friesner_utils,                  ONLY: give_scr_krylov_ref,&
                                             krylov_ref
  USE gsortho_utils,                   ONLY: gs_ortho
  USE hfx_drivers,                     ONLY: hfxpsi
  USE hpsi_utils,                      ONLY: give_scr_hpsi,&
                                             hpsi
  USE jrotation_utils,                 ONLY: my_set_orbdist
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE machine,                         ONLY: m_walltime
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_group,&
                                             mp_sum,&
                                             mp_sync
  USE parac,                           ONLY: parai,&
                                             paral
  USE prng_utils,                      ONLY: repprngu_vec
  USE randtowf_utils,                  ONLY: randtowf
  USE sort_utils,                      ONLY: sort2
  USE summat_utils,                    ONLY: give_scr_summat
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             cntr,&
                                             ncpw,&
                                             paraw
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: dist_friesner
  PUBLIC :: give_scr_dist_friesner

CONTAINS

  ! ==================================================================
  SUBROUTINE dist_friesner(ndiag,nel,nconv,nhpsi,nstate,&
       c0,cs,sc0,cscr,cgs,focc,vpot,psi,edav,jspin,trefine,tinfo)
    ! ==--------------------------------------------------------------==
    ! ==  Kohn-Sham Matrix Diagonalization by Krylov-Space Method     ==
    ! ==  W.T. Pollard and R.A. Friesner                              ==
    ! ==  J. Chem. Phys. 99, 6742 (1993)                              ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ndiag
    REAL(real_8)                             :: nel
    INTEGER                                  :: nconv, nhpsi, nstate
    COMPLEX(real_8)                          :: c0(ncpw%ngw,*), &
                                                cs(ncpw%ngw,*), &
                                                sc0(ncpw%ngw,*), &
                                                cscr(ncpw%ngw,*), &
                                                cgs(ncpw%ngw,*)
    REAL(real_8)                             :: focc(nstate), vpot(*)
    COMPLEX(real_8)                          :: psi(:)
    REAL(real_8)                             :: edav(*)
    INTEGER                                  :: jspin
    LOGICAL                                  :: trefine, tinfo

    CHARACTER(*), PARAMETER                  :: procedureN = 'dist_friesner'
    INTEGER, PARAMETER                       :: ispin = 1 
    REAL(real_8), PARAMETER                  :: big = 1.e33_real_8 , &
                                                edavmax = 300._real_8, &
                                                edavmin = -300._real_8
    CHARACTER(len=20)                        :: charspin
    INTEGER :: chunk_begin, chunk_begin_e, chunk_end, chunk_end_e, chunk_new, &
      chunk_new_e, i, iaux, icycle, ierr, ig, index1, index2, index4, ip, is, &
      isub, itest, j, lan_max, langrp, lanlist(parai%nproc), msglen, &
      nconvold, ncurr, iwork(5*ndiag), ifail(ndiag), nfound, ngw2, nhpsiold, &
      nleft, nolan, norb, norbx, lc_index(nstate), ntest, num_eig, nx, &
      old_send_cnt_e(parai%nproc), rank, send_cnt(parai%nproc), &
      send_cnt_e(parai%nproc), send_displ(parai%nproc), &
      send_displ_e(parai%nproc), start, start_col, sz
    INTEGER, ALLOCATABLE                     :: index0(:)
    LOGICAL                                  :: conv_flag, mem_flag, tcheck, &
                                                tconv, tdebug, tlsd2, ttest
    REAL(real_8) :: abtol, alpha, amu, aux, b2, b2max, b2min, beta, conv_tol, &
      d(ndiag), dnrm2, e(ndiag), fac, ld(ndiag), le(ndiag), local_alpha, &
      local_beta, nrm, sign, sm_old, teig(ndiag), tim, tim1, tim2, trotmax, &
      trotmin, vl, vu, work(5*ndiag)
    REAL(real_8), PARAMETER                  :: eps = 1.e-14_real_8

    REAL(real_8), ALLOCATABLE :: adav(:,:), alpha0(:), beta1(:), e_temp(:), &
      eig_vec(:), lanv(:), local_pw(:), pw(:), rmat(:), tmat(:), v0(:), &
      v_prev(:), vcur(:), w(:), xtmat(:,:), z(:,:), z_temp(:)
    REAL(real_8), EXTERNAL                   :: ddot

! ==================================================================

    rank     =  parai%me
    mem_flag = .TRUE.
    ! ==--------------------------------------------------------------==
    ! With LSD and 1 electron...
    IF (ndiag.EQ.0) RETURN
    ! ==--------------------------------------------------------------==
    CALL tiset('  DIST_FRIESNER',isub)
    ! ==--------------------------------------------------------------==
    ! TODO check stat
    ! TODO align for BG
    ALLOCATE(index0(ndiag),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(adav(ndiag, ndiag),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(alpha0(ndiag),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(beta1(ndiag),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    tim1=m_walltime()
    sign=-1._real_8
    IF (cntl%tlsd) sign=-2._real_8
    tlsd2=cntl%tlsd
    cntl%tlsd=.FALSE.
    tconv=.FALSE.
    tdebug=.FALSE.
    conv_flag=.FALSE.
    ngw2=ncpw%ngw*2
    nconv=0
    nconvold=0
    nhpsiold=0
    IF (fint1%ttrot) THEN
       trotmin=EXP(-edavmax*fint1%betap)
       trotmax=EXP(-edavmin*fint1%betap)
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (paral%parent.AND.tinfo) THEN
       IF (jspin.EQ.0) THEN
          charspin='<<<<<<<<<<<<<<<<<<<<'
       ELSEIF (jspin.EQ.1) THEN
          charspin='<<<<<<< ALPHA SPIN<<'
       ELSEIF (jspin.EQ.2) THEN
          charspin='<<<<<<<< BETA SPIN<<'
       ENDIF
       IF (paral%io_parent)&
            WRITE(6,'(/,1X,19("<"),A,A20)')&
            ' LANCZOS DIAGONALIZATION ',charspin
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ==================================================================
    ! DEFINE BLOCKS FOR DISTRIBUTED LANCZOS EIGENSOLVER
    ! ==================================================================
    CALL my_set_orbdist(ndiag,cnti%nstblk,norbx)
    norb=paraw%nwa12(parai%mepos,2)-paraw%nwa12(parai%mepos,1)+1
    norb=MAX(0,norb)
    sz = parai%nproc
    ! ==================================================================
    ! CREATE LANCZOS GROUP
    ! ==================================================================
    IF (paral%parent) THEN
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
    ! ==================================================================
    ! Subspace Matrix
    tcheck=.TRUE.
    ntest=0
    DO WHILE(tcheck)
       IF (fint1%ttrot) THEN
          CALL ehpsi(c0,cs,vpot,psi,ndiag)
       ELSE
          CALL hpsi(c0(:,1:ndiag),cs,sc0,vpot,psi,ndiag,1,ispin)
          CALL hfxpsi(cgs(:,1:nstate),c0(:,1:ndiag),cs(:,1:ndiag),focc,sign,psi,nstate,ndiag)
       ENDIF
       nhpsi=ndiag
       ! ================================================================
       ! DISTRIBUTED PROJECTED HAMILTONIAN IS CALCULATED HERE 
       ! ================================================================
       msglen=ndiag*norbx*8
       IF (mem_flag) THEN
          ALLOCATE(xtmat(ndiag,norbx),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(rmat(ndiag*norbx),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(tmat(ndiag*norbx),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
       DO ip=0,parai%nproc-1
          nx=paraw%nwa12(ip,2)-paraw%nwa12(ip,1)+1
          CALL zeroing(tmat)!,ndiag*norbx)
          IF (nx.GT.0) THEN
             CALL ovlap2(ncpw%ngw,ndiag,nx,tmat,c0,cs(1,paraw%nwa12(ip,1)),.TRUE.)
             CALL mp_sum(tmat,rmat,ndiag*norbx,parai%allgrp)
             IF (parai%me.EQ.ip) THEN
                CALL dcopy(ndiag*norbx,rmat(1),1,xtmat(1,1),1)
             ENDIF
          ENDIF
       ENDDO
       ! ================================================================
       ! WE NEED TO SYMMETRIZE OVERLAP MATRIX
       ! ================================================================
       ! CALL SYMM_DA(XTMAT,RMAT,TMAT,NDIAG,
       ! *              NWA12(0,1),NWA12(0,2),NORBX,MEPOS,NOLAN,LANGRP)
       ! call mp_sync(ALLGRP)
       ! ================================================================
       ! LANCZOS PARALLEL EIGENSOLVER: DIAGONALIZE XTMAT
       ! ONLY PROCESSORS THAT HAVE A NON-ZERO CHUNK OF THE RESTRICTED
       ! HAMILTONIAN IN XTMAT PARTICIPATE IN LANCZOS
       ! ================================================================
       IF ((paraw%nwa12(parai%me,2)-paraw%nwa12(parai%me,1)+1).NE.0) THEN
          sz = 0
          chunk_new = paraw%nwa12(parai%me,2)-paraw%nwa12(parai%me,1)+1
          DO ip=0, parai%nproc-1
             IF ((paraw%nwa12(ip,2)-paraw%nwa12(ip,1)+1).NE.0) sz = sz+1
          ENDDO
          nolan = sz
          lan_max = ndiag+1
          ! ==================
          ! ALLOCATE ONLY ONCE
          ! ==================
          IF (mem_flag) THEN
             ALLOCATE(vcur(ndiag),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             ALLOCATE(v0(ndiag),STAT=ierr)
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
          ENDIF

          CALL zeroing(pw)!,lan_max)
          CALL zeroing(lanv)!,chunk_new*lan_max)
          CALL zeroing(w)!,chunk_new)
          ! DETERMINE SEND AND RECEIVE BUFFERS
          msglen=8/2
          CALL my_concat(chunk_new, send_cnt(1), msglen, langrp)
          send_displ(1) = 0
          DO index4=2, sz
             send_displ(index4)=send_displ(index4-1)+send_cnt(index4-1)
          ENDDO
          ! DETERMINE A UNIT NORM INITIAL VECTOR V0
          IF (paral%parent) THEN
             CALL repprngu_vec(ndiag,v0) 
             nrm = dnrm2(ndiag, v0(1), 1)
             nrm = 1._real_8/nrm
             CALL dscal(ndiag, nrm, v0(1), 1)
          ENDIF
          ! BROADCAST INITIAL VECTOR TO ALL PROCESSORS IN LANCZOS GROUP
          CALL mp_bcast(v0, ndiag, parai%source, langrp)
          CALL dcopy(ndiag, v0(1), 1, vcur(1), 1)
          ! INITIALIZE V_PREV and BETA
          DO is=1, chunk_new
             v_prev(is) = 0._real_8
             w(is)=0._real_8
          ENDDO
          beta = 0._real_8
          chunk_begin = paraw%nwa12(parai%me,1)
          chunk_end   = paraw%nwa12(parai%me,2)
          ! INITIALIZE FIRST LANCZOS VECTOR
          CALL dcopy(chunk_new, vcur(chunk_begin), 1, lanv(1), 1)
          ! INITIALIZE INITIAL SUM OF WANTED EIGENVALUES
          sm_old = 0._real_8
          ! INITIALIZE CONVERGENCE TOLERANCE
          conv_tol = 1.e-8_real_8
          ! ==============================================================
          ! LANCZOS LOOP
          ! ==============================================================
          num_eig = ndiag
          DO index4 = 1, lan_max-1
             ! DISTRIBUTED MATRIX-VECTOR PRODUCT
             ! W: Holds the result
             CALL dgemv("T", ndiag,chunk_new,1._real_8,xtmat(1,1),ndiag,&
                  vcur(1), 1, 0._real_8, w(1), 1)
             ! DAXPY: W = W - BETA*V_(J-1)
             CALL daxpy(chunk_new, -beta, v_prev(1), 1, w(1), 1)
             ! CALCULATE ALPHA: INNER PRODUCT WITH CURRENT LANCZOS VECTOR
             local_alpha = ddot(chunk_new,w(1),1,&
                  lanv((index4-1)*chunk_new+1), 1)
             CALL mp_sum(local_alpha, alpha, langrp)
             ! DAXPY: W = W - ALPHA*V_J
             CALL daxpy(chunk_new,-alpha,lanv((index4-1)*chunk_new+1),1,&
                  w(1),1)
             ! ============================================================
             ! REORTHOGONALIZE AGAINST ALL PREVIOUS BASIS VECTORS
             ! ============================================================
             DO index1=1,index4
                local_pw(index1)=ddot(chunk_new,lanv((index1-1)*&
                     chunk_new+1),1, w(1), 1)
             ENDDO
             CALL mp_sum(local_pw, pw, index4, langrp)
             DO index1=1, index4
                CALL daxpy(chunk_new, -pw(index1),&
                     lanv((index1-1)*chunk_new+1), 1, w(1), 1)
             ENDDO
             ! ============================================================
             ! CALCULATE BETA: NORM OF VECTOR W
             local_beta = ddot(chunk_new, w(1), 1, w(1), 1)
             CALL mp_sum(local_beta, beta, langrp)
             beta = SQRT(beta)
             ! NORMALIZE NEXT LANCZOS VECTOR
             CALL dscal(chunk_new, 1/beta, w(1), 1)
             CALL dcopy(chunk_new, w(1), 1, lanv(index4*chunk_new+1), 1)
             d(index4) = alpha
             e(index4) = beta
             ! UPDATE AUXILIARY VECTORS
             CALL my_concatv(lanv((index4)*chunk_new+1),vcur(1),&
                  chunk_new,send_cnt(1),send_displ(1),langrp)
             CALL dcopy(chunk_new,lanv((index4-1)*chunk_new+1),1,&
                  v_prev(1),1)
             ! CHECK CONVERGENCE OF EIGENVALUES: OBSOLETE (FUTURE USE ONLY)
             CALL dcopy(index4, d(1), 1, ld(1), 1)
             CALL dcopy(index4-1, e(1), 1, le(1), 1)
             abtol = -1
          ENDDO
          ! =============================================================
          ! END OF LANCZOS LOOP
          ! =============================================================
          CALL mp_sync(langrp)
          IF (.NOT.conv_flag) index4 = index4-1
          chunk_begin_e = paraw%nwa12(parai%mepos,1)
          chunk_end_e   = paraw%nwa12(parai%mepos,2)
          chunk_new_e   = chunk_end_e - chunk_begin_e +1
          ! =======================
          ! DETERMINE SEND  BUFFERS
          ! =======================
          msglen=8/2
          CALL my_concat(chunk_new_e, send_cnt(1), msglen, langrp)
          send_displ(1) = 0
          DO index1=2, sz
             send_displ(index1) = send_displ(index1-1) +&
                  send_cnt(index1-1)
          ENDDO
          ! ALLOCATE MEMORY FOR EIGENVECTORS
          IF (mem_flag) THEN
             ALLOCATE(z(ndiag,(chunk_end_e-chunk_begin_e+1)*index4/ndiag),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
          ENDIF
          CALL  dstevx('V','I',index4,d(1),e(1),vl,vu,chunk_begin_e,&
               chunk_end_e, abtol, nfound, teig(1), z(1,1),&
               index4, work(1), iwork(1), ifail(1), ierr)

          CALL mp_sync(langrp)
          ! ==============================================================
          ! CALCULATION OF APPROXIMATE EIGENVECTORS
          ! EACH PROCESSOR WILL STORE COMPLETE EIGENVECTORS
          ! EIGENVECTORS ARE DISTRIBUTED ACROSS PROCESSORS
          ! ==============================================================
          IF (mem_flag) THEN
             ALLOCATE(z_temp(ndiag*(chunk_end_e-chunk_begin_e+1)*index4/ndiag),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
          ENDIF
          ! ALLOCATION OF MEMORY FOR APPROXIMATE EIGENVECTORS
          IF (mem_flag) THEN
             ALLOCATE(eig_vec(ndiag*chunk_new_e),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
          ENDIF
          ! ALLOCATE ADDITIONAL MEMORY FOR BOOK KEEPING
          msglen=8/2
          CALL my_concat(chunk_new, old_send_cnt_e(1), msglen, langrp)
          IF (mem_flag) THEN
             ALLOCATE(e_temp(2*old_send_cnt_e(parai%me+1)*&
                  MAX(chunk_new_e,send_cnt(1))),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
          ENDIF
          DO index1=1, nolan
             IF (rank.EQ.(index1-1) ) THEN
                CALL dcopy(index4*chunk_new_e, z(1,1), 1, z_temp(1), 1)
                DO index2 = 1, sz
                   send_cnt_e(index2) = old_send_cnt_e(index2) *&
                        send_cnt(rank+1)
                ENDDO
                send_displ_e(1) = 0
                DO index2=2, sz
                   send_displ_e(index2) = send_displ_e(index2-1) +&
                        send_cnt_e(index2-1)
                ENDDO
             ENDIF
             ! BROADCAST COUNTS AND DISPL. DATA FOR ALLGATHER
             CALL mp_bcast(send_cnt_e, sz, index1-1, langrp)
             CALL mp_bcast(send_displ_e, sz, index1-1, langrp)
             ! BROADCAST MY EIGENVECTORS OF THE TRIDIAGONAL MATRIX
             CALL mp_bcast(z_temp, send_cnt(index1)*index4,&
                  index1-1,langrp)
             ! MULTIPLY RECEIVED EIGENVECTORS WITH MY PART OF THE LANCZOS
             ! BASIS
             CALL dgemm("N","N",chunk_new,send_cnt(index1),index4,1._real_8,&
                  lanv(1),chunk_new,z_temp(1),index4,&
                  0._real_8, e_temp(1), chunk_new)
             ! BRING RESULTS BACK TO THE CORESSPONDING PROCESSOR
             CALL my_source_concatv(e_temp(1),eig_vec(1),&
                  chunk_new*send_cnt(index1),&
                  send_cnt_e(1), send_displ_e(1),&
                  index1-1, langrp)
          ENDDO
          ! ==================================
          ! SET CORRECT ORDER FOR EIGENVECTORS
          ! ==================================
          !$omp parallel do private(INDEX1)
          DO index1 = 1, ndiag
             lc_index(index1) = 0
          ENDDO
          start     = 1
          start_col = 1
          DO index1 = 1, sz
             DO index2 = 1, chunk_new_e
                CALL dcopy(old_send_cnt_e(index1), eig_vec(start_col), 1,&
                     z_temp(start+lc_index(1)), 1)
                start_col = start_col + old_send_cnt_e(index1)
                start     = start + ndiag
             ENDDO
             lc_index(1) = lc_index(1) + old_send_cnt_e(index1)
             start = 1
          ENDDO
          CALL dcopy(ndiag*chunk_new_e, z_temp(1), 1, eig_vec(1), 1)
          ! ======================
          ! NORMALIZE EIGENVECTORS
          ! ======================
          DO index1=1, chunk_new
             nrm = dnrm2(ndiag, eig_vec((index1-1)*ndiag+1), 1)
             nrm = 1/SQRT(nrm)
             CALL dscal(ndiag, nrm, eig_vec((index1-1)*ndiag+1), 1)
          ENDDO
          CALL mp_sync(langrp)
          ! ==============================================================
          ! END LANCZOS EIGENSOLVER ON LANGRP
          ! ==============================================================
       ELSE
          ! PROCESSORS THAT DID NOT PARTICIPATE IN LANCZOS
          ! WILL WAIT IN SUBSEQUENT COLLECTIVE
       ENDIF
       ! ==============================================================
       ! ALL PROCESSORS GET ALL EIGENVALUES
       ! ==============================================================
       IF (chunk_new.NE.0) THEN
          CALL my_concatv(teig(1),edav(1),chunk_new,&
               send_cnt(1),send_displ(1),langrp)
       ENDIF
       CALL mp_bcast(edav, ndiag, parai%source, parai%allgrp)
       ! ==============================================================
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
          IF (mem_flag) THEN
             ALLOCATE(z_temp(ndiag*(chunk_end_e-chunk_begin_e+1)*index4/ndiag),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
          ENDIF
       ENDIF
       ! ====================
       ! ROTATE WAVEFUNCTIONS
       ! ====================
       start_col = 0
       start = 1
       DO index1=1, nolan
          IF (rank.EQ.(index1-1) ) THEN
             ! INITIALIZE TEMP. SPACE TO ZERO
             CALL zeroing(z_temp(1:ndiag*chunk_new_e))!, ndiag*chunk_new_e)
             CALL dcopy(ndiag*chunk_new_e, eig_vec(1), 1, z_temp(1), 1)
          ENDIF

          ! CURRENT ROOT BROADCASTS ITS APPROX. EIGENVECTORS
          CALL mp_bcast(z_temp, send_cnt(index1)*ndiag, index1-1,&
               parai%allgrp)
          ! EACH PROCESSOR MULTIPLIES THE APPROXIMATE EIGVECS
          ! RECEIVED AND STORES THEM AT THE APPROPRIATE POSITION
          ! IN MATRIX CSCR
          CALL dgemm("N","N",2*ncpw%ngw,send_cnt(index1),ndiag,1._real_8,&
               c0(1,1),2*ncpw%ngw,z_temp(1),ndiag,&
               0._real_8, cscr(1,start), 2*ncpw%ngw)
          start = start + send_cnt(index1)
       ENDDO
       CALL dcopy(ngw2*ndiag,cscr(1,1),1,c0(1,1),1)
       ! Check if an eigenvalue is not equal to 0 (dim < NDIAG)
       tcheck=.FALSE.
       itest=0
       DO i=1,ndiag
          IF (fint1%ttrot) THEN
             ttest = (edav(i).GT.trotmax).OR.(edav(i).LT.trotmin)
          ELSE
             ttest = ABS(edav(i)).LT.eps
          ENDIF
          IF (ttest) THEN
             tcheck=.TRUE.
             IF (paral%io_parent)&
                  WRITE(6,*) 'DIST_FRIESNER| EIGENVECTOR',i,&
                  ' IS VERY BAD! '
             ! Not useful: Now we use H|C0>.
             ! CALL RANDTOWF(C0(1,I),1,0,0)
             itest=itest+1
             edav(i)=1.e33_real_8
          ENDIF
       ENDDO
       IF (tcheck) THEN
          ntest=ntest+1
          IF (ntest.GT.10.AND.paral%parent)&
               CALL stopgm(' DIST_FRIESNER', 'CAN''T FIND GUEST VECTORS',& 
               __LINE__,__FILE__)
          CALL sort2(edav,ndiag,index0)
          ! Put at the end randomized wavefunctions.
          CALL dcopy(ngw2*ndiag,c0,1,cs,1)
          DO i=1,ndiag
             j=index0(i)
             IF (edav(i).EQ.big) THEN
                CALL dcopy(ngw2,cs(1,j),1,c0(1,i),1)
             ELSE
                CALL dcopy(ngw2,cscr(1,j),1,c0(1,i),1)
             ENDIF
          ENDDO
          ! Orthogonalize.
          CALL gs_ortho(c0,ndiag-itest,c0(1,ndiag-itest+1),itest)
       ENDIF
       ! SET MEMORY ALLOCATION FLAG TO FALSE
       mem_flag = .FALSE.
    ENDDO
    ! ==================================================================
    CALL mp_sync(parai%allgrp)
    ! ====================
    ! ROTATE WAVEFUNCTIONS
    ! ====================
    start_col = 0
    start = 1
    DO index1=1, nolan
       IF (rank.EQ.(index1-1) ) THEN
          ! INITIALIZE TEMP. SPACE TO ZERO
          CALL zeroing(z_temp(1:ndiag*chunk_new_e))!, ndiag*chunk_new_e)
          CALL dcopy(ndiag*chunk_new_e, eig_vec(1), 1, z_temp(1), 1)
       ENDIF
       ! CURRENT ROOT BROADCASTS ITS APPROX. EIGENVECTORS
       CALL mp_bcast(z_temp, send_cnt(index1)*ndiag, index1-1,&
            parai%allgrp)
       ! EACH PROCESSOR MULTIPLIES THE APPROXIMATE EIGVECS
       ! RECEIVED AND STORES THEM AT THE APPROPRIATE POSITION
       ! IN MATRIX CSCR
       CALL dgemm("N","N",2*ncpw%ngw,send_cnt(index1),ndiag,1._real_8,&
            cs(1,1),2*ncpw%ngw,z_temp(1),ndiag,&
            0._real_8, cscr(1,start), 2*ncpw%ngw)
       start = start + send_cnt(index1)
    ENDDO
    CALL dcopy(ngw2*ndiag,cscr(1,1),1,cs(1,1),1)
    tim2=m_walltime()
    tim=(tim2-tim1)*0.001_real_8
    IF (paral%parent.AND.tinfo) THEN
       IF (paral%io_parent)&
            WRITE(6,'(A,T58,F8.2)')&
            ' >> TIME FOR INITIAL SUBSPACE DIAGONALIZATION:  ',tim
       IF (paral%io_parent)&
            WRITE(6,'(A,8X,A,8X,A)') ' >> CYCLE     NCONV',&
            'B2MAX','B2MIN     #HPSI      TIME'
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL mp_sync(parai%allgrp)
    CALL zeroing(index0)!,ndiag)
    DO icycle=1,cnti%n_fries
       tim1=m_walltime()
       b2max=0._real_8
       b2min=1.e30_real_8
       ! We calculate again the residue for all vectors
       nconv=0
       ncurr=nconv+1
       ! First Lanczos Step
       DO i=ncurr,ndiag
          !$omp parallel do private(IG)
          DO ig=1,ncpw%ngw
             cscr(ig,1)=cs(ig,i)-edav(i)*c0(ig,i)
          ENDDO
          b2=dotp(ncpw%ngw,cscr(:,1),cscr(:,1))
          CALL mp_sum(b2,parai%allgrp)
          beta1(i)=SQRT(b2)
          alpha0(i)=edav(i)
          ! IF (PARENT) WRITE(6,*) "RESIDUAL !!!", I, BETA, B2LIMIT
          IF (b2.LT.cntr%b2limit) THEN
             nconv=nconv+1
             index0(i)=nconv
             IF (b2.GT.b2max) b2max=b2
             IF (b2.LT.b2min) b2min=b2
          ELSE
             CALL dcopy(ngw2,cscr(1,1),1,cs(1,i),1)
             fac=1._real_8/beta1(i)
             CALL dscal(ngw2,fac,cs(1,i),1)
          ENDIF
       ENDDO
       ! Order states: converged first
       DO i=ncurr,ndiag
          j=index0(i)
          IF (j.NE.0.AND.j.NE.i) THEN
             CALL dswap(ngw2,c0(1,i),1,c0(1,j),1)
             CALL dswap(ngw2,cs(1,i),1,cs(1,j),1)
             index0(j)=j
             index0(i)=0
             aux = beta1(i)
             beta1(i) = beta1(j)
             beta1(j) = aux
             aux = alpha0(i)
             alpha0(i) = alpha0(j)
             alpha0(j) = aux
          ENDIF
       ENDDO
       ncurr=nconv+1
       ! Lanczos Refinement
       IF (icycle.EQ.1.AND.ncurr.GT.ndiag) THEN
          trefine=.FALSE.
       ELSE
          trefine=.TRUE.
       ENDIF
       CALL krylov_ref(ncurr,ndiag,nconv,nhpsi,nstate,&
            c0,cs,sc0,cscr,cgs,focc,sign,&
            vpot,psi,index0,alpha0,beta1,b2min,b2max)
       ! Order states : converged first
       DO i=ncurr,ndiag
          j=index0(i)
          IF (j.NE.0.AND.j.NE.i) THEN
             CALL dswap(ngw2,c0(1,i),1,c0(1,j),1)
             CALL dswap(ngw2,cs(1,i),1,cs(1,j),1)
             index0(j)=j
             index0(i)=0
             aux = beta1(i)
             beta1(i) = beta1(j)
             beta1(j) = aux
          ENDIF
       ENDDO
       ! Reorthogonalize
       IF (ncurr.LE.ndiag) THEN
          CALL gs_ortho(c0,ncurr-1,c0(1,ncurr),ndiag-ncurr+1)
       ENDIF
       ! Calculate new forces; only for states that entered Lanczos
       nleft=ndiag-ncurr+1

       IF (ncurr.LE.ndiag) THEN
          IF (fint1%ttrot) THEN
             CALL ehpsi(c0(1,ncurr),cs(1,ncurr),vpot,psi,nleft)
          ELSE
             CALL hpsi(c0(:,ncurr:ncurr+nleft-1),cs(1,ncurr),sc0(1,ncurr),&
                  vpot,psi,nleft,1,ispin)
             CALL hfxpsi(cgs(:,1:nstate),c0(:,ncurr:ncurr+nleft-1),cs(:,ncurr:ncurr+nleft-1),focc,sign,&
                  psi,nstate,nleft)
          ENDIF
          nhpsi=nhpsi+nleft
       ENDIF
       ! ===============================================================
       ! DISTRIBUTED PROJECTED HAMILTONIAN IS CALCULATED HERE 
       ! ===============================================================
       DO ip=0,parai%nproc-1
          nx=paraw%nwa12(ip,2)-paraw%nwa12(ip,1)+1
          CALL zeroing(tmat)!,ndiag*norbx)
          IF (nx.GT.0) THEN
             CALL ovlap2(ncpw%ngw,ndiag,nx,tmat,c0,cs(1,paraw%nwa12(ip,1)),.TRUE.)
             CALL mp_sum(tmat,rmat,ndiag*norbx,parai%allgrp)
             IF (parai%me.EQ.ip) THEN
                CALL dcopy(ndiag*norbx,rmat(1),1,xtmat(1,1),1)
             ENDIF
          ENDIF
       ENDDO
       ! ================================================================ 
       ! ====================================
       ! WE NEED TO SYMMETRIZE OVERLAP MATRIX
       ! ====================================
       ! CALL SYMM_DA(XTMAT,RMAT,TMAT,NDIAG,
       ! *              NWA12(0,1),NWA12(0,2),NORBX,MEPOS,NPROC,LANGRP)
       ! 
       ! ================================================================
       ! LANCZOS PARALLEL EIGENSOLVER
       ! ================================================================
       IF ((paraw%nwa12(parai%me,2)-paraw%nwa12(parai%me,1)+1).NE.0) THEN
          sz = 0
          chunk_new = paraw%nwa12(parai%me,2)-paraw%nwa12(parai%me,1)+1
          DO ip=0, parai%nproc-1
             IF ((paraw%nwa12(ip,2)-paraw%nwa12(ip,1)+1).NE.0) sz = sz+1
          ENDDO
          nolan = sz
          lan_max = ndiag+1
          CALL zeroing(pw)!,lan_max)
          CALL zeroing(lanv)!,chunk_new*lan_max)
          CALL zeroing(w)!,chunk_new)
          ! DETERMINE SEND AND RECEIVE BUFFERS
          msglen=8/2
          CALL my_concat(chunk_new, send_cnt(1), msglen, langrp)
          send_displ(1) = 0
          DO index4=2, sz
             send_displ(index4) = send_displ(index4-1) +&
                  send_cnt(index4-1)
          ENDDO
          ! DETERMINE A UNIT NORM INITIAL VECTOR V0
          IF (paral%parent) THEN
             CALL repprngu_vec(ndiag,v0) 
             nrm = dnrm2(ndiag, v0(1), 1)
             nrm = 1/nrm
             CALL dscal(ndiag, nrm, v0(1), 1)
          ENDIF
          ! BROADCAST INITIAL VECTOR TO ALL PROCESSORS
          CALL mp_bcast(v0, ndiag, parai%source, langrp)
          CALL dcopy(ndiag, v0(1), 1, vcur(1), 1)
          ! INITIALIZE V_PREV and BETA
          DO is=1, chunk_new
             v_prev(is) = 0._real_8
             w(is)=0._real_8
          ENDDO
          beta = 0._real_8
          chunk_begin = paraw%nwa12(parai%me,1)
          chunk_end   = paraw%nwa12(parai%me,2)
          ! INITIALIZE FIRST LANCZOS VECTOR
          CALL dcopy(chunk_new, vcur(chunk_begin), 1, lanv(1), 1)
          ! INITIALIZE INITIAL SUM OF WANTED EIGENVALUES
          sm_old = 0._real_8
          ! INITIALIZE CONVERGENCE TOLERANCE
          conv_tol = 1.e-8_real_8
          ! ============
          ! LANCZOS LOOP
          ! ============
          num_eig = ndiag
          DO index4 = 1, lan_max-1
             ! DISTRIBUTED MATRIX-VECTOR PRODUCT
             ! W: Holds the result
             CALL dgemv("T", ndiag,chunk_new,1._real_8,xtmat(1,1),ndiag,&
                  vcur(1), 1, 0._real_8, w(1), 1)
             ! DAXPY: W = W - BETA*V_(J-1)
             CALL daxpy(chunk_new, -beta, v_prev(1), 1, w(1), 1)
             ! CALCULATE ALPHA: INNER PRODUCT WITH CURRENT LANCZOS VECTOR
             local_alpha = ddot(chunk_new,w(1),1,&
                  lanv((index4-1)*chunk_new+1), 1)
             CALL mp_sum(local_alpha, alpha, langrp)
             ! DAXPY: W = W - ALPHA*V_J
             CALL daxpy(chunk_new,-alpha,lanv((index4-1)*chunk_new+1),1,&
                  w(1),1)
             ! ============================================================
             ! REORTHOGONALIZE AGAINST ALL PREVIOUS BASIS VECTORS
             ! ============================================================
             DO index1=1,index4
                local_pw(index1)=ddot(chunk_new,lanv((index1-1)*&
                     chunk_new+1),&
                     1, w(1), 1)
             ENDDO

             CALL mp_sum(local_pw, pw, index4, langrp)

             DO index1=1, index4
                CALL daxpy(chunk_new, -pw(index1),&
                     lanv((index1-1)*chunk_new+1), 1, w(1), 1)
             ENDDO
             ! ============================================================
             ! CALCULATE BETA: NORM OF VECTOR W
             local_beta = ddot(chunk_new, w(1), 1, w(1), 1)
             CALL mp_sum(local_beta, beta, langrp)
             beta = SQRT(beta)
             ! NORMALIZE NEXT LANCZOS VECTOR
             CALL dscal(chunk_new, 1/beta, w(1), 1)
             CALL dcopy(chunk_new, w(1), 1, lanv(index4*chunk_new+1), 1)
             d(index4) = alpha
             e(index4) = beta

             ! UPDATE AUXILIARY VECTORS
             CALL my_concatv(lanv((index4)*chunk_new+1),vcur(1),&
                  chunk_new, send_cnt(1),send_displ(1),langrp)
             CALL dcopy(chunk_new,lanv((index4-1)*chunk_new+1),1,&
                  v_prev(1),1)
             ! CHECK CONVERGENCE FOR EIGENVALUES: FUTURE USE ONLY
             CALL dcopy(index4, d(1), 1, ld(1), 1)
             CALL dcopy(index4-1, e(1), 1, le(1), 1)
             abtol = -1
          ENDDO
          ! ===================
          ! END OF LANCZOS LOOP
          ! ===================
          IF (.NOT.conv_flag) index4 = index4-1
          chunk_begin_e = paraw%nwa12(parai%mepos,1)
          chunk_end_e   = paraw%nwa12(parai%mepos,2)
          chunk_new_e   = chunk_end_e - chunk_begin_e +1
          ! =======================
          ! DETERMINE SEND  BUFFERS
          ! =======================
          msglen=8/2
          CALL my_concat(chunk_new_e, send_cnt(1), msglen, langrp)
          send_displ(1) = 0
          DO index1=2, sz
             send_displ(index1) = send_displ(index1-1) + send_cnt(index1-1)
          ENDDO
          CALL  dstevx('V','I',index4,d(1),e(1),vl,vu,chunk_begin_e,&
               chunk_end_e,&
               abtol, nfound, teig(1), z(1,1), index4, work(1),&
               iwork(1), ifail(1), ierr)
          ! ================================================================
          ! CALCULATION OF APPROXIMATE EIGENVECTORS
          ! EACH PROCESSOR WILL STORE COMPLETE EIGENVECTORS
          ! EIGENVECTORS ARE DISTRIBUTED ACROSS PROCESSORS
          ! ================================================================
          ! ALLOCATE ADDITIONAL MEMORY FOR BOOK KEEPING
          msglen=8/2
          CALL my_concat(chunk_new, old_send_cnt_e(1), msglen, langrp)
          DO index1=1, nolan
             IF (rank.EQ.(index1-1) ) THEN
                CALL dcopy(index4*chunk_new_e, z(1,1), 1, z_temp(1), 1)

                DO index2 = 1, sz
                   send_cnt_e(index2) = old_send_cnt_e(index2) *&
                        send_cnt(rank+1)
                ENDDO
                send_displ_e(1) = 0
                DO index2=2, sz
                   send_displ_e(index2) = send_displ_e(index2-1) +&
                        send_cnt_e(index2-1)
                ENDDO
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
                  0._real_8, e_temp(1), chunk_new)
             ! BRING RESULTS BACK TO THE CORESSPONDING PROCESSOR
             CALL my_source_concatv(e_temp(1),eig_vec(1),&
                  chunk_new*send_cnt(index1),&
                  send_cnt_e(1), send_displ_e(1),&
                  index1-1, langrp)
          ENDDO
          ! ================================================================
          DO index1 = 1, ndiag
             lc_index(index1) = 0
          ENDDO
          start     = 1
          start_col = 1
          DO index1 = 1, sz
             DO index2 = 1, chunk_new_e
                CALL dcopy(old_send_cnt_e(index1), eig_vec(start_col), 1,&
                     z_temp(start+lc_index(1)), 1)
                start_col = start_col + old_send_cnt_e(index1)
                start     = start + ndiag
             ENDDO
             lc_index(1) = lc_index(1) + old_send_cnt_e(index1)
             start = 1
          ENDDO
          CALL dcopy(ndiag*chunk_new_e, z_temp(1), 1, eig_vec(1), 1)
          ! ======================
          ! NORMALIZE EIGENVECTORS
          ! ======================
          DO index1=1, chunk_new
             nrm = dnrm2(ndiag, eig_vec((index1-1)*ndiag+1), 1)
             nrm = 1/SQRT(nrm)
             CALL dscal(ndiag, nrm, eig_vec((index1-1)*ndiag+1), 1)
          ENDDO
          CALL mp_sync(langrp)
          ! ================================================================
          ! END LANCZOS ITERATION ON LANGRP
          ! ================================================================
       ELSE
          ! PROCESSORS THAT DID NOT PARTICIPATE IN LANCZOS
          ! WILL WAIT IN SUBSEQUENT COLLECTIVE
       ENDIF
       ! ==============================================================
       ! ALL PROCESSORS GET ALL EIGENVALUES
       ! ==============================================================
       IF (chunk_new.NE.0) THEN
          CALL my_concatv(teig(1),edav(1),chunk_new,&
               send_cnt(1),send_displ(1),langrp)
       ENDIF
       CALL mp_bcast(edav, ndiag, parai%source, parai%allgrp)
       ! ==============================================================
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
       ENDIF
       ! ====================
       ! ROTATE WAVEFUNCTIONS
       ! ====================
       start_col = 0
       start = 1
       DO index1=1, nolan
          IF (rank.EQ.(index1-1) ) THEN
             ! INITIALIZE TEMP. SPACE TO ZERO
             CALL zeroing(z_temp(1:ndiag*chunk_new_e))!, ndiag*chunk_new_e)
             CALL dcopy(ndiag*chunk_new_e, eig_vec(1), 1, z_temp(1), 1)
          ENDIF
          ! CURRENT ROOT BROADCASTS ITS APPROX. EIGENVECTORS
          CALL mp_bcast(z_temp, send_cnt(index1)*ndiag, index1-1,&
               parai%allgrp)
          ! EACH PROCESSOR MULTIPLIES THE APPROXIMATE EIGVECS
          ! RECEIVED AND STORES THEM AT THE APPROPRIATE POSITION
          ! IN MATRIX CSCR
          CALL dgemm("N","N",2*ncpw%ngw,send_cnt(index1),ndiag,1._real_8,&
               c0(1,1),2*ncpw%ngw,z_temp(1),ndiag,&
               0._real_8, cscr(1,start), 2*ncpw%ngw)
          start = start + send_cnt(index1)
       ENDDO
       CALL dcopy(ngw2*ndiag,cscr(1,1),1,c0(1,1),1)
       start_col = 0
       start = 1
       DO index1=1, nolan
          IF (rank.EQ.(index1-1) ) THEN
             ! INITIALIZE TEMP. SPACE TO ZERO
             CALL zeroing(z_temp(1:ndiag*chunk_new_e))!, ndiag*chunk_new_e)
             CALL dcopy(ndiag*chunk_new_e, eig_vec(1), 1, z_temp(1), 1)
          ENDIF
          ! CURRENT ROOT BROADCASTS ITS APPROX. EIGENVECTORS
          CALL mp_bcast(z_temp, send_cnt(index1)*ndiag, index1-1,&
               parai%allgrp)
          ! EACH PROCESSOR MULTIPLIES THE APPROXIMATE EIGVECS
          ! RECEIVED AND STORES THEM AT THE APPROPRIATE POSITION
          ! IN MATRIX CSCR
          CALL dgemm("N","N",2*ncpw%ngw,send_cnt(index1),ndiag,1._real_8,&
               cs(1,1),2*ncpw%ngw,z_temp(1),ndiag,&
               0._real_8, cscr(1,start), 2*ncpw%ngw)
          start = start + send_cnt(index1)
       ENDDO
       CALL dcopy(ngw2*ndiag,cscr(1,1),1,cs(1,1),1)
       ! If one eigenvalue is equal to zero (like friesner_c).
       cnti%ntrans=nconv
       DO i=1,nconv
          IF (fint1%ttrot) THEN
             ttest = (edav(i).GT.trotmax).OR.(edav(i).LT.trotmin)
          ELSE
             ttest = ABS(edav(i)).LT.eps
          ENDIF
          IF (ttest) THEN
             nconv=nconv-1
             cnti%ntrans=cnti%ntrans+1
             IF (paral%io_parent)&
                  WRITE(6,*) ' DIST_FRIESNER| ','VERY BAD EIGENVALUE !!'
             IF (cnti%ntrans.GT.ndiag) THEN
                CALL randtowf(c0(:,i:i),1,0,0)
             ELSE
                CALL dswap(ngw2,c0(1,i),1,c0(1,cnti%ntrans),1)
                CALL dswap(ngw2,cs(1,i),1,cs(1,cnti%ntrans),1)
             ENDIF
          ENDIF
       ENDDO
       tim2=m_walltime()
       tim=(tim2-tim1)*0.001_real_8
       IF (paral%parent.AND.tinfo) THEN
          IF (paral%io_parent)&
               WRITE(6,'(4X,I5,5X,I5,2(1PE13.3),0PF10.2,0PF10.2)')&
               icycle,nconv,b2max,b2min,&
               REAL(nhpsi-nhpsiold,kind=real_8)/REAL(ndiag,kind=real_8),tim
       ENDIF
       nhpsiold=nhpsi
       ! Test if the number of eigenvalues is enough to stop.
       IF ((.NOT.fint1%tfral).AND.cntl%tfint.AND.&
            (nconv.NE.nconvold).AND.&
            (nconv.NE.ndiag)) THEN
          ! Spin or not.
          IF (tlsd2.AND.(nconv.GE.nel).OR.&
               (.NOT.tlsd2.AND.(nconv.GE.(nel+1)/2)) ) THEN
             tconv=.TRUE.
             nconvold=nconv
             CALL convfrie(ndiag,nconv,nel,edav,cntl%tlsd,amu,tconv)
          ENDIF
       ENDIF
       IF (tconv.OR.nconv.EQ.ndiag) THEN
          GOTO 200
       ENDIF
       ! End of Krylov Space Cycles
    ENDDO
200 CONTINUE
    IF (paral%parent.AND.((.NOT.tconv).AND.nconv.NE.ndiag).AND.tdebug) THEN
       IF (paral%io_parent)&
            WRITE(6,'(1X,64("!"))')
       IF (paral%io_parent)&
            WRITE(6,'(" !!",A,T64,"!!")')&
            ' DIST_FRIESNER| NOT ALL ROOTS ARE CONVERGED'
       IF (paral%io_parent)&
            WRITE(6,'(1X,64("!"))')
    ENDIF
    ! Order Roots with respect to energy
    CALL sort2(edav,ndiag,index0)
    CALL dcopy(ngw2*ndiag,c0(1,1),1,cs(1,1),1)
    IF (.NOT.fint1%ttrot) CALL dscal(ndiag,-1._real_8,edav(1),1)
    DO i = 1 , ndiag/2
       aux = edav(ndiag-i+1)
       edav(ndiag-i+1) = edav(i)
       edav(i) = aux
       iaux = index0(ndiag-i+1)
       index0(ndiag-i+1) = index0(i)
       index0(i) = iaux
    ENDDO
    DO i=1,ndiag
       j=index0(i)
       CALL dcopy(2*ncpw%ngw,cs(1,j),1,c0(1,i),1)
    ENDDO
    cntl%tlsd=tlsd2
    ! =================================================================
    ! FREE ALLOCATED MEMORY
    ! =================================================================
    IF (chunk_new.NE.0) THEN
       DEALLOCATE(v_prev,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(local_pw,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(pw,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(w,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(lanv,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(v0,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(vcur,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)

       DEALLOCATE(e_temp,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       ! C$$$         CALL FREEM(IP_EIG_VEC)
       DEALLOCATE(eig_vec,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(z_temp,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(z,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ELSE
       DEALLOCATE(z_temp,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! GLOBAL DEALLOCATIONS
    DEALLOCATE(tmat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rmat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(xtmat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    CALL tihalt('  DIST_FRIESNER',isub)
    ! ==--------------------------------------------------------------==
    DEALLOCATE(index0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(adav,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(alpha0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(beta1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)

    RETURN
  END SUBROUTINE dist_friesner
  ! ==================================================================
  SUBROUTINE give_scr_dist_friesner(lfriesner,tag,ndiag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lfriesner
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: ndiag

    INTEGER                                  :: ldiag, lhpsi, lkrylov_ref, &
                                                lscr, lsummat

    lscr=0
    ldiag=0
    lkrylov_ref=0
    lsummat=0
    lhpsi=0
    IF (.NOT.fint1%ttrot) THEN
       CALL give_scr_hpsi(lhpsi,tag,ndiag)
    ENDIF
    IF (.FALSE.) THEN
       CALL give_scr_summat(lsummat,tag,ndiag)
    ELSE
       lsummat = 0
    ENDIF
    CALL give_scr_krylov_ref(lkrylov_ref,tag,ndiag)
    IF (tkpts%tkpnt) THEN
       ldiag=((3*ndiag-2) + 5*(2*ndiag))! ZHEEV
       lscr=ndiag/2+1+2*ndiag*ndiag+2*ndiag! INDEX ADAV ALPHA0 BETA1
    ELSE
       ldiag=5*ndiag         ! DSYEV
       IF (.FALSE.) THEN
          lscr=ndiag/2+1+ndiag*ndiag+2*ndiag
       ELSE
          lscr = ndiag/2+1+ndiag+2*ndiag
       ENDIF                 ! INDEX ADAV ALPHA0 BETA1
    ENDIF
    ldiag=MAX(ldiag,2*ndiag)  ! CONVFRIE
    lfriesner=lscr+MAX(lhpsi,lsummat,lkrylov_ref,ldiag)+10
    tag='LSCR+MAX(HPSI,SUMHMAT,...)'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_dist_friesner
  ! ==================================================================

END MODULE dist_friesner_utils
