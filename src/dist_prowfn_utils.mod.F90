MODULE dist_prowfn_utils
  USE adat,                            ONLY: elem
  USE atwf,                            ONLY: atwf_mod,&
                                             atwp,&
                                             catom,&
                                             loadc_foc_array_size,&
                                             xsmat,&
                                             xxmat
  USE augchg_utils,                    ONLY: augchg
  USE cmaos_utils,                     ONLY: cmaos,&
                                             give_scr_cmaos,&
                                             give_scr_satch,&
                                             satch,&
                                             satch3,&
                                             satch4
  USE dotp_utils,                      ONLY: dotp
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_def,&
                                             fo_ufo
  USE forcep_utils,                    ONLY: rhoe_psi_size
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE jrotation_utils,                 ONLY: my_set_orbdist
  USE kinds,                           ONLY: real_8
  USE linalg_utils,                    ONLY: symm_da
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_group,&
                                             mp_recv,&
                                             mp_send,&
                                             mp_sum,&
                                             mp_sync
  USE parac,                           ONLY: parai,&
                                             paral
  USE prden,                           ONLY: numpr
  USE prmem_utils,                     ONLY: prmem
  USE prng_utils,                      ONLY: repprngu_vec
  USE prop,                            ONLY: prop1,&
                                             prop2,&
                                             prop3
  USE prowfn_utils,                    ONLY: mklabel
  USE pslo,                            ONLY: pslo_com
  USE rnlsm_utils,                     ONLY: give_scr_rnlsm,&
                                             rnlsm
  USE setbasis_utils,                  ONLY: loadc
  USE sfac,                            ONLY: fnl,&
                                             fnl2
  USE spin,                            ONLY: spin_mod
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             maxsys,&
                                             ncpw,&
                                             paraw
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: dist_prowfn
  PUBLIC :: dist_prtmat
  PUBLIC :: give_scr_dist_prowfn

CONTAINS

  ! ==================================================================
  SUBROUTINE dist_prowfn(c0,tau0,nstate)
    ! ==--------------------------------------------------------------==
    ! == PROJECT WAVEFUNCTIONS ON PSEUDO ATOMIC ORBITALS              ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(ncpw%ngw,*)
    REAL(real_8)                             :: tau0(*)
    INTEGER                                  :: nstate

    CHARACTER(*), PARAMETER                  :: procedureN = 'dist_prowfn'

    CHARACTER(len=15)                        :: label(1000)
    CHARACTER(len=30)                        :: tag
    COMPLEX(real_8), ALLOCATABLE             :: cscr(:,:), psi(:,:)
    INTEGER :: chunk_begin, chunk_begin_e, chunk_end, chunk_end_e, chunk_new, &
      chunk_new_e, i, i1, i2, i3max, i4max, ia, ia1, ia2, iao, iao1, iao2, &
      iaorb, iat, iat1, iat2, iat3, iat3m, iat4, iat4m, ib, ic, id, ierr, &
      ifail(atwp%nattot), ii, ijk, il_psi_1d, il_psi_2d, il_rhoe_1d, &
      il_rhoe_2d, index1, index2, index4, ip, is, is1, is2, isub, isw, &
      iwork(5*atwp%nattot), ixx, j, k, ki, kl, l, lan_max, langrp, &
      lanlist(parai%nproc), lscr, msglen, msweep, n3, n4, nao, nao1, nao2, &
      natst, ndd(0:parai%nproc-1,2), nfound, nolan, nomax, norb, norbx, &
      num_eig, numin, nx, old_send_cnt_e(parai%nproc)
    INTEGER :: send_cnt(parai%nproc), send_cnt_e(parai%nproc), &
      send_displ(parai%nproc), send_displ_e(parai%nproc), start, str, str1, &
      strc, sz
    LOGICAL                                  :: conv_flag, debug, ferror, flg
    REAL(real_8) :: abcm, abtol, alpha, beta, conv_tol, d(atwp%nattot), &
      dasum, dnrm2, e(atwp%nattot), eig(atwp%nattot), &
      foc(loadc_foc_array_size), ld(atwp%nattot), le(atwp%nattot), &
      local_alpha, local_beta, nrm, ql, qm, rlen, scal, sfc, sm, sm_old, &
      teig(atwp%nattot), unac, vl, vu, work(5*atwp%nattot)
    REAL(real_8), ALLOCATABLE :: abc(:), abcd(:), bab(:,:), comp(:), &
      cscr1(:,:), e_temp(:), eig_vec(:), lanv(:), local_pw(:), pw(:), &
      qa(:,:), rhoe(:,:), rmat(:), scr(:), small_xsmat(:), smat(:,:), &
      teig_vec(:), tmat(:), v0(:), v_prev(:), vcur(:), w(:), xtmat(:,:), &
      z(:,:), z_temp(:)
    REAL(real_8), EXTERNAL                   :: ddot

    CALL tiset('    DIST_PROWFN',isub)
    IF (MOD(atwp%nattot,parai%nproc).NE.0) CALL stopgm('DIST_PROWFN', &
         'NATTOT must be divisible by NPROC',& 
         __LINE__,__FILE__)
    conv_flag = .FALSE.
    nomax = MAX(prop2%numorb,atwp%nattot)
    ALLOCATE(comp(prop2%numorb),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(catom(ncpw%ngw,atwp%nattot),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(qa(ions1%nat,5),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(qa)!,5*ions1%nat)
    IF (prop1%dpan) THEN
       ALLOCATE(smat(atwp%nattot,atwp%nattot),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    n3=ions1%nat*ions1%nat*ions1%nat/6 + 1
    IF (prop1%dpan.AND.prop2%ncen.GE.3) THEN
       ALLOCATE(abc(n3),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    n4=ions1%nat*ions1%nat*ions1%nat*ions1%nat/24 + 1
    IF (prop1%dpan.AND.prop2%ncen.GE.4) THEN
       ALLOCATE(abcd(n4),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! Allocation of RHOE and PSI
    CALL rhoe_psi_size(il_psi_1d=il_psi_1d,il_psi_2d=il_psi_2d,&
         il_rhoe_1d=il_rhoe_1d,il_rhoe_2d=il_rhoe_2d)
    IF (prop1%locl.AND.(numpr.NE.0)) THEN
       ALLOCATE(psi(il_psi_1d,il_psi_2d),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(rhoe(il_rhoe_1d,il_rhoe_2d),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ELSEIF (pslo_com%tivan) THEN
       ALLOCATE(psi(il_psi_1d,il_psi_2d),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    debug=.FALSE.
    IF (atwp%nattot.LE.0) RETURN
    CALL my_set_orbdist(atwp%nattot,cnti%nstblk,norbx)
    norb=paraw%nwa12(parai%mepos,2)-paraw%nwa12(parai%mepos,1)+1
    norb=MAX(0,norb)
    sz = parai%nproc
    ALLOCATE(cscr(ncpw%ngw,atwp%nattot),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(cscr1(ncpw%ngw,2*atwp%nattot),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(xxmat(atwp%nattot,norbx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! =================================================================
    ! CREATE LANCZOS GROUP
    ! =================================================================
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
    ! =================================================================
    ! ==--------------------------------------------------------------==
    ! SCR ALLOCATION AND PARTITION (SCRATCH ARRAY).
    CALL give_scr_dist_prowfn(lscr,tag,norbx)
    ALLOCATE(scr(lscr),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    IF (paral%parent) CALL prmem('    DIST_PROWFN')
    ! ==--------------------------------------------------------------==
    ! Load AO in PW basis
    iaorb=1
    iat=0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          iat=iat+1
          CALL loadc(catom(1,iaorb),foc,ncpw%ngw,ncpw%ngw,atwp%nattot,SIZE(foc),&
               is,iat,natst)
          DO ixx=iaorb,iaorb+natst-1
             sfc=dotp(ncpw%ngw,catom(:,ixx),catom(:,ixx))
             CALL mp_sum(sfc,parai%allgrp)
             sfc=1._real_8/SQRT(sfc)
             CALL zdscal(ncpw%ngw,sfc,catom(1,ixx),1)
          ENDDO
          iaorb=iaorb+natst
       ENDDO
    ENDDO
    ! ..overlap matrix
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            CALL fileopen(41,'OVERLAP',fo_def+fo_ufo,ferror)
    ENDIF
    ALLOCATE(rmat(atwp%nattot * norb),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(tmat(atwp%nattot * norbx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    CALL zeroing(rmat)!,atwp%nattot*norb)
    ALLOCATE(xsmat(atwp%nattot,norbx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(xtmat(atwp%nattot,norbx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(xsmat)!,atwp%nattot*norbx)
    DO ip=0,parai%nproc-1
       nx=paraw%nwa12(ip,2)-paraw%nwa12(ip,1)+1
       CALL zeroing(tmat)!,atwp%nattot*norbx)
       IF (nx.GT.0) THEN
          CALL ovlap2(ncpw%ngw,atwp%nattot,nx,tmat,catom,catom(1,paraw%nwa12(ip,1)),&
               .TRUE.)
          CALL mp_sum(tmat,rmat,atwp%nattot*norbx,parai%allgrp)
          IF (parai%me.EQ.ip) THEN
             CALL dcopy(atwp%nattot*norbx,rmat(1),1,xsmat(1,1),1)
          ENDIF
          IF (paral%parent) THEN
             ! WRITE THE CURRENT COLUMNS CHUNK TO FILE
             DO ki=1, nx
                DO kl=1, atwp%nattot
                   IF (paral%io_parent)&
                        WRITE(41)  rmat(kl+(ki-1)*atwp%nattot)
                ENDDO
             ENDDO
          ENDIF
       ENDIF
    ENDDO
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            CALL fileclose(41)
    ENDIF
    CALL dcopy(atwp%nattot*norbx,xsmat(1,1),1,xtmat(1,1),1)
    ! ====================================================
    ! LANCZOS ITERATION TO DIAGONALIZE THE OVERLAP MATRIX
    ! ====================================================
    ! ==================================================================
    ! ONLY PROCESSORS THAT HAVE A NON-ZERO CHUNK OF THE RESTRICTED
    ! HAMILTONIAN PARTICIPATE IN LANCZOS
    ! ==================================================================
    IF ((paraw%nwa12(parai%me,2)-paraw%nwa12(parai%me,1)+1).NE.0) THEN
       ! ==================================================================
       ! START THE LANCZOS ITERATION WITH FULL GRAM-SCHIDT REORTH.
       ! ==================================================================
       ! --------------------------------------
       ! MEMORY ALLOCATIONS AND INITIALIZATIONS
       ! --------------------------------------
       ! EACH PROCESSOR ALLOCATES MEMORY FOR ROW-WISE DISTRIBUTED LANCZOS
       ! BASIS
       sz = 0
       chunk_new = paraw%nwa12(parai%me,2)-paraw%nwa12(parai%me,1)+1
       DO ip=0, parai%nproc-1
          IF ((paraw%nwa12(ip,2)-paraw%nwa12(ip,1)+1).NE.0) sz = sz+1
       ENDDO
       nolan = sz
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
       CALL zeroing(w)!,chunk_new)
       ! DETERMINE SEND AND RECEIVE BUFFERS
       msglen=8/2
       CALL my_concat(chunk_new, send_cnt(1), msglen, langrp)
       send_displ(1) = 0
       DO index4=2, sz
          send_displ(index4) = send_displ(index4-1) + send_cnt(index4-1)
       ENDDO
       ! DETERMINE A UNIT NORM INITIAL VECTOR V0
       IF (paral%parent) THEN
          CALL repprngu_vec(atwp%nattot,v0) 
          nrm = dnrm2(atwp%nattot, v0(1), 1)
          nrm = 1/nrm
          CALL dscal(atwp%nattot, nrm, v0(1), 1)
       ENDIF
       ! BROADCAST INITIAL VECTOR TO ALL PROCESSORS
       CALL mp_bcast(v0, atwp%nattot, parai%source, langrp)
       CALL dcopy(atwp%nattot, v0(1), 1, vcur(1), 1)
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
       ! ==================================================================
       ! LANCZOS LOOP
       ! ==================================================================
       num_eig = atwp%nattot
       DO index4 = 1, lan_max-1
          ! DISTRIBUTED MATRIX-VECTOR PRODUCT
          ! W: Holds the result
          CALL dgemv("T", atwp%nattot,chunk_new,1._real_8,xtmat(1,1),atwp%nattot,&
               vcur(1),1, 0._real_8, w(1), 1)

          ! DAXPY: W = W - BETA*V_(J-1)
          CALL daxpy(chunk_new, -beta, v_prev(1), 1, w(1), 1)
          ! CALCULATE ALPHA: INNER PRODUCT WITH CURRENT LANCZOS VECTOR
          local_alpha = ddot(chunk_new,w(1),1,&
               lanv((index4-1)*chunk_new+1), 1)
          CALL mp_sum(local_alpha, alpha, langrp)
          ! DAXPY: W = W - ALPHA*V_J
          CALL daxpy(chunk_new,-alpha,lanv((index4-1)*chunk_new+1),1,&
               w(1),1)
          ! ==================================================================
          ! REORTHOGONALIZE AGAINST ALL PREVIOUS BASIS VECTORS
          ! ==================================================================
          DO index1=1,index4
             local_pw(index1)=ddot(chunk_new,lanv((index1-1)*chunk_new+1),&
                  1, w(1), 1)
          ENDDO
          CALL mp_sum(local_pw, pw, index4, langrp)
          DO index1=1, index4
             CALL daxpy(chunk_new, -pw(index1),&
                  lanv((index1-1)*chunk_new+1), 1, w(1), 1)
          ENDDO
          ! ==================================================================
          ! CALCULATE BETA: NORM OF VECTOR W
          local_beta = ddot(chunk_new, w(1), 1, w(1), 1)
          CALL mp_sum(local_beta, beta, langrp)
          beta = SQRT(beta)
          ! NORMALIZE NEXT LANCZOS VECTOR
          CALL dscal(chunk_new, 1/beta, w(1), 1)
          CALL dcopy(chunk_new, w(1), 1, lanv(index4*chunk_new+1), 1)
          ! UPDATE TRIDIAGONAL MATRIX T (DIAGONAL AND SUP_DIAGONAL ENTRIES)
          d(index4) = alpha
          e(index4) = beta
          ! UPDATE AUXILIARY VECTORS
          CALL my_concatv(lanv((index4)*chunk_new+1),vcur(1),chunk_new,&
               send_cnt(1),send_displ(1),langrp)
          CALL dcopy(chunk_new,lanv((index4-1)*chunk_new+1),1,v_prev(1),1)
          ! CHECK CONVERGENCE FOR EIGENVALUES
          CALL dcopy(index4, d(1), 1, ld(1), 1)
          CALL dcopy(index4-1, e(1), 1, le(1), 1)
          abtol = -1
          ! ================================================================
          ! THIS TEST IS PRACTICALLY NEVER DONE: FOR FUTURE USE
          IF (index4.GT.num_eig*2) THEN
             !vw commented the following line as z isnt allocated at that point !
             !vw and added a stop
             !vw CALL  dstevx('N', 'I', index4, ld(1), le(1), vl,&
             !vw     vu, 1, num_eig, abtol,&
             !vw     nfound, eig(1), z(1,1), 1, work(1), iwork(1),&
             !vw     ifail(1), ierr)
             CALL stopgm(procedureN,'Z isnt allocated at that point',&
                  __LINE__,__FILE__)
             !vw
             sm = dasum(index4, eig(1), 1)
             IF ( (ABS(sm-sm_old)/sm_old).LE.conv_tol ) THEN
                conv_flag=.TRUE.
                GOTO 666
             ENDIF
             sm_old = sm
          ENDIF
          ! ================================================================
       ENDDO
666    CONTINUE
       ! ===================
       ! END OF LANCZOS LOOP
       ! ===================
       IF (.NOT.conv_flag) index4 = index4-1
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
       ! ALLOCATE MEMORY FOR EIGENVECTORS
       ALLOCATE(z(atwp%nattot,(chunk_end_e-chunk_begin_e+1)*index4/atwp%nattot),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL  dstevx('V','I',index4,d(1),e(1),vl,vu,chunk_begin_e,&
            chunk_end_e,&
            abtol, nfound, teig(1), z(1,1), index4, work(1),&
            iwork(1), ifail(1), ierr)
       DO index1=1, chunk_new_e
          teig(index1) = 1/SQRT(teig(index1))
       ENDDO
       ! ==================================================================
       ! ALL PROCESSORS GET ALL (SQUARED-RECIPROCAL) EIGENVALUES
       ! ==================================================================
       CALL my_concatv(teig(1),eig(1),chunk_new,&
            send_cnt(1),send_displ(1),langrp)
       ! ==================================================================
       ! CALCULATION OF APPROXIMATE EIGENVECTORS
       ! EACH PROCESSOR WILL STORE COMPLETE EIGENVECTORS
       ! EIGENVECTORS ARE DISTRIBUTED ACROSS PROCESSORS
       ! ==================================================================
       ALLOCATE(z_temp((chunk_end_e-chunk_begin_e+1)*index4),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ! ALLOCATION OF MEMORY FOR APPROXIMATE EIGENVECTORS
       ALLOCATE(eig_vec(atwp%nattot*norbx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(teig_vec(atwp%nattot*norbx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ! ALLOCATE ADDITIONAL MEMORY FOR BOOK KEEPING
       msglen=8/2
       CALL my_concat(chunk_new, old_send_cnt_e(1), msglen, langrp)
       ALLOCATE(e_temp(2*old_send_cnt_e(parai%me+1)*MAX(chunk_new_e,send_cnt(1))),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       strc = 1
       DO index1=1, nolan
          IF (parai%me.EQ.(index1-1) ) THEN
             CALL dcopy(atwp%nattot*chunk_new_e, z(1,1), 1, z_temp(1), 1)
             DO index2 = 1, sz
                send_cnt_e(index2) = old_send_cnt_e(index2) *&
                     send_cnt(parai%me+1)
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
          CALL mp_bcast(z_temp, send_cnt(index1)*atwp%nattot, index1-1,&
               langrp)
          ! MULTIPLY RECEIVED EIGENVECTORS WITH MY PART OF THE LANCZOS
          ! BASIS
          CALL dgemm("N","N",chunk_new,send_cnt(index1),atwp%nattot,1._real_8,&
               lanv(1),chunk_new,z_temp(1),atwp%nattot,&
               0._real_8, eig_vec(strc), chunk_new)
          ! NORMALIZE EIGENVECTORS TO UNITY
          DO index4 = 1, send_cnt(index1)
             alpha =ddot(chunk_new, eig_vec(strc+(index4-1)*chunk_new),&
                  1, eig_vec(strc+(index4-1)*chunk_new), 1)
             CALL mp_sum(alpha, scal, langrp)
             scal = 1/SQRT(scal)
             CALL dscal(chunk_new, scal,&
                  eig_vec(strc+(index4-1)*chunk_new),1)
          ENDDO
          strc = strc + chunk_new*send_cnt(index1)
       ENDDO
       ! ==================================================================
       ALLOCATE(small_xsmat(norbx*norbx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       str = 1
       DO index1 = 1, nolan
          ! ==========
          ! LOCAL COPY
          ! ==========
          IF (parai%me.EQ.(index1-1) ) THEN
             CALL dcopy(atwp%nattot*chunk_new, eig_vec(1), 1, z_temp(1), 1)
             ! ===========================
             ! SCALE WITH RECIPROCAL SQRTs
             ! ===========================
             DO index4=1, atwp%nattot
                CALL dscal(chunk_new, eig(index4),&
                     z_temp((index4-1)*chunk_new+1), 1)
             ENDDO
             DO index2 = 1, nolan
                send_cnt_e(index2) = old_send_cnt_e(index2)*chunk_new
             ENDDO
             send_displ_e(1) = 0
             DO index2=2, sz
                send_displ_e(index2) = send_displ_e(index2-1) +&
                     send_cnt_e(index2-1)
             ENDDO
          ENDIF
          ! =================
          ! ARRAYS FOR GATHER
          ! =================
          CALL mp_bcast(send_cnt_e, nolan, index1-1, langrp)
          CALL mp_bcast(send_displ_e, nolan, index1-1, langrp)
          ! =====================
          ! BROADCAST  LOCAL COPY
          ! =====================
          CALL mp_bcast(z_temp, atwp%nattot*send_cnt(index1),&
               index1-1, langrp)
          ! ===============================================
          ! MULTIPLY WITH LOCAL PART OF EIGENVECTORS MATRIX
          ! ===============================================
          CALL dgemm("N","T", send_cnt(index1), chunk_new,&
               atwp%nattot, 1._real_8, z_temp(1),&
               send_cnt(index1), eig_vec(1), chunk_new,&
               0._real_8, small_xsmat(1),&
               send_cnt(index1) )
          ! ===============================
          ! GATHER BACK LOCAL PART OF XSMAT 
          ! ===============================
          CALL my_source_concatv(small_xsmat(1),teig_vec(1),&
               chunk_new*send_cnt(index1),&
               send_cnt_e(1), send_displ_e(1),&
               index1-1, langrp)
       ENDDO
       DO index1=1, atwp%nattot
          CALL dcopy(chunk_new, teig_vec((index1-1)*chunk_new+1),&
               1, xsmat(index1,1), atwp%nattot)
       ENDDO
    ENDIF
    ! =================================================================  
    ! All processors work here since all of them
    ! hold a pice of the atomic orbitals
    ! =================================================================
    CALL mp_sync(parai%allgrp)
    ! Make sure that every processor knows SEND_CNT 
    CALL mp_bcast(send_cnt, parai%nproc, parai%source, parai%allgrp)
    ! Allocate memory Z_TEMP if not allocated so far
    IF ((parai%me+1).GT.nolan) THEN
       ALLOCATE(z_temp(atwp%nattot*(chunk_end_e-chunk_begin_e+1)*index4/atwp%nattot),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! Orthonormalize atomic orbitals
    CALL zeroing(cscr)!, atwp%nattot*ngw)
    strc = 1
    DO index1=1, nolan
       IF ((parai%me+1).EQ.index1) THEN
          CALL dcopy(atwp%nattot*chunk_new, xsmat(1,1), 1, z_temp(1), 1)
       ENDIF
       CALL mp_bcast(z_temp, atwp%nattot*send_cnt(index1),&
            index1-1, parai%allgrp)
       CALL dgemm('N','N',2*ncpw%ngw,send_cnt(index1),atwp%nattot,1.0_real_8,&
            catom(1,1),2*ncpw%ngw,z_temp(1),atwp%nattot,&
            0.0_real_8,cscr(1,strc),2*ncpw%ngw)
       strc = strc + send_cnt(index1)
    ENDDO
    CALL dcopy(2*atwp%nattot*ncpw%ngw,cscr(1,1),1,catom(1,1),1)
    ! =================================================================
    ! Overlap AO with wavefunctions
    ! =================================================================
    DO ip=0,parai%nproc-1
       nx=paraw%nwa12(ip,2)-paraw%nwa12(ip,1)+1
       CALL zeroing(z_temp(1:atwp%nattot*norbx))!,atwp%nattot*norbx)
       IF (nx.GT.0) THEN
          CALL ovlap2(ncpw%ngw,nx,prop2%numorb,z_temp,catom(1,paraw%nwa12(ip,1)),c0,&
               .TRUE.)
          CALL mp_sum(z_temp,rmat,prop2%numorb*norbx,parai%allgrp)
          IF (parai%me.EQ.ip) THEN
             DO index1=1, nx
                CALL dcopy(prop2%numorb,rmat(index1),nx,&
                     xxmat(1,index1),1)
             ENDDO
          ENDIF
       ENDIF
    ENDDO
    CALL mp_sync(parai%allgrp)

    DEALLOCATE(rmat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__) ! TODO check
    ! =================================================================
    ! =================================================================
    ! VANDERBILD CASE NOT IMPLEMENTED YET
    ! =================================================================
    IF (pslo_com%tivan) CALL stopgm('DIST_PROWFN','VANDERBILD NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    IF (pslo_com%tivan) THEN
       ! Calculate augmentation charges
       CALL augchg(fnl,crge%f,qa(1,4),prop2%numorb)
       ALLOCATE(fnl2(1,ions1%nat,prop2%numorb,maxsys%nhxs,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(fnl,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(fnl(1,ions1%nat,atwp%nattot,maxsys%nhxs,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL rnlsm(catom,atwp%nattot,1,1,.FALSE.)
       CALL zeroing(fnl2)!,ions1%nat*maxsys%nhxs*prop2%numorb)
       DO i=1,prop2%numorb
          DO j=1,atwp%nattot
             CALL daxpy(ions1%nat*maxsys%nhxs,xxmat(j,i),&
                  fnl(1,1,1,j,1),1,fnl2(1,1,1,i,1),1)
          ENDDO
       ENDDO
       CALL augchg(fnl2,crge%f,qa(1,5),prop2%numorb)
    ENDIF
    ! =================================================================
    DO i=1,prop2%numorb
       z_temp(1)=dotp(ncpw%ngw,c0(:,i),c0(:,i))
       CALL mp_sum(z_temp(1),rlen, parai%allgrp)
       ! ===============================================================
       ! COMP IS CALCULATED IN PARALLEL SINCE XXMAT IS DISTRIBUTED
       ! ROW-WISE TO ALL PROCESSORS THAT TAKE PART IN LANGRP
       ! ===============================================================
       IF (chunk_new.NE.0) THEN
          tmat(1)=ddot(chunk_new, xxmat(i,1), atwp%nattot, xxmat(i,1), atwp%nattot)
          CALL mp_sum(tmat(1), comp(i), langrp)
          comp(i)=comp(i)/rlen
       ENDIF
    ENDDO

    DEALLOCATE(tmat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__) ! TODO check
    ! ================================================================= 
    ! THE NEXT CHUNK OF THE CODE IS SUPPOSED TO BE EXECUTED ONLY BY
    ! THE PARENT PROCESS. BUT XXMAT IS DISTRIBUTED ROW-WISE NOW. SO ALL
    ! PROCESSORS INVOLVED IN THE LANGRP GROUP TAKE PART
    ! =================================================================
    IF (chunk_new.NE.0) THEN
       CALL zeroing(z_temp(1:norbx*atwp%nattot))!, norbx*atwp%nattot)
       DO index1=1, chunk_new
          DO index4=1, nstate
             z_temp(index1) = z_temp(index1) +&
                  crge%f(index4,1)*xxmat(index4,index1)**2
          ENDDO
       ENDDO
       ! PARENT GETS ARRAY Z_TEMP 
       CALL my_concatv(z_temp(1),z_temp(atwp%nattot+1),chunk_new,&
            send_cnt(1),send_displ(1), langrp)
       IF (paral%parent) THEN
          ! Lowdin Population Analysis
          iat=0
          iao=0
          DO is=1,ions1%nsp
             nao=atwf_mod%numaor(is)
             DO ia=1,ions0%na(is)
                iat=iat+1
                ! QA(.,4) is the true augmentation charge!
                qa(iat,2)=ions0%zv(is)-qa(iat,4)
                DO i=iao+1,iao+nao
                   qa(iat,2)=qa(iat,2)-z_temp(atwp%nattot+i)
                ENDDO
                iao=iao+nao
             ENDDO
          ENDDO
       ENDIF
    ENDIF
    ! =================================================================
    ! Rotate back to AO basis:
    ! RESULT IS OF SIZE (NATTOT x NUMORB). SO WE DISTRIBUTE IT 
    ! ROW-WISE
    ! =================================================================
    IF (chunk_new.NE.0) THEN
       strc=1
       DO index1=1, nstate
          CALL dgemv('N',atwp%nattot,chunk_new,1._real_8,xsmat(1,1),atwp%nattot,&
               xxmat(index1,1),atwp%nattot,0._real_8,scr(1),1)
          CALL mp_sum(scr,eig_vec,atwp%nattot,langrp)
          CALL dcopy(chunk_new, eig_vec(paraw%nwa12(parai%me,1)),&
               1, xxmat(index1,1), atwp%nattot)
       ENDDO
    ENDIF
    CALL mp_sync(parai%allgrp)
    ! ===============================================================
    ! MOVE OUTPUT TO FILE HERE: IMPORTANT FOR I/O
    ! ===============================================================
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            CALL fileopen(41,'WFNCOEF',fo_def+fo_ufo,ferror)
       IF (paral%io_parent)&
            WRITE(41) atwp%nattot,ions1%nsp,(ions0%zv(is),ions0%na(is),atwf_mod%numaor(is),is=1,ions1%nsp)
    ENDIF
    IF (chunk_new.NE.0) THEN
       DO index1=1, nolan
          IF ((parai%me+1).EQ.index1) THEN
             DO index4=1, prop2%numorb
                CALL dcopy(chunk_new,xxmat(index4,1),atwp%nattot,&
                     z_temp((index4-1)*chunk_new+1),1)
             ENDDO
          ENDIF
          CALL mp_bcast(z_temp, prop2%numorb*send_cnt(index1),&
               index1-1, langrp)
          IF (paral%parent) THEN
             DO kl=1, send_cnt(index1)
                DO ki=1, prop2%numorb
                   IF (paral%io_parent)&
                        WRITE(41) z_temp(kl+(ki-1)*send_cnt(index1))
                ENDDO
             ENDDO
          ENDIF
       ENDDO
    ENDIF
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            CALL fileclose(41)
    ENDIF
    CALL mp_sync(parai%allgrp)
    ! ===============================================================
    IF (paral%parent) THEN
       ! Print the wavefunctions
       CALL mklabel(label)
       IF (prop1%prto) THEN
          IF (paral%io_parent)&
               WRITE(6,'(/,A,/)') ' WAVEFUNCTIONS IN ATOMIC ORBITAL BASIS'
          IF (cntl%tlsd) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(21X,A)') ' ****** ALPHA SPIN ****** '
             CALL dist_prtmat(xxmat(1,1),atwp%nattot,spin_mod%nsup,label,comp,crge%f)
             IF (paral%io_parent)&
                  WRITE(6,'(/,21X,A)') ' ****** BETA  SPIN ****** '
             CALL dist_prtmat(xxmat(1,spin_mod%nsup+1),atwp%nattot,spin_mod%nsdown,label,&
                  comp(spin_mod%nsup+1),crge%f(spin_mod%nsup+1,1))
          ELSE
             CALL dist_prtmat(xxmat,atwp%nattot,prop2%numorb,label,comp,crge%f)
          ENDIF
       ENDIF
       IF (pslo_com%tivan) THEN
          ! ===============================================================
          ! THE TVAN CASE IS NOT IMPLEMENTED YET. THE CODE EXITS ABOVE
          ! ===============================================================       
          ! Print Augmentation Charges
          IF (paral%io_parent)&
               WRITE(6,'(/,A,/)')&
               ' VANDERBILT AUGMENTATION CHARGES '
          IF (paral%io_parent)&
               WRITE(6,'(A,A)') '       ATOM        FULL BASIS    ATOM BASIS'
          iat=0
          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                iat=iat+1
                IF (paral%io_parent)&
                     WRITE(6,'(I4,4X,A,9X,F10.3,4X,F10.3)')&
                     iat,elem%el(ions0%iatyp(is)),qa(iat,4),qa(iat,5)
             ENDDO
          ENDDO
       ENDIF
    ENDIF
    ! =================================================================
    IF (prop1%mpan.AND.(chunk_new.GT.0)) THEN
       ! Mulliken Population Analysis
       str = 1
       DO index1 = 1, nolan
          ! ==========
          ! LOCAL COPY
          ! ==========
          IF (parai%me.EQ.(index1-1) ) THEN
             DO index4=1, nstate
                CALL dcopy(chunk_new,xxmat(index4,1),atwp%nattot,&
                     z_temp((index4-1)*chunk_new+1),1)
             ENDDO
             ! =====================
             ! SCALE WITH OCCUPATION
             ! =====================
             DO index4=1, nstate
                CALL dscal(chunk_new, crge%f(index4,1),&
                     z_temp((index4-1)*chunk_new+1), 1)
             ENDDO
             DO index2 = 1, nolan
                send_cnt_e(index2) = old_send_cnt_e(index2)*chunk_new
             ENDDO
             send_displ_e(1) = 0
             DO index2=2, sz
                send_displ_e(index2) = send_displ_e(index2-1) +&
                     send_cnt_e(index2-1)
             ENDDO
          ENDIF
          ! =================
          ! ARRAYS FOR GATHER
          ! =================
          CALL mp_bcast(send_cnt_e, nolan, index1-1, langrp)
          CALL mp_bcast(send_displ_e, nolan, index1-1, langrp)

          ! =====================
          ! BROADCAST  LOCAL COPY
          ! =====================
          CALL mp_bcast(z_temp, nstate*send_cnt(index1),&
               index1-1, langrp)
          ! ===============================================
          ! MULTIPLY WITH LOCAL PART OF EIGENVECTORS MATRIX
          ! ===============================================
          CALL dgemm("N","N", send_cnt(index1), chunk_new,&
               nstate, 1._real_8, z_temp(1),&
               send_cnt(index1), xxmat(1,1), atwp%nattot,&
               0._real_8, small_xsmat(1),&
               send_cnt(index1) )

          ! ===============================
          ! GATHER BACK LOCAL PART OF XSMAT
          ! ===============================
          CALL my_source_concatv(small_xsmat(1),teig_vec(1),&
               chunk_new*send_cnt(index1),&
               send_cnt_e(1), send_displ_e(1),&
               index1-1, langrp)
       ENDDO
       DO index1=1, atwp%nattot
          CALL dcopy(chunk_new, teig_vec((index1-1)*chunk_new+1),&
               1, xsmat(index1,1), atwp%nattot)
       ENDDO
       CALL zeroing(scr(1:atwp%nattot))!,atwp%nattot)
       CALL zeroing(z_temp(1:atwp%nattot))!,atwp%nattot)
       DO j=1,chunk_new
          DO i=1,atwp%nattot
             z_temp(i)=z_temp(i)+xsmat(i,j)*xtmat(i,j)
          ENDDO
       ENDDO
       CALL mp_sum(z_temp, scr, atwp%nattot, langrp)
       ! ==============================================================        
       IF (paral%parent) THEN
          iat=0
          iao=0
          DO is=1,ions1%nsp
             nao=atwf_mod%numaor(is)
             DO ia=1,ions0%na(is)
                iat=iat+1
                ! QA(.,4) is the true augmentation charge!
                qa(iat,1)=ions0%zv(is)-qa(iat,4)
                DO i=iao+1,iao+nao
                   qa(iat,1)=qa(iat,1)-scr(i)
                ENDDO
                iao=iao+nao
             ENDDO
          ENDDO
       ENDIF
       ! ==============================================================
       ! Bond Orders
       CALL dcopy(atwp%nattot*norbx,xxmat(1,1),1,scr(1),1)
       ! CALL DGEMM('N','N',NATTOT,NATTOT,NATTOT,1.0_real_8,XSMAT(1,1),
       ! &               NATTOT,XTMAT(1,1),NATTOT,0.0_real_8,XXMAT(1,1),NATTOT)

       ! =================================================================
       str = 1
       DO index1 = 1, nolan
          ! ==========
          ! LOCAL COPY
          ! ==========
          IF (parai%me.EQ.(index1-1) ) THEN
             DO index4=1, atwp%nattot
                CALL dcopy(chunk_new,xsmat(index4,1),atwp%nattot,&
                     z_temp((index4-1)*chunk_new+1),1)
             ENDDO
             DO index2 = 1, nolan
                send_cnt_e(index2) = old_send_cnt_e(index2)*chunk_new
             ENDDO
             send_displ_e(1) = 0
             DO index2=2, sz
                send_displ_e(index2) = send_displ_e(index2-1) +&
                     send_cnt_e(index2-1)
             ENDDO
          ENDIF
          ! =================
          ! ARRAYS FOR GATHER
          ! =================
          CALL mp_bcast(send_cnt_e, nolan, index1-1, langrp)
          CALL mp_bcast(send_displ_e, nolan, index1-1, langrp)
          ! =====================
          ! BROADCAST  LOCAL COPY
          ! =====================
          CALL mp_bcast(z_temp, atwp%nattot*send_cnt(index1),&
               index1-1, langrp)
          ! ===============================================
          ! MULTIPLY WITH LOCAL PART OF EIGENVECTORS MATRIX
          ! ===============================================
          CALL dgemm("N","N", send_cnt(index1), chunk_new,&
               atwp%nattot, 1._real_8, z_temp(1),&
               send_cnt(index1), xtmat(1,1), atwp%nattot,&
               0._real_8, small_xsmat(1),&
               send_cnt(index1) )
          ! ===============================
          ! GATHER BACK LOCAL PART OF XSMAT
          ! ===============================
          CALL my_source_concatv(small_xsmat(1),teig_vec(1),&
               chunk_new*send_cnt(index1),&
               send_cnt_e(1), send_displ_e(1),&
               index1-1, langrp)
       ENDDO
       ! ===============================================================
       ! XXMAT HOLDS A CHUNK OF THE ROWS OF THE RESULT:
       ! 1ST COLUMN XXMAT(:,1) HOLDS THE FIRST ROW ETC
       ! =============================================================== 
       DO index1=1, atwp%nattot
          CALL dcopy(chunk_new, teig_vec((index1-1)*chunk_new+1),&
               1, xxmat(index1, 1), atwp%nattot)
       ENDDO
       ! ===============================================================
       ! ==============================================================
       ! OBSERVE THAT IN REALITY ONLY THE SUMS OF THE  ENTRIES ALONG 
       ! THE COLUMNS OF ARRAY BAB ARE NEEDED.
       ! WE ONLY STORE THESE.
       ! ==============================================================
       CALL zeroing(z_temp(1:atwp%nattot*norbx))!,atwp%nattot*norbx)
       ! ==============================================================
       ! FIRST WE NEED TO CALCULATE THE HADAMARD (SCALAR-ELEMENTWISE)
       ! PRODUCT OF XXMAT WITH ITS TRANSPOSE
       ! ============================================================== 
       strc = 0
       DO index1 = 1, nolan
          ! ==========
          ! LOCAL COPY
          ! ==========
          IF (parai%me.EQ.(index1-1) ) THEN
             DO index4=1, chunk_new
                CALL dcopy(atwp%nattot,xxmat(1,index4),1,&
                     z_temp((index4-1)*atwp%nattot+1),1)
             ENDDO

             DO index2 = 1, nolan
                send_cnt_e(index2) = old_send_cnt_e(index2)*chunk_new
             ENDDO
             send_displ_e(1) = 0
             DO index2=2, sz
                send_displ_e(index2) = send_displ_e(index2-1) +&
                     send_cnt_e(index2-1)
             ENDDO
          ENDIF
          ! =================
          ! ARRAYS FOR GATHER
          ! =================
          CALL mp_bcast(send_cnt_e, nolan, index1-1, langrp)
          CALL mp_bcast(send_displ_e, nolan, index1-1, langrp)
          ! =====================
          ! BROADCAST  LOCAL COPY
          ! =====================
          CALL mp_bcast(z_temp, atwp%nattot*send_cnt(index1),&
               index1-1, langrp)
          ! ===========================
          ! ELEMENT-WISE MULTIPLICATION
          ! ===========================
          str  = 1
          str1 = 1
          DO index2=(index1-1)*chunk_new+1, index1*chunk_new
             DO index4=1, send_cnt(index1)
                teig_vec(str+strc) = xxmat(index2,index4) *&
                     z_temp((str1-1)*atwp%nattot+index4+&
                     paraw%nwa12(parai%me,1)-1 )
                str = str+1
             ENDDO
             str1 = str1 + 1
          ENDDO

          strc = strc + chunk_new * send_cnt(index1)
       ENDDO
       CALL zeroing(z_temp(1:atwp%nattot*norbx))!,atwp%nattot*norbx)
       CALL zeroing(eig_vec(1:atwp%nattot))!,atwp%nattot)
       ! DETERMINE DISTRIBUTION OF BAB
       chunk_new_e = ions1%nat/nolan
       DO index1=1, nolan-1
          ndd(index1-1,1) = chunk_new_e*(index1-1)+1
          ndd(index1-1,2) = chunk_new_e*index1
       ENDDO
       ndd(nolan-1,1) = chunk_new_e*(nolan-1)+1
       ndd(nolan-1,2) = ions1%nat
       chunk_begin_e = parai%me*chunk_new_e+1
       chunk_end_e   = (parai%me+1)*chunk_new_e
       IF (parai%me.EQ.(nolan-1)) chunk_end_e = ions1%nat
       IF (parai%me.EQ.(nolan-1)) chunk_new_e = ions1%nat-(nolan-1)*chunk_new_e
       index2=0
       DO index1=1, nolan
          index2 = MAX(index2,ndd(index1-1,2)-ndd(index1-1,1)+1)
       ENDDO
       norbx = index2
       ! ALLOCATE MEMORY
       ALLOCATE(bab(ions1%nat,index2),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(bab)!,ions1%nat*norbx)
       IF (paral%parent) THEN
          IF (paral%io_parent)&
               CALL fileopen(41,'BAB_tmp',fo_def+fo_ufo,ferror)
       ENDIF
       CALL zeroing(eig_vec(1:ions1%nat))!,ions1%nat)
       flg = .FALSE.
       str = 1
       iat1=0
       iao1=0
       start = 1
       DO is1=1,ions1%nsp
          nao1=atwf_mod%numaor(is1)
          DO ia1=1,ions0%na(is1)
             iat1=iat1+1
             DO i1=iao1+1,iao1+nao1
                iat2=0
                iao2=0
                DO is2=1,ions1%nsp
                   nao2=atwf_mod%numaor(is2)
                   DO ia2=1,ions0%na(is2)
                      iat2=iat2+1
                      IF (iat2.LT.iat1) THEN
                         DO i2=iao2+1,iao2+nao2
                            IF ((i1.GE.paraw%nwa12(parai%me,1)).AND.&
                                 (i1.LE.paraw%nwa12(parai%me,2))) THEN

                               eig_vec(iat2) = eig_vec(iat2) +&
                                    teig_vec(chunk_new*(i2-1)+str)

                               flg = .TRUE.
                               eig_vec(iat1) = eig_vec(iat1) +&
                                    teig_vec(chunk_new*(i2-1)+str)

                               z_temp(iat2) = z_temp(iat2) +&
                                    teig_vec(chunk_new*(i2-1)+str)

                            ENDIF
                         ENDDO
                      ENDIF
                      iao2=iao2+nao2
                   ENDDO
                ENDDO
                IF ((flg).OR.(parai%me.EQ.0)) THEN
                   str = str+1
                ENDIF
             ENDDO
             iao1=iao1+nao1
             CALL mp_sum(z_temp, z_temp(ions1%nat+1:), ions1%nat, langrp)
             IF ((iat1.GE.chunk_begin_e).AND.(iat1.LE.chunk_end_e)) THEN
                CALL dcopy(ions1%nat, z_temp(ions1%nat+1), 1, bab(1,start), 1)
                start = start + 1
             ENDIF
             CALL zeroing(z_temp(1:2*ions1%nat))!,2*ions1%nat)
             CALL mp_sync(langrp)
          ENDDO
       ENDDO
       CALL mp_sum(eig_vec, qa(:,3), ions1%nat, langrp)
       CALL symm_da(bab,z_temp(1),z_temp(norbx*norbx+1),&
            ions1%nat,ndd(0,1),ndd(0,2),&
            norbx,parai%me,nolan,langrp)
       CALL dscal(ions1%nat*norbx, 2._real_8, bab(1,1), 1)
       IF (parai%me.EQ.0) THEN
          DO index1=1, chunk_new_e
             DO index2=1, ions1%nat
                IF (paral%io_parent)&
                     WRITE(41) bab(index2,index1)
             ENDDO
          ENDDO
       ENDIF
       DO index1=2, nolan
          IF (parai%me.EQ.(index1-1)) THEN
             CALL mp_send(bab,ions1%nat*norbx,0,1,parai%allgrp)
          ELSEIF (parai%me.EQ.0) THEN
             CALL mp_recv(z_temp,ions1%nat*norbx,index1-1,1,parai%allgrp)
             DO index2=1, ions1%nat*(ndd(index1-1,2)-ndd(index1-1,1)+1)
                IF (paral%io_parent)&
                     WRITE(41) z_temp(index2)
             ENDDO
          ENDIF
          CALL mp_sync(langrp)
       ENDDO
       IF ((paral%parent).AND.paral%io_parent)&
            CALL fileclose(41)
       ! ==============================================================
       ! Print Population Analysis
       IF (paral%parent) THEN
          IF (paral%io_parent)&
               WRITE(6,'(/,A,/)')&
               ' POPULATION ANALYSIS FROM PROJECTED WAVEFUNCTIONS'
          IF (paral%io_parent)&
               WRITE(6,'(A)')&
               '      ATOM          MULLIKEN        LOWDIN       VALENCE'
          iat=0
          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                iat=iat+1
                IF (paral%io_parent)&
                     WRITE(6,'(I4,4X,A,9X,F10.3,4X,F10.3,4X,F10.3)')&
                     iat,elem%el(ions0%iatyp(is)),qa(iat,1),qa(iat,2),qa(iat,3)
             ENDDO
          ENDDO
          WRITE(6,'(a,e16.8)') ' ChkSum(POP_MUL) = ',SUM(ABS(qa(:,1:3)))
          qm=0.0_real_8
          ql=0.0_real_8
          DO i=1,ions1%nat
             qm=qm+qa(i,1)
             ql=ql+qa(i,2)
          ENDDO
          IF (paral%io_parent)&
               WRITE(6,'(A,F10.3,4X,F10.3)') '  UNASSIGNED CHARGE',qm,ql
          ! Print Bond Orders
          IF (paral%io_parent)&
               WRITE(6,'(/,A)')&
               ' MAYER BOND ORDERS FROM PROJECTED WAVEFUNCTIONS'
          iat=0
          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                iat=iat+1
                IF (paral%io_parent)&
                     WRITE(label(iat),'(I4,1X,A2)') iat,elem%el(ions0%iatyp(is))
             ENDDO
          ENDDO
          ! ===========================================================
          ! OPEN AUXILLIARY FILE BAB_tmp THAT HOLDS BAB
          ! ===========================================================
          msweep=ions1%nat/8
          IF (MOD(ions1%nat,8).NE.0) msweep=msweep+1
          DO isw=1,msweep
             IF (paral%io_parent)&
                  CALL fileopen(41,'BAB_tmp',fo_def+fo_ufo,ferror)
             i1=(isw-1)*8+1
             i2=MIN(8*isw,ions1%nat)
             IF (paral%io_parent)&
                  WRITE(6,'(/,9(A8))') '        ',(label(ii),ii=i1,i2)
             DO ia=1,ions1%nat
                DO index1=1, i1-1
                   IF (paral%io_parent)&
                        READ(41) z_temp(1)
                ENDDO
                DO ii=i1, i2
                   IF (paral%io_parent)&
                        READ(41) bab(ii,1)
                ENDDO
                IF (paral%io_parent)&
                     WRITE(6,'(A8,8(F7.3,1X))')&
                     label(ia),(bab(ii,1),ii=i1,i2)
                DO index1=i2+1, ions1%nat
                   IF (paral%io_parent)&
                        READ(41) z_temp(1)
                ENDDO
             ENDDO
             IF (paral%io_parent)&
                  CALL fileclose(41)
          ENDDO
          IF (paral%io_parent)&
               WRITE(6,*)

          ! CALL FILECLOSE(41)
       ENDIF
       ! ==============================================================
    ENDIF
    ! ==================================================================
    ! END MULLIKEN POPULATION ANALYSIS
    ! ==================================================================
    IF (prop1%dpan) THEN
       CALL stopgm('DIST_PROWFN','DAVIDSON NOT IMPLEMENTED ',& 
            __LINE__,__FILE__)
       ! Davidson Population Analysis
       numin=0
       DO is=1,ions1%nsp
          numin=numin+ions0%na(is)*prop2%maos(is)
          IF (atwf_mod%numaor(is).LT.prop2%maos(is)) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(A,I4)')&
                  ' INCONSISTENT ATOMIC BASIS FOR SPECIES ',is
             CALL stopgm('DIST_PROWFN',' ',& 
                  __LINE__,__FILE__)
          ENDIF
       ENDDO
       ! Density Matrix
       DO i=1,atwp%nattot
          DO j=1,atwp%nattot
             xsmat(i,j)=0.0_real_8
             DO k=1,nstate
                xsmat(i,j)=xsmat(i,j)+crge%f(k,1)*xxmat(i,k)*xxmat(j,k)
             ENDDO
          ENDDO
       ENDDO
       unac=0._real_8
       DO iat=1,ions1%nat
          unac=unac+qa(iat,4)
       ENDDO
       ! Calculate modified atomic orbitals
       CALL dcopy(atwp%nattot*atwp%nattot,xtmat(1,1),1,smat(1,1),1)
       CALL cmaos(xsmat,smat,numin,unac)
       unac=unac/REAL(ions1%nat,kind=real_8)
       ! Calculate shared atomic charges
       CALL satch(xsmat,xtmat,smat,numin,bab)
       ! Calculate 3-center shared atomic charges
       IF (prop2%ncen.GE.3) CALL satch3(xsmat,xtmat,smat,numin,&
            bab,abc)
       ! Calculate 4-center shared atomic charges
       IF (prop2%ncen.GE.4) CALL satch4(xsmat,xtmat,smat,numin,&
            bab,abc,abcd)
       IF ((prop2%ncen.GE.5).AND.paral%io_parent)&
            WRITE(6,*) prop2%ncen,' CENTER TERMS NOT IMPLEMENTED '
       ! Charges
       iat=0
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             iat=iat+1
             ! QA(.,4) is the true augmentation charge!
             qa(iat,1)=ions0%zv(is)-qa(iat,4)-unac
          ENDDO
       ENDDO
       iat3=0
       iat4=0
       DO ia=1,ions1%nat
          qa(ia,1)=qa(ia,1)-bab(ia,ia)
          DO ib=ia+1,ions1%nat
             qa(ia,1)=qa(ia,1)+bab(ia,ib)/2._real_8
             qa(ib,1)=qa(ib,1)+bab(ia,ib)/2._real_8
             IF (prop2%ncen.GE.3) THEN
                DO ic=ib+1,ions1%nat
                   iat3=iat3+1
                   qa(ia,1)=qa(ia,1)-abc(iat3)/3._real_8
                   qa(ib,1)=qa(ib,1)-abc(iat3)/3._real_8
                   qa(ic,1)=qa(ic,1)-abc(iat3)/3._real_8
                   IF (prop2%ncen.GE.4) THEN
                      DO id=ic+1,ions1%nat
                         iat4=iat4+1
                         qa(ia,1)=qa(ia,1)+abcd(iat4)/4._real_8
                         qa(ib,1)=qa(ib,1)+abcd(iat4)/4._real_8
                         qa(ic,1)=qa(ic,1)+abcd(iat4)/4._real_8
                         qa(id,1)=qa(id,1)+abcd(iat4)/4._real_8
                      ENDDO
                   ENDIF
                ENDDO
             ENDIF
          ENDDO
       ENDDO
       ! Print Population Analysis
       IF (paral%io_parent)&
            WRITE(6,'(/,A,/)')&
            ' POPULATION ANALYSIS FROM PROJECTED WAVEFUNCTIONS'
       IF (paral%io_parent)&
            WRITE(6,'(A)')&
            '       ATOM          DAVIDSON        LOWDIN       '
       iat=0
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             iat=iat+1
             IF (paral%io_parent)&
                  WRITE(6,'(I4,4X,A,9X,F10.3,4X,F10.3)')&
                  iat,elem%el(ions0%iatyp(is)),qa(iat,1),qa(iat,2)
          ENDDO
       ENDDO
       qm=0.0_real_8
       ql=0.0_real_8
       DO i=1,ions1%nat
          qm=qm+qa(i,1)
          ql=ql+qa(i,2)
       ENDDO
       IF (paral%io_parent)&
            WRITE(6,'(A,F10.3,4X,F10.3)') '  UNASSIGNED CHARGE',qm,ql
       WRITE(6,'(a,e16.8)') ' ChkSum(POP_DAV) = ',SUM(ABS(qa(:,1:2)))
       ! Print Shared Electron Numbers
       IF (paral%io_parent)&
            WRITE(6,'(/,A)')&
            ' SHARED ELECTRON NUMBERS FROM PROJECTED WAVEFUNCTIONS'
       iat=0
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             iat=iat+1
             IF (paral%io_parent)&
                  WRITE(label(iat),'(I4,1X,A2)') iat,elem%el(ions0%iatyp(is))
          ENDDO
       ENDDO
       msweep=ions1%nat/8
       IF (MOD(ions1%nat,8).NE.0) msweep=msweep+1
       DO isw=1,msweep
          i1=(isw-1)*8+1
          i2=MIN(8*isw,ions1%nat)
          IF (paral%io_parent)&
               WRITE(6,'(/,9(A8))') '        ',(label(ii),ii=i1,i2)
          DO ia=1,ions1%nat
             IF (paral%io_parent)&
                  WRITE(6,'(A8,8(F7.3,1X))') label(ia),(bab(ia,ii),ii=i1,i2)
          ENDDO
       ENDDO
       IF (paral%io_parent)&
            WRITE(6,*)
       ! Print 3-Center Shared Electron Numbers
       IF (prop2%ncen.GE.3) THEN
          IF (paral%io_parent)&
               WRITE(6,'(/,A)')&
               ' 3-CENTER SHARED ELECTRON NUMBERS FROM PROJECTED WAVEFUNCTIONS'
          IF (paral%io_parent)&
               WRITE(6,'(A,F10.4)') ' CUTOFF FOR PRINTING IS :',prop3%cut3o
          ip=0
          abcm=0.0_real_8
          iat3=0
          DO i=1,ions1%nat
             DO j=i+1,ions1%nat
                DO k=j+1,ions1%nat
                   iat3=iat3+1
                   IF (ABS(abc(iat3)).GT.abcm) THEN
                      i3max=i+(j-1)*ions1%nat+(k-1)*ions1%nat*ions1%nat
                      iat3m=iat3
                      abcm=ABS(abc(iat3))
                   ENDIF
                   IF (abc(iat3).GT.prop3%cut3o) THEN
                      ip=ip+1
                      IF (paral%io_parent)&
                           WRITE(6,'(3A8,F10.4)')&
                           label(i),label(j),label(k),abc(iat3)
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
          IF (ip.EQ.0.AND.abcm.GT.0._real_8) THEN
             ijk=i3max
             i = MOD(ijk-1,ions1%nat)+1
             k = ijk/(ions1%nat*ions1%nat)+1
             j = (ijk-(k-1)*ions1%nat*ions1%nat)/ions1%nat+1
             IF (paral%io_parent)&
                  WRITE(6,'(3A8,F10.4)')&
                  label(i),label(j),label(k),abc(iat3m)
          ENDIF
       ENDIF
       ! Print 4-Center Shared Electron Numbers
       IF (prop2%ncen.GE.4) THEN
          IF (paral%io_parent)&
               WRITE(6,'(/,A)')&
               ' 4-CENTER SHARED ELECTRON NUMBERS FROM PROJECTED WAVEFUNCTIONS'
          IF (paral%io_parent)&
               WRITE(6,'(A,F10.4)') ' CUTOFF FOR PRINTING IS :',prop3%cut4o
          ip=0
          abcm=0.0_real_8
          iat4=0
          DO i=1,ions1%nat
             DO j=i+1,ions1%nat
                DO k=j+1,ions1%nat
                   DO l=k+1,ions1%nat
                      iat4=iat4+1
                      IF (ABS(abcd(iat4)).GT.abcm) THEN
                         i4max=i+(j-1)*ions1%nat+(k-1)*ions1%nat*ions1%nat+(l-1)*ions1%nat*ions1%nat*ions1%nat
                         iat4m=iat4
                         abcm=ABS(abcd(iat4))
                      ENDIF
                      IF (abcd(iat4).GT.prop3%cut4o) THEN
                         ip=ip+1
                         IF (paral%io_parent)&
                              WRITE(6,'(4A8,F10.4)')&
                              label(i),label(j),label(k),label(l),abcd(iat4)
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
          IF (ip.EQ.0.AND.abcm.GT.0.0_real_8) THEN
             ijk=i4max
             i = MOD(ijk-1,ions1%nat)+1
             k = ijk/(ions1%nat*ions1%nat)+1
             j = (ijk-(k-1)*ions1%nat*ions1%nat)/ions1%nat+1
             l = (ijk-(k-1)*ions1%nat*ions1%nat-(j-1)*ions1%nat*ions1%nat*ions1%nat)/ions1%nat+1
             IF (paral%io_parent)&
                  WRITE(6,'(4A8,F10.4)')&
                  label(i),label(j),label(k),label(l),abcd(iat4m)
          ENDIF
       ENDIF
    ENDIF
    ! ENDIF
3000 CONTINUE
    ! Clean up the memory
    DEALLOCATE(comp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(catom,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(cscr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(xsmat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(xxmat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(xtmat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (chunk_new.NE.0) THEN
       DEALLOCATE(bab,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    DEALLOCATE(qa,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(scr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (pslo_com%tivan) THEN
       DEALLOCATE(fnl2,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    CALL tihalt('    DIST_PROWFN',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE dist_prowfn
  ! ==================================================================
  SUBROUTINE dist_prtmat(xxmat,nattot,numorb,label,comp,f)
    ! CHANGED TO READ XXMAT DATA FROM FILE
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nattot
    REAL(real_8)                             :: xxmat(nattot,*)
    INTEGER                                  :: numorb
    CHARACTER(len=15)                        :: label(nattot)
    REAL(real_8)                             :: comp(*), f(*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'dist_prtmat'

    INTEGER                                  :: i1, i2, ia, ierr, ii, is, &
                                                isw, k, l, lnattot, msweep, &
                                                numaor(nattot)
    LOGICAL                                  :: ferror
    REAL(real_8)                             :: tmp
    REAL(real_8), ALLOCATABLE                :: lxxmat(:,:)

! Variables
! NEW VARIABLES
! ==--------------------------------------------------------------==

    !TK Bug fix: tmp only defined on io_parent...
    !io_parent only does something here anyway - move io_parent to the
    !very beginning
    IF(paral%io_parent)THEN
       msweep=numorb/8
       IF (MOD(numorb,8).NE.0) msweep=msweep+1
       ! ALLOCATE LOCAL MEMORY
       ALLOCATE(lxxmat(nattot,8),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       DO isw=1,msweep
          CALL fileopen(41,'WFNCOEF',fo_def+fo_ufo,ferror)
          READ(41) lnattot,ions1%nsp,(ions0%zv(is),ions0%na(is),numaor(is),is=1,ions1%nsp)
          i1=(isw-1)*8+1
          i2=MIN(8*isw,numorb)
          ! LOAD XXMAT(*,I1:I2) PART FROM FILE TO LXXMAT
          DO k=1, nattot
             is = 1
             DO l=1, numorb
                READ(41) tmp
                IF ((l.GE.i1).AND.(l.LE.i2)) THEN
                   lxxmat(k,is) = tmp
                   is = is +1
                ENDIF
             ENDDO
          ENDDO
          CALL fileclose(41)
          WRITE(6,'(/,A,8(I5,3X))') '      ORBITAL  ',(ii,ii=i1,i2)
          WRITE(6,'(A,8(F7.3,1X))')&
               '  COMPLETNESS  ',(comp(ii),ii=i1,i2)
          WRITE(6,'(A,8(F7.3,1X),/)')&
               '  OCCUPATION   ',(f(ii),ii=i1,i2)
          ! DO IA=1,NATTOT
          ! WRITE(6,'(A15,8(F7.3,1X))') LABEL(IA),
          ! *                                   (XXMAT(IA,II),II=I1,I2)
          DO ia=1,nattot
             WRITE(6,'(A15,8(F7.3,1X))') label(ia),&
                  (lxxmat(ia,ii),ii=1,8)
          ENDDO
       ENDDO
       ! FREE LOCAL MEMORY
       DEALLOCATE(lxxmat,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       ! ==--------------------------------------------------------------==
    END IF
    RETURN
  END SUBROUTINE dist_prtmat
  SUBROUTINE give_scr_dist_prowfn(lprowfn,tag,norbx)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lprowfn
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: norbx

    INTEGER                                  :: is, lcmaos, lrnlsm, lsatch, &
                                                lsummat, lwfnrho, nstate, &
                                                numin

    nstate=crge%n
    lrnlsm=0
    lcmaos=0
    lsatch=0
    lwfnrho=0
    lsummat=0
    ! CALL GIVE_SCR_SUMMAT(LSUMMAT,TAG,NATTOT)
    IF (pslo_com%tivan) CALL give_scr_rnlsm(lrnlsm,tag,atwp%nattot,.FALSE.)
    IF (prop1%dpan) THEN
       CALL stopgm('DIST_PROWFN','DISTRIBUTED DAVIDSON NOT IMPLEMENTED',& 
            __LINE__,__FILE__)
       ! Davidson Population Analysis
       numin=0
       DO is=1,ions1%nsp
          numin=numin+ions0%na(is)*prop2%maos(is)
       ENDDO
       CALL give_scr_cmaos(lcmaos,tag,numin)
       CALL give_scr_satch(lsatch,tag,numin)
    ENDIF
    ! PROWFN
    lprowfn=MAX(2*atwp%nattot+4*atwp%nattot,&
         atwp%nattot*prop2%numorb,&
         atwp%nattot*norbx,&
         lrnlsm,lsummat,lcmaos,lsatch,lwfnrho)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_dist_prowfn
  ! ==================================================================

END MODULE dist_prowfn_utils
