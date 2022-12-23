MODULE lr_tddft_utils
  USE canon_utils,                     ONLY: canon,&
                                             give_scr_canon
  USE corec_utils,                     ONLY: corec
  USE ddipo_utils,                     ONLY: give_scr_ddipo
  USE dftin_utils,                     ONLY: tdm_fun
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE forcedr_driver,                  ONLY: forcedr
  USE forcedr_utils,                   ONLY: give_scr_forcedr
  USE gettrans_utils,                  ONLY: gettrans
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE isos,                            ONLY: isos3
  USE kinds,                           ONLY: real_8
  USE linres,                          ONLY: &
       c1k, c_sh, ck, lr01, lr03, lrf4, lrsym, shfilen2, shlct, shsigma, &
       shsigmaf, td01, td03, tshi, tshl, v_sh, xfmqc
  USE lr_diag_utils,                   ONLY: give_scr_lr_diag,&
                                             lr_diag
  USE lr_tddft_drhoe,                  ONLY: print_ex_rho,&
                                             relax_rho
  USE lr_xcpot_utils,                  ONLY: lr_xcpot
  USE mm_dim_utils,                    ONLY: mm_dim
  USE mm_dimmod,                       ONLY: mm_go_mm,&
                                             mm_go_qm,&
                                             mm_revert
  USE mm_input,                        ONLY: lqmmm
  USE mm_qmmm_forcedr_utils,           ONLY: mm_qmmm_forcedr
  USE mm_rho_forcedr_utils,            ONLY: mm_rho_forcedr
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum
  USE nlcc,                            ONLY: corel,&
                                             roct
  USE ortho_utils,                     ONLY: give_scr_ortho
  USE ovlap_utils,                     ONLY: ovlap
  USE parac,                           ONLY: parai,&
                                             paral
  USE poin,                            ONLY: potr,&
                                             ptau,&
                                             rhoo,&
                                             rtau
  USE rhoofr_utils,                    ONLY: give_scr_rhoofr,&
                                             rhoofr
  USE rnlsm_utils,                     ONLY: give_scr_rnlsm
  USE ropt,                            ONLY: infi
  USE sh_tddft_utils,                  ONLY: &
       adjustphase, adjustphase_lr, ecoupl, getsigma, lz_trans, &
       normalize_c_sh, read_shcoeff, resize_f, td_spect, tmprd, tmpwr, &
       update_sh, write_shcoeff
  USE sh_utils,                        ONLY: do_tullytest,&
                                             propagate_c
  USE soc_types,                       ONLY: ene_sin,&
                                             ene_tri,&
                                             etot_isc_md,&
                                             md_on_sin,&
                                             md_on_tri,&
                                             nskip
  USE spin,                            ONLY: clsd,&
                                             spin_mod,&
                                             tdsp1
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             maxsys,&
                                             ncpw
  USE tauf,                            ONLY: tau,&
                                             vtau
  USE td_force_utils,                  ONLY: give_scr_td_force,&
                                             td_force
  USE td_mm_qmmm_forcedr_utils,        ONLY: td_mm_qmmm_forcedr
  USE td_os_utils,                     ONLY: td_os
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE tpot,                            ONLY: eigref,&
                                             eigval,&
                                             foccp,&
                                             mstate
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: lr_tddft
  PUBLIC :: give_scr_lr_tddft
  PUBLIC :: setspin

CONTAINS

  ! ==================================================================
  SUBROUTINE lr_tddft(c0,c1,c2,sc0,rhoe,psi,tau0,fion,eigv,&
       nstate,tfor,kprint)
    ! ==--------------------------------------------------------------==
    ! Qmmm
    COMPLEX(real_8)                          :: c0(:,:), c2(:,:), sc0(*)
    REAL(real_8)                             :: rhoe(:,:)
    COMPLEX(real_8)                          :: psi(:,:)
    REAL(real_8)                             :: tau0(:,:,:), fion(:,:,:), &
                                                eigv(*)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c1(ncpw%ngw,nstate)
    LOGICAL                                  :: tfor
    INTEGER                                  :: kprint

    CHARACTER(*), PARAMETER                  :: procedureN = 'lr_tddft'

    CHARACTER(len=10)                        :: orbital
    COMPLEX(real_8), ALLOCATABLE             :: cvirt(:,:), vtmp(:)
    INTEGER :: i, ia, ierr, is, isub, it, lddxc_1d, lddxc_2d, ms, nact, nap, &
      nbp, ncscr, nroot, ntddft, nvir, temp, transition(3,4,20)
    INTEGER, SAVE                            :: ifirst = 0
    LOGICAL                                  :: fexist, init_c1k, status, &
                                                statusdummy
    REAL(real_8), ALLOCATABLE                :: ddxc(:,:), edav(:), eigm(:), &
                                                eigs(:), ovlapm(:,:), &
                                                shprob(:)

! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    ! ==--------------------------------------------------------------==
    ! 
    IF (td03%treorder) THEN
       CALL stopgm("LR_TDDFT","NO REORDER",& 
            __LINE__,__FILE__)
    ENDIF
    ! ..TSH[
    init_c1k=.FALSE.
    nroot=MAX(td01%ns_mix,td01%ns_sin,td01%ns_tri,2)
    IF (tshl%tdtully.AND.ifirst.EQ.0) THEN
       ! done during the initialization
       tshl%s0_sh=.FALSE.
       IF (td01%fstate.EQ.0) tshl%s0_sh=.TRUE.
       tshl%tully_sh=.FALSE.
       ! this is done during the initialization of the forces 
       IF (paral%io_parent) INQUIRE(file=shfilen2,exist=fexist)
       CALL mp_bcast(fexist,parai%io_source,parai%cp_grp)
       IF (.NOT.fexist) THEN
          init_c1k=.TRUE.
       ELSE
          CALL tmprd(c1k,nstate,nroot,temp,shfilen2)
       ENDIF
    ENDIF

    IF (tshl%tdtully.AND.(.NOT.tshl%txfmqc).AND.ifirst.EQ.1) THEN
       ! done during the first cntl%md step
       IF (tshi%shstep.EQ.1) THEN
          ! INITIALIZE_SH
          tshl%s0_sh=.FALSE.
          IF (td01%fstate.EQ.0) tshl%s0_sh=.TRUE.
          tshl%tully_sh=.FALSE.
          CALL zeroing(c_sh)!,2*(nroot+1))
          CALL zeroing(v_sh)!,2*(nroot+1))
          CALL zeroing(shsigma)!,2*(nroot+1)*(nroot+1))
          CALL zeroing(shsigmaf)!,2*(nroot+1)*(nroot+1))
          c_sh(1,td01%fstate+1)=(1._real_8,0._real_8)
          c_sh(2,td01%fstate+1)=(1._real_8,0._real_8)
          ! FSTATE from input file &cntl%tddft
       ELSE
          tshl%s0_sh=.FALSE.
          tshl%tully_sh=.FALSE.
          CALL read_shcoeff
          tshi%shstep=tshi%shstep+1
       ENDIF
    ENDIF
    IF (tshl%tdtully.AND.tshl%txfmqc.AND.ifirst.EQ.0) THEN
       ! done during the first cntl%md step
       IF (tshi%shstep.EQ.0) THEN
          ! INITIALIZE_SH
          tshl%s0_sh=.FALSE.
          IF (td01%fstate.EQ.0) tshl%s0_sh=.TRUE.
          tshl%tully_sh=.FALSE. 
          CALL zeroing(c_sh)!,2*(nroot+1))
          CALL zeroing(v_sh)!,2*(nroot+1)) 
          CALL zeroing(shsigma)!,2*(nroot+1)*(nroot+1))
          c_sh(1,td01%fstate+1)=(1._real_8,0._real_8)
          c_sh(2,td01%fstate+1)=(1._real_8,0._real_8)
          !WRITE(6,*) 'FEDE', c_sh(1,td01%fstate+1)
          !STOP
          ! FSTATE from input file &cntl%tddft 
       ELSE
          tshl%s0_sh=.FALSE.
          tshl%tully_sh=.FALSE.
          CALL read_shcoeff
          tshi%shstep=tshi%shstep+1
       ENDIF
    ENDIF
    ! 
    IF (td03%tdlz.OR.tshl%tdtully) THEN
       IF (paral%qmnode) CALL mm_dim(mm_go_qm,statusdummy)
       nroot=MAX(td01%ns_mix,td01%ns_sin,td01%ns_tri,2)
       nvir=nroot
       ncscr=2*ncpw%ngw*nvir
       ntddft=nstate+nvir
       tshl%tully_sh=.FALSE.
       ALLOCATE(cvirt(ncpw%ngw,ncscr/ncpw%ngw),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(edav(ntddft),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(shprob((nroot+1)*(nroot+1)),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       IF (td03%sh_diag) THEN
          IF (ifirst.EQ.0) THEN
             CALL resize_f
          ENDIF
          CALL setspin("S",c0,cvirt,c2)
          CALL td_spect(c0,c2,rhoe,ntddft,cvirt,sc0,psi(:,1),edav)
          CALL setspin("R",c0,cvirt,c2)
          IF (cntl%tlsd) THEN
             nap=spin_mod%nsup-tdsp1%nupel
             nbp=spin_mod%nsdown-tdsp1%ndoel
             CALL dcopy(2*ncpw%ngw*nap,c0(1,tdsp1%nupel+1),1,cvirt(1,1),1)
             CALL dcopy(2*ncpw%ngw*nbp,c0(1,spin_mod%nsup+tdsp1%ndoel+1),1,cvirt(1,nap+1),1)
          ELSE
             CALL dcopy(2*ncpw%ngw*nvir,c0(1,crge%n+1),1,cvirt(1,1),1)
          ENDIF
       ENDIF
       ALLOCATE(ovlapm(tshi%nshs,tshi%nshs),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ! IF (PARENT) THEN
       ! write(6,'(A,3I7)') 'NSTATE,NROOT,NSHS',NSTATE,NROOT,NSHS
       ! ENDIF
    ENDIF
    ! ..TSH]
    ! calculate Kohn-Sham matrix
    ALLOCATE(eigm(nstate*nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    IF ( lqmmm%qmmm ) THEN
       CALL mm_dim(mm_go_qm,statusdummy)
       IF (cntl%tddft) THEN
          CALL td_mm_qmmm_forcedr(c0,c2,sc0,rhoe,psi,tau0,fion,&
               eigv,nstate,1,.FALSE.,.TRUE.,.TRUE.)
       ELSE
          CALL mm_qmmm_forcedr(c0,c2,sc0,rhoe,psi,tau0,fion,&
               eigv,nstate,1,.FALSE.,.FALSE.,.TRUE.)
       ENDIF
    ELSE
       CALL forcedr(c0(:,1:nstate),c2(:,1:nstate),sc0,rhoe,psi,tau0,fion,eigv,&
            nstate,1,.FALSE.,tfor)
    ENDIF
    ! transform to canonical orbitals and get eigenvalues
    CALL canon(c0,c2,crge%f,nstate,eigv)
    IF (lrf4%td_method.EQ.1) THEN
       orbital="GENERAL"
       CALL tdm_fun(lrf4%td_functional,1)
       IF (cntl%tpotential) THEN
          ! 
          ALLOCATE(eigval(nstate),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(eigref(nstate),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(foccp(nstate),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          CALL dcopy(nstate,crge%f,1,foccp,1)
          CALL dcopy(nstate,eigv,1,eigval,1)
          mstate=nstate
          IF (td03%molstat) THEN
             CALL stopgm("LR_TDDFT","MOLSTAT NOT IMPLEMENTED",& 
                  __LINE__,__FILE__)
          ENDIF
       ENDIF
       IF ( lqmmm%qmmm ) THEN
          CALL mm_dim(mm_go_qm,statusdummy)
          CALL mm_qmmm_forcedr(c0,c2,sc0,rhoe,psi,tau0,fion,&
               eigv,nstate,1,.FALSE.,.FALSE.,.TRUE.)
       ELSE
          CALL forcedr(c0(:,1:nstate),c2(:,1:nstate),sc0,rhoe,psi,tau0,fion,eigv,&
               nstate,1,.FALSE.,.FALSE.)
       ENDIF
       DO is=1,nstate
          CALL dscal(2*ncpw%ngw,1._real_8/crge%f(is,1),c2(1,is),1)
       ENDDO
       CALL ovlap(nstate,eigm,c0,c2)
       CALL mp_sum(eigm,nstate*nstate,parai%allgrp)
    ELSE
       orbital="CANON"
       CALL dcopy(nstate,eigv,1,eigm,1)
    ENDIF
    ! ..TSH[
    IF (td03%tdlz.OR.tshl%tdtully) THEN
       IF (cntl%tlsd) THEN
          nap=spin_mod%nsup-tdsp1%nupel
          nbp=spin_mod%nsdown-tdsp1%ndoel
          CALL dcopy(2*ncpw%ngw*nap,c0(1,tdsp1%nupel+1),1,cvirt(1,1),1)
          CALL dcopy(2*ncpw%ngw*nbp,c0(1,spin_mod%nsup+tdsp1%ndoel+1),1,cvirt(1,nap+1),1)
       ELSE
          CALL dcopy(2*ncpw%ngw*nvir,c0(1,crge%n+1),1,cvirt(1,1),1)
       ENDIF
    ENDIF
    ! ..TSH]
    ! store potential and density
    CALL dcopy(fpar%nnr1*clsd%nlsd,rhoe,1,potr,1)
    CALL rhoofr(c0(:,1:nstate),rhoe,psi(:,1),nstate)
    CALL dcopy(fpar%nnr1*clsd%nlsd,rhoe,1,rhoo,1)

    IF (corel%tinlc) THEN
       ALLOCATE(roct(fpar%nnr1*clsd%nlsd),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(roct)!,nnr1*clsd%nlsd)
       ALLOCATE(vtmp(ncpw%nhg),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       CALL corec(roct,vtmp,psi)
       DEALLOCATE(vtmp,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (isos3%ps_type.EQ.1) CALL stopgm('LR_TDDFT','HOCKNEY PS not impl.',& 
         __LINE__,__FILE__)
    ms=MAX(td01%ns_mix,td01%ns_sin,td01%ns_tri)
    ALLOCATE(eigs(ms),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! calculate and store analytic dmu/dn
    IF (lr03%txc_analytic) THEN
       lddxc_1d=fpar%nnr1
       lddxc_2d=(2*clsd%nlsd-1)
    ELSE
       lddxc_1d=1
       lddxc_2d=1
    ENDIF
    IF (lr01%lopti.GE.0) lr01%lopti=lr01%lopti-10
    ALLOCATE(ddxc(lddxc_1d,lddxc_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL lr_xcpot(ddxc,rhoo,.TRUE.)
    IF (cntl%tlsd) THEN
       lrsym=0
    ELSEIF (td01%ns_sin.NE.0) THEN
       lrsym=1
    ELSEIF (td01%ns_tri.NE.0) THEN
       lrsym=3
    ENDIF
    IF (cntl%ttau) THEN
       CALL dcopy(fpar%nnr1*clsd%nlsd,rtau,1,tau,1)
       CALL dcopy(fpar%nnr1*clsd%nlsd,ptau,1,vtau,1)
    ENDIF
    CALL lr_diag(c0,c1,c2,sc0,eigm,ddxc,rhoe,psi,&
         eigs,nstate,ms,orbital,kprint)
    CALL mp_bcast(eigs,SIZE(eigs),parai%io_source,parai%cp_grp)
    ! ==--------------------------------------------------------------==
    CALL td_os(nstate,c0,ms,eigs,c1)
    ! ==--------------------------------------------------------------==
    DEALLOCATE(ddxc,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    IF (td01%fstate.GT.0) THEN
       ener_com%etddft=eigs(td01%fstate)
    ELSE
       ener_com%etddft=0.0_real_8
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ..TSH[
    IF (td03%tdlz.OR.tshl%tdtully) THEN
       IF (tshl%isc) THEN 
          ! collect energies for energy derivatives in LZ formula
          IF (.NOT.tshl%isc_2ndlr.AND.(MOD(infi,nskip).EQ.0)) THEN
             IF (md_on_sin) THEN
                DO it=10,2,-1
                   DO is=1,nroot+1
                      ene_sin(it,is)=ene_sin(it-1,is)
                   ENDDO
                ENDDO
                DO is=1,nroot
                   ene_sin(1,is+1)=eigs(is)
                ENDDO
             ELSE
                DO it=10,2,-1
                   DO is=1,nroot+1
                      ene_tri(it,is)=ene_tri(it-1,is)
                   ENDDO
                ENDDO
                DO is=1,nroot
                   ene_tri(1,is+1)=eigs(is)
                ENDDO
             ENDIF
          ELSEIF(tshl%isc_2ndlr.AND.(MOD(infi,nskip).EQ.0)) THEN 
             IF (md_on_sin) THEN
                DO it=10,2,-1
                   DO is=1,nroot+1
                      ene_tri(it,is)=ene_tri(it-1,is)
                   ENDDO
                ENDDO
                DO is=1,nroot
                   ene_tri(1,is+1)=eigs(is)
                ENDDO
             ELSE
                DO it=10,2,-1
                   DO is=1,nroot+1
                      ene_sin(it,is)=ene_sin(it-1,is)
                   ENDDO
                ENDDO
                DO is=1,nroot
                   ene_sin(1,is+1)=eigs(is)
                ENDDO
             ENDIF
          ENDIF
          ! set ground state excitation energy to 0._real_8
          DO it=1,10
             ene_sin(it,1)=0._real_8
             ene_tri(it,1)=0._real_8
          ENDDO
       ENDIF
       !---------------------------------------------------------------------
       !if (paral%io_parent) then
       !  write(*,'(/,10F8.3)') (ene_sin(1,it),it=1,nroot)
       !  write(*,'(10F8.3)')   (ene_sin(2,it),it=1,nroot)
       !  write(*,'(10F8.3)')   (ene_tri(1,it),it=1,nroot)
       !  write(*,'(10F8.3,/)') (ene_tri(2,it),it=1,nroot)
       !endif
       !-------------------------------------------------------------------- 

       IF ((.NOT.tshl%isc) &
            .OR. (tshl%isc.AND.(.NOT.tshl%isc_2ndlr))) THEN
          IF (td03%tda.AND.td01%msubta.GT.0) THEN
             nact=td01%msubta
          ELSE
             IF (cntl%tlsd) THEN
                nact=tdsp1%nupel+tdsp1%ndoel
             ELSE
                nact=crge%n
             ENDIF
          ENDIF
          CALL zeroing(transition)!,3*4*20)
          CALL gettrans(transition,c1,cvirt,&
               nact,nvir,nap,nbp,nroot)
          IF (tshl%tdtully) THEN
             IF (tshi%shstep.EQ.0) CALL dcopy(2*ncpw%ngw*tshi%nshs,c0(1,1),1,ck(1,1),1)
             IF (init_c1k) CALL dcopy(2*ncpw%ngw*nstate*nroot,c1,1,c1k,1)
             DO is=1,nroot
                v_sh(1,is+1)=eigs(is)
                IF (tshl%txfmqc) THEN
                   xfmqc%eigv(is)=eigs(is)
                ENDIF
                ! if(parent) write(6,*) 'energies',V_SH(1,IS+1)
             ENDDO
             v_sh(1,1)=0._real_8
          ENDIF
          IF (ifirst.GT.0) THEN
             IF (td03%tdlz) THEN
                CALL lz_trans(transition,eigs,nact,nvir,nap,nbp,nroot)
             ELSEIF (tshl%tdtully.AND.tshi%shstep.GE.1) THEN
                IF (tshl%tdextpot) CALL ecoupl(nstate,c0,nroot,c1,eigs)
                IF (shlct%sh_lcontrol) CALL ecoupl(nstate,c0,nroot,c1,eigs)
                CALL zeroing(ovlapm)!,tshi%nshs*tshi%nshs)
                CALL adjustphase_lr(c1,c1k,nstate,nroot,ovlapm)
                CALL zeroing(ovlapm)!,tshi%nshs*tshi%nshs)
                CALL adjustphase(c0,ck,tshi%nshs,ovlapm)
                CALL getsigma(c0,c1,nstate,nroot)
                IF (tshl%tdextpot) CALL ecoupl(nstate,c0,nroot,c1,eigs)
                CALL propagate_c(nstate,nroot+1,shprob)
                CALL do_tullytest(nroot+1,shprob)
                CALL normalize_c_sh(nroot)
                CALL write_shcoeff(shprob,ener_com%etot,nroot)
                CALL tmpwr(c1,nstate,nroot,tshi%shstep,shfilen2)
                ! NEWinOLD 
                CALL update_sh(nroot)
                CALL dcopy(2*ncpw%ngw*nstate*nroot,c1,1,c1k,1)
             ENDIF
          ENDIF
          ifirst=ifirst+1
       ELSE ! tshl%isc.and.tshl%isc_2ndlr
          IF (td03%tda.AND.td01%msubta.GT.0) THEN
             nact=td01%msubta
          ELSE
             IF (cntl%tlsd) THEN
                nact=tdsp1%nupel+tdsp1%ndoel
             ELSE
                nact=crge%n
             ENDIF
          ENDIF
          CALL zeroing(transition)!,3*4*20)

          CALL gettrans(transition,c1,cvirt,nact,nvir,nap,nbp,nroot)
          !     CALL printtrans(transition,eigs,nact,nvir,nap,nbp,nroot)
          IF (tshl%isc.AND.md_on_sin) THEN
             IF(paral%io_parent) WRITE(6,'(/,1x,"ISC:",60("-"))') 
             DO is=1,nroot
                IF(paral%io_parent)  &
                     WRITE(6,'(A,I5,F15.7)') ' Triplet energies', is, eigs(is)
             ENDDO
             IF(paral%io_parent) WRITE(6,'(1x,"ISC:",60("-"),/)')
          ELSEIF (tshl%isc.AND.md_on_tri) THEN
             IF(paral%io_parent) WRITE(6,'(1x,"ISC:",60("-"))')
             DO is=1,nroot
                IF(paral%io_parent)  &
                     WRITE(6,'(A,I5,F15.7)') ' Singlet energies', is, eigs(is)
             ENDDO
             IF(paral%io_parent) WRITE(6,'(1x,"ISC:",60("-"),/)')
          ENDIF
       ENDIF
    ENDIF
    ! ..TSH]
    ! ==--------------------------------------------------------------==
    IF (tfor) THEN
       IF (tshl%isc) etot_isc_md=ener_com%etot
       ener_com%etot=ener_com%etot+ener_com%etddft
       ! ..TSH[
       IF (tshl%tdtully.AND.(tshl%s0_sh)) THEN
          IF (paral%io_parent) WRITE(6,*) 'RUNNING IN THE GROUND STATE'
          ! do nothing
       ELSE
          IF (tshl%txfmqc) THEN
             IF (paral%io_parent) WRITE(6,*) '...ADD EXCITED STATES FORCES'
             td01%store_fstate=td01%fstate
             ! cp fion into fion_state(:,:,:,1) for the  gs forces
             CALL dcopy(3*maxsys%nax*maxsys%nsx,fion(1,1,1),1,xfmqc%fion_state(1,1,1,1),1)
             DO i=1,tshi%nroot_sh
                !  azzero fion
                CALL zeroing(fion)
                td01%fstate=i
                CALL td_force(c0,c1,c2,sc0,eigm,rhoe,psi,eigs,&
                     nstate,ms,tau0,fion,orbital,kprint)
                ! add fion to fion_state(1) and put it into fion_state(fstate+1)
                DO is = 1,ions1%nsp
                   DO ia = 1,ions0%na(is)
                      xfmqc%fion_state(1,ia,is,td01%fstate+1)=&
                           xfmqc%fion_state(1,ia,is,1) + fion(1,ia,is) 
                      xfmqc%fion_state(2,ia,is,td01%fstate+1)=&
                           xfmqc%fion_state(2,ia,is,1) + fion(2,ia,is) 
                      xfmqc%fion_state(3,ia,is,td01%fstate+1)=&
                           xfmqc%fion_state(3,ia,is,1) + fion(3,ia,is) 
                   ENDDO
                ENDDO
             ENDDO
             ! ivano check this
             CALL dcopy(3*maxsys%nax*maxsys%nsx,xfmqc%fion_state(1,1,1,1),1,fion(1,1,1),1)
             td01%fstate=td01%store_fstate
          ELSE
             CALL td_force(c0,c1,c2,sc0,eigm,rhoe,psi,eigs,&
                  nstate,ms,tau0,fion,orbital,kprint)
             CALL relax_rho(c0,c1,c2,sc0,rhoo,rhoe,psi(:,1),nstate,ms)
             CALL print_ex_rho(rhoe,psi(:,1),tau0)
          ENDIF
       ENDIF
       ! ..TSH]
       ! ..TSH[
       IF (tshl%tdtully.AND.lqmmm%qmmm) THEN
          CALL mm_dim(mm_go_mm,status)
          CALL mm_rho_forcedr(rhoe,tau0,fion,nstate,.TRUE.)
          CALL mm_dim(mm_revert,status)
       ENDIF
       ! ..TSH]
    ENDIF
    IF ((lqmmm%qmmm .AND. cntl%tddft .AND. paral%parent).AND.paral%io_parent)&
         WRITE(41,'(3F15.9)') ener_com%etot-ener_com%etddft,ener_com%etot,ener_com%etddft
    ! ==--------------------------------------------------------------==
    IF ((lrf4%td_method.EQ.1).AND.(cntl%tpotential)) THEN
       CALL tdm_fun(lrf4%ks_functional,1)
       DEALLOCATE(eigval,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(eigref,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(foccp,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    DEALLOCATE(eigs,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(eigm,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ..TSH[
    IF (td03%tdlz.OR.tshl%tdtully) THEN
       DEALLOCATE(cvirt,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(edav,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    IF (tshl%tdtully) THEN
       DEALLOCATE(ovlapm,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(shprob,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! ..TSH]
    IF (corel%tinlc) DEALLOCATE(roct,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE lr_tddft
  ! ==================================================================
  SUBROUTINE give_scr_lr_tddft(lspectra,tfor,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lspectra
    LOGICAL                                  :: tfor
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: l_diag, l_for, lcanon, &
                                                lforces, llocal, llr_diag, &
                                                lortho, lrhoofr, lrnlsm, &
                                                nstate

    l_diag=0
    nstate=crge%n
    CALL give_scr_rhoofr(lrhoofr,tag)
    CALL give_scr_forcedr(lforces,tag,nstate,.FALSE.,tfor)
    CALL give_scr_canon(lcanon,tag,nstate)
    CALL give_scr_rnlsm(lrnlsm,tag,nstate,.FALSE.)
    CALL give_scr_lr_diag(llr_diag,tag)
    nstate=crge%n+MAX(td01%ns_tri,td01%ns_sin,td01%ns_mix)
    CALL give_scr_ortho(lortho,tag,nstate)
    IF (td03%tdlocal) THEN
       CALL give_scr_ddipo(llocal,tag)
       llocal=MAX(llocal,2*ncpw%nhg)
    ELSE
       llocal=0
    ENDIF
    l_for=0
    IF (tfor) THEN
       CALL give_scr_td_force(l_for,nstate,tag)
    ENDIF
    IF (corel%tinlc) THEN
       llocal=MAX(llocal,ncpw%nhg)
    ENDIF
    ! 
    lspectra=MAX(lrhoofr,lforces,lcanon,lrnlsm,&
         llr_diag,l_diag,lortho,llocal,l_for)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE give_scr_lr_tddft
  ! ==========================================================

  SUBROUTINE setspin(tag,c0,cvirt,cs)
    ! ==--------------------------------------------------------------==
    CHARACTER(len=1)                         :: tag
    COMPLEX(real_8)                          :: c0(ncpw%ngw,*), &
                                                cvirt(ncpw%ngw,*), &
                                                cs(ncpw%ngw,*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'setspin'

    INTEGER                                  :: i, isub, nap, nbp, nx

! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    IF (cntl%tlsd) THEN
       IF (tag.EQ."S") THEN
          nap=spin_mod%nsup-tdsp1%nupel
          nbp=spin_mod%nsdown-tdsp1%ndoel
          CALL dcopy(2*ncpw%ngw*nap,c0(1,tdsp1%nupel+1),1,cvirt(1,1),1)
          CALL dcopy(2*ncpw%ngw*nbp,c0(1,spin_mod%nsup+tdsp1%ndoel+1),1,cvirt(1,nap+1),1)
          CALL dcopy(2*ncpw%ngw*tdsp1%ndoel,c0(1,spin_mod%nsup+1),1,cs(1,1),1)
          CALL dcopy(2*ncpw%ngw*tdsp1%ndoel,cs(1,1),1,c0(1,tdsp1%nupel+1),1)
          !$omp parallel do private(I)
#ifdef __SR8000
          !poption parallel
#endif
#ifdef _vpp_
          !OCL NOALIAS
#endif
          DO i=1,spin_mod%nsup+spin_mod%nsdown
             crge%f(i,1)=0._real_8
          ENDDO
          spin_mod%nsup=tdsp1%nupel
          spin_mod%nsdown=tdsp1%ndoel
          !$omp parallel do private(I)
#ifdef __SR8000
          !poption parallel
#endif
#ifdef _vpp_
          !OCL NOALIAS
#endif
          DO i=1,tdsp1%nupel+tdsp1%ndoel
             crge%f(i,1)=1._real_8
          ENDDO
       ELSEIF (tag.EQ."R") THEN
          nx=crge%n-tdsp1%nupel-tdsp1%ndoel
          spin_mod%nsup=tdsp1%nupel+(nx+1)/2
          spin_mod%nsdown=tdsp1%ndoel+nx/2
          !$omp parallel do private(I)
#ifdef __SR8000
          !poption parallel
#endif
#ifdef _vpp_
          !OCL NOALIAS
#endif
          DO i=1,crge%n
             crge%f(i,1)=0._real_8
          ENDDO
          !$omp parallel do private(I)
#ifdef __SR8000
          !poption parallel
#endif
#ifdef _vpp_
          !OCL NOALIAS
#endif
          DO i=1,tdsp1%nupel
             crge%f(i,1)=1._real_8
          ENDDO
          !$omp parallel do private(I)
#ifdef __SR8000
          !poption parallel
#endif
#ifdef _vpp_
          !OCL NOALIAS
#endif
          DO i=spin_mod%nsup+1,spin_mod%nsup+tdsp1%ndoel
             crge%f(i,1)=1._real_8
          ENDDO
          nap=spin_mod%nsup-tdsp1%nupel
          nbp=spin_mod%nsdown-tdsp1%ndoel
          CALL dcopy(2*ncpw%ngw*tdsp1%ndoel,c0(1,tdsp1%nupel+1),1,cs(1,1),1)
          CALL dcopy(2*ncpw%ngw*tdsp1%ndoel,cs(1,1),1,c0(1,spin_mod%nsup+1),1)
          CALL dcopy(2*ncpw%ngw*nap,cvirt(1,1),1,c0(1,tdsp1%nupel+1),1)
          CALL dcopy(2*ncpw%ngw*nbp,cvirt(1,nap+1),1,c0(1,spin_mod%nsup+tdsp1%ndoel+1),1)
       ENDIF
    ELSE
       ! Nothing to do
    ENDIF
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE setspin
  ! ==================================================================

END MODULE lr_tddft_utils
