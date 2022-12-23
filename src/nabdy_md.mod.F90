MODULE nabdy_md
  USE andp,                            ONLY: rin0,&
                                             rmix,&
                                             rout0
  USE andr,                            ONLY: andr2
  USE anneal_utils,                    ONLY: anneal,&
                                             berendsen,&
                                             dampdyn,&
                                             tempramp
  USE atwf,                            ONLY: tmovr
  USE box_boundary_utils,              ONLY: box_boundary
  USE bs_forces_diag_utils,            ONLY: bs_forces_diag
  USE bsym,                            ONLY: bsclcs,&
                                             bsfac
  USE calc_alm_utils,                  ONLY: calc_alm
  USE cdft_utils,                      ONLY: cdft_reinit,&
                                             cdft_w,&
                                             init_cdft,&
                                             write_w
  USE cdftmod,                         ONLY: cdftcom,&
                                             wdiff
  USE cnst,                            ONLY: factem,&
                                             pi
  USE cnst_dyn,                        ONLY: ekincv,&
                                             vharm
  USE comvel_utils,                    ONLY: comvel
  USE comvelmod,                       ONLY: comvl
  USE coor,                            ONLY: fion,&
                                             tau0,&
                                             taup,&
                                             velp
  USE copot_utils,                     ONLY: copot
  USE cotr,                            ONLY: cotc0
  USE detdof_utils,                    ONLY: detdof
  USE dispp_utils,                     ONLY: dispp
  USE dynit_utils,                     ONLY: dynit
  USE ekinpp_utils,                    ONLY: ekinpp
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE extrap_utils,                    ONLY: extrap
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_app,&
                                             fo_nofp,&
                                             fo_verb
  USE finalp_utils,                    ONLY: finalp
  USE fint,                            ONLY: fint1
  USE forcep_utils,                    ONLY: rhoe_psi_size
  USE forces_diag_utils,               ONLY: forces_diag
  USE geofile_utils,                   ONLY: geofile
  USE gle_utils,                       ONLY: gle_init,&
                                             gle_step
  USE glemod,                          ONLY: glepar
  USE gsize_utils,                     ONLY: gsize
  USE hfxmod,                          ONLY: hfxc3
  USE initrun_driver,                  ONLY: initrun
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE isos,                            ONLY: isos1
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE linres,                          ONLY: td01
  USE localize_utils,                  ONLY: localize2
  USE lr_tddft_utils,                  ONLY: lr_tddft
  USE machine,                         ONLY: m_flush,&
                                             m_walltime
  USE md_driver,                       ONLY: give_scr_mddiag
  USE mddiag_interaction_p_utils,      ONLY: mddiag_interaction_p
  USE meta_multiple_walkers_utils,     ONLY: mw_filename
  USE mfep,                            ONLY: mfepi
  USE mm_dim_utils,                    ONLY: mm_dim
  USE mm_dimmod,                       ONLY: clsaabox,&
                                             cpat,&
                                             cpsp,&
                                             mm_go_mm,&
                                             mm_revert,&
                                             nat_grm
  USE mm_extrap,                       ONLY: cold, nnow, numcold
  USE moverho_utils,                   ONLY: moverho
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sync
  USE nabdy_ampli,                     ONLY: nabdy_action,&
                                             nabdy_initialize_action,&
                                             nabdy_keep_classical,&
                                             nabdy_repartition_ampli
  USE nabdy_forces,                    ONLY: nabdy_forces_calc
  USE nabdy_initialize,                ONLY: mom_to_vel,&
                                             nabdy_load,&
                                             nabdy_mem,&
                                             nabdy_rel,&
                                             vel_to_mom
  USE nabdy_types,                     ONLY: nabdyfric,&
                                             nabdyvar,&
                                             nabdyvec
  USE nlcc,                            ONLY: corel
  USE norm,                            ONLY: gnmax,&
                                             gnorm
  USE nose,                            ONLY: glib
  USE noseng_utils,                    ONLY: noseng
  USE nosepa_utils,                    ONLY: nosepa
  USE noseup_utils,                    ONLY: noseup
  USE nospinit_utils,                  ONLY: nospinit
  USE parac,                           ONLY: parai,&
                                             paral
  USE phfac_utils,                     ONLY: phfac
  USE pimd,                            ONLY: ipcurr
  USE poin,                            ONLY: potr,&
                                             rhoo
  USE posupi_utils,                    ONLY: posupi
  USE printave_utils,                  ONLY: paccc
  USE printp_utils,                    ONLY: printp,&
                                             printp2
  USE prng_utils,                      ONLY: repprngu
  USE proppt_utils,                    ONLY: propcal
  USE puttau_utils,                    ONLY: taucl
  USE rattle_utils,                    ONLY: rattle
  USE readsr_utils,                    ONLY: xstring
  USE resetac_utils,                   ONLY: resetac
  USE response_pmod,                   ONLY: dmbi
  USE rhopri_utils,                    ONLY: rhopri
  USE rinvel_utils,                    ONLY: rinvel,&
                                             rvscal
  USE rmas,                            ONLY: rmass
  USE ropt,                            ONLY: infi,&
                                             iteropt,&
                                             ropt_mod
  USE rotvel_utils,                    ONLY: rotvel
  USE rscvp_utils,                     ONLY: rscvp
  USE sample_utils,                    ONLY: sample_go,&
                                             sample_wait
  USE setbsstate_utils,                ONLY: setbsstate
  USE setirec_utils,                   ONLY: write_irec
  USE shake_utils,                     ONLY: cpmdshake
  USE soft,                            ONLY: soft_com
  USE spin,                            ONLY: clsd
  USE store_types,                     ONLY: &
       cprint, irec_ac, irec_nop1, irec_nop2, irec_nop3, irec_nop4, irec_vel, &
       restart1, rout1
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             cntr,&
                                             fpar,&
                                             maxsys,&
                                             nacc,&
                                             nkpt,&
                                             restf
  USE testex_utils,                    ONLY: testex
  USE teststore_utils,                 ONLY: teststore
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE totstr_utils,                    ONLY: totstr
  USE tpar,                            ONLY: dt_ions
  USE vdwcmod,                         ONLY: vdwl,&
                                             vdwwfl
  USE velupi_utils,                    ONLY: velupi
  USE wrener_utils,                    ONLY: wrprint_md
  USE wrgeo_utils,                     ONLY: wrgeo,&
                                             wrgeof
  USE wv30_utils,                      ONLY: zhwwf
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: nabdy_dyn
  PUBLIC :: nabdy_normalize

CONTAINS

  ! ==================================================================
  SUBROUTINE nabdy_dyn(c0,cm,c1,c2,sc0,vpp,gamx,gamy)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8), INTENT(inout) :: c0(nkpt%ngwk,crge%n,nkpt%nkpts), cm(:), &
      c1(*), c2(:,:), sc0(:)
    REAL(real_8), INTENT(inout)              :: vpp(:), gamx(:), gamy(:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'nabdy_dyn'

    CHARACTER(len=100)                       :: filen, ftmp
    CHARACTER(LEN=30)                        :: tag
    COMPLEX(real_8), ALLOCATABLE             :: c0mem(:,:,:), c2mem(:,:,:), &
                                                psi(:,:)
    INTEGER :: i1, i2, ia, ierr, ifcalc, il_psi, il_rhoe, irec(100), is, &
      isub, itemp, itraj, k, lenext, loopnfi, lscr, nfimax, nfimin, nmm, nnx, &
      nstate, nx, save_nfi
    INTEGER, SAVE                            :: ifirst = 0
    LOGICAL                                  :: ferror, fexist, qmmech, tstrng
    REAL(real_8) :: alfap, deltaeava, deltaeprime, disa, dummy(2), econs, &
      ekin1, ekin2, ekincp, ekinh1, ekinh2, ekinp, ekprime, enose, enosp, &
      etotprime, lmio(3), mm_ekin, mm_temp, natotene_diff, natotene_ref, &
      natotener(4), rmem, tcpu, temp1, temp2, tempp, time1, time2, totamplfe, &
      vcmio(4)
    REAL(real_8), ALLOCATABLE :: eigv(:,:), fionmem(:,:,:,:), rhoe(:,:), &
      rhoemem(:,:,:), rinp(:), rm1(:), scr(:), taui(:,:,:), tauio(:,:), &
      taur(:,:,:), velpmem(:,:,:,:)

! ==================================================================

    CALL tiset(procedureN,isub)
    IF (cntl%tddft.AND.cntl%tresponse) CALL stopgm("MDDIAG",&
         "cntl%tddft.AND.cntl%tresponse NOT POSSIBLE",& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    time1 =m_walltime()
    ALLOCATE(taup(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(fion(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(taui(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(taur(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(c0mem(nkpt%ngwk,crge%n,nabdyvar%ntrajbd),STAT=ierr)  
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(c2mem(nkpt%ngwk,crge%n,nabdyvar%ntrajbd),STAT=ierr)  
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(rhoemem(fpar%nnr1,clsd%nlsd,nabdyvar%ntrajbd),STAT=ierr)  
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(velpmem(3,maxsys%nax,maxsys%nsx,nabdyvar%ntrajbd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(fionmem(3,maxsys%nax,maxsys%nsx,nabdyvar%ntrajbd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(taup)!,SIZE(taup))
    CALL zeroing(fion)!,SIZE(fion))
    IF (comvl%tsubrot) THEN
       ALLOCATE(tauio(3,ions1%nat),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)! TODO check dimensions
    ENDIF
    ! Memory for densities
    IF (.NOT.cntl%bsymm) THEN
       nnx=fpar%nnr1*clsd%nlsd
       ALLOCATE(rin0(fpar%nnr1,clsd%nlsd),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(rout0(fpar%nnr1,clsd%nlsd),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(rmix(fpar%nnr1,clsd%nlsd),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(rm1(nnx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(rinp(nnx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)

       rhoo => rin0
       CALL zeroing(rin0)!,nnx)
       CALL zeroing(rout0)!,nnx)
       CALL zeroing(rmix)!,nnx)
    ENDIF
    nstate=crge%n
    nacc = 22
    iteropt%nfi  = 0
    ropt_mod%modens=.FALSE.
    ropt_mod%engpri=.FALSE.
    ropt_mod%calste=cntl%tpres
    ALLOCATE(eigv(crge%n,bsfac*nkpt%nkpts*clsd%nlsd*nstate/crge%n),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! Extrapolation
    IF (cntl%textrap) THEN
       lenext=2*nkpt%ngwk*nstate*nkpt%nkpts*cnti%mextra
       rmem = 16._real_8*lenext*1.e-6_real_8
       ALLOCATE(cold(nkpt%ngwk,crge%n,nkpt%nkpnt,lenext/(crge%n*nkpt%ngwk*nkpt%nkpnt)),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)

       IF (paral%io_parent)&
            WRITE(6,'(A,T51,F8.3,A)') ' MDDIAG| '&
            // 'EXTRAPOLATION WAVEFUNCTION HISTORY TAKES ',rmem,' MBYTES'
    ENDIF
    ! ==--------------------------------------------------------------==
    ! SCR ALLOCATION AND PARTITION (SCRATCH ARRAY).
    CALL rhoe_psi_size(il_rhoe,il_psi)
    ALLOCATE(rhoe(il_rhoe,1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    nmm=1
    IF (cntl%tddft) THEN
       ALLOCATE(rhoo(fpar%nnr1,il_rhoe/fpar%nnr1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(potr(il_rhoe,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       IF (td01%ns_tri.GT.0) nmm=2
    ENDIF
    il_psi=nmm*il_psi
    ALLOCATE(psi(il_psi,1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL give_scr_mddiag(lscr,tag)
    ALLOCATE(scr(lscr),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
99999 IF (cntl%tsampl) THEN
       CALL sample_wait
       IF (cnti%nomore.LT.0) GOTO 10000
    ENDIF
    restf%nfnow=1
    ! ==--------------------------------------------------------------==
    ! TIME STEP FUNCTIONS
    CALL dynit(ekincp,ekin1,ekin2,temp1,temp2,ekinh1,ekinh2)
    ekinp=0.0_real_8    ! McB: cf. META_EXT..()
    ! PARAMETERS FOR THE NOSE-HOOVER THERMOSTATS
    IF (cntl%tnosep.AND.paral%parent) CALL nosepa(1,1)
    ! Dont symmetrize density 
    cntl%tsymrho=.FALSE.
    ! ..Make sure TKFULL=.TRUE
    IF (tkpts%tkpnt.AND.(.NOT.tkpts%tkfull).AND.(.NOT.cntl%tmdeh).AND.paral%io_parent) THEN
       WRITE(6,*)&
            ' ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(6,*)&
            ' WARNING! USE KPOINTS FULL  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(6,*)&
            ' ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == INITIALISATION                                               ==
    ! ==--------------------------------------------------------------==
    IF (cntl%bsymm.AND.paral%io_parent) THEN
       WRITE(6,*)
       WRITE (6,*) 'BSYMM: BS WAVEFUNCTION INITIALIZATION'
    ENDIF
    IF (cntl%cdft)CALL init_cdft()
    CALL initrun(irec,c0(:,:,1:1),c2,sc0,rhoe,psi,eigv)
    IF (cntl%bsymm) THEN
       IF (paral%io_parent) WRITE (6,*) 'BSYMM: HS WAVEFUNCTION INITIALIZATION'
       bsclcs=2
       CALL setbsstate
       CALL initrun(irec,c0(:,:,2:),c2,sc0,rhoe,psi,eigv(1,2))
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Localize wavefunctions if HFX + screening
    IF (hfxc3%twscr.OR.(vdwl%vdwd.AND..NOT.vdwwfl%trwannc))&
         CALL localize2(tau0,c0,c2,sc0,crge%n)
    ! ==--------------------------------------------------------------==
    ! 
    CALL mp_bcast(taup,SIZE(taup),parai%io_source,parai%cp_grp)
    IF (paral%parent) CALL dcopy(3*maxsys%nax*maxsys%nsx,taup,1,taui,1)
    ! ==--------------------------------------------------------------==
    IF (cprint%iprint_step.EQ.0) cprint%iprint_step=cnti%nomore+1
    ! ==--------------------------------------------------------------==
    ! INITIALIZE VELOCITIES
    IF (paral%parent) CALL detdof(tau0,taur)

    IF (irec(irec_vel).EQ.0.AND..NOT.restart1%rgeo) THEN
       ener_com%ecnstr = 0.0_real_8
       ener_com%erestr = 0.0_real_8
       CALL rinvel(velp,c2,nstate)
       IF (paral%parent) CALL taucl(velp)
       IF (paral%parent) CALL rattle(tau0,velp)
       CALL rvscal(velp)
    ELSE
       IF (paral%parent) CALL taucl(velp)
       IF (paral%parent) CALL rattle(tau0,velp)
       IF (cntl%trescale) CALL rvscal(velp)
    ENDIF
    IF (cntl%quenchp) CALL zeroing(velp)!,3*maxsys%nax*maxsys%nsx)
    IF (cntl%trevers) THEN
       ! invert ionic velocities (useful for path sampling)
       CALL dscal(3*maxsys%nax*maxsys%nsx,-1._real_8,velp,1)
    ENDIF
    ! COMPUTE THE IONIC TEMPERATURE TEMPP
    IF (paral%parent) THEN
       CALL ekinpp(ekinp,velp)
       tempp=ekinp*factem*2._real_8/glib
    ENDIF
    call mp_bcast(tempp,parai%io_source,parai%cp_grp)
    ! RESET ACCUMULATORS
    IF (paral%parent)THEN
       IF (irec(irec_ac).EQ.0)CALL resetac(tau0,taui,iteropt%nfi)
       ! 
       itemp=irec(irec_nop1)+irec(irec_nop2)+irec(irec_nop3)&
            +irec(irec_nop4)
       IF (cntl%tnosep .AND. itemp.EQ.0) CALL nospinit(1)
    ENDIF
    ! 
    CALL write_irec(irec)
    ! ..-------------------------------------------------
    ! ..QM/MM coupling by Roethlisberger Group
    ! ..Initailize Molecular Mechanics subsystem:
    ! ..-------------------------------------------------
    ! ..
#if defined (__QMECHCOUPL)
    IF (paral%parent) THEN
       ! ---------------------------------------------------------
       ! CORE QM/MM initilization
       ! Please note this is the same routine used to initialize
       ! the QMMM functionality during a geometry optimization.
       ! cntl%md specific initialization is performed by the routine
       ! 'mm_cpmd_md_init'.
       ! ---------------------------------------------------------
       CALL mm_cpmd_init (cntl%tqmmech,tau0,ions0%na,ions1%nsp,maxsys%nax,&
            ions1%nat,restart1%rco,restart1%rvel,cntr%tempw)
       ! -------------
       ! set QMMM flag
       ! -------------
       IF (cntl%tqmmech) THEN
          CALL mm_run(qmmech)
       ELSE
          qmmech=.FALSE.
       ENDIF
       ! -----------------------------------
       ! initialize QM/MM Molecular Dynamics
       ! -----------------------------------
       IF (qmmech) THEN
          mm_ekin = 0.0_real_8
          mm_temp = 0.0_real_8
          CALL mm_cpmd_md_init
          CALL mm_cpmd_ext_pot_f77
          CALL mm_cpmd_update_links(tau0, ions0%na, ions1%nsp, maxsys%nax, ions1%nat)
       ENDIF
    ENDIF
#endif
    ! INITIALIZE FORCES
    nx = 1
    IF (cntl%tdiag) THEN
       IF (cntl%tlanc) nx=1
       IF (cntl%tdavi) nx=nkpt%ngwk*cnti%ndavv*nkpt%nkpnt+1
       IF (cntl%diis)  nx=((nkpt%ngwk*nstate+8)*cnti%mdiis*nkpt%nkpnt)/4
    ELSEIF (cntl%tsde) THEN
       nx=1
    ELSEIF (cntl%diis) THEN
       nx=(nkpt%ngwk*nstate+8)*cnti%mdiis*nkpt%nkpnt/2+4
    ELSEIF (cntl%pcg) THEN
       nx=1
    ENDIF

    ! INITIALIZE cntl%cdft
    IF (cntl%cdft)THEN
       CALL cdft_w(rhoe,tau0,dummy)
       IF (cntl%cdft_weight)CALL write_w(wdiff,"INIT")
    ENDIF

    IF (paral%io_parent) THEN
       WRITE(6,'(1X,64("="))')
       WRITE(6,'(1X,"==",T25,A,T64,"==")')&
            'FORCES INITIALIZATION'
       WRITE(6,'(1X,64("="))')
    ENDIF

    ! initialize nabdy
    CALL nabdy_mem
    !call nabdy_spline_mem
    CALL nabdy_load(tempp)
    CALL nabdy_initialize_action
    CALL nabdy_normalize(2)
    CALL dcopy(3*maxsys%nax*maxsys%nsx*nabdyvar%ntrajbd,nabdyvec%nacoor(1,1,1,1),1, &
         nabdyvec%nacoorp(1,1,1,1),1)
    CALL integrate_atom(1,1,100)

    DO itraj=1,nabdyvar%ntrajbd
       CALL dcopy(3*maxsys%nax*maxsys%nsx,nabdyvec%nacoor(1,1,1,itraj),1,tau0(1,1,1),1)
       CALL dcopy(3*maxsys%nax*maxsys%nsx,nabdyvec%namom(1,1,1,itraj),1,velp(1,1,1),1)
       CALL mom_to_vel(velp)
       CALL dcopy(3*maxsys%nax*maxsys%nsx,velp(1,1,1),1,velpmem(1,1,1,itraj),1)
       ifcalc=0
       CALL phfac(tau0)
       IF (paral%io_parent) WRITE(6,'(/,a,i20,/)') &
            '  loop over fluid elements (initialization): ', itraj
       IF(paral%parent) CALL wrgeo(tau0)
       IF(cntl%bsymm) THEN
          CALL bs_forces_diag(nstate,c0,c2,cm,sc0,cm(nx:),vpp,eigv,&
               rhoe,psi,&
               tau0,velp,taui,fion,ifcalc,&
               irec,.TRUE.,.TRUE.)
       ELSE
          CALL forces_diag(nstate,c0,c2,cm,sc0,cm(nx:),vpp,eigv,&
               rhoe,psi,&
               tau0,velp,taui,fion,ifcalc,&
               irec,.TRUE.,.TRUE.)
       ENDIF
       IF (cntl%tddft) THEN
          CALL lr_tddft(c0(:,:,1),c1,c2,sc0,rhoe,psi,tau0,fion,eigv,&
               nstate,.TRUE.,td01%ioutput)
       ELSE
          CALL dcopy(2*nkpt%ngwk*nstate,c0(1,1,1),1,c0mem(1,1,itraj),1)
          CALL dcopy(2*nkpt%ngwk*nstate,c2(1,1),1,c2mem(1,1,itraj),1)
          CALL dcopy(3*maxsys%nax*maxsys%nsx,fion(1,1,1),1,fionmem(1,1,1,itraj),1)
          CALL dcopy(fpar%nnr1*clsd%nlsd,rhoe(1,1),1,rhoemem(1,1,itraj),1)
       ENDIF
    ENDDO
    CALL nabdy_keep_classical(velpmem)  

    ! switch on info printing for the lin resp part
    IF (cntl%tresponse) dmbi%inter_pt_firstcall = .TRUE.

    IF (.NOT.cntl%bsymm) CALL dcopy(nnx,rin0,1,rm1(1),1)
    ! INITIALIZE THERMOSTATS
    IF (paral%io_parent) THEN
       CALL wrgeof(tau0,fion)
       ! FILEN=FPATH(IAPATH:IEPATH)//'ENERGIES'
       filen='ENERGIES'
       IF (cntl%tpmin) THEN
          ftmp='ENERGIES_'
          CALL mw_filename(ftmp,filen,ipcurr)
          CALL xstring(filen,i1,i2)
          ftmp=filen(i1:i2)//'.'
          CALL mw_filename(ftmp,filen,mfepi%istring)
       ENDIF
       CALL fileopen(3,filen,fo_app+fo_verb+fo_nofp,ferror)
    ENDIF
    IF (paral%io_parent) THEN
       WRITE(6,'(1X,64("="))')
       WRITE(6,'(1X,"==",T20,A,T64,"==")')&
            'END OF FORCES INITIALIZATION'
       WRITE(6,'(1X,64("="),/)')
    ENDIF

    ! INITIALIZE GLE THERMO
    CALL gle_init(tau0,velp,rmass%pma)

    CALL write_irec(irec)
    ! ==--------------------------------------------------------------==
    ! == END INITIALIZATION                                           ==
    ! ==--------------------------------------------------------------==
    IF (teststore(0).AND.cntl%tsampl)&
         CALL zhwwf(2,irec,c0,c2,nstate,eigv,tau0,velp,taui,iteropt%nfi)
    IF (paral%io_parent) THEN
       ! MEAN SQUARE DISPLACEMENT OF DIFFERENT IONIC SPECIES
       IF (paral%parent) CALL dispp(tau0,taui,disa)
       ! ENERGY OF THE NOSE THERMOSTATS
       CALL noseng(iteropt%nfi,velp,enose,enosp,dummy(1),1) 
       econs=ekinp+ener_com%etot+enose+enosp+ener_com%ecnstr
       IF (cntl%cdft)econs=econs+cdftcom%cdft_v(1)*(cdftcom%cdft_nc+cdftcom%vgrad(1))
       time2 =m_walltime()
       tcpu = (time2 - time1)*0.001_real_8
       WRITE(6,'(A,T50,F8.2,A8)') ' TIME FOR INITIALIZATION:',&
            tcpu,' SECONDS'
       WRITE(6,'(//,1X,64("="))')
       WRITE(6,'(1X,"=",T20,A,T65,"=")')&
            'MOLECULAR DYNAMICS SIMULATION'
       WRITE(6,'(1X,64("="))')
       ! CALL WRPRINT_MD(EIGV,F,AMU,NSTATE,TAU0,FION,
       ! &                  0._real_8,TEMPP,ETOT,ECONS,0._real_8,DISA,
       ! &                  TCPU,.FALSE.,NFI,0)
    ENDIF
    ! ==================================================================
    ! ==          THE BASIC LOOP FOR MOLECULAR DYNAMICS               ==
    ! ==                 USING VELOCITY VERLET                        ==
    ! ==================================================================
#if defined(USE_IBM_HPM)
    CALL hpm_start('MD LOOP')
#endif
    infi=0
    nfimin=iteropt%nfi+1
    nfimax=iteropt%nfi+cnti%nomore
    DO loopnfi=nfimin,nfimax

       ! FIRST UPDATE OF THE AMPLITUDES: dT/2 USING THE POSITIONS X(t) [TAU0=NACOOR]
       CALL nabdy_ampli_update(cntr%delt_ions/2.0,1)
       CALL nabdy_normalize(1)
       CALL mp_sync(parai%cp_grp)
       !GET THE AMPLITUDE FOR THE EVALUATION OF THE RHS OF THE ODE FOR THE AMPLITUDE
       !-----------------------------------------------------------------
       DO itraj=1,nabdyvar%ntrajbd
          !-----------------------------------------------------------------
          IF (paral%io_parent) THEN
             WRITE(6,'(/,a,i37,/)') '  loop over fluid elements: ', itraj
          ENDIF
          !IF (LOOPNFI.GT.NFIMIN) THEN
          CALL dcopy(2*nkpt%ngwk*nstate,c0mem(1,1,itraj),1,c0(1,1,1),1)
          CALL dcopy(2*nkpt%ngwk*nstate,c2mem(1,1,itraj),1,c2(1,1),1)
          !ENDIF
          CALL dcopy(3*maxsys%nax*maxsys%nsx,nabdyvec%nacoor(1,1,1,itraj),1,tau0(1,1,1),1)
          CALL dcopy(3*maxsys%nax*maxsys%nsx,velpmem(1,1,1,itraj),1,velp(1,1,1),1)
          CALL dcopy(3*maxsys%nax*maxsys%nsx,fionmem(1,1,1,itraj),1,fion(1,1,1),1)
          CALL dcopy(fpar%nnr1*clsd%nlsd,rhoemem(1,1,itraj),1,rhoe(1,1),1)
          IF(paral%io_parent) CALL wrgeo(tau0)

          time1=m_walltime()
          CALL mp_sync(parai%cp_grp)
          infi=infi+1
          iteropt%nfi=loopnfi
          comvl%subcom=comvl%tsubcom.AND.MOD(iteropt%nfi-1,comvl%ncomv).EQ.0
          comvl%subrot=comvl%tsubrot.AND.MOD(iteropt%nfi-1,comvl%nrotv).EQ.0
          tstrng=cntl%tpmin.AND.MOD(iteropt%nfi,cnti%nomore).EQ.0
          ! ANNEALING
          CALL anneal(velp,c2,nstate,scr)
          CALL berendsen(velp,c2,nstate,scr,0.0_real_8,0.0_real_8)
          ! UPDATE NOSE THERMOSTATS
          CALL noseup(velp,c2,nstate,1)
          ! FIRST HALF OF GLE EVOLUTION
          CALL gle_step(tau0,velp,rmass%pma)

          ! SUBTRACT CENTER OF MASS VELOCITY
          IF (paral%io_parent.AND.comvl%subcom) CALL comvel(velp,vcmio,.TRUE.)
          ! SUBTRACT ROTATION AROUND CENTER OF MASS
          IF (paral%io_parent.AND.comvl%subrot) CALL rotvel(tau0,velp,lmio,&
               tauio,.TRUE.)
          ! UPDATE VELOCITIES
          IF (paral%io_parent) CALL velupi(velp,fion,1)
#if defined (__QMECHCOUPL)
          IF (paral%io_parent .AND. qmmech) THEN
             CALL mm_cpmd_velup(cntr%delt_ions)
          ENDIF
#endif
          ! UPDATE POSITIONS
          ener_com%ecnstr = 0.0_real_8
          IF (paral%io_parent) THEN
             CALL posupi(tau0,taup,velp)
             IF (cotc0%mcnstr.NE.0) CALL cpmdshake(tau0,taup,velp)
#if defined (__QMECHCOUPL)
             IF (qmmech) THEN
                CALL mm_cpmd_update_links(taup, ions0%na, ions1%nsp, maxsys%nax, ions1%nat)
                CALL mm_cpmd_posup(cntr%delt_ions)
             ENDIF
#endif
          ENDIF
          CALL mp_bcast(taup,SIZE(taup),parai%io_source,parai%cp_grp)
          ! Reset swap file to update (option BLOCK and ALL)
          IF (tkpts%tkblock) CALL reskpt_swap
          CALL phfac(taup)
          IF (corel%tinlc) CALL copot(rhoe,psi,ropt_mod%calste)
          IF (tmovr) THEN
             CALL dcopy(nnx,rin0,1,rhoe(1,1),1)
             CALL moverho(rhoe,psi)
             CALL dcopy(nnx,rhoe(1,1),1,rin0,1)
          ENDIF
          IF (fint1%ttrot) THEN
             CALL calc_alm
          ENDIF
          ! CALCULATE THE FORCES
          ropt_mod%calste=cntl%tpres.AND.MOD(iteropt%nfi,cnti%npres).EQ.0
          IF (cntl%textrap) THEN
             ! Extrapolate wavefunctions
             CALL extrapwf(infi,c0,scr,cold,nnow,numcold,nstate,cnti%mextra)
          ENDIF
          ! mb - Wannier stuff for vdW-WC
          vdwwfl%twannup=vdwwfl%twannup.OR.(infi.EQ.1.AND..NOT.vdwwfl%trwannc)
          IF (vdwl%vdwd.AND.vdwwfl%twannup) THEN
             CALL localize2(taup,c0,c2,sc0,nstate)
          ENDIF
          IF (cntl%tlanc) nx=1
          IF (cntl%tdavi) nx=cnti%ndavv*nkpt%nkpnt+1
          ! RESPONSE calculation
          IF (cntl%tresponse) THEN
             ! save the cntl%md step number
             save_nfi=iteropt%nfi
             ! localisation at first iteration + PT
             CALL mddiag_interaction_p(nstate,c0,c2,cm,sc0,&
                  CM(NX:),VPP,EIGV,RHOE,PSI,&! SCR,LSCR,& !vw doesnt match procedure args
                  TAUP,VELP,TAUI,FION,IFCALC,IREC,.TRUE.,.FALSE.)
             ! sets the number of iteration and recover cntl%md step number
             ifcalc=iteropt%nfi
             iteropt%nfi=save_nfi

             ! switch off the info printing
             dmbi%inter_pt_firstcall=.FALSE.
          ELSE
             IF (cntl%bsymm) THEN
                CALL bs_forces_diag(nstate,c0,c2,cm,sc0,cm(nx:),vpp,eigv,&
                     rhoe,psi,&
                     taup,velp,taui,fion,ifcalc,&
                     irec,.TRUE.,.FALSE.)
             ELSE
                CALL forces_diag(nstate,c0,c2,cm,sc0,cm(nx:),vpp,eigv,&
                     rhoe,psi,&
                     taup,velp,taui,fion,ifcalc,&
                     irec,.TRUE.,.FALSE.)
             ENDIF
             IF (cntl%tddft) THEN
                CALL lr_tddft(c0(:,:,1),c1,c2,sc0,rhoe,psi,taup,fion,eigv,&
                     nstate,.TRUE.,td01%ioutput)
             ENDIF
          ENDIF
          ! PATH MINIMIZATION
          ! ==================================================================
          ! Damped Dynamics
          CALL dampdyn(velp,fion,cm,c2,nstate,scr(1),scr(10))
          ! ==================================================================
          IF(paral%parent) WRITE(6,*) '... adding quantum forces ...'
          CALL nabdy_forces_calc(itraj)
          ! ==================================================================

          IF (ropt_mod%calste) CALL totstr
#if defined (__QMECHCOUPL)
          ! ----------------------------------------------
          ! QMMM electrostatic coupling
          ! ----------------------------------------------
          IF (qmmech) THEN
             CALL mm_cpmd_ext_pot_f77
             CALL mm_cpmd_elstat(taup,ions0%na,ions1%nsp,maxsys%nax,ions1%nat,c0,scr)
          ENDIF
#endif
          ! FINAL UPDATE FOR VELOCITIES
          ener_com%ecnstr = 0.0_real_8
          IF (paral%io_parent) THEN
             CALL velupi(velp,fion,1)
#if defined (__QMECHCOUPL)
             IF (qmmech) THEN
                CALL mm_cpmd_velup(cntr%delt_ions)
             ENDIF
#endif
             IF (isos1%twall) CALL box_boundary(taup,velp)
             CALL rattle(taup,velp)
          ENDIF
          IF (paral%io_parent) CALL geofile(taup,velp,'WRITE')
          ! COMPUTE THE IONIC TEMPERATURE TEMPP
          IF (paral%io_parent) THEN
             CALL ekinpp(ekinp,velp)
             tempp=ekinp*factem*2._real_8/glib
#if defined (__QMECHCOUPL)
             ! -------------------------------------------------
             ! TEMPERTURE and KINETIC ENERGY of the MM subsystem
             ! MM_EKIN + EKINP = KE of whole QM/MM system
             ! -------------------------------------------------
             IF (qmmech) THEN
                CALL mm_cpmd_ekin(mm_ekin,mm_temp)
             ENDIF
#endif
          ENDIF
#if defined (__QMECHCOUPL)
          ! ---------------------------------------------------------------
          ! MM Temperature control - annealing schedule not yet implemented
          ! ---------------------------------------------------------------
          IF (paral%io_parent .AND. qmmech) THEN
             CALL mm_cpmd_temp_control(temp1,temp2,mm_temp,cntl%tcp)
             CALL mm_cpmd_ekin(mm_ekin,mm_temp)
          ENDIF
#endif
          ! IONIC TEMPERATURE CONTROL
          IF (paral%io_parent) CALL rscvp(temp1,temp2,tempp,velp)
          ! SUBTRACT ROTATION AROUND CENTER OF MASS
          IF (paral%io_parent.AND.comvl%subrot) CALL rotvel(tau0,velp,lmio,&
               tauio,.FALSE.)
          ! SUBTRACT CENTER OF MASS VELOCITY
          IF (paral%io_parent.AND.comvl%subcom) CALL comvel(velp,vcmio,.FALSE.)

          ! SECOND HALF OF GLE EVOLUTION
          CALL gle_step(tau0,velp,rmass%pma)

          ! UPDATE NOSE THERMOSTATS
          CALL noseup(velp,c2,nstate,1)
          CALL berendsen(velp,c2,nstate,scr,0.0_real_8,0.0_real_8)
          ! ANNEALING
          CALL anneal(velp,c2,nstate,scr)
          IF (paral%io_parent) THEN
             CALL ekinpp(ekinp,velp)
             tempp=ekinp*factem*2._real_8/glib
          ENDIF
          ! MEAN SQUARE DISPLACEMENT OF DIFFERENT IONIC SPECIES
          IF (paral%io_parent) CALL dispp(taup,taui,disa)
          ! ENERGY OF THE NOSE THERMOSTATS
          IF (paral%io_parent) CALL noseng(iteropt%nfi,velp,enose,enosp,dummy(1),1) 
          ! CALCULATE PROPERTIES DURING SIMULATION.
          cntl%caldip=cntl%tdipd.AND.MOD(iteropt%nfi-1,cnti%npdip).EQ.0
          IF (.NOT.cntl%bsymm)&
               CALL propcal(c0,c2,cm,sc0,taup,eigv,crge%f,ener_com%amu,&
               rhoe,psi,nstate,nkpt%nkpnt,iteropt%nfi,infi)
          ! PRINTOUT the evolution of the accumulators every time step
          IF (paral%io_parent) THEN
             nabdyvec%natotqp(itraj)=0._real_8
             DO is=1,ions1%nsp
                IF (ions0%iatyp(is).GT.nabdyvar%nabdy_zmax) CYCLE
                DO ia=1,ions0%na(is)
                   nabdyvec%natotqp(itraj)=nabdyvec%natotqp(itraj)+nabdyvec%naqp_pat(ia,is)
                ENDDO
             ENDDO
             nabdyvec%natotk(itraj)=ekinp
             nabdyvec%natotcl(itraj)=ener_com%etot
          ENDIF
          IF (paral%io_parent) THEN
             econs=ekinp+ener_com%etot+enose+enosp+ener_com%ecnstr+ekincv+vharm+glepar%egle
             IF (cntl%cdft)THEN
                econs=econs+cdftcom%cdft_v(1)*(cdftcom%cdft_nc+cdftcom%vgrad(1))
             ENDIF
#if defined (__QMECHCOUPL)
             IF (qmmech) THEN
                econs=econs + mm_ekin
                IF (paral%io_parent)&
                     WRITE (6,'(50x, "TEMP: ",f10.0)') mm_temp
             ENDIF
#endif
             time2=m_walltime()
             tcpu=(time2-time1)*0.001_real_8
             CALL wrprint_md(eigv,crge%f,ener_com%amu,nstate,taup,fion,&
                  0._real_8,tempp,ener_com%etot,econs,0._real_8,disa,&
                  tcpu,.FALSE.,iteropt%nfi,infi)
             ! UPDATE ACCUMULATORS
             CALL paccc(tempp,ener_com%etot,econs,enose,enosp,ener_com%ecnstr,ener_com%erestr,&
                  ener_com%ebogo,disa,tcpu,iteropt%nfi,1)
             ! STORE IONIC COORDINATES AND VELOCITIES FOR STATISTICS
             ropt_mod%movie=rout1%mout .AND. MOD(iteropt%nfi-1,cnti%imovie).EQ.0
             ropt_mod%rprint=rout1%rout .AND. MOD(iteropt%nfi-1,cnti%ntraj).EQ.0
             ropt_mod%txyz=rout1%xtout .AND. MOD(iteropt%nfi-1,cnti%ntraj).EQ.0
             ropt_mod%tdcd=rout1%dcout .AND. MOD(iteropt%nfi-1,cnti%ntraj).EQ.0
             CALL printp(taur,taup,velp)
             IF (cprint%twriteforcetrajectory) CALL printp2(taur,taup,velp,fion)
#if defined (__QMECHCOUPL)
             ! --------------------------------------------------------
             ! ..       Write to QMMM trajectory.  This trajectory file contains
             ! the atoms that make up the "real" system. (no dummies)
             ! --------------------------------------------------------
             IF (qmmech) THEN
                IF (rout1%rout .AND. MOD(iteropt%nfi-1,cnti%ntraj).EQ.0) THEN
                   CALL mm_cpmd_write_trajectory(iteropt%nfi)
                ENDIF
             ENDIF
#endif
             time2=m_walltime()
             tcpu=(time2-time1)*0.001_real_8
          ENDIF
          IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)
          IF (tstrng) THEN
             soft_com%exsoft=.TRUE.
             soft_com%exnomore=.TRUE.
          ELSE IF (iteropt%nfi.EQ.nfimax) THEN
             soft_com%exsoft=.TRUE.
          ENDIF
          ! periodic output of density/wavefunction etc.
          IF (rout1%rhoout.AND.(rout1%nrhoout.GT.0)) THEN
             IF (MOD(iteropt%nfi-1,rout1%nrhoout).EQ.0) THEN
                CALL rhopri(c0,tau0,rhoe,psi(:,1),nstate,nkpt%nkpnt)
             ENDIF
          ENDIF
          IF (teststore(iteropt%nfi).OR.soft_com%exsoft) THEN
             CALL zhwwf(2,irec,c0,c2,nstate,eigv,taup,velp,taui,iteropt%nfi)
          ENDIF

#if defined (__QMECHCOUPL)
          ! -----------------------
          ! Write QMMM restart file
          ! -----------------------
          IF (paral%io_parent.AND.qmmech) THEN
             IF (MOD(iteropt%nfi,store1%istore).EQ.0.OR.iteropt%nfi.EQ.nfimax.OR.soft_com%exsoft) THEN
                CALL mm_cpmd_write_restart
             ENDIF
          ENDIF
#endif
          ! temperature ramping
          CALL tempramp(temp1,temp2)
          ! UPDATE IONIC POSITIONS
          CALL dcopy(3*maxsys%nax*maxsys%nsx,taup(1,1,1),1,tau0(1,1,1),1)
          ! UPDATE DENSITY
          IF (.NOT.cntl%bsymm) THEN
             CALL extrap(nnx,andr2%alxmix,rm1,rin0,rinp)
             CALL dcopy(nnx,rin0,1,rm1(1),1)
             CALL dcopy(nnx,rinp(1),1,rin0,1)
          ENDIF
          ! NABDY[
          CALL dcopy(2*nkpt%ngwk*crge%n,c0(1,1,1),1,c0mem(1,1,itraj),1)
          CALL dcopy(2*nkpt%ngwk*crge%n,c2(1,1),1,c2mem(1,1,itraj),1)
          CALL dcopy(fpar%nnr1*clsd%nlsd,rhoe(1,1),1,rhoemem(1,1,itraj),1)
          CALL dcopy(3*maxsys%nax*maxsys%nsx,velp(1,1,1),1,velpmem(1,1,1,itraj),1)
          CALL dcopy(3*maxsys%nax*maxsys%nsx,fion(1,1,1),1,fionmem(1,1,1,itraj),1)
          CALL dcopy(3*maxsys%nax*maxsys%nsx,tau0(1,1,1),1,nabdyvec%nacoor(1,1,1,itraj),1)
          CALL dcopy(3*maxsys%nax*maxsys%nsx,velp(1,1,1),1,nabdyvec%namom(1,1,1,itraj),1)
          CALL vel_to_mom(nabdyvec%namom,itraj)

          !NABDY UPDATE OF THE ACTION, D ACTION, DD ACTION
          CALL nabdy_action(rhoe,velpmem,itraj,loopnfi)
          !=================================================================
          !.==     END OF FLUID ELEMENTS LOOP                              ==
          !=================================================================
       ENDDO ! loop over itraj

       ! NABDY]     
       ! FINAL ENERGIES
       IF (paral%io_parent) THEN
          !total amplitude per fluid element
          DO itraj=1,nabdyvar%ntrajbd
             nabdyvec%naampfe(itraj)=1._real_8
             DO is=1,ions1%nsp
                IF (ions0%iatyp(is).GT.nabdyvar%nabdy_zmax) CYCLE
                DO ia=1,ions0%na(is)
                   IF (nabdyvec%naampl(ia,is,itraj) .NE. 0.0_real_8) then
                     nabdyvec%naampfe(itraj)=nabdyvec%naampfe(itraj)*nabdyvec%naampl(ia,is,itraj)
                   ENDIF                   
                ENDDO
             ENDDO
          ENDDO

          totamplfe=0._real_8
          DO itraj=1,nabdyvar%ntrajbd
             totamplfe=totamplfe+nabdyvec%naampfe(itraj)**2
          ENDDO
          DO itraj=1,nabdyvar%ntrajbd
             nabdyvec%naampfe(itraj)=nabdyvec%naampfe(itraj)/SQRT(totamplfe)
          ENDDO
          !total enery per fluid element
          DO k=1,3
             natotener(k)=0._real_8
          ENDDO
          DO itraj=1,nabdyvar%ntrajbd
             natotener(1)=natotener(1)+nabdyvec%natotk(itraj)*nabdyvec%naampfe(itraj)**2
             natotener(2)=natotener(2)+nabdyvec%natotcl(itraj)*nabdyvec%naampfe(itraj)**2
             natotener(3)=natotener(3)+nabdyvec%natotqp(itraj)*nabdyvec%naampfe(itraj)**2
          ENDDO
          natotener(4)=natotener(1)+natotener(2)+natotener(3)

          fexist=.FALSE.
          INQUIRE(file='nabdy_energies',exist=fexist)
          IF (.NOT.fexist) THEN
             OPEN(unit=2111,file='nabdy_energies',status='new')
             OPEN(unit=2112,file='nabdy_fe_ampli',status='new')
          ELSE
             OPEN(unit=2111,file='nabdy_energies',status='unknown',position='append')
             OPEN(unit=2112,file='nabdy_fe_ampli',status='unknown',position='append')
          ENDIF
          WRITE(2111,'(i10,2x,4f15.5)') loopnfi,(natotener(k),k=1,4)
          CALL m_flush(2111)
          WRITE(2112,'(i10,2x,100f15.5)') loopnfi,(nabdyvec%naampfe(itraj),itraj=1,nabdyvar%ntrajbd)
          CALL m_flush(2112)
          CLOSE(2111)
          CLOSE(2112)

          !energy conservation 
          IF (ifirst.EQ.0) THEN
             natotene_ref=natotener(1)
             ifirst=ifirst+1
          ENDIF
          natotene_diff=natotene_ref-natotener(1)
          DO itraj=1,nabdyvar%ntrajbd
             ekprime=nabdyvec%natotk(itraj)*nabdyvec%naampfe(itraj)**2
             etotprime=nabdyvec%natotk(itraj)* nabdyvec%naampfe(itraj)**2 + &
                  nabdyvec%natotcl(itraj)*nabdyvec%naampfe(itraj)**2 + &
                  nabdyvec%natotqp(itraj)*nabdyvec%naampfe(itraj)**2
             !deltaeava=(natotene_diff/nabdyvar%ntrajbd)
             !deltaeava   = natotene_diff 
             !deltaeprime = deltaeava ! for the moment deltae(itraj) =< de_tot>
             !call ekinpp(naekinp,velpmem(1,1,1,itraj))
             !SQRT((naekinp-natotene_diff)/naekinp)
             !alfap=SQRT(MAX((deltaeprime/ekprime),-0.9_real_8)+1.0_real_8)
             alfap=SQRT((natotener(1)+natotene_diff)/natotener(1))
             WRITE(6,'(a,3f15.5)') 'alphap ',alfap
             IF (nabdyvar%scalep) THEN
              DO is=1,ions1%nsp
                IF (ions0%iatyp(is).GT.nabdyvar%nabdy_zmax) CYCLE
                DO ia=1,ions0%na(is)
                   velpmem(1,ia,is,itraj)=alfap*velpmem(1,ia,is,itraj)
                   velpmem(2,ia,is,itraj)=alfap*velpmem(2,ia,is,itraj)
                   velpmem(3,ia,is,itraj)=alfap*velpmem(3,ia,is,itraj)
                ENDDO
              ENDDO
             ENDIF
          ENDDO
          WRITE(6,'(a,3f15.5)') 'e_kin^ref,e_kin^now',natotene_ref,natotener(1)
       ENDIF
       CALL mp_bcast(velpmem,SIZE(velpmem),parai%source,parai%allgrp)
       CALL dcopy(3*maxsys%nax*maxsys%nsx*nabdyvar%ntrajbd,velpmem,1,nabdyvec%namom,1)
       DO itraj=1,nabdyvar%ntrajbd
          CALL vel_to_mom(nabdyvec%namom,itraj)
       ENDDO
       !FINAL UPDATE OF THE AMPLITUDES: dT 
       !COMPUTE THE AVERAGE POSITION (X(t)+X(t+dt))/2 [(TAU0+TAUP)/2]
       !AND USE THE AMPLITUDE COMPUTED FROM THE FIST UPDATE 
       CALL nabdy_ampli_update(cntr%delt_ions,2)

       IF (nabdyvar%tnafric) THEN
          CALL nabdy_friction_cm(velpmem,nstate)
       ENDIF

       IF (nabdyvar%tnarepart) THEN
         CALL check_small_ampli
         !CALL nabdy_repartition_ampli
       ENDIF

       CALL nabdy_keep_classical(velpmem)

       CALL nabdy_normalize(2)
       CALL mp_sync(parai%cp_grp)
       nabdyvar%nabdy_time=nabdyvar%nabdy_time+cntr%delt_ions
       CALL write_nabdy_traj

       IF (soft_com%exsoft) GOTO 100
       IF (cntl%cdft)CALL cdft_reinit()
       ! ==================================================================
       ! ==     END OF MAIN LOOP                                         ==
       ! ==================================================================
    ENDDO
    IF (paral%io_parent) THEN
       WRITE(6,'(1X,64("="))')
       WRITE(6,'(1X,"=",T17,A,T65,"=")')&
            'END OF MOLECULAR DYNAMICS SIMULATION'
       WRITE(6,'(1X,64("="),/,/)')
    ENDIF
    ! ==--------------------------------------------------------------==
100 CONTINUE

#if defined(USE_IBM_HPM)
    CALL hpm_stop('MD LOOP')
#endif

    IF (cntl%cdft)THEN
       IF (cntl%cdft_weight)CALL write_w(wdiff,"FINAL")
    ENDIF

    IF (rout1%rhoout.AND.(rout1%nrhoout.LE.0)) THEN
       CALL rhopri(c0,tau0,rhoe,psi(:,1),nstate,nkpt%nkpnt)
       IF (cntl%bsymm) THEN
          bsclcs=2
          CALL rhopri(c0(:,:,2:),tau0,rhoe,psi(:,1),nstate,nkpt%nkpnt)
          bsclcs=1
       ENDIF
    ENDIF
    ! Print accumulators.
    IF (paral%io_parent) CALL paccc(tempp,ener_com%etot,econs,enose,enosp,&
         ener_com%ecnstr,ener_com%erestr,ener_com%ebogo,disa,tcpu,iteropt%nfi,0)
    IF (paral%io_parent) CALL gsize(fion,gnmax,gnorm)
    IF (paral%io_parent) THEN
       IF (.NOT.cntl%bsymm) THEN
          CALL finalp(tau0,fion,velp,eigv)
       ELSE
          CALL wrgeof(tau0,fion)
          WRITE(6,'(A)') ' NUCLEAR GRADIENT:'
          WRITE(6,'(2(A,1PE15.5))') '    MAX. COMPONENT =',&
               gnmax,'         NORM =',gnorM
       ENDIF
    ENDIF
    IF (cntl%tsampl) THEN
       CALL sample_go
       GOTO 99999
    ENDIF
10000 CONTINUE
    ! ==--------------------------------------------------------------==
    DEALLOCATE(eigv,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(taur,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(taui,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(fion,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(taup,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (paral%io_parent) CALL fileclose(3)
    DEALLOCATE(rhoe,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(psi,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(scr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (cntl%textrap) DEALLOCATE(cold,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (comvl%tsubrot) DEALLOCATE(tauio,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (.NOT.cntl%bsymm) THEN
       DEALLOCATE(rinp,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(rm1,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(rmix,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(rout0,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(rin0,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    DEALLOCATE(c0mem,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(c2mem,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rhoemem,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(fionmem,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(velpmem,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__) 
    CALL nabdy_rel
    !call nebdy_spline_free
    !
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE nabdy_dyn

  ! ==================================================================
  SUBROUTINE nabdy_getnorm(icall)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: icall

    INTEGER                                  :: i, ia, is
    REAL(real_8)                             :: ampli

!
! The 3D gaussians assiciated to each fluid element has the form:
! phi(x)=a (1/((2pi)^3/2 * s^3) exp(-((x-xo)^2)/(2 s^2)).
! The norm of phi(x) in 3D is:
! N = a 
! For the phi^2(x)
! N = a^2 / (8 pi^(3/2) s^3)
!
! select one atom

    IF (paral%io_parent) THEN
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             ampli=0._real_8
             !    run over the fluuid elements
             IF(icall.EQ.1) THEN
                DO i=1,nabdyvar%ntrajbd
                   ampli=ampli+nabdyvec%naampl_temp(ia,is,i)**2 / &
                        (8._real_8 * pi**(1.5) * nabdyvec%naomega(ia,is,i)**(3.0))
                ENDDO
             ELSE
                DO i=1,nabdyvar%ntrajbd
                   ampli=ampli+nabdyvec%naampl(ia,is,i)**2 / &
                        (8._real_8 * pi**(1.5) * nabdyvec%naomega(ia,is,i)**(3.0))
                ENDDO
             ENDIF
             WRITE(6,'(A,2I6,A,F12.5)') 'Norm of atom ',is,ia,': ',ampli
          ENDDO
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE nabdy_getnorm
  ! ==================================================================
  SUBROUTINE nabdy_normalize(icall)
    ! ==--------------------------------------------------------------==
    !
    INTEGER                                  :: icall

    INTEGER                                  :: i, ia, is, j
    REAL(real_8)                             :: ampli, distij, ni, nj, nsinv, &
                                                sis

!
! The 3D gaussians assiciated to each fluid element has the form:
! phi(x)=a (1/((2pi)^3/2 * s^3) exp(-((x-xo)^2)/(2 s^2)).
! The norm of phi(x) in 3D is:
! N = ampl 
! For the phi^2(x)
! N = a^2 / (8 pi^(3/2) s^3)
!

    IF (paral%io_parent) THEN
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             ampli=0._real_8
             !      run over the fluid elements
             IF(icall.EQ.1) THEN
                DO i=1,nabdyvar%ntrajbd
                   ampli=ampli+nabdyvec%naampl_temp(ia,is,i)**2  / &
                        (8._real_8 * pi**(1.5) * nabdyvec%naomega(ia,is,i)**(3.0))
                   ni=1._real_8/((2._real_8*pi)**(1.5) * (nabdyvec%naomega(ia,is,i)**3))
                   sis=(1._real_8/SQRT(2._real_8))*nabdyvec%naomega(ia,is,i)
                   DO j=i+1,nabdyvar%ntrajbd
                      nj=ni
                      nsinv=(2._real_8*pi)**(1.5) * sis**3
                      distij=SQRT((nabdyvec%nacoor(1,ia,is,i)-nabdyvec%nacoor(1,ia,is,j))**2 + &
                           (nabdyvec%nacoor(2,ia,is,i)-nabdyvec%nacoor(2,ia,is,j))**2 + &
                           (nabdyvec%nacoor(3,ia,is,i)-nabdyvec%nacoor(3,ia,is,j))**2 &
                           )
                      ampli=ampli+ 2._real_8 * &
                           (nabdyvec%naampl_temp(ia,is,i)*nabdyvec%naampl_temp(ia,is,j)*ni*nj) * &
                           nsinv * &
                           EXP(-((1._real_8/SQRT(2._real_8))* distij**2)/(2._real_8*sis**2)  )
                   ENDDO
                ENDDO
                DO i=1,nabdyvar%ntrajbd
                   nabdyvec%naampl_temp(ia,is,i)=nabdyvec%naampl_temp(ia,is,i)/SQRT(ampli)
                ENDDO
             ELSE
                DO i=1,nabdyvar%ntrajbd
                   ampli=ampli+nabdyvec%naampl(ia,is,i)**2 / &
                        (8._real_8 * pi**(1.5) * nabdyvec%naomega(ia,is,i)**(3.0))
                   ni=1._real_8/((2._real_8*pi)**(1.5) * (nabdyvec%naomega(ia,is,i)**3))
                   sis=(1._real_8/SQRT(2._real_8))*nabdyvec%naomega(ia,is,i)
                   DO j=i+1,nabdyvar%ntrajbd
                      nj=ni
                      nsinv=(2._real_8*pi)**(1.5) * sis**3
                      distij=SQRT((nabdyvec%nacoor(1,ia,is,i)-nabdyvec%nacoor(1,ia,is,j))**2 + &
                           (nabdyvec%nacoor(2,ia,is,i)-nabdyvec%nacoor(2,ia,is,j))**2 + &
                           (nabdyvec%nacoor(3,ia,is,i)-nabdyvec%nacoor(3,ia,is,j))**2 &
                           )
                      ampli=ampli+ 2._real_8 * &
                           (nabdyvec%naampl(ia,is,i)*nabdyvec%naampl(ia,is,j)*ni*nj) * &
                           nsinv * &
                           EXP(-((1._real_8/SQRT(2._real_8))* distij**2)/(2._real_8*sis**2)  )
                   ENDDO
                ENDDO
                DO i=1,nabdyvar%ntrajbd
                   nabdyvec%naampl(ia,is,i)=nabdyvec%naampl(ia,is,i)/SQRT(ampli)
                ENDDO
             ENDIF
             WRITE(6,*) 'ampli in nabdy_normalize',ampli,icall
          ENDDO
       ENDDO
    ENDIF
    ! msglen = maxsys%nax*maxsys%nsx*nabdyvar%ntrajbd * 8
    ! write(*,*) 'msglen',msglen
    IF (icall.EQ.1) THEN
       CALL mp_bcast(nabdyvec%naampl_temp,SIZE(nabdyvec%naampl_temp),parai%source,parai%allgrp)
    ELSE
       CALL mp_bcast(nabdyvec%naampl,SIZE(nabdyvec%naampl_temp),parai%source,parai%allgrp)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE nabdy_normalize
  ! ==================================================================
  ! subroutine mom_to_vel(vect)
  ! Arguments
  ! real(real_8)    :: vect(3,maxsys%nax,maxsys%nsx)
  ! Variables
  ! integer         :: ia,is,k
  ! do is=1,ions1%nsp
  !  do ia=1,ions0%na(is)
  !   do k=1,3
  !     vect(k,ia,is)=vect(k,ia,is)/rmass%pma(is)
  !   enddo
  !  enddo
  ! enddo
  ! ==--------------------------------------------------------------==
  ! return
  ! end subroutine mom_to_vel
  ! ==================================================================
  ! subroutine vel_to_mom(vect,itraj)
  ! ==--------------------------------------------------------------==
  ! Arguments
  ! real(real_8)    :: vect(3,maxsys%nax,maxsys%nsx,*)
  ! integer         :: ITRAJ
  ! Variables
  ! integer         :: ia,is,k
  ! do is=1,ions1%nsp
  !  do ia=1,ions0%na(is)
  !   do k=1,3
  !     vect(k,ia,is,itraj)=vect(k,ia,is,itraj)*rmass%pma(is)
  !   enddo
  !  enddo
  ! enddo
  ! ==--------------------------------------------------------------==
  ! return
  ! end subroutine vel_to_mom
  ! ==================================================================
  SUBROUTINE nabdy_ampli_update(deltat,icall)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: deltat
    INTEGER                                  :: icall

    CHARACTER(*), PARAMETER :: procedureN = 'nabdy_ampli_update'

    INTEGER                                  :: i, ia, ierr, is, j, k, npoints
    REAL(real_8)                             :: ampli, ampli_mol(1000), &
                                                distij, ni, nj, nsinv, sis, &
                                                time_l
    REAL(real_8), ALLOCATABLE                :: dydx(:), wint(:,:,:), &
                                                xint(:,:), yterm(:)
    LOGICAL                                  :: debug

    debug=.FALSE.

    IF (paral%io_parent) THEN
       WRITE(6,*) '... entered in subr. NABDY_AMPLI_UPDATE'
       !
       ALLOCATE(wint(nabdyvar%ntrajbd,maxsys%nax,maxsys%nsx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(xint(nabdyvar%ntrajbd,nabdyvar%ntrajbd),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(yterm(nabdyvar%ntrajbd),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(dydx(nabdyvar%ntrajbd),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       !   COMPUTE FIRST TERM ON TH RHS
       !   Define NPOINTS for the integration of the first term
       time_l=nabdyvar%nabdy_time
       npoints=20
       IF (npoints.GT.100) THEN
          IF (paral%io_parent) WRITE(6,*)  &
               ' too many grid points. stop in nabdy_ampli_update'
       ENDIF

       DO k=1,nabdyvar%ntrajbd 
          IF (icall.EQ.1) THEN
             CALL wintegral(wint,nabdyvec%naampl,k,npoints)
          ELSE
             CALL wintegral(wint,nabdyvec%naampl_temp,k,npoints)
          ENDIF
       ENDDO
       !
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             CALL get_grid_shepard(ia,is,npoints) ! loads gridcoor(k,i) ,k=1,3; i=gridpoint index
             IF (icall.EQ.1) THEN
                DO k=1,nabdyvar%ntrajbd
                   yterm(k)=nabdyvec%naampl(ia,is,k)
                ENDDO
                DO k=1,nabdyvar%ntrajbd
                   CALL xintegral(xint,k,ia,is,npoints)
                ENDDO
             ELSE ! (icall.eq.2)
                DO k=1,nabdyvar%ntrajbd
                   yterm(k)=nabdyvec%naampl_temp(ia,is,k)
                ENDDO
                DO k=1,nabdyvar%ntrajbd
                   CALL xintegral(xint,k,ia,is,npoints)
                ENDDO
             ENDIF
             !      the starting amplitudes are always nabdyvec%naampl(ia,is,k) [icall =1 and 2]
             DO k=1,nabdyvar%ntrajbd
                yterm(k)=nabdyvec%naampl(ia,is,k)
             ENDDO

             !ivano
             IF (debug) WRITE(6,'(/,100(1X,E12.6),/)') (nabdyvec%ddnaact(ia,is,k),k=1,nabdyvar%ntrajbd)

             CALL naderivs(time_l,yterm,dydx,wint,xint,ia,is,nabdyvar%ntrajbd,icall) !(eq40)
             CALL rungesh(yterm,dydx,nabdyvar%ntrajbd,time_l,deltat,yterm,wint,xint, &
                  ia,is,1)
             !&                ia,is,1)  why '1'? set to 'icall' 18.10.2011 (ivano)
             !&                ia,is,icall)
             !&                no! set back to 1 because in rungesh the rhs needs to be updated
             !                 every time (and not kept equal to nabdyvec%naampl_temp when icall.eq.2
             IF (icall.EQ.1) THEN
                DO k=1,nabdyvar%ntrajbd
                   nabdyvec%naampl_temp(ia,is,k)=yterm(k)
                ENDDO
             ELSE
                DO k=1,nabdyvar%ntrajbd
                   nabdyvec%naampl(ia,is,k)=yterm(k)
                ENDDO
                !       if (paral%io_parent) write(6,*) (yterm(k),k=1,nabdyvar%ntrajbd)
             ENDIF
          ENDDO
       ENDDO
       !
       !   normalization
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             ampli=0._real_8
             !      run over the fluid elements
             IF(icall.EQ.1) THEN
                DO i=1,nabdyvar%ntrajbd
                   ampli=ampli+nabdyvec%naampl_temp(ia,is,i)**2 / &
                        (8._real_8 * pi**(1.5) * nabdyvec%naomega(ia,is,i)**(3.0))
                   ni=1._real_8/((2._real_8*pi)**(1.5) * (nabdyvec%naomega(ia,is,i)**3))
                   sis=(1._real_8/SQRT(2._real_8))*nabdyvec%naomega(ia,is,i)
                   DO j=i+1,nabdyvar%ntrajbd
                      nj=ni
                      nsinv=(2._real_8*pi)**(1.5) * sis**3
                      distij=SQRT((nabdyvec%nacoor(1,ia,is,i)-nabdyvec%nacoor(1,ia,is,j))**2 + &
                           (nabdyvec%nacoor(2,ia,is,i)-nabdyvec%nacoor(2,ia,is,j))**2 + & 
                           (nabdyvec%nacoor(3,ia,is,i)-nabdyvec%nacoor(3,ia,is,j))**2 &
                           )
                      ampli=ampli+ 2._real_8 * &
                           (nabdyvec%naampl_temp(ia,is,i)*nabdyvec%naampl_temp(ia,is,j)*ni*nj) * &
                           nsinv * &
                           EXP(-((1._real_8/SQRT(2._real_8))* distij**2)/(2._real_8*sis**2)  )
                   ENDDO
                ENDDO
                DO i=1,nabdyvar%ntrajbd
                   nabdyvec%naampl_temp(ia,is,i)=nabdyvec%naampl_temp(ia,is,i)/SQRT(ampli)
                ENDDO
             ELSE
                DO i=1,nabdyvar%ntrajbd
                   ampli=ampli+nabdyvec%naampl(ia,is,i)**2 / &
                        (8._real_8 * pi**(1.5) * nabdyvec%naomega(ia,is,i)**(3.0))
                   ni=1._real_8/((2._real_8*pi)**(1.5) * (nabdyvec%naomega(ia,is,i)**3))
                   sis=(1._real_8/SQRT(2._real_8))*nabdyvec%naomega(ia,is,i)
                   DO j=i+1,nabdyvar%ntrajbd
                      nj=ni
                      nsinv=(2._real_8*pi)**(1.5) * sis**3
                      distij=SQRT((nabdyvec%nacoor(1,ia,is,i)-nabdyvec%nacoor(1,ia,is,j))**2 + &
                           (nabdyvec%nacoor(2,ia,is,i)-nabdyvec%nacoor(2,ia,is,j))**2 + &
                           (nabdyvec%nacoor(3,ia,is,i)-nabdyvec%nacoor(3,ia,is,j))**2 &
                           )
                      ampli=ampli+ 2._real_8 * &
                           (nabdyvec%naampl(ia,is,i)*nabdyvec%naampl(ia,is,j)*ni*nj) * &
                           nsinv * &
                           EXP(-((1._real_8/SQRT(2._real_8))* distij**2)/(2._real_8*sis**2)  )
                   ENDDO
                ENDDO
                DO i=1,nabdyvar%ntrajbd
                   nabdyvec%naampl(ia,is,i)=nabdyvec%naampl(ia,is,i)/SQRT(ampli)
                ENDDO
             ENDIF
             IF (debug) WRITE(6,*) 'ampli in NABDY_AMPLI_UPDATE',ampli
          ENDDO
       ENDDO
       !
       !IF (icall.EQ.2) THEN
       !   Special output for the H2 molecule
       !   DO k=1,nabdyvar%ntrajbd
       !      ampli_mol(k)=nabdyvec%naampl(1,1,k)*nabdyvec%naampl(2,1,k)
       !   ENDDO
       !   OPEN(unit=666,file='aplitudes_molecular.dat')
       !   WRITE(666,'(F12.4,100(2X,E15.9))')  &
       !        time_l,(ampli_mol(k),k=1,nabdyvar%ntrajbd)
       !   CLOSE(666)
       !ENDIF
       !   ==------------------------------------------------------------==
       !   output
       IF (icall.EQ.1) THEN
          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                WRITE(668,'(1000(E15.9,2X))')  &
                     (nabdyvec%naampl_temp(ia,is,i),i=1,nabdyvar%ntrajbd)
             ENDDO
          ENDDO
       ELSE
          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                WRITE(669,'(1000(E15.9,2X))') (nabdyvec%naampl(ia,is,i),i=1,nabdyvar%ntrajbd)
             ENDDO
          ENDDO
       ENDIF
       !   ==------------------------------------------------------------==
       DEALLOCATE(wint,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(xint,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(yterm,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(dydx,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    !
    IF (icall.EQ.1) THEN
       CALL mp_bcast(nabdyvec%naampl_temp,SIZE(nabdyvec%naampl_temp),parai%io_source,parai%cp_grp)
    ELSE
       CALL mp_bcast(nabdyvec%naampl,SIZE(nabdyvec%naampl),parai%io_source,parai%cp_grp)
    ENDIF
    !
    IF (paral%io_parent) WRITE(6,*) '... exiting in subr. NABDY_AMPLI_UPDATE'
    CALL m_flush(6)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE nabdy_ampli_update
  ! ==================================================================

  ! UTILITIES FOR  NABDY_AMPLI_UPDATE

  ! ==================================================================
  SUBROUTINE get_grid_shepard(atno,spno,npoints)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: atno, spno, npoints

    INTEGER                                  :: i, k, l
    REAL(real_8)                             :: com(3), dist, maxl(3), &
                                                startl(3)

! Variables
!
! find the coordinates of the com of all trajectories 

    DO k=1,3
       com(k)=0._real_8
    ENDDO
    DO l=1,nabdyvar%ntrajbd
       DO k=1,3
          com(k)=com(k)+nabdyvec%nacoor(k,atno,spno,l)
       ENDDO
    ENDDO
    DO k=1,3
       com(k)=com(k)/nabdyvar%ntrajbd
       maxl(k)=0._real_8
    ENDDO
    DO l=1,nabdyvar%ntrajbd
       DO k=1,3
          dist=dabs(com(k)-nabdyvec%nacoor(k,atno,spno,l))
          IF (dist.GT.maxl(k)) maxl(k)=dist
       ENDDO
    ENDDO
    ! add 50% tolerance
    DO k=1,3
       maxl(k)=2.0*(maxl(k)+0.5_real_8*maxl(k))
       startl(k)=com(k)-0.5_real_8*(maxl(k))
    ENDDO
    ! compute grid points
    nabdyvar%naintdv=1._real_8
    DO k=1,3
       nabdyvec%ddelta(k)=maxl(k)/(npoints-1)
       nabdyvar%naintdv=nabdyvar%naintdv*nabdyvec%ddelta(k)
       DO i=1,npoints
          nabdyvec%gridshep(k,i)=startl(k)+(i-1)*nabdyvec%ddelta(k)
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE get_grid_shepard
  ! ==================================================================
  SUBROUTINE rungesh(yvec,dydx,nlimit,xval,displ,yres,wint,xint, &
       ia,is,icall)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nlimit
    REAL(real_8)                             :: dydx(nlimit), yvec(nlimit), &
                                                xval, displ, yres(nlimit), &
                                                wint(*), xint(nlimit,nlimit)
    INTEGER                                  :: ia, is, icall

    INTEGER, PARAMETER                       :: nmax = 1000  

    INTEGER                                  :: i
    REAL(real_8)                             :: dy_t(nmax), dym(nmax), inth2, &
                                                inth6, x_h, y_t(nmax)

    inth2=displ*0.5_real_8
    inth6=displ/6.0_real_8
    x_h=xval+inth2
    DO i=1,nlimit 
       y_t(i)=yvec(i)+inth2*dydx(i)
    ENDDO
    CALL naderivs(x_h,y_t,dy_t,wint,xint,ia,is,nlimit,icall) 
    DO i=1,nlimit
       y_t(i)=yvec(i)+inth2*dy_t(i)
    ENDDO
    CALL naderivs(x_h,y_t,dym,wint,xint,ia,is,nlimit,icall) 
    DO i=1,nlimit
       y_t(i)=yvec(i)+displ*dym(i)
       dym(i)=dy_t(i)+dym(i)
    ENDDO
    CALL naderivs(xval+displ,y_t,dy_t,wint,xint,ia,is,nlimit,icall) 
    DO i=1,nlimit 
       yres(i)=yvec(i)+inth6*(dydx(i)+dy_t(i)+2.*dym(i))
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rungesh
  ! ==================================================================
  SUBROUTINE naderivs(time_l,yterm,dydx,wint,xint,ia,is,nlimit,icall)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: time_l
    INTEGER                                  :: ia, is, nlimit
    REAL(real_8) :: xint(nlimit,nlimit), wint(nlimit,maxsys%nax,maxsys%nsx), &
      dydx(nlimit), yterm(nlimit)
    INTEGER                                  :: icall

    INTEGER                                  :: ia1, is1, j, k, l
    REAL(real_8)                             :: aa, aaj, comega, distij, nj, &
                                                nk, nsinv, sis

    DO k=1,nlimit
       dydx(k)=0.0_real_8
    ENDDO

    IF (icall.EQ.1) THEN
       DO k=1,nlimit
          DO l=1,nlimit
             DO is1=1,ions1%nsp
                IF (ions0%iatyp(is1).GT.nabdyvar%nabdy_zmax) CYCLE
                DO ia1=1,ions0%na(is1)
                   dydx(k)=dydx(k)-(yterm(k)/(2.0_real_8*rmass%pma(is1)))* &
                        wint(l,ia1,is1)*nabdyvec%ddnaact(ia1,is1,l)
                ENDDO
             ENDDO
             dydx(k)=dydx(k)-(1.0_real_8/(2.0_real_8*rmass%pma(is)))* &
                  (yterm(k)/nabdyvec%nanormalf(ia,is,k))  * &
                  xint(k,l)*nabdyvec%ddnaact(ia,is,l)
          ENDDO
       ENDDO

    ELSEIF(icall.EQ. 2) THEN

       DO k=1,nlimit
          DO l=1,nlimit
             DO is1=1,ions1%nsp
                IF (ions0%iatyp(is1).GT.nabdyvar%nabdy_zmax) CYCLE
                DO ia1=1,ions0%na(is1)
                   dydx(k)=dydx(k)-(nabdyvec%naampl_temp(ia,is,k)/(2.0_real_8*rmass%pma(is1)))* &
                        wint(l,ia1,is1)*nabdyvec%ddnaact(ia1,is1,l)
                ENDDO
             ENDDO
             dydx(k)=dydx(k)-(1.0_real_8/(2.0_real_8*rmass%pma(is)))* &
                  (nabdyvec%naampl_temp(ia,is,k)/nabdyvec%nanormalf(ia,is,k))* &
                  xint(k,l)*nabdyvec%ddnaact(ia,is,l)
          ENDDO
       ENDDO
    ENDIF
    !
    ! Jones normalization
    ! F= sum_i a_i^2/(8 pi^(3/2) sigma_i^3) =1
    ! F_{a_1} = a_1/(4 pi^(3/2) sigma_i^3)
    !
    GOTO 100

    comega=0._real_8
    DO k=1,nlimit
       IF (icall.EQ.1) THEN
          aa=yterm(k)
       ELSE
          aa=nabdyvec%naampl_temp(ia,is,k)
       ENDIF
       comega=comega+dydx(k)* &
            (aa/(4.0_real_8 * pi**(1.5) * nabdyvec%naomega(ia,is,k)**(3.0)))
       nk=1._real_8/((2._real_8*pi)**(1.5) * (nabdyvec%naomega(ia,is,k)**3))
       sis=(1._real_8/SQRT(2._real_8))*nabdyvec%naomega(ia,is,k)
       DO j=k+1,nabdyvar%ntrajbd
          IF (icall.EQ.1) THEN
             aaj=yterm(j)
          ELSE
             aaj=nabdyvec%naampl_temp(ia,is,j)
          ENDIF
          nj=nk
          nsinv=(2._real_8*pi)**(1.5) * sis**3
          distij=SQRT((nabdyvec%nacoor(1,ia,is,k)-nabdyvec%nacoor(1,ia,is,j))**2 + &
               (nabdyvec%nacoor(2,ia,is,k)-nabdyvec%nacoor(2,ia,is,j))**2 + &
               (nabdyvec%nacoor(3,ia,is,k)-nabdyvec%nacoor(3,ia,is,j))**2 &
               )
          comega=comega+ dydx(k)* 2._real_8 * &
               (aaj*nk*nj) * &
               nsinv * &
               EXP(-((1._real_8/SQRT(2._real_8))* distij**2)/(2._real_8*sis**2)  )
       ENDDO
    ENDDO
    comega=comega/2._real_8
    !
    DO k=1,nlimit
       IF (icall.EQ.1) THEN
          aa=yterm(k)
       ELSE
          aa=nabdyvec%naampl_temp(ia,is,k)
       ENDIF
       dydx(k)=dydx(k)-aa*comega
    ENDDO

100 CONTINUE

    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE naderivs
  ! ==================================================================
  SUBROUTINE wintegral(wint,ampli,l,npoints)
    ! ==--------------------------------------------------------------==
    REAL(real_8) :: wint(nabdyvar%ntrajbd,maxsys%nax,*), &
      ampli(maxsys%nax,maxsys%nsx,*)
    INTEGER                                  :: l, npoints

    REAL(real_8), PARAMETER                  :: p_shep = 2.0

    INTEGER                                  :: ia, is, ix, iy, iz, k, p
    LOGICAL                                  :: add_it
    REAL(real_8)                             :: distk, distk2, distl, distl2, &
                                                distp, distp2, minddelta, &
                                                rcoord(3), sumoverk, sump

    minddelta=1000.0
    DO k=1,3
       IF (nabdyvec%ddelta(k).LT.minddelta) minddelta=nabdyvec%ddelta(k)
    ENDDO
    !
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)

          wint(l,ia,is)=0._real_8

          DO ix=1,npoints
             DO iy=1,npoints
                DO iz=1,npoints
                   !
                   add_it=.TRUE.
                   !
                   rcoord(1)=nabdyvec%gridshep(1,ix)
                   rcoord(2)=nabdyvec%gridshep(2,iy)
                   rcoord(3)=nabdyvec%gridshep(3,iz)

                   sumoverk=0._real_8
                   DO k=1,nabdyvar%ntrajbd
                      distk2=(rcoord(1)-nabdyvec%nacoor(1,ia,is,k))**2 + &
                           (rcoord(2)-nabdyvec%nacoor(2,ia,is,k))**2 + &
                           (rcoord(3)-nabdyvec%nacoor(3,ia,is,k))**2
                      distk=SQRT(distk2)
                      !        --------------------------------------------------------
                      !        skip if the x,y,z point sits close to a grid point
                      IF (distk.LT.minddelta/2.0) add_it=.FALSE.
                      !        --------------------------------------------------------
                      sumoverk=sumoverk+ &
                           (ampli(ia,is,k)/nabdyvec%nanormalf(ia,is,k))* &
                           EXP(-distk2/(2.0_real_8*nabdyvec%naomega(ia,is,k)**2))
                   ENDDO
                   sumoverk=sumoverk**2

                   distl2=(rcoord(1)-nabdyvec%nacoor(1,ia,is,l))**2 + &
                        (rcoord(2)-nabdyvec%nacoor(2,ia,is,l))**2 + &
                        (rcoord(3)-nabdyvec%nacoor(3,ia,is,l))**2
                   distl=SQRT(distl2)
                   !
                   sump=0._real_8
                   DO p=1,nabdyvar%ntrajbd
                      distp2=(rcoord(1)-nabdyvec%nacoor(1,ia,is,p))**2 + &
                           (rcoord(2)-nabdyvec%nacoor(2,ia,is,p))**2 + &
                           (rcoord(3)-nabdyvec%nacoor(3,ia,is,p))**2
                      distp=SQRT(distp2)
                      sump=sump+1.0/((distp)**(p_shep))
                   ENDDO
                   !
                   IF (add_it) THEN
                      wint(l,ia,is)=wint(l,ia,is)+ &
                           sumoverk &
                           / &
                           ((distl**p_shep)*sump)* &
                           nabdyvar%naintdv
                   ENDIF
                   !
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE wintegral
  ! ==================================================================
  SUBROUTINE xintegral(xint,k,ia,is,npoints)
    ! ==--------------------------------------------------------------==
    REAL(real_8) :: xint(nabdyvar%ntrajbd,nabdyvar%ntrajbd)
    INTEGER                                  :: k, ia, is, npoints

    INTEGER                                  :: ix, iy, iz, l, l1, l2, p
    LOGICAL                                  :: debug
    REAL(real_8)                             :: distk2, distl2, distp2, &
                                                minddelta, rcoord(3), sump

    debug=.FALSE.

    DO l=1,nabdyvar%ntrajbd
       xint(k,l)=0._real_8
    ENDDO

    DO l=1,nabdyvar%ntrajbd

       DO ix=1,npoints
          DO iy=1,npoints
             DO iz=1,npoints

                rcoord(1)=nabdyvec%gridshep(1,ix)
                rcoord(2)=nabdyvec%gridshep(2,iy)
                rcoord(3)=nabdyvec%gridshep(3,iz)
                distl2=(rcoord(1)-nabdyvec%nacoor(1,ia,is,l))**2 + &
                     (rcoord(2)-nabdyvec%nacoor(2,ia,is,l))**2 + &
                     (rcoord(3)-nabdyvec%nacoor(3,ia,is,l))**2
                distk2=(rcoord(1)-nabdyvec%nacoor(1,ia,is,k))**2 + &
                     (rcoord(2)-nabdyvec%nacoor(2,ia,is,k))**2 + &
                     (rcoord(3)-nabdyvec%nacoor(3,ia,is,k))**2
                !        --------------------------------------------------------
                minddelta=1000.0
                DO l1=1,3
                   IF (nabdyvec%ddelta(l1).LT.minddelta) minddelta=nabdyvec%ddelta(l1)
                ENDDO
                IF ((SQRT(distl2).LT.minddelta)) CYCLE
                IF ((SQRT(distk2).LT.minddelta)) CYCLE
                !        --------------------------------------------------------

                sump=0._real_8
                DO p=1,nabdyvar%ntrajbd
                   distp2=(rcoord(1)-nabdyvec%nacoor(1,ia,is,p))**2 + &
                        (rcoord(2)-nabdyvec%nacoor(2,ia,is,p))**2 + &
                        (rcoord(3)-nabdyvec%nacoor(3,ia,is,p))**2
                   !         --------------------------------------------------------
                   minddelta=1000.0
                   DO l1=1,3
                      IF (nabdyvec%ddelta(l1).LT.minddelta) minddelta=nabdyvec%ddelta(l1)
                   ENDDO
                   IF ((SQRT(distp2).LT.minddelta)) CYCLE
                   !         --------------------------------------------------------
                   sump=sump+1.0/distp2
                ENDDO

                xint(k,l)=xint(k,l)+ &
                     ((EXP(-distk2/(2._real_8 * nabdyvec%naomega(ia,is,k)**2 ))) &
                     /(distl2*sump))*nabdyvar%naintdv

             ENDDO
          ENDDO
       ENDDO

    ENDDO
    IF (paral%io_parent.AND.debug) &
         WRITE(6,'(100E14.6)') (xint(k,l2),l2=1,nabdyvar%ntrajbd)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE xintegral
  ! ==================================================================
  SUBROUTINE write_nabdy_traj
    ! ==--------------------------------------------------------------==
    ! To be called in mdnabby.F because it needs the temperature
    ! If temp_inst.eq.0 then it takes tempw from system.h
    ! real(real_8)    ::     temp_inst
    ! Variables
    CHARACTER(*), PARAMETER :: procedureN = 'write_nabdy_traj'

    CHARACTER(LEN=14)                        :: filen1, filen4
    CHARACTER(LEN=28)                        :: filen2, filen3
    INTEGER                                  :: i, ia, ierr, is, j, k
    INTEGER, ALLOCATABLE                     :: gr_iat(:)
    LOGICAL                                  :: status
    REAL(real_8)                             :: mass
    REAL(real_8), ALLOCATABLE                :: ampli(:), coord(:,:), &
                                                sigma(:), veloc(:,:)

    CALL mm_dim(mm_go_mm,status)

    ALLOCATE(coord(3,nabdyvar%ntrajbd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(veloc(3,nabdyvar%ntrajbd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(ampli(nabdyvar%ntrajbd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(sigma(nabdyvar%ntrajbd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(gr_iat(ions1%nat+1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    ! from geofile.F
    i=0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          i=i+1
          gr_iat(nat_grm(i))=ions0%iatyp(is)
       ENDDO
    ENDDO
    !
    filen1="nabdy_geometry"
    filen2="nabdy_geometry.xyz"
    filen3="nabdy_gaussian.dat"
    filen4='nabdy_time.dat'

    IF (paral%io_parent) THEN

       OPEN(unit=777,file=filen1,status='replace')
       OPEN(unit=888,file=filen3,status='replace')
       OPEN(unit=999,file=filen4,status='old',position='append')
       DO j=1,nabdyvar%ntrajbd
          IF(cntl%tqmmm)THEN
             DO i=1,ions1%nat
                is=cpsp(i)
                mass=rmass%pma(is)
                WRITE(777,'(3f20.12,8x,3f20.12)') &
                     (coord(k,i),k=1,3),(veloc(k,i),k=1,3)
                WRITE(888,'(2e20.12)') ampli(i),sigma(i)
                DO k=1,3
                   nabdyvec%nacoor(k,cpat(i),cpsp(i),j)=coord(k,i)+clsaabox%mm_c_trans(k)
                   nabdyvec%namom(k,cpat(i),cpsp(i),j)=veloc(k,i)*mass
                END DO
                nabdyvec%naampl(cpat(i),cpsp(i),j)=ampli(i)
                nabdyvec%naomega(cpat(i),cpsp(i),j)=sigma(i)
             ENDDO
          ELSE
             DO is=1,ions1%nsp
                mass=rmass%pma(is)
                DO ia=1,ions0%na(is)
                   WRITE(777,'(3f20.12,8x,3f20.12)') &
                        (nabdyvec%nacoor(k,ia,is,j),k=1,3),(nabdyvec%namom(k,ia,is,j)/mass,k=1,3)
                   WRITE(888,'(2e20.12)') nabdyvec%naampl(ia,is,j),nabdyvec%naomega(ia,is,j)
                ENDDO
             ENDDO
          ENDIF
       ENDDO
       WRITE(999,'(F12.4)') nabdyvar%nabdy_time
       !
       CLOSE(777)
       CLOSE(888)
       CLOSE(999)
    ENDIF
    !
    CALL mm_dim(mm_revert,status)
    DEALLOCATE(coord,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(veloc,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(ampli,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(sigma,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(gr_iat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE write_nabdy_traj
  ! ==================================================================
  SUBROUTINE integrate_atom(atno,spno,npoints)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: atno, spno, npoints

    INTEGER                                  :: i, ix, iy, iz, k, l
    REAL(real_8)                             :: com(3), dist, dist2, &
                                                integral, maxl(3), rcoord(3), &
                                                startl(3)
    LOGICAL                                  :: debug

    debug=.FALSE.
! Variables
!
! find the coordinates of the com of all trajectories 

    DO k=1,3
       com(k)=0._real_8
    ENDDO
    DO l=1,nabdyvar%ntrajbd
       DO k=1,3
          com(k)=com(k)+nabdyvec%nacoor(k,atno,spno,l)
       ENDDO
    ENDDO
    DO k=1,3
       com(k)=com(k)/nabdyvar%ntrajbd
       maxl(k)=0.0
    ENDDO
    DO l=1,nabdyvar%ntrajbd
       DO k=1,3
          dist=dabs(com(k)-nabdyvec%nacoor(k,atno,spno,l))
          IF (dist.GT.maxl(k)) maxl(k)=dist
       ENDDO
    ENDDO
    ! add 100% tolerance
    DO k=1,3
       maxl(k)=2.0*(maxl(k)+1._real_8*maxl(k))
       startl(k)=com(k)-0.5_real_8*(maxl(k))
    ENDDO
    ! compute grid points
    nabdyvar%naintdv=1._real_8
    DO k=1,3
       nabdyvec%ddelta(k)=maxl(k)/(npoints-1)
       nabdyvar%naintdv=nabdyvar%naintdv*nabdyvec%ddelta(k)
       DO i=1,npoints
          nabdyvec%gridshep(k,i)=startl(k)+(i-1)*nabdyvec%ddelta(k)
       ENDDO
    ENDDO

    integral=0._real_8

    DO l=1,nabdyvar%ntrajbd

       DO ix=1,npoints
          DO iy=1,npoints
             DO iz=1,npoints
                rcoord(1)=nabdyvec%gridshep(1,ix)
                rcoord(2)=nabdyvec%gridshep(2,iy)
                rcoord(3)=nabdyvec%gridshep(3,iz)
                dist2=(rcoord(1)-nabdyvec%nacoor(1,atno,spno,l))**2 + &
                     (rcoord(2)-nabdyvec%nacoor(2,atno,spno,l))**2 + &
                     (rcoord(3)-nabdyvec%nacoor(3,atno,spno,l))**2
                integral=integral+( &
                     (nabdyvec%naampl(atno,spno,l)/nabdyvec%nanormalf(atno,spno,l))*  &
                     EXP(-dist2/(2._real_8 * nabdyvec%naomega(atno,spno,l)**2 )) &
                     )**2  &
                     *nabdyvar%naintdv   
             ENDDO
          ENDDO
       ENDDO

    ENDDO

    IF (debug) WRITE(6,'(A,2I6,F12.7)') 'integtal',atno,spno,integral
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE integrate_atom
  ! ==================================================================
  SUBROUTINE check_small_ampli
    ! ==--------------------------------------------------------------==
    ! Variables
    INTEGER                                  :: i, i_mem(1000), ia, &
                                                ia_mem(1000), icount, ind, &
                                                is, is_mem(1000), j, k
    REAL(real_8) :: a_mean, a_std, a_threshold, com(3), coord(3), dist, &
      dist_threshold, dx_variance, min_dist, x0, x1, z1, z2
    LOGICAL                                  :: success,try_once_more 

    a_threshold=1.0D-5
    dx_variance=0.05 ! au
    dist_threshold=0.05

    icount=0
    DO i=1,nabdyvar%ntrajbd
       DO is=1,ions1%nsp
          IF (rmass%pma(is).GT.nabdyvar%nabdy_zmax) CYCLE
          DO ia=1,ions0%na(is)
             IF (nabdyvec%naampl(ia,is,i).LT.a_threshold) THEN
                icount=icount+1
                ia_mem(icount)=ia
                is_mem(icount)=is
                i_mem(icount)=i
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    !
    IF (paral%io_parent) THEN
       WRITE(6,'(A,I5)') 'Number of redistributed FEs',icount
       !
       DO ind=1,icount
          !    compute the center of mass of the atomic distribution
          DO k=1,3
             com(k)=0.0
          ENDDO
          DO i=1,nabdyvar%ntrajbd
             DO k=1,3
                com(k)=com(k)+nabdyvec%nacoor(k,ia_mem(ind),is_mem(ind),i)
             ENDDO
          ENDDO
          DO k=1,3
             com(k)=com(k)/nabdyvar%ntrajbd
          ENDDO
          !    ------------------------------------------------------
          !   try to move the atom within a distance dx from the com
          !   iterate until success
 333      CONTINUE
          DO j=1,100
             DO k=1,3
                z1=repprngu()
                z2=repprngu()
                x0=com(k)
                x1=(-2._real_8*LOG(z1))**0.5 * COS(2._real_8*pi*z2)
                coord(k)=x0+dx_variance*x1
             ENDDO
             min_dist=1000.0
             DO i=1,nabdyvar%ntrajbd
                dist=(nabdyvec%nacoor(1,ia_mem(ind),is_mem(ind),i)-coord(1))**2 + &
                     (nabdyvec%nacoor(2,ia_mem(ind),is_mem(ind),i)-coord(2))**2 + &
                     (nabdyvec%nacoor(3,ia_mem(ind),is_mem(ind),i)-coord(3))**2
                IF (dist.LT.min_dist) min_dist=dist
             ENDDO
             min_dist=SQRT(min_dist)
 
             IF (min_dist.GT.dist_threshold) THEN
               ! succeeded
               success=.TRUE.
             ELSE
               IF (dist_threshold.LT.1.0D-3) THEN
                ! did not succeed and giveup
                success=.FALSE.
                try_once_more=.FALSE.
               ELSE
                ! try again with halph dist_threshold
                success=.FALSE.
                try_once_more=.TRUE.
                dist_threshold=dist_threshold/2.0_real_8
               ENDIF
             ENDIF
          ENDDO

          IF ((.NOT.success).AND.(try_once_more)) THEN
           GOTO 333
          ENDIF

          IF (success) THEN
            WRITE(6,*) '  .... one fluid element redistributed'  
            DO k=1,3
              nabdyvec%nacoor(k,ia_mem(ind),is_mem(ind),i_mem(ind))=coord(k)
            ENDDO
            !   ------------------------------------------------------
            !   assign an amplitude
            a_mean=0.0
            DO i=1,nabdyvar%ntrajbd
               a_mean=a_mean+nabdyvec%naampl(ia_mem(ind),is_mem(ind),i)
            ENDDO
            a_mean=a_mean/nabdyvar%ntrajbd
            a_std=0.0
            DO i=1,nabdyvar%ntrajbd
               a_std=a_std+(nabdyvec%naampl(ia_mem(ind),is_mem(ind),i)-a_mean)**2
            ENDDO
            a_std=SQRT(a_std/nabdyvar%ntrajbd)
            nabdyvec%naampl(ia_mem(ind),is_mem(ind),i_mem(ind))=a_mean+a_std
            !   ------------------------------------------------------
          ENDIF
       ENDDO
    ENDIF
    !
    ! distribute coordinates to all processors
    CALL mp_bcast(nabdyvec%nacoor,SIZE(nabdyvec%nacoor),parai%io_source,parai%cp_grp)
    CALL mp_bcast(nabdyvec%naampl,SIZE(nabdyvec%naampl),parai%io_source,parai%cp_grp)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE check_small_ampli
  ! ==================================================================
  SUBROUTINE nabdy_friction_cm(navelp,nstate)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: navelp(:,:,:,:)
    INTEGER                                  :: nstate

    CHARACTER(*), PARAMETER :: procedureN = 'nabdy_friction_cm'
    REAL(real_8), PARAMETER                  :: thresh = 1.0e-6_real_8  

    CHARACTER(LEN=16)                        :: filen
    INTEGER                                  :: ia, ierr, is, istart, it, k
    INTEGER, SAVE                            :: ifirst = 0
    LOGICAL                                  :: fexist
    REAL(real_8) :: dispp, dispp_r, dispx, dpmax, dpmin, dpmod, ekintmp, &
      lambda, lambda_p, scalp, scalp2, scalx, scalx2, tempcur, xx
    REAL(real_8), ALLOCATABLE                :: dp(:,:,:,:), dx(:,:,:,:), &
                                                meanp(:,:), meanx(:,:), &
                                                velpp(:,:,:)

    istart=0
    filen="nabdy_thermostat"

    ALLOCATE(velpp(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(meanx(3,maxsys%nax),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(meanp(3,maxsys%nax),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(dx(3,maxsys%nax,maxsys%nsx,nabdyvar%ntrajbd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(dp(3,maxsys%nax,maxsys%nsx,nabdyvar%ntrajbd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    ! ==--------------------------------------------------------------==
    DO is=1,ions1%nsp
       IF (((ions0%iatyp(is).GT.nabdyvar%nabdy_zmax)).OR.(ifirst.LT.istart)) CYCLE

       CALL zeroing(nabdyfric%nabdy_dx)!,maxsys%nsx*maxsys%nax)
       CALL zeroing(nabdyfric%nabdy_dx2)!,maxsys%nsx*maxsys%nax)
       CALL zeroing(nabdyfric%nabdy_dp)!,maxsys%nsx*maxsys%nax)
       CALL zeroing(nabdyfric%nabdy_dp2)!,maxsys%nsx*maxsys%nax)

       DO ia=1,ions0%na(is)
          DO k=1,3
             meanx(k,ia)=0._real_8
             meanp(k,ia)=0._real_8
          ENDDO
       ENDDO
       DO ia=1,ions0%na(is)
          DO it=1,nabdyvar%ntrajbd
             DO k=1,3 
                meanx(k,ia)=meanx(k,ia)+nabdyvec%nacoor(k,ia,is,it) 
                meanp(k,ia)=meanp(k,ia)+nabdyvec%namom(k,ia,is,it)
             ENDDO
          ENDDO
          DO k=1,3
             meanx(k,ia)=meanx(k,ia)/nabdyvar%ntrajbd
             meanp(k,ia)=meanp(k,ia)/nabdyvar%ntrajbd
          ENDDO
          DO k=1,3
             velpp(k,ia,is)=meanp(k,ia)/rmass%pma(is)
          ENDDO
       ENDDO

       DO ia=1,ions0%na(is)
          DO it=1,nabdyvar%ntrajbd
             DO k=1,3
                dx(k,ia,is,it)=nabdyvec%nacoor(k,ia,is,it)- meanx(k,ia)
                dp(k,ia,is,it)=nabdyvec%namom(k,ia,is,it) - meanp(k,ia)
             ENDDO
             scalx2=0._real_8
             scalp2=0._real_8
             DO k=1,3
                scalx2=scalx2+dx(k,ia,is,it)**2
                scalp2=scalp2+dp(k,ia,is,it)**2
             ENDDO
             scalx=SQRT(scalx2)
             scalp=SQRT(scalp2)
             nabdyfric%nabdy_dx(is,ia) =nabdyfric%nabdy_dx(is,ia)  + scalx /nabdyvar%ntrajbd
             nabdyfric%nabdy_dx2(is,ia)=nabdyfric%nabdy_dx2(is,ia) + scalx2/nabdyvar%ntrajbd
             nabdyfric%nabdy_dp(is,ia) =nabdyfric%nabdy_dp(is,ia)  + scalp /nabdyvar%ntrajbd
             nabdyfric%nabdy_dp2(is,ia)=nabdyfric%nabdy_dp2(is,ia) + scalp2/nabdyvar%ntrajbd
          ENDDO
       ENDDO
    ENDDO
    !
    IF (ifirst.EQ.istart) THEN
       IF (paral%io_parent) THEN
          INQUIRE(file=filen,exist=fexist)
          IF(.NOT.fexist) THEN
             OPEN(unit=777,file=filen,status='new')
             DO is=1,ions1%nsp
                DO ia=1,ions0%na(is)
                   nabdyfric%dispx_ref(is,ia)=nabdyfric%nabdy_dx2(is,ia)-nabdyfric%nabdy_dx(is,ia)**2
                   nabdyfric%dispp_ref(is,ia)=nabdyfric%nabdy_dp2(is,ia)-nabdyfric%nabdy_dp(is,ia)**2
                   WRITE(777,'(2f25.12)') nabdyfric%dispx_ref(is,ia),nabdyfric%dispp_ref(is,ia)
                ENDDO
             ENDDO
          ELSE
             OPEN(unit=777,file=filen,status='old')
             DO is=1,ions1%nsp
                DO ia=1,ions0%na(is)
                   READ(777,*) nabdyfric%dispx_ref(is,ia),nabdyfric%dispp_ref(is,ia)
                ENDDO
             ENDDO
          ENDIF
          CLOSE(777)
       ENDIF
       CALL mp_bcast(nabdyfric%dispx_ref,SIZE(nabdyfric%dispx_ref),parai%io_source,parai%cp_grp)
       CALL mp_bcast(nabdyfric%dispp_ref,SIZE(nabdyfric%dispp_ref),parai%io_source,parai%cp_grp)
    ENDIF
    !
    CALL ekinpp_qm(ekintmp,velpp)
    tempcur=ekintmp*factem*2.0_real_8/glib
    lambda_p=1.0_real_8
    lambda=1.0_real_8

    DO is=1,ions1%nsp
       IF (((ions0%iatyp(is).GT.nabdyvar%nabdy_zmax)).OR.(ifirst.LT.istart)) CYCLE
       !  NATAUBP=100.0
       !  TEMPW=300.0
       lambda_p=SQRT(1.0_real_8+0.5_real_8*dt_ions/nabdyfric%natempbp * &
            (nabdyfric%natempcm/tempcur-1.0_real_8))
       DO ia=1,ions0%na(is)
          !   SAKURAI p.57 for dispersion of a Gaussian Wavepacket in x and p
          dispx=nabdyfric%nabdy_dx2(is,ia)-nabdyfric%nabdy_dx(is,ia)**2
          dispp=nabdyfric%nabdy_dp2(is,ia)-nabdyfric%nabdy_dp(is,ia)**2
          dispp_r=nabdyfric%dispp_ref(is,ia)

          !   EMPW/TEMPCUR is replaced by (DISPP_target=1/(4*DISPX))/(DISPP_actual=DISPP)
          !   according to Eq. 1.7.41 in Sakurai p.58

          lambda=SQRT(1.0_real_8+0.5_real_8*dt_ions/(5.0*nabdyfric%natempbp)* &
               ((dispp_r/dispp)-1.0_real_8))

          IF (paral%io_parent) WRITE(*,'(a,10f12.3)')                            &
               'friction',nabdyfric%natempcm        , tempcur,             &
                                !                          nabdyfric%natempcm/tempcur, 0.25/(dispx*dispp),
               nabdyfric%natempcm/tempcur, dispp_r,dispp,       &
               lambda_p        , lambda
          DO k=1,3
             velpp(k,ia,is)=lambda_p*velpp(k,ia,is) 
          ENDDO
          dpmin=1.0e+10_real_8
          dpmax=0.0
          DO it=1,nabdyvar%ntrajbd
             dpmod=0._real_8
             DO k=1,3
                dpmod=dpmod+dp(k,ia,is,it)**2
             ENDDO
             dpmod=SQRT(dpmod)
             IF (dpmod.LT.dpmin) dpmin=dpmod
             IF (dpmod.GT.dpmax) dpmax=dpmod
          ENDDO

          xx=1.0_real_8
          !   LINEAR ---------------------------------------------------------
          !   XX=1 if DP=DPMIN
          !   XX=0 if DP=DPMAX
          !   XX   linearly interpolated in between
          !   with the y=0 axes tp DPMAX*1.5
          !   to avoid that dp at dpmax is multiplied by 0 I shift the cross 
          !        DPMAX=DPMAX*1.5
          !        XX=-(1._real_8/(DPMAX-DPMIN))*DPMOD+(DPMAX/(DPMAX-DPMIN))

          !   PARABOLIC ------------------------------------------------------
          !   SS=1.e-4_real_8
          !   AA=(SS-1.0-DPMAX**2+DPMIN**2)/(2._real_8*(DPMIN-DPMAX))
          !   BB=1.0-DPMIN**2+2._real_8*AA*DPMIN-AA**2
          !   XX=(DPMOD-AA)**2 + BB

          DO it=1,nabdyvar%ntrajbd
             nabdyvec%namom(1,ia,is,it)=velpp(1,ia,is)*rmass%pma(is)+lambda*xx        &
                  *dp(1,ia,is,it)   
             nabdyvec%namom(2,ia,is,it)=velpp(2,ia,is)*rmass%pma(is)+lambda*xx        &
                  *dp(2,ia,is,it)
             nabdyvec%namom(3,ia,is,it)=velpp(3,ia,is)*rmass%pma(is)+lambda*xx        &
                  *dp(3,ia,is,it)
             navelp(1,ia,is,it)=nabdyvec%namom(1,ia,is,it)/rmass%pma(is)
             navelp(2,ia,is,it)=nabdyvec%namom(2,ia,is,it)/rmass%pma(is)
             navelp(3,ia,is,it)=nabdyvec%namom(3,ia,is,it)/rmass%pma(is)
          ENDDO
       ENDDO
    ENDDO
    !
    ifirst=ifirst+1
    IF (ifirst.GT.1000) ifirst=1000
    ! endif
    !
    DEALLOCATE(velpp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(meanx,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(meanp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(dx,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(dp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE nabdy_friction_cm
  ! ==================================================================
  SUBROUTINE ekinpp_qm(ekinp,velp)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: ekinp, velp(:,:,:)

    INTEGER                                  :: j, k
    REAL(real_8)                             :: const

! Variables
! ==--------------------------------------------------------------==
! ==  calculate kinetic energy of the quantum nuclei              ==
! ==--------------------------------------------------------------==

    ekinp=0._real_8
    DO k=1,ions1%nsp
       IF ((ions0%iatyp(k).GT.nabdyvar%nabdy_zmax)) CYCLE
       const=0.5_real_8*rmass%pma(k)
       DO j=1,ions0%na(k)
          ekinp=ekinp+const*(velp(1,j,k)**2+ &
               velp(2,j,k)**2+ &
               velp(3,j,k)**2)
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE ekinpp_qm
  ! ==================================================================

END MODULE nabdy_md
