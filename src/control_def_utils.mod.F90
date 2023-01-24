MODULE control_def_utils
  USE andr,                            ONLY: andr2,&
                                             andr3
  USE atwf,                            ONLY: dmovrmix,&
                                             tmovr
  USE benc,                            ONLY: ibench
  USE broy,                            ONLY: broy1
  USE cdftmod,                         ONLY: &
       cdftci, cdftcom, cdftlog, cdftpred, cdftvgrs, cm_dir, czones, projlog, &
       sccomm, wg_n
  USE cnst_dyn,                        ONLY: lmeta
  USE comvelmod,                       ONLY: comvl
  USE conv,                            ONLY: nac,&
                                             nbc
  USE cotr,                            ONLY: sdpl
  USE cp_cuda_types,                   ONLY: cp_cuda_env
  USE fileopenmod,                     ONLY: fo_info
  USE fint,                            ONLY: fint1,&
                                             fint4,&
                                             fint5
  USE g_loc,                           ONLY: glocal,&
                                             gloci,&
                                             glocr
  USE glemod,                          ONLY: glepar,&
                                             tglepc
  USE ions,                            ONLY: coord_fdiff,&
                                             r_fdiff,&
                                             tref_fdiff
  USE isos,                            ONLY: isos1
  USE kpts,                            ONLY: tkpts
  USE linres,                          ONLY: tshl,&
                                             xfmqc
  USE machine,                         ONLY: m_getenv
  USE nort,                            ONLY: nort_com
  USE nose,                            ONLY: loct,&
                                             nosl,&
                                             tcafes,&
                                             tnosepc
  USE parac,                           ONLY: parai
  USE prden,                           ONLY: elfcb,&
                                             numpr
  USE prop,                            ONLY: prop5
  USE qspl,                            ONLY: qspl1,&
                                             qsrang
  USE readmod,                         ONLY: isrlength
  USE spin,                            ONLY: clsd,&
                                             lspin1,&
                                             lspin2,&
                                             lspin3
  USE store_types,                     ONLY: cprint,&
                                             iface1,&
                                             intfn,&
                                             iprint_max,&
                                             restart1,&
                                             rout1,&
                                             store1
  USE struc,                           ONLY: angle,&
                                             bond,&
                                             dihedral
  USE system,                          ONLY: &
       acc, cnti, cntl, cntr, cp_trace, group, locpot2, nacc, nacx, nssel, &
       restf
  USE time,                            ONLY: tname
  USE vdwcmod,                         ONLY: vdwl
  USE wann,                            ONLY: wan05,&
                                             wannc,&
                                             wanni,&
                                             wannl,&
                                             wannr
  USE xinr,                            ONLY: gnx_inr,&
                                             inr_integer,&
                                             inr_logical,&
                                             real_8,&
                                             rmixsd,&
                                             tol_inr,&
                                             tolx_inr
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: control_def

CONTAINS

  ! ==================================================================
  SUBROUTINE control_def
    ! ==--------------------------------------------------------------==
    ! ==  THIS ROUTINE SETS THE DEFAULT FOR VARIABLES IN CONTROL      ==
    ! ==--------------------------------------------------------------==
    ! ==  DEFAULTS                                                    ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: i

    nort_com%slimit=0._real_8
    nort_com%scond=0._real_8
    cntl%is_in_stream=.FALSE.
    cntl%is_out_stream=.FALSE.
    cntl%use_mpi_io=.FALSE.
    cntl%md=.FALSE.
    cntl%tmdbo=.FALSE.
    cntl%tmdfile=.FALSE.
    cntl%tprcp=.FALSE.
    cntl%geopt=.FALSE.
    cntl%wfopt=.FALSE.
    cntl%ksener=.FALSE.
    cntl%vibrat=.FALSE.
    cntl%tsdan=.FALSE.
    cntl%tspec=.FALSE.
    cntl%tsoc=.FALSE.
    cntl%tnabdy=.FALSE.
    cntl%tshop=.FALSE.
    cntl%tqmmm=.FALSE.
    cntl%tqmmech=.FALSE.
    cntl%tsampl=.FALSE.
    cntl%proper=.FALSE.
    cntl%use_mts=.false.
    restart1%restart=.FALSE.
    cntl%cnstmd=.TRUE.
    cntl%gaudyn=.FALSE.
    cntl%orbrot=.FALSE.
    cntl%nonort=.FALSE.
    cntl%tlsd=.FALSE.
    cntl%bsymm=.FALSE.
    clsd%nlsd=1
    clsd%nlsx=1
    cntl%tharm=.FALSE.
    cntl%tmass=.FALSE.
    cntl%trevers=.FALSE.
    cntl%fmatch=.FALSE.
    ! cntl%cdft OPTIONS
    cntl%cdft=.FALSE.
    cdftlog%newgrad=.FALSE.
    cdftlog%dekk=.FALSE.
    cntl%tsyscomb=.FALSE.
    cntl%tksham=.FALSE.
    cntl%tscombm=.FALSE.
    cntl%tscortho=.TRUE.
    cdftcom%cdft_v(1)=0.0_real_8
    cdftcom%cdft_v(2)=0.0_real_8
    cdftvgrs%vgstep=10.0_real_8
    cdftvgrs%maxvmov=0.1_real_8
    cdftcom%nother=0.0_real_8
    cdftlog%tczones=.FALSE.
    glepar%gle_mode=0              ! GLE
    glepar%gle_omega=0.0_real_8
    glepar%gle_ns=0
    glepar%gle_com=1
    tglepc=.TRUE. 
    cnti%iprng=123456
    czones(1,1)=10.0_real_8
    czones(1,2)=1.0e-3_real_8 ! outmost zone sloppy convergence
    czones(2,1)=0.3_real_8
    czones(2,2)=1.0e-4_real_8 ! medium zone average convergence
    czones(3,1)=0.1_real_8
    czones(3,2)=1.0e-5_real_8 ! innermost zone good convergence
    cntl%cdft_weight=.FALSE.
    cdftcom%wslice=0.5_real_8
    cdftci%wstep=1
    wg_n=0
    cdftlog%reswf=.FALSE.
    cdftlog%thda=.FALSE.
    cdftlog%thdaproj=.FALSE.
    cdftlog%ccor=.TRUE. ! cutoff correction per default on
    cdftcom%vccon=1.0e-5_real_8
    cdftcom%vcconu=1.0e-5_real_8
    cdftlog%tpred=.FALSE.
    cdftpred%predord=5
    cdftlog%tphio=.FALSE.
    cdftlog%tspinc=.FALSE.
    projlog%projnow=.FALSE.
    cdftlog%recw=.TRUE. ! cntl%cdft recalculate weights
    cntl%tipao=.FALSE.
    cdftlog%tvmirr=.FALSE.
    cm_dir=0
    sccomm%tsysk=-1
    sccomm%tsysl=-1
    ! Restart options
    restart1%rwf=.FALSE.
    restart1%rco=.FALSE.
    restart1%rcell=.FALSE.
    restart1%rac=.FALSE.
    restart1%rhe=.FALSE.
    restart1%rnoe=.FALSE.
    restart1%rnop=.FALSE.
    restart1%rnoc=.FALSE.
    restart1%rvel=.FALSE.
    restart1%rgeo=.FALSE.
    restart1%rlate=.FALSE.
    restart1%rvib=.FALSE.
    restart1%rpot=.FALSE.
    restart1%rrho =.FALSE.
    restart1%rkpt =.FALSE.
    restart1%rlr=.FALSE.
    restart1%rphes=.FALSE.
    restart1%rlscl=.FALSE.
    restart1%rdtl =.FALSE.
    restart1%rcon =.FALSE.
    restart1%rnon =.FALSE.
    restart1%rxtp =.FALSE.
    restart1%rprng =.FALSE.
    restart1%rgle =.FALSE.
    ! Store options
    store1%swf=.TRUE.
    store1%srho=.TRUE.
    store1%spot=.FALSE.
    ! Store options in another files
    rout1%rhoout=.FALSE.
    rout1%nrhoout=0
    rout1%rout=.TRUE.
    rout1%teband=.FALSE.
    cntl%tlowd=.FALSE.
    cntl%tlowdmat=.FALSE.
    cntl%quenchp=.FALSE.
    cntl%quenche=.FALSE.
    cntl%quenchc=.FALSE.
    cntl%quenchb=.FALSE.
    cntl%tberp=.FALSE.
    cntl%tbere=.FALSE.
    cntl%tberc=.FALSE.
    cntl%trane=.FALSE.
    cntl%tranp=.FALSE.
    cntl%tranc=.FALSE.
    cntl%thead=.TRUE.
    cntl%tc=.FALSE.
    cntl%tcp=.FALSE.
    cntl%annei=.FALSE.
    cntl%annee=.FALSE.
    cntl%annec=.FALSE.
    cntl%dampi=.FALSE.
    cntl%dampe=.FALSE.
    cntl%dampc=.FALSE.
    cntl%stopot=.FALSE.
    cntl%tchngw=.FALSE.
    cntl%tchngr=.FALSE.
    cntl%tsde=.FALSE.
    cntl%tsdp=.FALSE.
    sdpl%tcgp=.FALSE.
    sdpl%tsdpline=.FALSE.
    cntl%tsdc=.FALSE.
    cntl%diis=.FALSE.
    cntl%pcg=.FALSE.
    cntl%pcgmin=.FALSE.
    cntl%gdiis=.FALSE.
    cntl%bfgs=.FALSE.
    cntl%lbfgs=.FALSE.
    cntl%prfo=.FALSE.
    cntl%rfo=.FALSE.
    cntl%caldip=.FALSE.
    cntl%tnosee=.FALSE.
    cntl%tnosep=.FALSE.
    tcafes=.FALSE.
    nosl%tultra=.FALSE.
    loct%tloct=.FALSE.
    nosl%tmnose=.FALSE.
    tnosepc=.TRUE. 
    cntl%tnosec=.FALSE.
    cntl%timing=.TRUE.
    cntl%rcomp=.FALSE.
    cntl%wcomp=.FALSE.
#if defined (__NEC) 
    cntl%bigmem=.TRUE.
#else
    cntl%bigmem=.FALSE.
#endif
    cntl%fcnstr=.FALSE.
    cntl%krwfn=.FALSE.
    cntl%posver=.FALSE.
    cntl%initcm=.FALSE.
    cntl%finalcm=.FALSE.
    bond=.FALSE.
    angle=.FALSE.
    dihedral=.FALSE.
    cntl%tpres=.FALSE.
    isos1%tisos=.FALSE.
    isos1%tcent=.TRUE.
    qspl1%qspline=.FALSE.
    qspl1%initspl=.TRUE.
    cntl%simul=.FALSE.
    cntl%tmemchk=.FALSE.
    cntl%tpenn=.FALSE.
    cntl%tdiag=.FALSE.
    cntl%tdiagopt=.FALSE.
    cntl%tfint=.FALSE.
    fint1%tbogo=.FALSE.
    fint1%tfral=.FALSE.
    andr2%trand=.FALSE.
    cntl%tdavi=.FALSE.
    cntl%tlanc=.FALSE.
    cntl%tfrho_upw=.FALSE.
    fint1%ttrot=.FALSE.
    tmovr=.FALSE.
    cntl%textrap=.FALSE.
    cntl%tstrxtp=.FALSE.
    cntl%taspc=.TRUE.
    cntl%tdipd=.FALSE.
    cntl%tinter=.FALSE.
    elfcb%telf=.FALSE.
    cntl%tepot=.FALSE.
    cntl%tcc=.FALSE.
    cntl%texpot=.FALSE.
    cntl%texadd=.FALSE.
    cntl%tpath=.FALSE.
    cntl%tpimd=.FALSE.
    cntl%tpmin=.FALSE.
    cntl%trescale=.FALSE.
    comvl%tsubcom=.FALSE.
    comvl%tsubrot=.FALSE.
    cntl%tfrsblk=.FALSE.
    broy1%tgmix=.FALSE.
    broy1%tgbroy=.FALSE.
    broy1%tadbroy=.FALSE.
    cntl%tprec=.FALSE.
    wannl%twann=.FALSE.
    wannl%twpri=.FALSE.
    wannl%twdos=.FALSE.
    wannl%twmol=.FALSE.
    wannl%tsdens=.FALSE.
    cntl%tfdist=.FALSE.
    cntl%tsymrho=.TRUE.
    cntl%tresponse=.FALSE.
    cntl%thard=.FALSE.
    cntl%tfusi=.FALSE.
    cntl%tmerge=.FALSE.
    cntl%tnacvs=.FALSE.
    locpot2%tlpot=.FALSE.
    rout1%vout=.FALSE.
    rout1%acout=.FALSE.
    rout1%xtin=.FALSE.
    rout1%xtout=.FALSE.
    rout1%xgout=.FALSE.
    rout1%dcout=.FALSE.
    cprint%twritebintrajectory=.FALSE.
    cprint%minwriteatom=1
    cprint%maxwriteatom=1000000
    cprint%twriteforcetrajectory=.FALSE.
    ! Debugging
    cntl%tdebfor=.FALSE.
    store1%tdebio=.FALSE.
    store1%tdebacc=.FALSE.
    ! Tracing
    cp_trace%ttrace=.FALSE.
    cp_trace%ttrace_master_only=.FALSE.
    tname%trace_nbr_procedure=0
    tname%trace_procedure=REPEAT('X',LEN(tname%trace_procedure))
    tname%trace_max_depth=HUGE(0)
    tname%trace_max_calls=HUGE(0)
    ! Time-Dependent DFT
    cntl%tddft=.FALSE.
    cntl%tnogeocheck=.FALSE.
    cntl%tssel=.FALSE.
    ! EHR[
    cntl%cmplx_wf=.FALSE.
    cntl%tmdeh=.FALSE.
    cntl%tpspec=.FALSE.
    cntl%tpdist=.FALSE.
    cntl%tgaugep=.FALSE.
    cntl%tgaugef=.FALSE.
    cntl%cheby=.FALSE.
    cntl%cayley=.FALSE.
    cntl%ruku=.FALSE.
    ! EHR]     
    ! ROKS
    lspin2%tlse=.FALSE.
    lspin2%troks=.FALSE.
    lspin2%tros=.FALSE.
    lspin1%lsea=2._real_8
    lspin1%lseb=-1._real_8
    lspin3%mgab(1)=0.5_real_8
    lspin3%mgab(2)=0.5_real_8
    lspin3%mgab(3)=0.5_real_8
    lspin3%mgab(4)=0.5_real_8
    lspin3%mgab(5)=0.5_real_8
    lspin3%mgab(6)=0.5_real_8
    lspin3%mgab(7)=0.5_real_8
    lspin3%mgab(8)=0.5_real_8
    lspin3%mgab(9)=1.5_real_8
    lspin3%mgab(10)=-0.5_real_8
    lspin3%mgab(11)=0.5_real_8
    lspin3%mgab(12)=0.5_real_8
    ! Wfn Localization in G space
    glocal%tgloc=.FALSE.
    glocal%tglocprint=.FALSE.
    glocal%tg_kick = .FALSE.
    glocal%tg_complex = .TRUE.
    glocal%tg_real = .FALSE.
    glocal%tg_tune_st = .FALSE.
    glocal%tg_antisymm = .FALSE.
    glocal%tg_penalty = .FALSE.
    glocal%tg_linesearch = .FALSE.
    glocal%tg_read_matrix = .FALSE.
    glocal%tglocrealp = .FALSE.
    glocal%tgwannier = .FALSE.
    glocal%torbdist  = .FALSE.
    cntl%ttau      = .FALSE.
    cntl%tpotential= .FALSE.
    cntl%tr4a2a = .FALSE.
    ! Initialised here because Kohn-Sham energies used it.
    tkpts%tknoswap=.FALSE.
    inr_logical%tmixsd=.FALSE.
    inr_logical%tmixgdiis=.FALSE.
    cntl%tinr=.FALSE.
    inr_logical%inr_prec=.FALSE.
    inr_logical%inr_step=.FALSE.
    inr_logical%inr_verbose=.FALSE.
    inr_logical%inr_cont=.FALSE.
    cntl%twserial=.FALSE.
    lmeta%lsadpnt=.FALSE.   ! cmbike
    cntl%tdmal=.FALSE.
    cntl%tortho_new=.FALSE.        ! this set the new orthogonalization of the wavefunction
    ! C.Bekas and A.Curioni Comp. Phys. Comm. (2010)
    prop5%teprefg = .FALSE.
    prop5%tefg = .FALSE.
    cntl%tpcgfi = .FALSE.
    cntl%tpcgfic = .FALSE.
    cntl%gshellout = .FALSE.
    vdwl%dcacp = .FALSE.
    vdwl%vdwc  = .FALSE.
    vdwl%vdwd  = .FALSE.
    cntl%tcart=.FALSE.
    cntl%tshock = .FALSE.
    cntl%tsic=.FALSE.      ! cmb_ssic
    cp_cuda_env%use_blas = .FALSE.
    cp_cuda_env%use_fft = .FALSE.
    cp_cuda_env%fft_n_streams_per_device = 1
    cp_cuda_env%fft_n_devices_per_task = 1
    cp_cuda_env%blas_n_streams_per_device = 1
    cp_cuda_env%blas_n_devices_per_task = 1
    ! ==--------------------------------------------------------------==
    fo_info%fpath=' '
    CALL m_getenv('CPMD_FILEPATH',fo_info%fpath)
    IF ((INDEX(fo_info%fpath,' ')-1).EQ.0) fo_info%fpath='./'
    ! ==--------------------------------------------------------------==
    cnti%nomore=10000
    cnti%nomore_iter=10000
    cnti%nps=1
    cprint%tprint=.FALSE.
    cprint%iprint_step=0
    DO i=1,iprint_max
       cprint%iprint(i)=0
    ENDDO
    store1%istore=0
    store1%isctore=0
    cnti%ntgaus=1
    comvl%ncomv=10
    comvl%nrotv=10
    ! Timesteps for electrons and ions
    cntr%wnosp0 = 0.0_real_8
    cntr%wnose0 = 0.0_real_8
    cntr%wnosc0 = 0.0_real_8
    cntr%smf = 0.0_real_8
    cntr%tolini = 0.0_real_8
    cntr%toll = 0.0_real_8
    cntr%tolp = 0.0_real_8
    cntr%delt_ions=5._real_8
    cntr%delt_elec=5._real_8
    cntr%emass=400._real_8
    cntr%cmass=-1.0_real_8
    cntr%hthrs=0.5_real_8
    cntr%tolog=0._real_8
    cntr%tolad=0._real_8
    cntr%tolene=0._real_8
    cntr%tolfor=3._real_8
    cntr%tolrhofix=0.0_real_8
    cntr%toldetot=0._real_8
    cntr%tolkinc=-1.0_real_8
    cnti%nstcnv=0
    cntr%tolng=5.e-4_real_8
    cntr%tempw=0.0_real_8
    cntr%trampt=0.0_real_8
    cntr%trampr=0.0_real_8
    cnti%nperid=0
    cnti%maxit=30
    cnti%mdiis=10
    cnti%maxldiis=0
    cnti%minldiis=0
    cnti%mdiis_fr=4
    cnti%mgdiis=5
    cnti%nreset=0
    cntr%epsog=1.e-6_real_8
    cntr%amprc=0.001_real_8
    cntr%ampre=0.005_real_8
    cntr%amprp=0.05_real_8
    andr2%amprd=0.001_real_8
    cnti%insys=5
    cnti%npara=0
    cnti%nchp=3
    cnti%nchs=3
    cnti%nchb=3
    cnti%ncalls0=7
    cnti%nit0=1
    cntr%nedof0=6._real_8
    cnti%nkssta=0
    cnti%ndavv=0
    cnti%imovie=0
    rout1%mout=.FALSE.
    cnti%iproj=2
    cntr%epsdav=1.e-5_real_8
    cnti%ivdbrd=0
    cnti%ivdbwr=0
    ! Finite Difference options
    cntr%fdiff=1.e-2_real_8
    tref_fdiff=.FALSE.
    coord_fdiff(1) = 0._real_8
    coord_fdiff(2) = 0._real_8
    coord_fdiff(3) = 0._real_8
    r_fdiff = 0._real_8
    cntr%memsize=-1._real_8
    cnti%inwfun=2
    cnti%nrestf=1
    cnti%ntrans=1
    numpr=0
    cnti%ntraj=1
    cnti%ngxyz=1
    cnti%nvib=1
    cnti%npres=1
    group%nogrp=1
    cnti%nsplp=5000
    cnti%nsorder=0
    cnti%wcompb=1
    cnti%rcompb=1
    nacc=nacx
    CALL zeroing(acc)!,nacx)
    CALL zeroing(ibench)!,nbentr)
    andr2%ntabmix=1
    andr2%andrmix=0.2_real_8
    broy1%broymix=0.15_real_8
    broy1%w02broy=0.01_real_8
    broy1%ecutbroy=-1._real_8
    broy1%kermix=0.0_real_8
    andr2%alxmix=0.9_real_8
    fint1%betael=-1._real_8
    fint1%betap=0._real_8
    cntr%b2limit=1.e-18_real_8
    cnti%n_fries=-1
    cnti%nkry_max=6
    cnti%nkry_block=-1
    fint4%ntabtrot=10
    nssel=0
    cnti%nskip=0
    cnti%nsample=1
    cnti%mextra=-1
    cnti%naspc=-1
    cnti%nstblk=0! 16
    cnti%disortho_bsize=0
    ! .. 
    fint4%denstrot(1) = 1._real_8
    fint4%b2trot(1) = 1.e-7_real_8
    fint4%denstrot(2) = 0.1_real_8
    fint4%b2trot(2) = 1.e-8_real_8
    fint4%denstrot(3) = 0.07_real_8
    fint4%b2trot(3) = 1.e-10_real_8
    fint4%denstrot(4) = 0.05_real_8
    fint4%b2trot(4) = 1.e-11_real_8
    fint4%denstrot(5) = 0.02_real_8
    fint4%b2trot(5) = 1.e-12_real_8
    fint4%denstrot(6) = 0.01_real_8
    fint4%b2trot(6) = 1.e-13_real_8
    fint4%denstrot(7) = 0.007_real_8
    fint4%b2trot(7) = 1.e-14_real_8
    fint4%denstrot(8) = 0.004_real_8
    fint4%b2trot(8) = 1.e-15_real_8
    fint4%denstrot(9) = 0.002_real_8
    fint4%b2trot(9) = 1.e-16_real_8
    fint4%denstrot(10) = 0.0012_real_8
    fint4%b2trot(10) = 1.e-17_real_8
    ! ..
    cntr%dshift=0._real_8
    fint5%ntabbetap=1
    andr3%nrdiis=0
    andr3%nrdiismax=0
    andr3%ntabnrdiis=1
    dmovrmix=0.5_real_8
    cntr%tolcg=1._real_8
    cnti%npdip=1
    cnti%iftype=0
    elfcb%elfcut=0.001_real_8
    elfcb%elfeps=2.87e-5_real_8
    cntr%tolc=0._real_8
    cntr%ekinhr=0._real_8
    cntr%taubp=1._real_8
    cntr%taube=1._real_8
    cntr%taubc=1._real_8
    qsrang=1._real_8
    cnti%icmet=3
    broy1%nfrbroy=0
    broy1%ibreset=8
    broy1%ngbroy=0
    wanni%w_maxs=2000
    wanni%w_type=1
    wanni%w_opt=2
    wanni%ns_wann=1
    wanni%sw_orb=0
    wanni%sw_first=-1
    wanni%sw_last=-1
    wanni%sw_all=0
    wannr%sw_spread=0.0_real_8
    wan05%loc_npgrp = parai%cp_nproc
    wan05%loc_relocalize_every=1
    wan05%loc_relocalize=.FALSE.
    wan05%loc_recompute_dipole_matrices_every=1
    wan05%loc_relocalize_in_scf=.FALSE.
    wannr%w_step=0.1_real_8
    wannr%w_eps=1.e-8_real_8
    wannr%w_ran=0.0_real_8
    wannr%w_ref(1)=-1000000.0_real_8
    wannr%w_ref(2)=-1000000.0_real_8
    wannr%w_ref(3)=-1000000.0_real_8
    CALL zeroing(restf%nstepwr)!,maxrf)
    restf%nstepwr(1)=-1
    restf%nrcopy=0
    cntr%anneri=1._real_8
    cntr%annere=1._real_8
    cntr%annerc=1._real_8
    cntr%dampgi=0._real_8
    cntr%dampge=0._real_8
    cntr%dampgc=0._real_8
    rmixsd=5.e-3_real_8
    inr_integer%nreg=5
    gnx_inr(1)=5.e-3_real_8
    gnx_inr(2)=1.e-3_real_8
    gnx_inr(3)=5.e-4_real_8
    gnx_inr(4)=1.e-4_real_8
    gnx_inr(5)=5.e-5_real_8
    tolx_inr(1)=5.e-3_real_8
    tolx_inr(2)=1.e-3_real_8
    tolx_inr(3)=5.e-4_real_8
    tolx_inr(4)=1.e-4_real_8
    tolx_inr(5)=5.e-5_real_8
    tol_inr=tolx_inr(1)
    gloci%gloc_maxs=20000
    gloci%gloc_type=2
    gloci%gloc_opt=2
    gloci%gloc_const = 1
    gloci%gloc_orb=0
    gloci%gloc_first=-1
    gloci%gloc_last=-1
    gloci%gloc_all=0
    glocr%gloc_step=0.1_real_8
    glocr%gloc_eps=1.e-12_real_8
    glocr%gloc_ran=0.0_real_8
    gloci%gloc_init = 0
    gloci%gloc_maxit=20
    glocr%gepslag   = 0.000001_real_8
    glocr%wan_weight  = 0.0_real_8
    glocr%g2g_weight  = 1.0_real_8
    ! mb_ssic... SIC defaults
    cntr%asic=0.2_real_8
    cntr%bsic=0.0_real_8
    ! EXACT FACTORIZATION
    tshl%txfmqc=.FALSE.
    xfmqc%n_xftraj=1
    ! IPRINT
    CALL zeroing(cprint%iprint)!,iprint_max)
    ! 
    nac=-1
    nbc=-1
    ! ==--------------------------------------------------------------==
    wannc%nwanopt=0
    ! ==--------------------------------------------------------------==
    isrlength = 0
    ! ==--------------------------------------------------------------==
    iface1%intread = .FALSE.
    iface1%intwrite = .FALSE.
    intfn = "interface.bin"
    ! ==--------------------------------------------------------------==
    cntl%mimic = .FALSE.
    cntl%new_constraints = .FALSE.
    cntl%pbicgstab = .FALSE.
    cnti%shake_maxstep = 5000
    cnti%shake_cg_iter = 50
    cntl%anneal_dual = .FALSE.
    cntr%anneal_factors = 1.0_real_8
    RETURN
  END SUBROUTINE control_def
  ! ==================================================================

END MODULE control_def_utils
