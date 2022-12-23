MODULE mm_mdmain_utils
  USE anneal_utils,                    ONLY: anneal,&
                                             berendsen,&
                                             dampdyn,&
                                             tempramp
  USE bsym,                            ONLY: autocm,&
                                             bsclcs,&
                                             bsfac,&
                                             cnstwgt,&
                                             ekinc_bs,&
                                             ekinc_hs,&
                                             rtsasb
  USE cnst,                            ONLY: factem
  USE cnst_dyn,                        ONLY: ekincv,&
                                             fhills,&
                                             lmeta,&
                                             ltcglobal,&
                                             ncolvar,&
                                             rmeta,&
                                             vharm
  USE comvel_utils,                    ONLY: comvel
  USE comvelmod,                       ONLY: comvl
  USE coor,                            ONLY: fion,&
                                             tau0,&
                                             taup,&
                                             velp
  USE copot_utils,                     ONLY: copot,&
                                             give_scr_copot
  USE cotr,                            ONLY: cotc0,&
                                             lskcor,&
                                             ntcnst
  USE csize_utils,                     ONLY: csize
  USE ddipo_utils,                     ONLY: ddipo,&
                                             give_scr_ddipo
  USE deort_utils,                     ONLY: deort,&
                                             give_scr_deort
  USE detdof_utils,                    ONLY: detdof
  USE dispp_utils,                     ONLY: dispp
  USE dynit_utils,                     ONLY: dynit
  USE efld,                            ONLY: extf,&
                                             textfld
  USE ekinpp_utils,                    ONLY: ekinpp
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_app,&
                                             fo_verb
  USE finalp_utils,                    ONLY: finalp
  USE forcedr_driver,                  ONLY: forcedr
  USE forcedr_utils,                   ONLY: give_scr_forcedr
  USE forcep_utils,                    ONLY: rhoe_psi_size
  USE freqs_utils,                     ONLY: freqs
  USE geq0mod,                         ONLY: geq0
  USE gsize_utils,                     ONLY: gsize
  USE hubbardu,                        ONLY: hubbu
  USE initrun_driver,                  ONLY: initrun
  USE initrun_utils,                   ONLY: give_scr_initrun
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE localize_utils,                  ONLY: localize2
  USE machine,                         ONLY: m_flush,&
                                             m_walltime
  USE meta_colvar_inp_utils,           ONLY: colvar_structure
  USE meta_colvar_utils,               ONLY: meta_colvar
  USE meta_exl_mult_utils,             ONLY: meta_ext_mul
  USE meta_exlagr_methods,             ONLY: give_scr_meta_extlagr,&
                                             meta_extlagr
  USE meta_exlagr_utils,               ONLY: ekincv_global
  USE mm_dim_utils,                    ONLY: mm_dim
  USE mm_dimmod,                       ONLY: &
       clsaabox, gratom, mm_go_mm, mm_go_qm, mm_revert, mmdim, nam, naq, &
       solsolv
  USE mm_input,                        ONLY: clc,&
                                             g96_vel,&
                                             lqmmm,&
                                             rtr_l,&
                                             wp_i
  USE mm_parallel,                     ONLY: gparal
  USE mm_qmmm_forcedr_bs_utils,        ONLY: mm_qmmm_forcedr_bs
  USE mm_qmmm_forcedr_utils,           ONLY: mm_qmmm_forcedr
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum,&
                                             mp_sync
  USE nlcc,                            ONLY: corel
  USE norm,                            ONLY: cnorm,&
                                             gemax,&
                                             gnmax,&
                                             gnorm
  USE nose,                            ONLY: glib
  USE noseinit_utils,                  ONLY: noseinit
  USE noseng_utils,                    ONLY: noseng
  USE nosepa_utils,                    ONLY: nosepa
  USE noseup_utils,                    ONLY: noseup
  USE nospinit_utils,                  ONLY: nospinit
  USE ortho_utils,                     ONLY: give_scr_ortho
  USE parac,                           ONLY: parai,&
                                             paral
  USE phfac_utils,                     ONLY: phfac
  USE posupa_utils,                    ONLY: give_scr_posupa,&
                                             posupa
  USE posupi_utils,                    ONLY: posupi
  USE printave_utils,                  ONLY: pacca
  USE printp_utils,                    ONLY: printp,&
                                             printp2
  USE prmem_utils,                     ONLY: prmem
  USE proja_utils,                     ONLY: proja
  USE pslo,                            ONLY: pslo_com
  USE puttau_utils,                    ONLY: taucl
  USE quenbo_utils,                    ONLY: give_scr_quenbo,&
                                             quenbo
  USE rattle_utils,                    ONLY: rattle
  USE rekine_utils,                    ONLY: rekine
  USE resetac_utils,                   ONLY: resetac
  USE reshaper,                        ONLY: reshape_inplace
  USE rhopri_utils,                    ONLY: give_scr_rhopri,&
                                             rhopri
  USE rinvel_utils,                    ONLY: rinvel,&
                                             rvscal
  USE rmas,                            ONLY: rmass
  USE ropt,                            ONLY: infi,&
                                             iteropt,&
                                             ropt_mod
  USE rortv_utils,                     ONLY: give_scr_rortv,&
                                             rortv
  USE rotvel_utils,                    ONLY: rotvel
  USE rscve_utils,                     ONLY: rscve
  USE rscvp_utils,                     ONLY: rscvp
  USE sample_utils,                    ONLY: sample_go,&
                                             sample_wait
  USE setbsstate_utils,                ONLY: setbsstate
  USE setirec_utils,                   ONLY: write_irec
  USE shake_utils,                     ONLY: cpmdshake
  USE shop_ekinqm,                     ONLY: ekinqmsh
  USE soft,                            ONLY: soft_com
  USE spin,                            ONLY: spin_mod
  USE store_types,                     ONLY: &
       cprint, irec_ac, irec_noe, irec_nop1, irec_nop2, irec_nop3, irec_nop4, &
       irec_vel, irec_wf, restart1, rout1
  USE system,                          ONLY: &
       cnti, cntl, cntr, fpar, maxsys, nacc, ncpw, nkpt, restf
  USE testex_utils,                    ONLY: testex
  USE teststore_utils,                 ONLY: teststore
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE totstr_utils,                    ONLY: totstr
  USE tpar,                            ONLY: dt_ions
  USE tst2min_utils,                   ONLY: tst2min
  USE utils,                           ONLY: zclean,&
                                             zclean_k
  USE vdwcmod,                         ONLY: vdwl,&
                                             vdwwfl
  USE velupa_utils,                    ONLY: velupa
  USE velupi_utils,                    ONLY: velupi
  USE wann,                            ONLY: wannl
  USE wannier_print_utils,             ONLY: wannier_print
  USE wc_dos_utils,                    ONLY: wc_dos
  USE wrener_utils,                    ONLY: wrprint_md
  USE wrgeo_utils,                     ONLY: wrgeof
  USE wv30_utils,                      ONLY: zhwwf
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mm_mdmain
  PUBLIC :: give_scr_mm_mdmain
  PUBLIC :: mm_localt
  !public :: velocities

CONTAINS

  ! ==================================================================
  SUBROUTINE mm_mdmain(c0,cm,c2,sc0,gamx,gamy)
    ! ==--------------------------------------------------------------==

    ! Qmmm

    ! Qmmm locals
    COMPLEX(real_8) :: c0(ncpw%ngw,crge%n,nkpt%nkpnt), cm(ncpw%ngw,crge%n,*), &
      c2(ncpw%ngw,crge%n,*), sc0(ncpw%ngw,crge%n,*)
    REAL(real_8)                             :: gamx(crge%n*crge%n,*), &
                                                gamy(crge%n*crge%n,*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'mm_mdmain'

    CHARACTER(len=10)                        :: prch
    CHARACTER(len=100)                       :: filebs, filen
    CHARACTER(len=30)                        :: tag
    COMPLEX(real_8), ALLOCATABLE             :: psi(:,:)
    INTEGER :: i, ia, ierr, il_psi_1d, il_psi_2d, il_rhoe_1d, il_rhoe_2d, in, &
      irec(100), is, isub, itemp, j, k, lscr, nstate, out
    LOGICAL                                  :: ferror, lmetares, lquench, &
                                                oldstatus, resetcv, &
                                                statusdummy
    REAL(real_8) :: couplj, disa, dummy, econs, eham, ek_cv, ekin1, ekin2, &
      ekinc, ekincp, ekinh1, ekinh2, ekinp, enose, enosp, etotbs, etoths, ff, &
      jwn, lmio(3), scalbs, scalhs, spd_bs, spd_hs, spda_bs, spda_hs, tcpu, &
      temp1, temp2, tempp, time1, time2, vcmio(4)
    REAL(real_8), ALLOCATABLE :: center(:,:), eigm(:,:), eigv(:,:), &
      rhoe(:,:), scr(:), taui(:,:,:), tauio(:,:), taur(:,:,:), tsrc(:,:,:)
    REAL(real_8), DIMENSION(:, :, :), &
      POINTER                                :: scr_correct_shape

! Variables
! META DYNAMICS
! ==================================================================

    NULLIFY( scr_correct_shape )
#if defined (__GROMOS)
    CALL tiset(procedureN,isub)
    time1 =m_walltime()
    IF (paral%qmnode) THEN
       IF (textfld)THEN
          ALLOCATE(extf(fpar%kr1*fpar%kr2s*fpar%kr3s),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          CALL zeroing(extf)!,kr1*kr2s*kr3s)
       ENDIF
       nstate=crge%n
       nacc = 22
       ropt_mod%modens=.FALSE.
       ropt_mod%engpri=.FALSE.
       ropt_mod%calste=cntl%tpres
       IF (pslo_com%tivan) cntl%nonort=.TRUE.
       IF (cntl%nonort) THEN
          ! TODO why these arrays are allocated if cntl%nonort(=TIVAN)
          ! but this is not checked in children subroutines?
          ALLOCATE(eigv(crge%n,bsfac*nstate/crge%n),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(eigm(crge%n**2,bsfac*nstate**2/crge%n**2),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ELSE
          ! Avoid 'not allocated' runtime error
          ALLOCATE(eigv(1,1),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
          ALLOCATE(eigm(1,1),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
       ENDIF
       ! SCALING FACTORS FOR ELEC. KINETIC ENERGIES FOR BROKEN SYMMETRY
       IF (cntl%bsymm)THEN
          scalhs = -1.0_real_8*cnstwgt
          scalbs =  1.0_real_8+cnstwgt
       ENDIF
       ! 
       IF (cntl%tharm.AND.paral%parent) THEN
          ff=crge%f(1,1)
          DO i=1,nstate
             IF (ff.NE.crge%f(i,1)) THEN
                IF (paral%io_parent)&
                     WRITE(6,'(A,A)') ' HARMONIC REFERENCE SYSTEM INTEGRATION',&
                     ' ONLY POSSIBLE WITH EQUAL OCCUPATION NUMBERS'
                CALL stopgm('MM_MDMAIN','FATAL ERROR',& 
                     __LINE__,__FILE__)
             ENDIF
          ENDDO
       ENDIF
       ! ==--------------------------------------------------------------==
       ! SCR ALLOCATION AND PARTITION (SCRATCH ARRAY).
       CALL rhoe_psi_size(il_rhoe_1d=il_rhoe_1d, il_rhoe_2d=il_rhoe_2d,&
            il_psi_1d=il_psi_1d,il_psi_2d=il_psi_2d)
       ALLOCATE(rhoe(il_rhoe_1d, il_rhoe_2d),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(psi(il_psi_1d,il_psi_2d),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL give_scr_mm_mdmain(lscr,tag)
       ALLOCATE(scr(lscr),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF ! qmnode
    ! ==--------------------------------------------------------------==
99999 IF (cntl%tsampl) THEN
       CALL sample_wait
       IF (cnti%nomore.LT.0) GOTO 10000
    ENDIF
    restf%nfnow=1
    ! ==--------------------------------------------------------------==
    ! == INITIALIZATION                                               ==
    ! ==--------------------------------------------------------------==

    ! we may be able to save some memory by allocating more carefully
    ! on every node a different size. Could be a bad idea.

    CALL mm_dim(mm_go_mm,oldstatus)
    ALLOCATE(fion(3,maxsys%nax,maxsys%nsx),tsrc(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(taui(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(taur(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(taup(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(tauio(3,ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    IF (lmeta%lcolvardyn) THEN
       ALLOCATE(fhills(3,maxsys%nax,maxsys%nsx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(fhills)!,3*maxsys%nax*maxsys%nsx)
    ENDIF
    ! Initialize logical variable for Metadynamics
    lmetares=.FALSE.

    CALL zeroing(fion)!,3*maxsys%nax*maxsys%nsx)
    CALL zeroing(taui)!,3*maxsys%nax*maxsys%nsx)
    CALL zeroing(taup)!,3*maxsys%nax*maxsys%nsx)

    ! show sizes of the individual MPI-threads.
    ! use barriers so that the output is not garbled.
    CALL mp_sync(parai%qmmmgrp)
    IF (lqmmm%qmmm_verbose) THEN
       DO k=0,parai%nproc
          IF (gparal%mmparent) THEN
             IF (paral%io_parent)&
                  WRITE(prch,'(a10)') ' MMPARENT '
          ELSE
             IF (paral%io_parent)&
                  WRITE(prch,'(a5,i5)')' QM-N',parai%me
          ENDIF
          IF (k.EQ.parai%me) CALL prmem(prch)
          CALL mp_sync(parai%qmmmgrp)
       ENDDO
    ENDIF

    ! format for for BS messages
333 FORMAT (/,t2,a,/)
    IF (paral%qmnode ) THEN
       ! if (qmnode .and. .not.classical) then 
       IF (cntl%bsymm)THEN
          bsclcs=1
          IF (paral%io_parent)&
               WRITE(6,333) 'BROKEN SYMMETRY STATE IS INITIALIZING'
       ENDIF
       CALL initrun(irec,c0,cm,sc0,rhoe,psi,eigv)
       IF (cntl%bsymm)THEN
          bsclcs=2
          IF (paral%io_parent)&
               WRITE(6,333) 'HIGH SPIN STATE IS INITIALIZING'
          CALL setbsstate
          CALL initrun(irec,c0(:,:,2:),cm(1,1,2),sc0(1,1,2),rhoe,&
               psi,eigv(1,2))
          IF (paral%io_parent)&
               WRITE(6,333) 'BROKEN SYMMETRY INIT SUCCESSFUL'
       ENDIF
    ENDIF
    IF (paral%qmnode ) THEN
       ! ==--------------------------------------------------------------==
       ! TIME STEP FUNCTIONS
       CALL dynit(ekincp,ekin1,ekin2,temp1,temp2,ekinh1,ekinh2)
       ekinc=0.0_real_8  ! McB: cf. META_EXT..()
       ekinp=0.0_real_8  ! McB: cf. META_EXT..()
       ! PARAMETERS FOR THE NOSE-HOOVER THERMOSTATS
       IF ((cntl%tnosee.OR.cntl%tnosep).AND.paral%parent) CALL nosepa(1,1)

       ! Dont symmetrize density 
       cntl%tsymrho=.FALSE.
       CALL mp_bcast(taup,SIZE(taup),parai%source,parai%allgrp)

       irec(irec_wf)=1

       CALL dcopy(3*maxsys%nax*maxsys%nsx,taup,1,taui,1)! cmb-bugfix
    ENDIF
    IF (cntl%quenchb) THEN
       CALL mm_dim(mm_go_mm,statusdummy)
       CALL zeroing(cm(:,:,1:bsfac))!,ngw*nstate*bsfac)
       IF (.NOT.lqmmm%qmmm_reflex) THEN
          CALL mm_translate_qmmm(tau0,c0,cm,nstate)
          IF (cntl%bsymm) THEN
             CALL mm_translate_qmmm(tau0,c0(1,1,2),cm(1,1,2),nstate)
          ENDIF
       ENDIF
       CALL mm_dim(mm_go_qm,statusdummy)

       IF (cntl%bsymm)THEN
          bsclcs=1
          CALL setbsstate
          IF (paral%io_parent)&
               WRITE(6,333) 'QUENCHING BROKEN SYMMETRY STATE'
       ENDIF
       CALL quenbo(c0(:,:,1),c2,sc0,taur,rhoe,psi)
       IF (cntl%bsymm)THEN
          bsclcs=2
          CALL setbsstate
          IF (paral%io_parent)&
               WRITE(6,333) 'QUENCHING HIGH SPIN STATE'
          CALL quenbo(c0(:,:,2),c2(1,1,2),sc0(1,1,2),taur,&
               rhoe,psi)
       ENDIF
       CALL mm_dim(mm_go_mm,statusdummy)

       IF (paral%qmnode ) THEN
          CALL dcopy(3*maxsys%nax*maxsys%nsx,taup,1,taui,1)
       ENDIF
#if 0
       ! FIXME: AK 2005/06/17. the restart from QUENCH BO is broken.
       IF (paral%qmnode)&
            CALL ZHWWF(2,IREC,C0,CM,NSTATE,EIGV,TAU0,VELP,TAUI,iteropt%nfi)
#endif
       cntl%quenchb=.FALSE.
    ENDIF
    IF (paral%qmnode)THEN
       IF (pslo_com%tivan) THEN
          CALL mm_dim(mm_go_qm,statusdummy)
          IF (cntl%tlsd) THEN
             bsclcs=1
             IF (cntl%bsymm)CALL setbsstate
             CALL deort(ncpw%ngw,spin_mod%nsup,eigm,eigv,c0(1,1,1),sc0(1,1,1))
             CALL deort(ncpw%ngw,spin_mod%nsdown,eigm,eigv,c0(1,spin_mod%nsup+1,1),&
                  sc0(1,spin_mod%nsup+1,1))
             IF (cntl%bsymm)THEN
                bsclcs=2
                CALL setbsstate
                CALL deort(ncpw%ngw,spin_mod%nsup,eigm(1,2),eigv(1,2),c0(1,1,2),&
                     sc0(1,1,2))
                CALL deort(ncpw%ngw,spin_mod%nsdown,eigm(1,2),eigv(1,2),c0(1,spin_mod%nsup+1,2),&
                     sc0(1,spin_mod%nsup+1,2))
             ENDIF
          ELSE
             CALL deort(ncpw%ngw,nstate,eigm,eigv,c0,sc0)
          ENDIF
       ENDIF
       ! INITIALIZE VELOCITIES
       CALL mm_dim(mm_go_mm,statusdummy)
       IF (paral%parent) CALL detdof(tau0,tsrc)

       ! INITIALIZE METADYNAMICS VARIABLES used also for 
       ! SHOOTING from SADDLE POINT with RANDOM VELOCITIES
       IF (paral%parent .AND. (lmeta%lcolvardyn .OR. lmeta%lsadpnt)) THEN
          CALL colvar_structure(tau0,taur)
       ENDIF

       IF ((irec(irec_vel).EQ.0).AND.&
            (.NOT.restart1%rgeo).AND.(.NOT.rtr_l%restart_traj)) THEN
          ener_com%ecnstr = 0.0_real_8
          ener_com%erestr = 0.0_real_8
          CALL rinvel(velp,cm,crge%n)
          IF (cntl%bsymm) CALL zeroing(cm(:,:,1:bsfac))!,ngw*nstate*bsfac)
          IF (paral%parent) THEN
             CALL taucl(velp)
             CALL mm_solv_const('RATTL',dt_ions,tau0,velp,tau0)
             CALL rattle(tau0,velp)
          ENDIF
          ! RESCALE VELOCITIES ONLY IF NOT READ FROM GROMOS.G96 FILE
          IF ( g96_vel%ntx_vel.NE.1 ) CALL rvscal(velp)
       ELSE
          IF (paral%parent) THEN
             CALL taucl(velp)
             CALL mm_solv_const('RATTL',dt_ions,tau0,velp,tau0)
             CALL rattle(tau0,velp)
          ENDIF
          IF (cntl%trescale) CALL rvscal(velp)
       ENDIF
       IF (cntl%quenchp) CALL zeroing(velp)!,3*maxsys%nax*maxsys%nsx)
       IF (cntl%quenche) CALL zeroing(cm(:,:,1:bsfac))!,ngw*nstate*bsfac)
       IF (clc%classical) CALL zeroing(cm(:,:,1:bsfac))!,ngw*nstate*bsfac)
       IF (cntl%trevers) THEN
          ! invert electronic and ionic velocities (useful for path sampling)
          CALL dscal(2*ncpw%ngw*nstate*bsfac,-1._real_8,cm,1)
          CALL dscal(3*maxsys%nax*maxsys%nsx,-1._real_8,velp,1)
       ENDIF

       ! RESET ACCUMULATORS
       IF (paral%parent.AND.irec(irec_ac).EQ.0)&
            CALL resetac(tau0,taui,iteropt%nfi)
       ! INITIALIZE FORCES
       IF (paral%parent) THEN
          IF (paral%io_parent)&
               WRITE(6,'(1X,64("="))')
          IF (paral%io_parent)&
               WRITE(6,'(1X,"==",T25,A,T64,"==")') 'FORCES INITIALIZATION'
          IF (paral%io_parent)&
               WRITE(6,'(1X,64("="))')
       ENDIF
       CALL mm_dim(mm_go_qm,statusdummy)
       IF (tkpts%tkpnt) THEN
          IF (geq0) CALL zclean_k(c0,nstate,ncpw%ngw)
       ELSE
          IF (geq0) CALL zclean(c0,nstate,ncpw%ngw)
       ENDIF
       ! INITIALIZE WF CENTERS & SPREAD
       IF (vdwl%vdwd.AND..NOT.vdwwfl%trwannc) THEN
          CALL localize2(tau0,c0,c2,sc0,nstate)
       ENDIF
       ! Initialize Metadynamics contributions
       IF (lmeta%lcolvardyn .AND. lmeta%lextlagrange) THEN
          CALL mm_dim(mm_go_mm,statusdummy)
          lquench = .FALSE.
          lmetares= .FALSE.
          resetcv = .FALSE.
          ekinc=0._real_8
          disa=0._real_8
          IF (lmeta%tmulti) THEN
             CALL meta_ext_mul(tau0,velp,taur,&
                  lquench,lmetares,resetcv,ekinc,ekinp)
          ELSE

             CALL meta_extlagr(tau0,velp,taur,&
                  lquench,lmetares,resetcv,ekinc,ekinp)
          ENDIF
          CALL mm_dim(mm_revert,statusdummy)
       ENDIF

    ENDIF ! qmnode

    CALL mm_dim(mm_go_mm,statusdummy)
    IF ( gparal%mmparent ) THEN
       CALL mm_write_gromos_coord('CRD_INI.g96',tau0,velp,maxsys%nax,maxsys%nsx)
    ENDIF

    CALL mp_sync(parai%qmmmgrp)
    IF (cntl%bsymm)THEN
       CALL mm_qmmm_forcedr_bs(c0,c2,sc0,rhoe,psi,tau0,fion,&
            eigv,nstate,0,.FALSE.,.TRUE.,&
            spd_bs,spd_hs,spda_bs,spda_hs,etoths,etotbs,.FALSE.)
    ELSE
       CALL mm_qmmm_forcedr(c0,c2,sc0,rhoe,psi,tau0,fion,&
            eigv,nstate,0,.FALSE.,.TRUE.,.FALSE.)
    ENDIF

    IF (paral%qmnode) THEN

       ! Check orthogonality condition for wavefunction velocities
       CALL mm_dim(mm_go_qm,statusdummy)
       CALL freqs(crge%n,.TRUE.)
       IF (cntl%bsymm)THEN
          bsclcs=1
          CALL setbsstate
       ENDIF

       IF (.NOT.clc%classical) THEN
          CALL rortv(c0,cm,c2,sc0,gamy,nstate)
          IF (cntl%bsymm)THEN
             bsclcs=2
             CALL setbsstate
             CALL rortv(c0(1,1,2),cm(1,1,2),c2(1,1,2),sc0(1,1,2),&
                  gamy(1,2),nstate)
          ENDIF
       ENDIF

       CALL mm_dim(mm_go_mm,statusdummy)
       ! Initialize thermostats
       IF (paral%parent) THEN
          itemp=irec(irec_nop1)+irec(irec_nop2)+irec(irec_nop3)&
               +irec(irec_nop4)
          IF (cntl%tnosep .AND. itemp.EQ.0) CALL nospinit(1)
          IF (cntl%tnosee .AND. irec(irec_noe) .EQ.0) THEN
             CALL noseinit(1)
             IF (cntl%bsymm)CALL noseinit(2)
          ENDIF
          IF ( mmdim%natm.LE.500) CALL wrgeof(tau0,fion)
          filen='ENERGIES'
          IF (paral%io_parent)&
               CALL fileopen(3,filen,fo_app+fo_verb,ferror)
          ! OPEN OUTPUT DATA FILE FOR BROKEN SYMMETRY
          IF (cntl%bsymm) THEN
             filebs='BS_ENERG'
             IF (paral%io_parent)&
                  WRITE(6,'(/,2A)')&
                  '  BROKEN SYMMETRY DATA WILL BE WRITTEN IN ',filen
             IF (paral%io_parent)&
                  WRITE(6,'(/,A,/)')&
                  '  DATA IN BS_ENERG:'
             IF (paral%io_parent)&
                  WRITE(6,'(2A,/)')&
                  '  Step #, EKINC_BS, EKINC_HS, ETOT_LS, ETOT_BS,',&
                  ' ETOT_HS, Spin-density-BS, ..-HS, '
             IF (paral%io_parent)&
                  WRITE(6,'(2A,/)')&
                  '  Absolute-spin-density-BS, ..-HS,',&
                  ' Coupling constant J in a.u., J in cm-1 '
             ! 
             IF (paral%io_parent)&
                  CALL fileopen(277,filebs,fo_app+fo_verb,ferror)
          ENDIF
       ENDIF
       IF (paral%parent) THEN
          IF (paral%io_parent)&
               WRITE(6,'(1X,64("="))')
          IF (paral%io_parent)&
               WRITE(6,'(1X,"==",T20,A,T64,"==")')&
               'END OF FORCES INITIALIZATION'
          IF (paral%io_parent)&
               WRITE(6,'(1X,64("="),/)')
       ENDIF
       CALL write_irec(irec)
       ! ==--------------------------------------------------------------==
       ! == END INITIALIZATION                                           ==
       ! ==--------------------------------------------------------------==
       IF (teststore(0).AND.cntl%tsampl)&
            CALL zhwwf(2,irec,c0,cm,nstate,eigv,taup,velp,taui,iteropt%nfi)
       IF (paral%parent) THEN
          time2 =m_walltime()
          tcpu = (time2 - time1)*0.001_real_8
          IF (paral%io_parent)&
               WRITE(6,'(A,T50,F8.2,A8)') ' TIME FOR INITIALIZATION:',&
               tcpu,' SECONDS'
       ENDIF
    ENDIF ! qmnode
    ALLOCATE(center(4,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ==================================================================
    ! ==          THE BASIC LOOP FOR MOLECULAR DYNAMICS               ==
    ! ==                 USING VELOCITY VERLET                        ==
    ! ==================================================================
    DO infi=1,cnti%nomore
       CALL mm_dim(mm_go_mm,statusdummy)
       CALL mp_sync(parai%qmmmgrp)
       IF (paral%qmnode) THEN
          time1=m_walltime()
          iteropt%nfi=iteropt%nfi+1
          comvl%subcom=comvl%tsubcom.AND.MOD(iteropt%nfi-1,comvl%ncomv).EQ.0
          comvl%subrot=comvl%tsubrot.AND.MOD(iteropt%nfi-1,comvl%nrotv).EQ.0
          ropt_mod%prteig=cprint%tprint.AND.MOD(iteropt%nfi-1,cprint%iprint_step).EQ.0
          cntl%caldip=cntl%tdipd.AND.MOD(iteropt%nfi-1,cnti%npdip).EQ.0
          IF (hubbu%pfrqom.gt.0) THEN
             hubbu%tpom=MOD(iteropt%nfi-1,hubbu%pfrqom).EQ.0
          ELSE
             hubbu%tpom=.False.
          ENDIF
          IF (.NOT.paral%parent) ropt_mod%prteig=.FALSE.
          ropt_mod%engpri=cprint%tprint.AND.MOD(iteropt%nfi-1,cprint%iprint_step).EQ.0
          ! ANNEALING
          bsclcs=1
          CALL anneal(velp,cm(1,1,1),nstate,scr)
          CALL berendsen(velp,cm(1,1,1),nstate,scr,ekinc,0.0_real_8)
          IF (cntl%bsymm)THEN
             bsclcs=2
             CALL anneal(velp,cm(1,1,2),nstate,scr)
             CALL berendsen(velp,cm(1,1,2),nstate,scr,ekinc,0.0_real_8)
          ENDIF
          ! SUBTRACT CENTER OF MASS VELOCITY
          IF (paral%parent.AND.comvl%subcom) CALL comvel(velp,vcmio,.TRUE.)
          ! SUBTRACT ROTATION AROUND CENTER OF MASS
          IF (paral%parent.AND.comvl%subrot) CALL rotvel(tau0,velp,lmio,tauio,.TRUE.)
          ! UPDATE NOSE THERMOSTATS
          CALL noseup(velp,cm,nstate,1)
          IF (cntl%bsymm)THEN
             bsclcs=2
             CALL noseup(velp,cm(1,1,2),nstate,2)
          ENDIF
          ! UPDATE VELOCITIES
          IF (paral%parent) CALL velupi(velp,fion,1)
          IF (.NOT.clc%classical) THEN
             CALL velupa(c0,cm,c2,nstate,1)
             IF (cntl%bsymm)CALL velupa(c0(1,1,2),cm(1,1,2),c2(1,1,2),nstate,1)
          ENDIF
          ! UPDATE POSITIONS
          ener_com%ecnstr = 0.0_real_8
          ener_com%erestr = 0.0_real_8
          IF (paral%parent) THEN
             CALL posupi(tau0,taup,velp)
             CALL mm_solv_const('SHAKE',dt_ions,taup,velp,tau0)
             IF (cotc0%mcnstr.NE.0) THEN
                CALL cpmdshake(tau0,taup,velp)
             ENDIF
          ENDIF
       ENDIF! qmnode

       CALL mp_bcast(taup,SIZE(taup),parai%qmmmsource,parai%qmmmgrp)
       IF (.NOT.lqmmm%qmmm_reflex) THEN
          CALL mm_translate_qmmm(taup,c0,cm,nstate)
          IF (cntl%bsymm) THEN
             CALL mm_translate_qmmm (taup,c0(1,1,2),cm(1,1,2),nstate)
          ENDIF
       ENDIF
       IF (paral%qmnode) THEN
          CALL mm_dim(mm_go_qm,statusdummy)
          CALL phfac(taup)
          IF (corel%tinlc) CALL copot(rhoe,psi,ropt_mod%calste)
          IF (.NOT.clc%classical)THEN
             CALL posupa(c0,cm,c2,gamx,nstate)
          ENDIF
          ! ..Dipole moment
          IF ((cntl%caldip.OR.(vdwl%vdwd.AND.vdwwfl%twannup)).AND.(.NOT.cntl%bsymm)) THEN
             CALL ddipo(taup,c0(:,:,1),cm(:,:,1),c2(:,:,1),sc0,nstate,center)
             CALL wannier_print(iteropt%nfi,c0(:,:,1),taup,nstate,psi(:,1),center)
          ENDIF
          ropt_mod%calste=cntl%tpres.AND.MOD(iteropt%nfi-1,cnti%npres).EQ.0
          CALL mm_dim(mm_go_mm,statusdummy)
       ENDIF! qmnode

       IF (cntl%bsymm)THEN
          CALL mm_qmmm_forcedr_bs(c0,c2,sc0,rhoe,psi,taup,fion,&
               eigv,nstate,0,.FALSE.,.TRUE.,&
               spd_bs,spd_hs,spda_bs,spda_hs,etoths,etotbs,.TRUE.)
       ELSE
          CALL mm_qmmm_forcedr(c0,c2,sc0,rhoe,psi,taup,fion,&
               eigv,nstate,0,.FALSE.,.TRUE.,.TRUE.)
       ENDIF
       ! ==================================================================
       ! Damped Dynamics
       CALL dampdyn(velp,fion,cm(1,1,1),c2(1,1,1),nstate,&
            scr(1),scr(10))
       IF (cntl%bsymm)THEN
          bsclcs=2
          CALL dampdyn(velp,fion,cm(1,1,2),c2(1,1,2),nstate,&
               scr(1),scr(10))
       ENDIF
       ! ==================================================================
       ! Meta Dynamics of Collective Variables

       IF (paral%qmnode)THEN
          IF (lmeta%lcolvardyn) THEN
             lquench = .FALSE.
             lmetares= .FALSE.

             IF (lmeta%lextlagrange) THEN
                ! Metadynamics with Extended Lagrangian
                IF (lmeta%tmulti) THEN
                   CALL meta_ext_mul(taup,velp,taur,&
                        lquench,lmetares,resetcv,ekinc,ekinp)
                ELSE

                   CALL meta_extlagr(taup,velp,taur,&
                        lquench,lmetares,resetcv,ekinc,ekinp)
                ENDIF
             ELSE
                ! Time dependent potential applied directly on the Collective Variables
                CALL meta_colvar(taup,velp,fion,taur,&
                     lquench,lmetares,ekinc,ekinp)
             ENDIF
             IF ((rmeta%tolkin.GT.0.0_real_8 .AND. ekinc.GT.rmeta%tolkin)&
                  .OR. lquench) THEN
                ! Check for Quench on the BO surface
                IF (paral%parent) THEN
                   IF (paral%io_parent)&
                        WRITE(6,'(''TOLKIN ='',f16.8,''  EKINC ='',f16.8)')&
                        rmeta%tolkin,ekinc
                ENDIF
                CALL mm_dim(mm_go_qm,statusdummy)
                bsclcs=1
                IF (cntl%bsymm)CALL setbsstate
                cntl%quenchb=.TRUE.
                CALL quenbo(c0(:,:,1),c2,sc0,taur,rhoe,psi)
                IF (cntl%bsymm)THEN
                   bsclcs=2
                   CALL setbsstate
                   CALL quenbo(c0(:,:,2),c2(1,1,2),sc0(1,1,2),taur,rhoe,&
                        psi)
                ENDIF
                CALL mm_dim(mm_go_mm,statusdummy)
                CALL zeroing(cm(:,:,1:nkpt%ngwk*bsfac))!,nstate*nkpt%ngwk*bsfac)
                resetcv = .TRUE.
                cntl%quenchb=.FALSE.
             ENDIF
             IF (paral%parent)THEN
                !$omp parallel do private(IS,IA)
                DO is = 1,ions1%nsp
                   DO ia = 1,ions0%na(is)
                      fion(1,ia,is) = fion(1,ia,is) + fhills(1,ia,is)
                      fion(2,ia,is) = fion(2,ia,is) + fhills(2,ia,is)
                      fion(3,ia,is) = fion(3,ia,is) + fhills(3,ia,is)
                   ENDDO
                ENDDO
             ENDIF
          ENDIF

          ! ==================================================================
          ! From saddle point to minima

          IF (lmeta%lsadpnt) THEN
             ! Check on the predefined  known minima 
             CALL tst2min(taup,taur,iteropt%nfi)
             ! If one minimum is found the search is initialized as new (random vel.)
             IF (lmeta%lmdreinit) THEN
                restart1%restart = .TRUE.
                IF ((paral%parent).AND.paral%io_parent)&
                     CALL fileclose(3)
                GOTO   99999
             ENDIF
          ENDIF
       ENDIF

       ! ==================================================================

       IF (paral%qmnode) THEN
          IF (ropt_mod%calste) CALL totstr
          ! FINAL UPDATE FOR VELOCITIES
          ener_com%ecnstr = 0.0_real_8
          ener_com%erestr = 0.0_real_8
          IF (paral%parent) THEN
             CALL velupi(velp,fion,1)
             IF (lqmmm%qmmm_reflex) CALL mm_qm_boundary(taup,velp)! cmb
             CALL mm_solv_const('RATTL',dt_ions,taup,velp,tau0)
             CALL rattle(taup,velp)
          ENDIF
          bsclcs=1
          IF (cntl%bsymm)CALL setbsstate

          IF (.NOT.clc%classical) THEN
             CALL velupa(c0,cm,c2,nstate,1)
             IF (cntl%bsymm)THEN
                bsclcs=2
                CALL setbsstate
                CALL velupa(c0(1,1,2),cm(1,1,2),c2(1,1,2),nstate,1)
             ENDIF

             CALL mm_dim(mm_go_qm,statusdummy)
             IF (cntl%bsymm)THEN
                bsclcs=1
                CALL setbsstate
             ENDIF
             CALL rortv(c0,cm,c2,sc0,gamy,nstate)
             IF (cntl%bsymm)THEN
                bsclcs=2
                CALL setbsstate
                CALL rortv(c0(1,1,2),cm(1,1,2),c2(1,1,2),sc0(1,1,2),&
                     gamy(1,2),nstate)
             ENDIF
          ENDIF
          CALL mm_dim(mm_go_mm,statusdummy)
          ! COMPUTE THE IONIC TEMPERATURE TEMPP
          IF (paral%parent) THEN

             CALL ekinpp(ekinp,velp)
             IF (lmeta%lextlagrange.AND. ltcglobal) THEN
                CALL ekincv_global(ek_cv)
                tempp=(ek_cv+ekinp)*factem*2._real_8/(glib+REAL(ncolvar,kind=real_8))
             ELSE
                tempp=ekinp*factem*2._real_8/glib
             ENDIF
             ! McB  calculate local TEMP  a f t e r  corrections on VELP()!
          ENDIF
          ! IONIC TEMPERATURE CONTROL
          IF (paral%parent) CALL rscvp(temp1,temp2,tempp,velp)
          ! SUBTRACT ROTATION AROUND CENTER OF MASS
          IF (paral%parent.AND.comvl%subrot)CALL rotvel(tau0,velp,lmio,tauio,.FALSE.&
               )
          ! SUBTRACT CENTER OF MASS VELOCITY
          IF (paral%parent.AND.comvl%subcom)CALL comvel(velp,vcmio,.FALSE.)
          ! UPDATE NOSE THERMOSTATS
          IF (cntl%tnosee.OR.cntl%tc) THEN
             CALL rekine(cm,nstate,ekinc)
             IF (cntl%bsymm)THEN
                ekinc_bs=ekinc
                CALL rekine(cm(1,1,2),nstate,ekinc)
                ekinc_hs=ekinc
             ENDIF
          ENDIF
          bsclcs=1
          CALL noseup(velp,cm,nstate,1)
          IF (cntl%bsymm)THEN
             bsclcs=2
             CALL noseup(velp,cm(1,1,2),nstate,2)
          ENDIF
          ! ANNEALING
          bsclcs=1
          CALL anneal(velp,cm(1,1,1),nstate,scr)
          CALL berendsen(velp,cm(1,1,1),nstate,scr,ekinc,0.0_real_8)
          IF (cntl%bsymm)THEN
             bsclcs=2
             CALL anneal(velp,cm(1,1,2),nstate,scr)
             CALL berendsen(velp,cm(1,1,2),nstate,scr,ekinc,0.0_real_8)
          ENDIF
          IF (paral%parent) THEN
             CALL ekinpp(ekinp,velp)
             IF (lmeta%lextlagrange.AND. ltcglobal) THEN
                CALL ekincv_global(ek_cv)
                tempp=(ek_cv+ekinp)*factem*2._real_8/(glib+REAL(ncolvar,kind=real_8))
             ELSE
                tempp=ekinp*factem*2._real_8/glib
             ENDIF
             ! McB  now the temperatures fit to velocities written to TRAJECTORY!
             ! calculate local kinetic temperature in QM and MM subsystem
             ! and write to "QM_TEMP"
             IF (wp_i%nfi_lt.GT.0) THEN
                CALL mm_localt(velp,rmass%pma,factem,tempp,glib,iteropt%nfi,wp_i%nfi_lt)
             ENDIF
          ENDIF
          ! RESCALE ELECTRONIC VELOCITIES
          IF (cntl%tc.AND.(.NOT.cntl%bsymm))&
               CALL rscve(ekin1,ekin2,ekinc,cntr%ekinw,cm,nstate,ncpw%ngw)
          IF (cntl%bsymm.AND.cntl%tc)THEN
             CALL rscve(ekin1,ekin2,ekinc_bs,cntr%ekinw,cm,nstate,ncpw%ngw)
             CALL rscve(ekin1,ekin2,ekinc_hs,cntr%ekinw,cm(1,1,2),nstate,ncpw%ngw)
          ENDIF
          ! MEAN SQUARE DISPLACEMENT OF DIFFERENT IONIC SPECIES
          IF (paral%parent) CALL dispp(taup,taui,disa)
          ! KINETIC ENERGY OF THE ELECTRONS
          CALL rekine(cm,nstate,ekinc)
          IF (cntl%bsymm)THEN
             ekinc_bs=ekinc
             CALL rekine(cm(1,1,2),nstate,ekinc)
             ekinc_hs=ekinc
             ekinc=(scalhs*ekinc_hs)+(scalbs*ekinc_bs)
          ENDIF
          ! ENERGY OF THE NOSE THERMOSTATS
          IF (paral%parent) CALL noseng(iteropt%nfi,velp,enose,enosp,dummy,1)
          IF (paral%parent) THEN
             IF (paral%io_parent)&
                  CALL fileclose(3)
             IF (paral%io_parent)&
                  CALL fileopen(3,filen,fo_app,ferror)
             econs=ekinp+ener_com%etot+enose+enosp+ener_com%ecnstr+ener_com%erestr+ekincv+vharm
             eham=econs+ekinc
             time2=m_walltime()
             tcpu=(time2-time1)*0.001_real_8
             ! Printing data of Broken symmetry calculation
             IF (cntl%bsymm) THEN
                couplj = (etoths - etotbs) * rtsasb
                jwn=couplj*autocm
                IF (paral%io_parent)&
                     WRITE(277,'(I10,11F16.8)')&
                     iteropt%nfi,ekinc_bs,ekinc_hs,ener_com%etot,etotbs,etoths,spd_bs,&
                     spd_hs,spda_bs,spda_hs,couplj,jwn
                CALL m_flush(277)
             ENDIF
             ! PRINTOUT the evolution of the accumulators every time step
             CALL wrprint_md(eigv,crge%f,ener_com%amu,nstate,taup,fion,&
                  ekinc,tempp,ener_com%etot,econs,eham,disa,&
                  tcpu,.FALSE.,iteropt%nfi,infi)
             ! UPDATE ACCUMULATORS
             CALL pacca(ekinc,tempp,ener_com%etot,econs,eham,enose,enosp,ener_com%ecnstr,&
                  ener_com%erestr,disa,tcpu,iteropt%nfi,1)
             ! Store ionic coordinates and velocities for statistics
             ropt_mod%movie= rout1%mout  .AND. MOD(iteropt%nfi-1,cnti%imovie).EQ.0
             ropt_mod%rprint=rout1%rout  .AND. MOD(iteropt%nfi-1,cnti%ntraj ).EQ.0
             ropt_mod%txyz=  rout1%xtout .AND. MOD(iteropt%nfi-1,cnti%ntraj ).EQ.0
             ropt_mod%tdcd=  rout1%dcout .AND. MOD(iteropt%nfi-1,cnti%ntraj ).EQ.0
             IF (MOD(iteropt%nfi-1,cnti%ntraj).EQ.0) THEN
                IF (paral%io_parent)&
                     CALL fileopen(83,'MM_CELL_TRANS',fo_app,ferror)
                IF (paral%io_parent)&
                     WRITE(83,'(I10,3F15.10)') iteropt%nfi,(clsaabox%mm_c_trans(k),k=1,3)
                IF (paral%io_parent)&
                     CALL fileclose(83)
             ENDIF

             CALL reshape_inplace(scr, (/3,maxsys%nax,maxsys%nsx/), scr_correct_shape )
             CALL printp(scr_correct_shape,taup,velp)
             IF (cprint%twriteforcetrajectory) CALL printp2(taur,taup,velp,fion)
          ENDIF

          IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)
          IF (infi.EQ.cnti%nomore) THEN
             soft_com%exsoft=.TRUE.
             soft_com%exnomore=.TRUE.
          ENDIF

       ENDIF! qmnode

       ! STOP THE RUN IF THE USER HAS SET THE SIGNAL 30
       in=0
       out=0
       IF (soft_com%exsoft) in=1
       CALL mp_sum(in,out,parai%qmmmgrp)
       IF (out.NE.0) soft_com%exsoft=.TRUE.

       IF (paral%qmnode) THEN
          ! periodic output of density/wavefunction etc.
          IF (rout1%rhoout.AND.(rout1%nrhoout.GT.0)) THEN
             IF (MOD(iteropt%nfi-1,rout1%nrhoout).EQ.0) THEN
                CALL rhopri(c0,tau0,rhoe,psi(:,1),nstate,nkpt%nkpnt)
             ENDIF
          ENDIF
          IF (teststore(iteropt%nfi).OR.soft_com%exsoft.OR.lmetares) THEN
             CALL zhwwf(2,irec,c0,cm,nstate,eigv,taup,velp,taui,iteropt%nfi)
          ENDIF

          IF (soft_com%exsoft .AND.lmeta%lcolvardyn) THEN
             lmetares= .TRUE.

             IF (lmeta%lextlagrange) THEN
                ! Metadynamics with Extended Lagrangian
                IF (lmeta%tmulti) THEN
                   CALL meta_ext_mul(taup,velp,taur,&
                        lquench,lmetares,resetcv,ekinc,ekinp)
                ELSE
                   CALL meta_extlagr(taup,velp,taur,&
                        lquench,lmetares,resetcv,ekinc,ekinp)
                ENDIF
             ELSE
                ! Time dependent potential applied directly on the Collective Variables
                CALL meta_colvar(taup,velp,fion,taur,&
                     lquench,lmetares,ekinc,ekinp)
             ENDIF
          ENDIF
          ! temperature ramping
          CALL tempramp(temp1,temp2)
       ENDIF! qmnode

       IF (soft_com%exsoft) GOTO 100

       IF (paral%qmnode) THEN
          ! UPDATE IONIC POSITIONS
          CALL dcopy(3*maxsys%nax*maxsys%nsx,taup(1,1,1),1,tau0(1,1,1),1)
       ENDIF! qmnode
       ! ==================================================================
       ! ==     END OF MAIN LOOP                                         ==
       ! ==================================================================
    ENDDO
100 CONTINUE
    CALL mp_sync(parai%qmmmgrp)
    CALL mm_dim(mm_go_mm,statusdummy)
    IF ( gparal%mmparent ) THEN
       CALL mm_write_gromos_coord('CRD_FIN.g96',taup,velp,maxsys%nax,maxsys%nsx)
    ENDIF
    IF (paral%qmnode) THEN
       IF (wannl%twann) THEN
          CALL mm_dim(mm_go_qm,statusdummy)
          CALL ddipo(taup,c0(:,:,1),cm(:,:,1),c2(:,:,1),sc0,nstate,center)
          CALL forcedr(c0(:,:,1),c2(:,:,1),sc0(:,:,1),rhoe,psi,taup,fion,eigv,&
               nstate,1,.FALSE.,.TRUE.)
          CALL wc_dos(c0,c2,nstate,center)
          CALL mm_dim(mm_go_mm,statusdummy)
       ENDIF
       CALL mm_dim(mm_go_qm,statusdummy)
       IF (rout1%rhoout.AND.(rout1%nrhoout.LE.0)) THEN
          CALL rhopri(c0,tau0,rhoe,psi(:,1),nstate,nkpt%nkpnt)
       ENDIF
       CALL mm_dim(mm_go_mm,statusdummy)
       ! PRINT ACCUMULATOR
       IF (paral%parent) CALL pacca(ekinc,tempp,ener_com%etot,econs,eham,enose,enosp,&
            ener_com%ecnstr,ener_com%erestr,disa,tcpu,iteropt%nfi,2)
       bsclcs=1
       IF (cntl%bsymm)CALL setbsstate
       CALL mm_dim(mm_go_qm,statusdummy)
       CALL proja(c0,c2,sc0,nstate,cnti%iproj)
       IF (cntl%bsymm)THEN
          bsclcs=2
          CALL setbsstate
          CALL proja(c0(1,1,2),c2(1,1,2),sc0(1,1,2),&
               nstate,cnti%iproj)
       ENDIF
       bsclcs=1
       IF (cntl%bsymm)CALL setbsstate
       CALL csize(c2,crge%n,gemax,cnorm)
       IF (cntl%bsymm)THEN
          bsclcs=2
          CALL setbsstate
          CALL csize(c2(1,1,2),crge%n,gemax,cnorm)
       ENDIF
       IF (paral%parent) CALL gsize(fion,gnmax,gnorm)
       CALL mm_dim(mm_go_mm,statusdummy)
       IF (paral%parent) CALL finalp(tau0,fion,velp,eigv)
       ! 
    ENDIF                     ! qmnode
    DEALLOCATE(center,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    IF (cntl%tsampl) THEN
       CALL sample_go
       GOTO 99999
    ENDIF
10000 CONTINUE
    ! 
    IF (paral%qmnode) THEN
       IF (cntl%nonort) THEN
          DEALLOCATE(eigv,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
          DEALLOCATE(eigm,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
       ENDIF
       IF ((paral%parent).AND.paral%io_parent)&
            CALL fileclose(3)
       IF ((cntl%bsymm.AND.paral%parent).AND.paral%io_parent)&
            CALL fileclose(277)
       DEALLOCATE(rhoe,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(psi,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(scr,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL mm_dim(mm_revert,oldstatus)
    CALL mp_sync(parai%qmmmgrp)
    CALL tihalt(procedureN,isub)
#endif
    RETURN
  END SUBROUTINE mm_mdmain
  ! ==================================================================
  SUBROUTINE give_scr_mm_mdmain(lmdmain,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lmdmain
    CHARACTER(len=30)                        :: tag

    INTEGER :: lcopot, lddipo, ldeort, lforcedr, linitrun, lmtd, lortho, &
      lposupa, lquenbo, lrhopri, lrortv, nstate

    nstate=crge%n
    linitrun=0
    lcopot=0
    lortho=0
    lquenbo=0
    ldeort=0
    lrhopri=0
    lddipo=0
    CALL give_scr_initrun(linitrun,tag)
    IF (corel%tinlc) CALL give_scr_copot(lcopot,tag)
    IF (cntl%trane) CALL give_scr_ortho(lortho,tag,nstate)
    IF (cntl%quenchb) CALL give_scr_quenbo(lquenbo,tag)
    IF (pslo_com%tivan) CALL give_scr_deort(ldeort,tag,nstate)
    IF (cntl%tdipd) CALL give_scr_ddipo(lddipo,tag)
    CALL give_scr_forcedr(lforcedr,tag,nstate,.FALSE.,.TRUE.)
    CALL give_scr_rortv(lrortv,tag,nstate)
    CALL give_scr_posupa(lposupa,tag,nstate)
    IF (rout1%rhoout) CALL give_scr_rhopri(lrhopri,tag,nstate)
    CALL give_scr_meta_extlagr(lmtd,tag)
    lmdmain=MAX(lcopot,lortho,lquenbo,ldeort,lforcedr,&
         lrortv,lposupa,lrhopri,lddipo,linitrun,lmtd)
    IF (cntl%tqmmm)lmdmain=MAX(lmdmain,fpar%kr1*fpar%kr2s*fpar%kr3s)
    IF (cntl%tqmmm)lmdmain=MAX(lmdmain,maxsys%nax*maxsys%nsx*3)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_mm_mdmain
  ! ==================================================================
  SUBROUTINE mm_localt(velp,pma,factem,tempp,glib,nfi,NFI_lt)
    ! ==--------------------------------------------------------------==
    ! Compute kinetic temperature for QM and MM subsystem separately
    ! does not rely on CAFES being used and is more accurate!
    ! initially written by udo schmitt march 2003
    ! revised by mauro boero, axel kohlmeyer march/april 2005
    ! updated for separate solvent temperature and constraints, august 2005, AK.
    ! 
    REAL(real_8)                             :: velp(:,:,:), pma(*), factem, &
                                                tempp, glib
    INTEGER                                  :: nfi, NFI_lt

    INTEGER                                  :: atmidx, ia, iat, is, ndof
    LOGICAL                                  :: ferror
    REAL(real_8)                             :: const, ekincl, ekinmm, ekinp, &
                                                ekinqm, mmcnstr, mmdof, &
                                                qmcnstr, qmdof, tempcl, &
                                                tempmm, tempqm

! ==--------------------------------------------------------------==

#if defined (__GROMOS)
    IF (MOD(nfi-1,NFI_lt).NE.0) RETURN

    ! compute constraints contributions to DOFs.
    qmcnstr=0.0_real_8
    mmcnstr=0.0_real_8
    IF (cotc0%mcnstr.GT.0) THEN
       DO is=1,cotc0%mcnstr
          ! set prefactor depending on constraint type.
          IF (ntcnst(1,is).EQ.1) THEN
             const=0.5_real_8
          ELSEIF (ntcnst(1,is).EQ.2) THEN
             const=1.0_real_8/3.0_real_8
          ELSEIF (ntcnst(1,is).EQ.3) THEN
             const=0.25_real_8
          ELSEIF (ntcnst(1,is).EQ.4) THEN
             const=0.5_real_8
          ELSEIF (ntcnst(1,is).EQ.5) THEN
             const=0.25_real_8
          ELSEIF (ntcnst(1,is).EQ.6) THEN
             const=1.0_real_8
          ELSEIF (ntcnst(1,is).EQ.7) THEN
             const=1.0_real_8/3.0_real_8
          ELSEIF (ntcnst(1,is).EQ.8) THEN
             const=1.0_real_8
          ENDIF
          ! loop over list of atoms in constraint. unused entries are supposed to be 0.
          DO ia=2,6
             atmidx=ntcnst(ia,is)
             IF (atmidx.LE.0.OR.atmidx.GT.mmdim%natm) THEN
                const=0.0_real_8     ! do nothing
             ELSEIF (atmidx.LE.mmdim%natq) THEN
                qmcnstr=qmcnstr + const! atom is quantum
             ELSE IF (atmidx.LE.solsolv%nrpt) THEN
                mmcnstr=mmcnstr + const! atom is solute but not quantum
             ENDIF
          ENDDO
       ENDDO
    ENDIF

    ! QM atoms in solute
    iat=0
    ndof=0
    ekinqm=0.0_real_8
    DO is=1,mmdim%nspq
       const=0.5_real_8*pma(is)
       DO ia=1,NAq(is)
          iat=iat+1
          ndof=ndof+lskcor(1,iat)+lskcor(2,iat)+lskcor(3,iat)
          ekinqm=ekinqm+const*(velp(1,ia,is)*velp(1,ia,is)&
               +VELP(2,IA,IS)*VELP(2,IA,IS)&
               +VELP(3,IA,IS)*VELP(3,IA,IS))
       ENDDO
    ENDDO
    qmdof=REAL(ndof,kind=real_8)-qmcnstr
    IF (qmdof.GT.0.1_real_8) THEN
       tempqm=ekinqm*factem*2._real_8/qmdof
    ELSE
       tempqm=0.0_real_8
    ENDIF

    ! MM atoms in solute not counting the QM atoms
    ndof=0
    ekinmm=0.0_real_8
    DO is=mmdim%nspq+1,mmdim%nspm
       const=0.5_real_8*pma(is)
       DO ia=1,NAm(is)
          IF (gratom((is-1)*maxsys%nax+ia).LE.solsolv%nrpt) THEN
             iat=iat+1
             ndof=ndof+lskcor(1,iat)+lskcor(2,iat)+lskcor(3,iat)
             ekinmm=ekinmm+const*(velp(1,ia,is)*velp(1,ia,is)&
                  +VELP(2,IA,IS)*VELP(2,IA,IS)&
                  +VELP(3,IA,IS)*VELP(3,IA,IS))
          ENDIF
       ENDDO
    ENDDO
    mmdof=REAL(ndof,kind=real_8)-mmcnstr
    ! account for uncounted constraints in case we have no solvent
    IF (glib.LT.(qmdof+mmdof)) mmdof=glib-qmdof
    IF (mmdof.GT.0.1_real_8) THEN
       tempmm=ekinmm*factem*2._real_8/mmdof
    ELSE
       tempmm=0.0_real_8
    ENDIF

    ekinp=tempp/factem/2.0_real_8*glib
    ekincl=ekinp-ekinqm-ekinmm
    IF (glib.GT.(qmdof+mmdof)) THEN
       tempcl=ekincl*factem*2._real_8/(glib-qmdof-mmdof)
    ELSE
       tempcl=0.0_real_8
    ENDIF
    IF (tempcl.GT.0.0_real_8) THEN
       IF (paral%io_parent)&
            WRITE(6,'(/,A,F8.2,A,I4,1X,A,F8.2,A,I7,1X,A,F8.2,/)')&
            ' T(QM)=',tempqm,' DOF(QM)=',NINT(qmdof),&
            ' T(MM)=',tempmm,' DOF(MM)=',NINT(mmdof),&
            ' T(SOL)=',tempcl
    ELSE
       IF (paral%io_parent)&
            WRITE(6,'(/,A,F8.2,A,I4,1X,A,F8.2,A,I7,/)')&
            ' T(QM)=',tempqm,' DOF(QM)=',NINT(qmdof),&
            ' T(MM)=',tempmm,' DOF(MM)=',NINT(mmdof)
    ENDIF
    IF (paral%io_parent)&
         CALL fileopen(49,'QM_TEMP',fo_app,ferror)
    IF (paral%io_parent)&
         WRITE(49,'(I6,4F12.4)') nfi,tempqm,tempmm,tempcl,tempp
    IF (paral%io_parent)&
         CALL fileclose(49)
    ! ... store EKINQM for surface hopping
    ekinqmsh=ekinqm
#endif
    RETURN
  END SUBROUTINE mm_localt
  ! ==================================================================

END MODULE mm_mdmain_utils
