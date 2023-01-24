MODULE mdmain_utils
  USE anneal_utils,                    ONLY: anneal,&
                                             berendsen,&
                                             dampdyn,&
                                             tempramp
  use bicanonicalCpmd, only: bicanonicalTestExit,&
  biCanonicalEnsembleDo, getNameEnergiesTape, &
  bicanonicalCpmdConfig, PrintEnergies,&
  CpmdEnergiesGradients, SanityChecks
  USE box_boundary_utils,              ONLY: box_boundary
  USE bsym,                            ONLY: autocm,&
                                             bsclcs,&
                                             bsfac,&
                                             cnstwgt,&
                                             ekinc_bs,&
                                             ekinc_hs,&
                                             rtsasb
  USE bsympnt,                         ONLY: fnbs
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
  USE cotr,                            ONLY: cnsval,&
                                             cotc0,&
                                             ntcnst
  USE cppt,                            ONLY: inyh
  USE csize_utils,                     ONLY: csize
  USE ddipo_utils,                     ONLY: ddipo,&
                                             give_scr_ddipo
  USE deort_utils,                     ONLY: deort,&
                                             give_scr_deort
  USE detdof_utils,                    ONLY: detdof
  USE dispp_utils,                     ONLY: dispp
  USE dynit_utils,                     ONLY: dynit
  USE efld,                            ONLY: extf
  USE ekinpp_utils,                    ONLY: ekinpp
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: chrg,&
                                             ener_com
  USE epr_efg_utils,                   ONLY: epr_efg
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
  USE geofile_utils,                   ONLY: geofile
  USE geq0mod,                         ONLY: geq0
  USE gle_utils,                       ONLY: gle_init,&
                                             gle_step
  USE glemod,                          ONLY: glepar
  USE gsize_utils,                     ONLY: gsize
  USE hfxmod,                          ONLY: hfxc3
  USE hubbardu,                        ONLY: hubbu
  USE initrun_driver,                  ONLY: initrun
  USE initrun_utils,                   ONLY: give_scr_initrun
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE isos,                            ONLY: isos1
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE localize_utils,                  ONLY: localize2
  USE lsforce_utils,                   ONLY: lsforce
  USE machine,                         ONLY: m_walltime
  USE meta_colvar_inp_utils,           ONLY: colvar_structure
  USE meta_colvar_utils,               ONLY: meta_colvar,&
                                             meta_colvar_mw
  USE meta_exl_mult_utils,             ONLY: meta_ext_mul
  USE meta_exlagr_methods,             ONLY: give_scr_meta_extlagr,&
                                             meta_extlagr,&
                                             meta_extlagr_mw
  USE meta_exlagr_utils,               ONLY: ekincv_global
  USE meta_multiple_walkers_utils,     ONLY: mw_assign_filenames,&
                                             mw_filename
  USE mimic_wrapper,                   ONLY: mimic_ifc_collect_energies,&
                                             mimic_ifc_collect_forces,&
                                             mimic_ifc_sort_fragments,&
                                             mimic_ifc_init,&
                                             mimic_revert_dim,&
                                             mimic_save_dim,&
                                             mimic_ifc_send_coordinates,&
                                             mimic_sum_forces,&
                                             mimic_switch_dim,&
                                             mimic_update_coords,&
                                             mimic_control,&
                                             mimic_energy,&
                                             mimic_subsystem_temperatures,&
                                             mimic_subsystem_dof
  USE mfep,                            ONLY: mfepi
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sync
  USE mw,                              ONLY: mwi,&
                                             tmw
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
  USE pimd,                            ONLY: grandparent,&
                                             ipcurr,&
                                             np_low,&
                                             supergroup
  USE posupa_utils,                    ONLY: give_scr_posupa,&
                                             posupa
  USE posupi_utils,                    ONLY: posupi
  USE printave_utils,                  ONLY: pacca
  USE printp_utils,                    ONLY: printp,&
                                             printp2
  USE proja_utils,                     ONLY: proja
  USE pslo,                            ONLY: pslo_com
  USE puttau_utils,                    ONLY: taucl
  USE quenbo_utils,                    ONLY: give_scr_quenbo,&
                                             quenbo
  USE rattle_utils,                    ONLY: rattle
  USE readsr_utils,                    ONLY: xstring
  USE rekine_utils,                    ONLY: rekine
  USE resetac_utils,                   ONLY: resetac
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
  USE shake_utils,                     ONLY: cpmdshake,&
                                             init_constraints,&
                                             do_shake,&
                                             do_rattle
  USE soft,                            ONLY: soft_com
  USE spin,                            ONLY: spin_mod
  USE store_types,                     ONLY: &
       cprint, irec_ac, irec_noe, irec_nop1, irec_nop2, irec_nop3, irec_nop4, &
       irec_vel, restart1, rout1
  USE system,                          ONLY: &
       cnti, cntl, cntr, fpar, maxsys, nacc, ncpw, nkpt, restf
  USE testex_utils,                    ONLY: testex,&
                                             testex_mw
  USE teststore_utils,                 ONLY: teststore
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE totstr_utils,                    ONLY: totstr
  USE tst2min_utils,                   ONLY: tst2min
  USE utils,                           ONLY: zclean,&
                                             zclean_k
  USE vdwcmod,                         ONLY: trwanncx,&
                                             twannupx,&
                                             vdwl,&
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

  PUBLIC :: mdmain
  PUBLIC :: give_scr_mdmain

CONTAINS

  ! ==================================================================
  SUBROUTINE mdmain(c0,cm,c2,sc0,gamx,gamy)
    ! ==--------------------------------------------------------------==

    ! 
    COMPLEX(real_8)                          :: c0(:,:,:), &
                                                cm(ncpw%ngw,crge%n,bsfac), &
                                                c2(:,:,:), &
                                                sc0(ncpw%ngw,crge%n,*)
    REAL(real_8)                             :: gamx(crge%n*crge%n,*), &
                                                gamy(crge%n*crge%n,*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'mdmain'

    CHARACTER(len=100)                       :: filen, filenbs, ftmp
    CHARACTER(len=30)                        :: tag
    COMPLEX(real_8), ALLOCATABLE             :: psi(:,:)
    INTEGER :: i, i1, i2, ia, ierr, il_psi_1d, il_psi_2d, il_rhoe_1d, &
      il_rhoe_2d, ipwalk, ipx, irec(100), is, isub, itemp, lscr, nstate
    LOGICAL                                  :: ferror, ionode, lmetares, &
                                                lquench, resetcv, tstrng
    REAL(real_8) :: couplj, disa, dum1, dum2, dummy(10), econs, eham, ek_cv, &
      ekin1, ekin2, ekinc, ekincp, ekinh1, ekinh2, ekinp, enosbs, enose, &
      enoshs, enosp, etotbs, etoths, ff, jwn, lmio(3), scalbs, scalhs, &
      spd_bs, spd_hs, spda_bs, spda_hs, tcpu, temp1, temp2, tempp, time1, &
      time2, vcmio(4)
    REAL(real_8), ALLOCATABLE :: center(:,:), eigm(:,:), eigv(:,:), &
      rhoe(:,:), scr(:), taui(:,:,:), tauio(:,:), taur(:,:,:)

#if defined (__QMECHCOUPL)
    ! QM/MM
    LOGICAL :: qmmech
    REAL(real_8) :: mm_ekin
    REAL(real_8) :: mm_temp
#endif
    CALL tiset(procedureN,isub)
    IF (cntl%mimic) THEN
      CALL mimic_save_dim()
      CALL mimic_switch_dim(go_qm=.FALSE.)
    ENDIF
    time1 =m_walltime()
    ! Walker ID for multiple walker MTD. MTD part is in grandparent
    ipwalk=1
    ionode=paral%io_parent !bugfix
    filen='ENERGIES'
    filenbs='BS_ENERG'
    IF (tmw)THEN
       ipwalk=mwi%walker_id
       ionode=grandparent
       CALL mw_assign_filenames(ipwalk,filen,filenbs)
    else if (biCanonicalEnsembleDo) then
      filen = getNameEnergiesTape(bicanonicalCpmdConfig)
    ELSE IF (cntl%tpmin)THEN
       ftmp='ENERGIES_'
       CALL mw_filename(ftmp,filen,ipcurr)
       CALL xstring(filen,i1,i2)
       ftmp=filen(i1:i2)//'.'
       CALL mw_filename(ftmp,filen,mfepi%istring)
       ftmp='BS_ENERG_'
       CALL mw_filename(ftmp,filenbs,ipcurr)
       CALL xstring(filenbs,i1,i2)
       ftmp=filenbs(i1:i2)//'.'
       CALL mw_filename(ftmp,filenbs,mfepi%istring)
    ENDIF
    !
    IF (.NOT. cntl%mimic) THEN
       ALLOCATE(taup(3,maxsys%nax,maxsys%nsx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(taup)!, 3*maxsys%nax*maxsys%nsx)
       ALLOCATE(fion(3,maxsys%nax,maxsys%nsx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ELSE
        ALLOCATE(extf(fpar%kr1*fpar%kr2s*fpar%kr3s),STAT=ierr)
    END IF
    CALL zeroing(fion)!,SIZE(fion))
    ALLOCATE(taui(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(taur(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    IF (comvl%tsubrot)  THEN
       ALLOCATE(tauio(3,ions1%nat),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    IF (cntl%bsymm)    THEN
       ALLOCATE(fnbs(3*maxsys%nax*maxsys%nsx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    IF (lmeta%lcolvardyn) THEN
       ALLOCATE(fhills(3,maxsys%nax,maxsys%nsx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(fhills)!,3*maxsys%nax*maxsys%nsx)
    ENDIF
    ! Initialize logical variable for Metadynamics
    lmetares=.FALSE.
    nstate=crge%n
    nacc = 22
    ropt_mod%modens=.FALSE.
    ropt_mod%engpri=.FALSE.
    ropt_mod%calste=cntl%tpres
    IF (pslo_com%tivan) cntl%nonort=.TRUE.
    IF (cntl%nonort) THEN
       ALLOCATE(eigv(crge%n,bsfac),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(eigm(crge%n*crge%n,bsfac),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ELSE
       ! Avoid Fortran runtime 'not allocated' error
       ALLOCATE(eigv(1,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(eigm(1,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
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
             IF (paral%io_parent) THEN
                WRITE(6,'(A,A)') ' HARMONIC REFERENCE SYSTEM INTEGRATION',&
                     ' ONLY POSSIBLE WITH EQUAL OCCUPATION NUMBERS'
             ENDIF
             CALL stopgm('MDMAIN','FATAL ERROR',& 
                  __LINE__,__FILE__)
          ENDIF
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    ! SCR ALLOCATION AND PARTITION (SCRATCH ARRAY).
    CALL rhoe_psi_size(il_psi_1d=il_psi_1d,il_psi_2d=il_psi_2d,&
         il_rhoe_1d=il_rhoe_1d,il_rhoe_2d=il_rhoe_2d)
    ALLOCATE(rhoe(fpar%nnr1,il_rhoe_2d),STAT=ierr) !vw doenst work in parallel (should be equal to nnr1) !il_rhoe_1d
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(psi(il_psi_1d,il_psi_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL give_scr_mdmain(lscr,tag)
    ALLOCATE(scr(lscr),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
99999 IF (cntl%tsampl) THEN
       soft_com%exsoft=.FALSE.
       CALL sample_wait
       IF (cnti%nomore.LT.0) GOTO 10000
    ENDIF
    restf%nfnow=1
    IF (cntl%mimic) THEN
       CALL mimic_ifc_init(rhoe, extf, tau0)
       mimic_control%update_potential = .TRUE.
    END IF
    ! ==--------------------------------------------------------------==
    ! TIME STEP FUNCTIONS
    CALL dynit(ekincp,ekin1,ekin2,temp1,temp2,ekinh1,ekinh2)
    ekinc=0.0_real_8    ! McB: cf. META_EXT..()
    ekinp=0.0_real_8    ! McB: cf. META_EXT..()
    ! PARAMETERS FOR THE NOSE-HOOVER THERMOSTATS
    IF ((cntl%tnosee.OR.cntl%tnosep).AND.paral%parent) CALL nosepa(ipwalk,ipwalk)
    ! Dont symmetrize density 
    cntl%tsymrho=.FALSE.
    IF (cntl%mimic) THEN
       CALL mimic_switch_dim(go_qm=.TRUE.)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == INITIALIZATION                                               ==
    ! ==--------------------------------------------------------------==
    ! BROKEN SYMMETRY: INITIALIZE BROKEN SYMM. STATE
333 FORMAT (/,t2,a,/)
    IF (cntl%bsymm)THEN
       bsclcs=1
       IF (paral%io_parent)&
            WRITE(6,333) 'BROKEN SYMMETRY STATE IS INITIALIZING'
    ENDIF
    CALL initrun(irec,c0,cm,sc0,rhoe,psi,eigv)

    if (biCanonicalEnsembleDo .and. paral%parent) call SanityChecks(bicanonicalCpmdConfig)

    ! BROKEN SYMMETRY: INITIALIZE HIGH SPIN STATE
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
    ! 
    CALL mp_bcast(taup,SIZE(taup),parai%source,parai%allgrp)
    CALL dcopy(3*maxsys%nax*maxsys%nsx,taup,1,taui,1)
    ! INITIALIZE WF CENTERS & SPREAD
    IF (vdwl%vdwd) THEN
       IF (cntl%tpath) THEN
          ipx=ipcurr-np_low+1
       ELSE
          ipx=1
       ENDIF
       vdwwfl%trwannc=trwanncx(ipx)
    ENDIF
    IF (hfxc3%twscr.OR.(vdwl%vdwd.AND..NOT.vdwwfl%trwannc)) THEN
       CALL localize2(tau0,c0,c2,sc0,nstate)
    ENDIF

    IF (cntl%mimic) THEN
       CALL mimic_switch_dim(go_qm=.TRUE.)
    ENDIF
    ! 
    ! NN: BROKEN SYMMETRY: QUENCHING TO BO SURFACE OF BS STATE
    IF (cntl%quenchb) THEN
       IF (cntl%bsymm)THEN
          bsclcs=1
          CALL setbsstate
          IF (paral%io_parent)&
               WRITE(6,333) 'QUENCHING BROKEN SYMMETRY STATE'
       ENDIF
       CALL quenbo(c0(:,:,1),c2,sc0,taur,rhoe,psi)

       CALL zhwwf(2,irec,c0,cm,nstate,eigv,tau0,velp,taui,iteropt%nfi)
       ! NN: BROKEN SYMMETRY: QUENCHING TO BO SURFACE OF HS STATE
       IF (cntl%bsymm)THEN
          bsclcs=2
          CALL setbsstate
          IF (paral%io_parent)&
               WRITE(6,333) 'QUENCHING HIGH SPIN STATE'
          CALL quenbo(c0(:,:,2),c2(:,:,2),sc0(1,1,2),taur,&
               rhoe,psi)
          CALL zhwwf(2,irec,c0(:,:,2),cm(1,1,2),nstate,eigv(1,2),&
               tau0,velp,taui,iteropt%nfi)
       ENDIF
       cntl%quenchb=.FALSE.
       IF (tmw) CALL mp_sync(supergroup)
    ENDIF
    ! 
    IF (pslo_com%tivan) THEN
       IF (cntl%tlsd) THEN
          bsclcs=1
          IF (cntl%bsymm)CALL setbsstate
          CALL deort(ncpw%ngw,spin_mod%nsup,eigm,eigv,c0(:,1:spin_mod%nsup,1),sc0(1,1,1))
          CALL deort(ncpw%ngw,spin_mod%nsdown,eigm,eigv,c0(:,spin_mod%nsup+1:spin_mod%nsup+spin_mod%nsdown,1),&
               sc0(1,spin_mod%nsup+1,1))
          ! 
          IF (cntl%bsymm)THEN
             bsclcs=2
             CALL setbsstate
             CALL deort(ncpw%ngw,spin_mod%nsup,eigm(1,2),eigv(1,2),c0(:,1:spin_mod%nsup,2),&
                  sc0(1,1,2))
             CALL deort(ncpw%ngw,spin_mod%nsdown,eigm(1,2),eigv(1,2),c0(:,spin_mod%nsup+1:spin_mod%nsup+spin_mod%nsdown,2),&
                  sc0(1,spin_mod%nsup+1,2))
          ENDIF
       ELSE
          CALL deort(ncpw%ngw,nstate,eigm,eigv,c0,sc0)
       ENDIF
    ENDIF

    IF (cntl%mimic) THEN
       CALL mimic_switch_dim(go_qm=.FALSE.)
    ENDIF

    ! INITIALIZE VELOCITIES
    IF (paral%parent) CALL detdof(tau0,taur)

    ! INITIALIZE CONSTRAINTS
    IF (paral%io_parent.AND.cntl%new_constraints) THEN
       WRITE(6, '(1x,a)') "USING NEW CONSTRAINTS SOLVER"
       IF (cntl%pbicgstab) then
          WRITE(6, '(1x,a)') "WITH PRECONDITIONED BICONJUGATE GRADIENT STABILIZED"
       ELSE
          WRITE(6, '(1x,a)') "WITH PRECONDITIONED CONJUGATE GRADIENT"
       ENDIF
       WRITE(6, '(1x,a,t60,i6)') "MAX NUMBER OF SHAKE ITERATIONS:", &
                                  cnti%shake_maxstep
       WRITE(6, '(1x,a,t60,i6)') "MAX NUMBER OF PCG STEPS IN EACH SHAKE ITERATION:", &
                                  cnti%shake_cg_iter
       CALL init_constraints(ntcnst, cnsval, ions0%na, ions1%nsp)
       IF (.NOT.mimic_control%external_constraints) THEN
          WRITE(6, '(1x,a)') "EXTERNAL CONSTRAINTS DEFINED THROUGH MIMIC WILL BE IGNORED"
       ENDIF
       IF (cntl%mimic) THEN
          CALL mimic_subsystem_dof()
       END IF
    END IF

    ! INITIALIZE METADYNAMICS VARIABLES used also for 
    ! SHOOTING from SADDLE POINT with RANDOM VELOCITIES
    IF (ionode .AND. (lmeta%lcolvardyn .OR. lmeta%lsadpnt)) THEN
       CALL colvar_structure(tau0,taur)
    ENDIF

    IF (irec(irec_vel).EQ.0.AND..NOT.restart1%rgeo) THEN
       ener_com%ecnstr = 0.0_real_8
       ener_com%erestr = 0.0_real_8
       CALL rinvel(velp,cm,crge%n)
       IF (paral%parent) CALL taucl(velp)
       IF (paral%parent) THEN
          IF (cntl%new_constraints) THEN
             CALL do_rattle(tau0, velp)
          ELSE
             CALL rattle(tau0,velp)
          END IF
       END IF
       CALL rvscal(velp)
    ELSE
       IF (paral%parent) CALL taucl(velp)
       IF (paral%parent) THEN
          IF (cntl%new_constraints) THEN
             CALL do_rattle(tau0, velp)
          ELSE
             CALL rattle(tau0,velp)
          END IF
       END IF
       IF (cntl%trescale) CALL rvscal(velp)
    ENDIF
    ! FOR NO RESTART THE VELOCITIES OF THE COEFFICIENTS OF THE 
    ! BS_STATE AND HS_STATE ARE SET TO ZERO
    IF (cntl%bsymm.AND.irec(irec_vel).EQ.0)&
         CALL zeroing(cm(:,:,1:bsfac))!,ngw*nstate*bsfac)
    IF (cntl%quenchp) CALL zeroing(velp)!,3*maxsys%nax*maxsys%nsx)
    IF (cntl%quenche) CALL zeroing(cm(:,:,1:bsfac))!,ngw*nstate*bsfac)
    IF (cntl%trevers) THEN
       ! invert electronic and ionic velocities (useful for path sampling)
       CALL dscal(2*ncpw%ngw*nstate*bsfac,-1._real_8,cm,1)
       CALL dscal(3*maxsys%nax*maxsys%nsx,-1._real_8,velp,1)
    ENDIF
    ! >>>
    ! make sure the velocities are correctly replicated while using groups      
    CALL mp_bcast(velp,SIZE(velp),0,parai%cp_grp)
    CALL mp_bcast(cm,ncpw%ngw*nstate*bsfac,0,parai%cp_inter_grp)
    ! <<<
    ! RESET ACCUMULATORS
    IF (paral%parent.AND.irec(irec_ac).EQ.0)&
         CALL resetac(tau0,taui,iteropt%nfi)
    ! ..-------------------------------------------------
    ! ..QM/MM coupling by Roethlisberger Group
    ! ..Initailize Molecular Mechanics subsystem:
    ! ..-------------------------------------------------
    ! ..
    IF (cntl%mimic) THEN
      CALL mimic_switch_dim(go_qm=.TRUE.)
    ENDIF
#if defined (__QMECHCOUPL)
    IF (paral%parent) THEN
       ! ---------------------------------------------------------
       ! CORE QM/MM initilization
       ! Please note this is the same routine used to initialize
       ! the QMMM functionality during a geometry optimization.
       ! cntl%md specific initialization is performed by the routine
       ! 'mm_cpmd_md_init'.
       ! ---------------------------------------------------------
       CALL mm_cpmd_init (cntl%tqmmech,tau0,ions0%na,ions1%nsp,maxsys%nax,ions1%nat,restart1%rco,restart1%rvel,cntr%tempw)
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
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,'(1X,64("="))')
       IF (paral%io_parent)&
            WRITE(6,'(1X,"==",T25,A,T64,"==")') 'FORCES INITIALIZATION'
       IF (paral%io_parent)&
            WRITE(6,'(1X,64("="))')
    ENDIF
    IF (tkpts%tkpnt) THEN
       IF (geq0) CALL zclean_k(c0,nstate,ncpw%ngw)
    ELSE
       IF (geq0) CALL zclean(c0,nstate,ncpw%ngw)
    ENDIF
    ! FORCES FOR BROKEN SYMMETRY STATE 
    IF (cntl%bsymm)THEN
       bsclcs=1
       CALL setbsstate
    ENDIF
    IF (cntl%mimic) THEN
       CALL mimic_update_coords(tau0, c0, cm, nstate, ncpw%ngw, inyh)
       IF (paral%io_parent) THEN
          CALL mimic_ifc_send_coordinates()
          IF (mimic_control%tot_scf_energy) THEN
             CALL mimic_ifc_collect_energies()
          END IF
       END IF
    END IF
    CALL forcedr(c0(:,:,1),c2(:,:,1),sc0(:,:,1),rhoe,psi,&
         TAU0,FION,EIGV,NSTATE,1,.FALSE.,.TRUE.)
    ! STORE THE SPIN DENSITIES FOR PRINTING
    IF (cntl%bsymm) THEN
       spd_bs=chrg%csums
       spda_bs=chrg%csumsabs
       ! FORCES FOR BROKEN SYMMETRY STATE 
       IF (paral%parent)THEN
          etotbs=ener_com%etot
          ! STORING THE FORCES FOR BROKEN SYMMETRY STATE
          CALL dcopy(3*maxsys%nax*maxsys%nsx,fion,1,fnbs,1)
       ENDIF
       ! FORCES FOR HIGH SPIN STATE 
       bsclcs=2
       CALL setbsstate
       CALL forcedr(c0(:,:,2),c2(:,:,2),sc0(:,:,2),rhoe,psi,&
            tau0,fion,eigv(1,2),nstate,1,.FALSE.,.TRUE.)
       ! STORE THE SPIN DENSITIES FOR PRINTING
       spd_hs =chrg%csums
       spda_hs=chrg%csumsabs
       ! EVALUATE ENERGY AND FORCE FOR LOW SPIN STATE 
       IF (paral%parent)THEN
          etoths=ener_com%etot
          ener_com%etot = (1.0_real_8 + cnstwgt) * etotbs - cnstwgt * etoths
          CALL lsforce(fnbs,fion)
       ENDIF
    ENDIF
    ! Initialize Metadynamics contributions
    IF (lmeta%lcolvardyn .AND. lmeta%lextlagrange) THEN
       lquench = .FALSE.
       lmetares= .FALSE.
       resetcv = .FALSE.
       IF (lmeta%tmulti) THEN
          CALL meta_ext_mul(tau0,velp,taur,&
               lquench,lmetares,resetcv,ekinc,ekinp)
       ELSE IF (tmw)THEN
          CALL meta_extlagr_mw(tau0,velp,taur,&
               lquench,lmetares,resetcv,ekinc,ekinp)
       ELSE
          CALL meta_extlagr(tau0,velp,taur,&
               lquench,lmetares,resetcv,ekinc,ekinp)
       ENDIF
    ENDIF
    
    IF (cntl%mimic) THEN
       CALL mimic_sum_forces(fion)
    END IF

    if (biCanonicalEnsembleDo) then 
      call CpmdEnergiesGradients(bicanonicalCpmdConfig, ener_com%etot)
    end if 

    CALL freqs(crge%n,.TRUE.)
    ! Check orthogonality condition for wavefunction velocities
    ! SET UP BROKEN SYMMETRY STATE
    IF (cntl%bsymm)THEN
       bsclcs=1
       CALL setbsstate
    ENDIF
    CALL rortv(c0,cm,c2,sc0,gamy,nstate)
    ! SET UP HIGH SPIN STATE
    IF (cntl%bsymm)THEN
       bsclcs=2
       CALL setbsstate
       CALL rortv(c0(:,:,2),cm(1,1,2),c2(:,:,2),sc0(1,1,2),&
            gamy(1,2),nstate)
    ENDIF
    ! Initialize thermostats
    IF (paral%parent) THEN
       itemp=irec(irec_nop1)+irec(irec_nop2)+irec(irec_nop3)&
            +irec(irec_nop4)
       IF (cntl%tnosep .AND. itemp.EQ.0) CALL nospinit(ipwalk)
       IF (cntl%tnosee .AND. irec(irec_noe) .EQ.0) THEN
          ! THERMOSTAT FOR BROKEN SYMMETRY STATE WAVEFUNCTIONS
          CALL noseinit(ipwalk)
          ! THERMOSTAT FOR HIGH SPIN STATE WAVEFUNCTIONS
          IF (cntl%bsymm)CALL noseinit(2)
       ENDIF
       CALL wrgeof(tau0,fion)
       IF (paral%io_parent)&
            CALL fileopen(3,filen,fo_app+fo_verb,ferror)
       IF (cntl%bsymm) THEN
          ! OPEN OUTPUT DATA FILE FOR BROKEN SYMMETRY
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
               CALL fileopen(77,filenbs,fo_app+fo_verb,ferror)
       ENDIF
    ENDIF

    ! INITIALIZE GLE
    CALL gle_init(tau0,velp,rmass%pma)

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
    ! EVALUATE KINETIC ENERGY FOR BROKEN SYMMETRY STATE WAVEFUNCTIONS
    IF (cntl%bsymm)THEN
       CALL rekine(cm,nstate,ekinc)
       ekinc_bs=ekinc
       ! EVALUATE KINETIC ENERGY FOR HIGH SPIN STATE WAVEFUNCTIONS
       CALL rekine(cm(1,1,2),nstate,ekinc)
       ekinc_hs=ekinc
    ENDIF
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
    ! ==================================================================
    ! ==          THE BASIC LOOP FOR MOLECULAR DYNAMICS               ==
    ! ==                 USING VELOCITY VERLET                        ==
    ! ==================================================================
    ALLOCATE(center(4,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    DO infi=1,cnti%nomore
       CALL mp_sync(parai%cp_grp) ! allgrp-cp_grp bugfix
       time1=m_walltime()
       iteropt%nfi=iteropt%nfi+1
       ropt_mod%prteig=cprint%tprint.AND.MOD(iteropt%nfi-1,cprint%iprint_step).EQ.0
       cntl%caldip=cntl%tdipd.AND.MOD(iteropt%nfi-1,cnti%npdip).EQ.0
       comvl%subcom=comvl%tsubcom.AND.MOD(iteropt%nfi-1,comvl%ncomv).EQ.0
       comvl%subrot=comvl%tsubrot.AND.MOD(iteropt%nfi-1,comvl%nrotv).EQ.0
       IF (hubbu%pfrqom.gt.0) THEN
          hubbu%tpom=MOD(iteropt%nfi-1,hubbu%pfrqom).EQ.0
       ELSE
          hubbu%tpom=.False.
       ENDIF
       IF (.NOT.paral%parent) ropt_mod%prteig=.FALSE.
       ropt_mod%engpri=cprint%tprint.AND.MOD(iteropt%nfi-1,cprint%iprint_step).EQ.0
       tstrng=cntl%tpmin.AND.MOD(iteropt%nfi,cnti%nomore).EQ.0
       IF (cntl%mimic) THEN
         CALL mimic_switch_dim(go_qm=.FALSE.)
       ENDIF
       ! ANNEALING
       bsclcs=1
       CALL berendsen(velp,cm(1,1,1),nstate,dummy,ekinc,0.0_real_8)
       CALL anneal(velp,cm(1,1,1),nstate,dummy)
       IF (cntl%bsymm)THEN
          bsclcs=2
          CALL anneal(velp,cm(1,1,2),nstate,dummy)
          CALL berendsen(velp,cm(1,1,2),nstate,dummy,ekinc,0.0_real_8)
       ENDIF
       ! SUBTRACT CENTER OF MASS VELOCITY
       IF (paral%parent.AND.comvl%subcom) CALL comvel(velp,vcmio,.TRUE.)
       ! SUBTRACT ROTATION AROUND CENTER OF MASS
       IF (paral%parent.AND.comvl%subrot) CALL rotvel(tau0,velp,lmio,tauio,.TRUE.)
       ! UPDATE NOSE THERMOSTATS
       ! -FOR BROKEN SYMMETRY
       bsclcs=1
       CALL noseup(velp,cm,nstate,ipwalk)
       ! FIRST HALF OF GLE EVOLUTION
       IF (glepar%gle_mode.GT.0) CALL gle_step(tau0,velp,rmass%pma)

       ! -FOR HIGH SPIN 
       IF (cntl%bsymm)THEN
          bsclcs=2
          CALL noseup(velp,cm(1,1,2),nstate,2)
       ENDIF
       ! UPDATE VELOCITIES
       IF (paral%parent) CALL velupi(velp,fion,1)
       ! UPDATE BROKEN SYMMETRY WF VELOCITIES
       CALL velupa(c0,cm,c2,nstate,1)
       ! UPDATE HIGH SPIN WF VELOCITIES
       IF (cntl%bsymm)CALL velupa(c0(:,:,2),cm(1,1,2),c2(:,:,2),nstate,1)
#if defined (__QMECHCOUPL)
       IF (paral%parent .AND. qmmech) THEN
          CALL mm_cpmd_velup(cntr%delt_ions)
       ENDIF
#endif
       ! UPDATE POSITIONS
       ener_com%ecnstr = 0.0_real_8
       ener_com%erestr = 0.0_real_8
       IF (paral%parent) THEN
          CALL posupi(tau0,taup,velp)
       ENDIF

       IF (paral%parent) THEN
          IF (cotc0%mcnstr.NE.0) THEN
             IF (cntl%new_constraints) THEN
                CALL do_shake(tau0, taup, velp)
             ELSE 
                CALL cpmdshake(tau0,taup,velp)
             ENDIF
          ENDIF
#if defined (__QMECHCOUPL)
          IF (qmmech) THEN
             CALL mm_cpmd_update_links(taup, ions0%na, ions1%nsp, maxsys%nax, ions1%nat)
             CALL mm_cpmd_posup(cntr%delt_ions)
          ENDIF
#endif
       ENDIF

       CALL mp_bcast(taup,SIZE(taup),parai%io_source,parai%cp_grp) ! allgrp-cp_grp bugfix
       IF (cntl%mimic) THEN
          CALL mimic_update_coords(taup, c0, cm, nstate, ncpw%ngw, inyh)
          mimic_control%update_potential = .TRUE.
          IF (paral%io_parent) THEN
             CALL mimic_ifc_send_coordinates()
             IF (mimic_control%tot_scf_energy) THEN
                CALL mimic_ifc_collect_energies()
             END IF
          END IF
          IF (mimic_control%long_range_coupling.AND.(mod(iteropt%nfi-1,mimic_control%update_sorting).EQ.0)) THEN
             CALL mimic_ifc_sort_fragments()
          END IF
          CALL mimic_switch_dim(go_qm=.TRUE.)
       ENDIF
       CALL phfac(taup)
       IF (corel%tinlc) CALL copot(rhoe,psi,ropt_mod%calste)
       ! BROKEN SYMMETRY WF UPDATED USING VELOCITY VERLET
       bsclcs=1
       IF (cntl%bsymm)CALL setbsstate
       CALL posupa(c0,cm,c2,gamx,nstate)
       ! HIGH SPIN WF UPDATED USING VELOCITY VERLET
       IF (cntl%bsymm)THEN
          bsclcs=2
          CALL setbsstate
          CALL posupa(c0(:,:,2),cm(1,1,2),c2(:,:,2),gamx(1,2),&
               nstate)
       ENDIF
       ! ..Dipole moment
       ! !!! NOT ACTIVATED FOR BROKEN SYMMETRY !!!
       IF (vdwl%vdwd) THEN
          IF (cntl%tpath) THEN
             ipx=ipcurr-np_low+1
          ELSE
             ipx=1
          ENDIF
          vdwwfl%twannup=twannupx(ipx)
       ENDIF
       IF ((cntl%caldip.OR.(vdwl%vdwd.AND.vdwwfl%twannup)).AND.(.NOT.cntl%bsymm)) THEN
          CALL ddipo(taup,c0(:,:,1),cm(:,:,1),c2(:,:,1),sc0,nstate,center)
          CALL wannier_print(iteropt%nfi,c0(:,:,1),taup,nstate,psi(:,1),center)
       ENDIF
       ! CALCULATE THE FORCES
       ropt_mod%calste=cntl%tpres.AND.MOD(iteropt%nfi-1,cnti%npres).EQ.0
       ! FOR BROKEN SYMMETRY STATE
       bsclcs=1
       IF (cntl%bsymm)CALL setbsstate
       CALL forcedr(c0(:,:,1),c2(:,:,1),sc0(:,:,1),rhoe,psi,taup,fion,eigv,&
            nstate,1,.FALSE.,.TRUE.)
       IF (cntl%bsymm) THEN
          spd_bs=chrg%csums
          spda_bs=chrg%csumsabs
          IF (paral%parent)THEN
             etotbs=ener_com%etot
             CALL dcopy(3*maxsys%nax*maxsys%nsx,fion,1,fnbs,1)
          ENDIF
          ! FOR HIGH SPIN STATE
          bsclcs=2
          CALL setbsstate
          CALL forcedr(c0(:,:,2),c2(:,:,2),sc0(:,:,2),rhoe,psi,taup,&
               fion,eigv(1,2),nstate,1,.FALSE.,.TRUE.)
          spd_hs=chrg%csums
          spda_hs=chrg%csumsabs
          ! ENERGY AND FORCES FOR LOW SPIN STATE 
          IF (paral%parent)THEN
             etoths=ener_com%etot
             ener_com%etot = (1.0_real_8 + cnstwgt) * etotbs - cnstwgt * etoths
             CALL lsforce(fnbs,fion)
          ENDIF
       ENDIF
       ! ==================================================================
       ! Damped Dynamics
       CALL dampdyn(velp,fion,cm(1,1,1),c2(:,:,1),nstate,dummy,dummy)
       IF (cntl%bsymm)THEN
          bsclcs=2
          CALL dampdyn(velp,fion,cm(1,1,2),c2(:,:,2),nstate,dummy,dummy)
       ENDIF
       ! ==================================================================
       ! Meta Dynamics of Collective Variables

       IF (lmeta%lcolvardyn) THEN
          lquench = .FALSE.
          lmetares= .FALSE.

          IF (lmeta%lextlagrange) THEN
             ! Metadynamics with Extended Lagrangian
             IF (lmeta%tmulti) THEN
                CALL meta_ext_mul(taup,velp,taur,&
                     lquench,lmetares,resetcv,ekinc,ekinp)
             ELSE IF (tmw)THEN
                CALL meta_extlagr_mw(taup,velp,taur,&
                     lquench,lmetares,resetcv,ekinc,ekinp)
             ELSE
                CALL meta_extlagr(taup,velp,taur,&
                     lquench,lmetares,resetcv,ekinc,ekinp)
             ENDIF
          ELSE
             ! Time dipendent potential applied directly on the Collective Variables
             IF (tmw)THEN
                CALL meta_colvar_mw(taup,velp,fion,taur,&
                     lquench,lmetares,ekinc,ekinp)
             ELSE
                CALL meta_colvar(taup,velp,fion,taur,&
                     lquench,lmetares,ekinc,ekinp)
             END IF
          ENDIF

          IF ((rmeta%tolkin.GT.0.0_real_8 .AND. ekinc.GT.rmeta%tolkin)&
               .OR. lquench) THEN
             ! Check for Quench on the BO surface
             bsclcs=1
             IF (cntl%bsymm)CALL setbsstate
             cntl%quenchb=.TRUE.
             CALL quenbo(c0(:,:,1),c2,sc0,taur,rhoe,psi)
             IF (cntl%bsymm)THEN
                bsclcs=2
                CALL setbsstate
                CALL quenbo(c0(:,:,2),c2(:,:,2),sc0(1,1,2),taur,rhoe,&
                     psi)
             ENDIF
             CALL zeroing(cm(:,:,1:bsfac))!,nstate*nkpt%ngwk*bsfac)
             resetcv=.TRUE.
             cntl%quenchb=.FALSE.
          ENDIF
          IF (tmw) CALL mp_sync(supergroup)
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

       ! ==================================================================
       ! IF (lmeta%lcolvardyn .AND. paral%parent) THEN ! bugfix
       IF (lmeta%lcolvardyn .AND. paral%io_parent) THEN

          ! Additional Contribution to FION due to the Metadynamics
          ! (from coupling pot.if extended Lagrangian, directly from V(S,t) if not)
          DO is = 1,ions1%nsp
             DO ia = 1,ions0%na(is)
                fion(1,ia,is) = fion(1,ia,is) + fhills(1,ia,is)
                fion(2,ia,is) = fion(2,ia,is) + fhills(2,ia,is)
                fion(3,ia,is) = fion(3,ia,is) + fhills(3,ia,is)
             ENDDO
          ENDDO
       ENDIF

       IF (cntl%mimic) THEN
          CALL mimic_sum_forces(fion)
       END IF

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

       if (biCanonicalEnsembleDo) then 
         call CpmdEnergiesGradients(bicanonicalCpmdConfig, ener_com%etot)
       end if  

       IF (cntl%mimic) THEN
         CALL mimic_switch_dim(go_qm=.FALSE.)
       ENDIF

       ! FINAL UPDATE FOR VELOCITIES
       ener_com%ecnstr = 0.0_real_8
       ener_com%erestr = 0.0_real_8
       IF (paral%parent) THEN
          CALL velupi(velp,fion,1)
#if defined (__QMECHCOUPL)
          IF (qmmech) THEN
             CALL mm_cpmd_velup(cntr%delt_ions)
          ENDIF
#endif
          IF (isos1%twall) CALL box_boundary(taup,velp)
          IF (cotc0%mcnstr.NE.0) THEN
             IF (cntl%new_constraints) THEN
                CALL do_rattle(taup, velp)
             ELSE
                CALL rattle(taup,velp)
             END IF
          END IF
       ENDIF
       ! UPDATE BROKEN_SYMMETRY_WAVEFUNCTION_VELOCITIES
       bsclcs=1
       IF (cntl%bsymm)CALL setbsstate
       CALL velupa(c0,cm,c2,nstate,1)
       CALL rortv(c0,cm,c2,sc0,gamy,nstate)
       ! UPDATE HIGH_SPIN_WAVEFUNCTION_VELOCITIES
       IF (cntl%bsymm)THEN
          bsclcs=2
          CALL setbsstate
          CALL velupa(c0(:,:,2),cm(1,1,2),c2(:,:,2),nstate,1)
          CALL rortv(c0(:,:,2),cm(1,1,2),c2(:,:,2),sc0(1,1,2),&
               gamy(1,2),nstate)
       ENDIF

       IF (paral%parent.AND..NOT.cntl%mimic) CALL geofile(taup,velp,'WRITE')
       ! COMPUTE THE IONIC TEMPERATURE TEMPP
       IF (paral%parent) THEN
          CALL ekinpp(ekinp,velp)
          IF (lmeta%lextlagrange.AND. ltcglobal) THEN
             CALL ekincv_global(ek_cv)
             tempp=(ek_cv+ekinp)*factem*2._real_8/(glib+REAL(ncolvar,kind=real_8))
          ELSE
             tempp=ekinp*factem*2._real_8/glib
          ENDIF
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
       ! IONIC TEMPERATURE CONTROL
       IF (paral%parent) CALL rscvp(temp1,temp2,tempp,velp)
       ! SUBTRACT ROTATION AROUND CENTER OF MASS
       IF (paral%parent.AND.comvl%subrot) CALL rotvel(tau0,velp,lmio,tauio,.FALSE.)
       ! SUBTRACT CENTER OF MASS VELOCITY
       IF (paral%parent.AND.comvl%subcom) CALL comvel(velp,vcmio,.FALSE.)
#if defined (__QMECHCOUPL)
       ! ---------------------------------------------------------------
       ! MM Temperature control - annealing schedule not yet implemented
       ! ---------------------------------------------------------------
       IF (paral%parent .AND. qmmech) THEN
          CALL mm_cpmd_temp_control(temp1,temp2,mm_temp,cntl%tcp)
          CALL mm_cpmd_ekin(mm_ekin,mm_temp)
       ENDIF
#endif
       ! SECOND HALF OF GLE EVOLUTION
       IF (glepar%gle_mode.GT.0) CALL gle_step(tau0,velp,rmass%pma)

       ! UPDATE NOSE THERMOSTATS
       IF (cntl%tnosee.OR.cntl%tc) THEN
          CALL rekine(cm,nstate,ekinc)
          IF (cntl%bsymm)THEN
             ekinc_bs=ekinc
             CALL rekine(cm(1,1,2),nstate,ekinc)
             ekinc_hs=ekinc
          ENDIF
       ENDIF
       ! UPDATE FOR BROKEN SYMMETRY 
       bsclcs=1
       CALL noseup(velp,cm,nstate,ipwalk)
       CALL berendsen(velp,cm(1,1,1),nstate,dummy,ekinc,0.0_real_8)
       ! UPDATE FOR HIGH SPIN 
       IF (cntl%bsymm)THEN
          bsclcs=2
          CALL noseup(velp,cm(1,1,2),nstate,2)
          CALL berendsen(velp,cm(1,1,2),nstate,dummy,ekinc,0.0_real_8)
       ENDIF
       ! ANNEALING
       ! -FOR BROKEN SYMMETRY
       bsclcs=1
       CALL anneal(velp,cm,nstate,dummy)
       ! FOR HIGH SPIN
       IF (cntl%bsymm)THEN
          bsclcs=2
          CALL anneal(velp,cm(1,1,2),nstate,dummy)
       ENDIF
       IF (paral%parent) THEN
          CALL ekinpp(ekinp,velp)
          IF (lmeta%lextlagrange.AND. ltcglobal) THEN
             CALL ekincv_global(ek_cv)
             tempp=(ek_cv+ekinp)*factem*2._real_8/(glib+REAL(ncolvar,kind=real_8))
          ELSE
             tempp=ekinp*factem*2._real_8/glib
          ENDIF
          IF (cntl%mimic) then
             CALL mimic_subsystem_temperatures(velp)
          ENDIF
       ENDIF
       ! RESCALE ELECTRONIC VELOCITIES
       IF ((.NOT.cntl%bsymm).AND.cntl%tc)&
            CALL rscve(ekin1,ekin2,ekinc,cntr%ekinw,cm,nstate,ncpw%ngw)
       ! RESCALE BROKEN SYMMETRY AND HIGH SPIN ELECTRONIC VELOCITIES
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
          ! PROJECTED KINETIC ENERGY OF ELECTRONS
          ekinc=(scalhs*ekinc_hs)+(scalbs*ekinc_bs)
       ENDIF
       ! ENERGY OF THE NOSE THERMOSTATS
       IF (paral%parent) THEN
          bsclcs=1
          CALL noseng(iteropt%nfi,velp,enose,enosp,dummy(1),ipwalk)
          IF (cntl%bsymm)THEN
             enosbs=enose
             bsclcs=2
             CALL noseng(iteropt%nfi,velp,enose,dum1,dum2,2)
             enoshs=enose
             enose=(scalbs*enosbs)+(scalhs*enoshs)
          ENDIF
       ENDIF
       IF (paral%parent) THEN
          econs=ekinp+ener_com%etot+enose+enosp+ener_com%ecnstr+ener_com%erestr+ekincv+vharm+glepar%egle
          ! 
#if defined (__QMECHCOUPL)
          IF (qmmech) THEN
             econs=econs + mm_ekin
             IF (paral%io_parent)&
                  WRITE (6,'(50x, "TEMP: ",f10.0)') mm_temp
          ENDIF
#endif
          eham=econs+ekinc
          rmeta%eham_hill = eham
          time2=m_walltime()
          tcpu=(time2-time1)*0.001_real_8
          ! Printing data of Broken symmetry calculation
          IF (cntl%bsymm) THEN
             couplj = (etoths - etotbs) * rtsasb
             jwn=couplj*autocm
             IF (paral%io_parent)&
                  WRITE(77,'(I10,11F16.8)')&
                  iteropt%nfi,ekinc_bs,ekinc_hs,ener_com%etot,etotbs,etoths,spd_bs,&
                  spd_hs,spda_bs,spda_hs,couplj,jwn
          ENDIF
          ! PRINTOUT the evolution of the accumulators every time step
          CALL wrprint_md(eigv,crge%f,ener_com%amu,nstate,taup,fion,&
               ekinc,tempp,ener_com%etot,econs,eham,disa,&
               tcpu,.FALSE.,iteropt%nfi,infi)
          ! UPDATE ACCUMULATORS
          CALL pacca(ekinc,tempp,ener_com%etot,econs,eham,enose,enosp,ener_com%ecnstr,&
               ener_com%erestr,disa,tcpu,iteropt%nfi,1)
          ! Store ionic coordinates and velocities for statistics
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
       ENDIF
       if (biCanonicalEnsembleDo)&
           call PrintEnergies(bicanonicalCpmdConfig, iteropt%nfi)

       IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)
       IF (tmw) CALL testex_mw(soft_com%exsoft)
       if (biCanonicalEnsembleDo)&
         call bicanonicalTestExit(soft_com%exsoft)
       IF (tstrng) THEN
          soft_com%exsoft=.TRUE.
          soft_com%exnomore=.TRUE.
       ELSE IF (infi.EQ.cnti%nomore) THEN
          soft_com%exsoft=.TRUE.
       ENDIF
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
             ELSE IF (tmw)THEN
                CALL meta_extlagr_mw(taup,velp,taur,&
                     lquench,lmetares,resetcv,ekinc,ekinp)
             ELSE
                CALL meta_extlagr(taup,velp,taur,&
                     lquench,lmetares,resetcv,ekinc,ekinp)
             ENDIF
          ELSE
             ! Time dependent potential applied directly on the Collective Variables
             IF (tmw)THEN
                CALL meta_colvar_mw(taup,velp,fion,taur,&
                     lquench,lmetares,ekinc,ekinp)
             ELSE
                CALL meta_colvar(taup,velp,fion,taur,&
                     lquench,lmetares,ekinc,ekinp)
             END IF
          ENDIF
       ENDIF
#if defined (__QMECHCOUPL)
       ! -----------------------
       ! Write QMMM restart file
       ! -----------------------
       IF (paral%parent.AND.qmmech) THEN
          IF (MOD(iteropt%nfi,store1%istore).EQ.0.OR.infi.EQ.cnti%nomore.OR.soft_com%exsoft) THEN
             CALL mm_cpmd_write_restart
          ENDIF
       ENDIF
#endif
       ! Calculate EPR and EFG properties (not active for broken symmetry).
       IF (.NOT.cntl%bsymm) CALL epr_efg(rhoe,psi,iteropt%nfi)
       IF (cntl%mimic) THEN
          CALL mimic_switch_dim(go_qm=.FALSE.)
       ENDIF
       ! temperature ramping
       CALL tempramp(temp1,temp2)
       ! UPDATE IONIC POSITIONS
       CALL dcopy(3*maxsys%nax*maxsys%nsx,taup(1,1,1),1,tau0(1,1,1),1)
       IF (cntl%mimic) THEN
          CALL mimic_switch_dim(go_qm=.TRUE.)
       ENDIF
       ! STOP THE RUN IF THE USER HAS SET THE SIGNAL 30
       IF (soft_com%exsoft) GOTO 100
       ! ==================================================================
       ! ==     END OF MAIN LOOP                                         ==
       ! ==================================================================
    ENDDO
100 CONTINUE
    IF (cntl%mimic) THEN
       CALL mimic_update_coords(taup, c0, cm, nstate, ncpw%ngw, inyh)
       IF (paral%io_parent) THEN
          CALL mimic_ifc_send_coordinates()
       END IF
       IF (mimic_control%long_range_coupling) THEN
          CALL mimic_ifc_sort_fragments()
       END IF
    END IF
    IF (wannl%twann) THEN
       CALL ddipo(taup,c0(:,:,1),cm(:,:,1),c2(:,:,1),sc0,nstate,center)
       CALL forcedr(c0(:,:,1),c2(:,:,1),sc0(:,:,1),rhoe,psi,taup,fion,eigv,&
            nstate,1,.FALSE.,.TRUE.)
       CALL wc_dos(c0,c2,nstate,center)
    ENDIF
    DEALLOCATE(center,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    IF (rout1%rhoout.AND.(rout1%nrhoout.LE.0))THEN
       bsclcs=1
       CALL rhopri(c0,tau0,rhoe,psi(:,1),nstate,nkpt%nkpnt)
       IF (cntl%bsymm)THEN
          bsclcs=2
          CALL rhopri(c0(:,:,2:),tau0,rhoe,psi(:,1),nstate,nkpt%nkpnt)
       ENDIF
    ENDIF
    ! PRINT ACCUMULATOR
    IF (paral%parent) CALL pacca(ekinc,tempp,ener_com%etot,econs,eham,enose,enosp,&
         ener_com%ecnstr,ener_com%erestr,disa,tcpu,iteropt%nfi,2)
    bsclcs=1
    IF (cntl%bsymm)CALL setbsstate
    CALL proja(c0,c2,sc0,nstate,cnti%iproj)
    IF (cntl%bsymm)THEN
       bsclcs=2
       CALL setbsstate
       CALL proja(c0(:,:,2),c2(:,:,2),sc0(1,1,2),&
            nstate,cnti%iproj)
    ENDIF
    bsclcs=1
    IF (cntl%bsymm)CALL setbsstate
    CALL csize(c2,crge%n,gemax,cnorm)
    IF (cntl%bsymm)THEN
       bsclcs=2
       CALL setbsstate
       CALL csize(c2(:,:,2),crge%n,gemax,cnorm)
    ENDIF
    IF (paral%parent) CALL gsize(fion,gnmax,gnorm)
    IF (paral%parent) THEN
       IF (.NOT.cntl%bsymm)THEN
          CALL finalp(tau0,fion,velp,eigv)
       ELSE
          CALL wrgeof(tau0,fion)
       ENDIF
    ENDIF
    ! 
    IF (cntl%tsampl) THEN
       CALL sample_go
       GOTO 99999
    ENDIF
10000 CONTINUE
    ! 
    DEALLOCATE(taup,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(fion,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(taui,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(taur,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (cntl%nonort) THEN
       DEALLOCATE(eigv,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(eigm,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    IF (paral%parent.AND.paral%io_parent) CALL fileclose(3)
    DEALLOCATE(rhoe,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(psi,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(scr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (comvl%tsubrot) DEALLOCATE(tauio,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (cntl%bsymm) DEALLOCATE(fnbs,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (lmeta%lcolvardyn) THEN
       DEALLOCATE(fhills,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! NN: CLOSING THE FILE BS_ENG
    IF ((cntl%bsymm.AND.paral%parent).AND.paral%io_parent)&
         CALL fileclose(77)
    IF (cntl%mimic) THEN
       CALL mimic_revert_dim()
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
  END SUBROUTINE mdmain
  ! 
  ! ==================================================================
  SUBROUTINE give_scr_mdmain(lmdmain,tag)
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
    IF (cntl%tdipd.OR.vdwl%vdwd) CALL give_scr_ddipo(lddipo,tag)
    CALL give_scr_forcedr(lforcedr,tag,nstate,.FALSE.,.TRUE.)
    CALL give_scr_rortv(lrortv,tag,nstate)
    CALL give_scr_posupa(lposupa,tag,nstate)
    IF (rout1%rhoout) CALL give_scr_rhopri(lrhopri,tag,nstate)
    CALL give_scr_meta_extlagr(lmtd,tag)
    lmdmain=MAX(lcopot,lortho,lquenbo,ldeort,lforcedr,&
         lrortv,lposupa,lrhopri,lddipo,linitrun,lmtd)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_mdmain
  ! ==================================================================

END MODULE mdmain_utils
