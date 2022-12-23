MODULE pi_npt_cpmd_utils
  USE anneal_utils,                    ONLY: anneal,&
                                             berendsen
  USE cnst,                            ONLY: factem
  USE cnstpr_utils,                    ONLY: cnstpr
  USE comvel_utils,                    ONLY: comvel
  USE comvelmod,                       ONLY: comvl
  USE coor,                            ONLY: fion,&
                                             tau0,&
                                             taup,&
                                             velp
  USE copot_utils,                     ONLY: copot
  USE cotr,                            ONLY: cotc0
  USE ddipo_utils,                     ONLY: ddipo
  USE deort_utils,                     ONLY: deort
  USE detdof_utils,                    ONLY: detdof,&
                                             qmdof
  USE dispp_utils,                     ONLY: dispp
  USE dynit_utils,                     ONLY: dynit
  USE ekinpp_utils,                    ONLY: ekinpp,&
                                             s_ekinpp
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: chrg,&
                                             ener_com
  USE error_handling,                  ONLY: stopgm
  USE evirial_utils,                   ONLY: evirial
  USE fharm_utils,                     ONLY: fharm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_app,&
                                             fo_mark,&
                                             fo_verb
  USE filnmod,                         ONLY: filbod,&
                                             filn
  USE forcedr_driver,                  ONLY: forcedr
  USE forcep_utils,                    ONLY: rhoe_psi_size
  USE freqs_utils,                     ONLY: freqs
  USE geofile_utils,                   ONLY: geofile
  USE getcor_utils,                    ONLY: getcor
  USE getfnm_utils,                    ONLY: getfnm
  USE getfu_utils,                     ONLY: getfu
  USE getgyr_utils,                    ONLY: getgyration
  USE global_utils,                    ONLY: global
  USE ions,                            ONLY: ions1
  USE kinds,                           ONLY: real_8
  USE localize_utils,                  ONLY: localize2
  USE machine,                         ONLY: m_walltime
  USE metr,                            ONLY: eps,&
                                             metr_com,&
                                             veps
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sync
  USE newcell_utils,                   ONLY: newcell
  USE nlcc,                            ONLY: corel,&
                                             vnlcc,&
                                             vnlt
  USE nose,                            ONLY: etap1,&
                                             etap1dot,&
                                             etapm,&
                                             etapmdot,&
                                             glib,&
                                             loct,&
                                             nchx,&
                                             nosl,&
                                             ncdof,&
                                             gkt,&
                                             gkt1,&
                                             tnosepc
  USE noscinit_utils,                  ONLY: noscinit
  USE noseinit_utils,                  ONLY: noseinit
  USE noseng_utils,                    ONLY: noseng
  USE nosepa_utils,                    ONLY: nosepa
  USE noseup_utils,                    ONLY: noseup
  USE nospinit_utils,                  ONLY: nospinit
  USE ortho_utils,                     ONLY: ortho
  USE parac,                           ONLY: parai,&
                                             paral
  USE phfac_utils,                     ONLY: phfac
  USE pimd,                            ONLY: &
       csumgv, csumrv, egcv, ehar, eharv, ehtv, ekinv, enlv, epseuv, eselfv, &
       esrv, etotv, excv, fionks, grandparent, ipcurr, maxnp, np_high, &
       np_local, np_low, pimd1, pimd3, supergroup, parentgroup, supersource
  USE pi_md_utils,                     ONLY: pi_ener,&
                                             disten
  USE pinmtrans_utils,                 ONLY: pinmtrans
  USE pi_stress_utils,                 ONLY: pi_stress_vir,&
                                             wr_pi_stress
  USE pi_wf_utils,                     ONLY: initrep
  USE posupa_utils,                    ONLY: posupa
  USE posupi_utils,                    ONLY: posupi,&
                                             posupih,&
                                             posupih_iso
  USE prcp,                            ONLY: prcp_com,&
                                             prcpl
  USE prcpmd_utils,                    ONLY: give_scr_prcpmd
  USE printave_utils,                  ONLY: pacce
  USE printp_utils,                    ONLY: printp
  USE prmem_utils,                     ONLY: prmem
  USE prtgyr_utils,                    ONLY: prtgyr
  USE pslo,                            ONLY: pslo_com
  USE puttau_utils,                    ONLY: taucl
  USE quenbo_utils,                    ONLY: quenbo
  USE ranp_utils,                      ONLY: ranp
  USE rattle_utils,                    ONLY: rattle
  USE readsr_utils,                    ONLY: xstring
  USE rekine_utils,                    ONLY: rekine
  USE reshaper,                        ONLY: reshape_inplace
  USE rhopri_utils,                    ONLY: rhopri
  USE rinitwf_driver,                  ONLY: rinitwf
  USE rinvel_utils,                    ONLY: rinvel,&
                                             rvscal,&
                                             s_rinvel
  USE rnlsm_utils,                     ONLY: rnlsm
  USE ropt,                            ONLY: infi,&
                                             iteropt,&
                                             ropt_mod
  USE rortv_utils,                     ONLY: rortv
  USE rotvel_utils,                    ONLY: rotvel
  USE rscve_utils,                     ONLY: rscve
  USE rscvp_utils,                     ONLY: rscvp
  USE rv30_utils,                      ONLY: zhrwf
  USE setirec_utils,                   ONLY: read_irec,&
                                             write_irec
  USE setsc_utils,                     ONLY: ihmat
  USE shake_utils,                     ONLY: cpmdshake
  USE shock,                           ONLY: shock1
  USE soft,                            ONLY: soft_com
  USE spin,                            ONLY: spin_mod,&
                                             clsd
  USE stagetrans_utils,                ONLY: stagetrans
  USE store_types,                     ONLY: &
       cprint, iprint_coor, iprint_force, irec_nop1, irec_nop2, irec_nop3, &
       irec_nop4, irec_wf, restart1, rout1, irec_noc
  USE strs,                            ONLY: paiu
  USE system,                          ONLY: &
       acc, cnti, cntl, cntr, maxsys, nacc, nacx, ncpw, restf, parm
  USE testex_utils,                    ONLY: testex,&
                                             testex_mw
  USE teststore_utils,                 ONLY: teststore
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE totstr_utils,                    ONLY: totstr
  USE vdw_utils,                       ONLY: rwannier
  USE vdwcmod,                         ONLY: &
       nfragx, nwfcx, rwann, swann, tauref, trwanncx, twannupx, vdwl, vdwwfl
  USE velupa_utils,                    ONLY: velupa
  USE velupi_utils,                    ONLY: s_velupi,&
                                             velupc
  USE wannier_print_utils,             ONLY: wannier_print
  USE wr_temps_utils,                  ONLY: wr_temps
  USE wrener_utils,                    ONLY: wrener
  USE wrgeo_utils,                     ONLY: wrgeo,&
                                             wrgeof
  USE wv30_utils,                      ONLY: zhwwf
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: pi_npt_cpmd

CONTAINS

  ! ==================================================================
  SUBROUTINE pi_npt_cpmd(c0,cm,c2,sc0,gamx,gamy)
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8)                          :: c0(:,:,:), &
                                                cm(ncpw%ngw,crge%n,*), &
                                                c2(ncpw%ngw,crge%n,*), &
                                                sc0(ncpw%ngw,*)
    REAL(real_8)                             :: gamx(crge%n*crge%n,*), &
                                                gamy(crge%n*crge%n,*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'pi_npt_cpmd'

    CHARACTER(len=100)                       :: filen
    CHARACTER(len=12)                        :: cflbod, cipnum
    CHARACTER(len=30)                        :: tag
    COMPLEX(real_8), ALLOCATABLE             :: psi(:,:)
    INTEGER :: i, i1, i2, ierr, il_psi_1d, il_psi_2d, il_rhoe_1d, il_rhoe_2d, &
      ip, ipx, isub, itemp, k, l, lscr, n1, n2, npx, nrepl, nwfc, j
    INTEGER, ALLOCATABLE                     :: irec(:,:)
    LOGICAL                                  :: ferror, reset_gkt,&
                                                pitmnose(2), pitnosec(2)
    REAL(real_8) :: accus(nacx,maxnp), d1, d2, disa(maxnp), dummies(1), &
      dummy, econs(maxnp), econsa, eham(maxnp), ehama, ekin1, ekin2, ekina, &
      ekinc(maxnp), ekincp, ekinh1, ekinh2, ekinh, ekinp, enose(maxnp), &
      enosp(maxnp), equant, etota, ff, glib_s, lmio(3), qpkinp, qvkinp, rnp, &
      tcpu, temp1, temp2, tempa, tempp(maxnp), time1, time2, vcmio(4), &
      enosc(maxnp), temph, pipaiu(3,3)
    REAL(real_8), ALLOCATABLE :: center(:,:), eigm(:,:), eigv(:,:), &
      fstage(:,:,:,:), rhoe(:,:), scr(:), stagep(:,:,:,:), taui(:,:,:,:), &
      tauio(:,:), taur(:,:,:), vstage(:,:,:,:), paiuall(:,:,:)
    REAL(real_8), POINTER                    :: pifion(:,:,:,:), &
                                                pitaup(:,:,:,:), &
                                                pivelp(:,:,:,:)
    REAL(real_8), EXTERNAL                   :: ddot

! Variables
! NEW
! NEW
! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    restf%nfnow=1
    time1 =m_walltime()
    npx=pimd3%np_total
    rnp=REAL(npx,kind=real_8)
    ! TAUP, FION and VELP are defined in coor.mod as 3D allocatables
    ALLOCATE(taup(3,maxsys%nax,maxsys%nsx*npx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(fion(3,maxsys%nax,maxsys%nsx*npx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    IF(ALLOCATED(velp)) THEN
       DEALLOCATE(velp,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ALLOCATE(velp(3,maxsys%nax,maxsys%nsx*npx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(taup)!,SIZE(taup))
    CALL zeroing(fion)!,SIZE(fion))
    CALL zeroing(velp)!,SIZE(velp))
    ALLOCATE(taui(3,maxsys%nax,maxsys%nsx,npx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(taur(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(fionks(3,maxsys%nax,maxsys%nsx,npx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(fionks)
    IF (pimd1%tstage.OR.pimd1%tpinm) THEN
       ALLOCATE(stagep(3,maxsys%nax,maxsys%nsx,npx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(vstage(3,maxsys%nax,maxsys%nsx,npx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(fstage(3,maxsys%nax,maxsys%nsx,npx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    IF (comvl%tsubrot) THEN
       ALLOCATE(tauio(3,ions1%nat),stat=ierr)
       IF(ierr/=0) CALL stopgm(proceduren,'allocation problem',&
            __LINE__,__FILE__)! TODO check dimensions
    ENDIF
    IF (pslo_com%tivan.OR.cntl%nonort) THEN
       ALLOCATE(eigv(crge%n,np_local),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(eigm(crge%n*crge%n,np_local),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ELSE
       ! Avoid not allocated Fortran runtime error
       ALLOCATE(eigv(1,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    nacc = 30
    CALL zeroing(accus)!,nacx*maxnp)
    ropt_mod%modens=.FALSE.
    ropt_mod%engpri=.FALSE.
    ropt_mod%calste=.TRUE.
    ALLOCATE(paiuall(3,3,npx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    IF (cntl%tharm.AND.paral%parent) THEN
       ff=crge%f(1,1)
       DO i=1,crge%n
          IF (ff.NE.crge%f(i,1)) THEN
             ! DM1
             ! DM            WRITE(6,'(A,A)') ' HARMONIC REFERENCE SYSTEM INTEGRATION',
             IF (grandparent.AND.paral%io_parent)&
                  WRITE(6,'(A,A)') ' HARMONIC REFERENCE SYSTEM INTEGRATION',&
                  ' ONLY POSSIBLE WITH EQUAL OCCUPATION NUMBERS'
             CALL stopgm('MDMAIN',' ',& 
                  __LINE__,__FILE__)
          ENDIF
       ENDDO
    ENDIF
    ALLOCATE(irec(100,npx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    IF (corel%tinlc) THEN
       ALLOCATE(vnlt(ncpw%nhg*2),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)

       ALLOCATE(vnlcc(ncpw%nhg,clsd%nlsd),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! NEW  EQUIVALENCE FION <-> PIFION 
    CALL reshape_inplace(fion, (/3,maxsys%nax,maxsys%nsx, npx/), pifion)
    CALL reshape_inplace(taup, (/3,maxsys%nax,maxsys%nsx, npx/), pitaup)
    CALL reshape_inplace(velp, (/3,maxsys%nax,maxsys%nsx, npx/), pivelp)
    ! ==--------------------------------------------------------------==
    ! SCR ALLOCATION AND PARTITION (SCRATCH ARRAY).
    CALL rhoe_psi_size(il_psi_1d=il_psi_1d,il_psi_2d=il_psi_2d,&
         il_rhoe_1d=il_rhoe_1d,il_rhoe_2d=il_rhoe_2d)
    ALLOCATE(rhoe(il_rhoe_1d,il_rhoe_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(psi(il_psi_1d,il_psi_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL give_scr_prcpmd(lscr,tag)
    ALLOCATE(scr(lscr),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ..TIME STEP FUNCTIONS
    CALL dynit(ekincp,ekin1,ekin2,temp1,temp2,ekinh1,ekinh2)
    ! ..PARAMETERS FOR THE NOSE-HOOVER THERMOSTATS
    ncdof = 9
    IF (prcpl%tisot.OR.prcpl%tzflex) ncdof = 1
    IF ((cntl%tnosee.OR.cntl%tnosep.OR.cntl%tnosec).AND.paral%parent) CALL nosepa(np_low,np_high)
    ! thermostat
    reset_gkt=.FALSE.
    IF (ALLOCATED(gkt1)) reset_gkt=.TRUE.
    pitmnose(:)=nosl%tmnose
    IF (nosl%tmnose.AND.pimd1%tcentro) pitmnose(1)=.FALSE. !massive -> single for ip=1
    ! barostat
    pitnosec(:)=cntl%tnosec
    IF (cntl%tnosec) pitnosec(2)=.FALSE. !switch off for ip>1
    ! ..INITIALIZATION OF WAVEFUNCTION AND COORDINATS
    CALL zeroing(acc)!,nacc)
    DO ipcurr=np_low,np_high
       nosl%tmnose=pitmnose(MIN(ipcurr,2))
       cntl%tnosec=pitnosec(MIN(ipcurr,2))
       CALL read_irec(irec(:,ipcurr))
       IF (restart1%restart) THEN
          ! ..Construct filenames
          cflbod='RESTART_'
          IF (paral%io_parent)&
               WRITE(cipnum,'(I4)') ipcurr
          CALL xstring(cflbod,n1,n2)
          CALL xstring(cipnum,i1,i2)
          filbod=cflbod(n1:n2)//cipnum(i1:i2)//'.1'
          IF (restart1%rlate) THEN
             filn=filbod
          ELSE
             filn=cflbod(n1:n2)//cipnum(i1:i2)
          ENDIF
          ! ..Read restart file
          ipx=ipcurr-np_low+1
          CALL zhrwf(1,irec(:,ipcurr),c0(:,:,ipx:ipx),cm(:,:,ipx),crge%n,eigv(:,ipx),&
               pitaup(:,:,:,ipcurr),&
               pivelp(:,:,:,ipcurr),taui(:,:,:,ipcurr),iteropt%nfi)
          CALL dcopy(nacc,acc,1,accus(1,ipcurr),1)
          IF (restart1%rgeo) THEN
             IF (paral%io_parent) CALL geofile(pitaup(:,:,:,ipcurr),pivelp(:,:,:,ipcurr),'READ')
             CALL mp_bcast(pitaup(:,:,:,ipcurr),3*maxsys%nax*maxsys%nsx,parai%io_source,parai%cp_grp)
             CALL mp_bcast(pivelp(:,:,:,ipcurr),3*maxsys%nax*maxsys%nsx,parai%io_source,parai%cp_grp)
             CALL mp_bcast(metr_com%ht,SIZE(metr_com%ht),parai%io_source,parai%cp_grp)
             CALL mp_bcast(metr_com%htvel,SIZE(metr_com%htvel),parai%io_source,parai%cp_grp)
             restart1%rcell=.TRUE.; restart1%rvel=.TRUE.
          ENDIF
          ! ..Randomization of the atomic coordinates
          IF (cntl%tranp.AND.paral%parent) CALL ranp(pitaup(:,:,:,ipcurr))
          ! ..Initialization of wavefunction
          IF (irec(irec_wf,ipcurr).EQ.0) THEN
             CALL rinitwf(c0(:,:,ipx:ipx),cm,sc0,crge%n,pitaup(:,:,:,ipcurr),&
                  taur,rhoe,psi)
             cntl%quenchb=.TRUE.
          ENDIF
          CALL phfac(pitaup(:,:,:,ipcurr))
          IF (corel%tinlc) CALL copot(rhoe,psi,.TRUE.)
          ! Orthogonalize electronic wavefunctions
          IF (irec(irec_wf,ipcurr).EQ.0.OR.cntl%trane.OR.(pslo_com%tivan.AND..NOT.cntl%nonort)) THEN
             IF (pslo_com%tivan) CALL rnlsm(c0(:,:,ipx),crge%n,1,1,.FALSE.)
             CALL ortho(crge%n,c0(:,:,ipx),c2)
          ENDIF
       ENDIF
    ENDDO
    ! ..Update cell
    IF (restart1%rcell) THEN
       ! ..Make sure all nodes have exactly the same cell 
       IF (paral%io_parent) THEN
          CALL mp_bcast(metr_com%ht,SIZE(metr_com%ht),supersource,parentgroup)
       ENDIF
       CALL mp_bcast(metr_com%ht,SIZE(metr_com%ht),parai%io_source,parai%cp_grp)
       DO i=1,3
          parm%a1(i) = metr_com%ht(1,i)
          parm%a2(i) = metr_com%ht(2,i)
          parm%a3(i) = metr_com%ht(3,i)
       ENDDO
       CALL ihmat(metr_com%ht,metr_com%htm1,parm%omega)
       CALL newcell
    ENDIF
    ! ..Initialize replica coordinates if requested
    IF (.NOT.restart1%restart) THEN
       IF (pimd1%tinit) THEN
          CALL initrep(pitaup)
          cntl%quenchb=.TRUE.
       ELSE
          CALL stopgm(procedureN,'replica coordinates not read/generated',&
             __LINE__,__FILE__)
       ENDIF
    ENDIF
    ! ..Global distribution of ionic positions
    CALL global(pitaup,3*maxsys%nax*maxsys%nsx)
    CALL mp_bcast(pitaup,SIZE(pitaup),parai%io_source,parai%cp_grp)
    IF (paral%parent) THEN
       IF (pimd1%tstage) THEN
          CALL stagetrans(pitaup,stagep,0)
          CALL dcopy(3*maxsys%nax*maxsys%nsx,stagep,1,tau0,1)
       ELSEIF (pimd1%tpinm) THEN
          CALL pinmtrans(pitaup,stagep,0)
          CALL dcopy(3*maxsys%nax*maxsys%nsx,stagep,1,tau0,1)
       ENDIF
       ! ..Setup for constraints and degrees of freedom
       CALL detdof(tau0,taur)
       CALL qmdof
    ELSE
       glib=1._real_8
    ENDIF
    DO ipcurr=np_low,np_high
       ipx=ipcurr-np_low+1
       CALL phfac(pitaup(:,:,:,ipcurr))
       ! ..Read WF centers & spread from restart file
       IF (vdwl%vdwd) THEN
          IF (paral%io_parent) THEN
             nwfc=crge%n
             vdwwfl%trwannc=trwanncx(ipx)
             CALL rwannier(nwfc,tauref(:,:,:,ipx),rwann(:,:,ipx),swann(:,ipx),&
                  vdwwfl%trwannc)
             IF (.NOT.vdwwfl%trwannc) THEN
                CALL dcopy(3*maxsys%nax*maxsys%nsx,pitaup(1,1,1,ipcurr),1,&
                     tauref(1,1,1,ipx),1)
             ENDIF
          ENDIF
          CALL mp_bcast(tauref(:,:,:,ipx),3*maxsys%nax*maxsys%nsx,&
               parai%io_source,parai%cp_grp)
          CALL mp_bcast(rwann(:,:,ipx),3*nwfcx*nfragx,parai%io_source,parai%cp_grp)
          CALL mp_bcast(swann(:,ipx),nwfcx*nfragx,parai%io_source,parai%cp_grp)
          CALL mp_bcast(vdwwfl%trwannc,parai%io_source,parai%cp_grp)
          trwanncx(ipx)=vdwwfl%trwannc
       ENDIF
       ! ..Initialization of wavefunction
       IF (.NOT.restart1%restart) THEN
          CALL rinitwf(c0(:,:,ipx:ipx),cm,sc0,crge%n,pitaup(:,:,:,ipcurr),&
               taur,rhoe,psi)
       ENDIF
       ! ..Quenching to BO surface
       IF (cntl%quenchb) THEN
          CALL dcopy(3*maxsys%nax*maxsys%nsx,pitaup(1,1,1,ipcurr),1,tau0,1)
          CALL quenbo(c0(:,:,ipx),c2(1,1,ipx),sc0,taur,rhoe,psi)
       ENDIF
       IF (pslo_com%tivan) THEN
          IF (cntl%tlsd) THEN
             CALL deort(ncpw%ngw,spin_mod%nsup,eigm(1,ipx),eigv(1,ipx),&
                  c0(:,1:spin_mod%nsup,ipx),sc0(1,1))
             CALL deort(ncpw%ngw,spin_mod%nsdown,eigm(1,ipx),eigv(1,ipx),&
                  c0(:,spin_mod%nsup+1:spin_mod%nsup+spin_mod%nsdown,ipx),sc0(1,spin_mod%nsup+1))
          ELSE
             CALL deort(ncpw%ngw,crge%n,eigm(1,ipx),eigv(1,ipx),&
                  c0(:,:,ipx),sc0)
          ENDIF
       ENDIF
       ! ..Initialize velocities
       IF (.NOT.restart1%rvel) THEN
          ener_com%ecnstr = 0.0_real_8
          ener_com%erestr = 0.0_real_8
          IF (pimd1%tpinm.OR.pimd1%tstage) THEN
             CALL s_rinvel(vstage(:,:,:,ipcurr),cm(1,1,ipx),crge%n,ipcurr)
             IF (paral%parent.AND.ipcurr.EQ.1) CALL rattle(stagep(:,:,:,1),vstage(:,:,:,1))
             IF (ipcurr.EQ.1) CALL rvscal(vstage(:,:,:,ipcurr))
          ELSE
             CALL rinvel(pivelp(:,:,:,ipcurr),cm(1,1,ipx),crge%n)
          ENDIF
       ENDIF
       IF (cntl%quenchp) THEN
          IF (pimd1%tpinm.OR.pimd1%tstage) THEN
             CALL zeroing(vstage(:,:,:,ipcurr))!,3*maxsys%nax*maxsys%nsx)
          ELSE
             CALL zeroing(pivelp(:,:,:,ipcurr))!,3*maxsys%nax*maxsys%nsx)
          ENDIF
       ENDIF
       IF (cntl%quenche) CALL zeroing(cm(:,:,ipx))!,ngw*n)
    ENDDO
    ! ..Global distribution of ionic velocities
    IF (pimd1%tstage.OR.pimd1%tpinm) THEN
       IF (restart1%rvel) THEN
          CALL global(pivelp,3*maxsys%nax*maxsys%nsx)
          IF (paral%parent) THEN
             IF (pimd1%tstage) THEN
                CALL stagetrans(pivelp,vstage,0)
             ELSEIF (pimd1%tpinm) THEN
                CALL pinmtrans(pivelp,vstage,0)
             ENDIF
          ENDIF
       ELSE
          CALL global(vstage,3*maxsys%nax*maxsys%nsx)
          IF (paral%parent) THEN
             IF (pimd1%tstage) THEN
                CALL stagetrans(pivelp,vstage,1)
             ELSEIF (pimd1%tpinm) THEN
                CALL pinmtrans(pivelp,vstage,1)
             ENDIF
          ENDIF
       ENDIF
    ELSE
       CALL global(pivelp,3*maxsys%nax*maxsys%nsx)
    ENDIF
    IF (cntl%quenchc) CALL zeroing(metr_com%htvel)!,9)
    IF (paral%parent) THEN
       ! ..Reset accumulators
       IF (.NOT.restart1%rac) THEN
          iteropt%nfi=0
          CALL zeroing(accus)!,nacx*maxnp)
          CALL dcopy(3*maxsys%nax*maxsys%nsx*npx,pitaup,1,taui,1)
       ENDIF
    ENDIF
    ! ..Initialize forces
    DO ipcurr=np_low,np_high
       ipx=ipcurr-np_low+1
       CALL phfac(pitaup(:,:,:,ipcurr))
       IF (corel%tinlc) CALL copot(rhoe,psi,.TRUE.)
       ! INITIALIZE WF CENTERS & SPREAD
       IF (vdwl%vdwd) THEN
          IF (.NOT.trwanncx(ipx)) THEN
             CALL localize2(pitaup(:,:,:,ipcurr),c0(:,:,ipx),c2(:,:,ipx),sc0,crge%n)
          ENDIF
       ENDIF
       CALL forcedr(c0(:,:,ipx),c2(:,:,ipx),sc0,rhoe,psi,&
            pitaup(:,:,:,ipcurr),pifion(:,:,:,ipcurr),eigv(:,ipx),&
            crge%n,1,.FALSE.,.TRUE.)
       CALL dscal(3*maxsys%nax*maxsys%nsx,1._real_8/rnp,pifion(1,1,1,ipcurr),1)
       CALL dcopy(3*maxsys%nax*maxsys%nsx,pifion(1,1,1,ipcurr),1,fionks(1,1,1,ipcurr),1)
       CALL fharm(pitaup,pifion(:,:,:,ipcurr),ipcurr,.TRUE.)
       CALL pi_ener(ipcurr,0)
       IF (ipcurr.EQ.np_low) CALL freqs(crge%n,.FALSE.)
       ! ..Check orthogonality condition for wavefunction velocities
       CALL rortv(c0(:,:,ipx),cm(1,1,ipx),c2(1,1,ipx),sc0,&
            gamy(1,ipx),crge%n)
       CALL totstr
       CALL dcopy(9,paiu(1,1),1,paiuall(1,1,ipcurr),1)
       CALL dscal(9,1._real_8/rnp,paiuall(1,1,ipcurr),1)
    ENDDO
    IF (paral%parent) THEN
       CALL global(fionks,3*maxsys%nax*maxsys%nsx)
       CALL global(paiuall,9)
       IF (pimd1%tstage) THEN
          CALL getfu(fionks,fstage,stagep,1)
          CALL taucl(fstage(:,:,:,1))
       ELSEIF (pimd1%tpinm) THEN
          CALL getfnm(fionks,fstage,stagep,1)
          CALL taucl(fstage(:,:,:,1))
       ENDIF
       ! ..Initialize thermostats
       DO ipcurr=np_low,np_high
          nosl%tmnose=pitmnose(MIN(ipcurr,2))
          cntl%tnosec=pitnosec(MIN(ipcurr,2))
          itemp=irec(irec_nop1,ipcurr)+irec(irec_nop2,ipcurr)+&
                irec(irec_nop3,ipcurr)+irec(irec_nop4,ipcurr)
          IF (cntl%tnosep .AND. itemp.EQ.0) CALL nospinit(ipcurr)
          IF (cntl%tnosee .AND. .NOT.restart1%rnoe) CALL noseinit(ipcurr)
          IF (cntl%tnosec .AND. irec(irec_noc,ipcurr).EQ.0) CALL noscinit(ipcurr)
          ! DM1
          ! ..Centroid & Ring-polmer PIMD: Reset thermostat variables in case of 
          ! ..a centroid restart from a non-centroid thermosatted run
          IF (cntl%tnosep.AND.(pimd1%tcentro.OR.pimd1%tringp).AND.ipcurr.EQ.1) THEN
             nrepl=maxnp
             IF (nosl%tultra) THEN
                ! ..Not used with PICPMD 
                CALL stopgm(procedureN,'Ultra NHCs thermostat not supported',& 
                     __LINE__,__FILE__)
             ELSEIF (loct%tloct) THEN
                CALL stopgm(procedureN,'LOCAL TEMEPERATURE NOT IMPLEMENTED',& 
                     __LINE__,__FILE__)
             ELSEIF (nosl%tmnose) THEN
                IF (pimd1%tcentro) THEN
                   CALL stopgm(procedureN,'Massive NHCs thermostat not available for ip=1',&
                        __LINE__,__FILE__)
                ELSE
                   ! ..Switch off if requested
                   IF (.NOT.tnosepc) THEN
                      DO k=1,3*maxsys%nax*maxsys%nsx
                         DO l=1,nchx
                            etapm(k,l,1)=0._real_8
                            etapmdot(k,l,1)=0._real_8
                         ENDDO
                      ENDDO
                   ENDIF
                ENDIF
             ELSE
                ! ..Switch off if requested
                IF (.NOT.tnosepc) THEN
                   DO l=1,nchx
                      etap1(l,1)=0._real_8
                      etap1dot(l,1)=0._real_8
                   ENDDO
                ENDIF
             ENDIF
          ENDIF
          ! DM2
          IF (ipcurr.EQ.1) THEN
             ! Initialize Barostat
             shock1%vol0 = prcp_com%omega0
             IF (prcpl%tzflex) THEN
                eps = LOG(parm%omega/prcp_com%omega0)/3._real_8
                veps = metr_com%htvel(3,3)
             ELSE IF (prcpl%tisot) THEN
                eps = LOG(parm%omega/prcp_com%omega0)/3._real_8
                veps = (metr_com%htvel(1,1)+metr_com%htvel(2,2)+metr_com%htvel(3,3))/3._real_8
             ELSE IF (cntl%tshock) THEN
                shock1%vol0 = prcp_com%omega0
             ENDIF
          ENDIF
       ENDDO
    ENDIF
    ! stress tensor & cell forces
    IF (pimd1%tstage.OR.pimd1%tpinm) THEN
       pipaiu=pi_stress_vir(paiuall,stagep(:,:,:,1),vstage(:,:,:,1),pitaup,fionks,.FALSE.)
    ELSE
       CALL stopgm(procedureN,'not yet implemented',&
            __LINE__,__FILE__)
    ENDIF
    ! KINETIC ENERGY OF THE CELL
    ekinh=0.5_real_8*cntr%cmass*ddot(9,metr_com%htvel,1,metr_com%htvel,1)
    IF (grandparent) THEN
       IF (pimd3%levprt.GE.5) THEN
          DO ip=1,pimd3%np_total
             IF (paral%io_parent)&
                  WRITE(6,'(A,I3)') ' *** REPLICA :',ip
             CALL wrgeof(pitaup(:,:,:,ip),pifion(:,:,:,ip))
          ENDDO
       ENDIF
       CALL getgyration(pitaup,pivelp)!vw routine requies 2 dummy args, added pivelp
       CALL getcor(pitaup)
       CALL prtgyr
       filen='ENERGIES'
       IF (paral%io_parent)&
            CALL fileopen(3,filen,fo_app+fo_verb+fo_mark,ferror)
       IF (paral%io_parent)&
            WRITE(6,'(/,A)') '  CELL PARAMETERS '
       DO i=1,3
          IF (paral%io_parent)&
               WRITE(6,'(3(1X,F12.8))') (metr_com%ht(i,j),j=1,3)
       ENDDO
       IF (paral%io_parent)&
            WRITE(6,'(/,A)') '  CELL VELOCITY '
       DO i=1,3
          IF (paral%io_parent)&
               WRITE(6,'(3(1X,F12.8))') (metr_com%htvel(i,j),j=1,3)
       ENDDO
       IF (paral%io_parent)&
            WRITE(6,'(/,A)') '  CELL FORCES '
       DO i=1,3
          IF (paral%io_parent)&
               WRITE(6,'(3(1X,F12.8))') (metr_com%htfor(i,j),j=1,3)
       ENDDO
       ! ..END INITIALIZATION
       CALL prmem('     PI_NPT_CPMD')
       time2 =m_walltime()
       tcpu = (time2 - time1)*0.001_real_8
       IF (paral%io_parent)&
            WRITE(6,'(A,F8.2,A)') ' TIME FOR INITIALIZATION',&
            tcpu,' SECONDS'
    ENDIF
    CALL mp_sync(supergroup)
    ALLOCATE(center(4,crge%n),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ==================================================================
    ! ==          THE BASIC LOOP FOR MOLECULAR DYNAMICS               ==
    ! ==                 USING VELOCITY VERLET                        ==
    ! ==================================================================
    DO infi=1,cnti%nomore
       time1=m_walltime()
       iteropt%nfi=iteropt%nfi+1
       ropt_mod%prteig=MOD(iteropt%nfi-1,cprint%iprint_step).EQ.0
       IF (.NOT.paral%parent) ropt_mod%prteig=.FALSE.
       ropt_mod%engpri=MOD(iteropt%nfi-1,cprint%iprint_step).EQ.0
       comvl%subcom=comvl%tsubcom.AND.MOD(iteropt%nfi-1,comvl%ncomv).EQ.0
       comvl%subrot=comvl%tsubrot.AND.MOD(iteropt%nfi-1,comvl%nrotv).EQ.0
       ! >>>>Replica loop
       DO ipcurr=np_low,np_high
          ipx=ipcurr-np_low+1
          ! ANNEALING
          IF (pimd1%tstage.OR.pimd1%tpinm) THEN
             CALL anneal(vstage(:,:,:,ipcurr),cm(:,:,ipx),crge%n,metr_com%htvel)
             CALL berendsen(vstage(:,:,:,ipcurr),cm(:,:,ipx),crge%n,metr_com%htvel,ekinc(ipcurr),ekinh)
          ELSE
             CALL anneal(pivelp(:,:,:,ipcurr),cm(:,:,ipx),crge%n,metr_com%htvel)
             CALL berendsen(pivelp(:,:,:,ipcurr),cm(:,:,ipx),crge%n,metr_com%htvel,ekinc(ipcurr),ekinh)
          ENDIF
          ! ..UPDATE NOSE THERMOSTATS
          nosl%tmnose=pitmnose(MIN(ipcurr,2))
          cntl%tnosec=pitnosec(MIN(ipcurr,2))
          IF (reset_gkt.AND.ipcurr==1) THEN
             scr(1:2)=gkt(1:2,1); gkt(1:2,1)=gkt1(1:2,1)
          ENDIF
          IF (pimd1%tstage.OR.pimd1%tpinm) THEN
             CALL noseup(vstage(:,:,:,ipcurr),cm(1,1,ipx),crge%n,ipcurr)
          ELSE
             CALL stopgm(procedureN,'not yet implemented',&
                  __LINE__,__FILE__)
          ENDIF
          IF (reset_gkt.AND.ipcurr==1) gkt(1:2,1)=scr(1:2)
          IF ((pimd1%tcentro.OR.pimd1%tringp).AND.ipcurr.EQ.1) THEN
             IF (pimd1%tstage.OR.pimd1%tpinm) THEN
                ! SUBTRACT CENTER OF MASS VELOCITY
                IF (paral%parent.AND.comvl%subcom) THEN
                   CALL comvel(vstage(:,:,:,1),vcmio,.TRUE.)
                ENDIF
                ! SUBTRACT ROTATION AROUND CENTER OF MASS
                IF (paral%parent.AND.comvl%subrot) THEN
                   CALL rotvel(stagep(:,:,:,1),vstage(:,:,:,1),lmio,tauio,.TRUE.)
                ENDIF
             ENDIF
          ENDIF
          ! ..UPDATE VELOCITIES
          IF (paral%parent) THEN
             IF (pimd1%tstage.OR.pimd1%tpinm) THEN
                CALL s_velupi(vstage(:,:,:,ipcurr),fstage(:,:,:,ipcurr),1,ipcurr)
                IF (ipcurr.EQ.1) CALL velupc(1)
             ELSE
                CALL stopgm(procedureN,'not yet implemented',&
                     __LINE__,__FILE__)
             ENDIF
          ENDIF
          CALL velupa(c0(:,:,ipx),cm(1,1,ipx),c2(1,1,ipx),crge%n,1)
          ! ..UPDATE POSITIONS
          ener_com%ecnstr = 0.0_real_8
          ener_com%erestr = 0.0_real_8
          IF (paral%parent) THEN
             IF (pimd1%tstage.OR.pimd1%tpinm) THEN
                CALL dcopy(3*maxsys%nax*maxsys%nsx,stagep(1,1,1,ipcurr),1,tau0,1)
                IF (ipcurr.EQ.1) THEN
                   IF (prcpl%tisot) THEN
                      CALL posupih_iso(tau0,stagep(:,:,:,1),vstage(:,:,:,1))
                   ELSE
                      CALL posupih(tau0,stagep(:,:,:,1),vstage(:,:,:,1))
                   ENDIF
                   IF (cotc0%mcnstr.NE.0)&
                      CALL cpmdshake(tau0,stagep(:,:,:,1),vstage(:,:,:,1))
                ELSE
                   CALL posupi(tau0,stagep(:,:,:,ipcurr),vstage(:,:,:,ipcurr))
                ENDIF
             ELSE
                CALL stopgm(procedureN,'not yet implemented',&
                     __LINE__,__FILE__)
             ENDIF
          ENDIF
       ENDDO
       IF (paral%io_parent) THEN
          CALL mp_bcast(metr_com%ht,SIZE(metr_com%ht),supersource,parentgroup)
       ENDIF
       CALL mp_bcast(metr_com%ht,SIZE(metr_com%ht),parai%io_source,parai%cp_grp)
       DO i=1,3
          parm%a1(i) = metr_com%ht(1,i)
          parm%a2(i) = metr_com%ht(2,i)
          parm%a3(i) = metr_com%ht(3,i)
       ENDDO
       CALL ihmat(metr_com%ht,metr_com%htm1,parm%omega)
       CALL newcell
       ! <<<<End of replica loop
       IF (paral%parent) THEN
          IF (pimd1%tpinm) THEN
             CALL global(stagep,3*maxsys%nax*maxsys%nsx)
             CALL pinmtrans(pitaup,stagep,1)
          ELSEIF (pimd1%tstage) THEN
             CALL global(stagep,3*maxsys%nax*maxsys%nsx)
             CALL stagetrans(pitaup,stagep,1)
          ENDIF
       ENDIF
       ! ..Synchronization point! Distribute new replica coordinates to all PE
       CALL mp_sync(supergroup)
       CALL global(pitaup,3*maxsys%nax*maxsys%nsx)
       CALL mp_bcast(pitaup,SIZE(pitaup),parai%io_source,parai%cp_grp)
       ! >>>>Replica loop
       DO ipcurr=np_low,np_high
          ipx=ipcurr-np_low+1
          CALL phfac(pitaup(:,:,:,ipcurr))
          IF (corel%tinlc) THEN
             CALL copot(rhoe,psi,.TRUE.)
          ENDIF
          CALL posupa(c0(:,:,ipx),cm(1,1,ipx),c2(1,1,ipx),&
               gamx(1,ipx),crge%n)
          ! DM     *                GAMX,N,SCR,LSCR)
          ! Wannier stuff for vdW-WC
          IF (vdwl%vdwd) THEN
             IF (twannupx(ipx)) THEN
                CALL ddipo(pitaup(:,:,:,ipcurr),c0(:,:,ipx),cm(:,:,ipx),c2(:,:,ipx),&
                     sc0,crge%n,center)
                CALL wannier_print(iteropt%nfi,c0(:,:,ipx),pitaup(:,:,:,ipcurr),crge%n,&
                     psi(:,ipx),center)
             ENDIF
          ENDIF
          ! ..CALCULATE THE FORCES
          CALL forcedr(c0(:,:,ipx),c2(:,:,ipx),sc0,rhoe,psi,&
               pitaup(:,:,:,ipcurr),pifion(:,:,:,ipcurr),eigv(:,ipx),&
               crge%n,1,.FALSE.,.TRUE.)
          CALL dscal(3*maxsys%nax*maxsys%nsx,1._real_8/rnp,pifion(1,1,1,ipcurr),1)
          CALL dcopy(3*maxsys%nax*maxsys%nsx,pifion(1,1,1,ipcurr),1,&
               fionks(1,1,1,ipcurr),1)
          CALL fharm(pitaup,pifion(:,:,:,ipcurr),ipcurr,.TRUE.)
          CALL pi_ener(ipcurr,0)
          CALL totstr
          CALL dcopy(9,paiu(1,1),1,paiuall(1,1,ipcurr),1)
          CALL dscal(9,1._real_8/rnp,paiuall(1,1,ipcurr),1)
       ENDDO
       ! <<<<End of replica loop
       IF (paral%parent) THEN
          CALL global(fionks,3*maxsys%nax*maxsys%nsx)
          CALL global(paiuall,9)
          IF (pimd1%tpinm) THEN
             CALL getfnm(fionks,fstage,stagep,1)
             CALL taucl(fstage(:,:,:,1))
          ELSEIF (pimd1%tstage) THEN
             CALL getfu(fionks,fstage,stagep,1)
             CALL taucl(fstage(:,:,:,1))
          ENDIF
       ENDIF
       ! ..Update cell forces
       IF (pimd1%tstage.OR.pimd1%tpinm) THEN
          pipaiu=pi_stress_vir(paiuall,stagep(:,:,:,1),vstage(:,:,:,1),pitaup,fionks,.FALSE.)
       ELSE
          CALL stopgm(procedureN,'not yet implemented',&
               __LINE__,__FILE__)
       ENDIF
       ! >>>>Replica loop
       DO ipcurr=np_low,np_high
          ipx=ipcurr-np_low+1
          ! ..FINAL UPDATE FOR VELOCITIES
          ener_com%ecnstr = 0.0_real_8
          ener_com%erestr = 0.0_real_8
          IF (paral%parent) THEN
             IF (pimd1%tpinm.OR.pimd1%tstage) THEN
                ! IONS & CELL
                IF (ipcurr.EQ.1) CALL velupc(1)
                CALL s_velupi(vstage(:,:,:,ipcurr),fstage(:,:,:,ipcurr),1,ipcurr)
                IF (ipcurr.EQ.1) CALL rattle(stagep(:,:,:,1),vstage(:,:,:,1))
             ELSE
                CALL stopgm(procedureN,'not yet implemented',&
                     __LINE__,__FILE__)
             ENDIF
          ENDIF
          CALL velupa(c0(:,:,ipx),cm(1,1,ipx),c2(1,1,ipx),crge%n,1)
          CALL rortv(c0(:,:,ipx),cm(1,1,ipx),c2(1,1,ipx),&
               sc0,gamy(1,ipx),crge%n)
          ! ..IONIC TEMPERATURE CONTROL
          IF (paral%parent.AND.(cntl%tcp.OR.cntl%trescale)) THEN
             ! DMra          IF(PARENT.AND.(cntl%tcp.OR.ANNE)) THEN
             ! DM          IF(PARENT) THEN
             IF (pimd1%tpinm.OR.pimd1%tstage) THEN
                CALL s_ekinpp(ekinp,vstage(:,:,:,ipcurr),ipcurr)
                ! DM1
                ! ..GLIB is calculated in QMDOF for the first (staging) bead 
                ! ..that is possible constraint and is 3*NAT for all other beads
                ! DM2
                glib_s=3._real_8*ions1%nat
                IF (ipcurr.EQ.1) glib_s=glib
                tempp(ipcurr)=ekinp*factem*2._real_8/glib_s
                ! DM1ra
                ! DM              CALL RSCVP(TEMP1,TEMP2,TEMPP(IPCURR),VSTAGE(1,IPCURR))
                IF (.NOT.cntl%trescale) THEN
                   CALL rscvp(temp1,temp2,tempp(ipcurr),vstage(:,:,:,ipcurr))
                ELSE
                   ! ..Re-adjust ionic velocities  only after restart to desired temperature
                   IF (ipcurr.EQ.1) THEN
                      cntl%tcp=.TRUE.
                      CALL rscvp(0._real_8,0._real_8,tempp(ipcurr),vstage(:,:,:,ipcurr))
                      cntl%tcp=.FALSE.
                      cntl%trescale=.FALSE.
                      IF (paral%io_parent)&
                           WRITE(6,*)
                      IF (paral%io_parent)&
                           WRITE(6,'(A,F12.2,A)')&
                           ' RESTARTED VELOCITIES FOR REPLICA IP=1 RESCALED TO ',&
                           cntr%tempw,' KELVIN'
                   ENDIF
                ENDIF
                ! DM2ra
             ELSE
                CALL ekinpp(ekinp,pivelp(:,:,:,ipcurr))
                tempp(ipcurr)=ekinp*factem*2._real_8/glib
                ! DM1ra              CALL RSCVP(TEMP1,TEMP2,TEMPP(IPCURR),PIVELP(1,IPCURR))
                IF (.NOT.cntl%trescale) THEN
                   CALL rscvp(temp1,temp2,tempp(ipcurr),pivelp(:,:,:,ipcurr))
                ELSE
                   ! ..Re-adjust ionic velocities  only after restart to desired temperature
                   IF (ipcurr.EQ.1) THEN
                      cntl%tcp=.TRUE.
                      CALL rscvp(0._real_8,0._real_8,tempp(ipcurr),pivelp(:,:,:,ipcurr))
                      cntl%tcp=.FALSE.
                      cntl%trescale=.FALSE.
                      IF (paral%io_parent)&
                           WRITE(6,*)
                      IF (paral%io_parent)&
                           WRITE(6,'(A,F12.2,A)')&
                           ' RESTARTED VELOCITIES FOR REPLICA IP=1 RESCALED TO ',&
                           cntr%tempw,' KELVIN'
                   ENDIF
                ENDIF
                ! DM2ra
             ENDIF
          ENDIF
          IF ((pimd1%tcentro.OR.pimd1%tringp).AND.ipcurr.EQ.1) THEN
             IF (pimd1%tstage.OR.pimd1%tpinm) THEN
                ! SUBTRACT ROTATION AROUND CENTER OF MASS
                IF (paral%parent.AND.comvl%subrot) THEN
                   CALL rotvel(stagep(:,:,:,1),vstage(:,:,:,1),lmio,tauio,.FALSE.)
                ENDIF
                ! SUBTRACT CENTER OF MASS VELOCITY
                IF (paral%parent.AND.comvl%subcom) THEN
                   CALL comvel(vstage(:,:,:,1),vcmio,.FALSE.)
                ENDIF
             ENDIF
          ENDIF
          ! RESCALE KINETIC ENERGY OF THE CELL
          IF (paral%parent.AND.ipcurr.EQ.1) THEN
             ekinh=0.5_real_8*cntr%cmass*ddot(9,metr_com%htvel,1,metr_com%htvel,1)
             CALL rscvc(ekinh,metr_com%htvel,ekinh1,ekinh2)
          ENDIF
          ! ..UPDATE NOSE THERMOSTATS
          IF (cntl%tnosee.OR.cntl%tc) CALL rekine(cm(1,1,ipx),crge%n,ekinc(ipcurr))
          nosl%tmnose=pitmnose(MIN(ipcurr,2))
          cntl%tnosec=pitnosec(MIN(ipcurr,2))
          IF (reset_gkt.AND.ipcurr==1) THEN
             scr(1:2)=gkt(1:2,1); gkt(1:2,1)=gkt1(1:2,1)
          ENDIF
          IF (pimd1%tstage.OR.pimd1%tpinm) THEN
             CALL noseup(vstage(:,:,:,ipcurr),cm(1,1,ipx),crge%n,ipcurr)
          ELSE
             CALL stopgm(procedureN,'not yet implemented',&
                  __LINE__,__FILE__)
          ENDIF
          IF (reset_gkt.AND.ipcurr==1) gkt(1:2,1)=scr(1:2)
          ekinh=0.5_real_8*cntr%cmass*ddot(9,metr_com%htvel,1,metr_com%htvel,1)
          IF (pimd1%tpinm.OR.pimd1%tstage) THEN
             CALL berendsen(vstage(:,:,:,ipcurr),cm(:,:,ipx),crge%n,metr_com%htvel,ekinc(ipcurr),ekinh)
             ! ANNEALING
             CALL anneal(vstage(:,:,:,ipcurr),cm(:,:,ipx),crge%n,metr_com%htvel)
          ELSE
             CALL berendsen(pivelp(:,:,:,ipcurr),cm(:,:,ipx),crge%n,metr_com%htvel,ekinc(ipcurr),ekinh)
             ! ANNEALING
             CALL anneal(pivelp(:,:,:,ipcurr),cm(:,:,ipx),crge%n,metr_com%htvel)
          ENDIF
          IF (paral%parent) THEN
             IF (pimd1%tpinm.OR.pimd1%tstage) THEN
                CALL s_ekinpp(ekinp,vstage(:,:,:,ipcurr),ipcurr)
                glib_s=3._real_8*ions1%nat
                IF (ipcurr.EQ.1) glib_s=glib
                tempp(ipcurr)=ekinp*factem*2._real_8/glib_s
             ELSE
                CALL ekinpp(ekinp,pivelp(:,:,:,ipcurr))
                tempp(ipcurr)=ekinp*factem*2._real_8/glib
             ENDIF
          ENDIF
          ! ..RESCALE ELECTRONIC VELOCITIES
          IF (cntl%tc) CALL rscve(ekin1,ekin2,ekinc(ipcurr),&
               cntr%ekinw,cm(1,1,ipx),crge%n,ncpw%ngw)
          ! ..MEAN SQUARE DISPLACEMENT OF DIFFERENT IONIC SPECIES
          IF (paral%parent) CALL dispp(pitaup(:,:,:,ipcurr),taui(:,:,:,ipcurr),&
               disa(ipcurr))
          ! ..KINETIC ENERGY OF THE ELECTRONS
          CALL rekine(cm(1,1,ipx),crge%n,ekinc(ipcurr))
          IF (paral%parent) THEN
             ! ..ENERGY OF THE NOSE THERMOSTATS
             nosl%tmnose=pitmnose(MIN(ipcurr,2))
             cntl%tnosec=pitnosec(MIN(ipcurr,2))
             IF (reset_gkt.AND.ipcurr==1) THEN
                scr(1:2)=gkt(1:2,1); gkt(1:2,1)=gkt1(1:2,1)
             ENDIF
             IF (pimd1%tstage.OR.pimd1%tpinm) THEN
                CALL noseng(iteropt%nfi,vstage(:,:,:,ipcurr),enose(ipcurr),&
                     enosp(ipcurr),enosc(ipcurr),ipcurr)
             ELSE
                CALL stopgm(procedureN,'not yet implemented',&
                     __LINE__,__FILE__)
             ENDIF
             IF (reset_gkt.AND.ipcurr==1) gkt(1:2,1)=scr(1:2)
             ! DM1
             ! ..EFFECTIVELY CONSERVED AND HAMILTONIAN ENERGY
             ! ..NOTE: ECNSTR IS ALWAYS ZERO
             ! DM2
             econs(ipcurr)=ekinp+etotv(ipcurr)/rnp+&
                  enose(ipcurr)/rnp+enosp(ipcurr)+enosc(ipcurr)+&
                  eharv(ipcurr)+ener_com%ecnstr+ener_com%erestr
             eham(ipcurr)=econs(ipcurr)+ekinc(ipcurr)/rnp
          ENDIF
       ENDDO
       ! <<<<End of Replica Loop
       CALL mp_sync(supergroup)
       IF (paral%parent) THEN
          IF (pimd1%tpinm) THEN
             CALL global(vstage,3*maxsys%nax*maxsys%nsx)
             CALL pinmtrans(pivelp,vstage,1)
          ELSEIF (pimd1%tstage) THEN
             CALL global(vstage,3*maxsys%nax*maxsys%nsx)
             CALL stagetrans(pivelp,vstage,1)
          ENDIF
          DO ipcurr=np_low,np_high
             CALL geofile(pitaup(:,:,:,ipcurr),pivelp(:,:,:,ipcurr),'WRITE')
          ENDDO
       ENDIF
       CALL disten
       CALL global(tempp,1)
       CALL global(eham,1)
       CALL global(econs,1)
       CALL global(ekinc,1)
       CALL global(enose,1)
       CALL global(enosp,1)
       CALL global(enosc,1)
       CALL global(disa,1)
       ! DM C..PRINTOUT the evolution of the accumulators every time step
       IF (grandparent) THEN
          ! ..Radii of gyration
          CALL getgyration(pitaup,pivelp)!vw routine requies 2 dummy args, added pivelp
          ! ..Trotter correlation functions for positions and KS-forces
          CALL getcor(pitaup)
          ! ..extended output
          IF (ropt_mod%engpri) THEN
             IF (pimd3%levprt.GE.5) THEN
                DO ip=1,pimd3%np_total
                   IF (paral%io_parent)&
                        WRITE(6,'(1X,63("*"))')
                   IF (paral%io_parent)&
                        WRITE(6,'(A,I3,1X,46("*"))') ' *** REPLICA :',ip
                   IF (paral%io_parent)&
                        WRITE(6,'(1X,63("*"))')
                   CALL pi_ener(ip,1)
                   CALL wrener
                   IF (cprint%iprint(iprint_coor ).EQ.1) CALL wrgeo(pitaup(:,:,:,ip))
                   IF (cprint%iprint(iprint_force).EQ.1)&
                        CALL wrgeof(pitaup(:,:,:,ip),pifion(:,:,:,ip))
                ENDDO
             ENDIF
             IF (cprint%iprint(iprint_coor).EQ.1) CALL cnstpr
             CALL prtgyr
             IF (paral%io_parent)&
                  WRITE(6,'(/,A,A)')&
                  '      NFI   IP   EKINC  TEMPP     EHARM         EKS',&
                  '  ENOSE(C)  ENOSE(E)  ENOSE(P)    ECLASSIC        EHAM     DIS'
             DO ip=1,pimd3%np_total
                IF (paral%io_parent)&
                     WRITE(6,&
                     '(I9,I5,F8.5,F7.1,F10.6,F12.5,3F10.5,2F12.5,F8.2)')&
                     iteropt%nfi,ip,ekinc(ip),tempp(ip),eharv(ip),etotv(ip),&
                     enosc(ip),enose(ip),enosp(ip),econs(ip),eham(ip),disa(ip)
             ENDDO
          ENDIF
          ! ..Quantum kinetic energy of ions
          ! DM1
          ! ..Primitive (Barker) estimator
          ! DM C..Barker estimator
          ! DM2
          qpkinp=0.0_real_8
          DO ip=1,pimd3%np_total
             qpkinp=qpkinp+eharv(ip)
          ENDDO
          qpkinp=(glib+3._real_8*ions1%nat*(rnp-1))*cntr%tempw/(2._real_8*315777._real_8)-qpkinp
          ! ..Virial (Berne) estimator from KS force
          CALL evirial(pitaup,qvkinp)
          ! KINETIC ENERGY OF THE CELL
          ekinh=0.5_real_8*cntr%cmass*ddot(9,metr_com%htvel,1,metr_com%htvel,1)
          ! TEMPERATURE OF THE CELL
          temph = 2._real_8*ekinh*factem/REAL(ncdof,kind=real_8)
          ! ..System energies
          ekina=0.0_real_8
          etota=0.0_real_8
          econsa=0.0_real_8
          ehama=0.0_real_8
          tempa=0.0_real_8
          DO ip=1,pimd3%np_total
             ekina=ekina+ekinc(ip)/rnp
             etota=etota+etotv(ip)/rnp
             econsa=econsa+econs(ip)
             ehama=ehama+eham(ip)
             tempa=tempa+tempp(ip)/rnp
          ENDDO
          equant=etota+qpkinp
          econsa=econsa+prcp_com%druck*parm%omega
          ehama=ehama+prcp_com%druck*parm%omega+ekinh
          time2=m_walltime()
          tcpu=(time2-time1)*0.001_real_8
          IF ((ropt_mod%engpri).AND.paral%io_parent)&
               WRITE(6,'(/,A,A)')&
               '      NFI  EKINC/P    EKINH  TEMP  EKINP(P)  EKINP(V)       EKS/P',&
               '       EQUANT     ECLASSIC         EHAM     TCPU'
          IF (paral%io_parent)&
               WRITE(6,'(I9,2F9.5,F7.1,2F9.5,4F13.5,F9.2)')&
               iteropt%nfi,ekina,ekinh,tempa,qpkinp,qvkinp,etota,equant,econsa,ehama,tcpu
          IF (paral%io_parent)&
               WRITE(3,'(I9,2F15.10,F7.1,6F20.10,F9.2)')&
               iteropt%nfi,ekina,ekinh,tempa,qpkinp,qvkinp,etota,equant,econsa,ehama,tcpu
          ! ..Store ionic coordinates and velocities for statistics
          ropt_mod%movie=rout1%mout .AND. MOD(iteropt%nfi-1,cnti%imovie).EQ.0
          ropt_mod%rprint=rout1%rout .AND. MOD(iteropt%nfi-1,cnti%ntraj).EQ.0
          ropt_mod%txyz=rout1%xtout .AND. MOD(iteropt%nfi-1,cnti%ntraj).EQ.0
          DO ip=1,pimd3%np_total
             ! wr1
             ! wr       note: STAGEP and VSTAGE are already globally updated
             ! wr            CALL PRINTP(TAUR,PITAUP(1,IP),PIVELP(1,IP),
             ! wr     *                  STAGEP(1,IP),VSTAGE(1,IP))
             ! wr2
             CALL printp(taur,pitaup(:,:,:,ip),pivelp(:,:,:,ip))
          ENDDO
       ENDIF
       IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)
       CALL testex_mw(soft_com%exsoft)
       IF (infi.EQ.cnti%nomore) soft_com%exsoft=.TRUE.
       ! new
       nosl%tmnose=pitmnose(2)  ! recover original setup of tmnose
       IF (pimd1%tstage.OR.pimd1%tpinm) CALL wr_temps(iteropt%nfi,ekinc,tempp)
       ! new
       ! ..Recompute velocity-dependent part of pipaiu
       IF (pimd1%tstage.OR.pimd1%tpinm) THEN
          pipaiu=pi_stress_vir(paiuall,stagep(:,:,:,1),vstage(:,:,:,1),pitaup,fionks,.TRUE.)
       ELSE
          CALL stopgm(procedureN,'not yet implemented',&
               __LINE__,__FILE__)
       ENDIF
       IF (grandparent) CALL wr_pi_stress(pipaiu)
       IF (paral%parent) THEN
          DO ipcurr=np_low,np_high
             ! ..Update accumulators
             i=ipcurr
             CALL pacce(ekinc(i),etotv(i),econs(i),eham(i),tempp(i),&
                  eharv(i),enose(i),enosp(i),disa(i),&
                  tcpu,qpkinp,qvkinp,ekina,etota,econsa,&
                  ehama,tempa,equant,temph,parm%omega,accus(1,i),iteropt%nfi,1)
          ENDDO
       ENDIF
       CALL global(accus,nacx)
       DO ipcurr=np_low,np_high
          ! ..Construct filenames
          cflbod='RESTART_'
          IF (paral%io_parent)&
               WRITE(cipnum,'(I4)') ipcurr
          CALL xstring(cflbod,n1,n2)
          CALL xstring(cipnum,i1,i2)
          filbod=cflbod(n1:n2)//cipnum(i1:i2)//'.1'
          ! ..Write restart file
          ipx=ipcurr-np_low+1
          CALL dcopy(nacc,accus(1,ipcurr),1,acc,1)
          IF (teststore(infi).OR.soft_com%exsoft) THEN
             nosl%tmnose=pitmnose(MIN(ipcurr,2))
             cntl%tnosec=pitnosec(MIN(ipcurr,2))
             CALL write_irec(irec(:,ipcurr))
             CALL zhwwf(2,irec(:,ipcurr),c0(:,:,ipx),cm(:,:,ipx),crge%n,eigv(:,ipx),&
                  pitaup(:,:,:,ipcurr),pivelp(:,:,:,ipcurr),taui(:,:,:,ipcurr),iteropt%nfi)
          ENDIF
       ENDDO
       ! ..STOP THE RUN IF THE USER HAS SET THE SIGNAL 30
       IF (soft_com%exsoft) GOTO 100
       ! ==================================================================
       ! ==     END OF MAIN LOOP                                         ==
       ! ==================================================================
    ENDDO
100 CONTINUE
    ! ===================================================================
    ! AK: synchronisation point. risk of losing restarts on slow networks.
    CALL mp_sync(supergroup)
    ! ===================================================================
    ! ..PRINT ACCUMULATOR
    IF (rout1%rhoout) CALL rhopri(c0,tau0,rhoe,psi(:,1),crge%n,1)
    IF (grandparent) THEN
       ! Dummy arguments
       d1=0._real_8
       d2=0._real_8
       CALL pacce(d1,d1,d1,d1,d1,d1,d1,d1,d1,&
            d2,d2,d2,d2,d2,d2,d2,d2,d2,d2,d2,accus,iteropt%nfi,2)
    ENDIF
    IF (paral%parent) CALL prmem('     PI_MD')
    DEALLOCATE(center,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(eigv,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (pslo_com%tivan.OR.cntl%nonort) THEN
       DEALLOCATE(eigm,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    DEALLOCATE(fionks,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (pimd1%tstage.OR.pimd1%tpinm) THEN
       DEALLOCATE(stagep,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(vstage,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(fstage,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    IF (comvl%tsubrot) THEN
       DEALLOCATE(tauio,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    DEALLOCATE(rhoe,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(psi,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(scr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(paiuall,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(irec,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (corel%tinlc) THEN
       DEALLOCATE(vnlt,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(vnlcc,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    IF ((paral%parent).AND.paral%io_parent)&
         CALL fileclose(3)
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    RETURN
  END SUBROUTINE pi_npt_cpmd
  ! ==================================================================
  SUBROUTINE rscvc(ekinh,htvel,ekinh1,ekinh2)
    REAL(real_8)                             :: ekinh, htvel(3,3), ekinh1, &
                                                ekinh2

    INTEGER                                  :: k
    REAL(real_8)                             :: alfap

! Variables
! ==--------------------------------------------------------------==
! ==  Dynamical rescaling factor (ekinhr/ekinh), of cell kinetic  ==
! ==                                                   energy     ==
! ==--------------------------------------------------------------==

    IF (cntl%tcc) THEN
       IF (ekinh.GT.ekinh1.OR.ekinh.LT.ekinh2.AND.ekinh.NE.0._real_8) THEN
          alfap=SQRT(cntr%ekinhr/ekinh)
          DO k=1,3
             htvel(1,k)=alfap*htvel(1,k)
             htvel(2,k)=alfap*htvel(2,k)
             htvel(3,k)=alfap*htvel(3,k)
          ENDDO
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rscvc
  ! ==================================================================

END MODULE pi_npt_cpmd_utils
