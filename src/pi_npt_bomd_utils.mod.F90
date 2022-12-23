MODULE pi_npt_bomd_utils
  USE andp,                            ONLY: rin0,&
                                             rmix,&
                                             rout0
  USE anneal_utils,                    ONLY: anneal,&
                                             berendsen
  USE atwf,                            ONLY: tmovr
  USE calc_alm_utils,                  ONLY: calc_alm
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
  USE detdof_utils,                    ONLY: detdof,&
                                             qmdof
  USE dispp_utils,                     ONLY: dispp
  USE dynit_utils,                     ONLY: dynit
  USE ekinpp_utils,                    ONLY: ekinpp,&
                                             s_ekinpp
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: ener_com
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
  USE fint,                            ONLY: fint1
  USE forcep_utils,                    ONLY: rhoe_psi_size
  USE forces_diag_utils,               ONLY: forces_diag
  USE freqs_utils,                     ONLY: freqs
  USE geofile_utils,                   ONLY: geofile
  USE getcor_utils,                    ONLY: getcor
  USE getfnm_utils,                    ONLY: getfnm
  USE getfu_utils,                     ONLY: getfu
  USE getgyr_utils,                    ONLY: getgyration
  USE global_utils,                    ONLY: global
  USE ions,                            ONLY: ions1
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE localize_utils,                  ONLY: localize2
  USE machine,                         ONLY: m_walltime
  USE metr,                            ONLY: eps,&
                                             metr_com,&
                                             veps
  USE mm_extrap,                       ONLY: cold,&
                                             nnow,&
                                             numcold
  USE moverho_utils,                   ONLY: moverho
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sync
  USE newcell_utils,                   ONLY: newcell
  USE nlcc,                            ONLY: corel,&
                                             vnlcc,&
                                             vnlt
  USE noscinit_utils,                  ONLY: noscinit
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
  USE noseng_utils,                    ONLY: noseng
  USE nosepa_utils,                    ONLY: nosepa
  USE noseup_utils,                    ONLY: noseup
  USE nospinit_utils,                  ONLY: nospinit
  USE ortho_utils,                     ONLY: ortho
  USE parac,                           ONLY: parai,&
                                             paral
  USE phfac_utils,                     ONLY: phfac
  USE pi_md_utils,                     ONLY: disten,&
                                             pi_ener
  USE pimd,                            ONLY: &
       eharv, etotv, fionks, grandparent, ipcurr, maxnp, np_high, np_local, &
       np_low, pimd1, pimd3, supergroup, parentgroup, supersource
  USE pinmtrans_utils,                 ONLY: pinmtrans
  USE pi_stress_utils,                 ONLY: pi_stress_vir,&
                                             wr_pi_stress
  USE pi_wf_utils,                     ONLY: initrep
  USE poin,                            ONLY: rhoo
  USE posupi_utils,                    ONLY: posupi,&
                                             posupih,&
                                             posupih_iso
  USE prbomd_utils,                    ONLY: give_scr_prbomd
  USE prcp,                            ONLY: prcp_com,&
                                             prcpl
  USE printave_utils,                  ONLY: pacce
  USE printp_utils,                    ONLY: printp
  USE prmem_utils,                     ONLY: prmem
  USE prtgyr_utils,                    ONLY: prtgyr
  USE pslo,                            ONLY: pslo_com
  USE puttau_utils,                    ONLY: taucl
  USE ranp_utils,                      ONLY: ranp
  USE rattle_utils,                    ONLY: rattle
  USE readsr_utils,                    ONLY: xstring
  USE reshaper,                        ONLY: reshape_inplace
  USE rhoofr_c_utils,                  ONLY: rhoofr_c
  USE rhoofr_utils,                    ONLY: rhoofr
  USE rinitwf_driver,                  ONLY: rinitwf
  USE rinvel_utils,                    ONLY: rinvel,&
                                             rvscal,&
                                             s_rinvel
  USE rnlsm_utils,                     ONLY: rnlsm
  USE ropt,                            ONLY: infi,&
                                             iteropt,&
                                             ropt_mod
  USE rotvel_utils,                    ONLY: rotvel
  USE rrane_utils,                     ONLY: rrane
  USE rscvp_utils,                     ONLY: rscvp
  USE rv30_utils,                      ONLY: zhrwf
  USE setirec_utils,                   ONLY: read_irec,&
                                             write_irec
  USE setsc_utils,                     ONLY: ihmat
  USE shake_utils,                     ONLY: cpmdshake
  USE shock,                           ONLY: shock1
  USE soft,                            ONLY: soft_com
  USE spin,                            ONLY: clsd
  USE stagetrans_utils,                ONLY: stagetrans
  USE store_types,                     ONLY: &
       cprint, iprint_coor, iprint_force, irec_nop1, irec_nop2, irec_nop3, &
       irec_nop4, irec_rho, irec_wf, restart1, rout1, irec_noc
  USE strs,                            ONLY: paiu
  USE system,                          ONLY: &
       acc, cnti, cntl, cntr, fpar, maxsys, nacc, nacx, nkpt, restf, parm, ncpw
  USE testex_utils,                    ONLY: testex,&
                                             testex_mw
  USE teststore_utils,                 ONLY: teststore
  USE totstr_utils,                    ONLY: totstr
  USE vdw_utils,                       ONLY: rwannier
  USE vdwcmod,                         ONLY: &
       nfragx, nwfcx, rwann, swann, tauref, trwanncx, twannupx, vdwl, vdwwfl
  USE velupi_utils,                    ONLY: s_velupi,&
                                             velupc
  USE wr_temps_utils,                  ONLY: wr_temps
  USE wrener_utils,                    ONLY: wrener
  USE wrgeo_utils,                     ONLY: wrgeo,&
                                             wrgeof
  USE wv30_utils,                      ONLY: zhwwf
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: pi_npt_bomd

CONTAINS

  ! ==================================================================
  SUBROUTINE pi_npt_bomd(c0,cm,c2,sc0,vpp,gamx,gamy)
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8)                          :: c0(:,:,:), cm(*), c2(*), &
                                                sc0(*)
    REAL(real_8)                             :: vpp(*), &
                                                gamx(crge%n*crge%n,*), &
                                                gamy(crge%n*crge%n,*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'pi_npt_bomd'

    CHARACTER(len=100)                       :: filen
    CHARACTER(len=12)                        :: cflbod, cipnum
    CHARACTER(len=30)                        :: tag
    COMPLEX(real_8), ALLOCATABLE             :: coldall(:,:,:,:,:), psi(:,:)
    INTEGER :: i, i1, i2, ierr, ifcalc, ik, il_psi_1d, il_psi_2d, il_rhoe_1d, &
      il_rhoe_2d, ip, ipx, itemp, k, l, lenext, lscr, n1, n2, nnx, &
      npx, nrepl, nwfc, nx, j
    INTEGER, ALLOCATABLE                     :: nnowall(:), numcoldall(:), &
                                                irec(:,:)
    LOGICAL                                  :: ferror, reset_gkt, &
                                                pitmnose(2), pitnosec(2)
    REAL(real_8) :: accus(nacx,maxnp), d1, d2, disa(maxnp), dummies(1), &
      dummy, econs(maxnp), econsa, eham(maxnp), ehama, ekin1, ekin2, ekina, &
      ekinc(maxnp), ekincp, ekinh1, ekinh2, ekinh, ekinp, enosc(maxnp), &
      enose(maxnp), enosp(maxnp), equant, etota, glib_s, lmio(3), qpkinp, &
      qvkinp, rmem, rnp, tcpu, temp1, temp2, temph, tempa, tempp(maxnp), &
      time1, time2, vcmio(4), pipaiu(3,3)
    REAL(real_8), ALLOCATABLE :: eigv(:,:), fall(:,:), fstage(:,:,:,:), &
      rall(:,:), rhoe(:,:), rinp(:), rm1(:), scr(:), stagep(:,:,:,:), &
      taui(:,:,:,:), tauio(:,:), taur(:,:,:), vstage(:,:,:,:), &
      paiuall(:,:,:)
    REAL(real_8), POINTER                    :: pifion(:,:,:,:), &
                                                pitaup(:,:,:,:), &
                                                pivelp(:,:,:,:)
    REAL(real_8), EXTERNAL                   :: ddot

    accus=0.0_real_8
    disa=0.0_real_8
    econs=0.0_real_8
    eham=0.0_real_8
    ekinc=0.0_real_8
    enose=0.0_real_8
    enosp=0.0_real_8
    enosc=0.0_real_8
    tempp=0.0_real_8
    restf%nfnow=1
    time1 =m_walltime()
    npx=pimd3%np_total
    rnp=REAL(npx,kind=real_8)
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
    ALLOCATE(eigv(crge%n*nkpt%nkpts,np_local),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! Memory for occupation numbers
    ALLOCATE(fall(crge%n*nkpt%nkpts,np_local),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! Memory for densities
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
    IF (cntl%tdiag) THEN
       ALLOCATE(rall(fpar%nnr1*clsd%nlsd,nnx*np_local/(clsd%nlsd*fpar%nnr1)),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    nacc = 30
    iteropt%nfi  = 0
    CALL zeroing(accus)!,nacx*maxnp)
    CALL zeroing(ekinc)!,maxnp)
    ropt_mod%modens=.FALSE.
    ropt_mod%engpri=.FALSE.
    ropt_mod%calste=.TRUE.
    ALLOCATE(paiuall(3,3,npx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
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
    ! Extrapolation
    IF (cntl%textrap) THEN
       lenext=2*nkpt%ngwk*crge%n*nkpt%nkpts*cnti%mextra
       rmem = 16._real_8*lenext*1.e-6_real_8
       ALLOCATE(cold(nkpt%ngwk,crge%n,nkpt%nkpnt,lenext/(crge%n*nkpt%ngwk*nkpt%nkpnt)),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       call zeroing(cold)
       IF (np_local>1) THEN
          ALLOCATE(coldall(nkpt%ngwk,crge%n,nkpt%nkpnt,cnti%mextra,np_local),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          call zeroing(coldall)
          ALLOCATE(nnowall(np_local),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          call zeroing(nnowall)
          ALLOCATE(numcoldall(np_local),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          call zeroing(numcoldall)
          rmem=rmem+rmem*REAL(np_local,kind=real_8)
       ENDIF
       IF (paral%io_parent)&
            WRITE(6,'(A,T51,F8.3,A)') ' PI_NPT_BOMD| '&
            // 'EXTRAPOL WAVEFUNCTION HISTORY TAKES ',rmem,' MBYTES'
    ENDIF
    ! ==--------------------------------------------------------------==
    ! SCR ALLOCATION AND PARTITION (SCRATCH ARRAY).
    CALL rhoe_psi_size(il_rhoe_1d=il_rhoe_1d, il_rhoe_2d=il_rhoe_2d, &
         il_psi_1d=il_psi_1d,il_psi_2d=il_psi_2d)
    ALLOCATE(rhoe(il_rhoe_1d,il_rhoe_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(psi(il_psi_1d,il_psi_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL give_scr_prbomd(lscr,tag)
    ALLOCATE(scr(lscr),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ..TIME STEP FUNCTIONS
    CALL dynit(ekincp,ekin1,ekin2,temp1,temp2,ekinh1,ekinh2)
    ! ..PARAMETERS FOR THE NOSE-HOOVER THERMOSTATS
    ncdof = 9
    IF (prcpl%tisot.OR.prcpl%tzflex) ncdof = 1
    IF ((cntl%tnosep.OR.cntl%tnosec).AND.paral%parent) CALL nosepa(np_low,np_high)
    ! thermostat
    reset_gkt=.FALSE.
    IF (ALLOCATED(gkt1)) reset_gkt=.TRUE.
    pitmnose(:)=nosl%tmnose
    IF (nosl%tmnose.AND.pimd1%tcentro) pitmnose(1)=.FALSE. ! massive -> single for ip=1
    ! barostat
    pitnosec(:)=cntl%tnosec
    IF (cntl%tnosec) pitnosec(2)=.FALSE. ! switch off for ip>1
    ! Dont symmetrize density
    cntl%tsymrho=.FALSE.
    ! ..INITIALIZATION OF WAVEFUNCTION AND COORDINATS
    IF (cprint%iprint_step.EQ.0) cprint%iprint_step=cnti%nomore+1
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
          CALL zhrwf(1,irec(:,ipcurr),c0(:,:,ipx:ipx),c2,crge%n,eigv(:,ipx),&
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
          IF (irec(irec_rho,ipcurr).NE.0) THEN
             CALL dcopy(nnx,rhoo,1,rall(1,ipx),1)
          ENDIF
          ! ..Randomization of the atomic coordinates
          IF (cntl%tranp.AND.paral%parent) CALL ranp(pitaup(:,:,:,ipcurr))
          ! ..Initialization of wavefunction
          IF (irec(irec_wf,ipcurr).EQ.0) THEN
             CALL rinitwf(c0(:,:,ipx:ipx),c2,sc0,crge%n,pitaup(:,:,:,ipcurr),taur,rhoe,psi)
          ENDIF
          ! Randomize electronic wavefunctions around loaded ones
          IF (cntl%trane) CALL rrane(c0(:,:,ipx:ipx),c2,crge%n)
          CALL phfac(pitaup(:,:,:,ipcurr))
          IF (corel%tinlc) CALL copot(rhoe,psi,.TRUE.)
          ! Orthogonalize electronic wavefunctions
          IF (irec(irec_wf,ipcurr).EQ.0.OR.cntl%trane.OR.(pslo_com%tivan.AND..NOT.cntl%nonort)) THEN
             IF (pslo_com%tivan) CALL rnlsm(c0(:,:,ipx),crge%n,1,1,.FALSE.)
             CALL ortho(crge%n,c0(:,:,ipx),c2)
          ENDIF
          ! ..Store occupation numbers in fall
          CALL dcopy(crge%n*nkpt%nkpts,crge%f,1,fall(1,ipx),1)
          ! ..Store WF history in coldall
          IF (cntl%textrap.AND.np_local>1) THEN
             CALL zcopy(nkpt%ngwk*crge%n*nkpt%nkpnt*cnti%mextra,cold,1,&
                  coldall(1,1,1,1,ipx),1)
             nnowall(ipx)=nnow; numcoldall(ipx)=numcold
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
          CALL rinitwf(c0(:,:,ipx:ipx),c2,sc0,crge%n,pitaup(:,:,:,ipcurr),taur,rhoe,psi)
          IF (corel%tinlc) CALL copot(rhoe,psi,.TRUE.)
          IF (pslo_com%tivan) CALL rnlsm(c0(:,:,ipx),crge%n,1,1,.FALSE.)
          CALL ortho(crge%n,c0(:,:,ipx),c2)
       ENDIF
       ! ..Initialize velocities
       IF (.NOT.restart1%rvel) THEN
          ener_com%ecnstr = 0.0_real_8
          ener_com%erestr = 0.0_real_8
          IF (pimd1%tpinm.OR.pimd1%tstage) THEN
             CALL s_rinvel(vstage(:,:,:,ipcurr),c2,crge%n,ipcurr)
             IF (paral%parent.AND.ipcurr.EQ.1) CALL rattle(stagep(:,:,:,1),vstage(:,:,:,1))
             IF (ipcurr.EQ.1) CALL rvscal(vstage(:,:,:,ipcurr))
          ELSE
             CALL rinvel(pivelp(:,:,:,ipcurr),c2,crge%n)
          ENDIF
       ENDIF
       IF (cntl%quenchp) THEN
          IF (pimd1%tpinm.OR.pimd1%tstage) THEN
             CALL zeroing(vstage(:,:,:,ipcurr))!,3*maxsys%nax*maxsys%nsx)
          ELSE
             CALL zeroing(pivelp(:,:,:,ipcurr))!,3*maxsys%nax*maxsys%nsx)
          ENDIF
       ENDIF
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
    IF (cntl%tdiag) THEN
       IF (cntl%tlanc) nx=1
       IF (cntl%tdavi) nx=nkpt%ngwk*cnti%ndavv*nkpt%nkpnt+1
    ELSEIF (cntl%tsde) THEN
       nx=1
    ELSEIF (cntl%diis) THEN
       nx=(nkpt%ngwk*crge%n+8)*cnti%mdiis/2+4
    ELSEIF (cntl%pcg) THEN
       nx=1
    ENDIF
    IF (grandparent) THEN
       IF (paral%io_parent)&
            WRITE(6,'(1X,64("="))')
       IF (paral%io_parent)&
            WRITE(6,'(1X,"==",T25,A,T64,"==")')&
            'FORCES INITIALIZATION'
       IF (paral%io_parent)&
            WRITE(6,'(1X,64("="))')
    ENDIF
    ifcalc=0
    DO ipcurr=np_low,np_high
       ipx=ipcurr-np_low+1
       CALL phfac(pitaup(:,:,:,ipcurr))
       IF (corel%tinlc) CALL copot(rhoe,psi,.TRUE.)
       ! Initialization of density and potential
       ! for diagonalization schemes
       IF (cntl%tdiag) THEN
          IF (pslo_com%tivan) THEN
             CALL rnlsm(c0(:,:,ipx),crge%n,1,1,.FALSE.)
          ENDIF
          IF (irec(irec_rho,ipcurr).EQ.0) THEN
             ! ..Restore occupation numbers
             IF (restart1%restart) THEN
                CALL dcopy(crge%n*nkpt%nkpts,fall(1,ipx),1,crge%f,1)
             ENDIF
             IF (tkpts%tkpnt) THEN
                CALL rhoofr_c(c0(:,:,ipx),rhoe,psi(:,1),crge%n)
             ELSE
                CALL rhoofr(c0(:,:,ipx),rhoe,psi(:,1),crge%n)
             ENDIF
             CALL dcopy(nnx,rhoe,1,rin0,1)
          ELSE
             CALL dcopy(nnx,rall(1,ipx),1,rin0,1)
          ENDIF
          ! NON LOCAL PROJECTOR OVERLAP MATRIX
          IF (fint1%ttrot) CALL calc_alm
       ENDIF
       ! INITIALIZE WF CENTERS & SPREAD
       IF (vdwl%vdwd) THEN
          IF (.NOT.trwanncx(ipx)) THEN
             CALL localize2(pitaup(:,:,:,ipcurr),c0(:,:,ipx),c2,sc0,crge%n)
          ENDIF
       ENDIF
       CALL forces_diag(crge%n,c0(:,:,ipx:),c2,cm,sc0,cm(nx),vpp,eigv(:,ipx),&
            rhoe,psi,&
            pitaup(:,:,:,ipcurr),pivelp(:,:,:,ipcurr),taui(:,:,:,ipcurr),&
            pifion(:,:,:,ipcurr),&
            ifcalc,irec(:,ipcurr),.TRUE.,.FALSE.)
       CALL dscal(3*maxsys%nax*maxsys%nsx,1._real_8/rnp,pifion(1,1,1,ipcurr),1)
       CALL dcopy(3*maxsys%nax*maxsys%nsx,pifion(1,1,1,ipcurr),1,&
            fionks(1,1,1,ipcurr),1)
       CALL fharm(pitaup,pifion(:,:,:,ipcurr),ipcurr,.TRUE.)
       CALL pi_ener(ipcurr,0)
       IF (ipcurr.EQ.np_low) CALL freqs(crge%n,.FALSE.)
       CALL dcopy(crge%n*nkpt%nkpts,crge%f,1,fall(1,ipx),1)
       IF (cntl%tdiag) CALL dcopy(nnx,rin0,1,rall(1,ipx),1)
       CALL totstr
       CALL dcopy(9,paiu(1,1),1,paiuall(1,1,ipcurr),1)
       CALL dscal(9,1._real_8/rnp,paiuall(1,1,ipcurr),1)
    ENDDO
    ! sync all replica
    CALL mp_sync(supergroup)

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
          IF (cntl%tnosec .AND. irec(irec_noc,ipcurr).EQ.0) CALL noscinit(ipcurr)
          ! ..Centroid & Ring-polymer PIMD: Reset thermostat variables in case of
          ! ..a centroid restart from a non-centroid thermosatted run
          IF (cntl%tnosep.AND.(pimd1%tcentro.OR.pimd1%tringp).AND.ipcurr.EQ.1) THEN
             nrepl=maxnp
             IF (nosl%tultra) THEN
                ! ..Not used with PICPMD 
                CALL stopgm(procedureN,'Ultra NHCs thermostat not supported',&
                     __LINE__,__FILE__)
             ELSEIF (loct%tloct)  THEN
                CALL stopgm(procedureN,'LOCAL TEMPERATURE NOT IMPLEMENTED',& 
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
       CALL getgyration(pitaup,pivelp) !vw routine requies 2 dummy args, added pivelp
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
       CALL prmem('   PI_NPT_BOMD')
       time2 =m_walltime()
       tcpu = (time2 - time1)*0.001_real_8
       IF (paral%io_parent)&
            WRITE(6,'(A,F8.2,A)') ' TIME FOR INITIALIZATION',&
            tcpu,' SECONDS'
    ENDIF
    CALL mp_sync(supergroup)
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
             CALL anneal(vstage(:,:,:,ipcurr),c2,crge%n,metr_com%htvel)
             CALL berendsen(vstage(:,:,:,ipcurr),c2,crge%n,metr_com%htvel,0.0_real_8,ekinh)
          ELSE
             CALL anneal(pivelp(:,:,:,ipcurr),c2,crge%n,metr_com%htvel)
             CALL berendsen(pivelp(:,:,:,ipcurr),c2,crge%n,metr_com%htvel,0.0_real_8,ekinh)
          ENDIF
          ! ..UPDATE NOSE THERMOSTATS
          nosl%tmnose=pitmnose(MIN(ipcurr,2))
          cntl%tnosec=pitnosec(MIN(ipcurr,2))
          IF (reset_gkt.AND.ipcurr==1) THEN
             scr(1:2)=gkt(1:2,1); gkt(1:2,1)=gkt1(1:2,1)
          ENDIF
          IF (pimd1%tstage.OR.pimd1%tpinm) THEN
             CALL noseup(vstage(:,:,:,ipcurr),c2,crge%n,ipcurr)
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
                ! IONS & CELL
                CALL s_velupi(vstage(:,:,:,ipcurr),fstage(:,:,:,ipcurr),1,ipcurr)
                IF (ipcurr.EQ.1) CALL velupc(1)
             ELSE
                CALL stopgm(procedureN,'not yet implemented',&
                     __LINE__,__FILE__)
             ENDIF
          ENDIF
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
       ! <<<<End of replica loop
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
          IF (corel%tinlc) CALL copot(rhoe,psi,.TRUE.)
          IF (cntl%tdiag) CALL dcopy(nnx,rall(1,ipx),1,rin0,1)
          IF (tmovr) THEN
             CALL dcopy(nnx,rin0,1,rhoe,1)
             CALL moverho(rhoe,psi)
             CALL dcopy(nnx,rhoe,1,rin0,1)
          ENDIF
          IF (fint1%ttrot) CALL calc_alm
          ! Wannier stuff for vdW-WC
          IF (vdwl%vdwd) THEN
             twannupx(ipx)=twannupx(ipx).OR.(infi.EQ.1.AND..NOT.trwanncx(ipx))
             IF (twannupx(ipx)) THEN
                CALL localize2(pitaup(:,:,:,ipcurr),c0(:,:,ipx),c2,sc0,crge%n)
             ENDIF
          ENDIF
          IF (cntl%textrap) THEN
             ! Extrapolate wavefunctions
             IF (np_local>1) THEN
                CALL dcopy(crge%n*nkpt%nkpts,fall(1,ipx),1,crge%f,1)
                CALL zcopy(nkpt%ngwk*crge%n*nkpt%nkpnt*cnti%mextra,&
                     coldall(1,1,1,1,ipx),1,cold,1)
                nnow=nnowall(ipx); numcold=numcoldall(ipx)
             ENDIF
             CALL extrapwf(infi,c0(:,:,ipx),scr,cold,nnow,numcold,crge%n,cnti%mextra)
             IF (np_local>1) THEN
                CALL zcopy(nkpt%ngwk*crge%n*nkpt%nkpnt*cnti%mextra,&
                     cold,1,coldall(1,1,1,1,ipx),1)
                nnowall(ipx)=nnow; numcoldall(ipx)=numcold
             ENDIF
          ENDIF
          IF (cntl%tlanc) nx=1
          IF (cntl%tdavi) nx=cnti%ndavv*nkpt%nkpnt+1
          ! ..CALCULATE THE FORCES
          CALL forces_diag(crge%n,c0(:,:,ipx:),c2,cm,sc0,cm(nx),&
               vpp,eigv(:,ipx),rhoe,psi,&
               pitaup(:,:,:,ipcurr),pivelp(:,:,:,ipcurr),taui(:,:,:,1),&
               pifion(:,:,:,ipcurr),&
               ifcalc,irec(:,ipcurr),.TRUE.,.FALSE.)
          CALL dscal(3*maxsys%nax*maxsys%nsx,1._real_8/rnp,pifion(1,1,1,ipcurr),1)
          CALL dcopy(3*maxsys%nax*maxsys%nsx,pifion(1,1,1,ipcurr),1,&
               fionks(1,1,1,ipcurr),1)
          CALL fharm(pitaup,pifion(:,:,:,ipcurr),ipcurr,.TRUE.)
          IF (cntl%tdiag) CALL dcopy(nnx,rin0,1,rall(1,ipx),1)
          CALL dcopy(crge%n*nkpt%nkpts,crge%f,1,fall(1,ipx),1)
          CALL pi_ener(ipcurr,0)
          CALL totstr
          CALL dcopy(9,paiu(1,1),1,paiuall(1,1,ipcurr),1)
          CALL dscal(9,1._real_8/rnp,paiuall(1,1,ipcurr),1)
       ENDDO
       CALL mp_sync(supergroup)
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
          ! ..IONIC TEMPERATURE CONTROL
          IF (paral%parent.AND.(cntl%tcp.OR.cntl%trescale)) THEN
             IF (pimd1%tpinm.OR.pimd1%tstage) THEN
                CALL s_ekinpp(ekinp,vstage(:,:,:,ipcurr),ipcurr)
                ! ..GLIB is calculated in QMDOF for the first (staging) bead 
                ! ..that is possible constraint and is 3*NAT for all other beads
                glib_s=3._real_8*ions1%nat
                IF (ipcurr.EQ.1) glib_s=glib
                tempp(ipcurr)=ekinp*factem*2._real_8/glib_s
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
             ELSE
                CALL ekinpp(ekinp,pivelp(:,:,:,ipcurr))
                tempp(ipcurr)=ekinp*factem*2._real_8/glib
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
          nosl%tmnose=pitmnose(MIN(ipcurr,2))
          cntl%tnosec=pitnosec(MIN(ipcurr,2))
          IF (reset_gkt.AND.ipcurr==1) THEN
             scr(1:2)=gkt(1:2,1); gkt(1:2,1)=gkt1(1:2,1)
          ENDIF
          IF (pimd1%tstage.OR.pimd1%tpinm) THEN
             CALL noseup(vstage(:,:,:,ipcurr),c2,crge%n,ipcurr)
          ELSE
             CALL stopgm(procedureN,'not yet implemented',&
                  __LINE__,__FILE__)
          ENDIF
          IF (reset_gkt.AND.ipcurr==1) gkt(1:2,1)=scr(1:2)
          ekinh=0.5_real_8*cntr%cmass*ddot(9,metr_com%htvel,1,metr_com%htvel,1)
          IF (pimd1%tpinm.OR.pimd1%tstage) THEN
             CALL berendsen(vstage(:,:,:,ipcurr),c2,crge%n,metr_com%htvel,0.0_real_8,ekinh)
             ! ANNEALING
             CALL anneal(vstage(:,:,:,ipcurr),c2,crge%n,metr_com%htvel)
          ELSE
             CALL berendsen(pivelp(:,:,:,ipcurr),c2,crge%n,metr_com%htvel,0.0_real_8,ekinh)
             ! ANNEALING
             CALL anneal(pivelp(:,:,:,ipcurr),c2,crge%n,metr_com%htvel)
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
          ! ..MEAN SQUARE DISPLACEMENT OF DIFFERENT IONIC SPECIES
          IF (paral%parent) CALL dispp(pitaup(:,:,:,ipcurr),taui(:,:,:,ipcurr),&
               disa(ipcurr))
          IF (paral%parent) THEN
             ! ..ENERGY OF THE NOSE THERMOSTATS
             nosl%tmnose=pitmnose(MIN(ipcurr,2))
             cntl%tnosec=pitnosec(MIN(ipcurr,2))
             IF (reset_gkt.AND.ipcurr==1) THEN
                scr(1:2)=gkt(1:2,1); gkt(1:2,1)=gkt1(1:2,1)
             ENDIF
             IF (pimd1%tpinm.OR.pimd1%tstage) THEN
                CALL noseng(iteropt%nfi,vstage(:,:,:,ipcurr),dummy,enosp(ipcurr),enosc(ipcurr),ipcurr)
             ELSE
                CALL stopgm(procedureN,'not yet implemented',&
                     __LINE__,__FILE__)
             ENDIF
             IF (reset_gkt.AND.ipcurr==1) gkt(1:2,1)=scr(1:2)
             ! ..EFFECTIVELY CONSERVED AND HAMILTONIAN ENERGY
             ! ..NOTE: ECNSTR IS ALWAYS ZERO
             econs(ipcurr)=ekinp+etotv(ipcurr)/rnp+&
                  enosp(ipcurr)+enosc(ipcurr)+&
                  eharv(ipcurr)+ener_com%ecnstr+ener_com%erestr
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
       CALL global(enosp,1)
       CALL global(enosc,1)
       CALL global(disa,1)
       ! DM C..PRINTOUT the evolution of the accumulators every time step
       IF (grandparent) THEN
          ! ..Radii of gyration
          CALL getgyration(pitaup, pivelp) !vw routine requies 2 dummy args, added pivelp
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
                  '      NFI   IP  TEMPP     EHARM         EKS',&
                  '    ENOSE(C)    ENOSE(P)    ECLASSIC     DIS'
             DO ip=1,pimd3%np_total
                IF (paral%io_parent)&
                     WRITE(6,&
                     '(I9,I5,F7.1,F10.6,F12.5,2(2X,F10.5),F12.5, F8.2)')&
                     iteropt%nfi,ip,tempp(ip),eharv(ip),etotv(ip),&
                     enosc(ip),enosp(ip),econs(ip),disa(ip)
             ENDDO
          ENDIF
          ! ..Quantum kinetic energy of ions
          ! ..Primitive (Barker) estimator
          ! DM C..Barker estimator
          qpkinp=0.0_real_8
          DO ip=1,pimd3%np_total
             qpkinp=qpkinp+eharv(ip)
          ENDDO
          qpkinp=(glib+3._real_8*ions1%nat*(rnp-1))*cntr%tempw/(2._real_8*factem)-qpkinp
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
             etota=etota+etotv(ip)/rnp
             econsa=econsa+econs(ip)
             tempa=tempa+tempp(ip)/rnp
          ENDDO
          equant=etota+qpkinp
          econsa=econsa+prcp_com%druck*parm%omega
          ehama=econsa+ekinh
          time2=m_walltime()
          tcpu=(time2-time1)*0.001_real_8
          IF ((ropt_mod%engpri).AND.paral%io_parent)&
               WRITE(6,'(/,A,A)')&
               '      NFI    EKINH  TEMP  EKINP(P)  EKINP(V)       EKS/P',&
               '       EQUANT     ECLASSIC         EHAM     TCPU'
          IF (paral%io_parent)&
               WRITE(6,'(I9,F9.5,F7.1,2F9.5,4F13.5,F9.2)')&
               iteropt%nfi,ekinh,tempa,qpkinp,qvkinp,etota,equant,econsa,ehama,tcpu
          IF (paral%io_parent)&
               WRITE(3,'(I9,F9.5,F7.1,2F9.5,4F20.10,F9.2)')&
               iteropt%nfi,ekinh,tempa,qpkinp,qvkinp,etota,equant,econsa,ehama,tcpu
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

       nosl%tmnose=pitmnose(2)  ! recover original setup of tmnose
       IF (pimd1%tstage.OR.pimd1%tpinm) CALL wr_temps(iteropt%nfi,ekinc,tempp)
       ! ..Recompute velocity-dependent part of pipaiu
       IF (pimd1%tstage.OR.pimd1%tpinm) THEN
          pipaiu=pi_stress_vir(paiuall,stagep(:,:,:,1),vstage(:,:,:,1),pitaup,fionks,.TRUE.)
       ELSE
          CALL stopgm(procedureN,'not yet implemented',&
               __LINE__,__FILE__)
       ENDIF
       IF (grandparent) CALL wr_pi_stress(pipaiu)
 
       IF (grandparent) THEN
          DO ipcurr=np_low,np_high
             ! ..Update accumulators
             i=ipcurr
             CALL pacce(ekinc(i),etotv(i),econs(i),eham(i),tempp(i),&
                  eharv(i),enose(i),enosp(i),disa(i),&
                  tcpu,qpkinp,qvkinp,ekina,etota,econsa,&
                  ehama,tempa,equant,temph,parm%omega,accus(:,i),iteropt%nfi,1)
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
          IF (teststore(iteropt%nfi).OR.soft_com%exsoft) THEN
             IF (cntl%tdiag) CALL dcopy(nnx,rall(1,ipx),1,rhoe,1)
             CALL dcopy(crge%n*nkpt%nkpts,fall(1,ipx),1,crge%f,1)
             IF (cntl%textrap.AND.np_local>1) THEN
                CALL zcopy(nkpt%ngwk*crge%n*nkpt%nkpnt*cnti%mextra,&
                     coldall(1,1,1,1,ipx),1,cold,1)
                nnow=nnowall(ipx); numcold=numcoldall(ipx)
             ENDIF
             nosl%tmnose=pitmnose(MIN(ipcurr,2))
             cntl%tnosec=pitnosec(MIN(ipcurr,2))
             CALL write_irec(irec(:,ipcurr))
             CALL zhwwf(2,irec(:,ipcurr),c0(:,:,ipx),c2,crge%n,eigv(:,ipx),&
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
    ! NEW  What do we have to do with this call??
    ! NEW  IF(RHOOUT) CALL RHOPRI(C0,TAU0,RHOE,PSI,N,1)
    IF (grandparent) THEN
       ! Dummy arguments
       d1=0._real_8
       d2=0._real_8
       CALL pacce(d1,d1,d1,d1,d1,d1,d1,d1,d1,&
            d2,d2,d2,d2,d2,d2,d2,d2,d2,d2,d2,accus,iteropt%nfi,2)
    ENDIF
    IF (paral%parent) CALL prmem('   PI_NPT_BOMD')
    DEALLOCATE(eigv,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(fionks,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(fall,STAT=ierr)
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
    DEALLOCATE(rhoe,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(psi,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(scr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rin0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rout0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rmix,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rm1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rinp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (cntl%tdiag) THEN
       DEALLOCATE(rall,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    IF (cntl%textrap) THEN
       DEALLOCATE(cold,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       IF (np_local>1) THEN
          DEALLOCATE(coldall,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
          DEALLOCATE(nnowall,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
          DEALLOCATE(numcoldall,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
       ENDIF
    ENDIF
    IF (comvl%tsubrot) THEN
       DEALLOCATE(tauio,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
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
    RETURN
  END SUBROUTINE pi_npt_bomd
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

END MODULE pi_npt_bomd_utils
