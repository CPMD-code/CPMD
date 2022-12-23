MODULE pi_diag_utils
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
  USE glemod,                          ONLY: glepar,&
                                             glec,&
                                             gles,&
                                             glet,&
                                             glep
  USE gle_utils,                       ONLY: gle_init,&
                                             gle_step
  USE global_utils,                    ONLY: global
  USE ions,                            ONLY: ions1
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE localize_utils,                  ONLY: localize2
  USE machine,                         ONLY: m_walltime
  USE md_driver,                       ONLY: give_scr_mddiag
  USE mm_extrap,                       ONLY: cold,&
                                             nnow,&
                                             numcold
  USE moverho_utils,                   ONLY: moverho
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sync
! USE interface_utils,                 ONLY: get_external_forces
  USE mts_utils,                       ONLY: mts, set_mts_functional
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
                                             pi_ener,&
                                             gleback
  USE pimd,                            ONLY: &
       eharv, etotv, fionks, grandparent, ipcurr, maxnp, np_high, np_local, &
       np_low, pimd1, pimd3, supergroup, pma0s, pi_egle, pi_glec, pi_gles, &
       pi_glet, pi_glep
  USE pinmtrans_utils,                 ONLY: pinmtrans
  USE pi_stress_utils,                 ONLY: pi_stress_vir,&
                                             pi_stress_nmm,&
                                             pi_stress_pri,&
                                             wr_pi_stress
  USE pi_wf_utils,                     ONLY: initrep
  USE poin,                            ONLY: rhoo
  USE posupi_utils,                    ONLY: posupi
  USE printave_utils,                  ONLY: paccd
  USE printp_utils,                    ONLY: printp,&
                                             print_mts_forces
  USE prmem_utils,                     ONLY: prmem
  USE prng,                            ONLY: prng_com,&
                                             pi_prng_com
  USE prng_utils,                      ONLY: prnginit,&
                                             repprngu_vec
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
  USE rmas,                            ONLY: rmass
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
  USE shake_utils,                     ONLY: cpmdshake
  USE soft,                            ONLY: soft_com
  USE spin,                            ONLY: clsd
  USE stagetrans_utils,                ONLY: stagetrans
  USE store_types,                     ONLY: &
       cprint, iprint_coor, iprint_force, irec_nop1, irec_nop2, irec_nop3, &
       irec_nop4, irec_rho, irec_wf, irec_prng, irec_co, restart1, rout1
  USE strs,                            ONLY: paiu
  USE system,                          ONLY: &
       acc, cnti, cntl, cntr, fpar, maxsys, nacc, nacx, nkpt, restf, ncpw, iatpt
  USE testex_utils,                    ONLY: testex,&
                                             testex_mw
  USE teststore_utils,                 ONLY: teststore
  USE totstr_utils,                    ONLY: totstr
  USE vdw_utils,                       ONLY: rwannier
  USE vdwcmod,                         ONLY: &
       nfragx, nwfcx, rwann, swann, tauref, trwanncx, twannupx, vdwl, vdwwfl
  USE velupi_utils,                    ONLY: s_velupi,&
                                             velupi
  USE wr_temps_utils,                  ONLY: wr_temps
  USE wrener_utils,                    ONLY: wrener
  USE wrgeo_utils,                     ONLY: wrgeo,&
                                             wrgeof
  USE wv30_utils,                      ONLY: zhwwf
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: pi_diag

CONTAINS

  ! ==================================================================
  SUBROUTINE pi_diag(c0,cm,c2,sc0,vpp,gamx,gamy)
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8)                          :: c0(:,:,:), cm(*), c2(*), &
                                                sc0(*)
    REAL(real_8)                             :: vpp(*), &
                                                gamx(crge%n*crge%n,*), &
                                                gamy(crge%n*crge%n,*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'pi_diag'

    CHARACTER(len=100)                       :: filen
    CHARACTER(len=12)                        :: cflbod, cipnum
    CHARACTER(len=30)                        :: tag
    COMPLEX(real_8), ALLOCATABLE             :: coldall(:,:,:,:,:), &
                                                coldall_high(:,:,:,:,:), &
                                                psi(:,:)
    INTEGER :: i, i1, i2, ierr, ifcalc, ik, il_psi_1d, il_psi_2d, il_rhoe_1d, &
      il_rhoe_2d, ip, ipx, itemp, k, l, lenext, lscr, n1, n2, nnx, &
      npx, nrepl, nwfc, nx, iprng0
    INTEGER, ALLOCATABLE                     :: nnowall(:), numcoldall(:), &
                                                nnowall_high(:), numcoldall_high(:), &
                                                irec(:,:)
    LOGICAL                                  :: ferror, reset_gkt, &
                                                pitmnose(2)
    REAL(real_8) :: accus(nacx,maxnp), d1, d2, disa(maxnp), dummies(maxnp), &
      dummy, dummy2, econs(maxnp), econsa, eham(maxnp), ehama, ekin1, ekin2, ekina, &
      ekinc(maxnp), ekincp, ekinh1, ekinh2, ekinp, enose(maxnp), &
      enosp(maxnp), equant, etota, glib_s, lmio(3), qpkinp, qvkinp, rmem, &
      rnp, tcpu, temp1, temp2, tempa, tempp(maxnp), time1, time2, vcmio(4), &
      pipaiu(3,3)
    REAL(real_8), ALLOCATABLE :: eigv(:,:), fall(:,:), fstage(:,:,:,:), &
      rall(:,:), rhoe(:,:), rinp(:), rm1(:), scr(:), stagep(:,:,:,:), &
      taui(:,:,:,:), tauio(:,:), taur(:,:,:), vstage(:,:,:,:), paiuall(:,:,:)
    REAL(real_8), POINTER                    :: pifion(:,:,:,:), &
                                                pitaup(:,:,:,:), &
                                                pivelp(:,:,:,:)
    ! MTS[
    ! number of inner steps between two large steps and total number of large steps
    integer :: n_inner_steps, n_large_steps
    ! logical to know if the current step is a large step in the MTS scheme
    logical :: mts_large_step, mts_pure_dft
    ! high level ionic forces
    real(real_8), allocatable :: fion_high(:,:,:)
    ! high level wave-function parameters
    complex(real_8), allocatable :: c0_high(:,:,:)
    ! MTS]

    accus=0.0_real_8
    disa=0.0_real_8
    econs=0.0_real_8
    eham=0.0_real_8
    ekinc=0.0_real_8
    enose=0.0_real_8
    enosp=0.0_real_8
    tempp=0.0_real_8
    restf%nfnow=1
    time1 =m_walltime()

    mts_pure_dft = .false.
    if (cntl%use_mts) then

       ! allocate high level forces array
       allocate(fion_high(3,maxsys%nax,maxsys%nsx),stat=ierr)
       if(ierr/=0) call stopgm(proceduren,'allocation problem: fion_high',&
          __LINE__,__FILE__)
       call zeroing(fion_high)

       ! allocate high level WF param array only if needed
       mts_pure_dft = ( (mts%low_level == 'DFT') .and. (mts%high_level == 'DFT') )
       if (mts_pure_dft) then
          allocate( c0_high(size(c0,1),size(c0,2),size(c0,3)), stat=ierr )
          if(ierr/=0) call stopgm(proceduren,'allocation problem : c0_high',&
             __LINE__,__FILE__)
       end if

    end if

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
    IF (cntl%tpres) THEN
       ALLOCATE(paiuall(3,3,npx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ropt_mod%calste=cntl%tpres
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

       ! allocate array for high level WF extrapolation
       if (mts_pure_dft) then
          allocate(coldall_high(nkpt%ngwk,crge%n,nkpt%nkpnt,cnti%mextra,np_local),stat=ierr)
          if(ierr/=0) call stopgm(proceduren,'allocation problem: coldall_high',&
               __LINE__,__FILE__)
          call zeroing(coldall_high)
          allocate(nnowall_high(np_local),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
          call zeroing(nnowall_high)
          allocate(numcoldall_high(np_local),STAT=ierr)
          if(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
          call zeroing(numcoldall_high)
          rmem=rmem+rmem*REAL(np_local,kind=real_8)
       endif

       IF (paral%io_parent)&
            WRITE(6,'(A,T51,F8.3,A)') ' PI_DIAG| '&
            // 'EXTRAPOLATION WAVEFUNC. HISTORY TAKES ',rmem,' MBYTES'
    ENDIF

    ALLOCATE(pi_egle(npx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    IF (np_local>1) THEN
       ALLOCATE(pi_prng_com(np_local),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       IF (glepar%gle_mode>0) THEN
          ALLOCATE(pi_glec(np_local),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(pi_gles((glepar%gle_ns+1)*(glepar%gle_ns+1),np_local),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(pi_glet((glepar%gle_ns+1)*(glepar%gle_ns+1),np_local),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(pi_glep(3,maxsys%nax,maxsys%nsx,(glepar%gle_ns+1),np_local),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
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
    CALL give_scr_mddiag(lscr,tag)
    ALLOCATE(scr(lscr),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ..TIME STEP FUNCTIONS
    CALL dynit(ekincp,ekin1,ekin2,temp1,temp2,ekinh1,ekinh2)
    ! ..PARAMETERS FOR THE NOSE-HOOVER THERMOSTATS
    IF (cntl%tnosep.AND.paral%parent) CALL nosepa(np_low,np_high)
    ! thermostat
    reset_gkt=.FALSE.
    IF (ALLOCATED(gkt1)) reset_gkt=.TRUE.
    pitmnose(:)=nosl%tmnose
    IF (nosl%tmnose.AND.pimd1%tcentro) pitmnose(1)=.FALSE. ! massive -> single for ip=1
    ! Dont symmetrize density
    cntl%tsymrho=.FALSE.
    ! ..INITIALIZATION OF WAVEFUNCTION AND COORDINATS
    IF (cprint%iprint_step.EQ.0) cprint%iprint_step=cnti%nomore+1
    CALL zeroing(acc)!,nacc)
    DO ipcurr=np_low,np_high
       nosl%tmnose=pitmnose(MIN(ipcurr,2))
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
             restart1%rvel=.TRUE.
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
          IF (corel%tinlc) CALL copot(rhoe,psi,ropt_mod%calste)
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
          ! ..Store GLE-related variables
          CALL gleback(ipcurr,.FALSE.,0)
       ENDIF
    ENDDO
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
          IF (corel%tinlc) CALL copot(rhoe,psi,ropt_mod%calste)
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
       IF (corel%tinlc) CALL copot(rhoe,psi,ropt_mod%calste)
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
       n_inner_steps=0
       n_large_steps=0
       mts_large_step=.false.
       if (cntl%use_mts) then
          n_inner_steps=MOD(iteropt%nfi,mts%timestep_factor)
          n_large_steps=iteropt%nfi/mts%timestep_factor
          if (n_inner_steps==0) mts_large_step=.true.
          call get_mts_forces(mts_large_step,.true.,n_inner_steps,n_large_steps,&
             mts_pure_dft,fion_high,c0_high,coldall_high(:,:,:,:,ipx),&
             nnowall_high(ipx),numcoldall_high(ipx),infi,&
             crge%n,c0(:,:,ipx:),c2,cm,sc0,cm(nx),vpp,eigv(:,ipx),&
             rhoe,psi,&
             pitaup(:,:,:,ipcurr),pivelp(:,:,:,ipcurr),taui(:,:,:,ipcurr),&
             pifion(:,:,:,ipcurr),&
             ifcalc,irec(:,ipcurr),.true.,.true.)
       else
          CALL forces_diag(crge%n,c0(:,:,ipx:),c2,cm,sc0,cm(nx),vpp,eigv(:,ipx),&
               rhoe,psi,&
               pitaup(:,:,:,ipcurr),pivelp(:,:,:,ipcurr),taui(:,:,:,ipcurr),&
               pifion(:,:,:,ipcurr),&
               ifcalc,irec(:,ipcurr),.TRUE.,.TRUE.)
       endif
       CALL dscal(3*maxsys%nax*maxsys%nsx,1._real_8/rnp,pifion(1,1,1,ipcurr),1)
       CALL dcopy(3*maxsys%nax*maxsys%nsx,pifion(1,1,1,ipcurr),1,&
            fionks(1,1,1,ipcurr),1)
       CALL fharm(pitaup,pifion(:,:,:,ipcurr),ipcurr,.TRUE.)
       CALL pi_ener(ipcurr,0)
       IF (ipcurr.EQ.np_low) CALL freqs(crge%n,.FALSE.)
       CALL dcopy(crge%n*nkpt%nkpts,crge%f,1,fall(1,ipx),1)
       IF (cntl%tdiag) CALL dcopy(nnx,rin0,1,rall(1,ipx),1)
    ENDDO
    ! sync all replica
    CALL mp_sync(supergroup)

    IF (paral%parent) THEN
       CALL global(fionks,3*maxsys%nax*maxsys%nsx)
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
          itemp=irec(irec_nop1,ipcurr)+irec(irec_nop2,ipcurr)+irec(irec_nop3,ipcurr)&
               +irec(irec_nop4,ipcurr)
          IF (cntl%tnosep .AND. itemp.EQ.0) CALL nospinit(ipcurr)
          ! ..Centroid & Ring-polymer PIMD: Set thermostat variables to zero in case 
          ! ..of a centroid restart from a non-centroid thermosatted run
          IF (cntl%tnosep.AND.(pimd1%tcentro.OR.pimd1%tringp).AND.ipcurr.EQ.1) THEN
             nrepl=maxnp
             IF (nosl%tultra) THEN
                ! ..Not used with PICPMD 
             ELSEIF (loct%tloct)  THEN
                CALL stopgm('PI_DIAG','LOCAL TEMPERATURE NOT IMPLEMENTET',& 
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
       ENDDO
    ENDIF

    ! Initialize GLE thermostat
    iprng0=cnti%iprng
    CALL repprngu_vec(npx,dummies)
    DO ipcurr=np_low,np_high
       ipx=ipcurr-np_low+1
       CALL gleback(ipcurr,.FALSE.,1)
       IF (.NOT.restart1%rprng.OR.irec(irec_prng,ipcurr).EQ.0) THEN
          cnti%iprng=FLOOR(dummies(ipcurr)*iprng0)
          IF (paral%io_parent) WRITE(6,*) 'USING SEED ', cnti%iprng,&
             'TO REINIT. PSEUDO RANDOM NUMBER GEN.'
          CALL prnginit
       ENDIF
       IF (pimd1%tstage.OR.pimd1%tpinm) THEN
          CALL gle_init(stagep(:,:,:,ipcurr),vstage(:,:,:,ipcurr),pma0s(1,ipcurr))
       ELSE
          CALL gle_init(pitaup(:,:,:,ipcurr),pivelp(:,:,:,ipcurr),rmass%pma)
       ENDIF
       CALL gleback(ipcurr,.TRUE.,0)
    ENDDO

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
       ! ..END INITIALIZATION
       CALL prmem('   PI_DIAG')
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

       ! MTS time step counters
       n_inner_steps=n_inner_steps+1
       mts_large_step=.false.
       ! check if it is a large step
       if (n_inner_steps==mts%timestep_factor) then
          mts_large_step=.true.
          n_large_steps=n_large_steps+1
          n_inner_steps=0
       endif

       ropt_mod%prteig=MOD(iteropt%nfi-1,cprint%iprint_step).EQ.0
       IF (.NOT.paral%parent) ropt_mod%prteig=.FALSE.
       ropt_mod%engpri=MOD(iteropt%nfi-1,cprint%iprint_step).EQ.0
       comvl%subcom=comvl%tsubcom.AND.MOD(iteropt%nfi-1,comvl%ncomv).EQ.0
       comvl%subrot=comvl%tsubrot.AND.MOD(iteropt%nfi-1,comvl%nrotv).EQ.0
       ! >>>>Replica loop
       DO ipcurr=np_low,np_high
          ipx=ipcurr-np_low+1
          ! ANNEALING
          ! thermostats only if 1) standard dynamics 2) smaller MTS step
          if ( .not.cntl%use_mts .or. .not.mts_large_step ) then
             IF (pimd1%tstage.OR.pimd1%tpinm) THEN
                CALL anneal(vstage(:,:,:,ipcurr),c2,crge%n,scr)
                CALL berendsen(vstage(:,:,:,ipcurr),c2,crge%n,scr,0.0_real_8,0.0_real_8)
             ELSE
                CALL anneal(pivelp(:,:,:,ipcurr),c2,crge%n,scr)
                CALL berendsen(pivelp(:,:,:,ipcurr),c2,crge%n,scr,0.0_real_8,0.0_real_8)
             ENDIF
             ! ..UPDATE NOSE THERMOSTATS
             nosl%tmnose=pitmnose(MIN(ipcurr,2))
             IF (reset_gkt.AND.ipcurr==1) THEN
                scr(1:2)=gkt(1:2,1); gkt(1:2,1)=gkt1(1:2,1)
             ENDIF
             IF (pimd1%tstage.OR.pimd1%tpinm) THEN
                CALL noseup(vstage(:,:,:,ipcurr),c2,crge%n,ipcurr)
             ELSE
                CALL noseup(pivelp(:,:,:,ipcurr),c2,crge%n,ipcurr)
             ENDIF
             IF (reset_gkt.AND.ipcurr==1) gkt(1:2,1)=scr(1:2)
             ! FIRST HALF OF GLE EVOLUTION
             CALL gleback(ipcurr,.TRUE.,1)
             IF (pimd1%tstage.OR.pimd1%tpinm) THEN
                CALL gle_step(stagep(:,:,:,ipcurr),vstage(:,:,:,ipcurr),pma0s(1,ipcurr))
             ELSE
                CALL gle_step(pitaup(:,:,:,ipcurr),pivelp(:,:,:,ipcurr),rmass%pma)
             ENDIF
             CALL gleback(ipcurr,.TRUE.,0)
          endif
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
             ELSE
                CALL velupi(pivelp(:,:,:,ipcurr),pifion(:,:,:,ipcurr),1)
             ENDIF
          ENDIF
          ! ..UPDATE POSITIONS
          ener_com%ecnstr = 0.0_real_8
          ener_com%erestr = 0.0_real_8
          IF (paral%parent) THEN
             IF (pimd1%tstage.OR.pimd1%tpinm) THEN
                CALL dcopy(3*maxsys%nax*maxsys%nsx,stagep(1,1,1,ipcurr),1,tau0,1)
                CALL posupi(tau0,stagep(:,:,:,ipcurr),vstage(:,:,:,ipcurr))
                IF (cotc0%mcnstr.NE.0.AND.ipcurr.EQ.1)&
                     CALL cpmdshake(tau0,stagep(:,:,:,1),vstage(:,:,:,1))
             ELSE
                CALL dcopy(3*maxsys%nax*maxsys%nsx,pitaup(1,1,1,ipcurr),1,tau0,1)
                CALL posupi(tau0,pitaup(:,:,:,ipcurr),pivelp(:,:,:,ipcurr))
             ENDIF
          ENDIF
       ENDDO
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
          IF (corel%tinlc) CALL copot(rhoe,psi,ropt_mod%calste)
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
          ropt_mod%calste=cntl%tpres.AND.MOD(iteropt%nfi,cnti%npres).EQ.0
          if (cntl%use_mts) then
             call get_mts_forces(mts_large_step,.false.,n_inner_steps,n_large_steps,&
                mts_pure_dft,fion_high,c0_high,coldall_high(:,:,:,:,ipx),&
                nnowall_high(ipx),numcoldall_high(ipx),infi,&
                crge%n,c0(:,:,ipx:),c2,cm,sc0,cm(nx),vpp,eigv(:,ipx),&
                rhoe,psi,&
                pitaup(:,:,:,ipcurr),pivelp(:,:,:,ipcurr),taui(:,:,:,ipcurr),&
                pifion(:,:,:,ipcurr),&
                ifcalc,irec(:,ipcurr),.true.,.false.)
          else
             CALL forces_diag(crge%n,c0(:,:,ipx:),c2,cm,sc0,cm(nx),&
                  vpp,eigv(:,ipx),rhoe,psi,&
                  pitaup(:,:,:,ipcurr),pivelp(:,:,:,ipcurr),taui(:,:,:,ipcurr),&
                  pifion(:,:,:,ipcurr),&
                  ifcalc,irec(:,ipcurr),.TRUE.,.FALSE.)
          endif
          CALL dscal(3*maxsys%nax*maxsys%nsx,1._real_8/rnp,pifion(1,1,1,ipcurr),1)
          CALL dcopy(3*maxsys%nax*maxsys%nsx,pifion(1,1,1,ipcurr),1,&
               fionks(1,1,1,ipcurr),1)
          CALL fharm(pitaup,pifion(:,:,:,ipcurr),ipcurr,.TRUE.)
          IF (cntl%tdiag) CALL dcopy(nnx,rin0,1,rall(1,ipx),1)
          CALL dcopy(crge%n*nkpt%nkpts,crge%f,1,fall(1,ipx),1)
          CALL pi_ener(ipcurr,0)
          IF (ropt_mod%calste) THEN
             CALL totstr
             CALL dcopy(9,paiu(1,1),1,paiuall(1,1,ipcurr),1)
             CALL dscal(9,1._real_8/rnp,paiuall(1,1,ipcurr),1)
          ENDIF
       ENDDO
       CALL mp_sync(supergroup)
       ! <<<<End of replica loop
       IF (paral%parent) THEN
          CALL global(fionks,3*maxsys%nax*maxsys%nsx)
          IF (ropt_mod%calste) CALL global(paiuall,9)
          IF (pimd1%tpinm) THEN
             CALL getfnm(fionks,fstage,stagep,1)
             CALL taucl(fstage(:,:,:,1))
          ELSEIF (pimd1%tstage) THEN
             CALL getfu(fionks,fstage,stagep,1)
             CALL taucl(fstage(:,:,:,1))
          ENDIF
       ENDIF
       ! >>>>Replica loop
       DO ipcurr=np_low,np_high
          ipx=ipcurr-np_low+1
          ! ..FINAL UPDATE FOR VELOCITIES
          ener_com%ecnstr = 0.0_real_8
          ener_com%erestr = 0.0_real_8
          IF (paral%parent) THEN
             IF (pimd1%tpinm.OR.pimd1%tstage) THEN
                CALL s_velupi(vstage(:,:,:,ipcurr),fstage(:,:,:,ipcurr),1,ipcurr)
                IF (ipcurr.EQ.1) CALL rattle(stagep(:,:,:,1),vstage(:,:,:,1))
             ELSE
                CALL velupi(pivelp(:,:,:,ipcurr),pifion(:,:,:,ipcurr),1)
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
          
          ! thermostats only if 1) standard dynamics 2) smaller MTS step
          if ( .not.cntl%use_mts .or. .not.mts_large_step ) then
             ! SECOND HALF OF GLE EVOLUTION
             CALL gleback(ipcurr,.TRUE.,1)
             IF (pimd1%tstage.OR.pimd1%tpinm) THEN
                 CALL gle_step(stagep(:,:,:,ipcurr),vstage(:,:,:,ipcurr),pma0s(1,ipcurr))
             ELSE
                 CALL gle_step(pitaup(:,:,:,ipcurr),pivelp(:,:,:,ipcurr),rmass%pma)
             ENDIF
             CALL gleback(ipcurr,.TRUE.,0)
             ! ..UPDATE NOSE THERMOSTATS
             nosl%tmnose=pitmnose(MIN(ipcurr,2))
             IF (reset_gkt.AND.ipcurr==1) THEN
                scr(1:2)=gkt(1:2,1); gkt(1:2,1)=gkt1(1:2,1)
             ENDIF
             IF (pimd1%tstage.OR.pimd1%tpinm) THEN
                CALL noseup(vstage(:,:,:,ipcurr),c2,crge%n,ipcurr)
             ELSE
                CALL noseup(pivelp(:,:,:,ipcurr),c2,crge%n,ipcurr)
             ENDIF
             IF (reset_gkt.AND.ipcurr==1) gkt(1:2,1)=scr(1:2)
             IF (pimd1%tpinm.OR.pimd1%tstage) THEN
                CALL berendsen(vstage(:,:,:,ipcurr),c2,crge%n,scr,0.0_real_8,0.0_real_8)
                ! ANNEALING
                CALL anneal(vstage(:,:,:,ipcurr),c2,crge%n,scr)
             ELSE
                CALL berendsen(pivelp(:,:,:,ipcurr),c2,crge%n,scr,0.0_real_8,0.0_real_8)
                ! ANNEALING
                CALL anneal(pivelp(:,:,:,ipcurr),c2,crge%n,scr)
             ENDIF
          endif
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
             IF (reset_gkt.AND.ipcurr==1) THEN
                scr(1:2)=gkt(1:2,1); gkt(1:2,1)=gkt1(1:2,1)
             ENDIF
             IF (pimd1%tstage.OR.pimd1%tpinm) THEN
                CALL noseng(iteropt%nfi,vstage(:,:,:,ipcurr),dummy,enosp(ipcurr),&
                     dummy2,ipcurr)
             ELSE
                CALL noseng(iteropt%nfi,pivelp(:,:,:,ipcurr),dummy,enosp(ipcurr),&
                     dummy2,ipcurr)
             ENDIF
             IF (reset_gkt.AND.ipcurr==1) gkt(1:2,1)=scr(1:2)
             ! ..EFFECTIVELY CONSERVED AND HAMILTONIAN ENERGY
             ! ..NOTE: ECNSTR IS ALWAYS ZERO
             econs(ipcurr)=ekinp+etotv(ipcurr)/rnp+&
                  enosp(ipcurr)+&
                  eharv(ipcurr)+ener_com%ecnstr+ener_com%erestr+pi_egle(ipcurr)
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
       CALL global(pi_egle,1)
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
                  '    ENOSE(P)      EGLE    ECLASSIC     DIS'
             DO ip=1,pimd3%np_total
                IF (paral%io_parent)&
                     WRITE(6,&
                     '(I9,I5,F7.1,F10.6,F12.5,2X,F10.5,F10.5,F12.5, F8.2)')&
                     iteropt%nfi,ip,tempp(ip),eharv(ip),etotv(ip),&
                     enosp(ip),pi_egle(ip),econs(ip),disa(ip)
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
          time2=m_walltime()
          tcpu=(time2-time1)*0.001_real_8
          IF ((ropt_mod%engpri).AND.paral%io_parent)&
               WRITE(6,'(/,A,A)')&
               '      NFI   TEMP   EKINP(PRI)  EKINP(VIR)      EKS/P ',&
               '     EQUANT    ECLASSIC     TCPU'
          IF (paral%io_parent)&
               WRITE(6,'(I9,F7.1,5F12.5,F9.2)')&
               iteropt%nfi,tempa,qpkinp,qvkinp,etota,equant,econsa,tcpu
          IF (paral%io_parent)&
               WRITE(3,'(I9,F7.1,5F20.10,F9.2)')&
               iteropt%nfi,tempa,qpkinp,qvkinp,etota,equant,econsa,tcpu
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
       IF (ropt_mod%calste) THEN
          IF (pimd1%tstage.OR.pimd1%tpinm) THEN
             IF (cntl%tvirial) THEN
                pipaiu=pi_stress_vir(paiuall,stagep(:,:,:,1),vstage(:,:,:,1),pitaup,fionks,.TRUE.)
             ELSE
                pipaiu=pi_stress_nmm(paiuall,stagep,vstage,.TRUE.)
             ENDIF
          ELSE
             pipaiu=pi_stress_pri(paiuall,pitaup,pivelp,.TRUE.)
          ENDIF
          IF (grandparent) CALL wr_pi_stress(pipaiu)
       ENDIF
 
       IF (grandparent) THEN
          DO ipcurr=np_low,np_high
             ! ..Update accumulators
             i=ipcurr
             CALL paccd(ekinc(i),etotv(i),econs(i),eham(i),tempp(i),&
                  eharv(i),enose(i),enosp(i),disa(i),&
                  tcpu,qpkinp,qvkinp,ekina,etota,econsa,&
                  ehama,tempa,equant,accus(:,i),iteropt%nfi,1)
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
             CALL gleback(ipcurr,.FALSE.,1)
             nosl%tmnose=pitmnose(MIN(ipcurr,2))
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
       CALL paccd(d1,d1,d1,d1,d1,d1,d1,d1,d1,&
            d2,d2,d2,d2,d2,d2,d2,d2,d2,accus,iteropt%nfi,2)
    ENDIF
    IF (paral%parent) CALL prmem('   PI_DIAG')
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
    if (cntl%use_mts) then
       deallocate(fion_high,stat=ierr)
       if(ierr/=0) call stopgm(proceduren,'deallocation problem: fion_high',&
          __LINE__,__FILE__)
       if (mts_pure_dft) then
          deallocate(c0_high,stat=ierr)
          if(ierr/=0) call stopgm(proceduren,'deallocation problem: c0_high',&
             __LINE__,__FILE__)
          if (cntl%textrap) then
             DEALLOCATE(coldall_high,STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
                  __LINE__,__FILE__)
             DEALLOCATE(nnowall_high,STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
                  __LINE__,__FILE__)
             DEALLOCATE(numcoldall_high,STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
                  __LINE__,__FILE__)
          end if
       end if
    end if
    DEALLOCATE(pi_egle,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (np_local>1) THEN
       DEALLOCATE(pi_prng_com,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       IF (glepar%gle_mode>0) THEN
          DEALLOCATE(pi_glec,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
          DEALLOCATE(pi_gles,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
          DEALLOCATE(pi_glet,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
          DEALLOCATE(pi_glep,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
       ENDIF
    ENDIF
    IF (comvl%tsubrot) THEN
       DEALLOCATE(tauio,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    IF (cntl%tpres) THEN
       DEALLOCATE(paiuall,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
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
  END SUBROUTINE pi_diag
  ! ==================================================================


  ! Purpose: returns effective ionic forces in MTS algorithm
  !      The forces are either from the low level or from a 
  !      combination of high and low levels.
  !
  ! Author: Pablo Baudin
  ! Date: May 2018
  subroutine get_mts_forces(large_step,initialization,n_inner_steps, n_large_steps,&
        mts_pure_dft,fion_high,c0_high,cold_high,nnow_high,numcold_high,infi,&
        nstate,c0,c2,cr,csc0,cscr,vpp,eigv,&
        rhoe,psi,tau0,velp,taui,fion,ifcalc,irec,tfor,tinfo)

     implicit none

     logical, intent(in) :: large_step, initialization, mts_pure_dft
     integer, intent(inout) :: n_inner_steps
     integer, intent(in) :: n_large_steps
     logical :: tfor, tinfo
     integer :: nnow_high, numcold_high
     integer :: infi, nstate
     integer :: ifcalc, irec(:)
     real(real_8), allocatable :: fion_high(:,:,:)
     real(real_8) :: vpp(ncpw%ngw), eigv(*), rhoe(:,:)
     real(real_8) :: tau0(:,:,:), velp(:,:,:), taui(:,:,:), fion(:,:,:)
     complex(real_8), allocatable :: c0_high(:,:,:)
     complex(real_8) :: cold_high(:,:,:,:), c0(:,:,:), c2(nkpt%ngwk,nstate), cr(*)
     complex(real_8) :: csc0(nkpt%ngwk,nstate), cscr(*)
     complex(real_8) :: psi(:,:)

     character(*), parameter :: proceduren = 'get_mts_forces'
     character(len=100) :: title
     integer :: i, ia, is, x, N


     ! GET LOW LEVEL FORCES
     if (paral%io_parent) write(6,'(1x,3a)') 'MTS: LOW LEVEL FORCES: ', &
        mts%low_level, mts%low_dft_func
     select case(mts%low_level)
     case('EXTERNAL')
        CALL stopgm(procedureN,'LOW_LEVEL_FORCES EXTERNAL not available',&
           __LINE__,__FILE__)
     !  call get_external_forces('EXT_LOW_FORCES', tau0, fion)

     case('DFT')
        if (initialization .and. mts_pure_dft) then
           ! copy wf for extrapolation 
           call dcopy(2*size(c0,1)*size(c0,2)*size(c0,3),c0,1,c0_high,1)
        end if

        call set_mts_functional('LOW')

        call forces_diag(nstate,c0,c2,cr,csc0,cscr,vpp,eigv,&
           rhoe,psi,tau0,velp,taui,fion,ifcalc,irec,tfor,tinfo)

     case default
        if(paral%io_parent) print *, 'Low level forces not available with', mts%low_level
        call stopgm(proceduren,'wrong option for high level forces',&
           __LINE__,__FILE__)

     end select

     ! print low level forces to file
     if (mts%print_forces) then
        write(title,'(1x,a,i10,10x,a,i10)') 'STEP:',iteropt%nfi,'N_INNER_STEPS:',n_inner_steps
        call print_mts_forces(fion, title, 'LOW')
     end if

     if (large_step) then
        ! GET HIGH LEVEL FORCES
        if (paral%io_parent) write(6,'(1x,3a)') 'MTS: HIGH LEVEL FORCES: ', &
           mts%high_level, mts%high_dft_func
        select case(mts%high_level)
        case('EXTERNAL')
           CALL stopgm(procedureN,'HIGH_LEVEL_FORCES EXTERNAL not available',&
              __LINE__,__FILE__)
        !  call get_external_forces('EXT_HIGH_FORCES', tau0, fion_high)

        case('DFT')

           call set_mts_functional('HIGH')

           if (mts_pure_dft) then

              ! wf extrapolation
              if (.not.initialization .and. cntl%textrap) then
                 call extrapwf(infi,c0_high,cscr,cold_high,nnow_high,numcold_high,nstate,cnti%mextra)
              end if

              ! get forces
              call forces_diag(nstate,c0_high,c2,cr,csc0,cscr,vpp,eigv,&
                 rhoe,psi,tau0,velp,taui,fion_high,ifcalc,irec,tfor,tinfo)
           else

              ! In this case the wf extrap. is done in the main (outside) routine
              call forces_diag(nstate,c0,c2,cr,csc0,cscr,vpp,eigv,&
                 rhoe,psi,tau0,velp,taui,fion_high,ifcalc,irec,tfor,tinfo)
           end if

        case default
           if(paral%io_parent) print *, 'High level forces not available with', mts%high_level
           call stopgm(proceduren,'wrong option for high level forces',&
              __LINE__,__FILE__)

        end select

        ! print high level forces to file
        if (mts%print_forces) then
           write(title,'(1x,a,i10,10x,a,i10)') 'STEP:',iteropt%nfi,'N_LARGE_STEPS:',n_large_steps
           call print_mts_forces(fion_high, title, 'HIGH')
        end if

        ! get effective forces from difference between high and low level
        !
        !     F = F_low + (F_high - F_low) * N 
        !     F = F_high * N - F_low * (N - 1)
        !     
        !     where N is the MTS time-step factor
        ! 
        N = mts%timestep_factor
        do i=1,ions1%nat
           ia=iatpt(1,i)
           is=iatpt(2,i)
           do x=1,3
              fion(x,ia,is) = fion_high(x,ia,is) * N - fion(x,ia,is) * (N - 1)
           end do
        end do

        ! reinitialize inner counter
        n_inner_steps = 0

     end if

  end subroutine get_mts_forces

END MODULE pi_diag_utils
