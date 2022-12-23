MODULE prmdfile_utils
  USE andp,                            ONLY: rin0,&
                                             rmix,&
                                             rout0
  USE andr,                            ONLY: andr2
  USE anneal_utils,                    ONLY: anneal,&
                                             berendsen,&
                                             dampdyn,&
                                             tempramp
  USE atwf,                            ONLY: tmovr
  USE calc_alm_utils,                  ONLY: calc_alm,&
                                             give_scr_calc_alm
  USE cnst,                            ONLY: factem
  USE comvel_utils,                    ONLY: comvel
  USE comvelmod,                       ONLY: comvl
  USE coor,                            ONLY: fion,&
                                             tau0,&
                                             taup,&
                                             velp
  USE copot_utils,                     ONLY: copot,&
                                             give_scr_copot
  USE cotr,                            ONLY: cotc0
  USE ddipo_utils,                     ONLY: give_scr_ddipo
  USE detdof_utils,                    ONLY: detdof
  USE dispp_utils,                     ONLY: dispp
  USE do_perturbation_p_utils,         ONLY: do_perturbation
  USE dynit_utils,                     ONLY: dynit
  USE ekinpp_utils,                    ONLY: ekinpp
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: chrg,&
                                             ener_com
  USE error_handling,                  ONLY: stopgm
  USE extrap_utils,                    ONLY: extrap
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_app,&
                                             fo_old,&
                                             fo_verb
  USE finalp_utils,                    ONLY: finalp
  USE fint,                            ONLY: fint1
  USE forcep_utils,                    ONLY: rhoe_psi_size
  USE forces_diag_utils,               ONLY: forces_diag,&
                                             give_scr_forces_diag
  USE geofile_utils,                   ONLY: geofile
  USE gsize_utils,                     ONLY: gsize
  USE initrun_driver,                  ONLY: initrun
  USE initrun_utils,                   ONLY: give_scr_initrun
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE linres,                          ONLY: td01
  USE localize_utils,                  ONLY: localize2
  USE lr_tddft_utils,                  ONLY: give_scr_lr_tddft,&
                                             lr_tddft
  USE machine,                         ONLY: m_walltime
  USE metr,                            ONLY: metr_com
  USE mm_extrap,                       ONLY: cold, nnow, numcold
  USE moverho_utils,                   ONLY: give_scr_moverho,&
                                             moverho
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sync
  USE newcell_utils,                   ONLY: newcell
  USE nlcc,                            ONLY: corel
  USE norm,                            ONLY: gnmax,&
                                             gnorm
  USE nose,                            ONLY: glib,&
                                             ncdof
  USE parac,                           ONLY: parai,&
                                             paral
  USE phfac_utils,                     ONLY: phfac
  USE poin,                            ONLY: potr,&
                                             rhoo
  USE posupi_utils,                    ONLY: posuphf,&
                                             posupif
  USE prcp,                            ONLY: prcp_com
  USE printave_utils,                  ONLY: paccb
  USE printfor_utils,                  ONLY: printfor
  USE printp_utils,                    ONLY: printp
  USE proppt_utils,                    ONLY: give_scr_propcal,&
                                             propcal
  USE puttau_utils,                    ONLY: taucl
  USE rattle_utils,                    ONLY: rattle
  USE resetac_utils,                   ONLY: resetac
  USE rhopri_utils,                    ONLY: give_scr_rhopri,&
                                             rhopri
  USE rinvel_utils,                    ONLY: rinvel,&
                                             rvscal
  USE ropt,                            ONLY: infi,&
                                             iteropt,&
                                             ropt_mod
  USE rscvp_utils,                     ONLY: rscvp
  USE sample_utils,                    ONLY: sample_go,&
                                             sample_wait
  USE setirec_utils,                   ONLY: write_irec
  USE setsc_utils,                     ONLY: ihmat
  USE shake_utils,                     ONLY: cpmdshake
  USE soft,                            ONLY: soft_com
  USE spin,                            ONLY: clsd
  USE store_types,                     ONLY: cprint,&
                                             irec_ac,&
                                             irec_vel,&
                                             restart1,&
                                             rout1
  USE system,                          ONLY: &
       cnti, cntl, cntr, fpar, maxsys, nacc, ncpw, nkpt, parm, restf
  USE testex_utils,                    ONLY: testex
  USE teststore_utils,                 ONLY: teststore
  USE totstr_utils,                    ONLY: totstr
  USE vdwcmod,                         ONLY: vdwl,&
                                             vdwwfl
  USE velupi_utils,                    ONLY: velupc,&
                                             velupif
  USE wrener_utils,                    ONLY: wrprint_pr
  USE wv30_utils,                      ONLY: zhwwf
  USE zeroing_utils,                   ONLY: zeroing
#include "sizeof.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: prmdfile
  PUBLIC :: give_scr_prmdfile
  !public :: erback

CONTAINS

  ! ==================================================================
  SUBROUTINE prmdfile(c0,cm,c1,c2,sc0,vpp,gamx,gamy)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8) :: c0(nkpt%ngwk,crge%n,nkpt%nkpts), cm(:), c1(*), &
      c2(:,:), sc0(ncpw%ngw,*)
    REAL(real_8)                             :: vpp(*), gamx(*), gamy(*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'prmdfile'

    CHARACTER(len=100)                       :: filen
    CHARACTER(len=30)                        :: tag
    COMPLEX(real_8), ALLOCATABLE             :: psi(:,:)
    EXTERNAL                                 :: ddot
    INTEGER :: i, ierr, ifcalc, il_psi_1d, il_psi_2d, il_rhoe_1d, il_rhoe_2d, &
      irec(100), lenext, lscr, nmm, nnx, nstate, nx
    LOGICAL                                  :: ferror(2), qmmech
    REAL(real_8) :: ddot, disa, econs, eham, ekin1, ekin2, ekinc, ekincp, &
      ekinh, ekinh1, ekinh2, ekinp, enosc, enosp, mm_ekin, mm_temp, rmem, &
      tcpu, temp1, temp2, temph, tempp, time1, time2, vcmio(4)
    REAL(real_8), ALLOCATABLE :: eigv(:,:), rhoe(:,:), rinp(:), rm1(:), &
      save_fion(:,:,:), scr(:), taui(:,:,:), tauio(:,:), taur(:,:,:)

! Variables
! temporary alias for rhoo for allocation through memory90 
! QM/MM
! ==================================================================

    IF (cntl%tddft) CALL stopgm("MDFILE","TDDFT NOT YET SUPPORTED",&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    time1 =m_walltime()
    ALLOCATE(taup(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(fion(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(taup)!,SIZE(taup))
    CALL zeroing(fion)!,SIZE(fion))
    ALLOCATE(taui(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(taur(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(save_fion(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(tauio(3,maxsys%nax),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__) ! TODO check dimensions
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
    nstate=crge%n
    nacc = 22
    ncdof= 9
    iteropt%nfi  = 0
    ropt_mod%modens=.FALSE.
    ropt_mod%engpri=.FALSE.
    ropt_mod%calste=cntl%tpres
    ALLOCATE(eigv(crge%n,nkpt%nkpts),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! Extrapolation
    IF (cntl%textrap) THEN
       lenext=2*nkpt%ngwk*nstate*nkpt%nkpts*cnti%mextra
       rmem = 16._real_8*lenext*1.e-6_real_8
       ALLOCATE(cold(nkpt%ngwk,nstate,nkpt%nkpts,cnti%mextra),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       IF (paral%io_parent)&
            WRITE(6,'(A,T51,F8.3,A)') ' MDFILE| '&
            // 'EXTRAPOLATION WAVEFUNCTION HISTORY TAKES ',rmem,' MBYTES'
    ENDIF
    ! ==--------------------------------------------------------------==
    ! SCR ALLOCATION AND PARTITION (SCRATCH ARRAY).
    CALL rhoe_psi_size(il_psi_1d=il_psi_1d,il_psi_2d=il_psi_2d,&
         il_rhoe_1d=il_rhoe_1d,il_rhoe_2d=il_rhoe_2d)
    ALLOCATE(rhoe(il_rhoe_1d,il_rhoe_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    nmm=1
    IF (cntl%tddft) THEN
       ALLOCATE(rhoo(il_rhoe_1d,il_rhoe_2d),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)

       ALLOCATE(potr(il_rhoe_1d,il_rhoe_2d),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       IF (td01%ns_tri.GT.0) nmm=2
    ENDIF
    ALLOCATE(psi(il_psi_1d*nmm,il_psi_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL give_scr_prmdfile(lscr,tag)
    ALLOCATE(scr(lscr),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
99999 IF (cntl%tsampl) THEN
       CALL sample_wait
       IF (cnti%nomore<0) GOTO 10000
    ENDIF
    restf%nfnow=1
    ! ==--------------------------------------------------------------==
    ! TIME STEP FUNCTIONS
    CALL dynit(ekincp,ekin1,ekin2,temp1,temp2,ekinh1,ekinh2)
    ekinc=0.0_real_8
    ekinp=0.0_real_8
    ! Dont symmetrize density 
    cntl%tsymrho=.FALSE.
    ! ..Make sure TKFULL=.TRUE
    IF (tkpts%tkpnt.AND.(.NOT.tkpts%tkfull)) THEN
       IF (paral%io_parent) THEN
          WRITE(6,*)&
               ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          WRITE(6,*)&
               ' WARNING! USE KPOINTS FULL  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          WRITE(6,*)&
               ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == INITIALISATION                                               ==
    ! ==--------------------------------------------------------------==
    CALL initrun(irec,c0,c2,sc0,rhoe,psi,eigv)
    CALL mp_bcast(taup,SIZE(taup),parai%source,parai%allgrp)
    CALL dcopy(3*maxsys%nax*maxsys%nsx,taup(1,1,1),1,taui(1,1,1),1)
    ! ==--------------------------------------------------------------==
    IF (cprint%iprint_step.EQ.0) cprint%iprint_step=cnti%nomore+1
    ! ==--------------------------------------------------------------==
    ! INITIALIZE VELOCITIES
    IF (paral%parent) CALL detdof(tau0,taur)
    IF (irec(irec_vel)==0.AND..NOT.restart1%rgeo) THEN
       ener_com%ecnstr = 0.0_real_8
       ener_com%erestr = 0.0_real_8
       CALL rinvel(velp,c2,nstate)
       IF (paral%parent) CALL taucl(velp)
       IF (paral%parent) CALL rattle(tau0,velp)
       CALL rvscal(velp)
    ELSE
       IF (paral%parent) CALL taucl(velp)
       IF (paral%parent)  CALL rattle(tau0,velp)
       IF (cntl%trescale) CALL rvscal(velp)
    ENDIF
    IF (cntl%quenchp) CALL zeroing(velp)!,3*maxsys%nax*maxsys%nsx)
    ! COMPUTE THE IONIC TEMPERATURE TEMPP
    IF (paral%parent) THEN
       CALL ekinpp(ekinp,velp)
       tempp=ekinp*factem*2._real_8/glib
    ENDIF
    ! KINETIC ENERGY OF THE CELL
    IF (paral%parent) THEN
       ekinh=0.5_real_8*cntr%cmass*ddot(9,metr_com%htvel,1,metr_com%htvel,1)
    ENDIF
    ! RESET ACCUMULATORS
    IF (paral%parent.AND.irec(irec_ac).EQ.0) THEN
       CALL resetac(tau0,taui,iteropt%nfi)
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
    IF (cntl%tdiag) THEN
       IF (cntl%tlanc) nx=1
       IF (cntl%tdavi) nx=nkpt%ngwk*cnti%ndavv*nkpt%nkpnt+1
       IF (cntl%diis)  nx=((nkpt%ngwk*nstate+8)*cnti%mdiis*nkpt%nkpnt)/4
    ELSEIF (cntl%tsde) THEN
       nx=1
    ELSEIF (cntl%diis) THEN
       nx=(nkpt%ngwk*nstate+8)*cnti%mdiis/2+4
    ELSEIF (cntl%pcg) THEN
       nx=1
    ENDIF
    IF (paral%io_parent) THEN
       WRITE(6,'(1X,64("="))')
       WRITE(6,'(1X,"==",T25,A,T64,"==")') 'FORCES INITIALIZATION'
       WRITE(6,'(1X,64("="))')
    ENDIF
    IF (vdwl%vdwd.AND..NOT.vdwwfl%trwannc)&
         CALL localize2(TAU0,C0,C2,SC0,NSTATE)
    ifcalc=0
    CALL forces_diag(nstate,c0,c2,cm,sc0,cm(nx:),vpp,eigv,&
         rhoe,psi,&
         tau0,velp,taui,fion,ifcalc,&
         irec,.TRUE.,.TRUE.)
    IF (ropt_mod%calste) THEN
       CALL totstr
    ELSE
       CALL zeroing(metr_com%htfor)!,9)
    ENDIF
    ! 
    CALL dcopy(nnx,rin0,1,rm1(1),1)
    IF (paral%parent) THEN
       filen='ENERGIES'
       IF (paral%io_parent)&
            CALL fileopen(3,filen,fo_app+fo_verb,ferror(1))
    ENDIF
    IF (paral%io_parent) THEN
       WRITE(6,'(1X,64("="))')
       WRITE(6,'(1X,"==",T20,A,T64,"==")')'END OF FORCES INITIALIZATION'
       WRITE(6,'(1X,64("="),/)')
    ENDIF
    CALL write_irec(irec)
    ! ==--------------------------------------------------------------==
    ! == END INITIALIZATION                                           ==
    ! ==--------------------------------------------------------------==
    enosp=0.0_real_8
    enosc=0.0_real_8
    IF (teststore(0).AND.cntl%tsampl)&
         CALL zhwwf(2,irec,c0,c2,nstate,eigv,tau0,velp,taui,iteropt%nfi)
    IF (paral%parent) THEN
       ! MEAN SQUARE DISPLACEMENT OF DIFFERENT IONIC SPECIES
       IF (paral%parent) CALL dispp(tau0,taui,disa)
       econs=ekinp+ener_com%etot+enosp+enosc+ener_com%ecnstr+ener_com%erestr+prcp_com%druck*parm%omega
       eham=econs+ekinc+ekinh
       time2 =m_walltime()
       tcpu = (time2 - time1)*0.001_real_8
       IF (paral%io_parent) THEN
          WRITE(6,'(A,T50,F8.2,A8)') ' TIME FOR INITIALIZATION:',tcpu,' SECONDS'
          WRITE(6,'(//,1X,64("="))')
          WRITE(6,'(1X,"=",T20,A,T65,"=")')'MOLECULAR DYNAMICS SIMULATION'
          WRITE(6,'(1X,64("="))')
       ENDIF
    ENDIF
    ! OPEN PREVIOUSLY SAVED TRAJECTORIES
    IF (paral%io_parent) THEN
       IF (rout1%xtin) THEN
          CALL fileopen(64,'TRAJSAVED.xyz',fo_old,ferror(1))
          IF (ferror(1)) CALL stopgm('PRMDFILE',' TRAJSAVED.xyz NOT FOUND',&
               __LINE__,__FILE__)
       ELSE
          CALL fileopen(64,'TRAJSAVED',fo_old,ferror(1))
          IF (ferror(1)) CALL stopgm('PRMDFILE',' TRAJSAVED NOT FOUND',&
               __LINE__,__FILE__)
       ENDIF
       CALL fileopen(65,'CELLSAVED',fo_old,ferror(2))
       IF (ferror(2)) CALL stopgm('PRMDFILE',' CELLSAVED NOT FOUND',&
            __LINE__,__FILE__)
    ENDIF
    ! UPDATE NSKIP
    cnti%nskip=cnti%nskip+iteropt%nfi*cnti%nsample
    IF (iteropt%nfi>0.AND.paral%io_parent) THEN
       WRITE(6,'(/,A,T56,I10,/)') ' RESTART WITH NSKIP: ',cnti%nskip
    ENDIF
    ! ==================================================================
    ! ==          THE BASIC LOOP FOR MOLECULAR DYNAMICS               ==
    ! ==                 USING VELOCITY VERLET                        ==
    ! ==================================================================
    DO infi=1,cnti%nomore
       CALL mp_sync(parai%allgrp)
       time1=m_walltime()
       iteropt%nfi=iteropt%nfi+1
       comvl%subcom=comvl%tsubcom.AND.MOD(iteropt%nfi-1,comvl%ncomv)==0
       ! ANNEALING
       CALL anneal(velp,c2,nstate,metr_com%htvel)
       CALL berendsen(velp,c2,nstate,metr_com%htvel,0.0_real_8,ekinh)
       ! SUBTRACT CENTER OF MASS VELOCITY
       IF (paral%parent.AND.comvl%subcom) CALL comvel(velp,vcmio,.TRUE.)
       ! UPDATE VELOCITIES
       IF (paral%parent) CALL velupif(velp,fion,1)
#if defined (__QMECHCOUPL)
       IF (paral%parent .AND. qmmech) THEN
          CALL mm_cpmd_velup(cntr%delt_ions)
       ENDIF
#endif
       ! UPDATE POSITIONS
       ener_com%ecnstr = 0.0_real_8
       ener_com%erestr = 0.0_real_8

       CALL posupif(64,tau0,taup,velp,ferror(1))
       CALL posuphf(65,metr_com%ht,metr_com%htvel,ferror(2))
       IF (paral%io_parent) THEN
          WRITE(6,'(1X,64("*"))')
          IF (rout1%xtin) THEN
             WRITE(6,'(1X,"*",2X,A,I10,A,7X,"*")') 'CONFIGURATION ',&
                  cnti%nskip+infi*cnti%nsample,' READ FROM FILE TRAJSAVED.xyz'
          ELSE
             WRITE(6,'(1X,"*",2X,A,I10,A,11X,"*")') 'CONFIGURATION ',&
                  cnti%nskip+infi*cnti%nsample,' READ FROM FILE TRAJSAVED'
          ENDIF
          WRITE(6,'(1X,64("*"))')
       ENDIF
       IF (paral%parent) THEN
          IF (cotc0%mcnstr/=0) CALL cpmdshake(tau0,taup,velp)
#if defined (__QMECHCOUPL)
          IF (qmmech) THEN
             CALL mm_cpmd_update_links(taup, ions0%na, ions1%nsp, maxsys%nax, ions1%nat)
             CALL mm_cpmd_posup(cntr%delt_ions)
          ENDIF
#endif
       ENDIF
       CALL mp_bcast(ferror,SIZE(ferror),parai%source,parai%allgrp)
       IF (ferror(1).OR.ferror(2)) THEN
          IF (paral%io_parent) THEN
             WRITE(6,'(1X,64("!"))')
             WRITE(6,'(" !!",A,T64,"!!")')&
                  ' PRMDFILE| ENCOUNTERED END OF FILE '
             WRITE(6,'(" !!",A,I6,T64,"!!")')&
                  '         WHILE READING CONFIGURATION ',cnti%nskip+infi*cnti%nsample
             WRITE(6,'(1X,64("!"))')
          ENDIF
          GOTO 100
       ENDIF
       CALL mp_bcast(taup,SIZE(taup),parai%source,parai%allgrp)
       CALL mp_bcast(metr_com%ht,SIZE(metr_com%ht),parai%source,parai%allgrp)
       DO i=1,3
          parm%a1(i) = metr_com%ht(1,i)
          parm%a2(i) = metr_com%ht(2,i)
          parm%a3(i) = metr_com%ht(3,i)
       ENDDO
       CALL ihmat(metr_com%ht,metr_com%htm1,parm%omega)
       CALL newcell
       ! Reset swap file to update (option BLOCK and ALL)
       IF (tkpts%tkblock) CALL reskpt_swap
       CALL phfac(taup)
       IF (corel%tinlc) CALL copot(rhoe,psi,ropt_mod%calste)
       IF (tmovr) THEN
          CALL dcopy(nnx,rin0,1,rhoe,1)
          CALL moverho(rhoe,psi)
          CALL dcopy(nnx,rhoe,1,rin0,1)
       ENDIF
       IF (fint1%ttrot) THEN
          CALL calc_alm
       ENDIF
       ! CALCULATE THE FORCES
       IF (cntl%textrap) THEN
          ! Extrapolate wavefunctions
          CALL extrapwf(infi,c0,scr,cold,nnow,numcold,nstate,cnti%mextra)
       ENDIF
       IF (cntl%tlanc) nx=1
       IF (cntl%tdavi) nx=cnti%ndavv*nkpt%nkpnt+1
       ! RESPONSE calculation
       CALL forces_diag(nstate,c0,c2,cm,sc0,cm(nx:),vpp,eigv,&
            rhoe,psi,&
            taup,velp,taui,fion,ifcalc,&
            irec,.TRUE.,.FALSE.)
       IF (cntl%tddft) THEN
          CALL lr_tddft(c0(:,:,1),c1,c2,sc0,rhoe,psi,taup,fion,eigv,&
               nstate,.FALSE.,td01%ioutput)
       ENDIF
       IF (cntl%tresponse) THEN
          CALL erback(0)
          CALL dcopy(3*maxsys%nax*maxsys%nsx,fion,1,save_fion,1)
          CALL dcopy(3*maxsys%nax*maxsys%nsx,tau0(1,1,1),1,taur(1,1,1),1)
          CALL dcopy(3*maxsys%nax*maxsys%nsx,taup(1,1,1),1,tau0(1,1,1),1)
          CALL do_perturbation(c0,c2,nstate)
          CALL erback(1)
          CALL dcopy(3*maxsys%nax*maxsys%nsx,save_fion,1,fion,1)
          CALL dcopy(3*maxsys%nax*maxsys%nsx,taur(1,1,1),1,tau0(1,1,1),1)
       ENDIF
       IF (ropt_mod%calste) THEN
          CALL totstr
       ELSE
          CALL zeroing(metr_com%htfor)!,9)
       ENDIF
#if defined (__QMECHCOUPL)
       ! ----------------------------------------------
       ! QMMM electrostatic coupling
       ! ----------------------------------------------
       IF (qmmech) THEN
          CALL mm_cpmd_ext_pot_f77
          CALL mm_cpmd_elstat(taup,ions0%na,ions1%nsp,maxsys%nax,ions1%nat,c0,scr)
       ENDIF
#endif
       ! ==================================================================
       ! Damped Dynamics
       CALL dampdyn(velp,fion,cm,c2,nstate,metr_com%htvel,metr_com%htfor)
       ! ==================================================================
       ! FINAL UPDATE FOR VELOCITIES
       ener_com%ecnstr = 0.0_real_8
       ener_com%erestr = 0.0_real_8
       IF (paral%parent) THEN
          CALL velupif(velp,fion,1)
          CALL velupc(1)
#if defined (__QMECHCOUPL)
          IF (qmmech) THEN
             CALL mm_cpmd_velup(cntr%delt_ions)
          ENDIF
#endif
          CALL rattle(taup,velp)
       ENDIF
       IF (paral%parent) CALL geofile(taup,velp,'WRITE')
       ! IONIC TEMPERATURE CONTROL
       IF (paral%parent) THEN
          ! KINETIC ENERGY OF THE CELL
          ekinh=0.5_real_8*cntr%cmass*ddot(9,metr_com%htvel,1,metr_com%htvel,1)
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
          CALL rscvp(temp1,temp2,tempp,velp)
       ENDIF
#if defined (__QMECHCOUPL)
       ! ---------------------------------------------------------------
       ! MM Temperature control - annealing schedule not yet implemented
       ! ---------------------------------------------------------------
       IF (paral%parent .AND. qmmech) THEN
          CALL mm_cpmd_temp_control(temp1,temp2,mm_temp,cntl%tcp)
          CALL mm_cpmd_ekin(mm_ekin,mm_temp)
       ENDIF
#endif
       ! SUBTRACT CENTER OF MASS VELOCITY
       IF (paral%parent.AND.comvl%subcom) CALL comvel(velp,vcmio,.FALSE.)
       ! ANNEALING
       CALL berendsen(velp,c2,nstate,metr_com%htvel,0.0_real_8,ekinh)
       CALL anneal(velp,c2,nstate,metr_com%htvel)
       IF (paral%parent) THEN
          ! ..CELL TEMPERATURE CONTROL
          CALL rscvc(ekinh,metr_com%htvel,ekinh1,ekinh2)
          CALL ekinpp(ekinp,velp)
          tempp=ekinp*factem*2._real_8/glib
          ! KINETIC ENERGY OF THE CELL
          ekinh=0.5_real_8*cntr%cmass*ddot(9,metr_com%htvel,1,metr_com%htvel,1)
          ! TEMPERATURE OF THE CELL
          temph=2.0_real_8*ekinh*factem/REAL(ncdof,kind=real_8)
       ENDIF
       ! MEAN SQUARE DISPLACEMENT OF DIFFERENT IONIC SPECIES
       IF (paral%parent) CALL dispp(taup,taui,disa)
       ! CALCULATE PROPERTIES DURING SIMULATION.
       cntl%caldip=cntl%tdipd
       vdwwfl%twannup=vdwwfl%twannup.OR.(infi==1.AND..NOT.vdwwfl%trwannc)
       CALL propcal(c0,c2,cm,sc0,taup,eigv,crge%f,ener_com%amu,&
            rhoe,psi,nstate,nkpt%nkpnt,iteropt%nfi,infi)
       ! PRINTOUT the evolution of the accumulators every time step
       IF (paral%parent) THEN
          econs=ekinp+ener_com%etot+enosp+enosc+ener_com%ecnstr+ener_com%erestr+prcp_com%druck*parm%omega
          eham=econs+ekinc+ekinh
#if defined (__QMECHCOUPL)
          IF (qmmech) THEN
             econs=econs + mm_ekin
             IF (paral%io_parent)&
                  WRITE (6,'(50x, "TEMP: ",f10.0)') mm_temp
          ENDIF
#endif
          time2=m_walltime()
          tcpu=(time2-time1)*0.001_real_8
          ! PRINTOUT the evolution of the accumulators every time step
          CALL wrprint_pr(0.0_real_8,ekinh,tempp,ener_com%etot,econs,eham,&
               disa,taup,fion,tcpu,iteropt%nfi,infi)
          ! UPDATE ACCUMULATORS
          CALL paccb(0.0_real_8,tempp,ener_com%etot,econs,eham,disa,tcpu,temph,&
               parm%omega,0.0_real_8,enosp,iteropt%nfi,1)
          ! STORE IONIC COORDINATES AND VELOCITIES FOR STATISTICS
          ropt_mod%movie=rout1%mout .AND. MOD(iteropt%nfi-1,cnti%imovie).EQ.0
          ropt_mod%rprint=rout1%rout .AND. MOD(iteropt%nfi-1,cnti%ntraj).EQ.0
          ropt_mod%txyz=rout1%xtout .AND. MOD(iteropt%nfi-1,cnti%ntraj).EQ.0
          ropt_mod%tdcd=rout1%dcout .AND. MOD(iteropt%nfi-1,cnti%ntraj).EQ.0
          CALL printp(taur,taup,velp)
          CALL printfor (taup,fion)
#if defined (__QMECHCOUPL)
          ! --------------------------------------------------------
          ! ..       Write to QMMM trajectory.  This trajectory file contains
          ! the atoms that make up the "real" system. (no dummies)
          ! --------------------------------------------------------
          IF (qmmech) THEN
             IF (rout1%rout .AND. MOD(iteropt%nfi-1,cnti%ntraj)==0) THEN
                CALL mm_cpmd_write_trajectory(iteropt%nfi)
             ENDIF
          ENDIF
#endif
          time2=m_walltime()
          tcpu=(time2-time1)*0.001_real_8
       ENDIF
       IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)
       IF (infi==cnti%nomore) soft_com%exsoft=.TRUE.
       ! periodic output of density/wavefunction etc.
       IF (rout1%rhoout.AND.(rout1%nrhoout>0)) THEN
          IF (MOD(iteropt%nfi-1,rout1%nrhoout)==0) THEN
             CALL rhopri(c0,tau0,rhoe,psi(:,1),nstate,nkpt%nkpnt)
          ENDIF
       ENDIF
       IF (teststore(iteropt%nfi).OR.soft_com%exsoft)&
            CALL zhwwf(2,irec,c0,c2,nstate,eigv,taup,velp,taui,iteropt%nfi)
#if defined (__QMECHCOUPL)
       ! -----------------------
       ! Write QMMM restart file
       ! -----------------------
       IF (paral%parent.AND.qmmech) THEN
          IF (MOD(iteropt%nfi,store1%istore)==0.OR.infi==cnti%nomore.OR.soft_com%exsoft) THEN
             CALL mm_cpmd_write_restart
          ENDIF
       ENDIF
#endif
       ! temperature ramping
       CALL tempramp(temp1,temp2)
       ! STOP THE RUN IF THE USER HAS SET THE SIGNAL 30
       IF (soft_com%exsoft) GOTO 100
       ! UPDATE IONIC POSITIONS
       CALL dcopy(3*maxsys%nax*maxsys%nsx,taup(1,1,1),1,tau0(1,1,1),1)
       ! UPDATE DENSITY
       CALL extrap(nnx,andr2%alxmix,rm1,rin0,rinp)
       CALL dcopy(nnx,rin0,1,rm1(1),1)
       CALL dcopy(nnx,rinp(1),1,rin0,1)
       ! ==================================================================
       ! ==     END OF MAIN LOOP                                         ==
       ! ==================================================================
    ENDDO
    IF (paral%io_parent) THEN
       WRITE(6,'(1X,64("="))')
       WRITE(6,'(1X,"=",T17,A,T65,"=")')'END OF MOLECULAR DYNAMICS SIMULATION'
       WRITE(6,'(1X,64("="),/,/)')
    ENDIF
    ! ==--------------------------------------------------------------==
100 CONTINUE
    IF (rout1%rhoout.AND.(rout1%nrhoout<=0))&
         CALL rhopri(c0,tau0,rhoe,psi(:,1),nstate,nkpt%nkpnt)
    ! PRINT ACCUMULATORS.
    IF (paral%parent) CALL paccb(0.0_real_8,tempp,ener_com%etot,econs,eham,disa,tcpu,temph&
         ,parm%omega,0.0_real_8,enosp,iteropt%nfi,2)
    IF (paral%parent) CALL gsize(fion,gnmax,gnorm)
    IF (paral%parent) CALL finalp(tau0,fion,velp,eigv)
    IF (cntl%tsampl) THEN
       CALL sample_go
       GOTO 99999
    ENDIF
10000 CONTINUE
    ! ==--------------------------------------------------------------==
    DEALLOCATE(scr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(psi,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rhoe,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(eigv,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
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
    DEALLOCATE(tauio,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(save_fion,STAT=ierr)
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
    IF (paral%io_parent) CALL fileclose(64)
    IF (paral%io_parent) CALL fileclose(65)
    IF (cntl%textrap) THEN
       DEALLOCATE(cold,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    IF (cntl%tddft) THEN
       DEALLOCATE(potr,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(rhoo,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE prmdfile
  ! ==================================================================
  SUBROUTINE give_scr_prmdfile(lprmdfile,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lprmdfile
    CHARACTER(len=30)                        :: tag

    INTEGER :: lcalc_alm, lcopot, lddipo, lforces_diag, linitrun, lmoverho, &
      lpropcal, lrhopri, ltddft, nstate

    nstate=crge%n
    lcalc_alm=0
    lcopot=0
    lrhopri=0
    lmoverho=0
    ltddft=0
    lddipo=0
    CALL give_scr_initrun(linitrun,tag)
    IF (fint1%ttrot) CALL give_scr_calc_alm(lcalc_alm,tag)
    CALL give_scr_forces_diag(lforces_diag,tag,nstate,.TRUE.)
    IF (corel%tinlc) CALL give_scr_copot(lcopot,tag)
    IF (rout1%rhoout) CALL give_scr_rhopri(lrhopri,tag,nstate)
    CALL give_scr_propcal(lpropcal,tag,nstate)
    IF (tmovr) CALL give_scr_moverho(lmoverho,tag)
    IF (cntl%tddft) CALL give_scr_lr_tddft(ltddft,.TRUE.,tag)
    IF (vdwl%vdwd) CALL give_scr_ddipo(lddipo,tag)
    lprmdfile=MAX(linitrun,lcalc_alm,lforces_diag,lcopot,lrhopri,lpropcal,&
         lmoverho,lddipo)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_prmdfile
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
  SUBROUTINE erback(mode)
    ! ==--------------------------------------------------------------==
    ! == Save or restore values of energies                           ==
    ! == Mode 0 -> SAVE                                               ==
    ! ==      1 -> RESTORE                                            ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: mode

    INTEGER, SAVE                            :: bnfi, bnomore
    REAL(real_8), SAVE :: bcsumg, bcsumr, becnstr, begc, behee, behep, behii, &
      beht, bekin, benl, bepseu, berestr, beself, besr, betot, bexc, bvxc

    IF (mode.EQ.0) THEN
       betot=ener_com%etot
       bekin=ener_com%ekin
       bepseu=ener_com%epseu
       benl=ener_com%enl
       beht=ener_com%eht
       behep=ener_com%ehep
       behee=ener_com%ehee
       behii=ener_com%ehii
       bexc=ener_com%exc
       bvxc=ener_com%vxc
       beself=ener_com%eself
       besr=ener_com%esr
       begc=ener_com%egc
       becnstr=ener_com%ecnstr
       berestr=ener_com%erestr
       bcsumg=chrg%csumg
       bcsumr=chrg%csumr
       bnfi=iteropt%nfi
       bnomore=cnti%nomore
       cnti%nomore=cnti%nomore_iter
       ! BECONS=ECONS
       ! BENOSE=ENOSE
       ! BENOSP=ENOSP
       ! BECNSTR=ECNSTR
       ! BEHAM=EHAM
       ! BEKINC=EKINC
    ELSE
       ener_com%etot=betot
       ener_com%ekin=bekin
       ener_com%epseu=bepseu
       ener_com%enl=benl
       ener_com%eht=beht
       ener_com%ehep=behep
       ener_com%ehee=behee
       ener_com%ehii=behii
       ener_com%exc=bexc
       ener_com%vxc=bvxc
       ener_com%eself=beself
       ener_com%esr=besr
       ener_com%egc=begc
       ener_com%ecnstr=becnstr
       ener_com%erestr=berestr
       chrg%csumg=bcsumg
       chrg%csumr=bcsumr
       iteropt%nfi=bnfi
       cnti%nomore=bnomore
       ! ECONS=BECONS
       ! ENOSE=BENOSE
       ! ENOSP=BENOSP
       ! ECNSTR=BECNSTR
       ! EHAM=BEHAM
       ! EKINC=BEKINC
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE erback
  ! ==================================================================

END MODULE prmdfile_utils
