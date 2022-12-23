MODULE npt_md_utils
  USE andp,                            ONLY: rin0,&
                                             rmix,&
                                             rout0
  USE andr,                            ONLY: andr2
  USE anneal_utils,                    ONLY: anneal,&
                                             berendsen,&
                                             dampdyn,&
                                             tempramp
  USE cnst,                            ONLY: factem
  USE comvel_utils,                    ONLY: comvel
  USE comvelmod,                       ONLY: comvl
  USE coor,                            ONLY: fion,&
                                             tau0,&
                                             taup,&
                                             velp
  USE copot_utils,                     ONLY: copot,&
                                             give_scr_copot
  USE csize_utils,                     ONLY: csize
  USE ddipo_utils,                     ONLY: give_scr_ddipo
  USE deort_utils,                     ONLY: give_scr_deort
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
                                             fo_mark,&
                                             fo_verb
  USE finalp_utils,                    ONLY: finalp
  USE forcedr_driver,                  ONLY: forcedr
  USE forcedr_utils,                   ONLY: give_scr_forcedr
  USE forcep_utils,                    ONLY: rhoe_psi_size
  USE forces_diag_utils,               ONLY: forces_diag
  USE freqs_utils,                     ONLY: freqs
  USE geofile_utils,                   ONLY: geofile
  USE geq0mod,                         ONLY: geq0
  USE gsize_utils,                     ONLY: gsize
  USE initrun_driver,                  ONLY: initrun
  USE initrun_utils,                   ONLY: give_scr_initrun
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE localize_utils,                  ONLY: localize2
  USE machine,                         ONLY: m_walltime
  USE metr,                            ONLY: eps,&
                                             metr_com,&
                                             veps
  USE mm_extrap,                       ONLY: cold, nnow, numcold
  USE mp_interface,                    ONLY: mp_bcast
  USE newcell_utils,                   ONLY: give_scr_newcell,&
                                             newcell
  USE nlcc,                            ONLY: corel
  USE norm,                            ONLY: cnorm,&
                                             gemax,&
                                             gnmax,&
                                             gnorm
  USE noscinit_utils,                  ONLY: noscinit
  USE nose,                            ONLY: glib,&
                                             ncdof
  USE noseinit_utils,                  ONLY: noseinit
  USE noseng_utils,                    ONLY: noseng
  USE nosepa_utils,                    ONLY: nosepa
  USE noseup_utils,                    ONLY: noseup
  USE nospinit_utils,                  ONLY: nospinit
  USE ortho_utils,                     ONLY: give_scr_ortho
  USE parac,                           ONLY: parai,&
                                             paral
  USE phfac_utils,                     ONLY: phfac
  USE poin,                            ONLY: rhoo
  USE posupa_utils,                    ONLY: give_scr_posupa,&
                                             posupa
  USE posupi_utils,                    ONLY: posupih,&
                                             posupih_iso,&
                                             posupihshock
  USE prcp,                            ONLY: prcp_com,&
                                             prcpl
  USE printave_utils,                  ONLY: paccb
  USE printp_utils,                    ONLY: printp
  USE proja_utils,                     ONLY: proja
  USE pslo,                            ONLY: pslo_com
  USE puttau_utils,                    ONLY: taucl
  USE quenbo_utils,                    ONLY: give_scr_quenbo,&
                                             quenbo
  USE rattle_utils,                    ONLY: rattle
  USE rekine_utils,                    ONLY: rekine
  USE resetac_utils,                   ONLY: resetac
  USE rhopri_utils,                    ONLY: give_scr_rhopri,&
                                             rhopri
  USE rinvel_utils,                    ONLY: rinvel,&
                                             rvscal
  USE ropt,                            ONLY: infi,&
                                             iteropt,&
                                             ropt_mod
  USE rortv_utils,                     ONLY: give_scr_rortv,&
                                             rortv
  USE rscve_utils,                     ONLY: rscve
  USE rscvp_utils,                     ONLY: rscvp
  USE setirec_utils,                   ONLY: write_irec
  USE setsc_utils,                     ONLY: ihmat
  USE shock,                           ONLY: shock1
  USE soft,                            ONLY: soft_com
  USE spin,                            ONLY: clsd
  USE store_types,                     ONLY: &
       cprint, irec_ac, irec_noc, irec_noe, irec_nop1, irec_nop2, irec_nop3, &
       irec_vel, restart1, rout1
  USE system,                          ONLY: &
       cnti, cntl, cntr, fpar, maxsys, nacc, ncpw, nkpt, parm, restf
  USE testex_utils,                    ONLY: testex
  USE teststore_utils,                 ONLY: teststore
  USE totstr_utils,                    ONLY: totstr
  USE utils,                           ONLY: zclean,&
                                             zclean_k
  USE vdwcmod,                         ONLY: vdwl,&
                                             vdwwfl
  USE velupa_utils,                    ONLY: velupa
  USE velupi_utils,                    ONLY: velupc,&
                                             velupi,&
                                             velupidamp,&
                                             velupishock
  USE vepsup_utils,                    ONLY: vepsup,&
                                             vepsupdamp
  USE wrener_utils,                    ONLY: wrprint_pr
  USE wrgeo_utils,                     ONLY: wrgeof
  USE wv30_utils,                      ONLY: zhwwf
  USE zeroing_utils,                   ONLY: zeroing
#include "sizeof.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: npt_cpmd
  PUBLIC :: npt_bomd

CONTAINS


  ! ==================================================================
  SUBROUTINE npt_cpmd(c0,cm,c2,sc0,gamx,gamy)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(ncpw%ngw,crge%n,1), cm(:), &
                                                c2(:,:), sc0(ncpw%ngw,crge%n)
    REAL(real_8)                             :: gamx(*), gamy(*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'npt_cpmd'

    CHARACTER(len=100)                       :: filen
    CHARACTER(len=30)                        :: tag
    COMPLEX(real_8), ALLOCATABLE             :: psi(:,:)
    INTEGER                                  :: i, ierr, il_psi_1d, &
                                                il_psi_2d, il_rhoe_1d, &
                                                il_rhoe_2d, irec(100), itemp, &
                                                j, lscr, nstate
    LOGICAL                                  :: ferror
    REAL(real_8) :: disa, econs, eham, ekin1, ekin2, ekinc, ekincp, ekinh, &
      ekinh1, ekinh2, ekinp, enosc, enose, enosp, ff, tcpu, temp1, temp2, &
      temph, tempp, time1, time2, vcmio(4)
    REAL(real_8), ALLOCATABLE                :: eigm(:), eigv(:), rhoe(:,:), &
                                                scr(:), taui(:,:,:), &
                                                taur(:,:,:)
    REAL(real_8), EXTERNAL                   :: ddot

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
    CALL zeroing(taup)!,SIZE(taup))
    CALL zeroing(fion)!,SIZE(fion))
    CALL zeroing(taui)!,SIZE(taui))
    CALL zeroing(taur)!,SIZE(taur))
    nacc = 22
    nstate=crge%n
    ropt_mod%modens=.FALSE.
    ropt_mod%engpri=.FALSE.
    ropt_mod%calste=.TRUE.
    ! ==--------------------------------------------------------------==
    ! 06/2008 AK
    ! NPT FLEXIBLE CELL DOES NOT WORK WITH VANDERBILT PSEUDOPOTENTIALS.
    ! FOR STARTERS THE STRESS TENSOR IS BROKEN.
    ! IT WORKS IN QUANTUM ESPRESSO...
    ! 08/2017 MB-TI - The above statement is wrong. Fixed on 28/08/2017.
  ! IF (pslo_com%tivan) THEN
  !    CALL stopgm('NPT_CPMD','NPT MD IS BROKEN FOR VANDERBILT USPPS',& 
  !         __LINE__,__FILE__)
  ! ENDIF
    ! ==--------------------------------------------------------------==
    IF (pslo_com%tivan.OR.cntl%nonort) THEN
       ALLOCATE(eigv(nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(eigm(nstate*nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ELSE
       ! avoid fortran runtime error 'not allocated'
       ALLOCATE(eigv(1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(eigm(1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
    ENDIF
    IF (cntl%tharm.AND.paral%parent) THEN
       ff=crge%f(1,1)
       DO i=1,nstate
          IF (ff.NE.crge%f(i,1)) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(A,A)') ' HARMONIC REFERENCE SYSTEM INTEGRATION',&
                  ' ONLY POSSIBLE WITH EQUAL OCCUPATION NUMBERS'
             CALL stopgm('NPT_CPMD',' ',& 
                  __LINE__,__FILE__)
          ENDIF
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == INITIALIZATION                                               ==
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
    CALL give_scr_npt_md(lscr,tag)
    ALLOCATE(scr(lscr),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    restf%nfnow=1
    ! ==--------------------------------------------------------------==
    ! TIME STEP FUNCTIONS
    CALL dynit(ekincp,ekin1,ekin2,temp1,temp2,ekinh1,ekinh2)
    ekinc=0.0_real_8
    ekinp=0.0_real_8
    ! PARAMETERS FOR THE NOSE-HOOVER THERMOSTATS
    CALL initrun(irec,c0,cm,sc0,rhoe,psi,eigv)
    CALL mp_bcast(taup,SIZE(taup),parai%io_source,parai%cp_grp)
    CALL dcopy(3*maxsys%nax*maxsys%nsx,taup(1,1,1),1,taui(1,1,1),1)
    ! INITIALIZE WF CENTERS & SPREAD
    IF (vdwl%vdwd.AND..NOT.vdwwfl%trwannc) THEN
       CALL localize2(tau0,c0,c2,sc0,nstate)
    ENDIF
    ncdof = 1
    IF ((cntl%tnosep.OR.cntl%tnosee.OR.cntl%tnosec).AND.paral%parent) CALL nosepa(1,1)
    ! INITIALIZE VELOCITIES
    IF (paral%parent) CALL detdof(tau0,taur)
    IF (irec(irec_vel).EQ.0.AND..NOT.restart1%rgeo) THEN
       ener_com%ecnstr = 0.0_real_8
       ener_com%erestr = 0.0_real_8
       CALL rinvel(velp,cm,nstate)
       IF (paral%parent) CALL taucl(velp)
       IF (paral%parent) CALL rattle(tau0,velp)
       CALL rvscal(velp)
    ELSE
       IF (paral%parent) CALL taucl(velp)
       IF (paral%parent) CALL rattle(tau0,velp)
       IF (cntl%trescale) CALL rvscal(velp)
    ENDIF
    IF (cntl%quenchp) THEN
       CALL zeroing(velp)!,3*maxsys%nax*maxsys%nsx)
    ENDIF
    IF (cntl%quenche) CALL zeroing(cm)!,ngw*nstate)
    IF (cntl%quenchc) CALL zeroing(metr_com%htvel)!,9)
    ! >>>
    ! make sure the velocities are correctly replicated while using groups
    CALL mp_bcast(velp,SIZE(velp),parai%io_source,parai%cp_grp)
    CALL mp_bcast(cm,ncpw%ngw*nstate,parai%cp_inter_io_source,parai%cp_inter_grp)
    ! <<<
    ! RESET ACCUMULATORS
    IF (paral%parent.AND.irec(irec_ac).EQ.0)&
         CALL resetac(tau0,taui,iteropt%nfi)
    ! 
    IF (cntl%quenchb) THEN
       CALL quenbo(c0(:,:,1),c2,sc0,taur,rhoe,psi)
       CALL zhwwf(2,irec,c0,cm,nstate,eigv,tau0,velp,taui,iteropt%nfi)
       cntl%quenchb=.FALSE.
    ENDIF
    ! 
    ! INITIALIZE FORCES
    IF (tkpts%tkpnt) THEN
       IF (geq0) CALL zclean_k(c0,nstate,ncpw%ngw)
    ELSE
       IF (geq0) CALL zclean(c0,nstate,ncpw%ngw)
    ENDIF
    CALL forcedr(c0(:,:,1),c2,sc0,rhoe,psi,tau0,fion,eigv,&
         nstate,1,.FALSE.,.TRUE.)
    CALL totstr
    CALL freqs(nstate,.FALSE.)
    ! Check orthogonality condition for wavefunction velocities
    CALL rortv(c0,cm,c2,sc0,gamy,nstate)
    ! Initialize Barostat
    shock1%vol0 = prcp_com%omega0
    IF (prcpl%tzflex) THEN
       eps = LOG(parm%omega/prcp_com%omega0)/3._real_8
       veps = metr_com%htvel(3,3)
    ELSE IF (prcpl%tisot) THEN
       eps = LOG(parm%omega/prcp_com%omega0)/3._real_8
       veps = (metr_com%htvel(1,1)+metr_com%htvel(2,2)+metr_com%htvel(3,3))/3._real_8
    ELSE IF (cntl%tshock) THEN
       veps = metr_com%htvel(1,1)
    ENDIF
    ! Initialize thermostats
    IF (paral%parent) THEN
       itemp=irec(irec_nop1)+irec(irec_nop2)+irec(irec_nop3)
       IF (cntl%tnosep .AND. itemp.EQ.0) CALL nospinit(1)
       IF (cntl%tnosee .AND. irec(irec_noe ).EQ.0) CALL noseinit(1)
       IF (cntl%tnosec .AND. irec(irec_noc ).EQ.0) CALL noscinit(1)
       CALL wrgeof(tau0,fion)
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
    ENDIF
    CALL write_irec(irec)
    ! ==--------------------------------------------------------------==
    ! == END INITIALIZATION                                           ==
    ! ==--------------------------------------------------------------==
    IF (paral%parent) THEN
       time2 =m_walltime()
       tcpu = (time2 - time1)*0.001_real_8
       IF (paral%io_parent)&
            WRITE(6,'(A,F8.2,A)') ' TIME FOR INITIALIZATION',&
            tcpu,' SECONDS'
    ENDIF
    ! ==================================================================
    ! ==          THE BASIC LOOP FOR MOLECULAR DYNAMICS               ==
    ! ==                 USING VELOCITY VERLET                        ==
    ! ==================================================================
    DO infi=1,cnti%nomore
       time1=m_walltime()
       iteropt%nfi=iteropt%nfi+1
       ropt_mod%prteig=cprint%tprint.AND.MOD(iteropt%nfi-1,cprint%iprint_step).EQ.0
       IF (.NOT.paral%parent) ropt_mod%prteig=.FALSE.
       ropt_mod%engpri=cprint%tprint.AND.MOD(iteropt%nfi-1,cprint%iprint_step).EQ.0
       comvl%subcom=comvl%tsubcom.AND.MOD(iteropt%nfi-1,comvl%ncomv).EQ.0
       ! ANNEALING
       IF (.NOT.cntl%tshock) THEN
          CALL anneal(velp,cm,nstate,metr_com%htvel)
          CALL berendsen(velp,cm,nstate,metr_com%htvel,ekinc,ekinh)
          ! SUBTRACT CENTER OF MASS VELOCITY
          IF (paral%parent.AND.comvl%subcom) CALL comvel(velp,vcmio,.TRUE.)
       ENDIF
       ! UPDATE NOSE THERMOSTATS
       CALL noseup(velp,cm,nstate,1)
       ! UPDATE VELOCITIES
       IF (paral%parent) THEN
          IF (cntl%tshock) THEN
             ! FOR SHOCKS, APPLY DAMPING TO ION VELOCITIES
             CALL velupidamp(velp,1)
             ! FOR SHOCKS, APPLY DAMPING TO BAROSTAT VELOCITIES
             CALL vepsupdamp
             ! FOR SHOCKS, UPDATE BAROSTAT AND PARTICLE VELOCITIES
             CALL vepsup(velp)
             ! FOR SHOCKS, APPLY DAMPING TO BAROSTAT VELOCITIES
             CALL vepsupdamp
             ! FOR SHOCKS, UPDATE VELOCITIES
             CALL velupishock(velp,fion,1)
          ELSE
             ! IONS and CELL
             CALL velupi(velp,fion,1)
             CALL velupc(1)
          ENDIF
       ENDIF
       CALL velupa(c0,cm,c2,nstate,1)
       ! UPDATE POSITIONS
       IF (prcpl%tzflex) THEN
          ! AK: FIXME
          CALL stopgm('NPT_CPMD','NO NPT MD FOR Z-FLEXIBLE CELL',& 
               __LINE__,__FILE__)
          ! CALL POSUPIH_ZSC(TAU0,TAUP,VELP)
       ELSE IF (prcpl%tisot) THEN
          CALL posupih_iso(tau0,taup,velp)
       ELSE IF (cntl%tshock) THEN
          CALL posupihshock(tau0,taup,velp)
       ELSE
          CALL posupih(tau0,taup,velp)
       ENDIF
       CALL mp_bcast(metr_com%ht, SIZE(metr_com%ht),parai%io_source,parai%cp_grp)
       DO i=1,3
          parm%a1(i) = metr_com%ht(1,i)
          parm%a2(i) = metr_com%ht(2,i)
          parm%a3(i) = metr_com%ht(3,i)
       ENDDO
       CALL ihmat(metr_com%ht,metr_com%htm1,parm%omega)
       IF (paral%parent) THEN
          CALL dcopy(3*maxsys%nax*maxsys%nsx,taup(1,1,1),1,tau0(1,1,1),1)
       ENDIF
       CALL mp_bcast(taup,SIZE(taup),parai%io_source,parai%cp_grp)
       CALL newcell
       CALL phfac(taup)
       IF (corel%tinlc) CALL copot(rhoe,psi,.TRUE.)
       CALL posupa(c0,cm,c2,gamx,nstate)
       ! mb - Wannier stuff for vdW-WC
       IF (vdwl%vdwd.AND.vdwwfl%twannup) THEN
          CALL localize2(taup,c0,c2,sc0,nstate)
       ENDIF
       ! CALCULATE THE FORCES
       CALL forcedr(c0(:,:,1),c2,sc0,rhoe,psi,taup,fion,eigv,&
            nstate,1,.FALSE.,.TRUE.)
       CALL totstr
       ! ==================================================================
       ! Damped Dynamics
       CALL dampdyn(velp,fion,cm,c2,nstate,metr_com%htvel,metr_com%htfor)
       ! ==================================================================
       ! FINAL UPDATE FOR VELOCITIES
       IF (paral%parent) THEN
          IF (cntl%tshock) THEN
             CALL velupishock(velp,fion,1)
          ELSE
             CALL velupc(1)
             CALL velupi(velp,fion,1)
          ENDIF
       ENDIF
       CALL velupa(c0,cm,c2,nstate,1)
       CALL rortv(c0,cm,c2,sc0,gamy,nstate)
       IF (paral%parent) THEN
          CALL geofile(taup,velp,'WRITE')
          ! COMPUTE THE IONIC TEMPERATURE TEMPP
          CALL ekinpp(ekinp,velp)
          tempp=ekinp*factem*2._real_8/glib
          ! IONIC TEMPERATURE CONTROL
          CALL rscvp(temp1,temp2,tempp,velp)
       ENDIF
       ! UPDATE NOSE THERMOSTATS
       IF (cntl%tnosee.OR.cntl%tc) CALL rekine(cm,nstate,ekinc)
       IF (paral%parent) THEN
          IF (cntl%tshock) THEN
             ! FOR SHOCKS, APPLY DAMPING TO BAROSTAT VELOCITIES
             CALL vepsupdamp
             ! FOR SHOCKS, UPDATE BAROSTAT VELOCITIES
             CALL vepsup(velp)
             ! FOR SHOCKS, APPLY DAMPING TO BAROSTAT VELOCITIES
             CALL vepsupdamp
             ! FOR SHOCKS, APPLY DAMPING TO ION VELOCITIES
             CALL velupidamp(velp,1)
          ELSE
             ! IONS and CELL
          !  CALL velupi(velp,fion,1) !bugfix
          !  CALL velupc(1)
          ENDIF
       ENDIF
       ! THERMOSTAT UPDATE
       CALL noseup(velp,cm,nstate,1)
       IF (cntl%tshock) THEN
          metr_com%htvel(1,1) = veps
          metr_com%htvel(2,2) = 0._real_8
          metr_com%htvel(3,3) = 0._real_8
       ELSE
          ! SUBTRACT CENTER OF MASS VELOCITY
          IF (paral%parent.AND.comvl%subcom) CALL comvel(velp,vcmio,.FALSE.)
          CALL berendsen(velp,cm,nstate,metr_com%htvel,ekinc,ekinh)
          ! ANNEALING
          CALL anneal(velp,cm,nstate,metr_com%htvel)
       ENDIF
       IF (paral%parent) THEN
          CALL ekinpp(ekinp,velp)
          tempp=ekinp*factem*2._real_8/glib
       ENDIF
       ! RESCALE ELECTRONIC VELOCITIES
       IF (cntl%tc) CALL rscve(ekin1,ekin2,ekinc,cntr%ekinw,cm,nstate,ncpw%ngw)
       ! MEAN SQUARE DISPLACEMENT OF DIFFERENT IONIC SPECIES
       IF (paral%parent) CALL dispp(taup,taui,disa)
       ! KINETIC ENERGY OF THE ELECTRONS
       CALL rekine(cm,nstate,ekinc)
       ! KINETIC ENERGY OF THE CELL
       ekinh=0.5_real_8*cntr%cmass*ddot(9,metr_com%htvel,1,metr_com%htvel,1)
       IF (cntl%tshock) ekinh=0.5_real_8*cntr%cmass*veps*veps
       ! TEMPERATURE OF THE CELL
       temph = 2._real_8*ekinh*factem/REAL(ncdof,kind=real_8)
       ! ENERGY OF THE NOSE THERMOSTATS
       IF (paral%parent) CALL noseng(iteropt%nfi,velp,enose,enosp,enosc,1)
       IF (paral%parent) THEN
          IF (.NOT.cntl%tshock) THEN
             ! Conserved quantity for NPT ensembles
             econs=ekinp+ener_com%etot+enose+enosp+enosc+ener_com%ecnstr+ener_com%erestr+prcp_com%druck*parm%omega
             eham=econs+ekinc+ekinh
          ELSE
             ! Conserved quantity (HUGONIOT) for NPH-SHOCK ensemble
             shock1%eshock=-0.5_real_8*shock1%vshock*shock1%vshock*(1.0_real_8-parm%omega/shock1%vol0)&
                  *(1.0_real_8-parm%omega/shock1%vol0)-shock1%pshock*(shock1%vol0-parm%omega)
             econs=ekinp+ener_com%etot+enose+ener_com%ecnstr+ener_com%erestr+shock1%eshock
             eham=econs+ekinc+ekinh
          ENDIF
          time2=m_walltime()
          tcpu=(time2-time1)*0.001_real_8
          ! PRINTOUT the evolution of the accumulators every time step
          CALL wrprint_pr(ekinc,ekinh,tempp,ener_com%etot,econs,eham,&
               disa,taup,fion,tcpu,iteropt%nfi,infi)
          ! UPDATE ACCUMULATORS
          CALL paccb(ekinc,tempp,ener_com%etot,econs,eham,disa,tcpu,temph,&
               parm%omega,enose,enosp,iteropt%nfi,1)
          ! Store ionic coordinates and velocities for statistics
          ropt_mod%movie=rout1%mout .AND. MOD(iteropt%nfi-1,cnti%imovie).EQ.0
          ropt_mod%rprint=rout1%rout .AND. MOD(iteropt%nfi-1,cnti%ntraj).EQ.0
          ropt_mod%txyz=rout1%xtout .AND. MOD(iteropt%nfi-1,cnti%ntraj).EQ.0
          ropt_mod%tdcd=rout1%dcout .AND. MOD(iteropt%nfi-1,cnti%ntraj).EQ.0
          CALL printp(taur,taup,velp)
       ENDIF
       IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)
       IF (infi.EQ.cnti%nomore) THEN
          soft_com%exsoft=.TRUE.
          soft_com%exnomore=.TRUE.
       ENDIF
       ! periodic output of density/wavefunction etc.
       IF (rout1%rhoout.AND.(rout1%nrhoout.GT.0)) THEN
          IF (MOD(iteropt%nfi-1,rout1%nrhoout).EQ.0) THEN
             CALL rhopri(c0,tau0,rhoe,psi(:,1),nstate,nkpt%nkpnt)
          ENDIF
       ENDIF
       IF (teststore(iteropt%nfi).OR.soft_com%exsoft)&
            CALL zhwwf(2,irec,c0,cm,nstate,eigv,taup,velp,taui,iteropt%nfi)
       ! temperature ramping
       CALL tempramp(temp1,temp2)
       ! STOP THE RUN IF THE USER HAS SET THE SIGNAL 30
       IF (soft_com%exsoft) GOTO 100
       ! UPDATE IONIC POSITIONS
       CALL dcopy(3*maxsys%nax*maxsys%nsx,taup(1,1,1),1,tau0(1,1,1),1)
       ! ==================================================================
       ! ==     END OF MAIN LOOP                                         ==
       ! ==================================================================
    ENDDO
100 CONTINUE
    IF (rout1%rhoout.AND.(rout1%nrhoout.LE.0))&
         CALL rhopri(c0,tau0,rhoe,psi(:,1),nstate,nkpt%nkpnt)
    ! Print accumulators.
    IF (paral%parent) CALL paccb(ekinc,tempp,ener_com%etot,econs,eham,disa,tcpu,temph&
         ,parm%omega,enose,enosp,iteropt%nfi,2)
    CALL proja(c0,c2,sc0,nstate,cnti%iproj)
    CALL csize(c2,nstate,gemax,cnorm)
    IF (paral%parent) CALL gsize(fion,gnmax,gnorm)
    IF (paral%parent) CALL finalp(tau0,fion,velp,eigv)
    IF (pslo_com%tivan.OR.cntl%nonort) THEN
       DEALLOCATE(eigv,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(eigm,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    DEALLOCATE(taup,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(fion,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(velp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(taui,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(taur,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF ((paral%parent).AND.paral%io_parent)&
         CALL fileclose(3)
    DEALLOCATE(rhoe,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(psi,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(scr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE npt_cpmd
  ! ==================================================================
  SUBROUTINE give_scr_npt_md(lnpt_md,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lnpt_md
    CHARACTER(len=30)                        :: tag

    INTEGER :: lcopot, lddipo, ldeort, lforcedr, linitrun, lnewcell, lortho, &
      lposupa, lquenbo, lrhopri, lrortv, nstate

    nstate=crge%n
    lddipo=0
    lcopot=0
    lortho=0
    lquenbo=0
    ldeort=0
    lrhopri=0
    linitrun=0
    CALL give_scr_initrun(linitrun,tag)
    IF (corel%tinlc) CALL give_scr_copot(lcopot,tag)
    IF (cntl%trane) CALL give_scr_ortho(lortho,tag,nstate)
    IF (cntl%quenchb) CALL give_scr_quenbo(lquenbo,tag)
    IF (pslo_com%tivan) CALL give_scr_deort(ldeort,tag,nstate)
    CALL give_scr_forcedr(lforcedr,tag,nstate,.FALSE.,.TRUE.)
    CALL give_scr_rortv(lrortv,tag,nstate)
    CALL give_scr_newcell(lnewcell,tag)
    CALL give_scr_posupa(lposupa,tag,nstate)
    IF (rout1%rhoout) CALL give_scr_rhopri(lrhopri,tag,nstate)
    IF (vdwl%vdwd) CALL give_scr_ddipo(lddipo,tag)
    lnpt_md=MAX(lcopot,lortho,lquenbo,ldeort,lforcedr,&
         lrortv,lposupa,lrhopri,linitrun,&
         lnewcell,lddipo)
    lnpt_md=lnpt_md+10000
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_npt_md
  ! ==================================================================
  SUBROUTINE npt_bomd(c0,cm,c2,sc0,vpp,gamx,gamy)
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8)                          :: c0(ncpw%ngw,crge%n,1), cm(:), &
                                                c2(ncpw%ngw,crge%n), &
                                                sc0(ncpw%ngw,crge%n)
    REAL(real_8)                             :: vpp(*), gamx(*), gamy(*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'npt_bomd'

    CHARACTER(len=100)                       :: filen
    CHARACTER(len=30)                        :: tag
    COMPLEX(real_8), ALLOCATABLE             :: psi(:,:)
    INTEGER :: i, ierr, ifcalc, il_psi_1d, il_psi_2d, il_rhoe_1d, il_rhoe_2d, &
      irec(100), itemp, j, loopnfi, lscr, nfimax, nfimin, nnx, nstate, nx
    LOGICAL                                  :: ferror
    REAL(real_8) :: disa, econs, eham, ekin1, ekin2, ekinc, ekincp, ekinh, &
      ekinh1, ekinh2, ekinp, enosc, enose, enosp, tcpu, temp1, temp2, temph, &
      tempp, time1, time2, vcmio(4)
    REAL(real_8), ALLOCATABLE                :: eigv(:), rhoe(:,:), rinp(:), &
                                                rm1(:), scr(:), taui(:,:,:), &
                                                taur(:,:,:)
    REAL(real_8), EXTERNAL                   :: ddot

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
    CALL zeroing(taup)!,SIZE(taup))
    CALL zeroing(fion)!,SIZE(fion))
    CALL zeroing(taui)!,SIZE(taui))
    CALL zeroing(taur)!,SIZE(taur))
    ! ==--------------------------------------------------------------==
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
    CALL zeroing(rm1)
    ! ==--------------------------------------------------------------==
    nstate=crge%n
    nacc = 22
    iteropt%nfi  = 0
    ropt_mod%modens=.FALSE.
    ropt_mod%engpri=.FALSE.
    ropt_mod%calste=.TRUE.
    ! ==--------------------------------------------------------------==
    ! 06/2008 AK
    ! NPT FLEXIBLE CELL DOES NOT WORK WITH VANDERBILT PSEUDOPOTENTIALS.
    ! FOR STARTERS THE STRESS TENSOR IS BROKEN.
    ! IT WORKS IN QUANTUM ESPRESSO...
    ! 08/2017 MB-TI - The above statement is wrong. Fixed on 28/08/2017.
  ! IF (pslo_com%tivan) THEN
  !    CALL stopgm('NPT_BOMD','NPT MD IS BROKEN FOR VANDERBILT USPPS',& 
  !         __LINE__,__FILE__)
  ! ENDIF
    ! ==--------------------------------------------------------------==
    ALLOCATE(eigv(clsd%nlsd*nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! Extrapolation
    IF (cntl%textrap) THEN
       ALLOCATE(cold(nkpt%ngwk,crge%n,nkpt%nkpnt,cnti%mextra*nkpt%nkpts*nstate/(crge%n*nkpt%nkpnt)),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == INITIALIZATION                                               ==
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
    CALL give_scr_npt_md(lscr,tag)
    ALLOCATE(scr(lscr),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    restf%nfnow=1
    ! ==--------------------------------------------------------------==
    ! TIME STEP FUNCTIONS
    CALL dynit(ekincp,ekin1,ekin2,temp1,temp2,ekinh1,ekinh2)
    ekinc=0.0_real_8
    ekinp=0.0_real_8
    ! 
    CALL initrun(irec,c0,c2,sc0,rhoe,psi,eigv)
    CALL mp_bcast(taup,SIZE(taup),parai%io_source,parai%cp_grp)
    CALL dcopy(3*maxsys%nax*maxsys%nsx,taup(1,1,1),1,taui(1,1,1),1)
    ncdof = 1
    IF ((cntl%tnosep.OR.cntl%tnosee.OR.cntl%tnosec).AND.paral%parent) CALL nosepa(1,1)
    ! INITIALIZE VELOCITIES
    IF (paral%parent) CALL detdof(tau0,taur)
    IF (irec(irec_vel).EQ.0.AND..NOT.restart1%rgeo) THEN
       ener_com%ecnstr = 0.0_real_8
       ener_com%erestr = 0.0_real_8
       CALL rinvel(velp,cm,nstate)
       IF (paral%parent) CALL taucl(velp)
       IF (paral%parent) CALL rattle(tau0,velp)
       CALL rvscal(velp)
    ELSE
       IF (paral%parent) CALL taucl(velp)
       IF (paral%parent) CALL rattle(tau0,velp)
       IF (cntl%trescale) CALL rvscal(velp)
    ENDIF
    IF (cntl%quenchp) CALL zeroing(velp)!,3*maxsys%nax*maxsys%nsx)
    IF (cntl%quenchc) CALL zeroing(metr_com%htvel)!,9)
    IF (cntl%trevers) THEN
       ! invert ionic velocities (useful for path sampling)
       CALL dscal(3*maxsys%nax*maxsys%nsx,-1._real_8,velp,1)
    ENDIF
    ! COMPUTE THE IONIC TEMPERATURE TEMPP
    IF (paral%parent) THEN
       CALL ekinpp(ekinp,velp)
       tempp=ekinp*factem*2._real_8/glib
    ENDIF
    ! RESET ACCUMULATORS
    IF (paral%parent.AND.irec(irec_ac).EQ.0)&
         CALL resetac(tau0,taui,iteropt%nfi)
    ! INITIALIZE FORCES
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
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,'(1X,64("="))')
       IF (paral%io_parent)&
            WRITE(6,'(1X,"==",T25,A,T64,"==")')&
            'FORCES INITIALIZATION'
       IF (paral%io_parent)&
            WRITE(6,'(1X,64("="))')
    ENDIF
    ! INITIALIZE WF CENTERS & SPREAD
    IF(vdwl%vdwd.AND..NOT.vdwwfl%trwannc) THEN
       CALL localize2(tau0,c0,c2,sc0,nstate)
    ENDIF
    ifcalc=0
    CALL forces_diag(nstate,c0,c2,cm,sc0,cm(nx:),vpp,eigv,&
         rhoe,psi,&
         tau0,velp,taui,fion,ifcalc,&
         irec,.TRUE.,.TRUE.)
    CALL totstr
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
    ! Initialize thermostats
    IF (paral%parent) THEN
       itemp=irec(irec_nop1)+irec(irec_nop2)+irec(irec_nop3)
       IF (cntl%tnosep .AND. itemp.EQ.0) CALL nospinit(1)
       IF (cntl%tnosec .AND. irec(irec_noc ).EQ.0) CALL noscinit(1)
       CALL wrgeof(tau0,fion)
       filen='ENERGIES'
       IF (paral%io_parent)&
            CALL fileopen(3,filen,fo_app+fo_verb,ferror)
    ENDIF
    CALL write_irec(irec)
    ! ==--------------------------------------------------------------==
    ! == END INITIALIZATION                                           ==
    ! ==--------------------------------------------------------------==
    IF (paral%parent) THEN
       time2 =m_walltime()
       tcpu = (time2 - time1)*0.001_real_8
       IF (paral%io_parent)&
            WRITE(6,'(A,F8.2,A)') ' TIME FOR INITIALIZATION',&
            tcpu,' SECONDS'
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
    ENDIF
    ! ==================================================================
    ! ==          THE BASIC LOOP FOR MOLECULAR DYNAMICS               ==
    ! ==                 USING VELOCITY VERLET                        ==
    ! ==================================================================
    infi=0
    nfimin=iteropt%nfi+1
    nfimax=iteropt%nfi+cnti%nomore
    DO loopnfi=nfimin,nfimax
       time1=m_walltime()
       infi=infi+1
       iteropt%nfi=loopnfi
       comvl%subcom=comvl%tsubcom.AND.MOD(iteropt%nfi-1,comvl%ncomv).EQ.0
       ropt_mod%prteig=cprint%tprint.AND.MOD(iteropt%nfi-1,cprint%iprint_step).EQ.0
       IF (.NOT.paral%parent) ropt_mod%prteig=.FALSE.
       ropt_mod%engpri=cprint%tprint.AND.MOD(iteropt%nfi-1,cprint%iprint_step).EQ.0
       ! ANNEALING
       IF (.NOT.cntl%tshock) THEN
          CALL anneal(velp,cm,nstate,metr_com%htvel)
          CALL berendsen(velp,cm,nstate,metr_com%htvel,ekinc,ekinh)
          ! SUBTRACT CENTER OF MASS VELOCITY
          IF (paral%parent.AND.comvl%subcom) CALL comvel(velp,vcmio,.TRUE.)
       ENDIF
       ! UPDATE NOSE THERMOSTATS
       CALL noseup(velp,cm,nstate,1)
       IF (paral%parent) THEN
          IF (cntl%tshock) THEN
             ! FOR SHOCKS, APPLY DAMPING TO ION VELOCITIES
             CALL velupidamp(velp,1)
             ! FOR SHOCKS, APPLY DAMPING TO BAROSTAT VELOCITIES
             CALL vepsupdamp
             ! FOR SHOCKS, UPDATE BAROSTAT AND PARTICLE VELOCITIES
             CALL vepsup(velp)
             ! FOR SHOCKS, APPLY DAMPING TO BAROSTAT VELOCITIES
             CALL vepsupdamp
             ! FOR SHOCKS, UPDATE VELOCITIES
             CALL velupishock(velp,fion,1)
          ELSE
             ! IONS and CELL
             CALL velupi(velp,fion,1)
             CALL velupc(1)
          ENDIF
       ENDIF
       ! UPDATE POSITIONS
       IF (prcpl%tzflex) THEN
          ! AK: FIXME
          CALL stopgm('NPT_BOMD','NO NPT MD FOR Z-FLEXIBLE CELL',& 
               __LINE__,__FILE__)
          ! CALL POSUPIH_ZSC(TAU0,TAUP,VELP)
          ! UPDATE CELL
       ELSE IF (prcpl%tisot) THEN
          CALL posupih_iso(tau0,taup,velp)
       ELSE IF (cntl%tshock) THEN
          CALL posupihshock(tau0,taup,velp)
       ELSE
          CALL posupih(tau0,taup,velp)
       ENDIF
       CALL mp_bcast(metr_com%ht, SIZE(metr_com%ht),parai%io_source,parai%cp_grp)
       DO i=1,3
          parm%a1(i) = metr_com%ht(1,i)
          parm%a2(i) = metr_com%ht(2,i)
          parm%a3(i) = metr_com%ht(3,i)
       ENDDO
       CALL ihmat(metr_com%ht,metr_com%htm1,parm%omega)
       IF (paral%parent) THEN
          CALL dcopy(3*maxsys%nax*maxsys%nsx,taup(1,1,1),1,tau0(1,1,1),1)
       ENDIF
       CALL mp_bcast(taup,SIZE(taup),parai%io_source,parai%cp_grp)
       CALL newcell
       CALL phfac(taup)
       IF (corel%tinlc) CALL copot(rhoe,psi,.TRUE.)
       !mb - Wannier stuff for vdW-WC
       vdwwfl%twannup=vdwwfl%twannup.OR.(infi.EQ.1.AND..NOT.vdwwfl%trwannc)
       IF(vdwl%vdwd.AND.vdwwfl%twannup) THEN
          CALL localize2(taup,c0,c2,sc0,nstate)
       ENDIF
       ! CALCULATE THE FORCES
       IF (cntl%textrap) THEN
          ! Extrapolate wavefunctions
          CALL extrapwf(infi,c0,scr,cold,nnow,numcold,nstate,cnti%mextra)
       ENDIF
       ifcalc=0
       CALL forces_diag(nstate,c0,c2,cm,sc0,cm(nx:),vpp,eigv,&
            rhoe,psi,&
            taup,velp,taui,fion,ifcalc,&
            irec,.TRUE.,.FALSE.)
       CALL totstr
       ! ==================================================================
       ! Damped Dynamics
       CALL dampdyn(velp,fion,cm,c2,nstate,metr_com%htvel,metr_com%htfor)
       ! ==================================================================
       ! FINAL UPDATE FOR VELOCITIES
       IF (paral%parent) THEN
          IF (cntl%tshock) THEN
             CALL velupishock(velp,fion,1)
          ELSE
             CALL velupc(1)
             CALL velupi(velp,fion,1)
          ENDIF
          CALL geofile(taup,velp,'WRITE')
          ! COMPUTE THE IONIC TEMPERATURE TEMPP
          CALL ekinpp(ekinp,velp)
          tempp=ekinp*factem*2._real_8/glib
          ! IONIC TEMPERATURE CONTROL
          CALL rscvp(temp1,temp2,tempp,velp)
          IF (cntl%tshock) THEN
             ! FOR SHOCKS, APPLY DAMPING TO BAROSTAT VELOCITIES
             CALL vepsupdamp
             ! FOR SHOCKS, UPDATE BAROSTAT VELOCITIES
             CALL vepsup(velp)
             ! FOR SHOCKS, APPLY DAMPING TO BAROSTAT VELOCITIES
             CALL vepsupdamp
             ! FOR SHOCKS, APPLY DAMPING TO ION VELOCITIES
             CALL velupidamp(velp,1)
          ELSE
             ! IONS and CELL
          !  CALL velupi(velp,fion,1) !bugfix
          !  CALL velupc(1)
          ENDIF
       ENDIF
       ! THERMOSTAT UPDATE
       CALL noseup(velp,cm,nstate,1)
       IF (cntl%tshock) THEN
          metr_com%htvel(1,1) = veps
          metr_com%htvel(2,2) = 0._real_8
          metr_com%htvel(3,3) = 0._real_8
       ELSE
          ! SUBTRACT CENTER OF MASS VELOCITY
          IF (paral%parent.AND.comvl%subcom) CALL comvel(velp,vcmio,.FALSE.)
          CALL berendsen(velp,cm,nstate,metr_com%htvel,0.0_real_8,ekinh)
          ! ANNEALING
          CALL anneal(velp,cm,nstate,metr_com%htvel)
       ENDIF
       IF (paral%parent) THEN
          CALL ekinpp(ekinp,velp)
          tempp=ekinp*factem*2._real_8/glib
       ENDIF
       ! MEAN SQUARE DISPLACEMENT OF DIFFERENT IONIC SPECIES
       IF (paral%parent) CALL dispp(taup,taui,disa)
       ! KINETIC ENERGY OF THE CELL
       ekinh=0.5_real_8*cntr%cmass*ddot(9,metr_com%htvel,1,metr_com%htvel,1)
       IF (cntl%tshock) ekinh=0.5_real_8*cntr%cmass*veps*veps
       ! TEMPERATURE OF THE CELL
       temph = 2._real_8*ekinh*factem/REAL(ncdof,kind=real_8)
       ! ENERGY OF THE NOSE THERMOSTATS
       IF (paral%parent) CALL noseng(iteropt%nfi,velp,enose,enosp,enosc,1)
       IF (paral%parent) THEN
          IF (cntl%tshock) THEN
             ! Conserved quantity (HUGONIOT) for NPH-SHOCK ensemble
             shock1%eshock=-0.5_real_8*shock1%vshock*shock1%vshock*(1.0_real_8-parm%omega/shock1%vol0)&
                  *(1.0_real_8-parm%omega/shock1%vol0)-shock1%pshock*(shock1%vol0-parm%omega)
             econs=ekinp+ener_com%etot+ener_com%ecnstr+ener_com%erestr+shock1%eshock
          ELSE
             !eham=econs+ekinh ! econs is not initialized and eham not used....
             ! Conserved quantity for regular NPT ensembles
             econs=ekinp+ener_com%etot+enosc+enosp+ener_com%ecnstr+ener_com%erestr+prcp_com%druck*parm%omega
             eham=econs+ekinh
          ENDIF
          time2=m_walltime()
          tcpu=(time2-time1)*0.001_real_8
          ! PRINTOUT the evolution of the accumulators every time step
          CALL wrprint_pr(ekinc,ekinh,tempp,ener_com%etot,econs,eham,&
               disa,taup,fion,tcpu,iteropt%nfi,infi)
          ! UPDATE ACCUMULATORS
          CALL paccb(ekinc,tempp,ener_com%etot,econs,eham,disa,tcpu,temph,&
               parm%omega,0._real_8,enosp,iteropt%nfi,1)
          ! Store ionic coordinates and velocities for statistics
          ropt_mod%movie=rout1%mout .AND. MOD(iteropt%nfi-1,cnti%imovie).EQ.0
          ropt_mod%rprint=rout1%rout .AND. MOD(iteropt%nfi-1,cnti%ntraj).EQ.0
          ropt_mod%txyz=rout1%xtout .AND. MOD(iteropt%nfi-1,cnti%ntraj).EQ.0
          ropt_mod%tdcd=rout1%dcout .AND. MOD(iteropt%nfi-1,cnti%ntraj).EQ.0
          CALL printp(taur,taup,velp)
       ENDIF
       IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)
       IF (iteropt%nfi.EQ.nfimax) THEN
          soft_com%exsoft=.TRUE.
          soft_com%exnomore=.TRUE.
       ENDIF
       ! periodic output of density/wavefunction etc.
       IF (rout1%rhoout.AND.(rout1%nrhoout.GT.0)) THEN
          IF (MOD(iteropt%nfi-1,rout1%nrhoout).EQ.0) THEN
             CALL rhopri(c0,tau0,rhoe,psi(:,1),nstate,nkpt%nkpnt)
          ENDIF
       ENDIF
       IF (teststore(iteropt%nfi).OR.soft_com%exsoft)&
            CALL zhwwf(2,irec,c0,cm,nstate,eigv,taup,velp,taui,iteropt%nfi)
       ! temperature ramping
       CALL tempramp(temp1,temp2)
       ! STOP THE RUN IF THE USER HAS SET THE SIGNAL 30
       IF (soft_com%exsoft) GOTO 100
       ! UPDATE IONIC POSITIONS
       CALL dcopy(3*maxsys%nax*maxsys%nsx,taup(1,1,1),1,tau0(1,1,1),1)
       ! UPDATE DENSITY
       IF (.NOT.cntl%bsymm) THEN
          CALL extrap(nnx,andr2%alxmix,rm1,rin0,rinp)
          CALL dcopy(nnx,rin0,1,rm1(1),1)
          CALL dcopy(nnx,rinp(1),1,rin0,1)
       ENDIF
       ! ==================================================================
       ! ==     END OF MAIN LOOP                                         ==
       ! ==================================================================
    ENDDO
100 CONTINUE
    IF (rout1%rhoout.AND.(rout1%nrhoout.LE.0))&
         CALL rhopri(c0,tau0,rhoe,psi(:,1),nstate,nkpt%nkpnt)
    ! Print accumulators.
    IF (paral%parent) CALL paccb(ekinc,tempp,ener_com%etot,econs,eham,disa,tcpu,temph&
         ,parm%omega,enose,enosp,iteropt%nfi,2)
    CALL proja(c0,c2,sc0,nstate,cnti%iproj)
    CALL csize(c2,nstate,gemax,cnorm)
    IF (paral%parent) CALL gsize(fion,gnmax,gnorm)
    IF (paral%parent) CALL finalp(tau0,fion,velp,eigv)
    DEALLOCATE(taup,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(fion,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(velp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(taui,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(taur,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF ((paral%parent).AND.paral%io_parent)&
         CALL fileclose(3)
    DEALLOCATE(rhoe,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(psi,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(scr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(eigv,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (cntl%textrap) DEALLOCATE(cold,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE npt_bomd
  ! ==================================================================

END MODULE npt_md_utils
