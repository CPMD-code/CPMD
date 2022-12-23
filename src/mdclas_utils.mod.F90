MODULE mdclas_utils
  USE anneal_utils,                    ONLY: anneal,&
                                             berendsen,&
                                             dampdyn,&
                                             tempramp
  USE clas,                            ONLY: &
       clab, clas3, clas4, clas7, clasc, clasf, clasfold, clasv, defy, &
       delclasc, imovc, ndefo, ntrac
  USE clas_force_utils,                ONLY: check,&
                                             clas_force
  USE cnst,                            ONLY: factem,&
                                             fbohr,&
                                             pi,&
                                             scmass
  USE comvel_utils,                    ONLY: comvel
  USE comvelmod,                       ONLY: comvl
  USE coor,                            ONLY: fion,&
                                             tau0,&
                                             taup,&
                                             velp
  USE copot_utils,                     ONLY: copot
  USE cotr,                            ONLY: cotc0
  USE csize_utils,                     ONLY: csize
  USE ddipo_utils,                     ONLY: ddipo
  USE deort_utils,                     ONLY: deort
  USE detdof_utils,                    ONLY: detdof
  USE dispp_utils,                     ONLY: dispp
  USE dynit_utils,                     ONLY: dynit
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_app,&
                                             fo_mark,&
                                             fo_ufo,&
                                             fo_verb
  USE finalp_utils,                    ONLY: finalp
  USE forcedr_driver,                  ONLY: forcedr
  USE forcep_utils,                    ONLY: rhoe_psi_size
  USE freqs_utils,                     ONLY: freqs
  USE geofile_utils,                   ONLY: geofile
  USE geq0mod,                         ONLY: geq0
  USE gsize_utils,                     ONLY: gsize
  USE initrun_driver,                  ONLY: initrun
  USE ions,                            ONLY: ions0
  USE kinds,                           ONLY: real_4,&
                                             real_8
  USE machine,                         ONLY: m_walltime
  USE mdmain_utils,                    ONLY: give_scr_mdmain
  USE metr,                            ONLY: metr_com
  USE movi,                            ONLY: imtyp
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sync
  USE nlcc,                            ONLY: corel
  USE norm,                            ONLY: cnorm,&
                                             gemax,&
                                             gnmax,&
                                             gnorm
  USE noseinit_utils,                  ONLY: noseinit
  USE noseng_utils,                    ONLY: noseng
  USE nosepa_utils,                    ONLY: nosepa
  USE noseup_utils,                    ONLY: noseup
  USE nospinit_utils,                  ONLY: nospinit
  USE parac,                           ONLY: parai,&
                                             paral
  USE phfac_utils,                     ONLY: phfac
  USE posupa_utils,                    ONLY: posupa
  USE printave_utils,                  ONLY: pacca
  USE printp_utils,                    ONLY: printp
  USE prng_utils,                      ONLY: repprngu_vec
  USE proja_utils,                     ONLY: proja
  USE pslo,                            ONLY: pslo_com
  USE quenbo_utils,                    ONLY: quenbo
  USE rattle_utils,                    ONLY: rattle
  USE rekine_utils,                    ONLY: rekine
  USE resetac_utils,                   ONLY: resetac
  USE rhopri_utils,                    ONLY: rhopri
  USE ropt,                            ONLY: infi,&
                                             iteropt,&
                                             ropt_mod
  USE rortv_utils,                     ONLY: rortv
  USE rotvel_utils,                    ONLY: rotvel
  USE rscve_utils,                     ONLY: rscve
  USE setirec_utils,                   ONLY: write_irec
  USE shake_utils,                     ONLY: cpmdshake
  USE soft,                            ONLY: soft_com
  USE spin,                            ONLY: spin_mod
  USE store_types,                     ONLY: &
       cprint, irec_ac, irec_noe, irec_nop1, irec_nop2, irec_nop3, irec_nop4, &
       irec_vel, restart1, rout1, store1
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             cntr,&
                                             maxsys,&
                                             nacc,&
                                             ncpw,&
                                             nkpt
  USE temps,                           ONLY: tempcm,&
                                             tempqm
  USE testex_utils,                    ONLY: testex
  USE totstr_utils,                    ONLY: totstr
  USE tpar,                            ONLY: dt_ions
  USE utils,                           ONLY: zclean
  USE velupa_utils,                    ONLY: velupa
  USE wannier_print_utils,             ONLY: wannier_print
  USE wrener_utils,                    ONLY: wrprint_md
  USE wrgeo_utils,                     ONLY: wrgeof
  USE wv30_utils,                      ONLY: zhwwf
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mdclas
  !public :: q_coor
  !public :: add_qm_force
  !public :: init_velocities
  !public :: ekincl
  !public :: deformation
  !public :: temp_control
  !public :: wr_classic
  !public :: wr_cltra
  !public :: velocities

CONTAINS

  ! ==================================================================
  SUBROUTINE mdclas(c0,cm,c2,sc0,gamx,gamy)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8) :: c0(ncpw%ngw,crge%n,1), cm(ncpw%ngw,crge%n), &
      c2(ncpw%ngw,crge%n), sc0(ncpw%ngw,crge%n)
    REAL(real_8)                             :: gamx(*), gamy(*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'mdclas'

    CHARACTER(len=100)                       :: filen
    CHARACTER(len=30)                        :: tag
    COMPLEX(real_8), ALLOCATABLE             :: psi(:,:)
    INTEGER :: i, ia, ic, ierr, il_psi_1d, il_psi_2d, il_rhoe_1d, il_rhoe_2d, &
      irec(100), is, itemp, j, lscr, nc, nq, nstate
    LOGICAL                                  :: ferror, update
    REAL(real_8) :: alfap, disa, dummy, econs, eham, eint, ekin1, ekin2, &
      ekinc, ekincm, ekincp, ekinh1, ekinh2, ekinp, ekinqm, enose, enosp, &
      extwork, ff, lmio(3), tcpu, temp1, temp2, tempp, time1, time2, vcmio(4)
    REAL(real_8), ALLOCATABLE                :: center(:,:), eigm(:), &
                                                eigv(:), rhoe(:,:), scr(:), &
                                                taui(:,:,:), tauio(:,:), &
                                                taur(:,:,:)

! Variables
! ==================================================================

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
    IF (comvl%tsubrot)  THEN
       ALLOCATE(tauio(3,maxsys%nax),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)! TODO check dimensions
    ENDIF
    CALL zeroing(taup)!,SIZE(taup))
    CALL zeroing(fion)!,SIZE(fion))

    nstate=crge%n
    nacc = 22
    ropt_mod%modens=.FALSE.
    ropt_mod%engpri=.FALSE.
    ropt_mod%calste=cntl%tpres
    ropt_mod%calstc=cntl%tprec
    IF (pslo_com%tivan.OR.cntl%nonort) THEN
       ALLOCATE(eigv(nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(eigm(nstate*nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    IF (cntl%tharm.AND.paral%parent) THEN
       ff=crge%f(1,1)
       DO i=1,nstate
          IF (ff.NE.crge%f(i,1)) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(A,A)') ' HARMONIC REFERENCE SYSTEM INTEGRATION',&
                  ' ONLY POSSIBLE WITH EQUAL OCCUPATION NUMBERS'
             CALL stopgm('MDCLAS',' ',& 
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
    CALL give_scr_mdmain(lscr,tag)
    ALLOCATE(scr(lscr),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! TIME STEP FUNCTIONS
    CALL dynit(ekincp,ekin1,ekin2,temp1,temp2,ekinh1,ekinh2)
    ! PARAMETERS FOR THE NOSE-HOOVER THERMOSTATS
    IF ((cntl%tnosee.OR.cntl%tnosep).AND.paral%parent) CALL nosepa(1,1)
    ! ==--------------------------------------------------------------==
    ! == INITIALIZATION                                               ==
    ! ==--------------------------------------------------------------==
    CALL initrun(irec,c0,cm,sc0,rhoe,psi,eigv)
    CALL q_coor('GET',tau0,velp)
    ! ==--------------------------------------------------------------==
    ! INITIALIZE VELOCITIES
    IF (paral%parent) CALL detdof(tau0,taur)
    IF (irec(irec_vel).EQ.0.AND..NOT.restart1%rgeo) THEN
       CALL zeroing(cm)!,ngw*nstate)
       ener_com%ecnstr = 0.0_real_8
       ener_com%erestr = 0.0_real_8
       CALL init_velocities
       IF (cotc0%mcnstr.GT.0) THEN
          CALL q_coor('PUT',tau0,velp)
          IF (paral%parent) CALL rattle(tau0,velp)
          CALL q_coor('GET',tau0,velp)
       ENDIF
    ENDIF
    IF (cntl%quenchp) THEN
       CALL zeroing(velp)!,3*maxsys%nax*maxsys%nsx)
       CALL zeroing(clasv)!,3*clas3%nclatom)
    ENDIF
    IF (cntl%quenche) CALL zeroing(cm)!,ngw*nstate)
    ! RESET ACCUMULATORS
    IF (paral%parent.AND.irec(irec_ac).EQ.0) CALL resetac(tau0,taui,iteropt%nfi)
    IF (cntl%quenchb) THEN
       CALL quenbo(c0(:,:,1),c2,sc0,taur,rhoe,psi)
       CALL zhwwf(2,irec,c0,cm,nstate,eigv,tau0,velp,taui,iteropt%nfi)
       cntl%quenchb=.FALSE.
    ENDIF
    IF (pslo_com%tivan) THEN
       IF (cntl%tlsd) THEN
          CALL deort(ncpw%ngw,spin_mod%nsup,eigm,eigv,c0(:,1:spin_mod%nsup,1),sc0(1,1))
          CALL deort(ncpw%ngw,spin_mod%nsdown,eigm,eigv,&
               c0(:,spin_mod%nsup+1:spin_mod%nsup+spin_mod%nsdown,1),sc0(1,spin_mod%nsup+1))
       ELSE
          CALL deort(ncpw%ngw,nstate,eigm,eigv,c0,sc0)
       ENDIF
    ENDIF
    ! INITIALIZE FORCES
    IF (clas7%tfreeze) THEN
       CALL zeroing(fion)!,3*maxsys%nax*maxsys%nsx)
    ELSE
       CALL q_coor('PUT',tau0,velp)
       IF (geq0) CALL zclean(c0,nstate,ncpw%ngw)
       CALL forcedr(c0(:,:,1),c2,sc0,rhoe,psi,tau0,fion,eigv,&
            nstate,1,.FALSE.,.TRUE.)
       CALL freqs(crge%n,.TRUE.)
       ! Check orthogonality condition for wavefunction velocities
       CALL rortv(c0,cm,c2,sc0,gamy,nstate)
    ENDIF
    update=.TRUE.
    CALL clas_force(ener_com%eext,update,eint)
    ! rmdebug
    ! we save the classical force on quantum atoms to CLASFOLD
    DO ic=1,clas3%ncltyp
       is=clas3%is_qm(ic)
       IF (is.GT.0) THEN
          j=clas3%ncrang(1,ic)-1
          DO ia=1,ions0%na(is)
             clasfold(1,j+ia)=clasf(1,j+ia)
             clasfold(2,j+ia)=clasf(2,j+ia)
             clasfold(3,j+ia)=clasf(3,j+ia)
          ENDDO
       ENDIF
    ENDDO
    ! rmdebug
    CALL add_qm_force(fion)
    ! Initialize thermostats
    IF (paral%parent) THEN
       itemp=irec(irec_nop1)+irec(irec_nop2)+irec(irec_nop3)&
            +irec(irec_nop4)
       IF (cntl%tnosep .AND. itemp.EQ.0) CALL nospinit(1)
       IF (cntl%tnosee .AND. irec(irec_noe) .EQ.0) CALL noseinit(1)
       CALL wrgeof(tau0,fion)
       IF (clas7%twcoor.AND.paral%parent) CALL wr_classic
       filen='ENERGIES'
       IF (paral%io_parent)&
            CALL fileopen(3,filen,fo_app+fo_verb+fo_mark,ferror)
       filen='TEMPERATURES'
       IF (paral%io_parent)&
            CALL fileopen(15,filen,fo_app+fo_verb+fo_mark,ferror)
    ENDIF
    ! determine the number of classical and quantum atoms
    nc=0
    nq=0
    DO is=1,clas3%ncltyp
       DO ia=clas3%ncrang(1,is),clas3%ncrang(2,is)
          IF (clas3%is_qm(is).EQ.0) THEN
             nc=nc+1
          ELSE
             nq=nq+1
          ENDIF
       ENDDO
    ENDDO
    IF ((nc+nq).NE.clas3%nclatom) CALL stopgm('MDCLAS','(NC+NQ).NE.NCLATOM',& 
         __LINE__,__FILE__)
    CALL write_irec(irec)
    ! rmdebug
    extwork=0._real_8
    ! rmdebug
    ! ==--------------------------------------------------------------==
    ! == END INITIALIZATION                                           ==
    ! ==--------------------------------------------------------------==
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
       CALL mp_sync(parai%allgrp)
       time1=m_walltime()
       iteropt%nfi=iteropt%nfi+1
       ropt_mod%prteig=cprint%tprint.AND.MOD(iteropt%nfi-1,cprint%iprint_step).EQ.0
       cntl%caldip=cntl%tdipd.AND.MOD(iteropt%nfi-1,cnti%npdip).EQ.0
       comvl%subcom=comvl%tsubcom.AND.MOD(iteropt%nfi-1,comvl%ncomv).EQ.0
       comvl%subrot=comvl%tsubrot.AND.MOD(iteropt%nfi-1,comvl%nrotv).EQ.0
       IF (.NOT.paral%parent) ropt_mod%prteig=.FALSE.
       ropt_mod%engpri=cprint%tprint.AND.MOD(iteropt%nfi-1,cprint%iprint_step).EQ.0
       ! UPDATE NOSE THERMOSTATS
       IF (.NOT.clas7%tfreeze) THEN
          CALL q_coor('PUT',tau0,velp)
          CALL anneal(velp,cm,nstate,scr)
          CALL berendsen(velp,cm,nstate,scr,ekinc,0.0_real_8)
          ! SUBTRACT CENTER OF MASS VELOCITY
          IF (comvl%subcom) CALL comvel(velp,vcmio,.TRUE.)
          ! SUBTRACT ROTATION AROUND CENTER OF MASS
          IF (comvl%subrot) CALL rotvel(tau0,velp,lmio,tauio,.TRUE.)
          ! UPDATE NOSE THERMOSTATS
          CALL noseup(velp,cm,nstate,1)
          IF (cntl%annei) THEN
             alfap=cntr%anneri**(0.25_real_8)
             DO is=1,clas3%ncltyp
                DO ia=clas3%ncrang(1,is),clas3%ncrang(2,is)
                   clasv(1,ia)=alfap*clasv(1,ia)
                   clasv(2,ia)=alfap*clasv(2,ia)
                   clasv(3,ia)=alfap*clasv(3,ia)
                ENDDO
             ENDDO
          ENDIF
          CALL q_coor('GET',tau0,velp)
       ENDIF
       ! UPDATE VELOCITIES
       DO is=1,clas3%ncltyp
          IF (.NOT.clas7%tfreeze.OR.clas3%is_qm(is).EQ.0) THEN
             DO ia=clas3%ncrang(1,is),clas3%ncrang(2,is)
                clasv(1,ia)=clasv(1,ia)+clas4%cstep(is)*clasf(1,ia)
                clasv(2,ia)=clasv(2,ia)+clas4%cstep(is)*clasf(2,ia)
                clasv(3,ia)=clasv(3,ia)+clas4%cstep(is)*clasf(3,ia)
             ENDDO
          ENDIF
       ENDDO
       IF (.NOT.clas7%tfreeze) CALL velupa(c0,cm,c2,nstate,1)
       ! UPDATE POSITIONS
       CALL q_coor('PUT',tau0,velp)
       ener_com%ecnstr = 0.0_real_8
       DO is=1,clas3%ncltyp
          IF (.NOT.clas7%tfreeze.OR.clas3%is_qm(is).EQ.0) THEN
             DO ia=clas3%ncrang(1,is),clas3%ncrang(2,is)
                clasc(1,ia)=clasc(1,ia)+dt_ions*clasv(1,ia)
                clasc(2,ia)=clasc(2,ia)+dt_ions*clasv(2,ia)
                clasc(3,ia)=clasc(3,ia)+dt_ions*clasv(3,ia)
                ! 
                delclasc(1,ia)=dt_ions*clasv(1,ia)
                delclasc(2,ia)=dt_ions*clasv(2,ia)
                delclasc(3,ia)=dt_ions*clasv(3,ia)
                ! 
             ENDDO
          ENDIF
       ENDDO
       ropt_mod%calstc=cntl%tprec.AND.MOD(iteropt%nfi-1,cnti%nprec).EQ.0
       IF (clas7%tfreeze) THEN
          CALL zeroing(fion)!,3*maxsys%nax*maxsys%nsx)
       ELSE
          CALL q_coor('PUT',taup,velp)
          IF (cotc0%mcnstr.NE.0) THEN
             IF (paral%parent) CALL cpmdshake(tau0,taup,velp)
             CALL q_coor('GET',taup,velp)
          ENDIF
          CALL phfac(taup)
          IF (corel%tinlc) CALL copot(rhoe,psi,ropt_mod%calste)
          CALL posupa(c0,cm,c2,gamx,nstate)
          ! ..Dipole moment
          IF (cntl%caldip) THEN
             CALL ddipo(taup,c0(:,:,1),cm,c2,sc0,nstate,center)
             CALL wannier_print(iteropt%nfi,c0(:,:,1),taup,nstate,psi(:,1),center)
          ENDIF
          ! CALCULATE THE FORCES
          ropt_mod%calste=cntl%tpres.AND.MOD(iteropt%nfi-1,cnti%npres).EQ.0
          CALL forcedr(c0(:,:,1),c2,sc0,rhoe,psi,taup,fion,eigv,&
               nstate,1,.FALSE.,.TRUE.)
          IF (ropt_mod%calste) CALL totstr
       ENDIF
       ! we have to check whether the Verlet list need to be updated
       CALL check(update)
       ! 
       ! we calculate a contribution to external work on the cluster
       DO ic=1,clas3%ncltyp
          is=clas3%is_qm(ic)
          IF (is.GT.0) THEN
             j=clas3%ncrang(1,ic)-1
             DO ia=1,ions0%na(is)
                extwork=extwork+clasfold(1,j+ia)*delclasc(1,j+ia)/2._real_8
                extwork=extwork+clasfold(2,j+ia)*delclasc(2,j+ia)/2._real_8
                extwork=extwork+clasfold(3,j+ia)*delclasc(3,j+ia)/2._real_8
             ENDDO
          ENDIF
       ENDDO
       CALL clas_force(ener_com%eext,update,eint)
       ! 
       ! we calculate a contribution to external work on the cluster
       DO ic=1,clas3%ncltyp
          is=clas3%is_qm(ic)
          IF (is.GT.0) THEN
             j=clas3%ncrang(1,ic)-1
             DO ia=1,ions0%na(is)
                extwork=extwork+clasf(1,j+ia)*delclasc(1,j+ia)/2._real_8
                extwork=extwork+clasf(2,j+ia)*delclasc(2,j+ia)/2._real_8
                extwork=extwork+clasf(3,j+ia)*delclasc(3,j+ia)/2._real_8
             ENDDO
          ENDIF
       ENDDO
       ! WRITE(17,*)EXTWORK
       ! we save the classical force on quantum atoms to CLASFOLD
       DO ic=1,clas3%ncltyp
          is=clas3%is_qm(ic)
          IF (is.GT.0) THEN
             j=clas3%ncrang(1,ic)-1
             DO ia=1,ions0%na(is)
                clasfold(1,j+ia)=clasf(1,j+ia)
                clasfold(2,j+ia)=clasf(2,j+ia)
                clasfold(3,j+ia)=clasf(3,j+ia)
             ENDDO
          ENDIF
       ENDDO
       CALL add_qm_force(fion)
       ! ==================================================================
       ! Damped Dynamics
       CALL dampdyn(velp,fion,cm,c2,nstate,scr(1),scr(10))
       ! ==================================================================
       ! FINAL UPDATE FOR VELOCITIES
       DO is=1,clas3%ncltyp
          IF (.NOT.clas7%tfreeze.OR.clas3%is_qm(is).EQ.0) THEN
             DO ia=clas3%ncrang(1,is),clas3%ncrang(2,is)
                clasv(1,ia)=clasv(1,ia)+clas4%cstep(is)*clasf(1,ia)
                clasv(2,ia)=clasv(2,ia)+clas4%cstep(is)*clasf(2,ia)
                clasv(3,ia)=clasv(3,ia)+clas4%cstep(is)*clasf(3,ia)
             ENDDO
          ENDIF
       ENDDO
       IF (.NOT.clas7%tfreeze) THEN
          CALL velupa(c0,cm,c2,nstate,1)
          CALL q_coor('PUT',tau0,velp)
          IF (cotc0%mcnstr.NE.0) THEN
             CALL rattle(tau0,velp)
             CALL q_coor('GET',tau0,velp)
          ENDIF
          CALL rortv(c0,cm,c2,sc0,gamy,nstate)
       ENDIF
       IF (paral%parent) CALL geofile(tau0,velp,'WRITE')
       ! COMPUTE THE IONIC TEMPERATURE TEMPP
       CALL ekincl(ekinp,ekincm,ekinqm)
       tempp=ekinp*factem*2._real_8/(3*clas3%nclatom-3)
       ! IONIC TEMPERATURE CONTROL
       CALL temp_control(temp1,temp2,tempp)
       IF (clas7%tdefo.AND.MOD(iteropt%nfi,ndefo).EQ.0) THEN
          CALL deformation
       ENDIF
       ! UPDATE NOSE THERMOSTATS
       IF (cntl%tnosee.OR.cntl%tc) CALL rekine(cm,nstate,ekinc)
       IF (.NOT.clas7%tfreeze) THEN
          CALL q_coor('PUT',tau0,velp)
          CALL noseup(velp,cm,nstate,1)
          ! SUBTRACT ROTATION AROUND CENTER OF MASS
          IF (comvl%subrot) CALL rotvel(tau0,velp,lmio,tauio,.FALSE.)
          ! SUBTRACT CENTER OF MASS VELOCITY
          IF (comvl%subcom) CALL comvel(velp,vcmio,.FALSE.)
          CALL anneal(velp,cm,nstate,scr)
          CALL berendsen(velp,cm,nstate,scr,ekinc,0.0_real_8)
          IF (cntl%annei) THEN
             alfap=cntr%anneri**(0.25_real_8)
             DO is=1,clas3%ncltyp
                DO ia=clas3%ncrang(1,is),clas3%ncrang(2,is)
                   clasv(1,ia)=alfap*clasv(1,ia)
                   clasv(2,ia)=alfap*clasv(2,ia)
                   clasv(3,ia)=alfap*clasv(3,ia)
                ENDDO
             ENDDO
          ENDIF
          CALL q_coor('GET',tau0,velp)
          CALL ekincl(ekinp,ekincm,ekinqm)
          tempp=ekinp*factem*2._real_8/(3*clas3%nclatom-3)
          tempcm=ekincm*factem*2._real_8/(3*nc)
          tempqm=ekinqm*factem*2._real_8/(3*nq)
          ! RESCALE ELECTRONIC VELOCITIES
          IF (cntl%tc) CALL rscve(ekin1,ekin2,ekinc,cntr%ekinw,cm,nstate,ncpw%ngw)
          ! KINETIC ENERGY OF THE ELECTRONS
          CALL rekine(cm,nstate,ekinc)
          ! ENERGY OF THE NOSE THERMOSTATS
          IF (paral%parent) CALL noseng(iteropt%nfi,velp,enose,enosp,dummy,1)
       ELSE
          enose=0._real_8
          enosp=0._real_8
          ekinc=0._real_8
          CALL ekincl(ekinp,ekincm,ekinqm)
          tempp=ekinp*factem*2._real_8/(3*clas3%nclatom-3)
          tempcm=ekincm*factem*2._real_8/(3*nc)
          tempqm=ekinqm*factem*2._real_8/(3*nq)
       ENDIF
       ! MEAN SQUARE DISPLACEMENT OF DIFFERENT IONIC SPECIES
       IF (paral%parent) CALL dispp(tau0,taui,disa)
       IF (paral%parent) THEN
          econs=ekinp+ener_com%etot+enose+enosp+ener_com%ecnstr+ener_com%eext
          eham=econs+ekinc
          time2=m_walltime()
          tcpu=(time2-time1)*0.001_real_8
          ! PRINTOUT the evolution of the accumulators every time step
          IF (paral%io_parent)&
               WRITE(16,789)ekinqm,ekinqm+ener_com%etot,ener_com%eext&
               ,eint,extwork
789       FORMAT(5(f12.6,2x))
          IF (clas7%tfreeze) ropt_mod%engpri=.FALSE.
          CALL wrprint_md(eigv,crge%f,ener_com%amu,nstate,tau0,fion,&
               ekinc,tempp,ener_com%etot,econs,eham,disa,&
               tcpu,.FALSE.,iteropt%nfi,infi)
          ! UPDATE ACCUMULATORS
          CALL pacca(ekinc,tempp,ener_com%etot,econs,eham,enose,enosp,ener_com%ecnstr,&
               ener_com%erestr,disa,tcpu,iteropt%nfi,1)
          ! Store ionic coordinates and velocities for statistics
          ! MOVIE=MOD(NFI-1,IMOVIE).EQ.0
          ropt_mod%movie=clas7%tmovc.AND.MOD(iteropt%nfi-1,imovc).EQ.0
          ropt_mod%rprint=rout1%rout .AND. MOD(iteropt%nfi-1,cnti%ntraj).EQ.0
          ropt_mod%txyz=rout1%xtout .AND. MOD(iteropt%nfi-1,cnti%ntraj).EQ.0
          ropt_mod%tdcd=rout1%dcout .AND. MOD(iteropt%nfi-1,cnti%ntraj).EQ.0
          CALL printp(taur,tau0,velp)
          ! 
          clas7%cprint=clas7%tftra.AND.MOD(iteropt%nfi-1,ntrac).EQ.0
          IF (clas7%tftra.OR.clas7%tmovc) CALL wr_cltra
          ! IF(CPRINT)CALL WR_CLTRA
          ! 
       ENDIF
       IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)
       IF (infi.EQ.cnti%nomore) THEN
          soft_com%exsoft=.TRUE.
          soft_com%exnomore=.TRUE.
       ENDIF
       IF (MOD(iteropt%nfi,store1%istore).EQ.0.OR.infi.EQ.cnti%nomore.OR.soft_com%exsoft)&
            CALL zhwwf(2,irec,c0,cm,nstate,eigv,tau0,velp,taui,iteropt%nfi)
       ! temperature ramping
       CALL tempramp(temp1,temp2)
       ! STOP THE RUN IF THE USER HAS SET THE SIGNAL 30
       IF (soft_com%exsoft) GOTO 100
       ! ==================================================================
       ! ==     END OF MAIN LOOP                                         ==
       ! ==================================================================
    ENDDO
    DEALLOCATE(center,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
100 CONTINUE
    IF (rout1%rhoout.AND. .NOT.clas7%tfreeze)&
         CALL rhopri(c0,tau0,rhoe,psi(:,1),nstate,nkpt%nkpnt)
    ! PRINT ACCUMULATOR
    IF (paral%parent) CALL pacca(ekinc,tempp,ener_com%etot,econs,eham,enose,enosp,&
         ener_com%ecnstr,ener_com%erestr,disa,tcpu,iteropt%nfi,2)
    IF (.NOT.clas7%tfreeze) THEN
       CALL proja(c0,c2,sc0,nstate,cnti%iproj)
       CALL csize(c2,crge%n,gemax,cnorm)
    ENDIF
    IF (paral%parent) CALL gsize(fion,gnmax,gnorm)
    IF (paral%parent) CALL finalp(tau0,fion,velp,eigv)
    ! 
    IF (clas7%twcoor.AND.paral%parent)CALL wr_classic
    ! 
    IF (pslo_com%tivan.OR.cntl%nonort) THEN
       DEALLOCATE(eigv,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(eigm,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            CALL fileclose(3)
       IF (paral%io_parent)&
            CALL fileclose(15)
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
    IF (comvl%tsubrot) DEALLOCATE(tauio,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE mdclas
  ! ==================================================================
  SUBROUTINE q_coor(tag,tau0,velp)
    ! ==--------------------------------------------------------------==
    CHARACTER(len=3)                         :: tag
    REAL(real_8)                             :: tau0(:,:,:), velp(:,:,:)

    INTEGER                                  :: ia, ic, is, j

    CALL mp_bcast(tau0,SIZE(tau0),parai%source,parai%allgrp)
    CALL mp_bcast(velp,SIZE(velp),parai%source,parai%allgrp)
    IF (tag.EQ.'PUT') THEN
       DO ic=1,clas3%ncltyp
          is=clas3%is_qm(ic)
          IF (is.GT.0) THEN
             j=clas3%ncrang(1,ic)-1
             DO ia=1,ions0%na(is)
                tau0(1,ia,is)=clasc(1,j+ia)
                tau0(2,ia,is)=clasc(2,j+ia)
                tau0(3,ia,is)=clasc(3,j+ia)
                velp(1,ia,is)=clasv(1,j+ia)
                velp(2,ia,is)=clasv(2,j+ia)
                velp(3,ia,is)=clasv(3,j+ia)
             ENDDO
          ENDIF
       ENDDO
    ELSEIF (tag.EQ.'GET') THEN
       DO ic=1,clas3%ncltyp
          is=clas3%is_qm(ic)
          IF (is.GT.0) THEN
             j=clas3%ncrang(1,ic)-1
             DO ia=1,ions0%na(is)
                clasc(1,j+ia)=tau0(1,ia,is)
                clasc(2,j+ia)=tau0(2,ia,is)
                clasc(3,j+ia)=tau0(3,ia,is)
                clasv(1,j+ia)=velp(1,ia,is)
                clasv(2,j+ia)=velp(2,ia,is)
                clasv(3,j+ia)=velp(3,ia,is)
             ENDDO
          ENDIF
       ENDDO
    ELSE
       CALL stopgm('Q_COOR','WRONG TAG',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE q_coor
  ! ==================================================================
  SUBROUTINE add_qm_force(fion)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: fion(:,:,:)

    INTEGER                                  :: ia, ic, is, j

    CALL mp_bcast(fion,SIZE(fion),parai%source,parai%allgrp)
    DO ic=1,clas3%ncltyp
       is=clas3%is_qm(ic)
       IF (is.GT.0) THEN
          j=clas3%ncrang(1,ic)-1
          DO ia=1,ions0%na(is)
             clasf(1,j+ia)=clasf(1,j+ia)+fion(1,ia,is)
             clasf(2,j+ia)=clasf(2,j+ia)+fion(2,ia,is)
             clasf(3,j+ia)=clasf(3,j+ia)+fion(3,ia,is)
          ENDDO
       ENDIF
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE add_qm_force
  ! ==================================================================
  SUBROUTINE init_velocities
    ! ==--------------------------------------------------------------==
    ! Variables
    INTEGER                                  :: ia, is
    REAL(real_8)                             :: alfa1, alfa2, alfa3, ekincm, &
                                                ekinp, ekinqm, pma00, rnr(3), &
                                                sigma, tempp, tscal, vcm(3), &
                                                vscale

! ==--------------------------------------------------------------==
! Initialize Cell velocities

    CALL zeroing(metr_com%htvel)!,9)
    ! MAXWELL DISTRIBUTION FOR THE IONS
    IF (paral%parent) THEN
       IF (cntr%tempw.LT.1.e-5_real_8) GOTO 100
       CALL repprngu_vec(3,rnr)
       DO is=1,clas3%ncltyp
          sigma=SQRT(cntr%tempw/(clas4%clmas(is)*factem*scmass))
          IF (.NOT.clas7%tfreeze.OR.clas3%is_qm(is).EQ.0) THEN
             DO ia=clas3%ncrang(1,is),clas3%ncrang(2,is)
                CALL repprngu_vec(3,rnr)
                alfa1=2.0_real_8*pi*rnr(1)
                alfa2=2.0_real_8*pi*rnr(2)
                alfa3=2.0_real_8*pi*rnr(3)
                CALL repprngu_vec(3,rnr)
                clasv(1,ia)=SQRT(LOG(rnr(1))*(-2.0_real_8))*COS(alfa1)*sigma
                clasv(2,ia)=SQRT(LOG(rnr(2))*(-2.0_real_8))*COS(alfa2)*sigma
                clasv(3,ia)=SQRT(LOG(rnr(3))*(-2.0_real_8))*COS(alfa3)*sigma
             ENDDO
          ENDIF
       ENDDO
       ! SUBTRACT CENTER OF MASS VELOCITY
       vcm(1)=0.0_real_8
       vcm(2)=0.0_real_8
       vcm(3)=0.0_real_8
       pma00=0._real_8
       DO is=1,clas3%ncltyp
          IF (.NOT.clas7%tfreeze.OR.clas3%is_qm(is).EQ.0) THEN
             DO ia=clas3%ncrang(1,is),clas3%ncrang(2,is)
                vcm(1)=vcm(1)+clasv(1,ia)*clas4%clmas(is)
                vcm(2)=vcm(2)+clasv(2,ia)*clas4%clmas(is)
                vcm(3)=vcm(3)+clasv(3,ia)*clas4%clmas(is)
                pma00=pma00+clas4%clmas(is)
             ENDDO
          ENDIF
       ENDDO
       DO is=1,clas3%ncltyp
          IF (.NOT.clas7%tfreeze.OR.clas3%is_qm(is).EQ.0) THEN
             DO ia=clas3%ncrang(1,is),clas3%ncrang(2,is)
                clasv(1,ia)=clasv(1,ia)-vcm(1)/pma00
                clasv(2,ia)=clasv(2,ia)-vcm(2)/pma00
                clasv(3,ia)=clasv(3,ia)-vcm(3)/pma00
             ENDDO
          ENDIF
       ENDDO
       ! RESCALE VELOCITIES
       CALL ekincl(ekinp,ekincm,ekinqm)
       tempp=ekinp*factem*2._real_8/REAL(3*clas3%nclatom-3,kind=real_8)
       IF (tempp.GT.1.e-5_real_8) THEN
          tscal=cntr%tempw/tempp
          vscale=SQRT(tscal)
          DO ia=1,clas3%nclatom
             clasv(1,ia)=clasv(1,ia)*vscale
             clasv(2,ia)=clasv(2,ia)*vscale
             clasv(3,ia)=clasv(3,ia)*vscale
          ENDDO
       ENDIF
100    CONTINUE
    ENDIF
    CALL mp_bcast(clasv,SIZE(clasv),parai%source,parai%allgrp)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE init_velocities
  ! ==================================================================
  SUBROUTINE ekincl(ekinp,ekincm,ekinqm)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: ekinp, ekincm, ekinqm

    INTEGER                                  :: ia, is
    REAL(real_8)                             :: const, ekin

! Variables
! ==--------------------------------------------------------------==
! ==  CALCULATE KINETIC ENERGY OF THE IONS                        ==
! ==--------------------------------------------------------------==

    ekinp=0._real_8
    ekincm=0._real_8
    ekinqm=0._real_8
    DO is=1,clas3%ncltyp
       const=0.5_real_8*clas4%clmas(is)*scmass
       ekin=0._real_8
       DO ia=clas3%ncrang(1,is),clas3%ncrang(2,is)
          ekin=ekin+const*clasv(1,ia)*clasv(1,ia)
          ekin=ekin+const*clasv(2,ia)*clasv(2,ia)
          ekin=ekin+const*clasv(3,ia)*clasv(3,ia)
       ENDDO
       IF (clas3%is_qm(is).EQ.0) THEN
          ekincm=ekincm+ekin
       ELSE
          ekinqm=ekinqm+ekin
       ENDIF
    ENDDO
    ekinp=ekincm+ekinqm
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE ekincl
  ! ==================================================================
  SUBROUTINE deformation
    ! ==--------------------------------------------------------------==
    ! Variables
    INTEGER                                  :: ia, ic, is, j

! ==--------------------------------------------------------------==

    DO ic=1,clas3%ncltyp
       is=clas3%is_qm(ic)
       IF (is.GT.0) THEN
          j=clas3%ncrang(1,ic)-1
          DO ia=1,ions0%na(is)
             clasc(1,j+ia)=defy%xcen+defy%expx*(clasc(1,j+ia)-defy%xcen)
             clasc(2,j+ia)=defy%ycen+defy%expy*(clasc(2,j+ia)-defy%ycen)
             clasc(3,j+ia)=defy%zcen+defy%expz*(clasc(3,j+ia)-defy%zcen)
          ENDDO
       ENDIF
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE deformation
  ! ==================================================================
  SUBROUTINE temp_control(temp1,temp2,tempp)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: temp1, temp2, tempp

    INTEGER                                  :: ia, is
    REAL(real_8)                             :: alfap

! Variables
! ==--------------------------------------------------------------==
! ==  Dynamical rescaling factor (tempw/tempp), where tempp is    ==
! ==  calculated every step                                       ==
! ==--------------------------------------------------------------==

    IF (cntl%tcp) THEN
       IF (tempp.GT.temp1.OR.tempp.LT.temp2.AND.tempp.NE.0._real_8) THEN
          alfap=SQRT(cntr%tempw/tempp)
          DO is=1,clas3%ncltyp
             DO ia=clas3%ncrang(1,is),clas3%ncrang(2,is)
                clasv(1,ia)=alfap*clasv(1,ia)
                clasv(2,ia)=alfap*clasv(2,ia)
                clasv(3,ia)=alfap*clasv(3,ia)
             ENDDO
          ENDDO
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE temp_control
  ! ==================================================================
  SUBROUTINE wr_classic
    ! ==--------------------------------------------------------------==
    ! Variables
    INTEGER                                  :: ia, is, k

! ==--------------------------------------------------------------==

    IF (paral%io_parent)&
         WRITE(6,'(/,1X,64("*"))')
    IF (paral%io_parent)&
         WRITE(6,'(22X,A)') ' ATOMIC COORDINATES'
    DO is=1,clas3%ncltyp
       DO ia=clas3%ncrang(1,is),clas3%ncrang(2,is)
          IF (paral%io_parent)&
               WRITE(6,'(I5,3X,A4,3X,I5,3F15.6)')&
               ia,clab(is),clas3%is_qm(is),(clasc(k,ia),k=1,3)
       ENDDO
    ENDDO
    IF (paral%io_parent)&
         WRITE(6,'(1X,64("*"),/)')
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE wr_classic
  ! ==================================================================
  SUBROUTINE wr_cltra
    ! ==--------------------------------------------------------------==
    CHARACTER(len=100)                       :: fnmov, fntrj
    INTEGER                                  :: i, ia, iat, is, iz0
    INTEGER, SAVE                            :: ifim = 0, ifit = 0
    LOGICAL                                  :: ferror
    REAL(real_4)                             :: ccx, ccy, ccz, vcx, vcy, vcz

    fntrj='CLASSIC'
    ! mdebug classical trajectory written independently from 
    ! mdebug QM trajectory (===> better for freeze quantum cases) 
    IF (clas7%cprint) THEN
       IF (ifit.EQ.0) THEN
          IF (paral%io_parent)&
               CALL fileopen(4,fntrj,fo_app+fo_ufo+fo_verb,ferror)
          ifit=1
       ELSE
          IF (paral%io_parent)&
               CALL fileopen(4,fntrj,fo_app+fo_ufo,ferror)
       ENDIF
       DO is=1,clas3%ncltyp
          DO ia=clas3%ncrang(1,is),clas3%ncrang(2,is)
             ! mdebug            WRITE(4,'(I7,6(2X,F14.6))')
             ccx=clasc(1,ia)
             ccy=clasc(2,ia)
             ccz=clasc(3,ia)
             vcx=clasv(1,ia)
             vcy=clasv(2,ia)
             vcz=clasv(3,ia)
             ! mdebug     *                NFI,(CLASC(I,IA),I=1,3),(CLASV(I,IA),I=1,3)
             IF (paral%io_parent)&
                  WRITE(4)&
                  ccx,ccy,ccz,vcx,vcy,vcz
          ENDDO
       ENDDO
       IF (paral%io_parent)&
            CALL fileclose(4)
       ! mdebug
    ENDIF
    ! mdebug
    ! 
    fnmov='CLASSIC_MOVIE'
    IF (ropt_mod%movie)THEN
       IF (ifim.EQ.0) THEN
          IF (paral%io_parent)&
               CALL fileopen(4,fnmov,fo_app+fo_verb+fo_mark,ferror)
          ifim=1
       ELSE
          IF (paral%io_parent)&
               CALL fileopen(4,fnmov,fo_app,ferror)
       ENDIF
       DO is=1,clas3%ncltyp
          iz0=imtyp(is)
          IF (ions0%iatyp(is).EQ.0) THEN
             iat=6
          ELSE
             iat=ions0%iatyp(is)
          ENDIF
          DO ia=clas3%ncrang(1,is),clas3%ncrang(2,is)
             ! mdebug wrong movie format... 
             ! WRITE(4,'(I7,6(2X,F14.6))')
             ! mdebug
             IF (paral%io_parent)&
                  WRITE(4,'(3(2X,F12.4),2I4)')&
                  (clasc(i,ia)/fbohr,i=1,3),iat,iz0
          ENDDO
       ENDDO
       IF (paral%io_parent)&
            CALL fileclose(4)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE wr_cltra
  ! ==================================================================

END MODULE mdclas_utils
