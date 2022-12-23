MODULE mdshop_cp_utils
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
  USE cotr,                            ONLY: cotc0
  USE csize_utils,                     ONLY: csize
  USE ddipo_utils,                     ONLY: ddipo,&
                                             give_scr_ddipo
  USE deort_utils,                     ONLY: deort,&
                                             give_scr_deort
  USE detdof_utils,                    ONLY: detdof
  USE dispp_utils,                     ONLY: dispp
  USE dynit_utils,                     ONLY: dynit
  USE ekinpp_utils,                    ONLY: ekinpp
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_app,&
                                             fo_info,&
                                             fo_verb
  USE finalp_utils,                    ONLY: finalp
  USE forcedr_driver,                  ONLY: forcedr
  USE forcedr_utils,                   ONLY: give_scr_forcedr
  USE forcep_utils,                    ONLY: rhoe_psi_size
  USE geofile_utils,                   ONLY: geofile
  USE geq0mod,                         ONLY: geq0
  USE gsize_utils,                     ONLY: gsize
  USE hubbardu,                        ONLY: hubbu
  USE initrun_driver,                  ONLY: initrun
  USE initrun_utils,                   ONLY: give_scr_initrun
  USE kinds,                           ONLY: real_8
  USE machine,                         ONLY: m_walltime
  USE mp_interface,                    ONLY: mp_bcast,&
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
  USE printp_utils,                    ONLY: printp
  USE proja_utils,                     ONLY: proja
  USE pslo,                            ONLY: pslo_com
  USE puttau_utils,                    ONLY: taucl
  USE quenbo_utils,                    ONLY: give_scr_quenbo,&
                                             quenbo
  USE rattle_utils,                    ONLY: rattle
  USE rekine_utils,                    ONLY: rekine
  USE resetac_utils,                   ONLY: resetac
  USE rhopri_utils,                    ONLY: give_scr_rhopri
  USE rinvel_utils,                    ONLY: rinvel,&
                                             rvscal
  USE rk4ov_utils,                     ONLY: rk4ov_new,&
                                             rk4ov_old
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
  USE setirec_utils,                   ONLY: write_irec
  USE shake_utils,                     ONLY: cpmdshake
  USE shop,                            ONLY: sh02,&
                                             tlsd0,&
                                             tshopold
  USE shop_adds_utils,                 ONLY: decide,&
                                             s0_s1_overlap,&
                                             state_select,&
                                             write_shmd
  USE shop_rest,                       ONLY: prob1,&
                                             sh03
  USE soft,                            ONLY: soft_com
  USE spin,                            ONLY: clsd
  USE store_types,                     ONLY: &
       cprint, irec_ac, irec_noe, irec_nop1, irec_nop2, irec_nop3, irec_nop4, &
       irec_vel, restart1, rout1
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             cntr,&
                                             maxsys,&
                                             ncpw,&
                                             restf
  USE testex_utils,                    ONLY: testex
  USE teststore_utils,                 ONLY: teststore
  USE totstr_utils,                    ONLY: totstr
  USE utils,                           ONLY: zclean
  USE velupa_utils,                    ONLY: velupa
  USE velupi_utils,                    ONLY: velupi
  USE wann,                            ONLY: wannl
  USE wannier_print_utils,             ONLY: wannier_print
  USE wc_dos_utils,                    ONLY: wc_dos
  USE wrener_utils,                    ONLY: wrprint_md
  USE wrgeo_utils,                     ONLY: wrgeof
  USE wv30_utils,                      ONLY: zhwwf
  USE zeroing_utils,                   ONLY: zeroing
#include "sizeof.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mdshop_cp
  PUBLIC :: give_scr_mdshop_cp

CONTAINS

  ! ==================================================================
  SUBROUTINE mdshop_cp(c0,cm,c2,sc0,gamx,gamy)
    ! ==--------------------------------------------------------------==
    ! McB ... surface hopping stuff ...
    ! McB ............................................................
    COMPLEX(real_8) :: c0(ncpw%ngw,crge%n,1), cm(ncpw%ngw,crge%n), &
      c2(ncpw%ngw,crge%n), sc0(ncpw%ngw,crge%n)
    REAL(real_8)                             :: gamx(*), gamy(*)

    CHARACTER(len=100)                       :: filen
    CHARACTER(len=30)                        :: tag
    COMPLEX(real_8), ALLOCATABLE             :: psi(:,:), psi1(:,:)
    INTEGER :: i, ierr, il_psi_1d, il_psi_2d, il_rhoe_1d, il_rhoe_2d, &
      irec(100), itemp, j, lscr, ncoef, ns1, nstate, ntmp
    LOGICAL                                  :: ferror
    REAL(real_8) :: c12, c21, dcoupl(2,2), ddet, dei(2), delt, detold2, disa, &
      dtcoef, dummy, e(2), econs, eham, ei(sh02%nsurf), ekin1, ekin2, ekinc, &
      ekincp, ekinh1, ekinh2, ekinp, enose, enosp, lmio(3), tcpu, temp1, &
      temp2, tempp, time1, time2, vcmio(4)
    REAL(real_8), ALLOCATABLE :: center(:,:), eigm(:), eigv(:), eigv1(:), &
      fion0(:,:,:), fion1(:,:,:), rhoe(:,:), rhoe1(:,:), scr(:), taui(:,:,:), &
      tauio(:,:), taur(:,:,:)

! Variables
! .....BORN-OPPENHEIMER..................................................
! .....END BORN-OPPENHEIMER..............................................
! McB ... surface hopping stuff ...

    COMMON /ekin/ekinp
    LOGICAL :: lexist
    ! McB ............................................................
    CHARACTER(*),PARAMETER::procedureN='mdshop_cp'
    ! ==================================================================
    time1 =m_walltime()
    ALLOCATE(taup(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(taui(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(taur(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(fion(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(taup)!,SIZE(taup))
    CALL zeroing(fion)!,SIZE(fion))
    ALLOCATE(fion0(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(fion1(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    IF (comvl%tsubrot)  THEN
       ALLOCATE(tauio(3,maxsys%nax),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)! TODO check dimensions
    ENDIF

    CALL zeroing(fion)!,3*maxsys%nax*maxsys%nsx)
    CALL zeroing(taui)!,3*maxsys%nax*maxsys%nsx)
    CALL zeroing(taup)!,3*maxsys%nax*maxsys%nsx)

    CALL zeroing(sh03%coupl)!,2*2)
    CALL zeroing(sh03%couplold)!,2*2)
    CALL zeroing(fion0)!,3*maxsys%nax*maxsys%nsx)
    CALL zeroing(fion1)!,3*maxsys%nax*maxsys%nsx)
    DO i=1,2
       e(i)=0.0_real_8
       sh03%ec(i)=0.0_real_8
       sh03%eold(i)=0.0_real_8
    ENDDO
    sh03%det=0.0_real_8
    sh03%detold=0.0_real_8
    detold2=0.0_real_8
    sh03%isurf=0
    delt=cntr%delt_ions
    ncoef=INT(delt/0.04_real_8)+1
    dtcoef=delt/REAL(ncoef,kind=real_8)
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,*)'NCOEF=',ncoef,'DTCOEF=',dtcoef
       IF (paral%io_parent)&
            WRITE(6,*)'DT=',delt,'DTCOEF*NCOEF=',dtcoef*ncoef
    ENDIF
    tlsd0=cntl%tlsd
    crge%n=sh02%nst_s0+sh02%nst_s1
    nstate=crge%n
    clsd%nlsd=4
    clsd%nlsx=3
    sh02%nsttot=sh02%nst_s0 + sh02%nst_s1
    ns1=sh02%nst_s0+1
    iteropt%nfi=0
    ropt_mod%modens=.FALSE.
    ropt_mod%engpri=.FALSE.
    ropt_mod%calste=cntl%tpres
    ALLOCATE(eigv(nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(eigv1(nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    IF (pslo_com%tivan) cntl%nonort=.TRUE.
    IF (cntl%nonort) THEN
       ! CALL MEMORY(IP_EIGV,NSTATE,'EIGV')
       ALLOCATE(eigm(nstate*nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    IF (cntl%tharm.AND.paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,'(A,A)') ' HARMONIC REFERENCE SYSTEM INTEGRATION',&
            ' ONLY POSSIBLE WITH EQUAL OCCUPATION NUMBERS'
       CALL stopgm('MDSHOP_CP',' ',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! SCR ALLOCATION AND PARTITION (SCRATCH ARRAY).
    CALL rhoe_psi_size(il_psi_1d=il_psi_1d,il_psi_2d=il_psi_2d,&
         il_rhoe_1d=il_rhoe_1d,il_rhoe_2d=il_rhoe_2d)
    ALLOCATE(rhoe(il_rhoe_1d,il_rhoe_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(rhoe1(il_rhoe_1d,il_rhoe_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(psi(il_psi_1d,il_psi_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(psi1(il_psi_1d,il_psi_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL give_scr_mdshop_cp(lscr,tag)
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
    ! PARAMETERS FOR THE NOSE-HOOVER THERMOSTATS
    IF ((cntl%tnosee.OR.cntl%tnosep).AND.paral%parent) CALL nosepa(1,1)
    ! Dont symmetrize density 
    cntl%tsymrho=.FALSE.
    ! ==--------------------------------------------------------------==
    ! == INITIALISATION                                               ==
    ! ==--------------------------------------------------------------==
    ! read RESTART file:      
    CALL initrun(irec,c0,cm,sc0,rhoe,psi,eigv)
    CALL mp_bcast(taup,SIZE(taup),parai%source,parai%allgrp)
    CALL dcopy(3*maxsys%nax*maxsys%nsx,taup,1,taui,1)
    CALL phfac(tau0)
    IF (cntl%quenchb) THEN
       CALL state_select("S0")
       CALL quenbo(c0(:,1:ns1-1,1),c2,sc0,tau0,rhoe,psi) !vw added tau0
       CALL state_select("S1")
       CALL quenbo(c0(:,ns1:,1),c2,sc0,tau0,rhoe,psi) !vw added tau0
       ! McB... cf. elct.inc
       ntmp=crge%n
       crge%n=sh02%nsttot
       CALL zhwwf(2,irec,c0,cm,sh02%nsttot,eigv,tau0,velp,taui,iteropt%nfi)
       crge%n=ntmp
    ENDIF
    IF (pslo_com%tivan) THEN
       CALL deort(ncpw%ngw,sh02%nst_s0,eigm,eigv,c0(1,1,1),sc0)
       CALL deort(ncpw%ngw,sh02%nst_s1,eigm,eigv,c0(1,ns1,1),sc0)
    ENDIF
    ! INITIALIZE VELOCITIES
    IF (paral%parent) CALL detdof(tau0,taur)
    IF (irec(irec_vel).EQ.0.AND..NOT.restart1%rgeo) THEN
       ener_com%ecnstr = 0.0_real_8
       ener_com%erestr = 0.0_real_8
       CALL rinvel(velp,cm,sh02%nsttot)
       IF (paral%parent) CALL taucl(velp)
       IF (paral%parent) CALL rattle(tau0,velp)
       CALL rvscal(velp)
    ELSE
       IF (paral%parent) CALL taucl(velp)
       IF (paral%parent) CALL rattle(tau0,velp)
       IF (cntl%trescale) CALL rvscal(velp)
    ENDIF
    IF (cntl%quenchp) CALL zeroing(velp)!,3*maxsys%nax*maxsys%nsx)
    IF (cntl%quenche) CALL zeroing(cm)!,ngw*sh02%nsttot)
    IF (cntl%trevers) THEN
       ! invert electronic and ionic velocities (useful for path sampling)
       CALL dscal(2*ncpw%ngw*sh02%nsttot,-1._real_8,cm,1)
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
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,'(1X,64("="))')
       IF (paral%io_parent)&
            WRITE(6,'(1X,"==",T25,A,T64,"==")') 'FORCES INITIALIZATION'
       IF (paral%io_parent)&
            WRITE(6,'(1X,64("="))')
    ENDIF
    IF (geq0) CALL zclean(c0,sh02%nsttot,ncpw%ngw)

    ! McB ... surface hopping stuff ...
    CALL state_select("S0")
    CALL forcedr(c0(:,:,1),c2,sc0,rhoe,psi,&
         TAU0,FION0,EIGV,&
         sh02%nst_s0,1,.FALSE.,.TRUE.)
    e(1)=ener_com%etot

    ! Check orthogonality condition for wavefunction velocities (S0)
    CALL rortv(c0,cm,c2,sc0,gamy,sh02%nst_s0)

    CALL state_select("S1")
    CALL zeroing(rhoe)!,nnr1*clsd%nlsd)
    CALL zeroing(psi)!,nnr1*clsd%nlsd)
    CALL forcedr(c0(:,ns1:ns1+sh02%nst_s1-1,1),c2(:,ns1:ns1+sh02%nst_s1-1),sc0,rhoe,psi,&
         TAU0,FION1,EIGV,&
         sh02%nst_s1,1,.FALSE.,.TRUE.)

    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE (6,*) 'mdshop_bo: ',' EADDSH = ',sh02%eaddsh
    ENDIF

    e(2)=ener_com%etot
    e(2)=e(2)+sh02%eaddsh

    ! Check orthogonality condition for wavefunction velocities (S1)
    CALL rortv(c0(1,ns1,1),cm(1,ns1),c2(1,ns1),sc0,gamy,sh02%nst_s1)
    ! McB ............................................................

    ! Initialize thermostats
    IF (paral%parent) THEN
       itemp=irec(irec_nop1)+irec(irec_nop2)+irec(irec_nop3)&
            +irec(irec_nop4)
       IF (cntl%tnosep .AND. itemp.EQ.0) CALL nospinit(1)
       IF (cntl%tnosee .AND. irec(irec_noe) .EQ.0) CALL noseinit(1)
       CALL wrgeof(tau0,fion)
       filen=fo_info%fpath(fo_info%iapath:fo_info%iepath)//'ENERGIES'
       IF (paral%io_parent)&
            CALL fileopen(3,filen,fo_app+fo_verb,ferror)
       ! ... SURFACE HOPPING DATA ...
       filen=fo_info%fpath(fo_info%iapath:fo_info%iepath)//'ENERGIES_SH'
       IF (paral%io_parent)&
            CALL fileopen(32,filen,fo_app+fo_verb,ferror)
       CALL write_shmd(0,32,infi,tempp)
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

    ! McB ... surface hopping stuff ...
    IF (paral%parent) THEN
       ! ...     DEFINE INITIAL CONDITIONS FOR INTEGRATION OF STATE POPULATIONS
       IF (paral%io_parent)&
            INQUIRE(file='RESH',exist=lexist)
       IF (lexist) THEN
          IF (paral%io_parent)&
               OPEN(33,file='RESH')
          IF (paral%io_parent)&
               READ(33,*)sh03%isurf
          DO i=1,6
             IF (paral%io_parent)&
                  READ(33,*)sh03%pop(i)
          ENDDO
          IF (paral%io_parent)&
               CLOSE(33)
       ELSE
          IF (.NOT.sh03%tshopres) THEN
             ! ... standard setup ...
             sh03%isurf=2
             DO i=1,6
                sh03%pop(i)=0.0_real_8
             ENDDO
             sh03%pop(2)=1.0_real_8
          ENDIF
       ENDIF

       IF (paral%io_parent)&
            WRITE(6,'(/1X,''SURFACE HOPPING: INITIAL STATE ENERGIES'')')
       IF (paral%io_parent)&
            WRITE(6,'(5X,'' E(S0) = '',F20.10,''A.U.'')') e(1)
       IF (paral%io_parent)&
            WRITE(6,'(5X,'' E(S1) = '',F20.10,''A.U.'')') e(2)
       IF (paral%io_parent)&
            WRITE(6,'(5X,'' CURRENT STATE: '',I5/)') sh03%isurf
       IF (paral%io_parent)&
            WRITE(6,'(5X,''       POP(1) = '',F20.10)') sh03%pop(1)
       IF (paral%io_parent)&
            WRITE(6,'(5X,''       POP(2) = '',F20.10)') sh03%pop(2)
       IF (paral%io_parent)&
            WRITE(6,'(5X,''       POP(3) = '',F20.10)') sh03%pop(3)
       IF (paral%io_parent)&
            WRITE(6,'(5X,''       POP(4) = '',F20.10)') sh03%pop(4)
       IF (paral%io_parent)&
            WRITE(6,'(5X,''       POP(5) = '',F20.10)') sh03%pop(5)
       IF (paral%io_parent)&
            WRITE(6,'(5X,''       POP(6) = '',F20.10)') sh03%pop(6)
       IF (paral%io_parent)&
            WRITE(6,'(/1X,''SURFACE HOPPING: PREVIOUS STEP RESULTS'')')
       IF (paral%io_parent)&
            WRITE(6,'(5X,'' DET, DETOLD'',2F20.10)') sh03%det ,sh03%detold
       IF (paral%io_parent)&
            WRITE(6,'(5X,'' COUPL(,)'',   2F20.10)') sh03%coupl(1,1),sh03%coupl(1,2)
       IF (paral%io_parent)&
            WRITE(6,'(5X,'' COUPL(,)'',   2F20.10)') sh03%coupl(2,1),sh03%coupl(2,2)
       IF (paral%io_parent)&
            WRITE(6,'(5X,'' E1, E2  '',   2F20.10)') sh03%ec(1),sh03%ec(2)
    ENDIF
    CALL mp_bcast(sh03%pop,SIZE(sh03%pop),parai%source,parai%allgrp)
    CALL mp_bcast(sh03%isurf,parai%source,parai%allgrp)
    CALL decide(e,fion0,fion1,.FALSE.)
    ! McB ............................................................
    ! ==--------------------------------------------------------------==
    ! == END INITIALIZATION                                           ==
    ! ==--------------------------------------------------------------==
    IF ( teststore(0).AND.cntl%tsampl ) THEN
       ! McB... cf. elct.inc
       ntmp=crge%n
       crge%n=sh02%nsttot
       CALL zhwwf(2,irec,c0,cm,sh02%nsttot,eigv,taup,velp,taui,iteropt%nfi)
       crge%n=ntmp
    ENDIF
    IF (paral%parent) THEN
       time2 =m_walltime()
       tcpu = (time2 - time1)*0.001_real_8
       IF (paral%io_parent)&
            WRITE(6,'(A,T50,F8.2,A8)') ' TIME FOR INITIALIZATION:',&
            tcpu,' SECONDS'
       IF (paral%io_parent)&
            WRITE(6,'(//,1X,64("="))')
       IF (paral%io_parent)&
            WRITE(6,'(1X,"=",T20,A,T65,"=")')'MOLECULAR DYNAMICS SIMULATION'
       IF (paral%io_parent)&
            WRITE(6,'(1X,64("="))')
       IF ( tshopold ) THEN
          IF (paral%io_parent)&
               WRITE(6,'(/1X,''MDSHOP_CP:  USING OLD RK4OV! '')')
       ELSE
          IF (paral%io_parent)&
               WRITE(6,'(/1X,''MDSHOP_CP:  USING NEW RK4OV! '')')
       ENDIF
    ENDIF
    ! ==================================================================
    ! ==          THE BASIC LOOP FOR MOLECULAR DYNAMICS               ==
    ! ==                 USING VELOCITY VERLET                        ==
    ! ==================================================================

    ! NNNNN THIS SYNC SHOULD BE HERE !
    ALLOCATE(center(4,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    DO infi=1,cnti%nomore
       CALL mp_sync(parai%allgrp)
       time1=m_walltime()
       iteropt%nfi=iteropt%nfi+1

       ! McB ... surface hopping stuff ...
       ! ... store values at time (t-dt), (t-2dt) ...
       DO i=1,2
          DO j=1,2
             sh03%couplold(i,j)=sh03%coupl(i,j)
          ENDDO
       ENDDO
       detold2=sh03%detold
       sh03%detold=sh03%det
       ! McB ............................................................

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
       CALL anneal(velp,cm,sh02%nsttot,scr)
       CALL berendsen(velp,cm,nstate,scr,ekinc,0.0_real_8)
       ! UPDATE NOSE THERMOSTATS
       CALL noseup(velp,cm,sh02%nsttot,1)
       ! UPDATE VELOCITIES
       IF (paral%parent) CALL velupi(velp,fion,1)
       CALL velupa(c0,cm,c2,sh02%nsttot,1)
       ! SUBTRACT CENTER OF MASS VELOCITY
       IF (paral%parent.AND.comvl%subcom) CALL comvel(velp,vcmio,.TRUE.)
       ! SUBTRACT ROTATION AROUND CENTER OF MASS
       IF (paral%parent.AND.comvl%subrot) CALL rotvel(tau0,velp,lmio,tauio,.TRUE.)
       ! UPDATE POSITIONS
       ener_com%ecnstr = 0.0_real_8
       ener_com%erestr = 0.0_real_8
       IF (paral%parent) THEN
          CALL posupi(tau0,taup,velp)
          IF (cotc0%mcnstr.NE.0) THEN
             CALL cpmdshake(tau0,taup,velp)
          ENDIF
       ENDIF
       CALL mp_bcast(taup,SIZE(taup),parai%source,parai%allgrp)
       CALL phfac(taup)

       IF (corel%tinlc) CALL copot(rhoe,psi,ropt_mod%calste)

       ! McB ... surface hopping stuff ...
       ! ... select S0 ...
       CALL state_select("S0")
       CALL posupa(c0,cm,c2,gamx,sh02%nst_s0)
       ! ..Dipole moment
       IF (cntl%caldip) THEN
          CALL ddipo(taup,c0(:,:,1),cm,c2,sc0,sh02%nst_s0,center)
          CALL wannier_print(iteropt%nfi,c0(:,:,1),taup,sh02%nst_s0,psi(:,1),center)
       ENDIF
       ! CALCULATE THE FORCES
       ropt_mod%calste=cntl%tpres.AND.MOD(iteropt%nfi-1,cnti%npres).EQ.0
       CALL forcedr(c0(:,:,1),c2,sc0,rhoe,psi,&
            TAUP,FION0,EIGV,&
            sh02%nst_s0,1,.FALSE.,.TRUE.)
       e(1)=ener_com%etot
       IF (ropt_mod%calste) CALL totstr

       ! ... switch to S1 ...
       CALL state_select("S1")
       CALL posupa(c0(1,ns1,1),cm(1,ns1),c2(1,ns1),gamx,sh02%nst_s1)
       ! ..Dipole moment        
       IF (cntl%caldip) THEN
          CALL ddipo(taup,c0(:,ns1:,1),cm(:,ns1:),c2(:,ns1:),sc0,&
               sh02%nst_s1,center)
          CALL wannier_print(iteropt%nfi,c0(:,ns1:ns1+sh02%nst_s1-1,1),taup,sh02%nst_s1,psi(:,1),center)
       ENDIF
       ! CALCULATE THE FORCES
       ropt_mod%calste=cntl%tpres.AND.MOD(iteropt%nfi-1,cnti%npres).EQ.0
       CALL forcedr(c0(:,ns1:ns1+sh02%nst_s1-1,1),c2(:,ns1:ns1+sh02%nst_s1-1),sc0,rhoe,psi,taup,fion1,&
            eigv,sh02%nst_s1,1,.FALSE.,.TRUE.)
       e(2)=ener_com%etot
       e(2)=e(2)+sh02%eaddsh

       IF (ropt_mod%calste) CALL totstr

       ! -PARALLEL
       CALL decide(e,fion0,fion1,.FALSE.)
       ! -ENDPARALLEL        
       ! McB ............................................................

       ! ==================================================================
       ! Damped Dynamics
       CALL dampdyn(velp,fion,cm,c2,nstate,scr(1),scr(10))
       ! ==================================================================
       ! FINAL UPDATE FOR VELOCITIES
       ener_com%ecnstr = 0.0_real_8
       ener_com%erestr = 0.0_real_8
       IF (paral%parent) THEN
          CALL velupi(velp,fion,1)
          CALL rattle(taup,velp)
       ENDIF
       ! CALL VELUPA(C0,CM,C2,NSTTOT,1)  
       CALL velupa(c0,cm,c2,sh02%nsttot,1)
       CALL rortv(c0       ,cm       ,c2,&
            SC0,GAMY,sh02%nst_s0)
       CALL rortv(c0(1,ns1,1),cm(1,ns1),c2(1,ns1),&
            SC0,GAMY,NSTATE)
       IF (paral%parent) CALL geofile(taup,velp,'WRITE')
       ! COMPUTE THE IONIC TEMPERATURE TEMPP
       IF (paral%parent) THEN
          CALL ekinpp(ekinp,velp)
          tempp=ekinp*factem*2._real_8/glib
       ENDIF
       ! IONIC TEMPERATURE CONTROL
       IF (paral%parent) CALL rscvp(temp1,temp2,tempp,velp)
       ! SUBTRACT ROTATION AROUND CENTER OF MASS
       IF (paral%parent.AND.comvl%subrot) CALL rotvel(tau0,velp,lmio,tauio,.FALSE.)
       ! SUBTRACT CENTER OF MASS VELOCITY
       IF (paral%parent.AND.comvl%subcom) CALL comvel(velp,vcmio,.FALSE.)
       ! UPDATE NOSE THERMOSTATS
       IF (cntl%tnosee.OR.cntl%tc) CALL rekine(cm,sh02%nsttot,ekinc)
       CALL noseup(velp,cm,sh02%nsttot,1)
       CALL berendsen(velp,cm,nstate,scr,ekinc,0.0_real_8)
       ! ANNEALING
       CALL anneal(velp,cm,sh02%nsttot,scr)
       IF (paral%parent) THEN
          CALL ekinpp(ekinp,velp)
          tempp=ekinp*factem*2._real_8/glib
       ENDIF
       ! RESCALE ELECTRONIC VELOCITIES
       IF (cntl%tc) CALL rscve(ekin1,ekin2,ekinc,cntr%ekinw,cm,sh02%nsttot,ncpw%ngw)
       ! MEAN SQUARE DISPLACEMENT OF DIFFERENT IONIC SPECIES
       IF (paral%parent) CALL dispp(taup,taui,disa)
       ! KINETIC ENERGY OF THE ELECTRONS
       IF (sh03%isurf.EQ.1)CALL rekine(cm,crge%n,ekinc)
       IF (sh03%isurf.EQ.2)CALL rekine(cm(1,ns1),ns1,ekinc)
       ! ENERGY OF THE NOSE THERMOSTATS
       IF (paral%parent) CALL noseng(iteropt%nfi,velp,enose,enosp,dummy,1)
       IF (paral%parent) THEN
          econs=ekinp+ener_com%etot+enose+enosp+ener_com%ecnstr+ener_com%erestr
          eham=econs+ekinc
          time2=m_walltime()
          tcpu=(time2-time1)*0.001_real_8
          ! PRINTOUT the evolution of the accumulators every time step
          ! ND
          CALL wrprint_md(eigv,crge%f,ener_com%amu,nstate,taup,fion,&
               ekinc,tempp,ener_com%etot,econs,eham,disa,&
               tcpu,.FALSE.,iteropt%nfi,infi)
          ! /ND
          ! UPDATE ACCUMULATORS
          CALL pacca(ekinc,tempp,ener_com%etot,econs,eham,enose,enosp,ener_com%ecnstr,&
               ener_com%erestr,disa,tcpu,iteropt%nfi,1)
          ! Store ionic coordinates and velocities for statistics
          ropt_mod%movie=rout1%mout .AND. MOD(iteropt%nfi-1,cnti%imovie).EQ.0
          ropt_mod%rprint=rout1%rout .AND. MOD(iteropt%nfi-1,cnti%ntraj).EQ.0
          ropt_mod%txyz=rout1%xtout .AND. MOD(iteropt%nfi-1,cnti%ntraj).EQ.0
          ropt_mod%tdcd=rout1%dcout .AND. MOD(iteropt%nfi-1,cnti%ntraj).EQ.0
          CALL printp(taur,taup,velp)
       ENDIF

       ! McB ... surface hopping stuff ...
       ! ... compute DELTA(C0)/DELTA(t)=CM(t+dt/2) ...
       CALL s0_s1_overlap(c0,cm,sh03%det,c12,c21)
       sh03%coupl(1,2)=c12
       sh03%coupl(2,1)=c21
       DO i=1,2
          sh03%eold(i)=sh03%ec(i)
          sh03%ec(i)=e(i)
       ENDDO

       IF ( paral%parent ) THEN
          IF (paral%io_parent)&
               WRITE(6,*)'DET',sh03%det
          IF (paral%io_parent)&
               WRITE(6,*)'DELTA(DET) ',(sh03%det-detold2)/(2._real_8*delt)
          IF (paral%io_parent)&
               WRITE(6,*)'analytic:  ',sh03%couplold(2,1)+sh03%couplold(1,2),&
               ' delta:',(sh03%det-DETOLD2)/(2._real_8*DELT)-&
               (sh03%couplold(2,1)+sh03%couplold(1,2))
          IF (paral%io_parent)&
               WRITE(6,*)'D_21 + D_12',sh03%coupl(2,1)+sh03%coupl(1,2)
          IF (paral%io_parent)&
               WRITE(6,*)'D_21, D_12',sh03%coupl(2,1),sh03%coupl(1,2)
          IF (paral%io_parent)&
               WRITE(6,*)'DELTA(E)=',sh03%ec(2)-sh03%ec(1)
       ENDIF
       ! IF (INFI.GT.1) THEN
       IF ( (infi.GT.1).OR.sh03%tshopres ) THEN
          ddet=(sh03%det-sh03%detold)/REAL(ncoef,kind=real_8)
          dcoupl(1,2)=(sh03%coupl(1,2)-sh03%couplold(1,2))/REAL(ncoef,kind=real_8)
          dcoupl(2,1)=(sh03%coupl(2,1)-sh03%couplold(2,1))/REAL(ncoef,kind=real_8)
          dei(1)=(sh03%ec(1)-sh03%eold(1))/REAL(ncoef,kind=real_8)
          dei(2)=(sh03%ec(2)-sh03%eold(2))/REAL(ncoef,kind=real_8)
          ! \begin{interpolation loop for population integration} 
          DO i=1,ncoef
             ei(1)=sh03%eold(1)+dei(1)*REAL((i-1),kind=real_8)
             ei(2)=sh03%eold(2)+dei(2)*REAL((i-1),kind=real_8)
             sh03%coupl(1,2)=sh03%couplold(1,2)+dcoupl(1,2)*REAL((i-1),kind=real_8)
             sh03%coupl(2,1)=sh03%couplold(2,1)+dcoupl(2,1)*REAL((i-1),kind=real_8)
             sh03%det=sh03%detold+ddet*REAL((i-1),kind=real_8)
             IF ( paral%parent ) THEN
                ! WRITE(6,*)'INTERPOLATION',I,DET
                IF ( tshopold ) THEN
                   CALL rk4ov_old(dtcoef,sh03%pop,sh03%coupl,dcoupl,ei,dei,&
                        sh03%det,DDET)
                ELSE
                   CALL rk4ov_new(dtcoef,sh03%pop,sh03%coupl,dcoupl,ei,dei,&
                        sh03%det,DDET)
                ENDIF
             ENDIF
          ENDDO
          ei(1)=sh03%eold(1)+dei(1)*REAL(ncoef,kind=real_8)
          ei(2)=sh03%eold(2)+dei(2)*REAL(ncoef,kind=real_8)
          sh03%coupl(1,2)=sh03%couplold(1,2)+dcoupl(1,2)*REAL(ncoef,kind=real_8)
          sh03%coupl(2,1)=sh03%couplold(2,1)+dcoupl(2,1)*REAL(ncoef,kind=real_8)
          sh03%det=sh03%detold+ddet*REAL(ncoef,kind=real_8)
          ! \end{interpolation loop for population integration} 
       ELSE
          prob1%d11dot=0.0_real_8
          prob1%d22dot=0.0_real_8
          prob1%d1sq=1.0e-20_real_8
          prob1%d2sq=1.0_real_8
       ENDIF
       ! ... decide on which surface to propagate ...
       CALL mp_bcast_byte(prob1, size_in_bytes_of(prob1),parai%source,parai%allgrp)   ! /PROB1/
       CALL decide(e,fion0,fion1,.TRUE.)
       ! McB ............................................................
       IF (paral%parent) THEN
          CALL write_shmd(1,32,infi,tempp)
       ENDIF

       ! McB   ....................................
       IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)
       IF (infi.EQ.cnti%nomore) soft_com%exsoft=.TRUE.
       IF (teststore(iteropt%nfi).OR.soft_com%exsoft) THEN
          ! McB... cf. elct.inc
          ntmp=crge%n
          crge%n=sh02%nsttot
          CALL zhwwf(2,irec,c0,cm,sh02%nsttot,eigv,taup,velp,taui,iteropt%nfi)
          crge%n=ntmp
          ! McB ... surface hopping stuff ...
          ! ... write surface hopping restart file: 
          IF (paral%parent) THEN
             IF (paral%io_parent)&
                  OPEN(33,file='RESH')
             IF (paral%io_parent)&
                  WRITE(33,*)sh03%isurf
             DO i=1,6
                IF (paral%io_parent)&
                     WRITE(33,*)sh03%pop(i)
             ENDDO
             IF (paral%io_parent)&
                  CLOSE(33)
             IF (paral%io_parent)&
                  OPEN(33,file='COUPLSH')
             IF (paral%io_parent)&
                  WRITE(33,*) sh03%det,sh03%detold
             IF (paral%io_parent)&
                  WRITE(33,*) sh03%coupl(1,1), sh03%coupl(1,2)
             IF (paral%io_parent)&
                  WRITE(33,*) sh03%coupl(2,1), sh03%coupl(2,2)
             IF (paral%io_parent)&
                  WRITE(33,*) sh03%ec(1),sh03%ec(2)
             IF (paral%io_parent)&
                  CLOSE(33)
          ENDIF
          ! McB ............................................................
       ENDIF
       ! UPDATE IONIC POSITIONS
       CALL dcopy(3*maxsys%nax*maxsys%nsx,taup(1,1,1),1,tau0(1,1,1),1)
       ! temperature ramping
       CALL tempramp(temp1,temp2)
       ! STOP THE RUN IF THE USER HAS SET THE SIGNAL 30
       IF (soft_com%exsoft) GOTO 100

       ! ==================================================================
       ! ==     END OF MAIN LOOP                                         ==
       ! ==================================================================
    ENDDO
100 CONTINUE
    ! McB ... surface hopping stuff ...
    IF (paral%parent) CALL write_shmd(-1,32,infi,tempp)
    IF (wannl%twann) THEN
       CALL state_select("S0")
       CALL ddipo(taup,c0(:,:,1),cm,c2,sc0,sh02%nst_s0,center)
       CALL forcedr(c0(:,:,1),c2,sc0,rhoe,psi,taup,fion,eigv,&
            sh02%nst_s0,1,.FALSE.,.TRUE.)
       CALL wc_dos(c0,c2,sh02%nst_s0,center)
       CALL state_select("S1")
       CALL ddipo(taup,c0(:,ns1:,1),cm(:,ns1:),c2(:,ns1:),sc0,sh02%nst_s1,center)
       CALL forcedr(c0(:,ns1:ns1+sh02%nst_s1-1,1),c2(:,ns1:ns1+sh02%nst_s1-1),sc0,rhoe,psi,taup,fion,&
            eigv,&
            sh02%nst_s1,1,.FALSE.,.TRUE.)
       CALL wc_dos(c0,c2,sh02%nst_s1,center)
       ! McB... cf. elct.inc
       ntmp=crge%n
       crge%n=sh02%nsttot
       CALL zhwwf(2,irec,c0,cm,sh02%nsttot,eigv,taup,velp,taui,iteropt%nfi)
       crge%n=ntmp
    ENDIF
    DEALLOCATE(center,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! McB ............................................................
    ! PRINT ACCUMULATOR
    IF (paral%parent) CALL pacca(ekinc,tempp,ener_com%etot,econs,eham,enose,enosp,&
         ener_com%ecnstr,ener_com%erestr,disa,tcpu,iteropt%nfi,2)
    ! McB ... surface hopping stuff ...
    CALL proja(c0,c2,sc0,sh02%nst_s0,cnti%iproj)
    CALL proja(c0(1,ns1,1),c2(1,ns1),sc0,sh02%nst_s1,cnti%iproj)
    ! McB ............................................................
    CALL csize(c2,sh02%nsttot,gemax,cnorm)
    IF (paral%parent) CALL gsize(fion,gnmax,gnorm)
    IF (paral%parent) CALL finalp(tau0,fion,velp,eigv)
    ! 
    IF (cntl%tsampl) THEN
       CALL sample_go
       GOTO 99999
    ENDIF
10000 CONTINUE
    ! 
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
    IF ((paral%parent).AND.paral%io_parent)&
         CALL fileclose(32)
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
  END SUBROUTINE mdshop_cp
  ! ==================================================================
  SUBROUTINE give_scr_mdshop_cp(lmdshop,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lmdshop
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: lcopot, lddipo, ldeort, &
                                                lforcedr, linitrun, lortho, &
                                                lposupa, lquenbo, lrhopri, &
                                                lrortv, nstate

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
    lmdshop=MAX(lcopot,lortho,lquenbo,ldeort,lforcedr,&
         lrortv,lposupa,lrhopri,lddipo,linitrun)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_mdshop_cp

END MODULE mdshop_cp_utils
