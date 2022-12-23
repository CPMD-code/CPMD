MODULE mm_mdshop_cp_utils
  USE anneal_utils,                    ONLY: anneal,&
                                             berendsen,&
                                             dampdyn,&
                                             tempramp
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
  USE cotr,                            ONLY: cotc0
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
                                             fo_info,&
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
  USE machine,                         ONLY: m_walltime
  USE meta_colvar_inp_utils,           ONLY: colvar_structure
  USE meta_colvar_utils,               ONLY: meta_colvar
  USE meta_exl_mult_utils,             ONLY: meta_ext_mul
  USE meta_exlagr_methods,             ONLY: meta_extlagr
  USE meta_exlagr_utils,               ONLY: ekincv_global
  USE mm_dim_utils,                    ONLY: mm_dim
  USE mm_dimmod,                       ONLY: clsaabox,&
                                             mm_go_mm,&
                                             mm_go_qm,&
                                             mm_revert,&
                                             mmdim
  USE mm_input,                        ONLY: clc,&
                                             g96_vel,&
                                             lqmmm,&
                                             rtr_l,&
                                             wp_i
  USE mm_mdmain_utils,                 ONLY: mm_localt
  USE mm_parallel,                     ONLY: gparal
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
  USE rk4ov_utils,                     ONLY: rk4ov_new,&
                                             rk4ov_old
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
       irec_vel, irec_wf, restart1, rout1
  USE system,                          ONLY: &
       cnti, cntl, cntr, fpar, iatpt, maxsys, nacc, ncpw, nkpt, restf
  USE testex_utils,                    ONLY: testex
  USE teststore_utils,                 ONLY: teststore
  USE totstr_utils,                    ONLY: totstr
  USE tpar,                            ONLY: dt_ions
  USE tst2min_utils,                   ONLY: tst2min
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

  PUBLIC :: mm_mdshop_cp
  PUBLIC :: give_scr_mm_mdshop_cp

CONTAINS

  ! ==================================================================
  SUBROUTINE mm_mdshop_cp(c0,cm,c2,sc0,gamx,gamy)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8) :: c0(ncpw%ngw,crge%n,1), cm(ncpw%ngw,crge%n), &
      c2(ncpw%ngw,crge%n), sc0(ncpw%ngw,crge%n)
    REAL(real_8)                             :: gamx(*), gamy(*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'mm_mdshop_cp'

    CHARACTER(len=10)                        :: prch
    CHARACTER(len=100)                       :: filen, filensh
    CHARACTER(len=30)                        :: tag
    COMPLEX(real_8), ALLOCATABLE             :: psi(:,:)
    INTEGER :: i, ia, ierr, il_psi_1d, il_psi_2d, il_rhoe_1d, il_rhoe_2d, in, &
      irec(100), is, itemp, j, k, lscr, ncoef, ns1, nstate, ntmp, out
    LOGICAL                                  :: ferror, lexist, lmetares, &
                                                lquench, oldstatus, resetcv, &
                                                statusdummy
    REAL(real_8) :: c12, c21, DCOUPL(2,2), ddet, dei(2), delt, detold2, disa, &
      dt, dtcoef, dummy, e(2), econs, eham, ei(sh02%nsurf), ek_cv, ekin1, &
      ekin2, ekinc, ekincp, ekinh1, ekinh2, ekinp, enose, enosp, ff, lmio(3), &
      tcpu, temp1, temp2, tempp, time1, time2, vcmio(4)
    REAL(real_8), ALLOCATABLE :: center(:,:), eigm(:), eigv(:), fion0(:,:,:), &
      fion1(:,:,:), rhoe(:,:), scr(:), taui(:,:,:), tauio(:,:), taur(:,:,:), &
      tscr(:,:,:)
    REAL(real_8), DIMENSION(:, :, :), &
      POINTER                                :: scr_correct_shape

    COMMON /ekin/ekinp
    ! ==================================================================
    NULLIFY( scr_correct_shape )
#if defined (__GROMOS)
    time1 =m_walltime()
    CALL mm_dim(mm_go_mm,oldstatus)

    IF (paral%qmnode) THEN
       IF (textfld)THEN
          ALLOCATE(extf(fpar%kr1*fpar%kr2s*fpar%kr3s),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          CALL zeroing(extf)!,kr1*kr2s*kr3s)
       ENDIF
       nstate=crge%n
       nacc = 22
       ! McB ... surface hopping stuff ...
       delt=cntr%delt_ions
       sh03%detold=0.0_real_8
       detold2=0.0_real_8
       sh03%isurf=2
       ncoef=INT(delt/0.04_real_8)+1
       dtcoef=delt/REAL(ncoef,kind=real_8)
       IF (paral%parent) THEN
          IF (paral%io_parent)&
               WRITE(6,*)'NCOEF=',ncoef,'DTCOEF=',dtcoef
          IF (paral%io_parent)&
               WRITE(6,*)'DT=',delt,'DTCOEF*NCOEF=',dtcoef*ncoef
       ENDIF
       tlsd0=cntl%tlsd
       crge%n=sh02%nst_s1
       clsd%nlsd=4
       clsd%nlsx=3
       nstate=sh02%nst_s1
       sh02%nsttot=sh02%nst_s0 + sh02%nst_s1
       ns1=sh02%nst_s0+1
       iteropt%nfi=0
       ! .................................
       ropt_mod%modens=.FALSE.
       ropt_mod%engpri=.FALSE.
       ropt_mod%calste=cntl%tpres

       IF (cntl%tqmmm.AND.pslo_com%tivan) THEN
          CALL stopgm('MM_MDSHOP_CP','QMMM WITH VANDERBILT UNSUPPORTED',& 
               __LINE__,__FILE__)
       ENDIF

       ! McB ... surface hopping stuff ...
       ! CALL MEMORY(IP_EIGV,NSTATE,'EIGV')     ! BO only ?!
       ! CALL MEMORY(IP_EIGV1,NSTATE,'EIGV1')   ! BO only ?!
       ! .................................

       IF (pslo_com%tivan) cntl%nonort=.TRUE.
       IF (cntl%nonort) THEN
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
                CALL stopgm('MDMAIN',' ',& 
                     __LINE__,__FILE__)
             ENDIF
          END DO
       ENDIF
       ! ==--------------------------------------------------------------==
       ! SCR ALLOCATION AND PARTITION (SCRATCH ARRAY).
       CALL rhoe_psi_size(il_psi_1d=il_psi_1d,il_psi_2d=il_psi_2d,&
            il_rhoe_1d=il_rhoe_1d,il_rhoe_2d=il_rhoe_2d)
       ALLOCATE(rhoe(il_rhoe_1d,il_rhoe_2d),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ! McB  CALL MEMORY(IP_RHOE1,IL_RHOE,'RHOE1')
       ALLOCATE(psi(il_psi_1d,il_psi_2d),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ! McB  CALL MEMORY(IP_PSI1,IL_PSI,'PSI1')
       CALL give_scr_mm_mdshop_cp(lscr,tag)
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
    CALL mm_dim(mm_go_mm,statusdummy)

    ! we may be able to save some memory by allocating more carefully
    ! on every node a different size. Could be a bad idea.

    ALLOCATE(fion(3,maxsys%nax,maxsys%nsx),tscr(3,maxsys%nax,maxsys%nsx),STAT=ierr)
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
    IF (comvl%tsubrot)  THEN
       ALLOCATE(tauio(3,ions1%nat),STAT=ierr)
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

    CALL zeroing(fion)!,3*maxsys%nax*maxsys%nsx)
    CALL zeroing(taui)!,3*maxsys%nax*maxsys%nsx)
    CALL zeroing(taup)!,3*maxsys%nax*maxsys%nsx)

    ! McB ... surface hopping stuff ...
    ALLOCATE(fion0(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(fion1(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(sh03%coupl)!,2*2)
    CALL zeroing(sh03%couplold)!,2*2)
    CALL zeroing(fion0)!,3*maxsys%nax*maxsys%nsx)
    CALL zeroing(fion1)!,3*maxsys%nax*maxsys%nsx)
    DO i=1,2
       e(i)=0.0_real_8
       sh03%ec(i)=0.0_real_8
       sh03%eold(i)=0.0_real_8
    END DO
    sh03%det=0.0_real_8
    sh03%detold=0.0_real_8
    detold2=0.0_real_8
    sh03%isurf=0
    ! .................................

    ! AK 2005/05/06: FIXME: we need to loop over processes here,
    ! not have a garbled output. at the moment this is just
    ! disabled, because it is _very_ annoying.
    CALL mp_sync(parai%qmmmgrp)
#if 0
    IF (paral%io_parent) WRITE(prch,'(a6,i4)')'NODE: ',parai%me
    CALL mp_sync(parai%qmmmgrp)
#endif
    IF (paral%qmnode ) THEN

       CALL initrun(irec,c0,cm,sc0,rhoe,psi,eigv)

       ! time step functions
       CALL dynit(ekincp,ekin1,ekin2,temp1,temp2,ekinh1,ekinh2)
       ekinc=0.0_real_8  ! McB: cf. META_EXT..()
       ekinp=0.0_real_8  ! McB: cf. META_EXT..()
       ! parameters for the nose-hoover thermostats
       IF ((cntl%tnosee.OR.cntl%tnosep).AND.paral%parent) CALL nosepa(1,1)

       ! Dont symmetrize density 
       cntl%tsymrho=.FALSE.
       CALL mp_bcast(taup,SIZE(taup),parai%source,parai%allgrp)

       irec(irec_wf)=1

       CALL dcopy(3*maxsys%nax*maxsys%nsx,taup,1,taui,1)
    ENDIF
    IF (cntl%quenchb) THEN

       CALL mm_dim(mm_go_mm,statusdummy)
       CALL zeroing(cm(:,1:sh02%nsttot))!,ngw*sh02%nsttot)

       CALL state_select("S0")
       IF (.NOT.lqmmm%qmmm_reflex) THEN
          CALL mm_translate_qmmm(tau0,c0,cm,sh02%nst_s0)
       ENDIF
       CALL mm_dim(mm_go_qm,statusdummy)
       CALL quenbo(c0(:,:,1),c2,sc0,taur,rhoe,psi)
       CALL mm_dim(mm_go_mm,statusdummy)

       CALL state_select("S1")
       IF (.NOT.lqmmm%qmmm_reflex) THEN
          CALL mm_translate_qmmm(tau0,c0(1,ns1,1),cm,sh02%nst_s1)
       ENDIF
       CALL mm_dim(mm_go_qm,statusdummy)
       CALL quenbo(c0(:,ns1:,1),c2,sc0,taur,rhoe,psi)
       CALL mm_dim(mm_go_mm,statusdummy)
       IF ( paral%qmnode ) THEN
          ! McB    ... cf. elct.inc
          ntmp=crge%n
          crge%n=sh02%nsttot
          CALL zhwwf(2,irec,c0,cm,sh02%nsttot,eigv,tau0,velp,taui,iteropt%nfi)
          crge%n=ntmp
       ENDIF
    ENDIF
    IF (paral%qmnode ) THEN
       CALL dcopy(3*maxsys%nax*maxsys%nsx,taup,1,taui,1)
       IF (pslo_com%tivan) THEN
          CALL mm_dim(mm_go_qm,statusdummy)
          IF (cntl%tlsd) THEN
             ! McB  ... what about this in mdmain.F ?!?
             ! CALL DEORT(NGW,NSUP,EIGM,EIGV,C0(1,1),SC0(1,1))
             ! CALL DEORT(NGW,NSDOWN,EIGM,EIGV,C0(1,NSUP+1),SC0)
             ! *               SCR,LSCR)
             ! McB  ................................?!?
          ELSE
             CALL deort(ncpw%ngw,sh02%nst_s0,eigm,eigv,c0(1,1,1)  ,sc0)
             CALL deort(ncpw%ngw,sh02%nst_s1,eigm,eigv,c0(1,ns1,1),sc0)
          ENDIF
       ENDIF
       ! INITIALIZE VELOCITIES
       CALL mm_dim(mm_go_mm,statusdummy)
       IF (paral%parent) CALL detdof(tau0,tscr)

       ! INITIALIZE METADYNAMICS VARIABLES used also for 
       ! SHOOTING from SADDLE POINT with RANDOM VELOCITIES
       IF (paral%parent .AND. (lmeta%lcolvardyn .OR. lmeta%lsadpnt)) THEN
          CALL colvar_structure(tau0,taur)
       ENDIF

       IF ( (irec(irec_vel).EQ.0).AND.&
            (.NOT.restart1%rgeo).AND.(.NOT.rtr_l%restart_traj) ) THEN
          ener_com%ecnstr = 0.0_real_8
          ener_com%erestr = 0.0_real_8
          CALL rinvel(velp,cm,crge%n)
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
       IF (cntl%quenche) CALL zeroing(cm(:,1:sh02%nsttot))!,ngw*sh02%nsttot)
       IF (clc%classical) CALL zeroing(cm(:,1:sh02%nsttot))!,ngw*sh02%nsttot)
       ! RESET ACCUMULATORS
       IF (paral%parent.AND.irec(irec_ac).EQ.0) CALL resetac(tau0,taui,iteropt%nfi)
       ! INITIALIZE FORCES
       IF ( paral%parent ) THEN
          IF (paral%io_parent)&
               WRITE(6,'(1X,64("="))')
          IF (paral%io_parent)&
               WRITE(6,'(1X,"==",T25,A,T64,"==")') 'FORCES INITIALIZATION'
          IF (paral%io_parent)&
               WRITE(6,'(1X,64("="))')
       ENDIF
       CALL mm_dim(mm_go_qm,statusdummy)
       IF (tkpts%tkpnt) THEN
          ! McB  ... what about this in mdmain.F ?!?
          ! IF(GEQ0) CALL ZCLEAN_K(C0,NSTATE,NGW)
          ! McB  ............................... ?!?
       ELSE
          IF (geq0) CALL zclean(c0,sh02%nsttot,ncpw%ngw)
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

    ! McB ... surface hopping stuff ...
    ! ...   initial forces on s0: FION0
    CALL state_select("S0")
    CALL mm_qmmm_forcedr(c0,c2,sc0,rhoe,psi,&
         TAU0,FION0,EIGV,&
         sh02%nst_s0,0,.FALSE.,.TRUE.,.FALSE.)
    e(1)=ener_com%etot

    IF (paral%qmnode) THEN
       ! Check orthogonality condition for wavefunction velocities (S0)
       CALL mm_dim(mm_go_qm,statusdummy)
       CALL freqs(crge%n,.TRUE.)
       CALL rortv(c0,cm,c2,sc0,gamy,sh02%nst_s0)
       CALL mm_dim(mm_go_mm,statusdummy)
    ENDIF

    ! ...   initial forces on s1: FION1
    CALL state_select("S1")
    CALL zeroing(rhoe)!,il_rhoe)
    CALL zeroing(psi)!,SIZE(psi))     ! PSI is *16, but IL_PSI counts in *8
    CALL mm_qmmm_forcedr(c0(:,ns1:ns1+sh02%nst_s1-1,1),c2(:,ns1:ns1+sh02%nst_s1-1),sc0,rhoe,psi,&
         TAU0,FION1,EIGV,&
         sh02%nst_s1,0,.FALSE.,.TRUE.,.FALSE.)
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE (6,*) 'mdshop_cp: ',' EADDSH = ',sh02%eaddsh
    ENDIF

    e(2)=ener_com%etot
    e(2)=e(2)+sh02%eaddsh

    IF (paral%qmnode) THEN
       ! Check orthogonality condition for wavefunction velocities (S0)
       CALL mm_dim(mm_go_qm,statusdummy)
       CALL freqs(crge%n,.TRUE.)
       CALL rortv(c0(1,ns1,1),cm(1,ns1),c2(1,ns1),sc0,gamy,&
            sh02%nst_s1)
       CALL mm_dim(mm_go_mm,statusdummy)
    ENDIF
    ! McB ............................................................

    IF (paral%qmnode) THEN
       ! Initialize thermostats
       IF (paral%parent) THEN
          itemp=irec(irec_nop1)+irec(irec_nop2)+irec(irec_nop3)&
               +irec(irec_nop4)
          IF (cntl%tnosep .AND. itemp.EQ.0) CALL nospinit(1)
          IF (cntl%tnosee .AND. irec(irec_noe) .EQ.0) CALL noseinit(1)
          IF ( mmdim%natm.LE.500) CALL wrgeof(tau0,fion1)
          filen=fo_info%fpath(fo_info%iapath:fo_info%iepath)//'ENERGIES'
          IF (paral%io_parent)&
               CALL fileopen(3,filen,fo_app+fo_verb,ferror)
          ! ... SURFACE HOPPING DATA ...
          filensh=fo_info%fpath(fo_info%iapath:fo_info%iepath)//'ENERGIES_SH'
          IF (paral%io_parent)&
               CALL fileopen(32,filensh,fo_app+fo_verb,ferror)
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
          ! ...    DEFINE INITIAL CONDITIONS FOR INTEGRATION OF STATE POPULATIONS
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
             END DO
             IF (paral%io_parent)&
                  CLOSE(33)
          ELSE
             IF (.NOT.sh03%tshopres) THEN
                ! ... standard setup ...
                sh03%isurf=2
                DO i=1,6
                   sh03%pop(i)=0.0_real_8
                END DO
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
       ! .................................
       ! ==--------------------------------------------------------------==
       ! == END INITIALIZATION                                           ==
       ! ==--------------------------------------------------------------==
       IF (teststore(0).AND.cntl%tsampl) THEN
          ! McB    ... cf. elct.inc
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
                  WRITE(6,'(/1X,''MM_MDSHOP_CP:  USING OLD RK4OV! '')')
          ELSE
             IF (paral%io_parent)&
                  WRITE(6,'(/1X,''MM_MDSHOP_CP:  USING NEW RK4OV! '')')
          ENDIF
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
          ! McB ... surface hopping stuff ...
          ! ... store values at time (t-dt), (t-2dt) ...
          DO i=1,2
             DO j=1,2
                sh03%couplold(i,j)=sh03%coupl(i,j)
             END DO
          END DO
          detold2=sh03%detold
          sh03%detold=sh03%det
          ! McB ............................................................
          ! ANNEALING
          ! McB    ... currently  no broken symmetry with surface hopping ...
          ! BSCLCS=1
          CALL anneal(velp,cm,sh02%nsttot,scr)
          CALL berendsen(velp,cm,nstate,scr,ekinc,0.0_real_8)
          ! McB    ... currently  no broken symmetry with surface hopping ...
          ! IF(cntl%bsymm)THEN
          ! BSCLCS=2
          ! CALL ANNEAL(VELP,CM(1,1,2),NSTATE,SCR)
          ! endif
          ! McB    ..........................................................
          ! UPDATE NOSE THERMOSTATS
          CALL noseup(velp,cm,sh02%nsttot,1)
          ! SUBTRACT CENTER OF MASS VELOCITY
          IF (paral%parent.AND.comvl%subcom) CALL comvel(velp,vcmio,.TRUE.)
          ! SUBTRACT ROTATION AROUND CENTER OF MASS
          IF (paral%parent.AND.comvl%subrot) CALL rotvel(tau0,velp,lmio,tauio,.TRUE.)
          ! UPDATE VELOCITIES
          IF (paral%parent) CALL velupi(velp,fion,1)
          CALL velupa(c0,cm,c2,sh02%nsttot,1)
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
          CALL mm_translate_qmmm(taup,c0,cm,sh02%nsttot)
          ! mb     we must not translate the wavefunction
          ! mb     if we do not translate the atoms !
       ENDIF

       IF (paral%qmnode) THEN
          CALL mm_dim(mm_go_qm,statusdummy)
          CALL phfac(taup)
          IF (corel%tinlc) CALL copot(rhoe,psi,ropt_mod%calste)

          ! McB ... surface hopping stuff ...
          ! ... select S0 ...
          CALL state_select("S0")
          IF (.NOT.clc%classical)THEN
             CALL posupa(c0,cm,c2,gamx,sh02%nst_s0)
          ENDIF
          ! ..Dipole moment
          IF (cntl%caldip) THEN
             CALL ddipo(taup,c0(:,:,1),cm,c2,sc0,sh02%nst_s0,center)
             CALL wannier_print(iteropt%nfi,c0(:,:,1),taup,sh02%nst_s0,psi(:,1),center)
          ENDIF
          ropt_mod%calste=cntl%tpres.AND.MOD(iteropt%nfi-1,cnti%npres).EQ.0
          CALL mm_dim(mm_go_mm,statusdummy)
       ENDIF! qmnode

       CALL mm_qmmm_forcedr(c0,c2,sc0,rhoe,psi,&
            TAUP,FION0,EIGV,&
            sh02%nst_s0,0,.FALSE.,.TRUE.,.TRUE.)
       e(1)=ener_com%etot

       ! ... select S1 ...
       IF ( paral%qmnode ) THEN
          CALL mm_dim(mm_go_qm,statusdummy)
          CALL state_select("S1")
          IF ( .NOT.clc%classical ) THEN
             CALL posupa(c0(1,ns1,1),cm(1,ns1),c2(1,ns1),gamx,sh02%nst_s1)
          ENDIF
          ! ..Dipole moment
          IF (cntl%caldip.AND.(.NOT.cntl%bsymm)) THEN
             CALL ddipo(taup,c0(:,ns1:,1),cm(:,ns1:),c2(:,ns1:),sc0,&
                  sh02%nst_s1,center)
             CALL wannier_print(iteropt%nfi,c0(:,ns1:ns1+sh02%nst_s1-1,1),taup,sh02%nst_s1,psi(:,1),center)
          ENDIF
          ropt_mod%calste=cntl%tpres.AND.MOD(iteropt%nfi-1,cnti%npres).EQ.0
          CALL mm_dim(mm_go_mm,statusdummy)
       ENDIF! qmnode

       CALL mm_qmmm_forcedr(c0(1,ns1,1),c2(1,ns1),sc0,rhoe,psi,&
            TAUP,FION1,EIGV,&
            sh02%nst_s1,0,.FALSE.,.TRUE.,.TRUE.)
       e(2)=ener_com%etot
       e(2)=e(2)+sh02%eaddsh
       ! .................................
       ! ==================================================================
       ! McB ... skip this for a while ...
       go to 500
       ! McB ... 
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
             IF ((rmeta%tolkin.NE.-1.0_real_8 .AND. ekinc.GT.rmeta%tolkin)&
                  .OR. lquench) THEN
                ! Check for Quench on the BO surface
                IF (paral%parent) THEN
                   IF (paral%io_parent)&
                        WRITE(6,'(''TOLKIN ='',f16.8,''  EKINC ='',f16.8)')&
                        rmeta%tolkin,ekinc
                ENDIF
                CALL mm_dim(mm_go_qm,statusdummy)
                ! McB         ... currently  no broken symmetry with surface hopping ...
                ! BSCLCS=1
                ! IF(cntl%bsymm)CALL SETBSSTATE
                ! McB         ..........................................................
                CALL quenbo(c0(:,:,1),c2,sc0,taur,rhoe,psi)
                ! McB         ... currently  no broken symmetry with surface hopping ...
                ! IF(cntl%bsymm)THEN
                ! BSCLCS=2
                ! CALL SETBSSTATE
                ! CALL QUENBO(C0(1,1,2),C2(1,1,2),SC0(1,1,2),TAUR,RHOE,
                ! &                    PSI,SCR,LSCR)
                ! endif
                ! McB         ..........................................................
                CALL mm_dim(mm_go_mm,statusdummy)
                ! McB         ... currently  no broken symmetry with surface hopping ...
                ! CALL AZZERO(CM,2*NSTATE*NGWK*BSFAC)
                ! McB         ... currently  no broken symmetry with surface hopping ...
                resetcv = .TRUE.
             ENDIF
             IF (paral%parent)THEN
#if defined(__VECTOR) &&  ( !defined(__PRIMERGY) )
                !$omp parallel DO private(I,IS,IA)
                DO i=1,ions1%nat
                   ia=iatpt(1,i)
                   is=iatpt(2,i)
                   fion(1,ia,is) = fion(1,ia,is) + fhills(1,ia,is)
                   fion(2,ia,is) = fion(2,ia,is) + fhills(2,ia,is)
                   fion(3,ia,is) = fion(3,ia,is) + fhills(3,ia,is)
                END DO
#else
                ! ! !$OMP parallel do private(IS,IA) - fails on many compilers
                DO is=1,ions1%nsp
                   DO ia=1,ions0%na(is)
                      fion(1,ia,is) = fion(1,ia,is) + fhills(1,ia,is)
                      fion(2,ia,is) = fion(2,ia,is) + fhills(2,ia,is)
                      fion(3,ia,is) = fion(3,ia,is) + fhills(3,ia,is)
                   END DO
                END DO
#endif
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
                IF ((paral%parent).AND.paral%io_parent) CALL fileclose(3)
                GOTO   99999
             ENDIF
          ENDIF
       ENDIF

       ! ==================================================================
       ! McB
500    CONTINUE
       ! McB
       IF (paral%qmnode) THEN

          ! McB ... surface hopping stuff ...
          ! -PARALLEL
          CALL decide(e,fion0,fion1,.FALSE.)
          ! -ENDPARALLEL        
          ! McB ............................................................

          ! ==================================================================
          ! Damped Dynamics
          CALL dampdyn(velp,fion,cm,c2,nstate,scr(1),scr(10))
          ! ==================================================================
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
          CALL velupa(c0,cm,c2,sh02%nsttot,1)
          CALL mm_dim(mm_go_qm,statusdummy)
          CALL rortv(c0       ,cm       ,c2       ,&
               SC0,GAMY,sh02%nst_s0)
          CALL rortv(c0(1,ns1,1),cm(1,ns1),c2(1,ns1),&
               SC0,GAMY,sh02%nst_s1)
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
             ! calculate local kinetic temperature in QM and MM subsystem
             ! and write to "QM_TEMP"
             IF (wp_i%NFI_lt.GT.0) THEN
                CALL mm_localt(velp,rmass%pma,factem,tempp,glib,iteropt%nfi,wp_i%NFI_lt)
             ENDIF
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
             CALL rekine(cm,sh02%nsttot,ekinc)
          ENDIF
          CALL noseup(velp,cm,sh02%nsttot,1)
          CALL berendsen(velp,cm,nstate,scr,ekinc,0.0_real_8)
          ! ANNEALING
          CALL anneal(velp,cm,sh02%nsttot,scr)
          IF (paral%parent) THEN
             CALL ekinpp(ekinp,velp)
             IF (lmeta%lextlagrange.AND. ltcglobal) THEN
                CALL ekincv_global(ek_cv)
                tempp=(ek_cv+ekinp)*factem*2._real_8/(glib+REAL(ncolvar,kind=real_8))
             ELSE
                tempp=ekinp*factem*2._real_8/glib
             ENDIF
          ENDIF
          ! RESCALE ELECTRONIC VELOCITIES
          IF (cntl%tc) CALL rscve(ekin1,ekin2,ekinc,cntr%ekinw,cm,sh02%nsttot,ncpw%ngw)
          ! MEAN SQUARE DISPLACEMENT OF DIFFERENT IONIC SPECIES
          IF (paral%parent) CALL dispp(taup,taui,disa)
          ! KINETIC ENERGY OF THE ELECTRONS
          ! McB ... surface hopping stuff ...
          ! CALL REKINE(CM,NSTATE,EKINC)
          IF ( sh03%isurf.EQ.1 ) CALL rekine(cm       ,sh02%nst_s0,ekinc)
          IF ( sh03%isurf.EQ.2 ) CALL rekine(cm(1,ns1),sh02%nst_s1,ekinc)
          ! .................................
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
             ! PRINTOUT the evolution of the accumulators every time step
             CALL wrprint_md(eigv,crge%f,ener_com%amu,sh02%nsttot,taup,fion,&
                  EKINC,TEMPP,ener_com%etot,ECONS,EHAM,DISA,&
                  TCPU,.FALSE.,iteropt%nfi,infi)
             ! UPDATE ACCUMULATORS
             CALL pacca(ekinc,tempp,ener_com%etot,econs,eham,enose,enosp,ener_com%ecnstr,&
                  ener_com%erestr,disa,tcpu,iteropt%nfi,1)
             ! Store ionic coordinates and velocities for statistics
             ropt_mod%movie=rout1%mout .AND. MOD(iteropt%nfi-1,cnti%imovie).EQ.0
             ropt_mod%rprint=rout1%rout .AND. MOD(iteropt%nfi-1,cnti%ntraj).EQ.0
             ropt_mod%tdcd=rout1%dcout .AND. MOD(iteropt%nfi-1,cnti%ntraj).EQ.0
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

          ! McB   ... surface hopping stuff ...
          ! ... compute DELTA(C0)/DELTA(t)=CM(t+dt/2) ...
          CALL s0_s1_overlap(c0,cm,sh03%det,c12,c21)
          sh03%coupl(1,2)=c12
          sh03%coupl(2,1)=c21
          DO i=1,2
             sh03%eold(i)=sh03%ec(i)
             sh03%ec(i)=e(i)
          END DO

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
             END DO
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
          CALL mp_bcast_byte(prob1, size_in_bytes_of(prob1),parai%source,parai%allgrp)
          CALL decide(e,fion0,fion1,.TRUE.)
          ! McB ............................................................
          IF (paral%parent) THEN
             CALL write_shmd(1,32,infi,tempp)
          ENDIF

          IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)
          IF (infi.EQ.cnti%nomore) soft_com%exsoft=.TRUE.
          ! temperature ramping
          CALL tempramp(temp1,temp2)
       ENDIF! qmnode

       ! STOP THE RUN IF THE USER HAS SET THE SIGNAL 30
       in=0
       out=0
       IF (soft_com%exsoft) in=1
       CALL mp_sum(in,out,parai%qmmmgrp)
       IF (out.NE.0) soft_com%exsoft=.TRUE.

       IF (paral%qmnode) THEN
          IF (teststore(iteropt%nfi).OR.soft_com%exsoft.OR.lmetares) THEN
             ! McB      ... cf. elct.inc
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
                END DO
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
                ! Time dipendent potential applied directly on the Collective Variables
                CALL meta_colvar(taup,velp,fion,taur,&
                     lquench,lmetares,ekinc,ekinp)
             ENDIF
          ENDIF

       ENDIF! qmnode

       IF (paral%qmnode) THEN
          ! UPDATE IONIC POSITIONS
          CALL dcopy(3*maxsys%nax*maxsys%nsx,taup(1,1,1),1,tau0(1,1,1),1)
       ENDIF! qmnode

       IF (soft_com%exsoft) GOTO 100

       ! ==================================================================
       ! ==     END OF MAIN LOOP                                         ==
       ! ==================================================================
    END DO
100 CONTINUE
    CALL mm_dim(mm_go_mm,statusdummy)
    IF ( gparal%mmparent ) THEN
       CALL mm_write_gromos_coord('CRD_FIN.g96',taup,velp,maxsys%nax,maxsys%nsx)
    ENDIF
    IF (paral%qmnode) THEN
       ! McB ... surface hopping stuff ...
       IF (paral%parent) CALL write_shmd(-1,32,infi,tempp)
       IF (wannl%twann) THEN
          CALL mm_dim(mm_go_qm,statusdummy)
          ! McB ... surface hopping stuff ...
          CALL state_select("S0")
          CALL ddipo(taup,c0(:,:,1),cm,c2,sc0,sh02%nst_s0,center)
          CALL forcedr(c0(:,:,1),c2,sc0,rhoe,psi,&
               TAUP,FION0,EIGV,&
               sh02%nst_s0,1,.FALSE.,.TRUE.)
          CALL wc_dos(c0,c2,sh02%nst_s0,center)
          CALL state_select("S1")
          CALL ddipo(taup,c0(:,ns1:,1),cm(:,ns1:),c2(:,ns1:),&
               SC0,sh02%nst_s1,CENTER)
          CALL forcedr(c0(:,ns1:ns1+sh02%nst_s1-1,1),c2(:,ns1:ns1+sh02%nst_s1-1),sc0,rhoe,psi,&
               TAUP,FION1,EIGV,&
               sh02%nst_s1,1,.FALSE.,.TRUE.)
          CALL wc_dos(c0,c2,sh02%nst_s1,center)
          ! .................................
          CALL mm_dim(mm_go_mm,statusdummy)
          ! McB... cf. elct.inc
          ntmp=crge%n
          crge%n=sh02%nsttot
          CALL zhwwf(2,irec,c0,cm,sh02%nsttot,eigv,taup,velp,taui,iteropt%nfi)
          crge%n=ntmp
       ENDIF
       DEALLOCATE(center,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       CALL mm_dim(mm_go_qm,statusdummy)
       IF (rout1%rhoout) CALL rhopri(c0,tau0,rhoe,psi(:,1),sh02%nsttot,nkpt%nkpnt)
       CALL mm_dim(mm_go_mm,statusdummy)
       ! PRINT ACCUMULATOR
       IF (paral%parent) CALL pacca(ekinc,tempp,ener_com%etot,econs,eham,enose,enosp,&
            ener_com%ecnstr,ener_com%erestr,disa,tcpu,iteropt%nfi,2)
       ! McB ... surface hopping stuff ...
       CALL proja(c0       ,c2       ,sc0,sh02%nst_s0,cnti%iproj)
       CALL proja(c0(1,ns1,1),c2(1,ns1),sc0,sh02%nst_s1,cnti%iproj)
       ! .................................
       CALL csize(c2,sh02%nsttot,gemax,cnorm)
       IF (paral%parent) CALL gsize(fion,gnmax,gnorm)
       IF (paral%parent) CALL finalp(tau0,fion,velp,eigv)
       ! 
    ENDIF ! qmnode
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
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL mm_dim(mm_revert,oldstatus)
    CALL mp_sync(parai%qmmmgrp)
#endif
    RETURN
  END SUBROUTINE mm_mdshop_cp
  ! ==================================================================
  SUBROUTINE give_scr_mm_mdshop_cp(lmdshop,tag)
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
    IF (cntl%tqmmm)lmdshop=MAX(lmdshop,fpar%kr1*fpar%kr2s*fpar%kr3s)
    IF (cntl%tqmmm)lmdshop=MAX(lmdshop,maxsys%nax*maxsys%nsx*3)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_mm_mdshop_cp
  ! ==================================================================

END MODULE mm_mdshop_cp_utils
