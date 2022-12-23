MODULE prcpmd_utils
  USE anneal_utils,                    ONLY: anneal,&
                                             berendsen,&
                                             dampdyn,&
                                             tempramp
  USE cnst,                            ONLY: factem
  USE cnst_dyn,                        ONLY: dstrcv,&
                                             dstrmeta,&
                                             fhills,&
                                             fvbound,&
                                             lmeta,&
                                             ltcglobal,&
                                             ncolvar,&
                                             rmeta
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
                                             fo_verb
  USE finalp_utils,                    ONLY: finalp
  USE forcedr_driver,                  ONLY: forcedr
  USE forcedr_utils,                   ONLY: give_scr_forcedr
  USE forcep_utils,                    ONLY: rhoe_psi_size
  USE freqs_utils,                     ONLY: freqs
  USE geofile_utils,                   ONLY: geofile
  USE geq0mod,                         ONLY: geq0
  USE gsize_utils,                     ONLY: gsize
  USE initrun_driver,                  ONLY: initrun
  USE initrun_utils,                   ONLY: give_scr_initrun
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE localize_utils,                  ONLY: localize2
  USE machine,                         ONLY: m_walltime
  USE meta_cell_utils,                 ONLY: delrot_cell,&
                                             meta_cell
  USE meta_colvar_inp_utils,           ONLY: colvar_structure
  USE meta_colvar_utils,               ONLY: meta_colvar
  USE meta_exlagr_methods,             ONLY: meta_extlagr
  USE meta_exlagr_utils,               ONLY: ekincv_global,&
                                             meta_stress
  USE metr,                            ONLY: metr_com
  USE mp_interface,                    ONLY: mp_bcast
  USE newcell_utils,                   ONLY: give_scr_newcell,&
                                             newcell
  USE nlcc,                            ONLY: corel
  USE norm,                            ONLY: cnorm,&
                                             gemax,&
                                             gnmax,&
                                             gnorm
  USE nose,                            ONLY: glib,&
                                             ncdof
  USE ortho_utils,                     ONLY: give_scr_ortho
  USE parac,                           ONLY: parai,&
                                             paral
  USE phfac_utils,                     ONLY: phfac
  USE posupa_utils,                    ONLY: give_scr_posupa,&
                                             posupa
  USE posupi_utils,                    ONLY: posupi
  USE prcp,                            ONLY: prcp_com
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
  USE shake_utils,                     ONLY: cpmdshake
  USE soft,                            ONLY: soft_com
  USE spin,                            ONLY: spin_mod
  USE store_types,                     ONLY: cprint,&
                                             irec_ac,&
                                             irec_vel,&
                                             restart1,&
                                             rout1
  USE str2,                            ONLY: sfion,&
                                             stau0,&
                                             svelp
  USE system,                          ONLY: &
       cnti, cntl, cntr, maxsys, nacc, ncpw, nkpt, parm, restf
  USE testex_utils,                    ONLY: testex
  USE teststore_utils,                 ONLY: teststore
  USE totstr_utils,                    ONLY: totstr
  USE tpar,                            ONLY: dt_ions
  USE utils,                           ONLY: zclean,&
                                             zclean_k
  USE vdwcmod,                         ONLY: vdwl,&
                                             vdwwfl
  USE velupa_utils,                    ONLY: velupa
  USE velupi_utils,                    ONLY: c_to_s,&
                                             s_to_c,&
                                             velupc1,&
                                             velupc2
  USE wannier_print_utils,             ONLY: wannier_print
  USE wrener_utils,                    ONLY: wrprint_pr
  USE wrgeo_utils,                     ONLY: wrgeof
  USE wv30_utils,                      ONLY: zhwwf
  USE zeroing_utils,                   ONLY: zeroing
#include "sizeof.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: prcpmd
  PUBLIC :: give_scr_prcpmd
  !public :: rscvc

CONTAINS

  ! ==================================================================
  SUBROUTINE prcpmd(c0,cm,c2,sc0,gamx,gamy)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8) :: c0(ncpw%ngw,crge%n,1), cm(ncpw%ngw,crge%n), &
      c2(ncpw%ngw,crge%n), sc0(ncpw%ngw,crge%n)
    REAL(real_8)                             :: gamx(*), gamy(*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'prcpmd'

    CHARACTER(len=100)                       :: filen
    CHARACTER(len=30)                        :: tag
    COMPLEX(real_8), ALLOCATABLE             :: psi(:,:)
    INTEGER                                  :: i, ia, ierr, il_psi_1d, &
                                                il_psi_2d, il_rhoe_1d, &
                                                il_rhoe_2d, irec(100), is, j, &
                                                lscr, nstate
    LOGICAL                                  :: ferror, lmetares, lquench, &
                                                resetcv
    REAL(real_8) :: disa, econs, eham, ek_cv, ekin1, ekin2, ekinc, ekincp, &
      ekinh, ekinh1, ekinh2, ekinp, enose, enosp, ff, ha(3,3), hb(3,3), tcpu, &
      temp1, temp2, temph, tempp, time1, time2, vcmio(4)
    REAL(real_8), ALLOCATABLE                :: center(:,:), eigm(:), &
                                                eigv(:), rhoe(:,:), scr(:), &
                                                taui(:,:,:), taur(:,:,:)
    REAL(real_8), EXTERNAL                   :: ddot

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
    ALLOCATE(stau0(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(svelp(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(sfion(3,maxsys%nax,maxsys%nsx),STAT=ierr)
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
    nstate=crge%n
    nacc = 22
    ncdof=9
    ropt_mod%modens=.FALSE.
    ropt_mod%engpri=.FALSE.
    ropt_mod%calste=.TRUE.
    IF (pslo_com%tivan) cntl%nonort=.TRUE.
    IF (cntl%nonort) THEN
       ALLOCATE(eigv(crge%n),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(eigm(crge%n*crge%n),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ELSE
       ! Avoid 'not allocated' runtime error
       ALLOCATE(eigv(1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(eigm(1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
    ENDIF
    IF (cntl%tharm.AND.paral%parent) THEN
       ff=crge%f(1,1)
       DO i=1,crge%n
          IF (ff.NE.crge%f(i,1)) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(A,A)') ' HARMONIC REFERENCE SYSTEM INTEGRATION',&
                  ' ONLY POSSIBLE WITH EQUAL OCCUPATION NUMBERS'
             CALL stopgm('PRCPMD',' ',& 
                  __LINE__,__FILE__)
          ENDIF
       ENDDO
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
    CALL give_scr_prcpmd(lscr,tag)
    ALLOCATE(scr(lscr),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    restf%nfnow=1
    ! ==--------------------------------------------------------------==
    ! TIME STEP FUNCTIONS
    CALL dynit(ekincp,ekin1,ekin2,temp1,temp2,ekinh1,ekinh2)
    ekinc=0.0_real_8    ! McB: cf. META_EXT..()
    ekinp=0.0_real_8    ! McB: cf. META_EXT..()
    ! INITIALIZATION OF WAVEFUNCTION AND COORDINATES
    CALL initrun(irec,c0,cm,sc0,rhoe,psi,eigv)
    CALL mp_bcast(taup,SIZE(taup),parai%io_source,parai%cp_grp)
    CALL dcopy(3*maxsys%nax*maxsys%nsx,taup,1,taui,1)
    ! INITIALIZE WF CENTERS & SPREAD
    IF (vdwl%vdwd.AND..NOT.vdwwfl%trwannc) THEN
       CALL localize2(tau0,c0,c2,sc0,nstate)
    ENDIF
    ! QUENCHING TO BO SURFACE
    IF (cntl%quenchb) THEN
       CALL quenbo(c0(:,:,1),c2,sc0,taur,rhoe,psi)
#if 0
       ! AK: another broken restart (see QM/MM). why??
       CALL zhwwf(2,irec,c0,cm,n,eigv,tau0,velp,taui,iteropt%nfi)
#endif
       cntl%quenchb=.FALSE.
    ENDIF
    IF (pslo_com%tivan) THEN
       IF (cntl%tlsd) THEN
          CALL deort(ncpw%ngw,spin_mod%nsup,eigm,eigv,c0,sc0)
          CALL deort(ncpw%ngw,spin_mod%nsdown,eigm,eigv,c0(1,spin_mod%nsup+1,1),sc0(1,spin_mod%nsup+1))
       ELSE
          CALL deort(ncpw%ngw,crge%n,eigm,eigv,c0,sc0)
       ENDIF
    ENDIF
    ! INITIALIZE VELOCITIES
    ener_com%ecnstr = 0.0_real_8
    ener_com%erestr = 0.0_real_8
    IF (paral%parent) CALL detdof(tau0,taur)
    IF (irec(irec_vel).EQ.0.AND..NOT.restart1%rgeo) THEN
       ener_com%ecnstr = 0.0_real_8
       ener_com%erestr = 0.0_real_8
       CALL rinvel(velp,cm,crge%n)
       IF (paral%parent) CALL taucl(velp)
       IF (paral%parent) CALL rattle(tau0,velp)
       CALL rvscal(velp)
    ELSE
       IF (paral%parent) CALL taucl(velp)
       IF (paral%parent) CALL rattle(tau0,velp)
       IF (cntl%trescale) CALL rvscal(velp)
    ENDIF
    IF (cntl%quenchp) CALL zeroing(velp)!,3*maxsys%nax*maxsys%nsx)
    IF (cntl%quenche) CALL zeroing(cm)!,ngw*n)
    IF (cntl%quenchc) CALL zeroing(metr_com%htvel)!,9)

    ! INITIALIZE METADYNAMICS VARIABLES used also for 
    ! SHOOTING from SADDLE POINT with RANDOM VELOCITIES
    IF (paral%io_parent .AND. (lmeta%lcolvardyn .OR. lmeta%lsadpnt)) THEN
       CALL colvar_structure(tau0,taur)
    ENDIF
    IF (cntl%trevers) THEN
       ! invert electronic, ionic, and cell velocities (for path sampling)
       CALL dscal(2*ncpw%ngw*nstate,-1._real_8,cm,1)
       CALL dscal(3*maxsys%nax*maxsys%nsx,-1._real_8,velp,1)
       CALL dscal(9,-1._real_8,metr_com%htvel,1)
    ENDIF
    ! >>>
    ! make sure the velocities are correctly replicated while using groups
    CALL mp_bcast(velp,SIZE(velp),parai%io_source,parai%cp_grp)
    CALL mp_bcast(cm,ncpw%ngw*nstate,parai%cp_inter_io_source,parai%cp_inter_grp)
    ! <<<
    ! RESET ACCUMULATORS
    IF (paral%parent.AND.irec(irec_ac).EQ.0)&
         CALL resetac(tau0,taui,iteropt%nfi)
    ! INITIALIZE FORCES
    IF (paral%io_parent) THEN
       WRITE(6,'(1X,64("="))')
       WRITE(6,'(1X,"==",T25,A,T64,"==")') 'FORCES INITIALIZATION'
       WRITE(6,'(1X,64("="))')
    ENDIF
    IF (tkpts%tkpnt) THEN
       IF (geq0) CALL zclean_k(c0,crge%n,ncpw%ngw)
    ELSE
       IF (geq0) CALL zclean(c0,crge%n,ncpw%ngw)
    ENDIF
    CALL forcedr(c0(:,:,1),c2,sc0,rhoe,psi,tau0,fion,eigv,&
         crge%n,1,.FALSE.,.TRUE.)

    ! INITIALIZE STRESS CONTRIBUTION DUE TO METADYNAMICS
    CALL zeroing(dstrmeta)!,6)
    CALL zeroing(fvbound)!,6)
    CALL zeroing(dstrcv)!,6)


    ! INITIALIZE FORCES CONTRIBUTION DUE TO METADYNAMICS
    IF (lmeta%lmeta_cell) THEN
       lquench = .FALSE.
       lmetares= .FALSE.
       ! Time dipendent potential applied directly on the Collective Variables
       CALL meta_cell(velp,lquench,lmetares)
    ELSEIF (lmeta%lcolvardyn .AND. lmeta%lextlagrange) THEN
       lquench = .FALSE.
       lmetares= .FALSE.
       resetcv = .FALSE.
       CALL meta_extlagr(tau0,velp,taur,lquench,lmetares,&
            resetcv,ekinc,ekinp)
       CALL meta_stress(tau0)
    ENDIF

    CALL totstr
    CALL freqs(crge%n,.FALSE.)
    ! Check orthogonality condition for wavefunction velocities
    CALL rortv(c0,cm,c2,sc0,gamy,crge%n)

    IF (paral%io_parent) THEN
       CALL wrgeof(tau0,fion)
       filen='ENERGIES'
       CALL fileopen(3,filen,fo_app+fo_verb,ferror)
       WRITE(6,'(/,A)') '  CELL PARAMETERS '
       DO i=1,3
          WRITE(6,'(3(1X,F12.8))') (metr_com%ht(i,j),j=1,3)
       ENDDO
       WRITE(6,'(/,A)') '  CELL VELOCITY '
       DO i=1,3
          WRITE(6,'(3(1X,F12.8))') (metr_com%htvel(i,j),j=1,3)
       ENDDO
       WRITE(6,'(/,A)') '  CELL FORCES '
       DO i=1,3
          WRITE(6,'(3(1X,F12.8))') (metr_com%htfor(i,j),j=1,3)
       ENDDO
       WRITE(6,'(1X,64("="))')
       WRITE(6,'(1X,"==",T20,A,T64,"==")')&
            'END OF FORCES INITIALIZATION'
       WRITE(6,'(1X,64("="),/)')
    ENDIF
    CALL write_irec(irec)
    ! ==--------------------------------------------------------------==
    ! == END INITIALIZATION                                           ==
    ! ==--------------------------------------------------------------==
    enose=0.0_real_8
    enosp=0.0_real_8
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
    ALLOCATE(center(4,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    DO infi=1,cnti%nomore
       time1=m_walltime()
       iteropt%nfi=iteropt%nfi+1
       ropt_mod%prteig=cprint%tprint.AND.MOD(iteropt%nfi-1,cprint%iprint_step).EQ.0
       cntl%caldip=cntl%tdipd.AND.MOD(iteropt%nfi-1,cnti%npdip).EQ.0
       comvl%subcom=comvl%tsubcom.AND.MOD(iteropt%nfi-1,comvl%ncomv).EQ.0
       IF (.NOT.paral%parent) ropt_mod%prteig=.FALSE.
       ropt_mod%engpri=cprint%tprint.AND.MOD(iteropt%nfi-1,cprint%iprint_step).EQ.0
       ! ANNEALING
       CALL anneal(velp,cm,crge%n,metr_com%htvel)
       CALL berendsen(velp,cm,nstate,metr_com%htvel,ekinc,ekinh)
       ! SUBTRACT CENTER OF MASS VELOCITY
       IF (paral%parent.AND.comvl%subcom) CALL comvel(velp,vcmio,.TRUE.)
       CALL c_to_s(tau0,stau0)
       CALL c_to_s(velp,svelp)
       CALL c_to_s(fion,sfion)
       ! UPDATE VELOCITIES
       IF (paral%parent) CALL velupc1(svelp,sfion,velp)
       CALL velupa(c0,cm,c2,crge%n,1)
       ! UPDATE POSITIONS
       IF (lmeta%lmeta_cell) THEN
          CALL dcopy(9,metr_com%ht,1,ha,1)
       ENDIF
       ! UPDATE CELL
       CALL daxpy(9,dt_ions,metr_com%htvel,1,metr_com%ht,1)
       IF (lmeta%lmeta_cell) THEN
          CALL dcopy(9,metr_com%ht,1,hb,1)
          CALL delrot_cell(ha,hb)
          CALL dcopy(9,hb,1,metr_com%ht,1)
       ENDIF

       CALL mp_bcast_byte(metr_com, size_in_bytes_of(metr_com),parai%io_source,parai%cp_grp)
       DO i=1,3
          parm%a1(i) = metr_com%ht(1,i)
          parm%a2(i) = metr_com%ht(2,i)
          parm%a3(i) = metr_com%ht(3,i)
       ENDDO
       CALL ihmat(metr_com%ht,metr_com%htm1,parm%omega)
       IF (paral%parent) THEN
          CALL posupi(stau0,taup,svelp)
          CALL dcopy(3*maxsys%nax*maxsys%nsx,taup(1,1,1),1,stau0,1)
       ENDIF
       CALL mp_bcast(stau0,SIZE(stau0),parai%io_source,parai%cp_grp)
       CALL s_to_c(stau0,taup)
       CALL s_to_c(svelp,velp)
       CALL newcell
       ener_com%ecnstr = 0.0_real_8
       ener_com%erestr = 0.0_real_8
       IF (cotc0%mcnstr.NE.0) THEN
          IF(paral%parent) CALL cpmdshake(tau0,taup,velp)
          CALL mp_bcast(taup,SIZE(taup),parai%io_source,parai%cp_grp)
       ENDIF
       CALL phfac(taup)
       IF (corel%tinlc) CALL copot(rhoe,psi,.TRUE.)
       CALL posupa(c0,cm,c2,gamx,crge%n)
       ! ..Dipole moment
       IF (cntl%caldip.OR.(vdwl%vdwd.AND.vdwwfl%twannup)) THEN
          CALL ddipo(taup,c0(:,:,1),cm,c2,sc0,crge%n,center)
          CALL wannier_print(iteropt%nfi,c0(:,:,1),taup,crge%n,psi(:,1),center)
       ENDIF
       ! CALCULATE THE FORCES
       CALL forcedr(c0(:,:,1),c2,sc0,rhoe,psi,taup,fion,eigv,&
            crge%n,1,.FALSE.,.TRUE.)

       IF (lmeta%lmeta_cell) THEN
          lquench = .FALSE.
          lmetares= .FALSE.
          ! Time dipendent potential applied directly on the Collective Variables
          CALL meta_cell(velp,lquench,lmetares)

          IF ((rmeta%tolkin.NE.-1._real_8 .AND. ekinc.GT.rmeta%tolkin)&
               .OR. lquench) THEN
             ! Check for Quench on the BO surface
             cntl%quenchb=.TRUE.
             CALL quenbo(c0(:,:,1),c2,sc0,taur,rhoe,psi)
             CALL zeroing(cm)!,n*nkpt%ngwk)
          ENDIF

       ELSEIF (lmeta%lcolvardyn .AND. lmeta%lextlagrange) THEN
          lquench = .FALSE.
          lmetares= .FALSE.
          CALL meta_extlagr(taup,velp,taur,&
               lquench,lmetares,resetcv,ekinc,ekinp)
          CALL meta_stress(taup)
          IF ((rmeta%tolkin.NE.-1._real_8 .AND. ekinc.GT.rmeta%tolkin)&
               .OR. lquench) THEN
             ! Check for Quench on the BO surface
             CALL quenbo(c0(:,:,1),c2,sc0,taur,rhoe,psi)
             CALL zeroing(cm)!,n*nkpt%ngwk)
             resetcv = .TRUE.
             cntl%quenchb=.FALSE.
          ENDIF
       ELSEIF (cntr%tolkinc.NE.-1._real_8 .AND. ekinc.GT.cntr%tolkinc) THEN
          ! Check for Quench on the BO surface
          CALL quenbo(c0(:,:,1),c2,sc0,taur,rhoe,psi)
          CALL zeroing(cm)!,n*nkpt%ngwk)
          cntl%quenchb=.FALSE.
       ENDIF

       CALL totstr
       ! ==================================================================
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


       ! ==================================================================
       ! Damped Dynamics
       CALL dampdyn(velp,fion,cm,c2,nstate,metr_com%htvel,metr_com%htfor)
       ! ==================================================================
       CALL c_to_s(taup,stau0)
       CALL c_to_s(velp,svelp)
       CALL c_to_s(fion,sfion)

       ! FINAL UPDATE FOR VELOCITIES
       ener_com%ecnstr = 0.0_real_8
       ener_com%erestr = 0.0_real_8
       IF (paral%parent) THEN
          CALL velupc2(svelp,sfion,velp,taur)
          CALL s_to_c(svelp,velp)
          CALL rattle(taup,velp)
       ENDIF
       CALL velupa(c0,cm,c2,crge%n,1)
       CALL rortv(c0,cm,c2,sc0,gamy,crge%n)
       ! SUBTRACT CENTER OF MASS VELOCITY
       IF (paral%parent.AND.comvl%subcom) CALL comvel(velp,vcmio,.FALSE.)
       CALL berendsen(velp,cm,nstate,metr_com%htvel,ekinc,ekinh)
       ! ANNEALING
       CALL anneal(velp,cm,crge%n,metr_com%htvel)
       IF (paral%parent) CALL geofile(taup,velp,'WRITE')
       ! COMPUTE THE IONIC TEMPERATURE TEMPP
       IF (paral%parent) THEN
          CALL ekinpp(ekinp,velp)
          IF (lmeta%lextlagrange.AND. ltcglobal) THEN
             CALL ekincv_global(ek_cv)
             tempp=(ek_cv+ekinp)*factem*2._real_8/(glib+REAL(ncolvar,kind=real_8))
          ELSE
             tempp=ekinp*factem*2._real_8/glib
          ENDIF
          ! IONIC TEMPERATURE CONTROL
          CALL rscvp(temp1,temp2,tempp,velp)
       ENDIF
       ! RESCALE ELECTRONIC VELOCITIES
       IF (cntl%tc) CALL rscve(ekin1,ekin2,ekinc,cntr%ekinw,cm,crge%n,ncpw%ngw)
       ! MEAN SQUARE DISPLACEMENT OF DIFFERENT IONIC SPECIES
       IF (paral%parent) CALL dispp(taup,taui,disa)
       ! KINETIC ENERGY OF THE ELECTRONS
       CALL rekine(cm,crge%n,ekinc)
       ! KINETIC ENERGY OF THE CELL
       ekinh=0.5_real_8*cntr%cmass*ddot(9,metr_com%htvel,1,metr_com%htvel,1)
       IF (paral%parent) THEN
          ! ..CELL TEMPERATURE CONTROL
          CALL rscvc(ekinh,metr_com%htvel,ekinh1,ekinh2)
          CALL ekinpp(ekinp,velp)
          IF (lmeta%lextlagrange.AND. ltcglobal) THEN
             CALL ekincv_global(ek_cv)
             tempp=(ek_cv+ekinp)*factem*2._real_8/(glib+REAL(ncolvar,kind=real_8))
          ELSE
             tempp=ekinp*factem*2._real_8/glib
          ENDIF
          ekinh=0.5_real_8*cntr%cmass*ddot(9,metr_com%htvel,1,metr_com%htvel,1)
          temph = 2._real_8*ekinh*factem/REAL(ncdof,kind=real_8)
          econs=ekinp+ener_com%etot+ener_com%ecnstr+ener_com%erestr+prcp_com%druck*parm%omega
          eham=econs+ekinc+ekinh
          rmeta%eham_hill = eham
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
       IF (teststore(iteropt%nfi).OR.soft_com%exsoft.OR.lmetares) THEN
          CALL zhwwf(2,irec,c0,cm,crge%n,eigv,taup,velp,taui,iteropt%nfi)
       ENDIF
       IF (soft_com%exsoft .AND.lmeta%lcolvardyn) THEN
          lmetares= .TRUE.
          IF (lmeta%lextlagrange) THEN
             ! Metadynamics with Extended Lagrangian
             CALL meta_extlagr(taup,velp,taur,&
                  lquench,lmetares,resetcv,ekinc,ekinp)
          ELSE
             ! Time dipendent potential applied directly on the Collective Variables
             CALL meta_colvar(taup,velp,fion,taur,&
                  lquench,lmetares,ekinc,ekinp)
          ENDIF
       ENDIF
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
    DEALLOCATE(center,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
100 CONTINUE
    IF (rout1%rhoout.AND.(rout1%nrhoout.LE.0))&
         CALL rhopri(c0,tau0,rhoe,psi(:,1),crge%n,nkpt%nkpnt)
    ! PRINT ACCUMULATORS.
    IF (paral%parent) CALL paccb(ekinc,tempp,ener_com%etot,econs,eham,disa,tcpu,temph&
         ,parm%omega,enose,enosp,iteropt%nfi,2)
    CALL proja(c0,c2,sc0,crge%n,cnti%iproj)
    CALL csize(c2,crge%n,gemax,cnorm)
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
    DEALLOCATE(stau0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(svelp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(sfion,STAT=ierr)
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
  END SUBROUTINE prcpmd
  ! ==================================================================
  SUBROUTINE give_scr_prcpmd(lprcpmd,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lprcpmd
    CHARACTER(len=30)                        :: tag

    INTEGER :: lcopot, lddipo, ldeort, lforcedr, linitrun, lnewcell, lortho, &
      lposupa, lquenbo, lrhopri, lrortv, nstate

    nstate=crge%n
    lcopot=0
    lortho=0
    lquenbo=0
    ldeort=0
    lrhopri=0
    lddipo=0
    IF (corel%tinlc) CALL give_scr_copot(lcopot,tag)
    IF (cntl%trane) CALL give_scr_ortho(lortho,tag,nstate)
    IF (cntl%quenchb) CALL give_scr_quenbo(lquenbo,tag)
    IF (pslo_com%tivan) CALL give_scr_deort(ldeort,tag,nstate)
    IF (cntl%tdipd.OR.vdwl%vdwd) CALL give_scr_ddipo(lddipo,tag)
    CALL give_scr_forcedr(lforcedr,tag,nstate,.FALSE.,.TRUE.)
    CALL give_scr_rortv(lrortv,tag,nstate)
    CALL give_scr_newcell(lnewcell,tag)
    CALL give_scr_posupa(lposupa,tag,nstate)
    CALL give_scr_initrun(linitrun,tag)
    IF (rout1%rhoout) CALL give_scr_rhopri(lrhopri,tag,nstate)
    lprcpmd=MAX(lcopot,lortho,lquenbo,ldeort,lddipo,lforcedr,&
         lrortv,lposupa,lrhopri,lnewcell,linitrun)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_prcpmd
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


END MODULE prcpmd_utils
