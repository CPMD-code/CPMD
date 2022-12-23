MODULE mm_mdshop_bo_utils
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
  USE cnst_dyn,                        ONLY: fhills,&
                                             lmeta
  USE comvel_utils,                    ONLY: comvel
  USE comvelmod,                       ONLY: comvl
  USE coor,                            ONLY: fion,&
                                             tau0,&
                                             taup,&
                                             velp
  USE copot_utils,                     ONLY: copot,&
                                             give_scr_copot
  USE cotr,                            ONLY: cotc0
  USE detdof_utils,                    ONLY: detdof
  USE dispp_utils,                     ONLY: dispp
  USE dynit_utils,                     ONLY: dynit
  USE efld,                            ONLY: extf,&
                                             textfld
  USE ekinpp_utils,                    ONLY: ekinpp
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE extrap_utils,                    ONLY: extrap
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_app,&
                                             fo_info,&
                                             fo_verb
  USE finalp_utils,                    ONLY: finalp
  USE fint,                            ONLY: fint1
  USE forcep_utils,                    ONLY: rhoe_psi_size
  USE forces_diag_utils,               ONLY: give_scr_forces_diag
  USE geofile_utils,                   ONLY: geofile
  USE gsize_utils,                     ONLY: gsize
  USE hubbardu,                        ONLY: hubbu
  USE initrun_driver,                  ONLY: initrun
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE linres,                          ONLY: td01
  USE lr_tddft_utils,                  ONLY: give_scr_lr_tddft
  USE machine,                         ONLY: m_flush,&
                                             m_walltime
  USE mddiag_interaction_p_utils,      ONLY: mddiag_interaction_p
  USE mdshop_bo_utils,                 ONLY: give_scr_mdshop_bo,&
                                             mdshop_bo
  USE meta_colvar_inp_utils,           ONLY: colvar_structure
  USE meta_colvar_utils,               ONLY: meta_colvar
  USE meta_exl_mult_utils,             ONLY: meta_ext_mul
  USE meta_exlagr_methods,             ONLY: meta_extlagr
  USE mm_dim_utils,                    ONLY: mm_dim
  USE mm_dimmod,                       ONLY: clsaabox,&
                                             mm_go_mm,&
                                             mm_go_qm,&
                                             mm_revert,&
                                             mmdim
  USE mm_extrap,                       ONLY: cold, nnow, numcold
  USE mm_forces_diag_utils,            ONLY: mm_forces_diag
  USE mm_input,                        ONLY: cgrest_i,&
                                             clc,&
                                             g96_vel,&
                                             lqmmm,&
                                             rtr_l,&
                                             wp_i
  USE mm_mdmain_utils,                 ONLY: mm_localt
  USE mm_parallel,                     ONLY: gparal
  USE moverho_utils,                   ONLY: give_scr_moverho,&
                                             moverho
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum,&
                                             mp_sync
  USE nlcc,                            ONLY: corel
  USE norm,                            ONLY: gnmax,&
                                             gnorm
  USE nose,                            ONLY: glib
  USE noseng_utils,                    ONLY: noseng
  USE nosepa_utils,                    ONLY: nosepa
  USE noseup_utils,                    ONLY: noseup
  USE nospinit_utils,                  ONLY: nospinit
  USE parac,                           ONLY: parai,&
                                             paral
  USE phfac_utils,                     ONLY: phfac
  USE poin,                            ONLY: potr,&
                                             rhoo
  USE posupi_utils,                    ONLY: posupi
  USE printave_utils,                  ONLY: paccc
  USE printp_utils,                    ONLY: printp,&
                                             printp2
  USE prmem_utils,                     ONLY: prmem
  USE proppt_utils,                    ONLY: give_scr_propcal
  USE pslo,                            ONLY: pslo_com
  USE puttau_utils,                    ONLY: taucl
  USE rattle_utils,                    ONLY: rattle
  USE resetac_utils,                   ONLY: resetac
  USE response_pmod,                   ONLY: dmbi
  USE rhoofr_utils,                    ONLY: give_scr_rhoofr
  USE rhopri_utils,                    ONLY: give_scr_rhopri,&
                                             rhopri
  USE rinitwf_utils,                   ONLY: give_scr_rinitwf
  USE rinvel_utils,                    ONLY: rinvel,&
                                             rvscal
  USE rk4ov_utils,                     ONLY: rk4ov_new,&
                                             rk4ov_old
  USE rmas,                            ONLY: rmass
  USE rnlsm_utils,                     ONLY: give_scr_rnlsm
  USE ropt,                            ONLY: infi,&
                                             iteropt,&
                                             ropt_mod
  USE rotvel_utils,                    ONLY: rotvel
  USE rscvp_utils,                     ONLY: rscvp
  USE sample_utils,                    ONLY: sample_go,&
                                             sample_wait
  USE setirec_utils,                   ONLY: read_irec,&
                                             write_irec
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
  USE shop_rest_2,                     ONLY: c0old,&
                                             c0old2
  USE soft,                            ONLY: soft_com
  USE spin,                            ONLY: clsd
  USE store_types,                     ONLY: &
       cprint, irec_ac, irec_nop1, irec_nop2, irec_nop3, irec_nop4, irec_vel, &
       restart1, rout1
  USE system,                          ONLY: &
       cnti, cntl, cntr, fpar, maxsys, nacc, ncpw, nkpt, restf
  USE testex_utils,                    ONLY: testex
  USE teststore_utils,                 ONLY: teststore
  USE totstr_utils,                    ONLY: totstr
  USE tpar,                            ONLY: dt_ions
  USE tst2min_utils,                   ONLY: tst2min
  USE velupi_utils,                    ONLY: velupi
  USE wrener_utils,                    ONLY: wrprint_md
  USE wrgeo_utils,                     ONLY: wrgeof
  USE wv30_utils,                      ONLY: zhwwf
  USE zeroing_utils,                   ONLY: zeroing
#include "sizeof.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mm_mdshop_bo
  PUBLIC :: give_scr_mm_mdshop_bo

CONTAINS

  ! ==================================================================
  SUBROUTINE mm_mdshop_bo(c0,cm,c1,c2,sc0,vpp,gamx,gamy)
    ! ==--------------------------------------------------------------==
    ! McB ... surface hopping stuff ...
    ! McB ........................................................
    ! Qmmm


    ! Qmmm locals
    COMPLEX(real_8)                          :: c0(:,:,:), cm(ncpw%ngw,*), &
                                                c1(*), c2(nkpt%ngwk,*), &
                                                sc0(ncpw%ngw,*)
    REAL(real_8)                             :: vpp(*), gamx(*), gamy(*)

    CHARACTER(len=10)                        :: prch
    CHARACTER(len=100)                       :: filen, filensh
    CHARACTER(len=30)                        :: tag
    COMPLEX(real_8), ALLOCATABLE             :: psi(:,:), psi1(:,:)
    INTEGER :: cm_i, cm_j, i, ia, ifcalc, ig, il_psi_1d, il_psi_2d, &
      il_rhoe_1d, il_rhoe_2d, in, irec(100), is, istate, itemp, j, k, lscr, &
      ncoef, nmin, nmm, nnx, ns1, nstate, ntmp, nx, nxs0, nxs1, out, save_nfi
    LOGICAL                                  :: ferror, lmetares, oldstatus, &
                                                statusdummy
    REAL(real_8) :: c12, c21, dcoupl(2,2), ddet, dei(2), delt, detold2, disa, &
      dt, dtcoef, dummy, e(2), econs, ei(2), ekin1, ekin2, ekincp, ekinh1, &
      ekinh2, ekinp, enose, enosp, eold2(2), lmio(3), tcpu, temp1, temp2, &
      tempp, time1, time2, vcmio(4)
    REAL(real_8), ALLOCATABLE :: eigv(:), eigv1(:), fion0(:,:,:), &
      fion1(:,:,:), rhoe(:,:), rhoe1(:,:), rinp(:), rm1(:), scr(:), &
      taui(:,:,:), tauio(:,:), taur(:,:,:), vpp1(:)

! TODO refactor these arrays 
! Variables
! META DYNAMICS
! McB ... surface hopping stuff ...

    COMMON /ekin/ekinp
    LOGICAL :: lexist
    ! RHOO_1D is temporary alias for RHOO to allocate through MEMORY90
    ! TODO refactor RHOO 
    REAL(real_8), ALLOCATABLE :: rhoo_1d(:)
    ! McB ............................................................
    INTEGER :: ierr
    CHARACTER(*),PARAMETER :: procedureN='mm_mdshop_bo'

    ! ==================================================================
    IF (cntl%tddft.AND.cntl%tresponse) CALL stopgm("MDDIAG",&
         "cntl%tddft.AND.cntl%tresponse NOT POSSIBLE",& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
#if defined (__GROMOS)
    time1 =m_walltime()
    CALL mm_dim(mm_go_mm,oldstatus)
    ! 
    nstate=crge%n
    IF (paral%qmnode) THEN
       IF (textfld)THEN
          ALLOCATE(extf(fpar%kr1*fpar%kr2s*fpar%kr3s),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          CALL zeroing(extf)!,kr1*kr2s*kr3s)
       ENDIF
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

       ! McB ... surface hopping stuff ...
       ! Memory for surface hopping
       ALLOCATE(vpp1(nkpt%ngwk),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ! McB ............................................................

       nacc = 22
       iteropt%nfi  = 0
       ! MODENS=.FALSE.
       ! ENGPRI=.FALSE.
       ! CALSTE=cntl%tpres

       IF (cntl%tqmmm.AND.pslo_com%tivan) THEN
          CALL stopgm('MM_MDSHOP_BO','Qmmm WITH VANDERBILT UNSUPPORTED',& 
               __LINE__,__FILE__)
       ENDIF

       ! CALL MEMORY(IP_EIGV,NLSD*N*NKPTS,'EIGV')

       ! McB ... surface hopping stuff ...
       ! 
       ! Initialize variables for Surface Hopping
       CALL zeroing(sh03%coupl)!,2*2)
       CALL zeroing(sh03%couplold)!,2*2)
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
       ALLOCATE(eigv(nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)

       IF (cntl%tharm.AND.paral%parent) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,A)') ' HARMONIC REFERENCE SYSTEM INTEGRATION',&
               ' ONLY POSSIBLE WITH EQUAL OCCUPATION NUMBERS'
          CALL stopgm('MM_MDSHOP_BO',' ',& 
               __LINE__,__FILE__)
       ENDIF
       ! McB ............................................................
       ! ==------------------------------------------------------------==
       ! Extrapolation
       IF (cntl%textrap) THEN
          ALLOCATE(cold(nkpt%ngwk,crge%n,nkpt%nkpnt,cnti%mextra*nkpt%nkpts/nkpt%nkpnt),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
       ! ==------------------------------------------------------------==
       ! SCR ALLOCATION AND PARTITION (SCRATCH ARRAY).
       CALL rhoe_psi_size(il_psi_1d=il_psi_1d,il_psi_2d=il_psi_2d,&
            il_rhoe_1d=il_rhoe_1d,il_rhoe_2d=il_rhoe_2d)
       ALLOCATE(rhoe(il_rhoe_1d,il_rhoe_2d),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ! McB ... surface hopping stuff ...
       ALLOCATE(rhoe1(il_rhoe_1d,il_rhoe_2d),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ! McB ............................................................
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
       il_psi_1d=nmm*il_psi_1d
       ALLOCATE(psi(il_psi_1d,il_psi_2d),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ! McB ... surface hopping stuff ...
       ALLOCATE(psi1(il_psi_1d,il_psi_2d),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ! McB ............................................................
       IF ( lqmmm%qmmm ) THEN
          CALL give_scr_mm_mdshop_bo(lscr,tag)
       ELSE
          CALL give_scr_mdshop_bo(lscr,tag)
       ENDIF
       ALLOCATE(scr(lscr),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ELSE
       nmin=10
       ALLOCATE(eigv(nmin),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ! McB ... surface hopping stuff ...
       ALLOCATE(eigv1(nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ! McB ............................................................
       CALL give_scr_mdshop_bo(lscr,tag)
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
    CALL mm_dim(mm_go_mm,statusdummy)
    ALLOCATE(fion(3,maxsys%nax,maxsys%nsx),STAT=ierr)
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
    CALL zeroing(fion0)!,3*maxsys%nax*maxsys%nsx)
    CALL zeroing(fion1)!,3*maxsys%nax*maxsys%nsx)
    ! .................................

    ! show sizes of the individual MPI-threads.
    ! use barriers so that the output is not garbled.
    CALL mp_sync(parai%qmmmgrp)
    IF (lqmmm%qmmm_verbose) THEN
       DO k=0,parai%nproc
          IF (gparal%mmparent) THEN
             IF (paral%io_parent)&
                  WRITE(prch,'(a10)') ' MMPARENT '
          ELSE
             IF (paral%io_parent)&
                  WRITE(prch,'(a5,i5)')' QM-N',parai%me
          ENDIF
          IF (k.EQ.parai%me) CALL prmem(prch)
          CALL mp_sync(parai%qmmmgrp)
       ENDDO
    ENDIF
    ! 
    IF (paral%io_parent.AND.cgrest_i%n_cg.GE.1)&
         CALL fileopen(40,'restrain.out',fo_app,ferror)
    IF (paral%io_parent.AND.cntl%tqmmm.AND.cntl%tddft)&
         CALL fileopen(41,'ENERGIES_EX',fo_app,ferror)
    ! ==--------------------------------------------------------------==
    IF (paral%qmnode ) THEN
       IF (cprint%iprint_step.EQ.0) cprint%iprint_step=cnti%nomore+1
       ! ==--------------------------------------------------------------==
       ! Set IREC to restart file.
       CALL read_irec(irec)
       ! INITIALIZATION OF WAVEFUNCTION AND COORDINATES
       ! McB ... ! remember: mddiag.F reads velocities onto C2 ! ....
       CALL initrun(irec,c0,cm,sc0,rhoe,psi,eigv)
       ! mb-bug
       IF (paral%parent) CALL dcopy(3*maxsys%nax*maxsys%nsx,taup,1,taui,1)
       ! mb-bug
       ! ==--------------------------------------------------------------==
       ! TIME STEP FUNCTIONS
       CALL dynit(ekincp,ekin1,ekin2,temp1,temp2,ekinh1,ekinh2)
       ekinp=0.0_real_8  ! McB: cf. META_EXT..()
       ! PARAMETERS FOR THE NOSE-HOOVER THERMOSTATS
       IF (cntl%tnosep.AND.paral%parent) CALL nosepa(1,1)
       ! Dont symmetrize density
       cntl%tsymrho=.FALSE.
       ! ..Make sure TKFULL=.TRUE
       IF (tkpts%tkpnt.AND.(.NOT.tkpts%tkfull)) THEN
          IF (paral%io_parent)&
               WRITE(6,*)&
               ' ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          IF (paral%io_parent)&
               WRITE(6,*)&
               ' WARNING! USE KPOINTS FULL  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          IF (paral%io_parent)&
               WRITE(6,*)&
               ' ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       ENDIF

       ! ==--------------------------------------------------------------==
       ! 28/3/03 added
       CALL mm_dim(mm_go_qm,statusdummy)
       CALL phfac(tau0)
       IF (corel%tinlc) CALL copot(rhoe,psi,ropt_mod%calste)
       CALL mm_dim(mm_go_mm,statusdummy)
       CALL zeroing(cm(:,1:sh02%nsttot))!,ngw*sh02%nsttot)
       IF (.NOT.lqmmm%qmmm_reflex) THEN
          CALL mm_translate_qmmm(tau0,c0,cm,sh02%nsttot)
       ENDIF
       ! 28/3/03

       ! INITIALIZE VELOCITIES
       CALL mm_dim(mm_go_mm,statusdummy)
       IF (paral%parent) CALL detdof(tau0,taur)

       ! INITIALIZE METADYNAMICS VARIABLES used also for 
       ! SHOOTING from SADDLE POINT with RANDOM VELOCITIES
       IF (paral%parent .AND. (lmeta%lcolvardyn .OR. lmeta%lsadpnt)) THEN
          CALL colvar_structure(tau0,taur)
       ENDIF

       IF ( (irec(irec_vel).EQ.0).AND.&
            (.NOT.restart1%rgeo).AND.(.NOT.rtr_l%restart_traj) ) THEN
          ener_com%ecnstr = 0.0_real_8
          ener_com%erestr = 0.0_real_8
          CALL rinvel(velp,c2,crge%n)
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
       IF (clc%classical) CALL zeroing(cm(:,1:sh02%nsttot))!,ngw*sh02%nsttot)
       ! COMPUTE THE IONIC TEMPERATURE TEMPP
       IF (paral%parent) THEN
          CALL ekinpp(ekinp,velp)
          tempp=ekinp*factem*2._real_8/glib
       ENDIF
       ! RESET ACCUMULATORS
       IF (paral%parent.AND.irec(irec_ac).EQ.0)&
            CALL resetac(tau0,taui,iteropt%nfi)
       ! 
       CALL write_irec(irec)
       ! INITIALIZE FORCES
       IF (paral%parent) THEN
          IF (paral%io_parent)&
               WRITE(6,'(1X,64("="))')
          IF (paral%io_parent)&
               WRITE(6,'(1X,"==",T25,A,T64,"==")') 'FORCES INITIALIZATION'
          IF (paral%io_parent)&
               WRITE(6,'(1X,64("="))')
       ENDIF
       CALL mm_dim(mm_go_qm,statusdummy)
       IF (cntl%tdiag) THEN
          IF (cntl%tlanc) nx=1
          IF (cntl%tdavi) nx=nkpt%ngwk*cnti%ndavv*nkpt%nkpnt+1
          IF (cntl%diis)  nx=((nkpt%ngwk*crge%n+8)*cnti%mdiis*nkpt%nkpnt)/4
       ELSEIF (cntl%tsde) THEN
          nx=1
       ELSEIF (cntl%diis) THEN
          nx=(nkpt%ngwk*crge%n+8)*cnti%mdiis/2+4
       ELSEIF (cntl%pcg) THEN
          nx=1
       ENDIF
       nxs0=(nkpt%ngwk*sh02%nst_s0+8)*cnti%mdiis/2+4
       nxs1=(nkpt%ngwk*sh02%nst_s1+8)*cnti%mdiis/2+4
       ifcalc=0

       ! Initialize Metadynamics contributions
       IF (lmeta%lcolvardyn .AND. lmeta%lextlagrange) THEN
          CALL mm_dim(mm_go_mm,statusdummy)
          lmetares= .FALSE.
          disa=0._real_8
          IF (lmeta%tmulti) THEN
             CALL meta_ext_mul(tau0,velp,taur,&
                  .FALSE.,lmetares,.FALSE.,0._real_8,ekinp)
          ELSE

             CALL meta_extlagr(tau0,velp,taur,&
                  .FALSE.,lmetares,.FALSE.,0._real_8,ekinp)
          ENDIF
          CALL mm_dim(mm_revert,statusdummy)
       ENDIF
    ENDIF ! (qmnode)

    CALL mm_dim(mm_go_mm,statusdummy)
    IF ( gparal%mmparent ) THEN
       CALL mm_write_gromos_coord('CRD_INI.g96',tau0,velp,maxsys%nax,maxsys%nsx)
    ENDIF
    CALL mp_sync(parai%qmmmgrp)

    ! McB ... surface hopping stuff ...
    ! ...   initial forces on s0: FION0
    CALL state_select("S0")
    nx=nxs0
    ! ...  convert 1D index NX for CM0 to 2D indices for CM to get rid of CM0      
    cm_i = MOD(nx-1,ncpw%ngw) + 1
    cm_j = (nx-1)/ncpw%ngw + 1
    CALL mm_forces_diag(sh02%nst_s0,c0(:,:,1),c1,c2,cm,sc0,&
         CM(cm_i,cm_j),&
         VPP,EIGV,RHOE,PSI,&
         TAU0,VELP,TAUI,FION0,IFCALC,&
         IREC,.TRUE.,.TRUE.)
    e(1)=ener_com%etot

    ! ...   initial forces on s1 ...
    CALL state_select("S1")
    ifcalc=0
    nx=nxs1
    ! ...  convert 1D index NX for CM1 to 2D indices for CM to get rid of CM1
    cm_i = MOD(nx-1,ncpw%ngw) + 1
    cm_j = (nx-1)/ncpw%ngw + 1 + ns1
    CALL mm_forces_diag(sh02%nst_s1,c0(:,ns1:,1),c1,c2(1,ns1),cm,sc0,&
         CM(cm_i, cm_j),&
         VPP1,&
         EIGV1,RHOE1,PSI1,&
         TAU0,VELP,TAUI,FION1,IFCALC,&
         irec,.TRUE.,.TRUE.)

    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE (6,*) 'mdshop_bo: ',' EADDSH = ',sh02%eaddsh
    ENDIF

    e(2)=ener_com%etot
    e(2)=e(2)+sh02%eaddsh

    ! McB ............................................................
    ! 
    IF (paral%qmnode) THEN
       ! switch on info printing for the lin resp part
       IF (cntl%tresponse) dmbi%inter_pt_firstcall = .TRUE.

       CALL mm_dim(mm_go_mm,statusdummy)
       CALL dcopy(nnx,rin0,1,rm1,1)
       ! INITIALIZE THERMOSTATS
       IF (paral%parent) THEN
          itemp=irec(irec_nop1)+irec(irec_nop2)+irec(irec_nop3)&
               +irec(irec_nop4)
          IF (cntl%tnosep .AND. itemp.EQ.0) CALL nospinit(1)
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
          ! ...      DEFINE INITIAL CONDITIONS FOR INTEGRATION OF STATE POPULATIONS
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
                IF (paral%io_parent)&
                     WRITE(6,*) 'RESH: ',sh03%pop(i)
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
               WRITE(6,'(5X,'' COUPL(,)'',2F20.10)') sh03%coupl(1,1),sh03%coupl(1,2)
          IF (paral%io_parent)&
               WRITE(6,'(5X,'' COUPL(,)'',2F20.10)') sh03%coupl(2,1),sh03%coupl(2,2)
          IF (paral%io_parent)&
               WRITE(6,'(5X,'' E1, E2  '',2F20.10)') sh03%ec(1),sh03%ec(2)
          IF (paral%io_parent)&
               WRITE(6,'(5X,'' E1OLD, E2OLD '',2F20.10)') sh03%eold(1),sh03%eold(2)
       ENDIF
       CALL mp_bcast(sh03%pop,SIZE(sh03%pop),parai%source,parai%allgrp)
       CALL mp_bcast(sh03%isurf,parai%source,parai%allgrp)
       CALL decide(e,fion0,fion1,.FALSE.)
       ! McB ........................................................
       ! ==--------------------------------------------------------------==
       ! == END INITIALIZATION                                           ==
       ! ==--------------------------------------------------------------==

       ! McB ... surface hopping stuff ...
       ! zero out some arrays
       IF (.NOT.sh03%tshopres) CALL zeroing(c0old(:,1:sh02%nsttot))!,ngw*sh02%nsttot)
       CALL zeroing(c0old2(:,1:sh02%nsttot))!,ngw*sh02%nsttot)
       sh03%ec(1)=e(1)
       sh03%ec(2)=e(2)
       ! McB ........................................................
       CALL mm_dim(mm_go_mm,statusdummy)
       IF ( teststore(0).AND.cntl%tsampl ) THEN
          ! McB    ... cf. elct.inc
          ntmp=crge%n
          crge%n=sh02%nsttot
          CALL zhwwf(2,irec,c0,c2,sh02%nsttot,eigv,tau0,velp,taui,iteropt%nfi)
          crge%n=ntmp
       ENDIF
       IF (paral%parent) THEN
          ! MEAN SQUARE DISPLACEMENT OF DIFFERENT IONIC SPECIES
          IF (paral%parent) CALL dispp(tau0,taui,disa)
          ! ENERGY OF THE NOSE THERMOSTATS
          CALL noseng(iteropt%nfi,velp,enose,enosp,dummy,1)
          econs=ekinp+ener_com%etot+enose+enosp+ener_com%ecnstr+ener_com%erestr
          time2 =m_walltime()
          tcpu = (time2 - time1)*0.001_real_8
          IF (paral%io_parent)&
               WRITE(6,'(A,T50,F8.2,A8)') ' TIME FOR INITIALIZATION:',&
               tcpu,' SECONDS'
          IF (paral%io_parent)&
               WRITE(6,'(//,1X,64("="))')
          IF (paral%io_parent)&
               WRITE(6,'(1X,"=",T20,A,T65,"=")')&
               'MOLECULAR DYNAMICS SIMULATION'
          IF (paral%io_parent)&
               WRITE(6,'(1X,64("="))')
          ! McBdbg
          ! write (6,*) 'mm_mdshop_bo: ','>>> wrprint_md'
          ! CALL WRPRINT_MD(EIGV,F,AMU,N,TAU0,FION,
          ! &                  0._real_8,TEMPP,ETOT,ECONS,0._real_8,DISA,
          ! &                  TCPU,.FALSE.,NFI,0)
          ! write (6,*) 'mm_mdshop_bo: ','<<< wrprint_md'
          ! McBdbg
          IF ( tshopold ) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(/1X,''MM_MDSHOP_BO:  USING OLD RK4OV! '')')
          ELSE
             IF (paral%io_parent)&
                  WRITE(6,'(/1X,''MM_MDSHOP_BO:  USING NEW RK4OV! '')')
          ENDIF
       ENDIF
    ENDIF ! (qmnode)
    ! ==================================================================
    ! ==          THE BASIC LOOP FOR MOLECULAR DYNAMICS               ==
    ! ==                 USING VELOCITY VERLET                        ==
    ! ==================================================================
    DO infi=1,cnti%nomore
       CALL mm_dim(mm_go_mm,statusdummy)
       CALL mp_sync(parai%qmmmgrp)
       ! call mp_sync(ALLGRP)
       IF (paral%qmnode) THEN
          time1=m_walltime()
          iteropt%nfi=iteropt%nfi+1

          ! McB ... surface hopping stuff ...
          ! ... store values at time (t-dt), (t-2dt) ...
          DO istate=1,sh02%nsttot
             DO ig=1,ncpw%ngw
                c0old2(ig,istate)=c0old(ig,istate)
                c0old(ig,istate)=c0(ig,istate,1)
             ENDDO
          ENDDO
          DO i=1,2
             eold2(i)=sh03%eold(i)
             sh03%eold(i)=e(i)
             DO j=1,2
                sh03%couplold(i,j)=sh03%coupl(i,j)
             ENDDO
          ENDDO
          detold2=sh03%detold
          sh03%detold=sh03%det
          ! McB ............................................................

          comvl%subcom=comvl%tsubcom.AND.MOD(iteropt%nfi-1,comvl%ncomv).EQ.0
          comvl%subrot=comvl%tsubrot.AND.MOD(iteropt%nfi-1,comvl%nrotv).EQ.0
          IF (hubbu%pfrqom.gt.0) THEN
             hubbu%tpom=MOD(iteropt%nfi-1,hubbu%pfrqom).EQ.0
          ELSE
             hubbu%tpom=.False.
          ENDIF
          ! ANNEALING
          CALL anneal(velp,c2,sh02%nsttot,scr)
          CALL berendsen(velp,c2,nstate,scr,0.0_real_8,0.0_real_8)
          ! SUBTRACT CENTER OF MASS VELOCITY
          IF (paral%parent.AND.comvl%subcom) CALL comvel(velp,vcmio,.TRUE.)
          ! SUBTRACT ROTATION AROUND CENTER OF MASS
          IF (paral%parent.AND.comvl%subrot) CALL rotvel(tau0,velp,lmio,tauio,.TRUE.)
          ! UPDATE NOSE THERMOSTATS
          CALL noseup(velp,c2,sh02%nsttot,1)
          ! UPDATE VELOCITIES
          IF (paral%parent) CALL velupi(velp,fion,1)
          ! UPDATE POSITIONS
          ener_com%ecnstr = 0.0_real_8
          ener_com%erestr = 0.0_real_8
          IF (paral%parent) THEN
             CALL posupi(tau0,taup,velp)
             CALL mm_solv_const('SHAKE',dt_ions,taup,velp,tau0)
             IF (cotc0%mcnstr.NE.0) CALL cpmdshake(tau0,taup,velp)
          ENDIF
       ENDIF! (qmnode)
       ! 
       CALL mp_bcast(taup,SIZE(taup),parai%source,parai%qmmmgrp)
       CALL zeroing(cm(:,1:sh02%nsttot))!,nkpt%ngwk*sh02%nsttot)
       IF (.NOT.lqmmm%qmmm_reflex) THEN
          CALL mm_translate_qmmm(taup,c0,cm,sh02%nsttot)
       ENDIF
       IF (paral%qmnode) THEN
          CALL mm_dim(mm_go_qm,statusdummy)
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
          ropt_mod%calste=cntl%tpres.AND.MOD(infi,cnti%npres).EQ.0
          IF (cntl%textrap) THEN
             ! Extrapolate wavefunctions
             CALL extrapwf(infi,c0,scr,cold,nnow,numcold,sh02%nsttot,cnti%mextra)
          ENDIF
          CALL mm_dim(mm_go_mm,statusdummy)
       ENDIF! (qmnode)
       ! 
       CALL m_flush(6)
       IF (cntl%tlanc) nx=1
       IF (cntl%tdavi) nx=cnti%ndavv*nkpt%nkpnt+1
       CALL mp_sync(parai%qmmmgrp)

       ! RESPONSE calculation
       IF (cntl%tresponse) THEN
          ! McB
          CALL stopgm("MM_MDSHOP_BO",&
               "cntl%tshop.AND.cntl%tresponse NOT POSSIBLE",& 
               __LINE__,__FILE__)
          ! McB

          ! save the cntl%md step number
          save_nfi=iteropt%nfi
          ! localisation at first iteration + PT
          CALL mddiag_interaction_p(crge%n,c0,c2,cm,sc0,&
               CM(NX,1),VPP,EIGV,RHOE,PSI,&! SCR,LSCR,& !vw doesnt match procedure args
               TAUP,VELP,TAUI,FION,IFCALC,IREC,.TRUE.,.FALSE.)
          ! sets the number of iteration and recover cntl%md step number
          ifcalc=iteropt%nfi
          iteropt%nfi=save_nfi

          ! switch off the info printing
          dmbi%inter_pt_firstcall=.FALSE.
       ELSE

          ! McB ... surface hopping stuff ...
          ! ... select S0 ...
          CALL state_select("S0")
          nx=nxs0
          ! ...      convert 1D index for CM0 to 2D index for CM and get rid of CM0
          cm_i = MOD(nx-1,ncpw%ngw) + 1
          cm_j = (nx-1)/ncpw%ngw + 1
          CALL mm_forces_diag(sh02%nst_s0,c0(:,:,1),c1,c2,cm,&
               SC0,CM(cm_i, cm_j),VPP,EIGV,&
               RHOE,PSI,&
               TAUP,VELP,TAUI,FION0,IFCALC,&
               IREC,.TRUE.,.FALSE.)
          e(1)=ener_com%etot
          ! IF(CALSTE) CALL TOTSTR

          ! ... switch to S1 ...
          CALL state_select("S1")

          ! ... forces ...
          nx=nxs1
          ! ...      convert 1D index for CM0 to 2D index for CM and get rid of CM0
          cm_i = MOD(nx-1,ncpw%ngw) + 1
          cm_j = (nx-1)/ncpw%ngw + 1 + ns1
          CALL mm_forces_diag(sh02%nst_s1,c0(:,ns1:,1),c1,c2(1,ns1),cm,&
               SC0,CM(cm_i, cm_j),VPP1,EIGV1,&
               RHOE1,PSI1,&
               TAUP,VELP,TAUI,FION1,IFCALC,&
               IREC,.TRUE.,.FALSE.)
          e(2)=ener_com%etot
          e(2)=e(2)+sh02%eaddsh

          ! IF(CALSTE) CALL TOTSTR
          ! McB ............................................................

       ENDIF

       CALL mm_dim(mm_go_mm,statusdummy)

       ! ==================================================================
       ! Meta Dynamics of Collective Variables

       IF (paral%qmnode)THEN
          IF (lmeta%lcolvardyn) THEN
             lmetares= .FALSE.

             IF (lmeta%lextlagrange) THEN
                ! Metadynamics with Extended Lagrangian
                IF (lmeta%tmulti) THEN
                   CALL meta_ext_mul(taup,velp,taur,&
                        .FALSE.,lmetares,.FALSE.,0.0_real_8,ekinp)
                ELSE
                   CALL meta_extlagr(taup,velp,taur,&
                        .FALSE.,lmetares,.FALSE.,0.0_real_8,ekinp)
                ENDIF
             ELSE
                ! Time dipendent potential applied directly on the Collective Variables
                CALL meta_colvar(taup,velp,fion,taur,&
                     .FALSE.,lmetares,0.0_real_8,ekinp)
             ENDIF

             IF (paral%parent) THEN
                ! Additional Contribution to FION due to the Metadynamics
                ! (from coupling pot.if extended Lagrangian, directly from V(S,t) if not)
                !$omp parallel do private(IS,IA)
                DO is = 1,ions1%nsp
                   DO ia = 1,ions0%na(is)
                      fion(1,ia,is) = fion(1,ia,is) + fhills(1,ia,is)
                      fion(2,ia,is) = fion(2,ia,is) + fhills(2,ia,is)
                      fion(3,ia,is) = fion(3,ia,is) + fhills(3,ia,is)
                   ENDDO
                ENDDO
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
                IF ((paral%parent).AND.paral%io_parent)&
                     CALL fileclose(3)
                GOTO   99999
             ENDIF
          ENDIF
       ENDIF

       ! ==================================================================
       IF (paral%qmnode) THEN
          ! McB ... surface hopping stuff ...
          ! -PARALLEL
          CALL decide(e,fion0,fion1,.FALSE.)
          ! -ENDPARALLEL        
          ! McB ............................................................
          IF (ropt_mod%calste) CALL totstr
          ! ==================================================================
          ! Damped Dynamics
          CALL dampdyn(velp,fion,cm,c2,nstate,scr(1),scr(10))
          ! ==================================================================
          ! FINAL UPDATE FOR VELOCITIES
          ener_com%ecnstr = 0.0_real_8
          ener_com%erestr = 0.0_real_8
          IF (paral%parent) THEN
             CALL velupi(velp,fion,1)
             IF (lqmmm%qmmm_reflex) CALL mm_qm_boundary(taup,velp)! cmb
             CALL mm_solv_const('RATTL',dt_ions,taup,velp,tau0)
             CALL rattle(taup,velp)
          ENDIF
          IF (paral%parent) CALL geofile(taup,velp,'WRITE')
          ! COMPUTE THE IONIC TEMPERATURE TEMPP
          IF (paral%parent) THEN
             CALL ekinpp(ekinp,velp)
             tempp=ekinp*factem*2._real_8/glib
             ! calculate local kinetic temperature in QM and MM subsystem
             ! and write to "QM_TEMP"
             IF (wp_i%nfi_lt.GT.0) THEN
                CALL mm_localt(velp,rmass%pma,factem,tempp,glib,iteropt%nfi,wp_i%nfi_lt)
             ENDIF
          ENDIF
          ! IONIC TEMPERATURE CONTROL
          IF (paral%parent) CALL rscvp(temp1,temp2,tempp,velp)
          ! UPDATE NOSE THERMOSTATS
          CALL noseup(velp,c2,sh02%nsttot,1)
          ! SUBTRACT CENTER OF MASS VELOCITY
          IF (paral%parent.AND.comvl%subcom) CALL comvel(velp,vcmio,.FALSE.)
          ! SUBTRACT ROTATION AROUND CENTER OF MASS
          IF (paral%parent.AND.comvl%subrot) CALL rotvel(tau0,velp,lmio,tauio,.FALSE.&
               )
          CALL berendsen(velp,c2,nstate,scr,0.0_real_8,0.0_real_8)
          ! ANNEALING
          CALL anneal(velp,c2,sh02%nsttot,scr)
          IF (paral%parent) THEN
             CALL ekinpp(ekinp,velp)
             tempp=ekinp*factem*2._real_8/glib
          ENDIF
          ! MEAN SQUARE DISPLACEMENT OF DIFFERENT IONIC SPECIES
          IF (paral%parent) CALL dispp(taup,taui,disa)
          ! ENERGY OF THE NOSE THERMOSTATS
          IF (paral%parent) CALL noseng(iteropt%nfi,velp,enose,enosp,dummy,1)
          ! CALCULATE PROPERTIES DURING SIMULATION.
          cntl%caldip=cntl%tdipd
          ! McB
          ! call mm_dim(mm_go_qm,statusdummy)
          ! CALL PROPCAL(C0,C2,CM,SC0,TAUP,EIGV,F,AMU,
          ! &                RHOE,PSI,SCR,LSCR,N,NKPNT,NFI,INFI)
          ! call mm_dim(mm_go_mm,statusdummy)
          ! McB
          ! PRINTOUT the evolution of the accumulators every time step
          IF (paral%parent) THEN
             IF (paral%io_parent)&
                  CALL fileclose(3)
             IF (paral%io_parent)&
                  CALL fileopen(3,filen,fo_app,ferror)
             econs=ekinp+ener_com%etot+enose+enosp+ener_com%ecnstr+ener_com%erestr
             time2=m_walltime()
             tcpu=(time2-time1)*0.001_real_8
             CALL wrprint_md(eigv,crge%f,ener_com%amu,nstate,taup,fion,&
                  0._real_8,tempp,ener_com%etot,econs,0._real_8,disa,&
                  tcpu,.FALSE.,iteropt%nfi,infi)
             ! UPDATE ACCUMULATORS
             CALL paccc(tempp,ener_com%etot,econs,enose,enosp,ener_com%ecnstr,ener_com%erestr,&
                  ener_com%ebogo,disa,tcpu,iteropt%nfi,1)
             ! STORE IONIC COORDINATES AND VELOCITIES FOR STATISTICS
             ropt_mod%movie=rout1%mout .AND. MOD(iteropt%nfi-1,cnti%imovie).EQ.0
             ropt_mod%rprint=rout1%rout .AND. MOD(iteropt%nfi-1,cnti%ntraj).EQ.0
             ropt_mod%txyz=rout1%xtout .AND. MOD(iteropt%nfi-1,cnti%ntraj).EQ.0
             ropt_mod%tdcd=rout1%dcout .AND. MOD(iteropt%nfi-1,cnti%ntraj).EQ.0
             IF (MOD(iteropt%nfi-1,cnti%ntraj).EQ.0) THEN
                IF (paral%io_parent)&
                     CALL fileopen(83,'MM_CELL_TRANS',fo_app,ferror)
                IF (paral%io_parent)&
                     WRITE(83,'(I10,3F15.10)') iteropt%nfi,(clsaabox%mm_c_trans(k),k=1,3)
                IF (paral%io_parent)&
                     CALL fileclose(83)
             ENDIF
             CALL printp(taur,taup,velp)
             IF (cprint%twriteforcetrajectory) CALL printp2(taur,taup,velp,fion)
             time2=m_walltime()
             tcpu=(time2-time1)*0.001_real_8
          ENDIF
       ENDIF! qmnode

       IF (.NOT.soft_com%exsoft) THEN
          ! McB ... surface hopping stuff only if WFO succeeded ...
          ! ... compute DELTA(C0)/DELTA(t)=CM(t+dt/2) ...
          ! IF ( INFI.GT.2 ) THEN
          IF ( (infi.GT.1).OR.sh03%tshopres ) THEN
             DO ig=1,ncpw%ngw
                DO istate=1,sh02%nsttot
                   cm(ig,istate)=(c0(ig,istate,1)-c0old2(ig,istate))/&
                        (2.0_real_8*delt)
                ENDDO
             ENDDO
          ENDIF

          CALL s0_s1_overlap(c0old,cm,sh03%det,c12,c21)
          sh03%coupl(1,2)=c12
          sh03%coupl(2,1)=c21
          sh03%ec(1)=e(1)
          sh03%ec(2)=e(2)

          IF ( paral%parent ) THEN
             IF (paral%io_parent)&
                  WRITE(6,*)'DET',sh03%det
             IF (paral%io_parent)&
                  WRITE(6,*)'DELTA(DET) ',(sh03%det-detold2)/(2._real_8*delt)
             IF (paral%io_parent)&
                  WRITE(6,*)'analytic:  ',sh03%couplold(2,1)+sh03%couplold(1,2),&
                  ' delta:' ,(sh03%det-detold2)/(2._real_8*delt)-&
                  (sh03%couplold(2,1)+sh03%couplold(1,2))
             IF (paral%io_parent)&
                  WRITE(6,*)'D_21 + D_12',sh03%coupl(2,1)+sh03%coupl(1,2)
             IF (paral%io_parent)&
                  WRITE(6,*)'D_21, D_12',sh03%coupl(2,1),sh03%coupl(1,2)
             IF (paral%io_parent)&
                  WRITE(6,*)'DELTA(E)=',sh03%ec(2)-sh03%ec(1)
          ENDIF
          IF ( (infi.GT.2).OR.sh03%tshopres ) THEN
             ddet=(sh03%det-sh03%detold)/REAL(ncoef,kind=real_8)
             dcoupl(1,2)=(sh03%coupl(1,2)-sh03%couplold(1,2))/REAL(ncoef,kind=real_8)
             dcoupl(2,1)=(sh03%coupl(2,1)-sh03%couplold(2,1))/REAL(ncoef,kind=real_8)
             dei(1)=(sh03%eold(1)-eold2(1))/REAL(ncoef,kind=real_8)
             dei(2)=(sh03%eold(2)-eold2(2))/REAL(ncoef,kind=real_8)
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
                           sh03%det,ddet)
                   ELSE
                      CALL rk4ov_new(dtcoef,sh03%pop,sh03%coupl,dcoupl,ei,dei,&
                           sh03%det,ddet)
                   ENDIF
                ENDIF
             ENDDO
             ei(1)=eold2(1)+dei(1)*REAL(ncoef,kind=real_8)
             ei(2)=eold2(2)+dei(2)*REAL(ncoef,kind=real_8)
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
          CALL mp_bcast_byte(prob1, size_in_bytes_of(prob1),parai%source,parai%allgrp)! /PROB1/
          CALL decide(e,fion0,fion1,.TRUE.)
          ! McB    ............................................................
          IF (paral%parent) THEN
             CALL write_shmd(1,32,infi,tempp)
          ENDIF
       ENDIF

       IF (paral%qmnode) THEN
          IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)
          IF (infi.EQ.cnti%nomore) soft_com%exsoft=.TRUE.
          CALL mm_dim(mm_go_mm,statusdummy)
       ENDIF                  ! qmnode

       in=0
       out=0
       IF (soft_com%exsoft) in=1
       CALL mp_sum(in,out,parai%qmmmgrp)
       IF (out.NE.0) soft_com%exsoft=.TRUE.

       IF (paral%qmnode) THEN
          IF (teststore(iteropt%nfi).OR.soft_com%exsoft) THEN
             ! McB    ... cf. elct.inc
             ntmp=crge%n
             crge%n=sh02%nsttot
             CALL zhwwf(2,irec,c0,c2,sh02%nsttot,eigv,taup,velp,taui,iteropt%nfi)
             crge%n=ntmp
             ! McB    ... write surface hopping restart file: 
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

          IF (soft_com%exsoft .AND.lmeta%lcolvardyn) THEN
             lmetares= .TRUE.
             IF (lmeta%lextlagrange) THEN
                ! Metadynamics with Extended Lagrangian
                IF (lmeta%tmulti) THEN
                   CALL meta_ext_mul(taup,velp,taur,&
                        .FALSE.,lmetares,.FALSE.,0.0_real_8,ekinp)
                ELSE
                   CALL meta_extlagr(taup,velp,taur,&
                        .FALSE.,lmetares,.FALSE.,0.0_real_8,ekinp)
                ENDIF
             ELSE
                ! Time dependent potential applied directly on the Collective Variables
                CALL meta_colvar(taup,velp,fion,taur,&
                     .FALSE.,lmetares,0.0_real_8,ekinp)
             ENDIF
          ENDIF
          ! temperature ramping
          CALL tempramp(temp1,temp2)
       ENDIF                   ! qmnode

       IF (paral%qmnode) THEN
          ! UPDATE IONIC POSITIONS
          CALL dcopy(3*maxsys%nax*maxsys%nsx,taup(1,1,1),1,tau0(1,1,1),1)
          ! UPDATE DENSITY
          CALL extrap(nnx,andr2%alxmix,rm1,rin0,rinp)
          CALL dcopy(nnx,rin0,1,rm1,1)
          CALL dcopy(nnx,rinp,1,rin0,1)
       ENDIF! (qmnode)

       ! STOP THE RUN IF THE USER HAS SET THE SIGNAL 30
       IF (soft_com%exsoft) GOTO 100

       ! ==================================================================
       ! ==     END OF MAIN LOOP                                         ==
       ! ==================================================================
    ENDDO
    IF (paral%qmnode) THEN
       IF (paral%parent) THEN
          IF (paral%io_parent)&
               WRITE(6,'(1X,64("="))')
          IF (paral%io_parent)&
               WRITE(6,'(1X,"=",T17,A,T65,"=")')&
               'END OF MOLECULAR DYNAMICS SIMULATION'
          IF (paral%io_parent)&
               WRITE(6,'(1X,64("="),/,/)')
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
100 CONTINUE
    CALL mm_dim(mm_go_mm,statusdummy)
    IF ( gparal%mmparent ) THEN
       CALL mm_write_gromos_coord('CRD_FIN.g96',taup,velp,maxsys%nax,maxsys%nsx)
    ENDIF
    IF (paral%qmnode) THEN
       ! McB ... surface hopping stuff ...
       IF (paral%parent) CALL write_shmd(-1,32,infi,tempp)
       CALL mm_dim(mm_go_qm,statusdummy)
       IF (rout1%rhoout) CALL rhopri(c0,tau0,rhoe,psi(:,1),sh02%nsttot,nkpt%nkpnt)
       CALL mm_dim(mm_go_mm,statusdummy)
       ! Print accumulators.
       IF (paral%parent) CALL paccc(tempp,ener_com%etot,econs,enose,enosp,&
            ener_com%ecnstr,ener_com%erestr,ener_com%ebogo,disa,tcpu,iteropt%nfi,0)
       IF (paral%parent) CALL gsize(fion,gnmax,gnorm)
       IF (paral%parent) CALL finalp(tau0,fion,velp,eigv)
    ENDIF ! (qmnode)
    IF (cntl%tsampl) THEN
       CALL sample_go
       GOTO 99999
    ENDIF
10000 CONTINUE
    ! ==--------------------------------------------------------------==
    IF (paral%qmnode) THEN
       DEALLOCATE(eigv,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
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
       IF (textfld) THEN
          DEALLOCATE(extf,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
       ENDIF
       IF (cntl%textrap) THEN
          DEALLOCATE(cold,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
       ENDIF
    ELSE
       DEALLOCATE(eigv,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF ! (qmnode)
    ! ==--------------------------------------------------------------==
    CALL mm_dim(mm_revert,oldstatus)
    IF ((paral%parent.AND.cgrest_i%n_cg.GT.1).AND.paral%io_parent)&
         CALL fileclose(40)
    IF ((paral%parent.AND.cntl%tqmmm.AND.cntl%tddft).AND.paral%io_parent)&
         CALL fileclose(41)
    CALL mp_sync(parai%qmmmgrp)
#endif
    RETURN
  END SUBROUTINE mm_mdshop_bo
  ! ==================================================================
  SUBROUTINE give_scr_mm_mdshop_bo(lmmmdshopbo,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lmmmdshopbo
    CHARACTER(len=30)                        :: tag

    INTEGER :: lcalc_alm, lcopot, lforces_diag, lmoverho, lpropcal, lrhoofr, &
      lrhopri, lrinitwf, lrnlsm, ltddft, nstate

    nstate=crge%n
    lrnlsm=0
    lcalc_alm=0
    lcopot=0
    lrhopri=0
    lmoverho=0
    ltddft=0
    CALL give_scr_rinitwf(lrinitwf,tag,nstate)
    IF (pslo_com%tivan) CALL give_scr_rnlsm(lrnlsm,tag,nstate,.FALSE.)
    CALL give_scr_rhoofr(lrhoofr,tag)
    IF (fint1%ttrot) CALL give_scr_calc_alm(lcalc_alm,tag)
    CALL give_scr_forces_diag(lforces_diag,tag,nstate,.TRUE.)
    IF (corel%tinlc) CALL give_scr_copot(lcopot,tag)
    IF (rout1%rhoout) CALL give_scr_rhopri(lrhopri,tag,nstate)
    CALL give_scr_propcal(lpropcal,tag,nstate)
    IF (tmovr) CALL give_scr_moverho(lmoverho,tag)
    IF (cntl%tddft) CALL give_scr_lr_tddft(ltddft,.TRUE.,tag)
    lmmmdshopbo=MAX(lrinitwf,lrnlsm,lrhoofr,lforces_diag,ltddft,&
         lcopot,lcalc_alm,lrhopri,lpropcal,lmoverho)
    IF (cntl%tqmmm) lmmmdshopbo=MAX(lmmmdshopbo,fpar%kr1*fpar%kr2s*fpar%kr3s)
    IF (cntl%tqmmm) lmmmdshopbo=MAX(lmmmdshopbo,maxsys%nax*maxsys%nsx*3)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_mm_mdshop_bo
  ! ==================================================================

END MODULE mm_mdshop_bo_utils
