MODULE mm_mddiag_utils
  USE andp,                            ONLY: rin0,&
                                             rmix,&
                                             rout0
  USE andr,                            ONLY: andr2
  USE anneal_utils,                    ONLY: anneal,&
                                             berendsen,&
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
                                             fo_verb
  USE finalp_utils,                    ONLY: finalp
  USE fint,                            ONLY: fint1
  USE forcep_utils,                    ONLY: rhoe_psi_size
  USE forces_diag_utils,               ONLY: give_scr_forces_diag
  USE geofile_utils,                   ONLY: geofile
  USE gsize_utils,                     ONLY: gsize
  USE hubbardu,                        ONLY: hubbu
  USE initrun_driver,                  ONLY: initrun
  USE initrun_utils,                   ONLY: give_scr_initrun
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE linres,                          ONLY: &
       c_sh, ck, dt_sh, lrsym, shfilen1, shfilen2, shsigma, shsigmaf, td01, tempp_sh, &
       tshi, tshl, v_sh
  USE localize_utils,                  ONLY: localize2
  USE lr_tddft_utils,                  ONLY: give_scr_lr_tddft
  USE machine,                         ONLY: m_flush,&
                                             m_walltime
  USE md_driver,                       ONLY: give_scr_mddiag
  USE mddiag_interaction_p_utils,      ONLY: mddiag_interaction_p
  USE meta_colvar_inp_utils,           ONLY: colvar_structure
  USE meta_colvar_utils,               ONLY: meta_colvar
  USE meta_exl_mult_utils,             ONLY: meta_ext_mul
  USE meta_exlagr_methods,             ONLY: give_scr_meta_extlagr,&
                                             meta_extlagr
  USE mm_dim_utils,                    ONLY: mm_dim
  USE mm_dimmod,                       ONLY: clsaabox,&
                                             mm_go_mm,&
                                             mm_go_qm,&
                                             mm_revert
  USE mm_extrap,                       ONLY: cold, nnow, numcold
  USE mm_forces_diag_utils,            ONLY: mm_forces_diag
  USE mm_forces_prop_utils,            ONLY: mm_forces_prop
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
  USE proppt_utils,                    ONLY: give_scr_propcal,&
                                             propcal
  USE pslo,                            ONLY: pslo_com
  USE puttau_utils,                    ONLY: taucl
  USE rattle_utils,                    ONLY: rattle
  USE resetac_utils,                   ONLY: resetac
  USE response_pmod,                   ONLY: dmbi
  USE rhoofr_c_utils,                  ONLY: rhoofr_c
  USE rhoofr_utils,                    ONLY: give_scr_rhoofr
  USE rhopri_utils,                    ONLY: give_scr_rhopri,&
                                             rhopri
  USE rinitwf_utils,                   ONLY: give_scr_rinitwf
  USE rinvel_utils,                    ONLY: rinvel,&
                                             rvscal
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
  USE sh_tddft_utils,                  ONLY: rscvp_sh,&
                                             tmprd,&
                                             tmpwr
  USE shake_utils,                     ONLY: cpmdshake
  USE soc,                             ONLY: do_tsh_isc,&
                                             do_tsh_isc_lz,&
                                             isc_read_start,&
                                             isc_write_start
  USE soc_types,                       ONLY: &
       cs1, ct1, do_soc_because_sh, ene_sin, ene_tri, facthcm, md_on_sin, &
       md_on_tri, nskip, tausoc
  USE soft,                            ONLY: soft_com
  USE spin,                            ONLY: clsd
  USE store_types,                     ONLY: &
       cprint, irec_ac, irec_nop1, irec_nop2, irec_nop3, irec_nop4, irec_vel, &
       restart1, rout1, store1
  USE system,                          ONLY: &
       cnti, cntl, cntr, fpar, maxsys, nacc, ncpw, nkpt, restf
  USE td_input,                        ONLY: td_prop
  USE td_utils,                        ONLY: getnorm_k,&
                                             load_ex_states,&
                                             tmprd_prop,&
                                             tmpwr_prop
  USE testex_utils,                    ONLY: testex
  USE teststore_utils,                 ONLY: teststore
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE totstr_utils,                    ONLY: totstr
  USE tpar,                            ONLY: dt_ions
  USE tst2min_utils,                   ONLY: tst2min
  USE vdwcmod,                         ONLY: vdwl,&
                                             vdwwfl
  USE velupi_utils,                    ONLY: velupi
  USE wrener_utils,                    ONLY: wrprint_md
  USE wv30_utils,                      ONLY: zhwwf
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mm_mddiag
  PUBLIC :: give_scr_mm_mddiag

CONTAINS

  ! ==================================================================
  SUBROUTINE mm_mddiag(c0,cm,c1,c2,sc0,vpp,gamx,gamy)
    ! ==--------------------------------------------------------------==
    ! Qmmm
    ! Qmmm locals
    COMPLEX(real_8)                          :: c0(:,:,:), cm(:), c1(*), &
                                                c2(:,:,:), sc0(:)
    REAL(real_8)                             :: vpp(:), gamx(:), gamy(:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'mm_mddiag'

    CHARACTER(len=10)                        :: prch
    CHARACTER(len=100)                       :: filen
    CHARACTER(len=30)                        :: tag, tagw
    COMPLEX(real_8)                          :: soc_el(3)
    COMPLEX(real_8), ALLOCATABLE             :: psi(:,:)
    INTEGER :: i, ia, iat, ierr, ifcalc, il_psi_1d, il_psi_2d, il_rhoe_1d, &
      il_rhoe_2d, in, irec(100), is, isp, isub, itemp, j, j1, j2, j3, k, &
      lenext, loopnfi, lscr, nfimax, nfimin, nmin, nmm, nnat, nnx, nroot, &
      nstate, nx, out, save_nfi
    LOGICAL                                  :: ferror, fexist, lmetares, &
                                                oldstatus, statusdummy, &
                                                thalfmesh
    REAL(real_8) :: cubecenter(3), disa, dummy, econs, ekin1, ekin2, ekincp, &
      ekinh1, ekinh2, ekinp, enose, enosp, lmio(3), rmem, tcpu, temp1, temp2, &
      tempp, time1, time2, vcmio(4)
    REAL(real_8), ALLOCATABLE :: eigv(:,:), norms(:), rhoe(:,:), rhoes(:), &
      rinp(:), rm1(:), scr(:), soc_array(:), taui(:,:,:), tauio(:,:), &
      taur(:,:,:)

! Variables
! META DYNAMICS
! ..TSH[
! RHOO_1D is temporary alias for RHOO to allocate through MEMORY90
! ..TSH]
! 
! EHR[
! EHRENFEST
! complex(8) :: CTT(NGWK,*)
! POINTER    (IP_CTT,CTT)
! 

    CALL tiset(procedureN,isub)
    ! 
    thalfmesh=.FALSE.
    cubecenter(1) = 0._real_8
    cubecenter(2) = 0._real_8
    cubecenter(3) = 0._real_8
    nnat = 0
    DO isp=1,ions1%nsp
       DO iat=1,ions0%na(isp)
          cubecenter(1) = cubecenter(1) + tau0(1,iat,isp)
          cubecenter(2) = cubecenter(2) + tau0(2,iat,isp)
          cubecenter(3) = cubecenter(3) + tau0(3,iat,isp)
          nnat = nnat + 1
       ENDDO
    ENDDO
    cubecenter(1) =  cubecenter(1) / REAL(nnat,kind=real_8)
    cubecenter(2) =  cubecenter(2) / REAL(nnat,kind=real_8)
    cubecenter(3) =  cubecenter(3) / REAL(nnat,kind=real_8)
    ALLOCATE(rhoes(fpar%nnr1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! EHR]
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
       ! EHR[
       IF (cntl%tmdeh) THEN
          ! CALL MEMORY(IP_CTT,2*NGWK*N,'CTT')
          ALLOCATE(norms(nkpt%ngwk*crge%n),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ! CALL AZZERO(CTT,2*NGWK*N)
       ELSE
          ! CALL MEMORY(IP_CTT,1,'CTT')
       ENDIF
       ! EHR]
       IF (tshl%isc) THEN
          do_soc_because_sh=.FALSE.
          CALL isc_read_start
          IF (td01%ns_sin.NE.0) THEN
             md_on_sin=.TRUE.
             md_on_tri=.FALSE.
             lrsym=1
          ELSE
             md_on_sin=.FALSE.
             md_on_tri=.TRUE.
             lrsym=3
          ENDIF
          nroot=MAX(td01%ns_sin,td01%ns_tri)
          nstate=crge%n

          ALLOCATE(cs1(ncpw%ngw,nstate,nroot+1),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(ct1(ncpw%ngw,nstate,nroot+1),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(tausoc(3,maxsys%nax,maxsys%nsx),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(ene_tri(10,nroot+1),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(ene_sin(10,nroot+1),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(soc_array(nroot+1),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          CALL zeroing(ene_sin)!,size(ene_sin))
          CALL zeroing(ene_tri)!,size(ene_tri))
          fexist=.FALSE.
          IF (paral%io_parent) INQUIRE(file='ISC_SOC.dat',exist=fexist)
          IF (paral%io_parent .AND. fexist) THEN
             do_soc_because_sh=.FALSE.
             !call fileopen(50,'ISC_SOC.dat',fo_app,ferror) 
             OPEN(unit=50,file='ISC_SOC.dat',status='OLD',position='APPEND')
             BACKSPACE(50) 
             READ(50,*) j1,j2,j3,(soc_array(i),i=1,nroot) 
             DO i=1,nroot 
                soc_array(i)=soc_array(i)/facthcm 
             ENDDO
             CALL fileclose(50) 
          ELSE
             do_soc_because_sh=.TRUE. ! do soc calculation at the first step 
          ENDIF
          CALL mp_bcast(fexist,parai%io_source,parai%cp_grp)
          CALL mp_bcast(do_soc_because_sh,parai%io_source,parai%cp_grp)
          IF (fexist) CALL mp_bcast(soc_array,SIZE(soc_array),parai%io_source,parai%cp_grp)
          fexist=.FALSE.
          IF (paral%io_parent) INQUIRE(file="ISC_ENERGIES_ex.dat",exist=fexist)
          IF (fexist.AND.paral%io_parent) THEN
             OPEN(unit=50,file="ISC_ENERGIES_ex.dat",status='unknown',position='append')
             BACKSPACE(50)
             READ(50,*) k, &
                  (ene_sin(1,j),j=2,nroot+1),&
                  (ene_tri(1,j),j=2,nroot+1)
             CLOSE(50)
          ENDIF
          CALL mp_bcast(fexist,parai%io_source,parai%cp_grp)
          IF (fexist) CALL mp_bcast(ene_sin,SIZE(ene_sin),parai%io_source,parai%cp_grp)
          IF (fexist) CALL mp_bcast(ene_tri,SIZE(ene_tri),parai%io_source,parai%cp_grp)
       ENDIF

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
       nacc = 22
       iteropt%nfi  = 0
       ropt_mod%modens=.FALSE.
       ropt_mod%engpri=.FALSE.
       ropt_mod%calste=cntl%tpres

       IF (cntl%tqmmm.AND.pslo_com%tivan) THEN
          CALL stopgm('MM_MDDIAG','Qmmm WITH VANDERBILT UNSUPPORTED',& 
               __LINE__,__FILE__)
       ENDIF

       ALLOCATE(eigv(crge%n,nkpt%nkpts*clsd%nlsd),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ! ==------------------------------------------------------------==
       ! Extrapolation
       IF (cntl%textrap) THEN
          lenext=2*nkpt%ngwk*nstate*nkpt%nkpts*cnti%mextra
          rmem = 16._real_8*lenext*1.e-6_real_8
          ALLOCATE(cold(nkpt%ngwk,crge%n,nkpt%nkpnt,lenext/(crge%n*nkpt%ngwk*nkpt%nkpnt)),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          CALL zeroing(cold)
          IF (paral%parent) WRITE(6,'(A,T51,F8.3,A)') ' MM_MDDIAG| '//&
               ' EXTRAPOLATION WAVEFUNCTION HISTORY TAKES ',rmem,&
               ' MBYTES'
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
          IF (td01%ns_tri.GT.0.OR.tshl%isc) nmm=2
       ENDIF
       il_psi_1d=nmm*il_psi_1d
       ALLOCATE(psi(il_psi_1d,il_psi_2d),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       IF (lqmmm%qmmm) THEN
          CALL give_scr_mm_mddiag(lscr,tag)
       ELSE
          CALL give_scr_mddiag(lscr,tag)
       ENDIF
       ALLOCATE(scr(lscr),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ELSE
       nmin=10
       ALLOCATE(eigv(crge%n,nkpt%nkpts*clsd%nlsd),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL give_scr_mddiag(lscr,tag)
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
    ALLOCATE(tauio(3,ions1%nat),STAT=ierr)
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

    CALL zeroing(fion)!,3*maxsys%nax*maxsys%nsx)
    CALL zeroing(taui)!,3*maxsys%nax*maxsys%nsx)
    CALL zeroing(taup)!,3*maxsys%nax*maxsys%nsx)
    ! ..TSH[
    IF (tshl%tdtully) THEN
       tshi%nshs=nstate+tshi%nroot_sh
       ALLOCATE(shsigma(2,(tshi%nroot_sh+1),(tshi%nroot_sh+1)),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(shsigmaf(2,(tshi%nroot_sh+1),(tshi%nroot_sh+1)),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(v_sh(2,(tshi%nroot_sh+1)),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(c_sh(2,(tshi%nroot_sh+1)),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! ..TSH]
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
    IF (paral%io_parent.AND.cgrest_i%n_cg.GE.1)&
         CALL fileopen(40,'restrain.out',fo_app,ferror)
    IF (paral%io_parent.AND.cntl%tqmmm.AND.cntl%tddft)&
         CALL fileopen(41,'ENERGIES_EX',fo_app,ferror)
    ! ==--------------------------------------------------------------==
    CALL mm_dim(mm_go_mm,statusdummy)
    IF (paral%qmnode ) THEN
       IF (cprint%iprint_step.EQ.0) cprint%iprint_step=cnti%nomore+1
       ! ==--------------------------------------------------------------==
       ! Set IREC to restart file.
       CALL read_irec(irec)
       ! INITIALIZATION OF WAVEFUNCTION AND COORDINATES
       ! EHR[
       IF (cntl%tmdeh.AND.(td_prop%stextpot.OR.td_prop%td_extpot.OR.td_prop%tpointch))&
            CALL stopgm("MM_MDDIAG","EXTERNAL FIELD AND Qmmm NOT POSSIBLE",& 
            __LINE__,__FILE__)
       ! EHR]
       CALL initrun(irec,c0,c2,sc0,rhoe,psi,eigv)
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
       CALL zeroing(cm(1:ncpw%ngw*nstate))!,ngw*nstate)
       IF (.NOT.lqmmm%qmmm_reflex) THEN
          CALL mm_translate_qmmm(tau0,c0,cm,nstate)
       ENDIF
       ! 28/3/03
       ! EHR[
       IF (cntl%tmdeh) THEN
          IF (cntl%start_real) THEN
             ! LOAD WAVEFUNCTIONS FORM FILE DENSITY.-NM
             IF (td_prop%read_wf.EQ.1) THEN
                CALL load_ex_states(c0)
                CALL getnorm_k(c0,crge%n,norms)
                DO i=1,crge%n
                   CALL dscal(2*nkpt%ngwk,1._real_8/SQRT(norms(i)),c0(1,i,1),1)
                ENDDO
             ELSEIF (td_prop%read_wf.EQ.2) THEN
                ! keep the orbitals read in the RESTART
             ELSE
                ! READ FROM WAVEFUNCTIONS
                CALL tmprd_prop(c0,crge%n,itemp,'wavefunctions')
             ENDIF
          ELSE
             ! READ FROM FORT.120
             CALL tmprd_prop(c0,crge%n,itemp,'wavefunctions')
          ENDIF
          CALL rhoofr_c(c0,rhoe,psi(:,1),crge%n)
          CALL getnorm_k(c0,crge%n,norms)
          IF (paral%parent) THEN
             IF (paral%io_parent)&
                  WRITE(6,*) " NORM OF THE STARTING KS STATES"
             IF (paral%io_parent)&
                  WRITE(6,'(2X,63("-"))')
             IF (paral%io_parent)&
                  WRITE(6,'(2X,7(1X,F8.3))') (norms(i),i=1,crge%n)
             IF (paral%io_parent)&
                  WRITE(6,'(2X,63("-"))')
          ENDIF
       ENDIF
       ! EHR]
       ! ..TSH[
       IF (tshl%tdtully) THEN
          shfilen1="SH_WAVEFUNCTIONS"
          shfilen2="SH_LRWAVEFUNCTIONS"
          IF (paral%parent) THEN
             IF (paral%io_parent)&
                  INQUIRE(file=shfilen1,exist=fexist)
          ENDIF
          CALL mp_bcast(fexist,parai%source,parai%allgrp)
          IF (.NOT.fexist) THEN
             tshi%shstep=0
             ! CALL DCOPY(2*NGW*NSHS,C0(1,1,1),1,CK(1,1),1)
          ELSE
             CALL tmprd(ck,tshi%nshs,1,tshi%shstep,shfilen1)
             CALL mp_bcast(tshi%shstep,parai%source,parai%allgrp)
          ENDIF
       ENDIF
       ! ..TSH]

       ! INITIALIZE VELOCITIES
       CALL mm_dim(mm_go_mm,statusdummy)
       IF (paral%parent) CALL detdof(tau0,taur)

       ! INITIALIZE METADYNAMICS VARIABLES used also for 
       ! SHOOTING from SADDLE POINT with RANDOM VELOCITIES
       IF (paral%parent .AND. (lmeta%lcolvardyn .OR. lmeta%lsadpnt)) THEN
          CALL colvar_structure(tau0,taur)
       ENDIF

       IF ((irec(irec_vel).EQ.0).AND.&
            (.NOT.restart1%rgeo).AND.(.NOT.rtr_l%restart_traj)) THEN
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
       IF (clc%classical) CALL zeroing(cm(1:ncpw%ngw*nstate))!,ngw*nstate)
       ! COMPUTE THE IONIC TEMPERATURE TEMPP
       IF (paral%parent) THEN
          CALL ekinpp(ekinp,velp)
          tempp=ekinp*factem*2._real_8/glib
          ! ..TSH[
          tempp_sh=tempp
          ! ..TSH]
       ENDIF
       IF (cntl%trevers) THEN
          ! invert velocities (useful for path sampling)
          CALL dscal(3*maxsys%nax*maxsys%nsx,-1._real_8,velp,1)
       ENDIF
       ! RESET ACCUMULATORS
       IF (paral%parent.AND.irec(irec_ac).EQ.0)&
            CALL resetac(tau0,taui,iteropt%nfi)
       ! 
       CALL write_irec(irec)
       ! INITIALIZE FORCES
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
       IF (vdwl%vdwd.AND..NOT.vdwwfl%trwannc) THEN
          CALL localize2(tau0,c0,c2,sc0,nstate)
       ENDIF
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
    ! EHR[
    IF (cntl%tmdeh) THEN
       CALL mm_forces_prop(crge%n,c0,c1,c2(:,:,1),cm,sc0,cm(nx:),vpp,eigv,&
            rhoe,psi,&
            tau0,velp,taui,fion,ifcalc,&
            irec,.TRUE.,.TRUE.)
    ELSE
       CALL mm_forces_diag(crge%n,c0(:,:,1),c1,c2,cm,sc0,cm(nx:),vpp,eigv,&
            rhoe,psi, tau0,velp,taui,fion,ifcalc,&
            irec,.TRUE.,.TRUE.)
    ENDIF
    ! EHR]
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
          filen='ENERGIES'
          IF (paral%io_parent)&
               CALL fileopen(3,filen,fo_app+fo_verb,ferror)
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
       ! ==--------------------------------------------------------------==
       ! == END INITIALIZATION                                           ==
       ! ==--------------------------------------------------------------==
       CALL mm_dim(mm_go_mm,statusdummy)
       IF (teststore(0).AND.cntl%tsampl)&
            CALL zhwwf(2,irec,c0,c2,crge%n,eigv,tau0,velp,taui,iteropt%nfi)
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
          ! CALL WRPRINT_MD(EIGV,F,AMU,N,TAU0,FION,
          ! &                  0._real_8,TEMPP,ETOT,ECONS,0._real_8,DISA,
          ! &                  TCPU,.FALSE.,NFI,0)
       ENDIF
    ENDIF ! (qmnode)
    ! ==================================================================
    ! ==          THE BASIC LOOP FOR MOLECULAR DYNAMICS               ==
    ! ==                 USING VELOCITY VERLET                        ==
    ! ==================================================================
    infi=0
    nfimin=iteropt%nfi+1
    nfimax=iteropt%nfi+cnti%nomore
    DO loopnfi=nfimin,nfimax
       CALL mm_dim(mm_go_mm,statusdummy)
       CALL mp_sync(parai%qmmmgrp)
       ! call mp_sync(ALLGRP)
       infi=infi+1
       ! ..TSH[
       IF (tshl%tdtully) tshi%shstep=tshi%shstep+1
       ! ..TSH]
       iteropt%nfi=loopnfi
       IF (paral%qmnode) THEN
          time1=m_walltime()
          comvl%subcom=comvl%tsubcom.AND.MOD(iteropt%nfi-1,comvl%ncomv).EQ.0
          comvl%subrot=comvl%tsubrot.AND.MOD(iteropt%nfi-1,comvl%nrotv).EQ.0
          IF (hubbu%pfrqom.gt.0) THEN
             hubbu%tpom=MOD(iteropt%nfi-1,hubbu%pfrqom).EQ.0
          ELSE
             hubbu%tpom=.False.
          ENDIF
          ! ANNEALING
          CALL anneal(velp,c2,crge%n,scr)
          CALL berendsen(velp,c2,nstate,scr,0.0_real_8,0.0_real_8)
          ! SUBTRACT CENTER OF MASS VELOCITY
          IF (paral%parent.AND.comvl%subcom) CALL comvel(velp,vcmio,.TRUE.)
          ! SUBTRACT ROTATION AROUND CENTER OF MASS
          IF (paral%parent.AND.comvl%subrot) CALL rotvel(tau0,velp,lmio,tauio,.TRUE.)
          ! UPDATE NOSE THERMOSTATS
          CALL noseup(velp,c2,crge%n,1)
          ! UPDATE VELOCITIES
          IF (paral%parent) CALL velupi(velp,fion,1)
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
       ENDIF! (qmnode)
       ! 
       CALL mp_bcast(taup,SIZE(taup),parai%io_source,parai%qmmmgrp)
       CALL zeroing(cm(1:nkpt%ngwk*nstate))!,nkpt%ngwk*nstate)

       IF (.NOT.lqmmm%qmmm_reflex) THEN
          CALL mm_translate_qmmm(taup,c0,cm,nstate)
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
             CALL extrapwf(infi,c0,scr,cold,nnow,numcold,crge%n,cnti%mextra)
          ENDIF
          ! mb - Wannier stuff for vdW-WC
          vdwwfl%twannup=vdwwfl%twannup.OR.(infi.EQ.1.AND..NOT.vdwwfl%trwannc)
          IF (vdwl%vdwd.AND.vdwwfl%twannup) THEN
             CALL localize2(taup,c0,c2,sc0,nstate)
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
          ! save the cntl%md step number
          save_nfi=iteropt%nfi
          ! localisation at first iteration + PT
          CALL mddiag_interaction_p(crge%n,c0,c2,cm,sc0,&
               cm(nx:),vpp,eigv,rhoe,psi, &!SCR,LSCR,& !vw doesnt match procedure args
               taup,velp,taui,fion,ifcalc,irec,.TRUE.,.FALSE.)
          ! sets the number of iteration and recover cntl%md step number
          ifcalc=iteropt%nfi
          iteropt%nfi=save_nfi

          ! switch off the info printing
          dmbi%inter_pt_firstcall=.FALSE.
       ELSE
          ! EHR[
          IF (cntl%tmdeh) THEN
             CALL mm_forces_prop(crge%n,c0,c1,c2(:,:,1),cm,sc0,cm(nx:),vpp,eigv,&
                  rhoe,psi,&
                  taup,velp,taui,fion,ifcalc,&
                  irec,.TRUE.,.FALSE.)
          ELSE
             CALL mm_forces_diag(crge%n,c0(:,:,1),c1,c2,cm,sc0,cm(nx:),vpp,eigv,&
                  rhoe,psi, taup,velp,taui,fion,ifcalc,&
                  irec,.TRUE.,.FALSE.)
          ENDIF
          ! EHR]
       ENDIF
       CALL mm_dim(mm_go_mm,statusdummy)

       ! ..TSH[
       IF (tshl%tdtully.AND.tshl%tully_sh) THEN
          IF (paral%io_parent)&
               WRITE(6,*) ' velocity rescaling after SH'
          IF (paral%parent) CALL rscvp_sh(dt_sh,tempp,velp)
          tshl%tully_sh=.FALSE.
       ENDIF
       ! ..TSH]

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
          IF (paral%parent) CALL geofile(taup,velp,'WRITE')
          ! COMPUTE THE IONIC TEMPERATURE TEMPP
          tempp=0.0_real_8
          IF (paral%parent) THEN
             ! calculate local kinetic temperature in QM and MM subsystem
             ! and write to "QM_TEMP"
             CALL ekinpp(ekinp,velp)
             tempp=ekinp*factem*2._real_8/glib
             IF (wp_i%nfi_lt.GT.0) THEN
                CALL mm_localt(velp,rmass%pma,factem,tempp,glib,iteropt%nfi,wp_i%nfi_lt)
             ENDIF
          ENDIF
          ! ..TSH[
          tempp_sh=tempp
          ! ..TSH]
          ! IONIC TEMPERATURE CONTROL
          IF (paral%parent) CALL rscvp(temp1,temp2,tempp,velp)
          ! SUBTRACT ROTATION AROUND CENTER OF MASS
          IF (paral%parent.AND.comvl%subrot) CALL rotvel(tau0,velp,lmio,tauio,.FALSE.&
               )
          ! SUBTRACT CENTER OF MASS VELOCITY
          IF (paral%parent.AND.comvl%subcom) CALL comvel(velp,vcmio,.FALSE.)
          ! UPDATE NOSE THERMOSTATS
          CALL noseup(velp,c2,crge%n,1)
          CALL berendsen(velp,c2,nstate,scr,0.0_real_8,0.0_real_8)
          ! ANNEALING
          CALL anneal(velp,c2,crge%n,scr)
          IF (paral%parent) THEN
             CALL ekinpp(ekinp,velp)
             tempp=ekinp*factem*2._real_8/glib
          ENDIF
          ! MEAN SQUARE DISPLACEMENT OF DIFFERENT IONIC SPECIES
          IF (paral%parent) CALL dispp(taup,taui,disa)
          ! ENERGY OF THE NOSE THERMOSTATS
          IF (paral%parent) CALL noseng(iteropt%nfi,velp,enose,enosp,dummy,1)
          ! CALCULATE PROPERTIES DURING SIMULATION.
          cntl%caldip=cntl%tdipd.AND.MOD(iteropt%nfi-1,cnti%npdip).EQ.0
          CALL mm_dim(mm_go_qm,statusdummy)
          CALL propcal(c0,c2(:,:,1),cm,sc0,taup,eigv,crge%f,ener_com%amu,&
               rhoe,psi,crge%n,nkpt%nkpnt,iteropt%nfi,infi)
          CALL mm_dim(mm_go_mm,statusdummy)
          ! PRINTOUT the evolution of the accumulators every time step
          IF (paral%parent) THEN
             IF (paral%io_parent)&
                  CALL fileclose(3)
             IF (paral%io_parent)&
                  CALL fileopen(3,filen,fo_app,ferror)
             econs=ekinp+ener_com%etot+enose+enosp+ener_com%ecnstr+ener_com%erestr
             time2=m_walltime()
             tcpu=(time2-time1)*0.001_real_8
             CALL wrprint_md(eigv,crge%f,ener_com%amu,crge%n,taup,fion,&
                  0._real_8,tempp,ener_com%etot,econs,0._real_8,disa,&
                  tcpu,.FALSE.,iteropt%nfi,infi)
             !IF (tshl%isc) etot_isc_md=ener_com%etot
             ! UPDATE ACCUMULATORS
             CALL paccc(tempp,ener_com%etot,econs,enose,enosp,ener_com%ecnstr,ener_com%erestr,&
                  ener_com%ebogo,disa,tcpu,iteropt%nfi,1)
             ! STORE IONIC COORDINATES AND VELOCITIES FOR STATISTICS
             ropt_mod%movie=rout1%mout .AND. MOD(iteropt%nfi-1,cnti%imovie).EQ.0
             ropt_mod%rprint=rout1%rout .AND. MOD(iteropt%nfi-1,cnti%ntraj).EQ.0
             ropt_mod%txyz=rout1%xtout .AND. MOD(iteropt%nfi-1,cnti%ntraj).EQ.0
             ropt_mod%tdcd=rout1%dcout .AND. MOD(iteropt%nfi-1,cnti%ntraj).EQ.0
             ! EHR[
             IF (cntl%tmdeh) THEN
                IF (paral%io_parent)&
                     WRITE(200,'(9(F10.6,1X),I6)') tempp,ener_com%etot,econs,enose,&
                     enosp,ener_com%ecnstr,ener_com%ebogo,disa,tcpu,iteropt%nfi
                DO i=1,ions1%nsp
                   DO j=1,ions0%na(i)
                      IF (paral%io_parent)&
                           WRITE(201,'(I7,6(2X,F22.14))')&
                           iteropt%nfi,(taup(k,j,i),k=1,3),(velp(k,j,i),k=1,3)
                   ENDDO
                ENDDO
             ENDIF
             ! EHR]
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
          IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)
          IF (iteropt%nfi.EQ.nfimax) THEN
             soft_com%exsoft=.TRUE.
             soft_com%exnomore=.TRUE.
          ENDIF
          CALL mm_dim(mm_go_mm,statusdummy)
       ENDIF ! qmnode

       ! STOP THE RUN IF THE USER HAS SET THE SIGNAL 30
       in=0
       out=0
       IF (soft_com%exsoft) in=1
       CALL mp_sum(in,out,parai%qmmmgrp)
       IF (out.NE.0) soft_com%exsoft=.TRUE.

       IF (paral%qmnode) THEN
          ! periodic output of density/wavefunction etc.
          IF (rout1%rhoout.AND.(rout1%nrhoout.GT.0)) THEN
             IF (MOD(iteropt%nfi-1,rout1%nrhoout).EQ.0) THEN
                CALL rhopri(c0,tau0,rhoe,psi(:,1),nstate,nkpt%nkpnt)
             ENDIF
          ENDIF
          IF (teststore(iteropt%nfi).OR.soft_com%exsoft)&
               CALL zhwwf(2,irec,c0,c2,crge%n,eigv,taup,velp,taui,iteropt%nfi)

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
          ! EHR[
          IF (cntl%tmdeh.AND.store1%srho.AND.&
               ((infi.EQ.1).OR.(MOD(infi,store1%istore).EQ.0))) THEN
             CALL mm_dim(mm_go_qm,statusdummy)
             CALL rhopri(c0,tau0,rhoe,psi(:,1),crge%n,1)
             tagw="wavefunctions"
             CALL tmpwr_prop(c0,crge%n,infi,tagw)
             CALL mm_dim(mm_go_mm,statusdummy)
          ENDIF
          ! EHR]
       ENDIF! (qmnode)
       ! ..TSH[
       IF (tshl%tdtully) THEN
          CALL tmpwr(ck,tshi%nshs,1,tshi%shstep,shfilen1)
          CALL dcopy(2*ncpw%ngw*tshi%nshs,c0(1,1,1),1,ck(1,1),1)
       ENDIF
       ! ..TSH]
       IF (soft_com%exsoft) GOTO 100
       ! ISC[
       IF (tshl%isc) THEN
          CALL dcopy(3*maxsys%nax*maxsys%nsx,tau0,1,tausoc,1)
          IF (MOD(infi,nskip).EQ.0.OR.do_soc_because_sh) THEN
             !if (paral%qmnode)  write(*,*) 'etot',ener_com%etot
             CALL do_tsh_isc(c0(:,:,1),c1,c2(:,:,1),sc0,rhoe,psi,eigv,ener_com%etot,&
                  soc_el,soc_array,infi,nstate,nroot)
             CALL do_tsh_isc_lz(soc_array,cntr%delt_elec,nroot,nskip)
          ENDIF
          IF (do_soc_because_sh) do_soc_because_sh=.FALSE.
          CALL isc_write_start
       ENDIF
       ! ISC]
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
       CALL mm_dim(mm_go_qm,statusdummy)
       IF (rout1%rhoout.AND.(rout1%nrhoout.LE.0))&
            CALL rhopri(c0,tau0,rhoe,psi(:,1),crge%n,nkpt%nkpnt)
       CALL mm_dim(mm_go_mm,statusdummy)
       ! Print accumulators.
       IF (paral%parent) CALL paccc(tempp,ener_com%etot,econs,enose,enosp,ener_com%ecnstr,&
            ener_com%erestr,ener_com%ebogo,disa,tcpu,iteropt%nfi,0)
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
       IF ((paral%parent).AND.paral%io_parent)CALL fileclose(3)
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
       IF (comvl%tsubrot) THEN
          DEALLOCATE(tauio,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
       ENDIF
       ! ..TSH[
       IF (tshl%tdtully) THEN
          DEALLOCATE(shsigma,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
          DEALLOCATE(shsigmaf,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
          DEALLOCATE(v_sh,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
          DEALLOCATE(c_sh,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
       ENDIF
       ! ..TSH]
       ! ==--------------------------------------------------------------==
    ELSE
       DEALLOCATE(eigv,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF ! (qmnode)
    ! EHR[
    DEALLOCATE(rhoes,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! CALL FREEM(IP_CTT)
    IF (cntl%tmdeh) DEALLOCATE(norms,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! EHR]
    ! SOC[
    IF (tshl%isc) THEN
       DEALLOCATE(cs1,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(ct1,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(tausoc,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(ene_tri,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(ene_sin,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(soc_array,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! SOC]
    ! ==--------------------------------------------------------------==
    CALL mm_dim(mm_revert,oldstatus)
    IF ((paral%parent.AND.cgrest_i%n_cg.GT.1).AND.paral%io_parent)&
         CALL fileclose(40)
    IF ((paral%parent.AND.cntl%tqmmm.AND.cntl%tddft).AND.paral%io_parent)&
         CALL fileclose(41)
    CALL mp_sync(parai%qmmmgrp)
#endif
    CALL tihalt(procedureN,isub)
    RETURN
  END SUBROUTINE mm_mddiag
  ! ==================================================================
  SUBROUTINE give_scr_mm_mddiag(lmddiag,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lmddiag
    CHARACTER(len=30)                        :: tag

    INTEGER :: lcalc_alm, lcopot, lforces_diag, linitrun, lmoverho, lmtd, &
      lpropcal, lrhoofr, lrhopri, lrinitwf, lrnlsm, ltddft, nstate

! Variables
! real(8) :: ALM(*),AFNL(*),BILN(*)
! ==--------------------------------------------------------------==

    nstate=crge%n
    lrnlsm=0
    lcalc_alm=0
    lcopot=0
    lrhopri=0
    lmoverho=0
    ltddft=0
    linitrun=0
    CALL give_scr_initrun(linitrun,tag)
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
    CALL give_scr_meta_extlagr(lmtd,tag)
    lmddiag=MAX(lrinitwf,lrnlsm,lrhoofr,lforces_diag,ltddft,linitrun,&
         lcopot,lcalc_alm,lrhopri,lpropcal,lmoverho,lmtd)
    IF (cntl%tqmmm) lmddiag=MAX(lmddiag,fpar%kr1*fpar%kr2s*fpar%kr3s)
    IF (cntl%tqmmm) lmddiag=MAX(lmddiag,maxsys%nax*maxsys%nsx*3)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_mm_mddiag
  ! ==================================================================

END MODULE mm_mddiag_utils
