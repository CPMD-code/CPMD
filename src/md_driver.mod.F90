#if defined(__USE_IBM_HPM)
#define USE_IBM_HPM
#endif

MODULE md_driver
  USE andp,                            ONLY: rin0,&
                                             rmix,&
                                             rout0
  USE andr,                            ONLY: andr2
  USE anneal_utils,                    ONLY: anneal,&
                                             berendsen,&
                                             dampdyn,&
                                             tempramp
  USE atwf,                            ONLY: tmovr
  USE box_boundary_utils,              ONLY: box_boundary
  USE bs_forces_diag_utils,            ONLY: bs_forces_diag
  USE bsym,                            ONLY: bsclcs,&
                                             bsfac
  USE calc_alm_utils,                  ONLY: calc_alm,&
                                             give_scr_calc_alm
  USE cdft_utils,                      ONLY: cdft_reinit,&
                                             cdft_w,&
                                             init_cdft,&
                                             write_w
  USE cdftmod,                         ONLY: cdftcom,&
                                             wdiff
  USE cnst,                            ONLY: factem
  USE cnst_dyn,                        ONLY: ekincv,&
                                             fhills,&
                                             lmeta,&
                                             ltcglobal,&
                                             ncolvar,&
                                             vharm
  USE comvel_utils,                    ONLY: comvel
  USE comvelmod,                       ONLY: comvl
  USE coor,                            ONLY: fion,&
                                             tau0,&
                                             taup,&
                                             velp
  USE copot_utils,                     ONLY: copot,&
                                             give_scr_copot
  USE cotr,                            ONLY: cotc0,&
                                             cnsval,&
                                             ntcnst
  USE cppt,                            ONLY: inyh
  USE ddipo_utils,                     ONLY: give_scr_ddipo
  USE detdof_utils,                    ONLY: detdof
  USE dispp_utils,                     ONLY: dispp
  USE dynit_utils,                     ONLY: dynit
  USE efld,                            ONLY: extf
  USE ekinpp_utils,                    ONLY: ekinpp
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE extrap_utils,                    ONLY: extrap
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_app,&
                                             fo_nofp,&
                                             fo_verb
  USE finalp_utils,                    ONLY: finalp
  USE fint,                            ONLY: fint1
  USE forcep_utils,                    ONLY: rhoe_psi_size
  USE forces_diag_utils,               ONLY: forces_diag,&
                                             give_scr_forces_diag
  USE forces_prop_utils,               ONLY: forces_prop,&
                                             give_scr_forces_prop
  USE geofile_utils,                   ONLY: geofile
  USE gle_utils,                       ONLY: gle_init,&
                                             gle_step
  USE glemod,                          ONLY: glepar
  USE gsize_utils,                     ONLY: gsize
  USE hfxmod,                          ONLY: hfxc3
  USE hubbardu,                        ONLY: hubbu
  USE initrun_driver,                  ONLY: initrun
  USE initrun_utils,                   ONLY: give_scr_initrun
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE isos,                            ONLY: isos1
  USE kinds,                           ONLY: int_1,&
                                             int_2,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8
  USE kpts,                            ONLY: tkpts
  USE linres,                          ONLY: &
       c_sh, cdotlct, ck, dt_sh, lrsym, shfilen1, shfilen2, shfilen3, shlct, &
       shsigma, shsigmaf, td01, tempp_sh, tshf, tshi, tshl, v_sh, xfmqc
  USE localize_utils,                  ONLY: localize2
  USE lr_tddft_utils,                  ONLY: give_scr_lr_tddft,&
                                             lr_tddft
  USE machine,                         ONLY: m_walltime
  USE mddiag_interaction_p_utils,      ONLY: mddiag_interaction_p
  USE meta_colvar_inp_utils,           ONLY: colvar_structure
  USE meta_colvar_utils,               ONLY: meta_colvar,&
                                             meta_colvar_mw
  USE meta_exl_mult_utils,             ONLY: meta_ext_mul
  USE meta_exlagr_methods,             ONLY: give_scr_meta_extlagr,&
                                             meta_extlagr,&
                                             meta_extlagr_mw
  USE meta_exlagr_utils,               ONLY: ekincv_global
  USE meta_multiple_walkers_utils,     ONLY: mw_assign_filenames,&
                                             mw_filename
  USE mimic_wrapper,                   ONLY: mimic_ifc_collect_energies,&
                                             mimic_ifc_collect_forces,&
                                             mimic_ifc_sort_fragments,&
                                             mimic_ifc_init,&
                                             mimic_revert_dim,&
                                             mimic_save_dim,&
                                             mimic_ifc_send_coordinates,&
                                             mimic_sum_forces,&
                                             mimic_switch_dim,&
                                             mimic_update_coords,&
                                             mimic_control,&
                                             mimic_energy,&
                                             mimic_subsystem_temperatures,&
                                             mimic_subsystem_dof
  USE mfep,                            ONLY: mfepi
  USE mm_extrap,                       ONLY: cold , cold_high, &
                                             numcold, numcold_high, &
                                             nnow, nnow_high
  USE moverho_utils,                   ONLY: give_scr_moverho,&
                                             moverho
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum,&
                                             mp_sync
  USE interface_utils,                 ONLY: get_external_forces
  USE mts_utils,                       ONLY: mts, set_mts_functional
  USE mw,                              ONLY: mwi,&
                                             tmw
  USE nlcc,                            ONLY: corel
  USE norm,                            ONLY: gnmax,&
                                             gnorm
  USE nose,                            ONLY: glib
  USE noseng_utils,                    ONLY: noseng
  USE nosepa_utils,                    ONLY: nosepa
  USE noseup_utils,                    ONLY: noseup
  USE nospinit_utils,                  ONLY: nospinit
  USE ovlap_utils,                     ONLY: ovlap
  USE parac,                           ONLY: parai,&
                                             paral
  USE phfac_utils,                     ONLY: phfac
  USE pimd,                            ONLY: grandparent,&
                                             ipcurr,&
                                             np_low,&
                                             supergroup
  USE poin,                            ONLY: potr,&
                                             rhoo
  USE posupi_utils,                    ONLY: posupi
  USE printave_utils,                  ONLY: paccc
  USE printp_utils,                    ONLY: printp,&
                                             printp2,&
                                             print_mts_forces
  USE proppt_utils,                    ONLY: give_scr_propcal,&
                                             propcal
  USE pslo,                            ONLY: pslo_com
  USE puttau_utils,                    ONLY: taucl
  USE rattle_utils,                    ONLY: rattle
  USE readsr_utils,                    ONLY: xstring
  USE resetac_utils,                   ONLY: resetac
  USE response_pmod,                   ONLY: dmbi
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
  USE rotate_utils,                    ONLY: rotate
  USE rotvel_utils,                    ONLY: rotvel
  USE rscvp_utils,                     ONLY: rscvp
  USE sample_utils,                    ONLY: sample_go,&
                                             sample_wait
  USE setbsstate_utils,                ONLY: setbsstate
  USE setirec_utils,                   ONLY: write_irec
  USE sh_tddft_utils,                  ONLY: rscvp_sh,&
                                             sh_extpot,&
                                             sh_lct,&
                                             tmprd,&
                                             tmpwr,&
                                             xf_fill_all_buffers,&
                                             xf_force_components
  USE shake_utils,                     ONLY: cpmdshake,&
                                             init_constraints,&
                                             do_shake,&
                                             do_rattle
  USE soc,                             ONLY: do_tsh_isc,&
                                             do_tsh_isc_lz,&
                                             isc_read_start,&
                                             isc_write_start
  USE soc_types,                       ONLY: &
       cs1, ct1, do_soc_because_sh, ene_sin, ene_tri, facthcm, md_on_sin, &
       md_on_tri, nskip, tausoc
  USE soft,                            ONLY: soft_com
  USE spin,                            ONLY: clsd,&
                                             spin_mod
  USE store_types,                     ONLY: &
       cprint, irec_ac, irec_nop1, irec_nop2, irec_nop3, irec_nop4, irec_vel, &
       restart1, rout1, store1
  USE system,                          ONLY: &
       cnti, cntl, cntr, fpar, maxsys, nacc, ncpw, nkpt, restf, iatpt
  USE td_utils,                        ONLY: initialize_ehrenfest_dyn
  USE testex_utils,                    ONLY: testex,&
                                             testex_mw
  USE teststore_utils,                 ONLY: teststore
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE totstr_utils,                    ONLY: totstr
  USE vdwcmod,                         ONLY: trwanncx,&
                                             twannupx,&
                                             vdwl,&
                                             vdwwfl
  USE velupi_utils,                    ONLY: velupi
  USE wrener_utils,                    ONLY: wrprint_md
  USE wrgeo_utils,                     ONLY: wrgeof
  USE wv30_utils,                      ONLY: zhwwf
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mddiag
  PUBLIC :: give_scr_mddiag

CONTAINS

  ! ==================================================================
  SUBROUTINE mddiag(c0,cm,c1,c2,sc0,vpp,gamx,gamy)
    ! ==--------------------------------------------------------------==
    ! CB: BS case needs two wf (BSFAC=2)
    COMPLEX(real_8), INTENT(inout)           :: c0(:,:,:), cm(:), &
                                                c1(nkpt%ngwk,crge%n,*), &
                                                c2(:,:,:), sc0(:)
    REAL(real_8), INTENT(inout)              :: vpp(:), gamx(:), gamy(:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'mddiag'

    CHARACTER(len=100)                       :: filen, filew, ftmp
    CHARACTER(len=30)                        :: tag
    COMPLEX(real_8)                          :: soc_el(3)
    COMPLEX(real_8), ALLOCATABLE             :: ch(:,:,:), psi(:,:)
    INTEGER :: i, i1, i2, ia, ierr, ifcalc, il_psi_1d, il_psi_2d, il_rhoe_1d, &
      il_rhoe_2d, ipwalk, ipx, irec(100), is, isub, itemp, j, k, lenext, &
      loopnfi, lscr, nfimax, nfimin, nmm, nnx, nroot, nstate, nx, save_nfi, &
      temp
    LOGICAL                                  :: ferror, fexist, ionode, &
                                                lmetares, qmmech, tstrng
    REAL(real_8) :: disa, dummy(2), econs, ek_cv, ekin1, ekin2, ekincp, &
      ekinh1, ekinh2, ekinp, enose, enosp, lmio(3), mm_ekin, mm_temp, rmem, &
      tcpu, temp1, temp2, tempp, time1, time2, vcmio(4)
    REAL(real_8), ALLOCATABLE :: eigs(:), eigv(:,:), norms(:), rhoe(:,:), &
      rinp(:), rm1(:), scr(:), soc_array(:), taui(:,:,:), tauio(:,:), &
      taur(:,:,:)
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

    CALL tiset(procedureN,isub)
    IF (cntl%mimic) THEN
      CALL mimic_save_dim()
      CALL mimic_switch_dim(go_qm=.FALSE.)
    ENDIF
    IF (cntl%tddft.AND.cntl%tresponse) CALL stopgm("MDDIAG",&
         "cntl%tddft.AND.cntl%tresponse NOT POSSIBLE",& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==

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

    ! EHR[
    IF (cntl%tmdeh) THEN
       ALLOCATE(ch(nkpt%ngwk,crge%n,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(extf(fpar%kr1*fpar%kr2s*fpar%kr3s),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(norms(nkpt%ngwk*crge%n),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(ch)!,nkpt%ngwk*n)
       ! CALL SETPOSOP
    ELSE IF (cntl%mimic) THEN
       ALLOCATE(ch(1,1,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                               __LINE__,__FILE__)
       ALLOCATE(extf(fpar%kr1*fpar%kr2s*fpar%kr3s),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                               __LINE__,__FILE__)
    ELSE
       ALLOCATE(ch(1,1,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF

    ! Initialize file names
    ! Walker ID for multiple walker MTD. MTD part is in grandparent 
    ipwalk=1                                                        
    ionode=paral%io_parent
    filen='ENERGIES'
    IF (tmw)THEN                                                    
       ipwalk=mwi%walker_id                                         
       ionode=grandparent                                           
       CALL mw_assign_filenames(ipwalk,filen,ftmp)               
    ELSE IF (cntl%tpmin)THEN                                        
       ! FILEN=FPATH(IAPATH:IEPATH)//'ENERGIES'
       ftmp='ENERGIES_'                                             
       CALL mw_filename(ftmp,filen,ipcurr)                          
       CALL xstring(filen,i1,i2)                                    
       ftmp=filen(i1:i2)//'.'                                       
       CALL mw_filename(ftmp,filen,mfepi%istring)                   
    ENDIF

    ! EHR]
    ! SOC[
    IF (tshl%isc) THEN
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
       IF (paral%io_parent)&
            INQUIRE(file="ISC_SOC.dat",exist=fexist)
       IF (paral%io_parent .AND. fexist) THEN
          do_soc_because_sh=.FALSE.
          OPEN(unit=50,file="ISC_SOC.dat",status='unknown',position='append')
          BACKSPACE(50)
          READ(50,*) j,j,j,(soc_array(i),i=1,nroot)
          DO i=1,nroot
             soc_array(i)=soc_array(i)/facthcm
          ENDDO
          CLOSE(50)
       ELSE
          do_soc_because_sh=.TRUE. ! do soc calculation at the first step
       ENDIF
       CALL mp_bcast(do_soc_because_sh,parai%io_source,parai%cp_grp)
       CALL mp_bcast(fexist,parai%io_source,parai%cp_grp)
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
    ! SOC]
    IF (.NOT.cntl%mimic) THEN
       ALLOCATE(taup(3,maxsys%nax,maxsys%nsx),stat=ierr)
       IF(ierr/=0) CALL stopgm(proceduren,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(fion(3,maxsys%nax,maxsys%nsx),stat=ierr)
       IF(ierr/=0) CALL stopgm(proceduren,'allocation problem',&
            __LINE__,__FILE__)
    END IF
    ! [EXACT FACTORIZATION
    IF (tshl%txfmqc) THEN
       ekincv=0._real_8
       vharm=0._real_8
       ALLOCATE(xfmqc%fion_state(3,maxsys%nax,maxsys%nsx,tshi%nroot_sh+1),&
            stat=ierr)
       IF(ierr/=0) CALL stopgm(proceduren,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(xfmqc%fion_state)
       ALLOCATE(eigs(crge%n),STAT=ierr)
       IF(ierr/=0) CALL stopgm(proceduren,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(eigs)
       ALLOCATE(xfmqc%qmoment(3,maxsys%nax,maxsys%nsx,((tshi%nroot_sh+1)*tshi%nroot_sh)/2),&
            stat=ierr)
       IF(ierr/=0) CALL stopgm(proceduren,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(xfmqc%qmoment)
       ! initialize tauall ! TO BE DONE IN A PROPER WAY WHEN USING MW
       ALLOCATE(xfmqc%tauall(3,maxsys%nax,maxsys%nsx,mwi%nwalk),stat=ierr)
       IF(ierr/=0) CALL stopgm(proceduren,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(xfmqc%tauall)!
       !ALLOCATE(xfmqc%w_ij(xfmqc%n_xftraj,xfmqc%n_xftraj),STAT=ierr)
       !IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
       !     __LINE__,__FILE__)
       ![FEDE
       ALLOCATE(xfmqc%w_ij(3,maxsys%nax,maxsys%nsx,xfmqc%n_xftraj,xfmqc%n_xftraj),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       !FEDE]
       CALL zeroing(xfmqc%w_ij)
       ALLOCATE(xfmqc%k_li((tshi%nroot_sh+1),(tshi%nroot_sh+1)),STAT=ierr) !New - 17.03.16
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(xfmqc%k_li)
       ALLOCATE(xfmqc%f_nli(3,maxsys%nax,maxsys%nsx,tshi%nroot_sh+1),&
            STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(xfmqc%f_nli)
       ALLOCATE(xfmqc%r0_ni(3,maxsys%nax,maxsys%nsx,((tshi%nroot_sh+1)*(tshi%nroot_sh))/2),&
            STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(xfmqc%r0_ni)
       ALLOCATE(xfmqc%d0_ni(3,maxsys%nax,maxsys%nsx,((tshi%nroot_sh+1)*(tshi%nroot_sh))/2),&
            STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(xfmqc%d0_ni)
       ALLOCATE(xfmqc%bf_nli(3,maxsys%nax,maxsys%nsx,tshi%nroot_sh+1,&
            xfmqc%n_xftraj),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(xfmqc%bf_nli)
       ALLOCATE(xfmqc%bfion_state(3,maxsys%nax,maxsys%nsx,tshi%nroot_sh+1,&
            xfmqc%n_xftraj),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(xfmqc%bfion_state)
       ALLOCATE(xfmqc%bfa_ni(3,maxsys%nax,maxsys%nsx,xfmqc%n_xftraj),&
            stat=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(xfmqc%bfa_ni)
       ALLOCATE(xfmqc%fa_ni(3,maxsys%nax,maxsys%nsx),&
            stat=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(xfmqc%fa_ni)
       ! .. MIN[
       ALLOCATE(xfmqc%bw0_nkli(3,maxsys%nax,maxsys%nsx,((tshi%nroot_sh+1)*tshi%nroot_sh)/2,& ! MIN: upper triangle (a_ij, j>i)
            xfmqc%n_xftraj),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(xfmqc%bw0_nkli)
       ALLOCATE(xfmqc%brho_l(tshi%nroot_sh+1,xfmqc%n_xftraj),&
            stat=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(xfmqc%brho_l)
       ! .. MIN]
       ALLOCATE(xfmqc%cf(tshi%nroot_sh+1),&
            stat=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(xfmqc%cf)
       ALLOCATE(xfmqc%bcf(tshi%nroot_sh+1,xfmqc%n_xftraj),&
            stat=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ! initialize sigma ! TO BE FIXED
       ![FEDE
       ALLOCATE(xfmqc%sigma(3,maxsys%nax,maxsys%nsx,xfmqc%n_xftraj),&
            stat=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(xfmqc%sigma)
       !FEDE]
       !xfmqc%sigma=1._real_8
       ! initialize threshold ! TO BE FIXED
       xfmqc%threshold=0.01_real_8
       ! initialize tauall ! TO BE DONE IN A PROPER WAY WHEN USING MW
       ! initialize sigma ! TO BE FIXED
       !DO is = 1,ions1%nsp
       !   DO ia = 1,ions0%na(is)
       !      DO k=1,3
       !         xfmqc%tauall(k,ia,is,i)=tau0(k,ia,is)
       !      ENDDO
       !   ENDDO
       !ENDDO
       !DO i=1,tshi%nroot_sh+1
       !   DO is = 1,ions1%nsp
       !      DO ia = 1,ions0%na(is)
       !         DO k=1,3
       !            xfmqc%f_nli(k,ia,is,i)=0._real_8
       !         ENDDO
       !      ENDDO
       !   ENDDO
       !ENDDO
       !xfmqc%sigma=1._real_8
    ENDIF
    !] EXACT FACTORIZATION
    ALLOCATE(taui(3,maxsys%nax,maxsys%nsx),stat=ierr)
    IF(ierr/=0) CALL stopgm(proceduren,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(taur(3,maxsys%nax,maxsys%nsx),stat=ierr)
    IF(ierr/=0) CALL stopgm(proceduren,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(taup)!,SIZE(taup))
    CALL zeroing(fion)!,SIZE(fion))
    IF (comvl%tsubrot) THEN
       ALLOCATE(tauio(3,ions1%nat),stat=ierr)
       IF(ierr/=0) CALL stopgm(proceduren,'allocation problem',&
            __LINE__,__FILE__)! TODO check dimensions
    ENDIF
    ! Memory for densities
    IF (.NOT.cntl%bsymm) THEN
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
    ENDIF
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
    iteropt%nfi  = 0
    ropt_mod%modens=.FALSE.
    ropt_mod%engpri=.FALSE.
    ropt_mod%calste=cntl%tpres
    ALLOCATE(eigv(crge%n,bsfac*nkpt%nkpts*clsd%nlsd*nstate/crge%n),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ..TSH[
    IF (tshl%tdtully) THEN
       tshi%nshs=nstate+tshi%nroot_sh
       ALLOCATE(shsigma(2,(tshi%nroot_sh+1),(tshi%nroot_sh+1)),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(shsigmaf(2,(tshi%nroot_sh+1),(tshi%nroot_sh+1)),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(shsigma)
       CALL zeroing(shsigmaf)
       ALLOCATE(v_sh(2,(tshi%nroot_sh+1)),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(v_sh)
       ALLOCATE(c_sh(2,(tshi%nroot_sh+1)),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(c_sh)
       ALLOCATE(cdotlct(2,100),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       IF (tshl%txfmqc) THEN
          ALLOCATE(xfmqc%nacv(3,maxsys%nax,maxsys%nsx,100),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          CALL zeroing(xfmqc%nacv)
          ALLOCATE(xfmqc%eigv(tshi%nroot_sh+1),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          CALL zeroing(xfmqc%eigv)
       ENDIF
    ENDIF
    ! ..TSH]
    ! ==--------------------------------------------------------------==
    ! Extrapolation
    IF (cntl%textrap) THEN
       lenext=2*nkpt%ngwk*nstate*nkpt%nkpts*cnti%mextra
       rmem = 16._real_8*lenext*1.e-6_real_8
       ALLOCATE(cold(nkpt%ngwk,crge%n,nkpt%nkpnt,lenext/(crge%n*nkpt%ngwk*nkpt%nkpnt)),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       call zeroing(cold)

       ! allocate array for high level WF extrapolation
       if (mts_pure_dft) then
          allocate(cold_high(nkpt%ngwk,crge%n,nkpt%nkpnt,lenext/(crge%n*nkpt%ngwk*nkpt%nkpnt)),stat=ierr)
          if(ierr/=0) call stopgm(proceduren,'allocation problem: cold_high',&
               __LINE__,__FILE__)
          call zeroing(cold_high)
          numcold_high=0 
          nnow_high=0
          rmem=rmem*2._real_8
       endif

       IF (paral%io_parent)&
            WRITE(6,'(A,T51,F8.3,A)') ' MDDIAG| '&
            // 'EXTRAPOLATION WAVEFUNCTION HISTORY TAKES ',rmem,' MBYTES'
    ENDIF
    ! ==--------------------------------------------------------------==
    ! SCR ALLOCATION AND PARTITION (SCRATCH ARRAY).
    CALL rhoe_psi_size(il_psi_1d=il_psi_1d,il_psi_2d=il_psi_2d,&
         il_rhoe_1d=il_rhoe_1d, il_rhoe_2d=il_rhoe_2d)
    ALLOCATE(rhoe(il_rhoe_1d, il_rhoe_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    nmm=1
    IF (cntl%tddft) THEN
       ALLOCATE(rhoo(il_rhoe_1d, il_rhoe_2d),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(potr(il_rhoe_1d, il_rhoe_2d),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       IF (td01%ns_tri.GT.0.OR.tshl%isc) nmm=2
    ENDIF
    il_psi_1d=il_psi_1d*nmm
    ALLOCATE(psi(il_psi_1d,il_psi_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL give_scr_mddiag(lscr,tag)
    ALLOCATE(scr(lscr),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
99999 IF (cntl%tsampl) THEN
       CALL sample_wait
       IF (cnti%nomore.LT.0) GOTO 10000
    ENDIF
    restf%nfnow=1
    IF (cntl%mimic) THEN
       CALL mimic_ifc_init(rhoe, extf, tau0)
       mimic_control%update_potential = .TRUE.
    END IF
    ! ==--------------------------------------------------------------==
    ! TIME STEP FUNCTIONS
    CALL dynit(ekincp,ekin1,ekin2,temp1,temp2,ekinh1,ekinh2)
    ekinp=0.0_real_8    ! McB: cf. META_EXT..()
    ! PARAMETERS FOR THE NOSE-HOOVER THERMOSTATS
    IF (cntl%tnosep.AND.paral%parent) CALL nosepa(ipwalk,ipwalk)
    ! Dont symmetrize density 
    cntl%tsymrho=.FALSE.
    ! ..Make sure TKFULL=.TRUE
    IF (tkpts%tkpnt.AND.(.NOT.tkpts%tkfull).AND.(.NOT.cntl%tmdeh).AND.paral%io_parent) THEN
       WRITE(6,*)&
            ' ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(6,*)&
            ' WARNING! USE KPOINTS FULL  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(6,*)&
            ' ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    ENDIF
    IF (cntl%mimic) THEN
      CALL mimic_switch_dim(go_qm=.TRUE.)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == INITIALISATION                                               ==
    ! ==--------------------------------------------------------------==
    ! ..TSH[
    IF (tshl%tdtully) THEN
       ! ..XFMQC[
       IF (tshl%txfmqc.AND.tmw) THEN
          filew="SH_WAVEFUNCTIONS_"
          CALL mw_filename(filew,shfilen1,mwi%walker_id)
          filew="SH_LRWAVEFUNCTIONS_"
          CALL mw_filename(filew,shfilen2,mwi%walker_id)
          ! ..XFMQC]
       ELSE
          shfilen1="SH_WAVEFUNCTIONS"
          shfilen2="SH_LRWAVEFUNCTIONS"
       ENDIF
       IF (paral%io_parent) INQUIRE(file=shfilen1,exist=fexist)
       CALL mp_bcast(fexist,parai%io_source,parai%cp_grp)
       IF (.NOT.fexist) THEN
          tshi%shstep=0
          ! CALL DCOPY(2*NGW*NSHS,C0(1,1,1),1,CK(1,1),1)
       ELSE
          CALL tmprd(ck,tshi%nshs,1,tshi%shstep,shfilen1)
          CALL mp_bcast(tshi%shstep,parai%io_source,parai%cp_grp)
       ENDIF
       IF (tshl%tdextpot.OR.shlct%sh_lcontrol) THEN
          IF (tshl%tdextpot) THEN
             shfilen3="SH_EXTPT.dat"
          ELSEIF (shlct%sh_lcontrol) THEN
             shfilen3="SH_LCON.dat"
          ENDIF
          IF (paral%io_parent) THEN
             INQUIRE(file=shfilen3,exist=fexist)
             IF (.NOT.fexist) THEN
                !OPEN(unit=101,file=shfilen3,status='NEW')
                tshf%apot=0._real_8
                tshf%etpot=0._real_8
             ELSE
                OPEN(unit=101,file=shfilen3,status='UNKNOWN',&
                     position='APPEND')
                BACKSPACE(101)
                READ(101,'(I10,2X,2F15.6)') temp,tshf%apot,tshf%etpot
                CLOSE(101)
             ENDIF
          ENDIF
          CALL mp_bcast(tshf%apot,parai%io_source,parai%cp_grp)
          CALL mp_bcast(tshf%etpot,parai%io_source,parai%cp_grp)
       ENDIF
    ENDIF
    ! ..TSH]
    IF (cntl%bsymm.AND.paral%io_parent) THEN
       WRITE(6,*)
       WRITE (6,*) 'BSYMM: BS WAVEFUNCTION INITIALIZATION'
    ENDIF
    IF (cntl%cdft)CALL init_cdft()
    ! EHR[
    IF (cntl%tmdeh) THEN
       CALL initrun(irec,ch,c2,sc0,rhoe,psi,eigv)
    ELSE
       CALL initrun(irec,c0(:,:,1:1),c2,sc0,rhoe,psi,eigv)
    ENDIF
    ! EHR]
    IF (cntl%bsymm) THEN
       IF (paral%io_parent) WRITE (6,*) 'BSYMM: HS WAVEFUNCTION INITIALIZATION'
       bsclcs=2
       CALL setbsstate
       CALL initrun(irec,c0(:,:,2:2),c2,sc0,rhoe,psi,eigv(1,2))
    ENDIF
    IF (tshl%txfmqc) THEN ! MIN: iniailize buffers : tauall, bf_nli
       CALL xf_fill_all_buffers !MIN: call this subroutine in SH_TDDFT_UTILS
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Localize wavefunctions if HFX + screening
    IF (vdwl%vdwd) THEN
       IF (cntl%tpath) THEN
          ipx=ipcurr-np_low+1
       ELSE
          ipx=1
       ENDIF
       vdwwfl%trwannc=trwanncx(ipx)
    ENDIF
    IF (hfxc3%twscr.OR.(vdwl%vdwd.AND..NOT.vdwwfl%trwannc))&
         CALL localize2(tau0,c0,c2,sc0,crge%n)
    ! ==--------------------------------------------------------------==
    ! 
    IF (cntl%mimic) THEN
      CALL mimic_switch_dim(go_qm=.FALSE.)
    ENDIF
    CALL mp_bcast(taup,SIZE(taup),parai%io_source,parai%cp_grp)
    IF (paral%parent) CALL dcopy(3*maxsys%nax*maxsys%nsx,taup,1,taui,1)
    ! ==--------------------------------------------------------------==
    IF (cprint%iprint_step.EQ.0) cprint%iprint_step=cnti%nomore+1
    ! ==--------------------------------------------------------------==
    ! EHR[
    IF (cntl%tmdeh) THEN
       CALL initialize_ehrenfest_dyn(c0,ch,rhoe,psi,norms)
    ENDIF
    ! EHR]
    ! INITIALIZE VELOCITIES
    IF (paral%parent) CALL detdof(tau0,taur)
    IF (paral%io_parent.AND.cntl%new_constraints) THEN
       WRITE(6, '(1x,a)') "USING NEW CONSTRAINTS SOLVER"
       IF (cntl%pbicgstab) then
          WRITE(6, '(1x,a)') "WITH PRECONDITIONED BICONJUGATE GRADIENT STABILIZED"
       ELSE
          WRITE(6, '(1x,a)') "WITH PRECONDITIONED CONJUGATE GRADIENT"
       ENDIF
       WRITE(6, '(1x,a,t60,i6)') "MAX NUMBER OF SHAKE ITERATIONS:", &
                                  cnti%shake_maxstep
       WRITE(6, '(1x,a,t60,i6)') "MAX NUMBER OF PCG STEPS IN EACH SHAKE ITERATION:", &
                                  cnti%shake_cg_iter
       CALL init_constraints(ntcnst, cnsval, ions0%na, ions1%nsp)
       IF (.NOT.mimic_control%external_constraints) THEN
          WRITE(6, '(1x,a)') "EXTERNAL CONSTRAINTS DEFINED THROUGH MIMIC WILL BE IGNORED"
       ENDIF
       IF (cntl%mimic) THEN
          CALL mimic_subsystem_dof()
       END IF
    END IF

    ! INITIALIZE METADYNAMICS VARIABLES used also for 
    ! SHOOTING from SADDLE POINT with RANDOM VELOCITIES
    IF (tmw) CALL mp_sync(supergroup)
    IF (ionode .AND. (lmeta%lcolvardyn .OR. lmeta%lsadpnt)) THEN
       CALL colvar_structure(tau0,taur)
    ENDIF

    IF (irec(irec_vel).EQ.0.AND..NOT.restart1%rgeo) THEN
       ener_com%ecnstr = 0.0_real_8
       ener_com%erestr = 0.0_real_8
       CALL rinvel(velp,c2,nstate)
       IF (paral%parent) CALL taucl(velp)
       IF (paral%io_parent) THEN
          IF (cntl%new_constraints) THEN
             CALL do_rattle(tau0, velp)
          ELSE
             CALL rattle(tau0,velp)
          END IF
       END IF
       CALL rvscal(velp)
    ELSE
       IF (paral%parent) CALL taucl(velp)
       IF (paral%io_parent) THEN
          IF (cntl%new_constraints) THEN
             CALL do_rattle(tau0, velp)
          ELSE
             CALL rattle(tau0,velp)
          END IF
       END IF
       IF (cntl%trescale) CALL rvscal(velp)
    ENDIF
    IF (cntl%quenchp) CALL zeroing(velp)!,3*maxsys%nax*maxsys%nsx)
    IF (cntl%trevers) THEN
       ! invert ionic velocities (useful for path sampling)
       CALL dscal(3*maxsys%nax*maxsys%nsx,-1._real_8,velp,1)
    ENDIF
    ! COMPUTE THE IONIC TEMPERATURE TEMPP
    IF (paral%parent) THEN
       CALL ekinpp(ekinp,velp)
       tempp=ekinp*factem*2._real_8/glib
       ! ..TSH[
       tempp_sh=tempp
       ! ..TSH]
    ENDIF
    ! RESET ACCUMULATORS
    IF (paral%parent)THEN
       IF (irec(irec_ac).EQ.0)CALL resetac(tau0,taui,iteropt%nfi)
       ! 
       itemp=irec(irec_nop1)+irec(irec_nop2)+irec(irec_nop3)&
            +irec(irec_nop4)
       IF (cntl%tnosep .AND. itemp.EQ.0) CALL nospinit(ipwalk)
    ENDIF
    IF (cntl%mimic) THEN
      CALL mimic_switch_dim(go_qm=.TRUE.)
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
    nx = 1
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

    ! INITIALIZE cntl%cdft
    IF (cntl%cdft)THEN
       CALL cdft_w(rhoe,tau0,dummy)
       IF (cntl%cdft_weight)CALL write_w(wdiff,"INIT")
    ENDIF

    IF (paral%io_parent) THEN
       WRITE(6,'(1X,64("="))')
       WRITE(6,'(1X,"==",T25,A,T64,"==")')&
            'FORCES INITIALIZATION'
       WRITE(6,'(1X,64("="))')
    ENDIF
    IF (cntl%mimic) THEN
       CALL mimic_update_coords(tau0, c0, cm, nstate, ncpw%ngw, inyh)
       IF (paral%io_parent) THEN
          CALL mimic_ifc_send_coordinates()
          IF (mimic_control%tot_scf_energy) THEN
             CALL mimic_ifc_collect_energies()
          END IF
       END IF
    END IF
    ifcalc=0
    IF (cntl%bsymm) THEN
       CALL bs_forces_diag(nstate,c0,c2,cm,sc0,cm(nx:),vpp,eigv,&
          rhoe,psi,&
          tau0,velp,taui,fion,ifcalc,&
          irec,.TRUE.,.TRUE.)
    ELSE IF (cntl%tmdeh) THEN ! EHR[
       CALL forces_prop(nstate,c0,ch,c2,cm,sc0,cm(nx:),vpp,eigv,&
          rhoe,psi,&
          tau0,velp,taui,fion,ifcalc,&
          irec,.TRUE.,.TRUE.)
       ! EHR]

    ELSE IF (cntl%use_mts) then
       n_inner_steps=0
       n_large_steps=0
       call get_mts_forces(.true.,.true.,n_inner_steps,n_large_steps,&
          mts_pure_dft,fion_high,c0_high,cold_high,nnow_high,numcold_high,infi,&
          nstate,c0,c2,cm,sc0,cm(nx:),vpp,eigv,&
          rhoe,psi,tau0,velp,taui,fion,ifcalc,irec,.true.,.true.)

    else 
       CALL forces_diag(nstate,c0,c2,cm,sc0,cm(nx:),vpp,eigv,&
          rhoe,psi,tau0,velp,taui,fion,ifcalc,irec,.TRUE.,.TRUE.)

    ENDIF
    IF (cntl%tddft) THEN
       CALL lr_tddft(c0(:,:,1),c1,c2(:,:,1),sc0,rhoe,psi,tau0,fion,eigv,&
            nstate,.TRUE.,td01%ioutput)
    ENDIF
    ![ EXACT FACTORIZATION
    IF (tshl%tdtully.AND.tshl%txfmqc) THEN
       CALL xf_force_components(c0,c1,c2(:,:,1),sc0,rhoe,eigs)
    ENDIF
    !] EXACT FACTORIZATION 

    ! switch on info printing for the lin resp part
    IF (cntl%tresponse) dmbi%inter_pt_firstcall = .TRUE.

    ! Initialize Metadynamics contributions
    IF (lmeta%lcolvardyn .AND. lmeta%lextlagrange) THEN
       lmetares= .FALSE.
       ! Metadynamics with Extended Lagrangian
       IF (lmeta%tmulti) THEN
          CALL meta_ext_mul(tau0,velp,taur,&
               .FALSE.,lmetares,.FALSE.,0.0_real_8,ekinp)
       ELSE IF (tmw) THEN
          CALL meta_extlagr_mw(tau0,velp,taur,&
               .FALSE.,lmetares,.FALSE.,0.0_real_8,ekinp)
       ELSE
          CALL meta_extlagr(tau0,velp,taur,&
               .FALSE.,lmetares,.FALSE.,0.0_real_8,ekinp)
       ENDIF
    ENDIF

    IF (.NOT.cntl%bsymm) CALL dcopy(nnx,rin0,1,rm1(1),1)

    IF (cntl%mimic) THEN
       CALL mimic_sum_forces(fion)
    END IF

    IF (paral%io_parent) THEN
       CALL wrgeof(tau0,fion)
       CALL fileopen(3,filen,fo_app+fo_verb+fo_nofp,ferror)
    ENDIF
    IF (paral%io_parent) THEN
       WRITE(6,'(1X,64("="))')
       WRITE(6,'(1X,"==",T20,A,T64,"==")')&
            'END OF FORCES INITIALIZATION'
       WRITE(6,'(1X,64("="),/)')
    ENDIF

    ! INITIALIZE GLE THERMO
    CALL gle_init(tau0,velp,rmass%pma)

    CALL write_irec(irec)
    ! ==--------------------------------------------------------------==
    ! == END INITIALIZATION                                           ==
    ! ==--------------------------------------------------------------==
    IF (teststore(0).AND.cntl%tsampl)&
         CALL zhwwf(2,irec,c0,c2,nstate,eigv,tau0,velp,taui,iteropt%nfi)
    IF (paral%io_parent) THEN
       ! MEAN SQUARE DISPLACEMENT OF DIFFERENT IONIC SPECIES
       IF (paral%parent) CALL dispp(tau0,taui,disa)
       ! ENERGY OF THE NOSE THERMOSTATS
       CALL noseng(iteropt%nfi,velp,enose,enosp,dummy(1),1)
       IF (cntl%mimic) THEN
          mimic_energy%qm_energy = ener_com%etot
          ener_com%etot = ener_com%etot &
                          + ener_com%eext &
                          + mimic_energy%qmmm_energy &
                          + mimic_energy%mm_energy
       END IF
       econs=ekinp+ener_com%etot+enose+enosp+ener_com%ecnstr
       IF (cntl%cdft)econs=econs+cdftcom%cdft_v(1)*(cdftcom%cdft_nc+cdftcom%vgrad(1))
       time2 =m_walltime()
       tcpu = (time2 - time1)*0.001_real_8
       WRITE(6,'(A,T50,F8.2,A8)') ' TIME FOR INITIALIZATION:',&
            tcpu,' SECONDS'
       WRITE(6,'(//,1X,64("="))')
       WRITE(6,'(1X,"=",T20,A,T65,"=")')&
            'MOLECULAR DYNAMICS SIMULATION'
       WRITE(6,'(1X,64("="))')
       ! CALL WRPRINT_MD(EIGV,F,AMU,NSTATE,TAU0,FION,
       ! &                  0._real_8,TEMPP,ETOT,ECONS,0._real_8,DISA,
       ! &                  TCPU,.FALSE.,NFI,0)
    ENDIF
    ! ==================================================================
    ! ==          THE BASIC LOOP FOR MOLECULAR DYNAMICS               ==
    ! ==                 USING VELOCITY VERLET                        ==
    ! ==================================================================
#if defined(USE_IBM_HPM)
    CALL hpm_start('MD LOOP')
#endif
    infi=0
    n_inner_steps=0
    n_large_steps=0
    nfimin=iteropt%nfi+1
    nfimax=iteropt%nfi+cnti%nomore
    DO loopnfi=nfimin,nfimax
       time1=m_walltime()
       CALL mp_sync(parai%cp_grp)
       infi=infi+1

       ! MTS time step counters
       n_inner_steps=n_inner_steps+1
       mts_large_step=.false.
       ! check if it is a large step
       if (n_inner_steps==mts%timestep_factor) then
          mts_large_step=.true.
          n_large_steps=n_large_steps+1
          n_inner_steps=0
       endif

       ! ..TSH[
       IF (tshl%tdtully) tshi%shstep=tshi%shstep+1
       ! ..TSH]
       iteropt%nfi=loopnfi
       comvl%subcom=comvl%tsubcom.AND.MOD(iteropt%nfi-1,comvl%ncomv).EQ.0
       comvl%subrot=comvl%tsubrot.AND.MOD(iteropt%nfi-1,comvl%nrotv).EQ.0
       tstrng=cntl%tpmin.AND.MOD(iteropt%nfi,cnti%nomore).EQ.0
       IF (hubbu%pfrqom.gt.0) THEN
          hubbu%tpom=MOD(iteropt%nfi-1,hubbu%pfrqom).EQ.0
       ELSE
          hubbu%tpom=.False.
       ENDIF
       IF (cntl%mimic) THEN
          CALL mimic_switch_dim(go_qm=.FALSE.)
       ENDIF
       ! ANNEALING
       ! thermostats only if 1) standard dynamics 2) larger MTS step
       if ( .not.cntl%use_mts .or. mts_large_step ) then
          CALL anneal(velp,c2,nstate,scr)
          CALL berendsen(velp,c2,nstate,scr,0.0_real_8,0.0_real_8)
          ! UPDATE NOSE THERMOSTATS
          CALL noseup(velp,c2,nstate,ipwalk)
          ! FIRST HALF OF GLE EVOLUTION
          CALL gle_step(tau0,velp,rmass%pma)
       endif 
       ! SUBTRACT CENTER OF MASS VELOCITY
       IF (paral%io_parent.AND.comvl%subcom) CALL comvel(velp,vcmio,.TRUE.)
       ! SUBTRACT ROTATION AROUND CENTER OF MASS
       IF (paral%io_parent.AND.comvl%subrot) CALL rotvel(tau0,velp,lmio,&
            tauio,.TRUE.)
       ! UPDATE VELOCITIES
       IF (paral%io_parent) CALL velupi(velp,fion,1)
#if defined (__QMECHCOUPL)
       IF (paral%io_parent .AND. qmmech) THEN
          CALL mm_cpmd_velup(cntr%delt_ions)
       ENDIF
#endif
       ! UPDATE POSITIONS
       ener_com%ecnstr = 0.0_real_8
       IF (paral%io_parent) THEN
          CALL posupi(tau0,taup,velp)
          IF (cntl%new_constraints) THEN
             CALL do_shake(tau0, taup, velp)
          ELSE
             IF (cotc0%mcnstr.NE.0) CALL cpmdshake(tau0,taup,velp)
          END IF
#if defined (__QMECHCOUPL)
          IF (qmmech) THEN
             CALL mm_cpmd_update_links(taup, ions0%na, ions1%nsp, maxsys%nax, ions1%nat)
             CALL mm_cpmd_posup(cntr%delt_ions)
          ENDIF
#endif
       ENDIF
       CALL mp_bcast(taup,SIZE(taup),parai%io_source,parai%cp_grp)
       IF (cntl%mimic) THEN
          CALL mimic_update_coords(taup, c0, cm, nstate, ncpw%ngw, inyh)
          mimic_control%update_potential = .TRUE.
          IF (paral%io_parent) THEN
             CALL mimic_ifc_send_coordinates()
             IF (mimic_control%tot_scf_energy) THEN
                CALL mimic_ifc_collect_energies()
             END IF
          END IF
          IF (mimic_control%long_range_coupling.AND.(mod(iteropt%nfi-1,mimic_control%update_sorting).EQ.0)) THEN
             CALL mimic_ifc_sort_fragments()
          END IF
          CALL mimic_switch_dim(go_qm=.TRUE.)
       END IF
       ! Reset swap file to update (option BLOCK and ALL)
       IF (tkpts%tkblock) CALL reskpt_swap
       CALL phfac(taup)
       IF (corel%tinlc) CALL copot(rhoe,psi,ropt_mod%calste)
       IF (tmovr) THEN
          CALL dcopy(nnx,rin0,1,rhoe(1,1),1)
          CALL moverho(rhoe,psi)
          CALL dcopy(nnx,rhoe(1,1),1,rin0,1)
       ENDIF
       IF (fint1%ttrot) THEN
          CALL calc_alm
       ENDIF
       ! CALCULATE THE FORCES
       ropt_mod%calste=cntl%tpres.AND.MOD(iteropt%nfi,cnti%npres).EQ.0
       IF (cntl%textrap) THEN
          ! Extrapolate wavefunctions
          CALL extrapwf(infi,c0,scr,cold,nnow,numcold,nstate,cnti%mextra)
       ENDIF
       IF (cntl%tlanc) nx=1
       IF (cntl%tdavi) nx=cnti%ndavv*nkpt%nkpnt+1
       ! RESPONSE calculation
       IF (cntl%tresponse) THEN
          ! save the cntl%md step number
          save_nfi=iteropt%nfi
          ! localisation at first iteration + PT
          CALL mddiag_interaction_p(nstate,c0,c2,cm,sc0,&
               CM(NX:),VPP,EIGV,RHOE,PSI,&!SCR,LSCR,&
               TAUP,VELP,TAUI,FION,IFCALC,IREC,.TRUE.,.FALSE.)
          ! sets the number of iteration and recover cntl%md step number
          ifcalc=iteropt%nfi
          iteropt%nfi=save_nfi

          ! switch off the info printing
          dmbi%inter_pt_firstcall=.FALSE.
       ELSE 
          ifcalc = 0
          IF (cntl%bsymm) THEN
             CALL bs_forces_diag(nstate,c0,c2,cm,sc0,cm(nx:),vpp,eigv,&
                rhoe,psi,&
                taup,velp,taui,fion,ifcalc,&
                irec,.TRUE.,.FALSE.)
          ELSE IF (cntl%tmdeh) THEN
             CALL forces_prop(nstate,c0,ch,c2,cm,sc0,cm(nx:),vpp,eigv,&
                rhoe,psi,&
                tau0,velp,taui,fion,ifcalc,&
                irec,.TRUE.,.TRUE.)

          ELSE IF (cntl%use_mts) then
             call get_mts_forces(mts_large_step,.false.,n_inner_steps,n_large_steps,&
                mts_pure_dft,fion_high,c0_high,cold_high,nnow_high,numcold_high,infi,&
                nstate,c0,c2,cm,sc0,cm(nx:),vpp,eigv,&
                rhoe,psi,taup,velp,taui,fion,ifcalc,irec,.true.,.false.)

          ELSE
             CALL forces_diag(nstate,c0,c2,cm,sc0,cm(nx:),vpp,eigv,&
                rhoe,psi,taup,velp,taui,fion,ifcalc,irec,.TRUE.,.FALSE.)
          ENDIF
          IF (cntl%tddft) THEN
             CALL lr_tddft(c0(:,:,1),c1,c2(:,:,1),sc0,rhoe,psi,taup,fion,eigv,&
                  nstate,.TRUE.,td01%ioutput)
          ENDIF
       ENDIF
       ! ..TSH[
       IF (tshl%tdtully) THEN 
          ![ EXACT FACTORIZATION
          IF (tshl%txfmqc) THEN
             CALL xf_force_components(c0,c1,c2(:,:,1),sc0,rhoe,eigs)
          ENDIF
          !] EXACT FACTORIZATION

          IF (tshl%tully_sh.AND.(.NOT.tshl%sh_phex)) THEN
             IF (paral%io_parent)&
                  WRITE(6,*) ' velocity rescaling after SH'
             IF (paral%io_parent) CALL rscvp_sh(dt_sh,tempp,velp)
             tshl%tully_sh=.FALSE.
          ENDIF
          IF (shlct%sh_lcontrol) THEN
             IF(paral%io_parent) CALL sh_lct
             CALL mp_bcast(tshf%apot,parai%io_source,parai%cp_grp)
             CALL mp_bcast(tshf%etpot,parai%io_source,parai%cp_grp)
          ELSEIF (tshl%tdextpot) THEN
             IF (paral%io_parent) CALL sh_extpot
             CALL mp_bcast(tshf%apot,parai%io_source,parai%cp_grp)
             CALL mp_bcast(tshf%etpot,parai%io_source,parai%cp_grp)
          ENDIF
       ENDIF
       ! ..TSH]
       ! PATH MINIMIZATION
       ! ==================================================================
       ! Damped Dynamics
       CALL dampdyn(velp,fion,cm,c2,nstate,scr(1),scr(10))
       ! ==================================================================
       ! Meta Dynamics of Collective Variables

       IF (lmeta%lcolvardyn) THEN
          lmetares= .FALSE.

          IF (lmeta%lextlagrange) THEN
             ! Metadynamics with Extended Lagrangian
             IF (lmeta%tmulti) THEN
                CALL meta_ext_mul(taup,velp,taur,&
                     .FALSE.,lmetares,.FALSE.,0.0_real_8,ekinp)
             ELSE IF (tmw) THEN
                CALL meta_extlagr_mw(taup,velp,taur,&
                     .FALSE.,lmetares,.FALSE.,0.0_real_8,ekinp)
             ELSE
                CALL meta_extlagr(taup,velp,taur,&
                     .FALSE.,lmetares,.FALSE.,0.0_real_8,ekinp)
             ENDIF
          ELSE
             ! Time dipendent potential applied directly on the Collective Variables
             IF (tmw) THEN
                CALL meta_colvar_mw(taup,velp,fion,taur,&
                     .FALSE.,lmetares,0.0_real_8,ekinp)
             ELSE
                CALL meta_colvar(taup,velp,fion,taur,&
                     .FALSE.,lmetares,0.0_real_8,ekinp)
             END IF
          ENDIF

          IF (paral%io_parent) THEN
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
       IF (cntl%mimic) THEN
          CALL mimic_sum_forces(fion)
       END IF
       ! ==================================================================

       IF (ropt_mod%calste) CALL totstr
#if defined (__QMECHCOUPL)
       ! ----------------------------------------------
       ! QMMM electrostatic coupling
       ! ----------------------------------------------
       IF (qmmech) THEN
          CALL mm_cpmd_ext_pot_f77
          CALL mm_cpmd_elstat(taup,ions0%na,ions1%nsp,maxsys%nax,ions1%nat,c0,scr)
       ENDIF
#endif
       IF (cntl%mimic) THEN
          CALL mimic_switch_dim(go_qm=.FALSE.)
       ENDIF
       ! FINAL UPDATE FOR VELOCITIES
       ener_com%ecnstr = 0.0_real_8
       IF (paral%io_parent) THEN
          CALL velupi(velp,fion,1)
#if defined (__QMECHCOUPL)
          IF (qmmech) THEN
             CALL mm_cpmd_velup(cntr%delt_ions)
          ENDIF
#endif
          IF (isos1%twall) CALL box_boundary(taup,velp)
          IF (cntl%new_constraints) THEN
             CALL do_rattle(taup, velp)
          ELSE
             CALL rattle(taup,velp)
          END IF
       ENDIF
       IF (paral%parent.AND..NOT.cntl%mimic) CALL geofile(taup,velp,'WRITE')
       ! COMPUTE THE IONIC TEMPERATURE TEMPP
       IF (paral%io_parent) THEN
          CALL ekinpp(ekinp,velp)
          IF (lmeta%lextlagrange.AND. ltcglobal) THEN
             CALL ekincv_global(ek_cv)
             tempp=(ek_cv+ekinp)*factem*2._real_8/(glib+REAL(ncolvar,kind=real_8))
          ELSE
             tempp=ekinp*factem*2._real_8/glib
          ENDIF
          IF (cntl%mimic) THEN
             CALL mimic_subsystem_temperatures(velp)
          END IF
          ! ..TSH[
          tempp_sh=tempp
          ! ..TSH]
#if defined (__QMECHCOUPL)
          ! -------------------------------------------------
          ! TEMPERTURE and KINETIC ENERGY of the MM subsystem
          ! MM_EKIN + EKINP = KE of whole QM/MM system
          ! -------------------------------------------------
          IF (qmmech) THEN
             CALL mm_cpmd_ekin(mm_ekin,mm_temp)
          ENDIF
#endif
       ENDIF
#if defined (__QMECHCOUPL)
       ! ---------------------------------------------------------------
       ! MM Temperature control - annealing schedule not yet implemented
       ! ---------------------------------------------------------------
       IF (paral%io_parent .AND. qmmech) THEN
          CALL mm_cpmd_temp_control(temp1,temp2,mm_temp,cntl%tcp)
          CALL mm_cpmd_ekin(mm_ekin,mm_temp)
       ENDIF
#endif
       ! IONIC TEMPERATURE CONTROL
       IF (paral%io_parent) CALL rscvp(temp1,temp2,tempp,velp)
       ! SUBTRACT ROTATION AROUND CENTER OF MASS
       IF (paral%io_parent.AND.comvl%subrot) CALL rotvel(tau0,velp,lmio,&
            tauio,.FALSE.)
       ! SUBTRACT CENTER OF MASS VELOCITY
       IF (paral%io_parent.AND.comvl%subcom) CALL comvel(velp,vcmio,.FALSE.)

       ! thermostats only if 1) standard dynamics 2) larger MTS step
       if ( .not.cntl%use_mts .or. mts_large_step ) then
          ! SECOND HALF OF GLE EVOLUTION
          CALL gle_step(tau0,velp,rmass%pma)

          ! UPDATE NOSE THERMOSTATS
          CALL noseup(velp,c2,nstate,ipwalk)
          CALL berendsen(velp,c2,nstate,scr,0.0_real_8,0.0_real_8)
          ! ANNEALING
          CALL anneal(velp,c2,nstate,scr)
       endif
       IF (paral%io_parent) THEN
          CALL ekinpp(ekinp,velp)
          IF (lmeta%lextlagrange.AND. ltcglobal) THEN
             CALL ekincv_global(ek_cv)
             tempp=(ek_cv+ekinp)*factem*2._real_8/(glib+REAL(ncolvar,kind=real_8))
          ELSE
             tempp=ekinp*factem*2._real_8/glib
          ENDIF
       ENDIF
       ! MEAN SQUARE DISPLACEMENT OF DIFFERENT IONIC SPECIES
       IF (paral%io_parent) CALL dispp(taup,taui,disa)
       ! ENERGY OF THE NOSE THERMOSTATS
       IF (paral%io_parent) CALL noseng(iteropt%nfi,velp,enose,enosp,dummy(1),ipwalk)
       ! CALCULATE PROPERTIES DURING SIMULATION.
       cntl%caldip=cntl%tdipd.AND.MOD(iteropt%nfi-1,cnti%npdip).EQ.0
       ! mb - Wannier stuff for vdW-WC; call ddipo instead of localize2
       IF (vdwl%vdwd) THEN
          IF (cntl%tpath) THEN
             ipx=ipcurr-np_low+1
          ELSE
             ipx=1
          ENDIF
          vdwwfl%twannup=twannupx(ipx).OR.(infi.EQ.1.AND..NOT.vdwwfl%trwannc)
          twannupx(ipx)=vdwwfl%twannup
       ENDIF
       IF (.NOT.cntl%bsymm)&
            CALL propcal(c0,c2(:,:,1),cm,sc0,taup,eigv,crge%f,ener_com%amu,&
            rhoe,psi,nstate,nkpt%nkpnt,iteropt%nfi,infi)
       ! PRINTOUT the evolution of the accumulators every time step
       IF (paral%io_parent) THEN
          econs=ekinp+ener_com%etot+enose+enosp+ener_com%ecnstr+ekincv+vharm+glepar%egle 
          IF (cntl%cdft)THEN
             econs=econs+cdftcom%cdft_v(1)*(cdftcom%cdft_nc+cdftcom%vgrad(1))
          ENDIF
#if defined (__QMECHCOUPL)
          IF (qmmech) THEN
             econs=econs + mm_ekin
             IF (paral%io_parent)&
                  WRITE (6,'(50x, "TEMP: ",f10.0)') mm_temp
          ENDIF
#endif
          time2=m_walltime()
          tcpu=(time2-time1)*0.001_real_8
          CALL wrprint_md(eigv,crge%f,ener_com%amu,nstate,taup,fion,&
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
          CALL printp(taur,taup,velp)
          IF (cprint%twriteforcetrajectory) CALL printp2(taur,taup,velp,fion)
#if defined (__QMECHCOUPL)
          ! --------------------------------------------------------
          ! ..       Write to QMMM trajectory.  This trajectory file contains
          ! the atoms that make up the "real" system. (no dummies)
          ! --------------------------------------------------------
          IF (qmmech) THEN
             IF (rout1%rout .AND. MOD(iteropt%nfi-1,cnti%ntraj).EQ.0) THEN
                CALL mm_cpmd_write_trajectory(iteropt%nfi)
             ENDIF
          ENDIF
#endif
          time2=m_walltime()
          tcpu=(time2-time1)*0.001_real_8
       ENDIF
       IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)

       IF (tstrng) THEN
          soft_com%exsoft=.TRUE.
          soft_com%exnomore=.TRUE.
       ELSE IF (iteropt%nfi.EQ.nfimax) THEN
          soft_com%exsoft=.TRUE.
       ENDIF
       ! EHR[
       IF (cntl%tmdeh)&
            CALL zhwwf(2,irec,c0,c2,nstate,eigv,taup,velp,taui,iteropt%nfi)
       ! EHR]
       ! periodic output of density/wavefunction etc.
       IF (rout1%rhoout.AND.(rout1%nrhoout.GT.0)) THEN
          IF (MOD(iteropt%nfi-1,rout1%nrhoout).EQ.0) THEN
             CALL rhopri(c0,tau0,rhoe,psi(:,1),nstate,nkpt%nkpnt)
          ENDIF
       ENDIF
       IF (teststore(iteropt%nfi).OR.soft_com%exsoft.OR.lmetares) THEN
          CALL zhwwf(2,irec,c0,c2,nstate,eigv,taup,velp,taui,iteropt%nfi)
       ENDIF

       IF (soft_com%exsoft .AND.lmeta%lcolvardyn) THEN
          lmetares= .TRUE.
          IF (lmeta%lextlagrange) THEN
             ! Metadynamics with Extended Lagrangian
             IF (lmeta%tmulti) THEN
                CALL meta_ext_mul(taup,velp,taur,&
                     .FALSE.,lmetares,.FALSE.,0.0_real_8,ekinp)
             ELSE IF (tmw)THEN
                CALL meta_extlagr_mw(taup,velp,taur,&
                     .FALSE.,lmetares,.FALSE.,0.0_real_8,ekinp)
             ELSE
                CALL meta_extlagr(taup,velp,taur,&
                     .FALSE.,lmetares,.FALSE.,0.0_real_8,ekinp)
             ENDIF
          ELSE
             ! Time dependent potential applied directly on the Collective Variables
             IF (tmw) THEN
                CALL meta_colvar_mw(taup,velp,fion,taur,&
                     .FALSE.,lmetares,0.0_real_8,ekinp)
             ELSE
                CALL meta_colvar(taup,velp,fion,taur,&
                     .FALSE.,lmetares,0.0_real_8,ekinp)
             END IF
          ENDIF
       ENDIF
#if defined (__QMECHCOUPL)
       ! -----------------------
       ! Write QMMM restart file
       ! -----------------------
       IF (paral%io_parent.AND.qmmech) THEN
          IF (MOD(iteropt%nfi,store1%istore).EQ.0.OR.iteropt%nfi.EQ.nfimax.OR.soft_com%exsoft) THEN
             CALL mm_cpmd_write_restart
          ENDIF
       ENDIF
#endif
       IF (cntl%mimic) THEN
          CALL mimic_switch_dim(go_qm=.FALSE.)
       ENDIF
       ! temperature ramping
       CALL tempramp(temp1,temp2)
       ! UPDATE IONIC POSITIONS
       CALL dcopy(3*maxsys%nax*maxsys%nsx,taup(1,1,1),1,tau0(1,1,1),1)
       ! UPDATE DENSITY
       IF (.NOT.cntl%bsymm) THEN
          CALL extrap(nnx,andr2%alxmix,rm1,rin0,rinp)
          CALL dcopy(nnx,rin0,1,rm1(1),1)
          CALL dcopy(nnx,rinp(1),1,rin0,1)
       ENDIF
       IF (cntl%mimic) THEN
          CALL mimic_switch_dim(go_qm=.TRUE.)
       ENDIF
       ! STOP THE RUN IF THE USER HAS SET THE SIGNAL 30
       ! ..TSH[
       IF (tshl%tdtully) THEN
          CALL tmpwr(ck,tshi%nshs,1,tshi%shstep,shfilen1)
          CALL dcopy(2*ncpw%ngw*tshi%nshs,c0(1,1,1),1,ck(1,1),1)
          IF ((tshl%tdextpot.OR.shlct%sh_lcontrol).AND.paral%io_parent) THEN
             fexist=.FALSE.
             INQUIRE(file=shfilen3,exist=fexist)
             IF (.NOT.fexist) THEN
                OPEN(unit=101,file=shfilen3,status='NEW')
             ELSE
                OPEN(unit=101,file=shfilen3,status='UNKNOWN',&
                     position='APPEND')
             ENDIF
             WRITE(101,'(I10,2X,2F15.6)') tshi%shstep,tshf%apot,tshf%etpot
             CLOSE(101)
          ENDIF
       ENDIF
       ! ..TSH]
       ! EHR[
       IF (cntl%tmdeh.AND.store1%srho.AND.&
            ((infi.EQ.1).OR.(MOD(infi,store1%istore).EQ.0))) THEN
          CALL rhopri(c0,tau0,rhoe,psi(:,1),crge%n,1)
       ENDIF
       ! EHR]
       IF (soft_com%exsoft) GOTO 100
       IF (cntl%cdft)CALL cdft_reinit()
       ! ISC[
       IF (tshl%isc) THEN
          CALL dcopy(3*maxsys%nax*maxsys%nsx,tau0,1,tausoc,1)
          IF (MOD(infi,nskip).EQ.0.OR.do_soc_because_sh) THEN
             CALL do_tsh_isc(c0(:,:,1),c1,c2(:,:,1),sc0,rhoe,psi,eigv,ener_com%etot,& 
                  soc_el,soc_array,infi,nstate,nroot)
             CALL do_tsh_isc_lz(soc_array,cntr%delt_elec,nroot,nskip)
          ENDIF
          IF (do_soc_because_sh) do_soc_because_sh=.FALSE.
          CALL isc_write_start
       ENDIF
       ! ISC]
       ! [EXACT FACTORIZATION - here we fill the buffer for the coordinates
       !                        Each walker contributes with its coordinates
       IF (tshl%txfmqc) THEN
          CALL xf_fill_all_buffers !MIN: call this subroutine in SH_TDDFT_UTILS
       ENDIF
       ! ]EXACT FACTORIZATION
       ! ==================================================================
       ! ==     END OF MAIN LOOP                                         ==
       ! ==================================================================
    ENDDO
    IF (paral%io_parent) THEN
       WRITE(6,'(1X,64("="))')
       WRITE(6,'(1X,"=",T17,A,T65,"=")')&
            'END OF MOLECULAR DYNAMICS SIMULATION'
       WRITE(6,'(1X,64("="),/,/)')
    ENDIF
    ! ==--------------------------------------------------------------==
100 CONTINUE
    ! This is just to remove the EXIT file..
    IF (tmw) CALL testex_mw(soft_com%exsoft)

#if defined(USE_IBM_HPM)
    CALL hpm_stop('MD LOOP')
#endif

    IF (cntl%mimic) THEN
       CALL mimic_update_coords(taup, c0, cm, nstate, ncpw%ngw, inyh)
       IF (paral%io_parent) THEN
          CALL mimic_ifc_send_coordinates()
       END IF
       IF (mimic_control%long_range_coupling) THEN
          CALL mimic_ifc_sort_fragments()
       END IF
    END IF

    IF (cntl%cdft)THEN
       IF (cntl%cdft_weight)CALL write_w(wdiff,"FINAL")
    ENDIF

    !IF (.NOT.cntl%tmdeh.AND.rout1%rhoout.AND.(rout1%nrhoout.LE.0)) THEN
    IF (rout1%rhoout.AND.(rout1%nrhoout.LE.0)) THEN
       CALL rhopri(c0(:,:,:),tau0,rhoe,psi(:,1),nstate,nkpt%nkpnt)
       IF (cntl%bsymm) THEN
          bsclcs=2
          CALL rhopri(c0(:,:,2:),tau0,rhoe,psi(:,1),nstate,nkpt%nkpnt)
          bsclcs=1
       ENDIF
    ENDIF
    ! Print accumulators.
    IF (paral%io_parent) CALL paccc(tempp,ener_com%etot,econs,enose,enosp,&
         ener_com%ecnstr,ener_com%erestr,ener_com%ebogo,disa,tcpu,iteropt%nfi,0)
    IF (paral%io_parent) CALL gsize(fion,gnmax,gnorm)
    IF (paral%io_parent) THEN
       IF (.NOT.cntl%bsymm) THEN
          CALL finalp(tau0,fion,velp,eigv)
       ELSE
          CALL wrgeof(tau0,fion)
          WRITE(6,'(A)') ' NUCLEAR GRADIENT:'
          WRITE(6,'(2(A,1PE15.5))') '    MAX. COMPONENT =',&
               gnmax,'         NORM =',gnorM
       ENDIF
    ENDIF
    IF (cntl%tsampl) THEN
       CALL sample_go
       GOTO 99999
    ENDIF
10000 CONTINUE
    ! ==--------------------------------------------------------------==
    DEALLOCATE(eigv,STAT=ierr)
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
    IF (tshl%txfmqc) THEN
       DEALLOCATE(xfmqc%fion_state,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(eigs,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(xfmqc%qmoment,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(xfmqc%tauall,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(xfmqc%w_ij,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       ![FEDE
       DEALLOCATE(xfmqc%sigma,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       !FEDE]
       DEALLOCATE(xfmqc%k_li,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(xfmqc%f_nli,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(xfmqc%r0_ni,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(xfmqc%d0_ni,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(xfmqc%bf_nli,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(xfmqc%bfion_state,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(xfmqc%bfa_ni,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(xfmqc%fa_ni,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       ! .. MIN[
       DEALLOCATE(xfmqc%bw0_nkli,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(xfmqc%brho_l,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       ! .. MIN]
       DEALLOCATE(xfmqc%cf,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(xfmqc%bcf,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    DEALLOCATE(taup,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (paral%io_parent) CALL fileclose(3)
    DEALLOCATE(rhoe,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(psi,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(scr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (cntl%textrap) DEALLOCATE(cold,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    if (cntl%use_mts) then
       deallocate(fion_high,stat=ierr)
       if(ierr/=0) call stopgm(proceduren,'deallocation problem: fion_high',&
          __LINE__,__FILE__)
       if (mts_pure_dft) then
          deallocate(c0_high,stat=ierr)
          if(ierr/=0) call stopgm(proceduren,'deallocation problem: c0_high',&
             __LINE__,__FILE__)
          if (cntl%textrap) then 
             deallocate(cold_high,stat=ierr)
             if(ierr/=0) call stopgm(proceduren,'deallocation problem: cold_high',&
                __LINE__,__FILE__)
          end if
       end if
    endif 
    IF (comvl%tsubrot) DEALLOCATE(tauio,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (lmeta%lcolvardyn) DEALLOCATE(fhills,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (.NOT.cntl%bsymm) THEN
       !      DEALLOCATE(rinp,STAT=ierr)
       !      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
       !           __LINE__,__FILE__)
       !      DEALLOCATE(rm1,STAT=ierr)
       !      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
       !           __LINE__,__FILE__)
       !      DEALLOCATE(rmix,STAT=ierr)
       !      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
       !           __LINE__,__FILE__)
       !      DEALLOCATE(rout0,STAT=ierr)
       !      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
       !           __LINE__,__FILE__)
       !      DEALLOCATE(rin0,STAT=ierr)
       !      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
       !           __LINE__,__FILE__)
    ENDIF
    ! ..TSH[
    IF (tshl%tdtully) THEN
       DEALLOCATE(c_sh,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(v_sh,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(shsigma,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(shsigmaf,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(cdotlct,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       IF (tshl%txfmqc) THEN
          DEALLOCATE(xfmqc%eigv,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
          DEALLOCATE(xfmqc%nacv,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
       ENDIF
    ENDIF
    ! ..TSH]
    ! EHR[
    IF (cntl%tmdeh) THEN
       DEALLOCATE(norms,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(extf,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(ch,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
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
    IF (cntl%mimic) THEN
      CALL mimic_revert_dim()
    ENDIF
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE mddiag


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
     complex(real_8), allocatable :: c0_high(:,:,:), cold_high(:,:,:,:)
     complex(real_8) :: c0(:,:,:), c2(nkpt%ngwk,nstate), cr(*)
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
        call get_external_forces('EXT_LOW_FORCES', tau0, fion)

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
           call get_external_forces('EXT_HIGH_FORCES', tau0, fion_high)

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


  ! ==================================================================
  SUBROUTINE give_scr_mddiag(lmddiag,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lmddiag
    CHARACTER(len=30)                        :: tag

    INTEGER :: lcalc_alm, lcopot, lddipo, lforces_diag, linitrun, lmoverho, &
      lmtd, lpropcal, lrhoofr, lrhopri, lrinitwf, lrnlsm, ltddft, nstate

    nstate=crge%n
    lrnlsm=0
    lcalc_alm=0
    lcopot=0
    lrhopri=0
    lmoverho=0
    ltddft=0
    linitrun=0
    lddipo=0
    CALL give_scr_initrun(linitrun,tag)
    CALL give_scr_rinitwf(lrinitwf,tag,nstate)
    IF (pslo_com%tivan) CALL give_scr_rnlsm(lrnlsm,tag,nstate,.FALSE.)
    CALL give_scr_rhoofr(lrhoofr,tag)
    IF (fint1%ttrot) CALL give_scr_calc_alm(lcalc_alm,tag)
    IF (cntl%tmdeh) THEN
       CALL give_scr_forces_prop(lforces_diag,tag,nstate,.TRUE.)
    ELSE
       CALL give_scr_forces_diag(lforces_diag,tag,nstate,.TRUE.)
    ENDIF
    IF (corel%tinlc) CALL give_scr_copot(lcopot,tag)
    IF (rout1%rhoout) CALL give_scr_rhopri(lrhopri,tag,nstate)
    CALL give_scr_propcal(lpropcal,tag,nstate)
    IF (tmovr) CALL give_scr_moverho(lmoverho,tag)
    IF (cntl%tddft) CALL give_scr_lr_tddft(ltddft,.TRUE.,tag)
    CALL give_scr_meta_extlagr(lmtd,tag)
    IF (vdwl%vdwd) CALL give_scr_ddipo(lddipo,tag)
    lmddiag=MAX(lrinitwf,lrnlsm,lrhoofr,lforces_diag,ltddft,linitrun,&
         lcopot,lcalc_alm,lrhopri,lpropcal,lmoverho,lmtd,&
         lddipo)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE give_scr_mddiag

END MODULE md_driver

! ==================================================================
SUBROUTINE extrapwf(infi,c0,gam,cold,nnow,numcold,nstate,m)
  USE dotp_utils, ONLY: dotp
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE mp_interface, ONLY: mp_sum
  USE system , ONLY:cntl,ncpw,nkpt
  USE parac, ONLY : paral,parai
  USE kpts, ONLY : tkpts
  USE kpnt, ONLY : wk
  USE spin, ONLY : spin_mod
  USE elct, ONLY : crge
  USE geq0mod , ONLY:geq0
  USE pslo, ONLY : pslo_com
  USE legendre_p_utils, ONLY : d_binom
  USE ortho_utils, ONLY: ortho,preortho
  USE ovlap_utils, ONLY : ovlap
  USE rnlsm_utils, ONLY: rnlsm
  USE rotate_utils, ONLY: rotate
  USE utils, ONLY: zclean
  USE zeroing_utils,                   ONLY: zeroing
  IMPLICIT NONE
  INTEGER                                    :: nstate, infi, nnow, numcold, m
  COMPLEX(real_8)                            :: cold(nkpt%ngwk,nstate,nkpt%nkpnt,*)
  REAL(real_8)                               :: gam(nstate,*)
  COMPLEX(real_8)                            :: c0(nkpt%ngwk,nstate,*)

  INTEGER                                    :: i, ik, isub, ma, n1, nm, ierr
  REAL(real_8)                               :: fa, rsum, scalef
  REAL(real_8), EXTERNAL                     :: ddot
  complex(real_8), allocatable               :: c2(:,:)

! IF (pslo_com%tivan) CALL stopgm('EXTRAPWF','EXTRAPOLATION NOT (YET) ' //&
!      'SUPPORTED FOR VANDERBILT PSEUDOPOTENTIALS',& 
!      __LINE__,__FILE__)
  CALL tiset('  EXTRAPWF',isub)

  ALLOCATE(c2(nkpt%ngwk,nstate),stat=ierr)
  IF(ierr/=0) CALL stopgm('extrapwf','allocation problem : c2',&
       __LINE__,__FILE__)

  nnow=MOD(nnow,m)+1
  CALL dcopy(2*nkpt%ngwk*nstate*nkpt%nkpnt,c0,1,cold(1,1,1,nnow),1)

  numcold=MIN(numcold+1,m)

  ma=numcold
  IF (ma.GT.1) THEN
     CALL zeroing(c0(:,:,1:nkpt%nkpnt))!,nkpt%ngwk*nstate*nkpt%nkpnt)
     n1=nnow
     fa=-1.0_real_8
     DO i=1,ma
        nm=nnow-i+1
        IF (nm.LE.0) nm=nm+m
        ! construct extrapolation polynomial coefficient.
        IF (cntl%taspc) THEN
           fa=(-1._real_8)**(i+1)*REAL(i,kind=real_8)*d_binom(2*ma,ma-i)&
                /d_binom(2*ma-2,ma-1)
        ELSE
           fa=-1._real_8*fa*REAL(ma-i+1,kind=real_8)/REAL(i,kind=real_8)
        ENDIF
        IF (tkpts%tkpnt) THEN
           DO ik=1,nkpt%nkpnt
              CALL ovlap_c(nstate,gam,cold(:,:,ik,nm),cold(:,:,ik,n1))
              CALL mp_sum(gam,2*nstate*nstate,parai%allgrp)
              CALL rotate_c(CMPLX(fa,0._real_8,kind=real_8),cold(1,1,ik,nm),&
                   CMPLX(1._real_8,0._real_8,kind=real_8),c0(1,1,ik),gam,nstate)
           ENDDO
        ELSE
           CALL ovlap(nstate,gam,cold(:,:,1,nm),cold(:,:,1,n1))
           CALL mp_sum(gam,nstate*nstate,parai%allgrp)
           CALL rotate(fa,cold(:,:,1,nm),1._real_8,c0(:,:,1),gam,nstate,2*nkpt%ngwk,&
                cntl%tlsd,spin_mod%nsup,spin_mod%nsdown)
        ENDIF
     ENDDO
     ! orthogonalization: as done in k_forces
     DO ik=1,nkpt%nkpnt
        IF (cntl%nonort) THEN
           IF (geq0) CALL zclean(c0,nstate,ncpw%ngw)
        ELSE
           CALL preortho(c0(1,1,ik),nstate)
           IF (pslo_com%tivan) CALL rnlsm(c0(:,:,1),nstate,1,1,.FALSE.)
           CALL ortho(nstate,c0(:,:,ik),c2)
        ENDIF
     ENDDO
     IF (pslo_com%tivan) GOTO 1
     ! ..COUNT NUMBER OF ELECTRONS IN EXTRAPOLATED WFN AND RENORMALIZE
     rsum=0._real_8
     IF (tkpts%tkpnt) THEN
        DO ik=1,nkpt%nkpnt
           DO i=1,nstate
              IF (crge%f(i,ik).NE.0._real_8) THEN
                 rsum=rsum+&
                      crge%f(i,ik)*wk(ik) *&
                      ddot(nkpt%ngwk*2,c0(1,i, ik),1,c0(1,i,ik),1)
              ENDIF
           ENDDO
        ENDDO
     ELSE
        DO i=1,nstate
           IF (crge%f(i,1).NE.0._real_8) THEN
              rsum=rsum+crge%f(i,1)*dotp(ncpw%ngw,c0(:,i,1),c0(:,i,1))
           ENDIF
        ENDDO
     ENDIF
     CALL mp_sum(rsum,parai%allgrp)
     scalef=crge%nel/rsum
     CALL dscal(2*nkpt%ngwk*nstate*nkpt%nkpnt,scalef,c0,1)
  ENDIF
1 CONTINUE
  DEALLOCATE(c2,stat=ierr)
  IF(ierr/=0) call stopgm('extrapwf','deallocation problem : c2',&
       __LINE__,__FILE__)

  CALL tihalt('  EXTRAPWF',isub)
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE extrapwf
! ==================================================================
