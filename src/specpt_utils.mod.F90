MODULE specpt_utils
  USE adjmu_utils,                     ONLY: adjmu
  USE canon_utils,                     ONLY: canon,&
                                             give_scr_canon
  USE conv,                            ONLY: nac,&
                                             nbc
  USE coor,                            ONLY: fion,&
                                             tau0,&
                                             taup
  USE corec_utils,                     ONLY: corec
  USE cppt,                            ONLY: nzh
  USE davidson_utils,                  ONLY: davidson,&
                                             give_scr_davidson
  USE ddip,                            ONLY: lenbk
  USE ddipo_utils,                     ONLY: give_scr_ddipo
  USE dftin_utils,                     ONLY: tdm_fun
  USE dist_friesner_utils,             ONLY: dist_friesner,&
                                             give_scr_dist_friesner
  USE dynit_utils,                     ONLY: dynit
  USE efld,                            ONLY: extf,&
                                             textfld
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE fftmain_utils,                   ONLY: fwfftn
  USE fint,                            ONLY: fint1
  USE fnlalloc_utils,                  ONLY: fnlalloc,&
                                             fnldealloc
  USE forcedr_driver,                  ONLY: forcedr
  USE forcedr_utils,                   ONLY: give_scr_forcedr
  USE forcep_utils,                    ONLY: rhoe_psi_size
  USE friesner_utils,                  ONLY: friesner,&
                                             give_scr_friesner
  USE func,                            ONLY: func1
  USE gettrans_utils,                  ONLY: gettrans,&
                                             printtrans
  USE gsortho_utils,                   ONLY: gsortho
  USE initrun_driver,                  ONLY: initrun
  USE initrun_utils,                   ONLY: give_scr_initrun
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE isos,                            ONLY: isos3
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE ksdiag_utils,                    ONLY: ksdiag
  USE linres,                          ONLY: &
       clrwf, ireor, lr01, lr03, lractive, lrf4, lrsym, nlinr, nlinw, nolr, &
       nua, nub, td01, td03, urot, wcent
  USE localize_utils,                  ONLY: localize
  USE lr_diag_utils,                   ONLY: give_scr_lr_diag,&
                                             lr_diag
  USE lr_in_utils,                     ONLY: lr_in,&
                                             tddft_input
  USE lr_tddft_drhoe,                  ONLY: relax_rho
  USE lr_tddft_utils,                  ONLY: setspin
  USE lr_xcpot_utils,                  ONLY: lr_xcpot
  USE machine,                         ONLY: m_walltime
  USE mm_dim_utils,                    ONLY: mm_dim
  USE mm_dimmod,                       ONLY: mm_go_mm,&
                                             mm_go_qm
  USE mm_input,                        ONLY: cgrest_i,&
                                             lqmmm
  USE mm_qmmm_forcedr_utils,           ONLY: mm_qmmm_forcedr
  USE mols,                            ONLY: mstat,&
                                             numol
  USE molstates_utils,                 ONLY: molecular_states
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum,&
                                             mp_sync
  USE nlcc,                            ONLY: corel,&
                                             roct
  USE norm,                            ONLY: cnorm,&
                                             gemax
  USE ortho_utils,                     ONLY: give_scr_ortho
  USE ovlap_utils,                     ONLY: ovlap
  USE parac,                           ONLY: parai,&
                                             paral
  USE pbc_utils,                       ONLY: pbc3
  USE poin,                            ONLY: potr,&
                                             ptau,&
                                             rhoo,&
                                             rtau
  USE prmem_utils,                     ONLY: prmem
  USE pw_hfx_resp,                     ONLY: hfx_resp_finalize,&
                                             hfx_resp_init
  USE pw_hfx_resp_types,               ONLY: hfx_lin2,&
                                             hfx_resp_env
  USE readsr_utils,                    ONLY: xstring
  USE reshaper,                        ONLY: reshape_inplace
  USE rho1pri_utils,                   ONLY: rho1pri
  USE rhoofr_utils,                    ONLY: give_scr_rhoofr,&
                                             rhoofr
  USE rhopri_utils,                    ONLY: give_scr_rhopri,&
                                             rhopri
  USE rnlsm_utils,                     ONLY: give_scr_rnlsm
  USE ropt,                            ONLY: infw,&
                                             iteropt,&
                                             ropt_mod
  USE rv30_utils,                      ONLY: zhrwf
  USE setirec_utils,                   ONLY: read_irec,&
                                             write_irec
  USE soc_types,                       ONLY: c01,&
                                             cs1,&
                                             ct1,&
                                             socvar,&
                                             tausoc
  USE soft,                            ONLY: soft_com
  USE sort_utils,                      ONLY: sort2
  USE spin,                            ONLY: clsd,&
                                             spin_mod,&
                                             tdsp1
  USE stcop_utils,                     ONLY: stcop
  USE store_types,                     ONLY: restart1,&
                                             rout1
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             fpar,&
                                             maxsys,&
                                             ncpw,&
                                             parm
  USE tauf,                            ONLY: itaur,&
                                             tau,&
                                             vtau
  USE td_force_utils,                  ONLY: give_scr_td_force,&
                                             td_force
  USE td_nacvs_utils,                  ONLY: give_scr_vcoupling,&
                                             vcouplings
  USE td_os_utils,                     ONLY: td_os
  USE td_prop_utils,                   ONLY: tdcharge,&
                                             tddipo
  USE testex_utils,                    ONLY: testex
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE tpot,                            ONLY: c0y,&
                                             eigref,&
                                             eigval,&
                                             foccp,&
                                             mstate
  USE updwf_utils,                     ONLY: give_scr_updwf,&
                                             updwf
  USE utils,                           ONLY: nxxfun
  USE wrener_utils,                    ONLY: wreigen,&
                                             wrprint_wfopt
  USE wv30_utils,                      ONLY: zhwwf
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PUBLIC :: specpt

CONTAINS 
  ! ==================================================================
  SUBROUTINE specpt
    ! ==--------------------------------------------------------------==
    ! ==  ELECTRONIC SPECTRA                                          ==
    ! ==--------------------------------------------------------------==
    ! Variables
    CHARACTER(*), PARAMETER                  :: procedureN = 'specpt'

    COMPLEX(real_8), ALLOCATABLE             :: c2(:,:,:), sc0(:,:)
    COMPLEX(real_8), ALLOCATABLE, TARGET     :: c0(:,:,:)
    EXTERNAL                                 :: dasum
    INTEGER                                  :: i, ierr, isub, ms, nc0_1, &
                                                nc0_2, nc0_3, nc2_1, nc2_2, &
                                                nc2_3, no, nsc0_1, nsc0_2, &
                                                nstate, nx
    LOGICAL                                  :: statusdummy
    REAL(real_8)                             :: dasum, du1,du2,du3,du4,du5,du6,du7
    REAL(real_8), ALLOCATABLE                :: fo(:)

    CALL tiset(procedureN,isub)
    ! ==--------------------------------------------------------------==
    IF (tkpts%tkpnt) THEN
       CALL stopgm('SPECPT','K-Points not implemented',& 
            __LINE__,__FILE__)
    ENDIF
    ! ..input for linear response
    CALL lr_in
    CALL tddft_input
    CALL hfx_resp_init
    ! 
    !CSOC[
    IF (cntl%tsoc.AND.(.NOT.socvar%do_sing)) THEN
       td01%ns_tri=td01%ns_sin
       td01%ns_sin=0
    ENDIF
    IF (cntl%tsoc.AND.paral%io_parent) WRITE(6,'(A,2I5)') ' SOC: Number of Sates (sing,tripl): '&
         ,td01%ns_sin,td01%ns_tri
    !CSOC]
    du1=0._real_8
    du2=0._real_8
    du3=0._real_8
    du4=0._real_8
    du5=0._real_8
    du6=0._real_8
    du7=0._real_8
    CALL dynit(du1,du2,du3,du4,du5,du6,du7)
    ! 
    ms=MAX(td01%ns_mix,td01%ns_sin,td01%ns_tri,2)
    no=0
    DO i=1,crge%n
       IF (crge%f(i,1).GT.0.00001_real_8) no=no+1
    ENDDO
    CALL mp_bcast(no,parai%io_source,parai%cp_grp)
    CALL mp_bcast(ms,parai%io_source,parai%cp_grp)

    IF (no+ms.GT.crge%n) THEN
       ALLOCATE(fo(no+ms),STAT=ierr)
       !ALLOCATE(fo(crge%n),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL dcopy(crge%n,crge%f,1,fo,1)
       DEALLOCATE(crge%f,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       crge%n=no+ms
       ALLOCATE(crge%f(crge%n,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(crge%f)!,n)
       IF (cntl%tlsd) THEN
          tdsp1%nupel=NINT(dasum(spin_mod%nsup,fo,1))
          tdsp1%ndoel=NINT(dasum(spin_mod%nsdown,fo(spin_mod%nsup+1),1))
          nx=crge%n-spin_mod%nsup-spin_mod%nsdown
          spin_mod%nsup=spin_mod%nsup+(nx+1)/2
          spin_mod%nsdown=spin_mod%nsdown+nx/2
          DO i=1,tdsp1%nupel
             crge%f(i,1)=1._real_8
          ENDDO
          DO i=spin_mod%nsup+1,spin_mod%nsup+tdsp1%ndoel
             crge%f(i,1)=1._real_8
          ENDDO
       ELSE
          CALL dcopy(no,fo,1,crge%f,1)
       ENDIF
       DEALLOCATE(fo,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    nstate=crge%n
    IF (.NOT.cntl%tlanc) THEN
       IF (cnti%ndavv.EQ.0) THEN
          cnti%ndavv=2*crge%n
       ELSEIF (cnti%ndavv.LT.crge%n+1) THEN
          cnti%ndavv=2*crge%n
       ENDIF
       nc0_1=ncpw%ngw
       nc0_2=cnti%ndavv
       nc0_3=1
       nc2_1=ncpw%ngw
       nc2_2=cnti%ndavv
       nc2_3=1
    ELSE
       nc0_1=ncpw%ngw
       nc0_2=nstate
       nc0_3=1
       nc2_1=ncpw%ngw
       nc2_2=nstate
       nc2_3=1
    ENDIF
    IF (td03%tda.AND.td01%msubta.GT.0.AND..NOT.td03%tdacanon) THEN
       nsc0_1=ncpw%ngw
       nsc0_2=td01%msubta*clsd%nlsd
    ELSE
       nsc0_1=1
       nsc0_2=1
    ENDIF

    ALLOCATE(c0(nc0_1,nc0_2,nc0_3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(c0)
    ALLOCATE(c2(nc2_1,nc2_2,nc2_3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(sc0(nsc0_1,nsc0_2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL mm_dim(mm_go_mm,statusdummy)
    ALLOCATE(fion(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(taup(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(taup)!,SIZE(taup))
    CALL zeroing(fion)!,SIZE(fion))
    CALL mm_dim(mm_go_qm,statusdummy)
    ! 
    CALL spectra(c0,c2,sc0,no,tau0,fion)
    !CSOC[
    IF (cntl%tsoc) THEN
       IF (.NOT.socvar%do_sing) THEN
          CALL dcopy(2*ncpw%ngw*no,c0,1,c01,1)
       ENDIF
    ENDIF
    !CSOC]
    CALL hfx_resp_finalize
    ! 
    DEALLOCATE(c0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(c2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(sc0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(fion,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(taup,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE specpt
  ! ==================================================================
  SUBROUTINE spectra(c0,c2,sc0,no,tau0,fion)
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8), TARGET                  :: c0(:,:,:)
    COMPLEX(real_8)                          :: c2(:,:,:), sc0(:,:)
    INTEGER                                  :: no
    REAL(real_8)                             :: tau0(:,:,:), fion(:,:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'spectra'

    CHARACTER(len=10)                        :: orbital
    CHARACTER(len=30)                        :: tag
    COMPLEX(real_8), ALLOCATABLE             :: c_scr(:), cr(:,:), cscr(:,:), &
                                                cvirt(:,:), cx(:,:), gde(:), &
                                                pme(:), psi(:,:)
    COMPLEX(real_8), ALLOCATABLE, TARGET     :: c1(:,:,:,:)
    INTEGER :: ib, ierr, il_psi_1d, il_psi_2d, il_rhoe_1d, il_rhoe_2d, &
      irec(100), is, ispin, isub, j, k, kprint, lddxc_1d, lddxc_2d, loptib, &
      lscr, mtda, nact, nap, nbp, nconv, ncr, ncscr, ndiag, ngde, nhpsi, &
      norb, npme, nroot, ntp, nus, nvir, nvpp, nxx, transition(3,4,20)
    INTEGER, ALLOCATABLE                     :: order(:)
    LOGICAL                                  :: statusdummy, tadjmu, &
                                                tdasubopt, trefine, update_pot
    REAL(real_8)                             :: detot, dmom(3), en, ena, enb, &
                                                etoto, tcpu, thl(2), time0, &
                                                time1, time2, wk_fake(1)
    REAL(real_8), ALLOCATABLE :: ddxc(:,:), edav(:), eigt(:), eigv(:), &
      rhoe(:,:), rhot(:,:), scr(:), taux(:,:,:), velx(:,:,:), vpp(:)

    CALL tiset(procedureN,isub)
    ! ==================================================================
    cntl%tsymrho=.FALSE.
    lractive=.FALSE.
    !CSOC[
    IF (cntl%tsoc.AND.(.NOT.socvar%do_sing)) THEN
       cntl%tdavi=.TRUE.
    ENDIF
    !CSOC]
    ! 
    IF (lqmmm%qmmm)THEN
       IF (textfld)THEN
          ALLOCATE(extf(fpar%kr1*fpar%kr2s*fpar%kr3s),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          CALL zeroing(extf)!,kr1*kr2s*kr3s)
       ENDIF
    ENDIF
    ! 
    CALL mm_dim(mm_go_mm,statusdummy)
    ALLOCATE(taux(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(velx(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(eigv(crge%n),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(eigv)!,n)
    CALL mm_dim(mm_go_qm,statusdummy)
    ! set variable for Chelovky decomposition
    IF (hfx_resp_env%use_lin_lin) THEN
       hfx_resp_env%lin2_restart=.TRUE.
       norb=crge%n
       IF (cntl%tlsd) CALL stopgm(procedureN,'LSD not implemented',&
            __LINE__,__FILE__)
       ALLOCATE(hfx_lin2%C2(ncpw%ngw,norb),stat=ierr)
       IF (ierr.NE.0) CALL stopgm( procedureN, 'Allocation problem' ,&
            __LINE__,__FILE__)
       ALLOCATE(hfx_lin2%xi(ncpw%ngw,norb),stat=ierr)
       IF (ierr.NE.0) CALL stopgm( procedureN, 'Allocation problem' ,&
            __LINE__,__FILE__)
       ALLOCATE(hfx_lin2%m(norb,norb),stat=ierr)
       IF (ierr.NE.0) CALL stopgm( procedureN, 'Allocation problem' ,&
            __LINE__,__FILE__)
       ALLOCATE(hfx_lin2%l(norb,norb),stat=ierr)
       IF (ierr.NE.0) CALL stopgm( procedureN, 'Allocation problem' ,&
            __LINE__,__FILE__)
       CALL zeroing(hfx_lin2%c2)
       CALL zeroing(hfx_lin2%xi)
       WRITE(*,*) 'allocation in specpt'
    ENDIF
    ! ==--------------------------------------------------------------==
    loptib=lr01%lopti
    lr01%lopti=-1
    IF (lr03%txc_analytic) THEN
       lddxc_1d=fpar%nnr1
       lddxc_2d=(2*clsd%nlsd-1)
    ELSE
       lddxc_1d=1
       lddxc_2d=1
    ENDIF
    ALLOCATE(ddxc(lddxc_1d,lddxc_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    IF (isos3%ps_type.EQ.1) CALL stopgm('SPECTRA','HOCKNEY PS not impl.',& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! SCR ALLOCATION AND PARTITION (SCRATCH ARRAY).
    CALL rhoe_psi_size(il_psi_1d=il_psi_1d,il_psi_2d=il_psi_2d,&
         il_rhoe_1d=il_rhoe_1d,il_rhoe_2d=il_rhoe_2d)
    ntp=1
    IF (td01%ns_tri.GT.0) ntp=2
    ALLOCATE(rhoe(fpar%nnr1,2*il_rhoe_2d),STAT=ierr)!nnr1,2*il_rhoe/nnr1 !vw that is buggus il_rhoe_1d
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(psi(il_psi_1d,il_psi_2d*ntp),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(rhoe)
    CALL zeroing(psi)
    CALL give_scr_spectra(lscr,tag)
    ALLOCATE(scr(lscr),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(scr)
    ! 
    IF (lr03%txc_dd_ana) THEN
       ALLOCATE(rhoo(fpar%nnr1,ntp*il_rhoe_2d),STAT=ierr)!vw that is buggus il_rhoe_1d
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ELSE
       ALLOCATE(rhoo(fpar%nnr1,il_rhoe_2d),STAT=ierr)!vw that is buggus il_rhoe_1d
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    CALL zeroing(rhoo)
    ALLOCATE(potr(fpar%nnr1,il_rhoe_2d),STAT=ierr)!vw that is buggus il_rhoe_1d
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(potr)
    IF (cntl%ttau) THEN
       ALLOCATE(rtau(fpar%nnr1,il_rhoe_2d),STAT=ierr)!vw that is buggus il_rhoe_1d
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(rtau)
       ALLOCATE(ptau(fpar%nnr1,il_rhoe_2d),STAT=ierr)!vw that is buggus il_rhoe_1d
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(ptau)
    ENDIF
    ! 
    CALL fnlalloc(crge%n,.FALSE.,.FALSE.)
    ! 
    IF (td03%treorder) THEN
       ALLOCATE(order(crge%n),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       IF (td03%tdlocal)  THEN
          ALLOCATE(wcent(3,crge%n),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
    ENDIF
    IF (td03%molstat) THEN
       ALLOCATE(wcent(3,crge%n),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! 
    kprint=td01%ioutput
    ! ==--------------------------------------------------------------==
    ! == INITIALISATION                                               ==
    ! ==--------------------------------------------------------------==
    time0=m_walltime()
    ! TIME STEP FUNCTIONS
    ropt_mod%modens=.FALSE.
    ropt_mod%engpri=.FALSE.
    ropt_mod%calste=.FALSE.
    ! Set IREC to restart file.
    CALL read_irec(irec)
    CALL initrun(irec,c0,c2,sc0,rhoe,psi,eigv)
    CALL write_irec(irec)
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent) THEN
       CALL prmem('   SPECTRA')
       WRITE(6,'(1X,64("="))')
       WRITE(6,'(1X,"==",T25,A,T64,"==")')&
            '   REFERENCE POINT'
       WRITE(6,'(1X,64("="))')
    ENDIF
    ropt_mod%convwf=.FALSE.
    ropt_mod%sdiis=.TRUE.
    ropt_mod%spcg=.TRUE.
    ! ..allocate additional memory for optimisers
    IF (cntl%tsde) THEN
       npme = 1
       ngde = 1
       nvpp = 1
    ELSE IF (cntl%diis) THEN
       npme = (ncpw%ngw*crge%n+8)*cnti%mdiis/2
       ngde = ((ncpw%ngw*crge%n+8)*cnti%mdiis)/2
       nvpp = ncpw%ngw
    ELSE IF (cntl%pcg) THEN
       npme = ncpw%ngw*crge%n
       ngde = 1
       nvpp = ncpw%ngw
    ELSEIF (paral%io_parent) THEN
       WRITE(6,*) ' WRONG OPTION FOR WAVEFUNCTION OPTIMIZATION'
       CALL stopgm('SPECTRA',' ',& 
            __LINE__,__FILE__)
    ENDIF
    ALLOCATE(pme(npme),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(gde(ngde),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(vpp(nvpp),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ..virtual orbitals
    nvir=crge%n-no
    ncscr=2*ncpw%ngw*nvir
    ALLOCATE(cvirt(ncpw%ngw,ncscr/ncpw%ngw),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL setspin("S",c0,cvirt,c2)

    IF (lqmmm%qmmm) update_pot=.TRUE.

    etoto =1.e10_real_8
    DO infw=1,cnti%nomore
       time1=m_walltime()
       IF (lqmmm%qmmm .AND. infw.GT.1 ) THEN
          IF (cgrest_i%n_cg.EQ.0)THEN
             update_pot=.FALSE.
          ENDIF
       ENDIF
       ! ..UPDATE THE WAVEFUNCTIONS
       CALL updwf(c0(:,:,1),c2(:,:,1),sc0(:,:),tau0,taux,pme,gde,vpp,eigv,&
            rhoe,psi,no,.FALSE.,update_pot)
       ! ..PRINTOUT THE EVOLUTION OF THE ITERATIVE OPTIMIZATION
       IF (paral%io_parent) THEN
          detot=ener_com%etot-etoto
          IF (infw.EQ.1) detot=0.0_real_8
          time2=m_walltime()
          tcpu=(time2-time1)*0.001_real_8
          CALL wrprint_wfopt(eigv,crge%f,ener_com%amu,no,ener_com%etot,etoto,&
               tcpu,gemax,cnorm,thl,ropt_mod%engpri,infw,infw)
          etoto=ener_com%etot
       ENDIF
       IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)
       IF (ropt_mod%convwf.OR.soft_com%exsoft) GOTO 200
    ENDDO
    CALL setspin("R",c0,cvirt,c2)
    IF (crge%n.GT.no) THEN
       IF (cntl%tlsd) THEN
          CALL gsortho(c0,ncpw%ngw,tdsp1%nupel+1,spin_mod%nsup)
          CALL gsortho(c0(:,spin_mod%nsup+1,1),ncpw%ngw,tdsp1%ndoel+1,spin_mod%nsdown)
       ELSE
          CALL gsortho(c0,ncpw%ngw,no+1,crge%n)
       ENDIF
    ENDIF
    IF (paral%io_parent)WRITE(6,*) ' NO CONVERGENCE IN ',infw,' STEPS '
    CALL mm_dim(mm_go_mm,statusdummy)
    CALL zhwwf(2,irec,c0,c0,crge%n,eigv,tau0,tau0,tau0,infw)
    CALL mm_dim(mm_go_qm,statusdummy)
    ! ==--------------------------------------------------------------==
200 CONTINUE

    ! ..write restart for soft exit.
    IF (soft_com%exsoft) THEN
       CALL mm_dim(mm_go_mm,statusdummy)
       CALL zhwwf(2,irec,c0,c0,crge%n,eigv,tau0,tau0,tau0,infw)
       RETURN
    ENDIF
    CALL setspin("R",c0,cvirt,c2)
    ! ..release memory for ground state optimiser
    DEALLOCATE(vpp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(pme,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(gde,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ..orthogonalize additional states
    IF (crge%n.GT.no) THEN
       IF (cntl%tlsd) THEN
          CALL gsortho(c0,ncpw%ngw,tdsp1%nupel+1,spin_mod%nsup)
          CALL gsortho(c0(:,spin_mod%nsup+1,1),ncpw%ngw,tdsp1%ndoel+1,spin_mod%nsdown)
       ELSE
          CALL gsortho(c0,ncpw%ngw,no+1,crge%n)
       ENDIF
    ENDIF

    CALL mm_dim(mm_go_mm,statusdummy)
    CALL zhwwf(2,irec,c0,c0,crge%n,eigv,tau0,tau0,tau0,iteropt%nfi)
    CALL mm_dim(mm_go_qm,statusdummy)
    ! ..store potential and density
    DO ispin=1,clsd%nlsd
       CALL dcopy(fpar%nnr1,rhoe(1,ispin),1,potr(1,ispin),1)
    ENDDO
    CALL rhoofr(c0(:,:,1),rhoe,psi(:,1),crge%n)
    CALL dcopy(fpar%nnr1*clsd%nlsd,rhoe,1,rhoo,1)
    IF (cntl%ttau) THEN
       CALL dcopy(fpar%nnr1*clsd%nlsd,tau,1,rtau,1)
       CALL dcopy(fpar%nnr1*clsd%nlsd,vtau,1,ptau,1)
    ENDIF
    ! ..NLCC
    IF (corel%tinlc) THEN
       ALLOCATE(roct(fpar%nnr1*clsd%nlsd),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(roct)!,nnr1*clsd%nlsd)
       ALLOCATE(c_scr(ncpw%nhg),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL corec(roct,c_scr,psi)
       DEALLOCATE(c_scr,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    IF (cntl%tlsd) THEN
       ALLOCATE(rhot(fpar%nnr1,clsd%nlsd),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    time2=m_walltime()
    tcpu=(time2-time0)*0.001_real_8
    IF (paral%io_parent) THEN
       WRITE(6,'(/,A,T46,F12.3,A)') ' TIME FOR MINIMUM STRUCTURE :',&
            tcpu,' SECONDS'
    ENDIF
    IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)
    IF (soft_com%exsoft) RETURN
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent) THEN
       CALL prmem('   SPECTRA')
       WRITE(6,'(1X,64("="))')
       WRITE(6,'(1X,"==",T20,A,T64,"==")')&
            'END OF REFERENCE CALCULATION'
       WRITE(6,'(1X,"==",T19,A,T64,"==")')&
            'GENERATE INITIAL GUESS VECTORS'
       WRITE(6,'(1X,64("="),/)')
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (cntl%tlanc) THEN
       fint1%tfral=.TRUE.
       cntl%tlanc=.FALSE.
       ALLOCATE(edav(crge%n),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(edav)!,n)
       ncscr=2*ncpw%ngw*MAX(crge%n,cnti%nkry_max*cnti%nkry_block)
       ALLOCATE(cscr(ncpw%ngw,ncscr/ncpw%ngw),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       IF (func1%mhfx.NE.0) THEN
          ALLOCATE(cx(ncpw%ngw,ncscr/ncpw%ngw),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          CALL dcopy(2*ncpw%ngw*crge%n,c0,1,cx,1)
       ENDIF
       IF (cntl%tlsd) CALL dcopy(fpar%nnr1*clsd%nlsd,potr(1,1),1,rhot(1,1),1)
       IF (paral%io_parent) CALL prmem('   SPECTRA')
       IF (cntl%tlsd) THEN
          IF (cntl%tdmal) THEN
             CALL dist_friesner(spin_mod%nsup,REAL(tdsp1%nupel,kind=real_8),nconv,nhpsi,tdsp1%nupel,&
                  c0,c2,sc0,cscr,cx,crge%f,potr(1,1),psi(:,1),edav,&
                  0,trefine,.TRUE.)
          ELSE
             CALL friesner(spin_mod%nsup,REAL(tdsp1%nupel,kind=real_8),nconv,nhpsi,tdsp1%nupel,&
                  c0,c2,sc0,cscr,cx,crge%f(1:spin_mod%nsup,1),potr(1,1),psi(:,1),edav,&
                  0,trefine,.TRUE.)
          ENDIF
          IF (cntl%tlsd) CALL dcopy(fpar%nnr1*clsd%nlsd,rhot(1,1),1,potr(1,1),1)
          nac=nconv
          ib=spin_mod%nsup+1
          itaur=2
          IF (cntl%tdmal) THEN
             CALL dist_friesner(spin_mod%nsdown,REAL(tdsp1%ndoel,kind=real_8),nconv,nhpsi,tdsp1%ndoel,&
                  c0(:,ib,1),c2,sc0,cscr,cx(:,ib),crge%f(ib,1),potr(1,2),psi(:,1),&
                  edav(ib),0,trefine,.TRUE.)
          ELSE
             CALL friesner(spin_mod%nsdown,REAL(tdsp1%ndoel,kind=real_8),nconv,nhpsi,tdsp1%ndoel,&
                  c0(:,ib,1),c2,sc0,cscr,cx(:,ib),crge%f(ib:ib+spin_mod%nsdown-1,1),potr(1,2),&
                  psi(:,1),edav(ib),0, trefine,.TRUE.)
          ENDIF
          IF (cntl%tlsd) CALL dcopy(fpar%nnr1*clsd%nlsd,rhot(1,1),1,potr(1,1),1)
          itaur=1
          nbc=nconv
          ndiag=crge%n
       ELSE
          ndiag=crge%n
          IF (cntl%tdmal) THEN
             CALL dist_friesner(ndiag,crge%nel,nconv,nhpsi,crge%n,&
                  c0,c2,sc0,cscr,cx,crge%f,potr,psi(:,1),edav,&
                  0,trefine,.TRUE.)
          ELSE
             CALL friesner(ndiag,crge%nel,nconv,nhpsi,crge%n,&
                  c0,c2,sc0,cscr,cx,crge%f(1:ndiag,1),potr,psi(:,1),edav,&
                  0,trefine,.TRUE.)
          ENDIF
          nac=nconv
       ENDIF
       DEALLOCATE(cscr,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       IF (func1%mhfx.NE.0) DEALLOCATE(cx,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ELSEIF (cntl%tdavi) THEN
       !ivano
       IF (hfx_resp_env%use_lin_lin) THEN                                                               
          hfx_resp_env%switch_lin2=.TRUE.                                                                
          hfx_resp_env%use_lin_lin=.FALSE.                                                               
       ENDIF
       ALLOCATE(edav(2*cnti%ndavv),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(edav)!,2*cnti%ndavv)
       ncscr=2*ncpw%ngw*cnti%ndavv
       ALLOCATE(cscr(ncpw%ngw,ncscr/ncpw%ngw),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ncr=2*ncpw%ngw*cnti%ndavv
       ALLOCATE(cr(ncpw%ngw,ncr/ncpw%ngw),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(cx(ncpw%ngw,ncscr/ncpw%ngw),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(vpp(ncpw%ngw),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL ksdiag(vpp)
       IF (paral%io_parent) CALL prmem('   SPECTRA')
       IF (cntl%tlsd) THEN
          CALL dcopy(2*ncpw%ngw*spin_mod%nsup,c0,1,cx,1)
          CALL dcopy(fpar%nnr1*clsd%nlsd,potr(1,1),1,rhot(1,1),1)
          CALL davidson(spin_mod%nsup,cnti%ndavv,cx,c2,cr,sc0,cscr,potr(:,1),&
               c0(:,1:spin_mod%nsup,1),crge%f(1:spin_mod%nsup,1),tdsp1%nupel,psi(:,1),edav,vpp)
          CALL dcopy(fpar%nnr1*clsd%nlsd,rhot(1,1),1,potr(1,1),1)
          CALL dcopy(2*ncpw%ngw*spin_mod%nsup,cx,1,c0,1)
          ib = spin_mod%nsup+1
          itaur=2
          CALL dcopy(2*ncpw%ngw*spin_mod%nsdown,c0(1,ib,1),1,cx,1)
          CALL davidson(spin_mod%nsdown,cnti%ndavv,cx,c2,cr,sc0,cscr,potr(:,2),&
               c0(:,ib:ib+spin_mod%nsdown-1,1),crge%f(ib:ib+spin_mod%nsdown-1,1),tdsp1%ndoel,&
               psi(:,1),edav(cnti%ndavv+1),vpp)
          itaur=1
          CALL dcopy(2*ncpw%ngw*spin_mod%nsdown,cx,1,c0(1,ib,1),1)
          CALL dcopy(fpar%nnr1*clsd%nlsd,rhot(1,1),1,potr(1,1),1)
          DO is=1,spin_mod%nsdown
             edav(spin_mod%nsup+is)=edav(cnti%ndavv+is)
          ENDDO
       ELSE
          ndiag=crge%n
          CALL dcopy(2*ncpw%ngw*ndiag,c0,1,cx,1)
          CALL davidson(ndiag,cnti%ndavv,cx,c2,cr,sc0,cscr,potr(:,1),&
               c0(:,1:ndiag,1),crge%f(1:no,1),no,psi(:,1),edav,vpp)
          CALL dcopy(2*ncpw%ngw*ndiag,cx,1,c0,1)
       ENDIF
       DEALLOCATE(cx,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(vpp,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(cr,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(cscr,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       cntl%tdavi=.FALSE.
    ENDIF
    IF (cntl%tlsd) THEN
       DEALLOCATE(rhot,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF

    CALL mp_sync(parai%cp_grp)
    CALL mm_dim(mm_go_mm,statusdummy)
    CALL zhwwf(2,irec,c0,c0,crge%n,edav,tau0,tau0,tau0,iteropt%nfi)
    CALL mm_dim(mm_go_qm,statusdummy)
    IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)
    IF (soft_com%exsoft) RETURN
    IF (paral%parent) THEN
       tadjmu=.TRUE.
       wk_fake(:)=-1.0_real_8!vw need to set this guy, it is referenced in adjmu
       CALL adjmu(crge%n,1,crge%nel,fint1%betael,edav,wk_fake,cntl%tlsd,ener_com%amu,tadjmu)
       CALL wreigen(edav,crge%f,ener_com%amu,crge%n)
       IF (paral%io_parent) THEN
          WRITE(6,'(1X,64("="))')
          WRITE(6,'(1X,"==",T20,A,T64,"==")')&
               'END OF STATE INITIALIZATION'
          WRITE(6,'(1X,64("="),/)')
       ENDIF
    ENDIF
    ! ..store virtual states
    IF (cntl%tlsd) THEN
       nap=spin_mod%nsup-tdsp1%nupel
       nbp=spin_mod%nsdown-tdsp1%ndoel
       CALL dcopy(2*ncpw%ngw*nap,c0(1,tdsp1%nupel+1,1),1,cvirt(1,1),1)
       CALL dcopy(2*ncpw%ngw*nbp,c0(1,spin_mod%nsup+tdsp1%ndoel+1,1),1,cvirt(1,nap+1),1)
    ELSE
       CALL dcopy(2*ncpw%ngw*nvir,c0(1,no+1,1),1,cvirt(1,1),1)
    ENDIF
    ! ..make some choices on the form of the occupied orbitals
    IF (td03%treorder) THEN
       DO is=1,no
          order(is)=is
       ENDDO
       DO is=1,td01%nreor
          order(no-is+1)=ireor(is)
       ENDDO
       IF (paral%io_parent) THEN
          ! Check if the ordering is reasonable
          DO is=1,no
             DO j=is+1,no
                IF (order(is).EQ.order(j)) CALL stopgm("SPECTRA",&
                     "Illegal reordering of states",& 
                     __LINE__,__FILE__)
             ENDDO
          ENDDO
       ENDIF
    ENDIF
    IF (td03%tdlocal) THEN
       orbital="GENERAL"
       time1=m_walltime()
       CALL setspin("S",c0,cvirt,c2)
       lenbk=nxxfun(no)
       nxx=MAX(2*lenbk*parai%nproc,2*ncpw%ngw*no)
       ALLOCATE(cx(ncpw%ngw,ncscr/ncpw%ngw),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(cr(ncpw%ngw,ncr/ncpw%ngw),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL localize(tau0,c0,cx,cr,no)
       DEALLOCATE(cr,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(cx,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       time2=m_walltime()
       tcpu=(time2-time1)*0.001_real_8
       IF (paral%io_parent) THEN
          WRITE(6,'(/,A,T46,F12.3,A)')&
               ' TIME FOR WAVEFUNCTION LOCALIZATION :',tcpu,' SECONDS'
       ENDIF
       DEALLOCATE(eigv,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(eigv(no*no),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(eigv)!,no*no)
       IF (td03%treorder.AND.td01%zatom.NE.0) THEN
          CALL worder(no,tau0,order)
       ENDIF
       CALL setspin("R",c0,cvirt,c2)
    ELSE
       orbital="CANON"
    ENDIF
    ! ==--------------------------------------------------------------==
    tdasubopt=.FALSE.
    IF (td03%tda.AND.td01%msubta.GT.0) THEN
       nact=td01%msubta
       IF (.NOT.td03%tdacanon) tdasubopt=.TRUE.
    ELSE
       IF (cntl%tlsd) THEN
          nact=tdsp1%nupel+tdsp1%ndoel
       ELSE
          nact=no
       ENDIF
    ENDIF
    ! ..calculate KS-Matrix
    IF (lrf4%td_method.EQ.1) THEN
       IF (INDEX(orbital,"CANON").NE.0) THEN
          ! transform to canonical orbitals and get eigenvalues
          CALL setspin("S",c0,cvirt,c2)
          IF (lqmmm%qmmm)THEN
             CALL mm_dim(mm_go_qm,statusdummy)
             CALL mm_qmmm_forcedr(c0,c2,sc0,rhoe,psi,tau0,fion,&
                  eigv,no,1,.FALSE.,.FALSE.,.TRUE.)
          ELSE
             CALL forcedr(c0(:,:,1),c2(:,:,1),sc0,rhoe,psi,tau0,fion,eigv,&
                  no,1,.FALSE.,.FALSE.)
          ENDIF
          CALL canon(c0,c2,crge%f,no,eigv)
          IF (td03%treorder) THEN
             ! This seems not to be right for LSD -> empty states?
             DO is=1,no
                CALL dcopy(2*ncpw%ngw,c0(1,order(is),1),1,c2(1,is,1),1)
                eigv(is)=edav(order(is))
             ENDDO
             CALL dcopy(2*ncpw%ngw*no,c2,1,c0,1)
             CALL dcopy(no,eigv,1,edav,1)
          ENDIF
          CALL dcopy(crge%n,edav,1,eigv,1)
          IF (td03%molstat) THEN
             ncscr=2*ncpw%ngw*no
             ALLOCATE(c0y(ncscr,1,1),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             CALL dcopy(ncscr,c0,1,c0y,1)
             CALL molecular_states(c0y,eigv,no,tau0)
          ENDIF
          CALL tdm_fun(lrf4%td_functional,1)
          IF (cntl%tpotential) THEN
             ALLOCATE(eigval(no),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             ALLOCATE(eigref(no),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             ALLOCATE(foccp(no),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             CALL dcopy(no,crge%f,1,foccp,1)
             mstate = no
             IF (td03%molstat) THEN
                IF (cntl%tlsd) THEN
                   CALL stopgm("SPECTRA","cntl%tlsd MOLSTAT",& 
                        __LINE__,__FILE__)
                ELSE
                   CALL dcopy(no,eigv,1,eigval,1)
                   is=0
                   DO j=1,numol
                      en=eigval(is+mstat(j))
                      DO k=is+1,is+mstat(j)
                         eigref(k)=en
                      ENDDO
                      is=is+mstat(j)
                   ENDDO
                ENDIF
             ELSE
                IF (cntl%tlsd) THEN
                   CALL dcopy(tdsp1%nupel,eigv,1,eigval,1)
                   CALL dcopy(tdsp1%ndoel,eigv(tdsp1%nupel+nap+1),1,eigval(tdsp1%nupel+1),1)
                ELSE
                   CALL dcopy(no,eigv,1,eigval,1)
                ENDIF
                c0y => c0
                IF (cntl%tlsd) THEN
                   DO is=1,tdsp1%nupel+nap
                      IF (foccp(is).GT.1.e-4_real_8) ena=eigv(is)
                   ENDDO
                   DO is=1,tdsp1%nupel+nap
                      eigref(is)=ena
                   ENDDO
                   DO is=tdsp1%nupel+nap+1,mstate
                      IF (foccp(is).GT.1.e-4_real_8) enb=eigv(is)
                   ENDDO
                   DO is=tdsp1%nupel+nap+1,mstate
                      eigref(is)=enb
                   ENDDO
                ELSE
                   DO is=1,mstate
                      IF (foccp(is).GT.1.e-4_real_8) en=eigv(is)
                   ENDDO
                   DO is=1,mstate
                      eigref(is)=en
                   ENDDO
                ENDIF
             ENDIF
          ENDIF
          orbital="GENERAL"
          DEALLOCATE(eigv,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(eigv(no*no),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          CALL zeroing(eigv)!,no*no)
          DEALLOCATE(scr,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
          CALL give_scr_spectra(lscr,tag)
          ALLOCATE(scr(lscr),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          IF (lqmmm%qmmm)THEN
             CALL mm_dim(mm_go_qm,statusdummy)
             CALL mm_qmmm_forcedr(c0,c2,sc0,rhoe,psi,tau0,fion,&
                  eigv,no,1,.FALSE.,.FALSE.,.FALSE.)
          ELSE
             CALL forcedr(c0(:,:,1),c2(:,:,1),sc0,rhoe,psi,tau0,fion,eigv,&
                  no,1,.FALSE.,.TRUE.)
          ENDIF
          CALL dcopy(fpar%nnr1*clsd%nlsd,rhoe,1,potr,1)
          DO is=1,no
             CALL dscal(2*ncpw%ngw,1._real_8/crge%f(is,1),c2(1,is,1),1)
          ENDDO
          CALL ovlap(no,eigv,c0(:,:,1),c2(:,:,1))
          CALL mp_sum(eigv,no*no,parai%allgrp)
          CALL setspin("R",c0,cvirt,c2)
          IF (td03%tda.AND.td01%msubta.GT.0) THEN
             ! Reorder KS Matrix
             CALL stopgm("SPECTRA","TDMETHOD-MSUBTA",& 
                  __LINE__,__FILE__)
          ENDIF
          IF (td03%molstat) THEN
             DEALLOCATE(c0y,STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
                  __LINE__,__FILE__)
             c0y => c0
          ENDIF
          IF (cntl%tpotential) THEN
             DEALLOCATE(eigval,STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             DEALLOCATE(eigref,STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             DEALLOCATE(foccp,STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
          END IF
       ELSE
          CALL stopgm("SPECTRA","TDMETHOD-GENERAL",& 
               __LINE__,__FILE__)
       ENDIF
    ELSEIF (INDEX(orbital,"CANON").NE.0) THEN
       ! transform to canonical orbitals and get eigenvalues
       CALL setspin("S",c0,cvirt,c2)
       IF (lqmmm%qmmm)THEN
          CALL mm_dim(mm_go_qm,statusdummy)
          CALL mm_qmmm_forcedr(c0,c2,sc0,rhoe,psi,tau0,fion,&
               eigv,no,1,.FALSE.,.FALSE.,.TRUE.)
       ELSE
          CALL forcedr(c0(:,:,1),c2(:,:,1),sc0,rhoe,psi,tau0,fion,eigv,&
               no,1,.FALSE.,.TRUE.)
       ENDIF
       CALL canon(c0,c2,crge%f,no,eigv)
       IF (td03%treorder) THEN
          ! This seems not to be right for LSD -> empty states?
          DO is=1,no
             CALL dcopy(2*ncpw%ngw,c0(1,order(is),1),1,c2(1,is,1),1)
             eigv(is)=edav(order(is))
          ENDDO
          CALL dcopy(2*ncpw%ngw*no,c2,1,c0,1)
          CALL dcopy(no,eigv,1,edav,1)
       ENDIF
       CALL dcopy(crge%n,edav,1,eigv,1)
       CALL setspin("R",c0,cvirt,c2)
    ELSE
       IF (td03%treorder) THEN
          DO is=1,no
             CALL dcopy(2*ncpw%ngw,c0(1,order(is),1),1,c2(1,is,1),1)
          ENDDO
          CALL dcopy(2*ncpw%ngw*no,c2,1,c0,1)
       ENDIF
       CALL setspin("S",c0,cvirt,c2)
       IF (lqmmm%qmmm)THEN
          CALL mm_dim(mm_go_qm,statusdummy)
          CALL mm_qmmm_forcedr(c0,c2,sc0,rhoe,psi,tau0,fion,&
               eigv,no,1,.FALSE.,.FALSE.,.TRUE.)
       ELSE
          CALL forcedr(c0(:,:,1),c2(:,:,1),sc0,rhoe,psi,tau0,fion,eigv,&
               no,1,.FALSE.,.TRUE.)
       ENDIF
       DO is=1,no
          CALL dscal(2*ncpw%ngw,1._real_8/crge%f(is,1),c2(1,is,1),1)
       ENDDO
       CALL ovlap(no,eigv,c0(:,:,1),c2(:,:,1))
       CALL mp_sum(eigv,no*no,parai%allgrp)
       IF (td03%tda.AND.td01%msubta.GT.0) THEN
          ! Reorder KS Matrix
          DO is=1,td01%msubta
             j=(no-td01%msubta+1) + (no-td01%msubta+is-1)*no
             k=(is-1)*td01%msubta+1
             CALL dcopy(td01%msubta,eigv(j),1,eigv(k),1)
          ENDDO
       ENDIF
       CALL setspin("R",c0,cvirt,c2)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == END INITIALIZATION                                           ==
    ! ==--------------------------------------------------------------==
    IF (cntl%tlsd) THEN
       crge%n=tdsp1%nupel+tdsp1%ndoel
    ELSE
       crge%n=no
    ENDIF
    nroot=MAX(td01%ns_sin,td01%ns_tri,td01%ns_mix)
    ALLOCATE(eigt(nroot),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(eigt)!,nroot)
    IF (paral%io_parent) CALL prmem('   SPECTRA')
    lractive=.TRUE.
    ! ..LR-cntl%tddft EQUATIONS AND TAMM-DANCOFF APPROXIMATION
    IF (td03%tda.AND.td01%ldiag.EQ.1) THEN
       mtda=1
    ELSE
       mtda=2
    ENDIF
    IF (paral%io_parent) THEN
       IF (cntl%tlsd) WRITE(6,'(1X,"==",T28,A,T64,"==")') 'MIXED STATES'
       IF (td01%ns_sin.NE.0) WRITE(6,'(1X,"==",T26,A,T64,"==")')&
            'SINGLET STATES'
       IF (td01%ns_tri.NE.0) WRITE(6,'(1X,"==",T26,A,T64,"==")')&
            'TRIPLET STATES'
    ENDIF
    ! ..calculate and store analytic dmu/dn
    CALL lr_xcpot(ddxc,rhoo,.TRUE.)
    ! 
    IF (td03%tda.AND.td01%msubta.GT.0) THEN
       nolr=clsd%nlsd*td01%msubta
    ELSE
       IF (cntl%tlsd) THEN
          nolr=tdsp1%nupel+tdsp1%ndoel
       ELSE
          nolr=no
       ENDIF
    ENDIF
    nus=clsd%nlsd*nact*crge%n*nroot+1
    ALLOCATE(urot(nus,1,1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(urot)!,nus)
    ! 
    ALLOCATE(c1(ncpw%ngw,nolr,nroot,mtda),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    nua=crge%n
    nub=nact
    nlinw=nroot
    nlinr=0
    clrwf => c1
    IF (restart1%rlr) THEN
       IF (cntl%tlsd) THEN
          crge%n=nolr+nap+nbp
       ELSE
          crge%n=no+nvir
       ENDIF
       CALL read_irec(irec)
       CALL zhrwf(1,irec,c0,c2,crge%n,edav,tau0,velx,taux,iteropt%nfi)
       CALL write_irec(irec)
       IF (cntl%tlsd) THEN
          crge%n=tdsp1%nupel+tdsp1%ndoel
       ELSE
          crge%n=no
       ENDIF
    ENDIF

    ! >>>>QDF
    ! vw this should be moved inside ZHRWF... but so many pointer carbages...
    CALL mp_bcast(taux,SIZE(taux),parai%io_source,parai%cp_grp) ! causes memory corruption while inside zhrwf
    CALL mp_bcast(c1,nroot*ncpw%ngw*nolr*mtda,&
         parai%cp_inter_io_source,parai%cp_inter_grp)
    ! <<<

    IF (cntl%tlsd) THEN
       lrsym=0
       CALL setspin("S",c0,cvirt,c2)
       CALL stcop(td01%ns_mix,c0,c2,crge%n,eigv,edav,cvirt,nap,nbp,c1,nact,&
            ddxc,psi,rhoe,orbital)
       DO is=tdsp1%nupel+1,crge%n
          eigv(is)=edav(nap+is)
       ENDDO
    ELSEIF (td01%ns_sin.NE.0) THEN
       lrsym=1
       CALL stcop(td01%ns_sin,c0,c2,crge%n,eigv,edav,cvirt,nvir,0,c1,nact,&
            ddxc,psi,rhoe,orbital)
    ELSEIF (td01%ns_tri.NE.0) THEN
       lrsym=3
       CALL stcop(td01%ns_tri,c0,c2,crge%n,eigv,edav,cvirt,nvir,0,c1,nact,&
            ddxc,psi,rhoe,orbital)
    ENDIF
    IF (paral%io_parent) CALL prmem('   SPECTRA')
    IF (cntl%ttau) THEN
       CALL dcopy(fpar%nnr1*clsd%nlsd,rtau,1,tau,1)
       CALL dcopy(fpar%nnr1*clsd%nlsd,ptau,1,vtau,1)
    ENDIF
    CALL lr_diag(c0(:,:,1),c1(:,:,:,1),c2(:,:,1),sc0,eigv,ddxc,rhoe,psi,&
         eigt,crge%n,nroot,orbital,kprint)
    CALL zeroing(transition)!,3*4*20)
    CALL gettrans(transition,c1,cvirt,&
         nact,nvir,nap,nbp,nroot)
    CALL printtrans(transition,eigt,nact,nvir,nap,nbp,nroot)
    CALL td_os(nact,c0,nroot,eigt,c1)
    IF (rout1%rhoout) CALL rho1pri(c0,c1,tau0,crge%n,nroot,"S")
    IF (cntl%tlsd) CALL setspin("R",c0,cvirt,c2)
    !CSOC[ 
    IF (cntl%tsoc) THEN
       IF (socvar%do_sing) THEN
          CALL dcopy(nroot*2*ncpw%ngw*nolr*mtda,c1,1,cs1,1)
       ELSE
          CALL dcopy(nroot*2*ncpw%ngw*nolr*mtda,c1,1,ct1,1)
          CALL dcopy(3*maxsys%nax*maxsys%nsx,tau0,1,tausoc,1)
       ENDIF
    ENDIF
    !CSOC]
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent) CALL prmem('   SPECTRA')
    IF (cntl%tlsd) THEN
       crge%n=nolr+nap+nbp
    ELSE
       crge%n=no+nvir
    ENDIF
    IF (cntl%tnacvs) THEN
       CALL give_scr_vcoupling(SIZE(psi),tag)
       IF (SIZE(psi).GT.lscr) THEN
          DEALLOCATE(scr,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(scr(SIZE(psi)),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          lscr=SIZE(scr)
       ENDIF
       CALL vcouplings(c0(:,:,1),c1,c2(:,:,1),sc0,rhoe,eigv,eigt,&
            crge%n,no,nvir,nact,nroot)
    ENDIF
    CALL mm_dim(mm_go_mm,statusdummy)
    CALL zhwwf(2,irec,c0,c0,crge%n,eigv,tau0,tau0,tau0,iteropt%nfi)
    CALL mm_dim(mm_go_qm,statusdummy)
    IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)
    IF (soft_com%exsoft) RETURN
    IF (lrf4%td_method.EQ.1) THEN
       CALL tdm_fun(lrf4%ks_functional,1)
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (td03%tprop) THEN
       IF (lrf4%td_method.EQ.1) THEN
          CALL stopgm("TD_METHOD","NOT IMPLEMENTED",& 
               __LINE__,__FILE__)
       ENDIF
       IF (paral%io_parent) THEN
          WRITE(6,*)
          WRITE(6,'(1X,64("="))')
          WRITE(6,'(1X,"==",T13,A,T64,"==")')&
               '      EXCITED STATE PROPERTY CALCULATION    '
          WRITE(6,'(1X,64("="))')
       ENDIF
       IF (cntl%tlsd) CALL setspin("S",c0,cvirt,c2)
       IF (cntl%tlsd) THEN
          crge%n=tdsp1%nupel+tdsp1%ndoel
       ELSE
          crge%n=no
       ENDIF
       lr01%lopti=loptib-10
       ! Print ground state density
       IF (paral%io_parent) THEN
          WRITE(6,'(" >>>>>> ",A,T63,I3)') 'Process ground state'
       ENDIF
       CALL rhopri(c0,tau0,rhoe,psi(:,1),crge%n,1)
       CALL tdcharge(rhoo(:,1),tau0,psi(:,1))
       CALL tddipo(dmom,rhoo(:,1),tau0,psi(:,1))
       IF (td01%npstate.EQ.0) THEN
          CALL dcopy(3*maxsys%nax*maxsys%nsx,fion,1,taux,1)
          DO is=1,nroot
             td01%fstate=is
             IF (paral%io_parent) THEN
                WRITE(6,'(" >>>>>> ",A,T63,I3)') 'Process state ',td01%fstate
             ENDIF
             CALL dcopy(3*maxsys%nax*maxsys%nsx,taux,1,fion,1)
             CALL td_force(c0(:,:,1),c1,c2(:,:,1),sc0,eigv,rhoe,psi,eigt,&
                  crge%n,nroot,tau0,fion,orbital,1)
             ! calculate relaxed density
             CALL relax_rho(c0,c1,c2,sc0,rhoo,rhoe,psi(:,1),crge%n,nroot)
             CALL relax_pri(rhoe,tau0,scr,psi(:,1),td01%fstate)
             CALL tdcharge(rhoe(:,1),tau0,psi(:,1))
             CALL tddipo(dmom,rhoe(:,1),tau0,psi(:,1))
          ENDDO
       ELSE
          td01%fstate=MIN(td01%npstate,nroot)
          IF (td01%fstate.NE.td01%npstate) CALL stopgm("SPECPT",&
               "INVALID STATE FOR PROPERTY CALCULATION",& 
               __LINE__,__FILE__)
          IF (paral%io_parent) THEN
             WRITE(6,'(" >>>>>> ",A,T63,I3)') 'Process state ',td01%fstate
          ENDIF
          CALL td_force(c0(:,:,1),c1,c2(:,:,1),sc0,eigv,rhoe,psi,eigt,&
               crge%n,nroot,tau0,fion,orbital,1)
          ! calculate relaxed density
          CALL relax_rho(c0,c1,c2,sc0,rhoo,rhoe,psi(:,1),crge%n,nroot)
          CALL relax_pri(rhoe,tau0,scr,psi(:,1),td01%fstate)
          CALL tdcharge(rhoe(:,1),tau0,psi(:,1))
          CALL tddipo(dmom,rhoe(:,1),tau0,psi(:,1))
       ENDIF
       IF (cntl%tlsd) CALL setspin("R",c0,cvirt,c2)
    ENDIF
    ! ==--------------------------------------------------------------==
    DEALLOCATE(ddxc,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(scr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(c1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(urot,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    CALL fnldealloc(.FALSE.,.FALSE.)
    DEALLOCATE(edav,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(taux,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(velx,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rhoe,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(psi,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(eigt,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(cvirt,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (td03%treorder) DEALLOCATE(order,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (td03%treorder.AND.td03%tdlocal) DEALLOCATE(wcent,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (td03%molstat) THEN
       DEALLOCATE(wcent,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    IF (corel%tinlc) THEN
       DEALLOCATE(roct,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    IF (cntl%ttau) THEN
       DEALLOCATE(rtau,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(ptau,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    DEALLOCATE(potr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rhoo,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (lqmmm%qmmm.AND.textfld) THEN
       DEALLOCATE(extf,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    IF (hfx_resp_env%use_lin_lin) THEN
       DEALLOCATE(hfx_lin2%c2, stat=ierr)
       IF (ierr.NE.0) CALL stopgm( procedureN, 'Deallocation problem' ,&
            __LINE__,__FILE__)
       DEALLOCATE(hfx_lin2%m, stat=ierr)
       IF (ierr.NE.0) CALL stopgm( procedureN, 'Deallocation problem' ,&
            __LINE__,__FILE__)
       DEALLOCATE(hfx_lin2%l, stat=ierr)
       IF (ierr.NE.0) CALL stopgm( procedureN, 'Deallocation problem' ,&
            __LINE__,__FILE__)
       DEALLOCATE(hfx_lin2%xi, stat=ierr)
       IF (ierr.NE.0) CALL stopgm( procedureN, 'Deallocation problem' ,&
            __LINE__,__FILE__)
       WRITE(*,*) 'de allocation in specpt'
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE spectra
  ! ==================================================================
  SUBROUTINE give_scr_spectra(lspectra,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lspectra
    CHARACTER(len=30)                        :: tag

    INTEGER :: l_diag, lcanon, lforces, linitrun, llocal, llr_diag, lortho, &
      lrhoofr, lrhopri, lrnlsm, ltd_force, lupdwf, nstate

    nstate=crge%n
    CALL give_scr_initrun(linitrun,tag)
    CALL give_scr_updwf(lupdwf,tag,nstate,.FALSE.)
    CALL give_scr_rhoofr(lrhoofr,tag)
    CALL give_scr_forcedr(lforces,tag,nstate,.FALSE.,.FALSE.)
    CALL give_scr_canon(lcanon,tag,nstate)
    CALL give_scr_rnlsm(lrnlsm,tag,nstate,.FALSE.)
    CALL give_scr_lr_diag(llr_diag,tag)
    nstate=crge%n+MAX(td01%ns_tri,td01%ns_sin,td01%ns_mix)
    IF (cntl%tlanc) THEN
       IF (cntl%tdmal) THEN
          CALL give_scr_dist_friesner(l_diag,tag,nstate)
       ELSE
          CALL give_scr_friesner(l_diag,tag,nstate)
       ENDIF
    ELSE
       CALL give_scr_davidson(l_diag,tag,nstate,cnti%ndavv)
    ENDIF
    CALL give_scr_ortho(lortho,tag,nstate)
    IF (td03%tdlocal) THEN
       CALL give_scr_ddipo(llocal,tag)
       llocal=MAX(llocal,2*ncpw%nhg)
    ELSE
       llocal=0
    ENDIF
    IF (corel%tinlc) THEN
       llocal=MAX(llocal,ncpw%nhg)
    ENDIF
    lrhopri=0
    ltd_force=0
    IF (td03%tprop) THEN
       CALL give_scr_rhopri(lrhopri,tag,nstate)
       lrhopri=MAX(lrhopri,2*ncpw%nhg)
       CALL give_scr_td_force(ltd_force,nstate,tag)
    ENDIF
    ! 
    lspectra=MAX(linitrun,lupdwf,lrhoofr,lforces,lcanon,lrnlsm,&
         llr_diag,l_diag,lortho,llocal,lrhopri,ltd_force)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_spectra
  ! ==================================================================
  SUBROUTINE worder(nstate,tau0,order)
    ! Variables
    INTEGER                                  :: nstate
    REAL(real_8)                             :: tau0(:,:,:)
    INTEGER                                  :: order(*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'worder'

    INTEGER                                  :: ia, iat, ierr, iref, is, iz
    REAL(real_8)                             :: dd, dx, dy, dz, dx_, dy_, dz_
    REAL(real_8), ALLOCATABLE                :: distwr(:), ref(:,:)

! ==--------------------------------------------------------------==

    IF (paral%io_parent) THEN
       ALLOCATE(distwr(nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(ref(3,td01%zatom),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       iat=0
       iref=0
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             iat=iat+1
             DO iz=1,td01%zatom
                IF (ireor(iz).EQ.iat) THEN
                   iref=iref+1
                   ref(1,iref)=tau0(1,ia,is)
                   ref(2,iref)=tau0(2,ia,is)
                   ref(3,iref)=tau0(3,ia,is)
                   GOTO 2000
                ENDIF
             ENDDO
2000         CONTINUE
          ENDDO
       ENDDO
       DO is=1,nstate
          distwr(is)=1.e30_real_8
          DO ia=1,td01%zatom
             dx_=wcent(1,is)-ref(1,ia)
             dy_=wcent(2,is)-ref(2,ia)
             dz_=wcent(3,is)-ref(3,ia)
             CALL pbc3(dx_,dy_,dz_,dx,dy,dz,1,parm%apbc,parm%ibrav)
             dd=(dx*dx+dy*dy+dz*dz)
             distwr(is)=MIN(dd,distwr(is))
          ENDDO
       ENDDO
       CALL dscal(nstate,-1._real_8,distwr,1)
       CALL sort2 ( distwr, nstate, order )
       DEALLOCATE(ref,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(distwr,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    CALL mp_bcast(order,nstate,parai%io_source,parai%cp_grp)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE worder
  ! ==================================================================
  SUBROUTINE relax_pri(rhoe,tau0,rg_in,psi,fstate)
    REAL(real_8)                             :: rhoe(:,:), tau0(*)
    REAL(real_8), TARGET                     :: rg_in(:)
    COMPLEX(real_8)                          :: psi(:)
    INTEGER                                  :: fstate

    CHARACTER(len=40)                        :: cflbod, cipnum, filen
    COMPLEX(real_8), DIMENSION(:), POINTER   :: rg
    INTEGER                                  :: i1, i2, ig, ir, n1, n2, newdim

    newdim=SIZE(rg_in)/2
    CALL reshape_inplace(rg_in, (/newdim/), rg)
    cflbod='RELAXDEN_'
    IF (paral%io_parent)&
         WRITE(cipnum,'(I4)') fstate
    CALL xstring(cflbod,n1,n2)
    CALL xstring(cipnum,i1,i2)
    filen=cflbod(n1:n2)//cipnum(i1:i2)
    !$omp parallel do private (IR)
#ifdef __SR8000
    !poption parallel
#endif
    DO ir=1,fpar%nnr1
       psi(ir)=CMPLX(rhoe(ir,1),0._real_8,kind=real_8)
    ENDDO
    CALL  fwfftn(psi,.FALSE.,parai%allgrp)
    DO ig=1,ncpw%nhg
       rg(ig) = psi(nzh(ig))
    ENDDO
    CALL densto(rg,tau0,filen)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE relax_pri
  ! ==================================================================
END MODULE specpt_utils
