#include "cpmd_global.h"

MODULE rwfopt_utils
  USE andp
  USE cdft_utils,                      ONLY: cdft_adarrays,&
                                             cdft_reinit,&
                                             cdft_w,&
                                             init_cdft,&
                                             vupdate,&
                                             write_w
  USE cdftmod,                         ONLY: &
       cdftci, cdftcom, cdfthda, cdftlog, cdftpi, sccomm, wd, wdiff, wgaussl
  USE cnst
  USE coor
  USE cplngs_utils,                    ONLY: cplngs,&
                                             give_scr_cplngs,&
                                             prcplngs
  USE cplngsmod
  USE ddipo_utils,                     ONLY: give_scr_ddipo
  USE detdof_utils,                    ONLY: detdof
  USE dynit_utils,                     ONLY: dynit
  USE efld
  USE ehrenfest_utils,                 ONLY: ehrenfest
  USE elct
  USE enbandpri_utils,                 ONLY: enbandpri
  USE ener
  USE epr_efg_utils,                   ONLY: epr_efg
  USE error_handling,                  ONLY: stopgm
  USE fft_maxfft,                      ONLY: maxfft
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE finalp_utils,                    ONLY: finalp
  USE forcep_utils,                    ONLY: rhoe_psi_size
  USE forces_diag_utils,               ONLY: updiag
  USE geofile_utils,                   ONLY: geofile
  USE gsize_utils,                     ONLY: gsize
  USE hfxmod
  USE hubbardu,                        ONLY: hubbu
  USE hubbardu_utils,                  ONLY: hubbardUcorrection
  USE initrun_driver,                  ONLY: initrun
  USE initrun_utils,                   ONLY: give_scr_initrun
  USE ions
  USE isos
  USE k_updwf_utils,                   ONLY: give_scr_kupdwf,&
                                             k_updwf
  USE kinds,                           ONLY: real_8
  USE kpts
  USE linres
  USE localize_utils,                  ONLY: localize,&
                                             localize2
  USE locpot
  USE lscal
  USE machine,                         ONLY: m_walltime
  USE mm_dim_utils,                    ONLY: mm_dim
  USE mm_dimmod
  USE mm_input
  USE mm_qmmm_forcedr_utils,           ONLY: mm_qmmm_forcedr
  USE mp_interface,                    ONLY: mp_sum
  USE norm
  USE parac,                           ONLY: parai,&
                                             paral
  USE poin
  USE prmem_utils,                     ONLY: prmem
  USE rhoofr_c_utils,                  ONLY: rhoofr_c
  USE rhoofr_utils,                    ONLY: rhoofr
  USE rhopri_utils,                    ONLY: give_scr_rhopri,&
                                             rhopri
  USE rinitwf_driver,                  ONLY: rinitwf
  USE ropt
  USE rscpot_utils,                    ONLY: give_scr_rscpot,&
                                             rscpot
  USE rswfmod
  USE rv30_utils,                      ONLY: zhrwf
  USE setirec_utils,                   ONLY: read_irec,&
                                             write_irec
  USE spin
  USE store_types
  USE syscomb_utils,                   ONLY: syscomb,&
                                             write_ksham
  USE system,                          ONLY: &
       cnti, cntl, cntr, fpar, locpot2, maxsys, nacc, ncpw, nkpt, parm, spar
  USE td_input
  USE td_utils,                        ONLY: getnorm_k,&
                                             load_ex_states,&
                                             tmprd_prop
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE totstr_utils,                    ONLY: totstr
  USE transme_utils,                   ONLY: transme
  USE updrho_utils,                    ONLY: give_scr_updrho,&
                                             updrho
  USE updwf_utils,                     ONLY: give_scr_updwf,&
                                             updwf
  USE utils,                           ONLY: icopy
  USE vdwcmod
  USE wann
  USE wrener_utils,                    ONLY: wreigen,&
                                             wrener,&
                                             wrprint_wfopt
  USE wv30_utils,                      ONLY: zhwwf
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rwfopt
  PUBLIC :: give_scr_rwfopt
  !public :: lseprof

CONTAINS

  ! ==================================================================
  SUBROUTINE rwfopt(c0,c2,sc0,pme,gde,vpp,eigv)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8), INTENT(inout)           :: c0(:,:,:), c2(:,:), sc0(:,:,:)
    COMPLEX(real_8)                          :: pme(:), gde(:)
    REAL(real_8)                             :: vpp(:), &
                                                eigv(crge%n,nkpt%nkpts)

    CHARACTER(*), PARAMETER                  :: procedureN = 'rwfopt'
    CHARACTER(len=11), PARAMETER             :: filename = 'LOCPOT.cube'

    CHARACTER(len=30)                        :: filen, tag
    COMPLEX(real_8), ALLOCATABLE             :: c0ini(:,:,:), csave(:), &
                                                rsave(:), rsave2(:),c2u(:,:)

#if defined(_HAS_CUDA)
    COMPLEX(real_8), DIMENSION(:,:), POINTER __CONTIGUOUS :: psi
#else
    COMPLEX(real_8), ALLOCATABLE             :: psi(:,:)
#endif

    COMPLEX(real_8), ALLOCATABLE, TARGET     :: cf_4d(:,:,:,:)
    INTEGER :: cdfti, i, iaa, ierr, il_psi_1d, il_psi_2d, il_rhoe_1d, &
         il_rhoe_2d, int_update_pot, irec(100), iss, isub, ixx, iyy, izz, ldim, &
         lscr, nhpsi, nnx, nrxyz1s, nrxyz2s, &
         recompute_dipole_matrices_every_tmp, w_opt_tmp
    LOGICAL                                  :: convtest, error, oldstatus, &
         statusdummy, tfor, twscc, &
         update_pot
    REAL(real_8) :: chrg_(2), delx1, delx2, delx3, dely1, dely2, dely3, &
         delz1, delz2, delz3, ekin1, ekin2, ekincp, ekinh1, ekinh2, esave, &
         etot0, n2save(2), nsave(2), tcpu, temp1, temp2, thl(2), time1, time2, &
         vsave(2)
    REAL(real_8), ALLOCATABLE                :: norms(:), rhoe(:,:), &
         rhoini(:,:), scr(:), &
         tscr(:,:,:), vpotx3s(:,:,:), &
         wdsave(:)
    REAL(real_8), POINTER                    :: vpotx3(:,:,:)

    CALL tiset(procedureN,isub)
    time1 =m_walltime()
    ! ==--------------------------------------------------------------==
    ! SCR creation
    CALL rhoe_psi_size(il_psi_1d=il_psi_1d,il_psi_2d=il_psi_2d,&
         il_rhoe_1d=il_rhoe_1d,il_rhoe_2d=il_rhoe_2d)
    ALLOCATE(rhoe(fpar%nnr1,il_rhoe_2d),STAT=ierr) !vw doenst work in parallel with il_rhoe_1d
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(rhoe)
    ! EHR[
    IF (cntl%tpspec.OR.cntl%tpdist) THEN
       il_psi_1d=MAX(il_psi_1d,2*maxfft)
       ALLOCATE(norms(nkpt%ngwk*crge%n),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ! CALL SETPOSOP
    ENDIF
    ! EHR]
    !#if defined(_HAS_CUDA)
    !    CALL cuda_alloc_host(psi,[il_psi_1d,il_psi_2d])
    !#else
    ALLOCATE(psi(il_psi_1d,il_psi_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', __LINE__,__FILE__)
    !#endif
    CALL zeroing(psi)
    CALL give_scr_rwfopt(lscr,tag)
    ALLOCATE(scr(lscr),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    tfor = (cprint%iprint(iprint_force).EQ.1)
    IF (lqmmm%qmmm) THEN
       IF (textfld) THEN
          ALLOCATE(extf(fpar%kr1*fpar%kr2s*fpar%kr3s),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          CALL zeroing(extf)!,kr1*kr2s*kr3s)
       ENDIF
    ENDIF

    CALL mm_dim(mm_go_mm,oldstatus)
    ! CB: TAUP,FION,VELP handled by BS_WFO
    IF (.NOT.cntl%bsymm) THEN
       ALLOCATE(taup(3,maxsys%nax,maxsys%nsx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(fion(3,maxsys%nax,maxsys%nsx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(taup)!,3*maxsys%nax*maxsys%nsx)
       CALL zeroing(fion)!,3*maxsys%nax*maxsys%nsx)
    ENDIF
    ALLOCATE(tscr(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    !
    ! Non-adiabatic couplings
    IF (tcpl) THEN
       cntl%tsymrho=.FALSE.
       ALLOCATE(cf_4d(nkpt%ngwk,crge%n,nkpt%nkpnt,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(c0ini(nkpt%ngwk,crge%n,nkpt%nkpnt),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       IF (cntl%tdiag.OR.tcpllr)  THEN
          ALLOCATE(rhoini(fpar%nnr1,clsd%nlsd),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
       IF (tcpllr) THEN
          clrwf => cf_4d

          IF (isos3%ps_type.EQ.1) THEN
             CALL stopgm('COUPLINGS','HOCKNEY PS not impl.',&
                  __LINE__,__FILE__)
          ENDIF
       ENDIF
       IF (tkpts%tkpnt) THEN
          ALLOCATE(cplion(3,maxsys%nax,maxsys%nsx,nsurf),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(cplcfg(3,maxsys%nax,maxsys%nsx),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ELSE
          ALLOCATE(cplion(3,maxsys%nax,maxsys%nsx,nsurf),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(cplcfg(3,maxsys%nax,maxsys%nsx),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
       CALL zeroing(cf_4d)!,SIZE(cf_4d))
       CALL zeroing(c0ini)!,SIZE(c0ini))
    ENDIF
    !
    CALL mm_dim(mm_go_qm,statusdummy)
    nacc = 7
    IF (cntl%tdiag.OR.cntl%tdiagopt) THEN
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
       rhoo => rin0
    ELSEIF (tkpts%tkpnt) THEN
       nnx=fpar%nnr1*clsd%nlsd
       ALLOCATE(rin0(fpar%nnr1,clsd%nlsd),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(rmix(fpar%nnr1,clsd%nlsd),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    !FIXME error checking
    IF(cntl%thubb) ALLOCATE(c2u(ncpw%ngw,crge%n))
    iteropt%nfi = 0
    ener_com%ecnstr=0.0_real_8
    ener_com%erestr=0.0_real_8
    ! TIME STEP FUNCTIONS
    CALL dynit(ekincp,ekin1,ekin2,temp1,temp2,ekinh1,ekinh2)
    ropt_mod%convwf=.FALSE.
    ropt_mod%sdiis=.TRUE.
    ropt_mod%spcg=.TRUE.
    ropt_mod%modens=.FALSE.
    ropt_mod%engpri=.FALSE.
    ropt_mod%calste=.FALSE.
    ! ==--------------------------------------------------------------==
    ! == INITIALIZATION                                               ==
    ! ==--------------------------------------------------------------==
    IF (cntl%cdft)CALL init_cdft()
    IF (cdftlog%thda)THEN
       ldim  = (maxstates+1)/2 * lwdim
       ALLOCATE(csave(ncpw%ngw*crge%n+4),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(rsave(ldim),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(rsave2(ldim),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ELSE
       cdfthda%hdafirst=.FALSE.
    ENDIF
    IF (cntl%tsyscomb)THEN
       CALL syscomb(c0,c2,sccomm%n_s0,sccomm%n_s1,sccomm%n_s0up,sccomm%n_s1up)
       GOTO 150
    ENDIF
    CALL initrun(irec,c0,c2,sc0,rhoe,psi,eigv)
    IF (cntl%tksham)THEN
       CALL write_ksham(c0,c2,sc0,rhoe,psi,eigv)
       GOTO 150
    ENDIF
    ! McB  this is not necessary in standard runs
    ! McB  (maybe for cntl%cdft, but it is WRONG! for QM/MM (cf. setsys.F)!)
    ! IF(TCLUST.AND.TCENT.AND..NOT.TCLAS)CALL QM_CENTRE
    CALL write_irec(irec)
    IF (paral%parent) THEN
       ! McB    DETDOF needs to be called in full system dimensions
       CALL mm_dim(mm_go_mm,statusdummy)
       CALL detdof(tau0,tscr)
       CALL mm_dim(mm_go_qm,statusdummy)
    ENDIF
    IF (cntl%cdft)THEN
       twscc=.FALSE.
       IF (hfxc3%twscr) THEN
          hfxc3%twscr=.FALSE.
          twscc=.TRUE.
       ENDIF
       CALL rhoofr(c0(:,:,1),rhoe,psi(:,1),crge%n)
       IF (twscc) hfxc3%twscr=.TRUE.
       CALL cdft_w(rhoe,tau0,chrg_)
       IF (cdftlog%tauto)THEN
          cdftcom%cdft_nc=chrg_(2)-chrg_(1)
          IF (paral%io_parent)THEN
             WRITE(6,'(1X,64("-"),/)')
             WRITE(6,'(1X,A)') 'CDFT NC AUTO'
             WRITE(6,'(2X,A,T54,F10.3)') "Setting NC to",cdftcom%cdft_nc
          ENDIF
          IF (cdftlog%thda)THEN
             IF (.NOT.wgaussl%thdas)THEN
                cdftcom%nother=-cdftcom%cdft_nc
             ELSE
                cdftcom%nother=cdftcom%cdft_nc
             ENDIF
             IF (paral%io_parent)THEN
                WRITE(6,'(2X,A,T54,F10.3)') "and other state to",cdftcom%nother
             ENDIF
          ENDIF
          IF (paral%io_parent)THEN
             WRITE(6,'(2X,A,T54,E10.2E1)') "DONOR: ",chrg_(2)
             WRITE(6,'(2X,A,T54,E10.2E1,/)') "ACCEPTOR: ",chrg_(1)
             WRITE(6,'(1X,64("-"),/)')
          ENDIF
       ENDIF
       IF (cntl%cdft_weight)THEN
          IF (cdftlog%thda)THEN
             IF (wgaussl%thdas.OR.wgaussl%thdawm)THEN
                CALL write_w(wdiff,"HDA-STATE1")
             ELSE
                CALL write_w(wdiff,"HDA")
             ENDIF
          ELSE
             CALL write_w(wdiff,"WFOPT")
          ENDIF
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == END INITIALIZATION                                           ==
    ! ==--------------------------------------------------------------==
    etot0=0._real_8
    ener_com%ebogo=0._real_8
    IF (paral%io_parent) THEN
       time2 =m_walltime()
       tcpu = (time2 - time1)*0.001_real_8
       WRITE(6,'(A,T50,F8.2,A8)')&
            ' TIME FOR WAVEFUNCTION INITIALIZATION:',tcpu,' SECONDS'
       CALL prmem('    RWFOPT')
    ENDIF
    DO ! cntl%cdft HDA
       ! ME... IF statement introduced by M. Eichinger for interface mode
       ! EHR[
       IF ((cntl%tinter.AND.restart1%restart.AND.restart1%rwf) .OR.&
            cntl%tpspec.OR.cntl%tpdist ) THEN
          ! initialize external field
          IF (lqmmm%qmmm) THEN
             update_pot=.TRUE.
             int_update_pot=1
             CALL mm_qmmm_forcedr(c0,c2,sc0,rhoe,psi,tau0,fion,&
                  eigv,crge%n,0,.FALSE.,update_pot,.FALSE.)
          ENDIF
          ! do not optimize
          GOTO 100
          ! EHR]
       ELSE
          update_pot=.TRUE.
          int_update_pot=1
          ! initialize external field
          IF (cntl%tdiag .AND. lqmmm%qmmm)&
               CALL mm_qmmm_forcedr(c0,c2,sc0,rhoe,psi,tau0,fion,&
               eigv,crge%n,0,.FALSE.,update_pot,.FALSE.)
          ! ==--------------------------------------------------------------==
          ! ==      THE BASIC LOOP FOR WAVEFUNCTION OPTIMIZATION (cntl%cdft)       ==
          ! ==--------------------------------------------------------------==
          IF (cntl%cdft) THEN
             DO cdfti=1,cdftci%cdft_end
                DO infi=1,cnti%nomore_iter
                   time1=m_walltime()
                   iteropt%nfi=iteropt%nfi+1
                   IF (infi.GT.1) THEN
                      int_update_pot=0
                      update_pot=.FALSE.
                   ENDIF
                   CALL mm_dim(mm_go_qm,statusdummy)
                   IF (cntl%tdiag) THEN
                      CALL updrho(c0,c2,gde,sc0,pme,vpp,tau0,fion,eigv,&
                           rhoe,psi,&
                           crge%n,.FALSE.,.TRUE.,.FALSE.,infi,thl,nhpsi)
                      IF(cntl%thubb)THEN
                         hubbu%tpom=.true.
                         IF(hubbu%debug)  THEN
                             IF (paral%io_parent) write(6,*) procedureN, "| starting hubbardUcorrection"
                         ENDIF
                         CALL hubbardUcorrection(c0(:,:,1),c2u,tau0,fion,crge%n,psi(:,1),.FALSE.,0)
                      ENDIF
                   ELSE
                      IF (tkpts%tkpnt) THEN
                         CALL k_updwf(c0,c2,sc0,tau0,fion,pme,gde,vpp,eigv,&
                              rhoe,psi,crge%n,.FALSE.,infi)
                      ELSE
                         IF (cntl%tdipd.AND.(wan05%loc_relocalize_in_scf.OR.infi==1)) &
                              CALL localize2(tau0,c0,c2,sc0,crge%n)
                         ! UPDATE THE WAVEFUNCTIONS
                         CALL updwf(c0(:,:,1),c2,sc0(:,:,1),tau0,fion,pme,gde,vpp,eigv,&
                              rhoe,psi,crge%n,.FALSE.,update_pot)
                         IF (cntl%tdiagopt.AND.(cnti%nperid.NE.0.OR.cnti%nrperi.NE.0)) THEN
                            CALL updiag(c0,c2,gde,sc0,pme,vpp,tau0,fion,eigv,&
                                 rhoe,psi,crge%n,iteropt%nfi,.FALSE.,.TRUE.,&
                                 time1,etot0)
                         ENDIF
                      ENDIF
                   ENDIF
                   CALL mm_dim(mm_go_mm,statusdummy)
                   IF (paral%io_parent) THEN
                      ropt_mod%engpri=MOD(infi-1,cprint%iprint_step).EQ.0
                   ELSE
                      ropt_mod%engpri=.FALSE.
                   ENDIF

                   ! PRINTOUT THE EVOLUTION OF THE ITERATIVE OPTIMIZATION
                   IF (paral%io_parent) THEN
                      time2=m_walltime()
                      tcpu=(time2-time1)*0.001_real_8
                      ! VENER=CDFT_V*(CDFT_NC+VGRAD)
                      CALL wrprint_wfopt(eigv,crge%f,ener_com%amu,crge%n,ener_com%etot,etot0,tcpu,&
                           gemax,cnorm,thl,.FALSE.,bsnfi+iteropt%nfi,2)
                   ENDIF
                   ! PERIODICALLY WRITE THE RESTART FILE
                   ! CB:   In BS case RESTART is written in BS_WFO
                   IF (.NOT.cntl%bsymm.AND.(MOD(iteropt%nfi,store1%istore).EQ.0.OR.&
                        infi.EQ.cnti%nomore_iter.OR.ropt_mod%convwf).AND.MOD(cdfti,store1%isctore).EQ.0)&
                        THEN
                      CALL mm_dim(mm_go_mm,statusdummy)
                      CALL zhwwf(2,irec,c0,c2,crge%n,eigv,tau0,velp,taup,iteropt%nfi)
                      CALL mm_dim(mm_revert,statusdummy)
                   ENDIF
                   IF (ropt_mod%convwf) THEN
                      GOTO 101
                   ENDIF
                ENDDO
101             CONTINUE
                ropt_mod%convwf=.FALSE.
                CALL vupdate(c0,psi,rhoe,convtest)
                ropt_mod%sdiis=.TRUE.
                ropt_mod%spcg=.TRUE.
                IF (paral%io_parent)CALL wrener
                IF (convtest) THEN
                   IF (.NOT.cntl%bsymm.AND.MOD(cdfti,store1%isctore).NE.0)THEN
                      CALL mm_dim(mm_go_mm,statusdummy)
                      CALL zhwwf(2,irec,c0,c2,crge%n,eigv,tau0,velp,taup,iteropt%nfi)
                      CALL mm_dim(mm_revert,statusdummy)
                   ENDIF
                   ropt_mod%convwf=.TRUE.
                   GOTO 100
                ENDIF
                IF (cdftlog%reswf) THEN
                   IF (paral%io_parent) WRITE(6,*) "Re-initialising Wavefunction"
                   CALL rinitwf(c0,c2,sc0,crge%n,tau0,taup,rhoe,psi)
                ENDIF
             ENDDO
          ELSE
             ! ==--------------------------------------------------------------==
             ! ==     NORMAL LOOP                                              ==
             ! ==--------------------------------------------------------------==
             DO infi=1,cnti%nomore_iter
                time1=m_walltime()
                iteropt%nfi=iteropt%nfi+1
                IF (infi.GT.1) THEN
                   int_update_pot=0
                   update_pot=.FALSE.
                ENDIF
                CALL mm_dim(mm_go_qm,statusdummy)
                IF (cntl%tdiag) THEN
                   CALL updrho(c0,c2,gde,sc0,pme,vpp,tau0,fion,eigv,&
                        rhoe,psi,&
                        crge%n,.FALSE.,.TRUE.,.FALSE.,infi,thl,nhpsi)
                ELSE
                   IF (tkpts%tkpnt) THEN
                      CALL k_updwf(c0,c2,sc0,tau0,fion,pme,gde,vpp,eigv,&
                           rhoe,psi,crge%n,.FALSE.,infi)
                   ELSE
                      IF (cntl%tdipd.AND.(wan05%loc_relocalize_in_scf.OR.infi==1))&
                           CALL localize2(tau0,c0,c2,sc0,crge%n)
                      ! UPDATE THE WAVEFUNCTIONS
                      CALL updwf(c0(:,:,1),c2,sc0(:,:,1),tau0,fion,pme,gde,vpp,eigv,&
                           rhoe,psi,crge%n,.FALSE.,update_pot)
                      IF (cntl%tdiagopt.AND.(cnti%nperid.NE.0.OR.cnti%nrperi.NE.0)) THEN
                         CALL updiag(c0,c2,gde,sc0,pme,vpp,tau0,fion,eigv,&
                              rhoe,psi,crge%n,iteropt%nfi,.FALSE.,.TRUE.,&
                              time1,etot0)
                      ENDIF
                   ENDIF
                ENDIF
                CALL mm_dim(mm_go_mm,statusdummy)
                IF (paral%io_parent) THEN
                   ropt_mod%engpri=MOD(infi-1,cprint%iprint_step).EQ.0
                ELSE
                   ropt_mod%engpri=.FALSE.
                ENDIF
                ! PRINTOUT THE EVOLUTION OF THE ITERATIVE OPTIMIZATION
                IF (paral%io_parent) THEN
                   time2=m_walltime()
                   tcpu=(time2-time1)*0.001_real_8
                   CALL wrprint_wfopt(eigv,crge%f,ener_com%amu,crge%n,ener_com%etot,etot0,tcpu,&
                        gemax,cnorm,thl,ropt_mod%convwf,bsnfi+iteropt%nfi,infi)
                ENDIF
                ! PERIODICALLY WRITE THE RESTART FILE
                ! CB:   In BS case RESTART is written in BS_WFO
!!! CUR
                IF(ropt_mod%convwf.OR.infi.EQ.cnti%nomore_iter) THEN
                   IF(cntl%tdipd.OR.hfxc3%twscr) THEN
                      w_opt_tmp=wanni%w_opt
                      wanni%w_opt=2
                      recompute_dipole_matrices_every_tmp=wan05%loc_recompute_dipole_matrices_every
                      wan05%loc_recompute_dipole_matrices_every=1
                      CALL localize2(tau0,c0,c2,sc0,crge%n)
                      wanni%w_opt=w_opt_tmp
                      wan05%loc_recompute_dipole_matrices_every=recompute_dipole_matrices_every_tmp
                   ENDIF
                ENDIF
!!! CUR
                IF (.NOT.cntl%bsymm.AND.&
                     (MOD(infi,store1%istore).EQ.0.OR.infi.EQ.cnti%nomore_iter.OR.ropt_mod%convwf))THEN
                   CALL mm_dim(mm_go_mm,statusdummy)
                   CALL zhwwf(2,irec,c0,c2,crge%n,eigv,tau0,velp,taup,iteropt%nfi)
                   CALL mm_dim(mm_revert,statusdummy)
                ENDIF
                IF (ropt_mod%convwf) THEN
                   GOTO 100
                ENDIF
             ENDDO
          ENDIF
       ENDIF
       ! ==--------------------------------------------------------------==
       ! ==     END OF MAIN LOOP                                         ==
       ! ==--------------------------------------------------------------==
       IF (paral%io_parent.AND.gemax.GT.cntr%tolog.AND.&
            .NOT.cntl%ksener.AND..NOT.tkpts%tknoswap.AND..NOT.cntl%bsymm) THEN
          WRITE(6,'(1X,64("!"))')
          WRITE(6,'(" !!",A,T64,"!!")')&
               ' RWFOPT| THE MAXIMUM NUMBER OF STEP IS REACHED'
          WRITE(6,'(" !!",A,F10.6,A,T64,"!!")')&
               '         BUT NO CONVERGENCE (DRHOMAX=',gemax,')'
          WRITE(6,'(1X,64("!"))')
       ENDIF
       IF (tkpts%tknoswap) THEN
          ! No calculation of ionic forces.
          ropt_mod%convwf=.FALSE.
       ENDIF
100    CONTINUE
       ! kk-mb === print local potential (start) ===
       IF (locpot2%tlpot) THEN
          IF (lg_vpotx3a) THEN
             vpotx3 => vpotx3a
          ELSE IF (lg_vpotx3b) THEN
             vpotx3 => vpotx3b
          ENDIF
          nrxyz1s=spar%nr1s*spar%nr2s*spar%nr3s
          ALLOCATE(vpotx3s(spar%nr1s,spar%nr2s,spar%nr3s),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          CALL zeroing(vpotx3s)!,nrxyz1s)
          !msglenn=8*nrxyz1s
          CALL mp_sum(vpotx3,vpotx3s,nrxyz1s,parai%allgrp)
          !msglenn=1*8/irat
          CALL mp_sum(nrxyz2,nrxyz2s,parai%allgrp)
          IF (paral%io_parent) THEN
             IF (nrxyz2s .NE. nrxyz1s) THEN
                WRITE(6,*) ' *** WRONG NON-ZERO DATA POINTSs ***'
                WRITE(6,*) '   DATA POINTS:',nrxyz2S
                WRITE(6,*) '   NR1S*NR2S*NR3S:',nrxyz1S
                CALL stopgm('RWFOPT','PROBLEM IN LOCPOT',&
                     __LINE__,__FILE__)
             ENDIF
             WRITE(6,'(A)')&
                  ' LOCAL POTENTIAL = Delta_Vloc(r)+Vxc(r)+VHt(r)+Vcore(r)'
             WRITE(6,'(A,A11)')&
                  ' WRITTEN TO FILE ',FILENAME
             CALL fileopen(12,filename,64,error)
             IF (error) CALL stopgm('RWFOPT','PROBLEM IN LOCPOT',&
                  __LINE__,__FILE__)
             WRITE(12,'(A)')&
                  'LOCAL POTENTIAL = Delta_Vloc(r)+Vxc(r)+VHt(r)+Vcore(r)'
             WRITE(12,'(A)') 'created by CPMD'
             WRITE(12,'(I5)') ions1%nat
             delx1=parm%a1(1)/REAL(spar%nr1s,kind=real_8)
             delx2=0._real_8
             delx3=0._real_8
             dely1=0._real_8
             dely2=parm%a2(2)/REAL(spar%nr2s,kind=real_8)
             dely3=0._real_8
             delz1=0._real_8
             delz2=0._real_8
             delz3=parm%a3(3)/REAL(spar%nr3s,kind=real_8)
             WRITE(12,500) spar%nr1s,delx1,delx2,delx3
             WRITE(12,500) spar%nr2s,dely1,dely2,dely3
             WRITE(12,500) spar%nr3s,delz1,delz2,delz3
             DO iss=1,ions1%nsp
                DO iaa=1,ions0%na(iss)
                   WRITE(12,'(I5,4F12.6)') ions0%iatyp(iss),REAL(ions0%iatyp(iss),kind=real_8),&
                        TAU0(1,IAA,ISS),&
                        TAU0(2,IAA,ISS),&
                        TAU0(3,IAA,ISS)
                ENDDO
             ENDDO
             DO ixx=1,spar%nr1s
                DO iyy=1,spar%nr2s
                   WRITE(12,'(6E13.5)') (vpotx3s(ixx,iyy,izz),izz=1,spar%nr3s)
                ENDDO
             ENDDO
500          FORMAT(i5,3f12.6)
             CALL fileclose(12)
          ENDIF
       ENDIF
       ! kk-mb === print local potential (end)  ===
       ! ==--------------------------------------------------------------==
       ! Diagonalize optimized wavefunction if required
       IF (cntl%tdiagopt.AND..NOT.tcpl) THEN
          IF (paral%io_parent)&
               WRITE(6,'(1X,A)') 'DIAGONALIZING OPTIMIZED WAVEFUNCTION'
          CALL updiag(c0,c2,gde,sc0,pme,vpp,tau0,fion,eigv,&
               rhoe,psi,crge%n,iteropt%nfi,tfor,.TRUE.,&
               time1,etot0)
          IF (paral%io_parent) THEN
             time2=m_walltime()
             tcpu=(time2-time1)*0.001_real_8
             CALL wrprint_wfopt(eigv,crge%f,ener_com%amu,crge%n,ener_com%etot,etot0,tcpu,&
                  gemax,cnorm,thl,ropt_mod%convwf,iteropt%nfi,infi)
             IF (cprint%iprint(iprint_eigen).EQ.1)&
                  CALL wreigen(eigv,crge%f,ener_com%amu,crge%n)
          ENDIF
       ENDIF
       CALL mm_dim(mm_go_qm,statusdummy)
       IF (rout1%rhoout.AND.(.NOT.cntl%bsymm.OR.ropt_mod%convwf)) THEN
          CALL rhopri(c0,tau0,rhoe,psi(:,1),crge%n,nkpt%nkpnt)
       ENDIF
       ! EHR[
       ! ==--------------------------------------------------------------==
       ! ==     COMPUTE SPECTRA USING PROPAGATION OF PERTURBED WFs
       ! ==--------------------------------------------------------------==
       IF (cntl%tpdist) THEN
          CALL mm_dim(mm_go_qm,statusdummy)
          textfld=.TRUE.
          ALLOCATE(extf(fpar%kr1*fpar%kr2s*fpar%kr3s),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          CALL zeroing(extf)!,kr1*kr2s*kr3s)
          IF (.NOT.(cntl%tsde)) THEN
             IF (paral%io_parent) THEN
                WRITE(6,*) 'STOP IN rwfopt.F.'
                WRITE(6,*) 'USE STEEPEST DESCENT FOR WF OPTIMIZATION.'
             ENDIF
             STOP
          ENDIF
          IF (lqmmm%qmmm) THEN
             IF (paral%io_parent) THEN
                WRITE(6,*) 'STOP IN rwfopt.F.'
                WRITE(6,*) 'TPDIST AND QMMM NOT COMPATIBLE'
             ENDIF
             STOP
          ENDIF
       ENDIF
       IF (cntl%tpspec.OR.cntl%tpdist) THEN
          CALL mm_dim(mm_go_qm,statusdummy)
          IF (cntl%ruku) THEN
             IF (paral%io_parent) WRITE(6,*)&
                  'NO RUNGE-KUTTA WITH SPECTRA CALCULATION'
             cntl%ruku=.FALSE.
          ENDIF
          IF (cntl%start_real) THEN
             ! LOAD WAVEFUNCTIONS FORM FILE DENSITY.-NM
             IF (td_prop%read_wf.EQ.1) THEN
                CALL load_ex_states(c0)
                ! CALL GETNORM_K(C0,N,NORMS)
                ! DO I=1,N
                ! CALL DSCAL(2*NGWK,1._real_8/SQRT(NORMS(I)),C0(1,I,1),1)
                ! ENDDO
             ELSEIF (td_prop%read_wf.EQ.2) THEN
                IF (paral%io_parent)&
                     WRITE(6,*) 'WAVEFUNCTIONS FROM RESTART FILE'
                ! TAKE PROP. WFS FROM THE OPT. WFS
             ENDIF
          ELSE
             ! READ FROM FORT.120
             filen="wavefunctions"
             CALL tmprd_prop(c0,crge%n,i,filen)
             IF (paral%io_parent)&
                  WRITE(6,*) 'WAVEFUNCTIONS READ FROM THE FILE WAVEFUNCTIONS'
          ENDIF
          CALL rhoofr_c(c0,rhoe,psi(:,1),crge%n)
          CALL getnorm_k(c0,crge%n,norms)
          IF (paral%io_parent) THEN
             WRITE(6,*) " NORM OF THE STARTING KS STATES"
             WRITE(6,'(2X,63("-"))')
             WRITE(6,'(2X,7(1X,F8.3))') (norms(i),i=1,crge%n)
             WRITE(6,'(2X,63("-"))')
          ENDIF
          CALL ehrenfest(c0,c2,rhoe,psi(:,1),sc0,eigv)
          ropt_mod%convwf=.TRUE.
          GOTO 150
       ENDIF
       ! EHR]
       ! ==--------------------------------------------------------------==
       ! cntl%cdft calculate the transition matrix
       IF (cdfthda%hdafirst)THEN
          IF (cntl%tdipd) CALL localize(tau0,c0,c2,sc0,crge%n)
          CALL dcopy(2*ncpw%ngw*crge%n,c0,1,csave,1)
          CALL zcopy(ldim,rswf,1,rsave,1)
          esave=ener_com%etot
          vsave=cdftcom%cdft_v
          nsave(1)=cdftcom%cdft_nc+cdftcom%vgrad(1)
          nsave(2)=cdftcom%cdft_ns+cdftcom%vgrad(2)
          cdfthda%hdafirst=.FALSE.
          ropt_mod%convwf=.FALSE.
          cdftcom%cdft_nc=cdftcom%nother
          IF (wgaussl%thdas)THEN
             CALL icopy(cdftpi%naccr,cdftpi%cdft_a,1,cdftpi%spd,1)
             CALL icopy(cdftpi%naccr,cdftpi%cdft_d,1,cdftpi%cdft_a,1)
             CALL icopy(cdftpi%naccr,cdftpi%spd,1,cdftpi%cdft_d,1)
             CALL cdft_adarrays
          ELSE IF (wgaussl%thdawm)THEN
             ALLOCATE(wdsave(fpar%nnr1),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             CALL dcopy(fpar%nnr1,wdiff,1,wdsave,1)
             CALL icopy(cdftpi%naccr,cdftpi%cdft_a(cdftpi%naccr+1),1,cdftpi%cdft_a,1)
             CALL icopy(cdftpi%naccr,cdftpi%spa(cdftpi%naccr+1),1,cdftpi%spa,1)
             CALL icopy(cdftpi%naccr,cdftpi%sna(cdftpi%naccr+1),1,cdftpi%sna,1)
             CALL icopy(cdftpi%ndon,cdftpi%cdft_d(cdftpi%ndon+1),1,cdftpi%cdft_d,1)
             CALL icopy(cdftpi%ndon,cdftpi%spd(cdftpi%ndon+1),1,cdftpi%spd,1)
             CALL icopy(cdftpi%ndon,cdftpi%snd(cdftpi%ndon+1),1,cdftpi%snd,1)
          ELSE
             cdftcom%cdft_v=-cdftcom%cdft_v
          ENDIF

          CALL cdft_reinit()
          ropt_mod%sdiis=.TRUE.
          ropt_mod%spcg=.TRUE.
          IF (paral%io_parent)THEN
             WRITE(6,*) "HDA: Switching NC to ",cdftcom%nother,&
                  "to calculate second state"
          ENDIF
          IF (.NOT.cdftlog%rcdft)THEN
             CALL rinitwf(c0,c2,sc0,crge%n,tau0,taup,rhoe,psi)
          ELSE
             IF (cdfthda%hdaresb.AND..NOT.restart1%rnon)THEN
                restart1%rco=.FALSE.
                CALL read_irec(irec)
                CALL zhrwf(1,irec,c0,c2,crge%n,eigv,tau0,velp,taup,iteropt%nfi)
             ELSE
                CALL rinitwf(c0,c2,sc0,crge%n,tau0,taup,rhoe,psi)
             ENDIF
          ENDIF
          CALL cdft_w(rhoe,tau0,chrg_)
          IF (cntl%cdft_weight.AND.cdftlog%thda.AND.(wgaussl%thdas.OR.wgaussl%thdawm))THEN
             CALL write_w(wdiff,"HDA-STATE2")
          ENDIF
       ELSE
          EXIT
       ENDIF
    ENDDO ! cntl%cdft HDA
    IF (cdftlog%thda)THEN
       IF (cntl%tdipd) CALL localize(tau0,c0,c2,sc0,crge%n)
       CALL zcopy(ldim,rswf,1,rsave2,1)
       CALL zeroing(wd)!,nnr1)
       CALL dcopy(fpar%nnr1,wdiff,1,wd,1)

       IF (wgaussl%thdas)THEN
          CALL icopy(cdftpi%naccr,cdftpi%cdft_a,1,cdftpi%spd,1)
          CALL icopy(cdftpi%naccr,cdftpi%cdft_d,1,cdftpi%cdft_a,1)
          CALL icopy(cdftpi%naccr,cdftpi%spd,1,cdftpi%cdft_d,1)
          CALL cdft_adarrays
          cdftlog%recw=.TRUE.
          CALL cdft_w(rhoe,tau0,chrg_)
       ELSE IF (wgaussl%thdawm)THEN
          CALL dcopy(fpar%nnr1,wdsave,1,wdiff,1)
       ENDIF

       n2save(1)=cdftcom%cdft_nc+cdftcom%vgrad(1)
       n2save(2)=cdftcom%cdft_ns+cdftcom%vgrad(2)
       CALL transme(csave,c0,rsave,rsave2,esave,ener_com%etot,vsave,cdftcom%cdft_v,&
            nsave,n2save)
       DEALLOCATE(csave,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(rsave,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(rsave2,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (rout1%teband) THEN
       CALL enbandpri(eigv,crge%f,ener_com%amu,crge%n,nkpt%nkpts)
    ENDIF
    IF (tkpts%tonlydiag) THEN
       GOTO 150
    ENDIF
    ! Calculate ionic forces
    IF (ropt_mod%convwf.AND.(tfor.OR.cntl%tpres.OR.vdwl%vdwd)) THEN
       IF (vdwl%vdwd) CALL localize(tau0,c0,c2,sc0,crge%n)
       IF (cntl%tdiag) THEN
          CALL updrho(c0,c2,gde,sc0,pme,vpp,tau0,fion,eigv,&
               rhoe,psi,&
               crge%n,.TRUE.,.FALSE.,cntl%tpres,infi,thl,nhpsi)
       ELSE
          ropt_mod%calste=cntl%tpres
          IF (tkpts%tkpnt) THEN
             CALL k_updwf(c0,c2,sc0,tau0,fion,pme,gde,vpp,eigv,&
                  rhoe,psi,crge%n,.TRUE.,int_update_pot)
          ELSE
             update_pot=.TRUE.
             int_update_pot=1
             CALL updwf(c0(:,:,1),c2,sc0(:,:,1),tau0,fion,pme,gde,vpp,eigv,&
                  rhoe,psi,crge%n,.TRUE.,update_pot)
          ENDIF
       ENDIF
       ! Calculate norm of nuclear gradient
       CALL gsize(fion,gnmax,gnorm)
       IF (cntl%tpres) THEN
          CALL totstr
       ENDIF
    ELSE
       gnmax=0.0_real_8
       gnorm=0.0_real_8
    ENDIF
150 CONTINUE
    !
    IF (lspin2%teprof.AND.lspin2%tlse) THEN
       CALL lseprof(c0,rhoe,psi,tau0,fion,crge%n)
    ENDIF
    ! EHR[
    IF (cntl%tpspec) THEN
       IF (rout1%rhoout) CALL rhopri(c0,tau0,rhoe,psi(:,1),crge%n,nkpt%nkpnt)
       CALL mm_dim(mm_go_mm,statusdummy)
       CALL zhwwf(2,irec,c0,c2,crge%n,eigv,tau0,velp,taup,iteropt%nfi)
    ENDIF
    ! EHR]
    !
    ! Calculate electric field gradient and EPR anisotropic tensors
    CALL epr_efg(rhoe,psi,0)

    CALL mm_dim(mm_go_mm,statusdummy)
    ! ==--------------------------------------------------------------==
    ! Calculate non-adiabatic couplings if required
    IF (tcpl) THEN
       IF (paral%parent) CALL finalp(tau0,fion,tau0,eigv)
       CALL cplngs(c0,c2,cf_4d,c0ini,gde,sc0,pme,vpp,eigv(:,1),&
            rhoe,rhoini,psi,tau0,fion,infi,&
            irec,.TRUE.,tfor,cplion,cplcfg,crge%n)
       CALL prcplngs(cplion,cplcfg,.FALSE.)
       CALL prcplngs(cplion,cplcfg,.TRUE.)
       IF (tfor.AND.tcplfd.OR.tcpllr) THEN
          IF (tcplfd) irec(irec_phes) = 1
          CALL zhwwf(2,irec,c0,cf_4d,crge%n,eigv,tau0,velp,taup,iteropt%nfi)
       ENDIF
       DEALLOCATE(cplion,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(cplcfg,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       IF (cntl%tdiag.OR.tcpllr) DEALLOCATE(rhoini,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       IF (.NOT.tallat) THEN
          DEALLOCATE(iatfd,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
       ENDIF
       IF (tfor.AND.tcplfd.AND.paral%parent) THEN
          DEALLOCATE(hesscr,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
       ENDIF
       DEALLOCATE(c0ini,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(cf_4d,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    IF (paral%parent) THEN
       CALL prmem('    RWFOPT')
       IF (.NOT.cntl%bsymm.OR.ropt_mod%convwf) CALL finalp(tau0,fion,tau0,eigv)
       CALL geofile(tau0,fion,'WRITE')
    ENDIF
    ! CB: TAUP,FION handled by BS_WFO
    IF (.NOT.cntl%bsymm) THEN
       DEALLOCATE(taup,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(fion,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    IF (cntl%tdiag) THEN
       DEALLOCATE(rin0,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(rout0,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(rmix,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ELSEIF (tkpts%tkpnt) THEN
       DEALLOCATE(rin0,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(rmix,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    DEALLOCATE(rhoe,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    !#if defined(_HAS_CUDA)
    !    CALL cuda_dealloc_host(psi)
    !#else
    DEALLOCATE(psi,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    !#endif
    DEALLOCATE(scr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! EHR[
    IF (cntl%tpdist) THEN
       DEALLOCATE(extf,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    IF (cntl%tpspec.OR.cntl%tpdist) THEN
       DEALLOCATE(norms,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! EHR]
    !FIXME error checking
    IF(cntl%thubb)DEALLOCATE(c2u)
    ! ==--------------------------------------------------------------==
    CALL mm_dim(mm_revert,oldstatus)
    CALL tihalt(procedureN,isub)
    RETURN
  END SUBROUTINE rwfopt
  ! ==================================================================
  SUBROUTINE give_scr_rwfopt(lrwfopt,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lrwfopt
    CHARACTER(len=30)                        :: tag

    CHARACTER(*), PARAMETER                  :: procedureN = 'give_scr_rwfopt'

    INTEGER                                  :: isub, lcplngs, ldipd, &
                                                linitrun, lprof, lrhopri, &
                                                lupd, lupdrho, nstate

    CALL tiset(procedureN,isub)
    CALL give_scr_initrun(linitrun,tag)
    nstate=crge%n
    lupd=0
    lrhopri=0
    lupdrho=0
    lcplngs=0
    IF (cntl%tdiag) THEN
       CALL give_scr_updrho(lupd,tag,nstate,.TRUE.,cntl%tpres)
    ELSE
       IF (tkpts%tkpnt.AND.(.NOT.(cntl%tpspec.OR.cntl%tpdist))) THEN
          CALL give_scr_kupdwf(lupd,tag,nstate,.TRUE.)
       ELSE
          CALL give_scr_updwf(lupd,tag,nstate,.TRUE.)
       ENDIF
       IF (cntl%tdiagopt) THEN
          CALL give_scr_updrho(lupdrho,tag,nstate,.TRUE.,cntl%tpres)
          lupd = MAX(lupdrho,lupd)
       ENDIF
    ENDIF
    lprof=0
    IF (lspin2%tlse.AND.lspin2%teprof) THEN
       lprof=4*ncpw%ngw+100
    ENDIF
    IF (rout1%rhoout) THEN
       CALL give_scr_rhopri(lrhopri,tag,nstate)
    ENDIF
    IF (cntl%tdipd.OR.vdwl%vdwd) THEN
       CALL give_scr_ddipo(ldipd,tag)
    ELSE
       ldipd=0
    ENDIF
    IF (tcpl) CALL give_scr_cplngs(lcplngs,tag)
    lrwfopt=MAX(linitrun,lupd,lrhopri,ldipd,lcplngs)+lprof+100
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_rwfopt
  ! ==================================================================
  SUBROUTINE lseprof(c0,rhoe,psi,tau0,fion,nstate)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rhoe(:,:)
    COMPLEX(real_8)                          :: psi(:,:)
    REAL(real_8)                             :: tau0(:,:,:), fion(:,:,:)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(ncpw%ngw,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'lseprof'

    CHARACTER(len=30)                        :: tag
    COMPLEX(real_8), ALLOCATABLE             :: ca(:), cb(:), sca(:)
    INTEGER                                  :: ia, ierr, lsca
    LOGICAL                                  :: debug, tccc
    REAL(real_8)                             :: ang, angd, cosa, deee, ezero, &
                                                sina, tpi

    IF (lspin2%tcas22) THEN
       lspin2%troks=.TRUE.
       lspin2%tcas22=.FALSE.
       tccc=.TRUE.
    ELSE
       tccc=.FALSE.
    ENDIF
    tpi=0.5_real_8*fpi
    debug=.FALSE.

    ALLOCATE(ca(ncpw%ngw),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(cb(ncpw%ngw),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    CALL give_scr_rscpot(lsca, tag, .FALSE.)
    ALLOCATE(sca(lsca),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    CALL dcopy(2*ncpw%ngw,c0(1,clsd%ialpha),1,ca,1)
    CALL dcopy(2*ncpw%ngw,c0(1,clsd%ibeta),1,cb,1)
    IF (paral%io_parent) THEN
       WRITE(6,'(/,1X,64("*"))')
       WRITE(6,*)
       WRITE(6,'(T20,A,T60,A)') '  Rotation Angle ','Energy'
    ENDIF
    lspin3%rotab=lspin3%rotab/tpi*360._real_8
    DO ia=0,9
       ang=ia*5._real_8
       angd=ang/360._real_8*tpi
       cosa=COS(angd)
       sina=SIN(angd)
       CALL dcopy(2*ncpw%ngw,ca,1,c0(1,clsd%ialpha),1)
       CALL dscal(2*ncpw%ngw,cosa,c0(1,clsd%ialpha),1)
       CALL daxpy(2*ncpw%ngw,sina,cb,1,c0(1,clsd%ialpha),1)
       CALL dcopy(2*ncpw%ngw,cb,1,c0(1,clsd%ibeta),1)
       CALL dscal(2*ncpw%ngw,cosa,c0(1,clsd%ibeta),1)
       CALL daxpy(2*ncpw%ngw,-sina,ca,1,c0(1,clsd%ibeta),1)

       CALL rscpot(c0,tau0,fion,rhoe,psi,&
            .FALSE.,.FALSE.,nstate,1)

       IF (ia.EQ.0) THEN
          ezero=ener_com%etot
       ENDIF
       deee=ener_com%etot-ezero
       IF (paral%io_parent) THEN
          WRITE(6,'(T28,F8.0,T50,F16.8)') lspin3%rotab+ang,deeE
       ENDIF
    ENDDO
    IF (paral%io_parent) THEN
       WRITE(6,'(1X,64("*"))')
       WRITE(6,'(A,T50,F16.8)') "  LAGRANGE MULTIPLIERS : A-B ",lspin3%hablse
       WRITE(6,'(A,T50,F16.8)') "  LAGRANGE MULTIPLIERS : A-O ",lspin3%haolse
       WRITE(6,'(A,T50,F16.8)') "  LAGRANGE MULTIPLIERS : B-O ",lspin3%hbolse
       WRITE(6,'(1X,64("*"))')
       WRITE(6,*)
    ENDIF

    DEALLOCATE(ca,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(cb,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(sca,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)

    IF (tccc) THEN
       lspin2%tcas22=.TRUE.
       lspin2%troks=.FALSE.
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE lseprof
  ! ==================================================================

END MODULE rwfopt_utils
