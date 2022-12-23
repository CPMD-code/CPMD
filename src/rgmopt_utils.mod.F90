MODULE rgmopt_utils
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

!!$public :: rgmopt
!!$public :: give_scr_rgmopt
!!$
!!$contains

END MODULE rgmopt_utils


! ==================================================================
SUBROUTINE rgmopt(c0,c1,c2,cm,sc0,pme,gde,vpp,eigv)
  ! ==--------------------------------------------------------------==
  ! == GEOMETRY OPTIMISATION                                        ==
  ! == XPAR contains TAU0 non fixed coordinates                     ==
  ! ==--------------------------------------------------------------==
#include "sizeof.h"
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE machine, ONLY: m_walltime
  USE mp_interface, ONLY: mp_bcast
  USE setsc_utils, ONLY : ihmat
  USE dynit_utils, ONLY : dynit
  USE symtrz_utils, ONLY : symvec, symmat, give_scr_symvec, give_scr_symmat, give_scr_symrho, symrho
  USE chksym_utils, ONLY : updatsym
  USE lsforce_utils, ONLY : lsforce
  USE mm_qmmm_forcedr_utils, ONLY : mm_qmmm_forcedr
  USE wrccfl_utils, ONLY : wrccfl
  USE sdcell_utils, ONLY : sdcell
  USE puttau_utils, ONLY : puttau, gettau, taucl
  USE rrfo_utils, ONLY : rrfo, give_scr_rrfo
  USE hessup_utils, ONLY : hessup, hespow
  USE hessout_utils, ONLY : hessout
  USE hessin_utils, ONLY : hessin
  USE detdof_utils, ONLY : detdof
  USE sdion_utils, ONLY : sdion, give_scr_sdion
  USE adapttol_utils, ONLY : tol_init, tol_chk_force, tol_chk_cnvgrad, tol_chk_cnvener, tol_det_grad, tol_det_ener, tol_init_shared
  USE cnstfc_utils, ONLY : cnstfc
  USE rgdiis_utils, ONLY : rgdiis, give_scr_rgdiis
  USE cnstpr_utils, ONLY : cnstpr
  USE rbfgs_utils, ONLY : rbfgs, give_scr_rbfgs
  USE inr_dr_utils, ONLY : give_scr_inr,inr_dr
  USE rlbfgs_utils, ONLY : rlbfgs, rprfo, pvibana, stack_destroy, give_work_lbfgs, stack_init !, stack_do_sched
  USE parac, ONLY : paral,parai
  USE atwf , ONLY:tmovr
  USE spin , ONLY:clsd
  USE ener , ONLY:ener_com
  USE cnst , ONLY:au_kb
  USE elct , ONLY:crge
  USE tpar , ONLY:dt2bye,dt2_elec
  USE pslo , ONLY:pslo_com
  USE ions , ONLY:ions1
  USE soft , ONLY:soft_com
  USE norm , ONLY:cnorm,gemax,gnmax,gnorm
  USE ropt , ONLY:infi,infw,iteropt,ropt_mod
  USE coor , ONLY:fion,tau0,taup,velp
  USE sfac , ONLY:ddfnl
  USE cotr , ONLY:cotc0,dtm,hess
  USE nlcc , ONLY:corel
  USE andr , ONLY:andr2
  USE andp , ONLY:rin0,rmix,rout0
  USE nlps , ONLY:imagp,ndfnl
  USE fint , ONLY:fint1
  USE poin , ONLY:potr,rhoo
  USE kpts , ONLY:tkpts
  USE kpnt , ONLY:kdrho
  USE store_types , ONLY:cprint,iprint_force,restart1,rout1,store1
  USE metr , ONLY:metr_com
  USE symm , ONLY:symmi
  USE xinr , ONLY:correct,direct,inr_integer,inr_logical,tol_inr,zprec
  USE response_pmod , ONLY:dfnl00,fnl00,rho0,vofrho0
  USE implhv , ONLY:rs_v,sd0
  USE linres , ONLY:td01
  USE bsym , ONLY:bsclcs,bsfac,cnstwgt,resthswf
  USE bsympnt , ONLY:fnbs
  USE mm_dimmod , ONLY:mm_go_mm,mm_go_qm,mm_revert
  USE mm_input , ONLY:clc,lqmmm
  USE lscal , ONLY:deaccp,hesscr,ielstk,lvlhes,lwlbfgs,lwprfo,nsvib,st_put,wlbfgs
  USE efld , ONLY:extf,textfld
  USE cdftmod , ONLY:cdftci,cdftlog,wdiff
  USE system , ONLY:cnti,cntl,cntr,fpar,maxsys,nacc,ncpw,nkpt,parap,parm
  USE vdwcmod , ONLY:vdwl,vdwwfl
  USE hfxmod , ONLY:hfxc3
  USE wann , ONLY:wan05
  USE rinitwf_driver, ONLY : rinitwf
  USE setbsstate_utils, ONLY : setbsstate
  USE mm_dim_utils, ONLY : mm_dim
  USE bs_forces_diag_utils, ONLY : bs_forces_diag
  USE newcell_utils, ONLY : newcell, give_scr_newcell
  USE totstr_utils, ONLY : totstr, dstre
  USE dum2_utils, ONLY : dumpr
  USE moverho_utils, ONLY : moverho,give_scr_moverho
  USE copot_utils, ONLY : give_scr_copot, copot
  USE updrho_utils, ONLY : give_scr_updrho,updrho
  USE forces_diag_utils, ONLY : forces_diag,updiag,dens_for_diag,give_scr_forces_diag
  USE calc_alm_utils, ONLY : calc_alm,calc_k_alm,give_scr_calc_alm
  USE extrap_utils, ONLY : extrap
  USE localize_utils, ONLY : localize,localize2,wc_print
  USE lr_tddft_utils, ONLY : give_scr_lr_tddft, lr_tddft
  USE updwf_utils, ONLY : give_scr_updwf, updwf
  USE cdft_utils, ONLY :  wcdft_restart, cdft_forces, cdft_finalize, write_w, vupdate, cdft_w, cdft_reinit, init_cdft
  USE k_updwf_utils, ONLY : k_updwf
  USE ortho_utils, ONLY : ortho,give_scr_ortho
  USE rhopri_utils, ONLY : rhopri 
  USE finalp_utils, ONLY : finalp
  USE phfac_utils, ONLY : phfac
  USE wrener_utils, ONLY : wrprint_wfopt,wrener,wrprint_geopt,wrcell
  USE forcep_utils, ONLY : rhoe_psi_size
  USE wrener_utils, ONLY : wreigen,wrprint_wfopt
  USE wrgeo_utils, ONLY : wrgeof
  USE gsize_utils, ONLY : gnodim
  USE hesele_utils, ONLY : hesele
  USE geofile_utils, ONLY : geofile
  USE setirec_utils, ONLY : write_irec
  USE wv30_utils, ONLY : zhwwf
  USE initrun_driver, ONLY : initrun
  USE initrun_utils, ONLY : give_scr_initrun
  USE testex_utils, ONLY : testex
  USE rnlsm_utils, ONLY : rnlsm
  USE forcedr_driver, ONLY : forcedr
  USE zeroing_utils,                   ONLY: zeroing
  IMPLICIT NONE
  ! Arguments
  COMPLEX(real_8) :: c0(nkpt%ngwk,crge%n,nkpt%nkpnt*bsfac),c2(nkpt%ngwk,crge%n,bsfac),c1(*),&
       cm(*),sc0(nkpt%ngwk,crge%n,bsfac),pme(*),gde(*)
  REAL(real_8) :: vpp(ncpw%ngw),eigv(crge%n,nkpt%nkpts*bsfac)
  ! Variables
  LOGICAL :: tfor,cvbswf,cvhswf
  INTEGER :: irec(100),nstate,nnx,nhpsi,nx,ifcalc,i,ig,j,&
       lscr,nmm,linr

  REAL(real_8), ALLOCATABLE :: tscr(:,:,:)
  REAL(real_8), ALLOCATABLE :: xpar(:)
  REAL(real_8), ALLOCATABLE :: dxpar(:)
  REAL(real_8), ALLOCATABLE :: ypar(:)
  REAL(real_8), ALLOCATABLE :: dypar(:)
  REAL(real_8), ALLOCATABLE :: gold(:)
  REAL(real_8), ALLOCATABLE :: told(:)
  REAL(real_8), ALLOCATABLE :: rm1(:)
  REAL(real_8), ALLOCATABLE :: rinp(:)

  COMPLEX(real_8), ALLOCATABLE :: psi(:,:)

  REAL(real_8), ALLOCATABLE :: rhoe(:,:)
  REAL(real_8), ALLOCATABLE :: scr(:)
  ! 
  CHARACTER (len=30) :: tag
  REAL(real_8) :: dummy(2),time1,time2,time11,time22,tcpu,ekincp,&
       ekin1,ekin2,etoto,etotg,detot,etotbs,etoths,&
       temp1,temp2,cnmax,pf1,dstress,dsmax,thl(2),&
       dx,dxmax,ekinh1,ekinh2
  ! ... inr...
  COMPLEX(real_8), ALLOCATABLE :: eirop(:)
  COMPLEX(real_8), ALLOCATABLE :: eivps(:)
  COMPLEX(real_8), ALLOCATABLE :: c11(:,:)

  REAL(real_8), ALLOCATABLE :: z11(:,:)
  REAL(real_8), ALLOCATABLE :: drhoe(:,:)

  ! real(8) :: ddfnl(imagp,nat,maxsys%nhxs,3,3,n,nkpnt)
  REAL(real_8), ALLOCATABLE :: tmp(:)
  REAL(real_8), ALLOCATABLE :: tras(:,:)

  INTEGER :: ldfnl



  ! pointer (ip_ddfnl,ddfnl)
  ! 
  LOGICAL :: status,statusdummy,update_pot,qmmmrest
  ! Linear scaling
  REAL(real_8) :: genvmx,wcmpsv
  LOGICAL :: lret,lcvmsk,ltolad,lpexit
  INTEGER :: lrlbfgs,lrprfo


  ! rhoo_1d is alias for rhoo for allocation through memory90 
  REAL(real_8), ALLOCATABLE :: rhoo_1d(:)

  ! integer ia, is, k,iat
  INTEGER :: cdfti
  LOGICAL :: convtest
  INTEGER :: isub,ierr,il_rhoe_1d, il_rhoe_2d, il_psi_1d, il_psi_2d
  CHARACTER(*),PARAMETER :: procedureN='RGMOPT'
  ! ==================================================================
  CALL tiset(procedureN,isub)
  ! Check some options
  IF (cntl%tddft.AND.cntl%tdiag) THEN
     CALL stopgm("RGMOPT","TDDFT.AND.TDIAG NOT POSSIBLE",& 
          __LINE__,__FILE__)
  ENDIF
  ! ==--------------------------------------------------------------==
  time1 =m_walltime()
  ! cntl%tsymrho=.FALSE.           ! do not symmetrize density
  nstate=crge%n
  CALL mm_dim(mm_go_mm,status)
  ! ==--------------------------------------------------------------==
  ! ALLOCATIONS
  ALLOCATE(fion(3,maxsys%nax,maxsys%nsx),STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
       __LINE__,__FILE__)
  ALLOCATE(tscr(3,maxsys%nax,maxsys%nsx),STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
       __LINE__,__FILE__)
  ALLOCATE(taup(3,maxsys%nax,maxsys%nsx),STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
       __LINE__,__FILE__)
  CALL zeroing(tscr)!,3*maxsys%nax*maxsys%nsx)
  CALL zeroing(taup)!,3*maxsys%nax*maxsys%nsx)
  CALL zeroing(velp)!,3*maxsys%nax*maxsys%nsx)
  CALL zeroing(fion)!,3*maxsys%nax*maxsys%nsx)
  qmmmrest = .TRUE.
  IF (lqmmm%qmmm)THEN
     IF (textfld)THEN
        ALLOCATE(extf(fpar%kr1*fpar%kr2s*fpar%kr3s),STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
             __LINE__,__FILE__)
        CALL zeroing(extf)!,kr1*kr2s*kr3s)
     ENDIF
     IF (.NOT.paral%qmnode) qmmmrest = .FALSE.
  ENDIF
  nacc = 7
  iteropt%nfi  = 0
  IF (cntl%tdiag.AND..NOT.cntl%bsymm) THEN
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
  ELSEIF (tkpts%tkpnt .AND. .NOT. cntl%tdiag .AND. .NOT. cntl%bsymm) THEN
     nnx=fpar%nnr1*clsd%nlsd
     ALLOCATE(rin0(fpar%nnr1,clsd%nlsd),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     ALLOCATE(rmix(fpar%nnr1,clsd%nlsd),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
  ENDIF
  CALL rhoe_psi_size(il_rhoe_1d=il_rhoe_1d, il_rhoe_2d=il_rhoe_2d,&
       il_psi_1d=il_psi_1d, il_psi_2d=il_psi_2d)
  ALLOCATE(rhoe(il_rhoe_1d, il_rhoe_2d),STAT=ierr)
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
     IF (td01%ns_tri.GT.0) THEN
        nmm=2
     ENDIF
  ENDIF
  IF (cntl%bsymm) THEN
     ALLOCATE(fnbs(3*maxsys%nax*maxsys%nsx),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
  ENDIF
  il_psi_1d = nmm* il_psi_1d
  ALLOCATE(psi(il_psi_1d, il_psi_2d),STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
       __LINE__,__FILE__)
  ! ==--------------------------------------------------------------==
  ! TIME STEP FUNCTIONS
  CALL dynit(ekincp,ekin1,ekin2,temp1,temp2,ekinh1,ekinh2)
  ropt_mod%convge=.FALSE.
  ropt_mod%modens=.FALSE.
  ropt_mod%engpri=.FALSE.
  ropt_mod%calste=.FALSE.
  ! DFNL
  ndfnl=parap%nst12(parai%mepos,2)-parap%nst12(parai%mepos,1)+1
  ldfnl=imagp*3*ions1%nat*maxsys%nhxs*ndfnl*nkpt%nkpnt
  IF (ldfnl.LE.0) THEN
     ldfnl=1
  ENDIF
  ! ==--------------------------------------------------------------==
  IF (cntl%tinr) THEN
     ALLOCATE(rho0(il_rhoe_1d,il_rhoe_2d),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     CALL zeroing(rho0)!,il_rhoe)

     ALLOCATE(drhoe(il_rhoe_1d,il_rhoe_2d),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     CALL zeroing(drhoe)!,il_rhoe)
     ! potentials and related:
     ALLOCATE(VofRho0(il_rhoe_1d,il_rhoe_2d),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     CALL zeroing(VofRho0)!,il_rhoe)

     ALLOCATE(eirop(ncpw%nhg),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     CALL zeroing(eirop)!,nhg)

     ALLOCATE(eivps(ncpw%nhg),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     CALL zeroing(eivps)!,nhg)

     ALLOCATE(c11(ncpw%ngw,crge%n*nkpt%ngwk/ncpw%ngw),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     CALL zeroing(c11)!,nkpt%ngwk*n)

     ! others & scratch:
     ALLOCATE(z11(crge%n,crge%n),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     CALL zeroing(z11)!,nstate*nstate)

     ALLOCATE(rs_v(3*ions1%nat,2*ions1%nat/ions1%nat),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     CALL zeroing(rs_v)!, 2*3*ions1%nat)

     ALLOCATE(tmp(3*ions1%nat),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     CALL zeroing(tmp)!,3*ions1%nat)

     ALLOCATE(fnl00(imagp,nstate,ions1%nat,maxsys%nhxs,nkpt%nkpnt),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     CALL zeroing(fnl00)!,imagp*nstate*ions1%nat*maxsys%nhxs*nkpt%nkpnt)

     IF (ldfnl.GT.0) THEN
        ALLOCATE(dfnl00(imagp,3,ions1%nat,maxsys%nhxs,ndfnl,nkpt%nkpnt),STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
             __LINE__,__FILE__)
     ELSE
        ALLOCATE(dfnl00(1,1,1,1,1,1),STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
             __LINE__,__FILE__)
     ENDIF
     CALL zeroing(dfnl00)!,ldfnl)

     ALLOCATE(ddfnl(3*ldfnl,1,1,1),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     CALL zeroing(ddfnl)!,3*ldfnl)

     ALLOCATE(sd0(3*ions1%nat,3*ions1%nat),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     CALL zeroing(sd0)!,9*ions1%nat*ions1%nat)

     ALLOCATE(tras(3*ions1%nat,3*ions1%nat/ions1%nat),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     CALL zeroing(tras)!,9*ions1%nat)

     ALLOCATE(direct(3*ions1%nat),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     CALL zeroing(direct)!,3*ions1%nat)
     ALLOCATE(correct(3*ions1%nat),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     CALL zeroing(correct)!,3*ions1%nat)
     ALLOCATE(zprec(3*ions1%nat),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     CALL zeroing(zprec)!,3*ions1%nat)

  ENDIF
  ! ==--------------------------------------------------------------==
  ! == INITIALIZATION                                               ==
  ! ==--------------------------------------------------------------==
  ! We allocate only for initrun
  IF (cntl%bsymm.AND.paral%parent) THEN
     IF (paral%io_parent) THEN
        WRITE(6,*)
        WRITE(6,*) 'BSYMM: BS WAVEFUNCTION INITIALIZATION'
     ENDIF
  ENDIF
  IF (cntl%cdft)CALL init_cdft()
  CALL initrun(irec,c0,c2,sc0,rhoe,psi,eigv)
  IF (cntl%bsymm) THEN
     IF (paral%io_parent)&
          WRITE(6,*) 'BSYMM: HS WAVEFUNCTION INITIALIZATION'
     bsclcs=2
     CALL setbsstate
     CALL initrun(irec,c0(:,:,2:),c2(1,1,2),sc0(1,1,2),rhoe,psi,&
          eigv(1,2))
  ENDIF
  CALL write_irec(irec)
  ! ==--------------------------------------------------------------==
  IF (cntl%lbfgs.OR.cntl%prfo) THEN
     CALL stack_init(ielstk,2*nkpt%ngwk*nstate,1,st_put)
  ENDIF
  IF (paral%parent) THEN
     ! Determination of NODIM
     CALL detdof(tau0,tscr)
  ENDIF
  CALL mp_bcast(cotc0%nodim,parai%source,parai%allgrp)
  IF (cntl%tinr) THEN
     CALL mp_bcast_byte(inr_integer, size_in_bytes_of(inr_integer),parai%source,parai%allgrp)
     CALL mp_bcast(tol_inr,parai%source,parai%allgrp)
  ENDIF
  ! ==--------------------------------------------------------------==
  ! UPDATE SYMMETRIES (POINT GROUP)
  CALL updatsym(tau0)
  ! We have all information to do that (specially NODIM)
  CALL give_scr_rgmopt(lscr,tag)
  IF (cntl%tinr) THEN
     CALL give_scr_inr(linr,tag,nstate)
     lscr=MAX(lscr,linr)
  ENDIF
  ALLOCATE(scr(lscr),STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
       __LINE__,__FILE__)
  ! ==--------------------------------------------------------------==
  IF (cntl%tinr .OR. paral%parent) THEN ! Don t understand this [sebast]
     ALLOCATE(xpar(cotc0%nodim),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     ALLOCATE(dxpar(cotc0%nodim),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     ALLOCATE(ypar(cotc0%nodim),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     ALLOCATE(dypar(cotc0%nodim),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
  ENDIF
  IF ((cntl%gdiis.OR.cntl%tinr).AND.paral%parent) THEN
     ALLOCATE(gold(cotc0%nodim*cnti%mgdiis),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     ALLOCATE(told(cotc0%nodim*cnti%mgdiis),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
  ENDIF
  IF (paral%parent) THEN
     IF (cntl%tsdp) THEN
        cnti%npara=-1 !vw is that safe? shall it be bcast??
     ENDIF
     IF (cntl%lbfgs .AND. lwlbfgs.EQ.0) THEN
        CALL give_work_lbfgs(cotc0%nodim,lrlbfgs)! L-cntl%bfgs HESSIAN
        ALLOCATE(wlbfgs(lrlbfgs),STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
             __LINE__,__FILE__)
     ELSE IF (cntl%prfo .AND. lwprfo.EQ.0) THEN
        IF (lwlbfgs.EQ.0) THEN
           IF (lrlbfgs.GT.0) THEN
              ALLOCATE(wlbfgs(lrlbfgs),STAT=ierr)
              IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                   __LINE__,__FILE__)
           ENDIF
        ENDIF
     ELSE
        IF (.NOT.(cntl%tsdp.OR.cntl%lbfgs.OR.cntl%prfo)) THEN
           CALL hessin(tau0,restart1%rhe)
        ENDIF
     ENDIF
     CALL puttau(tau0,xpar)
     CALL dcopy(cotc0%nodim,xpar(1),1,ypar(1),1)
     CALL wrgeof(tau0,fion)
     IF (cntl%tprcp) THEN
        CALL wrcell
     ENDIF
  ENDIF
  ! ==--------------------------------------------------------------==
  ! == INITIALISATION OF DENSITY FOR DIAGONALISATION SCHEMES        ==
  ! ==--------------------------------------------------------------==
  ifcalc=0
  IF (cntl%tdavi) THEN
     nx=nkpt%ngwk*cnti%ndavv*nkpt%nkpnt+1
  ELSEIF (cntl%tdiag .AND. cntl%diis) THEN
     nx=((nkpt%ngwk*crge%n+8)*cnti%mdiis*nkpt%nkpnt)/4
  ELSE
     nx=1
  ENDIF
  IF (cntl%cdft)THEN
     CALL cdft_w(rhoe,tau0,dummy)
     IF (cntl%cdft_weight)CALL write_w(wdiff,"INIT")
  ENDIF
  IF (cntl%tdiag) THEN
     IF (paral%parent) THEN
        IF (paral%io_parent)&
             WRITE(6,'(1X,64("="))')
        IF (paral%io_parent)&
             WRITE(6,'(1X,"==",T25,A,T64,"==")')&
             'FORCES INITIALIZATION'
        IF (paral%io_parent)&
             WRITE(6,'(1X,64("="))')
     ENDIF
     ropt_mod%calste=cntl%tpres.OR.cntl%tprcp
     IF (cntl%bsymm) THEN
        CALL bs_forces_diag(crge%n,c0,c2,cm,sc0,cm(nx),vpp,eigv,&
             rhoe,psi,&
             tau0,tau0,tau0,fion,ifcalc,&
             irec,.TRUE.,.TRUE.)
     ELSE
        CALL forces_diag(crge%n,c0,c2,cm,sc0,cm(nx),vpp,eigv,&
             rhoe,psi,&
             tau0,tau0,tau0,fion,ifcalc,&
             irec,.TRUE.,.TRUE.)
     ENDIF
     IF (cntl%tprcp) THEN
        CALL dstre(metr_com%htfor,dstress,dsmax)
     ENDIF
     IF (ropt_mod%calste) THEN
        CALL totstr
     ENDIF
     IF (paral%parent) THEN
        IF (paral%io_parent)&
             WRITE(6,'(1X,64("="))')
        IF (paral%io_parent)&
             WRITE(6,'(1X,"==",T20,A,T64,"==")')&
             'END OF FORCES INITIALIZATION'
        IF (paral%io_parent)&
             WRITE(6,'(1X,64("="))')
     ENDIF
     IF (.NOT.cntl%bsymm) CALL dcopy(nnx,rin0,1,rm1,1)
  ENDIF
  ! ..-------------------------------------------------
  ! ..QM/MM mechanical coupling by Roethlisberger Group
  ! ..Initailize Molecular Mechanics subsystem:
  ! ..  (read MM_INPUT,initialize variables,
  ! ..       set up MM force field)
  ! ..-------------------------------------------------
  ! ..
#if defined (__QMECHCOUPL)
  IF (paral%parent) THEN
     CALL mm_cpmd_init (cntl%tqmmech ,tau0, ions0%na, ions1%nsp,&
          maxsys%nax, ions1%nat, restart1%rco, restart1%rvel, cntr%tempw)
  ENDIF
#endif
  ! ==--------------------------------------------------------------==
  ! == INITIALIZATION OF CONVERGENCE CRITERIA                       ==
  ! ==--------------------------------------------------------------==
  ! ELECTRONS
  CALL tol_init(ltolad,deaccp,genvmx)
  ! IONS
  lcvmsk = (.NOT.cntl%prfo)
  lpexit = .FALSE.
  ! ==--------------------------------------------------------------==
  ! == END OF INITIALIZATION                                        ==
  ! ==--------------------------------------------------------------==
  etoto=0.0_real_8
  etotg=0.0_real_8
  IF (paral%parent) THEN
     time2 =m_walltime()
     tcpu = (time2 - time1)*0.001_real_8
     IF (paral%io_parent)&
          WRITE(6,'(A,T50,F8.2,A8)') ' CPU TIME FOR INITIALIZATION',&
          tcpu,' SECONDS'
     IF (paral%io_parent)&
          WRITE(6,'(//,1X,64("="))')
     IF (paral%io_parent)&
          WRITE(6,'(1X,"=",T20,A,T65,"=")')&
          ' GEOMETRY OPTIMIZATION'
     IF (paral%io_parent)&
          WRITE(6,'(1X,64("="))')
     cnmax=0._real_8
     dxmax=0._real_8
     CALL wrprint_geopt(eigv,crge%f,ener_com%amu,nstate,tau0,fion,&
          ener_com%etot,etotg,tcpu,gemax,gnmax,dxmax,gnorm,cnmax,&
          dstress,dsmax,.FALSE.,iteropt%nfi,infi)
  ENDIF
  ropt_mod%calste=.FALSE.
  ! ==================================================================
  ! ==      THE BASIC LOOP FOR GEOMETRY OPTIMIZATION                ==
  ! ==================================================================
  IF (cntl%tsdp.AND.cntl%tsde) THEN
     IF (cntl%tddft) THEN
        CALL stopgm("RGMOPT","cntl%tddft NOT POSSIBLE",& 
             __LINE__,__FILE__)
     ENDIF
     ! ==--------------------------------------------------------------==
     ! == COMBINED STEEPEST DESCENT                                    ==
     ! ==--------------------------------------------------------------==
     IF (tkpts%tkpnt) THEN
        CALL stopgm('RGMOPT','K POINTS NOT IMPLEMENTED',& 
             __LINE__,__FILE__)
     ENDIF
     IF (lqmmm%qmmm) THEN
        CALL stopgm('RGMOPT','Qmmm + steepest not tested',& 
             __LINE__,__FILE__)
     ENDIF
     ropt_mod%calste=cntl%tprcp
     etoto=0.0_real_8
     DO infi=1,cnti%nomore
        time1=m_walltime()
        iteropt%nfi=iteropt%nfi+1
        ifcalc=infi
        ! ..       -------------------------------------------------
        ! ..       QM/MM mechanical coupling by Roethlisberger Group
        ! ..       OPTIMIZE QM/MM subsystem with QM atoms frozen
        ! ..       -------------------------------------------------
#if defined (__QMECHCOUPL)
        IF (paral%parent.AND.cntl%tqmmech) THEN
           CALL mm_cpmd_geopt (tau0, ions0%na, ions1%nsp, maxsys%nax, ions1%nat)
        ENDIF
#endif
        ! AK 2005/05/06 the following line was deleted out in the cntl%bsymm
        ! version, but it should not interfere. 
        IF (cntl%tdipd.OR.vdwl%vdwd) THEN
           CALL localize2(tau0,c0,c2,sc0,nstate)
           vdwwfl%twannup=.TRUE.
        ENDIF
        ener_com%ecnstr=0.0_real_8
        ener_com%erestr=0.0_real_8
        ! Initialize forces
        IF (lqmmm%qmmm)THEN
           CALL mm_qmmm_forcedr(c0,c2,sc0,rhoe,psi,tau0,fion,&
                eigv,nstate,0,.FALSE.,.TRUE.,.TRUE.)
        ELSE
           CALL forcedr(c0(:,:,1),c2(:,:,1),sc0(:,:,1),rhoe,psi,tau0,fion,eigv,&
                nstate,1,.TRUE.,.TRUE.)
        ENDIF
        IF (ropt_mod%calste) THEN
           CALL totstr
        ENDIF
        IF (paral%parent) THEN
           CALL dscal(3*maxsys%nax*maxsys%nsx,-1.0_real_8,fion(1,1,1),1)
           CALL puttau(fion,dxpar)
           CALL cnstfc(tau0,tscr,dxpar,cnmax,.FALSE.)
           ! CHECK GRADIENT
           CALL gnodim(dxpar,cotc0%nodim,gnmax,gnorm)
           IF (cntl%tprcp) THEN
              CALL dstre(metr_com%htfor,dstress,dsmax)
              ropt_mod%convge=gnmax.LT.cntr%tolng.AND.dstress.LT.cntr%tolcg
           ELSE
              ropt_mod%convge=gnmax.LT.cntr%tolng
           ENDIF
        ENDIF
        IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)
        CALL testex(ropt_mod%convge)! check if any node is converged
        IF (ropt_mod%convge .OR. soft_com%exsoft) GOTO 200
        ! UPDATE THE NUCLEAR POSITIONS
        IF (paral%parent) THEN
           IF (MOD(ifcalc-1,cprint%iprint_step).EQ.0) THEN
              IF (cprint%iprint(iprint_force).EQ.1) THEN
                 CALL wrgeof(tau0,fion)
              ENDIF
#if defined (__QMECHCOUPL)
              IF (cntl%tqmmech) THEN
                 CALL mm_cpmd_wrgeo
              ENDIF
#endif
           ENDIF
           CALL puttau(tau0,xpar)
           CALL sdion(xpar,dxpar)
           CALL gettau(tau0,xpar)
        ENDIF
        CALL mp_bcast(tau0,SIZE(tau0),parai%source,parai%allgrp)
        ! UPDATE SYMMETRIES (POINT GROUP)
        CALL updatsym(tau0)
        ! UPDATE WAVEFUNCTION 
        IF (cntl%prec) THEN
           CALL mm_dim(mm_go_qm,statusdummy)
           CALL hesele(dt2bye,vpp)
           CALL mm_dim(mm_revert,statusdummy)
           pf1=dt2_elec*cntr%hthrs/cntr%emass
           DO i=1,nstate
              DO ig=1,nkpt%ngwk
                 c0(ig,i,1)=c0(ig,i,1)+vpp(ig)*pf1*c2(ig,i,1)
              ENDDO
           ENDDO
        ELSE
           CALL daxpy(2*nkpt%ngwk*nstate,dt2bye,c2(1,1,1),1,c0(1,1,1),1)
        ENDIF
        IF (qmmmrest.AND.MOD(ifcalc,store1%isctore).EQ.0) THEN
           CALL zhwwf(1,irec,c0,c2,crge%n,eigv,tau0,velp,tau0,iteropt%nfi)
        ENDIF
        ! Move the cell
        IF (cntl%tsdc) THEN
           CALL sdcell(tau0)
        ENDIF
        CALL mp_bcast(tau0,SIZE(tau0),parai%source,parai%allgrp)
        IF (cntl%tprcp) THEN
           CALL mp_bcast_byte(metr_com, size_in_bytes_of(metr_com),parai%source,parai%allgrp)
           DO i=1,3
              parm%a1(i) = metr_com%ht(1,i)
              parm%a2(i) = metr_com%ht(2,i)
              parm%a3(i) = metr_com%ht(3,i)
           ENDDO
           CALL ihmat(metr_com%ht,metr_com%htm1,parm%omega)
           CALL newcell
        ENDIF
        ! REORTHOGONALIZE WAVEFUNCTIONS
        CALL mm_dim(mm_go_qm,statusdummy)
        CALL phfac(tau0)
        IF (corel%tinlc) THEN
           CALL copot(rhoe,psi,cntl%tprcp)
        ENDIF
        IF (.NOT.cntl%nonort) THEN
           IF (pslo_com%tivan) THEN
              CALL rnlsm(c0(:,:,1),nstate,1,1,.FALSE.)
           ENDIF
           CALL ortho(nstate,c0(:,:,1),c2)
        ENDIF
        CALL mm_dim(mm_revert,statusdummy)
        ! Printout the evolution of the iterative optimization
        IF (paral%parent) THEN
           detot=ener_com%etot-etoto
           IF (infi.EQ.1) detot=0.0_real_8
           time2=m_walltime()
           tcpu=(time2-time1)*0.001_real_8
           IF (cntl%tprcp) THEN
              IF (paral%io_parent)&
                   WRITE(6,*) '  TOTAL STRESS - P (kB): Step:', iteropt%nfi
              DO i=1,3
                 IF (paral%io_parent)&
                      WRITE(6,'(6X,3(F18.8))')&
                      ((metr_com%htfp(i,j)/parm%omega)*au_kb,j=1,3)
              ENDDO
           ENDIF
           IF (MOD(ifcalc-1,cprint%iprint_step).EQ.0) THEN
              CALL wrener
              IF (cntl%tsdc) THEN
                 IF (paral%io_parent)&
                      WRITE(6,'(A,A)') ' NFI    GNMAX    GEMAX   HFOR  ',&
                      '         ETOT         DETOT      TCPU'
              ELSE
                 IF (paral%io_parent)&
                      WRITE(6,'(A)')&
                      ' NFI    GNMAX    GEMAX          ETOT         DETOT      TCPU'
              ENDIF
           ENDIF
           IF (cntl%tprcp) THEN
              CALL dstre(metr_com%htfor,dstress,dsmax)
              IF (paral%io_parent)&
                   WRITE(6,'(I4,F9.6,F9.6,F9.6,F14.5,F14.6,F10.2) ')&
                   iteropt%nfi,gnmax,gemax,dstress,ener_com%etot,detot,tcpu
           ELSE
              IF (paral%io_parent)&
                   WRITE(6,'(I4,F9.6,F9.6,F14.5,F14.6,F10.2) ')&
                   iteropt%nfi,gnmax,gemax,ener_com%etot,detot,tcpu
           ENDIF
           etoto=ener_com%etot
        ENDIF
        IF (ifcalc.GE.cnti%nomore_iter) GOTO 200
     ENDDO
  ELSEIF (cntl%simul) THEN
     CALL stopgm('RGMOPT','SIMUL NOT WORKING',& 
          __LINE__,__FILE__)
  ELSE
     ! ==--------------------------------------------------------------==
     ! == CONVENTIONAL GEOMETRY OPTIMIZATION                           ==
     ! == SETUP ATOMIC BASIS TO EXPAND DENSITY                         ==
     ! ==--------------------------------------------------------------==
     CALL dcopy(3*maxsys%nax*maxsys%nsx,tau0(1,1,1),1,taup(1,1,1),1)
     etoto=0.0_real_8
     etotg=0.0_real_8
     ifcalc=0
     DO infi=1,cnti%nomore
        time1=m_walltime()
        ener_com%ecnstr=0.0_real_8
        ener_com%erestr=0.0_real_8
        iteropt%nfi=iteropt%nfi+1
        ropt_mod%convwf=.FALSE.
        cvbswf=.FALSE.
        cvhswf=.FALSE.
        ropt_mod%sdiis=.TRUE.
        ropt_mod%spcg=.TRUE.
        IF (cntl%cdft)CALL cdft_reinit()
        ! ..       -------------------------------------------------
        ! ..       QM/MM mechanical coupling by Roethlisberger Group
        ! ..       OPTIMIZE QM/MM subsystem with QM atoms frozen
        ! ..       -------------------------------------------------
#if defined (__QMECHCOUPL)
        IF (paral%parent.AND.cntl%tqmmech) THEN
           CALL mm_cpmd_geopt (tau0, ions0%na, ions1%nsp, maxsys%nax, ions1%nat)
        ENDIF
#endif
        ! Preoptimization of the wavefunction
        ropt_mod%calste=(cntl%tpres.AND.MOD(infi,cnti%npres).EQ.0).OR.cntl%tprcp
        IF (cntl%tdiag) THEN
           IF (infi.EQ.1) THEN
              GOTO 100
           ENDIF
           IF (cntl%bsymm) THEN
              CALL bs_forces_diag(crge%n,c0,c2,cm,sc0,cm(nx),vpp,eigv,&
                   rhoe,psi,&
                   tau0,tau0,tau0,fion,ifcalc,&
                   irec,.TRUE.,.FALSE.)
           ELSE
              CALL forces_diag(crge%n,c0,c2,cm,sc0,cm(nx),vpp,eigv,&
                   rhoe,psi,&
                   tau0,tau0,tau0,fion,ifcalc,&
                   irec,.TRUE.,.FALSE.)
           ENDIF
           IF (soft_com%exsoft.OR.ifcalc.EQ.cnti%nomore_iter) THEN
              GOTO 200
           ENDIF
        ELSE
           ! ==--------------------------------------------------------==
           IF (cntl%cdft) THEN
              IF (cdftlog%recw)CALL cdft_w(rhoe,tau0,dummy)
              DO cdfti=1,cdftci%cdft_end
                 update_pot=.TRUE.
                 DO infw=1,cnti%nomore_iter
                    IF (infw.GT.1)update_pot=.FALSE.
                    time11=m_walltime()
                    ifcalc=ifcalc+1
                    ! Check if forces need to be calculated
                    CALL tol_chk_force(tfor,ropt_mod%calste,gemax,infw)
                    IF (clc%classical)tfor=.TRUE.
                    IF (cntl%tdipd.OR.vdwl%vdwd) THEN
                       IF (tfor.OR.hfxc3%twscr) THEN
                          CALL localize2(tau0,c0,c2,sc0,nstate)
                          vdwwfl%twannup=.TRUE.
                       ENDIF
                    ENDIF
                    ! Update the wavefunctions
                    IF (tkpts%tkpnt) THEN
                       IF (kdrho.LT.10._real_8*cntr%tolog) THEN
                          ropt_mod%calste = cntl%tprcp
                          tfor   = .TRUE.
                       ENDIF
                       CALL k_updwf(c0,c2,sc0,tau0,fion,pme,gde,vpp,eigv,&
                            rhoe,psi,nstate,tfor,infw)
                    ELSE
                       IF (.NOT. cntl%bsymm) THEN
                          CALL updwf(c0(:,:,1),c2(:,:,1),sc0(:,:,1),tau0,fion,pme,gde,vpp,eigv,&
                               rhoe,psi,nstate,tfor,update_pot)
                       ELSE
                          IF (.NOT.cvhswf) THEN
                             ! HS STATE
                             IF (infw.EQ.1) THEN
                                bsclcs=2
                                CALL setbsstate
                             ENDIF
                             CALL updwf(c0(:,:,2),c2(:,:,2),sc0(:,:,2),tau0,&
                                  fion,pme,gde,vpp,eigv(1,2),rhoe,&
                                  psi,nstate,tfor,update_pot)
                             cvhswf=ropt_mod%convwf
                          ELSEIF (.NOT.cvbswf) THEN
                             ! BS STATE
                             IF (bsclcs.EQ.2) THEN
                                IF (paral%parent) THEN
                                   IF (paral%io_parent)&
                                        WRITE(6,*) 'BSYMM: HS STATE CONVERGED. ',&
                                        'SWITCHING TO BS STATE.'
                                   etoto=0.0_real_8
                                   IF ((infi.EQ.1).AND.&
                                        (.NOT.resthswf))THEN
                                      IF (paral%io_parent)&
                                           WRITE(6,'(/,T2,2A,/)')&
                                           '!!!! WARNING : COPYING C0, EIGV & ',&
                                           'C2 OF HS TO BS !!!!'
                                      CALL dcopy(2*nkpt%ngwk*crge%n,c0(1,1,2),1,c0,1)
                                      CALL dcopy(crge%n,eigv(1,2),1,eigv,1)
                                      CALL dcopy(2*nkpt%ngwk*crge%n,c2(1,1,2),1,c2,1)
                                   ENDIF
                                ENDIF
                                bsclcs=1
                                CALL setbsstate
                                CALL updwf(c0(:,:,1),c2(:,:,1),sc0(:,:,1),tau0,fion,pme,gde,&
                                     vpp,eigv,rhoe,psi,nstate,&
                                     tfor,update_pot)
                                IF (paral%parent) THEN
                                   CALL wrener
                                   IF (paral%io_parent)&
                                        WRITE(6,'(2A)')&
                                        ' NFI      GEMAX       CNORM',&
                                        '           ETOT        DETOT      TCPU'
                                ENDIF
                             ELSE
                                CALL updwf(c0(:,:,1),c2(:,:,1),sc0(:,:,1),tau0,fion,pme,gde,&
                                     vpp,eigv,rhoe,psi,nstate,&
                                     tfor,update_pot)
                             ENDIF
                             cvbswf=ropt_mod%convwf
                          ENDIF
                          ropt_mod%convwf=cvhswf.AND.cvbswf
                       ENDIF
                    ENDIF
                    ropt_mod%modens=.FALSE.
                    ! Check for convergence in terms of energy change
                    CALL tol_chk_cnvener(ropt_mod%convwf,ener_com%etot,etoto)
                    ! Printout the evolution of the iterative optimization
                    IF (paral%parent) THEN
                       time22=m_walltime()
                       tcpu=(time22-time11)*0.001_real_8
                       CALL wrprint_wfopt(eigv(1,bsclcs),crge%f,ener_com%amu,nstate,&
                            ener_com%etot,etoto,tcpu,gemax,cnorm,dummy,.FALSE.,infw,-1)
                    ENDIF
                    CALL mp_bcast_byte(ener_com,size_in_bytes_of(ener_com),parai%source,parai%allgrp)
                    CALL mp_bcast(etoto,parai%source,parai%allgrp)
                    ! Periodically write the restart file
                    IF (MOD(ifcalc,store1%istore).EQ.0) THEN
                       CALL zhwwf(2,irec,c0,c2,nstate,eigv,&
                            tau0,tau0,tau0,iteropt%nfi)
                    ENDIF
                    ! End of wf optimization
                    IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)
                    IF (soft_com%exsoft.OR.ifcalc.EQ.cnti%nomore_iter) THEN
                       GOTO 200
                    ENDIF
                    IF (clc%classical) THEN
                       GOTO 101
                    ENDIF

                    IF (ropt_mod%convwf.AND.tfor) THEN
                       GOTO 101
                    ENDIF
                 ENDDO
101              CONTINUE
                 ropt_mod%convwf=.FALSE.
                 CALL vupdate(c0,psi,rhoe,convtest)
                 ropt_mod%sdiis=.TRUE.
                 ropt_mod%spcg=.TRUE.
                 tfor=.TRUE.
                 IF (convtest) THEN
                    GOTO 100
                 ENDIF
                 IF (cdftlog%reswf) THEN
                    IF (paral%io_parent)&
                         WRITE(6,*) "Re-initialising Wavefunction"
                    CALL rinitwf(c0,c2,sc0,crge%n,tau0,tau0,rhoe,psi)
                 ENDIF
              ENDDO
           ELSE
              update_pot=.TRUE.
              DO infw=1,cnti%nomore_iter
                 IF (infw.GT.1)update_pot=.FALSE.
                 time11=m_walltime()
                 ifcalc=ifcalc+1
                 ! Check if forces need to be calculated
                 CALL tol_chk_force(tfor,ropt_mod%calste,gemax,infw)
                 IF (clc%classical)tfor=.TRUE.
                 IF ((cntl%tdipd.AND.(wan05%loc_relocalize_in_scf.OR.infw==1)) &
                      .OR.vdwl%vdwd) THEN
                    IF (tfor.OR.hfxc3%twscr) THEN
                       CALL localize2(tau0,c0,c2,sc0,nstate)
                       vdwwfl%twannup=.TRUE.
                    ENDIF
                 ENDIF
                 ! Update the wavefunctions
                 IF (tkpts%tkpnt) THEN
                    IF (kdrho.LT.10._real_8*cntr%tolog) THEN
                       ropt_mod%calste = cntl%tprcp
                       tfor   = .TRUE.
                    ENDIF
                    CALL k_updwf(c0,c2,sc0,tau0,fion,pme,gde,vpp,eigv,&
                         rhoe,psi,nstate,tfor,infw)
                 ELSE
                    IF (.NOT. cntl%bsymm) THEN
                       CALL updwf(c0(:,:,1),c2(:,:,1),sc0(:,:,1),tau0,fion,pme,gde,vpp,eigv,&
                            rhoe,psi,nstate,tfor,update_pot)
                    ELSE
                       IF (.NOT.cvhswf) THEN
                          ! HS STATE
                          IF (infw.EQ.1) THEN
                             bsclcs=2
                             CALL setbsstate
                          ENDIF
                          CALL updwf(c0(:,:,2),c2(:,:,2),sc0(:,:,2),tau0,&
                               fion,pme,gde,vpp,eigv(1,2),rhoe,&
                               psi,nstate,tfor,update_pot)
                          cvhswf=ropt_mod%convwf
                       ELSEIF (.NOT.cvbswf) THEN
                          ! BS STATE
                          IF (bsclcs.EQ.2) THEN
                             IF (paral%parent) THEN
                                IF (paral%io_parent)&
                                     WRITE(6,*) 'BSYMM: HS STATE CONVERGED. ',&
                                     'SWITCHING TO BS STATE.'
                                etoto=0.0_real_8
                                IF ((infi.EQ.1).AND.&
                                     (.NOT.resthswf))THEN
                                   IF (paral%io_parent)&
                                        WRITE(6,'(/,T2,2A,/)')&
                                        '!!!! WARNING : COPYING C0, EIGV & ',&
                                        'C2 OF HS TO BS !!!!'
                                   CALL dcopy(2*nkpt%ngwk*crge%n,c0(1,1,2),1,c0,1)
                                   CALL dcopy(crge%n,eigv(1,2),1,eigv,1)
                                   CALL dcopy(2*nkpt%ngwk*crge%n,c2(1,1,2),1,c2,1)
                                ENDIF
                             ENDIF
                             bsclcs=1
                             CALL setbsstate
                             CALL updwf(c0(:,:,1),c2(:,:,1),sc0(:,:,1),tau0,fion,pme,gde,&
                                  vpp,eigv,rhoe,psi,nstate,&
                                  tfor,update_pot)
                             IF (paral%parent) THEN
                                CALL wrener
                                IF (paral%io_parent)&
                                     WRITE(6,'(2A)')&
                                     ' NFI      GEMAX       CNORM',&
                                     '           ETOT        DETOT      TCPU'
                             ENDIF
                          ELSE
                             CALL updwf(c0(:,:,1),c2(:,:,1),sc0(:,:,1),tau0,fion,pme,gde,&
                                  vpp,eigv,rhoe,psi,nstate,&
                                  tfor,update_pot)
                          ENDIF
                          cvbswf=ropt_mod%convwf
                       ENDIF
                       ropt_mod%convwf=cvhswf.AND.cvbswf
                    ENDIF
                 ENDIF
                 ropt_mod%modens=.FALSE.
                 ! Check for convergence in terms of energy change
                 CALL tol_chk_cnvener(ropt_mod%convwf,ener_com%etot,etoto)
                 ! Printout the evolution of the iterative optimization
                 IF (paral%parent) THEN
                    time22=m_walltime()
                    tcpu=(time22-time11)*0.001_real_8
                    CALL wrprint_wfopt(eigv(1,bsclcs),crge%f,ener_com%amu,nstate,&
                         ener_com%etot,etoto,tcpu,gemax,cnorm,dummy,.FALSE.,infw,-1)
                 ENDIF
                 CALL mp_bcast_byte(ener_com,size_in_bytes_of(ener_com),parai%source,parai%allgrp)
                 CALL mp_bcast(etoto,parai%source,parai%allgrp)
                 ! Periodically write the restart file
                 IF (qmmmrest.AND.(MOD(ifcalc,store1%isctore).EQ.0)) THEN
                    CALL zhwwf(2,irec,c0,c2,nstate,eigv,&
                         tau0,tau0,tau0,iteropt%nfi)
                 ENDIF
                 ! End of wf optimization
                 IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)
                 IF (soft_com%exsoft.OR.ifcalc.EQ.cnti%nomore_iter) THEN
                    GOTO 200
                 ENDIF
                 IF (clc%classical) THEN
                    GOTO 100
                 ENDIF

                 IF (ropt_mod%convwf.AND.tfor) THEN
                    GOTO 100
                 ENDIF
              ENDDO
           ENDIF
        ENDIF
        ! ==----------------------------------------------------------==
100     CONTINUE
        IF (cntl%bsymm) THEN
           ! BS energy calculation 
           ! BS state
           bsclcs=1
           CALL setbsstate
           CALL updwf(c0(:,:,1),c2(:,:,1),sc0(:,:,1),tau0,fion,pme,gde,vpp,eigv,&
                rhoe,psi,nstate,tfor,update_pot)
           IF (paral%parent) THEN
              IF (infi.EQ.cnti%nomore) THEN
                 IF (paral%io_parent) THEN
                    WRITE(6,*)
                    WRITE(6,*) 'BSYMM: BS STATE ENERGIES:'
                 END IF
                 CALL wrener
              ENDIF
              etotbs=ener_com%etot
              CALL dcopy(3*maxsys%nax*maxsys%nsx,fion,1,fnbs,1)
           ENDIF
           ! HS state      
           bsclcs=2
           CALL setbsstate
           CALL updwf(c0(:,:,2),c2(:,:,2),sc0(:,:,2),tau0,fion,pme,&
                gde,vpp,eigv(1,2),rhoe,psi,&
                nstate,tfor,update_pot)
           IF (paral%parent) THEN
              IF (infi.EQ.cnti%nomore) THEN
                 IF (paral%io_parent) THEN
                    WRITE(6,*)
                    WRITE(6,*) 'BSYMM: HS STATE ENERGIES:'
                 END IF
                 CALL wrener
              ENDIF
              etoths=ener_com%etot
              ! LS state      
              ener_com%etot = (1 + cnstwgt) * etotbs - cnstwgt * etoths
              CALL wrccfl(ener_com%etot,etotbs,etoths)
              CALL lsforce(fnbs,fion)
           ENDIF
        ENDIF
        ! if(cntl%prfo) then       ! use CONVERGENCE CALFOR 1.0 instead
        ! CALL FORCEDR(C0,C2,SC0,RHOE,PSI,TAU0,FION,EIGV,
        ! &                           N,1,.TRUE.,.TRUE.)
        ! endif
        IF (lqmmm%qmmm)THEN
           CALL updwf(c0(:,:,1),c2(:,:,1),sc0(:,:,1),tau0,fion,pme,gde,vpp,eigv,&
                rhoe,psi,nstate,tfor,.TRUE.)
        ENDIF
        IF (cntl%tddft) THEN
           CALL lr_tddft(c0(:,:,1),c1,c2(:,:,1),sc0,rhoe,psi,tau0,fion,eigv,&
                nstate,.TRUE.,td01%ioutput)
        ENDIF
        IF (qmmmrest.AND.&
             (store1%isctore.EQ.0 .OR. store1%isctore.GT.cnti%nomore_iter) ) THEN
           CALL zhwwf(2,irec,c0,c2,nstate,eigv,tau0,tau0,tau0,iteropt%nfi)
        ENDIF
        ! ==----------------------------------------------------------==
        ! Update stress tensor.
        IF (ropt_mod%calste) THEN
           CALL totstr
        ENDIF
        IF (paral%parent) THEN
           detot=ener_com%etot-etotg
           IF (infi.NE.1) THEN
              ! Decrease step for very large drops in energy..AA
              IF (detot.GT.1.e-6_real_8.OR.detot.LE.-5.e-3_real_8) THEN
                 CALL dscal(cotc0%nodim,0.5_real_8,dtm(1),1)
              ELSE
                 ! CALL DSCAL(NODIM,1.25_real_8,DTM(1),1)
                 ! Use a more conservative upward scaling factor..AA
                 CALL dscal(cotc0%nodim,1.10_real_8,dtm(1),1)
              ENDIF
           ENDIF
           CALL dscal(3*maxsys%nax*maxsys%nsx,-1.0_real_8,fion(1,1,1),1)
           ! UPDATE NUCLEAR HESSIAN
           CALL puttau(fion,dxpar)

           IF (cntl%rfo) THEN
              ! ..USE POWELL UPDATE FOR TRANSITION STATE SEARCH
              IF (infi.NE.1) THEN
                 IF (cnti%nsorder.EQ.0) THEN
                    CALL hessup(xpar,ypar,dxpar,dypar)
                 ELSE
                    CALL hespow(xpar,ypar,dxpar,dypar)
                 ENDIF
              ENDIF
              CALL cnstfc(tau0,tscr,dxpar,cnmax,.TRUE.)
           ELSE
              CALL cnstfc(tau0,tscr,dxpar,cnmax,.FALSE.)
              IF (.NOT.(cntl%tsdp.OR.cntl%lbfgs.OR.cntl%prfo).AND.infi.NE.1.AND.&
                   (.NOT.cntl%tinr.OR.inr_logical%tmixgdiis)) THEN
                 CALL hessup(xpar,ypar,dxpar,dypar)
              ENDIF
           ENDIF
           IF (.NOT.(cntl%tsdp.OR.cntl%lbfgs.OR.cntl%prfo)) THEN
              ! Symmetrize Hessian
              IF ( symmi%indpg.NE.0.AND.(symmi%nrot.NE.1.OR.symmi%ntvec.GT.1) ) THEN
                 IF (cotc0%nodim.EQ.3*ions1%nat) THEN
                    CALL symmat(hess,0)
                 ENDIF
              ENDIF
              ! AK: hessian for typical QM/MM is huge. don't write it.
              IF (.NOT.cntl%tqmmm) CALL hessout
           ENDIF
           CALL dcopy(cotc0%nodim, xpar(1),1, ypar(1),1)
           CALL dcopy(cotc0%nodim,dxpar(1),1,dypar(1),1)
           ! CHECK GRADIENT
           CALL gnodim(dxpar,cotc0%nodim,gnmax,gnorm)
           IF (cntl%tprcp) THEN
              ! Check cell gradient (stress)
              CALL dstre(metr_com%htfor,dstress,dsmax)
              ropt_mod%convge=gnmax.LT.cntr%tolng.AND.dstress.LT.cntr%tolcg
           ELSE
              IF (gnmax.LT.cntr%tolng) ropt_mod%convge=.TRUE.
           ENDIF
           ! Printout the evolution of the iterative optimization
           time2=m_walltime()
           tcpu=(time2-time1)*0.001_real_8
           CALL wrprint_geopt(eigv,crge%f,ener_com%amu,nstate,tau0,fion,&
                ener_com%etot,etotg,tcpu,gemax,gnmax,dxmax,gnorm,cnmax,&
                dstress,dsmax,ropt_mod%convge,ifcalc,infi)
           ropt_mod%convge = (ropt_mod%convge.AND.lcvmsk)
        ENDIF  ! IF (PARENT)
        ! ==----------------------------------------------------------==
        ! CHECK IF IONS CONVERGED
        IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)
        CALL testex(ropt_mod%convge)
        IF (ropt_mod%convge.OR.soft_com%exsoft) GOTO 200
        IF (ifcalc.GE.cnti%nomore_iter) THEN
           IF (paral%parent) THEN
              IF (paral%io_parent)&
                   WRITE(6,'(1X,64("!"))')
              IF (paral%io_parent)&
                   WRITE(6,'(" !!",A,A,I5,A,T64,"!!")')&
                   ' RGMOPT| THE MAXIMUM NUMBER OF ITERATIONS IS REACHED',&
                   ' (',ifcalc,')'
              IF (paral%io_parent)&
                   WRITE(6,'(1X,64("!"))')
           ENDIF
           GOTO 200
        ENDIF
        IF (infi.EQ.cnti%nomore) THEN
           IF (paral%parent) THEN
              IF (paral%io_parent)&
                   WRITE(6,'(1X,64("!"))')
              IF (paral%io_parent)&
                   WRITE(6,'(" !!",A,A,I5,A,T64,"!!")')&
                   ' RGMOPT| THE MAX. NUMBER OF GEOSTEPS IS REACHED',&
                   ' (',infi,')'
              IF (paral%io_parent)&
                   WRITE(6,'(1X,64("!"))')
           ENDIF
           GOTO 200
        ENDIF
        ! ==----------------------------------------------------------==
        IF ( qmmmrest.AND.(.NOT.(MOD(ifcalc,store1%isctore).EQ.0)).AND.&
             MOD(iteropt%nfi,store1%istore).EQ.0) THEN
           wcmpsv = cnti%wcompb
           cnti%wcompb = MIN(cnti%wcompb*4, 8)
           CALL zhwwf(1,irec,c0,c2,nstate,eigv,tau0,tau0,tau0,iteropt%nfi)
           cnti%wcompb = INT(wcmpsv)
        ENDIF
        ! ==----------------------------------------------------------==
        ! UPDATE THE NUCLEAR POSITIONS
        CALL dcopy(3*maxsys%nax*maxsys%nsx,tau0(1,1,1),1,taup(1,1,1),1)
        IF (cntl%tinr) THEN
           CALL puttau(tau0,xpar)
           CALL inr_dr(c0,c11,nstate,xpar,dxpar,psi,rhoe,&
                tras,tmp,eirop,eivps,told,gold)
           IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)
           IF (soft_com%exsoft) THEN
              GOTO 200
           ENDIF
        ENDIF
        ! SAVE THE WAVEFUNCTION (ONLY AT INIT)
        IF (cntl%lbfgs.OR.cntl%prfo) THEN
           CALL stack_do_sched(ielstk,c0,2*nkpt%ngwk*nstate)
        ENDIF
        IF (paral%parent) THEN
           IF (cntl%tsdp) THEN
              CALL sdion(xpar,dxpar)
           ELSEIF (cntl%gdiis) THEN
              CALL rgdiis(xpar,dxpar,told,gold)
           ELSEIF (cntl%bfgs) THEN
              CALL rbfgs(xpar,dxpar)
           ELSEIF (cntl%lbfgs) THEN
              lret = .FALSE.
              CALL rlbfgs(xpar,dxpar,cotc0%nodim,ener_com%etot+ener_com%ecnstr+ener_com%erestr,lret)
           ELSEIF (cntl%prfo) THEN
              CALL rprfo(xpar,dxpar,ener_com%etot,&
                   ltolad,genvmx,lcvmsk,lpexit)
              IF (lpexit) THEN
                 ropt_mod%convge = .TRUE.
              ENDIF
           ELSEIF (cntl%rfo) THEN
              CALL rrfo(xpar,dxpar)
           ENDIF
           CALL gettau(tau0,xpar)
           dxmax=0.0_real_8
           ! Check convergence on atomic positions (for the next it.)
           ! (DXMAX .LT. 0.01*TOLNG)
           DO i=1,cotc0%nodim
              dx=ABS(xpar(i)-ypar(i))
              dxmax=MAX(dx,dxmax)
           ENDDO
           IF (dxmax.LT.1.e-2_real_8*cntr%tolng.AND..NOT.(cntl%tsdp.OR.cntl%prfo))&
                ropt_mod%convge=.TRUE.
           CALL geofile(tau0,fion,'WRITE')
        ENDIF              ! IF (PARENT)
        ! CONVERGENCE CRITERIA FOR THE WAVEFUNCTION GRADIENT
        CALL tol_det_grad(ltolad,gnmax,genvmx)
        ! CONVERGENCE CRITERIA BASED ON ENERGY CHANGE
        CALL tol_det_ener(ltolad,infi,detot,deaccp)
        ! IF SCHEDULED, STORE OR RESTORE THE WAVEFUNCTION NOW
        IF (cntl%lbfgs.OR.cntl%prfo) THEN
           CALL stack_do_sched(ielstk,c0,2*nkpt%ngwk*nstate)
        ENDIF
        ! MOVE THE CELL
        IF (cntl%tprcp) THEN
           CALL totstr
           IF (cntl%tsdc) THEN
              CALL sdcell(tau0)
           ENDIF
        ENDIF
        ! Broadcast new nuclear positions
        CALL mp_bcast(tau0,SIZE(tau0),parai%source,parai%allgrp)
        ! UPDATE SYMMETRIES (POINT GROUP)
        CALL updatsym(tau0)
        IF (cntl%tprcp) THEN
           CALL mp_bcast_byte(metr_com, size_in_bytes_of(metr_com),parai%source,parai%allgrp)
           DO i=1,3
              parm%a1(i) = metr_com%ht(1,i)
              parm%a2(i) = metr_com%ht(2,i)
              parm%a3(i) = metr_com%ht(3,i)
           ENDDO
           CALL ihmat(metr_com%ht,metr_com%htm1,parm%omega)
           CALL newcell
        ENDIF
        ! Reset swap file to update (option BLOCK and ALL)
        IF (tkpts%tkblock) THEN
           CALL reskpt_swap
        ENDIF
        ! Initialize form factors
        ! Reorthogonalize wavefunctions
        CALL mm_dim(mm_go_qm,statusdummy)
        CALL phfac(tau0)
        IF (cntl%tdiag) THEN
           IF (tmovr) THEN
              CALL dcopy(nnx,rin0,1,rhoe,1)
              CALL moverho(rhoe,psi)
              CALL dcopy(nnx,rhoe,1,rin0,1)
           ENDIF
           ! Update density
           IF (.NOT.cntl%bsymm) THEN
              CALL extrap(nnx,andr2%alxmix,rm1,rin0,rinp)
              CALL dcopy(nnx,rin0,1, rm1,1)
              CALL dcopy(nnx,rinp,1,rin0,1)
           ENDIF
        ENDIF
        IF (fint1%ttrot) THEN
           CALL calc_alm
        ENDIF
        IF (corel%tinlc) THEN
           CALL copot(rhoe,psi,cntl%tsdc)
        ENDIF
        IF (.NOT.cntl%nonort) THEN
           IF (.NOT.cntl%bsymm)THEN
              IF (pslo_com%tivan) THEN
                 CALL rnlsm(c0(:,:,1),nstate,1,1,.FALSE.)
              ENDIF
              CALL ortho(nstate,c0(:,:,1),c2)
           ELSE
              bsclcs=1
              CALL setbsstate
              IF (pslo_com%tivan) THEN
                 CALL rnlsm(c0(:,:,1),nstate,1,1,.FALSE.)
              ENDIF
              CALL ortho(nstate,c0(:,:,1),c2)
              bsclcs=2
              CALL setbsstate
              IF (pslo_com%tivan)CALL rnlsm(c0(:,:,2),nstate,&
                   1,1,.FALSE.)
              CALL ortho(nstate,c0(:,:,2),c2(1,1,2))
           ENDIF
        ENDIF
        CALL mm_dim(mm_revert,statusdummy)
     ENDDO
     ! ==================================================================
     ! ==     END OF MAIN LOOP                                         ==
     ! ==================================================================
  ENDIF
200 CONTINUE
  IF (cntl%cdft)THEN
     IF (cntl%cdft_weight)CALL write_w(wdiff,"FINAL")
  ENDIF
  IF (paral%parent) THEN
     IF (paral%io_parent)&
          WRITE(6,'(1X,64("="))')
     IF (paral%io_parent)&
          WRITE(6,'(1X,"=",T17,A,T65,"=")')&
          'END OF GEOMETRY OPTIMIZATION'
     IF (paral%io_parent)&
          WRITE(6,'(1X,64("="),/,/)')
  ENDIF
  ! Store wavefunction and positions
  IF (qmmmrest)&
       CALL zhwwf(2,irec,c0,c2,nstate,eigv,tau0,tau0,tau0,iteropt%nfi)
  IF (rout1%rhoout.AND..NOT.cntl%bsymm) THEN
     CALL rhopri(c0,tau0,rhoe,psi(:,1),nstate,nkpt%nkpnt)
  ENDIF
  IF (rout1%rhoout.AND.cntl%bsymm) THEN
     bsclcs=1
     CALL setbsstate
     CALL rhopri(c0,tau0,rhoe,psi(:,1),crge%n,nkpt%nkpnt)
     bsclcs=2
     CALL setbsstate
     CALL rhopri(c0(:,:,2:),tau0,rhoe,psi(:,1),crge%n,nkpt%nkpnt)
     bsclcs=1
     CALL setbsstate
  ENDIF
  IF (.NOT.ropt_mod%calste.AND.cntl%tpres) THEN
     ! Calculate stress tensor
     IF (cntl%tdiag) THEN
        CALL zeroing(fion)!,3*maxsys%nax*maxsys%nsx)
        CALL updrho(c0,c2,cm,sc0,cm(nx),vpp,tau0,fion,eigv,&
             rhoe,psi,&
             crge%n,.TRUE.,.FALSE.,cntl%tpres,infi,thl,nhpsi)
     ELSE
        ropt_mod%calste=cntl%tpres
        CALL updwf(c0(:,:,1),c2(:,:,1),sc0(:,:,1),tau0,fion,pme,gde,vpp,eigv,&
             rhoe,psi,crge%n,.TRUE.,.TRUE.)
     ENDIF
     CALL totstr
  ENDIF
  IF (paral%parent) THEN
     CALL dumpr
     CALL cnstpr
     ! (T.D.) in order to have same GNORM in results and good FION
     ! CALL TAUCL(FION)
     ! CALL GSIZE(FION,GNMAX,GNORM)
     CALL puttau(fion,dxpar)
     CALL gnodim(dxpar,cotc0%nodim,gnmax,gnorm)
     ! Give final results
     IF (.NOT.cntl%bsymm) THEN
        CALL finalp(tau0,fion,velp,eigv)
     ELSE
        CALL wrgeof(tau0,fion)
        IF (paral%io_parent)&
             WRITE(6,'(A)') ' NUCLEAR GRADIENT: '
        IF (paral%io_parent)&
             WRITE(6,'(2(A,1PE15.5))') '    MAX. COMPONENT =',&
             gnmax,'         NORM =',gnorm
     ENDIF
     CALL geofile(tau0,fion,'WRITE')
     ! Do vibrational analysis
     IF (cntl%prfo.AND.nsvib.NE.0.AND.lvlhes) THEN
        CALL pvibana
     ENDIF
  ENDIF
  ! Deallocation
  IF (cntl%tinr.OR.paral%parent) THEN
     DEALLOCATE(xpar,STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
          __LINE__,__FILE__)
     DEALLOCATE(dxpar,STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
          __LINE__,__FILE__)
     DEALLOCATE(ypar,STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
          __LINE__,__FILE__)
     DEALLOCATE(dypar,STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
          __LINE__,__FILE__)
  ENDIF
  IF (paral%parent) THEN
     DEALLOCATE(dtm,STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
          __LINE__,__FILE__)
     IF (cntl%gdiis.OR.cntl%tinr) THEN
        DEALLOCATE(gold,STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
             __LINE__,__FILE__)
        DEALLOCATE(told,STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
             __LINE__,__FILE__)
     ENDIF
     IF (cntl%lbfgs) THEN
        DEALLOCATE(wlbfgs,STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
             __LINE__,__FILE__)
     ENDIF
     IF (cntl%prfo)THEN
        DEALLOCATE(hesscr,STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
             __LINE__,__FILE__)
        IF (lrlbfgs.GT.0) THEN
           DEALLOCATE(wlbfgs,STAT=ierr)
           IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
                __LINE__,__FILE__)
        ENDIF
     ENDIF
     IF (cntl%lbfgs.OR.cntl%prfo) THEN
        CALL stack_destroy(ielstk)
     ENDIF
  ENDIF
  DEALLOCATE(fion,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
       __LINE__,__FILE__)
  DEALLOCATE(tscr,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
       __LINE__,__FILE__)
  DEALLOCATE(taup,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
       __LINE__,__FILE__)
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
     DEALLOCATE(rinp,STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
          __LINE__,__FILE__)
     DEALLOCATE(rm1,STAT=ierr)
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
  DEALLOCATE(psi,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
       __LINE__,__FILE__)
  DEALLOCATE(scr,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
       __LINE__,__FILE__)
  IF (cntl%bsymm) DEALLOCATE(fnbs,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
       __LINE__,__FILE__)
  IF (cntl%tddft) THEN
     DEALLOCATE(potr,STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
          __LINE__,__FILE__)
     DEALLOCATE(rhoo,STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
          __LINE__,__FILE__)
  ENDIF
  ! ==--------------------------------------------------------------==
  IF (cntl%tinr) THEN
     DEALLOCATE(rho0,STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
          __LINE__,__FILE__)
     DEALLOCATE(drhoe,STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
          __LINE__,__FILE__)

     ! potentials and related:
     DEALLOCATE(VofRho0,STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
          __LINE__,__FILE__)
     DEALLOCATE(eirop,STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
          __LINE__,__FILE__)
     DEALLOCATE(eivps,STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
          __LINE__,__FILE__)

     ! others & scratch:
     DEALLOCATE(z11,STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
          __LINE__,__FILE__)
     DEALLOCATE(rs_v,STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
          __LINE__,__FILE__)
     DEALLOCATE(tmp,STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
          __LINE__,__FILE__)
     DEALLOCATE(fnl00,STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
          __LINE__,__FILE__)
     DEALLOCATE(dfnl00,STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
          __LINE__,__FILE__)
     DEALLOCATE(ddfnl,STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
          __LINE__,__FILE__)
     DEALLOCATE(sd0,STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
          __LINE__,__FILE__)
     DEALLOCATE(tras,STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
          __LINE__,__FILE__)
     DEALLOCATE(c11,STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
          __LINE__,__FILE__)
     DEALLOCATE(direct,STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
          __LINE__,__FILE__)
     DEALLOCATE(correct,STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
          __LINE__,__FILE__)
     DEALLOCATE(zprec,STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
          __LINE__,__FILE__)
  ENDIF

  ! ==--------------------------------------------------------------==
  CALL tihalt(procedureN,isub)
  RETURN
END SUBROUTINE rgmopt
! ==================================================================
SUBROUTINE give_scr_rgmopt(lrgmopt,tag)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE newcell_utils, ONLY : newcell, give_scr_newcell
  USE parac, ONLY : paral,parai
  USE elct , ONLY:crge
  USE nlcc , ONLY:corel
  USE store_types , ONLY:restart1,rout1
  USE pslo , ONLY:pslo_com
  USE fint , ONLY:fint1
  USE atwf , ONLY:tmovr
  USE system , ONLY:cntl
  USE vdwcmod , ONLY:vdwl
  USE rgdiis_utils, ONLY : rgdiis, give_scr_rgdiis
  USE moverho_utils, ONLY : moverho,give_scr_moverho
  USE rbfgs_utils, ONLY : rbfgs, give_scr_rbfgs
  USE rlbfgs_utils, ONLY : give_scr_rlbfgs, give_scr_rprfo
  USE rrfo_utils, ONLY : rrfo, give_scr_rrfo
  USE sdion_utils, ONLY : sdion, give_scr_sdion
  USE copot_utils, ONLY : give_scr_copot, copot
  USE ddipo_utils, ONLY : setdip, give_scr_ddipo
  USE forces_diag_utils, ONLY : forces_diag,updiag,dens_for_diag,give_scr_forces_diag
  USE calc_alm_utils, ONLY : calc_alm,calc_k_alm,give_scr_calc_alm
  USE inr_dr_utils, ONLY : give_scr_inr,inr_dr
  USE lr_tddft_utils, ONLY : give_scr_lr_tddft, lr_tddft
  USE updwf_utils, ONLY : give_scr_updwf, updwf
  USE ortho_utils, ONLY : ortho,give_scr_ortho
  USE rhopri_utils, ONLY : rhopri ,give_scr_rhopri
  USE rhoofr_utils, ONLY : give_scr_rhoofr,rhoofr
  USE initrun_utils, ONLY : give_scr_initrun
  USE rnlsm_utils, ONLY : rnlsm, give_scr_rnlsm
  USE forcedr_utils, ONLY : give_scr_forcedr
  IMPLICIT NONE
  INTEGER                                    :: lrgmopt
  CHARACTER(len=30)                          :: tag

  INTEGER :: lcalc_alm, lcopot, ldeort, ldipd, lforces, lforces_diag, &
      linitrun, lmoverho, lnewcell, lortho, lposupa, lprepv, lrbfgs, lrgdiis, &
      lrhoofr, lrhopri, lrinr, lrlbfgs, lrnlsm, lrortv, lrprfo, lrrfo, &
      lsdion, ltddft, lupdwf, nstate

  nstate=crge%n
  lcopot=0
  lortho=0
  lrnlsm=0
  lrhoofr=0
  lcalc_alm=0
  lforces_diag=0
  lforces=0
  lsdion=0
  lrgdiis=0
  lrbfgs=0
  lrlbfgs=0
  lrprfo=0
  lrrfo=0
  lrinr=0
  lnewcell=0
  ldeort=0
  lupdwf=0
  lprepv=0
  lposupa=0
  lrortv=0
  lforces_diag=0
  lrhopri=0
  lmoverho=0
  ltddft=0
  CALL give_scr_initrun(linitrun,tag)
  IF (restart1%restart) THEN
     IF (corel%tinlc) THEN
        CALL give_scr_copot(lcopot,tag)
     ENDIF
     CALL give_scr_ortho(lortho,tag,nstate)
  ENDIF
  IF (cntl%tdiag) THEN
     IF (pslo_com%tivan) THEN
        CALL give_scr_rnlsm(lrnlsm,tag,nstate,.FALSE.)
     ENDIF
     CALL give_scr_rhoofr(lrhoofr,tag)
     IF (fint1%ttrot) THEN
        CALL give_scr_calc_alm(lcalc_alm,tag)
     ENDIF
     CALL give_scr_forces_diag(lforces_diag,tag,nstate,.TRUE.)
  ENDIF
  IF (cntl%tsdp.AND.cntl%tsde) THEN
     CALL give_scr_forcedr(lforces,tag,nstate,.TRUE.,.TRUE.)
     IF (paral%parent) THEN
        CALL give_scr_sdion(lsdion,tag)
     ENDIF
     IF (cntl%tprcp) THEN
        CALL give_scr_newcell(lnewcell,tag)
     ENDIF
     IF (corel%tinlc) THEN
        CALL give_scr_copot(lcopot,tag)
     ENDIF
     IF (.NOT.cntl%nonort) THEN
        CALL give_scr_ortho(lortho,tag,nstate)
     ENDIF
  ELSE
     IF (cntl%tprcp) THEN
        CALL give_scr_newcell(lnewcell,tag)
     ENDIF
     IF (cntl%tdiag) THEN
        CALL give_scr_forces_diag(lforces_diag,tag,nstate,.TRUE.)
        IF (tmovr) THEN
           CALL give_scr_moverho(lmoverho,tag)
        ENDIF
     ELSE
        CALL give_scr_updwf(lupdwf,tag,nstate,.TRUE.)
     ENDIF
     IF (paral%parent) THEN
        ! HESSUP needs 3*NODIM which is small than the other routines
        IF (cntl%tsdp) THEN
           CALL give_scr_sdion(lsdion,tag)
        ELSEIF (cntl%gdiis) THEN
           CALL give_scr_rgdiis(lrgdiis,tag)
        ELSEIF (cntl%bfgs) THEN
           CALL give_scr_rbfgs(lrbfgs,tag)
        ELSEIF (cntl%lbfgs) THEN
           CALL give_scr_rlbfgs(lrlbfgs,tag)
        ELSEIF (cntl%prfo) THEN
           CALL give_scr_rprfo(lrprfo,tag)
        ELSEIF (cntl%rfo) THEN
           CALL give_scr_rrfo(lrrfo,tag)
        ELSEIF (cntl%tinr) THEN
           CALL give_scr_inr(lrinr,tag,nstate)
        ENDIF
     ENDIF
     IF (fint1%ttrot) THEN
        CALL give_scr_calc_alm(lcalc_alm,tag)
     ENDIF
     IF (.NOT.cntl%nonort) THEN
        CALL give_scr_ortho(lortho,tag,nstate)
     ENDIF
  ENDIF
  IF (cntl%tddft) THEN
     CALL give_scr_lr_tddft(ltddft,.TRUE.,tag)
  ENDIF
  IF (rout1%rhoout) THEN
     CALL give_scr_rhopri(lrhopri,tag,nstate)
  ENDIF
  IF (cntl%tdipd.OR.vdwl%vdwd) THEN
     CALL give_scr_ddipo(ldipd,tag)
  ELSE
     ldipd=0
  ENDIF
  lrgmopt=MAX(linitrun,lcopot,lortho,lrnlsm,lrhoofr,&
       lcalc_alm,lforces_diag,lforces,&
       lsdion,lrgdiis,lrbfgs,lrrfo,lrinr,lrlbfgs,lrprfo,&
       lnewcell,ldeort,lupdwf,lprepv,lposupa,lrortv,&
       ldipd,lforces_diag,lrhopri,lmoverho,ltddft)
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE give_scr_rgmopt
! ==================================================================

