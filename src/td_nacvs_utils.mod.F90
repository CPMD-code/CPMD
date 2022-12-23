MODULE td_nacvs_utils
  USE canon_utils,                     ONLY: canon
  USE cell,                            ONLY: cell_com
  USE cnst,                            ONLY: factem,&
                                             fbohr
  USE constr_utils,                    ONLY: vecprod
  USE coor,                            ONLY: tau0
  USE cplngs_utils,                    ONLY: cplconfig,&
                                             get_eind,&
                                             give_scr_get_eind,&
                                             prcplngs
  USE cplngsmod,                       ONLY: &
       csurf, f_high, f_low, f_med, fc_high, fc_low, fc_med, nsurf, nvect, &
       tallat, talldof, tolcpl, tspecv
  USE cppt,                            ONLY: gk,&
                                             inyh,&
                                             nzh,&
                                             rhops,&
                                             scg,&
                                             vps
  USE dotp_utils,                      ONLY: dotp
  USE eicalc_utils,                    ONLY: eicalc,&
                                             eicalc1
  USE eind_ii_utils,                   ONLY: eind_ii
  USE eind_loc_utils,                  ONLY: eind_loc
  USE eind_nl_utils,                   ONLY: eind_nl
  USE ekinpp_utils,                    ONLY: ekinpp
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE fftmain_utils,                   ONLY: fwfftn
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_def
  USE fnlalloc_utils,                  ONLY: fnl_set,&
                                             fnlalloc
  USE forcedr_driver,                  ONLY: forcedr
  USE forcedr_utils,                   ONLY: give_scr_forcedr
  USE forcep_utils,                    ONLY: rhoe_psi_size
  USE geq0mod,                         ONLY: geq0
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE ksdiag_utils,                    ONLY: ksdiag
  USE linres,                          ONLY: lr01,&
                                             lr02,&
                                             lr03,&
                                             nolr,&
                                             tshl,&
                                             xfmqc
  USE lr_tddft_drhoe,                  ONLY: print_ex_rho,&
                                             relax_rho
  USE lr_xcpot_utils,                  ONLY: lr_xcpot
  USE mp_interface,                    ONLY: mp_sum
  USE nl_res_utils,                    ONLY: give_scr_nl_res,&
                                             nl_res
  USE nlps,                            ONLY: imagp,&
                                             ndfnl,&
                                             nghtol,&
                                             nlm,&
                                             nlps_com,&
                                             wsg
  USE nose,                            ONLY: glib
  USE opt_lr_utils,                    ONLY: give_scr_opt_lr,&
                                             opt_lr
  USE parac,                           ONLY: parai,&
                                             paral
  USE poin,                            ONLY: potr,&
                                             rhoo
  USE pslo,                            ONLY: pslo_com
  USE rho1ofr_utils,                   ONLY: rho1ofr,&
                                             rhoabofr
  USE rhoofr_utils,                    ONLY: give_scr_rhoofr,&
                                             rhoofr
  USE rmas,                            ONLY: rmass
  USE rnlsm1_utils,                    ONLY: rnlsm1
  USE rnlsm2_utils,                    ONLY: rnlsm2
  USE rnlsm_2d_utils,                  ONLY: give_scr_rnlsm_2d,&
                                             rnlsm_2d
  USE rnlsm_utils,                     ONLY: give_scr_rnlsm,&
                                             rnlsm
  USE ropt,                            ONLY: iteropt,&
                                             ropt_mod
  USE sfac,                            ONLY: ddfnl,&
                                             dfnl,&
                                             ei1,&
                                             ei2,&
                                             ei3,&
                                             eigrb,&
                                             fnl
  USE sgpp,                            ONLY: sgpp1,&
                                             sgpp2
  USE spin,                            ONLY: clsd,&
                                             lspin2,&
                                             spin_mod,&
                                             tdsp1
  USE symm,                            ONLY: symmi
  USE symtrz_utils,                    ONLY: give_scr_symvec
  USE system,                          ONLY: &
       cnti, cntl, fpar, ipept, maxsys, ncpw, nkpt, parap, parm, spar
  USE td_force_utils,                  ONLY: td_force
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE v1ofrho1_utils,                  ONLY: give_scr_v1ofrho1,&
                                             v1ofrho1
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: vcouplings
  PUBLIC :: vcouplings_lr
  PUBLIC :: give_scr_vcoupling

CONTAINS

  ! ==================================================================
  SUBROUTINE vcouplings(c0,c1,c2,sc0,rhoe,eigv,eigo,&
       n,no,nvir,nact,nroot)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(:,:), c2(:,:), sc0(*)
    REAL(real_8), TARGET                     :: rhoe(:,:)
    REAL(real_8)                             :: eigv(:), eigo(*)
    INTEGER                                  :: n, no, nvir, nact
    COMPLEX(real_8)                          :: c1(ncpw%ngw,nact,*)
    INTEGER                                  :: nroot

    CHARACTER(*), PARAMETER                  :: procedureN = 'vcouplings'

    INTEGER                                  :: ia, iaa, ic, ierr, iroot, is, &
                                                j, k, nc, ndim
    LOGICAL                                  :: ttrinv
    REAL(real_8), ALLOCATABLE                :: fion_na(:,:,:), sder(:,:), &
                                                vcloca(:,:,:,:), &
                                                vclocb(:,:,:,:), vibe(:)

!(nnr1,clsd%nlsd)
! NAT
! NAT
! 
! Variables
! 

    ALLOCATE(vcloca(3,maxsys%nax,maxsys%nsx,100),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(fion_na(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    IF (cntl%tlsd) THEN
       ALLOCATE(vclocb(3,maxsys%nax,maxsys%nsx,100),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ELSE
       ALLOCATE(vclocb(3,maxsys%nax,maxsys%nsx,100),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ALLOCATE(vibe(3*ions1%nat*ions1%nsp*nroot*nroot),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(sder(3*ions1%nat,ions1%nsp*nroot*nroot),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! 
    ! ----------------------------------------------------------------
    ttrinv=.TRUE.
    ! ----------------------------------------------------------------
    CALL zeroing(vcloca)!,3*maxsys%nax*maxsys%nsx*100)
    IF (cntl%tlsd) CALL zeroing(vclocb)!,3*maxsys%nax*maxsys%nsx*100)
    ! ----------------------------------------------------------------
    CALL vcouplings_lr(c0,c1,c2,eigv,eigo,sc0,vcloca,vclocb,&
         nact,n,no,nroot,rhoe,fion_na)
    ! ----------------------------------------------------------------
    ! WRITE OUTPUT FILE FOR MOLEKEL
    ndim=3*ions1%nat
    CALL zeroing(sder)!,3*ions1%nat*ions1%nsp*nroot*nroot)
    CALL zeroing(vibe)!,3*ions1%nat*ions1%nsp*nroot*nroot)
    ! REFINE THE NACS
    IF (paral%io_parent)&
         WRITE(6,'(" ==",60("-"),"==")')
    IF (ttrinv)THEN
       CALL nacstrans(vcloca,vclocb,nroot)
    ENDIF
    ! J is the index that runs over the set of NACs.
    ! it starts with 4 for FORMAT reasons. Next set will be 5
    DO iroot=1,NINT((nroot*(nroot+1._real_8))/2._real_8)
       j=4+(iroot-1)
       nc=1
       ! NC is the number of non-adiabatic couplings pairs
       ! computed between couple of states
       DO ic=1,nc
          ia=0
          DO is=1,ions1%nsp
             DO iaa=1,ions0%na(is)
                DO k=1,3
                   sder(k+ia,j)=1.0*vcloca(k,iaa,is,iroot)
                ENDDO
                ia=ia+3
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    IF (tshl%txfmqc) THEN
       DO iroot=1,NINT((nroot*(nroot+1._real_8))/2._real_8)
          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                DO k=1,3
                   xfmqc%nacv(k,ia,is,iroot)=vcloca(k,ia,is,iroot)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    CALL writenacs(vibe,sder,ndim,tau0,nroot,0)
    IF (cntl%tlsd) THEN
       DO iroot=1,NINT((nroot*(nroot+1._real_8))/2._real_8)
          j=4+(iroot-1)
          nc=1
          ! NC is the number of non-adiabatic couplings pairs
          ! computed between couple of states
          DO ic=1,nc
             ia=0
             DO is=1,ions1%nsp
                DO iaa=1,ions0%na(is)
                   DO k=1,3
                      sder(k+ia,j)=1.0*vclocb(k,iaa,is,iroot)
                   ENDDO
                   ia=ia+3
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       CALL writenacs(vibe,sder,ndim,tau0,nroot,1)
    ENDIF
    ! 
    DEALLOCATE(sder,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(vibe,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(vcloca,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(vclocb,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(fion_na,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE vcouplings
  ! ==================================================================
  SUBROUTINE give_scr_vcoupling(lcplngs,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lcplngs
    CHARACTER(len=30)                        :: tag

    INTEGER :: lcplsub, lddxc_1d, lddxc_2d, ldfnl, ldone, leirop, leivps, &
      lforcedr, lget_eind, LH1NL, LNL_RES, LOPT_LR, LRHOE, lrhoofr, lrnlsm, &
      LRNLSM_2D, lsymvec, LV1OFRHO1, LVECT, NSTATE

! NHG
! N
! NLSD
! TCPLFD,TCPLLR,TALLDOF,NSURF
! NAT
! NDFNL
! NHG
! N
! NLSD
! TCPLFD,TCPLLR,TALLDOF,NSURF
! NAT
! NDFNL
! ==--------------------------------------------------------------==

    tag     = 'VCOUPLING'
    leivps = 2*ncpw%nhg
    leirop = 2*ncpw%nhg
    ldfnl  = 6*ions1%nat*maxsys%nhxs*ndfnl
    lh1nl  = 2*nkpt%ngwk*crge%n
    lddxc_1d = fpar%nnr1
    lddxc_2d = (2*clsd%nlsd-1)
    CALL rhoe_psi_size(lrhoe)
    lvect = 3*ions1%nat*nvect
    ldone = 3*ions1%nat/2
    lcplngs = 2*leivps + 2*leirop + ldfnl + lh1nl + 2*lrhoe +&
         LDDXC_1d*LDDXC_2d + NSURF + LVECT
    lcplngs = lcplngs + ldone + nsurf
    ! 
    nstate=crge%n
    CALL give_scr_forcedr (lforcedr,tag,nstate,.FALSE.,.FALSE.)
    CALL give_scr_rnlsm(lrnlsm,tag,nstate,.TRUE.)
    CALL give_scr_symvec(lsymvec,tag)
    lcplsub=MAX(lforcedr,lrnlsm,lsymvec)
    CALL give_scr_rhoofr(lrhoofr,tag)
    CALL give_scr_v1ofrho1(lv1ofrho1,tag)
    CALL give_scr_rnlsm_2d(lrnlsm_2d,tag,nstate)
    CALL give_scr_nl_res(lnl_res,nstate,tag)
    CALL give_scr_opt_lr(lopt_lr,"PHONON",tag)
    CALL give_scr_get_eind(lget_eind,tag)
    lcplsub=MAX(lcplsub,lrhoofr,lrnlsm_2d,lnl_res,lopt_lr,&
         LGET_EIND,LV1OFRHO1)

    lcplngs=lcplngs+lcplsub+100
    ! ==--------------------------------------------------------------==
  END SUBROUTINE give_scr_vcoupling
  ! ==================================================================
  SUBROUTINE vcouplings_lr(c0,clrtd,cr,eigv,eigo,sc0,vcloca,vclocb,&
       nact,nstate,no,nroot,rhoe,fion_na)
    ! ==--------------------------------------------------------------==
    ! ==        ANALYTIC CALCULATION OF INTERSURFACE COUPLINGS        ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(:,:), CR(:,:)
    REAL(real_8)                             :: EIGV(:), EIGO(*), &
                                                vcloca(:,:,:,:), &
                                                vclocb(:,:,:,:)
    INTEGER                                  :: nact
    COMPLEX(real_8)                          :: clrtd(ncpw%ngw,nact,*)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: SC0(nkpt%ngwk,nstate)
    INTEGER                                  :: no, nroot
    REAL(real_8), TARGET                     :: rhoe(:,:)
    REAL(real_8)                             :: fion_na(:,:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'vcouplings_lr'

    CHARACTER(len=10)                        :: TYPE
    CHARACTER(len=5)                         :: orbtyp
    COMPLEX(real_8), ALLOCATABLE             :: c1(:,:), c2(:,:), cscr(:), &
                                                eirop(:), eirop1(:), &
                                                eivps(:), eivps1(:), h1nl(:), &
                                                psi(:,:)
    INTEGER :: i, ia, iat, ic, icc, ici, id, ierr, ig, il_psi_1d, il_psi_2d, &
      imax, ir, is1, is2, isp, isub, iv, ivv, k, kmax, lddxc_1d, lddxc_2d, &
      ldfnl, ldone, lh1nl, loptib, lrhoe, lvect, n3, nap, natoms, nbp, ncscr, &
      nn, nolrb, nvpp
    INTEGER, ALLOCATABLE                     :: idone(:)
    LOGICAL                                  :: debug, dodirect, dolinres, &
                                                lconv, lprint, lprint2, &
                                                RELAXRHOE, tfdens
    REAL(real_8)                             :: dasum, ddot, deinv, e2(5), &
                                                eind, fc, fcomp, fcompa, &
                                                fcompb, fcsum, rlen, rmax, &
                                                tol_lr_bak, trace, vp
    REAL(real_8), ALLOCATABLE :: ddxc(:,:), deinvv(:), dfnlc(:,:,:,:,:,:), &
      eigs(:), eigt(:), fcsrf(:), fnlc(:,:,:,:,:), rhoov(:), tauc(:,:,:), &
      vect(:,:,:), vpp(:)

!(nnr1,clsd%nlsd)
! TAU0
! SPCG,IESR
! DDFNL
! NDFNL
! NLSD,NSUP
! F
! RHOO, POTR
! TAU0
! SPCG,IESR
! DDFNL
! NDFNL
! NLSD,NSUP
! F
! RHOO, POTR
! Variables
! TODO refactor these arrays
! 
! ==--------------------------------------------------------------==

    cntl%tsymrho=.FALSE.
    talldof=.TRUE.
    tallat=.TRUE.
    nvect=10
    symmi%nrot=1
    dodirect=.TRUE.
    IF (tshl%nacv_direct_only) THEN
       dolinres=.FALSE.
    ELSE
       dolinres=.TRUE.
    ENDIF
    relaxrhoe=.FALSE.
    lr01%lopti=3    ! cntl%diis
    ! LOPTI=2    ! cntl%pcg

    debug = .FALSE.
    IF (paral%parent.AND.debug) THEN
       IF (paral%io_parent)&
            WRITE(6,*) 'NACT,NSTATE,NO,NROOT'
       IF (paral%io_parent)&
            WRITE(6,*)  nact,nstate,no,nroot
    ENDIF
    ! ----------------------------------------------------------------==
    ! ------  WORKING AREA ---------------------------------------------
    ! C    cntl%vibrat=.TRUE.
    ! C    cntl%tsdin=.FALSE.
    ! C    cntl%tsdan=.TRUE.
    ! C    NVIB=1
    ! C    DLRS=0.1_real_8
    ! C    LOPTI=2    ! cntl%pcg
    ! ------  WORKING AREA ---------------------------------------------
    ! ----------------------------------------------------------------==
    IF (cntl%tlsd) THEN
       nap=spin_mod%nsup-tdsp1%nupel
       nbp=spin_mod%nsdown-tdsp1%ndoel
    ENDIF
    IF (cntl%tsde) THEN
       nvpp = 1
    ELSE IF (cntl%diis) THEN
       nvpp = ncpw%ngw
    ELSE IF (cntl%pcg) THEN
       nvpp = ncpw%ngw
    ENDIF
    ncscr=2*ncpw%ngw*MAX(crge%n,cnti%nkry_max*cnti%nkry_block)
    ! ==--------------------------------------------------------------==
    lprint  = .TRUE.          ! Progress report
    lprint2 = .FALSE.         ! Diagnostics
    lprint  = (lprint  .AND. paral%parent)
    lprint2 = (lprint2 .AND. paral%parent)
    CALL tiset ('   LRCPLNG',isub)

    ldfnl  = 6*ions1%nat*maxsys%nhxs*ndfnl
    lh1nl  = 2*nkpt%ngwk*nstate
    IF (lr03%txc_analytic) THEN
       lddxc_1d = fpar%nnr1
       lddxc_2d = (2*clsd%nlsd-1)
    ELSEIF (lr01%lopti.EQ.0 .OR. lr01%lopti.EQ.2) THEN
       lddxc_1d = fpar%nnr1
       lddxc_2d = (2*clsd%nlsd-1)
    ELSE
       lddxc_1d = 1
       lddxc_2d = 1
    ENDIF
    CALL rhoe_psi_size(lrhoe,il_psi_1d=il_psi_1d,il_psi_2d=il_psi_2d)
    IF (talldof) THEN
       lvect = 3*ions1%nat
    ELSE
       lvect = 3*ions1%nat*nvect
    ENDIF
    ldone = 3*ions1%nat/2+1

    ALLOCATE(ddfnl(ions1%nat,maxsys%nhxs,6,ndfnl),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__) ! TODO check this
    !   ALLOCATE(potr(lrhoe,1),STAT=ierr)
    !   IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
    !        __LINE__,__FILE__) ! TODO check this

    ALLOCATE(eivps(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(eirop(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(eivps1(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(eirop1(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(h1nl(lh1nl),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(rhoov(lrhoe),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(ddxc(lddxc_1d,lddxc_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(deinvv(nroot),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(vect(3, ions1%nat, nvect),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    IF (.NOT.talldof) THEN
       ALLOCATE(idone(3*ions1%nat),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(fcsrf(nroot),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
    ENDIF

    ALLOCATE(c1(nkpt%ngwk,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(c2(nkpt%ngwk,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(psi(il_psi_1d,il_psi_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(vpp(nvpp),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(cscr(ncscr),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(tauc(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(eigt(nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(eigs(nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(fnlc(imagp,ions1%nat,maxsys%nhxs,nstate,nkpt%nkpnt),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(dfnlc(imagp,ions1%nat,maxsys%nhxs,3,ndfnl,nkpt%nkpnt),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! IF (TSPECV) CALL CPLRDVEC (VECT, NAT, NVECT, TFDENS)
    ! WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
    nn=no
    ! WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
    ! 
    CALL zeroing(c2)!,ngw*nstate)
    ! ==--------------------------------------------------------------==
    CALL zeroing(eigt)!,nstate)
    ! EIGV ARE TAKEN FROM THE ARGUMENTS: cntl%tddft values 
    CALL forcedr(c0,c2,sc0,rhoe,psi,tau0,fion_na,eigt,&
         nn,1,.FALSE.,.TRUE.)
    CALL canon(c0,c2,crge%f,nn,eigt)
    orbtyp="CANON"
    ! recalculate FNL and DFNL for transformed orbitals.
    ! One could also rotate the old ones, but this is easier...
    CALL rnlsm(c0,nn,1,1,.TRUE.)
    !   CALL dcopy(nnr1*clsd%nlsd,rhoe,1,potr,1)
    ! RHOE must be computed using RELAXRHO for cntl%tddft
    CALL rhoofr(c0(:,1:nn),rhoe,psi(:,1),nn)
    CALL dcopy(fpar%nnr1*clsd%nlsd,rhoe,1,rhoo,1)

    ! ==--------------------------------------------------------------==
    ! == direct term (explicit derivative term)                        ==
    ! ==--------------------------------------------------------------==
    IF (dodirect) THEN
       ici=0
       DO ic = 0,nroot
          DO icc = ic+1, nroot
             ici=ici+1
             DO id=1,nn
                is1=id
                IF (cntl%tlsd) THEN
                   IF (is1.LE.tdsp1%nupel) THEN
                      is2=is1
                   ELSE
                      is2=nap+is1
                   ENDIF
                ELSE
                   is2=is1
                ENDIF
                CALL zeroing(rhoov(1:fpar%nnr1))!,nnr1)
                CALL zeroing(psi)!,maxfft*group%nogrp)
                IF (ic.EQ.0) THEN
                   CALL rhoabofr(1,c0(:,is1:is1),clrtd(:,is2:is2,icc),rhoov,psi(:,1))
                ELSE
                   CALL rhoabofr(1,clrtd(:,is1:is1,ic),clrtd(:,is2:is2,icc),rhoov,psi(:,1))
                ENDIF
                ! Transform density to G space
                !$omp parallel do private(IR)
                DO ir = 1,fpar%nnr1
                   psi(ir,1)=CMPLX(rhoov(ir),0._real_8,kind=real_8)
                ENDDO
                CALL  fwfftn(psi(:,1),.FALSE.,parai%allgrp)
                IF (cntl%tlsd) THEN
                   IF (is1.LE.tdsp1%nupel) THEN
                      CALL tdrhofc(vcloca(:,:,:,ici),rhoov,psi(:,1))
                   ELSE
                      CALL tdrhofc(vclocb(:,:,:,ici),rhoov,psi(:,1))
                   ENDIF
                ELSE
                   CALL tdrhofc(vcloca(:,:,:,ici),rhoov,psi(:,1))
                ENDIF
             ENDDO
          ENDDO
       ENDDO
       ! NACS contribution from the nonlocal PP
       ici=0
       DO ic=0,nroot
          DO icc = ic+1, nroot
             ici=ici+1
             CALL zeroing(fnl)
             CALL zeroing(dfnl)
             CALL rnlsm1(clrtd(1,1,icc),nn,1)
             CALL rnlsm2(clrtd(1,1,icc),nn,1,1)
             CALL dcopy(imagp*ions1%nat*maxsys%nhxs*nn*nkpt%nkpnt,fnl,1,fnlc,1)
             CALL dcopy(imagp*ions1%nat*maxsys%nhxs*3*ndfnl*nkpt%nkpnt,dfnl,1,dfnlc,1)
             IF (ic.EQ.0) THEN
                CALL rnlsm1(c0,nn,1)
                CALL rnlsm2(c0,nn,1,1)
             ELSE
                CALL rnlsm1(clrtd(1,1,ic),nn,1)
                CALL rnlsm2(clrtd(1,1,ic),nn,1,1)
             ENDIF
             CALL drhopp(vcloca,vclocb,fnlc,dfnlc,nn,nap,ici)
          ENDDO
       ENDDO
       ! 
       ici=0
       DO ic=0,nroot
          DO icc = ic+1, nroot
             ici=ici+1
             IF (ic.EQ.0) THEN
                deinv=1._real_8/eigo(icc)
             ELSE
                deinv=1._real_8/(eigo(icc)-eigo(ic))
             ENDIF
             CALL dscal(3*maxsys%nax*maxsys%nsx,-deinv,vcloca(1,1,1,ici),1)
             IF (cntl%tlsd) THEN
                CALL dscal(3*maxsys%nax*maxsys%nsx,-deinv,vclocb(1,1,1,ici),1)
             ENDIF
          ENDDO
       ENDDO

       CALL mp_sum(vcloca,3*maxsys%nax*maxsys%nsx*100,parai%allgrp)
       IF (cntl%tlsd) CALL mp_sum(vclocb,3*maxsys%nax*maxsys%nsx*100,parai%allgrp)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ==  Linear response for dependence of KS operator on density    ==
    ! ==--------------------------------------------------------------==
    IF (dolinres) THEN
       ! Recalculate density, but re-use POTR, FNL, DFNL, EIVPS, and EIROP
       ! IF recalculate FNL, DFNL, EIVPS, and EIROP THEN do next call
       !      CALL dcopy(nnr1*clsd%nlsd,rhoe,1,potr,1)
       ! RHOE must be computed using RELAXRHO for cntl%tddft
       CALL rhoofr(c0(:,1:nn),rhoe,psi(:,1),nn)
       CALL dcopy(fpar%nnr1*clsd%nlsd,rhoe,1,rhoo,1)
       IF (relaxrhoe) THEN
          loptib=lr01%lopti
          lr01%lopti=loptib-10
          nolrb=nolr
          nolr=nact
          CALL td_force(c0,clrtd,cr,sc0,eigv,rhoe,psi,eigs,&
               no,nroot,tau0,fion_na,'CANON',1)
          CALL relax_rho(c0,clrtd,cr,sc0,rhoo,rhoe,psi(:,1),no,nroot)
          lr01%lopti=loptib
          nolr=nolrb
          CALL dcopy(fpar%nnr1*clsd%nlsd,rhoe,1,rhoo,1)
       ENDIF
       ! Non-local pseudopotential, calculate DDFNL - allocated above
       CALL rnlsm_2d(c0,nn)
       ! Save FNL
       CALL fnl_set("SAVE")
       CALL fnlalloc(nn,.TRUE.,.FALSE.)
       ! Local pseudopotential, use structure factors calculated before
       CALL eicalc(eivps,eirop)
       ! Calculate and store analytic dmu/dn
       CALL lr_xcpot(ddxc,rhoo,.FALSE.)
       ! Preconditioning for LR optimization - cf. SDLINRES
       cntl%prec=.TRUE.
       CALL ksdiag(vpp)
       trace=dasum(nn,eigt,1)/REAL(nn,kind=real_8)
       !$omp parallel do private(IG,VP)
       DO ig=1,ncpw%ngw
          vp=vpp(ig)+trace
          vpp(ig)=ABS(vp/(vp**2+lr02%lr_hthrs**2))
       ENDDO
       ! ==-------------------------------------------------------------==
       ! ==  Unit vector along constant density coupling vector         ==
       ! ==  Result in TAUC and then VECT                               ==
       ! ==-------------------------------------------------------------==
       GOTO 1234
       CALL cplconfig(vcloca,tauc)
       rlen=0.0_real_8
       iat=0
       DO isp=1,ions1%nsp
          DO ia=1,ions0%na(isp)
             iat=iat + 1
             IF (tfdens.OR..NOT.tspecv) THEN
                DO k=1,3
                   vect(k,iat,1)=tauc(k,ia,isp)
                   rlen=rlen+vect(k,iat,1)**2
                ENDDO
             ENDIF
          ENDDO            ! IA = 1,NA(ISP)
       ENDDO               ! ISP = 1,NSP
       natoms=iat
       n3    =3*natoms
       IF (tfdens.OR..NOT.tspecv) THEN
          rlen=1.0_real_8/SQRT(rlen)
          DO iat=1,natoms
             vect(1,iat,1)=rlen*vect(1,iat,1)
             vect(2,iat,1)=rlen*vect(2,iat,1)
             vect(3,iat,1)=rlen*vect(3,iat,1)
          ENDDO
       ENDIF
1234   CONTINUE
       ! ==--------------------------------------------------------------==
       ! ==  Brute-force loop over all DOF of ions                       ==
       ! ==--------------------------------------------------------------==
       CALL zeroing(c1)!,ngw*nstate)
       ! 
       IF (talldof) THEN
          ! WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
          nn=no
          ! WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
          ! 
          ! ==  Main finite difference loop over the dimensions             ==
          ! 
          symmi%nrot=1
          iat = 0
          DO isp = 1,ions1%nsp
             DO ia = 1,ions0%na(isp)
                iat = iat + 1
                DO k = 1,3
                   IF ((lprint).AND.paral%io_parent)&
                        WRITE(6,'(A,3I3)')&
                        ' TYPE, ATOM, DIMENSION: ', isp, ia, k
                   ! 
                   ! ==  Linear response of nuclear displacement                     ==
                   ! 
                   ! The following lines stem from SDLINRES
                   eind=0.0_real_8
                   CALL eind_ii(tau0,isp,iat,k,eind,iteropt%iesr)
                   CALL eind_loc(eind,isp,iat,k,rhoo,psi(:,1),eirop)
                   ! Calculate constant force from NL-PP
                   CALL fnl_set("SWITCH")
                   CALL nl_res(isp,iat,k,h1nl,nn)
                   CALL eind_nl(eind,isp,iat,k,nn)
                   CALL fnl_set("SWITCH")
                   ! Calculate the constant part of the local potential
                   CALL eicalc1(k,isp,iat,eivps1,eirop1)
                   ! Calculate first-order wavefunction
                   ropt_mod%sdiis=.TRUE.
                   ropt_mod%spcg=.TRUE.
                   CALL opt_lr(c0,c1,c2,sc0,eind,eigv,rhoe,&
                        h1nl,eivps1,eirop1,ddxc,vpp,cscr,cr,&
                        psi,nn,"PHONON",orbtyp)
                   ! End SDLINRES - calculate first order density
                   CALL rho1ofr(c0,c1,crge%f(:,1),rhoe,psi(:,1),nn)
                   ! Hartree and XC first-order terms
                   CALL v1ofrho1(e2,rhoe,ddxc,psi)
                   ! CALL DSCAL(NNR1,-1.0_real_8,RHOE(1,1),1)
                   ! IF(cntl%tlsd) CALL DSCAL(NNR1,-1.0_real_8,RHOE(1,2),1)
                   ! CALL AZZERO(C2,2*NGWK*N)
                   ! CALL VPSIMT(C0,C2,F,RHOE,PSI,NN,NLSD,.TRUE.)
                   CALL print_ex_rho(rhoe,psi(:,1),tau0)
                   ! Coupling component for each pair of surfaces
                   ! IC=0 is the GS surface
                   ici=0
                   DO ic=0,nroot
                      DO icc=ic+1,nroot
                         ici=ici+1
                         ! if (parent) write(6,*)'omega ', ICI,EIGO(ICI)
                         DO id=1,nn
                            is1=id
                            IF (cntl%tlsd) THEN
                               IF (is1.LE.tdsp1%nupel) THEN
                                  is2=is1
                               ELSE
                                  is2=nap+is1
                               ENDIF
                            ELSE
                               is2=is1
                            ENDIF
                            CALL zeroing(rhoov)!,nnr1)
                            IF (ic.EQ.0) THEN
                               CALL rhoabofr(1,c0(:,is1:is1),clrtd(:,is2:is2,icc),rhoov,psi(:,1))
                               deinv = 1._real_8/eigo(icc)
                            ELSE
                               CALL rhoabofr(1,clrtd(:,is1:is1,ic),clrtd(:,is2:is2,icc),&
                                    rhoov,psi(:,1))
                               ! CALL RHOABOFR(1,C0(1,IS1),C0(1,IS2),RHOOV,PSI)
                               deinv = 1._real_8/(eigo(icc)-eigo(ic))
                            ENDIF
                            ! --------------------------------------------------------------------
                            IF (debug.AND.ic.GT.0) THEN
                               IF (paral%io_parent)&
                                    WRITE(6,*) '------------'
                               IF (paral%io_parent)&
                                    WRITE(6,*) "NUPEL,NAP,NLSD"
                               IF (paral%io_parent)&
                                    WRITE(6,*)  tdsp1%nupel,nap,clsd%nlsd
                               IF (paral%io_parent)&
                                    WRITE(6,*) 'NN,IC,ID,IS1,IS2'
                               IF (paral%io_parent)&
                                    WRITE(6,*)  nn,ic,id,is1,is2
                               IF (paral%io_parent)&
                                    WRITE(6,*) '------------'
                               fcompa=dotp(ncpw%ngw,clrtd(:,is2,ic),clrtd(:,is2,ic))
                               CALL mp_sum(fcompa,parai%allgrp)
                               IF (paral%io_parent)&
                                    WRITE(6,*)'module CLRTD',is2
                               IF (paral%io_parent)&
                                    WRITE(6,*)'module CLRTD',fcompa
                               fcompa=dotp(ncpw%ngw,c0(:,is1),clrtd(:,is2,ic))
                               CALL mp_sum(fcompa,parai%allgrp)
                               IF (paral%io_parent)&
                                    WRITE(6,*) 'scalar product C0 CLRTD',is1,is2
                               IF (paral%io_parent)&
                                    WRITE(6,*) 'scalar product C0 CLRTD',fcompa
                               fcompa=0._real_8
                               !$omp parallel do private(IR) reduction(+:FCOMPA) schedule(static)
                               DO ir=1,fpar%nnr1
                                  fcompa=fcompa+rhoov(ir)
                               ENDDO
                               fcompa=fcompa*parm%omega/(spar%nr1s*spar%nr2s*spar%nr3s)
                               CALL mp_sum(fcompa,parai%allgrp)
                               IF (paral%io_parent)&
                                    WRITE(6,*) 'int RHOOV',fcompa
                               fcompa=0._real_8
                               !$omp parallel do private(IR) reduction(+:FCOMPA) schedule(static)
                               DO ir=1,fpar%nnr1
                                  fcompa=fcompa+rhoov(ir)**2
                               ENDDO
                               fcompa=fcompa*parm%omega/(spar%nr1s*spar%nr2s*spar%nr3s)
                               CALL mp_sum(fcompa,parai%allgrp)
                               IF (paral%io_parent)&
                                    WRITE(6,*) 'int RHOOV2',fcompa

                               fcompa=0._real_8
                               !$omp parallel do private(IR) reduction(+:FCOMPA) schedule(static)
                               DO ir=1,fpar%nnr1
                                  fcompa=fcompa+rhoe(ir,1)
                               ENDDO
                               fcompa=fcompa*parm%omega/(spar%nr1s*spar%nr2s*spar%nr3s)
                               CALL mp_sum(fcompa,parai%allgrp)
                               IF (paral%io_parent)&
                                    WRITE(6,*) 'int DPOTalpha',fcompa
                               fcompa=0._real_8
                               !$omp parallel do private(IR) reduction(+:FCOMPA) schedule(static)
                               DO ir=1,fpar%nnr1
                                  fcompa=fcompa+rhoe(ir,2)
                               ENDDO
                               fcompa=fcompa*parm%omega/(spar%nr1s*spar%nr2s*spar%nr3s)
                               CALL mp_sum(fcompa,parai%allgrp)
                               IF (paral%io_parent)&
                                    WRITE(6,*) 'int DPOTbeta',fcompa
                            ENDIF
                            ! --------------------------------------------------------------------
                            IF (lspin2%tlse) THEN
                               fcompa=0.5_real_8*(&
                                    DDOT(fpar%nnr1,RHOOV,1,RHOE(1,2),1) +&
                                    DDOT(fpar%nnr1,RHOOV,1,RHOE(1,3),1))
                               fcompa=(fcompa*parm%omega/REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)) * deinv
                               CALL mp_sum(fcompa,parai%allgrp)
                            ELSEIF (cntl%tlsd) THEN
                               IF (is1.LE.tdsp1%nupel) THEN
                                  ! FCOMPA=DOTP(NGW,CLRTD(1,IS2,ICI),C2(1,IS1))
                                  ! FCOMPA=DOTP(NGW,C0(1,IS2),C2(1,IS1))/ EIGO(ICI)
                                  fcompa=ddot(fpar%nnr1,rhoov,1,rhoe(1,1),1)
                                  fcompa=(fcompa*parm%omega/REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)) * deinv
                                  CALL mp_sum(fcompa,parai%allgrp)
                                  fcompb=0.0_real_8
                               ELSE
                                  ! FCOMPB=DOTP(NGW,CLRTD(1,IS2,ICI),C2(1,IS1))
                                  ! FCOMPB=DOTP(NGW,C0(1,IS2),C2(1,IS1))/ EIGO(ICI)
                                  fcompb=ddot(fpar%nnr1,rhoov,1,rhoe(1,2),1)
                                  fcompb=(fcompb*parm%omega/REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)) * deinv
                                  CALL mp_sum(fcompb,parai%allgrp)
                                  fcompa=0.0_real_8
                               ENDIF
                            ELSE
                               ! FCOMPA=DOTP(NGW,CLRTD(1,IS2,ICI),C2(1,IS1))
                               ! FCOMPA=FCOMPA / EIGO(ICI) ! * DEINV
                               fcompa=ddot(fpar%nnr1,rhoov,1,rhoe(1,1),1)
                               fcompa=(fcompa*parm%omega/REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)) * deinv
                               CALL mp_sum(fcompa,parai%allgrp)
                            ENDIF


                            IF (lprint.AND.paral%parent) THEN
                               IF (paral%io_parent)&
                                    WRITE(6,'(A,2F12.5)') ' FIXED DENSITY AND RESPONSE: ',&
                                    VCLOCA(K,IA,ISP,ICI), FCOMPA
                               IF (paral%io_parent)&
                                    WRITE(6,'(A,F12.5)') ' TOTAL COUPLING:             ',&
                                    VCLOCA(K,IA,ISP,ICI) + FCOMPA
                            ENDIF
                            vcloca(k,ia,isp,ici)=vcloca(k,ia,isp,ici)+fcompa
                            IF (cntl%tlsd) THEN
                               vclocb(k,ia,isp,ici)=vclocb(k,ia,isp,ici)+fcompb
                            ENDIF
                            ! 
                            ! FCOMP  = DOTP (NGW, C0(1,IS1), C1(1,IS2))
                            ! FCOMP2 = DOTP (NGW, C0(1,IS2), C1(1,IS1))
                            ! call mp_sum(FCOMP,1,allgrp)
                            ! call mp_sum(FCOMP2,1,allgrp)
                            ! IF (LPRINT) WRITE(6,'(A,2F12.5)')
                            !               'COUPLINGS FROM OVERLAP: ', FCOMP, FCOMP
                         ENDDO
                      ENDDO     ! ICC = IC+1,NROOT
                   ENDDO       ! IC = 1,NROOT
                   IF (.FALSE.) THEN
                      DO i=1,nstate
                         CALL zeroing(rhoov)!,nnr1)
                         CALL rhoabofr (1,c0(:,i:i),c0(:,i:i),rhoov,psi(:,1))
                         fcomp=ddot(fpar%nnr1,rhoov,1,rhoe(1,1),1)
                         CALL mp_sum(fcomp,parai%allgrp)
                         fcomp=fcomp*parm%omega/(spar%nr1s*spar%nr2s*spar%nr3s)
                         IF (paral%io_parent)&
                              WRITE(6,*)'State ', i, ': ', fcomp
                      ENDDO
                   ENDIF
                ENDDO          ! K = 1,3
             ENDDO              ! IA = 1,NA(ISP)
          ENDDO                  ! ISP = 1,NSP
          ! ==--------------------------------------------------------------==
          ! ==  Iterative scheme for contributions of response of density   ==
          ! ==--------------------------------------------------------------==
       ELSE                     ! IF (TALLDOF)
          ! 
          ! ==  Response of density along the direction at constant density ==
          ! 
          IF ((lprint).AND.paral%io_parent)&
               WRITE(6,'(A,A)') ' LINEAR RESPONSE ALONG ',&
               'CONSTANT DENSITY COUPLING VECTOR'
          iv=1
          ! Clear history of searched components
          !$omp parallel do private(K)
          DO k=1,3*ions1%nat
             idone(k) = 0
          ENDDO
          ! 
          ! ==  Iterative search of the direction                           ==
          ! 
10        CONTINUE           ! DO ... UNTIL (IV.GE.NVECT .OR. LCONV)
          IF (lprint) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(A,I3,A)') ' VECTOR ', iv, ':'
             !vw commented out following line as natoms is not defined !
             !vwIF (paral%io_parent)&
             !vw     WRITE(6,'(2(2X,3F8.4))')&
             !vw     ((VECT(K,IAT,IV), K=1,3), IAT=1,NATOMS)
          ENDIF
          ! Calculate constant part of local and non-local interactions
          CALL get_eind (vect(1,1,iv), eind, eivps1, eirop1, h1nl,&
               TAU0, iteropt%iesr, RHOO, PSI, NSTATE, EIROP)
          ! Calculate first-order wavefunction
          ropt_mod%sdiis = .TRUE.
          ropt_mod%spcg  = .TRUE.
          CALL opt_lr (c0, c1, c2, sc0, eind, eigv, rhoe,&
               H1NL, EIVPS1, EIROP1, DDXC, VPP, CSCR, CR,&
               PSI, NSTATE, "PHONON", ORBTYP)
          ! Calculate first order density
          CALL rho1ofr (c0, c1, crge%f(:,1), rhoe, psi(:,1), nstate)
          ! Hartree and XC first-order terms
          CALL v1ofrho1(e2, rhoe, ddxc, psi)
          ! Coupling component for each pair of surfaces
          fcsum = 0.0_real_8
          DO ic = 1,nroot
             ! CC            IS1 = ISURF(1,IC)
             ! CC            IS2 = ISURF(2,IC)
             IF (ic.GT.1) THEN
                CALL zeroing(rhoov)!, nnr1)
                CALL rhoabofr (1, c0(:,is1:is1), c0(:,is2:is2), rhoov, psi(:,1))
                deinv = deinvv(ic)
             ENDIF
             IF (lspin2%tlse) THEN
                fcomp = 0.5_real_8 * (&
                     DDOT (fpar%nnr1, RHOOV, 1, RHOE(1,2), 1) +&
                     DDOT (fpar%nnr1, RHOOV, 1, RHOE(1,3), 1))
             ELSE IF (cntl%tlsd.AND.is2.GT.spin_mod%nsup) THEN
                fcomp = ddot (fpar%nnr1, rhoov, 1, rhoe(1,2), 1)
             ELSE
                fcomp = ddot (fpar%nnr1, rhoov, 1, rhoe(1,1), 1)
             ENDIF
             fcomp = fcomp * parm%omega/(spar%nr1s*spar%nr2s*spar%nr3s) * deinv
             CALL mp_sum(fcomp,parai%allgrp)
             fcsrf(ic) = fcomp
             fcsum = fcsum + csurf(ic)*fcomp
             IF ((lprint).AND.paral%io_parent)&
                  WRITE(6,'(A,I3,A,F12.5)')&
                  ' RESPONSE ALONG VECTOR ', IV, ': ', FCOMP
          ENDDO            ! IC = 1,ROOT
          ! 
          ! ==  Response along search vector large, increase accuracy       ==
          ! 
          IF (iv.EQ.1 .OR. ABS(fcsum).GT.fc_low) THEN
             TYPE = "NOINIT"
             tol_lr_bak = lr02%tol_lr
             IF (iv.EQ.1 .OR. ABS(fcsum).GT.fc_high) THEN
                lr02%tol_lr = lr02%tol_lr / f_high
             ELSE IF (ABS(fcsum).GT.fc_med) THEN
                lr02%tol_lr = lr02%tol_lr / f_med
             ELSE
                lr02%tol_lr = lr02%tol_lr / f_low
             ENDIF
             IF ((lprint).AND.paral%io_parent)&
                  WRITE(6,'(A,F12.5)')&
                  ' INCREASED ACCURACY: ', lr02%tol_lr
             CALL opt_lr (c0, c1, c2, sc0, eind, eigv, rhoe,&
                  H1NL, EIVPS1, EIROP1, DDXC, VPP, CSCR, CR,&
                  PSI, NSTATE, TYPE, ORBTYP)
             CALL rho1ofr (c0, c1, crge%f(:,1), rhoe, psi(:,1), nstate)
             ! Hartree and XC first-order terms
             CALL v1ofrho1(e2, rhoe, ddxc, psi)
             ! Coupling component for each pair of surfaces
             DO ic = 1,nroot
                ! CC              IS1 = ISURF(1,IC)
                ! CC              IS2 = ISURF(2,IC)
                IF (ic.GT.1) THEN
                   CALL zeroing(rhoov)!, nnr1)
                   CALL rhoabofr (1, c0(:,is1:is1), c0(:,is2:is2), rhoov, psi(:,1))
                   deinv = deinvv(ic)
                ENDIF
                IF (lspin2%tlse) THEN
                   fcomp = 0.5_real_8 * (&
                        DDOT (fpar%nnr1, RHOOV, 1, RHOE(1,2), 1) +&
                        DDOT (fpar%nnr1, RHOOV, 1, RHOE(1,3), 1))
                ELSE IF (cntl%tlsd.AND.is2.GT.spin_mod%nsup) THEN
                   fcomp = ddot (fpar%nnr1, rhoov, 1, rhoe(1,2), 1)
                ELSE
                   fcomp = ddot (fpar%nnr1, rhoov, 1, rhoe(1,1), 1)
                ENDIF
                fcomp = fcomp * parm%omega/(spar%nr1s*spar%nr2s*spar%nr3s) * deinv
                CALL mp_sum(fcomp,parai%allgrp)
                fcsrf(ic) = fcomp
                IF ((lprint).AND.paral%io_parent)&
                     WRITE(6,'(A,F12.5)')&
                     ' MORE ACCURATE RESPONSE: ', FCOMP
             ENDDO        ! IC = 1,NROOT
             lr02%tol_lr = tol_lr_bak
          ENDIF            ! IF (IV.EQ.1 .OR. FCSUM.GT.FC_LOW)
          ! 
          ! ==  Store response                                              ==
          ! 
          DO ic = 1,nroot
             iat = 0
             DO isp = 1,ions1%nsp
                DO ia = 1,ions0%na(isp)
                   iat = iat + 1
                   DO k = 1,3
                      fc = fcsrf(ic) * vect(k,iat,iv)
                      vclocb(k,ia,isp,ic) = vclocb(k,ia,isp,ic) + fc
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
          ! 
          ! ==  Not converged, find new search direction                    ==
          ! 
          lconv = (ABS(fcomp).LT.tolcpl)
          IF (iv.LT.nvect .AND. .NOT.lconv) THEN
             IF (lprint) THEN
                IF (paral%io_parent)&
                     WRITE(6,'(A,I3)') ' COUPLINGS, ITERATION: ', iv
                CALL prcplngs (vclocb, tauc, .FALSE.)
             ENDIF
             iv = iv + 1
             ! Skip generation of new vector if specified
             IF (.NOT.tspecv) THEN
                ! Look for maximum component not yet searched
                rmax = 0.0_real_8
                imax = 0
                kmax = 0
                CALL zeroing(vect(:,:,iv))!, n3)
                DO iat = 1,natoms
                   DO k = 1,3
                      IF (idone(k+3*(iat-1)).EQ.0) THEN
                         IF (ABS(vect(k,iat,iv-1)).GT.rmax) THEN
                            rmax = ABS(vect(k,iat,iv-1))
                            imax = iat
                            kmax = k
                         ENDIF
                      ENDIF
                   ENDDO
                ENDDO
                IF ((lprint).AND.paral%io_parent)&
                     WRITE(6,*) 'COMPONENT: ', imax, kmax
                ! Select and mark new component as searched
                idone(kmax+3*(imax-1)) = iv
                vect(kmax,imax,iv) = 1.0_real_8
             ENDIF        ! (.NOT.TSPECV)
             ! Orthogonalize new vector with respect to all previous ones
             rlen = SQRT(ddot(n3, vect(1,1,iv), 1, vect(1,1,iv), 1))
             CALL dscal (n3, 1.0_real_8/rlen, vect(1,1,iv), 1)
             DO ivv = 1,iv-1
                rlen = ddot (n3, vect(1,1,ivv), 1, vect(1,1,iv), 1)
                CALL daxpy (n3, -rlen, vect(1,1,ivv), 1, vect(1,1,iv), 1)
                rlen = SQRT(ddot(n3, vect(1,1,iv), 1, vect(1,1,iv), 1))
                CALL dscal (n3, 1.0_real_8/rlen, vect(1,1,iv), 1)
             ENDDO
             ! Report singular transformation
             DO ivv = 1,iv-1
                rlen = ddot (n3, vect(1,1,ivv), 1, vect(1,1,iv), 1)
                IF ((lprint2).AND.paral%io_parent)&
                     WRITE(6,'(A,I3,A,F10.6)')&
                     ' OVERLAP WITH VECTOR ', IVV, ': ', RLEN
                IF (paral%parent .AND. ABS(rlen).GT.1.0e-6_real_8) THEN
                   IF (paral%io_parent)&
                        WRITE(6,'(A,I3,A,I3,A,E10.4)')&
                        ' WARNING: OVERLAP BETWEEN VECTOR ', IVV, ' AND ',&
                        IV, ': ', RLEN
                ENDIF
             ENDDO
             GOTO 10
             ! 
             ! ==  Done with iterative scheme for one or the other reason      ==
             ! 
          ELSE IF (lconv) THEN
             IF ((lprint).AND.paral%io_parent)&
                  WRITE(6,'(A,A,I3,A)') ' ITERATIVE SCHEME FOR ',&
                  'COUPLINGS CONVERGED IN ', IV, ' CYCLES'
          ELSE IF (iv.GE.nvect .AND. tolcpl.GT.0.0_real_8) THEN
             IF ((lprint).AND.paral%io_parent)&
                  WRITE(6,'(A)')&
                  ' WARNING: ITERATIVE SCHEME FOR COUPLINGS NOT CONVERGED'
          ELSE IF (iv.GE.nvect) THEN
             IF ((lprint.AND.tolcpl.GT.0.0_real_8).AND.paral%io_parent)&
                  WRITE(6,'(A,I3)')&
                  ' ITERATIVE SCHEME FOR COUPLINGS DONE, CYCLE ', IV
          ENDIF
          CONTINUE           ! enddo ... UNTIL (IV.GE.NVECT .OR. LCONV)
       ENDIF                ! IF (TALLDOF) ... ELSE ...
    ENDIF  ! IF (DOLINRES) THEN ... ENDIF

    ! ==--------------------------------------------------------------==
    DEALLOCATE(ddfnl,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    !   DEALLOCATE(potr,STAT=ierr)
    !   IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
    !        __LINE__,__FILE__)
    DEALLOCATE(eivps,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(eirop,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(eivps1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(eirop1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(h1nl,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(rhoov,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(ddxc,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(deinvv,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(vect,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    IF (.NOT.talldof) THEN
       DEALLOCATE(idone,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(fcsrf,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    ENDIF

    DEALLOCATE(c1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(c2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(psi,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(vpp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(cscr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(tauc,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(eigt,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(eigs,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(fnlc,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(dfnlc,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    CALL tihalt('   LRCPLNG',isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE vcouplings_lr
  ! ==================================================================
  SUBROUTINE drhopp(vcloca,vclocb,fnlc,dfnlc,nn,nap,ic)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==  THE NON-LOCAL POTENTIAL CONTRIBUTION TO THE FORCE ON THE    ==
    ! ==  IONIC DEGREES OF FREEDOM                                    ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: vcloca(:,:,:,:), &
                                                vclocb(:,:,:,:), &
                                                fnlc(:,:,:,:,:), &
                                                dfnlc(:,:,:,:,:,:)
    INTEGER                                  :: nn, nap, ic

    CHARACTER(*), PARAMETER                  :: procedureN = 'drhopp'

    INTEGER                                  :: ia, iatl, id, ierr, ii, is, &
                                                is1, is2, isa, isa0, isub, &
                                                iv, jv, k, ki, kj, l, l2, li, &
                                                lj
    REAL(real_8)                             :: tt
    REAL(real_8), ALLOCATABLE                :: dfab(:,:,:)

! NLM
! NLM
! Variables
! ==--------------------------------------------------------------==
! If no non-local components -> return.

    IF (nlm.EQ.0) RETURN
    ! ==--------------------------------------------------------------==
    CALL tiset('    DRHOPP',isub)
    ALLOCATE(dfab(ions1%nat,maxsys%nhxs,2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! DO IC=1,NROOT
    DO id=1,nn
       is1=id
       IF (cntl%tlsd) THEN
          IF (is1.LE.tdsp1%nupel) THEN
             is2=is1
          ELSE
             is2=nap+is1
          ENDIF
       ELSE
          is2=is1
       ENDIF
       DO k=1,3
          CALL zeroing(dfab)!,2*ions1%nat*maxsys%nhxs)
          IF (is1.GE.parap%nst12(parai%mepos,1).AND.&
               is1.LE.parap%nst12(parai%mepos,2) ) THEN
             ii=is1-parap%nst12(parai%mepos,1)+1
             CALL dcopy(maxsys%nhxs*ions1%nat,dfnl(1,1,1,k,ii,1),1,dfab(1,1,1),1)
          ENDIF
          IF (is2.GE.parap%nst12(parai%mepos,1).AND.&
               is2.LE.parap%nst12(parai%mepos,2) ) THEN
             ii=is2-parap%nst12(parai%mepos,1)+1
             CALL dcopy(maxsys%nhxs*ions1%nat,dfnlc(1,1,1,k,ii,1),1,dfab(1,1,2),1)
          ENDIF
          CALL mp_sum(dfab,2*ions1%nat*maxsys%nhxs,parai%allgrp)
          isa0=0
          DO is=1,ions1%nsp
             IF (pslo_com%tvan(is)) THEN
                ! Vanderbild pp
                CALL stopgm("DRHOPP","VDB NOT IMPLEMENTED",&
                     __LINE__,__FILE__)
             ELSEIF (sgpp1%tsgp(is)) THEN
                ! Stefan Goedecker pp
                CALL stopgm("DRHOPP","SG NOT YET WORKING",&
                     __LINE__,__FILE__)
                DO iv=1,nlps_com%ngh(is)
                   l=nghtol(iv,is)+1
                   ki=sgpp2%lfval(iv,is)
                   li=sgpp2%lpval(iv,is)
                   DO jv=1,nlps_com%ngh(is)
                      l2=nghtol(jv,is)+1
                      lj=sgpp2%lpval(jv,is)
                      IF (l2.EQ.l.AND.li.EQ.lj) THEN
                         kj=sgpp2%lfval(jv,is)
                         tt=sgpp2%hlsg(ki,kj,l,is)
                         DO ia=1,ions0%na(is)
                            isa=isa0+ia
                            IF (cntl%tfdist) THEN
                               iatl=isa-ipept(1,parai%mepos)+1
                               IF (cntl%tlsd) THEN
                                  IF (is1.LE.tdsp1%nupel) THEN
                                     vcloca(k,ia,is,ic)=vcloca(k,ia,is,ic)&
                                          -0.5_real_8*tt*(&
                                          dfab(isa,iv,1)*fnlc(1,iatl,jv,is2,1)&
                                          +dfab(isa,iv,2)*fnl(1,iatl,jv,is1,1)&
                                          +dfab(isa,jv,1)*fnlc(1,iatl,iv,is2,1)&
                                          +dfab(isa,jv,2)*fnl(1,iatl,iv,is1,1))
                                  ELSE
                                     vclocb(k,ia,is,ic)=vclocb(k,ia,is,ic)&
                                          -0.5_real_8*tt*(&
                                          dfab(isa,iv,1)*fnlc(1,iatl,jv,is2,1)&
                                          +dfab(isa,iv,2)*fnl(1,iatl,jv,is1,1)&
                                          +dfab(isa,jv,1)*fnlc(1,iatl,iv,is2,1)&
                                          +dfab(isa,jv,2)*fnl(1,iatl,iv,is1,1))
                                  ENDIF
                               ELSE
                                  vcloca(k,ia,is,ic)=vcloca(k,ia,is,ic)&
                                       -0.5_real_8*tt*(&
                                       dfab(isa,iv,1)*fnlc(1,iatl,jv,is2,1)&
                                       +dfab(isa,iv,2)*fnl(1,iatl,jv,is1,1)&
                                       +dfab(isa,jv,1)*fnlc(1,iatl,iv,is2,1)&
                                       +dfab(isa,jv,2)*fnl(1,iatl,iv,is1,1))
                               ENDIF
                            ELSE
                               IF (cntl%tlsd) THEN
                                  IF (is1.LE.tdsp1%nupel) THEN
                                     vcloca(k,ia,is,ic)=vcloca(k,ia,is,ic)-tt*(&
                                          dfab(isa,iv,1)*fnlc(1,isa,jv,is2,1)&
                                          +dfab(isa,iv,2)*fnl(1,isa,jv,is1,1)&
                                          +dfab(isa,jv,1)*fnlc(1,isa,iv,is2,1)&
                                          +dfab(isa,jv,2)*fnl(1,isa,iv,is1,1))
                                  ELSE
                                     vclocb(k,ia,is,ic)=vclocb(k,ia,is,ic)-tt*(&
                                          dfab(isa,iv,1)*fnlc(1,isa,jv,is2,1)&
                                          +dfab(isa,iv,2)*fnl(1,isa,jv,is1,1)&
                                          +dfab(isa,jv,1)*fnlc(1,isa,iv,is2,1)&
                                          +dfab(isa,jv,2)*fnl(1,isa,iv,is1,1))
                                  ENDIF
                               ELSE
                                  vcloca(k,ia,is,ic)=vcloca(k,ia,is,ic)-tt*(&
                                       dfab(isa,iv,1)*fnlc(1,isa,jv,is2,1)&
                                       +dfab(isa,iv,2)*fnl(1,isa,jv,is1,1)&
                                       +dfab(isa,jv,1)*fnlc(1,isa,iv,is2,1)&
                                       +dfab(isa,jv,2)*fnl(1,isa,iv,is1,1))
                               ENDIF
                            ENDIF
                         ENDDO
                      ENDIF
                   ENDDO
                ENDDO
             ELSE
                ! Other pp (numeric)
                DO iv=1,nlps_com%ngh(is)
                   DO ia=1,ions0%na(is)
                      isa=isa0+ia
                      IF (isa.GE.ipept(1,parai%mepos).AND.&
                           isa.LE.ipept(2,parai%mepos)) THEN
                         IF (cntl%tfdist) THEN
                            iatl=isa-ipept(1,parai%mepos)+1
                            IF (cntl%tlsd) THEN
                               IF (is1.LE.tdsp1%nupel) THEN
                                  vcloca(k,ia,is,ic)=vcloca(k,ia,is,ic)&
                                       -wsg(is,iv)*&
                                       (dfab(isa,iv,1)*fnlc(1,iatl,iv,is2,1)&
                                       +fnl(1,iatl,iv,is1,1)*dfab(isa,iv,2))
                               ELSE
                                  vclocb(k,ia,is,ic)=vclocb(k,ia,is,ic)&
                                       -wsg(is,iv)*&
                                       (dfab(isa,iv,1)*fnlc(1,iatl,iv,is2,1)&
                                       +fnl(1,iatl,iv,is1,1)*dfab(isa,iv,2))
                               ENDIF
                            ELSE
                               vcloca(k,ia,is,ic)=vcloca(k,ia,is,ic)&
                                    -wsg(is,iv)*&
                                    (dfab(isa,iv,1)*fnlc(1,iatl,iv,is2,1)&
                                    +fnl(1,iatl,iv,is1,1)*dfab(isa,iv,2))
                            ENDIF
                         ELSE
                            IF (cntl%tlsd) THEN
                               IF (is1.LE.tdsp1%nupel) THEN
                                  vcloca(k,ia,is,ic)=vcloca(k,ia,is,ic)&
                                       -wsg(is,iv)*&
                                       (dfab(isa,iv,1)*fnlc(1,isa,iv,is2,1)&
                                       +fnl(1,isa,iv,is1,1)*dfab(isa,iv,2))
                               ELSE
                                  vclocb(k,ia,is,ic)=vclocb(k,ia,is,ic)&
                                       -wsg(is,iv)*&
                                       (dfab(isa,iv,1)*fnlc(1,isa,iv,is2,1)&
                                       +fnl(1,isa,iv,is1,1)*dfab(isa,iv,2))
                               ENDIF
                            ELSE
                               vcloca(k,ia,is,ic)=vcloca(k,ia,is,ic)&
                                    -wsg(is,iv)*&
                                    (dfab(isa,iv,1)*fnlc(1,isa,iv,is2,1)&
                                    +fnl(1,isa,iv,is1,1)*dfab(isa,iv,2))
                            ENDIF
                         ENDIF
                      ENDIF
                   ENDDO
                ENDDO
             ENDIF
             isa0 = isa0 + ions0%na(is)
          ENDDO
       ENDDO
    ENDDO
    ! ENDDO
    DEALLOCATE(dfab,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    CALL tihalt('    DRHOPP',isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE drhopp
  ! ==================================================================
  SUBROUTINE tdrhofc(fion,rho,psi)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: fion(:,:,:), rho(:)
    COMPLEX(real_8)                          :: psi(:)

    COMPLEX(real_8)                          :: ei123, gx, gy, gz, rhets, &
                                                txx, tyy, tzz, vcgs
    INTEGER                                  :: ia, ig, ig1, is, isa
    REAL(real_8)                             :: omtp

! Variables
! &           EI3(NATX,(2*NR3S-1))
! ==--------------------------------------------------------------==
! DO IR=1,NNR1
! PSI(IR)=CMPLX(RHO(IR),0._real_8)
! ENDDO
! CALL FWFFT(PSI)
! 

    omtp=2._real_8*parm%omega*parm%tpiba
    ig1=1
    IF (geq0) ig1=2
    DO ig=ig1,ncpw%nhg
       rhets=CONJG(psi(nzh(ig)))
       gx=CMPLX(0._real_8,gk(1,ig),kind=real_8)
       gy=CMPLX(0._real_8,gk(2,ig),kind=real_8)
       gz=CMPLX(0._real_8,gk(3,ig),kind=real_8)
       vcgs=scg(ig)*rhets
       isa=0
       DO is=1,ions1%nsp
          txx=(rhops(is,ig)*vcgs+vps(is,ig)*rhets)*gx
          tyy=(rhops(is,ig)*vcgs+vps(is,ig)*rhets)*gy
          tzz=(rhops(is,ig)*vcgs+vps(is,ig)*rhets)*gz
          DO ia=1,ions0%na(is)
             isa=isa+1
             IF (cntl%bigmem) THEN
                ei123=eigrb(ig,isa)
             ELSE
                ei123=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
                     ei3(isa,inyh(3,ig))
             ENDIF
             fion(1,ia,is)=fion(1,ia,is)+REAL(ei123*txx)*omtp
             fion(2,ia,is)=fion(2,ia,is)+REAL(ei123*tyy)*omtp
             fion(3,ia,is)=fion(3,ia,is)+REAL(ei123*tzz)*omtp
          ENDDO
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
  END SUBROUTINE tdrhofc
  ! ==================================================================
  SUBROUTINE  nacstrans(vcloca,vclocb,nroot)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: vcloca(:,:,:,:), &
                                                vclocb(:,:,:,:)
    INTEGER                                  :: nroot

    CHARACTER(*), PARAMETER                  :: procedureN = 'nacstrans'

    INTEGER                                  :: ia, ierr, iroot, is, k, niter
    REAL(real_8)                             :: ca(3,100), cb(3,100)
    REAL(real_8), ALLOCATABLE                :: vhelpa(:,:,:), vhelpb(:,:,:)

! Variables
! 
! 

    ALLOCATE(vhelpa(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    IF (cntl%tlsd) THEN
       ALLOCATE(vhelpb(3,maxsys%nax,maxsys%nsx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ELSE
       ALLOCATE(vhelpb(3,maxsys%nax,maxsys%nsx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! 
    niter= (nroot*(nroot-1))/2
    DO iroot=1,niter
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             DO k=1,3
                vhelpa(k,ia,is)=vcloca(k,ia,is,iroot)
             ENDDO
          ENDDO
       ENDDO
       IF (cntl%tlsd) THEN
          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                DO k=1,3
                   vhelpb(k,ia,is)=vclocb(k,ia,is,iroot)
                ENDDO
             ENDDO
          ENDDO
       ENDIF
       ! 
       DO k=1,3
          ca(k,iroot)=0.0
          cb(k,iroot)=0.0
       ENDDO
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             DO k=1,3
                ca(k,iroot)=ca(k,iroot)+vhelpa(k,ia,is)
             ENDDO
          ENDDO
       ENDDO
       IF (cntl%tlsd) THEN
          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                DO k=1,3
                   cb(k,iroot)=cb(k,iroot)+vhelpb(k,ia,is)
                ENDDO
             ENDDO
          ENDDO
       ENDIF
       ! 
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             DO k=1,3
                vcloca(k,ia,is,iroot)=vhelpa(k,ia,is)-ca(k,iroot)/ions1%nat
             ENDDO
             IF (cntl%tlsd) THEN
                DO k=1,3
                   vclocb(k,ia,is,iroot)=vhelpb(k,ia,is)-cb(k,iroot)/ions1%nat
                ENDDO
             ENDIF
          ENDDO
       ENDDO
       ! 
       ! CALL ROTVEC(TAU0,VCLOCA(1,1,1,IROOT))
       ! IF (cntl%tlsd) CALL ROTVEC(TAU0,VCLOCB(1,1,1,IROOT))
    ENDDO
    ! 
  END SUBROUTINE nacstrans
  ! ==================================================================
  SUBROUTINE rotvec(tau0,velp)
    ! ==--------------------------------------------------------------==
    ! REMOVE FIX ROTATIONS
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:), velp(:,:,:)

    INTEGER                                  :: ia, iat, ierr, is
    REAL(real_8) :: cm(3), ekinp1, evl(3), evm(3,3), lm(3), ltens(3,3), &
      pscr(maxsys%nax*maxsys%nsx), pscrt, r1, r2, small, tempp1, &
      tscr(3,maxsys%nax*ions1%nsp), vscr(3,maxsys%nax*maxsys%nsx), vtmp(3), &
      wk(9), wx(3)

    PARAMETER (small=1.0e-10)  ! small angular momentum where we stop considering it
    ! ==--------------------------------------------------------------==
    ! COMPUTE THE IONIC TEMPERATURE TEMPP1
    CALL ekinpp(ekinp1,velp)
    tempp1=ekinp1*factem*2.0_real_8/glib

    ! GET CENTER OF MASS AND COPY TO HELPER ARRAY
    CALL zeroing(cm)!,3)
    iat=0
    pscrt=0.0_real_8
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          iat=iat+1
          r1=rmass%pma0(is)
          pscrt=pscrt+r1
          pscr(iat)=r1
          r2=tau0(1,ia,is)
          tscr(1,iat)=r2
          cm(1)=cm(1)+r2*r1
          r2=tau0(2,ia,is)
          tscr(2,iat)=r2
          cm(2)=cm(2)+r2*r1
          r2=tau0(3,ia,is)
          tscr(3,iat)=r2
          cm(3)=cm(3)+r2*r1
          vscr(1,iat)=velp(1,ia,is)
          vscr(2,iat)=velp(2,ia,is)
          vscr(3,iat)=velp(3,ia,is)
       ENDDO
    ENDDO
    cm(1)=cm(1)/pscrt
    cm(2)=cm(2)/pscrt
    cm(3)=cm(3)/pscrt

    ! translate C.O.M to origin.
    !$omp parallel do private(IAT)
    DO iat=1,ions1%nat
       tscr(1,iat)=tscr(1,iat)-cm(1)
       tscr(2,iat)=tscr(2,iat)-cm(2)
       tscr(3,iat)=tscr(3,iat)-cm(3)
    ENDDO

    ! get total angular momentum
    CALL zeroing(lm)!,3)
    DO iat=1,ions1%nat
       CALL vecprod(tscr(1,iat),vscr(1,iat),vtmp)
       r1=pscr(iat)
       lm(1)=lm(1)+vtmp(1)*r1
       lm(2)=lm(2)+vtmp(2)*r1
       lm(3)=lm(3)+vtmp(3)*r1
    ENDDO

    ! get moment of inertia
    CALL zeroing(ltens)!,9)
    DO iat=1,ions1%nat
       ltens(1,1)=ltens(1,1)+pscr(iat)*(tscr(2,iat)**2+tscr(3,iat)**2)
       ltens(1,2)=ltens(1,2)-pscr(iat)*(tscr(1,iat)*tscr(2,iat))
       ltens(1,3)=ltens(1,3)-pscr(iat)*(tscr(1,iat)*tscr(3,iat))

       ltens(2,1)=ltens(2,1)-pscr(iat)*(tscr(1,iat)*tscr(2,iat))
       ltens(2,2)=ltens(2,2)+pscr(iat)*(tscr(1,iat)**2+tscr(3,iat)**2)
       ltens(2,3)=ltens(2,3)-pscr(iat)*(tscr(2,iat)*tscr(3,iat))

       ltens(3,1)=ltens(3,1)-pscr(iat)*(tscr(1,iat)*tscr(3,iat))
       ltens(3,2)=ltens(3,2)-pscr(iat)*(tscr(2,iat)*tscr(3,iat))
       ltens(3,3)=ltens(3,3)+pscr(iat)*(tscr(1,iat)**2+tscr(2,iat)**2)
    ENDDO

    ierr=0
    CALL dcopy(9,ltens,1,evm,1)
    CALL dsyev('V','U',3,evm,3,evl,wk,9,ierr)

    ! get angular velocity in body frame. ignore if moment of inertia is small.
    CALL dgemv('T',3,3,1.0_real_8,evm,3,lm,1,0.0_real_8,wx,1)
    IF (ABS(evl(1)).GT.small) THEN
       wx(1)=wx(1)/evl(1)
    ELSE
       wx(1)=0.0_real_8
    ENDIF
    IF (ABS(evl(2)).GT.small) THEN
       wx(2)=wx(2)/evl(2)
    ELSE
       wx(2)=0.0_real_8
    ENDIF
    IF (ABS(evl(3)).GT.small) THEN
       wx(3)=wx(3)/evl(3)
    ELSE
       wx(3)=0.0_real_8
    ENDIF
    ! transform back to lab frame and subtract
    CALL dgemv('N',3,3,1.0_real_8,evm,3,wx,1,0.0_real_8,lm,1)
    iat=0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          iat=iat+1
          CALL vecprod(lm,tscr(1,iat),vtmp)
          velp(1,ia,is)=velp(1,ia,is)-vtmp(1)
          velp(2,ia,is)=velp(2,ia,is)-vtmp(2)
          velp(3,ia,is)=velp(3,ia,is)-vtmp(3)
       ENDDO
    ENDDO

    ! ==--------------------------------------------------------------==
  END SUBROUTINE rotvec
  ! ==================================================================
  SUBROUTINE writenacs(vibe,sder,ndim,tau0,nroot,switch)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: vibe(:), sder(:,:)
    INTEGER                                  :: ndim
    REAL(real_8)                             :: tau0(:,:,:)
    INTEGER                                  :: nroot, switch

    CHARACTER(len=10)                        :: a
    CHARACTER(len=100)                       :: file1
    CHARACTER(len=17)                        :: fr, rm(5)
    CHARACTER(len=2)                         :: typ
    CHARACTER(len=7)                         :: ck(3)
    INTEGER                                  :: at, he, i, ia, ii, is, iunit, &
                                                j, k, l, m, n, t(ions1%nat)
    LOGICAL                                  :: ferror, lw
    REAL(real_8)                             :: bx, by, bz, gz, &
                                                xc(ions1%nat), yc(ions1%nat), &
                                                zc(ions1%nat)

    file1='NACV.log'
    IF (cntl%tlsd) THEN
       IF (switch.EQ.0) THEN
          file1='NACVa.log'
       ELSE
          file1='NACVb.log'
       ENDIF
    ENDIF
    ! write gaussianformated output into VIB1.log and VIB2.log files,
    ! which are readable in molden/molekel to visualise the vibrations.
    fr=' Frequencies --  '
    rm(1)=' Red. masses --  '
    rm(2)=' Frc consts  --  '
    rm(3)=' IR Inten    --  '
    rm(4)=' Raman Activ --  '
    rm(5)=' Depolar     --  '
    typ='?A'
    a=' Atom AN  '
    ck(1)='      X'
    ck(2)='      Y'
    ck(3)='      Z'
    at=0
    gz=0._real_8
    lw=.FALSE.
    ferror=.FALSE.
    ! ---- cell dimensions---------------------------
    bx=cell_com%celldm(1)
    by=cell_com%celldm(2)
    bz=cell_com%celldm(3)
    IF (paral%io_parent)&
         CALL fileopen(15,file1,fo_def,ferror)
    IF (ferror) GOTO 999
    ! ---- read coordinates and atomtyps ------------
    k=1
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          xc(k)=tau0(1,ia,is)
          yc(k)=tau0(2,ia,is)
          zc(k)=tau0(3,ia,is)
          t(k)=ions0%iatyp(is)
          IF (xc(k).LT.0) lw=.TRUE.
          IF (yc(k).LT.0) lw=.TRUE.
          IF (zc(k).LT.0) lw=.TRUE.
          k=k+1
       ENDDO
    ENDDO
    ! ---- write some lines needed for molden
    iunit=15
    IF (paral%io_parent)&
         WRITE (iunit,*)'Entering Gaussian System'
    IF (paral%io_parent)&
         WRITE (iunit,*)'this file is generated from the WRITENACS',&
         ' subroutine in the file td_nacvs.F'
    IF (paral%io_parent)&
         WRITE (iunit,*)'Please note, that this is a "faked" output;'
    IF (paral%io_parent)&
         WRITE (iunit,*)'there are no intensities computed in CPMD.'
    IF (paral%io_parent)&
         WRITE (iunit,*)'Standard orientation:'
    IF (paral%io_parent)&
         WRITE (iunit,*)'---------------------------------------',&
         '------------------------------'
    IF (paral%io_parent)&
         WRITE (iunit,'(A,2(5X,A),14X,A)')&
         'Center','Atomic','Atomic','Coordinates (Angstroms)'
    IF (paral%io_parent)&
         WRITE (iunit,'(2(A,5X),1X,A,3X,3(11X,A))')&
         'Number','Number','Type','X','Y','Z'
    IF (paral%io_parent)&
         WRITE (iunit,*)'---------------------------------------',&
         '------------------------------'
    ! ---- make sure that the center of mass is the origin or move the atoms
    ! ---- from the center of the box to the origin and printout
    IF (lw) THEN
       DO i=1,ions1%nat
          xc(i)=(xc(i))/fbohr
          yc(i)=(yc(i))/fbohr
          zc(i)=(zc(i))/fbohr
          IF (paral%io_parent)&
               WRITE (15,22)  i, t(i), at, xc(i), yc(i), zc(i)
       ENDDO
    ELSE
       DO i=1,ions1%nat
          xc(i)=(xc(i)-bx/2)/fbohr
          yc(i)=(yc(i)-bx*by/2)/fbohr
          zc(i)=(zc(i)-bx*bz/2)/fbohr
          IF (paral%io_parent)&
               WRITE (15,22)  i, t(i), at, xc(i), yc(i), zc(i)
       END DO
    ENDIF
    ! ---- write some lines for molden -----------------------------
    iunit=15
    IF (paral%io_parent)&
         WRITE(iunit,*)'--------------------------------------------',&
         '-------------------------'
    IF (paral%io_parent)&
         WRITE(iunit,*)'      basis functions          primitive ',&
         'gaussians'
    IF (paral%io_parent)&
         WRITE(iunit,*)'      alpha electrons          beta electrons'
    IF (paral%io_parent)&
         WRITE(iunit,*)'********************************************',&
         '**************************'
    IF (paral%io_parent)&
         WRITE(iunit,*)
    IF (paral%io_parent)&
         WRITE(iunit,*)'Harmonic frequencies (cm**-1), IR intensities ',&
         '(KM/Mole),'
    IF (paral%io_parent)&
         WRITE(iunit,*)'Raman scattering activities (A**4/AMU), Raman ',&
         'depolarization ratios,'
    IF (paral%io_parent)&
         WRITE(iunit,*)'reduced masses (AMU), force constants ',&
         '(mDyne/A) and normal coordinates:'
    ! ---- write eigenvalues and eigenvectors in both files
    DO  ii=1,NINT(nroot*(nroot+1.0)/2.0),3
       i=ii+3
       IF (paral%io_parent)&
            WRITE(15,23) i-3, i-2, i-1
       IF (paral%io_parent)&
            WRITE(15,24) typ, typ, typ
       IF (paral%io_parent)&
            WRITE(15,25) fr, (vibe(l),l=i,i+2)
       DO n=1,5
          IF (paral%io_parent)&
               WRITE(15,25) rm(n), gz, gz, gz
       ENDDO
       IF (paral%io_parent)&
            WRITE(15,26) a,(ck(n),n=1,3),(ck(n),n=1,3),(ck(n),n=1,3)
       DO j=1,ions1%nat*3,3
          he=(j-1)/3+1
          IF (paral%io_parent)&
               WRITE(15,27) he,t(he),&
               (sder(j,m),sder(j+1,m),sder(j+2,m),m=i,i+2)
       END DO
    END DO
    IF (paral%io_parent)&
         WRITE(15,*) 'Normal termination of Gaussian 98.'
    IF (paral%io_parent)&
         CALL fileclose(15)
22  FORMAT(i5,i11,i14,4x,3(3x,f9.6))
23  FORMAT(i22,2i23)
24  FORMAT(20x,a2,2(21x,a2))
25  FORMAT(a17,f9.4,2f23.4)
    ! 26   FORMAT(a8,2X,3a17,2(2x,3a17))
26  FORMAT(a8,3a7,2(2x,3a7))
27  FORMAT(2i4,3(f19.8,2f17.8))
    GOTO 888
999 CONTINUE
    IF (paral%io_parent)&
         WRITE(6,*) 'COULD NOT OPEN FILE VIB[1,2].log! '
888 CONTINUE
    ! ==--------------------------------------------------------------==
  END SUBROUTINE writenacs
  ! ==================================================================

END MODULE td_nacvs_utils
