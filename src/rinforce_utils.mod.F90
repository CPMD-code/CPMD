MODULE rinforce_utils
  USE aainit_utils,                    ONLY: aainit
  USE aavan,                           ONLY: indv,&
                                             lpl
  USE atom,                            ONLY: gnl,&
                                             rps,&
                                             rw,&
                                             vr
  USE cnst,                            ONLY: pi
  USE cppt,                            ONLY: &
       gk, hg, qrad, rhops, tshel, tshels, twnl, vps, ylmb
  USE cvan,                            ONLY: deeq,&
                                             dvan,&
                                             nelev,&
                                             qq
  USE dpot,                            ONLY: dpot_mod
  USE eam,                             ONLY: tieam
  USE eam_pot_utils,                   ONLY: eamin
  USE ehpsi_utils,                     ONLY: eg2
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE fint,                            ONLY: fint1
  USE fitpack_utils,                   ONLY: curv2
  USE formf_utils,                     ONLY: formfn,&
                                             formfv,&
                                             formsg,&
                                             formup
  USE geq0mod,                         ONLY: geq0
  USE gvec,                            ONLY: gvec_com
  USE ions,                            ONLY: al,&
                                             ions0,&
                                             ions1,&
                                             rcl
  USE isos,                            ONLY: isos1
  USE kinds,                           ONLY: real_8
  USE kpclean_utils,                   ONLY: r_clean
  USE kpnt,                            ONLY: hgkm,&
                                             hgkp,&
                                             rk
  USE kpts,                            ONLY: kpts_com,&
                                             tkpts
  USE mm_input,                        ONLY: lqmmm
  USE mm_ion_dens,                     ONLY: mm_raggio,&
                                             mm_rhops
  USE mp_interface,                    ONLY: mp_sum
  USE nlcc,                            ONLY: corecg,&
                                             corel,&
                                             drhoc,&
                                             rcgrid,&
                                             rhoc
  USE nlccset_utils,                   ONLY: nlccset
  USE nlps,                            ONLY: nghcom,&
                                             nghtol,&
                                             nlps_com
  USE parac,                           ONLY: parai,&
                                             paral
  USE prmem_utils,                     ONLY: prmem
  USE pslo,                            ONLY: pslo_com
  USE qrada_s_utils,                   ONLY: qrada_s
  USE qspl,                            ONLY: &
       ggng, ggnh, nqdim, nsplpa, nsplpe, nsplpo, qspl1, qsrang, twns, voo
  USE ragg,                            ONLY: raggio
  USE rnlin_utils,                     ONLY: nlin
  USE rnlset_utils,                    ONLY: rnlset
  USE sgpp,                            ONLY: sgpp1,&
                                             sgpp2
  USE sphe,                            ONLY: gcutwmax
  USE str2,                            ONLY: qrada
  USE system,                          ONLY: &
       cnti, cntl, kpbeg, lx, maxsys, mx, nbrx, ncpw, nhx, nkpbl, nkpt, nlx, &
       parm
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: numcpus
  USE vdbinit_utils,                   ONLY: qinit,&
                                             vdbinit
  USE vdbp,                            ONLY: &
       betar, dion, ncpr1, qfunc, qqq, qrl, r, rab, rsatom, rscore, ru, rucore
  USE ylmr2_utils,                     ONLY: ylmr2
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rinforce
  PUBLIC :: testspline
  PUBLIC :: putps
  PUBLIC :: putwnl
  PUBLIC :: give_scr_putwnl
  PUBLIC :: calc_twnl
  PUBLIC :: setspline

CONTAINS

  ! ==================================================================
  SUBROUTINE rinforce
    ! ==--------------------------------------------------------------==
    ! Variables
    CHARACTER(*), PARAMETER                  :: procedureN = 'rinforce'

    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: ierr, im, is, istep, isub, &
                                                iv, jv, lmaxv, lp, lpmax, &
                                                lscr, lval, mrscr, nylmb
    REAL(real_8)                             :: aa, eself, pf, pub
    REAL(real_8), ALLOCATABLE                :: rs1(:), rs2(:), rs3(:)

! ==--------------------------------------------------------------==
! ==  ALLOCATION OF MEMORY                                        ==
! ==--------------------------------------------------------------==

    CALL tiset('  RINFORCE',isub)
    CALL setspline
    ALLOCATE(vps(maxsys%nsx,ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__) ! FIXME deallocate missing
    ALLOCATE(rhops(maxsys%nsx,ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__) ! FIXME deallocate missing
    IF (lqmmm%qmmm)  THEN
       ALLOCATE(mm_RHOPS(maxsys%nsx,ncpw%nhg),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    IF (pslo_com%tivan) THEN
       ALLOCATE(indv(nhx,maxsys%nsx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)! FIXME deallocate missing
       nqdim    = ncpw%nhgl
       IF (qspl1%qspline) nqdim=2*nsplpo
       ALLOCATE(qrad(nqdim,nbrx,nbrx,lx,maxsys%nsx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)! FIXME deallocate missing
    ELSE
       ! TODO what to do with INDV and QRAD?? 
    ENDIF
    IF (corel%tinlc) THEN
       ALLOCATE(rhoc(ncpw%nhg,maxsys%nsx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(rhoc)!,maxsys%nsx*nhg)
       IF (cntl%tpres.OR.cntl%tprcp) THEN
          ALLOCATE(drhoc(ncpw%nhg,maxsys%nsx,6),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          CALL zeroing(drhoc)!,6*maxsys%nsx*nhg)
       ENDIF
    ENDIF
    ! SCRATCH SPACE ALLOCATION
    CALL numcpus(parai%ncpus)
    mrscr   = MAX(8*maxsys%mmaxx,parai%ncpus*maxsys%mmaxx,nsplpo,nhx*nhx)
    ALLOCATE(rs1(mrscr),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(rs2(mrscr),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(rs3(mrscr),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! SPACE ALLOCATION FOR SPLINE
    ALLOCATE(voo(nsplpo,2,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(twns(nsplpo,2,maxsys%nhxs,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(voo)!,2*nsplpo*maxsys%nsx)
    CALL zeroing(twns)!,2*nsplpo*maxsys%nsx*maxsys%nhxs)
    ! INITIALIZATION FOR VANDERBILT PSEUDOPOTENTIALS
    IF (pslo_com%tivan) THEN
       lmaxv=0
       DO is=1,ions1%nsp
          IF (pslo_com%tvan(is)) THEN
             lval=ncpr1%nvales(is)
             IF (lval.GT.3) CALL stopgm('RINFORCE',' L>3, L=',& 
                  __LINE__,__FILE__)
             IF (lval.GT.lmaxv) lmaxv=lval
          ENDIF
       ENDDO
       ! INITIALIZE CLEBSCH-GORDON COEFFICIENTS
       CALL aainit(lmaxv)
       ALLOCATE(nelev(maxsys%nsx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(qq(maxsys%nhxs,maxsys%nhxs,maxsys%nsx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(dvan(maxsys%nhxs,maxsys%nhxs,maxsys%nsx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(deeq(ions1%nat,maxsys%nhxs,maxsys%nhxs,2),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(deeq)
    ENDIF
    ! ==--------------------------------------------------------------==
    DO is=1,ions1%nsp
       IF (pslo_com%tvan(is)) THEN
          CALL formfv(is,rs1,rs2)
       ELSEIF (sgpp1%tsgp(is)) THEN
          IF (pslo_com%tnum(is)) THEN
             CALL formup(is,rs1,rs2)
          ELSE
             CALL formsg(is,rs1,rs2)
          ENDIF
       ELSE
          IF (ions0%igau(is).GE.0) CALL formfn(is,rs1,rs2)
       ENDIF
    ENDDO
    CALL putps
    DO is=1,ions1%nsp
       IF (pslo_com%tvan(is)) THEN
          CALL vdbinit(is,rs1,rs2)
       ELSE
          CALL rnlset(is,rs1,rs2,rs3)
          CALL nlin(is,rs1,rs2)
       ENDIF
    ENDDO
    IF (pslo_com%tivan) CALL qinit(rs1,rs2)
    pub = 0._real_8
    IF (geq0) THEN
       !$omp parallel do private(IS) reduction(+:PUB)
#ifdef __SR8000
       !poption parallel, tlocal(IS), psum(PUB)
#endif
       DO is=1,ions1%nsp
          pub = pub + REAL(ions0%na(is),kind=real_8)*vps(is,1)
       ENDDO
    ENDIF
    aa = pub
    CALL mp_sum(aa, pub, parai%allgrp)

    ! ==--------------------------------------------------------------==
    ! ==  SELF-ENERGY (NOTHING TO DO WITH MANY-BODY SELF-ENERGY)      ==
    ! ==--------------------------------------------------------------==
    ener_com%eself=0._real_8
    eself=0.0_real_8
    !$omp parallel do private(IS) reduction(+:ESELF)
#ifdef __SR8000
    !poption parallel, tlocal(IS), psum(ESELF)
#endif
    DO is=1,ions1%nsp
       eself=eself+REAL(ions0%na(is),kind=real_8)*ions0%zv(is)*ions0%zv(is)/raggio(is)
    ENDDO
    ener_com%eself=eself
    ener_com%eself=ener_com%eself/SQRT(2._real_8*pi)
    ! Correction for charged systems
    IF (.NOT.isos1%tclust) THEN
       pf=crge%charge
       ener_com%eself=ener_com%eself-pf*pub
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Calculation of Maximal value of L
    maxsys%lpmax=0
    IF (pslo_com%tivan) THEN
       lpmax=0
       !$omp parallel do private(IM,IV,JV,LP) reduction(max:lpmax)
       DO im=1,mx
          DO iv=1,nlx
             DO jv=iv,nlx
                lp=lpl(iv,jv,im)
                lpmax=MAX(lpmax,lp)
             ENDDO
          ENDDO
       ENDDO
       maxsys%lpmax=lpmax
    ENDIF
    DO is=1,ions1%nsp
       IF (nlps_com%ngh(is).GT.0) THEN
          IF (dpot_mod%tkb(is)) THEN
             iv=nlps_com%ngh(is)
             maxsys%lpmax=MAX(maxsys%lpmax,nghcom(iv,is))
          ELSEIF (sgpp1%tsgp(is)) THEN
             iv=nlps_com%ngh(is)
             maxsys%lpmax=MAX(maxsys%lpmax,sgpp2%lpval(iv,is))
          ELSEIF (.NOT.pslo_com%tvan(is)) THEN
             istep=NINT(nlps_com%rmaxn(is))
             DO iv=1,nlps_com%ngh(is)
                lp=(iv-1)/istep+1
                maxsys%lpmax=MAX(maxsys%lpmax,lp)
             ENDDO
          ENDIF
       ENDIF
    ENDDO
    IF (pslo_com%tivan) THEN
       nylmb=nkpt%nhgk*maxsys%lpmax*nkpt%nkpnt
       ALLOCATE(ylmb(nkpt%nhgk,maxsys%lpmax,nkpt%nkpnt),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(ylmb)!,nylmb)
    ELSE
       nylmb=0
    ENDIF
    ALLOCATE(twnl(nkpt%ngwk,maxsys%nhxs,maxsys%nsx,kpts_com%nkptall),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__) ! FIXME deallocate missing

    CALL zeroing(twnl)!,maxsys%nsx*maxsys%nhxs*nkpt%ngwk*kpts_com%nkptall)
    CALL give_scr_putwnl(lscr,tag)
    CALL putwnl
    ! ==--------------------------------------------------------------==
    ! ==  NONLINEAR CORE CORRECTION                                   ==
    ! ==--------------------------------------------------------------==
    IF (corel%tinlc) CALL nlccset
    ! ==--------------------------------------------------------------==
    ! ==  FINITE TEMPERATURE ELECTRONS                                ==
    ! ==--------------------------------------------------------------==
    IF (fint1%ttrot) CALL eg2
    ! ==--------------------------------------------------------------==
    ! ==  EAM POTENTIALS                                              ==
    ! ==--------------------------------------------------------------==
    tieam=.FALSE.
    DO is=1,ions1%nsp
       tieam=tieam.OR.dpot_mod%team(is)
    ENDDO
    IF (tieam) CALL eamin
    ! ==--------------------------------------------------------------==
    ! ==  INITIALIZE DERIVATIVE ARRAYS FOR STRESS CALCULATION         ==
    ! ==--------------------------------------------------------------==
    IF (cntl%tprcp.OR.cntl%tpres) THEN
       IF (pslo_com%tivan) THEN
          nqdim    = ncpw%nhgl
          IF (qspl1%qspline) nqdim=2*nsplpo
          ALLOCATE(qrada(nqdim,nbrx,nbrx,lx,maxsys%nsx),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ELSE
          ! TODO what to do with QRADA
       ENDIF
       IF (pslo_com%tivan) CALL qrada_s(rs1,rs2)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ==  DEALLOCATE ALL MEMORY NOT LONGER NEEDED                     ==
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent) CALL prmem('  RINFORCE')
    DEALLOCATE(rs1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rs2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rs3,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rcgrid,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(corecg,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(al,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! FIXME vw this free causes a crash with gcc-mpi-4 
    ! call free90(bl)
    DEALLOCATE(gnl,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rcl,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rps,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rw,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(vr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (pslo_com%tivan) THEN
       DEALLOCATE(rscore,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(dion,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(betar,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(qqq,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(qfunc,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(qrl,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(rucore,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(ru,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(rsatom,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(r,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(rab,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    CALL tihalt('  RINFORCE',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rinforce
  ! ==================================================================
  SUBROUTINE putps
    ! ==--------------------------------------------------------------==
    ! == RHOPS: Calculates the Gaussian charge distributions which    ==
    ! == replaced the ionic point charges (smeared charge density)    ==
    ! == in G space(see ener.inc and cppt.inc)                        ==
    ! == VPS: local pseudopotential per species in G space            ==
    ! ==--------------------------------------------------------------==
    ! Variables
    CHARACTER(*), PARAMETER                  :: procedureN = 'putps'

    INTEGER                                  :: ig, is, isub
    REAL(real_8)                             :: emax, mm_R2MAX, qmax, r2max, &
                                                vol

! ==--------------------------------------------------------------==

    CALL tiset('     PUTPS',isub)
    vol=1._real_8/parm%omega
    DO is=1,ions1%nsp
       r2max=raggio(is)*raggio(is)
       !$omp parallel do private(IG,QMAX,EMAX)
#ifdef __SR8000
       !poption parallel, tlocal(IG,QMAX,EMAX)
#endif 
       DO ig=1,ncpw%nhg
          qmax=0.25_real_8*r2max*hg(ig)*parm%tpiba2
          emax=EXP(-qmax)
          rhops(is,ig)=-ions0%zv(is)*emax*vol
          vps(is,ig)=vol*curv2(hg(ig),nsplpo,ggnh(1),voo(1,1,is),&
               voo(1,2,is),0.0_real_8)
       ENDDO
       IF (lqmmm%qmmm) THEN
          mm_R2MAX=mm_RAGGIO(is)*mm_raggio(is)
          !$omp parallel do private(IG,QMAX,EMAX)
#ifdef __SR8000
          !poption parallel, tlocal(IG,QMAX,EMAX)
#endif 
          DO ig=1,ncpw%nhg
             qmax=0.25_real_8*mm_R2MAX*hg(ig)*parm%tpiba2
             emax=EXP(-qmax)
             mm_RHOPS(is,ig)=-ions0%zv(is)*emax*vol
          ENDDO
       ENDIF
       IF (geq0) rhops(is,1)=-ions0%zv(is)*vol
       IF (geq0.AND.lqmmm%qmmm) mm_RHOPS(is,1)=-ions0%zv(is)*vol
    ENDDO
    CALL tihalt('     PUTPS',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE putps
  ! ==================================================================
  SUBROUTINE putwnl
    ! ==--------------------------------------------------------------==
    ! == CALCULATE TWNL(1:NGW,1:NGH(IS),1:NSP,1:NKPNT) [cppt.inc]     ==
    ! ==        Non-Local projectors array                            ==
    ! ==        for each G-components (Kleinman-Bylander form)        ==
    ! == FOR TIVAN calculate also YLMB                                ==
    ! ==--------------------------------------------------------------==
    ! Variables
    CHARACTER(*), PARAMETER                  :: procedureN = 'putwnl'

    INTEGER                                  :: ierr, ig, ik, ikk, ikpt, &
                                                ikylmb, is, istep, isub, iv, &
                                                kbeg, kend, kinc, l, lp
    REAL(real_8)                             :: tw, vol, xx
    REAL(real_8), ALLOCATABLE                :: gkrk(:,:)
    REAL(real_8), EXTERNAL                   :: dasum

    CALL tiset('    PUTWNL',isub)
    IF (pslo_com%tivan) THEN
       ! YLMB is already allocated in RINFORCE
    ELSE
       ALLOCATE(ylmb(nkpt%nhgk, MAX(maxsys%lpmax,1), 1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       IF (tkpts%tkpnt) THEN
          ALLOCATE(gkrk(3, ncpw%nhg),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    vol=1._real_8/SQRT(parm%omega)
    CALL inq_swap(kbeg,kend,kinc)
    DO ikpt=kbeg,kend,kinc
       IF (tkpts%tkblock) CALL rkpt_swap(gkrk,1,ikpt,'HGKP HGKM MASKGW')
       DO ik=1,nkpbl(ikpt)
          ikk=kpbeg(ikpt)+ik
          IF (tkpts%tkpnt) THEN
             IF (pslo_com%tivan) THEN
                ikylmb=ik
             ELSE
                ikylmb=1
             ENDIF
             DO lp=1,maxsys%lpmax
                !$omp parallel do private(IG)
#ifdef __SR8000
                !poption parallel
#endif
                DO ig = 1 , ncpw%nhg
                   gkrk(1,ig)=gk(1,ig)+rk(1,ikk)
                   gkrk(2,ig)=gk(2,ig)+rk(2,ikk)
                   gkrk(3,ig)=gk(3,ig)+rk(3,ikk)
                ENDDO

                CALL ylmr2(lp,ncpw%nhg,hgkp(1,ik),gkrk,ylmb(1,lp,ikylmb))

                !$omp parallel do private(IG)
#ifdef __SR8000
                !poption parallel
#endif
                DO ig = 1 , ncpw%nhg
                   gkrk(1,ig)=-gk(1,ig)+rk(1,ikk)
                   gkrk(2,ig)=-gk(2,ig)+rk(2,ikk)
                   gkrk(3,ig)=-gk(3,ig)+rk(3,ikk)
                ENDDO

                CALL ylmr2(lp,ncpw%nhg,hgkm(1,ik),gkrk,ylmb(1+ncpw%nhg,lp,ikylmb))

             ENDDO
          ELSE
             DO lp=1,maxsys%lpmax
                CALL ylmr2(lp,ncpw%nhg,hg,gk,ylmb(1,lp,1))
             ENDDO
          ENDIF
          DO is=1,ions1%nsp
             ! Shell structure code disabled
             tshel(is)=.FALSE.
             tshels=.FALSE.
             DO iv=1,nlps_com%ngh(is)
                IF (pslo_com%tvan(is)) THEN
                   istep=ncpr1%nvales(is)*ncpr1%nvales(is)
                   l=nghtol(iv,is)+1
                   lp=1+MOD(iv-1,istep)
                ELSEIF (dpot_mod%tkb(is)) THEN
                   l=nghtol(iv,is)+1
                   lp=nghcom(iv,is)
                ELSEIF (sgpp1%tsgp(is)) THEN
                   lp=sgpp2%lpval(iv,is)
                ELSE
                   istep=NINT(nlps_com%rmaxn(is))
                   l=nghtol(iv,is)+1
                   lp=(iv-1)/istep+1
                ENDIF
                xx=dasum(nsplpo,twns(1,1,iv,is),1)
                IF (xx.GT.1.e-12_real_8) THEN
                   IF (tkpts%tkpnt) THEN
                      DO ig=1,ncpw%ngw
                         tw=curv2(hgkp(ig,ik),nsplpo,ggng(1),twns(1,1,iv,is),&
                              twns(1,2,iv,is),0.0_real_8)
                         twnl(ig,iv,is,ik)=ylmb(ig,lp,ikylmb)*tw*vol
                      ENDDO
                      DO ig=1,ncpw%ngw
                         tw=curv2(hgkm(ig,ik),nsplpo,ggng(1),twns(1,1,iv,is),&
                              twns(1,2,iv,is),0.0_real_8)
                         twnl(ig+ncpw%ngw,iv,is,ik)=ylmb(ig+ncpw%nhg,lp,ikylmb)*tw*vol
                      ENDDO
                      CALL r_clean(twnl(1,iv,is,ik),1,ik)
                   ELSE
                      DO ig=1,ncpw%ngw
                         tw=curv2(hg(ig),nsplpo,ggng(1),twns(1,1,iv,is),&
                              twns(1,2,iv,is),0.0_real_8)
                         twnl(ig,iv,is,1)=ylmb(ig,lp,1)*tw*vol
                      ENDDO
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDDO                 ! END OF 1,NKPNT
       IF (tkpts%tkblock) THEN
          IF (pslo_com%tivan) THEN
             CALL wkpt_swap(gkrk,1,ikpt,'TWNL YLMB')
          ELSE
             CALL wkpt_swap(gkrk,1,ikpt,'TWNL')
          ENDIF
       ENDIF
    ENDDO                     ! END OF NBLKP
    ! ==--------------------------------------------------------------==
    IF (.NOT.pslo_com%tivan) THEN
       DEALLOCATE(ylmb,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       IF (tkpts%tkpnt) THEN
          DEALLOCATE(gkrk,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
               __LINE__,__FILE__)
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL tihalt('    PUTWNL',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE putwnl
  ! ==================================================================
  SUBROUTINE calc_twnl(ikpt)
    ! ==--------------------------------------------------------------==
    ! == CALCULATE TWNL(1:NGW,1:NGH(IS),1:NSP,1:NKPNT) [cppt.inc]     ==
    ! ==        Non-Local projectors array                            ==
    ! ==        for each G-components (Kleinman-Bylander form)        ==
    ! == FOR TIVAN calculate also YLMB                                ==
    ! ==--------------------------------------------------------------==
    ! == IF TIVAN=.TRUE. YLMB is used (need allocation)               ==
    ! == otherwise  allocate temporarily YLMB and GKRK                ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ikpt

    CHARACTER(*), PARAMETER                  :: procedureN = 'calc_twnl'

    INTEGER                                  :: ierr, ig, ik, ikk, ikylmb, &
                                                is, istep, isub, iv, l, lp
    REAL(real_8)                             :: tw, vol, xx
    REAL(real_8), ALLOCATABLE                :: gkrk(:,:)
    REAL(real_8), EXTERNAL                   :: dasum

    CALL tiset(' CALC_TWNL',isub)
    IF (pslo_com%tivan) THEN
       ! YLMB is already allocated and is used in this routine.
       IF (tkpts%tkpnt)  THEN
          ALLOCATE(gkrk(3,ncpw%nhg),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
    ELSE
       ALLOCATE(ylmb(nkpt%nhgk, MAX(maxsys%lpmax,1),1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       !call memory90(ylmb, (/ nhgk, max(maxsys%lpmax,1),1/), 'YLMB')
       IF (tkpts%tkpnt)  THEN
          ALLOCATE(gkrk(3,ncpw%nhg),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    vol=1._real_8/SQRT(parm%omega)
    DO ik=1,nkpbl(ikpt)
       ikk=kpbeg(ikpt)+ik
       IF (tkpts%tkpnt) THEN
          IF (pslo_com%tivan) THEN
             ikylmb=ik
          ELSE
             ikylmb=1
          ENDIF
          DO lp=1,maxsys%lpmax

             !$omp parallel do private(IG)
#ifdef __SR8000
             !poption parallel
#endif
             DO ig = 1 , ncpw%nhg
                gkrk(1,ig)=gk(1,ig)+rk(1,ikk)
                gkrk(2,ig)=gk(2,ig)+rk(2,ikk)
                gkrk(3,ig)=gk(3,ig)+rk(3,ikk)
             ENDDO

             CALL ylmr2(lp,ncpw%nhg,hgkp(1,ik),gkrk,ylmb(1,lp,ikylmb))

             !$omp parallel do private(IG)
#ifdef __SR8000
             !poption parallel
#endif
             DO ig = 1 , ncpw%nhg
                gkrk(1,ig)=-gk(1,ig)+rk(1,ikk)
                gkrk(2,ig)=-gk(2,ig)+rk(2,ikk)
                gkrk(3,ig)=-gk(3,ig)+rk(3,ikk)
             ENDDO

             CALL ylmr2(lp,ncpw%nhg,hgkm(1,ik),gkrk,ylmb(1+ncpw%nhg,lp,ikylmb))
          ENDDO
       ELSE
          DO lp=1,maxsys%lpmax
             CALL ylmr2(lp,ncpw%nhg,hg,gk,ylmb(1,lp,1))
          ENDDO
       ENDIF
       DO is=1,ions1%nsp
          ! Shell structure code disabled
          tshel(is)=.FALSE.
          tshels=.FALSE.
          DO iv=1,nlps_com%ngh(is)
             IF (pslo_com%tvan(is)) THEN
                istep=ncpr1%nvales(is)*ncpr1%nvales(is)
                l=nghtol(iv,is)+1
                lp=1+MOD(iv-1,istep)
             ELSEIF (dpot_mod%tkb(is)) THEN
                l=nghtol(iv,is)+1
                lp=nghcom(iv,is)
             ELSEIF (sgpp1%tsgp(is)) THEN
                lp=sgpp2%lpval(iv,is)
             ELSE
                istep=NINT(nlps_com%rmaxn(is))
                l=nghtol(iv,is)+1
                lp=(iv-1)/istep+1
             ENDIF
             xx=dasum(nsplpo,twns(1,1,iv,is),1)
             IF (xx.GT.1.e-12_real_8) THEN
                IF (tkpts%tkpnt) THEN
                   DO ig=1,ncpw%ngw
                      tw=curv2(hgkp(ig,ik),nsplpo,ggng(1),twns(1,1,iv,is),&
                           twns(1,2,iv,is),0.0_real_8)
                      twnl(ig,iv,is,ik)=ylmb(ig,lp,ikylmb)*tw*vol
                   ENDDO
                   DO ig=1,ncpw%ngw
                      tw=curv2(hgkm(ig,ik),nsplpo,ggng(1),twns(1,1,iv,is),&
                           twns(1,2,iv,is),0.0_real_8)
                      twnl(ig+ncpw%ngw,iv,is,ik)=ylmb(ig+ncpw%nhg,lp,ikylmb)*tw*vol
                   ENDDO
                   CALL r_clean(twnl(1,iv,is,ik),1,ik)
                ELSE
                   DO ig=1,ncpw%ngw
                      tw=curv2(hg(ig),nsplpo,ggng(1),twns(1,1,iv,is),&
                           twns(1,2,iv,is),0.0_real_8)
                      twnl(ig,iv,is,1)=ylmb(ig,lp,1)*tw*vol
                   ENDDO
                ENDIF
             ENDIF
          ENDDO
       ENDDO
    ENDDO                     ! END OF 1,NKPNT
    ! ==--------------------------------------------------------------==
    IF (pslo_com%tivan) THEN
       ! YLMB is already allocated and is used in this routine.
       IF (tkpts%tkpnt) DEALLOCATE(gkrk,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ELSE
       DEALLOCATE(ylmb,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       IF (tkpts%tkpnt) DEALLOCATE(gkrk,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    CALL tihalt(' CALC_TWNL',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE calc_twnl
  ! ==================================================================
  SUBROUTINE give_scr_putwnl(lputwnl,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lputwnl
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: lylmb

    IF (tkpts%tkpnt) THEN
       lputwnl=3*ncpw%nhg
    ELSE
       lputwnl=0
    ENDIF
    IF (pslo_com%tivan) THEN
       lylmb=0
    ELSE
       lylmb=nkpt%nhgk*maxsys%lpmax
    ENDIF
    lputwnl=MAX(1,lputwnl+lylmb)
    tag='3*NHG+NHGK*maxsys%lpmax'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_putwnl
  ! ==================================================================
  SUBROUTINE setspline
    ! ==--------------------------------------------------------------==
    ! Variables

    CHARACTER(*), PARAMETER                  :: procedureN = 'setspline'
    INTEGER, PARAMETER                       :: nadd = 10 

    INTEGER                                  :: i, ierr, il, nn
    REAL(real_8)                             :: dgl, xsaim, xsnow

! ==--------------------------------------------------------------==

    nsplpo=cnti%nsplp
    IF (nsplpo.GT.maxsys%mmaxx-nadd) THEN
       IF (paral%io_parent)&
            WRITE(6,*) ' SETSPLINE| MAXIMUM NUMBER OF SPLINE POINTS IS ',&
            maxsys%mmaxx-nadd
       CALL stopgm('SETSPLINE','TOO MANY SPLINE POINTS REQUESTED',& 
            __LINE__,__FILE__)
    ENDIF
    nn=nsplpo+nadd+10
    ALLOCATE(ggnh(nn),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(ggng(nn),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    dgl=gvec_com%gcut/REAL(nsplpo-1,kind=real_8)
    DO il=1,nsplpo
       ggnh(il)=(il-1)*dgl
    ENDDO
    IF (qsrang.GT.1) THEN
       dgl=(qsrang-1._real_8)*gvec_com%gcut/REAL(nadd,kind=real_8)
       DO il=nsplpo+1,nsplpo+nadd
          ggnh(il)=ggnh(il-1)+dgl
       ENDDO
    ENDIF
    dgl=(gcutwmax)/REAL(nsplpo-1,kind=real_8)
    DO il=1,nsplpo
       ggng(il)=(il-1)*dgl
    ENDDO
    IF (qsrang.GT.1) THEN
       dgl=(qsrang-1._real_8)*(gcutwmax)/REAL(nadd,kind=real_8)
       DO il=nsplpo+1,nsplpo+nadd
          ggng(il)=ggng(il-1)+dgl
       ENDDO
    ENDIF
    IF (qsrang.GT.1) nsplpo=nsplpo+nadd
    ! ==--------------------------------------------------------------==
    ! DISTRIBUTE SPLINE POINTS
    xsnow=0._real_8
    DO i=parai%nproc,1,-1
       xsaim = xsnow + REAL(nsplpo,kind=real_8)/REAL(parai%nproc,kind=real_8)
       IF (i.EQ.parai%mepos+1) THEN
          nsplpa=NINT(xsnow)+1
          nsplpe=NINT(xsaim)
          IF (NINT(xsaim).GT.nsplpo) nsplpe=nsplpo
          IF (i.EQ.1) nsplpe=nsplpo
       ENDIF
       xsnow = xsaim
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE setspline
  ! ==================================================================
  SUBROUTINE testspline
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: warn

! ==--------------------------------------------------------------==

    warn=0._real_8
    IF (ggnh(nsplpo).LT.hg(ncpw%nhg)) warn=1._real_8
    CALL mp_sum(warn,parai%allgrp)
    IF (warn.GT.0.5_real_8.AND.paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,'(A,A)') ' WARNING| SPLINE REGION SMALLER THAN',&
            ' MAXIMUM G-VALUE'
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE testspline
  ! ==================================================================

END MODULE rinforce_utils
