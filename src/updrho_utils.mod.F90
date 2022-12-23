MODULE updrho_utils
  USE adjmu_utils,                     ONLY: adjmu,&
                                             enert,&
                                             occup
  USE andp,                            ONLY: rin0,&
                                             rmix,&
                                             rout0
  USE bogol_utils,                     ONLY: bogol,&
                                             give_scr_bogol
  USE broy,                            ONLY: broy1
  USE conv,                            ONLY: nac,&
                                             nbc
  USE davidson_utils,                  ONLY: davidson,&
                                             give_scr_davidson
  USE dist_friesner_utils,             ONLY: dist_friesner,&
                                             give_scr_dist_friesner
  USE ehpsi_utils,                     ONLY: change_b2l,&
                                             change_betap,&
                                             epot
  USE elct,                            ONLY: crge
  USE elct2,                           ONLY: tfixo
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE fint,                            ONLY: fint1,&
                                             fint4
  USE friesner_c_p_utils,              ONLY: frie_c_p
  USE friesner_c_utils,                ONLY: friesner_c
  USE friesner_utils,                  ONLY: friesner,&
                                             give_scr_friesner
  USE frsblk_c_utils,                  ONLY: frsblk_c,&
                                             give_scr_frsblk_c
  USE frsblk_utils,                    ONLY: frsblk,&
                                             give_scr_frsblk
  USE func,                            ONLY: func1
  USE hfx_drivers,                     ONLY: hfx
  USE hpsi_utils,                      ONLY: give_scr_hpsi,&
                                             hpsi
  USE hubbardu,                        ONLY: c2u0,hubbu
  USE hubbardu_utils,                  ONLY: hubbardUcorrection,&
                                             add_hubbardu
  USE k_diis_rhofix_utils,             ONLY: k_diis_rhofix
  USE k_forces_utils,                  ONLY: give_scr_kforces
  USE kdp,                             ONLY: &
       akdp, akdp2, auxdiag, bkdp, brkdp, ckdp, ekdp, fkdp, pkdp, rauxdiag, &
       wkdp, xkdp, xlambda, xlambda0
  USE kdp_diag_utils,                  ONLY: kdp_diag
  USE kdp_prep_utils,                  ONLY: kdp_prep
  USE kdp_rho_utils,                   ONLY: kdp_rho
  USE kdpc,                            ONLY: bmix,&
                                             nkdp,&
                                             tkdp
  USE kinds,                           ONLY: real_8
  USE kpnt,                            ONLY: wk
  USE kpts,                            ONLY: tkpts
  USE ksdiag_utils,                    ONLY: ksdiag
  USE mixing_g_utils,                  ONLY: give_scr_mixing,&
                                             mixing_g
  USE mixing_r_utils,                  ONLY: mixing_r
  USE mp_interface,                    ONLY: mp_max,&
                                             mp_sum
  USE norm,                            ONLY: cnorm,&
                                             gemax
  USE parac,                           ONLY: parai,&
                                             paral
  USE pcgrad_driver,                   ONLY: eback
  USE pslo,                            ONLY: pslo_com
  USE ptheory_utils,                   ONLY: give_scr_ptheory
  USE rhoofr_c_utils,                  ONLY: rhoofr_c
  USE rhoofr_utils,                    ONLY: give_scr_rhoofr,&
                                             rhoofr
  USE rnlfor_utils,                    ONLY: rnlfor
  USE rnlrh_utils,                     ONLY: rnlrh
  USE rnlsm_utils,                     ONLY: give_scr_rnlsm,&
                                             rnlsm
  USE ropt,                            ONLY: iteropt,&
                                             ropt_mod
  USE rpiiint_utils,                   ONLY: rpiiint
  USE soft,                            ONLY: soft_com
  USE spin,                            ONLY: amudo,&
                                             amuup,&
                                             clsd,&
                                             spin_mod,&
                                             tdsp1
  USE stress_utils,                    ONLY: give_scr_stress,&
                                             stress
  USE symtrz_utils,                    ONLY: give_scr_symvec,&
                                             symvec
  USE system,                          ONLY: &
       cnti, cntl, cntr, fpar, kpbeg, maxsys, ncpw, nkpbl, nkpt, parm
  USE tauf,                            ONLY: itaur
  USE testex_utils,                    ONLY: testex
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: dgive
  USE vbeta_utils,                     ONLY: vbeta
  USE vdw_utils,                       ONLY: vdw,&
                                             vdw_wf
  USE vdwcmod,                         ONLY: &
       empvdwi, empvdwr, idvdw, ivdw, jvdw, vdwbe, vdwi, vdwl, vdwr, vdwrm, &
       vdwst
  USE vofrho_utils,                    ONLY: give_scr_vofrho,&
                                             vofrho
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: updrho
  PUBLIC :: give_scr_updrho

CONTAINS

  ! ==================================================================
  SUBROUTINE updrho(c0,c2,cr,sc0,cscr,vpp,tau0,fion,eigv,rhoe,&
       psi,nstate,tfor,tinfo,tstress,nfr,thl,nhpsi)
    ! ==--------------------------------------------------------------==
    ! ==               UPDATES THE DENSITY                            ==
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8)                          :: cr(*), sc0(*), cscr(*)
    REAL(real_8)                             :: vpp(:), tau0(:,:,:), &
                                                fion(:,:,:), rhoe(:,:)
    COMPLEX(real_8)                          :: psi(:,:)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: eigv(nstate,nkpt%nkpts)
    COMPLEX(real_8) :: c2(nkpt%ngwk,nstate), c0(nkpt%ngwk,nstate,nkpt%nkpnt)
    LOGICAL                                  :: tfor, tinfo, tstress
    INTEGER                                  :: nfr
    REAL(real_8)                             :: thl(2)
    INTEGER                                  :: nhpsi

    CHARACTER(*), PARAMETER                  :: procedureN = 'updrho'
    REAL(real_8), PARAMETER                  :: betaelmax = 1.e33_real_8 , &
                                                tollpar = 1.e-6_real_8

    COMPLEX(real_8), ALLOCATABLE             :: cgs(:,:)
    INTEGER :: i, ib, idamax, ierr, ig, ik, ikind, ikk, ikpt, imax, isub, &
      kbeg, kend, kinc, nconv, nkpoint, nnx, nnxs, noccup, nwfc
    INTEGER, SAVE                            :: ifirst = 0, netot
    LOGICAL                                  :: tadjmu, trefine
    REAL(real_8)                             :: detot, eband1, eband2, eeig1, &
                                                eeig2, ehfx, entropy1, &
                                                entropy2, vhfx
    REAL(real_8), ALLOCATABLE                :: edav(:), focc(:,:), &
                                                rhot(:,:), we(:,:)
    REAL(real_8), EXTERNAL                   :: dasum
    REAL(real_8), SAVE                       :: ddrho, etot0

    CALL tiset(procedureN,isub)
    IF (cntl%tlsd .OR. func1%mhfx.NE.0) THEN
       ALLOCATE(cgs(nkpt%ngwk,nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ELSE
       ! avoid 'not allocated' fortran runtime error
       ALLOCATE(cgs(1,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
    ENDIF
    IF (cntl%tlsd) THEN
       ALLOCATE(rhot(fpar%nnr1,clsd%nlsd),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (ifirst.EQ.0) THEN
       ddrho = 1.0_real_8
       ifirst = 1
       etot0=0.0_real_8
       netot=0
    ELSE IF (nfr.EQ.1) THEN
       ddrho = ddrho*1.5e+01_real_8! Change to looser wftol for first loop
       etot0=0.0_real_8
       netot=0
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == ALLOCATE ARRAYS FOR KDP ROUTINES                             ==
    ! ==--------------------------------------------------------------==
    IF (tkdp.AND.ifirst.EQ.0) THEN
       ifirst=1
       ALLOCATE(pkdp(3,nstate,nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(ekdp(nstate,nkdp),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(fkdp(nstate,nkdp),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(xlambda0(nstate,nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(ckdp(nkpt%ngwk,nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(akdp(nstate,nstate,nkdp),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(bkdp(nstate,nstate,nkdp),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(brkdp(nstate,nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(xlambda(nstate,nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(akdp2(nstate,nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(auxdiag(2*nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(rauxdiag(3*nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    nnx=fpar%nnr1*clsd%nlsd
    ! Copy density to scratch space used in VOFRHO
    CALL dcopy(nnx,rin0(1,1),1,rhoe,1)
    CALL zeroing(fion)!,3*maxsys%nax*maxsys%nsx)
    CALL rpiiint(ener_com%esr,tau0,fion,iteropt%iesr,tfor)
    ! ==--------------------------------------------------------------==
    IF (.NOT.cntl%tlanc)  THEN
       ALLOCATE(edav(cnti%ndavv),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    IF (tkpts%tkpnt.AND. cntl%tdavi) CALL stopgm("UPDRHO",&
         "DAVIDSON AND K-POINTS NOT IMPLEMENTED",& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! == Calculate the Stress Tensor                                  ==
    ! ==--------------------------------------------------------------==
    IF (tstress) CALL stress(c0,tau0,crge%f,psi(:,1),nstate)
    ! Calculate Potential
    IF (tfor) CALL eback(0)
    CALL vofrho(tau0,fion,rhoe,psi,tfor,tstress)
    IF (tfor) CALL eback(1)
    ! ==--------------------------------------------------------------==
    ! == Force on the ions from nonlocal pseudopotential              ==
    ! ==--------------------------------------------------------------==
    IF (tfor) THEN
       CALL inq_swap(kbeg,kend,kinc)
       DO ikpt=kbeg,kend,kinc
          nkpoint=nkpbl(ikpt)
          IF (tkpts%tkblock) CALL rkpt_swap(c0,nstate,ikpt,&
               'HGKP HGKM MASKGW EIGKR TWNL C0')
          DO ikind=1,nkpoint
             CALL rnlsm(c0(:,:,ikind),nstate,&
                  ikpt,ikind,.TRUE.)
             IF (.NOT.fint1%tbogo) THEN
                CALL hpsi(c0(:,:,ikind),c2,sc0,&
                     rhoe,psi(:,1),nstate,ikind,clsd%nlsd)
                ! ,.TRUE.,TSTRESS
             ENDIF
          ENDDO
          ! ENL: Non-local PP energy
          CALL rnlrh(ener_com%enl,nstate,nkpt%nkpts)
          CALL rnlfor(fion,crge%f(1,kpbeg(ikpt)+1),wk(kpbeg(ikpt)+1),&
               nstate,nkpoint)
          hubbu%ehub=0.D0
          if(cntl%thubb) then 
            if(.not.allocated(c2u0)) allocate(c2u0(ncpw%ngw,nstate))
            call hubbardUcorrection(c0(:,:,1),c2u0,tau0,fion,nstate,psi(:,1),tfor,0)
            call mp_sum(hubbu%ehub,parai%allgrp)
            call add_hubbardu(c2,c2u0,nstate)
          endif
       ENDDO
       CALL mp_sum(fion,3*maxsys%nax*maxsys%nsx,parai%allgrp)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == Energy and force on the ions from van der Waals interaction  ==
    ! ==--------------------------------------------------------------==
    vdwr%evdw=0._real_8
    CALL zeroing(vdwr%devdw)!,6)
    IF (vdwl%vdwc) CALL vdw(tau0,empvdwi%nvdw,idvdw,ivdw,jvdw,vdwst,vdwrm,vdwbe,&
         empvdwr%vdweps,empvdwr%s6grim,vdwi%nxvdw,vdwi%nyvdw,vdwi%nzvdw,vdwr%evdw,fion,vdwr%devdw)
    ! ==--------------------------------------------------------------==
    ! == Energy and force on the ions for Wannier-based van der Waals ==
    ! ==--------------------------------------------------------------==
    IF (vdwl%vdwd) THEN
       vdwr%evdw=0._real_8
       CALL zeroing(vdwr%devdw)!,6)
       nwfc=nstate
       IF (tfor) CALL vdw_wf(tau0,fion,nwfc)
    ENDIF
    ! Symmetrize forces (if needed).
    ! With special k-points (and point group specified)
    ! we symmetrize RHO.
    IF (cntl%tsymrho) CALL symvec(fion)
    ! ==--------------------------------------------------------------==
    IF (ropt_mod%convwf.AND..NOT.cntl%tdiagopt) GOTO 666
    ! ==--------------------------------------------------------------==
    IF (fint1%ttrot) THEN
       CALL epot(rhoe)
       IF (cntl%tdavi) THEN
          CALL ksdiag(vpp)
          !$omp parallel do private(IG)
          DO ig=1,ncpw%ngw
             vpp(ig)=EXP(fint1%betap*vpp(ig))
          ENDDO
       ENDIF
    ELSE
       IF (cntl%tfint) CALL dcopy(nstate*nkpt%nkpts,1._real_8,0,crge%f(1,1),1)
       IF (cntl%tdavi) CALL ksdiag(vpp)
    ENDIF
    ener_com%eeig=0._real_8
    ! ==--------------------------------------------------------------==
    ! == Loop over k points                                           ==
    ! ==--------------------------------------------------------------==
    IF ((cntl%diis.OR.cntl%pcg).AND..NOT.cntl%tdiagopt) THEN
       CALL k_diis_rhofix(c0,c2,sc0,tau0,fion,&
            cscr,cr,vpp,eigv,rhoe,psi,&
            nstate,tfor,tinfo,ddrho)
    ELSE
       CALL inq_swap(kbeg,kend,kinc)
       DO ikpt=kbeg,kend,kinc
          nkpoint=nkpbl(ikpt)
          IF (tkpts%tkblock) THEN
             IF (fint1%ttrot) THEN
                CALL rkpt_swap(c0,nstate,ikpt,&
                     'HGKP HGKM MASKGW TWNL EMK2 EIGKR ALM C0')
             ELSE
                CALL rkpt_swap(c0,nstate,ikpt,&
                     'HGKP HGKM MASKGW TWNL EIGKR C0')
             ENDIF
          ENDIF
          DO ikind = 1, nkpoint
             ikk=kpbeg(ikpt)+ikind
             IF (fint1%ttrot) CALL vbeta(rhoe,psi,ikind)
             IF (cntl%tlsd) THEN
                IF (cntl%tlanc) THEN
                   CALL dcopy(nnx,rhoe(1,1),1,rhot(1,1),1)
                   IF (tkpts%tkpnt) THEN
                      CALL friesner_c(spin_mod%nsup,REAL(spin_mod%nsup,kind=real_8),nconv,nhpsi,&
                           c0(:,:,ikind),c2,sc0,cscr,&
                           rhoe(:,1),psi(:,1),eigv(1,ikk),&
                           1,ikind,ikk,trefine,tinfo)
                   ELSE
                      IF (func1%mhfx.NE.0) CALL dcopy(2*ncpw%ngw*nstate,c0,1,cgs,1)
                      IF (cntl%tdmal) THEN
                         CALL dist_friesner(spin_mod%nsup,REAL(spin_mod%nsup,kind=real_8),&
                              nconv,nhpsi,spin_mod%nsup,&
                              c0,c2,sc0,cscr,cgs,crge%f,&
                              rhoe(:,1),psi(:,1),eigv,&
                              1,trefine,tinfo)
                      ELSE
                         CALL friesner(spin_mod%nsup,REAL(spin_mod%nsup,kind=real_8),&
                              nconv,nhpsi,spin_mod%nsup,&
                              c0,c2,sc0,cscr,cgs,crge%f(1:spin_mod%nsup,1),&
                              rhoe(:,1),psi(:,1),eigv,&
                              1,trefine,tinfo)
                      ENDIF
                      nac=nconv
                      CALL dcopy(nnx,rhot(1,1),1,rhoe,1)
                   ENDIF
                ELSEIF (cntl%tdavi)THEN
                   CALL dcopy(2*ncpw%ngw*nstate,c0,1,cgs,1)
                   CALL dcopy(nnx,rhoe,1,rhot(1,1),1)
                   CALL davidson(spin_mod%nsup,cnti%ndavv,c0(:,:,1),c2,cr,sc0,cscr,&
                        rhoe(:,1),cgs,crge%f(1:spin_mod%nsup,1),spin_mod%nsup,psi(:,1),edav,vpp)
                   CALL dcopy(2*ncpw%ngw*spin_mod%nsup,c0,1,cgs,1)
                   CALL dcopy(spin_mod%nsup,edav,1,eigv(1,ikk),1)
                   CALL dcopy(nnx,rhot(1,1),1,rhoe,1)
                ENDIF
                ib = spin_mod%nsup+1
                itaur=2
                IF (cntl%tlanc) THEN
                   IF (tkpts%tkpnt) THEN
                      CALL friesner_c(spin_mod%nsdown,REAL(spin_mod%nsdown,kind=real_8),nconv,nhpsi,&
                           c0(:,ib:,ikind),c2(1,ib),sc0,cscr,&
                           rhoe(:,2),psi(:,1),eigv(ib,ikk),&
                           2,ikind,ikk,trefine,tinfo)
                   ELSE
                      IF (cntl%tdmal) THEN
                         CALL dist_friesner(spin_mod%nsdown,&
                              REAL(spin_mod%nsdown,kind=real_8),nconv,nhpsi,spin_mod%nsdown,&
                              c0(:,ib:,1),c2(1,ib),sc0,cscr,&
                              cgs(1,ib),crge%f(ib,1),&
                              rhoe(:,2),psi(:,1),eigv(ib,1),&
                              2,trefine,tinfo)
                      ELSE
                         CALL friesner(spin_mod%nsdown,&
                              REAL(spin_mod%nsdown,kind=real_8),nconv,nhpsi,spin_mod%nsdown,&
                              c0(:,ib:,1),c2(1,ib),sc0,cscr,&
                              cgs(1,ib),crge%f(ib:ib+spin_mod%nsdown-1,1),&
                              rhoe(:,2),psi(:,1),eigv(ib,1),&
                              2,trefine,tinfo)
                      ENDIF
                      nbc=nconv
                      CALL dcopy(nnx,rhot(1,1),1,rhoe,1)
                   ENDIF
                ELSEIF (cntl%tdavi)THEN
                   CALL dcopy(2*ncpw%ngw*spin_mod%nsdown,cgs(1,ib),1,c0,1)
                   CALL davidson(spin_mod%nsdown,cnti%ndavv,c0(:,:,1),c2,cr,sc0,cscr,rhoe(:,2),&
                        cgs(:,ib:ib+spin_mod%nsdown-1),crge%f(ib:ib+spin_mod%nsdown-1,1),spin_mod%nsdown,&
                        psi(:,1),edav,vpp)
                   CALL dcopy(2*ncpw%ngw*spin_mod%nsdown,c0,1,cgs(1,ib),1)
                   CALL dcopy(spin_mod%nsdown,edav,1,eigv(ib,ikk),1)
                   CALL dcopy(2*ncpw%ngw*nstate,cgs,1,c0,1)
                   CALL dcopy(nnx,rhot(1,1),1,rhoe,1)
                ENDIF
                itaur=1
             ELSE
                IF (cntl%tlanc) THEN
                   IF (tkpts%tkpnt) THEN
                      IF (.NOT.cntl%tfrsblk) THEN
                         IF (tkpts%tonlydiag) THEN
                            CALL frie_c_p(nstate,crge%nel,nconv,nhpsi,&
                                 c0(:,:,ikind),c2,sc0,cscr,rhoe(:,1),psi(:,1),&
                                 eigv(1,ikk),&
                                 0,ikind,ikk,trefine,tinfo)
                         ELSE
                            CALL friesner_c(nstate,crge%nel,nconv,nhpsi,&
                                 c0(:,:,ikind),c2,sc0,cscr,rhoe(:,1),psi(:,1),&
                                 eigv(1,ikk),&
                                 0,ikind,ikk,trefine,tinfo)
                         ENDIF
                      ELSE
                         CALL frsblk_c(nstate,crge%nel,nconv,nhpsi,&
                              c0(:,:,ikind),c2,sc0,cscr,rhoe(:,1),psi(:,1),&
                              eigv(1,ikk),&
                              0,ikind,ikk,trefine,tinfo)
                      ENDIF
                   ELSE
                      IF (.NOT.cntl%tfrsblk) THEN
                         IF (func1%mhfx.NE.0) CALL dcopy(2*ncpw%ngw*nstate,c0,1,cgs,1)
                         IF (cntl%tdmal) THEN
                            CALL dist_friesner(nstate,crge%nel,nconv,nhpsi,nstate,&
                                 c0,c2,sc0,cscr,cgs,crge%f,rhoe(:,1),psi(:,1),eigv,&
                                 0,trefine,tinfo)
                         ELSE
                            CALL friesner(nstate,crge%nel,nconv,nhpsi,nstate,&
                                 c0,c2,sc0,cscr,cgs,crge%f(:,1),rhoe(:,1),psi(:,1),eigv,&
                                 0,trefine,tinfo)
                         ENDIF
                      ELSE
                         CALL frsblk(nstate,crge%nel,nconv,nhpsi,&
                              c0,c2,sc0,cscr,rhoe(:,1),psi(:,1),eigv,&
                              0,trefine,tinfo)
                      ENDIF
                      nac=nconv
                   ENDIF
                ELSE
                   ! CGS is only needed for HF
                   IF (func1%mhfx.NE.0) CALL dcopy(2*ncpw%ngw*nstate,c0,1,cgs,1)
                   CALL davidson(nstate,cnti%ndavv,c0(:,:,1),c2,cr,sc0,cscr,rhoe(:,1),&
                        cgs,crge%f(:,1),nstate,psi(:,1),eigv,vpp)
                ENDIF
             ENDIF
             ! Check memory
          ENDDO               ! LOOP OVER IKIND
          IF (tkpts%tkblock) THEN
             CALL wkpt_swap(c0,nstate,ikpt,'C0')
          ENDIF
       ENDDO                   ! LOOP OVER IKPT
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == End loop over kpoints                                        ==
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    ! == Solve K.P eigenvalue equations                               ==
    ! ==--------------------------------------------------------------==
    IF (tkdp.AND..NOT.cntl%tlsd) THEN
       ! 
       ! 1) Compute < c_i | p | c_j > and <c_i| H |c_j> matrix elements
       CALL kdp_prep(nkdp,c0,nstate,pkdp,xlambda0,ckdp,rhoe,psi,&
            fint1%betap,eigv)

       ! 2) Solve the eigenvalue equations: eigenvalues in ekdp and
       ! eigenvectors in akdp.
       CALL kdp_diag(nstate,pkdp,xlambda0,xkdp,nkdp,ekdp,akdp,&
            noccup,auxdiag,rauxdiag)

       ! 3) Adjust chemical potential and find occupation numbers.
       CALL adjmu(nstate,nkdp,crge%nel,fint1%betael,ekdp,wkdp,cntl%tlsd,ener_com%amu,tadjmu)
       CALL occup(ekdp,wkdp,ener_com%amu,fint1%betael,fkdp,nstate,nkdp,crge%nel,cntl%tlsd)
       ! 
       ! 4) Calculate the band energy
       CALL enert(ekdp,wkdp,fkdp,ener_com%amu,fint1%betael,& !vw swap fkdp,ener_com%amu
            nstate,nstate,nkdp,crge%nel,cntl%tlsd,&
            ener_com%eeig,ener_com%eband,ener_com%entropy)
       ! 
       ! 5) Calculate output density
       CALL kdp_rho(nkdp,c0,ckdp,rhoe,psi,nstate,&
            akdp,noccup,bmix,fkdp,wkdp,bkdp)
       ! 
       GOTO 1111
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    ! Calculate occupation numbers
    IF (fint1%ttrot) THEN
       DO ikind=1,nkpt%nkpts
          DO i=1,nstate
             eigv(i,ikind)=-LOG(eigv(i,ikind))/fint1%betap
          ENDDO
       ENDDO
    ENDIF
    tadjmu=.FALSE.
    IF (tfixo) THEN
       ener_com%amu=eigv(nstate,1)
    ELSEIF (cntl%tlsd.AND.tdsp1%tfxsp) THEN
       ! ..for fixed spin, get amu and f of each spin species
       ! ..prepare eigenvalues in WE array
       ! ..alpha spin
       ALLOCATE(we(spin_mod%nsup,nkpt%nkpts),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(focc(spin_mod%nsup,nkpt%nkpts),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       DO ik=1,nkpt%nkpts
          CALL dcopy(spin_mod%nsup,eigv(1,ik),1,we(1,ik),1)
       ENDDO
       CALL adjmu(spin_mod%nsup,nkpt%nkpts,REAL(tdsp1%nupel,kind=real_8),fint1%betael,we,wk,&
            cntl%tlsd,amuup,tadjmu)
       CALL occup(we,wk,amuup,fint1%betael,focc,spin_mod%nsup,nkpt%nkpts,&
            REAL(tdsp1%nupel,kind=real_8),cntl%tlsd)
       DO ik=1,nkpt%nkpts
          CALL dcopy(spin_mod%nsup,focc(1,ik),1,crge%f(1,ik),1)
       ENDDO
       DEALLOCATE(we,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(focc,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       ! ..beta spin
       ALLOCATE(we(spin_mod%nsdown,nkpt%nkpts),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(focc(spin_mod%nsdown,nkpt%nkpts),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       DO ik=1,nkpt%nkpts
          CALL dcopy(spin_mod%nsdown,eigv(1+spin_mod%nsup,ik),1,we(1,ik),1)
       ENDDO
       CALL adjmu(spin_mod%nsdown,nkpt%nkpts,REAL(tdsp1%ndoel,kind=real_8),fint1%betael,we,wk,&
            cntl%tlsd,amudo,tadjmu)
       CALL occup(we,wk,amudo,fint1%betael,focc,&
            spin_mod%nsdown,nkpt%nkpts,REAL(tdsp1%ndoel,kind=real_8),cntl%tlsd)
       DO ik=1,nkpt%nkpts
          CALL dcopy(spin_mod%nsdown,focc(1,ik),1,crge%f(spin_mod%nsup+1,ik),1)
       ENDDO
       DEALLOCATE(focc,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(we,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    ELSE
       CALL adjmu(nstate,nkpt%nkpts,crge%nel,fint1%betael,eigv,wk,cntl%tlsd,ener_com%amu,tadjmu)
       CALL occup(eigv,wk,ener_com%amu,fint1%betael,crge%f,nstate,nkpt%nkpts,crge%nel,cntl%tlsd)
    ENDIF

    ! Check if the number of states is sufficient.
    IF (paral%parent.AND.cntl%tfint.AND.fint1%betael.LT.betaelmax) THEN
       DO ikind=1,nkpt%nkpts
          IF (cntl%tlsd) THEN
             IF (crge%f(spin_mod%nsup,ikind).GT.tollpar) THEN
                IF (paral%io_parent)&
                     WRITE(6,'(1X,64("!"))')
                IF (paral%io_parent)&
                     WRITE(6,'(" !!",A,T64,"!!")')&
                     ' UPDRHO| THE NUMBER OF ALPHA STATES IS TOO SMALL'
                IF (paral%io_parent)&
                     WRITE(6,'(" !!",A,I4,A,I3,A,1PE15.5,T64,"!!")')&
                     '         F(',spin_mod%nsup,',',ikind,') =',crge%f(spin_mod%nsup,ikind)
                IF (paral%io_parent)&
                     WRITE(6,'(1X,64("!"))')
             ENDIF
             IF (crge%f(nstate,ikind).GT.tollpar) THEN
                IF (paral%io_parent)&
                     WRITE(6,'(1X,64("!"))')
                IF (paral%io_parent)&
                     WRITE(6,'(" !!",A,T64,"!!")')&
                     ' UPDRHO| THE NUMBER OF BETA STATES IS TOO SMALL'
                IF (paral%io_parent)&
                     WRITE(6,'(" !!",A,I4,A,I3,A,1PE15.5,T64,"!!")')&
                     '         F(',spin_mod%nsdown,',',ikind,') =',&
                     crge%f(nstate,ikind)
                IF (paral%io_parent)&
                     WRITE(6,'(1X,64("!"))')
             ENDIF
          ELSE
             IF (crge%f(nstate,ikind).GT.tollpar) THEN
                IF (paral%io_parent)&
                     WRITE(6,'(1X,64("!"))')
                IF (paral%io_parent)&
                     WRITE(6,'(" !!",A,T64,"!!")')&
                     ' UPDRHO| THE NUMBER OF STATES IS TOO SMALL'
                IF (paral%io_parent)&
                     WRITE(6,'(" !!",A,I4,A,I3,A,1PE15.5,T64,"!!")')&
                     '         F(',nstate,',',ikind,') =',&
                     crge%f(nstate,ikind)
                IF (paral%io_parent)&
                     WRITE(6,'(1X,64("!"))')
             ENDIF
          ENDIF
       ENDDO
    ENDIF
    ! Calculate free electron energy, band energy and kt*entropy.
    IF (cntl%tlsd.AND.tdsp1%tfxsp) THEN
       CALL enert(eigv,wk,crge%f,amuup,fint1%betael,&
            spin_mod%nsup,nstate,nkpt%nkpts,REAL(tdsp1%nupel,kind=real_8),.TRUE.,&
            eeig1,eband1,entropy1)
       CALL enert(eigv(1+spin_mod%nsup:,1),wk,crge%f(1+spin_mod%nsup:,1),amudo,fint1%betael,&
            spin_mod%nsdown,nstate,nkpt%nkpts,REAL(tdsp1%ndoel,kind=real_8),.TRUE.,&
            eeig2,eband2,entropy2)
       ener_com%eeig=eeig1+eeig2
       ener_com%eband=eband1+eband2
       ener_com%entropy=entropy1+entropy2
    ELSE
       ! Calculate free band energy.
       CALL enert(eigv,wk,crge%f,ener_com%amu,fint1%betael,nstate,nstate,nkpt%nkpts,crge%nel,cntl%tlsd,&
            ener_com%eeig,ener_com%eband,ener_com%entropy)
       ! IF Fixed occupation numbers EEIG is not correct.
       IF (tfixo) ener_com%eeig = ener_com%eband
    ENDIF
    ! Calculate output density
    IF (pslo_com%tivan) THEN
       CALL rnlsm(c0(:,:,1),nstate,1,1,.FALSE.)
    ENDIF
    IF (tkpts%tkpnt) THEN
       CALL rhoofr_c(c0,rhoe,psi(:,1),nstate)
    ELSE
       CALL rhoofr(c0(:,:,1),rhoe,psi(:,1),nstate)
    ENDIF
    ! 
    ! From k.p
1111 CONTINUE
    CALL dcopy(nnx,rhoe(1,1),1,rout0(1,1),1)
    ! 
    ! Total Energy summation
    ! Sum over EKIN (cal. in RHOOFR), EPSEU, ENL (unused)
    ! EHT, EHEI, EHEE, EHEP, EXC, VXC, EGC )
    CALL mp_sum(ener_com%ekin,parai%allgrp)
    CALL mp_sum(ener_com%epseu,parai%allgrp)
    CALL mp_sum(ener_com%enl,parai%allgrp)
    CALL mp_sum(ener_com%eht,parai%allgrp)
    CALL mp_sum(ener_com%ehep,parai%allgrp)
    CALL mp_sum(ener_com%ehee,parai%allgrp)
    CALL mp_sum(ener_com%ehii,parai%allgrp)
    CALL mp_sum(ener_com%exc,parai%allgrp)
    CALL mp_sum(ener_com%vxc,parai%allgrp)
    CALL mp_sum(ener_com%egc,parai%allgrp)

    ! Hartree-Fock contribution to energy
    CALL hfx(c0(:,:,1),c2,crge%f(:,1),psi(:,1),nstate,ehfx,vhfx,.TRUE.)
    ener_com%exc=ener_com%exc-ehfx
    ! If you want to have non-local pp contribution
    ! CALL RNLRH(ENL,NSTATE,NKPTS)
    ener_com%etot=ener_com%eeig+(ener_com%exc-ener_com%vxc)+ener_com%eht&
         +vdwr%evdw ! Empirical van der Waals correction
    ! Check total energy difference 
    detot=ener_com%etot-etot0
    IF (ABS(detot).LT.cntr%toldetot) THEN
       netot=netot+1
    ELSE
       netot=0
    ENDIF
    etot0=ener_com%etot
    ! Check for convergence
    CALL dcopy(nnx,rout0(1,1),1,rmix(1,1),1)
    CALL daxpy(nnx,-1.0_real_8,rin0(1,1),1,rmix(1,1),1)
    imax=idamax(nnx,rmix(1,1),1)
    gemax=ABS(dgive(rmix(1,1),imax))
    cnorm=dasum(nnx,rmix(1,1),1)
    CALL mp_sum(cnorm,parai%allgrp)
    CALL mp_max(gemax,parai%allgrp)
    nnxs=fpar%kr1s*fpar%kr2s*fpar%kr3s
    gemax=SQRT(parm%omega*gemax/nnxs)
    cnorm=parm%omega*cnorm/nnxs
    ddrho=gemax
    ! Check if soft exit
    IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)
    IF (gemax.LT.cntr%tolog .OR. soft_com%exsoft .OR. cntl%ksener .OR. netot.GE.2) THEN
       ropt_mod%convwf=.TRUE.
       ! Check if B2MAX is not too small (T.D.)
       IF (cntl%tlanc) THEN
          IF (nfr.EQ.1.AND.fint4%ntabtrot.GT.1.AND.cntr%b2limit.GT.1.e-8_real_8) THEN
             IF (.NOT.soft_com%exsoft) THEN
                IF (paral%parent) THEN
                   IF (paral%io_parent)&
                        WRITE(6,'(1X,64("!"))')
                   IF (paral%io_parent)&
                        WRITE(6,'(" !!",A,T64,"!!")')&
                        ' UPDRHO| B2MAX IS VERY HIGH'
                   IF (paral%io_parent)&
                        WRITE(6,'(" !!",A,T64,"!!")')&
                        ' UPDRHO| DOUBT ON CONVERGENCE. WE TRY ONCE AGAIN'
                   IF (paral%io_parent)&
                        WRITE(6,'(1X,64("!"))')
                ENDIF
                ropt_mod%convwf=.FALSE.
             ENDIF
          ENDIF
       ENDIF
       ! Calculate Bogoliubov correction
       IF (fint1%tbogo.AND.fint1%ttrot) THEN
          CALL dcopy(nnx,rin0(1,1),1,rhoe(1,1),1)
          CALL eback(0)
          CALL vofrho(tau0,fion,rhoe,psi,tfor,.FALSE.)
          CALL eback(1)
          CALL bogol(c0,c2,sc0,eigv,crge%f,rhoe,psi(:,1),&
               nstate,ener_com%ebogo)
          ener_com%etot=ener_com%etot+ener_com%ebogo
       ENDIF
    ENDIF
    IF (tkpts%tonlydiag) THEN
       ropt_mod%convwf=.TRUE.
    ENDIF
    ! Change b2limit if needed
    IF (cntl%tlanc) THEN
       CALL change_b2l(gemax,trefine)
    ENDIF
    ! Change Betap if needed
    IF (fint1%ttrot) THEN
       CALL change_betap(gemax,trefine)
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (.NOT.cntl%tlanc) THEN
       DEALLOCATE(edav,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    IF (.NOT.cntl%tlanc.AND.cntl%tlsd) THEN
       DEALLOCATE(cgs,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ELSE
       DEALLOCATE(cgs,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    ENDIF
    IF (cntl%tlsd) THEN
       DEALLOCATE(rhot,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == Density mixing                                               ==
    ! ==--------------------------------------------------------------==
    IF (broy1%tgmix) THEN
       CALL mixing_g(nfr,gemax,rin0,rout0,rmix,psi,thl)
    ELSE
       CALL mixing_r(nfr,gemax,rin0,rout0,rmix,thl)
    ENDIF
    CALL dcopy(nnx,rmix(1,1),1,rin0(1,1),1)
    ! ==--------------------------------------------------------------==
666 CONTINUE
    CALL tihalt(procedureN,isub)
    RETURN
  END SUBROUTINE updrho
  ! ==================================================================
  SUBROUTINE give_scr_updrho(lupdrho,tag,nstate,tfor,tstress)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lupdrho
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate
    LOGICAL                                  :: tfor, tstress

    INTEGER                                  :: lbogol, lhpsi, lmixing, &
                                                lptheory, lrhofix, lrhoofr, &
                                                lrnlsm, lscrdiag, lstress, &
                                                lsymvec, lvofrho

    lscrdiag=0
    CALL give_scr_vofrho(lvofrho,tag)
    IF (tfor) THEN
       CALL give_scr_rnlsm(lrnlsm,tag,nstate,tfor)
       CALL give_scr_symvec(lsymvec,tag)
    ELSE
       lrnlsm=0
       lsymvec=0
    ENDIF
    IF (tstress) THEN
       CALL give_scr_stress(lstress,tag)
    ELSE
       lstress=0
    ENDIF
    lrhofix=0
    IF (cntl%tlanc) THEN
       IF (.NOT.cntl%tfrsblk) THEN
          IF (cntl%tdmal) THEN
             CALL give_scr_dist_friesner(lscrdiag,tag,nstate)
          ELSE
             CALL give_scr_friesner(lscrdiag,tag,nstate)
          ENDIF
       ELSE
          IF (tkpts%tkpnt) THEN
             CALL give_scr_frsblk_c(lscrdiag,tag,nstate)
          ELSE
             CALL give_scr_frsblk(lscrdiag,tag,nstate)
          ENDIF
       ENDIF
    ELSEIF (cntl%tdavi) THEN
       CALL give_scr_davidson(lscrdiag,tag,nstate,cnti%ndavv)
    ELSEIF (cntl%diis.OR.cntl%pcg) THEN
       CALL give_scr_kforces(lrhofix,tag,nstate,.TRUE.,tfor)
    ENDIF
    IF (pslo_com%tivan) CALL give_scr_rnlsm(lrnlsm,tag,nstate,.FALSE.)
    CALL give_scr_rhoofr(lrhoofr,tag)
    lhpsi=0
    IF (fint1%tbogo) THEN
       CALL give_scr_bogol(lbogol,tag,nstate)
    ELSE
       lbogol=0
       IF (tfor)  CALL give_scr_hpsi(lhpsi,tag,nstate)
    ENDIF
    CALL give_scr_ptheory(lptheory,tag,nstate)
    CALL give_scr_mixing(lmixing,tag)
    lupdrho=MAX(lvofrho,lrnlsm,lsymvec,lstress,lrhofix,&
         lscrdiag,lrhoofr,lbogol,lhpsi,lptheory,lmixing)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_updrho
  ! ==================================================================

END MODULE updrho_utils
