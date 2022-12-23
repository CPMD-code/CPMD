MODULE cplngs_utils
  USE andp,                            ONLY: rin0
  USE canon_utils,                     ONLY: canon
  USE coor,                            ONLY: fion,&
                                             tau0
  USE cplngsmod,                       ONLY: &
       csurf, eps_c, f_high, f_low, f_med, fc_high, fc_low, fc_med, iatfd, &
       isurf, natfd, nsurf, nvect, tallat, talldof, tcpl, tcplfd, tcpllr, &
       tolcpl, tspecv
  USE cppt,                            ONLY: gk,&
                                             inyh,&
                                             nzh,&
                                             rhops,&
                                             scg,&
                                             vps
  USE eicalc_utils,                    ONLY: eicalc,&
                                             eicalc1
  USE eind_ii_utils,                   ONLY: eind_ii
  USE eind_loc_utils,                  ONLY: eind_loc
  USE eind_nl_utils,                   ONLY: eind_nl
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE fftmain_utils,                   ONLY: fwfftn
  USE fnlalloc_utils,                  ONLY: fnl_set,&
                                             fnlalloc
  USE forcedr_driver,                  ONLY: forcedr
  USE forcep_utils,                    ONLY: rhoe_psi_size
  USE forces_diag_utils,               ONLY: forces_diag
  USE geq0mod,                         ONLY: geq0
  USE hpsi_utils,                      ONLY: hpsi
  USE inscan_utils,                    ONLY: inscan
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE ksdiag_utils,                    ONLY: ksdiag
  USE linres,                          ONLY: lr01,&
                                             lr02,&
                                             lr03
  USE lr_xcpot_utils,                  ONLY: lr_xcpot
  USE lscal,                           ONLY: hesscr,&
                                             lvlhes,&
                                             nvar
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum
  USE nl_res_utils,                    ONLY: give_scr_nl_res,&
                                             nl_res
  USE nlps,                            ONLY: ndfnl,&
                                             nghtol,&
                                             nlm,&
                                             nlps_com,&
                                             wsg
  USE opt_lr_utils,                    ONLY: give_scr_opt_lr,&
                                             opt_lr
  USE ovlap_utils,                     ONLY: ovlap
  USE parac,                           ONLY: parai,&
                                             paral
  USE phfac_utils,                     ONLY: phfac
  USE poin,                            ONLY: potr,&
                                             rhoo
  USE pslo,                            ONLY: pslo_com
  USE rho1ofr_utils,                   ONLY: rho1ofr,&
                                             rhoabofr
  USE rhoofr_utils,                    ONLY: give_scr_rhoofr,&
                                             rhoofr
  USE rnlsm1_utils,                    ONLY: rnlsm1
  USE rnlsm2_utils,                    ONLY: rnlsm2
  USE rnlsm_2d_utils,                  ONLY: give_scr_rnlsm_2d,&
                                             rnlsm_2d
  USE rnlsm_utils,                     ONLY: give_scr_rnlsm
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
                                             spin_mod
  USE summat_utils,                    ONLY: give_scr_summat,&
                                             summat
  USE symtrz_utils,                    ONLY: give_scr_symvec,&
                                             symvec
  USE system,                          ONLY: &
       cnti, cntl, fpar, ipept, maxsys, ncpw, nkpt, parap, parm, spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE v1ofrho1_utils,                  ONLY: give_scr_v1ofrho1,&
                                             v1ofrho1
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: cplngs

  PUBLIC :: cplana
  PUBLIC :: get_eind
  PUBLIC :: cplconfig
  PUBLIC :: prcplngs
  !public :: cpl_wreigen
  PUBLIC :: cplngs_init
  PUBLIC :: cpl_para
  !public :: rnlcpl
  PUBLIC :: give_scr_cplngs
  !public :: give_scr_cplsub
  PUBLIC :: give_scr_get_eind
  !public :: switch_st
  !public :: tdrhof_dis
  !public :: cplrdvec

CONTAINS

  ! ==================================================================
  SUBROUTINE cplngs (c0, c2, cf, c0ini, cr, sc0, cscr, vpp, eigv,&
       RHOE, RHOINI, PSI, &
       TAU0, FION, IFCALC, IREC, TINFO, TFOR,&
       CPLVEC, CPLCONF, NSTATE)
    ! ==--------------------------------------------------------------==
    ! ==      DRIVER FOR CALCULATION OF INTERSURFACE COUPLINGS        ==
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8)                          :: c0(*), c2(*), cf(*), &
                                                C0INI(*), cr(*), sc0(*), &
                                                cscr(*)
    REAL(real_8)                             :: vpp(:), eigv(:), RHOE(:,:), &
                                                RHOINI(*)
    COMPLEX(real_8)                          :: psi(:,:)
    REAL(real_8)                             :: tau0(:,:,:), fion(:,:,:)
    INTEGER                                  :: ifcalc, irec(:)
    LOGICAL                                  :: tinfo, tfor
    REAL(real_8)                             :: CPLVEC(:,:,:,:), &
                                                CPLCONF(:,:,:)
    INTEGER                                  :: nstate

    INTEGER                                  :: isub

    IF (.NOT.tcpl) RETURN
    CALL tiset('    CPLNGS',isub)
    IF (tcplfd) THEN
       CALL cplfd (c0, c2, cf, c0ini, cr, sc0, cscr, vpp, eigv,&
            RHOE, RHOINI, PSI, &
            TAU0, FION, IFCALC, IREC, TINFO, TFOR,&
            CPLVEC, NSTATE)
    ELSE
       CALL cplana (c0, cf, c2, cr, sc0, cscr, vpp, rhoini, cplvec,&
            CPLCONF, NSTATE, RHOE, PSI, EIGV)
    ENDIF
    CALL cplconfig (cplvec, cplconf)
    CALL tihalt('    CPLNGS',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE cplngs
  ! ==================================================================
  SUBROUTINE cplana (c0, c1, c2, cr, sc0, cscr, vpp, rhoini, cplvec,&
       CPLCONF, NSTATE, RHOE, PSI, EIGV)
    ! ==--------------------------------------------------------------==
    ! ==        ANALYTIC CALCULATION OF INTERSURFACE COUPLINGS        ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: CR(*), cscr(*)
    REAL(real_8)                             :: vpp(:)
    REAL(real_8), TARGET                     :: rhoini(fpar%nnr1,clsd%nlsd)
    REAL(real_8)                             :: cplvec(:,:,:,:), &
                                                cplconf(:,:,:)
    INTEGER                                  :: nstate
    COMPLEX(real_8) :: SC0(nkpt%ngwk,nstate), c2(nkpt%ngwk,nstate), &
      c1(nkpt%ngwk,nstate), c0(nkpt%ngwk,nstate)
    REAL(real_8), TARGET                     :: rhoe(fpar%nnr1,clsd%nlsd)
    COMPLEX(real_8)                          :: PSI(:,:)
    REAL(real_8)                             :: eigv(:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'cplana'

    CHARACTER(len=10)                        :: TYPE
    CHARACTER(len=5)                         :: orbtyp
    COMPLEX(real_8), ALLOCATABLE             :: eirop(:), eirop1(:), &
                                                eivps(:), eivps1(:), h1nl(:)
    EXTERNAL                                 :: dasum, ddot, dotp
    INTEGER :: i, ia, iact, iat, ic, ierr, ig, il_rhoe_1d, il_rhoe_2d, imax, &
      ir, is1, is2, isp, ISUB, IV, IVV, K, KMAX, LDDXC_1d, LDDXC_2d, LDFNL, &
      LDONE, LEIROP, LEIVPS, LH1NL, LVECT, N3, NATOMS
    INTEGER, ALLOCATABLE                     :: idone(:)
    LOGICAL                                  :: debug, lconv, lprint, &
                                                lprint2, tfdens, tksmat
    REAL(real_8) :: c, dasum, ddot, deinv, dotp, e2(5), eexc, egnd, eind, fc, &
      fcomp, fcomp2, fcsum, rlen, rmax, tol_lr_bak, trace, v(3), vc(3), vp
    REAL(real_8), ALLOCATABLE                :: ddxc(:,:), deinvv(:), &
                                                fcsrf(:), vect(:,:,:)
    REAL(real_8), POINTER                    :: drhoe(:,:), rhoov(:,:)

! DRHOE(NNR1,NLSD)
! RHOOV(:,1)
! ==--------------------------------------------------------------==

    lprint  = .TRUE.          ! Progress report
    lprint2 = .FALSE.         ! Diagnostics
    lprint  = (lprint  .AND. paral%parent)
    lprint2 = (lprint2 .AND. paral%parent)
    CALL tiset ('    CPLANA', isub)

    ! For linear response - space for locally stored arrays
    IF (tcpllr) THEN
       ! scratch will be allocated inside CPLANA

       leivps = 2*ncpw%nhg
       leirop = 2*ncpw%nhg
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
       CALL rhoe_psi_size (il_rhoe_1d=il_rhoe_1d,il_rhoe_2d=il_rhoe_2d)
       IF (talldof) THEN
          lvect = 3*ions1%nat
       ELSE
          lvect = 3*ions1%nat*nvect
       ENDIF
       ldone = 3*ions1%nat/2+1
       ! 
       debug = .FALSE.

       ! TODO check stat
       ALLOCATE(eivps(leivps),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(eirop(leirop),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(eivps1(leivps),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(eirop1(leirop),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(h1nl(lh1nl),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(rhoov(il_rhoe_1d,il_rhoe_2d),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(ddxc(lddxc_1d,lddxc_2d),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(deinvv(nsurf),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)

       IF (talldof) THEN
          ALLOCATE(vect(3, ions1%nat, 1),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
       ELSE
          ALLOCATE(vect(3, ions1%nat, nvect),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
       ENDIF

       ! TODO understand why memory90 instead of allocate?
       ALLOCATE(ddfnl(ions1%nat,maxsys%nhxs,6,ndfnl),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(potr(il_rhoe_1d,il_rhoe_2d),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)

       IF (.NOT.talldof) THEN
          ! TODO check stat
          ALLOCATE(idone(ldone),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
          ALLOCATE(fcsrf(nsurf),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
       ENDIF

       rhoo => rhoini
       drhoe => rhoe
    ELSE                      ! IF (TCPLLR) THEN
       ! scratch was allocated outside CPLANA,
       ! do not need all these arrays, will exit early

       rhoov => rhoe
    ENDIF
    IF (tspecv) CALL cplrdvec (vect, ions1%nat, nvect, tfdens)
    ! ==--------------------------------------------------------------==
    ! ==  Terms at constant density                                   ==
    ! ==--------------------------------------------------------------==
    CALL zeroing(cplvec)!, 3*maxsys%nax*maxsys%nsx*nsurf)
    tksmat = lspin2%tlse
    IF (tksmat) THEN
       ! Excitation between lowest two singlet states
       IF (lspin2%tlse) THEN
          CALL switch_st ('SET', crge%f, nstate, lspin2%tlse, lspin2%tlsets)
          CALL switch_st ('GND', crge%f, nstate, lspin2%tlse, lspin2%tlsets)
          CALL forcedr (c0, c2, sc0, rhoe, psi, tau0, fion, eigv,&
               NSTATE, 1, .FALSE., .FALSE.)
          egnd = ener_com%etot
          IF ((lprint).AND.paral%io_parent)&
               WRITE(6,110) 'ENERGY OF STATE S0:', egnd
          CALL switch_st ('EXC', crge%f, nstate, lspin2%tlse, lspin2%tlsets)
          CALL forcedr (c0, c2, sc0, rhoe, psi, tau0, fion, eigv,&
               NSTATE, 1, .FALSE., .FALSE.)
          eexc = ener_com%etot
          IF ((lprint).AND.paral%io_parent)&
               WRITE(6,110) 'ENERGY OF STATE S1:', eexc
          CALL switch_st ('DEF', crge%f, nstate, lspin2%tlse, lspin2%tlsets)
          IF ((lprint).AND.paral%io_parent)&
               WRITE(6,110)&
               'EXCITATION ENERGY FROM S0 TO S1:', EEXC-EGND
       ENDIF
       ! Diagonal KS matrix elements, swapped states
       DO i = 1,nstate-2
          CALL dcopy (ncpw%ngw, c0(1,i), 1, c1(1,i), 1)
       ENDDO
       CALL dcopy (ncpw%ngw, c0(1,nstate-1), 1, c1(1,nstate), 1)
       CALL dcopy (ncpw%ngw, c0(1,nstate), 1, c1(1,nstate-1), 1)
       CALL forcedr (c1, c2, sc0, rhoe, psi, tau0, fion, eigv,&
            NSTATE, 1, .FALSE., .FALSE.)
       IF (lprint.AND.paral%io_parent) WRITE (6,110)&
            'ENERGY OF STATE TS (SWAPPED ORBITALS):', ener_com%etot
       CALL hpsi (c1, c2, sc0, rhoe, psi(:,1), nstate,&
            1, clsd%nlsd)
       CALL ovlap (nstate, eigv, c1, c2)
       CALL summat (eigv, nstate)
       IF (.FALSE..AND.paral%parent) THEN
          IF (paral%io_parent)&
               WRITE(6,*)'KS MATRIX (SWAPPED ORBITALS)'
          DO i = 1,nstate
             PRINT '(10F12.8)', (eigv(k+(i-1)*nstate), k=1,nstate)
          ENDDO
       ENDIF
       DO ic = 1,nsurf
          is1 = isurf(1,ic)
          is2 = isurf(2,ic)
          eexc = eigv(is2+(is2-1)*nstate)
          egnd = eigv(is1+(is1-1)*nstate)
          deinvv(ic) = eexc-egnd
          IF ((lprint).AND.paral%io_parent)&
               WRITE(6,110)&
               'EXCITATION ENERGY (SWAPPED ORBITALS):', EEXC-EGND
       ENDDO
       ! Diagonal KS matrix elements, right order
       CALL forcedr (c0, c2, sc0, rhoe, psi, tau0, fion, eigv,&
            NSTATE, 1, .FALSE., .FALSE.)
       IF (lprint.AND.paral%io_parent)&
            WRITE (6,110) 'ENERGY OF STATE TS:', ener_com%etot

       CALL hpsi (c0, c2, sc0, rhoe, psi(:,1), nstate,&
            1, clsd%nlsd)
       CALL ovlap (nstate, eigv, c0, c2)
       CALL summat (eigv, nstate)
       DO ic = 1,nsurf
          is1 = isurf(1,ic)
          is2 = isurf(2,ic)
          eexc = eigv(is1+(is1-1)*nstate)
          egnd = eigv(is2+(is2-1)*nstate)
          deinvv(ic) = 0.5_real_8 * (deinvv(ic) + eexc-egnd)
          IF ((lprint).AND.paral%io_parent)&
               WRITE(6,110)&
               'EXCITATION ENERGY FROM KS ORBITALS:', EEXC-EGND
       ENDDO
       IF (.FALSE..AND.paral%parent) THEN
          IF (paral%io_parent)&
               WRITE(6,*)'KS MATRIX'
          DO i = 1,nstate
             PRINT '(10F12.8)', (eigv(k+(i-1)*nstate), k=1,nstate)
          ENDDO
       ENDIF
       orbtyp="GENER"
    ELSE
       ! Other excitations
       CALL forcedr (c0, c2, sc0, rhoe, psi, tau0, fion, eigv,&
            NSTATE, 1, .FALSE., .FALSE.)
       CALL canon (c0, c2, crge%f, crge%n, eigv)
       orbtyp="CANON"
    ENDIF
    ! Overlap density and local forces
    DO ic = 1,nsurf
       is1 = isurf(1,ic)
       is2 = isurf(2,ic)
       CALL zeroing(rhoov(:,1))
       CALL rhoabofr (1, c0(:,is1:is1), c0(:,is2:is2), rhoov(:,1), psi(:,1))
       ! Transform density to G space
       !$omp parallel do private(IR)
       DO ir = 1,fpar%nnr1
          psi(ir,1) = CMPLX(rhoov(ir,1),0._real_8,kind=real_8)
       ENDDO
       CALL  fwfftn (psi(:,1),.FALSE.,parai%allgrp)
       CALL tdrhof_dis (cplvec(:,:,:,ic), rhoov, psi)
    ENDDO
    ! Forces due to nonlocal PP
    CALL rnlsm1 (c0, nstate, 1)
    CALL rnlsm2 (c0, nstate, 1, 1)
    CALL rnlcpl (cplvec)
    CALL mp_sum(cplvec,3*maxsys%nax*maxsys%nsx*nsurf,parai%allgrp)
    DO ic = 1,nsurf
       is1 = isurf(1,ic)
       is2 = isurf(2,ic)
       IF (tksmat) THEN
          deinv = 1.0_real_8 / deinvv(ic)
          IF ((lprint).AND.paral%io_parent)&
               WRITE(6,110)&
               'EIGENVALUE DIFFERENCE: ', DEINVV(IC)
       ELSE
          deinv = 1.0_real_8 / (eigv(is2)-eigv(is1))
          IF ((lprint).AND.paral%io_parent)&
               WRITE(6,120) 'EXCITATION ENERGY (STATES ',&
               IS1, ' - ', IS2, '):', EIGV(IS2)-EIGV(IS1)
       ENDIF
       IF (tcpllr) deinvv(ic) = deinv
       CALL dscal (3*maxsys%nax*maxsys%nsx, -deinv, cplvec(1,1,1,ic), 1)
       IF (cntl%tsymrho.AND.paral%parent)&
            CALL SYMVEC (CPLVEC(:,:,:,IC))
    ENDDO
    ! ==--------------------------------------------------------------==
    ! ==  Linear response for dependence of KS operator on density    ==
    ! ==--------------------------------------------------------------==
    IF (.NOT.tcpllr) THEN
       CALL tihalt ('    CPLANA', isub)
       RETURN
    ENDIF
    IF (lprint) THEN
       IF (paral%io_parent)&
            WRITE(6,'(A)') ' COUPLINGS AT CONSTANT DENSITY'
       CALL prcplngs (cplvec, cplconf, .FALSE.)
    ENDIF
    ! Recalculate density, but re-use POTR, FNL, DFNL, EIVPS, and EIROP
    CALL dcopy (fpar%nnr1*clsd%nlsd, rhoe, 1, potr, 1)
    CALL rhoofr(c0(:,1:nstate), rhoe, psi(:,1), nstate)
    CALL dcopy (fpar%nnr1*clsd%nlsd, rhoe, 1, rhoo, 1)
    ! Non-local pseudopotential, calculate DDFNL - allocated above
    CALL rnlsm_2d (c0, nstate)
    ! Save FNL
    CALL fnl_set ("SAVE")
    CALL fnlalloc (nstate, .TRUE., .FALSE.)
    ! Local pseudopotential, use structure factors calculated before
    CALL eicalc (eivps, eirop)
    ! Calculate and store analytic dmu/dn
    CALL lr_xcpot (ddxc, rhoo, .FALSE.)
    ! Preconditioning for LR optimization - cf. SDLINRES
    cntl%prec = .TRUE.
    CALL ksdiag (vpp)
    IF (tksmat) THEN
       trace = dasum (nstate, eigv, nstate+1) / REAL(nstate,kind=real_8)
    ELSE
       trace = dasum (nstate, eigv, 1) / REAL(nstate,kind=real_8)
    ENDIF
    !$omp parallel do private(IG,VP) shared(TRACE)
    DO ig = 1,ncpw%ngw
       vp = vpp(ig) + trace
       vpp(ig) = ABS(vp/(vp**2+lr02%lr_hthrs**2))
    ENDDO
    ! ==--------------------------------------------------------------==
    ! ==  Unit vector along constant density coupling vector          ==
    ! ==--------------------------------------------------------------==
    CALL cplconfig (cplvec, cplconf)
    rlen = 0.0_real_8
    iat = 0
    DO isp = 1,ions1%nsp
       DO ia = 1,ions0%na(isp)
          iat = iat + 1
          IF (tfdens.OR..NOT.tspecv) THEN
             DO k = 1,3
                vect(k,iat,1) = cplconf(k,ia,isp)
                rlen = rlen + vect(k,iat,1)**2
             ENDDO
          ENDIF
       ENDDO                ! IA = 1,NA(ISP)
    ENDDO                    ! ISP = 1,NSP
    natoms = iat
    n3     = 3*natoms
    IF (tfdens.OR..NOT.tspecv) THEN
       rlen = 1.0_real_8 / SQRT(rlen)
       DO iat = 1,natoms
          vect(1,iat,1) = rlen*vect(1,iat,1)
          vect(2,iat,1) = rlen*vect(2,iat,1)
          vect(3,iat,1) = rlen*vect(3,iat,1)
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ==  Brute-force loop over all DOF of ions                       ==
    ! ==--------------------------------------------------------------==
    IF (talldof) THEN
       ! 
       ! ==  Main finite difference loop over the dimensions             ==
       ! 
       iat = 0
       iact = 1
       DO isp = 1,ions1%nsp
          DO ia = 1,ions0%na(isp)
             iat = iat + 1
             IF (tallat.OR.iat.EQ.iatfd(iact)) THEN
                iact = iact + 1
             ELSE
                DO ic = 1,nsurf
                   cplvec(k,ia,isp,ic) = 0.0_real_8
                ENDDO
                GOTO 9
             ENDIF
             DO k = 1,3
                IF ((lprint).AND.paral%io_parent)&
                     WRITE(6,'(A,3I3)')&
                     ' TYPE, ATOM, DIMENSION: ', ISP, IA, K
                ! 
                ! ==  Linear response of nuclear displacement                     ==
                ! 
                ! The following lines stem from SDLINRES
                eind = 0.0_real_8
                CALL eind_ii (tau0, isp, ia, k, eind, iteropt%iesr)
                CALL eind_loc (eind, isp, iat, k, rhoo, psi(:,1), eirop)
                ! Calculate constant force from NL-PP
                CALL fnl_set ("SWITCH")
                CALL nl_res (isp, iat, k, h1nl, nstate)
                CALL eind_nl (eind, isp, iat, k, nstate)
                CALL fnl_set ("SWITCH")
                ! Calculate the constant part of the local potential
                CALL eicalc1 (k, isp, iat, eivps1, eirop1)
                ! Calculate first-order wavefunction
                ropt_mod%sdiis = .TRUE.
                ropt_mod%spcg  = .TRUE.
                CALL opt_lr (c0, c1, c2, sc0, eind, eigv, drhoe,&
                     H1NL, EIVPS1, EIROP1, DDXC, VPP, CSCR, CR,&
                     PSI, NSTATE, "PHONON", ORBTYP)
                ! End SDLINRES - calculate first order density
                CALL rho1ofr (c0, c1, crge%f(:,1), drhoe, psi(:,1), nstate)
                ! Hartree and XC first-order terms
                CALL v1ofrho1 (e2, drhoe, ddxc, psi)
                ! Coupling component for each pair of surfaces
                DO ic = 1,nsurf
                   is1 = isurf(1,ic)
                   is2 = isurf(2,ic)
                   IF (ic.GT.1) THEN
                      CALL zeroing(rhoov(:,1))!, nnr1)
                      CALL rhoabofr (1, c0(:,is1:is1), c0(:,is2:is2), rhoov(:,1), psi(:,1))
                      deinv = deinvv(ic)
                   ENDIF
                   IF (lspin2%tlse) THEN
                      fcomp = 0.5_real_8 * (&
                           DDOT (fpar%nnr1, RHOOV, 1, DRHOE(1,2), 1) +&
                           DDOT (fpar%nnr1, RHOOV, 1, DRHOE(1,3), 1))
                   ELSE IF (cntl%tlsd.AND.is2.GT.spin_mod%nsup) THEN
                      fcomp = ddot (fpar%nnr1, rhoov, 1, drhoe(1,2), 1)
                   ELSE
                      fcomp = ddot (fpar%nnr1, rhoov, 1, drhoe(1,1), 1)
                   ENDIF
                   fcomp = fcomp * parm%omega/(spar%nr1s*spar%nr2s*spar%nr3s) * deinv
                   CALL mp_sum(fcomp,parai%allgrp)
                   IF (lprint) THEN
                      IF (paral%io_parent)&
                           WRITE(6,'(A,2F12.5)') ' FIXED DENSITY AND RESPONSE: ',&
                           CPLVEC(K,IA,ISP,IC), FCOMP
                      IF (paral%io_parent)&
                           WRITE(6,'(A,F12.5)') ' TOTAL COUPLING:             ',&
                           CPLVEC(K,IA,ISP,IC) + FCOMP
                   ENDIF
                   cplvec(k,ia,isp,ic) = cplvec(k,ia,isp,ic) + fcomp
                   ! 
                   fcomp  = dotp (ncpw%ngw, c0(1,is1), c1(1,is2))
                   fcomp2 = dotp (ncpw%ngw, c0(1,is2), c1(1,is1))
                   CALL mp_sum(fcomp,parai%allgrp)
                   CALL mp_sum(fcomp2,parai%allgrp)
                   IF ((lprint).AND.paral%io_parent)&
                        WRITE(6,'(A,2F12.5)')&
                        'COUPLINGS FROM OVERLAP: ', FCOMP, FCOMP
                ENDDO    ! IC = 1,NSURF
                IF (.FALSE.) THEN
                   DO i = 1,nstate
                      CALL zeroing(rhoov(:,1))!, nnr1)
                      CALL rhoabofr (1, c0(:,i:i), c0(:,i:i), rhoov(:,1), psi(:,1))
                      fcomp = ddot (fpar%nnr1, rhoov, 1, drhoe(1,1), 1)
                      CALL mp_sum(fcomp,parai%allgrp)
                      fcomp = fcomp * parm%omega/(spar%nr1s*spar%nr2s*spar%nr3s)
                      IF (paral%io_parent)&
                           WRITE(6,*)'State ', i, ': ', fcomp
                   ENDDO
                ENDIF
             ENDDO        ! K = 1,3
9            CONTINUE      ! jump here if atom is not interesting
          ENDDO            ! IA = 1,NA(ISP)
       ENDDO                ! ISP = 1,NSP
       ! ==--------------------------------------------------------------==
       ! ==  Iterative scheme for contributions of response of density   ==
       ! ==--------------------------------------------------------------==
    ELSE                      ! IF (TALLDOF)
       ! 
       ! ==  Response of density along the direction at constant density ==
       ! 
       IF ((lprint).AND.paral%io_parent)&
            WRITE(6,'(A,A)') ' LINEAR RESPONSE ALONG ',&
            'CONSTANT DENSITY COUPLING VECTOR'
       iv = 1
       ! Clear history of searched components
       !$omp parallel do private(K)
       DO k = 1,3*ions1%nat
          idone(k) = 0
       ENDDO
       ! 
       ! ==  Iterative search of the direction                           ==
       ! 
10     CONTINUE           ! DO ... UNTIL (IV.GE.NVECT .OR. LCONV)
       IF (lprint) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,I3,A)') ' VECTOR ', iv, ':'
          IF (paral%io_parent)&
               WRITE(6,'(2(2X,3F8.4))')&
               ((VECT(K,IAT,IV), K=1,3), IAT=1,NATOMS)
       ENDIF
       ! Calculate constant part of local and non-local interactions
       CALL get_eind(vect(1,1,iv), eind, eivps1, eirop1, h1nl,&
            TAU0, iteropt%iesr, RHOO, PSI, NSTATE, EIROP)
       ! Calculate first-order wavefunction
       ropt_mod%sdiis = .TRUE.
       ropt_mod%spcg  = .TRUE.
       CALL opt_lr (c0, c1, c2, sc0, eind, eigv, drhoe,&
            H1NL, EIVPS1, EIROP1, DDXC, VPP, CSCR, CR,&
            PSI, NSTATE, "PHONON", ORBTYP)
       ! Calculate first order density
       CALL rho1ofr (c0, c1, crge%f(:,1), drhoe, psi(:,1), nstate)
       ! Hartree and XC first-order terms
       CALL v1ofrho1 (e2, drhoe, ddxc, psi)
       ! Coupling component for each pair of surfaces
       fcsum = 0.0_real_8
       DO ic = 1,nsurf
          is1 = isurf(1,ic)
          is2 = isurf(2,ic)
          IF (ic.GT.1) THEN
             CALL zeroing(rhoov(:,1))!, nnr1)
             CALL rhoabofr (1, c0(:,is1:is1), c0(:,is2:is2), rhoov(:,1), psi(:,1))
             deinv = deinvv(ic)
          ENDIF
          IF (lspin2%tlse) THEN
             fcomp = 0.5_real_8 * (&
                  DDOT (fpar%nnr1, RHOOV, 1, DRHOE(1,2), 1) +&
                  DDOT (fpar%nnr1, RHOOV, 1, DRHOE(1,3), 1))
          ELSE IF (cntl%tlsd.AND.is2.GT.spin_mod%nsup) THEN
             fcomp = ddot (fpar%nnr1, rhoov, 1, drhoe(1,2), 1)
          ELSE
             fcomp = ddot (fpar%nnr1, rhoov, 1, drhoe(1,1), 1)
          ENDIF
          fcomp = fcomp * parm%omega/(spar%nr1s*spar%nr2s*spar%nr3s) * deinv
          CALL mp_sum(fcomp,parai%allgrp)
          fcsrf(ic) = fcomp
          fcsum = fcsum + csurf(ic)*fcomp
          IF ((lprint).AND.paral%io_parent)&
               WRITE(6,'(A,I3,A,F12.5)')&
               ' RESPONSE ALONG VECTOR ', IV, ': ', FCOMP
       ENDDO                ! IC = 1,NSURF
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
          CALL opt_lr (c0, c1, c2, sc0, eind, eigv, drhoe,&
               H1NL, EIVPS1, EIROP1, DDXC, VPP, CSCR, CR,&
               PSI, NSTATE, TYPE, ORBTYP)
          CALL rho1ofr (c0, c1, crge%f(:,1), drhoe, psi(:,1), nstate)
          ! Hartree and XC first-order terms
          CALL v1ofrho1 (e2, drhoe, ddxc, psi)
          ! Coupling component for each pair of surfaces
          DO ic = 1,nsurf
             is1 = isurf(1,ic)
             is2 = isurf(2,ic)
             IF (ic.GT.1) THEN
                CALL zeroing(rhoov(:,1))!, nnr1)
                CALL rhoabofr (1, c0(:,is1:is1), c0(:,is2:is2), rhoov(:,1), psi(:,1))
                deinv = deinvv(ic)
             ENDIF
             IF (lspin2%tlse) THEN
                fcomp = 0.5_real_8 * (&
                     DDOT (fpar%nnr1, RHOOV, 1, DRHOE(1,2), 1) +&
                     DDOT (fpar%nnr1, RHOOV, 1, DRHOE(1,3), 1))
             ELSE IF (cntl%tlsd.AND.is2.GT.spin_mod%nsup) THEN
                fcomp = ddot (fpar%nnr1, rhoov, 1, drhoe(1,2), 1)
             ELSE
                fcomp = ddot (fpar%nnr1, rhoov, 1, drhoe(1,1), 1)
             ENDIF
             fcomp = fcomp * parm%omega/(spar%nr1s*spar%nr2s*spar%nr3s) * deinv
             CALL mp_sum(fcomp,parai%allgrp)
             fcsrf(ic) = fcomp
             IF ((lprint).AND.paral%io_parent)&
                  WRITE(6,'(A,F12.5)')&
                  ' MORE ACCURATE RESPONSE: ', FCOMP
          ENDDO            ! IC = 1,NSURF
          lr02%tol_lr = tol_lr_bak
       ENDIF                ! IF (IV.EQ.1 .OR. FCSUM.GT.FC_LOW)
       ! 
       ! ==  Store response                                              ==
       ! 
       DO ic = 1,nsurf
          iat = 0
          DO isp = 1,ions1%nsp
             DO ia = 1,ions0%na(isp)
                iat = iat + 1
                DO k = 1,3
                   fc = fcsrf(ic) * vect(k,iat,iv)
                   cplvec(k,ia,isp,ic) = cplvec(k,ia,isp,ic) + fc
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
             CALL prcplngs (cplvec, cplconf, .FALSE.)
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
          ENDIF            ! (.NOT.TSPECV)
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
       CONTINUE              ! enddo ... UNTIL (IV.GE.NVECT .OR. LCONV)
    ENDIF                    ! IF (TALLDOF) ... ELSE ...
    ! ==--------------------------------------------------------------==
    ! ==  Brute force and iterative scheme: translational invariance  ==
    ! ==--------------------------------------------------------------==
    IF (.TRUE.) THEN
       CALL cplconfig (cplvec, cplconf)
       IF (lprint) THEN
          IF (paral%io_parent)&
               WRITE(6,*) 'COUPLINGS, NOT CORRECTED FOR TRANSLATIONAL ',&
               'INVARIANCE:'
          CALL prcplngs (cplvec, cplconf, .FALSE.)
       ENDIF
       ! 
       ! ==  Analytically minimize deviation from invariance (for C)     ==
       ! 
       DO k = 1,3
          v(k) = 0.0_real_8
          vc(k) = 0.0_real_8
       ENDDO
       iat = 0
       DO isp = 1,ions1%nsp
          DO ia = 1,ions0%na(isp)
             iat = iat + 1
             v(1) = v(1) + cplconf(1,ia,isp)
             v(2) = v(2) + cplconf(2,ia,isp)
             v(3) = v(3) + cplconf(3,ia,isp)
             vc(1) = vc(1) + vect(1,iat,1)
             vc(2) = vc(2) + vect(2,iat,1)
             vc(3) = vc(3) + vect(3,iat,1)
          ENDDO
       ENDDO
       c = -ddot (3, vc, 1, v, 1) / ddot (3, vc, 1, vc, 1)
       IF ((lprint).AND.paral%io_parent)&
            WRITE(6,'(A,F9.6)')&
            ' COEFFICIENT ENSURING TRANSLATIONAL INVARIANCE: ', C
       ! 
       ! == Apply correction ==
       ! 
       DO ic = 1,nsurf
          iat = 0
          DO isp = 1,ions1%nsp
             DO ia = 1,ions0%na(isp)
                iat = iat + 1
                DO k = 1,3
                   cplvec(k,ia,isp,ic) = cplvec(k,ia,isp,ic) +&
                        C*VECT(K,IAT,1)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       IF (lprint) THEN
          IF (paral%io_parent)&
               WRITE(6,*) 'COUPLINGS, CORRECTED ALONG CONSTANT DENSITY ',&
               'TERM:'
          CALL prcplngs (cplvec, cplconf, .FALSE.)
       ENDIF
       ! 
       ! ==  Undo correction                                             ==
       ! 
       DO ic = 1,nsurf
          iat = 0
          DO isp = 1,ions1%nsp
             DO ia = 1,ions0%na(isp)
                iat = iat + 1
                DO k = 1,3
                   cplvec(k,ia,isp,ic) = cplvec(k,ia,isp,ic) -&
                        C*VECT(K,IAT,1)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==

    ! if we are here TCPLLR == .true.

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

    DEALLOCATE(ddfnl,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(potr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    IF (.NOT.talldof) THEN
       ! TODO check stat
       DEALLOCATE(idone,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(fcsrf,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    ENDIF

    CALL tihalt ('    CPLANA', isub)
    ! ==--------------------------------------------------------------==
110 FORMAT (1x,a,t41,f20.8,' A.U.')
120 FORMAT (1x,a,i3,a,i3,a,t41,f20.8, ' A.U.')
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE cplana
  ! ==================================================================
  SUBROUTINE get_eind (vect, eind, eivps1, eirop1, h1nl,&
       TAU0, IESR, RHOO, PSI, NSTATE, EIROP)
    ! ==--------------------------------------------------------------==
    ! ==  Linear response along unit vector VECT                      ==
    ! ==  Input:  VECT                                                ==
    ! ==  Output: EIND, EIVPS1, EIROP1, H1NL                          ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: vect(3,ions1%nat), eind
    COMPLEX(real_8)                          :: eivps1(ncpw%nhg), &
                                                eirop1(ncpw%nhg)
    REAL(real_8)                             :: tau0(:,:,:)
    INTEGER                                  :: iesr
    REAL(real_8)                             :: rhoo(*)
    COMPLEX(real_8)                          :: psi(:,:)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: h1nl(nkpt%ngwk*nstate), &
                                                EIROP(ncpw%nhg)

    CHARACTER(*), PARAMETER                  :: procedureN = 'get_eind'

    COMPLEX(real_8), ALLOCATABLE             :: eirop1_k(:), eivps1_k(:), &
                                                h1nl_k(:)
    INTEGER                                  :: ia, iat, ierr, isp, isub, k
    LOGICAL                                  :: debug
    REAL(real_8)                             :: comp_k, eind_k

! ==--------------------------------------------------------------==

    debug = .FALSE.
    CALL tiset ('  GET_EIND', isub)

    ALLOCATE(eivps1_k(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(eirop1_k(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(h1nl_k(nkpt%ngwk*nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    eind = 0.0_real_8
    CALL zeroing(eivps1)!, nhg)
    CALL zeroing(eirop1)!, nhg)
    CALL zeroing(h1nl)!, nkpt%ngwk*nstate)
    iat = 0
    DO isp = 1,ions1%nsp
       DO ia = 1,ions0%na(isp)
          iat = iat + 1
          DO k = 1,3
             eind_k = 0.0_real_8
             CALL eind_ii (tau0, isp, ia, k, eind_k, iesr)
             CALL eind_loc (eind_k, isp, iat, k, rhoo, psi(:,1), eirop)
             ! 
             CALL fnl_set ("SWITCH")
             CALL nl_res (isp, iat, k, h1nl_k, nstate)
             CALL eind_nl (eind, isp, iat, k, nstate)
             CALL fnl_set ("SWITCH")
             ! 
             CALL eicalc1 (k, isp, iat, eivps1_k, eirop1_k)
             ! 
             comp_k = vect(k,iat)
             eind = eind + comp_k*eind_k
             CALL daxpy (2*ncpw%nhg, comp_k, eivps1_k, 1, eivps1, 1)
             CALL daxpy (2*ncpw%nhg, comp_k, eirop1_k, 1, eirop1, 1)
             CALL daxpy (2*nkpt%ngwk*nstate, comp_k, h1nl_k, 1, h1nl, 1)
          ENDDO
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    DEALLOCATE(eivps1_k,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(eirop1_k,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(h1nl_k,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL tihalt ('  GET_EIND', isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE get_eind
  ! ==================================================================
  SUBROUTINE cplconfig (cplvec, cplconf)
    ! ==--------------------------------------------------------------==
    ! ==        ANALYTIC CALCULATION OF INTERSURFACE COUPLINGS        ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: cplvec(:,:,:,:), &
                                                cplconf(:,:,:)

    INTEGER                                  :: ia, ic, isp, k
    REAL(real_8)                             :: r

    DO isp = 1,ions1%nsp
       DO ia = 1,ions0%na(isp)
          DO k = 1,3
             r = 0.0_real_8
             DO ic = 1,nsurf
                r = r + csurf(ic)*cplvec(k,ia,isp,ic)
             ENDDO
             cplconf(k,ia,isp) = r
          ENDDO
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE cplconfig
  ! ==================================================================
  SUBROUTINE prcplngs (cplvec, cplconf, ttrans)
    ! ==--------------------------------------------------------------==
    ! ==        PRINT NON-ADIABATIC COUPLING VECTOR TO STDOUT         ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: cplvec(:,:,:,:), &
                                                cplconf(:,:,:)
    LOGICAL                                  :: ttrans

    LOGICAL, PARAMETER                       :: tfor = .TRUE.

    INTEGER                                  :: i, ia, iat, ic, isp, k
    REAL(real_8) :: cplwrk(3,maxsys%nax,maxsys%nsx,nsurf), v(3)

    IF (.NOT. tcpl) RETURN
    IF (.NOT. paral%parent) RETURN
    ! 
    ! ==  Copy to CPLWRK applying or not applying uniform correction  ==
    ! 
    DO ic = 1,nsurf
       DO k = 1,3
          v(k) = 0.0_real_8
       ENDDO
       IF (ttrans) THEN      ! Uniform correction for translation
          iat = 0
          DO isp = 1,ions1%nsp
             DO ia = 1,ions0%na(isp)
                iat = iat + 1
                v(1) = v(1) + cplvec(1,ia,isp,ic)
                v(2) = v(2) + cplvec(2,ia,isp,ic)
                v(3) = v(3) + cplvec(3,ia,isp,ic)
             ENDDO
          ENDDO
          DO k = 1,3
             v(k) = v(k) / iat
          ENDDO
       ENDIF                ! IF (TTRANS)
       DO isp = 1,ions1%nsp
          DO ia = 1,ions0%na(isp)
             cplwrk(1,ia,isp,ic) = cplvec(1,ia,isp,ic) - v(1)
             cplwrk(2,ia,isp,ic) = cplvec(2,ia,isp,ic) - v(2)
             cplwrk(3,ia,isp,ic) = cplvec(3,ia,isp,ic) - v(3)
          ENDDO
       ENDDO
    ENDDO
    ! 
    CALL cplconfig (cplwrk, cplconf)
    ! 
    IF ((ttrans).AND.paral%io_parent)&
         WRITE(6,*) 'PRINTING COUPLINGS APPLYING ',&
         'TRANSLATIONAL INVARIANCE USING CONSTANT SHIFT'
    DO ic = 1,nsurf
       IF ((nsurf.GT.1).AND.paral%io_parent)&
            WRITE(6,*) 'PAIR ', ic
       IF (paral%io_parent)&
            WRITE(6,*) 'COUPLING BETWEEN KS STATES ',&
            ISURF(1,IC), ' AND ', ISURF(2,IC)
       i = 0
       DO isp = 1,ions1%nsp
          DO ia = 1,ions0%na(isp)
             i = i + 1
             IF (paral%io_parent)&
                  WRITE(6,'(T3,I4,3F12.6)')&
                  I,(CPLWRK(K,IA,ISP,IC), K=1,3)
          ENDDO
       ENDDO
    ENDDO
    IF (nsurf.GT.1) THEN
       IF (paral%io_parent)&
            WRITE(6,*) 'COUPLING BETWEEN CONFIGURATIONS'
       i = 0
       DO isp = 1,ions1%nsp
          DO ia = 1,ions0%na(isp)
             i = i + 1
             IF (paral%io_parent)&
                  WRITE(6,'(T3,I4,3F12.6)')&
                  I,(CPLCONF(K,IA,ISP), K=1,3)
          ENDDO
       ENDDO
    ENDIF
    ! 
    CALL cplconfig (cplvec, cplconf)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE prcplngs
  ! ==================================================================
  SUBROUTINE cpl_wreigen (eigv, nstate)
    ! ==--------------------------------------------------------------==
    ! ==       PRINT ADIABATIC KS ENERGIES FROM DIAGONALIZATION       ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    REAL(real_8)                             :: eigv(nstate,nkpt%nkpnt)

    INTEGER                                  :: i, is1, is2

    IF (.NOT.paral%parent) RETURN
    IF (tcpl) THEN
       IF (paral%io_parent)&
            WRITE(6,'(A,A)') ' ENERGIES OF KS STATES RELEVANT FOR ',&
            'COUPLING VECTOR GIVEN IN AU'
       DO i = 1,nsurf
          is1 = isurf(1,i)
          is2 = isurf(2,i)
          IF (paral%io_parent)&
               WRITE(6,fmt=10) is1, eigv(is1,1)
          IF (paral%io_parent)&
               WRITE(6,fmt=10) is2, eigv(is2,1)
       ENDDO
    ELSE
       IF (paral%io_parent)&
            WRITE(6,'(A)') ' ENERGIES OF ALL KS STATES GIVEN IN AU'
       DO is1 = 1,crge%n
          IF (paral%io_parent)&
               WRITE(6,fmt=10) is1, eigv(is1,1)
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
10  FORMAT (3x,'KS ENERGY OF STATE ', i3, ': ', f12.7)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE cpl_wreigen
  ! ==================================================================
  SUBROUTINE cplngs_init
    ! ==--------------------------------------------------------------==
    ! == Block data: init variables for couplings                     ==
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    tcpl       = .FALSE.  ! Don't calculate couplings by default
    tcplfd     = .FALSE.  ! Analytic instead of finite diff.
    tallat     = .TRUE.   ! Do FD couplings with all atoms
    tcpllr     = .FALSE.  ! Don't use linear response for cpl.
    talldof    = .FALSE.  ! Don't do brute force LR over all DOF
    tspecv     = .FALSE.  ! Don't specify vectors for LR
    nsurf      = 0        ! Coupling between no KS orbitals
    natfd      = 0        ! Not required if (TALLAT)
    nvect      = 5        ! Maximum number of vectors
    eps_c      = 5.0e-3_real_8   ! Finite difference displacement
    tolcpl     = 0.0_real_8    ! Tolerance for iterative LR couplings
    fc_low     = 0.05_real_8   ! Threshold for increased LR accuracy
    fc_med     = 0.1_real_8    ! Threshold for increased LR accuracy
    fc_high    = 0.5_real_8    ! Threshold for increased LR accuracy
    f_low      = 2.0_real_8    ! Increased LR accuracy
    f_med      = 4.0_real_8    ! Increased LR accuracy
    f_high     = 8.0_real_8    ! Increased LR accuracy
    ! ==--------------------------------------------------------------==
  END SUBROUTINE cplngs_init
  ! ==================================================================
  SUBROUTINE cpl_para
    ! ==--------------------------------------------------------------==
    ! ==  SETUP AND DISTRIBUTION OF PARAMETERS FOR COUPLING (SYSIN)   ==
    ! ==--------------------------------------------------------------==
    ! Variables
    CHARACTER(*), PARAMETER                  :: procedureN = 'cpl_para'

    INTEGER                                  :: i, ierr
    REAL(real_8)                             :: c

! ==--------------------------------------------------------------==
! 
! == Distribution of parameters for non-adiabatic couplings       ==
! 

    CALL mp_bcast (tcpl, parai%source, parai%cp_grp)
    CALL mp_bcast (tcplfd, parai%source, parai%cp_grp)
    CALL mp_bcast (nsurf, parai%source, parai%cp_grp)
    CALL mp_bcast (eps_c,  parai%source, parai%cp_grp)
    CALL mp_bcast (tolcpl,  parai%source, parai%cp_grp)
    CALL mp_bcast (fc_low,  parai%source, parai%cp_grp)
    CALL mp_bcast (fc_med,  parai%source, parai%cp_grp)
    CALL mp_bcast (fc_high,  parai%source, parai%cp_grp)
    CALL mp_bcast (f_low,  parai%source, parai%cp_grp)
    CALL mp_bcast (f_med,  parai%source, parai%cp_grp)
    CALL mp_bcast (f_high,  parai%source, parai%cp_grp)
    CALL mp_bcast (tcplfd, parai%source, parai%cp_grp)
    CALL mp_bcast (tallat, parai%source, parai%cp_grp)
    CALL mp_bcast (tcpllr, parai%source, parai%cp_grp)
    CALL mp_bcast (talldof, parai%source, parai%cp_grp)
    CALL mp_bcast (tspecv, parai%source, parai%cp_grp)
    CALL mp_bcast (natfd, parai%source, parai%cp_grp)
    CALL mp_bcast (nvect, parai%source, parai%cp_grp)
    IF (tcpl) THEN
       IF (paral%io_parent) THEN
          IF (tcplfd) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(A)') ' CALCULATE NON-ADIABATIC COUPLINGS'
             IF (paral%io_parent)&
                  WRITE(6,'(A,T56,1PE10.4)')&
                  ' FINITE-DIFFERENCE DISPLACEMENT', ABS(EPS_C)
             IF ((eps_c.LT.0.0_real_8).AND.paral%io_parent)&
                  WRITE(6,'(A)')&
                  ' PRODUCT ANSATZ FOR FINITE-DIFFERENCE COUPLINGS'
             IF (.NOT.tallat) THEN
                IF (paral%io_parent)&
                     WRITE(6,'(A,A,T63,I3)') ' FINITE-DIFFERENCE ONLY ',&
                     ' FOR A SUBSET OF ATOMS: ', NATFD
                IF (paral%io_parent)&
                     WRITE(6,'(2X,16I4)') (iatfd(i), i=1,natfd)
             ENDIF
          ELSE
             IF (paral%io_parent)&
                  WRITE(6,'(A)')&
                  ' CALCULATE ANALYTIC NON-ADIABATIC COUPLINGS'
             IF (tcpllr) THEN
                IF (paral%io_parent)&
                     WRITE(6,'(A,A)') ' USING LINEAR RESPONSE FOR ',&
                     'DENSITY DEPENDENCE OF KS-OPERATOR'
                IF (talldof) THEN
                   IF (paral%io_parent)&
                        WRITE(6,'(A,A)') ' LINEAR RESPONSE ALONG CARTESIAN ',&
                        'COORDINATES'
                   IF (.NOT.tallat) THEN
                      IF (paral%io_parent)&
                           WRITE(6,'(A,A,T63,I3)') ' FINITE-DIFFERENCE ONLY ',&
                           ' FOR A SUBSET OF ATOMS: ', NATFD
                      IF (paral%io_parent)&
                           WRITE(6,'(2X,16I4)') (iatfd(i), i=1,natfd)
                   ENDIF
                ELSE
                   IF (paral%io_parent)&
                        WRITE(6,'(A)') ' ITERATIVE LINEAR-RESPONSE SCHEME'
                   IF (paral%io_parent)&
                        WRITE(6,'(A,T62,I4)') ' MAXIMUM NUMBER OF VECTORS: ',&
                        NVECT
                   IF ((tspecv).AND.paral%io_parent)&
                        WRITE(6,'(A)') ' VECTORS ARE SPECIFIED'
                   IF (paral%io_parent)&
                        WRITE(6,'(A,T59,F7.5)') ' TOLERANCE FOR COUPLINGS: ',&
                        TOLCPL
                   IF (fc_high.NE.0.5_real_8 .AND. f_high.NE.8.0_real_8) THEN
                      IF (paral%io_parent)&
                           WRITE(6,'(A,F7.5,A,F7.4)')&
                           ' LR TOLERANCE DIVIDED BY ', F_LOW,&
                           ' IF COMPONENT EXCEEDS ', FC_LOW
                      IF (paral%io_parent)&
                           WRITE(6,'(A,F7.5,A,F7.4)')&
                           ' LR TOLERANCE DIVIDED BY ', F_MED,&
                           ' IF COMPONENT EXCEEDS ', FC_MED
                      IF (paral%io_parent)&
                           WRITE(6,'(A,F7.5,A,F7.4)')&
                           ' LR TOLERANCE DIVIDED BY ', F_HIGH,&
                           ' IF COMPONENT EXCEEDS ', FC_HIGH
                   ENDIF
                ENDIF
             ELSE
                IF (paral%io_parent)&
                     WRITE(6,'(A)') ' NEGLECTING RESPONSE OF DENSITY'
             ENDIF
          ENDIF
          IF (paral%io_parent)&
               WRITE(6,'(A,I1,A)') ' CONTRIBUTIONS FROM ', nsurf,&
               ' PAIRS OF KS STATES'
          c = 0.0_real_8
          DO i = 1,nsurf
             c = c + csurf(i)
             IF (paral%io_parent)&
                  WRITE(6,'(A,I3,A,I3,A,T59,F7.5)') ' STATES ',&
                  ISURF(1,I), ' - ', ISURF(2,I), ': ', CSURF(I)
          ENDDO
          IF (paral%io_parent)&
               WRITE(6,'(A,T59,F7.5)') ' SUM OF COEFFICIENTS: ', c
          IF (paral%io_parent)&
               WRITE(6,*)
       ELSE
          ALLOCATE(isurf(2,nsurf),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(csurf(nsurf),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          IF (.NOT.tallat)  THEN
             ALLOCATE(iatfd(natfd),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
          ENDIF
       ENDIF
       CALL mp_bcast (isurf, 2*nsurf, parai%source, parai%cp_grp)
       CALL mp_bcast (csurf, nsurf, parai%source, parai%cp_grp)
       CALL mp_bcast (iatfd, natfd, parai%source, parai%cp_grp)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE cpl_para
  ! ==================================================================
  SUBROUTINE rnlcpl(cplvec)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==  THE NON-LOCAL POTENTIAL CONTRIBUTION TO THE FORCE ON THE    ==
    ! ==  IONIC DEGREES OF FREEDOM                                    ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: cplvec(:,:,:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'rnlcpl'

    INTEGER                                  :: ia, iatl, ic, ierr, ii, is, &
                                                is1, is2, isa, isa0, isub, &
                                                iv, jv, k, ki, kj, l, l2, li, &
                                                lj
    REAL(real_8)                             :: tt
    REAL(real_8), ALLOCATABLE                :: dfab(:,:,:)

! Variables
! ==--------------------------------------------------------------==
! If no non-local components -> return.

    IF (nlm.EQ.0) RETURN
    ! ==--------------------------------------------------------------==

    CALL tiset('    RNLCPL',isub)
    ALLOCATE(dfab(ions1%nat,maxsys%nhxs,2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    DO ic=1,nsurf
       is1=isurf(1,ic)
       is2=isurf(2,ic)
       DO k=1,3
          CALL zeroing(dfab)!,2*ions1%nat*maxsys%nhxs)
          IF (is1.GE.parap%nst12(parai%mepos,1).AND.&
               is1.LE.parap%nst12(parai%mepos,2) ) THEN
             ii=is1-parap%nst12(parai%mepos,1)+1
             CALL dcopy(maxsys%nhxs*ions1%nat,dfnl(1,1,1,k,ii,1),1,dfab(1,1,1),1)
          ENDIF
          IF (is2.GE.parap%nst12(parai%mepos,1).AND.is2.LE.parap%nst12(parai%mepos,2) ) THEN
             ii=is2-parap%nst12(parai%mepos,1)+1
             CALL dcopy(maxsys%nhxs*ions1%nat,dfnl(1,1,1,k,ii,1),1,dfab(1,1,2),1)
          ENDIF
          CALL mp_sum(dfab,2*ions1%nat*maxsys%nhxs,parai%allgrp)
          isa0=0
          DO is=1,ions1%nsp
             IF (pslo_com%tvan(is)) THEN
                ! Vanderbild pp
                CALL stopgm("RNLCPL","VDB NOT IMPLEMENTED",& 
                     __LINE__,__FILE__)
             ELSEIF (sgpp1%tsgp(is)) THEN
                ! Stefan Goedecker pp
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
                               cplvec(k,ia,is,ic)=cplvec(k,ia,is,ic)-0.5_real_8*tt*(&
                                    dfab(isa,iv,1)*fnl(1,iatl,jv,is2,1)&
                                    +dfab(isa,iv,2)*fnl(1,iatl,jv,is1,1)&
                                    +dfab(isa,jv,1)*fnl(1,iatl,iv,is2,1)&
                                    +dfab(isa,jv,2)*fnl(1,iatl,iv,is1,1))
                            ELSE
                               cplvec(k,ia,is,ic)=cplvec(k,ia,is,ic)-tt*(&
                                    dfab(isa,iv,1)*fnl(1,isa,jv,is2,1)&
                                    +dfab(isa,iv,2)*fnl(1,isa,jv,is1,1)&
                                    +dfab(isa,jv,1)*fnl(1,isa,iv,is2,1)&
                                    +dfab(isa,jv,2)*fnl(1,isa,iv,is1,1))
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
                            cplvec(k,ia,is,ic)=cplvec(k,ia,is,ic)-wsg(is,iv)*&
                                 (dfab(isa,iv,1)*fnl(1,iatl,iv,is2,1)+&
                                 fnl(1,iatl,iv,is1,1)*dfab(isa,iv, 2))
                         ELSE
                            cplvec(k,ia,is,ic)=cplvec(k,ia,is,ic)-wsg(is,iv)*&
                                 (dfab(isa,iv,1)*fnl(1,isa,iv,is2,1)+&
                                 fnl(1,isa,iv,is1,1)*dfab(isa,iv,2))
                         ENDIF
                      ENDIF
                   ENDDO
                ENDDO
             ENDIF
             isa0 = isa0 + ions0%na(is)
          ENDDO
       ENDDO
    ENDDO
    DEALLOCATE(dfab,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    CALL tihalt('    RNLCPL',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rnlcpl
  ! ==================================================================
  SUBROUTINE give_scr_cplngs (lcplngs, tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lcplngs
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: lcplsub, lddxc, ldfnl, ldone, &
                                                leirop, leivps, LH1NL, LRHOE, &
                                                LVECT

! ==--------------------------------------------------------------==

    IF (tcplfd.OR..NOT.tcpllr) THEN
       tag     ='1' ! All subroutines are called by RWFOPT too
       lcplngs = 1
    ELSE
       tag     = 'CPLNGS'
       leivps = 2*ncpw%nhg
       leirop = 2*ncpw%nhg
       ldfnl  = 6*ions1%nat*maxsys%nhxs*ndfnl
       lh1nl  = 2*nkpt%ngwk*crge%n
       IF (lr03%txc_analytic) THEN
          lddxc = fpar%nnr1*(2*clsd%nlsd-1)
       ELSEIF (lr01%lopti.EQ.0 .OR. lr01%lopti.EQ.2) THEN
          lddxc = fpar%nnr1*(2*clsd%nlsd-1)
       ELSE
          lddxc = 1
       ENDIF
       CALL rhoe_psi_size (lrhoe)
       IF (talldof) THEN
          lvect = 3*ions1%nat
       ELSE
          lvect = 3*ions1%nat*nvect
       ENDIF
       ldone = 3*ions1%nat/2
       ! 
       lcplngs = 2*leivps + 2*leirop + ldfnl + lh1nl + 2*lrhoe +&
            LDDXC + NSURF + LVECT
       IF (.NOT.talldof) lcplngs = lcplngs + ldone + nsurf
       ! 
       CALL give_scr_cplsub (lcplsub, tag)
       lcplngs = lcplngs + lcplsub + 100
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_cplngs
  ! ==================================================================
  SUBROUTINE give_scr_cplsub (lcplsub, tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lcplsub
    CHARACTER(len=30)                        :: tag

    CHARACTER(*), PARAMETER                  :: procedureN = 'give_scr_cplsub'

    INTEGER :: lforcedr, lget_eind, LNL_RES, LOPT_LR, lrhoofr, lrnlsm, &
      LRNLSM_2D, LSUMMAT, lsymvec, LV1OFRHO1, NSTATE

! ==--------------------------------------------------------------==

    nstate = crge%n
    IF (tcplfd.OR..NOT.tcpllr) THEN
       lcplsub = 0         ! CPLFD: same scratch space as RWFOPT
    ELSE
       CALL stopgm (procedureN, 'fix me',&
            __LINE__,__FILE__)
       !vw this is buggus     call give_scr_forcedr (lforcedr,tag,nstate,.false.)
       lforcedr = 0!vw need to set that to something
       CALL give_scr_rnlsm (lrnlsm, tag, nstate, .TRUE.)
       CALL give_scr_symvec (lsymvec, tag)
       lcplsub = MAX (lforcedr, lrnlsm, lsymvec)
       CALL give_scr_rhoofr (lrhoofr, tag)
       CALL give_scr_v1ofrho1 (lv1ofrho1,tag)
       CALL give_scr_rnlsm_2d (lrnlsm_2d, tag, nstate)
       CALL give_scr_nl_res (lnl_res, nstate, tag)
       CALL give_scr_opt_lr (lopt_lr, "PHONON", tag)
       IF (talldof) THEN
          lget_eind = 0
       ELSE
          CALL give_scr_get_eind (lget_eind, tag)
       ENDIF
       IF (lspin2%tlse) THEN
          CALL give_scr_summat (lsummat, tag, nstate)
       ELSE
          lsummat = 0
       ENDIF
       lcplsub = MAX (lcplsub, lrhoofr, lrnlsm_2d, lnl_res, lopt_lr,&
            LGET_EIND, LSUMMAT, LV1OFRHO1)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_cplsub
  ! ==================================================================
  SUBROUTINE give_scr_get_eind (lget_eind, tag)
    ! ==--------------------------------------------------------------==
    ! Arguments and local variables
    INTEGER                                  :: lget_eind
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: lnl_res

! ==--------------------------------------------------------------==

    CALL give_scr_nl_res (lnl_res, crge%n, tag)
    lget_eind = 4*ncpw%nhg + 2*nkpt%ngwk*crge%n + 100
    tag = '4*NHG+2*NGWK*N+100+LNL_RES'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_get_eind
  ! ==================================================================
  SUBROUTINE switch_st (state, f, nstate, tlse, tlsets)
    ! ==--------------------------------------------------------------==
    CHARACTER(len=3)                         :: state
    INTEGER                                  :: nstate
    REAL(real_8)                             :: f(nstate)
    LOGICAL                                  :: tlse, tlsets

    LOGICAL, SAVE                            :: tlse_def, tlsets_def
    REAL(real_8), SAVE                       :: f_def(2)

! ==--------------------------------------------------------------==

    IF (state.EQ.'SET') THEN
       tlse_def = tlse
       tlsets_def = tlsets
       f_def(nstate-1) = f(nstate-1)
       f_def(nstate) = f(nstate)
    ELSE IF (state.EQ.'DEF') THEN
       tlse = tlse_def
       tlsets = tlsets_def
       f(nstate-1) = f_def(nstate-1)
       f(nstate) = f_def(nstate)
    ELSE IF (state.EQ.'GND') THEN
       tlse = .FALSE.
       tlsets = .FALSE.
       f(nstate-1) = 2.0_real_8
       f(nstate) = 0.0_real_8
    ELSE IF (state.EQ.'EXC') THEN
       tlse = .TRUE.
       tlsets = .FALSE.
       f(nstate-1) = 1.0_real_8
       f(nstate) = 1.0_real_8
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE switch_st
  ! ==================================================================
  SUBROUTINE tdrhof_dis(fion,rho,psi)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: fion(:,:,:), rho(*)
    COMPLEX(real_8)                          :: psi(:,:)

    COMPLEX(real_8)                          :: ei123, gx, gy, gz, rhets, &
                                                txx, tyy, tzz, vcgs
    INTEGER                                  :: ia, ig, ig1, ir, is, isa
    REAL(real_8)                             :: omtp

! Variables
! &           EI3(NATX,(2*NR3S-1))
! ==--------------------------------------------------------------==

    !$omp parallel do private(IR)
    DO ir=1,fpar%nnr1
       psi(ir,1)=CMPLX(rho(ir),0._real_8,kind=real_8)
    ENDDO
    CALL  fwfftn(psi(:,1),.FALSE.,parai%allgrp)
    ! 
    omtp=2._real_8*parm%omega*parm%tpiba
    ig1=1
    IF (geq0) ig1=2
    DO ig=ig1,ncpw%nhg
       rhets=CONJG(psi(nzh(ig),1))
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
    RETURN
  END SUBROUTINE tdrhof_dis
  ! ==================================================================
  SUBROUTINE cplrdvec (vect, nat, nvec, tfdens)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nat, nvec
    REAL(real_8)                             :: vect(3,nat,nvec)
    LOGICAL                                  :: tfdens

    CHARACTER(len=80)                        :: line
    INTEGER                                  :: iat, ierr, ifirst, iunit, &
                                                ivect, k
    REAL(real_8)                             :: rlen

! ==--------------------------------------------------------------==

    iunit = 5
    IF (paral%parent) THEN
       ierr = inscan (iunit,'&VECTORS')
       ivect = 0
       IF (paral%io_parent)&
            READ(iunit,*,END=99) line
       tfdens = (INDEX(line,'FIXED').NE.0)
       ifirst = 1
       IF (tfdens) ifirst = 2
       DO ivect = ifirst,nvec
          IF (paral%io_parent)&
               READ(iunit,*,END=99,err=98)&
               ((VECT(K,IAT,IVECT), K=1,3), IAT=1,NAT)
          rlen = 0.0_real_8
          DO iat = 1,nat
             DO k = 1,3
                rlen = rlen + vect(k,iat,ivect)**2
             ENDDO
          ENDDO
          rlen = 1.0_real_8 / SQRT(rlen)
          DO iat = 1,nat
             vect(1,iat,ivect) = vect(1,iat,ivect)*rlen
             vect(2,iat,ivect) = vect(2,iat,ivect)*rlen
             vect(3,iat,ivect) = vect(3,iat,ivect)*rlen
          ENDDO
       ENDDO
    ENDIF
    CALL mp_bcast (vect, 3*nat*nvec, parai%source, parai%cp_grp)
    CALL mp_bcast (tfdens, parai%source, parai%cp_grp)
    ! ==--------------------------------------------------------------==
    RETURN
98  IF (paral%io_parent)&
         WRITE(6,*) 'ERROR READING RESPONSE VECTOR ', ivect
    CALL stopgm ('CPLRDVEC', 'ERROR IN INPUT',& 
         __LINE__,__FILE__)
99  IF (paral%io_parent)&
         WRITE(6,*) 'END OF FILE REACHED READING RESPONSE VECTOR ', ivect
    CALL stopgm ('CPLRDVEC', 'END OF FILE',& 
         __LINE__,__FILE__)
  END SUBROUTINE cplrdvec
  ! ==================================================================

  SUBROUTINE cplfd (c0, c2, cf, c0ini, cr, sc0, cscr, vpp, eigv,&
       RHOE, RHOINI, PSI, &
       TAU0, FION, IFCALC, IREC, TINFO, TFOR,&
       CPLVEC, NSTATE)
    ! ==--------------------------------------------------------------==
    ! ==   FINITE DIFFERENCE CALCULATION OF INTERSURFACE COUPLINGS    ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c2(nkpt%ngwk,*), CR(*), &
                                                SC0(*), CSCR(*)
    REAL(real_8)                             :: vpp(:), eigv(*), RHOE(:,:), &
                                                RHOINI(fpar%nnr1,*)
    COMPLEX(real_8)                          :: PSI(:,:)
    REAL(real_8)                             :: tau0(:,:,:), fion(:,:,:)
    INTEGER                                  :: ifcalc, irec(:)
    LOGICAL                                  :: tinfo, tfor
    REAL(real_8)                             :: CPLVEC(:,:,:,:)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: C0INI(nkpt%ngwk,NSTATE,*), &
                                                CF(nkpt%ngwk,NSTATE,1), &
                                                c0(nkpt%ngwk,nstate,1)

    INTEGER                                  :: ia, iact, ic, idf, ik, is1, &
                                                is2, isa, isp, JA, JACT, JDF, &
                                                JK, JSA, JSP, K, NC0, NRHO
    LOGICAL                                  :: lprint
    REAL(real_8)                             :: eps_c2n, etotf, fcomp, r02, &
                                                R10, R11, R12, r1m2, r1p2, &
                                                R22, r2m2, r2p2

! Variables
! TODO refactor common block

    COMMON     /cplscr/ r02, r1m2, r1p2, r2m2, r2p2, r10, r11, r12,&
         R22

    COMPLEX(real_8) :: z12

    REAL(real_8), ALLOCATABLE :: oldf(:,:,:)

    ! Externals
    REAL(real_8), EXTERNAL :: dotp
    COMPLEX(real_8), EXTERNAL :: zdotc
    INTEGER :: ierr
    CHARACTER(*),PARAMETER :: procedureN='cplfd'
    ! ==--------------------------------------------------------------==
    lprint  = .TRUE.
    CALL forces_diag(nstate, c0, c2, cr, sc0, cscr, vpp, eigv,&
         RHOE, PSI,&
         TAU0, TAU0, TAU0, FION, IFCALC,&
         IREC, TFOR, TINFO)
    IF (paral%parent.AND.lprint) THEN
       IF (paral%io_parent)&
            WRITE(6,*)&
            'FINITE-DIFFERENCE COUPLING VECTOR, INITIAL STATE'
       CALL cpl_wreigen (eigv, nstate)
       IF (paral%io_parent)&
            WRITE(6,'(A,F12.5)') ' INITIAL POINT, ENERGY: ', ener_com%etot
    ENDIF
    IF (cntl%tdiag.AND.cntl%tdavi) THEN
       nc0     = 2*nkpt%ngwk*cnti%ndavv*nkpt%nkpnt + 8
    ELSE
       nc0     = 2*nkpt%ngwk*nstate*nkpt%nkpnt + 8
    ENDIF
    nrho    = fpar%nnr1*clsd%nlsd
    eps_c2n = 1.0_real_8 / (2.0_real_8 * eps_c)
    ! 
    ! ==  Save the initial converged wavefuntion and density          ==
    ! 
    CALL dcopy (nc0, c0, 1, c0ini, 1)
    IF (cntl%tdiag) CALL dcopy (nrho, rin0, 1, rhoini, 1)
    IF (paral%parent) THEN
       IF (tfor) THEN
          IF (tallat) THEN
             nvar = 3*ions1%nat
          ELSE
             nvar = 3*natfd
          ENDIF

          ALLOCATE(hesscr(nvar,nvar),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(oldf(3,maxsys%nax,maxsys%nsx),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)

          CALL zeroing(hesscr)!, nvar*nvar)
          IF (paral%io_parent)&
               WRITE(6,'(A,A,I3,A)') ' CALCULATING FINITE DIFFERENCE ',&
               'HESSIAN FOR ', NVAR, ' DOF'
       ENDIF
    ENDIF
    ! 
    ! ==  Product (constant) density if EPS_C<0                       ==
    ! 
    IF (eps_c.LT.0.0_real_8) THEN
       eps_c = -eps_c
       eps_c2n = -eps_c2n
       cnti%nomore = 0
    ENDIF
    ! 
    ! ==  Main finite difference loop over the dimensions             ==
    ! 
    isa = 0
    idf = 0
    iact = 1
    DO isp = 1,ions1%nsp
       DO ia = 1,ions0%na(isp)
          isa = isa + 1
          IF (tallat.OR.isa.EQ.iatfd(iact)) THEN
             iact = iact + 1
          ELSE
             GOTO 9
          ENDIF
          DO k = 1,3
             idf = idf + 1
             IF (paral%parent) THEN
                IF ((lprint).AND.paral%io_parent)&
                     WRITE(6,'(A,3I3)')&
                     ' TYPE, ATOM, DIMENSION: ', ISP, IA, K
             ENDIF
             ! 
             ! ==  Finite difference step backward                              ==
             ! 
             tau0(k,ia,isp) = tau0(k,ia,isp) - eps_c
             CALL phfac (tau0)
             CALL dcopy (nc0, c0ini, 1, cf, 1)
             IF (cntl%tdiag) CALL dcopy (nrho, rhoini, 1, rin0, 1)
             CALL forces_diag(nstate, cf, c2, cr, sc0, cscr, vpp, eigv,&
                  RHOE, PSI,&
                  TAU0, TAU0, TAU0, FION, IFCALC,&
                  IREC, TFOR, TINFO)
             etotf = ener_com%etot
             IF (paral%parent) THEN
                IF (tfor) CALL dcopy (3*maxsys%nsx*maxsys%nax, fion, 1, oldf, 1)
                IF ((lprint).AND.paral%io_parent)&
                     WRITE(6,'(A,F12.5)')&
                     ' DONE BACKWARD, ENERGY: ', ener_com%etot
             ENDIF
             ! 
             ! ==  Two finite difference steps forward                        ==
             ! 
             tau0(k,ia,isp) = tau0(k,ia,isp) + 2.0_real_8*eps_c
             CALL phfac (tau0)
             CALL dcopy (nc0, c0ini, 1, c0, 1)
             IF (cntl%tdiag) CALL dcopy (nrho, rhoini, 1, rin0, 1)
             CALL forces_diag(nstate, c0, c2, cr, sc0, cscr, vpp, eigv,&
                  RHOE, PSI,&
                  TAU0, TAU0, TAU0, FION, IFCALC,&
                  IREC, TFOR, TINFO)
             IF (paral%parent) THEN
                IF (tfor) THEN
                   CALL daxpy (3*maxsys%nax*maxsys%nsx, -1.0_real_8, fion, 1, oldf, 1)
                   CALL dscal (3*maxsys%nax*maxsys%nsx, 0.5_real_8*eps_c2n, oldf, 1)
                   jsa = 0
                   jdf = 0
                   jact = 1
                   DO jsp = 1,ions1%nsp
                      DO ja = 1,ions0%na(jsp)
                         jsa = jsa + 1
                         IF (tallat.OR.jsa.EQ.iatfd(jact)) THEN
                            jact = jact + 1
                            DO jk = 1,3
                               jdf = jdf + 1
                               hesscr(idf+nvar*(jdf-1),1) &
                                    =HESSCR(IDF+NVAR*(JDF-1),1) + OLDF(JK,JA,JSP&
                                    )
                               hesscr(jdf+nvar*(idf-1),1) &
                                    =HESSCR(JDF+NVAR*(IDF-1),1) + OLDF(JK,JA,JSP&
                                    )
                            ENDDO
                         ENDIF
                      ENDDO
                   ENDDO
                ENDIF
                IF ((lprint).AND.paral%io_parent)&
                     WRITE(6,'(A,F12.5)')&
                     ' DONE FORWARD, ENERGY: ', ener_com%etot
             ENDIF
             ! 
             DO ic = 1,nsurf
                is1 = isurf(1,ic)
                is2 = isurf(2,ic)
                IF (tkpts%tkpnt) THEN
                   DO ik = 1,nkpt%nkpnt
                      z12 = zdotc (nkpt%ngwk, c0(1,is1,ik), 1, cf(1,is2,ik), 1)
                      CALL mp_sum(z12,parai%allgrp)
                   ENDDO
                ELSE
                   r02 = dotp (ncpw%ngw, c0ini(1,is2,1), cf(1,is2,1))
                   r10 = dotp (ncpw%ngw, c0(1,is1,1), c0ini(1,is1,1))
                   r11 = dotp (ncpw%ngw, c0(1,is1,1), cf(1,is1,1))
                   r12 = dotp (ncpw%ngw, c0(1,is1,1), cf(1,is2,1))
                   r22 = dotp (ncpw%ngw, c0(1,is2,1), cf(1,is2,1))
                   r1p2 = dotp (ncpw%ngw, c0ini(1,is1,1), cf(1,is2,1))
                   r2p2 = dotp (ncpw%ngw, c0ini(1,is2,1), cf(1,is2,1))
                   r1m2 = dotp (ncpw%ngw, c0ini(1,is1,1), c0(1,is2,1))
                   r2m2 = dotp (ncpw%ngw, c0ini(1,is2,1), c0(1,is2,1))
                   CALL mp_sum(r02,parai%allgrp)
                   CALL mp_sum(r10,parai%allgrp)
                   CALL mp_sum(r11,parai%allgrp)
                   CALL mp_sum(r12,parai%allgrp)
                   CALL mp_sum(r22,parai%allgrp)
                   CALL mp_sum(r1p2,parai%allgrp)
                   CALL mp_sum(r2p2,parai%allgrp)
                   CALL mp_sum(r1m2,parai%allgrp)
                   CALL mp_sum(r2m2,parai%allgrp)

                   IF (r2p2.LT.0.0_real_8) r1p2 = -r1p2
                   IF (r2m2.LT.0.0_real_8) r1m2 = -r1m2
                   IF (paral%parent.AND.lprint) THEN
                      IF (paral%io_parent)&
                           WRITE(6,'(A,I3,A,I3,A)')&
                           ' DOT-PRODUCTS FOR STATES a = ', IS1,&
                           ' AND b = ', IS2, ':'
                      IF (paral%io_parent)&
                           WRITE(6,'(2(A,F12.8,A,F12.8,/),A,F12.8)')&
                           '   <a|a''> = ', R11, '   <b|b''>  = ', R22,&
                           '   <a|a0> = ', R10, '   <b0|b''> = ', R02,&
                           '   <a|b''> = ', R12
                      IF (paral%io_parent)&
                           WRITE(6,'(A,F12.8,A,F12.8)')&
                           '   <a0|b> = ', R1M2, '   <a0|b''> = ', R1P2
                   ENDIF
                   r12 = r12 * eps_c2n
                   IF (r11*r22.LT.0.0_real_8) r12 = -r12
                   fcomp = -(etotf-ener_com%etot)*eps_c2n
                   cplvec(k,ia,isp,ic) = (r1p2 - r1m2)* eps_c2n
                   IF (paral%parent) THEN
                      IF (lprint) THEN
                         IF (paral%io_parent)&
                              WRITE(6,'(A,1PE12.5)')&
                              ' COUPLING (SCALED):     ', CPLVEC(K,IA,ISP,IC)
                         IF (paral%io_parent)&
                              WRITE(6,'(A,1PE12.5)')&
                              ' COUPLING FROM <a|b''>:  ', R12
                         IF (paral%io_parent)&
                              WRITE(6,'(A,1PE12.5)')&
                              ' F/D FORCE COMPONENT:   ', FCOMP
                      ENDIF
                   ENDIF
                   IF (.NOT.tfor) fion(k,ia,isp) = fcomp
                ENDIF
             ENDDO        ! IC = 1,NSURF
             ! 
             ! ==  Undo finite difference step                                 ==
             ! 
             tau0(k,ia,isp) = tau0(k,ia,isp) - eps_c
          ENDDO            ! K = 1,3
9         CONTINUE          ! jump here if atom is not interesting
       ENDDO                ! IA = 1,NA(ISP)
    ENDDO                    ! ISP = 1,NSP
    IF (tfor) THEN
       lvlhes = .TRUE.
       IF (paral%parent) THEN
          DEALLOCATE(oldf,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE cplfd

  ! ==================================================================

END MODULE cplngs_utils
