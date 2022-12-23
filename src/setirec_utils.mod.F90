MODULE setirec_utils
  USE clas,                            ONLY: tclas
  USE cotr,                            ONLY: cotc0
  USE glemod,                          ONLY: glepar
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE mp_interface,                    ONLY: mp_bcast
  USE nose,                            ONLY: loct,&
                                             nosl
  USE parac,                           ONLY: parai
  USE shop_rest,                       ONLY: sh03
  USE store_types,                     ONLY: &
       irec_ac, irec_cell, irec_clas, irec_co, irec_cons, irec_ctrans, &
       irec_cut, irec_eigv, irec_gle, irec_he, irec_ico, irec_info, irec_kpt, &
       irec_lbfgs, irec_lrwf, irec_lrwv, irec_nas, irec_noc, irec_noe, &
       irec_nop1, irec_nop2, irec_nop3, irec_nop4, irec_nose, irec_occ, &
       irec_phes, irec_pot, irec_prfo, irec_prng, irec_pwv, irec_rdtl, &
       irec_rho, irec_sdp, irec_shop, irec_shopbo, irec_vel, irec_wf, &
       irec_xtrp, restart1, store1
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             cntr
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: read_irec
  PUBLIC :: write_irec
  PUBLIC :: isetone

CONTAINS

  ! ==================================================================
  SUBROUTINE read_irec(irec)
    ! ==--------------------------------------------------------------==
    ! ==                     GENERAL FILE FORMAT (FOR IREC ALSO)      ==
    ! ==--------------------------------------------------------------==
    ! ==  Section  1: Header                                          ==
    ! ==  Section  2: Symmetry and Cell info                          ==
    ! ==  Section  3: Number of atomic species and atoms per species  ==
    ! ==  Section  4: Atomic coordinates                              ==
    ! ==  Section  5: Atomic velocities                               ==
    ! ==  Section  6: Initial atomic coordinates                      ==
    ! ==  Section  7: Cutoff, # of electrons, grid                    ==
    ! ==  Section  8: States, dual, plane waves                       ==
    ! ==  Section  9: PW Coefficients                                 ==
    ! ==  Section 10: PW Velocities                                   ==
    ! ==  Section 11: Accumulators                                    ==
    ! ==  Section 12: Nose thermostats general info                   ==
    ! ==  Section 13: Nose thermostats electrons                      ==
    ! ==  Section 14: Nose thermostats ions                           ==
    ! ==  Section 15: Nose thermostats ions (ULTRA)                   ==
    ! ==  Section 16: Nose thermostats ions (MASSIVE)                 ==
    ! ==  Section 17: Nose thermostats cell                           ==
    ! ==  Section 18: Potential                                       ==
    ! ==  Section 19: PW single states                                ==
    ! ==  Section 20: H Matrices (cell parameters)                    ==
    ! ==  Section 21: K-Points                                        ==
    ! ==  Section 22: Electron density                                ==
    ! ==  Section 23: Occupation numbers                              ==
    ! ==  Section 24: Fermi energy and eigenvalues                    ==
    ! ==  Section 25: Classical Particles (coordinates and velocities)==
    ! ==  Section 26: LinRes PW Coefficients                          ==
    ! ==  Section 27: LinRes PW Velocities                            ==
    ! ==  Section 28: Partial Hessian (for microiterative TS search)  ==
    ! ==  Section 29: P-cntl%rfo status and work arrays                    ==
    ! ==  Section 30: L-cntl%bfgs history and status                       ==
    ! ==  Section 31: Adaptive tolerance status                       ==
    ! ==  Section 32: Constraints values                              ==
    ! ==  Section 33: Cell translation vector in QM/MM runs           ==
    ! ==  Section 34: Local temperature parameters                    ==
    ! ==  Section 35: Wavefunction history for extrapolation          ==
    ! ==  Section 36: Surface Hopping I                               ==
    ! ==  Section 37: Surface Hopping II                              ==
    ! ==  Section 38: Pseudo random number generator                  ==
    ! ==  Section 39: (generalized) Langevin equation                 ==
    ! ==  Section 40 - 99 : empty                                     ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: irec(:)

! ==--------------------------------------------------------------==

    CALL zeroing(irec)!,100)
    irec(irec_info) = 1
    irec(irec_cell) = 1
    irec(irec_nas)  = 1
    irec(irec_co)   = isetone(restart1%rco)
    irec(irec_vel)  = isetone(restart1%rvel)
    irec(irec_ico)  = 1
    irec(irec_cut)  = 1
    irec(irec_sdp)  = 1
    irec(irec_wf)   = isetone(restart1%rwf)
    ! By default, PW velocities are not read for diagonalization
    ! and BOMD.
    irec(irec_pwv)  = isetone(restart1%rwf.AND. .NOT.(cntl%tmdbo.OR.cntl%geopt.OR.&
         cntl%tdiag.OR.cntl%wfopt))
    irec(irec_ac)   = isetone(restart1%rac)
    irec(irec_nose) = isetone((restart1%rnoe.OR.restart1%rnop.OR.restart1%rnoc).AND.cntl%md)
    irec(irec_noe)  = isetone(restart1%rnoe.AND.cntl%md)
    irec(irec_nop1) = isetone(restart1%rnop.AND.cntl%md.AND..NOT.&
         (nosl%tultra.OR.nosl%tmnose.OR.loct%tloct))
    irec(irec_nop2) = isetone(restart1%rnop.AND.cntl%md.AND.nosl%tultra)
    irec(irec_nop3) = isetone(restart1%rnop.AND.cntl%md.AND.nosl%tmnose)
    irec(irec_nop4) = isetone(restart1%rnop.AND.cntl%md.AND.loct%tloct)
    irec(irec_noc)  = isetone(restart1%rnoc.AND.cntl%md.AND.cntl%tnosec)
    irec(irec_pot)  = isetone(restart1%rpot.AND.store1%spot)
    irec(irec_he)   = isetone(restart1%rcell.AND.cntl%tprcp)
    irec(irec_kpt)  = isetone(restart1%rkpt)
    irec(irec_rho)  = isetone(restart1%rrho.AND.cntl%tdiag)
    irec(irec_occ)  = isetone(restart1%rocc)
    irec(irec_eigv) = isetone(restart1%reigv)
    irec(irec_clas) = isetone(tclas.AND.restart1%rco)
    !CSOC[
    irec(irec_lrwf) = isetone((cntl%tddft.OR.cntl%tspec.OR.cntl%tsoc).AND.restart1%rlr)
    !CSOC]
    irec(irec_lrwv) = 0
    irec(irec_phes) = isetone(restart1%rphes)
    irec(irec_lbfgs)= isetone(restart1%rlscl.AND.cntl%lbfgs)
    irec(irec_prfo) = isetone(restart1%rlscl.AND.cntl%prfo)
    irec(irec_rdtl) = isetone(restart1%rdtl)
    irec(irec_cons) = isetone(restart1%rcon)
    irec(irec_ctrans) = isetone(cntl%tqmmm)
    irec(irec_xtrp) = isetone(cntl%textrap.AND.restart1%rxtp)
    irec(irec_shop) = isetone(sh03%tshopres.AND.(.NOT.cntl%tfusi))
    irec(irec_shopbo) = isetone(sh03%tshopres.AND.cntl%tmdbo.AND.(.NOT.cntl%tfusi))
    irec(irec_prng) = isetone(restart1%rprng)
    irec(irec_gle) = isetone(restart1%rgle.AND.cntl%md.AND.(glepar%gle_mode.NE.0))
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE read_irec
  ! ==================================================================
  SUBROUTINE write_irec(irec)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: irec(:)

! ==--------------------------------------------------------------==

    CALL zeroing(irec)!,100)
    irec(irec_info) = 1
    irec(irec_cell) = 1
    irec(irec_nas)  = 1
    irec(irec_co)   = 1
    irec(irec_vel)  = isetone(.NOT.cntl%ksener)
    irec(irec_cut)  = 1
    irec(irec_sdp)  = 1
    irec(irec_ico)  = 1
    irec(irec_cut)  = 1
    irec(irec_sdp)  = 1
    irec(irec_wf)   = isetone(store1%swf)
    IF (tkpts%tonlydiag) irec(irec_wf)   = 0
    ! By default, PW velocities are not stored for diagonalization
    ! and BOMD.
    irec(irec_pwv)  = isetone(.NOT.(cntl%tmdbo.OR.cntl%tdiag.OR.cntl%geopt.OR.&
         cntl%wfopt))
    ! By default, accumulators are stored for cntl%md and GEOMETRY.
    irec(irec_ac)   = isetone((cntl%geopt.OR.cntl%md).AND..NOT.store1%tdebacc)
    ! Nose information
    irec(irec_nose) = isetone(cntl%tnosee.OR.cntl%tnosep.OR.cntl%tnosec)
    irec(irec_noe)  = isetone(cntl%tnosee)
    irec(irec_nop1) = isetone(cntl%tnosep.AND.&
         .NOT.(nosl%tultra.OR.nosl%tmnose.OR.loct%tloct))
    irec(irec_nop2) = isetone(cntl%tnosep.AND.nosl%tultra)
    irec(irec_nop3) = isetone(cntl%tnosep.AND.nosl%tmnose)
    irec(irec_nop4) = isetone(cntl%tnosep.AND.loct%tloct)
    irec(irec_noc)  = isetone(cntl%tnosec)
    ! By default, the potential are stored for Kohn-Sham calculation
    irec(irec_pot)  = isetone(store1%spot)
    ! By default, the cell information are stored
    irec(irec_he)   = 1
    ! K-points
    irec(irec_kpt)  = isetone(tkpts%tkpnt)
    ! By default, density is stored for diagonalization
    irec(irec_rho)  = isetone(store1%srho)
    ! By default, the occupation numbers are stored for diagonalization
    irec(irec_occ)  = isetone(cntl%tdiag)
    ! By default, the fermi energy and eigenvalues are stored for diag.
    irec(irec_eigv) = isetone(cntl%tdiag)
    irec(irec_clas) = isetone(tclas)
    ! cntl%tddft
    !CSOC[
    irec(irec_lrwf) = isetone(cntl%tddft.OR.cntl%tspec.OR.cntl%tsoc)
    !CSOC]
    irec(irec_lrwv) = 0
    ! The partial Hessian is written whenever it is valid
    irec(irec_phes) = isetone(cntl%prfo)
    ! The linear scaling optimizers should always dump their status
    irec(irec_lbfgs)= isetone(cntl%lbfgs)
    irec(irec_prfo) = isetone(cntl%prfo)
    irec(irec_rdtl) = isetone(cntr%tolene.NE.0.0_real_8.OR.cntr%tolad.NE.0.0_real_8)
    irec(irec_cons) = isetone(cotc0%mcnstr.NE.0)
    irec(irec_ctrans) = isetone(cntl%tqmmm)
    irec(irec_xtrp) = isetone(cntl%textrap.AND.cntl%tstrxtp)
    irec(irec_shop) = isetone(cntl%tshop.AND.(.NOT.cntl%tfusi))
    irec(irec_shopbo) = isetone(cntl%tshop.AND.cntl%tmdbo.AND.(.NOT.cntl%tfusi))
    irec(irec_prng) = isetone(cnti%iprng.NE.0)
    irec(irec_gle) = isetone(glepar%gle_mode.NE.0)

    CALL mp_bcast(irec,SIZE(irec),parai%source,parai%allgrp)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE write_irec
  ! ==================================================================
  FUNCTION isetone(zlog)
    ! ==--------------------------------------------------------------==
    LOGICAL                                  :: zlog
    INTEGER                                  :: isetone

! ==--------------------------------------------------------------==

    IF (zlog) THEN
       isetone = 1
    ELSE
       isetone = 0
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION isetone
  ! ==================================================================

END MODULE setirec_utils
