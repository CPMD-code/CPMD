MODULE spin
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==================================================================
  ! == FOR LOCAL SPIN DENSITY                                       ==
  ! ==--------------------------------------------------------------==
  ! == NLSD = 2 IF TLSD, OTHERWISE 1                                ==
  ! == NLSX = 3 IF TLSD, OTHERWISE 1                                ==
  ! == NLSD = 4 IF TLSE                                             ==
  ! == IALPHA = alpha f=1 orbital                                   ==
  ! == IBETA  = beta  f=1 orbital                                   ==
  ! == TDIAG options:                                               ==
  ! == TFXSP: fixed spin ( multiplicity) calculation                ==
  ! == AMUUP AMUDO: chemical potentials for up and down electrons   ==
  ! == NUPEL NDOEL: number of up and down electronc ( TDIAG)        ==
  ! ==--------------------------------------------------------------==

  TYPE :: clsd_t
     INTEGER :: nlsd = HUGE(0)
     INTEGER :: nlsx = HUGE(0)
     INTEGER :: ialpha = HUGE(0)
     INTEGER :: ibeta = HUGE(0)
  END TYPE clsd_t
  TYPE(clsd_t), SAVE :: clsd
  ! ==--------------------------------------------------------------==

  TYPE :: spin_mod_t
     INTEGER :: nsup = HUGE(0)
     INTEGER :: nsdown = HUGE(0)
     INTEGER :: nspin = HUGE(0)
  END TYPE spin_mod_t
  TYPE(spin_mod_t), SAVE :: spin_mod
  ! ==--------------------------------------------------------------==

  TYPE :: lspin1_t
     REAL(real_8) :: lsea = HUGE(0.0_real_8)
     REAL(real_8) :: lseb = HUGE(0.0_real_8)
  END TYPE lspin1_t
  TYPE(lspin1_t), SAVE :: lspin1

  TYPE :: lspin2_t
     LOGICAL :: tlse = .FALSE.
     LOGICAL :: troks = .FALSE.
     LOGICAL :: tross = .FALSE.
     LOGICAL :: tcas22 = .FALSE.
     LOGICAL :: troot = .FALSE.
     LOGICAL :: teprof = .FALSE.
     LOGICAL :: tpenal = .FALSE.
     LOGICAL :: tmgoe = .FALSE.
     LOGICAL :: TLSETS = .FALSE.
     LOGICAL :: TROS = .FALSE.
  END TYPE lspin2_t
  TYPE(lspin2_t), SAVE :: lspin2

  TYPE :: lspin3_t
     REAL(real_8) :: hablse = HUGE(0.0_real_8)
     REAL(real_8) :: haolse = HUGE(0.0_real_8)
     REAL(real_8) :: hbolse = HUGE(0.0_real_8)
     REAL(real_8) :: rotab = HUGE(0.0_real_8)
     REAL(real_8) :: mgmab = HUGE(0.0_real_8)
     REAL(real_8) :: mgmba = HUGE(0.0_real_8)
     REAL(real_8) :: mgab(12) = HUGE(0.0_real_8)
  END TYPE lspin3_t
  TYPE(lspin3_t), SAVE :: lspin3

  TYPE :: lspin4_t
     REAL(real_8) :: apenal = HUGE(0.0_real_8)
  END TYPE lspin4_t
  TYPE(lspin4_t), SAVE :: lspin4
  ! ==--------------------------------------------------------------==
  ! .. TDIAG variables for fixed spin ( two amus )


  REAL(real_8), SAVE :: amuup=0.0_real_8,amudo=0.0_real_8
  TYPE :: tdsp1_t
     INTEGER :: nupel=HUGE(0)
     INTEGER :: ndoel=HUGE(0)
     LOGICAL :: tfxsp=.TRUE.
  END TYPE tdsp1_t
  TYPE(tdsp1_t), SAVE :: tdsp1
  TYPE :: tdsp_t
     REAL(real_8) :: amuup=HUGE(0.0_real_8)
     REAL(real_8) :: amudo=HUGE(0.0_real_8)
  END TYPE tdsp_t
  TYPE(tdsp_t), SAVE :: tdsp
  ! ==================================================================

END MODULE spin
