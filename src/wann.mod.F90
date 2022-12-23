MODULE wann
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==================================================================
  ! == INCLUDE FILE FOR WANNIER FUNCTIONS                           ==
  ! ==--------------------------------------------------------------==
  ! == TWANN : calculate Wannier orbitals                           ==
  ! == TWPRI : print Wannier orbitals                               ==
  ! == TSDENS: print square of Wannier orbitals                     ==
  ! == TWDOS: print wannier density of states                       ==
  ! == TWMOL: print molecular orbitals                              ==
  ! ==--------------------------------------------------------------==
  INTEGER, PARAMETER :: nwannl=5 

  ! ==--------------------------------------------------------------==
  ! == W_MAXS : Maximum number of optimization steps                ==
  ! == W_TYPE : Type of optimization functional                     ==
  ! ==          1 Vanderbilt                                        == 
  ! ==          2 RESTA                                             ==
  ! == W_OPT  : Optimisation procedure                              ==
  ! ==          1 steepest descent                                  ==
  ! ==          2 Jacobi Rotation                                   ==
  ! == NS_WANN: sampling rate for Wannier function output           ==
  ! == SW_FIRST : first Wannier orbital to be printed               ==
  ! == SW_LAST  : last Wannier orbital to be printed                ==
  ! == SW_ALL   : all Wannier orbital to be printed                 ==
  ! == SW_ORB   : number of orbitals to be printed                  ==
  ! == SW_LIST  : list of orbitals to be printed                    ==
  ! ==--------------------------------------------------------------==
  INTEGER, PARAMETER :: nwanni=8 
  INTEGER, ALLOCATABLE, SAVE :: sw_list(:)

  ! ==--------------------------------------------------------------==
  ! == W_STEP : step size for optimisation                          ==
  ! == W_EPS  : convergence criteria                                ==
  ! == W_RAN  : randomization amplitude                             ==
  ! == W_REF  : reference point for Wannier centres                 ==
  ! == SW_SPREAD: spread threshold for orbitals to be printed       ==
  ! ==--------------------------------------------------------------==
  INTEGER, PARAMETER :: nwannr=7 

  ! ==--------------------------------------------------------------==
  ! == NWANOPT : number of operators                                ==
  ! == IOW     : operators 4-6                                      ==
  ! == WWEI    : weights                                            ==
  ! == WQUADI  : weights for quadrupoles (RV)                       ==
  ! ==--------------------------------------------------------------==
  TYPE :: wannc_t
     INTEGER :: nwanopt=HUGE(0)
     INTEGER :: iow(2,3)=HUGE(0)
     REAL(real_8) :: wwei(6)=HUGE(0.0_real_8)
     REAL(real_8) :: wquadi(6,6)=HUGE(0.0_real_8)
  END TYPE wannc_t
  TYPE(wannc_t), SAVE :: wannc
  ! ==================================================================

  LOGICAL, SAVE :: hmat_spread=.FALSE.
  REAL(real_8), SAVE :: minspr=HUGE(0.0_real_8)
  ! ==================================================================

  TYPE :: wan05_t
     INTEGER :: loc_npgrp=HUGE(0)
     INTEGER :: loc_relocalize_every=HUGE(0)
     INTEGER :: loc_recompute_dipole_matrices_every=HUGE(0)
     LOGICAL :: loc_relocalize=.FALSE.
     LOGICAL :: loc_relocalize_in_scf=.FALSE.
  END TYPE wan05_t
  TYPE(wan05_t), SAVE :: wan05
  ! ==================================================================
  TYPE :: wanni_t
     INTEGER :: w_maxs=HUGE(0)
     INTEGER :: w_type=HUGE(0)
     INTEGER :: w_opt=HUGE(0)
     INTEGER :: ns_wann=HUGE(0)
     INTEGER :: sw_orb=HUGE(0)
     INTEGER :: sw_first=HUGE(0)
     INTEGER :: sw_last=HUGE(0)
     INTEGER :: sw_all=HUGE(0)
  END TYPE wanni_t
  TYPE(wanni_t), SAVE :: wanni
  TYPE :: wannl_t
     LOGICAL :: twann=.FALSE.
     LOGICAL :: twpri=.FALSE.
     LOGICAL :: tsdens=.FALSE.
     LOGICAL :: twdos=.FALSE.
     LOGICAL :: twmol=.FALSE.
  END TYPE wannl_t
  TYPE(wannl_t), SAVE :: wannl
  TYPE :: wannr_t
     REAL(real_8) :: w_step=HUGE(0.0_real_8)
     REAL(real_8) :: w_eps=HUGE(0.0_real_8)
     REAL(real_8) :: w_ran=HUGE(0.0_real_8)
     REAL(real_8) :: w_ref(3)=HUGE(0.0_real_8)
     REAL(real_8) :: sw_spread=HUGE(0.0_real_8)
  END TYPE wannr_t
  TYPE(wannr_t), SAVE :: wannr

END MODULE wann
