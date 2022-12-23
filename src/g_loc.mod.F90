MODULE g_loc
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==================================================================
  ! == INCLUDE FILE FOR FUNCTIONS LOCALIZED IN RECIPROCAL SPACE     ==
  ! ==--------------------------------------------------------------==
  ! == TGLOC     : calculate  orbitals localized in G space         ==
  ! == TGLOCPRINT  : print orbitals loc. in G space                 ==
  ! == TG_KICK   : starting random unitary transformation on C0     ==
  ! == TG_COMPLEX: starting C0 already complex in real space        ==
  ! == TG_REAL   : starting C0 real in real space, to be expanded   ==
  ! == TG_TUNE_ST: tuning the step on the size of the gradient      ==
  ! == TG_ANTISYMM: wfn forced  not to be symmetric wrt  G = 0      ==
  ! == TG_LINESEARCH: tuning of the time step by 3 points method    ==
  ! == TG_READ_MATRIX: read unitary  transformation in file MATRIX  ==
  ! == TGLOCREALP: print orbitals in real space                     ==
  ! == TGWANNIER:  localization in real space                       ==
  ! ==--------------------------------------------------------------==
  INTEGER, PARAMETER :: lgloclog = 13 

  ! ==--------------------------------------------------------------==
  ! == GLOC_MAXS : Maximum number of optimization steps             ==
  ! == GLOC_TYPE : Type of optimization functional                  ==
  ! ==          1   SUM_(n,i) -<n|Gi|n>^2                           == 
  ! ==          2   SUM_(n,i) 1-|<n|exp{iGi*r_0}|n>|^2              ==
  ! == GLOC_OPT  : Optimisation procedure                           ==
  ! ==          1 U matrix integration C0 update at the very end    ==
  ! ==          2 U kept as identity + D_U, C0 update each step     ==
  ! == GLOC_CONST: How to impose the unitaruty constraint           ==
  ! ==          1 Lagrangian multipliers                            ==
  ! ==          2 Exact orthogonalization by S^(-1/2)               ==
  ! ==          3 Approximation U = 1+iA, + orthogonalization of C0 ==
  ! == GLOC_FIRST : first orbital to be printed                     ==
  ! == GLOC_LAST  : last  orbital to be printed                     ==
  ! == GLOC_ALL   : all   orbital to be printed                     ==
  ! == GLOC_ORB : number of orbitals to be printed                  ==
  ! == GLOC_LIST: list of orbitals to be printed                    ==
  ! == GLOC_MAXIT: Maximum number of U constraint iteration         ==
  ! == GLOC_INIT : First step -1, in case of restart                ==
  ! ==--------------------------------------------------------------==
  INTEGER, ALLOCATABLE :: gloc_list(:)

  ! ==--------------------------------------------------------------==
  ! == GLOC_STEP : step size for optimisation                       ==
  ! == GLOC_EPS  : convergence criteria                             ==
  ! == GLOC_RAN  : randomization amplitude                          ==
  ! == GEPSLAG   : convergence criteria for Lagrangian multipliers  ==
  ! == WAN_WEIGHT : weight of R spread functional (wannier) (def=0) ==
  ! == G2-G_WEIGHT: weight of G spread functional      (def=1)      ==
  ! ==--------------------------------------------------------------==
  INTEGER, PARAMETER :: lglocpar = 6 
  REAL(real_8) :: gepslag
  ! ==--------------------------------------------------------------==
  COMPLEX(real_8), ALLOCATABLE :: lag_mul(:)

  REAL(real_8), ALLOCATABLE :: omega_n(:)
  REAL(real_8) ::shiftkx,shiftky,shiftkz,wan_mem,g2g_mem
  INTEGER :: nstep,gifirst
  INTEGER, PARAMETER :: lglocre = 11 


  INTEGER, ALLOCATABLE :: ist_list(:)
  INTEGER, ALLOCATABLE :: ind_st(:)


  ! ==================================================================
  CHARACTER(len=20) :: filspr,filcen,filpen,filmat

  ! ==================================================================

  TYPE :: glocal_t
     LOGICAL :: tgloc
     LOGICAL :: tglocprint
     LOGICAL :: tg_kick
     LOGICAL :: tg_complex
     LOGICAL :: tg_real
     LOGICAL :: tg_tune_st
     LOGICAL :: tg_antisymm
     LOGICAL :: tg_penalty
     LOGICAL :: tg_linesearch
     LOGICAL :: tg_read_matrix
     LOGICAL :: tglocrealp
     LOGICAL :: tgwannier
     LOGICAL :: torbdist
  END TYPE glocal_t
  TYPE(glocal_t) :: glocal
  TYPE :: gloci_t
     INTEGER :: gloc_maxs
     INTEGER :: gloc_type
     INTEGER :: gloc_opt
     INTEGER :: gloc_const
     INTEGER :: gloc_orb
     INTEGER :: gloc_first
     INTEGER :: gloc_last
     INTEGER :: gloc_all
     INTEGER :: gloc_maxit
     INTEGER :: gloc_init
  END TYPE gloci_t
  TYPE(gloci_t) :: gloci
  TYPE :: gloc_re_t
     REAL(real_8) :: ofun
     REAL(real_8) :: omega_tot
     REAL(real_8) :: gmax
     REAL(real_8) :: ofun0
     REAL(real_8) :: ggnorm
     REAL(real_8) :: dif_fun
     REAL(real_8) :: xyzfun
     REAL(real_8) :: mass_r
     REAL(real_8) :: mass_g
  END TYPE gloc_re_t
  TYPE(gloc_re_t) :: gloc_re
  TYPE :: glocr_t
     REAL(real_8) :: gloc_step
     REAL(real_8) :: gloc_eps
     REAL(real_8) :: gloc_ran
     REAL(real_8) :: gepslag
     REAL(real_8) :: wan_weight
     REAL(real_8) :: g2g_weight
  END TYPE glocr_t
  TYPE(glocr_t) :: glocr
  TYPE :: indstate_t
     INTEGER :: ist_first
     INTEGER :: ist_last
     INTEGER :: nst_list
  END TYPE indstate_t
  TYPE(indstate_t) :: indstate
  TYPE :: lostate_t
     LOGICAL :: state_all
     LOGICAL :: state_range
     LOGICAL :: state_list
  END TYPE lostate_t
  TYPE(lostate_t) :: lostate

END MODULE g_loc
