MODULE lscal
  USE kinds,                           ONLY: real_8
  USE nvarmod,                         ONLY: nvar

  IMPLICIT NONE

  ! ==================================================================
  ! == INPUT PARAMETERS FOR LINEAR SCALING OPTIMIZATION             ==
  ! ==--------------------------------------------------------------==
  ! == NREM     Number of remembered steps for L-BFGS               ==
  ! ==          Default: 0 -> MIN(40,NDIM) at first step            ==
  ! == NTRUST   Type of trust radius / line search                  ==
  ! ==          Default: 0                                          ==
  ! ==          0: HDLCopt trust radius algorithm                   ==
  ! ==          1: Line search                                      ==
  ! ==          2: CPMD optimizer trust radius algorithm            ==
  ! == NRESTT   Periodic restart every NRESTT performed steps       ==
  ! ==          Default: 0 (never)                                  ==
  ! == NTRSTR   Number of performed TR steps before line search     ==
  ! ==          Default: 0 (keep trust radius or line search)       ==
  ! == TRUSTR   Maximum and initial trust radius                    ==
  ! ==          Default: 0.5                                        ==
  ! == See also BLOCK DATA LBFGS_INIT                               ==
  ! ==--------------------------------------------------------------==
  INTEGER :: nrem, ntrust, nrestt, ntrstr
  REAL(real_8) :: step_bmb_l, step_ini_l, step_max_l, step_min_l,&
       trustr
  ! 
  ! ==================================================================
  ! == CONSTANTS FOR LINEAR SCALING OPTIMIZATION                    ==
  ! ==--------------------------------------------------------------==
  ! == BETAP_L, BETA_L: beta' and beta for Wolfe conditions         ==
  ! ==--------------------------------------------------------------==
  REAL(real_8), PARAMETER :: betap_l = 1.0e-4_real_8, beta_l = 0.9_real_8 
  ! 
  ! ==================================================================
  ! == SHARED INFORMATION FOR LINEAR SCALING OPTIMIZATION           ==
  ! ==--------------------------------------------------------------==
  LOGICAL :: lfevmn, linils
  INTEGER :: ibnd_l, igpnt, info_l, ipnt_l, iprnt, ispnt,&
       iter_l, iwpnt, jump, lwlbfgs, nfev_l,nfevmn, nres_l, nstep_l,&
       nstried_l, npt_l
  REAL(real_8) :: deaccp, etotmn, oldene, scl1_l, scl2_l, step_l,&
       sz_l, sz_lmn
  REAL(real_8), ALLOCATABLE :: wlbfgs(:)

  ! ==================================================================
  ! 
  ! ==================================================================
  ! == INPUT PARAMETERS FOR TRANSITION STATE SEARCH                 ==
  ! ==--------------------------------------------------------------==
  ! == MODE       Hessian eigenmode to be followed                  ==
  ! == MODELK     0: Mode following, 1: mode locked                 ==
  ! == NCORE      Number of atoms in the reaction core              ==
  ! == NSMAXP     Maximum number of P-RFO steps                     ==
  ! == NVAR       Number of variables in the core                   ==
  ! == NSVIB      Vibrational analysis is done every NSVIB steps    ==
  ! == ICORE      Atom sequence numbers in the reaction core        ==
  ! == CUTEIG     Cutoff for Hessian eigenvalues                    ==
  ! == DEMIN_P    Minimum energy change to nevertheless accept step ==
  ! == OMIN       Minimum overlap between Hessian modes             ==
  ! == RMAX_P     Maximum ratio actual / predicted energy change    ==
  ! == RMIN_P     Minimum ratio actual / predicted energy change    ==
  ! == STEP_INI_P Initial P-RFO trust radius                        ==
  ! == STEP_MAX_P Maximum P-RFO trust radius                        ==
  ! == STEP_MIN_P Minimum P-RFO trust radius                        ==
  ! == TOLENV     Gradient tolerance for the environment            ==
  ! == TRUSTP     User set P-RFO trust radius                       ==
  ! == EPS_H      Step for finite-difference Hessian                ==
  ! == M_HESS     Type of initial Hessian (0: INIHES, 1: HESSIN)    ==
  ! ==--------------------------------------------------------------==
  INTEGER :: mode, modelk, ncore, nsmaxp, nsvib, m_hess


  INTEGER, ALLOCATABLE :: icore(:)

  REAL(real_8) :: cuteig, demin_p, omin, rmax_p, rmin_p, step_ini_p,&
       step_max_p, step_min_p, tolenv, trustp, eps_h
  ! 
  ! ==================================================================
  ! == CONSTANTS FOR TRANSITION STATE SEARCH                        ==
  ! == FSTEP_P: Step for determining the minimizing lambda in P-RFO ==
  ! ==--------------------------------------------------------------==
  REAL(real_8), PARAMETER :: fstep_p = 5.0e-2_real_8 
  ! 
  ! ==================================================================
  ! == SHARED INFORMATION FOR TRANSITION STATE SEARCH               ==
  ! ==--------------------------------------------------------------==
  LOGICAL :: lallhs, llock, lmicro, lnomap, lprhes, lrstrt,&
       lundop, lvlhes
  INTEGER :: lmap_p, lwprfo, modcur, nshess, nstep_p, nstried_p
  REAL(real_8) :: depred, oprven, ostep_p, osz_p, prvene,step_p,&
       sz_p

  INTEGER, ALLOCATABLE :: map_p(:)

  REAL(real_8), ALLOCATABLE :: wprfo(:)

  ! TODO refactor these arrays 
  REAL(real_8), ALLOCATABLE :: hesscr(:,:)
  REAL(real_8), ALLOCATABLE :: oldhss(:,:)
  REAL(real_8), ALLOCATABLE :: vmode(:)
  REAL(real_8), ALLOCATABLE :: gx(:)
  REAL(real_8), ALLOCATABLE :: oldgx(:)
  REAL(real_8), ALLOCATABLE :: oldg(:)
  REAL(real_8), ALLOCATABLE :: ooldg(:)
  REAL(real_8), ALLOCATABLE :: eigval(:)
  REAL(real_8), ALLOCATABLE :: step(:)
  REAL(real_8), ALLOCATABLE :: oldeig(:)
  REAL(real_8), ALLOCATABLE :: eigvec(:,:)
  REAL(real_8), ALLOCATABLE :: oldev(:,:)


  ! 
  ! ==================================================================
  ! == SHARED INFORMATION FOR STACKS                                ==
  ! ==--------------------------------------------------------------==
  INTEGER :: ielstk, nstack,st_nop, st_pop, st_push, st_get,&
       st_put
  PARAMETER       (ielstk =  1, nstack  = 1,st_nop  =  0,st_pop  = -&
       2, st_push = 2,st_get  = -1, st_put  = 1)
  INTEGER :: iopnxt(nstack), istlen(nstack), iszstk(nstack),&
       itpstk(nstack)
  REAL(real_8), ALLOCATABLE :: elestk(:)

  ! 
  ! ==================================================================
  ! == SHARED INFORMATION FOR ADAPTIVE TOLERANCE                    ==
  ! ==--------------------------------------------------------------==
  LOGICAL :: tnofor
  INTEGER :: nrllst
  REAL(real_8) :: demin, gnmin, tolmin, tolrel

  TYPE :: adtolr_t
     REAL(real_8) :: delast
     REAL(real_8) :: demin
     REAL(real_8) :: gnmin
     REAL(real_8) :: tolmin
     REAL(real_8) :: tolrel
  END TYPE adtolr_t
  TYPE(adtolr_t) :: adtolr

END MODULE lscal
