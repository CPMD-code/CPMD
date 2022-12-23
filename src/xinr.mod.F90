MODULE xinr
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! =====================================================
  ! == TOL:   TOLERANCE IN IMPLICIT NEWTON-RAPHSON     ==
  ! == NREG = number of gradien regions in INR         ==
  ! == TOLX_INR = convergence thresholds for various   ==
  ! ==            regions                              ==
  ! == GNX_INR = gradient boundaries for the regions   ==
  ! == ITMAX_INR = dimension of INR subspace           ==
  ! == DIRECT = directions in the inner CG-INR cycle   ==
  ! == CORRECT = residuals in the inner CG-INR cycle   ==
  ! == ZPREC = preconditioned directions in CG-INR cyc.==
  ! == HESS_1 = empirical hessian for preconditioning  ==
  ! == INR_STEP -> INR ALTERNATIVE STEP (obsolete?)    ==
  ! == TMIXSD -> mixed procedure, Steepest Descent     ==
  ! == TMIXGDIIS -> mixed procedure, GDIIS             ==
  ! == INR_CONT -> continuation of a broken CG-INR cyc.==
  ! == INR_PREC -> preconditioning of INR              ==
  ! == INR_VERBOSE -> verbosity in INR                 ==
  ! == RMIXSD = gradient for starting INR in mixed opt.==
  ! =====================================================
  INTEGER, PARAMETER :: maxreg=5 
  REAL(real_8) :: tol_inr,tolx_inr(maxreg),gnx_inr(maxreg)

  REAL(real_8), ALLOCATABLE :: direct(:)
  REAL(real_8), ALLOCATABLE :: correct(:)
  REAL(real_8), ALLOCATABLE :: zprec(:)
  REAL(real_8), ALLOCATABLE :: hess_1(:,:)




  REAL(real_8) :: rmixsd

  TYPE :: inr_integer_t
     INTEGER :: itmax_inr
     INTEGER :: nreg
  END TYPE inr_integer_t
  TYPE(inr_integer_t) :: inr_integer
  TYPE :: inr_logical_t
     LOGICAL :: tmixsd
     LOGICAL :: tmixgdiis
     LOGICAL :: inr_step
     LOGICAL :: inr_cont
     LOGICAL :: inr_prec
     LOGICAL :: inr_verbose
  END TYPE inr_logical_t
  TYPE(inr_logical_t) :: inr_logical

END MODULE xinr
