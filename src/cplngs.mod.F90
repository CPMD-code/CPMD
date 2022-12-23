MODULE cplngsmod
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==================================================================
  ! == INPUT PARAMETERS FOR NONADIABATIC COUPLINGS                  ==
  ! ==--------------------------------------------------------------==
  ! == TCPL     Calculate couplings between surfaces                ==
  ! == TCPLFD   Use finite difference for couplings                 ==
  ! == TALLAT   Do FD couplings with all atoms                      ==
  ! == TCPLLR   Use linear response for analytic couplings          ==
  ! == TALLDOF  Brute force instead of iterative response density   ==
  ! == TSPECV   Specify vectors                                     ==
  ! == NSURF    Number of KS orbital pairs for couplings            ==
  ! == ISURF    KS orbital numbers (below)                          ==
  ! == NATFD    Number of atoms involved in FD couplings            ==
  ! == NVECT    Maximum number of vectors for iterative LR scheme   ==
  ! == IATFD    Atoms involved in FD couplings                      ==
  ! == EPS_C    Finite difference displacement                      ==
  ! == TOLCPL   Tolerance for iterative LR scheme of couplings      ==
  ! == FC_LOW   Threshold increased LR accuracy (FC_MED, FC_HIGH)   ==
  ! == F_LOW    Factor for increased LR accuracy (F_MED, F_HIGH)    ==
  ! ==--------------------------------------------------------------==
  LOGICAL :: tcpl, tcplfd, tallat, tcpllr, talldof,tspecv
  INTEGER :: nsurf, natfd, nvect
  REAL(real_8) :: eps_c, tolcpl,fc_low, fc_med, fc_high, f_low,&
       f_med, f_high
  ! ==================================================================
  ! == DYNAMIC ALLOCATION OF ARRAYS RELATED TO COUPLINGS            ==
  ! ==--------------------------------------------------------------==
  ! == ISURF    (see above, also CSURF)                             ==
  ! == IATFD    (see above)                                         ==
  ! == CPLION   Non-adiabatic coupling vector between KS states     ==
  ! == CPLCFG   Non-adiabatic coupling vector between configs       ==
  ! ==--------------------------------------------------------------==
  INTEGER, ALLOCATABLE :: isurf(:,:)
  INTEGER, ALLOCATABLE :: iatfd(:)

  REAL(real_8), ALLOCATABLE :: csurf(:)

  REAL(real_8), ALLOCATABLE :: cplion(:,:,:,:)
  REAL(real_8), ALLOCATABLE :: cplcfg(:,:,:)


  ! ==================================================================

END MODULE cplngsmod
