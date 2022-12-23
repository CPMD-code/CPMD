MODULE norm
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==================================================================
  ! == INCLUDE FILE FOR CONVERGENCE (ELECTRONIC DENSITY)            ==
  ! ==================================================================
  ! == GEMAX  Max. change in wavefunction components (dE/dC(G))     ==
  ! ==        (or for TDIAG, max. change for density)               ==
  ! == CNORM  Norm of (DE/DC(G)) of for TDIAG norm of DRHO          ==
  ! == GNMAX  Max. force components                                 ==
  ! == GNORM  Norm of Forces                                        ==
  ! ==--------------------------------------------------------------==
  REAL(real_8), SAVE :: gemax=0.0_real_8,&
       cnorm=0.0_real_8,&
       gnmax=0.0_real_8,&
       gnorm=0.0_real_8
  ! ==================================================================

END MODULE norm
