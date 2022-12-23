MODULE tauf
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==--------------------------------------------------------------==
  ! == TAU(NNR1,NLSD)                                               ==
  ! == VTAU(NNR1,NLSD)                                              ==
  ! ==--------------------------------------------------------------==
  INTEGER :: itaur

  REAL(real_8), ALLOCATABLE, TARGET :: tau(:,:)
  REAL(real_8), ALLOCATABLE, TARGET :: vtau(:,:)

  ! ==================================================================

END MODULE tauf
