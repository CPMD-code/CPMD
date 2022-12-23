MODULE locpot
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  REAL(real_8), ALLOCATABLE, TARGET :: vpotx3a(:,:,:)

  REAL(real_8), ALLOCATABLE, TARGET :: vpotx3b(:,:,:)

  LOGICAL :: lg_vpotx3a,lg_vpotx3b
  INTEGER :: iclpot,nrxyz2


END MODULE locpot
