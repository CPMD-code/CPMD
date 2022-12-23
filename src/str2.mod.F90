MODULE str2
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==================================================================
  ! == INCLUDE FILE for STRESS CALCULATIONS                         ==
  ! ==================================================================
  ! == BECS(1:IMAGP,1:NAT,1:NGH(IS),1:6,1:NDFNL,1:NKPNT)            ==
  ! ==     THE DERIVATIVE OF FNL WITH RESPECT TO H (NLSM1_S)        ==
  ! == BETAA(1:NGWK,1:NGH(IS),1:6,1:NSP,1:NKPNT)                    ==
  ! ==     THE DERIVATIVE OF TWNL (calculated in PUTBET)            ==
  ! ==--------------------------------------------------------------==
  ! == GAGK(NHG,6) GK(ALPHA(KK),IG)*GK(BETA(KK),IG) (KK=1,6)        ==
  ! == GAGKP(NGW,6,NKPNT) (RK+GK)_ALPHA*(RK+GK)_BETA                ==
  ! == GAGKM(NGW,6,NKPNT) (RK-GK)_ALPHA*(RK-GK)_BETA                ==
  ! ==--------------------------------------------------------------==
  COMPLEX(real_8), ALLOCATABLE :: drhovg(:,:)

  ! dimension DQG(NNR1,*)      
  COMPLEX(real_8), ALLOCATABLE, TARGET :: dqg(:,:)

  REAL(real_8), ALLOCATABLE :: becs(:,:,:,:,:,:)
  REAL(real_8), ALLOCATABLE :: betaa(:,:,:,:,:)
  REAL(real_8), ALLOCATABLE :: qrada(:,:,:,:,:)
  REAL(real_8), ALLOCATABLE :: gagk(:,:)
  REAL(real_8), ALLOCATABLE :: sfion(:,:,:)
  REAL(real_8), ALLOCATABLE :: stau0(:,:,:)
  REAL(real_8), ALLOCATABLE :: svelp(:,:,:)
  REAL(real_8), ALLOCATABLE :: gagkp(:,:,:)
  REAL(real_8), ALLOCATABLE :: gagkm(:,:,:)

  ! ==================================================================

END MODULE str2
