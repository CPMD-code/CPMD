MODULE kpnt
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==================================================================
  ! ==   DYNAMIC ALLOCATION OF K-POINT ARRAYS                       ==
  ! ==   RK(3,NKPTS)           : K-POINTS ARRAYS                    ==
  ! ==   WK(NKPTS)             : WEIGHT                             ==
  ! ==   HGKP(NHG,NKPNT)       : (RK+G)**2                          ==
  ! ==   HGKM(NHG,NKPNT)       : (RK-G)**2                          ==
  ! ==   EIKR(NKPTS,NAT)       : PHASE FACTOR WITH RK     (COMPLEX) ==
  ! ==   EIGKR(NGWK,NAT,NKPNT) : PHASE FACTOR WITH (G+KR) (COMPLEX) ==
  ! ==   KGEMAX(NKPNT)         : MAX GRADIENT PER K POINT 
  ! ==   KCNORM(NKPNT)         : GRADIENT NORM  
  ! ==================================================================

  REAL(real_8), ALLOCATABLE :: rk(:,:) ! RK(3,NKPNT,NBLKP)

  REAL(real_8), ALLOCATABLE :: wk(:)
  REAL(real_8), ALLOCATABLE :: hgkp(:,:)
  REAL(real_8), ALLOCATABLE :: hgkm(:,:)
  REAL(real_8), ALLOCATABLE :: kgemax(:)
  REAL(real_8), ALLOCATABLE :: kcnorm(:)

  REAL(real_8) :: kdrho,kdrhon

  COMPLEX(real_8), ALLOCATABLE :: eikr(:,:)
  COMPLEX(real_8), POINTER :: eigkr(:,:,:)

  ! ==================================================================

END MODULE kpnt
