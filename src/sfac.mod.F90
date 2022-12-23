MODULE sfac
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==--------------------------------------------------------------==
  ! == FNL(IMAGP,NAT,NHXS,NSTATE,NKPNT)                             ==
  ! == DFNL(IMAGP,NAT,NHXS,3,NDFNL,NKPNT)                           ==
  ! ==--------------------------------------------------------------==


  REAL(real_8), POINTER :: fnl(:,:,:,:,:)
  REAL(real_8), POINTER :: fnl2(:,:,:,:,:)

  REAL(real_8), POINTER :: dfnl(:,:,:,:,:,:)
  REAL(real_8), ALLOCATABLE :: ddfnl(:,:,:,:)

  ! NOTE: not clear what happens with this var... 
  REAL(real_8), ALLOCATABLE :: FNLGP(:,:,:,:,:) ! ! FNLGP(IMAGP,NAT,NHXS,LDF2,NKPOINT)

  ! (:,:,1) because EIGR is used only if NKPNT=1. in K-points EIGKR is used instead 
  COMPLEX(real_8), ALLOCATABLE, TARGET :: eigr(:,:,:)

  COMPLEX(real_8), ALLOCATABLE :: eigrb(:,:)

  LOGICAL :: tfnl2
  INTEGER :: ldf1,ldf2
  ! ==================================================================
  INTEGER :: natx
  COMPLEX(real_8), ALLOCATABLE :: ei1(:,:)
  COMPLEX(real_8), ALLOCATABLE :: ei2(:,:)
  COMPLEX(real_8), ALLOCATABLE :: ei3(:,:)


  ! ==================================================================

END MODULE sfac
