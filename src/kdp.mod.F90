MODULE kdp
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==================================================================
  ! ==   DYNAMIC ALLOCATION OF K.P ARRAYS                           ==
  ! ==   WKDP(NKDP)              : WEIGHT OF THE K-POINTS           ==
  ! ==   XKDP(NKDP)              : K-POINTS                         ==
  ! ==   PKDP(3,NSTATE,NSTATE)   : <c0|grad|c0>                     ==
  ! ==   EKDP(NSTATE,NKDP)       : EIGENVALUES                      ==
  ! ==   FKDP(NKDP)              : OCCUPATIONS                      ==
  ! ==   XLAMBDA0(NSTATE,NSTATE) : <c0|g^2|c0>                      ==
  ! ==   CAUX(NGW,NSTATE)        : WFN-LIKE AUXILIARY ARRAY         ==
  ! ==   AKDP(NSTATE,NSTATE,NKDP): TRANSFORMATION MATRICES          ==
  ! ==================================================================
  REAL(real_8), ALLOCATABLE :: wkdp(:)
  REAL(real_8), ALLOCATABLE :: xkdp(:,:)
  REAL(real_8), ALLOCATABLE :: pkdp(:,:,:)
  REAL(real_8), ALLOCATABLE :: ekdp(:,:)
  REAL(real_8), ALLOCATABLE :: fkdp(:,:)
  REAL(real_8), ALLOCATABLE :: rauxdiag(:)

  COMPLEX(real_8), ALLOCATABLE :: xlambda0(:,:)
  COMPLEX(real_8), ALLOCATABLE :: ckdp(:,:)
  COMPLEX(real_8), ALLOCATABLE :: akdp(:,:,:)
  COMPLEX(real_8), ALLOCATABLE :: bkdp(:,:,:)
  COMPLEX(real_8), ALLOCATABLE :: brkdp(:,:)
  COMPLEX(real_8), ALLOCATABLE :: xlambda(:,:)
  COMPLEX(real_8), ALLOCATABLE :: akdp2(:,:)
  COMPLEX(real_8), ALLOCATABLE :: auxdiag(:)

  ! ==================================================================

END MODULE kdp
