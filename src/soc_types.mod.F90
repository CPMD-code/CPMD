MODULE soc_types
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  LOGICAL                   :: md_on_sin, &
       md_on_tri
  LOGICAL                   :: do_soc_because_sh
  INTEGER                   :: nskip,socskip
  REAL(real_8)              :: etot_isc_md
  REAL(real_8), ALLOCATABLE :: tausoc(:,:,:)
  REAL(real_8), ALLOCATABLE :: norm_s(:)
  REAL(real_8), ALLOCATABLE :: norm_t(:)
  REAL(real_8), ALLOCATABLE :: ene_tri(:,:)
  REAL(real_8), ALLOCATABLE :: ene_sin(:,:)
  REAL(real_8), ALLOCATABLE :: soc_mem(:,:)
  COMPLEX(real_8), ALLOCATABLE :: c01(:,:)
  COMPLEX(real_8), ALLOCATABLE :: cs1(:,:,:)
  COMPLEX(real_8), ALLOCATABLE :: ct1(:,:,:)

  REAL(real_8),PARAMETER    ::facthcm=219474.53638012137444
  ! 27.21138386*8065.541154

  TYPE :: socvar_t
     LOGICAL      :: do_sing
     REAL(real_8) :: rr_soc(3)
  END TYPE socvar_t
  TYPE(socvar_t) :: socvar
  ! ==--------------------------------------------------------------==
END MODULE soc_types
