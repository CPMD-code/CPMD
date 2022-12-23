MODULE nabdy_types
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  INTEGER, PARAMETER ::  NEXP=4           ! max no. of spline exponent:0,1,2,...NEXP-1
  INTEGER, PARAMETER ::  NGRIDP=5         ! grid points for the atomic integrals
  !  ==================================================================
  !  == NABDY GLOBAL VARIABLES   nabdy.inc                           ==
  !  ==================================================================
  TYPE nabdyvec_t
     REAL(real_8), ALLOCATABLE :: NACOOR(:,:,:,:)  !(3,NAX,NSX,NTRAJBD) coordinates of the NABDY fluid elements
     REAL(real_8), ALLOCATABLE :: NACOORP(:,:,:,:) !(3,NAX,NSX,NTRAJBD) ! previous coordinates 
     REAL(real_8), ALLOCATABLE :: NAMOM(:,:,:,:)   !(3,NAX,NSX,NTRAJBD)  ! momenta of the NABDY fluid elements
     REAL(real_8), ALLOCATABLE :: NAFOR(:,:,:,:)   !(3,NAX,NSX,NTRAJBD)  ! forces acting on the NABDY fluid elements
     REAL(real_8), ALLOCATABLE :: NAAMPL(:,:,:)    !(NAX,NSX,NTRAJBD)   ! aplitude of each gaussian
     REAL(real_8), ALLOCATABLE :: NAOMEGA(:,:,:)     !(NAX,NSX,NTRAJBD)  ! width of each gaussian
     REAL(real_8), ALLOCATABLE :: NAAMPL_TEMP(:,:,:)!(NAX,NSX,NTRAJBD)   !
     REAL(real_8), ALLOCATABLE :: NAACT(:,:,:)     !(NAX,NSX,NTRAJBD)    ! action at the atomic positions
     REAL(real_8), ALLOCATABLE :: DDNAACT(:,:,:)     !(NAX,NSX,NTRAJBD)
     REAL(real_8), ALLOCATABLE :: NANORMALF(:,:,:)     !(NAX,NSX,NTRAJBD)! normalization factor 
     REAL(real_8), ALLOCATABLE :: GRIDSHEP(:,:)      !(3,*)             ! Grid for Shepard integral
     REAL(real_8), ALLOCATABLE :: DDELTA(:)
     REAL(real_8), ALLOCATABLE :: NAQP_pat(:,:)    !(NAX,NSX)     !
     REAL(real_8), ALLOCATABLE :: NATOTQP(:),NATOTCL(:),NATOTK(:),NAAMPFE(:) !(NTRAJBD)

     REAL(real_8), ALLOCATABLE :: past_namom(:,:,:,:)
     REAL(real_8), ALLOCATABLE :: past_naact(:,:,:)
     REAL(real_8), ALLOCATABLE :: past_dds(:,:,:)
  END TYPE nabdyvec_t
  TYPE(nabdyvec_t) :: nabdyvec

  TYPE nabdyvar_t
     REAL(real_8)  :: NAINTDV,NABDY_TIME,NAFRITC,NASOFTENING
     INTEGER       :: NABDY_ZMAX = 1, &
                      NTRAJBD = 1     
     LOGICAL       :: TNAREPART        = .FALSE., &
                      SCALEP           = .FALSE., &
                      TNAFRIC          = .FALSE., &
                      NAFORCE_SCREENED = .FALSE.
  END TYPE nabdyvar_t
  TYPE(nabdyvar_t) :: nabdyvar
  !     ==--------------------------------------------------------------==
  !     for the spline algorithm
  !     ==--------------------------------------------------------------==
  TYPE nasplin_t
     REAL(real_8) :: GRIDCOOR(3,NGRIDP)      ! coordinated of the atomic grid points 
     REAL(real_8) :: DELTA                   ! grid spacing
     REAL(real_8) :: WSPL(3,NGRIDP)          ! weight functions in eq. 3
     REAL(real_8) :: SCOEFF(NGRIDP,NGRIDP,NGRIDP,NEXP,NEXP,NEXP) ! coeff. of the spline
     REAL(real_8) :: FS(NGRIDP,NGRIDP,NGRIDP)       ! eigenvectors f_{n1,n2} in eq 4 (and 7)
     REAL(real_8) :: X(NEXP*(NGRIDP-1),NEXP*(NGRIDP-1)),&
          XINV(NEXP*(NGRIDP-1),NEXP*(NGRIDP-1)),&
          Y(NEXP*(NGRIDP-1),NEXP*(NGRIDP-1)),&
          YINV(NEXP*(NGRIDP-1),NEXP*(NGRIDP-1)),&
          Z(NEXP*(NGRIDP-1),NEXP*(NGRIDP-1)),&
          ZINV(NEXP*(NGRIDP-1),NEXP*(NGRIDP-1))
  END TYPE nasplin_t
  TYPE (nasplin_t) :: nasplin
  !     ==--------------------------------------------------------------==
  !     for the friction algorithm
  !     ==--------------------------------------------------------------==
  TYPE nabdyfric_t
     REAL(real_8), ALLOCATABLE :: NABDY_MEAN(:)
     REAL(real_8) :: NATEMPCM,NATEMPBP
     REAL(real_8), ALLOCATABLE, DIMENSION(:,:) :: NABDY_DX,NABDY_DX2,NABDY_DP,NABDY_DP2
     REAL(real_8), ALLOCATABLE :: dispx_ref(:,:),dispp_ref(:,:)
  END TYPE nabdyfric_t
  TYPE (nabdyfric_t) :: nabdyfric

END MODULE nabdy_types
