MODULE func
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  !
  REAL(real_8) :: f0hf,f0hfm,alhfx
  !
  REAL(real_8) :: ashcroft_coulomb_rcut
  !
  !
  ! mfxcx
  INTEGER, PARAMETER :: mfxcx_is_skipped = 0
  INTEGER, PARAMETER :: mfxcx_is_slaterx = 1
  !
  ! mfxcc
  INTEGER, PARAMETER :: mfxcc_is_skipped = 0 
  INTEGER, PARAMETER :: mfxcc_is_pz      = 1
  INTEGER, PARAMETER :: mfxcc_is_vwn     = 2
  INTEGER, PARAMETER :: mfxcc_is_lyp     = 3
  INTEGER, PARAMETER :: mfxcc_is_pw      = 4
  INTEGER, PARAMETER :: mfxcc_is_wigner  = 5
  INTEGER, PARAMETER :: mfxcc_is_hedin   = 6
  INTEGER, PARAMETER :: mfxcc_is_obpz    = 7
  INTEGER, PARAMETER :: mfxcc_is_obpw    = 8
  INTEGER, PARAMETER :: mfxcc_is_pade    = 9
  !
  ! mgcx
  INTEGER, PARAMETER :: mgcx_is_skipped          = 0 
  INTEGER, PARAMETER :: mgcx_is_becke88          = 1
  INTEGER, PARAMETER :: mgcx_is_ggax             = 2
  INTEGER, PARAMETER :: mgcx_is_pbex             = 3
  INTEGER, PARAMETER :: mgcx_is_revpbex          = 4
  INTEGER, PARAMETER :: mgcx_is_hcth             = 5
  INTEGER, PARAMETER :: mgcx_is_optx             = 6
  INTEGER, PARAMETER :: mgcx_is_ox               = 7
  INTEGER, PARAMETER :: mgcx_is_ox_hybrid        = 8
  INTEGER, PARAMETER :: mgcx_is_pbesx            = 9
  INTEGER, PARAMETER :: mgcx_is_dfr_zpbex        = 10
  INTEGER, PARAMETER :: mgcx_is_dfr_xpbex        = 11
  INTEGER, PARAMETER :: mgcx_is_dfr_xpbex_hybrid = 12
  !
  ! mgcc
  INTEGER, PARAMETER :: mgcc_is_skipped          = 0 
  INTEGER, PARAMETER :: mgcc_is_perdew86         = 1
  INTEGER, PARAMETER :: mgcc_is_lyp              = 2
  INTEGER, PARAMETER :: mgcc_is_ggac             = 3
  INTEGER, PARAMETER :: mgcc_is_pbec             = 4
  INTEGER, PARAMETER :: mgcc_is_hse              = 5
  INTEGER, PARAMETER :: mgcc_is_optc             = 6
  INTEGER, PARAMETER :: mgcc_is_pbesc            = 7
  INTEGER, PARAMETER :: mgcc_is_dfr_zpbec        = 8
  !
  ! mgcc
  INTEGER, PARAMETER :: mtau_is_skipped          = 0
  INTEGER, PARAMETER :: mtau_is_tpss             = 1
  !
  !
  INTEGER, PARAMETER :: mhfx_is_skipped          = 0
  INTEGER, PARAMETER :: mhfx_is_hfx              = 1
  INTEGER, PARAMETER :: mhfx_is_hartree          = 2
  !
  ! mgcsrx
  INTEGER, PARAMETER :: msrx_is_skipped          = 0
  INTEGER, PARAMETER :: msrx_is_exp              = 1
  INTEGER, PARAMETER :: msrx_is_erfc             = 2
  INTEGER, PARAMETER :: msrx_is_ashcroft         = 3
  INTEGER, PARAMETER :: msrx_is_CAM              = 4
  !
  ! mgcsrx
  INTEGER, PARAMETER :: mgcsrx_is_skipped        = 0
  INTEGER, PARAMETER :: mgcsrx_is_hse            = 1
  !
  !
  TYPE :: func1_t
     INTEGER :: mfxcx
     INTEGER :: mfxcc
     INTEGER :: mgcx
     INTEGER :: mgcc
     INTEGER :: mhfx
     INTEGER :: mtau
     INTEGER :: msrx
     INTEGER :: mgcsrx
  END TYPE func1_t
  TYPE(func1_t) :: func1
  TYPE :: func2_t
     REAL(real_8) :: salpha
     REAL(real_8) :: bbeta
     REAL(real_8) :: betapp
     REAL(real_8) :: srxa
     REAL(real_8) :: cam_alpha
     REAL(real_8) :: cam_beta
  END TYPE func2_t
  TYPE(func2_t) :: func2
  TYPE :: func3_t
     REAL(real_8) :: pxlda
     REAL(real_8) :: pxgc
     REAL(real_8) :: pclda
     REAL(real_8) :: pcgc
     REAL(real_8) :: phfx
  END TYPE func3_t
  TYPE(func3_t) :: func3

END MODULE func
