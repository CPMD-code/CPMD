MODULE hfxmod
  USE kinds,                           ONLY: real_8
! ==================================================================

  IMPLICIT NONE
  INTEGER, SAVE :: ipoolhfx = HUGE(0)
  ! ==================================================================

  TYPE :: hfxc2_t
     REAL(real_8) :: hfxwfe = HUGE(0.0_real_8)
     REAL(real_8) :: hfxdee = HUGE(0.0_real_8)
  END TYPE hfxc2_t
  TYPE(hfxc2_t), SAVE :: hfxc2
  ! ==================================================================

  TYPE :: hfxc3_t
     LOGICAL :: twscr=.FALSE.
     LOGICAL :: twfc=.FALSE.
     LOGICAL :: twft=.FALSE.
     LOGICAL :: tdiaw=.FALSE.
     LOGICAL :: keep_radius_fixed=.FALSE.
     LOGICAL :: use_new_hfx=.FALSE.
  END TYPE hfxc3_t
  TYPE(hfxc3_t), SAVE :: hfxc3

  TYPE :: hfxc4_t
     REAL(real_8) :: dwfc = HUGE(0.0_real_8)
     REAL(real_8) :: dwf_integral_thresh = HUGE(0.0_real_8)
  END TYPE hfxc4_t
  TYPE(hfxc4_t), SAVE :: hfxc4

  TYPE :: hfxc5_t
     INTEGER :: hfx_block_size = HUGE(0)
     INTEGER :: hfx_distribution = HUGE(0)
     INTEGER :: recomp_two_int_list_every = HUGE(0)
  END TYPE hfxc5_t
  TYPE(hfxc5_t), SAVE :: hfxc5
  ! ==================================================================
  REAL(real_8), ALLOCATABLE, SAVE :: wcentx(:,:) ! (4,*)
  ! ==================================================================

END MODULE hfxmod
