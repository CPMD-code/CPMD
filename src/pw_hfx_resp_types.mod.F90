MODULE pw_hfx_resp_types

  USE kinds,                           ONLY: real_8

  IMPLICIT NONE
  PRIVATE 

  TYPE , PUBLIC :: pw_hfx_resp_t
     LOGICAL :: init = .FALSE.
     ! use symmetry to save space (ib>=ia)
     ! check for unrestricted case (lsd)
     REAL(real_8), DIMENSION(:,:,:), ALLOCATABLE :: v_ab
     LOGICAL :: is_set = .FALSE.
  END TYPE pw_hfx_resp_t
  TYPE(pw_hfx_resp_t), PUBLIC, SAVE :: hfx_resp

  TYPE , PUBLIC :: pw_hfx_lin2_t
     LOGICAL :: init = .FALSE.
     COMPLEX(real_8), ALLOCATABLE, &
          DIMENSION(:, :)                        :: c2, &
          xi
     REAL(real_8), ALLOCATABLE, &
          DIMENSION(:, :)                        :: m, l
  END TYPE pw_hfx_lin2_t
  TYPE(pw_hfx_lin2_t), PUBLIC, SAVE :: hfx_lin2


  TYPE , PUBLIC :: pw_hfx_resp_env_t
     LOGICAL :: init = .FALSE.
     LOGICAL :: store_v_ab = .FALSE.
     LOGICAL :: write_v_ab = .FALSE.
     LOGICAL :: use_lin_lin = .FALSE.
     LOGICAL :: lin2_restart = .FALSE.
     LOGICAL :: switch_lin2 = .FALSE.
     REAL(real_8) :: sigma_eps = 0.00001_real_8
     INTEGER :: interval_restart_lin2 = 4
  END TYPE pw_hfx_resp_env_t
  TYPE(pw_hfx_resp_env_t), PUBLIC, SAVE :: hfx_resp_env

END MODULE pw_hfx_resp_types
