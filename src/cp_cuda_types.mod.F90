MODULE cp_cuda_types

  IMPLICIT NONE

  PRIVATE

  TYPE, PUBLIC :: cp_cuda_env_t
     LOGICAL :: use_fft = .FALSE.
     LOGICAL :: use_blas = .FALSE.
     INTEGER :: fft_n_streams_per_device = 1
     INTEGER :: fft_n_devices_per_task = HUGE(0)
     INTEGER :: blas_n_streams_per_device = 1
     INTEGER :: blas_n_devices_per_task = HUGE(0)
  END TYPE cp_cuda_env_t

  TYPE(cp_cuda_env_t), SAVE, PUBLIC :: cp_cuda_env

  TYPE, PUBLIC :: cp_cuda_devices_t

     !list of DEVICE id's used by the current task
     ! bounds (1:n_devices_per_task)
     INTEGER, DIMENSION( : ), ALLOCATABLE :: ids

  END TYPE cp_cuda_devices_t

  TYPE(cp_cuda_devices_t), TARGET, SAVE, PUBLIC :: cp_cuda_devices_fft
  TYPE(cp_cuda_devices_t), TARGET, SAVE, PUBLIC :: cp_cuda_devices_blas

END MODULE cp_cuda_types
