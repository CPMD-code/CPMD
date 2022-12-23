MODULE cp_cuvpsi_types

  USE cuda_types,                      ONLY: cuda_memory_t

  IMPLICIT NONE

  PRIVATE


  TYPE, PUBLIC :: cp_cuvpsi_device_t
     LOGICAL :: init = .FALSE.

     !device id
     INTEGER :: id = HUGE(0)

     TYPE( cuda_memory_t ) :: vpotdg_d

  END TYPE cp_cuvpsi_device_t


  TYPE, PUBLIC :: cp_cuvpsi_t
     LOGICAL :: init = .FALSE.

     TYPE( cp_cuvpsi_device_t ), DIMENSION(:), ALLOCATABLE :: devices

  END TYPE cp_cuvpsi_t

END MODULE cp_cuvpsi_types
