MODULE cufft_types

  USE cuda_types,                      ONLY: cuda_device_null
  USE cufft_interfaces,                ONLY: cufftHandle_t

  IMPLICIT NONE

  PRIVATE

  TYPE, PUBLIC :: cufft_plan_t
     LOGICAL :: init = .FALSE.
     INTEGER :: device = cuda_device_null
     TYPE( cufftHandle_t ) :: h
  END TYPE cufft_plan_t

END MODULE cufft_types
