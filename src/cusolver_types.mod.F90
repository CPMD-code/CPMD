MODULE cusolver_types

  USE cuda_types,                      ONLY: cuda_device_null
  USE cusolver_interfaces,             ONLY: cusolverDnHandle_t

  IMPLICIT NONE

  PRIVATE

  TYPE, PUBLIC :: cusolver_handle_t
     LOGICAL :: init = .FALSE.
     INTEGER :: device = cuda_device_null
     TYPE( cusolverDnHandle_t ) :: h
  END TYPE cusolver_handle_t

END MODULE cusolver_types
