MODULE cublas_types

  USE cublas_interfaces,               ONLY: cublasHandle_t
  USE cuda_types,                      ONLY: cuda_device_null

  IMPLICIT NONE

  PRIVATE

  TYPE, PUBLIC :: cublas_handle_t
     LOGICAL :: init = .FALSE.
     INTEGER :: device = cuda_device_null
     TYPE( cublasHandle_t ) :: h
  END TYPE cublas_handle_t


END MODULE cublas_types
