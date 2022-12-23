MODULE cuda_types

  USE cuda_interfaces,                 ONLY: cudaEvent_t,&
                                             cudaStream_t
  USE kinds,                           ONLY: int_8

  USE, INTRINSIC :: ISO_C_BINDING,  ONLY: C_PTR, C_NULL_PTR

  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER, PUBLIC :: cuda_device_null = -HUGE( 0 )

  TYPE, PUBLIC :: cuda_stream_t
     LOGICAL :: init = .FALSE.
     INTEGER :: device = cuda_device_null
     TYPE( cudaStream_t ) :: s
  END TYPE cuda_stream_t

  TYPE, PUBLIC :: cuda_event_t
     LOGICAL :: init = .FALSE.
     TYPE( cudaEvent_t ) :: e
  END TYPE cuda_event_t

  TYPE, PUBLIC :: cuda_memory_t
     LOGICAL :: init = .FALSE.
     INTEGER :: device = cuda_device_null
     INTEGER(int_8) :: n_bytes = 0_int_8
     TYPE( C_PTR ) :: ptr = C_NULL_PTR
  END TYPE cuda_memory_t

END MODULE cuda_types
