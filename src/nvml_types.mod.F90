MODULE nvml_types

  USE nvml_interfaces,                 ONLY: nvmlDevice_t

  IMPLICIT NONE

  PRIVATE

  TYPE, PUBLIC :: nvml_device_t
     LOGICAL :: init = .FALSE.
     TYPE( nvmlDevice_t ) :: d
  END TYPE nvml_device_t

END MODULE nvml_types
