MODULE array_utils

  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: int_4,&
                                             real_8

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: array_alloc
  PUBLIC :: array_dealloc
  PUBLIC :: array_realloc
  PUBLIC :: array_ensure_size

  INTERFACE array_alloc
     MODULE PROCEDURE i4_array_alloc_r1
     MODULE PROCEDURE i4_array_alloc_r2
     MODULE PROCEDURE l_array_alloc_r1
     MODULE PROCEDURE l_array_alloc_r2
  END INTERFACE array_alloc

  INTERFACE array_dealloc
     MODULE PROCEDURE i4_array_dealloc_r1
     MODULE PROCEDURE i4_array_dealloc_r2
     MODULE PROCEDURE l_array_dealloc_r1
     MODULE PROCEDURE l_array_dealloc_r2
  END INTERFACE array_dealloc

  INTERFACE array_realloc
     MODULE PROCEDURE i4_array_realloc_r1
     MODULE PROCEDURE i4_array_realloc_r2
     MODULE PROCEDURE l_array_realloc_r1
     MODULE PROCEDURE l_array_realloc_r2
  END INTERFACE array_realloc

  INTERFACE array_ensure_size
     MODULE PROCEDURE i4_array_ensure_size_r1
     MODULE PROCEDURE i4_array_ensure_size_r2
     MODULE PROCEDURE l_array_ensure_size_r1
     MODULE PROCEDURE l_array_ensure_size_r2
  END INTERFACE array_ensure_size

  REAL(real_8), PARAMETER, PRIVATE :: array_default_growing_factor = 1.3_real_8

CONTAINS

  !
  ! include file for the interfaces
  !
#include "array_utils.inc"


END MODULE array_utils
