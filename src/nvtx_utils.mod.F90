#include "cpmd_global.h"

MODULE nvtx_utils

  USE, INTRINSIC :: ISO_C_BINDING,  ONLY: C_INT, C_SIZE_T, &
       C_PTR,&
       C_DOUBLE, C_FLOAT, &
       C_DOUBLE_COMPLEX,&
       C_CHAR, C_NULL_CHAR, C_NULL_PTR
  USE nvtx_interfaces, ONLY: nvtxRangePushA, nvtxRangePop


  IMPLICIT NONE

  PRIVATE

  PUBLIC :: nvtx_range_push_a
  PUBLIC :: nvtx_range_pop

CONTAINS

  SUBROUTINE nvtx_range_push_a ( name )
    CHARACTER(*), INTENT(IN)                 :: name

#if defined(_HAS_CUDA)

    CALL nvtxRangePushA ( name//C_NULL_CHAR )

#endif

  END SUBROUTINE nvtx_range_push_a

  SUBROUTINE nvtx_range_pop ( )

#if defined(_HAS_CUDA)

    CALL nvtxRangePop ( )

#endif

  END SUBROUTINE nvtx_range_pop

END MODULE nvtx_utils
