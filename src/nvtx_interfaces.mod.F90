MODULE nvtx_interfaces

  USE, INTRINSIC :: ISO_C_BINDING,  ONLY: C_INT, C_SIZE_T, &
       C_PTR,&
       C_DOUBLE, C_FLOAT, &
       C_DOUBLE_COMPLEX,&
       C_CHAR, C_NULL_CHAR, C_NULL_PTR

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: nvtxRangePushA
  PUBLIC :: nvtxRangePop

  INTERFACE

     SUBROUTINE nvtxRangePushA( name ) BIND( C, NAME = 'nvtxRangePushA' )
       IMPORT :: C_CHAR
       CHARACTER(C_CHAR) :: name(*)
     END SUBROUTINE nvtxRangePushA

     SUBROUTINE nvtxRangePop( ) BIND( C, NAME = 'nvtxRangePop' )
     END SUBROUTINE nvtxRangePop

  END INTERFACE

END MODULE nvtx_interfaces
