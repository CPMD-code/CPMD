#include "cpmd_global.h"

MODULE cufft_interfaces

  USE cuda_interfaces,                 ONLY: cudaStream_t

  USE, INTRINSIC :: ISO_C_BINDING,  ONLY: C_INT,&
       C_PTR,&
       C_DOUBLE,&
       C_DOUBLE_COMPLEX,&
       C_CHAR, C_NULL_CHAR, C_NULL_PTR, C_BOOL
  USE kinds,                           ONLY: real_8

#include "sizeof.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: cufftPlanMany
  PUBLIC :: cufftDestroy
  PUBLIC :: cufftExecZ2Z
  PUBLIC :: cufftSetStream
  PUBLIC :: cufftGetVersion


  !cufftResult
  PUBLIC :: CUFFT_SUCCESS

  !cufftType
  PUBLIC :: CUFFT_Z2Z


  ENUM, BIND( C ) !:: cufftResult
     ENUMERATOR :: CUFFT_SUCCESS = 0
  END ENUM

  ENUM, BIND( C ) !:: cufftType
     ENUMERATOR :: CUFFT_R2C = 42  !z'2a'    ! Real to complex (interleaved)
     ENUMERATOR :: CUFFT_C2R = 44  !z'2c'    ! Complex (interleaved) to real
     ENUMERATOR :: CUFFT_C2C = 41  !z'29'    ! Complex to complex (interleaved)
     ENUMERATOR :: CUFFT_D2Z = 106 !z'6a'    ! Double to double-complex (interleaved)
     ENUMERATOR :: CUFFT_Z2D = 108 !z'6c'    ! Double-complex (interleaved) to double
     ENUMERATOR :: CUFFT_Z2Z = 105 !z'69'    ! Double-complex to double-complex (interleaved)
  END ENUM

  INTEGER, PUBLIC, PARAMETER :: CUFFT_FORWARD = -1
  INTEGER, PUBLIC, PARAMETER :: CUFFT_INVERSE = 1


  TYPE, PUBLIC, BIND( C ) :: cufftHandle_t
     INTEGER( C_INT ) :: cufftHandle
  END TYPE cufftHandle_t


  INTERFACE

     !cufftResult CUFFTAPI cufftPlanMany(cufftHandle *plan,
     !                                   int rank,
     !                                   int *n,
     !                                   int *inembed, int istride, int idist,
     !                                   int *onembed, int ostride, int odist,
     !                                   cufftType type,
     !                                   int batch);
     INTEGER( KIND( CUFFT_SUCCESS ) ) FUNCTION cufftPlanMany ( plan, rank, n, inembed, istride, idist,  &
          & onembed, ostride, odist, TYPE, batch ) BIND( C, NAME="cufftPlanMany" )
       IMPORT :: C_INT, CUFFT_SUCCESS, CUFFT_R2C, cufftHandle_t, C_PTR
       IMPLICIT NONE
       TYPE( cufftHandle_t ) :: plan
       INTEGER( C_INT ), VALUE :: rank
       INTEGER( C_INT ), DIMENSION(*) :: n
       INTEGER( C_INT ), DIMENSION(*) :: inembed
       INTEGER( C_INT ), VALUE :: istride
       INTEGER( C_INT ), VALUE :: idist
       INTEGER( C_INT ), DIMENSION(*) :: onembed
       INTEGER( C_INT ), VALUE :: ostride
       INTEGER( C_INT ), VALUE :: odist
       INTEGER( KIND( CUFFT_R2C ) ), VALUE :: TYPE
       INTEGER( C_INT ), VALUE :: batch
     END FUNCTION cufftPlanMany

     !cufftDestroy(cufftHandle plan);
     INTEGER( KIND( CUFFT_SUCCESS ) ) FUNCTION cufftDestroy ( plan ) BIND( C, NAME="cufftDestroy" )
       IMPORT :: cufftHandle_t, CUFFT_SUCCESS
       IMPLICIT NONE
       TYPE( cufftHandle_t ), VALUE :: plan
     END FUNCTION cufftDestroy

     !cufftResult cufftExecZ2Z(cufftHandle plan, cufftDoubleComplex *idata, cufftDoubleComplex *odata, int direction);
     INTEGER( KIND( CUFFT_SUCCESS ) ) FUNCTION cufftExecZ2Z ( plan, idata, odata, direction ) BIND( C, NAME="cufftExecZ2Z" )
       IMPORT :: cuffthandle_t, CUFFT_SUCCESS, C_PTR, C_INT
       IMPLICIT NONE
       TYPE( cufftHandle_t ), VALUE :: plan
       TYPE( C_PTR ), VALUE :: idata
       TYPE( C_PTR ), VALUE :: odata
       INTEGER( C_INT ), VALUE :: direction
     END FUNCTION cufftExecZ2Z

     !cufftResult cufftSetStream(cufftHandle plan, cudaStream_t stream);
     INTEGER( KIND( CUFFT_SUCCESS ) ) FUNCTION cufftSetStream ( plan, stream ) BIND( C, NAME='cufftSetStream' )
       IMPORT :: CUFFT_SUCCESS, cufftHandle_t, cudaStream_t
       IMPLICIT NONE
       TYPE( cufftHandle_t ), VALUE :: plan
       TYPE( cudaStream_t ), VALUE :: stream
     END FUNCTION cufftSetStream

     !cufftResult CUFFTAPI cufftGetVersion(int *version);
     INTEGER( KIND( CUFFT_SUCCESS ) ) FUNCTION cufftGetVersion ( version )  BIND( C, NAME='cufftGetVersion' )
       IMPORT :: CUFFT_SUCCESS, C_INT
       IMPLICIT NONE
       INTEGER(C_INT) :: version
     END FUNCTION cufftGetVersion

  END INTERFACE

END MODULE cufft_interfaces
