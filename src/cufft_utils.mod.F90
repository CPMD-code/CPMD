#include "cpmd_global.h"

MODULE cufft_utils

  USE cuda_types,                      ONLY: cuda_device_null,&
                                             cuda_memory_t,&
                                             cuda_stream_t
  USE cuda_utils,                      ONLY: cuda_check_device,&
                                             cuda_set_device
  USE cuda_interfaces,                 ONLY: cudaMemGetInfo
  USE cufft_interfaces,                ONLY: CUFFT_SUCCESS,&
                                             cufftDestroy,&
                                             cufftExecZ2Z,&
                                             cufftGetVersion,&
                                             cufftPlanMany,&
                                             cufftSetStream
  USE cufft_types,                     ONLY: cufft_plan_t
  USE error_handling,                  ONLY: stopgm
  USE string_utils,                    ONLY: int2str

  USE, INTRINSIC :: ISO_C_BINDING,  ONLY: C_INT,&
       C_PTR,&
       C_DOUBLE,&
       C_DOUBLE_COMPLEX,&
       C_CHAR, C_NULL_CHAR, C_NULL_PTR, C_BOOL, C_SIZE_T
  USE kinds,                           ONLY: real_8, int_8


  IMPLICIT NONE

  PRIVATE

  PUBLIC :: cufft_plan_create
  PUBLIC :: cufft_plan_destroy
  PUBLIC :: cufft_execz2z
  PUBLIC :: cufft_set_stream
  PUBLIC :: cufft_get_version

CONTAINS


  SUBROUTINE cufft_plan_create ( plan, rank, n, inembed, istride, idist, onembed, ostride, odist, TYPE, batch, device )
    TYPE(cufft_plan_t)                       :: plan
    INTEGER, INTENT(IN)                      :: rank, n, inembed, istride, &
                                                idist, onembed, ostride, &
                                                odist, TYPE, batch, device

    CHARACTER(*), PARAMETER :: procedureN = 'cufft_plan_create'

    INTEGER(C_INT) :: c_batch, c_idist, c_inembed, c_istride, c_n, c_odist, &
      c_onembed, c_ostride, c_rank, c_status, c_type, c_status_info
    INTEGER(C_SIZE_T) :: c_free, c_total

#if defined(_HAS_CUDA)

    CALL cuda_check_device ( device, procedureN )

    IF( plan%init ) CALL stopgm(procedureN,'plan already created',&
         __LINE__,__FILE__)

    c_rank = INT( rank, C_INT )
    c_n = INT( n, C_INT )
    c_inembed = INT( inembed, C_INT )
    c_istride = INT( istride, C_INT )
    c_idist = INT( idist, C_INT )
    c_onembed = INT( onembed, C_INT )
    c_ostride = INT( ostride, C_INT )
    c_odist = INT( odist, C_INT )
    c_type = INT( TYPE, C_INT )
    c_batch = INT( batch, C_INT )

    CALL cuda_set_device ( device )
    c_status = cufftPlanMany ( plan%h, c_rank, [c_n], [c_inembed], c_istride, c_idist, [c_onembed], c_ostride, &
         & c_odist, c_type, c_batch )

    IF( c_status /= cufft_success ) THEN
       c_status_info = cudaMemGetInfo ( c_free, c_total )
       CALL stopgm(procedureN,"cufft error: "//TRIM(int2str( INT( c_status ) ))//'. '//&
            'Memory free='//TRIM(int2str(INT(c_free,int_8)))//', total='//TRIM(int2str(INT(c_total,int_8))),&
            __LINE__,__FILE__)
    ENDIF

    plan%device = device
    plan%init = .TRUE.

#else

    CALL stopgm(procedureN,"no cuda available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE cufft_plan_create


  SUBROUTINE cufft_plan_destroy ( plan )
    TYPE(cufft_plan_t)                       :: plan

    CHARACTER(*), PARAMETER :: procedureN = 'cufft_plan_destroy'

    INTEGER(C_INT)                           :: c_status

#if defined(_HAS_CUDA)

    IF( .NOT.plan%init ) CALL stopgm(procedureN,'plan not created',&
         __LINE__,__FILE__)

    CALL cuda_set_device ( plan%device )
    c_status = cufftDestroy ( plan%h )
    IF( c_status /= cufft_success ) CALL stopgm(procedureN,"cufft error: "//TRIM(int2str( INT( c_status ) )),&
         __LINE__,__FILE__)

    plan%device = cuda_device_null
    plan%init = .FALSE.

#else

    CALL stopgm(procedureN,"no cuda available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE cufft_plan_destroy


  SUBROUTINE cufft_set_stream ( plan, stream )
    TYPE(cufft_plan_t)                       :: plan
    TYPE(cuda_stream_t)                      :: stream

    CHARACTER(*), PARAMETER :: procedureN = 'cufft_set_stream'

    INTEGER(C_INT)                           :: c_status

#if defined(_HAS_CUDA)

    IF( .NOT.plan%init ) CALL stopgm(procedureN,'plan not created',&
         __LINE__,__FILE__)

    IF( .NOT.stream%init ) CALL stopgm(procedureN,'stream not created',&
         __LINE__,__FILE__)

    IF( plan%device /= stream%device ) CALL stopgm(procedureN,'stream and plan dont share the same device',&
         __LINE__,__FILE__)

    CALL cuda_check_device ( plan%device, procedureN )

    c_status = cufftSetStream ( plan%h, stream%s )
    IF( c_status /= cufft_success ) CALL stopgm(procedureN,"cufft error: "//TRIM(int2str( INT( c_status ) )),&
         __LINE__,__FILE__)

#else

    CALL stopgm(procedureN,"no cuda available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE cufft_set_stream


  SUBROUTINE cufft_get_version ( version )
    INTEGER, INTENT(OUT)                     :: version

    CHARACTER(*), PARAMETER :: procedureN = 'cufft_get_version'

    INTEGER(C_INT)                           :: c_status, c_version

    version = 0

#if defined(_HAS_CUDA)

    c_status = cufftGetVersion ( c_version )
    IF( c_status /= cufft_success ) CALL stopgm(procedureN,"cufft error: "//TRIM(int2str( INT( c_status ) )),&
         __LINE__,__FILE__)

    version = INT( c_version )

#else

    CALL stopgm(procedureN,"no cuda available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE cufft_get_version

  SUBROUTINE cufft_execz2z ( plan, i_devc, o_devc, direction )
    TYPE(cufft_plan_t), INTENT(IN)           :: plan
    TYPE(cuda_memory_t)                      :: i_devc, o_devc
    INTEGER, INTENT(IN)                      :: direction

    CHARACTER(*), PARAMETER                  :: procedureN = 'cufft_execz2z'

    INTEGER(C_INT)                           :: c_direction, c_status

#if defined(_HAS_CUDA)

    IF( .NOT.plan%init ) CALL stopgm(procedureN,'plan not created',&
         __LINE__,__FILE__)
    IF( .NOT.i_devc%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)
    IF( .NOT.o_devc%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)
    IF( plan%device /= i_devc%device .OR. plan%device /= o_devc%device ) &
         CALL stopgm(procedureN,'plan and memory dont share the same device',&
         __LINE__,__FILE__)
    CALL cuda_check_device ( plan%device, procedureN )

    c_direction = INT( direction, C_INT )

    CALL cuda_set_device ( plan%device )
    c_status = cufftExecZ2Z ( plan%h, i_devc%ptr, o_devc%ptr, c_direction )
    IF( c_status /= cufft_success ) CALL stopgm(procedureN,"cufft error: "//TRIM(int2str( INT( c_status ) )),&
         __LINE__,__FILE__)

#else

    CALL stopgm(procedureN,"no cuda available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE cufft_execz2z

END MODULE cufft_utils
