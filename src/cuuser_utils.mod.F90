#include "cpmd_global.h"

MODULE cuuser_utils

  USE cuda_types,                      ONLY: cuda_memory_t,&
                                             cuda_stream_t
  USE cuda_utils,                      ONLY: cuda_set_device
  USE cuuser_interfaces,               ONLY: &
       CuUserDensitySum, CuUserEicalc, CuUserGetZ, CuUserIdentityDouble, &
       CuUserPack_x2y, CuUserPack_y2x, CuUserPhaseN, CuUserPointwise_CxR, &
       CuUserPutZ, CuUserSetBlock2Zero, CuUserSetBlock2ZeroScale, &
       CuUserSetPsi1StateG, CuUserSetPsi2StatesG, CuUserUnpack_x2y, &
       CuUserUnpack_y2x
  USE error_handling,                  ONLY: stopgm

  USE, INTRINSIC :: ISO_C_BINDING,  ONLY: C_INT,&
       C_PTR,&
       C_DOUBLE,&
       C_DOUBLE_COMPLEX,&
       C_CHAR, C_NULL_CHAR, C_NULL_PTR, C_BOOL
  USE kinds,                           ONLY: real_8


  IMPLICIT NONE

  PRIVATE

  PUBLIC :: CuUser_pack_x2y
  PUBLIC :: CuUser_unpack_x2y
  PUBLIC :: CuUser_pack_y2x
  PUBLIC :: CuUser_unpack_y2x
  PUBLIC :: CuUser_phasen
  PUBLIC :: CuUser_putz
  PUBLIC :: CuUser_getz
  PUBLIC :: CuUser_setblock2zero
  PUBLIC :: CuUser_setblock2zero_scale
  PUBLIC :: CuUser_identity
  PUBLIC :: CuUser_density_sum
  PUBLIC :: CuUser_pointwise_cxr
  PUBLIC :: CuUser_setpsi_1_state_g
  PUBLIC :: CuUser_setpsi_2_states_g
  PUBLIC :: CuUser_eicalc

CONTAINS


  SUBROUTINE CuUser_pack_x2y ( xf,yf,nrays,lda,jrxpl,sp5,maxfft,mproc,tr4a2a, stream )
    TYPE(cuda_memory_t)                      :: xf, yf
    INTEGER, INTENT(IN)                      :: nrays, lda
    TYPE(cuda_memory_t)                      :: jrxpl, sp5
    INTEGER, INTENT(IN)                      :: maxfft, mproc
    LOGICAL, INTENT(IN)                      :: tr4a2a
    TYPE(cuda_stream_t), INTENT(IN)          :: stream

    CHARACTER(*), PARAMETER                  :: procedureN = 'CuUser_pack_x2y'

    INTEGER(C_INT)                           :: c_lda, c_maxfft, c_mproc, &
                                                c_nrays
    LOGICAL(C_BOOL)                          :: c_tr4a2a

#if defined(_HAS_CUDA)

    IF( .NOT.xf%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT.yf%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT.jrxpl%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT.sp5%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT.stream%init ) CALL stopgm(procedureN,'stream not created',&
         __LINE__,__FILE__)

    IF(  xf%device /= yf%device .OR. &
         xf%device /= jrxpl%device .OR. &
         xf%device /= sp5%device .OR. &
         xf%device /= stream%device ) CALL stopgm(procedureN,'devices not consistent',&
         __LINE__,__FILE__)


    c_nrays = INT( nrays, C_INT )
    c_lda = INT( lda, C_INT )
    c_maxfft =INT( maxfft, C_INT )
    c_mproc = INT( mproc, C_INT )
    c_tr4a2a = LOGICAL( tr4a2a, C_BOOL )

    CALL cuda_set_device ( stream%device )
    CALL CuUserPack_x2y ( xf%ptr, yf%ptr, c_nrays, c_lda, jrxpl%ptr, sp5%ptr, c_maxfft, c_mproc, c_tr4a2a, stream%s )

#else

    CALL stopgm(procedureN,"no cuda available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE CuUser_pack_x2y


  SUBROUTINE CuUser_unpack_x2y ( xf,yf,m,lr1,lda,msp,lmsp,sp8,maxfft,mproc,tr4a2a, stream )
    TYPE(cuda_memory_t)                      :: xf, yf
    INTEGER, INTENT(IN)                      :: m, lr1, lda
    TYPE(cuda_memory_t)                      :: msp
    INTEGER, INTENT(IN)                      :: lmsp
    TYPE(cuda_memory_t)                      :: sp8
    INTEGER, INTENT(IN)                      :: maxfft, mproc
    LOGICAL, INTENT(IN)                      :: tr4a2a
    TYPE(cuda_stream_t), INTENT(IN)          :: stream

    CHARACTER(*), PARAMETER :: procedureN = 'CuUser_unpack_x2y'

    INTEGER(C_INT)                           :: c_lda, c_lmsp, c_lr1, c_m, &
                                                c_maxfft, c_mproc
    LOGICAL(C_BOOL)                          :: c_tr4a2a

#if defined(_HAS_CUDA)

    IF( .NOT.xf%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT.yf%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT.sp8%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT.msp%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT.stream%init ) CALL stopgm(procedureN,'stream not created',&
         __LINE__,__FILE__)

    IF(  xf%device /= yf%device .OR. &
         xf%device /= msp%device .OR. &
         xf%device /= sp8%device .OR. &
         xf%device /= stream%device ) CALL stopgm(procedureN,'devices not consistent',&
         __LINE__,__FILE__)

    c_m = INT( m, C_INT )
    c_lr1 = INT( lr1, C_INT )
    c_lda = INT( lda, C_INT )
    c_lmsp = INT( lmsp, C_INT )
    c_maxfft =INT( maxfft, C_INT )
    c_mproc = INT( mproc, C_INT )
    c_tr4a2a = LOGICAL( tr4a2a, C_BOOL )

    CALL cuda_set_device ( stream%device )
    CALL CuUserUnpack_x2y ( xf%ptr, yf%ptr, c_m, c_lr1, c_lda, msp%ptr, c_lmsp, sp8%ptr, c_maxfft, c_mproc, c_tr4a2a, stream%s )

#else

    CALL stopgm(procedureN,"no cuda available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE CuUser_unpack_x2y


  SUBROUTINE CuUser_pack_y2x ( xf, yf, m, lr1, lda, msp, lmsp, sp8, maxfft, mproc, tr4a2a, stream )
    TYPE(cuda_memory_t)                      :: xf, yf
    INTEGER, INTENT(IN)                      :: m, lr1, lda
    TYPE(cuda_memory_t)                      :: msp
    INTEGER, INTENT(IN)                      :: lmsp
    TYPE(cuda_memory_t)                      :: sp8
    INTEGER, INTENT(IN)                      :: maxfft, mproc
    LOGICAL, INTENT(IN)                      :: tr4a2a
    TYPE(cuda_stream_t), INTENT(IN)          :: stream

    CHARACTER(*), PARAMETER                  :: procedureN = 'CuUser_pack_y2x'

    INTEGER(C_INT)                           :: c_lda, c_lmsp, c_lr1, c_m, &
                                                c_maxfft, c_mproc
    LOGICAL(C_BOOL)                          :: c_tr4a2a

#if defined(_HAS_CUDA)

    IF( .NOT.xf%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT.yf%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT.msp%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT.sp8%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT.stream%init ) CALL stopgm(procedureN,'stream not created',&
         __LINE__,__FILE__)

    IF(  xf%device /= yf%device .OR. &
         xf%device /= msp%device .OR. &
         xf%device /= sp8%device .OR. &
         xf%device /= stream%device ) CALL stopgm(procedureN,'devices not consistent',&
         __LINE__,__FILE__)

    c_m = INT( m, C_INT )
    c_lr1 = INT( lr1, C_INT )
    c_lda = INT( lda, C_INT )
    c_lmsp = INT( lmsp, C_INT )
    c_maxfft =INT( maxfft, C_INT )
    c_mproc = INT( mproc, C_INT )
    c_tr4a2a = LOGICAL( tr4a2a, C_BOOL )

    CALL cuda_set_device ( stream%device )
    CALL CuUserPack_y2x ( xf%ptr, yf%ptr, c_m, c_lr1, c_lda, msp%ptr, c_lmsp, sp8%ptr, c_maxfft, c_mproc, c_tr4a2a, stream%s )

#else

    CALL stopgm(procedureN,"no cuda available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE CuUser_pack_y2x


  SUBROUTINE CuUser_unpack_y2x ( xf, yf, m, nrays, lda, jrxpl, sp5, maxfft, mproc, tr4a2a, stream )
    TYPE(cuda_memory_t)                      :: xf, yf
    INTEGER, INTENT(IN)                      :: m, nrays, lda
    TYPE(cuda_memory_t)                      :: jrxpl, sp5
    INTEGER, INTENT(IN)                      :: maxfft, mproc
    LOGICAL, INTENT(IN)                      :: tr4a2a
    TYPE(cuda_stream_t), INTENT(IN)          :: stream

    CHARACTER(*), PARAMETER :: procedureN = 'CuUser_unpack_y2x'

    INTEGER(C_INT)                           :: c_lda, c_m, c_maxfft, &
                                                c_mproc, c_nrays
    LOGICAL(C_BOOL)                          :: c_tr4a2a

#if defined(_HAS_CUDA)

    IF( .NOT.xf%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT.yf%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT.sp5%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT.jrxpl%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT.stream%init ) CALL stopgm(procedureN,'stream not created',&
         __LINE__,__FILE__)

    IF(  xf%device /= yf%device .OR. &
         xf%device /= jrxpl%device .OR. &
         xf%device /= sp5%device .OR. &
         xf%device /= stream%device ) CALL stopgm(procedureN,'devices not consistent',&
         __LINE__,__FILE__)

    c_m = INT( m, C_INT )
    c_nrays = INT( nrays, C_INT )
    c_lda = INT( lda, C_INT )
    c_maxfft =INT( maxfft, C_INT )
    c_mproc = INT( mproc, C_INT )
    c_tr4a2a = LOGICAL( tr4a2a, C_BOOL )

    CALL cuda_set_device ( stream%device )
    CALL CuUserUnpack_y2x ( xf%ptr, yf%ptr, c_m, c_nrays, c_lda, jrxpl%ptr, sp5%ptr, c_maxfft, c_mproc, c_tr4a2a, stream%s )

#else

    CALL stopgm(procedureN,"no cuda available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE CuUser_unpack_y2x

  SUBROUTINE CuUser_phasen ( f, kr1, kr2s, kr3s, n1u, n1o, nr2s, nr3s, stream )
    TYPE(cuda_memory_t)                      :: f
    INTEGER, INTENT(IN)                      :: kr1, kr2s, kr3s, n1u, n1o, &
                                                nr2s, nr3s
    TYPE(cuda_stream_t), INTENT(IN)          :: stream

    CHARACTER(*), PARAMETER                  :: procedureN = 'CuUser_phasen'

    INTEGER(C_INT)                           :: c_kr1, c_kr2s, c_kr3s, c_n1o, &
                                                c_n1u, c_nr2s, c_nr3s

#if defined(_HAS_CUDA)

    IF( .NOT.f%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT.stream%init ) CALL stopgm(procedureN,'stream not created',&
         __LINE__,__FILE__)

    IF( f%device /= stream%device ) CALL stopgm(procedureN,'devices not consistent',&
         __LINE__,__FILE__)

    c_kr1 = INT( kr1, C_INT )
    c_kr2s = INT( kr2s, C_INT )
    c_kr3s = INT( kr3s, C_INT )
    c_n1u = INT( n1u, C_INT )
    c_n1o = INT( n1o, C_INT )
    c_nr2s = INT( nr2s, C_INT )
    c_nr3s = INT( nr3s, C_INT )

    CALL cuda_set_device ( stream%device )
    CALL CuUserPhaseN ( f%ptr, c_kr1, c_kr2s, c_kr3s, c_n1u, c_n1o, c_nr2s, c_nr3s, stream%s )

#else

    CALL stopgm(procedureN,"no cuda available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE CuUser_phasen


  SUBROUTINE CuUser_setblock2zero ( odata, trans, N, M, ldbx, ldby, stream )
    TYPE(cuda_memory_t)                      :: odata
    CHARACTER(len=1), INTENT(IN)             :: trans
    INTEGER, INTENT(IN)                      :: n, m, ldbx, ldby
    TYPE(cuda_stream_t), INTENT(IN)          :: stream

    CHARACTER(*), PARAMETER :: procedureN = 'CuUser_setblock2zero'

    CHARACTER(KIND=C_CHAR, LEN=2)            :: c_trans
    INTEGER(C_INT)                           :: c_ldbx, c_ldby, c_m, c_n

#if defined(_HAS_CUDA)

    IF( .NOT.odata%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT.stream%init ) CALL stopgm(procedureN,'stream not created',&
         __LINE__,__FILE__)

    IF( odata%device /= stream%device ) CALL stopgm(procedureN,'devices not consistent',&
         __LINE__,__FILE__)

    c_trans = trans//C_NULL_CHAR
    c_n = INT( n, C_INT )
    c_m = INT( m, C_INT )
    c_ldbx = INT( ldbx, C_INT )
    c_ldby = INT( ldby, C_INT )

    CALL cuda_set_device ( stream%device )
    CALL CuUserSetBlock2Zero ( odata%ptr, c_trans, c_n, c_m, c_ldbx, c_ldby, stream%s )

#else

    CALL stopgm(procedureN,"no cuda available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE CuUser_setblock2zero


  SUBROUTINE CuUser_setblock2zero_scale ( odata, scale, trans, N, M, ldbx, ldby, stream )
    TYPE(cuda_memory_t)                      :: odata
    REAL(real_8), INTENT(IN)                 :: scale
    CHARACTER(len=1), INTENT(IN)             :: trans
    INTEGER, INTENT(IN)                      :: n, m, ldbx, ldby
    TYPE(cuda_stream_t), INTENT(IN)          :: stream

    CHARACTER(*), PARAMETER :: procedureN = 'CuUser_setblock2zero_scale'

    CHARACTER(KIND=C_CHAR, LEN=2)            :: c_trans
    INTEGER(C_INT)                           :: c_ldbx, c_ldby, c_m, c_n
    REAL(C_DOUBLE)                           :: c_scale

#if defined(_HAS_CUDA)

    IF( .NOT.odata%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT.stream%init ) CALL stopgm(procedureN,'stream not created',&
         __LINE__,__FILE__)

    IF( odata%device /= stream%device ) CALL stopgm(procedureN,'devices not consistent',&
         __LINE__,__FILE__)

    c_trans = trans//C_NULL_CHAR
    c_scale = REAL( scale, C_DOUBLE )
    c_n = INT( n, C_INT )
    c_m = INT( m, C_INT )
    c_ldbx = INT( ldbx, C_INT )
    c_ldby = INT( ldby, C_INT )

    CALL cuda_set_device ( stream%device )
    CALL CuUserSetBlock2ZeroScale ( odata%ptr, c_scale, c_trans, c_n, c_m, c_ldbx, c_ldby, stream%s )

#else

    CALL stopgm(procedureN,"no cuda available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE CuUser_setblock2zero_scale


  SUBROUTINE CuUser_putz ( a, b, krmin, krmax, kr, m, stream )
    ! isize in unit of COMPLEX(real_8)
    TYPE(cuda_memory_t)                      :: a, b
    INTEGER, INTENT(IN)                      :: krmin, krmax, kr, m
    TYPE(cuda_stream_t), INTENT(IN)          :: stream

    CHARACTER(*), PARAMETER                  :: procedureN = 'CuUser_putz'

    INTEGER(C_INT)                           :: c_kr, c_krmax, c_krmin, c_m

#if defined(_HAS_CUDA)

    IF( .NOT.a%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT.b%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT.stream%init ) CALL stopgm(procedureN,'stream not created',&
         __LINE__,__FILE__)

    IF( a%device /= stream%device ) CALL stopgm(procedureN,'devices not consistent',&
         __LINE__,__FILE__)

    c_krmin = INT( krmin, C_INT )
    c_krmax = INT( krmax, C_INT)
    c_kr    = INT( kr, C_INT)
    c_m     = INT( m, C_INT)
    !vw should check that complex(real_8) is the same as complex(c_double_complex)

    CALL cuda_set_device ( stream%device )
    CALL CuUserPutZ ( a%ptr, b%ptr, c_krmin, c_krmax, c_kr, c_m, stream%s )

#else

    CALL stopgm(procedureN,"no cuda available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE CuUser_putz


  SUBROUTINE CuUser_getz ( a, b, krmin, krmax, kr, m, stream )
    ! isize in unit of COMPLEX(real_8)
    TYPE(cuda_memory_t)                      :: a, b
    INTEGER, INTENT(IN)                      :: krmin, krmax, kr, m
    TYPE(cuda_stream_t), INTENT(IN)          :: stream

    CHARACTER(*), PARAMETER                  :: procedureN = 'CuUser_getz'

    INTEGER(C_INT)                           :: c_kr, c_krmax, c_krmin, c_m

#if defined(_HAS_CUDA)

    IF( .NOT.a%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT.b%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT.stream%init ) CALL stopgm(procedureN,'stream not created',&
         __LINE__,__FILE__)

    IF(  a%device /= b%device .OR. &
         a%device /= stream%device ) CALL stopgm(procedureN,'devices not consistent',&
         __LINE__,__FILE__)

    c_krmin = INT( krmin, C_INT )
    c_krmax = INT( krmax, C_INT)
    c_kr    = INT( kr, C_INT)
    c_m     = INT( m, C_INT)
    !vw should check that complex(real_8) is the same as complex(c_double_complex)

    CALL cuda_set_device ( stream%device )
    CALL CuUserGetZ ( a%ptr, b%ptr, c_krmin, c_krmax, c_kr, c_m, stream%s )

#else

    CALL stopgm(procedureN,"no cuda available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE CuUser_getz


  SUBROUTINE CuUser_identity ( a, m, stream )
    TYPE(cuda_memory_t)                      :: a
    INTEGER, INTENT(IN)                      :: m
    TYPE(cuda_stream_t), INTENT(IN)          :: stream

    CHARACTER(*), PARAMETER                  :: procedureN = 'CuUser_identity'

    INTEGER(C_INT)                           :: c_m

#if defined(_HAS_CUDA)

    IF( .NOT.a%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT.stream%init ) CALL stopgm(procedureN,'stream not created',&
         __LINE__,__FILE__)

    IF(  a%device /= stream%device ) CALL stopgm(procedureN,'devices not consistent',&
         __LINE__,__FILE__)

    c_m = INT( m, C_INT)

    CALL cuda_set_device ( stream%device )
    CALL CuUserIdentityDouble ( a%ptr, c_m, stream%s )

#else

    CALL stopgm(procedureN,"no cuda available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE CuUser_identity

  SUBROUTINE CuUser_density_sum ( alpha_real, alpha_imag, psi, rho, n, stream )
    REAL(real_8), INTENT(IN)                 :: alpha_real, alpha_imag
    TYPE(cuda_memory_t)                      :: psi, rho
    INTEGER, INTENT(IN)                      :: n
    TYPE(cuda_stream_t), INTENT(IN)          :: stream

    CHARACTER(*), PARAMETER :: procedureN = 'CuUser_density_sum'

    INTEGER(C_INT)                           :: c_n
    REAL(C_DOUBLE)                           :: c_alpha_imag, c_alpha_real

#if defined(_HAS_CUDA)

    IF( .NOT.stream%init ) CALL stopgm(procedureN,'stream not created',&
         __LINE__,__FILE__)

    IF( .NOT.psi%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT.rho%init ) CALL stopgm(procedureN,'stream not created',&
         __LINE__,__FILE__)

    IF( psi%device /= stream%device ) CALL stopgm(procedureN,'devices not consistent',&
         __LINE__,__FILE__)

    c_alpha_real = REAL( alpha_real, C_DOUBLE )
    c_alpha_imag = REAL( alpha_imag, C_DOUBLE )
    c_n = INT( n, C_INT )

    CALL cuda_set_device ( stream%device )
    CALL CuUserDensitySum ( c_alpha_real, c_alpha_imag, psi%ptr, rho%ptr, c_n, stream%s )

#else

    CALL stopgm(procedureN,"no cuda available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE CuUser_density_sum


  SUBROUTINE CuUser_pointwise_cxr ( xf, yf, n, stream )
    TYPE(cuda_memory_t)                      :: xf, yf
    INTEGER, INTENT(IN)                      :: n
    TYPE(cuda_stream_t), INTENT(IN)          :: stream

    CHARACTER(*), PARAMETER :: procedureN = 'CuUser_pointwise_cxr'

    INTEGER(C_INT)                           :: c_n

#if defined(_HAS_CUDA)

    IF( .NOT.stream%init ) CALL stopgm(procedureN,'stream not created',&
         __LINE__,__FILE__)

    IF( .NOT.xf%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT.yf%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( xf%device /= stream%device ) CALL stopgm(procedureN,'devices not consistent',&
         __LINE__,__FILE__)

    c_n = INT( n, C_INT )

    CALL cuda_set_device ( stream%device )
    CALL CuUserPointwise_CxR ( xf%ptr, yf%ptr, c_n, stream%s )

#else

    CALL stopgm(procedureN,"no cuda available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE CuUser_pointwise_cxr


  SUBROUTINE CuUser_setpsi_1_state_g ( alpha, c1, psi, jgw, nzfs, inzs, geq0, stream )
    COMPLEX(real_8), INTENT(IN)              :: alpha
    TYPE(cuda_memory_t), INTENT(IN)          :: c1
    TYPE(cuda_memory_t), INTENT(INOUT)       :: psi
    INTEGER, INTENT(IN)                      :: jgw
    TYPE(cuda_memory_t), INTENT(IN)          :: nzfs, inzs
    LOGICAL, INTENT(IN)                      :: geq0
    TYPE(cuda_stream_t), INTENT(IN)          :: stream

    CHARACTER(*), PARAMETER :: procedureN = 'CuUser_setpsi_1_state_g'

    INTEGER(C_INT)                           :: c_jgw
    LOGICAL(C_BOOL)                          :: c_geq0

#if defined(_HAS_CUDA)

    IF( .NOT.stream%init ) CALL stopgm(procedureN,'stream not created',&
         __LINE__,__FILE__)

    IF( .NOT.c1%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT.psi%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT.nzfs%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT.inzs%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( c1%device /= stream%device ) CALL stopgm(procedureN,'devices not consistent',&
         __LINE__,__FILE__)

    c_jgw = INT( jgw, C_INT )
    c_geq0 = LOGICAL( geq0, C_BOOL )

    CALL cuda_set_device ( stream%device )
    CALL CuUserSetPsi1StateG ( REAL(alpha), AIMAG(alpha), c1%ptr, psi%ptr, c_jgw, nzfs%ptr, inzs%ptr, c_geq0, stream%s )

#else

    CALL stopgm(procedureN,"no cuda available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE CuUser_setpsi_1_state_g


  SUBROUTINE CuUser_setpsi_2_states_g ( c1, c2, psi, jgw, nzfs, inzs, geq0, stream )
    TYPE(cuda_memory_t), INTENT(IN)          :: c1, c2
    TYPE(cuda_memory_t), INTENT(INOUT)       :: psi
    INTEGER, INTENT(IN)                      :: jgw
    TYPE(cuda_memory_t), INTENT(IN)          :: nzfs, inzs
    LOGICAL, INTENT(IN)                      :: geq0
    TYPE(cuda_stream_t), INTENT(IN)          :: stream

    CHARACTER(*), PARAMETER :: procedureN = 'CuUser_setpsi_2_states_g'

    INTEGER(C_INT)                           :: c_jgw
    LOGICAL(C_BOOL)                          :: c_geq0

#if defined(_HAS_CUDA)

    IF( .NOT.stream%init ) CALL stopgm(procedureN,'stream not created',&
         __LINE__,__FILE__)

    IF( .NOT.c1%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT.c2%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT.psi%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT.nzfs%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT.inzs%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( c1%device /= stream%device ) CALL stopgm(procedureN,'devices not consistent',&
         __LINE__,__FILE__)

    IF( c2%device /= stream%device ) CALL stopgm(procedureN,'devices not consistent',&
         __LINE__,__FILE__)

    c_jgw = INT( jgw, C_INT )
    c_geq0 = LOGICAL( geq0, C_BOOL )

    CALL cuda_set_device ( stream%device )
    CALL CuUserSetPsi2StatesG ( c1%ptr, c2%ptr, psi%ptr, c_jgw, nzfs%ptr, inzs%ptr, c_geq0, stream%s )

#else

    CALL stopgm(procedureN,"no cuda available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE CuUser_setpsi_2_states_g


  SUBROUTINE CuUser_eicalc ( eivps, eirop, nhg, nat, iatpt, inyh, ei1, ei2, ei3, vps, rhops, stream )
    TYPE(cuda_memory_t), INTENT(OUT)         :: eivps, eirop
    INTEGER, INTENT(IN)                      :: nhg, nat
    TYPE(cuda_memory_t), INTENT(IN)          :: iatpt, inyh, ei1, ei2, ei3, &
                                                vps, rhops
    TYPE(cuda_stream_t), INTENT(IN)          :: stream

    CHARACTER(*), PARAMETER                  :: procedureN = 'CuUser_eicalc'

    INTEGER(C_INT)                           :: c_nat, c_nhg

! Complex Real_8
! Integers
! Complex Real_8
! Real_8

#if defined(_HAS_CUDA)

    IF( .NOT.stream%init ) CALL stopgm(procedureN,'stream not created',&
         __LINE__,__FILE__)

    IF( .NOT.eivps%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT.eirop%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT.iatpt%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT.inyh%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT.ei1%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT.ei2%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT.ei3%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT.vps%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( .NOT.rhops%init ) CALL stopgm(procedureN,'memory not allocated',&
         __LINE__,__FILE__)

    IF( eivps%device /= stream%device ) CALL stopgm(procedureN,'devices not consistent',&
         __LINE__,__FILE__)

    c_nhg = INT( nhg, C_INT )
    c_nat = INT( nat, C_INT )

    CALL cuda_set_device ( stream%device )
    CALL CuUserEicalc ( eivps%ptr, eirop%ptr, c_nhg, c_nat, iatpt%ptr, inyh%ptr, ei1%ptr,&
         ei2%ptr, ei3%ptr, vps%ptr, rhops%ptr, stream%s )

#else

    CALL stopgm(procedureN,"no cuda available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE CuUser_eicalc

END MODULE cuuser_utils
