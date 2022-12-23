#include "cpmd_global.h"

MODULE cuuser_interfaces

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

  PUBLIC :: CuUserPack_x2y
  PUBLIC :: CuUserUnpack_x2y
  PUBLIC :: CuUserPack_y2x
  PUBLIC :: CuUserUnpack_y2x
  PUBLIC :: CuUserPhaseN
  PUBLIC :: CuUserPutZ
  PUBLIC :: CuUserGetZ
  PUBLIC :: CuUserSetBlock2Zero
  PUBLIC :: CuUserSetBlock2ZeroScale
  PUBLIC :: CuUserIdentityDouble
  PUBLIC :: CuUserDScal
  PUBLIC :: CuUserDensitySum
  PUBLIC :: CuUserPointwise_CxR
  PUBLIC :: CuUserSetPsi1StateG
  PUBLIC :: CuUserSetPsi2StatesG
  PUBLIC :: CuUserEicalc

  INTERFACE

     !void CuUserUnpack_x2y ( void *xf, void *yf, int m, int lr1, int lda, void *msp, int lmsp, void *sp8, int maxfft, int mproc, Bool tr4a2a, cudaStream_t stream )
     SUBROUTINE CuUserUnpack_x2y ( xf, yf, m, lr1, lda, msp, lmsp, sp8, maxfft, &
          & mproc, tr4a2a, stream ) BIND( C, name='CuUser_C_Unpack_x2y' )
       IMPORT :: C_INT, C_PTR, C_BOOL, cudaStream_t
       IMPLICIT NONE
       TYPE( C_PTR ), VALUE :: xf, yf
       INTEGER( C_INT ), VALUE :: m, lr1, lda, lmsp, maxfft, mproc
       TYPE( C_PTR ), VALUE :: msp, sp8
       LOGICAL( C_BOOL ), VALUE :: tr4a2a
       TYPE(cudaStream_t), VALUE :: stream
     END SUBROUTINE CuUserUnpack_x2y

     !void CuUserPack_x2y ( void *xf, void *yf, int nrays, int lda, void *jrxpl, void *sp5, int maxfft, int mproc, Bool tr4a2a, cudaStream_t stream  )
     SUBROUTINE CuUserPack_x2y ( xf, yf, nrays, lda, jrxpl, sp5, maxfft, mproc, &
          & tr4a2a, stream ) BIND( C, name='CuUser_C_Pack_x2y' )
       IMPORT :: C_INT, C_PTR, C_BOOL, cudastream_t
       IMPLICIT NONE
       TYPE( C_PTR ), VALUE :: xf, yf
       INTEGER( C_INT ), VALUE :: nrays, lda, maxfft, mproc
       TYPE( C_PTR ), VALUE :: jrxpl, sp5
       LOGICAL( C_BOOL ), VALUE :: tr4a2a
       TYPE(cudaStream_t), VALUE :: stream
     END SUBROUTINE CuUserPack_x2y

     !void CuUserPack_y2x ( void *xf, void *yf, int m, int lr1, int lda, void *msp, int lmsp, void *sp8, int maxfft, int mproc, Bool tr4a2a, cudaStream_t stream  )
     SUBROUTINE CuUserPack_y2x ( xf, yf, m, lr1, lda, msp, lmsp, sp8, maxfft, &
          & mproc, tr4a2a, stream ) BIND( C, name='CuUser_C_Pack_y2x' )
       IMPORT :: C_INT, C_PTR, C_BOOL, cudaStream_t
       IMPLICIT NONE
       TYPE( C_PTR ), VALUE :: xf, yf
       INTEGER( C_INT ), VALUE :: m, lr1, lda, lmsp,  maxfft, mproc
       TYPE( C_PTR ), VALUE :: msp, sp8
       LOGICAL( C_BOOL ), VALUE :: tr4a2a
       TYPE(cudaStream_t), VALUE :: stream
     END SUBROUTINE CuUserPack_y2x

     !void CuUserUnpack_y2x ( void *xf, void *yf, int m, int nrays, int lda, void *jrxpl, void *sp5, int maxfft, int mproc, Bool tr4a2a, cudaStream_t stream  )
     SUBROUTINE CuUserUnpack_y2x ( xf, yf, m, nrays, lda, jrxpl, sp5, maxfft, mproc, &
          & tr4a2a, stream ) BIND( C, name='CuUser_C_Unpack_y2x' )
       IMPORT :: C_INT, C_PTR, C_BOOL, cudaStream_t
       IMPLICIT NONE
       TYPE( C_PTR ), VALUE :: xf, yf
       INTEGER( C_INT ), VALUE ::  m, nrays, lda, maxfft, mproc
       TYPE( C_PTR ), VALUE :: jrxpl, sp5
       LOGICAL( C_BOOL ), VALUE :: tr4a2a
       TYPE(cudaStream_t), VALUE :: stream
     END SUBROUTINE CuUserUnpack_y2x

     !void CuUserPhaseN ( void *f, int kr1, int kr2s, int kr3s, int n1u, int n1o, int nr2s, int nr3s, cudaStream_t stream  )
     SUBROUTINE CuUserPhaseN ( f, kr1, kr2s, kr3s, n1u, n1o, nr2s, nr3s, stream ) BIND( C, NAME='CuUser_C_PhaseN' )
       IMPORT :: C_INT, C_PTR, cudaStream_t
       IMPLICIT NONE
       TYPE( C_PTR ), VALUE :: f
       INTEGER ( C_INT ), VALUE :: kr1, kr2s, kr3s, n1u, n1o, nr2s, nr3s
       TYPE(cudaStream_t), VALUE :: stream
     END SUBROUTINE CuUserPhaseN

     !void CuUserPutZ( void *a, void *b, int krmin, int krmax, int kr, int m, cudaStream_t stream  )
     SUBROUTINE CuUserPutZ ( a, b, krmin, krmax, kr, m, stream ) BIND( C, NAME='CuUser_C_PutZ' )
       IMPORT :: C_INT, C_PTR, cudaStream_t
       IMPLICIT NONE
       TYPE( C_PTR ), VALUE :: a, b
       INTEGER( C_INT ), VALUE :: krmin, krmax, kr, m
       TYPE(cudaStream_t), VALUE :: stream
     END SUBROUTINE CuUserPutZ

     !void CuUserGetZ( void *a, void *b, int krmin, int krmax, int kr, int m, cudaStream_t stream  )
     SUBROUTINE CuUserGetZ ( a, b, krmin, krmax, kr, m, stream ) BIND( C, NAME='CuUser_C_GetZ' )
       IMPORT :: C_INT, C_PTR, cudaStream_t
       IMPLICIT NONE
       TYPE( C_PTR ), VALUE :: a, b
       INTEGER( C_INT ), VALUE :: krmin, krmax, kr, m
       TYPE(cudaStream_t), VALUE :: stream
     END SUBROUTINE CuUserGetZ

     !void CuUserSetBlockToZero (void *a, char *trans, int N, int M, int LDBX, int LDBY, cudaStream_t stream  )
     SUBROUTINE CuUserSetBlock2Zero ( odata, trans, n, m, ldbx, ldby, stream ) BIND( C, NAME="CuUser_C_SetBlock2Zero" )
       IMPORT :: C_INT, C_PTR, C_CHAR, cudaStream_t
       IMPLICIT NONE
       TYPE( C_PTR ), VALUE :: odata
       CHARACTER ( C_CHAR ), DIMENSION( * ) :: trans
       INTEGER ( C_INT ), VALUE :: n, m, ldbx, ldby
       TYPE(cudaStream_t), VALUE :: stream
     END SUBROUTINE CuUserSetBlock2Zero

     !void CuUserSetBlockToZeroScale (void *a, double scale, char *trans, int N, int M, int LDBX, int LDBY, cudaStream_t stream  )
     SUBROUTINE CuUserSetBlock2ZeroScale ( odata, scale, trans, n, m, ldbx, ldby, stream ) &
          BIND( C, NAME="CuUser_C_SetBlock2ZeroScale" )
       IMPORT :: C_INT, C_DOUBLE, C_PTR, C_CHAR, cudaStream_t
       IMPLICIT NONE
       TYPE( C_PTR ), VALUE :: odata
       REAL( C_DOUBLE ), VALUE :: scale
       CHARACTER ( C_CHAR ), DIMENSION( * ) :: trans
       INTEGER ( C_INT ), VALUE :: n, m, ldbx, ldby
       TYPE(cudaStream_t), VALUE :: stream
     END SUBROUTINE CuUserSetBlock2ZeroScale

     !void CuUser_C_IdentityDouble( double *a, int m, cudaStream_t stream )
     SUBROUTINE CuUserIdentityDouble ( a, m, stream ) BIND( C, NAME='CuUser_C_IdentityDouble' )
       IMPORT :: C_INT, C_PTR, cudaStream_t
       IMPLICIT NONE
       TYPE( C_PTR ), VALUE :: a
       INTEGER ( C_INT ), VALUE :: m
       TYPE(cudaStream_t), VALUE :: stream
     END SUBROUTINE CuUserIdentityDouble

     !void CuUser_C_DScal( double *a, int m, double scale, cudaStream_t stream )
     SUBROUTINE CuUserDScal ( a, m, scale, stream ) BIND( C, NAME='CuUser_C_DScal' )
       IMPORT :: C_INT, C_DOUBLE, C_PTR, cudaStream_t
       IMPLICIT NONE
       TYPE( C_PTR ), VALUE :: a
       INTEGER ( C_INT ), VALUE :: m
       REAL( C_DOUBLE ), VALUE :: scale
       TYPE(cudaStream_t), VALUE :: stream
     END SUBROUTINE CuUserDScal

     !void CuUser_C_Build_Density_Sum ( double alpha_real, double alpha_imag, void *psi, void *rho, int n, cudaStream_t stream )
     SUBROUTINE CuUserDensitySum ( alpha_real, alpha_imag, psi, rho, n, stream ) BIND( C, name='CuUser_C_Build_Density_Sum' )
       IMPORT :: C_INT, C_DOUBLE, C_PTR, C_BOOL, cudastream_t
       IMPLICIT NONE
       REAL( C_DOUBLE ), VALUE :: alpha_real, alpha_imag
       TYPE( C_PTR ), VALUE :: psi, rho
       INTEGER( C_INT ), VALUE :: n
       TYPE(cudaStream_t), VALUE :: stream
     END SUBROUTINE CuUserDensitySum

     !void CuUser_C_Build_Pointwise_CxR ( void *xf, void *yf, int n, cudaStream_t stream )
     SUBROUTINE CuUserPointwise_CxR ( xf, yf, n, stream ) BIND( C, name='CuUser_C_Build_Pointwise_CxR' )
       IMPORT :: C_INT, C_PTR, C_BOOL, cudastream_t
       IMPLICIT NONE
       TYPE( C_PTR ), VALUE :: xf, yf
       INTEGER( C_INT ), VALUE :: n
       TYPE(cudaStream_t), VALUE :: stream
     END SUBROUTINE CuUserPointwise_CxR

     !void CuUserSetPsi1StateG ( double alpha_real, double alpha_imag, void *c1, void *psi, int jgw, void *nzfs, void *inzs, Bool geq0, cudaStream_t stream )
     SUBROUTINE CuUserSetPsi1StateG ( alpha_real, alpha_imag, c1, psi, jgw, nzfs, inzs, geq0, stream )&
          BIND( C, name='CuUser_C_Set_Psi_1_Stage_G' )
       IMPORT :: C_DOUBLE, C_INT, C_PTR, C_BOOL, cudastream_t
       IMPLICIT NONE
       REAL( C_DOUBLE ), VALUE :: alpha_real, alpha_imag
       TYPE( C_PTR ), VALUE :: c1, psi
       INTEGER( C_INT ), VALUE :: jgw
       TYPE( C_PTR ), VALUE :: nzfs, inzs
       LOGICAL( C_BOOL ), VALUE :: geq0
       TYPE(cudaStream_t), VALUE :: stream
     END SUBROUTINE CuUserSetPsi1StateG

     !void CuUserSetPsi2StatesG ( void *c1, void *c2, void *psi, int jgw, void *nzfs, void *inzs, Bool geq0, cudaStream_t stream )
     SUBROUTINE CuUserSetPsi2StatesG ( c1, c2, psi, jgw, nzfs, inzs, geq0, stream )&
          BIND( C, name='CuUser_C_Set_Psi_2_Stages_G' )
       IMPORT :: C_INT, C_PTR, C_BOOL, cudastream_t
       IMPLICIT NONE
       TYPE( C_PTR ), VALUE :: c1, c2, psi
       INTEGER( C_INT ), VALUE :: jgw
       TYPE( C_PTR ), VALUE :: nzfs, inzs
       LOGICAL( C_BOOL ), VALUE :: geq0
       TYPE(cudaStream_t), VALUE :: stream
     END SUBROUTINE CuUserSetPsi2StatesG

     !void CuUser_C_CuUserEicalc ( void *eivps, void *eirop, int nhg, int nat, void *iatpt, void *inyh, void *ei1, void *ei2, void *ei3, void *vps, void *rhops, cudaStream_t stream )
     SUBROUTINE CuUserEicalc ( eivps, eirop, nhg, nat, iatpt, inyh, ei1, ei2, ei3, vps, rhops, stream )&
          BIND( C, name='CuUser_C_Eicalc' )
       IMPORT :: C_INT, C_PTR, cudastream_t
       IMPLICIT NONE
       TYPE( C_PTR ), VALUE :: eivps, eirop
       INTEGER( C_INT ), VALUE :: nhg, nat
       TYPE( C_PTR ), VALUE :: iatpt, inyh
       TYPE( C_PTR ), VALUE :: ei1, ei2, ei3
       TYPE( C_PTR ), VALUE :: vps, rhops
       TYPE(cudaStream_t), VALUE :: stream
     END SUBROUTINE CuUserEicalc

  END INTERFACE

END MODULE cuuser_interfaces
