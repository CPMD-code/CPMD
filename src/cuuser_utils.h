#if defined(__HAS_CUDA)
#define __CUDA
#endif

#if defined(__CUDA)

#ifndef __CUUSER_C_UTILS_HEADER___
#define __CUUSER_C_UTILS_HEADER___

#include <stdbool.h>
#include <cstdio>

#include <cuda.h>
#include <cuComplex.h>

// Max number of threads per block
// This should be GPU dependent
#define maxTxB 1024

// Pseudo-Optimal number of threads per block
// This should be GPU and function dependent
#define optTxB  128

#ifdef __cplusplus
extern "C" {
#endif

  __global__ void CuUser_Kernel_SetBlock2Zero (cuDoubleComplex *a, int N, int M, int LDBX, int LDBY );
  __global__ void CuUser_Kernel_SetBlock2ZeroScale (cuDoubleComplex *a, double scale, int N, int M, int LDBX, int LDBY );

  __global__ void CuUser_Kernel_Identity_double( double *devMatrix, int m );

  __global__ void CuUser_Kernel_MatMov ( int n, int m, cuDoubleComplex *a, int lda, cuDoubleComplex *b, int ldb);
  __global__ void CuUser_Kernel_Zeroing ( cuDoubleComplex *a, int size);

  __global__ void CuUser_Kernel_PhaseN ( cuDoubleComplex *f, int kr1, int kr2s, int kr3s, int n1u, int n1o, int nr2s, int nr3s );

  __global__ void CuUser_Kernel_Unpack_x2y_8 ( cuDoubleComplex *xf, cuDoubleComplex *yf, int m, int lr1, int lda, int *msp, int lmsp, int *sp8, int maxfft, int mproc );
  __global__ void CuUser_Kernel_Unpack_x2y_4 ( cuComplex *xf, cuDoubleComplex *yf, int m, int lr1, int lda, int *msp, int lmsp, int *sp8, int maxfft, int mproc );

  __global__ void CuUser_Kernel_Unpack_y2x_8 ( cuDoubleComplex *xf, cuDoubleComplex *yf, int m, int nrays, int lda, int *jrxpl, int *sp5, int maxfft, int mproc );
  __global__ void CuUser_Kernel_Unpack_y2x_4 ( cuDoubleComplex *xf, cuComplex *yf, int m, int nrays, int lda, int *jrxpl, int *sp5, int maxfft, int mproc );

  __global__ void CuUser_Kernel_Pack_x2y_8 ( cuDoubleComplex *xf, cuDoubleComplex *yf, int nrays, int lda, int *jrxpl, int *sp5, int maxfft, int mproc );
  __global__ void CuUser_Kernel_Pack_x2y_4 ( cuDoubleComplex *xf, cuComplex *yf, int nrays, int lda, int *jrxpl, int *sp5, int maxfft, int mproc );

  __global__ void CuUser_Kernel_Pack_y2x_8 ( cuDoubleComplex *xf, cuDoubleComplex *yf, int m, int lr1, int lda, int *msp, int lmsp, int *sp8, int maxfft, int mproc );
  __global__ void CuUser_Kernel_Pack_y2x_4 ( cuComplex *xf, cuDoubleComplex *yf, int m, int lr1, int lda, int *msp, int lmsp, int *sp8, int maxfft, int mproc );

  __global__ void CuUser_Kernel_Build_Density_Sum ( double alpha_real, double alpha_imag, cuDoubleComplex *psi, double *rho, int n );

  __global__ void CuUser_Kernel_Pointwise_CxR ( cuDoubleComplex *xf_p, double *yf_p, int n );

  __global__ void CuUser_Kernel_Set_Psi_1_Stage_G ( double alpha_real, double alpha_imag, cuDoubleComplex * c1_p, cuDoubleComplex * psi_p, int jgw, int *nzfs_p, int *inzs_p, _Bool geq0 );
  __global__ void CuUser_Kernel_Set_Psi_2_Stages_G ( cuDoubleComplex * c1_p, cuDoubleComplex * c2_p, cuDoubleComplex * psi_p, int jgw, int *nzfs_p, int *inzs_p, _Bool geq0 );

  __global__ void CuUser_Kernel_Eicalc ( cuDoubleComplex * eivps_p, cuDoubleComplex * eirop_p, int nhg, int nat, int *iatpt_p, int *inyh_p, cuDoubleComplex *ei1_p, cuDoubleComplex *ei2_p, cuDoubleComplex *ei3_p, double *vps_p, double *rhops_p );



  void CuUser_C_SetBlock2Zero (void *a, char *trans, int N, int M, int LDBX, int LDBY, cudaStream_t stream );
  void CuUser_C_SetBlock2ZeroScale (void *a, double scale, char *trans, int N, int M, int LDBX, int LDBY, cudaStream_t stream );

  void CuUser_C_IdentityDouble( void *a, int m, cudaStream_t stream );

  void CuUser_C_PutZ( void *a, void *b, int krmin, int krmax, int kr, int m, cudaStream_t stream );

  void CuUser_C_GetZ( void *a, void *b, int krmin, int krmax, int kr, int m, cudaStream_t stream );

  void CuUser_C_PhaseN ( void *f, int kr1, int kr2s, int kr3s, int n1u, int n1o, int nr2s, int nr3s, cudaStream_t stream );

  void CuUser_C_Unpack_x2y ( void *xf, void *yf, int m, int lr1, int lda, void *msp, int lmsp, void *sp8, int maxfft, int mproc, _Bool tr4a2a, cudaStream_t stream );

  void CuUser_C_Unpack_y2x ( void *xf, void *yf, int m, int nrays, int lda, void *jrxpl, void *sp5, int maxfft, int mproc, _Bool tr4a2a, cudaStream_t stream );

  void CuUser_C_Pack_x2y ( void *xf, void *yf, int nrays, int lda, void *jrxpl, void *sp5, int maxfft, int mproc, _Bool tr4a2a, cudaStream_t stream );

  void CuUser_C_Pack_y2x ( void *xf, void *yf, int m, int lr1, int lda, void *msp, int lmsp, void *sp8, int maxfft, int mproc, _Bool tr4a2a, cudaStream_t stream );

  void CuUser_C_Build_Density_Sum ( double alpha_real, double alpha_imag, void *psi, void *rho, int n, cudaStream_t stream );

  void CuUser_C_Build_Pointwise_CxR ( void *xf, void *yf, int n, cudaStream_t stream );

  void CuUser_C_Set_Psi_1_Stage_G ( double alpha_real, double alpha_imag, void *c1, void *psi, int jgw, void *nzfs, void *inzs, _Bool geq0, cudaStream_t stream );

  void CuUser_C_Set_Psi_2_Stages_G ( void *c1, void *c2, void *psi, int jgw, void *nzfs, void *inzs, _Bool geq0, cudaStream_t stream );

  void CuUser_C_Eicalc ( void *eivps, void *eirop, int nhg, int nat, void *iatpt, void *inyh, void *ei1, void *ei2, void *ei3, void *vps, void *rhops, cudaStream_t stream );

#ifdef __cplusplus
}
#endif

#endif //__CUUSER_C_UTILS_HEADER___
#endif //__CUDA
