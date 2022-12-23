#include "cuuser_utils.h"

#if defined(__CUDA)

#ifndef __CUUSER_C_UTILS_SOURCES___
#define __CUUSER_C_UTILS_SOURCES___

void CuUser_C_Synchronize()
{
    cudaError_t cudaerr = cudaDeviceSynchronize();
// TODO warning: comparison between ‘cudaError_t {aka enum cudaError}’ and ‘enum cudaError_enum’ [-Wenum-compare]
//    if (cudaerr != CUDA_SUCCESS)
//        printf("kernel launch failed with error \"%s\".\n",cudaGetErrorString(cudaerr));

    printf("DEBUG CUDA kernel - Synchronized, remove once optimized\n");
}

void CuUser_C_SetBlock2Zero (void *a, char *trans, int N, int M, int LDBX, int LDBY, cudaStream_t stream )
{
    cuDoubleComplex *a_p = (cuDoubleComplex *) a;

    dim3 dimBlock(1,optTxB,1);  // TODO acm: probably this can be better blocked
    dim3 dimGrid(int(ceil(float(LDBX)/dimBlock.x)),int(ceil(float(LDBY)/dimBlock.y)),1);

    if ( *trans == 'n' || *trans == 'N' )
        CuUser_Kernel_SetBlock2Zero<<<dimGrid,dimBlock,0,stream>>> ( a_p, N, M, LDBX, LDBY);
    else
        CuUser_Kernel_SetBlock2Zero<<<dimGrid,dimBlock,0,stream>>> ( a_p, M, N, LDBX, LDBY);
}

void CuUser_C_SetBlock2ZeroScale (void *a, double scale, char *trans, int N, int M, int LDBX, int LDBY, cudaStream_t stream )
{
    cuDoubleComplex *a_p = (cuDoubleComplex *) a;

    dim3 dimBlock(1,maxTxB,1);  // TODO acm: probably this can be better blocked
    dim3 dimGrid(int(ceil(float(LDBX)/dimBlock.x)),int(ceil(float(LDBY)/dimBlock.y)),1);

    if ( *trans == 'n' || *trans == 'N' )
      CuUser_Kernel_SetBlock2ZeroScale<<<dimGrid,dimBlock,0,stream>>> ( a_p, scale, N, M, LDBX, LDBY);
    else
      CuUser_Kernel_SetBlock2ZeroScale<<<dimGrid,dimBlock,0,stream>>> ( a_p, scale, M, N, LDBX, LDBY);
}

void CuUser_C_IdentityDouble( void *a, int m, cudaStream_t stream )
{
    double *a_p = (double *) a;

    dim3 dimBlock(1,optTxB,1);  // TODO acm: probably this can be better blocked
    dim3 dimGrid(int(ceil(float(m)/dimBlock.x)),int(ceil(float(m)/dimBlock.y)),1);

    CuUser_Kernel_Identity_double<<<dimGrid,dimBlock,0,stream>>> ( a_p, m );
}

void CuUser_C_PutZ( void *a, void *b, int krmin, int krmax, int kr, int m, cudaStream_t stream )
{
    // TODO: Same blocking for the two kernels would probably improve performance
    cuDoubleComplex *a_p = (cuDoubleComplex *) a;
    cuDoubleComplex *b_p = (cuDoubleComplex *) b;

    const int size = m * kr;

    dim3 dimBlock(optTxB,1,1);
    dim3 dimGrid(int(ceil(float(size)/dimBlock.x)),1,1);

    CuUser_Kernel_Zeroing<<<dimGrid,dimBlock,0,stream>>> ( b_p, size );

    const int n = krmax - krmin + 1;

    dim3 dimBlock2(1,optTxB,1);  // TODO acm: probably this can be better blocked
    dim3 dimGrid2(int(ceil(float(n)/dimBlock2.x)),int(ceil(float(m)/dimBlock2.y)),1);

    CuUser_Kernel_MatMov<<<dimGrid2,dimBlock2,0,stream>>> ( n, m, a_p, n, &b_p[krmin-1], kr );
}

void CuUser_C_GetZ( void *a, void *b, int krmin, int krmax, int kr, int m, cudaStream_t stream )
{
    cuDoubleComplex *a_p = (cuDoubleComplex *) a;
    cuDoubleComplex *b_p = (cuDoubleComplex *) b;

    const int n = krmax - krmin + 1;

    dim3 dimBlock(1,optTxB,1);  // TODO acm: probably this can be better blocked
    dim3 dimGrid(int(ceil(float(n)/dimBlock.x)),int(ceil(float(m)/dimBlock.y)),1);

    CuUser_Kernel_MatMov<<<dimGrid,dimBlock,0,stream>>> ( n, m, &a_p[krmin-1], kr, b_p, n );
}

void CuUser_C_PhaseN ( void *f, int kr1, int kr2s, int kr3s, int n1u, int n1o, int nr2s, int nr3s, cudaStream_t stream )
{
    cuDoubleComplex *f_p = (cuDoubleComplex *) f;

    dim3 dimBlock(1,optTxB,1);
    dim3 dimGrid(int(ceil(float(n1o - n1u + 1)/dimBlock.x)),int(ceil(float(nr2s)/dimBlock.y)),int(ceil(float(nr3s)/dimBlock.z)));

    CuUser_Kernel_PhaseN<<<dimGrid,dimBlock,0,stream>>> ( f_p, kr1, kr2s, kr3s, n1u, n1o, nr2s, nr3s );
}

void CuUser_C_Unpack_x2y ( void *xf, void *yf, int m, int lr1, int lda, void *msp, int lmsp, void *sp8, int maxfft, int mproc, _Bool tr4a2a, cudaStream_t stream )
{
    cuDoubleComplex *yf_p = (cuDoubleComplex *) yf;

    // TODO acm: is this necessary?
    int *msp_p = (int *) msp;
    int *sp8_p = (int *) sp8;

    // TODO: Same blocking for the two kernels would probably improve performance
    dim3 dimBlock2(optTxB,1,1);
    dim3 dimGrid2(int(ceil(float(maxfft)/dimBlock2.x)),1,1);

    CuUser_Kernel_Zeroing<<<dimGrid2,dimBlock2,0,stream>>> ( yf_p, maxfft );

    dim3 dimBlock(1,optTxB,1);
    dim3 dimGrid(int(ceil(float(mproc)/dimBlock.x)),int(ceil(float(maxfft)/dimBlock.y)),int(ceil(float(lr1)/dimBlock.z)));

    if ( tr4a2a )
    {
        cuComplex *xf4_p = (cuComplex *) xf;
        CuUser_Kernel_Unpack_x2y_4<<<dimGrid,dimBlock,0,stream>>> ( xf4_p, yf_p, m, lr1, lda, msp_p, lmsp, sp8_p, maxfft, mproc );
    }
    else
    {
        cuDoubleComplex *xf8_p = (cuDoubleComplex *) xf;
        CuUser_Kernel_Unpack_x2y_8<<<dimGrid,dimBlock,0,stream>>> ( xf8_p, yf_p, m, lr1, lda, msp_p, lmsp, sp8_p, maxfft, mproc );
    }
}

void CuUser_C_Unpack_y2x ( void *xf, void *yf, int m, int nrays, int lda, void *jrxpl, void *sp5, int maxfft, int mproc, _Bool tr4a2a, cudaStream_t stream )
{
    cuDoubleComplex *xf_p = (cuDoubleComplex *) xf;

    // TODO acm: is this necessary?
    int *jrxpl_p = (int *) jrxpl;
    int *sp5_p = (int *) sp5;

    dim3 dimBlock(1,optTxB,1);  // TODO acm: probably this can be better blocked
    dim3 dimGrid(int(ceil(float(mproc)/dimBlock.x)),int(ceil(float(maxfft)/dimBlock.y)),1);

    if ( tr4a2a )
    {
        cuComplex *yf4_p = (cuComplex *) yf;
        CuUser_Kernel_Unpack_y2x_4<<<dimGrid,dimBlock,0,stream>>> ( xf_p, yf4_p, m, nrays, lda, jrxpl_p, sp5_p, maxfft, mproc );
    }
    else
    {
        cuDoubleComplex *yf8_p = (cuDoubleComplex *) yf;
        CuUser_Kernel_Unpack_y2x_8<<<dimGrid,dimBlock,0,stream>>> ( xf_p, yf8_p, m, nrays, lda, jrxpl_p, sp5_p, maxfft, mproc );
    }
}

void CuUser_C_Pack_x2y ( void *xf, void *yf, int nrays, int lda, void *jrxpl, void *sp5, int maxfft, int mproc, _Bool tr4a2a, cudaStream_t stream )
{
    cuDoubleComplex *xf_p  = (cuDoubleComplex *) xf;

    int *jrxpl_p = (int *) jrxpl;
    int *sp5_p   = (int *) sp5;

    dim3 dimBlock(1,optTxB,1);  // TODO acm: probably this can be better blocked
    dim3 dimGrid(int(ceil(float(mproc)/dimBlock.x)),int(ceil(float(maxfft)/dimBlock.y)),1);

    if ( tr4a2a )
    {
        cuComplex *yf4_p = (cuComplex *) yf;
        CuUser_Kernel_Pack_x2y_4<<<dimGrid,dimBlock,0,stream>>> ( xf_p, yf4_p, nrays, lda, jrxpl_p, sp5_p, maxfft, mproc );
    }
    else
    {
        cuDoubleComplex *yf8_p = (cuDoubleComplex *) yf;
        CuUser_Kernel_Pack_x2y_8<<<dimGrid,dimBlock,0,stream>>> ( xf_p, yf8_p, nrays, lda, jrxpl_p, sp5_p, maxfft, mproc );
    }
}

void CuUser_C_Pack_y2x ( void *xf, void *yf, int m, int lr1, int lda, void *msp, int lmsp, void *sp8, int maxfft, int mproc, _Bool tr4a2a, cudaStream_t stream )
{
    cuDoubleComplex *yf_p = (cuDoubleComplex *) yf;

    int *msp_p = (int *) msp;
    int *sp8_p = (int *) sp8;

    dim3 dimBlock(1,optTxB,1);
    dim3 dimGrid(int(ceil(float(mproc)/dimBlock.x)),int(ceil(float(maxfft)/dimBlock.y)),int(ceil(float(lr1)/dimBlock.z)));

    if ( tr4a2a )
    {
        cuComplex *xf4_p = (cuComplex *) xf;
        CuUser_Kernel_Pack_y2x_4<<<dimGrid,dimBlock,0,stream>>> ( xf4_p, yf_p, m, lr1, lda, msp_p, lmsp, sp8_p, maxfft, mproc );
    }
    else
    {
        cuDoubleComplex *xf8_p = (cuDoubleComplex *) xf;
        CuUser_Kernel_Pack_y2x_8<<<dimGrid,dimBlock,0,stream>>> ( xf8_p, yf_p, m, lr1, lda, msp_p, lmsp, sp8_p, maxfft, mproc );
    }
}

void CuUser_C_Build_Density_Sum ( double alpha_real, double alpha_imag, void *psi, void *rho, int n, cudaStream_t stream )
{
    cuDoubleComplex *psi_p = (cuDoubleComplex *) psi;

    double *rho_p = (double *) rho;

    dim3 dimBlock(optTxB,1,1);
    dim3 dimGrid(int(ceil(float(n)/dimBlock.x)),1,1);

    CuUser_Kernel_Build_Density_Sum<<<dimGrid,dimBlock,0,stream>>> ( alpha_real, alpha_imag, psi_p, rho_p, n );
}

void CuUser_C_Build_Pointwise_CxR ( void *xf, void *yf, int n, cudaStream_t stream )
{
    cuDoubleComplex *xf_p = (cuDoubleComplex *) xf;

    double *yf_p = (double *) yf;

    dim3 dimBlock(optTxB,1,1);
    dim3 dimGrid(int(ceil(float(n)/dimBlock.x)),1,1);

    CuUser_Kernel_Pointwise_CxR<<<dimGrid,dimBlock,0,stream>>> ( xf_p, yf_p, n );
}

void CuUser_C_Set_Psi_1_Stage_G ( double alpha_real, double alpha_imag, void *c1, void *psi, int jgw, void *nzfs, void *inzs, _Bool geq0, cudaStream_t stream )
{
    cuDoubleComplex *c1_p  = (cuDoubleComplex *) c1;
    cuDoubleComplex *psi_p = (cuDoubleComplex *) psi;

    int *nzfs_p = (int *) nzfs;
    int *inzs_p = (int *) inzs;

    dim3 dimBlock(optTxB,1,1);
    dim3 dimGrid(int(ceil(float(jgw)/dimBlock.x)),1,1);

    CuUser_Kernel_Set_Psi_1_Stage_G<<<dimGrid,dimBlock,0,stream>>> ( alpha_real, alpha_imag, c1_p, psi_p, jgw, nzfs_p, inzs_p, geq0 );
}

void CuUser_C_Set_Psi_2_Stages_G ( void *c1, void *c2, void *psi, int jgw, void *nzfs, void *inzs, _Bool geq0, cudaStream_t stream )
{
    cuDoubleComplex *c1_p  = (cuDoubleComplex *) c1;
    cuDoubleComplex *c2_p  = (cuDoubleComplex *) c2;
    cuDoubleComplex *psi_p = (cuDoubleComplex *) psi;

    int *nzfs_p = (int *) nzfs;
    int *inzs_p = (int *) inzs;

    dim3 dimBlock(optTxB,1,1);
    dim3 dimGrid(int(ceil(float(jgw)/dimBlock.x)),1,1);

    CuUser_Kernel_Set_Psi_2_Stages_G<<<dimGrid,dimBlock,0,stream>>> ( c1_p, c2_p, psi_p, jgw, nzfs_p, inzs_p, geq0 );
}

void CuUser_C_Eicalc ( void *eivps, void *eirop, int nhg, int nat, void *iatpt, void *inyh, void *ei1, void *ei2, void *ei3, void *vps, void *rhops, cudaStream_t stream )
{
    cuDoubleComplex *eivps_p  = (cuDoubleComplex *) eivps;
    cuDoubleComplex *eirop_p  = (cuDoubleComplex *) eirop;

    int *iatpt_p = (int *) iatpt;
    int *inyh_p  = (int *) inyh;

    cuDoubleComplex *ei1_p  = (cuDoubleComplex *) ei1;
    cuDoubleComplex *ei2_p  = (cuDoubleComplex *) ei2;
    cuDoubleComplex *ei3_p  = (cuDoubleComplex *) ei3;

    double *vps_p   = (double *) vps;
    double *rhops_p = (double *) rhops;

    dim3 dimBlock(optTxB,1,1); // TODO acm: probably this can be better blocked
    dim3 dimGrid(int(ceil(float(nhg)/dimBlock.x)),int(ceil(float(nat)/dimBlock.y)),1);

    CuUser_Kernel_Eicalc<<<dimGrid,dimBlock,0,stream>>> ( eivps_p, eirop_p, nhg, nat, iatpt_p, inyh_p, ei1_p, ei2_p, ei3_p, vps_p, rhops_p);
}

#endif // __CUUSER_C_UTILS_SOURCES___
#endif // __CUDA
