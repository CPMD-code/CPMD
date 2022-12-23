#include "cuuser_utils.h"

#if defined(__CUDA)

#ifndef __CUUSER_C_UTILS_KERNELS___
#define __CUUSER_C_UTILS_KERNELS___

__global__ void CuUser_Kernel_SetBlock2Zero (cuDoubleComplex *a, int N, int M, int LDBX, int LDBY )
{
    const int row = threadIdx.x + blockDim.x * blockIdx.x;
    const int col = threadIdx.y + blockDim.y * blockIdx.y;

    const int id = col*LDBX+row;

    if( col >= M && col < LDBY && row < LDBX)
    {
        a[id] = make_cuDoubleComplex(0.0, 0.0);

        //a[id].x = 0.0;
        //a[id].y = 0.0;
    }

    if( row >= N && row < LDBX && col < M)
    {
        a[id] = make_cuDoubleComplex(0.0, 0.0);

        //a[id].x = 0.0;
        //a[id].y = 0.0;
    }
}

__global__ void CuUser_Kernel_SetBlock2ZeroScale (cuDoubleComplex *a, double scale, int N, int M, int LDBX, int LDBY )
{
    const int row = threadIdx.x + blockDim.x * blockIdx.x;
    const int col = threadIdx.y + blockDim.y * blockIdx.y;

    const int id = col*LDBX+row;

    if( col >= M && col < LDBY && row < LDBX)
    {
        a[id] = make_cuDoubleComplex(0.0, 0.0);

//        a[id].x = 0.0;
//        a[id].y = 0.0;
    }

    if( row >= N && row < LDBX && col < M)
    {
        a[id] = make_cuDoubleComplex(0.0, 0.0);

//        a[id].x = 0.0;
//        a[id].y = 0.0;
    }

    if( col < LDBY && row < LDBX)
    {
        a[id].x = a[id].x * scale;
        a[id].y = a[id].y * scale;
    }
}

__global__ void CuUser_Kernel_Identity_double( double *devMatrix, int m )
{
    const int row = threadIdx.x + blockDim.x * blockIdx.x;
    const int col = threadIdx.y + blockDim.y * blockIdx.y;

    const int id = col*m+row;

    if(row < m && col < m)
    {
        if(row == col)
            devMatrix[id] = 1.0;
        else
            devMatrix[id] = 0.0;
    }
}

__global__ void CuUser_Kernel_MatMov ( int n, int m, cuDoubleComplex *a, int LDA, cuDoubleComplex *b, int LDB )
{
    const int row = threadIdx.x + blockDim.x * blockIdx.x;
    const int col = threadIdx.y + blockDim.y * blockIdx.y;

    const int ida = col*LDA+row;
    const int idb = col*LDB+row;

    if( col < m && row < n)
        b[idb] = a[ida];
}

__global__ void CuUser_Kernel_Zeroing ( cuDoubleComplex *a, int size)
{
    const int id = threadIdx.x + blockDim.x * blockIdx.x;

    if (id < size)
    {
        a[id] = make_cuDoubleComplex(0.0, 0.0);
//        a[id].x = 0.0;
//        a[id].y = 0.0;
    }
}

__global__ void CuUser_Kernel_PhaseN ( cuDoubleComplex *f, int kr1, int kr2s, int kr3s, int n1u, int n1o, int nr2s, int nr3s )
{
    const int i   = threadIdx.x + blockDim.x * blockIdx.x;
    const int j   = threadIdx.y + blockDim.y * blockIdx.y;
    const int k   = threadIdx.z + blockDim.z * blockIdx.z;

    const int imax = n1o - n1u + 1;
    const double pf[2] = {1.0, -1.0};

    if( i < imax && j < nr2s && k < nr3s )
    {
        const int ijk  = ( k + j + i + n1u + 3 ) % 2;
        const int iijk = k * kr2s * kr1 + j * kr1 + i;

        f[iijk].x = f[iijk].x * pf[ijk];
        f[iijk].y = f[iijk].y * pf[ijk];
    }
}

__global__ void CuUser_Kernel_Unpack_x2y_8 ( cuDoubleComplex *xf, cuDoubleComplex *yf, int m, int lr1, int lda, int *msp, int lmsp, int *sp8, int maxfft, int mproc )
{
    const int ip  = threadIdx.x + blockDim.x * blockIdx.x;
    const int k   = threadIdx.y + blockDim.y * blockIdx.y;
    const int i   = threadIdx.z + blockDim.z * blockIdx.z;

    const int mxrp = sp8[ip];

    if( ip < mproc && k < mxrp && i < lr1)
    {
        const int idx = k + ip * lda + i * mxrp;
        const int idy = i * m + msp[k + ip * lmsp] - 1;

        yf[idy] = xf[idx];
    }
}

__global__ void CuUser_Kernel_Unpack_x2y_4 ( cuComplex *xf, cuDoubleComplex *yf, int m, int lr1, int lda, int *msp, int lmsp, int *sp8, int maxfft, int mproc )
{
    const int ip  = threadIdx.x + blockDim.x * blockIdx.x;
    const int k   = threadIdx.y + blockDim.y * blockIdx.y;
    const int i   = threadIdx.z + blockDim.z * blockIdx.z;

    const int mxrp = sp8[ip];

    if( ip < mproc && k < mxrp && i < lr1)
    {
        const int idx = k + ip * lda + i * mxrp;
        const int idy = i * m + msp[k + ip * lmsp] - 1;

//        yf[idy] = static_cast<cuDoubleComplex> (xf[idx]);

        yf[idy].x = xf[idx].x;
        yf[idy].y = xf[idx].y;
    }
}

__global__ void CuUser_Kernel_Unpack_y2x_8 ( cuDoubleComplex *xf, cuDoubleComplex *yf, int m, int nrays, int lda, int *jrxpl, int *sp5, int maxfft, int mproc )
{
    const int ip  = threadIdx.x + blockDim.x * blockIdx.x;
    const int k   = threadIdx.y + blockDim.y * blockIdx.y;

    const int nrx = nrays *  sp5[ip];
    const int nrs = nrays * (jrxpl[ip] - 1);

    if( ip < mproc && k < nrx)
    {
        const int idx = k + nrs;
        const int idy = k + ip * lda;

        xf[idx] = yf[idy];
    }
}

__global__ void CuUser_Kernel_Unpack_y2x_4 ( cuDoubleComplex *xf, cuComplex *yf, int m, int nrays, int lda, int *jrxpl, int *sp5, int maxfft, int mproc )
{
    const int ip  = threadIdx.x + blockDim.x * blockIdx.x;
    const int k   = threadIdx.y + blockDim.y * blockIdx.y;

    const int nrx = nrays *  sp5[ip];
    const int nrs = nrays * (jrxpl[ip] - 1);

    if( ip < mproc && k < nrx)
    {
        const int idx = k + nrs;
        const int idy = k + ip * lda;

//        xf[idx] = yf[idy];

        xf[idx].x = yf[idy].x;
        xf[idx].y = yf[idy].y;
    }
}

__global__ void CuUser_Kernel_Pack_x2y_8 ( cuDoubleComplex *xf, cuDoubleComplex *yf, int nrays, int lda, int *jrxpl, int *sp5, int maxfft, int mproc )
{
    const int ip  = threadIdx.x + blockDim.x * blockIdx.x;
    const int k   = threadIdx.y + blockDim.y * blockIdx.y;

    const int nrx = nrays *  sp5[ip];
    const int nrs = nrays * (jrxpl[ip] - 1);

    if( ip < mproc && k < nrx)
    {
        const int idx = k + nrs;
        const int idy = k + ip * lda;

        yf[idy] = xf[idx];
    }
}

__global__ void CuUser_Kernel_Pack_x2y_4 ( cuDoubleComplex *xf, cuComplex *yf, int nrays, int lda, int *jrxpl, int *sp5, int maxfft, int mproc )
{
    const int ip  = threadIdx.x + blockDim.x * blockIdx.x;
    const int k   = threadIdx.y + blockDim.y * blockIdx.y;

    const int nrx = nrays *  sp5[ip];
    const int nrs = nrays * (jrxpl[ip] - 1);

    if( ip < mproc && k < nrx)
    {
        const int idx = k + nrs;
        const int idy = k + ip * lda;

//        yf[idy] = xf[idx];

        yf[idy].x = xf[idx].x;
        yf[idy].y = xf[idx].y;
    }
}

__global__ void CuUser_Kernel_Pack_y2x_8 ( cuDoubleComplex *xf, cuDoubleComplex *yf, int m, int lr1, int lda, int *msp, int lmsp, int *sp8, int maxfft, int mproc )
{
    const int ip  = threadIdx.x + blockDim.x * blockIdx.x;
    const int k   = threadIdx.y + blockDim.y * blockIdx.y;
    const int i   = threadIdx.z + blockDim.z * blockIdx.z;

    const int mxrp = sp8[ip];

    if( ip < mproc && k < mxrp && i < lr1)
    {
        const int idx = k + ip * lda + i * mxrp;
        const int idy = i * m + msp[k + ip * lmsp] - 1;

        xf[idx] = yf[idy];
    }
}

__global__ void CuUser_Kernel_Pack_y2x_4 ( cuComplex *xf, cuDoubleComplex *yf, int m, int lr1, int lda, int *msp, int lmsp, int *sp8, int maxfft, int mproc )
{
    const int ip  = threadIdx.x + blockDim.x * blockIdx.x;
    const int k   = threadIdx.y + blockDim.y * blockIdx.y;
    const int i   = threadIdx.z + blockDim.z * blockIdx.z;

    const int mxrp = sp8[ip];

    if( ip < mproc && k < mxrp && i < lr1)
    {
        const int idx = k + ip * lda + i * mxrp;
        const int idy = i * m + msp[k + ip * lmsp] - 1;

//        xf[idx] = yf[idy];

        xf[idx].x = yf[idy].x;
        xf[idx].y = yf[idy].y;
    }
}

__global__ void CuUser_Kernel_Build_Density_Sum ( double alpha_real, double alpha_imag, cuDoubleComplex *psi, double *rho, int n )
{
    const int idx = threadIdx.x + blockDim.x * blockIdx.x;

    if( idx < n )
    {
        rho[idx] += alpha_real * psi[idx].x * psi[idx].x + alpha_imag * psi[idx].y * psi[idx].y;
    }
}

__global__ void CuUser_Kernel_Pointwise_CxR ( cuDoubleComplex *xf_p, double *yf_p, int n )
{
    const int idx = threadIdx.x + blockDim.x * blockIdx.x;

    if( idx < n )
    {
        xf_p[idx].x *= yf_p[idx];
        xf_p[idx].y *= yf_p[idx];
    }
}

__global__ void CuUser_Kernel_Set_Psi_1_Stage_G ( double alpha_real, double alpha_imag, cuDoubleComplex * c1_p, cuDoubleComplex * psi_p, int jgw, int *nzfs_p, int *inzs_p, bool geq0 )
{
    const int idx = threadIdx.x + blockDim.x * blockIdx.x;

    // TODO acm: for alpha = (1, 0) we can have a special case, with less operations.
    if( idx < jgw )
    {
        psi_p[nzfs_p[idx]-1].x =  alpha_real * c1_p[idx].x - alpha_imag * c1_p[idx].y;
        psi_p[nzfs_p[idx]-1].y =  alpha_real * c1_p[idx].y + alpha_imag * c1_p[idx].x;

        psi_p[inzs_p[idx]-1].x =  alpha_real * c1_p[idx].x + alpha_imag * c1_p[idx].y;
        psi_p[inzs_p[idx]-1].y = -alpha_real * c1_p[idx].y + alpha_imag * c1_p[idx].x;
    }

    if( idx == 0 )
        if( geq0 )
        {
            psi_p[nzfs_p[0]-1].x = alpha_real * c1_p[0].x - alpha_imag * c1_p[0].y;
            psi_p[nzfs_p[0]-1].y = alpha_real * c1_p[0].y + alpha_imag * c1_p[0].x;
        }
}

__global__ void CuUser_Kernel_Set_Psi_2_Stages_G ( cuDoubleComplex * c1_p, cuDoubleComplex * c2_p, cuDoubleComplex * psi_p, int jgw, int *nzfs_p, int *inzs_p, bool geq0 )
{
    const int idx = threadIdx.x + blockDim.x * blockIdx.x;

    if( idx < jgw )
    {
        psi_p[nzfs_p[idx]-1].x =  c1_p[idx].x - c2_p[idx].y;
        psi_p[nzfs_p[idx]-1].y =  c1_p[idx].y + c2_p[idx].x;

        psi_p[inzs_p[idx]-1].x =  c1_p[idx].x + c2_p[idx].y;
        psi_p[inzs_p[idx]-1].y = -c1_p[idx].y + c2_p[idx].x;
    }

    if( idx == 0 )
        if( geq0 )
        {
            psi_p[nzfs_p[0]-1].x = c1_p[0].x - c2_p[0].y;
            psi_p[nzfs_p[0]-1].y = c1_p[0].y + c2_p[0].x;
        }
}

__global__ void CuUser_Kernel_Eicalc ( cuDoubleComplex * eivps_p, cuDoubleComplex * eirop_p, int nhg, int nat, int *iatpt_p, int *inyh_p, cuDoubleComplex *ei1_p, cuDoubleComplex *ei2_p, cuDoubleComplex *ei3_p, double *vps_p, double *rhops_p )
{
    const int idx = threadIdx.x + blockDim.x * blockIdx.x;
    const int idy = threadIdx.y + blockDim.y * blockIdx.y;

    if( idx < nhg  && idy < nat )
    {
        const int inyh0 = inyh_p[idx];
        const int inyh1 = inyh_p[idx+nhg];
        const int inyh2 = inyh_p[idx+nhg*2];

        const cuDoubleComplex ei123 = cuCmul(cuCmul(ei1_p[inyh0+idy*nat], ei2_p[inyh1+idy*nat]), ei3_p[inyh2+idy*nat]); // acm: idy*nat maybe wrong!

        //const int ia = iatpt_p[idy];
        const int is = iatpt_p[idy+nat];

        eivps_p[idx] = cuCadd(eivps_p[idx], make_cuDoubleComplex(ei123.x * vps_p[idx+is*nhg],   ei123.y * vps_p[idx+is*nhg])); // vps_p[idx+is*nhg] maybe wrong
        eirop_p[idx] = cuCadd(eirop_p[idx], make_cuDoubleComplex(ei123.x * rhops_p[idx+is*nhg], ei123.y * rhops_p[idx+is*nhg]));
    }
}

#endif // __CUUSER_C_UTILS_KERNELS___
#endif // __CUDA
