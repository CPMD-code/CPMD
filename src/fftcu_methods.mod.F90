MODULE fftcu_methods

  USE cp_cufft_types,                  ONLY: cp_cufft_plans_t
  USE cp_cufft_utils,                  ONLY: cp_cufft_get_plan
  USE cublas_types,                    ONLY: cublas_handle_t
  USE cuda_types,                      ONLY: cuda_memory_t,&
                                             cuda_stream_t
  USE cuda_utils,                      ONLY: cuda_memcpy_async_device_to_host,&
                                             cuda_memcpy_async_host_to_device,&
                                             cuda_stream_synchronize
  USE cufft_types,                     ONLY: cufft_plan_t
  USE cuuser_utils,                    ONLY: CuUser_getz,&
                                             CuUser_pack_x2y,&
                                             CuUser_pack_y2x,&
                                             CuUser_phasen,&
                                             CuUser_putz,&
                                             CuUser_unpack_x2y,&
                                             CuUser_unpack_y2x
  USE fft,                             ONLY: &
       lfrm, lmsq, lr1, lr1m, lr1s, lr2s, lr3s, lrxpl, lsrm, mfrays, msqf, &
       msqs, msrays, qr1, qr1s, qr2s, qr3max, qr3min, qr3s, sp8, sp9
  USE fft_maxfft,                      ONLY: maxfftn
  USE fftutil_utils,                   ONLY: pack_y2x,&
                                             unpack_x2y
  USE kinds,                           ONLY: real_8
  USE mltfft_utils,                    ONLY: mltfft_cuda

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: fftcu_inv_sprs_1
  PUBLIC :: fftcu_inv_sprs_2
  PUBLIC :: fftcu_inv_full_1
  PUBLIC :: fftcu_inv_full_2
  PUBLIC :: fftcu_frw_sprs_1
  PUBLIC :: fftcu_frw_sprs_2
  PUBLIC :: fftcu_frw_full_1
  PUBLIC :: fftcu_frw_full_2

  LOGICAL, PARAMETER :: use_cpu_unpack_x2y = .TRUE.
  LOGICAL, PARAMETER :: use_cpu_pack_y2x = .TRUE.

CONTAINS


  SUBROUTINE fftcu_inv_sprs_1 ( f, xf, yf, nproc, tr4a2a, t1_d, t2_d, sp5_d, lrxpl_d, &
       & copy_to_device, plans_d, blas_handle, stream )
    COMPLEX(real_8)                          :: f(:), xf(:), yf(:)
    INTEGER, INTENT(IN)                      :: nproc
    LOGICAL, INTENT(IN)                      :: tr4a2a
    TYPE(cuda_memory_t), INTENT(IN)          :: t1_d, t2_d, sp5_d, lrxpl_d
    LOGICAL, INTENT(IN)                      :: copy_to_device
    TYPE(cp_cufft_plans_t), INTENT(INOUT)    :: plans_d
    TYPE(cublas_handle_t), INTENT(IN)        :: blas_handle
    TYPE(cuda_stream_t), INTENT(IN)          :: stream

    INTEGER                                  :: isign, lda, m, mm
    REAL(real_8)                             :: scale
    TYPE(cufft_plan_t), POINTER              :: plan_p

    NULLIFY( plan_p )
    isign = -1
    scale=1._real_8

    m=msrays
    lda=lsrm*lr1m
    mm=qr2s*(qr3max-qr3min+1)

    IF( copy_to_device ) CALL cuda_memcpy_async_host_to_device ( f(1:qr1s*m), t1_d, stream )

    CALL cp_cufft_get_plan ( 'N', 'T', qr1s, m, lr1s, m, plan_p, plans_d, stream )

    CALL mltfft_cuda('N','T',t1_d,qr1s,m,t2_d,m,qr1s,lr1s,m,isign,scale, plan_p, blas_handle, stream )

    CALL CuUser_pack_x2y(t2_d,t1_d,msrays,lda,lrxpl_d,sp5_d,maxfftn,nproc,tr4a2a, stream )

    IF( tr4a2a ) THEN
       CALL cuda_memcpy_async_device_to_host ( t1_d, yf(1:(maxfftn+1)/2), stream )
    ELSE
       CALL cuda_memcpy_async_device_to_host ( t1_d, yf(1:maxfftn), stream )
    ENDIF

  END SUBROUTINE fftcu_inv_sprs_1

  SUBROUTINE fftcu_inv_sprs_2 ( f, xf, yf, nproc, tr4a2a, t1_d, t2_d, sp9_d, msqs_d, &
       & copy_to_host, plans_d, blas_handle, stream )
    COMPLEX(real_8)                          :: f(:), xf(:), yf(:)
    INTEGER, INTENT(IN)                      :: nproc
    LOGICAL, INTENT(IN)                      :: tr4a2a
    TYPE(cuda_memory_t), INTENT(IN)          :: t1_d, t2_d, sp9_d, msqs_d
    LOGICAL, INTENT(IN)                      :: copy_to_host
    TYPE(cp_cufft_plans_t), INTENT(INOUT)    :: plans_d
    TYPE(cublas_handle_t), INTENT(IN)        :: blas_handle
    TYPE(cuda_stream_t), INTENT(IN)          :: stream

    INTEGER                                  :: isign, lda, m, mm
    REAL(real_8)                             :: scale
    TYPE(cufft_plan_t), POINTER              :: plan_p

    NULLIFY( plan_p )
    isign = -1
    scale=1._real_8

    lda=lsrm*lr1m
    mm=qr2s*(qr3max-qr3min+1)

    IF( use_cpu_unpack_x2y ) THEN
       CALL unpack_x2y(xf,yf,mm,lr1,lda,msqs,lmsq,sp9,maxfftn,nproc,tr4a2a)
       CALL cuda_memcpy_async_host_to_device ( yf(1:maxfftn), t2_d, stream )
    ELSE
       IF( tr4a2a ) THEN
          CALL cuda_memcpy_async_host_to_device ( xf(1:(maxfftn+1)/2), t1_d, stream )
       ELSE
          CALL cuda_memcpy_async_host_to_device ( xf(1:maxfftn), t1_d, stream )
       ENDIF

       CALL CuUser_unpack_x2y(t1_d,t2_d,mm,lr1,lda,msqs_d,lmsq,sp9_d,maxfftn,nproc,tr4a2a, stream )
    ENDIF

    m=(qr3max-qr3min+1)*qr1
    CALL cp_cufft_get_plan ( 'N', 'T', qr2s, m, lr2s, m, plan_p, plans_d, stream )
    CALL mltfft_cuda('N','T',t2_d,qr2s,m,t1_d,m,qr2s,lr2s,m,isign,scale, plan_p, blas_handle, stream )
    m=qr1*qr2s
    CALL CuUser_putz(t1_d,t2_d,qr3min,qr3max,qr3s,m, stream )
    CALL cp_cufft_get_plan ( 'N', 'T', qr3s, m, lr3s, m, plan_p, plans_d, stream )
    CALL mltfft_cuda('N','T',t2_d,qr3s,m,t1_d,m,qr3s,lr3s,m,isign,scale, plan_p, blas_handle, stream )

    IF( copy_to_host ) CALL cuda_memcpy_async_device_to_host ( t1_d, f(1:m*qr3s), stream )

  END SUBROUTINE fftcu_inv_sprs_2

  SUBROUTINE fftcu_inv_full_1 ( f, xf, yf, nproc, tr4a2a, t1_d, t2_d, sp5_d, lrxpl_d, &
       & copy_to_device, plans_d, blas_handle, stream )
    COMPLEX(real_8)                          :: f(:), xf(:), yf(:)
    INTEGER, INTENT(IN)                      :: nproc
    LOGICAL, INTENT(IN)                      :: tr4a2a
    TYPE(cuda_memory_t), INTENT(IN)          :: t1_d, t2_d, sp5_d, lrxpl_d
    LOGICAL, INTENT(IN)                      :: copy_to_device
    TYPE(cp_cufft_plans_t), INTENT(INOUT)    :: plans_d
    TYPE(cublas_handle_t), INTENT(IN)        :: blas_handle
    TYPE(cuda_stream_t), INTENT(IN)          :: stream

    INTEGER                                  :: isign, lda, m, mm
    REAL(real_8)                             :: scale
    TYPE(cufft_plan_t), POINTER              :: plan_p

    NULLIFY( plan_p )
    isign = -1
    scale=1._real_8

    m=mfrays
    lda=lfrm*lr1m
    mm=qr2s*qr3s

    IF( copy_to_device ) CALL cuda_memcpy_async_host_to_device ( f(1:qr1s*m), t1_d, stream )

    CALL cp_cufft_get_plan ( 'N', 'T', qr1s, m, lr1s, m, plan_p, plans_d, stream )
    CALL mltfft_cuda('N','T',t1_d,qr1s,m,t2_d,m,qr1s,lr1s,m,isign,scale, plan_p, blas_handle, stream )
    CALL CuUser_pack_x2y(t2_d,t1_d,mfrays,lda,lrxpl_d,sp5_d,maxfftn,nproc,tr4a2a, stream )

    IF( tr4a2a ) THEN
       CALL cuda_memcpy_async_device_to_host ( t1_d, yf(1:(maxfftn+1)/2), stream )
    ELSE
       CALL cuda_memcpy_async_device_to_host ( t1_d, yf(1:maxfftn), stream )
    ENDIF

  END SUBROUTINE fftcu_inv_full_1

  SUBROUTINE fftcu_inv_full_2 ( f, xf, yf, mepos, nproc, tr4a2a, t1_d, t2_d, sp8_d, msqf_d, &
       & copy_to_host, plans_d, blas_handle, stream )
    COMPLEX(real_8)                          :: f(:), xf(:), yf(:)
    INTEGER, INTENT(IN)                      :: mepos, nproc
    LOGICAL, INTENT(IN)                      :: tr4a2a
    TYPE(cuda_memory_t), INTENT(IN)          :: t1_d, t2_d, sp8_d, msqf_d
    LOGICAL, INTENT(IN)                      :: copy_to_host
    TYPE(cp_cufft_plans_t), INTENT(INOUT)    :: plans_d
    TYPE(cublas_handle_t), INTENT(IN)        :: blas_handle
    TYPE(cuda_stream_t), INTENT(IN)          :: stream

    INTEGER                                  :: isign, lda, m, mm, n1o, n1u
    REAL(real_8)                             :: scale
    TYPE(cufft_plan_t), POINTER              :: plan_p

    NULLIFY( plan_p )
    isign = -1
    scale=1._real_8

    lda=lfrm*lr1m
    mm=qr2s*qr3s

    IF( use_cpu_unpack_x2y ) THEN
       CALL unpack_x2y(xf,yf,mm,lr1,lda,msqf,lmsq,sp8,maxfftn,nproc,tr4a2a)
       CALL cuda_memcpy_async_host_to_device ( yf(1:maxfftn), t1_d, stream )
    ELSE
       IF( tr4a2a ) THEN
          CALL cuda_memcpy_async_host_to_device ( xf(1:(maxfftn+1)/2), t2_d, stream )
       ELSE
          CALL cuda_memcpy_async_host_to_device ( xf(1:maxfftn), t2_d, stream )
       ENDIF

       CALL CuUser_unpack_x2y ( t2_d,t1_d,mm,lr1,lda,msqf_d,lmsq,sp8_d,maxfftn,nproc,tr4a2a, stream )
    ENDIF

    m=qr1*qr3s
    CALL cp_cufft_get_plan ( 'N', 'T', qr2s, m, lr2s, m, plan_p, plans_d, stream )
    CALL mltfft_cuda('N','T',t1_d,qr2s,m,t2_d,m,qr2s,lr2s,m,isign,scale, plan_p, blas_handle, stream )

    m=qr1*qr2s
    CALL cp_cufft_get_plan ( 'N', 'T', qr3s, m, lr3s, m, plan_p, plans_d, stream )
    CALL mltfft_cuda('N','T',t2_d,qr3s,m,t1_d,m,qr3s,lr3s,m,isign,scale, plan_p, blas_handle, stream )

    n1u=lrxpl(mepos,1)
    n1o=lrxpl(mepos,2)
    CALL CuUser_phasen (t1_d,qr1,qr2s,qr3s,n1u,n1o,lr2s,lr3s, stream )

    IF( copy_to_host ) CALL cuda_memcpy_async_device_to_host ( t1_d, f(1:qr3s*m), stream )

  END SUBROUTINE fftcu_inv_full_2


  SUBROUTINE fftcu_frw_sprs_1 ( f, xf, yf, nproc, tr4a2a, t1_d, t2_d, sp9_d, msqs_d, &
       & copy_to_device, plans_d, blas_handle, stream )
    COMPLEX(real_8)                          :: f(:), xf(:), yf(:)
    INTEGER, INTENT(IN)                      :: nproc
    LOGICAL, INTENT(IN)                      :: tr4a2a
    TYPE(cuda_memory_t), INTENT(IN)          :: t1_d, t2_d, sp9_d, msqs_d
    LOGICAL, INTENT(IN)                      :: copy_to_device
    TYPE(cp_cufft_plans_t), INTENT(INOUT)    :: plans_d
    TYPE(cublas_handle_t), INTENT(IN)        :: blas_handle
    TYPE(cuda_stream_t), INTENT(IN)          :: stream

    INTEGER                                  :: isign, lda, m, mm
    REAL(real_8)                             :: scale
    TYPE(cufft_plan_t), POINTER              :: plan_p

    NULLIFY( plan_p )
    isign = 1
    scale=1._real_8

    m=qr1*qr2s

    IF( copy_to_device ) CALL cuda_memcpy_async_host_to_device ( f(1:m*qr3s), t1_d, stream )

    CALL cp_cufft_get_plan ( 'T', 'N', m, qr3s, lr3s, m, plan_p, plans_d, stream )
    CALL mltfft_cuda('T','N',t1_d,m,qr3s,t2_d,qr3s,m,lr3s,m,isign,scale, plan_p, blas_handle, stream )
    CALL CuUser_getz(t2_d,t1_d,qr3min,qr3max,qr3s,m, stream )
    m=(qr3max-qr3min+1)*qr1
    CALL cp_cufft_get_plan ( 'T', 'N', m, qr2s, lr2s, m, plan_p, plans_d, stream )
    CALL mltfft_cuda('T','N',t1_d,m,qr2s,t2_d,qr2s,m,lr2s,m,isign,scale, plan_p, blas_handle, stream )
    lda=lsrm*lr1m
    mm=qr2s*(qr3max-qr3min+1)

    IF( use_cpu_pack_y2x ) THEN
       CALL cuda_memcpy_async_device_to_host ( t2_d, yf(1:maxfftn), stream )
       CALL cuda_stream_synchronize( stream )
       CALL pack_y2x(xf,yf,mm,lr1,lda,msqs,lmsq,sp9,maxfftn,nproc,tr4a2a)
    ELSE
       CALL CuUser_pack_y2x(t1_d,t2_d,mm,lr1,lda,msqs_d,lmsq,sp9_d,maxfftn,nproc,tr4a2a, stream )

       IF( tr4a2a ) THEN
          CALL cuda_memcpy_async_device_to_host ( t1_d, xf(1:(maxfftn+1)/2), stream )
       ELSE
          CALL cuda_memcpy_async_device_to_host ( t1_d, xf(1:maxfftn), stream )
       ENDIF
    ENDIF

  END SUBROUTINE fftcu_frw_sprs_1

  SUBROUTINE fftcu_frw_sprs_2 ( f, xf, yf, nproc, tr4a2a, t1_d, t2_d, sp5_d, lrxpl_d, &
       & copy_to_host, plans_d, blas_handle, stream )
    COMPLEX(real_8)                          :: f(:), xf(:), yf(:)
    INTEGER, INTENT(IN)                      :: nproc
    LOGICAL, INTENT(IN)                      :: tr4a2a
    TYPE(cuda_memory_t), INTENT(IN)          :: t1_d, t2_d, sp5_d, lrxpl_d
    LOGICAL, INTENT(IN)                      :: copy_to_host
    TYPE(cp_cufft_plans_t), INTENT(INOUT)    :: plans_d
    TYPE(cublas_handle_t), INTENT(IN)        :: blas_handle
    TYPE(cuda_stream_t), INTENT(IN)          :: stream

    INTEGER                                  :: isign, lda, m, mm
    REAL(real_8)                             :: scale
    TYPE(cufft_plan_t), POINTER              :: plan_p

    NULLIFY( plan_p )
    isign = 1
    scale=1._real_8/( REAL(lr1s,kind=real_8)*REAL(lr2s,kind=real_8)*REAL(lr3s,kind=real_8) )

    lda=lsrm*lr1m
    mm=qr2s*(qr3max-qr3min+1)

    IF( tr4a2a ) THEN
       CALL cuda_memcpy_async_host_to_device ( yf(1:(maxfftn+1)/2), t1_d, stream )
    ELSE
       CALL cuda_memcpy_async_host_to_device ( yf(1:maxfftn), t1_d, stream )
    ENDIF

    CALL CuUser_unpack_y2x(t2_d,t1_d,mm,msrays,lda,lrxpl_d,sp5_d,maxfftn,nproc,tr4a2a, stream )
    m=msrays
    CALL cp_cufft_get_plan ( 'T', 'N', m, qr1s, lr1s, m, plan_p, plans_d, stream )
    CALL mltfft_cuda('T','N',t2_d,m,qr1s,t1_d,qr1s,m,lr1s,m,isign,scale, plan_p, blas_handle, stream )

    IF( copy_to_host ) CALL cuda_memcpy_async_device_to_host ( t1_d, f(1:qr1s*m), stream )

  END SUBROUTINE fftcu_frw_sprs_2

  SUBROUTINE fftcu_frw_full_1 ( f, xf, yf, mepos, nproc, tr4a2a, t1_d, t2_d, sp8_d, msqf_d, &
       & copy_to_device, plans_d, blas_handle, stream )
    COMPLEX(real_8)                          :: f(:), xf(:), yf(:)
    INTEGER, INTENT(IN)                      :: mepos, nproc
    LOGICAL, INTENT(IN)                      :: tr4a2a
    TYPE(cuda_memory_t), INTENT(IN)          :: t1_d, t2_d, sp8_d, msqf_d
    LOGICAL, INTENT(IN)                      :: copy_to_device
    TYPE(cp_cufft_plans_t), INTENT(INOUT)    :: plans_d
    TYPE(cublas_handle_t), INTENT(IN)        :: blas_handle
    TYPE(cuda_stream_t), INTENT(IN)          :: stream

    INTEGER                                  :: isign, lda, m, mm, n1o, n1u
    REAL(real_8)                             :: scale
    TYPE(cufft_plan_t), POINTER              :: plan_p

    NULLIFY( plan_p )
    isign = 1
    scale=1._real_8

    m=qr1*qr2s

    IF( copy_to_device ) CALL cuda_memcpy_async_host_to_device ( f(1:m*qr3s), t1_d, stream )

    n1u=lrxpl(mepos,1)
    n1o=lrxpl(mepos,2)
    CALL CuUser_phasen(t1_d,qr1,qr2s,qr3s,n1u,n1o,lr2s,lr3s, stream )

    CALL cp_cufft_get_plan ( 'T', 'N', m, qr3s, lr3s, m, plan_p, plans_d, stream )
    CALL mltfft_cuda('T','N',t1_d,m,qr3s,t2_d,qr3s,m,lr3s,m,isign,scale, plan_p, blas_handle, stream )
    m=qr1*qr3s

    CALL cp_cufft_get_plan ( 'T', 'N', m, qr2s, lr2s, m, plan_p, plans_d, stream )
    CALL mltfft_cuda('T','N',t2_d,m,qr2s,t1_d,qr2s,m,lr2s,m,isign,scale, plan_p, blas_handle, stream )
    lda=lfrm*lr1m
    mm=qr2s*qr3s

    IF( use_cpu_pack_y2x ) THEN
       CALL cuda_memcpy_async_device_to_host ( t1_d, yf(1:maxfftn), stream )
       CALL cuda_stream_synchronize( stream )
       CALL pack_y2x(xf,yf,mm,lr1,lda,msqf,lmsq,sp8,maxfftn,nproc,tr4a2a)
    ELSE
       CALL CuUser_pack_y2x(t2_d,t1_d,mm,lr1,lda,msqf_d,lmsq,sp8_d,maxfftn,nproc,tr4a2a, stream )

       IF( tr4a2a ) THEN
          CALL cuda_memcpy_async_device_to_host ( t2_d, xf(1:(maxfftn+1)/2), stream )
       ELSE
          CALL cuda_memcpy_async_device_to_host ( t2_d, xf(1:maxfftn), stream )
       ENDIF
    ENDIF

  END SUBROUTINE fftcu_frw_full_1

  SUBROUTINE fftcu_frw_full_2 ( f, xf, yf, nproc, tr4a2a, t1_d, t2_d, sp5_d, lrxpl_d, &
       & copy_to_host, plans_d, blas_handle, stream )
    COMPLEX(real_8)                          :: f(:), xf(:), yf(:)
    INTEGER, INTENT(IN)                      :: nproc
    LOGICAL, INTENT(IN)                      :: tr4a2a
    TYPE(cuda_memory_t), INTENT(IN)          :: t1_d, t2_d, sp5_d, lrxpl_d
    LOGICAL, INTENT(IN)                      :: copy_to_host
    TYPE(cp_cufft_plans_t), INTENT(INOUT)    :: plans_d
    TYPE(cublas_handle_t), INTENT(IN)        :: blas_handle
    TYPE(cuda_stream_t), INTENT(IN)          :: stream

    INTEGER                                  :: isign, lda, m, mm
    REAL(real_8)                             :: scale
    TYPE(cufft_plan_t), POINTER              :: plan_p

    NULLIFY( plan_p )
    isign = 1
    scale=1._real_8/( REAL(lr1s,kind=real_8)*REAL(lr2s,kind=real_8)*REAL(lr3s,kind=real_8) )

    lda=lfrm*lr1m
    mm=qr2s*qr3s
    m=mfrays

    IF( tr4a2a ) THEN
       CALL cuda_memcpy_async_host_to_device ( yf(1:(maxfftn+1)/2), t1_d, stream )
    ELSE
       CALL cuda_memcpy_async_host_to_device ( yf(1:maxfftn), t1_d, stream )
    ENDIF

    CALL CuUser_unpack_y2x(t2_d,t1_d,mm,mfrays,lda,lrxpl_d,sp5_d,maxfftn,nproc,tr4a2a, stream )
    CALL cp_cufft_get_plan ( 'T', 'N', m, qr1s, lr1s, m, plan_p, plans_d, stream )
    CALL mltfft_cuda('T','N',t2_d,m,qr1s,t1_d,qr1s,m,lr1s,m,isign,scale, plan_p, blas_handle, stream )

    IF( copy_to_host ) CALL cuda_memcpy_async_device_to_host ( t1_d, f(1:qr1s*m), stream )

  END SUBROUTINE fftcu_frw_full_2


END MODULE fftcu_methods
