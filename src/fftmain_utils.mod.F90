#include "cpmd_global.h"

MODULE fftmain_utils

  USE cp_cuda_types,                   ONLY: cp_cuda_env
  USE cp_cufft_types,                  ONLY: cp_cufft,&
                                             cp_cufft_device_get_ptrs,&
                                             cp_cufft_plans_t,&
                                             cp_cufft_stream_get_ptrs
  USE cp_cufft_utils,                  ONLY: cp_cufft_get_plan
  USE cublas_types,                    ONLY: cublas_handle_t
  USE cuda_types,                      ONLY: cuda_memory_t,&
                                             cuda_stream_t
  USE cuda_utils,                      ONLY: cuda_memcpy_async_device_to_host,&
                                             cuda_memcpy_async_host_to_device,&
                                             cuda_stream_synchronize
  USE cufft_types,                     ONLY: cufft_plan_t
  USE error_handling,                  ONLY: stopgm
  USE fft,                             ONLY: &
       lfrm, lmsq, lr1, lr1m, lr1s, lr2s, lr3s, lrxpl, lsrm, mfrays, msqf, &
       msqs, msrays, qr1, qr1s, qr2s, qr3max, qr3min, qr3s, sp5, sp8, sp9, &
       xf, yf
  USE fft_maxfft,                      ONLY: maxfftn
  USE fftcu_methods,                   ONLY: fftcu_frw_full_1,&
                                             fftcu_frw_full_2,&
                                             fftcu_frw_sprs_1,&
                                             fftcu_frw_sprs_2,&
                                             fftcu_inv_full_1,&
                                             fftcu_inv_full_2,&
                                             fftcu_inv_sprs_1,&
                                             fftcu_inv_sprs_2
  USE fftutil_utils,                   ONLY: fft_comm,&
                                             getz,&
                                             pack_x2y,&
                                             pack_y2x,&
                                             phasen,&
                                             putz,&
                                             unpack_x2y,&
                                             unpack_y2x
  USE kinds,                           ONLY: real_8
  USE mltfft_utils,                    ONLY: mltfft_cuda,&
                                             mltfft_default,&
                                             mltfft_essl,&
                                             mltfft_fftw,&
                                             mltfft_hp
  USE parac,                           ONLY: parai
  USE system,                          ONLY: cntl
  USE thread_view_types,               ONLY: thread_view_t
  USE timer,                           ONLY: tihalt,&
                                             tiset

  !$ USE omp_lib, ONLY: omp_in_parallel, omp_get_thread_num


  USE, INTRINSIC :: ISO_C_BINDING,     ONLY: C_PTR, C_NULL_PTR

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mltfft
  PUBLIC :: invfftn
  PUBLIC :: fwfftn
  !public :: fftnew


CONTAINS


  ! ==================================================================
  SUBROUTINE fftnew(isign,f,sparse,comm)
    ! ==--------------------------------------------------------------==
    INTEGER, INTENT(IN)                      :: isign
    COMPLEX(real_8)                          :: f(:)
    LOGICAL, INTENT(IN)                      :: sparse
    INTEGER, INTENT(IN)                      :: comm

    CHARACTER(*), PARAMETER                  :: procedureN = 'fftnew'

    COMPLEX(real_8), DIMENSION(:), &
      POINTER __CONTIGUOUS                   :: xf_ptr, yf_ptr
    INTEGER                                  :: lda, m, mm, n1o, n1u
    REAL(real_8)                             :: scale

    xf_ptr => xf(:,1)
    yf_ptr => yf(:,1)

    LDA=HUGE(0);MM=HUGE(0);N1U=HUGE(0);N1O=HUGE(0);M=HUGE(0)
    scale=HUGE(0.0_real_8)
    IF (isign.EQ.-1) THEN
       scale=1._real_8
       IF (sparse) THEN
          m=msrays
          CALL mltfft('N','T',f,qr1s,m,xf_ptr,m,qr1s,lr1s,m,isign,scale )
          lda=lsrm*lr1m
          mm=qr2s*(qr3max-qr3min+1)
          CALL pack_x2y(xf_ptr,yf_ptr,msrays,lda,lrxpl,sp5,maxfftn,parai%nproc,cntl%tr4a2a)
          CALL fft_comm(yf_ptr,xf_ptr,lda,cntl%tr4a2a,comm)
          CALL unpack_x2y(xf_ptr,yf_ptr,mm,lr1,lda,msqs,lmsq,sp9,maxfftn,parai%nproc,cntl%tr4a2a)
          m=(qr3max-qr3min+1)*qr1
          CALL mltfft('N','T',yf_ptr,qr2s,m,xf_ptr,m,qr2s,lr2s,m,isign,scale )
          m=qr1*qr2s
          CALL putz(xf_ptr,yf_ptr,qr3min,qr3max,qr3s,m)
          CALL mltfft('N','T',yf_ptr,qr3s,m,f,m,qr3s,lr3s,m,isign,scale )
       ELSE
          m=mfrays
          CALL mltfft('N','T',f,qr1s,m,xf_ptr,m,qr1s,lr1s,m,isign,scale )
          lda=lfrm*lr1m
          mm=qr2s*qr3s
          CALL pack_x2y(xf_ptr,yf_ptr,mfrays,lda,lrxpl,sp5,maxfftn,parai%nproc,cntl%tr4a2a)
          CALL fft_comm(yf_ptr,xf_ptr,lda,cntl%tr4a2a,comm)
          CALL unpack_x2y(xf_ptr,yf_ptr,mm,lr1,lda,msqf,lmsq,sp8,maxfftn,parai%nproc,cntl%tr4a2a)
          m=qr1*qr3s
          CALL mltfft('N','T',yf_ptr,qr2s,m,xf_ptr,m,qr2s,lr2s,m,isign,scale )
          m=qr1*qr2s
          CALL mltfft('N','T',xf_ptr,qr3s,m,f,m,qr3s,lr3s,m,isign,scale )
          n1u=lrxpl(parai%mepos,1)
          n1o=lrxpl(parai%mepos,2)
          CALL phasen(f,qr1,qr2s,qr3s,n1u,n1o,lr2s,lr3s)
       ENDIF
    ELSE
       IF (sparse) THEN
          scale=1._real_8
          m=qr1*qr2s
          CALL mltfft('T','N',f,m,qr3s,xf_ptr,qr3s,m,lr3s,m,isign,scale )
          CALL getz(xf_ptr,f,qr3min,qr3max,qr3s,m)
          m=(qr3max-qr3min+1)*qr1
          CALL mltfft('T','N',f,m,qr2s,yf_ptr,qr2s,m,lr2s,m,isign,scale )
          lda=lsrm*lr1m
          mm=qr2s*(qr3max-qr3min+1)
          CALL pack_y2x(xf_ptr,yf_ptr,mm,lr1,lda,msqs,lmsq,sp9,maxfftn,parai%nproc,cntl%tr4a2a)
          CALL fft_comm(xf_ptr,yf_ptr,lda,cntl%tr4a2a,comm)
          CALL unpack_y2x(xf_ptr,yf_ptr,mm,msrays,lda,lrxpl,sp5,maxfftn,parai%nproc,cntl%tr4a2a)
          scale=1._real_8/REAL(lr1s*lr2s*lr3s,kind=real_8)
          m=msrays
          CALL mltfft('T','N',xf_ptr,m,qr1s,f,qr1s,m,lr1s,m,isign,scale )
       ELSE
          scale=1._real_8
          n1u=lrxpl(parai%mepos,1)
          n1o=lrxpl(parai%mepos,2)
          CALL phasen(f,qr1,qr2s,qr3s,n1u,n1o,lr2s,lr3s)
          m=qr1*qr2s
          CALL mltfft('T','N',f,m,qr3s,xf_ptr,qr3s,m,lr3s,m,isign,scale )
          m=qr1*qr3s
          CALL mltfft('T','N',xf_ptr,m,qr2s,yf_ptr,qr2s,m,lr2s,m,isign,scale )
          lda=lfrm*lr1m
          mm=qr2s*qr3s
          CALL pack_y2x(xf_ptr,yf_ptr,mm,lr1,lda,msqf,lmsq,sp8,maxfftn,parai%nproc,cntl%tr4a2a)
          CALL fft_comm(xf_ptr,yf_ptr,lda,cntl%tr4a2a,comm)
          CALL unpack_y2x(xf_ptr,yf_ptr,mm,mfrays,lda,lrxpl,sp5,maxfftn,parai%nproc,cntl%tr4a2a)
          scale=1._real_8/REAL(lr1s*lr2s*lr3s,kind=real_8)
          m=mfrays
          CALL mltfft('T','N',xf_ptr,m,qr1s,f,qr1s,m,lr1s,m,isign,scale )
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
  END SUBROUTINE fftnew


  ! ==================================================================
  SUBROUTINE fftnew_cuda(isign,f,sparse,comm,thread_view, copy_data_to_device, copy_data_to_host )
    ! ==--------------------------------------------------------------==
    INTEGER, INTENT(IN)                      :: isign
    COMPLEX(real_8)                          :: f(:)
    LOGICAL, INTENT(IN)                      :: sparse
    INTEGER, INTENT(IN)                      :: comm
    TYPE(thread_view_t), INTENT(IN), &
      OPTIONAL                               :: thread_view
    LOGICAL, INTENT(IN), OPTIONAL            :: copy_data_to_device, &
                                                copy_data_to_host

    CHARACTER(*), PARAMETER                  :: procedureN = 'fftnew_cuda'

    COMPLEX(real_8), POINTER __CONTIGUOUS    :: xf_ptr(:), yf_ptr(:)
    INTEGER                                  :: device_idx, host_buff_ptr, &
                                                lda, stream_idx
    LOGICAL                                  :: copy_to_device, copy_to_host
    REAL(real_8)                             :: scale
    TYPE(cp_cufft_plans_t), POINTER          :: plans_d
    TYPE(cublas_handle_t), POINTER           :: blas_handle_p
    TYPE(cuda_memory_t), POINTER             :: lrxpl_d, msqf_d, msqs_d, &
                                                sp5_d, sp8_d, sp9_d, t1_d, &
                                                t2_d
    TYPE(cuda_stream_t), POINTER             :: stream_p

    copy_to_host = .TRUE.
    copy_to_device = .TRUE.
    IF( PRESENT( copy_data_to_host ) ) copy_to_host = copy_data_to_host
    IF( PRESENT( copy_data_to_device ) ) copy_to_device = copy_data_to_device

    !vw if thread view is available, the proben the device/stream
    !vw else NEED TO CLEAN THAT

    IF( PRESENT( thread_view ) ) THEN
       device_idx = thread_view%device_idx
       stream_idx = thread_view%stream_idx - 1 !vw FIX that -1
       host_buff_ptr = thread_view%host_buff_ptr
    ELSE
       device_idx = 1
       stream_idx = 0
       !$ stream_idx = omp_get_thread_num()
       host_buff_ptr = stream_idx + 1 !vw FIX that +1
    ENDIF

    xf_ptr => xf(:, host_buff_ptr )
    yf_ptr => yf(:, host_buff_ptr )

    CALL cp_cufft_stream_get_ptrs ( cp_cufft, device_idx, stream_idx, stream=stream_p, &
         t1_d=t1_d, t2_d=t2_d, plans=plans_d, blas_handle=blas_handle_p )
    CALL cp_cufft_device_get_ptrs ( cp_cufft, device_idx, sp5_d=sp5_d, sp8_d=sp8_d, &
         sp9_d=sp9_d, msqs_d=msqs_d, msqf_d=msqf_d, lrxpl_d=lrxpl_d ) 

    LDA=HUGE(0)
    scale=HUGE(0.0_real_8)
    IF (isign.EQ.-1) THEN
       scale=1._real_8
       IF (sparse) THEN
          CALL fftcu_inv_sprs_1 ( f, xf_ptr, yf_ptr, parai%nproc, cntl%tr4a2a,&
               t1_d, t2_d, sp5_d, lrxpl_d, copy_to_device, plans_d,&
               blas_handle_p, stream_p )

          CALL cuda_stream_synchronize(stream_p)
          lda=lsrm*lr1m
          CALL fft_comm(yf_ptr,xf_ptr,lda,cntl%tr4a2a,comm)

          CALL fftcu_inv_sprs_2 ( f, xf_ptr, yf_ptr, parai%nproc, cntl%tr4a2a,&
               t1_d, t2_d, sp9_d, msqs_d, copy_to_host, plans_d,&
               blas_handle_p, stream_p )

       ELSE
          CALL fftcu_inv_full_1 ( f, xf_ptr, yf_ptr, parai%nproc, cntl%tr4a2a,&
               t1_d, t2_d, sp5_d, lrxpl_d, copy_to_device, plans_d,&
               blas_handle_p, stream_p )

          CALL cuda_stream_synchronize(stream_p)
          lda=lfrm*lr1m
          CALL fft_comm(yf_ptr,xf_ptr,lda,cntl%tr4a2a,comm)

          CALL fftcu_inv_full_2 ( f, xf_ptr, yf_ptr, parai%mepos, parai%nproc, cntl%tr4a2a,&
               t1_d, t2_d, sp8_d, msqf_d, copy_to_host, plans_d,&
               blas_handle_p, stream_p )

       ENDIF
    ELSE
       IF (sparse) THEN
          CALL fftcu_frw_sprs_1 ( f, xf_ptr, yf_ptr, parai%nproc, cntl%tr4a2a,&
               t1_d, t2_d, sp9_d, msqs_d, copy_to_device, plans_d,&
               blas_handle_p, stream_p )

          CALL cuda_stream_synchronize(stream_p)
          lda=lsrm*lr1m
          CALL fft_comm(xf_ptr,yf_ptr,lda,cntl%tr4a2a,comm)

          CALL fftcu_frw_sprs_2 ( f, xf_ptr, yf_ptr, parai%nproc, cntl%tr4a2a,&
               t1_d, t2_d, sp5_d, lrxpl_d, copy_to_host, plans_d,&
               blas_handle_p, stream_p )

       ELSE
          CALL fftcu_frw_full_1 ( f, xf_ptr, yf_ptr, parai%mepos, parai%nproc, cntl%tr4a2a,&
               t1_d, t2_d, sp8_d, msqf_d, copy_to_device, plans_d,&
               blas_handle_p, stream_p )

          CALL cuda_stream_synchronize(stream_p)
          lda=lfrm*lr1m
          CALL fft_comm(xf_ptr,yf_ptr,lda,cntl%tr4a2a,comm)

          CALL fftcu_frw_full_2 ( f, xf_ptr, yf_ptr, parai%nproc, cntl%tr4a2a,&
               t1_d, t2_d, sp5_d, lrxpl_d, copy_to_host, plans_d,&
               blas_handle_p, stream_p )

       ENDIF
    ENDIF

    IF( copy_to_host ) CALL cuda_stream_synchronize(stream_p)

    ! ==--------------------------------------------------------------==
  END SUBROUTINE fftnew_cuda

  ! ==================================================================
  SUBROUTINE mltfft(transa,transb,a,ldax,lday,b,ldbx,ldby,n,m,isign,scale,thread_view)
    ! ==--------------------------------------------------------------==
    CHARACTER(len=1)                         :: transa, transb
    COMPLEX(real_8)                          :: a(:)
    INTEGER                                  :: ldax, lday
    COMPLEX(real_8)                          :: b(:)
    INTEGER                                  :: ldbx, ldby, n, m, isign
    REAL(real_8)                             :: scale
    TYPE(thread_view_t), INTENT(IN), &
      OPTIONAL                               :: thread_view

    CHARACTER(*), PARAMETER                  :: procedureN = 'mltfft'

    INTEGER                                  :: device_idx, host_buff_ptr, &
                                                stream_idx
    TYPE(cp_cufft_plans_t), POINTER          :: plans_d
    TYPE(cublas_handle_t), POINTER           :: blas_handle_p
    TYPE(cuda_memory_t), POINTER             :: t1_d, t2_d
    TYPE(cuda_stream_t), POINTER             :: stream_p
    TYPE(cufft_plan_t), POINTER              :: plan_p

! ==--------------------------------------------------------------==
!vw for the moment we use CPU FFT 
!IF( cp_cuda_env%use_fft ) THEN

    IF( .FALSE. ) THEN

       PRINT *,'START DEBUG ----------------------------------------'

       IF( PRESENT( thread_view ) ) THEN
          device_idx = thread_view%device_idx
          stream_idx = thread_view%stream_idx - 1 !vw FIX that -1
          host_buff_ptr = thread_view%host_buff_ptr
       ELSE
          device_idx = 1
          stream_idx = 0
          !$ IF( omp_in_parallel() ) THEN
          !$    stream_idx = omp_get_thread_num()
          !    !$    WRITE(6,*) 'FFTNEW_CUDA - Stream n.: ',stream_idx ! acm: this is just for debug
          !$ ENDIF
          host_buff_ptr = stream_idx + 1 !vw FIX that +1
       ENDIF

       t1_d => cp_cufft%devices(device_idx)%streams(stream_idx)%t1
       t2_d => cp_cufft%devices(device_idx)%streams(stream_idx)%t2
       plans_d => cp_cufft%devices(device_idx)%streams(stream_idx)%plans
       stream_p => cp_cufft%devices(device_idx)%streams(stream_idx)%stream
       blas_handle_p => cp_cufft%devices(device_idx)%streams(stream_idx)%blas_handle

       PRINT *,'cuda_memcpy_async_host_to_device'

       CALL cuda_memcpy_async_host_to_device ( a(1:maxfftn), t1_d, stream_p )!vw check if maxfftn is correct

       PRINT *,'cp_cufft_get_plan'
       CALL cp_cufft_get_plan ( transa,transb,     ldax,lday,               n,m, plan_p, plans_d, stream_p )

       PRINT *,'mltfft_cuda'
       CALL mltfft_cuda(        transa,transb, t1_d,ldax,lday,t2_d,ldbx,ldby,n,m,isign,scale, plan_p, blas_handle_p, stream_p )

       PRINT *,'cuda_memcpy_async_device_to_host'
       WRITE (6,*) 'yf_d%init',t2_d%init
       WRITE (6,*) 'yf_d%idevice',t2_d%device
       WRITE (6,*) 'yf_d%n_bytes',t2_d%n_bytes
       WRITE (6,*) 'stream_p%init',stream_p%init
       WRITE (6,*) 'stream_p%device',stream_p%device

       CALL cuda_memcpy_async_device_to_host ( t2_d, b(1:maxfftn), stream_p )!vw check if maxfftn is correct

       PRINT *,'cuda_stream_synchronize'
       CALL cuda_stream_synchronize(stream_p)


       PRINT *,'DONE DEBUG ----------------------------------------'
    ELSE

#if defined(__HAS_FFT_DEFAULT)
       CALL mltfft_default(transa,transb,a,ldax,lday,b,ldbx,ldby,n,m,isign,scale)
#elif defined(__HAS_FFT_ESSL)
       CALL mltfft_essl(transa,transb,a,ldax,lday,b,ldbx,ldby,n,m,isign,scale)
#elif defined(__HAS_FFT_HP)
       CALL mltfft_hp(transa,transb,a,ldax,lday,b,ldbx,ldby,n,m,isign,scale)
#elif defined(__HAS_FFT_FFTW3)
       CALL mltfft_fftw(transa,transb,a,ldax,lday,b,ldbx,ldby,n,m,isign,scale,.FALSE.)
#else
       CALL stopgm(procedureN,"MLTFFT ROUTINE NOT AVAILABLE",&
            __LINE__,__FILE__)
#endif
    ENDIF
    ! ==--------------------------------------------------------------==
  END SUBROUTINE mltfft
  ! ==================================================================
  ! NEW FFT CODE
  ! ==================================================================
  SUBROUTINE invfftn(f,sparse,comm,thread_view, copy_data_to_device, copy_data_to_host )
    ! ==--------------------------------------------------------------==
    ! == COMPUTES THE INVERSE FOURIER TRANSFORM OF A COMPLEX          ==
    ! == FUNCTION F. THE FOURIER TRANSFORM IS                         ==
    ! == RETURNED IN F (THE INPUT F IS OVERWRITTEN).                  ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: f(:)
    LOGICAL                                  :: sparse
    INTEGER, INTENT(IN)                      :: comm
    TYPE(thread_view_t), INTENT(IN), &
      OPTIONAL                               :: thread_view
    LOGICAL, INTENT(IN), OPTIONAL            :: copy_data_to_device, &
                                                copy_data_to_host

    CHARACTER(*), PARAMETER                  :: procedureN = 'invfftn'

    INTEGER                                  :: isign, isub

    CALL tiset(procedureN,isub)
    isign=-1
    IF( cp_cuda_env%use_fft ) THEN
       CALL fftnew_cuda(isign,f,sparse, comm, thread_view=thread_view, &
            & copy_data_to_device=copy_data_to_device, copy_data_to_host=copy_data_to_host )
    ELSE
       CALL fftnew(isign,f,sparse, parai%allgrp )
    ENDIF
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE invfftn
  ! ==================================================================
  SUBROUTINE fwfftn(f,sparse,comm,thread_view, copy_data_to_device, copy_data_to_host )
    ! ==--------------------------------------------------------------==
    ! == COMPUTES THE FORWARD FOURIER TRANSFORM OF A COMPLEX          ==
    ! == FUNCTION F. THE FOURIER TRANSFORM IS                         ==
    ! == RETURNED IN F IN OUTPUT (THE INPUT F IS OVERWRITTEN).        ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: f(:)
    LOGICAL                                  :: sparse
    INTEGER, INTENT(IN)                      :: comm
    TYPE(thread_view_t), INTENT(IN), &
      OPTIONAL                               :: thread_view
    LOGICAL, INTENT(IN), OPTIONAL            :: copy_data_to_device, &
                                                copy_data_to_host

    CHARACTER(*), PARAMETER                  :: procedureN = 'fwfftn'

    INTEGER                                  :: isign, isub

    CALL tiset(procedureN,isub)
    isign=1
    IF( cp_cuda_env%use_fft ) THEN
       CALL fftnew_cuda(isign,f,sparse, comm, thread_view=thread_view, &
            & copy_data_to_device=copy_data_to_device, copy_data_to_host=copy_data_to_host )
    ELSE
       CALL fftnew(isign,f,sparse, parai%allgrp )
    ENDIF
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE fwfftn
  ! ==================================================================

END MODULE fftmain_utils
