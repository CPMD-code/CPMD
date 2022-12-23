#include "cpmd_global.h"

MODULE fft

  USE, INTRINSIC :: ISO_C_BINDING,     ONLY: C_PTR, C_NULL_PTR
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE


  ! ==================================================================
  INTEGER, ALLOCATABLE :: mg(:,:)
  INTEGER, ALLOCATABLE, TARGET :: ms(:,:)
  INTEGER, ALLOCATABLE :: mz(:)

  INTEGER :: naux1,naux2

#if defined(_HAS_CUDA)
  COMPLEX(real_8), POINTER __CONTIGUOUS :: xf(:,:)
  COMPLEX(real_8), POINTER __CONTIGUOUS :: yf(:,:)
#else
  COMPLEX(real_8), ALLOCATABLE, TARGET :: xf(:,:)
  COMPLEX(real_8), ALLOCATABLE, TARGET :: yf(:,:)
#endif

  ! ==================================================================
  INTEGER :: nr1m,nr2m,nr3m,kr1m,kr2m,kr3m,nhrm,ngrm,&
       kr2min,kr2max,kr3min,kr3max,maxrpt
  INTEGER, DIMENSION(:), ALLOCATABLE :: mxy !(maxcpu)
  ! ==================================================================
  ! ==   GATHER/SCATTER ARRAYS                                      ==
  ! ==================================================================
  INTEGER, ALLOCATABLE :: msp(:,:,:)

  ! ==================================================================
  ! NEW GENERAL PARALLEL FFT CODE
  ! ==================================================================
  INTEGER :: msrays,mfrays,llr1
  INTEGER :: qr1s,qr2s,qr3s,qr1,qr2,qr3
  INTEGER :: lr1s,lr2s,lr3s,lr1,lr2,lr3
  INTEGER :: qr2max,qr2min,qr3max,qr3min
  INTEGER :: lsrm,lfrm,lr1m,lmsq
  INTEGER :: jgw,jgws,jhg,jhgs
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: lrxpl !(0:maxcpu,2)
  INTEGER, DIMENSION(:), ALLOCATABLE :: sp5,sp8,sp9!(0:maxcpu)
  INTEGER, ALLOCATABLE :: msqs(:)
  INTEGER, ALLOCATABLE :: msqf(:)

  INTEGER, POINTER :: nzff(:)
  INTEGER, POINTER :: nzfs(:)
  INTEGER, POINTER :: inzf(:)
  INTEGER, POINTER :: inzs(:)
  INTEGER, POINTER :: inzh(:,:)

  ! ==================================================================
  ! POOL
  ! ==================================================================
  INTEGER, PARAMETER :: fftpoolsize=3
  INTEGER :: fftpool
  INTEGER :: fpoolv(28,fftpoolsize)
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: lrxpool!(0:maxcpu,2,fftpoolsize)
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: spm!(9,0:maxcpu,fftpoolsize)
  INTEGER :: lmsqmax,lnzf,lnzs
  INTEGER, ALLOCATABLE :: msqspool(:,:,:)
  INTEGER, ALLOCATABLE :: msqfpool(:,:,:)

  INTEGER, ALLOCATABLE, TARGET :: nzffp(:,:)
  INTEGER, ALLOCATABLE, TARGET :: nzfsp(:,:)
  INTEGER, ALLOCATABLE, TARGET :: inzfp(:,:)
  INTEGER, ALLOCATABLE, TARGET :: inzsp(:,:)
  INTEGER, ALLOCATABLE, TARGET :: inzhp(:,:,:)

  ! ==================================================================

END MODULE fft

