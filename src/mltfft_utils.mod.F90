MODULE mltfft_utils

  USE cublas_types,                    ONLY: cublas_handle_t
  USE cublas_utils,                    ONLY: cublas_zdscal
  USE cuda_types,                      ONLY: cuda_memory_t,&
                                             cuda_stream_t
  USE cufft_interfaces,                ONLY: cufft_forward,&
                                             cufft_inverse
  USE cufft_types,                     ONLY: cufft_plan_t
  USE cufft_utils,                     ONLY: cufft_execz2z
  USE cuuser_utils,                    ONLY: CuUser_setblock2zero,&
                                             CuUser_setblock2zero_scale
  USE error_handling,                  ONLY: stopgm
  USE gfft_utils,                      ONLY: ctrig
  USE kinds,                           ONLY: int_4,&
                                             int_8,&
                                             real_8
  USE utils,                           ONLY: numcpus

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_NULL_CHAR,&
       C_INT, C_CHAR,&
       C_PTR,&
       C_LONG_DOUBLE_COMPLEX
  USE timer,                           ONLY: tihalt,&
       tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mltfft_default
  PUBLIC :: mltfft_essl
  PUBLIC :: mltfft_hp
  PUBLIC :: mltfft_fftw
  PUBLIC :: mltfft_cuda

CONTAINS

  ! ==================================================================
  SUBROUTINE mltfft_default(transa,transb,a,ldax,lday,b,ldbx,ldby,n,m,isign,scale)
    CHARACTER(len=1)                         :: transa, transb
    INTEGER                                  :: ldax, lday
    COMPLEX(real_8)                          :: a(ldax,lday)
    INTEGER                                  :: ldbx, ldby
    COMPLEX(real_8)                          :: b(ldbx,ldby)
    INTEGER                                  :: n, m, isign
    REAL(real_8)                             :: scale

    CHARACTER(*), PARAMETER                  :: procedureN = 'mltfft_default'

    INTEGER :: after(20), before(20), i, ic, icpu, ierr, il, inzee, ir, isig, &
      itr, j, length, lot, mh, ml, mm, ncpus, nfft, now(20)
    LOGICAL                                  :: tscal
    REAL(real_8)                             :: trig(2,1024)

#if defined(__NEC)
    INTEGER, PARAMETER :: ncache=1024*1024
#elif defined(__PRIMEHPC)
    INTEGER, PARAMETER :: ncache=1024*32
#elif defined(__VECTOR)
#ifdef __SR8000
    INTEGER, PARAMETER :: ncache=1024*8   ! cmb - for SR8K
#else
    INTEGER, PARAMETER :: ncache=1024*128
#endif
#elif defined(__SGI)
    INTEGER, PARAMETER :: ncache=1024*4
#elif defined(__IBM)
#ifdef __SR11KIBM
    INTEGER, PARAMETER :: ncache=1024*32  ! cmb - for SR11K
#else
    INTEGER, PARAMETER :: ncache=1024*10
#endif
#elif defined(__OSX) || defined(__OSX_IFC)
    INTEGER, PARAMETER :: ncache=1024*10
#elif defined(__i386)
    INTEGER, PARAMETER :: ncache=1024*10
#elif defined(__x86_64)
#ifdef __PACS
    INTEGER, PARAMETER :: ncache=1024*8   ! cmb - for PACS-CS
#else
    INTEGER, PARAMETER :: ncache=1024*6
#endif
#elif defined(__HP)
    INTEGER, PARAMETER :: ncache=1024*64
#elif defined(__ALTIX)
    INTEGER, PARAMETER :: ncache=1024*16  ! Altix3700/3900
#elif defined(_FFT_NCACHE)
    INTEGER, PARAMETER :: ncache=_fft_ncache  ! for convenient overrides in makefile
#else
    INTEGER, PARAMETER :: ncache=1024*2
#endif
#if defined(__SR8000) || defined(__SR11000)
    REAL(real_8) :: z(2,ncache/4+1,2,maxcpu/8)
#else
    REAL(real_8), ALLOCATABLE, SAVE :: z(:,:,:,:)
#endif
    INTEGER, SAVE :: ifirst=0

    ! ==-------------------------------------------------------------==
    CALL numcpus(ncpus)
    IF (ifirst.EQ.0) THEN
       ifirst=1
       length=4*(ncache/4+1)*ncpus
#if defined(__SR8000) || defined(__SR11000)
       ! mb     do not call MEMORY
#else
       ALLOCATE(z(2,ncache/4+1,2,ncpus),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
#endif
    ENDIF
    isig=-isign
    tscal=(ABS(scale-1._real_8).GT.1.e-12_real_8)
    CALL ctrig(n,trig,after,before,now,isig,ic)
    lot=ncache/(4*n)
    lot=lot-MOD(lot+1,2)
    lot=MAX(1,lot)
    IF (4*n*lot.GT.ncache) CALL stopgm('MLFFT_DEFAULT',' NCACHE ',&
         __LINE__,__FILE__)
    mm=MAX(m/ncpus,1)
#if defined(__SR11000) || defined(__BG)
    !$omp parallel do private(nfft,icpu,ml,mh,itr,inzee,i,j,il,ir) &
    !$omp  shared(trig,now,after,before,tscal,scale, &
    !$omp         ldax,lday,ldbx,ldby,isig,n,a,b,lot,mm,m,ncpus, &
    !$omp         transa,transb,ic)
#else
    !$omp parallel do private(nfft,icpu,ml,mh,itr,inzee,i) &
    !$omp  shared(trig,now,after,before,tscal,scale, &
    !$omp         ldax,lday,ldbx,ldby,isig,n,a,b,lot,mm,m,ncpus, &
    !$omp         transa,transb,ic)
#endif
    DO icpu=1,ncpus
       ml=(icpu-1)*mm+1
       mh=MIN(icpu*mm,m)
       IF (icpu.EQ.ncpus) mh=m
#ifdef __SR8000
       ! mb-warning: something might be screwed up on SR11000-J1
       !poption parallel
       !poption tlocal(NFFT,Z)
#endif
       DO itr=ml,mh,lot
          nfft=MIN(mh-itr+1,lot)
          IF (transa.EQ.'N'.OR.transa.EQ.'n') THEN
             CALL fftpre(nfft,nfft,ldax,lot,n,a(1,itr),z(1,1,1,icpu),&
                  trig,now(1),after(1),before(1),isig)
          ELSE
             CALL fftstp(ldax,nfft,n,lot,n,a(itr,1),z(1,1,1,icpu),&
                  trig,now(1),after(1),before(1),isig)
          ENDIF
          ! IF(TSCAL) CALL DSCAL(2*LOT*N,SCALE,Z(1,1,1,ICPU),1)
          IF (tscal) THEN
             IF (lot.EQ.nfft) THEN
#if defined(__SR8000) || defined(__SR11000) || defined(__BG)
                DO j=1,lot*n
                   z(1,j,1,icpu)=scale*z(1,j,1,icpu)
                   z(2,j,1,icpu)=scale*z(2,j,1,icpu)
                ENDDO
#else
                CALL dscal(2*lot*n,scale,z(1,1,1,icpu),1)
#endif
             ELSE
#if defined(__SR8000) || defined(__SR11000) || defined(__BG)
                DO i=1,n
                   DO j=1,nfft
                      il=lot*(i-1)+j
                      ir=il
                      z(1,il,1,icpu)=scale*z(1,ir,1,icpu)
                      z(2,il,1,icpu)=scale*z(2,ir,1,icpu)
                   ENDDO
                ENDDO
#else
                DO i=1,n
                   CALL dscal(2*nfft,scale,z(1,lot*(i-1)+1,1,icpu),1)
                ENDDO
#endif
             ENDIF
          ENDIF
          IF (ic.EQ.1) THEN
             IF (transb.EQ.'N'.OR.transb.EQ.'n') THEN
#if defined(__EREG)
                CALL cmplx_transpose(nfft,n,z(1,1,1,icpu),lot,&
                     b(itr,1),ldbx)
#elif defined(__SR8000) || defined(__SR11000) || defined(__BG)
                DO i=1,nfft
                   DO j=1,n
                      il=ldbx*(i-1)+j
                      ir= lot*(j-1)+i
                      b(il,itr)=CMPLX(z(1,ir,1,icpu),z(2,i,ir,icpu),kind=real_8)
                   ENDDO
                ENDDO
#else
                CALL zgetmo(z(1,1,1,icpu),lot,nfft,n,b(1,itr),ldbx)
#endif
             ELSE
#if defined(__SR8000) || defined(__SR11000) || defined(__BG)
                DO i=1,nfft
                   DO j=1,n
                      il=ldbx*(i-1)+j
                      ir= lot*(i-1)+j
                      b(itr,il)=CMPLX(z(1,ir,1,icpu),z(2,i,ir,icpu),kind=real_8)
                   ENDDO
                ENDDO
#else
                CALL matmov(nfft,n,z(1,1,1,icpu),lot,b(itr,1),ldbx)
#endif
             ENDIF
          ELSE
             inzee=1
             DO i=2,ic-1
                CALL fftstp(lot,nfft,n,lot,n,z(1,1,inzee,icpu),&
                     z(1,1,3-inzee,icpu),trig,now(i),after(i),&
                     before(i),isig)
                inzee=3-inzee
             ENDDO
             IF (transb.EQ.'N'.OR.transb.EQ.'n') THEN
                CALL fftrot(lot,nfft,n,nfft,ldbx,z(1,1,inzee,icpu),&
                     b(1,itr),trig,now(ic),after(ic),before(ic),isig)
             ELSE
                CALL fftstp(lot,nfft,n,ldbx,n,z(1,1,inzee,icpu),&
                     b(itr,1),trig,now(ic),after(ic),before(ic),isig)
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    IF (transb.EQ.'N'.OR.transb.EQ.'n') THEN
       !$omp parallel do private (I,J)
       DO j=m+1,ldby
          DO i=1,ldbx
             b(i,j)=CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
          ENDDO
       ENDDO
       !$omp parallel do private (I,J)
       DO j=1,m
          DO i=n+1,ldbx
             b(i,j)=CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
          ENDDO
       ENDDO
    ELSE
       !$omp parallel do private (I,J)
       DO j=n+1,ldby
          DO i=1,ldbx
             b(i,j)=CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
          ENDDO
       ENDDO
       !$omp parallel do private (I,J)
       DO j=1,n
          DO i=m+1,ldbx
             b(i,j)=CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
          ENDDO
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE mltfft_default
  ! ==================================================================
  SUBROUTINE mltfft_essl(transa,transb,a,ldax,lday,b,ldbx,ldby,n,m,&
       isign,scale)
    ! ==--------------------------------------------------------------==
    CHARACTER(len=1)                         :: transa, transb
    INTEGER                                  :: ldax, lday
    COMPLEX(real_8)                          :: a(ldax,lday)
    INTEGER                                  :: ldbx, ldby
    COMPLEX(real_8)                          :: b(ldbx,ldby)
    INTEGER                                  :: n, m, isign
    REAL(real_8)                             :: scale

    INTEGER, PARAMETER                       :: naux1 = 100000, naux2 = 100000

    INTEGER                                  :: i, j
    REAL(real_8)                             :: aux1(naux1), aux2(naux2)

#if defined(__HAS_FFT_ESSL)
    IF (transa.EQ.'N'.OR.transa.EQ.'n') THEN
       IF (transb.EQ.'N'.OR.transb.EQ.'n') THEN
          CALL dcft(1,a,1,ldax,b,1,ldbx,n,m,isign,scale,&
               aux1,naux1,aux2,naux2)
          CALL dcft(0,a,1,ldax,b,1,ldbx,n,m,isign,scale,&
               aux1,naux1,aux2,naux2)
       ELSE
          CALL dcft(1,a,1,ldax,b,ldbx,1,n,m,isign,scale,&
               aux1,naux1,aux2,naux2)
          CALL dcft(0,a,1,ldax,b,ldbx,1,n,m,isign,scale,&
               aux1,naux1,aux2,naux2)
       ENDIF
    ELSE
       IF (transb.EQ.'N'.OR.transb.EQ.'n') THEN
          CALL dcft(1,a,ldax,1,b,1,ldbx,n,m,isign,scale,&
               aux1,naux1,aux2,naux2)
          CALL dcft(0,a,ldax,1,b,1,ldbx,n,m,isign,scale,&
               aux1,naux1,aux2,naux2)
       ELSE
          CALL dcft(1,a,ldax,1,b,ldbx,1,n,m,isign,scale,&
               aux1,naux1,aux2,naux2)
          CALL dcft(0,a,ldax,1,b,ldbx,1,n,m,isign,scale,&
               aux1,naux1,aux2,naux2)
       ENDIF
    ENDIF
#endif
    ! ==--------------------------------------------------------------==
    IF (transb.EQ.'N'.OR.transb.EQ.'n') THEN
       DO j=m+1,ldby
          DO i=1,ldbx
             b(i,j)=CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
          ENDDO
       ENDDO
       DO j=1,m
          DO i=n+1,ldbx
             b(i,j)=CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
          ENDDO
       ENDDO
    ELSE
       DO j=n+1,ldby
          DO i=1,ldbx
             b(i,j)=CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
          ENDDO
       ENDDO
       DO j=1,n
          DO i=m+1,ldbx
             b(i,j)=CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
          ENDDO
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE mltfft_essl
  ! ==================================================================
  SUBROUTINE mltfft_hp(transa,transb,a,ldax,lday,b,ldbx,ldby,n,m,&
       isign,scale)
    ! ==--------------------------------------------------------------==
    CHARACTER(len=1)                         :: transa, transb
    INTEGER                                  :: ldax, lday
    COMPLEX(real_8)                          :: a(ldax,lday)
    INTEGER                                  :: ldbx, ldby
    COMPLEX(real_8)                          :: b(ldbx,ldby)
    INTEGER                                  :: n, m, isign
    REAL(real_8)                             :: scale

    INTEGER(int_4)                           :: i, ier, iopt, j
    REAL(real_8)                             :: ascale

! ==--------------------------------------------------------------==

#if defined(__HAS_FFT_HP)
    IF (transa.EQ.'N'.OR.transa.EQ.'n') THEN
       IF (transb.EQ.'N'.OR.transb.EQ.'n') THEN
          ! CALL DCFT(1,A,1,LDAX,B,1,LDBX,N,M,ISIGN,SCALE,
          ! *            AUX1,NAUX1,AUX2,NAUX2)
          ! CALL DCFT(0,A,1,LDAX,B,1,LDBX,N,M,ISIGN,SCALE,
          ! *            AUX1,NAUX1,AUX2,NAUX2)
          ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          CALL zgecpy ('N', n, m, a, ldax, b, ldbx)
          IF (isign .GE. 0) THEN
             iopt = -1
          ELSE
             iopt = 1
          ENDIF
          CALL zffts (b, n, 1, m, ldbx, iopt, ier)
          IF (ier .NE. 0) THEN
             CALL stopgm('MLFFT_HP','ERROR IN CALL TO ZFFTS',&
                  __LINE__,__FILE__)
          ENDIF
          IF (isign .GE. 0) THEN
             ascale = scale*REAL(n)
          ELSE
             ascale = scale
          ENDIF
          IF (ascale .NE. 1.0_real_8) THEN
             DO i = 1,m
                CALL zdscal(n, ascale, b(1,i), 1)
             ENDDO
          ENDIF
          ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       ELSE
          ! CALL DCFT(1,A,1,LDAX,B,LDBX,1,N,M,ISIGN,SCALE,
          ! *            AUX1,NAUX1,AUX2,NAUX2)
          ! CALL DCFT(0,A,1,LDAX,B,LDBX,1,N,M,ISIGN,SCALE,
          ! *            AUX1,NAUX1,AUX2,NAUX2)
          ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          CALL zgecpy ('T', n, m, a, ldax, b, ldbx)
          IF (isign .GE. 0) THEN
             iopt = -1
          ELSE
             iopt = 1
          ENDIF
          CALL zffts (b, n, ldbx, m, 1, iopt, ier)
          IF (ier .NE. 0) THEN
             CALL stopgm('MLFFT_HP','ERROR IN CALL TO ZFFTS',&
                  __LINE__,__FILE__)
          ENDIF
          IF (isign .GE. 0) THEN
             ascale = scale*REAL(n)
          ELSE
             ascale = scale
          ENDIF
          IF (ascale .NE. 1.0_real_8) THEN
             DO i = 1,n
                CALL zdscal(m, ascale, b(1,i), 1)
             ENDDO
          ENDIF
          ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       ENDIF
    ELSE
       IF (transb.EQ.'N'.OR.transb.EQ.'n') THEN
          ! CALL DCFT(1,A,LDAX,1,B,1,LDBX,N,M,ISIGN,SCALE,
          ! *            AUX1,NAUX1,AUX2,NAUX2)
          ! CALL DCFT(0,A,LDAX,1,B,1,LDBX,N,M,ISIGN,SCALE,
          ! *            AUX1,NAUX1,AUX2,NAUX2)
          ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          CALL zgecpy ('T', m, n, a, ldax, b, ldbx)
          IF (isign .GE. 0) THEN
             iopt = -1
          ELSE
             iopt = 1
          ENDIF
          CALL zffts (b, n, 1, m, ldbx, iopt, ier)
          IF (ier .NE. 0) THEN
             CALL stopgm('MLFFT_HP','ERROR IN CALL TO ZFFTS',&
                  __LINE__,__FILE__)
          ENDIF
          IF (isign .GE. 0) THEN
             ascale = scale*REAL(n)
          ELSE
             ascale = scale
          ENDIF
          IF (ascale .NE. 1.0_real_8) THEN
             DO i = 1,m
                CALL zdscal(n, ascale, b(1,i), 1)
             ENDDO
          ENDIF
          ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       ELSE
          ! CALL DCFT(1,A,LDAX,1,B,LDBX,1,N,M,ISIGN,SCALE,
          ! *            AUX1,NAUX1,AUX2,NAUX2)
          ! CALL DCFT(0,A,LDAX,1,B,LDBX,1,N,M,ISIGN,SCALE,
          ! *            AUX1,NAUX1,AUX2,NAUX2)
          ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          CALL zgecpy ('N', m, n, a, ldax, b, ldbx)
          IF (isign .GE. 0) THEN
             iopt = -1
          ELSE
             iopt = 1
          ENDIF
          CALL zffts (b, n, ldbx, m, 1, iopt, ier)
          IF (ier .NE. 0) THEN
             CALL stopgm('MLFFT_HP','ERROR IN CALL TO ZFFTS',&
                  __LINE__,__FILE__)
          ENDIF
          IF (isign .GE. 0) THEN
             ascale = scale*REAL(n)
          ELSE
             ascale = scale
          ENDIF
          IF (ascale .NE. 1.0_real_8) THEN
             DO i = 1,n
                CALL zdscal(m, ascale, b(1,i), 1)
             ENDDO
          ENDIF
          ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       ENDIF
    ENDIF
#endif
    ! ==--------------------------------------------------------------==
    IF (transb.EQ.'N'.OR.transb.EQ.'n') THEN
       DO j=m+1,ldby
          DO i=1,ldbx
             b(i,j)=CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
          ENDDO
       ENDDO
       DO j=1,m
          DO i=n+1,ldbx
             b(i,j)=CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
          ENDDO
       ENDDO
    ELSE
       DO j=n+1,ldby
          DO i=1,ldbx
             b(i,j)=CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
          ENDDO
       ENDDO
       DO j=1,n
          DO i=m+1,ldbx
             b(i,j)=CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
          ENDDO
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE mltfft_hp
  ! ==================================================================
  SUBROUTINE mltfft_fftw(transa,transb,a,ldax,lday,b,ldbx,ldby,n,m,&
       isign,scale,use_wisdom)
    ! ==--------------------------------------------------------------==
    CHARACTER(len=1)                         :: transa, transb
    INTEGER                                  :: ldax
    COMPLEX(real_8)                          :: a(ldax,*)
    INTEGER                                  :: lday, ldbx
    COMPLEX(real_8)                          :: b(ldbx,*)
    INTEGER                                  :: ldby, n, m, isign
    REAL(real_8)                             :: scale
    LOGICAL                                  :: use_wisdom

    INTEGER, PARAMETER                       :: fftw_b = 1, fftw_es = 64, &
                                                fftw_f = -1 , fftw_me = 0, &
                                                fftw_pa = 32

    INTEGER                                  :: fftw_dir, fftw_flags, i, j
    LOGICAL                                  :: tscal

    TYPE(C_PTR) :: plan
#if defined(__HAS_FFT_FFTW3)
    tscal=(ABS(scale-1._real_8).GT.1.e-12_real_8)
    IF (isign.EQ.1) THEN
       fftw_dir=fftw_f
    ELSE
       fftw_dir=fftw_b
    ENDIF
    IF (use_wisdom) THEN
       fftw_flags=fftw_me
    ELSE
       fftw_flags=fftw_es
    ENDIF

    IF (transa.EQ.'N'.OR.transa.EQ.'n') THEN
       IF (transb.EQ.'N'.OR.transb.EQ.'n') THEN
          CALL dfftw_plan_many_dft(plan,1,n,m,a,n,1,ldax,b,n,1,ldbx,&
               fftw_dir,fftw_flags)
          CALL dfftw_execute_dft(plan,a,b)
          IF (tscal) THEN
             !$omp parallel do private(I,J)
             DO i = 1,m
                DO j = 1,n
                   b(j,i)=scale*b(j,i)
                ENDDO
             ENDDO
          ENDIF
       ELSE
          CALL dfftw_plan_many_dft(plan,1,n,m,a,n,1,ldax,b,n,ldbx,1,&
               fftw_dir,fftw_flags)
          CALL dfftw_execute_dft(plan,a,b)
          IF (tscal) THEN
             !$omp parallel do private(I,J)
             DO i = 1,n
                DO j = 1,m
                   b(j,i)=scale*b(j,i)
                ENDDO
             ENDDO
          ENDIF
       ENDIF
    ELSE
       IF (transb.EQ.'N'.OR.transb.EQ.'n') THEN
          CALL dfftw_plan_many_dft(plan,1,n,m,a,n,ldax,1,b,n,1,ldbx,&
               fftw_dir,fftw_flags)
          CALL dfftw_execute_dft(plan,a,b)
          IF (tscal) THEN
             !$omp parallel do private(I,J)
             DO i = 1,m
                DO j = 1,n
                   b(j,i)=scale*b(j,i)
                ENDDO
             ENDDO
          ENDIF
       ELSE
          CALL dfftw_plan_many_dft(plan,1,n,m,a,n,ldax,1,b,n,ldbx,1,&
               fftw_dir,fftw_flags)
          CALL dfftw_execute_dft(plan,a,b)
          IF (tscal) THEN
             !$omp parallel do private(I,J)
             DO i = 1,n
                DO j = 1,m
                   b(j,i)=scale*b(j,i)
                ENDDO
             ENDDO
          ENDIF
       ENDIF
    ENDIF
    CALL dfftw_destroy_plan(plan)
#endif
    IF (transb.EQ.'N'.OR.transb.EQ.'n') THEN
       !$omp parallel do private (I,J)
       DO j=m+1,ldby
          DO i=1,ldbx
             b(i,j)=CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
          ENDDO
       ENDDO
       !$omp parallel do private (I,J)
       DO j=1,m
          DO i=n+1,ldbx
             b(i,j)=CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
          ENDDO
       ENDDO
    ELSE
       !$omp parallel do private (I,J)
       DO j=n+1,ldby
          DO i=1,ldbx
             b(i,j)=CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
          ENDDO
       ENDDO
       !$omp parallel do private (I,J)
       DO j=1,n
          DO i=m+1,ldbx
             b(i,j)=CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
          ENDDO
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE mltfft_fftw
  ! ==================================================================

  SUBROUTINE mltfft_cuda(transa,transb,a,ldax,lday,b,ldbx,ldby,n,m,isign,scale, plan, blas_handle, stream)
    !     ==--------------------------------------------------------------==
    CHARACTER(len=1)                         :: transa, transb
    TYPE(cuda_memory_t)                      :: a
    INTEGER                                  :: ldax, lday
    TYPE(cuda_memory_t)                      :: b
    INTEGER                                  :: ldbx, ldby, n, m, isign
    REAL(real_8)                             :: scale
    TYPE(cufft_plan_t), INTENT(IN)           :: plan
    TYPE(cublas_handle_t), INTENT(IN)        :: blas_handle
    TYPE(cuda_stream_t), INTENT(IN)          :: stream

    CHARACTER(*), PARAMETER                  :: procedureN = 'mltfft_cuda'

    INTEGER                                  :: cufft_dir
    LOGICAL                                  :: tscal

!     ==--------------------------------------------------------------==

    IF( plan%device /= blas_handle%device .OR. &
         plan%device /= stream%device ) CALL stopgm(procedureN,'devices are not consistent',&
         __LINE__,__FILE__)


    IF (isign.EQ.1) THEN
       cufft_dir=cufft_forward
    ELSE
       cufft_dir=cufft_inverse
    ENDIF

    CALL cufft_execz2z(plan, a, b, cufft_dir)

    tscal=(ABS(scale-1._real_8).GT.1.e-12_real_8)
    IF( .TRUE. ) THEN
       IF(tscal) CALL cublas_zdscal( blas_handle, ldax*lday, scale, b, 1); !acm: should it be ldbx*ldby ?
       CALL CuUser_setblock2zero ( b, transb, n, m, ldbx, ldby, stream )
    ELSE
       IF(tscal) THEN
          CALL CuUser_setblock2zero_scale ( b, scale, transb, n, m, ldbx, ldby, stream )
       ELSE
          CALL CuUser_setblock2zero ( b, transb, n, m, ldbx, ldby, stream )
       ENDIF
    ENDIF

  END SUBROUTINE mltfft_cuda

END MODULE mltfft_utils
