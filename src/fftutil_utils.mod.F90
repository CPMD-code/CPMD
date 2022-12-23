#if defined(__FFT_HASNT_THREADED_COPIES)
#define HASNT_THREADED_COPIES .TRUE.
#else
#define HASNT_THREADED_COPIES .FALSE.
#endif

#if defined(__FFT_HAS_LOW_LEVEL_TIMERS)
#define HAS_LOW_LEVEL_TIMERS .TRUE.
#else
#define HAS_LOW_LEVEL_TIMERS .FALSE.
#endif

#if defined(__FFT_HAS_SPECIAL_COPY)
#define HAS_SPECIAL_COPY __FFT_HAS_SPECIAL_COPY
#else
#define HAS_SPECIAL_COPY 0
#endif

#if defined(__FFT_HAS_OMP_COLLAPSE)
#define __COLLAPSE2 collapse(2)
#else
#define __COLLAPSE2
#endif


MODULE fftutil_utils
  USE fft_maxfft,                      ONLY: maxfft
  USE kinds,                           ONLY: real_4,&
                                             real_8
  USE mp_interface,                    ONLY: mp_all2all
  USE parac,                           ONLY: parai
  USE reshaper,                        ONLY: type_cast
  USE system,                          ONLY: fpar,&
                                             parap,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: zgthr_no_omp,&
                                             zsctr_no_omp
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: phase
  PUBLIC :: putz
  PUBLIC :: getz
  PUBLIC :: unpack_y2x
  PUBLIC :: pack_y2x
  PUBLIC :: fft_comm
  PUBLIC :: pack_x2y
  PUBLIC :: unpack_x2y
  PUBLIC :: phasen

CONTAINS


  ! ==================================================================
  SUBROUTINE phase(f)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8) :: f(fpar%kr1,fpar%kr2s,fpar%kr3s)

    INTEGER                                  :: ii, ijk1, ijk2, isub, j, k

! ==--------------------------------------------------------------==

    IF (HAS_LOW_LEVEL_TIMERS) CALL tiset('     PHASE',isub)
    ijk2=parap%nrxpl(parai%mepos,2)-parap%nrxpl(parai%mepos,1)+1 ! jh-mb
    !$omp parallel do private (J,K,II,IJK1) shared(IJK2) __COLLAPSE2
#ifdef __SR8000
    !poption tlocal(K,J,II,IJK1)
#endif
    DO k=1,spar%nr3s
       DO j=1,spar%nr2s
          ii=k+j+parap%nrxpl(parai%mepos,1)
          ijk1=MOD(ii,2)+1
          DO ii=ijk1,ijk2,2
             f(ii,j,k)=-f(ii,j,k)
          ENDDO
       ENDDO
    ENDDO
    IF (HAS_LOW_LEVEL_TIMERS) CALL tihalt('     PHASE',isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE phase
  ! ==================================================================
  SUBROUTINE putz(a,b,krmin,krmax,kr,m)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: krmin, krmax
    COMPLEX(real_8)                          :: a(krmax-krmin+1,*)
    INTEGER                                  :: kr, m
    COMPLEX(real_8)                          :: b(kr,m)

    CHARACTER(*), PARAMETER                  :: procedureN = 'putz'

    INTEGER                                  :: isub, n

    IF (HAS_LOW_LEVEL_TIMERS) CALL tiset(procedureN,isub)
    n=krmax-krmin+1
    CALL zeroing(b)!,m*kr)
    CALL matmov(n,m,a,n,b(krmin,1),kr)
    IF (HAS_LOW_LEVEL_TIMERS) CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE putz
  ! ==================================================================
  SUBROUTINE getz(a,b,krmin,krmax,kr,m)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: krmin, krmax
    COMPLEX(real_8)                          :: b(krmax-krmin+1,*)
    INTEGER                                  :: kr
    COMPLEX(real_8)                          :: a(kr,*)
    INTEGER                                  :: m

    CHARACTER(*), PARAMETER                  :: procedureN = 'getz'

    INTEGER                                  :: isub, n

! ==--------------------------------------------------------------==

    IF (HAS_LOW_LEVEL_TIMERS) CALL tiset(procedureN,isub)
    n=krmax-krmin+1
    CALL matmov(n,m,a(krmin,1),kr,b,n)
    IF (HAS_LOW_LEVEL_TIMERS) CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE getz
  ! ==================================================================
  ! ==================================================================
  ! CODE FOR NEW FFT ROUTINES
  ! ==================================================================


  SUBROUTINE unpack_y2x(xf,yf,m,nrays,lda,jrxpl,sp5,maxfft,mproc,&
       tr4a2a)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: xf(*), yf(*)
    INTEGER                                  :: m, nrays, lda, maxfft, mproc, &
                                                sp5(0:mproc-1), &
                                                jrxpl(0:mproc-1)
    LOGICAL                                  :: tr4a2a

    COMPLEX(real_4), POINTER                 :: yf4(:)
    INTEGER                                  :: ip, ipp, isub1, k, nrs, nrx

    !$    INTEGER   max_threads
    !$    INTEGER, EXTERNAL :: omp_get_max_threads
    CHARACTER(*),PARAMETER :: procedureN='UNPACK_Y2X'
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    ! ..Pack the data for sending
    IF (HAS_LOW_LEVEL_TIMERS) CALL tiset(procedureN,isub1)
    !$    IF(HASNT_THREADED_COPIES) THEN
    !$       max_threads = omp_get_max_threads()
    !$       call omp_set_num_threads(1)
    !$    ENDIF
    IF (tr4a2a) THEN
       CALL type_cast(yf, maxfft, yf4)
       SELECT CASE(HAS_SPECIAL_COPY)
       CASE default
          DO ip=0,mproc-1
             nrx = sp5(ip)*nrays
             nrs = (jrxpl(ip)-1)*nrays + 1
             ipp = ip*lda + 1
             CALL dcopy_s(2*nrx,yf4(ipp),1,xf(nrs),1)
          ENDDO
       CASE(1)
          ! that seems to work better on the P
          !$omp parallel do shared(SP5,NRAYS,JRXPL,LDA) &
          !$omp             private(NRX,NRS,IPP,K)
          DO ip=0,mproc-1
             nrx = sp5(ip)*nrays
             nrs = (jrxpl(ip)-1)*nrays
             ipp = ip*lda
             DO k=1,nrx
                xf(nrs+k)=yf4(ipp+k)
             ENDDO
          ENDDO
       END SELECT
    ELSE
       SELECT CASE(HAS_SPECIAL_COPY)
       CASE default
          DO ip=0,mproc-1
             nrx = sp5(ip)*nrays
             nrs = (jrxpl(ip)-1)*nrays + 1
             ipp = ip*lda + 1
             CALL dcopy(2*nrx,yf(ipp),1,xf(nrs),1)
          ENDDO
       CASE(1)
          ! that seems to work better on the P
          !$omp parallel do shared(SP5,NRAYS,JRXPL,LDA) &
          !$omp             private(NRX,NRS,IPP,K)
          DO ip=0,mproc-1
             nrx = sp5(ip)*nrays
             nrs = (jrxpl(ip)-1)*nrays
             ipp = ip*lda
             DO k=1,nrx
                xf(nrs+k)=yf(ipp+k)
             ENDDO
          ENDDO
       END SELECT
    ENDIF
    !$    IF(HASNT_THREADED_COPIES) call omp_set_num_threads(max_threads)
    IF (HAS_LOW_LEVEL_TIMERS) CALL tihalt(procedureN,isub1)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE unpack_y2x
  ! ==================================================================
  SUBROUTINE pack_y2x(xf,yf,m,lr1,lda,msp,lmsp,sp8,maxfft,mproc,&
       tr4a2a)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: xf(*), yf(*)
    INTEGER                                  :: m, lr1, lda, lmsp, &
                                                msp(lmsp,*), maxfft, mproc, &
                                                sp8(0:mproc-1)
    LOGICAL                                  :: tr4a2a

    COMPLEX(real_4), POINTER                 :: xf4(:)
    INTEGER                                  :: i, ii, ip, isub1, jj, k, &
                                                mxrp, nrx

    !$    INTEGER   max_threads
    !$    INTEGER, EXTERNAL :: omp_get_max_threads
    CHARACTER(*),PARAMETER :: procedureN='PACK_Y2X'
    ! ==--------------------------------------------------------------==
    IF (HAS_LOW_LEVEL_TIMERS) CALL tiset(procedureN,isub1)
    ! ..Pack the data for sending
    ! IF(HAS_LOW_LEVEL_TIMERS) CALL TISET(procedureN,ISUB1)
    !$    IF(HASNT_THREADED_COPIES) THEN
    !$       max_threads = omp_get_max_threads()
    !$       call omp_set_num_threads(1)
    !$    ENDIF
    IF (tr4a2a) THEN
       CALL type_cast(xf, maxfft, xf4)
       SELECT CASE(HAS_SPECIAL_COPY)
       CASE default
          DO ip=0,mproc-1
             mxrp = sp8(ip)
             DO i=1,lr1
                ii = ip*lda + (i-1)*mxrp + 1
                jj = (i-1)*m + 1
                CALL cgthr_z(mxrp,yf(jj),xf4(ii),msp(1,ip+1))
             ENDDO
          ENDDO
       CASE(1)
          ! that seems to work better on the P
          !$omp parallel do shared(SP8,LR1,M,LDA) &
          !$omp             private(MXRP,I,II,JJ,K)
          DO ip=0,mproc-1
             mxrp = sp8(ip)
             DO i=1,lr1
                ii = ip*lda + (i-1)*mxrp
                jj = (i-1)*m
                DO k=1,mxrp
                   xf4(ii+k) = yf(jj+msp(k,ip+1))
                ENDDO
             ENDDO
          ENDDO
       END SELECT
    ELSE
       SELECT CASE(HAS_SPECIAL_COPY)
       CASE default

          !$omp parallel do shared(SP8,NRX,M,LDA) &
          !$omp             private(MXRP,I,II,JJ)
          DO ip=0,mproc-1
             mxrp = sp8(ip)
             DO i=1,lr1
                ii = ip*lda + (i-1)*mxrp + 1
                jj = (i-1)*m + 1
                CALL zgthr_no_omp(mxrp,yf(jj),xf(ii),msp(1,ip+1))
             ENDDO
          ENDDO
       CASE(1)
          ! that seems to work better on the P
          !$omp parallel do shared(SP8,LR1,M,LDA) &
          !$omp             private(MXRP,I,II,JJ,K)
          DO ip=0,mproc-1
             mxrp = sp8(ip)
             DO i=1,lr1
                ii = ip*lda + (i-1)*mxrp
                jj = (i-1)*m
                DO k=1,mxrp
                   xf(ii+k) = yf(jj+msp(k,ip+1))
                ENDDO
             ENDDO
          ENDDO
       END SELECT
    ENDIF
    !$    IF(HASNT_THREADED_COPIES) call omp_set_num_threads(max_threads)
    IF (HAS_LOW_LEVEL_TIMERS) CALL tihalt(procedureN,isub1)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE pack_y2x
  ! ==================================================================
  SUBROUTINE fft_comm(xf,yf,lda,tr4a2a, comm )
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8), TARGET                  :: xf(*), yf(*)
    INTEGER, INTENT(IN)                      :: lda
    LOGICAL, INTENT(IN)                      :: tr4a2a
    INTEGER, INTENT(IN)                      :: comm

    CHARACTER(*), PARAMETER                  :: procedureN = 'fft_comm'

    COMPLEX(real_4), POINTER                 :: xf4(:), yf4(:)
    INTEGER                                  :: isub1

! Variables
! ==--------------------------------------------------------------==
! ..All to all communication

    IF (HAS_LOW_LEVEL_TIMERS) CALL tiset(procedureN,isub1)
    IF (tr4a2a) THEN
       ! is that needed on P?
       CALL type_cast(xf, maxfft, xf4)
       CALL type_cast(yf, maxfft, yf4)       
       CALL mp_all2all( xf4, yf4, lda, comm )
    ELSE
       CALL mp_all2all( xf, yf, lda, comm )
    ENDIF
    IF (HAS_LOW_LEVEL_TIMERS) CALL tihalt(procedureN,isub1)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE fft_comm
  ! ==================================================================
  SUBROUTINE pack_x2y(xf,yf,nrays,lda,jrxpl,sp5,&
       maxfft,mproc,tr4a2a)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: xf(*), yf(*)
    INTEGER                                  :: nrays, lda, maxfft, mproc, &
                                                sp5(0:mproc-1), &
                                                jrxpl(0:mproc-1)
    LOGICAL                                  :: tr4a2a

    COMPLEX(real_4), POINTER                 :: yf4(:)
    INTEGER                                  :: ip, ipp, isub1, k, nrs, nrx

    !$    INTEGER   max_threads
    !$    INTEGER, EXTERNAL :: omp_get_max_threads
    CHARACTER(*),PARAMETER :: procedureN='PACK_X2Y'
    ! ==--------------------------------------------------------------==
    IF (HAS_LOW_LEVEL_TIMERS) CALL tiset(procedureN,isub1)
    ! ..Prepare data for sending
    !$    IF(HASNT_THREADED_COPIES) THEN
    !$       max_threads = omp_get_max_threads()
    !$       call omp_set_num_threads(1)
    !$    ENDIF
    IF (tr4a2a) THEN
       CALL type_cast(yf, maxfft, yf4)
       SELECT CASE(HAS_SPECIAL_COPY)
       CASE default
          DO ip=0,mproc-1
             nrx = sp5(ip)*nrays
             nrs = (jrxpl(ip)-1)*nrays + 1
             ipp = ip*lda + 1
             CALL scopy_d(2*nrx,xf(nrs),1,yf4(ipp),1)
          ENDDO
       CASE(1)
          ! that seems to work better on the P
          !$omp parallel do shared(SP5,NRAYS,JRXPL,LDA) &
          !$omp             private(NRX,NRS,IPP,K)
          DO ip=0,mproc-1
             nrx = sp5(ip)*nrays
             nrs = (jrxpl(ip)-1)*nrays
             ipp = ip*lda
             DO k=1,nrx
                yf4(ipp+k)=xf(nrs+k)
             ENDDO
          ENDDO
       END SELECT
    ELSE
       SELECT CASE(HAS_SPECIAL_COPY)
       CASE default
          DO ip=0,mproc-1
             nrx = sp5(ip)*nrays
             nrs = (jrxpl(ip)-1)*nrays + 1
             ipp = ip*lda + 1
             CALL dcopy(2*nrx,xf(nrs),1,yf(ipp),1)
          ENDDO
       CASE(1)
          ! that seems to work better on the P
          !$omp parallel do shared(SP5,NRAYS,JRXPL,LDA) &
          !$omp             private(NRX,NRS,IPP,K)
          DO ip=0,mproc-1
             nrx = sp5(ip)*nrays
             nrs = (jrxpl(ip)-1)*nrays
             ipp = ip*lda
             DO k=1,nrx
                yf(ipp+k)=xf(nrs+k)
             ENDDO
          ENDDO
       END SELECT
    ENDIF
    !$    IF(HASNT_THREADED_COPIES) call omp_set_num_threads(max_threads)
    IF (HAS_LOW_LEVEL_TIMERS) CALL tihalt(procedureN,isub1)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE pack_x2y
  ! ==================================================================
  SUBROUTINE unpack_x2y(xf,yf,m,lr1,lda,msp,lmsp,sp8,maxfft,mproc,tr4a2a)
    ! ==--------------------------------------------------------------==
    ! include 'parac.inc'
    COMPLEX(real_8)                          :: xf(*)
    INTEGER                                  :: m, lr1, lda, lmsp, &
                                                msp(lmsp,*), maxfft
    COMPLEX(real_8)                          :: yf(maxfft)
    INTEGER                                  :: mproc, sp8(0:mproc-1)
    LOGICAL                                  :: tr4a2a

    COMPLEX(real_4), POINTER                 :: xf4(:)
    INTEGER                                  :: i, ii, ip, isub1, jj, k, mxrp

    !$    INTEGER   max_threads
    !$    INTEGER, EXTERNAL :: omp_get_max_threads
    CHARACTER(*),PARAMETER :: procedureN='UNPACK_X2Y'
    ! ==--------------------------------------------------------------==
    IF (HAS_LOW_LEVEL_TIMERS) CALL tiset(procedureN,isub1)
    ! ..Unpacking the data
    CALL zeroing(yf)!,maxfft)
    !$    IF(HASNT_THREADED_COPIES) THEN
    !$       max_threads = omp_get_max_threads()
    !$       call omp_set_num_threads(1)
    !$    ENDIF
    IF (tr4a2a) THEN
       CALL type_cast(xf, maxfft, xf4)
       SELECT CASE(HAS_SPECIAL_COPY)
       CASE default
          DO ip=0,mproc-1
             mxrp = sp8(ip)
             DO i=1,lr1
                ii = ip*lda + (i-1)*mxrp + 1
                jj = (i-1)*m + 1
                CALL zsctr_c(mxrp,xf4(ii),msp(1,ip+1),yf(jj))
             ENDDO
          ENDDO
       CASE(1)
          ! that seems to work better on the P
          !$omp parallel do shared(SP8,LR1,M,LDA) &
          !$omp             private(MXRP,I,II,JJ,K)
          DO ip=0,mproc-1
             mxrp = sp8(ip)
             DO i=1,lr1
                ii = ip*lda + (i-1)*mxrp
                jj = (i-1)*m
                DO k=1,mxrp
                   yf(jj+msp(k,ip+1)) = xf4(ii+k)
                ENDDO
             ENDDO
          ENDDO
       END SELECT
    ELSE
       SELECT CASE(HAS_SPECIAL_COPY)
       CASE default
          !$omp parallel do shared(SP8,LR1,M,LDA) &
          !$omp             private(MXRP,I,II,JJ)
          DO ip=0,mproc-1
             mxrp = sp8(ip)
             DO i=1,lr1
                ii = ip*lda + (i-1)*mxrp + 1
                jj = (i-1)*m + 1
                CALL zsctr_no_omp(mxrp,xf(ii),msp(1,ip+1),yf(jj))
             ENDDO
          ENDDO
       CASE(1)
          ! that seems to work better on the P
          !$omp parallel do shared(SP8,LR1,M,LDA) &
          !$omp             private(MXRP,I,II,JJ,K)
          DO ip=0,mproc-1
             mxrp = sp8(ip)
             DO i=1,lr1
                ii = ip*lda + (i-1)*mxrp
                jj = (i-1)*m
                DO k=1,mxrp
                   yf(jj+msp(k,ip+1)) = xf(ii+k)
                ENDDO
             ENDDO
          ENDDO
       END SELECT
    ENDIF
    !$    IF(HASNT_THREADED_COPIES) call omp_set_num_threads(max_threads)
    IF (HAS_LOW_LEVEL_TIMERS) CALL tihalt(procedureN,isub1)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE unpack_x2y
  ! ==================================================================
  SUBROUTINE phasen(f,kr1,kr2s,kr3s,n1u,n1o,nr2s,nr3s)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: kr1, kr2s, kr3s
    COMPLEX(real_8)                          :: f(kr1,kr2s,kr3s)
    INTEGER                                  :: n1u, n1o, nr2s, nr3s

    INTEGER                                  :: i, ii, ijk, isub, j, k
    REAL(real_8), DIMENSION(2)               :: pf = (/1._real_8,-1._real_8/)

    IF (HAS_LOW_LEVEL_TIMERS) CALL tiset('     PHASE',isub)
    !$omp parallel do default(none) __COLLAPSE2 &
    !$omp             private(K,J,I,II,IJK) &
    !$omp             shared(F,PF,NR3S,NR2S,N1U,N1O)
    DO k=1,nr3s
       DO j=1,nr2s
          DO i=n1u,n1o
             ii=i-n1u+1
             ijk=MOD(k+j+i+1,2)+1
             f(ii,j,k)=f(ii,j,k)*pf(ijk)
          ENDDO
       ENDDO
    ENDDO
    IF (HAS_LOW_LEVEL_TIMERS) CALL tihalt('     PHASE',isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE phasen
  ! ==================================================================
END MODULE fftutil_utils
