MODULE vbeta_utils
  USE cnst,                            ONLY: uimag
  USE cppt,                            ONLY: indzs,&
                                             nzhs,&
                                             twnl
  USE fft_maxfft,                      ONLY: maxfftn
  USE fftmain_utils,                   ONLY: fwfftn,&
                                             invfftn
  USE fint,                            ONLY: cnl,&
                                             nlptr
  USE geq0mod,                         ONLY: geq0
  USE kinds,                           ONLY: real_8
  USE kpnt,                            ONLY: eigkr
  USE kpts,                            ONLY: tkpts
  USE nlps,                            ONLY: nghtol,&
                                             nlm
  USE parac,                           ONLY: parai
  USE sfac,                            ONLY: eigr
  USE system,                          ONLY: fpar,&
                                             ncpw
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: vbeta

CONTAINS

  ! ==================================================================
  SUBROUTINE vbeta(rhoe,psi,ikind)
    ! ==--------------------------------------------------------------==
    ! ==  APPLIES e^(-beta/(2P)Vlocal) to NL-PP projectors            ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rhoe(fpar%nnr1)
    COMPLEX(real_8)                          :: psi(maxfftn)
    INTEGER                                  :: ikind

    COMPLEX(real_8)                          :: cp1, cp2, fm, fp
    INTEGER                                  :: i, i1, i2, ifft, ig, ir, is1, &
                                                is2, isa1, isa2, isub, iv1, &
                                                iv2, njump

    IF (nlm.EQ.0) RETURN
    CALL tiset('     VBETA',isub)
    ! ..prepare the potential
    !$omp parallel do private(IR)
    DO ir = 1 , fpar%nnr1
       rhoe(ir)=SQRT(rhoe(ir))
    ENDDO
    ifft=1
    ! kpt  We do not use the trick for FFT.
    IF (tkpts%tkpnt) THEN
       njump=1
    ELSE
       njump=2
    ENDIF
    DO i = 1, nlm, njump
       i1 = i
       i2 = i+1
       CALL zeroing(psi)!,maxfft)
       IF ((i2.LE.nlm).AND.(.NOT.tkpts%tkpnt)) THEN
          is1=nlptr(1,i1)
          isa1=nlptr(2,i1)
          iv1=nlptr(3,i1)
          is2=nlptr(1,i2)
          isa2=nlptr(2,i2)
          iv2=nlptr(3,i2)
          cp1=(0.0_real_8,-1.0_real_8)**nghtol(iv1,is1)
          cp2=(0.0_real_8,-1.0_real_8)**nghtol(iv2,is2)
          !CDIR NODEP
#ifdef __SR8000
          !poption parallel
#endif
#if defined(__VECTOR)
          !$omp parallel do private(IG,FP,FM)
#else
          !$omp parallel do private(IG,FP,FM) schedule(static)
#endif
          DO ig = 1,ncpw%ngw
             fp=cp1*twnl(ig,iv1,is1,ikind)*eigr(ig,isa1,1)
             fm=cp2*twnl(ig,iv2,is2,ikind)*eigr(ig,isa2,1)
             psi(nzhs(ig)) = fp+uimag*fm
             psi(indzs(ig)) = CONJG(fp)+uimag*CONJG(fm)
          ENDDO
          IF (geq0) psi(nzhs(1))&
               =cp1*twnl(1,iv1,is1,ikind)*eigr(1,isa1,1)+&
               uimag*cp2*twnl(1,iv2,is2,ikind)*eigr(1,isa2,1)
       ELSE
          is1=nlptr(1,i1)
          isa1=nlptr(2,i1)
          iv1=nlptr(3,i1)
          cp1=(0.0_real_8,-1.0_real_8)**nghtol(iv1,is1)
          IF (tkpts%tkpnt) THEN
             !CDIR NODEP
#ifdef __SR8000
             !poption parallel
#endif
#if defined(__VECTOR)
             !$omp parallel do private(IG)
#else
             !$omp parallel do private(IG) schedule(static)
#endif
             DO ig = 1,ncpw%ngw
                psi(nzhs(ig))  = cp1*twnl(ig    ,iv1,is1,ikind)*&
                     eigkr(ig    ,isa1,ikind)
                psi(indzs(ig)) = cp1*twnl(ig+ncpw%ngw,iv1,is1,ikind)*&
                     eigkr(ig+ncpw%ngw,isa1,ikind)
             ENDDO
             IF (geq0) psi(nzhs(1))&
                  =cp1*twnl(1,iv1,is1,ikind)*eigkr(1,isa1,ikind)
          ELSE
             !CDIR NODEP
#ifdef __SR8000
             !poption parallel
#endif
#if defined(__VECTOR)
             !$omp parallel do private(IG,FP)
#else
             !$omp parallel do private(IG,FP) schedule(static)
#endif
             DO ig = 1,ncpw%ngw
                fp=cp1*twnl(ig,iv1,is1,ikind)*eigr(ig,isa1,1)
                psi(nzhs(ig)) = fp
                psi(indzs(ig)) = CONJG(fp)
             ENDDO
             IF (geq0) psi(nzhs(1))&
                  =cp1*twnl(1,iv1,is1,ikind)*eigr(1,isa1,1)
          ENDIF
       ENDIF
       ! ==------------------------------------------------------------==
       ! == INVERSE FFT PSI TO OBTAIN PSI ON REAL SPACE MESH           ==
       ! ==------------------------------------------------------------==
       CALL  invfftn(psi,.TRUE.,parai%allgrp)
       ! ==------------------------------------------------------------==
       ! == APPLY THE EXPONENTIATED POTENTIAL TO PSI                   ==
       ! == APPLY e^(-b/2P V) IN REAL SPACE, AND RETURN ANSWER IN PSI. ==
       ! == ARRAY POT CONTAINS e^(-b/P V(r)).                          ==
       ! ==------------------------------------------------------------==
       !$omp parallel do private(IR)
       DO ir = 1 , fpar%nnr1
          psi(ir) = psi(ir)*rhoe(ir)
       ENDDO
       ! ==------------------------------------------------------------==
       ! == TRANSFORM BACK TO FOURIER SPACE                            ==
       ! ==------------------------------------------------------------==
       CALL  fwfftn(psi,.TRUE.,parai%allgrp)
       IF ((i2.LE.nlm).AND.(.NOT.tkpts%tkpnt)) THEN
          !CDIR NODEP
#ifdef __SR8000
          !poption parallel
#endif
#if defined(__VECTOR)
          !$omp parallel do private(IG,FP,FM)
#else
          !$omp parallel do private(IG,FP,FM) schedule(static)
#endif
          DO ig = 1 , ncpw%ngw
             fp=(psi(nzhs(ig))+psi(indzs(ig)))*0.5_real_8
             fm=(psi(nzhs(ig))-psi(indzs(ig)))*0.5_real_8
             cnl(ig,i1) = CMPLX(REAL(fp),AIMAG(fm),kind=real_8)
             cnl(ig,i2) = CMPLX(AIMAG(fp),-REAL(fm),kind=real_8)
          ENDDO
          IF (geq0) THEN
             cnl(1,i1) = CMPLX(REAL(cnl(1,i1)),0._real_8,kind=real_8)
             cnl(1,i2) = CMPLX(REAL(cnl(1,i2)),0._real_8,kind=real_8)
          ENDIF
       ELSE
          IF (tkpts%tkpnt) THEN
#if defined(__VECTOR)
             !$omp parallel do private(IG)
#else
             !$omp parallel do private(IG) schedule(static)
#endif
             DO ig = 1, ncpw%ngw
                cnl(ig    ,i1) = psi(nzhs(ig))
                cnl(ig+ncpw%ngw,i1) = psi(indzs(ig))
             ENDDO
             IF (geq0) cnl(1+ncpw%ngw,i1) = CMPLX(0._real_8,0._real_8,kind=real_8)
          ELSE
#if defined(__VECTOR)
             !$omp parallel do private(IG,FP,FM)
#else
             !$omp parallel do private(IG,FP,FM) schedule(static)
#endif
             DO ig = 1, ncpw%ngw
                fp=(psi(nzhs(ig))+psi(indzs(ig)))*0.5_real_8
                fm=(psi(nzhs(ig))-psi(indzs(ig)))*0.5_real_8
                cnl(ig,i1) = CMPLX(REAL(fp),AIMAG(fm),kind=real_8)
             ENDDO
             IF (geq0) cnl(1,i1) = CMPLX(REAL(cnl(1,i1)),0._real_8,kind=real_8)
          ENDIF
       ENDIF
    ENDDO
    ! ..Return to the old potential
    !$omp parallel do private(IR)
    DO ir = 1 , fpar%nnr1
       rhoe(ir) = rhoe(ir)*rhoe(ir)
    ENDDO
    CALL tihalt('     VBETA',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE vbeta
  ! ==================================================================

END MODULE vbeta_utils
