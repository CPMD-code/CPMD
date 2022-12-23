MODULE xcener_utils
  USE fft_maxfft,                      ONLY: maxfft
  USE func,                            ONLY: func1,&
                                             func2,&
                                             mfxcc_is_lyp,&
                                             mfxcc_is_pade,&
                                             mfxcc_is_skipped,&
                                             mfxcx_is_skipped,&
                                             mfxcx_is_slaterx
  USE functionals_utils,               ONLY: xc
  USE kinds,                           ONLY: real_8
  USE lsd_func_utils,                  ONLY: xc_lsd
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             parm,&
                                             spar
  USE tbxc,                            ONLY: eexc,&
                                             tabx,&
                                             toldcode,&
                                             vvxc
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: xcener

CONTAINS

  ! ==================================================================
  SUBROUTINE xcener(sxc,vxc,rhoval,rhoe,v)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! == EXC:     EXCHANGE AND CORRELATION ENERGY                     ==
    ! == VXC:     ENERGY FROM EXC POTENTIAL (Int{n exc(n)})           ==
    ! == RHOVAL:  VALENCE CHARGE DENSITY                              ==
    ! == RHOE:    TOTAL CHARGE DENSITY (IF TINLC CORE+VALENCE)        ==
    ! == V:        IN -> POTENTIAL IN REAL SPACE                      ==
    ! ==          OUT -> POTENTIAL + EXC POTENTIAL IN REAL SPACE      ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: sxc, vxc, rhoval(:,:), &
                                                rhoe(:,:)
    COMPLEX(real_8)                          :: v(:,:)

!rhoe(nnr1,clsd%nlsd)
!v(maxfft,clsd%nlsd)
! ==--------------------------------------------------------------==

    IF (toldcode) THEN
       CALL xcener_old(sxc,vxc,rhoval,rhoe,v)
    ELSE
       CALL xcener_new(sxc,vxc,rhoval,rhoe,v)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE xcener
  ! ==================================================================
  SUBROUTINE xcener_old(sxc,vxc,rhoval,rhoe,v)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! == EXC:     EXCHANGE AND CORRELATION ENERGY                     ==
    ! == VXC:     ENERGY FROM EXC POTENTIAL (Int{n exc(n)})           ==
    ! == RHOVAL:  VALENCE CHARGE DENSITY                              ==
    ! == RHOE:    TOTAL CHARGE DENSITY                                ==
    ! == V:        IN -> POTENTIAL IN REAL SPACE                      ==
    ! ==          OUT -> POTENTIAL + EXC POTENTIAL IN REAL SPACE      ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: sxc, vxc, rhoval(:,:), &
                                                rhoe(:,:)
    COMPLEX(real_8)                          :: v(:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'xcener_old'

    INTEGER                                  :: i, i1, ir, isub
    REAL(real_8)                             :: dd, ec, ee, eta, ex, exc, &
                                                exc1, ratio, roe, vc, vca, &
                                                vcb, vx, vxa, vxb, vxc1, &
                                                vxca, vxcb

    CALL tiset(procedureN,isub)
    sxc=0._real_8
    vxc=0._real_8
    IF (cntl%tlsd) CALL dcopy(2*fpar%nnr1,v(1,1),1,v(1,2),1)
    IF (func1%mfxcc /= mfxcc_is_skipped.OR.func1%mfxcx /= mfxcx_is_skipped) THEN
       IF (tabx%narray.LT.2.AND..NOT.cntl%tlsd) THEN
#ifdef __SR8000
          !poption parallel
          !poption tlocal(EX,EC,VX,VC,VXC1,EXC1)
#endif
          !$omp parallel do private(I,ROE,EX,EC,VX,VC,VXC1,EXC1) &
          !$omp reduction(+:SXC,VXC)
          DO i=1,fpar%nnr1
             roe=MAX(rhoe(i,1),0.0_real_8)
             CALL xc(roe,ex,ec,vx,vc)
             ! hybrid weights are now in subroutine XC
             vxc1=vx+vc
             exc1=ex+ec
             v(i,1)=v(i,1)+CMPLX(vxc1,0.0_real_8,kind=real_8)
             sxc=sxc+exc1*roe
             vxc=vxc+vxc1*MAX(rhoval(i,1),0.0_real_8)
          ENDDO
       ELSE
          IF (.NOT.cntl%tlsd) THEN
#ifdef __SR8000
             !poption parallel
             !poption tlocal(EX,EC,VX,VC,VXC1,EXC1)
             !poption tlocal(ROE)
             !voption pvfunc(3)
#endif
             !$omp parallel do private(I,ROE,RATIO,I1,DD,EE) &
             !$omp private(VXC1,EXC1,EX,VX,EC,VC) &
             !$omp reduction(+:SXC,VXC)
             DO i=1,fpar%nnr1
                roe=MAX(rhoe(i,1),0.0_real_8)
                ratio=(1._real_8/tabx%ddro)*roe
                i1=ratio
                IF (i1.LT.tabx%narray) THEN
                   dd=ratio-REAL(i1,kind=real_8)
                   ee=1._real_8-dd
                   vxc1=vvxc(i1)*ee+vvxc(i1+1)*dd
                   exc1=eexc(i1)*ee+eexc(i1+1)*dd
                ELSE
                   CALL xc(roe,ex,ec,vx,vc)
                   ! hybrid weights are now in subroutine XC
                   vxc1=vx+vc
                   exc1=ex+ec
                ENDIF
                v(i,1)=v(i,1)+CMPLX(vxc1,0.0_real_8,kind=real_8)
                sxc=sxc+exc1*roe
                vxc=vxc+vxc1*MAX(rhoval(i,1),0.0_real_8)
             ENDDO
          ELSE
             ! ==----------------------------------------------------------==
             ! ==  INPUT : RHOE(1..NNR1)  ALPHA + BETA DENSITY             ==
             ! ==          RHOE(NNR1+1..2*NNR1)  BETA DENSITY             ==
             ! ==          V(1..NNR1)   LOCAL POTENTIAL                    ==
             ! ==  OUTPUT  V(1..NNR1)   ALPHA POTENTIAL                    ==
             ! ==          V(NNR1+1..2*NNR1)   BETA POTENTIAL              ==
             ! ==----------------------------------------------------------==
#ifdef __SR8000
             !poption parallel
             !poption tlocal(EX,EC,VXA,VCA,VXB,VCB,VXCA,VXCB,EXC)
#endif
             !$omp parallel do private(IR,ROE,ETA,EX,EC,VXA,VCA,VXCA) &
             !$omp private(VXB,VCB,VXCB,EXC) &
             !$omp reduction(+:SXC,VXC)
             DO ir=1,fpar%nnr1
                roe = MAX(rhoe(ir,1),1.0e-15_real_8)
                eta = (rhoe(ir,1)-2._real_8*rhoe(ir,2))/roe
                IF (ABS(eta).GT.1._real_8) eta=SIGN(1.0_real_8,eta)
                CALL xc_lsd(roe,eta,ex,ec,vxa,vca,vxb,vcb)
                ! hybrid weights are now in subroutine xc_lsd
                vxca=vxa+vca
                vxcb=vxb+vcb
                exc=ex+ec
                v(ir,1)=v(ir,1)+CMPLX(vxca,0.0_real_8,kind=real_8)
                v(ir,2)=v(ir,2)+CMPLX(vxcb,0.0_real_8,kind=real_8)
                sxc=sxc+exc*roe
                ! .aa
                vxc=vxc+(vxca*(rhoval(ir,1)-rhoval(ir,2))+&
                     vxcb*rhoval(ir,2))
                ! .aa
             ENDDO
          ENDIF
       ENDIF
    ENDIF
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE xcener_old
  ! ==================================================================
  SUBROUTINE xcener_new(sxc,vxc,rhoval,rhoe,v)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! == EXC:     EXCHANGE AND CORRELATION ENERGY                     ==
    ! == VXC:     ENERGY FROM EXC POTENTIAL (Int{n exc(n)})           ==
    ! == RHOVAL:  VALENCE CHARGE DENSITY                              ==
    ! == RHOE:    TOTAL CHARGE DENSITY                                ==
    ! == V:        IN -> POTENTIAL IN REAL SPACE                      ==
    ! ==          OUT -> POTENTIAL + EXC POTENTIAL IN REAL SPACE      ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: sxc, vxc, rhoval(:,:), &
                                                rhoe(:,:)
    COMPLEX(real_8)                          :: v(:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'xcener_new'

    INTEGER                                  :: ir, isub, ntab
    REAL(real_8)                             :: flops

    CALL tiset(procedureN,isub)
    sxc=0._real_8
    vxc=0._real_8
    ntab=0
    IF (.NOT.cntl%tlsd) THEN
       IF (func1%mfxcc == mfxcc_is_skipped.AND.func1%mfxcx == mfxcx_is_slaterx) THEN
          !$omp parallel do private(IR) reduction(+:VXC)
          DO ir=1,fpar%nnr1
             vxc=vxc+(-ABS(rhoval(ir,1))*REAL(v(ir,1)))
          ENDDO
          CALL xalpu(func2%salpha,parm%nr1,spar%nr2s,spar%nr3s,fpar%kr1,fpar%kr2s,fpar%kr3s,&
               rhoe,v,sxc,flops)
          !$omp parallel do private(IR) reduction(+:VXC)
          DO ir=1,fpar%nnr1
             vxc=vxc+ABS(rhoval(ir,1))*REAL(v(ir,1))
          ENDDO
       ELSEIF (func1%mfxcc == mfxcc_is_lyp.AND.func1%mfxcx == mfxcx_is_slaterx) THEN
          !$omp parallel do private(IR) reduction(+:VXC)
          DO ir=1,fpar%nnr1
             vxc=vxc+(-ABS(rhoval(ir,1))*REAL(v(ir,1)))
          ENDDO
          CALL lypuu(parm%nr1,spar%nr2s,spar%nr3s,fpar%kr1,fpar%kr2s,fpar%kr3s,rhoe,v,sxc,flops)
          !$omp parallel do private(IR) reduction(+:VXC)
          DO ir=1,fpar%nnr1
             vxc=vxc+ABS(rhoval(ir,1))*REAL(v(ir,1))
          ENDDO
       ELSEIF (func1%mfxcc == mfxcc_is_pade) THEN
          !$omp parallel do private(IR) reduction(+:VXC)
          DO ir=1,fpar%nnr1
             vxc=vxc+(-ABS(rhoval(ir,1))*REAL(v(ir,1)))
          ENDDO
          CALL mikeu(parm%nr1,spar%nr2s,spar%nr3s,fpar%kr1,fpar%kr2s,fpar%kr3s,rhoe,v,sxc,flops)
          !$omp parallel do private(IR) reduction(+:VXC)
          DO ir=1,fpar%nnr1
             vxc=vxc+ABS(rhoval(ir,1))*REAL(v(ir,1))
          ENDDO
       ELSE
          CALL xcener_old(sxc,vxc,rhoval,rhoe,v)
       ENDIF
    ELSE
       ! ==------------------------------------------------------------==
       ! ==  INPUT : RHOE(1..NNR1)  ALPHA + BETA DENSITY               ==
       ! ==          RHOE(NNR1+1..2*NNR1)  BETA DENSITY               ==
       ! ==          V(1..NNR1)   LOCAL POTENTIAL                      ==
       ! ==  OUTPUT  V(1..NNR1,1)  ALPHA POTENTIAL                     ==
       ! ==          V(1..NNR1,2)   BETA POTENTIAL                     ==
       ! ==------------------------------------------------------------==
       CALL dcopy(2*fpar%nnr1,v(1,1),1,v(1,2),1)
       IF (func1%mfxcc == mfxcc_is_skipped.AND.func1%mfxcx == mfxcx_is_slaterx) THEN
          !$omp parallel do private(IR) reduction(+:VXC)
          DO ir=1,fpar%nnr1
             vxc=vxc+(-ABS(rhoval(ir,1))*REAL(v(ir,1)))
          ENDDO
          CALL xalpsp(func2%salpha,parm%nr1,spar%nr2s,spar%nr3s,fpar%kr1,fpar%kr2s,fpar%kr3s,&
               rhoe,v(1,1),v(1,2),sxc,flops)
          !$omp parallel do private(IR) reduction(+:VXC)
          DO ir=1,fpar%nnr1
             vxc=vxc+(ABS(rhoval(ir,1)-rhoval(ir,2))*REAL(v(ir,1))&
                  +ABS(rhoval(ir,2))*REAL(v(ir,2)))
          ENDDO
       ELSEIF (func1%mfxcc == mfxcc_is_lyp.AND.func1%mfxcx == mfxcx_is_slaterx) THEN
          !$omp parallel do private(IR) reduction(+:VXC)
          DO ir=1,fpar%nnr1
             vxc=vxc+(-ABS(rhoval(ir,1))*REAL(v(ir,1)))
          ENDDO
          CALL slypsp(parm%nr1,spar%nr2s,spar%nr3s,fpar%kr1,fpar%kr2s,fpar%kr3s,rhoe,v(1,1),v(1,2),&
               func2%salpha,sxc,flops)
          !$omp parallel do private(IR) reduction(+:VXC)
          DO ir=1,fpar%nnr1
             vxc=vxc+(ABS(rhoval(ir,1)-rhoval(ir,2))*REAL(v(ir,1))&
                  +ABS(rhoval(ir,2))*REAL(v(ir,2)))
          ENDDO
       ELSEIF (func1%mfxcc == mfxcc_is_pade) THEN
          !$omp parallel do private(IR) reduction(+:VXC)
          DO ir=1,fpar%nnr1
             vxc=vxc+(-ABS(rhoval(ir,1))*REAL(v(ir,1)))
          ENDDO
          CALL mikesp(parm%nr1,spar%nr2s,spar%nr3s,fpar%kr1,fpar%kr2s,fpar%kr3s,rhoe,v(1,1),v(1,2),&
               sxc,flops)
          !$omp parallel do private(IR) reduction(+:VXC)
          DO ir=1,fpar%nnr1
             vxc=vxc+(ABS(rhoval(ir,1)-rhoval(ir,2))*REAL(v(ir,1))&
                  +ABS(rhoval(ir,2))*REAL(v(ir,2)))
          ENDDO
       ELSE
          CALL xcener_old(sxc,vxc,rhoval,rhoe,v)
       ENDIF
    ENDIF
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE xcener_new
  ! ==================================================================

END MODULE xcener_utils
