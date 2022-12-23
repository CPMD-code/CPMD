MODULE epr_dv0_utils
  USE cnst,                            ONLY: uimag
  USE cppt,                            ONLY: gk,&
                                             indz,&
                                             nzh
  USE fft_maxfft,                      ONLY: maxfft
  USE fftmain_utils,                   ONLY: fwfftn,&
                                             invfftn
  USE fftnew_utils,                    ONLY: setfftn
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: parai
  USE system,                          ONLY: fpar,&
                                             ncpw,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: zgthr
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: dv0

CONTAINS

  ! ==================================================================
  SUBROUTINE dv0(v0,vt,grad,vt2)
    ! ==--------------------------------------------------------------==
    ! ==   CALCULATION OF |nabla.V0|                                  ==
    ! ==   ON INPUT : V0 : V0 IN REAL SPACE                           ==
    ! ==              VT   : UNDEFINED                                ==
    ! ==              GRAD : UNDEFINED                                ==
    ! ==              VT2  : UNDeFINED                                ==
    ! ==   ON OUTPUT: V0  : V0 IN REAL SPACE (SMOOTH)                 ==
    ! ==              VT   : UNDEFINED                                ==
    ! ==              GRAD : GRADIENT OF V0 IN REAL SPACE             ==
    ! ==              VT2  : V0 IN G SPACE                            ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: v0(fpar%nnr1)
    COMPLEX(real_8)                          :: vt(maxfft)
    REAL(real_8)                             :: grad(fpar%nnr1,3)
    COMPLEX(real_8)                          :: vt2(ncpw%nhg)

    INTEGER                                  :: ig, ir, isub

! ==--------------------------------------------------------------==
! ==  TRANSFORM DENSITY TO G SPACE                                ==
! ==--------------------------------------------------------------==

    CALL tiset('    DIFFV0',isub)
    CALL setfftn(0)
    !$omp parallel do private(IR)
    DO ir=1,fpar%nnr1
       vt(ir) = CMPLX(v0(ir),0.0_real_8,kind=real_8)
    ENDDO
    CALL  fwfftn(vt,.FALSE.,parai%allgrp)
    CALL zgthr(ncpw%nhg,vt,vt2,nzh)
    ! ==--------------------------------------------------------------==
    ! ==  FFT OF RHO AND NABLA(X)*RHOE                                ==
    ! ==--------------------------------------------------------------==
    CALL zeroing(vt)!,maxfft)
    !$omp parallel do private(IG)
#ifdef __NEC
    !CDIR NODEP
#endif
#ifdef __SR8000
    !poption parallel
#endif
#ifdef _vpp_
    !OCL NOVREC(VT)
#endif
    DO ig=1,ncpw%nhg
       vt(nzh(ig))=vt2(ig)-parm%tpiba*gk(1,ig)*vt2(ig)
       vt(indz(ig))=CONJG(vt2(ig)+parm%tpiba*gk(1,ig)*vt2(ig))
    ENDDO
    CALL  invfftn(vt,.FALSE.,parai%allgrp)
    !$omp parallel do private(IR)
    DO ir=1,fpar%nnr1
       v0(ir)=REAL(vt(ir))
       grad(ir,1)=AIMAG(vt(ir))
    ENDDO
    ! ==--------------------------------------------------------------==
    ! ==  FFT OF NABLA(Y)*RHO AND NABLA(Z)*RHOE                       ==
    ! ==--------------------------------------------------------------==
    CALL zeroing(vt)!,maxfft)
#ifdef __ES  
    !CDIR NODEP(VT)
    DO ig=1,nhg
       vt(nzh(ig))=parm%tpiba*(uimag*gk(2,ig)-gk(3,ig))*vt2(ig)
       vt(indz(ig))=parm%tpiba*(-uimag*gk(2,ig)+gk(3,ig))*CONJG(vt2(ig))
    ENDDO
#else 
    !$omp parallel do private(IG)
#ifdef __SR8000
    !poption parallel
#endif
#ifdef _vpp_
    !OCL NOVREC(VT)
#endif
    DO ig=1,ncpw%nhg
       vt(nzh(ig))=parm%tpiba*(uimag*gk(2,ig)-gk(3,ig))*vt2(ig)
       vt(indz(ig))=parm%tpiba*(-uimag*gk(2,ig)+gk(3,ig))*CONJG(vt2(ig))
    ENDDO
#endif 
    CALL  invfftn(vt,.FALSE.,parai%allgrp)
    !$omp parallel do private(IR)
    DO ir=1,fpar%nnr1
       grad(ir,2)=REAL(vt(ir))
       grad(ir,3)=AIMAG(vt(ir))
    ENDDO
    CALL tihalt('    DIFFV0',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE dv0
  ! ==================================================================

END MODULE epr_dv0_utils
