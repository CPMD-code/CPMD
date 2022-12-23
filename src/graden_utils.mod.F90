#include "cpmd_global.h"

MODULE graden_utils
  USE cnst,                            ONLY: uimag
  USE cppt,                            ONLY: gk,&
                                             hg,&
                                             indz,&
                                             nzh
  USE fft_maxfft,                      ONLY: maxfftn
  USE fftmain_utils,                   ONLY: fwfftn,&
                                             invfftn
  USE fftnew_utils,                    ONLY: setfftn
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_max
  USE nvtx_utils
  USE parac,                           ONLY: parai
  USE state_utils,                     ONLY: copy_to_re
  USE system,                          ONLY: cntl,&
                                             cntr,&
                                             fpar,&
                                             ncpw,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: zgthr
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: graden

CONTAINS

  ! ==================================================================
  SUBROUTINE graden(rhoe,v,grad,vtmp)
    ! ==--------------------------------------------------------------==
    ! ==   SMOOTHING OF THE DENSITY AND CALCULATION OF |nabla.RHO|    ==
    ! ==   ON INPUT : RHOE : DENSITY IN REAL SPACE                    ==
    ! ==              V    : UNDEFINED                                ==
    ! ==              GRAD : UNDEFINED                                ==
    ! ==              VTMP : UNDEFINED                                ==
    ! ==   ON OUTPUT: RHOE : DENSITY IN REAL SPACE (SMOOTH)           ==
    ! ==              V    : UNDEFINED                                ==
    ! ==              GRAD : (GRADIENT OF RHO)^2 IN REAL SPACE        ==
    ! ==              VTMP : DENSITY IN G SPACE (SMOOTH)              ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rhoe(fpar%nnr1)
    COMPLEX(real_8)                          :: v(maxfftn)
    REAL(real_8)                             :: grad(fpar%nnr1,4)
    COMPLEX(real_8)                          :: vtmp(ncpw%nhg)

    CHARACTER(*), PARAMETER                  :: procedureN = 'graden'

    INTEGER                                  :: ig, ir, isub
    REAL(real_8)                             :: eg, gcs, gmax, smfac

! Variables
! ==--------------------------------------------------------------==
! ==  TRANSFORM DENSITY TO G SPACE                                ==
! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    __NVTX_TIMER_START ( procedureN )

    CALL setfftn(0)

    CALL zeroing ( v )
    CALL copy_to_re ( fpar%nnr1, rhoe, v )
!!$omp parallel do private(IR)
    !DO ir=1,fpar%nnr1
    !   v(ir) = CMPLX(rhoe(ir),0.0_real_8,kind=real_8)
    !ENDDO
    CALL  fwfftn(v,.FALSE.,parai%allgrp)
    ! ==--------------------------------------------------------------==
    ! ==  SMOOTHING                                                   ==
    ! ==--------------------------------------------------------------==
    IF (cntl%tsmooth) THEN
       gmax=hg(ncpw%nhg)
       CALL mp_max(gmax,parai%allgrp)
       gcs=cntr%smf*gmax
       !$omp parallel do private(IG,EG,SMFAC)
       DO ig=1,ncpw%nhg
          eg=(hg(ig)-gcs)/(cntr%sdelta*gmax)
          smfac=1.0_real_8/(1.0_real_8+EXP(eg))
          vtmp(ig)=v(nzh(ig))*smfac
       ENDDO
    ELSE
       CALL zgthr(ncpw%nhg,v,vtmp,nzh)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ==  FFT OF RHO AND NABLA(X)*RHOE                                ==
    ! ==--------------------------------------------------------------==
    CALL zeroing(v)!,maxfft)
    !ocl novrec(v)
    !$omp parallel do private(IG)
#ifdef __NEC
    !CDIR NODEP
#endif
#ifdef __SR8000
    !poption parallel
#endif
    DO ig=1,ncpw%nhg
       v(nzh(ig))=vtmp(ig)-parm%tpiba*gk(1,ig)*vtmp(ig)
       v(indz(ig))=CONJG(vtmp(ig)+parm%tpiba*gk(1,ig)*vtmp(ig))
    ENDDO
    CALL  invfftn(v,.FALSE.,parai%allgrp)
    !$omp parallel do private(IR)
    DO ir=1,fpar%nnr1
       rhoe(ir)=REAL(v(ir))
       grad(ir,1)=AIMAG(v(ir))*AIMAG(v(ir))
       grad(ir,2)=AIMAG(v(ir))
    ENDDO
    ! ==--------------------------------------------------------------==
    ! ==  FFT OF NABLA(Y)*RHO AND NABLA(Z)*RHOE                       ==
    ! ==--------------------------------------------------------------==
    CALL zeroing(v)!,maxfft)
#ifdef __ES  
    !CDIR NODEP(V)
    DO ig=1,nhg
       v(nzh(ig))=parm%tpiba*(uimag*gk(2,ig)-gk(3,ig))*vtmp(ig)
       v(indz(ig))=parm%tpiba*(-uimag*gk(2,ig)+gk(3,ig))*CONJG(vtmp(ig))
    ENDDO
#else 
    !ocl novrec(v)
    !$omp parallel do private(IG)
#ifdef __SR8000
    !poption parallel
#endif
    DO ig=1,ncpw%nhg
       v(nzh(ig))=parm%tpiba*(uimag*gk(2,ig)-gk(3,ig))*vtmp(ig)
       v(indz(ig))=parm%tpiba*(-uimag*gk(2,ig)+gk(3,ig))*CONJG(vtmp(ig))
    ENDDO
#endif 
    CALL  invfftn(v,.FALSE.,parai%allgrp)
    !$omp parallel do private(IR)
    DO ir=1,fpar%nnr1
       grad(ir,1)=grad(ir,1)+REAL(v(ir)*CONJG(v(ir)))
       grad(ir,3)=REAL(v(ir))
       grad(ir,4)=AIMAG(v(ir))
    ENDDO

    __NVTX_TIMER_STOP
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE graden
  ! ==================================================================

END MODULE graden_utils
