
#include "cpmd_global.h"

MODULE gcener_utils
  USE cnst,                            ONLY: uimag
  USE cp_xc_driver,                    ONLY: cp_xc_compute
  USE cppt,                            ONLY: gk,&
                                             hg,&
                                             indz,&
                                             nzh
  USE error_handling,                  ONLY: stopgm
  USE fftmain_utils,                   ONLY: fwfftn,&
                                             invfftn
  USE fftnew_utils,                    ONLY: setfftn
  USE func,                            ONLY: &
       func1, func2, func3, mgcc_is_ggac, mgcc_is_lyp, mgcc_is_pbec, &
       mgcc_is_perdew86, mgcc_is_skipped, mgcsrx_is_skipped, mgcx_is_becke88, &
       mgcx_is_ggax, mgcx_is_pbex, mgcx_is_revpbex
  USE functionals_utils,               ONLY: gcxc,&
                                             pbexsr,&
                                             pbexsr_lsd
  USE gcxctbl_utils,                   ONLY: &
       gcgga, gcpbe, gcrevpbe, gcspbe, gcsrevpbe, gcsxlyp, gcsxonly, gcsxp86, &
       gcxlyp, gcxonly, gcxp86
  USE kinds,                           ONLY: real_8
  USE lsd_func_utils,                  ONLY: gc_lsd
  USE metafun_utils,                   ONLY: taufun,&
                                             taufuns
  USE mp_interface,                    ONLY: mp_max
  USE nvtx_utils
  USE parac,                           ONLY: parai
  USE strs,                            ONLY: alpha,&
                                             beta,&
                                             degc
  USE system,                          ONLY: cntl,&
                                             cntr,&
                                             fpar,&
                                             ncpw,&
                                             parm,&
                                             spar
  USE tauf,                            ONLY: tau,&
                                             vtau
  USE tbxc,                            ONLY: toldcode
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: gcener
  PUBLIC :: gclsd
  !public :: gcold
  !public :: gclsdold

CONTAINS

  ! ==================================================================
  SUBROUTINE gcener(sgcx,sgcc,rhoe,v,vtemp,vtmp,grad,tstress)
    ! ==--------------------------------------------------------------==
    ! ==   CALCULATION OF THE GRADIENT CORRECTION TO THE EXCHANGE     ==
    ! ==   AND CORRELATION FUNCTIONAL                                 ==
    ! ==                                                              ==
    ! ==   ON INPUT : RHOE : DENSITY IN REAL SPACE                    ==
    ! ==              V    : UNDEFINED                                ==
    ! ==              GRAD : (NABLA*RHOE)^2,RHOx,RHOy,RHOz            ==
    ! ==              VTEMP: POTENTIAL IN G SPACE                     ==
    ! ==              VTMP : DENSITY IN G SPACE                       ==
    ! ==   ON OUTPUT: RHOE : DENSITY IN REAL SPACE                    ==
    ! ==              V    : TOTAL POTENTIAL IN REAL SPACE            ==
    ! ==              GRAD : (NABLA*RHOE)^2,RHOx,RHOy,RHOz            ==
    ! ==              VTEMP: UNDEFINED                                ==
    ! ==              VTMP : UNDEFINED                                ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: sgcx, sgcc, rhoe(:,:)
    COMPLEX(real_8)                          :: v(:,:), vtemp(*)
    REAL(real_8), TARGET __CONTIGUOUS        :: vtmp(:)
    REAL(real_8)                             :: grad(:,:)
    LOGICAL                                  :: tstress, is_present = .false.

    CHARACTER(*), PARAMETER                  :: procedureN = 'gcener'

    COMPLEX(real_8)                          :: fa, fb, vg1, vg2, vin, vnz
    INTEGER                                  :: i1, i2, i3, ia, ib, ierr, ig, &
                                                ii3, ir, isub, kk
    REAL(real_8)                             :: eg, flops, gcs, gmax, smfac, &
                                                vfac, vfac2, vfac3
    REAL(real_8), DIMENSION(:, :), &
      POINTER __CONTIGUOUS                   :: vtau_2d_p, vtmp_2d_p

!
! Temporary fix if tau is not allocated (to be nicified in the future)
!

    IF (.NOT.cntl%tgcx .AND. .NOT.cntl%tgcc) RETURN
    CALL tiset(procedureN,isub)
    __NVTX_TIMER_START ( procedureN )

    CALL setfftn(0)
    ! ==--------------------------------------------------------------==
    ! ==  V1 = dF/dn stored in V(*)                                   ==
    ! ==  V2 = (1/Nabla.n) dF/d(Nabla.n) stored in VTMP(*)            ==
    ! ==--------------------------------------------------------------==
    CALL zeroing(v(:,1))!,nnr1)
    CALL zeroing(vtmp)!,nnr1)
    sgcx = 0.0_real_8
    sgcc = 0.0_real_8
    flops=0.0_real_8

    IF( cntl%use_xc_driver ) THEN
       vtmp_2d_p(1:SIZE(vtmp,1),1:1) => vtmp
       !
       ! Temporary fix for metafunctionals (if tau is not needed, allocate
       ! one array for tau and vtau, rather than two)
       !
       IF (cntl%ttau) THEN
          vtau_2d_p(1:SIZE(vtau,1),1:1) => vtau(:,1)
       ELSE 
          IF (ALLOCATED(tau)) THEN
             is_present = .true.
          ELSE
          ALLOCATE(tau(fpar%nnr1,1), stat=ierr)
          IF (ierr /= 0) CALL stopgm(procedureN,'allocation problem',& 
               __LINE__,__FILE__)
          ENDIF
          CALL zeroing(tau)
          vtau_2d_p(1:fpar%nnr1,1:1) => tau(:,1)
       END IF
       !
       ! End of temporary fix
       !
       CALL cp_xc_compute( fpar%nnr1, rhoe, grad, tau, &
            sgcx, sgcc, v, vtmp_2d_p, vtau_2d_p )
       !
       ! Temporary fix for metafunctionals
       !
       IF (ASSOCIATED(vtau_2d_p)) NULLIFY(vtau_2d_p)
       IF (.NOT. cntl%ttau) THEN
          IF (.NOT. is_present) THEN
          DEALLOCATE(tau, stat=ierr)
          IF (ierr /= 0) CALL stopgm(procedureN,'deallocation problem',& 
               __LINE__,__FILE__)
          ENDIF
       END IF
       !
       ! End of temporary fix
       !
    ELSE
       IF (toldcode) THEN
          CALL gcold(sgcx,sgcc,rhoe(:,1),v(:,1),vtmp,grad,flops)
       ELSE
          IF (func1%mgcx == mgcx_is_becke88) THEN
             IF (func1%mgcc == mgcc_is_skipped) THEN
                CALL gcxonly(func2%bbeta,sgcx,sgcc,rhoe,v,vtmp,grad,flops)
             ELSEIF (func1%mgcc == mgcc_is_perdew86) THEN
                CALL gcxp86(func2%bbeta,sgcx,sgcc,rhoe,v,vtmp,grad,flops)
             ELSEIF (func1%mgcc == mgcc_is_lyp) THEN
                CALL gcxlyp(func2%bbeta,sgcx,sgcc,rhoe,v,vtmp,grad,flops)
             ELSE
                CALL stopgm(procedureN,'NOT IMPLEMENTED',& 
                     __LINE__,__FILE__)
             ENDIF
          ELSEIF (func1%mgcx == mgcx_is_ggax.AND.func1%mgcc == mgcc_is_ggac) THEN
             CALL gcgga(sgcx,sgcc,rhoe,v,vtmp,grad,flops)
          ELSEIF (func1%mgcx == mgcx_is_pbex.AND.func1%mgcc == mgcc_is_pbec) THEN
             CALL gcpbe(sgcx,sgcc,rhoe,v,vtmp,grad,flops)
          ELSEIF (func1%mgcx == mgcx_is_revpbex.AND.func1%mgcc == mgcc_is_pbec) THEN
             CALL gcrevpbe(sgcx,sgcc,rhoe,v,vtmp,grad,flops)
          ELSE
             CALL gcold(sgcx,sgcc,rhoe(:,1),v(:,1),vtmp,grad,flops)
          ENDIF
          ! Set to 0.0 VTMP(1:KR1,1:KR2S,KR3S),VTMP(1:NR1,KR2S,:), and
          ! VTMP(KR1,1:KR2S,1:KR3S).
          DO i3=spar%nr3s+1,fpar%kr3s
             !          CALL azzero(vtmp((i3-1)*kr1*kr2s+1),kr1*kr2s)
             CALL zeroing(vtmp((i3-1)*fpar%kr1*fpar%kr2s+1:i3*fpar%kr1*fpar%kr2s) )! ,kr1*kr2s)
          ENDDO
          DO i3=1,spar%nr3s
             DO i2=spar%nr2s+1,fpar%kr2s
                !             CALL azzero(vtmp((i3-1)*kr1*kr2s+(i2-1)*kr1+1),kr1)
                CALL zeroing(vtmp((i3-1)*fpar%kr1*fpar%kr2s+(i2-1)*fpar%kr1+1:(i3-1)*fpar%kr1*fpar%kr2s+i2*fpar%kr1))!,kr1)
             ENDDO
             ii3 = (i3-1)*fpar%kr1*fpar%kr2s
             DO i2=1,spar%nr2s
                DO i1=parm%nr1+1,fpar%kr1
                   vtmp(ii3+(i2-1)*fpar%kr1+i1)=0._real_8
                ENDDO
             ENDDO
          ENDDO
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (tstress) THEN
       ! Contribution of dF/d(delta.n) to the stress tensor
       !$omp parallel do private(KK,IR,IA,IB)
       DO kk=1,6
          ia=alpha(kk)+1
          ib=beta(kk)+1
          degc(kk)=0.0_real_8
          DO ir=1,fpar%nnr1
             degc(kk)=degc(kk)+vtmp(ir)*grad(ir,ia)*grad(ir,ib)
          ENDDO
          degc(kk)=degc(kk)*parm%omega/(spar%nr1s*spar%nr2s*spar%nr3s)
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    ! FFT of V1 and V2*RHOx to G-Space
    !$omp parallel do private(IR)
    DO ir=1,fpar%nnr1
       v(ir,1) = v(ir,1) + CMPLX(0.0_real_8,vtmp(ir)*grad(ir,2),kind=real_8)
    ENDDO
    CALL  fwfftn(v(:,1),.FALSE.,parai%allgrp)
    gcs=cntr%smf*hg(ncpw%nhg)
    IF (cntl%tsmooth) THEN
       gmax=hg(ncpw%nhg)
       CALL mp_max(gmax,parai%allgrp)
       gcs=cntr%smf*gmax
       DO ig=1,ncpw%nhg
          eg = (hg(ig)-gcs)/(cntr%sdelta*gmax)
          smfac=1.0_real_8/(1.0_real_8+EXP(eg))
          vfac=smfac*parm%tpiba*gk(1,ig)
          ! !        FA  = V(NZH(IG)) + V(INDZ(IG))
          ! !        FB  = V(NZH(IG)) - V(INDZ(IG))
          vnz = v(nzh(ig),1)
          vin = v(indz(ig),1)
          fa  = vnz + vin
          fb  = vnz - vin
          vg1 = 0.5_real_8*CMPLX(REAL(fa),AIMAG(fb),kind=real_8)
          vg2 = 0.5_real_8*CMPLX(AIMAG(fa),-REAL(fb),kind=real_8)
          vtemp(ig) = vtemp(ig) + smfac*vg1 - vfac*uimag*vg2
       ENDDO
    ELSE
       !$omp parallel do private(IG,VFAC,vnz,vin,FA,FB,VG1,VG2)
       DO ig=1,ncpw%nhg
          vfac=parm%tpiba*gk(1,ig)
          ! !        FA  = V(NZH(IG)) + V(INDZ(IG))
          ! !        FB  = V(NZH(IG)) - V(INDZ(IG))
          vnz = v(nzh(ig),1)
          vin = v(indz(ig),1)
          fa  = vnz + vin
          fb  = vnz - vin
          vg1 = 0.5_real_8*CMPLX(REAL(fa),AIMAG(fb),kind=real_8)
          vg2 = 0.5_real_8*CMPLX(AIMAG(fa),-REAL(fb),kind=real_8)
          vtemp(ig) = vtemp(ig) + vg1 - vfac*uimag*vg2
       ENDDO
    ENDIF
    ! FFT of V2*RHOy and V2*RHOz to G-Space
    CALL zeroing(v(:,1))!,nnr1)
    !$omp parallel do private(IR)
    DO ir=1,fpar%nnr1
       v(ir,1) = CMPLX(vtmp(ir)*grad(ir,3),vtmp(ir)*grad(ir,4),kind=real_8)
    ENDDO
    CALL  fwfftn(v(:,1),.FALSE.,parai%allgrp)
    gcs=cntr%smf*hg(ncpw%nhg)
    IF (cntl%tsmooth) THEN
       gmax=hg(ncpw%nhg)
       CALL mp_max(gmax,parai%allgrp)
       gcs=cntr%smf*gmax
       DO ig=1,ncpw%nhg
          eg = (hg(ig)-gcs)/(cntr%sdelta*gmax)
          smfac=1.0_real_8/(1.0_real_8+EXP(eg))
          vfac2=smfac*parm%tpiba*gk(2,ig)
          vfac3=smfac*parm%tpiba*gk(3,ig)
          ! !        FA  = V(NZH(IG)) + V(INDZ(IG))
          ! !        FB  = V(NZH(IG)) - V(INDZ(IG))
          vnz = v(nzh(ig),1)
          vin = v(indz(ig),1)
          fa  = vnz + vin
          fb  = vnz - vin
          vg1 = 0.5_real_8*CMPLX(REAL(fa),AIMAG(fb),kind=real_8)
          vg2 = 0.5_real_8*CMPLX(AIMAG(fa),-REAL(fb),kind=real_8)
          vtemp(ig) = vtemp(ig) - vfac2*uimag*vg1 - vfac3*uimag*vg2
       ENDDO
    ELSE
       !$omp parallel do private(IG,VFAC2,VFAC3,vnz,vin,FA,FB,VG1,VG2)
       DO ig=1,ncpw%nhg
          vfac2=parm%tpiba*gk(2,ig)
          vfac3=parm%tpiba*gk(3,ig)
          ! !        FA  = V(NZH(IG)) + V(INDZ(IG))
          ! !        FB  = V(NZH(IG)) - V(INDZ(IG))
          vnz = v(nzh(ig),1)
          vin = v(indz(ig),1)
          fa  = vnz + vin
          fb  = vnz - vin
          vg1 = 0.5_real_8*CMPLX(REAL(fa),AIMAG(fb),kind=real_8)
          vg2 = 0.5_real_8*CMPLX(AIMAG(fa),-REAL(fb),kind=real_8)
          vtemp(ig) = vtemp(ig) - vfac2*uimag*vg1 - vfac3*uimag*vg2
       ENDDO
    ENDIF
    CALL zeroing(v(:,1))!,maxfft)
    !CDIR NODEP
    !ocl novrec
    !$omp parallel do private(IG)
    DO ig=1,ncpw%nhg
       v(nzh(ig),1) = vtemp(ig)
       v(indz(ig),1) = CONJG(vtemp(ig))
    ENDDO
    CALL invfftn(v(:,1),.FALSE.,parai%allgrp)

    __NVTX_TIMER_STOP
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE gcener
  ! ==================================================================
  SUBROUTINE gclsd(sgcx,sgcc,rhoe,v,vtemp,vtmp,grad,tstress)
    ! ==--------------------------------------------------------------==
    ! ==   CALCULATION OF THE GRADIENT CORRECTION TO THE EXCHANGE     ==
    ! ==   AND CORRELATION FUNCTIONAL FOR LSD CASE                    ==
    ! ==                                                              ==
    ! ==   ON INPUT : RHOE : DENSITY IN REAL SPACE                    ==
    ! ==              V    : UNDEFINED                                ==
    ! ==              GRAD : (NABLA*RHOE)^2,RHOx,RHOy,RHOz ALPHA      ==
    ! ==                     (NABLA*RHOE)^2,RHOx,RHOy,RHOz BETA       ==
    ! ==              VTEMP: POTENTIAL IN G SPACE                     ==
    ! ==              VTMP : DENSITY IN G SPACE                       ==
    ! ==   ON OUTPUT: RHOE : DENSITY IN REAL SPACE                    ==
    ! ==              V    : TOTAL POTENTIAL IN REAL SPACE            ==
    ! ==              GRAD : (NABLA*RHOE)^2,RHOx,RHOy,RHOz ALPHA      ==
    ! ==                     (NABLA*RHOE)^2,RHOx,RHOy,RHOz BETA       ==
    ! ==              VTEMP: UNDEFINED                                ==
    ! ==              VTMP : UNDEFINED                                ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: sgcx, sgcc, rhoe(:,:)
    COMPLEX(real_8)                          :: v(:,:), vtemp(:,:)
    REAL(real_8)                             :: vtmp(:,:), grad(:,:)
    LOGICAL                                  :: tstress, is_present = .false.

    CHARACTER(*), PARAMETER                  :: procedureN = 'gclsd'

    COMPLEX(real_8)                          :: fa, fb, vg1, vg2, vin, vnz
    INTEGER                                  :: ia, ib, ierr, ig, ir, isub, kk
    REAL(real_8)                             :: eg, flops, gcs, smfac, vfac, &
                                                vfac2, vfac3
    REAL(real_8), DIMENSION(:, :), &
      POINTER __CONTIGUOUS                   :: vtau_2d_p

!
! Temporary fix if tau is not allocated (to be nicified in the future)
!

    IF (.NOT.cntl%tgcx .AND. .NOT.cntl%tgcc) RETURN
    CALL tiset(procedureN,isub)
    ! ==--------------------------------------------------------------==
    ! ==  V1 = dF/dn stored in V(*)                                   ==
    ! ==  V2 = (1/Nabla.n) dF/d(Nabla.n) stored in VTMP(*)            ==
    ! ==--------------------------------------------------------------==
    CALL zeroing(v)!,2*maxfftn)
    CALL zeroing(vtmp)!,4*nhg)
    sgcx = 0.0_real_8
    sgcc = 0.0_real_8
    flops=0.0_real_8
    ! ==--------------------------------------------------------------==
    ! ==  POTENTIALS                                                  ==
    ! ==--------------------------------------------------------------==
    IF( cntl%use_xc_driver ) THEN
       !
       ! Temporary fix for metafunctionals (if tau is not needed, allocate
       ! one array for tau and vtau, rather than two)
       !
       IF (cntl%ttau) THEN
          vtau_2d_p(1:SIZE(vtau,1),1:2) => vtau(:,1:2)
       ELSE 
          IF (ALLOCATED(tau)) THEN
             is_present = .true.
          ELSE
          ALLOCATE(tau(fpar%nnr1,2), stat=ierr)
          IF (ierr /= 0) CALL stopgm(procedureN,'allocation problem',& 
               __LINE__,__FILE__)
          ENDIF
          CALL zeroing(tau)
          vtau_2d_p(1:fpar%nnr1,1:2) => tau(:,1:2)
       END IF
       !
       ! End of temporary fix
       !
       CALL cp_xc_compute( fpar%nnr1, rhoe, grad, tau, &
            sgcx, sgcc, v, vtmp, vtau_2d_p )
       !
       ! Temporary fix for metafunctionals
       !
       IF (ASSOCIATED(vtau_2d_p)) NULLIFY(vtau_2d_p)
       IF (.NOT. cntl%ttau) THEN
          IF (.NOT. is_present) THEN
          DEALLOCATE(tau, stat=ierr)
          IF (ierr /= 0) CALL stopgm(procedureN,'deallocation problem',& 
               __LINE__,__FILE__)
          ENDIF
       END IF
       !
       ! End of temporary fix
       !
    ELSE
       IF (toldcode) THEN
          CALL gclsdold(sgcx,sgcc,rhoe,v,vtmp,grad,flops)
       ELSE
          IF (func1%mgcx == mgcx_is_becke88) THEN
             IF (func1%mgcc == mgcc_is_skipped) THEN
                CALL gcsxonly(func2%bbeta,sgcx,sgcc,rhoe,v(:,1),v(:,2),&
                     vtmp,grad,flops)
             ELSEIF (func1%mgcc == mgcc_is_perdew86) THEN
                CALL gcsxp86(func2%bbeta,sgcx,sgcc,rhoe,v(:,1),v(:,2),&
                     vtmp,grad,flops)
             ELSEIF (func1%mgcc == mgcc_is_lyp) THEN
                CALL gcsxlyp(func2%bbeta,sgcx,sgcc,rhoe,v(:,1),v(:,2),&
                     vtmp,grad,flops)
             ELSE
                CALL stopgm(procedureN,'NOT IMPLEMENTED',& 
                     __LINE__,__FILE__)
             ENDIF
          ELSEIF (func1%mgcx == mgcx_is_pbex .AND. func1%mgcc == mgcc_is_pbec) THEN
             CALL gcspbe(sgcx,sgcc,rhoe,v(:,1),v(:,2),vtmp,grad,flops)
          ELSEIF (func1%mgcx == mgcx_is_revpbex .AND. func1%mgcc == mgcc_is_pbec) THEN
             CALL gcsrevpbe(sgcx,sgcc,rhoe,v(:,1),v(:,2),vtmp,grad,flops)
          ELSE
             ! this is also somewhat broken for some combinations of functionals.
             CALL gclsdold(sgcx,sgcc,rhoe,v,vtmp,grad,flops)
          ENDIF
       ENDIF
    ENDIF

    IF (tstress) THEN
       ! Contribution of dF/d(delta.n) to the stress tensor
       !$omp parallel do private(KK,IR,IA,IB)
       DO kk=1,6
          ia=alpha(kk)+1
          ib=beta(kk)+1
          degc(kk)=0.0_real_8
          DO ir=1,fpar%nnr1
             degc(kk)=degc(kk)+(vtmp(ir,1)+vtmp(ir,3))*&
                  grad(ir,ia)*grad(ir,ib)
          ENDDO
       ENDDO
       DO kk=1,6
          ia=alpha(kk)+5
          ib=beta(kk)+5
          DO ir=1,fpar%nnr1
             degc(kk)=degc(kk)+(vtmp(ir,2)+vtmp(ir,3))*&
                  grad(ir,ia)*grad(ir,ib)
          ENDDO
          degc(kk)=degc(kk)*parm%omega/(spar%nr1s*spar%nr2s*spar%nr3s)
       ENDDO
    ENDIF
    ! HACK_TEST
    ! ==--------------------------------------------------------------==
    ! ==  ALPHA SPIN                                                  ==
    ! ==--------------------------------------------------------------==
    ! FFT of V1 and V2*RHOx to G-Space
    !$omp parallel do private(IR)
    DO ir=1,fpar%nnr1
       v(ir,1) = v(ir,1) + CMPLX(0.0_real_8,vtmp(ir,1)*grad(ir,2),kind=real_8)&
            + CMPLX(0.0_real_8,vtmp(ir,3)*grad(ir,6),kind=real_8)
    ENDDO
    CALL  fwfftn(v(:,1),.FALSE.,parai%allgrp)
    gcs=cntr%smf*hg(ncpw%nhg)
    IF (cntl%tsmooth) THEN
       gcs=cntr%smf*hg(ncpw%nhg)
       DO ig=1,ncpw%nhg
          eg = (hg(ig)-gcs)/cntr%sdelta
          smfac=1.0_real_8/(1.0_real_8+EXP(eg))
          vfac=smfac*parm%tpiba*gk(1,ig)
          ! !        FA  = V(NZH(IG),1) + V(INDZ(IG),1)
          ! !        FB  = V(NZH(IG),1) - V(INDZ(IG),1)
          vnz = v(nzh(ig),1)
          vin = v(indz(ig),1)
          fa  = vnz + vin
          fb  = vnz - vin
          vg1 = 0.5_real_8*CMPLX(REAL(fa),AIMAG(fb),kind=real_8)
          vg2 = 0.5_real_8*CMPLX(AIMAG(fa),-REAL(fb),kind=real_8)
          vtemp(ig,1) = vtemp(ig,1) + smfac*vg1 - vfac*uimag*vg2
       ENDDO
    ELSE
       !$omp parallel do private(IG,VFAC,vnz,vin,FA,FB,VG1,VG2)
       DO ig=1,ncpw%nhg
          vfac=parm%tpiba*gk(1,ig)
          ! !        FA  = V(NZH(IG),1) + V(INDZ(IG),1)
          ! !        FB  = V(NZH(IG),1) - V(INDZ(IG),1)
          vnz = v(nzh(ig),1)
          vin = v(indz(ig),1)
          fa  = vnz + vin
          fb  = vnz - vin
          vg1 = 0.5_real_8*CMPLX(REAL(fa),AIMAG(fb),kind=real_8)
          vg2 = 0.5_real_8*CMPLX(AIMAG(fa),-REAL(fb),kind=real_8)
          vtemp(ig,1) = vtemp(ig,1) + vg1 - vfac*uimag*vg2
       ENDDO
    ENDIF
    ! FFT of V2*RHOy and V2*RHOz to G-Space
    CALL zeroing(v(:,1))!,nnr1)
    !$omp parallel do private(IR)
    DO ir=1,fpar%nnr1
       v(ir,1) = CMPLX(vtmp(ir,1)*grad(ir,3),vtmp(ir,1)*grad(ir,4),kind=real_8)&
            +CMPLX(vtmp(ir,3)*grad(ir,7),vtmp(ir,3)*grad(ir,8),kind=real_8)
    ENDDO
    CALL  fwfftn(v(:,1),.FALSE.,parai%allgrp)
    gcs=cntr%smf*hg(ncpw%nhg)
    IF (cntl%tsmooth) THEN
       gcs=cntr%smf*hg(ncpw%nhg)
       DO ig=1,ncpw%nhg
          eg = (hg(ig)-gcs)/cntr%sdelta
          smfac=1.0_real_8/(1.0_real_8+EXP(eg))
          vfac2=smfac*parm%tpiba*gk(2,ig)
          vfac3=smfac*parm%tpiba*gk(3,ig)
          ! !        FA  = V(NZH(IG),1) + V(INDZ(IG),1)
          ! !        FB  = V(NZH(IG),1) - V(INDZ(IG),1)
          vnz = v(nzh(ig),1)
          vin = v(indz(ig),1)
          fa  = vnz + vin
          fb  = vnz - vin
          vg1 = 0.5_real_8*CMPLX(REAL(fa),AIMAG(fb),kind=real_8)
          vg2 = 0.5_real_8*CMPLX(AIMAG(fa),-REAL(fb),kind=real_8)
          vtemp(ig,1) = vtemp(ig,1) - vfac2*uimag*vg1 - vfac3*uimag*vg2
       ENDDO
    ELSE
       !$omp parallel do private(IG,VFAC2,VFAC3,vnz,vin,FA,FB,VG1,VG2)
       DO ig=1,ncpw%nhg
          vfac2=parm%tpiba*gk(2,ig)
          vfac3=parm%tpiba*gk(3,ig)
          ! !        FA  = V(NZH(IG),1) + V(INDZ(IG),1)
          ! !        FB  = V(NZH(IG),1) - V(INDZ(IG),1)
          vnz = v(nzh(ig),1)
          vin = v(indz(ig),1)
          fa  = vnz + vin
          fb  = vnz - vin
          vg1 = 0.5_real_8*CMPLX(REAL(fa),AIMAG(fb),kind=real_8)
          vg2 = 0.5_real_8*CMPLX(AIMAG(fa),-REAL(fb),kind=real_8)
          vtemp(ig,1) = vtemp(ig,1) - vfac2*uimag*vg1 - vfac3*uimag*vg2
       ENDDO
    ENDIF
    CALL zeroing(v(:,1))!,nnr1)
    !CDIR NODEP
    !ocl novrec
#if defined(__SR11000) || defined(__PRIMERGY) || defined(__PRIMEHPC)
    !$omp parallel do private(IG)
#endif
    DO ig=1,ncpw%nhg
       v(nzh(ig),1) = vtemp(ig,1)
       v(indz(ig),1) = CONJG(vtemp(ig,1))
    ENDDO
    CALL  invfftn(v(:,1),.FALSE.,parai%allgrp)
    ! ==--------------------------------------------------------------==
    ! ==  BETA SPIN                                                   ==
    ! ==--------------------------------------------------------------==
    ! FFT of V1 and V2*RHOx to G-Space
    !$omp parallel do private(IR)
    DO ir=1,fpar%nnr1
       v(ir,2) = v(ir,2) + CMPLX(0.0_real_8,vtmp(ir,2)*grad(ir,6),kind=real_8)&
            + CMPLX(0.0_real_8,vtmp(ir,3)*grad(ir,2),kind=real_8)
    ENDDO
    CALL  fwfftn(v(:,2),.FALSE.,parai%allgrp)
    gcs=cntr%smf*hg(ncpw%nhg)
    IF (cntl%tsmooth) THEN
       gcs=cntr%smf*hg(ncpw%nhg)
       DO ig=1,ncpw%nhg
          eg = (hg(ig)-gcs)/cntr%sdelta
          smfac=1.0_real_8/(1.0_real_8+EXP(eg))
          vfac=smfac*parm%tpiba*gk(1,ig)
          ! !        FA  = V(NZH(IG),2) + V(INDZ(IG),2)
          ! !        FB  = V(NZH(IG),2) - V(INDZ(IG),2)
          vnz = v(nzh(ig),2)
          vin = v(indz(ig),2)
          fa  = vnz + vin
          fb  = vnz - vin
          vg1 = 0.5_real_8*CMPLX(REAL(fa),AIMAG(fb),kind=real_8)
          vg2 = 0.5_real_8*CMPLX(AIMAG(fa),-REAL(fb),kind=real_8)
          vtemp(ig,2) = vtemp(ig,2) + smfac*vg1 - vfac*uimag*vg2
       ENDDO
    ELSE
       !$omp parallel do private(IG,VFAC,vnz,vin,FA,FB,VG1,VG2)
       DO ig=1,ncpw%nhg
          vfac=parm%tpiba*gk(1,ig)
          ! !        FA  = V(NZH(IG),2) + V(INDZ(IG),2)
          ! !        FB  = V(NZH(IG),2) - V(INDZ(IG),2)
          vnz = v(nzh(ig),2)
          vin = v(indz(ig),2)
          fa  = vnz + vin
          fb  = vnz - vin
          vg1 = 0.5_real_8*CMPLX(REAL(fa),AIMAG(fb),kind=real_8)
          vg2 = 0.5_real_8*CMPLX(AIMAG(fa),-REAL(fb),kind=real_8)
          vtemp(ig,2) = vtemp(ig,2) + vg1 - vfac*uimag*vg2
       ENDDO
    ENDIF
    ! FFT of V2*RHOy and V2*RHOz to G-Space
    CALL zeroing(v(:,2))!,nnr1)
    !$omp parallel do private(IR)
    DO ir=1,fpar%nnr1
       v(ir,2) = CMPLX(vtmp(ir,2)*grad(ir,7),vtmp(ir,2)*grad(ir,8),kind=real_8)&
            +CMPLX(vtmp(ir,3)*grad(ir,3),vtmp(ir,3)*grad(ir,4),kind=real_8)
    ENDDO
    CALL  fwfftn(v(:,2),.FALSE.,parai%allgrp)
    gcs=cntr%smf*hg(ncpw%nhg)
    IF (cntl%tsmooth) THEN
       gcs=cntr%smf*hg(ncpw%nhg)
       DO ig=1,ncpw%nhg
          eg = (hg(ig)-gcs)/cntr%sdelta
          smfac=1.0_real_8/(1.0_real_8+EXP(eg))
          vfac2=smfac*parm%tpiba*gk(2,ig)
          vfac3=smfac*parm%tpiba*gk(3,ig)
          ! !        FA  = V(NZH(IG),2) + V(INDZ(IG),2)
          ! !        FB  = V(NZH(IG),2) - V(INDZ(IG),2)
          vnz = v(nzh(ig),2)
          vin = v(indz(ig),2)
          fa  = vnz + vin
          fb  = vnz - vin
          vg1 = 0.5_real_8*CMPLX(REAL(fa),AIMAG(fb),kind=real_8)
          vg2 = 0.5_real_8*CMPLX(AIMAG(fa),-REAL(fb),kind=real_8)
          vtemp(ig,2) = vtemp(ig,2) - vfac2*uimag*vg1 - vfac3*uimag*vg2
       ENDDO
    ELSE
       !$omp parallel do private(IG,VFAC2,VFAC3,vnz,vin,FA,FB,VG1,VG2)
       DO ig=1,ncpw%nhg
          vfac2=parm%tpiba*gk(2,ig)
          vfac3=parm%tpiba*gk(3,ig)
          ! !        FA  = V(NZH(IG),2) + V(INDZ(IG),2)
          ! !        FB  = V(NZH(IG),2) - V(INDZ(IG),2)
          vnz = v(nzh(ig),2)
          vin = v(indz(ig),2)
          fa  = vnz + vin
          fb  = vnz - vin
          vg1 = 0.5_real_8*CMPLX(REAL(fa),AIMAG(fb),kind=real_8)
          vg2 = 0.5_real_8*CMPLX(AIMAG(fa),-REAL(fb),kind=real_8)
          vtemp(ig,2) = vtemp(ig,2) - vfac2*uimag*vg1 - vfac3*uimag*vg2
       ENDDO
    ENDIF
    CALL zeroing(v(:,2))!,nnr1)
    !CDIR NODEP
    !ocl novrec
#if defined(__SR11000) || defined(__PRIMERGY) || defined(__PRIMEHPC)
    !$omp parallel do private(IG)
#endif
    DO ig=1,ncpw%nhg
       v(nzh(ig),2) = vtemp(ig,2)
       v(indz(ig),2) = CONJG(vtemp(ig,2))
    ENDDO
    CALL  invfftn(v(:,2),.FALSE.,parai%allgrp)
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE gclsd
  ! ==================================================================
  SUBROUTINE gcold(sgcx,sgcc,rhoe,v,vtmp,grad,flops)
    ! ==--------------------------------------------------------------==

    REAL(real_8)                             :: sgcx, sgcc, rhoe(:)
    COMPLEX(real_8)                          :: v(:)
    REAL(real_8)                             :: vtmp(:), grad(:,:), flops

    INTEGER                                  :: ir
    REAL(real_8)                             :: dsmoo, grho, rho, sc, smoo, &
                                                sx, sxsr, tar, texp, v1, v1c, &
                                                v1x, v1xsr, v2, v2c, v2x, &
                                                v2xsr, vtt

    IF (.NOT.cntl%ttau) THEN
       !$omp parallel do private(IR,RHO,SMOO,DSMOO,TEXP,GRHO,SX,SC) &
       !$omp             private(V1X,V2X,V1C,V2C,V1,V2,TAR) &
       !$omp             private(SXSR,V1XSR,V2XSR) &
       !$omp             shared(GRAD,V,VTMP,fpar,RHOE) &
       !$omp             reduction(+:SGCX,SGCC)
#ifdef __SR8000
       !poption parallel
       !poption tlocal(SX,SC,V1X,V2X,V1C,V2C)
#endif
#if defined (__SR11000)
       !poption psum(SGCX,SGCC)
#endif
       DO ir=1,fpar%nnr1
          rho   = MAX(rhoe(ir),0.0_real_8)
          IF (rho.GT.0.1_real_8*cntr%gceps) THEN
             IF (rho.GT.4.0_real_8*cntr%gceps) THEN
                smoo  = 1.0_real_8
                dsmoo = 0.0_real_8
             ELSE
                texp=EXP(3.0_real_8*(1.0_real_8-rho/cntr%gceps))
                smoo=1.0_real_8/(1.0_real_8+texp)
                dsmoo=3.0_real_8/cntr%gceps * texp*smoo*smoo
             ENDIF
             grho    = grad(ir,1)
             CALL gcxc(rho,grho,sx,sc,v1x,v2x,v1c,v2c)
             IF (cntl%thybrid) THEN
                ! wPBE stabbing, adding short-range GC
                IF (func1%mgcsrx /= mgcsrx_is_skipped) THEN
                   CALL pbexsr(rho,grho,sxsr,v1xsr,v2xsr,func2%srxa)
                   SX=SX*func3%pxgc-func3%phfx*SXSR;
                   V1X=V1X*func3%pxgc-func3%phfx*V1XSR;
                   V2X=V2X*func3%pxgc-func3%phfx*V2XSR;
                ELSE
                   sx=sx*func3%pxgc
                   v1x=v1x*func3%pxgc
                   v2x=v2x*func3%pxgc
                ENDIF
                sc=sc*func3%pcgc
                v1c=v1c*func3%pcgc
                v2c=v2c*func3%pcgc
             ENDIF
             sgcx  = sgcx + smoo*sx
             sgcc  = sgcc + smoo*sc
             v1    = dsmoo*(sx+sc) + smoo*(v1x+v1c)
             v2    = smoo*(v2x+v2c)
          ELSE
             v1    = 0.0_real_8
             v2    = 0.0_real_8
          ENDIF
          v(ir) = v1
          vtmp(ir) = v2
       ENDDO
    ELSE
       DO ir=1,fpar%nnr1
          rho   = MAX(rhoe(ir),0.0_real_8)
          IF (rho.GT.0.1_real_8*cntr%gceps) THEN
             IF (rho.GT.4.0_real_8*cntr%gceps) THEN
                smoo  = 1.0_real_8
                dsmoo = 0.0_real_8
             ELSE
                texp=EXP(3.0_real_8*(1.0_real_8-rho/cntr%gceps))
                smoo=1.0_real_8/(1.0_real_8+texp)
                dsmoo=3.0_real_8/cntr%gceps * texp*smoo*smoo
             ENDIF
             grho    = grad(ir,1)
             tar=tau(ir,1)
             CALL taufun(rho,grho,tar,sx,sc,v1x,v2x,v1c,v2c,vtt)
             sgcx  = sgcx + smoo*sx
             sgcc  = sgcc + smoo*sc
             v1    = dsmoo*(sx+sc) + smoo*(v1x+v1c)
             v2    = smoo*(v2x+v2c)
          ELSE
             v1    = 0.0_real_8
             v2    = 0.0_real_8
             vtt   = 0.0_real_8
          ENDIF
          v(ir) = v1
          vtmp(ir) = v2
          vtau(ir,1) = vtt
       ENDDO
    ENDIF
    flops=0._real_8
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gcold
  ! ==================================================================
  SUBROUTINE gclsdold(sgcx,sgcc,rhoe,v,vtmp,grad,flops)
    ! ==--------------------------------------------------------------==

    REAL(real_8)                             :: sgcx, sgcc, rhoe(:,:)
    COMPLEX(real_8)                          :: v(:,:)
    REAL(real_8)                             :: vtmp(:,:), grad(:,:), flops

    CHARACTER(*), PARAMETER                  :: procedureN = 'gclsdold'

    INTEGER                                  :: ir, isub
    REAL(real_8) :: dsmoo, grhoa, grhoab, grhob, rho, rhoa, rhob, sc, smoo, &
      sx, sxsr, tara, tarb, texp, v1a, v1b, v1ca, v1cb, v1xa, v1xasr, v1xb, &
      v1xbsr, v2a, v2ab, v2b, v2ca, v2cab, v2cb, v2xa, v2xab, v2xabsr, &
      v2xasr, v2xb, v2xbsr, vtta, vttb

    CALL tiset(procedureN,isub)
    flops=0._real_8
    IF (.NOT.cntl%ttau) THEN
       !$omp parallel do private(IR,RHO,RHOA,RHOB,SMOO,DSMOO,TEXP,SX,SC) &
       !$omp       private(GRHOA,GRHOB,GRHOAB,V1XA,V1XB,V2XA,V2XB,V1CA,V2CA) &
       !$omp       private(V1CB,V2CB,V2XAB,V2CAB,V1A,V2A,V1B,V2B,V2AB) &
       !$omp       private(TARA,TARB,VTTA,VTTB) &
       !$omp       private(SXSR,V1XASR,V2XASR,V1XBSR,V2XBSR,V2XABSR) &
       !$omp       shared(GRAD,V,VTMP,fpar,RHOE) &
       !$omp       reduction(+:SGCX,SGCC)
#ifdef __SR8000
       !poption parallel
       !poption tlocal(SX,SC)
       !poption tlocal(V1XA,V2XA,V1XB,V2XB,V1CA,V2CA,V1CB,V2CB)
       !poption tlocal(V2XAB,V2CAB)
#endif
#if defined (__SR11000)
       !poption psum(SGCX,SGCC)
#endif
       DO ir=1,fpar%nnr1
          rhoa   = MAX(rhoe(ir,1),0.0_real_8)
          rhob   = MAX(rhoe(ir,2),0.0_real_8)
          rho    = rhoa+rhob
          IF (rho.GT.0.1_real_8*cntr%gceps) THEN
             IF (rho.GT.4.0_real_8*cntr%gceps) THEN
                smoo  = 1.0_real_8
                dsmoo = 0.0_real_8
             ELSE
                texp=EXP(3.0_real_8*(1.0_real_8-rho/cntr%gceps))
                smoo=1.0_real_8/(1.0_real_8+texp)
                dsmoo=3.0_real_8/cntr%gceps * texp*smoo*smoo
             ENDIF
             IF (rhoa.LT.0.1_real_8*cntr%gceps) rhoa=0.0_real_8
             IF (rhob.LT.0.1_real_8*cntr%gceps) rhob=0.0_real_8
             grhoa   = grad(ir,1)
             grhob   = grad(ir,5)
             grhoab  = grad(ir,2)*grad(ir,6)+grad(ir,3)*grad(ir,7)+&
                  grad(ir,4)*grad(ir,8)
             CALL gc_lsd(rhoa,rhob,grhoa,grhob,grhoab,sx,sc,v1xa,v2xa,&
                  v1xb,v2xb,v1ca,v2ca,v1cb,v2cb,v2xab,v2cab)
             IF (cntl%thybrid) THEN
                ! wPBE stabbing, adding short-range GC
                IF (func1%mgcsrx /= mgcsrx_is_skipped) THEN
                   CALL pbexsr_lsd(rhoa,rhob,grhoa,grhob,sxsr,&
                        v1xasr,v2xasr,v1xbsr,v2xbsr,v2xabsr,func2%srxa)
                   sx=sx*func3%pxgc - func3%phfx*sxsr
                   v1xa=v1xa*func3%pxgc - func3%phfx*v1xasr
                   v2xa=v2xa*func3%pxgc - func3%phfx*v2xasr
                   v1xb=v1xb*func3%pxgc - func3%phfx*v1xbsr
                   v2xb=v2xb*func3%pxgc - func3%phfx*v2xbsr
                   v2xab=v2xab*func3%pxgc - func3%phfx*v2xabsr
                ELSE
                   sx=sx*func3%pxgc
                   v1xa=v1xa*func3%pxgc
                   v2xa=v2xa*func3%pxgc
                   v1xb=v1xb*func3%pxgc
                   v2xb=v2xb*func3%pxgc
                   v2xab=v2xab*func3%pxgc
                ENDIF
                sc=sc*func3%pcgc
                v1ca=v1ca*func3%pcgc
                v2ca=v2ca*func3%pcgc
                v1cb=v1cb*func3%pcgc
                v2cb=v2cb*func3%pcgc
                v2cab=v2cab*func3%pcgc
             ENDIF
             sgcx  = sgcx + smoo*sx
             sgcc  = sgcc + smoo*sc
             v1a   = dsmoo*(sx+sc) + smoo*(v1xa+v1ca)
             v2a   = smoo*(v2xa+v2ca)
             v1b   = dsmoo*(sx+sc) + smoo*(v1xb+v1cb)
             v2b   = smoo*(v2xb+v2cb)
             v2ab  = smoo*(v2xab+v2cab)
          ELSE
             v1a   = 0.0_real_8
             v1b   = 0.0_real_8
             v2a   = 0.0_real_8
             v2b   = 0.0_real_8
             v2ab  = 0.0_real_8
             vtta  = 0.0_real_8
             vttb  = 0.0_real_8
          ENDIF
          v(ir,1) = v1a
          v(ir,2) = v1b
          vtmp(ir,1) = v2a
          vtmp(ir,2) = v2b
          vtmp(ir,3) = v2ab
       ENDDO
    ELSE
       DO ir=1,fpar%nnr1
          rhoa   = MAX(rhoe(ir,1),0.0_real_8)
          rhob   = MAX(rhoe(ir,2),0.0_real_8)
          rho    = rhoa+rhob
          IF (rho.GT.0.1_real_8*cntr%gceps) THEN
             IF (rho.GT.4.0_real_8*cntr%gceps) THEN
                smoo  = 1.0_real_8
                dsmoo = 0.0_real_8
             ELSE
                texp=EXP(3.0_real_8*(1.0_real_8-rho/cntr%gceps))
                smoo=1.0_real_8/(1.0_real_8+texp)
                dsmoo=3.0_real_8/cntr%gceps * texp*smoo*smoo
             ENDIF
             IF (rhoa.LT.0.1_real_8*cntr%gceps) rhoa=0.0_real_8
             IF (rhob.LT.0.1_real_8*cntr%gceps) rhob=0.0_real_8
             grhoa   = grad(ir,1)
             grhob   = grad(ir,5)
             grhoab  = grad(ir,2)*grad(ir,6)+grad(ir,3)*grad(ir,7)+&
                  grad(ir,4)*grad(ir,8)
             tara = tau(ir,1)
             tarb = tau(ir,2)
             CALL taufuns(rhoa,rhob,grhoa,grhob,grhoab,tara,tarb,&
                  sx,sc,v1xa,v2xa,&
                  v1xb,v2xb,v1ca,v2ca,v1cb,v2cb,v2xab,v2cab,&
                  vtta,vttb)
             sgcx  = sgcx + smoo*sx
             sgcc  = sgcc + smoo*sc
             v1a   = dsmoo*(sx+sc) + smoo*(v1xa+v1ca)
             v2a   = smoo*(v2xa+v2ca)
             v1b   = dsmoo*(sx+sc) + smoo*(v1xb+v1cb)
             v2b   = smoo*(v2xb+v2cb)
             v2ab  = smoo*(v2xab+v2cab)
          ELSE
             v1a   = 0.0_real_8
             v1b   = 0.0_real_8
             v2a   = 0.0_real_8
             v2b   = 0.0_real_8
             v2ab  = 0.0_real_8
             vtta  = 0.0_real_8
             vttb  = 0.0_real_8
          ENDIF
          v(ir,1) = v1a
          v(ir,2) = v1b
          vtmp(ir,1) = v2a
          vtmp(ir,2) = v2b
          vtmp(ir,3) = v2ab
          vtau(ir,1) = vtta
          vtau(ir,2) = vttb
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE gclsdold
  ! ==================================================================

END MODULE gcener_utils
