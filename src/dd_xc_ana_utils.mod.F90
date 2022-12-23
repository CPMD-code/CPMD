MODULE dd_xc_ana_utils
  USE cnst,                            ONLY: uimag
  USE cppt,                            ONLY: gk,&
                                             hg,&
                                             indz,&
                                             nzh
  USE cp_xc_utils,                     ONLY: cp_xc
  USE cp_dxc_driver,                   ONLY: cp_dxc_compute
  USE dd_functionals_utils,            ONLY: &
       b88_x, lyp88_c, lyp88_c_loc, p86_c, p86_c_loc, pbe96_c, pbe96_c_loc, &
       pbe96_x, slater_x
  USE error_handling,                  ONLY: stopgm
  USE fft_maxfft,                      ONLY: maxfft
  USE fftmain_utils,                   ONLY: fwfftn,&
                                             invfftn
  USE fftnew_utils,                    ONLY: setfftn
  USE func,                            ONLY: func1,&
                                             mgcc_is_lyp,&
                                             mgcc_is_pbec,&
                                             mgcc_is_perdew86,&
                                             mgcc_is_skipped,&
                                             mgcx_is_becke88,&
                                             mgcx_is_pbex
  USE graden_utils,                    ONLY: graden
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_max
  USE parac,                           ONLY: parai
  USE spin,                            ONLY: clsd,&
                                             lspin2
  USE switch_functionals_utils,        ONLY: switch_functionals
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

  PUBLIC :: dd_xc_ana
  !public :: graden_comp
  PUBLIC :: give_scr_dd_xc_ana
  !public :: fftmult

CONTAINS

  ! ==================================================================
  SUBROUTINE dd_xc_ana(rhoe,drhoe,rhox,v,dxc,vtmp,vtemp,grad,switch)
    ! ==--------------------------------------------------------------==
    ! Calculate dVxc/dn*n1 using analytical second derivatives
    ! ==--------------------------------------------------------------==
    ! == INPUT: rhoe = density (a+b,b)                                ==
    ! ==        drhoe = LR density (a+b,b)                            ==
    ! ==        v = undef. (size: psi)                                ==
    ! ==        dxc = undef.                                          ==
    ! ==        vtmp = undef.                                         ==
    ! ==        vtemp = undef.                                        ==
    ! ==        grad = undef.                                         ==
    ! == OUTPUT: rhoe = density (a+b,b)                               ==
    ! ==        drhoe = LR density (a+b,b)                            ==
    ! ==        v = undef. (size: psi)                                ==
    ! ==        dxc = dVxc/dn*n1                                      ==
    ! ==        vtmp = undef.                                         ==
    ! ==        vtemp = undef.                                        ==
    ! ==        grad = (|grad|**2,x,y,z)                              ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rhoe(:,:), drhoe(:,:), &
                                                rhox(:,:)
    COMPLEX(real_8)                          :: v(:,:)
    REAL(real_8)                             :: dxc(:,:)
    COMPLEX(real_8)                          :: vtmp(:), vtemp(:)
    REAL(real_8)                             :: grad(:,:)
    LOGICAL                                  :: switch

    CHARACTER(*), PARAMETER                  :: procedureN = 'dd_xc_ana'
    REAL(real_8), PARAMETER :: &
      a = -0.930525736349100025002010218071667_real_8, hf = 0.5_real_8, &
      qt = 0.25_real_8, tw = 2.0_real_8

    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: ierr, ig, ir, ispin, isub, &
                                                ldd_xc_ana
    LOGICAL                                  :: debug
    REAL(real_8) :: d2fdg2, d2fdn2, d2fdndg, dfdg, dfdga, dfdn, ga, gaga, &
      gagb, ganb, gask, granb, lda, loca, locb, naga, nagb, nagrb, nana, &
      nanb, nask, p86a, p86b, sgrad, sgradb, sk, skaa, skab, skalar, skba, &
      skbb, skgb, sknb, sksk
    REAL(real_8), ALLOCATABLE                :: dgrad(:,:), tmp1(:), tmp2(:), &
                                                tmp3(:), tmp4(:), tmp5(:)
    REAL(real_8), DIMENSION(:,:,:), &
                  ALLOCATABLE                :: dxc_tmp

!(nnr1,clsd%nlsd)
!maxfft,clsd%nlsd
!(nnr1,*)
!nnr1,4*clsd%nlsd
! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    ! ==--------------------------------------------------------------==
    CALL setfftn(0)
    CALL give_scr_dd_xc_ana(ldd_xc_ana,tag)
    debug=.FALSE.
    !
    IF (associated(cp_xc)) THEN
       !
       ! If xc_driver is not used, cp_xc is not necessarily associated
       !
       IF (cp_xc%use_libxc) CALL stopgm(procedureN,'Analytical derivatives are not available with libxc', &
                                     __LINE__, __FILE__)
    ENDIF
    !
    ! Allocation of scratch space for gradients (if necessary)
    !
    IF (cntl%tgc .and. .not. cntl%use_xc_driver) THEN
       ALLOCATE(tmp1(fpar%nnr1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       CALL dcopy(fpar%nnr1, rhox(1,1), 1, tmp1, 1)
       ALLOCATE(tmp2(fpar%nnr1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       IF (cntl%tlsd) THEN
          CALL dcopy(fpar%nnr1, rhox(1,2), 1, tmp2, 1)
       ENDIF
       ALLOCATE(tmp3(fpar%nnr1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(tmp4(fpar%nnr1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(tmp5(fpar%nnr1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(dgrad(fpar%nnr1, 3*clsd%nlsd),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)

       CALL zeroing(tmp1)!,nnr1)
       CALL zeroing(tmp2)!,nnr1)
       CALL zeroing(tmp3)!,nnr1)
       CALL zeroing(tmp4)!,nnr1)
       CALL zeroing(tmp5)!,nnr1)
       CALL zeroing(dgrad)!,nnr1*3*clsd%nlsd)
    ELSE IF (cntl%tgc .AND. cntl%use_xc_driver) THEN
       ALLOCATE(dxc_tmp(fpar%nnr1,5,clsd%nlsd),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(dgrad(fpar%nnr1, 3*clsd%nlsd),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       CALL zeroing(dxc_tmp)
       CALL zeroing(dgrad)
    ENDIF
    ! Calculation of gradients (if necessary)
    IF (.NOT.cntl%tlsd.AND.cntl%tgc) THEN
       CALL graden(rhoe,v,grad,vtmp)
       CALL graden_comp(drhoe,v,dgrad,vtmp)
    ELSEIF (cntl%tlsd.AND..NOT.cntl%tgc) THEN
       DO ir=1,fpar%nnr1
          rhoe(ir,1)=rhoe(ir,1)-rhoe(ir,2)
          drhoe(ir,1)=drhoe(ir,1)-drhoe(ir,2)
       ENDDO
    ELSEIF (cntl%tlsd.AND.cntl%tgc) THEN
       DO ir=1,fpar%nnr1
          rhoe(ir,1)=rhoe(ir,1)-rhoe(ir,2)
          drhoe(ir,1)=drhoe(ir,1)-drhoe(ir,2)
       ENDDO
       CALL graden(rhoe,v,grad,vtmp)
       CALL graden(rhoe(:,2),v,grad(:,5),vtmp)
       CALL graden_comp(drhoe,v,dgrad,vtmp)
       CALL graden_comp(drhoe(:,2),v,dgrad(:,4),vtmp)
    ENDIF

    IF (switch) CALL switch_functionals

    IF (cntl%tlsd) THEN
       CALL zeroing(v)!,2*maxfft)
    ELSE
       CALL zeroing(v)!,maxfft)
    ENDIF
    IF (cntl%tgc) CALL zeroing(vtemp)!,nhg)

    IF (cntl%use_xc_driver) THEN
      IF (.not. cntl%tlsd) THEN
         ! no spin polarisation --> a-dens = b-dens = 1/2 * dens      
         CALL dscal(fpar%nnr1,hf,rhoe,1)
         CALL dscal(fpar%nnr1,hf,drhoe,1)
         IF (cntl%tgc) THEN
            CALL dscal(fpar%nnr1,qt,grad(1,1),1)
            CALL dscal(3*fpar%nnr1,hf,grad(1,2),1)
            CALL dscal(3*fpar%nnr1,hf,dgrad(1,1),1)
         ENDIF
         CALL cp_dxc_compute( fpar%nnr1, rhoe, drhoe, grad, dgrad, dxc_tmp)
         !
         ! CALCULATION OF DIVERGENCE IN RECIPROCAL SPACE:
         ! real space               G-space
         ! dot(nabla,gradient) --> dot(i*G-Vector,gradient)
         ! tmp2: terms to be multiplied with gradient of a-density
         ! tmp3: terms to be multiplied with gradient of a-linres-density
         ! tmp4: terms to be multiplied with gradient of b-density
         ! tmp5: terms to be multiplied with gradient of b-linres-density
         !
         ! FFTMULT:calculates dot product of i*G-Vector with gradient 
         DO ir=1,fpar%nnr1
            v(ir,1)=CMPLX(dxc_tmp(ir,2,1)*grad(ir,2),dxc_tmp(ir,2,1)*grad(ir,3),kind=real_8)
         ENDDO
         CALL  fwfftn(v(:,1),.FALSE.,parai%allgrp)
         CALL fftmult(v,vtemp,1,2)

         CALL zeroing(v(:,1))!,maxfft)
         DO ir=1,fpar%nnr1
            v(ir,1)=CMPLX(dxc_tmp(ir,2,1)*grad(ir,4),dxc_tmp(ir,3,1)*dgrad(ir,1),kind=real_8)
         ENDDO
         CALL  fwfftn(v(:,1),.FALSE.,parai%allgrp)
         CALL fftmult(v,vtemp,3,1)

         CALL zeroing(v(:,1))!,maxfft)
         DO ir=1,fpar%nnr1
            v(ir,1)=CMPLX(dxc_tmp(ir,3,1)*dgrad(ir,2),dxc_tmp(ir,3,1)*dgrad(ir,3),kind=real_8)
         ENDDO
         CALL  fwfftn(v(:,1),.FALSE.,parai%allgrp)
         CALL fftmult(v,vtemp,2,3)

         CALL zeroing(v(:,1))!,maxfft)
         DO ir=1,fpar%nnr1
            v(ir,1)=CMPLX(dxc_tmp(ir,4,1)*grad(ir,2),dxc_tmp(ir,4,1)*grad(ir,3),kind=real_8)
         ENDDO
         CALL  fwfftn(v(:,1),.FALSE.,parai%allgrp)
         CALL fftmult(v,vtemp,1,2)

         CALL zeroing(v(:,1))!,maxfft)
         DO ir=1,fpar%nnr1
            v(ir,1)=CMPLX(dxc_tmp(ir,4,1)*grad(ir,4),dxc_tmp(ir,5,1)*dgrad(ir,1),kind=real_8)
         ENDDO
         CALL  fwfftn(v(:,1),.FALSE.,parai%allgrp)
         CALL fftmult(v,vtemp,3,1)

         CALL zeroing(v(:,1))!,maxfft)
         DO ir=1,fpar%nnr1
            v(ir,1)=CMPLX(dxc_tmp(ir,5,1)*dgrad(ir,2),dxc_tmp(ir,5,1)*dgrad(ir,3),kind=real_8)
         ENDDO
         CALL  fwfftn(v(:,1),.FALSE.,parai%allgrp)
         CALL fftmult(v,vtemp,2,3)

         CALL zeroing(v(:,1))!,maxfft)
         DO ig=1,ncpw%nhg
            v(nzh(ig),1) = vtemp(ig)
            v(indz(ig),1) = CONJG(vtemp(ig))
         ENDDO

         CALL  invfftn(v(:,1),.FALSE.,parai%allgrp)
         DO ir=1,fpar%nnr1
            v(ir,1)=v(ir,1)+CMPLX(dxc_tmp(ir,1,1),0._real_8,kind=real_8)
         ENDDO
         CALL dscal(fpar%nnr1,tw,rhoe,1)
         CALL dscal(fpar%nnr1,tw,drhoe,1)
         IF (cntl%tgc) THEN
            CALL dscal(fpar%nnr1,tw*tw,grad(1,1),1)
            CALL dscal(3*fpar%nnr1,tw,grad(1,2),1)
            CALL dscal(3*fpar%nnr1,tw,dgrad(1,1),1)
         ENDIF
      !
      ! Open-shell
      !
      ELSE
         CALL cp_dxc_compute( fpar%nnr1, rhoe, drhoe, grad, dgrad, dxc_tmp)
         !
         ! Alpha spin
         !
         DO ir=1,fpar%nnr1
            v(ir,1)=CMPLX(dxc_tmp(ir,2,1)*grad(ir,2),dxc_tmp(ir,2,1)*grad(ir,3),kind=real_8)
         ENDDO
         CALL  fwfftn(v(:,1),.FALSE.,parai%allgrp)
         CALL fftmult(v,vtemp,1,2)

         CALL zeroing(v(:,1))!,maxfft)
         DO ir=1,fpar%nnr1
            v(ir,1)=CMPLX(dxc_tmp(ir,2,1)*grad(ir,4),dxc_tmp(ir,3,1)*dgrad(ir,1),kind=real_8)
         ENDDO
         CALL  fwfftn(v(:,1),.FALSE.,parai%allgrp)
         CALL fftmult(v,vtemp,3,1)

         CALL zeroing(v(:,1))!,maxfft)
         DO ir=1,fpar%nnr1
            v(ir,1)=CMPLX(dxc_tmp(ir,3,1)*dgrad(ir,2),dxc_tmp(ir,3,1)*dgrad(ir,3),kind=real_8)
         ENDDO
         CALL  fwfftn(v(:,1),.FALSE.,parai%allgrp)
         CALL fftmult(v,vtemp,2,3)

         CALL zeroing(v(:,1))!,maxfft)
         DO ir=1,fpar%nnr1
            v(ir,1)=CMPLX(dxc_tmp(ir,4,1)*grad(ir,6),dxc_tmp(ir,4,1)*grad(ir,7),kind=real_8)
         ENDDO
         CALL  fwfftn(v(:,1),.FALSE.,parai%allgrp)
         CALL fftmult(v,vtemp,1,2)

         CALL zeroing(v(:,1))!,maxfft)
         DO ir=1,fpar%nnr1
            v(ir,1)=CMPLX(dxc_tmp(ir,4,1)*grad(ir,8),dxc_tmp(ir,5,1)*dgrad(ir,4),kind=real_8)
         ENDDO
         CALL  fwfftn(v(:,1),.FALSE.,parai%allgrp)
         CALL fftmult(v,vtemp,3,1)

         CALL zeroing(v(:,1))!,maxfft)
         DO ir=1,fpar%nnr1
            v(ir,1)=CMPLX(dxc_tmp(ir,5,1)*dgrad(ir,5),dxc_tmp(ir,5,1)*dgrad(ir,6),kind=real_8)
         ENDDO
         CALL  fwfftn(v(:,1),.FALSE.,parai%allgrp)
         CALL fftmult(v,vtemp,2,3)

         CALL zeroing(v(:,1))!,maxfft)
         DO ig=1,ncpw%nhg
            v(nzh(ig),1) = vtemp(ig)
            v(indz(ig),1) = CONJG(vtemp(ig))
         ENDDO
         CALL  invfftn(v(:,1),.FALSE.,parai%allgrp)
         DO ir=1,fpar%nnr1
            v(ir,1)=v(ir,1)+CMPLX(dxc_tmp(ir,1,1),0._real_8,kind=real_8)
         ENDDO

         !
         ! Beta Spin
         !
         CALL zeroing(vtemp)!,nhg)
         !
         DO ir=1,fpar%nnr1
            v(ir,2)=CMPLX(dxc_tmp(ir,2,2)*grad(ir,6),dxc_tmp(ir,2,2)*grad(ir,7),kind=real_8)
         ENDDO
         CALL  fwfftn(v(:,2),.FALSE.,parai%allgrp)
         CALL fftmult(v(:,2),vtemp,1,2)

         CALL zeroing(v(:,2))!,maxfft)
         DO ir=1,fpar%nnr1
            v(ir,2)=CMPLX(dxc_tmp(ir,2,2)*grad(ir,8),dxc_tmp(ir,3,2)*dgrad(ir,4),kind=real_8)
         ENDDO
         CALL  fwfftn(v(:,2),.FALSE.,parai%allgrp)
         CALL fftmult(v(:,2),vtemp,3,1)

         CALL zeroing(v(:,2))!,maxfft)
         DO ir=1,fpar%nnr1
            v(ir,2)=CMPLX(dxc_tmp(ir,3,2)*dgrad(ir,5),dxc_tmp(ir,3,2)*dgrad(ir,6),kind=real_8)
         ENDDO
         CALL  fwfftn(v(:,2),.FALSE.,parai%allgrp)
         CALL fftmult(v(:,2),vtemp,2,3)

         CALL zeroing(v(:,2))!,maxfft)
         DO ir=1,fpar%nnr1
            v(ir,2)=CMPLX(dxc_tmp(ir,4,2)*grad(ir,2),dxc_tmp(ir,4,2)*grad(ir,3),kind=real_8)
         ENDDO
         CALL  fwfftn(v(:,2),.FALSE.,parai%allgrp)
         CALL fftmult(v(:,2),vtemp,1,2)

         CALL zeroing(v(:,2))!,maxfft)
         DO ir=1,fpar%nnr1
            v(ir,2)=CMPLX(dxc_tmp(ir,4,2)*grad(ir,4),dxc_tmp(ir,5,2)*dgrad(ir,1),kind=real_8)
         ENDDO
         CALL  fwfftn(v(:,2),.FALSE.,parai%allgrp)
         CALL fftmult(v(:,2),vtemp,3,1)

         CALL zeroing(v(:,2))!,maxfft)
         DO ir=1,fpar%nnr1
            v(ir,2)=CMPLX(dxc_tmp(ir,5,2)*dgrad(ir,2),dxc_tmp(ir,5,2)*dgrad(ir,3),kind=real_8)
         ENDDO
         CALL  fwfftn(v(:,2),.FALSE.,parai%allgrp)
         CALL fftmult(v(:,2),vtemp,2,3)

         CALL zeroing(v(:,2))!,maxfft)
         DO ig=1,ncpw%nhg
            v(nzh(ig),2) = vtemp(ig)
            v(indz(ig),2) = CONJG(vtemp(ig))
         ENDDO
         CALL  invfftn(v(:,2),.FALSE.,parai%allgrp)

         DO ir=1,fpar%nnr1
            v(ir,2)=v(ir,2)+CMPLX(dxc_tmp(ir,1,2),0._real_8,kind=real_8)
         ENDDO
      ENDIF
    ELSE
      IF (.NOT.cntl%tlsd) THEN
         ! no spin polarisation --> a-dens = b-dens = 1/2 * dens      
         CALL dscal(fpar%nnr1,hf,rhoe,1)
         CALL dscal(fpar%nnr1,hf,drhoe,1)
         IF (cntl%tgc) THEN
            CALL dscal(fpar%nnr1,qt,grad(1,1),1)
            CALL dscal(3*fpar%nnr1,hf,grad(1,2),1)
            CALL dscal(3*fpar%nnr1,hf,dgrad(1,1),1)
         ENDIF
         ! ==--------------------------------------------------------------==
         ! ==  SLATER LDA only  -  not spin polarised                      ==
         ! ==--------------------------------------------------------------==
         IF (.NOT.cntl%tgc) THEN
            DO ir=1,fpar%nnr1
               IF (rhoe(ir,1).LE.0._real_8) THEN
                  v(ir,1)=0._real_8
               ELSE
                  CALL slater_x(rhoe(ir,1),lda)
                  v(ir,1)=lda*drhoe(ir,1)
               ENDIF
            ENDDO
            GOTO 99
            ! ==--------------------------------------------------------------==
            ! ==  BECKE88 ONLY - not spin polarised                           ==
            ! ==--------------------------------------------------------------==
            ! no spin polarisation --> a-dens = b-dens = 1/2 * dens      
         ELSEIF (func1%mgcx == mgcx_is_becke88 .AND. func1%mgcc == mgcc_is_skipped) THEN
            DO ir=1,fpar%nnr1
               IF (rhoe(ir,1).LE.0._real_8) THEN
               ELSEIF(rhoe(ir,1).GT.0._real_8.AND.(rhoe(ir,1).LT.cntr%gceps.OR.&
                    grad(IR,1).LE.0._real_8)) THEN                           
                  CALL slater_x(rhoe(ir,1),lda)
                  tmp1(ir)=lda*drhoe(ir,1)
               ELSE
                  CALL b88_x(rhoe(ir,1),grad(ir,1),dfdn,dfdg,&
                       d2fdn2,d2fdndg,d2fdg2,lda)
                  skaa=grad(ir,2)*dgrad(ir,1)+grad(ir,3)*dgrad(ir,2)&
                       +grad(ir,4)*dgrad(ir,3)
                  sgrad=SQRT(grad(ir,1))
                  skalar=skaa/sgrad
                  tmp1(ir)=(d2fdn2+lda)*drhoe(ir,1)+d2fdndg*skalar
                  tmp2(ir)=1._real_8/sgrad*(d2fdndg*drhoe(ir,1)+&
                       d2fdg2*skalar-dfdg*skalar/sgrad)
                  tmp3(ir)=dfdg/sgrad
               ENDIF
            ENDDO
            ! ==--------------------------------------------------------------==
            ! ==  BLYP  -  not spin polarised                                 ==
            ! ==--------------------------------------------------------------==
            ! no spin polarisation --> a-dens = b-dens = 1/2 * dens
         ELSEIF (func1%mgcx == mgcx_is_becke88 .AND. func1%mgcc == mgcc_is_lyp) THEN
            DO ir=1,fpar%nnr1
               IF (rhoe(ir,1).GT.1.e-24_real_8.AND.(rhoe(ir,1).LT.cntr%gceps.OR.&
                    grad(IR,1).LE.1.e-48_real_8)) THEN       
                  CALL lyp88_c_loc(rhoe(ir,1),rhoe(ir,1),loca,locb)
                  tmp1(ir)=(loca+locb)*drhoe(ir,1)
               ELSEIF (rhoe(ir,1).GT.1.e-24_real_8) THEN
                  CALL lyp88_c(rhoe(ir,1),rhoe(ir,1),&
                       grad(ir,1),grad(ir,2),grad(ir,3),grad(ir,4),&
                       grad(IR,1),grad(IR,2),grad(IR,3),grad(IR,4),&
                       dfdga,nana,nanb,naga,granb,gagb,nagb,nagrb,ganb)
                  CALL b88_x(rhoe(ir,1),grad(ir,1),dfdn,dfdg,&
                       d2fdn2,d2fdndg,d2fdg2,lda)
                  skaa=grad(ir,2)*dgrad(ir,1)+grad(ir,3)*dgrad(ir,2)&
                       +grad(ir,4)*dgrad(ir,3)
                  sgrad=SQRT(grad(ir,1))
                  skalar=skaa/sgrad
                  tmp1(ir)=(nana+d2fdn2+lda)*drhoe(ir,1)&
                       +nanb*drhoe(ir,1)+naga*2._real_8*skaa+nagb*skaa&
                       +nagrb*2._real_8*skaa+nagb*skaa+d2fdndg*skalar
                  tmp2(ir)=2._real_8*(naga*drhoe(ir,1)+granb*drhoe(ir,1))&
                       +1._real_8/sgrad*(d2fdndg*drhoe(ir,1)&
                       +skalar*(d2fdg2-dfdg/sgrad))
                  tmp3(ir)=dfdga*2._real_8+dfdg/sgrad
                  tmp4(ir)=nagb*drhoe(ir,1)+ganb*drhoe(ir,1)
                  tmp5(ir)=gagb
               ENDIF
            ENDDO
            ! ==--------------------------------------------------------------==
            ! ==  BP  -  not spin polarised                                   ==
            ! ==--------------------------------------------------------------==
            ! no spin polarisation --> a-dens = b-dens = 1/2 * dens
         ELSEIF (func1%mgcx == mgcx_is_becke88 .AND. func1%mgcc == mgcc_is_perdew86) THEN
            DO ir=1,fpar%nnr1
               IF (rhoe(ir,1).GT.1.e-24_real_8.AND.(rhoe(ir,1).LT.cntr%gceps.OR.&
                    grad(IR,1).LE.1.e-48_real_8)) THEN
                  CALL slater_x(rhoe(ir,1),lda)
                  CALL p86_c_loc(rhoe(ir,1),rhoe(ir,1),p86a,p86b)
                  tmp1(ir)=(lda+p86a)*drhoe(ir,1)+p86b*drhoe(ir,1)
               ELSEIF (rhoe(ir,1).GT.1.e-24_real_8) THEN
                  CALL p86_c(rhoe(ir,1),rhoe(ir,1),grad(ir,1),&
                       grad(IR,2),grad(IR,3),grad(IR,4),&
                       grad(IR,1),grad(IR,2),grad(IR,3),&
                       grad(ir,4),nana,naga,gaga,nask,gask,sk,ga,sksk,&
                       gagb,nagb,skgb,sknb,ganb,nanb)
                  CALL b88_x(rhoe(ir,1),grad(ir,1),dfdn,dfdg,d2fdn2,&
                       d2fdndg,d2fdg2,LDA)
                  skaa=(grad(ir,2)*dgrad(ir,1)+grad(ir,3)*dgrad(ir,2)&
                       +grad(ir,4)*dgrad(ir,3))
                  skab=skaa
                  skba=skaa
                  skbb=skaa
                  sgrad=SQRT(grad(ir,1))
                  sgradb=sgrad
                  skalar=skaa/sgrad
                  tmp1(ir)=(nana+d2fdn2+lda)*drhoe(ir,1)+&
                       nanb*drhoe(IR,1)&
                       +naga*skalar+nask*(skba+skab)+nagb*skbb/sgradb&
                       +d2fdndg*skalar
                  tmp2(ir)=(naga*drhoe(ir,1)+ganb*drhoe(ir,1)+gask*&
                       (skba+skab)+gaga*skalar+gagb*skbb/sgradb&
                       -ga*skalar/sgrad)/sgrad&
                       +(d2fdndg*drhoe(IR,1)+skalar*&
                       (d2fdg2-dfdg/sgrad))/sgrad
                  tmp3(ir)=(ga+dfdg)/sgrad
                  tmp4(ir)=nask*drhoe(ir,1)+sknb*drhoe(ir,1)&
                       +gask*skalar+skgb*skbb/sgradb+sksk*(skba+skab)
                  tmp5(ir)=sk
               ENDIF
            ENDDO
            ! ==--------------------------------------------------------------==
            ! ==  PBE  -  not spin polarised
            ! ==--------------------------------------------------------------==
            ! no spin polarisation --> a-dens = b-dens = 1/2 * dens
         ELSEIF (func1%mgcx == mgcx_is_pbex .AND. func1%mgcc == mgcc_is_pbec) THEN
            DO ir=1,fpar%nnr1
               IF (rhoe(ir,1).GT.1.e-24_real_8.AND.(rhoe(ir,1).LT.cntr%gceps.OR.&
                    grad(IR,1).LE.1.e-48_real_8)) THEN  
                  CALL slater_x(rhoe(ir,1),lda)
                  CALL pbe96_c_loc(rhoe(ir,1),rhoe(ir,1),nana,nanb)
                  tmp1(ir)=(lda+nana)*drhoe(ir,1)+nanb*drhoe(ir,1)
               ELSEIF (rhoe(ir,1).GT.1.e-24_real_8) THEN
                  CALL pbe96_x(rhoe(ir,1),grad(ir,1),dfdg,d2fdn2,&
                       d2fdndg,d2fdg2)
                  CALL pbe96_c(rhoe(ir,1),rhoe(ir,1),grad(ir,1),&
                       grad(IR,2),grad(IR,3),grad(IR,4),&
                       grad(IR,1),grad(IR,2),grad(IR,3),&
                       grad(ir,4),nana,naga,gaga,nask,gask,sk,ga,sksk,gagb,&
                       nagb,skgb,sknb,ganb,nanb)
                  skaa=(grad(ir,2)*dgrad(ir,1)+grad(ir,3)*dgrad(ir,2)&
                       +grad(ir,4)*dgrad(ir,3))
                  skab=skaa
                  skba=skaa
                  skbb=skaa
                  tmp1(ir)=(nana+d2fdn2)*drhoe(ir,1)+nanb*&
                       drhoe(IR,1)+naga*2._real_8*skaa+nask*(skba+skab)+&
                       nagb*2._real_8*skbb+d2fdndg*2._real_8*skaa
                  tmp2(ir)=2._real_8*(naga*drhoe(ir,1)+ganb*drhoe(ir,1)+&
                       gask*(skba+skab)+gaga*2._real_8*skaa+gagb*2._real_8*skbb&
                       +d2fdndg*drhoe(IR,1)+2._real_8*skaa*d2fdg2)
                  tmp3(ir)=2._real_8*(ga+dfdg)
                  tmp4(ir)=nask*drhoe(ir,1)+sknb*drhoe(ir,1)+&
                       gask*2._real_8*skaa+skgb*2._real_8*skbb+sksk*(skba+skab)
                  tmp5(ir)=sk
               ENDIF
            ENDDO
         ELSE
            CALL stopgm(procedureN,'Functional not implemented',& 
                 __LINE__,__FILE__)
         ENDIF

         ! ==--------------------------------------------------------------==
         ! CALCULATION OF DIVERGENCE IN RECIPROCAL SPACE:
         ! real space               G-space
         ! dot(nabla,gradient) --> dot(i*G-Vector,gradient)
         ! tmp2: terms to be multiplied with gradient of a-density
         ! tmp3: terms to be multiplied with gradient of a-linres-density
         ! tmp4: terms to be multiplied with gradient of b-density
         ! tmp5: terms to be multiplied with gradient of b-linres-density
         ! ==--------------------------------------------------------------==
         ! FFTMULT:calculates dot product of i*G-Vector with gradient 
         DO ir=1,fpar%nnr1
            v(ir,1)=CMPLX(tmp2(ir)*grad(ir,2),tmp2(ir)*grad(ir,3),kind=real_8)
         ENDDO
         CALL  fwfftn(v(:,1),.FALSE.,parai%allgrp)
         CALL fftmult(v,vtemp,1,2)

         CALL zeroing(v(:,1))!,maxfft)
         DO ir=1,fpar%nnr1
            v(ir,1)=CMPLX(tmp2(ir)*grad(ir,4),tmp3(ir)*dgrad(ir,1),kind=real_8)
         ENDDO
         CALL  fwfftn(v(:,1),.FALSE.,parai%allgrp)
         CALL fftmult(v,vtemp,3,1)

         CALL zeroing(v(:,1))!,maxfft)
         DO ir=1,fpar%nnr1
            v(ir,1)=CMPLX(tmp3(ir)*dgrad(ir,2),tmp3(ir)*dgrad(ir,3),kind=real_8)
         ENDDO
         CALL  fwfftn(v(:,1),.FALSE.,parai%allgrp)
         CALL fftmult(v,vtemp,2,3)

         CALL zeroing(v(:,1))!,maxfft)
         DO ir=1,fpar%nnr1
            v(ir,1)=CMPLX(tmp4(ir)*grad(ir,2),tmp4(ir)*grad(ir,3),kind=real_8)
         ENDDO
         CALL  fwfftn(v(:,1),.FALSE.,parai%allgrp)
         CALL fftmult(v,vtemp,1,2)

         CALL zeroing(v(:,1))!,maxfft)
         DO ir=1,fpar%nnr1
            v(ir,1)=CMPLX(tmp4(ir)*grad(ir,4),tmp5(ir)*dgrad(ir,1),kind=real_8)
         ENDDO
         CALL  fwfftn(v(:,1),.FALSE.,parai%allgrp)
         CALL fftmult(v,vtemp,3,1)

         CALL zeroing(v(:,1))!,maxfft)
         DO ir=1,fpar%nnr1
            v(ir,1)=CMPLX(tmp5(ir)*dgrad(ir,2),tmp5(ir)*dgrad(ir,3),kind=real_8)
         ENDDO
         CALL  fwfftn(v(:,1),.FALSE.,parai%allgrp)
         CALL fftmult(v,vtemp,2,3)

         CALL zeroing(v(:,1))!,maxfft)
         DO ig=1,ncpw%nhg
            v(nzh(ig),1) = vtemp(ig)
            v(indz(ig),1) = CONJG(vtemp(ig))
         ENDDO

         CALL  invfftn(v(:,1),.FALSE.,parai%allgrp)
         DO ir=1,fpar%nnr1
            v(ir,1)=v(ir,1)+CMPLX(tmp1(ir),0._real_8,kind=real_8)
         ENDDO
         CALL dscal(fpar%nnr1,tw,rhoe,1)
         CALL dscal(fpar%nnr1,tw,drhoe,1)
         IF (cntl%tgc) THEN
            CALL dscal(fpar%nnr1,tw*tw,grad(1,1),1)
            CALL dscal(3*fpar%nnr1,tw,grad(1,2),1)
            CALL dscal(3*fpar%nnr1,tw,dgrad(1,1),1)
         ENDIF

      ELSE
         ! SPIN POLARISED
         ! ==--------------------------------------------------------------==
         ! ==  SLATER LDA only
         ! ==--------------------------------------------------------------==
         IF (.NOT.cntl%tgc) THEN
            DO ir=1,fpar%nnr1
               IF (rhoe(ir,1).GT.0._real_8) THEN
                  CALL slater_x(rhoe(ir,1),lda)
                  v(ir,1)=lda*drhoe(ir,1)
               ENDIF
               IF (rhoe(ir,2).GT.0._real_8) THEN
                  CALL slater_x(rhoe(ir,2),lda)
                  v(ir,2)=lda*drhoe(ir,2)
               ENDIF
            ENDDO
            GOTO 99
            ! Alpha Spin
            ! ==--------------------------------------------------------------==
            ! ==  BECKE88 ONLY: individual terms
            ! ==--------------------------------------------------------------==
         ELSEIF (func1%mgcx == mgcx_is_becke88 .AND. func1%mgcc == mgcc_is_skipped) THEN
            DO ir=1,fpar%nnr1
               IF (rhoe(ir,1).LE.0._real_8) THEN
               ELSEIF(rhoe(ir,1).GT.0._real_8.AND.(rhoe(ir,1).LT.cntr%gceps.OR.&
                    grad(IR,1).LE.0._real_8)) THEN              
                  CALL slater_x(rhoe(ir,1),lda)
                  tmp1(ir)=lda*drhoe(ir,1)
               ELSE
                  CALL b88_x(rhoe(ir,1),grad(ir,1),dfdn,dfdg,&
                       d2fdn2,d2fdndg,d2fdg2,lda)
                  skaa=grad(ir,2)*dgrad(ir,1)+grad(ir,3)*dgrad(ir,2)&
                       +grad(ir,4)*dgrad(ir,3)
                  sgrad=SQRT(grad(ir,1))
                  skalar=skaa/sgrad
                  tmp1(ir)=(d2fdn2+lda)*drhoe(ir,1)+d2fdndg*skalar
                  tmp2(ir)=1._real_8/sgrad*(d2fdndg*drhoe(ir,1)+d2fdg2*skalar&
                       -dfdg*skalar/sgrad)
                  tmp3(ir)=dfdg/sgrad
               ENDIF
            ENDDO

            ! ==--------------------------------------------------------------==
            ! ==  BLYP: individual terms
            ! ==--------------------------------------------------------------==
         ELSEIF (func1%mgcx == mgcx_is_becke88 .AND. func1%mgcc == mgcc_is_lyp) THEN
            DO ir=1,fpar%nnr1
               IF (rhoe(ir,1).GT.1.e-24_real_8.AND.(rhoe(ir,1).LT.cntr%gceps.OR.&
                    grad(IR,1).LE.1.e-48_real_8)) THEN       
                  CALL lyp88_c_loc(rhoe(ir,1),rhoe(ir,2),loca,locb)
                  tmp1(ir)=loca*drhoe(ir,1)+locb*drhoe(ir,2)
               ELSEIF (rhoe(ir,1).GT.1.e-24_real_8) THEN
                  CALL lyp88_c(rhoe(ir,1),rhoe(ir,2),grad(ir,1),grad(ir,2)&
                       ,grad(ir,3),grad(ir,4),grad(ir,5),grad(ir,6),grad(ir,7)&
                       ,grad(ir,8),dfdga,nana,nanb,naga,granb,gagb,nagb,&
                       nagrb,ganb)
                  CALL b88_x(rhoe(ir,1),grad(ir,1),dfdn,dfdg,d2fdn2,d2fdndg&
                       ,d2fdg2,LDA)
                  skaa=grad(ir,2)*dgrad(ir,1)+grad(ir,3)*dgrad(ir,2)&
                       +grad(ir,4)*dgrad(ir,3)
                  skab=grad(ir,2)*dgrad(ir,4)+grad(ir,3)*dgrad(ir,5)&
                       +grad(ir,4)*dgrad(ir,6)
                  skba=grad(ir,6)*dgrad(ir,1)+grad(ir,7)*dgrad(ir,2)&
                       +grad(ir,8)*dgrad(ir,3)
                  skbb=grad(ir,6)*dgrad(ir,4)+grad(ir,7)*dgrad(ir,5)&
                       +grad(ir,8)*dgrad(ir,6)
                  sgrad=SQRT(grad(ir,1))
                  skalar=skaa/sgrad
                  tmp1(ir)=(nana+d2fdn2+lda)*drhoe(ir,1)&
                       +nanb*drhoe(ir,2)+naga*2._real_8*skaa+nagb*skab&
                       +nagrb*2._real_8*skbb+nagb*skba+d2fdndg*skalar
                  tmp2(ir)=2._real_8*(naga*drhoe(ir,1)+granb*drhoe(ir,2))&
                       +1._real_8/sgrad*(d2fdndg*drhoe(ir,1)&
                       +skalar*(d2fdg2-dfdg/sgrad))
                  tmp3(ir)=dfdga*2._real_8+dfdg/sgrad
                  tmp4(ir)=nagb*drhoe(ir,1)+ganb*drhoe(ir,2)
                  tmp5(ir)=gagb
               ENDIF
            ENDDO
            ! ==--------------------------------------------------------------==
            ! ==  BP
            ! ==--------------------------------------------------------------==
         ELSEIF (func1%mgcx == mgcx_is_becke88 .AND. func1%mgcc == mgcc_is_perdew86) THEN
            DO ir=1,fpar%nnr1
               IF (rhoe(ir,1).GT.1.e-24_real_8.AND.(rhoe(ir,1).LT.cntr%gceps.OR.&
                    grad(IR,1).LE.1.e-48_real_8)) THEN  
                  CALL slater_x(rhoe(ir,1),lda)
                  CALL p86_c_loc(rhoe(ir,1),rhoe(ir,2),p86a,p86b)
                  tmp1(ir)=(lda+p86a)*drhoe(ir,1)+p86b*drhoe(ir,2)
               ELSEIF (rhoe(ir,1).GT.1.e-24_real_8) THEN
                  CALL p86_c(rhoe(ir,1),rhoe(ir,2),grad(ir,1),grad(ir,2)&
                       ,grad(ir,3),grad(ir,4),grad(ir,5),grad(ir,6),grad(ir,7)&
                       ,grad(ir,8),nana,naga,gaga,nask,gask,sk,ga,sksk,gagb&
                       ,nagb,skgb,sknb,ganb,nanb)
                  CALL b88_x(rhoe(ir,1),grad(ir,1),dfdn,dfdg,d2fdn2,d2fdndg&
                       ,d2fdg2,LDA)
                  skaa=grad(ir,2)*dgrad(ir,1)+grad(ir,3)*dgrad(ir,2)&
                       +grad(ir,4)*dgrad(ir,3)
                  skab=grad(ir,2)*dgrad(ir,4)+grad(ir,3)*dgrad(ir,5)&
                       +grad(ir,4)*dgrad(ir,6)
                  skba=grad(ir,6)*dgrad(ir,1)+grad(ir,7)*dgrad(ir,2)&
                       +grad(ir,8)*dgrad(ir,3)
                  skbb=grad(ir,6)*dgrad(ir,4)+grad(ir,7)*dgrad(ir,5)&
                       +grad(ir,8)*dgrad(ir,6)
                  sgrad=SQRT(grad(ir,1))
                  sgradb=SQRT(grad(ir,5))
                  skalar=skaa/sgrad
                  tmp1(ir)=(nana+d2fdn2+lda)*drhoe(ir,1)+nanb*drhoe(ir,2)&
                       +naga*skalar+nask*(skba+skab)+nagb*skbb/sgradb&
                       +d2fdndg*skalar
                  tmp2(ir)=(naga*drhoe(ir,1)+ganb*drhoe(ir,2)+gask*&
                       (skba+skab)+gaga*skalar+gagb*skbb/sgradb&
                       -ga*skalar/sgrad)/sgrad&
                       +1._real_8/sgrad*(d2fdndg*drhoe(IR,1)+skalar*&
                       (d2fdg2-dfdg/sgrad))
                  tmp3(ir)=(ga+dfdg)/sgrad
                  tmp4(ir)=nask*drhoe(ir,1)+sknb*drhoe(ir,2)+gask*skalar+&
                       skgb*skbb/sgradb+sksk*(skba+skab)
                  tmp5(ir)=sk
               ENDIF
            ENDDO
            ! ==--------------------------------------------------------------==
            ! ==  PBE
            ! ==--------------------------------------------------------------==
         ELSEIF (func1%mgcx == mgcx_is_pbex .AND. func1%mgcc == mgcc_is_pbec) THEN
            DO ir=1,fpar%nnr1
               IF (rhoe(ir,1).GT.1.e-24_real_8.AND.(rhoe(ir,1).LT.cntr%gceps.OR.&
                    grad(IR,1).LE.1.e-48_real_8)) THEN  
                  CALL slater_x(rhoe(ir,1),lda)
                  CALL pbe96_c_loc(rhoe(ir,1),rhoe(ir,2),nana,nanb)
                  tmp1(ir)=(lda+nana)*drhoe(ir,1)+nanb*drhoe(ir,2)
               ELSEIF (rhoe(ir,1).GT.1.e-24_real_8) THEN
                  CALL pbe96_x(rhoe(ir,1),grad(ir,1),dfdg,d2fdn2,d2fdndg,&
                       d2fdg2)
                  CALL pbe96_c(rhoe(ir,1),rhoe(ir,2),grad(ir,1),grad(ir,2)&
                       ,grad(ir,3),grad(ir,4),grad(ir,5),grad(ir,6),grad(ir,7)&
                       ,grad(ir,8),nana,naga,gaga,nask,gask,sk,ga,sksk,gagb&
                       ,nagb,skgb,sknb,ganb,nanb)
                  skaa=grad(ir,2)*dgrad(ir,1)+grad(ir,3)*dgrad(ir,2)&
                       +grad(ir,4)*dgrad(ir,3)
                  skab=grad(ir,2)*dgrad(ir,4)+grad(ir,3)*dgrad(ir,5)&
                       +grad(ir,4)*dgrad(ir,6)
                  skba=grad(ir,6)*dgrad(ir,1)+grad(ir,7)*dgrad(ir,2)&
                       +grad(ir,8)*dgrad(ir,3)
                  skbb=grad(ir,6)*dgrad(ir,4)+grad(ir,7)*dgrad(ir,5)&
                       +grad(ir,8)*dgrad(ir,6)
                  tmp1(ir)=(nana+d2fdn2)*drhoe(ir,1)+nanb*drhoe(ir,2)&
                       +naga*2._real_8*skaa+nask*(skba+skab)+nagb*2._real_8*skbb&
                       +d2fdndg*2._real_8*skaa
                  tmp2(ir)=2._real_8*(naga*drhoe(ir,1)+ganb*drhoe(ir,2)+gask*&
                       (skba+skab)+gaga*2._real_8*skaa+gagb*2._real_8*skbb&
                       +d2fdndg*drhoe(IR,1)+2._real_8*skaa*d2fdg2)
                  tmp3(ir)=2._real_8*(ga+dfdg)
                  tmp4(ir)=nask*drhoe(ir,1)+sknb*drhoe(ir,2)+gask*2._real_8*skaa+&
                       skgb*2._real_8*skbb+sksk*(skba+skab)
                  tmp5(ir)=sk
               ENDIF
            ENDDO
         ELSE
            CALL stopgm(procedureN,'Functional not implemented',& 
                 __LINE__,__FILE__)
         ENDIF

         ! ==--------------------------------------------------------------==
         ! CALCULATION OF DIVERGENCE IN RECIPROCAL SPACE:
         ! real space               G-space
         ! dot(nabla,gradient) --> dot(i*G-Vector,gradient)
         ! tmp2: terms to be multiplied with gradient of a-density
         ! tmp3: terms to be multiplied with gradient of a-linres-density
         ! tmp4: terms to be multiplied with gradient of b-density
         ! tmp5: terms to be multiplied with gradient of b-linres-density
         ! ==--------------------------------------------------------------==
         ! FFTMULT:calculates dot product ofi*G-Vector with gradient 
         DO ir=1,fpar%nnr1
            v(ir,1)=CMPLX(tmp2(ir)*grad(ir,2),tmp2(ir)*grad(ir,3),kind=real_8)
         ENDDO
         CALL  fwfftn(v(:,1),.FALSE.,parai%allgrp)
         CALL fftmult(v,vtemp,1,2)

         CALL zeroing(v(:,1))!,maxfft)
         DO ir=1,fpar%nnr1
            v(ir,1)=CMPLX(tmp2(ir)*grad(ir,4),tmp3(ir)*dgrad(ir,1),kind=real_8)
         ENDDO
         CALL  fwfftn(v(:,1),.FALSE.,parai%allgrp)
         CALL fftmult(v,vtemp,3,1)

         CALL zeroing(v(:,1))!,maxfft)
         DO ir=1,fpar%nnr1
            v(ir,1)=CMPLX(tmp3(ir)*dgrad(ir,2),tmp3(ir)*dgrad(ir,3),kind=real_8)
         ENDDO
         CALL  fwfftn(v(:,1),.FALSE.,parai%allgrp)
         CALL fftmult(v,vtemp,2,3)

         CALL zeroing(v(:,1))!,maxfft)
         DO ir=1,fpar%nnr1
            v(ir,1)=CMPLX(tmp4(ir)*grad(ir,6),tmp4(ir)*grad(ir,7),kind=real_8)
         ENDDO
         CALL  fwfftn(v(:,1),.FALSE.,parai%allgrp)
         CALL fftmult(v,vtemp,1,2)

         CALL zeroing(v(:,1))!,maxfft)
         DO ir=1,fpar%nnr1
            v(ir,1)=CMPLX(tmp4(ir)*grad(ir,8),tmp5(ir)*dgrad(ir,4),kind=real_8)
         ENDDO
         CALL  fwfftn(v(:,1),.FALSE.,parai%allgrp)
         CALL fftmult(v,vtemp,3,1)

         CALL zeroing(v(:,1))!,maxfft)
         DO ir=1,fpar%nnr1
            v(ir,1)=CMPLX(tmp5(ir)*dgrad(ir,5),tmp5(ir)*dgrad(ir,6),kind=real_8)
         ENDDO
         CALL  fwfftn(v(:,1),.FALSE.,parai%allgrp)
         CALL fftmult(v,vtemp,2,3)

         CALL zeroing(v(:,1))!,maxfft)
         DO ig=1,ncpw%nhg
            v(nzh(ig),1) = vtemp(ig)
            v(indz(ig),1) = CONJG(vtemp(ig))
         ENDDO
         CALL  invfftn(v(:,1),.FALSE.,parai%allgrp)
         DO ir=1,fpar%nnr1
            v(ir,1)=v(ir,1)+CMPLX(tmp1(ir),0._real_8,kind=real_8)
         ENDDO

         ! Beta Spin
         CALL zeroing(vtemp)!,nhg)
         CALL zeroing(tmp1)!,nnr1)
         CALL zeroing(tmp2)!,nnr1)
         CALL zeroing(tmp3)!,nnr1)
         CALL zeroing(tmp4)!,nnr1)
         CALL zeroing(tmp5)!,nnr1)

         ! ==--------------------------------------------------------------==
         ! ==  BECKE88 ONLY: individual terms
         ! ==--------------------------------------------------------------==
         IF (func1%mgcx == mgcx_is_becke88 .AND. func1%mgcc == mgcc_is_skipped) THEN
            DO ir=1,fpar%nnr1
               IF (rhoe(ir,2).LE.0._real_8) THEN
               ELSEIF(rhoe(ir,2).GT.0._real_8.AND.(rhoe(ir,2).LT.cntr%gceps.OR.&
                    grad(IR,5).LE.0._real_8)) THEN              
                  CALL slater_x(rhoe(ir,2),lda)
                  tmp1(ir)=lda*drhoe(ir,2)
               ELSE
                  CALL b88_x(rhoe(ir,2),grad(ir,5),dfdn,dfdg,&
                       d2fdn2,d2fdndg,d2fdg2,lda)
                  skaa=grad(ir,6)*dgrad(ir,4)+grad(ir,7)*dgrad(ir,5)&
                       +grad(ir,8)*dgrad(ir,6)
                  sgrad=SQRT(grad(ir,5))
                  skalar=skaa/sgrad
                  tmp1(ir)=(d2fdn2+lda)*drhoe(ir,2)+d2fdndg*skalar
                  tmp2(ir)=1._real_8/sgrad*(d2fdndg*drhoe(ir,2)+d2fdg2*skalar&
                       -dfdg*skalar/sgrad)
                  tmp3(ir)=dfdg/sgrad
               ENDIF
            ENDDO

            ! ==--------------------------------------------------------------==
            ! ==  BLYP: individual terms
            ! ==--------------------------------------------------------------==
         ELSEIF (func1%mgcx == mgcx_is_becke88 .AND. func1%mgcc == mgcc_is_lyp) THEN
            DO ir=1,fpar%nnr1
               IF (rhoe(ir,2).GT.1.e-24_real_8.AND.(rhoe(ir,2).LT.cntr%gceps.OR.&
                    grad(IR,5).LE.1.e-48_real_8)) THEN       
                  CALL lyp88_c_loc(rhoe(ir,2),rhoe(ir,1),loca,locb)
                  tmp1(ir)=loca*drhoe(ir,2)+locb*drhoe(ir,1)
               ELSEIF (rhoe(ir,2).GT.1.e-24_real_8) THEN
                  CALL lyp88_c(rhoe(ir,2),rhoe(ir,1),grad(ir,5),grad(ir,6)&
                       ,grad(ir,7),grad(ir,8),grad(ir,1),grad(ir,2),grad(ir,3)&
                       ,grad(ir,4),dfdga,nana,nanb,naga,granb,gagb,nagb,nagrb&
                       ,ganb)
                  CALL b88_x(rhoe(ir,2),grad(ir,5),dfdn,dfdg,d2fdn2,&
                       d2fdndg,d2fdg2,lda)
                  skaa=grad(ir,6)*dgrad(ir,4)+grad(ir,7)*dgrad(ir,5)&
                       +grad(ir,8)*dgrad(ir,6)
                  skab=grad(ir,6)*dgrad(ir,1)+grad(ir,7)*dgrad(ir,2)&
                       +grad(ir,8)*dgrad(ir,3)
                  skba=grad(ir,2)*dgrad(ir,4)+grad(ir,3)*dgrad(ir,5)&
                       +grad(ir,4)*dgrad(ir,6)
                  skbb=grad(ir,2)*dgrad(ir,1)+grad(ir,3)*dgrad(ir,2)&
                       +grad(ir,4)*dgrad(ir,3)
                  sgrad=SQRT(grad(ir,5))
                  skalar=skaa/sgrad
                  tmp1(ir)=(nana+d2fdn2+lda)*drhoe(ir,2)&
                       +nanb*drhoe(ir,1)+naga*2._real_8*skaa+nagb*skab&
                       +nagrb*2._real_8*skbb+nagb*skba+d2fdndg*skalar
                  tmp2(ir)=2._real_8*(naga*drhoe(ir,2)+granb*drhoe(ir,1))&
                       +1._real_8/sgrad*(d2fdndg*drhoe(ir,2)&
                       +skalar*(d2fdg2-dfdg/sgrad))
                  tmp3(ir)=dfdga*2._real_8+dfdg/sgrad
                  tmp4(ir)=nagb*drhoe(ir,2)+ganb*drhoe(ir,1)
                  tmp5(ir)=gagb
               ENDIF
            ENDDO
            ! ==--------------------------------------------------------------==
            ! ==  BP
            ! ==--------------------------------------------------------------==
         ELSEIF (func1%mgcx == mgcx_is_becke88 .AND. func1%mgcc == mgcc_is_perdew86) THEN
            DO ir=1,fpar%nnr1
               IF (rhoe(ir,2).GT.1.e-24_real_8.AND.(rhoe(ir,2).LT.cntr%gceps.OR.&
                    grad(IR,5).LE.1.e-48_real_8)) THEN
                  CALL  slater_x(rhoe(ir,2),lda)
                  CALL p86_c_loc(rhoe(ir,2),rhoe(ir,1),p86a,p86b)
                  tmp1(ir)=(lda+p86a)*drhoe(ir,2)+p86b*drhoe(ir,1)
               ELSEIF (rhoe(ir,2).GT.1.e-24_real_8) THEN
                  CALL p86_c(rhoe(ir,2),rhoe(ir,1),grad(ir,5),grad(ir,6)&
                       ,grad(ir,7),grad(ir,8),grad(ir,1),grad(ir,2),grad(ir,3)&
                       ,grad(ir,4),nana,naga,gaga,nask,gask,sk,ga,sksk,gagb&
                       ,nagb,skgb,sknb,ganb,nanb)
                  CALL b88_x(rhoe(ir,2),grad(ir,5),dfdn,dfdg,d2fdn2,d2fdndg&
                       ,d2fdg2,LDA)
                  skaa=grad(ir,6)*dgrad(ir,4)+grad(ir,7)*dgrad(ir,5)&
                       +grad(ir,8)*dgrad(ir,6)
                  skab=grad(ir,6)*dgrad(ir,1)+grad(ir,7)*dgrad(ir,2)&
                       +grad(ir,8)*dgrad(ir,3)
                  skba=grad(ir,2)*dgrad(ir,4)+grad(ir,3)*dgrad(ir,5)&
                       +grad(ir,4)*dgrad(ir,6)
                  skbb=grad(ir,2)*dgrad(ir,1)+grad(ir,3)*dgrad(ir,2)&
                       +grad(ir,4)*dgrad(ir,3)
                  sgrad=SQRT(grad(ir,5))
                  sgradb=SQRT(grad(ir,1))
                  skalar=skaa/sgrad
                  tmp1(ir)=(nana+d2fdn2+lda)*drhoe(ir,2)+nanb*drhoe(ir,1)&
                       +naga*skalar+nask*(skba+skab)+nagb*skbb/sgradb&
                       +d2fdndg*skalar
                  tmp2(ir)=(naga*drhoe(ir,2)+ganb*drhoe(ir,1)+gask*&
                       (skba+skab)+gaga*skalar+gagb*skbb/sgradb&
                       -ga*skalar/sgrad)/sgrad&
                       +1._real_8/sgrad*(d2fdndg*drhoe(IR,2)+skalar*&
                       (d2fdg2-dfdg/sgrad))
                  tmp3(ir)=(ga+dfdg)/sgrad
                  tmp4(ir)=nask*drhoe(ir,2)+sknb*drhoe(ir,1)+gask*skalar+&
                       skgb*skbb/sgradb+sksk*(skba+skab)
                  tmp5(ir)=sk
               ENDIF
            ENDDO
            ! ==--------------------------------------------------------------==
            ! ==  PBE
            ! ==--------------------------------------------------------------==
         ELSEIF (func1%mgcx == mgcx_is_pbex .AND. func1%mgcc == mgcc_is_pbec) THEN
            DO ir=1,fpar%nnr1
               IF (rhoe(ir,2).GT.1.e-24_real_8.AND.(rhoe(ir,2).LT.cntr%gceps.OR.&
                    grad(IR,5).LE.1.e-48_real_8)) THEN  
                  CALL slater_x(rhoe(ir,2),lda)
                  CALL pbe96_c_loc(rhoe(ir,2),rhoe(ir,1),nana,nanb)
                  tmp1(ir)=(lda+nana)*drhoe(ir,2)+nanb*drhoe(ir,1)
               ELSEIF (rhoe(ir,2).GT.1.e-24_real_8) THEN
                  CALL pbe96_x(rhoe(ir,2),grad(ir,5),dfdg,d2fdn2,d2fdndg,&
                       d2fdg2)
                  CALL pbe96_c(rhoe(ir,2),rhoe(ir,1),grad(ir,5),grad(ir,6)&
                       ,grad(ir,7),grad(ir,8),grad(ir,1),grad(ir,2),grad(ir,3)&
                       ,grad(ir,4),nana,naga,gaga,nask,gask,sk,ga,sksk,gagb&
                       ,nagb,skgb,sknb,ganb,nanb)
                  skaa=grad(ir,6)*dgrad(ir,4)+grad(ir,7)*dgrad(ir,5)&
                       +grad(ir,8)*dgrad(ir,6)
                  skab=grad(ir,6)*dgrad(ir,1)+grad(ir,7)*dgrad(ir,2)&
                       +grad(ir,8)*dgrad(ir,3)
                  skba=grad(ir,2)*dgrad(ir,4)+grad(ir,3)*dgrad(ir,5)&
                       +grad(ir,4)*dgrad(ir,6)
                  skbb=grad(ir,2)*dgrad(ir,1)+grad(ir,3)*dgrad(ir,2)&
                       +grad(ir,4)*dgrad(ir,3)
                  tmp1(ir)=(nana+d2fdn2)*drhoe(ir,2)+nanb*drhoe(ir,1)&
                       +naga*2._real_8*skaa+nask*(skba+skab)+nagb*2._real_8*skbb&
                       +d2fdndg*2._real_8*skaa
                  tmp2(ir)=2._real_8*(naga*drhoe(ir,2)+ganb*drhoe(ir,1)+gask*&
                       (skba+skab)+gaga*2._real_8*skaa+gagb*2._real_8*skbb&
                       +d2fdndg*drhoe(IR,2)+2._real_8*skaa*d2fdg2)
                  tmp3(ir)=2._real_8*(ga+dfdg)
                  tmp4(ir)=nask*drhoe(ir,2)+sknb*drhoe(ir,1)+gask*2._real_8*skaa+&
                       skgb*2._real_8*skbb+sksk*(skba+skab)
                  tmp5(ir)=sk
               ENDIF
            ENDDO
         ELSE
            CALL stopgm(procedureN,'Functional not implemented',& 
                 __LINE__,__FILE__)
         ENDIF

         ! ==--------------------------------------------------------------==
         ! CALCULATION OF DIVERGENCE IN RECIPROCAL SPACE:
         ! real space               G-space
         ! dot(nabla,gradient) --> dot(G-Vector,gradient)
         ! tmp2: terms to be multiplied with gradient of a-density
         ! tmp3: terms to be multiplied with gradient of a-linres-density
         ! tmp4: terms to be multiplied with gradient of b-density
         ! tmp5: terms to be multiplied with gradient of b-linres-density
         ! ==--------------------------------------------------------------==
         ! FFTMULT:calculates dot product of G-Vector with gradient 
         DO ir=1,fpar%nnr1
            v(ir,2)=CMPLX(tmp2(ir)*grad(ir,6),tmp2(ir)*grad(ir,7),kind=real_8)
         ENDDO
         CALL  fwfftn(v(:,2),.FALSE.,parai%allgrp)
         CALL fftmult(v(:,2),vtemp,1,2)

         CALL zeroing(v(:,2))!,maxfft)
         DO ir=1,fpar%nnr1
            v(ir,2)=CMPLX(tmp2(ir)*grad(ir,8),tmp3(ir)*dgrad(ir,4),kind=real_8)
         ENDDO
         CALL  fwfftn(v(:,2),.FALSE.,parai%allgrp)
         CALL fftmult(v(:,2),vtemp,3,1)

         CALL zeroing(v(:,2))!,maxfft)
         DO ir=1,fpar%nnr1
            v(ir,2)=CMPLX(tmp3(ir)*dgrad(ir,5),tmp3(ir)*dgrad(ir,6),kind=real_8)
         ENDDO
         CALL  fwfftn(v(:,2),.FALSE.,parai%allgrp)
         CALL fftmult(v(:,2),vtemp,2,3)

         CALL zeroing(v(:,2))!,maxfft)
         DO ir=1,fpar%nnr1
            v(ir,2)=CMPLX(tmp4(ir)*grad(ir,2),tmp4(ir)*grad(ir,3),kind=real_8)
         ENDDO
         CALL  fwfftn(v(:,2),.FALSE.,parai%allgrp)
         CALL fftmult(v(:,2),vtemp,1,2)

         CALL zeroing(v(:,2))!,maxfft)
         DO ir=1,fpar%nnr1
            v(ir,2)=CMPLX(tmp4(ir)*grad(ir,4),tmp5(ir)*dgrad(ir,1),kind=real_8)
         ENDDO
         CALL  fwfftn(v(:,2),.FALSE.,parai%allgrp)
         CALL fftmult(v(:,2),vtemp,3,1)

         CALL zeroing(v(:,2))!,maxfft)
         DO ir=1,fpar%nnr1
            v(ir,2)=CMPLX(tmp5(ir)*dgrad(ir,2),tmp5(ir)*dgrad(ir,3),kind=real_8)
         ENDDO
         CALL  fwfftn(v(:,2),.FALSE.,parai%allgrp)
         CALL fftmult(v(:,2),vtemp,2,3)

         CALL zeroing(v(:,2))!,maxfft)
         DO ig=1,ncpw%nhg
            v(nzh(ig),2) = vtemp(ig)
            v(indz(ig),2) = CONJG(vtemp(ig))
         ENDDO
         CALL  invfftn(v(:,2),.FALSE.,parai%allgrp)

         DO ir=1,fpar%nnr1
            v(ir,2)=v(ir,2)+CMPLX(tmp1(ir),0._real_8,kind=real_8)
         ENDDO
      ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
99  CONTINUE
    IF (cntl%tlsd) THEN
       DO ir=1,fpar%nnr1
          rhoe(ir,1)=rhoe(ir,1)+rhoe(ir,2)
          drhoe(ir,1)=drhoe(ir,1)+drhoe(ir,2)
       ENDDO
    ENDIF
    DO ispin=1,clsd%nlsd
       DO ir=1,fpar%nnr1
          dxc(ir,ispin) = REAL(v(ir,ispin))
       ENDDO
    ENDDO
    IF (cntl%tgc .and. cntl%use_xc_driver) THEN
       DEALLOCATE(dxc_tmp,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    ELSE IF (cntl%tgc .and. .not. cntl%use_xc_driver) THEN
       DEALLOCATE(tmp1,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(tmp2,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(tmp3,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(tmp4,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(tmp5,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(dgrad,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE dd_xc_ana
  ! ==================================================================
  SUBROUTINE graden_comp(rhoe,v,grad,vtmp)
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
    COMPLEX(real_8)                          :: v(maxfft)
    REAL(real_8)                             :: grad(fpar%nnr1,3)
    COMPLEX(real_8)                          :: vtmp(ncpw%nhg)

    INTEGER                                  :: ig, ir, isub
    REAL(real_8)                             :: eg, gcs, gmax, smfac

! Variables
! ==--------------------------------------------------------------==
! ==  TRANSFORM DENSITY TO G SPACE                                ==
! ==--------------------------------------------------------------==

    CALL zeroing(v)!,maxfft)

    CALL tiset('GRADEN_COM',isub)
    DO ir=1,fpar%nnr1
       v(ir) = CMPLX(rhoe(ir),0.0_real_8,kind=real_8)
    ENDDO
    CALL  fwfftn(v,.FALSE.,parai%allgrp)
    ! ==--------------------------------------------------------------==
    ! ==  SMOOTHING                                                   ==
    ! ==--------------------------------------------------------------==
    IF (cntl%tsmooth) THEN
       gmax=hg(ncpw%nhg)
       CALL mp_max(gmax,parai%allgrp)
       gcs=cntr%smf*gmax
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
    DO ig=1,ncpw%nhg
       v(nzh(ig))=vtmp(ig)-parm%tpiba*gk(1,ig)*vtmp(ig)
       v(indz(ig))=CONJG(vtmp(ig)+parm%tpiba*gk(1,ig)*vtmp(ig))
    ENDDO
    CALL  invfftn(v,.FALSE.,parai%allgrp)
    DO ir=1,fpar%nnr1
       rhoe(ir)=REAL(v(ir))
       grad(ir,1)=AIMAG(v(ir))
    ENDDO
    ! ==--------------------------------------------------------------==
    ! ==  FFT OF NABLA(Y)*RHO AND NABLA(Z)*RHOE                       ==
    ! ==--------------------------------------------------------------==
    CALL zeroing(v)!,maxfft)
    DO ig=1,ncpw%nhg
       v(nzh(ig))=parm%tpiba*(uimag*gk(2,ig)-gk(3,ig))*vtmp(ig)
       v(indz(ig))=parm%tpiba*(-uimag*gk(2,ig)+gk(3,ig))*CONJG(vtmp(ig))
    ENDDO
    CALL  invfftn(v,.FALSE.,parai%allgrp)
    DO ir=1,fpar%nnr1
       grad(ir,2)=REAL(v(ir))
       grad(ir,3)=AIMAG(v(ir))
    ENDDO
    CALL tihalt('GRADEN_COM',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE graden_comp
  ! ==================================================================
  SUBROUTINE give_scr_dd_xc_ana(ldd_xc_ana,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ldd_xc_ana
    CHARACTER(len=*)                         :: tag

! ==--------------------------------------------------------------==

    IF (.NOT.(cntl%tlsd.OR.lspin2%tlse).AND.cntl%tgc) THEN
       ldd_xc_ana=5*fpar%nnr1+3*fpar%nnr1-1*fpar%nnr1
    ELSEIF (cntl%tlsd.AND..NOT.cntl%tgc) THEN
       ldd_xc_ana=1
    ELSEIF ((cntl%tlsd.OR.lspin2%tlse).AND.cntl%tgc) THEN
       ldd_xc_ana=5*fpar%nnr1+6*fpar%nnr1-2*fpar%nnr1
    ELSE
       ldd_xc_ana=1
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_dd_xc_ana
  ! =================================================================
  SUBROUTINE fftmult(v,vtemp,coord1,coord2)
    COMPLEX(real_8)                          :: v(maxfft), vtemp(ncpw%nhg)
    INTEGER                                  :: coord1, coord2

    COMPLEX(real_8)                          :: fa, fb, vg1, vg2
    INTEGER                                  :: ig
    REAL(real_8)                             :: vfac1, vfac2

! ==--------------------------------------------------------------==

    DO ig=1,ncpw%nhg
       vfac1=parm%tpiba*gk(coord1,ig)
       vfac2=parm%tpiba*gk(coord2,ig)
       fa  = v(nzh(ig)) + v(indz(ig))
       fb  = v(nzh(ig)) - v(indz(ig))
       vg1 = 0.5_real_8*CMPLX(REAL(fa),AIMAG(fb),kind=real_8)
       vg2 = 0.5_real_8*CMPLX(AIMAG(fa),-REAL(fb),kind=real_8)
       vtemp(ig) = vtemp(ig) - vfac1*uimag*vg1 - vfac2*uimag*vg2
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE fftmult
  ! =================================================================     

END MODULE dd_xc_ana_utils
