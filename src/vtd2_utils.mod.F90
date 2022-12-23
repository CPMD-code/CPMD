MODULE vtd2_utils
  USE dd_xc_utils,                     ONLY: vxc_calc
  USE density_functionals_utils,       ONLY: pade_lda
  USE error_handling,                  ONLY: stopgm
  USE func,                            ONLY: func1,&
                                             mfxcc_is_pade,&
                                             mfxcc_is_skipped,&
                                             mfxcx_is_skipped
  USE kinds,                           ONLY: real_8
  USE linres,                          ONLY: lr02,&
                                             lr03,&
                                             lrsym,&
                                             td01,&
                                             tshl
  USE poin,                            ONLY: rhoo
  USE spin,                            ONLY: clsd
  USE switch_functionals_utils,        ONLY: switch_functionals
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             ncpw
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: vtd2
  PUBLIC :: vtd2t
  PUBLIC :: give_scr_vtd2

CONTAINS

  ! ==================================================================
  SUBROUTINE vtd2(vp2,rhod,psi,swfun)
    ! ==--------------------------------------------------------------==
    ! 
    REAL(real_8)                             :: vp2(:,:), rhod(:,:)
    COMPLEX(real_8)                          :: psi(:,:)
    LOGICAL                                  :: swfun

    CHARACTER(*), PARAMETER                  :: procedureN = 'vtd2'
    REAL(real_8), PARAMETER                  :: small = 1.e-14_real_8 

    COMPLEX(real_8), ALLOCATABLE             :: vtemp(:,:), vtmp(:)
    INTEGER                                  :: ierr, ir, isub, ndim
    REAL(real_8)                             :: oeps
    REAL(real_8), ALLOCATABLE                :: grad(:,:), rhox(:,:)

! ==--------------------------------------------------------------==

    CALL tiset('      VTD2',isub)
    ! ==--------------------------------------------------------------==
    ! 
    IF (swfun) CALL switch_functionals
    ! 
    IF (.NOT.lr03%txc_analytic) THEN
       ALLOCATE(rhox(fpar%nnr1, clsd%nlsd),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       IF (cntl%tgc) THEN
          ndim=clsd%nlsd
          IF(td01%ns_tri.GT.0.OR.tshl%isc) ndim=2

          ALLOCATE(vtemp(ncpw%nhg,ndim),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
          ALLOCATE(vtmp(MAX(ncpw%nhg, fpar%nnr1*ndim)),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
          ALLOCATE(grad(fpar%nnr1, ndim * 4),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
       ENDIF
    ENDIF
    ! 
    IF (cntl%tlsd) THEN
       IF (lr03%txc_analytic) THEN
          CALL stopgm('VTD2','ANALYTIC XC NOT AVAILABLE WITH LSD',& 
               __LINE__,__FILE__)
       ELSE
          IF (func1%mfxcc == mfxcc_is_skipped.AND.func1%mfxcx == mfxcx_is_skipped.AND..NOT.cntl%tgc) THEN
             CALL zeroing(vp2)!,2*nnr1)
          ELSE
             IF (.NOT.ALLOCATED(vtemp)) CALL stopgm(procedureN,'vtemp not allocated ',&
                  __LINE__,__FILE__)
             ! 5-point formula for finite differences
             ! R0
             CALL vxc_calc(rhoo,psi,vtemp,vtmp,grad,is_alpha_beta_dens=.TRUE.)
             DO ir=1,fpar%nnr1
                vp2(ir,1)=-30._real_8*REAL(psi(ir,1))
                vp2(ir,2)=-30._real_8*REAL(psi(ir,2))
             ENDDO
             ! R0+2*RD
             CALL dcopy(2*fpar%nnr1,rhoo,1,rhox,1)
             CALL daxpy(2*fpar%nnr1,2._real_8*lr02%xc_eps,rhod,1,rhox,1)
             CALL vxc_calc(rhox,psi,vtemp,vtmp,grad,is_alpha_beta_dens=.TRUE.)
             DO ir=1,fpar%nnr1
                vp2(ir,1)=vp2(ir,1)-REAL(psi(ir,1))
                vp2(ir,2)=vp2(ir,2)-REAL(psi(ir,2))
             ENDDO
             ! R0+RD
             CALL dcopy(2*fpar%nnr1,rhoo,1,rhox,1)
             CALL daxpy(2*fpar%nnr1,lr02%xc_eps,rhod,1,rhox,1)
             CALL vxc_calc(rhox,psi,vtemp,vtmp,grad,is_alpha_beta_dens=.TRUE.)
             DO ir=1,fpar%nnr1
                vp2(ir,1)=vp2(ir,1)+16._real_8*REAL(psi(ir,1))
                vp2(ir,2)=vp2(ir,2)+16._real_8*REAL(psi(ir,2))
             ENDDO
             ! R0-RD
             CALL dcopy(2*fpar%nnr1,rhoo,1,rhox,1)
             CALL daxpy(2*fpar%nnr1,-lr02%xc_eps,rhod,1,rhox,1)
             CALL vxc_calc(rhox,psi,vtemp,vtmp,grad,is_alpha_beta_dens=.TRUE.)
             DO ir=1,fpar%nnr1
                vp2(ir,1)=vp2(ir,1)+16._real_8*REAL(psi(ir,1))
                vp2(ir,2)=vp2(ir,2)+16._real_8*REAL(psi(ir,2))
             ENDDO
             ! R0-2*RD
             CALL dcopy(2*fpar%nnr1,rhoo,1,rhox,1)
             CALL daxpy(2*fpar%nnr1,-2._real_8*lr02%xc_eps,rhod,1,rhox,1)
             CALL vxc_calc(rhox,psi,vtemp,vtmp,grad,is_alpha_beta_dens=.TRUE.)
             DO ir=1,fpar%nnr1
                vp2(ir,1)=vp2(ir,1)-REAL(psi(ir,1))
                vp2(ir,2)=vp2(ir,2)-REAL(psi(ir,2))
             ENDDO
             oeps=1._real_8/(12._real_8*lr02%xc_eps*lr02%xc_eps)
             CALL dscal(2*fpar%nnr1,oeps,vp2,1)
             CALL dscal(2*fpar%nnr1,0.25_real_8,vp2,1)
          ENDIF
       ENDIF
    ELSE
       IF (lr03%txc_analytic) THEN
          IF (func1%mfxcc == mfxcc_is_skipped.AND.func1%mfxcx == mfxcx_is_skipped.AND..NOT.cntl%tgc) THEN
             CALL stopgm('VTD2','this functional not available',& 
                  __LINE__,__FILE__)
          ELSEIF (func1%mfxcc == mfxcc_is_pade.AND.func1%mfxcx == mfxcx_is_skipped.AND..NOT.cntl%tgc) THEN
             DO ir=1,fpar%nnr1
                IF (rhoo(ir,1).GE.small)THEN
                   vp2(ir,1)=pade_lda(rhoo(ir,1),3)*rhod(ir,1)*rhod(ir,1)
                ELSE
                   vp2(ir,1)=0._real_8
                ENDIF
             ENDDO
          ELSE
             CALL stopgm('VTD2','this functional not available',& 
                  __LINE__,__FILE__)
          ENDIF
          IF (cntl%thybrid) THEN
             CALL stopgm('VTD2','this functional not available',& 
                  __LINE__,__FILE__)
          ENDIF
       ELSE
          IF (func1%mfxcc == mfxcc_is_skipped.AND.func1%mfxcx == mfxcx_is_skipped.AND..NOT.cntl%tgc) THEN
             CALL zeroing(vp2)!,nnr1)
          ELSE
             IF (.NOT.ALLOCATED(vtemp)) CALL stopgm(procedureN,'vtemp not allocated ',&
                  __LINE__,__FILE__)
             ! 5-point formula for finite differences
             ! R0
             CALL vxc_calc(rhoo,psi,vtemp,vtmp,grad)
             DO ir=1,fpar%nnr1
                vp2(ir,1)=-30._real_8*REAL(psi(ir,1))
             ENDDO
             ! R0+2*RD
             CALL dcopy(fpar%nnr1,rhoo,1,rhox,1)
             CALL daxpy(fpar%nnr1,2._real_8*lr02%xc_eps,rhod,1,rhox,1)
             CALL vxc_calc(rhox,psi,vtemp,vtmp,grad)
             DO ir=1,fpar%nnr1
                vp2(ir,1)=vp2(ir,1)-REAL(psi(ir,1))
             ENDDO
             ! R0+RD
             CALL dcopy(fpar%nnr1,rhoo,1,rhox,1)
             CALL daxpy(fpar%nnr1,lr02%xc_eps,rhod,1,rhox,1)
             CALL vxc_calc(rhox,psi,vtemp,vtmp,grad)
             DO ir=1,fpar%nnr1
                vp2(ir,1)=vp2(ir,1)+16._real_8*REAL(psi(ir,1))
             ENDDO
             ! R0-RD
             CALL dcopy(fpar%nnr1,rhoo,1,rhox,1)
             CALL daxpy(fpar%nnr1,-lr02%xc_eps,rhod,1,rhox,1)
             CALL vxc_calc(rhox,psi,vtemp,vtmp,grad)
             DO ir=1,fpar%nnr1
                vp2(ir,1)=vp2(ir,1)+16._real_8*REAL(psi(ir,1))
             ENDDO
             ! R0-2*RD
             CALL dcopy(fpar%nnr1,rhoo,1,rhox,1)
             CALL daxpy(fpar%nnr1,-2._real_8*lr02%xc_eps,rhod,1,rhox,1)
             CALL vxc_calc(rhox,psi,vtemp,vtmp,grad)
             DO ir=1,fpar%nnr1
                vp2(ir,1)=vp2(ir,1)-REAL(psi(ir,1))
             ENDDO
             oeps=1._real_8/(12._real_8*lr02%xc_eps*lr02%xc_eps)
             CALL dscal(fpar%nnr1,oeps,vp2,1)
          ENDIF
       ENDIF
    ENDIF
    ! 
    IF (.NOT. lr03%txc_analytic) THEN
       DEALLOCATE(rhox,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       IF (cntl%tgc) THEN
          DEALLOCATE(vtmp,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
               __LINE__,__FILE__)
          DEALLOCATE(vtemp,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
               __LINE__,__FILE__)
          DEALLOCATE(grad,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
               __LINE__,__FILE__)
       ENDIF
    ENDIF
    ! 
    IF (swfun) CALL switch_functionals
    ! ==--------------------------------------------------------------==
    CALL tihalt('      VTD2',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE vtd2
  ! ==================================================================
  SUBROUTINE vtd2t(vp2,rhod,psi,swfun)
    ! ==--------------------------------------------------------------==
    ! 
    REAL(real_8)                             :: vp2(fpar%nnr1), &
                                                rhod(fpar%nnr1)
    COMPLEX(real_8)                          :: psi(:,:)
    LOGICAL                                  :: swfun

    CHARACTER(*), PARAMETER                  :: procedureN = 'vtd2t'
    REAL(real_8), PARAMETER                  :: small = 1.e-14_real_8 

    COMPLEX(real_8), ALLOCATABLE             :: vtemp(:,:), vtmp(:)
    INTEGER                                  :: ierr, ir, isub, ndim
    REAL(real_8)                             :: oeps
    REAL(real_8), ALLOCATABLE                :: grad(:,:), rhox(:,:)

    CALL tiset('      VTD2',isub)
    ! ==--------------------------------------------------------------==
    ! 
    IF (swfun) CALL switch_functionals
    ! 
    cntl%tlsd=.TRUE.
    ! 
    IF (.NOT. lr03%txc_analytic) THEN
       ALLOCATE(rhox(fpar%nnr1, 2),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       IF (cntl%tgc) THEN
          ndim=clsd%nlsd
          IF(td01%ns_tri.GT.0.OR.tshl%isc) ndim=2

          ALLOCATE(vtmp(MAX(ncpw%nhg, fpar%nnr1*ndim)),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
          ALLOCATE(vtemp(ncpw%nhg,ndim),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
          ALLOCATE(grad(fpar%nnr1,ndim*4),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
       ENDIF
    ENDIF

    ! 
    IF (lr03%txc_analytic) THEN
       CALL stopgm('VTD2','this functional not available',& 
            __LINE__,__FILE__)
    ELSE
       IF (func1%mfxcc == mfxcc_is_skipped.AND.func1%mfxcx == mfxcx_is_skipped.AND..NOT.cntl%tgc) THEN
          CALL zeroing(vp2)!,nnr1)
       ELSE
          ! 5-point formula for finite differences
          ! R0
          CALL dcopy(fpar%nnr1,rhoo,1,rhox(1,1),1)
          CALL dcopy(fpar%nnr1,rhoo,1,rhox(1,2),1)
          CALL dscal(fpar%nnr1,0.5_real_8,rhox(1,2),1)
          CALL vxc_calc(rhox,psi,vtemp,vtmp,grad)
          DO ir=1,fpar%nnr1
             vp2(ir)=-30._real_8*REAL(psi(ir,1))
          ENDDO
          ! R0+2*RD
          CALL dcopy(fpar%nnr1,rhoo,1,rhox(1,1),1)
          CALL dcopy(fpar%nnr1,rhoo,1,rhox(1,2),1)
          CALL daxpy(fpar%nnr1,2._real_8*lr02%xc_eps,rhod,1,rhox(1,2),1)
          CALL dscal(fpar%nnr1,0.5_real_8,rhox(1,2),1)
          CALL vxc_calc(rhox,psi,vtemp,vtmp,grad)
          DO ir=1,fpar%nnr1
             vp2(ir)=vp2(ir)-REAL(psi(ir,1))
          ENDDO
          ! R0+RD
          CALL dcopy(fpar%nnr1,rhoo,1,rhox(1,1),1)
          CALL dcopy(fpar%nnr1,rhoo,1,rhox(1,2),1)
          CALL daxpy(fpar%nnr1,lr02%xc_eps,rhod,1,rhox(1,2),1)
          CALL dscal(fpar%nnr1,0.5_real_8,rhox(1,2),1)
          CALL vxc_calc(rhox,psi,vtemp,vtmp,grad)
          DO ir=1,fpar%nnr1
             vp2(ir)=vp2(ir)+16._real_8*REAL(psi(ir,1))
          ENDDO
          ! R0-RD
          CALL dcopy(fpar%nnr1,rhoo,1,rhox(1,1),1)
          CALL dcopy(fpar%nnr1,rhoo,1,rhox(1,2),1)
          CALL daxpy(fpar%nnr1,-lr02%xc_eps,rhod,1,rhox(1,2),1)
          CALL dscal(fpar%nnr1,0.5_real_8,rhox(1,2),1)
          CALL vxc_calc(rhox,psi,vtemp,vtmp,grad)
          DO ir=1,fpar%nnr1
             vp2(ir)=vp2(ir)+16._real_8*REAL(psi(ir,1))
          ENDDO
          ! R0-2*RD
          CALL dcopy(fpar%nnr1,rhoo,1,rhox(1,1),1)
          CALL dcopy(fpar%nnr1,rhoo,1,rhox(1,2),1)
          CALL daxpy(fpar%nnr1,-2._real_8*lr02%xc_eps,rhod,1,rhox(1,2),1)
          CALL dscal(fpar%nnr1,0.5_real_8,rhox(1,2),1)
          CALL vxc_calc(rhox,psi,vtemp,vtmp,grad)
          DO ir=1,fpar%nnr1
             vp2(ir)=vp2(ir)-REAL(psi(ir,1))
          ENDDO
          oeps=1._real_8/(12._real_8*lr02%xc_eps*lr02%xc_eps)
          CALL dscal(fpar%nnr1,oeps,vp2,1)
       ENDIF
    ENDIF
    ! 
    IF (.NOT. lr03%txc_analytic) THEN
       DEALLOCATE(rhox,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       IF (cntl%tgc) THEN
          DEALLOCATE(vtmp,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
               __LINE__,__FILE__)
          DEALLOCATE(vtemp,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
               __LINE__,__FILE__)
          DEALLOCATE(grad,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
               __LINE__,__FILE__)
       ENDIF
    ENDIF
    ! 
    IF (swfun) CALL switch_functionals
    ! 
    cntl%tlsd=.FALSE.
    ! ==--------------------------------------------------------------==
    CALL tihalt('      VTD2',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE vtd2t
  ! ==================================================================
  SUBROUTINE give_scr_vtd2(lvtd2,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lvtd2
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: lgrad, lrhox, lvtemp, lvtmp, &
                                                nm

! ==--------------------------------------------------------------==

    IF (.NOT.lr03%txc_analytic) THEN
       IF (lrsym.EQ.3) THEN
          nm=2
       ELSE
          nm=clsd%nlsd
       ENDIF
       lrhox=fpar%nnr1*nm
       IF (cntl%tgc) THEN
          lgrad=fpar%nnr1*nm*4
          lvtemp=2*ncpw%nhg
          lvtmp=MAX(2*ncpw%nhg,fpar%nnr1*(2*nm-1))
       ELSE
          lgrad=1
          lvtemp=1
          lvtmp=1
       ENDIF
       lvtd2=lrhox+lgrad+lvtemp+lvtmp
    ELSE
       lvtd2=0
    ENDIF
    lvtd2=lvtd2+200
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_vtd2
  ! ==================================================================

END MODULE vtd2_utils
