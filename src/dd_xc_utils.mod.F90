MODULE dd_xc_utils
  USE cppt,                            ONLY: nzh
  USE dd_xc_ana_utils,                 ONLY: dd_xc_ana
  USE error_handling,                  ONLY: stopgm
  USE fftmain_utils,                   ONLY: fwfftn
  USE gcener_utils,                    ONLY: gcener,&
                                             gclsd
  USE graden_utils,                    ONLY: graden
  USE kinds,                           ONLY: real_8
  USE linres,                          ONLY: lr01,&
                                             lr02,&
                                             lr03
  USE nlcc,                            ONLY: corel,&
                                             roct
  USE parac,                           ONLY: parai
  USE spin,                            ONLY: clsd
  USE switch_functionals_utils,        ONLY: switch_functionals
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             ncpw
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: zgthr
  USE xcener_utils,                    ONLY: xcener
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: dd_xc
  PUBLIC :: dd_xct
  PUBLIC :: vxc_calc

CONTAINS

  ! ==================================================================
  SUBROUTINE dd_xc(rhoe,drhoe,rhox,psi,dxc,vtmp,vtemp,grad,switch)
    ! ==--------------------------------------------------------------==
    ! == Calculate dVxc/dn*n1 by finite differences                   ==
    ! ==--------------------------------------------------------------==
    ! == INPUT: rhoe = density                                        ==
    ! ==        drhoe = LR density                                    ==
    ! ==        rhox = undef.                                         ==
    ! ==        psi = undef.                                          ==
    ! ==        dxc = undef.                                          ==
    ! ==        vtmp = undef.                                         ==
    ! ==        vtemp = undef.                                        ==
    ! ==        grad = undef.                                         ==
    ! == OUTPUT: rhoe = density                                       ==
    ! ==        drhoe = LR density                                    ==
    ! ==        rhox = undef.                                         ==
    ! ==        psi = undef.                                          ==
    ! ==        dxc = dVxc/dn*n1                                      ==
    ! ==        vtmp = potential in G space                           ==
    ! ==        vtemp = undef.                                        ==
    ! ==        grad = undef.                                         ==
    ! ==--------------------------------------------------------------==
    CHARACTER(len=*), PARAMETER              :: procedureN = 'dd_xc'
    REAL(real_8)                             :: rhoe(:,:), drhoe(:,:), &
                                                rhox(:,:)
    COMPLEX(real_8)                          :: psi(:,:)
    REAL(real_8)                             :: dxc(:,:)
    COMPLEX(real_8)                          :: vtmp(:), vtemp(:,:)
    REAL(real_8)                             :: grad(:,:)
    LOGICAL                                  :: switch

    INTEGER                                  :: ir, isub
    REAL(real_8)                             :: oeps

! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)

    IF (lr03%txc_dd_ana) THEN
       CALL dd_xc_ana(rhoe,drhoe,rhox,psi,dxc,vtmp,vtemp(:,1),grad,switch)
       CALL tihalt(procedureN,isub)
       RETURN
    ENDIF
    ! ==--------------------------------------------------------------==
    ! 
    IF (switch) CALL switch_functionals
    ! 
    IF (cntl%ttau) THEN
       CALL stopgm('DD_XC','META FUNCTIONALS NOT IMPLENTED',& 
            __LINE__,__FILE__)
    ENDIF
    ! 
    IF (.NOT. cntl%tlsd) THEN
       ! ..NON-spin polarised
       IF (lr01%ndpoint.EQ.2) THEN
          ! ..2-point formula
          CALL dcopy(fpar%nnr1,rhoe,1,rhox,1)
          CALL daxpy(fpar%nnr1,lr02%xc_eps,drhoe,1,rhox,1)
          CALL vxc_calc(rhox,psi,vtemp,vtmp,grad)
          DO ir=1,fpar%nnr1
             dxc(ir,1)=REAL(psi(ir,1))
          ENDDO
          CALL dcopy(fpar%nnr1,rhoe,1,rhox,1)
          CALL daxpy(fpar%nnr1,-lr02%xc_eps,drhoe,1,rhox,1)
          CALL vxc_calc(rhox,psi,vtemp,vtmp,grad)
          DO ir=1,fpar%nnr1
             dxc(ir,1)=dxc(ir,1)-REAL(psi(ir,1))
          ENDDO
          oeps=0.5_real_8/lr02%xc_eps
          CALL dscal(fpar%nnr1,oeps,dxc,1)
       ELSEIF (lr01%ndpoint.EQ.4) THEN
          ! ..4-point formula
          CALL dcopy(fpar%nnr1,rhoe,1,rhox,1)
          CALL daxpy(fpar%nnr1,2._real_8*lr02%xc_eps,drhoe,1,rhox,1)
          CALL vxc_calc(rhox,psi,vtemp,vtmp,grad)
          DO ir=1,fpar%nnr1
             dxc(ir,1)=REAL(psi(ir,1))
          ENDDO
          CALL dcopy(fpar%nnr1,rhoe,1,rhox,1)
          CALL daxpy(fpar%nnr1,lr02%xc_eps,drhoe,1,rhox,1)
          CALL vxc_calc(rhox,psi,vtemp,vtmp,grad)
          DO ir=1,fpar%nnr1
             dxc(ir,1)=dxc(ir,1)+8._real_8*REAL(psi(ir,1))
          ENDDO
          CALL dcopy(fpar%nnr1,rhoe,1,rhox,1)
          CALL daxpy(fpar%nnr1,-lr02%xc_eps,drhoe,1,rhox,1)
          CALL vxc_calc(rhox,psi,vtemp,vtmp,grad)
          DO ir=1,fpar%nnr1
             dxc(ir,1)=dxc(ir,1)-8._real_8*REAL(psi(ir,1))
          ENDDO
          CALL dcopy(fpar%nnr1,rhoe,1,rhox,1)
          CALL daxpy(fpar%nnr1,-2._real_8*lr02%xc_eps,drhoe,1,rhox,1)
          CALL vxc_calc(rhox,psi,vtemp,vtmp,grad)
          DO ir=1,fpar%nnr1
             dxc(ir,1)=dxc(ir,1)-REAL(psi(ir,1))
          ENDDO
          oeps=1._real_8/(12._real_8*lr02%xc_eps)
          CALL dscal(fpar%nnr1,oeps,dxc,1)
       ELSE
          CALL stopgm('DD_XC','N Point formula not supported',& 
               __LINE__,__FILE__)
       ENDIF
    ELSE
       ! ..spin polarised
       IF (lr01%ndpoint.EQ.2) THEN
          ! ..2-point formula
          CALL dcopy(2*fpar%nnr1,rhoe,1,rhox,1)
          CALL daxpy(2*fpar%nnr1,lr02%xc_eps,drhoe,1,rhox,1)
          CALL vxc_calc(rhox,psi,vtemp,vtmp,grad)
          DO ir=1,fpar%nnr1
             dxc(ir,1)=REAL(psi(ir,1))
             dxc(ir,2)=REAL(psi(ir,2))
          ENDDO
          CALL dcopy(2*fpar%nnr1,rhoe,1,rhox,1)
          CALL daxpy(2*fpar%nnr1,-lr02%xc_eps,drhoe,1,rhox,1)
          CALL vxc_calc(rhox,psi,vtemp,vtmp,grad)
          DO ir=1,fpar%nnr1
             dxc(ir,1)=dxc(ir,1)-REAL(psi(ir,1))
             dxc(ir,2)=dxc(ir,2)-REAL(psi(ir,2))
          ENDDO
          oeps=0.5_real_8/lr02%xc_eps
          CALL dscal(fpar%nnr1,oeps,dxc(1,1),1)
          CALL dscal(fpar%nnr1,oeps,dxc(1,2),1)
       ELSEIF (lr01%ndpoint.EQ.4) THEN
          ! ..4-point formula
          CALL dcopy(2*fpar%nnr1,rhoe,1,rhox,1)
          CALL daxpy(2*fpar%nnr1,2._real_8*lr02%xc_eps,drhoe,1,rhox,1)
          CALL vxc_calc(rhox,psi,vtemp,vtmp,grad)
          DO ir=1,fpar%nnr1
             dxc(ir,1)=REAL(psi(ir,1))
             dxc(ir,2)=REAL(psi(ir,2))
          ENDDO
          CALL dcopy(2*fpar%nnr1,rhoe,1,rhox,1)
          CALL daxpy(2*fpar%nnr1,lr02%xc_eps,drhoe,1,rhox,1)
          CALL vxc_calc(rhox,psi,vtemp,vtmp,grad)
          DO ir=1,fpar%nnr1
             dxc(ir,1)=dxc(ir,1)+8._real_8*REAL(psi(ir,1))
             dxc(ir,2)=dxc(ir,2)+8._real_8*REAL(psi(ir,2))
          ENDDO
          CALL dcopy(2*fpar%nnr1,rhoe,1,rhox,1)
          CALL daxpy(2*fpar%nnr1,-lr02%xc_eps,drhoe,1,rhox,1)
          CALL vxc_calc(rhox,psi,vtemp,vtmp,grad)
          DO ir=1,fpar%nnr1
             dxc(ir,1)=dxc(ir,1)-8._real_8*REAL(psi(ir,1))
             dxc(ir,2)=dxc(ir,2)-8._real_8*REAL(psi(ir,2))
          ENDDO
          CALL dcopy(2*fpar%nnr1,rhoe,1,rhox,1)
          CALL daxpy(2*fpar%nnr1,-2._real_8*lr02%xc_eps,drhoe,1,rhox,1)
          CALL vxc_calc(rhox,psi,vtemp,vtmp,grad)
          DO ir=1,fpar%nnr1
             dxc(ir,1)=dxc(ir,1)-REAL(psi(ir,1))
             dxc(ir,2)=dxc(ir,2)-REAL(psi(ir,2))
          ENDDO
          oeps=1._real_8/(12._real_8*lr02%xc_eps)
          CALL dscal(fpar%nnr1,oeps,dxc(1,1),1)
          CALL dscal(fpar%nnr1,oeps,dxc(1,2),1)
       ELSE
          CALL stopgm(procedureN,'N Point formula not supported',& 
               __LINE__,__FILE__)
       ENDIF
    ENDIF
    ! 
    IF (switch) CALL switch_functionals
    ! 
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE dd_xc
  ! ==================================================================
  SUBROUTINE dd_xct(rhoe,drhoe,rhox,psi,dxc,vtmp,vtemp,grad,switch)
    ! ==--------------------------------------------------------------==
    ! == Calculate dVxc/dn*n1 by finite differences                   ==
    ! == cntl%tddft triplet case                                           ==
    ! ==--------------------------------------------------------------==
    ! == INPUT: rhoe = density                                        ==
    ! ==        drhoe = LR density                                    ==
    ! ==        rhox = undef.                                         ==
    ! ==        psi = undef.                                          ==
    ! ==        dxc = undef.                                          ==
    ! ==        vtmp = undef.                                         ==
    ! ==        vtemp = undef.                                        ==
    ! ==        grad = undef.                                         ==
    ! == OUTPUT: rhoe = density                                       ==
    ! ==        drhoe = LR density                                    ==
    ! ==        rhox = undef.                                         ==
    ! ==        psi = undef.                                          ==
    ! ==        dxc = dVxc/dn*n1                                      ==
    ! ==        vtmp = potential in G space                           ==
    ! ==        vtemp = undef.                                        ==
    ! ==        grad = undef.                                         ==
    ! ==--------------------------------------------------------------==
    CHARACTER(len=*), PARAMETER              :: procedureN = 'dd_xct'
    REAL(real_8)                             :: rhoe(:,:), drhoe(:,:), &
                                                rhox(:,:)
    COMPLEX(real_8)                          :: psi(:,:)
    REAL(real_8)                             :: dxc(:,:)
    COMPLEX(real_8)                          :: vtmp(:), vtemp(:,:)
    REAL(real_8)                             :: grad(:,:)
    LOGICAL                                  :: switch

    INTEGER                                  :: ir, isub, nlsdsv
    LOGICAL                                  :: tlsdsv
    REAL(real_8)                             :: oeps

! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    tlsdsv=cntl%tlsd
    nlsdsv=clsd%nlsd
    cntl%tlsd=.TRUE.
    clsd%nlsd=2
    IF (lr03%txc_dd_ana) THEN
       CALL dcopy(fpar%nnr1,rhoe(1,1),1,rhoe(1,2),1)
       CALL dscal(fpar%nnr1,0.5_real_8,rhoe(1,2),1)
       CALL dcopy(fpar%nnr1,drhoe(1,1),1,drhoe(1,2),1)
       CALL zeroing(drhoe(:,1))!,nnr1)
       CALL dd_xc_ana(rhoe,drhoe,rhox,psi,dxc,vtmp,vtemp(:,1),grad,&
            switch)
       CALL dscal(fpar%nnr1,-0.5_real_8,dxc(1,1),1)
    ELSE
       ! 
       IF (switch) CALL switch_functionals
       ! 
       IF (cntl%ttau) THEN
          CALL stopgm('DD_XC','META FUNCTIONALS NOT IMPLENTED',& 
               __LINE__,__FILE__)
       ENDIF
       ! 
       IF (lr01%ndpoint.EQ.2) THEN
          ! ..2-point formula
          CALL dcopy(fpar%nnr1,rhoe,1,rhox(1,1),1)
          CALL dcopy(fpar%nnr1,rhoe,1,rhox(1,2),1)
          CALL daxpy(fpar%nnr1,2._real_8*lr02%xc_eps,drhoe,1,rhox(1,2),1)
          CALL dscal(fpar%nnr1,0.5_real_8,rhox(1,2),1)
          CALL vxc_calc(rhox,psi,vtemp,vtmp,grad)
          DO ir=1,fpar%nnr1
             dxc(ir,1)=REAL(psi(ir,1))
          ENDDO
          CALL dcopy(fpar%nnr1,rhoe,1,rhox(1,1),1)
          CALL dcopy(fpar%nnr1,rhoe,1,rhox(1,2),1)
          CALL daxpy(fpar%nnr1,-2._real_8*lr02%xc_eps,drhoe,1,rhox(1,2),1)
          CALL dscal(fpar%nnr1,0.5_real_8,rhox(1,2),1)
          CALL vxc_calc(rhox,psi,vtemp,vtmp,grad)
          DO ir=1,fpar%nnr1
             dxc(ir,1)=dxc(ir,1)-REAL(psi(ir,1))
          ENDDO
          oeps=0.5_real_8/lr02%xc_eps
          CALL dscal(fpar%nnr1,-0.5_real_8*oeps,dxc(1,1),1)
       ELSEIF (lr01%ndpoint.EQ.4) THEN
          ! ..4-point formula
          CALL dcopy(fpar%nnr1,rhoe,1,rhox(1,1),1)
          CALL dcopy(fpar%nnr1,rhoe,1,rhox(1,2),1)
          CALL daxpy(fpar%nnr1,4._real_8*lr02%xc_eps,drhoe,1,rhox(1,2),1)
          CALL dscal(fpar%nnr1,0.5_real_8,rhox(1,2),1)
          CALL vxc_calc(rhox,psi,vtemp,vtmp,grad)
          DO ir=1,fpar%nnr1
             dxc(ir,1)=REAL(psi(ir,1))
          ENDDO
          CALL dcopy(fpar%nnr1,rhoe,1,rhox(1,1),1)
          CALL dcopy(fpar%nnr1,rhoe,1,rhox(1,2),1)
          CALL daxpy(fpar%nnr1,2._real_8*lr02%xc_eps,drhoe,1,rhox(1,2),1)
          CALL dscal(fpar%nnr1,0.5_real_8,rhox(1,2),1)
          CALL vxc_calc(rhox,psi,vtemp,vtmp,grad)
          DO ir=1,fpar%nnr1
             dxc(ir,1)=dxc(ir,1)+8._real_8*REAL(psi(ir,1))
          ENDDO
          CALL dcopy(fpar%nnr1,rhoe,1,rhox(1,1),1)
          CALL dcopy(fpar%nnr1,rhoe,1,rhox(1,2),1)
          CALL daxpy(fpar%nnr1,-2._real_8*lr02%xc_eps,drhoe,1,rhox(1,2),1)
          CALL dscal(fpar%nnr1,0.5_real_8,rhox(1,2),1)
          CALL vxc_calc(rhox,psi,vtemp,vtmp,grad)
          DO ir=1,fpar%nnr1
             dxc(ir,1)=dxc(ir,1)-8._real_8*REAL(psi(ir,1))
          ENDDO
          CALL dcopy(fpar%nnr1,rhoe,1,rhox(1,1),1)
          CALL dcopy(fpar%nnr1,rhoe,1,rhox(1,2),1)
          CALL daxpy(fpar%nnr1,-4._real_8*lr02%xc_eps,drhoe,1,rhox(1,2),1)
          CALL dscal(fpar%nnr1,0.5_real_8,rhox(1,2),1)
          CALL vxc_calc(rhox,psi,vtemp,vtmp,grad)
          DO ir=1,fpar%nnr1
             dxc(ir,1)=dxc(ir,1)-REAL(psi(ir,1))
          ENDDO
          oeps=1._real_8/(12._real_8*lr02%xc_eps)
          CALL dscal(fpar%nnr1,-0.5_real_8*oeps,dxc(1,1),1)
       ELSE
          CALL stopgm('DD_XC','N Point formula not supported',& 
               __LINE__,__FILE__)
       ENDIF
       ! 
       IF (switch) CALL switch_functionals
       ! 
    ENDIF
    cntl%tlsd=tlsdsv
    clsd%nlsd=nlsdsv
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE dd_xct
  ! ==================================================================
  SUBROUTINE vxc_calc(rhoe,v,vtemp,vtmp,grad,is_alpha_beta_dens)
    ! ==--------------------------------------------------------------==
    ! == COMPUTES XC POTENTIAL                                        ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rhoe(:,:)
    COMPLEX(real_8)                          :: v(:,:), vtemp(:,:), vtmp(:)
    REAL(real_8)                             :: grad(:,:)
    LOGICAL, INTENT(in), OPTIONAL            :: is_alpha_beta_dens

    CHARACTER(*), PARAMETER                  :: procedureN = 'vxc_calc'

    INTEGER                                  :: ierr, ir, isub
    LOGICAL                                  :: my_is_alpha_beta_dens
    REAL(real_8)                             :: sgcc, sgcx, sxc, vxc
    REAL(real_8), ALLOCATABLE                :: rvtmp(:,:)

! Variables
! ==--------------------------------------------------------------==
! == ADD CORE CHARGE TO RHOE                                      ==
! == about  my_is_alpha_beta_dens                                 ==
! ==  my_is_alpha_beta_dens is used as optional parameter to      ==
! ==  skip the subtraction of rhoe(:,2) from rhoe(:,1)            ==
! ==  when rhoe(:,1) contains the alpha density instaed of        ==
! ==  the total density. This occurs for instance when computing  ==
! ==  the TDDFT forces with LSD                                   ==
! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    my_is_alpha_beta_dens=.FALSE.
    IF (PRESENT(is_alpha_beta_dens)) my_is_alpha_beta_dens=is_alpha_beta_dens
    ! ==--------------------------------------------------------------==
    IF (cntl%ttau) THEN
       CALL stopgm('VXC_CALC','META FUNCTIONALS NOT IMPLENTED',& 
            __LINE__,__FILE__)
    ENDIF
    IF (corel%tinlc) THEN
       CALL daxpy(clsd%nlsd*fpar%nnr1,1._real_8,roct,1,rhoe,1)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == COMPUTE EXCHANGE AND CORRELATION ENERGY (EXC)                ==
    ! ==--------------------------------------------------------------==
    CALL zeroing(v)!,maxfft*clsd%nlsd)
    CALL xcener(sxc,vxc,rhoe,rhoe,v)
    IF (cntl%tgc) THEN
       ! ..CALCULATE THE GRADIENT OF THE DENSITY
       IF (cntl%tlsd) THEN
          CALL  fwfftn(v(:,1),.FALSE.,parai%allgrp)
          CALL zgthr(ncpw%nhg,v(:,1),vtemp(:,1),nzh)
          CALL  fwfftn(v(:,2),.FALSE.,parai%allgrp)
          CALL zgthr(ncpw%nhg,v(:,2),vtemp(:,2),nzh)
          IF (.NOT.my_is_alpha_beta_dens) THEN
             DO ir=1,fpar%nnr1
                rhoe(ir,1)=rhoe(ir,1)-rhoe(ir,2)
             ENDDO
          ENDIF
          !write(*,*) 'max rhoe',sum(abs(rhoe(:,1))), sum(abs(rhoe(:,2)))
          CALL graden(rhoe(:,1),v,grad(:,1:4),vtmp)
          CALL graden(rhoe(:,2),v,grad(:,5:8),vtmp)
          !write(*,*) 'max grad',sum(abs(grad(:,1:4))),sum(abs(grad(:,5:8))),shape(grad)
       ELSE
          CALL  fwfftn(v(:,1),.FALSE.,parai%allgrp)
          CALL zgthr(ncpw%nhg,v,vtemp,nzh)
          CALL graden(rhoe,v,grad,vtmp)
          !write(*,*) 'max rhoe',sum(abs(rhoe(:,1)))
          !write(*,*) 'max grad',sum(abs(grad(:,1:4))),shape(grad)
       ENDIF
       ! ..GRADIENT CORRECTION TO THE EXCHANGE ENERGY
       ALLOCATE(rvtmp(fpar%nnr1,4*clsd%nlsx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       IF (cntl%tlsd) THEN
          CALL gclsd(sgcx,sgcc,rhoe,v,vtemp,rvtmp,grad,.FALSE.)
       ELSE
          CALL gcener(sgcx,sgcc,rhoe,v,vtemp,rvtmp(:,1),grad,.FALSE.)
       ENDIF
       DEALLOCATE(rvtmp,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE vxc_calc
  ! ==================================================================

END MODULE dd_xc_utils
