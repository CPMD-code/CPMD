MODULE v1ofrho1_utils
  USE cppt,                            ONLY: indz,&
                                             nzh,&
                                             scg
  USE dd_xc_ana_utils,                 ONLY: give_scr_dd_xc_ana
  USE dd_xc_utils,                     ONLY: dd_xc,&
                                             dd_xct
  USE error_handling,                  ONLY: stopgm
  USE fftmain_utils,                   ONLY: fwfftn,&
                                             invfftn
  USE fftnew_utils,                    ONLY: setfftn
  USE geq0mod,                         ONLY: geq0
  USE isos,                            ONLY: isos1,&
                                             isos3
  USE kinds,                           ONLY: real_8
  USE linres,                          ONLY: lr03,&
                                             lrf2,&
                                             lrsym,&
                                             td01,&
                                             td03,&
                                             tshl
  USE parac,                           ONLY: parai,&
                                             paral
  USE poin,                            ONLY: rhoo
  USE spin,                            ONLY: clsd,&
                                             lspin1,&
                                             lspin2
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             ncpw,&
                                             parm,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: v1ofrho1
  PUBLIC :: give_scr_v1ofrho1
  PUBLIC :: vtdofrho1
  PUBLIC :: give_scr_vtdofrho1

CONTAINS

  ! ==================================================================
  SUBROUTINE v1ofrho1(e2,drhoe,ddxc,psi)
    ! ==--------------------------------------------------------------==
    ! calculate the implicit part of the first order potential, i.e.
    ! the part induced by the response of the density wrt the bare
    ! perturbation.
    ! input:  first order density in
    ! drhoe(r) = <phi1|r><r|phi0> + cc
    ! output: first order potential (implicit part) in
    ! drhoe(r)
    ! ==--------------------------------------------------------------==
    ! 
    REAL(real_8)                             :: e2(:), drhoe(:,:), ddxc(:,:)
    COMPLEX(real_8)                          :: psi(:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'v1ofrho1'

    COMPLEX(real_8), ALLOCATABLE             :: v1(:), vtemp(:,:), vtmp(:)
    INTEGER                                  :: ierr, ig, il_grad_1d, &
                                                il_grad_2d, il_vtemp_1d, &
                                                il_vtemp_2d, il_vtmp, ir, &
                                                isub, nlsds
    REAL(real_8)                             :: dra, drb, dvma, dvmb, dvta, &
                                                dvtb, eht2, exc2, exc2m, &
                                                exc2t, tt
    REAL(real_8), ALLOCATABLE                :: dxc(:,:), grad(:,:), rhox(:,:)

! ==--------------------------------------------------------------==

    CALL tiset('   V1OFRHO1',isub)
    CALL setfftn(0)
    ! ==--------------------------------------------------------------==
    IF (lspin2%tlse) THEN
       nlsds=2         ! Number of spin densities treated together
    ELSE
       nlsds=clsd%nlsd
    ENDIF
    ALLOCATE(v1(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(dxc(fpar%nnr1, clsd%nlsd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(rhox(fpar%nnr1 , nlsds),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! 
    IF (cntl%tgc) THEN
       il_grad_1d=fpar%nnr1
       il_grad_2d=nlsds*4
       il_vtemp_1d=ncpw%nhg
       il_vtemp_2d=nlsds
       il_vtmp=MAX(2*ncpw%nhg,fpar%nnr1*2*nlsds)
    ELSE
       il_grad_1d=1
       il_grad_2d=1
       il_vtemp_1d=1
       il_vtemp_2d=1
       il_vtmp=1
    ENDIF
    ALLOCATE(vtemp(il_vtemp_1d,il_vtemp_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(vtmp(il_vtmp/2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(grad(il_grad_1d,il_grad_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! 
    IF (isos1%tclust.AND.isos3%ps_type.EQ.1)&
         CALL stopgm('V1OFRHO1','HOCKNEY NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    ! ..transform n1 to g-space
    CALL zeroing(psi(:,1))!,maxfft)
    DO ir=1,fpar%nnr1
       psi(ir,1)=CMPLX(drhoe(ir,1),0._real_8,kind=real_8)
    ENDDO
    CALL  fwfftn(psi(:,1),.FALSE.,parai%allgrp)
    ! ==--------------------------------------------------------------==
    ! == This is the Hartree Potential                                ==
    ! ==--------------------------------------------------------------==
    eht2=0._real_8
    DO ig=1,ncpw%nhg
       v1(ig) = psi(nzh(ig),1) * scg(ig)
       eht2 = eht2 + REAL(CONJG(v1(ig))*psi(nzh(ig),1))
    ENDDO
    eht2 = eht2*parm%omega
    ! ==--------------------------------------------------------------==
    ! == XC contributions                                             ==
    ! ==--------------------------------------------------------------==
    IF (lspin2%tlse) THEN
       ! Switch to LSD for mixed and triplet states
       lspin2%tlse=.FALSE.
       cntl%tlsd=.TRUE.
       clsd%nlsd=2
       IF (lr03%txc_analytic) THEN
          IF (paral%io_parent)&
               WRITE(6,*) 'V1OFRHO1: use XC_DD_ANALYTIC or numeric with LSE'
          CALL stopgm('V1OFRHO1','XC_ANALYTIC not implemented',& 
               __LINE__,__FILE__)
       ENDIF
       ! Mixed state
       CALL dd_xc(rhoo(:,1:2),drhoe(:,1:2),rhox,psi,dxc(:,1:2),&
            VTMP,VTEMP,GRAD,.FALSE.)
       exc2m=0._real_8
       DO ir=1,fpar%nnr1
          drb=drhoe(ir,2)
          dra=drhoe(ir,1)-drhoe(ir,2)
          exc2m=exc2m+dxc(ir,1)*dra+dxc(ir,2)*drb
       ENDDO
       exc2m=0.5_real_8*exc2m*parm%omega/REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)
       ! Triplet state
       CALL dd_xc(rhoo(:,3:4),drhoe(:,3:4),rhox,psi,dxc(:,3:4),&
            VTMP,VTEMP,GRAD,.FALSE.)
       exc2t=0._real_8
       DO ir=1,fpar%nnr1
          drb=drhoe(ir,4)
          dra=drhoe(ir,3)-drhoe(ir,4)
          exc2t=exc2t+dxc(ir,3)*dra+dxc(ir,4)*drb
       ENDDO
       exc2t=0.5_real_8*exc2t*parm%omega/REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)
       ! Join and switch back to LSE
       exc2 = lspin1%lsea*exc2m + lspin1%lseb*exc2t
       IF (lspin2%tlsets) THEN
          tt=1.0_real_8/1.5_real_8
          DO ir=1,fpar%nnr1
             dvma=dxc(ir,1)
             dvmb=dxc(ir,2)
             dvta=dxc(ir,3)
             dvtb=dxc(ir,4)
             dxc(ir,1)=0.5_real_8*(lspin1%lsea*(dvma+dvmb)+lspin1%lseb*(dvta+dvtb))
             dxc(ir,2)=tt*(lspin1%lsea*(dvma+0.5_real_8*dvmb)+lspin1%lseb*(0.5_real_8*dvta+dvtb))
             dxc(ir,3)=lspin1%lsea*dvmb+lspin1%lseb*dvtb
          ENDDO
       ELSE
          DO ir=1,fpar%nnr1
             dvma=dxc(ir,1)
             dvmb=dxc(ir,2)
             dvta=dxc(ir,3)
             dvtb=dxc(ir,4)
             dxc(ir,1)=0.5_real_8*(lspin1%lsea*(dvma+dvmb)+lspin1%lseb*(dvta+dvtb))
             dxc(ir,2)=lspin1%lsea*dvma+lspin1%lseb*dvtb
             dxc(ir,3)=lspin1%lsea*dvmb+lspin1%lseb*dvtb
          ENDDO
       ENDIF
       lspin2%tlse=.TRUE.
       cntl%tlsd=.FALSE.
       clsd%nlsd=4
    ELSEIF (cntl%tlsd) THEN
       IF (lr03%txc_analytic) THEN
          DO ir=1,fpar%nnr1
             drb=drhoe(ir,2)
             dra=drhoe(ir,1)-drhoe(ir,2)
             dxc(ir,1)=ddxc(ir,1)*dra+ddxc(ir,3)*drb
             dxc(ir,2)=ddxc(ir,2)*drb+ddxc(ir,3)*dra
          ENDDO
       ELSE
          CALL dd_xc(rhoo,drhoe,rhox,psi,dxc,vtmp,vtemp,grad,&
               .FALSE.)
       ENDIF
       exc2=0._real_8
       DO ir=1,fpar%nnr1
          drb=drhoe(ir,2)
          dra=drhoe(ir,1)-drhoe(ir,2)
          exc2=exc2+dxc(ir,1)*dra+dxc(ir,2)*drb
       ENDDO
       exc2=0.5_real_8*exc2*parm%omega/REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)
    ELSE
       IF (lr03%txc_analytic) THEN
          DO ir=1,fpar%nnr1
             dxc(ir,1)=ddxc(ir,1)*drhoe(ir,1)
          ENDDO
       ELSE
          CALL dd_xc(rhoo,drhoe,rhox,psi,dxc,vtmp,vtemp,grad,&
               .FALSE.)
       ENDIF
       exc2=0._real_8
       DO ir=1,fpar%nnr1
          exc2=exc2+dxc(ir,1)*drhoe(ir,1)
       ENDDO
       exc2=0.5_real_8*exc2*parm%omega/REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL zeroing(psi(:,1))!,maxfft)
    !CDIR NODEP
    DO ig=1,ncpw%nhg
       psi(indz(ig),1) = CONJG(v1(ig))
       psi(nzh(ig),1)  = v1(ig)
    ENDDO
    IF (geq0)   psi(nzh(1),1) = v1(1)
    CALL  invfftn(psi(:,1),.FALSE.,parai%allgrp)
    ! ==--------------------------------------------------------------==
    ! == Total potential in drhoe                                     ==
    ! ==--------------------------------------------------------------==
    IF (lspin2%tlse) THEN
       DO ir=1,fpar%nnr1
          drhoe(ir,1)=REAL(psi(ir,1))+dxc(ir,1)
          drhoe(ir,2)=REAL(psi(ir,1))+dxc(ir,2)
          drhoe(ir,3)=REAL(psi(ir,1))+dxc(ir,3)
       ENDDO
    ELSEIF (cntl%tlsd) THEN
       DO ir=1,fpar%nnr1
          drhoe(ir,1)=REAL(psi(ir,1))+dxc(ir,1)
          drhoe(ir,2)=REAL(psi(ir,1))+dxc(ir,2)
       ENDDO
    ELSE
       DO ir=1,fpar%nnr1
          drhoe(ir,1)=REAL(psi(ir,1))+dxc(ir,1)
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    e2(1) = eht2
    e2(2) = exc2
    ! ==--------------------------------------------------------------==
    DEALLOCATE(v1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(dxc,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(rhox,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(vtemp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(vtmp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(grad,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL tihalt('   V1OFRHO1',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE v1ofrho1
  ! ==================================================================
  SUBROUTINE give_scr_v1ofrho1(lv1ofrho1,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lv1ofrho1
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: il_dxc, il_grad, il_rhox, &
                                                il_v1, il_vtemp, il_vtmp, &
                                                ildd_xc_ana, NLSDS

! ==--------------------------------------------------------------==

    IF (lspin2%tlse) THEN
       nlsds=2
    ELSE
       nlsds=clsd%nlsd
    ENDIF
    IF (cntl%tgc) THEN
       il_grad=fpar%nnr1*nlsds*4
       il_vtemp=2*ncpw%nhg*nlsds
       il_vtmp=MAX(2*ncpw%nhg,fpar%nnr1*(2*nlsds-1))
    ELSE
       il_grad=1
       il_vtemp=1
       il_vtmp=1
    ENDIF
    IF (lr03%txc_dd_ana) THEN
       CALL give_scr_dd_xc_ana(ildd_xc_ana,tag)
    ELSE
       ildd_xc_ana=1
    ENDIF
    il_dxc=fpar%nnr1*clsd%nlsd
    il_rhox=fpar%nnr1*nlsds
    il_v1=2*ncpw%nhg
    lv1ofrho1=il_grad+il_vtemp+il_vtmp+il_dxc+il_v1+il_rhox+&
         ILDD_XC_ANA
    lv1ofrho1=lv1ofrho1+200
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_v1ofrho1
  ! ==================================================================
  SUBROUTINE vtdofrho1(e2,drhoe,ddxc,psi,swfun)
    ! ==--------------------------------------------------------------==
    ! 
    REAL(real_8)                             :: e2(:), drhoe(:,:), ddxc(:,:)
    COMPLEX(real_8)                          :: psi(:,:)
    LOGICAL                                  :: swfun

    CHARACTER(*), PARAMETER                  :: procedureN = 'vtdofrho1'

    COMPLEX(real_8), ALLOCATABLE             :: v1(:), vtemp(:,:), vtmp(:)
    INTEGER                                  :: ierr, ig, il_grad_1d, &
                                                il_grad_2d, il_vtemp_1d, &
                                                il_vtemp_2d, il_vtmp, ir, &
                                                isub, nmu
    REAL(real_8)                             :: dra, drb, eht2, exc2
    REAL(real_8), ALLOCATABLE                :: dxc(:,:), grad(:,:), rhox(:,:)

! 
! ==--------------------------------------------------------------==

    CALL tiset(' VTDOFRHO1',isub)
    ! ==--------------------------------------------------------------==
    eht2=0._real_8
    exc2=0._real_8
    ! 
    IF (cntl%tlsd.OR.lrsym.EQ.3) THEN
       nmu=2
    ELSE
       nmu=1
    ENDIF
    ALLOCATE(v1(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(dxc(fpar%nnr1, nmu),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(rhox(fpar%nnr1, nmu),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    IF ((cntl%tgc.OR.lrf2%td_tgc).AND..NOT.lr03%txc_analytic) THEN
       il_grad_1d=fpar%nnr1
       il_grad_2d=nmu*4
       il_vtemp_1d=ncpw%nhg
       il_vtemp_2d=nmu
       il_vtmp=MAX(2*ncpw%nhg,fpar%nnr1*(2*nmu-1))
    ELSE
       il_grad_1d=1
       il_grad_2d=1
       il_vtemp_1d=1
       il_vtemp_2d=1
       il_vtmp=1
    ENDIF
    ALLOCATE(vtemp(il_vtemp_1d,il_vtemp_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(vtmp(il_vtmp/2+1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(grad(il_grad_1d,il_grad_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! 
    IF (isos1%tclust.AND.isos3%ps_type.EQ.1)&
         CALL stopgm('VTDOFRHO1','HOCKNEY NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    ! ..first: Hartree (only if not triplet)
    IF (lrsym.NE.3) THEN
       ! ..transform n1 to g-space
       CALL zeroing(psi)!,maxfft)
       DO ir=1,fpar%nnr1
          psi(ir,1)=CMPLX(drhoe(ir,1),0._real_8,kind=real_8)
       ENDDO
       CALL  fwfftn(psi(:,1),.FALSE.,parai%allgrp)
       ! ==--------------------------------------------------------------==
       ! == This is the Hartree Potential                                ==
       ! ==--------------------------------------------------------------==
       DO ig=1,ncpw%nhg
          v1(ig) = psi(nzh(ig),1) * scg(ig)
          eht2 = eht2 + REAL(CONJG(v1(ig))*psi(nzh(ig),1))
       ENDDO
    ENDIF
    eht2 = eht2*parm%omega
    ! ==--------------------------------------------------------------==
    ! == XC contributions                                             ==
    ! ==--------------------------------------------------------------==
    IF (cntl%tlsd) THEN
       IF (lr03%txc_analytic) THEN
          DO ir=1,fpar%nnr1
             drb=drhoe(ir,2)
             dra=drhoe(ir,1)-drhoe(ir,2)
             dxc(ir,1)=ddxc(ir,1)*dra+ddxc(ir,3)*drb
             dxc(ir,2)=ddxc(ir,2)*drb+ddxc(ir,3)*dra
          ENDDO
       ELSE
          CALL dd_xc(rhoo,drhoe,rhox,psi,dxc,vtmp,vtemp,grad,&
               swfun)
       ENDIF
       DO ir=1,fpar%nnr1
          drb=drhoe(ir,2)
          dra=drhoe(ir,1)-drhoe(ir,2)
          exc2=exc2+dxc(ir,1)*dra+dxc(ir,2)*drb
       ENDDO
       exc2=0.5_real_8*exc2*parm%omega/REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)
    ELSEIF (lrsym.EQ.1) THEN
       ! ..SINGLET
       IF (lr03%txc_analytic) THEN
          DO ir=1,fpar%nnr1
             dxc(ir,1)=ddxc(ir,1)*drhoe(ir,1)
          ENDDO
       ELSE
          CALL dd_xc(rhoo,drhoe,rhox,psi,dxc,vtmp,vtemp,grad,&
               swfun)
       ENDIF
       DO ir=1,fpar%nnr1
          exc2=exc2+dxc(ir,1)*drhoe(ir,1)
       ENDDO
       exc2=0.5_real_8*exc2*parm%omega/REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)
    ELSEIF (lrsym.EQ.3) THEN
       ! ..TRIPLET
       IF (lr03%txc_analytic) THEN
          DO ir=1,fpar%nnr1
             dxc(ir,1)=ddxc(ir,1)*drhoe(ir,1)
          ENDDO
       ELSE
          CALL dd_xct(rhoo,drhoe,rhox,psi,dxc,vtmp,vtemp,grad,swfun)
       ENDIF
       DO ir=1,fpar%nnr1
          exc2=exc2+dxc(ir,1)*drhoe(ir,1)
       ENDDO
       exc2=0.5_real_8*exc2*parm%omega/REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)
    ELSE
       CALL stopgm('VTDOFRHO1','WRONG LRSYM',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL zeroing(psi(:,1))!,maxfft)
    IF (lrsym.NE.3) THEN
       !CDIR NODEP
       DO ig=1,ncpw%nhg
          psi(indz(ig),1) = CONJG(v1(ig))
          psi(nzh(ig),1)  = v1(ig)
       ENDDO
       IF (geq0)   psi(nzh(1),1) = v1(1)
       CALL  invfftn(psi(:,1),.FALSE.,parai%allgrp)
    ENDIF
    e2(1) = eht2
    e2(2) = exc2

    ! ==--------------------------------------------------------------==
    ! == Total potential in drhoe                                     ==
    ! ==--------------------------------------------------------------==
    IF (cntl%tlsd) THEN
       DO ir=1,fpar%nnr1
          drhoe(ir,1)=0.5_real_8*(REAL(psi(ir,1))+dxc(ir,1))
          drhoe(ir,2)=0.5_real_8*(REAL(psi(ir,1))+dxc(ir,2))
       ENDDO
       IF (.NOT.td03%tda) CALL dscal(fpar%nnr1,2._real_8,drhoe(1,1),1)
       IF (.NOT.td03%tda) CALL dscal(fpar%nnr1,2._real_8,drhoe(1,2),1)
    ELSE
       DO ir=1,fpar%nnr1
          drhoe(ir,1)=REAL(psi(ir,1))+dxc(ir,1)
       ENDDO
       IF (.NOT.td03%tda) CALL dscal(fpar%nnr1,2._real_8,drhoe(1,1),1)
    ENDIF

    DEALLOCATE(v1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(dxc,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(rhox,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(vtemp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(vtmp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(grad,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL tihalt(' VTDOFRHO1',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE vtdofrho1
  ! ==================================================================
  SUBROUTINE give_scr_vtdofrho1(lvtdofrho1,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lvtdofrho1
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: il_dxc, il_grad, il_rhox, &
                                                il_v1, il_vtemp, il_vtmp, &
                                                LDD_XC_ANA, nmu
    LOGICAL                                  :: tlsdsv

! ==--------------------------------------------------------------==

    IF (cntl%tlsd.OR.td01%ns_tri.GT.0) THEN
       nmu=2
    ELSE
       nmu=1
    ENDIF
    IF (tshl%isc) nmu=2
    il_dxc=fpar%nnr1*nmu
    il_rhox=fpar%nnr1*nmu
    il_v1=2*ncpw%nhg
    IF ((cntl%tgc.OR.lrf2%td_tgc).AND..NOT.lr03%txc_analytic) THEN
       il_grad=fpar%nnr1*nmu*4
       il_vtemp=2*ncpw%nhg*nmu
       il_vtmp=MAX(2*ncpw%nhg,fpar%nnr1*(2*nmu-1))
    ELSE
       il_grad=1
       il_vtemp=1
       il_vtmp=1
    ENDIF
    IF (lr03%txc_dd_ana) THEN
       IF (td01%ns_tri.GT.0.OR.tshl%isc) THEN
          tlsdsv=cntl%tlsd
          cntl%tlsd=.TRUE.
          CALL give_scr_dd_xc_ana(ldd_xc_ana,tag)
          cntl%tlsd=tlsdsv
       ELSE
          CALL give_scr_dd_xc_ana(ldd_xc_ana,tag)
       ENDIF
    ELSE
       ldd_xc_ana=1
    ENDIF
    lvtdofrho1=il_grad+il_vtemp+il_vtmp+il_dxc+il_v1+il_rhox+&
         LDD_XC_ANA
    lvtdofrho1=lvtdofrho1+200
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_vtdofrho1
  ! ==================================================================


END MODULE v1ofrho1_utils
