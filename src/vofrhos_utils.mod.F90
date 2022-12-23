MODULE vofrhos_utils
  USE cp_xc_utils,                     ONLY: cp_xc
  USE cnst,                            ONLY: fpi
  USE cofor_utils,                     ONLY: cofor
  USE corec_utils,                     ONLY: corec
  USE cppt,                            ONLY: hg,&
                                             indz,&
                                             nzh,&
                                             scg
  USE efld,                            ONLY: extf,&
                                             textfld
  USE ener,                            ONLY: ener_c,&
                                             ener_com,&
                                             ener_d
  USE error_handling,                  ONLY: stopgm
  USE extpotmod,                       ONLY: extpot
  USE fftmain_utils,                   ONLY: fwfftn,&
                                             invfftn
  USE fftnew_utils,                    ONLY: setfftn
  USE func,                            ONLY: func1,&
                                             mfxcc_is_pade,&
                                             mfxcc_is_skipped,&
                                             mfxcx_is_skipped,&
                                             mgcc_is_skipped,&
                                             mgcx_is_skipped
  USE gcener_utils,                    ONLY: gclsd
  USE geq0mod,                         ONLY: geq0
  USE graden_utils,                    ONLY: graden
  USE isos,                            ONLY: isos1,&
                                             isos3
  USE kinds,                           ONLY: real_8
  USE nlcc,                            ONLY: corel,&
                                             corer,&
                                             vnlt
  USE parac,                           ONLY: parai
  USE potfor_utils,                    ONLY: potabfor
  USE pslo,                            ONLY: pslo_com
  USE spin,                            ONLY: clsd,&
                                             lspin1,&
                                             lspin2
  USE str2,                            ONLY: dqg
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             maxsys,&
                                             ncpw,&
                                             parm,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: zgthr
  USE vofrhoc_utils,                   ONLY: potab
  USE xcener_utils,                    ONLY: xcener
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: vofrhos
  PUBLIC :: give_scr_vofrhos
  PUBLIC :: vofrhor
  PUBLIC :: give_scr_vofrhor

CONTAINS

  ! ==================================================================
  SUBROUTINE vofrhos(fion,rhoe,v,vtemp,tfor,tstress)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==  THE ONE-PARTICLE POTENTIAL V IN REAL SPACE                  ==
    ! ==  THE TOTAL ENERGY ETOT                                       ==
    ! ==  THE FORCES FION ACTING ON THE IONS                          ==
    ! ==  SUBROUTINE FOR LOW SPIN EXCITATION                          ==
    ! ==--------------------------------------------------------------==
    ! == RHOE:  in electronic density in real space                   ==
    ! ==       out total electronic potential in real space           ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: fion(:,:,:)
    REAL(real_8), TARGET                     :: rhoe(:,:)
    COMPLEX(real_8)                          :: v(:,:), vtemp(:,:)
    LOGICAL                                  :: tfor, tstress

    CHARACTER(*), PARAMETER                  :: procedureN = 'vofrhos'

    COMPLEX(real_8)                          :: epotab
    COMPLEX(real_8), ALLOCATABLE             :: vtmp(:,:)
    INTEGER                                  :: ierr, ig, il_grad, il_rhoval, &
                                                il_vtmp, ir, isub, length, &
                                                nnrs
    REAL(real_8)                             :: od, sgcc, sgcc1, sgcc2, sgcx, &
                                                sgcx1, sgcx2, sxc, sxc1, &
                                                sxc2, vcc, vgc, vgc1, vgc2, &
                                                vma, vmb, vta, vtb, vxc1, vxc2
    REAL(real_8), ALLOCATABLE                :: grad(:,:), rvtmp(:,:)
    REAL(real_8), POINTER                    :: rhoval(:,:)

    CALL tiset(procedureN,isub)
    CALL setfftn(0)
    IF (tstress) THEN
       CALL stopgm(procedureN,'STRESS NOT IMPLEMENTED',& 
            __LINE__,__FILE__)
    ENDIF
    IF (cntl%ttau) THEN
       CALL stopgm(procedureN,'META FUNCTIONALS NOT IMPLENTED',& 
            __LINE__,__FILE__)
    ENDIF
    ! SCR partition.
    il_vtmp = MAX(fpar%nnr1,ncpw%nhg*2)*clsd%nlsx      ! GCENER and COREC
    IF (cntl%tgc) THEN
       il_grad=fpar%nnr1*8   ! GRADEN
    ELSE
       il_grad=0
    ENDIF
    IF (corel%tinlc) THEN
       il_rhoval=fpar%nnr1*clsd%nlsd   ! XCENER
    ELSE
       il_rhoval=0
    ENDIF
    length = il_vtmp+MAX(il_grad,il_rhoval)
    IF (pslo_com%tivan) length=MAX(length,ncpw%nhg*4+maxsys%nax) ! NEWD
    ALLOCATE(vtmp(ncpw%nhg,il_vtmp/ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(rvtmp(fpar%nnr1,clsd%nlsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    IF (cntl%tgc) ALLOCATE(grad(fpar%nnr1,8),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    IF (corel%tinlc) THEN
       ALLOCATE(rhoval(fpar%nnr1,clsd%nlsd),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
    ELSE
       rhoval => rhoe
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == ADD CORE CHARGE TO RHOE                                      ==
    ! ==--------------------------------------------------------------==
    IF (corel%tinlc) THEN
       CALL dcopy(fpar%nnr1*clsd%nlsd,rhoe(1,1),1,rhoval(1,1),1)
       CALL corec(rhoe,vtmp,v)
       CALL dcopy(2*ncpw%nhg,vtemp(1,1),1,vnlt(1),1)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == VTEMP (Potential in G-Space) -FFT-> V(R)                     ==
    ! ==--------------------------------------------------------------==
    CALL zeroing(v)!,maxfft*clsd%nlsd)
    !CDIR NODEP
    DO ig=1,ncpw%nhg
       v(indz(ig),1) = CONJG(vtemp(ig,1))
       v(nzh(ig),1)  = vtemp(ig,1)
    ENDDO
    IF (geq0) THEN
       v(nzh(1),1) = vtemp(1,1)
    ENDIF
    CALL  invfftn(v(:,1),.FALSE.,parai%allgrp)
    ! ==--------------------------------------------------------------==
    ! == ADD EXTERNAL POTENTIAL TO V                                  ==
    ! ==--------------------------------------------------------------==
    IF (textfld) THEN
       ! ..potential from classical interface
       !$omp parallel do private(IR)
       DO ir=1,fpar%nnr1
          v(ir,1)=v(ir,1)+CMPLX(extf(ir),0._real_8,kind=real_8)
       ENDDO
    ENDIF
    IF (cntl%texpot) THEN
       ! ..static external potential
       !$omp parallel do private(IR)
       DO ir=1,fpar%nnr1
          v(ir,1)=v(ir,1)+CMPLX(extpot(ir),0._real_8,kind=real_8)
       ENDDO
    ENDIF
    ! SAVE THE GROUND STATE DENSITY
    IF (lspin2%tcas22) THEN
       ! Warning DQG is complex (address (1,2))
       CALL dcopy(fpar%nnr1,rhoe(1,3),1,dqg(1,2),1)
       CALL daxpy(fpar%nnr1,-1._real_8,rhoe(1,4),1,dqg(1,2),1)
       CALL dscal(fpar%nnr1,2._real_8,dqg(1,2),1)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == MIXED STATE                                                  ==
    ! ==--------------------------------------------------------------==
    ! == COMPUTE EXCHANGE AND CORRELATION ENERGY (EXC)                ==
    ! ==--------------------------------------------------------------==
    CALL dcopy(2*fpar%nnr1,v(1,1),1,dqg(1,1),1)
    CALL zeroing(v)!,nnr1)
    sxc1  = 0.0_real_8
    vxc1  = 0.0_real_8
    sgcx1 = 0.0_real_8
    sgcc1 = 0.0_real_8
    vgc1  = 0.0_real_8
    IF (ABS(lspin1%lsea) .GT. 1.e-6_real_8) THEN
       cntl%tlsd=.TRUE.
       IF (ASSOCIATED(cp_xc)) cp_xc%polarized = .NOT. cp_xc%polarized
       CALL xcener(sxc1,vxc1,rhoval,rhoe(:,1:2),v)
       cntl%tlsd=.FALSE.
       IF (cntl%tgc) THEN
          ! CALCULATE THE GRADIENT OF THE DENSITY
          CALL fwfftn(v(:,1),.FALSE.,parai%allgrp)
          CALL zgthr(ncpw%nhg,v(:,1),vtemp(:,1),nzh)
          CALL fwfftn(v(:,2),.FALSE.,parai%allgrp)
          CALL zgthr(ncpw%nhg,v(:,2),vtemp(:,2),nzh)
          DO ir=1,fpar%nnr1
             rhoe(ir,1)=rhoe(ir,1)-rhoe(ir,2)
          ENDDO
          CALL graden(rhoe(:,1),v,grad(1,1),vtmp)
          CALL graden(rhoe(:,2),v,grad(1,5),vtmp)
          ! GRADIENT CORRECTION TO THE EXCHANGE ENERGY (EGCX)
          CALL gclsd(sgcx1,sgcc1,rhoe,v,vtemp,rvtmp,grad,.FALSE.)
       ENDIF
       IF (ASSOCIATED(cp_xc)) cp_xc%polarized = .NOT. cp_xc%polarized
    ENDIF
    ! V CONTAINS THE MIXED STATE VXC IN R-SPACE, MOVE IT TO RHOE
    !$omp parallel do private(IR)
    DO ir=1,fpar%nnr1
       rhoe(ir,1)=REAL(v(ir,1))
       rhoe(ir,2)=REAL(v(ir,2))
    ENDDO
    ! ==--------------------------------------------------------------==
    ! == TRIPLET STATE                                                ==
    ! ==--------------------------------------------------------------==
    ! == COMPUTE EXCHANGE AND CORRELATION ENERGY (EXC)                ==
    ! ==--------------------------------------------------------------==
    CALL zeroing(v)!,nnr1)
    sxc2  = 0.0_real_8
    vxc2  = 0.0_real_8
    sgcx2 = 0.0_real_8
    sgcc2 = 0.0_real_8
    vgc2  = 0.0_real_8
    IF (ABS(lspin1%lseb) .GT. 1.e-6_real_8) THEN
       cntl%tlsd=.TRUE.
       IF (ASSOCIATED(cp_xc)) cp_xc%polarized = .NOT. cp_xc%polarized
       CALL xcener(sxc2,vxc2,rhoval,rhoe(:,3:4),v)
       cntl%tlsd=.FALSE.
       IF (cntl%tgc) THEN
          ! CALCULATE THE GRADIENT OF THE DENSITY
          CALL fwfftn(v(:,1),.FALSE.,parai%allgrp)
          CALL zgthr(ncpw%nhg,v(:,1),vtemp(:,1),nzh)
          CALL fwfftn(v(:,2),.FALSE.,parai%allgrp)
          CALL zgthr(ncpw%nhg,v(:,2),vtemp(:,2),nzh)
          !$omp parallel do private(IR)
          DO ir=1,fpar%nnr1
             rhoe(ir,3)=rhoe(ir,3)-rhoe(ir,4)
          ENDDO
          CALL graden(rhoe(:,3),v,grad(1,1),vtmp)
          CALL graden(rhoe(:,4),v,grad(1,5),vtmp)
          ! GRADIENT CORRECTION TO THE EXCHANGE ENERGY (EGCX)
          CALL gclsd(sgcx2,sgcc2,rhoe(:,3:),v,vtemp,rvtmp,grad,.FALSE.)
       ENDIF
       IF (ASSOCIATED(cp_xc)) cp_xc%polarized = .NOT. cp_xc%polarized
    ENDIF
    ! V CONTAINS THE TRIPLET STATE VXC IN R-SPACE, MOVE IT TO RHOE
    !$omp parallel do private(IR)
    DO ir=1,fpar%nnr1
       rhoe(ir,3)=REAL(v(ir,1))
       rhoe(ir,4)=REAL(v(ir,2))
    ENDDO
    IF (lspin2%troot) THEN
       !$omp parallel do private(IR,VMA,VMB,VTA,VTB,VCC)
       DO ir=1,fpar%nnr1
          vma=rhoe(ir,1)
          vmb=rhoe(ir,2)
          vta=rhoe(ir,3)
          vtb=rhoe(ir,4)
          vcc=REAL(dqg(ir,1))
          rhoe(ir,1)=vcc+0.5_real_8*(lspin1%lsea*vma+lspin1%lsea*vmb+lspin1%lseb*vta+lspin1%lseb*vtb)
          rhoe(ir,2)=0.5_real_8*(-lspin1%lsea*vma+lspin1%lsea*vmb+lspin1%lseb*vta-lspin1%lseb*vtb)
          rhoe(ir,3)=0.5_real_8*(lspin1%lsea*vma-lspin1%lsea*vmb+lspin1%lseb*vta-lspin1%lseb*vtb)
       ENDDO
    ELSE IF (lspin2%tlsets) THEN
       !$omp parallel do private(IR,VMA,VMB,VTA,VTB,VCC) 
       DO ir=1,fpar%nnr1
          vma=rhoe(ir,1)
          vmb=rhoe(ir,2)
          vta=rhoe(ir,3)
          vtb=rhoe(ir,4)
          vcc=REAL(dqg(ir,1))*(lspin1%lsea+lspin1%lseb)
          rhoe(ir,1)=vcc+0.5_real_8*(lspin1%lsea*(vma+vmb)+lspin1%lseb*(vta+vtb))
          rhoe(ir,2)=vcc+&
               (lspin1%lsea*(VMA+0.5_real_8*VMB)+lspin1%lseb*(0.5_real_8*VTA+VTB))/1.5_real_8
          rhoe(ir,3)=vcc+lspin1%lsea*vmb+lspin1%lseb*vtb
       ENDDO
    ELSE
       !$omp parallel do private(IR,VMA,VMB,VTA,VTB,VCC) 
       DO ir=1,fpar%nnr1
          vma=rhoe(ir,1)
          vmb=rhoe(ir,2)
          vta=rhoe(ir,3)
          vtb=rhoe(ir,4)
          vcc=REAL(dqg(ir,1))*(lspin1%lsea+lspin1%lseb)
          rhoe(ir,1)=vcc+0.5_real_8*(lspin1%lsea*(vma+vmb)+lspin1%lseb*(vta+vtb))
          rhoe(ir,2)=vcc+lspin1%lsea*vma+lspin1%lseb*vtb
          rhoe(ir,3)=vcc+lspin1%lsea*vmb+lspin1%lseb*vtb
       ENDDO
    ENDIF
    IF (lspin2%tcas22) THEN
       CALL dcopy(fpar%nnr1,dqg(1,2),1,rhoe(1,4),1)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == ADD THE PIECES                                               ==
    ! ==--------------------------------------------------------------==
    sxc=lspin1%lsea*sxc1+lspin1%lseb*sxc2
    ener_com%vxc=lspin1%lsea*vxc1+lspin1%lseb*vxc2
    vgc=lspin1%lsea*vgc1+lspin1%lseb*vgc2
    sgcx=lspin1%lsea*sgcx1+lspin1%lseb*sgcx2
    sgcc=lspin1%lsea*sgcc1+lspin1%lseb*sgcc2
    ! ==--------------------------------------------------------------==
    ! == CORRECT ENERGY FOR CORE CHARGE                               ==
    ! ==--------------------------------------------------------------==
    IF (corel%tinlc) THEN
       sxc=sxc-corer%excco
       sgcx=sgcx-corer%egcxco
       sgcc=sgcc-corer%egccco
    ENDIF
    IF (pslo_com%tivan.OR.(corel%tinlc.AND.tfor)) THEN
       ! ==--------------------------------------------------------------==
       ! == TRANSF. TOTAL POTENTIAL IN G-SPACE                           ==
       ! == PUT THE TOTAL POTENTIAL IN G-SPACE INTO VTEMP                ==
       ! ==--------------------------------------------------------------==
       CALL  fwfftn(v(:,1),.FALSE.,parai%allgrp)
       CALL zgthr(ncpw%nhg,v,vtemp,nzh)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == FORCE ON IONS DUE TO NLCC                                    ==
    ! ==--------------------------------------------------------------==
    IF (corel%tinlc.AND.tfor) CALL cofor(fion,vtemp)
    ! ==--------------------------------------------------------------==
    IF (lspin2%tpenal) THEN
       !$omp parallel do private(IR)
       DO ir=1,fpar%nnr1
          v(ir,1)=CMPLX(rhoe(ir,5),0._real_8,kind=real_8)
       ENDDO
       CALL  fwfftn(v(:,1),.FALSE.,parai%allgrp)
       CALL zgthr(ncpw%nhg,v,vtemp,nzh)
       IF (tfor) CALL potabfor(fion,vtemp)
       CALL potab(epotab,vtemp)
       ener_c%eht_ab=REAL(epotab)*parm%omega
       CALL zeroing(v)!,maxfft)
       DO ig=1,ncpw%nhg
          v(indz(ig),1) = CONJG(vtemp(ig,1))
          v(nzh(ig),1)  = vtemp(ig,1)
       ENDDO
       IF (geq0) THEN
          v(nzh(1),1) = vtemp(1,1)
       ENDIF
       CALL  invfftn(v(:,1),.FALSE.,parai%allgrp)
       !$omp parallel do private(IR)
       DO ir=1,fpar%nnr1
          rhoe(ir,5)=REAL(v(ir,1))
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == SINGLE PIECES OF THE ENERGY:                                 ==
    ! ==--------------------------------------------------------------==
    nnrs = spar%nr1s*spar%nr2s*spar%nr3s
    od = parm%omega/REAL(nnrs,kind=real_8)
    ener_com%egc=(sgcc+sgcx)*od
    ener_com%exc=sxc*od + ener_com%egc
    ener_com%vxc=(ener_com%vxc + vgc)*od
    ! Triplet energy
    IF (corel%tinlc) THEN
       ener_d%etot_t = (sxc2-corer%excco+sgcx2-corer%egcxco+sgcc2-corer%egccco)*od
    ELSE
       ener_d%etot_t = (sxc2+sgcx2+sgcc2)*od
    ENDIF
    !
    !
    DEALLOCATE(rvtmp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)

    ! 
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE vofrhos
  ! ==================================================================
  SUBROUTINE give_scr_vofrhos(lvofrhos,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lvofrhos
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: ltgc, ltinlc

    IF (cntl%tgc) THEN
       ltgc=fpar%nnr1*clsd%nlsd*4      ! GRADEN
    ELSE
       ltgc=0
    ENDIF
    IF (corel%tinlc) THEN            ! XCENER
       ltinlc=fpar%nnr1*clsd%nlsd
    ELSE
       ltinlc=0
    ENDIF
    lvofrhos = MAX(fpar%nnr1,ncpw%nhg*2)*clsd%nlsx+MAX(ltgc,ltinlc) ! GCENER and COREC
    IF (pslo_com%tivan) lvofrhos=MAX(lvofrhos,ncpw%nhg*4+maxsys%nax) ! NEWD
    tag='MAX(NNR1,NHG*2)*NLSX+.........'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_vofrhos
  ! ==================================================================
  SUBROUTINE vofrhor(fion,rhoe,v,vtemp,tfor,tstress)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==  THE ONE-PARTICLE POTENTIAL V IN REAL SPACE                  ==
    ! ==  THE TOTAL ENERGY ETOT                                       ==
    ! ==  THE FORCES FION ACTING ON THE IONS                          ==
    ! ==  SUBROUTINE FOR LOW SPIN EXCITATION (ROSS TYPE)              ==
    ! ==--------------------------------------------------------------==
    ! == RHOE:  in electronic density in real space                   ==
    ! ==       out total electronic potential in real space           ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: fion(:,:,:)
    REAL(real_8), TARGET                     :: rhoe(:,:)
    COMPLEX(real_8)                          :: v(:,:), vtemp(:,:)
    LOGICAL                                  :: tfor, tstress

    CHARACTER(*), PARAMETER                  :: procedureN = 'vofrhor'

    COMPLEX(real_8)                          :: epotab, scgy
    COMPLEX(real_8), ALLOCATABLE             :: vtmp(:,:)
    INTEGER                                  :: ierr, ig, ig1, il_grad, &
                                                il_rhoval, il_vtmp, ir, isub, &
                                                length, mfxcc_1, mfxcx_1, &
                                                mgcc_1, mgcx_1, nnrs
    REAL(real_8)                             :: ek, g2, sgcc, sgcc1, sgcc2, &
                                                sgcx, sgcx1, sgcx2, sxc, &
                                                sxc1, sxc2, vgc, vgc1, vgc2, &
                                                vp1, vp2, vp3, vxc1, vxc2
    REAL(real_8), ALLOCATABLE                :: grad(:,:), rvtmp(:,:)
    REAL(real_8), POINTER                    :: rhoval(:,:)

    CALL tiset('   VOFRHOR',isub)
    IF (tstress) THEN
       CALL stopgm('VOFRHOR','STRESS NOT IMPLEMENTED',& 
            __LINE__,__FILE__)
    ENDIF
    IF (cntl%ttau) THEN
       CALL stopgm('VOFRHOR','META FUNCTIONALS NOT IMPLENTED',& 
            __LINE__,__FILE__)
    ENDIF
    IF (pslo_com%tivan) THEN
       CALL stopgm('VOFRHOR','VDB PP NOT IMPLEMENTED',& 
            __LINE__,__FILE__)
    ENDIF
    ! SCR partition.
    il_vtmp = MAX(fpar%nnr1,ncpw%nhg*2)*clsd%nlsx      ! GCENER and COREC
    IF (cntl%tgc) THEN
       il_grad=fpar%nnr1*8   ! GRADEN
    ELSE
       il_grad=0
    ENDIF
    IF (corel%tinlc) THEN
       il_rhoval=fpar%nnr1*clsd%nlsd   ! XCENER
    ELSE
       il_rhoval=0
    ENDIF
    length = il_vtmp+MAX(il_grad,il_rhoval)
    ALLOCATE(vtmp(ncpw%nhg,il_vtmp/ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(rvtmp(fpar%nnr1,clsd%nlsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    IF (cntl%tgc) THEN
       ALLOCATE(grad(fpar%nnr1,8),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
    ENDIF
    IF (corel%tinlc) THEN
       ALLOCATE(rhoval(fpar%nnr1,clsd%nlsd),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
    ELSE
       rhoval => rhoe
    ENDIF
    ! ==--------------------------------------------------------------==
    mfxcx_1 = func1%mfxcx
    mfxcc_1 = func1%mfxcc
    mgcx_1  = func1%mgcx
    mgcc_1  = func1%mgcc
    IF (func1%mfxcc == mfxcc_is_pade) THEN
       CALL stopgm("ROSS METHOD"," PADE FUNCTIONAL NOT POSSIBLE",& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == ADD CORE CHARGE TO RHOE                                      ==
    ! ==--------------------------------------------------------------==
    IF (corel%tinlc) THEN
       CALL dcopy(fpar%nnr1*clsd%nlsd,rhoe(1,1),1,rhoval(1,1),1)
       CALL corec(rhoe,vtmp,v)
       CALL dcopy(2*ncpw%nhg,vtemp(1,1),1,vnlt(1),1)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == VTEMP (Potential in G-Space) -FFT-> V(R)                     ==
    ! ==--------------------------------------------------------------==
    CALL zeroing(v)!,maxfft*clsd%nlsd)
    DO ig=1,ncpw%nhg
       v(indz(ig),1) = CONJG(vtemp(ig,1))
       v(nzh(ig),1)  = vtemp(ig,1)
    ENDDO
    IF (geq0) THEN
       v(nzh(1),1) = vtemp(1,1)
    ENDIF
    CALL  invfftn(v(:,1),.FALSE.,parai%allgrp)
    ! ==--------------------------------------------------------------==
    ! == ADD EXTERNAL POTENTIAL TO V                                  ==
    ! ==--------------------------------------------------------------==
    IF (textfld) THEN
       ! ..potential from classical interface
       !$omp parallel do private(IR)
       DO ir=1,fpar%nnr1
          v(ir,1)=v(ir,1)+CMPLX(extf(ir),0._real_8,kind=real_8)
       ENDDO
    ENDIF
    IF (cntl%texpot) THEN
       ! ..static external potential
       !$omp parallel do private(IR)
       DO ir=1,fpar%nnr1
          v(ir,1)=v(ir,1)+CMPLX(extpot(ir),0._real_8,kind=real_8)
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == MIXED STATE                                                  ==
    ! ==--------------------------------------------------------------==
    func1%mfxcx = mfxcx_is_skipped
    func1%mfxcc = mfxcc_1
    func1%mgcx  = mgcx_is_skipped
    func1%mgcc  = mgcc_1
    cntl%tgcx=.FALSE.
    IF (func1%mgcc /= mgcc_is_skipped) THEN
       cntl%tgcc=.TRUE.
    ELSE
       cntl%tgcc=.FALSE.
    ENDIF
    cntl%tgc=cntl%tgcc.OR.cntl%tgcx
    ! ==--------------------------------------------------------------==
    ! == COMPUTE CORRELATION ENERGY (EC)                              ==
    ! ==--------------------------------------------------------------==
    CALL dcopy(2*fpar%nnr1,v(1,1),1,dqg(1,1),1)
    cntl%tlsd=.TRUE.
    IF (ASSOCIATED(cp_xc)) cp_xc%polarized = .NOT. cp_xc%polarized
    CALL xcener(sxc1,vxc1,rhoval,rhoe(:,1:2),v)
    cntl%tlsd=.FALSE.
    sgcx1 = 0.0_real_8
    sgcc1 = 0.0_real_8
    vgc1  = 0.0_real_8
    IF (cntl%tgc) THEN
       ! CALCULATE THE GRADIENT OF THE DENSITY
       CALL fwfftn(v(:,1),.FALSE.,parai%allgrp)
       CALL zgthr(ncpw%nhg,v(:,1),vtemp(:,1),nzh)
       CALL fwfftn(v(:,2),.FALSE.,parai%allgrp)
       CALL zgthr(ncpw%nhg,v(:,2),vtemp(:,2),nzh)
       !$omp parallel do private(IR)
       DO ir=1,fpar%nnr1
          rhoe(ir,1)=rhoe(ir,1)-rhoe(ir,2)
       ENDDO
       CALL graden(rhoe(:,1),v,grad(1,1),vtmp)
       CALL graden(rhoe(:,2),v,grad(1,5),vtmp)
       ! GRADIENT CORRECTION TO THE EXCHANGE ENERGY (EGCX)
       CALL gclsd(sgcx1,sgcc1,rhoe,v,vtemp,rvtmp,grad,.FALSE.)
    ENDIF
    IF (ASSOCIATED(cp_xc)) cp_xc%polarized = .NOT. cp_xc%polarized
    ! V CONTAINS THE TOTAL POTENTIAL IN R-SPACE, MOVE IT TO RHOE
    !$omp parallel do private(IR)
    DO ir=1,fpar%nnr1
       rhoe(ir,1)=REAL(v(ir,1))
       rhoe(ir,2)=REAL(v(ir,2))
    ENDDO
    ! ==--------------------------------------------------------------==
    ! == TRIPLET STATE                                                ==
    ! ==--------------------------------------------------------------==
    func1%mfxcx = mfxcx_1
    func1%mfxcc = mfxcc_is_skipped
    func1%mgcx  = mgcx_1
    func1%mgcc  = 0
    cntl%tgcc=.FALSE.
    IF (func1%mgcx /= mgcx_is_skipped) THEN
       cntl%tgcx=.TRUE.
    ELSE
       cntl%tgcx=.FALSE.
    ENDIF
    cntl%tgc=cntl%tgcc.OR.cntl%tgcx
    ! ==--------------------------------------------------------------==
    ! == COMPUTE EXCHANGE ENERGY (EX)                                 ==
    ! ==--------------------------------------------------------------==
    CALL dcopy(2*fpar%nnr1,dqg(1,1),1,v(1,1),1)
    cntl%tlsd=.TRUE.
    IF (ASSOCIATED(cp_xc)) cp_xc%polarized = .NOT. cp_xc%polarized
    CALL xcener(sxc2,vxc2,rhoval,rhoe(:,3:4),v)
    cntl%tlsd=.FALSE.
    sgcx2 = 0.0_real_8
    sgcc2 = 0.0_real_8
    vgc2  = 0.0_real_8
    IF (cntl%tgc) THEN
       ! CALCULATE THE GRADIENT OF THE DENSITY
       CALL fwfftn(v(:,1),.FALSE.,parai%allgrp)
       CALL zgthr(ncpw%nhg,v(:,1),vtemp(:,1),nzh)
       CALL fwfftn(v(:,2),.FALSE.,parai%allgrp)
       CALL zgthr(ncpw%nhg,v(:,2),vtemp(:,2),nzh)
       DO ir=1,fpar%nnr1
          rhoe(ir,3)=rhoe(ir,3)-rhoe(ir,4)
       ENDDO
       CALL graden(rhoe(:,3),v,grad(1,1),vtmp)
       CALL graden(rhoe(:,4),v,grad(1,5),vtmp)
       ! GRADIENT CORRECTION TO THE EXCHANGE ENERGY (EGCX)
       CALL gclsd(sgcx2,sgcc2,rhoe(:,3:),v,vtemp,rvtmp,grad,.FALSE.)
    ENDIF
    IF (ASSOCIATED(cp_xc)) cp_xc%polarized = .NOT. cp_xc%polarized
    ! V CONTAINS THE TOTAL POTENTIAL IN R-SPACE, MOVE IT TO RHOE
    !$omp parallel do private(IR)
    DO ir=1,fpar%nnr1
       rhoe(ir,3)=REAL(v(ir,1))
       rhoe(ir,4)=REAL(v(ir,2))
    ENDDO
    ! ==--------------------------------------------------------------==
    func1%mfxcx = mfxcx_1
    func1%mfxcc = mfxcc_1
    func1%mgcx  = mgcx_1
    func1%mgcc  = mgcc_1
    IF (func1%mgcx /= mgcx_is_skipped) THEN
       cntl%tgcx=.TRUE.
    ELSE
       cntl%tgcx=.FALSE.
    ENDIF
    IF (func1%mgcc /= mgcc_is_skipped) THEN
       cntl%tgcc=.TRUE.
    ELSE
       cntl%tgcc=.FALSE.
    ENDIF
    cntl%tgc=cntl%tgcc.OR.cntl%tgcx
    ! ==--------------------------------------------------------------==
    ! == ADD THE PIECES                                               ==
    ! ==--------------------------------------------------------------==
    sxc=sxc1+sxc2
    ener_com%vxc=vxc1+vxc2
    vgc=vgc1+vgc2
    sgcx=sgcx1+sgcx2
    sgcc=sgcc1+sgcc2
    ! ==--------------------------------------------------------------==
    ! == CORRECT ENERGY FOR CORE CHARGE                               ==
    ! ==--------------------------------------------------------------==
    IF (corel%tinlc) THEN
       sxc=sxc-corer%excco
       sgcx=sgcx-corer%egcxco
       sgcc=sgcc-corer%egccco
    ENDIF
    IF (corel%tinlc.AND.tfor) THEN
       ! ==--------------------------------------------------------------==
       ! == FORCE ON IONS DUE TO NLCC                                    ==
       ! ==--------------------------------------------------------------==
       CALL stopgm("VOFRHOR","DOES THIS WORK?",& 
            __LINE__,__FILE__)
       CALL  fwfftn(v(:,1),.FALSE.,parai%allgrp)
       CALL zgthr(ncpw%nhg,v,vtemp,nzh)
       CALL cofor(fion,vtemp)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! There are always two potentials that contribute to the total
    ! potential for each state 
    ! occ = 0.5 * ( V1 + V2 + V3 + V4 )
    ! alpha = V1 + V4
    ! beta  = V2 + V4
    ! -> add 0.5 * local potential to the XC parts, this sums always to 
    ! one total contribution (1-0.5)
    !$omp parallel do private(IR)
    DO ir=1,fpar%nnr1
       rhoe(ir,1)=rhoe(ir,1)-0.5_real_8*REAL(dqg(ir,1))
       rhoe(ir,2)=rhoe(ir,2)-0.5_real_8*REAL(dqg(ir,1))
       rhoe(ir,3)=rhoe(ir,3)-0.5_real_8*REAL(dqg(ir,1))
       rhoe(ir,4)=rhoe(ir,4)-0.5_real_8*REAL(dqg(ir,1))
    ENDDO
    ! Rearange the potentials
    !$omp parallel do private(IR,VP1,VP2,VP3)
    DO ir=1,fpar%nnr1
       vp1 = 0.5_real_8*(rhoe(ir,1)+rhoe(ir,2)+rhoe(ir,3)+rhoe(ir,4))
       vp2 = rhoe(ir,1)+rhoe(ir,4)
       vp3 = rhoe(ir,2)+rhoe(ir,4)
       rhoe(ir,1) = vp1
       rhoe(ir,2) = vp2
       rhoe(ir,3) = vp3
    ENDDO
    ! ==--------------------------------------------------------------==
    ! == EXCHANGE POTENTIAL FROM A-B STATES                           ==
    ! ==--------------------------------------------------------------==
    ! TRANSFORM THE DENSITY TO G SPACE
    CALL zeroing(v(:,1))!,maxfft)
    !$omp parallel do private(IR)
    DO ir=1,fpar%nnr1
       v(ir,1)=CMPLX(rhoe(ir,5),0._real_8,kind=real_8)
    ENDDO
    CALL  fwfftn(v(:,1),.FALSE.,parai%allgrp)
    ek=0._real_8
    IF (isos1%tclust) THEN
       IF (isos3%ps_type.EQ.1) THEN
          CALL stopgm('VOFRHOR','HOCKNEY PS not impl.',& 
               __LINE__,__FILE__)
       ELSE
          DO ig=1,ncpw%nhg
             vtemp(ig,1) = 2._real_8*v(nzh(ig),1) * scg(ig)
             ek = ek + 2._real_8 * REAL(CONJG(vtemp(ig,1))*v(nzh(ig),1))
          ENDDO
       ENDIF
    ELSE
       ig1=1
       IF (geq0) ig1=2
       DO ig=ig1,ncpw%nhg
          g2=parm%tpiba2*hg(ig)
          scgy=CMPLX(fpi/g2,0._real_8,kind=real_8)
          vtemp(ig,1) = 2._real_8*v(nzh(ig),1) * scgy
          ek = ek + 2._real_8 * REAL(CONJG(vtemp(ig,1))*v(nzh(ig),1))
       ENDDO
       IF (geq0) vtemp(1,1)=0._real_8
    ENDIF
    CALL zeroing(v(:,1))!,maxfft)
    DO ig=1,ncpw%nhg
       v(indz(ig),1) = CONJG(vtemp(ig,1))
       v(nzh(ig),1)  = vtemp(ig,1)
    ENDDO
    IF (geq0) THEN
       v(nzh(1),1) = vtemp(1,1)
    ENDIF
    CALL  invfftn(v(:,1),.FALSE.,parai%allgrp)
    !$omp parallel do private(IR)
    DO ir=1,fpar%nnr1
       rhoe(ir,4)=REAL(v(ir,1))
    ENDDO
    ! ==--------------------------------------------------------------==
    IF (lspin2%tpenal) THEN
       !$omp parallel do private(IR)
       DO ir=1,fpar%nnr1
          v(ir,1)=CMPLX(rhoe(ir,5),0._real_8,kind=real_8)
       ENDDO
       CALL  fwfftn(v(:,1),.FALSE.,parai%allgrp)
       CALL zgthr(ncpw%nhg,v,vtemp,nzh)
       IF (tfor) THEN
          CALL potabfor(fion,vtemp)
       ENDIF
       CALL potab(epotab,vtemp)
       ener_c%eht_ab=REAL(epotab)*parm%omega
       CALL zeroing(v(:,1))!,maxfft)
       DO ig=1,ncpw%nhg
          v(indz(ig),1) = CONJG(vtemp(ig,1))
          v(nzh(ig),1)  = vtemp(ig,1)
       ENDDO
       IF (geq0) THEN
          v(nzh(1),1) = vtemp(1,1)
       ENDIF
       CALL  invfftn(v(:,1),.FALSE.,parai%allgrp)
       !$omp parallel do private(IR)
       DO ir=1,fpar%nnr1
          rhoe(ir,5)=REAL(v(ir,1))
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == SINGLE PIECES OF THE ENERGY:                                 ==
    ! ==--------------------------------------------------------------==
    nnrs = spar%nr1s*spar%nr2s*spar%nr3s
    ener_com%egc=(sgcc+sgcx)*parm%omega/REAL(nnrs,kind=real_8)
    ener_com%exc=sxc*parm%omega/REAL(nnrs,kind=real_8) + ener_com%egc + ek*parm%omega
    ! 
    !
    DEALLOCATE(rvtmp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    !
    CALL tihalt('   VOFRHOR',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE vofrhor
  ! ==================================================================
  SUBROUTINE give_scr_vofrhor(lvofrhor,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lvofrhor
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: ltgc, ltinlc

    IF (cntl%tgc) THEN
       ltgc=fpar%nnr1*clsd%nlsd*4      ! GRADEN
    ELSE
       ltgc=0
    ENDIF
    IF (corel%tinlc) THEN            ! XCENER
       ltinlc=fpar%nnr1*clsd%nlsd
    ELSE
       ltinlc=0
    ENDIF
    lvofrhor = MAX(fpar%nnr1,ncpw%nhg*2)*clsd%nlsx+MAX(ltgc,ltinlc) ! GCENER and COREC
    IF (pslo_com%tivan) lvofrhor=MAX(lvofrhor,ncpw%nhg*4+maxsys%nax) ! NEWD
    tag='MAX(NNR1,NHG*2)*NLSX+.........'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_vofrhor
  ! ==================================================================

END MODULE vofrhos_utils
