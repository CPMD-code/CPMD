MODULE vpsi_lse_utils
  USE cnst,                            ONLY: uimag
  USE cppt,                            ONLY: hg,&
                                             indzs,&
                                             nzhs
  USE dotp_utils,                      ONLY: dotp
  USE ener,                            ONLY: ener_c,&
                                             ener_d
  USE error_handling,                  ONLY: stopgm
  USE fft_maxfft,                      ONLY: maxfft
  USE fftmain_utils,                   ONLY: fwfftn,&
                                             invfftn
  USE geq0mod,                         ONLY: geq0
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai
  USE spin,                            ONLY: clsd,&
                                             lspin2
  USE system,                          ONLY: fpar,&
                                             group,&
                                             ncpw,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: vpsi_lse

CONTAINS

  ! ==================================================================
  SUBROUTINE vpsi_lse(c0,c2,f,vpot,psi,nstate,tekin)
    ! ==================================================================
    COMPLEX(real_8)                          :: c0(ncpw%ngw,*), c2(ncpw%ngw,*)
    REAL(real_8)                             :: f(*), vpot(fpar%nnr1,*)
    COMPLEX(real_8)                          :: psi(:)
    INTEGER                                  :: nstate
    LOGICAL                                  :: tekin

    CHARACTER(*), PARAMETER                  :: procedureN = 'vpsi_lse'

    COMPLEX(real_8)                          :: fm, fp
    COMPLEX(real_8), ALLOCATABLE             :: psiab(:), vaa(:), vbb(:)
    INTEGER                                  :: i, ierr, ig, ir, isub
    REAL(real_8)                             :: ca, ff, fh, fia, fib, g2, sa, &
                                                tab
    REAL(real_8), ALLOCATABLE                :: kvab(:,:)

    CALL tiset('  VPSI_LSE',isub)
    IF (group%nogrp.GT.1) CALL stopgm("VPSI","GROUP CODE MISSING",& 
         __LINE__,__FILE__)
    ALLOCATE(psiab(maxfft),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL zeroing(psiab)!,maxfft)
#ifdef __SR8000
    !poption parallel
#endif
#ifndef __SR11000
    !$omp parallel do private(IG)
#endif
    DO ig=1,ncpw%ngw
       psiab(nzhs(ig))=c0(ig,clsd%ialpha)+uimag*c0(ig,clsd%ibeta)
       psiab(indzs(ig))=CONJG(c0(ig,clsd%ialpha))+&
            uimag*CONJG(c0(ig,clsd%ibeta))
    ENDDO
    IF (geq0) psiab(nzhs(1))=c0(1,clsd%ialpha)+uimag*c0(1,clsd%ibeta)
    CALL  invfftn(psiab,.TRUE.,parai%allgrp)
    IF (lspin2%troot) THEN
#ifdef __SR8000
       !poption parallel
#endif
       !$omp parallel do private(IR)
       DO ir=1,fpar%nnr1
          psi(ir)=vpot(ir,1)*psiab(ir)
       ENDDO
    ELSE
#ifdef __SR8000
       !poption parallel
#endif
       !$omp parallel do private(IR)
       DO ir=1,fpar%nnr1
          psi(ir)=vpot(ir,2)*REAL(psiab(ir))&
               +uimag*vpot(ir,3)*AIMAG(psiab(ir))
       ENDDO
    ENDIF
    CALL  fwfftn(psi,.TRUE.,parai%allgrp)
    fh=0.5_real_8
    IF (lspin2%tlsets) THEN
       fia=f(clsd%ialpha)
       fib=f(clsd%ibeta)
       IF (tekin) THEN
#ifdef __SR8000
          !poption parallel
          !voption indep(NZHS,INDZS,C2)
#endif
          DO ig=1,ncpw%ngw
             ff=-fh*parm%tpiba2*hg(ig)
             fp=psi(nzhs(ig))+psi(indzs(ig))
             fm=psi(nzhs(ig))-psi(indzs(ig))
             c2(ig,clsd%ialpha)=fia*(ff*c0(ig,clsd%ialpha)&
                  -fh*CMPLX(REAL(fp),AIMAG(fm),kind=real_8))
             c2(ig,clsd%ibeta) =fib*(ff*c0(ig,clsd%ibeta)&
                  -fh*CMPLX(AIMAG(fp),-REAL(fm),kind=real_8))
          ENDDO
       ELSE
#ifdef __SR8000
          !poption parallel
          !voption indep(NZHS,INDZS,C2)
#endif
          DO ig=1,ncpw%ngw
             fp=psi(nzhs(ig))+psi(indzs(ig))
             fm=psi(nzhs(ig))-psi(indzs(ig))
             c2(ig,clsd%ialpha)=c2(ig,clsd%ialpha)&
                  -fia*fh*CMPLX(REAL(fp),AIMAG(fm),kind=real_8)
             c2(ig,clsd%ibeta) =c2(ig,clsd%ibeta)&
                  -fib*fh*CMPLX(AIMAG(fp),-REAL(fm),kind=real_8)
          ENDDO
       ENDIF
    ELSE
       IF (tekin) THEN
#ifdef __SR8000
          !poption parallel
          !voption indep(NZHS,INDZS,C2)
#endif
          DO ig=1,ncpw%ngw
             ff=-fh*parm%tpiba2*hg(ig)
             fp=psi(nzhs(ig))+psi(indzs(ig))
             fm=psi(nzhs(ig))-psi(indzs(ig))
             c2(ig,clsd%ialpha)=ff*c0(ig,clsd%ialpha)&
                  -fh*CMPLX(REAL(fp),AIMAG(fm),kind=real_8)
             c2(ig,clsd%ibeta) =ff*c0(ig,clsd%ibeta)&
                  -fh*CMPLX(AIMAG(fp),-REAL(fm),kind=real_8)
          ENDDO
       ELSE
#ifdef __SR8000
          !poption parallel
          !voption indep(NZHS,INDZS,C2)
#endif
          DO ig=1,ncpw%ngw
             fp=psi(nzhs(ig))+psi(indzs(ig))
             fm=psi(nzhs(ig))-psi(indzs(ig))
             c2(ig,clsd%ialpha)=c2(ig,clsd%ialpha)&
                  -fh*CMPLX(REAL(fp),AIMAG(fm),kind=real_8)
             c2(ig,clsd%ibeta) =c2(ig,clsd%ibeta)&
                  -fh*CMPLX(AIMAG(fp),-REAL(fm),kind=real_8)
          ENDDO
       ENDIF
    ENDIF
    IF (lspin2%troot) THEN
       ALLOCATE(vaa(ncpw%ngw),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(vbb(ncpw%ngw),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(kvab(2,nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
#ifdef __SR8000
       !poption parallel
#endif
       !$omp parallel do private(IR)
       DO ir=1,fpar%nnr1
          psi(ir)=vpot(ir,2)*REAL(psiab(ir))&
               +uimag*vpot(ir,3)*AIMAG(psiab(ir))
       ENDDO
       CALL  fwfftn(psi,.TRUE.,parai%allgrp)
#ifdef __SR8000
       !poption parallel
       !voption indep(NZHS,INDZS,VAA,VBB)
#endif
       !$omp parallel do private(IG,FP,FM)
       DO ig=1,ncpw%ngw
          fp=-0.5_real_8*(psi(nzhs(ig))+psi(indzs(ig)))
          fm=-0.5_real_8*(psi(nzhs(ig))-psi(indzs(ig)))
          vaa(ig) = CMPLX(REAL(fp),AIMAG(fm),kind=real_8)
          vbb(ig) = CMPLX(AIMAG(fp),-REAL(fm),kind=real_8)
       ENDDO
       CALL daxpy(2*ncpw%ngw,-1.0_real_8,vaa,1,c2(1,clsd%ialpha),1)
       CALL daxpy(2*ncpw%ngw,-1.0_real_8,vbb,1,c2(1,clsd%ibeta),1)
       DO i=1,nstate
          kvab(1,i)=dotp(ncpw%ngw,vaa,c0(:,i))
          kvab(2,i)=dotp(ncpw%ngw,vbb,c0(:,i))
       ENDDO
       CALL mp_sum(kvab,2*nstate,parai%allgrp)
       DO i=1,clsd%ialpha-1
          tab=f(i)*kvab(1,i)
          CALL daxpy(2*ncpw%ngw,tab,c0(1,clsd%ialpha),1,c2(1,i),1)
          tab=f(i)*kvab(2,i)
          CALL daxpy(2*ncpw%ngw,tab,c0(1,clsd%ibeta),1,c2(1,i),1)
          tab=2._real_8*kvab(1,i)
          CALL daxpy(2*ncpw%ngw,tab,c0(1,i),1,c2(1,clsd%ialpha),1)
          tab=2._real_8*kvab(2,i)
          CALL daxpy(2*ncpw%ngw,tab,c0(1,i),1,c2(1,clsd%ibeta),1)
       ENDDO
       tab=0.5_real_8*(kvab(1,clsd%ibeta)-kvab(2,clsd%ialpha))
       CALL daxpy(2*ncpw%ngw,tab,c0(1,clsd%ibeta),1,c2(1,clsd%ialpha),1)
       tab=0.5_real_8*(kvab(2,clsd%ialpha)-kvab(1,clsd%ibeta))
       CALL daxpy(2*ncpw%ngw,tab,c0(1,clsd%ialpha),1,c2(1,clsd%ibeta),1)
       DEALLOCATE(vaa,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(vbb,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(kvab,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ELSEIF (lspin2%tross) THEN
       ! Special code for exchange contribution of singly occupied states
       ! for the open shell (biradical) singlet method
#ifdef __SR8000
       !poption parallel
#endif
       !$omp parallel do private(IR)
       DO ir=1,fpar%nnr1
          psi(ir)=psiab(ir)*vpot(ir,4)
       ENDDO
       CALL  fwfftn(psi,.TRUE.,parai%allgrp)
       ! In the following loop states are interchanged
       ! a contributes to b and b to a
#ifdef __SR8000
       !poption parallel
#endif
       DO ig=1,ncpw%ngw
          fp=0.5_real_8*(psi(nzhs(ig))+psi(indzs(ig)))
          fm=0.5_real_8*(psi(nzhs(ig))-psi(indzs(ig)))
          c2(ig,clsd%ialpha)= c2(ig,clsd%ialpha)-CMPLX(AIMAG(fp),-REAL(fm),kind=real_8)
          c2(ig,clsd%ibeta) = c2(ig,clsd%ibeta) -CMPLX(REAL(fp),AIMAG(fm),kind=real_8)
       ENDDO
    ELSEIF (lspin2%tcas22) THEN
       ! Special code for the contribution to the forces from the
       ! Kinetic energy in the CAS22 method
       ca=COS(ener_d%casang)
       sa=SIN(ener_d%casang)
       DO ig=1,ncpw%ngw
          g2=0.5_real_8*parm%tpiba2*hg(ig)
          c2(ig,clsd%ialpha)=c2(ig,clsd%ialpha)+&
               g2*(sa*ca*c0(ig,clsd%ibeta)-sa*sa*c0(ig,clsd%ialpha))
          c2(ig,clsd%ibeta) =c2(ig,clsd%ibeta) +&
               g2*(sa*ca*c0(ig,clsd%ialpha)+sa*sa*c0(ig,clsd%ibeta))
       ENDDO
       ! now the contribution |a> = V*|b>
       !CDIR NODEP
#ifdef __SR8000
       !poption parallel
#endif
       DO ig=1,ncpw%ngw
          psi(nzhs(ig))=c0(ig,clsd%ialpha)+uimag*c0(ig,clsd%ibeta)
          psi(indzs(ig))=CONJG(c0(ig,clsd%ialpha))+&
               uimag*CONJG(c0(ig,clsd%ibeta))
       ENDDO
       IF (geq0) psi(nzhs(1))=c0(1,clsd%ialpha)+uimag*c0(1,clsd%ibeta)
       CALL  invfftn(psi,.TRUE.,parai%allgrp)
#ifdef __SR8000
       !voption indep(PSI,VPOT)
#endif
       !$omp parallel do private(IR)
       DO ir=1,fpar%nnr1
          psi(ir)=psi(ir)*vpot(ir,4)
       ENDDO
       CALL  fwfftn(psi,.TRUE.,parai%allgrp)
       ! In the following loop states are interchanged
       ! a contributes to b and b to a
#ifdef __SR8000
       !poption parallel
#endif
       DO ig=1,ncpw%ngw
          fp=0.5_real_8*(psi(nzhs(ig))+psi(indzs(ig)))
          fm=0.5_real_8*(psi(nzhs(ig))-psi(indzs(ig)))
          c2(ig,clsd%ialpha)= c2(ig,clsd%ialpha)-CMPLX(AIMAG(fp),-REAL(fm),kind=real_8)
          c2(ig,clsd%ibeta) = c2(ig,clsd%ibeta) -CMPLX(REAL(fp),AIMAG(fm),kind=real_8)
       ENDDO
    ENDIF
    IF (lspin2%tpenal) THEN
       ! Contribution from Penalty function
       !$omp parallel do private(IR)
       DO ir=1,fpar%nnr1
          psi(ir)=psiab(ir)*vpot(ir,5)
       ENDDO
       CALL  fwfftn(psi,.TRUE.,parai%allgrp)
       ! In the following loop states are interchanged
       ! a contributes to b and b to a
#ifdef __SR8000
       !poption parallel
#endif
       DO ig=1,ncpw%ngw
          g2=0.25_real_8*parm%tpiba2*hg(ig)
          fp=0.5_real_8*(psi(nzhs(ig))+psi(indzs(ig)))
          fm=0.5_real_8*(psi(nzhs(ig))-psi(indzs(ig)))
          c2(ig,clsd%ialpha)= c2(ig,clsd%ialpha)-&
               2._real_8*ener_c%etot_ab*(g2*c0(ig,clsd%ibeta)+CMPLX(AIMAG(fp),-REAL(fm),kind=real_8))
          c2(ig,clsd%ibeta) = c2(ig,clsd%ibeta) -&
               2._real_8*ener_c%etot_ab*(g2*c0(ig,clsd%ialpha)+CMPLX(REAL(fp),AIMAG(fm),kind=real_8))
       ENDDO
    ENDIF
    ! 
    DEALLOCATE(psiab,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! 
    CALL tihalt('  VPSI_LSE',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE vpsi_lse
  ! ==================================================================

END MODULE vpsi_lse_utils
