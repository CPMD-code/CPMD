MODULE lr_pcg_utils
  USE cppt,                            ONLY: nzh,&
                                             scg
  USE dotp_utils,                      ONLY: dotp
  USE elct,                            ONLY: crge
  USE fftmain_utils,                   ONLY: fwfftn
  USE hpsi_utils,                      ONLY: give_scr_hpsi,&
                                             hpsi
  USE kinds,                           ONLY: real_8
  USE linres,                          ONLY: lr02,&
                                             lr03
  USE lr_ortho_utils,                  ONLY: give_scr_lr_ortho,&
                                             lr_ortho
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai
  USE poin,                            ONLY: potr
  USE rho1ofr_utils,                   ONLY: rho1ofr
  USE spin,                            ONLY: clsd
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

  PUBLIC :: lr_pcg
  PUBLIC :: give_scr_lr_pcg
  PUBLIC :: quadsr
  PUBLIC :: give_scr_quadsr
  !public :: edofrhod

CONTAINS

  ! ==================================================================
  SUBROUTINE lr_pcg(c1,c2,vpp,hnm1,c0,drhoe,ddxc,eigv,sc0,psi,&
       nstate,dinit)
    ! ==--------------------------------------------------------------==
    ! ==  PRECONDITIONED CONJUGATE GRADIENT OPTIMIZATION              ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: vpp(:), drhoe(:,:), ddxc(*), &
                                                eigv(*)
    COMPLEX(real_8)                          :: sc0(*), psi(:,:)
    INTEGER                                  :: nstate
    COMPLEX(real_8) :: c0(ncpw%ngw,nstate), hnm1(ncpw%ngw,nstate), &
      c2(ncpw%ngw,nstate), c1(ncpw%ngw,nstate)
    LOGICAL                                  :: dinit

    INTEGER                                  :: i, ig, isub
    LOGICAL                                  :: qsearch
    REAL(real_8)                             :: a0, gamma, ggnorm
    REAL(real_8), DIMENSION(2), SAVE         :: ghist

    CALL tiset('    LR_PCG',isub)
    ! ==--------------------------------------------------------------==
    ! ..preconditioning
    DO i=1,nstate
       DO ig=1,ncpw%ngw
          c2(ig,i)=vpp(ig)*c2(ig,i)
       ENDDO
    ENDDO
    ! ..norm
    ggnorm=0.0_real_8
    DO i=1,nstate
       ggnorm=ggnorm+dotp(ncpw%ngw,c2(:,i),c2(:,i))
    ENDDO
    CALL mp_sum(ggnorm,parai%allgrp)
    qsearch=.TRUE.
    IF (ggnorm.LT.lr02%tol_qs.AND..NOT.lr03%txc_analytic) qsearch=.FALSE.
    IF (dinit) THEN
       ! ==--------------------------------------------------------------==
       ! ==  INITIAL CALL                                                ==
       ! ==--------------------------------------------------------------==
       ghist(1)=ggnorm
       CALL dcopy(2*ncpw%ngw*nstate,c2,1,hnm1,1)
       CALL lr_ortho(nstate,c0,hnm1)
       IF (qsearch) THEN
          CALL quadsr(a0,c0,vpp,hnm1,c2,drhoe,ddxc,eigv,sc0,psi(:,1),&
               nstate)
          CALL daxpy(2*ncpw%ngw*nstate,a0,hnm1,1,c1,1)
       ELSE
          CALL daxpy(2*ncpw%ngw*nstate,lr02%dlrs,hnm1,1,c1,1)
       ENDIF
    ELSE
       ! ==--------------------------------------------------------------==
       ! ==  CONJUGATE GRADIENT                                          ==
       ! ==--------------------------------------------------------------==
       ghist(2)=ggnorm
       gamma=ghist(2)/ghist(1)
       CALL dscal(2*ncpw%ngw*nstate,gamma,hnm1,1)
       CALL daxpy(2*ncpw%ngw*nstate,1._real_8,c2,1,hnm1,1)
       CALL lr_ortho(nstate,c0,hnm1)
       IF (qsearch)THEN
          CALL quadsr(a0,c0,vpp,hnm1,c2,drhoe,ddxc,eigv,sc0,psi(:,1),&
               nstate)
          CALL daxpy(2*ncpw%ngw*nstate,a0,hnm1,1,c1,1)
       ELSE
          CALL daxpy(2*ncpw%ngw*nstate,lr02%dlrs,hnm1,1,c1,1)
       ENDIF
       ghist(1)=ghist(2)
    ENDIF
    dinit=.FALSE.
    CALL tihalt('    LR_PCG',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE lr_pcg
  ! ==================================================================
  SUBROUTINE give_scr_lr_pcg(llr_pcg,tag,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: llr_pcg
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate

    INTEGER                                  :: llr_ortho, lquadsr

! ==--------------------------------------------------------------==

    CALL give_scr_quadsr(lquadsr,tag,nstate)
    CALL give_scr_lr_ortho(llr_ortho,nstate)
    llr_pcg=MAX(lquadsr,llr_ortho)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_lr_pcg
  ! ==================================================================
  SUBROUTINE quadsr(a0,c0,vpp,dc,c2,drhoe,ddxc,eigv,sc0,psi,&
       nstate)
    ! ==--------------------------------------------------------------==
    ! ==  Quadratic line search                                       ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: a0, vpp(:), drhoe(:,:), &
                                                ddxc(*), eigv(*)
    COMPLEX(real_8)                          :: sc0(*), psi(:)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c2(ncpw%ngw,nstate), &
                                                dc(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)

    INTEGER                                  :: i, ig
    REAL(real_8)                             :: eq(2), zlin, zquad

! ==--------------------------------------------------------------==
! ..de_preconditioning

    DO i=1,nstate
       DO ig=1,ncpw%ngw
          c2(ig,i)=c2(ig,i)/vpp(ig)
       ENDDO
    ENDDO
    ! ..linear terms
    zlin=0._real_8
    DO i=1,nstate
       zlin=zlin-dotp(ncpw%ngw,dc(:,i),c2(:,i))
    ENDDO
    ! ..delta force
    CALL hpsi(dc,c2,sc0,potr,psi,nstate,1,clsd%nlsd)
    DO i=1,nstate
       CALL dscal(2*ncpw%ngw,crge%f(i,1),c2(1,i),1)
    ENDDO
    DO i=1,nstate
       CALL daxpy(2*ncpw%ngw,eigv(i),dc(1,i),1,c2(1,i),1)
    ENDDO
    ! ..XC and Hartree kernel
    CALL rho1ofr(c0,dc,crge%f(:,1),drhoe,psi,nstate)
    CALL edofrhod(eq,drhoe,ddxc,psi)
    ! ..quadratic terms
    zquad=-eq(1)-eq(2)
    DO i=1,nstate
       zquad=zquad+dotp(ncpw%ngw,dc(:,i),c2(:,i))
    ENDDO
    CALL mp_sum(zquad,parai%allgrp)
    CALL mp_sum(zlin,parai%allgrp)
    IF (ABS(zquad).LT.1.e-12_real_8) THEN
       a0=lr02%dlrs
    ELSE
       a0=zlin/zquad
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE quadsr
  ! ==================================================================
  SUBROUTINE give_scr_quadsr(lquadsr,tag,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lquadsr
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate

    INTEGER                                  :: lhpsi

! ==--------------------------------------------------------------==

    CALL give_scr_hpsi(lhpsi,tag,nstate)
    lquadsr=lhpsi
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_quadsr
  ! ==================================================================
  SUBROUTINE edofrhod(e2,drhoe,ddxc,psi)
    ! ==--------------------------------------------------------------==
    ! calculate the quadratic contribution of the XC and Hartree Kernel
    ! ==--------------------------------------------------------------==
    ! 
    REAL(real_8)                             :: e2(2), drhoe(:,:), &
                                                ddxc(fpar%nnr1,clsd%nlsd)
    COMPLEX(real_8)                          :: psi(:)

    COMPLEX(real_8)                          :: v1
    INTEGER                                  :: ig, ir
    REAL(real_8)                             :: eht2, exc2

! 
! ==--------------------------------------------------------------==
! ..transform dn1 to g-space

    CALL zeroing(psi)!,maxfft)
    !$omp parallel do private(IR)
    DO ir=1,fpar%nnr1
       psi(ir)=CMPLX(drhoe(ir,1),0._real_8,kind=real_8)
    ENDDO
    CALL  fwfftn(psi,.FALSE.,parai%allgrp)
    ! ==--------------------------------------------------------------==
    ! == This is the Hartree Potential                                ==
    ! ==--------------------------------------------------------------==
    eht2=0._real_8
    DO ig=1,ncpw%nhg
       v1 = psi(nzh(ig)) * scg(ig)
       eht2 = eht2 + REAL(CONJG(v1)*psi(nzh(ig)))
    ENDDO
    eht2 = eht2*parm%omega
    ! ==--------------------------------------------------------------==
    ! == XC contributions (we make a LDA approximation here)          ==
    ! == DDXC holds at least a reasonable approximation to W(d)       ==
    ! ==--------------------------------------------------------------==
    exc2=0._real_8
    DO ir=1,fpar%nnr1
       exc2=exc2+ddxc(ir,1)*drhoe(ir,1)*drhoe(ir,1)
    ENDDO
    IF (cntl%tlsd) THEN
       DO ir=1,fpar%nnr1
          exc2=exc2+ddxc(ir,2)*drhoe(ir,2)*drhoe(ir,2)
       ENDDO
    ENDIF
    exc2=0.5_real_8*exc2*parm%omega/REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)
    ! ==--------------------------------------------------------------==
    e2(1) = eht2
    e2(2) = exc2
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE edofrhod
  ! ==================================================================

END MODULE lr_pcg_utils
