MODULE elf_utils
  !!use densto_utils, only : densto
  USE bsym,                            ONLY: bsclcs
  USE cnst,                            ONLY: pi,&
                                             uimag
  USE cppt,                            ONLY: gk,&
                                             indz,&
                                             indzs,&
                                             nzh,&
                                             nzhs
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE fft_maxfft,                      ONLY: maxfft
  USE fftmain_utils,                   ONLY: fwfftn,&
                                             invfftn
  USE fftnew_utils,                    ONLY: setfftn
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE meta_multiple_walkers_utils,     ONLY: mw_filename
  USE mp_interface,                    ONLY: mp_max
  USE mw,                              ONLY: mwi,&
                                             tmw
  USE parac,                           ONLY: parai,&
                                             paral
  USE phfac_utils,                     ONLY: phfac
  USE pimd,                            ONLY: ipcurr
  USE prden,                           ONLY: elfcb
  USE pslo,                            ONLY: pslo_com
  USE readsr_utils,                    ONLY: xstring
  USE rhoofr_c_utils,                  ONLY: rhoofr_c
  USE rhoofr_utils,                    ONLY: give_scr_rhoofr,&
                                             rhoofr
  USE rnlsm_utils,                     ONLY: give_scr_rnlsm,&
                                             rnlsm
  USE spin,                            ONLY: clsd
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             ncpw,&
                                             nkpt,&
                                             parm
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: elf
  PUBLIC :: give_scr_elf

CONTAINS

  ! ==================================================================
  SUBROUTINE elf(c0,tau0,rhoe,v,nstate,nkpoint)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==  THE NORMALIZED ELECTRON DENSITY RHOE IN REAL SPACE          ==
    ! ==     AND WRITES IT TO THE FILE rho.movie                      ==
    ! ==  THE ELECTRON LOCALIZATION FUNCTION:                         ==
    ! ==     ELF=1/(1+( (D+ELFEPS) /D_h)^2)                           ==
    ! ==     D  =(1/2) sum_i |nabla psi_i|^2 -                        ==
    ! ==         (1/8) {|nabla rho|^2} / rho                          ==
    ! ==     D_h=(3/10) (3 pi^2)^2/3 rho^5/3                          ==
    ! ==     ELFEPS=2.87e-5_real_8 (-> set ELF=0.5 for D=0 and the van der   ==
    ! ==     Waals density rho=0.001 to avoid artifacts for small     ==
    ! ==     densities rho)                                           ==
    ! ==     (ELF is set to zero for ELF<ELFCUT)                      ==
    ! ==                                                              ==
    ! ==     psi_i is the KS-spin(!)orbital, rho is the total         ==
    ! ==     density, and i=1, ..., number of electrons               ==
    ! ==     (this definition also applies to spin polarized systems  ==
    ! ==     and yields a kind of averaged ELF value)                 ==
    ! ==     see                                                      ==
    ! ==     A. D. Becke and K. E. Edgecombe, JCP 92 (1990) 5397      ==
    ! ==     B. Silvi and A. Savin, Nature 371 (1994) 683             ==
    ! ==     (note that the exponent of 3pi^2 is wrong in this paper) ==
    ! ==     A. Savin, private communication (for epsilon)            ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:), &
                                                rhoe(fpar%nnr1,clsd%nlsd)
    COMPLEX(real_8)                          :: v(maxfft)
    INTEGER                                  :: nstate, nkpoint
    COMPLEX(real_8)                          :: c0(nkpt%ngwk,nstate,nkpoint)

    CHARACTER(*), PARAMETER                  :: procedureN = 'elf'

    CHARACTER(len=12)                        :: cflbod, cipnum
    CHARACTER(len=15)                        :: filen
    COMPLEX(real_8), ALLOCATABLE             :: vtemp(:)
    INTEGER                                  :: i, i1, i2, ierr, ig, ikind, &
                                                ir, k, n1, n2, space
    REAL(real_8)                             :: aelf, coefe, coefs, elfmax, &
                                                elfmin, rhom
    REAL(real_8), ALLOCATABLE                :: elff(:)

    CALL setfftn(0)
    IF (paral%io_parent) WRITE(6,'(/,A)') ' CALCULATE ELF FUNCTION '
    ! ==--------------------------------------------------------------==
    ALLOCATE(elff(fpar%nnr1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    ! TODO check stat 
    ALLOCATE(vtemp(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    CALL phfac(tau0)
    DO ikind=1,nkpoint
       CALL rnlsm(c0(:,:,ikind),nstate,1,ikind,.FALSE.)
    ENDDO
    IF (tkpts%tkpnt) THEN
       CALL rhoofr_c(c0,rhoe,v,nstate)
    ELSE
       CALL rhoofr  (c0(:,:,1),rhoe,v,nstate)
    ENDIF
    ! ..VTEMP=density in G-space
    !$omp parallel do private(IR) schedule(static)
    DO ir=1,fpar%nnr1
       v(ir)=CMPLX(rhoe(ir,1),0._real_8,kind=real_8)
    ENDDO
    CALL  fwfftn(v,.FALSE.,parai%allgrp)
    !$omp parallel do private(IG) schedule(static)
    DO ig=1,ncpw%nhg
       vtemp(ig) = v(nzh(ig))
    ENDDO
    ! ..calculate |nabla*rho|^2
    CALL zeroing(v)!,maxfft)
    !CDIR NODEP
    DO ig=1,ncpw%nhg
       v(nzh(ig))=parm%tpiba*uimag*gk(1,ig)*vtemp(ig)
       v(indz(ig))=-parm%tpiba*uimag*gk(1,ig)*CONJG(vtemp(ig))
    ENDDO
    CALL  invfftn(v,.FALSE.,parai%allgrp)
    !$omp parallel do private(IR) schedule(static)
    DO ir=1,fpar%nnr1
       elff(ir)=REAL(v(ir))*REAL(v(ir))
    ENDDO
    CALL zeroing(v)!,maxfft)
    !CDIR NODEP
    DO ig=1,ncpw%nhg
       v(nzh(ig))=parm%tpiba*(uimag*gk(2,ig)-gk(3,ig))*vtemp(ig)
       v(indz(ig))=parm%tpiba*(-uimag*gk(2,ig)+gk(3,ig))*CONJG(vtemp(ig))
    ENDDO
    CALL  invfftn(v,.FALSE.,parai%allgrp)
    !$omp parallel do private(IR) schedule(static)
    DO ir=1,fpar%nnr1
       elff(ir)=elff(ir)+REAL(v(ir)*CONJG(v(ir)))
    ENDDO
    ! ..calculate -{|nabla rho|^2 / 8*rho} / D_h
    coefe=pi**2
    coefe=(3._real_8*coefe)**(2._real_8/3._real_8)
    coefe=0.3_real_8*coefe
    DO ir=1,fpar%nnr1
       rhom=MAX(rhoe(ir,1),0._real_8)
       IF (rhom.GT.elfcb%elfcut) THEN
          elff(ir)=-elff(ir)/(8._real_8*coefe*rhoe(ir,1)**2.6666666667_real_8)
       ELSE
          elff(ir)=0._real_8
       ENDIF
    ENDDO
    ! ..calculate 0.5*SUM_i{|nabla*phi_i|^2 / D_h }
    DO ikind=1,nkpoint
       DO i=1,nstate
          coefs=0.5_real_8*crge%f(i,ikind)/parm%omega
          DO k=1,3
             CALL zeroing(v)!,maxfft)
             IF (tkpts%tkpnt) THEN
                CALL stopgm('elf','not programmed for k-points',& 
                     __LINE__,__FILE__)
                !CDIR NODEP
                DO ig=1,ncpw%ngw
                   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                   ! here we have to add GK for the K-POINTS
                   ! does this also need symmetrization?????
                   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                   v(nzhs(ig))=parm%tpiba*uimag*gk(k,ig)*c0(ig,i,ikind)
                   v(indzs(ig))=CONJG(parm%tpiba*uimag*gk(k,ig)*c0(ig,i,ikind))
                ENDDO
             ELSE
                !CDIR NODEP
                DO ig=1,ncpw%ngw
                   v(nzhs(ig))=parm%tpiba*uimag*gk(k,ig)*c0(ig,i,ikind)
                   v(indzs(ig))=CONJG(parm%tpiba*uimag*gk(k,ig)*c0(ig,i,ikind))
                ENDDO
             ENDIF
             CALL  invfftn(v,.TRUE.,parai%allgrp)
             ! ..here to add Vanderbilt part of nabla*psi
             IF (pslo_com%tivan) THEN
                CALL stopgm('ELF','VANDERBILT NOT IMPLEMENTED',& 
                     __LINE__,__FILE__)
             ENDIF
             DO ir=1,fpar%nnr1
                rhom=MAX(rhoe(ir,1),0._real_8)
                IF (rhom.GT.elfcb%elfcut) THEN
                   elff(ir)=elff(ir)+coefs*REAL(v(ir))**2/&
                        (coefe*rhoe(ir,1)**1.6666666667_real_8)
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    ! ..calculate 1/(1+(D/D_h)^2)
    DO ir=1,fpar%nnr1
       rhom=MAX(rhoe(ir,1),0._real_8)
       IF (rhom.GT.elfcb%elfcut) THEN
          aelf=elff(ir)+elfcb%elfeps/(coefe*rhoe(ir,1)**1.6666666667_real_8)
          elff(ir)=1._real_8/(1._real_8+aelf**2)
       ENDIF
    ENDDO
    filen='ELF'
    ! ..VTEMP=ELF in G-space
    !$omp parallel do private(IR) schedule(static)
    DO ir=1,fpar%nnr1
       v(ir)=CMPLX(elff(ir),0._real_8,kind=real_8)
    ENDDO
    CALL  fwfftn(v,.FALSE.,parai%allgrp)
    !$omp parallel do private(IG) schedule(static)
    DO ig=1,ncpw%nhg
       vtemp(ig) = v(nzh(ig))
    ENDDO
    IF (cntl%tpath) THEN
       IF (ipcurr.EQ.0) THEN
          filen='ELF'
       ELSE
          cflbod='ELF_'
          IF (paral%io_parent)&
               WRITE(cipnum,'(I4)') ipcurr
          CALL xstring(cflbod,n1,n2)
          CALL xstring(cipnum,i1,i2)
          filen=cflbod(n1:n2)//cipnum(i1:i2)
       ENDIF
    ELSEIF (tmw) THEN
       cflbod='ELF_'
       CALL mw_filename(cflbod,filen,mwi%walker_id)
    ELSE
       filen='ELF'
    ENDIF
    ! ==--------------------------------------------------------------==
    ! CB: Filenames for the BS states
    IF (cntl%bsymm) THEN
       space=INDEX(filen,' ')
       IF (bsclcs.EQ.1) THEN
          IF (space.GT.10) THEN
             filen='DENSITY.bs'
          ELSE
             filen(space:space+2)='.bs'
          ENDIF
       ELSE
          IF (space.GT.10) THEN
             filen='DENSITY.hs'
          ELSE
             filen(space:space+2)='.hs'
          ENDIF
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL densto(vtemp,tau0,filen)
    elfmin=1.e30_real_8
    elfmax=-1.e30_real_8
    DO ir=1,fpar%nnr1
       elfmax=MAX(elfmax,elff(ir))
       elfmin=MIN(elfmin,elff(ir))
    ENDDO
    CALL mp_max(elfmax,parai%allgrp)
    elfmin = -elfmin
    CALL mp_max(elfmin,parai%allgrp)
    elfmin = -elfmin
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,'(A,T50,F15.8)')  ' ELF: MINIMUM OF ELF   = ',elfmin
       IF (paral%io_parent)&
            WRITE(6,'(A,T52,F13.8)')  ' ELF: MAXIMUM OF ELF   = ',elfmax
    ENDIF
    DEALLOCATE(elff,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    DEALLOCATE(vtemp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE elf
  ! ==================================================================
  SUBROUTINE give_scr_elf(lelf,tag,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lelf
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate

    INTEGER                                  :: lrhoofr, lrnlsm

! Variables
! ==--------------------------------------------------------------==
! VTEMP(2*NHG)

    lelf=2*ncpw%nhg
    CALL give_scr_rnlsm(lrnlsm,tag,nstate,.FALSE.)
    CALL give_scr_rhoofr(lrhoofr,tag)
    lelf=MAX(lelf,lrnlsm,lrhoofr)
    tag='MAX(LELF,LRNLSM,LRHOOFR)'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_elf
  ! ==================================================================

END MODULE elf_utils
