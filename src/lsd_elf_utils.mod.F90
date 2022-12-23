MODULE lsd_elf_utils
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
  USE spin,                            ONLY: spin_mod
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             ncpw,&
                                             nkpt,&
                                             parm
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: lsd_elf
  PUBLIC :: give_scr_lsd_elf

CONTAINS

  ! ==================================================================
  SUBROUTINE lsd_elf(c0,tau0,rhoe,v,nstate,nkpoint)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! %% THE ELECTRON LOCALIZATION FUNCTION FOR THE UNRESTRICTED CASE %%
    ! %% PLUS FOR ALPHA AND BETA SPIN INDIVIDUALLY                    %%
    ! %% use keywords: LSD and ELF an the 3 fcts are calculated       %%
    ! %% automatically                                                %%
    ! %% formulas based on: M. Kohut and A.Savin, INT. J. Quant.      %%
    ! %% Chem,60, 875-882 (1996), eq. (2) for alpha and beta spin,    %%
    ! %% eq (9) for the spin unrestricted case                        %%
    ! %% more information in the elf.F file                           %%
    ! %% 17/03/00 last check: Kirchner                                %%
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:), rhoe(fpar%nnr1,2)
    COMPLEX(real_8)                          :: v(maxfft)
    INTEGER                                  :: nstate, nkpoint
    COMPLEX(real_8)                          :: c0(nkpt%ngwk,nstate,nkpoint)

    CHARACTER(*), PARAMETER                  :: procedureN = 'lsd_elf'

    CHARACTER(len=12)                        :: cflbod, cflboda, cflbodb, &
                                                cipnum, cipnuma, cipnumb
    CHARACTER(len=15)                        :: filen, filena, filenb
    COMPLEX(real_8), ALLOCATABLE             :: va(:), vb(:), vtemp(:), &
                                                vtempa(:), vtempb(:)
    INTEGER                                  :: i, i1, i2, ierr, ig, ikind, &
                                                ir, k, n1, n2, space
    REAL(real_8) :: aelf, aelfa, aelfb, coefe, coefe23, coefe53, coefs, &
      elfmax, elfmaxa, elfmaxb, elfmin, elfmina, elfminb, p53, rhom
    REAL(real_8), ALLOCATABLE                :: elff(:), elffa(:), elffb(:)

    CALL setfftn(0)
    IF (paral%io_parent) WRITE(6,'(/,A)') ' CALCULATE LSD_ELF FUNCTION '
    ! ==--------------------------------------------------------------==
    ! Allocations of arrays
    ALLOCATE(elff(fpar%nnr1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(elffa(fpar%nnr1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(elffb(fpar%nnr1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(va(maxfft),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(vb(maxfft),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(vtemp(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(vtempa(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(vtempb(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL phfac(tau0)
    ! 
    DO ikind=1,nkpoint
       CALL rnlsm(c0(:,:,ikind),nstate,1,ikind,.FALSE.)
    ENDDO
    IF (tkpts%tkpnt) THEN
       CALL rhoofr_c(c0,rhoe,v,nstate)
    ELSE
       CALL rhoofr(c0(:,:,1),rhoe,v,nstate)
    ENDIF
    ! alpha spin: rho_alpha(IR)=(RHOE(IR,1)-RHOE(IR,2)
    DO ir=1,fpar%nnr1
       va(ir)=CMPLX((rhoe(ir,1)-rhoe(ir,2)),0._real_8,kind=real_8)
    ENDDO
    CALL  fwfftn(va,.FALSE.,parai%allgrp)
    DO ig=1,ncpw%nhg
       vtempa(ig) = va(nzh(ig))
    ENDDO
    ! beta spin: rho_beta(IR)=RHOE(IR,2)
    DO ir=1,fpar%nnr1
       vb(ir)=CMPLX(rhoe(ir,2),0._real_8,kind=real_8)
    ENDDO
    CALL  fwfftn(vb,.FALSE.,parai%allgrp)
    DO ig=1,ncpw%nhg
       vtempb(ig) = vb(nzh(ig))
    ENDDO
    ! .....calculate |nabla*rhoa|^2 and |nabla*rhob|^2
    CALL zeroing(va)!,maxfft)
    CALL zeroing(vb)!,maxfft)
    ! |nabla*rhoa| and |nabla*rhob|
    !CDIR NODEP
    DO ig=1,ncpw%nhg
       va(nzh(ig))=parm%tpiba*uimag*gk(1,ig)*vtempa(ig)
       vb(nzh(ig))=parm%tpiba*uimag*gk(1,ig)*vtempb(ig)
       va(indz(ig))=-parm%tpiba*uimag*gk(1,ig)*CONJG(vtempa(ig))
       vb(indz(ig))=-parm%tpiba*uimag*gk(1,ig)*CONJG(vtempb(ig))
    ENDDO
    CALL  invfftn(va,.FALSE.,parai%allgrp)
    CALL  invfftn(vb,.FALSE.,parai%allgrp)
    ! |nabla*rhoa|^2 and |nabla*rhob|^2
    DO ir=1,fpar%nnr1
       elffa(ir)=REAL(va(ir))*REAL(va(ir))
       elffb(ir)=REAL(vb(ir))*REAL(vb(ir))
    ENDDO

    CALL zeroing(va)!,maxfft)
    CALL zeroing(vb)!,maxfft)
    ! double FFT
    !CDIR NODEP
    DO ig=1,ncpw%nhg
       va(nzh(ig))=parm%tpiba*(uimag*gk(2,ig)-gk(3,ig))*vtempa(ig)
       vb(nzh(ig))=parm%tpiba*(uimag*gk(2,ig)-gk(3,ig))*vtempb(ig)
       va(indz(ig))=parm%tpiba*(-uimag*gk(2,ig)+gk(3,ig))*&
            CONJG(vtempa(ig))
       vb(indz(ig))=parm%tpiba*(-uimag*gk(2,ig)+gk(3,ig))*&
            CONJG(vtempb(ig))
    ENDDO
    CALL  invfftn(va,.FALSE.,parai%allgrp)
    CALL  invfftn(vb,.FALSE.,parai%allgrp)
    DO ir=1,fpar%nnr1
       elffa(ir)=elffa(ir)+REAL(va(ir)*CONJG(va(ir)))
       elffb(ir)=elffb(ir)+REAL(vb(ir)*CONJG(vb(ir)))
    ENDDO
    ! .....end calculate |nabla*rhoa|^2 and |nabla*rhob|^2
    ! 
    ! calculate C_F
    p53=5._real_8/3._real_8
    coefe=pi**2
    coefe=(3._real_8*coefe)**(2._real_8/3._real_8)
    coefe=0.3_real_8*coefe
    coefe23 = coefe*2**(2._real_8/3._real_8)
    coefe53 = coefe*2**p53
    ! 1/8 (|nabla*rhoa|^2/rhoa + |nabla*rhob|^2/rhob )/denominator
    DO ir=1,fpar%nnr1
       rhom=MAX(rhoe(ir,1),0._real_8)
       IF (rhom.GT.elfcb%elfcut) THEN
          elff(ir)=-elffa(ir)/(8._real_8*(rhoe(ir,1)-rhoe(ir,2)))
          elff(ir)=elff(ir)-elffb(ir)/(8._real_8*rhoe(ir,2))
          elff(ir) = elff(ir)/(coefe23*((rhoe(ir,1)&
               -RHOE(IR,2))**P53+RHOE(IR,2)**P53))
       ELSE
          elff(ir)=0._real_8
       ENDIF
       ! alpha part: 1/4*|nabla*rhoa|^2/rhoa/denominator_a
       rhom=MAX((rhoe(ir,1)-rhoe(ir,2)),0._real_8)
       IF (rhom.GT.(elfcb%elfcut)) THEN
          elffa(ir)=-elffa(ir)/(4._real_8*(rhoe(ir,1)-rhoe(ir,2)))
          elffa(ir)=elffa(ir)/(coefe53*&
               ((RHOE(IR,1)-RHOE(IR,2))**P53))
       ELSE
          elffa(ir)=0._real_8
       ENDIF
       ! beta part: 1/4*|nabla*rhoa|^2/rhoa/denominator_b
       rhom=MAX(rhoe(ir,2),0._real_8)
       IF (rhom.GT.elfcb%elfcut) THEN
          elffb(ir)=-elffb(ir)/(4._real_8*rhoe(ir,2))
          elffb(ir)=elffb(ir)/(coefe53*rhoe(ir,2)**p53)
       ELSE
          elffb(ir)=0._real_8
       ENDIF
    ENDDO
    ! 
    ! calculate 0.5*SUM_i{|nabla*phi_i|^2 / D_h }
    DO ikind=1,nkpoint
       DO i=1,nstate
          coefs=crge%f(i,ikind)/parm%omega
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
                IF (i .LE. spin_mod%nsup) THEN
                   rhom=MAX((rhoe(ir,1)-rhoe(ir,2)),0._real_8)
                   IF (rhom.GT.elfcb%elfcut) THEN
                      elffa(ir)=elffa(ir)+coefs*REAL(v(ir))**2&
                           /(COEFE53*(RHOE(IR,1)-RHOE(IR,2))**P53)
                   ENDIF
                ELSEIF (i .GT. spin_mod%nsup) THEN
                   rhom=MAX(rhoe(ir,2),0._real_8)
                   IF (rhom.GT.elfcb%elfcut) THEN
                      elffb(ir)=elffb(ir)+coefs*REAL(v(ir))**2&
                           /(COEFE53*RHOE(IR,2)**P53)
                   ENDIF
                ENDIF
                rhom=MAX(rhoe(ir,1),0._real_8)
                IF (rhom.GT.elfcb%elfcut) THEN
                   elff(ir)=elff(ir)+0.5_real_8*coefs*REAL(v(ir))**2/&
                        (COEFE23*((RHOE(IR,1)-RHOE(IR,2))**P53&
                        +RHOE(IR,2)**P53))
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    ! ..calculate 1/(1+(D/D_h)^2)
    DO ir=1,fpar%nnr1
       rhom=MAX(rhoe(ir,1),0._real_8)
       IF (rhom.GT.elfcb%elfcut) THEN
          aelf=elff(ir)+elfcb%elfeps/(coefe23*((rhoe(ir,1)&
               -RHOE(IR,2))**P53+RHOE(IR,2)**P53))
          elff(ir)=1._real_8/(1._real_8+aelf**2)
       ENDIF
       ! ..alpha elf
       rhom=MAX((rhoe(ir,1)-rhoe(ir,2)),0._real_8)
       IF (rhom.GT.elfcb%elfcut) THEN
          aelfa=elffa(ir)+elfcb%elfeps/(coefe53*&
               (RHOE(IR,1)-RHOE(IR,2))**P53)
          elffa(ir)=1._real_8/(1._real_8+aelfa**2)
       ENDIF
       ! ..beta elf
       rhom=MAX(rhoe(ir,2),0._real_8)
       IF (rhom.GT.elfcb%elfcut) THEN
          aelfb=elffb(ir)+elfcb%elfeps/(coefe53*rhoe(ir,2)**p53)
          elffb(ir)=1._real_8/(1._real_8+aelfb**2)
       ENDIF
    ENDDO

    filen='LSD_ELF'
    filena='ELF_ALPHA'
    filenb='ELF_BETA'
    ! ..   VTEMP=ELF in G-space
    DO ir=1,fpar%nnr1
       v(ir)=CMPLX(elff(ir),0._real_8,kind=real_8)
       va(ir)=CMPLX(elffa(ir),0._real_8,kind=real_8)
       vb(ir)=CMPLX(elffb(ir),0._real_8,kind=real_8)
    ENDDO

    CALL  fwfftn(v,.FALSE.,parai%allgrp)
    CALL  fwfftn(va,.FALSE.,parai%allgrp)
    CALL  fwfftn(vb,.FALSE.,parai%allgrp)

    DO ig=1,ncpw%nhg
       vtemp(ig) = v(nzh(ig))
       vtempa(ig) = va(nzh(ig))
       vtempb(ig) = vb(nzh(ig))
    ENDDO

    IF (cntl%tpath) THEN
       IF (ipcurr.EQ.0) THEN
          filen='LSD_ELF'
          filena='ELF_ALPHA'
          filenb='ELF_BETA'
       ELSE
          cflbod='LSD_ELF_'
          cflboda='ELF_ALPHA_'
          cflbodb='ELF_BETA_'
          IF (paral%io_parent)&
               WRITE(cipnum,'(I4)') ipcurr
          IF (paral%io_parent)&
               WRITE(cipnuma,'(I4)') ipcurr
          IF (paral%io_parent)&
               WRITE(cipnumb,'(I4)') ipcurr
          CALL xstring(cflbod,n1,n2)
          CALL xstring(cipnum,i1,i2)
          filen=cflbod(n1:n2)//cipnum(i1:i2)
          CALL xstring(cflboda,n1,n2)
          CALL xstring(cipnuma,i1,i2)
          filena=cflboda(n1:n2)//cipnuma(i1:i2)
          CALL xstring(cflbodb,n1,n2)
          CALL xstring(cipnumb,i1,i2)
          filenb=cflbodb(n1:n2)//cipnumb(i1:i2)
       ENDIF
    ELSEIF (tmw) THEN
       cflbod='LSD_ELF_'
       CALL mw_filename(cflbod,filen,mwi%walker_id)
       cflboda='ELF_ALPHA_'
       CALL mw_filename(cflboda,filena,mwi%walker_id)
       cflbodb='ELF_BETA_'
       CALL mw_filename(cflbodb,filenb,mwi%walker_id)
    ELSE
       filen='LSD_ELF'
       filena='ELF_ALPHA'
       filenb='ELF_BETA'
    ENDIF

    ! CB: Filenames for the BS states
    IF (cntl%bsymm) THEN
       space=INDEX(filen,' ')
       IF (bsclcs.EQ.1) THEN
          IF (space.GT.10) THEN
             filen='LSD_ELF.bs'
          ELSE
             filen(space:space+2)='.bs'
          ENDIF
       ELSE
          IF (space.GT.10) THEN
             filen='LSD_ELF.hs'
          ELSE
             filen(space:space+2)='.hs'
          ENDIF
       ENDIF
       space=INDEX(filen,'LSD_ELF')
       filena=filen
       filena(space:space+6)='ELF_ALP'
       filenb=filen
       filenb(space:space+6)='ELF_BET'
    ENDIF

    CALL densto(vtemp,tau0,filen)
    CALL densto(vtempa,tau0,filena)
    CALL densto(vtempb,tau0,filenb)
    elfmin=1.e30_real_8
    elfmax=-1.e30_real_8
    DO ir=1,fpar%nnr1
       elfmax=MAX(elfmax,elff(ir))
       elfmin=MIN(elfmin,elff(ir))
       elfmaxa=MAX(elfmax,elffa(ir))
       elfmina=MIN(elfmin,elffa(ir))
       elfmaxb=MAX(elfmax,elffb(ir))
       elfminb=MIN(elfmin,elffb(ir))
    ENDDO

    CALL mp_max(elfmax,parai%allgrp)
    elfmin = -elfmin
    CALL mp_max(elfmin,parai%allgrp)
    elfmin = -elfmin

    CALL mp_max(elfmaxa,parai%allgrp)
    elfmina = -elfmina
    CALL mp_max(elfmina,parai%allgrp)
    elfmina = -elfmina

    CALL mp_max(elfmaxb,parai%allgrp)
    elfminb = -elfminb
    CALL mp_max(elfminb,parai%allgrp)
    elfminb = -elfminb

    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,'(A,T50,F15.8)')&
            'LSD_ELF: MINIMUM OF LSD_ELF  = ',ELFMIN
       IF (paral%io_parent)&
            WRITE(6,'(A,T52,F13.8)')&
            'LSD_ELF: MAXIMUM OF LSD_ELF  = ',ELFMAX
       IF (paral%io_parent)&
            WRITE(6,'(A,T50,F15.8)')&
            'ELFA: MINIMUM OF ELF_ALPHA  = ',ELFMINA
       IF (paral%io_parent)&
            WRITE(6,'(A,T52,F13.8)')&
            'ELFA: MAXIMUM OF ELF_ALPHA  = ',ELFMAXA
       IF (paral%io_parent)&
            WRITE(6,'(A,T50,F15.8)')&
            'ELFB: MINIMUM OF ELF_BETA  = ',ELFMINB
       IF (paral%io_parent)&
            WRITE(6,'(A,T52,F13.8)')&
            'ELFB: MAXIMUM OF ELF_BETA  = ',ELFMAXB
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Deallocations.
    DEALLOCATE(elff,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(elffa,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(elffb,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(va,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(vb,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(vtemp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(vtempa,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(vtempb,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE lsd_elf
  ! ==================================================================
  SUBROUTINE give_scr_lsd_elf(lelf,tag,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lelf
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate

    INTEGER                                  :: lrhoofr, lrnlsm

    lelf=2*ncpw%nhg
    CALL give_scr_rnlsm(lrnlsm,tag,nstate,.FALSE.)
    CALL give_scr_rhoofr(lrhoofr,tag)
    lelf=MAX(lelf,lrnlsm,lrhoofr)
    tag='MAX(LELF,LRNLSM,LRHOOFR)'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_lsd_elf
  ! ==================================================================

END MODULE lsd_elf_utils
