MODULE rhopri_utils
  USE bsym,                            ONLY: bsclcs
  USE cdftmod,                         ONLY: cdfthda,&
                                             cdftlog
  USE cppt,                            ONLY: nzh
  USE eicalc_utils,                    ONLY: eicalc
  USE elf_utils,                       ONLY: elf,&
                                             give_scr_elf
  USE elstpo_utils,                    ONLY: elstpo
  USE ener,                            ONLY: chrg,&
                                             ener_com
  USE error_handling,                  ONLY: stopgm
  USE fftmain_utils,                   ONLY: fwfftn
  USE hip_utils,                       ONLY: give_qphi
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE lsd_elf_utils,                   ONLY: give_scr_lsd_elf,&
                                             lsd_elf
  USE meta_multiple_walkers_utils,     ONLY: mw_filename
  USE mw,                              ONLY: mwi,&
                                             tmw
  USE noforce_utils,                   ONLY: give_scr_noforce
  USE norhoe_utils,                    ONLY: norhoe
  USE parac,                           ONLY: parai,&
                                             paral
  USE phfac_utils,                     ONLY: phfac
  USE pimd,                            ONLY: ipcurr
  USE prden,                           ONLY: elfcb,&
                                             mwfn,&
                                             numpr
  USE readsr_utils,                    ONLY: xstring
  USE response_pmod,                   ONLY: response1
  USE rhoofr_c_utils,                  ONLY: rhoofr_c
  USE rhoofr_utils,                    ONLY: give_scr_rhoofr,&
                                             rhoofr
  USE rnlsm_utils,                     ONLY: give_scr_rnlsm,&
                                             rnlsm
  USE ropt,                            ONLY: iteropt
  USE spin,                            ONLY: clsd,&
                                             spin_mod
  USE store_types,                     ONLY: rout1
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             ncpw,&
                                             nkpt
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: unitmx
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rhopri
  PUBLIC :: give_scr_rhopri
  PUBLIC :: rhopri_eh
  !public :: currpri

CONTAINS

  ! ==================================================================
  SUBROUTINE rhopri(c0,tau0,rhoe,psi,nstate,nkpoint)
    ! ==--------------------------------------------------------------==
    ! ==  CALCULATE THE ELECTRONIC DENSITY IN G-SPACE AND DUMP        ==
    ! ==  ALL DATA TO FILE DENSITY                                    ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8), TARGET                  :: c0(:,:,:)
    REAL(real_8)                             :: tau0(:,:,:), rhoe(:,:)
    COMPLEX(real_8)                          :: psi(:)
    INTEGER                                  :: nstate, nkpoint

    CHARACTER(*), PARAMETER                  :: procedureN = 'rhopri'

    CHARACTER(len=14)                        :: cflbod, cipnum
    CHARACTER(len=30)                        :: filen
    COMPLEX(real_8), ALLOCATABLE             :: eirop(:), eivps(:), qphi(:), &
                                                vtemp(:)
    COMPLEX(real_8), ALLOCATABLE, SAVE       :: sc0_scr(:,:)
    INTEGER                                  :: i, i1, i2, ierr, ig, ii, ik, &
                                                ikind, il_qphi, is, ispst, &
                                                isub, n1, n2, nlsdbk, space
    INTEGER, SAVE                            :: icall = 0, icon = 0
    LOGICAL                                  :: ldowrt, tlsdbk
    REAL(real_8)                             :: ekin_bak, eta, svsumg, svsumr
    REAL(real_8), ALLOCATABLE, SAVE          :: eigv_scr(:), xmat_1(:,:), &
                                                xmat_2(:,:)

!c0(nkpt%ngwk,nstate,nkpoint)

    CALL tiset(procedureN,isub)
    ! ==--------------------------------------------------------------==
    CALL phfac(tau0)
    DO ikind=1,nkpoint
       CALL rnlsm(c0(:,:,ikind),nstate,1,ikind,.FALSE.)
    ENDDO
    ! calling RHOOFR will overwrite EKIN so we need to make
    ! a backup. the rest of COMMON /ENER/ is untouched. AK+NN.
    ekin_bak=ener_com%ekin
    IF (tkpts%tkpnt) THEN
       CALL rhoofr_c(c0,rhoe,psi,nstate)
    ELSE
       IF (cntl%nonort)THEN
          IF (icon.EQ.0)THEN
             ALLOCATE(sc0_scr(ncpw%ngw,nstate),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             ALLOCATE(eigv_scr(nstate),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             ALLOCATE(xmat_1(nstate,nstate),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             ALLOCATE(xmat_2(nstate,nstate),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             CALL unitmx(xmat_1,nstate)
          ENDIF
          icon=1
          CALL norhoe(c0,sc0_scr,eigv_scr,rhoe,psi,xmat_1,xmat_2,&
               nstate)
       ELSE
          CALL rhoofr(c0(:,:,1),rhoe,psi,nstate)
       ENDIF
    ENDIF
    ! 
    ener_com%ekin=ekin_bak
    ALLOCATE(vtemp(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    !$omp parallel do private(I)
    DO i=1,fpar%nnr1
       psi(i)=CMPLX(rhoe(i,1),0._real_8,kind=real_8)
    ENDDO
    !$omp end parallel do
    CALL  fwfftn(psi,.FALSE.,parai%allgrp)
    !$omp parallel do private(IG)
    DO ig=1,ncpw%nhg
       vtemp(ig) = psi(nzh(ig))
    ENDDO
    !$omp end parallel do
    IF (cntl%tpath) THEN
       IF (ipcurr.EQ.0) THEN
          filen='DENSITY'
       ELSE
          cflbod='DENSITY_'
          IF (paral%io_parent)&
               WRITE(cipnum,'(I4)') ipcurr
          CALL xstring(cflbod,n1,n2)
          CALL xstring(cipnum,i1,i2)
          filen=cflbod(n1:n2)//cipnum(i1:i2)
       ENDIF
    ELSEIF (response1%tdummyatom_ref) THEN
       filen='rhoofr_ref.dat'
       ! EHR[
    ELSEIF (cntl%tmdeh) THEN
       cflbod='DENSITY_'
       IF (paral%io_parent)&
            WRITE(cipnum,'(I4)') icall
       CALL xstring(cflbod,n1,n2)
       CALL xstring(cipnum,i1,i2)
       filen=cflbod(n1:n2)//cipnum(i1:i2)
       icall=icall+1
       ! EHR]
    ELSEIF (tmw) THEN
       cflbod='DENSITY_'
       CALL mw_filename(cflbod,filen,mwi%walker_id)
    ELSE
       IF (cdftlog%thda)THEN
          IF (cdfthda%hdafirst)THEN
             filen='DENSITY_1'
          ELSE
             filen='DENSITY_2'
          ENDIF
       ELSE
          filen='DENSITY'
       ENDIF
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
    CALL xstring(filen,i1,i2)
    IF (rout1%nrhoout.GT.0) THEN
       IF (paral%io_parent)&
            WRITE(cipnum,'(I14)') iteropt%nfi
       CALL xstring(cipnum,n1,n2)
       filen=filen(i1:i2) // '.' // cipnum(n1:n2)
       CALL xstring(filen,i1,i2)
    ENDIF
    CALL densto(vtemp,tau0,filen(i1:i2))
    ! 
    IF (cntl%tlsd) THEN
       !$omp parallel do private(I,ETA)
       DO i=1,fpar%nnr1
          eta=rhoe(i,1)-2._real_8*rhoe(i,2)
          psi(i)=CMPLX(eta,0._real_8,kind=real_8)
       ENDDO
       !$omp end parallel do
       CALL  fwfftn(psi,.FALSE.,parai%allgrp)
       !$omp parallel do private(IG)
       DO ig=1,ncpw%nhg
          vtemp(ig) = psi(nzh(ig))
       ENDDO
       !$omp end parallel do
       IF (cntl%tpath) THEN
          IF (ipcurr.EQ.0) THEN
             filen='SPINDEN'
          ELSE
             cflbod='SPINDEN_'
             IF (paral%io_parent)&
                  WRITE(cipnum,'(I4)') ipcurr
             CALL xstring(cflbod,n1,n2)
             CALL xstring(cipnum,i1,i2)
             filen=cflbod(n1:n2)//cipnum(i1:i2)
          ENDIF
          ! EHR[
       ELSEIF (cntl%tmdeh) THEN
          cflbod='SPINDEN_'
          IF (paral%io_parent)&
               WRITE(cipnum,'(I4)') icall
          CALL xstring(cflbod,n1,n2)
          CALL xstring(cipnum,i1,i2)
          filen=cflbod(n1:n2)//cipnum(i1:i2)
          icall=icall+1
          ! EHR]
       ELSEIF (tmw) THEN
          cflbod='SPINDEN_'
          CALL mw_filename(cflbod,filen,mwi%walker_id)
       ELSE
          IF (cdftlog%thda)THEN
             IF (cdfthda%hdafirst)THEN
                filen='SPINDEN_1'
             ELSE
                filen='SPINDEN_2'
             ENDIF
          ELSE
             filen='SPINDEN'
          ENDIF
       ENDIF
       ! ==--------------------------------------------------------------==
       ! CB: Filenames for the BS states
       IF (cntl%bsymm) THEN
          space=INDEX(filen,' ')
          IF (bsclcs.EQ.1) THEN
             IF (space.GT.10) THEN
                filen='SPINDEN.bs'
             ELSE
                filen(space:space+2)='.bs'
             ENDIF
          ELSE
             IF (space.GT.10) THEN
                filen='SPINDEN.hs'
             ELSE
                filen(space:space+2)='.hs'
             ENDIF
          ENDIF
       ENDIF
       ! ==--------------------------------------------------------------==
       CALL xstring(filen,i1,i2)
       IF (rout1%nrhoout.GT.0) THEN
          CALL xstring(cipnum,n1,n2)
          IF (paral%io_parent)&
               WRITE(cipnum,'(I14)') iteropt%nfi
          filen=filen(i1:i2) // '.' // cipnum(n1:n2)
          CALL xstring(filen,i1,i2)
       ENDIF
       CALL densto(vtemp,tau0,filen(i1:i2))
    ENDIF
    ! 
    ! Save CSUMG/CSUMR from all electron calc.
    svsumg=chrg%csumg
    svsumr=chrg%csumr
    IF (numpr.GE.1) THEN
       DO is=1,clsd%nlsd
          DO ii=1,numpr
             ispst=ABS(mwfn(ii))
             IF (is.EQ.2) ispst=ispst+spin_mod%nsup
             CALL zeroing(vtemp)!,nhg)
             ! ..Store density associated with wavefunction ISPST
             DO ik=1,nkpoint
                IF (mwfn(ii).GT.0.AND.ispst.LE.nstate) THEN
                   nlsdbk=clsd%nlsd
                   tlsdbk=cntl%tlsd
                   clsd%nlsd=1
                   cntl%tlsd=.FALSE.
                   IF (tkpts%tkpnt) THEN
                      CALL rhoofr_c(c0(:,ispst,ik),rhoe,psi,1)
                   ELSE
                      CALL rhoofr(c0(:,ispst:ispst,1),rhoe,psi,1)
                   ENDIF
                   clsd%nlsd=nlsdbk
                   cntl%tlsd=tlsdbk
                   !$omp parallel do private(I)
                   DO i=1,fpar%nnr1
                      psi(i)=CMPLX(rhoe(i,1),0._real_8,kind=real_8)
                   ENDDO
                   !$omp end parallel do
                   CALL  fwfftn(psi,.FALSE.,parai%allgrp)
                   !$omp parallel do private(IG)
                   DO ig=1,ncpw%nhg
                      vtemp(ig) = psi(nzh(ig))
                   ENDDO
                   !$omp end parallel do
                ELSEIF (ispst.LE.nstate) THEN
                   ! ..Store directly the wavefunction
                   IF (tkpts%tkpnt) THEN
                      IF (paral%io_parent)&
                           WRITE(6,*)&
                           '!!!! WRITING OUT COMPLEX WAVFUNCTIONS (I.E. KPOINTS) NOT'&
                           ,'IMPLEMENTED'
                      GOTO 100
                   ELSE
                      !$omp parallel do private(IG)
                      DO ig=1,nkpt%ngwk
                         vtemp(ig) = c0(ig,ABS(ispst),1)
                      ENDDO
                      !$omp end parallel do
                   ENDIF
                ENDIF
                IF ( mwfn(ii) .LT. 0 ) THEN
                   IF (cntl%tlsd) THEN
                      IF (is.EQ.1) THEN
                         cflbod='WAVEFUNCTION.A'
                         ldowrt=(ABS(mwfn(ii)).LE.spin_mod%nsup)
                      ELSE
                         cflbod='WAVEFUNCTION.B'
                         ldowrt=(ABS(mwfn(ii)).LE.spin_mod%nsdown)
                      ENDIF
                   ELSE
                      cflbod='WAVEFUNCTION.'
                      ldowrt=(ispst.LE.nstate)
                   ENDIF
                ELSE
                   IF (cntl%tlsd) THEN
                      IF (is.EQ.1) THEN
                         cflbod='DENSITY.A'
                         ldowrt=(ABS(mwfn(ii)).LE.spin_mod%nsup)
                      ELSE
                         cflbod='DENSITY.B'
                         ldowrt=(ABS(mwfn(ii)).LE.spin_mod%nsdown)
                      ENDIF
                   ELSE
                      cflbod='DENSITY.'
                      ldowrt=(ispst.LE.nstate)
                   ENDIF
                ENDIF
                IF (ldowrt) THEN
                   IF (paral%io_parent)&
                        WRITE(cipnum,'(I4)') ABS(mwfn(ii))
                   IF (tmw) THEN
                      CALL mw_filename(cflbod,filen,mwi%walker_id)
                      cflbod=filen
                   ENDIF
                   CALL xstring(cflbod,n1,n2)
                   CALL xstring(cipnum,i1,i2)
                   filen=cflbod(n1:n2)//cipnum(i1:i2)
                   IF (tkpts%tkpnt) THEN
                      n1=LEN(filen)
                      IF (n1.GT.14) n1=14
                      cflbod=filen(1:n1)
                      IF (paral%io_parent)&
                           WRITE(cipnum,'(I4)') ik
                      CALL xstring(cflbod,n1,n2)
                      CALL xstring(cipnum,i1,i2)
                      filen=cflbod(n1:n2)//'.KP.'//cipnum(i1:i2)
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
                   CALL xstring(filen,i1,i2)
                   IF (rout1%nrhoout.GT.0) THEN
                      IF (paral%io_parent)&
                           WRITE(cipnum,'(I14)') iteropt%nfi
                      CALL xstring(cipnum,n1,n2)
                      filen=filen(i1:i2) // '.' // cipnum(n1:n2)
                      CALL xstring(filen,i1,i2)
                   ENDIF
                   CALL densto(vtemp,tau0,filen(i1:i2))
                ENDIF! LDOWRT
             ENDDO! IK
          ENDDO ! II
       ENDDO     ! IS
    ENDIF
    ! Restore CSUMG/CSUMR from all electron calc.
    chrg%csumg=svsumg
    chrg%csumr=svsumr
100 CONTINUE
    ! 
    ! ..Electrostatic potential
    IF (cntl%tepot) THEN
       ALLOCATE(eivps(ncpw%nhg),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(eirop(ncpw%nhg),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)

       CALL give_qphi(il_qphi)
       ALLOCATE(qphi(il_qphi),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)

       CALL eicalc(eivps,eirop)
       !$omp parallel do private(I)
       DO i=1,fpar%nnr1
          psi(i)=CMPLX(rhoe(i,1),0._real_8,kind=real_8)
       ENDDO
       !$omp end parallel do
       CALL  fwfftn(psi,.FALSE.,parai%allgrp)
       !$omp parallel do private(IG)
       DO ig=1,ncpw%nhg
          vtemp(ig) = psi(nzh(ig))
       ENDDO
       !$omp end parallel do
       CALL elstpo(vtemp,eirop,eivps,psi,qphi)
       IF (cntl%tpath) THEN
          IF (ipcurr.EQ.0) THEN
             filen='ELPOT'
          ELSE
             cflbod='ELPOT_'
             IF (paral%io_parent)&
                  WRITE(cipnum,'(I4)') ipcurr
             CALL xstring(cflbod,n1,n2)
             CALL xstring(cipnum,i1,i2)
             filen=cflbod(n1:n2)//cipnum(i1:i2)
          ENDIF
       ELSEIF (tmw) THEN
          cflbod='ELPOT_'
          CALL mw_filename(cflbod,filen,mwi%walker_id)
       ELSE
          filen='ELPOT'
       ENDIF
       ! ==--------------------------------------------------------------==
       ! CB: Filenames for the BS states
       IF (cntl%bsymm) THEN
          space=INDEX(filen,' ')
          IF (bsclcs.EQ.1) THEN
             IF (space.GT.10) THEN
                filen='ELPOT.bs'
             ELSE
                filen(space:space+2)='.bs'
             ENDIF
          ELSE
             IF (space.GT.10) THEN
                filen='ELPOT.hs'
             ELSE
                filen(space:space+2)='.hs'
             ENDIF
          ENDIF
       ENDIF
       ! ==--------------------------------------------------------------==
       CALL xstring(filen,i1,i2)
       IF (rout1%nrhoout.GT.0) THEN
          IF (paral%io_parent)&
               WRITE(cipnum,'(I14)') iteropt%nfi
          CALL xstring(cipnum,n1,n2)
          filen=filen(i1:i2) // '.' // cipnum(n1:n2)
          CALL xstring(filen,i1,i2)
       ENDIF
       CALL densto(psi,tau0,filen(i1:i2))

       DEALLOCATE(eirop,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(eivps,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(qphi,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    ENDIF
    ! ..ELF
    IF (elfcb%telf) CALL elf(c0,tau0,rhoe,psi,nstate,nkpoint)
    ! ..ELF for the unrestricted case
    IF (elfcb%telf.AND.cntl%tlsd)&
         CALL lsd_elf(C0,TAU0,RHOE,PSI,NSTATE,NKPOINT)

    DEALLOCATE(vtemp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rhopri
  ! ==================================================================
  SUBROUTINE give_scr_rhopri(lrhopri,tag,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lrhopri
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate

    INTEGER :: il_auxc, il_ddia, il_gam, il_qphi = 0, il_smat, lelf = 0, &
      lnorho = 0, lrhoofr = 0, lrnlsm = 0

    CALL give_scr_rnlsm(lrnlsm,tag,nstate,.FALSE.)
    CALL give_scr_rhoofr(lrhoofr,tag)
    IF (cntl%tepot) THEN
       ! VTEMP(2*NHG),EIVPS(2*NHG),EIROP(2*NHG)
       lrhopri=6*ncpw%nhg
    ELSE
       ! VTEMP(2*NHG)
       lrhopri=2*ncpw%nhg
    ENDIF
    ! QPHI
    CALL give_qphi(il_qphi)
    lrhopri=lrhopri+il_qphi
    IF (cntl%nonort) CALL give_scr_noforce(lnorho,il_gam,il_auxc,&
         il_smat,il_ddia,tag,nstate,.FALSE.)
    IF (elfcb%telf) CALL give_scr_elf(lelf,tag,nstate)
    IF (elfcb%telf.AND.cntl%tlsd) CALL give_scr_lsd_elf(lelf,tag,nstate)
    lrhopri=MAX(lrnlsm,lrhoofr,lrhopri,lelf,lnorho)
    tag='MAX(LRNLSM,LRHOOFR,LRHOPRI...)'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_rhopri
  ! ==================================================================
  SUBROUTINE rhopri_eh(c0,tau0,rhoe,psi,nstate,nkpoint)
    ! ==--------------------------------------------------------------==
    ! ==  DUMP THE DENSITY IN DDENSITY_ICALL                          ==
    ! ==  IN: RHOE in real space                                      ==
    ! ==  ICALL is a local incremental number (SAVE STATEMENT)        ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:), rhoe(:,:)
    COMPLEX(real_8)                          :: psi(:)
    INTEGER                                  :: nstate, nkpoint
    COMPLEX(real_8)                          :: c0(nkpt%ngwk,nstate,nkpoint)

    CHARACTER(*), PARAMETER                  :: procedureN = 'rhopri_eh'

    CHARACTER(len=12)                        :: cflbod, cipnum
    CHARACTER(len=15)                        :: filen
    COMPLEX(real_8), ALLOCATABLE             :: vtemp(:)
    INTEGER                                  :: i, i1, i2, ierr, ig, isub, &
                                                n1, n2
    INTEGER, SAVE                            :: icall = 0

! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    ! ==--------------------------------------------------------------==
    ALLOCATE(vtemp(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    DO i=1,fpar%nnr1
       psi(i)=CMPLX(rhoe(i,1),0._real_8,kind=real_8)
    ENDDO
    CALL  fwfftn(psi,.FALSE.,parai%allgrp)
    DO ig=1,ncpw%nhg
       vtemp(ig) = psi(nzh(ig))
    ENDDO
    ! 
    cflbod='DDENSITY_'
    IF (paral%io_parent)&
         WRITE(cipnum,'(I4)') icall
    CALL xstring(cflbod,n1,n2)
    CALL xstring(cipnum,i1,i2)
    filen=cflbod(n1:n2)//cipnum(i1:i2)
    CALL densto(vtemp,tau0,filen)
    icall=icall+1
    CALL tihalt(procedureN,isub)

    DEALLOCATE(vtemp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rhopri_eh
  ! ==================================================================
  SUBROUTINE currpri(tau0,rhoe,psi,nstate,nkpoint,icoor)
    ! ==--------------------------------------------------------------==
    ! ==  DUMP THE CURRENT IN RHOE(NNR1,I) TO THE FILE CURRENT_       ==
    ! ==  THE CURRENT IS IN RHOE (in real space)                      ==
    ! ==  THE VARIABLE ICORR controls the output (component=x,y,...)  ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:), rhoe(:,:)
    COMPLEX(real_8)                          :: psi(:)
    INTEGER                                  :: nstate, nkpoint, icoor

    CHARACTER(*), PARAMETER                  :: procedureN = 'currpri'

    CHARACTER(len=12)                        :: cflbod, cipnum
    CHARACTER(len=15)                        :: filen
    COMPLEX(real_8), ALLOCATABLE             :: vtemp(:)
    INTEGER                                  :: i, i1, i2, ierr, ig, isub, &
                                                n1, n2
    INTEGER, SAVE                            :: icall = 0

! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    ! ==--------------------------------------------------------------==
    ALLOCATE(vtemp(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    DO i=1,fpar%nnr1
       psi(i)=CMPLX(rhoe(i,icoor),0._real_8,kind=real_8)
    ENDDO
    CALL  fwfftn(psi,.FALSE.,parai%allgrp)
    DO ig=1,ncpw%nhg
       vtemp(ig) = psi(nzh(ig))
    ENDDO
    ! 
    IF (icoor.EQ.0) THEN
       cflbod='CURRENT_'
    ELSEIF (icoor.EQ.1) THEN
       cflbod='CURRENTX_'
    ELSEIF (icoor.EQ.2) THEN
       cflbod='CURRENTY_'
    ELSEIF (icoor.EQ.3) THEN
       cflbod='CURRENTZ_'
    ELSEIF (icoor.EQ.4) THEN
       cflbod='CURRENTD_'
    ELSEIF (icoor.EQ.5) THEN
       cflbod='CURRENTT_'
    ENDIF
    IF (paral%io_parent)&
         WRITE(cipnum,'(I4)') icall
    CALL xstring(cflbod,n1,n2)
    CALL xstring(cipnum,i1,i2)
    filen=cflbod(n1:n2)//cipnum(i1:i2)
    CALL densto(vtemp,tau0,filen)
    IF (icoor.EQ.0) icall=icall+1
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE currpri
  ! ==================================================================
  SUBROUTINE field_pr(tau0,rhoe,psi,nstate)
    ! ==--------------------------------------------------------------==
    ! ==  DUMP the external firld FIELD in RHOE (real space) to the   ==
    ! ==  FILE FIELD_ICALL                                            ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:), rhoe(:)
    COMPLEX(real_8)                          :: psi(:)
    INTEGER                                  :: nstate

    CHARACTER(*), PARAMETER                  :: procedureN = 'field_pr'

    CHARACTER(len=12)                        :: cflbod, cipnum
    CHARACTER(len=15)                        :: filen
    COMPLEX(real_8), ALLOCATABLE             :: vtemp(:)
    INTEGER                                  :: i, i1, i2, ierr, ig, isub, &
                                                n1, n2
    INTEGER, SAVE                            :: icall = 0

! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    ! ==--------------------------------------------------------------==
    ALLOCATE(vtemp(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    DO i=1,fpar%nnr1
       psi(i)=CMPLX(rhoe(i),0._real_8,kind=real_8)
    ENDDO
    CALL  fwfftn(psi,.FALSE.,parai%allgrp)
    DO ig=1,ncpw%nhg
       vtemp(ig) = psi(nzh(ig))
    ENDDO
    ! 
    cflbod='FIELD_'
    IF (paral%io_parent)&
         WRITE(cipnum,'(I4)') icall
    CALL xstring(cflbod,n1,n2)
    CALL xstring(cipnum,i1,i2)
    filen=cflbod(n1:n2)//cipnum(i1:i2)
    CALL densto(vtemp,tau0,filen)
    icall=icall+1
    CALL tihalt(procedureN,isub)

    DEALLOCATE(vtemp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE field_pr
  ! ==================================================================

END MODULE rhopri_utils
