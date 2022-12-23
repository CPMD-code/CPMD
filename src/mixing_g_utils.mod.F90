MODULE mixing_g_utils
  USE anderson_utils,                  ONLY: anderson_c,&
                                             change_xmix
  USE andr,                            ONLY: andr2,&
                                             andr3
  USE broy,                            ONLY: broy1
  USE broyden_utils,                   ONLY: broyden,&
                                             give_scr_broyden,&
                                             kerker
  USE cppt,                            ONLY: hg,&
                                             indz,&
                                             nzh
  USE error_handling,                  ONLY: stopgm
  USE fftmain_utils,                   ONLY: fwfftn,&
                                             invfftn
  USE fftnew_utils,                    ONLY: setfftn
  USE geq0mod,                         ONLY: geq0
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: parai
  USE rhodiis_utils,                   ONLY: change_nrdiis
  USE spin,                            ONLY: clsd
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             ncpw
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mixing_g
  PUBLIC :: give_scr_mixing

CONTAINS

  ! ==================================================================
  SUBROUTINE mixing_g(nfr,gemax,rin0,rout0,rmix,v,thl)
    ! ==--------------------------------------------------------------==
    ! == Routine to perform density mixing in g-space                 ==
    ! == RIN0 is input density in real space                          ==
    ! == ROUT0 is output density in real space                        ==
    ! == RMIX is mixed density  in real space                         ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nfr
    REAL(real_8) :: gemax, rin0(fpar%nnr1,clsd%nlsd), &
      rout0(fpar%nnr1,clsd%nlsd), rmix(fpar%nnr1,clsd%nlsd)
    COMPLEX(real_8)                          :: v(fpar%nnr1)
    REAL(real_8)                             :: thl(2)

    CHARACTER(*), PARAMETER                  :: procedureN = 'mixing_g'

    COMPLEX(real_8), ALLOCATABLE             :: rgin0(:,:), rgout0(:,:)
    COMPLEX(real_8), ALLOCATABLE, SAVE       :: rginm1(:,:), rgmix(:,:), &
                                                rgoutm1(:,:)
    INTEGER                                  :: ierr, ig, ilsd, ir, itb, nngx
    INTEGER, SAVE                            :: ifirst = 0

    CALL setfftn(0)
    IF (andr3%nrdiismax.GT.0) CALL stopgm('MIXING_G',&
         'DIIS MIXING NOT PROGRAMMED',& 
         __LINE__,__FILE__)
    nngx=ncpw%nhg*clsd%nlsd
    ! ==--------------------------------------------------------------==
    ! Test SCR length
    ALLOCATE(rgin0(ncpw%nhg, clsd%nlsd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(rgout0(ncpw%nhg, clsd%nlsd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    IF (ifirst.EQ.0) THEN
       ifirst=1
       IF (andr3%nrdiismax.GT.0) THEN
          ! CALL MEMORY(IP_RGITER, 2*NRDIISMAX*NNGX,'RGITER')
          ! CALL MEMORY(IP_DRGITER,2*NRDIISMAX*NNGX,'DRGITER')
          ! We can use Anderson mixing and after cntl%diis mixing.
          ! IP_RGINM1=IP_RGITER
          ! IP_RGOUTM1=IP_DRGITER
          ! CALL MEMORY(IP_DMAT,NRDIISMAX*NRDIISMAX*NLSD)
       ELSE
          ALLOCATE(rginm1(ncpw%nhg,clsd%nlsd),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(rgoutm1(ncpw%nhg,clsd%nlsd),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(rgmix(ncpw%nhg,clsd%nlsd),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          CALL zeroing(rginm1)!,nngx)
          CALL zeroing(rgoutm1)!,nngx)
       ENDIF
    ENDIF
    ! =================================================================
    ! Transform RIN0 and ROUT0 to g-space ( RGIN0 and RGOUT0)
    DO ilsd=1,clsd%nlsd
       !$omp parallel do private(IR)
       DO ir=1,fpar%nnr1
          v(ir)=CMPLX(rin0(ir,ilsd),0._real_8,kind=real_8)
       ENDDO
       CALL  fwfftn(v,.FALSE.,parai%allgrp)
       !$omp parallel do private(IG)
       DO ig=1,ncpw%nhg
          rgin0(ig,ilsd)=v(nzh(ig))
       ENDDO
       !$omp parallel do private(IR)
       DO ir=1,fpar%nnr1
          v(ir)=CMPLX(rout0(ir,ilsd),0._real_8,kind=real_8)
       ENDDO
       CALL  fwfftn(v,.FALSE.,parai%allgrp)
       !$omp parallel do private(IG)
       DO ig=1,ncpw%nhg
          rgout0(ig,ilsd)=v(nzh(ig))
       ENDDO
    ENDDO
    ! Initialization for density mixing (needed in cntl%md).
    ! IF(NFR.EQ.1) THEN
    ! IF(NRDIISMAX.GT.0) THEN
    ! INRDIIS=1
    ! CALL RHODIIS(0,RMIX(1,1),RIN0(1,1),ROUT0(1,1),RGITER(1,1),
    ! &         DRGITER(1,1),ANDRMIX,DMAT(1,1),0,1,THL(1))
    ! IF(cntl%tlsd) THEN
    ! CALL RHODIIS(0,RMIX(1,2),RIN0(1,2),ROUT0(1,2),RGITER(1,2),
    ! &           DRGITER(1,2),ANDRMIX,DMAT(1,1),0,2,THL(1))
    ! ENDIF
    ! ENDIF
    ! ENDIF
    ! Change parameter if needed.
    CALL change_xmix(andr2%andrmix,gemax)
    IF (andr3%nrdiismax.GT.0) CALL change_nrdiis(gemax)
    IF (cntl%tlsd) THEN
       IF (andr3%nrdiis.GT.0) THEN
          ! CALL RHODIIS(NNR1,RMIX(1,1),RIN0(1,1),ROUT0(1,1),
          ! *         RGITER(1,1),DRGITER(1,1),ANDRMIX,DMAT(1,1),
          ! &         NRDIIS,1,THL(1))
          ! CALL RHODIIS(NNR1,RMIX(1,2),RIN0(1,2),ROUT0(1,2),
          ! *         RGITER(1,2),DRGITER(1,2),ANDRMIX,DMAT(1,2),
          ! &         NRDIIS,2,THL(2))
       ELSE IF (broy1%tgbroy.AND.(nfr.GT.broy1%nfrbroy)) THEN
          itb=nfr-broy1%nfrbroy
          itb=MOD(itb-1,broy1%ibreset)+1
          CALL broyden(broy1%broymix,itb,rginm1(1,1),&
               rgin0(1,1),rgout0(1,1),rgmix(1,1),ncpw%nhg,broy1%ngbroy,1,hg)
          thl(1)=0._real_8
          CALL broyden(broy1%broymix,itb,rginm1(1,2),&
               rgin0(1,2),rgout0(1,2),rgmix(1,2),ncpw%nhg,broy1%ngbroy,2,hg)
          thl(2)=0._real_8
       ELSE IF (broy1%kermix.NE.0._real_8) THEN
          CALL kerker(broy1%broymix,broy1%kermix,rgin0(1,1),rgout0(1,1),&
               rgmix(1,1),ncpw%nhg,hg)
          thl(1)=0._real_8
          CALL kerker(broy1%broymix,broy1%kermix,rgin0(1,2),rgout0(1,2),&
               rgmix(1,2),ncpw%nhg,hg)
          thl(2)=0._real_8
       ELSE
          CALL anderson_c(andr2%andrmix,nfr,rginm1(:,1),rgoutm1(:,1),rgin0(:,1),&
               rgout0(:,1),rgmix(:,1),2*ncpw%nhg,thl(1))
          CALL anderson_c(andr2%andrmix,nfr,rginm1(:,2),rgoutm1(:,2),rgin0(:,2),&
               rgout0(:,2),rgmix(:,2),2*ncpw%nhg,thl(2))
       ENDIF
       CALL dcopy(2*nngx,rgin0,1,rginm1,1)
       CALL dcopy(2*nngx,rgout0,1,rgoutm1,1)
    ELSE
       IF (andr3%nrdiis.GT.0) THEN
          ! CALL RHODIIS(NNR1,RMIX,RIN0,ROUT0,RGITER,DRGITER,
          ! *                 ANDRMIX,DMAT,NRDIIS,1,THL)
          ! THL(1)=0._real_8
       ELSE
          IF (broy1%tgbroy.AND.(nfr.GT.broy1%nfrbroy)) THEN
             itb=nfr-broy1%nfrbroy
             ! IF(ITB.GT.IBRESET)ITB=ITB-IBRESET
             itb=MOD(itb-1,broy1%ibreset)+1
             CALL broyden(broy1%broymix,itb,rginm1,&
                  rgin0,rgout0,rgmix,ncpw%nhg,broy1%ngbroy,1,hg)
             thl(1)=0._real_8
          ELSE IF (broy1%kermix.NE.0._real_8) THEN
             CALL kerker(broy1%broymix,broy1%kermix,rgin0,rgout0,rgmix,ncpw%nhg,hg)
             thl(1)=0._real_8
          ELSE
             CALL anderson_c(andr2%andrmix,nfr,rginm1,rgoutm1,rgin0,rgout0,&
                  rgmix,2*ncpw%nhg,thl(1))
          ENDIF
          CALL dcopy(2*nngx,rgin0,1,rginm1,1)
          CALL dcopy(2*nngx,rgout0,1,rgoutm1,1)
       ENDIF
    ENDIF
    ! ..Transform RGMIX to to real-space
    DO ilsd=1,clsd%nlsd
       CALL zeroing(v)!,nnr1)
       !$omp parallel do private(IG)
       !CDIR NODEP
       DO ig=1,ncpw%nhg
          v(nzh(ig))=rgmix(ig,ilsd)
          v(indz(ig))=CONJG(rgmix(ig,ilsd))
       ENDDO
       IF (geq0) v(nzh(1))=rgmix(1,ilsd)
       CALL  invfftn(v,.FALSE.,parai%allgrp)
       !$omp parallel do private(IR)
       DO ir=1,fpar%nnr1
          rmix(ir,ilsd)=REAL(v(ir))
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    DEALLOCATE(rgin0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(rgout0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE mixing_g
  ! ==================================================================
  SUBROUTINE give_scr_mixing(lmixing,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lmixing
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: lbroyden

! ==--------------------------------------------------------------==

    IF (broy1%tgmix) THEN
       ! Mixing in G-space
       IF (broy1%tgbroy) THEN
          CALL give_scr_broyden(lbroyden,tag,ncpw%nhg)
       ELSE
          lbroyden=0
       ENDIF
       lmixing=lbroyden+2*2*ncpw%nhg*clsd%nlsd! RGIN0 RGOUT0
    ELSE
       ! Mixing in R-space
       lmixing=(andr3%nrdiismax+1)*(andr3%nrdiismax+1+1)! DM VB
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_mixing
  ! ==================================================================

END MODULE mixing_g_utils
