MODULE secder_utils
  USE adat,                            ONLY: elem
  USE calc_alm_utils,                  ONLY: give_scr_calc_alm
  USE cell,                            ONLY: cell_com
  USE cnst,                            ONLY: fbohr
  USE copot_utils,                     ONLY: give_scr_copot
  USE ddipo_utils,                     ONLY: give_scr_ddipo
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_app,&
                                             fo_def,&
                                             fo_mark
  USE fint,                            ONLY: fint1
  USE forcedr_utils,                   ONLY: give_scr_forcedr
  USE initrun_utils,                   ONLY: give_scr_initrun
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE isos,                            ONLY: isos1
  USE kddipo_utils,                    ONLY: give_scr_kddipo
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE lr_tddft_utils,                  ONLY: give_scr_lr_tddft
  USE nlcc,                            ONLY: corel
  USE ortho_utils,                     ONLY: give_scr_ortho
  USE parac,                           ONLY: paral
  USE prop,                            ONLY: prop1
  USE pslo,                            ONLY: pslo_com
  USE rhoofr_utils,                    ONLY: give_scr_rhoofr
  USE rmas,                            ONLY: rmass
  USE rnlsm_utils,                     ONLY: give_scr_rnlsm
  USE store_types,                     ONLY: restart1
  USE symm,                            ONLY: symmi
  USE symtrz_utils,                    ONLY: give_scr_symmat
  USE system,                          ONLY: cnti,&
                                             cntl
  USE updrho_utils,                    ONLY: give_scr_updrho
  USE updwf_utils,                     ONLY: give_scr_updwf
  USE utils,                           ONLY: unitmx
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: purgeh
  PUBLIC :: molvib
  PUBLIC :: vibeig
  PUBLIC :: purged
  PUBLIC :: mldgau
  !public :: settransd
  !public :: setrotatd
  !public :: purgeh
  !public :: settransh
  !public :: setrotath
  PUBLIC :: give_scr_secder
  PUBLIC :: writeaclimax

CONTAINS

  ! ==================================================================
  SUBROUTINE molvib(tau0,sder)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:), sder(:,:)

    CHARACTER(len=100)                       :: filen
    INTEGER                                  :: i, ia, is, j, k
    LOGICAL                                  :: ferror

    filen='MOLVIB'
    IF (paral%io_parent)&
         CALL fileopen(21,filen,fo_def,ferror)
    IF (paral%io_parent)&
         WRITE(21,'(A)') ' &CART'
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          IF (paral%io_parent)&
               WRITE(21,'(I6,3F15.8,F15.6)') ions0%iatyp(is),(tau0(k,ia,is),k=1,3),&
               rmass%pma0(is)
       ENDDO
    ENDDO
    IF (paral%io_parent)&
         WRITE(21,'(A)') ' &END '
    IF (paral%io_parent)&
         WRITE(21,'(A)') ' &FCON '
    DO i=1,3*ions1%nat
       IF (paral%io_parent)&
            WRITE(21,*) (sder(i,j),j=1,3*ions1%nat)
    ENDDO
    IF (paral%io_parent)&
         WRITE(21,'(A)') ' &END '
    IF (paral%io_parent)&
         CALL fileclose(21)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE molvib
  ! ==================================================================
  SUBROUTINE vibeig(vibe,sder,ndim,rw)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: vibe(:), sder(:,:)
    INTEGER                                  :: ndim
    LOGICAL                                  :: rw

    CHARACTER(len=100)                       :: filen
    INTEGER                                  :: i, j, k, ncol, nn
    LOGICAL                                  :: ferror

!(ndim,ndim)

    ncol=8
    filen='VIBEIGVEC'
    IF (rw) THEN
       IF (paral%io_parent)&
            CALL fileopen(21,filen,fo_def,ferror)
    ELSE
       IF (paral%io_parent)&
            CALL fileopen(21,filen,fo_app+fo_mark,ferror)
    ENDIF
    DO i=1,ndim,ncol
       nn=MIN(ncol,ndim-i+1)
       IF (paral%io_parent)&
            WRITE(21,'(14I12)') (j,j=i,i+nn-1)
       IF (paral%io_parent)&
            WRITE(21,'(14(F12.3))') (vibe(j),j=i,i+nn-1)
       IF (paral%io_parent)&
            WRITE(21,*)
       DO k=1,ndim
          IF (paral%io_parent)&
               WRITE(21,'(14F12.6)') (sder(k,j),j=i,i+nn-1)
       ENDDO
    ENDDO
    IF (paral%io_parent)&
         CALL fileclose(21)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE vibeig
  ! ==================================================================
  SUBROUTINE purged(sder,scrd,tau0)
    ! ==--------------------------------------------------------------==
    ! == PROJECT OUT TRANSLATIONS (AND ROTATIONS) FROM THE DYNAMICAL  ==
    ! == MATRIX (MASS WEIGHTED HESSIAN)                               ==
    ! ==                                                              ==
    ! == USE SUBROUTINE PURGEH() FOR THE HESSIAN (NOT MASS WEIGHTED)  ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: sder(:,:), scrd(:,:), &
                                                tau0(:,:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'purged'

    INTEGER                                  :: ia, ierr, ii, is, nd
    REAL(real_8)                             :: cx, cy, cz, mi(3)
    REAL(real_8), ALLOCATABLE                :: cmc(:,:), proj(:,:), &
                                                rot(:,:), tras(:,:)

    ALLOCATE(tras(3*ions1%nat,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(rot(3*ions1%nat,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(proj(3*ions1%nat,3*ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(cmc(3,ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ...translation
    CALL settransd(tras,ions1%nat)
    ! ...rotation
    IF (isos1%tisos.OR.isos1%tclust) THEN
       cx=0._real_8
       cy=0._real_8
       cz=0._real_8
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             cx=cx+tau0(1,ia,is)*rmass%pma0(is)/rmass%pmat0
             cy=cy+tau0(2,ia,is)*rmass%pma0(is)/rmass%pmat0
             cz=cz+tau0(3,ia,is)*rmass%pma0(is)/rmass%pmat0
          ENDDO
       ENDDO
       ii=0
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             ii=ii+1
             cmc(1,ii)=tau0(1,ia,is)-cx
             cmc(2,ii)=tau0(2,ia,is)-cy
             cmc(3,ii)=tau0(3,ia,is)-cz
          ENDDO
       ENDDO
       CALL setrotatd(rot,cmc,ions1%nat,mi)
    ELSE
       mi(1)=0.0_real_8
       mi(2)=0.0_real_8
       mi(3)=0.0_real_8
       CALL zeroing(rot)!,9*ions1%nat)
    ENDIF
    nd=3*ions1%nat
    ! ...set up projectors
    CALL unitmx(proj,nd)
    CALL dsyr('U',nd,-1._real_8,tras(1,1),1,proj,nd)
    CALL dsyr('U',nd,-1._real_8,tras(1,2),1,proj,nd)
    CALL dsyr('U',nd,-1._real_8,tras(1,3),1,proj,nd)
    ! tu   project out rotations only for non-zero moments of inertia
    DO ii=1,3
       IF (mi(ii) .GT. 1.0e-2_real_8) THEN
          CALL dsyr('U',nd,-1._real_8,rot(1,ii),1,proj,nd)
       ENDIF
    ENDDO
    ! ...apply projectors
    CALL dsymm('L','U',nd,nd,1._real_8,proj,nd,sder,nd,0._real_8,scrd,nd)
    CALL dsymm('R','U',nd,nd,1._real_8,proj,nd,scrd,nd,0._real_8,sder,nd)
    DEALLOCATE(tras,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rot,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(proj,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(cmc,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE purged
  ! ==================================================================
  SUBROUTINE settransd(tras,n)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: n
    REAL(real_8)                             :: tras(3,n,3)

    INTEGER                                  :: i, ia, is
    REAL(real_8)                             :: x

! variables
! ==--------------------------------------------------------------==

    CALL zeroing(tras)!,9*n)
    ! tu   components of TRAS: sqrt(atomic mass)
    ! tu   normalization of TRAS: scale by 1/sqrt(total mass)
    i=0
    DO is=1,ions1%nsp
       x=SQRT(rmass%pma0(is)/rmass%pmat0)
       DO ia=1,ions0%na(is)
          i=i+1
          tras(1,i,1)=x
          tras(2,i,2)=x
          tras(3,i,3)=x
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE settransd
  ! ==================================================================
  SUBROUTINE setrotatd(rot,cmc,n,mi)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: n
    REAL(real_8)                             :: cmc(3,n), rot(3,n,3), mi(3)

    INTEGER                                  :: i, ia, info, is, j
    REAL(real_8)                             :: aux(9), mit(3,3), x

! Variables
! ==--------------------------------------------------------------==
! tu   prepare the moments of inertia tensor

    CALL zeroing(mit)!,9)
    i=0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          i=i+1
          mit(1,1)=mit(1,1)+rmass%pma0(is)*(cmc(2,i)**2+cmc(3,i)**2)
          mit(2,1)=mit(2,1)-rmass%pma0(is)*cmc(1,i)*cmc(2,i)
          mit(3,1)=mit(3,1)-rmass%pma0(is)*cmc(1,i)*cmc(3,i)
          mit(2,2)=mit(2,2)+rmass%pma0(is)*(cmc(1,i)**2+cmc(3,i)**2)
          mit(3,2)=mit(3,2)-rmass%pma0(is)*cmc(2,i)*cmc(3,i)
          mit(3,3)=mit(3,3)+rmass%pma0(is)*(cmc(1,i)**2+cmc(2,i)**2)
       ENDDO
    ENDDO
    mit(1,2)=mit(2,1)
    mit(1,3)=mit(3,1)
    mit(2,3)=mit(3,2)
    ! tu   diagonalize the moment of inertia tensor
    CALL dsyev('V','U',3,mit,3,mi,aux,9,info)
    IF (info.NE.0) CALL stopgm('SETROTAT','DSYEV INFO',& 
         __LINE__,__FILE__)
    ! tu   mass weight the coordinates CMC
    i=0
    DO is=1,ions1%nsp
       x=SQRT(rmass%pma0(is))
       DO ia=1,ions0%na(is)
          i=i+1
          cmc(1,i)=x*cmc(1,i)
          cmc(2,i)=x*cmc(2,i)
          cmc(3,i)=x*cmc(3,i)
       ENDDO
    ENDDO
    ! tu   components of ROT: cross product of coordinates and principal
    ! tu   axes of inertia (for non-zero moments of inertia only!)
    ! tu   normalization of ROT: scale by 1/sqrt(moment of inertia)          
    DO j=1,3
       IF (mi(j) .GT. 1.0e-2_real_8) THEN
          x=1.0_real_8/SQRT(mi(j))
          DO i=1,n
             rot(1,i,j)=x*(mit(2,j)*cmc(3,i)-mit(3,j)*cmc(2,i))
             rot(2,i,j)=x*(mit(3,j)*cmc(1,i)-mit(1,j)*cmc(3,i))
             rot(3,i,j)=x*(mit(1,j)*cmc(2,i)-mit(2,j)*cmc(1,i))
          ENDDO
       ENDIF
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE setrotatd
  ! ==================================================================
  SUBROUTINE purgeh(sder,scrd,tau0)
    ! ==--------------------------------------------------------------==
    ! == PROJECT OUT TRANSLATIONS (AND ROTATIONS) FROM THE HESSIAN    ==
    ! == MATRIX (NOT MASS WEIGHTED)                                   ==
    ! ==                                                              ==
    ! == USE SUBROUTINE PURGED() FOR THE DYNAMICAL MATRIX (MASS       ==
    ! == WEIGHTED HESSIAN)                                            ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: sder(:,:), scrd(:,:), &
                                                tau0(:,:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'purgeh'

    INTEGER                                  :: ia, ierr, ii, is, nd
    REAL(real_8)                             :: cx, cy, cz
    REAL(real_8), ALLOCATABLE                :: cmc(:,:), proj(:,:), &
                                                rot(:,:), tras(:,:)

! Variables
! (3*nat,3)
! (3*nat,3)
! (3*nat,3*nat)
! (3,nat)
! ==--------------------------------------------------------------==

    ALLOCATE(tras(3*ions1%nat,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(rot(3*ions1%nat,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(proj(3*ions1%nat,3*ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(cmc(3,ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ...translation
    CALL settransh(tras,ions1%nat)
    ! ...rotation
    IF (isos1%tisos.OR.isos1%tclust) THEN
       cx=0._real_8
       cy=0._real_8
       cz=0._real_8
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             cx=cx+tau0(1,ia,is)*rmass%pma0(is)/rmass%pmat0
             cy=cy+tau0(2,ia,is)*rmass%pma0(is)/rmass%pmat0
             cz=cz+tau0(3,ia,is)*rmass%pma0(is)/rmass%pmat0
          ENDDO
       ENDDO
       ii=0
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             ii=ii+1
             cmc(1,ii)=tau0(1,ia,is)-cx
             cmc(2,ii)=tau0(2,ia,is)-cy
             cmc(3,ii)=tau0(3,ia,is)-cz
          ENDDO
       ENDDO
       CALL setrotath(rot,tras,cmc,ions1%nat)
    ELSE
       CALL zeroing(rot)!,9*ions1%nat)
    ENDIF
    nd=3*ions1%nat
    ! ...set up projectors
    CALL unitmx(proj,nd)
    CALL dsyr('U',nd,-1._real_8,tras(1,1),1,proj,nd)
    CALL dsyr('U',nd,-1._real_8,tras(1,2),1,proj,nd)
    CALL dsyr('U',nd,-1._real_8,tras(1,3),1,proj,nd)
    CALL dsyr('U',nd,-1._real_8,rot(1,1),1,proj,nd)
    CALL dsyr('U',nd,-1._real_8,rot(1,2),1,proj,nd)
    CALL dsyr('U',nd,-1._real_8,rot(1,3),1,proj,nd)
    ! ...apply projectors
    CALL dsymm('L','U',nd,nd,1._real_8,proj,nd,sder,nd,0._real_8,scrd,nd)
    CALL dsymm('R','U',nd,nd,1._real_8,proj,nd,scrd,nd,0._real_8,sder,nd)
    DEALLOCATE(tras,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rot,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(proj,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(cmc,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE purgeh
  ! ==================================================================
  SUBROUTINE settransh(tras,nat)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nat
    REAL(real_8)                             :: tras(3,nat,3)

    INTEGER                                  :: i
    REAL(real_8)                             :: sc

! variables
! ==--------------------------------------------------------------==

    CALL zeroing(tras)!,9*nat)
    !$omp parallel do private(I)
    DO i=1,nat
       tras(1,i,1)=1._real_8
       tras(2,i,2)=1._real_8
       tras(3,i,3)=1._real_8
    ENDDO
    sc=1._real_8/SQRT(REAL(nat,kind=real_8))
    CALL dscal(9*nat,sc,tras,1)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE settransh
  ! ==================================================================
  SUBROUTINE setrotath(rot,tras,cmc,nat)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nat
    REAL(real_8)                             :: cmc(3,nat), tras(3,nat,3), &
                                                rot(3,nat,3)

    INTEGER                                  :: i, j
    REAL(real_8)                             :: sc
    REAL(real_8), EXTERNAL                   :: ddot

    !$omp parallel do private(I)
    DO i=1,nat
       rot(1,i,1)=     0._real_8
       rot(2,i,1)= cmc(3,i)
       rot(3,i,1)=-cmc(2,i)
       rot(1,i,2)=-cmc(3,i)
       rot(2,i,2)=     0._real_8
       rot(3,i,2)= cmc(1,i)
       rot(1,i,3)= cmc(2,i)
       rot(2,i,3)=-cmc(1,i)
       rot(3,i,3)=     0._real_8
    ENDDO
    DO i=1,3
       DO j=1,3
          sc=ddot(3*nat,rot(1,1,i),1,tras(1,1,j),1)
          CALL daxpy(3*nat,-sc,tras(1,1,j),1,rot(1,1,i),1)
       ENDDO
       DO j=1,i-1
          sc=ddot(3*nat,rot(1,1,i),1,rot(1,1,j),1)
          CALL daxpy(3*nat,-sc,rot(1,1,j),1,rot(1,1,i),1)
       ENDDO
       sc=ddot(3*nat,rot(1,1,i),1,rot(1,1,i),1)
       IF (sc.GT.1.e-10_real_8) THEN
          sc=1._real_8/SQRT(sc)
          CALL dscal(3*nat,sc,rot(1,1,i),1)
       ELSE
          CALL zeroing(rot(:,:,i))!,3*nat)
       ENDIF
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE setrotath
  ! ==================================================================
  SUBROUTINE give_scr_secder(lsecder,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lsecder
    CHARACTER(len=30)                        :: tag

    INTEGER :: lcalc_alm, lcopot, lddipo, lforces, linitrun, lortho, lrhoofr, &
      lrnlsm, lsymmat, ltddft, lupdate, nstate

    nstate=crge%n
    lcopot=0
    lortho=0
    lrnlsm=0
    lrhoofr=0
    lcalc_alm=0
    lforces=0
    lsymmat=0
    ltddft=0
    lddipo=0
    CALL give_scr_initrun(linitrun,tag)
    IF (restart1%restart)THEN
       IF (corel%tinlc) CALL give_scr_copot(lcopot,tag)
       CALL give_scr_ortho(lortho,tag,nstate)
    ENDIF
    IF (cntl%tdiag) THEN
       IF (pslo_com%tivan) CALL give_scr_rnlsm(lrnlsm,tag,nstate,.FALSE.)
       CALL give_scr_rhoofr(lrhoofr,tag)
       IF (fint1%ttrot) CALL give_scr_calc_alm(lcalc_alm,tag)
       CALL give_scr_updrho(lupdate,tag,nstate,.TRUE.,cntl%tpres)
    ELSE
       CALL give_scr_updwf(lupdate,tag,nstate,.FALSE.)
       CALL give_scr_forcedr(lforces,tag,nstate,.TRUE.,.TRUE.)
    ENDIF
    IF (cntl%tddft) CALL give_scr_lr_tddft(ltddft,.TRUE.,tag)
    IF (.NOT.restart1%rvib) CALL give_scr_ortho(lortho,tag,nstate)
    IF (.NOT.(symmi%indpg.EQ.0.OR.symmi%nrot.EQ.1))&
         CALL give_scr_symmat(lsymmat,tag)
    IF (prop1%dberry) THEN
       IF (tkpts%tkpnt) THEN
          CALL give_scr_kddipo(lddipo,nstate,tag)
       ELSE
          CALL give_scr_ddipo(lddipo,tag)
       ENDIF
    ENDIF
    lsecder=MAX(9*ions1%nat*ions1%nat+9*ions1%nat,&
         linitrun,lcopot,lortho,lrnlsm,lrhoofr,lcalc_alm,&
         lupdate,lforces,lsymmat,ltddft,lddipo)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_secder
  ! ==================================================================
  SUBROUTINE mldgau(vibe,sder,ndim,tau0)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: vibe(:), sder(:,:)
    INTEGER                                  :: ndim
    REAL(real_8)                             :: tau0(:,:,:)

    CHARACTER(len=100)                       :: file1, file2
    CHARACTER(len=17)                        :: fr, rm(5)
    CHARACTER(len=2)                         :: typ
    CHARACTER(len=7)                         :: ck(3)
    CHARACTER(len=8)                         :: a
    INTEGER                                  :: at, he, i, ia, is, iunit, j, &
                                                k, l, m, n, t(ions1%nat)
    LOGICAL                                  :: ferror, lw
    REAL(real_8)                             :: bx, by, bz, gz, &
                                                xc(ions1%nat), yc(ions1%nat), &
                                                zc(ions1%nat)

    file1='VIB1.log'
    IF (cnti%nvib.EQ.2) file2='VIB2.log'
    ! write gaussianformated output into VIB1.log and VIB2.log files,
    ! which are readable in molden/molekel to visualise the vibrations.
    fr=' Frequencies --  '
    rm(1)=' Red. masses --  '
    rm(2)=' Frc consts  --  '
    rm(3)=' IR Inten    --  '
    rm(4)=' Raman Activ --  '
    rm(5)=' Depolar     --  '
    typ='?A'
    a=' Atom AN'
    ck(1)='      X'
    ck(2)='      Y'
    ck(3)='      Z'
    at=0
    gz=0._real_8
    lw=.FALSE.
    ferror=.FALSE.
    ! ---- cell dimensions---------------------------
    bx=cell_com%celldm(1)
    by=cell_com%celldm(2)
    bz=cell_com%celldm(3)
    IF (paral%io_parent)&
         CALL fileopen(15,file1,fo_def,ferror)
    IF (ferror) GOTO 999
    IF ((cnti%nvib.EQ.2).AND.paral%io_parent)&
         CALL fileopen(17,file2,fo_def,ferror)
    ! ---- read coordinates and atomtyps ------------
    k=1
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          xc(k)=tau0(1,ia,is)
          yc(k)=tau0(2,ia,is)
          zc(k)=tau0(3,ia,is)
          t(k)=ions0%iatyp(is)
          IF (xc(k).LT.0) lw=.TRUE.
          IF (yc(k).LT.0) lw=.TRUE.
          IF (zc(k).LT.0) lw=.TRUE.
          k=k+1
       ENDDO
    ENDDO
    ! ---- write some lines needed for molden
    DO iunit=15,17,2
       IF (paral%io_parent)&
            WRITE (iunit,*)'Entering Gaussian System'
       IF (paral%io_parent)&
            WRITE (iunit,*)'this file is generated from the MLDGAU',&
            ' subroutine in the file secder.F'
       IF (paral%io_parent)&
            WRITE (iunit,*)'Please note, that this is a "faked" output;'
       IF (paral%io_parent)&
            WRITE (iunit,*)'there are no intensities computed in CPMD.'
       IF (paral%io_parent)&
            WRITE (iunit,*)'Standard orientation:'
       IF (paral%io_parent)&
            WRITE (iunit,*)'---------------------------------------',&
            '------------------------------'
       IF (paral%io_parent)&
            WRITE (iunit,'(A,2(5X,A),14X,A)')&
            'Center','Atomic','Atomic','Coordinates (Angstroms)'
       IF (paral%io_parent)&
            WRITE (iunit,'(2(A,5X),1X,A,3X,3(11X,A))')&
            'Number','Number','Type','X','Y','Z'
       IF (paral%io_parent)&
            WRITE (iunit,*)'---------------------------------------',&
            '------------------------------'
       IF (cnti%nvib.NE.2) GOTO 55
    ENDDO
55  CONTINUE
    ! ---- make sure that the center of mass is the origin or move the atoms
    ! ---- from the center of the box to the origin and printout
    IF (lw) THEN
       DO i=1,ions1%nat
          xc(i)=(xc(i))/fbohr
          yc(i)=(yc(i))/fbohr
          zc(i)=(zc(i))/fbohr
          IF (paral%io_parent)&
               WRITE (15,22)  i, t(i), at, xc(i), yc(i), zc(i)
          IF ((cnti%nvib.EQ.2).AND.paral%io_parent)&
               WRITE (17,22)  i, t(i), at, xc(i), yc(i), zc(i)
       ENDDO
    ELSE
       DO i=1,ions1%nat
          xc(i)=(xc(i)-bx/2)/fbohr
          yc(i)=(yc(i)-bx*by/2)/fbohr
          zc(i)=(zc(i)-bx*bz/2)/fbohr
          IF (paral%io_parent)&
               WRITE (15,22)  i, t(i), at, xc(i), yc(i), zc(i)
          IF ((cnti%nvib.EQ.2).AND.paral%io_parent)&
               WRITE (17,22)  i, t(i), at, xc(i), yc(i), zc(i)
       END DO
    ENDIF
    ! ---- write some lines for molden -----------------------------
    DO iunit=15,17,2
       IF (paral%io_parent)&
            WRITE(iunit,*)'--------------------------------------------',&
            '-------------------------'
       IF (paral%io_parent)&
            WRITE(iunit,*)'      basis functions          primitive ',&
            'gaussians'
       IF (paral%io_parent)&
            WRITE(iunit,*)'      alpha electrons          beta electrons'
       IF (paral%io_parent)&
            WRITE(iunit,*)'********************************************',&
            '**************************'
       IF (paral%io_parent)&
            WRITE(iunit,*)
       IF (paral%io_parent)&
            WRITE(iunit,*)'Harmonic frequencies (cm**-1), IR intensities ',&
            '(KM/Mole),'
       IF (paral%io_parent)&
            WRITE(iunit,*)'Raman scattering activities (A**4/AMU), Raman ',&
            'depolarization ratios,'
       IF (paral%io_parent)&
            WRITE(iunit,*)'reduced masses (AMU), force constants ',&
            '(mDyne/A) and normal coordinates:'
       IF (cnti%nvib.NE.2) GOTO 56
    ENDDO
56  CONTINUE
    ! ---- write eigenvalues and eigenvectors in both files
    DO i=4,(ndim),3
       IF (paral%io_parent)&
            WRITE(15,23) i-3, i-2, i-1
       IF (paral%io_parent)&
            WRITE(15,24) typ, typ, typ
       IF (paral%io_parent)&
            WRITE(15,25) fr, (vibe(l),l=i,i+2)
       DO n=1,5
          IF (paral%io_parent)&
               WRITE(15,25) rm(n), gz, gz, gz
       ENDDO
       IF (paral%io_parent)&
            WRITE(15,26) a,(ck(n),n=1,3),(ck(n),n=1,3),(ck(n),n=1,3)
       DO j=1,ndim,3
          he=(j-1)/3+1
          IF (paral%io_parent)&
               WRITE(15,27) he,t(he),&
               (sder(j,m),sder(j+1,m),sder(j+2,m),m=i,i+2)
       END DO
    END DO
    IF (paral%io_parent)&
         WRITE(15,*) 'Normal termination of Gaussian 98.'
    IF (cnti%nvib.EQ.2) THEN
       DO i=1,(ndim-3),3
          IF (paral%io_parent)&
               WRITE(17,23) i, i+1, i+2
          IF (paral%io_parent)&
               WRITE(17,24) typ, typ, typ
          IF (paral%io_parent)&
               WRITE(17,25) fr, (vibe(l),l=i,i+2)
          DO n=1,5
             IF (paral%io_parent)&
                  WRITE(17,25) rm(n), gz, gz, gz
          ENDDO
          IF (paral%io_parent)&
               WRITE(17,26) a,(ck(n),n=1,3),(ck(n),n=1,3),(ck(n),n=1,3)
          DO j=1,ndim,3
             he=(j-1)/3+1
             IF (paral%io_parent)&
                  WRITE(17,27) he,t(he),&
                  (sder(j,m),sder(j+1,m),sder(j+2,m),m=i,i+2)
          END DO
       END DO
       IF (paral%io_parent)&
            WRITE(17,*) 'Normal termination of Gaussian 98.'
    ENDIF
    IF (paral%io_parent)&
         CALL fileclose(15)
    IF ((cnti%nvib.EQ.2).AND.paral%io_parent)&
         CALL fileclose(17)
22  FORMAT(i5,i11,i14,4x,3(3x,f9.6))
23  FORMAT(i22,2i23)
24  FORMAT(20x,a2,2(21x,a2))
25  FORMAT(a17,f9.4,2f23.4)
26  FORMAT(a8,3a7,2(2x,3a7))
27  FORMAT(2i4,3(f9.2,2f7.2))
    GOTO 888
999 CONTINUE
    IF (paral%io_parent)&
         WRITE(6,*) 'COULD NOT OPEN FILE VIB[1,2].log! '
888 CONTINUE
    ! ==--------------------------------------------------------------==
    RETURN
  END  SUBROUTINE mldgau
  ! ==================================================================
  SUBROUTINE writeaclimax(vibe,sder,ndim,tau0)
    ! ==--------------------------------------------------------------==
    IMPLICIT NONE
    ! Arguments
    INTEGER :: ndim
    REAL(real_8) :: vibe(:),sder(:,:),tau0(:,:,:)
    ! Variables
    LOGICAL :: ferror
    CHARACTER (len=100) :: file
    INTEGER :: i,j,k,l,n,t(ions1%nat),is,ia,he,m
    ! ==--------------------------------------------------------------==
    file='VIB.aclimax'
    ! write aClimax formated output into cpmd.aclimax file
    ! which are readable in aClimax to generate INS spectrum
    ! INS is inelastic neutron scattering
    ! convert to windows format befoe using aClimax for windows
    ferror=.FALSE.
    IF (paral%io_parent)&
         CALL fileopen(17,file,fo_def,ferror)
    IF (ferror) GOTO 999
    ! ---- read atomtypes ------------
    k=1
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          t(k)=ions0%iatyp(is)
          k=k+1
       ENDDO
    ENDDO
    ! ---- write aClimax comment header
    IF (paral%io_parent)&
         WRITE(17,*)'! This file is generated from the WRITEACLIMAX',&
         ' subroutine in the file secder.F'
    IF (paral%io_parent)&
         WRITE(17,*)'! #  Type  Number      Mass      ',&
         'Coord        Coord        Coord'
    IF (paral%io_parent)&
         WRITE(17,*)
    IF (paral%io_parent)&
         WRITE(17,*)'BEGIN ATOMS'
    i=1
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          IF (paral%io_parent)&
               WRITE(17,22) i,elem%el(ions0%iatyp(is)),ions0%iatyp(is),rmass%pma0(is),&
               (tau0(k,ia,is),k=1,3)
       ENDDO
    ENDDO
    IF (paral%io_parent)&
         WRITE(17,*)'END ATOMS'
    IF (paral%io_parent)&
         WRITE(17,*)
    ! ---- write eigenvalues and eigenvectors in both files
    IF (paral%io_parent)&
         WRITE(17,*)'BEGIN FREQUENCIES'
    DO i=4,ndim
       IF (paral%io_parent)&
            WRITE(17,*)'! Mode#    Frequency'
       IF (paral%io_parent)&
            WRITE(17,'(I7,4X,f7.2)') (i-3),(vibe(i))
       IF (paral%io_parent)&
            WRITE(17,*)
       IF (paral%io_parent)&
            WRITE(17,*)'! Displacments'
       IF (paral%io_parent)&
            WRITE(17,*)
       DO j=1,ndim,3
          he=(j-1)/3+1
          IF (paral%io_parent)&
               WRITE(17,27) he,elem%el(t(he)),t(he),&
               (sder(m,i),m=j,j+2)
       ENDDO
       IF (paral%io_parent)&
            WRITE(17,*)
    ENDDO
    IF (paral%io_parent)&
         WRITE(17,*)'END FREQUENCIES'
    ! ---- close file and exit subroutine
    IF (paral%io_parent)&
         CALL fileclose(17)
22  FORMAT(i3,a4,i6,f15.6,3f15.8)
27  FORMAT(i4,a4,i4,3f12.6)
    GOTO 888
999 CONTINUE
    IF (paral%io_parent)&
         WRITE(6,*) 'COULD NOT OPEN FILE cpmd.aclimax! '
888 CONTINUE
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE writeaclimax
  ! ==================================================================

END MODULE secder_utils
