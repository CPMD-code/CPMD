MODULE stcop_utils
  USE dotp_utils,                      ONLY: dotp
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE hfx_drivers,                     ONLY: hfxpsi
  USE kinds,                           ONLY: real_8
  USE linres,                          ONLY: &
       lrsym, nlinr, nlinw, nua, nub, td01, td02, td03, urot
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum
  USE parac,                           ONLY: parai,&
                                             paral
  USE prng_utils,                      ONLY: repprngu_vec_cmplx
  USE rho1ofr_utils,                   ONLY: rho1ofr
  USE spin,                            ONLY: spin_mod
  USE system,                          ONLY: cntl,&
                                             ncpw
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: zclean
  USE v1ofrho1_utils,                  ONLY: vtdofrho1
  USE vpsi_utils,                      ONLY: vpsimt
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: stcop
  !public :: stcop1
  !public :: edia
  !public :: stcop2

CONTAINS

  ! ==================================================================
  SUBROUTINE stcop(ns,c0,c2,nstate,eigv,edav,cv,nva,nvb,c1,nact,&
       ddxc,psi,rhoe,orbital)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ns
    COMPLEX(real_8)                          :: c0(ncpw%ngw,*), c2(ncpw%ngw,*)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: eigv(*), edav(*)
    COMPLEX(real_8)                          :: cv(ncpw%ngw,*)
    INTEGER                                  :: nva, nvb, nact
    COMPLEX(real_8)                          :: c1(ncpw%ngw,nact,nlinw,*)
    REAL(real_8)                             :: ddxc(:,:)
    COMPLEX(real_8)                          :: psi(:,:)
    REAL(real_8)                             :: rhoe(:,:)
    CHARACTER(len=*)                         :: orbital

    CHARACTER(*), PARAMETER                  :: procedureN = 'stcop'

    INTEGER                                  :: isub

! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    IF (INDEX(orbital,"CANON").NE.0) THEN
       CALL stcop1(ns,c0,c2,nstate,eigv,edav,cv,nva,nvb,c1,nact,ddxc,&
            psi,rhoe)
    ELSE
       CALL stcop2(ns,c0,c2,nstate,eigv,edav,cv,nva,nvb,c1,nact,ddxc,&
            psi,rhoe)
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    RETURN
  END SUBROUTINE stcop
  ! ==================================================================
  SUBROUTINE stcop1(ns,c0,c2,nstate,eigv,edav,cv,nva,nvb,c1,nact,&
       ddxc,psi,rhoe)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ns
    COMPLEX(real_8)                          :: c0(ncpw%ngw,*), c2(ncpw%ngw,*)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: eigv(*), edav(*)
    COMPLEX(real_8)                          :: cv(ncpw%ngw,*)
    INTEGER                                  :: nva, nvb, nact
    COMPLEX(real_8)                          :: c1(ncpw%ngw,nact,nlinw,*)
    REAL(real_8)                             :: ddxc(:,:)
    COMPLEX(real_8)                          :: psi(:,:)
    REAL(real_8)                             :: rhoe(:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'stcop1'

    COMPLEX(real_8), ALLOCATABLE             :: cs(:,:), ct(:,:)
    COMPLEX(real_8), ALLOCATABLE, &
      DIMENSION(:)                           :: zrandom
    INTEGER                                  :: i, ierr, ij, is, is1, is2, &
                                                isub, ix, ix1, ix2, iy, iy1, &
                                                iy2, j, laux, no, nov, nv
    LOGICAL                                  :: debug, ostda
    REAL(real_8)                             :: dev, devmax, e2(5), focc(1)
    REAL(real_8), ALLOCATABLE                :: amat(:,:), aux(:), deig(:)

! ==--------------------------------------------------------------==
! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)

    debug=.FALSE.
    IF (td03%tda.AND.td01%msubta.GT.0.AND..NOT.td03%tdacanon) THEN
       ostda=.TRUE.
    ELSE
       ostda=.FALSE.
    ENDIF
    nv=nva+nvb
    IF (cntl%tlsd) THEN
       IF (paral%io_parent) WRITE(6,'(A,T56,I10)')&
            " NUMBER OF STATES TO BE INITIALIZED ",ns-nlinr
       DO is=nlinr+1,ns,2
          is1=is
          is2=is+1
          CALL zeroing(c1(:,:,is1,1))!,ngw*nact)
          IF (.NOT.td03%tda) CALL zeroing(c1(:,:,is1,2))!,ngw*nact)
          IF (is2.LE.ns) THEN
             CALL zeroing(c1(:,:,is2,1))!,ngw*nact)
             IF (.NOT.td03%tda) CALL zeroing(c1(:,:,is2,2))!,ngw*nact)
          ENDIF
          iy1=MAX(nact/2-(is-1)/nva,1)
          iy2=spin_mod%nsup+MAX(nact/2-(is-1)/nvb,1)
          ix1=MOD(is,nva)
          IF (ix1.LE.0) ix1=ix1+nva
          ix1=MIN(ix1,nva)
          ix2=MOD(is,nvb)
          IF (ix2.LE.0) ix2=ix2+nvb
          ix2=nva+MIN(ix2,nvb)
          CALL dcopy(2*ncpw%ngw,cv(1,ix1),1,c1(1,iy1,is1,1),1)
          CALL dcopy(2*ncpw%ngw,cv(1,ix2),1,c1(1,iy2,is1,1),1)
          IF (is2.LE.ns) THEN
             CALL dcopy(2*ncpw%ngw,cv(1,ix1),1,c1(1,iy1,is2,1),1)
             CALL dcopy(2*ncpw%ngw,cv(1,ix2),1,c1(1,iy2,is2,1),1)
             CALL dscal(2*ncpw%ngw,-1._real_8,c1(1,iy2,is2,1),1)
          ENDIF
          IF (td02%tdrandom.GT.0._real_8) THEN
             ALLOCATE(zrandom(ncpw%ngw),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
                  __LINE__,__FILE__)
             DO iy=1,nact
                !call vrand(2*ngw,c1(1,iy,is1,1),tdrandom)
                CALL repprngu_vec_cmplx(ncpw%ngw,zrandom)
                c1(:,iy,is1,1)=c1(:,iy,is1,1)+zrandom(:)*td02%tdrandom
                !if(.not.tda) call vrand(2*ngw,c1(1,iy,is1,2),tdrandom)
                IF(.NOT.td03%tda) THEN
                   CALL repprngu_vec_cmplx(ncpw%ngw,zrandom)
                   c1(:,iy,is1,1)=c1(:,iy,is1,1)+zrandom(:)*td02%tdrandom
                ENDIF
             ENDDO
             CALL zclean(c1(1,1,is1,1),nact,ncpw%ngw)
             IF (.NOT.td03%tda) CALL zclean(c1(1,1,is1,2),nact,ncpw%ngw)
             IF (is2.LE.ns) THEN
                DO iy=1,nact
                   !call vrand(2*ngw,c1(1,iy,is2,1),tdrandom)
                   CALL repprngu_vec_cmplx(ncpw%ngw,zrandom)
                   c1(:,iy,is2,1)=c1(:,iy,is2,1)+zrandom(:)*td02%tdrandom
                   !if(.not.tda) call vrand(2*ngw,c1(1,iy,is2,2),tdrandom)
                   IF(.NOT.td03%tda) THEN
                      CALL repprngu_vec_cmplx(ncpw%ngw,zrandom)
                      c1(:,iy,is2,1)=c1(:,iy,is2,1)+zrandom(:)*td02%tdrandom
                   ENDIF
                ENDDO
                CALL zclean(c1(1,1,is2,1),nact,ncpw%ngw)
                IF (.NOT.td03%tda) CALL zclean(c1(1,1,is2,2),nact,ncpw%ngw)
             ENDIF
             DEALLOCATE(zrandom,STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
                  __LINE__,__FILE__)
          ENDIF
       ENDDO
    ELSE
       IF (nlinr.LT.ns) THEN
          ! initialize states
          IF (paral%io_parent) WRITE(6,'(A,T56,I10)')&
               " NUMBER OF STATES TO BE INITIALIZED ",ns-nlinr
          dev=eigv(nstate+nv)-eigv(nstate+1)
          no=0
          DO i=1,nstate
             IF (eigv(i).GT.eigv(nstate)-dev) no=no+1
          ENDDO
          no=MIN(no,nv)
          no=MAX(no,2)
          no=MIN(no,nstate)
          IF (td03%tdacanon) no=MIN(no,nact)
          nov=no*nv
          IF (paral%io_parent) WRITE(6,'(A,T56,I10)')&
               " TOTAL NUMBER OF TEST VECTORS ",nov
          ALLOCATE(ct(ncpw%ngw,nov),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          CALL zeroing(ct)
          ALLOCATE(cs(ncpw%ngw,nov),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          CALL zeroing(cs)

          DO i=1,no
             DO j=1,nv
                ij=(i-1)*nv+j
                CALL dcopy(2*ncpw%ngw,cv(1,j),1,ct(1,ij),1)
                ix=nstate+1-i
                CALL dcopy(2*ncpw%ngw,cv(1,j),1,cs(1,ij),1)
                dev=eigv(nstate+j)-eigv(ix)
                CALL dscal(2*ncpw%ngw,-dev,cs(1,ij),1)
                IF (lrsym.EQ.1) THEN
                   focc=1._real_8
                   CALL rho1ofr(c0(1,ix),ct(1,ij),focc,rhoe,psi(:,1),1)
                   CALL vtdofrho1(e2,rhoe,ddxc,psi,.TRUE.)
                   CALL vpsimt(c0(:,ix:ix),cs(:,ij:ij),focc,rhoe,psi(:,1),1,1,.TRUE.)
                   CALL hfxpsi(c0(:,1:nstate),c0(:,ix:ix),cs(:,ij:ij),crge%f(:,1),-1._real_8,psi(:,1),nstate,1)
                ENDIF
             ENDDO
          ENDDO

          ALLOCATE(amat(nov,nov),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
          ALLOCATE(deig(nov),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)

          laux=20*nov
          ALLOCATE(aux(laux),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)

          CALL edia(ct,cs,no,nv,amat,deig,aux,laux)

          DO is=nlinr+1,ns
             CALL zeroing(c1(:,:,is,1))!,ngw*nact)
             IF (.NOT.td03%tda) CALL zeroing(c1(:,:,is,2))!,ngw*nact)
             devmax=0._real_8
             ix1=0
             DO i=1,nov
                ix=(is-1)*nov+i
                !dev=amat(ix)
                dev=amat(i,is)
                IF (ABS(dev).GT.devmax) THEN
                   devmax=dev
                   ix1=i
                ENDIF
                IF (ABS(dev).GT.1.e-12_real_8) THEN
                   iy=nstate-(i-1)/nv
                   iy=MOD(iy-1,nact)+1
                   CALL daxpy(2*ncpw%ngw,dev,ct(1,i),1,c1(1,iy,is,1),1)
                ENDIF
             ENDDO
             IF (td02%tdrandom.GT.0._real_8) THEN
                ALLOCATE(zrandom(ncpw%ngw),STAT=ierr)
                IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
                     __LINE__,__FILE__)
                DO iy=1,nact
                   !call vrand(2*ngw,c1(1,iy,is,1),tdrandom)
                   CALL repprngu_vec_cmplx(ncpw%ngw,zrandom)
                   c1(:,iy,is,1)=c1(:,iy,is,1)+zrandom(:)*td02%tdrandom
                   !if(.not.tda) call vrand(2*ngw,c1(1,iy,is,2),tdrandom)
                   IF(.NOT.td03%tda) THEN
                      CALL repprngu_vec_cmplx(ncpw%ngw,zrandom)
                      c1(:,iy,is,2)=c1(:,iy,is,2)+zrandom(:)*td02%tdrandom
                   ENDIF
                   IF(ostda) THEN
                      !call vrand(nua-nub,urot(1,iy,is),10._real_8*tdrandom)
                      CALL repprngu_vec_cmplx(nua-nub,zrandom)
                      urot(1:nua-nub,iy,is)=urot(1:nua-nub,iy,is) &
                           +zrandom(1:nua-nub)*10._real_8*td02%tdrandom
                      DO i=1,nua-nub
                         urot(i,iy,is)=urot(i,iy,is)*REAL(i,kind=real_8)/REAL(nua,kind=real_8)
                      ENDDO
                   ENDIF
                ENDDO
                CALL zclean(c1(1,1,is,1),nact,ncpw%ngw)
                IF (.NOT.td03%tda) CALL zclean(c1(1,1,is,2),nact,ncpw%ngw)
                IF (ostda) THEN
                   CALL mp_bcast(urot(:,:,is),nua*nub,parai%source,parai%allgrp)
                ENDIF
                DEALLOCATE(zrandom,STAT=ierr)
                IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
                     __LINE__,__FILE__)
             ENDIF
          ENDDO
          DEALLOCATE(amat,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
               __LINE__,__FILE__)
          DEALLOCATE(deig,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
               __LINE__,__FILE__)
          DEALLOCATE(aux,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
               __LINE__,__FILE__)
          DEALLOCATE(ct,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
          DEALLOCATE(cs,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    RETURN
  END SUBROUTINE stcop1
  ! ==================================================================
  SUBROUTINE edia(ct,cs,no,nv,amat,deig,aux,laux)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: ct(:,:), cs(:,:)
    INTEGER                                  :: no, nv
    REAL(real_8)                             :: amat(:,:), deig(:), aux(:)
    INTEGER                                  :: laux

    INTEGER                                  :: i, info, io, j, jo, nov

    nov=no*nv
    CALL zeroing(amat)!,nov*nov)
    DO i=1,nov
       io=i/nv+1
       DO j=1,i
          jo=j/nv+1
          IF (io.EQ.jo) THEN
             amat(i,j) = -dotp(ncpw%ngw,cs(:,i),ct(:,j))
             amat(j,i) = amat(i,j)
          ENDIF
       ENDDO
    ENDDO
    CALL mp_sum(amat,nov*nov,parai%allgrp)
    IF (paral%parent) THEN
       CALL dsyev('V','L',nov,amat,nov,deig,aux,laux,info)
       IF (info.NE.0) CALL stopgm('EDIA','INFO ! = 0 AFTER DSYEV',& 
            __LINE__,__FILE__)
    ENDIF
    CALL mp_bcast(amat,nov*nov,parai%source,parai%allgrp)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE edia
  ! ==================================================================
  SUBROUTINE stcop2(ns,c0,c2,nstate,eigv,edav,cv,nva,nvb,c1,nact,&
       ddxc,psi,rhoe)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ns
    COMPLEX(real_8)                          :: c0(ncpw%ngw,*), c2(ncpw%ngw,*)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: eigv(*), edav(*)
    COMPLEX(real_8)                          :: cv(ncpw%ngw,*)
    INTEGER                                  :: nva, nvb, nact
    COMPLEX(real_8)                          :: c1(ncpw%ngw,nact,nlinw,*)
    REAL(real_8)                             :: ddxc(:,:)
    COMPLEX(real_8)                          :: psi(:,:)
    REAL(real_8)                             :: rhoe(:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'stcop2'

    COMPLEX(real_8), ALLOCATABLE, &
      DIMENSION(:)                           :: zrandom
    INTEGER                                  :: i, ierr, is, is1, is2, ix, &
                                                ix1, ix2, iy, iy1, iy2, nv
    LOGICAL                                  :: debug, ostda

! ==--------------------------------------------------------------==
! ==--------------------------------------------------------------==

    debug=.FALSE.
    IF (td03%tda.AND.td01%msubta.GT.0.AND..NOT.td03%tdacanon) THEN
       ostda=.TRUE.
    ELSE
       ostda=.FALSE.
    ENDIF
    nv=nva+nvb
    IF (cntl%tlsd) THEN
       IF (paral%io_parent)&
            WRITE(6,'(A,T56,I10)')&
            " NUMBER OF STATES TO BE INITIALIZED ",ns-nlinr
       DO is=nlinr+1,ns,2
          is1=is
          is2=is+1
          CALL zeroing(c1(:,:,is1,1))!,ngw*nact)
          IF (.NOT.td03%tda) CALL zeroing(c1(:,:,is1,2))!,ngw*nact)
          IF (is2.LE.ns) THEN
             CALL zeroing(c1(:,:,is2,1))!,ngw*nact)
             IF (.NOT.td03%tda) CALL zeroing(c1(:,:,is2,2))!,ngw*nact)
          ENDIF
          iy1=MAX(nact/2-(is-1)/nva,1)
          iy2=spin_mod%nsup+MAX(nact/2-(is-1)/nvb,1)
          ix1=MOD(is,nva)
          IF (ix1.LE.0) ix1=ix1+nva
          ix1=MIN(ix1,nva)
          ix2=MOD(is,nvb)
          IF (ix2.LE.0) ix2=ix2+nvb
          ix2=nva+MIN(ix2,nvb)
          CALL dcopy(2*ncpw%ngw,cv(1,ix1),1,c1(1,iy1,is1,1),1)
          CALL dcopy(2*ncpw%ngw,cv(1,ix2),1,c1(1,iy2,is1,1),1)
          IF (is2.LE.ns) THEN
             CALL dcopy(2*ncpw%ngw,cv(1,ix1),1,c1(1,iy1,is2,1),1)
             CALL dcopy(2*ncpw%ngw,cv(1,ix2),1,c1(1,iy2,is2,1),1)
             CALL dscal(2*ncpw%ngw,-1._real_8,c1(1,iy2,is2,1),1)
          ENDIF
          IF (td02%tdrandom.GT.0._real_8) THEN
             ALLOCATE(zrandom(ncpw%ngw),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
                  __LINE__,__FILE__)
             DO iy=1,nact
                !call vrand(2*ngw,c1(1,iy,is1,1),tdrandom)
                CALL repprngu_vec_cmplx(ncpw%ngw,zrandom)
                c1(:,iy,is1,1)=c1(:,iy,is1,1)+zrandom(:)*td02%tdrandom
                !if(.not.tda) call vrand(2*ngw,c1(1,iy,is1,2),tdrandom)
                IF(.NOT.td03%tda) THEN
                   CALL repprngu_vec_cmplx(ncpw%ngw,zrandom)
                   c1(:,iy,is1,2)=c1(:,iy,is1,2)+zrandom(:)*td02%tdrandom
                ENDIF
             ENDDO
             CALL zclean(c1(1,1,is1,1),nact,ncpw%ngw)
             IF (.NOT.td03%tda) CALL zclean(c1(1,1,is1,2),nact,ncpw%ngw)
             IF (is2.LE.ns) THEN
                DO iy=1,nact
                   !call vrand(2*ngw,c1(1,iy,is2,1),tdrandom)
                   CALL repprngu_vec_cmplx(ncpw%ngw,zrandom)
                   c1(:,iy,is2,1)=c1(:,iy,is2,1)+zrandom(:)*td02%tdrandom
                   !if(.not.tda) call vrand(2*ngw,c1(1,iy,is2,2),tdrandom)
                   IF(.NOT.td03%tda) THEN
                      CALL repprngu_vec_cmplx(ncpw%ngw,zrandom)
                      c1(:,iy,is2,2)=c1(:,iy,is2,2)+zrandom(:)*td02%tdrandom
                   ENDIF
                ENDDO
                CALL zclean(c1(1,1,is2,1),nact,ncpw%ngw)
                IF (.NOT.td03%tda) CALL zclean(c1(1,1,is2,2),nact,ncpw%ngw)
             ENDIF
             DEALLOCATE(zrandom,STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
                  __LINE__,__FILE__)
          ENDIF
       ENDDO
    ELSEIF (nlinr.LT.ns) THEN
       ! initialize states
       IF (paral%io_parent)&
            WRITE(6,'(A,T56,I10)')&
            " NUMBER OF STATES TO BE INITIALIZED ",ns-nlinr
       DO is=nlinr+1,ns
          CALL zeroing(c1(:,:,is,1))!,ngw*nact)
          IF (.NOT.td03%tda) CALL zeroing(c1(:,:,is,2))!,ngw*nact)
          iy=MAX(nact-(is-1)/nv,1)
          ix=MOD(is,nv)
          IF (ix.LE.0) ix=ix+nv
          ix=MIN(ix,nv)
          CALL dcopy(2*ncpw%ngw,cv(1,ix),1,c1(1,iy,is,1),1)
          IF (td02%tdrandom.GT.0._real_8) THEN
             ALLOCATE(zrandom(ncpw%ngw),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
                  __LINE__,__FILE__)
             DO iy=1,nact
                !call vrand(2*ngw,c1(1,iy,is,1),tdrandom)
                CALL repprngu_vec_cmplx(ncpw%ngw,zrandom)
                c1(:,iy,is,1)=c1(:,iy,is,1)+zrandom(:)*td02%tdrandom
                !if(.not.tda) call vrand(2*ngw,c1(1,iy,is,2),tdrandom)
                IF(.NOT.td03%tda) THEN
                   CALL repprngu_vec_cmplx(ncpw%ngw,zrandom)
                   c1(:,iy,is,2)=c1(:,iy,is,2)+zrandom(:)*td02%tdrandom
                ENDIF
                IF(ostda) THEN
                   !call vrand(nua-nub,urot(1,iy,is),10._real_8*tdrandom)
                   CALL repprngu_vec_cmplx(nua-nub,zrandom)
                   urot(1:nua-nub,iy,is)=urot(1:nua-nub,iy,is)+ &
                        zrandom(:)*10._real_8*td02%tdrandom
                   DO i=1,nua-nub
                      urot(i,iy,is)=urot(i,iy,is)*REAL(i,kind=real_8)/REAL(nua,kind=real_8)
                   ENDDO
                ENDIF
             ENDDO
             CALL zclean(c1(1,1,is,1),nact,ncpw%ngw)
             IF (.NOT.td03%tda) CALL zclean(c1(1,1,is,2),nact,ncpw%ngw)
             IF (ostda) THEN
                CALL mp_bcast(urot(:,:,is),nua*nub,parai%source,parai%allgrp)
             ENDIF
             DEALLOCATE(zrandom,STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
                  __LINE__,__FILE__)
          ENDIF
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE stcop2
  ! ==================================================================

END MODULE stcop_utils
