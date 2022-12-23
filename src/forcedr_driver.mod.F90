MODULE forcedr_driver
  USE bsym,                            ONLY: bsclcs
  USE cnstfc_utils,                    ONLY: restfc
  USE cotr,                            ONLY: cotr007,&
                                             gsrate,&
                                             resval,&
                                             resval_dest
  USE error_handling,                  ONLY: stopgm
  USE forces_driver,                   ONLY: forces
  USE kinds,                           ONLY: real_8
  USE mm_dim_utils,                    ONLY: mm_dim
  USE mm_dimmod,                       ONLY: mm_go_mm,&
                                             mm_go_qm,&
                                             mm_revert
  USE mp_interface,                    ONLY: mp_bcast
  USE noforce_utils,                   ONLY: noforce
  USE parac,                           ONLY: parai,&
                                             paral
  USE pimd,                            ONLY: ipcurr,&
                                             np_local,&
                                             np_low
  USE reshaper,                        ONLY: reshape_inplace
  USE system,                          ONLY: cntl,&
                                             maxsys,&
                                             ncpw
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE tpar,                            ONLY: dt_ions
  USE utils,                           ONLY: unitmx
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: forcedr

CONTAINS

  ! ==================================================================
  SUBROUTINE forcedr(c0,c2,sc0,rhoe,psi,tau0,fion,eigv,&
       nstate,nkpoint,lproj,tfor)
    ! ==--------------------------------------------------------------==
    ! ==             DRIVER FOR FORCE CALCULATION                     ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8), TARGET                  :: c0(:,:)
    COMPLEX(real_8)                          :: c2(:,:)
    REAL(real_8)                             :: rhoe(:,:)
    COMPLEX(real_8)                          :: psi(:,:)
    REAL(real_8)                             :: tau0(:,:,:), fion(:,:,:), &
                                                eigv(*)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: sc0(ncpw%ngw,nstate)
    INTEGER                                  :: nkpoint
    LOGICAL                                  :: lproj, tfor

    CHARACTER(*), PARAMETER                  :: procedureN = 'forcedr'

    COMPLEX(real_8), POINTER                 :: c0_ptr(:,:,:)
    INTEGER                                  :: i, ierr, ipx, isub, numx
    INTEGER, SAVE                            :: icon = 0
    LOGICAL                                  :: oldstatus, statusdummy
    REAL(real_8), ALLOCATABLE, SAVE          :: xmat1(:,:), xmat2(:)

    CALL tiset(procedureN,isub)
    CALL mm_dim(mm_go_mm,oldstatus)
    CALL zeroing(fion)!,3*maxsys%nax*maxsys%nsx)
    CALL mm_dim(mm_go_qm,statusdummy)
    IF (icon.EQ.0) THEN
       IF (cntl%nonort) THEN
          IF (paral%parent) THEN
             IF (cntl%tpath) THEN
                numx=np_local
             ELSE
                numx=1
                ! For cntl%bsymm: we have two wavefunctions
                IF (cntl%bsymm)numx=2
             ENDIF
             ALLOCATE(xmat1(nstate*nstate,numx),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             ALLOCATE(xmat2(nstate*nstate),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             DO i=1,numx
                CALL unitmx(xmat1(1,i),nstate)
             ENDDO
             icon=1
          ELSE
             ! To avoid Fortran runtime 'not allocated' error
             IF(ALLOCATED(xmat1)) THEN
                DEALLOCATE(xmat1,STAT=ierr)
                IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
                     __LINE__,__FILE__)
             ENDIF
             ALLOCATE(xmat1(1,1),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)! TODO XMAT1 is not used in NOFORCE
             IF(ALLOCATED(xmat2)) THEN
                DEALLOCATE(xmat2,STAT=ierr)
                IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
                     __LINE__,__FILE__)
             ENDIF
             ALLOCATE(xmat2(1),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
          ENDIF
       ENDIF
    ENDIF
    IF (cntl%nonort) THEN
       IF (cntl%tpath) THEN
          ipx=ipcurr-np_low+1
          ! For cntl%bsymm: BSCLCS=1 for BS and BSCLCS==2 for HS
       ELSE IF (cntl%bsymm) THEN
          ipx=bsclcs
       ELSE
          ipx=1
       ENDIF
       CALL noforce(c0,c2,sc0,tau0,fion,eigv,rhoe,psi,&
            xmat1(1,ipx),xmat2,nstate,tfor)
    ELSE
       NULLIFY(c0_ptr)
       CALL reshape_inplace(c0, (/SIZE(c0,1),SIZE(c0,2),1/), c0_ptr)
       CALL forces(c0_ptr,c2,tau0,fion,rhoe,psi,&
            nstate,nkpoint,lproj,tfor)
    ENDIF
    ! Forces from restraints (may be in the MM part as well)
    IF ( tfor ) THEN
       CALL mm_dim(mm_go_mm,statusdummy)
       IF (paral%parent) THEN
          DO i=1,cotr007%mrestr
             resval(i)=resval(i)+gsrate(i)*dt_ions
             IF (resval_dest(i).NE.-999._real_8) THEN
                IF (gsrate(i).GT.0._real_8.AND.resval(i).GT.resval_dest(i))&
                     resval(i)=resval_dest(i) ! increase
                IF (gsrate(i).LT.0._real_8.AND.resval(i).LT.resval_dest(i))&
                     resval(i)=resval_dest(i) ! decrease
             ENDIF
          ENDDO
          CALL restfc(tau0,fion)
       ENDIF
       !CALL mp_bcast(fion,3*maxsys%nax*maxsys%nsx,parai%source,parai%allgrp)
       CALL mp_bcast(fion,3*maxsys%nax*maxsys%nsx,parai%io_source,parai%cp_grp)
    ENDIF
    CALL mm_dim(mm_revert,oldstatus)
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE forcedr

END MODULE forcedr_driver
