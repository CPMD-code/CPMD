MODULE deort_utils
  USE dotp_utils,                      ONLY: dotp
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: dspevy

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: deort
  PUBLIC :: give_scr_deort

CONTAINS

  ! ==================================================================
  SUBROUTINE deort(ngw,nstate,rmat,reig,c0,sc0)
    ! ==--------------------------------------------------------------==
    ! ==    DEORTHOGONALIZE WAVEFUNCTIONS FOR VANDERBILT PP           ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ngw, nstate
    REAL(real_8)                             :: rmat(*), reig(*)
    COMPLEX(real_8)                          :: c0(ngw,*), sc0(ngw,*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'deort'

    INTEGER                                  :: i, ierr, ij, iopt, isub, j, k
    REAL(real_8)                             :: serr
    REAL(real_8), ALLOCATABLE                :: work(:), z(:), zz(:)

    CALL tiset('     DEORT',isub)
    serr=0.0_real_8
    k=0
    DO i=1,nstate
       DO j=i,nstate
          k=k+1
          rmat(k)=dotp(ngw,c0(:,i),c0(:,j))
       ENDDO
    ENDDO
    CALL mp_sum(rmat,k,parai%allgrp)
    DO i=1,k
       serr=serr+rmat(i)
    ENDDO
    serr=serr-REAL(nstate,kind=real_8)
    IF (ABS(serr).GT.1.e-8_real_8) THEN
       iopt=1
       ALLOCATE(z(nstate*nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(zz(nstate*nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(work(3*nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       CALL dspevy(iopt,rmat,reig,z,nstate,nstate,work,3*nstate)
       DO i=1,nstate
          reig(i)=1.0_real_8/SQRT(reig(i))
       ENDDO
       ij=1
       DO i=1,nstate
          DO j=1,nstate
             zz(ij)=reig(i)*z(ij)
             ij=ij+1
          ENDDO
       ENDDO
       CALL dgemm('N','T',nstate,nstate,nstate,1.0_real_8,zz,nstate,&
            z,nstate,0.0_real_8,rmat,nstate)
       CALL dgemm('N','N',2*ngw,nstate,nstate,1.0_real_8,c0,2*ngw,&
            rmat,nstate,0.0_real_8,sc0,2*ngw)
       CALL dcopy(2*ngw*nstate,sc0,1,c0,1)
       DEALLOCATE(z,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(work,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(zz,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    ENDIF
    CALL tihalt('     DEORT',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE deort
  ! ==================================================================
  SUBROUTINE give_scr_deort(ldeort,tag,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ldeort
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate

! ==--------------------------------------------------------------==

    ldeort=nstate*nstate+3*nstate
    tag  ='NSTATE*NSTATE+3*NSTATE'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_deort
  ! ==================================================================

END MODULE deort_utils
