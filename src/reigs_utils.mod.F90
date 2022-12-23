MODULE reigs_utils
  USE cnst,                            ONLY: ry
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: paral
  USE ropt,                            ONLY: ropt_mod
  USE spin,                            ONLY: spin_mod
  USE system,                          ONLY: cntl
  USE timer,                           ONLY: tiset
  USE utils,                           ONLY: dspevy

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: reigs

CONTAINS

  ! ==================================================================
  SUBROUTINE reigs(nstate,gam,f)
    ! ==--------------------------------------------------------------==
    ! ==  COMPUTES THE EIGENVALUES OF THE REAL SYMMETRIC MATRIX GAM   ==
    ! ==  AND UPDATES THE DENSITY OF STATES DOFS ACCORDINGLY.         ==
    ! ==  THE EIGENVALUES ARE PRINTED OUT IN ELECTRON VOLTS.          ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    REAL(real_8)                             :: gam(nstate,nstate), f(nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'reigs'

    INTEGER                                  :: i, ierr, iopt, j, k, nleft
    REAL(real_8)                             :: dummy(1), ffi
    REAL(real_8), ALLOCATABLE                :: ap(:), fm1(:), wr(:)

! Variables
! ==--------------------------------------------------------------==
! vw      CALL TISET('     REIGS',ISUB)

    ALLOCATE(fm1(3*nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(wr(nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(ap(nstate*(nstate+1)/2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    IF (.NOT.cntl%tlsd) THEN
       k=0
       DO j=1,nstate
          DO i=j,nstate
             k=k+1
             ap(k)=gam(i,j)
          ENDDO
       ENDDO
       iopt=0
       CALL dspevy(iopt,ap,wr,dummy,nstate,nstate,fm1,3*nstate)
       DO i=1,nstate
          ffi=f(i)
          IF (ffi.LT.1.e-5_real_8) ffi=1._real_8
          wr(i)=wr(i)/ffi
       ENDDO
       IF (ropt_mod%engpri) THEN
          IF (paral%io_parent)&
               WRITE(6,'(/," EIGENVALUES (EV) AND OCCUPATION:")')
          DO i=1,nstate,2
             nleft=MIN(nstate-i,1)
             IF (paral%io_parent)&
                  WRITE(6,'(3X,2(I3,F15.7,4X,F6.3,6X))')&
                  (i+j,wr(i+j)*(2*ry),f(i+j),j=0,nleft)
          ENDDO
       ENDIF
    ELSE
       k=0
       DO j=1,spin_mod%nsup
          DO i=j,spin_mod%nsup
             k=k+1
             ap(k)=gam(i,j)
          ENDDO
       ENDDO
       iopt=0
       CALL dspevy(iopt,ap,wr,dummy,spin_mod%nsup,spin_mod%nsup,fm1,3*spin_mod%nsup)
       k=0
       DO j=spin_mod%nsup+1,nstate
          DO i=j,nstate
             k=k+1
             ap(k)=gam(i,j)
          ENDDO
       ENDDO
       iopt=0
       CALL dspevy(iopt,ap,wr(spin_mod%nsup+1),dummy,spin_mod%nsdown,spin_mod%nsdown,fm1,3*spin_mod%nsdown)
       DO i=1,nstate
          ffi=f(i)
          IF (ffi.LT.1.e-5_real_8) ffi=1._real_8
          wr(i)=wr(i)/ffi
       ENDDO
       IF (ropt_mod%engpri) THEN
          IF (paral%io_parent)&
               WRITE(6,'(/," EIGENVALUES (EV) AND OCCUPATION:")')
          IF (paral%io_parent)&
               WRITE(6,'(" ALPHA SPIN ")')
          DO i=1,spin_mod%nsup,2
             nleft=MIN(spin_mod%nsup-i,1)
             IF (paral%io_parent)&
                  WRITE(6,'(3X,2(I3,F15.7,4X,F6.3,6X))')&
                  (i+j,wr(i+j)*(2*ry),f(i+j),j=0,nleft)
          ENDDO
          IF (paral%io_parent)&
               WRITE(6,'(" BETA SPIN ")')
          DO i=spin_mod%nsup+1,nstate,2
             nleft=MIN(nstate-i,1)
             IF (paral%io_parent)&
                  WRITE(6,'(3X,2(I3,F15.7,4X,F6.3,6X))')&
                  (i+j,wr(i+j)*(2*ry),f(i+j),j=0,nleft)
          ENDDO
       ENDIF
    ENDIF
    DEALLOCATE(fm1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(wr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(ap,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! CALL TIHALT('     REIGS',ISUB)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE reigs
  ! ==================================================================



END MODULE reigs_utils
