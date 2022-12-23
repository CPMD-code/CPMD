MODULE lowdin_utils
  USE csmat_utils,                     ONLY: csmat
  USE error_handling,                  ONLY: stopgm
  USE geq0mod,                         ONLY: geq0
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: paral
  USE reshaper,                        ONLY: reshape_inplace
  USE sfac,                            ONLY: fnl
  USE spin,                            ONLY: spin_mod
  USE system,                          ONLY: cntl,&
                                             ncpw
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: dspevy
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: lowdin
  PUBLIC :: give_scr_lowdin

CONTAINS

  ! ==================================================================
  SUBROUTINE lowdin(cp,cn,nstate)
    ! ==--------------------------------------------------------------==
    ! ==  LOWDIN ORTHOGONALIZATION                                    ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: cp(ncpw%ngw,*), cn(ncpw%ngw,*)
    INTEGER                                  :: nstate

    CHARACTER(*), PARAMETER                  :: procedureN = 'lowdin'

    INTEGER                                  :: i, ierr, ij, isub, j, k, kk
    REAL(real_8)                             :: fac
    REAL(real_8), ALLOCATABLE                :: a(:), aux(:), w(:)
    REAL(real_8), ALLOCATABLE, TARGET        :: z(:,:)
    REAL(real_8), POINTER                    :: z2(:)

    CALL tiset('    LOWDIN',isub)
    ! SCR partition
    ALLOCATE(a(nstate*nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(w(nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(z(nstate, nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(aux(3*nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL reshape_inplace(z, (/nstate*nstate/), z2)
    ! ==--------------------------------------------------------------==
    CALL csmat(z,cp,fnl,nstate,1)
    IF (.NOT.cntl%tlsd) THEN
       k=0
       DO j=1,nstate
          DO i=j,nstate
             k=k+1
             a(k)=z(i,j)
          ENDDO
       ENDDO
       CALL dspevy(1,a,w,z,nstate,nstate,aux,3*nstate)
       DO i=1,nstate
          w(i)=1._real_8/SQRT(w(i))
       ENDDO
       CALL zeroing(a)!,n2)
       DO k=1,nstate
          ij=1
          DO j=1,nstate
             fac=w(k)*z(j,k)
             CALL daxpy(nstate,fac,z(1,k),1,a(ij),1)
             ij=ij+nstate
          ENDDO
       ENDDO
       IF (ncpw%ngw.GT.0) THEN
          CALL dgemm('N','N',2*ncpw%ngw,nstate,nstate,1.0_real_8,cp(1,1),2*ncpw%ngw,&
               a(1),nstate,0.0_real_8,cn(1,1),2*ncpw%ngw)
          CALL dcopy(2*nstate*ncpw%ngw,cn(1,1),1,cp(1,1),1)
          IF (geq0) THEN
             DO i=1,nstate
                cp(1,i)=CMPLX(REAL(cp(1,i)),0.0_real_8,kind=real_8)
             ENDDO
          ENDIF
       ENDIF
       IF (paral%io_parent.AND.cntl%tlowdmat)THEN
          WRITE(6,'(1X,"Writing out Lowdin transformation matrix")')
          OPEN(41,file="LOWDIN_A",status='UNKNOWN')
          REWIND(unit=41)
          DO i=1,nstate
             WRITE(41,*)(a(i*nstate+j),j=1,nstate)
          ENDDO
          CLOSE(41)
       ENDIF
    ELSE
       ! ALPHA SPIN
       IF (spin_mod%nsup.GT.0) THEN
          k=0
          DO j=1,spin_mod%nsup
             DO i=j,spin_mod%nsup
                ij=(j-1)*nstate+i
                k=k+1
                a(k)=z2(ij)
             ENDDO
          ENDDO
          CALL dspevy(1,a,w,z2,spin_mod%nsup,spin_mod%nsup,aux,3*spin_mod%nsup)
          DO i=1,spin_mod%nsup
             w(i)=1._real_8/SQRT(w(i))
          ENDDO
          CALL zeroing(a(1:spin_mod%nsup**2))!,spin_mod%nsup*spin_mod%nsup)
          DO k=1,spin_mod%nsup
             ij=1
             DO j=1,spin_mod%nsup
                kk=(k-1)*spin_mod%nsup
                fac=w(k)*z2(kk+j)
                CALL daxpy(spin_mod%nsup,fac,z2(kk+1),1,a(ij),1)
                ij=ij+spin_mod%nsup
             ENDDO
          ENDDO
          IF (ncpw%ngw.GT.0) THEN
             CALL dgemm('N','N',2*ncpw%ngw,spin_mod%nsup,spin_mod%nsup,1.0_real_8,cp(1,1),2*ncpw%ngw,a(1),&
                  spin_mod%nsup,0.0_real_8,cn(1,1),2*ncpw%ngw)
             CALL dcopy(2*spin_mod%nsup*ncpw%ngw,cn(1,1),1,cp(1,1),1)
          ENDIF
          IF (paral%io_parent.AND.cntl%tlowdmat)THEN
             WRITE(6,'(1X,"Writing out Lowdin transformation matrix")')
             OPEN(41,file="LOWDIN_A",status='UNKNOWN')
             REWIND(unit=41)
             DO i=0,spin_mod%nsup-1
                WRITE(41,*)(a(i*spin_mod%nsup+j),j=1,spin_mod%nsup),(0.0,j=1,spin_mod%nsdown)
             ENDDO
          ENDIF
       ENDIF
       ! BETA SPIN
       IF (spin_mod%nsdown.GT.0) THEN
          k=0
          DO j=1,spin_mod%nsdown
             DO i=j,spin_mod%nsdown
                ij=spin_mod%nsup*nstate+(j-1)*nstate+spin_mod%nsup+i
                k=k+1
                a(k)=z2(ij)
             ENDDO
          ENDDO
          CALL dspevy(1,a,w,z2,spin_mod%nsdown,spin_mod%nsdown,aux,3*spin_mod%nsdown)
          DO i=1,spin_mod%nsdown
             w(i)=1._real_8/SQRT(w(i))
          ENDDO
          CALL zeroing(a(1:spin_mod%nsdown**2))!,spin_mod%nsdown*spin_mod%nsdown)
          DO k=1,spin_mod%nsdown
             ij=1
             DO j=1,spin_mod%nsdown
                kk=(k-1)*spin_mod%nsdown
                fac=w(k)*z2(kk+j)
                CALL daxpy(spin_mod%nsdown,fac,z2(kk+1),1,a(ij),1)
                ij=ij+spin_mod%nsdown
             ENDDO
          ENDDO
          IF (ncpw%ngw.GT.0) THEN
             CALL dgemm('N','N',2*ncpw%ngw,spin_mod%nsdown,spin_mod%nsdown,1.0_real_8,cp(1,spin_mod%nsup+1),&
                  2*ncpw%ngw,a(1),spin_mod%nsdown,0.0_real_8,cn(1,1),2*ncpw%ngw)
             CALL dcopy(2*spin_mod%nsdown*ncpw%ngw,cn(1,1),1,cp(1,spin_mod%nsup+1),1)
          ENDIF
          IF (paral%io_parent.AND.cntl%tlowdmat)THEN
             DO i=0,spin_mod%nsdown-1
                WRITE(41,*)(0.0,j=1,spin_mod%nsup),(a(i*spin_mod%nsdown+j),j=1,spin_mod%nsdown)
             ENDDO
             CLOSE(41)
          ENDIF
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    DEALLOCATE(a,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(w,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(z,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(aux,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL tihalt('    LOWDIN',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE lowdin
  ! ==================================================================
  SUBROUTINE give_scr_lowdin(llowdin,tag,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: llowdin
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate

! ==--------------------------------------------------------------==

    llowdin=2*nstate*nstate+4*nstate+10
    tag='2*NSTATE*NSTATE+4*NSTATE'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_lowdin
  ! ==================================================================

END MODULE lowdin_utils
