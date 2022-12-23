MODULE summat_utils
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: int_1,&
                                             int_2,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8
  USE mp_interface,                    ONLY: mp_sum
  USE nlps,                            ONLY: imagp
  USE parac,                           ONLY: parai,&
                                             paral
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: summat
  PUBLIC :: summat_parent
  !!public :: sumhmat
  PUBLIC :: give_scr_summat

CONTAINS

  ! ==================================================================
  SUBROUTINE summat(a,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    REAL(real_8), TARGET                     :: a(nstate,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'summat'

    INTEGER                                  :: i, ierr, isub, j, k, n2, ntr
    REAL(real_8), ALLOCATABLE                :: aux(:,:)

! Variables
! ==--------------------------------------------------------------==
! == GLOBAL SUMMATION OF A SYMMETRIC MATRIX                       ==
! ==--------------------------------------------------------------==

    IF (parai%nproc.LE.1) RETURN
    CALL tiset('    SUMMAT',isub)
    n2 = nstate*nstate
    ntr = (nstate*(nstate+1))/2
    ALLOCATE(aux(ntr,1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    !$omp parallel do private(I,J,K) schedule(static,1)
#ifdef __SR8000
    !poption parallel, tlocal(I,J,K), cyclic 
#endif
    DO i=1,nstate
       DO j=i,nstate
          k=i+(j*(j-1))/2
          aux(k,1)=a(i,j)
       ENDDO
    ENDDO
    CALL mp_sum(aux,a,ntr,parai%allgrp)
    CALL dcopy(ntr,a(1,1),1,aux,1)
    !$omp parallel do private(I,J,K) schedule(static,1)
#ifdef __SR8000
    !poption parallel, cyclic 
#endif 
    DO i=1,nstate
       DO j=i,nstate
          k=i+(j*(j-1))/2
          a(i,j)=aux(k,1)
          a(j,i)=aux(k,1)
       ENDDO
    ENDDO
    DEALLOCATE(aux,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    CALL tihalt('    SUMMAT',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE summat
  ! ==================================================================
  SUBROUTINE summat_parent(a,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    REAL(real_8)                             :: a(nstate,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'summat_parent'

    INTEGER                                  :: i, ierr, isub, j, k, n2, ntr
    REAL(real_8), ALLOCATABLE                :: aux(:,:)

! Variables
! ==--------------------------------------------------------------==
! == GLOBAL SUMMATION OF A SYMMETRIC MATRIX                       ==
! ==--------------------------------------------------------------==

    IF (parai%nproc.LE.1) RETURN
    CALL tiset('    SUMMAT',isub)
    n2 = nstate*nstate
    ntr = (nstate*(nstate+1))/2
    ALLOCATE(aux(ntr,1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    !$omp parallel do private(I,J,K) schedule(static,1)
#ifdef __SR8000
    !poption parallel, tlocal(I,J,K), cyclic 
#endif 
    DO i=1,nstate
       DO j=i,nstate
          k=i+(j*(j-1))/2
          aux(k,1)=a(i,j)
       ENDDO
    ENDDO
    CALL mp_sum(aux,a,ntr,parai%source,parai%allgrp)
    IF (paral%parent) THEN
       CALL dcopy(ntr,a,1,aux,1)
       !$omp parallel do private(I,J,K) schedule(static,1)
#ifdef __SR8000
       !poption parallel, cyclic 
#endif 
       DO i=1,nstate
          DO j=i,nstate
             k=i+(j*(j-1))/2
             a(i,j)=aux(k,1)
             a(j,i)=aux(k,1)
          ENDDO
       ENDDO
    ENDIF
    DEALLOCATE(aux,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    CALL tihalt('    SUMMAT',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE summat_parent
  ! ==================================================================
  SUBROUTINE give_scr_summat(lsummat,tag,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lsummat
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate

! ==--------------------------------------------------------------==

    IF (parai%nproc.LE.1) THEN
       lsummat=0
    ELSE
       lsummat=imagp*(nstate*(nstate+1))/2
       tag   ='IMAGP*(NSTATE*(NSTATE+1))/2'
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_summat
  ! ==================================================================

END MODULE summat_utils

SUBROUTINE sumhmat(a,nstate)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE mp_interface, ONLY: mp_sum
  USE parac, ONLY : paral,parai
  IMPLICIT NONE
  INTEGER                                    :: nstate
  COMPLEX(real_8)                            :: a(nstate,nstate)

  CHARACTER(*), PARAMETER                    :: procedureN = 'sumhmat'

  COMPLEX(real_8), ALLOCATABLE               :: aux(:,:)
  INTEGER                                    :: i, ierr, isub, j, k, n2, ntr, &
                                                ntr2

! Variables
! ==--------------------------------------------------------------==
! == GLOBAL SUMMATION OF A HERMITIAN MATRIX                       ==
! ==--------------------------------------------------------------==

  IF (parai%nproc.LE.1) RETURN
  CALL tiset('   SUMHMAT',isub)
  n2 = nstate*nstate
  ntr = (nstate*(nstate+1))/2
  ntr2 = ntr*2
  ALLOCATE(aux(ntr2,1),STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
       __LINE__,__FILE__)
  k=0
  DO i=1,nstate
     DO j=i,nstate
        k=k+1
        aux(k,1)=a(i,j)
     ENDDO
  ENDDO
  CALL mp_sum(aux,a,ntr,parai%allgrp)
  CALL dcopy(ntr2,a(1,1),1,aux,1)
  k=0
  DO i=1,nstate
     DO j=i,nstate
        k=k+1
        a(i,j)=aux(k,1)
        IF (i.NE.j) a(j,i)=CMPLX(REAL(aux(k,1)),-AIMAG(aux(k,1)),kind=real_8)
     ENDDO
  ENDDO
  DEALLOCATE(aux,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
       __LINE__,__FILE__)
  CALL tihalt('   SUMHMAT',isub)
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE sumhmat
! ==================================================================
