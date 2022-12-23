MODULE jacobi_utils
  USE kinds,                           ONLY: real_8
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: jacobi

CONTAINS

  ! ==================================================================
  SUBROUTINE jacobi(nm,n,a,w,eivr,ierr)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nm, n
    REAL(real_8)                             :: a(nm,n), w(n), eivr(nm,n)
    INTEGER                                  :: ierr

    INTEGER                                  :: info, isub
    REAL(real_8)                             :: work(10)

    CALL tiset('    JACOBI',isub)
    IF (nm.LE.2) THEN
       CALL dsyev('V','U',n,a,nm,w,work,10,info)
    ELSE
       CALL dsyev('V','U',n,a,nm,w,eivr,nm*n,info)
    ENDIF
    CALL dcopy(nm*n,a(1,1),1,eivr(1,1),1)
    CALL tihalt('    JACOBI',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE jacobi
  ! ==================================================================

END MODULE jacobi_utils
