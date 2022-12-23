MODULE jacobi_c_utils
  USE kinds,                           ONLY: real_8
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: jacobi_c

CONTAINS

  ! ==================================================================
  SUBROUTINE jacobi_c(nm,n,a,w,eivr,scr,ierr)
    ! ==--------------------------------------------------------------==
    ! == VERSION FOR HERMITIAN MATRIX                                 ==
    ! == USE LAPACK ROUTINE                                           ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nm, n
    COMPLEX(real_8)                          :: a(nm,n)
    REAL(real_8)                             :: w(n)
    COMPLEX(real_8)                          :: eivr(nm,n)
    REAL(real_8)                             :: scr( (3*n-2) +2*(2*n) )
    INTEGER                                  :: ierr

    INTEGER                                  :: isub

! ==--------------------------------------------------------------==

    CALL tiset('  JACOBI_C',isub)
    CALL zheev('V','U',n,a(1,1),n,w(1),scr(3*n-2+1),2*n,scr(1),ierr)
    CALL dcopy(n*n*2,a(1,1),1,eivr(1,1),1)
    CALL tihalt('  JACOBI_C',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE jacobi_c
  ! ==================================================================


END MODULE jacobi_c_utils
