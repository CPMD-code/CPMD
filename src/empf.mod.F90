MODULE empf
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==================================================================
  ! == USED BY FSTART TO CALCULATE EMPIRICAL FORCE CONSTANT MATRIX  ==
  ! ==================================================================
  ! ==--------------------------------------------------------------==
  REAL(real_8), ALLOCATABLE :: zan(:)
  REAL(real_8), ALLOCATABLE :: c(:,:)

  ! ==--------------------------------------------------------------==
  INTEGER, ALLOCATABLE :: ibind(:)

  INTEGER :: naty
  ! ==--------------------------------------------------------------==
  ! == PARAMETERS USED BY DISCO ROUTINES                            ==
  ! ==--------------------------------------------------------------==
  INTEGER, PARAMETER :: ir=5,iw=6 
  ! ==--------------------------------------------------------------==
  REAL(real_8), PARAMETER  :: zero=0._real_8,one=1._real_8,two=2._real_8,&
       three=3._real_8,four=4._real_8,five=5._real_8,big=1.e30_real_8,&
       dp12=0.5_real_8,dp34=0.75_real_8,dp32=1.5_real_8
  ! ==================================================================
END MODULE empf
