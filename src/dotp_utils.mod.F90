MODULE dotp_utils

  USE error_handling,                  ONLY: stopgm
  USE geq0mod,                         ONLY: geq0
  USE kinds,                           ONLY: int_1,&
                                             int_2,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: dotp

  INTERFACE dotp
     MODULE PROCEDURE dotp_c
  END INTERFACE dotp

CONTAINS


  FUNCTION dotp_c(n,a,b) RESULT(dotp)
    INTEGER, INTENT(IN)                      :: n
    COMPLEX(real_8), INTENT(IN)              :: a(:), b(:)
    REAL(real_8)                             :: dotp

    CHARACTER(*), PARAMETER                  :: procedureN = 'dotp_c'

    REAL(real_8), EXTERNAL                   :: ddot

! ==--------------------------------------------------------------==

    dotp = 0._real_8
    IF (n.EQ.0) RETURN

    IF (n>SIZE(a).OR.n>SIZE(b)) CALL stopgm(procedureN,&
         'wrong dimension',&
         __LINE__,__FILE__)

    IF (geq0) THEN
       dotp = REAL(a(1),KIND=real_8) * REAL(b(1),KIND=real_8)
    ELSE
       dotp=2.0_real_8*( REAL(a(1),KIND=real_8)*REAL(b(1),KIND=real_8) + AIMAG(a(1))*AIMAG(b(1)))
    ENDIF
    IF (n.GT.1) THEN
       dotp=dotp+2.0_real_8*ddot(2*n-2,a(2),1,b(2),1)
    ENDIF
  END FUNCTION dotp_c

END MODULE dotp_utils

! ==================================================================
FUNCTION dotp(n,a,b)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE geq0mod , ONLY:geq0
  IMPLICIT NONE
  INTEGER                                    :: n
  REAL(real_8)                               :: a(2*n), b(2*n), dotp

#if defined (_vpp_) || defined (__SR8000) || defined (__ES)
  INTEGER                                    :: ii
  REAL(real_8)                               :: ddo
#else
  REAL(real_8), EXTERNAL                     :: ddot
#endif
  ! ==--------------------------------------------------------------==
  IF (n.EQ.0) THEN
     dotp=0._real_8
     RETURN
  ENDIF
  IF (geq0) THEN
     dotp=a(1)*b(1)
  ELSE
     dotp=2.0_real_8*(a(1)*b(1)+a(2)*b(2))
  ENDIF
  IF (n.GT.1) THEN
#if defined (_vpp_) || defined (__SR8000)  || defined (__ES)
     ddo=0._real_8
     !$omp parallel do private(II) reduction(+:DDO)
#ifdef __SR11000
     !poption parallel, tlocal(II), psum(DDO)
#endif
     DO ii=3,2*n
        ddo=ddo+a(ii)*b(ii)
     ENDDO
     dotp=dotp+ddo+ddo
#else
     dotp=dotp+2.0_real_8*ddot(2*n-2,a(3),1,b(3),1)
#endif
  ENDIF
  ! ==--------------------------------------------------------------==
  RETURN
END FUNCTION dotp
! ==================================================================
