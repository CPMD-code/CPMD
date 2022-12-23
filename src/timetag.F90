! ==================================================================
SUBROUTINE timetag
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE parac, ONLY : paral
  IMPLICIT NONE

#if defined(_vpp_)
#elif defined(__NEC)
#elif defined(__HP)
#else
  IF (paral%io_parent)&
       WRITE(6,'(18X,5A,/)') '***  ',&
       __DATE__ , ' -- ',&
       __TIME__&
       ,'  ***'
#endif
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE timetag
! ==================================================================
