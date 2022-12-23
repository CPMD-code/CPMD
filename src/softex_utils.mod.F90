MODULE softex_utils
  USE machine,                         ONLY: m_datum
  USE parac,                           ONLY: paral
  USE soft,                            ONLY: soft_com

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: softex

CONTAINS

  ! ==================================================================
  SUBROUTINE softex(signum)
    ! ==--------------------------------------------------------------==
    ! ==  THIS ROUTINE IS CALLED WHEN THE USER SPECIFIES SIGNAL 30    ==
    ! ==  BY THE COMMAND "kill -30 PID"                              ==
    ! ==--------------------------------------------------------------==
    INTEGER, INTENT(in)                      :: signum

    CHARACTER(len=26)                        :: datx

    CALL m_datum(datx)
    soft_com%exsoft=.TRUE.
    IF (paral%io_parent)&
         WRITE(6,'(1X,64("*"))')
    IF (paral%io_parent)&
         WRITE(6,'(1X,"*",62X,"*")')
    IF (paral%io_parent)&
         WRITE(6,'(1X,"*",13X,A,14X,"*")')&
         'CPMD RECEIVED THE SOFT EXIT REQUEST'
    IF (paral%io_parent)&
         WRITE(6,'(A,A,A,A)')&
         ' *     ',' THE COMMAND WAS ISSUED AT ',datx(1:24),'      *'
    IF (paral%io_parent)&
         WRITE(6,'(1X,"*",62X,"*")')
    IF (paral%io_parent)&
         WRITE(6,'(1X,64("*"))')
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE softex
  ! ==================================================================

END MODULE softex_utils
