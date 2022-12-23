MODULE envj
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==================================================================
  ! == ENVIRONMENT                                                  ==
  ! ==================================================================
  ! == USER   User name                                             ==
  ! == HNAME  Computer name                                         ==
  ! == CURDIR Current Directory                                     ==
  ! == TMPDIR Temporary Directory                                   ==
  ! ==--------------------------------------------------------------==
  ! == MY_PID Process ID Number                                     ==
  ! == MY_UID User ID Number                                        ==
  ! ==--------------------------------------------------------------==
  ! == TJLIMIT Job or CPU limit time (seconds)                      ==
  ! ==--------------------------------------------------------------==
  CHARACTER (len=10) :: user=''
  CHARACTER (len=64) :: hname=''
  CHARACTER (len=255) :: curdir='', tmpdir=''
  INTEGER :: my_pid,my_uid
  REAL(real_8) :: tjlimit
  ! ==================================================================

END MODULE envj
