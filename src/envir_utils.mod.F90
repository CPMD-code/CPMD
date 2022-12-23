MODULE envir_utils
  USE envj,                            ONLY: curdir,&
                                             hname,&
                                             my_pid,&
                                             my_uid,&
                                             real_8,&
                                             tjlimit,&
                                             tmpdir,&
                                             user
  USE machine,                         ONLY: m_getcwd,&
                                             m_getenv,&
                                             m_getlog,&
                                             m_getpid,&
                                             m_getuid,&
                                             m_hostname

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: envir

CONTAINS

  ! ==================================================================
  SUBROUTINE envir
    ! ==--------------------------------------------------------------==
    ! ==  LOOK FOR THE ENVIRONMENT THIS JOB IS RUNNING                ==
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    ! ==  LOOK FOR THE ENVIRONMENT THIS JOB IS RUNNING                ==
    ! ==--------------------------------------------------------------==
    CALL m_getpid(my_pid)
    CALL m_getlog(user)
    CALL m_getuid(my_uid)
    CALL m_hostname(hname)
    CALL m_getcwd(curdir)
    tjlimit=0._real_8
    CALL m_getenv('TMPDIR',tmpdir)
    IF (tmpdir(1:1).EQ.' ') tmpdir=curdir
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE envir
  ! ==================================================================

END MODULE envir_utils
