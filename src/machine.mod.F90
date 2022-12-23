#include "cpmd_global.h"

MODULE machine
  USE kinds,                           ONLY: real_8

  USE, INTRINSIC :: ISO_C_BINDING,  ONLY: C_INT, C_NULL_CHAR, C_CHAR, C_PTR, C_NULL_PTR, C_ASSOCIATED, C_F_POINTER

  !$ use omp_lib, only : omp_get_wtime



  IMPLICIT NONE

  PRIVATE


  PUBLIC :: m_getenv
  PUBLIC :: m_datum
  PUBLIC :: m_walltime
  PUBLIC :: m_cputime
  PUBLIC :: m_flush
  PUBLIC :: m_getcwd
  PUBLIC :: m_getarg
  PUBLIC :: m_iargc
  PUBLIC :: m_system
  PUBLIC :: m_getpid
  PUBLIC :: m_getuid
  PUBLIC :: m_hostname
  PUBLIC :: m_sleep
  PUBLIC :: m_signal
  PUBLIC :: m_getlog
  PUBLIC :: m_compiler_version
  PUBLIC :: m_compiler_options
  PUBLIC :: m_getpagesize


  INTERFACE

     FUNCTION getcwd(buf, buflen) BIND(C,name="getcwd") RESULT(reslt)
       IMPORT :: C_CHAR, C_INT, C_PTR
       CHARACTER(C_CHAR), DIMENSION(*) :: buf
       INTEGER(C_INT), VALUE           :: buflen
       TYPE(C_PTR)                     :: reslt
     END FUNCTION getcwd

     FUNCTION getpid() BIND(C,name="getpid") RESULT(pid)
       IMPORT :: C_INT
       INTEGER(C_INT)                  :: pid
     END FUNCTION getpid

     FUNCTION getuid() BIND(C,name="getuid") RESULT(uid)
       IMPORT :: C_INT
       INTEGER(C_INT)                  :: uid
     END FUNCTION getuid

     FUNCTION  gethostname(buf, buflen) BIND(C,name="gethostname") RESULT(errno)
       IMPORT :: C_CHAR, C_INT
       CHARACTER(C_CHAR), DIMENSION(*)     :: buf
       INTEGER(C_INT), VALUE               :: buflen
       INTEGER(C_INT)                      :: errno
     END FUNCTION gethostname

     FUNCTION sleep(seconds) BIND(C,NAME="sleep") RESULT(errno)
       IMPORT :: C_INT
       INTEGER(C_INT), VALUE               :: seconds
       INTEGER(C_INT)                      :: errno
     END FUNCTION sleep

     INTEGER(C_INT) FUNCTION getpagesize() BIND(C,NAME="getpagesize")
       IMPORT :: C_INT
       IMPLICIT NONE
     END FUNCTION getpagesize

  END INTERFACE

  ABSTRACT INTERFACE
     SUBROUTINE signal_interface(signum)
       INTEGER, INTENT(in)                      :: signum

     END SUBROUTINE signal_interface
  END INTERFACE


CONTAINS


  SUBROUTINE m_compiler_version ( compiler )

#if defined(_HASNT_F08_ISO_FORTRAN_ENV)

    CHARACTER( * ), INTENT( OUT ) :: compiler

    compiler = 'unknown'

#else

    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : compiler_version

    CHARACTER( * ), INTENT( OUT ) :: compiler

    compiler = COMPILER_VERSION ( )

#endif

  END SUBROUTINE m_compiler_version


  SUBROUTINE m_compiler_options ( options )

#if defined(_HASNT_F08_ISO_FORTRAN_ENV)

    CHARACTER( * ), INTENT( OUT ) :: options

    options = 'unknown'

#else

    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : compiler_options

    CHARACTER( * ), INTENT( OUT ) :: options

    options = COMPILER_OPTIONS ( )

#endif

  END SUBROUTINE m_compiler_options


  SUBROUTINE m_getenv(name,val)
    CHARACTER(*), INTENT(in)                 :: name
    CHARACTER(*), INTENT(out)                :: val

    val = ''
    CALL GET_ENVIRONMENT_VARIABLE(name,val)

  END SUBROUTINE m_getenv


  SUBROUTINE m_datum(date_time)
    CHARACTER(*), INTENT(out)                :: date_time

    CHARACTER(10)                            :: time
    CHARACTER(8)                             :: date

    CALL DATE_AND_TIME(date=date,time=time)
    date_time=date(1:4)//"-"//date(5:6)//"-"//date(7:8)//" "//time(1:2)//":"//time(3:4)//":"//time(5:10)

  END SUBROUTINE m_datum


  ! For the time being this subroutine MUST return the time in milliseconds
  REAL(real_8) FUNCTION m_walltime()
    INTEGER :: count
    LOGICAL, SAVE :: first_time = .TRUE.
    INTEGER, SAVE :: count_max, count_rate, ncycles, previous_count

    !$ if(.false.) then
    IF ( first_time ) THEN
       CALL SYSTEM_CLOCK( count_rate=count_rate, count_max=count_max )
       ncycles = 0
       previous_count = 0
       first_time = .FALSE.
    ENDIF

    CALL SYSTEM_CLOCK( count=count )

    IF ( count < previous_count) ncycles = ncycles + 1

    previous_count = count

    m_walltime = ( REAL( count ,kind=real_8) + REAL( ncycles ,kind=real_8) * REAL( count_max ,kind=real_8) )&
         / REAL( count_rate ,kind=real_8) * 1000._real_8

    !$ else
    !$   m_walltime = omp_get_wtime() * 1000._real_8
    !$ endif

  END FUNCTION m_walltime

  ! This function returns the time in SECONDS: if there is any problem, it is
  ! entirely to be blamed the compiler! no hacks here!
  FUNCTION m_cputime()
    REAL(real_8)                             :: m_cputime

    CALL CPU_TIME(m_cputime)

  END FUNCTION m_cputime


  SUBROUTINE m_flush(iunit)
    INTEGER                                  :: iunit

    FLUSH(iunit)

  END SUBROUTINE m_flush


  SUBROUTINE m_getcwd(curdir)
    CHARACTER(len=*), INTENT(out)            :: curdir

    CHARACTER(LEN(curdir)+1,C_CHAR), TARGET  :: tmp
    TYPE(C_PTR)                              :: c_stat
    INTEGER(C_INT) :: l

    l = LEN(tmp,KIND=C_INT)!vw bypass pgi bug
    c_stat = getcwd(tmp, l)
    IF(.NOT. C_ASSOCIATED(c_stat)) STOP 12345
    curdir = ''
    curdir = tmp(1:INDEX(tmp,C_NULL_CHAR)-1)

  END SUBROUTINE m_getcwd


  SUBROUTINE m_getarg(m,string)
    INTEGER, INTENT(in)                      :: m
    CHARACTER(len=*), INTENT(out)            :: string

    INTEGER                                  :: istat

    CALL GET_COMMAND_ARGUMENT(m, string, STATUS=istat)
    IF( istat /= 0 ) STOP 12345
    string = TRIM(string)

  END SUBROUTINE m_getarg


  FUNCTION m_iargc() RESULT (ic)
    INTEGER                                  :: ic

    ic = COMMAND_ARGUMENT_COUNT()

  END FUNCTION m_iargc


  SUBROUTINE m_system(command)
    CHARACTER(len=*)                         :: command

#if defined(_HASNT_F03_EXECUTE_COMMAND_LINE)

    CALL system(command)

#else

    CALL EXECUTE_COMMAND_LINE(command)

#endif

  END SUBROUTINE m_system


  SUBROUTINE m_getpid(pid)
    INTEGER, INTENT(OUT)                     :: pid

    pid = getpid()

  END SUBROUTINE m_getpid


  SUBROUTINE m_getuid(uid)
    INTEGER, INTENT(OUT)                     :: uid

    uid = getuid( )

  END SUBROUTINE m_getuid


  SUBROUTINE m_hostname(hname)
    CHARACTER(len=*), INTENT(OUT)            :: hname

    CHARACTER(LEN(hname)+1,C_CHAR)           :: tmp
    INTEGER                                  :: istat
    INTEGER(C_INT) :: l

    l = LEN(tmp,KIND=C_INT)!vw bypass pgi bug
    istat = gethostname(tmp, l)
    IF(istat /= 0) STOP 12346

    hname = ''
    hname = tmp(1:INDEX(tmp,C_NULL_CHAR)-1)

  END SUBROUTINE m_hostname


  SUBROUTINE m_sleep(seconds)
    INTEGER, INTENT(IN)                      :: seconds

    INTEGER                                  :: istat

    istat = sleep(INT(seconds,C_INT))

  END SUBROUTINE m_sleep


  SUBROUTINE m_signal(number,handler)
    INTEGER, INTENT(in)                      :: number
    PROCEDURE(signal_interface)              :: handler

  END SUBROUTINE m_signal


  SUBROUTINE m_getlog(string)
    CHARACTER(len=*)                         :: string

    string='xxx'

  END SUBROUTINE m_getlog


  FUNCTION m_getpagesize() RESULT(reslt)
    INTEGER :: reslt

    reslt = getpagesize()

  END FUNCTION m_getpagesize


END MODULE machine
