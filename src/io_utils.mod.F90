#include "cpmd_global.h"

MODULE io_utils
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: int_8,&
       real_8
  USE string_utils,                    ONLY: int2str

#ifdef __PARALLEL
  USE mpi
#endif


  IMPLICIT NONE

  PRIVATE

  PUBLIC :: file_close
  PUBLIC :: file_open
  PUBLIC :: file_seek
  PUBLIC :: file_get_position
  PUBLIC :: file_write_l
  PUBLIC :: file_write_char
  PUBLIC :: file_write_i
  PUBLIC :: file_write_d
  PUBLIC :: file_write_z
  PUBLIC :: file_read_l
  PUBLIC :: file_read_char
  PUBLIC :: file_read_i
  PUBLIC :: file_read_d
  PUBLIC :: file_read_z

CONTAINS

#if defined(__HASNT_F03_FEATURES)
#define HASNT_F03_FEATURES
#endif

  SUBROUTINE file_close(unit,use_mpi)
    INTEGER                                  :: unit
    LOGICAL                                  :: use_mpi

    CHARACTER(*), PARAMETER                  :: procedureN = 'file_close'

    CHARACTER(256)                           :: file_name, iomsg
    INTEGER                                  :: ierr, iostat

    iostat=0
    iomsg=''
#if defined(__PARALLEL)
    IF (use_mpi) THEN
       CALL mpi_file_close(unit,ierr)
    ELSE
#endif
#if ! defined(HASNT_F03_FEATURES)
       INQUIRE(unit=unit,name=file_name,&
            iomsg=iomsg,&
            iostat=iostat)
#else
       INQUIRE(unit=unit,name=file_name,&
            iostat=iostat)
#endif
       IF (iostat /= 0 ) CALL stopgm(procedureN,&
            'Problem inquiring FILE = '//&
            ' IOMSG = '//TRIM(iomsg)//&
            ' IOSTAT = '//TRIM(int2str(iostat)),& 
            __LINE__,__FILE__)

#if ! defined(HASNT_F03_FEATURES)
       CLOSE(unit=unit,&
            iomsg=iomsg,&
            iostat=iostat)
#else
       CLOSE(unit=unit,&
            iostat=iostat)
#endif
       IF (iostat /= 0 ) CALL stopgm(procedureN,&
            'Problem closing FILE = '//TRIM(file_name)//&
            ' IOMSG = '//TRIM(iomsg)//&
            ' IOSTAT = '//TRIM(int2str(iostat)),& 
            __LINE__,__FILE__)

#if defined(__PARALLEL)
    ENDIF
#endif
    RETURN
  END SUBROUTINE file_close
  SUBROUTINE file_open(file_name,form,status,action,&
       is_stream,use_mpi,grp,unit)
    CHARACTER(*)                             :: file_name, form, status, &
                                                action
    LOGICAL                                  :: is_stream, use_mpi
    INTEGER                                  :: grp, unit

    CHARACTER(*), PARAMETER                  :: procedureN = 'file_open'

    CHARACTER(256)                           :: iomsg
    INTEGER                                  :: iostat

#if defined(__PARALLEL) && ! defined(_HASNT_BF_STREAM_IO)
    !    INCLUDE 'mpif.h'
    INTEGER :: ierr,mode
#endif
    iostat=0
    iomsg=''
    IF(.FALSE.) THEN
       WRITE(6,*) procedureN,' opening file: NAME=',TRIM(FILE_NAME), &
            ' FORM=',TRIM(FORM),' STATUS=',TRIM(STATUS), & 
            ' ACTION=',TRIM(ACTION),' IS_STREAM=',IS_STREAM, &
            ' USE_MPI=',USE_MPI,' GRP=',GRP
    ENDIF
#if ! defined(_HASNT_BF_STREAM_IO)
    IF (is_stream) THEN
#if defined(__PARALLEL)
       IF (use_mpi) THEN
          SELECT CASE(action)
          CASE('READ' ); mode=MPI_MODE_RDONLY
          CASE('WRITE'); mode=MPI_MODE_WRONLY
          CASE('READWRITE'); mode=MPI_MODE_RDWR
          CASE default
             CALL stopgm(procedureN,'invalid ACTION',& 
                  __LINE__,__FILE__)
          END SELECT
          CALL mpi_file_open(grp,file_name,mode,MPI_INFO_NULL,&
               unit,ierr)
       ELSE
#endif

#if ! defined(HASNT_F03_FEATURES)
          OPEN(unit=unit,file=file_name,form=form,status=status,&
               action=action,&
               access='stream',&
               iomsg=iomsg,&
               iostat=iostat)
#else
          OPEN(unit=unit,file=file_name,form=form,status=status,&
               action=action,&
               access='stream',&
               iostat=iostat)
#endif

#if defined(__PARALLEL)
       ENDIF
#endif
    ELSE
#endif
#if ! defined(HASNT_F03_FEATURES)
       OPEN(unit=unit,file=file_name,form=form,status=status,&
            action=action,&
            iomsg=iomsg,&
            iostat=iostat)
#else
       OPEN(unit=unit,file=file_name,form=form,status=status,&
            action=action,&
            iostat=iostat)
#endif

#if ! defined(_HASNT_BF_STREAM_IO)
    ENDIF
#endif
    IF (iostat /= 0 ) CALL stopgm(procedureN,&
         'Problem opening FILE = '//TRIM(file_name)//&
         ' IOMSG = '//TRIM(iomsg)//&
         ' IOSTAT = '//TRIM(int2str(iostat)),& 
         __LINE__,__FILE__)
    ! write(6,*) 'opening file: UNIT=',unit
    RETURN
  END SUBROUTINE file_open
  SUBROUTINE file_seek(unit,use_mpi,offset)
    INTEGER                                  :: unit
    LOGICAL                                  :: use_mpi
    INTEGER(int_8)                           :: offset

#if defined(__PARALLEL) && ! defined(_HASNT_BF_STREAM_IO)
    !    INCLUDE 'mpif.h'
    INTEGER :: ierr
    INTEGER(kind=MPI_OFFSET_KIND) :: offset_tmp
#endif
#if ! defined(_HASNT_BF_STREAM_IO)
    IF (use_mpi) THEN
#if defined(__PARALLEL)
       offset_tmp = offset
       CALL mpi_file_seek(unit,offset_tmp,MPI_SEEK_SET,ierr)
#else
       READ(unit,pos=offset)
#endif
    ELSE
       READ(unit,pos=offset)
    ENDIF
#else
    offset = -1
#endif
  END SUBROUTINE file_seek
  SUBROUTINE file_get_position(unit,use_mpi,offset)
    INTEGER                                  :: unit
    LOGICAL                                  :: use_mpi
    INTEGER(int_8)                           :: offset

#if defined(__PARALLEL) && ! defined(_HASNT_BF_STREAM_IO)
    !    INCLUDE 'mpif.h'
    INTEGER :: ierr
    INTEGER(kind=MPI_OFFSET_KIND) :: offset_tmp
#endif
#if ! defined(_HASNT_BF_STREAM_IO)
    IF (use_mpi) THEN
#if defined(__PARALLEL)
       CALL mpi_file_get_position(unit,offset_tmp,ierr)
       offset = offset_tmp
#else
       INQUIRE(unit,pos=offset)
#endif
    ELSE
       INQUIRE(unit,pos=offset)
    ENDIF
#else
    offset = -1
#endif
  END SUBROUTINE file_get_position
  SUBROUTINE file_write_l(unit,n,DATA,use_mpi)
    INTEGER                                  :: unit, n
    LOGICAL                                  :: DATA(*), use_mpi

#if defined(__PARALLEL) && ! defined(_HASNT_BF_STREAM_IO)
    !    INCLUDE 'mpif.h'
    INTEGER :: ierr
#endif
#if ! defined(_HASNT_BF_STREAM_IO)
    IF (use_mpi) THEN
#if defined(__PARALLEL)
       CALL mpi_file_write(unit,DATA,n,MPI_LOGICAL,&
            MPI_STATUS_IGNORE,ierr)
#else
       WRITE(unit) DATA(1:n)
#endif
    ELSE
       WRITE(unit) DATA(1:n)
    ENDIF
#else
    WRITE(unit) DATA(1:n)
#endif
  END SUBROUTINE file_write_l
  SUBROUTINE file_write_char(unit,n,DATA,use_mpi)
    INTEGER                                  :: unit, n
    CHARACTER(1)                             :: DATA(*)
    LOGICAL                                  :: use_mpi

#if defined(__PARALLEL) && ! defined(_HASNT_BF_STREAM_IO)
    !    INCLUDE 'mpif.h'
    INTEGER :: ierr
#endif
#if ! defined(_HASNT_BF_STREAM_IO)
    IF (use_mpi) THEN
#if defined(__PARALLEL)
       CALL mpi_file_write(unit,DATA,n,MPI_CHARACTER,&
            MPI_STATUS_IGNORE,ierr)
#else
       WRITE(unit) DATA(1:n)
#endif
    ELSE
       WRITE(unit) DATA(1:n)
    ENDIF
#else
    WRITE(unit) DATA(1:n)
#endif
  END SUBROUTINE file_write_char
  SUBROUTINE file_write_i(unit,n,DATA,use_mpi)
    INTEGER                                  :: unit, n, DATA(*)
    LOGICAL                                  :: use_mpi

#if defined(__PARALLEL) && ! defined(_HASNT_BF_STREAM_IO)
    !    INCLUDE 'mpif.h'
    INTEGER :: ierr
#endif
#if ! defined(_HASNT_BF_STREAM_IO)
    IF (use_mpi) THEN
#if defined(__PARALLEL)
       CALL mpi_file_write(unit,DATA,n,MPI_INTEGER,&
            MPI_STATUS_IGNORE,ierr)
#else
       WRITE(unit) DATA(1:n)
#endif
    ELSE
       WRITE(unit) DATA(1:n)
    ENDIF
#else
    WRITE(unit) DATA(1:n)
#endif
  END SUBROUTINE file_write_i
  SUBROUTINE file_write_d(unit,n,DATA,use_mpi)
    INTEGER                                  :: unit, n
    REAL(real_8)                             :: DATA(*)
    LOGICAL                                  :: use_mpi

#if defined(__PARALLEL) && ! defined(_HASNT_BF_STREAM_IO)
    !    INCLUDE 'mpif.h'
    INTEGER :: ierr
#endif
#if ! defined(_HASNT_BF_STREAM_IO)
    IF (use_mpi) THEN
#if defined(__PARALLEL)
       CALL mpi_file_write(unit,DATA,n,MPI_DOUBLE_PRECISION,&
            MPI_STATUS_IGNORE,ierr)
#else
       WRITE(unit) DATA(1:n)
#endif
    ELSE
       WRITE(unit) DATA(1:n)
    ENDIF
#else
    WRITE(unit) DATA(1:n)
#endif
  END SUBROUTINE file_write_d
  SUBROUTINE file_write_z(unit,n,DATA,use_mpi)
    INTEGER                                  :: unit, n
    COMPLEX(real_8)                          :: DATA(*)
    LOGICAL                                  :: use_mpi

#if defined(__PARALLEL) && ! defined(_HASNT_BF_STREAM_IO)
    !    INCLUDE 'mpif.h'
    INTEGER :: ierr
#endif
#if ! defined(_HASNT_BF_STREAM_IO)
    IF (use_mpi) THEN
#if defined(__PARALLEL)
       CALL mpi_file_write(unit,DATA,n,MPI_DOUBLE_COMPLEX,&
            MPI_STATUS_IGNORE,ierr)
#else
       WRITE(unit) DATA(1:n)
#endif
    ELSE
       WRITE(unit) DATA(1:n)
    ENDIF
#else
    WRITE(unit) DATA(1:n)
#endif
  END SUBROUTINE file_write_z
  SUBROUTINE file_read_l(unit,n,DATA,use_mpi)
    INTEGER                                  :: unit, n
    LOGICAL                                  :: DATA(*), use_mpi

#if defined(__PARALLEL) && ! defined(_HASNT_BF_STREAM_IO)
    !    INCLUDE 'mpif.h'
    INTEGER :: ierr
#endif
#if ! defined(_HASNT_BF_STREAM_IO)
    IF (use_mpi) THEN
#if defined(__PARALLEL)
       CALL mpi_file_read(unit,DATA,n,MPI_LOGICAL,&
            MPI_STATUS_IGNORE,ierr)
#else
       READ(unit) DATA(1:n)
#endif
    ELSE
       READ(unit) DATA(1:n)
    ENDIF
#else
    READ(unit) DATA(1:n)
#endif
  END SUBROUTINE file_read_l
  SUBROUTINE file_read_char(unit,n,DATA,use_mpi)
    INTEGER                                  :: unit, n
    CHARACTER(1)                             :: DATA(*)
    LOGICAL                                  :: use_mpi

#if defined(__PARALLEL) && ! defined(_HASNT_BF_STREAM_IO)
    !    INCLUDE 'mpif.h'
    INTEGER :: ierr
#endif
#if ! defined(_HASNT_BF_STREAM_IO)
    IF (use_mpi) THEN
#if defined(__PARALLEL)
       CALL mpi_file_read(unit,DATA,n,MPI_CHARACTER,&
            MPI_STATUS_IGNORE,ierr)
#else
       READ(unit) DATA(1:n)
#endif
    ELSE
       READ(unit) DATA(1:n)
    ENDIF
#else
    READ(unit) DATA(1:n)
#endif
  END SUBROUTINE file_read_char
  SUBROUTINE file_read_i(unit,n,DATA,use_mpi)
    INTEGER                                  :: unit, n, DATA(*)
    LOGICAL                                  :: use_mpi

#if defined(__PARALLEL) && ! defined(_HASNT_BF_STREAM_IO)
    !    INCLUDE 'mpif.h'
    INTEGER :: ierr
#endif
#if ! defined(_HASNT_BF_STREAM_IO)
    IF (use_mpi) THEN
#if defined(__PARALLEL)
       CALL mpi_file_read(unit,DATA,n,MPI_INTEGER,&
            MPI_STATUS_IGNORE,ierr)
#else
       READ(unit) DATA(1:n)
#endif
    ELSE
       READ(unit) DATA(1:n)
    ENDIF
#else
    READ(unit) DATA(1:n)
#endif
  END SUBROUTINE file_read_i
  SUBROUTINE file_read_d(unit,n,DATA,use_mpi)
    INTEGER                                  :: unit, n
    REAL(real_8)                             :: DATA(*)
    LOGICAL                                  :: use_mpi

#if defined(__PARALLEL) && ! defined(_HASNT_BF_STREAM_IO)
    !    INCLUDE 'mpif.h'
    INTEGER :: ierr
#endif
#if ! defined(_HASNT_BF_STREAM_IO)
    IF (use_mpi) THEN
#if defined(__PARALLEL)
       CALL mpi_file_read(unit,DATA,n,MPI_DOUBLE_PRECISION,&
            MPI_STATUS_IGNORE,ierr)
#else
       READ(unit) DATA(1:n)
#endif
    ELSE
       READ(unit) DATA(1:n)
    ENDIF
#else
    READ(unit) DATA(1:n)
#endif
  END SUBROUTINE file_read_d
  SUBROUTINE file_read_z(unit,n,DATA,use_mpi)
    INTEGER                                  :: unit, n
    COMPLEX(real_8)                          :: DATA(*)
    LOGICAL                                  :: use_mpi

#if defined(__PARALLEL) && ! defined(_HASNT_BF_STREAM_IO)
    !    INCLUDE 'mpif.h'
    INTEGER :: ierr
#endif
#if ! defined(_HASNT_BF_STREAM_IO)
    IF (use_mpi) THEN
#if defined(__PARALLEL)
       CALL mpi_file_read(unit,DATA,n,MPI_DOUBLE_COMPLEX,&
            MPI_STATUS_IGNORE,ierr)
#else
       READ(unit) DATA(1:n)
#endif
    ELSE
       READ(unit) DATA(1:n)
    ENDIF
#else
    READ(unit) DATA(1:n)
#endif
  END SUBROUTINE file_read_z

END MODULE io_utils
