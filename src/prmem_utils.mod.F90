MODULE prmem_utils
  USE envj,                            ONLY: my_pid
  USE kinds,                           ONLY: real_8
  USE machine,                         ONLY: m_flush,&
                                             m_system
  USE parac,                           ONLY: paral

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: prmem

CONTAINS

  ! ==================================================================
  SUBROUTINE prmem(subr)
    ! ==--------------------------------------------------------------==
    CHARACTER(*)                             :: subr

! Variables

#if defined(__NO_MEMSIZE) || defined(__SR8000) || defined(__ES) || defined(__WINNT) || defined(__BG)
    IF (paral%io_parent.AND..FALSE.)&
         WRITE(6,'(A4,A10,A,10X,A)')&
         ' ***',SUBR,'| SIZE OF THE PROGRAM IS NOT AVAILABLE','***' 
#elif defined(__NEC)
    CHARACTER (len=140) :: string
    ! ==--------------------------------------------------------------==
#if defined(__ES)
    IF (paral%io_parent.AND..FALSE.)&
         WRITE(6,'(A,A,A)')&
         ' ***',SUBR,'| SIZE OF THE PROGRAM IS NOT AVAILABLE' 
#else
    IF (paral%io_parent) THEN
       WRITE(6,'(A,A,A,$)')&
            ' ***',subr,'| SIZE OF THE PROGRAM IS '
       CALL m_flush(6)
       WRITE(string,'(A,I10,A,A)') 'ps -lp ',my_pid,&
            ' | grep -v PID |  cut -c49-55 | ',&
            ' awk '' { print $1/25 " MBytes ***"}'''
       CALL m_system(string)
       WRITE(6,*)
    ENDIF
#endif
#elif defined(_vpp_)
    IF (paral%io_parent) THEN
       WRITE(6,'(A,A,A,$)')&
            ' ***',subr,'| SIZE OF THE PROGRAM IS '
       CALL m_flush(6)
       CALL printmemsize
       WRITE(6,'(A)')&
            ' kBYTES ***'
    ENDIF
#elif defined(__HP)
    REAL(real_8) :: rdata,vdata,rstack,vstack
    ! ==--------------------------------------------------------------==
    CALL memme(vdata,vstack)
    IF (paral%io_parent)&
         WRITE(6,'(A,A,A,F7.2,A,A,F7.2,A)')&
         ' ***',subr,'| PROGSIZE DATA=',vdata,' MByte',&
         '  /  STACK=',vstack,' MByte'
#elif defined(__SR8000)
    ! ==--------------------------------------------------------------==
    ! call to the shell is prohibitively expensive on sr8000
    IF (paral%io_parent.AND..FALSE.)&
         WRITE(6,'(A,A,A)') ' ***',subr,&
         '| SIZE OF THE PROGRAM IS NOT AVAILABLE ***'
    ! WRITE(6,'(A,A,A,$)')
    ! *    ' ***',SUBR,'| SIZE OF THE PROGRAM IS '
    ! WRITE(STRING,'(A,I10,A,A)') 'ps -o rss -p ',MY_PID,
    ! *  ' | sed -e "/RSS/d" | sed -e "s/K.*$/ KBytes ***/" ',
    ! *  ' | sed -e "s/M.*$/ MBytes ***/" '
    ! CALL SYSTEM(STRING)
#elif defined(__OSX) || defined(__OSX_IFC)
    CHARACTER (len=140) :: string
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent) THEN
       WRITE(6,'(A,A,A,$)')&
            ' ***',subr,'| SIZE OF THE PROGRAM IS '
       WRITE(string,'(A,I12,A)') 'ps -ovsz -p ',my_pid,&
            ' | sed -e "/VSZ/d" | sed -e "s/$/ KBytes ***/" '
       CALL m_system(string)
    ENDIF
#elif defined(__IBM)  && ! defined(__PWRLinux)
    CHARACTER (len=140) :: string
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent) THEN
       WRITE(6,'(A,A,A,$)')&
            ' ***',subr,'| SIZE OF THE PROGRAM IS '
       WRITE(string,'(A,I12,A)') 'ps -Fvsz -p ',my_pid,&
            ' | sed -e "/VSZ/d" | sed -e "s/$/ KBytes ***/" '
       CALL m_system(string)
    ENDIF
#else
    ! unknown platforms
    IF (paral%io_parent.AND..FALSE.)&
         WRITE(6,'(A4,A10,A,10X,A)')&
         ' ***',SUBR,'| SIZE OF THE PROGRAM IS NOT AVAILABLE','***' 
#endif
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE prmem
  ! ==================================================================
END MODULE prmem_utils
