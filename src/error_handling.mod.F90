MODULE error_handling
  USE parac,                           ONLY: parai

  !$ USE omp_lib, ONLY: omp_get_thread_num, omp_get_level

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: stopgm
CONTAINS
  ! ==================================================================
  SUBROUTINE stopgm(a,b,line,file)
    ! ==--------------------------------------------------------------==
    CHARACTER(len=*)                         :: a, b
    INTEGER                                  :: line
    CHARACTER(len=*)                         :: file

    CHARACTER(*), PARAMETER                  :: file_name_base = 'LocalError'
    INTEGER, PARAMETER                       :: file_unit = 666

    CHARACTER(100)                           :: buff, file_name
    EXTERNAL                                 :: tistopgm
    INTEGER                                  :: i_level, i_thread, nc

! ==--------------------------------------------------------------==

    !$omp master
    CALL end_swap
    !$omp end master

    i_thread = 0
    i_level = 0
    !$ i_thread = omp_get_thread_num( )
    !$ i_level = omp_get_level( )
    file_name=file_name_base
    WRITE(buff,'(i0,A,i0,A,i0)') parai%cp_me,'-',i_thread,'-',i_level
    file_name=TRIM(file_name)//'-'//TRIM(ADJUSTL(buff))//'.log'
    OPEN(unit=file_unit,file=file_name,action='write')
    WRITE(file_unit,'(5(A,I0))') ' process id''s: ',&
         parai%cp_me,', ',parai%me,', ',parai%cp_inter_me,', ',i_thread,', ',i_level
    WRITE(file_unit,'(A,A)')&
         ' process stops in file: ',TRIM(ADJUSTL(file))
    WRITE(file_unit,'(A,I0)')&
         '               at line: ',line
    WRITE(file_unit,'(A,A)')&
         '               in procedure: ',TRIM(ADJUSTL(a))
    WRITE(file_unit,'(A,A)') ' error message: ',TRIM(ADJUSTL(b))
    CALL tistopgm(file_unit)
    CLOSE(file_unit)

    nc=999
    CALL my_stopall(nc)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE stopgm
END MODULE error_handling
