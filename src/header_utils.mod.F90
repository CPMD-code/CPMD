MODULE header_utils
  USE envj,                            ONLY: curdir,&
                                             hname,&
                                             my_pid,&
                                             real_8,&
                                             tjlimit,&
                                             tmpdir,&
                                             user
  USE parac,                           ONLY: paral
  USE readsr_utils,                    ONLY: xstring

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: header

CONTAINS

  ! ==================================================================
  SUBROUTINE header(filename)
    ! ==--------------------------------------------------------------==
    ! ==  Writes a header to the standard output                      ==
    ! ==--------------------------------------------------------------==
    CHARACTER(len=*)                         :: filename

    CHARACTER(len=9)                         :: fformat
    INTEGER                                  :: i1, i2

! ==--------------------------------------------------------------==

    IF (paral%io_parent) THEN
       WRITE(6,'(/)')
       WRITE(6,'(12X,A)') '   ******  ******    ****  ****  ******   '
       WRITE(6,'(12X,A)') '  *******  *******   **********  *******  '
       WRITE(6,'(12X,A)') ' ***       **   ***  ** **** **  **   *** '
       WRITE(6,'(12X,A)') ' **        **   ***  **  **  **  **    ** '
       WRITE(6,'(12X,A)') ' **        *******   **      **  **    ** '
       WRITE(6,'(12X,A)') ' ***       ******    **      **  **   *** '
       WRITE(6,'(12X,A)') '  *******  **        **      **  *******  '
       WRITE(6,'(12X,A)') '   ******  **        **      **  ******   '
       WRITE(6,'(/,23X,A,A,/)') '   VERSION 4.3-', SVN_REV
#if defined(__GROMOS)
       WRITE(6,'(/,16X,A,/)') 'COMPILED WITH GROMOS-AMBER QM/MM SUPPORT'
#endif
       WRITE(6,'(14X,A)') '              COPYRIGHT'
       WRITE(6,'(14X,A)') '        IBM RESEARCH DIVISION'
       WRITE(6,'(14X,A)') '  MPI FESTKOERPERFORSCHUNG STUTTGART'
       WRITE(6,'(/,14X,A)')  '         The CPMD consortium'
       WRITE(6,'(14X,A)') '    Home Page: http://www.cpmd.org'
       WRITE(6,'(14X,A)') ' Mailing List: cpmd-list@cpmd.org'
       WRITE(6,'(14X,A)') '       E-mail: cpmd@cpmd.org'
       WRITE(6,'(/)')
       CALL timetag
       ! ==--------------------------------------------------------------==
       CALL xstring(filename,i1,i2)
       WRITE(fformat,'(A,I2,A)') '(A,T',MAX(21,65-(i2-i1)),',A)'
       WRITE(6,fformat) ' THE INPUT FILE IS: ',ADJUSTL(TRIM(filename))
       CALL xstring(hname,i1,i2)
       WRITE(fformat,'(A,I2,A)') '(A,T',MAX(20,65-(i2-i1)),',A)'
       WRITE(6,fformat) ' THIS JOB RUNS ON: ',ADJUSTL(TRIM(hname))
       WRITE(6,'(A)')   ' THE CURRENT DIRECTORY IS: '
       CALL xstring(curdir,i1,i2)
       WRITE(fformat,'(A,I2,A)') '(T',MAX(2,65-(i2-i1)),',A)'
       WRITE(6,fformat) ADJUSTL(TRIM(curdir))
       WRITE(6,'(A)')   ' THE TEMPORARY DIRECTORY IS: '
       CALL xstring(tmpdir,i1,i2)
       WRITE(fformat,'(A,I2,A)') '(T',MAX(2,65-(i2-i1)),',A)'
       WRITE(6,fformat) ADJUSTL(TRIM(tmpdir))
       WRITE(6,'(A,T50,I16)')  ' THE PROCESS ID IS: ',my_piD
       CALL xstring(user,i1,i2)
       IF (i2.NE.0) THEN
          WRITE(fformat,'(A,I2,A)') '(A,T',MAX(28,65-(i2-i1)),',A)'
          WRITE(6,fformat) ' THE JOB WAS SUBMITTED BY: ',&
               ADJUSTL(TRIM(user))
       ENDIF
       IF (tjlimit.NE.0._real_8) THEN
          WRITE(6,'(A,T48,F10.0,A8)')&
               ' THE JOB TIME LIMIT IS:',tjlimit,' SECONDS'
       ENDIF
       WRITE(6,*)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE header
  ! ==================================================================

END MODULE header_utils
