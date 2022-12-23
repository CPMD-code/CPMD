MODULE velocitinp_utils
  USE coor,                            ONLY: lvelini,&
                                             velp
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions0
  USE kinds,                           ONLY: long_string_length
  USE parac,                           ONLY: paral
  USE readsr_utils,                    ONLY: readsi
  USE system,                          ONLY: maxsys

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: velocitinp

CONTAINS

  ! ==================================================================
  SUBROUTINE velocitinp(iunit)
    ! ==--------------------------------------------------------------==
    ! ==  READS VELOCITIES INPUT FOR INITIAL STEP (Section &ATOMS)    ==
    ! ==--------------------------------------------------------------==
    ! ==                                                              ==
    ! ==  VELOCITIES                                                  ==
    ! ==    nvel ia ib ic .... iz                                     ==
    ! ==      vx vy vz                                                ==
    ! ==      ....                                                    ==
    ! ==  END VELOCITIES                                              ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: iunit

    CHARACTER(*), PARAMETER                  :: procedureN = 'velocitinp'

    CHARACTER(:), ALLOCATABLE                :: long_line
    CHARACTER(LEN=long_string_length)        :: line
    INTEGER                                  :: i, ia, ierr, inin, iout, is, &
                                                j, nvel
    INTEGER, ALLOCATABLE                     :: lvel(:)
    LOGICAL                                  :: erread

    ALLOCATE(lvel(maxsys%nax*maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
10  CONTINUE
    long_line=""
    IF (paral%io_parent) THEN
       DO WHILE (.TRUE.)
          READ(iunit,ADVANCE="NO",EOR=5,ERR=20,END=20,fmt='(A)')line
          long_line=long_line(1:LEN_TRIM(long_line))//TRIM(line)
       END DO
5      long_line=long_line(1:LEN_TRIM(long_line))//TRIM(line)
    END IF
    IF (INDEX(long_line,'END').NE.0 .AND. INDEX(long_line,'VELOC').NE.0)&
         GOTO 30
    CALL readsi(long_line,1,iout,nvel,erread)
    IF (erread) GOTO 20
    inin=iout
    DO i=1,nvel
       CALL readsi(long_line,inin,iout,lvel(i),erread)
       IF (erread) GOTO 20
       inin=iout
    ENDDO
    DO i=1,nvel
       CALL igivenumbers(lvel(i),ia,is)
       IF (is.EQ.0)&
            CALL stopgm('VELOCITINP', 'BAD NUMBER OF ATOM',& 
            __LINE__,__FILE__)
       IF (paral%io_parent)&
            READ(iunit,err=20,END=20,fmt=*) (velp(j,ia,is),j=1,3)
       lvelini(0,is)=.TRUE.
       lvelini(ia,is)=.TRUE.
    ENDDO
    GOTO 10
30  CONTINUE
    DEALLOCATE(lvel,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
20  CONTINUE
    IF (paral%io_parent)&
         WRITE(6,*) ' VELOCITINP: ERROR WHILE READING ON UNIT ',iunit
    CALL stopgm('VELOCITINP',' ',& 
         __LINE__,__FILE__)
  END SUBROUTINE velocitinp
  ! ==================================================================
  SUBROUTINE igivenumbers(number,ia,is)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: number, ia, is

    INTEGER                                  :: num

! ==--------------------------------------------------------------==

    num=number
    IF (num.LE.0) THEN
       is=0
       RETURN
    ENDIF
    DO is=1,maxsys%nsx
       IF (num.LE.ions0%na(is)) THEN
          ia=num
          RETURN
       ELSE
          num=num-ions0%na(is)
       ENDIF
    ENDDO
    is=0
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE igivenumbers
  ! ==--------------------------------------------------------------==

END MODULE velocitinp_utils
