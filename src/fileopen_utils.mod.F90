MODULE fileopen_utils
  USE error_handling,                  ONLY: stopgm
  USE fileopenmod,                     ONLY: &
       fo_app, fo_debug, fo_info, fo_mark, fo_nochk, fo_nofp, fo_scratch, &
       fo_stat_max, fo_ufo, fo_verb
  USE parac,                           ONLY: paral
  USE readsr_utils,                    ONLY: xstring

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_fileopen
  PUBLIC :: fileopen
  PUBLIC :: fileclose

CONTAINS

  ! ==================================================================
  ! INITIALIZE ARRAYS FOR THE FILEOPEN CALLS.
  SUBROUTINE init_fileopen
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: i

    fo_info%fo_tdebug=.FALSE.
    fo_info%iapath=1
    fo_info%iepath=2
    fo_info%fpath='./                            '
    DO i=1,fo_stat_max
       fo_info%fo_stat(i)=.FALSE.
    ENDDO
    ! SOME DEFAULT UNITS ARE ALREADY FLAGGED OPEN.
    fo_info%fo_stat(5)=.TRUE.
    fo_info%fo_stat(6)=.TRUE.
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE init_fileopen
  ! ==================================================================
  ! see fileopen.inc for documentation of the API.
  ! ==================================================================
  SUBROUTINE fileopen(iunit,filen,fo_mode,error)
    INTEGER                                  :: iunit
    CHARACTER(len=*), OPTIONAL               :: filen
    INTEGER                                  :: fo_mode
    LOGICAL                                  :: error

    CHARACTER(len=12)                        :: fo_form
    CHARACTER(len=255)                       :: fo_flags, fullpath
    CHARACTER(len=8), DIMENSION(0:10) :: fo_mode_string = (/'FO_DEF  ',&
      'FO_NEW  ','FO_OLD  ','FO_APP  ','FO_UFO  ','FO_VERB ','FO_MARK ',&
      'FO_NOCHK','FO_NOFP ','FO_DEBUG','FO_SCRA '/)
    INTEGER                                  :: i, ia, ia1, ia2, ie, ie1, &
                                                ie2, openmode
    LOGICAL                                  :: fexist, tfodbg, tfoscra, &
                                                tfovrb

! ==--------------------------------------------------------------==

    error=.FALSE.
    fexist=.FALSE.
    tfovrb=.FALSE.
    tfoscra=.FALSE.
    tfodbg=fo_info%fo_tdebug
    ! ========================================================================

222 FORMAT(1X,A,I4)

    ! SCRATCH FILE
    IF (MOD(fo_mode/fo_scratch,2).NE.0) tfoscra=.TRUE.
    IF (tfoscra) THEN
       IF (.NOT.PRESENT(filen)) THEN
          ! we should never get here
          CALL stopgm('FILEOPEN','FILENAME CAN BE OMITTED ONLY FOR SCRATCH FILES!',& 
               __LINE__,__FILE__)
       END IF
    END IF
    CALL xstring(filen,ia1,ie1)

    ! BE VERY VERBOSE, IF REQUESTED
    IF (MOD(fo_mode/fo_debug,2).NE.0) tfodbg=.TRUE.
    IF (tfodbg) THEN
       WRITE(fo_flags,'(A)') fo_mode_string(MOD(fo_mode,4))
       CALL xstring(fo_flags,ia,ie)
       DO i=4,10
          IF (MOD(fo_mode/(2**(i-2)),2).NE.0) THEN
             WRITE(fo_flags(ie+1:255),'(2A)') '+',fo_mode_string(i)
             CALL xstring(fo_flags,ia,ie)
          ENDIF
       ENDDO
       WRITE(6,'(A,I4,5A)') ' FILEOPEN: OPEN UNIT=',iunit,&
            ' FLAGS=',fo_flags(1:ie),' FILENAME=',filen(ia1:ie1),'"'
    ENDIF

    ! SCRATCH FILE
    IF (tfoscra.AND.tfodbg) WRITE(6,*) 'FILEOPEN: SCRATCH FLAG'

    ! BE VERBOSE, IF REQUESTED
    IF (MOD(fo_mode/fo_verb,2).NE.0) THEN
       IF (tfodbg) WRITE(6,*) 'FILEOPEN: VERBOSE FLAG'
       tfovrb=.TRUE.
    ENDIF

    ! do not use FILEPATH keyword (or CPMD_FILEPATH environment variable)
    IF (MOD(fo_mode/fo_nofp,2).NE.0) THEN
       IF (tfodbg) WRITE(6,*) 'FILEOPEN: DO NOT USE FPATH'
       fullpath=filen(ia1:ie1)
    ELSE
       fullpath=fo_info%fpath(fo_info%iapath:fo_info%iepath)//filen(ia1:ie1)
    ENDIF
    CALL xstring(fullpath,ia2,ie2)
    IF (tfodbg) WRITE(6,*) 'FILEOPEN: FULLPATH=',fullpath(ia2:ie2)

    ! CHECK WHETHER THE CHANNEL IS ALREADY OPEN
    IF (MOD(fo_mode/fo_nochk,2).EQ.0) THEN
       IF (tfodbg) WRITE(6,222) 'FILEOPEN: CHECK FOR OPEN UNIT:',iunit
       IF (iunit.GT.0.AND.iunit.LT.fo_stat_max) THEN
          IF (fo_info%fo_stat(iunit)) THEN
             WRITE(6,'(A,I4,2A)') ' FILEOPEN: UNIT ',iunit,&
                  ' ALREADY IN USE WHEN OPENING FILE ',filen(ia1:ie1)
             error=.TRUE.
             RETURN
          ENDIF
       ENDIF
    ELSE
       IF (tfodbg) WRITE(6,*) 'FILEOPEN: UNIT OPEN CHECK SUPPRESSED'
    ENDIF

    ! CHECK WHETHER THE FILE ALREADY EXISTS AND ADJUST THE OPEN MODE IF NOT.
    openmode=MOD(fo_mode,4)
    IF (tfodbg) WRITE(6,222) 'FILEOPEN: OPENMODE:',openmode
    IF ((MOD(fo_mode/fo_nochk,2).EQ.0).AND.(openmode.EQ.fo_app)) THEN
       INQUIRE(file=fullpath(ia2:ie2),exist=fexist)
       IF (.NOT.fexist) openmode=0
       IF (tfodbg) WRITE(6,222) 'FILEOPEN: OPENMODE NOW:',openmode
    ENDIF

    ! SET FORMATTED/UNFORMATTED
    fo_form='FORMATTED'
    IF (MOD(fo_mode/fo_ufo,2).NE.0) fo_form='UNFORMATTED'
    ! NOW OPEN FILE
    IF (tfoscra) THEN
       IF (tfodbg) WRITE(6,*) 'FILEOPEN: OPEN FILE WITH STATUS=SCRATCH'
       OPEN(unit=iunit,status='SCRATCH',form=fo_form,err=999)
    ELSE
       IF (openmode.EQ.0) THEN
          IF (tfodbg) WRITE(6,*) 'FILEOPEN: OPEN FILE WITH STATUS=UNKNOWN'
          OPEN(unit=iunit,file=fullpath(ia2:ie2),status='UNKNOWN',&
               form=fo_form,err=999)
       ELSEIF (openmode.EQ.1) THEN
          IF (tfodbg) WRITE(6,*) 'FILEOPEN: OPEN FILE WITH STATUS=NEW'
          OPEN(unit=iunit,file=fullpath(ia2:ie2),status='NEW',&
               form=fo_form,err=999)
       ELSEIF (openmode.EQ.2) THEN
          IF (tfodbg) WRITE(6,*) 'FILEOPEN: OPEN FILE WITH STATUS=OLD'
          OPEN(unit=iunit,file=fullpath(ia2:ie2),status='OLD',&
               form=fo_form,err=999)
       ELSEIF (openmode.EQ.3) THEN
          IF (tfodbg) WRITE(6,*) 'FILEOPEN: OPEN FILE FOR APPENDING'
          OPEN(unit=iunit,file=fullpath(ia2:ie2),status='OLD',&
               position='APPEND',&
               form=fo_form,err=999)
          IF (tfovrb) WRITE(6,'(3A)') ' FILE ',filen(ia1:ie1),&
               ' EXISTS, NEW DATA WILL BE APPENDED'
       ELSE
          ! we should never get here
          CALL stopgm('FILEOPEN','UNKNOWN FILE ACCESS MODE',& 
               __LINE__,__FILE__)
       ENDIF
    END IF
    ! RECORD THAT WE HAVE OPENED THIS UNIT
    IF (iunit.LT.fo_stat_max) THEN
       fo_info%fo_stat(iunit)=.TRUE.
    ENDIF

    ! WRITE NEW DATA MARK
    IF (openmode.EQ.3.AND.MOD(fo_mode/fo_mark,2).NE.0) THEN
       IF (tfodbg) WRITE(6,*) 'FILEOPEN: WRITE NEW DATA MARK'
       WRITE(iunit,'(A)') '   <<<<<<  NEW DATA  >>>>>>'
    ENDIF
    IF (tfodbg) WRITE(6,*) 'FILEOPEN: --------------------'
    RETURN

    ! FILE OPEN ERROR
999 CONTINUE
    IF (tfodbg) WRITE(6,*) 'FILEOPEN: ERROR WHEN OPENING FILE'
    error=.TRUE.
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE fileopen
  ! ==================================================================


  ! ==================================================================
  SUBROUTINE fileclose(iunit)
    INTEGER                                  :: iunit

    IF (paral%io_parent) CLOSE(iunit)
    ! RECORD THAT WE HAVE CLOSED THIS UNIT
    IF (paral%io_parent.AND.fo_info%fo_tdebug)&
         WRITE(6,'(A,I4)') ' FILECLOSE: CLOSING UNIT:',iunit
    IF (iunit.LT.fo_stat_max) THEN
       fo_info%fo_stat(iunit)=.FALSE.
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE fileclose
  ! ==================================================================

END MODULE fileopen_utils
