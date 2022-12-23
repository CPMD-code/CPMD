MODULE inscan_utils
  USE error_handling,                  ONLY: stopgm
  USE readsr_utils,                    ONLY: xstring
  USE store_types,                     ONLY: store1

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: inscan

CONTAINS

  ! ==================================================================
  FUNCTION inscan(iunit,label)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: iunit
    CHARACTER(len=*)                         :: label
    INTEGER                                  :: inscan

    CHARACTER(len=80)                        :: line
    INTEGER                                  :: ia, ie, ifind, iline
    LOGICAL                                  :: exists

! Variables
! ==--------------------------------------------------------------==
! ==  Scans unit 'iunit' for the section header 'label'           ==
! ==--------------------------------------------------------------==

    INQUIRE(unit=iunit,exist=exists)
    IF (.NOT.exists) THEN
       WRITE(6,'(A,I3,A)') ' INSCAN: UNIT',iunit,' DOES NOT EXIST'
       CALL stopgm('INSCAN','TRYING TO READ FROM A CLOSED FILE',& 
            __LINE__,__FILE__)
    ENDIF
    inscan=1
    iline=0
    REWIND(iunit)
    ! ==--------------------------------------------------------------==
    CALL xstring(label,ia,ie)
10  CONTINUE
    READ(iunit,END=20,err=20,fmt='(A80)') line
    iline=iline+1
    ifind=INDEX(line,label(ia:ie))
    IF (ifind.NE.0) THEN
       IF (store1%tdebio)WRITE(6,'(A,A,T40,A,I3,A,I5)')&
            ' INSCAN| FOUND SECTION: ',label(ia:ie),&
            ' FOR UNIT',iunit,' ON LINE:',iline
       inscan=0
       GOTO 20
    ENDIF
    GOTO 10
20  CONTINUE
    IF (store1%tdebio.AND.(inscan.NE.0)) WRITE(6,'(A,A,T40,A,I3)')&
         ' INSCAN|    NO SECTION: ',label(ia:ie),' FOR UNIT',iunit
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION inscan
  ! ==================================================================

END MODULE inscan_utils
