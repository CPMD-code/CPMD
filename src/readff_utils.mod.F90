MODULE readff_utils
  USE clas,                            ONLY: clab,&
                                             clas3,&
                                             pote1,&
                                             pr12,&
                                             rcr12
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: paral
  USE readsr_utils,                    ONLY: readsr
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: readff
  !public :: get_tag
  !public :: cl_typ

CONTAINS

  ! ==================================================================
  SUBROUTINE readff(iunit)
    ! ==--------------------------------------------------------------==
    ! ==  Reads Classical Force Fields                                ==
    ! ==--------------------------------------------------------------==
    ! ==                                                              ==
    ! ==                                                              ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: iunit

    CHARACTER(*), PARAMETER                  :: procedureN = 'readff'

    CHARACTER(len=4)                         :: tag
    CHARACTER(len=80)                        :: line
    INTEGER                                  :: ia1, ia2, ierr, iout, ipos
    LOGICAL                                  :: erread
    REAL(real_8)                             :: para

! 
! ==--------------------------------------------------------------==

10  READ(iunit,err=20,END=20,fmt='(A)') line
    IF (INDEX(line,'END').NE.0 .AND. INDEX(line,'FORCE').NE.0)&
         GOTO 30
    ipos=1
    CALL get_tag(tag,4,line,ipos,80)
    IF (tag(1:3).EQ.'R12') THEN
       ! ..R12 potential
       IF (pote1%initr12.EQ.0) THEN
          pote1%t_r12=.TRUE.
          pote1%initr12=1
          ALLOCATE(pr12(clas3%ncltyp,clas3%ncltyp),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          CALL zeroing(pr12)!,clas3%ncltyp*clas3%ncltyp)
          ALLOCATE(rcr12(clas3%ncltyp,clas3%ncltyp),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          CALL zeroing(rcr12)!,clas3%ncltyp*clas3%ncltyp)
       ENDIF
       CALL get_tag(tag,4,line,ipos,80)
       ia1=cl_typ(tag)
       CALL get_tag(tag,4,line,ipos,80)
       ia2=cl_typ(tag)
       CALL readsr(line,ipos,iout,para,erread)
       pr12(ia1,ia2)=para
       pr12(ia2,ia1)=para
       ipos=iout
       CALL readsr(line,ipos,iout,para,erread)
       rcr12(ia1,ia2)=para
       rcr12(ia2,ia1)=para
    ELSEIF (tag(1:2).EQ.'LJ') THEN
       ! ..LJ potential
    ENDIF
    GOTO 10
20  CONTINUE
    IF (paral%io_parent)&
         WRITE(6,*) ' ERROR WHILE READING FORCE FIELD '
    CALL stopgm('READFF',' ',& 
         __LINE__,__FILE__)
30  CONTINUE
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE readff
  ! ==================================================================
  SUBROUTINE get_tag(tag,len,line,ipos,iend)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: len
    CHARACTER(len=*)                         :: tag, line
    INTEGER                                  :: ipos, iend

    INTEGER                                  :: i, l

    tag='           '
    DO i=ipos,iend
       IF (line(i:i).NE.' ') GOTO 10
    ENDDO
    ipos=iend+1
    RETURN
10  CONTINUE
    ipos=i
    l=0
    DO i=ipos,iend
       IF (line(i:i).NE.' ') THEN
          l=l+1
          IF (l.GT.len) RETURN
          tag(l:l)=line(i:i)
       ELSE
          ipos=i
          RETURN
       ENDIF
    ENDDO
    ipos=iend+1
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE get_tag
  ! ==================================================================
  FUNCTION cl_typ(tag)
    CHARACTER(len=4)                         :: tag
    INTEGER                                  :: cl_typ

    INTEGER                                  :: i

! ==--------------------------------------------------------------==

    cl_typ=0
    DO i=1,clas3%ncltyp
       IF (tag.EQ.clab(i)) THEN
          cl_typ=i
          RETURN
       ENDIF
    ENDDO
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent)&
         WRITE(6,*) 'READFF| Unknown atom label ',tag
    IF (paral%io_parent)&
         WRITE(6,'(A4)') tag
    DO i=1,clas3%ncltyp
       IF (paral%io_parent)&
            WRITE(6,'(A4)') clab(i)
    ENDDO
    CALL stopgm('CL_TYP','ERROR',& 
         __LINE__,__FILE__)
  END FUNCTION cl_typ
  ! ==================================================================

END MODULE readff_utils
