MODULE restart_p_utils
  USE coor,                            ONLY: tau0
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_def,&
                                             fo_old,&
                                             fo_ufo
  USE kinds,                           ONLY: int_8,&
                                             real_8
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum
  USE parac,                           ONLY: parai,&
                                             paral
  USE readsr_utils,                    ONLY: xstring
  USE response_pmod,                   ONLY: lower_left,&
                                             response1,&
                                             response_read,&
                                             response_write
  USE system,                          ONLY: maxsys,&
                                             ncpw,&
                                             parm,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
!!use wr30wfn_utils, only : rd30wfn
!!use wr30wfn_utils, only : wr30wfn
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: restart_p
  PUBLIC :: restart_nmr
  PUBLIC :: restart_epr

CONTAINS

  ! ==================================================================
  SUBROUTINE restart_p(c1,tag,nstate,ireadwrite)
    ! ==--------------------------------------------------------------==
    ! == Read (1) or write (2) one wavefunction in Version 3.0 format ==
    ! == into a file: RESTART.tag.                                    ==
    ! ==                                                              ==
    ! == IF THE FILE DOES NOT EXIST (for the read case), then the     ==
    ! == wavefunction is just put to zero. No error is issued.        ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c1(ncpw%ngw,*)
    CHARACTER(len=*)                         :: tag
    INTEGER                                  :: nstate, ireadwrite

    CHARACTER(len=128)                       :: filename
    INTEGER                                  :: controlbyte, fileunit, ia, &
                                                ie, ierror, info, ir, &
                                                irecord, isub, ngws_, nstate_
    INTEGER(int_8)                           :: fpos
    LOGICAL                                  :: error
    REAL(real_8)                             :: a1_(3), a2_(3), a3_(3)

    CALL tiset(' RESTART_P',isub)

    ! Construct the filename
    fileunit = 2602
    CALL xstring(tag,ia,ie)
    filename = 'RESTART.'//tag(ia:ie)
    CALL xstring(filename,ia,ie)

    ! ==--------------------------------------------------------------==
    IF (ireadwrite .EQ. response_write) THEN ! WRITE mode
       IF (paral%parent) THEN
          error=.FALSE.
          IF (paral%io_parent)&
               CALL fileopen(fileunit,filename(ia:ie),fo_def+fo_ufo,error)
          IF (paral%io_parent)&
               REWIND(fileunit,err=401)
          GOTO 402
401       CONTINUE
          error=.TRUE.
402       CONTINUE
          IF (error) CALL stopgm('RESTART_P',&
               'ERROR OPENING RESTART FILE ('&
               //filename(ia:ie)//') for writing.',& 
               __LINE__,__FILE__)

          controlbyte = 260273
          IF (paral%io_parent)&
               WRITE (fileunit) controlbyte
          IF (paral%io_parent)&
               WRITE (fileunit) parm%a1,parm%a2,parm%a3
          IF (paral%io_parent)&
               WRITE (fileunit) nstate,spar%ngws
       ENDIF                 ! parent
       ! ==--------------------------------------------------------------==
       ierror=0
       CALL wr30wfn(fileunit,ierror,nstate,c1,tau0,'NIL',irecord)
       IF (ierror.NE.0 .AND. paral%parent) THEN
          IF (paral%io_parent)&
               WRITE (6,'(A,I4,A)') 'RESTART_P: ERROR ',ierror,&
               ' WRITING RESTART_P file! File is corrupt!'
       ENDIF
       ! ==--------------------------------------------------------------==
       IF (paral%parent) THEN
          CALL xstring(filename,ia,ie)
          IF (paral%io_parent)&
               WRITE (6,'(A)')&
               'RESPONSE wavefunction written to RESTART file '&
               //filename(ia:ie)//'.'
          IF (paral%io_parent)&
               CALL fileclose(fileunit)
       ENDIF
       ! ==--------------------------------------------------------------==
    ELSEIF (ireadwrite .EQ. response_read) THEN ! READ mode
       ! ==--------------------------------------------------------------==
       IF (paral%parent) THEN
          error=.FALSE.
          IF (paral%io_parent)&
               CALL fileopen(fileunit,filename(ia:ie),fo_old+fo_ufo,error)
          IF (error) CALL stopgm('RESTART_P',&
               'RESTART FILE NOT FOUND: '//filename(ia:ie)//'.',& 
               __LINE__,__FILE__)
          ! if the file _does_ exist: read it
          IF (paral%io_parent)&
               REWIND(fileunit,err=501)
          error=.FALSE.
          GOTO 502
501       CONTINUE
          error=.TRUE.
502       CONTINUE
          IF (error) THEN
             CALL stopgm('RESTART_P',&
                  'ERROR OPENING RESTART FILE ('&
                  //filename(ia:ie)//') for reading.',& 
                  __LINE__,__FILE__)
          ENDIF
          ! ==--------------------------------------------------------------==
          ! GENERAL PART of restart file:
          IF (paral%io_parent)&
               READ(fileunit) controlbyte
          IF (paral%io_parent)&
               READ(fileunit) a1_,a2_,a3_
          IF (paral%io_parent)&
               READ(fileunit) nstate_,ngws_
          error=.FALSE.
          IF (controlbyte .NE. 260273)&
               error = .TRUE.
          DO ir=1,3
             IF ((parm%a1(ir) .NE. a1_(ir))&
                  .OR. (parm%a2(ir) .NE. a2_(ir))&
                  .OR. (parm%a3(ir) .NE. a3_(ir)))&
                  error = .TRUE.
          ENDDO
          IF ((nstate_.NE.nstate).OR.(spar%ngws.NE.ngws_))&
               error = .TRUE.
          IF (error) THEN
             CALL stopgm('RESTART_P',&
                  'ERROR READING RESTART FILE ('&
                  //filename(ia:ie)//'): wrong parameters.',& 
                  __LINE__,__FILE__)
          ENDIF
          ! ==--------------------------------------------------------------==
          IF (paral%io_parent)&
               READ(fileunit) irecord,fpos ! dummy read.
       ENDIF                 ! parent
       CALL rd30wfn(fileunit,c1,nstate,tau0,info,&
            .FALSE.,nstate,1,'NIL',fpos)
       ! ==--------------------------------------------------------------==
       IF (paral%parent) THEN
          CALL xstring(filename,ia,ie)
          IF (paral%io_parent)&
               WRITE (6,'(A)')&
               'RESPONSE wavefunction read from RESTART file '&
               //filename(ia:ie)//'.'
          IF (paral%io_parent)&
               CALL fileclose(fileunit)
       ENDIF
       response1%t_initialguess_c1 = .TRUE.
    ELSE                      ! read-write mode
       CALL stopgm('RESTART_P','BAD RESTART MODE.',& 
            __LINE__,__FILE__)
    ENDIF                     ! read - write
    ! ==--------------------------------------------------------------==
    CALL tihalt(' RESTART_P',isub)
    RETURN
  END SUBROUTINE restart_p
  ! ==================================================================





  ! ==================================================================
  ! The SAVE/RESTORE combination needs a little explanation for 
  ! use in parallel. When several processors keep only their parts of
  ! the shift and chi matrices, a cntl%proper save/restore is not easy.
  ! Instead, in this case, a global sum is done at each save.
  ! 
  ! When restoring, only the parent reads the file, whereas the other
  ! processors assume to have initially empty (zero) shift
  ! and chi matrices.

  ! ==================================================================
  SUBROUTINE restart_nmr(simple_done,full_done,&
       shift_sumup,chi_sumup,&
       shift_saved,chi_saved,nstate,ireadwrite)

    ! Arguments:
    LOGICAL                                  :: simple_done(6)
    REAL(real_8) :: shift_sumup(3,3,maxsys%nax,maxsys%nsx), chi_sumup(3,3), &
      shift_saved(3,3,maxsys%nax,maxsys%nsx), chi_saved(3,3)
    INTEGER                                  :: nstate
    LOGICAL                                  :: full_done(3,nstate)
    INTEGER                                  :: ireadwrite

    CHARACTER(*), PARAMETER                  :: procedureN = 'restart_nmr'

    CHARACTER(len=128)                       :: filename
    INTEGER                                  :: cbyte, controlbyte, fileunit, &
                                                ia, ie, ierr, ir, is, ll_(3), &
                                                ngws_, nstate_
    LOGICAL                                  :: error, errors
    REAL(real_8)                             :: a1_(3), a2_(3), a3_(3)
    REAL(real_8), ALLOCATABLE                :: scr(:,:), scr2(:,:)

    fileunit=2604
    controlbyte=19732604
    filename='RESTART.NMR'
    CALL xstring(filename,ia,ie)
    ALLOCATE(scr(3*3*maxsys%nax*maxsys%nsx,2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(scr2(3*3,2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    IF (ireadwrite .EQ. response_write) THEN ! WRITE MODE
       ! ==--------------------------------------------------------------==
       ! Preparation:
       CALL dcopy(3*3*maxsys%nax*maxsys%nsx,shift_saved,1,scr(1,1),1)
       CALL dcopy(3*3*maxsys%nax*maxsys%nsx,shift_sumup,1,scr(1,2),1)
       CALL dcopy(3*3,chi_saved,1,scr2(1,1),1)
       CALL dcopy(3*3,chi_sumup,1,scr2(1,2),1)
       CALL mp_sum(scr2,2*3*3,parai%allgrp)
       CALL mp_sum(scr,2*3*3*maxsys%nax*maxsys%nsx,parai%allgrp)

       IF (paral%parent) THEN
          error = .FALSE.
          IF (paral%io_parent)&
               CALL fileopen(fileunit,filename(ia:ie),fo_def+fo_ufo,error)
          IF (paral%io_parent)&
               REWIND(fileunit)
          IF (paral%io_parent)&
               WRITE (fileunit) controlbyte
          IF (paral%io_parent)&
               WRITE (fileunit) parm%a1,parm%a2,parm%a3
          IF (paral%io_parent)&
               WRITE (fileunit) nstate,spar%ngws
          DO is=1,nstate
             IF (paral%io_parent)&
                  WRITE (fileunit) (lower_left(ir,is),ir=1,3)
          ENDDO
          IF (paral%io_parent)&
               WRITE (fileunit) simple_done
          IF (paral%io_parent)&
               WRITE (fileunit) full_done
          IF (paral%io_parent)&
               WRITE (fileunit) scr
          IF (paral%io_parent)&
               WRITE (fileunit) scr2
          IF (paral%io_parent)&
               CALL fileclose(fileunit)
       ENDIF                 ! Parent
       ! ==--------------------------------------------------------------==
    ELSEIF (ireadwrite .EQ. response_read) THEN ! READ MODE
       ! ==--------------------------------------------------------------==
       ! OPEN the file
       IF (paral%parent) THEN
          error = .FALSE.
          IF (paral%io_parent)&
               CALL fileopen(fileunit,filename(ia:ie),fo_old+fo_ufo,error)
          IF (error) THEN
             CALL stopgm('RESTART_P',&
                  'ERROR OPENING RESTART FILE ('&
                  //filename(ia:ie)//'.',&
                  __LINE__,__FILE__)
          ENDIF
       ENDIF
       ! ==--------------------------------------------------------------==
       ! READ the file and check control parameters
       IF (paral%parent) THEN
          IF (paral%io_parent)&
               REWIND(fileunit)
          errors = .FALSE.
          IF (paral%io_parent)&
               READ(fileunit) cbyte
          IF (paral%io_parent)&
               READ(fileunit) a1_,a2_,a3_
          IF (paral%io_parent)&
               READ(fileunit) nstate_,ngws_
          error =  (cbyte .NE. controlbyte)
          errors = errors .OR. error
          IF ((error).AND.paral%io_parent)&
               WRITE(6,*) 'CONTROL BYTE CHANGED.'
          error =  (nstate_ .NE. nstate)
          errors = errors .OR. error
          IF ((error).AND.paral%io_parent)&
               WRITE(6,*) 'NSTATE CHANGED.'
          error =  (ngws_ .NE. spar%ngws)
          errors = errors .OR. error
          IF ((error).AND.paral%io_parent)&
               WRITE(6,*) 'NGWS CHANGED.'
          DO ir=1,3
             error = (ABS(parm%a1(ir)-a1_(ir)).GE.1.e-9_real_8)&
                  .OR. (ABS(parm%a2(ir)-a2_(ir)).GE.1.e-9_real_8)&
                  .OR. (ABS(parm%a3(ir)-a3_(ir)).GE.1.e-9_real_8)
             errors = errors .OR. error
             IF ((error).AND.paral%io_parent)&
                  WRITE(6,*) 'UNIT CELL CHANGED.'
          ENDDO
          DO is=1,nstate
             IF (paral%io_parent)&
                  READ(fileunit) (ll_(ir),ir=1,3)
             DO ir=1,3
                error = (ll_(ir) .NE. lower_left(ir,is))
                errors = errors .OR. error
                IF ((error).AND.paral%io_parent)&
                     WRITE(6,*) 'LLC CHANGED.'
             ENDDO
          ENDDO
       ENDIF                 ! parent
       CALL mp_bcast(errors,parai%source,parai%allgrp)
       IF (errors) THEN
          CALL stopgm('RESTART_P',&
               'ERROR READING RESTART FILE ('&
               //filename(ia:ie)//': parameters changed! !',&
               __LINE__,__FILE__)
       ENDIF
       ! ==--------------------------------------------------------------==
       ! READ the values calculated so far
       IF (paral%parent) THEN
          IF (paral%io_parent)&
               READ(fileunit) simple_done
          IF (paral%io_parent)&
               READ(fileunit) full_done
          IF (paral%io_parent)&
               READ(fileunit) scr
          IF (paral%io_parent)&
               READ(fileunit) scr2
          IF (paral%io_parent)&
               CALL fileclose(fileunit)
          CALL dcopy(3*3*maxsys%nax*maxsys%nsx,scr(1,1),1,shift_saved,1)
          CALL dcopy(3*3*maxsys%nax*maxsys%nsx,scr(1,2),1,shift_sumup,1)
          CALL dcopy(3*3,scr2(1,1),1,chi_saved,1)
          CALL dcopy(3*3,scr2(1,2),1,chi_sumup,1)
       ELSE
          CALL zeroing(shift_saved)!,3*3*maxsys%nax*maxsys%nsx)
          CALL zeroing(shift_sumup)!,3*3*maxsys%nax*maxsys%nsx)
          CALL zeroing(chi_saved)!,3*3)
          CALL zeroing(chi_sumup)!,3*3)
       ENDIF
       CALL mp_bcast(simple_done,SIZE(simple_done),parai%source,parai%allgrp)
       CALL mp_bcast(full_done,SIZE(full_done),parai%source,parai%allgrp)
       IF (paral%parent) THEN
          IF (paral%io_parent)&
               WRITE (6,'(A)')&
               'NMR RESTART DATA read from file '&
               //filename(ia:ie)//'.'
       ENDIF
    ELSE                      ! READ/WRITE MODES.
       CALL stopgm('RESTART_NMR','BAD RESTART MODE.',& 
            __LINE__,__FILE__)
    ENDIF                     ! read - write

    DEALLOCATE(scr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(scr2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    RETURN
  END SUBROUTINE restart_nmr
  ! ==================================================================


  ! ==================================================================
  SUBROUTINE restart_epr(simple_done,nstate,ireadwrite)
    ! Arguments:
    LOGICAL                                  :: simple_done(9)
    INTEGER                                  :: nstate, ireadwrite

    CHARACTER(len=128)                       :: filename
    INTEGER                                  :: cbyte, controlbyte, fileunit, &
                                                ia, ie, ir, is, ll_(3), &
                                                ngws_, nstate_
    LOGICAL                                  :: error, errors
    REAL(real_8)                             :: a1_(3), a2_(3), a3_(3)

    fileunit=2604
    controlbyte=19732604
    filename='RESTART.EPR'
    CALL xstring(filename,ia,ie)

    IF (ireadwrite .EQ. response_write) THEN ! WRITE MODE
       ! ==--------------------------------------------------------------==
       IF (paral%parent) THEN
          IF (paral%io_parent)&
               OPEN(unit=fileunit,file=filename,&
               status='UNKNOWN',form='UNFORMATTED')
          IF (paral%io_parent)&
               REWIND(fileunit)
          IF (paral%io_parent)&
               WRITE (fileunit) controlbyte
          IF (paral%io_parent)&
               WRITE (fileunit) parm%a1,parm%a2,parm%a3
          IF (paral%io_parent)&
               WRITE (fileunit) nstate,spar%ngws
          DO is=1,nstate
             IF (paral%io_parent)&
                  WRITE (fileunit) (lower_left(ir,is),ir=1,3)
          ENDDO
          IF (paral%io_parent)&
               WRITE (fileunit) simple_done
          IF (paral%io_parent)&
               CLOSE(fileunit)
       ENDIF                 ! Parent
       ! ==--------------------------------------------------------------==
    ELSEIF (ireadwrite .EQ. response_read) THEN ! READ MODE
       ! ==--------------------------------------------------------------==
       ! OPEN the file
       IF (paral%parent) THEN
          error = .FALSE.
          IF (paral%io_parent)&
               OPEN(unit=fileunit,file=filename,&
               status='OLD',form='UNFORMATTED')
       ENDIF
       CALL mp_bcast(error,parai%source,parai%allgrp)
       IF (error) THEN
          CALL stopgm('RESTART_P',&
               'ERROR OPENING RESTART FILE ('&
               //filename(ia:ie)//'.',&
               __LINE__,__FILE__)
       ENDIF
       ! ==--------------------------------------------------------------==
       ! READ the file and check control parameters
       IF (paral%parent) THEN
          errors = .FALSE.
          IF (paral%io_parent)&
               READ(fileunit) cbyte
          IF (paral%io_parent)&
               READ(fileunit) a1_,a2_,a3_
          IF (paral%io_parent)&
               READ(fileunit) nstate_,ngws_
          error =  (cbyte .NE. controlbyte)
          errors = errors .OR. error
          IF ((error).AND.paral%io_parent)&
               WRITE(6,*) 'CONTROL BYTE CHANGED.'
          error =  (nstate_ .NE. nstate)
          errors = errors .OR. error
          IF ((error).AND.paral%io_parent)&
               WRITE(6,*) 'NSTATE CHANGED.'
          error =  (ngws_ .NE. spar%ngws)
          errors = errors .OR. error
          IF ((error).AND.paral%io_parent)&
               WRITE(6,*) 'NGWS CHANGED.'
          DO ir=1,3
             error = (ABS(parm%a1(ir)-a1_(ir)).GE.1.e-9_real_8)&
                  .OR. (ABS(parm%a2(ir)-a2_(ir)).GE.1.e-9_real_8)&
                  .OR. (ABS(parm%a3(ir)-a3_(ir)).GE.1.e-9_real_8)
             errors = errors .OR. error
             IF ((error).AND.paral%io_parent)&
                  WRITE(6,*) 'UNIT CELL CHANGED.'
          ENDDO
          DO is=1,nstate
             IF (paral%io_parent)&
                  READ(fileunit) (ll_(ir),ir=1,3)
             DO ir=1,3
                error = (ll_(ir) .NE. lower_left(ir,is))
                errors = errors .OR. error
                IF ((error).AND.paral%io_parent)&
                     WRITE(6,*) 'LLC CHANGED.'
             ENDDO
          ENDDO
       ENDIF             ! parent
       CALL mp_bcast(errors,parai%source,parai%allgrp)
       IF (errors) THEN
          CALL stopgm('RESTART_P',&
               'ERROR READING RESTART FILE ('&
               //filename(ia:ie)//': parameters changed! !',&
               __LINE__,__FILE__)
       ENDIF
       ! ==--------------------------------------------------------------==
       ! READ the values calculated so far
       IF (paral%parent) THEN
          IF (paral%io_parent)&
               READ(fileunit) simple_done
          IF (paral%io_parent)&
               CLOSE(fileunit)
       ENDIF
       CALL mp_bcast(simple_done,SIZE(simple_done),parai%source,parai%allgrp)
       IF (paral%parent) THEN
          CALL xstring(filename,ia,ie)
          IF (paral%io_parent)&
               WRITE (6,'(A)')&
               'EPR RESTART DATA read from file '&
               //filename(ia:ie)//'.'
          IF (paral%io_parent)&
               CLOSE(fileunit)
       ENDIF
    ELSE                  ! READ/WRITE MODES.
       CALL stopgm('RESTART_EPR','BAD RESTART MODE.',& 
            __LINE__,__FILE__)
    ENDIF                 ! read - write

    RETURN
  END SUBROUTINE restart_epr
  ! ==================================================================

END MODULE restart_p_utils
