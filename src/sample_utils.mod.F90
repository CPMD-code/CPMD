MODULE sample_utils
  USE error_handling,                  ONLY: stopgm
  USE fileopenmod,                     ONLY: fo_info
  USE filnmod,                         ONLY: filbod
  USE machine,                         ONLY: m_sleep
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sync
  USE parac,                           ONLY: parai,&
                                             paral
  USE readsr_utils,                    ONLY: xstring
  USE store_types,                     ONLY: restart1
  USE system,                          ONLY: cnti,&
                                             cntr,&
                                             maxrf,&
                                             restf
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: sample_wait
  PUBLIC :: sample_go

CONTAINS

  ! ==================================================================
  SUBROUTINE sample_wait
    ! ==--------------------------------------------------------------==
    ! ==    STRUCTURE OF THE PSAMPLE FILE                             ==
    ! ==                                                              ==
    ! ==    NOMORE    number of CP steps                              ==
    ! ==    DELT_ELEC CP time step (+ or -)                           ==
    ! ==    FILENAME  full name of the restart file                   ==
    ! ==    FBODY     name for new restart files (FBODY.1 FBODY.2     ==
    ! ==              will be the names of the new restart files)     ==
    ! ==    NRESTF    number of new restart files                     ==
    ! ==    S1 S2 ... time steps at which the new restart files are   ==
    ! ==              written                                         ==
    ! ==                                                              ==
    ! ==--------------------------------------------------------------==
    CHARACTER(len=1), DIMENSION(10), PARAMETER :: &
      filnum = (/'1','2','3','4','5','6','7','8','9','0'/)
    INTEGER, PARAMETER                       :: iunit = 99 

    CHARACTER(len=80)                        :: fbody, filen, filename
    INTEGER                                  :: i, i1, i2, ia, icalls, ie, &
                                                isub, nw
    LOGICAL                                  :: test

    CALL tiset('SAMPLEWAIT',isub)
    IF (paral%parent) THEN
       filen=fo_info%fpath(fo_info%iapath:fo_info%iepath)//'PSAMPLE'
       IF (paral%io_parent)&
            WRITE(6,*) 'SAMPLE_WAIT| WAIT FOR CONTINUE-FILE ', filen
1      CONTINUE
       IF (paral%io_parent)&
            INQUIRE(file=filen,exist=test)
       IF (test) THEN
          IF (paral%io_parent)&
               OPEN(iunit,file=filen)
          IF (paral%io_parent)&
               REWIND(iunit)
          IF (paral%io_parent)&
               READ(iunit,*) cnti%nomore
          IF (paral%io_parent)&
               READ(iunit,*) cntr%delt_elec
          IF (paral%io_parent)&
               READ(iunit,*) filename
          IF (paral%io_parent)&
               READ(iunit,*) fbody
          IF (paral%io_parent)&
               READ(iunit,*) cnti%nrestf
          IF (cnti%nrestf.GT.maxrf)&
               CALL stopgm('SAMPLE_WAIT','TOO MANY RESTART FILES',& 
               __LINE__,__FILE__)
          IF (paral%io_parent)&
               READ(iunit,*) (restf%nstepwr(i),i=1,cnti%nrestf)
          IF (paral%io_parent)&
               CLOSE(iunit)
          IF (paral%io_parent)&
               WRITE(6,'(A)') ' SAMPLE_WAIT| CONTINUE CALCULATION'
          IF (paral%io_parent)&
               WRITE(6,'(A,T56,I10)')&
               ' SAMPLE_WAIT| NUMBER OF TIME STEPS ',cnti%nomore
          IF (paral%io_parent)&
               WRITE(6,'(A,T56,F10.1)')&
               ' SAMPLE_WAIT| NEW TIME STEP IS ',cntr%delt_elec
          CALL xstring(filename,ia,ie)
          IF (paral%io_parent)&
               WRITE(6,'(A,T50,A)')&
               ' SAMPLE_WAIT| RESTART FILE NAME',filename(ia:ie)
          CALL xstring(fbody,i1,i2)
          IF (paral%io_parent)&
               WRITE(6,'(A,T56,I10)')&
               ' SAMPLE_WAIT| NUMBER OF NEW RESTART FILES',cnti%nrestf
          DO i=1,cnti%nrestf
             IF (paral%io_parent)&
                  WRITE(6,'(A,T20,A,A,T56,I10)')&
                  ' SAMPLE_WAIT| FILE ',fbody(i1:i2)//filnum(i),&
                  ' WRITTEN AT STEP',restf%nstepwr(i)
          ENDDO
          ! ..prepare restart files
          restart1%rlate=.TRUE.
          nw=2
          IF (paral%io_parent)&
               OPEN(unit=nw,file='LATEST',status='UNKNOWN')
          IF (paral%io_parent)&
               REWIND(nw)
          IF (paral%io_parent)&
               WRITE(nw,'(A)') filename(ia:ie)
          icalls=0
          IF (paral%io_parent)&
               WRITE(nw,*) icalls
          IF (paral%io_parent)&
               CLOSE(nw)
          filbod=fbody(i1:i2)
          GOTO 3
       ENDIF
       IF (paral%io_parent)&
            INQUIRE(file='EXIT',exist=test)
       IF (test) cnti%nomore=-1
       IF (test) GOTO 3
       CALL m_sleep(1)
       GOTO 1
3      CONTINUE
    ENDIF
    CALL mp_sync(parai%allgrp)
    CALL mp_bcast(cnti%nomore,parai%source,parai%allgrp)
    CALL mp_bcast(cntr%delt_elec,parai%source,parai%allgrp)
    cntr%delt_ions=cntr%delt_elec
    CALL mp_bcast(cnti%nrestf,parai%source,parai%allgrp)
    CALL mp_bcast(restf%nstepwr,cnti%nrestf,parai%source,parai%allgrp)
    CALL tihalt('SAMPLEWAIT',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE sample_wait
  ! ==================================================================
  SUBROUTINE sample_go
    ! ==--------------------------------------------------------------==
    CHARACTER(len=80)                        :: filen, line
    INTEGER                                  :: ia, ie
    LOGICAL                                  :: test

! ==--------------------------------------------------------------==

    IF (paral%parent) THEN
       filen=fo_info%fpath(fo_info%iapath:fo_info%iepath)//'PSAMPLE'
       IF (paral%io_parent)&
            INQUIRE(file=filen,exist=test)
       IF (test) THEN
          CALL xstring(filen,ia,ie)
          IF (paral%io_parent)&
               WRITE(line,'(A,A)') 'rm -f ',filen(ia:ie)
          ! TODO: fix this: 'use system' and 'call system' conflicts 
          PRINT *, 'sample.F: SYSTEM call should be used here'
          STOP
       ENDIF
       IF (paral%io_parent)&
            WRITE(6,*) ' SAMPLE_GO| PASS CONTROL BACK TO SAMPLING PROGRAM'
    ENDIF
    CALL mp_sync(parai%allgrp)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE sample_go
  ! ==================================================================




END MODULE sample_utils
