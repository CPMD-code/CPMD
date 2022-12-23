MODULE timer
  USE envj,                            ONLY: tjlimit
  USE kinds,                           ONLY: int_1,&
                                             int_2,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8
  USE machine,                         ONLY: m_cputime,&
                                             m_datum,&
                                             m_flush,&
                                             m_walltime
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_max,&
                                             mp_sum
  USE parac,                           ONLY: parai,&
                                             paral
  USE system,                          ONLY: cp_trace
  USE time,                            ONLY: &
       ctsta, ctsum, incl_wtime, lentim, ncall, ntimer, ntimx, time_in, &
       timeloop, timestop, tjsloop, tjstart, tname, wtsta, wtsum
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: tiset
  PUBLIC :: tihalt
  PUBLIC :: tipri
  PUBLIC :: tistart
  PUBLIC :: tilimex
  PUBLIC :: tilimit
  PUBLIC :: ttimp

CONTAINS

  ! ==================================================================
  SUBROUTINE tiset(sname,nstim)
    ! ==--------------------------------------------------------------==
    ! ==  START THE CPU AND WALL CLOCK TIME FOR ROUTINE SNAME         ==
    ! ==--------------------------------------------------------------==
    CHARACTER(*)                             :: sname
    INTEGER                                  :: nstim

    CHARACTER(*), PARAMETER                  :: procedureN = 'tiset'

    INTEGER                                  :: i
    REAL(real_8)                             :: ct, wt

    !$omp master

    tname%trace_depth = tname%trace_depth + 1
    IF (tname%trace_depth.GT.ntimx) THEN
       WRITE(*,*) "wrong trace_depth..."
       CALL m_flush(6)
       STOP
    ENDIF
    IF (LEN(ADJUSTL(sname))>lentim) THEN
       WRITE(*,*) "routine name ("//ADJUSTL(sname)//") too long!"
       CALL m_flush(6)
       STOP
    ENDIF
    tname%trace_names(tname%trace_depth) = ADJUSTL(sname)
    ! <<<<

    DO i=1,ntimer
       IF (sname.EQ.tname%tnam(i)) THEN
          nstim=i
          GOTO 100
       ENDIF
    ENDDO
    ntimer=ntimer+1
    IF (ntimer.GT.ntimx) THEN
       IF (paral%io_parent)&
            WRITE(6,*) ' PARAMETER NTIMX IN INCLUDE FILE time.inc TOO ',&
            'SMALL'
       CALL m_flush(6)
       STOP
       !       call stopgm('TISET',' PARAMETER',& 
       !            __LINE__,__FILE__)
    ENDIF
    nstim=ntimer
    tname%tnam(nstim)=sname
100 CONTINUE

    ct = m_cputime()
    wt = m_walltime()
    ctsta(nstim) = ct
    wtsta(nstim) = wt
    time_in(nstim,1) = ct
    time_in(nstim,2) = wt

    ! >>>> the trace
    ! we trace if needed
    IF (cp_trace%ttrace) THEN
       IF (tname%trace_depth.LE.tname%trace_max_depth.AND.&
            ncall(nstim)+1.LE.tname%trace_max_calls) THEN
          IF (cp_trace%ttrace_master_only) THEN
             IF (parai%cp_me.EQ.0) THEN
                IF (trace_this_procedure(&
                     tname%trace_names(tname%trace_depth))) THEN
                   WRITE(6,10) parai%cp_me,&
                        ncall(nstim)+1,&
                        incl_wtime(nstim),&
                        tname%trace_depth,&
                        REPEAT('..',tname%trace_depth),&
                        TRIM(tname%trace_names(tname%trace_depth))
                ENDIF
             ENDIF
          ELSE
             IF (trace_this_procedure(&
                  tname%trace_names(tname%trace_depth))) THEN
                WRITE(6,10) parai%cp_me,&
                     ncall(nstim)+1,&
                     incl_wtime(nstim),&
                     tname%trace_depth,&
                     REPEAT('..',tname%trace_depth),&
                     TRIM(tname%trace_names(tname%trace_depth))
             ENDIF
          ENDIF
          CALL m_flush(6)
10        FORMAT('#tr: ',i4,' calls: ',i6,&
               ' time: ',7x,f8.2,1x,i2,a,' >> ',2x,a)
       ENDIF
    ENDIF

    !$omp end master

  END SUBROUTINE tiset

  ! ==================================================================
  SUBROUTINE tihalt(sname,nstim)
    ! ==--------------------------------------------------------------==
    ! ==  STOP THE CPU AND WALL CLOCK TIME FOR ROUTINE SNAME          ==
    ! ==--------------------------------------------------------------==
    CHARACTER(*)                             :: sname
    INTEGER                                  :: nstim

    CHARACTER(*), PARAMETER                  :: procedureN = 'tihalt'

    CHARACTER(len=lentim)                    :: SNAME_tmp
    INTEGER                                  :: i
    REAL(real_8)                             :: ct, dct, dwt, wt

    !$omp master

    SNAME_tmp=sname
    IF (nstim.EQ.0) THEN
       DO i=1,ntimer
          IF (SNAME_tmp.EQ.tname%tnam(i)) THEN
             nstim=i
             GOTO 100
          ENDIF
       ENDDO
       WRITE(6,*) ' NO SUCH TIME HANDLE ',SNAME_tmp
       CALL m_flush(6)
       STOP
       !       call stopgm('TIHALT',' HANDLE   ',& 
       !            __LINE__,__FILE__)
100    CONTINUE
    ELSE
       IF (SNAME_tmp.NE.tname%tnam(nstim)) THEN
          WRITE(6,*) ' INCONSISTENT HANDLE IN TIHALT ',SNAME_tmp,&
               tname%tnam(nstim)
          CALL m_flush(6)
          STOP
          !          call stopgm('TIHALT',' HANDLE   ',& 
          !               __LINE__,__FILE__)
       ENDIF
    ENDIF

    ct = m_cputime()
    wt = m_walltime()
    dct=ct-ctsta(nstim)
    dwt=wt-wtsta(nstim)
    ctsum(nstim)=ctsum(nstim)+dct
    wtsum(nstim)=wtsum(nstim)+dwt
    incl_wtime(nstim) = incl_wtime(nstim) +&
         wt - time_in(nstim,2)

    ncall(nstim)=ncall(nstim)+1
    DO i=1,ntimer
       ctsta(i)=ctsta(i)+dct
       wtsta(i)=wtsta(i)+dwt
    ENDDO
    ! >>>> the trace
    ! we trace if needed
    IF (cp_trace%ttrace) THEN
       IF (tname%trace_depth.LE.tname%trace_max_depth.AND.&
            ncall(nstim).LE.tname%trace_max_calls) THEN
          IF (cp_trace%ttrace_master_only) THEN
             IF (parai%cp_me.EQ.0) THEN
                IF (trace_this_procedure(&
                     tname%trace_names(tname%trace_depth))) THEN
                   WRITE(6,10) parai%cp_me,&
                        ncall(nstim),&
                        wt-time_in(nstim,2),&
                        incl_wtime(nstim),&
                        tname%trace_depth,&
                        REPEAT('..',tname%trace_depth),&
                        TRIM(tname%trace_names(tname%trace_depth))
                ENDIF
             ENDIF
          ELSE
             IF (trace_this_procedure(&
                  tname%trace_names(tname%trace_depth))) THEN
                WRITE(6,10) parai%cp_me,&
                     ncall(nstim),&
                     wt-time_in(nstim,2),&
                     incl_wtime(nstim),&
                     tname%trace_depth,&
                     REPEAT('..',tname%trace_depth),&
                     TRIM(tname%trace_names(tname%trace_depth))
             ENDIF
          ENDIF
          CALL m_flush(6)
10        FORMAT('#tr: ',i4,' calls: ',i6,&
               ' time: ',f7.2,f8.2,1x,i2,a,' << ',2x,a)
       ENDIF
    ENDIF
    tname%trace_depth = tname%trace_depth - 1
    IF (tname%trace_depth.LT.0) THEN
       WRITE(*,*) "wrong trace_depth..."
       CALL m_flush(6)
       STOP
       !       call stopgm(procedureN,' wrong trace_depth...',& 
       !            __LINE__,__FILE__)
    ENDIF

    !$omp end master

  END SUBROUTINE tihalt

  ! ==================================================================
  LOGICAL FUNCTION trace_this_procedure(procedure_name)
    IMPLICIT NONE
    ! ==--------------------------------------------------------------==
    ! ==  decide to trace this procedure or not                       ==
    ! ==--------------------------------------------------------------==
    ! Arguments
    CHARACTER(*) :: procedure_name
    ! Variables
    INTEGER :: i
    ! ==--------------------------------------------------------------==
    IF (tname%trace_nbr_procedure.EQ.0)  THEN
       trace_this_procedure = .TRUE.
    ELSE
       trace_this_procedure = .FALSE.
       DO i=1,tname%trace_nbr_procedure
          IF (TRIM(procedure_name).EQ.TRIM(tname%trace_procedure(i))) THEN
             trace_this_procedure = .TRUE.
             EXIT
          ENDIF
       ENDDO
    ENDIF
  END FUNCTION trace_this_procedure

  ! ==================================================================
  SUBROUTINE tipri
    ! ==--------------------------------------------------------------==
    ! ==  PRINT TIMING INFO                                           ==
    ! ==--------------------------------------------------------------==
    ! Variables
    CHARACTER(*), PARAMETER                  :: procedureN = 'tipri'

    CHARACTER(lentim)                        :: timerN(ntimx)
    INTEGER                                  :: found, i, ierr, imax, j, &
                                                NTIMER_master
    REAL(real_8)                             :: cpmax, cpthrs, totalct, &
                                                totalwt
    REAL(real_8), ALLOCATABLE, DIMENSION(:)  :: avrg_excl_time, &
                                                avrg_incl_time, excl_time, &
                                                incl_time, max_excl_time, &
                                                max_incl_time

    IF (paral%io_parent)WRITE(6,'(/,/,1X,80("*"))')
    IF (paral%io_parent)WRITE(6,'(" *",78X,"*")')
    IF (paral%io_parent)WRITE(6,'(" *",35X,A,35X,"*")') ' TIMING '
    IF (paral%io_parent)WRITE(6,'(" *",78X,"*")')
    IF (paral%io_parent)WRITE(6,'(1X,80("*"))')
    !IF (paral%io_parent)WRITE(6,'(T2,A,T19,A,T37,A,T56,A)') 'SUBROUTINE',&
    IF (paral%io_parent)WRITE(6,'(T2,A,T35,A,T53,A,T72,A)') 'SUBROUTINE',&
         'CALLS','SELF TIME','TOTAL TIME'
    !IF (paral%io_parent)WRITE (6,'(T26,4A10)')&
    IF (paral%io_parent)WRITE (6,'(T42,4A10)')&
         'AVERAGE','MAXIMUM','AVERAGE','MAXIMUM'

    ! ==--------------------------------------------------------------==
    ! the io node drives which procedures are used for the average
    ALLOCATE(excl_time(ntimx),incl_time(ntimx),&
         avrg_excl_time(ntimx),avrg_incl_time(ntimx),&
         max_excl_time(ntimx),max_incl_time(ntimx),stat=ierr)
    IF (ierr.NE.0) STOP "Allocation error" !call stopgm(procedureN,'Allocation error',& 
    !         __LINE__,__FILE__)
    excl_time(:)=0.0_real_8; incl_time(:)=0.0_real_8
    NTIMER_master = ntimer
    timerN(:)=tname%tnam(:)
    CALL mp_bcast(NTIMER_master,parai%io_source,parai%cp_grp)
    CALL mp_bcast(timerN,ntimx,parai%io_source,parai%cp_grp)
    DO i=1,NTIMER_master
       found=0
       DO j=1,ntimer
          IF (timerN(i).EQ.tname%tnam(j)) THEN
             found=j
             EXIT
          ENDIF
       ENDDO
       IF (found.GT.0) THEN
          excl_time(i) = wtsum(j)*1.e-3_real_8
          incl_time(i) = incl_wtime(j)*1.e-3_real_8
       ENDIF
    ENDDO
    ! ==--------------------------------------------------------------==

    avrg_excl_time(:) = excl_time(:)
    max_excl_time(:)  = excl_time(:)
    avrg_incl_time(:) = incl_time(:)
    max_incl_time(:)  = incl_time(:)
    CALL mp_max(max_excl_time,ntimx,parai%cp_grp)
    CALL mp_max(max_incl_time,ntimx,parai%cp_grp)
    CALL mp_sum(avrg_excl_time,ntimx,parai%cp_grp)
    CALL mp_sum(avrg_incl_time,ntimx,parai%cp_grp)
    avrg_excl_time=avrg_excl_time/REAL(parai%cp_nproc,kind=real_8)
    avrg_incl_time=avrg_incl_time/REAL(parai%cp_nproc,kind=real_8)

    totalct=0._real_8
    totalwt=0._real_8
    cpthrs=0._real_8
    DO j=1,ntimx
       imax=0
       cpmax=-0.1_real_8
       DO i=1,ntimx
          IF (max_incl_time(i).GT.cpmax) THEN
             imax=i
             cpmax=max_incl_time(i)
          ENDIF
       ENDDO
       IF (j.EQ.1) THEN
          cpthrs=cpmax/1000._real_8
       ENDIF
       IF (cpmax.LT.cpthrs) THEN
          GOTO 100
       ENDIF
       i=imax
       !       IF (paral%io_parent)WRITE(6,'(T2,A14,T18,I6,T26,4F10.2)')&
       IF (paral%io_parent)WRITE(6,'(T2,A30,T34,I6,T42,4F10.2)')&            
            ADJUSTL(tname%tnam(i)),ncall(i),&
            avrg_excl_time(i),max_excl_time(i),&
            avrg_incl_time(i),max_incl_time(i)
       totalct=totalct+ctsum(i)
       totalwt=totalwt+wtsum(i)*1.e-3_real_8
       max_incl_time(i)=-1.0_real_8
    ENDDO
100 CONTINUE
    IF (paral%io_parent)WRITE(6,'(1X,80("*"),/)')
    DEALLOCATE(avrg_excl_time,avrg_incl_time,excl_time,incl_time,&
         max_excl_time,max_incl_time,stat=ierr)
    IF (ierr.NE.0)STOP "Deallocation error"

  END SUBROUTINE tipri

  ! ==================================================================
  SUBROUTINE tistart(t_cputime, t_walltime)
    ! ==--------------------------------------------------------------==
    ! == CALLED ONCE AT THE BEGINNING IN ORDER TO START THE TIMER     ==
    ! == AND DO SOME INITIALISATION                                   ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: t_cputime, t_walltime

    INTEGER                                  :: i
    INTEGER, SAVE                            :: ifirst = 0

    IF (ifirst.EQ.0) THEN

       ! ==--------------------------------------------------------------==
       ! == INITIALIZATION OF ARRAYS FOR SUBROUTINE TIMING               ==
       ! ==--------------------------------------------------------------==
       ntimer=0
       DO i=1,ntimx
          tname%tnam(i) = ' '
       ENDDO
       CALL zeroing(ctsum)!,ntimx)
       CALL zeroing(ctsta)!,ntimx)
       CALL zeroing(wtsum)!,ntimx)
       CALL zeroing(wtsta)!,ntimx)
       CALL zeroing(ctsum)!,ntimx)
       CALL zeroing(ncall)!,ntimx)
       CALL zeroing(incl_wtime)!,ntimx)

       ! >>>> init trace
       tname%trace_depth = 0
       tname%trace_names=REPEAT('X',LEN(tname%trace_names))
       ! <<<<

       t_walltime=m_walltime()
       t_cputime=m_cputime()
       ifirst=1
       tjsloop=t_walltime
       tjstart=t_walltime
    ELSE
       WRITE(*,*) "CALLED ONCE"
       CALL m_flush(6)
       STOP
    ENDIF
  END SUBROUTINE tistart

  ! ==================================================================
  SUBROUTINE tilimit(tstop)
    ! ==--------------------------------------------------------------==
    ! == .TRUE. IF THE JOB TIME IS TOO SHORT FOR ANOTHER LOOP         ==
    ! == A LOOP IS DEFINED BY 2 CALLS TO TILIMIT (TESTEX)             ==
    ! ==--------------------------------------------------------------==
    LOGICAL                                  :: tstop

    REAL(real_8)                             :: tjeloop, tremain

    IF (timestop) THEN
       tstop=.TRUE.
       RETURN
    ENDIF
    tstop=.FALSE.
    timestop=.FALSE.
    IF (tjlimit.EQ.0) THEN
       RETURN
    ENDIF
    ! End of time for the last loop
    tjeloop=m_walltime()

    ! Total time for the last loop
    timeloop=(tjeloop-tjsloop)*1.e-3_real_8

    ! Remaining job time 
    tremain=tjlimit-(tjeloop-tjstart)*1.e-3_real_8
    ! We add 10%.
    IF (tremain.LT.(1.1_real_8*timeloop)) THEN
       tstop=.TRUE.
       timestop=.TRUE.
    ENDIF
    tjsloop=tjeloop
  END SUBROUTINE tilimit

  ! ==================================================================
  SUBROUTINE tilimex
    ! ==--------------------------------------------------------------==
    ! == EXIT IF THE JOB TIME LIMIT IS REACHED                        ==
    ! ==--------------------------------------------------------------==
    ! Variables
    CHARACTER(len=26)                        :: datx
    REAL(real_8)                             :: time1, tremain

    IF (.NOT.timestop) THEN
       RETURN
    ENDIF
    ! Remaining job time
    time1=m_walltime()
    tremain=tjlimit-(time1-tjstart)*1.e-3_real_8
    CALL m_datum(datx)
    IF (paral%io_parent)&
         WRITE(6,'(1X,64("*"))')
    IF (paral%io_parent)&
         WRITE(6,'(1X,"*",62X,"*")')
    IF (paral%io_parent)&
         WRITE(6,'(1X,"*",12X,A,T65,"*")')&
         'JOB LIMIT TIME EXCEEDED FOR A NEW LOOP'
    IF (paral%io_parent)&
         WRITE(6,'(1X,"*",62X,"*")')
    IF (paral%io_parent)&
         WRITE(6,'(1X,64("*"))')
    IF (paral%io_parent)&
         WRITE(6,'(1X,"*",A,T43,F13.3," SECONDS *")')&
         ' THE JOB TIME LIMIT IS:',tjlimit
    IF (paral%io_parent)&
         WRITE(6,'(1X,"*",A,T43,F13.3," SECONDS *")')&
         ' THE TIME OF THE LAST LOOP IS:',timeloop
    IF (paral%io_parent)&
         WRITE(6,'(1X,"*",A,T43,F13.3," SECONDS *")')&
         ' THE REMAINING JOB TIME IS:',tremain
    IF (paral%io_parent)&
         WRITE(6,'(A,A,A,A)')&
         ' *     ',' THE COMMAND WAS ISSUED AT ',datx(1:24),'      *'
    IF (paral%io_parent)&
         WRITE(6,'(1X,"*",62X,"*")')
    IF (paral%io_parent)&
         WRITE(6,'(1X,64("*"))')
  END SUBROUTINE tilimex

  ! ==================================================================
  SUBROUTINE tanalyze(label,routine)
    CHARACTER(len=*)                         :: label, routine

    INTEGER                                  :: i, nc, nr
    INTEGER, SAVE                            :: ifirst = 0, nclast(ntimx)
    REAL(real_8)                             :: ct, wt
    REAL(real_8), SAVE                       :: ctlast(ntimx), wtlast(ntimx)

    IF (ifirst.EQ.0) THEN
       DO i=1,ntimx
          ctlast(i)=0._real_8
          wtlast(i)=0._real_8
          nclast(i)=0
       ENDDO
       ifirst=1
    ENDIF
    DO i=1,ntimer
       IF (INDEX(tname%tnam(i),routine).NE.0) THEN
          nr=i
          GOTO 100
       ENDIF
    ENDDO
    IF (paral%io_parent)&
         WRITE(6,'(A6,A,A,A)') label," >>>>>>>>>>",routine," NOT FOUND"
    RETURN
100 CONTINUE
    ct=ctsum(nr)-ctlast(nr)
    wt=wtsum(nr)-wtlast(nr)
    nc=ncall(nr)-nclast(nr)
    IF (paral%io_parent)&
         WRITE(6,'(A6,A,A10,A,I8,A,F10.2,A,F10.2)') label,&
         " >>>>>>>>>>",routine,"   CALLS=",nc,"  CPU-TIME=",ct,&
         "  WALL-TIME=",wt
    ctlast(nr)=ctsum(nr)
    wtlast(nr)=wtsum(nr)
    nclast(nr)=ncall(nr)
  END SUBROUTINE tanalyze

  ! ==================================================================
  SUBROUTINE ttimp(tcpu,twck)
    REAL(real_8)                             :: tcpu, twck

    CHARACTER(len=40)                        :: string1, string2
    INTEGER                                  :: nth1, nth2, ntm1, ntm2
    REAL(real_8)                             :: t1, t2, ts1, ts2

    nth1 = INT(tcpu/3600._real_8)
    nth2 = INT(twck/3600._real_8)
    t1   = tcpu-REAL(nth1,kind=real_8)*3600._real_8
    t2   = twck-REAL(nth2,kind=real_8)*3600._real_8
    ntm1 = INT(t1/60._real_8)
    ntm2 = INT(t2/60._real_8)
    ts1  = t1-REAL(ntm1,kind=real_8)*60._real_8
    ts2  = t2-REAL(ntm2,kind=real_8)*60._real_8
    IF (paral%io_parent)&
         WRITE(string1,'(I4,A,I2,A,F5.2,A)') nth1,' HOURS ',ntm1,&
         ' MINUTES ',ts1,' SECONDS'
    IF (paral%io_parent)&
         WRITE(string2,'(I4,A,I2,A,F5.2,A)') nth2,' HOURS ',ntm2,&
         ' MINUTES ',ts2,' SECONDS'
    IF (paral%io_parent)&
         WRITE(6,'(A,A)') '       CPU TIME : ',string1
    IF (paral%io_parent)&
         WRITE(6,'(A,A)') '   ELAPSED TIME : ',string2
  END SUBROUTINE ttimp
  ! ==================================================================

END MODULE timer


!vw there is circular dependency between the error_handling and timer modules.
!vw so we keep this guy outside.
SUBROUTINE tistopgm(file_unit)
  ! ==--------------------------------------------------------------==
  ! == DISPLAY THE STACK OF CALLS IF STOPGM IS CALLED               ==
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE machine, ONLY: m_flush
  USE time , ONLY:tname
  IMPLICIT NONE
  INTEGER                                    :: file_unit

  INTEGER                                    :: i

  IF (tname%trace_depth.GT.0) THEN
     WRITE(file_unit,'(A)') ' call stack:'
     DO i=tname%trace_depth,1,-1
        WRITE(file_unit,'(10X,I3,2X,A)') i,TRIM(tname%trace_names(i))
     ENDDO
     CALL m_flush(file_unit)
  ENDIF
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE tistopgm
! ==================================================================
