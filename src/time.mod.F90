MODULE time
  USE kinds,                           ONLY: real_8
! ==================================================================
! == INCLUDE FILE FOR TIMING USED BY TISET AND TIHALT ROUTINES    ==
! ==================================================================
! == NTIMX  Maximum number of timers                              ==
! == LENTIM Maximum length of names                               ==
! ==--------------------------------------------------------------==

  IMPLICIT NONE
  INTEGER, PARAMETER :: ntimx=200 
  INTEGER, PARAMETER :: lentim=63
  ! ==--------------------------------------------------------------==
  ! == TNAM(1:NTIMX) Name of routines                               ==
  ! ==--------------------------------------------------------------==



  TYPE :: tname_t
     CHARACTER(len=lentim) :: tnam(ntimx) = ''
     CHARACTER(len=lentim) :: trace_names(ntimx) = ''
     INTEGER :: trace_depth = HUGE(0)
     INTEGER :: trace_nbr_procedure = HUGE(0)
     CHARACTER(len=lentim) :: trace_procedure(ntimx) = ''
     INTEGER :: trace_max_depth = HUGE(0)
     INTEGER :: trace_max_calls = HUGE(0)
  END TYPE tname_t
  TYPE(tname_t), SAVE :: tname
  ! ==================================================================
  ! == CTSUM(1:NTIMX) Sum of CPU Time for each routine              ==
  ! == CTSTA(1:NTIMX) CPU Time Start                                ==
  ! == WTSUM(1:NTIMX) Sum of Elapsed Time                           ==
  ! == WTSTA(1:NTIMX) Elapsed Time Start                            ==
  ! ==--------------------------------------------------------------==
  REAL(real_8), SAVE :: ctsum(ntimx)=HUGE(0.0_real_8),ctsta(ntimx)=HUGE(0.0_real_8),&
       wtsum(ntimx)=HUGE(0.0_real_8),wtsta(ntimx)=HUGE(0.0_real_8),&
       time_in(ntimx,2)=HUGE(0.0_real_8),incl_wtime(ntimx)=HUGE(0.0_real_8)
  TYPE :: timbu_t
     REAL(real_8) :: ctsum(ntimx)=HUGE(0.0_real_8)
     REAL(real_8) :: ctsta(ntimx)=HUGE(0.0_real_8)
     REAL(real_8) :: wtsum(ntimx)=HUGE(0.0_real_8)
     REAL(real_8) :: wtsta(ntimx)=HUGE(0.0_real_8)
     REAL(real_8) :: time_in(ntimx,2)=HUGE(0.0_real_8)
     REAL(real_8) :: incl_wtime(ntimx)=HUGE(0.0_real_8)
  END TYPE timbu_t
  TYPE(timbu_t), SAVE :: timbu
  ! ==================================================================
  ! == NTIMER         Number of timers                              ==
  ! == NCALL(1:NTIMX) Number of calls for each routine              ==
  ! ==--------------------------------------------------------------==
  INTEGER, SAVE :: ntimer,ncall(ntimx)=HUGE(0)
  TYPE :: timin_t
     INTEGER :: ntimer=HUGE(0)
     INTEGER :: ncall(ntimx)=HUGE(0)
  END TYPE timin_t
  TYPE(timin_t), SAVE :: timin
  ! ==================================================================
  ! == TJSTART  Time start for CPMD                                 ==
  ! == TJSLOOP  Time of the last call for TILIMIT                   ==
  ! ==          Give time for each loop                             ==
  ! == TIMELOOP CPU time of the last loop i.e.                      ==
  ! ==          between 2 calls of TILIMIT                          ==
  ! == TIMESTOP .TRUE. if the Job time limit is almost exceeded     ==
  ! ==--------------------------------------------------------------==
  REAL(real_8), SAVE :: tjstart=HUGE(0.0_real_8),tjsloop=HUGE(0.0_real_8),&
       timeloop=HUGE(0.0_real_8)
  LOGICAL, SAVE :: timestop=.FALSE.
  TYPE :: timli_t
     REAL(real_8) :: tjstart=HUGE(0.0_real_8)
     REAL(real_8) :: tjsloop=HUGE(0.0_real_8)
     REAL(real_8) :: timeloop=HUGE(0.0_real_8)
     LOGICAL :: timestop=.FALSE.
  END TYPE timli_t
  TYPE(timli_t), SAVE :: timli
  ! ==================================================================
END MODULE time
