MODULE pm_cntl_utils
  USE cotr,                            ONLY: cotc0,&
                                             cotr007,&
                                             lskcor
  USE error_handling,                  ONLY: stopgm
  USE inscan_utils,                    ONLY: inscan
  USE ions,                            ONLY: ions1
  USE kinds,                           ONLY: real_8
  USE mfep,                            ONLY: alpha_pmin,&
                                             factor,&
                                             irest,&
                                             mfepi,&
                                             mfepl,&
                                             tolds
  USE parac,                           ONLY: parai,&
                                             paral
  USE pimd,                            ONLY: maxnp,&
                                             pc_groups,&
                                             pimd1,&
                                             pimd2,&
                                             pimd3,&
                                             repfname
  USE readsr_utils,                    ONLY: xstring
  USE store_types,                     ONLY: restart1
  USE system,                          ONLY: cnti,&
                                             cntl
#include "sizeof.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: pm_cntl

CONTAINS

  ! ==================================================================
  SUBROUTINE pm_cntl
    ! ==--------------------------------------------------------------==
    ! ==  THIS ROUTINE READS THE SECTION &PATH &END ON INPUT UNIT     ==
    ! ==--------------------------------------------------------------==
    ! ==  THE INPUT HAS TO BE IN THE FOLLOWING FORMAT                 ==
    ! ==     &PATH                                                    ==
    ! ==        Options                                               ==
    ! ==     &END                                                     ==
    ! ==--------------------------------------------------------------==
    ! ==  LIST OF OPTIONS                                             ==
    ! ==                                                              ==
    ! ==    REPLICA NUMBER                                            ==
    ! ==      np_total                                                ==
    ! ==    CLASSICAL TEST                                            == 
    ! ==    INITIALIZATION                                            ==
    ! ==    READ REPLICAS                                             ==
    ! ==       filename                                               ==
    ! ==    GENERATE REPLICAS                                         ==
    ! ==    PROCESSOR GROUPS                                          ==
    ! ==       pc_groups                                              ==
    ! ==    OUTPUT [GROUPS,PARENT,ALL]                                ==
    ! ==    PRINT LEVEL                                               ==
    ! ==       levprt                                                 ==
    ! ==    NLOOP                                                     ==
    ! ==      nloop                                                   ==
    ! ==    FACTOR                                                    ==
    ! ==      factor                                                  ==
    ! ==    ALPHA                                                     ==
    ! ==      alpha_pmin                                              ==
    ! ==    NEQUI                                                     ==
    ! ==      nequi                                                   ==
    ! ==    NPREVIOUS [LATEST]                                        ==
    ! ==       nloopold                                               ==
    ! ==    PROJOUT                                                   ==
    ! ==    CONVERGENCE                                               ==
    ! ==       tolds                                                  ==
    ! ==    COLLECTIVE VARIABLES                                      ==
    ! ==       ncvsp                                                  ==
    ! ==       irest(1) irest(2) ... irest(ncvsp)                     ==
    ! ==--------------------------------------------------------------==
    ! Variables
    CHARACTER(*), PARAMETER                  :: procedureN = 'pm_cntl'

    CHARACTER(len=80)                        :: line
    INTEGER                                  :: i, i1, i2, ia, ierr, ire, &
                                                iunit, maxir, minir, nfix
    LOGICAL                                  :: tcvsp

! ==--------------------------------------------------------------==
! ==  DEFAULTS                                                    ==
! ==--------------------------------------------------------------==

    IF (.NOT.paral%io_parent) GOTO 9999
    iunit = 5
    pimd1%rrep=.FALSE.
    pimd1%repread=.FALSE.
    pimd1%tread_cent = .TRUE.
    pimd1%tstage=.FALSE.
    pimd1%tpinm=.FALSE.
    pimd1%tinit=.FALSE.
    pimd1%tcentro=.FALSE.
    pimd1%testpi=.FALSE.
    pimd3%np_total = 1
    pimd2%tempb = 500._real_8
    pimd2%wmass = 1._real_8
    pimd2%facstage=1._real_8
    pc_groups=0
    repfname=' '
    pimd3%levprt=0
    pimd3%loutfn=1
    mfepi%nloop=1
    mfepi%nloopold=0
    mfepi%nequi=0
    factor=1.0_real_8
    alpha_pmin=0.2_real_8
    mfepl%projout=.FALSE.
    mfepl%rlate_string=.FALSE.
    tolds=0.0_real_8
    tcvsp=.FALSE.
    mfepi%ncvsp=cotr007%mrestr
    ! ==--------------------------------------------------------------==
    ierr=inscan(iunit,'&PATH')
    IF (ierr.NE.0) GOTO 100
    ! ==--------------------------------------------------------------==
10  CONTINUE
    IF (paral%io_parent)&
         READ(iunit,err=99,END=99,fmt='(A80)') line
    IF (INDEX(line,'&END').NE.0) THEN
       GOTO 100
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (INDEX(line,'REPL').NE.0 .AND. INDEX(line,'NUM').NE.0) THEN
       ! Number of replica NP_TOTAL
       IF (paral%io_parent)&
            READ(5,err=99,END=99,fmt=*) pimd3%np_total
       GOTO 10
    ELSE IF&
         (INDEX(line,'CLASS').NE.0.AND.INDEX(line,'TEST').NE.0) THEN
       pimd1%testpi=.TRUE.
       GOTO 10
    ELSE IF&
         (INDEX(line,'INITIALIZA').NE.0) THEN
       ! Initialization run: supply full start configuration either by
       ! (i)  generation -> with RREP=.TRUE.    or by
       ! (ii) reading it -> with REPREAD=.TRUE.
       ! is default if no restart of coordinates
       pimd1%tinit=.TRUE.
       GOTO 10
    ELSE IF&
         (INDEX(line,'GENERATE').NE.0) THEN
       ! Generate replica coordinates
       pimd1%rrep=.TRUE.
       GOTO 10
    ELSE IF&
         (INDEX(line,'READ').NE.0.AND.INDEX(line,'REPL').NE.0) THEN
       ! Read replica coordinates from external file 
       ! to be named 'repfname' in the loop format 
       ! "DO IP, DO IS, DO IA ... (K=1,3)" (see 'rread.F')
       ! with an empty line before EVERY set of replicas (also the first!)
       pimd1%repread=.TRUE.
       IF (paral%io_parent)&
            READ(iunit,err=99,END=99,fmt='(A80)') repfname
       GOTO 10
    ELSE IF&
         (INDEX(line,'PROC').NE.0.AND.INDEX(line,'GROUP').NE.0) THEN
       ! Number of processor goups 
       IF (paral%io_parent)&
            READ(5,err=99,END=99,fmt=*) pc_groups
       GOTO 10
    ELSE IF&
         (INDEX(line,'PRIN').NE.0.AND.INDEX(line,'LEVEL').NE.0) THEN
       ! Printing level: minimal output for LEVPRT < 5 
       IF (paral%io_parent)&
            READ(5,err=99,END=99,fmt=*) pimd3%levprt
       GOTO 10
    ELSE IF&
         (INDEX(line,'OUTPUT').NE.0) THEN
       ! Output files for each processor, group or only grandparent
       IF (INDEX(line,'PARENT').NE.0) pimd3%loutfn=0
       IF (INDEX(line,'GROUP').NE.0) pimd3%loutfn=1
       IF (INDEX(line,'ALL').NE.0) pimd3%loutfn=2
       GOTO 10
    ELSE IF&
         (INDEX(line,'NLOOP').NE.0) THEN
       ! Number of loop for minimum free energy path
       READ(5,err=99,END=99,fmt=*) mfepi%nloop
       GOTO 10
    ELSE IF&
         (INDEX(line,'FACTOR').NE.0) THEN
       ! Factor for path propagation
       READ(5,err=99,END=99,fmt=*) factor
       GOTO 10
    ELSE IF&
         (INDEX(line,'ALPHA').NE.0) THEN
       ! Factor for path smoothing
       READ(5,err=99,END=99,fmt=*) alpha_pmin
       GOTO 10
    ELSE IF&
         (INDEX(line,'NEQUI').NE.0) THEN
       ! Number of equilibration steps to be skipped 
       ! for minimum free energy path
       READ(5,err=99,END=99,fmt=*) mfepi%nequi
       GOTO 10
    ELSE IF&
         (INDEX(line,'NPREVIOUS').NE.0) THEN
       ! Loop number of last run read from LATEST_STRING file
       ! or specified below
       IF (INDEX(line,'LATEST').NE.0) THEN
          mfepl%rlate_string=.TRUE.
       ELSE
          READ(5,err=99,END=99,fmt=*) mfepi%nloopold
       ENDIF
       GOTO 10
    ELSE IF&
         (INDEX(line,'PROJOUT').NE.0) THEN
       ! Forces projected out
       mfepl%projout=.TRUE.
       GOTO 10
    ELSE IF&
         (INDEX(line,'CONVERG').NE.0) THEN
       ! Convergence criterion on string positions
       READ(5,err=99,END=99,fmt=*) tolds
       GOTO 10
    ELSE IF(INDEX(line,'COLLE').NE.0.AND.&
         INDEX(line,'VAR').NE.0) THEN
       ! Confined collective variables space
       READ(5,err=99,END=99,fmt=*) mfepi%ncvsp
       ALLOCATE(irest(mfepi%ncvsp),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       READ(5,err=99,END=99,fmt=*) (irest(i),i=1,mfepi%ncvsp)
       tcvsp=.TRUE.
       GOTO 10
    ELSE
       ! Dummy line
       GOTO 10
    ENDIF
    ! ==--------------------------------------------------------------==
99  CONTINUE
    IF (paral%io_parent)&
         WRITE(6,*) ' PM_CNTL: ERROR IN READING INPUT FILE'
    CALL stopgm('PM_CNTL',' ',& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
100 CONTINUE
    ! ==--------------------------------------------------------------==
    ! ==  Test of options                                             ==
    ! ==--------------------------------------------------------------==
    IF (pimd3%np_total.GT.maxnp) THEN
       IF (paral%io_parent)&
            WRITE(6,'(A,I4,A,A,I4)')&
            ' SELECTED REPLICA NUMBER ',pimd3%np_total,' LARGER THAN',&
            ' CURRENT MAXIMUM VALUE MAXNP = ',maxnp
       CALL stopgm('PM_CNTL',' REPLICA NUMBER LARGER THAN MAXNP ',& 
            __LINE__,__FILE__)
    ENDIF
    IF (cntl%geopt.OR.cntl%md) THEN
    ELSEIF (cntl%wfopt) THEN
       pc_groups = 1
       IF (.NOT.restart1%restart) THEN
          pimd1%tinit=.TRUE.
       ENDIF
       IF (pimd1%tinit) THEN
          IF (.NOT.(pimd1%rrep.OR.pimd1%repread)) pimd1%rrep=.TRUE.
       ELSE
          restart1%rco=.TRUE.
       ENDIF
    ELSE
       IF (paral%io_parent)&
            WRITE(6,*) ' ONLY MD, AND GEOMETRY OR WAVEFUNCTION ',&
            'OPTIMIZATION IMPLEMENTED'
       CALL stopgm('PM_CNTL',' WRONG OPTION ',& 
            __LINE__,__FILE__)
    ENDIF
    nfix=0
    DO ia=1,ions1%nat
       IF (lskcor(1,ia).EQ.0) THEN
          nfix=nfix+1
       ENDIF
       IF (lskcor(2,ia).EQ.0) THEN
          nfix=nfix+1
       ENDIF
       IF (lskcor(3,ia).EQ.0) THEN
          nfix=nfix+1
       ENDIF
    ENDDO
    nfix=nfix+cotc0%mcnstr
    IF (cntl%md.AND.(cnti%nomore-mfepi%nequi*cnti%ntraj).LE.0) THEN
       WRITE(6,'(A,I6,A)') ' REQUESTED MAXSTEP OF ',cnti%nomore,&
            ' IS TOO SMALL'
       CALL stopgm('PM_CNTL',' WRONG OPTION ',& 
            __LINE__,__FILE__)
    ENDIF
    IF (tcvsp) THEN
       IF (mfepi%ncvsp.LT.1.OR.mfepi%ncvsp.GT.cotr007%mrestr)&
            CALL stopgm('PM_CNTL',' WRONG NUMBER OF NCVSP ',& 
            __LINE__,__FILE__)
    ELSE
       ALLOCATE(irest(mfepi%ncvsp),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       DO i=1,mfepi%ncvsp
          irest(i)=i
       ENDDO
    ENDIF
    minir=1
    maxir=cotr007%mrestr
    DO i=1,mfepi%ncvsp
       ire=irest(i)
       minir=MIN(ire,minir)
       maxir=MAX(ire,maxir)
    ENDDO
    IF (minir.LT.1.OR.maxir.GT.cotr007%mrestr)&
         CALL stopgm('PM_CNTL',' WRONG NUMBER OF IREST ',& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! ==  WRITE INFO TO OUTPUT                                        ==
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent)&
         WRITE(6,*)
    IF (paral%io_parent)&
         WRITE(6,'(1X,19(">"),A,15("<"))') ' PATH MINIMISATION PARAMETERS '
    IF (pimd1%testpi) THEN
       pimd3%np_total=1
       IF (paral%io_parent)&
            WRITE(6,'(1X,4("! "),A)')&
            ' THIS IS A TEST WITH NP=1 FOR COMPARISON WITH CLASSICAL RUNS'
       IF (paral%io_parent)&
            WRITE(6,'(1X,4("! "),A)')&
            ' PRIMITIVE PROPAGATOR AND LEVY INITIAL GUESS USED'
    ENDIF
    IF (paral%io_parent)&
         WRITE(6,'(1X,4(">"),A,I4)')&
         ' NUMBER OF REPLICA:',pimd3%np_total
    IF (pimd1%tinit) THEN
       IF (paral%io_parent)&
            WRITE(6,'(1X,4(">"),A)')&
            ' INITIALIZATION RUN'
    ENDIF
    IF (pimd1%rrep) THEN
       IF (paral%io_parent)&
            WRITE(6,'(1X,4(">"),A)')&
            ' GENERATE REPLICA COORDINATES '
    ENDIF
    CALL xstring(repfname,i1,i2)
    IF (pimd1%repread) THEN
       IF (paral%io_parent)&
            WRITE(6,'(1X,4(">"),A)')&
            ' READ REPLICA COORDINATES FROM FILE:'
       IF (paral%io_parent)&
            WRITE(6,'(1X,4(">"),A,A,A)') '"',repfname(i1:i2),'"'
    ENDIF
    WRITE(6,'(1X,4(">"),A,I4)')&
         ' PRINT LEVEL:',pimd3%levprt
    IF (tolds.GT.0.0_real_8)&
         WRITE(6,'(1X,4(">"),A,1PE12.6)')&
         ' CONVERGENCE CRITERIA FOR STRING POSITIONS: ',tolds
    IF (mfepi%ncvsp.LT.cotr007%mrestr) THEN
       WRITE(6,'(1X,4(">"),A,I3)')&
            ' DIMENSION OF CONFINED CV SPACE: ',mfepi%ncvsp
       WRITE(6,'(1X,4(">"),A)')&
            '   RESTRAINT NUMBER: '
       WRITE(6,'(1X,4(">"),T10,17I3)') (irest(i),i=1,mfepi%ncvsp)
    ENDIF
    WRITE(6,'(1X,4(">"))')
    ! ==--------------------------------------------------------------==
9999 CONTINUE
    ! ==--------------------------------------------------------------==
    ! ==  BROADCAST INFO TO ALL OTHER PROCESSORS                      ==
    ! ==--------------------------------------------------------------==
    CALL mp_bcast_byte(pimd1, size_in_bytes_of(pimd1),parai%io_source,parai%cp_grp)
    CALL mp_bcast_byte(pimd2, size_in_bytes_of(pimd2),parai%io_source,parai%cp_grp)
    CALL mp_bcast_byte(pimd3, size_in_bytes_of(pimd3),parai%io_source,parai%cp_grp)
    CALL mp_bcast_byte(mfepi, size_in_bytes_of(mfepi),parai%io_source,parai%cp_grp)
    CALL mp_bcast_byte(mfepl, size_in_bytes_of(mfepl),parai%io_source,parai%cp_grp)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE pm_cntl
  ! ==================================================================

END MODULE pm_cntl_utils
