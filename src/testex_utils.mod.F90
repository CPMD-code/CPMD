MODULE testex_utils
  USE error_handling,                  ONLY: stopgm
  USE fileopenmod,                     ONLY: fo_stat_max
  USE kinds,                           ONLY: default_path_length
  USE mp_interface,                    ONLY: mp_sum
  USE mw,                              ONLY: tmw
  USE parac,                           ONLY: parai,&
                                             paral
  USE pimd,                            ONLY: grandparent,&
                                             pc_groups,&
                                             supergroup
  USE softex_utils,                    ONLY: softex
  USE system,                          ONLY: cntl
  USE timer,                           ONLY: tilimex,&
                                             tilimit

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: testex
  PUBLIC :: testex_mw

CONTAINS

  ! ==================================================================
  SUBROUTINE testex(sexit)
    ! ==--------------------------------------------------------------==
    ! ==  TEST FOR FILE EXIT                                          ==
    ! ==  TEST FOR A SOFT EXIT CALL OVER N PROCESSORS                 ==
    ! ==  TEST IF THE JOB TIME LIMIT IS ENOUGH                        ==
    ! ==  CHECK MEMORY                                                ==
    ! ==--------------------------------------------------------------==
    LOGICAL                                  :: sexit

    CHARACTER(*), PARAMETER                  :: procedureN = 'testex'

    CHARACTER(LEN=default_path_length)       :: filen
    INTEGER                                  :: in, mpigroup, out, unit
    LOGICAL                                  :: exist, opened, t_io_node, &
                                                test1
    LOGICAL, SAVE                            :: test2 = .FALSE.

#ifdef __ES 
    CHARACTER(len=130) :: exit_fpath ! PATH TO EXIT FILE ON ES
    COMMON/fpath/ exit_fpath
#endif 
    ! ==--------------------------------------------------------------==
    ! AK: 2007-07-02 
    ! with processor groups we cannot sync across the SUPERGROUP 
    ! as this will create a deadlock when using BO-PIMD or MW MTD.
    IF (cntl%tpath) THEN
       IF (pc_groups.GT.1) THEN
          !  t_io_node = paral%parent
          !  mpigroup = parai%allgrp
          t_io_node = paral%io_parent
          mpigroup = parai%cp_grp
       ELSE
          t_io_node = grandparent
          mpigroup = supergroup
       ENDIF
    ELSE
       ! Also multiple walkers should check only internally and not on the
       ! supergroup! TL
       t_io_node = paral%io_parent
       mpigroup = parai%cp_grp
    ENDIF
    IF (t_io_node) THEN
       filen='EXIT'
#ifdef __ES
       filen=TRIM(exit_fpath)//TRIM(filen)
#endif
       INQUIRE(file=TRIM(filen),exist=test1)
       IF (test1.OR.test2) THEN
          sexit=.TRUE.
          CALL softex(0)
       ENDIF
       ! AK: 2007-07-02
       ! do not remove the EXIT file here for PIMD and MW MTD, so that _all_
       ! processor groups can see it before the first one deletes it.
       ! we set PC_GROUPS to 1 at the end of PI_MDPT and enter here again.
       ! for MW MTD the EXIT file is removed in TESTEX_MW instead. 
       IF ((.NOT.(cntl%tpath.OR.tmw)).OR.(cntl%tpath.AND.(pc_groups.EQ.1))) THEN
          IF (test1) THEN
             DO unit=0,fo_stat_max
                INQUIRE(unit=unit,opened=opened,exist=exist)
                IF (.NOT.opened.AND.exist) EXIT
             ENDDO
             IF (unit>fo_stat_max) CALL stopgm(procedureN,&
                  'cannot find a valid unit',& 
                  __LINE__,__FILE__)
             OPEN(unit=unit,file='EXIT',status='OLD')
             CLOSE(unit=unit,status='DELETE')
          ENDIF
       ENDIF
       ! Check job limit time
       CALL tilimit(test1)
       IF (test1) THEN
          sexit=.TRUE.
          CALL tilimex
       ENDIF
    ENDIF
    in=0
    out=0
    IF (sexit) THEN
       in=1
    ENDIF
    CALL mp_sum(in,out,mpigroup)
    sexit = .FALSE.
    IF (out.NE.0) THEN
       sexit=.TRUE.
    ENDIF
    IF (sexit) THEN
       test2=.TRUE.
    ENDIF
  END SUBROUTINE testex
  ! ==================================================================
  SUBROUTINE testex_mw(sexit)
    ! ==--------------------------------------------------------------==
    ! == Reset SEXIT=.TRUE. if any distributed SEXIT is set to .TRUE. == 
    ! ==--------------------------------------------------------------==
    LOGICAL                                  :: sexit

    CHARACTER(*), PARAMETER                  :: procedureN = 'testex_mw'

    CHARACTER(LEN=default_path_length)       :: filen
    INTEGER                                  :: in, out, unit
    LOGICAL                                  :: exist, opened, test1

#ifdef __ES
    CHARACTER(len=130) :: exit_fpath ! PATH TO EXIT FILE ON ES
    COMMON/fpath/ exit_fpath
#endif
    ! ==--------------------------------------------------------------==
    ! Only one proc should ever arrive at this point
    IF (grandparent) THEN
       filen='EXIT'
#ifdef __ES
       filen=TRIM(exit_fpath)//TRIM(filen)
#endif
       INQUIRE(file=TRIM(filen),exist=test1)
       IF (test1) THEN
          DO unit=0, fo_stat_max
             INQUIRE(unit=unit,opened=opened,exist=exist)
             IF (.NOT.opened.AND.exist) EXIT
          ENDDO
          IF (unit>fo_stat_max) CALL stopgm(procedureN,&
               'cannot find a valid unit',& 
               __LINE__,__FILE__)
          OPEN(unit=unit,file='EXIT',status='OLD')
          CLOSE(unit=unit,status='DELETE')
       ENDIF
    ENDIF
    in=0
    out=0
    IF (sexit) in=1
    CALL mp_sum(in,out,supergroup)
    sexit = .FALSE.
    IF (out.NE.0) sexit=.TRUE.
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE testex_mw
  ! ==================================================================

END MODULE testex_utils
