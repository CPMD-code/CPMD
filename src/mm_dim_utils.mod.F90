MODULE mm_dim_utils
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE mm_dimmod,                       ONLY: mm_go_mm,&
                                             mm_go_qm,&
                                             mm_revert,&
                                             mm_stat,&
                                             mmdim,&
                                             nam,&
                                             naq
  USE mm_input,                        ONLY: clc,&
                                             lqmmm
  USE parac,                           ONLY: paral
  USE system,                          ONLY: cntl,&
                                             maxsys
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mm_dim
  !public :: mm_dim_query

CONTAINS

  SUBROUTINE mm_dim(ch,stat)
    ! this routine modifies NA(*) and NSP in ions.inc for QM and MM ordering
    ! available options for 'ch': mm_go_qm, mm_go_mm, mm_revert
    ! on exit 'stat' contains the previous state (= value of mm_stat)
    ! mm_stat=.true. stands for QM, mm_stat=.false. for MM dimensions.
    ! for mm_revert the value of 'stat' determines the action.
    INTEGER                                  :: ch
    LOGICAL                                  :: stat

    CHARACTER(*), PARAMETER                  :: procedureN = 'mm_dim'
    CHARACTER(len=9), DIMENSION(3), PARAMETER :: &
      mm_go_flag = (/'mm_go_qm ','mm_go_mm ','mm_revert'/)

    INTEGER                                  :: act, is, isub

    CALL tiset(procedureN,isub)

    IF (lqmmm%qmmm_verbose.AND.paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,*) 'MM_DIM| called with: ', mm_go_flag(ch)
       IF (mm_stat) THEN
          IF (paral%io_parent)&
               WRITE(6,*) 'MM_DIM| old status : QM'
       ELSE
          IF (paral%io_parent)&
               WRITE(6,*) 'MM_DIM| old status : MM'
       ENDIF
    ENDIF

    ! nothing to do unless we have a QM/MM calculation
    IF (.NOT.cntl%tqmmm)GOTO 666
    IF (.NOT.lqmmm%qmmm)GOTO 666

    ! copy the action flag to a local variable in 
    ! case ch points to a read-only location
    act=ch

    ! revert to previous state
    IF (act.EQ.mm_revert) THEN
       IF (stat) THEN
          act=mm_go_qm
       ELSE
          act=mm_go_mm
       ENDIF
    ENDIF

    ! record previous state
    stat=mm_stat

    IF (act.EQ.mm_go_qm)THEN
       mm_stat=.TRUE.
       maxsys%nsx=mmdim%NSPq
       ions1%nsp=mmdim%NSPq
       ions1%nat=0
       DO is=1,mmdim%NSPq
          ions0%na(is)=NAq(is)
          ions1%nat=ions1%nat+ions0%na(is)
       ENDDO
       IF (clc%classical) THEN
          maxsys%nsx=0
          ions1%nsp=0
          ions1%nat=0
          !$omp parallel do private(is)
          DO is=1,mmdim%NSPq
             ions0%na(is)=0
          ENDDO
       ENDIF
    ELSE IF (act.EQ.mm_go_mm) THEN
       mm_stat=.FALSE.
       maxsys%nsx=mmdim%NSPm
       ions1%nsp=mmdim%NSPm
       ions1%nat=0
       DO is=1,mmdim%NSPm
          ions0%na(is)=NAm(is)
          ions1%nat=ions1%nat+ions0%na(is)
       ENDDO
    ELSE
       IF (paral%io_parent)&
            WRITE(6,*) 'selector =', act
       CALL stopgm('MM_DIM','illegal flag in mm_dim',& 
            __LINE__,__FILE__)
    ENDIF

    IF (lqmmm%qmmm_verbose.AND.paral%parent) THEN
       IF (mm_stat) THEN
          IF (paral%io_parent)&
               WRITE(6,*) 'MM_DIM| new status : QM'
       ELSE
          IF (paral%io_parent)&
               WRITE(6,*) 'MM_DIM| new status : MM'
       ENDIF
    ENDIF

666 CONTINUE
    CALL tihalt(procedureN,isub)
    RETURN
  END SUBROUTINE mm_dim

  ! write current state to output.
  SUBROUTINE mm_dim_query(tag)

    CHARACTER(len=*)                         :: tag

    IF (paral%parent) THEN
       IF (mm_stat) THEN
          IF (paral%io_parent)&
               WRITE(6,*) 'MM_DIM| current status : QM| ', tag
       ELSE
          IF (paral%io_parent)&
               WRITE(6,*) 'MM_DIM| current status : MM| ', tag
       ENDIF
    ENDIF
    RETURN
  END SUBROUTINE mm_dim_query

END MODULE mm_dim_utils
