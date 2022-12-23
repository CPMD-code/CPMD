MODULE system_utils

  USE error_handling,                  ONLY: stopgm
  USE system,                          ONLY: ipept,&
                                             parap,&
                                             parap_t,&
                                             paraw,&
                                             paraw_t

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: alloc_system

CONTAINS

  SUBROUTINE alloc_system( max_nproc )
    INTEGER, INTENT(IN)                      :: max_nproc

    CHARACTER(*), PARAMETER                  :: procedureN = 'alloc_system'

    INTEGER                                  :: ierr

    ALLOCATE( ipept(2,0:max_nproc), &
         &    STAT=ierr )
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    ipept = HUGE(0)

    ! some more types allocations
    CALL alloc_parap( parap, max_nproc )
    CALL alloc_paraw( paraw, max_nproc )

  END SUBROUTINE alloc_system


  SUBROUTINE alloc_parap( parap, max_nproc )
    TYPE(parap_t), INTENT(INOUT)             :: parap
    INTEGER, INTENT(IN)                      :: max_nproc

    CHARACTER(*), PARAMETER                  :: procedureN = 'alloc_parap'

    INTEGER                                  :: ierr

    ALLOCATE( parap%nrxpl (0:max_nproc,2), &
         &    parap%nrzpl (0:max_nproc,2), &
         &    parap%sparm (9,0:max_nproc), &
         &    parap%nst12 (0:max_nproc,2), &
         &    parap%pgroup(0:max_nproc)  , &
         &    parap%nlink (0:max_nproc)  , &
         &    STAT=ierr )
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    parap%nrxpl  = HUGE(0)
    parap%nrzpl  = HUGE(0)
    parap%sparm  = HUGE(0)
    parap%nst12  = HUGE(0)
    parap%pgroup = HUGE(0)
    parap%nlink  = HUGE(0)

  END SUBROUTINE alloc_parap


  SUBROUTINE alloc_paraw( paraw, max_nproc )
    TYPE(paraw_t), INTENT(INOUT)             :: paraw
    INTEGER, INTENT(IN)                      :: max_nproc

    CHARACTER(*), PARAMETER                  :: procedureN = 'alloc_paraw'

    INTEGER                                  :: ierr

    ALLOCATE( paraw%nwa12(0:max_nproc,2), &
         &    STAT=ierr )
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    paraw%nwa12 = HUGE(0)

  END SUBROUTINE alloc_paraw

END MODULE system_utils
