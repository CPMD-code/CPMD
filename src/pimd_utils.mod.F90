MODULE pimd_utils

  USE error_handling,                  ONLY: stopgm
  USE pimd,                            ONLY: pc_grp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: alloc_pimd

CONTAINS

  SUBROUTINE alloc_pimd( max_nproc )
    INTEGER, INTENT(IN)                      :: max_nproc

    CHARACTER(*), PARAMETER                  :: procedureN = 'alloc_pimd'

    INTEGER                                  :: ierr

    ALLOCATE( pc_grp( max_nproc ), &
         &    STAT=ierr )
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    pc_grp = HUGE(0)

  END SUBROUTINE alloc_pimd

END MODULE pimd_utils
