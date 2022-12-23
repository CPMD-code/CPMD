MODULE vdwcmod_utils

  USE error_handling,                  ONLY: stopgm
  USE vdwcmod,                         ONLY: npt12

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: alloc_vdwcmod

CONTAINS

  SUBROUTINE alloc_vdwcmod( max_nproc )
    INTEGER, INTENT(IN)                      :: max_nproc

    CHARACTER(*), PARAMETER                  :: procedureN = 'alloc_vdwcmod'

    INTEGER                                  :: ierr

    ALLOCATE( npt12(0:max_nproc,2), &
         &    STAT=ierr )
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    npt12 = HUGE(0)

  END SUBROUTINE alloc_vdwcmod

END MODULE vdwcmod_utils
