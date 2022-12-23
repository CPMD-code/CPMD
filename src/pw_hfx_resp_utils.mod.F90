MODULE pw_hfx_resp_utils
  USE error_handling,                  ONLY: stopgm
  USE pw_hfx_resp_types,               ONLY: hfx_resp_env,&
                                             pw_hfx_resp_t

  IMPLICIT NONE

  PRIVATE 

  PUBLIC :: pw_hfx_resp_create, &
       pw_hfx_resp_destroy

CONTAINS

  SUBROUTINE pw_hfx_resp_create(this, nnr1, nstate)
    TYPE(pw_hfx_resp_t), INTENT(inout)       :: this
    INTEGER, INTENT(in)                      :: nnr1, nstate

    CHARACTER(*), PARAMETER :: procedureN = 'pw_hfx_resp_create'

    INTEGER                                  :: ierr

    IF( this%init ) CALL stopgm(procedureN,'already initialized',&
         __LINE__,__FILE__)
    this%init=.TRUE.

    IF (.NOT.hfx_resp_env%write_v_ab) THEN
       ALLOCATE(this%v_ab(nnr1,nstate,nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
    ENDIF
  END SUBROUTINE pw_hfx_resp_create

  SUBROUTINE pw_hfx_resp_destroy(this)
    TYPE(pw_hfx_resp_t), INTENT(inout)       :: this

    CHARACTER(*), PARAMETER :: procedureN = 'pw_hfx_resp_destroy'

    INTEGER                                  :: ierr

    IF( .NOT.this%init ) CALL stopgm(procedureN,'not initialized',&
         __LINE__,__FILE__)
    this%init=.FALSE.

    IF (.NOT.hfx_resp_env%write_v_ab) THEN 
       DEALLOCATE(this%v_ab,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    ENDIF
  END SUBROUTINE pw_hfx_resp_destroy

END MODULE pw_hfx_resp_utils

