MODULE shock
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==================================================================
  ! == SHOCK WAVE PARAMETERS                                        ==
  ! ==================================================================

  ! ==================================================================
  TYPE :: shock1_t
     REAL(real_8) :: vshock
     REAL(real_8) :: pshock
     REAL(real_8) :: vol0
     REAL(real_8) :: eshock
  END TYPE shock1_t
  TYPE(shock1_t) :: shock1

END MODULE shock
