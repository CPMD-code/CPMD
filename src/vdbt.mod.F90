MODULE vdbt
  USE system,                          ONLY: maxsp

  IMPLICIT NONE

  ! ==================================================================
  ! == VANDERBILT PSEUDOPOTENTIAL                                   ==
  ! ==================================================================
  INTEGER :: itmax(maxsp)
  CHARACTER (len=66) :: vdbti(60,maxsp)
  ! ==================================================================

END MODULE vdbt
