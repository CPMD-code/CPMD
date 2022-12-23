MODULE benc
  IMPLICIT NONE

  ! ==================================================================
  ! ==                 BENCHMARK OPTIONS                            ==
  ! ==                                                              ==
  ! == IBENCH(1)  : 0  default, no special actions                  ==
  ! ==              1  skip writing the restart files               ==
  ! ==                                                              ==
  ! ==================================================================
  INTEGER, PARAMETER ::   nbentr=10

  INTEGER :: ibench(nbentr)
  ! ==================================================================

END MODULE benc
