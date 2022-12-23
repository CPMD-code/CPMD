MODULE broy
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==================================================================
  ! == BROYDEN MIXING PARAMETERS                                    ==
  ! ==--------------------------------------------------------------==
  ! == NFRBROY  BROYDEN MIXING STARTS AFTER NFRBROY STEPS           ==
  ! == NGBROY   NUMBER OF PLANE WAVES FOR BROYDEN MIXING            ==
  ! == IBRESET  BROYDEN MIXING RESET AFTER IBRESET STEPS            ==
  ! == BROYMIX  BROYDEN MIXING PARAMETER                            ==
  ! == ECUTBROY BROYDEN MIXING CUTOFF                               ==
  ! == W02BROY  BROYDEN MIXING W02                                  ==
  ! == KERMIX   KERKER  MIXING PARAMETER                            ==
  ! == TGMIX    .TRUE.  MIXING IN G-SPACE                           ==
  ! == TGBROY   .TRUE.  BROYDEN MIXING                              ==
  ! == TADBROY  .TRUE.  ADAPTIVE WEIGHTING                          ==
  ! ==--------------------------------------------------------------==

  TYPE :: broy1_t
     REAL(real_8) :: broymix
     REAL(real_8) :: ecutbroy
     REAL(real_8) :: w02broy
     INTEGER :: nfrbroy
     INTEGER :: ngbroy
     REAL(real_8) :: kermix
     INTEGER :: ibreset
     LOGICAL :: tgmix
     LOGICAL :: tgbroy
     LOGICAL :: tadbroy
  END TYPE broy1_t
  TYPE(broy1_t) :: broy1

END MODULE broy
