MODULE dpot
  USE system,                          ONLY: maxsp

  IMPLICIT NONE

  ! ==================================================================
  ! == INFORMATION ABOUT PSEUDOPOTENTIALS                           ==
  ! ==================================================================
  ! == TKB(MAXSP)  True if Kleinman-Bylander separable form         ==
  ! == LMAX(MAXSP) Gives maximum angular momentum quantum number    ==
  ! == LSKIP(MAXSP) Gives number of unused ang. mom. qu. number     ==
  ! == LLOC(MAXSP)  Gives ang. mom. qu. number which is considered  ==
  ! ==              as local part                                   ==
  ! ==--------------------------------------------------------------==


  TYPE :: dpot_mod_t
     LOGICAL :: tkb(maxsp)
     LOGICAL :: team(maxsp)
     INTEGER :: lmax(maxsp)
     INTEGER :: lskip(maxsp)
     INTEGER :: lloc(maxsp)
  END TYPE dpot_mod_t
  TYPE(dpot_mod_t) :: dpot_mod
  ! ==================================================================

END MODULE dpot
