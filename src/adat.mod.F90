MODULE adat
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==================================================================
  ! == Atomic data (see atoms.F)                                    ==
  ! ==--------------------------------------------------------------==
  ! == EL  Symbol of atomic elements                                ==
  ! ==--------------------------------------------------------------==
  TYPE :: elem_t
     CHARACTER(len=2) :: el(99)
  END TYPE elem_t
  TYPE(elem_t) :: elem
  ! ==--------------------------------------------------------------==
  ! == COVRAD Covalent radius                                       ==
  ! == ATWT   Atomic Weight                                         ==
  ! == DEFRAG                                                       ==
  ! ==--------------------------------------------------------------==
  REAL(real_8) :: covrad(99),atwt(99),defrag(99),mm_core_raggio(99)
  TYPE :: adat_mod_t
     REAL(real_8) :: covrad(99)
     REAL(real_8) :: atwt(99)
     REAL(real_8) :: defrag(99)
     REAL(real_8) :: mm_core_raggio(99)
  END TYPE adat_mod_t
  TYPE(adat_mod_t) :: adat_mod
  ! ==--------------------------------------------------------------==
  ! == NELCON  Electronic configuration of atom                     ==
  ! ==--------------------------------------------------------------==
  INTEGER :: nelcon(4,7,99)
  TYPE :: atco_t
     INTEGER :: nelcon(4,7,99)
  END TYPE atco_t
  TYPE(atco_t) :: atco
  ! ==================================================================

END MODULE adat
