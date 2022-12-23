MODULE elct
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==================================================================
  ! == ELECTRONIC STATES                                            ==
  ! ==================================================================
  ! == CHARGE: CHARGE OF THE SYSTEM                                 ==
  ! == F(N,NKPNT): OCCUPATION NUMBER                                ==
  ! == NEL: TOTAL NUMBER OF ELECTRONS (now real(8))                 ==
  ! == N: NUMBER OF STATES                                          ==
  ! ==--------------------------------------------------------------==
  TYPE :: crge_t
     REAL(real_8) :: charge
     REAL(real_8) :: nel
     REAL(real_8), ALLOCATABLE :: f(:,:)
     INTEGER :: n
  END TYPE crge_t
  TYPE(crge_t), SAVE :: crge

END MODULE elct
