! 
MODULE shop_rest
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! surface hopping data for RESTART 

  REAL(real_8) :: prob(2),zz

  TYPE :: prob1_t
     REAL(real_8) :: d1sq
     REAL(real_8) :: d2sq
     REAL(real_8) :: d11dot
     REAL(real_8) :: d22dot
  END TYPE prob1_t
  TYPE(prob1_t) :: prob1
  TYPE :: sh03_t
     REAL(real_8) :: pop(6)
     REAL(real_8) :: coupl(2,2)
     REAL(real_8) :: couplold(2,2)
     REAL(real_8) :: det
     REAL(real_8) :: detold
     REAL(real_8) :: ec(2)
     REAL(real_8) :: eold(2)
     INTEGER :: isurf
     INTEGER :: ioldsurf
     INTEGER :: ishrest
     LOGICAL :: tshopres
  END TYPE sh03_t
  TYPE(sh03_t) :: sh03

END MODULE shop_rest
