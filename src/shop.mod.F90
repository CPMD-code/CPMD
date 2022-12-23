MODULE shop
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! surface hopping.
  CHARACTER(len=20) :: s0_filn,s1_filn
  LOGICAL :: tlsd0,tshopold
  REAL(real_8), ALLOCATABLE :: fs0(:)
  REAL(real_8), ALLOCATABLE :: fs1(:)





  TYPE :: sh02_t
     REAL(real_8) :: eaddsh
     INTEGER :: nst_s0
     INTEGER :: nst_s1
     INTEGER :: nsttot
     INTEGER :: nsurf
     INTEGER :: nelb2
  END TYPE sh02_t
  TYPE(sh02_t) :: sh02

END MODULE shop
