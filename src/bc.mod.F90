MODULE bc
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==================================================================

  TYPE :: bc_com_t
     REAL(real_8) :: v1crys(4)
     REAL(real_8) :: v2crys(4)
     REAL(real_8) :: v3crys(4)
  END TYPE bc_com_t
  TYPE(bc_com_t) :: bc_com
  ! ==================================================================

END MODULE bc
