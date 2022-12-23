MODULE cell
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==================================================================
  ! == CELLDM: PARAMETERS FOR THE PERIODIC SUPERCELL                ==
  ! ==         Dx Dy Dz alpha(yz) beta(zx) gamma(xy)                ==
  ! ==--------------------------------------------------------------==
  ! ==================================================================
  TYPE :: cell_com_t
     REAL(real_8) :: celldm(6)=HUGE(0.0_real_8)
     REAL(real_8) :: ar1(3)=HUGE(0.0_real_8)
     REAL(real_8) :: ar2(3)=HUGE(0.0_real_8)
     REAL(real_8) :: ar3(3)=HUGE(0.0_real_8)
     REAL(real_8) :: volcel=HUGE(0.0_real_8)
  END TYPE cell_com_t
  TYPE(cell_com_t), SAVE :: cell_com
  TYPE :: lcell_t
     LOGICAL :: tcellvectors = .FALSE.
     LOGICAL :: tcellrefvec = .FALSE.
  END TYPE lcell_t
  TYPE(lcell_t), SAVE :: lcell

END MODULE cell
