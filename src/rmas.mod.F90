MODULE rmas
  USE kinds,                           ONLY: real_8
  USE system,                          ONLY: maxsp

  IMPLICIT NONE

  ! ==================================================================
  ! == INFORMATION ABOUT ATOMIC MASS                                ==
  ! ==================================================================
  TYPE :: rmass_t
     REAL(real_8) :: pmat0 = HUGE(0.0_real_8)
     REAL(real_8) :: pmatot = HUGE(0.0_real_8)
     REAL(real_8) :: pma(maxsp) = 0.0_real_8
     REAL(real_8) :: pma0(maxsp) = 0.0_real_8 !vw this is checked as == 0.0_real_8 in many different places !
  END TYPE rmass_t
  TYPE(rmass_t), SAVE :: rmass

END MODULE rmas
