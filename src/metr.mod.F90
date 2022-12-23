MODULE metr
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==================================================================
  ! ==  METRIC TENSOR OF SUPERCELL                                  ==
  ! ==================================================================
  ! == HT    = (A1 A2 A3): metric tensor of supercell               ==
  ! == HTM1  = HT^-1 = (B1 B2 B3)                                   ==
  ! == HTVEL Cell velocities                                        ==
  ! == HTFOR Ionic forces on cell in A1,A2,A3 basis                 ==
  ! ==       (stress tensor - required stress tensor)               ==
  ! == HTFP  Ionic forces on cell in cartesian coordinates          ==
  ! ==       (stress tensor - required stress tensor)               ==
  ! == HT0   Metric tensor of the reference supercell               ==
  ! == HTM10 = HT0^_1                                               ==
  ! ==================================================================

  ! ==================================================================
  REAL(real_8), SAVE :: eps=0.0_real_8,veps=0.0_real_8
  ! ==================================================================
  TYPE :: metr_com_t
     REAL(real_8) :: ht(3,3)=0.0_real_8
     REAL(real_8) :: htm1(3,3)=0.0_real_8
     REAL(real_8) :: htvel(3,3)=0.0_real_8
     REAL(real_8) :: htfor(3,3)=0.0_real_8
     REAL(real_8) :: htfp(3,3)=0.0_real_8
     REAL(real_8) :: ht0(3,3)=0.0_real_8
     REAL(real_8) :: htm10(3,3)=0.0_real_8
  END TYPE metr_com_t
  TYPE(metr_com_t), SAVE :: metr_com

END MODULE metr
