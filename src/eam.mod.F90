MODULE eam
  USE kinds,                           ONLY: real_8
  USE system,                          ONLY: maxsp

  IMPLICIT NONE

  LOGICAL :: tieam,team_cross(maxsp,maxsp)

  TYPE :: eam2_t
     REAL(real_8) :: eamre(maxsp)
     REAL(real_8) :: eamfe(maxsp)
     REAL(real_8) :: eampe(maxsp)
     REAL(real_8) :: eama(maxsp)
     REAL(real_8) :: eamb(maxsp)
     REAL(real_8) :: eamc(maxsp)
     REAL(real_8) :: eamrho(maxsp)
     REAL(real_8) :: eambpe(maxsp)
     REAL(real_8) :: eamec(maxsp)
     REAL(real_8) :: eamcut(maxsp)
     REAL(real_8) :: eamsmooth(maxsp)
  END TYPE eam2_t
  TYPE(eam2_t) :: eam2

END MODULE eam
