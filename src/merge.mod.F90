MODULE mergemod
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE




  TYPE :: merge01_t
     REAL(real_8) :: mshift(3)
     CHARACTER(len=20) :: mfiln1
     CHARACTER(len=20) :: mfiln2
  END TYPE merge01_t
  TYPE(merge01_t) :: merge01
  TYPE :: merge02_t
     INTEGER :: mnst1
     INTEGER :: mnst2
     INTEGER :: mortho
  END TYPE merge02_t
  TYPE(merge02_t) :: merge02
END MODULE mergemod
