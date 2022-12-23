MODULE isos
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==================================================================
  ! == ISOLATED SYSTEM (CLUSTER,SURFACE,POLYMER)                    ==
  ! ==================================================================
  TYPE :: isos1_t
     LOGICAL :: tclust
     LOGICAL :: tisos
     LOGICAL :: tcent
     LOGICAL :: toned
     LOGICAL :: ttwod
     LOGICAL :: twall
  END TYPE isos1_t
  TYPE(isos1_t) :: isos1
  TYPE :: isos2_t
     REAL(real_8) :: rcerf
     REAL(real_8) :: alphal
     REAL(real_8) :: duals
     REAL(real_8) :: wall_skin
  END TYPE isos2_t
  TYPE(isos2_t) :: isos2
  TYPE :: isos3_t
     INTEGER :: ps_type
     INTEGER :: snormal
  END TYPE isos3_t
  TYPE(isos3_t) :: isos3

END MODULE isos
