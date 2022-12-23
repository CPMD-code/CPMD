MODULE aavan
  USE kinds,                           ONLY: real_8
  USE system,                          ONLY: mx,&
                                             nlx

  IMPLICIT NONE

  ! ==================================================================
  ! == Used by aainit                                               ==
  ! ==--------------------------------------------------------------==
  REAL(real_8) :: ap(25,nlx,nlx)
  INTEGER, ALLOCATABLE :: indv(:,:)
  INTEGER :: lpx(nlx,nlx),lpl(nlx,nlx,mx)

  TYPE :: aavan_mod_t
     REAL(real_8) :: ap(25,nlx,nlx)
     INTEGER :: lpx(nlx,nlx)
     INTEGER :: lpl(nlx,nlx,mx)
  END TYPE aavan_mod_t
  TYPE(aavan_mod_t) :: aavan_mod
  ! ==================================================================

END MODULE aavan
