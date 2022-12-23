MODULE mm_cpmd_add_MM_forces_f77_utils
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE system,                          ONLY: maxsys

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mm_cpmd_add_mm_forces_f77

CONTAINS

  ! ==================================================================
  SUBROUTINE mm_cpmd_add_mm_forces_f77 (fion)
    ! ==================================================================
    ! This is simply a wrapper subroutine to pass variables 
    ! from the common blocks into the fortran90 subroutine, for
    ! which I do not want to have any common blocks.  Okay this
    ! is really snobbish, and probably unnecessary since
    ! my programming is all crap anyways....
    ! ==================================================================
    ! ==================================================================
    REAL(real_8)                             :: fion(:,:,:)

#if defined (__QMECHCOUPL)

    CALL mm_cpmd_add_mm_forces (fion, ions0%na, ions1%nsp, maxsys%nax)

#endif

    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE mm_cpmd_add_mm_forces_f77

END MODULE mm_cpmd_add_MM_forces_f77_utils
