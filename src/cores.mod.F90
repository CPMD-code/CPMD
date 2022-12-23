MODULE cores
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==================================================================
  ! == DYNAMICAL ALLOCATION OF ARRAYS RELATED TO CORE SPECTRA STUFF ==
  ! ==================================================================



  ! ==--------------------------------------------------------------==
  ! ==================================================================
  REAL(real_8), ALLOCATABLE :: core_atwfr(:)
  REAL(real_8), ALLOCATABLE :: core_cat(:,:)

  COMPLEX(real_8), ALLOCATABLE :: core_c0(:,:)


  ! ==================================================================

  TYPE :: coresi_t
     INTEGER :: core_atom
     INTEGER :: core_nqsto
     INTEGER :: core_lshell
  END TYPE coresi_t
  TYPE(coresi_t) :: coresi
  TYPE :: coresl_t
     LOGICAL :: tcores
  END TYPE coresl_t
  TYPE(coresl_t) :: coresl
  TYPE :: coresr_t
     REAL(real_8) :: core_level
     REAL(real_8) :: core_stoexp
  END TYPE coresr_t
  TYPE(coresr_t) :: coresr

END MODULE cores
