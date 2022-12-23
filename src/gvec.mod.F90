MODULE gvec
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==================================================================
  ! == DEFINE G-VECTOR SETS                                         ==
  ! ==================================================================
  ! == EPSGX: PRECISION NUMBER                                      ==
  ! == EPSG:  USE FOR THE G-VECTOR REORDERING                       ==
  ! ==--------------------------------------------------------------==
  REAL(real_8), PARAMETER :: epsgx=1.e-12_real_8 
  REAL(real_8) :: epsg
  ! ==--------------------------------------------------------------==
  ! == B1(3),B2(3),B3(3): RECIPROCAL BASIS VECTOR                   ==
  ! == GCUT: DENSITY CUTOFF ENERGY                                  ==
  ! == GCUTW: PLANE WAVE CUTOFF ENERGY                              ==
  ! ==--------------------------------------------------------------==


  ! ==--------------------------------------------------------------==
  ! == TNHGWISH: TRUE IF THE NUMBER OF PW COMPONENTS FOR DENSITY    ==
  ! ==        IS REQUIRED                                           ==
  ! == NHGWISH: IF NON-ZERO, NUMBER OF PW COMPONENTS                ==
  ! ==          FOR THE DENSITY REQUIRED                            ==
  ! ==--------------------------------------------------------------==
  INTEGER :: nhgwish
  LOGICAL :: tnhgwish
  ! ==================================================================
  TYPE :: gvec_com_t
     REAL(real_8) :: b1(3)
     REAL(real_8) :: b2(3)
     REAL(real_8) :: b3(3)
     REAL(real_8) :: gcut
     REAL(real_8) :: gcutw
  END TYPE gvec_com_t
  TYPE(gvec_com_t) :: gvec_com

END MODULE gvec
