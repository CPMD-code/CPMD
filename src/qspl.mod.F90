MODULE qspl
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==================================================================
  ! == INCLUDE FILE FOR NON-LOCAL PSEUDOPOTENTIAL                   ==
  ! ==================================================================
  ! == QSPLINE=.TRUE. Use SPLINE INTERPOLATION OF Q-FUNCTIONS       ==
  ! ==                For Vanderbilt pseudopotentials               ==
  ! == INITSPL=.TRUE. Always TRUE                                   ==
  ! ==         SPLINE INTERPOLATION IN G-SPACE FOR PP FUNCTIONS     ==
  ! ==--------------------------------------------------------------==

  ! ==================================================================
  ! == NSPLPO  Number of points for the spline interpolation        ==
  ! == NSPLPA  First spline point calculated on this cpu            ==
  ! == NSPLPE  Last spline point calculated on this cpu             ==
  ! == NQDIM   Number of points for the spline int. for Q-function  ==
  ! == GGNG(1:NSPLPO) Point grid from 0 to GCUTW+GCUTKA             ==
  ! == GGNH(1:NSPLPO) Point grid from 0 to GCUT                     ==
  ! == VOO(1:NSPLPO,1:2,1:NSP) Grid values for local part of pp     ==
  ! == TWNS(1:NSPLPO,1:2,1:NGH(IS),1:NSP) Grid of values            ==
  ! ==     in G-space of non-local potential part for each species  ==
  ! ==     and each l,m components                                  ==
  ! ==     index 2 is the derivatives (use with spline interpolat.) ==
  ! ==--------------------------------------------------------------==
  INTEGER :: nsplpo,nqdim,nsplpa,nsplpe
  REAL(real_8), ALLOCATABLE :: voo(:,:,:)
  REAL(real_8), ALLOCATABLE :: twns(:,:,:,:)
  REAL(real_8), ALLOCATABLE :: ggnh(:)
  REAL(real_8), ALLOCATABLE :: ggng(:)

  ! ==================================================================
  ! == QSRANG Range for Q-functions used with spline interpolation  ==
  ! ==--------------------------------------------------------------==
  REAL(real_8) :: qsrang
  ! ==================================================================
  TYPE :: qspl1_t
     LOGICAL :: qspline
     LOGICAL :: initspl
  END TYPE qspl1_t
  TYPE(qspl1_t) :: qspl1

END MODULE qspl
