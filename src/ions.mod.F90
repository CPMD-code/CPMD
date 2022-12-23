MODULE ions
  USE kinds,                           ONLY: real_8
  USE system,                          ONLY: maxsp

  IMPLICIT NONE

  ! ==================================================================
  ! == PSEUDOPOTENTIALS AND IONS                                    ==
  ! == MAXGAU: MAXIMUM GAUSSIAN NUMBER FOR PP                       ==
  ! ==================================================================
  INTEGER, PARAMETER :: maxgau=3
  ! ==================================================================
  TYPE :: ions0_t
     REAL(real_8) :: zv(maxsp)
     REAL(real_8) :: wrc(maxsp,2)
     REAL(real_8) :: rc(maxsp,3)
     INTEGER :: iatyp(maxsp)
     INTEGER :: na(maxsp)
     INTEGER :: igau(maxsp)
     INTEGER :: nambl(maxsp)
  END TYPE ions0_t
  TYPE(ions0_t) :: ions0
  ! ==================================================================
  ! == AL(MAXGAU,NSX,LMAXX)                                         ==
  ! == BL(MAXGAU,NSX,LMAXX)                                         ==
  ! == RCL(MAXGAU,NSX,LMAXX)                                        ==
  ! ==--------------------------------------------------------------==
  REAL(real_8), ALLOCATABLE :: al(:,:,:)
  REAL(real_8), ALLOCATABLE :: bl(:,:,:)
  REAL(real_8), ALLOCATABLE, TARGET :: RCL(:,:,:) ! ! TARGET = DEBUG

  TYPE :: ions1_t
     INTEGER :: nat
     INTEGER :: nsp
     INTEGER :: nspnl
     INTEGER :: nsanl
  END TYPE ions1_t
  TYPE(ions1_t) :: ions1
  ! ==================================================================
  ! == For finite difference: option only for some atoms            ==
  ! == IF TREF_FDIFF == TRUE.:                                      ==
  ! == COORD_FDIFF gives the reference center                       ==
  ! == R_FDIFF gives the radius                                     ==
  ! == We calculate the finite diffrence only for atoms inside      ==
  ! == a sphere of radius "R_FDIFF" centered at the "COORD_FDIFF".  ==
  ! ==--------------------------------------------------------------==
  LOGICAL :: tref_fdiff
  REAL(real_8) :: coord_fdiff(3)
  REAL(real_8) :: r_fdiff
  ! ==================================================================

END MODULE ions
