MODULE condu
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==================================================================
  ! == CONDUCTIVITY:                                                ==
  ! ==    TCONDUCT if true calculates conductivity                  ==
  ! ==    NCONDUCT is the dimension of arrays                       ==
  ! ==    ICONDUCT is the number of steps                           ==
  ! ==             between two calculations of conductivity         ==
  ! ==    CONDSTEP histogram bin width in eV.                       ==
  ! ==    FILECOND Name of the file to store arrays                 ==
  ! ==    IUCOND   Unit of the file                                 ==
  ! ==                                                              ==
  ! ==    CONDUCT (NCONDUCT)  conductivity array                    ==
  ! ==    CONDUCT2(NCONDUCT) square conductivity array              ==
  ! ==    CONDAVG (NCONDUCT)  average conductivity array            ==
  ! ==    CONDAVG2(NCONDUCT) average square conductivity array      ==
  ! ==================================================================
  INTEGER, PARAMETER :: iucond=63 
  ! ==--------------------------------------------------------------==




  ! ==--------------------------------------------------------------==
  INTEGER, ALLOCATABLE :: normcon(:)
  REAL(real_8), ALLOCATABLE :: conduct(:)
  REAL(real_8), ALLOCATABLE :: conduct2(:)
  REAL(real_8), ALLOCATABLE :: condavg(:)
  REAL(real_8), ALLOCATABLE :: condavg2(:)

  ! ==================================================================
  TYPE :: condpa_t
     REAL(real_8) :: condstep
     INTEGER :: nconduct
     INTEGER :: iconduct
     INTEGER :: nfconduct
     LOGICAL :: tconduct
     CHARACTER(len=100) :: filecond
  END TYPE condpa_t
  TYPE(condpa_t) :: condpa

END MODULE condu
