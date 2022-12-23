MODULE andr
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==================================================================
  ! == ANDERSON AND ALEXANDER MIXING INCLUDE FILE                   ==
  ! ==================================================================
  ! == MAXMIX: MAXIMAL NUMBER OF CHOICE FOR ANDRMIX                 ==
  ! ==--------------------------------------------------------------==
  INTEGER, PARAMETER :: maxmix=10
  ! ==--------------------------------------------------------------==
  ! == NTABMIX: NUMBER OF CHOICE FOR ANDRMIX                        ==
  ! == ANDRMIX: ANDERSOM MIXING PARAMETER (BY DEFAULT=0.2)          ==
  ! == DENSMIX(MAXMIX): DENSITY THRESHOLD ARRAY                     ==
  ! == ANDMIX(MAXMIX): ANDRMIX ARRAY                                ==
  ! ==--------------------------------------------------------------==



  TYPE :: andr2_t
     REAL(real_8) :: andrmix
     REAL(real_8) :: amprd
     REAL(real_8) :: alxmix
     REAL(real_8) :: densmix(maxmix)
     REAL(real_8) :: andmix(maxmix)
     INTEGER :: ntabmix
     LOGICAL :: trand
  END TYPE andr2_t
  TYPE(andr2_t) :: andr2
  ! ==--------------------------------------------------------------==
  ! == NRDIISMAX: MAXIMUM NUMBER OF DENSITY ARRAYS FOR DIIS MIXING  ==
  ! == NRDIIS: CURRENT NUMBER OF DENSITY ARRAYS FOR DIIS MIXING     ==
  ! ==--------------------------------------------------------------==


  TYPE :: andr3_t
     REAL(real_8) :: densnrdiis(maxmix)
     INTEGER :: nrdiismax
     INTEGER :: nrdiis
     INTEGER :: ntabnrdiis
     INTEGER :: inrdiis
     INTEGER :: tabnrdiis(maxmix)
  END TYPE andr3_t
  TYPE(andr3_t) :: andr3

END MODULE andr
