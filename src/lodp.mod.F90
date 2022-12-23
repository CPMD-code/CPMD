MODULE lodp
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==================================================================
  ! == COMMON BLOCK FOR LOCAL DIPOLE MOMENTS                        ==
  ! ==--------------------------------------------------------------==
  INTEGER, PARAMETER :: nmaxld=20 
  REAL(real_8) :: rcc(3,nmaxld),dmomlo(3,nmaxld),xminld(nmaxld),yminld(&
       nmaxld),zminld(nmaxld),xmaxld(nmaxld),ymaxld(nmaxld),zmaxld(&
       nmaxld),chrld(nmaxld)
  ! ==================================================================
  INTEGER :: numbld
  ! ==================================================================
  ! == COMMON BLOCK FOR EXCITED STATE DIPOLE MOMENTS                ==
  ! ==================================================================
  INTEGER, ALLOCATABLE :: nsdip(:,:)

  REAL(real_8), ALLOCATABLE :: focc(:)
  REAL(real_8), ALLOCATABLE :: exd(:,:,:)
  REAL(real_8), ALLOCATABLE :: trmom(:,:,:)


  ! ==================================================================
  TYPE :: extd_t
     INTEGER :: nexdip
     INTEGER :: ntrmom
  END TYPE extd_t
  TYPE(extd_t) :: extd

END MODULE lodp
