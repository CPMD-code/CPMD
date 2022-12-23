MODULE dg
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! CASPUR 02/05/2004
  ! ==================================================================
  ! TDG = .TRUE. Use the DOUBLEGRID

  ! ==================================================================
  ! Identifier of the FFT grid (returned by ADDFFTNSET into DG_INIT)
  INTEGER :: ipooldg
  ! ==================================================================
  ! Wavefunction and density cutoff
  REAL(real_8) :: ecutdg
  ! ==================================================================

  TYPE :: edgcomm_t
     REAL(real_8) :: ecutwdg
     REAL(real_8) :: ecutdg
  END TYPE edgcomm_t
  TYPE(edgcomm_t) :: edgcomm
  TYPE :: tdgcomm_t
     LOGICAL :: tdg
  END TYPE tdgcomm_t
  TYPE(tdgcomm_t) :: tdgcomm

END MODULE dg
