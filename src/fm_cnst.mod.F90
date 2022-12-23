! ==================================================================
! == INCLUDE FILE for forcematching, unit conversion factors      ==
! ==--------------------------------------------------------------==
! == KJM_AU  kJ/mol to ATOMIC UNITS                               ==
! == NM_BOHR nanometer to BOHR
! ==================================================================
MODULE fm_cnst
  USE cnst,                            ONLY: au_kjm,&
                                             fbohr
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE
  ! 
  REAL(real_8), PARAMETER :: kjm_au = 1._real_8 / au_kjm 
  REAL(real_8), PARAMETER :: nm_bohr = 10._real_8 * fbohr 
  REAL(real_8), PARAMETER :: bohr_nm = 1._real_8 / nm_bohr 

  ! conversion factors for bond force constants
  REAL(real_8), PARAMETER :: kb_grm2au = kjm_au / nm_bohr**4 
  REAL(real_8), PARAMETER :: kb_au2grm = 1._real_8 / kb_grm2au 
  REAL(real_8), PARAMETER :: kb_amb2au = kjm_au / nm_bohr**2 
  REAL(real_8), PARAMETER :: kb_au2amb = 1._real_8 / kb_amb2au 

END MODULE fm_cnst
! ==================================================================
