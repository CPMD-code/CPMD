MODULE cnst
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==================================================================
  ! == INCLUDE FILE OF NUMERIC CONSTANTS                            ==
  ! == last update: 2008-03-25 by akohlmey@cmm.chem.upenn.edu       ==
  ! == CODATA 2006, Version 5.1, http://physics.nist.gov/constants  ==
  ! ==--------------------------------------------------------------==
  ! == UIMAG = CMPLX(0._real_8,1._real_8)                          ==
  ! == PI      PI NUMBER                                            ==
  ! == FPI     4*PI                                                 ==
  ! == RY      RYDBERG IN ELECTRON-VOLT                             ==
  ! == EVTOKEL ELECTRON-VOL IN KELVIN                               ==
  ! == FACTEM  1 HARTREE IN KELVIN                                  ==
  ! == SCMASS  ATOMIC MASS UNIT IN ATOMIC UNITS                     ==
  ! == FBOHR   ANGSTROM TO ATOMIC UNITS (BOHR)                      ==
  ! == AU_K    ATOMIC UNITS TO KBAR                                 ==
  ! == KB_AU   KBAR TO ATOMIC UNITS                                 ==
  ! == AU_FS   ATOMIC UNITS TO FEMTOSECONDS                         ==
  ! == AU_DEB  ATOMIC UNITS TO DEBYE                                ==
  ! == AU_KJM  ATOMIC UNITS TO kJ/mol                               ==
  ! == AU_KCM  ATOMIC UNITS TO kcal/mol                             ==
  ! == KBOLTZ  kB IN ATOMIC UNITS                                   ==
  ! ==================================================================
  COMPLEX(real_8), PARAMETER ::  uimag=(0._real_8,1._real_8)  
  ! NEW UNITS
  REAL(real_8), PARAMETER :: pi=3.141592653589793_real_8, fpi=4.0_real_8*pi,&
       ry=13.60569193_real_8, evtokel=11604.505_real_8, factem=2._real_8*ry*evtokel,&
       scmass=1822.888485_real_8, fbohr=1._real_8/0.529177210859_real_8,&
       au_deb=2.5417462289_real_8, au_kb=294210.1080_real_8, kb_au=1._real_8/au_kb,&
       au_fs=2.418884326505e-2_real_8, kboltz=0.316681534e-5_real_8,&
       au_kjm=2.62549962505e3_real_8, au_kcm=6.275094706142e2_real_8
  CHARACTER(len=11), PARAMETER :: unit_txt='CODATA 2006'
  ! ==================================================================

END MODULE cnst
