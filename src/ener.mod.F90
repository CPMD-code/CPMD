MODULE ener
  USE kinds,                           ONLY: real_8
  USE system,                          ONLY: maxsp

  IMPLICIT NONE

  ! ==================================================================
  ! == ENERGY AND CHARGE QUANTITIES                                 ==
  ! ==================================================================
  ! == ETOT:  Total energy                                          ==
  ! == EKIN:  Kinetic energy                                        ==
  ! == EHT:   Electrostatic energy                                  ==
  ! ==  If TDIAG EHII-EHEE+ESR-ESELF   (ELEC2)                      ==
  ! ==  else     EHEP+ESR-ESELF        (ELEC1)                      ==
  ! == EHEP   (e- pseudocharge) density Hartree energy              ==
  ! == EHEE   electron-electron Hartree energy                      ==
  ! == EHII   ion-ion energy                                        ==
  ! == EPSEU: Contribution of local pseudopotential part            ==
  ! == ENL:   Non-local part contribution of pseudopotential        ==
  ! == EEIG:  Band energy (free energy functional)                  ==
  ! == EBAND: band energy (sum of eigen. x occup. numbers)          ==
  ! == ENTROPY: KT*Entropy                                          ==
  ! == AMU:   Chemical potential (i.e. Fermi energy)                ==
  ! == EBOGO: Bogoliubov correction energy (free energy)            ==
  ! == EXC:   Exchange-correlation energy                           ==
  ! == VXC:   Contribution of exchange-correlation potential (-)    ==
  ! == EGC:   Gradient correction energy                            ==
  ! == EEXT:  Energy of external potential(for interface calc.)     ==
  ! == ETDDFT:Excitation energy from TDDFT theory                   ==
  ! == EHSIC: Self interaction energy (Hartree term only)           ==
  ! ==--------------------------------------------------------------==
  ! == ION-ION INTERACTION:                                         ==
  ! == The ion-ion interaction energy, the Hartree and the local    ==
  ! == pseudo potential energy are diverging due to the presence of ==
  ! == slowly decaying Coulomb forces                               ==
  ! == To calculate ion-ion interaction, we replace ionic point     ==
  ! == charges by Gaussian charge distributions determined by raggio==
  ! == (in RHOPS in G space)                                        ==
  ! == E_{Ion-Ion} = Hartree term due to RHOPS -ESELF - ESR         ==
  ! == ESELF: Self-energy (Sum_i Zv(i)^2 /Raggio(i) )               ==
  ! ==        calculated in RINFORCE once                           ==
  ! == ESR:   Residual part of the ion-ion interaction              ==
  ! ==        (nearest neighbor energy)                             ==
  ! ==--------------------------------------------------------------==
  ! == Global summation:                                            ==
  ! == ETOT,EKIN,EPSEU,ENL,EHT,EHEP,EHEE,EHII,EXC,VXC,EGC,ESR       ==
  ! ==--------------------------------------------------------------==
  TYPE :: ener_com_t
     REAL(real_8) :: etot = 0.0_real_8
     REAL(real_8) :: ekin = 0.0_real_8
     REAL(real_8) :: epseu = 0.0_real_8
     REAL(real_8) :: enl = 0.0_real_8
     REAL(real_8) :: eht = 0.0_real_8
     REAL(real_8) :: ehep = 0.0_real_8
     REAL(real_8) :: ehee = 0.0_real_8
     REAL(real_8) :: ehii = 0.0_real_8
     REAL(real_8) :: exc = 0.0_real_8
     REAL(real_8) :: vxc = 0.0_real_8
     REAL(real_8) :: egc = 0.0_real_8
     REAL(real_8) :: esr = 0.0_real_8
     REAL(real_8) :: eeig = 0.0_real_8
     REAL(real_8) :: eband = 0.0_real_8
     REAL(real_8) :: entropy = 0.0_real_8
     REAL(real_8) :: eself = 0.0_real_8
     REAL(real_8) :: ecnstr = 0.0_real_8
     REAL(real_8) :: amu = 0.0_real_8
     REAL(real_8) :: ebogo = 0.0_real_8
     REAL(real_8) :: eext = 0.0_real_8
     REAL(real_8) :: etddft = 0.0_real_8
     REAL(real_8) :: ehsic = 0.0_real_8
     REAL(real_8) :: erestr = 0.0_real_8
     REAL(real_8) :: eefield = 0.0_real_8
  END TYPE ener_com_t
  TYPE(ener_com_t), SAVE :: ener_com
  ! ==--------------------------------------------------------------==
  ! == CAS22 Energies                                               ==
  ! ==--------------------------------------------------------------==
  TYPE :: ener_c_t
     REAL(real_8) :: etot_a = 0.0_real_8
     REAL(real_8) :: ekin_a = 0.0_real_8
     REAL(real_8) :: epseu_a = 0.0_real_8
     REAL(real_8) :: enl_a = 0.0_real_8
     REAL(real_8) :: eht_a = 0.0_real_8
     REAL(real_8) :: exc_a = 0.0_real_8
     REAL(real_8) :: eext_a = 0.0_real_8
     REAL(real_8) :: etot_2 = 0.0_real_8
     REAL(real_8) :: ekin_2 = 0.0_real_8
     REAL(real_8) :: epseu_2 = 0.0_real_8
     REAL(real_8) :: enl_2 = 0.0_real_8
     REAL(real_8) :: eht_2 = 0.0_real_8
     REAL(real_8) :: exc_2 = 0.0_real_8
     REAL(real_8) :: eext_2 = 0.0_real_8
     REAL(real_8) :: etot_ab = 0.0_real_8
     REAL(real_8) :: ekin_ab = 0.0_real_8
     REAL(real_8) :: epseu_ab = 0.0_real_8
     REAL(real_8) :: enl_ab = 0.0_real_8
     REAL(real_8) :: eht_ab = 0.0_real_8
     REAL(real_8) :: exc_ab = 0.0_real_8
     REAL(real_8) :: eext_ab = 0.0_real_8
  END TYPE ener_c_t
  TYPE(ener_c_t), SAVE :: ener_c

  TYPE :: ener_d_t
     REAL(real_8) :: etot_b = 0.0_real_8
     REAL(real_8) :: ecas = 0.0_real_8
     REAL(real_8) :: etot_t = 0.0_real_8
     REAL(real_8) :: casang = 0.0_real_8
  END TYPE ener_d_t
  TYPE(ener_d_t), SAVE :: ener_d
  ! ==--------------------------------------------------------------==
  ! == CSUMG: Total sum in reciprocal space                         ==
  ! == CSUMR: Total sum in real space                               ==
  ! == CSUMS: Integral of spin polarized density                    ==
  ! == CSUMSABS: Integral of abs value of spin polarized density    ==
  ! ==--------------------------------------------------------------==
  TYPE :: chrg_t
     REAL(real_8) :: csumg = 0.0_real_8
     REAL(real_8) :: csumr = 0.0_real_8
     REAL(real_8) :: csums = 0.0_real_8
     REAL(real_8) :: csumsabs = 0.0_real_8
     REAL(real_8) :: vdbchg(maxsp) = 0.0_real_8
  END TYPE chrg_t
  TYPE(chrg_t), SAVE :: chrg
  ! ==--------------------------------------------------------------==
  ! == TENERGY_OK                                                   ==
  ! ==--------------------------------------------------------------==
  LOGICAL, SAVE :: tenergy_ok = .FALSE.
  ! ==================================================================

END MODULE ener
