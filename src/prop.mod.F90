MODULE prop
  USE kinds,                           ONLY: real_8
  USE system,                          ONLY: maxsp

  IMPLICIT NONE

  ! ==================================================================
  ! == CALCULATED PORPERTIES                                        ==
  ! ==--------------------------------------------------------------==
  ! == .TRUE. if the corresponding property has to be calculated    ==
  ! == PWFN    Project wavefunctions on Atomic Orbital              ==
  ! == LOCL    Localized Orbitals                                   ==
  ! == PRTO    Print the wavefunctions in atomic orbital basis      ==
  ! == OPTS    Optimizs Slater exponents                            ==
  ! == LDIP    Dipole moments                                       ==
  ! == LRSDIP  Calculate dipole moments in real space               ==
  ! == DBERRY  Berry phase dipole moment                            ==
  ! == LOCD    Local dipole                                         ==
  ! == LEXT    Excited  dipole                                      ==
  ! == LTRM    Transition moments                                   ==
  ! == MPAN    Mulliken population analysis                         ==
  ! == DPAN    Davidson population analysis                         ==
  ! == CEIG    Need to read the eigenvalues in RESTART file         ==
  ! == ESPC    Atomic charges                                       ==
  ! == TRAMAN  Raman phonons calculation                            ==
  ! == PERTURB Perturbation theory (ramon, phonon)                  ==
  ! == PHONON  Phonon calculation                                   ==
  ! == ELDISP  Electron dispersion function                         ==
  ! == GLOCL   Localization in reciprocal space                     ==
  ! == PYLM    Projection on Spherical Harmonics                    ==
  ! ==--------------------------------------------------------------==
  TYPE :: prop1_t
     LOGICAL :: pwfn
     LOGICAL :: wfr1
     LOGICAL :: wfr2
     LOGICAL :: wfr3
     LOGICAL :: cor1
     LOGICAL :: cor2
     LOGICAL :: abas
     LOGICAL :: locl
     LOGICAL :: prto
     LOGICAL :: opts
     LOGICAL :: ldip
     LOGICAL :: locd
     LOGICAL :: lext
     LOGICAL :: mpan
     LOGICAL :: dpan
     LOGICAL :: ceig
     LOGICAL :: espc
     LOGICAL :: eldisp
     LOGICAL :: glocl
     LOGICAL :: dberry
     LOGICAL :: ltrm
     LOGICAL :: lrsdip
     LOGICAL :: tavgp
     LOGICAL :: pylm
  END TYPE prop1_t
  TYPE(prop1_t) :: prop1
  ! ==--------------------------------------------------------------==
  ! == Population analysis                                          ==
  ! == NUMORB                                                       ==
  ! == NCEN  Highest n-center terms included                        ==
  ! == MAOS  For Davidson PAN                                       ==
  ! ==--------------------------------------------------------------==
  TYPE :: prop2_t
     INTEGER :: numorb
     INTEGER :: ncen
     INTEGER :: maos(maxsp)
  END TYPE prop2_t
  TYPE(prop2_t) :: prop2
  ! ==--------------------------------------------------------------==
  ! == For population analysis                                      ==
  ! == CUT3O 3-center cutoff                                        ==
  ! == CUT4O 4-center cutoff                                        ==
  ! ==--------------------------------------------------------------==
  TYPE :: prop3_t
     REAL(real_8) :: cut3o
     REAL(real_8) :: cut4o
  END TYPE prop3_t
  TYPE(prop3_t) :: prop3
  ! ==--------------------------------------------------------------==
  ! == For plotting CUBE files                                      ==
  ! == TCUBEFILE_ORB  = write KS-orbitals in cubefiles              ==
  ! == TCUBEFILE_DENS = write the electron density as cubefile      ==
  ! == TCUBEFILE_POT  = write the potential as cubefile             ==
  ! == THALFMESH      = write cube with half the points per         ==
  ! ==                  direction -> 1/8th of size                  ==
  ! == TCUBECENTER    = cubefile center is set from input           ==
  ! == CUBECIN        = cube center as read in (in a.u.)            ==
  ! == NUM_CUBEORB    = number of orbitals to write                 ==
  ! == ICUBEORB       = list of orbitals to write                   ==
  ! ==--------------------------------------------------------------==
  TYPE :: prop4_t
     REAL(real_8) :: cubecin(3)
     LOGICAL :: tcubefile_orb
     LOGICAL :: tcubefile_dens
     LOGICAL :: tcubefile_pot
     LOGICAL :: thalfmesh
     LOGICAL :: tcubecenter
     INTEGER :: num_cubeorb
  END TYPE prop4_t
  TYPE(prop4_t) :: prop4
  INTEGER, ALLOCATABLE :: icubeorb(:)

  REAL(real_8), ALLOCATABLE :: z_11(:)

  ! ==================================================================
  ! FOR COMPUTING EFG TENSOR AND EPR ANISOTROPIC HF TENSOR
  ! ==================================================================
  TYPE :: prop5_t
     LOGICAL :: teprefg
     LOGICAL :: tefg
  END TYPE prop5_t
  TYPE(prop5_t) :: prop5
  REAL(real_8), ALLOCATABLE :: rho_save_epr(:)
  REAL(real_8), ALLOCATABLE :: rho_save_efg(:)

  ! ==================================================================
  ! FOR COMPUTING AVERAGED ATOMIC ELECTROSTATIC POTENTIAL
  ! ==================================================================
  REAL(real_8) :: rcut

  ! ==================================================================
  ! ==================================================================
  ! FOR COMPUTING PROJECTION ON SPHERICAL HARMONICS
  ! ==================================================================
  TYPE :: prop7_t
     REAL(real_8) :: centylm(3)
     REAL(real_8) :: rylmax
     REAL(real_8) :: spr_min
     REAL(real_8) :: spr_max
     INTEGER :: numylm
     INTEGER :: nylmax
     LOGICAL :: radial
     LOGICAL :: spread_ham
  END TYPE prop7_t
  TYPE(prop7_t) :: prop7
END MODULE prop
