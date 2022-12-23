MODULE store_types
  IMPLICIT NONE

  ! ==================================================================
  ! == FILE IN WHICH WE HAVE FLAGS                                  ==
  ! == FOR RESTART, PRINT AND STORE OPTIONS                         ==
  ! ==================================================================
  ! == MRESTART: Number of parameters (for parallel)                ==
  ! == RESTART: True if restart options                             ==
  ! == RWF:   Wavefunctions                                         ==
  ! == RCO:   Coordinates                                           ==
  ! == RVEL:  Velocities                                            ==
  ! == RAC:   Accumulators                                          ==
  ! == RHE:   Hessian matrix                                        ==
  ! == RNOE:  NOSEE                                                 ==
  ! == RNOP:  NOSEP                                                 ==
  ! == RNOC:  NOSEC                                                 ==
  ! == RGEO:  GEOFILE                                               ==
  ! == RLATE: Latest restart                                        ==
  ! == RVIB:  Vibrational properties                                ==
  ! == RCELL: Cell                                                  ==
  ! == RPOT:  Potential                                             ==
  ! == RRHO:  Density                                               ==
  ! == RKPT:  Kpoints                                               ==
  ! == ROCC:  Occupation numbers                                    ==
  ! == ROH:   Orbital Hardness                                      ==
  ! == RCON:  Constraints                                           ==
  ! == RXTP:  Extrapolation                                         ==
  ! == RPRNG: Pseudo-random number generator                        ==
  ! == RGLE:  (generalized) Langevin equation                       ==
  ! == RNON:  NOT RESTART                                           ==
  ! ==--------------------------------------------------------------==
  ! ==                     GENERAL FILE FORMAT (FOR IREC ALSO)      ==
  ! ==--------------------------------------------------------------==
  ! ==  Section  1: Header                                          ==
  ! ==  Section  2: Symmetry and Cell info                          ==
  ! ==  Section  3: Number of atomic species and atoms per species  ==
  ! ==  Section  4: Atomic coordinates                              ==
  ! ==  Section  5: Atomic velocities                               ==
  ! ==  Section  6: Initial atomic coordinates                      ==
  ! ==  Section  7: Cutoff, # of electrons, grid                    ==
  ! ==  Section  8: States, dual, plane waves                       ==
  ! ==  Section  9: PW Coefficients                                 ==
  ! ==  Section 10: PW Velocities                                   ==
  ! ==  Section 11: Accumulators                                    ==
  ! ==  Section 12: Nose thermostats general info                   ==
  ! ==  Section 13: Nose thermostats electrons                      ==
  ! ==  Section 14: Nose thermostats ions                           ==
  ! ==  Section 15: Nose thermostats ions (ULTRA)                   ==
  ! ==  Section 16: Nose thermostats ions (MASSIVE)                 ==
  ! ==  Section 17: Nose thermostats cell                           ==
  ! ==  Section 18: Potential                                       ==
  ! ==  Section 19: PW single states                                ==
  ! ==  Section 20: H Matrices (cell parameters)                    ==
  ! ==  Section 21: K-Points                                        ==
  ! ==  Section 22: Electron density                                ==
  ! ==  Section 23: Occupation numbers                              ==
  ! ==  Section 24: Fermi energy and eigenvalues                    ==
  ! ==  Section 25: Classical Particles (coordinates and velocities)==
  ! ==  Section 26: LinRes PW Coefficients                          ==
  ! ==  Section 27: LinRes PW Velocities                            ==
  ! ==  Section 28: Partial Hessian (for microiterative TS search)  ==
  ! ==  Section 29: L-BFGS history and status                       ==
  ! ==  Section 30: P-RFO status                                    ==
  ! ==  Section 31: Adaptive tolerance status                       ==
  ! ==  Section 32: Constraints                                     ==
  ! ==  Section 33: Cell translation vector in QM/MM runs           ==
  ! ==  Section 34: Local temperature parameters                    ==
  ! ==  Section 35: Surface Hopping I                               ==
  ! ==  Section 36: Surface Hopping II                              ==
  ! ==  Section 37: Pseudo-random number generator                  ==
  ! ==  Section 38: (generalized) Langevin equation                 ==
  ! ==  Section 39 - 98 : empty                                     ==
  ! ==  Section 99: Wavefunction history for extrapolation          ==
  ! ==--------------------------------------------------------------==
  INTEGER, PARAMETER :: irec_info = 1 ,   irec_cell = 2 ,   irec_nas  = 3 ,&
       irec_co   = 4 ,   irec_vel  = 5 ,   irec_ico  = 6 ,&
       irec_cut  = 7 ,   irec_sdp  = 8 ,   irec_wf   = 9 ,&
       irec_pwv  = 10,   irec_ac   = 11,   irec_nose = 12,&
       irec_noe  = 13,   irec_nop1 = 14,   irec_nop2 = 15,&
       irec_nop3 = 16,   irec_noc  = 17,   irec_pot  = 18,&
       irec_pws  = 19,   irec_he   = 20,   irec_kpt  = 21,&
       irec_rho  = 22,   irec_occ  = 23,   irec_eigv = 24,&
       irec_clas = 25,   irec_lrwf = 26,   irec_lrwv = 27,&
       irec_phes = 28,   irec_prfo = 29,   irec_lbfgs= 30,&
       irec_rdtl = 31,   irec_cons = 32,   irec_ctrans = 33,&
       irec_nop4 = 34,   irec_shop = 35,   irec_shopbo = 36,&
       irec_prng = 37,   irec_gle = 38,  irec_xtrp = 99
  ! ==================================================================
  INTEGER, PARAMETER :: mrestart=26 
  ! ==================================================================
  ! == STORE OPTIONS                                                ==
  ! == ISTORE: Store at each ISTORE step                            ==
  ! == ISCTORE: Store at each ISCTORE self-consistent iteration     ==
  ! == SWF:   Wavefunctions                                         ==
  ! == SRHO:  Density                                               ==
  ! == SPOT:  Potential                                             ==
  ! == TDEBIO: Debug RESTART read and write                         ==
  ! ==================================================================
  INTEGER, PARAMETER :: mstore=7 


  ! ==================================================================
  ! == STORE IN ANOTHER FILE                                        ==
  ! == MROUT: Number of parameters (parallel)                       ==
  ! == ROUT: Trajectories                                           ==
  ! == MOUT: Trajectories                                           ==
  ! == VOUT: vibrational frequencies and eigenvectors - Gaussian    ==
  ! == ACOUT: vibrational frequencies and eigenvectors - aClimax    ==
  ! == XTOUT: Trajectories in xyz-format                            ==
  ! == XGOUT: Trajectories in xyz-format of geometry optimization   ==
  ! == DCOUT: Trajectories in dcd-format                            ==
  ! == RHOOUT: DENSITY FILE                                         ==
  ! == TEBAND: ENERGYBANDS FILE                                     ==
  ! ==================================================================
  INTEGER, PARAMETER :: mrout=11 

  LOGICAL :: trajsmall
  INTEGER :: trajsmalln
  ! ==================================================================
  ! == PRINT OPTIONS                                                ==
  ! == TPRINT: True if PRINT OPTION                                 ==
  ! == IPRINT_STEP: PRINT EACH IPRINT_STEP STEP                     ==
  ! == IPRINT:                                                      ==
  ! ==================================================================
  INTEGER, PARAMETER :: iprint_max     = 33,   iprint_info    =  1,&
       iprint_eigen   =  2,   iprint_coor    =  3,&
       iprint_force   =  4,   iprint_wann    =  5,&
       iprint_ekin    =  6,   iprint_eeig    =  7,&
       iprint_eband   =  8,   iprint_entropy =  9,&
       iprint_eht     = 10,   iprint_elec1   = 11,&
       iprint_elec2   = 12,   iprint_ehep    = 13,&
       iprint_ehee    = 14,   iprint_ehii    = 15,&
       iprint_eself   = 16,   iprint_esr     = 17,&
       iprint_epseu   = 18,   iprint_enl     = 19,&
       iprint_exc     = 20,   iprint_vxc     = 21,&
       iprint_egc     = 22,   iprint_ebogo   = 23,&
       iprint_etot1   = 24,   iprint_etot2   = 25,&
       iprint_ecas    = 26,   iprint_epen    = 27,&
       iprint_etddft  = 28,   iprint_lscal   = 29,&
       iprint_eext    = 30,   iprint_ehsic   = 31,&
       iprint_evdw    = 32,   iprint_ehub    = 33
  ! ==--------------------------------------------------------------==
  LOGICAL :: twritebintrajectory
  LOGICAL :: twriteforcetrajectory

  INTEGER :: maxwriteatom

  ! ==================================================================
  ! == INTERFACE FILE TO OTHER PW CODES                             ==
  ! == INTREAD : Read Interface file                                ==
  ! == INTWRITE : Write Interface file                              ==
  ! == INTFN : Filename of Interface file                           ==
  ! ==================================================================

  CHARACTER(len=40) :: intfn

  TYPE :: cprint_t
     LOGICAL :: tprint
     INTEGER :: iprint_step
     INTEGER :: iprint(iprint_max)
     LOGICAL :: twritebintrajectory
     INTEGER :: MINWRITEATOM
     INTEGER :: MAXWRITEATOM
     LOGICAL :: TWRITEFORCETRAJECTORY
  END TYPE cprint_t
  TYPE(cprint_t) :: cprint
  TYPE :: iface1_t
     LOGICAL :: intread
     LOGICAL :: intwrite
  END TYPE iface1_t
  TYPE(iface1_t) :: iface1
  TYPE :: restart1_t
     LOGICAL :: restart
     LOGICAL :: rwf
     LOGICAL :: rco
     LOGICAL :: rvel
     LOGICAL :: rac
     LOGICAL :: rhe
     LOGICAL :: rnoe
     LOGICAL :: rnop
     LOGICAL :: rnoc
     LOGICAL :: rgeo
     LOGICAL :: rlate
     LOGICAL :: rvib
     LOGICAL :: rcell
     LOGICAL :: rpot
     LOGICAL :: rkpt
     LOGICAL :: rrho
     LOGICAL :: rocc
     LOGICAL :: reigv
     LOGICAL :: roh
     LOGICAL :: rlr
     LOGICAL :: rphes
     LOGICAL :: rlscl
     LOGICAL :: rdtl
     LOGICAL :: rcon
     LOGICAL :: rnon
     LOGICAL :: rxtp
     LOGICAL :: rprng
     LOGICAL :: rgle
  END TYPE restart1_t
  TYPE(restart1_t) :: restart1
  TYPE :: rout1_t
     LOGICAL :: rout
     LOGICAL :: mout
     LOGICAL :: vout
     LOGICAL :: acout
     LOGICAL :: xtin
     LOGICAL :: dcout
     LOGICAL :: xtout
     LOGICAL :: xgout
     LOGICAL :: rhoout
     LOGICAL :: teband
     INTEGER :: nrhoout
  END TYPE rout1_t
  TYPE(rout1_t) :: rout1
  TYPE :: store1_t
     INTEGER :: istore
     INTEGER :: isctore
     LOGICAL :: swf
     LOGICAL :: srho
     LOGICAL :: spot
     LOGICAL :: tdebio
     LOGICAL :: tdebacc
  END TYPE store1_t
  TYPE(store1_t) :: store1

END MODULE store_types
