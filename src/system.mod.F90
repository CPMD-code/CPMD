MODULE system
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  PRIVATE

  ! ==================================================================
  ! == MAXSP  : MAXIMUM NUMBER OF ATOMIC SPECIES                    ==
  ! == NHX    : MAXIMUM NUMBER OF GAUSS-HERMIT POINTS               ==
  ! ==        : MAXIMUM NUMBER OF BETA FUNCTIONS (VANDERBILT)       ==
  ! ==        : MAXIMUM NUMBER OF PROJECTORS FOR NON-LOCAL PP       ==
  ! == NACX   : SIZE OF THE ACCUMULATOR ARRAY                       ==
  ! == MAXDIS : MAX. NUMBER OF VECTORS IN DIIS                      ==
  ! == MXGDIS : MAX. NUMBER OF VECTORS IN GDIIS                     ==
  ! == NBRX   : NUMBER OF DISTINCT RADIAL BETA FUNCTIONS            ==
  ! == LX     : 2*LMAX-1 WHERE LMAX=3 IF S,P AND D FUNCTIONS        ==
  ! == NLX    : COMBINED ANGULAR MOMENTUM FOR LLi (VANDERBILT)      ==
  ! == MX     : 2*LX-1                            (VANDERBILT)      ==
  ! == LIX    : MAX. LLi                          (VANDERBILT)      ==
  ! == MIX    : 2*LIX-1                           (VANDERBILT)      ==
  ! == LMAXX  : MAXIMUM NUMBER OF ANGULAR MOMENTUM                  ==
  ! == LM1X   : LMAXX-1                                             ==
  ! ==================================================================
  INTEGER, PARAMETER, PUBLIC :: maxsp=140
  INTEGER, PARAMETER, PUBLIC :: nhx=208
  INTEGER, PARAMETER, PUBLIC :: nacx=50
  INTEGER, PARAMETER, PUBLIC :: maxdis=20
  INTEGER, PARAMETER, PUBLIC :: mxgdis=20
  INTEGER, PARAMETER, PUBLIC :: nbrx=6
  INTEGER, PARAMETER, PUBLIC :: lx=5
  INTEGER, PARAMETER, PUBLIC :: nlx=9
  INTEGER, PARAMETER, PUBLIC :: mx=2*lx-1
  INTEGER, PARAMETER, PUBLIC :: lix=3
  INTEGER, PARAMETER, PUBLIC :: mix=lix*2-1
  INTEGER, PARAMETER, PUBLIC :: lmaxx=4
  INTEGER, PARAMETER, PUBLIC :: lm1x=lmaxx-1
  ! ==================================================================
  ! == MAXSYS : Maximum values for system                           ==
  ! ==          NSX   Number of atomic species                      ==
  ! ==          NAX   Maximum number of atoms per species           ==
  ! ==          NCORX (NSX*(NSX+1))/2 (detsp.F)                     ==
  ! ==          NHXS  Maximum number for all species as NHX         ==
  ! ==          LPMAX Maximum value of angular momentum             ==
  ! ==          MMAXX Maximum number of spline points               ==
  ! ==--------------------------------------------------------------==
  TYPE, PUBLIC :: maxsys_t
     INTEGER :: nsx = HUGE(0)
     INTEGER :: nax = HUGE(0)
     INTEGER :: nsax = HUGE(0)
     INTEGER :: ncorx = HUGE(0)
     INTEGER :: nhxs = HUGE(0)
     INTEGER :: lpmax = HUGE(0)
     INTEGER :: mmaxx = HUGE(0)
  END TYPE maxsys_t
  TYPE(maxsys_t), SAVE, PUBLIC :: maxsys
  ! ==================================================================
  ! == SPAR   : Global array dimensions (parallel work)             ==
  ! ==          NHGS = sum of NHG for all processors                ==
  ! ==          NHGLS= sum of NHGL     --                           ==
  ! ==          NHGKS= sum of NHGK     --                           ==
  ! ==          NGWS = sum of NGW      --                           ==
  ! ==          NGWKS= sum of NGWKS    --                           ==
  ! ==          NGWLS= sum of NGWL     --                           ==
  ! ==          NR1S,NR2S,NR3S total mesh dimension                 ==
  ! ==================================================================
  TYPE, PUBLIC :: spar_t
     INTEGER :: nhgs = HUGE(0)
     INTEGER :: nhgls = HUGE(0)
     INTEGER :: nhgks = HUGE(0)
     INTEGER :: ngws = HUGE(0)
     INTEGER :: ngwks = HUGE(0)
     INTEGER :: ngwls = HUGE(0)
     INTEGER :: nldxs = HUGE(0)
     INTEGER :: nr1s = HUGE(0)
     INTEGER :: nr2s = HUGE(0)
     INTEGER :: nr3s = HUGE(0)
  END TYPE spar_t
  TYPE(spar_t), SAVE, PUBLIC :: spar
  ! ==--------------------------------------------------------------==
  ! == NCPW   : Number of plane waves                               ==
  ! ==          NHG  number of G-components for electronic density  ==
  ! ==          NHGL number of G-shells for electronic density      ==
  ! ==          NGW  number of G-components for wavefunctions       ==
  ! ==          NGWL number of G-shells for wavefunctions           ==
  ! ==--------------------------------------------------------------==
  TYPE, PUBLIC :: ncpw_t
     INTEGER :: nhg = HUGE(0)
     INTEGER :: nhgl = HUGE(0)
     INTEGER :: ngw = HUGE(0)
     INTEGER :: ngwl = HUGE(0)
  END TYPE ncpw_t
  TYPE(ncpw_t), SAVE, PUBLIC :: ncpw
  ! ==--------------------------------------------------------------==
  ! == NKPT   : Used with K-points options                          ==
  ! ==          NHGK = 2 * NHG (no inversion symmetry)              ==
  ! ==          NGWK = 2 * NGW (no inversion symmetry)              ==
  ! ==          NKPNT - number of k points  (in memory)             ==
  ! ==          NKPTS - number of total k points >= NKPNT           ==
  ! ==          NBLKP - number of k points block                    ==
  ! ==          NKPBL(NBLKP) - number of k points for each block    ==
  ! ==          KPBEG(NBLKP) - number - 1 of first k-pt in block    ==
  ! ==--------------------------------------------------------------==
  TYPE :: nkpt_t
     INTEGER :: nhgk = HUGE(0)
     INTEGER :: ngwk = HUGE(0)
     INTEGER :: nkpnt = HUGE(0)
     INTEGER :: nkpts = HUGE(0)
     INTEGER :: nblkp = HUGE(0)
  END TYPE nkpt_t
  TYPE(nkpt_t), SAVE, PUBLIC :: nkpt

  INTEGER, ALLOCATABLE, PUBLIC :: nkpbl(:)
  INTEGER, ALLOCATABLE, PUBLIC :: kpbeg(:)
  ! ==--------------------------------------------------------------==
  ! == PARM   : Information about supercell and mesh                ==
  ! ==          ALAT lattice parameter                              ==
  ! ==          A1(3), A2(3), A3(3) crystal vector basis            ==
  ! ==          OMEGA Supercell volume                              ==
  ! ==          TPIBA=2._real_8*PI/ALAT                                  ==
  ! ==          TPIBA2=TPIBA*TPIBA                                  ==
  ! ==          APBC(4) used by PBC routine                         ==
  ! ==                  for periodic boundary condition             ==
  ! ==          IBRAV Determines shape and symmetry of the supercell==
  ! ==              =0 Isolated system (cubic cell)                 ==
  ! ==               1 Simple cubic                                 ==
  ! ==               2 Face centred cubic                           ==
  ! ==               3 Body centred cubic                           ==
  ! ==               4 Hexagonal cell                               ==
  ! ==               5 Trigonal or rhombohedral cell                ==
  ! ==               6 Tetragonal                                   ==
  ! ==               7 Body tetragonal                              ==
  ! ==               8 Orthorhombic                                 ==
  ! ==              12 Monoclinic                                   ==
  ! ==              14 Triclinic                                    ==
  ! ==          NR1, NR2, NR3 mesh dimension for each processor     ==
  ! ==--------------------------------------------------------------==
  TYPE, PUBLIC :: parm_t
     REAL(real_8) :: alat = 0.0_real_8
     REAL(real_8) :: a1(3) = 0.0_real_8
     REAL(real_8) :: a2(3) = 0.0_real_8
     REAL(real_8) :: a3(3) = 0.0_real_8
     REAL(real_8) :: omega = 0.0_real_8
     REAL(real_8) :: tpiba = 0.0_real_8
     REAL(real_8) :: tpiba2 = 0.0_real_8
     REAL(real_8) :: apbc(4) = 0.0_real_8
     INTEGER :: ibrav = HUGE(0)
     INTEGER :: nr1 = HUGE(0)
     INTEGER :: nr2 = HUGE(0)
     INTEGER :: nr3 = HUGE(0)
  END TYPE parm_t
  TYPE(parm_t), SAVE, PUBLIC :: parm
  ! ==--------------------------------------------------------------==
  ! == FPAR   : Leading dimensions for some arrays, especially FFT  ==
  ! ==          KR1, KR2, KR3 FFT mesh dimension for each proc.     ==
  ! ==          KR1S, KR2S, KR3S FFT mesh for all proc.             ==
  ! ==          KRX Used in FFTPRP routine with groups              ==
  ! ==          KRY, KRZ UNUSED                                     ==
  ! ==          NNR1 Number of real components for density          ==
  ! ==               for each processor                             ==
  ! ==--------------------------------------------------------------==
  TYPE :: fpar_t
     INTEGER :: kr1 = HUGE(0)
     INTEGER :: kr2 = HUGE(0)
     INTEGER :: kr3 = HUGE(0)
     INTEGER :: kr1s = HUGE(0)
     INTEGER :: kr2s = HUGE(0)
     INTEGER :: kr3s = HUGE(0)
     INTEGER :: krx = HUGE(0)
     INTEGER :: kry = HUGE(0)
     INTEGER :: krz = HUGE(0)
     INTEGER :: nnr1 = HUGE(0)
  END TYPE fpar_t
  TYPE(fpar_t), SAVE, PUBLIC :: fpar
  ! ==--------------------------------------------------------------==
  ! == ACCU   : Accumulator quantities                              ==
  ! ==          NACC number of quantities calculated per self-c iter==
  ! ==          ACC(1:NACC) values of calculated quantities         ==
  ! ==--------------------------------------------------------------==
  INTEGER, SAVE, PUBLIC :: nacc = HUGE(0)
  REAL(real_8), SAVE, PUBLIC :: acc(nacx)  = 0.0_real_8 !vw never set to zero 
  ! ==================================================================
  ! == ALL THE INFO FOR PARALLEL                                    ==
  ! ==================================================================
  ! ==--------------------------------------------------------------==
  ! == MAXCPU  : Maximum number of CPU
  ! ==--------------------------------------------------------------==
  !INTEGER, PARAMETER :: maxcpu =2**16
  ! ==--------------------------------------------------------------==
  ! == NRXPL   : NR1 for each processor                             ==
  ! == NRZPL   : Number of z-planes for each processor              ==
  ! == SPARM   : Information about dimension of distributed data    ==
  ! == NST12   : index of the first and last distributed states     ==
  ! ==         : per proc                                           ==
  ! == PGROUP  : Pointer position --> PE number                     ==
  ! == NLINK   : Pointer PE number --> position - 1                 ==
  ! == NWA12   : index of distributed states (dynamical)            ==
  ! ==--------------------------------------------------------------==
  TYPE, PUBLIC :: parap_t
     INTEGER, ALLOCATABLE, DIMENSION(:,:) :: nrxpl !(0:maxcpu,2)
     INTEGER, ALLOCATABLE, DIMENSION(:,:) :: nrzpl !(0:maxcpu,2)
     INTEGER, ALLOCATABLE, DIMENSION(:,:) :: sparm !(9,0:maxcpu)
     INTEGER, ALLOCATABLE, DIMENSION(:,:) :: nst12 !(0:maxcpu,2)
     INTEGER, ALLOCATABLE, DIMENSION(:)   :: pgroup!(0:maxcpu)
     INTEGER, ALLOCATABLE, DIMENSION(:)   :: nlink !(0:maxcpu)
  END TYPE parap_t
  TYPE(parap_t), SAVE, PUBLIC :: parap

  TYPE, PUBLIC :: paraw_t
     INTEGER, ALLOCATABLE, DIMENSION(:,:) :: nwa12 !(0:maxcpu,2)
  END TYPE paraw_t
  TYPE(paraw_t), SAVE, PUBLIC :: paraw
  ! ==================================================================
  ! == MAPPING G VECTORS-->PROCESSES                                ==
  ! ==================================================================
  ! == SUBDIVIDE THE PROCESSORS IN GROUPS (see GROUPS routine)      ==
  ! == MAPGP(1:NHG) : Use for distribution of wavefunction          ==
  ! == MAXGRP=MAXCPU : Maximum taskgroups (for distributed FFT)     ==
  ! == NOGRP  : Number of groups                                    ==
  ! == NPGRP  : Number of processors per groups                     ==
  ! == MEOGRP : Index of the group for the processor                ==
  ! == MEPGRP : Index of the processor in the group                 ==
  ! == MPEN   : Dimensions used for groups (see FFTPRP routine)     ==
  ! == MPENM  : ''                  ''          ''                  ==
  ! == NOLIST(1:MAXGRP) : List of processors in the orbital group   ==
  ! == NPLIST(1:MAXGRP) : List of processors in the plane wave group==
  ! ==--------------------------------------------------------------==
  INTEGER, ALLOCATABLE, SAVE, PUBLIC :: mapgp(:)

  INTEGER, PARAMETER, PUBLIC :: maxgrp=1 !maxcpu 
  TYPE, PUBLIC :: group_t
     INTEGER :: nogrp = HUGE(0)
     INTEGER :: npgrp = HUGE(0)
     INTEGER :: meogrp = HUGE(0)
     INTEGER :: mepgrp = HUGE(0)
     INTEGER :: mpen = HUGE(0)
     INTEGER :: mpenm = HUGE(0)
     INTEGER :: nolist(maxgrp) = HUGE(0)
     INTEGER :: nplist(maxgrp) = HUGE(0)
  END TYPE group_t
  TYPE(group_t), SAVE, PUBLIC :: group

  ! ==================================================================
  ! == MAPPING ATOMS-->PROCESSES                                    ==
  ! ==================================================================
  ! == NATPE           =IPEPT(2,MEPOS)-IPEPT(1,MEPOS)+1             ==
  ! ==                 Number of atoms per processes                ==
  ! == NORBPE          =NST12(MEPOS,2)-NST12(MEPOS,1)+1             ==
  ! ==                 Number of orbitals per processes             ==
  ! == IPEPT(2,0:MAXCPU) Indexes of atoms per processes             ==
  ! == IATPT(2,NAT)                                                 ==
  ! == IATPE(NAT)                                                   ==
  ! ==--------------------------------------------------------------==
  INTEGER, SAVE, PUBLIC :: natpe = HUGE(0),norbpe = HUGE(0)
  INTEGER, ALLOCATABLE, DIMENSION(:,:), SAVE, PUBLIC :: ipept !(2,0:maxcpu)
  INTEGER, ALLOCATABLE, SAVE, PUBLIC :: iatpt(:,:)
  INTEGER, ALLOCATABLE, SAVE, PUBLIC :: iatpe(:)
  ! ==================================================================
  ! == A LOT OF ITEMS IN THE COMMON BLOCKS IS NEEDED FOR PARALLEL   ==
  ! ==================================================================
  ! == CNTL1  : Use for the broadcast (dummy variable)              ==
  ! == MD     : Molecular dynamics option                           ==
  ! == GEOPT  : Geometry optimization                               ==
  ! == WFOPT  : Wavefunction optimization                           ==
  ! == TLOWD  : Lowdin orthogonalization instead of Gram-Schmidt    ==
  ! == TLOWDMAT : Lowdin transformation matrix output               ==
  ! == TRANE  : Randomize wavefunctions                             ==
  ! == TRANP  : Randomize ionic coordinates                         ==
  ! == TC     : Temperature control for electrons                   ==
  ! == TCP    : Temperature control for ions                        ==
  ! == TBERP  : Berendsen for ions                                  == 
  ! == TBERE  : Berendsen for electrons                             == 
  ! == TBERC  : Berendsen for cell                                  == 
  ! == ANNEI  : Simulated annealing for ions                        ==
  ! == ANNEE  : Simulated annealing for electrons                   ==
  ! == ANNEC  : Simulated annealing for cell                        ==
  ! == DAMPI  : Damped dynamics for ions                            ==
  ! == DAMPE  : Damped dynamics for electron                        ==
  ! == DAMPC  : Damped dynamics for cell                            ==
  ! == QUENCHP: Quench ionic velocities                             ==
  ! == QUENCHE: Quench electronic velocities (Car-Parrinello)       ==
  ! == QUENCHB: Quench system to the Born-Oppenheimer surface       ==
  ! == CNSTMD : Iterative orthogonalization (Car-Parrinello)        ==
  ! == TSDE   : Steepest descent (electronic coordinates)           ==
  ! == TSDE_P : Steepest descent in DFTP
  ! == TSDP   : Steepest descent (ionic coordinates)                ==
  ! == DIIS   : Wavefunction optimization by DIIS                   ==
  ! == PREC   : Preconditioning for Steepest descent (electronic)   ==
  ! == PREC_P : Preconditioning for DFTP
  ! == TSCALE : Use Scaled ionic coordinates                        ==
  ! == BOHR   : Use atomic units for ionic coordinates (bohr)       ==
  ! == PCG    : Wavefunction optimization by Precond. Conj. Gradient==
  ! == PCG_P  : Wavefun. optimization by Precond.Conj.Gradient DFTP ==
  ! == GDIIS  : Geometry optimization by GDIIS                      ==
  ! == RFO    : Geometry optimization by Rational Function Approx.  ==
  ! == BFGS   : Geometry optimization by quasi-Newton update        ==
  ! == LBFGS  : Geometry optimization low-memory BFGS implicit upd. ==
  ! == PRFO   : Use partitioned rational function optimizer (P-RFO) ==
  ! == TSDIN  : Second derivatives from input                       ==
  ! == TGC    : Gradient correction                                 ==
  ! == TGCX   : Gradient correction for exchange part               ==
  ! == TGCC   : Gradient correction for correlation part            ==
  ! == TTAU   : Tau dependent functional                            ==
  ! == THYBRID: Hybrid functionals                                  ==
  ! == TSMOOTH: Smoothing of the density                            ==
  ! == CALDIP : Calculate Dipole dynamics                           ==
  ! == TNOSEE : Nose-Hoover thermostats for electrons dynamics      ==
  ! == TNOSEP : Nose-Hoover thermostats for ions dynamics           ==
  ! == TNOSES : One thermostat chain per species                    ==
  ! == KSENER : Calculate Kohn-Sham energies for the final potential==
  ! == WCOMP  : Wavefunction in compressed form                     ==
  ! == NONORT : Use nonorthogonal orbital                           ==
  ! == BIGMEM : Big memory option                                   ==
  ! == THARM  : Harmonic reference system integration               ==
  ! == TMASS  : Scaled electron masses                              ==
  ! == KRWFN  : Store real space representation of wavefunctions    ==
  ! == PCGMIN : Perform quadratic line search                       ==
  ! == PCGMIN_P : Perform quadratic line search  in DFPT            ==
  ! == PROPER : Properties calculation                              ==
  ! == VIBRAT : Vibrational analysis                                ==
  ! == TSDAN  : Linear Response for second derivatives              ==
  ! == TLSD   : Local Spin Density Approximation                    ==
  ! == TPRES  : Stress tensor calculation                           ==
  ! == TVIRIAL : Use virial estimator of stress tensor              ==
  ! == SIMUL  : Combined geometry/wavefunction scheme               ==
  ! == TPRCP  : Parrinello-Rahman-Car-Parrinello                    ==
  ! == TNOSEC : Nose-Hoover thermostats for cell dynamics           ==
  ! == QUENCHC: Quench electronic velocities (Parrinello-Rahman)    ==
  ! == TSDC   : Cell dynamics by Steepest descent                   ==
  ! == TRANC  : Randomize cell                                      ==
  ! == TMEMCHK: Memory checking                                     ==
  ! == TPENN  : NPT dynamics                                        ==
  ! == TFINT  : Free Energy Functional (Alavi''s method)            ==
  ! == TDIAG  : Diagonalisation scheme (Lanczos or Davidson)        ==
  ! == TDIAGOPT : Optimize then diagonalize                         ==
  ! == TDAVI  : Davidson diagonalisation                            ==
  ! == TLANC  : Lanczos diagonalisation                             ==
  ! == TFRHO_UPW: Update Wavefunction with Rho fixed (DIIS + mixing)==
  ! == TFRSBLK: New Block Lanczos diagonalisation                   ==
  ! == TDIPD  : Dipole Dynamics                                     ==
  ! == TINTER : Interface to classical MD program                   ==
  ! == TCC    : Temperature Control for the cell                    ==
  ! == TEPOT  : Electrostatic Potential                             ==
  ! == TEXPOT : External Potential                                  ==
  ! == TEXADD : Add energy and forces of external potential         ==
  ! == TPATH  : Path Integrals (PIMD or PATH for saddle point)      ==
  ! == TPIMD  : Path Integral Molecular Dynamics                    ==
  ! == TPMIN  : Path Minimisation (Saddle point determination)      ==
  ! == TRESCALE : Re-adjust ionic velocities only after restart to  ==
  ! ==          desired temperature TEMPW                           ==
  ! == TSYMRHO: SYMMTRIZE DENSITY                                   ==
  ! == THARD  : Calculate orbital hardness matrix                   ==
  ! == TSPEC  : Calculate spectra with TDDFT                        ==
  ! == TSHOP  : Surface hopping dynamics S0/S1                      ==
  ! == TFUSI  : Fuse two Restart files for S0/S1                    ==
  ! == TMERGE : Merge two RESTART files and shift the wavefunction  ==
  ! == TSEP   : Separate one Restart into two files for S0/S1       ==
  ! == TREVERS: Invert ionic and electronic velocities              ==
  ! == TQMMM  : Use the QM/MM code (full Bio code)                  ==
  ! == TQMMECH: Use the QM/MM code (simple coupling only)           ==
  ! == TDEBFOR: Run option for force debugging                      ==
  ! == TDDFT  : Use TDDFT                                           ==
  ! == TNOGEOCHECK : Dont check geometry for close atoms            ==
  ! == TSSEL  : Print structure only for selected atoms             ==
  ! == TWSERIAL : Calculate Wannier centers with serial code        ==
  ! == TFIELD : Apply an external static electric field             ==
  ! == TPOTENTIAL : Indicates a potential only functional           ==
  ! == TR4A2A : Use single precision in ALLTOALL communication      ==
  ! == TEXTRAP: Use wavefunction extrapolation in BOMD              ==
  ! == TSTRXTP: Store wavefunction history in restart file          ==
  ! == TASPC  : Use the ASPC interpolator for extrapolation         ==
  ! == TDMAL  : Use distributed linear algebra routines             ==
  ! == TORTHO_NEW: Use new ortho for wavefunctions                  ==
  ! ==             C.Bekas and A.Curioni Comp. Phys. Comm. (2010)   ==
  ! == TPCGFI : Use PCG MINIMIZE for first Wfopt (EGO/GMX interface)==
  ! == TPCGFIC: Use PCG MINIMIZE for first Wfopt (CDFT)             ==
  ! == TCART  : Use cartesian taskgroups                            ==
  ! == GSHELLOUT : Write file with Gshell information               ==
  ! == TSHOCK : Shock wave simulation                               ==
  ! == TSIC   : SS-SIC self-intercation correction                  ==
  ! == TNACVS : TDDFT based non adiabatic coupling vectors          ==
  ! == CDFT   : Use charge constraints                              ==
  ! == CDFT_NEWT   : Use Newton Optimiser for CDFT                  ==
  ! == CDFT_DEKK   : Use Dekker Optimiser for CDFT                  ==
  ! == CDFT_WEIGHT : Write out a slice of the weight                ==
  ! == CDFT_WF     : Write out the full weight                      ==
  ! == TSYSCOMB    : Combine two reference systems                  ==
  ! == TSCOMBM     : Mirror the resulting WF around COM             ==
  ! == TSCORTHO    : Orthogonalise the resulting WF                 ==
  ! == TKSHAM      : Write out Kohn Sham Hamiltonian                ==
  ! == TIPAO       : Initialise with the primitive AO method        ==
  ! == TLPOT  : Store local potential                               ==
  ! == TSOC   : do spin-orbit coupling calculation                  ==
  ! == FMATCH : Use to run force matching in QM/MM                  ==
  ! == IS_IN_STREAM : use stream input                              ==
  ! == IS_OUT_STREAM : use stream output                            ==
  ! == USE_MPI_IO : use mpi for the parallel io                     ==
  ! == USE_MTS: use multiple time step algorithm in MD              ==
  ! ==--------------------------------------------------------------==
  TYPE, PUBLIC :: cntl_t
     LOGICAL :: md = .FALSE.
     LOGICAL :: tmdbo = .FALSE.
     LOGICAL :: tmdfile = .FALSE.
     LOGICAL :: geopt = .FALSE.
     LOGICAL :: wfopt = .FALSE.
     LOGICAL :: tlowd = .FALSE.
     LOGICAL :: tlowdmat = .FALSE.
     LOGICAL :: trane = .FALSE.
     LOGICAL :: tranp = .FALSE.
     LOGICAL :: thead = .FALSE.
     LOGICAL :: tc = .FALSE.
     LOGICAL :: tcp = .FALSE.
     LOGICAL :: annei = .FALSE.
     LOGICAL :: annee = .FALSE.
     LOGICAL :: annec = .FALSE.
     LOGICAL :: dampi = .FALSE.
     LOGICAL :: dampe = .FALSE.
     LOGICAL :: dampc = .FALSE.
     LOGICAL :: quenchp = .FALSE.
     LOGICAL :: quenche = .FALSE.
     LOGICAL :: quenchb = .FALSE.
     LOGICAL :: cnstmd = .FALSE.
     LOGICAL :: stopot = .FALSE.
     LOGICAL :: tchngw = .FALSE.
     LOGICAL :: tsde = .FALSE.
     LOGICAL :: tsdp = .FALSE.
     LOGICAL :: diis  = .FALSE.
     LOGICAL :: prec = .FALSE.
     LOGICAL :: tscale = .FALSE.
     LOGICAL :: bohr = .FALSE.
     LOGICAL :: tchngr = .FALSE.
     LOGICAL :: pcg = .FALSE.
     LOGICAL :: gdiis = .FALSE.
     LOGICAL :: rfo = .FALSE.
     LOGICAL :: bfgs = .FALSE.
     LOGICAL :: tgc = .FALSE.
     LOGICAL :: ttau = .FALSE.
     LOGICAL :: tgcx  = .FALSE.
     LOGICAL :: tgcc = .FALSE.
     LOGICAL :: tsmooth = .FALSE.
     LOGICAL :: tberp = .FALSE.
     LOGICAL :: tbere = .FALSE.
     LOGICAL :: tberc = .FALSE.
     LOGICAL :: caldip = .FALSE.
     LOGICAL :: tnosee = .FALSE.
     LOGICAL :: tnosep = .FALSE.
     LOGICAL :: tnoses  = .FALSE.
     LOGICAL :: timing = .FALSE.
     LOGICAL :: ksener = .FALSE.
     LOGICAL :: orbrot = .FALSE.
     LOGICAL :: rcomp = .FALSE.
     LOGICAL :: wcomp = .FALSE.
     LOGICAL :: nonort = .FALSE.
     LOGICAL :: bigmem = .FALSE.
     LOGICAL :: gaudyn = .FALSE.
     LOGICAL :: tharm  = .FALSE.
     LOGICAL :: tmass = .FALSE.
     LOGICAL :: fcnstr = .FALSE.
     LOGICAL :: krwfn = .FALSE.
     LOGICAL :: pcgmin = .FALSE.
     LOGICAL :: posver = .FALSE.
     LOGICAL :: proper = .FALSE.
     LOGICAL :: vibrat = .FALSE.
     LOGICAL :: initcm  = .FALSE.
     LOGICAL :: finalcm = .FALSE.
     LOGICAL :: tlsd = .FALSE.
     LOGICAL :: tpres = .FALSE.
     LOGICAL :: tvirial = .FALSE.
     LOGICAL :: simul = .FALSE.
     LOGICAL :: tprcp = .FALSE.
     LOGICAL :: tnosec = .FALSE.
     LOGICAL :: quenchc = .FALSE.
     LOGICAL :: tsdc = .FALSE.
     LOGICAL :: tranc = .FALSE.
     LOGICAL :: tmemchk = .FALSE.
     LOGICAL :: tpenn = .FALSE.
     LOGICAL :: tfint = .FALSE.
     LOGICAL :: tdiag = .FALSE.
     LOGICAL :: tdiagopt = .FALSE.
     LOGICAL :: tdavi = .FALSE.
     LOGICAL :: tlanc = .FALSE.
     LOGICAL :: tfrsblk = .FALSE.
     LOGICAL :: tdipd  = .FALSE.
     LOGICAL :: tinter = .FALSE.
     LOGICAL :: tcc = .FALSE.
     LOGICAL :: tepot = .FALSE.
     LOGICAL :: texpot = .FALSE.
     LOGICAL :: texadd = .FALSE.
     LOGICAL :: tpath = .FALSE.
     LOGICAL :: tpimd = .FALSE.
     LOGICAL :: tpmin = .FALSE.
     LOGICAL :: trescale  = .FALSE.
     LOGICAL :: tsampl = .FALSE.
     LOGICAL :: tprec = .FALSE.
     LOGICAL :: tfdist = .FALSE.
     LOGICAL :: tsymrho = .FALSE.
     LOGICAL :: tresponse = .FALSE.
     LOGICAL :: tsdan = .FALSE.
     LOGICAL :: thard = .FALSE.
     LOGICAL :: tspec  = .FALSE.
     LOGICAL :: tshop = .FALSE.
     LOGICAL :: tfusi = .FALSE.
     LOGICAL :: tmerge = .FALSE.
     LOGICAL :: trevers = .FALSE.
     LOGICAL :: tqmmm = .FALSE.
     LOGICAL :: tqmmech = .FALSE.
     LOGICAL :: tinr = .FALSE.
     LOGICAL :: thybrid  = .FALSE.
     LOGICAL :: tdebfor = .FALSE.
     LOGICAL :: tddft = .FALSE.
     LOGICAL :: lbfgs = .FALSE.
     LOGICAL :: prfo = .FALSE.
     LOGICAL :: tsdin = .FALSE.
     LOGICAL :: tnogeocheck = .FALSE.
     LOGICAL :: tssel = .FALSE.
     LOGICAL :: tfrho_upw  = .FALSE.
     LOGICAL :: twserial = .FALSE.
     LOGICAL :: bsymm = .FALSE.
     LOGICAL :: tfield = .FALSE.
     LOGICAL :: tpotential = .FALSE.
     LOGICAL :: tr4a2a = .FALSE.
     LOGICAL :: textrap = .FALSE.
     LOGICAL :: tdmal  = .FALSE.
     LOGICAL :: tortho_new = .FALSE.
     LOGICAL :: tsep = .FALSE.
     LOGICAL :: tpcgfi = .FALSE.
     LOGICAL :: tpcgfic = .FALSE.
     LOGICAL :: tcart = .FALSE.
     LOGICAL :: gshellout = .FALSE.
     LOGICAL :: tshock = .FALSE.
     LOGICAL :: tsic  = .FALSE.
     LOGICAL :: tstrxtp = .FALSE.
     LOGICAL :: taspc = .FALSE.
     LOGICAL :: tnacvs = .FALSE.
     LOGICAL :: cdft = .FALSE.
     LOGICAL :: cdft_newt = .FALSE.
     LOGICAL :: cdft_dekk  = .FALSE.
     LOGICAL :: cdft_weight = .FALSE.
     LOGICAL :: cdft_wf = .FALSE.
     LOGICAL :: tsyscomb = .FALSE.
     LOGICAL :: tscombm = .FALSE.
     LOGICAL :: tscortho = .FALSE.
     LOGICAL :: tipao = .FALSE.
     LOGICAL :: tksham  = .FALSE.
     LOGICAL :: fmatch = .FALSE.
     LOGICAL :: cmplx_wf = .FALSE.
     LOGICAL :: start_real = .FALSE.
     LOGICAL :: tmdeh = .FALSE.
     LOGICAL :: cheby = .FALSE.
     LOGICAL :: cayley = .FALSE.
     LOGICAL :: ruku = .FALSE.
     LOGICAL :: tgaugep  = .FALSE.
     LOGICAL :: tgaugef = .FALSE.
     LOGICAL :: tpdist = .FALSE.
     LOGICAL :: tpspec = .FALSE.
     LOGICAL :: is_in_stream = .FALSE.
     LOGICAL :: is_out_stream = .FALSE.
     LOGICAL :: use_mpi_io  = .FALSE.
     LOGICAL :: tsoc = .FALSE.
     LOGICAL :: tnabdy = .FALSE.
     LOGICAL :: use_xc_driver = .FALSE.
     LOGICAL :: div_analytical = .FALSE.
     LOGICAL :: thubb  = .FALSE.
     LOGICAL :: use_mts = .FALSE.
     LOGICAL :: use_scaled_hfx = .FALSE.
  END TYPE cntl_t
  TYPE(cntl_t), SAVE, PUBLIC :: cntl
  ! ==================================================================
  ! == CNTI1  : Use for the broadcast (dummy variable)              ==
  ! == ISOCS,JSOCT : states coupled in SOC calculation              ==
  ! == NOMORE : Maximum number of steps                             ==
  ! == NOMORE_ITER: Maximum number of iteration for SC              ==
  ! ==              Used only with TDIAG or TMDBO                   ==
  ! == MAXIT  : Parameters for the RATTLE step                      ==
  ! == MDIIS  : Maximum number of vectors retained for DIIS (Wave.) ==
  ! == INSYS  : Logic unit for properties (by default=5)            ==
  ! == MGDIIS : Size of GDIIS matrix (geometry optimization)        ==
  ! == NPARA  : Initial Hessian (Unit, Disco or Schlegel)           ==
  ! == NCHP   : Number of Nose thermostats (ions)                   ==
  ! == NCHB   : Number of Nose thermostats (cell)                   ==
  ! == NCHS   : Number of Nose thermostats (electrons)              ==
  ! == NCALLS0: Number of Yoshida-Suzuki steps (Nose)               ==
  ! == NIT0   : Number of integration cycles (Nose)                 ==
  ! == NKSSTA : Number of Kohn-Sham eigenvalues                     ==
  ! == IMOVIE : Write movie file every imovie steps                 ==
  ! == IPROJ  : Electronic gradient projection to be used           ==
  ! == INWFUN : Wavefunction initialization (1=random, 2=atomic)    ==
  ! == NRESTF : Number of different restart files                   ==
  ! == NTRANS : Use in Lanczos scheme                               ==
  ! == NTRAJ  : Trajectories are saved on file every NTRAJ          ==
  ! == NPRES  : Stress tensor calculation each NPRES                ==
  ! == NSPLP  : Number of points for spline                         ==
  ! == NSORDER: Saddle point order for geo. opt. by rational func.  ==
  ! == NVIB   : Write VIB[1,2].log file for vibrational analysis    ==
  ! == NGXYZ  : Write xyz file every NGXYZ steps for geo. opt.      ==
  ! == RCOMPB : Give index of read wavefunction compression scheme  ==
  ! == WCOMPB : Give index of written wave. compression scheme      ==
  ! == NDAVV  : Davidson parameter (n. of iteratation for one diag. ==
  ! == N_FRIES: Lanczos para. (n. of iteratation for one diag.)     ==
  ! == NKRY_MAX: Lanczos para.(Krylov subspace dimension)           ==
  ! == NKRY_BLOCK: Lanczos para.(Block dimension)                   ==
  ! == NPDIP:   Store dipole moments each NPDIP steps               ==
  ! == IFTYPE : Interface to classical MD program (EGO)             ==
  ! == ICMET  : Interface to classical MD program                   ==
  ! == NRESET : WF steps between ODIIS resets on poor progress      ==
  ! == NSTCNV : WF steps until WF convergence criteria are relaxed  ==
  ! == NPERID : WF steps until WF is diagonalized                   ==
  ! == NRPERI : DIIS resets until periodic WF diagonalization       ==
  ! == MDIIS_FR: Maximum # of vectors retained for DIIS (rho fix)   ==
  ! == MAXLDIIS: Maximum # of cycles at fixed rho                   ==
  ! == MINLDIIS: Minimum # of cycles at fixed rho                   ==
  ! == NSKIP: # of configurations to skip from TRAJSAVED            ==
  ! == NSAMPLE: # of conf. to sample TRAJSAVED                      ==
  ! == MEXTRA : Extrapolation order for BOMD                        ==
  ! == NASPC  : Number of corrector steps when using ASPC           ==
  ! == NSTBLK : CPU block size for dist. linalg                     ==
  ! == IPRNG   : seed for pseudo-random number generator            ==
  ! ==--------------------------------------------------------------==
  TYPE, PUBLIC :: cnti_t
     INTEGER :: nomore = HUGE(0)
     INTEGER :: nomore_iter = HUGE(0)
     INTEGER :: maxit = HUGE(0)
     INTEGER :: nps = HUGE(0)
     INTEGER :: mdiis = HUGE(0)
     INTEGER :: insys = HUGE(0)
     INTEGER :: mgdiis = HUGE(0)
     INTEGER :: npara = HUGE(0)
     INTEGER :: nchp = HUGE(0)
     INTEGER :: nchb = HUGE(0)
     INTEGER :: nchs = HUGE(0)
     INTEGER :: ncalls0 = HUGE(0)
     INTEGER :: nit0 = HUGE(0)
     INTEGER :: kssta = HUGE(0) !vw not initialized at all
     INTEGER :: nkssta = HUGE(0)
     INTEGER :: imovie = HUGE(0)
     INTEGER :: iproj = HUGE(0)
     INTEGER :: ivdbrd = HUGE(0)
     INTEGER :: ivdbwr = HUGE(0)
     INTEGER :: ntgaus = HUGE(0)
     INTEGER :: inwfun = HUGE(0)
     INTEGER :: nrestf = HUGE(0)
     INTEGER :: ntrans = HUGE(0)
     INTEGER :: ntraj = HUGE(0)
     INTEGER :: npres = HUGE(0)
     INTEGER :: nsplp = HUGE(0)
     INTEGER :: nsorder = HUGE(0)
     INTEGER :: nvib = HUGE(0)
     INTEGER :: ngxyz = HUGE(0)
     INTEGER :: rcompb = HUGE(0)
     INTEGER :: wcompb = HUGE(0)
     INTEGER :: ndavv = HUGE(0)
     INTEGER :: n_fries = HUGE(0)
     INTEGER :: nkry_max = HUGE(0)
     INTEGER :: nkry_block = HUGE(0)
     INTEGER :: npdip = HUGE(0)
     INTEGER :: iftype = HUGE(0)
     INTEGER :: icmet = HUGE(0)
     INTEGER :: nprec = HUGE(0) !vw not initialized at all
     INTEGER :: nreset = HUGE(0)
     INTEGER :: nstcnv = HUGE(0)
     INTEGER :: mdiis_fr = HUGE(0)
     INTEGER :: maxldiis = HUGE(0)
     INTEGER :: minldiis = HUGE(0)
     INTEGER :: nperid = HUGE(0) !vw not initialized at all
     INTEGER :: nrperi = HUGE(0) !vw not initialized at all
     INTEGER :: nskip = HUGE(0)
     INTEGER :: nsample = HUGE(0)
     INTEGER :: mextra = HUGE(0)
     INTEGER :: naspc = HUGE(0)
     INTEGER :: nstblk = HUGE(0)
     INTEGER :: iprng = HUGE(0)
     INTEGER :: lfit = HUGE(0) !vw not initialized at all
     INTEGER :: isocs = HUGE(0) !vw not initialized at all
     INTEGER :: jsoct = HUGE(0) !vw not initialized at all
     INTEGER :: disortho_bsize = HUGE(0)
  END TYPE cnti_t
  TYPE(cnti_t), SAVE, PUBLIC :: cnti
  ! ==================================================================
  ! == CNTR1  : Use for the broadcast (dummy variable)              ==
  ! == DELT_ELEC : Time step for electrons                          ==
  ! == DELT_IONS : Time step for ions                               ==
  ! == EMASS  : Electronic mass (Car-Parrinello)                    ==
  ! == EPSOG  : Parameters for the RATTLE step                      ==
  ! == EKINW  : Electron dynamics with rescaling of velocities      ==
  ! ==          Average kinetic energy                              ==
  ! == TOLL   : Tolerance (electron dynamics)                       ==
  ! == AMPRE  : Randomize wavefunction parameter                    ==
  ! == TEMPW  : Ion dynamics with rescaling of velocities           ==
  ! ==          desired temperature                                 ==
  ! == TOLP   : Tolerance (ion dynamics)                            ==
  ! == TAUBP  : Berendsen time constant for ions                    ==
  ! == TAUBE  : Berendsen time constant for electrons               ==
  ! == TAUBC  : Berendsen time constant for cell                    ==
  ! == TRAMPT : Temperature ramping target value                    ==
  ! == TRAMPR : Temperature ramping rate (in kelvin / a.u.)         ==
  ! == ANNERI : Simulated annealing parameter (ions)                ==
  ! == ANNERE : Simulated annealing parameter (electrons)           ==
  ! == ANNERC : Simulated annealing parameter (cell)                ==
  ! == DAMPGI : Damped dynamics friction parameter (ions)           ==
  ! == DAMPGE : Damped dynamics friction parameter (electrons)      ==
  ! == DAMPGC : Damped dyanmics friction parameter (cell)           ==
  ! == HTHRS  : Hamiltonian cutoff                                  ==
  ! == HTHRS  : Hamiltonian cutoff in DFTP                          ==
  ! == TOLOG  : Convergence orbital tolerance                       ==
  ! == TOLNG  : Convergence geometry tolerance                      ==
  ! == TOLAD  : Adaptive convergence geometry tolerance: ratio      ==
  ! == TOLENE : Adaptive convergence energy step tolerance: ratio   ==
  ! == TOLFOR : Relaxed convergence criteria for force calculation  ==
  ! == TOLINI : Initial convergence criterion for orbitals          ==
  ! == TOLDETOT : Convergence orbital criteria based on total energy==
  ! == ECUT   : Cutoff energy                                       ==
  ! == SDELTA : Smooth density parameter                            ==
  ! == SMF    : ''      ''      ''                                  ==
  ! == GCEPS  : Gradient correction cutoff                          ==
  ! == WNOSE0 : Characteristic frequency (electrons -- Nose)        ==
  ! == WNOSP0 : Characteristic frequency (ions -- Nose)             ==
  ! == EPSDAV : Davidson parameter                                  ==
  ! == AMPRP  : Randomize ionic positions parameter                 ==
  ! == FDIFF  : Finite difference step                              ==
  ! == CMASS  : Fictitious MD cell mass                             ==
  ! == AMPRC  : Randomize initial cell parameters, amplitude        ==
  ! == TEMPC  : Target cell temperature(Kelvin -- Nose)             ==
  ! == WNOSC0 : Characteristic frequency (cell -- Nose)             ==
  ! == B2LIMIT: Lanczos parameter                                   ==
  ! == EKINHR : Cell dynamics with rescaling of velocities          ==
  ! ==          Average kinetic energy                              ==
  ! == TOLC   : Tolerance (cell dynamics)                           ==
  ! == TOLKINC: Tolerance of elec. kinetic energy for cell dynamics ==
  ! == TOLCG  : Convergence cell tolerance                          ==
  ! == NEDOF0 : Scaling for elec. DOF (Nose)                        ==
  ! == DSHIFT : Shift of potential in Davidson diagonalisation      ==
  ! == TOLRHOFIX: Convergence tolerance when diis with fix rho      ==
  ! == MEMSIZE: Memory reserved for real space wavefunction         ==
  ! == ASIC   : alpha parameter for Hartree SIC correction          ==
  ! == BSIC   : beta parameter for XC SIC correction                ==
  ! == NOSPT0 : Temperature at which Nose velocities are initialized==
  ! ==--------------------------------------------------------------==
  TYPE, PUBLIC :: cntr_t
     REAL(real_8) :: delt_elec = HUGE(0.0_real_8)
     REAL(real_8) :: delt_ions = HUGE(0.0_real_8)
     REAL(real_8) :: emass = HUGE(0.0_real_8)
     REAL(real_8) :: trampt = HUGE(0.0_real_8)
     REAL(real_8) :: trampr = HUGE(0.0_real_8)
     REAL(real_8) :: epsog = HUGE(0.0_real_8)
     REAL(real_8) :: ekinw = HUGE(0.0_real_8)
     REAL(real_8) :: toll = HUGE(0.0_real_8)
     REAL(real_8) :: ampre = HUGE(0.0_real_8)
     REAL(real_8) :: tempw = HUGE(0.0_real_8)
     REAL(real_8) :: tolp = HUGE(0.0_real_8)
     REAL(real_8) :: taubp = HUGE(0.0_real_8)
     REAL(real_8) :: taube = HUGE(0.0_real_8)
     REAL(real_8) :: taubc = HUGE(0.0_real_8)
     REAL(real_8) :: anneri = HUGE(0.0_real_8)
     REAL(real_8) :: hthrs = HUGE(0.0_real_8)
     REAL(real_8) :: tolog = HUGE(0.0_real_8)
     REAL(real_8) :: tolng = HUGE(0.0_real_8)
     REAL(real_8) :: ecut = -HUGE(0.0_real_8)
     REAL(real_8) :: sdelta = HUGE(0.0_real_8) !vw not initialized at all
     REAL(real_8) :: smf = HUGE(0.0_real_8)
     REAL(real_8) :: gceps = HUGE(0.0_real_8) !vw not initialized at all
     REAL(real_8) :: wnose0 = HUGE(0.0_real_8) !vw not initialized at all
     REAL(real_8) :: wnosp0 = HUGE(0.0_real_8) !vw not initialized at all
     REAL(real_8) :: epsdav = HUGE(0.0_real_8)
     REAL(real_8) :: amprp = HUGE(0.0_real_8)
     REAL(real_8) :: fdiff = HUGE(0.0_real_8)
     REAL(real_8) :: cmass = HUGE(0.0_real_8)
     REAL(real_8) :: amprc = HUGE(0.0_real_8)
     REAL(real_8) :: tempc = HUGE(0.0_real_8) !vw not initialized at all
     REAL(real_8) :: wnosc0 = HUGE(0.0_real_8) !vw not initialized at all
     REAL(real_8) :: b2limit = HUGE(0.0_real_8)
     REAL(real_8) :: ekinhr = HUGE(0.0_real_8)
     REAL(real_8) :: tolc = HUGE(0.0_real_8)
     REAL(real_8) :: tolkinc = HUGE(0.0_real_8)
     REAL(real_8) :: annerc = HUGE(0.0_real_8)
     REAL(real_8) :: tolcg = HUGE(0.0_real_8)
     REAL(real_8) :: annere = HUGE(0.0_real_8)
     REAL(real_8) :: nedof0 = HUGE(0.0_real_8)
     REAL(real_8) :: dshift = HUGE(0.0_real_8)
     REAL(real_8) :: tolad = HUGE(0.0_real_8)
     REAL(real_8) :: tolene = HUGE(0.0_real_8)
     REAL(real_8) :: tolfor = HUGE(0.0_real_8)
     REAL(real_8) :: tolrhofix = HUGE(0.0_real_8)
     REAL(real_8) :: tolini = HUGE(0.0_real_8)
     REAL(real_8) :: toldetot = HUGE(0.0_real_8)
     REAL(real_8) :: memsize = HUGE(0.0_real_8)
     REAL(real_8) :: asic = HUGE(0.0_real_8)
     REAL(real_8) :: bsic = HUGE(0.0_real_8)
     REAL(real_8) :: nospt0 = HUGE(0.0_real_8) !vw not initialized at all
     REAL(real_8) :: dampgi = HUGE(0.0_real_8)
     REAL(real_8) :: dampge = HUGE(0.0_real_8)
     REAL(real_8) :: dampgc = HUGE(0.0_real_8)
     REAL(real_8) :: gfreq = HUGE(0.0_real_8) !vw not initialized at all
  END TYPE cntr_t
  TYPE(cntr_t), SAVE, PUBLIC :: cntr
  ! strings
  TYPE, PUBLIC :: cnts_t
     CHARACTER(len=255)   :: inputfile = ''
  END TYPE cnts_t
  TYPE(cnts_t), SAVE, PUBLIC :: cnts

  ! ==================================================================
  ! == FILEPATH  OPTION (OTHERWISE="./")                            ==
  ! == DEFINITION AND COMMON BLOCK MOVED TO FILE fileopen.inc       ==
  ! ==================================================================
  ! == DUAL OPTION: SIZE OF THE DENSITY MESH VERSUS                 ==
  ! ==              THE WAVEFUNCTION CUTOFF (DEFAULT=4)             ==
  ! ==      DUAL=.TRUE. IF USE THIS OPTION                          ==
  ! ==--------------------------------------------------------------==
  TYPE, PUBLIC :: dual00_t
     REAL(real_8) :: cdual = HUGE(0.0_real_8)
     LOGICAL :: dual
  END TYPE dual00_t
  TYPE(dual00_t), SAVE, PUBLIC :: dual00
  ! ==================================================================
  INTEGER, PARAMETER, PUBLIC :: maxrf = 10
  TYPE, PUBLIC :: restf_t
     INTEGER :: nstepwr(maxrf) = HUGE(0)
     INTEGER :: nfnow = HUGE(0)
     INTEGER :: nrcopy = HUGE(0)
  END TYPE restf_t
  TYPE(restf_t), SAVE, PUBLIC :: restf
  ! ==================================================================
  INTEGER, SAVE, PUBLIC :: nssel = HUGE(0)
  INTEGER, ALLOCATABLE, SAVE, PUBLIC :: nassel(:)
  ! ==================================================================

  ! ==================================================================
  ! == Control to dump local potential
  ! ==================================================================
  TYPE, PUBLIC :: locpot2_t
     LOGICAL :: tlpot = .FALSE.
  END TYPE locpot2_t
  TYPE(locpot2_t), SAVE, PUBLIC :: locpot2
  ! ==================================================================
  ! ==================================================================
  ! == Control tracing
  ! ==================================================================
  TYPE, PUBLIC :: cp_trace_t
     LOGICAL :: ttrace = .FALSE.
     LOGICAL :: ttrace_master_only = .FALSE.
  END TYPE cp_trace_t
  TYPE(cp_trace_t), SAVE, PUBLIC :: cp_trace
  ! ==================================================================

END MODULE system
