#include "cpmd_global.h"

MODULE control_utils
  USE andr,                            ONLY: andr2,&
                                             andr3,&
                                             maxmix
  USE atwf,                            ONLY: dmovrmix,&
                                             tmovr
  USE benc,                            ONLY: ibench,&
                                             nbentr
  USE broy,                            ONLY: broy1
  USE cdftmod,                         ONLY: &
       cdftci, cdftcom, cdftlog, cdftpred, cdftvgrs, cm_dir, cm_dr, czones, &
       sccomm
  USE comvelmod,                       ONLY: comvl
  USE control_bcast_utils,             ONLY: control_bcast
  USE control_def_utils,               ONLY: control_def
  USE control_pri_utils,               ONLY: control_pri
  USE control_test_utils,              ONLY: control_test
  USE cotr,                            ONLY: sdpl
  USE cp_cuda_types,                   ONLY: cp_cuda_env
  USE envj,                            ONLY: tjlimit
  USE error_handling,                  ONLY: stopgm
  USE fileopenmod,                     ONLY: fo_info
  USE fint,                            ONLY: fint1,&
                                             fint4,&
                                             fint5,&
                                             maxbetap,&
                                             maxtrot
  USE g_loc,                           ONLY: gloc_list,&
                                             glocal,&
                                             gloci,&
                                             glocr
  USE glemod,                          ONLY: &
       gle_cp_ns, gle_cpmd, gle_cust, gle_opt, gle_opt_ns, gle_smart, &
       gle_smart_ns, gle_white, glepar, tglepc
  USE header_utils,                    ONLY: header
  USE hubbardu,                        ONLY: hubbu
  USE inscan_utils,                    ONLY: inscan
  USE ions,                            ONLY: coord_fdiff,&
                                             ions1,&
                                             r_fdiff,&
                                             tref_fdiff
  USE isos,                            ONLY: isos1
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE linres,                          ONLY: tshl,&
                                             xfmqc
  USE lscal,                           ONLY: &
       eps_h, icore, lnomap, lprhes, m_hess, mode, modelk, ncore, nrem, &
       nrestt, nsmaxp, nsvib, ntrstr, ntrust, nvar, omin, tolenv, trustp, &
       trustr
  USE machine,                         ONLY: m_getarg,&
                                             m_iargc,&
                                             m_sleep
  USE mergemod,                        ONLY: merge01,&
                                             merge02
  USE mm_input,                        ONLY: clc,&
                                             lqmmm
  USE mm_parallel,                     ONLY: gparal
  USE nabdy_types,                     ONLY: nabdyfric,&
                                             nabdyvar
  USE nort,                            ONLY: nort_com
  USE nose,                            ONLY: &
       cafesini, cafesinr, lctrng, loct, loctpin, loctt0, ncafesgrp, nosl, &
       tcafes, tnosepc
  USE para_global,                     ONLY: para_buff_size,&
                                             para_stack_buff_size,&
                                             para_use_mpi_in_place
  USE parac,                           ONLY: parai,&
                                             paral
  USE prden,                           ONLY: elfcb,&
                                             mwfn,&
                                             numpr
  USE qspl,                            ONLY: qspl1,&
                                             qsrang
  USE readsr_utils,                    ONLY: index_of_delimiter,&
                                             input_string_len,&
                                             keyword_contains,&
                                             readsi,&
                                             readsr,&
                                             xstring
  USE rlbfgs_utils,                    ONLY: lscal_init
  USE shop,                            ONLY: s0_filn,&
                                             s1_filn,&
                                             sh02,&
                                             tshopold
  USE shop_rest,                       ONLY: sh03
  USE spin,                            ONLY: clsd,&
                                             lspin1,&
                                             lspin2,&
                                             lspin3
  USE store_types,                     ONLY: &
       cprint, iface1, intfn, iprint_coor, iprint_eband, iprint_ebogo, &
       iprint_ecas, iprint_eeig, iprint_egc, iprint_ehee, iprint_ehep, &
       iprint_ehii, iprint_ehsic, iprint_eht, iprint_eigen, iprint_ekin, &
       iprint_elec1, iprint_elec2, iprint_enl, iprint_entropy, iprint_epen, &
       iprint_epseu, iprint_eself, iprint_esr, iprint_etddft, iprint_etot1, &
       iprint_etot2, iprint_exc, iprint_force, iprint_info, iprint_lscal, &
       iprint_vxc, iprint_wann, iprint_ehub, restart1, rout1, store1, &
       trajsmall, trajsmalln
  USE struc,                           ONLY: angle,&
                                             bond,&
                                             dihedral
  USE system,                          ONLY: &
       cnti, cntl, cntr, cnts, cp_trace, group, locpot2, maxrf, maxsys, &
       nassel, nssel, restf
  USE time,                            ONLY: tname
  use bicanonicalInputReader, only: New, Delete, ReadInput,&
     bicanonicalInputReaderType
  use bicanonicalCpmd, only: bicanonicalCpmdInputConfig
  USE vdwcmod,                         ONLY: vdwl
  USE wann,                            ONLY: sw_list,&
                                             wan05,&
                                             wanni,&
                                             wannl,&
                                             wannr
  USE xinr,                            ONLY: gnx_inr,&
                                             inr_integer,&
                                             inr_logical,&
                                             maxreg,&
                                             rmixsd,&
                                             tolx_inr

  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER :: output_unit = 6
  INTEGER, PARAMETER :: maxinfo = 30
  INTEGER, PARAMETER :: max_unknown_lines = 30

  PUBLIC :: control
  PUBLIC :: get_input_name

CONTAINS

  ! ==================================================================
  SUBROUTINE control
    ! ==--------------------------------------------------------------==
    ! ==  THIS ROUTINE READS THE INPUT FILE                           ==
    ! ==--------------------------------------------------------------==
    ! ==  THE INPUT HAS TO BE IN THE FOLLOWING FORMAT                 ==
    ! ==     &CPMD                                                    ==
    ! ==        Options                                               ==
    ! ==     &END                                                     ==
    ! ==--------------------------------------------------------------==
    ! ==  LIST OF OPTIONS                                             ==
    ! ==                                                              ==
    ! ==    MOLECULAR DYNAMICS [CP,BO,PT,                             ==
    ! ==                        FILE XYZ [NSKIP=N] [NSAMPLE=N],       ==
    ! ==                        CLASSICAL]                            ==
    ! ==    PARRINELLO-RAHMAN [NPT,SHOCK,TOLKINC=tolkinc]             ==
    ! ==    OPTIMIZE {GEOMETRY,WAVEFUNCTION,COMBINED} [XYZ] [SAMPLE]  ==
    ! ==      ngxyz                                                   ==
    ! ==    FREE ENERGY FUNCTIONAL                                    ==
    ! ==    INTERFACE {EGO,GMX} [MULLIKEN,LOWDIN,ESP,HIRSHFELD]       ==
    ! ==    PATH SAMPLING                                             ==
    ! ==    KOHN-SHAM ENERGIES {OFF} [NOWAVEFUNCTION]                 ==
    ! ==      nkssta                                                  ==
    ! ==    VIBRATIONAL ANALYSIS [LR,FD,IN] [GAUSS]                   ==
    ! ==      nvib                                                    ==
    ! ==    PROPERTIES                                                ==
    ! ==    FORCEMATCH                                                ==
    ! ==    DEBUG FORCES                                              ==
    ! ==   DEBUG [MEMORY] [RESTART] [FILEOPEN] [NOACC] [IO]           ==
    ! ==    ORBITAL HARDNESS [LR,FD]                                  ==
    ! ==    ELECTRONIC SPECTRA                                        ==
    ! ==    PATH INTEGRALS                                            ==
    ! ==    PATH MINIMISATION                                         ==
    ! ==    SURFACE HOPPING                                           ==
    ! ==    QMMM [QMMMEASY]                                           ==
    ! ==    LINEAR RESPONSE                                           ==
    ! ==    NONORTHOGONAL ORBITALS {OFF}                              ==
    ! ==       slimit                                                 ==
    ! ==    LSD                                                       ==
    ! ==    LOCAL SPIN DENSITY                                        ==
    ! ==    HARMONIC REFERENCE SYSTEM {OFF}                           ==
    ! ==    SCALED MASSES {OFF}                                       ==
    ! ==    STEEPEST DESCENT [ELECTRONS,IONS,CELL]                    ==
    ! ==                     [NOPRECONDITIONING] [LINE]               ==
    ! ==      or {cntl%tsde,cntl%tsdp} [LINE]                                   ==
    ! ==      or {cntl%tsdc}                                               ==
    ! ==    CONJUGATE GRADIENT [ELECTRONS,IONS]                       ==
    ! ==                     [MINIMIZE,NOPRECONDITIONING]             ==
    ! ==      or cntl%pcg  [MINIMIZE,NOPRECONDITIONING]                    ==
    ! ==      or TCGP                                                 ==
    ! ==    ODIIS [NOPRECONDITIONING,NO_RESET=nreset]                 ==
    ! ==      mdiis                                                   ==
    ! ==    HAMILTONIAN CUTOFF                                        ==
    ! ==      hthrs                                                   ==
    ! ==    cntl%gdiis                                                     ==
    ! ==      mgdiis                                                  ==
    ! ==    cntl%bfgs                                                      ==
    ! ==    cntl%lbfgs [NREM,NTRUST,NRESTT,NTRSTR]                         ==
    ! ==      nrem|ntrust|nrestt|ntrstr                               ==
    ! ==    cntl%prfo [NVAR,MODE,MDLOCK,TOLENV,TRUSTP,OMIN,CORE,NSMAXP,    ==
    ! ==          NSVIB,PRJHES,DISPL,HESSTYPE]                        ==
    ! ==      nvar|mode|mdlock|tolenv|trustp|omin|core|nsmaxp|m_hess  ==
    ! ==    cntl%prfo CORE=natcor                                          ==
    ! ==      <atoms_list>                                            ==
    ! ==    HESSCORE                                                  ==
    ! ==    cntl%rfo [ORDER=nsorder]                                       ==
    ! ==    IMPLICIT NEWTON RAPHSON []                                ==
    ! ==        [cntl%prec,CONTINUE,VERBOSE,ALTERNATIVE STEP]              ==
    ! ==      itmax_inr                                               ==
    ! ==    IMPLICIT NEWTON RAPHSON PARAMETERS N=nreg                 ==
    ! ==      gnx_inr(i)   tolx_inr(i)                                ==
    ! ==          ...                                                 ==
    ! ==    MIXDIIS                                                   ==
    ! ==      rmixsd (USED BY INR)                                    ==
    ! ==    MIXSD                                                     ==
    ! ==      rmixsd (USED BY INR)                                    ==
    ! ==    RESTART [WAVEFUNCTION,COORDINATES,DENSITY,ACCUMULATORS,   ==
    ! ==             HESSIAN,VELOCITIES,NOSEE,NOSEP,GEOFILE,KPOINTS,  ==
    ! ==             VIBANALYSIS,CELL,NOSEC,POTENTIAL,OCCUPATION,     ==
    ! ==             ORBHARDNESS,PHESSIAN,LSSTAT,ADPTTL,EXTRAP,GLE,   ==
    ! ==             PRNG,LATEST]                                     ==
    ! ==    INTFILE {READ,WRITE,FILENAME}                             ==
    ! ==      intfn                                                   ==
    ! ==    FILE FUSION                                               ==
    ! ==      s0_filn                                                 ==
    ! ==      s1_filn                                                 ==
    ! ==    FILE MERGE                                                ==
    ! ==      mfiln1                                                  ==
    ! ==      mfiln2                                                  ==
    ! ==      mnst1 mnst2 mortho                                      ==
    ! ==      mshift(1) mshift(2) mshift(3)                           ==
    ! ==    INITIALIZE WAVEFUNCTION {RANDOM,ATOMS}                    ==
    ! ==    MAXSTEP                                                   ==
    ! ==      nomore                                                  ==
    ! ==    MAXITER                                                   ==
    ! ==      nomore_iter                                             ==
    ! ==    MAXRUNTIME                                                ==
    ! ==      tjlimit                                                 ==
    ! ==    ORTHOGONALIZATION {LOWDIN,GRAM-SCHMIDT}                   ==
    ! ==    RATTLE                                                    ==
    ! ==      [ maxit  epsog ]                                        ==
    ! ==    QUENCH [IONS, ELECTRONS, CELL, BO]                        ==
    ! ==    BOGOLIUBOV CORRECTION {OFF}                               ==
    ! ==    TROTTER FACTORISATION OFF                                 ==
    ! ==    TROTTER FACTOR                                            ==
    ! ==      betap                                                   ==
    ! ==    TROTTER FACTOR N=ntabbetap                                ==
    ! ==      densbetap(1)         tabbetap(1)                        ==
    ! ==          ...                                                 ==
    ! ==      densbetap(ntabbetap) tabbetap(ntabbetap)                ==
    ! ==    LANCZOS DIAGONALIZATION [ALL,NEW,OPT,PART,PERI,RESET]     ==
    ! ==    LANCZOS PARAMETER [ALL]                                   ==
    ! ==      n_fries nkry_max nkry_block b2limit                     ==
    ! ==    LANCZOS PARAMETER   N=ntabtrot [ALL]                      ==
    ! ==      n_fries nkry_max  b2limit                               ==
    ! ==     denstrot(2)        b2trot(2)                             ==
    ! ==          ...                                                 ==
    ! ==     denstrot(ntabtrot) b2trot(ntabtrot)                      ==
    ! ==    DAVIDSON DIAGONALIZATION                                  ==
    ! ==    DAVIDSON PARAMETER                                        ==
    ! ==      ndavv epsdav n_fries                                    ==
    ! ==    FIXRHO UPWFN                                              ==
    ! ==    FIXRHO VECT                                               ==
    ! ==       mdiis_fr                                               ==
    ! ==    FIXRHO LOOP                                               ==
    ! ==       minldiis  maxldiis                                     ==
    ! ==    FIXRHO WFTOL                                              ==
    ! ==       tolrhofix                                              ==
    ! ==    ANDERSON MIXING [G-SPACE]                                 ==
    ! ==     andrmix                                                  ==
    ! ==    ANDERSON MIXING N=ntabmix                                 ==
    ! ==     densmix(1)       andmix(1)                               ==
    ! ==          ...                                                 ==
    ! ==     densmix(ntabmix) andmix(ntabmix)                         ==
    ! ==    cntl%diis MIXING                                               ==
    ! ==     nrdiis                                                   ==
    ! ==    cntl%diis MIXING N=ntabnrdiis                                  ==
    ! ==     densnrdiis(1)          nrdiis(1)                         ==
    ! ==          ...                                                 ==
    ! ==     densnrdiis(ntabnrdiis) nrdiis(ntabmix)                   ==
    ! ==    BROYDEN MIXING                                            ==
    ! ==     broymix,ecutbroy,w02broy,nfrbroy,ibreset,kermix          ==
    ! ==    BROYDEN MIXING                                            ==
    ! ==     DEFAULT                                                  ==
    ! ==    BROYDEN MIXING                                            ==
    ! ==     [BROYMIX=broymix] [ECUTBROY=ecutbroy] [W02BROY=w02broy]  ==
    ! ==     [NFRBROY=nfrbroy] [IBRESET=ibreset] [KERMIX=kermix]      ==
    ! ==    ALEXANDER MIXING                                          ==
    ! ==     alxmix                                                   ==
    ! ==    MOVERHO                                                   ==
    ! ==     dmovrmix                                                 ==
    ! ==    EXTRAPOLATE WFN {ASPC,POLY} [STORE] [CSTEPS=naspc]        ==
    ! ==     mextra                                                   ==
    ! ==    RANDOMIZE [WAVEFUNCTION,COORDINATES,CELL,DENSITY]         ==
    ! ==      [ ampre amprp amprc amprd ]                             ==
    ! ==    TIMESTEP                                                  ==
    ! ==      delt                                                    ==
    ! ==    TIMESTEP ELECTRONS                                        ==
    ! ==      delt_elec                                               ==
    ! ==    TIMESTEP IONS                                             ==
    ! ==      delt_ions                                               ==
    ! ==    EMASS                                                     ==
    ! ==      emass                                                   ==
    ! ==    CMASS                                                     ==
    ! ==      cmass                                                   ==
    ! ==    TEMPERATURE  ELECTRON                                     ==
    ! ==      { betael }                                              ==
    ! ==    TEMPERATURE [RAMP]                                        ==
    ! ==      { tempw trampt trampr }                                 ==
    ! ==    RESCALE OLD VELOCITIES                                    ==
    ! ==    REVERSE VELOCITIES                                        ==
    ! ==    SUBTRACT COMVEL                                           ==
    ! ==      ncomv                                                   ==
    ! ==    SUBTRACT ROTVEL                                           ==
    ! ==      nrotv                                                   ==
    ! ==    PRNGSEED                                                  ==
    ! ==      iprng                                                   ==
    ! ==    LANGEVIN [WHITE,CPMD,OPTIMAL,SMART|CUSTOM] [MOVECM]       ==
    ! ==      [ omega0, custom_ns ]                                   ==
    ! ==    TEMPCONTROL {IONS,ELECTRONS,CELL}                         ==
    ! ==      { tempw tolp , ekinw  toll , ekinhr tolc}               ==
    ! ==    BERENDSEN {IONS,ELECTRONS,CELL}                           ==
    ! ==      { tempw taubp, ekinw taube,  ekinhr taubc}              ==
    ! ==    NOSE [IONS,ELECTRONS,CELL] [ULTRA,MASSIVE,CAFES]          ==
    ! ==      { tempw wnosep , ekinw wnosee , tempc wnosec}           ==
    ! ==    NOSE PARAMETERS                                           ==
    ! ==      nchp nche nchb nedof ncalls nit                         ==
    ! ==    PRINT [INFO,EIGENVALUES,COORDINATES,FORCES,WANNIER]       ==
    ! ==      iprint_step                                             ==
    ! ==    PRINT {ON,OFF} [INFO,EIGENVALUES,COORDINATES,FORCES,      ==
    ! ==                    WANNIER]                                  ==
    ! ==    PRINT ENERGY {ON,OFF} [EKIN,ELECTROSTATIC,ELEC1,ELEC2,    ==
    ! ==                    ESR,ESELF,                                ==
    ! ==                    EFREE,EBAND,ENTROPY,EPSEU,EHEP,EHEE,EHII, ==
    ! ==                    ENL,EXC,VXC,EGC,EBOGO,ETOT1,ETOT2,ECAS,   ==
    ! ==                    ETDDFT,EPEN]                              ==
    ! ==    STRUCTURE [BONDS,ANGLES,DIHEDRALS,SELECT]                 ==
    ! ==     natom                                                    ==
    ! ==     n1, n2, n3, ... nn                                       ==
    ! ==    ISOLATED MOLECULE                                         ==
    ! ==    CENTER MOLECULE [OFF]                                     ==
    ! ==    STORE [WAVEFUNCTION,DENSITY,POTENTIAL]                    ==
    ! ==      istore [SC=isctore]                                     ==
    ! ==    STORE OFF [WAVEFUNCTION,DENSITY,POTENTIAL]                ==
    ! ==    RESTFILE [SAMPLE]                                         ==
    ! ==      nrestf                                                  ==
    ! ==      [s1 s2 ...]                                             ==
    ! ==    CONVERGENCE {ORBITALS,GEOMETRY,CELL}                      ==
    ! ==      tolog   tolng   tolcg                                   ==
    ! ==    CONVERGENCE {ADAPT,ENERGY,CALFOR,RELAX,INITIAL}           ==
    ! ==      tolad,tolene,tolfor,nstcnv,initial                      ==
    ! ==    ANNEALING [IONS,ELECTRONS,CELL]                           ==
    ! ==      {anneri,annere,annerc}                                  ==
    ! ==    DAMPING [IONS,ELECTRONS,CELL]                             ==
    ! ==      {dampgi,dampge,dampgc}                                  ==
    ! ==    HESSIAN {DISCO,SCHLEGEL,UNIT} [PARTIAL]                   ==
    ! ==    PROJECT {NONE,DIAGONAL,FULL}                              ==
    ! ==    STRESS TENSOR [VIRIAL]                                    ==
    ! ==      npres                                                   ==
    ! ==    CLASSTRESS                                                ==
    ! ==      nprec                                                   ==
    ! ==    RHOOUT  [BANDS,SAMPLE=nrhooout]                           ==
    ! ==      numbands                                                ==
    ! ==      n1 n2 ... nx                                            ==
    ! ==    ELF [PARAMETER]                                           ==
    ! ==      elfcut elfeps                                           ==
    ! ==    ENERGYBANDS                                               ==
    ! ==    EXTERNAL POTENTIAL [ADD]                                  ==
    ! ==    ELECTROSTATIC POTENTIAL  [SAMPLE=nrhooout]                ==
    ! ==    DIPOLE DYNAMICS [SAMPLE,WANNIER]                          ==
    ! ==      npdip                                                   ==
    ! ==    WANNIER PARAMETER                                         ==
    ! ==      w_step w_eps w_ran w_maxs                               ==
    ! ==    WANNIER OPTIMIZATION {SD,JACOBI}                          ==
    ! ==    WANNIER SERIAL                                            ==
    ! ==    WANNIER TYPE {VANDERBILT,RESTA}                           ==
    ! ==    WANNIER REFERENCE                                         ==
    ! ==      w_ref(1..3)                                             ==
    ! ==    WANNIER DOS                                               ==
    ! ==    WANNIER MOLECULAR                                         ==
    ! ==    WANNIER WFNOUT [ALL,PARTIAL,LIST,DENSITY,SPREAD]          ==
    ! ==      sw_orb                                                  ==
    ! ==      sw_list(1) ...                                          ==
    ! ==    TRAJECTORY [XYZ] [SAMPLE] {OFF}                           ==
    ! ==      ntraj                                                   ==
    ! ==    FINITE DIFFERENCES                                        ==
    ! ==     fdiff COORD=coord_fdiff(1..3) RADIUS=r_fdiff             ==
    ! ==    MOVIE [SAMPLE] {OFF}                                      ==
    ! ==      imovie                                                  ==
    ! ==    MEMORY {SMALL,BIG}                                        ==
    ! ==    COMPRESS WRITEnn                                          ==
    ! ==    FILEPATH                                                  ==
    ! ==      fpath                                                   ==
    ! ==    TASKGROUPS {MAXIMUM,MINIMUM,CARTESIAN}                    ==
    ! ==      nogrp                                                   ==
    ! ==    DISTRIBUTE FNL [ON,OFF]                                   ==
    ! ==    SPLINE [POINTS QFUNCTION INIT RANGE]                      ==
    ! ==      nsplp qsrang                                            ==
    ! ==    REAL SPACE WFN [KEEP, SIZE]                               ==
    ! ==      memsize                                                 ==
    ! ==    BENCHMARK                                                 ==
    ! ==     b1 b2 b3 b4 b5 b6 b7 b8 b9 b0                            ==
    ! ==    CHECK MEMORY                                              ==
    ! ==    MIRROR                                                    ==
    ! ==    SHIFT POTENTIAL                                           ==
    ! ==     dshift                                                   ==
    ! ==    GLOCALIZATION PARAMETER                                   ==
    ! ==      gloc_step gloc_eps gloc_ran gloc_maxs                   ==
    ! ==    GLOCALIZATION OPTIMIZATION                                ==
    ! ==        {SD INTEGRATION,ALMOST IDENTITY}                      ==
    ! ==    GFUNCTIONAL TYPE {GSPREAD,ZETANUMBER}                     ==
    ! ==    UNITARITY CONSTRAINT                                      ==
    ! ==       {TGLAGRANGE, TGORTHO, TGAPPROX}                        ==
    ! ==        if tglagrange: G_LAGRANGE PARAMETERS                  ==
    ! ==         epslag , gloc_maxit                                  ==
    ! ==    RANDOM KICK                                               == 
    ! ==       gloc_kick                                              ==
    ! ==    G_CONJUGATE GRADIENT (only if SD INTEGRATION )            ==
    ! ==    GLOCALIZE STARTING WAFEFUNCTION                           ==
    ! ==      {complex, real}                                         ==
    ! ==    G_STEP MODULATION                                         ==
    ! ==       parameters                                             ==
    ! ==    GLOCALAZED WFNOUT [ALL,PARTIAL,LIST,DENSITY]              ==
    ! ==      tglocprint                                              ==
    ! ==      gloc_orb                                                ==
    ! ==      gloc_list(1) ...                                        ==
    ! ==    NOGEOCHECK                                                ==
    ! ==    ALLTOALL [SINGLE,DOUBLE]                                  ==
    ! ==    DISTRIBUTED LINALG [ON,OFF]                               ==
    ! ==    BLOCKSIZE STATES                                          ==
    ! ==     nstblk                                                   ==
    ! ==    GSHELL                                                    ==
    ! ==    VDW CORRECTION                                            ==
    ! ==    VDW WANNIER                                               ==
    ! ==    BICANONICAL ENSEMBLE                                      ==
    ! ==--------------------------------------------------------------==
    CHARACTER(len=*), PARAMETER              :: procedureN = 'control'

    CHARACTER(len=30)                        :: fformat
    CHARACTER(len=input_string_len)          :: line, previous_line, next_words, &
                                                error_message, infomsg(maxinfo), &
                                                unknown(max_unknown_lines)
    INTEGER                                  :: first, last, keep_first, ierr, &
                                                iunit, nbr_unknown_lines, &
                                                nbr_info_lines, counter
    INTEGER                                  :: i, iflag
    LOGICAL                                  :: something_went_wrong, go_on_reading, &
                                                wait_for_qmstart, is_there
    LOGICAL                                  :: erread, mirror, test, store_or_not, &
                                                tsrho, tsyssab
    REAL(real_8)                             :: var
    type (BicanonicalInputReaderType) :: bicanonicalInput

    !
    ! Block data with variable size arrays
    !
    CALL lscal_init()
    !
    ! Defaults that are not in control_def, part I
    ! QM/MM code init
    !
    lqmmm%qmmm    = .false.
    paral%qmnode  = .true.
    gparal%mmnode = .false.
    !
    IF (paral%io_parent) THEN
       iunit = 5
       ! 
       ! Test whether the extension of the input file is .run
       ! if this is the case wait for the file QMSTART before
       ! reading in the input file (needed for the interface mode).
       !
       IF (INDEX(cnts%inputfile,'.run') /= 0) THEN
          wait_for_qmstart = .true.
          DO WHILE(wait_for_qmstart)
             WRITE(output_unit,*) ' CPMD is waiting for QMSTART.'
             INQUIRE(file='QMSTART',exist=is_there)
             IF (is_there) THEN
                WRITE(output_unit,*) ' CPMD detected QMSTART.'
                OPEN(iunit,file='QMSTART',status='UNKNOWN')
                CLOSE(iunit,status='DELETE')
                wait_for_qmstart = .false.
             ENDIF
             CALL m_sleep(1)
          END DO
       ENDIF
       !
       ! Input is opened here and remains opened for subsequent routines
       !
       OPEN(unit=iunit,file=trim(adjustl(cnts%inputfile)),iostat=ierr,status='OLD')
       !
       IF (ierr /= 0) CALL stopgm(procedureN,"Input file not found / unreadable: "//trim(adjustl(cnts%inputfile)),&
                                  __LINE__, __FILE__)
       !
       ! Copy the info section
       !
       nbr_info_lines = 0
       ierr = inscan(iunit,'&INFO')
       IF (ierr == 0) THEN
          go_on_reading        = .true.
          DO WHILE (go_on_reading)
             READ(iunit,'(A)',iostat=ierr) line
             IF ( keyword_contains(line,'&END') .OR. &
                 nbr_info_lines >= maxinfo .OR. &
                 ierr /= 0) THEN
                go_on_reading        = .false.
             ELSE
                nbr_info_lines = nbr_info_lines + 1
                infomsg(nbr_info_lines) = line
             ENDIF
          ENDDO
       ENDIF
       !
       ! Variables for reading
       !
       nbr_unknown_lines = 0
       line          = ' '
       previous_line = ' '
       error_message        = ' '
       !
       ! Defaults that are not in control_def, part II
       !
       mirror        = .FALSE.
       tsrho         = .FALSE.
       clc%classical = .FALSE.
       trajsmall     = .FALSE.
       trajsmalln    = 1
       CALL control_def()
       !
       ! Search for Control section
       !     
       ierr = inscan(iunit,'&CPMD')
       !
       IF (ierr == 0) THEN
          !
          ! Main loop
          !
          go_on_reading        = .true.
          something_went_wrong = .false.
          !
          DO WHILE(go_on_reading)
             !
             ! Read a line and store the old one
             !
             previous_line = line
             READ(iunit,'(A80)',iostat=ierr) line
             IF (ierr /= 0) THEN
                something_went_wrong = .TRUE.
                go_on_reading        = .FALSE.
             ELSEIF ( keyword_contains(line,'&END') ) THEN
                go_on_reading        = .FALSE.
             ELSE IF ( keyword_contains(line,'BICANONICAL') .AND.  keyword_contains(line,'ENSEMBLE')) THEN
                call New(bicanonicalInput, iunit)
                call ReadInput(bicanonicalInput, bicanonicalCpmdInputConfig)
                call Delete(bicanonicalInput)
             ELSEIF ( keyword_contains(line,'CDFT') ) THEN
                cntl%cdft=.TRUE.
                IF ( keyword_contains(line,'SPIN'))  cdftlog%tspinc = .TRUE.
                IF ( keyword_contains(line,'PCGFI')) cntl%tpcgfic   = .TRUE.
                IF ( keyword_contains(line,'ALL') ) THEN
                   cdftlog%tspinc=.FALSE.
                   cdftlog%tcall=.TRUE.
                ENDIF
                IF ( keyword_contains(line,'NEWTON') ) THEN
                   cntl%cdft_newt=.TRUE.
                ELSE IF ( keyword_contains(line,'DEKKER')) THEN
                   cntl%cdft_dekk=.TRUE.
                ENDIF
                IF ( keyword_contains(line,'HDA') ) THEN
                   cdftlog%thda=.TRUE.
                   cntl%krwfn=.TRUE.
                   IF ( keyword_contains(line,'PROJECT'))cdftlog%thdaproj=.TRUE.
                ENDIF
                IF ( keyword_contains(line,'AUTO'))cdftlog%tauto=.TRUE.
                IF ( keyword_contains(line,'NOCCOR'))cdftlog%ccor=.FALSE.
                IF ( keyword_contains(line,'RESWF'))cdftlog%reswf=.TRUE.
                IF ( keyword_contains(line,'PHIOUT'))cdftlog%tphio=.TRUE.
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                first=1
                CALL readsr(line,first,last,cdftcom%cdft_nc,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                first=last
                IF (cdftlog%tcall) THEN
                   CALL readsr(line,first,last,cdftcom%cdft_ns,erread)
                   IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                   first=last
                ENDIF
                CALL readsr(line,first,last,cdftcom%cdft_v(1),erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                first=last
                IF (cdftlog%tcall) THEN
                   CALL readsr(line,first,last,cdftcom%cdft_v(2),erread)
                   IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                   first=last
                ENDIF
                CALL readsi(line,first,last,cdftci%cdft_end,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                IF (cdftlog%thda) THEN
                   previous_line = line
                   READ(iunit,'(A)',iostat=ierr) line
                   first=1
                   CALL readsr(line,first,last,cdftcom%nother,erread)
                   IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                   first=last
                   IF (cdftlog%tcall) THEN
                      CALL readsr(line,first,last,cdftcom%nsother,erread)
                      IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                   ENDIF
                ENDIF
             ELSEIF ( keyword_contains(line,'VGFACTOR') ) THEN
                READ(iunit,'(A)',iostat=ierr) line
                first=1
                CALL readsr(line,first,last,cdftvgrs%vgstep,erread)
                IF (erread.OR.cdftvgrs%vgstep<=0.0_real_8) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                first=last
                CALL readsr(line,first,last,cdftvgrs%maxvmov,erread)
                IF (erread.OR.cdftvgrs%maxvmov<=0.0_real_8) cdftvgrs%maxvmov=0.1_real_8
             ELSEIF ( keyword_contains(line,'VMIRROR') ) THEN
                cdftlog%tvmirr=.TRUE.
             ELSEIF ( keyword_contains(line,'COMBINE', and='SYSTEMS') ) THEN
                cntl%tsyscomb=.TRUE.
                IF ( keyword_contains(line,'REFLECT'))cntl%tscombm=.TRUE.
                IF ( keyword_contains(line,'NONORTH'))cntl%tscortho=.FALSE.
                IF ( keyword_contains(line,'SAB'))tsyssab=.TRUE.
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                first=1
                CALL readsi(line,first,last,sccomm%n_s0,erread)
                first=last
                CALL readsi(line,first,last,sccomm%n_s1,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                first=last
                CALL readsi(line,first,last,sccomm%n_s0up,erread)
                first=last
                CALL readsi(line,first,last,sccomm%n_s1up,erread)
                IF (cntl%tscombm) THEN
                   first=last
                   CALL readsi(line,first,last,cm_dir,erread)
                   IF (cm_dir==4) THEN
                      first=last
                      CALL readsr(line,first,last,cm_dr,erread)
                   ENDIF
                ENDIF
                IF (tsyssab) THEN
                   previous_line = line
                   READ(iunit,'(A)',iostat=ierr) line
                   first=1
                   CALL readsi(line,first,last,sccomm%tsysk,erread)
                   first=last
                   CALL readsi(line,first,last,sccomm%tsysl,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ENDIF
                erread=.FALSE.
             ELSEIF ( keyword_contains(line,'KSHAM') ) THEN
                cntl%tksham=.TRUE.
                IF ( keyword_contains(line,'MATRIX')) cntl%tlowdmat=.TRUE.
                IF ( keyword_contains(line,'ROUT') ) THEN
                   rout1%rhoout=.TRUE.
                   numpr=2
                   ALLOCATE(mwfn(numpr),STAT=ierr)
                   IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                        __LINE__,__FILE__)
                ENDIF
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                first=1
                CALL readsi(line,first,last,sccomm%n_s0up,erread)
                first=last
                CALL readsi(line,first,last,sccomm%n_s1up,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ! Convergence Zones for cntl%cdft
             ELSEIF ( keyword_contains(line,'CZONES') ) THEN
                cdftlog%tczones=.TRUE.
                IF ( keyword_contains(line,'SET') ) THEN
                   previous_line = line
                   READ(iunit,'(A)',iostat=ierr) line
                   first=1
                   CALL readsr(line,first,last,czones(1,2),erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   first=last
                   CALL readsr(line,first,last,czones(2,2),erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   first=last
                   CALL readsr(line,first,last,czones(3,2),erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   previous_line = line
                   READ(iunit,'(A)',iostat=ierr) line
                   first=1
                   CALL readsr(line,first,last,czones(2,1),erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   first=last
                   CALL readsr(line,first,last,czones(3,1),erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ENDIF
             ! cntl%cdft write out weight slice
             ELSEIF ( keyword_contains(line,'WOUT') ) THEN
                cntl%cdft_weight=.TRUE.
                IF ( keyword_contains(line,'FULL') ) THEN
                   cntl%cdft_wf=.TRUE.
                ELSE
                   cntl%cdft_wf=.FALSE.
                   previous_line = line
                   READ(iunit,'(A)',iostat=ierr) line
                   first=1
                   CALL readsr(line,first,last,cdftcom%wslice,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   first=last
                   CALL readsi(line,first,last,cdftci%wstep,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ENDIF
             ! EXACT FACT
             ELSEIF ( keyword_contains(line,'XFMQC') ) THEN
                tshl%txfmqc=.TRUE.
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                CALL readsi(line ,1,last,xfmqc%n_xftraj,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF ( keyword_contains(line,'MOLECULAR',and='DYNAMICS')) THEN
                ! Molecular Dynamics
                cntl%md=.TRUE.
                IF ( keyword_contains(line,'BO')) cntl%tmdbo=.TRUE.
                IF ( keyword_contains(line,'FILE') ) THEN
                   IF ( keyword_contains(line,'XYZ')) rout1%xtin=.TRUE.
                   cntl%tmdfile=.TRUE.
                   IF ( keyword_contains(line,'NSKIP',cut_at='=') ) THEN
                      first = index_of_delimiter(line,'NSKIP','=') 
                      CALL readsi(line,first,last,cnti%nskip,erread)
                      write(*,*) 'SKIP', first, line(first:), cnti%nskip
                      IF (erread) THEN
                         error_message        = "ERROR WHILE READING VALUE (NSKIP)"
                         something_went_wrong = .true.
                         go_on_reading        = .false.
                      ENDIF
                   ENDIF
                   IF ( keyword_contains(line,'NSAMPLE',cut_at='=') ) THEN
                      first = index_of_delimiter(line,'NSAMPLE','=')
                      CALL readsi(line,first,last,cnti%nsample,erread)
                      write(*,*) 'SAMPLE', first, line(first:), cnti%nsample
                      IF (erread) THEN
                         error_message        = "ERROR WHILE READING VALUE (NSAMPLE)"
                         something_went_wrong = .true.
                         go_on_reading        = .false.
                      ENDIF
                   ENDIF
                ENDIF
                IF (cntl%tmdbo.AND.keyword_contains(line,'PT')) cntl%tresponse=.TRUE.
                IF ( keyword_contains(line,'CP')) cntl%tmdbo=.FALSE.
                IF ( keyword_contains(line,'CLASSICAL')) clc%classical=.TRUE.
                ! EHR[
                IF ( keyword_contains(line,'EH',alias='EHRENFEST') ) THEN
                   cntl%tmdbo=.FALSE.
                   cntl%tmdeh=.TRUE.
                   cntl%cmplx_wf=.TRUE.
                ENDIF
                ! EHR]
                ! NABDY[
                IF ( keyword_contains(line,'BD') ) THEN
                   cntl%tmdbo=.FALSE.
                   cntl%tnabdy=.TRUE.
                   !Number of trajectories
                   previous_line = line
                   READ(iunit,'(A)',iostat=ierr) line
                   CALL readsi(line ,1,last,nabdyvar%ntrajbd,erread)
                   IF(erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ENDIF
                ! NABDY]
             ELSEIF ( keyword_contains(line,'PARRINELLO',and='RAHMAN',alias='PARRINELLO-RAHMAN')) THEN
                cntl%tprcp=.TRUE.
                IF ( keyword_contains(line,'NPT')) cntl%tpenn=.TRUE.
                IF ( keyword_contains(line,'SHOCK') ) THEN
                   cntl%tpenn=.TRUE.
                   cntl%tshock=.TRUE.
                ENDIF
                IF ( keyword_contains(line,'TOLKINC',cut_at='=') ) THEN
                   first = index_of_delimiter(line,'TOLKINC','=')
                   CALL readsr(line,first,last,cntr%tolkinc,erread)
                   IF(erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ENDIF
             ELSEIF ( keyword_contains(line,'OPTIMIZE',and='GEOMETRY') .OR. &
                      keyword_contains(line,'GEOMETRY',and='OPTIMIZATION') ) THEN
                ! Optimisation of Ionic Coordinates
                cntl%geopt=.TRUE.
                cprint%iprint(iprint_force) = 1
                IF ( keyword_contains(line,'CLASSICAL')) clc%classical=.TRUE.
                IF ( keyword_contains(line,'XYZ') ) THEN
                   rout1%xgout=.TRUE.
                   IF ( keyword_contains(line,'SAMPLE') ) THEN
                      previous_line = line
                      READ(iunit,'(A)',iostat=ierr) line
                      CALL readsi(line,1,last,cnti%ngxyz,erread)
                      IF (erread) THEN
                         error_message        = "ERROR WHILE READING VALUE"
                         something_went_wrong = .true.
                         go_on_reading        = .false.
                      ENDIF
                   ENDIF
                ENDIF
             ELSEIF ( keyword_contains(line,'OPTIMIZE',and='COMBINED') .OR. &
                      keyword_contains(line,'COMBINED',and='OPTIMIZATION') ) THEN
                ! Optimisation of Ionic Coordinates (combined scheme)
                cntl%geopt=.TRUE.
                cntl%simul=.TRUE.
                cntl%tharm=.TRUE.
                cntl%tmass=.TRUE.
                IF ( keyword_contains(line,'XYZ') ) THEN
                   rout1%xgout=.TRUE.
                   IF ( keyword_contains(line,'SAMPLE') ) THEN
                      previous_line = line
                      READ(iunit,'(A)',iostat=ierr) line
                      CALL readsi(line,1,last,cnti%ngxyz,erread)
                      IF (erread) THEN
                         error_message        = "ERROR WHILE READING VALUE"
                         something_went_wrong = .true.
                         go_on_reading        = .false.
                      ENDIF
                   ENDIF
                ENDIF
             ELSEIF ( keyword_contains(line,'OPTIMIZE',and='WAVEFUNCTION') .OR. &
                      keyword_contains(line,'WAVEFUNCTION',and='OPTIMIZATION') ) THEN
                ! Optimisation of Wavefunction
                cntl%wfopt=.TRUE.
             ELSEIF ( keyword_contains(line,'CHEBY') ) THEN
                cntl%cheby=.TRUE.
                cntl%cayley=.FALSE.
             ELSEIF ( keyword_contains(line,'CAYLEY') ) THEN
                cntl%cayley=.TRUE.
                cntl%cheby=.FALSE.
             ELSEIF ( keyword_contains(line,'RUNGE',and='KUTTA',alias='RUNGE-KUTTA') ) THEN
                cntl%ruku=.TRUE.
             ELSEIF ( keyword_contains(line,'FORCEMATCH') ) THEN
                ! Forcematching
                cntl%fmatch=.TRUE.
             ELSEIF ( keyword_contains(line,'DEBUG') ) THEN
                ! Debug force calculation
                IF ( keyword_contains(line,'FORCES') ) THEN
                   cntl%tdebfor=.TRUE.
                ENDIF
                ! Debug FILEOPEN processing
                IF ( keyword_contains(line,'FILE') ) THEN
                   fo_info%fo_tdebug=.TRUE.
                ENDIF
                ! Debug RESTART file read or write
                IF ( keyword_contains(line,'IO') ) THEN
                   store1%tdebio=.TRUE.
                ENDIF
                ! Do not read/write accumulator information from/to the RESTART file
                ! This avoids putting cntl%timing information to the RESTART and makes 
                ! RESTART files identical for identical runs.
                IF ( keyword_contains(line,'NOACC') ) THEN
                   store1%tdebacc=.TRUE.
                ENDIF
                IF ( keyword_contains(line,'HUBBARD') ) THEN
                   hubbu%debug=.TRUE.
                ENDIF
             ELSEIF ( keyword_contains(line,'KOHN-SHAM',and='ENERGIES',alias='KS-EIGENVALUES') ) THEN
                ! Calculate Kohn-Sham energies for the final potential
                IF ( keyword_contains(line,'OFF') ) THEN
                   cntl%ksener=.FALSE.
                ELSE
                   cntl%ksener=.TRUE.
                   vdwl%vdwc=.FALSE.  ! Empirical vdW affect the number of electrons
                ENDIF
                IF (cntl%ksener) THEN
                   IF ( keyword_contains(line,'NOWAVEFUNCTION')) tkpts%tknoswap=.TRUE.
                   previous_line = line
                   READ(iunit,'(A)',iostat=ierr) line
                   CALL readsi(line,1,last,cnti%nkssta,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   IF (cnti%nkssta<0) THEN
                      error_message        = 'NKSSTA HAS TO BE GREATER OR EQUAL TO ZERO'
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ENDIF
             ELSEIF ( keyword_contains(line,'SURFACE',and='HOPPING') ) THEN
                ! Surface Hopping molecular dynamics                 
                cntl%tshop=.TRUE.
                tshopold=.TRUE.
                sh03%tshopres=.FALSE.
                sh02%eaddsh=0.0_real_8
                IF ( keyword_contains(line,'NEW') )     tshopold=.FALSE.
                IF ( keyword_contains(line,'RESTART') ) sh03%tshopres=.TRUE.
                IF ( keyword_contains(line,'SHIFT',cut_at='=') ) THEN
                   first = index_of_delimiter(line,'SHIFT','=')
                   CALL readsr(line,first,last,sh02%eaddsh,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ENDIF
                IF ( keyword_contains(line,'NEW') ) tshopold=.FALSE.
                sh02%nsurf=2
             ELSEIF ( keyword_contains(line,'ROKS') ) THEN
                ! RESTRICTED OPEN SHELL KOHN SHAM
                lspin2%tros=.TRUE.
                clsd%nlsd=4
                clsd%nlsx=3
                ! SINGLET/TRIPLET: PARAMETERS FOR LOW SPIN EXCITATION
                IF ( keyword_contains(line,'TRIPLET') ) THEN
                   lspin1%lsea=0._real_8
                   lspin1%lseb=1._real_8
                ENDIF
                ! SINGLET IS DEFAULT
                IF ( keyword_contains(line,'GOEDECKER',but_not='LOCALIZED') ) THEN
                   ! (UNMODIFIED) GOEDECKER
                   lspin3%mgab(9) = 0.5_real_8
                   lspin3%mgab(10) = 0.5_real_8
                ELSEIF ( keyword_contains(line,'EXPERT') ) THEN
                   READ(iunit, fmt=*, iostat=ierr) lspin3%mgab(1), lspin3%mgab(2),&
                        lspin3%mgab(3), lspin3%mgab(4), lspin3%mgab(5), lspin3%mgab(6), lspin3%mgab(7),&
                        lspin3%mgab(8), lspin3%mgab(9), lspin3%mgab(10), lspin3%mgab(11), lspin3%mgab(12)
                   ! MGAB(1) = A_AC, MGAB(2) = B_AC, MGAB(3) = A_BC,
                   ! MGAB(4) = B_BC, MGAB(5) = A_CA, MGAB(6) = B_CA,
                   ! MGAB(7) = A_CB, MGAB(8) = B_CB, MGAB(9) = A_BA,
                   ! MGAB(10)= B_BA, MGAB(11)= A_AB, MGAB(12)= B_AB
                   IF (ierr /= 0) THEN
                       error_message        = "ERROR WHILE READING EXPERT VALUES"
                       something_went_wrong = .true.
                       go_on_reading        = .false.
                   ENDIF
                ELSEIF ( keyword_contains(line,'LOCALIZED',but_not='DELOCALIZED') ) THEN
                   ! ANTISYMMETRIZED GOEDECKER (LOCALIZED VARIANT)
                   lspin3%mgab(9) = 0.5_real_8
                   lspin3%mgab(10) = -0.5_real_8
                   lspin3%mgab(12) = -0.5_real_8
                   ! MODIFIED GOEDECKER (DELOCALIZED VARIANT) IS DEFAULT
                ENDIF
             ELSEIF ( keyword_contains(line,'PATH',and='SAMPLING') ) THEN
                ! MC sampling of reaction path                       
                cntl%tsampl=.TRUE.
             ELSEIF ( keyword_contains(line,'FREE',and='ENERGY',alias='FREE-ENERGY') ) THEN
                ! Free Energy functional (by default trotter and bogoliubov).
                cntl%tfint=.TRUE.
                cntl%tdiag=.TRUE.
                cntl%tlanc=.TRUE.
                fint1%ttrot=.TRUE.
                fint1%tbogo=.TRUE.
                fint1%betap=0.001_real_8
                ! Wfn Opt. via cntl%diis at fixed RHO 
                ! (recommended for metals with virtual states)
             ELSEIF ( keyword_contains(line,'FIXRHO',and='UPWFN') ) THEN
                IF ( keyword_contains(line,'PCG') ) THEN
                   cntl%pcg=.TRUE.
                   cntl%diis=.FALSE.
                ELSE
                   cntl%pcg=.FALSE.
                   cntl%diis =.TRUE.
                ENDIF
                cntl%tdiag=.TRUE.
                cntl%tlanc=.FALSE.
                cntl%tdavi=.FALSE.
                cntl%diis =.TRUE.
                cntl%prec=.TRUE.
                cntl%tfrho_upw=.TRUE.
                cntl%tfint=.TRUE.
                fint1%ttrot=.FALSE.
                fint1%tbogo=.FALSE.
                fint1%betap=0.001_real_8
                ! Number of vectors for cntl%diis storage
             ELSEIF ( keyword_contains(line,'FIXRHO',and='VECT')) THEN
                READ(iunit,'(A)',iostat=ierr) line
                CALL readsi(line,1,last,cnti%mdiis_fr,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                ! Min and Max Number of cntl%diis Cycles at FIXRHO
             ELSEIF ( keyword_contains(line,'FIXRHO',and='LOOP')) THEN
                READ(iunit,'(A)',iostat=ierr) line
                CALL readsi(line,1,last,cnti%minldiis,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                first = last+1
                CALL readsi(line,first,last,cnti%maxldiis,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                ! Tolerance parameter for cntl%diis convergence of wfn at FIXRHO
             ELSEIF ( keyword_contains(line,'FIXRHO',and='WFTOL')) THEN
                READ(iunit,'(A)',iostat=ierr) line
                CALL readsr(line,1,last,cntr%tolrhofix,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                ! Bogoliubov correction
             ELSEIF ( keyword_contains(line,'BOGOLIUBOV',and='CORRECTION') ) THEN
                IF ( keyword_contains(line,'OFF') ) THEN
                   fint1%tbogo=.FALSE.
                ELSE
                   fint1%tbogo=.TRUE.
                ENDIF
             ELSEIF ( keyword_contains(line,'VIBRATIONAL',and='ANALYSIS') ) THEN
                ! Vibrational Analysis          
                cntl%vibrat=.TRUE.
                cntl%tsdin=.FALSE.
                mode=0
     
                IF ( keyword_contains(line,'FD')) cntl%tsdan=.FALSE.
                IF ( keyword_contains(line,'LR')) cntl%tsdan=.TRUE.
                IF ( keyword_contains(line,'IN')) cntl%tsdin=.TRUE.
                IF ( keyword_contains(line,'ACLIMAX')) rout1%acout=.TRUE.
                IF ( keyword_contains(line,'GAUSS')) rout1%vout=.TRUE.
                IF ( keyword_contains(line,'SAMPLE') ) THEN
                   previous_line = line
                   READ(iunit,'(A)',iostat=ierr) line
                   CALL readsi(line,1,last,cnti%nvib,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ENDIF
                IF ( keyword_contains(line,'MODE',cut_at='=') ) THEN
                   first = index_of_delimiter(line,'MODE','=')
                   CALL readsi(line,first,last,mode,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ENDIF
             ELSEIF ( keyword_contains(line,'ELECTRONIC',and='SPECTRA') ) THEN
                ! Electronic spectra using cntl%tddft
                cntl%tspec=.TRUE.
                ! EHR[
                !CSOC[
             ELSEIF ( keyword_contains(line,'SPIN-ORBIT',and='COUPLING') ) THEN
                cntl%tsoc=.TRUE.
                READ(iunit,'(A)',iostat=ierr) line
                first=1
                CALL readsi(line,first,last,cnti%isocs,erread)
                IF(erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                first=last
                CALL readsi(line,first,last,cnti%jsoct,erread)
                IF(erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                !CSOC]
             ELSEIF ( keyword_contains(line,'PROPAGATION',and='SPECTRA') ) THEN
                ! Electronic spectra using real time propagation DFT
                cntl%cmplx_wf=.TRUE.
                cntl%tpspec=.TRUE.
                cntl%tsde=.TRUE.
                isos1%tcent=.FALSE.
             ELSEIF&
                  (keyword_contains(line,'PROPAGATION',and='DISTRUB') ) THEN
                ! Electronic spectra using real time propagation DFT
                cntl%cmplx_wf=.TRUE.
                cntl%tpdist=.TRUE.
                cntl%tsde=.TRUE.
             ELSEIF ( keyword_contains(line,'GAUGEPULSE') ) THEN
                ! Constant external electronic potential added as gauge field
                cntl%tgaugep=.TRUE.
             ELSEIF ( keyword_contains(line,'GAUGEFIELD') ) THEN
                READ(iunit,'(A)',iostat=ierr) line
                CALL readsr(line ,1,last,cntr%gfreq,erread)
                ! Constant external electronic potential added as gauge field
                cntl%tgaugef=.TRUE.
                cntl%tgaugep=.TRUE.
                ! EHR]
             ELSEIF ( keyword_contains(line,'NACV') ) THEN
                ! Non adiabatic coupling vectors using cntl%tddft
                cntl%tnacvs=.TRUE.
             ELSEIF ( keyword_contains(line,'ORBITAL',and='HARDNESS')) THEN
                ! Orbital Hardness Matrix
                cntl%thard=.TRUE.
                IF ( keyword_contains(line,'FD')) cntl%tsdan=.FALSE.
                IF ( keyword_contains(line,'LR')) cntl%tsdan=.TRUE.
             ELSEIF ( keyword_contains(line,'PROPERTIES',alias='PROPERTY') ) THEN
                ! Properties                    
                cntl%proper=.TRUE.
             ELSEIF ( keyword_contains(line,'PATH',and='INTEGRAL') .OR. &
                      keyword_contains(line,'PATH',and='INTEGRALS') ) THEN
                ! Path Integrals                
                cntl%tpath=.TRUE.
                cntl%tpimd=.TRUE.
             ELSEIF ( keyword_contains(line,'PATH',and='MINIMIZATION') ) THEN
                ! Path Minimisation (Saddle point determination)
                cntl%tpath=.TRUE.
                cntl%tpmin=.TRUE.
             ELSEIF ( keyword_contains(line,'PRNGSEED') ) THEN
                READ(iunit,'(A)',iostat=ierr) line
                first=1
                CALL readsi(line,first,last,cnti%iprng,erread)
             ELSEIF ( keyword_contains(line,'LANGEVIN') ) THEN
                IF ( keyword_contains(line,'MOVECM',alias='MOVECOM')) glepar%gle_com=0
                IF ( keyword_contains(line,'CENTROIDOFF')) tglepc=.FALSE.
                IF ( keyword_contains(line,'WHITE') ) THEN
                   glepar%gle_mode=gle_white
                   glepar%gle_ns=0
                ELSEIF ( keyword_contains(line,'CPMD') ) THEN
                   glepar%gle_mode=gle_cpmd
                   glepar%gle_ns=gle_cp_ns
                ELSEIF ( keyword_contains(line,'OPTIMAL') ) THEN
                   glepar%gle_mode=gle_opt
                   glepar%gle_ns=gle_opt_ns
                ELSEIF ( keyword_contains(line,'SMART') ) THEN
                   glepar%gle_mode=gle_smart
                   glepar%gle_ns=gle_smart_ns          
                ELSEIF ( keyword_contains(line,'CUSTOM') ) THEN
                   glepar%gle_mode=gle_cust
                ENDIF
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                IF (glepar%gle_mode/=gle_cust) THEN
                   first=1
                   CALL readsr(line ,first,last,glepar%gle_omega,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ELSE
                   first=1
                   CALL readsi(line ,first,last,glepar%gle_ns,erread)
                   IF (erread .OR. glepar%gle_ns<0) THEN 
                      error_message        = "SIZE_NS FOR LANGEVIN CUSTOM NOT SPECIFIED"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ENDIF
             ELSEIF ( keyword_contains(line,'QMMM') ) THEN
                ! QMMM
                IF ( keyword_contains(line,'QMMMEASY') ) THEN
                   cntl%tqmmech=.TRUE.
                ELSE
                   cntl%tqmmm=.TRUE.
                ENDIF
             ELSEIF ( keyword_contains(line,'INTERFACE') ) THEN
                ! Interface to classical cntl%md program
                cntl%tinter=.TRUE.
                ! original EGO interface
                IF ( keyword_contains(line,'EGO') ) THEN
                   cnti%iftype=1
                   ! biswas/2005. modified EGO interface for Gromacs
                ELSEIF ( keyword_contains(line,'GMX') ) THEN
                   cnti%iftype=2
                ELSEIF ( keyword_contains(line,'IPHIGENIE') ) THEN
                   cnti%iftype=3
                ELSE
                   error_message        = "UNKNOWN MM INTERFACE TYPE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                IF ( keyword_contains(line,'MULLIKEN')) cnti%icmet=1
                IF ( keyword_contains(line,'HIRSHFELD')) cnti%icmet=2
                IF ( keyword_contains(line,'ESP')) cnti%icmet=3
                IF ( keyword_contains(line,'LOWDIN')) cnti%icmet=4
                IF ( keyword_contains(line,'PCGFIRST')) cntl%tpcgfi=.TRUE.
             ELSEIF ( keyword_contains(line,'TROTTER',and='FACTOR',cut_at='=') .OR. &
                      keyword_contains(line,'TROTTER',and='FACTORIZATION',cut_at='=') ) THEN
                ! Trotter factor for free energy functional
                IF ( keyword_contains(line,'OFF') ) THEN
                   fint1%ttrot=.FALSE.
                   fint1%tbogo=.FALSE.
                ENDIF
                fint1%ttrot=.TRUE.
                first = index_of_delimiter(line,'FACTOR','=')
                ! If equal sign is present, read advanced parameters
                IF (first > 1) THEN
                   CALL readsi(line,first,last,fint5%ntabbetap,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   IF (fint5%ntabbetap<=0.OR.fint5%ntabbetap>maxbetap) THEN
                      WRITE(output_unit,'(A,I3)')&
                           ' CONTROL| NUMBER OF TROTTER FACTORS : ',fint5%ntabbetap
                      WRITE(output_unit,'(A,I3)')&
                           ' CONTROL| BETWEEN 0 AND',maxbetaP
                      CALL stopgm(procedureN,'WRONG TROTTER FACTOR NUMBER',& 
                           __LINE__,__FILE__)
                   ENDIF
                   previous_line = line
                   DO i=1,fint5%ntabbetap
                      READ(iunit,'(A)',iostat=ierr) line
                      first=1
                      CALL readsr(line,first,last,fint5%densbetap(i),erread)
                      IF (erread) THEN
                         error_message        = "ERROR WHILE READING VALUE"
                         something_went_wrong = .true.
                         go_on_reading        = .false.
                      ENDIF
                      first=last
                      CALL readsr(line,first,last,fint5%tabbetap(i), erread)
                      IF (erread) THEN
                         error_message        = "ERROR WHILE READING VALUE"
                         something_went_wrong = .true.
                         go_on_reading        = .false.
                      ENDIF
                   ENDDO
                   fint1%betap=fint5%tabbetap(1)
                ! No equal sign? Do as described in the manual
                ELSE
                   previous_line = line
                   READ(iunit,'(A)',iostat=ierr) line
                   CALL readsr(line,1,last,fint1%betap,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ENDIF
             ELSEIF ( keyword_contains(line,'LINEAR',and='RESPONSE') ) THEN
                ! Linear response calculation   
                cntl%tresponse  = .TRUE.
                cntl%wfopt=.TRUE.
             ELSEIF ( keyword_contains(line,'HARMONIC',and='REFERENCE') ) THEN
                ! Switch on/off harmonic reference system
                IF ( keyword_contains(line,'OFF') ) THEN
                   cntl%tharm=.FALSE.
                ELSE
                   cntl%tharm=.TRUE.
                ENDIF
             ELSEIF ( keyword_contains(line,'SCALED',and='MASSES') ) THEN
                ! Switch on/off scaling of electron masses
                IF ( keyword_contains(line,'OFF') ) THEN
                   cntl%tmass=.FALSE.
                ELSE
                   cntl%tmass=.TRUE.
                ENDIF
             ELSEIF ( keyword_contains(line,'TDDFT') ) THEN
                cntl%tddft=.TRUE.
             ELSEIF ( keyword_contains(line,'LSD',alias='LSDA') .OR. &
                    keyword_contains(line,'LOCAL',and='SPIN') ) THEN
                ! Local Spin Density Approximation
                cntl%tlsd=.TRUE.
                clsd%nlsd=2
                clsd%nlsx=3
             ELSEIF ( keyword_contains(line,'SSIC') ) THEN ! cmb_ssic
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                first=1
                CALL readsr(line,first,last,cntr%asic,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                ! mb          IA=IE ! xc not yet implemented
                ! mb          CALL READSR(LINE,IA,IE,bsic,ERREAD)
                ! mb          IF(ERREAD) GOTO 99
                cntl%tlsd=.TRUE.
                clsd%nlsd=2
                clsd%nlsx=3
                cntl%tsic=.TRUE.
             ELSEIF ( keyword_contains(line,'NONORTHOGONAL',and='ORBITALS') ) THEN
                ! Use nonorthogonal orbital
                IF ( keyword_contains(line,'OFF') ) THEN
                   cntl%nonort=.FALSE.
                   cntl%cnstmd=.TRUE.
                ELSE
                   cntl%nonort=.TRUE.
                   cntl%cnstmd=.FALSE.
                   previous_line = line
                   READ(iunit,'(A)',iostat=ierr) line
                   CALL readsr(line,1,last,nort_com%slimit,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ENDIF
             ELSEIF ( keyword_contains(line,'LANCZOS',and='DIAGONALIZATION') ) THEN
                ! Lanczos diagonalization (Friesner type)
                cntl%tdiag=.TRUE.
                cntl%tlanc=.TRUE.
                IF ( keyword_contains(line,'ALL')) fint1%tfral=.TRUE.
                IF ( keyword_contains(line,'NEW')) cntl%tfrsblk=.TRUE.
                IF ( keyword_contains(line,'OPT') ) THEN
                   cntl%tdiag=.FALSE.
                   cntl%tdiagopt=.TRUE.
                ELSEIF ( keyword_contains(line,'PERIODIC',cut_at='=') ) THEN
                   cntl%tdiag=.FALSE.
                   cntl%tdiagopt=.TRUE.
                   IF (cntr%b2limit<=1.0e-16_real_8) cntr%b2limit=1.0e-8_real_8
                   first = index_of_delimiter(line,'PERIODIC','=')
                   IF (first > 1) THEN
                      CALL readsi(line,first,last,cnti%nperid,erread)
                      IF (erread) THEN
                         error_message        = "ERROR WHILE READING VALUE"
                         something_went_wrong = .true.
                         go_on_reading        = .false.
                      ENDIF
                   ELSE
                      error_message        = "PERIODIC DIAGONALIZATION: NO FREQUENCY SPECIFIED"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ELSEIF ( keyword_contains(line,'RESET',cut_at='=') ) THEN
                   cntl%tdiag=.FALSE.
                   cntl%tdiagopt=.TRUE.
                   IF (cntr%b2limit<=1.0e-16_real_8) cntr%b2limit=1.0e-8_real_8
                   first = index_of_delimiter(line,'RESET','=')
                   CALL readsi(line,first,last,cnti%nrperi,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ENDIF
             ELSEIF ( keyword_contains(line,'LANCZOS',and='PARAMETER',cut_at='=') .OR. &
                     keyword_contains(line,'LANCZOS',and='PARAMETERS',cut_at='=')) THEN
                ! Parameter for Lanczos diagonalization (Friesner type)
                IF ( keyword_contains(line,'ALL',cut_at='=')) fint1%tfral=.TRUE.
                keep_first = max(index_of_delimiter(line,'N','='), &
                                    index_of_delimiter(line,'n','='))
                IF (keep_first > 1) THEN
                   CALL readsi(line,keep_first,last,fint4%ntabtrot,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   IF (fint4%ntabtrot<=0.OR.fint4%ntabtrot>maxtrot) THEN
                      WRITE(output_unit,'(A,I3)')&
                           ' CONTROL| NUMBER OF BETATOL : ',fint4%ntabtrot
                      WRITE(output_unit,'(A,I3)')&
                           ' CONTROL| BETWEEN 1 AND',maxtroT
                      CALL stopgm('CONTROL','WRONG BETATOL NUMBER',& 
                           __LINE__,__FILE__)
                   ENDIF
                ELSE
                   fint4%ntabtrot=1
                ENDIF
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                first=1
                CALL readsi(line,first,last,cnti%n_fries,   erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                first=last
                CALL readsi(line,first,last,cnti%nkry_max,  erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                first=last
                CALL readsi(line,first,last,cnti%nkry_block,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                first=last
                CALL readsr(line,first,last,cntr%b2limit,   erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                IF (keep_first > 1) THEN
                   fint4%b2trot(1)=cntr%b2limit
                   DO i=2,fint4%ntabtrot
                      READ(iunit,'(A)',iostat=ierr) line
                      first=1
                      CALL readsr(line,first,last,fint4%denstrot(i),erread)
                      IF (erread) THEN
                         error_message        = "ERROR WHILE READING VALUE"
                         something_went_wrong = .true.
                         go_on_reading        = .false.
                      ENDIF
                      first=last
                      CALL readsr(line,first,last,fint4%b2trot(i),  erread)
                      IF (erread) THEN
                         error_message        = "ERROR WHILE READING VALUE"
                         something_went_wrong = .true.
                         go_on_reading        = .false.
                      ENDIF
                   ENDDO
                   IF (fint4%ntabtrot==1) THEN
                      fint4%denstrot(1)=1._real_8
                   ELSE
                      fint4%denstrot(1)=2._real_8*fint4%denstrot(2)
                   ENDIF
                ENDIF
             ELSEIF ( keyword_contains(line,'DAVIDSON',and='DIAGONALIZATION') ) THEN
                ! Davidson diagonalization 
                cntl%tdiag=.TRUE.
                cntl%tdavi=.TRUE.
             ELSEIF ( keyword_contains(line,'DAVIDSON',and='PARAMETER') .OR. &
                     keyword_contains(line,'DAVIDSON',and='PARAMETERS') ) THEN
                ! Parameter for Davidson diagonalization
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                first=1
                CALL readsi(line,first,last,cnti%ndavv, erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                first=last
                CALL readsr(line,first,last,cntr%epsdav, erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                first=last
                CALL readsi(line,first,last,cnti%n_fries,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF ( keyword_contains(line,'ALEXANDER',and='MIXING',but_not='ANDERSON') ) THEN
                ! Alexander mixing parameter
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                CALL readsr(line,1,last,andr2%alxmix, erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF ( keyword_contains(line,'ANDERSON',and='MIXING',but_not='ALEXANDER',cut_at='=') ) THEN
                ! Anderson mixing parameter
                IF ( keyword_contains(line,'G-SPACE')) broy1%tgmix=.TRUE.
                first = max(index_of_delimiter(line,'N','='), index_of_delimiter(line,'n','='))
                IF (first > 1) THEN
                   CALL readsi(line,first,last,andr2%ntabmix,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   IF (andr2%ntabmix<=0.OR.andr2%ntabmix>maxmix) THEN
                      WRITE(output_unit,'(A,I3)')&
                           ' CONTROL| NUMBER OF MIXING : ',andr2%ntabmix
                      WRITE(output_unit,'(A,I3)')&
                           ' CONTROL| BETWEEN 0 AND',maxmiX
                      CALL stopgm('CONTROL','WRONG MIXING NUMBER',& 
                           __LINE__,__FILE__)
                   ENDIF
                   previous_line = line
                   DO i=1,andr2%ntabmix
                      READ(iunit,'(A)',iostat=ierr) line
                      first=1
                      CALL readsr(line,first,last,andr2%densmix(i),erread)
                      IF (erread) THEN
                         error_message        = "ERROR WHILE READING VALUE"
                         something_went_wrong = .true.
                         go_on_reading        = .false.
                      ENDIF
                      first=last
                      CALL readsr(line,first,last,andr2%andmix(i), erread)
                      IF (erread) THEN
                         error_message        = "ERROR WHILE READING VALUE"
                         something_went_wrong = .true.
                         go_on_reading        = .false.
                      ENDIF
                   ENDDO
                   andr2%andrmix=andr2%densmix(1)
                ELSE
                   previous_line = line
                   READ(iunit,'(A)',iostat=ierr) line
                   CALL readsr(line,1,last,andr2%andrmix,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ENDIF
             ELSEIF ( keyword_contains(line,'BROYDEN',and='MIXING') ) THEN
                ! Broyden mixing parameter
                broy1%tgmix=.TRUE.
                broy1%tgbroy=.TRUE.
                IF ( keyword_contains(line,'ADAPT',cut_at='=')) broy1%tadbroy=.TRUE.
                !
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                IF (ierr /=0) THEN
                   error_message        = "ERROR WHILE READING LINE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                !
                ! Special mixing parameters using = or standard format?
                !
                IF (INDEX(line,'=') /= 0) THEN
                   IF ( keyword_contains(line,'MIX',alias='BROYMIX',cut_at='=') ) THEN
                      first = max(index_of_delimiter(line,'MIX','='), &
                                  index_of_delimiter(line,'BROYMIX','='))
                      CALL readsr(line,first,last,broy1%broymix,erread)
                      IF (erread) THEN
                         error_message        = "ERROR WHILE READING VALUE"
                         something_went_wrong = .true.
                         go_on_reading        = .false.
                      ENDIF
                   ENDIF
                   IF ( keyword_contains(line,'CUTBROY',cut_at='=') ) THEN
                      first = index_of_delimiter(line,'CUTBROY','=')
                      CALL readsr(line,first,last,broy1%ecutbroy,erread)
                      IF (erread) THEN
                         error_message        = "ERROR WHILE READING VALUE"
                         something_went_wrong = .true.
                         go_on_reading        = .false.
                      ENDIF
                   ENDIF
                   IF ( keyword_contains(line,'W02BROY',alias='W02',cut_at='=') ) THEN
                      first = max(index_of_delimiter(line,'W02BROY','='), &
                                  index_of_delimiter(line,'W02','='))
                      CALL readsr(line,first,last,broy1%w02broy,erread)
                      IF (erread) THEN
                         error_message        = "ERROR WHILE READING VALUE"
                         something_went_wrong = .true.
                         go_on_reading        = .false.
                      ENDIF
                   ENDIF
                   IF ( keyword_contains(line,'NFRBROY',cut_at='=') ) THEN
                      first = index_of_delimiter(line,'NFRBROY','=')
                      CALL readsi(line,first,last,broy1%nfrbroy,erread)
                      IF (erread) THEN
                         error_message        = "ERROR WHILE READING VALUE"
                         something_went_wrong = .true.
                         go_on_reading        = .false.
                      ENDIF
                   ENDIF
                   IF ( keyword_contains(line,'RESET',alias='IBRESET',cut_at='=') ) THEN
                      first = max(index_of_delimiter(line,'RESET','='), &
                                  index_of_delimiter(line,'IBRESET','='))
                      CALL readsi(line,first,last,broy1%ibreset,erread)
                      IF (erread) THEN
                         error_message        = "ERROR WHILE READING VALUE"
                         something_went_wrong = .true.
                         go_on_reading        = .false.
                      ENDIF
                   ENDIF
                   IF ( keyword_contains(line,'KERMIX',cut_at='=') ) THEN
                      first = index_of_delimiter(line,'KERMIX','=')
                      CALL readsr(line,first,last,broy1%kermix,erread)
                      IF (erread) THEN
                         error_message        = "ERROR WHILE READING VALUE"
                         something_went_wrong = .true.
                         go_on_reading        = .false.
                      ENDIF
                   ENDIF
                ELSE
                   ! Use default
                   READ(line ,*,iostat=ierr)&
                        broy1%broymix,broy1%ecutbroy,broy1%w02broy,broy1%nfrbroy,broy1%ibreset,broy1%kermix
                   first=1
                   CALL readsr(line,first,last,broy1%broymix, erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   first=last
                   CALL readsr(line,first,last,broy1%ecutbroy,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   first=last
                   CALL readsr(line,first,last,broy1%w02broy, erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   first=last
                   CALL readsi(line,first,last,broy1%nfrbroy, erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   first=last
                   CALL readsi(line,first,last,broy1%ibreset, erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   IF (broy1%ibreset==0) broy1%ibreset=50
                   first=last
                   CALL readsr(line,first,last,broy1%kermix,  erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ENDIF
             ELSEIF ( keyword_contains(line,'DIIS',and='MIXING',cut_at='=') ) THEN 
                ! cntl%diis mixing for density  
                first = max(index_of_delimiter(line,'N','='), index_of_delimiter(line,'n','='))
                IF (first > 1) THEN
                   CALL readsi(line,first,last,andr3%ntabnrdiis,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   IF (andr3%ntabnrdiis<=0.OR.andr3%ntabnrdiis>maxmix) THEN
                      WRITE(output_unit,'(A,I3)')&
                           ' CONTROL| NUMBER OF DIIS MIXING : ',andr3%ntabnrdiis
                      WRITE(output_unit,'(A,I3)')&
                           ' CONTROL| BETWEEN 0 AND',maxmiX
                      CALL stopgm('CONTROL','WRONG DIIS MIXING NUMBER',& 
                           __LINE__,__FILE__)
                   ENDIF
                   counter = 0
                   previous_line = line
                   DO i=1,andr3%ntabnrdiis
                      READ(iunit,'(A)',iostat=ierr) line
                      first=1
                      CALL readsr(line,first,last,andr3%densnrdiis(i),erread)
                      IF (erread) THEN
                         error_message        = "ERROR WHILE READING VALUE"
                         something_went_wrong = .true.
                         go_on_reading        = .false.
                      ENDIF
                      first=last
                      CALL readsi(line,first,last,andr3%tabnrdiis(i), erread)
                      IF (erread) THEN
                         error_message        = "ERROR WHILE READING VALUE"
                         something_went_wrong = .true.
                         go_on_reading        = .false.
                      ENDIF
                      IF (andr3%tabnrdiis(i)<counter) CALL stopgm('CONTROL',&
                           'THE DIIS MIXING NUMBERS HAVE TO BE INCREASED',& 
                           __LINE__,__FILE__)
                      counter=andr3%tabnrdiis(i)
                   ENDDO
                   andr3%nrdiismax=andr3%tabnrdiis(andr3%ntabnrdiis)
                ELSE
                   previous_line = line
                   READ(iunit,'(A)',iostat=ierr) line
                   CALL readsi(line ,1,last,andr3%nrdiis,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   andr3%nrdiismax=andr3%nrdiis
                   andr3%tabnrdiis(1)=andr3%nrdiis
                ENDIF
             ELSEIF&
                  (keyword_contains(line,'MOVERHO') ) THEN
                ! Moverho 
                tmovr=.TRUE.
                andr2%alxmix=0._real_8
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                CALL readsr(line,1,last,dmovrmix,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                ! Wavefunction extrapolation (polynomial or aspc) for BOMD.
             ELSEIF ( keyword_contains(line,'EXTRAPOLATE',and='WFN') ) THEN
                cntl%textrap=.TRUE.
                IF ( keyword_contains(line,'STORE')) cntl%tstrxtp=.TRUE.
                IF ( keyword_contains(line,'POLY',alias='POLYNOMIAL')) cntl%taspc=.FALSE.
                IF ( keyword_contains(line,'ASPC')) cntl%taspc=.TRUE.
                IF ( keyword_contains(line,'CSTEPS',cut_at='=') ) THEN
                   first = index_of_delimiter(line,'CSTEPS','=')
                   CALL readsi(line,first,last,cnti%naspc,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ENDIF
                READ(iunit,'(A)',iostat=ierr) line
                CALL readsi(line ,1,last,cnti%mextra,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF ( keyword_contains(line,'EXTRAPOLATE',and='CONSTRAINT') ) THEN
                cdftlog%tpred=.TRUE.
                READ(iunit,'(A)',iostat=ierr) line
                CALL readsi(line ,1,last,cdftpred%predord,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF ( keyword_contains(line,'TSDE') ) THEN
                ! Electronic dynamics by Steepest descent
                cntl%tsde=.TRUE.
                cntl%prec=.TRUE.
                IF ( keyword_contains(line,'NOPREC',alias='NOPRECONDITIONING')) cntl%prec=.FALSE.
             ELSEIF ( keyword_contains(line,'TSDP') ) THEN
                ! Ionic dynamics by Steepest descent
                cntl%tsdp=.TRUE.
                ! Steepest descent with line minimisation
                IF ( keyword_contains(line,'LINE')) sdpl%tsdpline=.TRUE.
                ! Ionic minimization by Conjugate Gradient
             ELSEIF ( keyword_contains(line,'TCGP') ) THEN
                sdpl%tcgp=.TRUE.
                cntl%tsdp=.TRUE.
             ELSEIF ( keyword_contains(line,'TSDC') ) THEN
                ! Cell dynamics by Steepest descent
                cntl%tsdc=.TRUE.
             ELSEIF ( keyword_contains(line,'STEEPEST',and='DESCENT') ) THEN
                IF ( keyword_contains(line,'ELECTRONS',alias='ELECTRON') ) THEN
                   cntl%tsde=.TRUE.
                   cntl%prec=.TRUE.
                   IF ( keyword_contains(line,'NOPREC',alias='NOPRECONDITIONING')) cntl%prec=.FALSE.
                ENDIF
                IF ( keyword_contains(line,'IONS',alias='ION') ) THEN
                   cntl%tsdp=.TRUE.
                   IF ( keyword_contains(line,'LINE')) sdpl%tsdpline=.TRUE.
                ENDIF
                IF ( keyword_contains(line,'CELL')) cntl%tsdc=.TRUE.
                ! Conjugate Gradient
                ! PCG (electrons-only keyword alias)
             ELSEIF ( keyword_contains(line,'CONJUGATE',and='GRADIENT',alias='PCG') ) THEN
                IF ( keyword_contains(line,'ELECTRONS',alias='PCG') ) THEN
                   cntl%pcg=.TRUE.
                   cntl%prec=.TRUE.
                   IF ( keyword_contains(line,'NOPREC',alias='NOPRECONDITIONING')) cntl%prec=.FALSE.
                   IF ( keyword_contains(line,'MINIMIZE')) cntl%pcgmin=.TRUE.
                ENDIF
                IF ( keyword_contains(line,'IONS',alias='ION',but_not='PCG') ) THEN
                   sdpl%tcgp=.TRUE.
                   cntl%tsdp=.TRUE.
                ENDIF
             ELSEIF ( keyword_contains(line,'ODIIS') ) THEN
                ! Wavefunction optimization by cntl%diis
                cntl%diis=.TRUE.
                cntl%prec=.TRUE.
                IF ( keyword_contains(line,'NOPREC',alias='NOPRECONDITIONING')) cntl%prec=.FALSE.
                IF ( keyword_contains(line,'NO_RES',alias='NO_RESET',cut_at='=') ) THEN
                   first = max(index_of_delimiter(line,'NO_RESET','='), &
                                  index_of_delimiter(line,'NO_RES','='))
                   CALL readsi(line,first,last,cnti%nreset,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ENDIF
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                CALL readsi(line,1,last,cnti%mdiis,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF ( keyword_contains(line,'HAMILTONIAN',and='CUTOFF',alias='HTHRS') ) THEN
                ! Threshold for Hessian custore_or_not in preconditioning matrix for cntl%diis
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                CALL readsr(line ,1,last,cntr%hthrs,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF ( keyword_contains(line,'GDIIS') ) THEN
                ! Geometry optimization by cntl%gdiis/cntl%bfgs
                cntl%gdiis=.TRUE.
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                CALL readsi(line ,1,last,cnti%mgdiis,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF ( keyword_contains(line,'LBFGS',alias='L-BFGS') ) THEN
                ! Geometry optimization by low-memory cntl%bfgs
                cntl%lbfgs=.TRUE.
                IF ( keyword_contains(line,'NREM') ) THEN
                   READ(iunit,'(A)',iostat=ierr) line
                   CALL readsi(line,1,last,nrem,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ELSEIF ( keyword_contains(line,'NTRUST') ) THEN
                   READ(iunit,'(A)',iostat=ierr) line
                   CALL readsi(line,1,last,ntrust,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ELSEIF ( keyword_contains(line,'NRESTT') ) THEN
                   READ(iunit,'(A)',iostat=ierr) line
                   CALL readsi(line,1,last,nrestt,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ELSEIF ( keyword_contains(line,'NTRSTR') ) THEN
                   READ(iunit,'(A)',iostat=ierr) line
                   CALL readsi(line,1,last,ntrstr,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ELSEIF ( keyword_contains(line,'TRUSTR') ) THEN
                   READ(iunit,'(A)',iostat=ierr) line
                   CALL readsr(line,1,last,trustr,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ENDIF
             ELSEIF ( keyword_contains(line,'PRFO') ) THEN
                ! Transition state search by partitioned cntl%rfo
                cntl%prfo=.TRUE.
                IF ( keyword_contains(line,'NVAR') ) THEN
                   READ(iunit,'(A)',iostat=ierr) line
                   CALL readsi(line,1,last,nvar,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ELSEIF ( keyword_contains(line,'MODE') ) THEN
                   READ(iunit,'(A)',iostat=ierr) line
                   CALL readsi(line,1,last,mode,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ELSEIF ( keyword_contains(line,'MDLOCK') ) THEN
                   READ(iunit,'(A)',iostat=ierr) line
                   CALL readsi(line,1,last,modelk,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ELSEIF ( keyword_contains(line,'TOLENV') ) THEN
                   READ(iunit,'(A)',iostat=ierr) line
                   CALL readsr(line,1,last,tolenv,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ELSEIF ( keyword_contains(line,'TRUSTP') ) THEN
                   READ(iunit,'(A)',iostat=ierr) line
                   CALL readsr(line,1,last,trustp,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ELSEIF ( keyword_contains(line,'OMIN') ) THEN
                   READ(iunit,'(A)',iostat=ierr) line
                   CALL readsr(line,1,last,omin,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ELSEIF ( keyword_contains(line,'NSVIB') ) THEN
                   READ(iunit,'(A)',iostat=ierr) line
                   CALL readsi(line,1,last,nsvib,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ELSEIF ( keyword_contains(line,'CORE',cut_at='=') ) THEN
                   first = index_of_delimiter(line,'CORE','=')
                   CALL readsi(line,first,last,ncore,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   lnomap = .FALSE.
                   ALLOCATE(icore(ncore),STAT=ierr)
                   IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                        __LINE__,__FILE__)
                   READ(iunit,'(A)',iostat=ierr) line
                   first=1
                   DO i = 1,ncore
                      CALL readsi(line,first,last,icore(i),erread)
                      IF (erread) THEN
                         error_message        = "ERROR WHILE READING VALUE"
                         something_went_wrong = .true.
                         go_on_reading        = .false.
                      ENDIF
                      first=last
                   ENDDO
                   IF (nvar==0) nvar = 3*ncore
                ELSEIF ( keyword_contains(line,'NSMAXP') ) THEN
                   READ(iunit,'(A)',iostat=ierr) line
                   CALL readsi(line,1,last,nsmaxp,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ELSEIF ( keyword_contains(line,'PRJHES') ) THEN
                   lprhes=.TRUE.
                ELSEIF ( keyword_contains(line,'DISPL') ) THEN
                   READ(iunit,'(A)',iostat=ierr) line
                   CALL readsr(line,1,last,eps_h,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ELSEIF ( keyword_contains(line,'HESSTYPE') ) THEN
                   READ(iunit,'(A)',iostat=ierr) line
                   CALL readsi(line,1,last,m_hess,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ENDIF
             ELSEIF ( keyword_contains(line,'HESSCORE') ) THEN
                ! Relaxation of environment and calculation of partial Hessian
                cntl%prfo=.TRUE.
                nsmaxp=0
             ELSEIF ( keyword_contains(line,'BFGS') ) THEN
                ! Geometry optimization by quasi-Newton update
                cntl%bfgs=.TRUE.
             ELSEIF ( keyword_contains(line,'RFO') ) THEN
                ! Geometry optimization by Rational Function Approximation
                cntl%rfo=.TRUE.
                IF ( keyword_contains(line,'ORDER',cut_at='=') ) THEN
                   first = index_of_delimiter(line,'ORDER','=')
                   CALL readsi(line,first,last,cnti%nsorder,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ENDIF
             ELSEIF ( keyword_contains(line,'INR',alias='TURBOCICCIO',cut_at='=') .AND. &
                     keyword_contains(line,'PARAMETER',alias='PARAMETERS',cut_at='=') ) THEN
                ! Param. for INR geometry optim.: nreg; gnx_inr(nreg); tolx_inr(nreg) 
                first = max(index_of_delimiter(line,'N','='), index_of_delimiter(line,'n','='))
                IF (first > 1) THEN
                   CALL readsi(line,first,last,inr_integer%nreg,erread)
                   IF (inr_integer%nreg<=0.OR.inr_integer%nreg>maxreg) THEN
                      WRITE(output_unit,'(A,I3)')&
                           ' CONTROL| NUMBER OF GRADIENT REGIONS : ',inr_integer%nreg
                      WRITE(output_unit,'(A,I3)')&
                           ' CONTROL| BETWEEN 1 AND',maxreG
                      CALL stopgm('CONTROL','WRONG NREG NUMBER',& 
                           __LINE__,__FILE__)
                   ENDIF
                ENDIF
                DO i=1,inr_integer%nreg
                   READ(iunit,'(A)',iostat=ierr) line
                   CALL readsr(line,1,last,gnx_inr(i),erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   CALL readsr(line,last+1,last,tolx_inr(i),erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ENDDO
             ELSEIF ( keyword_contains(line,'IMPLICIT',and='NEWTON',alias='NEWTON-RAPHSON') .OR. &
                     keyword_contains(line,'TURBOCICCIO') ) THEN
                cntl%tinr=.TRUE.
                ! Geometry optimization by Implicit Newton-Raphson
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                CALL readsi(line,1,last,inr_integer%itmax_inr,erread)
                IF ( keyword_contains(previous_line,'CONTINUE',alias='CONT')) inr_logical%inr_cont=.TRUE.
                IF ( keyword_contains(previous_line,'PRECONDITIONED',alias='PREC') ) THEN
                   inr_logical%inr_prec=.TRUE.
                ENDIF
                IF ( keyword_contains(previous_line,'ALTERNATIVE',and='STEP')) inr_logical%inr_step=.TRUE.
                IF ( keyword_contains(previous_line,'VERBOSE')) inr_logical%inr_verbose=.TRUE.
             ELSEIF ( keyword_contains(line,'MIXSD') ) THEN
                inr_logical%tmixsd=.TRUE.
                READ(iunit,'(A)',iostat=ierr) line
                CALL readsr(line,1,last,rmixsd,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF ( keyword_contains(line,'MIXDIIS') ) THEN
                inr_logical%tmixgdiis=.TRUE.
                READ(iunit,'(A)',iostat=ierr) line
                CALL readsr(line,1,last,rmixsd,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF ( keyword_contains(line,'RESTART') ) THEN
                ! Restart
                restart1%restart=.TRUE.
                IF ( keyword_contains(line,'WAVEFUNCTION')) restart1%rwf  =.TRUE.
                IF ( keyword_contains(line,'COORDINATES')) restart1%rco  =.TRUE.
                IF ( keyword_contains(line,'VELOCITIES')) restart1%rvel =.TRUE.
                IF ( keyword_contains(line,'ACCUMULATORS')) restart1%rac  =.TRUE.
                IF ( keyword_contains(line,'HESSIAN')) restart1%rhe  =.TRUE.
                IF ( keyword_contains(line,'NOSEE')) restart1%rnoe =.TRUE.
                IF ( keyword_contains(line,'NOSEP')) restart1%rnop =.TRUE.
                IF ( keyword_contains(line,'NOSEC')) restart1%rnoc =.TRUE.
                IF ( keyword_contains(line,'GEOFILE')) restart1%rgeo =.TRUE.
                IF ( keyword_contains(line,'LATEST')) restart1%rlate=.TRUE.
                IF ( keyword_contains(line,'VIBANALYSIS',alias='VIBANA')) restart1%rvib =.TRUE.
                IF ( keyword_contains(line,'CELL')) restart1%rcell=.TRUE.
                IF ( keyword_contains(line,'POTENTIAL')) restart1%rpot =.TRUE.
                IF ( keyword_contains(line,'DENSITY')) restart1%rrho =.TRUE.
                IF ( keyword_contains(line,'KPOINTS')) restart1%rkpt =.TRUE.
                IF ( keyword_contains(line,'OCCUPATION')) restart1%rocc =.TRUE.
                IF ( keyword_contains(line,'ORBHA',alias='HARDNESS')) restart1%roh  =.TRUE.
                IF ( keyword_contains(line,'LINRES')) restart1%rlr  =.TRUE.
                IF ( keyword_contains(line,'PHES')) restart1%rphes=.TRUE.
                IF ( keyword_contains(line,'LSSTAT')) restart1%rlscl=.TRUE.
                IF ( keyword_contains(line,'ADPTTL')) restart1%rdtl =.TRUE.
                IF ( keyword_contains(line,'CONSTRAINTS')) restart1%rcon =.TRUE.
                IF ( keyword_contains(line,'EXTRAP',alias='EXTRAPOLATION')) restart1%rxtp =.TRUE.
                IF ( keyword_contains(line,'NOREST',alias='NORESTART')) restart1%rnon =.TRUE.
                IF ( keyword_contains(line,'CVDFT')) cdftlog%rcdft=.TRUE.
                IF ( keyword_contains(line,'PRNG')) restart1%rprng=.TRUE.
                IF ( keyword_contains(line,'GLE')) restart1%rgle= .TRUE.
     
                IF ( keyword_contains(line,'ALL') ) THEN
                   restart1%rwf  =.TRUE.
                   restart1%rco  =.TRUE.
                   restart1%rvel =.TRUE.
                   restart1%rac  =.TRUE.
                   restart1%rhe  =.TRUE.
                   restart1%rnoe =.TRUE.
                   restart1%rnop =.TRUE.
                   restart1%rnoc =.TRUE.
                   restart1%rgeo =.TRUE.
                   restart1%rlate=.TRUE.
                   restart1%rvib =.TRUE.
                   restart1%rcell=.TRUE.
                   restart1%rpot =.TRUE.
                   restart1%rrho =.TRUE.
                   restart1%rkpt =.TRUE.
                   restart1%rocc =.TRUE.
                   restart1%roh  =.TRUE.
                   restart1%rlr  =.TRUE.
                   restart1%rcon =.TRUE.
                   cdftlog%rcdft=.TRUE.
                   restart1%rxtp =.TRUE.
                   restart1%rprng=.TRUE.
                   restart1%rgle =.TRUE.
                ENDIF
             ELSEIF ( keyword_contains(line,'INTFILE') ) THEN
                IF ( keyword_contains(line,'READ')) iface1%intread=.TRUE.
                IF ( keyword_contains(line,'WRITE')) iface1%intwrite=.TRUE.
                IF ( keyword_contains(line,'FILENAME') ) THEN
                   READ(iunit,'(A)',iostat=ierr) intfn
                   IF (ierr /= 0) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ENDIF
             ELSEIF ( keyword_contains(line,'MAXSTEP',alias='MAXSTEPS') ) THEN
                ! Maximum number of steps
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                CALL readsi(line ,1,last,cnti%nomore,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF ( keyword_contains(line,'MAXIT',alias='MAXITER') ) THEN
                ! Maximum number of iterations (for diag or cntl%md BO)
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                CALL readsi(line ,1,last,cnti%nomore_iter,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF ( keyword_contains(line,'MAXRUN',alias='MAXRUNTIME') ) THEN
                ! Maximum time for the job
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                CALL readsr(line ,1,last,var,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                IF (tjlimit<=0._real_8) THEN
                   tjlimit=var
                ELSE
                   tjlimit=MIN(tjlimit,var)
                ENDIF
             ELSEIF ( keyword_contains(line,'INITIALIZE',and='WAVEFUNCTION') .OR. &
                      keyword_contains(line,'INITIAL',and='WAVEFUNCTION') ) THEN
                ! Wavefunction initialization
                IF ( keyword_contains(line,'RANDOM',alias='RANDOMIZED')) cnti%inwfun=1
                IF ( keyword_contains(line,'ATOMIC',alias='ATOMS') ) THEN
                   cnti%inwfun=2
                   IF ( keyword_contains(line,'PRIMITIVE')) cntl%tipao=.TRUE.
                ENDIF
                IF ( keyword_contains(line,'SIMPLE')) cnti%inwfun=3
             ELSEIF ( keyword_contains(line,'RATTLE') ) THEN
                ! Parameters for the RATTLE step
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                first=1
                CALL readsi(line ,first,last,cnti%maxit,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                first=last
                CALL readsr(line ,first,last,cntr%epsog,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF ( keyword_contains(line,'ORTHOGONALIZATION',alias='ORTHOGONALIZE') ) THEN
                ! Method for orthogonalization  
                IF ( keyword_contains(line,'GRAM',and='SCHMIDT',alias='GRAM-SCHMIDT')) cntl%tlowd=.FALSE.
                IF ( keyword_contains(line,'LOWDIN') ) THEN
                   cntl%tlowd=.TRUE.
                   IF ( keyword_contains(line,'MATRIX'))&
                        cntl%tlowdmat=.TRUE.
                ENDIF
             ELSEIF ( keyword_contains(line,'QUENCH') ) THEN
                ! Quench velocities
                IF ( keyword_contains(line,'IONS')) cntl%quenchp=.TRUE.
                IF ( keyword_contains(line,'ELECTRONS')) cntl%quenche=.TRUE.
                IF ( keyword_contains(line,'BO')) cntl%quenchb=.TRUE.
                IF ( keyword_contains(line,'CELL')) cntl%quenchc=.TRUE.
             ELSEIF ( keyword_contains(line,'RANDOMIZE') ) THEN
                ! Randomize ionic positions and/or wavefunction
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                IF ( keyword_contains(previous_line,'WAVEFUNCTION') ) THEN
                   cntl%trane=.TRUE.
                   CALL readsr(line,1,last,cntr%ampre, erread)
                ELSEIF ( keyword_contains(previous_line,'COORDINATES') ) THEN
                   cntl%tranp=.TRUE.
                   CALL readsr(line,1,last,cntr%amprp, erread)
                ELSEIF ( keyword_contains(previous_line,'CELL') ) THEN
                   cntl%tranc=.TRUE.
                   CALL readsr(line,1,last,cntr%amprc, erread)
                ELSEIF ( keyword_contains(previous_line,'DENSITY') ) THEN
                   andr2%trand=.TRUE.
                   CALL readsr(line,1,last,andr2%amprd, erread)
                ELSE
                   error_message        = "UNKNOWN OPTION FOR RANDOMIZATION"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF ( keyword_contains(line,'USE_MTS') ) THEN
                cntl%use_mts=.true.
             ELSEIF ( keyword_contains(line,'TIMESTEP') ) THEN
                ! Time step for electrons and ions
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                IF (.not. keyword_contains(previous_line,'IONS',alias='ION') .and. &
                    .not. keyword_contains(previous_line,'ELECTRONS',alias='ELECTRON') ) THEN
                   CALL readsr(line ,1,last,cntr%delt_elec,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   cntr%delt_ions = cntr%delt_elec
                ELSEIF ( keyword_contains(previous_line,'IONS',alias='ION') .and. &
                         keyword_contains(previous_line,'ELECTRONS',alias='ELECTRON') ) THEN
                   CALL readsr(line ,1,last,cntr%delt_elec,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   cntr%delt_ions = cntr%delt_elec
                ELSEIF ( keyword_contains(previous_line,'ELECTRONS',alias='ELECTRON') ) THEN
                   CALL readsr(line ,1,last,cntr%delt_elec,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ELSEIF ( keyword_contains(previous_line,'IONS',alias='ION') ) THEN
                   CALL readsr(line ,1,last,cntr%delt_ions,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ENDIF
             ELSEIF ( keyword_contains(line,'EMASS') ) THEN
                ! Electron mass
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                CALL readsr(line ,1,last,cntr%emass,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF ( keyword_contains(line,'CMASS') ) THEN
                ! Cell mass
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                CALL readsr(line ,1,last,cntr%cmass,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                !NABDY[
             ELSEIF ( keyword_contains(line,'NABDY_ZMAX') ) THEN
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                CALL readsi(line ,1,last,nabdyvar%nabdy_zmax,erread)
                IF(erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                !NABDY]
             ELSEIF ( keyword_contains(line,'TEMPERATURE') ) THEN
                ! Temperature 
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                first=1
                IF ( keyword_contains(previous_line,'ELECTRONS',alias='ELECTRON') ) THEN
                   CALL readsr(line,first,last,fint1%betael,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ELSE
                   CALL readsr(line,first,last,cntr%tempw,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   first=last
                   IF ( keyword_contains(previous_line,'RAMP') ) THEN
                      CALL readsr(line,first,last,cntr%trampt,erread)
                      IF (erread) THEN
                         error_message        = "ERROR WHILE READING VALUE"
                         something_went_wrong = .true.
                         go_on_reading        = .false.
                      ENDIF
                      first=last
                      CALL readsr(line,first,last,cntr%trampr,erread)
                      IF (erread) THEN
                         error_message        = "ERROR WHILE READING VALUE"
                         something_went_wrong = .true.
                         go_on_reading        = .false.
                      ENDIF
                   ELSE
                      cntr%trampt=cntr%tempw
                      cntr%trampr=0.0_real_8
                   ENDIF
                ENDIF
             ELSEIF ( keyword_contains(line,'RESCALE',and='VELOCITIES') ) THEN
                ! Re-adjust ionic velocities only after restart to
                ! desired temperature TEMPW
                cntl%trescale=.TRUE.
             ELSEIF ( keyword_contains(line,'SUBTRACT') ) THEN
                !    Subtract the center of mass velocity each NCOMV steps
                !    and/or subtract the rotational velocity each NROTV steps
                IF ((keyword_contains(line,'COMVEL')).AND.&
                     (keyword_contains(line,'ROTVEL')) ) THEN
                   !    do both
                   comvl%tsubcom=.TRUE.
                   comvl%tsubrot=.TRUE.
                   READ(iunit,'(A)',iostat=ierr) line
                   CALL readsi(line,1,last,comvl%ncomv,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   comvl%nrotv=comvl%ncomv
                ELSEIF ( keyword_contains(line,'COMVEL') ) THEN
                   !    only COM
                   comvl%tsubcom=.TRUE.
                   READ(iunit,'(A)',iostat=ierr) line
                   CALL readsi(line,1,last,comvl%ncomv,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ELSEIF ( keyword_contains(line,'ROTVEL') ) THEN
                   comvl%tsubrot=.TRUE.
                   !    only ROT
                   READ(iunit,'(A)',iostat=ierr) line
                   CALL readsi(line,1,last,comvl%nrotv,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ENDIF
             ELSEIF ( keyword_contains(line,'TEMPCONTROL') ) THEN
                ! Temperature control
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                IF ( keyword_contains(previous_line,'IONS') ) THEN
                   ! For the ions
                   cntl%tcp=.TRUE.
                   first=1
                   CALL readsr(line,first,last,cntr%tempw,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   first=last
                   CALL readsr(line,first,last,cntr%tolp, erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ELSEIF ( keyword_contains(previous_line,'ELECTRONS') ) THEN
                   ! For the electrons
                   cntl%tc=.TRUE.
                   first=1
                   CALL readsr(line,first,last,cntr%ekinw,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   first=last
                   CALL readsr(line,first,last,cntr%toll, erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ELSEIF ( keyword_contains(previous_line,'CELL') ) THEN
                   ! For the cell
                   cntl%tcc=.TRUE.
                   first=1
                   CALL readsr(line,first,last,cntr%ekinhr,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   first=last
                   CALL readsr(line,first,last,cntr%tolc, erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ELSE
                   CALL stopgm('CONTROL','UNKNOWN TEMPERATURE CONTROL',& 
                        __LINE__,__FILE__)
                ENDIF
                !CNABDY[
             ELSEIF ( keyword_contains(line,'NABDY_SOFT') ) THEN
                nabdyvar%naforce_screened=.TRUE.
                READ(iunit,'(A)',iostat=ierr) line
                first=1
                CALL readsr(line,first,last,nabdyvar%nasoftening,erread)
                IF(erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF ( keyword_contains(line,'NABDY_REDISTR_AMPLI') ) THEN
                nabdyvar%tnarepart=.TRUE.
             ELSEIF ( keyword_contains(line,'NABDY_SCALEP') ) THEN
                nabdyvar%scalep=.TRUE.
             ELSEIF ( keyword_contains(line,'NABDY_THERMO') ) THEN
                nabdyvar%tnafric=.TRUE.
                READ(iunit,'(A)',iostat=ierr) line
                first=1
                CALL readsr(line,first,last,nabdyfric%natempcm,erread)
                IF(erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                first=last
                CALL readsr(line,first,last,nabdyfric%natempbp,erread)
                IF(erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                !CNABDY]
             ELSEIF ( keyword_contains(line,'BERENDSEN') ) THEN
                ! Temperature control
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                IF ( keyword_contains(previous_line,'IONS') ) THEN
                   ! For the ions
                   cntl%tberp=.TRUE.
                   first=1
                   CALL readsr(line,first,last,cntr%tempw,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   first=last
                   CALL readsr(line,first,last,cntr%taubp, erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ELSEIF ( keyword_contains(previous_line,'ELECTRONS') ) THEN
                   ! For the electrons
                   cntl%tbere=.TRUE.
                   first=1
                   CALL readsr(line,first,last,cntr%ekinw,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   first=last
                   CALL readsr(line,first,last,cntr%taube, erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ELSEIF ( keyword_contains(previous_line,'CELL') ) THEN
                   ! For the cell
                   cntl%tberc=.TRUE.
                   first=1
                   CALL readsr(line,first,last,cntr%ekinhr,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   first=last
                   CALL readsr(line,first,last,cntr%taubc, erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ELSE
                   error_message        = "UNKNOWN BERENDSEN FLAG"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF ( keyword_contains(line,'NOSE',and='IONS') ) THEN
                ! Nose thermostats for ions
                cntl%tnosep=.TRUE.
                IF ( keyword_contains(line,'ULTRA')) nosl%tultra=.TRUE.
                IF ( keyword_contains(line,'MASSIVE')) nosl%tmnose=.TRUE.
                IF ( keyword_contains(line,'CENTROIDOFF')) tnosepc=.FALSE. 
                IF ( keyword_contains(line,'CAFES') ) THEN
                   tcafes=.TRUE.
                   nosl%tmnose=.TRUE.
                   nosl%tultra=.FALSE.
                   READ(iunit,'(A)',iostat=ierr) line
                   first=1
                   CALL readsi(line,first,last,ncafesgrp,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   IF (ncafesgrp<=1) THEN
                      WRITE(output_unit,*) 'CONTROL| WARNING! NCAFESGRP SHOULD BE > 1'
                   ENDIF
                   ALLOCATE(cafesini(2,ncafesgrp),STAT=ierr)
                   IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                        __LINE__,__FILE__)
                   ALLOCATE(cafesinr(2,ncafesgrp),STAT=ierr)
                   IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                        __LINE__,__FILE__)
                   DO i=1,ncafesgrp
                      READ(iunit,'(A)',iostat=ierr) line
                      first=1
                      CALL readsi(line,first,last,cafesini(1,i),erread)
                      IF (erread) THEN
                         error_message        = "ERROR WHILE READING VALUE"
                         something_went_wrong = .true.
                         go_on_reading        = .false.
                      ENDIF
                      first=last
                      CALL readsi(line,first,last,cafesini(2,i),erread)
                      IF (erread) THEN
                         error_message        = "ERROR WHILE READING VALUE"
                         something_went_wrong = .true.
                         go_on_reading        = .false.
                      ENDIF
                      first=last
                      CALL readsr(line,first,last,cafesinr(1,i),erread)
                      IF (erread) THEN
                         error_message        = "ERROR WHILE READING VALUE"
                         something_went_wrong = .true.
                         go_on_reading        = .false.
                      ENDIF
                      first=last
                      CALL readsr(line,first,last,cafesinr(2,i),erread)
                      IF (erread) THEN
                         error_message        = "ERROR WHILE READING VALUE"
                         something_went_wrong = .true.
                         go_on_reading        = .false.
                      ENDIF
                   ENDDO
                   DO first=1,ions1%nat
                      iflag=0
                      DO last=1,ncafesgrp
                         IF ((cafesini(1,last)<=first)&
                              .AND. (cafesini(2,last)>=first) ) THEN
                            iflag=iflag+1
                         ENDIF
                      ENDDO
                      IF (iflag/=1) THEN
                         WRITE(output_unit,*) 'ATOM ',first,' IS LISTED IN ',&
                              iflag,'CAFES GROUPS'
                         CALL stopgm('CONTROL','CAFES INPUT ERROR',& 
                              __LINE__,__FILE__)
                      ENDIF
                   ENDDO
                   ! (GM) local thermostats
                ELSEIF ( keyword_contains(line,'LOCAL') ) THEN
                   loct%tloct=.TRUE.
                   tcafes=.FALSE.
                   nosl%tmnose=.FALSE.
                   nosl%tultra=.FALSE.
                   previous_line = line
                   READ(iunit,'(A)',iostat=ierr) line
                   first=1
                   ! (GM) NUMBER OF LOCAL THEROMOSTATS
                   CALL readsi(line,first,last,loct%nloct,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   ALLOCATE(loctpin(2,loct%nloct),STAT=ierr)
                   IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                        __LINE__,__FILE__)
                   ALLOCATE(loctt0(loct%nloct),STAT=ierr)
                   IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                        __LINE__,__FILE__)
                   ! (GM) PARAMETERS FOR EACH THERMOSTAT
                   DO i=1,loct%nloct
                      READ(iunit,'(A)',iostat=ierr) line
                      first=1
                      CALL readsr(line,first,last,loctpin(1,i),erread)
                      IF (erread) THEN
                         error_message        = "ERROR WHILE READING VALUE"
                         something_went_wrong = .true.
                         go_on_reading        = .false.
                      ENDIF
                      first=last
                      CALL readsr(line,first,last,loctpin(2,i),erread)
                      IF (erread) THEN
                         error_message        = "ERROR WHILE READING VALUE"
                         something_went_wrong = .true.
                         go_on_reading        = .false.
                      ENDIF
                      IF ( keyword_contains(previous_line,'T0') ) THEN
                         first=last
                         CALL readsr(line,first,last,loctt0(i),erread)
                         IF (erread) THEN
                            error_message        = "ERROR WHILE READING VALUE"
                            something_went_wrong = .true.
                            go_on_reading        = .false.
                         ENDIF
                      ELSE
                         loctt0(i)=loctpin(1,i)
                      ENDIF
                   ENDDO
                   cntr%tempw=loctpin(1,1)
                   ! (GM) NEXT LINE: NUMBER OF LOCAL THERMOSTAT RANGES
                   READ(iunit,'(A)',iostat=ierr) line
                   first=1
                   CALL readsi(line,first,last,loct%nlocrng,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   ALLOCATE(lctrng(3,loct%nlocrng),STAT=ierr)
                   IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                        __LINE__,__FILE__)
                   ! (GM) READ IN THE RANGES:  <THERMOSTAT NO>  <START ATOM>  <END ATOM>
                   DO i=1,loct%nlocrng
                      READ(iunit,'(A)',iostat=ierr) line
                      first=1
                      CALL readsi(line,first,last,lctrng(1,i),erread)
                      IF (erread) THEN
                         error_message        = "ERROR WHILE READING VALUE"
                         something_went_wrong = .true.
                         go_on_reading        = .false.
                      ENDIF
                      first=last
                      CALL readsi(line,first,last,lctrng(2,i),erread)
                      IF (erread) THEN
                         error_message        = "ERROR WHILE READING VALUE"
                         something_went_wrong = .true.
                         go_on_reading        = .false.
                      ENDIF
                      first=last
                      CALL readsi(line,first,last,lctrng(3,i),erread)
                      IF (erread) THEN
                         error_message        = "ERROR WHILE READING VALUE"
                         something_went_wrong = .true.
                         go_on_reading        = .false.
                      ENDIF
                   ENDDO
                ELSE
                   previous_line = line
                   READ(iunit,'(A)',iostat=ierr) line
                   first=1
                   CALL readsr(line,first,last,cntr%tempw, erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   first=last
                   CALL readsr(line,first,last,cntr%wnosp0,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   IF ( keyword_contains(previous_line,'T0') ) THEN
                      first=last
                      CALL readsr(line,first,last,cntr%nospt0,erread)
                   ELSE
                      cntr%nospt0=cntr%tempw
                   ENDIF
                ENDIF
             ELSEIF ( keyword_contains(line,'NOSE',and='ELECTRONS') ) THEN
                ! Nose thermostats for ions
                cntl%tnosee=.TRUE.
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                first=1
                CALL readsr(line,first,last,cntr%ekinw, erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                first=last
                CALL readsr(line,first,last,cntr%wnose0,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF ( keyword_contains(line,'NOSE',and='CELL') ) THEN
                ! Nose thermostats for cell
                cntl%tnosec=.TRUE.
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                first=1
                CALL readsr(line,first,last,cntr%tempc, erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                first=last
                CALL readsr(line,first,last,cntr%wnosc0,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF ( keyword_contains(line,'NOSE',and='PARAMETER') .OR. &
                    keyword_contains(line,'NOSE',and='PARAMETERS') ) THEN
                ! Nose thermostats explicit form
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                first=1
                CALL readsi(line,first,last,cnti%nchp,   erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                first=last
                CALL readsi(line,first,last,cnti%nchs,   erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                first=last
                CALL readsi(line,first,last,cnti%nchb,   erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                first=last
                CALL readsr(line,first,last,cntr%nedof0, erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                first=last
                CALL readsi(line,first,last,cnti%ncalls0,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                first=last
                CALL readsi(line,first,last,cnti%nit0,   erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF ( keyword_contains(line,'USE_IN_STREAM') ) THEN
                cntl%is_in_stream=.TRUE.
             ELSEIF ( keyword_contains(line,'USE_OUT_STREAM') ) THEN
                cntl%is_out_stream=.TRUE.
             ELSEIF ( keyword_contains(line,'USE_CUBLAS') ) THEN
                cp_cuda_env%use_blas = .TRUE.
             ELSEIF ( keyword_contains(line,'USE_CUFFT') ) THEN
                cp_cuda_env%use_fft = .TRUE.
             ELSEIF ( keyword_contains(line,'BLAS_N_STREAMS_PER_DEVICE') ) THEN
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                first=1
                CALL readsi(line ,first,last,cp_cuda_env%blas_n_streams_per_device,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF ( keyword_contains(line,'BLAS_N_DEVICES_PER_TASK') ) THEN
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                first=1
                CALL readsi(line ,first,last,cp_cuda_env%blas_n_devices_per_task,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF ( keyword_contains(line,'FFT_N_STREAMS_PER_DEVICE') ) THEN
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                first=1
                CALL readsi(line ,first,last,cp_cuda_env%fft_n_streams_per_device,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF ( keyword_contains(line,'FFT_N_DEVICES_PER_TASK') ) THEN
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                first=1
                CALL readsi(line ,first,last,cp_cuda_env%fft_n_devices_per_task,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF ( keyword_contains(line,'N_STREAMS') ) THEN
                WRITE(output_unit,*) 'The N_STREAMS keyword is deprecated'
                CALL stopgm(procedureN,'The N_STREAMS keyword is deprecated',& 
                     __LINE__,__FILE__)
             ELSEIF ( keyword_contains(line,'USE_MPI_IO') ) THEN
                cntl%use_mpi_io=.TRUE.
             ELSEIF ( keyword_contains(line,'TRACE') ) THEN
                ! Tracing
                cp_trace%ttrace = .FALSE.
                cp_trace%ttrace_master_only = .FALSE.
                IF ( keyword_contains(line,'ALL') ) THEN
                   cp_trace%ttrace = .TRUE.
                ELSEIF ( keyword_contains(line,'MASTER') ) THEN
                   cp_trace%ttrace = .TRUE.
                   cp_trace%ttrace_master_only = .TRUE.
                ENDIF
             ELSEIF ( keyword_contains(line,'TRACE_PROCEDURE') ) THEN
                ! Tracing selected procedure
                cp_trace%ttrace = .TRUE.
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                tname%trace_nbr_procedure=tname%trace_nbr_procedure+1
                tname%trace_procedure(tname%trace_nbr_procedure)=TRIM(ADJUSTL(line))
             ELSEIF ( keyword_contains(line,'TRACE_MAX_DEPTH') ) THEN
                ! Tracing select max depth
                cp_trace%ttrace = .TRUE.
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                first=1
                CALL readsi(line ,first,last,tname%trace_max_depth,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF ( keyword_contains(line,'TRACE_MAX_CALLS') ) THEN
                ! Tracing select max calls
                cp_trace%ttrace = .TRUE.
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                first=1
                CALL readsi(line ,first,last,tname%trace_max_calls,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF ( keyword_contains(line,'PRINT') ) THEN
                ! Print intermediate results
                cprint%tprint=.TRUE.
                IF ( keyword_contains(line,'OFF') ) THEN
                   iflag=-1
                ELSE
                   iflag=1
                ENDIF
                IF ( keyword_contains(line,'INFO'))  cprint%iprint(iprint_info)  = iflaG
                IF ( keyword_contains(line,'EIGENVALUES')) cprint%iprint(iprint_eigen) = iflaG
                IF ( keyword_contains(line,'COORDINATES'))  cprint%iprint(iprint_coor)  = iflaG
                IF ( keyword_contains(line,'FORCES'))  cprint%iprint(iprint_force) = iflaG
                IF ( keyword_contains(line,'WANNIER'))  cprint%iprint(iprint_wann)  = iflaG
                IF ( keyword_contains(line,'LSCAL')) cprint%iprint(iprint_lscal) = iflaG
                IF ( keyword_contains(line,'ALL') ) THEN
                   cprint%iprint(iprint_info)  = iflag
                   cprint%iprint(iprint_eigen) = iflag
                   cprint%iprint(iprint_coor)  = iflag
                   cprint%iprint(iprint_force) = iflag
                   cprint%iprint(iprint_wann)  = iflag
                ENDIF
                IF ( keyword_contains(line,'ENERGY') ) THEN
                   IF ( keyword_contains(line,'EKIN'))  cprint%iprint(iprint_ekin)  = iflaG
                   IF ( keyword_contains(line,'EFREE')) cprint%iprint(iprint_eeig)  = iflaG
                   IF ( keyword_contains(line,'EBAND')) cprint%iprint(iprint_eband) = iflaG
                   IF ( keyword_contains(line,'ENTROPY'))&
                        cprint%iprint(iprint_entropy) = iflag
                   IF ( keyword_contains(line,'ELECTROSTATIC'))  THEN
                      cprint%iprint(iprint_eht)   = iflag
                      IF ( keyword_contains(line,'ELEC1'))&
                           cprint%iprint(iprint_elec1)   = iflag
                      IF ( keyword_contains(line,'ELEC2'))&
                           cprint%iprint(iprint_elec2)   = iflag
                   ENDIF
                   IF ( keyword_contains(line,'EHEP'))  cprint%iprint(iprint_ehep)  = iflaG
                   IF ( keyword_contains(line,'EHEE'))  cprint%iprint(iprint_ehee)  = iflaG
                   IF ( keyword_contains(line,'EHII'))  cprint%iprint(iprint_ehii)  = iflaG
                   IF ( keyword_contains(line,'ESR'))   cprint%iprint(iprint_esr)   = iflag
                   IF ( keyword_contains(line,'ESELF')) cprint%iprint(iprint_eself) = iflaG
                   IF ( keyword_contains(line,'EPSEU')) cprint%iprint(iprint_epseu) = iflaG
                   IF ( keyword_contains(line,'ENL'))   cprint%iprint(iprint_enl)   = iflag
                   IF ( keyword_contains(line,'EXCH'))  cprint%iprint(iprint_exc)   = iflaG
                   IF ( keyword_contains(line,'VXC'))   cprint%iprint(iprint_vxc)   = iflag
                   IF ( keyword_contains(line,'EGC'))   cprint%iprint(iprint_egc)   = iflag
                   IF ( keyword_contains(line,'EBOGO')) cprint%iprint(iprint_ebogo) = iflaG
                   IF ( keyword_contains(line,'ETOT1')) cprint%iprint(iprint_etot1) = iflaG
                   IF ( keyword_contains(line,'ETOT2')) cprint%iprint(iprint_etot2) = iflaG
                   IF ( keyword_contains(line,'ECAS'))  cprint%iprint(iprint_ecas)  = iflaG
                   IF ( keyword_contains(line,'EPEN'))  cprint%iprint(iprint_epen)  = iflaG
                   IF ( keyword_contains(line,'ETDDFT'))cprint%iprint(iprint_etddft)= iflaG
                   IF ( keyword_contains(line,'EHUB'))  cprint%iprint(iprint_ehub)  = iflaG
                ENDIF
                IF (cntl%tsic) cprint%iprint(iprint_ehsic) = 1! cmb_ssic
                IF ( .not. keyword_contains(line,'ON') .AND. .not. keyword_contains(line,'OFF') ) THEN
                   previous_line = line
                   READ(iunit,'(A)',iostat=ierr) line
                   CALL readsi(line,1,last,cprint%iprint_step,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ENDIF
             ELSEIF ( keyword_contains(line,'STRUCTURE') ) THEN
                ! Print structure parameters
                IF ( keyword_contains(line,'BONDS')) bond=.TRUE.
                IF ( keyword_contains(line,'ANGLES')) angle=.TRUE.
                IF ( keyword_contains(line,'DIHEDRALS')) dihedral=.TRUE.
                IF ( keyword_contains(line,'SELECT') ) THEN
                   cntl%tssel=.TRUE.
                   READ(iunit,*,iostat=ierr) nssel
                   IF (ierr /= 0) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   ALLOCATE(nassel(nssel),STAT=ierr)
                   IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                        __LINE__,__FILE__)
                   READ(iunit,*,iostat=ierr) (nassel(i),i=1,nssel)
                   IF (ierr /= 0) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ELSE
                   cntl%tssel=.FALSE.
                   nssel=0
                ENDIF
             ELSEIF ( keyword_contains(line,'STORE') ) THEN
                ! Store intermediate results
                previous_line = line
                IF ( keyword_contains(line,'OFF') ) THEN
                   store_or_not=.FALSE.
                ELSE
                   store_or_not=.TRUE.
                ENDIF
                IF ( keyword_contains(line,'POTENTIAL')) store1%spot=.TRUE..AND.store_or_not
                IF ( keyword_contains(line,'WAVEFUNCUNCTIONS',alias='WAVEFUNCTION')) store1%swf=.TRUE..AND.store_or_not
                IF ( keyword_contains(line,'DENSITY') ) THEN
                   tsrho       = .TRUE..AND.store_or_not
                   store1%srho = .TRUE..AND.store_or_not
                ENDIF
                IF (store_or_not) THEN
                   READ(iunit,'(A)',iostat=ierr) line
                   CALL readsi(line,1,last,store1%istore,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   IF ( keyword_contains(line,'SC',cut_at='=') ) THEN
                      first = index_of_delimiter(line,'SC','=')
                      CALL readsi(line,first,last,store1%isctore,erread)
                      IF (erread) THEN
                         error_message        = "ERROR WHILE READING VALUE"
                         something_went_wrong = .true.
                         go_on_reading        = .false.
                      ENDIF
                   ENDIF
                ENDIF
             ELSEIF ( keyword_contains(line,'RESTFILE',alias='RESTFILES') ) THEN
                ! Number of different restart files
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                CALL readsi(line,1,last,cnti%nrestf,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                IF ( keyword_contains(previous_line,'SAMPLE') ) THEN
                   IF (cnti%nrestf > maxrf)&
                        CALL stopgm('CONTROL','TOO MANY RESTART FILES',& 
                        __LINE__,__FILE__)
                   ! We can use more than one line.
                   READ(iunit,*,iostat=ierr) (restf%nstepwr(i),i=1,cnti%nrestf)
                   IF (ierr /= 0) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ENDIF
             ELSEIF ( keyword_contains(line,'CONVERGENCE') ) THEN
                ! Convergence criteria
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                IF ( keyword_contains(previous_line,'ORBITAL',alias='ORBITALS') .OR. &
                     keyword_contains(previous_line,'WAVEFUNCTION') ) THEN
                   CALL readsr(line,1,last,cntr%tolog,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   first = last
                   CALL readsr(line,first,last,cntr%toldetot,erread)
                   !
                   ! In case no toldetot is specified, nothing happens, so we
                   ! set erread to false.
                   !
                   erread=.FALSE.
                ELSEIF ( keyword_contains(previous_line,'GEOMETRY') ) THEN
                   CALL readsr(line,1,last,cntr%tolng,erread)
                ELSEIF ( keyword_contains(previous_line,'CELL') ) THEN
                   CALL readsr(line,1,last,cntr%tolcg,erread)
                ELSEIF ( keyword_contains(previous_line,'ADAPT') ) THEN
                   CALL readsr(line,1,last,cntr%tolad,erread)
                ELSEIF ( keyword_contains(previous_line,'ENERGY') ) THEN
                   CALL readsr(line,1,last,cntr%tolene,erread)
                ELSEIF ( keyword_contains(previous_line,'CALFOR') ) THEN
                   CALL readsr(line,1,last,cntr%tolfor,erread)
                ELSEIF ( keyword_contains(previous_line,'RELAX') ) THEN
                   CALL readsi(line,1,last,cnti%nstcnv,erread)
                ELSEIF ( keyword_contains(previous_line,'RHOFIX') ) THEN
                   CALL readsr(line,1,last,cntr%tolrhofix,erread)
                ELSEIF ( keyword_contains(previous_line,'INITIAL') ) THEN
                   CALL readsr(line,1,last,cntr%tolini,erread)
                ELSEIF ( keyword_contains(previous_line,'CONSTRAINT') ) THEN
                   CALL readsr(line,1,last,cdftcom%vccon,erread)
                   first=last
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   CALL readsr(line,first,last,cdftcom%vcconu,erread)
                   IF (erread) cdftcom%vcconu=cdftcom%vccon
                   erread=.FALSE.
                ELSE
                   first=1
                   CALL readsr(line,first,last,cntr%tolog,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   first=last
                   CALL readsr(line,first,last,cntr%tolng,erread)
                   !
                   ! In case the read goes wrong, we want a crash -
                   ! otherwise, tolng will be set to 0.0 !
                   !
                   ! erread=.FALSE.
                ENDIF
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF ( keyword_contains(line,'ANNEAL',alias='ANNEALING') ) THEN
                ! Simulated Annealing
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                IF ( keyword_contains(previous_line,'IONS') ) THEN
                    cntl%annei=.TRUE.
                    IF ( keyword_contains(previous_line,'DUAL') ) THEN
                       cntl%anneal_dual = .TRUE.
                       cntr%anneal_factors = 1.0_real_8
                       first = 1
                       DO i = 1, 2
                          CALL readsr(line,first,last,cntr%anneal_factors(i),erread)
                          first = last
                          IF (erread) THEN
                             error_message        = "ERROR WHILE READING VALUE"
                             something_went_wrong = .true.
                             go_on_reading        = .false.
                          ENDIF
                       END DO
                   ELSE
                      CALL readsr(line,1,last,cntr%anneri,erread)
                   END IF
                ELSEIF ( keyword_contains(previous_line,'ELECTRONS',alias='ELECTRON') ) THEN
                   cntl%annee=.TRUE.
                   CALL readsr(line,1,last,cntr%annere,erread)
                ELSEIF ( keyword_contains(previous_line,'CELL') ) THEN
                   cntl%annec=.TRUE.
                   CALL readsr(line,1,last,cntr%annerc,erread)
                ELSE
                   cntl%annei=.TRUE.
                   cntl%annee=.TRUE.
                   cntl%annec=.TRUE.
                   first=1
                   CALL readsr(line,first,last,cntr%anneri,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   first=last
                   CALL readsr(line,first,last,cntr%annere,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   first=last
                   CALL readsr(line,first,last,cntr%annerc,erread)
                ENDIF
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF ( keyword_contains(line,'DAMPING',alias='DAMP') ) THEN
                ! Damped Dynamics
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                IF ( keyword_contains(previous_line,'IONS') ) THEN
                   cntl%dampi=.TRUE.
                   CALL readsr(line,1,last,cntr%dampgi,erread)
                ELSEIF ( keyword_contains(previous_line,'ELECTRONS') ) THEN
                   cntl%dampe=.TRUE.
                   CALL readsr(line,1,last,cntr%dampge,erread)
                ELSEIF ( keyword_contains(previous_line,'CELL') ) THEN
                   cntl%dampc=.TRUE.
                   CALL readsr(line,1,last,cntr%dampgc,erread)
                ELSE
                   cntl%dampi=.TRUE.
                   cntl%dampe=.TRUE.
                   cntl%dampc=.TRUE.
                   first=1
                   CALL readsr(line,first,last,cntr%dampgi,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   first=last
                   CALL readsr(line,first,last,cntr%dampge,erread)
                   IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                   first=last
                   CALL readsr(line,first,last,cntr%dampgc,erread)
                ENDIF
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF ( keyword_contains(line,'MOLECULE',and='ISOLATED') ) THEN
                ! Isolated molecule
                isos1%tisos=.TRUE.
             ELSEIF ( keyword_contains(line,'MOLECULE',and='CENTER') ) THEN
                ! Put COM to center of box
                IF ( keyword_contains(line,'OFF') ) THEN
                   isos1%tcent=.FALSE.
                ELSE
                   isos1%tcent=.TRUE.
                ENDIF
             ELSEIF ( keyword_contains(line,'EXTERNAL',and='POTENTIAL') ) THEN
                cntl%texpot=.TRUE.
                IF ( keyword_contains(line,'ADD')) cntl%texadd=.TRUE.
             ELSEIF ( keyword_contains(line,'ELECTROSTATIC',and='POTENTIAL') ) THEN
                cntl%tepot=.TRUE.
                rout1%rhoout=.TRUE.
                IF ( keyword_contains(line,'SAMPLE',cut_at='=') ) THEN
                   first = index_of_delimiter(line,'SAMPLE','=')
                   CALL readsi(line,first,last,rout1%nrhoout,erread)
                   IF (erread .or. rout1%nrhoout<0) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ENDIF
             ELSEIF ( keyword_contains(line,'DIPOLE',and='DYNAMICS') ) THEN
                ! Store Dipole moments
                cntl%tdipd=.TRUE.
                IF ( keyword_contains(line,'WANNIER')) wannl%twann=.TRUE.
                IF ( keyword_contains(line,'SAMPLE') ) THEN
                   previous_line = line
                   READ(iunit,'(A)',iostat=ierr) line
                   CALL readsi(line,1,last,cnti%npdip,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ENDIF
             ELSEIF ( keyword_contains(line,'WANNIER',and='PARAMETER') .OR. &
                      keyword_contains(line,'WANNIER',and='PARAMETERS') ) THEN
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                first=1
                CALL readsr(line,first,last,wannr%w_step,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                first=last
                CALL readsr(line,first,last,wannr%w_eps, erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                first=last
                CALL readsr(line,first,last,wannr%w_ran, erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                first=last
                CALL readsi(line,first,last,wanni%w_maxs,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF ( keyword_contains(line,'WANNIER',and='OPTIMIZATION') ) THEN
                IF ( keyword_contains(line,'SD')) wanni%w_opt=1  
                IF ( keyword_contains(line,'JACOBI')) wanni%w_opt=2  
                IF ( keyword_contains(line,'SVD')) wanni%w_opt=3
             ELSEIF ( keyword_contains(line,'WANNIER',and='NPROC') ) THEN
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                first=1
                CALL readsi(line,first,last,wan05%loc_npgrp,erread)
                IF(erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                IF(wan05%loc_npgrp>parai%nproc.OR.wan05%loc_npgrp<=0) THEN
                   WRITE(output_unit,*)'WARNING: LOC_NPGRP must be at most NPROC, ' &
                        // 'reset LOC_NPGRP = NPROC'
                   wan05%loc_npgrp=parai%nproc
                ENDIF
             ELSEIF ( keyword_contains(line,'WANNIER',and='RELOCALIZE_IN_SCF') ) THEN
                wan05%loc_relocalize_in_scf=.TRUE.
             ELSEIF ( keyword_contains(line,'WANNIER',and='RECOMPUTE_DIPOLE_MATRICES_EVERY') ) THEN
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                first=1
                CALL readsi(line,first,last,wan05%loc_recompute_dipole_matrices_every,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF ( keyword_contains(line,'WANNIER',and='RELOCALIZE_EVERY') ) THEN
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                first=1
                wan05%loc_relocalize=.TRUE.
                CALL readsi(line,first,last,wan05%loc_relocalize_every,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF ( keyword_contains(line,'WANNIER',and='TYPE')) THEN
                IF ( keyword_contains(line,'VANDERBILT')) wanni%w_type=1  
                IF ( keyword_contains(line,'RESTA')) wanni%w_type=2  
             ELSEIF ( keyword_contains(line,'WANNIER',and='REFERENCE') ) THEN
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                first=1
                DO i=1,3
                   CALL readsr(line,first,last,wannr%w_ref(i),erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   first=last
                ENDDO
             ELSEIF ( keyword_contains(line,'WANNIER',and='SERIAL')) THEN
                cntl%twserial=.TRUE.
             ELSEIF ( keyword_contains(line,'WANNIER',and='DOS')) THEN
                wannl%twdos=.TRUE.
             ELSEIF ( keyword_contains(line,'WANNIER',and='MOLECULAR)')) THEN
                wannl%twmol=.TRUE.
             ELSEIF ( keyword_contains(line,'WANNIER',and='WFNOUT')) THEN
                wannl%twpri=.TRUE.
                IF ( keyword_contains(line,'DENSITY')) wannl%tsdens=.TRUE.
                IF ( keyword_contains(line,'ALL') ) THEN
                   wanni%sw_all=1
                ELSE
                   previous_line = line
                   READ(iunit,'(A)',iostat=ierr) line
                   IF ( keyword_contains(previous_line,'PARTIAL') ) THEN
                      first=1
                      CALL readsi(line,first,last,wanni%sw_first,erread)
                      IF (erread) THEN
                         error_message        = "ERROR WHILE READING VALUE"
                         something_went_wrong = .true.
                         go_on_reading        = .false.
                      ENDIF
                      first=last
                      CALL readsi(line,first,last,wanni%sw_last,erread)
                      IF (erread) THEN
                         error_message        = "ERROR WHILE READING VALUE"
                         something_went_wrong = .true.
                         go_on_reading        = .false.
                      ENDIF
                   ELSEIF ( keyword_contains(previous_line,'LIST') ) THEN
                      CALL readsi(line,1,last,wanni%sw_orb,erread)
                      IF (erread) THEN
                         error_message        = "ERROR WHILE READING VALUE"
                         something_went_wrong = .true.
                         go_on_reading        = .false.
                      ENDIF
                      ALLOCATE(sw_list(wanni%sw_orb),STAT=ierr)
                      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                           __LINE__,__FILE__)
                      ! We can use more than one line
                      READ(iunit,*,iostat=ierr) (sw_list(i),i=1,wanni%sw_orb)
                   ELSEIF ( keyword_contains(previous_line,'SPREAD') ) THEN
                      CALL readsr(line,1,last,wannr%sw_spread,erread)
                      IF (erread) THEN
                         error_message        = "ERROR WHILE READING VALUE"
                         something_went_wrong = .true.
                         go_on_reading        = .false.
                      ENDIF
                   ELSE
                      error_message        = "ERROR IN WANNIER PARAMETERS"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ENDIF
             ELSEIF ( keyword_contains(line,'PARA_USE_MPI_IN_PLACE') ) THEN
                para_use_mpi_in_place=.TRUE.
             ELSEIF ( keyword_contains(line,'PARA_BUFF_SIZE') ) THEN
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                first=1
                CALL readsi(line,first,last,para_buff_size,erread)
                IF(erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF ( keyword_contains(line,'PARA_STACK_BUFF_SIZE') ) THEN
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                first=1
                CALL readsi(line,first,last,para_buff_size,erread)
                IF(erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF ( keyword_contains(line,'RHOOUT') ) THEN
                ! Store density
                rout1%rhoout=.TRUE.
                IF ( keyword_contains(line,'SAMPLE') ) THEN
                   first = index_of_delimiter(line,'SAMPLE','=')
                   CALL readsi(line,first,last,rout1%nrhoout,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   IF (rout1%nrhoout<0) THEN
                      error_message        = "WRONG VALUE FOR RHOOUT SAMPLE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ENDIF
                IF ( keyword_contains(line,'BANDS') ) THEN
                   previous_line = line
                   READ(iunit,'(A)',iostat=ierr) line
                   CALL readsi(line,1,last,numpr,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   IF (numpr<=0) THEN
                      error_message        = "WRONG VALUE FOR RHOOUT BANDS"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   ALLOCATE(mwfn(numpr),STAT=ierr)
                   IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                        __LINE__,__FILE__)
                   ! We can use more than one line
                   READ(iunit,*,iostat=ierr) (mwfn(i),i=1,numpr )
                   IF (ierr /= 0) THEN
                      error_message        = "ERROR WHILE READING RHOOUT BANDS"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   rout1%rhoout=.TRUE.
                ENDIF
             ELSEIF ( keyword_contains(line,'ELF') ) THEN
                elfcb%telf=.TRUE.
                rout1%rhoout=.TRUE.
                IF ( keyword_contains(line,'PARAMETER') ) THEN
                   previous_line = line
                   READ(iunit,'(A)',iostat=ierr) line
                   first=1
                   CALL readsr(line,first,last,elfcb%elfcut,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   first=last
                   CALL readsr(line,first,last,elfcb%elfeps,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ENDIF
             ELSEIF ( keyword_contains(line,'ENERGYBAND',alias='ENERGYBANDS') ) THEN
                rout1%teband=.TRUE.
             ELSEIF ( keyword_contains(line,'HESSIAN') ) THEN
                ! Initial hessian
                IF ( keyword_contains(line,'UNIT')) cnti%npara=-1
                IF ( keyword_contains(line,'DISCO')) cnti%npara=0
                IF ( keyword_contains(line,'SCHLEGEL')) cnti%npara=1
                IF ( keyword_contains(line,'PARTIAL')) cnti%npara=cnti%npara+3
             ELSEIF ( keyword_contains(line,'PROJECT') ) THEN
                ! Electronic gradient projection to be used
                IF ( keyword_contains(line,'NONE')) cnti%iproj=0
                IF ( keyword_contains(line,'DIAGONAL')) cnti%iproj=1
                IF ( keyword_contains(line,'FULL')) cnti%iproj=2
             ELSEIF ( keyword_contains(line,'TRAJECTORY',alias='TRAJ') ) THEN
                ! Store trajectory
                rout1%rout=.TRUE.
                IF ( keyword_contains(line,'XYZ')) rout1%xtout=.TRUE.
                IF ( keyword_contains(line,'DCD')) rout1%dcout=.TRUE.
                IF ( keyword_contains(line,'OFF') ) THEN
                   rout1%rout=.FALSE.
                   rout1%xtout=.FALSE.
                   rout1%dcout=.FALSE.
                ENDIF
                IF ( keyword_contains(line,'BINARY',alias='BIN') ) THEN
                   cprint%twritebintrajectory=.TRUE.
                ENDIF
                IF ( keyword_contains(line,'FORCES') ) THEN
                   cprint%twriteforcetrajectory=.TRUE.
                ENDIF
                IF ( keyword_contains(line,'SMALL')) trajsmall=.TRUE.
                IF ( keyword_contains(line,'RANGE') ) THEN
                   previous_line = line
                   READ(iunit,*,iostat=ierr)&
                        cprint%minwriteatom, cprint%maxwriteatom
                ENDIF
                IF ( keyword_contains(line,'SAMPLE') ) THEN
                   previous_line = line
                   READ(iunit,'(A)',iostat=ierr) line
                   CALL readsi(line,1,last,cnti%ntraj,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   ! turn off trajectory writing for NTRAJ=0
                   ! and disable TRAJECTORY file, but not TRAJEC.xyz for NTRAJ < 0
                   IF (cnti%ntraj==0) THEN
                      cnti%ntraj=1
                      rout1%rout=.FALSE.
                      rout1%xtout=.FALSE.
                      rout1%dcout=.FALSE.
                   ELSEIF (cnti%ntraj<0) THEN
                      rout1%rout=.FALSE.
                      cnti%ntraj=-cnti%ntraj
                   ENDIF
                   IF (trajsmall) THEN
                      READ(iunit,'(A)',iostat=ierr) line
                      CALL readsi(line,1,last,trajsmalln,erread)
                   ENDIF
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ENDIF
             ELSEIF ( keyword_contains(line,'MOVIE') ) THEN
                ! Write movie file every imovie steps    
                rout1%mout=.TRUE.
                IF ( keyword_contains(line,'OFF')) rout1%mout=.FALSE.
                IF ( keyword_contains(line,'SAMPLE') ) THEN
                   previous_line = line
                   READ(iunit,'(A)',iostat=ierr) line
                   CALL readsi(line,1,last,cnti%imovie,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ENDIF
             ELSEIF ( keyword_contains(line,'STRESS',and='TENSOR') ) THEN
                ! Stress tensor calculation
                cntl%tpres=.TRUE.
                IF ( keyword_contains(line,'VIRIAL') ) cntl%tvirial=.TRUE.
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                CALL readsi(line,1,last,cnti%npres,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF ( keyword_contains(line,'CLASSTRESS') ) THEN
                cntl%tprec=.TRUE.
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                CALL readsi(line,1,last,cnti%nprec,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF ( keyword_contains(line,'COMPRESS') ) THEN
                ! Wavefunction in compressed form          
                cntl%wcomp=.TRUE.
                cnti%wcompb=2
                IF ( keyword_contains(line,'WRITEAO64')) cnti%wcompb=-1
                IF ( keyword_contains(line,'WRITEAO32')) cnti%wcompb=-2
                IF ( keyword_contains(line,'WRITEAO16')) cnti%wcompb=-4
                IF ( keyword_contains(line,'WRITEAO8' )) cnti%wcompb=-8
                IF ( keyword_contains(line,'WRITE64')) cnti%wcompb=1
                IF ( keyword_contains(line,'WRITE32')) cnti%wcompb=2
                IF ( keyword_contains(line,'WRITE16')) cnti%wcompb=4
                IF ( keyword_contains(line,'WRITE8' )) cnti%wcompb=8
             ELSEIF ( keyword_contains(line,'MEMORY') ) THEN
                ! Use small memory version of Vanderbilt code
                IF ( keyword_contains(line,'SMALL')) cntl%bigmem=.FALSE.
                IF ( keyword_contains(line,'BIG')) cntl%bigmem=.TRUE.
                ! Check memory frequently
                IF ( keyword_contains(line,'CHECK')) cntl%tmemchk=.TRUE.
             ELSEIF ( keyword_contains(line,'REALSPACE',and='WFN') .OR. &
                      keyword_contains(line,'SPACE',and='WFN') ) THEN
                ! Info on real space wavefunctions
                IF ( keyword_contains(line,'KEEP') ) THEN
                   cntl%krwfn=.TRUE.
                ELSE
                   cntl%krwfn=.FALSE.
                ENDIF
                IF ( keyword_contains(line,'SIZE') ) THEN
                   READ(iunit,'(A)',iostat=ierr) line
                   CALL readsr(line,1,last,cntr%memsize,erread)
                ELSE
                   cntr%memsize=-1._real_8
                ENDIF
             ELSEIF ( keyword_contains(line,'SPLINE') ) THEN
                IF ( keyword_contains(line,'POINTS') ) THEN
                   ! Number of spline points
                   previous_line = line
                   READ(iunit,'(A)',iostat=ierr) line
                   CALL readsi(line,1,last,cnti%nsplp,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ENDIF
                IF ( keyword_contains(line,'QFUNCTION')) qspl1%qspline=.TRUE.
                IF ( keyword_contains(line,'INIT')) qspl1%initspl=.TRUE.
                IF ( keyword_contains(line,'RANGE') ) THEN
                   previous_line = line
                   READ(iunit,'(A)',iostat=ierr) line
                   CALL readsr(line,1,last,qsrang,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ENDIF
             ELSEIF ( keyword_contains(line,'FINITE',and='DIFFERENCE') .OR. &
                      keyword_contains(line,'FINITE',and='DIFFERENCES')) THEN
                ! Step length for finite differences
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                CALL readsr(line,1,last,cntr%fdiff,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                IF ( keyword_contains(line,'COORD',cut_at='=') ) THEN
                   first = index_of_delimiter(line,'COORD','=')
                   tref_fdiff=.TRUE.
                   ! Position of the finite differences sphere.
                   DO i=1,3
                      CALL readsr(line,first,last,coord_fdiff(i),erread)
                      IF (erread) THEN
                         error_message        = "ERROR WHILE READING VALUE"
                         something_went_wrong = .true.
                         go_on_reading        = .false.
                      ENDIF
                      first=last+1
                   ENDDO
                   IF ( keyword_contains(line,'RADIUS',cut_at='=') ) THEN
                      first = index_of_delimiter(line,'RADIUS','=')
                      CALL readsr(line,first,last,r_fdiff,erread)
                      IF (erread) THEN
                         error_message        = "ERROR WHILE READING VALUE"
                         something_went_wrong = .true.
                         go_on_reading        = .false.
                      ENDIF
                   ELSE
                      error_message        = 'FINITE DIFFERENCE WITH ATOM NEEDS RADIUS'
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ENDIF
             ELSEIF ( keyword_contains(line,'TASKGROUP',alias='TASKGROUPS') ) THEN
                IF ( keyword_contains(line,'MAXI') ) THEN
                   CALL stopgm('CONTROL','TASKGROUPS/MAXI has been removed',& 
                        __LINE__,__FILE__)
                ELSEIF ( keyword_contains(line,'MINI') ) THEN
                   CALL stopgm('CONTROL','TASKGROUPS/MINI has been removed',& 
                        __LINE__,__FILE__)
                ELSE
                   IF ( keyword_contains(line,'CARTESIAN') ) THEN
                      cntl%tcart=.TRUE.
                   ENDIF
                   previous_line = line
                   READ(iunit,'(A)',iostat=ierr) line
                   CALL readsi(line,1,last,group%nogrp,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ENDIF
             ELSEIF ( keyword_contains(line,'CP_GROUPS') ) THEN
                ! fake read; CP_GROUPS are initialised at the very beginning
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
             ELSEIF ( keyword_contains(line,'DISTRIBUTE',and='FNL') .OR. &
                      keyword_contains(line,'DISTRIBUTED',and='FNL') ) THEN
                ! storage form of fnl
                IF ( keyword_contains(line,'ON')) cntl%tfdist=.TRUE.
                IF ( keyword_contains(line,'OFF')) cntl%tfdist=.FALSE.
             ELSEIF ( keyword_contains(line,'FILEPATH') ) THEN
                ! Path to the restart files (all the line)
                READ(iunit,'(A)',iostat=ierr) fo_info%fpath
             ELSEIF ( keyword_contains(line,'FILE',and='FUSION')) THEN
                ! File fusion
                cntl%tfusi=.TRUE.
                READ(iunit,'(A)',iostat=ierr) s0_filn
                READ(iunit,'(A)',iostat=ierr) s1_filn
             ELSEIF ( keyword_contains(line,'FILE',and='MERGE')) THEN
                ! File merge
                cntl%tmerge=.TRUE.
                READ(iunit,'(A)',iostat=ierr) merge01%mfiln1
                READ(iunit,'(A)',iostat=ierr) merge01%mfiln2
                READ(iunit,'(A)',iostat=ierr) line
                first=1
                CALL readsi(line,first,last,merge02%mnst1,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                first=last
                CALL readsi(line,first,last,merge02%mnst2,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                first=last
                CALL readsi(line,first,last,merge02%mortho,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                READ(iunit,'(A)',iostat=ierr) line
                first=1
                CALL readsr(line,first,last,merge01%mshift(1),erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                first=last
                CALL readsr(line,first,last,merge01%mshift(2),erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                first=last
                CALL readsr(line,first,last,merge01%mshift(3),erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF ( keyword_contains(line,'FILE',and='SEPARATION')) THEN
                ! File separation
                cntl%tsep=.TRUE.
                READ(iunit,'(A)',iostat=ierr) s0_filn
                READ(iunit,'(A)',iostat=ierr) s1_filn
             ELSEIF ( keyword_contains(line,'SHIFT',and='POTENTIAL')) THEN
                ! Shift potential
                READ(iunit,'(A)',iostat=ierr) cntr%dshift
                ! Reverse Velocities (allow old 'INVERSE' and new 'REVERSE')
             ELSEIF ( keyword_contains(line,'REVERSE',and='VELOCITIES')) THEN
                cntl%trevers=.TRUE.
             ELSEIF ( keyword_contains(line,'BENCHMARK') ) THEN
                ! Info on Benchmark
                ! We can use more than one line.
                READ(iunit,*,iostat=ierr) (ibench(i),i=1,nbentr)
             ELSEIF ( keyword_contains(line,'MIRROR') ) THEN
                ! MIRROR input file on output file
                mirror = .TRUE.
             !
             ! TODO: Already covered in ODIIS NO_RESET, remove NO_RESET as a standalone keyword
             !
             ELSEIF ( keyword_contains(line,'NO_RESET') ) THEN
                ! Number of steps between ODIIS resets if poor convergence
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                CALL readsi(line,1,last,cnti%nreset,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                ! LOCALIZAZION of WF in RECIPROCAL SPACE
             ELSEIF ( keyword_contains(line,'GLOCALIZATION',and='PARAMETER') .OR. &
                      keyword_contains(line,'GLOCALIZATION',and='PARAMETERS')) THEN
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                first=1
                CALL readsr(line,first,last,glocr%gloc_step,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                first=last
                CALL readsr(line,first,last,glocr%gloc_eps, erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                first=last
                CALL readsr(line,first,last,glocr%gloc_ran, erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                first=last
                CALL readsi(line,first,last,gloci%gloc_maxs,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                first=last
                CALL readsi(line,first,last,gloci%gloc_init,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF ( keyword_contains(line,'GLOCALIZATION',and='OPTIMIZATION')) THEN
                glocal%tgloc = .TRUE.
                IF ( keyword_contains(line,'SD') .AND. keyword_contains(line,'INTEGRATION'))&
                     gloci%gloc_opt=1
                IF ( keyword_contains(line,'ALMO') .AND. keyword_contains(line,'IDENTITY'))&
                     gloci%gloc_opt=2
             ELSEIF ( keyword_contains(line,'GFUNCTIONAL',and='TYPE')) THEN
                IF ( keyword_contains(line,'GSPREAD')) gloci%gloc_type=1  
                IF ( keyword_contains(line,'ZETANUMBER')) gloci%gloc_type=2
             ELSEIF ( keyword_contains(line,'SPREAD',and='RSPACE',cut_at='=')) THEN
                glocal%tgwannier = .TRUE.
                first = index_of_delimiter(line,'RSPACE','=')
                IF (first > 1) THEN
                   CALL readsr(line,first,last,glocr%g2g_weight,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   first=last
                   CALL readsr(line,first,last,glocr%wan_weight,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   ! IF(WAN_WEIGHT+G2G_WEIGHT > 1.0_real_8) 
                   ! *         CALL STOPGM('CONTROL','GLOC: SUM OF WEIGHTS > 1')
                ENDIF
                ! G2G_WEIGHT = 1._real_8
                ! WAN_WEIGHT = 0._real_8
             !
             ELSEIF ( keyword_contains(line,'PIPPO',cut_at='=') ) THEN
                ! ELSEIF (INDEX(LINE,'UNITARY') .AND.
                ! *           INDEX(LINE,'CONST') ) THEN
                IF ( keyword_contains(line,'LAGRANGE',cut_at='=') ) THEN
                   gloci%gloc_const = 1
                   first = index_of_delimiter(line,'LAGRANGE','=')
                   IF (first > 1) THEN
                      CALL readsr(line,first,last,glocr%gepslag,erread)
                      IF (erread) THEN
                         error_message        = "ERROR WHILE READING VALUE"
                         something_went_wrong = .true.
                         go_on_reading        = .false.
                      ENDIF
                      first=last
                      CALL readsi(line,first,last,gloci%gloc_maxit,erread)
                      IF (erread) THEN
                         error_message        = "ERROR WHILE READING VALUE"
                         something_went_wrong = .true.
                         go_on_reading        = .false.
                      ENDIF
                   ENDIF
                ELSEIF ( keyword_contains(line,'ORTHO') ) THEN
                   gloci%gloc_const = 2
                ELSEIF ( keyword_contains(line,'APPROX') ) THEN
                   gloci%gloc_const = 3
                ENDIF
             !
             ELSEIF ( keyword_contains(line,'STEP',and='TUNING')) THEN
                glocal%tg_linesearch = .TRUE.
             ELSEIF ( keyword_contains(line,'G_ANTISYM') ) THEN
                glocal%tg_antisymm = .TRUE.
                IF ( keyword_contains(line,'PENALTY')) glocal%tg_penalty = .TRUE.
             ELSEIF ( keyword_contains(line,'G_KICK') ) THEN
                glocal%tg_kick = .TRUE.
             ELSEIF ( keyword_contains(line,'G_COMPLEX') ) THEN
                glocal%tg_complex = .TRUE.
             ELSEIF ( keyword_contains(line,'G_REAL') ) THEN
                glocal%tg_real = .TRUE.
             ELSEIF ( keyword_contains(line,'READ',and='MATRIX') ) THEN
                glocal%tg_read_matrix = .TRUE.
             ELSEIF ( keyword_contains(line,'G_STEP') .AND.&
                  keyword_contains(line,'TUNE') ) THEN
                glocal%tg_tune_st = .TRUE.
             ELSEIF ( keyword_contains(line,'GLOC',and='WFNOUT') ) THEN
                glocal%tglocprint = .TRUE.
                IF ( keyword_contains(line,'REAL')) glocal%tglocrealp = .TRUE.
                IF ( keyword_contains(line,'ALL')) gloci%gloc_all=1
                previous_line = line
                IF ( keyword_contains(line,'PART',alias='PARTIAL') ) THEN
                   READ(iunit,'(A)',iostat=ierr) line
                   first=1
                   CALL readsi(line,first,last,gloci%gloc_first,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   first=last
                   CALL readsi(line,first,last,gloci%gloc_last,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ENDIF
                IF ( keyword_contains(previous_line,'LIST') ) THEN
                   READ(iunit,'(A)',iostat=ierr) line
                   CALL readsi(line,1,last,gloci%gloc_orb,erread)
                   IF (erread) THEN
                      error_message        = "ERROR WHILE READING VALUE"
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   ALLOCATE(gloc_list(gloci%gloc_orb),STAT=ierr)
                   IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                        __LINE__,__FILE__)
                   ! We can use more than one line
                   READ(iunit,*,iostat=ierr) (gloc_list(i),i=1,gloci%gloc_orb)
                ENDIF
             ELSEIF ( keyword_contains(line,'NO_GEO_CHECK',alias='NOGEOCHECK') ) THEN
                ! Dont check geometry for close atoms
                cntl%tnogeocheck=.TRUE.
             ELSEIF ( keyword_contains(line,'BROKEN') ) THEN
                ! Use Noodleman (broken symmetry) method for spin coupled systems
                cntl%bsymm=.TRUE.
                ! Distributed Linear Algebra
             ELSEIF ( keyword_contains(line,'DISTRIBUTED',and='LINALG',but_not='NEWORTHO') .OR. &
                      keyword_contains(line,'DIST',and='LINALG',but_not='NEWORTHO') ) THEN
                IF ( keyword_contains(line,'ON') ) THEN
                   cntl%tdmal=.TRUE.
                ELSEIF ( keyword_contains(line,'OFF') ) THEN
                   cntl%tdmal=.FALSE.
                ENDIF
                ! The new otrhogonalization of the wavefunction C.Bekas and A.Curioni Comp. Phys. Comm. (2010)
             ELSEIF ( keyword_contains(line,'LINALG',and='NEWORTHO')) THEN
                IF ( keyword_contains(line,'ON') ) THEN
                   cntl%tortho_new=.TRUE.
                ELSEIF ( keyword_contains(line,'OFF') ) THEN
                   cntl%tortho_new=.FALSE.
                ENDIF
                ! Block size of states distribution for CUDA newortho
             ELSEIF ( keyword_contains(line,'DISORTHO_BSIZE') ) THEN
                READ(iunit,'(A)',iostat=ierr) line
                CALL readsi(line,1,last,cnti%disortho_bsize,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                ! Block size of states distribution
             ELSEIF ( keyword_contains(line,'BLOCKSIZE',and='STATES')) THEN
                READ(iunit,'(A)',iostat=ierr) line
                CALL readsi(line,1,last,cnti%nstblk,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF ( keyword_contains(line,'ALLTOALL') ) THEN
                ! Use single/double precision in ALLTOALL communicaton
                IF ( keyword_contains(line,'SINGLE')) cntl%tr4a2a=.TRUE.
                IF ( keyword_contains(line,'DOUBLE')) cntl%tr4a2a=.FALSE.
             ELSEIF ( keyword_contains(line,'GSHELL') ) THEN
                ! Write file with GSHELL information
                cntl%gshellout=.TRUE.
                ! Empirical VDW correction
             ELSEIF ( keyword_contains(line,'VDW',and='CORRECTION') ) THEN
                IF ( keyword_contains(line,'OFF') ) THEN
                   vdwl%vdwc=.FALSE.
                ELSE
                   IF (.NOT.cntl%ksener) THEN
                      vdwl%vdwc=.TRUE. 
                   ELSE
                      CALL stopgm(procedureN,'vdWc incompatible with KS', &
                                  __LINE__,__FILE__)
                   ENDIF
                ENDIF
                ! Wannier VDW correction
             ELSEIF ( keyword_contains(line,'VDW',and='WANNIER') ) THEN
                IF ( keyword_contains(line,'OFF') ) THEN
                   vdwl%vdwd=.FALSE.
                ELSE
                   vdwl%vdwd=.TRUE.
                ENDIF
                ! Printing local potential
             ELSEIF ( keyword_contains(line,'LOCAL',and='POTENTIAL') ) THEN
                locpot2%tlpot=.TRUE.! cmb-kk -Printing local potential
             ELSEIF ( keyword_contains(line,'DCACP') ) THEN
                vdwl%dcacp = .TRUE.
             ELSEIF ( keyword_contains(line,'MIMIC') ) THEN
#ifdef __MIMIC
                cntl%mimic = .TRUE.
#else
                something_went_wrong = .true.
                error_message        = 'MIMIC IS REQUESTED BUT CPMD IS NOT COMPILED WITH IT!'
#endif
             ELSEIF ( keyword_contains(line,'NEW', and='CONSTRAINTS')) THEN
                cntl%new_constraints = .TRUE.
                IF ( keyword_contains(line,'PBICGSTAB') ) cntl%pbicgstab = .TRUE. !Use pbicgstab (default is pcg)
             ELSEIF ( keyword_contains(line,'SHAKE_MAXSTEP',alias='SHAKE_MAXSTEPS')) THEN
                ! Maximum number of shake iterations (cnti%shake_maxstep)
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                CALL readsi(line ,1,last,cnti%shake_maxstep,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF ( keyword_contains(line,'SHAKE_CG_ITER')) THEN
                ! Maximum number of pcg steps in each shake iteration (cnti%shake_cg_iter)
                previous_line = line
                READ(iunit,'(A)',iostat=ierr) line
                CALL readsi(line ,1,last,cnti%shake_cg_iter,erread)
                IF (erread) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSE
                ! Unknown keyword
                IF (' '/=line) THEN
                   IF (nbr_unknown_lines < max_unknown_lines) THEN
                      nbr_unknown_lines = nbr_unknown_lines+1
                      unknown(nbr_unknown_lines)=line
                   ELSE
                      DO i=1,max_unknown_lines-1
                         unknown(i) = unknown(i+1)
                      ENDDO
                      unknown(nbr_unknown_lines)=line
                   ENDIF
                ENDIF
             ENDIF
          END DO
          !
          ! End of read-loop
          !
       ELSE
          something_went_wrong = .true.
          error_message        = 'MISSING &CPMD SECTION - SECTION MANDATORY'
       ENDIF
       !
       IF (something_went_wrong) THEN
           WRITE(output_unit,'(/,1X,64("!"))')
           WRITE(output_unit,'(1X, A, 1X, A)') 'ERROR:', 'PROBLEM WHILE READING &CPMD SECTION:'
           WRITE(output_unit,'(8X, A)') TRIM(ADJUSTL(error_message))
           IF (line /= ' ' .or. previous_line /= ' ') THEN
              WRITE(output_unit,'(8X, A)') 'THE LAST TWO LINES READ WITHIN THE SECTION WERE:'
              WRITE(output_unit,'(/,1X, A)') TRIM(ADJUSTL(previous_line))
              WRITE(output_unit,'(1X, A)') TRIM(ADJUSTL(line))
           ENDIF
           WRITE(output_unit,'(1X,64("!"))')
           CALL stopgm(procedureN,'Error while reading &CPMD section, cf. output file',&
                __LINE__,__FILE__)
       ENDIF
  
       ! FileOpen settings
       !
       CALL set_fo_info()
       !
       ! Write the CPMD header here
       !
       CALL header(cnts%inputfile)
       !
       ! Print &INFO and MIRROR sections
       !
       CALL info_pri()
       CALL mirror_pri()
       !
       ! Check options and write to output
       ! 
       CALL control_test(tsrho)
       CALL control_pri(unknown,nbr_unknown_lines,fformat)
       
    ENDIF

    CALL control_bcast()
    
    ! NSPLP is broadcasted (CONTROL_BCAST).
    ! Set maxsys%mmaxx (maximum spline points). +20 to have NADD (setspline).
    ! 999 was the old value (cpmd 3.0)
    maxsys%mmaxx=MAX(cnti%nsplp+20,999)

    CONTAINS

    ! ==--------------------------------------------------------------==
    SUBROUTINE set_fo_info()

       CALL xstring(fo_info%fpath,fo_info%iapath,fo_info%iepath)
       IF (fo_info%fpath(fo_info%iepath:fo_info%iepath)/='/') THEN
          fo_info%iepath=fo_info%iepath+1
          fo_info%fpath(fo_info%iepath:fo_info%iepath)='/'
       ENDIF

    END SUBROUTINE set_fo_info
    ! ==--------------------------------------------------------------==
    SUBROUTINE info_pri()

       INTEGER :: i

       IF (nbr_info_lines>0) THEN
          WRITE(output_unit,'(/,1X,78("*"))')
          WRITE(output_unit,'(" * ",10("INFO - "),"INFO *")')
          WRITE(output_unit,'(1X,78("*"))')
          DO i=1,nbr_info_lines,1
             WRITE(output_unit,'(" * ",A74," *")') infomsg(i)
          ENDDO
          WRITE(output_unit,'(1X,78("*"))')
       ENDIF
       
    END SUBROUTINE info_pri
    ! ==--------------------------------------------------------------==
    SUBROUTINE mirror_pri()

       INTEGER :: ierr
       LOGICAL :: go_on_reading

       IF (mirror) THEN
          REWIND(iunit)
          WRITE(output_unit,'(/,1X,78("*"))')
          WRITE(output_unit,'(" ** ",T33,A,T77," ** ")') ' INPUT FILE'
          WRITE(output_unit,'(1X,78("*"))')
          go_on_reading = .true.
          DO WHILE(go_on_reading)
             READ(iunit,'(A80)',iostat=ierr) line
             IF (ierr /= 0) THEN
                go_on_reading = .false.
             ELSE
                WRITE(output_unit,'(A4,A72,A4)')' ** ',line,' **'
             ENDIF
          ENDDO
          WRITE(output_unit,'(1X,78("*"))')
          WRITE(output_unit,'(1X,78("*"))')
          REWIND(iunit)
       ENDIF

    END SUBROUTINE mirror_pri
    ! ==--------------------------------------------------------------==
  END SUBROUTINE control
  ! ==================================================================
  SUBROUTINE get_input_name(filename)
    ! ==--------------------------------------------------------------==
    CHARACTER(len=*), PARAMETER              :: procedureN='get_input_name'

    CHARACTER(len=*), INTENT(inout)          :: filename

    INTEGER                                  :: icarg

    icarg=m_iargc()
    IF (icarg<1) THEN
       IF (paral%io_parent)&
            WRITE(output_unit,*) ' NO INPUT FILE NAME SPECIFIED '
       CALL stopgm(procedureN,'No input file name specified',& 
            __LINE__,__FILE__)
    ELSE
       CALL m_getarg(1,filename)
    ENDIF
    ! ==--------------------------------------------------------------==
  END SUBROUTINE get_input_name
  !  ==================================================================
END MODULE control_utils
