MODULE respin_p_utils
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE inscan_utils,                    ONLY: inscan
  USE kinds,                           ONLY: real_8
  USE kpnt,                            ONLY: rk,&
                                             wk
  USE kpts,                            ONLY: kpts_com,&
                                             tkpts
  USE mp_interface,                    ONLY: mp_bcast
  USE parac,                           ONLY: parai,&
                                             paral
  USE readsr_utils,                    ONLY: readsi,&
                                             readsr,&
                                             xstring
  USE response_pmod,                   ONLY: &
       bdxvecx, bdxvecy, bdxvecz, btheta, bvecx, bvecy, bvecz, deltarvalues, &
       dmbi, dmbr, eig1, eig2, epr_options, iepr_csgt, iepr_default, &
       iepr_llc, iepr_novirt, iepr_wc, inmr_csgt, inmr_default, inmr_llc, &
       inmr_novirt, inmr_wc, kpertpar, lancphon, nf, nmr_options, nmr_para, &
       numf, ownpotvalue, response1, response2, tweight, voa_options, wghtf
  USE store_types,                     ONLY: restart1
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             nkpt
#include "sizeof.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: respin_p

CONTAINS

  SUBROUTINE respin_p
    ! ==--------------------------------------------------------------==
    ! ==  THIS ROUTINE READS THE SECTION &RESP &END ON UNIT IUNIT     ==
    ! ==--------------------------------------------------------------==
    ! ==  THE INPUT HAS TO BE IN THE FOLLOWING FORMAT                 ==
    ! ==     &RESP                                                    ==
    ! ==        Options                                               ==
    ! ==     &END                                                     ==
    ! ==--------------------------------------------------------------==
    ! ==  LIST OF OPTIONS                                             ==
    ! ==                                                              ==
    ! ==  KEEPREALSPACE                                               ==
    ! ==    like the standard CPMD option, this keeps the C0 ground   ==
    ! ==    state wavefunctions in the direct space representation    ==
    ! ==    during the calculation. Can save a lot of time, but is    ==
    ! ==    incredibly memory intensive.                             ==
    ! ==                                                              ==
    ! ==                                                              ==
    ! ==  CONVERGENCE                                                 ==
    ! ==    [0.00001]                                                 ==
    ! ==    convergence criterium on the gradient dE/dpsi             ==
    ! ==                                                              ==
    ! ==  NOOPT                                                       ==
    ! ==    Do not perform a ground state wfn optimization. Be sure   ==
    ! ==    the restarted wfn is at the BO-surface.                   ==
    ! ==                                                              ==
    ! ==  HTHRS or HAMiltonian CUToff                                 ==
    ! ==    [0.5]                                                     ==
    ! ==    The value where the preconditioner (the approximate       ==
    ! ==    Greens function (V_kin+V_coul-eps_KS)^-1) is truncated.   ==
    ! ==                                                              ==
    ! ==  POLAK                                                       ==
    ! ==    Uses the Polak-Ribiere formula for the conjugate          ==
    ! ==    gradient algorithm. Can be safer in the convergence.      ==
    ! ==                                                              ==
    ! ==  CG-ANALYTIC                                                 ==
    ! ==    The number of steps for which the step length in the      ==
    ! ==    conjugate gradient optimization is calculated assuming    ==
    ! ==    a quadratic functional E(2) (quadratic in C1). Default    ==
    ! ==    is 99, and 3 for NMR.                                     ==
    ! ==                                                              ==
    ! ==  CG-FACTOR                                                   ==
    ! ==    The analytic length calculation yields in general         ==
    ! ==    a result that is slightly too large. This factor          ==
    ! ==    is used to correct for that deficiency. Default is        ==
    ! ==    0.8.                                                      ==
    ! ==                                                              ==
    ! ==                                                              ==
    ! ==  TIGHTPREC                                                   ==
    ! ==    Uses a harder preconditioner. For experts: The            ==
    ! ==    Hamiltonian is approximated by the kinetic energy,        ==
    ! ==    the G-diagonal Coulomb potential and the KS-energies.     ==
    ! ==    The number obtained this way must not be close to zero.   ==
    ! ==    This is achieved by smoothing it with                     ==
    ! ==      x -> f(x) = sqrt(x^2 + eps^2) [default]                 ==
    ! ==    or                                                        ==
    ! ==      x -> f(x) = (x^2 + eps^2)/x   [this option]             ==
    ! ==    The HARD option conserves the sign of the approximate     ==
    ! ==    Hamiltonian whereas the default formula does never        ==
    ! ==                                                              ==
    ! ==    diverge.                                                  ==
    ! ==                                                              ==
    ! ==                                                              ==
    ! ==    RAMAN                                                     ==
    ! ==                                                              ==
    ! ==    HARDNESS                                                  ==
    ! ==                                                              ==
    ! ==    PHONON                                                    ==
    ! ==      lambda                                                  ==
    ! ==                                                              ==
    ! ==    LANCZOS  [CONTINUE,DETAILS]                               ==
    ! ==    lanczos_dim  iterations   conv_threshold                  ==
    ! ==    lanczos_dim= dimension of the vibrational d.o.f.          ==
    ! ==    iterations = no. of iterations desired for this run       ==
    ! ==    conv_threshold = threshold for convergence on eigenvectors==
    ! ==    CONTINUE = argument for continuig Lanczos diagonalization ==
    ! ==               from a previous run                            ==
    ! ==               (reads file LANCZOS_CONTINUE)                  ==
    ! ==    DETAILS  = argument for verbosity. prints a lot of stuff  ==
    ! ==                                                              ==
    ! ==    EIGEN                                                     ==
    ! ==      lanczos-dimemsion num_converged_states conv_threshold   ==
    ! ==                                                              ==
    ! ==      shift_constant                                          ==
    ! ==                                                              ==
    ! ==    FUKUI [N=NF] [COEFFICIENTS]                               ==
    ! ==      n [weight]                                              ==
    ! ==      ... (NF lines)                                          ==
    ! ==                                                              ==
    ! ==    NMR [options, see response_p.inc]                         ==
    ! ==       CURRENT [TCURRENT]:    Calculate the current density   ==
    ! ==         and the shift tensor using the current density       ==
    ! ==       VERBOSE [TVERBOSE]:    Print additional (technical)    ==
    ! ==         information                                          ==
    ! ==       PSI0 [TPRINTWFNS]:     Write the Wannier states to     ==
    ! ==         files named psi0-##.dens                             ==
    ! ==       RHO0 [TPRINTRHO ]:     Write the orbital densities     ==
    ! ==         to files named rho0-##.dens                          ==
    ! ==       SPLIT [TSPLITSTATES]:  Calculate the contribution of 
    ! ==         each electronic state (to the shift tensor) separately
    ! ==       ISO [TISOONLY]:        Calculate only diagonal values 
    ! ==         of the shifttensor (S_ii)
    ! ==       LLC [!TCOMMONLLC]:     Special option. EXPERIMENTAL.
    ! ==       NOVIRT [TNOVIRTUAL]:   Implies NOLOC, and uses the same 
    ! ==         virtual cell (equal to the physical cell) for all the
    ! ==         orbitals. (ONLY FOR MOLECULES! )
    ! ==       NOLOC [!TLOCALIZE]:    Do not use Wannier functions but 
    ! ==         canonical orbitals (ONLY FOR MOLECULES! Use of NOVIRT
    ! ==         instead of NOLOC is recommended.)
    ! ==       FAST [TFAST]:          Omit some time-consuming operations,
    ! ==         e.g. the loss in information when doing FFTs R->G
    ! ==       WANNIERCENTER [TWC]:   The center of the wavefunction-specific
    ! ==         box is put where psi0^2 reaches its maximum. Default is
    ! ==         that the LLC is placed where psi0^2 -> min.
    ! ==       NOSMOOTH [!TSMOOTHING]: By default, an (exp - x^2) 
    ! ==         is applied around the LLC. Switched off by this option.
    ! ==       TRIPACK/SIXPACK:       Distribute the perturbation 
    ! ==         calculations for p_i and L_i (the momentum and 
    ! ==         angular momentum operators applied to Phi_0)
    ! ==         on several proc groups.
    ! ==           TRIPACK: nmr_threads = 3
    ! ==           SIXPACK: nmr_threads = 6
    ! ==       IGLO:                  Do an IGLO-type calculation.
    ! ==                              (no longer supported)           ==
    ! ==       FULL:                 Calculate the errors due to the  ==
    ! ==         neglected current contribution Delta J (see paper)   ==
    ! ==       RESTART                                                ==
    ! ==         restart a previous (or interrupted) calculation      ==
    ! ==       ISOLATED                                               ==
    ! ==         equivalent to OVERLAP / 0.1                          ==
    ! ==                                                              ==
    ! ==    EPR                                                       ==
    ! ==      ...with the same options as for the NMR keyword.        ==
    ! ==      PLEASE NOTE that the EPR keyword has to be preceded     ==
    ! ==      by at least one space, due to overlap with other        ==
    ! ==      keywords (stateprec).                                   ==
    ! ==                                                              ==
    ! ==                                                              ==
    ! ==                                                              ==
    ! ==    SIGMA                                                     ==
    ! ==      [0.1] Enter a value for the spread of the Gaussian      ==
    ! ==      for the calculation of the virtual boxes (see also      ==
    ! ==      nmr_util_p.F, routine calc_lower_left). The LLC         ==
    ! ==      of the virtual cell is put where the integral of the    ==
    ! ==      orbital density over the Gaussian exp -(x-x0)^2/2sigma  ==
    ! ==      is smallest. Should be used with care (may affect the   ==
    ! ==      results if chosen incorrectly).                         ==
    ! ==                                                              ==
    ! ==    DISCARD [OFF,PARTIAL,TOTAL,LINEAR]                        ==
    ! ==      keyword for discarding trivial modes in vibrational     ==
    ! ==      analysis (boht PHONONS and LANCZOS)                     ==
    ! ==      OFF = argument for performing no projection             ==
    ! ==      PARTIAL = argument for projecting out only tranlsations ==
    ! ==                (this is the default)                         ==
    ! ==      TOTAL = argument for projecting both rotations and      ==
    ! ==              translations                                    ==
    ! ==      LINEAR = argument for projecting rotations in a linear  ==
    ! ==               molecule (not implemented yet)                 ==
    ! ==                                                              ==
    ! [marcella]
    ! ==     
    ! ==    KPERT[options, see response_p.inc]
    ! ==        sampling of the BZ                                    ==
    ! ==        MONKHORSTPACK [TMONKP]                                ==
    ! ==        NK1  NK2  NK3                                         ==
    ! ==        otherwise read K-POINT coordinates,                   ==
    ! ==                  {SCALE} for scalingwrt basis vectors in BZ  ==
    ! ==                  NKPTS                                       ==
    ! ==                  RK(1,IK) RK(2,IK) RK(3,IK) WK(IK)           ==
    ! ==                  ........                                    ==
    ! ==       Type of calculation                                    ==
    ! ==       OPTIMIZATION RESPONSE WFN:                             ==
    ! ==       WRITE_C1 {RESTART_CU1} C1 are written in RESTA_P       == 
    ! ==       BUILD_C00{TK_PREP_DIAG} a standard restart file is     ==
    ! ==                prepared where 2*NSTATE*NKPTS wfn are stored  ==
    ! ==                for each  k-point:
    ! ==                1 : NSTATE are ground state wfn : C0          ==
    ! ==                NSTATE+1 : 2*NSTATE are                       ==
    ! ==                          i(kx c1x + ky c1y + kz c1z)         ==
    ! ==                this RESTART should be used by standard CPMD  ==
    ! ==                with K-POINT option (same k vectors)          ==
    ! ==                for the calculation of the KS energies by a   ==
    ! ==                single diagonalization of the H matrix        ==
    ! ==                (option ONLYDIAG[TONLYDIAG] see sysin.F)      ==
    ! ==                {READ_C1}[TKREAD_C1] C1 are not calculated    ==
    ! ==                                    but read in RESTA_P       ==
    ! ==       NORESTART only calculation of C1 and E2                ==
    ! ==       HAMILTONIAN[TPERT_HAMILT]{READ_C1} C0 and C1 are used  ==
    ! ==                  for H(k) construction with the second order ==
    ! ==                  approximation, and KS energies are computed ==
    ! ==                  {{C0},{i(kxC1x+kyc1y+kxc1z)}} 2*NSTATE      ==
    ! [anatole]                                                       ==
    ! ==                                                              ==
    ! ==       DUMMY_OPT for the linking dummyatom optimization       ==
    ! ==          DUMPBD parabolic descent optimizer, don't trust it  ==
    ! ==          DUMSD  steepest descent, more trustworthy           ==
    ! ==          DUMN   Newton optimizer, not trustworthy at all     ==
    ! ==          DUMMY_STEP asks for step-sizes for the steepest     ==
    ! ==                     descent for each of the 5 sigma of a     ==
    ! ==                     1996 SG-PP of the 2nd row                ==
    ! ==                                                              ==
    ! ==       DUMMY_REF will print out a reference electron density  ==
    ! ==                 and coordinates                              ==
    ! ==                                                              ==
    ! ==       DUMMY_TEC will just print out coordinates              ==
    ! ==                                                              ==
    ! ==    VOA [AT MD POSITiON_FORM DEBUG VERBOSE AMATRIX CURRENT]   ==
    ! ==       ATOMLIST 0001 0002 ...                                 ==
    ! ==                                                              ==
    ! ==--------------------------------------------------------------==
    ! Variables
    CHARACTER(*), PARAMETER                  :: procedureN = 'respin_p'

    CHARACTER(len=80)                        :: line, newline
    INTEGER                                  :: i, i1, ia, ie, ierr, ik, &
                                                iout, iresponse, is, iunit, &
                                                nstate
    LOGICAL                                  :: erread, read_more

! ==--------------------------------------------------------------==
! ==  DEFAULTS                                                    ==
! ==--------------------------------------------------------------==

    iunit=cnti%insys
    nstate=crge%n
    response1%response_running = .FALSE. ! this is for WRENER only.
    response1%tkeeprealspacewfn = .FALSE.
    response1%t_restart_c1 = .FALSE.
    response1%t_initialguess_c1 = .FALSE.

    ! anatoles variables:
    response1%tdummyatom=.FALSE.
    response1%tdummyatom_ref=.FALSE.
    response1%tdummyatom_tec=.FALSE.
    response1%toacp4f = .FALSE.
    response1%tdummy_step=.FALSE.
    response1%tdumpbd=.FALSE.
    response1%tdumn=.FALSE.
    response1%tdumsd=.FALSE.

    ! TYPE of the perturbation:
    response1%traman=.FALSE.
    response1%phonon=.FALSE.
    response1%teigensystem=.FALSE.
    response1%thardness=.FALSE.
    response1%tnmr=.FALSE.
    response1%tepr=.FALSE.
    response1%pr_energy=.TRUE.
    response1%tkpert=.FALSE.
    response1%tfukui=.FALSE.
    response1%tvoa=.FALSE.

    ! INTERACTION routines, initialization
    response1%tinteraction = .FALSE.
    dmbi%cutoff_trick = .FALSE.
    dmbi%cutoff_restr = .FALSE.
    dmbi%tlinearscaling = .FALSE.
    ! Default number of SCF iterations
    dmbi%bptscfiter=1
    ! initialising the Wannier rotation scalar variables (DMB)
    dmbi%blmax=4
    dmbi%nmol=1
    ! initialisation flags
    dmbi%in_btheta=.FALSE.
    dmbi%in_bvec=.FALSE.
    dmbi%in_bdxvec=.FALSE.
    ! By default no save, no load, full rotations, 
    ! single reference file, no gaussian representation
    dmbi%wann_save=.FALSE.
    dmbi%wann_load=.FALSE.
    dmbi%wann_allrot=.TRUE.
    dmbi%wann_multi=.FALSE.
    dmbi%wann_gauss=.FALSE.
    ! by default only one localisation cycle if WC prediction is on...
    dmbi%wc_pred_step=1000000

    ! d2Exc/dn1 dn2 calculation using LDA (FALSE by default)
    dmbi%fastdexc=.FALSE.

    ! by default there is a non-zero force inducing a response
    dmbi%tnoresponse=.FALSE.

    ! by default the tolerance on the SCF-PT cycles is 1E-6
    dmbr%scf_tol=1.e-6_real_8

    ! By default, W0 is not saved in a RESTART file
    dmbi%save_localised_w0=.FALSE.

    ! By default use wannier functions for W0 and not atomic wf
    dmbi%tatomicwavefunction=.FALSE.

    ! By default use wannier density and not atomic one
    dmbi%tatomicrho=.FALSE.

    ! By default the Wannier orbitals are orthogonalised
    dmbi%torthog_wannier=.TRUE.

    ! By default full PT calculation
    dmbi%tsimple_model=.FALSE.
    dmbi%tinter_ortho = .FALSE.

    ! BornOppenheimer dynamics.. 
    ! By default no predictions of the WC position
    dmbi%wcpredict=.FALSE.
    dmbi%wc_sample=99

    ! Maximum number of states saved at the end of PT-INTERACTION
    dmbi%max_num_disp=10

    ! use parallel version by default
    dmbi%tmanno = .FALSE.
    ! end of initialisation (DMB)

    ! wfn optimization method:
    response1%tsde_p=.FALSE.
    response1%prec_p=.TRUE.
    response1%pcg_p=.TRUE.
    response1%pcgmin_p=.TRUE.
    response1%tpolak_p=.FALSE.
    response2%hthrs_p=0.5_real_8
    response2%tolog_p=0.00001_real_8
    response2%cg_factor=0.8_real_8
    response1%preconditioner_p=1
    response1%cg_analytic=99
    response1%opt1st = .FALSE.
    response1%mdiis_p = 10

    ! Lanczos Phonon option:
    response1%tlanphon=.FALSE.
    response1%tlanph_cont=.FALSE.
    lancphon%details=.FALSE.

    ! ... discarding defaults
    response1%projout=.TRUE.
    response1%rotout=.FALSE.



    response1%teigens=.FALSE.


    ! NMR special options:
    nmr_options%tcurrent=.FALSE.
    nmr_options%tverbose=.FALSE.
    nmr_options%tprintwfns=.FALSE.
    nmr_options%tprintrho=.FALSE.
    nmr_options%tlocalize=.TRUE.
    nmr_options%tfast=.FALSE.
    nmr_options%tsmoothing=.TRUE.
    deltarvalues%sigma=0.1_real_8
    nmr_para%nmr_threads=1
    nmr_options%inmr_method=inmr_default
    nmr_options%inmr_virtual=inmr_default
    nmr_options%tnmr_full = .FALSE.
    nmr_options%tnmr_only_error = .FALSE.
    nmr_options%tnmr_overlaps= .FALSE.
    nmr_options%tforce_xyz(1)= .FALSE.
    nmr_options%tforce_xyz(2)= .FALSE.
    nmr_options%tforce_xyz(3)= .FALSE.
    nmr_para%nmr_superparent   = paral%parent
    nmr_para%nmr_supergroup    = parai%cp_grp
    nmr_para%nmr_total_nproc   = parai%cp_nproc
    nmr_para%nmr_supersource   = parai%source

    ! EPR special options:
    epr_options%teverbose=.FALSE.
    epr_options%teprintwfns=.FALSE.
    epr_options%teprintrho=.FALSE.
    epr_options%telocalize=.TRUE.
    epr_options%tefast=.FALSE.
    epr_options%tesmoothing=.TRUE.
    epr_options%iepr_method=iepr_default
    epr_options%iepr_virtual=iepr_default
    epr_options%tepr_full = .FALSE.
    epr_options%tepr_only_error = .FALSE.
    epr_options%tepr_overlaps= .FALSE.
    epr_options%teforce_xyz(1)= .FALSE.
    epr_options%teforce_xyz(2)= .FALSE.
    epr_options%teforce_xyz(3)= .FALSE.
    epr_options%tepr_smart = .FALSE.
    epr_options%tepr_hyp = .FALSE.
    epr_options%tepr_ownpot = .FALSE.

    ! KPERT special option
    kpertpar%restart_cu1 = .FALSE.
    kpertpar%tk_prep_diag =.FALSE.
    kpertpar%tkread_c1 =.FALSE.
    kpertpar%tpert_hamilt =.FALSE.

    ! Num Orbitals for Fukui
    numf=0
    erread = .FALSE.

    IF (paral%parent) THEN
       ! ==------------------------------------------------------------==
       read_more = (inscan(iunit,'&RESP') .EQ. 0)
       ! read more [data from input file] if &RESP was found.

       DO WHILE (read_more)

          IF (paral%io_parent)&
               READ(iunit,err=99,END=99,fmt='(A80)') line
          read_more = (INDEX(line,'&END').EQ.0)
          ! read more if &END is not part of the current line.

          IF (INDEX(line,'RESTART').NE.0) THEN
             ! A RESTART implies that C0 is not re-optimized!
             response1%t_restart_c1 = .TRUE.
             cntl%wfopt=.FALSE.
             restart1%restart=.TRUE.
             restart1%rwf=.TRUE.
          ENDIF

          ! anatoles calls:
          IF (INDEX(line,'OACP DENSITY').NE.0) THEN
             response1%tdummyatom   = .TRUE.
             IF (paral%io_parent)&
                  READ(iunit,err=99,END=99,fmt='(A80)') newline
             CALL xstring(newline,ia,ie)
             CALL readsi(newline,ia,ie,response2%tdumnum,erread)
             IF (erread) GOTO 99
             ia=ie
             CALL readsr(newline,ia,ie,response2%dumcrit,erread)
             IF (erread) GOTO 99

             IF (paral%io_parent)&
                  READ(iunit,err=99,END=99,fmt='(A80)') newline
             CALL xstring(newline,ia,ie)
             CALL readsr(newline,ia,ie,response2%dumstep1,erread)
             IF (erread) GOTO 99
             ia=ie
             CALL readsr(newline,ia,ie,response2%dumstep2,erread)
             IF (erread) GOTO 99
             ia=ie
             CALL readsr(newline,ia,ie,response2%dumstep3,erread)
             IF (erread) GOTO 99
             ia=ie
             CALL readsr(newline,ia,ie,response2%dumstep4,erread)
             IF (erread) GOTO 99
             ia=ie
             CALL readsr(newline,ia,ie,response2%dumstep5,erread)
             IF (erread) GOTO 99

             IF (paral%io_parent)&
                  READ(iunit,err=99,END=99,fmt='(A80)') newline
             CALL xstring(newline,ia,ie)
             CALL readsr(newline,ia,ie,response2%rad_dum,erread)
             IF (erread) GOTO 99
             ia=ie
             CALL readsr(newline,ia,ie,response2%rad_norm,erread)
             IF (erread) GOTO 99
          ENDIF
          IF (INDEX(line,'OACP REF_DENSITY').NE.0)&
               response1%tdummyatom_ref = .TRUE.
          IF (INDEX(line,'OACP FORCE').NE.0) THEN
             response1%toacp4f = .TRUE.
             IF (paral%io_parent)&
                  READ(iunit,err=99,END=99,fmt=*)&
                  response2%tdumnum,response2%dumcrit,response2%tfdimen
             IF (paral%io_parent)&
                  READ(iunit,err=99,END=99,fmt=*)&
                  response2%ang_mom,response2%dumstep4,response2%dumstep5
          ENDIF

          IF (INDEX(line,'RAMAN').NE.0)       response1%traman       = .TRUE.
          IF (INDEX(line,'PHONON').NE.0)      response1%phonon       = .TRUE.
          IF (INDEX(line,'EIGEN').NE.0) THEN
             response1%teigensystem = .TRUE.
             IF (paral%io_parent)&
                  READ(iunit,err=99,END=99,fmt=*)&
                  eig1%lanczos_dim,eig1%num_lanczos_states,eig2%conv_threshold,&
                  eig2%eigenvalue_shift
          ENDIF
          IF (INDEX(line,'HARDNESS').NE.0) THEN
             response1%thardness    = .TRUE.
             IF (paral%io_parent)&
                  READ(iunit,err=99,END=99,fmt=*)&
                  eig1%lanczos_dim,eig1%num_lanczos_states,eig2%conv_threshold,&
                  eig2%eigenvalue_shift
          ENDIF

          IF (INDEX(line,'FUKUI').NE.0)  THEN
             response1%tfukui    = .TRUE.
             tweight   = .TRUE.
             IF (INDEX(line,'COEF').NE.0) tweight=.FALSE.
             i1=INDEX(line,'=')
             IF (i1.NE.0) THEN
                CALL readsi(line,i1+1,iout,numf,erread)
             ELSE
                numf=1
             ENDIF
             IF (erread) GOTO 99
             ALLOCATE(nf(numf),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             ALLOCATE(wghtf(numf+1),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             DO i=1,numf
                IF (paral%io_parent)&
                     READ(iunit,err=99,END=99,fmt='(A80)') newline
                CALL readsi(newline,1,iout,nf(i),erread)
                IF (erread) GOTO 99
                IF (nf(i).GE.0)&
                     CALL stopgm('RESPIN_P','NF must be negative', &
                     __LINE__,__FILE__)
                CALL readsr(newline,iout,i1,wghtf(i),erread)
                IF (erread) wghtf(i)=1.0_real_8
                ! READ(IUNIT,ERR=99,END=99,FMT=*) NF(I),WGHTF(I)
             ENDDO
             wghtf(numf+1)=1.0_real_8
          ENDIF
          IF (INDEX(line,'INTERACT').NE.0)    response1%tinteraction = .TRUE.
          IF (INDEX(line,'KPERT').NE.0)      THEN
             IF (paral%io_parent)&
                  WRITE(6,'("************* KPERT RESPONSE ********")')
             IF (paral%io_parent)&
                  WRITE(6,'(">>>>>>>>>>>>>>>   READ THE KPERT OPTIONS")')
             response1%tkpert       = .TRUE.
             IF (INDEX(line,'MONKH') .NE. 0) THEN
                tkpts%tmonkp = .TRUE.
                IF (paral%io_parent)&
                     READ(iunit,err=99,END=99,fmt=*) kpts_com%nk1,kpts_com%nk2,kpts_com%nk3
                IF (paral%io_parent)&
                     WRITE(6,&
                     '("MONKHORST_PACK MESH ",i4,"X",i4,"X",i4)')&
                     kpts_com%nk1,kpts_com%nk2,kpts_com%nk3
             ELSE
                IF (paral%io_parent)&
                     WRITE(6,'("CHOSEN SET OF K VECTORS")')
                IF (INDEX(line,'SCALE') .NE. 0) THEN
                   tkpts%tkscale =.TRUE.
                ENDIF
                IF (paral%io_parent)&
                     READ(iunit,err=99,END=99,fmt=*) nkpt%nkpts
                ! 
                IF (paral%io_parent)&
                     WRITE(6,'("NKPTS=",i5)') nkpt%nkpts
                ! 
                IF (nkpt%nkpts.LE.0)&
                     CALL stopgm&
                     ('  RESPIN','WRONG NUMBER OF K POINTS',& 
                     __LINE__,__FILE__)
                ALLOCATE(rk(3,nkpt%nkpts),STAT=ierr)
                IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                     __LINE__,__FILE__)
                ALLOCATE(wk(nkpt%nkpts),STAT=ierr)
                IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                     __LINE__,__FILE__)
                DO ik=1,nkpt%nkpts
                   IF (paral%io_parent)&
                        READ(iunit,err=99,END=99,fmt='(A)') line
                   ia=1
                   DO i=1,3
                      CALL readsr(line,ia,ie,rk(i,ik),erread)
                      IF (erread) THEN
                         IF (paral%io_parent)&
                              WRITE(6,&
                              '(" KPOINTS: RK(",I1,",",I4,")=")')&
                              i,ik
                         IF (paral%io_parent)&
                              WRITE(6,*) 'ERROR READING KPOINTS INPUT'
                         GOTO 99
                      ENDIF
                      ia=ie
                   ENDDO
                   CALL readsr(line,ia,ie,wk(ik),erread)
                   IF (erread) THEN
                      IF (paral%io_parent)&
                           WRITE(6,'(" KPOINTS: WEIGHT OF KPOINT ",I4)')&
                           ik
                      IF (paral%io_parent)&
                           WRITE(6,*) 'ERROR READING KPOINTS INPUT'
                      GOTO 99
                   ENDIF
                ENDDO
             ENDIF
75           CONTINUE
             IF (paral%io_parent)&
                  READ(iunit,err=99,END=99,fmt='(A80)') line
             IF (INDEX(line,'WRITE_C1') .NE. 0) THEN
                IF (paral%io_parent)&
                     WRITE(6,&
                     '(">>>>>>>>>> AT THE END, WRITE THE 3 C1 SETS",'&
                     //'" ON 3 SEPARATED RESTART FILES")')
                kpertpar%restart_cu1 = .TRUE.
                GOTO 75
             ELSEIF (INDEX(line,'BUILD_C00') .NE. 0) THEN
                IF (paral%io_parent)&
                     WRITE(6,&
                     '(">>>>>>>>>> USE THE C0 AND THE C1 TO BUILD THE C00")')
                kpertpar%tk_prep_diag =.TRUE.
                IF (INDEX(line,'READ_C1').NE. 0) THEN
                   IF (paral%io_parent)&
                        WRITE(6,'(">>>>>>>>>> THE C1 ARE READ FROM RESTA_P")')
                   kpertpar%tkread_c1 = .TRUE.
                ENDIF
                GOTO 75
             ELSEIF (INDEX(line,'HAMILTONIAN') .NE. 0) THEN
                IF (paral%io_parent)&
                     WRITE(6,&
                     '(">>>>>>>>>> BUILD H BY USING THE PERT. THEORY")')
                IF (paral%io_parent)&
                     WRITE(6,&
                     '(">>>>>>>>>> H IS APPROXIMATED TO SECOND ORDER IN K")')
                kpertpar%tpert_hamilt =.TRUE.
                IF (INDEX(line,'READ_C1').NE. 0) THEN
                   IF (paral%io_parent)&
                        WRITE(6,&
                        '(">>>>>>>>>> THE C1 ARE READ FROM RESTA_P")')
                   kpertpar%tkread_c1 = .TRUE.
                ENDIF
                GOTO 75
             ELSEIF (INDEX(line,'NORESTART') .NE. 0) THEN
                kpertpar%restart_cu1 = .FALSE.
                kpertpar%tk_prep_diag =.FALSE.
                IF (paral%io_parent)&
                     WRITE(6,'(">>>>>>>>>> NO RESTART FILE IS WRITTEN")')
                GOTO 75
             ELSE
                GOTO 76
             ENDIF
76           CONTINUE
             IF (kpertpar%tk_prep_diag .AND. kpertpar%tpert_hamilt) CALL stopgm&
                  ('  RESPIN',' 2 MANY METHODS',& 
                  __LINE__,__FILE__)
          ENDIF
          IF (INDEX(line,'LANCZOS').NE.0)      THEN
             response1%tlanphon=.TRUE.
             IF (INDEX(line,'CONT').NE.0) response1%tlanph_cont=.TRUE.
             IF (paral%io_parent)&
                  READ(iunit,err=99,END=99,fmt=*) eig1%lanczos_dim,&
                  lancphon%nlan_st,eig2%conv_threshold
             eig1%lanczos_dim=eig1%lanczos_dim+1
          ENDIF
          ! ... VERBOSITY IN LANCZOS ROUTINE
          IF (INDEX(line,'DETAILS').NE.0) lancphon%details=.TRUE.
          ! ... DISCARDING TRANSLATIONAL AND/OR ROTATIONAL DEGREES OF FREEDOM
          IF (INDEX(line,'DISCARD').NE.0) THEN
             IF (INDEX(line,'OFF').NE.0) THEN
                response1%projout=.FALSE.
                response1%rotout=.FALSE.
             ELSEIF (INDEX(line,'PARTIAL').NE.0) THEN
                response1%projout=.TRUE.
                response1%rotout=.FALSE.
             ELSEIF (INDEX(line,'TOTAL').NE.0) THEN
                response1%projout=.TRUE.
                response1%rotout=.TRUE.
             ELSEIF (INDEX(line,'LINEAR').NE.0) THEN
                IF (paral%io_parent)&
                     WRITE(6,'(A)') ' LINEAR MOLECULE NOT IMPLEMENTED YET '
                CALL stopgm('RESPIN_P',' KEYWORD LINEAR NOT'//&
                     '  IMPLEMENTED YET ',& 
                     __LINE__,__FILE__)
             ELSE
                IF (paral%io_parent)&
                     WRITE(6,'(A)') ' KEYWORD DISCARD NEEDS AN ARGUMENT '
                CALL stopgm('RESPIN_P',' KEYWORD DISCARD NEEDS AN'&
                     //'ARGUMENT ',& 
                     __LINE__,__FILE__)
             ENDIF
          ENDIF


          ! ==--  NMR  ---------------------------------------------------==
          IF (INDEX(line,'NMR ').NE.0)         THEN
             response1%tnmr         = .TRUE.
             ! IF(INDEX(LINE,' IGLO').NE.0)     INMR_METHOD  = INMR_IGLO
             ! No longer supported.
             IF (INDEX(line,' CSGT').NE.0)     nmr_options%inmr_method  = inmr_csgt
             IF (INDEX(line,' FAST').NE.0)     nmr_options%tfast        = .TRUE.
             IF (INDEX(line,' NOLOC').NE.0)    nmr_options%tlocalize    = .FALSE.
             IF (INDEX(line,' NOVIRT').NE.0)&
                  nmr_options%inmr_virtual = inmr_novirt
             IF (INDEX(line,' WANNIER').NE.0)  nmr_options%inmr_virtual = inmr_wc
             IF (INDEX(line,' LLC').NE.0)      nmr_options%inmr_virtual = inmr_llc
             IF (INDEX(line,' FORCEX').NE.0)   nmr_options%tforce_xyz(1)= .TRUE.
             IF (INDEX(line,' FORCEY').NE.0)   nmr_options%tforce_xyz(2)= .TRUE.
             IF (INDEX(line,' FORCEZ').NE.0)   nmr_options%tforce_xyz(3)= .TRUE.
             IF (INDEX(line,' CURRENT').NE.0)  nmr_options%tcurrent     = .TRUE.
             IF (INDEX(line,' VERBOSE').NE.0)  nmr_options%tverbose     = .TRUE.
             IF (INDEX(line,' PSI0').NE.0)     nmr_options%tprintwfns   = .TRUE.
             IF (INDEX(line,' RHO0').NE.0)     nmr_options%tprintrho    = .TRUE.
             IF (INDEX(line,' NOSMOOTH').NE.0) nmr_options%tsmoothing   = .FALSE.
             IF (INDEX(line,' TRIPACK').NE.0)  nmr_para%nmr_threads  = 3
             IF (INDEX(line,' SIXPACK').NE.0)  nmr_para%nmr_threads  = 6
             response1%cg_analytic=3
             IF (nmr_para%nmr_threads .NE. 1) THEN
                cntl%wfopt=.FALSE.
                restart1%restart=.TRUE.
                restart1%rwf=.TRUE.
             ENDIF
          ENDIF           ! NMR-line
          IF (INDEX(line,'ISOLAT').NE.0) THEN
             nmr_options%tnmr_overlaps= .TRUE.
             deltarvalues%overlap_threashold=0.1_real_8
             nmr_options%tnmr_full = .FALSE.
          ENDIF
          IF (INDEX(line,'OVERLAP').NE.0) THEN
             nmr_options%tnmr_overlaps= .TRUE.
             IF (paral%io_parent)&
                  READ(iunit,err=99,END=99,fmt='(A80)') newline
             CALL readsr(newline ,1,iout,deltarvalues%overlap_threashold,erread)
             IF (erread)&
                  CALL stopgm('RESPIN','ERROR reading ovl threashold',& 
                  __LINE__,__FILE__)
          ENDIF
          IF (INDEX(line,'FULL').NE.0)    THEN
             nmr_options%tnmr_full   = .TRUE.
          ENDIF
          ! NMR END

          ! ==--  VOA  ---------------------------------------------------==
          IF (INDEX(line,'VOA ').NE.0) THEN
             response1%tvoa = .TRUE.
             voa_options%atomlistline = ' '
             IF (INDEX(line,' DEBUG').NE.0) voa_options%tdebug = .TRUE.
             IF (INDEX(line,' VERBOSE').NE.0) voa_options%tverbose = .TRUE.
             IF (INDEX(line,' MD').NE.0) voa_options%tmd = .TRUE.
             IF (INDEX(line,' AT').NE.0) voa_options%tat = .TRUE.
             IF (INDEX(line,' POSITION_FORM').NE.0) voa_options%tpf = .TRUE.
             IF (INDEX(line,' CURRENT').NE.0) voa_options%tcurrent = .TRUE.
             IF (INDEX(line,' DENSITY').NE.0) voa_options%tdensity = .TRUE.
             IF (INDEX(line,' AMATRIX').NE.0) voa_options%tamat = .TRUE.
          ENDIF ! VOA-line
          IF(voa_options%tdebug) voa_options%tverbose = .TRUE.
          IF (INDEX(line,' ATOMLIST').NE.0) THEN
             voa_options%tatomlist = .TRUE.
             voa_options%atomlistline=TRIM(line)
          ENDIF
          ! VOA END
          ! ==--  EPR  ---------------------------------------------------==
          IF (INDEX(line,' EPR').NE.0)         THEN
             response1%tepr         = .TRUE.
             IF (INDEX(line,' CSGT').NE.0)     &
                  epr_options%iepr_method   = iepr_csgt
             IF (INDEX(line,' FAST').NE.0)     epr_options%tefast        = .TRUE.
             IF (INDEX(line,' NOLOC').NE.0)    epr_options%telocalize    = .FALSE.
             IF (INDEX(line,' NOVIRT').NE.0)&
                  epr_options%iepr_virtual = iepr_novirt
             IF (INDEX(line,' WANNIER').NE.0)  epr_options%iepr_virtual  = iepr_wc
             IF (INDEX(line,' LLC').NE.0)      epr_options%iepr_virtual  = iepr_llc
             IF (INDEX(line,' FORCEX').NE.0)   epr_options%teforce_xyz(1)= .TRUE.
             IF (INDEX(line,' FORCEY').NE.0)   epr_options%teforce_xyz(2)= .TRUE.
             IF (INDEX(line,' FORCEZ').NE.0)   epr_options%teforce_xyz(3)= .TRUE.
             IF (INDEX(line,' VERBOSE').NE.0)  epr_options%teverbose     = .TRUE.
             IF (INDEX(line,' PSI0').NE.0)     epr_options%teprintwfns   = .TRUE.
             IF (INDEX(line,' RHO0').NE.0)     epr_options%teprintrho    = .TRUE.
             IF (INDEX(line,' NOSMOOTH').NE.0) epr_options%tesmoothing   = .FALSE.
             IF (INDEX(line,' HYP').NE.0)      epr_options%tepr_hyp   = .TRUE.
             response1%cg_analytic=99
             ! Synchronize NMR/EPR options
             nmr_options%inmr_method=epr_options%iepr_method
             nmr_options%tfast=epr_options%tefast
             nmr_options%tlocalize=epr_options%telocalize
             nmr_options%inmr_virtual=epr_options%iepr_virtual
             nmr_options%tforce_xyz(1)=epr_options%teforce_xyz(1)
             nmr_options%tforce_xyz(2)=epr_options%teforce_xyz(2)
             nmr_options%tforce_xyz(3)=epr_options%teforce_xyz(3)
             nmr_options%tverbose=epr_options%teverbose
             nmr_options%tprintwfns=epr_options%teprintwfns
             nmr_options%tprintrho=epr_options%teprintrho
             nmr_options%tsmoothing=epr_options%tesmoothing
          ENDIF           ! EPR-line
          IF (INDEX(line,'OWNPOT').NE.0) THEN
             epr_options%tepr_ownpot = .TRUE.
             IF (paral%io_parent)&
                  READ(iunit,err=99,END=99,fmt='(A80)') newline
             CALL readsr(newline ,1,iout,ownpotvalue%epr_ownpot_value,erread)
             IF (erread)&
                  CALL stopgm('RESPIN','ERROR reading potential smoothing value',& 
                  __LINE__,__FILE__)
             IF (ownpotvalue%epr_ownpot_value.LE.0.0_real_8)&
                  CALL stopgm('RESPIN','ERROR reading potential smoothing value (must be > zero)',& 
                  __LINE__,__FILE__)
          ENDIF
          IF (INDEX(line,'SMART').NE.0) THEN
             epr_options%tepr_smart = .TRUE.
             IF (paral%io_parent)&
                  READ(iunit,err=99,END=99,fmt='(A80)') newline
             CALL readsr(newline ,1,iout,ownpotvalue%epr_smart_value,erread)
             IF (erread)&
                  CALL stopgm('RESPIN','ERROR reading threshold value in FULL calculation',& 
                  __LINE__,__FILE__)
             IF ((ownpotvalue%epr_smart_value.LT.0.0_real_8).OR.(ownpotvalue%epr_smart_value.GT.0.1_real_8))&
                  CALL stopgm('RESPIN','ERROR reading threshold value in FULL calculation (must be >= 0 and <= 1)',& 
                  __LINE__,__FILE__)
          ENDIF

          ! These are already read in the NMR section
          epr_options%tepr_overlaps = nmr_options%tnmr_overlaps
          epr_options%tepr_full = nmr_options%tnmr_full
          ! EPR END


          ! ==------------------------------------------------------------==
          IF (INDEX(line,'DETAILS').NE.0) lancphon%details=.TRUE.


          ! ==--  INTERACTION  -------------------------------------------==
          IF (response1%tinteraction) THEN
             IF (INDEX(line,' LINEARSCALING').NE.0) THEN
                ! We impose only that the response orbitals be orthogonal 
                ! to their "own" ground state orbital, i.e. <psi (1)_i | psi(0)_i> = 0 
                ! but not for other psi(0)_j.
                dmbi%tlinearscaling = .TRUE.
             ENDIF

             IF (INDEX(line,' ORTHOGONAL').NE.0) THEN
                ! Then, we want the perturbation orbitals to be orthogonal 
                ! between themselves.
                dmbi%tinter_ortho = .TRUE.
             ENDIF

             IF (INDEX(line,'LMAX').NE.0) THEN
                ! reads in the maximum angular momentum used for rotations of Wannier functions
                IF (paral%io_parent)&
                     READ(iunit,err=99,END=99,fmt=*) dmbi%blmax
             ENDIF

             IF (INDEX(line,'NMOL').NE.0) THEN
                ! reads in the number of molecules present in the system
                IF (paral%io_parent)&
                     READ(iunit,err=99,END=99,fmt=*) dmbi%nmol
             ENDIF

             IF (INDEX(line,'NSAVESTATE').NE.0) THEN
                ! reads in the number of states to display and save at the end of the PT-INTERACTION calculation
                IF (paral%io_parent)&
                     READ(iunit,err=99,END=99,fmt=*) dmbi%max_num_disp
             ENDIF

             IF (INDEX(line,'SELFCONSISTENT').NE.0) THEN
                ! reads in the number of self-consistent iteration for the optimisation of W1
                IF (paral%io_parent)&
                     READ(iunit,err=99,END=99,fmt=*) dmbi%bptscfiter
             ENDIF
             IF (INDEX(line,'SCFPTTOL').NE.0) THEN
                ! reads in the convergence tolerance on the self-consistent iteration for th
                IF (paral%io_parent)&
                     READ(iunit,err=99,END=99,fmt=*) dmbr%scf_tol
             ENDIF

             IF (INDEX(line,'THETA').NE.0) THEN
                ! reads in the rotation angle used for rotations of Wannier functions
                ! Allocating memory
                ALLOCATE(btheta(dmbi%nmol),STAT=ierr)
                IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                     __LINE__,__FILE__)
                ! Reading in the different data
                DO i=1,dmbi%nmol
                   IF (paral%io_parent)&
                        READ(iunit,err=99,END=99,fmt=*) btheta(i)
                ENDDO
                ! Setting flag to initialised!
                dmbi%in_btheta=.TRUE.
             ENDIF

             IF (INDEX(line,'VEC').NE.0) THEN
                ! reads in the rotation axis used for rotations of Wannier functions
                ! Allocating memory
                ALLOCATE(BVECx(dmbi%nmol),STAT=ierr)
                IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                     __LINE__,__FILE__)
                ALLOCATE(BVECy(dmbi%nmol),STAT=ierr)
                IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                     __LINE__,__FILE__)
                ALLOCATE(BVECz(dmbi%nmol),STAT=ierr)
                IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                     __LINE__,__FILE__)
                ! Reading in the different data
                DO i=1,dmbi%nmol
                   IF (paral%io_parent)&
                        READ(iunit,err=99,END=99,fmt=*) BVECx(i),BVECy(i),&
                        BVECz(i)
                ENDDO
                ! Setting flag to initialised!
                dmbi%in_bvec=.TRUE.
             ENDIF

             IF (INDEX(line,'DX').NE.0) THEN
                ! reads in the displacement vector for the wannier function 
                ! Allocating memory
                ALLOCATE(BDXVECx(dmbi%nmol),STAT=ierr)
                IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                     __LINE__,__FILE__)
                ALLOCATE(BDXVECy(dmbi%nmol),STAT=ierr)
                IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                     __LINE__,__FILE__)
                ALLOCATE(BDXVECz(dmbi%nmol),STAT=ierr)
                IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                     __LINE__,__FILE__)
                ! Reading in the different data
                DO i=1,dmbi%nmol
                   IF (paral%io_parent)&
                        READ(iunit,err=99,END=99,fmt=*) BDXVECx(i),&
                        BDXVECy(i),BDXVECz(i)
                ENDDO
                ! Setting flag to initialised!
                dmbi%in_bdxvec=.TRUE.
             ENDIF

             IF (INDEX(line,'WLOAD').NE.0) THEN
                ! READING the wannier funtions FROM file
                dmbi%wann_load=.TRUE.
                dmbi%wann_save=.FALSE.
             ENDIF

             IF (INDEX(line,'WSAVE').NE.0) THEN
                ! SAVING the wannier funtions TO a file
                dmbi%wann_load=.FALSE.
                dmbi%wann_save=.TRUE.
             ENDIF

             IF (INDEX(line,'NOROT').NE.0) THEN
                ! rotation by zero degrees is made without projection
                dmbi%wann_allrot=.FALSE.
             ENDIF

             IF (INDEX(line,'WGAUSSIAN').NE.0) THEN
                ! Using a Gaussian representation for the input/output of Wannier functions 
                dmbi%wann_gauss=.TRUE.
             ENDIF

             IF (INDEX(line,'MULTIREF').NE.0) THEN
                ! Allowing a reference wannier file for each molecule
                dmbi%wann_multi=.TRUE.
             ENDIF

             IF (INDEX(line,'TRICK').NE.0) THEN
                ! Cutoff trick (lower cutoff for PT wavefunction)    
                dmbi%cutoff_trick=.TRUE.
             ENDIF

             IF (INDEX(line,'RESTRAINED').NE.0) THEN
                ! Cutoff is lowered for the perturbation calculation 
                dmbi%cutoff_restr=.TRUE.
             ENDIF

             ! Disables the localisation of the optimised wf 
             IF (INDEX(line,' NOLOC').NE.0) nmr_options%tlocalize = .FALSE.

             ! Enables the Wannier centre prediction routine (for BOdyn ONLY...)
             IF (INDEX(line,' WCPREDICT').NE.0) THEN
                dmbi%wcpredict=.TRUE.
                ! initialises the sample counter
                dmbi%wc_sample=1
                ! reads the number of cntl%md steps between localisation cycles 
                IF (paral%io_parent)&
                     READ(iunit,err=99,END=99,fmt=*) dmbi%wc_pred_step
             ENDIF

             ! Uses LDA to calculate the d2Exc/dn1 dn2 in PT-GGA calculations
             IF (INDEX(line,' FASTDEXC').NE.0) dmbi%fastdexc = .TRUE.

             ! Saves the orthogonalised W0 in a RESTART file and quits
             IF (INDEX(line,' SAVELOCW0').NE.0)&
                  dmbi%save_localised_w0 = .TRUE.

             ! Use the atomic wavefunction for W0                     
             IF (INDEX(line,' ATOMICWAVEFUNCTION').NE.0)&
                  dmbi%tatomicwavefunction = .TRUE.
             ! Use the atomic density with Wannier W0
             IF (INDEX(line,' ATOMICRHO').NE.0)&
                  dmbi%tatomicrho = .TRUE.
             ! Use non orthogonal Wannier combination
             IF (INDEX(line,' NONORTHOGWANNIER').NE.0)&
                  dmbi%torthog_wannier = .FALSE.
             ! Use simplified interaction model
             IF (INDEX(line,' SIMPLEMODEL').NE.0)&
                  dmbi%tsimple_model = .TRUE.
             ! Use OLD MANNO PT-INTER code (only works in serial)
             IF (INDEX(line,' MANNOPT').NE.0)&
                  dmbi%tmanno = .TRUE.
          ENDIF           ! INTERACTION-parameters

          ! ==--  GENERAL  -----------------------------------------------==
          ! ..   REALSPACE representation of the ground state orbitals:
          IF (INDEX(line,'REALSPACE').NE.0) THEN
             response1%tkeeprealspacewfn = .TRUE.
          ENDIF
          ! ..   Polak-Ribieres-Variant of the conjugate gradients:
          IF (INDEX(line,'POLAK').NE.0)  THEN
             response1%tpolak_p=.TRUE.
          ENDIF
          ! ..   Additional optimization control options:
          IF (INDEX(line,'NOOPT').NE.0)  THEN
             cntl%wfopt=.FALSE.
             restart1%restart=.TRUE.
             restart1%rwf=.TRUE.
          ENDIF
          ! Number of steps for which the CG-length is calculated analytically
          IF (INDEX(line,'CG-ANALYTIC').NE.0)  THEN
             IF (paral%io_parent)&
                  READ(iunit,err=99,END=99,fmt='(A80)') line
             CALL readsi(line ,1,iout,response1%cg_analytic,erread)
             IF (erread) GOTO 99
          ENDIF
          ! Factor the last analytic length is multiplied with,
          ! for istep > CG-ANALYTIC
          IF (INDEX(line,'CG-FACTOR').NE.0)  THEN
             IF (paral%io_parent)&
                  READ(iunit,err=99,END=99,fmt='(A80)') line
             CALL readsr(line ,1,iout,response2%cg_factor,erread)
             IF (erread) GOTO 99
          ENDIF
          ! Preconditioner type (divergence smoothing formula)
          IF (INDEX(line,'TIGHTPREC').NE.0) THEN
             response1%preconditioner_p=2
             ! State dependent preconditioner 
             ! DS variant
          ELSEIF (INDEX(line,'STATEPREC').NE.0) THEN
             response1%preconditioner_p=3
             ! AP variant
          ELSEIF (INDEX(line,'TIGHTSTATEPREC').NE.0) THEN
             response1%preconditioner_p=4
          ENDIF
          IF (INDEX(line,'NOMIN').NE.0)  response1%pcgmin_p=.FALSE.

          ! Used by no-cntl%pcg, for a first step with cntl%pcg
          IF (INDEX(line,'PCG_FIRST').NE.0)  response1%opt1st=.TRUE.


          ! ..   cntl%diis optimization (preconditioning active)
          IF ( INDEX(line,'DIIS').NE.0) THEN
             response1%pcg_p=.FALSE.
             response1%tsde_p=.FALSE.
             response1%diis_p=.TRUE.
          ENDIF
          IF (INDEX(line,'MDIIS').NE.0)THEN
             response1%pcg_p=.FALSE.
             response1%tsde_p=.FALSE.
             response1%diis_p=.TRUE.
             IF (paral%io_parent)&
                  READ(iunit,err=99,END=99,fmt='(A80)') line
             CALL readsi(line,1,ie,response1%mdiis_p,erread)
             IF (erread) GOTO 99
          ENDIF
          ! ..   Steepest descent:
          IF ( INDEX(line,'STEE').NE.0 .AND.&
               INDEX(line,'DESC').NE.0) THEN
             response1%pcg_p=.FALSE.
             response1%tsde_p=.TRUE.
             CALL stopgm('RESPIN_P',&
                  'WHO THE HELL is sill using steepest descent?',& 
                  __LINE__,__FILE__)
          ENDIF
          IF (INDEX(line,'PCG_P').NE.0) THEN
             ! ..   Precond. Conjugate Gradients
             response1%tsde_p=.FALSE.
             response1%pcg_p=.TRUE.
             response1%diis_p=.FALSE.
          ENDIF
          ! Convergence threashold (max gradient)
          IF (INDEX(line,'CONVERGE').NE.0) THEN
             IF (paral%io_parent)&
                  READ(iunit,err=99,END=99,fmt='(A80)') line
             CALL readsr(line ,1,iout,response2%tolog_p,erread)
             IF (erread) GOTO 99
          ENDIF
          IF ((INDEX(line,'HTHRS').NE.0) .OR.&
               (INDEX(line,'HAM')*INDEX(line,'CUT').NE. 0)) THEN
             ! ..   Threshold for Hessian cutoff in preconditioning matrix
             IF (paral%io_parent)&
                  READ(iunit,err=99,END=99,fmt='(A80)') line
             CALL readsr(line ,1,iout,response2%hthrs_p,erread)
             IF (erread) GOTO 99
          ENDIF
          IF (INDEX(line,'SIGMA').NE.0) THEN
             ! ..   The spread of the Gaussian in the virtual cell calculation
             IF (paral%io_parent)&
                  READ(iunit,err=99,END=99,fmt='(A80)') line
             CALL readsr(line ,1,iout,deltarvalues%sigma,erread)
             IF (erread) GOTO 99
          ENDIF
          ! ==------------------------------------------------------------==
       ENDDO                 ! while READ_MORE
       ! ==------------------------------------------------------------==
       ! ==------------------------------------------------------------==
       ! If we are loading a saved wavefunction:
       ! no optimisations and restart (use different file format..)
       IF ((dmbi%wann_load.OR.dmbi%tatomicwavefunction)&
            .AND. response1%tinteraction) THEN
          cntl%wfopt = .FALSE.
          restart1%restart = .FALSE.
          restart1%rwf = .FALSE.
       ENDIF
       IF (paral%parent.AND.dmbi%tmanno.AND.(parai%cp_nproc.GT.1)) THEN
          IF (paral%io_parent)&
               WRITE(6,*) 'Old MANNOPT code is not parallel...'
          IF (paral%io_parent)&
               WRITE(6,*) 'Switching to PARA-PT instead'
          dmbi%tmanno=.FALSE.
       ENDIF
       ! ==------------------------------------------------------------==
    ENDIF                     ! parent

    CALL mp_bcast_byte(response1,size_in_bytes_of(response1),&
         parai%io_source,parai%cp_grp)
    CALL mp_bcast_byte(response2, size_in_bytes_of(response2),parai%io_source,parai%cp_grp)

    CALL mp_bcast_byte(eig1,size_in_bytes_of(eig1),parai%io_source,parai%cp_grp)
    CALL mp_bcast_byte(eig2,size_in_bytes_of(eig2),parai%io_source,parai%cp_grp)
    CALL mp_bcast_byte(nmr_options,size_in_bytes_of(nmr_options), &
         parai%io_source,parai%cp_grp)
    CALL mp_bcast_byte(epr_options,size_in_bytes_of(epr_options),&
         parai%io_source,parai%cp_grp)
    CALL mp_bcast(nmr_para%nmr_threads,parai%io_source,parai%cp_grp)
    CALL mp_bcast_byte(voa_options,size_in_bytes_of(voa_options), &
         parai%io_source,parai%cp_grp)

    CALL mp_bcast_byte(deltarvalues,size_in_bytes_of(deltarvalues),parai%io_source,parai%cp_grp)
    CALL mp_bcast_byte(ownpotvalue,size_in_bytes_of(ownpotvalue),parai%io_source,parai%cp_grp)
    CALL mp_bcast_byte(kpertpar,size_in_bytes_of(kpertpar),parai%io_source,parai%cp_grp)
    CALL mp_bcast_byte(dmbi,size_in_bytes_of(dmbi),parai%io_source,parai%cp_grp)
    CALL mp_bcast_byte(dmbr,size_in_bytes_of(dmbr),parai%io_source,parai%cp_grp)

    CALL mp_bcast(cntl%wfopt,parai%io_source,parai%cp_grp)
    CALL mp_bcast_byte(restart1,size_in_bytes_of(restart1),parai%io_source,parai%cp_grp)
    CALL mp_bcast_byte(lancphon,size_in_bytes_of(lancphon),parai%io_source,parai%cp_grp)
    CALL mp_bcast(numf,parai%source,parai%cp_grp)
    CALL mp_bcast(tweight,parai%source,parai%cp_grp)
    IF (numf.NE.0) THEN
       IF (.NOT.paral%parent)  THEN
          ALLOCATE(nf(numf),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
       IF (.NOT.paral%parent)  THEN
          ALLOCATE(wghtf(numf+1),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
       CALL mp_bcast(nf,SIZE(nf),parai%source,parai%cp_grp)
       CALL mp_bcast(wghtf,SIZE(wghtf),parai%source,parai%cp_grp)
    ENDIF
    ! ==------------------------------------------------------------==
    IF (.NOT. paral%parent) THEN
       IF (response1%tinteraction) THEN
          IF (dmbi%in_btheta) THEN
             ALLOCATE(btheta(dmbi%nmol),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             CALL mp_bcast(btheta,SIZE(btheta),parai%io_source,parai%cp_grp)
          ENDIF
          IF (dmbi%in_bvec) THEN
             ALLOCATE(BVECx(dmbi%nmol),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             ALLOCATE(BVECy(dmbi%nmol),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             ALLOCATE(BVECz(dmbi%nmol),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             CALL mp_bcast(BVECx,SIZE(BVECx),parai%io_source,parai%cp_grp)
             CALL mp_bcast(BVECy,SIZE(BVECy),parai%io_source,parai%cp_grp)
             CALL mp_bcast(BVECz,SIZE(BVECz),parai%io_source,parai%cp_grp)
          ENDIF
          IF (dmbi%in_bdxvec) THEN
             ALLOCATE(BDXVECx(dmbi%nmol),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             ALLOCATE(BDXVECy(dmbi%nmol),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             ALLOCATE(BDXVECz(dmbi%nmol),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             CALL mp_bcast(BDXVECx,SIZE(BDXVECx),parai%io_source,parai%cp_grp)
             CALL mp_bcast(BDXVECy,SIZE(BDXVECy),parai%io_source,parai%cp_grp)
             CALL mp_bcast(BDXVECz,SIZE(BDXVECz),parai%io_source,parai%cp_grp)
          ENDIF
       ENDIF
    ENDIF
    ! ==------------------------------------------------------------==
    ! ==------------------------------------------------------------==
    ! Are all options consistent? [Only one response at a time...!]
    iresponse=0
    IF (response1%tnmr) iresponse=iresponse+1
    IF (response1%tepr) iresponse=iresponse+1
    IF (response1%traman) iresponse=iresponse+1
    IF (response1%teigensystem) iresponse=iresponse+1
    IF (response1%phonon) iresponse=iresponse+1
    IF (response1%thardness) iresponse=iresponse+1
    IF (response1%tfukui) iresponse=iresponse+1
    IF (response1%tinteraction) iresponse=iresponse+1
    IF (response1%tlanphon) iresponse=iresponse+1
    IF (response1%tkpert) iresponse=iresponse+1
    IF (response1%tdummyatom) iresponse=iresponse+1
    IF (response1%tvoa) iresponse=iresponse+1

    IF ((.NOT.cntl%tinr).AND.(iresponse.NE.1)) THEN
       IF (paral%parent) THEN
          IF (paral%io_parent)&
               WRITE (6,*) 'PLEASE CHOOSE EXACTLY ONE OF THE FOLLOWING ',&
               '&RESP OPTIONS:'
          IF (paral%io_parent)&
               WRITE (6,*) '       NMR'
          IF (paral%io_parent)&
               WRITE (6,*) '       EPR'
          IF (paral%io_parent)&
               WRITE (6,*) '       RAMAN'
          IF (paral%io_parent)&
               WRITE (6,*) '       PHONON'
          IF (paral%io_parent)&
               WRITE (6,*) '       EIGENSYSTEM'
          IF (paral%io_parent)&
               WRITE (6,*) '       HARDNESS'
          IF (paral%io_parent)&
               WRITE (6,*) '       FUKUI'
          IF (paral%io_parent)&
               WRITE (6,*) '       INTERACTION'
          IF (paral%io_parent)&
               WRITE (6,*) '       LANCZOS'
          IF (paral%io_parent)&
               WRITE (6,*) '       KPERT'
          IF (paral%io_parent)&
               WRITE (6,*) '       OAPC'
          IF (paral%io_parent)&
               WRITE (6,*) '       VOA'
       ENDIF
       CALL stopgm('RESPIN','REQUESTED FLAGS INCONSISTENT'&
            // 'WITH LINEAR RESPONSE CODE',& 
            __LINE__,__FILE__)
    ENDIF

    IF (response1%tkpert)THEN
       tkpts%tsymkp=.FALSE.
       tkpts%tkfull=.FALSE.
       tkpts%tkall=.FALSE.
       tkpts%tkblock= .FALSE.
       nkpt%nkpnt=1
       nkpt%nblkp=1
    ENDIF
    IF (response1%tkpert .AND. tkpts%tkpnt) THEN
       CALL stopgm&
            ('  RESPIN',&
            'STANDARD K-POINT AND KPERT ARE NOT CONSISTENT',& 
            __LINE__,__FILE__)
    ENDIF


    RETURN
    ! ==--------------------------------------------------------------==
99  CONTINUE
    CALL stopgm('RESPIN',' ERROR IN READING INPUT FILE',& 
         __LINE__,__FILE__)
  END SUBROUTINE respin_p
  ! ==--------------------------------------------------------------==

END MODULE respin_p_utils
