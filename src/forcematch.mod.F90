! =======================================================================
! include file for forcematching
! 
! =======================================================================
MODULE forcematch
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! logicals
  ! copy of FMATCH in system.h for files which do not include system.h
  LOGICAL :: fm_fmatch
  ! compute a whole reference trajectory (i.e. do MD), not yet
  ! implemented
  LOGICAL :: fm_compute_traj
  ! compute forces, fields, etc. via single points from a reference
  ! TRAJECTORY
  LOGICAL :: fm_compute_sp
  ! read in forces from external file
  LOGICAL :: fm_read_forces
  ! read in covalent forces from external file
  LOGICAL :: fm_read_covforces
  ! equivalencies
  LOGICAL :: fm_equivalence
  ! fit charges
  LOGICAL :: fm_fit_charges
  ! fit charges only
  LOGICAL :: fm_charges_only
  ! fit force constants only (not equilibrium bond lengths/angles)
  LOGICAL :: fm_fit_fc_only
  ! fit parameters for bonds,angles,dihedrals,and improper dihedrals
  LOGICAL :: fm_fit_bonds,fm_fit_angles,fm_fit_diheds,fm_fit_improps
  ! re-initialize the wavefunction for each step
  LOGICAL :: fm_rinit_wf
  ! restart of the single point calculations
  ! restart from the last frame index written to the FM_REF files
  LOGICAL :: fm_sp_restart
  ! compute RMS of total forces
  LOGICAL :: fm_compute_rms

  TYPE :: fm_l_t
     LOGICAL :: fm_fmatch
     LOGICAL :: fm_compute_traj
     LOGICAL :: fm_compute_sp
     LOGICAL :: fm_read_forces
     LOGICAL :: fm_read_covforces
     LOGICAL :: fm_equivalence
     LOGICAL :: fm_fit_charges
     LOGICAL :: fm_charges_only
     LOGICAL :: fm_fit_fc_only
     LOGICAL :: fm_fit_bonds
     LOGICAL :: fm_fit_angles
     LOGICAL :: fm_fit_diheds
     LOGICAL :: fm_fit_improps
     LOGICAL :: fm_rinit_wf
     LOGICAL :: fm_sp_restart
     LOGICAL :: fm_compute_rms
  END TYPE fm_l_t
  TYPE(fm_l_t) :: fm_l

  ! strings
  ! filename with reference forces, etc.
  CHARACTER(len=80) :: fm_reference_filename
  ! filenames for qfit
  ! Manuel: I renamed these variables
  ! character*80      fm_qm_ref_filename, fm_nn_ref_filename
  ! filename for the reference Hirschfeld charges      
  CHARACTER(len=80) :: fm_chj_ref_filename
  ! Manuel: filename for the reference pip file
  CHARACTER(len=80) :: fm_pip_ref_filename
  ! Manuel: filename for the reference forces
  CHARACTER(len=80) :: fm_forces_ref_filename
  ! Manuel: reference trajectory file name from which 
  ! geos for single points are extracted
  CHARACTER(len=80) :: fm_ref_traj_filename
  ! file where covalent QM forces are written to (to fit force
  ! constants, etc.)
  CHARACTER(len=80) :: fm_covforces_ref_filename
  ! name of optimized topology file
  CHARACTER(len=80) :: fm_outtopo_filename

  TYPE :: fm_s_t
     CHARACTER(len=80) :: fm_reference_filename
     CHARACTER(len=80) :: fm_chj_ref_filename
     CHARACTER(len=80) :: fm_ref_traj_filename
     CHARACTER(len=80) :: fm_outtopo_filename
     CHARACTER(len=80) :: fm_pip_ref_filename
     CHARACTER(len=80) :: fm_forces_ref_filename
     CHARACTER(len=80) :: fm_covforces_ref_filename
  END TYPE fm_s_t
  TYPE(fm_s_t) :: fm_s

  ! integers
  ! total nr. of frames (that are used for fitting)
  INTEGER :: fm_nframes
  ! total nr. of QM atoms (that we fit on..)
  INTEGER :: fm_nrqm
  ! stride to apply to reference trajectory when extracting
  ! the geos for the single points
  INTEGER :: fm_ref_traj_stride
  ! keep the nr. of frame (that is in TRAJECTORY_REF!) available for
  ! consistent output
  INTEGER :: fm_iframe
  ! equivalencies
  ! NOTE: fm_equiv is read from input file and contains 
  ! GROMOS atom numbers!
  INTEGER, PARAMETER :: fm_maxeq = 1000 
  INTEGER :: fm_nreq
  INTEGER :: fm_equiv(2, fm_maxeq)

  INTEGER, POINTER :: grm_equiv(:)
  INTEGER, POINTER :: qm_equiv(:)

  ! index  when extracting a frame from the TRAJECTORY_REF
  ! file. this will be printed in the FM_REF_PIP, FM_REF_CHJ 
  ! and the FM_REF_FORCES files.
  INTEGER :: fm_frame_indx
  ! maximum number of allowed fixed charges
  INTEGER, PARAMETER :: fm_max_q_fixers = 1000 
  INTEGER, PARAMETER :: fm_max_qw=1000 
  ! cpmd atom indexes for which charges should stay fixed
  INTEGER :: fm_fixed_q_indx(fm_max_q_fixers)
  INTEGER :: fm_fixed_q_grm_indx(fm_max_q_fixers)
  ! actual number of fixed charges
  INTEGER :: fm_n_fix_q
  ! number of weights for the charge fitting given by user
  INTEGER :: fm_n_qw
  ! index array for the weights on charges specified by the user      
  INTEGER :: fm_wq_indx(fm_max_qw)
  INTEGER :: fm_wq_grm_indx(fm_max_qw)


  TYPE :: fm_i_t
     INTEGER :: fm_nframes
     INTEGER :: fm_nrqm
     INTEGER :: fm_ref_traj_stride
     INTEGER :: fm_nreq
     INTEGER :: fm_equiv(2, fm_maxeq)
     INTEGER, POINTER :: grm_equiv(:)
     INTEGER, POINTER :: qm_equiv(:)
     INTEGER :: fm_frame_indx
     INTEGER :: fm_fixed_q_indx(fm_max_q_fixers)
     INTEGER :: fm_n_fix_q
     INTEGER :: fm_n_qw
     INTEGER :: fm_wq_indx(fm_max_qw)
  END TYPE fm_i_t
  TYPE(fm_i_t) :: fm_i
  ! =======================================================================
  ! specific for charge fitting
  ! logicals
  ! Find optimal weights
  LOGICAL :: fm_optimize_weights
  ! MANUEL
  ! fix some of the charges
  LOGICAL :: fm_fix_q

  TYPE :: fm_qfit_l_t
     LOGICAL :: fm_optimize_weights
     LOGICAL :: fm_fix_q
  END TYPE fm_qfit_l_t
  TYPE(fm_qfit_l_t) :: fm_qfit_l

  ! integers
  ! MM atoms excluded from charge-fitting
  INTEGER, PARAMETER :: fm_max_nexclude = 100 
  INTEGER :: fm_nexclude
  INTEGER :: fm_exclude(fm_max_nexclude)


  ! reals
  ! weight for potential
  REAL(real_8) :: fm_wv
  ! weight for field
  REAL(real_8) :: fm_wf
  ! MANUEL (array) weights for individual charge restraints
  REAL(real_8) :: fm_wq(fm_max_qw)
  REAL(real_8) :: fm_wq_unique(fm_max_qw)
  ! weight for the atoms for which no individual weight is
  ! specified in the input
  REAL(real_8) :: fm_wq_general
  ! weight for total charge
  REAL(real_8) :: fm_wtot
  ! Manuel target value for the charges to be fixed
  REAL(real_8) :: fm_fixed_q_trgt(fm_max_q_fixers)

  TYPE :: fm_qfit_r_t
     REAL(real_8) :: fm_wv
     REAL(real_8) :: fm_wf
     REAL(real_8) :: fm_wq(fm_max_qw)
     REAL(real_8) :: fm_wq_unique(fm_max_qw)
     REAL(real_8) :: fm_wtot
     REAL(real_8) :: fm_fixed_q_trgt(fm_max_q_fixers)
  END TYPE fm_qfit_r_t
  TYPE(fm_qfit_r_t) :: fm_qfit_r
  ! =======================================================================
  ! topology related stuff and stuff for fitting covalent forces
  ! logicals
  ! amber or gromos functional form? (copy of qmmm_amber in
  ! mm_input.inc)
  LOGICAL :: fm_amber

  TYPE :: fm_topo_l_t
     LOGICAL :: fm_amber
  END TYPE fm_topo_l_t
  TYPE(fm_topo_l_t) :: fm_topo_l

  ! integers
  ! save exclusion lists
  INTEGER, POINTER :: fm_save_INE(:)
  INTEGER, POINTER :: fm_save_KNE(:)
  INTEGER, POINTER :: fm_save_JSNE(:)
  INTEGER, POINTER :: fm_save_INE14(:)
  INTEGER, POINTER :: fm_save_KNE14(:)
  INTEGER, POINTER :: fm_save_JSNE14(:)
  ! bond, etc. types to be optimized. fm_nbonopt is the number of bond
  ! types to optimize. fm_bonopt(1:fm_nbonopt) contains the numbers of 
  ! optimized bond types. Analogous for ang,imp,dih.
  INTEGER :: fm_nbonopt, fm_nangopt
  INTEGER :: fm_nimpopt, fm_ndihopt
  INTEGER, POINTER :: fm_bonopt(:), fm_angopt(:)
  INTEGER, POINTER :: fm_impopt(:), fm_dihopt(:)

  ! number of individual bonds, etc. that are optimized
  INTEGER :: fm_nbon, fm_nang, fm_nimp, fm_ndih

  ! individual bonds, etc. that belong to types that are optimized
  ! fm_ib(:, :) will be dimensioned fm_ib(0:2, MAX_NR_OF_BONDS), and
  ! for bond ib holds the type in fm_ib(0, ib), and the CPMD numbers
  ! of the atoms forming the bond in fm_ib(1,ib) and fm_ib(2,ib)
  ! analogous for the others, e.g. for dihedrals type is in
  ! fm_ip(0,ip), and the 4 atoms in fm_ip(1:4,ip)
  ! _ib, _it, _iq, _ip hold bonds, angles, impropers, dihedrals, resp.
  INTEGER, POINTER :: fm_ib(:, :), fm_it(:, :)
  INTEGER, POINTER :: fm_iq(:, :), fm_ip(:, :)

  ! max nr. of iterations for Levenberg-Marquardt
  INTEGER :: fm_max_kfit_iterations

  TYPE :: fm_topo_i_t
     INTEGER, POINTER :: fm_save_INE(:)
     INTEGER, POINTER :: fm_save_KNE(:)
     INTEGER, POINTER :: fm_save_JSNE(:)
     INTEGER, POINTER :: fm_save_INE14(:)
     INTEGER, POINTER :: fm_save_KNE14(:)
     INTEGER, POINTER :: fm_save_JSNE14(:)
     INTEGER :: fm_nbonopt
     INTEGER :: fm_nangopt
     INTEGER :: fm_nimpopt
     INTEGER :: fm_ndihopt
     INTEGER, POINTER :: fm_bonopt(:)
     INTEGER, POINTER :: fm_angopt(:)
     INTEGER, POINTER :: fm_impopt(:)
     INTEGER, POINTER :: fm_dihopt(:)
     INTEGER :: fm_nbon
     INTEGER :: fm_nang
     INTEGER :: fm_nimp
     INTEGER :: fm_ndih
     INTEGER, POINTER :: fm_ib(:, :)
     INTEGER, POINTER :: fm_it(:, :)
     INTEGER, POINTER :: fm_iq(:, :)
     INTEGER, POINTER :: fm_ip(:, :)
     INTEGER :: fm_max_kfit_iterations
  END TYPE fm_topo_i_t
  TYPE(fm_topo_i_t) :: fm_topo_i

  ! reals
  ! coordinates of QM atoms (dimensions are 3, nrqm, nframes)
  REAL(real_8), POINTER :: fm_rqm(:, :, :)
  ! reference forces
  REAL(real_8), POINTER :: fm_fqm_ref(:, :, :)
  ! store force constants of QM bonds, etc. (they are set to 0 for
  ! QM/MM, but we want to keep the original as an initial guess for
  ! fitting
  REAL(real_8), POINTER :: fm_kbon(:), fm_kang(:)
  REAL(real_8), POINTER :: fm_kimp(:), fm_kdih(:)

  ! force constants and equilibrium values for fitting (we use our own
  ! arrays going from 1..NR_OF_OPTIMIZED_TYPES
  REAL(real_8), POINTER :: fm_kb(:), fm_b0(:)
  REAL(real_8), POINTER :: fm_kt(:), fm_t0(:)
  REAL(real_8), POINTER :: fm_kq(:), fm_q0(:)
  REAL(real_8), POINTER :: fm_kp(:), fm_p0(:), fm_pn(:)

  TYPE :: fm_topo_r_t
     REAL(real_8), POINTER :: fm_rqm(:, :, :)
     REAL(real_8), POINTER :: fm_fqm_ref(:, :, :)
     REAL(real_8), POINTER :: fm_kbon(:)
     REAL(real_8), POINTER :: fm_kang(:)
     REAL(real_8), POINTER :: fm_kimp(:)
     REAL(real_8), POINTER :: fm_kdih(:)
     REAL(real_8), POINTER :: fm_kb(:)
     REAL(real_8), POINTER :: fm_b0(:)
     REAL(real_8), POINTER :: fm_kt(:)
     REAL(real_8), POINTER :: fm_t0(:)
     REAL(real_8), POINTER :: fm_kq(:)
     REAL(real_8), POINTER :: fm_q0(:)
     REAL(real_8), POINTER :: fm_kp(:)
     REAL(real_8), POINTER :: fm_p0(:)
     REAL(real_8), POINTER :: fm_pn(:)
  END TYPE fm_topo_r_t
  TYPE(fm_topo_r_t) :: fm_topo_r
  ! =======================================================================
  ! capping
  ! logicals
  ! are there capping hydrogens?
  LOGICAL :: fm_capping
  LOGICAL, POINTER :: fm_iscap(:)

  TYPE :: fm_cap_l_t
     LOGICAL, POINTER :: fm_iscap(:)
     LOGICAL :: fm_capping
  END TYPE fm_cap_l_t
  TYPE(fm_cap_l_t) :: fm_cap_l
  ! integers
  ! max. nr of capping hydrogens (this should be equal to 
  ! the parameter npax which is defined in MM_Interface/mm_cap_H.F
  INTEGER, PARAMETER :: fm_maxcap = 200 
  ! actual nr of capping hydrogens
  INTEGER :: fm_ncap
  ! In array fm_cap, atom fm_cap(1, i) is the nr of the hydrogen atom
  ! that caps atom fm_cap(2, i), initially with GROMOS numbers
  INTEGER :: fm_cap(2, fm_maxcap)

  TYPE :: fm_cap_i_t
     INTEGER :: fm_ncap
     INTEGER :: fm_cap(2, fm_maxcap)
  END TYPE fm_cap_i_t
  TYPE(fm_cap_i_t) :: fm_cap_i

END MODULE forcematch

! =======================================================================

