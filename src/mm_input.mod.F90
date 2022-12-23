MODULE mm_input
  USE kinds,                           ONLY: real_8
! CP species and gromos atom numbers of the QM atoms (read in mm_detsp)

  IMPLICIT NONE

  TYPE :: gqmmm_t
     INTEGER, ALLOCATABLE :: gr_atom(:,:)
     INTEGER, ALLOCATABLE :: gr_nasp(:)
     LOGICAL, ALLOCATABLE :: is_dummy(:)
  END TYPE gqmmm_t
  TYPE(gqmmm_t), SAVE :: gqmmm
  TYPE :: igqmmm_t
     INTEGER :: sp_H_added
  END TYPE igqmmm_t
  TYPE(igqmmm_t), SAVE :: igqmmm
  ! QMMM block in standard input (read in mm_read_qmmm_input)
  ! Gromos configuration files
  TYPE :: cqmmm_t
     CHARACTER(len=120) :: fileconf
     CHARACTER(len=120) :: filegrin
     CHARACTER(len=120) :: filetopo
  END TYPE cqmmm_t
  TYPE(cqmmm_t), SAVE :: cqmmm
  TYPE :: lqmmm_t
     LOGICAL :: qmmm
     LOGICAL :: qmmm_time
     LOGICAL :: qmmm_verbose
     LOGICAL :: qmmm_reflex
     LOGICAL :: qmmm_dcd
  END TYPE lqmmm_t
  TYPE(lqmmm_t), SAVE :: lqmmm
  TYPE :: rqmmm_t
     REAL(real_8) :: Rcut_el
     REAL(real_8) :: Rcut_small
     REAL(real_8) :: Rcut_esp
     REAL(real_8) :: box_toll
     REAL(real_8) :: box_wall
  END TYPE rqmmm_t
  TYPE(rqmmm_t), SAVE :: rqmmm
  TYPE :: iqmmm_t
     INTEGER :: coupl_model ! electrostatic coupling model
     INTEGER :: nt_nlist  ! the electrostatic coupling nn list
     ! is updated every nt_nlist steps
     INTEGER :: nt_sample ! sampling of trajectory of the interacting atoms
  END TYPE iqmmm_t
  TYPE(iqmmm_t), SAVE :: iqmmm
  ! carbons to be hydrogenized
  INTEGER, PARAMETER :: N_carbon_max=100
  TYPE :: addh_t
     INTEGER :: igdummy
     INTEGER :: N_carbon
     INTEGER :: ig_carbon(N_carbon_max)
     INTEGER :: n_added_H
  END TYPE addh_t
  TYPE(addh_t), SAVE :: addh
  ! Cap carbons
  TYPE :: capp_t
     LOGICAL :: cap_H
  END TYPE capp_t
  TYPE(capp_t), SAVE :: capp
  ! Fully classical MD
  TYPE :: clc_t
     LOGICAL :: classical
  END TYPE clc_t
  TYPE(clc_t), SAVE :: clc
  ! Flexible solvent (water)
  INTEGER, PARAMETER :: nwq_max=100
  TYPE :: solqmmm_t
     LOGICAL :: tflexsol
     LOGICAL :: all_water
  END TYPE solqmmm_t
  TYPE(solqmmm_t), SAVE :: solqmmm
  TYPE :: solqmmi_t
     INTEGER :: solv_bond(20)
     INTEGER :: nwq
     INTEGER :: ig_wq(nwq_max)
  END TYPE solqmmi_t
  TYPE(solqmmi_t), SAVE :: solqmmi
  ! smearing function for electrostatic interaction
  TYPE :: n_sml_t
     LOGICAL :: no_hirshfeld
  END TYPE n_sml_t
  TYPE(n_sml_t), SAVE :: n_sml
  TYPE :: n_sm_t
     INTEGER :: n_smear
  END TYPE n_sm_t
  TYPE(n_sm_t), SAVE :: n_sm
  TYPE :: r_esp_t
     REAL(real_8) :: esp_weight
  END TYPE r_esp_t
  TYPE(r_esp_t), SAVE :: r_esp
  ! max number of interacting atoms. MAXNN in input.
  TYPE :: mne_t
     INTEGER :: MAXNAT_el
  END TYPE mne_t
  TYPE(mne_t), SAVE :: mne
  ! write potential
  TYPE :: wp_l_t
     LOGICAL :: write_potential
     LOGICAL :: write_density
     LOGICAL :: write_orbs
  END TYPE wp_l_t
  TYPE(wp_l_t), SAVE :: wp_l
  TYPE :: wp_i_t
     INTEGER :: n_stride
     INTEGER :: NFI_wp
     INTEGER :: NFI_wd
     INTEGER :: NFI_lt
  END TYPE wp_i_t
  TYPE(wp_i_t), SAVE :: wp_i
  TYPE :: wp_c_t
     CHARACTER(len=30) :: cubename_dens
     CHARACTER(len=30) :: cubename_pot
  END TYPE wp_c_t
  TYPE(wp_c_t), SAVE :: wp_c
  ! write orbitals
  TYPE :: wp_p_t
     INTEGER :: n_print_wf
     INTEGER :: i_wf_pr(999)
  END TYPE wp_p_t
  TYPE(wp_p_t), SAVE :: wp_p
  ! amber force field
  TYPE :: agr_t
     LOGICAL :: qmmm_amber
  END TYPE agr_t
  TYPE(agr_t), SAVE :: agr
  TYPE :: ragr_t
     REAL(real_8) :: scale_14_amber
  END TYPE ragr_t
  TYPE(ragr_t), SAVE :: ragr
  ! split stuff
  TYPE :: sppph_t
     LOGICAL :: mm_split
  END TYPE sppph_t
  TYPE(sppph_t), SAVE :: sppph
  ! charge group restrain
  INTEGER, PARAMETER :: max_N_CG=100
  TYPE :: cgrest_i_t
     INTEGER :: n_cg
     INTEGER :: atom_qm_cg(max_N_CG)
  END TYPE cgrest_i_t
  TYPE(cgrest_i_t), SAVE :: cgrest_i
  TYPE :: cgrest_r_t
     REAL(real_8) :: q_rest
     REAL(real_8) :: lambda
  END TYPE cgrest_r_t
  TYPE(cgrest_r_t), SAVE :: cgrest_r
  ! restart from trajectory
  TYPE :: rtr_l_t
     LOGICAL :: restart_traj
  END TYPE rtr_l_t
  TYPE(rtr_l_t), SAVE :: rtr_l
  TYPE :: rtr_i_t
     INTEGER :: iframe_restart
  END TYPE rtr_i_t
  TYPE(rtr_i_t), SAVE :: rtr_i
  TYPE :: rtr_c_t
     CHARACTER(len=40) :: file_traj_restart
  END TYPE rtr_c_t
  TYPE(rtr_c_t), SAVE :: rtr_c
  TYPE :: rtr_r_t
     REAL(real_8) :: rev_vel
  END TYPE rtr_r_t
  TYPE(rtr_r_t), SAVE :: rtr_r
  ! EXCLUSION STUFF
  INTEGER, PARAMETER :: NCe_max=10000
  TYPE :: excl_comm_t
     INTEGER :: atom_qm_excl(NCe_max)
     INTEGER :: atom_mm_excl(NCe_max)
     INTEGER :: NCe
  END TYPE excl_comm_t
  TYPE(excl_comm_t), SAVE :: excl_comm
  TYPE :: excl_comm_l_t
     LOGICAL :: gromos_excl
     LOGICAL :: excl_mech
  END TYPE excl_comm_l_t
  TYPE(excl_comm_l_t), SAVE :: excl_comm_l

  ! Write out EQM and QM/MM related stuff
  TYPE :: eqm_r_t
     REAL(real_8) :: eqm=0.0_real_8
     REAL(real_8) :: eqmmm=0.0_real_8
     REAL(real_8) :: eqmmm_cl=0.0_real_8
     REAL(real_8) :: emm=0.0_real_8
     REAL(real_8) :: eqmmm0=0.0_real_8
     REAL(real_8) :: eext0=0.0_real_8
     REAL(real_8) :: eexcl=0.0_real_8
  END TYPE eqm_r_t
  TYPE(eqm_r_t), SAVE :: eqm_r

  ! GROMOS VELOCITY BLOCK
  TYPE :: g96_vel_t
     INTEGER :: NTX_vel
  END TYPE g96_vel_t
  TYPE(g96_vel_t), SAVE :: g96_vel
END MODULE mm_input
