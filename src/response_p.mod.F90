MODULE response_pmod
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==================================================================
  ! == global response parameters & variables                       ==
  ! ==================================================================
  ! t_initial_guess_c1: if true, then use the input c1 as initial
  ! guess. otherwise, take c1=0.
  ! restart_c1: read c1 set
  INTEGER,PARAMETER :: lresponse=39
  TYPE :: response1_t
     LOGICAL :: response_running
     LOGICAL :: tdummyatom
     LOGICAL :: tdummyatom_ref
     LOGICAL :: tdummy_step
     LOGICAL :: tdumpbd
     LOGICAL :: tdumn
     LOGICAL :: tdumsd
     LOGICAL :: tdummyatom_tec
     LOGICAL :: toacp4f
     LOGICAL :: traman
     LOGICAL :: phonon
     LOGICAL :: teigensystem
     LOGICAL :: tnmr
     LOGICAL :: thardness
     LOGICAL :: tkpert
     LOGICAL :: tfukui
     LOGICAL :: tvoa
     LOGICAL :: tsde_p
     LOGICAL :: prec_p
     LOGICAL :: diis_p
     LOGICAL :: pcg_p
     LOGICAL :: pcgmin_p
     LOGICAL :: opt1st
     LOGICAL :: tpolak_p
     LOGICAL :: pr_energy
     LOGICAL :: tlanphon
     LOGICAL :: tlanph_cont
     LOGICAL :: teigens
     LOGICAL :: tinteraction
     LOGICAL :: projout
     LOGICAL :: rotout
     LOGICAL :: trho1
     LOGICAL :: tnonlocal
     LOGICAL :: tkeeprealspacewfn
     LOGICAL :: t_restart_c1
     LOGICAL :: t_initialguess_c1
     INTEGER :: cg_analytic
     INTEGER :: mdiis_p
     INTEGER :: preconditioner_p
     LOGICAL :: tepr
  END TYPE response1_t
  TYPE(response1_t), SAVE :: response1

  INTEGER,PARAMETER :: rresponse=15

  TYPE :: response2_t
     REAL(real_8) :: hthrs_p
     REAL(real_8) :: tolog_p
     REAL(real_8) :: cg_factor
     REAL(real_8) :: dumstep1
     REAL(real_8) :: dumstep2
     REAL(real_8) :: dumstep3
     REAL(real_8) :: dumstep4
     REAL(real_8) :: dumstep5
     REAL(real_8) :: dumcrit
     INTEGER :: tdumnum
     INTEGER :: tdumqm
     REAL(real_8) :: rad_dum
     REAL(real_8) :: rad_norm
     INTEGER :: tfdimen
     INTEGER :: ang_mom
  END TYPE response2_t
  TYPE(response2_t), SAVE :: response2

  INTEGER,PARAMETER :: response_write=314,response_read=315
  ! these two are constant definitions for the interface 
  ! to the restart routines
  ! ==================================================================
  ! == energy and charge quantities                                 ==
  ! ==================================================================
  TYPE :: ener1_t
     REAL(real_8) :: eloc1=0.0_real_8
     REAL(real_8) :: enl1=0.0_real_8
     REAL(real_8) :: eht1=0.0_real_8
     REAL(real_8) :: exc1=0.0_real_8
     REAL(real_8) :: elag1=0.0_real_8
     REAL(real_8) :: eh0=0.0_real_8
  END TYPE ener1_t
  TYPE(ener1_t), SAVE :: ener1
  ! ==--------------------------------------------------------------==
  ! where:
  ! eloc1:   <1|  h(1)^local |0> + cc  =  eivps(r) n1(r)
  ! enl1:    <1| h(1)^direct |0> + cc  =  2 * <1| h1nl>
  ! which is calculated via the (d)fnl(00) for phonons
  ! enmr:    <1| h(1)^direct |0> + cc  for option nmr
  ! eht1:    1/2  n1(r') n1(r)  /  |r-r'|  +  v1_pp[gaussians](r) n1(r)
  ! where this exists only for phonons: ^^^^^^^^^^^^^^^^^^^^^^^^^
  ! exc1:    1/2  d^2 e_xc / dn^2  n1(r) n1(r')
  ! elag1:   - <1|  eps_ks  |1>  (the lagrange multipliers)
  ! eh0:    <1| h0 |1>
  ! eraman:  ?
  ! ==================================================================
  ! ==  data derived from the ground state wavefunctions that is
  ! ==  used during the perturbation calculation.
  ! ==================================================================
  REAL(real_8), ALLOCATABLE, SAVE :: vofrho0(:,:)

  REAL(real_8), POINTER, SAVE :: rho0(:,:) ! rho0(nnr1,*)
  REAL(real_8), POINTER, SAVE :: vxc0(:,:) ! vxc0(nnr1,*)

  COMPLEX(real_8), ALLOCATABLE, SAVE :: c0real(:,:)
  ! ==================================================================



  ! ==================================================================
  ! interaction parameters and variables
  ! ==================================================================
  INTEGER,PARAMETER :: linteraction=44

  REAL(real_8), ALLOCATABLE, SAVE :: wc_array1(:,:)
  REAL(real_8), ALLOCATABLE, SAVE :: wc_array2(:,:)
  REAL(real_8), ALLOCATABLE, SAVE :: wc_velocity(:,:)

  REAL(real_8), ALLOCATABLE, SAVE :: btheta(:)
  REAL(real_8), ALLOCATABLE, SAVE :: bvecx(:)
  REAL(real_8), ALLOCATABLE, SAVE :: bvecy(:)
  REAL(real_8), ALLOCATABLE, SAVE :: bvecz(:)

  REAL(real_8), ALLOCATABLE, SAVE :: bdxvecx(:)
  REAL(real_8), ALLOCATABLE, SAVE :: bdxvecy(:)
  REAL(real_8), ALLOCATABLE, SAVE :: bdxvecz(:)

  REAL(real_8), ALLOCATABLE, SAVE :: emono(:)

  INTEGER, ALLOCATABLE, SAVE :: statemin(:)
  INTEGER, ALLOCATABLE, SAVE :: statemax(:)


  TYPE :: dmbi_t
     INTEGER :: interaction_first
     INTEGER :: wc_sample
     INTEGER :: blmax
     INTEGER :: nmol
     INTEGER :: ngw_zero
     INTEGER :: bptscfiter
     LOGICAL :: wann_load
     LOGICAL :: wann_save
     LOGICAL :: wann_allrot
     LOGICAL :: wann_multi
     LOGICAL :: in_btheta
     LOGICAL :: in_bvec
     LOGICAL :: in_bdxvec
     LOGICAL :: cutoff_trick
     LOGICAL :: cutoff_restr
     LOGICAL :: tinter_ortho
     LOGICAL :: wann_gauss
     LOGICAL :: tlinearscaling
     LOGICAL :: inter_pt_firstcall
     LOGICAL :: wcpredict
     LOGICAL :: fastdexc
     LOGICAL :: tnoresponse
     LOGICAL :: save_localised_w0
     LOGICAL :: tatomicwavefunction
     LOGICAL :: tmanno
     LOGICAL :: tatomicrho
     LOGICAL :: torthog_wannier
     LOGICAL :: tsimple_model
     INTEGER :: wc_pred_step
     INTEGER :: max_num_disp
  END TYPE dmbi_t
  TYPE(dmbi_t), SAVE :: dmbi

  TYPE :: dmbr_t
     REAL(real_8) :: scf_tol
  END TYPE dmbr_t
  TYPE(dmbr_t), SAVE :: dmbr

  ! ==================================================================
  ! simple_model (option of interaction)
  ! parameters and variables
  ! ==================================================================
  REAL(real_8), ALLOCATABLE, SAVE :: s_star(:,:)

  ! ==================================================================
  ! ==  phonon                                                      ==
  ! ==================================================================

  REAL(real_8), ALLOCATABLE, SAVE :: fnl00(:,:,:,:,:)
  REAL(real_8), ALLOCATABLE, SAVE :: dfnl00(:,:,:,:,:,:)

  ! ==================================================================
  ! ==  lanczos/phonon                                              ==
  ! ==================================================================
  TYPE :: lancphon_t
     LOGICAL :: details
     INTEGER :: nlan_st
  END TYPE lancphon_t
  TYPE(lancphon_t), SAVE :: lancphon

  ! ==================================================================
  ! == kpert                                                        ==
  ! ================================================================== 
  ! rrk  storage array for k points
  ! wwk  storage array for k weights
  ! restart_cu1 logical variable for writing the cu1 on resta_p
  ! tk_prep_diag logical variable for building the c00 ={c0}+{cu1}
  ! tkread_c1 logical variable for reading the cu1 from resta_p
  ! tpert_hamilt for calculation of the h matrix dim 2*nstate
  ! dtwnl_dk(1:3,1:2*ngw,1:ngh(is),1:nsp)  first derivatives 
  ! of the non local projectors with respect to kx ky and kz,
  ! for each g-component (kleiman-bylander form)
  ! ddtwnl_ddk(1:2*ngw,1:ngh(is),1:nsp,1:3,1:3) second derivatives
  ! of the non local projectors with respect to kxkx kxky 
  ! kxkz,kykx, kyky, kykz, kzkx, kzky, kzkz
  ! twnl_p, twnl_m
  ! dfnl_dk,ddfnl_ddk
  ! pp projectors in g+ki and g-ki (where ki is considered in kx 
  ! or ky or kz direction)
  ! twnl_p  (1:2*ngw,1:ngh(is),1:nsp,1:3)  |
  ! twnl_m  (1:2*ngw,1:ngh(is),1:nsp,1:3)  | initialized in putwnl_kpert
  ! ylmb_p  (1:2*ngh,lpmax,1:3)            | form factors coefficients in
  ! ylmb_m  (1:2*ngh,lpmax,1:3)            | g+ki  and g-ki
  ! e2_nl_c0 (3,3) correction to the perturbed energy due to the second
  ! order terms related to the nlpp and which 
  ! do not depend on the  response wfn

  TYPE :: kpertpar_t
     LOGICAL :: norestart
     LOGICAL :: restart_cu1
     LOGICAL :: tk_prep_diag
     LOGICAL :: tkread_c1
     LOGICAL :: tpert_hamilt
  END TYPE kpertpar_t
  TYPE(kpertpar_t), SAVE :: kpertpar

  REAL(real_8), ALLOCATABLE, SAVE :: rrk(:,:)
  REAL(real_8), ALLOCATABLE, SAVE :: wwk(:)


  REAL(real_8), ALLOCATABLE, SAVE :: dtwnl_dk(:,:,:,:)
  REAL(real_8), ALLOCATABLE, SAVE :: ddtwnl_ddk(:,:,:,:,:)


  REAL(real_8), ALLOCATABLE, SAVE :: twnl_p(:,:,:,:)
  REAL(real_8), ALLOCATABLE, SAVE :: twnl_m(:,:,:,:)


  REAL(real_8), ALLOCATABLE, SAVE :: dfnl_dk(:,:,:,:,:)
  REAL(real_8), ALLOCATABLE, SAVE :: ddfnl_ddk(:,:,:,:,:,:)

  REAL(real_8), ALLOCATABLE, SAVE :: e2_nl_c0(:,:)

  ! ==================================================================
  ! ==   dynamic allocation of  arrays for kpert                    ==
  ! ==                                                              ==
  ! ==                                                              ==
  ! ==================================================================
  ! ==   pkpert(3,nstate,nstate)   : <c0|grad|c0>                   ==
  ! ==   ekpert(2*nstate,nkpts)    : eigenvalues                    ==
  ! ==   fkpert(2*nstate,nkpts)    : occupations                    ==
  ! ==   hamilk(2*nstate,2*nstate) : matrix elements depending on k ==
  ! ==================================================================
  REAL(real_8), ALLOCATABLE, SAVE :: ekpert(:,:)
  REAL(real_8), ALLOCATABLE, SAVE :: fkpert(:,:)

  COMPLEX(real_8), ALLOCATABLE, SAVE :: hamilk(:,:)

  ! 
  ! ==================================================================
  ! ==     hgkp_pk(nhg,nkpts) hgkm_pk(nhg.nkpts)
  ! ==     tk_pk(2*ngw,nhxs,nsx)


  ! ===================================================================
  ! == matrix not depending on k vector
  COMPLEX(real_8), ALLOCATABLE, SAVE :: h1_00(:,:,:)

  COMPLEX(real_8), ALLOCATABLE, SAVE :: h2_00(:,:,:,:)

  COMPLEX(real_8), ALLOCATABLE, SAVE :: h1_10(:,:,:,:)

  REAL(real_8), ALLOCATABLE, SAVE :: h0_11(:,:,:,:)





  ! ==================================================================
  ! == selected eigenvalues                                         ==
  ! ==================================================================
  INTEGER, SAVE :: neigens,iter_max_eigens,ikrylov_max
  INTEGER, PARAMETER :: neigens_max=20
  REAL(real_8), SAVE :: tol_eigens
  REAL(real_8),PARAMETER :: tol_direct=1.e-3_real_8
  TYPE :: eigens1_t
     INTEGER :: neigens
     INTEGER :: iter_max_eigens
     REAL(real_8) :: tol_eigens
     INTEGER :: ikrylov_max
  END TYPE eigens1_t
  TYPE(eigens1_t), SAVE :: eigens1

  ! ==================================================================
  ! == eigensystem                                                  ==
  ! ==================================================================
  TYPE :: eig1_t
     INTEGER :: lanczos_dim
     INTEGER :: num_lanczos_states
  END TYPE eig1_t
  TYPE(eig1_t), SAVE :: eig1
  TYPE :: eig2_t
     REAL(real_8) :: conv_threshold
     REAL(real_8) :: eigenvalue_shift
  END TYPE eig2_t
  TYPE(eig2_t), SAVE :: eig2



  ! ==================================================================
  ! ==  nmr                                                         ==
  ! ==================================================================
  ! == inmr_method:  0  determine default                           ==
  ! ==               1  csgt                                        ==
  ! == inmr_virtual  0  determine default                           ==
  ! ==               3  llc where int[rho_i] = min                  ==
  ! ==               4  center = wannier centers from localize      ==
  ! ==               5  no virtual cells, llc = (1,1,1)             ==
  ! ==                                                              ==
  ! ==                                                              ==
  ! ==                                                              ==
  ! == tnmr_full     if the full calculation is desired.            ==
  ! ==                                                              ==
  ! ==                                                              ==
  TYPE :: nmr_options_t
     LOGICAL :: t_nmr_first
     LOGICAL :: tcurrent
     LOGICAL :: tverbose
     LOGICAL :: tprintwfns
     LOGICAL :: tprintrho
     LOGICAL :: tlocalize
     LOGICAL :: tfast
     LOGICAL :: tsmoothing
     INTEGER :: inmr_method
     INTEGER :: inmr_virtual
     LOGICAL :: tnmr_full
     LOGICAL :: tnmr_overlaps
     LOGICAL :: tforce_xyz(3)
     LOGICAL :: tnmr_only_error
  END TYPE nmr_options_t
  TYPE(nmr_options_t), SAVE :: nmr_options

  INTEGER, PARAMETER ::inmr_default=0,inmr_iglo=1,inmr_csgt=2,&
       inmr_llc=3,inmr_wc=4,inmr_novirt=5
  ! ==================================================================
  ! ==  nmr superparallel variables                                 ==
  ! ==================================================================
  TYPE :: nmr_para_t
     INTEGER :: nmr_supergroup
     INTEGER :: nmr_threads
     INTEGER :: nmr_mygroup
     INTEGER :: nmr_total_nproc
     INTEGER :: parents(0:7)
     INTEGER :: nmr_supersource
     LOGICAL :: nmr_superparent
  END TYPE nmr_para_t
  TYPE(nmr_para_t), SAVE :: nmr_para

  ! ==================================================================
  ! ==  internal nmr variables             -------------------------==
  ! ==  used only by the nmr_*_p routines  -------------------------==
  REAL(real_8), ALLOCATABLE, SAVE :: wanniercenters(:,:)
  REAL(real_8), SAVE :: firstshelldistance

  INTEGER, ALLOCATABLE, SAVE :: lower_left(:,:)
  INTEGER, ALLOCATABLE, SAVE :: lower_left_value(:,:)
  INTEGER, SAVE :: llc_ind(3),llc_val(3),nris(3)

  TYPE :: origin_t
     INTEGER :: llc_ind(3)
     INTEGER :: llc_val(3)
     INTEGER :: nris(3)
  END TYPE origin_t
  TYPE(origin_t), SAVE :: origin
  TYPE :: deltarvalues_t
     REAL(real_8) :: sigma
     REAL(real_8) :: overlap_threashold
  END TYPE deltarvalues_t
  TYPE(deltarvalues_t), SAVE :: deltarvalues
  ! the lower_left(3,*) contains the optimal starting coordinates for
  ! the cell for the wavefunctions, i.e. the box on whose borders
  ! the wannier function psi_i vanishes. 
  ! for the calculation of these borders,
  ! a smearing gaussian of spread "sigma" is applied.
  ! ==--------------------------------------------------------------==
  REAL(real_8), ALLOCATABLE, SAVE :: shift_matrix(:,:,:,:)
  REAL(real_8), ALLOCATABLE, SAVE :: chi(:,:)
  REAL(real_8), SAVE :: shift_factor
  ! factor is alpha^2 * ppm
  REAL(real_8), SAVE :: chi_factor,chi_si_to_shift_ppm,chi_si_to_ppmcgs,&
       chi_si_iso=0.0_real_8

  TYPE :: shifts_t
     REAL(real_8) :: shift_factor
     REAL(real_8) :: chi_factor
     REAL(real_8) :: chi_si_to_shift_ppm
     REAL(real_8) :: chi_si_to_ppmcgs
     REAL(real_8) :: chi_si_iso
  END TYPE shifts_t
  TYPE(shifts_t), SAVE :: shifts
  CHARACTER(len=44), SAVE :: timetag

  ! ==================================================================
  ! ==  epr                                                         ==
  ! ==================================================================
  ! == iepr_method:  0  determine default                           ==
  ! ==               1  csgt                                        ==
  ! == iepr_virtual  0  determine default                           ==
  ! ==               3  llc where int[rho_i] = min                  ==
  ! ==               4  center = wannier centers from localize      ==
  ! ==               5  no virtual cells, llc = (1,1,1)             ==
  ! ==                                                              ==
  ! ==                                                              ==
  ! ==                                                              ==
  ! == tepr_full     if the full calculation is desired.            ==
  ! ==                                                              ==
  ! ==                                                              ==

  TYPE :: epr_options_t
     LOGICAL :: t_epr_first
     LOGICAL :: teverbose
     LOGICAL :: teprintwfns
     LOGICAL :: teprintrho
     LOGICAL :: telocalize
     LOGICAL :: tefast
     LOGICAL :: tesmoothing
     INTEGER :: iepr_method
     INTEGER :: iepr_virtual
     LOGICAL :: tepr_full
     LOGICAL :: tepr_overlaps
     LOGICAL :: teforce_xyz(3)
     LOGICAL :: tepr_only_error
     LOGICAL :: tepr_smart
     LOGICAL :: tepr_hyp
     LOGICAL :: tepr_ownpot
  END TYPE epr_options_t
  TYPE(epr_options_t), SAVE :: epr_options

  INTEGER, PARAMETER :: iepr_default=0,iepr_iglo=1,iepr_csgt=2,&
       iepr_llc=3,iepr_wc=4,iepr_novirt=5


  ! ==================================================================
  ! ==  internal epr variables             -------------------------==
  ! ==  used only by the epr_*_p routines  -------------------------==
  TYPE :: ownpotvalue_t
     REAL(real_8) :: epr_ownpot_value
     REAL(real_8) :: epr_smart_value
  END TYPE ownpotvalue_t
  TYPE(ownpotvalue_t), SAVE :: ownpotvalue

  ! ==================================================================
  ! ==  fukui calculation variables                                 ==
  ! ==================================================================
  LOGICAL, SAVE :: tweight
  INTEGER, ALLOCATABLE, SAVE :: nf(:)
  INTEGER, SAVE :: numf
  REAL(real_8), ALLOCATABLE, SAVE :: wghtf(:)


  ! ==================================================================
  ! ==  voa calculation variables                                   ==
  ! ==================================================================
  TYPE :: voa_options_t
     LOGICAL :: tat
     LOGICAL :: tmd
     LOGICAL :: tpf
     LOGICAL :: tcurrent
     LOGICAL :: tdensity
     LOGICAL :: tdebug
     LOGICAL :: tverbose
     LOGICAL :: tatomlist
     LOGICAL :: tamat
     LOGICAL :: tlag
     CHARACTER(len=1000) :: atomlistline
  END TYPE voa_options_t
  TYPE(voa_options_t), SAVE :: voa_options

  TYPE :: voa_data_t
     LOGICAL :: timagpert
     LOGICAL, ALLOCATABLE :: atomlist(:)
     CHARACTER(len=6) :: fmt = 'E18.10'
     REAL(real_8) :: box_center(3)
     REAL(real_8) :: p_gs(3)
     REAL(real_8) :: w_eps = 1.e-4_real_8
     REAL(real_8), ALLOCATABLE :: rotmat(:,:)
     REAL(real_8), ALLOCATABLE :: apt_pf(:,:,:)
     REAL(real_8), ALLOCATABLE :: apt_vf(:,:,:)
     REAL(real_8), ALLOCATABLE :: aat_cog(:,:,:)
     REAL(real_8), ALLOCATABLE :: aat_dog(:,:,:)
     REAL(real_8), ALLOCATABLE :: amat_nl(:,:,:)
     REAL(real_8), ALLOCATABLE :: p_vf(:,:)
     REAL(real_8), ALLOCATABLE :: m_vf(:,:)
     REAL(real_8), ALLOCATABLE :: r_wc(:,:)
     COMPLEX(real_8), ALLOCATABLE :: rc0(:,:)
     COMPLEX(real_8), ALLOCATABLE :: rc1(:,:)
     COMPLEX(real_8), ALLOCATABLE :: h1psi0(:,:)
     COMPLEX(real_8), ALLOCATABLE :: vnlc0(:,:)
     COMPLEX(real_8), ALLOCATABLE :: pic0(:,:,:)
     COMPLEX(real_8), ALLOCATABLE :: ric0(:,:,:)
     COMPLEX(real_8), ALLOCATABLE :: vnlric0(:,:,:)
     COMPLEX(real_8), ALLOCATABLE :: rxpic0(:,:,:)
     COMPLEX(real_8), ALLOCATABLE :: c1src(:,:,:,:)
  END TYPE voa_data_t
  TYPE(voa_data_t), SAVE :: voa_data

END MODULE response_pmod
