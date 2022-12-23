! ==================================================================
! == ELECTRONIC PROPAGATION AND SPECTRA                           ==
! ==================================================================
MODULE td_input
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE
  REAL(real_8) :: gpot,gpotv(3),genl,gampl(3)
  INTEGER :: curr_step,gndir
  INTEGER :: itermax
  LOGICAL ::tcurr
  LOGICAL :: ttransg,ppcorr,tberry
  INTEGER :: niterations

  COMPLEX(real_8), ALLOCATABLE :: cz(:)
  COMPLEX(real_8), ALLOCATABLE :: cz1(:)
  COMPLEX(real_8), ALLOCATABLE :: cz2(:)
  COMPLEX(real_8), ALLOCATABLE :: czm(:)
  COMPLEX(real_8), ALLOCATABLE :: czp(:)
  COMPLEX(real_8), ALLOCATABLE :: bei1(:,:)
  COMPLEX(real_8), ALLOCATABLE :: bei2(:,:)
  COMPLEX(real_8), ALLOCATABLE :: bei3(:,:)

  REAL(real_8), ALLOCATABLE :: rpos(:)
  REAL(real_8), ALLOCATABLE :: cpos(:)
  REAL(real_8), ALLOCATABLE :: spos(:)
  REAL(real_8), ALLOCATABLE :: f_bohmtraj(:,:,:)
  REAL(real_8), ALLOCATABLE :: bpot(:,:,:)

  TYPE :: gaugefield_t
     REAL(real_8) :: gswtt
     REAL(real_8) :: gdfreq
     REAL(real_8) :: muij
     INTEGER :: gimax
     INTEGER :: nperiod
     LOGICAL :: pi_pulse
     LOGICAL :: tkgauge
     INTEGER :: kgauge
  END TYPE gaugefield_t
  TYPE(gaugefield_t), SAVE, PUBLIC :: gaugefield
  TYPE :: pointcharge_t
     REAL(real_8) :: pcht0
     REAL(real_8) :: pchdt
     REAL(real_8) :: pchint
     REAL(real_8) :: pchfx
     REAL(real_8) :: pchfy
     REAL(real_8) :: pchfz
     INTEGER :: pchwr
  END TYPE pointcharge_t
  TYPE(pointcharge_t), SAVE, PUBLIC :: pointcharge
  TYPE :: td_prop_t
     REAL(real_8) :: tintrvll
     REAL(real_8) :: tt
     REAL(real_8) :: tt2
     REAL(real_8) :: ttime
     REAL(real_8) :: epsil
     REAL(real_8) :: ampl
     REAL(real_8) :: tdfreq
     INTEGER :: nzeros
     INTEGER :: n_cycles
     INTEGER :: pertdir
     INTEGER :: pert_type
     INTEGER :: read_wf
     INTEGER :: ionized_state
     LOGICAL :: stextpot
     LOGICAL :: td_extpot
     LOGICAL :: tpointch
     LOGICAL :: ionize_from_state
     LOGICAL :: ionize_from_rdm
     LOGICAL :: td_fix_spin_dens
     LOGICAL :: do_not_normalize
  END TYPE td_prop_t
  TYPE(td_prop_t), SAVE, PUBLIC :: td_prop
  TYPE :: masklog_t
     LOGICAL :: tmask    
  END TYPE masklog_t
  TYPE(masklog_t), SAVE, PUBLIC :: masklog
  TYPE :: maskreal_t
     REAL(real_8) :: maskpar1
     REAL(real_8) :: maskpar2
     REAL(real_8) :: maskpar3
  END TYPE maskreal_t
  TYPE(maskreal_t), SAVE, PUBLIC :: maskreal
  TYPE :: bohm_t
     REAL(real_8) :: gw_par1
     REAL(real_8) :: gw_par2
     REAL(real_8) :: gw_p1
     REAL(real_8) :: gw_p2
     REAL(real_8) :: gw_p3
     REAL(real_8) :: gw_r1
     REAL(real_8) :: gw_r2
     REAL(real_8) :: gw_r3
     INTEGER :: nf_bohmtraj
     INTEGER :: eh_restart_n
     LOGICAL :: tgaussian_wp
     LOGICAL :: tplanewave_wp
     LOGICAL :: tbohmtraj
     LOGICAL :: restart_bohm
     LOGICAL :: eh_restart_dens
  END TYPE bohm_t
  TYPE(bohm_t), PUBLIC, SAVE :: bohm
END MODULE td_input
! ==================================================================
