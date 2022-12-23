MODULE cnst_dyn
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! =================================================================
  ! == INCLUDE FILE OF CONSTRAINTS DYNAMICS                         ==
  ! ==================================================================
  ! ==--------------------------------------------------------------==
  REAL(real_8), ALLOCATABLE, SAVE :: cscl_fac(:,:)
  REAL(real_8), ALLOCATABLE, TARGET, SAVE :: cnst_val(:,:)
  REAL(real_8), ALLOCATABLE, SAVE :: cscl_val(:,:)

  REAL(real_8), ALLOCATABLE, SAVE :: hllw_val(:,:)
  REAL(real_8), ALLOCATABLE, SAVE :: hllh_val(:,:)

  REAL(real_8), ALLOCATABLE, SAVE :: vbound(:,:)

  INTEGER, ALLOCATABLE, SAVE :: ibound(:)

  LOGICAL, ALLOCATABLE, SAVE :: tad_scf(:)

  CHARACTER(len=20), SAVE :: fmtdres
  ! ==--------------------------------------------------------------==
  ! superseded by AU_KCM in cnst.inc.  AK:2008/03/21
  ! real(8) :: KcFAC
  ! PARAMETER :: KcFAC = 627.5132632235_real_8 
  ! ==================================================================
  ! ==  COLLECTIVE VARIABLES DYNAMICS                               ==
  ! ==================================================================
  INTEGER, PARAMETER :: maxlcv=3 
  ! ==--------------------------------------------------------------==
  INTEGER, SAVE :: ncolvar,inter_hill,inter_hill_max
  INTEGER, ALLOCATABLE, SAVE :: tycvar(:)
  INTEGER, ALLOCATABLE, SAVE :: atcvar(:,:)
  INTEGER, ALLOCATABLE, SAVE :: iangcv(:)


  INTEGER, ALLOCATABLE, SAVE :: iatdlmn(:)
  INTEGER, ALLOCATABLE, SAVE :: specindex(:)


  REAL(real_8), ALLOCATABLE, SAVE :: iqatdlmn(:)
  REAL(real_8), ALLOCATABLE, SAVE :: fhills(:,:,:)
  REAL(real_8), ALLOCATABLE, SAVE :: dof_incv(:,:,:)
  REAL(real_8), ALLOCATABLE, SAVE :: cv_ist(:)
  REAL(real_8), ALLOCATABLE, SAVE :: det_colvar(:,:)
  REAL(real_8), ALLOCATABLE, SAVE :: cvpar(:,:)
  REAL(real_8), ALLOCATABLE, TARGET, SAVE :: cv_path(:,:)

  ! ==--------------------------------------------------------------==
  REAL(real_8), SAVE :: toll_avcv
  ! ==--------------------------------------------------------------==
  LOGICAL, SAVE :: allocall
  ! ==================================================================
  ! ==  EXTENDED LAGRANGIAN DYNAMICS                                ==
  ! ==================================================================
  REAL(real_8), ALLOCATABLE, SAVE :: cv_dyn(:)
  REAL(real_8), ALLOCATABLE, SAVE :: cv_vel(:)
  REAL(real_8), ALLOCATABLE, SAVE :: cv_dyn_0(:)
  REAL(real_8), ALLOCATABLE, SAVE :: cv_mass(:)
  REAL(real_8), ALLOCATABLE, SAVE :: kharm(:)
  REAL(real_8), ALLOCATABLE, SAVE :: ra(:)
  REAL(real_8), ALLOCATABLE, SAVE :: ekincv_walk(:)
  REAL(real_8), ALLOCATABLE, SAVE :: vharm_walk(:)

  LOGICAL, ALLOCATABLE, SAVE :: initial_value(:)

  LOGICAL, ALLOCATABLE, SAVE :: skiphill(:)



  LOGICAL, SAVE :: lcvtc,&
       lchekharm,&
       lkfix,&
       cv_associated,&
       ltcglobal,&
       tcvscale,&
       tcvlangevin
  ! ==--------------------------------------------------------------==
  REAL(real_8), SAVE :: cv_temp0,&
       cv_temp,&
       cv_dtemp,&
       ekincv=HUGE(0.0_real_8),&
       vharm=HUGE(0.0_real_8),&
       cv_langevintemp,&
       cv_langamma
  ! ==================================================================
  ! ==  FROM TRANSITION STATE TO MINIMA                             ==
  ! ==================================================================
  INTEGER, PARAMETER :: max_pippo=10 
  ! ==--------------------------------------------------------------==
  INTEGER, SAVE :: ncvmin,imincheck,max_minchk,max_search,itol_type
  ! ==--------------------------------------------------------------==
  REAL(real_8), ALLOCATABLE, SAVE :: cv_min_tol(:)
  REAL(real_8), ALLOCATABLE, SAVE :: cv_min2(:,:)
  REAL(real_8), ALLOCATABLE, SAVE :: cv_tol_ext(:,:)

  ! ==================================================================
  ! ==  HYDROGEN BOND CHAIN                                         ==
  ! ==================================================================
  INTEGER, SAVE :: nshort ,nlong,lupd
  REAL(real_8), SAVE :: rcw,rch,bbeta,optdist
  ! ==================================================================
  ! ==  META DYNAMICS OF THE CELL                                   ==
  ! ==================================================================
  REAL(real_8), ALLOCATABLE, SAVE :: ht_ist(:)
  REAL(real_8), ALLOCATABLE, SAVE :: det_ht(:,:)
  REAL(real_8), ALLOCATABLE, TARGET, SAVE :: ht_path(:,:)

  CHARACTER(len=10), SAVE :: ht_name(6)
  ! ==--------------------------------------------------------------==
  LOGICAL, SAVE :: lfullc,lisoc,l2dc,tvolbound
  ! ==================================================================
  ! ==  MULTI -    META DYNAMICS                                    ==
  ! ==================================================================
  INTEGER, SAVE :: nsubsys
  ! ==--------------------------------------------------------------==
  INTEGER, ALLOCATABLE, SAVE :: ncvsys(:)

  LOGICAL, ALLOCATABLE, SAVE :: lha_sys(:)

  REAL(real_8), ALLOCATABLE, SAVE :: hwm(:)
  REAL(real_8), ALLOCATABLE, SAVE :: hhm(:)
  REAL(real_8), ALLOCATABLE, SAVE :: hlhm(:)
  REAL(real_8), ALLOCATABLE, SAVE :: hthm(:)
  REAL(real_8), ALLOCATABLE, SAVE :: hvm0(:)


  ! ==================================================================
  ! ==   SPIN DENSITY LOCALIZATION                                  ==
  ! ==================================================================
  REAL(real_8), SAVE :: rcc(3),en_harm_locs
  ! ==--------------------------------------------------------------==
  REAL(real_8), ALLOCATABLE, SAVE :: rcc0(:,:)

  INTEGER, ALLOCATABLE, SAVE :: icv_spin(:)

  ! ==================================================================
  ! ==  META DYNAMICS WITH VARIABLE CELL                            ==
  ! ==================================================================
  REAL(real_8), SAVE :: dstrmeta(6)=0.0_real_8,&
       fvbound(6)=0.0_real_8,&
       dstrcv(6)=0.0_real_8
  ! ==--------------------------------------------------------------==
  INTEGER, ALLOCATABLE, SAVE :: icv_cell(:)

  REAL(real_8), ALLOCATABLE, SAVE :: det_celvar(:,:,:)

  ! ==================================================================
  ! ==  RMSD_AB                                                     ==
  ! ==================================================================
  LOGICAL, SAVE :: trmsd_ab
  INTEGER, SAVE :: nrmsd_ab
  ! ==--------------------------------------------------------------==
  INTEGER, ALLOCATABLE, SAVE :: icv_rmsd_ab(:)
  INTEGER, SAVE :: maxrmsdat
  INTEGER, PARAMETER :: max_nnvar=200 
  CHARACTER(len=50), SAVE :: file_str_ab(max_nnvar)
  ! ==================================================================
  ! ==   COORGROUP                                                  ==
  ! ==================================================================
  INTEGER, PARAMETER :: max_natcnga=10 
  INTEGER, PARAMETER :: max_natcngb=50 
  INTEGER, ALLOCATABLE, SAVE :: iatcnga(:)
  INTEGER, ALLOCATABLE, SAVE :: iatcngb(:)
  INTEGER, ALLOCATABLE, SAVE :: natcngb(:)

  REAL(real_8), ALLOCATABLE, SAVE :: rccnga(:)

  ! ==================================================================
  TYPE :: imeta_t
     INTEGER :: i_meta_max
     INTEGER :: i_meta_res
     INTEGER :: i_cnst_min
     INTEGER :: i_cnst_max
     INTEGER :: icheck
     INTEGER :: nstep_relax
     INTEGER :: wcv_freq
     INTEGER :: st_freq
     INTEGER :: tr_freq
     INTEGER :: qw_freq
  END TYPE imeta_t
  TYPE(imeta_t), SAVE :: imeta
  TYPE :: mdcellr_t
     REAL(real_8) :: fh_cell(3,3)
     REAL(real_8) :: volmin
     REAL(real_8) :: volmax
     REAL(real_8) :: vbarrier
  END TYPE mdcellr_t
  TYPE(mdcellr_t), SAVE :: mdcellr
  TYPE :: rmeta_t
     REAL(real_8) :: hllw=0.0_real_8
     REAL(real_8) :: hllh=0.0_real_8
     REAL(real_8) :: toll_fcnst=0.0_real_8
     REAL(real_8) :: htop=0.0_real_8
     REAL(real_8) :: hlow=0.0_real_8
     REAL(real_8) :: gausspot=0.0_real_8
     REAL(real_8) :: eham_hill=0.0_real_8
     REAL(real_8) :: tolkin=0.0_real_8
     REAL(real_8) :: expup=0.0_real_8
     REAL(real_8) :: expdown=0.0_real_8
     REAL(real_8) :: rshift=0.0_real_8
     REAL(real_8) :: fboost=0.0_real_8
     REAL(real_8) :: hvol0=0.0_real_8
  END TYPE rmeta_t
  TYPE(rmeta_t), SAVE :: rmeta
  TYPE :: lmeta_t
     LOGICAL :: lcnstdyn
     LOGICAL :: lcolvardyn
     LOGICAL :: meta_restart
     LOGICAL :: lhills
     LOGICAL :: hlore
     LOGICAL :: hratio
     LOGICAL :: hshift
     LOGICAL :: sphere
     LOGICAL :: ltune_cscl
     LOGICAL :: ltune_hh
     LOGICAL :: doublequench
     LOGICAL :: lextlagrange
     LOGICAL :: lsadpnt
     LOGICAL :: lmdreinit
     LOGICAL :: lmeta_cell
     LOGICAL :: lmc_ext
     LOGICAL :: tmulti
     LOGICAL :: tcvanalysis
     LOGICAL :: tcvmonitor
     LOGICAL :: ttrajovrr
     LOGICAL :: thillvol
     LOGICAL :: tlocalizespin
     LOGICAL :: tcvcell
     LOGICAL :: tresfile
     LOGICAL :: qmmmorder
     LOGICAL :: skiphill_mw=.FALSE.
     LOGICAL :: randwalk
  END TYPE lmeta_t
  TYPE(lmeta_t), SAVE :: lmeta

END MODULE cnst_dyn
