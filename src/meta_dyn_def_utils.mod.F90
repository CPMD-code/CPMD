MODULE meta_dyn_def_utils
  USE cnst_dyn,                        ONLY: &
       bbeta, cv_associated, cv_langamma, cv_langevintemp, cv_temp, cv_temp0, &
       ekincv, fmtdres, imeta, inter_hill, inter_hill_max, l2dc, lcvtc, &
       lfullc, lisoc, lkfix, lmeta, ltcglobal, lupd, maxrmsdat, mdcellr, &
       ncolvar, nlong, nshort, nsubsys, optdist, rch, rcw, rmeta, &
       tcvlangevin, tcvscale, toll_avcv, tvolbound, vharm
  USE kinds,                           ONLY: real_8
  USE mw,                              ONLY: mwi,&
                                             tmw
  USE system,                          ONLY: cntr

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: meta_dyn_def
  !public :: s

CONTAINS

  SUBROUTINE meta_dyn_def
    ! ==================================================================
    ! ==      SETS DEFAULT FOR VARIABLES OF META DYNAMICS             ==
    ! ==================================================================



    lmeta%lcnstdyn     = .FALSE.
    lmeta%lcolvardyn   = .FALSE.
    lmeta%lhills       = .TRUE.
    lmeta%hlore        = .FALSE.
    lmeta%hshift       = .FALSE.
    lmeta%hratio       = .FALSE.
    lmeta%sphere       = .FALSE.
    lmeta%meta_restart = .FALSE.
    lmeta%ltune_cscl   = .FALSE.
    lmeta%ltune_hh     = .FALSE.
    lmeta%thillvol     = .FALSE.
    lmeta%tresfile     = .FALSE.
    lmeta%qmmmorder    = .TRUE.

    imeta%i_meta_max   = 100
    imeta%i_meta_res   = 0
    imeta%st_freq      = 50
    imeta%tr_freq      = 1
    imeta%qw_freq      = 10000

    fmtdres      = 'MTD_RESTART'

    imeta%icheck       = 20
    rmeta%hllw         = 0.1_real_8
    rmeta%hllh         = 0.001_real_8
    rmeta%htop         = 0.016_real_8
    rmeta%hlow         = 0.0001_real_8
    rmeta%expup        = 2.0_real_8
    rmeta%expdown      = 16._real_8
    rmeta%rshift       = 2.0_real_8
    rmeta%fboost       = 1.0_real_8

    rmeta%tolkin       = -1.0_real_8

    ! Constraints only
    imeta%i_cnst_min   = 50
    imeta%i_cnst_max   = 400
    rmeta%toll_fcnst   = 3.e-3_real_8
    imeta%nstep_relax  = 6
    lmeta%doublequench = .TRUE.

    ! Collective Variables only
    inter_hill_max = 300
    ncolvar     = 0
    inter_hill  = 100
    toll_avcv   = 0.001_real_8

    ! Extended Lagrangean
    lmeta%lextlagrange=.FALSE.
    lcvtc=.FALSE.
    ltcglobal=.FALSE.
    tcvscale=.FALSE.
    tcvlangevin=.FALSE.
    lkfix=.TRUE.
    cv_associated =.FALSE.
    cv_temp0=cntr%tempw
    cv_temp =0.0_real_8
    cv_langevintemp=0.0_real_8
    cv_langamma=0.0_real_8
    ekincv  =0.0_real_8
    vharm   =0.0_real_8

    ! Hydrogen Bond Chain
    Rcw     = 6.4_real_8
    Rch     = 3.0_real_8
    optdist = 5.5_real_8
    bbeta   = 100._real_8
    nshort  = 1
    nlong   = 5
    lupd    = 40

    ! RMDS_SEQ
    maxrmsdat = 0

    ! Cell Metadynamics
    lmeta%lmeta_cell = .FALSE.
    lmeta%lmc_ext    = .FALSE.
    lfullc     = .TRUE.
    lisoc      = .FALSE.
    l2dc       = .FALSE.
    tvolbound  = .FALSE.
    mdcellr%vbarrier = 0.1_real_8

    ! Multi-Metadynamics
    nsubsys = 1
    lmeta%tmulti  = .FALSE.

    ! Spin Localization
    lmeta%tlocalizespin=.FALSE.

    ! Analysis
    lmeta%tcvanalysis = .FALSE.

    ! use trajectory interval from METASTORE (.true.) or &CPMD (.false.)
    lmeta%ttrajovrr = .TRUE.

    ! Monitor CV along cntl%md trajectory
    lmeta%tcvmonitor = .FALSE.
    imeta%wcv_freq   = 50

    ! Metadynamics with CV as HT functions
    lmeta%tcvcell = .FALSE.

    ! Multiple walker metadynamics
    tmw=.FALSE.
    mwi%nwalk=1
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE meta_dyn_def
  ! ==================================================================

END MODULE meta_dyn_def_utils
