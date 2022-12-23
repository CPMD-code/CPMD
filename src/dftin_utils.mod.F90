#ifdef __SR8000
!option MP(P(0)), LANGLVL(SAVE(0))
#endif
! ==================================================================

MODULE dftin_utils
  USE cp_gga_correlation_utils,        ONLY: cp_gga_c_param
  USE cp_gga_exchange_utils,           ONLY: cp_gga_x_param
  USE cp_xc_utils,                     ONLY: cp_xc,&
                                             cp_xc_mts_low_func_env,&
                                             cp_xc_mts_low_func,&
                                             cp_xc_functional_env,&
                                             cp_xc_functional,&
                                             cp_xc_kernel,&
                                             cp_xc_kernel_env,&
                                             cp_xc_max_nbr_funcs
  USE hubbardu,                        ONLY: hubbu,maxuatm,maxluatm 
  USE ener,                            ONLY: tenergy_ok
  USE error_handling,                  ONLY: stopgm
  USE func,                            ONLY: &
       ashcroft_coulomb_rcut, func1, func2, func3, &
       mfxcc_is_hedin, mfxcc_is_lyp, mfxcc_is_obpw, &
       mfxcc_is_obpz, mfxcc_is_pade, mfxcc_is_pw, mfxcc_is_pz, &
       mfxcc_is_skipped, mfxcc_is_vwn, mfxcc_is_wigner, mfxcx_is_skipped, &
       mfxcx_is_slaterx, mgcc_is_dfr_zpbec, mgcc_is_ggac, mgcc_is_hse, &
       mgcc_is_lyp, mgcc_is_optc, mgcc_is_pbec, mgcc_is_pbesc, &
       mgcc_is_perdew86, mgcc_is_skipped, mgcsrx_is_hse, mgcsrx_is_skipped, &
       mgcx_is_becke88, mgcx_is_dfr_xpbex, mgcx_is_dfr_xpbex_hybrid, &
       mgcx_is_dfr_zpbex, mgcx_is_ggax, mgcx_is_hcth, mgcx_is_optx, &
       mgcx_is_ox, mgcx_is_ox_hybrid, mgcx_is_pbesx, mgcx_is_pbex, &
       mgcx_is_revpbex, mgcx_is_skipped, mhfx_is_hartree, mhfx_is_hfx, &
       mhfx_is_skipped, msrx_is_CAM, msrx_is_ashcroft, msrx_is_erfc, &
       msrx_is_exp, msrx_is_skipped, mtau_is_skipped, mtau_is_tpss
  USE hfxmod,                          ONLY: hfxc3,&
                                             hfxc4,&
                                             hfxc5
  USE initclust_utils,                 ONLY: hf_init
  USE inscan_utils,                    ONLY: inscan
  USE kinds,                           ONLY: real_8
  USE linres,                          ONLY: lrf1,&
                                             lrf2,&
                                             lrf3,&
                                             lrf4
  USE mp_interface,                    ONLY: mp_bcast
  USE parac,                           ONLY: parai,&
                                             paral
  USE pw_hfx_input_cnst,               ONLY: hfx_dist_block_cyclic,&
                                             hfx_dist_dynamic
  USE readsr_utils,                    ONLY: input_string_len, & 
                                             keyword_contains, &
                                             readsi,&
                                             readsr,&
                                             xstring
  USE store_types,                     ONLY: cprint,iprint_ehub
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             cntr
  USE tbxc,                            ONLY: tabx,&
                                             toldcode
  USE vdwcmod,                         ONLY: empvdwc
  USE wann,                            ONLY: wannl
  USE zeroing_utils,                   ONLY: zeroing
  USE mts_utils,                       ONLY: mts

#include "sizeof.h"

  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER :: output_unit = 6

  PUBLIC :: dftin
  PUBLIC :: tdm_fun

CONTAINS

  ! ==================================================================
  SUBROUTINE dftin
    ! ==--------------------------------------------------------------==
    ! ==  THIS ROUTINE READS THE SECTION &DFT &END ON UNIT IUNIT      ==
    ! ==--------------------------------------------------------------==
    ! ==  THE INPUT HAS TO BE IN THE FOLLOWING FORMAT                 ==
    ! ==     &DFT                                                     ==
    ! ==        Options                                               ==
    ! ==     &END                                                     ==
    ! ==--------------------------------------------------------------==
    ! ==  LIST OF OPTIONS                                             ==
    ! ==                                                              ==
    ! ==    OLDCODE/NEWCODE                                           ==
    ! ==    LDA CORRELATION functionals                               ==
    ! ==    SLATER                                                    ==
    ! ==      alpha                                                   ==
    ! ==    GRADIENT CORRECTION functionals                           ==
    ! ==    FUNCTIONAL {LDA,BONLY,BP,BLYP,GGA/PW91,PBE,REVPBE,HCTH,   ==
    ! ==                OPTX,OLYP,XLYP}                               ==
    ! ==    FUNCTIONAL {B3LYP,OLDB3LYP,MODB3LYP,X3LYP,OLDX3LYP,       ==
    ! ==                MODX3LYP,B1LYP,OLDB1LYP,MODB1LYP,             ==
    ! ==                PBE0,PBE1W,PBES,REVPBE0}                      ==
    ! ==    FUNCTIONAL {TPSS}                                         ==
    ! ==    MTS_HIGH_FUNC (alias for FUNCTIONAL)                      ==
    ! ==    MTS_LOW_FUNC (same as for FUNCTIONAL)                     ==
    ! ==    HARTREE                                                   ==
    ! ==    HARTREE-FOCK                                              ==
    ! ==    ACM0                                                      ==
    ! ==    ACM1                                                      ==
    ! ==      a1                                                      ==
    ! ==    ACM3                                                      ==
    ! ==      a1 a2 a3                                                ==
    ! ==    LR KERNEL {NONE,LDA,...}                                  ==
    ! ==    REFUNCT {LDA,BP,BLYP,PBE,OLYP}                            ==
    ! ==    GC-CUTOFF                                                 ==
    ! ==      gceps                                                   ==
    ! ==    SMOOTH                                                    ==
    ! ==      mf sdelta                                               ==
    ! ==    [NO] EXCHANGE CORRELATION TABLE                           ==
    ! ==      narray rmaxxc                                           ==
    ! ==    BECKE BETA                                                ==
    ! ==      beta                                                    ==
    ! ==    SCREENED EXCHANGE [SELF] {EXP,ERFC}                       ==
    ! ==      srxa                                                    ==
    ! ==      srxa                                                    ==
    ! ==    WANNIER SCREENING [WFC,DENSITY,DIAG]                      ==
    ! ==      dwfc dwfmax                                             ==
    ! ==                                                              ==
    ! ==    LDA functionals : NO,PZ,VWN,LYP,PW,WIGNER,HEDIN,          ==
    ! ==                      OBPZ,OBPW,TETER                         ==
    ! ==    Correlation functionals : PZ,VWN,LYP,PW,HCTH              ==
    ! ==    Gradient correction functionals :                         ==
    ! ==       EXCH, CORREL, BECKE88, GGAX, PERDEW86, LYP, GGAC,      ==
    ! ==       PBEX, PBEC , HCTH/120, OPTX, OLYP,PBESC,PBESX          ==
    ! ==                                                              ==
    ! ==--------------------------------------------------------------==
    CHARACTER(*), PARAMETER                  :: procedureN = 'dftin'
    INTEGER, PARAMETER                       :: max_unknown_lines = 30 

    CHARACTER(len=input_string_len)          :: line, previous_line, &
                                                next_words, error_message, &
                                                xc_functional, xc_kernel, &
                                                xc_mts_low_func, unknown(max_unknown_lines) 
    INTEGER                                  :: ia, ie
    INTEGER                                  :: i, first, last, ierr, iunit, nbr_unknown_lines, &
                                                j, tr_c, tr_gc, tr_gx, tr_x
    INTEGER                                  :: cp_xc_functional_n_scale_control, cp_xc_functional_n_lib_control
    INTEGER                                  :: cp_xc_kernel_n_scale_control, cp_xc_kernel_n_lib_control
    INTEGER                                  :: cp_xc_mts_low_func_n_scale_control, cp_xc_mts_low_func_n_lib_control
    LOGICAL                                  :: cp_xc_functional_has_scales, cp_xc_functional_has_library
    LOGICAL                                  :: cp_xc_kernel_has_scales, cp_xc_kernel_has_library
    LOGICAL                                  :: cp_xc_mts_low_func_has_scales, cp_xc_mts_low_func_has_library
    LOGICAL                                  :: erread, force_cpfunc, force_div_analytical
    LOGICAL                                  :: force_div_numerical, force_libxc, force_new_hfx
    LOGICAL                                  :: overwrite_pcgc, overwrite_pclda, overwrite_phfx 
    LOGICAL                                  :: overwrite_pxgc, overwrite_pxlda
    LOGICAL                                  :: tr_code, treff, something_went_wrong, go_on_reading
    LOGICAL                                  :: get_functional, get_kernel, get_mts_low_func
    LOGICAL                                  :: functional_is_defined
    LOGICAL                                  :: set_via_gc_or_sx, set_via_acm_keyword, need_dft
    REAL(real_8)                             :: a1par, a2par, a3par, &
                                                new_pcgc, new_pclda, &
                                                new_phfx, new_pxgc, &
                                                new_pxlda, scalei


    !
    ! The read loop is only accessed by the io_parent, therefore, the
    ! io_parent checks within the loop have been removed (redundancy).
    !
    IF (paral%io_parent) THEN
       iunit = 5
       !
       ! Variables for reading
       !
       nbr_unknown_lines     = 0
       cp_xc_functional_has_scales      = .FALSE.
       cp_xc_functional_has_library     = .FALSE.
       cp_xc_functional_n_scale_control = 0
       cp_xc_functional_n_lib_control   = 0
       cp_xc_kernel_has_scales      = .FALSE.
       cp_xc_kernel_has_library     = .FALSE.
       cp_xc_kernel_n_scale_control = 0
       cp_xc_kernel_n_lib_control   = 0
       cp_xc_mts_low_func_has_scales      = .FALSE.
       cp_xc_mts_low_func_has_library     = .FALSE.
       cp_xc_mts_low_func_n_scale_control = 0
       cp_xc_mts_low_func_n_lib_control   = 0
       force_div_analytical  = .FALSE.
       force_div_numerical   = .FALSE.
       force_cpfunc          = .FALSE.
       force_libxc           = .FALSE.
       xc_functional         = ' '
       get_functional        = .FALSE.
       functional_is_defined = .FALSE.
       xc_kernel             = ' '
       get_kernel            = .FALSE.
       xc_mts_low_func       = ' '
       get_mts_low_func      = .FALSE.
       set_via_gc_or_sx      = .FALSE.
       set_via_acm_keyword   = .FALSE.
       line                  = ' '
       previous_line         = ' '
       error_message         = ' '
       !
       ! Defaults
       !
       func1%mfxcx  = -1 ! needed for GRADIENT CORRECTION keyword  
       func1%mfxcc  = -1 ! needed for GRADIENT CORRECTION keyword  
       func1%mgcx   = mgcx_is_skipped
       func1%mgcc   = mgcc_is_skipped
       func1%mhfx   = mhfx_is_skipped
       func1%mtau   = mtau_is_skipped
       func1%msrx   = msrx_is_skipped
       func1%mgcsrx = mgcsrx_is_skipped
       func2%salpha = 2.0_real_8/3.0_real_8
       func2%bbeta  = 0.0042_real_8
       func2%betapp = 0.0_real_8
       func2%srxa   = 0.0_real_8
       func3%pxlda  = 1.0_real_8
       func3%pxgc   = 1.0_real_8
       func3%pclda  = 1.0_real_8
       func3%pcgc   = 1.0_real_8
       func3%phfx   = 0.0_real_8
       new_phfx     = 0.0_real_8
       new_pxgc     = 0.0_real_8
       new_pcgc     = 0.0_real_8
       new_pxlda    = 0.0_real_8
       new_pclda    = 0.0_real_8
       cntl%use_xc_driver  = .FALSE.
       cntl%tgc            = .FALSE.
       cntl%tgcx           = .FALSE.
       cntl%tgcc           = .FALSE.
       cntl%ttau           = .FALSE.
       cntl%thybrid        = .FALSE.
       cntl%tsmooth        = .FALSE.
       cntl%div_analytical = .TRUE.
       cntl%use_scaled_hfx = .FALSE.
       cntl%thubb          = .FALSE.
       cntr%gceps          = 1.0e-8_real_8
       toldcode            = .FALSE.
       overwrite_phfx      = .FALSE.
       overwrite_pxgc      = .FALSE.
       overwrite_pcgc      = .FALSE.
       overwrite_pxlda     = .FALSE.
       overwrite_pclda     = .FALSE.
       ashcroft_coulomb_rcut = 6.0_real_8
       hfxc3%twscr             = .FALSE.
       hfxc3%twft              = .FALSE.
       hfxc3%twfc              = .FALSE.
       hfxc3%tdiaw             = .FALSE.
       hfxc3%keep_radius_fixed = .FALSE.
       hfxc3%use_new_hfx       = .TRUE.
       hfxc4%dwfc              = 1.e10_real_8
       hfxc4%dwf_integral_thresh = 0.0_real_8
       hfxc5%hfx_distribution  = hfx_dist_block_cyclic
       hfxc5%hfx_block_size    = 2
       hfxc5%recomp_two_int_list_every = HUGE(0)
       force_new_hfx = .FALSE.
       !
       treff = .FALSE.
       !
       lrf1%td_x  = -1
       lrf1%td_c  = -1
       lrf1%td_gx = -1
       lrf1%td_gc = -1
       lrf1%td_hf = -1
       lrf1%td_mtau = -1
       lrf2%td_ttau=.FALSE.
       lrf2%td_tgc=.FALSE.
       lrf2%td_tgcx=.FALSE.
       lrf2%td_tgcc=.FALSE.
       lrf2%td_hybrid=.FALSE.
       lrf2%td_code=.FALSE.
       lrf3%tdpxlda=1._real_8
       lrf3%tdpxgc=1._real_8
       lrf3%tdpclda=1._real_8
       lrf3%tdpcgc=1._real_8
       lrf3%tdphfx=0._real_8
       lrf4%td_functional=0
       !
       empvdwc%dft_func=" "
       !
       tenergy_ok=.TRUE.
       !
       tabx%narray=0
       tabx%rmaxxc=2.0_real_8
       tabx%rmaxbx=100.0_real_8
       !
       need_dft = .true.
       ! dft section is not always needed with the MTS scheme
       if (cntl%use_mts) need_dft = mts%low_level=='DFT' .or. mts%high_level=='DFT'
       !
       ! If the DFT section is not there, we simply move on.
       !
       ierr = inscan(iunit,'&DFT')
       !
       IF (ierr == 0) THEN
          !
          ! Main loop
          !
          go_on_reading = .true.
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
             ELSEIF (keyword_contains(line,'&END')) THEN
                go_on_reading = .FALSE.
             ELSEIF ( keyword_contains(line,'USE_CP',alias='CP_LIBRARY_ONLY') .OR. &
                      keyword_contains(line,'INTERNAL_ONLY') ) THEN
                cntl%use_xc_driver = .TRUE.
                cntl%tgc           = .TRUE.
                cntl%tgcx          = .TRUE.
                cntl%tgcc          = .TRUE.
      
                cp_xc_functional_env%via_libxc(:) = .FALSE.
                force_cpfunc = .TRUE.
                force_libxc  = .FALSE.
      
                func1%mfxcx = mfxcx_is_skipped
                func1%mgcx  = mgcx_is_skipped
                func1%mfxcc = mfxcc_is_skipped
                func1%mgcc  = mgcc_is_skipped
             ELSEIF ( keyword_contains(line,'USE_LIBXC',alias='LIBXC_ONLY') .OR. &
                      keyword_contains(line,'LIBXC_ONLY') ) THEN
                cntl%use_xc_driver = .TRUE.
                cntl%tgc           = .TRUE.
                cntl%tgcx          = .TRUE.
                cntl%tgcc          = .TRUE.
      
                cp_xc_functional_env%via_libxc(:) = .TRUE.
                force_cpfunc = .FALSE.
                force_libxc  = .TRUE.
      
                func1%mfxcx = mfxcx_is_skipped
                func1%mgcx  = mgcx_is_skipped
                func1%mfxcc = mfxcc_is_skipped
                func1%mgcc  = mgcc_is_skipped
             ELSEIF ( keyword_contains(line,'XC_DRIVER',alias='USE_XC_DRIVER') .OR. &
                      keyword_contains(line,'USE_DRIVER') ) THEN
                cntl%use_xc_driver = .TRUE.
                cntl%tgc           = .TRUE.
                cntl%tgcx          = .TRUE.
                cntl%tgcc          = .TRUE.
      
                func1%mfxcx = mfxcx_is_skipped
                func1%mgcx  = mgcx_is_skipped
                func1%mfxcc = mfxcc_is_skipped
                func1%mgcc  = mgcc_is_skipped
             ELSEIF (keyword_contains(line,'LIBRARY')) THEN
                CALL xstring(line,first,last)
                next_words = line(last+1:)
                i = 1
                library_loop: DO
                   CALL xstring(next_words,first,last)
                   IF (last - first < 0) EXIT library_loop
                   IF (i > cp_xc_max_nbr_funcs) CALL stopgm(procedureN,' LIBRARY: Request exceeds cp_xc_max_nbr_funcs, &
                                                                                  &increase its value in cp_xc_utils',&
                                                             __LINE__,__FILE__)
                   IF (keyword_contains(next_words(first:last),'CP',alias='INTERNAL')) THEN
                      IF (force_libxc) CALL stopgm(procedureN,' You may not choose a CP functional when forcing libxc ',&
                           __LINE__,__FILE__)
                      cp_xc_functional_env%via_libxc(i) = .FALSE.
                   ELSEIF (keyword_contains(next_words(first:last),'LXC',alias='LIBXC')) THEN
                      IF (force_cpfunc) CALL stopgm(procedureN,&
                           ' You may not choose a libxc functional when forcing CP/INTERNAL ',&
                           __LINE__,__FILE__)
                      cp_xc_functional_env%via_libxc(i) = .TRUE.
                   ELSE
                      CALL stopgm(procedureN,&
                           ' Unknown xc library (CP (INTERNAL) or LIBXC): '//TRIM(ADJUSTL(next_words(ia:ie))) ,&
                           __LINE__,__FILE__)
                   ENDIF
                   next_words = next_words(last+1:)
                   i = i+1
                ENDDO library_loop
                cp_xc_functional_n_lib_control = i-1
                cp_xc_functional_has_library   = .TRUE.

             ELSEIF (keyword_contains(line,'KERNEL_LIBRARY')) THEN
                CALL xstring(line,first,last)
                next_words = line(last+1:)
                i = 1
                kernel_library_loop: DO
                   CALL xstring(next_words,first,last)
                   IF (last - first < 0) EXIT kernel_library_loop
                   IF (i > cp_xc_max_nbr_funcs) CALL stopgm(procedureN,' LIBRARY: Request exceeds cp_xc_max_nbr_funcs, &
                                                                                  &increase its value in cp_xc_utils',&
                                                             __LINE__,__FILE__)
                   IF (keyword_contains(next_words(first:last),'CP',alias='INTERNAL')) THEN
                      IF (force_libxc) CALL stopgm(procedureN,' You may not choose a CP kernel when forcing libxc ',&
                           __LINE__,__FILE__)
                      cp_xc_kernel_env%via_libxc(i) = .FALSE.
                   ELSEIF (keyword_contains(next_words(first:last),'LXC',alias='LIBXC')) THEN
                      IF (force_cpfunc) CALL stopgm(procedureN,' You may not choose a libxc kernel when forcing CP/INTERNAL ',&
                           __LINE__,__FILE__)
                      cp_xc_kernel_env%via_libxc(i) = .TRUE.
                   ELSE
                      CALL stopgm(procedureN,' Unknown xc library (CP (INTERNAL) or LIBXC): '//TRIM(ADJUSTL(next_words(ia:ie))) ,&
                           __LINE__,__FILE__)
                   ENDIF
                   next_words = next_words(last+1:)
                   i = i+1
                ENDDO kernel_library_loop
                cp_xc_kernel_n_lib_control = i-1
                cp_xc_kernel_has_library   = .TRUE.

             ELSEIF (keyword_contains(line,'MTS_LOW_FUNC_LIBRARY')) THEN
                CALL xstring(line,first,last)
                next_words = line(last+1:)
                i = 1
                mts_low_func_library_loop: DO
                   CALL xstring(next_words,first,last)
                   IF (last - first < 0) EXIT mts_low_func_library_loop
                   IF (i > cp_xc_max_nbr_funcs) CALL stopgm(procedureN,' LIBRARY: Request exceeds cp_xc_max_nbr_funcs, &
                                                                                  &increase its value in cp_xc_utils',&
                                                             __LINE__,__FILE__)
                   IF (keyword_contains(next_words(first:last),'CP',alias='INTERNAL')) THEN
                      IF (force_libxc) CALL stopgm(procedureN,' You may not choose a CP MTS low functional when forcing libxc ',&
                           __LINE__,__FILE__)
                      cp_xc_mts_low_func_env%via_libxc(i) = .FALSE.
                   ELSEIF (keyword_contains(next_words(first:last),'LXC',alias='LIBXC')) THEN
                      IF (force_cpfunc) CALL stopgm(procedureN,&
                           ' You may not choose a libxc MTS low functional when forcing CP/INTERNAL ',&
                           __LINE__,__FILE__)
                      cp_xc_mts_low_func_env%via_libxc(i) = .TRUE.
                   ELSE
                      CALL stopgm(procedureN,' Unknown xc library (CP (INTERNAL) or LIBXC): '//TRIM(ADJUSTL(next_words(ia:ie))) ,&
                           __LINE__,__FILE__)
                   ENDIF
                   next_words = next_words(last+1:)
                   i = i+1
                ENDDO mts_low_func_library_loop
                cp_xc_mts_low_func_n_lib_control = i-1
                cp_xc_mts_low_func_has_library   = .TRUE.

             ELSEIF (keyword_contains(line,'SCALES',alias='SCALE',but_not='HFX')) THEN
                CALL xstring(line,first,last)
                next_words = line(last+1:)
                i=1
                scale_loop: DO 
                   erread=.FALSE.
                   CALL readsr(next_words,first,last,scalei,erread)
                   IF (erread) EXIT scale_loop
                   IF (i > cp_xc_max_nbr_funcs) CALL stopgm(procedureN,' SCALES: Request exceeds cp_xc_max_nbr_funcs, &
                                                                         &increase its value in cp_xc_utils',&
                                                             __LINE__,__FILE__)
                   cp_xc_functional_env%scales(i) = scalei
                   first = last
                   i = i+1
                ENDDO scale_loop
                cp_xc_functional_n_scale_control = i-1
                cp_xc_functional_has_scales     = .TRUE.

             ELSEIF (keyword_contains(line,'KERNEL_SCALES',alias='KERNEL_SCALE',but_not='HFX')) THEN
                CALL xstring(line,first,last)
                next_words = line(last+1:)
                i=1
                kernel_scale_loop: DO 
                   erread=.FALSE.
                   CALL readsr(next_words,first,last,scalei,erread)
                   IF (erread) EXIT kernel_scale_loop
                   IF (i > cp_xc_max_nbr_funcs) CALL stopgm(procedureN,' SCALES: Request exceeds cp_xc_max_nbr_funcs, &
                                                                         &increase its value in cp_xc_utils',&
                                                             __LINE__,__FILE__)
                   cp_xc_kernel_env%scales(i) = scalei
                   first = last
                   i = i+1
                ENDDO kernel_scale_loop
                cp_xc_kernel_n_scale_control = i-1
                cp_xc_kernel_has_scales     = .TRUE.

             ELSEIF (keyword_contains(line,'MTS_LOW_FUNC_SCALES',alias='MTS_LOW_FUNC_SCALE',but_not='HFX')) THEN
                CALL xstring(line,first,last)
                next_words = line(last+1:)
                i=1
                mts_low_func_scale_loop: DO 
                   erread=.FALSE.
                   CALL readsr(next_words,first,last,scalei,erread)
                   IF (erread) EXIT mts_low_func_scale_loop
                   IF (i > cp_xc_max_nbr_funcs) CALL stopgm(procedureN,' SCALES: Request exceeds cp_xc_max_nbr_funcs, &
                                                                         &increase its value in cp_xc_utils',&
                                                             __LINE__,__FILE__)
                   cp_xc_mts_low_func_env%scales(i) = scalei
                   first = last
                   i = i+1
                ENDDO mts_low_func_scale_loop
                cp_xc_mts_low_func_n_scale_control = i-1
                cp_xc_mts_low_func_has_scales     = .TRUE.

             ELSEIF (keyword_contains(line,'COMPATIBILITY',alias='OLD_DEFINITIONS')) THEN
                cp_xc_functional_env%use_compatibility_mode = .TRUE.
             ELSEIF(keyword_contains(line,'ANALYTICAL_DIV')) THEN
                force_div_analytical = .TRUE.
                force_div_numerical  = .FALSE.
             ELSEIF(keyword_contains(line,'NUMERICAL_DIV')) THEN
                force_div_analytical = .FALSE.
                force_div_numerical  = .TRUE.
             ELSEIF(keyword_contains(line,'PBE_FLEX_KAPPA')) THEN
                READ(iunit,*,iostat=ierr) cp_gga_x_param%pbe%kappa
                cp_gga_x_param%init = .true.
             ELSEIF(keyword_contains(line,'PBE_FLEX_MU')) THEN
                READ(iunit,*,iostat=ierr) cp_gga_x_param%pbe%mu
                cp_gga_x_param%init = .true.
             ELSEIF(keyword_contains(line,'PBE_FLEX_BETA')) THEN
                READ(iunit,*,iostat=ierr) cp_gga_c_param%pbe%pbe_beta
                cp_gga_c_param%init = .true.
             ELSEIF(keyword_contains(line,'PBE_FLEX_GAMMA')) THEN
                READ(iunit,*,iostat=ierr) cp_gga_c_param%pbe%pbe_gamma
                cp_gga_c_param%init = .true.
             ELSEIF(keyword_contains(line,'PBE_FLEX_UEG_CORRELATION')) THEN
                CALL xstring(line,first,last)
                next_words=line(last+1:)
                CALL xstring(next_words,first,last)
                cp_gga_c_param%pbe%lda_c = 'CP_'//TRIM(ADJUSTL(next_words(first:last)))
             ELSEIF (keyword_contains(line,'GRADIENT',and='CORRECTION')) THEN
                functional_is_defined = .TRUE.
                set_via_gc_or_sx      = .TRUE.
                cntl%tgc=.TRUE.
                IF (keyword_contains(line,'EXCH'))     cntl%tgcx  = .TRUE.
                IF (keyword_contains(line,'CORREL'))   cntl%tgcc  = .TRUE.
                IF (keyword_contains(line,'BECKE88'))  func1%mgcx = 1
                IF (keyword_contains(line,'PBEX'))     func1%mgcx=3
                IF (keyword_contains(line,'REVPBEX'))  func1%mgcx=4
                IF (keyword_contains(line,'PBESX'))    func1%mgcx=9
                IF (keyword_contains(line,'PERDEW86')) func1%mgcc=1
                IF (keyword_contains(line,'LYP'))      func1%mgcc=2
                IF (keyword_contains(line,'PBESC'))    func1%mgcc=7
                IF (keyword_contains(line,'GGAX',alias='PW91X'))   func1%mgcx=2
                IF (keyword_contains(line,'GGAC',alias='PW91C'))   func1%mgcc=3
                IF (keyword_contains(line,'PBEC',alias='REVPBEC')) func1%mgcc=4
                IF (keyword_contains(line,'HCTH')) THEN
                   toldcode=.TRUE.
                   func1%mfxcx=0
                   func1%mfxcc=0
                   func1%mgcx=5
                   func1%mgcc=5
                ENDIF
                IF (keyword_contains(line,'OPTX')) THEN
                   toldcode=.TRUE.
                   func1%mfxcx=0
                   func1%mgcx=6
                ENDIF
                IF (.NOT.cntl%tgcx.AND..NOT.cntl%tgcc.AND.func1%mgcx == 0.AND.func1%mgcc == 0) THEN
                   cntl%tgcx=.TRUE.
                   cntl%tgcc=.TRUE.
                ENDIF
                IF (cntl%tgcx.AND.func1%mgcx == 0) func1%mgcx=1     ! Set defaults only if not already set
                IF (cntl%tgcc.AND.func1%mgcc == 0) func1%mgcc=1
                IF (func1%mgcx > 0) cntl%tgcx=.TRUE.
                IF (func1%mgcc > 0) cntl%tgcc=.TRUE.

             ELSEIF(keyword_contains(line,'HARTREE',and='FOCK',alias='HARTREE-FOCK')) THEN
                ! TODO: Stand-alone functionality based on xc_driver
                functional_is_defined = .TRUE.
                cntl%tgc=.FALSE.
                cntl%tgcx=.FALSE.
                cntl%tgcc=.FALSE.
                func1%mfxcx = mfxcx_is_skipped
                func1%mfxcc=mfxcc_is_skipped
                func1%mgcx = mgcx_is_skipped
                func1%mgcc = mgcc_is_skipped
                func1%mhfx = mhfx_is_hfx
                func2%salpha=2._real_8/3._real_8
             ELSEIF (keyword_contains(line,'HARTREE',but_not='FOCK')) THEN
                ! TODO: Stand-alone functionality based on xc_driver
                functional_is_defined = .TRUE.
                cntl%tgc=.FALSE.
                cntl%tgcx=.FALSE.
                cntl%tgcc=.FALSE.
                func1%mfxcx = mfxcx_is_skipped
                func1%mfxcc=mfxcc_is_skipped
                func1%mgcx = mgcx_is_skipped
                func1%mgcc = mgcc_is_skipped
                func1%mhfx = mhfx_is_hartree
                func2%salpha=2._real_8/3._real_8
             ELSEIF (keyword_contains(line,'SLATER')) THEN
                ! TODO: Stand-alone functionality based on xc_driver
                ! ..Slater exchange
                IF (keyword_contains(line,'NO')) THEN
                   functional_is_defined = .TRUE.
                   set_via_gc_or_sx      = .TRUE.
                   func1%mfxcx = mfxcx_is_skipped
                ELSEIF (keyword_contains(line,'ALPHA')) THEN
                   ! Hack as long as we have the old code and xc_driver in parallel:
                   ! use of "alpha" only changes the value of salpha
                   READ(iunit,*,iostat=ierr) func2%salpha
                   IF (ierr /= 0) THEN
                      error_message        = 'COULD NOT READ SLATER ALPHA'
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ELSE
                   READ(iunit,*,iostat=ierr) func2%salpha
                   IF (ierr /= 0) THEN
                      error_message        = 'COULD NOT READ SLATER ALPHA'
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                   functional_is_defined = .TRUE.
                   set_via_gc_or_sx      = .TRUE.
                   func1%mfxcx = mfxcx_is_slaterx
                ENDIF
             ELSEIF(keyword_contains(line,'LDA',and='CORRELATION')) THEN
                functional_is_defined = .TRUE.
                set_via_gc_or_sx      = .TRUE.
                ! ..LDA Correlation functionals
                IF (keyword_contains(line,'NO')) THEN
                   func1%mfxcc = mfxcc_is_skipped
                ELSEIF (keyword_contains(line,'PZ')) THEN
                   func1%mfxcc = mfxcc_is_pz
                ELSEIF (keyword_contains(line,'VWN')) THEN
                   func1%mfxcc = mfxcc_is_vwn
                ELSEIF (keyword_contains(line,'LYP')) THEN
                   func1%mfxcc = mfxcc_is_lyp
                ELSEIF (keyword_contains(line,'PW')) THEN
                   func1%mfxcc = mfxcc_is_pw
                ELSEIF (keyword_contains(line,'WIGNER')) THEN
                   func1%mfxcc = mfxcc_is_wigner
                ELSEIF (keyword_contains(line,'HEDIN')) THEN
                   func1%mfxcc = mfxcc_is_hedin
                ELSEIF (keyword_contains(line,'OBPZ')) THEN
                   func1%mfxcc = mfxcc_is_obpz
                ELSEIF (keyword_contains(line,'OBPW')) THEN
                   func1%mfxcc = mfxcc_is_obpw
                ELSEIF (keyword_contains(line,'TETER')) THEN
                   func1%mfxcc = mfxcc_is_pade
                   func1%mfxcx = mfxcx_is_skipped
                ELSEIF (keyword_contains(line,'HCTH')) THEN
                   toldcode = .TRUE.
                   func1%mfxcx = mfxcx_is_skipped
                   func1%mfxcc = mfxcc_is_skipped
                   func1%mgcx  = mgcx_is_hcth
                   func1%mgcc  = mgcc_is_hse
                ENDIF
             ELSEIF (keyword_contains(line,'GC-CUTOFF')) THEN
                ! ..Cutoff for the density for GC functionals
                READ(iunit,*,iostat=ierr) cntr%gceps
                IF (ierr /= 0) THEN
                   error_message        = 'COULD NOT READ GC-CUTOFF VALUE'
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF (keyword_contains(line,'HFX_BLOCK_SIZE')) THEN
                READ(iunit,*,iostat=ierr) hfxc5%hfx_block_size
                IF (ierr /= 0) THEN
                   error_message        = 'COULD NOT READ HFX BLOCK SIZE'
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF (keyword_contains(line,'SCALED',and='EXCHANGE',alias='SCEX')) THEN
                cntl%use_scaled_hfx = .true.
             ELSEIF (keyword_contains(line,'HFX_DISTRIBUTION')) THEN
                IF (keyword_contains(line,'BLOCK_CYCLIC')) THEN
                   hfxc5%hfx_distribution=hfx_dist_block_cyclic
                ENDIF
                IF (keyword_contains(line,'DYNAMIC')) THEN
                   hfxc5%hfx_distribution=hfx_dist_dynamic
                ENDIF
                !
                ! Combined range separation GGA & HFX
             ELSEIF(keyword_contains(line,'RANGE',and='SEPARATION') .OR. &
                    keyword_contains(line,'CAM',alias='ATTENUATION',but_not='SCREENED')) THEN
                cntl%use_xc_driver      = .TRUE.
                cntl%div_analytical     = .TRUE.
                cp_xc_functional_env%overwrite_hfx = .TRUE.
                cp_xc_functional_env%get_hfx       = .TRUE.
                cp_xc_functional_env%hfx_operator  = 'CAM'
                cp_xc_functional_env%get_CAM_GGA   = .TRUE.
                READ(iunit,*,iostat=ierr) cp_xc_functional_env%hfx_screening, cp_xc_functional_env%hfx_constant, &
                                          cp_xc_functional_env%hfx_attenuated
                IF (ierr /= 0) THEN
                   error_message        = 'COULD NOT READ CAM PARAMETERS'
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                cp_xc_kernel_env%overwrite_hfx  = cp_xc_functional_env%overwrite_hfx     
                cp_xc_kernel_env%get_hfx        = cp_xc_functional_env%get_hfx           
                cp_xc_kernel_env%hfx_operator   = cp_xc_functional_env%hfx_operator      
                cp_xc_kernel_env%get_CAM_GGA    = cp_xc_functional_env%get_CAM_GGA       
                cp_xc_kernel_env%hfx_screening  = cp_xc_functional_env%hfx_screening 
                cp_xc_kernel_env%hfx_constant   = cp_xc_functional_env%hfx_constant  
                cp_xc_kernel_env%hfx_attenuated = cp_xc_functional_env%hfx_attenuated
                !
                cp_xc_mts_low_func_env%overwrite_hfx  = cp_xc_functional_env%overwrite_hfx     
                cp_xc_mts_low_func_env%get_hfx        = cp_xc_functional_env%get_hfx           
                cp_xc_mts_low_func_env%hfx_operator   = cp_xc_functional_env%hfx_operator      
                cp_xc_mts_low_func_env%get_CAM_GGA    = cp_xc_functional_env%get_CAM_GGA       
                cp_xc_mts_low_func_env%hfx_screening  = cp_xc_functional_env%hfx_screening 
                cp_xc_mts_low_func_env%hfx_constant   = cp_xc_functional_env%hfx_constant  
                cp_xc_mts_low_func_env%hfx_attenuated = cp_xc_functional_env%hfx_attenuated
                !
                ! Screening for HFX _ONLY_
             ELSEIF(keyword_contains(line,'SCREENED',and='EXCHANGE')) THEN
                ! ..Short range exact exchange
                func1%msrx              = msrx_is_exp
                cp_xc_functional_env%overwrite_hfx = .TRUE.
                cp_xc_functional_env%get_hfx       = .TRUE.
                cp_xc_functional_env%hfx_operator  = 'exp'
                IF (keyword_contains(line,'ASHCROFT')) THEN
                   func1%msrx = msrx_is_ashcroft
                   cp_xc_functional_env%hfx_operator = 'Ashcroft'
                ELSEIF (keyword_contains(line,'EXP')) THEN
                   func1%msrx = msrx_is_exp
                   cp_xc_functional_env%hfx_operator = 'exp'
                ELSEIF (keyword_contains(line,'ERFC')) THEN
                   func1%msrx = msrx_is_erfc
                   cp_xc_functional_env%hfx_operator = 'erfc'
                ELSEIF (keyword_contains(line,'CAM',alias='LC')) THEN
                   func1%msrx             = msrx_is_CAM
                   cp_xc_functional_env%hfx_operator = 'CAM'
                   cntl%div_analytical    = .TRUE.
                ENDIF
                !
                SELECT CASE(cp_xc_functional_env%hfx_operator)
                CASE("Ashcroft","ashcroft","ASHCROFT")
                   ! READ THE RADIUS HERE
                   READ(iunit,*,iostat=ierr) ashcroft_coulomb_rcut
                   IF (ierr /= 0) THEN
                      error_message        = 'COULD NOT READ ASHCROFT PARAMETERS'
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                CASE("CAM")
                   READ(iunit,*,iostat=ierr) cp_xc_functional_env%hfx_screening, cp_xc_functional_env%hfx_constant, &
                                             cp_xc_functional_env%hfx_attenuated
                   IF (ierr /= 0) THEN
                      error_message        = 'COULD NOT READ CAM PARAMETERS'
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                CASE DEFAULT
                   READ(iunit,*,iostat=ierr) cp_xc_functional_env%hfx_screening
                   IF (ierr /= 0) THEN
                      error_message        = 'COULD NOT READ PARAMETERS FOR SCREENED EXCHANGE'
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                END SELECT
                func2%srxa      = cp_xc_functional_env%hfx_screening
                func2%cam_alpha = cp_xc_functional_env%hfx_constant
                func2%cam_beta  = cp_xc_functional_env%hfx_attenuated
                !
                cp_xc_kernel_env%overwrite_hfx  = cp_xc_functional_env%overwrite_hfx     
                cp_xc_kernel_env%get_hfx        = cp_xc_functional_env%get_hfx           
                cp_xc_kernel_env%hfx_operator   = cp_xc_functional_env%hfx_operator      
                cp_xc_kernel_env%get_CAM_GGA    = cp_xc_functional_env%get_CAM_GGA       
                cp_xc_kernel_env%hfx_screening  = cp_xc_functional_env%hfx_screening 
                cp_xc_kernel_env%hfx_constant   = cp_xc_functional_env%hfx_constant  
                cp_xc_kernel_env%hfx_attenuated = cp_xc_functional_env%hfx_attenuated
                !
                cp_xc_mts_low_func_env%overwrite_hfx  = cp_xc_functional_env%overwrite_hfx     
                cp_xc_mts_low_func_env%get_hfx        = cp_xc_functional_env%get_hfx           
                cp_xc_mts_low_func_env%hfx_operator   = cp_xc_functional_env%hfx_operator      
                cp_xc_mts_low_func_env%get_CAM_GGA    = cp_xc_functional_env%get_CAM_GGA       
                cp_xc_mts_low_func_env%hfx_screening  = cp_xc_functional_env%hfx_screening 
                cp_xc_mts_low_func_env%hfx_constant   = cp_xc_functional_env%hfx_constant  
                cp_xc_mts_low_func_env%hfx_attenuated = cp_xc_functional_env%hfx_attenuated
             ELSEIF(keyword_contains(line,'BECKE',and='BETA')) THEN
                ! ..Parameter value for Becke exchange
                READ(iunit,*,iostat=ierr) func2%bbeta
                IF (ierr /= 0) THEN
                   error_message        = 'COULD NOT READ BECKE BETA'
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                func2%betapp=func2%bbeta
             ELSEIF (keyword_contains(line,'SMOOTHING',alias='SMOOTH')) THEN
                ! ..Smoothing of the density
                cntl%tsmooth=.TRUE.
                READ(iunit,*,iostat=ierr) cntr%smf,cntr%sdelta
                IF (ierr /= 0) THEN
                   error_message        = 'COULD NOT READ DENSITY SMOOTHING'
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF(keyword_contains(line,'EXCHANGE',and='TABLE')) THEN
                ! ..Tables for exchange and correlation
                IF (keyword_contains(line,'NO')) THEN
                   tabx%narray=0
                   tabx%rmaxxc=1.0_real_8
                ELSE
                   READ(iunit,*,iostat=ierr) tabx%narray,tabx%rmaxxc
                   IF (ierr /= 0) THEN
                      error_message        = 'COULD NOT READ EXCHANGE-CORRELATION TABLE'
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ENDIF
             ELSEIF (keyword_contains(line,'OLDCODE')) THEN
                toldcode=.TRUE.
             ELSEIF (keyword_contains(line,'NEWCODE')) THEN
                toldcode=.FALSE.
             ELSEIF (keyword_contains(line,'HFX_SCALE',alias='HFX_SCALES')) THEN
                cp_xc_functional_env%overwrite_hfx = .TRUE.
                cp_xc_functional_env%get_hfx       = .TRUE.
                CALL xstring(line,first,last)
                next_words = line(last+1:)
                READ(next_words,*, iostat=ierr) cp_xc_functional_env%hfx_constant
                IF (ierr /= 0) THEN
                   error_message        = 'COULD NOT READ HFX SCALING VALUE'
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF (keyword_contains(line,'KERNEL_HFX_SCALE',alias='KERNEL_HFX_SCALES')) THEN
                cp_xc_kernel_env%overwrite_hfx = .TRUE.
                cp_xc_kernel_env%get_hfx       = .TRUE.
                CALL xstring(line,first,last)
                next_words = line(last+1:)
                READ(next_words,*, iostat=ierr) cp_xc_kernel_env%hfx_constant
                IF (ierr /= 0) THEN
                   error_message        = 'COULD NOT READ HFX SCALING VALUE'
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF (keyword_contains(line,'MTS_LOW_FUNC_HFX_SCALE',alias='MTS_LOW_FUNC_HFX_SCALES')) THEN
                cp_xc_mts_low_func_env%overwrite_hfx = .TRUE.
                cp_xc_mts_low_func_env%get_hfx       = .TRUE.
                CALL xstring(line,first,last)
                next_words = line(last+1:)
                READ(next_words,*, iostat=ierr) cp_xc_mts_low_func_env%hfx_constant
                IF (ierr /= 0) THEN
                   error_message        = 'COULD NOT READ HFX SCALING VALUE'
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF (keyword_contains(line,'PHFX')) THEN
                ! ..For controlling only the amount of Hartree exchange
                cp_xc_functional_env%overwrite_hfx = .TRUE.
                cp_xc_functional_env%get_hfx       = .TRUE.
                READ(iunit,*, iostat=ierr) cp_xc_functional_env%hfx_constant
                overwrite_phfx = cp_xc_functional_env%overwrite_hfx
                new_phfx       = cp_xc_functional_env%hfx_constant
                IF (ierr /= 0) THEN
                   error_message        = 'COULD NOT READ HFX PERCENTAGE'
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF (keyword_contains(line,'PXGC')) THEN
                overwrite_pxgc=.TRUE.
                READ(iunit,*,iostat=ierr) new_pxgc
                IF (ierr /= 0) THEN
                   error_message        = 'COULD NOT READ GCX PERCENTAGE'
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF (keyword_contains(line,'PCGC')) THEN
                overwrite_pcgc=.TRUE.
                READ(iunit,*,iostat=ierr) new_pcgc
                IF (ierr /= 0) THEN
                   error_message        = 'COULD NOT READ GCC PERCENTAGE'
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF (keyword_contains(line,'PXLDA')) THEN
                overwrite_pxlda=.TRUE.
                READ(iunit,*,iostat=ierr) new_pxlda
                IF (ierr /= 0) THEN
                   error_message        = 'COULD NOT READ LDAX PERCENTAGE'
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF (keyword_contains(line,'PCLDA')) THEN
                overwrite_pclda=.TRUE.
                READ(iunit,*,iostat=ierr) new_pclda
                IF (ierr /= 0) THEN
                   error_message        = 'COULD NOT READ LDAC PERCENTAGE'
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF (keyword_contains(line,'ACM0')) THEN
                set_via_acm_keyword = .TRUE.
                cntl%thybrid = .TRUE.
                func3%pxlda  = 0.75_real_8
                func3%pxgc   = 0.75_real_8
                func3%pclda  = 1.0_real_8
                func3%pcgc   = 1.0_real_8
                func3%phfx   = 0.25_real_8
                func1%mhfx   = mhfx_is_hfx
             ELSEIF (keyword_contains(line,'ACM1')) THEN
                set_via_acm_keyword = .TRUE.
                READ(iunit,*,iostat=ierr) a1par
                IF (ierr /= 0) THEN
                   error_message        = 'COULD NOT READ ACM1 PARAMETER'
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                cntl%thybrid = .TRUE.
                func3%pxlda  = 1.0_real_8-a1par
                func3%pxgc   = 1.0_real_8-a1par
                func3%pclda  = 1.0_real_8
                func3%pcgc   = 1.0_real_8
                func3%phfx   = a1par
                func1%mhfx   = mhfx_is_hfx
             ELSEIF (keyword_contains(line,'ACM3')) THEN
                set_via_acm_keyword = .TRUE.
                READ(iunit,*,iostat=ierr) a1par,a2par,a3par
                IF (ierr /= 0) THEN
                   WRITE(output_unit,*) 'ERORR: COULD NOT READ ACM3 PARAMETERS'
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
                cntl%thybrid = .TRUE.
                func3%pxlda  = 1.0_real_8-a1par
                func3%pxgc   = 1.0_real_8-a2par
                func3%pclda  = 1.0_real_8
                func3%pcgc   = 1.0_real_8-a3par
                func3%phfx   = a1par
                func1%mhfx   = mhfx_is_hfx
             ELSEIF (keyword_contains(line,'FUNCTIONAL',alias='MTS_HIGH_FUNC')) THEN
                functional_is_defined = .true.
                xc_functional  = line
                get_functional = .true.
             ELSEIF (keyword_contains(line,'MTS_LOW_FUNC')) THEN
                functional_is_defined = .true.
                xc_mts_low_func  = line
                get_mts_low_func = .true.
             ELSEIF (keyword_contains(line,'REFUNCT')) THEN
                treff=.TRUE.
                tr_x=0
                tr_c=9
                tr_gx=0
                tr_gc=0
                tr_code=.FALSE.
                IF (keyword_contains(line,'LDA')) THEN
                   tr_x=0
                   tr_c=9
                   tr_gx=0
                   tr_gc=0
                   tr_code=.FALSE.
                ELSEIF (keyword_contains(line,'PBE')) THEN
                   tr_x=0
                   tr_c=9
                   tr_gx=3
                   tr_gc=4
                   tr_code=.FALSE.
                ELSEIF (keyword_contains(line,'BLYP')) THEN
                   tr_x=1
                   tr_c=3
                   tr_gx=1
                   tr_gc=2
                   tr_code=.FALSE.
                ELSEIF (keyword_contains(line,'OLYP')) THEN
                   tr_x=0
                   tr_c=3
                   tr_gx=6
                   tr_gc=2
                   tr_code=.TRUE.
                ELSEIF (keyword_contains(line,'BP',alias='BP86')) THEN
                   tr_x=0
                   tr_c=9
                   tr_gx=1
                   tr_gc=1
                   tr_code=.FALSE.
                ELSE
                   CALL stopgm(procedureN,'Reference functional not available',& 
                        __LINE__,__FILE__)
                ENDIF
             ELSEIF (keyword_contains(line,'LR',and='KERNEL',alias='LR_KERNEL') .OR. &
                     keyword_contains(line,'XC_KERNEL')) THEN
                xc_kernel  = line
                get_kernel = .true.
             ELSEIF(keyword_contains(line,'FORCE_NEW_HFX')) THEN
                force_new_hfx = .TRUE.
             ELSEIF(keyword_contains(line,'HFX',and='SCREENING',alias='HFX-SCREENING')) THEN
                hfxc3%twscr=.TRUE.
                hfxc3%use_new_hfx=.FALSE.
                IF (keyword_contains(line,'RADIUS')) THEN
                   hfxc3%twfc=.TRUE.
                   READ(iunit,*,iostat=ierr) hfxc4%dwfc
                   IF (ierr /= 0) THEN
                      error_message        = 'COULD NOT READ HFX SCREENING RADIUS'
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ELSEIF (keyword_contains(line,'DIAG')) THEN
                   hfxc3%tdiaw=.TRUE.
                ELSEIF (keyword_contains(line,'KEEP_RADIUS_FIXED')) THEN
                   hfxc3%keep_radius_fixed=.TRUE.
                ELSEIF (keyword_contains(line,'EPS_INT')) THEN
                   hfxc3%twft=.TRUE.
                   READ(iunit,*,iostat=ierr) hfxc4%dwf_integral_thresh
                   IF (ierr /= 0) THEN
                      error_message        = 'COULD NOT READ HFX EPS_INT THRESHOLD'
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ELSEIF (keyword_contains(line,'RECOMPUTE_TWO_INT_LIST_EVERY')) THEN
                   READ(iunit,*,iostat=ierr) hfxc5%recomp_two_int_list_every
                   IF (ierr /= 0) THEN
                      error_message        = 'COULD NOT READ HFX RECOMPUTE FREQUENCY'
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF
                ELSE
                   CALL stopgm(procedureN,'need to specify something...',& 
                        __LINE__,__FILE__)
                ENDIF
                IF (hfxc3%tdiaw) THEN
                   hfxc3%twft=.FALSE.
                   hfxc3%twfc=.FALSE.
                ENDIF
             ELSEIF(keyword_contains(LINE,'HUBBARD'))THEN
                   ! DFT+U
                   cntl%thubb=.TRUE.
                   cprint%iprint(iprint_ehub) = 1
                   hubbu%pfrqom=1
                   hubbu%uverb=.FALSE.
                   ! verbose output else occupation matrix in file OCCMAT
                   IF(INDEX(LINE,'VERB') /= 0)hubbu%uverb=.TRUE.
                   ! frequency of writing the occupation matrix
                   I=INDEX(LINE,'OCCMAT=')
                   IF(I /= 0)THEN
                     I=I+7
                     CALL READSI(LINE,I,J,hubbu%pfrqom,ERREAD)
                     IF(ERREAD)CALL STOPGM(procedureN,&
                               'ERROR READING OCCMAT= KEYWORD',&
__LINE__,__FILE__)
                   ENDIF
                   ! type of projection (treatment of projectors)
                   hubbu%portho=.false.
                   hubbu%pnorm=.false.
                   IF((INDEX(LINE,'NORM') /= 0)&
                     .AND.(INDEX(LINE,'ORTHO') /= 0))THEN
                     hubbu%portho=.true.
                     hubbu%pnorm=.true.
                   ELSEIF(INDEX(LINE,'NORM') /= 0)THEN
                     hubbu%pnorm=.true.
                   ELSEIF(INDEX(LINE,'ORTHO') /= 0)THEN
                     hubbu%portho=.true.
                   ENDIF
                   ! number of U atoms
                   I=INDEX(LINE,'NUATM=')
                   IF(I /= 0)THEN
                     I=I+6
                     CALL READSI(LINE,I,J,hubbu%NUATM,ERREAD)
                     IF(ERREAD)CALL STOPGM(procedureN,&
                               'ERROR READING DFTIN= KEYWORD', &
__LINE__,__FILE__)
                   ELSE
                      WRITE(output_unit,*)'!!! WARNING! DEFAULT NUATM=1 IS USED !!'
                      hubbu%NUATM=1
                   ENDIF
                   IF(hubbu%NUATM > MAXUATM)&
                      CALL STOPGM(procedureN,'Increase MAXUATM in hubbardu',&
 __LINE__,__FILE__)
                   DO I=1,hubbu%NUATM
                      !    reading U-atom number, U parameter, L.R. parameter alpha, 
                      !    number of different angular momentum for projectors
                      IF(cntl%BSYMM)THEN
                         ! N. Siemer: Hubbard U + Broken Symmetry  was implemented in the Bochumer
                         ! version r1887. Upon reimplementation in CPMD 4 it was not
                         ! checked, therefore STOPGM.
                         ! contact: niklas.siemer@theochem.rub.de
                         CALL STOPGM(procedureN,&
                         "Broken Symmetry + Hubbard U not implemented",&
 __LINE__,__FILE__)
                         READ(IUNIT,*)&
                             hubbu%uatm(I),hubbu%u(I),hubbu%hs(I),hubbu%a(I),hubbu%nl(I)
                      ELSE
                         READ(IUNIT,*)&
                             hubbu%UATM(I),hubbu%U(I),hubbu%A(I),hubbu%NL(I)
                      ENDIF
                      IF(hubbu%NL(I) > MAXLUATM)&
                         CALL STOPGM(procedureN,'Increase MAXLUATM in hubbardu', &
__LINE__,__FILE__)
                      DO J=1,hubbu%NL(I)
                         ! reading the shell number of the projector and angular momentum l
                         READ(IUNIT,*)hubbu%S(I,J),hubbu%L(I,J)
                      ENDDO
                   ENDDO
             ELSE
                ! Unknown Keyword. store and continue
                IF (' ' /= line) THEN
                   IF (nbr_unknown_lines < max_unknown_lines) THEN
                      nbr_unknown_lines = nbr_unknown_lines+1
                      unknown(nbr_unknown_lines) = line
                   ELSE
                      DO i=1,max_unknown_lines-1
                         unknown(i) = unknown(i+1)
                      ENDDO
                      unknown(nbr_unknown_lines) = line
                   ENDIF
                ENDIF
             ENDIF
          END DO
          !
          ! End of read-loop
          !
          IF (.not. functional_is_defined) THEN
             something_went_wrong = .true.
             error_message        = 'NO FUNCTIONAL HAS BEEN DEFINED - DEFAULTS NO LONGER APPLY'
          ELSE IF (set_via_gc_or_sx .and. cntl%use_xc_driver) THEN
             something_went_wrong = .true.
             error_message        = 'XC_DRIVER: XC FUNCTIONAL CAN ONLY BE SET VIA KEYWORD "FUNCTIONAL"'
          ELSE IF (set_via_acm_keyword .and. cntl%use_xc_driver) THEN
             something_went_wrong = .true.
             error_message        = 'XC_DRIVER: ACM HAS TO BE SET UP MANUALLY VIA KEYWORD "SCALES"'
          END IF
       ELSE IF (need_dft) then
           something_went_wrong = .true.
           error_message        = 'MISSING &DFT SECTION - SECTION MANDATORY'
       ENDIF
       !
       IF (something_went_wrong) THEN
           WRITE(output_unit,'(/,1X,64("!"))')
           WRITE(output_unit,'(1X, A, 1X, A)') 'ERROR:', 'PROBLEM WHILE READING &DFT SECTION:' 
           WRITE(output_unit,'(8X, A)') TRIM(ADJUSTL(error_message))
           IF (line /= ' ' .or. previous_line /= ' ') THEN
              WRITE(output_unit,'(8X, A)') 'THE LAST TWO LINES READ WITHIN THE SECTION WERE:'
              WRITE(output_unit,'(/,1X, A)') TRIM(ADJUSTL(previous_line))
              WRITE(output_unit,'(1X, A)') TRIM(ADJUSTL(line))
           END IF
           WRITE(output_unit,'(1X,64("!"))')
           CALL stopgm(procedureN,'Error while reading &DFT section, cf. output file',& 
                __LINE__,__FILE__)
       ENDIF
       !
       ! Now that we know which functional to use, set it
       !
       if (get_functional) CALL set_xc_functional(xc_functional)
       if (get_kernel)     CALL set_xc_kernel(xc_kernel)
       if (get_mts_low_func)  CALL set_xc_mts_low_func(xc_mts_low_func)

       ! 
       ! MTS: set names for high and low level functionals in mts structure
       ! 
       if (cntl%use_mts) then
          if (mts%low_level=='DFT') mts%low_dft_func = xc_mts_low_func
          if (mts%high_level=='DFT') mts%high_dft_func = xc_functional
       end if
       !
       ! Check for conflicting options and do some variable magic
       !
       CALL check_options()

       !
       ! Output
       !
       IF (cntl%use_xc_driver) THEN

          ! Sanity checks
          IF (cp_xc_functional_env%init) THEN
             IF ((cp_xc_functional_has_scales) .and. &
                (cp_xc_functional_n_scale_control /= cp_xc_functional_env%n_funcs)) CALL stopgm(procedureN,'Number '&
                &//'of scaling values and functionals does not match',& 
                __LINE__,__FILE__)
             IF ((cp_xc_functional_has_library) .and. &
                (cp_xc_functional_n_lib_control /= cp_xc_functional_env%n_funcs)) CALL stopgm(procedureN,'Number ' &
                //'of library entries and functionals does not match',& 
                __LINE__,__FILE__)
             CALL cp_xc_functional_env%report( func2%salpha, cntr%gceps, output_unit )
          ENDIF

          IF (cp_xc_kernel_env%init) THEN
             IF ((cp_xc_kernel_has_scales) .and. &
                 (cp_xc_kernel_n_scale_control /= cp_xc_kernel_env%n_funcs)) CALL stopgm(procedureN,'Number '&
                 &//'of scaling values and kernels does not match',& 
                 __LINE__,__FILE__)
             IF ((cp_xc_kernel_has_library) .and. &
                 (cp_xc_kernel_n_lib_control /= cp_xc_kernel_env%n_funcs)) CALL stopgm(procedureN,'Number ' &
                 //'of library entries and kernels does not match',& 
                 __LINE__,__FILE__)

             CALL cp_xc_kernel_env%report( func2%salpha, cntr%gceps, output_unit ) 
          ENDIF

          IF (cp_xc_mts_low_func_env%init) THEN
             IF ((cp_xc_mts_low_func_has_scales) .and. &
                 (cp_xc_mts_low_func_n_scale_control /= cp_xc_mts_low_func_env%n_funcs)) CALL stopgm(procedureN,'Number '&
                 &//'of scaling values and MTS low functionals does not match',& 
                 __LINE__,__FILE__)
             IF ((cp_xc_mts_low_func_has_library) .and. &
                 (cp_xc_mts_low_func_n_lib_control /= cp_xc_mts_low_func_env%n_funcs)) CALL stopgm(procedureN,'Number ' &
                 //'of library entries and MTS low functionals does not match',& 
                 __LINE__,__FILE__)

             CALL cp_xc_mts_low_func_env%report( func2%salpha, cntr%gceps, output_unit ) 
          ENDIF

       ELSE
          CALL newcode_oldcode_report()
       ENDIF

       CALL hfx_report()

       IF (cntl%tsmooth) WRITE(output_unit,'(A,F6.4,A,F6.4)') ' SMOOTHING OF DENSITY: ALPHA=', cntr%smf,'  BETA=',cntr%sdelta

       IF (.not. cntl%use_xc_driver) CALL tddft_newcode_oldcode_report()

       IF (nbr_unknown_lines /= 0) THEN
          WRITE(output_unit,'(/,1X,64("="))')
          WRITE(output_unit,'(1X,A,14X,A,14X,A)') '= ','UNKNOWN KEYWORDS IN SECTION &DFT','='
          DO i=1,nbr_unknown_lines
             previous_line = unknown(i)
             CALL xstring(previous_line,first,last)
             WRITE(output_unit,'(A,A)') ' = ',previous_line(first:last)
          ENDDO
          WRITE(output_unit,'(1X,64("="),/)')
       ENDIF

    ENDIF

    CALL broadcast_dftin()

    !
    ! Create xc driver environment
    !
    IF( cntl%use_xc_driver ) THEN

       IF (cp_xc_functional_env%init) THEN
          CALL cp_xc_functional%create( cp_xc_functional_env )
       ELSE
          IF (cntl%use_mts .and. mts%high_level/='DFT') THEN
             !
             ! If the functional environment is not set, we set functional = mts_low_func
             !
             IF (cp_xc_mts_low_func_env%init) THEN
                CALL cp_xc_functional%create( cp_xc_mts_low_func_env )
             ELSE
                CALL stopgm(procedureN,'cp_xc_mts_low_func_env has not been initialised',&
                   __LINE__,__FILE__)
             ENDIF
          ELSE

             CALL stopgm(procedureN,'cp_xc_functional_env has not been initialised',&
                __LINE__,__FILE__)
          END IF
       ENDIF
       !
       IF (cp_xc_kernel_env%init) THEN
          CALL cp_xc_kernel%create( cp_xc_kernel_env )
       ELSE
          !
          ! If the kernel environment is not set, we set kernel = functional
          !
          CALL cp_xc_kernel%create( cp_xc_functional_env )
       ENDIF
       !
       if (cntl%use_mts .and. mts%low_level == 'DFT') then
          IF (cp_xc_mts_low_func_env%init) THEN
             CALL cp_xc_mts_low_func%create( cp_xc_mts_low_func_env )
          ELSE
             CALL stopgm(procedureN,'cp_xc_mts_low_func_env has not been initialised',&
                __LINE__,__FILE__)
          ENDIF
       end if
       !
       ! Copy back values from cp_xc types into ex-COMMONs
       !
       cp_xc => cp_xc_functional
       CALL cp_xc%get( tgc=cntl%tgc, tgcx=cntl%tgcx, tgcc=cntl%tgcc, ttau=cntl%ttau, &
                       thybrid=cntl%thybrid, mhfx=func1%mhfx, phfx=func3%phfx, &
                       msrx=func1%msrx, srxa=func2%srxa, cam_alpha=func2%cam_alpha, cam_beta=func2%cam_beta )
    ENDIF

    CONTAINS

    ! ==--------------------------------------------------------------==
    SUBROUTINE set_xc_functional(line)
       !
       ! TODO: Once we get rid of the old xc code, move to cp_xc_utils


       CHARACTER(len=*), INTENT(in) :: line
       
       IF( cntl%use_xc_driver ) THEN

          cp_xc_functional_env%init = .TRUE.
          cp_xc_functional_env%polarized = cntl%tlsd
          CALL xstring(line,first,last)
          next_words = line(last+1:)
          i=1
          DO
             CALL xstring(next_words,first,last)
             IF (last-first < 0) EXIT
             IF (i > cp_xc_max_nbr_funcs) CALL stopgm(procedureN,' FUNCTIONAL: Request exceeds cp_xc_max_nbr_funcs, &
                                                                  &increase its value in cp_xc_utils',&
                                                      __LINE__,__FILE__)
             cp_xc_functional_env%funcs(i) = next_words(first:last)
             next_words = next_words(last+1:)
             i = i+1
          ENDDO
          cp_xc_functional_env%n_funcs = i-1
          !
          ! these are always 0 for xc_driver functionals
          !
          func1%mfxcc = mfxcc_is_skipped
          func1%mfxcx = mfxcx_is_skipped

          !
          ! Add HFX for hybrids
          !
          CALL cp_xc_functional_env%set_hybrid_defaults(default_divergence=cntl%div_analytical)

          !
          ! Postprocess the VDW string for Grimme
          !
          IF (cp_xc_functional_env%n_funcs == 1) THEN
             IF     (cp_xc_functional_env%funcs(1) == 'GGA_XC_BLYP') THEN
                empvdwc%dft_func='BLYP'
             ELSEIF (cp_xc_functional_env%funcs(1) == 'HYB_GGA_XC_B3LYP') THEN
                empvdwc%dft_func='B3LYP'
             ELSEIF (cp_xc_functional_env%funcs(1) == 'GGA_XC_PBE') THEN
                empvdwc%dft_func='PBE'
             ELSEIF (cp_xc_functional_env%funcs(1) == 'HYB_GGA_XC_PBE0') THEN
                empvdwc%dft_func='PBE0'
             ELSEIF (cp_xc_functional_env%funcs(1) == 'GGA_XC_BP' .OR. &
                     cp_xc_functional_env%funcs(1) == 'GGA_XC_BP86') THEN
                empvdwc%dft_func='BP86'
             ELSEIF (cp_xc_functional_env%funcs(1) == 'MGGA_XC_TPSS') THEN
                empvdwc%dft_func='TPSS'
             ELSEIF (cp_xc_functional_env%funcs(1) == 'GGA_XC_PBE_R' .OR. &
                     cp_xc_functional_env%funcs(1) == 'GGA_XC_REVPBE') THEN
                empvdwc%dft_func='REVPBE'
             ELSEIF (cp_xc_functional_env%funcs(1) == 'HYB_GGA_XC_PBE0_R' .OR. &
                     cp_xc_functional_env%funcs(1) == 'HYB_GGA_XC_REVPBE0') THEN
                empvdwc%dft_func='REVPBE0'
             ELSEIF (cp_xc_functional_env%funcs(1) == 'HYB_GGA_XC_CAM_B3LYP') THEN
                empvdwc%dft_func='CAM-B3LYP'
             ELSEIF (cp_xc_functional_env%funcs(1) == 'GGA_XC_HCTH_120') THEN
                empvdwc%dft_func='HCTH'
             ELSEIF (cp_xc_functional_env%funcs(1) == 'GGA_XC_OPBE') THEN
                empvdwc%dft_func='OPBE'
             ENDIF
          ELSEIF (cp_xc_functional_env%n_funcs == 2) THEN
             IF     ( any(cp_xc_functional_env%funcs(:) == 'GGA_X_B88') .AND. &
                      any(cp_xc_functional_env%funcs(:) == 'GGA_C_LYP') ) THEN
                empvdwc%dft_func='BLYP'
             ELSEIF ( any(cp_xc_functional_env%funcs(:) == 'GGA_X_PBE') .AND. &
                      any(cp_xc_functional_env%funcs(:) == 'GGA_C_PBE') ) THEN
                empvdwc%dft_func='PBE'
             ELSEIF ( any(cp_xc_functional_env%funcs(:) == 'GGA_X_B88') .AND. &
                      any(cp_xc_functional_env%funcs(:) == 'GGA_C_P86') ) THEN
                empvdwc%dft_func='BP86'
             ELSEIF ( any(cp_xc_functional_env%funcs(:) == 'MGGA_X_TPSS') .AND. &
                      any(cp_xc_functional_env%funcs(:) == 'MGGA_C_TPSS') ) THEN
                empvdwc%dft_func='TPSS'
             ELSEIF ( any(cp_xc_functional_env%funcs(:) == 'GGA_X_PBE_R') .AND. &
                      any(cp_xc_functional_env%funcs(:) == 'GGA_C_PBE') ) THEN
                empvdwc%dft_func='REVPBE'
             ELSEIF ( any(cp_xc_functional_env%funcs(:) == 'GGA_X_PW91') .AND. &
                      any(cp_xc_functional_env%funcs(:) == 'GGA_C_PW91') ) THEN
                empvdwc%dft_func='PW91'
             ELSEIF ( any(cp_xc_functional_env%funcs(:) == 'GGA_X_PBE_SOL') .AND. &
                      any(cp_xc_functional_env%funcs(:) == 'GGA_C_PBE_SOL') ) THEN
                empvdwc%dft_func='PBES'
             ELSEIF ( any(cp_xc_functional_env%funcs(:) == 'GGA_X_OPTX') .AND. &
                      any(cp_xc_functional_env%funcs(:) == 'GGA_C_LYP') ) THEN
                empvdwc%dft_func='OLYP'
             ELSEIF ( any(cp_xc_functional_env%funcs(:) == 'GGA_X_OPTX') .AND. &
                      any(cp_xc_functional_env%funcs(:) == 'GGA_C_PBE') ) THEN
                empvdwc%dft_func='OPBE'
             ENDIF
          ENDIF

       ELSE
          !
          ! Oldcode & newcode
          ! 
          IF (keyword_contains(line,'NONE')) THEN
             cntl%tgc=.FALSE.
             cntl%tgcx=.FALSE.
             cntl%tgcc=.FALSE.
             func1%mfxcx = mfxcx_is_skipped
             func1%mfxcc=mfxcc_is_skipped
             func1%mgcx = mgcx_is_skipped
             func1%mgcc = mgcc_is_skipped
             func1%mhfx = mhfx_is_skipped
             func2%salpha=2._real_8/3._real_8
          ELSEIF (keyword_contains(line,'SONLY')) THEN
             cntl%tgc=.FALSE.
             cntl%tgcx=.FALSE.
             cntl%tgcc=.FALSE.
             func1%mfxcx = mfxcx_is_slaterx
             func1%mfxcc=mfxcc_is_skipped
             func1%mgcx = mgcx_is_skipped
             func1%mgcc = mgcc_is_skipped
             func1%mhfx = mhfx_is_skipped
             func2%salpha=2._real_8/3._real_8
          ELSEIF (keyword_contains(line,'LDA')) THEN
             cntl%tgc=.FALSE.
             cntl%tgcx=.FALSE.
             cntl%tgcc=.FALSE.
             func1%mfxcx = mfxcx_is_skipped
             func1%mfxcc=mfxcc_is_pade
             func1%mgcx = mgcx_is_skipped
             func1%mgcc = mgcc_is_skipped
             func1%mhfx = mhfx_is_skipped
             func2%salpha=2._real_8/3._real_8
          ELSEIF (keyword_contains(line,'BONLY')) THEN
             cntl%tgc=.TRUE.
             cntl%tgcx=.TRUE.
             cntl%tgcc=.FALSE.
             func1%mfxcx = mfxcx_is_skipped
             func2%salpha=2._real_8/3._real_8
             func1%mfxcc=mfxcc_is_pade
             func1%mgcx = mgcx_is_becke88
             func1%mgcc = mgcc_is_skipped
             func1%mhfx = mhfx_is_skipped
          ELSEIF (keyword_contains(line,'BP',alias='BP86')) THEN
             empvdwc%dft_func='BP86'
             cntl%tgc=.TRUE.
             cntl%tgcx=.TRUE.
             cntl%tgcc=.TRUE.
             func1%mfxcx = mfxcx_is_skipped
             func2%salpha=2._real_8/3._real_8
             func1%mfxcc=mfxcc_is_pade
             func1%mgcx = mgcx_is_becke88
             func1%mgcc = mgcc_is_perdew86
             func1%mhfx = mhfx_is_skipped
          ELSEIF (keyword_contains(line,'BLYP')) THEN
             cntl%tgc=.TRUE.
             cntl%tgcx=.TRUE.
             cntl%tgcc=.TRUE.
             func1%mfxcx = mfxcx_is_slaterx
             empvdwc%dft_func='BLYP'
             func2%salpha=2._real_8/3._real_8
             func1%mfxcc=mfxcc_is_lyp
             func1%mgcx = mgcx_is_becke88
             func1%mgcc = mgcc_is_lyp
             func1%mhfx = mhfx_is_skipped
          ELSEIF (keyword_contains(line,'XLYP')) THEN
             toldcode=.TRUE.
             cntl%tgc=.TRUE.
             cntl%tgcx=.TRUE.
             cntl%tgcc=.TRUE.
             func1%mfxcx = mfxcx_is_slaterx
             func2%salpha=2._real_8/3._real_8
             func1%mfxcc=mfxcc_is_lyp
             func1%mgcx = mgcx_is_ox
             func1%mgcc = mgcc_is_lyp
             func1%mhfx = mhfx_is_skipped
          ELSEIF(keyword_contains(line,'GGA',alias='PW91')) THEN
             empvdwc%dft_func='PW91'
             cntl%tgc=.TRUE.
             cntl%tgcx=.TRUE.
             cntl%tgcc=.TRUE.
             func1%mfxcx = mfxcx_is_skipped
             func2%salpha=2._real_8/3._real_8
             func1%mfxcc=mfxcc_is_pade
             func1%mgcx = mgcx_is_ggax
             func1%mgcc = mgcc_is_ggac
             func1%mhfx = mhfx_is_skipped
          ELSEIF (keyword_contains(line,'PBE1W')) THEN
             toldcode=.TRUE.
             cntl%tgc=.TRUE.
             cntl%tgcx=.TRUE.
             cntl%tgcc=.TRUE.
             func1%mfxcx = mfxcx_is_slaterx
             func2%salpha=2._real_8/3._real_8
             func1%mfxcc=mfxcc_is_pz
             func1%mgcx = mgcx_is_pbex
             func1%mgcc = mgcc_is_optc
             func1%mhfx = mhfx_is_skipped
          ELSEIF (keyword_contains(line,'REVPBE0')) THEN
             empvdwc%dft_func='REVPBE0'
             toldcode=.TRUE.
             cntl%tgc=.TRUE.
             cntl%tgcx=.TRUE.
             cntl%tgcc=.TRUE.
             cntl%thybrid=.TRUE.
             func1%mfxcx = mfxcx_is_slaterx
             func2%salpha=2.0_real_8/3.0_real_8
             func1%mfxcc=mfxcc_is_vwn
             func1%mgcx = mgcx_is_revpbex
             func1%mgcc = mgcc_is_pbec
             func1%mhfx = mhfx_is_hfx
             func3%pxlda=0.75_real_8
             func3%pxgc=0.75_real_8
             func3%pclda=1._real_8
             func3%pcgc=1._real_8
             func3%phfx=0.25_real_8
          ELSEIF (keyword_contains(line,'REVPBE')) THEN
             empvdwc%dft_func='REVPBE'
             cntl%tgc=.TRUE.
             cntl%tgcx=.TRUE.
             cntl%tgcc=.TRUE.
             func1%mfxcx = mfxcx_is_skipped
             func2%salpha=2._real_8/3._real_8
             func1%mfxcc=mfxcc_is_pade
             func1%mgcx = mgcx_is_revpbex
             func1%mgcc = mgcc_is_pbec
             func1%mhfx = mhfx_is_skipped
          ELSEIF (keyword_contains(line,'ZPBE0')) THEN
             toldcode=.TRUE.
             cntl%tgc=.TRUE.
             cntl%tgcx=.TRUE.
             cntl%tgcc=.TRUE.
             cntl%thybrid=.TRUE.
       
             func1%mfxcx = mfxcx_is_skipped ! 1 (1=LDAx l-part)
             func2%salpha=0.0_real_8! 2.0_real_8/3.0_real_8
             func1%mgcx = mgcx_is_dfr_zpbex ! 10 (10=ZPBEx+PBEx)
             func1%mfxcc=mfxcc_is_skipped ! 2 (2=VWNc l-part, 9=PADExc l-part)
             func1%mgcc = mgcc_is_dfr_zpbec  ! 4 (4=PBEc g-part, 8=PBEc external)
             func1%msrx = msrx_is_ashcroft  ! ASHCROFT EXCHANGE
       
             func1%mhfx = mhfx_is_hfx
             func3%pxlda=0.75_real_8
             func3%pxgc=0.75_real_8
             func3%pclda=1._real_8
             func3%pcgc=1._real_8
             func3%phfx=0.25_real_8
          ELSEIF (keyword_contains(line,'XPBE0')) THEN
             toldcode=.TRUE.
             cntl%tgc=.TRUE.
             cntl%tgcx=.TRUE.
             cntl%tgcc=.TRUE.
             cntl%thybrid=.TRUE.
       
             func1%mfxcx = mfxcx_is_skipped ! 1 (1=LDAx l-part)
             func2%salpha=0.0_real_8! 2.0_real_8/3.0_real_8
             func1%mgcx = mgcx_is_dfr_xpbex_hybrid ! 12 (10=XPBEx+PBEx)
             func1%mfxcc=mfxcc_is_skipped ! 2 (2=VWNc l-part, 9=PADExc l-part)
             func1%mgcc = mgcc_is_dfr_zpbec  ! 4 (4=PBEc g-part, 8=PBEc external)
             func1%msrx = msrx_is_ashcroft  ! ASHCROFT EXCHANGE
       
             func1%mhfx = mhfx_is_hfx
             func3%pxlda=0.75_real_8
             func3%pxgc=0.75_real_8
             func3%pclda=1._real_8
             func3%pcgc=1._real_8
             func3%phfx=0.25_real_8
          ELSEIF (keyword_contains(line,'PBEoriginal')) THEN
       
             CALL stopgm(procedureN,'PBEoriginal isnt supported anymore',&
                  __LINE__,__FILE__)
       
             cntl%tgc=.TRUE.
             cntl%tgcx=.TRUE.
             cntl%tgcc=.TRUE.
             func1%mfxcx = mfxcx_is_skipped
             !                                  !vw old       ! with DFRepository.F
             func2%salpha=2.0_real_8/3.0_real_8 !vw 2._real_8/3._real_8 ! 0.0_real_8
             func1%mfxcc=mfxcc_is_pade            !vw 9         ! 0
             func1%mgcx = mgcx_is_pbex             !vw 3         ! 11
             func1%mgcc = mgcc_is_pbec             !vw 4         ! 8
             func1%mhfx = mhfx_is_skipped
          ELSEIF (keyword_contains(line,'PBE0')) THEN
             empvdwc%dft_func='PBE0'
             toldcode=.TRUE.
             cntl%tgc=.TRUE.
             cntl%tgcx=.TRUE.
             cntl%tgcc=.TRUE.
             cntl%thybrid=.TRUE.
             !vw old       ! with DFRepository.F 
             func1%mfxcx = mfxcx_is_slaterx              !vw 1         ! 0
             func2%salpha=2.0_real_8/3.0_real_8   !vw 2._real_8/3._real_8 ! 0.0_real_8
             func1%mfxcc=mfxcc_is_vwn              !vw 2         ! 0
             func1%mgcx = mgcx_is_pbex               !vw 3         ! 11
             func1%mgcc = mgcc_is_pbec               !vw 4         ! 8                  
             func1%mhfx = mhfx_is_hfx
             func3%pxlda=0.75_real_8
             func3%pxgc=0.75_real_8
             func3%pclda=1._real_8
             func3%pcgc=1._real_8
             func3%phfx=0.25_real_8
          ELSEIF (keyword_contains(line,'HSE06')) THEN
             empvdwc%dft_func='HSE06'
             toldcode=.TRUE.
             cntl%tgc=.TRUE.
             cntl%tgcx=.TRUE.
             cntl%tgcc=.TRUE.
             cntl%thybrid=.TRUE.
             func1%mfxcx = mfxcx_is_slaterx
             func2%salpha=2._real_8/3._real_8
             func1%mfxcc=mfxcc_is_pz
             func1%mgcx = mgcx_is_pbex! Use regular PBEX
             func1%mgcc = mgcc_is_pbec! PBE-C
             func1%mhfx = mhfx_is_hfx
             func3%pxlda=1._real_8
             func3%pxgc=1._real_8
             func3%pclda=1._real_8
             func3%pcgc=1._real_8
             func3%phfx=0.25_real_8   ! Default alpha, adjust with keyword ACM1
             func1%msrx = msrx_is_erfc        ! Turn on erfc screened exchange
             func2%srxa=0.106_real_8  ! Default omega, adjust with keyword SCREENED EXCHANGE
             func1%mgcsrx = mgcsrx_is_hse      ! Take into account the short-range GGA-exchange
          ELSEIF (keyword_contains(line,'PBES')) THEN
             empvdwc%dft_func='PBES'
             cntl%tgc=.TRUE.
             cntl%tgcx=.TRUE.
             cntl%tgcc=.TRUE.
             func1%mfxcx = mfxcx_is_skipped
             func2%salpha=2._real_8/3._real_8
             func1%mfxcc= mfxcc_is_pade
             func1%mgcx = mgcx_is_pbesx
             func1%mgcc = mgcc_is_pbesc
             func1%mhfx = mhfx_is_skipped
          ELSEIF (keyword_contains(line,'PBE')) THEN
             empvdwc%dft_func='PBE'
             cntl%tgc=.TRUE.
             cntl%tgcx=.TRUE.
             cntl%tgcc=.TRUE.
             func1%mfxcx = mfxcx_is_skipped
             func2%salpha=2._real_8/3._real_8
             func1%mfxcc=mfxcc_is_pade
             func1%mgcx = mgcx_is_pbex
             func1%mgcc = mgcc_is_pbec
             func1%mhfx = mhfx_is_skipped
          ELSEIF (keyword_contains(line,'HCTH')) THEN
             empvdwc%dft_func='HCTH'
             toldcode=.TRUE.
             cntl%tgc=.TRUE.
             cntl%tgcx=.TRUE.
             cntl%tgcc=.TRUE.
             func1%mfxcx = mfxcx_is_skipped
             func1%mfxcc=mfxcc_is_skipped
             func1%mgcx = mgcx_is_hcth
             func1%mgcc = mgcc_is_hse
             func1%mhfx = mhfx_is_skipped
          ELSEIF (keyword_contains(line,'OPTX')) THEN
             toldcode=.TRUE.
             cntl%tgc=.TRUE.
             cntl%tgcx=.TRUE.
             cntl%tgcc=.TRUE.
             func1%mfxcx = mfxcx_is_skipped
             func1%mfxcc=mfxcc_is_skipped
             func1%mgcx = mgcx_is_optx
             func1%mgcc = mgcc_is_optc
             func1%mhfx = mhfx_is_skipped
          ELSEIF (keyword_contains(line,'OLYP')) THEN
             empvdwc%dft_func='OLYP'
             toldcode=.TRUE.
             cntl%tgc=.TRUE.
             cntl%tgcx=.TRUE.
             cntl%tgcc=.TRUE.
             func1%mfxcx = mfxcx_is_skipped
             func1%mfxcc=mfxcc_is_lyp
             func1%mgcx = mgcx_is_optx
             func1%mgcc = mgcc_is_lyp
             func1%mhfx = mhfx_is_skipped
          ELSEIF (keyword_contains(line,'OLDX3LYP')) THEN
             WRITE(output_unit,*) 'DFTIN| Using X3LYP as defined in CPMD 3.17 and prior.'
             WRITE(output_unit,*) '       The LDA part of LYP is used instead of VWN or PZ.'
             toldcode=.TRUE.
             cntl%tgc=.TRUE.
             cntl%tgcx=.TRUE.
             cntl%tgcc=.TRUE.
             cntl%thybrid=.TRUE.
             func1%mfxcx = mfxcx_is_slaterx
             func2%salpha=2._real_8/3._real_8
             func1%mfxcc=mfxcc_is_lyp ! this is not forced anymore
             func1%mgcx = mgcx_is_ox_hybrid
             func1%mgcc = mgcc_is_lyp
             func1%mhfx = mhfx_is_hfx
             func3%pxlda=0.782_real_8
             func3%pxgc=1._real_8
             func3%pclda=1._real_8 ! stay consistent with old versions
             func3%pcgc=0.871_real_8
             func3%phfx=0.218_real_8
          ELSEIF (keyword_contains(line,'MODX3LYP')) THEN
             WRITE(output_unit,*) 'DFTIN| Using X3LYP with 100% PZ instead of 19% VWN:'
             WRITE(output_unit,*) '       J. Phys. Chem. B 2006 110 (8) 3685/3691'
             toldcode=.TRUE.
             cntl%tgc=.TRUE.
             cntl%tgcx=.TRUE.
             cntl%tgcc=.TRUE.
             cntl%thybrid=.TRUE.
             func1%mfxcx = mfxcx_is_slaterx
             func2%salpha=2._real_8/3._real_8
             func1%mfxcc=mfxcc_is_pz ! PZ, not VWN
             func1%mgcx = mgcx_is_ox_hybrid
             func1%mgcc = mgcc_is_lyp
             func1%mhfx = mhfx_is_hfx
             func3%pxlda=0.782_real_8
             func3%pxgc=1._real_8
             func3%pclda=1._real_8
             func3%pcgc=0.871_real_8
             func3%phfx=0.218_real_8
          ELSEIF (keyword_contains(line,'X3LYP')) THEN
             WRITE(output_unit,*) 'DFTIN| Using the standard definition of X3LYP.'
             toldcode=.TRUE.
             cntl%tgc=.TRUE.
             cntl%tgcx=.TRUE.
             cntl%tgcc=.TRUE.
             cntl%thybrid=.TRUE.
             func1%mfxcx = mfxcx_is_slaterx
             func2%salpha=2._real_8/3._real_8
             func1%mfxcc=mfxcc_is_vwn ! VWN rather than PZ
             func1%mgcx = mgcx_is_ox_hybrid
             func1%mgcc = mgcc_is_lyp
             func1%mhfx = mhfx_is_hfx
             func3%pxlda=0.782_real_8
             func3%pxgc=1._real_8
             func3%pclda=0.19_real_8 ! this is not 1.0!
             func3%pcgc=0.871_real_8
             func3%phfx=0.218_real_8
          ELSEIF (keyword_contains(line,'OLDB3LYP')) THEN
             WRITE(output_unit,*) 'DFTIN| Using B3LYP as defined in CPMD 3.17 and prior.'
             WRITE(output_unit,*) '       The LDA part of LYP is used instead of VWN or PZ.'
             empvdwc%dft_func='B3LYP'
             toldcode=.TRUE.
             cntl%tgc=.TRUE.
             cntl%tgcx=.TRUE.
             cntl%tgcc=.TRUE.
             cntl%thybrid=.TRUE.
             func1%mfxcx = mfxcx_is_slaterx
             func2%salpha=2._real_8/3._real_8
             func1%mfxcc= mfxcc_is_lyp      ! this is to force LYP correlation. in versions <4, this was forced later on
             func1%mgcx = mgcx_is_becke88
             func1%mgcc = mgcc_is_lyp
             func1%mhfx = mhfx_is_hfx
             func3%pxlda=0.8_real_8
             func3%pxgc=0.72_real_8
             func3%pclda=1._real_8
             func3%pcgc=0.81_real_8
             func3%phfx=0.2_real_8
          ELSEIF (keyword_contains(line,'MODB3LYP')) THEN
             WRITE(output_unit,*) 'DFTIN| Using B3LYP with 100% PZ instead of 19% VWN:'
             WRITE(output_unit,*) '       J. Phys. Chem. B 2006 110 (8) 3685/3691'
             empvdwc%dft_func='B3LYP'
             toldcode=.TRUE.
             cntl%tgc=.TRUE.
             cntl%tgcx=.TRUE.
             cntl%tgcc=.TRUE.
             cntl%thybrid=.TRUE.
             func1%mfxcx = mfxcx_is_slaterx
             func2%salpha=2._real_8/3._real_8
             func1%mfxcc=mfxcc_is_pz ! PZ, not VWN
             func1%mgcx = mgcx_is_becke88
             func1%mgcc = mgcc_is_lyp
             func1%mhfx = mhfx_is_hfx
             func3%pxlda=0.8_real_8
             func3%pxgc=0.72_real_8
             func3%pclda=1._real_8
             func3%pcgc=0.81_real_8
             func3%phfx=0.2_real_8
          ELSEIF (keyword_contains(line,'B3LYP')) THEN
             WRITE(output_unit,*) 'DFTIN| Using the standard definition of B3LYP:'
             WRITE(output_unit,*) '       J. Phys. Chem. 1994 98 (45), 11623-11627'
             empvdwc%dft_func='B3LYP'
             toldcode=.TRUE.
             cntl%tgc=.TRUE.
             cntl%tgcx=.TRUE.
             cntl%tgcc=.TRUE.
             cntl%thybrid=.TRUE.
             func1%mfxcx = mfxcx_is_slaterx
             func2%salpha=2._real_8/3._real_8
             func1%mfxcc=mfxcc_is_vwn   ! vwn correlation as defined by Frisch and coworkers
             func1%mgcx = mgcx_is_becke88
             func1%mgcc = mgcc_is_lyp
             func1%mhfx = mhfx_is_hfx
             func3%pxlda=0.8_real_8
             func3%pxgc=0.72_real_8
             func3%pclda=0.19_real_8 ! this is not 1.0
             func3%pcgc=0.81_real_8
             func3%phfx=0.2_real_8
          ELSEIF (keyword_contains(line,'OLDB1LYP')) THEN
             WRITE(output_unit,*) 'DFTIN| Using B1LYP as defined in CPMD 3.17 and prior.'
             toldcode=.TRUE.
             cntl%tgc=.TRUE.
             cntl%tgcx=.TRUE.
             cntl%tgcc=.TRUE.
             cntl%thybrid=.TRUE.
             func1%mfxcx = mfxcx_is_slaterx
             func2%salpha=2._real_8/3._real_8
             func1%mfxcc=mfxcc_is_lyp ! this is not forced anymore
             func1%mgcx = mgcx_is_becke88
             func1%mgcc = mgcc_is_lyp
             func1%mhfx = mhfx_is_hfx
             func3%pxlda=0.75_real_8
             func3%pxgc=0.75_real_8
             func3%pclda=1._real_8
             func3%pcgc=1._real_8
             func3%phfx=0.25_real_8
          ELSEIF (keyword_contains(line,'MODB1LYP')) THEN
             WRITE(output_unit,*) 'DFTIN| Using B1LYP with 100% PZ+LYP instead of 100% VWN:'
             WRITE(output_unit,*) '       J. Phys. Chem. B 2006 110 (8) 3685/3691'
             toldcode=.TRUE.
             cntl%tgc=.TRUE.
             cntl%tgcx=.TRUE.
             cntl%tgcc=.TRUE.
             cntl%thybrid=.TRUE.
             func1%mfxcx = mfxcx_is_slaterx
             func2%salpha=2._real_8/3._real_8
             func1%mfxcc=mfxcc_is_pz ! PZ, not VWN
             func1%mgcx = mgcx_is_becke88
             func1%mgcc = mgcc_is_lyp
             func1%mhfx = mhfx_is_hfx
             func3%pxlda=0.75_real_8
             func3%pxgc=0.75_real_8
             func3%pclda=1._real_8
             func3%pcgc=1._real_8
             func3%phfx=0.25_real_8
          ELSEIF (keyword_contains(line,'B1LYP')) THEN
             WRITE(output_unit,*) 'DFTIN| Using the standard definition of B1LYP:'
             WRITE(output_unit,*) '       Chem. Phys. Lett. 1997, 274, 242-250 eq. (9)'
             toldcode=.TRUE.
             cntl%tgc=.TRUE.
             cntl%tgcx=.TRUE.
             cntl%tgcc=.TRUE.
             cntl%thybrid=.TRUE.
             func1%mfxcx = mfxcx_is_slaterx
             func2%salpha=2._real_8/3._real_8
             func1%mfxcc=mfxcc_is_lyp ! use LYP only, Chemical Physics Letters Volume 274, Issues 1-3, 1 August 1997, Pages 242-250
             func1%mgcx = mgcx_is_becke88
             func1%mgcc = mgcc_is_lyp
             func1%mhfx = mhfx_is_hfx
             func3%pxlda=0.75_real_8
             func3%pxgc=0.75_real_8
             func3%pclda=1._real_8
             func3%pcgc=1._real_8
             func3%phfx=0.25_real_8
          ELSEIF (keyword_contains(line,'TPSS')) THEN
             empvdwc%dft_func='TPSS'
             toldcode=.TRUE.
             cntl%tgc=.TRUE.
             cntl%tgcx=.TRUE.
             cntl%tgcc=.TRUE.
             cntl%ttau=.TRUE.
             func1%mfxcx = mfxcx_is_skipped
             func1%mfxcc=mfxcc_is_skipped
             func1%mgcx = mgcx_is_skipped
             func1%mgcc = mgcc_is_skipped
             func1%mtau = mtau_is_tpss
             func1%mhfx = mhfx_is_skipped
          ELSEIF (keyword_contains(line,'SAOP')) THEN
             lrf4%td_functional=10
             CALL tdm_fun(lrf4%td_functional,1)
             tenergy_ok=.FALSE.
          ELSEIF (keyword_contains(line,'LB94')) THEN
             lrf4%td_functional=11
             CALL tdm_fun(lrf4%td_functional,1)
             tenergy_ok=.FALSE.
          ELSEIF (keyword_contains(line,'GLLB')) THEN
             lrf4%td_functional=12
             CALL tdm_fun(lrf4%td_functional,1)
             tenergy_ok=.FALSE.
          ELSE
             CALL stopgm(procedureN,'Functional not available',& 
                  __LINE__,__FILE__)
          ENDIF

          !
          ! Does not apply to the new driver
          !
          IF (overwrite_pxgc ) func3%pxgc=new_pxgc
          IF (overwrite_pcgc ) func3%pcgc=new_pcgc
          IF (overwrite_pxlda) func3%pxlda=new_pxlda
          IF (overwrite_pclda) func3%pclda=new_pclda
   
       ENDIF

       !
       ! General overwriting
       !
       IF (overwrite_phfx ) THEN
          func3%phfx=new_phfx
       ENDIF

    END SUBROUTINE set_xc_functional
    ! ==--------------------------------------------------------------==
    SUBROUTINE set_xc_mts_low_func(line)
       !
       ! TODO: Once we get rid of the old xc code, move to cp_xc_utils


       CHARACTER(len=*), INTENT(in) :: line
       
       IF( cntl%use_xc_driver ) THEN

          cp_xc_mts_low_func_env%set_name  = 'MTS LOW FUNCTIONAL EXCHANGE-CORRELATION'
          cp_xc_mts_low_func_env%init = .TRUE.
          cp_xc_mts_low_func_env%polarized = cntl%tlsd
          CALL xstring(line,first,last)
          next_words = line(last+1:)
          i=1
          DO
             CALL xstring(next_words,first,last)
             IF (last-first < 0) EXIT
             IF (i > cp_xc_max_nbr_funcs) CALL stopgm(procedureN,' FUNCTIONAL: Request exceeds cp_xc_max_nbr_funcs, &
                                                                  &increase its value in cp_xc_utils',&
                                                      __LINE__,__FILE__)
             cp_xc_mts_low_func_env%funcs(i) = next_words(first:last)
             next_words = next_words(last+1:)
             i = i+1
          ENDDO
          cp_xc_mts_low_func_env%n_funcs = i-1
          !
          ! these are always 0 for xc_driver functionals
          !
          func1%mfxcc = mfxcc_is_skipped
          func1%mfxcx = mfxcx_is_skipped

          !
          ! Add HFX for hybrids
          !
          CALL cp_xc_mts_low_func_env%set_hybrid_defaults(default_divergence=cntl%div_analytical)

          !
          ! Postprocess the VDW string for Grimme
          !
          IF (cp_xc_mts_low_func_env%n_funcs == 1) THEN
             IF     (cp_xc_mts_low_func_env%funcs(1) == 'GGA_XC_BLYP') THEN
                empvdwc%dft_func='BLYP'
             ELSEIF (cp_xc_mts_low_func_env%funcs(1) == 'HYB_GGA_XC_B3LYP') THEN
                empvdwc%dft_func='B3LYP'
             ELSEIF (cp_xc_mts_low_func_env%funcs(1) == 'GGA_XC_PBE') THEN
                empvdwc%dft_func='PBE'
             ELSEIF (cp_xc_mts_low_func_env%funcs(1) == 'HYB_GGA_XC_PBE0') THEN
                empvdwc%dft_func='PBE0'
             ELSEIF (cp_xc_mts_low_func_env%funcs(1) == 'GGA_XC_BP' .OR. &
                     cp_xc_mts_low_func_env%funcs(1) == 'GGA_XC_BP86') THEN
                empvdwc%dft_func='BP86'
             ELSEIF (cp_xc_mts_low_func_env%funcs(1) == 'MGGA_XC_TPSS') THEN
                empvdwc%dft_func='TPSS'
             ELSEIF (cp_xc_mts_low_func_env%funcs(1) == 'GGA_XC_PBE_R' .OR. &
                     cp_xc_mts_low_func_env%funcs(1) == 'GGA_XC_REVPBE') THEN
                empvdwc%dft_func='REVPBE'
             ELSEIF (cp_xc_mts_low_func_env%funcs(1) == 'HYB_GGA_XC_PBE0_R' .OR. &
                     cp_xc_mts_low_func_env%funcs(1) == 'HYB_GGA_XC_REVPBE0') THEN
                empvdwc%dft_func='REVPBE0'
             ELSEIF (cp_xc_mts_low_func_env%funcs(1) == 'HYB_GGA_XC_CAM_B3LYP') THEN
                empvdwc%dft_func='CAM-B3LYP'
             ELSEIF (cp_xc_mts_low_func_env%funcs(1) == 'GGA_XC_HCTH_120') THEN
                empvdwc%dft_func='HCTH'
             ELSEIF (cp_xc_mts_low_func_env%funcs(1) == 'GGA_XC_OPBE') THEN
                empvdwc%dft_func='OPBE'
             ENDIF
          ELSEIF (cp_xc_mts_low_func_env%n_funcs == 2) THEN
             IF     ( any(cp_xc_mts_low_func_env%funcs(:) == 'GGA_X_B88') .AND. &
                      any(cp_xc_mts_low_func_env%funcs(:) == 'GGA_C_LYP') ) THEN
                empvdwc%dft_func='BLYP'
             ELSEIF ( any(cp_xc_mts_low_func_env%funcs(:) == 'GGA_X_PBE') .AND. &
                      any(cp_xc_mts_low_func_env%funcs(:) == 'GGA_C_PBE') ) THEN
                empvdwc%dft_func='PBE'
             ELSEIF ( any(cp_xc_mts_low_func_env%funcs(:) == 'GGA_X_B88') .AND. &
                      any(cp_xc_mts_low_func_env%funcs(:) == 'GGA_C_P86') ) THEN
                empvdwc%dft_func='BP86'
             ELSEIF ( any(cp_xc_mts_low_func_env%funcs(:) == 'MGGA_X_TPSS') .AND. &
                      any(cp_xc_mts_low_func_env%funcs(:) == 'MGGA_C_TPSS') ) THEN
                empvdwc%dft_func='TPSS'
             ELSEIF ( any(cp_xc_mts_low_func_env%funcs(:) == 'GGA_X_PBE_R') .AND. &
                      any(cp_xc_mts_low_func_env%funcs(:) == 'GGA_C_PBE') ) THEN
                empvdwc%dft_func='REVPBE'
             ELSEIF ( any(cp_xc_mts_low_func_env%funcs(:) == 'GGA_X_PW91') .AND. &
                      any(cp_xc_mts_low_func_env%funcs(:) == 'GGA_C_PW91') ) THEN
                empvdwc%dft_func='PW91'
             ELSEIF ( any(cp_xc_mts_low_func_env%funcs(:) == 'GGA_X_PBE_SOL') .AND. &
                      any(cp_xc_mts_low_func_env%funcs(:) == 'GGA_C_PBE_SOL') ) THEN
                empvdwc%dft_func='PBES'
             ELSEIF ( any(cp_xc_mts_low_func_env%funcs(:) == 'GGA_X_OPTX') .AND. &
                      any(cp_xc_mts_low_func_env%funcs(:) == 'GGA_C_LYP') ) THEN
                empvdwc%dft_func='OLYP'
             ELSEIF ( any(cp_xc_mts_low_func_env%funcs(:) == 'GGA_X_OPTX') .AND. &
                      any(cp_xc_mts_low_func_env%funcs(:) == 'GGA_C_PBE') ) THEN
                empvdwc%dft_func='OPBE'
             ENDIF
          ENDIF

       ELSE
          CALL stopgm(procedureN,' MTS: low functional has to be set up with XC_DRIVER',&
                                                      __LINE__,__FILE__)
       ENDIF

    END SUBROUTINE set_xc_mts_low_func
    ! ==--------------------------------------------------------------==
    SUBROUTINE set_xc_kernel(line)

       CHARACTER(len=*), INTENT(in) :: line
       
       IF( cntl%use_xc_driver ) THEN

          cp_xc_kernel_env%set_name  = 'LR EXCHANGE-CORRELATION KERNEL'
          cp_xc_kernel_env%init      = .TRUE.
          cp_xc_kernel_env%polarized = cntl%tlsd
          CALL xstring(line,first,last)
          next_words = line(last+1:)
          i=1
          DO
             CALL xstring(next_words,first,last)
             IF (last-first < 0) EXIT
             IF (i > cp_xc_max_nbr_funcs) CALL stopgm(procedureN,' FUNCTIONAL: Request exceeds cp_xc_kernel_max_nbr_funcs, &
                                                                  &increase its value in cp_xc_utils',&
                                                      __LINE__,__FILE__)
             cp_xc_kernel_env%funcs(i) = next_words(first:last)
             next_words = next_words(last+1:)
             i = i+1
          ENDDO
          cp_xc_kernel_env%n_funcs = i-1

          CALL cp_xc_kernel_env%set_hybrid_defaults(default_divergence=cntl%div_analytical)

       ELSE
          IF (keyword_contains(line,'NONE')) THEN
             lrf2%td_tgc=.FALSE.
             lrf2%td_tgcx=.FALSE.
             lrf2%td_tgcc=.FALSE.
             lrf1%td_x=0
             lrf1%td_c=0
             lrf1%td_gx=0
             lrf1%td_gc=0
             lrf1%td_hf=0
             lrf1%td_mtau=0
             lrf2%td_code=.FALSE.
          ELSEIF (keyword_contains(line,'PBE0')) THEN
             lrf2%td_code=.TRUE.
             lrf2%td_tgc=.TRUE.
             lrf2%td_tgcx=.TRUE.
             lrf2%td_tgcc=.TRUE.
             lrf2%td_hybrid=.TRUE.
             lrf1%td_x=1
             lrf1%td_c=2
             lrf1%td_gx=3
             lrf1%td_gc=4
             lrf1%td_hf=1
             lrf3%tdpxlda=0.75_real_8
             lrf3%tdpxgc=0.75_real_8
             lrf3%tdpclda=1._real_8
             lrf3%tdpcgc=1._real_8
             lrf3%tdphfx=0.25_real_8
          ELSEIF (keyword_contains(line,'PBE')) THEN
             lrf2%td_tgc=.TRUE.
             lrf2%td_tgcx=.TRUE.
             lrf2%td_tgcc=.TRUE.
             lrf1%td_x=0
             lrf1%td_c=9
             lrf1%td_gx=3
             lrf1%td_gc=4
             lrf1%td_hf=0
             lrf1%td_mtau=0
             lrf2%td_code=.FALSE.
          ELSEIF (keyword_contains(line,'LDA')) THEN
             lrf2%td_tgc=.FALSE.
             lrf2%td_tgcx=.FALSE.
             lrf2%td_tgcc=.FALSE.
             lrf1%td_x=0
             lrf1%td_c=9
             lrf1%td_gx=0
             lrf1%td_gc=0
             lrf1%td_hf=0
             lrf1%td_mtau=0
             lrf2%td_code=.FALSE.
          ELSE
             CALL stopgm(procedureN,'TDDFT kernel functional not available',& 
                  __LINE__,__FILE__)
          ENDIF
       ENDIF

    END SUBROUTINE set_xc_kernel
    ! ==--------------------------------------------------------------==
    SUBROUTINE check_options()
       ! 
       ! Test of options
       !
       IF (.not. cntl%use_xc_driver) THEN
          IF (cntl%thybrid) THEN
             toldcode=.TRUE.
             tabx%narray=0
             tabx%rmaxxc=1.0_real_8
          ENDIF
          ! 
          ! Defaulting to PADE is only needed for GRADIENT CORRECTION
          !
          IF (func1%mfxcx < 0) THEN
             func1%mfxcx = mfxcx_is_skipped
             IF (toldcode) func1%mfxcx = mfxcx_is_slaterx
          ENDIF
          IF (func1%mfxcc < 0) THEN
             func1%mfxcc=mfxcc_is_pade
             IF (toldcode) func1%mfxcc=mfxcc_is_pz
          ENDIF
          !  this was changed to avoid overwriting in case of .thybrid.
          ! (it is still overwriting if no hybrid is used, maybe this can be omitted completely?)
          IF (func1%mgcc == mgcc_is_lyp .AND. .NOT. cntl%thybrid) func1%mfxcc=mfxcc_is_lyp
          ! endthis
          IF (func1%mfxcc == mfxcc_is_pade.AND.func1%mfxcx /= mfxcx_is_skipped) CALL stopgm(procedureN,'FUNCTIONALS',& 
               __LINE__,__FILE__)
          IF (toldcode) THEN
             IF (func1%mfxcc == mfxcc_is_pade.AND.(func1%mgcc == mgcc_is_ggac .OR. func1%mgcc == mgcc_is_pbec)) THEN
                func1%mfxcc = mfxcc_is_pz
                func1%mfxcx = mfxcx_is_slaterx
             ENDIF
          ENDIF
          ! this is probably obsolete due to the overwrite statement above - now includes .not. thybrid
          IF (func1%mgcc == mgcc_is_lyp.AND.func1%mfxcc /= mfxcc_is_lyp .AND. .NOT. cntl%thybrid) THEN
             WRITE(output_unit,'(A)')&
                  ' WARNING: INCONSISTENT DEFINITION OF LYP FUNCTIONAL '
          ENDIF
          ! endthis
          IF ( func1%mfxcx == mfxcx_is_skipped .AND. func1%mfxcc == mfxcc_is_lyp .AND. func1%mgcc == mgcc_is_lyp .AND.&
               .NOT.toldcode ) THEN
             WRITE(output_unit,'(A)') ' IF YOU REALLY WANT TO USE THIS FUNCTIONAL'
             WRITE(output_unit,'(A)') ' (LYP WITHOUT SLATER EXCHANGE) USE THE'
             WRITE(output_unit,'(A)') ' OLDCODE OPTION'
             WRITE(output_unit,'(A)') ' BUT MOST LIKELY THIS IS AN INPUT ERROR'
             WRITE(output_unit,'(A)') ' SEE THE MANUAL ON HOW TO SPECIFY'
             WRITE(output_unit,'(A)') ' THE LYP FUNCTIONAL'
             CALL stopgm(procedureN,'FUNCTIONALS',& 
                  __LINE__,__FILE__)
          ENDIF
    
          IF (func1%mgcx == mgcx_is_dfr_zpbex.AND.func1%mgcc == mgcc_is_dfr_zpbec) THEN
             IF (func1%msrx /= msrx_is_ashcroft) CALL stopgm(procedureN,&
                  'wrong screening for the exact exchange',& 
                  __LINE__,__FILE__)
             IF (.NOT.cntl%thybrid.OR.func1%mhfx /= 1) CALL stopgm(procedureN,&
                  'exact exchange need to be turned on',& 
                  __LINE__,__FILE__)
          ENDIF
          IF (tabx%narray == 1) tabx%narray=0
          !
          ! cntl%tddft KERNEL
          IF (lrf1%td_x < 0 .OR. lrf1%td_c < 0 .OR. lrf1%td_gx < 0 .OR. lrf1%td_gc < 0) THEN
             lrf1%td_x = func1%mfxcx
             lrf1%td_c = func1%mfxcc
             lrf1%td_gx = func1%mgcx
             lrf1%td_gc = func1%mgcc
             lrf1%td_hf = func1%mhfx
             lrf2%td_tgc = cntl%tgc
             lrf2%td_tgcx = cntl%tgcx
             lrf2%td_tgcc = cntl%tgcc
             lrf2%td_code = toldcode
             lrf2%td_hybrid = cntl%thybrid
             lrf3%tdpxlda=func3%pxlda
             lrf3%tdpxgc=func3%pxgc
             lrf3%tdpclda=func3%pclda
             lrf3%tdpcgc=func3%pcgc
             lrf3%tdphfx=func3%phfx
          ENDIF
       ENDIF
       !
       ! HFX check
       IF (hfxc3%twscr.AND.( (.NOT. cntl%use_xc_driver .AND. func1%mhfx == 0 ) .OR. &
                             (cntl%use_xc_driver .AND. .NOT. cp_xc_functional_env%get_hfx) ) ) hfxc3%twscr=.FALSE.
       !
       ! force new hfx module even for thresholding
       IF(force_new_hfx) hfxc3%use_new_hfx=.TRUE.
       !
       ! force analytical or numerical computation of divergence term
       IF (force_div_analytical) THEN
          IF (force_div_numerical) CALL stopgm(procedureN,'Contradictory options (DIV_ANALYTICAL)',& 
               __LINE__,__FILE__)
          cntl%div_analytical = .TRUE.
       ENDIF
       IF (force_div_numerical) THEN
          IF (force_div_analytical) CALL stopgm(procedureN,'Contradictory options (DIV_NUMERICAL)',& 
               __LINE__,__FILE__)
          cntl%div_analytical = .FALSE.
       ENDIF
       !
       ! MTS check
       if (cntl%use_mts) then
          if (.not. cntl%use_xc_driver .and. need_dft) then
             call stopgm(procedureN,'MTS only comptatible with XC_DRIVER in DFT section',&
                __LINE__,__FILE__)
          end if

          if ( (mts%high_level=='DFT') .and. (.not. cp_xc_functional_env%init) ) then
             write(output_unit,*) 'ERROR: add FUNCTIONAL or MTS_HIGH_FUNC keyword to DFT section'
             call stopgm(procedureN,'MTS set up requires a functional for the high level forces',&
                __LINE__,__FILE__)
          end if

          if ( (mts%low_level=='DFT') .and. (.not. cp_xc_mts_low_func_env%init) ) then
             write(output_unit,*) 'ERROR: add MTS_LOW_FUNC keyword to DFT section'
             call stopgm(procedureN,'MTS set up requires a functional for the low level forces',&
                __LINE__,__FILE__)
          end if
       end if

    END SUBROUTINE check_options
    ! ==--------------------------------------------------------------==
    SUBROUTINE hfx_report()

       IF ( func1%mhfx /= 0 .OR. cp_xc_functional_env%get_hfx ) THEN
          WRITE(output_unit,'(A)') 
          IF (hfxc3%twscr) THEN
             IF (hfxc3%twfc) THEN
                WRITE(output_unit,'(A,T56,F10.3)')&
                     ' HFX SCREENING:  RADIUS FOR WANNIER CENTERS THRESHOLDING ',hfxc4%dwfc
                WRITE(output_unit,'(A,T65,L1)')&
                     ' HFX SCREENING:  KEEP RADIUS FIXED ',hfxc3%keep_radius_fixed
             ENDIF
             IF (hfxc3%twft) THEN
                WRITE(output_unit,'(A,T56,G10.3)')&
                     ' HFX SCREENING:  INTEGRAL THRESHOLD ',hfxc4%dwf_integral_thresh
                WRITE(output_unit,'(A,T56,I10)')&
                     ' HFX SCREENING:  RECOMPUTE INTEGRAL EVERY ',hfxc5%recomp_two_int_list_every
             ENDIF
             IF (hfxc3%tdiaw) WRITE(output_unit,'(A)') ' DIAGONAL SCREENING IN HFX '
          ENDIF
   
          WRITE(output_unit,'(A,T65,L1)') ' USE NEW HFX MODULE ',hfxc3%use_new_hfx
          IF (hfxc3%use_new_hfx) THEN
             SELECT CASE(hfxc5%hfx_distribution)
             CASE(hfx_dist_block_cyclic)
                WRITE(output_unit,'(A)') ' USE BLOCK CYCLIC DISTRIBUTION'
             CASE(hfx_dist_dynamic)
                WRITE(output_unit,'(A)') ' USE DYNAMIC DISTRIBUTION'
             END SELECT
             WRITE(output_unit,'(A,T62,I4)') ' BLOCK SIZE ',hfxc5%hfx_block_size
          ENDIF
          IF (cntl%use_scaled_hfx) THEN
             WRITE(output_unit,'(1X,A)') 'USE COORDINATE-SCALED EXACT EXCHANGE (ScEX) FOR ISOLATED SYSTEMS'
             WRITE(output_unit,'(1X,A)') 'AND THE TUCKERMAN-MARTYNA POISSON SOLVER. PLEASE CITE:'
             WRITE(output_unit,'(1X,A)') '                               M.P. Bircher and U. Rothlisberger'
             WRITE(output_unit,'(1X,A)') '                J. Phys. Chem. Lett., 2018, 9 (14), pp 3886–3890'
             WRITE(output_unit,'(1X,A)') '                                DOI: 10.1021/acs.jpclett.8b01620'
          ENDIF
       ENDIF

    END SUBROUTINE hfx_report
    ! ==--------------------------------------------------------------==
    SUBROUTINE newcode_oldcode_report()

       WRITE(output_unit,'(/,1x,64("!"))')
       WRITE(output_unit,'(1x,A)') 'WARNING!'
       WRITE(output_unit,'(7x,A)') 'NEWCODE AND OLDCODE WILL BE DEPRECATED IN A FUTURE RELEASE'
       WRITE(output_unit,'(7x,A)') 'PLEASE SWITCH TO THE NEW XC DRIVER AND/OR LIBXC'
       WRITE(output_unit,'(7x,A)') '(KEYWORD: XC_DRIVER IN &DFT)'
       WRITE(output_unit,'(1x,64("!"))')
       WRITE(output_unit,*)
       WRITE(output_unit,'(A)') ' EXCHANGE CORRELATION FUNCTIONALS '
       IF (lrf4%td_functional /= 0) THEN
          WRITE(output_unit,'(A)') ' EXCHANGE CORRELATION POTENTIAL METHOD'
          WRITE(output_unit,'(A)')&
               ' WARNING: NO CORRESPONDING ENERGY FUNCTIONAL AVAILABLE'
          IF (lrf4%td_functional == 10) THEN
             WRITE(output_unit,'(A)') ' SAOP POTENTIAL (modified)'
          ELSEIF (lrf4%td_functional == 11) THEN
             WRITE(output_unit,'(A)') ' LB94 POTENTIAL (modified)'
          ELSEIF (lrf4%td_functional == 12) THEN
             WRITE(output_unit,'(A)') ' GLLB POTENTIAL (modified)'
          ENDIF
       ELSEIF (cntl%thybrid) THEN
          WRITE(output_unit,'(A)') ' HYBRID FUNCTIONAL '
          IF (func1%mfxcx == mfxcx_is_skipped) THEN
             WRITE(output_unit,'(A,T62,A4)') '    LDA EXCHANGE:','NONE'
          ELSEIF (func1%mfxcx == mfxcx_is_slaterx) THEN
             WRITE(output_unit,'(A,T40,A,T58,F8.2)')&
                  '    LDA EXCHANGE:','SLATER',func3%pxlda
          ENDIF
          IF (func1%mfxcc == mfxcc_is_skipped) THEN
             WRITE(output_unit,'(A,T62,A4)') '    LDA CORRELATION:','NONE'
          ELSEIF (func1%mfxcc == mfxcc_is_pz) THEN
             WRITE(output_unit,'(A,T40,A,T58,F8.2)')&
                  '    LDA CORRELATION:','PERDEW & ZUNGER',func3%pclda
          ELSEIF (func1%mfxcc == mfxcc_is_vwn) THEN
             WRITE(output_unit,'(A,T40,A,T58,F8.2)')&
                  '    LDA CORRELATION:','VWN',func3%pclda
             IF (cntl%tlsd) THEN
                CALL stopgm(procedureN,'LSD_VWN N/A, you may want to switch to USE_XC_DRIVER',& 
                     __LINE__,__FILE__)
             ENDIF
          ELSEIF (func1%mfxcc == mfxcc_is_lyp) THEN
             WRITE(output_unit,'(A,T40,A,T58,F8.2)')&
                  '    LDA CORRELATION:','LYP',func3%pclda
          ELSEIF (func1%mfxcc == mfxcc_is_pw) THEN
             WRITE(output_unit,'(A,T40,A,T58,F8.2)')&
                  '    LDA CORRELATION:','PERDEW & WANG',func3%pclda
             IF (cntl%tlsd) THEN
                CALL stopgm(procedureN,'LSD_PW N/A, you may want to switch to USE_XC_DRIVER',& 
                     __LINE__,__FILE__)
             ENDIF
          ELSEIF (func1%mfxcc == mfxcc_is_pade) THEN
             WRITE(output_unit,'(A,T40,A,T58,F8.2)')&
                  '    LDA CORRELATION:','PADE',func3%pclda
          ENDIF
          ! this part was changed to report the LDA part of LYP separately (including the weirdest possible combinations)
          IF (func1%mgcc==2)  WRITE(output_unit,'(A,T40,A,T58,F8.2)')&
               '                    ','LYP',func3%pcgc
          ! endthis
          IF (func1%mgcx == mgcx_is_becke88) THEN
             WRITE(output_unit,'(A,T40,A,T58,F8.2)')&
                  '    GC EXCHANGE:','BECKE88',func3%pxgc
          ELSEIF (func1%mgcx == mgcx_is_ggax) THEN
             WRITE(output_unit,'(A,T40,A,T58,F8.2)')&
                  '    GC EXCHANGE:','PERDEW-WANG',func3%pxgc
          ELSEIF (func1%mgcx == mgcx_is_pbex) THEN
             WRITE(output_unit,'(A,T40,A,T58,F8.2)')&
                  '    GC EXCHANGE:','PBE',func3%pxgc
          ELSEIF (func1%mgcx == mgcx_is_revpbex) THEN
             WRITE(output_unit,'(A,T40,A,T58,F8.2)')&
                  '    GC EXCHANGE:','revPBE',func3%pxgc
          ELSEIF (func1%mgcx == mgcx_is_hcth) THEN
             WRITE(output_unit,'(A,T40,A,T58,F8.2)')&
                  '    GC EXCHANGE:','HCTH',func3%pxgc
          ELSEIF (func1%mgcx == mgcx_is_optx) THEN
             WRITE(output_unit,'(A,T40,A,T58,F8.2)')&
                  '    GC EXCHANGE:','OPTX',func3%pxgc
          ELSEIF (func1%mgcx == mgcx_is_ox) THEN
             WRITE(output_unit,'(A,T40,A,T58,F8.2)')&
                  '    GC EXCHANGE:','X(0.722 B88 + 0.347 PW91)',func3%pxgc
          ELSEIF (func1%mgcx == mgcx_is_ox_hybrid) THEN
             WRITE(output_unit,'(A,T40,A,T58,F8.2)')&
                  '    GC EXCHANGE:','X(0.542 B88 + 0.167 PW91)',func3%pxgc
          ELSEIF (func1%mgcx == mgcx_is_dfr_zpbex) THEN
             WRITE(output_unit,'(A,T40,A,T58,F8.2)')&
                  '       EXCHANGE:','PBE',func3%pxgc
             WRITE(output_unit,'(A,T40,A,T58,F8.2)')&
                  '                ','ZPBE',func3%phfx
          ELSEIF (func1%mgcx == mgcx_is_dfr_xpbex) THEN
             WRITE(output_unit,'(A,T40,A,T58,F8.2)')&
                  '       EXCHANGE:','PBE',func3%pxgc
          ELSEIF (func1%mgcx == mgcx_is_dfr_xpbex_hybrid) THEN
             WRITE(output_unit,'(A,T40,A,T58,F8.2)')&
                  '       EXCHANGE:','PBE',func3%pxgc
             WRITE(output_unit,'(A,T40,A,T58,F8.2)')&
                  '                ','XPBE',func3%phfx
          ENDIF
          IF (func1%mgcc == mgcc_is_perdew86) THEN
             WRITE(output_unit,'(A,T40,A,T58,F8.2)')&
                  '    GC CORRELATION:','PERDEW86',func3%pcgc
          ELSEIF (func1%mgcc == mgcc_is_lyp) THEN
             WRITE(output_unit,'(A,T40,A,T58,F8.2)')&
                  '    GC CORRELATION:','LYP',func3%pcgc
          ELSEIF (func1%mgcc == mgcc_is_ggac) THEN
             WRITE(output_unit,'(A,T40,A,T58,F8.2)')&
                  '    GC CORRELATION:','PERDEW-WANG',func3%pcgc
          ELSEIF (func1%mgcc == mgcc_is_pbec) THEN
             WRITE(output_unit,'(A,T40,A,T58,F8.2)')&
                  '    GC CORRELATION:','PBE',func3%pcgc
          ELSEIF (func1%mgcc == mgcc_is_hse) THEN
             WRITE(output_unit,'(A,T40,A,T58,F8.2)')&
                  '    GC EXCHANGE:','HCTH',func3%pcgc
          ELSEIF (func1%mgcc == mgcc_is_optc) THEN
             WRITE(output_unit,'(A,T40,A,T58,F8.2)')&
                  '    GC CORRELATION:','PBE1W',func3%pcgc
          ELSEIF (func1%mgcc == mgcc_is_pbesc) THEN
             WRITE(output_unit,'(A,T40,A,T58,F8.2)')&
                  '    GC CORRELATION:','PBES',func3%pcgc
          ELSEIF (func1%mgcc == mgcc_is_dfr_zpbec) THEN
             WRITE(output_unit,'(A,T40,A,T58,F8.2)')&
                  '       CORRELATION:','PBE',func3%pcgc
          ENDIF
          IF (func1%mhfx == 1) THEN
             WRITE(output_unit,'(A,T58,F8.2)') '    HARTREE-FOCK EXCHANGE:',func3%phfx
             IF (func1%mfxcx == mfxcx_is_slaterx.AND.func1%mfxcc == mfxcc_is_pz.AND.func1%mgcx == mgcx_is_pbex.AND.&
                  func1%mgcc == mgcc_is_pbec.AND.func1%msrx == msrx_is_erfc.AND.&
                  (ABS(func2%srxa-0.106_real_8)) < 1.e-10_real_8) THEN
                WRITE(output_unit,'(A,T58)') '    HYBRID HSE06 FUNCTIONAL'
             ELSEIF(func1%mfxcx == mfxcx_is_slaterx.AND.func1%mfxcc == mfxcc_is_vwn.AND.func1%mgcx == mgcx_is_pbex.AND.&
                  func1%mgcc == mgcc_is_pbec) THEN
                WRITE(output_unit,'(A,T58)') '    HYBRID PBE0 FUNCTIONAL'
             ENDIF
             IF (func1%mgcx == mgcx_is_dfr_zpbex.AND.func1%mgcc == mgcc_is_dfr_zpbec) THEN
                WRITE(output_unit,'(A,T58)') '    HYBRID ZPBE0 FUNCTIONAL'
             ENDIF
             IF (func1%mgcx == mgcx_is_dfr_xpbex_hybrid.AND.func1%mgcc == mgcc_is_dfr_zpbec) THEN
                WRITE(output_unit,'(A,T58)') '    HYBRID XPBE0 FUNCTIONAL'
             ENDIF
          ELSEIF (func1%mhfx == 2) THEN
             WRITE(output_unit,'(A,T58,F8.2)') '    HARTREE CORRECTION:',func3%phfx
          ENDIF
       ELSEIF (func1%mhfx == 0) THEN
          IF (func1%mfxcx == mfxcx_is_skipped) THEN
             WRITE(output_unit,'(A,T62,A4)') '    LDA EXCHANGE:','NONE'
          ELSEIF (func1%mfxcx == mfxcx_is_slaterx) THEN
             WRITE(output_unit,'(A,T42,A15,F8.5,A1)')&
                  '    LDA EXCHANGE:','SLATER (ALPHA =',func2%salpha,')'
          ENDIF
          IF (func1%mfxcc == mfxcc_is_skipped) THEN
             WRITE(output_unit,'(A,T62,A4)') '    LDA CORRELATION:','NONE'
          ELSEIF (func1%mfxcc == mfxcc_is_pz) THEN
             WRITE(output_unit,'(A,T51,A15,/,A)')&
                  '    LDA CORRELATION:','PERDEW & ZUNGER',&
                  '       [J.P. PERDEW AND A ZUNGER, PRB 235048 (1981)]'
          ELSEIF (func1%mfxcc == mfxcc_is_vwn) THEN
             WRITE(output_unit,'(A,T46,A20,/,2A)')&
                  '    LDA CORRELATION:','VOSKO, WILK & NUSAIR',&
                  '       [S.H. VOSKO, L. WILK, AND M. NUSAIR,',&
                  ' CAN. J. PHYS. 581200 (1980)]'
             IF (cntl%tlsd) THEN
                CALL stopgm(procedureN,'LSD_VWN N/A, you may want to switch to USE_XC_DRIVER',& 
                     __LINE__,__FILE__)
             ENDIF
          ELSEIF (func1%mfxcc == mfxcc_is_lyp) THEN
             WRITE(output_unit,'(A,T50,A16,/,2A)')&
                  '    LDA CORRELATION:','LEE, YANG & PARR',&
                  '       [C.L. LEE, W. YANG, AND R.G. PARR, PRB 37785',&
                  ' (1988)]'
          ELSEIF (func1%mfxcc == mfxcc_is_pw) THEN
             WRITE(output_unit,'(A,T53,A13,/,A)')&
                  '    LDA CORRELATION:','PERDEW & WANG',&
                  '       [J.P. PERDEW AND Y. WANG, PRB 4513244 (1992)]'
             IF (cntl%tlsd) THEN
                CALL stopgm(procedureN,'LSD_PW N/A, you may want to switch to USE_XC_DRIVER',& 
                     __LINE__,__FILE__)
             ENDIF
          ELSEIF (func1%mfxcc == mfxcc_is_pade) THEN
             WRITE(output_unit,'(A,/,A)') '    LDA XC THROUGH PADE APPROXIMATION',&
                  '    S.GOEDECKER, J.HUTTER, M.TETER PRB 541703 (1996)'
             !
             ! mfxcc is neither coded nor set - commented it out
             !
             ! ELSEIF (func1%mfxcc == 10) THEN
             !    WRITE(output_unit,'(A)') ' +++++++++++++ WARNING LYP+VWN'
             !    IF (cntl%tlsd) THEN
             !       CALL stopgm('LSD_VWN','N/A, you may want to switch to USE_XC_DRIVER',& 
             !            __LINE__,__FILE__)
             !    ENDIF
          ENDIF
          IF (cntl%tgcc .OR. cntl%tgcx) THEN
             WRITE(output_unit,'(A,/,A,T51,1PE15.5)')&
                  '    GRADIENT CORRECTED FUNCTIONAL',&
                  '    DENSITY THRESHOLD: ',cntr%gceps
          ENDIF
          IF (cntl%tgcx) THEN
             IF (func1%mgcx == mgcx_is_becke88) THEN
                WRITE(output_unit,'(A,/,A,/,A,T56,F10.6)')&
                     '    EXCHANGE ENERGY',&
                     '       [A.D. BECKE, PHYS. REV. A 38, 3098 (1988)]',&
                     '       PARAMETER BETA:',func2%bbeta
             ELSEIF (func1%mgcx == mgcx_is_ggax) THEN
                WRITE(output_unit,'(A,/,A,A)') '    EXCHANGE ENERGY',&
                     '       [GGA: J.P. PERDEW ET AL. PHYS. REV. B 46, 6671',&
                     ' (1992)]'
             ELSEIF (func1%mgcx == mgcx_is_pbex) THEN
                WRITE(output_unit,'(A,/,A)') '    EXCHANGE ENERGY',&
                     '       [PBE: J.P. PERDEW ET AL. PRL 77, 3865 (1996)]'
             ELSEIF (func1%mgcx == mgcx_is_revpbex) THEN
                WRITE(output_unit,'(A,A)') '    EXCHANGE ENERGY    ',&
                     ' [revPBE: Y. ZHANG ET AL. PRL 80, 890 (1998)]'
             ELSEIF (func1%mgcx == mgcx_is_hcth) THEN
                WRITE(output_unit,'(A,A)') '    HCTH/120 XC FUNCTIONAL ',&
                     ' [HCTH: N.C. HANDY ET AL. JCP 109, 6264 (1998)]'
             ELSEIF (func1%mgcx == mgcx_is_optx) THEN
                WRITE(output_unit,'(A,A)') '    OPTX XC FUNCTIONAL ',&
                     ' [OPTX: N.C. HANDY ET AL. JCP 116, 5411 (2002)]'
             ELSEIF (func1%mgcx == mgcx_is_ox .OR. func1%mgcx == mgcx_is_ox_hybrid) THEN
                WRITE(output_unit,'(A,A)') '    EXCHANGE ENERGY    ',&
                     ' [XIN XU and GODDARD PNAS 101, 2673 (2004)]'
             ELSEIF (func1%mgcx == mgcx_is_pbesx) THEN
                WRITE(output_unit,'(A,A)') '    PBEsol  XC FUNCTIONAL ',&
                     ' [PBEsol: J.P. PERDEW et al (2007) ]'
             ELSEIF (func1%mgcx == mgcx_is_dfr_zpbex) THEN
                WRITE(output_unit,'(A,/,A)') '    EXCHANGE ENERGY',&
                     '       [PBE: J.P. PERDEW ET AL. PRL 77, 3865 (1996)]'
                WRITE(output_unit,'(A,/,A)') '                   ',&
                     '       [ZPBE: ... on its way]'
             ELSEIF (func1%mgcx == mgcx_is_dfr_xpbex_hybrid) THEN
                WRITE(output_unit,'(A,/,A)') '    EXCHANGE ENERGY',&
                     '       [PBE: J.P. PERDEW ET AL. PRL 77, 3865 (1996)]'
                WRITE(output_unit,'(A,/,A)') '                   ',&
                     '       [XPBE: ... on its way]'
             ENDIF
          ENDIF
          IF (cntl%tgcc) THEN
             IF (func1%mgcc == mgcc_is_perdew86) THEN
                WRITE(output_unit,'(A,/,A)') '    CORRELATION ENERGY',&
                     '       [J.P. PERDEW, PHYS. REV. B 33, 8822 (1986)]'
             ELSEIF (func1%mgcc == mgcc_is_lyp) THEN
                WRITE(output_unit,'(A,/,A)') '    CORRELATION ENERGY',&
                     '       [LYP: C.L. LEE ET AL. PHYS. REV. B 37, 785 (1988)]'
             ELSEIF (func1%mgcc == mgcc_is_ggac) THEN
                WRITE(output_unit,'(A,/,A)') '    CORRELATION ENERGY ',&
                     '       [GGA: J.P. PERDEW ET AL. PHYS. REV. B 46, 6671',&
                     ' (1992)]'
                IF (cntl%tlsd) THEN
                   CALL stopgm(procedureN,'LSD_GGAC NYI, you may want to link to libxc',& 
                        __LINE__,__FILE__)
                ENDIF
             ELSEIF (func1%mgcc == mgcc_is_pbec) THEN
                WRITE(output_unit,'(A,/,A)') '    CORRELATION ENERGY ',&
                     '       [PBE: J.P. PERDEW ET AL. PRL 77, 3865 (1996)]'
             ELSEIF (func1%mgcc == mgcc_is_hse) THEN
                WRITE(output_unit,'(A,A)') '                 '
             ELSEIF (func1%mgcc == mgcc_is_pbesc) THEN
                WRITE(output_unit,'(A,A)') '    PBEsol  XC FUNCTIONAL ',&
                     ' [PBEsol: J.P. PERDEW et al (2007) ]'
             ELSEIF (func1%mgcc == mgcc_is_dfr_zpbec) THEN
                WRITE(output_unit,'(A,/,A)') '    CORRELATION ENERGY ',&
                     '       [GGA: J.P. PERDEW ET AL. PHYS. REV. B 46, 6671',&
                     ' (1992)]'
             ENDIF
          ENDIF
          IF (cntl%ttau) THEN
             IF (func1%mtau == mtau_is_tpss) THEN
                WRITE(output_unit,'(A)') '    META FUNCTIONAL TPSS'
             ENDIF
          ENDIF
       ELSEIF (func1%mhfx == mhfx_is_hartree) THEN
          WRITE(output_unit,'(T37,A)') 'ORTHOGONALIZED HARTREE METHOD'
       ELSEIF (func1%mhfx == mhfx_is_hfx) THEN
          IF (MAX(lrf1%td_x,lrf1%td_c,lrf1%td_gx,lrf1%td_gc) == 0) THEN
             WRITE(output_unit,'(T47,A)') 'HARTREE-FOCK METHOD'
          ENDIF
       ENDIF
       IF (func1%msrx /= msrx_is_skipped.AND.func1%mhfx /= 0) THEN
          IF (func1%msrx == msrx_is_exp) THEN
             WRITE(output_unit,'(T25,A)')&
                  'SCREENED EXACT EXCHANGE USING EXP(-A*R)/R'
          ELSEIF (func1%msrx == msrx_is_erfc) THEN
             WRITE(output_unit,'(T25,A)')&
                  'SCREENED EXACT EXCHANGE USING ERFC(A*R)/R'
          ELSEIF (func1%msrx == msrx_is_ashcroft) THEN
             WRITE(output_unit,'(T25,A)')&
                  'SCREENED EXACT EXCHANGE USING ASHCROFT'
          ELSEIF (func1%msrx == msrx_is_CAM) THEN
             WRITE(output_unit,'(T25,A)')&
                  'SCREENED EXACT EXCHANGE FOR CAM-B3LYP'
          ELSE
             CALL stopgm(procedureN,'Unknown tag for range separated HFX',& 
                  __LINE__,__FILE__)
          ENDIF
          IF (func1%msrx == msrx_is_ashcroft) THEN
             WRITE(output_unit,'(T25,A,T57,F9.3)')&
                  'SCREENING PARAMETER [BOHR]:',ashcroft_coulomb_rcuT
          ELSE
             WRITE(output_unit,'(T24,A,T57,F9.3)')&
                  ' SCREENING PARAMETER [BOHR^-1]:',func2%srxa
          ENDIF
       ENDIF
       IF (func1%mgcsrx /= mgcsrx_is_skipped) THEN
          WRITE(output_unit,'(T24,A,T57,F9.3)')' USING COMPENSATING PBE-SR'
       ENDIF

    END SUBROUTINE newcode_oldcode_report
    ! ==--------------------------------------------------------------==
    SUBROUTINE tddft_newcode_oldcode_report()

       WRITE(output_unit,*)
       IF (lrf1%td_x /= func1%mfxcx .OR. lrf1%td_c /= func1%mfxcc .OR.&
            lrf1%td_gx /= func1%mgcx .OR. lrf1%td_gc /= func1%mgcc .OR.lrf1%td_hf /= func1%mhfx) THEN
          WRITE(output_unit,'(A)') ' FUNCTIONAL FOR XC KERNEL IN LINEAR RESPONSE'
          IF (MAX(lrf1%td_x,lrf1%td_c,lrf1%td_gx,lrf1%td_gc,lrf1%td_hf) == 0) THEN
             WRITE(output_unit,'(T60,A)') 'NONE'
          ENDIF
          IF (lrf1%td_c == 9) THEN
             WRITE(output_unit,'(A,/,A)') '    LDA XC THROUGH PADE APPROXIMATION',&
                  '    S.GOEDECKER, J.HUTTER, M.TETER PRB 541703 (1996)'
          ENDIF
          IF (lrf1%td_gx == 3) THEN
             WRITE(output_unit,'(A,/,A)') '    EXCHANGE ENERGY',&
                  '       [PBE: J.P. PERDEW ET AL. PRL 77, 3865 (1996)]'
          ENDIF
          IF (lrf1%td_gc == 4) THEN
             WRITE(output_unit,'(A,/,A)') '    CORRELATION ENERGY',&
                  '       [PBE: J.P. PERDEW ET AL. PRL 77, 3865 (1996)]'
          ENDIF
          WRITE(output_unit,*)
       ENDIF
       IF (treff) THEN
          IF (tr_x == 0.AND.tr_c == 9.AND.tr_gx == 0.AND.tr_gc == 0) THEN
             WRITE(output_unit,'(A,T62,A,/)') ' REFERENCE FUNCTIONAL: ',' LDA'
             ELSEIF&
                  (tr_x == 0.AND.tr_c == 9.AND.tr_gx == 3.AND.tr_gc == 4) THEN
             WRITE(output_unit,'(A,T62,A,/)') ' REFERENCE FUNCTIONAL: ',' PBE'
             ELSEIF&
                  (tr_x == 1.AND.tr_c == 3.AND.tr_gx == 1.AND.tr_gc == 2) THEN
             WRITE(output_unit,'(A,T62,A,/)') ' REFERENCE FUNCTIONAL: ',' BLYP'
             ELSEIF&
                  (tr_x == 0.AND.tr_c == 3.AND.tr_gx == 6.AND.tr_gc == 2) THEN
             WRITE(output_unit,'(A,T62,A,/)') ' REFERENCE FUNCTIONAL: ',' OLYP'
             ELSEIF&
                  (tr_x == 0.AND.tr_c == 9.AND.tr_gx == 1.AND.tr_gc == 1) THEN
             WRITE(output_unit,'(A,T62,A,/)') ' REFERENCE FUNCTIONAL: ',' BP'
          ELSE
             CALL stopgm(procedureN,'Reference functional n/a',& 
                  __LINE__,__FILE__)
          ENDIF
       ENDIF

    END SUBROUTINE tddft_newcode_oldcode_report
    ! ==--------------------------------------------------------------==
    SUBROUTINE broadcast_dftin()

       ! ..CNTL
       CALL mp_bcast_byte(cntl,size_in_bytes_of(cntl),parai%io_source,parai%cp_grp)
       ! ..CNTI
       CALL mp_bcast_byte(cnti,size_in_bytes_of(cnti),parai%io_source,parai%cp_grp)
       ! ..CNTR
       CALL mp_bcast_byte(cntr,size_in_bytes_of(cntr),parai%io_source,parai%cp_grp)
       ! ..TABX
       CALL mp_bcast(toldcode,parai%io_source,parai%cp_grp)
       CALL mp_bcast_byte(tabx, size_in_bytes_of(tabx),parai%io_source,parai%cp_grp)
       CALL mp_bcast(tabx%narray,parai%io_source,parai%cp_grp)
       ! ..FUNC
       CALL mp_bcast_byte(func2, size_in_bytes_of(func2),parai%io_source,parai%cp_grp)
       CALL mp_bcast_byte(func1, size_in_bytes_of(func1),parai%io_source,parai%cp_grp)
       CALL mp_bcast_byte(func3, size_in_bytes_of(func3),parai%io_source,parai%cp_grp)
       CALL mp_bcast(ashcroft_coulomb_rcut,parai%io_source,parai%cp_grp)
       ! ..LINRES
       CALL mp_bcast_byte(lrf1, size_in_bytes_of(lrf1),parai%io_source,parai%cp_grp)
       CALL mp_bcast_byte(lrf2, size_in_bytes_of(lrf2),parai%io_source,parai%cp_grp)
       CALL mp_bcast_byte(lrf3, size_in_bytes_of(lrf3),parai%io_source,parai%cp_grp)
       CALL mp_bcast(lrf4%td_functional,parai%io_source,parai%cp_grp)
       ! ..ENERGY
       CALL mp_bcast(tenergy_ok,parai%io_source,parai%cp_grp)
       ! ..HFX 3/4
 
       CALL mp_bcast_byte(hfxc3, size_in_bytes_of(hfxc3),parai%io_source,parai%cp_grp)
       CALL mp_bcast_byte(hfxc5, size_in_bytes_of(hfxc5),parai%io_source,parai%cp_grp)
       CALL mp_bcast_byte(hfxc4, size_in_bytes_of(hfxc4),parai%io_source,parai%cp_grp)
 
       IF (hfxc3%twscr) THEN
          cntl%tdipd=.TRUE.
          wannl%twann=.TRUE.
       ENDIF
 
       ! xc_driver
       CALL mp_bcast_byte(cp_xc_functional_env, size_in_bytes_of(cp_xc_functional_env),&
            parai%io_source,parai%cp_grp)
       CALL mp_bcast_byte(cp_xc_kernel_env, size_in_bytes_of(cp_xc_kernel_env),&
            parai%io_source,parai%cp_grp)
       CALL mp_bcast_byte(cp_xc_mts_low_func_env, size_in_bytes_of(cp_xc_mts_low_func_env),&
            parai%io_source,parai%cp_grp)
       CALL mp_bcast_byte(cp_gga_x_param, size_in_bytes_of(cp_gga_x_param),parai%io_source,parai%cp_grp)
       CALL mp_bcast_byte(cp_gga_c_param, size_in_bytes_of(cp_gga_c_param),parai%io_source,parai%cp_grp)

    END SUBROUTINE broadcast_dftin
    ! ==--------------------------------------------------------------==
  END SUBROUTINE dftin
  ! ==================================================================
  SUBROUTINE tdm_fun(fun,itag)
    INTEGER                                  :: fun, itag

    CHARACTER(*), PARAMETER                  :: procedureN = 'tdm_fun'

! ==--------------------------------------------------------------==

    IF (itag > 0) THEN
       toldcode=.FALSE.
       cntl%tpotential=.FALSE.
       func2%salpha=2._real_8/3._real_8
       IF (fun == 0) THEN
          cntl%tgc=.FALSE.
          cntl%tgcx=.FALSE.
          cntl%tgcc=.FALSE.
          func1%mfxcx = mfxcx_is_skipped
          func1%mfxcc=mfxcc_is_skipped
          func1%mgcx = mgcx_is_skipped
          func1%mgcc = mgcc_is_skipped
       ELSEIF (fun == 1) THEN
          cntl%tgc=.FALSE.
          cntl%tgcx=.FALSE.
          cntl%tgcc=.FALSE.
          func1%mfxcx = mfxcx_is_skipped
          func1%mfxcc=mfxcc_is_pade
          func1%mgcx = mgcx_is_skipped
          func1%mgcc = mgcc_is_skipped
          func1%mhfx = mhfx_is_skipped
       ELSEIF (fun == 2) THEN
          cntl%tgc=.TRUE.
          cntl%tgcx=.TRUE.
          cntl%tgcc=.TRUE.
          func1%mfxcx = mfxcx_is_skipped
          func1%mfxcc=mfxcc_is_pade
          func1%mgcx = mgcx_is_becke88
          func1%mgcc = mgcc_is_perdew86
          func1%mhfx = mhfx_is_skipped
       ELSEIF (fun == 3) THEN
          cntl%tgc=.TRUE.
          cntl%tgcx=.TRUE.
          cntl%tgcc=.TRUE.
          func1%mfxcx = mfxcx_is_skipped
          func1%mfxcc=mfxcc_is_pade
          func1%mgcx = mgcx_is_pbex
          func1%mgcc = mgcc_is_pbec
          func1%mhfx = mhfx_is_skipped
       ELSEIF (fun == 4) THEN
          cntl%tgc=.TRUE.
          cntl%tgcx=.TRUE.
          cntl%tgcc=.TRUE.
          func1%mfxcx = mfxcx_is_slaterx
          func1%mfxcc=mfxcc_is_lyp
          func1%mgcx = mgcx_is_becke88
          func1%mgcc = mgcc_is_lyp
          func1%mhfx = mhfx_is_skipped
       ELSEIF (fun == 5) THEN
          toldcode=.TRUE.
          cntl%tgc=.TRUE.
          cntl%tgcx=.TRUE.
          cntl%tgcc=.TRUE.
          func1%mfxcx = mfxcx_is_skipped
          func1%mfxcc=mfxcc_is_lyp
          func1%mgcx = mgcx_is_optx
          func1%mgcc = mgcc_is_lyp
          func1%mhfx = mhfx_is_skipped
       ELSEIF (fun == 6) THEN
          toldcode=.TRUE.
          cntl%tgc=.TRUE.
          cntl%tgcx=.TRUE.
          cntl%tgcc=.TRUE.
          func1%mfxcx = mfxcx_is_skipped
          func1%mfxcc=mfxcc_is_skipped
          func1%mgcx = mgcx_is_hcth
          func1%mgcc = mgcc_is_hse
          func1%mhfx = mhfx_is_skipped
       ELSEIF (fun == 7) THEN
          toldcode=.TRUE.
          cntl%tgc=.TRUE.
          cntl%tgcx=.TRUE.
          cntl%tgcc=.TRUE.
          cntl%ttau=.TRUE.
          func1%mfxcx = mfxcx_is_skipped
          func1%mfxcc=mfxcc_is_skipped
          func1%mgcx = mgcx_is_skipped
          func1%mgcc = mgcc_is_skipped
          func1%mtau = mtau_is_tpss
          func1%mhfx = mhfx_is_skipped
       ELSEIF (fun == 10 .OR. fun == 11 .OR. fun == 12) THEN
          cntl%tpotential=.TRUE.
          toldcode=.TRUE.
          cntl%tgc=.TRUE.
          cntl%tgcx=.TRUE.
          cntl%tgcc=.TRUE.
          cntl%ttau=.FALSE.
          func1%mfxcx = mfxcx_is_skipped
          func1%mfxcc=mfxcc_is_skipped
          func1%mgcx = mgcx_is_skipped
          func1%mgcc = mgcc_is_skipped
          func1%mtau = mtau_is_skipped
          func1%mhfx = mhfx_is_skipped
       ELSEIF (fun == 20) THEN
          toldcode=.TRUE.
          cntl%tgc=.TRUE.
          cntl%tgcx=.TRUE.
          cntl%tgcc=.TRUE.
          cntl%thybrid=.TRUE.
          func1%mfxcx = mfxcx_is_slaterx
          func1%mfxcc=mfxcc_is_vwn
          func1%mgcx = mgcx_is_pbex
          func1%mgcc = mgcc_is_pbec
          func1%mhfx = mhfx_is_hfx
          func3%pxlda=0.75_real_8
          func3%pxgc=0.75_real_8
          func3%pclda=1._real_8
          func3%pcgc=1._real_8
          func3%phfx=0.25_real_8
          CALL hf_init(.TRUE.)
       ELSEIF (fun == 21) THEN
          toldcode=.TRUE.
          cntl%tgc=.TRUE.
          cntl%tgcx=.TRUE.
          cntl%tgcc=.TRUE.
          cntl%thybrid=.TRUE.
          func1%mfxcx = mfxcx_is_slaterx
          func1%mfxcc=mfxcc_is_lyp
          func1%mgcx = mgcx_is_becke88
          func1%mgcc = mgcc_is_lyp
          func1%mhfx = mhfx_is_hfx
          func3%pxlda=0.8_real_8
          func3%pxgc=0.72_real_8
          func3%pclda=1._real_8
          func3%pcgc=0.81_real_8
          func3%phfx=0.2_real_8
          CALL hf_init(.TRUE.)
       ELSEIF (fun == 30) THEN
          cntl%tgc=.FALSE.
          cntl%tgcx=.FALSE.
          cntl%tgcc=.FALSE.
          func1%mfxcx = mfxcx_is_skipped
          func1%mfxcc=mfxcc_is_skipped
          func1%mgcx = mgcx_is_skipped
          func1%mgcc = mgcc_is_skipped
          func1%mhfx = mhfx_is_hfx
          CALL hf_init(.TRUE.)
       ELSE
          CALL stopgm(procedureN,"Functional not supported(tag=1)",& 
               __LINE__,__FILE__)
       ENDIF
    ELSEIF (itag < 0) THEN
       fun=-1
       ! Hybrid functionals
       IF (cntl%thybrid) THEN
          fun=20
       ELSEIF (func1%mhfx==mhfx_is_hfx) THEN
          ! Hartree Fock
          fun=30
       ELSE
          ! NONE
          IF (func1%mfxcx == mfxcx_is_skipped .AND. func1%mfxcc == mfxcc_is_skipped .AND.&
               func1%mgcx == mgcx_is_skipped .AND. func1%mgcc == mgcc_is_skipped) fun = 0
          ! LDA
          IF (func1%mfxcx == mfxcx_is_skipped .AND. func1%mfxcc == mfxcc_is_pade .AND.&
               func1%mgcx == mgcx_is_skipped .AND. func1%mgcc == mgcc_is_skipped) fun = 1
          ! BP
          IF (func1%mfxcx == mfxcx_is_skipped .AND. func1%mfxcc == mfxcc_is_pade .AND.&
               func1%mgcx == mgcx_is_becke88 .AND. func1%mgcc == mgcc_is_perdew86) fun = 2
          ! PBE
          IF (func1%mfxcx == mfxcx_is_skipped .AND. func1%mfxcc == mfxcc_is_pade .AND.&
               func1%mgcx == mgcx_is_pbex .AND. func1%mgcc == mgcc_is_pbec) fun = 3
          ! BLYP
          IF (func1%mfxcx == mfxcx_is_slaterx .AND. func1%mfxcc == mfxcc_is_lyp .AND.&
               func1%mgcx == mgcx_is_becke88 .AND. func1%mgcc == mgcc_is_lyp) fun = 4
          ! OLYP
          IF (func1%mfxcx == mfxcx_is_skipped .AND. func1%mfxcc == mfxcc_is_lyp .AND.&
               func1%mgcx == mgcx_is_optx .AND. func1%mgcc == mgcc_is_lyp) fun = 5
          ! HCTH
          IF (func1%mfxcx == mfxcx_is_skipped .AND. func1%mfxcc == mfxcc_is_skipped .AND.&
               func1%mgcx == mgcx_is_hcth .AND. func1%mgcc == mgcc_is_hse) fun = 6
          ! TPSS
          IF (func1%mfxcx == mfxcx_is_slaterx .AND. func1%mfxcc == mfxcc_is_lyp .AND.&
               func1%mgcx == mgcx_is_becke88 .AND. func1%mgcc == mgcc_is_lyp .AND. func1%mtau == mtau_is_tpss) fun = 7
       ENDIF
       IF (fun == -1) THEN
          CALL stopgm(procedureN,"Functional not supported(tag<0)",& 
               __LINE__,__FILE__)
       ENDIF
    ELSE
       CALL stopgm(procedureN,"Tag not supported",& 
            __LINE__,__FILE__)
    ENDIF

  END SUBROUTINE tdm_fun
  ! ==================================================================
END MODULE dftin_utils
