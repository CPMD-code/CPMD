MODULE lr_in_utils
  USE adat,                            ONLY: elem
  USE coor,                            ONLY: tau0
  USE dftin_utils,                     ONLY: tdm_fun
  USE error_handling,                  ONLY: stopgm
  USE inscan_utils,                    ONLY: inscan
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE linres,                          ONLY: &
       ireor, lr01, lr02, lr03, lrf4, shlct, td01, td02, td03, tshf, tshi, &
       tshl
  USE mols,                            ONLY: msnum,&
                                             mstat,&
                                             numol
  USE mp_interface,                    ONLY: mp_bcast
  USE parac,                           ONLY: parai,&
                                             paral
  USE pw_hfx_resp_types,               ONLY: hfx_resp_env
  USE readsr_utils,                    ONLY: input_string_len,&
                                             keyword_contains,&
                                             xstring
  USE soc_types,                       ONLY: nskip,&
                                             socskip
  USE spin,                            ONLY: lspin2
  USE switch_functionals_utils,        ONLY: switch_functionals
  USE system,                          ONLY: cnti,&
                                             cntl
  USE zeroing_utils,                   ONLY: zeroing
#include "sizeof.h"

  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER :: output_unit = 6

  PUBLIC :: lr_in
  PUBLIC :: tddft_input
  !public :: ms_in

CONTAINS

#ifdef __SR8000
  !option MP(P(0)), LANGLVL(SAVE(0))
#endif
  ! ==================================================================
  SUBROUTINE lr_in
    ! ==--------------------------------------------------------------==
    ! ==  THIS ROUTINE READS THE SECTION &LINRES &END ON UNIT IUNIT   ==
    ! ==--------------------------------------------------------------==
    ! ==  THE INPUT HAS TO BE IN THE FOLLOWING FORMAT                 ==
    ! ==     &LINRES                                                  ==
    ! ==        Options                                               ==
    ! ==     &END                                                     ==
    ! ==--------------------------------------------------------------==
    ! ==  LIST OF OPTIONS                                             ==
    ! ==                                                              ==
    ! ==  XC_ANALYTIC       :use analytic second and third            ==
    ! ==                     derivatives of the XC functional         ==
    ! ==                     (only available for some LDA             ==
    ! ==                     functionals)                             ==
    ! ==  XC_DD_ANALYTIC    :use analytic second derivatives of some  ==
    ! ==                     LDA and gradient corrected XC            ==
    ! ==                     functionals (third derivatives either    ==
    ! ==                     LDA or numeric (XC_ANALYTIC)             ==
    ! ==  DIFF FORMULA      :Number of points used in finite diff     ==
    ! ==    ndpoint                                                   ==
    ! ==  XC_EPS            :finite difference parameter for XC       ==
    ! ==    eps              derivative                               ==
    ! ==  MAXSTEP           :Max. number of optimization steps        ==
    ! ==    nmaxlr                                                    ==
    ! ==  HTHRS             :Threshold for Hessian                    ==
    ! ==    lr_hthrs                                                  ==
    ! ==  OPTIMIZER [SD,DIIS,PCG,AUTO] :Optimizer for LR equations    ==
    ! ==  THAUTO            :Threshold to switch from cntl%pcg to ODIIS    ==
    ! ==    thauto(1) thauto(2)                                       ==
    ! ==  ZDIIS             :size of ZDIIS buffer                     ==
    ! ==   mldiis                                                     ==
    ! ==  STEPLENGTH        :Step length for SD and cntl%pcg               ==
    ! ==    dlrs                                                      ==
    ! ==  CONVERGENCE       :Conv. criteria for LR optimisation       ==
    ! ==    tol_lr                                                    ==
    ! ==  QS_LIMIT          :Tollerance above which we use QS search  ==
    ! ==    tol_qs                                                    ==
    ! ==  GAUGE [PARA|GEN|ALL] :Gauge of response functions           ==
    ! ==--------------------------------------------------------------==
    CHARACTER(*), PARAMETER                  :: procedureN = 'lr_in'
    INTEGER, PARAMETER                       :: max_unknown_lines = 30

    CHARACTER(len=input_string_len)          :: line, previous_line, &
                                                error_message, unknown(max_unknown_lines)
    INTEGER                                  :: nbr_unknown_lines
    INTEGER                                  :: i, ierr, igauge, iunit, first, last
    LOGICAL                                  :: go_on_reading, something_went_wrong


    IF (paral%io_parent) THEN
       iunit = 5
       !
       ! Defaults
       !
       lr03%txc_analytic=.FALSE.
       lr03%txc_dd_ana=.FALSE.
       igauge=0
       lr02%xc_eps=0.0005_real_8
       lr01%nmaxlr=1000
       lr02%lr_hthrs=0.5_real_8
       lr01%lopti=0
       lr02%dlrs=0.1_real_8
       lr02%tol_lr=1.e-5_real_8
       lr02%tol_qs=1.e-4_real_8
       lr02%thauto(1)=1.e-1_real_8
       lr02%thauto(2)=1.e-3_real_8
       lr01%ndpoint=2
       lr01%mldiis=4
       !
       line              = ' '
       previous_line     = ' '
       error_message     = ' '
       nbr_unknown_lines = 0
       ! 
       ierr=inscan(iunit,'&LINRES')
       !
       IF (ierr == 0) THEN
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
             ELSEIF (keyword_contains(line,'&END')) THEN
                go_on_reading = .FALSE.
             ELSEIF (keyword_contains(line,'XC_ANALYTIC')) THEN
                ! ..analytic dmu/dn
                lr03%txc_analytic=.TRUE.
             ELSEIF (keyword_contains(line,'XC_DD_ANALYTIC')) THEN
                ! ..analytic dmu/dn (gradient corrected)
                lr03%txc_dd_ana=.TRUE.
             ELSEIF (keyword_contains(line,'XC_EPS')) THEN
                ! ..finite difference parameter for dmu/dn
                READ(iunit,*,iostat=ierr) lr02%xc_eps
                IF (ierr /= 0) THEN
                   error_message        = 'ERROR WHILE READING VALUE'
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF (keyword_contains(line,'MAXSTEP')) THEN
                ! ..max. number of optimisation steps
                READ(iunit,*,iostat=ierr) lr01%nmaxlr
                IF (ierr /= 0) THEN
                   error_message        = 'ERROR WHILE READING VALUE'
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF (keyword_contains(line,'HTHRS')) THEN
                ! ..Hessian threshold
                READ(iunit,*,iostat=ierr) lr02%lr_hthrs
                IF (ierr /= 0) THEN
                   error_message        = 'ERROR WHILE READING VALUE'
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF (keyword_contains(line,'ZDIIS')) THEN
                ! ..Size of ZDIIS buffer
                READ(iunit,*,iostat=ierr) lr01%mldiis
                IF (ierr /= 0) THEN
                   error_message        = 'ERROR WHILE READING VALUE'
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF (keyword_contains(line,'OPTIMIZER',alias='OPTIMIZATION')) THEN
                ! ..Optimizer
                IF (keyword_contains(line,'SD')) lr01%lopti=1
                IF (keyword_contains(line,'PCG')) lr01%lopti=2
                IF (keyword_contains(line,'DIIS')) lr01%lopti=3
                IF (keyword_contains(line,'AUTO')) lr01%lopti=0
             ELSEIF (keyword_contains(line,'GAUGE')) THEN
                ! ..Gauge of response functions
                IF (keyword_contains(line,'PARA')) igauge=1
                IF (keyword_contains(line,'GEN')) igauge=2
                IF (keyword_contains(line,'ALL')) igauge=3
             ELSEIF (keyword_contains(line,'STEPLENGTH')) THEN
                ! ..Step length      
                READ(iunit,*,iostat=ierr) lr02%dlrs
                IF (ierr /= 0) THEN
                   error_message        = 'ERROR WHILE READING VALUE'
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF (keyword_contains(line,'CONVERGENCE')) THEN
                ! ..Convergence criterion
                READ(iunit,*,iostat=ierr) lr02%tol_lr
                IF (ierr /= 0) THEN
                   error_message        = 'ERROR WHILE READING VALUE'
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF (keyword_contains(line,'QS_LIMIT')) THEN
                ! ..TOLERANCE FOR QS IN cntl%pcg
                READ(iunit,*,iostat=ierr) lr02%tol_qs
                IF (ierr /= 0) THEN
                   error_message        = 'ERROR WHILE READING VALUE'
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF (keyword_contains(line,'THAUTO')) THEN
                ! ..TOLERANCE FOR SWITCH FROM cntl%pcg TO ODIIS
                READ(iunit,*,iostat=ierr) lr02%thauto(1),lr02%thauto(2)
                IF (ierr /= 0) THEN
                   error_message        = 'ERROR WHILE READING VALUE'
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF(keyword_contains(line,'DIFF',and='FORMULA')) THEN
                ! ..Number of finite difference points
                READ(iunit,*,iostat=ierr) lr01%ndpoint
                IF (ierr /= 0) THEN
                   error_message        = 'ERROR WHILE READING VALUE'
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSEIF(keyword_contains(line,'TAMM-DANCOFF') .OR.&
                    keyword_contains(line,'DIAGONALIZER')) THEN
                CALL stopgm(procedureN,"Confusing input; detected keywords that belong to &TDDFT",& 
                     __LINE__,__FILE__)
             ELSE
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
          ENDDO
       ELSE
          !
          ! &LINRES is not there... (this is legit)
          !
          something_went_wrong = .false.
       ENDIF
       !
       IF (something_went_wrong) THEN
           WRITE(output_unit,'(/,1X,64("!"))')
           WRITE(output_unit,'(1X, A, 1X, A)') 'ERROR:', 'PROBLEM WHILE READING &LINRES SECTION:' 
           WRITE(output_unit,'(8X, A)') TRIM(ADJUSTL(error_message))
           IF (line /= ' ' .or. previous_line /= ' ') THEN
              WRITE(output_unit,'(8X, A)') 'THE LAST TWO LINES READ WITHIN THE SECTION WERE:'
              WRITE(output_unit,'(/,1X, A)') TRIM(ADJUSTL(previous_line))
              WRITE(output_unit,'(1X, A)') TRIM(ADJUSTL(line))
           END IF
           WRITE(output_unit,'(1X,64("!"))')
           CALL stopgm(procedureN,'Error while reading &LINRES section, cf. output file',& 
                __LINE__,__FILE__)
       ENDIF
       !
       IF (nbr_unknown_lines /= 0) THEN
          WRITE(output_unit,'(/,1X,64("="))')
          WRITE(output_unit,'(1X,A,13X,A,12X,A)') '= ','UNKNOWN KEYWORDS IN SECTION &LINRES','='
          DO i=1,nbr_unknown_lines
             previous_line = unknown(i)
             CALL xstring(previous_line,first,last)
             WRITE(output_unit,'(A,A)') ' = ',previous_line(first:last)
          ENDDO
          WRITE(output_unit,'(1X,64("="),/)')
       ENDIF
       
       CALL determine_gauge()
       CALL lr_in_report()

    ENDIF

    lr02%rea1=0._real_8
    lr01%int1=0
    lr03%log1=.FALSE.

    CALL broadcast_lr_in()

    CONTAINS

    ! ==--------------------------------------------------------------==
    SUBROUTINE determine_gauge()

       !
       ! Determine gauge
       !
       lr03%tgauge_all=.FALSE.
       IF (igauge.EQ.0) THEN      ! Default
          lr03%tpara_gauge=(.NOT.lspin2%tlse)
       ELSE IF (igauge.EQ.1) THEN ! Parallel transport
          lr03%tpara_gauge=.TRUE.
          IF (lspin2%tlse) THEN
             IF (paral%io_parent)&
                  WRITE(output_unit,*) 'WARNING! Enforcing parallel-transport gauge ',&
                  'with ROKS'
             IF (paral%io_parent)&
                  WRITE(output_unit,*) 'Results may be wrong! '
          ENDIF
       ELSE IF (igauge.EQ.2) THEN ! General between blocks
          lr03%tpara_gauge=.FALSE.
          lr03%tgauge_all=(.NOT.lspin2%tlse)
          IF (.NOT.lspin2%tlse) THEN
             IF (paral%io_parent)&
                  WRITE(output_unit,*) 'WARNING! Enforcing non-parallel-transport gauge'
             IF (paral%io_parent)&
                  WRITE(output_unit,*) 'Convergence may be slower! '
          ENDIF
       ELSE IF (igauge.EQ.3) THEN ! General for all orbitals
          lr03%tpara_gauge=.FALSE.
          lr03%tgauge_all=.TRUE.
          IF (paral%io_parent)&
               WRITE(output_unit,*) 'WARNING! Enforcing non-parallel-transport gauge ',&
               'for all orbitals'
          IF (paral%io_parent)&
               WRITE(output_unit,*) 'Convergence may be slower! '
       ENDIF

    END SUBROUTINE determine_gauge
    ! ==--------------------------------------------------------------==
    SUBROUTINE lr_in_report()

       !
       ! Write info to output
       !
       WRITE(output_unit,'(/,A,A)')&
       ' *********************** LINEAR RESPONSE',' ************************'
       IF (lr03%txc_dd_ana) THEN
          IF (lr03%txc_analytic) THEN
             WRITE(output_unit,'(A)')&
             ' Analytic second derivatives of GC XC functionals'
             WRITE(output_unit,'(A)')&
             ' Analytic third derivatives of LDA XC functionals'
          ELSE
             WRITE(output_unit,'(A)')&
             ' Analytic second derivatives of GC XC functionals'
             IF (cntl%tddft.AND.(cntl%md.OR.cntl%geopt.OR.cntl%tdebfor.OR.cntl%vibrat)) THEN
                WRITE(output_unit,'(A)')&
                ' Numeric third derivatives of GC XC functionals'
                WRITE(output_unit,'(A,T56,E10.3)')&
                ' Step size for numeric dmu/dn :',lr02%xc_eps
                WRITE(output_unit,'(A,T56,I10)')&
                ' Number of calculations for dmu/dn :',lr01%ndpoint
             ENDIF
          ENDIF
       ELSE
          IF (lr03%txc_analytic) THEN
             WRITE(output_unit,'(A)')&
             ' Analytic second derivatives of LDA XC functionals'
          ELSE
             WRITE(output_unit,'(A,T56,E10.3)')&
             ' Step size for numeric dmu/dn :',lr02%xc_eps
             WRITE(output_unit,'(A,T56,I10)')&
             ' Number of calculations for dmu/dn :',lr01%ndpoint
          ENDIF
       ENDIF
       WRITE(output_unit,'(A,T58,I8)')&
       ' Maximum number of optimisation steps:',lr01%nmaxlr
       WRITE(output_unit,'(A,T58,F8.4)')&
       ' Threshold for Hessian (Preconditioner)',lr02%lr_hthrs
       IF (lr01%lopti.EQ.0) THEN
          WRITE(output_unit,'(A,T56,A10)')&
          ' Optimizer for LR equations ',' AUTOMATIC'
          WRITE(output_unit,'(T20,A,T60,I6)') ' Size of ODIIS buffer ',cnti%mdiis
          WRITE(output_unit,'(T20,A,T60,I6)') ' Size of ZDIIS buffer ',lr01%mldiis
          WRITE(output_unit,'(T20,A,T56,E10.4)')&
          ' SWITCH FROM PCG TO ODIIS AT ',lr02%thauto(1)
          WRITE(output_unit,'(T20,A,T56,E10.4)')&
          ' Switch to full preconditioning at ',lr02%thauto(2)
       ENDIF
       IF (lr01%lopti.EQ.1)&
            WRITE(output_unit,'(A,T56,A10)')&
            ' Optimizer for LR equations ','        SD'
       IF (lr01%lopti.EQ.2)&
            WRITE(output_unit,'(A,T56,A10)')&
            ' Optimizer for LR equations ','       PCG'
       IF (lr01%lopti.EQ.2.AND..NOT.lr03%txc_analytic)&
            WRITE(output_unit,'(A,T56,E10.4)')&
            ' Tolerance for quadratic search ',lr02%tol_qs
       IF (lr01%lopti.EQ.3)&
       WRITE(output_unit,'(A,T56,A10)')&
            ' Optimizer for LR equations ','      DIIS'
       WRITE(output_unit,'(A,T56,F10.4)')&
            ' Step length',lr02%dlrs
       WRITE(output_unit,'(A,T56,E10.4)') ' Convergence criteria ',lr02%tol_lr
       IF (lr03%tpara_gauge) THEN
          WRITE(output_unit,*) 'Parallel-transport gauge for response functions'
       ELSE
          IF (lr03%tgauge_all) THEN
             WRITE(output_unit,*) 'Generic gauge for all response functions'
          ELSE
             WRITE(output_unit,*) 'Generic gauge for involved response functions'
          ENDIF
       ENDIF
       WRITE(output_unit,'(A,A)')  ' ***************************************',&
                         '*************************'

    END SUBROUTINE lr_in_report
    ! ==--------------------------------------------------------------==
    SUBROUTINE broadcast_lr_in()

    CALL mp_bcast_byte(lr01,size_in_bytes_of(lr01),parai%io_source,parai%cp_grp)
    CALL mp_bcast_byte(lr02,size_in_bytes_of(lr02),parai%io_source,parai%cp_grp)
    CALL mp_bcast_byte(lr03,size_in_bytes_of(lr03),parai%io_source,parai%cp_grp)

    END SUBROUTINE broadcast_lr_in
    ! ==--------------------------------------------------------------==
  END SUBROUTINE lr_in
  ! ==================================================================
#ifdef __SR8000
  !option MP(P(0)), LANGLVL(SAVE(0))
#endif
  ! ==================================================================
  SUBROUTINE tddft_input
    ! ==--------------------------------------------------------------==
    ! ==  THIS ROUTINE READS THE SECTION &TDDFT &END ON UNIT IUNIT    ==
    ! ==--------------------------------------------------------------==
    ! ==  THE INPUT HAS TO BE IN THE FOLLOWING FORMAT                 ==
    ! ==     &TDDFT                                                   ==
    ! ==        Options                                               ==
    ! ==     &END                                                     ==
    ! ==--------------------------------------------------------------==
    ! ==  LIST OF OPTIONS                                             ==
    ! ==                                                              ==
    ! ==  TD_METHOD_A  [functional]                                   ==
    ! ==  DIAGONALIZER [DAVIDSON,.......] :Diagonalizer for Spectra   ==
    ! ==  DAVIDSON PARAMETER                                          ==
    ! ==   ndavmax,epstdav,ndavspace                                  ==
    ! ==  DAVIDSON RDIIS                                              ==
    ! ==   nrdiismax,nrestdmax,rdiistin                               ==
    ! ==  PCG PARAMETER                                               ==
    ! ==   pcgiter,epstdav                                            ==
    ! ==   pcgiter,epstdav,pcgstep                                    ==
    ! ==  ROTATION PARAMETER                                          ==
    ! ==   rotit,epsrot,rotstep                                       ==
    ! ==  STATES [MIXED,SINGLET,TRIPLET]                              ==
    ! ==    ns                                                        ==
    ! ==  LOCALIZATION                                                ==
    ! ==  MOLECULAR STATES                                            ==
    ! ==  REORDER                                                     ==
    ! ==   nreor                                                      ==
    ! ==    i1 i2 ... in                                              ==
    ! ==  REORDER LOCAL                                               ==
    ! ==   zatom                                                      ==
    ! ==    i1 i2 ... in                                              ==
    ! ==  RANDOMIZE                                                   ==
    ! ==   tdrandom                                                   ==
    ! ==  LR-TDDFT                                                    ==
    ! ==  TAMM-DANCOFF [SUBSPACE,OPTIMIZE]                            ==
    ! ==   [msubta]                                                   ==
    ! ==  FORCE STATE                                                 ==
    ! ==   fstate                                                     ==
    ! ==  OSCILLATOR STRENGTH options: compute oscillator strengths,  ==
    ! ==     options:QP,BF,BFH,BERRY   if QP only for KS states       ==
    ! ==                               (i.e. no cntl%tddft). BF uses  ==
    ! ==                               Bessel function integration for==
    ! ==                               the transition dipoles. BFH    ==
    ! ==                               is the Hutter-Bessel option.   ==
    ! ==                               lberry uses generalised Berry  ==
    ! ==                               phase formula.                 ==
    ! ==  PROPERTY  [STATE]                                           ==
    ! ==   npstate                                                    ==
    ! ==  OUTPUT {FULL,NORMAL,MINIMAL}                                ==
    ! ==  LANDAU-ZENER [DIAG]                                         ==
    ! ==  T-SHTDDFT                                                   ==
    ! ==  EXTPOT                                                      ==
    ! ==   aampl,adir,afreq,apara1                                    ==
    ! ==  OUTPUT {FULL,NORMAL,MINIMAL}                                ==
    ! ==  LOCALCONTROL                                                ==
    ! ==  sh_lct_seed_t, sh_lct_seed_m, sh_lct_l, adir, lct_state     ==
    ! ==  STORE_V_AB                                                  ==
    ! ==  store the vab potential in the calc of the kernel           ==
    ! ==  USE_LIN_LIN                                                 ==
    ! ==  uses the Lin Lin Scheme for the speed-up of exact exchange  ==
    ! ==--------------------------------------------------------------==
    CHARACTER(*), PARAMETER                  :: procedureN = 'tddft_input'
    INTEGER, PARAMETER                       :: max_unknown_lines = 30

    CHARACTER(len=input_string_len)          :: line, previous_line, &
                                                error_message, &
                                                unknown(max_unknown_lines)
    INTEGER                                  :: i, ierr, iunit, nn, nss, first, last
    INTEGER                                  :: nbr_unknown_lines
    LOGICAL                                  :: go_on_reading, something_went_wrong


    IF (paral%io_parent) THEN
       iunit = 5
       !
       ! Defaults
       !
       td03%tda=.TRUE.
       td03%tdacanon=.FALSE.
       td01%ns_mix=0
       td01%ns_sin=0
       td01%ns_tri=0
       td01%ldiag=1
       td01%ndavmax=100
       td02%epstdav=-1._real_8
       td01%ndavspace=10
       td02%tdrandom=-1._real_8
       td01%msubta=0
       td03%tdpcgmin=.FALSE.
       td01%tdpcgiter=100
       td02%tdpcgstep=0.5_real_8
       td02%epsrot=1.e-6_real_8
       td01%rotit=50
       td02%rotstep=0.5_real_8
       td01%fstate=1
       td03%tdlocal=.FALSE.
       td03%treorder=.FALSE.
       td03%molstat=.FALSE.
       td03%tprop=.FALSE.
       td01%npstate=0
       td01%nreor=0
       td01%zatom=0
       td01%ntdiismax=20
       td01%nrestdmax=3
       td02%rdiistin=1.e-3_real_8
       td03%los    = .FALSE.
       td03%lqp    = .FALSE.
       td03%lbf    = .FALSE.
       td03%lbfh   = .FALSE.
       td03%lberry = .FALSE.
       lrf4%td_method=0
       lrf4%ks_functional=0
       lrf4%td_functional=0
       lrf4%kxc_functional=0
       td01%ioutput=0
       ! ..TSH[
       td03%tdlz=.FALSE.
       shlct%sh_lcontrol=.FALSE.
       ! [EXACT FACTORIZATION
       IF (tshl%txfmqc) THEN
          tshl%tdtully=.TRUE.
       ELSE
          tshl%tdtully=.FALSE.
       ENDIF
       ! ]EXACT FACTORIZATION
       tshl%isc=.FALSE.
       tshl%isc_2ndlr=.FALSE.
       td03%sh_diag=.FALSE.
       tshl%tdextpot=.FALSE.
       shlct%sh_lct_seed_t=100._real_8
       shlct%sh_lct_l=100._real_8
       shlct%sh_lct_seed_m=1._real_8
       shlct%lct_state=1
       nskip=1
       socskip=1
       tshl%nacv_direct_only=.TRUE.
       ! ..TSH]
       !
       line              = ' '
       previous_line     = ' '
       error_message     = ' '
       nbr_unknown_lines = 0
       !
       ierr=inscan(iunit,'&TDDFT')
       !
       IF (ierr == 0) THEN
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
             ELSEIF (keyword_contains(line,'&END')) THEN
                go_on_reading = .FALSE.
             ELSEIF ( keyword_contains(line,'DIAGONALIZER',alias='DIAGONALIZATION') ) THEN
                ! ..Diagonalizer
                IF ( keyword_contains(line,'NONHERMIT',and='DAVIDSON',alias='NONHERMITIAN') ) THEN
                   td01%ldiag=2
                ELSEIF ( keyword_contains(line,'DAVIDSON') ) THEN
                   td01%ldiag=1
                ELSEIF ( keyword_contains(line,'PCG') ) THEN
                   IF ( keyword_contains(line,'MINIMIZE') ) td03%tdpcgmin=.TRUE.
                   td01%ldiag=3
                ENDIF
             ELSE IF( keyword_contains(line,'DAVIDSON',and='PARAMETER') .OR. &
                      keyword_contains(line,'DAVIDSON',and='PARAMETERS') ) THEN
                ! ..Davidson parameters
                READ(iunit,*,iostat=ierr) td01%ndavmax,td02%epstdav,td01%ndavspace
                IF (ierr /= 0) THEN
                   error_message        = 'ERROR WHILE READING VALUE'
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSE IF( keyword_contains(line,'DAVIDSON',and='RDIIS') ) THEN
                ! ..Davidson (RDIIS) parameters
                READ(iunit,*,iostat=ierr) td01%ntdiismax,td01%nrestdmax,td02%rdiistin
                IF (ierr /= 0) THEN
                   error_message        = 'ERROR WHILE READING VALUE'
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF
             ELSE IF( keyword_contains(line,'PCG',and='PARAMETER') .OR. &
                      keyword_contains(line,'PCG',and='PARAMETERS') ) THEN
                ! ..cntl%pcg parameters
                IF (cntl%pcgmin) THEN
                   READ(iunit,*,iostat=ierr) td01%tdpcgiter,td02%epstdav
                ELSE
                   READ(iunit,*,iostat=ierr) td01%tdpcgiter,td02%epstdav,td02%tdpcgstep
                ENDIF
                IF (ierr /= 0) THEN
                   error_message        = 'ERROR WHILE READING VALUE'
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF 
             ELSE IF( keyword_contains(line,'MOLECULAR',and='STATES') ) THEN
                ! ..Use molecular states for statistical potential mixing
                td03%molstat=.TRUE.
             ELSE IF( keyword_contains(line,'FORCE',and='STATE') ) THEN
                ! ..State to calculate force from
                READ(iunit,*,iostat=ierr) td01%fstate
                IF (ierr /= 0) THEN
                   error_message        = 'ERROR WHILE READING VALUE'
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF 
             ELSE IF ( keyword_contains(line,'STATE',alias='STATES') ) THEN
                READ(iunit,*,iostat=ierr) nss
                IF (ierr /= 0) THEN
                   error_message        = 'ERROR WHILE READING VALUE'
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF 
                IF ( keyword_contains(line,'MIXED') ) td01%ns_mix=nsS
                IF ( keyword_contains(line,'SINGLET') ) td01%ns_sin=nsS
                IF ( keyword_contains(line,'TRIPLET') ) td01%ns_tri=nsS
             ELSE IF( keyword_contains(line,'ROTATION',and='PARAMETER') .OR. &
                      keyword_contains(line,'ROTATION',and='PARAMETERS') ) THEN
                ! ..Subspace rotation parameters
                READ(iunit,*,iostat=ierr) td01%rotit,td02%epsrot,td02%rotstep
                IF (ierr /= 0) THEN
                   error_message        = 'ERROR WHILE READING VALUE'
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF 
             ELSE IF ( keyword_contains(line,'RANDOMIZE') ) THEN
                ! ..Randomize parameter
                READ(iunit,*,iostat=ierr) td02%tdrandom
                IF (ierr /= 0) THEN
                   error_message        = 'ERROR WHILE READING VALUE'
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF 
             ELSE IF ( keyword_contains(line,'REORDER') ) THEN
                ! ..Reorder states
                td03%treorder=.TRUE.
                IF ( keyword_contains(line,'LOCAL') ) THEN
                   td03%tdlocal=.TRUE.
                   READ(iunit,*,iostat=ierr) td01%zatom
                   IF (ierr /= 0) THEN
                      error_message        = 'ERROR WHILE READING VALUE'
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF 
                   ALLOCATE(ireor(td01%zatom),STAT=ierr)
                   IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                        __LINE__,__FILE__)
                   READ(iunit,*,iostat=ierr) (ireor(i),i=1,td01%zatom)
                   IF (ierr /= 0) THEN
                      error_message        = 'ERROR WHILE READING VALUE'
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF 
                ELSE
                   READ(iunit,*,iostat=ierr) td01%nreor
                   IF (ierr /= 0) THEN
                      error_message        = 'ERROR WHILE READING VALUE'
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF 
                   ALLOCATE(ireor(td01%nreor),STAT=ierr)
                   IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                        __LINE__,__FILE__)
                   READ(iunit,*,iostat=ierr) (ireor(i),i=1,td01%nreor)
                   IF (ierr /= 0) THEN
                      error_message        = 'ERROR WHILE READING VALUE'
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF 
                ENDIF
             ELSE IF ( keyword_contains(line,'STORE_V_AB') ) THEN
                ! ..store vab potential in hk tddft kernel
                hfx_resp_env%store_v_ab=.TRUE.
             ELSE IF ( keyword_contains(line,'WRITE_V_AB') ) THEN
                ! ..store vab potential in hk tddft kernel
                hfx_resp_env%store_v_ab=.TRUE.
                hfx_resp_env%write_v_ab=.TRUE.
             ELSE IF ( keyword_contains(line,'USE_LIN_LIN') ) THEN
                ! .. use Cholevski decomposition for LR-TDDFT with EX
                hfx_resp_env%use_lin_lin=.TRUE.
             ELSE IF ( keyword_contains(line,'LOCAL') ) THEN
                ! ..localize states
                td03%tdlocal=.TRUE.
             ELSE IF ( keyword_contains(line,'PROPERTY') ) THEN
                ! ..Property calculation
                td03%tprop=.TRUE.
                IF ( keyword_contains(line,'STATE') ) THEN
                   READ(iunit,*,iostat=ierr) td01%npstate
                   IF (ierr /= 0) THEN
                      error_message        = 'ERROR WHILE READING VALUE'
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF 
                ENDIF
             ELSE IF ( keyword_contains(line,'TD_METHOD_A') ) THEN
                ! ..Advanced TD Methods
                lrf4%td_method=1
                IF ( keyword_contains(line,'LDA') ) THEN
                   lrf4%td_functional=1
                ELSEIF ( keyword_contains(line,'BP') ) THEN
                   lrf4%td_functional=2
                ELSEIF ( keyword_contains(line,'PBE0') ) THEN
                   lrf4%td_functional=20
                ELSEIF ( keyword_contains(line,'PBE') ) THEN
                   lrf4%td_functional=3
                ELSEIF ( keyword_contains(line,'BLYP') ) THEN
                   lrf4%td_functional=4
                ELSEIF ( keyword_contains(line,'OLYP') ) THEN
                   lrf4%td_functional=5
                ELSEIF ( keyword_contains(line,'HCTH') ) THEN
                   lrf4%td_functional=6
                ELSEIF ( keyword_contains(line,'SAOP') ) THEN
                   lrf4%td_functional=10
                ELSEIF ( keyword_contains(line,'LB94') ) THEN
                   lrf4%td_functional=11
                ELSEIF ( keyword_contains(line,'GLLB') ) THEN
                   lrf4%td_functional=12
                ELSEIF ( keyword_contains(line,'B3LYP') ) THEN
                   lrf4%td_functional=21
                ENDIF
             ELSE IF ( keyword_contains(line,'LR-TDDFT') ) THEN
                ! ..Full LR-cntl%tddft
                td03%tda=.FALSE.
             ELSE IF( keyword_contains(line,'TAMM',and='DANCOFF',alias='TAMM-DANCOFF') ) THEN
                ! ..Tamm-Dancoff Approximation
                td03%tda=.TRUE.
                IF ( keyword_contains(line,'OPTIMIZE') ) THEN
                   td03%tdacanon=.FALSE.
                ELSE
                   td03%tdacanon=.TRUE.
                ENDIF
                IF ( keyword_contains(line,'SUBSPACE') ) THEN
                   READ(iunit,*,iostat=ierr) td01%msubta
                   IF (ierr /= 0) THEN
                      error_message        = 'ERROR WHILE READING VALUE'
                      something_went_wrong = .true.
                      go_on_reading        = .false.
                   ENDIF 
                ELSE
                   td01%msubta=-1
                ENDIF
             ELSE IF( keyword_contains(line,'OSCILLATOR',and='STRENGTH') ) THEN
                td03%los = .TRUE.
                IF (( keyword_contains(line,'QP') )) td03%lqp = .TRUE.
                IF (( keyword_contains(line,'BF') )) td03%lbf = .TRUE.
                IF (( keyword_contains(line,'BFH') )) THEN
                   td03%lbfh = .TRUE.
                   td03%lbf  = .FALSE.
                ENDIF
                IF (( keyword_contains(line,'BERRY') )) THEN
                   td03%lberry = .TRUE.
                   td03%lbfh   = .FALSE.
                   td03%lbf    = .FALSE.
                   td03%lqp    = .FALSE.
                ENDIF
             ELSE IF ( keyword_contains(line,'OUTPUT') ) THEN
                IF ( keyword_contains(line,'NORMAL') ) td01%ioutput=0
                IF ( keyword_contains(line,'FULL') ) td01%ioutput=1
                IF ( keyword_contains(line,'MINIMAL') ) td01%ioutput=-1
             ELSE IF ( keyword_contains(line,'LZ-SHTDDFT') ) THEN
                ! LZ Landau-Zener surface hopping sheme
                td03%tdlz=.TRUE.
                IF ( keyword_contains(line,'DIAGONAL',alias='DIAGONALIZE') ) td03%sh_diag=.TRUE.
             ELSE IF ( keyword_contains(line,'T-SHTDDFT') ) THEN
                ! ..TSH[
                tshl%tdtully=.TRUE.
                td03%sh_diag=.FALSE.
             ELSE IF ( keyword_contains(line,'NACV_FULL') ) THEN
                tshl%nacv_direct_only=.FALSE.
             ELSE IF ( keyword_contains(line,'ISC-SHTDDFT') ) THEN
                READ(iunit,*,iostat=ierr) nskip, socskip
                IF (ierr /= 0) THEN
                   error_message        = 'ERROR WHILE READING VALUE'
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF 
                tshl%tdtully=.TRUE.
                td03%sh_diag=.FALSE.
                tshl%isc=.TRUE.
             ELSE IF ( keyword_contains(line,'EXTPOT') ) THEN
                tshl%tdextpot=.TRUE.
                READ(iunit,*,iostat=ierr) tshf%aampl, tshi%adir, tshf%afreq, tshf%apara1
                IF (ierr /= 0) THEN
                   error_message        = 'ERROR WHILE READING VALUE'
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF 
             ELSE IF ( keyword_contains(line,'SH_LCT') ) THEN
                shlct%sh_lcontrol=.TRUE.
                READ(iunit,*,iostat=ierr) shlct%sh_lct_seed_t, shlct%sh_lct_seed_m, shlct%sh_lct_l ,&
                                          shlct%adir, shlct%lct_state 
                IF (ierr /= 0) THEN
                   error_message        = 'ERROR WHILE READING VALUE'
                   something_went_wrong = .true.
                   go_on_reading        = .false.
                ENDIF 
                ! ..TSH]
             ELSE
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
          ENDDO
       ELSE
          !
          ! &TDDFT is not there... (this is legit)
          !
          something_went_wrong = .false.
       ENDIF
       !
       IF (something_went_wrong) THEN
           WRITE(output_unit,'(/,1X,64("!"))')
           WRITE(output_unit,'(1X, A, 1X, A)') 'ERROR:', 'PROBLEM WHILE READING &TDDFT SECTION:' 
           WRITE(output_unit,'(8X, A)') TRIM(ADJUSTL(error_message))
           IF (line /= ' ' .or. previous_line /= ' ') THEN
              WRITE(output_unit,'(8X, A)') 'THE LAST TWO LINES READ WITHIN THE SECTION WERE:'
              WRITE(output_unit,'(/,1X, A)') TRIM(ADJUSTL(previous_line))
              WRITE(output_unit,'(1X, A)') TRIM(ADJUSTL(line))
           END IF
           WRITE(output_unit,'(1X,64("!"))')
           CALL stopgm(procedureN,'Error while reading &TDDFT section, cf. output file',& 
                __LINE__,__FILE__)
       ENDIF
       !
       IF (nbr_unknown_lines /= 0) THEN
          WRITE(output_unit,'(/,1X,64("="))')
          WRITE(output_unit,'(1X,A,13X,A,13X,A)') '= ','UNKNOWN KEYWORDS IN SECTION &TDDFT','='
          DO i=1,nbr_unknown_lines
             previous_line = unknown(i)
             CALL xstring(previous_line,first,last)
             WRITE(output_unit,'(A,A)') ' = ',previous_line(first:last)
          ENDDO
          WRITE(output_unit,'(1X,64("="),/)')
       ENDIF
   
       CALL adjust_parameters()
       CALL tddft_report()

    ENDIF

    td02%rea2 = 0._real_8
    td01%int2 = 0
    td03%log2 = .FALSE.

    CALL broadcast_tddft()

    IF (td03%molstat) CALL ms_in

    CONTAINS

    ! ==--------------------------------------------------------------==
    SUBROUTINE adjust_parameters()

       IF (td01%ns_mix+td01%ns_sin+td01%ns_tri.EQ.0) THEN
          IF (cntl%tlsd) THEN
             td01%ns_mix=1
          ELSE
             td01%ns_sin=1
          ENDIF
       ENDIF
       IF (cntl%tlsd.AND.td01%ns_mix.EQ.0) THEN
          td01%ns_mix=1
          td01%ns_sin=0
          td01%ns_tri=0
       ENDIF
       IF (.NOT.cntl%tlsd.AND.td01%ns_sin+td01%ns_tri.EQ.0) THEN
          td01%ns_mix=0
          td01%ns_sin=1
          td01%ns_tri=0
       ENDIF
       IF (td03%tda .AND. td01%msubta.GT.0) THEN
          IF (.NOT.td03%tdacanon) td01%ldiag=3
       ENDIF
       IF (.NOT.td03%tda .AND. td01%ldiag.EQ.1) td01%ldiag=2
       IF (td02%epstdav.LT.0._real_8) THEN
          IF (td01%ldiag.EQ.1) td02%epstdav=1.e-8_real_8
          IF (td01%ldiag.EQ.2) td02%epstdav=1.e-6_real_8
          IF (td01%ldiag.EQ.3) td02%epstdav=1.e-7_real_8
       ENDIF
       IF (lrf4%td_method.NE.0) THEN
          IF (cntl%use_xc_driver) CALL stopgm(procedureN,'TD_METHOD_A is not yet supported by the new driver', &
                                                         __LINE__, __FILE__)
          CALL tdm_fun(lrf4%ks_functional,-1)
          CALL switch_functionals
          CALL tdm_fun(lrf4%kxc_functional,-1)
          CALL switch_functionals
       ENDIF
       ! ..TSH[
       tshi%nroot_sh=MAX(td01%ns_sin,td01%ns_tri,td01%ns_mix,2)
       ! ..TSH]

    END SUBROUTINE adjust_parameters
    ! ==--------------------------------------------------------------==
    SUBROUTINE tddft_report()
       !
       ! Write to output
       ! 
       WRITE(output_unit,'(/,A,A)')&
            ' ***************************  TDDFT ',&
            ' ****************************'
       IF (lrf4%td_method == 1) THEN
          WRITE(output_unit,'(A)') ' Use Time-Dependent DFT Perturbation Method A'
          IF (lrf4%td_functional == 1) THEN
             WRITE(output_unit,'(T34,A,T58,A)') ' STATE FUNCTIONAL','     LDA'
          ELSEIF (lrf4%td_functional == 2) THEN
             WRITE(output_unit,'(T34,A,T58,A)') ' STATE FUNCTIONAL','      BP'
          ELSEIF (lrf4%td_functional == 3) THEN
             WRITE(output_unit,'(T34,A,T58,A)') ' STATE FUNCTIONAL','     PBE'
          ELSEIF (lrf4%td_functional == 4) THEN
             WRITE(output_unit,'(T34,A,T58,A)') ' STATE FUNCTIONAL','    BLYP'
          ELSEIF (lrf4%td_functional == 5) THEN
             WRITE(output_unit,'(T34,A,T58,A)') ' STATE FUNCTIONAL','    OLYP'
          ELSEIF (lrf4%td_functional == 6) THEN
             WRITE(output_unit,'(T34,A,T58,A)') ' STATE FUNCTIONAL','    HCTH'
          ELSEIF (lrf4%td_functional == 10) THEN
             WRITE(output_unit,'(T34,A,T58,A)') ' STATE FUNCTIONAL','    SAOP'
             IF (td03%molstat) THEN
                WRITE(output_unit,'(T45,A)') 'WITH MOLECULAR STATES'
             ENDIF
          ELSEIF (lrf4%td_functional == 11) THEN
             WRITE(output_unit,'(T34,A,T58,A)') ' STATE FUNCTIONAL','    LB94'
          ELSEIF (lrf4%td_functional == 12) THEN
             WRITE(output_unit,'(T34,A,T58,A)') ' STATE FUNCTIONAL','    GLLB'
          ELSEIF (lrf4%td_functional == 20) THEN
             WRITE(output_unit,'(T34,A,T58,A)') ' STATE FUNCTIONAL','    PBE0'
          ELSEIF (lrf4%td_functional == 21) THEN
             WRITE(output_unit,'(T34,A,T58,A)') ' STATE FUNCTIONAL','   B3LYP'
          ELSE
             CALL stopgm(procedureN,"TD_METHOD_A: State functional not supported.",& 
                  __LINE__,__FILE__)
          ENDIF
          IF (lrf4%kxc_functional == 0) THEN
             WRITE(output_unit,'(T34,A,T58,A)') ' KERNEL FUNCTIONAL','    NONE'
          ELSEIF (lrf4%kxc_functional == 1) THEN
             WRITE(output_unit,'(T34,A,T58,A)') ' KERNEL FUNCTIONAL','     LDA'
          ELSEIF (lrf4%kxc_functional == 2) THEN
             WRITE(output_unit,'(T34,A,T58,A)') ' KERNEL FUNCTIONAL','      BP'
          ELSEIF (lrf4%kxc_functional == 3) THEN
             WRITE(output_unit,'(T34,A,T58,A)') ' KERNEL FUNCTIONAL','     PBE'
          ELSEIF (lrf4%kxc_functional == 4) THEN
             WRITE(output_unit,'(T34,A,T58,A)') ' KERNEL FUNCTIONAL','    BLYP'
          ELSEIF (lrf4%kxc_functional == 5) THEN
             WRITE(output_unit,'(T34,A,T58,A)') ' KERNEL FUNCTIONAL','    OLYP'
          ELSEIF (lrf4%kxc_functional == 6) THEN
             WRITE(output_unit,'(T34,A,T58,A)') ' KERNEL FUNCTIONAL','    HCTH'
          ELSEIF (lrf4%kxc_functional == 10) THEN
             WRITE(output_unit,'(T34,A,T58,A)') ' KERNEL FUNCTIONAL','    SAOP'
          ELSEIF (lrf4%kxc_functional == 20) THEN
             WRITE(output_unit,'(T34,A,T58,A)') ' KERNEL FUNCTIONAL','    PBE0'
          ELSEIF (lrf4%kxc_functional == 21) THEN
             WRITE(output_unit,'(T34,A,T58,A)') ' KERNEL FUNCTIONAL','   B3LYP'
          ELSE
             CALL stopgm(procedureN,"TD_METHOD_A: Kernel not supported.",& 
                  __LINE__,__FILE__)
          ENDIF
       ENDIF
       IF (lr03%txc_dd_ana) THEN
          WRITE(output_unit,'(A)')&
               ' Analytic second derivatives of GC XC functionals'
       ELSEIF (lr03%txc_analytic) THEN
          WRITE(output_unit,'(A)')&
               ' Analytic second derivatives of LDA XC functionals'
       ELSE
          WRITE(output_unit,'(A,T56,E10.3)')&
               ' Step size for numeric dmu/dn :',lr02%xc_eps
          WRITE(output_unit,'(A,T56,I10)')&
               ' Number of calculations for dmu/dn :',lr01%ndpoint
       ENDIF
       IF (td03%tda) THEN
          IF (td01%msubta > 0) THEN
             WRITE(output_unit,'(A,T56,I10)')&
                  ' Tamm-Dancoff Approximation (Subspace) ',td01%msubta
             IF (td03%tdacanon) THEN
                WRITE(output_unit,'(T36,A)') ' Fixed Subspace Approximation'
             ELSE
                WRITE(output_unit,'(T36,A)') ' Optimized Subspace Approximation'
             ENDIF
          ELSE
             WRITE(output_unit,'(A)') ' Tamm-Dancoff Approximation '
          ENDIF
       ENDIF
       IF (td01%ldiag == 1) THEN
          WRITE(output_unit,'(A,T56,A10)')&
               ' Diagonalization Method ','  DAVIDSON'
          WRITE(output_unit,'(T10,A,T56,I10)')&
               ' Max. number of iterations',td01%ndavmax
          WRITE(output_unit,'(T10,A,T56,E10.3)') ' Convergence criteria',td02%epstdav
          WRITE(output_unit,'(T10,A,T56,I10)') ' Max. size of Davidson matrix',&
               td01%ndavspace
          WRITE(output_unit,'(T10,A)') ' RDIIS Parameters'
          WRITE(output_unit,'(T10,A,I4,T30,A,I4,T46,A,E10.3)')&
               ' NTDIISMAX=',td01%ntdiismax,' NRESTDMAX=',td01%nrestdmax,&
               ' RDIISTIN=',td02%rdiistin
       ELSEIF (td01%ldiag == 2) THEN
          WRITE(output_unit,'(A,T42,A24)')&
               ' Diagonalization Method ','  NON-HERMITIAN DAVIDSON'
          WRITE(output_unit,'(T10,A,T56,I10)')&
               ' Max. number of iterations',td01%ndavmax
          WRITE(output_unit,'(T10,A,T56,E10.3)') ' Convergence criteria',td02%epstdav
          WRITE(output_unit,'(T10,A,T56,I10)') ' Max. size of Davidson matrix',&
               td01%ndavspace
       ELSEIF (td01%ldiag == 3) THEN
          WRITE(output_unit,'(A,T42,A24)')&
               ' Diagonalization Method ','    CONJUGATED GRADIENTS'
          WRITE(output_unit,'(T10,A,T56,I10)')&
               ' Max. number of iterations',td01%tdpcgiter
          IF (td03%tdpcgmin) THEN
             WRITE(output_unit,'(T10,A)') ' Line minimization'
          ELSE
             WRITE(output_unit,'(T10,A,T56,E10.3)') ' Fixed step size',td02%tdpcgstep
          ENDIF
          WRITE(output_unit,'(T10,A,T56,E10.3)') ' Convergence criteria',td02%epstdav
       ENDIF
       IF (td03%tda.AND.td01%msubta > 0.AND..NOT.td03%tdacanon) THEN
          WRITE(output_unit,'(A)') ' Subspace Rotation by Exponentials'
          WRITE(output_unit,'(T10,A,T56,I10)') ' Max. number of iterations',td01%rotit
          WRITE(output_unit,'(T10,A,T56,E10.3)') ' Convergence criteria',td02%epsrot
          WRITE(output_unit,'(T10,A,T56,E10.3)') ' Fixed step size',td02%rotstep
       ENDIF
       IF (td01%ns_sin > 0) &
          WRITE(output_unit,'(A,T56,I10)') ' Number of Singlet States ',td01%ns_sin
       IF (td01%ns_tri > 0) &
          WRITE(output_unit,'(A,T56,I10)') ' Number of Triplet States ',td01%ns_tri
       IF (td01%ns_mix > 0) &
          WRITE(output_unit,'(A,T56,I10)') ' Number of Mixed States ',td01%ns_mix
       WRITE(output_unit,'(A,T56,I10)') ' Forces calculated for state ',td01%fstate
       IF (td02%tdrandom > 0._real_8) THEN
          WRITE(output_unit,'(A,T51,E15.6)')&
               ' Randomize Initial Guess (Amplitude)',td02%tdrandom
       ENDIF
       IF (td03%tdlocal) THEN
          WRITE(output_unit,'(A)') ' Localize ground state orbitals'
       ENDIF
       IF (td03%treorder) THEN
          WRITE(output_unit,'(A)') ' Reorder ground state orbitals'
          IF (td01%zatom /= 0) THEN
             WRITE(output_unit,'(T40,A)') "w.r.t. distance to atoms"
             WRITE(output_unit,'(T20,10I5)') (ireor(i),i=1,td01%zatom)
          ELSE
             WRITE(output_unit,'(T20,10I5)') (ireor(i),i=1,td01%nreor)
          ENDIF
       ENDIF
       IF (td03%los) THEN
          IF (.NOT.td03%lqp) THEN
             WRITE(output_unit,'(A)') ' COMPUTE TDDFT OSCILLATOR STRENGTHS'
             IF ((td03%lberry))&
                WRITE(output_unit,'(A)') '      Use Berry phase formula'
          ELSE
             WRITE(output_unit,'(A)')&
                  ' COMPUTE ONLY KOHN-SHAM OSCILLATOR STRENGTHS (NO TDDFT)'
          ENDIF
       ENDIF
       IF (td03%tprop) THEN
          IF (td01%npstate == 0) THEN
             WRITE(output_unit,'(A)')&
                  ' Compute excited state properties for all states'
          ELSE
             WRITE(output_unit,'(A,T63,I3)')&
                  ' Compute excited state properties for state',td01%npstate
          ENDIF
       ENDIF
       IF (td03%tdlz) THEN
          WRITE(output_unit,'(A)')&
               ' LANDAU-ZENER SURFACE HOPPING TDDFT DYNAMICS'
       ENDIF
       ! ..TSH[
       IF (tshl%tdtully) THEN
          WRITE(output_unit,'(A)') ' TDDFT FEWEST SWITCHES: TULLYS SURFACE HOPPING'
          ! ..XFMQC[    MIN       JUST FOR MORE INFORMATION IN OUTPUT 
          IF (tshl%txfmqc .AND. paral%io_parent) &
          WRITE(output_unit,'(A)') ' EXACT FACTORIZATION CORRECTION'
          ! ..XFMQC]
       ENDIF
       ! ..TSH]
       WRITE(output_unit,'(A,A)')&
            ' ***************************************',&
            '*************************'

    END SUBROUTINE tddft_report
    ! ==--------------------------------------------------------------==
    SUBROUTINE broadcast_tddft()

       CALL mp_bcast_byte(td02,size_in_bytes_of(td02),parai%io_source,parai%cp_grp)
       CALL mp_bcast_byte(td01,size_in_bytes_of(td01),parai%io_source,parai%cp_grp)
       CALL mp_bcast_byte(td03,size_in_bytes_of(td03),parai%io_source,parai%cp_grp)
       CALL mp_bcast_byte(hfx_resp_env,size_in_bytes_of(hfx_resp_env),parai%io_source,parai%cp_grp)

       IF (td03%treorder) THEN
          IF (.NOT.paral%parent) THEN
             nn=MAX(td01%zatom,td01%nreor)
             ALLOCATE(ireor(nn),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
          ENDIF
          CALL mp_bcast(ireor,SIZE(ireor),parai%io_source,parai%cp_grp)
       ENDIF
       CALL mp_bcast_byte(lrf4, size_in_bytes_of(lrf4),parai%io_source,parai%cp_grp)
       ! ..TSH[
       CALL mp_bcast(tshl%tdtully,parai%io_source,parai%cp_grp)
       CALL mp_bcast(tshl%nacv_direct_only,parai%io_source,parai%cp_grp)
       CALL mp_bcast(tshl%isc,parai%io_source,parai%cp_grp)
       CALL mp_bcast(tshl%isc_2ndlr,parai%io_source,parai%cp_grp)
       CALL mp_bcast(tshl%tdextpot,parai%io_source,parai%cp_grp)
       CALL mp_bcast(tshi%nroot_sh,parai%io_source,parai%cp_grp)
       CALL mp_bcast(tshi%adir,parai%io_source,parai%cp_grp)
       CALL mp_bcast(shlct%sh_lcontrol,parai%io_source,parai%cp_grp)
       CALL mp_bcast(shlct%adir,parai%io_source,parai%cp_grp)
       CALL mp_bcast(nskip,parai%io_source,parai%cp_grp)
       CALL mp_bcast(socskip,parai%io_source,parai%cp_grp)

    END SUBROUTINE broadcast_tddft
    ! ==--------------------------------------------------------------==
  END SUBROUTINE tddft_input
  ! ==================================================================
#ifdef __SR8000
  !option MP(P(0)), LANGLVL(SAVE(0))
#endif
  ! ==================================================================
  SUBROUTINE ms_in
    ! ==--------------------------------------------------------------==
    ! ==  THIS ROUTINE READS THE SECTION &MOLSTATES &END ON UNIT IUNIT==
    ! ==--------------------------------------------------------------==
    ! ==  THE INPUT HAS TO BE IN THE FOLLOWING FORMAT                 ==
    ! ==     &MOLSTATES                                               ==
    ! ==        Options                                               ==
    ! ==     &END                                                     ==
    ! ==--------------------------------------------------------------==
    ! ==  LIST OF OPTIONS                                             ==
    ! ==                                                              ==
    ! ==  MOLECULES                                                   ==
    ! ==  n1 n2 n3 ... (all atom)                                     ==
    ! ==  STATES                                                      ==
    ! ==  na nb nc ... (all molecules)                                ==
    ! ==                                                              ==
    ! ==--------------------------------------------------------------==
    CHARACTER(*), PARAMETER                  :: procedureN = 'ms_in'

    INTEGER, PARAMETER                       :: max_unknown_lines = 30

    CHARACTER(len=input_string_len)          :: line, previous_line, &
                                                error_message, unknown(max_unknown_lines)
    INTEGER                                  :: nbr_unknown_lines, first, last
    INTEGER                                  :: i, ia, iat, ierr, is, iunit, &
                                                k, len
    LOGICAL                                  :: go_on_reading, something_went_wrong


    len=ions1%nat
    ALLOCATE(msnum(len),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(mstat(len),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(msnum)
    CALL zeroing(mstat)
    numol=0

    !
    IF (paral%io_parent) THEN
       iunit = 5
       !
       line              = ' '
       previous_line     = ' '
       error_message     = ' '
       nbr_unknown_lines = 0
       !
       ierr=inscan(iunit,'&MOLSTATES')
       !
       IF (ierr == 0) THEN
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
             ELSEIF (keyword_contains(line,'&END')) THEN
                go_on_reading = .FALSE.
             ELSEIF (keyword_contains(line,'MOLECULES')) THEN
                READ(iunit,*,iostat=ierr) (msnum(i),i=1,ions1%nat)
                IF (ierr /= 0) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .TRUE.
                   go_on_reading        = .FALSE.
                ENDIF
                DO i=1,ions1%nat
                   IF (numol.LT.msnum(i)) numol=msnum(i)
                ENDDO
             ELSEIF (keyword_contains(line,'STATES')) THEN
                IF (numol <= 0) CALL stopgm(procedureN,"Molecules not defined",& 
                     __LINE__,__FILE__)
                READ(iunit,*,iostat=ierr) (mstat(i),i=1,numol)
                IF (ierr /= 0) THEN
                   error_message        = "ERROR WHILE READING VALUE"
                   something_went_wrong = .TRUE.
                   go_on_reading        = .FALSE.
                ENDIF
             ELSE
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
          ENDDO
       ELSE
          !
          ! &MOLSTATES is not there... (this is legit)
          !
          something_went_wrong = .false.
       ENDIF
       !
       IF (something_went_wrong) THEN
           WRITE(output_unit,'(/,1X,64("!"))')
           WRITE(output_unit,'(1X, A, 1X, A)') 'ERROR:', 'PROBLEM WHILE READING &MOLSTATES SECTION:' 
           WRITE(output_unit,'(8X, A)') TRIM(ADJUSTL(error_message))
           IF (line /= ' ' .or. previous_line /= ' ') THEN
              WRITE(output_unit,'(8X, A)') 'THE LAST TWO LINES READ WITHIN THE SECTION WERE:'
              WRITE(output_unit,'(/,1X, A)') TRIM(ADJUSTL(previous_line))
              WRITE(output_unit,'(1X, A)') TRIM(ADJUSTL(line))
           END IF
           WRITE(output_unit,'(1X,64("!"))')
           CALL stopgm(procedureN,'Error while reading &MOLSTATES section, cf. output file',& 
                __LINE__,__FILE__)
       ENDIF
       !
       IF (nbr_unknown_lines /= 0) THEN
          WRITE(output_unit,'(/,1X,64("="))')
          WRITE(output_unit,'(1X,A,11X,A,11X,A)') '= ','UNKNOWN KEYWORDS IN SECTION &MOLSTATES','='
          DO i=1,nbr_unknown_lines
             previous_line = unknown(i)
             CALL xstring(previous_line,first,last)
             WRITE(output_unit,'(A,A)') ' = ',previous_line(first:last)
          ENDDO
          WRITE(output_unit,'(1X,64("="),/)')
       ENDIF

       CALL report_ms_in()

    ENDIF

    CALL broadcast_ms_in()


    CONTAINS

    ! ==--------------------------------------------------------------==
    SUBROUTINE report_ms_in()

       WRITE(output_unit,'(/,A,A)')&
          ' *********************** MOLECULAR STATES',&
          ' ***********************'
       DO i=1,numol
          WRITE(output_unit,'(A,I5,T52,A,I5)') '  MOLECULE :',i,'STATES :',mstat(i)
          iat=0
          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                iat=iat+1
                IF (i.EQ.msnum(iat)) THEN
                   WRITE(output_unit,'(T6,I4,4X,A2,4X,3(F15.6))')&
                        iat,elem%el(ions0%iatyp(is)),(tau0(k,ia,is),k=1,3)
                ENDIF
             ENDDO
          ENDDO
       ENDDO
       WRITE(output_unit,'(A,A)')&
            ' ***************************************',&
            '*************************'

    END SUBROUTINE report_ms_in
    ! ==--------------------------------------------------------------==
    SUBROUTINE broadcast_ms_in()

       CALL mp_bcast(numol,parai%io_source,parai%cp_grp)
       CALL mp_bcast(msnum,ions1%nat,parai%io_source,parai%cp_grp)
       CALL mp_bcast(mstat,numol,parai%io_source,parai%cp_grp)

    END SUBROUTINE broadcast_ms_in
    ! ==--------------------------------------------------------------==
  END SUBROUTINE ms_in
  ! ==================================================================

END MODULE lr_in_utils
