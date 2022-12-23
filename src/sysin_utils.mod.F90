#ifdef __SR8000
!option MP(P(0)), LANGLVL(SAVE(0))
#endif

MODULE sysin_utils
  USE cdftmod,                         ONLY: cdft_mat,&
                                             cdftlog,&
                                             cdftpi,&
                                             wcfix,&
                                             wg_n,&
                                             wg_sigma,&
                                             wgaussl
  USE cell,                            ONLY: cell_com,&
                                             lcell
  USE clas,                            ONLY: clas8
  USE cnst,                            ONLY: au_fs,&
                                             fbohr,&
                                             kb_au,&
                                             pi
  USE cplngs_utils,                    ONLY: cplngs_init
  USE cplngsmod,                       ONLY: &
       csurf, eps_c, f_high, f_low, f_med, fc_high, fc_low, fc_med, iatfd, &
       isurf, natfd, nsurf, nvect, tallat, talldof, tcpl, tcplfd, tcpllr, &
       tolcpl, tspecv
  USE ddip,                            ONLY: aifield,&
                                             lenbk
  USE dg,                              ONLY: edgcomm,&
                                             tdgcomm
  USE elct,                            ONLY: crge
  USE elct2,                           ONLY: tfixo
  USE error_handling,                  ONLY: stopgm
  USE fcas,                            ONLY: init_flag
  USE gvec,                            ONLY: nhgwish,&
                                             tnhgwish
  USE hfxmod,                          ONLY: hfxc2,&
                                             ipoolhfx
  USE inscan_utils,                    ONLY: inscan
  USE isos,                            ONLY: isos1,&
                                             isos2,&
                                             isos3
  USE kdpc,                            ONLY: bmix,&
                                             tkdp
  USE kinds,                           ONLY: real_8
  USE kpnt,                            ONLY: rk,&
                                             wk
  USE kpts,                            ONLY: kpts_com,&
                                             tkpts,&
                                             wvk0
  USE mp_interface,                    ONLY: mp_bcast
  USE parac,                           ONLY: parai,&
                                             paral
  USE prcp,                            ONLY: prcp_com,&
                                             prcpl
  USE readsr_utils,                    ONLY: input_string_len,&
                                             index_of_delimiter,&
                                             keyword_contains,&
                                             readsi,&
                                             readsr,&
                                             xstring
  USE ropt,                            ONLY: iteropt,&
                                             ropt_mod
  USE shock,                           ONLY: shock1
  USE sphe,                            ONLY: tsphere
  USE spin,                            ONLY: clsd,&
                                             lspin1,&
                                             lspin2,&
                                             lspin3,&
                                             lspin4,&
                                             spin_mod
  USE symm,                            ONLY: stag,&
                                             sxscale,&
                                             symmi,&
                                             symmr,&
                                             symmt,&
                                             syscale,&
                                             szscale,&
                                             tcartesian
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             cntr,&
                                             dual00,&
                                             nkpt,&
                                             parm
  USE zeroing_utils,                   ONLY: zeroing
#include "sizeof.h"

  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER :: output_unit = 6

  PUBLIC :: sysin

CONTAINS

  ! ==================================================================
  SUBROUTINE sysin
    ! ==--------------------------------------------------------------==
    ! ==  THIS ROUTINE READS THE SECTION &SYSTEM FROM IUNIT           ==
    ! ==--------------------------------------------------------------==
    ! ==  DESCRIPTION OF INPUT                                        ==
    ! ==                                                              ==
    ! ==  &SYSTEM                                                     ==
    ! ==     OPTIONS                                                  ==
    ! ==  &END                                                        ==
    ! ==                                                              ==
    ! ==  LIST OF OPTIONS                                             ==
    ! ==    SYMMETRY                                                  ==
    ! ==      [ibrav] [name]                                          ==
    ! ==    CHECK SYMMETRY [OFF]                                      ==
    ! ==      [boxeps]                                                ==
    ! ==    POINT GROUP [MOLECULE] [DELTA=deltasym]                   ==
    ! ==      [indpg] [NAME=name] [AUTO]                              ==
    ! ==    CELL [ABSOLUTE] [DEGREE] [VECTORS]                        ==
    ! ==      celldm(1.. 6)                                          ==
    ! ==    REFERENCE CELL [ABSOLUTE] [DEGREE]                        ==
    ! ==      cellrf(1.. 6)                                          ==
    ! ==    CLASSICAL CELL [ABSOLUTE] [DEGREE]                        ==
    ! ==      cellcl(1.. 6)                                          ==
    ! ==    ISOTROPIC CELL                                            ==
    ! ==    ZFLEXIBLE CELL                                            ==
    ! ==    POISSON [HOCKNEY,TUCKERMAN,MORTENSEN,PARAMETER]           ==
    ! ==       coulomb smoothing (Hockney)                            ==
    ! ==       alpha*L (Tuckerman)                                    ==
    ! ==    BOX WALLS                                                 ==
    ! ==      wall_skin                                               ==
    ! ==    LOW SPIN EXCITATION [ROKS,ROSS,ROOTHAAN,CAS22,PENALTY]    ==
    ! ==      apenal                                                  ==
    ! ==    LOW SPIN EXCITATION LSETS                                 ==
    ! ==    LSE PARAMETERS                                            ==
    ! ==      lsea lseb       (see control_def for default values)    ==
    ! ==    MODIFIED GOEDECKER [PARAMETERS]                           ==
    ! ==      mgmab mgmba                                             ==
    ! ==    ENERGY PROFILE                                            ==
    ! ==    EXTERNAL FIELD                                            ==
    ! ==      ex ey ez                                                ==
    ! ==    COUPLINGS [FD,FD=eps_c,PROD=eps_c]                        ==
    ! ==    COUPLINGS NAT=nat                                         ==
    ! ==       iatfd(1),...                                           ==
    ! ==    COUPLINGS NSURF=nsurf                                     ==
    ! ==       kssta(1,1)   kssta(2,1)   coeff(1)                     ==
    ! ==       kssta(1,...) kssta(2,...) coeff(...)                   ==
    ! ==    COUPLINGS LINRES [TOL=<tol>] [NVECT=<nvect> [SPECIFY]]    ==
    ! ==    COUPLINGS LINRES THRESHOLDS                               ==
    ! ==      fc_thr_low  f_low                                       ==
    ! ==      fc_thr_med  f_med                                       ==
    ! ==      fc_thr_high f_high                                      ==
    ! ==    CHARGE                                                    ==
    ! ==      charge                                                  ==
    ! ==    MULTIPLICITY                                              ==
    ! ==      nspin (2S+1)                                            ==
    ! ==    NSUP                                                      ==
    ! ==      nsup                                                    ==
    ! ==    STATES                                                    ==
    ! ==      n                                                       ==
    ! ==    OCCUPATION [FIXED]                                        ==
    ! ==      f(i),i=1,n                                              ==
    ! ==    [SPHERICAL,NOSPHERICAL] CUTOFF                            ==
    ! ==      ecut                                                    ==
    ! ==    DENSITY CUTOFF                                            ==
    ! ==      ecutd                                                   ==
    ! ==    DUAL                                                      ==
    ! ==      cdual                                                   ==
    ! ==    CONSTANT CUTOFF                                           ==
    ! ==      akin skin eckin                                         ==
    ! ==    DENSITY CUTOFF NUMBER                                     ==
    ! ==      nhgwish                                                 ==
    ! ==    KPOINT [SCALED] [BLOCK=nkpnt] [ALL] [CALCULATED] [NOSWAP] ==
    ! ==           [TONLYDIAG]                                        ==
    ! ==      nkpts                                                   ==
    ! ==      rk(1.. 3) wk (1)                                       ==
    ! ==      rk(1.. 3) wk (nkpnt)                                   ==
    ! ==    KPOINT MONKHORST-PACK [SYMMETRIZE] [FULL] [KDP] ...       ==
    ! ==                          [BLOCK=nkpnt] [ALL]                 ==
    ! ==      nk1 nk2 nk3 [SHIFT=wvk0(1) wvk0(2) wvk0(3)]             ==
    ! ==    KPOINT BANDS [SCALED]                                     ==
    ! ==      nkpbd rk(1..3) rk(1..3)                                 ==
    ! ==      .                                                       ==
    ! ==      .                                                       ==
    ! ==      0     0  0  0  0  0  0                                  ==
    ! ==    MESH                                                      ==
    ! ==      nr1 nr2 nr3                                             ==
    ! ==    SCALE CARTESIAN [S=sascale]                               ==
    ! ==                    [SX=sxscale] [SY=syscale] [SZ=szscale]    ==
    ! ==    DOUBLE GRID [ON,OFF]                                      ==
    ! ==    SYMMETRIZE COORDINATES                                    ==
    ! ==    ANGSTROM                                                  ==
    ! ==    TESR                                                      ==
    ! ==     iesr                                                     ==
    ! ==    SURFACE [XY, YZ, ZX]                                      ==
    ! ==    POLYMER                                                   ==
    ! ==    CLUSTER                                                   ==
    ! ==    PRESSURE                                                  ==
    ! ==     druck                                                    ==
    ! ==    SHOCK VELOCITY                                            ==
    ! ==     vshock                                                   ==
    ! ==    STRESS TENSOR                                             ==
    ! ==     t11 t12 t13                                              ==
    ! ==     t21 t22 t23                                              ==
    ! ==     t31 t32 t33                                              ==
    ! ==    HFX CUTOFF                                                ==
    ! ==     hfxwfe hfxdee                                            ==
    ! ==    DONOR                                                     ==
    ! ==     NDONOR n1 n2 n3 n4                                       ==
    ! ==    ACCEPTOR                                                  ==
    ! ==     NACCEPTOR n1 n2 n3 n4                                    ==
    ! ==    WGAUSS [NSP or 1]                                         ==
    ! ==      sigma1                                                  ==
    ! ==      ...                                                     ==
    ! ==      sigmaN                                                  ==
    ! ==                                                              ==
    ! ==--------------------------------------------------------------==
    CHARACTER(*), PARAMETER                  :: procedureN = 'sysin'
    CHARACTER(len=3), DIMENSION(32), PARAMETER :: pgrd = (/'c1 ','ci ','c2 ',&
      'c1h','c2h','c3 ','c3i','d3 ','c3v','d3 ','c4 ','s4 ','c4h','d4 ','c4v',&
      'd2d','d4h','c6 ','c3h','c6h','d6 ','c6v','d3h','d6h','d2 ','c2v','d2h',&
      't  ','th ','o  ','td ','oh '/)
    CHARACTER(len=5), DIMENSION(32), PARAMETER :: pgrp = (/'    1','  <1>',&
      '    2','    m','  2/m','    3','  <3>','   32','   3m',' <3>m','    4',&
      '  <4>','  4/m','  422','  4mm','<4>2m','4/mmm','    6','  <6>','  6/m',&
      '  622','  6mm','<6>m2','6/mmm','  222','  mm2','  mmm','   23','   m3',&
      '  432','<4>3m','  m3m'/)
    INTEGER, PARAMETER                       :: max_unknown_lines = 30 

    CHARACTER(len=10000)                     :: linelong
    CHARACTER(len=input_string_len)          :: line, error_message, previous_line, &
                                                unknown(max_unknown_lines)
    INTEGER                                  :: first, last, ierr, iunit, nbr_unknown_lines
    INTEGER                                  :: i, iequal, ik, in1, in2, &
                                                ip1, ip2, nkpbd
    LOGICAL                                  :: erread, tpressure, tstrtensor, &
                                                go_on_reading, something_went_wrong
    REAL(real_8)                             :: boxeps, ecutd, rkf(3), &
                                                rki(3), xa
    REAL(real_8), ALLOCATABLE                :: rkt(:,:)


    init_flag = 0
    !
    ! Don't rely on the compiler to initialize LENBK properly.
    lenbk = 0
    !
    ! Block data
    !
    CALL cplngs_init()
    
    IF (paral%io_parent) THEN
       iunit = 5
       !
       ! Variables for reading
       !
       nbr_unknown_lines = 0
       line              = ' '
       previous_line     = ' '
       error_message     = ' '
       !
       ! Defaults (not in types)
       !
       boxeps = 1.0e-4_real_8
       ecutd = -1._real_8
       ipoolhfx = 0
       nhgwish = 0
       tcpl = .FALSE. ! see CPLNGS_INIT
       tfixo = .FALSE.
       tkdp = .FALSE.
       tnhgwish = .FALSE.
       tpressure = .FALSE.
       tsphere = .TRUE.
       tstrtensor = .FALSE.
       !
       ! Defaults (in types), now in alphabetical order
       !
       cell_com%ar1(:) = 0.0_real_8
       cell_com%ar2(:) = 0.0_real_8
       cell_com%ar3(:) = 0.0_real_8
       cntl%bohr = .TRUE.
       cntl%tfield = .FALSE.
       cntl%tscale = .FALSE.
       crge%charge = 0.0_real_8
       crge%n = 0
       dual00%cdual = 4.0_real_8
       dual00%dual = .FALSE.
       hfxc2%hfxdee  =  0.0_real_8
       hfxc2%hfxwfe  =  0.0_real_8
       isos1%tclust = .FALSE.
       isos1%toned = .FALSE.
       isos1%ttwod = .FALSE.
       isos1%twall = .FALSE.    ! reflecting walls
       isos2%alphal = 7._real_8
       isos2%rcerf = 0._real_8
       isos2%wall_skin = 7.0_real_8
       isos3%ps_type = 0
       isos3%snormal = 3
       iteropt%iesr = 0
       lcell%tcellrefvec  =  .FALSE.
       lcell%tcellvectors = .FALSE.
       lspin2%tcas22 = .FALSE.
       lspin2%teprof = .FALSE.
       lspin2%tlsets = .FALSE.
       lspin2%tmgoe  =  .FALSE.
       lspin2%tpenal = .FALSE.
       lspin2%troot = .FALSE.
       lspin2%tross = .FALSE.
       lspin3%mgmab  =  -0.5_real_8
       lspin3%mgmba  =   0.5_real_8
       lspin4%apenal = 1._real_8
       nkpt%nblkp = 1
       nkpt%nkpnt = 1
       nkpt%nkpts = 1
       parm%ibrav = -1
       parm%nr1 = 0
       parm%nr2 = 0
       parm%nr3 = 0
       prcp_com%akin = 0.0_real_8
       prcp_com%druck = 0.0_real_8
       prcp_com%eckin = 0.0_real_8
       prcp_com%gckin = 0.0_real_8
       prcp_com%skin = 0.0_real_8
       prcpl%tisot = .FALSE.
       prcpl%tzflex = .FALSE.
       ropt_mod%tesr = .FALSE.
       shock1%eshock  =  0._real_8
       shock1%pshock  =  0._real_8
       shock1%vol0  =  0._real_8
       shock1%vshock  =  0._real_8
       spin_mod%nspin = 0
       spin_mod%nsup = 0
       symmi%indpg = 0
       symmi%naxis = 0
       symmi%nrot = 0
       symmr%deltasym = 1.e-6_real_8
       symmt%tmsym = .FALSE.
       symmt%tpgauto = .FALSE.
       symmt%tsymm = .FALSE.
       tdgcomm%tdg = .FALSE.
       tkpts%tkall = .FALSE.
       tkpts%tkbcalc = tkpts%tknoswap ! TKNOSWAP is initialized in control.F (KOHN-SHAM EIGENVALUES)
       tkpts%tkblock = tkpts%tknoswap ! TKNOSWAP is initialized in control.F (KOHN-SHAM EIGENVALUES)
       tkpts%tkfull = .FALSE.
       tkpts%tkpnt = .FALSE.
       tkpts%tkscale = .FALSE.
       tkpts%tmonkp = .FALSE.
       tkpts%tonlydiag  = .FALSE.
       tkpts%tsymkp = .FALSE.
       wcfix%wcut = -1.0_real_8
       wgaussl%thdas = .FALSE.
       wgaussl%thdawm = .FALSE.
       wgaussl%twgauss = .FALSE.
       !
       CALL zeroing(cell_com%celldm)
       CALL zeroing(clas8%cellcl)
       CALL zeroing(prcp_com%cellrf)
       CALL zeroing(prcp_com%stens)
       CALL zeroing(symmr%ftau)
       CALL zeroing(wvk0)
       !
       ! Is the &SYSTEM section there?
       !
       ierr = inscan(iunit,'&SYSTEM')
       !
       IF (ierr == 0) THEN
          !
          ! Main loop
          !
          go_on_reading = .true.
          something_went_wrong = .false.
          !
          DO WHILE (go_on_reading)
             !
             previous_line = line
             READ(iunit,'(A80)',iostat=ierr) line
             IF (ierr /= 0) THEN
                something_went_wrong = .TRUE.
                go_on_reading        = .FALSE.
             ELSEIF ( keyword_contains(line,'&END') ) THEN
                go_on_reading = .FALSE.
             !
             ! cntl%cdft
             !
             ELSEIF ( keyword_contains(line,'DONOR') ) THEN
                first = INDEX(line,'DONOR') + 6
                CALL readsi(line,first,last,cdftpi%ndon,erread)
                IF (cdftpi%ndon > cdft_mat)  CALL stopgm(procedureN, 'Too many atoms in cdft donor group', &
                                                         __LINE__, __FILE__)
                IF ( keyword_contains(line,'WMULT') ) wgaussl%thdawm=.TRUE.
                IF (cdftpi%ndon > 0)THEN
                   previous_line = line
                   READ(iunit,'(A10000)',iostat=ierr) linelong
                   IF (ierr /= 0) THEN
                      error_message = 'ERROR WHILE READING LINE'
                      something_went_wrong = .true.
                      go_on_reading = .false.
                   ENDIF
                   first = 1
                ENDIF
                DO i=1,cdftpi%ndon
                   CALL readsi(linelong,first,last,cdftpi%cdft_d(i),erread)
                   IF (erread) THEN
                      error_message = 'ERROR WHILE READING VALUE'
                      something_went_wrong = .true.
                      go_on_reading = .false.
                   ENDIF
                   first = last 
                ENDDO
                IF (wgaussl%thdawm)THEN
                   previous_line = line
                   READ(iunit,'(A500)',iostat=ierr) linelong
                   IF (ierr /= 0) THEN
                      error_message = 'ERROR WHILE READING LINE'
                      something_went_wrong = .true.
                      go_on_reading = .false.
                   ENDIF
                   first = 1
                   DO i=1,cdftpi%ndon
                      CALL readsi(linelong,first,last,cdftpi%cdft_d(i+cdftpi%ndon),erread)
                      IF (erread) THEN
                         error_message = 'ERROR WHILE READING VALUE'
                         something_went_wrong = .true.
                         go_on_reading = .false.
                      ENDIF
                      first = last
                   ENDDO
                ENDIF
             ELSEIF ( keyword_contains(line,'ACCEPTOR') ) THEN
                first = INDEX(line,'ACCEPTOR') + 9
                CALL readsi(line,first,last,cdftpi%naccr,erread)
                IF (erread) THEN
                   error_message = 'ERROR WHILE READING VALUE'
                   something_went_wrong = .true.
                   go_on_reading = .false.
                ENDIF
                IF (cdftpi%naccr > cdft_mat) CALL stopgm(procedureN,'Too many atoms in cdft acceptor group',&
                                                         __LINE__, __FILE__)
                IF ( keyword_contains(line,'HDAS') ) wgaussl%thdas=.TRUE.
                IF ( keyword_contains(line,'WMULT') ) wgaussl%thdawm=.TRUE.
                IF (wgaussl%thdas.AND.wgaussl%thdawm) CALL stopgm(procedureN, 'HDASINGLE and WMULT are mutually exclusive', &
                                                                  __LINE__, __FILE__)
                previous_line=line
                READ(iunit,'(A10000)',iostat=ierr) linelong
                first = 1
                DO i=1,cdftpi%naccr
                   CALL readsi(linelong,first,last,cdftpi%cdft_a(i),erread)
                   IF (erread) THEN
                      error_message = 'ERROR WHILE READING VALUE'
                      something_went_wrong = .true.
                      go_on_reading = .false.
                   ENDIF
                   first = last 
                ENDDO
                IF (wgaussl%thdas)THEN
                   previous_line=line
                   READ(iunit,'(A10000)',iostat=ierr) linelong
                   IF (ierr /= 0) THEN
                      error_message = 'ERROR WHILE READING LINE'
                      something_went_wrong = .true.
                      go_on_reading = .false.
                   ENDIF
                   first = 1
                   DO i=1,cdftpi%naccr
                      CALL readsi(linelong,first,last,cdftpi%cdft_d(i),erread)
                      IF (erread) THEN
                         error_message = 'ERROR WHILE READING VALUE'
                         something_went_wrong = .true.
                         go_on_reading = .false.
                      ENDIF
                      first = last 
                   ENDDO
                ELSEIF (wgaussl%thdawm)THEN
                   previous_line=line
                   READ(iunit,'(A500)',iostat=ierr) linelong
                   IF (ierr /= 0) THEN
                      error_message = 'ERROR WHILE READING LINE'
                      something_went_wrong = .true.
                      go_on_reading = .false.
                   ENDIF
                   first = 1
                   DO i=1,cdftpi%naccr
                      CALL readsi(linelong,first,last,cdftpi%cdft_a(i+cdftpi%naccr),erread)
                      IF (erread) THEN
                         error_message = 'ERROR WHILE READING VALUE'
                         something_went_wrong = .true.
                         go_on_reading = .false.
                      ENDIF
                      first = last 
                   ENDDO
                ENDIF
             ELSEIF ( keyword_contains(line,'WGAUSS') ) THEN
                wgaussl%twgauss=.TRUE.
                first = INDEX(line,'WGAUSS') + 7
                CALL readsi(line, first,last,wg_n,erread)
                IF (erread) THEN
                   error_message = 'ERROR WHILE READING VALUE'
                   something_went_wrong = .true.
                   go_on_reading = .false.
                ENDIF
                DO i=1,wg_n
                   previous_line = line
                   READ(iunit,'(A80)',iostat=ierr) line
                   CALL readsr(line,1,last,wg_sigma(i),erread)
                   IF (erread) THEN
                      error_message = 'ERROR WHILE READING VALUE'
                      something_went_wrong = .true.
                      go_on_reading = .false.
                   ENDIF
                ENDDO
             ELSEIF ( keyword_contains(line,'WCUT') ) THEN
                first = INDEX(line,'WCUT') + 5
                CALL readsr(line, first,last,wcfix%wcut,erread)
                IF (erread) THEN
                   error_message = 'ERROR WHILE READING VALUE'
                   something_went_wrong = .true.
                   go_on_reading = .false.
                ENDIF
             !
             ! End of CDFT - start of normal SYSTEM input
             !
             ELSEIF ( keyword_contains(line,'CHECK',and='SYMMETRY') ) THEN
                IF ( keyword_contains(line,'OFF') ) THEN
                   boxeps=-1.0_real_8
                ENDIF
                previous_line=line
                READ(iunit,'(A80)',iostat=ierr) line
                CALL readsr(line,1,last,boxeps,erread)
                IF (erread) THEN
                   error_message = 'ERROR READING SYMMETRY PRECISION'
                   something_went_wrong = .TRUE.
                   go_on_reading = .FALSE.
                ENDIF
             !
             ! Unit cell symmetry
             !
             ELSEIF ( keyword_contains(line,'SYMMETRY') ) THEN
                previous_line = line
                READ(iunit,'(A80)',iostat=ierr) line
                IF ( keyword_contains(line,'ISOLATED') ) THEN
                   parm%ibrav=0
                ELSEIF ( keyword_contains(line,'CUBIC') ) THEN
                   IF ( keyword_contains(line,'FACE') ) THEN
                      parm%ibrav=2
                   ELSEIF ( keyword_contains(line,'BODY') ) THEN
                      parm%ibrav=3
                   ELSE
                      parm%ibrav=1
                   ENDIF
                ELSEIF ( keyword_contains(line,'FCC') ) THEN
                   parm%ibrav=2
                ELSEIF ( keyword_contains(line,'BCC') ) THEN
                   parm%ibrav=3
                ELSEIF ( keyword_contains(line,'HEXAGONAL') ) THEN
                   parm%ibrav=4
                ELSEIF( keyword_contains(line,'TRIGONAL',alias='RHOMBOHEDRAL') ) THEN
                   parm%ibrav=5
                ELSEIF ( keyword_contains(line,'TETRAGONAL') ) THEN
                   IF ( keyword_contains(line,'BODY') ) THEN
                      parm%ibrav=7
                   ELSE
                      parm%ibrav=6
                   ENDIF
                ELSEIF ( keyword_contains(line,'BCT') ) THEN
                   parm%ibrav=7
                ELSEIF ( keyword_contains(line,'ORTHORHOMBIC') ) THEN
                   parm%ibrav=8
                ELSEIF ( keyword_contains(line,'MONOCLINIC') ) THEN
                   parm%ibrav=12
                ELSEIF ( keyword_contains(line,'TRICLINIC') ) THEN
                   parm%ibrav=14
                ELSE
                   CALL readsi(line,1,last,parm%ibrav,erread)
                   IF (erread .or. (parm%ibrav < 0) .or. (parm%ibrav > 14)) THEN
                      error_message = 'ERROR READING SYMMETRY NUMBER'
                      something_went_wrong = .TRUE.
                      go_on_reading = .FALSE.
                   ENDIF
                ENDIF
             !
             ! Point group symmetry
             !
             ELSEIF( keyword_contains(line,'POINT',and='GROUP') ) THEN
                IF ( keyword_contains(line,'DELTA',cut_at='=') ) THEN
                   first = index_of_delimiter(line,'DELTA','=')
                   CALL readsr(line,first,last,symmr%deltasym,erread)
                   IF (erread) THEN
                      error_message = 'ERROR READING DELTA VALUE'
                      something_went_wrong = .TRUE.
                      go_on_reading = .FALSE.
                   ENDIF
                   IF (symmr%deltasym <= 0._real_8) symmr%deltasym=1.e-6_real_8
                ENDIF
                IF ( keyword_contains(line,'MOLECULE',alias='MOLECULAR') ) THEN
                   symmt%tmsym=.TRUE.
                   previous_line=line
                   READ(iunit,'(A)',iostat=ierr) line
                   IF (ierr /= 0) THEN
                      error_message = 'ERROR WHILE READING POINT GROUP INFORMATION'
                      something_went_wrong = .TRUE.
                      go_on_reading = .FALSE.
                   ENDIF
                   CALL xstring(line,first,last)
                   READ(line(first:),'(A3)',iostat=ierr) stag
                   IF (ierr /= 0) THEN
                      error_message = 'ERROR WHILE READING LINE'
                      something_went_wrong = .true.
                      go_on_reading = .false.
                   ENDIF
                   CALL readsr(line,last+1,last,xa,erread)
                   IF (erread) THEN
                      error_message = 'ERROR READING THE ORDER OF PRINCIPLE AXIS'
                      something_went_wrong = .TRUE.
                      go_on_reading = .FALSE.
                   ENDIF
                   symmi%naxis=NINT(xa)
                   symmi%indpg=100
                ELSE
                   previous_line=line
                   READ(iunit,'(A80)',iostat=ierr) line
                   IF (ierr /= 0) THEN
                      error_message = 'ERROR WHILE READING LINE'
                      something_went_wrong = .true.
                      go_on_reading = .false.
                   ENDIF
                   first = index_of_delimiter(line,'NAME','=')
                   IF (first > 0) THEN
                      CALL xstring(line(first+1:80),in1,in2)
                      in1=iequal+in1
                      in2=iequal+in2
                      DO i=32,1,-1
                         CALL xstring(pgrp(i),ip1,ip2)
                         IF (line(in1:in2) == pgrp(i)(ip1:ip2)) THEN
                            symmi%indpg=i
                         ENDIF
                      ENDDO
                      DO i=32,1,-1
                         CALL xstring(pgrd(i),ip1,ip2)
                         IF (line(in1:in2) == pgrd(i)(ip1:ip2)) THEN
                            symmi%indpg=i
                         ENDIF
                      ENDDO
                      CALL stopgm(procedureN,'Invalid point group',& 
                           __LINE__,__FILE__)
                   ELSEIF ( keyword_contains(line,'AUTO') ) THEN
                      symmt%tpgauto=.TRUE.
                   ELSE
                      CALL readsi(line,1,last,symmi%indpg,erread)
                      IF (erread) THEN
                         error_message = 'ERROR READING THE POINT GROUP NUMBER'
                         something_went_wrong = .true.
                         go_on_reading = .false.
                      ENDIF
                   ENDIF
                ENDIF
                !
             ELSEIF ( keyword_contains(line,'CELL') ) THEN
                !
                IF ( keyword_contains(line,'REFERENCE') ) THEN
                   IF ( keyword_contains(line,'VECTOR',alias='VECTORS') ) THEN
                      READ(iunit,*,iostat=ierr) cell_com%ar1, cell_com%ar2, cell_com%ar3
                      IF (ierr /= 0) THEN
                         error_message = 'ERROR WHILE READING LINE'
                         something_went_wrong = .true.
                         go_on_reading = .false.
                      ENDIF
                      lcell%tcellrefvec=.TRUE.
                   ELSE
                      READ(iunit,*,iostat=ierr) (prcp_com%cellrf(i),i=1,6)
                      IF (ierr /= 0) THEN
                         error_message = 'ERROR WHILE READING LINE'
                         something_went_wrong = .true.
                         go_on_reading = .false.
                      ENDIF
                      IF ( keyword_contains(line,'ABSOLUTE',alias='ABS') ) THEN
                         prcp_com%cellrf(2:3) = prcp_com%cellrf(2:3)/prcp_com%cellrf(1)
                      ENDIF
                   ENDIF
                   IF ( keyword_contains(line,'DEGREE',alias='DEGREES') ) THEN
                      prcp_com%cellrf(4:6) = COS(prcp_com%cellrf(4:6)*pi/180._real_8)
                   ENDIF
                ELSEIF ( keyword_contains(line,'CLASSIC') ) THEN
                   READ(iunit,*,iostat=ierr) (clas8%cellcl(i),i=1,6)
                   IF (ierr /= 0) THEN
                      error_message = 'ERROR WHILE READING LINE'
                      something_went_wrong = .true.
                      go_on_reading = .false.
                   ENDIF
                   IF ( keyword_contains(line,'ABSOLUTE',alias='ABS') ) THEN
                      clas8%cellcl(2:3) = clas8%cellcl(2:3)/clas8%cellcl(1)
                   ENDIF
                   IF ( keyword_contains(line,'DEGREE',alias='DEGREES') ) THEN
                      clas8%cellcl(4:6)=COS(clas8%cellcl(4:6)*pi/180._real_8)
                   ENDIF
                ELSEIF ( keyword_contains(line,'ISOTROPIC') ) THEN
                   prcpl%tisot=.TRUE.
                ELSEIF ( keyword_contains(line,'ZFLEXIBLE') ) THEN
                   prcpl%tzflex=.TRUE.
                ELSE
                   IF ( keyword_contains(line,'VECTOR',alias='VECTORS') ) THEN
                      READ(iunit,*,iostat=ierr) parm%a1,parm%a2,parm%a3
                      IF (ierr /= 0) THEN
                         error_message = 'ERROR WHILE READING LINE'
                         something_went_wrong = .true.
                         go_on_reading = .false.
                      ENDIF
                      lcell%tcellvectors=.TRUE.
                   ELSE
                      READ(iunit,*,iostat=ierr) (cell_com%celldm(i),i=1,6)
                      IF (ierr /= 0) THEN
                         error_message = 'ERROR WHILE READING LINE'
                         something_went_wrong = .true.
                         go_on_reading = .false.
                      ENDIF
                      IF ( keyword_contains(line,'ABSOLUTE',alias='ABS') ) THEN
                         cell_com%celldm(2:3)=cell_com%celldm(2:3)/cell_com%celldm(1)
                      ENDIF
                      IF ( keyword_contains(line,'DEGREE',alias='DEGREES') ) THEN
                         cell_com%celldm(4:6)=COS(cell_com%celldm(4:6)*pi/180._real_8)
                      ENDIF
                   ENDIF
                ENDIF
             !
             ELSEIF ( keyword_contains(line,'CHARGE') ) THEN
                READ(iunit,*,iostat=ierr) crge%charge
                IF (ierr /= 0) THEN
                   error_message = 'ERROR WHILE READING LINE'
                   something_went_wrong = .true.
                   go_on_reading = .false.
                ENDIF
             !
             ELSEIF ( keyword_contains(line,'LOW-SPIN',and='EXCITATION')  .OR. & 
                      keyword_contains(line,'SPIN',and='EXCITATION') ) THEN
                IF ( cntl%tshop ) THEN
                   ! McB     ... surface hopping stuff ...
                   WRITE (6,'(/A/A/)') ' WARNING: IGNORING LSE SPECIFICATIONS',&
                                       '          USING SURFACE HOPPING DEFAULTS INSTEAD'
                   ! McB     
                   lspin2%troks=.TRUE.
                   clsd%nlsd=4
                   clsd%nlsx=3
                ENDIF
                ! Low spin excitation
                IF ( keyword_contains(line,'ROOT') ) THEN
                   lspin2%troot=.TRUE.
                   clsd%nlsd=4
                   clsd%nlsx=3
                ELSEIF ( keyword_contains(line,'ROSS') ) THEN
                   lspin2%tross=.TRUE.
                   clsd%nlsd=5
                   clsd%nlsx=3
                ELSEIF ( keyword_contains(line,'CAS22') ) THEN
                   lspin2%tcas22=.TRUE.
                   clsd%nlsd=7
                   clsd%nlsx=3
                ELSE
                   ! McB_0705  ... default: RESTRICTED OPEN SHELL ...
                   lspin2%troks=.TRUE.
                   clsd%nlsd=4
                   clsd%nlsx=3
                ENDIF
                IF ( keyword_contains(line,'PENALTY') ) THEN
                   lspin2%tpenal=.TRUE.
                   READ(iunit,*,iostat=ierr) lspin4%apenal
                   clsd%nlsd=MAX(clsd%nlsd,5)
                   clsd%nlsx=3
                ENDIF
                IF ( keyword_contains(line,'LSETS') ) THEN
                   lspin2%tlsets=.TRUE.
                ENDIF
                lspin2%tlse=.TRUE.
                ! SMG    modified goedecker
             ELSEIF ( keyword_contains(line,'MODIFIED',and='GOEDECKER') ) THEN
                lspin2%tmgoe = .TRUE.
                IF ( keyword_contains(line,'PARAMETERS',alias='PARAMETER') ) THEN
                   READ(iunit,*,iostat=ierr) lspin3%mgmab,lspin3%mgmba
                   IF (ierr /= 0) THEN
                      error_message = 'ERROR WHILE READING LINE'
                      something_went_wrong = .true.
                      go_on_reading = .false.
                   ENDIF
                ENDIF
             ELSEIF ( keyword_contains(line,'LSE',and='PARAMETERS') .OR. & 
                      keyword_contains(line,'LSE',and='PARAMETER') ) THEN
                ! Parameters for low spin excitation
                READ(iunit,*,iostat=ierr) lspin1%lsea,lspin1%lseb
                IF (ierr /= 0) THEN
                   error_message = 'ERROR WHILE READING LINE'
                   something_went_wrong = .true.
                   go_on_reading = .false.
                ENDIF
             ELSEIF ( keyword_contains(line,'ENERGY', and='PROFILE') ) THEN
                ! Requests energy profile in lse calculation
                lspin2%teprof=.TRUE.
             ELSEIF ( keyword_contains(line,'COUPLING',alias='COUPLINGS') ) THEN
                tcpl=.TRUE.
                IF ( keyword_contains(line,'FD') ) THEN
                   tcplfd=.TRUE.
                   first = index_of_delimiter(line,'FD','=')
                   IF (first > 0) THEN
                      CALL readsr(line,first,last,eps_c,erread)
                      IF (erread) THEN
                         error_message = 'COUPLINGS: ERROR READING F.D. DISPLACEMENT'
                         something_went_wrong = .true.
                         go_on_reading = .false.
                      ENDIF
                   ENDIF
                ENDIF
                IF ( keyword_contains(line,'PROD') ) THEN
                   tcplfd=.TRUE.
                   first = index_of_delimiter(line,'PROD','=')
                   CALL readsr(line,first,last,eps_c,erread)
                   IF (erread) THEN
                      error_message = 'COUPLINGS: ERROR READING F.D. DISPLACEMENT'
                      something_went_wrong = .true.
                      go_on_reading = .false.
                   ENDIF
                   eps_c = -eps_c
                ENDIF
                IF ( keyword_contains(line,'LINRES') ) THEN
                   tcplfd=.FALSE.
                   tcpllr=.TRUE.
                   first = index_of_delimiter(line,'TOL','=')
                   IF (first > 0) THEN
                      CALL readsr(line,first,last,tolcpl,erread)
                      IF (erread) THEN
                         error_message = 'COUPLINGS: ERROR READING TOLERANCE'
                         something_went_wrong = .true.
                         go_on_reading = .false.
                      ENDIF
                   ENDIF
                   first = index_of_delimiter(line,'NVECT','=')
                   IF (first > 0) THEN
                      CALL readsi(line,first,last,nvect,erread)
                      IF (erread) THEN
                         error_message = 'COUPLINGS: ERROR READING NUMBER OF VECTORS'
                         something_went_wrong = .true.
                         go_on_reading = .false.
                      ENDIF
                      IF ( keyword_contains(line,'SPECIFY') ) tspecv = .TRUE.
                   ENDIF
                   IF ( keyword_contains(line,'THRESHOLDS',alias='THR') ) THEN
                      READ(iunit,*,iostat=ierr) fc_low, f_low
                      IF (ierr /= 0) THEN
                         error_message = 'ERROR WHILE READING LINE'
                         something_went_wrong = .true.
                         go_on_reading = .false.
                      ENDIF
                      READ(iunit,*,iostat=ierr) fc_med, f_med
                      IF (ierr /= 0) THEN
                         error_message = 'ERROR WHILE READING LINE'
                         something_went_wrong = .true.
                         go_on_reading = .false.
                      ENDIF
                      READ(iunit,*,iostat=ierr) fc_high, f_high
                      IF (ierr /= 0) THEN
                         error_message = 'ERROR WHILE READING LINE'
                         something_went_wrong = .true.
                         go_on_reading = .false.
                      ENDIF
                   ENDIF
                   IF ( keyword_contains(line,'BRUTE',and='FORCE') ) talldof = .TRUE.
                ENDIF
                IF ( keyword_contains(line,'NSURF') ) THEN
                   first = index_of_delimiter(line,'NSURF','=')
                   IF (first > 0) THEN
                      CALL readsi(line,first,last,nsurf,erread)
                      IF (erread) THEN
                         error_message = 'COUPLINGS: ERROR READING KS STATES'
                         something_went_wrong = .true.
                         go_on_reading = .false.
                      ENDIF
                      ALLOCATE(isurf(2,nsurf),STAT=ierr)
                      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                           __LINE__,__FILE__)
                      ALLOCATE(csurf(nsurf),STAT=ierr)
                      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                           __LINE__,__FILE__)
                      DO i = 1,nsurf
                         READ(iunit,*,iostat=ierr) isurf(1,i), isurf(2,i), csurf(i)
                      ENDDO
                   ENDIF
                ENDIF
                IF ( keyword_contains(line,'NAT') ) THEN
                   first = index_of_delimiter(line,'NAT','=')
                   IF (first > 0) THEN
                      CALL readsi(line,first,last,natfd,erread)
                      IF (erread) THEN
                         error_message = 'COUPLINGS: ERROR READING ATOMS FOR FD'
                         something_went_wrong = .true.
                         go_on_reading = .false.
                      ENDIF
                      tallat=.FALSE.
                      ALLOCATE(iatfd(natfd),STAT=ierr)
                      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                           __LINE__,__FILE__)
                      READ(iunit,*,iostat=ierr) (iatfd(i), i=1,natfd)
                   ENDIF
                ENDIF
             ELSEIF ( keyword_contains(line,'EXTERNAL',and='FIELD') ) THEN
                ! Applied external electric field
                cntl%tfield=.TRUE.
                READ(iunit,*,iostat=ierr) aifield
                IF (ierr /= 0) THEN
                   error_message = 'ERROR WHILE READING LINE'
                   something_went_wrong = .true.
                   go_on_reading = .false.
                ENDIF
             ELSEIF ( keyword_contains(line,'STATES',alias='STATE') ) THEN
                !CSOC[
                IF (cntl%tspec .OR. cntl%tsoc) THEN
                   CALL stopgm(procedureN,'No states with electronic spectra',& 
                        __LINE__,__FILE__)
                ENDIF
                !CSOC]
                ! Number of states to be optimized
                READ(iunit,*,iostat=ierr) crge%n
                IF (ierr /= 0) THEN
                   error_message = 'ERROR WHILE READING LINE'
                   something_went_wrong = .true.
                   go_on_reading = .false.
                ENDIF
                ! NKSSTA is the number of unoccupied states
                ! (parameter of KOHN-SHAM ENERGIES)
                IF (cntl%tlsd) THEN
                   in1=crge%n+2*cnti%nkssta
                ELSE
                   in1=crge%n+cnti%nkssta
                ENDIF
                ALLOCATE(crge%f(in1,1),STAT=ierr)
                IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                     __LINE__,__FILE__)
                DO i=1,in1
                   crge%f(i,1)=-1._real_8
                ENDDO
             ELSEIF ( keyword_contains(line,'OCCUPATION') ) THEN
                ! Occupation numbers
                IF (crge%n == 0) CALL stopgm(procedureN,&
                     'Specify the number of states before the occupation numbers',& 
                     __LINE__,__FILE__)
                IF ( keyword_contains(line,'FIXED') ) tfixo=.TRUE.
                READ(iunit,*,iostat=ierr) (crge%f(i,1),i=1,crge%n)
                IF (ierr /= 0) THEN
                   error_message = 'ERROR WHILE READING LINE'
                   something_went_wrong = .true.
                   go_on_reading = .false.
                ENDIF
             !
             ELSEIF ( keyword_contains(line,'MULTIPLICITY') ) THEN
                ! 2*S+1
                READ(iunit,*,iostat=ierr) spin_mod%nspin
                IF (ierr /= 0) THEN
                   error_message = 'ERROR WHILE READING LINE'
                   something_went_wrong = .true.
                   go_on_reading = .false.
                ENDIF
             ELSEIF ( keyword_contains(line,'NSUP') ) THEN
                ! Number of alpha spin states
                READ(iunit,*,iostat=ierr) spin_mod%nsup
                IF (ierr /= 0) THEN
                   error_message = 'ERROR WHILE READING LINE'
                   something_went_wrong = .true.
                   go_on_reading = .false.
                ENDIF
             !
             ELSEIF ( keyword_contains(line,'CUTOFF') ) THEN
                ! Cutoff for PW     
                IF ( keyword_contains(line,'CONSTANT') ) THEN
                   READ(iunit,*,iostat=ierr) prcp_com%akin,prcp_com%skin,prcp_com%eckin
                   IF (ierr /= 0) THEN
                      error_message = 'ERROR WHILE READING LINE'
                      something_went_wrong = .true.
                      go_on_reading = .false.
                   ENDIF
                ELSEIF ( keyword_contains(line,'HFX') ) THEN
                   READ(iunit,*,iostat=ierr) hfxc2%hfxwfe,hfxc2%hfxdee
                   IF (ierr /= 0) THEN
                      error_message = 'ERROR WHILE READING LINE'
                      something_went_wrong = .true.
                      go_on_reading = .false.
                   ENDIF
                ELSEIF ( keyword_contains(line,'DENSITY') ) THEN
                   IF ( keyword_contains(line,'NUMBER') ) THEN
                      READ(iunit,*,iostat=ierr) nhgwish
                      IF (ierr /= 0) THEN
                         error_message = 'ERROR WHILE READING LINE'
                         something_went_wrong = .true.
                         go_on_reading = .false.
                      ENDIF
                      tnhgwish=.TRUE.
                   ELSE
                      READ(iunit,*,iostat=ierr) ecutd
                      IF (ierr /= 0) THEN
                         error_message = 'ERROR WHILE READING LINE'
                         something_went_wrong = .true.
                         go_on_reading = .false.
                      ENDIF
                   ENDIF
                ELSE
                   IF ( keyword_contains(line,'NOSPHERICAL',alias='NONSPHERICAL') ) THEN
                      tsphere=.FALSE.
                   ELSEIF ( keyword_contains(line,'SPHERICAL',alias='SPHERIC') ) THEN
                      tsphere=.TRUE.
                   ENDIF
                   READ(iunit,*,iostat=ierr) cntr%ecut
                   IF (ierr /= 0) THEN
                      error_message = 'ERROR WHILE READING LINE'
                      something_went_wrong = .true.
                      go_on_reading = .false.
                   ENDIF
                ENDIF
             ELSEIF ( keyword_contains(line,'POISSON',alias='POISSON-SOLVER') ) THEN
                ! Poisson solver
                IF ( keyword_contains(line,'HOCKNEY') ) THEN
                   isos3%ps_type=1
                   IF ( keyword_contains(line,'PARAMETER',alias='PARAMETERS') ) THEN
                      READ(iunit,*,iostat=ierr) isos2%rcerf
                      IF (ierr /= 0) THEN
                         error_message = 'ERROR WHILE READING LINE'
                         something_went_wrong = .true.
                         go_on_reading = .false.
                      ENDIF
                   ENDIF
                ELSEIF ( keyword_contains(line,'TUCKERMAN',alias='TUCKERMAN-MARTYNA') ) THEN
                   isos3%ps_type=2
                   IF ( keyword_contains(line,'PARAMETER',alias='PARAMETERS') ) THEN
                      READ(iunit,*,iostat=ierr) isos2%alphal
                      IF (ierr /= 0) THEN
                         error_message = 'ERROR WHILE READING LINE'
                         something_went_wrong = .true.
                         go_on_reading = .false.
                      ENDIF
                   ENDIF
                ELSEIF ( keyword_contains(line,'MORTENSEN') ) THEN
                   isos3%ps_type=3
                !
                ! No more defaults!
                !
                ! ELSE
                !    isos3%ps_type=1
                ENDIF
             ELSEIF(( keyword_contains(line,'BOX') ).AND.&
                  ( keyword_contains(line,'WALLS') )) THEN
                READ(iunit,*,iostat=ierr) isos2%wall_skin
                IF (ierr /= 0) THEN
                   error_message = 'ERROR WHILE READING LINE'
                   something_went_wrong = .true.
                   go_on_reading = .false.
                ENDIF
                isos1%twall=.TRUE.
             ELSEIF ( keyword_contains(line,'DUAL') ) THEN
                ! ..Density Cutoff in units of wavefunction cutoff
                READ(iunit,*,iostat=ierr) dual00%cdual
                IF (ierr /= 0) THEN
                   error_message = 'ERROR WHILE READING LINE'
                   something_went_wrong = .true.
                   go_on_reading = .false.
                ENDIF
                dual00%dual=.TRUE.
             ELSEIF ( keyword_contains(line,'PRESSURE') ) THEN
                ! PRESSURE          
                READ(iunit,*,iostat=ierr) prcp_com%druck
                IF (ierr /= 0) THEN
                   error_message = 'ERROR WHILE READING LINE'
                   something_went_wrong = .true.
                   go_on_reading = .false.
                ENDIF
                tpressure=.TRUE.
             ELSEIF ( keyword_contains(line,'STRESS',and='TENSOR') ) THEN
                ! STRESS TENSOR     
                READ(iunit,*,iostat=ierr) (prcp_com%stens(i,1),i=1,9)
                IF (ierr /= 0) THEN
                   error_message = 'ERROR WHILE READING LINE'
                   something_went_wrong = .true.
                   go_on_reading = .false.
                ENDIF
                tstrtensor=.TRUE.
             ! 
             ! K-points
             !
             ELSEIF ( keyword_contains(line,'KPOINT',alias='KPOINTS') ) THEN
                tkpts%tkpnt=.TRUE.
                IF ( keyword_contains(line,'KDP') ) THEN
                   tkdp=.TRUE.
                   tkpts%tkpnt=.FALSE.
                   bmix=1._real_8
                ENDIF
                IF ( keyword_contains(line,'ONLYDIAG',alias='ONLYDIAGONAL') ) THEN
                   tkpts%tonlydiag = .TRUE.
                ENDIF
                IF ( keyword_contains(line,'BLOCK') ) THEN
                   tkpts%tkblock=.TRUE.
                   tkpts%tkbcalc=.TRUE.
                   IF (.NOT.cntl%tlanc) CALL stopgm(procedureN,'Swap files only with lanczos',& 
                                                    __LINE__,__FILE__)
                   IF (parai%cp_nproc > 1) CALL stopgm(procedureN,'No swap files with parallel',& 
                                                       __LINE__,__FILE__)
                   first = index_of_delimiter(line,'BLOCK','=')
                   CALL readsi(line,first,last,nkpt%nkpnt,erread)
                   IF (erread) THEN
                      error_message = 'KPOINTS, BLOCK OPTION: ERROR READING THE NUMBER'
                      something_went_wrong = .true.
                      go_on_reading = .false.
                   ENDIF
                   IF (nkpt%nkpnt <= 0) CALL stopgm(procedureN,'Wrong number of k-points block',& 
                                                    __LINE__,__FILE__)
                   IF ( keyword_contains(line,'ALL') ) THEN
                      tkpts%tkall=.TRUE.
                      tkpts%tkbcalc=.FALSE.
                   ELSEIF ( keyword_contains(line,'CALCULATED') ) THEN
                      tkpts%tkall=.FALSE.
                      tkpts%tkbcalc=.TRUE.
                   ENDIF
                   IF ( keyword_contains(line,'NOSWAP') ) THEN
                      tkpts%tknoswap=.TRUE.
                   ENDIF
                ENDIF
                IF ( keyword_contains(line,'MONKHORST',alias='MONKHORST-PACK') ) THEN
                   ! Monkhorst-pack k-points
                   tkpts%tmonkp=.TRUE.
                   IF ( keyword_contains(line,'SYMMETRIZE',alias='SYMMETRIZED') ) tkpts%tsymkp=.TRUE.
                   IF ( keyword_contains(line,'FULL') ) tkpts%tkfull=.TRUE.
                   ! Read NK1, NK2, and NK3
                   READ(iunit,'(A)',iostat=ierr) line
                   IF (ierr /= 0) THEN
                      error_message = 'ERROR WHILE READING LINE'
                      something_went_wrong = .true.
                      go_on_reading = .false.
                   ENDIF
                   first = 1
                   CALL readsi(line,first,last,kpts_com%nk1,erread)
                   IF (erread) THEN
                      error_message = 'KPOINTS, MONKHORST-PACK OPTION: ERROR READING NK1'
                      something_went_wrong = .true.
                      go_on_reading = .false.
                   ENDIF
                   first = last 
                   CALL readsi(line,first,last,kpts_com%nk2,erread)
                   IF (erread) THEN
                      error_message = 'KPOINTS, MONKHORST-PACK OPTION: ERROR READING NK2'
                      something_went_wrong = .true.
                      go_on_reading = .false.
                   ENDIF
                   first = last 
                   CALL readsi(line,first,last,kpts_com%nk3,erread)
                   IF (erread) THEN
                      error_message = 'KPOINTS, MONKHORST-PACK OPTION: ERROR READING NK3'
                      something_went_wrong = .true.
                      go_on_reading = .false.
                   ENDIF
                   IF (kpts_com%nk1 <= 0 .OR. kpts_com%nk2 <= 0 .OR. kpts_com%nk3 <= 0) THEN
                      WRITE(output_unit,'(T3,"NK1=",I5," NK2=",I5," NK3=",I5)') kpts_com%nk1,kpts_com%nk2,kpts_com%nk3
                      CALL stopgm(procedureN, 'Bad monkhorst-pack mesh value, cf. output file',& 
                           __LINE__,__FILE__)
                   ENDIF
                   ! MacDonald shift
                   IF ( keyword_contains(line,'SHIFT',cut_at='=') ) THEN
                      first = index_of_delimiter(line,'SHIFT','=')
                      DO i=1,3
                         CALL readsr(line,first,last,wvk0(i),erread)
                         IF (erread) THEN
                            error_message = 'KPOINTS, ERROR READING MacDONALD SHIFT'
                            something_went_wrong = .true.
                            go_on_reading = .false.
                         ENDIF
                         first = last 
                      ENDDO
                   ENDIF
                ELSEIF ( keyword_contains(line,'BAND',alias='BANDS') ) THEN
                   IF ( keyword_contains(line,'SCALE',alias='SCALED') ) tkpts%tkscale=.TRUE.
                   nkpbd=1
                   nkpt%nkpts=0
                   DO WHILE (nkpbd /= 0)
                      READ(iunit,*,iostat=ierr) nkpbd, (rki(i),i=1,3), (rkf(i),i=1,3)
                      IF (ierr /= 0) THEN
                         error_message = 'ERROR WHILE READING LINE'
                         something_went_wrong = .true.
                         go_on_reading = .false.
                      ENDIF
                      IF (nkpbd /= 0) THEN
                         IF (nkpbd <= 1) CALL stopgm(procedureN,'Wrong number of k-points per band', &
                                                     __LINE__,__FILE__)
                         ALLOCATE(rkt(3,nkpt%nkpts+nkpbd),STAT=ierr)
                         IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                                                 __LINE__,__FILE__)
                         IF (nkpt%nkpts /= 0) CALL dswap(3*nkpt%nkpts,rk,1,rkt,1)
                         DO ik=1,nkpbd
                            DO i=1,3
                               rkt(i,nkpt%nkpts+ik)=rki(i)&
                                    +REAL(ik-1,kind=real_8)*(rkf(i)-rki(i))/REAL(nkpbd-1,kind=real_8)
                            ENDDO
                         ENDDO
                         IF (nkpt%nkpts /= 0) DEALLOCATE(rk,STAT=ierr)
                         IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
                              __LINE__,__FILE__)
                         nkpt%nkpts=nkpt%nkpts+nkpbd
                         ALLOCATE(rk(3,nkpt%nkpts),STAT=ierr)
                         IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                              __LINE__,__FILE__)
                         CALL dswap(3*nkpt%nkpts,rkt,1,rk,1)
                         DEALLOCATE(rkt,STAT=ierr)
                         IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
                              __LINE__,__FILE__)
                      ENDIF
                   ENDDO
                   IF (nkpt%nkpts <= 0)&
                        CALL stopgm('  SYSIN','WRONG NUMBER OF K POINTS',& 
                        __LINE__,__FILE__)
                   ALLOCATE(wk(nkpt%nkpts),STAT=ierr)
                   IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                        __LINE__,__FILE__)
                   DO i=1,nkpt%nkpts
                      wk(i)=1._real_8
                   ENDDO
                ELSE
                   IF ( keyword_contains(line,'SCALE',alias='SCALED') ) tkpts%tkscale=.TRUE.
                   READ(iunit,*,iostat=ierr) nkpt%nkpts
                   IF (ierr /= 0) THEN
                      error_message = 'ERROR WHILE READING LINE'
                      something_went_wrong = .true.
                      go_on_reading = .false.
                   ENDIF
                   IF (nkpt%nkpts <= 0)&
                      CALL stopgm(procedureN,'Wrong number of k-points',& 
                                    __LINE__,__FILE__)
                   ALLOCATE(rk(3,nkpt%nkpts),STAT=ierr)
                   IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                        __LINE__,__FILE__)
                   ALLOCATE(wk(nkpt%nkpts),STAT=ierr)
                   IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                        __LINE__,__FILE__)
                   DO ik=1,nkpt%nkpts
                      READ(iunit,'(A)',iostat=ierr) line
                      IF (ierr /= 0) THEN
                         error_message = 'ERROR WHILE READING LINE'
                         something_went_wrong = .true.
                         go_on_reading = .false.
                      ENDIF
                      first = 1
                      DO i=1,3
                         CALL readsr(line,first,last,rk(i,ik),erread)
                         IF (erread) THEN
                            error_message = 'ERROR WHILE READING VALUE' 
                            something_went_wrong = .true.
                            go_on_reading = .false.
                         ENDIF
                         first = last 
                      ENDDO
                      CALL readsr(line,first,last,wk(ik),erread)
                      IF (erread) THEN
                         error_message = 'ERROR WHILE READING VALUE' 
                         something_went_wrong = .true.
                         go_on_reading = .false.
                      ENDIF
                   ENDDO
                ENDIF
             !
             ! End of k-points
             !
             ELSEIF( keyword_contains(line,'MESH',but_not='MONKHOSRT') ) THEN
                ! Real space mesh for wavefunction
                     READ(iunit,*,iostat=ierr) parm%nr1,parm%nr2,parm%nr3
             ELSEIF ( keyword_contains(line,'SCALE') ) THEN
                ! Scale atomic coordinates with unit cell dimensions
                cntl%tscale=.TRUE.
                sxscale=1._real_8
                syscale=1._real_8
                szscale=1._real_8
                IF ( keyword_contains(line,'CARTESIAN') ) THEN
                   tcartesian=.TRUE.
                ELSE
                   tcartesian=.FALSE.
                ENDIF
                IF ( keyword_contains(line,'S') ) THEN
                   first = index_of_delimiter(line,'S','=')
                   CALL readsr(line,first,last,sxscale,erread)
                   IF (erread .OR. (sxscale <= 0._real_8)) THEN
                      error_message = 'SCALE OPTION: ERROR READING THE NUMBER  S='
                      something_went_wrong = .true.
                      go_on_reading = .false.
                   ENDIF
                   syscale=sxscale
                   szscale=sxscale
                ENDIF
                IF ( keyword_contains(line,'SX') ) THEN
                   first = index_of_delimiter(line,'SX','=')
                   CALL readsr(line,first,last,sxscale,erread)
                   IF (erread .OR. (sxscale <= 0._real_8)) THEN
                      error_message = 'SCALE OPTION: ERROR READING THE NUMBER SX='
                      something_went_wrong = .true.
                      go_on_reading = .false.
                   ENDIF
                ENDIF
                IF ( keyword_contains(line,'SY') ) THEN
                   first = index_of_delimiter(line,'SY','=')
                   CALL readsr(line,first,last,syscale,erread)
                   IF (erread .OR. (syscale <= 0._real_8)) THEN
                      error_message = 'SCALE OPTION: ERROR READING THE NUMBER SY='
                      something_went_wrong = .true.
                      go_on_reading = .false.
                   ENDIF
                ENDIF
                IF ( keyword_contains(line,'SZ') ) THEN
                   first = index_of_delimiter(line,'SZ','=')
                   CALL readsr(line,first,last,szscale,erread)
                   IF (erread .OR. (szscale <= 0._real_8)) THEN
                      error_message = 'SCALE OPTION: ERROR READING THE NUMBER SZ='
                      something_went_wrong = .true.
                      go_on_reading = .false.
                   ENDIF
                ENDIF
             ELSEIF( keyword_contains(line,'DOUBLE',and='GRID') ) THEN
                ! Double Grid
                IF ( keyword_contains(line,'ON') ) tdgcomm%tdg=.TRUE.
                IF ( keyword_contains(line,'OFF') ) tdgcomm%tdg=.FALSE.
             ELSEIF( keyword_contains(line,'SYMMETRIZE',and='COORDINATES') ) THEN
                ! Symmetrize atomic coordinates
                symmt%tsymm=.TRUE.
             ELSEIF ( keyword_contains(line,'ANGSTROM',alias='ANGSTROMS') ) THEN
                ! Units of input are angstroms
                IF (cntl%tqmmm) THEN
                   CALL stopgm(procedureN,'ANGSTROM not available with QM/MM',& 
                               __LINE__,__FILE__)
                ENDIF
                cntl%bohr=.FALSE.
             ELSEIF ( keyword_contains(line,'TESR') ) THEN
                ! Extended real space integration for image charges
                ropt_mod%tesr=.TRUE.
                READ(iunit,*,iostat=ierr) iteropt%iesr
             ELSEIF ( keyword_contains(line,'SURFACE') ) THEN
                ! System is periodic in two dimensions [x and y in default]
                isos1%tclust=.TRUE.
                isos1%ttwod=.TRUE.
                IF ( keyword_contains(line,'XY') ) isos3%snormal=3
                IF ( keyword_contains(line,'YZ') ) isos3%snormal=1
                IF ( keyword_contains(line,'ZX') ) isos3%snormal=2
             ELSEIF ( keyword_contains(line,'POLYMER') ) THEN
                ! System is periodic in one dimension [x]
                isos1%tclust=.TRUE.
                isos1%toned=.TRUE.
             ELSEIF ( keyword_contains(line,'CLUSTER') ) THEN
                ! System is not periodic
                isos1%tclust=.TRUE.
             ELSEIF ( keyword_contains(line,'SHOCK',and='VELOCITY') ) THEN
                ! Velocity of shock wave in m/sec
                READ(iunit,*,iostat=ierr) shock1%vshock
                ! Convert from m/sec -> a/fs -> cntl%bohr/fs -> au
                shock1%vshock = shock1%vshock * 1.e-5_real_8 * fbohr * au_fs
             ELSE
                ! Unknown Keyword. store and continue
                IF (' ' /= line) THEN
                   IF (nbr_unknown_lines < max_unknown_lines) THEN
                      nbr_unknown_lines=nbr_unknown_lines+1
                      unknown(nbr_unknown_lines)=line
                   ELSE
                      DO i=1,max_unknown_lines-1
                         unknown(i)=unknown(i+1)
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
          error_message        = 'MISSING &SYSTEM SECTION - SECTION MANDATORY'
       ENDIF
       !
       IF (something_went_wrong) THEN
          WRITE(output_unit,'(/,1X,64("!"))')
          WRITE(output_unit,'(1X, A, 1X, A)') 'ERROR:', 'PROBLEM WHILE READING &SYSTEM SECTION:'
          WRITE(output_unit,'(8X, A)') TRIM(ADJUSTL(error_message))
          IF (line /= ' ' .or. previous_line /= ' ') THEN
             WRITE(output_unit,'(8X, A)') 'THE LAST TWO LINES READ WITHIN THE SECTION WERE:'
             WRITE(output_unit,'(/,1X, A)') TRIM(ADJUSTL(previous_line))
             WRITE(output_unit,'(1X, A)') TRIM(ADJUSTL(line))
          END IF
          WRITE(output_unit,'(1X,64("!"))')
          CALL stopgm(procedureN,'Error while reading &SYSTEM section, cf. output file',&
               __LINE__,__FILE__)
       ENDIF
       !
       !
       ! EHR[
       IF (cntl%cmplx_wf) THEN
          tkpts%tkpnt=.TRUE.
          nkpt%nkpts=1
          ALLOCATE(rk(3,nkpt%nkpts),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(wk(nkpt%nkpts),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          rk(1,1)=0._real_8
          rk(2,1)=0._real_8
          rk(3,1)=0._real_8
          wk(1)=1._real_8
          ! all bcast in setsys.F
       ENDIF
       ! EHR]
       !
       !
       CALL check_options()

       IF (nbr_unknown_lines /= 0) THEN
          WRITE(output_unit,'(/,1X,64("="))')
          WRITE(output_unit,'(1X,A,14X,A,12X,A)') '= ', 'UNKNOWN KEYWORDS IN SECTION &SYSTEM','='
          DO i=1,nbr_unknown_lines
             previous_line=unknown(i)
             CALL xstring(previous_line,first,last)
             WRITE(output_unit,'(A,A)') ' = ',previous_line(first:last)
          ENDDO
          WRITE(output_unit,'(1X,64("="),/)')
       ENDIF

    ENDIF

    CALL broadcast_sysin()


    CONTAINS

    ! ==--------------------------------------------------------------==
    SUBROUTINE check_options()

       ! mb   Prevent the use of non-isolated cluster options in QMMM
       IF (cntl%tqmmm) THEN
          IF (parm%ibrav /= 0) THEN
             WRITE(output_unit,'(A)') ' ERROR: SYMMETRY IBRAV =',parm%ibrav,' NOT ALLOWED IN QMMM'
             WRITE(output_unit,'(A)') '        A VALUE OF 0 (ISOLATED) IS MANDATORY'
             CALL stopgm(procedureN,'Cell dimensions incompatible with QM/MM', &
                         __LINE__, __FILE__)
          ENDIF
          !
          ! MPB: This warning may only cause confusion, since the use of "CLUSTER"
          !      is discouraged in the manual and not mentioned a single time in the QM/MM
          !      manual - and iv param%ibrav /= 0, isos1%tclust is set to true later in this
          !      routine anyway (!)
          !
          ! IF (.NOT. isos1%tclust) THEN! this is redundant but safe
          !    WRITE(output_unit,'(A)') ' WARNING: ISOLATED CLUSTER OPTION NEEDED IN QMMM'
          !    WRITE(output_unit,'(A)') '          NOW TCLUST WILL BE SET TO .TRUE.'
          !    parm%ibrav = 0
          !    isos1%tclust = .TRUE.
          ! ENDIF
          !
       ENDIF

       ! 
       ! AK: sanity check. ibrav=-2 is used for implicit triclinic symmetry
       ! with given cell vectors. symmetry should not be used in that case.
       !
       IF (lcell%tcellvectors) THEN
          IF (parm%ibrav /= -1) THEN
             CALL stopgm(procedureN,'SYMMETRY cannot be used with CELL VECTORS',& 
                  __LINE__,__FILE__)
          ELSE
             parm%ibrav = -2
          ENDIF
       ENDIF
       !
       IF (parm%ibrav == -1) THEN
          CALL stopgm(procedureN,'No symmetry specified.',& 
                      __LINE__,__FILE__)
       ENDIF
       IF (.NOT.lcell%tcellvectors .AND. cell_com%celldm(1) < 1.e-5_real_8) THEN
          CALL stopgm(procedureN,'Supercell dimensions are not specified; or the lattice constant is nought.',& 
               __LINE__,__FILE__)
       ENDIF

       ! AK: hmmm. parrinello-rahman or steepest descend cell need
       ! triclinic supercell symmetry (even for isotropic cell)
       ! and will reset ibrav. better do it here and tell people
       ! about it now, so that the output is correct.
       !
       ! MPB: Changed warning to stopgm
       IF ((cntl%tprcp .OR. cntl%tsdc).AND.(parm%ibrav /= 14)&
            .AND.(.NOT. lcell%tcellvectors)) THEN
          WRITE(output_unit,*) 'ERROR: USE SYMMETRY 14 (TRICLINIC) WITH VARIABLE CELLS'
          CALL stopgm(procedurEN,'Variable cell requires triclinic symmetry (applies to isotropic cells, too.)', &
                      __LINE__, __FILE__ )
       ENDIF

       ! AK: Ceck whether cell dimensions in input and symmetry do match
       CALL chkcsymm(boxeps)

       !
       ! LSE
       !
       IF (lspin2%tlse)THEN
          IF (cntl%tlsd) THEN
             CALL stopgm(procedureN,'LSDA not available with LSE option (NYI)',& 
                  __LINE__,__FILE__)
          ENDIF
          IF (cntl%tdiag) THEN
             CALL stopgm(procedureN,'Diagonalization not available with LSE (NYI)',& 
                  __LINE__,__FILE__)
          ENDIF
          IF (cntl%tpres) THEN
             CALL stopgm(procedureN,'Stress not available with LSE (NYI) ',& 
                  __LINE__,__FILE__)
          ENDIF
       ENDIF
       IF (.NOT.cntl%tlsd .AND. spin_mod%nspin /= 0) THEN
          WRITE(output_unit,*) ' MULTIPLICITY ',spin_mod%nspin,' INCOMPATIBLE WITH LDA'
          CALL stopgm(procedureN,'Multiplicity requires local spin density approximation (LSDA)',& 
               __LINE__,__FILE__)
       ENDIF

       !
       ! CB Some checks for broken symmetry
       !
       IF (cntl%bsymm) THEN
          IF (crge%n /= 0 .OR. spin_mod%nsup /= 0) THEN
             WRITE(output_unit,*) ' ERROR: FOR BROKEN SYMMETRY ARBITRARY STATE NUMBERS &
                         &ARE NOT SUPPORTED.'
             CALL stopgm(procedureN,'No arbitrary state numbers with broken symmetry',& 
                  __LINE__,__FILE__)
          ENDIF
          IF (spin_mod%nspin == 0) THEN
             WRITE(output_unit,*) ' FOR BROKEN SYMMETRY THE MULTIPLICITY &
                         &OF THE HIGHSPIN WAVEFUNCTION IS NEEDED.'
             CALL stopgm(procedureN,'No spin-restricted systems with broken symmetry',& 
                  __LINE__,__FILE__)
          ENDIF
          IF (MOD(spin_mod%nspin,2) == 0) THEN
             WRITE(output_unit,*) ' FOR BROKEN SYMMETRY THE MULTIPLICITY &
                         &OF THE HIGHSPIN WAVEFUNCTION MUST BE ODD.'
             CALL stopgm(procedureN,'No even spins with broken symmetry',&
                  __LINE__,__FILE__)
          ENDIF
       ENDIF
       IF (isos1%toned.AND.isos1%ttwod) THEN
          WRITE(output_unit,*) ' PERIODICITY IN ONE OR TWO DIMENSIONS?'
          CALL stopgm(procedureN,'Conflicting keywords, cf. output file',& 
               __LINE__,__FILE__)
       ENDIF
       IF (symmt%tpgauto.AND.symmt%tsymm) THEN
          WRITE(output_unit,*)' IMPOSSIBLE TO DETERMINE POINT GROUP AND &
                     &TO SYMMETRIZE COORDINATES'
          CALL stopgm(procedureN,'Conflicting keywords, cf. output file',& 
               __LINE__,__FILE__)
       ENDIF
       IF (tfixo.AND..NOT.cntl%tdiag) THEN
          WRITE(output_unit,*) ' FIXED OCCUPATION NUMBERS ONLY IN COMBINATION WITH &
                      &DIAGONALIZATION'
          CALL stopgm(procedureN,'Conflicting keywords, cf. output file',& 
               __LINE__,__FILE__)
       ENDIF
       IF (tfixo.AND.tkpts%tkpnt) THEN
          WRITE(output_unit,*) ' FIXED OCCUPATION NUMBERS IN COMBINATION WITH &
                      &KPOINTS NOT IMPLEMENTED'
          CALL stopgm(procedureN,'Conflicting keywords, cf. output file',& 
               __LINE__,__FILE__)
       ENDIF
       IF (.NOT.lspin2%tlse.AND.lspin2%tmgoe) THEN
          WRITE(output_unit,*) ' MODIFIED GOEDECKER ONLY AVAILABLE FOR ROKS '
          CALL stopgm(procedureN,'Conflicting keywords, cf. output file ',& 
               __LINE__,__FILE__)
       ENDIF
       !
       ! Overwriting...
       !
       IF (cntl%tshop) THEN
          ! ND     cntl%tlsd=.FALSE.
          lspin2%tlse = .TRUE.
       ENDIF
       !
       ! Check unsupported symmetries
       !
       IF (parm%ibrav == 9) THEN
          CALL stopgm(procedureN,'Symmetry type not supported',& 
               __LINE__,__FILE__)
       ELSEIF (parm%ibrav == 10) THEN
          CALL stopgm(procedureN,'Symmetry type not supported',& 
               __LINE__,__FILE__)
       ELSEIF (parm%ibrav == 11) THEN
          CALL stopgm(procedureN,'Symmetry type not supported',& 
               __LINE__,__FILE__)
       ELSEIF (parm%ibrav == 13) THEN
          CALL stopgm(procedureN,'Symmetry type not supported',& 
               __LINE__,__FILE__)
       ELSEIF (parm%ibrav < -2) THEN
          CALL stopgm(procedureN,'Illegal symmetry value',& 
               __LINE__,__FILE__)
       ELSEIF (parm%ibrav > 14) THEN
          CALL stopgm(procedureN,'Illegal symmetry value',& 
               __LINE__,__FILE__)
       ENDIF
       ! ==--------------------------------------------------------------==
       IF (cntl%tprcp .OR. cntl%tsdc) THEN
          IF (parm%ibrav == 0) THEN
             CALL stopgm(procedureN,'Variable cell and SYMMETRY 0 are incompatible ',& 
                  __LINE__,__FILE__)
          ENDIF
          IF (tstrtensor.AND.tpressure) THEN
             CALL stopgm(procedureN,'Stress tensor and pressure may not be specified simultaneously',& 
                  __LINE__,__FILE__)
          ENDIF
          IF (tpressure) THEN
             CALL zeroing(prcp_com%stens)
             DO i=1,3
                prcp_com%stens(i,i)=prcp_com%druck
             ENDDO
          ENDIF
          prcp_com%druck=prcp_com%druck*kb_au
          shock1%pshock=prcp_com%druck
          CALL dscal(9,kb_au,prcp_com%stens,1)
       ENDIF
       !
       ! Set symmetry to simple cubic or orthorhombic for isolated systems
       !
       IF (parm%ibrav == 0) THEN
          isos1%tclust=.TRUE.
          !
          ! qm/mm, surface, polymer need
          ! isos1%tclust, but are not isolated
          IF (.NOT.(cntl%tqmmm .OR. isos1%ttwod .OR. isos1%toned)) isos1%tisos=.TRUE.
          !
          ! Test if we are cubic.
          !
          IF ((ABS(cell_com%celldm(2)-1.0_real_8) + ABS(cell_com%celldm(3) - 1._real_8)) < 1.e-5_real_8) THEN
             parm%ibrav         = 1
             cell_com%celldm(2) = 1.0_real_8
             cell_com%celldm(3) = 1.0_real_8
             cell_com%celldm(4) = 0.0_real_8
             cell_com%celldm(5) = 0.0_real_8
             cell_com%celldm(6) = 0.0_real_8
          ELSE
             parm%ibrav         = 8
          ENDIF
       ENDIF
       !
       ! Default cell dimensions
       !
       IF (parm%ibrav == 1) THEN
          cell_com%celldm(2:3) = 1.0_real_8
          cell_com%celldm(4:6) = 0.0_real_8
       ELSEIF (parm%ibrav == 6) THEN
          cell_com%celldm(2)   = 1.0_real_8
          cell_com%celldm(4:6) = 0.0_real_8
       ELSEIF (parm%ibrav == 8) THEN
          cell_com%celldm(4:6) = 0.0_real_8
       ENDIF
       !
       IF ( (cell_com%celldm(4) < -1._real_8 .OR. cell_com%celldm(4) > 1._real_8) .OR. &
            (cell_com%celldm(5) < -1._real_8 .OR. cell_com%celldm(5) > 1._real_8) .OR. &
            (cell_com%celldm(6) < -1._real_8 .OR. cell_com%celldm(6) > 1._real_8) ) THEN
          CALL stopgm(procedureN,'Wrong angle cosines',& 
               __LINE__,__FILE__)
       ENDIF
       ! 
       IF (isos1%tclust .AND. isos3%ps_type == 0) THEN
          !
          ! No more defaults!
          !
          ! IF (isos1%toned .OR. isos1%ttwod) THEN
          !    isos3%ps_type = 3
          ! ELSE
          !    isos3%ps_type = 1
          ! ENDIF
          CALL stopgm(procedureN,'A POISSON SOLVER has to be specified for an ISOLATED/SURFACE/POLYMER system', &
                                 __LINE__, __FILE__)
       ENDIF
       !
       !
       IF (isos1%tclust) THEN
          IF (isos3%ps_type == 1) THEN
             IF (isos1%toned) THEN
                WRITE(output_unit,*) 'WARNING! POLYMER OPTION NOT TESTED WITH HOCKNEY POISSON SOLVER'
             ELSEIF (isos1%ttwod) THEN
                WRITE(output_unit,*) 'SURFACE WITH HOCKNEY POISSON SOLVER: &
                                      &INFLUENCE FUNCTION IN RECIPROCAL SPACE'
             ENDIF
          ELSEIF (isos3%ps_type == 2) THEN
             IF (isos1%toned .OR. isos1%ttwod) THEN
                CALL stopgm(procedureN,'SURFACE and POLYMER options are not available &
                                       &with the Tuckerman-Martyna poisson solver.',& 
                     __LINE__,__FILE__)
             ENDIF
          ELSEIF (isos3%ps_type == 3) THEN
             IF (isos1%toned) THEN
                IF (ABS(cell_com%celldm(2) - cell_com%celldm(3)) > 1.e-10_real_8) THEN
                   CALL stopgm(procedureN,'Cell dimensions b and c have to be equal &
                                           &for the Mortensen 1D Poisson solver',& 
                        __LINE__,__FILE__)
                ENDIF
             ELSEIF (isos1%ttwod) THEN
             ELSE
                IF (ABS(cell_com%celldm(2)*cell_com%celldm(3)-1._real_8) > 1.e-10_real_8) THEN
                        WRITE(output_unit,*) 'Only cubic boxes allowed for Mortensen',&
                        ' 0d(=ISOLATED SYSTEM) Poisson solver'
                   CALL stopgm(procedureN,'Mortensen 1D Poisson solver allows only for &
                                           &cubic boxes',&
                        __LINE__,__FILE__)
                ENDIF
             ENDIF
          ENDIF
       ENDIF

       !
       IF (lspin2%tcas22.AND.lspin2%tpenal) CALL stopgm("SYSIN","INCOMPATIBLE OPTIONS",& 
            __LINE__,__FILE__)
   
       !
       ! Double grid and cutoff
       !
       IF (cntr%ecut <= 0.0_real_8) THEN
          CALL stopgm(procedureN,'No plane wave cutoff specified (keyword CUTOFF) ',& 
               __LINE__,__FILE__)
       ENDIF
       !
       edgcomm%ecutwdg=cntr%ecut
       IF (tdgcomm%tdg) THEN
          edgcomm%ecutdg = 4._real_8*cntr%ecut
          IF (.NOT.dual00%dual) THEN
             IF (ecutd > 0._real_8) THEN
                dual00%cdual = ecutd/cntr%ecut
                dual00%dual  = .TRUE.
             ENDIF
          ELSE
             IF ( ecutd > 0._real_8 .and. &
                 (ABS(dual00%cdual-ecutd/cntr%ecut) > 1.e-5_real_8) ) CALL stopgm(procedureN, &
                        'Inconsistent dual/density cutoff',& 
                        __LINE__,__FILE__)
          ENDIF
          IF (dual00%cdual < 4._real_8) CALL stopgm(procedureN,'Inconsistent density cutoff',& 
                  __LINE__,__FILE__)
       ELSE
          IF (.NOT. dual00%dual) THEN
             IF (ecutd > 0._real_8) THEN
                dual00%cdual = ecutd/cntr%ecut
                dual00%dual  = .TRUE.
             ENDIF
          ELSE
             IF ( ecutd > 0._real_8 .and. &
                 (ABS(dual00%cdual-ecutd/cntr%ecut) > 1.e-5_real_8) ) CALL stopgm(procedureN, &
                        'Inconsistent dual/density cutoff',& 
                        __LINE__,__FILE__)
          ENDIF
          edgcomm%ecutdg=dual00%cdual*cntr%ecut
       ENDIF
       !
       IF (cdftlog%thda)THEN
          IF (cdftpi%ndon == 0 .AND. .NOT. wgaussl%thdas)THEN
             CALL stopgm(procedureN,'No donor atoms in HDA calculation; set HDAS flag.',& 
                  __LINE__,__FILE__)
          ENDIF
       ENDIF
       !
       ! External electric field
       ! 
       IF (cntl%tfield) THEN 
          WRITE(output_unit,'(A,T25,3E12.4,A,/)') ' APPLY ELECTRIC FIELD: ', (AIFIELD(I),I=1,3),' A.U.'
       ENDIF

    END SUBROUTINE check_options
    ! ==--------------------------------------------------------------==
    SUBROUTINE broadcast_sysin()

       ! ..DUAL
       CALL mp_bcast_byte(dual00, size_in_bytes_of(dual00),parai%io_source,parai%cp_grp)
       ! TESR
       CALL mp_bcast(ropt_mod%tesr,parai%io_source,parai%cp_grp)
       CALL mp_bcast(iteropt%iesr,parai%io_source,parai%cp_grp)
       ! ISOS
       CALL mp_bcast_byte(isos1, size_in_bytes_of(isos1),parai%io_source,parai%cp_grp)
       CALL mp_bcast_byte(isos2, size_in_bytes_of(isos2),parai%io_source,parai%cp_grp)
       CALL mp_bcast_byte(isos3, size_in_bytes_of(isos3),parai%io_source,parai%cp_grp)
       ! TSPHERE
       CALL mp_bcast(tsphere,parai%io_source,parai%cp_grp)
       ! CLSD (NLSD)
       CALL mp_bcast(clsd%nlsd,parai%io_source,parai%cp_grp)
       CALL mp_bcast(clsd%nlsx,parai%io_source,parai%cp_grp)
       ! SPIN2 (TLSE)
       CALL mp_bcast_byte(lspin1, size_in_bytes_of(lspin1),parai%io_source,parai%cp_grp)
       CALL mp_bcast_byte(lspin2, size_in_bytes_of(lspin2),parai%io_source,parai%cp_grp)
       CALL mp_bcast(lspin3%mgmab,parai%io_source,parai%cp_grp)
       CALL mp_bcast(lspin3%mgmba,parai%io_source,parai%cp_grp)
 
       CALL mp_bcast_byte(lspin4, size_in_bytes_of(lspin4),parai%io_source,parai%cp_grp)
       ! TFIXO
       CALL mp_bcast(tfixo,parai%io_source,parai%cp_grp)
       ! TKDP 
       CALL mp_bcast(tkdp,parai%io_source,parai%cp_grp)
       ! NK1,NK2,NK3
       CALL mp_bcast(kpts_com%nk1,parai%io_source,parai%cp_grp)
       CALL mp_bcast(kpts_com%nk2,parai%io_source,parai%cp_grp)
       CALL mp_bcast(kpts_com%nk3,parai%io_source,parai%cp_grp)
       CALL mp_bcast(cntl%tfield,parai%io_source,parai%cp_grp)
       CALL mp_bcast(aifield,3,parai%io_source,parai%cp_grp)
       ! HFX
       CALL mp_bcast_byte(hfxc2, size_in_bytes_of(hfxc2),parai%io_source,parai%cp_grp)
       ! DOUBLEGRID
       CALL mp_bcast_byte(tdgcomm, size_in_bytes_of(tdgcomm),parai%io_source,parai%cp_grp)
       CALL mp_bcast_byte(edgcomm, size_in_bytes_of(edgcomm),parai%io_source,parai%cp_grp)
       ! cntl%tshock
       CALL mp_bcast_byte(shock1, size_in_bytes_of(shock1),parai%io_source,parai%cp_grp)

    END SUBROUTINE broadcast_sysin
    ! ==--------------------------------------------------------------==
  END SUBROUTINE sysin
  ! ==================================================================

  ! ==================================================================
  ! CHECK THE INPUT FOR THE CELL DIMENSIONS AGAINST THE SYMMETRY
  ! MAKE SURE THAT THE PARAMETERS IN THE INPUT FILE ARE VALID.
  ! 01/2004 <axel.kohlmeyer@theochem.ruhr-uni-bochum.de>
  ! ==--------------------------------------------------------------==
  SUBROUTINE chkcsymm(boxeps)
    REAL(real_8)                             :: boxeps

    INTEGER                                  :: i
    LOGICAL                                  :: cellok

! VARIABLES
! ARE THE TESTS DISABLED?

    IF (boxeps < 0.0_real_8) RETURN

    ! WE ASSUME THE BEST AND STOP ONLY ON KNOWN ERRORS.
    cellok=.TRUE.
    IF (parm%ibrav == 0) THEN
       IF ((ABS(cell_com%celldm(4)) > boxeps)&
             .OR. (ABS(cell_com%celldm(5)) > BOXEPS)&
             .OR. (ABS(cell_com%celldm(6)) > BOXEPS)) THEN
               WRITE(0,'(A,A)') 'SUPERCELL IS NOT ORTHORHOMBIC AS',&
               ' REQUIRED FOR ISOLATED SYSTEMS'
               WRITE(0,'(A,A)') 'EXPECTED FORM:',&
               ' a b c 0.00.0 0.0'
          cellok=.FALSE.
       ENDIF
    ELSEIF (parm%ibrav == 1) THEN
       IF ((ABS(1.0_real_8-cell_com%celldm(2)) > boxeps)&
             .OR. (ABS(1.0_real_8-cell_com%celldm(3)) > BOXEPS)&
             .OR. (ABS(cell_com%celldm(4)) > BOXEPS)&
             .OR. (ABS(cell_com%celldm(5)) > BOXEPS)&
             .OR. (ABS(cell_com%celldm(6)) > BOXEPS)) THEN
               WRITE(0,*) 'SUPERCELL IS NOT SIMPLE CUBIC'
               WRITE(0,'(A,A)') 'EXPECTED FORM:',&
               ' a 1.01.0 0.00.0 0.0'
          cellok=.FALSE.
       ENDIF
    ELSEIF (parm%ibrav == 2) THEN
       IF ((ABS(1.0_real_8-cell_com%celldm(2)) > boxeps)&
             .OR. (ABS(1.0_real_8-cell_com%celldm(3)) > BOXEPS)&
             .OR. (ABS(cell_com%celldm(4)) > BOXEPS)&
             .OR. (ABS(cell_com%celldm(5)) > BOXEPS)&
             .OR. (ABS(cell_com%celldm(6)) > BOXEPS)) THEN
               WRITE(0,*) 'SUPERCELL IS NOT FACE CENTERED CUBIC'
               WRITE(0,'(A,A)') 'EXPECTED FORM:',&
               ' a 1.01.0 0.00.0 0.0'
          cellok=.FALSE.
       ENDIF
    ELSEIF (parm%ibrav == 3) THEN
       IF ((ABS(1.0_real_8-cell_com%celldm(2)) > boxeps)&
             .OR. (ABS(1.0_real_8-cell_com%celldm(3)) > BOXEPS)&
             .OR. (ABS(cell_com%celldm(4)) > BOXEPS)&
             .OR. (ABS(cell_com%celldm(5)) > BOXEPS)&
             .OR. (ABS(cell_com%celldm(6)) > BOXEPS)) THEN
               WRITE(0,*) 'SUPERCELL IS NOT BODY CENTERED CUBIC'
               WRITE(0,'(A,A)') 'EXPECTED FORM:',&
               ' a 1.01.0 0.00.0 0.0'
          cellok=.FALSE.
       ENDIF
    ELSEIF (parm%ibrav == 4) THEN
       IF ((ABS(1.0_real_8-cell_com%celldm(2)) > boxeps)&
             .OR. (ABS(cell_com%celldm(4)) > BOXEPS)&
             .OR. (ABS(cell_com%celldm(5)) > BOXEPS)&
             .OR. (ABS(cell_com%celldm(6)+0.5_real_8) > BOXEPS)) THEN
               WRITE(0,*) 'SUPERCELL IS NOT HEXAGONAL'
               WRITE(0,'(A,A)') 'EXPECTED FORM:',&
               ' a 1.0 cpa 0.00.0 -0.5'
          cellok=.FALSE.
       ENDIF
    ELSEIF (parm%ibrav == 5) THEN
       IF ((ABS(1.0_real_8-cell_com%celldm(2)) > boxeps)&
             .OR. (ABS(1.0_real_8-cell_com%celldm(3)) > BOXEPS)&
             .OR. (ABS(cell_com%celldm(5)-cell_com%celldm(4)) > BOXEPS)&
             .OR. (ABS(cell_com%celldm(6)-cell_com%celldm(4)) > BOXEPS)) THEN
               WRITE(0,*) 'SUPERCELL IS NOT RHOMBOHEDRAL'
               WRITE(0,'(A,A)') 'EXPECTED FORM:',&
               ' a 1.01.0 d d d'
          cellok=.FALSE.
       ENDIF
    ELSEIF (parm%ibrav == 6) THEN
       IF ((ABS(1.0_real_8-cell_com%celldm(2)) > boxeps)&
             .OR. (ABS(cell_com%celldm(4)) > BOXEPS)&
             .OR. (ABS(cell_com%celldm(5)) > BOXEPS)&
             .OR. (ABS(cell_com%celldm(6)) > BOXEPS)) THEN
               WRITE(0,*) 'SUPERCELL IS NOT TETRAGONAL'
               WRITE(0,'(A,A)') 'EXPECTED FORM:',&
               ' a 1.0 cpa 0.00.0 0.0'
          cellok=.FALSE.
       ENDIF
    ELSEIF (parm%ibrav == 7) THEN
       IF ((ABS(1.0_real_8-cell_com%celldm(2)) > boxeps)&
             .OR. (ABS(cell_com%celldm(4)) > BOXEPS)&
             .OR. (ABS(cell_com%celldm(5)) > BOXEPS)&
             .OR. (ABS(cell_com%celldm(6)) > BOXEPS)) THEN
               WRITE(0,*) 'SUPERCELL IS NOT BODY CENTERED TETRAGONAL'
               WRITE(0,'(A,A)') 'EXPECTED FORM:',&
               ' a 1.0 c 0.00.0 0.0'
          cellok=.FALSE.
       ENDIF
    ELSEIF (parm%ibrav == 8) THEN
       IF ((ABS(cell_com%celldm(4)) > boxeps)&
             .OR. (ABS(cell_com%celldm(5)) > BOXEPS)&
             .OR. (ABS(cell_com%celldm(6)) > BOXEPS)) THEN
               WRITE(0,*) 'SUPERCELL IS NOT ORTHORHOMBIC'
               WRITE(0,'(A,A)') 'EXPECTED FORM:',&
               ' a b c 0.00.0 0.0'
          cellok=.FALSE.
       ENDIF
    ELSEIF (parm%ibrav == 12) THEN
       IF ((ABS(cell_com%celldm(5)) > boxeps)&
             .OR. (ABS(cell_com%celldm(6)) > BOXEPS)) THEN
               WRITE(0,*) 'SUPERCELL IS NOT MONOCLINIC'
               WRITE(0,'(A,A)') 'EXPECTED FORM:',&
               ' a b c d 0.00.0'
          cellok=.FALSE.
       ENDIF
    ELSEIF (parm%ibrav == 14) THEN
       ! No sanity check possible with triclinic system
       cellok=.TRUE.
    ELSEIF (parm%ibrav == -2) THEN
       ! The cell was read from vectors -> implicit ibrav=14 (will be set in setsc.f)
       ! no sanity check possible with triclinic system
       cellok=.TRUE.
    ELSE
            WRITE(0,'(A,I2)') ' UNKNOWN SUPERCELL SYMMETRY IBRAV=', parm%ibrav
            WRITE(0,'(A,A)')  ' IF THIS IS A VALID SYMMETRY, PLEASE',&
            'ADD A SUITABLE TEST TO CHKCSYMM'
    ENDIF
    IF (.NOT.cellok) THEN
            WRITE(0,'(A,1PE10.3)') ' CHECK PRECISION: ', boxeps
            WRITE(0,'(A,6F10.5)' ) ' CELL DIMENSIONS: ', (cell_com%celldm(i),i=1,6)
       CALL stopgm('CHKCSYMM','SUPERCELL DOES NOT MATCH SYMMETRY',& 
            __LINE__,__FILE__)
    ENDIF

    ! ==--------------------------------------------------------------==
  END SUBROUTINE chkcsymm
  ! ==================================================================

END MODULE sysin_utils
