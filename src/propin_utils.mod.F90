MODULE propin_utils
  USE cnst,                            ONLY: fbohr
  USE condu,                           ONLY: condpa
  USE cores,                           ONLY: coresi,&
                                             coresl,&
                                             coresr
  USE error_handling,                  ONLY: stopgm
  USE g_loc,                           ONLY: indstate,&
                                             ist_list,&
                                             lostate
  USE inscan_utils,                    ONLY: inscan
  USE ldosmod,                         ONLY: cldos
  USE lodp,                            ONLY: &
       exd, extd, focc, nmaxld, nsdip, numbld, trmom, xmaxld, xminld, ymaxld, &
       yminld, zmaxld, zminld
  USE mp_interface,                    ONLY: mp_bcast
  USE parac,                           ONLY: parai,&
                                             paral
  USE pola,                            ONLY: ipolarise,&
                                             tpolarb,&
                                             tzeff
  USE prop,                            ONLY: icubeorb,&
                                             prop1,&
                                             prop2,&
                                             prop3,&
                                             prop4,&
                                             prop5,&
                                             prop7,&
                                             rcut
  USE readsr_utils,                    ONLY: input_string_len,&
                                             index_of_delimiter,&
                                             keyword_contains,&
                                             readsi,&
                                             readsr,&
                                             xstring
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             maxsys
  USE wann,                            ONLY: hmat_spread,&
                                             minspr,&
                                             real_8
  USE zeroing_utils,                   ONLY: zeroing
#include "sizeof.h"

  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER :: output_unit = 6

  PUBLIC :: propin

CONTAINS

  ! ==================================================================
  SUBROUTINE propin(need_property)
    ! ==--------------------------------------------------------------==
    ! ==  THIS ROUTINE READS THE SECTION &PROP &END ON UNIT IUNIT     ==
    ! ==  Input:                                                      ==
    ! ==  need_property  =.TRUE.  section reads at the beginning of   == 
    ! ==                          CPMD, read only                     ==
    ! ==                  .FALSE. initialize too                      ==
    ! ==  Output:         .TRUE. a property has to be calculated      ==
    ! ==--------------------------------------------------------------==
    ! ==  THE INPUT HAS TO BE IN THE FOLLOWING FORMAT                 ==
    ! ==     &PROP                                                    ==
    ! ==        Options                                               ==
    ! ==     &END                                                     ==
    ! ==--------------------------------------------------------------==
    ! ==  LIST OF OPTIONS                                             ==
    ! ==                                                              ==
    ! ==   PROJECT WAVEFUNCTIONS                                      ==
    ! ==   PROJECT HARMONICS [RADIAL] [SPREAD]                        ==
    ! ==      NUMBER OF HARMONICS, CENTER                             == 
    ! ==      [RMAX, NMAX]                                            ==
    ! ==      [SPREAD_MIN, SPREAD_MAX]                                ==
    ! ==--------------------------------------------------------------==
    ! ==    POPULATION ANALYSIS [MULLIKEN,DAVIDSON][n-CENTER]         ==
    ! ==     no1 ... nox                                              ==
    ! ==    n-CENTER CUTOFF                                           ==
    ! ==       cutoff                                                 ==
    ! ==    CHARGES [PARA, HYPER]                                     ==
    ! ==    OPTIMIZE SLATER EXPONENTS                                 ==
    ! ==    LOCALIZE                                                  ==
    ! ==    NOPRINT {ORBITALS}                                        ==
    ! ==    DIPOLE MOMENT {BERRY,RS}                                  ==
    ! IF   ==    LOCAL DIPOLE                                              ==
    ! IF   ==      numbld                                                  ==
    ! IF   ==      xminld(1) yminld(1) zminld(1)                           ==
    ! IF   ==      xmaxld(1) ymaxld(1) zmaxld(1)                           ==
    ! IF   ==      ...                                                     ==
    ! IF   ==      xminld(numbld) yminld(numbld) zminld(numbld)            ==
    ! IF   ==      xmaxld(numbld) ymaxld(numbld) zmaxld(numbld)            ==
    ! ==    EXCITED DIPOLE                                            ==
    ! ==      nexdip numorb                                           ==
    ! ==      n1 n2 n3 ... nx                                         ==
    ! ==      ...                                                     ==
    ! ==      n1 n2 n3 ... nx                                         ==
    ! ==    TRANSTION MOMENT                                          ==
    ! ==      ntrmom                                                  ==
    ! ==      ns1 ns2                                                 ==
    ! ==      ...                                                     ==
    ! ==    CONDUCTIVITY                                              ==
    ! ==      iconduct [STEP=condstep]                                ==
    ! ==    POLARISABILITY                                            ==
    ! ==      ipolarise                                               ==
    ! ==    CORE SPECTRA                                              ==
    ! ==      core_atom core_level                                    ==
    ! ==      core_nqsto core_lshell core_stoexp                      ==
    ! ==    LDOS                                                      ==
    ! ==    CUBECENTER                                                ==
    ! ==      cubecin(1) cubecin(2),cubecin(3)                        ==
    ! ==    CUBEFILE {POTENTIAL|ORBITALS|DENSITY} [HALFMESH]          ==
    ! ==      numorb               [for ORBITALS]                     ==
    ! ==      n1 n2 n3 n4                                             ==
    ! ==    GLOCALIZE                                                 ==
    ! ==    EPR                                                       ==
    ! ==    EFG                                                       ==
    ! ==    AVERAGED POTENTIAL                                        ==       
    ! ==                                                              ==
    ! ==--------------------------------------------------------------==
    LOGICAL                                  :: need_property

    CHARACTER(*), PARAMETER                  :: procedureN = 'propin'
    INTEGER, PARAMETER                       :: max_unknown_lines = 30 

    CHARACTER(len=input_string_len)          :: line, error_message, previous_line, &
                                                unknown(max_unknown_lines)
    INTEGER                                  :: i, ii, j, first, last, ierr, iunit, &
                                                nbr_unknown_lines
    LOGICAL                                  :: erread, is_at_start, &
                                                go_on_reading, something_went_wrong

    !
    is_at_start = need_property
    !
    IF (paral%io_parent) THEN
       iunit = cnti%insys
       !
       ! Variables for reading
       !
       nbr_unknown_lines = 0
       line              = ' '
       previous_line     = ' '
       error_message     = ' '
       !
       ! Defaults
       !
       iunit = cnti%insys
       prop1%pwfn = .FALSE.
       prop1%wfr1 = .FALSE.
       prop1%wfr2 = .FALSE.
       prop1%wfr3 = .FALSE.
       prop1%cor1 = .FALSE.
       prop1%cor2 = .FALSE.
       prop1%abas = .FALSE.
       prop1%ceig = .FALSE.
       prop1%locl = .FALSE.
       prop1%ldip = .FALSE.
       prop1%locd = .FALSE.
       prop1%lext = .FALSE.
       prop1%ltrm = .FALSE.
       prop1%lrsdip = .FALSE.
       prop1%prto = .TRUE.
       prop1%opts = .FALSE.
       prop1%mpan = .FALSE.
       prop1%dpan = .FALSE.
       need_property = .FALSE.
       prop1%espc = .FALSE.
       prop1%eldisp = .FALSE.
       prop1%glocl  = .FALSE.
       prop1%dberry  = .FALSE.
       lostate%state_all   = .FALSE.
       lostate%state_range = .FALSE.
       lostate%state_list  = .FALSE.
       coresl%tcores = .FALSE. ! cmb-ike
       prop1%tavgp = .FALSE.
       prop1%pylm = .FALSE.
       !
       ! Options which can be calculated during simulations
       ! (whatever that means?)
       !
       condpa%nfconduct = 0
       condpa%condstep = 0.5_real_8
       condpa%nconduct = 40
       IF (is_at_start) THEN
          condpa%tconduct=.FALSE.
          tpolarb =.FALSE.
          condpa%iconduct=1
          ipolarise=1
       ENDIF
       !
       prop4%tcubefile_pot  = .FALSE.
       prop4%tcubefile_orb  = .FALSE.
       prop4%tcubefile_dens = .FALSE.
       prop4%thalfmesh      = .FALSE.
       prop4%tcubecenter    = .FALSE.
       prop5%teprefg        = .FALSE.
       prop5%tefg           = .FALSE.
       !
       prop2%ncen=2
       extd%nexdip=0
       extd%ntrmom=0
       prop3%cut3o=0.05_real_8
       prop3%cut4o=0.01_real_8
       prop2%numorb=0
       numbld=0
       CALL zeroing(prop4%cubecin)!,3)
       !$omp parallel do private(I)
       DO i=1,nmaxld
          xminld(i)=0
          yminld(i)=0
          zminld(i)=0
          xmaxld(i)=0
          ymaxld(i)=0
          zmaxld(i)=0
       ENDDO
       !
       ! Is the &PROP(ERTIES) section somewhere?
       !
       ierr = inscan(iunit,'&PROP')
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
             ELSEIF ( keyword_contains(line,'PROJECT',and='WAVEFUNCTION') .OR. &
                      keyword_contains(line,'PROJECT',and='WAVEFUNCTIONS') ) THEN
                ! Project wavefunctions
                prop1%pwfn=.TRUE.
             ELSEIF ( keyword_contains(line,'PROJECT',and='HARMONICS') ) THEN
                ! Project wavefunction on spherical harmonics
                prop1%pylm=.TRUE.
                prop7%radial=.FALSE.
                prop7%spread_ham=.FALSE.
                IF ( keyword_contains(line,'RADIAL') ) prop7%radial = .TRUE.
                IF ( keyword_contains(line,'SPREAD') ) prop7%spread_ham = .TRUE.
                !
                READ(iunit,*,iostat=ierr) prop7%numylm, prop7%centylm(1:3)
                IF (ierr /= 0) THEN
                   error_message = 'ERROR WHILE READING LINE'
                   something_went_wrong = .true.
                   go_on_reading = .false.
                ENDIF
                IF (prop7%radial) THEN
                   READ(iunit,*,iostat=ierr) prop7%rylmax, prop7%nylmax
                   IF (ierr /= 0) THEN
                      error_message = 'ERROR WHILE READING LINE'
                      something_went_wrong = .true.
                      go_on_reading = .false.
                   ENDIF
                ENDIF
                IF (prop7%spread_ham)THEN
                   READ(iunit,*,iostat=ierr) prop7%spr_min, prop7%spr_max
                   IF (ierr /= 0) THEN
                      error_message = 'ERROR WHILE READING LINE'
                      something_went_wrong = .true.
                      go_on_reading = .false.
                   ENDIF
                ENDIF
             ELSEIF ( keyword_contains(line,'CHARGES') ) THEN
                cnti%lfit=0
                IF ( keyword_contains(line,'PARA') ) THEN
                   cnti%lfit=1
                ELSEIF ( keyword_contains(line,'HYPER') ) THEN
                   cnti%lfit=2
                ENDIF
                ! ..calculate atomic charges
                prop1%espc=.TRUE.
             ELSEIF ( keyword_contains(line,'AVERAGED',and='POTENTIAL') ) THEN
                ! ..calculate averaged potential
                READ(iunit,*,iostat=ierr) rcut
                IF (ierr /= 0) THEN
                   error_message = 'ERROR WHILE READING LINE'
                   something_went_wrong = .true.
                   go_on_reading = .false.
                ENDIF
                prop1%tavgp=.TRUE.
             ELSEIF ( keyword_contains(line,'DIPOLE') ) THEN
                ! Dipole moments
                IF ( keyword_contains(line,'MOMENT') ) THEN
                   prop1%ldip=.TRUE.
                   IF ( keyword_contains(line,'RS') ) prop1%lrsdip=.TRUE.
                   IF ( keyword_contains(line,'BERRY') ) THEN
                      prop1%ldip  =.FALSE.
                      prop1%dberry=.TRUE.
                   ENDIF
                ELSEIF ( keyword_contains(line,'BERRY') ) THEN
                   prop1%dberry=.TRUE.
                ELSEIF ( keyword_contains(line,'LOCAL') ) THEN
                   prop1%locd=.TRUE.
                   READ(iunit,*,iostat=ierr) numbld
                   IF (ierr /= 0) THEN
                      error_message = 'ERROR WHILE READING LINE'
                      something_went_wrong = .true.
                      go_on_reading = .false.
                   ENDIF
                   IF (numbld > nmaxld) CALL stopgm(procedureN,'Number of local dipoles exceeds maximum', &
                                                    __LINE__, __FILE__)
                   DO i=1,numbld
                      READ(iunit,*,iostat=ierr) xminld(i), yminld(i), zminld(i), xmaxld(i), ymaxld(i), zmaxld(i)
                   ENDDO
                   !
                   ! NB: This can be buggy if &PROP is read before &SYSTEM, make sure this is always the case!
                   !
                   IF (.NOT.cntl%bohr) THEN
                      xminld(:) = xminld(:)*fbohr
                      yminld(:) = yminld(:)*fbohr
                      zminld(:) = zminld(:)*fbohr
                      xmaxld(:) = xmaxld(:)*fbohr
                      ymaxld(:) = ymaxld(:)*fbohr
                      zmaxld(:) = zmaxld(:)*fbohr
                   ENDIF
                ENDIF
                IF ( keyword_contains(line,'EXCITED') ) THEN
                   prop1%lext=.TRUE.
                   READ(iunit,*,iostat=ierr) extd%nexdip,prop2%numorb
                   IF (ierr /= 0) THEN
                      error_message = 'ERROR WHILE READING LINE'
                      something_went_wrong = .true.
                      go_on_reading = .false.
                   ENDIF
                   ALLOCATE(focc(extd%nexdip*prop2%numorb),STAT=ierr)
                   IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                        __LINE__,__FILE__)
                   ALLOCATE(exd(3,nmaxld,extd%nexdip),STAT=ierr)
                   IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                        __LINE__,__FILE__)
                   DO i=1,extd%nexdip
                      ii=(i-1)*prop2%numorb
                      READ(iunit,*,iostat=ierr) (focc(j+ii),j=1,prop2%numorb)
                      IF (ierr /= 0) THEN
                         error_message = 'ERROR WHILE READING LINE'
                         something_went_wrong = .true.
                         go_on_reading = .false.
                      ENDIF
                   ENDDO
                ENDIF
             ELSEIF ( keyword_contains(line,'TRANSITION',and='MOMENT') ) THEN
                prop1%ltrm=.TRUE.
                READ(iunit,*,iostat=ierr) extd%ntrmom
                IF (ierr /= 0) THEN
                   error_message = 'ERROR WHILE READING LINE'
                   something_went_wrong = .true.
                   go_on_reading = .false.
                ENDIF
                ALLOCATE(nsdip(2,extd%ntrmom),STAT=ierr)
                IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                     __LINE__,__FILE__)
                ALLOCATE(trmom(3,nmaxld,extd%ntrmom),STAT=ierr)
                IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                     __LINE__,__FILE__)
                DO i=1,extd%ntrmom
                   READ(iunit,*,iostat=ierr) (nsdip(j,i),j=1,2)
                   IF (ierr /= 0) THEN
                      error_message = 'ERROR WHILE READING LINE'
                      something_went_wrong = .true.
                      go_on_reading = .false.
                   ENDIF
                ENDDO
             ELSEIF ( keyword_contains(line,'POPULATION',and='ANALYSIS') ) THEN
                ! Population analysis 
                prop1%pwfn=.TRUE.
                IF ( keyword_contains(line,'MULLIKEN') ) prop1%mpan=.TRUE.
                IF ( keyword_contains(line,'DAVIDSON') ) prop1%dpan=.TRUE.
                IF ( keyword_contains(line,'2-CENTER') ) prop2%ncen=2
                IF ( keyword_contains(line,'3-CENTER') ) prop2%ncen=3
                IF ( keyword_contains(line,'4-CENTER') ) prop2%ncen=4
                IF (prop1%dpan) READ(iunit,*,iostat=ierr) (prop2%maos(ii),ii=1,maxsys%nsx)
                IF (ierr /= 0) THEN
                   error_message = 'ERROR WHILE READING LINE'
                   something_went_wrong = .true.
                   go_on_reading = .false.
                ENDIF
             ELSEIF ( keyword_contains(line,'3-CENTER',and='CUTOFF') ) THEN
                READ(iunit,*,iostat=ierr) prop3%cut3o
                IF (ierr /= 0) THEN
                   error_message = 'ERROR WHILE READING LINE'
                   something_went_wrong = .true.
                   go_on_reading = .false.
                ENDIF
             ELSEIF ( keyword_contains(line,'4-CENTER',and='CUTOFF') ) THEN
                READ(iunit,*,iostat=ierr) prop3%cut4o
                IF (ierr /= 0) THEN
                   error_message = 'ERROR WHILE READING LINE'
                   something_went_wrong = .true.
                   go_on_reading = .false.
                ENDIF
             ELSEIF ( keyword_contains(line,'OPTIMIZE',and='SLATER') ) THEN
                ! Optimizs Slater exponents
                prop1%opts=.TRUE.
             ELSEIF ( keyword_contains(line,'LOCALIZE') ) THEN
                ! Localize Orbitals
                prop1%locl=.TRUE.
                hmat_spread=.FALSE.
                IF ( keyword_contains(line,'WAN_SPREAD') ) THEN
                   hmat_spread=.TRUE.
                   READ(iunit,*,iostat=ierr) minspr
                ENDIF
             ELSEIF ( keyword_contains(line,'NOPRINT') ) THEN
                ! Print flags 
                IF ( keyword_contains(line,'ORBITAL',alias='ORBITALS') ) prop1%prto=.FALSE.
             !
             ! Options which can be calculated during simulations
             !
             ELSEIF ( keyword_contains(line,'CONDUCTIVITY') ) THEN
                IF (is_at_start) THEN
                   condpa%tconduct = .TRUE.
                ELSE
                   condpa%tconduct = .NOT.condpa%tconduct! If already calculated -> False
                ENDIF
                ! Calculates conductivity at each iconduct steps.
                READ(iunit,'(A)',iostat=ierr) line
                IF (ierr /= 0) THEN
                   error_message = 'ERROR WHILE READING LINE'
                   something_went_wrong = .true.
                   go_on_reading = .false.
                ENDIF
                first = 1
                CALL readsi(line,first,last,condpa%iconduct,erread)
                first = index_of_delimiter(line,'STEP','=')
                IF (first > 0) THEN
                   CALL readsr(line,first,last,condpa%condstep,erread)
                   IF (condpa%condstep <= 0.0_real_8) condpa%condstep=0.5_real_8
                ENDIF
                first = index_of_delimiter(line,'N','=')
                IF (first > 0) THEN
                   CALL readsi(line,first,last,condpa%nconduct,erread)
                   IF (condpa%nconduct <= 0) condpa%nconduct=40
                ENDIF
                ! We need the eigenvalues (RESTART file)
                prop1%ceig=.TRUE.
                !
             ELSEIF ( keyword_contains(line,'POLARIZABILITY',alias='POLARIZABILITIES') .OR. &
                      keyword_contains(line,'POLARISABILITY',alias='POLARISABILITIES') ) THEN
                ! Calculates polarisability every ipolarise steps 
                READ(iunit,*,iostat=ierr) ipolarise
                IF (ierr /= 0) THEN
                   error_message = 'ERROR WHILE READING LINE'
                   something_went_wrong = .true.
                   go_on_reading = .false.
                ENDIF
                IF (is_at_start) THEN
                   tpolarb = .TRUE.
                   tzeff   = .TRUE.
                ELSE
                   tpolarb = .NOT.tpolarb! If already calculated -> False
                   tzeff   = .TRUE.
                ENDIF
                prop1%ceig=.TRUE.
             !
             ! mb-ike
             !
             ELSEIF ( keyword_contains(line,'CORE',and='SPECTRA') ) THEN
                ! Calculate core optical spectra
                READ(iunit,*,iostat=ierr) coresi%core_atom,coresr%core_level
                IF (ierr /= 0) THEN
                   error_message = 'ERROR WHILE READING LINE'
                   something_went_wrong = .true.
                   go_on_reading = .false.
                ENDIF
                READ(iunit,*,iostat=ierr) coresi%core_nqsto,coresi%core_lshell,coresr%core_stoexp
                IF (ierr /= 0) THEN
                   error_message = 'ERROR WHILE READING LINE'
                   something_went_wrong = .true.
                   go_on_reading = .false.
                ENDIF
                coresl%tcores=.TRUE.
                prop1%ceig=.TRUE.
                !
             ELSEIF ( keyword_contains(line,'LDOS') ) THEN
                ! Calculates LDOS every steps
                READ(iunit,*,iostat=ierr) cldos%nlayer,cldos%ldosn
                IF (ierr /= 0) THEN
                   error_message = 'ERROR WHILE READING LINE'
                   something_went_wrong = .true.
                   go_on_reading = .false.
                ENDIF
                cldos%tldos=.TRUE.
                prop1%ceig=.TRUE.
                !
             ELSEIF ( keyword_contains(line,'CUBECENTER') ) THEN
                ! Sets the center of the cubefile (default is geometric center).
                prop4%tcubecenter=.TRUE.
                READ(iunit,*,iostat=ierr) prop4%cubecin(1:3)
                IF (ierr /= 0) THEN
                   error_message = 'ERROR WHILE READING LINE'
                   something_went_wrong = .true.
                   go_on_reading = .false.
                ENDIF
                ! Angstrom?
                IF (.NOT.cntl%bohr) prop4%cubecin(:) = prop4%cubecin(:)*fbohr
                !
             ELSEIF ( keyword_contains(line,'CUBEFILE') ) THEN
                ! Plots potential, wfns or density as .cube file
                IF ( keyword_contains(line,'HALFMESH') )  prop4%thalfmesh = .TRUE.
                IF ( keyword_contains(line,'POTENTIAL') ) prop4%tcubefile_pot = .TRUE.
                IF ( keyword_contains(line,'DENSITY') )   prop4%tcubefile_dens = .TRUE.
                IF ( keyword_contains(line,'ORBITALS',alias='ORBITAL') .OR. &
                     keyword_contains(line,'WAVEFUNCTION') ) THEN
                   prop4%tcubefile_orb = .TRUE.
                   READ(iunit,*,iostat=ierr) prop4%num_cubeorb
                   IF (ierr /= 0) THEN
                      error_message = 'ERROR WHILE READING LINE'
                      something_went_wrong = .true.
                      go_on_reading = .false.
                   ENDIF
                   IF(ALLOCATED(icubeorb)) THEN
                      DEALLOCATE(icubeorb,STAT=ierr)
                      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
                           __LINE__,__FILE__)
                   ENDIF
                   ALLOCATE(icubeorb(prop4%num_cubeorb),STAT=ierr)
                   IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                        __LINE__,__FILE__)
                   CALL zeroing(icubeorb)
                   READ(iunit,*,iostat=ierr) (icubeorb(ii),ii=1,prop4%num_cubeorb)
                   IF (ierr /= 0) THEN
                      error_message = 'ERROR WHILE READING LINE'
                      something_went_wrong = .true.
                      go_on_reading = .false.
                   ENDIF
                ENDIF
                !
             ELSEIF ( keyword_contains(line,'EPR') ) THEN
                IF (cntl%tlsd) prop5%teprefg=.TRUE.
                !
             ELSEIF ( keyword_contains(line,'EFG') ) THEN
                prop5%tefg=.TRUE.
                !
             ELSEIF ( keyword_contains(line,'G_SPACE_LOC') ) THEN
                ! Localize Orbitals in Reciprocal space
                prop1%glocl=.TRUE.
                ! IF( keyword_contains(LINE,'CONJ_GRAD') ) THEN
                ! TCGRAD = .TRUE.
                ! ENDIF
                IF ( keyword_contains(line,'ALL') ) THEN
                   lostate%state_all   = .TRUE.
                   lostate%state_range = .FALSE.
                   lostate%state_list  = .FALSE.
                ELSEIF ( keyword_contains(line,'PART') ) THEN
                   lostate%state_all   = .FALSE.
                   lostate%state_range = .TRUE.
                   lostate%state_list  = .FALSE.
                   previous_line=line
                   READ(iunit,fmt='(A)',iostat=ierr) line
                   IF (ierr /= 0) THEN
                      error_message = 'ERROR WHILE READING LINE'
                      something_went_wrong = .true.
                      go_on_reading = .false.
                   ENDIF
                   first = 1
                   CALL readsi(line,first,last,indstate%ist_first,erread)
                   IF (erread) THEN
                      error_message = 'ERROR WHILE READING VALUE'
                      something_went_wrong = .true.
                      go_on_reading = .false.
                   ENDIF
                   first = last
                   CALL readsi(line,first,last,indstate%ist_last,erread)
                   IF (erread) THEN
                      error_message = 'ERROR WHILE READING VALUE'
                      something_went_wrong = .true.
                      go_on_reading = .false.
                   ENDIF
                ELSEIF ( keyword_contains(line,'LIST') )  THEN
                   lostate%state_all   = .FALSE.
                   lostate%state_range = .FALSE.
                   lostate%state_list  = .TRUE.
                   previous_line=line
                   READ(iunit,fmt='(A)',iostat=ierr) line
                   IF (ierr /= 0) THEN
                      error_message = 'ERROR WHILE READING LINE'
                      something_went_wrong = .true.
                      go_on_reading = .false.
                   ENDIF
                   CALL readsi(line,1,last,indstate%nst_list,erread)
                   IF (erread) THEN
                      error_message = 'ERROR WHILE READING VALUE'
                      something_went_wrong = .true.
                      go_on_reading = .false.
                   ENDIF
                   ALLOCATE(ist_list(indstate%nst_list),STAT=ierr)
                   IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                        __LINE__,__FILE__)
                   ! We can use more than one line
                   READ(iunit,*,iostat=ierr) (ist_list(i),i=1,indstate%nst_list)
                   IF (ierr /= 0) THEN
                      error_message = 'ERROR WHILE READING LINE'
                      something_went_wrong = .true.
                      go_on_reading = .false.
                   ENDIF
                ELSE
                   lostate%state_all   = .TRUE.
                   lostate%state_range = .FALSE.
                   lostate%state_list  = .FALSE.
                ENDIF
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
                      unknown(nbr_unknown_lines)=line
                   ENDIF
                ENDIF
                !
                ! Replace old GOTO 10/20 with a cycle (if property section is
                ! empty or cannot be read, we don't need need_property (was: TINFO)
                !
                CYCLE
             ENDIF
             !
             need_property = .true.
             !
          ENDDO
          !
          ! End of read-loop
          !
       ELSE
          !
          ! &PROP is not there... (this is legit)
          !
          something_went_wrong = .false.
       ENDIF
       !
       IF (something_went_wrong) THEN
           WRITE(output_unit,'(/,1X,64("!"))')
           WRITE(output_unit,'(1X, A, 1X, A)') 'ERROR:', 'PROBLEM WHILE READING &PROP SECTION:'
           WRITE(output_unit,'(8X, A)') TRIM(ADJUSTL(error_message))
           IF (line /= ' ' .or. previous_line /= ' ') THEN
              WRITE(output_unit,'(8X, A)') 'THE LAST TWO LINES READ WITHIN THE SECTION WERE:'
              WRITE(output_unit,'(/,1X, A)') TRIM(ADJUSTL(previous_line))
              WRITE(output_unit,'(1X, A)') TRIM(ADJUSTL(line))
           END IF
           WRITE(output_unit,'(1X,64("!"))')
           CALL stopgm(procedureN,'Error while reading &PROP section, cf. output file',&
                __LINE__,__FILE__)
       ENDIF
       ! 
       ! Stop if is_a_start=.true. (used to be TSTART) was replaced by an if-not
       !
       IF (.not. is_at_start) CALL propin_report()

    ENDIF

    CALL broadcast_propin()

    CONTAINS

    ! ==--------------------------------------------------------------==
    SUBROUTINE propin_report()

       WRITE(output_unit,*)
       WRITE(output_unit,'(A,/)') ' ACTIVE FLAGS FOR PROPERTIES RUN:'
       IF (prop1%ceig) WRITE(output_unit,'(A)')&
            ' READ EIGENVALUES FROM FILE RESTART'
       IF (prop1%pwfn) WRITE(output_unit,'(A)')&
            ' PROJECT WAVEFUNCTIONS ON ATOMIC PSEUDO-WAVEFUNCTIONS'
       IF (prop1%pylm) WRITE(output_unit,'(A)')&
            ' PROJECT WAVEFUNCTIONS ON SPHERICAL HARMONICS'
       IF (prop1%mpan) WRITE(output_unit,'(A)')&
            ' MULLIKEN POPULATION ANALYSIS AND MAYER BOND-ORDERS'
       IF (prop1%dpan) THEN
          WRITE(output_unit,'(A,/,A)') ' DAVIDSON POPULATION ANALYSIS',&
          '           C.Erhardt and R.Ahlrichs, TCA 68:231 (1985)'
          WRITE(output_unit,'(A,I3)')&
          '           HIGHEST N-CENTER TERMS INCLUDED : ',prop2%ncen
       ENDIF
       IF (prop1%espc) WRITE(output_unit,'(A,/,A,/,A)')&
            ' CALCULATE ATOMIC CHARGES FROM ',&
            '     REAL SPACE INTEGRATION AND ',&
            '     DERIVED FROM ELECTROSTATIC POTENTIAL '
       IF (prop1%espc.AND.cnti%lfit == 0)&
            WRITE(output_unit,'(A)') ' STANDARD METHOD '
       IF (prop1%espc.AND.cnti%lfit == 1)&
            WRITE(output_unit,'(A)') ' PARABOLIC RESTRAINT '
       IF (prop1%espc.AND.cnti%lfit == 2)&
            WRITE(output_unit,'(A)') ' HYPERBOLIC RESTRAINT '
       IF (prop1%locl)&
            WRITE(output_unit,'(A)') ' LOCALIZE ORBITALS'
       IF (prop1%glocl)&
            WRITE(output_unit,'(A)') ' LOCALIZE ORBITALS IN G-SPACE'
       IF (prop1%ldip) THEN
          IF (prop1%lrsdip) THEN
             WRITE(output_unit,'(A)') ' CALCULATE DIPOLE MOMENT IN REAL SPACE'
          ELSE
             WRITE(output_unit,'(A)') ' CALCULATE DIPOLE MOMENT '
          ENDIF
       ENDIF
       IF (prop1%dberry)&
            WRITE(output_unit,'(A,A)') ' CALCULATE DIPOLE MOMENT ',&
            'USING BERRY PHASE FORMULA'
       IF (prop1%locd)&
            WRITE(output_unit,'(A)') ' CALCULATE LOCAL DIPOLE MOMENTS '
       IF (prop1%lext)&
            WRITE(output_unit,'(A)') ' CALCULATE EXCITED STATE DIPOLE MOMENTS'
       IF (prop1%ltrm)&
            WRITE(output_unit,'(A)') ' CALCULATE TRANSITION MOMENTS'
       IF (condpa%tconduct)&
            WRITE(output_unit,'(A)') ' CALCULATE CONDUCTIVITY'
       IF (tpolarb)&
            WRITE(output_unit,'(A)') ' CALCULATE POLARISABILITY'
       IF (coresl%tcores) THEN! cmb-ike
          IF (paral%io_parent)&
               WRITE(output_unit,'(A,/,A,/,A,A,/,A,A)')&
               ' CALCULATE X-RAY SPECTRA AFTER: ',&
               '    M. Cavalleri et al. JCP 121, 10065 (2004)',&
               '    WARNING: IN PRINCIPLE THE METHOD IS VALID',&
               ' ONLY FOR K-EDGE OF ',&
               '    THE FIRST ROW ELEMENTS AND REQUIRES',&
               ' SPECIAL PSEUDOPOTENTIALS ! '
       ENDIF! cmb-ike
       IF (prop4%tcubecenter) THEN
          WRITE(output_unit,'(A)') ' SET CENTER OF CUBEFILES FROM INPUT'
       ENDIF
       IF (prop4%tcubefile_orb) THEN
          WRITE (6,'(A)') ' PLOTTING CUBEFILES (ORBITALS)'
          WRITE (6,'(15(1X,I0))') (icubeorb(ii),ii=1,prop4%num_cubeorb)
       ENDIF
       IF (prop4%tcubefile_dens)&
            WRITE (6,'(A)') ' PLOTTING CUBEFILES (DENSITY)'
       IF (prop4%tcubefile_pot)&
            WRITE (6,'(A)') ' PLOTTING CUBEFILES (POTENTIAL)'
       IF (prop1%tavgp)&
            WRITE(output_unit,'(A,/,A,1F5.2)')&
            ' CALCULATE AVERAGED ELECTROSTATIC POTENTIAL',&
            '     AROUND ATOMIC CORES RCUT= ',rcut
       WRITE(output_unit,*)
       IF (nbr_unknown_lines /= 0) THEN
          WRITE(output_unit,'(/,1X,64("="))')
          WRITE(output_unit,'(1X,A,14X,A,14X,A)') '= ', 'UNKNOWN KEYWORDS IN SECTION &PROP','='
          DO i=1,nbr_unknown_lines
             previous_line = unknown(i)
             CALL xstring(previous_line,first,last)
             WRITE(output_unit,'(A,A)') ' = ',previous_line(first:last)
          ENDDO
          WRITE(output_unit,'(1X,64("="),/)')
       ENDIF

    END SUBROUTINE propin_report
    ! ==--------------------------------------------------------------==
    SUBROUTINE broadcast_propin()

       CALL mp_bcast_byte(prop1,size_in_bytes_of(prop1),parai%io_source,parai%cp_grp)
       CALL mp_bcast_byte(extd, size_in_bytes_of(extd),parai%io_source,parai%cp_grp)
       CALL mp_bcast(extd%ntrmom,parai%io_source,parai%cp_grp)
       CALL mp_bcast(need_property,parai%io_source,parai%cp_grp)
       CALL mp_bcast_byte(prop2, size_in_bytes_of(prop2),parai%io_source,parai%cp_grp)
       IF (extd%nexdip /= 0) THEN
          IF (.NOT.paral%parent) THEN
             ALLOCATE(focc(extd%nexdip*prop2%numorb),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             ALLOCATE(exd(3,nmaxld,extd%nexdip),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
          ENDIF
          CALL mp_bcast(focc,SIZE(focc),parai%io_source,parai%cp_grp)
       ENDIF
       CALL mp_bcast_byte(prop2, size_in_bytes_of(prop2),parai%io_source,parai%cp_grp)
       CALL mp_bcast_byte(prop3,size_in_bytes_of(prop3),parai%io_source,parai%cp_grp)
       CALL mp_bcast(cldos%tldos,parai%io_source,parai%cp_grp)
       CALL mp_bcast_byte(cldos, size_in_bytes_of(cldos),parai%io_source,parai%cp_grp)
       CALL mp_bcast(cldos%ldosn,parai%io_source,parai%cp_grp)
       CALL mp_bcast(xminld,SIZE(xminld),parai%io_source,parai%cp_grp)
       CALL mp_bcast(numbld,parai%io_source,parai%cp_grp)
       ! 
       ! Maybe we should also broadcast the entries in condu.inc and pola.inc ?
       ! 
       CALL mp_bcast(prop1%ceig,parai%io_source,parai%cp_grp)
       CALL mp_bcast(condpa%nconduct,parai%io_source,parai%cp_grp)
       CALL mp_bcast(condpa%tconduct,parai%io_source,parai%cp_grp)
       CALL mp_bcast(condpa%condstep,parai%io_source,parai%cp_grp)
       CALL mp_bcast(tpolarb,parai%io_source,parai%cp_grp)
       CALL mp_bcast(tzeff,parai%io_source,parai%cp_grp)
       ! -ike
       CALL mp_bcast_byte(coresl, size_in_bytes_of(coresl),parai%io_source,parai%cp_grp)
       CALL mp_bcast_byte(coresi, size_in_bytes_of(coresi),parai%io_source,parai%cp_grp)
       CALL mp_bcast_byte(coresr, size_in_bytes_of(coresr),parai%io_source,parai%cp_grp)

       CALL mp_bcast_byte(prop4, size_in_bytes_of(prop4), &
            parai%io_source, parai%cp_grp)

       CALL mp_bcast_byte(prop5, size_in_bytes_of(prop5), &
            parai%io_source, parai%cp_grp)

       CALL mp_bcast(rcut,parai%io_source,parai%cp_grp)
       IF (prop4%tcubefile_orb) THEN
          IF (.NOT. paral%parent) THEN
             IF(ALLOCATED(icubeorb)) THEN
                DEALLOCATE(icubeorb,STAT=ierr)
                IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
                     __LINE__,__FILE__)
             ENDIF
             ALLOCATE(icubeorb(prop4%num_cubeorb),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             CALL zeroing(  icubeorb)!,prop4%num_cubeorb)
          ENDIF
          CALL mp_bcast(icubeorb,SIZE(icubeorb),parai%io_source,parai%cp_grp)
       ENDIF
       CALL mp_bcast(prop1%glocl,parai%io_source,parai%cp_grp)
       CALL mp_bcast_byte(lostate, size_in_bytes_of(lostate),parai%io_source,parai%cp_grp)
       CALL mp_bcast_byte(indstate, size_in_bytes_of(indstate),parai%io_source,parai%cp_grp)
       IF (lostate%state_list) THEN
          IF (.NOT.paral%parent) THEN
             ALLOCATE(ist_list(indstate%nst_list),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
          ENDIF
          CALL mp_bcast(ist_list,SIZE(ist_list),parai%io_source,parai%cp_grp)
       ENDIF
       IF (extd%ntrmom /= 0) THEN
          IF (.NOT.paral%parent) THEN
             ALLOCATE(nsdip(2,extd%ntrmom),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             ALLOCATE(trmom(3,nmaxld,extd%ntrmom),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
          ENDIF
          CALL mp_bcast(nsdip,SIZE(nsdip),parai%io_source,parai%cp_grp)
       ENDIF
       ! PROP7
       CALL mp_bcast_byte(prop7, size_in_bytes_of(prop7), &
            parai%io_source, parai%cp_grp)
       ! WANNIER SPEAD ORDERED HAMILTONIAN
       CALL mp_bcast(minspr,parai%io_source,parai%cp_grp)
       CALL mp_bcast(hmat_spread,parai%io_source,parai%cp_grp)

    END SUBROUTINE broadcast_propin
    ! ==--------------------------------------------------------------==
  END SUBROUTINE propin
  ! ==================================================================
END MODULE propin_utils
