MODULE pi_cntl_utils
  USE cotr,                            ONLY: cotc0,&
                                             lskcor
  USE error_handling,                  ONLY: stopgm
  USE glemod,                          ONLY: glepar
  USE inscan_utils,                    ONLY: inscan
  USE ions,                            ONLY: ions1
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: parai,&
                                             paral
  USE pimd,                            ONLY: maxnp,&
                                             pc_groups,&
                                             pimd1,&
                                             pimd2,&
                                             pimd3,&
                                             repfname
  USE readsr_utils,                    ONLY: input_string_len,&
                                             keyword_contains,&
                                             xstring
  USE store_types,                     ONLY: restart1
  USE system,                          ONLY: cntl,&
                                             cntr
#include "sizeof.h"

  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER :: output_unit = 6

  PUBLIC :: pi_cntl

CONTAINS

  ! ==================================================================
  SUBROUTINE pi_cntl
    ! ==--------------------------------------------------------------==
    ! ==  THIS ROUTINE READS THE SECTION &PIMD &END ON INPUT UNIT     ==
    ! ==--------------------------------------------------------------==
    ! ==  THE INPUT HAS TO BE IN THE FOLLOWING FORMAT                 ==
    ! ==     &PIMD                                                    ==
    ! ==        Options                                               ==
    ! ==     &END                                                     ==
    ! ==--------------------------------------------------------------==
    ! ==  LIST OF OPTIONS                                             ==
    ! ==                                                              ==
    ! ==    TROTTER DIMENSION                                         ==
    ! ==      np_total                                                ==
    ! ==    CENTROID DYNAMICS                                         ==
    ! ==    RING-POLYMER DYNAMICS                                     ==
    ! ==    CLASSICAL TEST                                            == 
    ! ==    STAGING                                                   ==
    ! ==      facstage                                                ==
    ! ==    NORMAL MODES                                              ==
    ! ==      facstage                                                ==
    ! ==    DEBROGLIE [CENTROID]                                      ==
    ! ==      tempb                                                   ==
    ! ==    FACMASS                                                   ==
    ! ==      wmass                                                   ==
    ! ==    INITIALIZATION                                            ==
    ! ==    READ REPLICAS                                             ==
    ! ==       filename                                               ==
    ! ==    GENERATE REPLICAS                                         ==
    ! ==    PROCESSOR GROUPS                                          ==
    ! ==       pc_groups                                              ==
    ! ==    OUTPUT [ALL,GROUPS,PARENT]                                ==
    ! ==    PRINT LEVEL                                               ==
    ! ==       levprt                                                 ==
    ! ==    GLE_LAMBDA                                                ==
    ! ==       gle_lambda                                             ==
    ! ==--------------------------------------------------------------==
    ! Variables
    CHARACTER(*), PARAMETER                  :: procedureN = 'pi_cntl'
    INTEGER, PARAMETER                       :: max_unknown_lines = 30

    CHARACTER(len=input_string_len)          :: line, error_message, previous_line, &
                                                unknown(max_unknown_lines)
    INTEGER                                  :: first, last, i, ierr, iunit, nbr_unknown_lines

    LOGICAL                                  :: go_on_reading, something_went_wrong

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
       ! Defaults
       !
       pimd1%rrep = .FALSE.
       pimd1%repread = .FALSE.
       pimd1%tread_cent  =  .TRUE.
       pimd1%tstage = .FALSE.
       pimd1%tpinm = .FALSE.
       pimd1%tinit = .FALSE.
       pimd1%tcentro = .FALSE.
       pimd1%tringp = .FALSE.
       pimd1%testpi = .FALSE.
       pimd3%np_total  =  1
       pimd2%tempb  =  500._real_8
       pimd2%wmass  =  1._real_8
       pimd2%facstage = 1._real_8
       pimd2%gle_lambda = 0.5_real_8
       pc_groups = 0
       repfname = ' '
       pimd3%levprt = 0
       pimd3%loutfn = 0
       !
       ! Is the &PIMD section there?
       !
       ierr = inscan(iunit,'&PIMD')
       !
       IF (ierr == 0) THEN
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
             ELSEIF ( keyword_contains(line,'TROTTER',and='DIMENSION') ) THEN
                ! Treat the nuclear coordinates in path integral representation
                ! with Trotter dimension NP_TOTAL
                READ(iunit,*,iostat=ierr) pimd3%np_total
                IF (ierr /= 0) THEN
                   error_message        = 'ERROR WHILE READING LINE'
                   something_went_wrong = .TRUE.
                   go_on_reading        = .FALSE.
                ENDIF
             ELSEIF ( keyword_contains(line,'CENTROID',and='DYNAMICS',alias='CMD') ) THEN
                ! do centroid path integral cntl%md and
                ! DM C do centroid path integral cntl%md (only together with normal modes) and
                ! do not thermostat centroid mode, i.e., first Trotter slice
                pimd1%tcentro=.TRUE.
             ELSEIF ( keyword_contains(line,'RING-POLYMER',and='DYNAMICS',alias='RPMD') ) THEN
                ! do ring-polymer path integral md 
                pimd1%tringp=.TRUE.
             ELSEIF ( keyword_contains(line,'CLASSICAL',and='TEST') ) THEN
                pimd1%testpi=.TRUE.
             ELSEIF ( keyword_contains(line,'STAGING') ) THEN
                ! use the staging representation for the path integrals
                pimd1%tstage=.TRUE.
                ! DM1 delete following two lines:
                ! the ficticious staging mass is (FACSTAGE*WMASS)-times heavier
                ! than the corresponding real physical mass of the same species
                ! 
                ! the fictitious non-staging masses are FACSTAGE time lighter than
                ! the staging mass (which has the physical mass as fictitious mass)
                ! DM2
                READ(iunit,*,iostat=ierr) pimd2%facstage
                IF (ierr /= 0) THEN
                   error_message        = 'ERROR WHILE READING LINE'
                   something_went_wrong = .TRUE.
                   go_on_reading        = .FALSE.
                ENDIF
             ELSEIF ( keyword_contains(line,'NORMAL',and='MODES') ) THEN
                ! use the normal mode representation for the path integrals
                pimd1%tpinm=.TRUE.
                ! DM1 delete following two lines:
                ! the ficticious normal mode mass is (FACSTAGE*WMASS)-times heavier
                ! than the corresponding real physical mass of the same species
                ! 
                ! the fictitious non-centroid masses are FACSTAGE time lighter than
                ! the centroid mass (which has the physical mass as fictitious mass)
                ! DM2
                READ(iunit,*,iostat=ierr) pimd2%facstage
                IF (ierr /= 0) THEN
                   error_message        = 'ERROR WHILE READING LINE'
                   something_went_wrong = .TRUE.
                   go_on_reading        = .FALSE.
                ENDIF
             ELSEIF ( keyword_contains(line,'DEBROGLIE') ) THEN
                ! DM1 C use TEMPB as the reference temperature (in Kelvin) to generate the
                ! DM C De Broglie displaced NP-1 other particles
                ! create free particle paths using a Gaussian Levy walk
                ! at temperature TEMPB (in Kelvin) according to the masses of the species
                ! starting from the first bead from the classical configuration
                ! that is obtained from input 
                ! DM2
                READ(iunit,*,iostat=ierr) pimd2%tempb
                IF (ierr /= 0) THEN
                   error_message        = 'ERROR WHILE READING LINE'
                   something_went_wrong = .TRUE.
                   go_on_reading        = .FALSE.
                ENDIF
                ! DM1
                pimd1%tread_cent = .FALSE.
                ! DM2
                ! DM1 C Classical coordinates are read in as centroid positions and
                ! DM C normal modes are sampled.
                ! create free particle paths using a Gaussian normal mode sampling
                ! at temperature TEMPB (in Kelvin) according to the masses of the species
                ! where the classical configuration yields the centroid positions
                ! DM2
                IF ( keyword_contains(line,'CENTROID',alias='CENTROIDS') ) pimd1%tread_cent = .TRUE.
             ELSEIF ( keyword_contains(line,'FACMASS') ) THEN
                ! DM1: delete the following 3 lines:
                ! MASS: ficticious ionic masses are M*P*(\beta *\hbar / WMASS)**2
                ! MASS: no specification of WMASS -> ficticious ionic masses are set to M
                ! MASS: (give WMASS in atomic units)
                ! DM2
                ! ficticious ionic masses are WMASS*M
                ! no specification of WMASS -> WMASS=1.0 (choice of Berne et al.)
                READ(iunit,*,iostat=ierr) pimd2%wmass
                IF (ierr /= 0) THEN
                   error_message        = 'ERROR WHILE READING LINE'
                   something_went_wrong = .TRUE.
                   go_on_reading        = .FALSE.
                ENDIF
             ELSEIF ( keyword_contains(line,'INITIALIZATION') ) THEN
                ! DM1 initialization run 
                ! initialization run: supply full quantum start configuration either by
                ! (i)  generation -> with RREP=.TRUE.    or by
                ! (ii) reading it -> with REPREAD=.TRUE.
                ! is default if no restart of coordinates
                ! DM2
                pimd1%tinit=.TRUE.
             ELSEIF ( keyword_contains(line,'GENERATE',and='REPLICAS') ) THEN
                ! DM     *    (INDEX(LINE,'REGENERATE').NE.0) THEN
                ! DM1 C generate replica coordinates
                ! generate replica coordinates as specified by 'DEBROGLIE' keyword
                ! DM2
                pimd1%rrep=.TRUE.
             ELSEIF ( keyword_contains(line,'READ',and='REPLICAS') ) THEN
                ! DM1 read replica coordinates from external file
                ! read replica coordinates from external file to be named 'repfname'  
                ! in the loop format "DO IP, DO IS, DO IA ... (K=1,3)" (see 'rread.F')
                ! with an empty line before EVERY set of replicas (also the first!)
                ! DM2
                pimd1%repread=.TRUE.
                READ(iunit,'(A80)',iostat=ierr) repfname
                IF (ierr /= 0) THEN
                   error_message        = 'ERROR WHILE READING LINE'
                   something_went_wrong = .TRUE.
                   go_on_reading        = .FALSE.
                ENDIF
             ELSEIF ( keyword_contains(line,'PROCESSOR',and='GROUPS') ) THEN
                ! number of processor goups 
                READ(iunit,*,iostat=ierr) pc_groups
                IF (ierr /= 0) THEN
                   error_message        = 'ERROR WHILE READING LINE'
                   something_went_wrong = .TRUE.
                   go_on_reading        = .FALSE.
                ENDIF
             ELSEIF ( keyword_contains(line,'PRINT',and='LEVEL') ) THEN
                ! DM1 printing level 
                ! printing level: minimal output for LEVPRT < 5 
                ! DM2
                READ(iunit,*,iostat=ierr) pimd3%levprt
                IF (ierr /= 0) THEN
                   error_message        = 'ERROR WHILE READING LINE'
                   something_went_wrong = .TRUE.
                   go_on_reading        = .FALSE.
                ENDIF
             ELSEIF ( keyword_contains(line,'OUTPUT') ) THEN
                ! output files for each processor, group or only grandparent
                IF ( keyword_contains(line,'PARENT') ) pimd3%loutfn=0
                IF ( keyword_contains(line,'GROUPS',alias='GROUP') ) pimd3%loutfn=1
                IF ( keyword_contains(line,'ALL') ) pimd3%loutfn=2
             ELSEIF ( keyword_contains(line,'GLE_LAMBDA') ) THEN
                READ(iunit,*,iostat=ierr) pimd2%gle_lambda
                IF (ierr /= 0) THEN
                   error_message        = 'ERROR WHILE READING LINE'
                   something_went_wrong = .TRUE.
                   go_on_reading        = .FALSE.
                ENDIF
             ELSE
                ! Dummy line
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
          ENDDO
       ENDIF
       !
       IF (something_went_wrong) THEN
          WRITE(output_unit,'(/,1X,64("!"))')
          WRITE(output_unit,'(1X, A, 1X, A)') 'ERROR:', 'PROBLEM WHILE READING &PIMD SECTION:'
          WRITE(output_unit,'(8X, A)') TRIM(ADJUSTL(error_message))
          IF (line /= ' ' .or. previous_line /= ' ') THEN
             WRITE(output_unit,'(8X, A)') 'THE LAST TWO LINES READ WITHIN THE SECTION WERE:'
             WRITE(output_unit,'(/,1X, A)') TRIM(ADJUSTL(previous_line))
             WRITE(output_unit,'(1X, A)') TRIM(ADJUSTL(line))
          END IF
          WRITE(output_unit,'(1X,64("!"))')
          CALL stopgm(procedureN,'Error while reading &PIMD section, cf. output file',&
               __LINE__,__FILE__)
       ENDIF
       !
       CALL check_options()
       CALL pi_cntl_report()

       IF (nbr_unknown_lines /= 0) THEN
          WRITE(output_unit,'(/,1X,64("="))')
          WRITE(output_unit,'(1X,A,14X,A,14X,A)') '= ','UNKNOWN KEYWORDS IN SECTION &PIMD','='
          DO i=1,nbr_unknown_lines
             previous_line = unknown(i)
             CALL xstring(previous_line,first,last)
             WRITE(output_unit,'(A,A)') ' = ',previous_line(first:last)
          ENDDO
          WRITE(output_unit,'(1X,64("="),/)')
       ENDIF

    END IF

    CALL broadcast_pi_cntl()


    CONTAINS

    ! ==--------------------------------------------------------------==
    SUBROUTINE check_options()

       INTEGER                                  :: ia, max_pc_groups, nfix

       !
       ! Test of options
       !
       ! DM1
       IF (pimd3%np_total.GT.maxnp) THEN
          WRITE(output_unit,'(A,I4,A,A,I4)')&
               ' SELECTED TROTTER NUMBER ',pimd3%np_total,' LARGER THAN',&
               ' CURRENT MAXIMUM VALUE MAXNP = ',maxnp
          CALL stopgm(procedureN,'Trotter number larger than maxnp.',& 
               __LINE__,__FILE__)
       ENDIF
     ! IF (cntl%tpres) CALL stopgm(procedureN,'Stress calculations not possible with path integrals.', &
     !                __LINE__, __FILE__)
       ! DM2
       IF (.NOT.pimd1%testpi) THEN
          IF (pimd3%np_total == 1) THEN
             IF (pimd1%tstage.OR.pimd1%tpinm) THEN
                IF (paral%io_parent)&
                     WRITE(output_unit,'(A)') ' STAGING OR NORMAL OPTION NEED NP_TOTAL > 1'
                CALL stopgm(procedureN,' WRONG OPTION ',& 
                     __LINE__,__FILE__)
             ENDIF
          ENDIF
       ENDIF
       IF (cntl%md) THEN
       !  pimd1%rrep=.FALSE.
       !  pimd1%repread=.FALSE.
          IF (cntr%tempw.LE.1.0e-5_real_8) CALL stopgm(procedureN,&
               'PIMD: You have to set temperature or use a thermostat',& 
               __LINE__,__FILE__)
          IF ((cntl%trane).AND.paral%io_parent) CALL stopgm(procedureN,&
               'RANDOMIZE WAVEFUNCTIONS not implemented with PIMD.',&
               __LINE__,__FILE__)
       !  IF (cntl%quenchb) CALL stopgm(procedureN,&
       !       'QUENCH BO  not implemented with PIMD.',& 
       !       __LINE__,__FILE__)
          ! check selected number of pc_groups not negative 
          IF (pc_groups < 0) CALL stopgm(procedureN,'Negative number for PI PC_GROUPS',&
               __LINE__,__FILE__)
          ! check selected number of pc_groups not exceed the upper limit
          max_pc_groups = MIN(pimd3%np_total,parai%cp_nproc)
          IF (pc_groups>max_pc_groups) pc_groups=max_pc_groups
          IF (cntl%tprcp) THEN
             IF (.NOT.(pimd1%tpinm.OR.pimd1%tstage)) CALL stopgm(procedureN,&
                'NPT PIMD needs staging or normal mode representation',&
                __LINE__,__FILE__)
          ENDIF
          IF (.NOT.restart1%restart) THEN
             pimd1%tinit=.TRUE.
             IF (.NOT.cntl%tmdbo) cntl%quenchb=.TRUE.
             IF (.NOT.(pimd1%rrep.OR.pimd1%repread)) pimd1%rrep=.TRUE.
          ENDIF
       ELSEIF (cntl%wfopt) THEN
          pc_groups = 1
          IF (.NOT.restart1%restart) THEN
             pimd1%tinit=.TRUE.
          ENDIF
          IF (pimd1%tinit) THEN
             IF (.NOT.(pimd1%rrep.OR.pimd1%repread)) pimd1%rrep=.TRUE.
          ELSE
             restart1%rco=.TRUE.
          ENDIF
       ELSE
          CALL stopgm(procedureN,'Path integrals: Only MD or wavefunction optimization are available.',& 
               __LINE__,__FILE__)
       ENDIF
       nfix=0
       DO ia=1,ions1%nat
          IF (lskcor(1,ia) == 0) nfix=nfix+1
          IF (lskcor(2,ia) == 0) nfix=nfix+1
          IF (lskcor(3,ia) == 0) nfix=nfix+1
       ENDDO
       nfix=nfix+cotc0%mcnstr
       IF (nfix > 0 .AND.(.NOT.(pimd1%tpinm.OR.pimd1%tstage))) CALL stopgm(procedureN,&
               'CONSTRAINTS only available with STAGING or NORMAL MODE representation',& 
               __LINE__,__FILE__)
       ! DM1
       IF (pimd1%tcentro) THEN
          IF ( .NOT.pimd1%tpinm ) WRITE(output_unit,'(/,A,A,/)')&
               ' WARNING: ONLY NORMAL MODE PROPAGATOR YIELDS',&
               ' CANONICAL DENSITY DISTRIBUTION: PLEASE VERIFY YOUR INPUT !'
          IF (pimd2%wmass /= 1.0_real_8) CALL stopgm(procedureN,&
                'WMASS must be unity for CENTROID path integral MD.',&
                 __LINE__,__FILE__)
       ELSE IF (pimd2%facstage /= 1.0_real_8) THEN
          WRITE(output_unit,'(/,A,/)')&
                ' WARNING: FACSTAGE SHOULD BE UNITY, PLEASE VERIFY YOUR INPUT'
       ENDIF
       ! DM2
       IF (pimd1%tringp) THEN
          IF (pimd1%tstage) CALL stopgm(procedureN,&
             'STAGING representation not yet supported for RING-POLYMER DYNAMICS',&
             __LINE__,__FILE__)
       ENDIF
       IF (pimd2%gle_lambda < 0.0_real_8) CALL stopgm(procedureN,&
          'Negative GLE_LAMBDA not acceptable',&
          __LINE__,__FILE__)

    END SUBROUTINE check_options
    ! ==--------------------------------------------------------------==
    SUBROUTINE pi_cntl_report()

       WRITE(output_unit,*)
       WRITE(output_unit,'(1X,19(">"),1x,A,1x,19("<"))') 'PATH-INTEGRAL PARAMETERS'
       ! DM1
       IF (pimd1%testpi) THEN
          pimd3%np_total=1
          pimd1%tread_cent=.FALSE.
          pimd1%tstage=.FALSE.
          pimd1%tpinm=.FALSE.
          pc_groups=1
          WRITE(output_unit,'(A)')&
               ' WARNING: THIS IS A TEST WITH NP=1 FOR COMPARISON WITH CLASSICAL RUNS'
          WRITE(output_unit,'(A)')&
               '          PRIMITIVE PROPAGATOR AND LEVY INITIAL GUESS USED'
       ENDIF
       IF (pimd1%tcentro) THEN
          WRITE(output_unit,'(A)')&
               ' ADIABATIC QUASICLASSICAL CENTROID DYNAMICS'
          IF (pimd1%tpinm) WRITE(output_unit,'(A)')&
               ' YIELDING THE CANONICAL BOLTZMANN DISTRIBUTION'
          IF (pimd1%tstage) WRITE(output_unit,'(A)')&
               ' YIELDING THE WIGNER DISTRIBUTION'
       ENDIF
       IF (pimd1%tringp) THEN
          WRITE(output_unit,'(A)')&
               ' QUASICLASSICAL RING-POLYMER DYNAMICS'
       ENDIF
       ! DM2
       WRITE(output_unit,'(A,7x,I4)') ' TROTTER DIMENSION:',pimd3%np_total
       IF (pimd1%tinit) WRITE(output_unit,'(A)') ' INITIALIZATION RUN'
       !
       ! DM1
       ! IF(RREP) WRITE(output_unit,'(A)') ' >>>> GENERATE REPLICA COORDINATES '
       IF (pimd1%rrep) THEN
          WRITE(output_unit,'(A)')&
               ' GENERATE QUANTUM FREE PARTICLE REPLICA COORDINATES '
          IF (pimd1%tread_cent) THEN
             WRITE(output_unit,'(A)')&
                  ' USING THE CLASSICAL POSITIONS AS CENTROIDS'
          ELSE
             WRITE(output_unit,'(A)')&
                  ' USING THE GAUSSIAN LEVY WALK IN REAL SPACE'
          ENDIF
          WRITE(output_unit,'(A,2x,F12.2,A)')&
               ' AT A TEMPERATURE OF:',pimd2%tempb,'K'
       ENDIF
       ! DM2
       CALL xstring(repfname,first,last)
       IF (pimd1%repread) WRITE(output_unit,'(A,A)')&
            ' READ REPLICA COORDINATES FROM FILE ',trim(adjustl(repfname(first:last)))
       ! DM1
       IF (pimd1%tstage) THEN
          WRITE(output_unit,'(A)') ' USE STAGING VARIABLES '
       ELSEIF (pimd1%tpinm) THEN
          WRITE(output_unit,'(A)') ' USE NORMAL MODE VARIABLES '
       ELSE
          WRITE(output_unit,'(A)') ' USE PRIMITIVE VARIABLES '
       ENDIF
       ! DM2
       IF ((pimd1%tpinm .OR. pimd1%tstage).AND.(.NOT.pimd1%tringp)) THEN
          ! DM        WRITE(output_unit,'(A,G12.6)') ' >>>> FACMASS : ',WMASS
          WRITE(output_unit,'(A,4x,G12.6)') ' SCALING FACTOR FACMASS:',pimd2%wmass
          ! DM        WRITE(output_unit,'(A,G12.6)') ' >>>> FACSTAGE: ',FACSTAGE
          WRITE(output_unit,'(A,3x,G12.6)') ' MASS DISPARITY FACSTAGE:',pimd2%facstage
       ENDIF
       IF (glepar%gle_mode>0) THEN
          WRITE(output_unit,'(A,2x,G12.6)') ' SCALING FACTOR GLE_LAMBDA:',pimd2%gle_lambda
       ENDIF
       WRITE(output_unit,'(A,13x,I4)') ' PRINT LEVEL:',pimd3%levprt

    END SUBROUTINE pi_cntl_report
    ! ==--------------------------------------------------------------==
    SUBROUTINE broadcast_pi_cntl

       ! PIMD 
       CALL mp_bcast_byte(pimd1, size_in_bytes_of(pimd1),parai%io_source,parai%cp_grp)
       CALL mp_bcast_byte(pimd2, size_in_bytes_of(pimd2),parai%io_source,parai%cp_grp)
       CALL mp_bcast_byte(pimd3, size_in_bytes_of(pimd3),parai%io_source,parai%cp_grp)

    END SUBROUTINE broadcast_pi_cntl
    ! ==--------------------------------------------------------------==
  END SUBROUTINE pi_cntl
  ! ==================================================================
END MODULE pi_cntl_utils
