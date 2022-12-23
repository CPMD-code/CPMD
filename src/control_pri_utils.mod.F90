MODULE control_pri_utils
  USE andr,                            ONLY: andr2,&
                                             andr3
  USE atwf,                            ONLY: dmovrmix,&
                                             tmovr
  USE broy,                            ONLY: broy1
  USE comvelmod,                       ONLY: comvl
  USE cotr,                            ONLY: sdpl
  USE envj,                            ONLY: tjlimit
  USE fileopenmod,                     ONLY: fo_info
  USE fint,                            ONLY: fint1,&
                                             fint4
  USE glemod,                          ONLY: gle_cpmd,&
                                             gle_cust,&
                                             gle_omega_def,&
                                             gle_opt,&
                                             gle_smart,&
                                             gle_white,&
                                             glepar
  USE isos,                            ONLY: isos1
  USE kpts,                            ONLY: tkpts
  USE lscal,                           ONLY: mode,&
                                             nsvib
  USE nort,                            ONLY: nort_com
  USE nose,                            ONLY: loct,&
                                             loctpin,&
                                             nosl,&
                                             tcafes
  USE parac,                           ONLY: paral
  USE qspl,                            ONLY: qspl1
  USE readsr_utils,                    ONLY: xstring
  USE shop,                            ONLY: s0_filn,&
                                             s1_filn
  USE store_types,                     ONLY: cprint,&
                                             iface1,&
                                             intfn,&
                                             restart1,&
                                             rout1,&
                                             store1,&
                                             trajsmall,&
                                             trajsmalln
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             cntr,&
                                             group,&
                                             restf
  USE vdwcmod,                         ONLY: vdwl
  USE wann,                            ONLY: sw_list,&
                                             wanni,&
                                             wannl,&
                                             wannr
  USE xinr,                            ONLY: gnx_inr,&
                                             inr_integer,&
                                             inr_logical,&
                                             real_8,&
                                             rmixsd,&
                                             tolx_inr

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: control_pri

CONTAINS

  ! ==================================================================
  SUBROUTINE control_pri(unknown,iunk,fformat)
    ! ==--------------------------------------------------------------==
    ! ==  Prints info on variables of control file                    ==
    ! ==--------------------------------------------------------------==

    ! INPUT
    CHARACTER(len=60)                        :: unknown(30)
    INTEGER                                  :: iunk
    CHARACTER(len=30)                        :: fformat

    INTEGER                                  :: i, ia, ie

! Variables
! INTEGER  IT !! unused variable
! ==--------------------------------------------------------------==
! ==  WRITE INFO TO OUTPUT                                        ==
! ==--------------------------------------------------------------==

    IF (.NOT. paral%io_parent) RETURN
    ! Debug options
    IF (cntl%tdebfor.OR.store1%tdebio) THEN
       WRITE(6,'(/,A)') ' DEBUG OPTIONS:'
       IF (cntl%tdebfor) THEN
          WRITE(6,'(4X,A)') 'DEBUG FORCES'
       ENDIF
       IF (store1%tdebio) THEN
          WRITE(6,'(4X,A)') 'DEBUG RESTART (READ AND WRITE)'
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (cntl%tpath) THEN
       IF (cntl%tpimd) THEN
          WRITE(6,'(/,A,/)') ' PATH INTEGRAL CALCULATION'
       ELSEIF (cntl%tpmin) THEN
          WRITE(6,'(/,A,/)') ' PATH MINIMISATION CALCULATION'
       ENDIF
    ENDIF
    IF (cntl%md) THEN
       IF (cntl%tprcp) THEN
          IF (cntl%tmdbo) THEN
             WRITE(6,'(/,A,A,/)') ' PARRINELLO-RAHMAN VARIABLE CELL',&
                  ' BORN-OPPENHEIMER MOLECULAR DYNAMICS '
          ELSEIF(cntl%tmdfile) THEN
             WRITE(6,'(/,A,A)') ' PARRINELLO-RAHMAN VARIABLE CELL',&
                  ' MOLECULAR DYNAMICS FROM FILES' 
             WRITE(6,'(T30,A,T61,I5)')   '       NSKIP   = ',cnti%nskip
             WRITE(6,'(T30,A,T61,I5,/)') '       NSAMPLE = ',cnti%nsample

          ELSE
             WRITE(6,'(/,A,A,/)') ' PARRINELLO-RAHMAN VARIABLE CELL',&
                  ' CAR-PARRINELLO MOLECULAR DYNAMICS '
          ENDIF
          IF (cntl%tpenn)  WRITE(6,'(A,/)')&
               ' GENERATE NPT-ENSEMBLE (PENN VERSION)'
          IF (cntl%tshock)  WRITE(6,'(A,/)')&
               ' SHOCK WAVE SIMULATION'
       ELSEIF (cntl%tmdbo) THEN
          WRITE(6,'(/,A,/)') ' BORN-OPPENHEIMER MOLECULAR DYNAMICS'
       ELSEIF (cntl%tmdfile) THEN
          WRITE(6,'(/,A,T56,A)')&
               ' MOLECULAR DYNAMICS FROM FILE',' TRAJSAVED'
          WRITE(6,'(T30,A,T61,I5)')   '       NSKIP   = ',cnti%nskip
          WRITE(6,'(T30,A,T61,I5,/)') '       NSAMPLE = ',cnti%nsample
       ELSEIF (cntl%tshop) THEN
          WRITE(6,'(/,A,/)') ' SURFACE HOPPING MOLECULAR DYNAMICS'
       ELSE
          WRITE(6,'(/,A,/)') ' CAR-PARRINELLO MOLECULAR DYNAMICS'
       ENDIF
    ELSEIF (cntl%geopt) THEN
       WRITE(6,'(/,A,/)') ' OPTIMIZATION OF IONIC POSITIONS'
    ELSEIF (cntl%wfopt) THEN
       WRITE(6,'(/,A,/)') ' SINGLE POINT DENSITY OPTIMIZATION'
       IF (cntl%tresponse)&
            WRITE(6,'(A,/)') ' WITH LINEAR RESPONSE CALCULATION '

    ELSEIF (cntl%tinter) THEN
       WRITE(6,'(/,A,/)') ' INTERFACE TO CLASSICAL MD PROGRAM'
    ELSEIF (cntl%tsampl) THEN
       WRITE(6,'(/,A,/)') ' REACTION PATH SAMPLING PROGRAM'
    ENDIF
    IF (cntl%ksener) THEN
       WRITE(6,'(/,A)')&
            ' EXACT DIAGONALIZATION OF KOHN-SHAM MATRIX '
       IF (tkpts%tknoswap) THEN
          WRITE(6,'(T28,A,/,T37,A)')&
               '[ONLY EIGENVALUES VALUES ARE RELEVANT]',&
               '[WAVEFUNCTIONS ARE NOT SAVED]'
       ENDIF
       WRITE(6,*)
    ENDIF
    IF (cntl%vibrat.AND.cntl%tsdin) THEN
       WRITE(6,'(/,A)')&
            ' PERFORM A VIBRATIONAL ANALYSIS FROM A GIVEN HESSIAN'
       IF (mode.NE.0) THEN
          WRITE(6,'(A,T63,I3,/)') ' PRINT HESSIAN EIGENMODE ', modE
       ELSE
          WRITE(6,'(/)')
       ENDIF
    ELSEIF (cntl%vibrat.AND..NOT.cntl%tsdan) THEN
       WRITE(6,'(/,A,/)')&
            ' PERFORM A VIBRATIONAL ANALYSIS BY FINITE DIFFERENCES'
    ELSEIF (cntl%vibrat.AND.cntl%tsdan) THEN
       WRITE(6,'(/,A,/)')&
            ' PERFORM A VIBRATIONAL ANALYSIS BY LINEAR RESPONSE'
    ENDIF
    IF (cntl%proper) THEN
       WRITE(6,'(/,A,/)') ' CALCULATE SOME PROPERTIES'
    ENDIF
    IF (cntl%thard) THEN
       WRITE(6,'(/,A)') ' CALCULATE ORBITAL HARDNESS MATRIX '
       IF (cntl%tsdan) THEN
          WRITE(6,'(A)') ' BY LINEAR RESPONSE ALGORITHM'
       ELSE
          WRITE(6,'(A)') ' BY FINITE DIFFERENCES'
       ENDIF
    ENDIF
    IF (cntl%tspec) THEN
       WRITE(6,'(/,A)') ' CALCULATE ELECTRONIC SPECTRA '
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (cntl%tfint) THEN
       WRITE(6,'(/,A)', advance='no') ' FREE ENERGY FUNCTIONAL'
       IF (fint1%ttrot) THEN
          IF (fint1%tbogo) THEN
             WRITE(6,'(/,T14,A)', advance='no')&
                  'WITH TROTTER FACTORISATION AND BOGOLIUBOV CORRECTION'
          ELSE
             WRITE(6,'(A)', advance='no') ' WITH TROTTER FACTORISATION'
          ENDIF
       ENDIF
       WRITE(6,'(/,1X,A,A,/)')&
            'A.ALAVI, J.KOHANOFF, M.PARRINELLO, D.FRENKEL',&
            '  PRL 732599 (1994)'
    ENDIF
    IF (cntl%tlsd) THEN
       WRITE(6,'(/,A)') ' LOCAL SPIN DENSITY APPROXIMATION'
    ENDIF
    IF (cntl%tddft) THEN
       WRITE(6,'(/,A)') ' LINEAR RESPONSE TO TIME-DEPENDENT DFT'
    ENDIF
    IF (cntl%nonort) THEN
       WRITE(6,'(/,A,A)') ' USE NONORTHOGONAL ORBITALS WITH',&
            ' NORM CONSTRAINTS'
       WRITE(6,'(A,1PE14.8,/)') ' OVERLAP MATRIX LIMIT: ',nort_com%slimit
    ENDIF
    IF (cntl%tharm) THEN
       WRITE(6,'(A)') ' HARMONIC REFERENCE SYSTEM INTEGRATION'
    ENDIF
    IF (trajsmall) THEN
       WRITE(6,*) ' WRITING ONLY A small TRAJECTORY'
       WRITE(6,*) ' meaning only ',trajsmalln,' atoms per frame'
    ENDIF
    IF (cnti%iprng.GT.0) THEN
       WRITE(6,*) 'USING SEED ', cnti%iprng,&
            'TO INIT. PSEUDO RANDOM NUMBER GEN.'
    ENDIF
    IF (glepar%gle_mode.GT.0) THEN
       WRITE(6,*) 'GENERALIZED LANGEVIN THERMOSTAT LOADED.'
       IF (glepar%gle_mode.NE.gle_cust) THEN
          IF (glepar%gle_omega .EQ. 0.0_real_8) THEN
             WRITE(6,*) '  USING DEFAULT OMEGA0 ', gle_omega_def,&
                  ' CM^-1'
             glepar%gle_omega=gle_omega_def
          ELSE
             WRITE(6,*) '  OMEGA0 SET TO', glepar%gle_omega, ' CM^-1'
          ENDIF
       ENDIF

       IF (glepar%gle_mode .EQ. gle_white) THEN
          WRITE(6,*) '  USING WHITE NOISE'
       ELSEIF (glepar%gle_mode .EQ. gle_opt) THEN
          WRITE(6,*) '  USING OPTIMAL SAMPLING, omega_min=',&
               0.0001_real_8*glepar%gle_omega, "CM^-1; omega_max=",&
               1.0_real_8*glepar%gle_omega, "CM^-1"
       ELSEIF (glepar%gle_mode .EQ. gle_cpmd)  THEN
          WRITE(6,*) '  USING STOCHASTIC CPMD VERSION omega_cut=',&
               glepar%gle_omega, "CM^-1"
       ELSEIF (glepar%gle_mode .EQ. gle_smart) THEN
          WRITE(6,*) '  USING SMART SAMPLING, tau_opt=',&
               1.e3_real_8/(glepar%gle_omega*0.029979246_real_8), "ps; omega_max=",&
               1.0_real_8*glepar%gle_omega, "CM^-1"             
       ELSEIF (glepar%gle_mode .EQ. gle_cust) THEN
          WRITE(6,*) '  USING CUSTOM MATRIX WITH ns=', glepar%gle_ns
          WRITE(6,*)&
               '  PARAMETERS WILL BE READ FROM ./GLE-A, [./GLE-C]'
       ENDIF
    ENDIF
    IF (cntl%tmass) THEN
       WRITE(6,'(A)') ' SCALED ELECTRON MASSES'
    ENDIF
    WRITE(fformat,'(A,I2,A)') '(A,T',MAX(29,65-(fo_info%iepath-fo_info%iapath)),',A)'
    WRITE(6,fformat) ' PATH TO THE RESTART FILES: ',&
         fo_info%fpath(fo_info%iapath:fo_info%iepath)
    IF (restart1%restart) THEN
       IF (restart1%rwf)   WRITE(6,'(A)') ' RESTART WITH OLD ORBITALS'
       IF (restart1%rrho)  WRITE(6,'(A)') ' RESTART WITH OLD DENSITY'
       IF (restart1%rocc)  WRITE(6,'(A)') ' RESTART WITH OLD OCCUPATION NUMBERS'
       IF (restart1%rkpt)  WRITE(6,'(A)') ' RESTART WITH OLD KPOINTS'
       IF (restart1%rco)   WRITE(6,'(A)') ' RESTART WITH OLD ION POSITIONS'
       IF (restart1%rvel)  WRITE(6,'(A)') ' RESTART WITH OLD VELOCITIES'
       IF (restart1%rac)   WRITE(6,'(A)') ' RESTART WITH OLD ACCUMULATORS'
       IF (restart1%rnoe)  WRITE(6,'(A)')&
            ' RESTART WITH OLD ELECTRON THERMOSTAT'
       IF (restart1%rnop)  WRITE(6,'(A)') ' RESTART WITH OLD ION THERMOSTAT'
       IF (restart1%rnoc)  WRITE(6,'(A)') ' RESTART WITH OLD CELL THERMOSTAT'
       IF (restart1%rgeo)  WRITE(6,'(A)') ' RESTART FROM FILE GEOMETRY'
       IF (restart1%rlate) WRITE(6,'(A)') ' RESTART WITH LATEST RESTART FILE'
       IF (restart1%rvib)  WRITE(6,'(A)')&
            ' RESTART FINITE DIFFERENCE CALCULATION'
       IF (restart1%rcell) WRITE(6,'(A)') ' RESTART WITH OLD MD CELL'
       IF (restart1%rpot)  WRITE(6,'(A)') ' RESTART WITH OLD POTENTIAL'
       IF (restart1%roh)   WRITE(6,'(A)') ' RESTART ORBITAL HARDNESS CALCULATION '
       IF (restart1%rlr)   WRITE(6,'(A)')&
            ' RESTART WITH OLD LINEAR RESPONSE VECTORS'
       IF (restart1%rphes) WRITE(6,'(A)') ' RESTART WITH OLD PARTIAL HESSIAN'
       IF (restart1%rlscl) WRITE(6,'(A)')&
            ' RESTART LINEAR SCALING OPTIMIZERS FROM OLD STATUS'
       IF (restart1%rdtl)  WRITE(6,'(A)')&
            ' RESTART WITH OLD CONVERGENCE CRITERIA'
       IF (restart1%rcon)  WRITE(6,'(A)') ' RESTART WITH OLD CONSTRAINT VALUES'
       IF (restart1%rxtp)  WRITE(6,'(A)') ' RESTART WITH WAVEFUNCTION HISTORY'
       IF (restart1%rprng)  WRITE(6,'(A)')&
            ' RESTART WITH PSEUDO RANDOM NUMBER GENERATOR'
       IF (restart1%rgle)  WRITE(6,'(A)')&
            ' RESTART WITH GENERALIZED LANGEVIN MOMENTA'
    ENDIF
    IF (cntl%tmdbo.OR.cntl%tmdfile.OR.(.NOT.cntl%md)) THEN
       IF (cntl%tlowd) THEN
          WRITE(6,'(A)') ' LOWDIN (SYMMETRIC) ORTHOGONALIZATION'
       ELSE
          WRITE(6,'(A)') ' GRAM-SCHMIDT ORTHOGONALIZATION'
       ENDIF
    ELSE
       IF (.NOT.cntl%nonort) THEN
          IF (cntl%cnstmd) THEN
             WRITE(6,'(A,/,A,T61,I5,/,A,T54,1PE12.2)')&
                  ' ITERATIVE ORTHOGONALIZATION',&
                  '    MAXIT:',cnti%maxit,&
                  '    EPS:',cntr%epsog
          ENDIF
       ENDIF
    ENDIF
    IF (tjlimit.GT.0._real_8) THEN
       WRITE(6,'(A,T54,F12.2)') ' MAXIMUM ELAPSED (RUN) TIME: ',tjlimit
    ENDIF
    WRITE(6,'(A,T50,I10,A6)')&
         ' MAXIMUM NUMBER OF STEPS: ',cnti%nomore,' STEPS'
    WRITE(6,'(A,T50,I10,A6)') ' MAXIMUM NUMBER OF ITERATIONS FOR SC:',&
         cnti%nomore_iter,' STEPS'
    IF (cntl%vibrat) THEN
       WRITE(6,'(A,1PE12.4)')&
            ' STEP SIZE FOR FINITE DIFFERENCES (BOHR):',cntr%fdiff
    ENDIF
    IF (cntl%krwfn) THEN
       WRITE(6,'(A)')&
            ' STORE REAL SPACE REPRESENTATION OF WAVEFUNCTIONS'
    ENDIF
    IF (cntl%wcomp) THEN
       WRITE(6,'(A)')&
            ' WRITE WAVEFUNCTIONS IN COMPRESSED FORM TO FILE'
       IF (cnti%wcompb.LT.0) THEN
          WRITE(6,'(A)') '   USE ATOMIC ORBITAL PROJECTION'
       ELSE
          WRITE(6,'(A,I3)') '   COMPRESSION FACTOR IS ',cnti%wcompb
       ENDIF
    ENDIF
    WRITE(6,'(A,T50,I10,A6)')&
         ' PRINT INTERMEDIATE RESULTS EVERY ',cprint%iprint_step,' STEPS'
    WRITE(6,'(A,T50,I10,A6)')&
         ' STORE INTERMEDIATE RESULTS EVERY ',store1%istore,' STEPS'
    IF (cntl%md.OR.cntl%tprcp.OR.cntl%geopt) THEN
       WRITE(6,'(A,T38,I6,A22)')&
            ' STORE INTERMEDIATE RESULTS EVERY ',store1%isctore,&
            ' SELF-CONSISTENT STEPS'
    ENDIF
    IF (.NOT.store1%swf) THEN
       WRITE(6,'(A)')&
            ' DO NOT STORE WAVEFUNCTIONS IN RESTART FILE'
    ENDIF
    IF (store1%srho) THEN
       WRITE(6,'(A)')&
            ' STORE ELECTRONIC DENSITY IN RESTART FILE'
    ENDIF
    WRITE(6,'(A,T56,I10)')&
         ' NUMBER OF DISTINCT RESTART FILES:',cnti%nrestf
    IF (restf%nstepwr(1).GE.0) THEN
       WRITE(6,'(A)') '   RESTART FILES WRITTEN AT STEPS'
       WRITE(6,'(10(8X,8I8))') (restf%nstepwr(i),i=1,cnti%nrestf)
    ENDIF
    IF (cntl%md.AND.cntl%tpres) THEN
       WRITE(6,'(A,I3,A)')&
            ' EVALUATE THE STRESS TENSOR OF THE SYSTEM EVERY ',&
            cnti%npres,' STEPS'
    ENDIF
    IF (cntl%md.AND.cntl%tprec) THEN
       WRITE(6,'(A,I3,A)')&
            ' EVALUATE THE CLASSICAL STRESS TENSOR OF SYSTEM EVERY   ',&
            cnti%nprec,' STEPS'
    ENDIF
    IF (isos1%tisos) THEN
       WRITE(6,'(A)')&
            ' TEMPERATURE IS CALCULATED ASSUMING AN ISOLATED MOLECULE '
    ELSE
       WRITE(6,'(A)')&
            ' TEMPERATURE IS CALCULATED ASSUMING EXTENDED BULK BEHAVIOR '
    ENDIF
    IF (rout1%rhoout) THEN
       IF (rout1%nrhoout.GT.0)&
            WRITE(6,'(A,I6,A)') ' STORE ELECTRON DENSITY EVERY ',&
            rout1%nrhoout,' STEPS DURING THE RUN'
       WRITE(6,'(A)')&
            ' STORE ELECTRON DENSITY AT THE END OF THE RUN'
    ENDIF
    IF (cntl%tepot) THEN
       WRITE(6,'(A)') ' PRINT ELECTROSTATIC POTENTIAL'
    ENDIF
    IF (cntl%tsic) THEN             ! cmb_ssic
       WRITE(6,'(/,A)') ' SELF INTERACTION CORRECTION AFTER'
       WRITE(6,'(A)')&
            '   J. VANDEVONDELE AND M. SPRIK, PCCP 7, 1363 (2005)'
       WRITE(6,'(A,T54,F12.4)')'   ALPHA =',cntr%asic
       WRITE(6,'(A,T54,F12.4,/)')'   BETA  =',cntr%bsic
    ENDIF                     ! cmb_ssic
    IF (.NOT.cntl%tdiag) THEN
       WRITE(6,'(A,T54,F12.4)')&
            ' FICTITIOUS ELECTRON MASS: ',cntr%emass
    ENDIF
    IF (cntl%tsdc.OR.cntl%tprcp) THEN
       IF (cntr%cmass.LE.0.0) THEN
          WRITE(6,'(A,T57,A)') ' FICTITIOUS MD CELL MASS: ','AUTOMATIC'
       ELSE
       !  WRITE(6,'(A,T54,F12.4)') ' FICTITIOUS MD CELL MASS: ',cntr%cmass
          WRITE(6,'(A,T54,E12.6)') ' FICTITIOUS MD CELL MASS: ',cntr%cmass
       ENDIF
    ENDIF
    ! Time step for ions and electrons
    WRITE(6,'(A,T54,F12.4)') ' TIME STEP FOR ELECTRONS: ',cntr%delt_elec
    WRITE(6,'(A,T54,F12.4)') ' TIME STEP FOR IONS: ',cntr%delt_ions
    IF (cntl%trane) THEN
       WRITE(6,'(A,T54,F12.4)')&
            ' RANDOMIZE INITIAL WAVEFUNCTION,  AMPLITUDE=',cntr%ampre
    ENDIF
    IF (cntl%tranp) THEN
       WRITE(6,'(A,T54,F12.4)')&
            ' RANDOMIZE INITIAL IONIC POSITIONS,  AMPLITUDE=',cntr%amprp
    ENDIF
    IF (cntl%tranc) THEN
       WRITE(6,'(A,T54,F12.4)')&
            ' RANDOMIZE INITIAL CELL PARAMETERS,  AMPLITUDE=',cntr%amprc
    ENDIF
    IF (andr2%trand) THEN
       WRITE(6,'(A,T54,F12.4)')&
            ' RANDOMIZE INITIAL DENSITY,  AMPLITUDE=',andr2%amprd
    ENDIF
    IF (group%nogrp.NE.1) THEN
       WRITE(6,'(A,T62,I4)')&
            ' NUMBER OF TASKGROUPS:',group%nogrp
       IF (cntl%tcart)  THEN
          WRITE(6,'(A)') ' USING CARTESIAN COMMUNICATORS '
       ENDIF
    ENDIF
    IF (rout1%acout) THEN
       WRITE(6,'(A,/,A)')&
            ' VIBRATIONAL FREQUENCIES AND EIGENVECTORS ARE',&
            '  SAVED IN a-CLIMAX FORMATTED FILE VIB.aclimax'
    ENDIF
    IF (rout1%vout) THEN
       IF (cnti%nvib.EQ.2) THEN
          WRITE(6,'(A,/,A)')&
               ' ALL VIBRATIONAL FREQUENCIES AND EIGENVECTORS ARE',&
               '  SAVED IN GAUSSIAN-STYLE FORMATTED FILES'
       ELSE
          WRITE(6,'(A,/,A)')&
               ' VIBRATIONAL FREQUENCIES AND EIGENVECTORS ARE',&
               '  SAVED IN THE GAUSSIAN-STYLE FORMATTED FILE VIB1.log'
       ENDIF
    ENDIF
    IF (cntl%md) THEN
       IF (cntl%tshop.AND.cntl%tfusi) THEN
          WRITE(6,'(A)')&
               ' COMBINE RESTART FILES FOR S0 AND S1 STATES '
          CALL xstring(s0_filn,ia,ie)
          WRITE(6,'(A,T20,A)') '     S0 STATE:',s0_filn(ia:ie)
          CALL xstring(s1_filn,ia,ie)
          WRITE(6,'(A,T20,A)') '     S1 STATE:',s1_filn(ia:ie)
       ENDIF
       IF (cntl%tshop.AND.cntl%tsep) THEN
          WRITE(6,'(A)')&
               ' SEPARATE RESTART FILE FOR S0 AND S1 STATES '
          CALL xstring(s0_filn,ia,ie)
          WRITE(6,'(A,T20,A)') '     S0 STATE:',s0_filn(ia:ie)
          CALL xstring(s1_filn,ia,ie)
          WRITE(6,'(A,T20,A)') '     S1 STATE:',s1_filn(ia:ie)
       ENDIF
       IF (cntl%tdiag) THEN
          WRITE(6,'(A,T60,F6.2)')&
               ' ALEXANDER MIXING PARAMETER:',&
               andr2%alxmix
       ENDIF
       IF (cntl%quenchp) THEN
          WRITE(6,'(A)')&
               ' RESET INITIAL VELOCITIES OF THE IONS TO ZERO'
       ENDIF
       IF (cntl%quenche) THEN
          WRITE(6,'(A)')&
               ' RESET INITIAL VELOCITIES OF THE ELECTRONS TO ZERO'
       ENDIF
       IF (cntl%quenchb) THEN
          WRITE(6,'(A)')&
               ' QUENCH SYSTEM TO THE BORN-OPPENHEIMER SURFACE'
       ENDIF
       IF (cntl%quenchc) THEN
          WRITE(6,'(A)')&
               ' QUENCH INITIAL VELOCITIES OF THE CELL TO ZERO'
       ENDIF
       IF (cntl%trevers) THEN
          WRITE(6,'(A)')&
               ' IONIC AND ELECTRONIC VELOCITIES HAVE BEEN INVERTED'
       ENDIF
       IF (rout1%rout) THEN
          IF (cnti%ntraj.EQ.1) THEN
             WRITE(6,'(A)') ' TRAJECTORIES ARE SAVED ON FILE'
          ELSE
             WRITE(6,'(A,T56,I4,A6)')&
                  ' TRAJECTORIES ARE SAVED ON FILE EVERY',cnti%ntraj,' STEPS'
          ENDIF
       ENDIF
       IF (rout1%xtout) THEN
          IF (cnti%ntraj.EQ.1) THEN
             WRITE(6,'(A)') ' TRAJEC.xyz IS SAVED ON FILE'
          ELSE
             WRITE(6,'(A,T56,I4,A6)')&
                  ' TRAJEC.xyz IS SAVED ON FILE EVERY',cnti%ntraj,' STEPS'
          ENDIF
       ENDIF
       IF (rout1%dcout) THEN
          IF (cnti%ntraj.EQ.1) THEN
             WRITE(6,'(A)') ' TRAJEC.dcd IS SAVED ON FILE'
          ELSE
             WRITE(6,'(A,T56,I4,A6)')&
                  ' TRAJEC.dcd IS SAVED ON FILE EVERY',cnti%ntraj,' STEPS'
          ENDIF
       ENDIF
       IF (rout1%mout) THEN
          WRITE(6,'(A,I0,A)') ' TRAJECTORIES ARE SAVED '// &
               'EVERY ',cnti%imovie,' STEPS IN MOVIE FORMAT'
       ENDIF
       IF (cntl%annei) THEN
          WRITE(6,'(A,T46,A8,F12.6)')&
               ' SIMULATED ANNEALING OF IONS WITH ','ANNERI =',cntr%anneri
       ENDIF
       IF (cntl%annee) THEN
          WRITE(6,'(A,T46,A8,F12.6)')&
               ' SIMULATED ANNEALING OF ELECTRONS WITH ','ANNERE =',cntr%annere
       ENDIF
       IF (cntl%annec) THEN
          WRITE(6,'(A,T46,A8,F12.6)')&
               ' SIMULATED ANNEALING OF CELL WITH ','ANNERC =',cntr%annerc
       ENDIF
       IF (cntl%dampi) THEN
          WRITE(6,'(A,T46,A8,F12.6)')&
               ' DAMPED DYNAMICS OF IONS WITH ','DAMPGI =',cntr%dampgi
       ENDIF
       IF (cntl%dampe) THEN
          WRITE(6,'(A,T46,A8,F12.6)')&
               ' DAMPED DYNAMICS OF ELECTRONS WITH ','DAMPGE =',cntr%dampge
       ENDIF
       IF (cntl%dampc) THEN
          WRITE(6,'(A,T46,A8,F12.6)')&
               ' DAMPED DYNAMICS OF CELL WITH ','DAMPGC =',cntr%dampgc
       ENDIF
       IF (cntl%tr4a2a) THEN
          WRITE(6,'(A)')&
               ' SINGLE PRECISION IS USED IN SOME ALLTOALL COMMUNICATION '
       ENDIF
       IF (cntl%tc) THEN
          WRITE(6,'(A,/,A,T54,1PE12.6,/,A,T54,1PE12.6)')&
               ' ELECTRON DYNAMICS WITH RESCALING OF VELOCITIES',&
               '    AVERAGE KINETIK ENERGY(A.U.)',cntr%ekinw,&
               '    TOLERANCE',cntr%toll
          IF (cntl%trescale) THEN
             WRITE(6,'(A)')&
                  '   ***** BUT RESCALE THEM ONLY AFTER THE RESTART TO'
             WRITE(6,'(A,T47,F12.2,A)')&
                  '   ***** SPECIFIED TEMPERATURE: ',cntr%tempw,' KELVIN'
          ENDIF
          IF (cntl%tbere) THEN
             WRITE(6,'(A,/,A,T54,1PE12.6,/,A,T54,1PE12.6)')&
                  ' ELECTRONS DYNAMICS WITH BERENDSEN-STYLE THERMOSTAT',&
                  '    AVERAGE KINETIK ENERGY(A.U.)',cntr%ekinw,&
                  '    CHARACTERISTIC TIME(A.U.):',cntr%taube
          ENDIF
          IF (cntl%tprcp.AND.cntl%tcc) THEN
             WRITE(6,'(A,/,A,T54,1PE12.6,A,/,A,T54,1PE12.6)')&
                  ' CELL DYNAMICS WITH RESCALING OF VELOCITIES',&
                  '    EKINH =',cntr%ekinhr,' AU',&
                  '    TOLERANCE=',cntr%tolc
          ENDIF
       ELSEIF (cntl%tnosee) THEN
          WRITE(6,'(A,T21,A)') ' ELECTRON DYNAMICS:',&
               'TEMPERATURE CONTROL (NOSE-HOOVER THERMOSTATS)'
          WRITE(6,'(A,T54,1PE12.6,/,A,T54,0PF12.2)')&
               '    TARGET KINETIK ENERGY(A.U.):',cntr%ekinw,&
               '    CHARACTERISTIC FREQUENCY(CM**-1):',cntr%wnose0
       ELSEIF (cntl%tbere) THEN
          WRITE(6,'(A,/,A,T54,1PE12.6,/,A,T54,1PE12.6)')&
               ' ELECTRONS DYNAMICS WITH BERENDSEN-STYLE THERMOSTAT',&
               '    AVERAGE KINETIK ENERGY(A.U.)',cntr%ekinw,&
               '    CHARACTERISTIC TIME(A.U.):',cntr%taube
       ELSE
          WRITE(6,'(A,T21,A)')' ELECTRON DYNAMICS:',&
               'THE TEMPERATURE IS NOT CONTROLLED'
       ENDIF
       IF (cntl%tcp) THEN
          WRITE(6,'(A,/,A,T54,1PE12.6,/,A,T54,1PE12.6)')&
               ' ION DYNAMICS WITH RESCALING OF VELOCITIES',&
               '    TEMPERATURE(KELVIN):',cntr%tempw,&
               '    TOLERANCE:',cntr%tolp
          IF (cntl%tberp) THEN
             WRITE(6,'(A,/,A,T54,1PE12.6,/,A,T54,1PE12.6)')&
                  ' ION DYNAMICS WITH BERENDSEN-STYLE THERMOSTAT',&
                  '    TEMPERATURE(KELVIN):',cntr%tempw,&
                  '    CHARACTERISTIC TIME(A.U.):',cntr%taubp
          ENDIF
       ELSEIF (cntl%tnosep) THEN
          IF (tcafes) THEN
             WRITE(6,*) '********** CAFES *************'
             WRITE(6,*) 'CANONICAL, ADIABATIC FREE ENERGY SAMPLING'
          ELSE
             WRITE(6,'(A,T21,A)')' ION DYNAMICS:',&
                  'TEMPERATURE CONTROL (NOSE-HOOVER THERMOSTATS)'
          ENDIF
          IF (nosl%tultra) THEN
             WRITE(6,'(A)') ' ONE THERMOSTAT CHAIN PER SPECIES'
          ENDIF
          IF (nosl%tmnose) THEN
             WRITE(6,'(A)') ' ONE THERMOSTAT CHAIN PER DEGREE OF FREEDOM'
          ENDIF
          IF (tcafes) THEN
             WRITE(6,*) 'TEMPERATURES FOR THE IONS TO BE DEFINED...'
             WRITE(6,*) 'FREQUENCIES OF THE THERMOSTATS TOO...'
          ELSEIF (loct%tloct) THEN
             WRITE(6,'(A)') ' THERMOSTAT FOR SUBGROUPS'
             DO i=1,loct%nloct
                WRITE(6,'(1X,I5,'':  '',A,T54,1PE12.6,/,1X,I5,'':  '',A,T54,0PF12.2)')&
                     i,'    TARGET TEMPERATURE(KELVIN):',LOCTPIN(1,i),&
                     i,'    CHARACTERISTIC FREQUENCY(CM**-1):',LOCTPIN(2,i)
             ENDDO
             WRITE(6,'(/,A,F10.3,/)')&
                  ' WARNING: INITIAL TEMPERATURE TAKEN AS ',cntr%tempw

          ELSE
             WRITE(6,'(A,T54,1PE12.6,/,A,T54,0PF12.2)')&
                  '    TARGET TEMPERATURE(KELVIN):',cntr%tempw,&
                  '    CHARACTERISTIC FREQUENCY(CM**-1):',cntr%wnosp0
          ENDIF
       ELSEIF (cntl%tberp) THEN
          WRITE(6,'(A,/,A,T54,1PE12.6,/,A,T54,1PE12.6)')&
               ' ION DYNAMICS WITH BERENDSEN-STYLE THERMOSTAT',&
               '    TEMPERATURE(KELVIN):',cntr%tempw,&
               '    CHARACTERISTIC TIME(A.U.):',cntr%taubp
       ELSE
          WRITE(6,'(A,T21,A)')' ION DYNAMICS:',&
               'THE TEMPERATURE IS NOT CONTROLLED'
       ENDIF
       IF (cntl%tnosec.AND.cntl%tprcp) THEN
          WRITE(6,'(A,T21,A)') ' CELL DYNAMICS:',&
               'TEMPERATURE CONTROL (NOSE-HOOVER THERMOSTATS)'
          WRITE(6,'(A,T54,1PE12.6,/,A,T54,0PF12.2)')&
               '    TARGET CELL TEMPERATURE(KELVIN):',cntr%tempc,&
               '    CHARACTERISTIC FREQUENCY(CM**-1):',cntr%wnosc0
       ELSEIF (cntl%tberc.AND.cntl%tprcp) THEN
          WRITE(6,'(A,/,A,T54,1PE12.6,/,A,T54,1PE12.6)')&
               ' CELL DYNAMICS WITH BERENDSEN-STYLE THERMOSTAT',&
               '    AVERAGE KINETIK ENERGY(A.U.)',cntr%ekinhr,&
               '    CHARACTERISTIC TIME(A.U.):',cntr%taubc
       ELSEIF (cntl%tcc.AND.cntl%tprcp) THEN
          WRITE(6,'(A,/,A,T54,1PE12.6,A,/,A,T54,1PE12.6)')&
               ' CELL DYNAMICS WITH RESCALING OF VELOCITIES',&
               '    EKINH =',cntr%ekinhr,' AU',&
               '    TOLERANCE=',cntr%tolc
       ELSEIF (cntl%tprcp) THEN
          WRITE(6,'(A,T21,A)') ' CELL DYNAMICS:',&
               'THE TEMPERATURE IS NOT CONTROLLED'
       ENDIF
       IF (cntl%tnosee.OR.cntl%tnosep.OR.cntl%tnosec) THEN
          WRITE(6,'(A)') ' NOSE PARAMETERS'
          WRITE(6,'(A,T63,I3)')&
               '    NUMBER OF THERMOSTATS (IONS)      :',cnti%nchp
          WRITE(6,'(A,T63,I3)')               &
               '    NUMBER OF THERMOSTATS (ELECTRONS) :',cnti%nchs
          WRITE(6,'(A,T63,I3)')               &
               '    NUMBER OF THERMOSTATS (CELL)      :',cnti%nchb
          WRITE(6,'(A,T60,F6.2)')               &
               '    SCALING FOR ELEC. DOF             :',cntr%nedof0
          WRITE(6,'(A,T63,I3)')               &
               '    NUMBER OF YOSHIDA-SUZUKI STEPS    :',cnti%ncalls0
          WRITE(6,'(A,T63,I3)')               &
               '    NUMBER OF INTEGRATION CYCLES (NIT):',cnti%nit0
       ENDIF
       IF ((cntl%tcp.OR.cntl%tnosep.OR.cntl%tberp).AND.(cntr%trampr.GT.0.0_real_8)) THEN
          WRITE(6,'(A,/,A,T54,1PE12.6,/,A,T54,1PE12.6)')&
               ' USING LINEAR TEMPERATURE RAMP ',&
               '    TARGET REFERENCE TEMPERATURE (KELVIN)',cntr%trampt,&
               '    RAMPING RATE (KELVIN/A.U.):',cntr%trampr
       ENDIF
       IF (comvl%tsubcom) THEN
          WRITE(6,'(A,T50,I10,A6)')&
               ' SUBTRACT CENTER OF MASS VELOCITY EVERY ',comvl%ncomv,' STEPS'
       ENDIF
       IF (comvl%tsubrot) THEN
          WRITE(6,'(A,A,T50,I10,A6)') ' SUBTRACT ROTATION AROUND',&
               ' CENTER OF MASS EVERY ',comvl%nrotv,' STEPS'
       ENDIF
       IF (cntl%tdipd) THEN
          WRITE(6,'(A)') ' DIPOLE MOMENT CALCULATION'
          WRITE(6,'(A,T50,I10,A6)')&
               '     STORE DIPOLE MOMENTS EVERY ',cnti%npdip,' STEPS'
          IF (wannl%twann) THEN
             WRITE(6,'(A)') '     WANNIER FUNCTION DYNAMICS'
             IF (wanni%w_type.EQ.1)&
                  WRITE(6,'(A)')   '         FUNCTIONAL: VANDERBILT |M|^2'
             IF (wanni%w_type.EQ.2)&
                  WRITE(6,'(A)')   '         FUNCTIONAL: RESTA Log|M|^2'
             IF (wanni%w_opt.EQ.1)&
                  WRITE(6,'(A)') '         OPTIMISATION: STEEPEST DESCENT'
             IF (wanni%w_opt.EQ.2)&
                  WRITE(6,'(A)') '         OPTIMISATION: JACOBI ROTATION'
             WRITE(6,'(A,T56,1PE10.4)')&
                  '         CONVERGENCE CRITERIA:',wannr%w_eps
             WRITE(6,'(A,T56,I10)')&
                  '         MAXIMUM # ITERATIONS:',wanni%w_maxs
             WRITE(6,'(A,T56,1PE10.4)')&
                  '         OPTIMISATION STEP   :',wannr%w_step
             IF (wannr%w_ref(1).LT.-999999.9_real_8) THEN
                WRITE(6,'(A)') '         REF POINT: AUTOMATIC'
             ELSE
                WRITE(6,'(A,3F11.5)')&
                     '         REF POINT[input units]:',(wannr%w_ref(i),i=1,3)
             ENDIF
             IF (cntl%twserial) THEN
                WRITE(6,'(A,3F11.5)')&
                     '         CALCULATE WANNIER CENTERS USING SERIAL CODE'
             ENDIF
             IF (wannl%twpri) THEN
                WRITE(6,'(A,T56,I4,A)')&
                     '     PRINT WANNIER FUNCTIONS EVERY ',&
                     wanni%ns_wann*cnti%npdip,' STEPS'
                IF (wanni%sw_all.EQ.1) THEN
                   WRITE(6,'(T20,A)') 'PRINT ALL WANNIER FUNCTIONS'
                ENDIF
                IF (ABS(wannr%sw_spread).GT.0.001_real_8) THEN
                   IF (wannr%sw_spread.LT.0.0_real_8) THEN
                      WRITE(6,'(T20,A,A,F6.2)') 'PRINT WANNIER FUNCTIONS ',&
                           'WITH A SPREAD SMALLER THAN ', -wannr%sw_spread
                   ELSE
                      WRITE(6,'(T20,A,A,F6.2)') 'PRINT WANNIER FUNCTIONS ',&
                           'WITH A SPREAD LARGER THAN ', wannr%sw_spread
                   ENDIF
                ENDIF
                IF (wanni%sw_first.GT.0) THEN
                   WRITE(6,'(T20,A,I4,A,I4)') 'PRINT WANNIER FUNCTIONS ',&
                        wanni%sw_first,' TROUGH ',wanni%sw_last
                ENDIF
                IF (wanni%sw_orb.GT.0) THEN
                   WRITE(6,'(T20,10I5)') (sw_list(i),i=1,wanni%sw_orb)
                ENDIF
             ENDIF
          ENDIF
       ENDIF
    ELSE
       WRITE(6,'(A,T54,1PE12.4)')&
            ' CONVERGENCE CRITERIA FOR WAVEFUNCTION OPTIMIZATION: ',&
            cntr%tolog
       IF (cntr%toldetot.GT.0._real_8) THEN
          WRITE(6,'(A,T54,1PE12.4)')&
               '    TOTAL ENERGY DIFFERENCE: ', cntr%toldetot
       ENDIF
       IF (cntr%tolad.GT.0._real_8) THEN
          WRITE(6,'(A,/,A,T54,1PE12.4)')&
               ' CONVERGENCE CRITERIA FOR WAVEFUNCTION RELATIVE TO ',&
               ' GRADIENT ON IONS: ', cntr%tolad
       ENDIF
       IF (cntr%tolene.GT.0._real_8) THEN
          WRITE(6,'(A,/,A,T54,1PE12.4)')&
               ' CONVERGENCE CRITERIA FOR RELAVTIVE WAVEFUNCTION ',&
               ' ENERGY CHANGE: ', cntr%tolene
       ENDIF
       IF (cntr%tolfor.NE.3._real_8) THEN
          WRITE(6,'(A,/,A,T54,1PE12.4)')&
               ' FORCE ON IONS CALCULATED IF CONVERGENCE IS REACHED ',&
               ' RELAXED BY: ', cntr%tolfor
       ENDIF
       IF (cnti%nstcnv.NE.0) THEN
          WRITE(6,'(A,/,A,T63,1I3)')&
               ' NUMBER OF STEPS UNTIL ELECTRONIC CONVERGENCE CRITERIA ',&
               ' ARE RELAXED: ', cnti%nstcnv
       ENDIF
       IF (cntl%diis) THEN
          IF (cntl%prec) THEN
             WRITE(6,'(A,/,A,T58,F8.4)')&
                  ' WAVEFUNCTION OPTIMIZATION BY PRECONDITIONED DIIS',&
                  ' THRESHOLD FOR THE WF-HESSIAN IS ',cntr%hthrs
          ELSE
             WRITE(6,'(A)')&
                  ' WAVEFUNCTION OPTIMIZATION BY DIIS'
          ENDIF
          WRITE(6,'(A,A,T63,I3)') ' MAXIMUM NUMBER OF VECTORS',&
               ' RETAINED FOR DIIS: ',cnti%mdiis
          WRITE(6,'(A,A,T63,I3)') ' STEPS UNTIL DIIS RESET ON POOR',&
               ' PROGRESS: ', cnti%nreset
       ENDIF
       IF (cntl%pcg) THEN
          IF (cntl%prec) THEN
             WRITE(6,'(A,/,A,T58,F8.4)')&
                  ' WAVEFUNCTION OPTIMIZATION BY PRECONDITIONED CG',&
                  ' THRESHOLD FOR THE HESSIAN IS ',cntr%hthrs
          ELSE
             WRITE(6,'(A)')&
                  ' WAVEFUNCTION OPTIMIZATION BY CG'
          ENDIF
          IF (cntl%pcgmin) WRITE(6,'(A)') ' PERFORM QUADRATIC LINE SEARCH'
       ENDIF
       IF (cntl%tsde) THEN
          WRITE(6,'(A)')&
               ' WAVEFUNCTION OPTIMIZATION BY STEEPEST DESCENT'
          IF (cntl%prec) THEN
             WRITE(6,'(A,/,A,T58,F8.4)')&
                  '           USING PRECONDITIONED MASSES',&
                  ' THRESHOLD FOR THE WF-HESSIAN IS ',cntr%hthrs
          ENDIF
       ENDIF
       IF (cntl%tsdc) THEN
          WRITE(6,'(A)')&
               ' CELL OPTIMIZATION BY STEEPEST DESCENT'
       ENDIF
       IF (.NOT.cntl%tdiag) THEN
          IF (cnti%iproj.EQ.0) WRITE(6,'(A)')&
               ' RAW ELECTRONIC GRADIENT (=H*PSI) IS USED'
          IF (cnti%iproj.EQ.1) WRITE(6,'(A)')&
               ' DIAGONALLY PROJECTED ELECTRONIC GRADIENT IS USED'
          IF (cnti%iproj.EQ.2) WRITE(6,'(A)')&
               ' FULL ELECTRONIC GRADIENT IS USED '
       ENDIF
    ENDIF
    IF (cntl%tdiagopt) THEN
       WRITE(6,'(A)')' DIAGONALIZE OPTIMIZED WAVEFUNCTION'
       IF (cnti%nperid.GT.0) WRITE(6,'(A,T63,I3)')&
            ' STEPS UNTIL WF DIAGONALIZATION: ',cnti%nperid
       IF (cnti%nrperi.GT.0) WRITE(6,'(A,T63,I3)')&
            ' DIIS RESETS UNTIL WF DIAGONALIZATION: ',cnti%nrperi
    ENDIF
    IF (cntl%tdiag.OR.cntl%tdiagopt) THEN
       IF (cntl%tlanc) THEN
          WRITE(6,'(A)')&
               ' LANCZOS DIAGONALIZATION (KRYLOV SUBSPACE)'
          WRITE(6,'(2(A,T63,I3,/))', advance='no')&
               '    MAX. FRIESNER ITERATIONS ',cnti%n_fries,&
               '    MAX. KRYLOV SUBSPACE     ',cnti%nkry_max
          IF (cnti%nkry_block.LE.0) THEN
             WRITE(6,'(A,T34,A32)')&
                  '    MAX. KRYLOV BLOCK SIZE',&
                  'FIXED LATER (WAITING FOR NSTATE)'
          ELSE
             WRITE(6,'(A,T63,I3)')&
                  '    MAX. KRYLOV BLOCK SIZE',cnti%nkry_block
          ENDIF
          WRITE(6,'(A,T54,1PE12.4)')&
               '    MAX. BETA^2              ',cntr%b2limit
          IF (fint4%ntabtrot.NE.1) THEN
             WRITE(6,'(A,T64,I2)')&
                  '    MAX. BETA^2 CHOICE NUMBER:',fint4%ntabtrot
             DO i=2,fint4%ntabtrot
                WRITE(6,'(A,1PE12.4,T42,A12,1PE12.4)')&
                     '    BELOW DRHOMAX=',fint4%denstrot(i),&
                     'MAX. BETA^2=',fint4%b2trot(i)
             ENDDO
          ENDIF
       ELSEIF (cntl%tdavi) THEN
          WRITE(6,'(A)') ' DAVIDSON DIAGONALIZATION'
          WRITE(6,'(A,T63,I3,/,A,T54,1PE12.4)')&
               '    MAX. DAVIDSON ITERATIONS ',cnti%n_fries,&
               '    CONVERGENCE CRITERIA     ',cntr%epsdav
       ELSEIF (cntl%tfrho_upw) THEN
          IF (cntl%diis) THEN
             WRITE(6,'(A)') ' WFN OPT. WITH DIIS AT FIXED RHO'
             WRITE(6,'(10x,A,I5)') '   MAX # of DIIS Cycles = ',cnti%maxldiis
             WRITE(6,'(10x,A)')&
                  '  DENSITY UPDATED AFTER EVERY DIIS LOOP'
          ELSE IF (cntl%pcg) THEN
             WRITE(6,'(A)') ' WFN OPT. WITH PCG AT FIXED RHO'
             WRITE(6,'(10x,A,I5)') '   MAX # OF PCG CYCLES = ',cnti%maxldiis
             WRITE(6,'(10x,A)')&
                  '  DENSITY UPDATED AFTER EVERY PCG LOOP'
          ENDIF
       ENDIF
       IF (cntl%md.OR.cntl%geopt.OR.cntl%wfopt.OR.cntl%tdebfor) THEN
          IF (andr3%nrdiismax.NE.0) THEN
             WRITE(6,'(A,T62,I4)')&
                  ' DIIS ON DENSITIES: MAX. SIZE',andr3%nrdiismax
          ENDIF
          IF (andr3%ntabnrdiis.NE.1) THEN
             WRITE(6,'(A,T64,I2)')&
                  ' DIIS MIXING NUMBER:',andr3%tabnrdiis(1)
             DO i=2,andr3%ntabnrdiis
                WRITE(6,'(4X,A,1PE12.4,T52,A10,I4)')&
                     'BELOW DRHOMAX=',andr3%densnrdiis(i),&
                     'MAX. SIZE=',andr3%tabnrdiis(i)
             ENDDO
          ENDIF
          WRITE(6,'(A,T54,1PE12.4)')&
               ' ANDERSON MIXING PARAMETER:',andr2%andrmix
          IF (andr2%ntabmix.NE.1) THEN
             WRITE(6,'(A,T64,I2)')&
                  ' ANDERSON MIXING NUMBER:',andr2%ntabmix
             DO i=2,andr2%ntabmix
                WRITE(6,'(4X,A,D12.4,T47,A7,1PE12.4)')&
                     'BELOW DRHOMAX=',andr2%densmix(i),'MIXING=',andr2%andmix(i)
             ENDDO
          ENDIF
          IF (broy1%tgmix) THEN
             WRITE(6,'(A)')' G-SPACE ANDERSON MIXING'
          ENDIF
          IF (broy1%tgbroy) THEN
             WRITE(6,'(A,T54,1PE12.4)')&
                  ' BROYDEN MIXING PARAMETER [BROYMIX]',broy1%broymix
             IF (broy1%ecutbroy.LE.0._real_8) THEN
                WRITE(6,'(A,T39,A27)')&
                     ' BROYDEN CUTOFF [ECUTBROY]',&
                     'EQUAL TO THE DENSITY CUTOFF'
             ELSE
                WRITE(6,'(A,T54,1PE12.4)')&
                     ' BROYDEN CUTOFF [ECUTBROY]',broy1%ecutbroy
             ENDIF
             WRITE(6,'(A,T57,I3,A)')&
                  ' BROYDEN MIXING STARTS [NFRBROY] AFTER ',&
                  broy1%nfrbroy, ' STEPS'
             WRITE(6,'(A,T57,I3,A)')&
                  ' BROYDEN MIXING RESET [IBRESET] AFTER ',&
                  broy1%ibreset, ' STEPS'
             WRITE(6,'(A,T54,1PE12.4)')&
                  ' BROYDEN MIXING W02',broy1%w02broy
             WRITE(6,'(A,T54,1PE12.4)')&
                  ' KERKER MIXING PARAMETER [KERMIX]',broy1%kermix
          ENDIF
          IF (tmovr) THEN
             WRITE(6,'(A,T54,F12.4)')&
                  ' MOVE DENSITY WITH ATOMIC MOTION, MIXING PARAMETER',&
                  dmovrmiX
          ELSE
             WRITE(6,'(A,T54,F12.4)')&
                  ' ALEXANDER MIXING:',andr2%alxmix
          ENDIF
       ENDIF
    ENDIF
    IF ((cntl%tmdbo.OR.cntl%tmdfile).AND.cntl%textrap) THEN
       IF (cntl%taspc) THEN
          WRITE(6,'(A,T60,I6)') ' ASPC WAVEFUNCTION '&
               // 'EXTRAPOLATION (MAX. HISTORY)',cnti%mextra
          IF (cnti%naspc.GT.0) WRITE(6,'(A,T60,I6)')&
               ' NUMBER OF ASPC CORRECTOR STEPS: ',cnti%naspc
       ELSE
          WRITE(6,'(A,T60,I6)') ' POLYNOMIAL WAVEFUNCTION '&
               // 'EXTRAPOLATION (MAX. HISTORY)',cnti%mextra
       ENDIF
    ENDIF
    IF (cntl%geopt) THEN
       WRITE(6,'(A,T54,1PE12.6)')&
            ' CONVERGENCE CRITERIA FOR GEOMETRY OPTIMIZATION:',&
            cntr%tolng
       IF (cntl%tsdp) THEN
          IF (sdpl%tcgp) THEN
             WRITE(6,'(A)')&
                  ' GEOMETRY OPTIMIZATION BY CONJUGATE GRADIENT'
          ELSE
             WRITE(6,'(A)')&
                  ' GEOMETRY OPTIMIZATION BY STEEPEST DESCENT'
             IF (sdpl%tsdpline) WRITE(6,'(A)')&
                  ' WITH MINIMISATION ALONG LINES'
          ENDIF
       ELSEIF (cntl%bfgs) THEN
          WRITE(6,'(A)')&
               ' GEOMETRY OPTIMIZATION BY QUASI-NEWTON UPDATE'
       ELSEIF (cntl%lbfgs) THEN
          WRITE(6,'(A)')&
               ' GEOMETRY OPTIMIZATION BY LOW-MEMORY BFGS'
       ELSEIF (cntl%prfo) THEN
          WRITE(6,'(A)')&
               ' TRANSITION STATE SEARCH BY PARTITIONED RFO'
          IF (nsvib.GT.0) WRITE(6,'(A,T53,I7,A6)')&
               ' VIBRATIONAL ANALYSIS EVERY ', nsvib, ' STEPS'
          IF (nsvib.EQ.-1) WRITE(6,'(A)')&
               ' VIBRATIONAL ANALYSIS DONE AT END ONLY'
       ELSEIF (cntl%rfo) THEN
          WRITE(6,'(A)')&
               ' GEOMETRY OPTIMIZATION BY RATIONAL FUNCTION OPTIMIZATION'
          IF (cnti%nsorder.NE.0) THEN
             WRITE(6,'(A,T62,I4)')&
                  '     SADDLE POINT ORDER:',cnti%nsorder
          ENDIF
       ELSEIF (cntl%gdiis) THEN
          WRITE(6,'(A,/,A,T62,I4)')&
               ' GEOMETRY OPTIMIZATION BY GDIIS/BFGS      ',&
               '   SIZE OF GDIIS MATRIX: ',cnti%mgdiis
       ELSEIF (cntl%tinr) THEN
          WRITE(6,'(A)')&
               ' GEOMETRY OPTIMIZATION BY IMPLICIT NEWTON-RAPHSON '&
               //' (USES DFPT) '
          IF (inr_logical%inr_prec) THEN
             WRITE(6,'(A,I4)')&
                  ' PRECONDITIONED --  MAXIMUM NUMBER OF INTERNAL '//&
                  'INR ITERATIONS', inr_integer%itmax_inr
          ELSE
             WRITE(6,'(A,19X,I4)')&
                  ' MAXIMUM NUMBER OF INTERNAL INR ITERATIONS',&
                  inr_integer%itmax_inr
          ENDIF
          WRITE(6,'(a,35x,i4)') ' NUMBER OF GRADIENT ZONES ',inr_integer%nreg
          DO i=1,inr_integer%nreg
             WRITE(6,'(A,E10.4,3X,A,E10.4)') ' FOR GNMAX BIGGER THAN ',&
                  gnx_inr(i),' RESIDUAL CONV. AT ',tolx_inr(i)
          ENDDO
          IF (inr_logical%tmixsd) THEN
             WRITE(6,'(1X,A)')'INITIAL OPTIMIZATION WITH IONIC STEEPEST'&
                  //' DESCENT'
             WRITE(6,'(A,T54,E12.6)')&
                  ' INR STARTING AFTER GNMAX = ',rmixsd
          ELSEIF (inr_logical%tmixgdiis) THEN
             WRITE(6,'(1X,A)') 'INITIAL OPTIMIZATION WITH GDIIS'
             WRITE(6,'(A,T54,E12.6)')&
                  ' INR STARTING AFTER GNMAX = ',rmixsd
          ENDIF
       ENDIF
       IF (cntl%tsdc) THEN
          WRITE(6,'(A)')&
               ' CELL OPTIMIZATION BY STEEPEST DESCENT'
       ENDIF
       IF (cntl%simul) THEN
          WRITE(6,'(A)')&
               ' COMBINED GEOMETRY/WAVEFUNCTION SCHEME '
       ENDIF
       IF (rout1%xgout) THEN
          IF (cnti%ngxyz.EQ.1) THEN
             WRITE(6,'(A,A)') 'GEOMETRY OPTIMIZATION IS SAVED ON FILE ',&
                  'GEO_OPT.xyz'
          ELSE
             WRITE(6,'(A,T56,I4,A6)')&
                  ' GEOMETRY OPTIMIZATION IS SAVED ON FILE EVERY',cnti%ngxyz,' STEPS'
          ENDIF
       ENDIF
    ENDIF
    IF (cntl%geopt.OR.(cntl%vibrat.AND.cntl%tsdin)) THEN
       IF (restart1%rhe) THEN
          WRITE(6,'(A)') ' RESTART WITH OLD HESSIAN'
       ELSE
          IF (cnti%npara.EQ.-1) THEN
             WRITE(6,'(A)')&
                  ' INITIAL HESSIAN IS UNIT MATRIX '
          ELSEIF (cnti%npara.EQ.0) THEN
             WRITE(6,'(A)')&
                  ' EMPIRICAL INITIAL HESSIAN (DISCO PARAMETRISATION) '
          ELSEIF (cnti%npara.EQ.1) THEN
             WRITE(6,'(A)')&
                  ' EMPIRICAL INITIAL HESSIAN (SCHLEGEL PARAMETRISATION) '
          ELSEIF (cnti%npara.EQ.2) THEN
             WRITE(6,'(A)')&
                  ' INITIAL HESSIAN FROM UNIT AND TS SEARCH PARTIAL HESSIAN'
          ELSEIF (cnti%npara.EQ.3) THEN
             WRITE(6,'(A)')&
                  ' INITIAL HESSIAN FROM DISCO AND TS SEARCH PARTIAL HESSIAN'
          ELSEIF (cnti%npara.EQ.4) THEN
             WRITE(6,'(A)')&
                  ' INITIAL HESSIAN FROM SCHLEGEL AND TS SEARCH PARTIAL HESSIAN'
          ENDIF
       ENDIF
    ENDIF
    IF (vdwl%vdwc) THEN
       WRITE(6,'(A)') ' USE EMPIRICAL VAN DER WAALS CORRECTION'
    ENDIF
    IF (vdwl%vdwd) THEN
       WRITE(6,'(A,A)') ' USE VAN DER WAALS CORRECTION BASED ON',&
            ' WANNIER FUNCTIONS/CENTERS'
    ENDIF
    IF (cntl%tdmal) THEN
       WRITE(6,'(A)') ' USE PARALLEL MATRIX ALGEBRA CODE'
    ENDIF
    IF (qspl1%initspl) THEN
       WRITE(6,'(A,A)')&
            ' SPLINE INTERPOLATION IN G-SPACE ',&
            'FOR PSEUDOPOTENTIAL FUNCTIONS'
    ENDIF
    IF (qspl1%qspline) THEN
       WRITE(6,'(A)')&
            ' SPLINE INTERPOLATION OF Q-FUNCTIONS '
    ENDIF
    IF (qspl1%initspl.OR.qspl1%qspline) THEN
       WRITE(6,'(A,T56,I10)')&
            '    NUMBER OF SPLINE POINTS: ',cnti%nsplp
    ENDIF
    IF (cntl%tfdist) THEN
       WRITE(6,'(A)')&
            ' FNL ARRAY DISTRIBUTED OVER PROCESSORS'
    ENDIF
    IF (iface1%intread) THEN
       WRITE(6,'(A)')&
            ' READ INFORMATION FROM INTERFACE FILE '
    ENDIF
    IF (iface1%intwrite) THEN
       WRITE(6,'(A)')&
            ' WRITE INFORMATION TO INTERFACE FILE '
    ENDIF
    IF (iface1%intwrite.OR.iface1%intread) THEN
       WRITE(6,'(A,/,A)')&
            ' INTERFACE FILE NAME ',intfN
    ENDIF
    IF (iunk.NE.0) THEN
       WRITE(6,'(/,1X,64("="))')
       WRITE(6,'(1X,A,14X,A,14X,A)') '= ',&
            'UNKNOWN KEYWORDS IN SECTION &CPMD','='
       DO i=1,iunk
          WRITE(6,'(1X,A,A,A)') '= ',unknown(i),' ='
       ENDDO
       WRITE(6,'(1X,64("="),/)')
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE control_pri
  ! ==================================================================

END MODULE control_pri_utils
