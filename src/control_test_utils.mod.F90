#include "cpmd_global.h"

MODULE control_test_utils
  USE broy,                            ONLY: broy1
  USE cdftmod,                         ONLY: cdftci,&
       cdftlog
  USE comvelmod,                       ONLY: comvl
  USE cp_cuda_types,                   ONLY: cp_cuda_env
  USE error_handling,                  ONLY: stopgm
  USE fileopenmod,                     ONLY: fo_info
  USE fint,                            ONLY: fint1,&
       fint4
  USE kpts,                            ONLY: tkpts
  USE linres,                          ONLY: tshl,&
       xfmqc
  USE lscal,                           ONLY: nsmaxp
  USE mm_input,                        ONLY: clc,&
       lqmmm
  USE mw,                              ONLY: mwi,&
       tmw                                                             
  USE parac,                           ONLY: parai
  USE qspl,                            ONLY: qspl1
  USE spin,                            ONLY: lspin2
  USE store_types,                     ONLY: cprint,&
       restart1,&
       store1
  USE system,                          ONLY: cnti,&
       cntl,&
       cntr,&
       group,&
       maxdis,&
       mxgdis
  USE vdwcmod,                         ONLY: vdwl
  USE wann,                            ONLY: real_8,&
       wannl
  !$ USE omp_lib, ONLY: omp_get_max_threads

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: control_test

CONTAINS

  ! ==================================================================
  SUBROUTINE control_test(tsrho)
    ! ==--------------------------------------------------------------==
    ! ==  Test of variables read in control routine                   ==
    ! ==--------------------------------------------------------------==
    LOGICAL                                  :: tsrho

    CHARACTER(len=*), PARAMETER              :: procedureN = 'control_test'

    INTEGER                                  :: it
    LOGICAL                                  :: exist

    ! Variables
    ! ==--------------------------------------------------------------==
    ! ==  Test of options                                             ==
    ! ==--------------------------------------------------------------==

    IF (cntl%tmdfile.AND.cntl%wfopt) cntl%wfopt=.FALSE.
    it=0
    IF (cntl%md) THEN
       it=it+1
    ENDIF
    IF (cntl%wfopt) THEN
       it=it+1
    ENDIF
    IF (cntl%geopt) THEN
       it=it+1
    ENDIF
    IF (cntl%vibrat) THEN
       it=it+1
    ENDIF
    IF (cntl%ksener) THEN
       it=it+1
    ENDIF
    IF (cntl%proper) THEN
       it=it+1
    ENDIF
    IF (cntl%tinter) THEN
       it=it+1
    ENDIF
    IF (cntl%thard) THEN
       it=it+1
    ENDIF
    IF (cntl%tspec) THEN
       it=it+1
    ENDIF
    IF (cntl%tsoc) THEN
       it=it+1
    ENDIF
    ! EHR[
    IF (cntl%tpdist) THEN
       it=it+1
    ENDIF
    IF (cntl%tpspec) THEN
       it=it+1
    ENDIF
    ! EHR]
    IF (cntl%tdebfor) THEN
       it=it+1
    ENDIF
    IF (cntl%fmatch) THEN
       it=it+1
    ENDIF
    IF (it.LT.1) THEN
       WRITE(6,'(A,/A,/,A,/,A,/,A,/,A,/,A,/,A,/,A,/,A,/,A,/,A,/,A)')&
            ' CHOOSING ONE OF THE FOLLOWING OPTIONS IS COMPULSORY ',&
            '                                  MOLECULAR DYNAMICS ',&
            '                               GEOMETRY OPTIMIZATION ',&
            '                           WAVEFUNCTION OPTIMIZATION ',&
            '                                  KOHN-SHAM ENERGIES ',&
            '                                VIBRATIONAL ANALYSIS ',&
            '                                          PROPERTIES ',&
            '                                    ORBITAL HARDNESS ',&
            '                                  ELECTRONIC SPECTRA ',&
            '                                     LINEAR RESPONSE ',&
            '                                        DEBUG FORCES ',&
            '                                 CLASSICAL INTERFACE ',&
            '                                          FORCEMATCH ',&
            '                                 SPIN-ORBIT COUPLING '
       CALL stopgm('CONTROL','RUNOPTIONS ',& 
            __LINE__,__FILE__)
    ELSEIF (it.GT.1) THEN
       WRITE(6,'(A,/A,/,A,/,A,/,A,/,A,/,A,/,A,/,A,/,A,/,A,/,A,/,A)')&
            '       THE FOLLOWING OPTIONS ARE MUTUALLY  EXCLUSIVE ',&
            '                                  MOLECULAR DYNAMICS ',&
            '                               GEOMETRY OPTIMIZATION ',&
            '                           WAVEFUNCTION OPTIMIZATION ',&
            '                                  KOHN-SHAM ENERGIES ',&
            '                                VIBRATIONAL ANALYSIS ',&
            '                                          PROPERTIES ',&
            '                                    ORBITAL HARDNESS ',&
            '                                  ELECTRONIC SPECTRA ',&
            '                                        DEBUG FORCES ',&
            '                                     LINEAR RESPONSE ',&
            '                                 CLASSICAL INTERFACE ',&
            '                                          FORCEMATCH ',&
            '                                 SPIN-ORBIT COUPLING '
       CALL stopgm('CONTROL','RUNOPTIONS ',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Path calculations (integral or minimisation)
    IF (cntl%tpath) THEN
       IF (cnti%nrestf.GT.1) THEN
          WRITE(6,'(/,1X,8("WARNING! "),/,A,/,1X,8("WARNING!"))')&
               ' WARNING! RESETTING RESTFILE TO 1 FOR Path-INTEGRALS'
          cnti%nrestf=1
       ENDIF
       ! Path integral
       IF (cntl%tpimd) THEN
          IF (cntl%geopt) CALL stopgm('CONTROL_TEST','GEOMETRY OPTIMIZATION '&
               // 'WITH Path INTEGRALS IS NOT SUPPORTED',& 
               __LINE__,__FILE__)
          IF (cntl%ksener) CALL stopgm('CONTROL_TEST','KOHN-SHAM ENERGIES '&
               // 'WITH Path INTEGRALS IS NOT SUPPORTED',& 
               __LINE__,__FILE__)
          IF (cntl%vibrat) CALL stopgm('CONTROL_TEST','VIBRATIONAL ANALYSIS '&
               // 'WITH Path INTEGRALS IS NOT SUPPORTED',& 
               __LINE__,__FILE__)
          IF (cntl%tinter) CALL stopgm('CONTROL_TEST','CLASSICAL INTERFACE '&
               // 'WITH Path INTEGRALS IS NOT SUPPORTED',& 
               __LINE__,__FILE__)
          IF (cntl%tmdfile) CALL stopgm('CONTROL_TEST','FILE MD '&
               // 'WITH Path INTEGRALS IS NOT SUPPORTED',& 
               __LINE__,__FILE__)
          IF (cntl%tshop) CALL stopgm('CONTROL_TEST','SURFACE HOPPING '&
               // 'WITH Path INTEGRALS IS NOT SUPPORTED',& 
               __LINE__,__FILE__)
          IF (cntl%bsymm) CALL stopgm('CONTROL_TEST','BROKEN SYMMETRY '&
               // 'WITH Path INTEGRALS IS NOT SUPPORTED',& 
               __LINE__,__FILE__)
          IF (cntl%tqmmm) CALL stopgm('CONTROL_TEST','QM/MM '&
               // 'WITH Path INTEGRALS IS NOT SUPPORTED',& 
               __LINE__,__FILE__)
       ELSEIF (cntl%tpmin) THEN
          ! Path minimisation
          IF (cntl%ksener) CALL stopgm('CONTROL_TEST','KOHN-SHAM ENERGIES '&
               // 'WITH Path MINIMIZATION IS NOT SUPPORTED',& 
               __LINE__,__FILE__)
          IF (cntl%vibrat) CALL stopgm('CONTROL_TEST','VIBRATIONAL ANALYSIS '&
               // 'WITH Path MINIMIZATION IS NOT SUPPORTED',& 
               __LINE__,__FILE__)
          IF (cntl%tinter) CALL stopgm('CONTROL_TEST','CLASSICAL INTERFACE '&
               // 'WITH Path MINIMIZATION IS NOT SUPPORTED',& 
               __LINE__,__FILE__)
          IF (cntl%tmdfile) CALL stopgm('CONTROL_TEST','FILE MD '&
               // 'WITH Path MINIMIZATION IS NOT SUPPORTED',& 
               __LINE__,__FILE__)
          IF (cntl%tshop) CALL stopgm('CONTROL_TEST','SURFACE HOPPING '&
               // 'WITH Path MINIMIZATION IS NOT SUPPORTED',& 
               __LINE__,__FILE__)
          IF (cntl%bsymm) CALL stopgm('CONTROL_TEST','BROKEN SYMMETRY '&
               // 'WITH Path MINIMIZATION IS NOT SUPPORTED',& 
               __LINE__,__FILE__)
          IF (cntl%tqmmm) CALL stopgm('CONTROL_TEST','QM/MM '&
               // 'WITH Path MINIMIZATION IS NOT SUPPORTED',& 
               __LINE__,__FILE__)
       ELSE
          ! It should never be possible to go here.
          CALL stopgm('CONTROL','ERROR IN Path OPTIONS',& 
               __LINE__,__FILE__)
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    ! _FM[
    IF (cntl%fmatch.AND..NOT.cntl%tqmmm)THEN
       ! Forcematching only implemented with QM/MM
       CALL stopgm('CONTROL_TEST', 'FORCEMATCHING IMPLEMENTED IN'&
            // 'QMMM',& 
            __LINE__,__FILE__)
    ENDIF
    ! _FM]
    ! ==--------------------------------------------------------------==
    IF (cntl%wfopt)THEN
       ! Backward compatibility. Allow use of MAXSTEP instead of MAXITER
       ! for regular Wavefunction Optimization
       IF (cnti%nomore.NE.10000.AND.cnti%nomore_iter.EQ.10000) THEN
          cnti%nomore_iter=cnti%nomore
       ENDIF
    ENDIF
    ! Kohn-Sham Energies
    IF (cntl%ksener) THEN
       ! By default, we use optimisation of wavefunctions with
       ! diagonalisation scheme (Lanczos algorithm).
       cntl%wfopt=.TRUE.
       cntl%tdiag=.TRUE.
       vdwl%vdwc=.FALSE. ! By default, empirical vdW is OFF
       IF (cntl%cdft)cdftci%cdft_end=1
       IF (.NOT.cntl%tdavi) THEN
          cntl%tlanc=.TRUE.
          ! All eigenvalues are accuratly calculated.
          fint1%tfral=.TRUE.
          IF (fint4%ntabtrot.NE.1) THEN
             fint4%ntabtrot=1
             IF (fint1%ttrot) THEN
                cntr%b2limit=1.e-12_real_8
             ELSE
                cntr%b2limit=1.e-8_real_8
             ENDIF
          ENDIF
          ! No self-consistency (1 iteration), but need an RESTART file.
          cnti%nomore_iter=1
          IF (.NOT.restart1%restart) THEN
             CALL stopgm('CONTROL',&
                  'KOHN-SHAM ENERGIES OPTION NEEDS TO USE RESTART OPTION',& 
                  __LINE__,__FILE__)
          ELSE IF (.NOT.(restart1%rwf.OR.restart1%rpot.OR.restart1%rrho)) THEN
             WRITE(6,'(A,/,A,/,A)')&
                  ' CONTROL! KOHN-SHAM ENERGIES OPTION:',&
                  ' CONTROL! RESTART OPTION NEEDS ONE OF THE ITEMS',&
                  ' CONTROL! WAVEFUNCTIONS, DENSITY OR POTENTIAL'
             CALL stopgm('CONTROL',&
                  'SPECIFY AT LEAST ONE ITEM',& 
                  __LINE__,__FILE__)
          ELSE IF (restart1%rwf.AND.tkpts%tknoswap) THEN
             WRITE(6,'(A,A,/,A,/,A)')&
                  ' CONTROL! KOHN-SHAM ENERGIES ',&
                  'WITH THE OPTION NOWAVEFUNCTION:',&
                  ' CONTROL! WAVEFUNCTIONS ARE NOT USED',&
                  ' CONTROL! YOU NEED A DENSITY OR A POTENTIAL TO RESTART'
             CALL stopgm('CONTROL',&
                  'DO NOT SPECIFY WAVEFUNCTIONS IN RESTART OPTION',& 
                  __LINE__,__FILE__)
          ELSE IF (cntl%cdft.AND..NOT.cdftlog%rcdft) THEN
             WRITE(6,'(A,/,A,/,A)')&
                  'WARNING! KOHN-SHAM ENERGIES',&
                  'Not Restarting from CDFT_RESTART',&
                  'make sure initial Vc is correctly set'
          ENDIF
          IF (tkpts%tknoswap) THEN
             ! The wavefunctions are not stored.
             store1%swf=.FALSE.
             ! The density are not stored.
             store1%srho=.FALSE.
             tsrho=.FALSE.
             ! The potential are not stored.
             store1%spot=.FALSE.
          ENDIF
          ! Set the default for number of Lanczos/Friesner iterations
          IF (cnti%n_fries.EQ.-1) cnti%n_fries=200
       ELSE
          cntl%tlanc=.FALSE.
          ! No self-consistency (1 iteration), but need an RESTART file.
          cnti%nomore_iter=1
          IF (.NOT.restart1%restart .OR. .NOT.restart1%rwf) THEN
             CALL stopgm('CONTROL',&
                  'KOHN-SHAM ENERGIES NEEDS RESTART WAVEUNCTION',& 
                  __LINE__,__FILE__)
          ENDIF
          ! Set the default for number of Davidson iterations
          IF (cnti%n_fries.EQ.-1) THEN
             cnti%n_fries=200
          ENDIF
       ENDIF
       !CSOC[
    ELSEIF (cntl%tspec.OR.cntl%tsoc) THEN
       !CSOC]
       ! Electronic Spectra
       IF (cnti%n_fries.EQ.-1) THEN
          cnti%n_fries=200
       ENDIF
    ELSEIF (cntl%tddft.AND.cntl%tmdbo) THEN
       cnti%n_fries=500
       cntr%b2limit=1.0e-8_real_8
       cnti%nkry_max=8
    ELSEIF (cntl%tdiagopt) THEN
       IF (cnti%n_fries.EQ.-1) THEN
          cnti%n_fries=50
       ENDIF
    ELSE
       IF (cnti%n_fries.EQ.-1) THEN
          cnti%n_fries=1
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (cprint%iprint_step.EQ.0) THEN
       cprint%iprint_step=MAX(cnti%nomore,10000)+1
    ENDIF
    IF (store1%istore.EQ.0) THEN
       store1%istore =MAX(cnti%nomore,10000)+1
    ENDIF
    IF (store1%isctore.EQ.0) THEN
       store1%isctore=cnti%nomore_iter+1
    ENDIF
    IF (cnti%imovie.EQ.0) THEN
       cnti%imovie=cnti%nomore+1
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (cntl%md) THEN
       ! if we have no restart we needed to zero out electron velocities.
       IF (restart1%rnon) cntl%quenche=.TRUE.
       ! 
       IF (cntr%trampr.LT.0.0_real_8)&
            WRITE(6,'(A,A)') ' WARNING! NEGATIVE ',&
            'TEMPERATURE RAMP RATE SPECIFIED. RAMP IGNORED'
       ! 
       IF (lqmmm%qmmm.AND.clc%classical.AND.cntl%tmdbo) cnti%nomore_iter=1

       IF (cntl%tnosee) THEN
          IF (cntl%tc) THEN
             WRITE(6,'(A,A)') ' WARNING! BOTH NOSE ELECTRONS AND ',&
                  'TEMPCONTROL ELECTRONS WERE SPECIFIED'
             WRITE(6,'(A)') ' WARNING! TURNING OFF TEMPCONTROL ELECTRONS'
          ENDIF
          cntl%tc=.FALSE.
          IF (cntl%tbere) THEN
             WRITE(6,'(A,A)') ' WARNING! BOTH NOSE ELECTRONS AND ',&
                  'BERENDSEN ELECTRONS WERE SPECIFIED'
             WRITE(6,'(A)') ' WARNING! TURNING OFF BERENDSEN ELECTRONS'
          ENDIF
          cntl%tbere=.FALSE.
          IF (cntl%tmdbo) THEN
             WRITE(6,'(A,A)') ' WARNING! NOSE ELECTRONS WAS ',&
                  'SPECIFIED WITH BO DYNAMICS'
             WRITE(6,'(A)') ' WARNING! TURNING OFF NOSE ELECTRONS'
             cntl%tnosee=.FALSE.
          ENDIF
       ENDIF
       IF (cntl%tnosep) THEN
          IF (cntl%tcp) THEN
             WRITE(6,'(A,A)')' WARNING! BOTH NOSE IONS AND ',&
                  'TEMPCONTROL IONS WERE SPECIFIED'
             WRITE(6,'(A)')' WARNING! TURNING OFF TEMPCONTROL IONS'
          ENDIF
          cntl%tcp=.FALSE.
          IF (cntl%tberp) THEN
             WRITE(6,'(A,A)') ' WARNING! BOTH NOSE IONS AND ',&
                  'BERENDSEN IONS WERE SPECIFIED'
             WRITE(6,'(A)') ' WARNING! TURNING OFF BERENDSEN IONS'
          ENDIF
          cntl%tberp=.FALSE.
       ENDIF
       IF (cntl%tnosec) THEN
          IF (cntl%tcc) THEN
             WRITE(6,'(A,A)')' WARNING! BOTH NOSE CELL AND ',&
                  'TEMPCONTROL CELL WERE SPECIFIED'
             WRITE(6,'(A)')' WARNING! TURNING OFF TEMPCONTROL CELL'
          ENDIF
          cntl%tcc=.FALSE.
          IF (cntl%tberc) THEN
             WRITE(6,'(A,A)') ' WARNING! BOTH NOSE CELL AND ',&
                  'BERENDSEN CELL WERE SPECIFIED'
             WRITE(6,'(A)') ' WARNING! TURNING OFF BERENDSEN CELL'
          ENDIF
          cntl%tberc=.FALSE.
       ENDIF
       IF (cntl%tdiag.AND.(.NOT.cntl%tmdbo)) THEN
          WRITE(6,'(A)')' WARNING! DIAGONALIZATION REQUESTED.'
          WRITE(6,'(A)')' WARNING! SWITCHING TO BORN-OPPENHEIMER MD'
          cntl%tmdbo=.TRUE.
       ENDIF
       IF (wannl%twmol.AND.(.NOT.cntl%tmdbo))&
            CALL stopgm('CONTROL_TEST','MOLECULAR ORBITALS NOT '&
            // 'NOT COMPATIBLE WITH CP DYNAMICS',& 
            __LINE__,__FILE__)
       IF (cntl%tshock.AND.(cntl%tberc.OR.cntl%tberp.OR.cntl%tbere))&
            CALL stopgm('CONTROL_TEST','BERENDSEN THERMOSTAT NOT '&
            // 'NOT COMPATIBLE WITH SHOCK.',& 
            __LINE__,__FILE__)
       IF (cntl%tshock.AND.(cntl%annec.OR.cntl%annei.OR.cntl%annee))&
            CALL stopgm('CONTROL_TEST','ANNEALING NOT '&
            // 'NOT COMPATIBLE WITH SHOCK.',& 
            __LINE__,__FILE__)
    ELSE
       cntl%tshop=.FALSE.
    ENDIF
    IF (cntl%tnosee.OR.cntl%tnosep.OR.cntl%tnosec) THEN
       ! ===============================================================
       ! T.D. in NOSEPA we use only DELT_ELEC
       ! A.K. ... but with BO-dynamics it doesn't matter
       ! and the wavefunction optimizers use DT2BYE 
       ! to tune the algorithm and the typical values
       ! for BO-dynamics may be too large.
       IF ((.NOT.cntl%tmdbo).AND.(cntr%delt_elec.NE.cntr%delt_ions)) THEN
          WRITE(6,*)&
               'THE TIMESTEPS FOR IONS AND ELECTRONS HAVE TO BE EQUAL'
          WRITE(6,*)&
               'WHEN USING NOSE-HOOVER THERMOSTATS'
          CALL stopgm('CONTROL','TIMESTEPS ARE NOT EQUAL',& 
               __LINE__,__FILE__)
       ENDIF
    ENDIF

    IF (cntl%tbere.AND.(cntr%taube.LE.0.0_real_8)) CALL stopgm('CONTROL',&
         'BERENDSEN CHARACTERISTIC TIME <=0.0 FOR ELECTRONS',& 
         __LINE__,__FILE__)
    IF (cntl%tberp.AND.(cntr%taubp.LE.0.0_real_8)) CALL stopgm('CONTROL',&
         'BERENDSEN CHARACTERISTIC TIME <=0.0 FOR IONS',& 
         __LINE__,__FILE__)
    IF (cntl%tberc.AND.(cntr%taubc.LE.0.0_real_8)) CALL stopgm('CONTROL',&
         'BERENDSEN CHARACTERISTIC TIME <=0.0 FOR CELL',& 
         __LINE__,__FILE__)
    IF (.NOT.(cntl%tmdbo.OR.cntl%tmdfile.OR.cnti%iftype.EQ.3)) THEN
       IF (cntl%textrap)&
            WRITE(6,*) 'WARNING! WAVEFUNCTION EXTRAPOLATION '&
            // 'ONLY SUPPORTED FOR BO-MD. DISABLING...'
       cntl%textrap=.FALSE.
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (cntl%wfopt .OR. cntl%geopt .OR. cntl%vibrat .OR. cntl%tinter .OR. cntl%tmdbo .OR.&
         cntl%fmatch .OR. cntl%tmdfile .OR. cntl%thard .OR. cntl%tspec .OR. cntl%tpspec .OR.&
         cntl%tdebfor.OR.cntl%tsoc) THEN
       IF (cntl%diis .AND. .NOT. cntl%tfrho_upw) THEN
          cnti%mdiis=MIN(cnti%mdiis,maxdis)
          cntl%pcg=.FALSE.
          cntl%tsde=.FALSE.
          IF (cnti%nreset.GE.0.AND.cnti%nreset.LT.cnti%mdiis) THEN
             cnti%nreset=cnti%mdiis
          ENDIF
       ELSEIF (cntl%pcg) THEN
          cntl%tsde=.FALSE.
       ELSEIF (cntl%tsde) THEN
          ! Nothing to be done.
       ELSE
          cntl%diis=.TRUE.
          cntl%prec=.TRUE.
          cnti%mdiis=MIN(cnti%mdiis,maxdis)
          IF (cnti%nreset.GE.0.AND.cnti%nreset.LT.cnti%mdiis) THEN
             cnti%nreset=cnti%mdiis
          ENDIF
       ENDIF
    ENDIF
    ! change the default value of emass (unless already specified differently)
    IF (cntl%cmplx_wf.AND.(cntr%emass.EQ.400)) cntr%emass=1._real_8
    ! ==--------------------------------------------------------------==
    !CSOC[
    IF (cntl%tspec.OR.cntl%tsoc) THEN
       cntl%tdiag=.FALSE.
       IF (.NOT.cntl%tdavi .AND. .NOT.cntl%tlanc) THEN
          cntl%tdavi=.TRUE.
       ENDIF
    ENDIF
    !CSOC]
    ! ==--------------------------------------------------------------==
    ! Diagonalisation scheme
    IF (cntl%tdiag.AND.cntl%wfopt)THEN
       ! MAXSTEP is also used for MAXITER
       IF (cnti%nomore.NE.10000.AND.cnti%nomore_iter.EQ.10000) THEN
          cnti%nomore_iter=cnti%nomore
       ENDIF
    ENDIF
    IF (cntl%tdiag .AND. .NOT.cntl%tfrho_upw) THEN
       cntl%diis=.FALSE.
       cntl%pcg=.FALSE.
       cntl%tsde=.FALSE.
       ! Only Gram-Schmidt works
       cntl%tlowd=.FALSE.
    ELSEIF (cntl%tfrho_upw) THEN
       IF (cntl%diis) THEN
          cntl%pcg=.FALSE.
          cntl%tsde=.FALSE.
       ELSE IF (cntl%pcg) THEN
          cntl%tsde=.FALSE.
       ENDIF
       cnti%mdiis=MIN(cnti%mdiis_fr,maxdis)
       IF (cnti%maxldiis .EQ. 0) cnti%maxldiis = 20
       IF (cnti%minldiis .EQ. 0) cnti%minldiis = 4
       IF (cnti%nreset.GE.0.AND.cnti%nreset.LT.cnti%mdiis) cnti%nreset=cnti%mdiis
       IF (cntr%tolrhofix .EQ.0.0_real_8) cntr%tolrhofix=1.e-7_real_8
       IF (fint1%betael .EQ. 0.0_real_8) cntl%tfint = .FALSE.
    ELSE
       fint1%betael=0._real_8
    ENDIF
    IF (.NOT.cntl%tdiag) THEN
       ! No density mixing
       broy1%tgmix=.FALSE.
       ! No Broyden density mixing
       broy1%tgbroy=.FALSE.
    ENDIF
    ! ==--------------------------------------------------------------==
    ! GEometry OPTimisation
    IF (cntl%geopt) THEN
       IF (cntl%gdiis) THEN
          cntl%rfo=.FALSE.
          cntl%bfgs=.FALSE.
          cntl%lbfgs=.FALSE.
          cntl%tsdp=.FALSE.
          cntl%prfo=.FALSE.
          cnti%mgdiis=MIN(cnti%mgdiis,mxgdis)
          cntl%tinr=.FALSE.
       ELSEIF (cntl%rfo) THEN
          cntl%bfgs=.FALSE.
          cntl%lbfgs=.FALSE.
          cntl%tsdp=.FALSE.
          cntl%tinr=.FALSE.
          cntl%prfo=.FALSE.
       ELSEIF (cntl%bfgs) THEN
          cntl%lbfgs=.FALSE.
          cntl%tsdp=.FALSE.
          cntl%prfo=.FALSE.
       ELSEIF (cntl%prfo) THEN
          cntl%tsdp=.FALSE.
          cntl%lbfgs=.FALSE.
          IF (nsmaxp.EQ.-1) nsmaxp=cnti%nomore
       ELSEIF (cntl%lbfgs) THEN
          cntl%tsdp=.FALSE.
          cntl%tinr=.FALSE.
       ELSEIF (cntl%tsdp) THEN
          cntl%tinr=.FALSE.
       ELSEIF (cntl%tinr) THEN
          cntl%gdiis=.FALSE.
          cntl%tdiag=.FALSE.
          cntl%simul=.FALSE.
       ELSE
          ! AK: QM/MM defaults to cntl%lbfgs for performance reasons
          IF (cntl%tqmmm) THEN
             cntl%lbfgs=.TRUE.
          ELSE
             IF (cntl%tdiag) THEN
                cntl%gdiis=.TRUE.
             ELSE
                cntl%gdiis=.TRUE.
             ENDIF
          ENDIF
       ENDIF
       IF (cntl%tsdc) cntl%tprcp=.TRUE.
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (cntl%tprcp) THEN
       IF (.NOT.qspl1%qspline) THEN
          WRITE(6,*) 'NOTE: VARIABLE CELL ENFORCES "SPLINE QFUNCTION"'
          qspl1%qspline=.TRUE.
       ENDIF
       IF (comvl%tsubrot) THEN
          WRITE(6,*) 'NOTE: NO "SUBTRACT ROTVEL" WITH VARIABLE CELL'
          comvl%tsubrot=.FALSE.
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (group%nogrp.LT.1.OR.group%nogrp.GT.parai%cp_nproc) CALL stopgm('control_test',&
         'The number of groups should be > 0 and < nproc+1 ',& 
         __LINE__,__FILE__)
    IF (cntl%tfint) THEN
       cntl%tdavi=.FALSE.
       IF (fint1%betael.LT.0._real_8) fint1%betael=1000._real_8
       IF (.NOT.fint1%ttrot) fint1%tbogo=.FALSE.
    ELSE
       fint1%betael=0._real_8
    ENDIF
    IF (cntl%tdavi.AND.cntl%tlanc) cntl%tlanc=.FALSE.
    IF (cntl%tlanc) THEN
       IF (cntr%b2limit.LE.0._real_8) THEN
          IF (fint1%ttrot) THEN
             cntr%b2limit=1.e-12_real_8
          ELSE
             cntr%b2limit=1.e-8_real_8
          ENDIF
       ENDIF
    ENDIF
    IF (cntr%tolog.EQ.0._real_8) THEN
       IF (cntl%tdiag) THEN
          cntr%tolog=1.e-3_real_8
       ELSE
          cntr%tolog=1.e-5_real_8
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Store density.
    IF (tsrho) THEN
       store1%srho=.TRUE.
    ELSE IF (.NOT.cntl%tdiag) THEN
       store1%srho=.FALSE.
    ENDIF
    ! Check for distributed linear algebra
    !    IF (cntl%tdmal.AND.(parai%cp_nproc.LE.1)) THEN
    !       WRITE(6,*) ' YOU HAVE A SINGLE TASK! !!',&
    !            ' TURNING OFF DISTRIBUTED LINEAR ALGEBRA'
    !       cntl%tdmal=.FALSE.
    !       cnti%nstblk=1
    !    ENDIF
    ! ==--------------------------------------------------------------==
    ! EXACT FACTORIZATION
    IF (tshl%txfmqc) THEN
       tmw=.TRUE.
       mwi%nwalk=xfmqc%n_xftraj
       IF (mwi%nwalk > parai%cp_nproc) THEN
          WRITE(6,*)'WARNING: NWALK must be at most CP_NPROC, ' &
               // 'reset NWALK = CP_NPROC'
          mwi%nwalk = parai%cp_nproc
          xfmqc%n_xftraj=mwi%nwalk
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    ! CB Some (probably insufficient) checks for BS
    IF (cntl%bsymm) THEN
       WRITE(6,*) ' ! WARNING: PARALLEL BS CODE IS UNDER CONSTRUCTION!'
       IF (.NOT.(cntl%wfopt.OR.(cntl%md.AND.cntl%tmdbo))) THEN
          WRITE(6,*) ' CURRENTLY ONLY WF OPTIMIZATION',&
               ' AND BO MOLECULAR DYNAMICS',    &
               ' ARE IMPLEMENTED FOR BROKEN SYMMETRY.'
          CALL stopgm('CONTROL','BROKEN SYMMETRY & UNSUPPORTED',& 
               __LINE__,__FILE__)
       ENDIF
       IF (cntl%tpath.OR.lspin2%tlse.OR.tkpts%tkpnt.OR.cntl%tdiag) THEN
          WRITE(6,*) ' PI, K-POINTS, LSE, AND KS DIAGONALIZATION',&
               ' ARE NOT IMPLEMENTED FOR BROKEN SYMMETRY.'
          CALL stopgm('CONTROL','BROKEN SYMMETRY & UNSUPPORTED',& 
               __LINE__,__FILE__)
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    ! check that the restart path is valid 
    exist = dir_exists( TRIM(ADJUSTL(fo_info%fpath)) )
    IF (.NOT.exist) CALL stopgm(procedureN,' THE FILEPATH '//&
         TRIM(ADJUSTL(fo_info%fpath))//' IS NOT VALID',& 
         __LINE__,__FILE__)

    ! ==--------------------------------------------------------------==
    ! check that cpmd was compiled with the right cpp if needed

#if defined(_HASNT_BF_STREAM_IO)
    ! reset the variable if we didnt compile with the streaming option
    IF (cntl%is_in_stream) THEN
       WRITE(6,*) 'CONTROL! cannot use the stream'&
            //' IO because CPMD was compiled with the'&
            //' -D__HASNT_BF_STREAM_IO'
       WRITE(6,*) 'CONTROL! reset USE_IN_STREAM to FALSE'
    ENDIF
    cntl%is_in_stream=.FALSE.
#endif

    ! ==--------------------------------------------------------------==
    ! check that cpmd was compiled with the right cpp if needed

#if defined(_HASNT_BF_STREAM_IO)
    ! reset the variable if we didnt compile with the streaming option
    IF (cntl%is_out_stream) THEN
       WRITE(6,*) 'CONTROL! cannot use the stream'&
            //' IO because CPMD was compiled with the'&
            //' -D__HASNT_BF_STREAM_IO'
       WRITE(6,*) 'CONTROL! reset USE_OUT_STREAM to FALSE'
    ENDIF
    cntl%is_out_stream=.FALSE.
#endif
    ! ==--------------------------------------------------------------==
    ! check that cpmd was compiled with the right cpp if needed

#if ! defined(_HAS_CUDA)
    ! reset the variable if we didnt compile with the streaming option
    IF (cp_cuda_env%use_fft) THEN
       WRITE(6,*) 'CONTROL| cannot use cufft'&
            //' because CPMD was not compiled with the'&
            //' -D__HAS_CUDA'
       WRITE(6,*) 'CONTROL| reset USE_FFT to FALSE'
    ENDIF
    cp_cuda_env%use_fft=.FALSE.
#endif
    ! ==--------------------------------------------------------------==
    ! check that cpmd was compiled with the right cpp if needed

#if ! defined(_HAS_CUDA)
    ! reset the variable if we didnt compile with the streaming option
    IF (cp_cuda_env%use_blas) THEN
       WRITE(6,*) 'CONTROL| cannot use cublas'&
            //' because CPMD was not compiled with the'&
            //' -D__HAS_CUDA'
       WRITE(6,*) 'CONTROL| reset USE_BLAS to FALSE'
    ENDIF
    cp_cuda_env%use_blas=.FALSE.
#endif
    ! ==--------------------------------------------------------------==
    ! check that the number of streams and devices per task are valid

    CALL cntr_test_cuda_env ( cp_cuda_env%use_fft, cp_cuda_env%fft_n_devices_per_task, &
         & cp_cuda_env%fft_n_streams_per_device, 'FFT' )

    CALL cntr_test_cuda_env ( cp_cuda_env%use_blas, cp_cuda_env%blas_n_devices_per_task, &
         & cp_cuda_env%blas_n_streams_per_device, 'BLAS' )

    ! ==--------------------------------------------------------------==

  END SUBROUTINE control_test
  ! ==================================================================


  SUBROUTINE cntr_test_cuda_env ( using, n_devices_per_task, n_streams_per_device, name_str )
    LOGICAL, INTENT(INOUT) :: using
    INTEGER, INTENT(INOUT) :: n_devices_per_task, n_streams_per_device
    CHARACTER(*), INTENT(IN) :: name_str

#if defined(_HAS_CUDA)

    IF( using .AND. n_devices_per_task <= 0 ) THEN
       WRITE(6,*) 'CONTROL| USE_'//name_str//' is set to TRUE, but the '//name_str//'_N_DEVICES_PER_TASK <= 0'
       WRITE(6,*) 'CONTROL| reset USE_'//name_str//' to FALSE'
       using = .FALSE.
    ENDIF

    IF( using .AND. n_streams_per_device <= 0 ) THEN
       WRITE(6,*) 'CONTROL| USE_'//name_str//' is set to TRUE, but the '//name_str//'_N_STREAMS_PER_DEVICE <= 0'
       WRITE(6,*) 'CONTROL| reset USE_'//name_str//' to FALSE'
       using = .FALSE.
    ENDIF


    IF( using ) THEN
       !$ IF( .TRUE. ) THEN
       !$    IF( n_streams_per_device > omp_get_max_threads() ) THEN
       !$      WRITE(6,*) 'CONTROL| '//name_str//'_N_STREAMS_PER_DEVICE ', n_streams_per_device,&
       !$           & ' is greater than the avalaible number of threads',&
       !$           & omp_get_max_threads()
       !$      WRITE(6,*) 'CONTROL| '//name_str//'_N_STREAMS_PER_DEVICE is reset to the number of threads'
       !$      n_streams_per_device = omp_get_max_threads()
       !$    ENDIF
       !$ ELSE
       IF( n_streams_per_device > 1 ) THEN
          WRITE(6,*) 'CONTROL| '//name_str//'_N_STREAMS_PER_DEVICE ', n_streams_per_device,&
               & ' is greater than one and CPMD was not compiled with -omp'
          WRITE(6,*) 'CONTROL| reset '//name_str//'_N_STREAMS_PER_DEVICE to one'
          n_streams_per_device = 1
       ENDIF
       !$ ENDIF
    ELSE
       IF( n_streams_per_device > 1 ) THEN
          WRITE(6,*) 'CONTROL| '//name_str//'_N_STREAMS_PER_DEVICE ', n_streams_per_device,&
               & ' is greater than one and USE_'//name_str//'=F'
          WRITE(6,*) 'CONTROL| reset '//name_str//'_N_STREAMS_PER_DEVICE to one'
          n_streams_per_device = 1
       ENDIF
       IF ( n_devices_per_task > 0 ) THEN
          WRITE(6,*) 'CONTROL| '//name_str//'_N_DEVICES_PER_TASK ', n_devices_per_task,&
               & ' is greater than zero and USE_'//name_str//'=F'
          WRITE(6,*) 'CONTROL| reset '//name_str//'_N_DEVICES_PER_TASK to zero'
          n_devices_per_task = 0
       ENDIF
    ENDIF
#else
    IF ( n_devices_per_task > 0 ) THEN
       WRITE(6,*) 'CONTROL| '//name_str//'_N_DEVICES_PER_TASK ', n_devices_per_task,&
            & ' is greater than zero and CPMD was not compiled with',&
            & ' -D__HAS_CUDA'
       WRITE(6,*) 'CONTROL| reset '//name_str//'_N_DEVICES_PER_TASK to zero'
       n_devices_per_task = 0
    ENDIF

    IF ( n_streams_per_device > 1 ) THEN
       WRITE(6,*) 'CONTROL| '//name_str//'_N_STREAMS_PER_DEVICE ', n_streams_per_device,&
            & ' is greater than one and CPMD was not compiled with',&
            & ' -D__HAS_CUDA'
       WRITE(6,*) 'CONTROL| reset '//name_str//'_N_STREAMS_PER_DEVICE to one'
       n_streams_per_device = 1
    ENDIF
#endif

  END SUBROUTINE cntr_test_cuda_env


  LOGICAL FUNCTION dir_exists( Path )
    IMPLICIT NONE
    CHARACTER(*) :: Path

    INTEGER, PARAMETER :: ntrial = 10, unit_beg=10, unit_end=99
    INTEGER, PARAMETER :: nrand=20
    INTEGER :: i, irand, unit, ir, iostat
    CHARACTER(128) :: file, id
    CHARACTER(128+LEN(Path)) :: path_file
    CHARACTER(len=*),PARAMETER :: procedureN='DIR_EXISTS'
    LOGICAL :: exist, opened
    REAL(real_8) :: r

    dir_exists = .FALSE.
    ! find a openable unit
    DO unit = unit_beg, unit_end
       INQUIRE(unit=unit,opened=opened,exist=exist)
       IF (.NOT.opened.AND.exist) EXIT
    ENDDO
    IF (unit>unit_end) CALL stopgm(procedureN,' Cannot find a '//&
         'valid unit',& 
         __LINE__,__FILE__)
    DO i = 1,ntrial
       file='fake_file_'
       DO irand = 1, nrand
          CALL RANDOM_NUMBER(r)
          ir = ANINT ( 10._real_8  * r - 0.5_real_8 )
          WRITE(id,'(I1)') ir
          file=TRIM(file)//TRIM(id)
       ENDDO
       path_file = TRIM(Path)//TRIM(file)
       INQUIRE(file=path_file,exist=exist)
       IF (.NOT.exist) EXIT
    ENDDO
    IF (i>ntrial) CALL stopgm(procedureN,' Cannot find a '//&
         'valid file name',& 
         __LINE__,__FILE__)
    OPEN(unit=unit,file=path_file,iostat=iostat)
    dir_exists = iostat == 0
    CLOSE(unit=unit,status='delete',iostat=iostat)
  END FUNCTION dir_exists

END MODULE control_test_utils
