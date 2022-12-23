! ==================================================================
PROGRAM cpmd_stuttgart

#ifdef __HAS_IEEE_EXCEPTIONS
  USE, INTRINSIC :: ieee_exceptions, ONLY: ieee_set_halting_mode
  USE, INTRINSIC :: ieee_exceptions, ONLY: ieee_usual
#endif

  IMPLICIT NONE

#ifdef __HAS_IEEE_EXCEPTIONS
  CALL ieee_set_halting_mode(ieee_usual,.TRUE.)
#endif

  CALL cpmd

  STOP  ! for F_trace option
END PROGRAM cpmd_stuttgart
! ==================================================================
! ==================================================================
SUBROUTINE cpmd
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt, tilimex, tipri, tistart, ttimp
  USE machine, ONLY: m_signal, m_datum, m_walltime, m_cputime
  USE mp_interface, ONLY : mp_end
  USE prng_utils, ONLY : prnginit, repprngu, repprngu_vec, repprngu_vec_cmplx
  USE control_utils, ONLY : control,get_input_name
  USE dftin_utils, ONLY : dftin
  USE sysin_utils, ONLY : sysin
  USE setsc_utils, ONLY : setsc
  USE detsp_utils, ONLY : detsp
  USE mm_init_utils, ONLY : mm_init
  USE read_prop_utils, ONLY : read_prop
  USE ratom_utils, ONLY : ratom
  USE vdwin_utils, ONLY : vdwin
  USE propin_utils, ONLY : propin
  USE respin_p_utils, ONLY : respin_p
  USE setsys_utils, ONLY : setsys
  USE genxc_utils, ONLY : genxc
  USE numpw_utils, ONLY : numpw
  USE pi_cntl_utils, ONLY : pi_cntl
  USE pi_init_utils, ONLY : pi_init
  USE meta_multiple_walkers_utils, ONLY : mw_init
  USE bicanonicalCpmd, ONLY :  bicanonicalCpmdConfig,&
    bicanonicalCpmdInputConfig, New, Delete,& 
    biCanonicalEnsembleDo
  USE bicanonicalConfig, ONLY : New, Delete
  USE nmr_para_p_utils, ONLY : nmr_para_p
  USE rinit_utils, ONLY : rinit
  USE rinforce_utils, ONLY : rinforce
  USE fftprp_utils, ONLY : fft_init, fft_finalize
  USE initclust_utils, ONLY : initclust
  USE dginit_utils, ONLY : dg_init
  USE nosalloc_utils, ONLY : nosalloc
  USE exterp_utils, ONLY : exterp
  USE setbasis_utils, ONLY : setbasis
  USE dqgalloc_utils, ONLY : dqgalloc
  USE gle_utils, ONLY : gle_alloc
  USE pi_wf_utils, ONLY : pi_wf
  USE pi_mdpt_utils, ONLY : pi_mdpt
  USE pi_prpt_utils, ONLY : pi_prpt
  USE pm_wf_utils, ONLY : pm_wf
  USE pm_gmopts_utils, ONLY : pm_gmopts
  USE pm_mdpt_utils, ONLY : pm_mdpt
  USE prpt_utils, ONLY : prpt
  USE mdpt_utils, ONLY : mdpt
  USE gmopts_utils, ONLY : gmopts
  USE vdw_wf_alloc_utils, ONLY : vdw_wf_alloc
  USE wfopts_utils, ONLY : wfopts
  USE interpt_utils, ONLY : interpt
  USE secdpt_utils, ONLY : secdpt
  USE proppt_utils, ONLY : proppt
  USE orbhard_utils, ONLY : orbhard
  USE response_p_utils, ONLY : response_p
  USE prep_forcematch_utils, ONLY : prep_forcematch
  USE prmem_utils, ONLY : prmem
  USE system, ONLY : cntl,cnts
  USE parac, ONLY : paral,parai
  USE pimd, ONLY : supergroup
  USE mw, ONLY : tmw
  USE specpt_utils, ONLY: specpt
  USE ropt, ONLY : init_pinf_pointers
  USE soc, ONLY : soc_calc
  USE set_cp_grp_utils, ONLY : finalize_cp_grp
  USE pm_init_utils, ONLY : pm_init
  USE startpa_utils, ONLY : startpa
  USE envir_utils, ONLY : envir
  USE pm_cntl_utils, ONLY : pm_cntl
  USE setcnst_utils, ONLY : setcnst
  USE softex_utils, ONLY : softex
  USE fileopen_utils, ONLY : init_fileopen
  USE linres, ONLY : tshl
  USE lr_in_utils, ONLY: lr_in,&
       tddft_input
  USE cp_cuda_utils, ONLY: cp_cuda_init, cp_cuda_finalize
  USE ortho_utils, ONLY: ortho_init, ortho_finalize
  USE mts_utils, ONLY: read_mts_input

  IMPLICIT NONE
  CHARACTER(*), PARAMETER                    :: procedureN = 'cpmd'

  CHARACTER(len=26)                          :: datx
  INTEGER                                    :: isub
  LOGICAL                                    :: tinfo
  REAL(real_8)                               :: tcpu, time1, time2, twck, &
                                                wclk1, wclk2

  CALL m_signal(24,SOFTEX)
  CALL m_signal(30,SOFTEX)
  CALL m_signal(1,SOFTEX)      ! SIGHUP(1), SIGUSR1 (10) and SIGUSR(12) can also
  CALL m_signal(10,SOFTEX)     ! be caught. VERY useful for queuing systems...
  CALL m_signal(12,SOFTEX)
  CALL tistart(time1,wclk1)
  CALL init_fileopen        ! initialize here to avoid Bus Error
  CALL get_input_name(cnts%inputfile)
  CALL startpa

  call New(bicanonicalCpmdConfig)
  call New(bicanonicalCpmdInputConfig)

  tinfo=.TRUE.

  ! infi and infw are now parts of iteropts_t but
  ! but iteropts%infi cannot be DO-loop iterator (we have plenty...)
  CALL init_pinf_pointers()

  CALL m_datum (datx)
  IF (paral%io_parent) WRITE(6,'(A,A)') ' PROGRAM CPMD STARTED AT: ',datx
  ! ==--------------------------------------------------------------==
#ifndef __ES
  CALL init_fileopen        ! initialize here outside ES
#endif

  CALL envir
  CALL setcnst
  ! READ CONTROL PARAMETERS

  CALL control


  ! we need to start the timer after the call to control
  CALL tiset(procedureN, isub)

  ! READ MTS INFORMATION
  CALL read_mts_input

  ! READ INFORMATION ON XC FUNCTIONAL
  CALL dftin

  ! READ UNIT CELL DEFINITIONS
  CALL sysin

  ! SET UP SUPERCELL
  CALL setsc

  ! init cuda/cublas/cufft
  CALL cp_cuda_init ( )


  ! READ ATOMS, PSEUDOPOTENTIALS, COORDINATES
  CALL detsp                ! detect number of species.
  CALL mm_init              ! set up some QM/MM stuff
  ! also needed for non-qmmm runs, so we can 
  ! reduce use of #if defined(__GROMOS)/#endif 

  ! EHR[
  IF (cntl%cmplx_wf) THEN
     CALL read_prop
  ENDIF
  ! EHR]
  IF ( paral%qmnode ) THEN
     ! READ ATOMIC INPUT
     CALL ratom
     ! READ VDW PARAMETERS
     CALL vdwin
     ! READ PROPERTIES
     CALL propin(tinfo)

     ! linear response or implicit newton raphson (docs?).
     IF (cntl%tresponse.OR.cntl%tinr) CALL respin_p

     ! SET UP SYSTEM
     CALL setsys

    call New(bicanonicalCpmdConfig, bicanonicalCpmdInputConfig)

    ! GENERATE EXCHANGE AND CORRELATION TABULATION VS. DENSITY
     CALL genxc
     ! CALCULATE NUMBER OF PW
     CALL numpw

     IF (cntl%tpath) THEN
        IF (cntl%tpimd) THEN
           ! PATH INTEGRAL INPUT
           CALL pi_cntl
           CALL pi_init
        ELSEIF (cntl%tpmin) THEN
           ! PATH MINIMISATION INPUT
           CALL pm_cntl
           CALL pm_init
        ENDIF
     ENDIF
     ! EXACT FACTORIZATION
     IF (tshl%txfmqc) THEN
        IF (cntl%tddft) THEN
           CALL lr_in
           CALL tddft_input
           cntl%tsymrho=.FALSE.
        ENDIF
     ENDIF
     ! MULTIPLE WALKER METADYNAMICS INITIALIZATION
     IF (tmw)CALL mw_init
     IF (cntl%tresponse) THEN
        CALL nmr_para_p
     ENDIF
     ! INITIALIZE G-VECTORS AND POINTERS
     CALL rinit
     ! FORM FACTORS
     CALL rinforce
     ! INITIALIZE FFT PARAMETERS
     CALL fft_init ( )
     ! INITIALIZE ORTHO
     CALL ortho_init ( )
     ! CLUSTER BOUNDARY CONDITIONS 
     CALL initclust
     ! INITIALIZE FFT DATA FOR THE DOUBLEGRID
     CALL dg_init
     ! NOSE ARRAYS
     CALL nosalloc
     ! EXTERNAL POTENTIAL
     CALL exterp
     ! SET ATOMIC BASIS (READ BASIS SECTION)
     CALL setbasis
     ! ALLOCATE DQG ARRAY
     CALL dqgalloc
     ! PRNG
     CALL prnginit
     ! GLE
     CALL gle_alloc
     ! vdW-WF
     CALL vdw_wf_alloc
  ENDIF
  wclk2 =m_walltime()
  tcpu  = (wclk2 - wclk1)*1.0e-3_real_8
  IF (paral%parent) THEN
     IF (paral%io_parent)&
          WRITE(6,'(/,A,T46,F12.2,A8,/)')&
          ' INITIALIZATION TIME:',tcpu,' SECONDS'
  ENDIF

  IF (cntl%tpath) THEN
     IF (cntl%tpimd) THEN
        ! ==----------------------------------------------------------==
        ! ==  PATH INTEGRAL (REPLICA)                                 ==
        ! ==----------------------------------------------------------==
        IF (cntl%wfopt) THEN
           ! ==--------------------------------------------------------==
           ! ==  SINGLE POINT WAVEFUNCTION OPTIMIZATION                ==
           ! ==--------------------------------------------------------==
           CALL pi_wf
        ELSEIF (cntl%md) THEN
           ! ==--------------------------------------------------------==
           ! ==  MOLECULAR DYNAMICS                                    ==
           ! ==--------------------------------------------------------==
           IF (cntl%tprcp) THEN
              CALL pi_prpt
           ELSE
              CALL pi_mdpt
           ENDIF
        ENDIF
     ELSEIF (cntl%tpmin) THEN
        ! ==----------------------------------------------------------==
        ! ==  PATH MINIMISATION (REPLICA)                             ==
        ! ==----------------------------------------------------------==
        IF (cntl%wfopt) THEN
           ! ==--------------------------------------------------------==
           ! ==  SINGLE POINT WAVEFUNCTION OPTIMIZATION                ==
           ! ==--------------------------------------------------------==
           CALL pm_wf
        ELSEIF (cntl%geopt) THEN
           ! ==--------------------------------------------------------==
           ! ==  GEOMETRY OPTIMIZATION                                 ==
           ! ==--------------------------------------------------------==
           CALL pm_gmopts
        ELSEIF (cntl%md) THEN
           ! ==--------------------------------------------------------==
           ! ==  MOLECULAR DYNAMICS                                    ==
           ! ==--------------------------------------------------------==
           CALL pm_mdpt
        ENDIF
     ENDIF
  ELSE
     ! ==------------------------------------------------------------==
     ! ==  NO REPLICA                                                ==
     ! ==------------------------------------------------------------==
     IF (cntl%md) THEN
        ! ==----------------------------------------------------------==
        ! ==  MOLECULAR DYNAMICS                                      ==
        ! ==----------------------------------------------------------==
        IF (cntl%tprcp) THEN
           CALL prpt
        ELSE
           CALL mdpt
        ENDIF
     ELSEIF (cntl%geopt.OR.cntl%tdebfor) THEN
        ! ==----------------------------------------------------------==
        ! ==  GEOMETRY OPTIMIZATION OR FORCE DEBUGGING                ==
        ! ==----------------------------------------------------------==
        CALL gmopts
        ! EHR[
     ELSEIF (cntl%tpspec.OR.cntl%tpdist) THEN
        ! ==----------------------------------------------------------==
        ! ==  ELECTRONIC SPECTRA COMPUTED WITH TD-KS EQUATIONS        ==
        ! ==  PROPAGATION OF A DENSITY PERTURBATION                   ==
        ! ==----------------------------------------------------------==
        CALL wfopts
        ! EHR]
     ELSEIF (cntl%wfopt) THEN
        ! ==----------------------------------------------------------==
        ! ==  SINGLE POINT WAVEFUNCTION OPTIMIZATION OR               ==
        ! ==  DIAGONALIZATION OF KOHN-SHAM MATRIX                     ==
        ! ==----------------------------------------------------------==
        CALL wfopts
     ELSEIF (cntl%tinter) THEN
        ! ==----------------------------------------------------------==
        ! ==  INTERFACE TO CLASSICAL cntl%md CODE                          ==
        ! ==----------------------------------------------------------==
        CALL interpt
     ELSEIF (cntl%vibrat) THEN
        ! ==----------------------------------------------------------==
        ! ==  VIBRATIONAL ANALYSIS                                    ==
        ! ==----------------------------------------------------------==
        CALL secdpt
     ELSEIF (cntl%proper) THEN
        ! ==----------------------------------------------------------==
        ! ==  CALCULATE SOME PROPERTIES                               ==
        ! ==----------------------------------------------------------==
        CALL proppt
     ELSEIF (cntl%thard) THEN
        ! ==----------------------------------------------------------==
        ! ==  ORBITAL HARDNESS                                        ==
        ! ==----------------------------------------------------------==
        CALL orbhard
     ELSEIF (cntl%tspec) THEN
        ! ==----------------------------------------------------------==
        ! ==  ELECTRONIC SPECTRA (cntl%tddft)                              ==
        ! ==----------------------------------------------------------==
        CALL specpt
        !CSOC[
     ELSEIF(cntl%tsoc) THEN
        ! ==----------------------------------------------------------==
        ! ==  SPIN-ORBIT COUPLING CALCULATION (TDDFT)                 ==
        ! ==----------------------------------------------------------==
        CALL soc_calc
        !CSOC]
     ELSEIF (cntl%tresponse) THEN
        ! ==----------------------------------------------------------==
        ! ==  CALCULATE SOME OTHER PROPERTIES                         ==
        ! ==----------------------------------------------------------==
        CALL response_p
        ! ==----------------------------------------------------------==
        ! ==  TOPOLOGY UPDATE VIA FORCE MATCHING                      ==
        ! ==----------------------------------------------------------==
     ELSEIF (cntl%fmatch) THEN
        CALL prep_forcematch
     ENDIF
  ENDIF
  ! ==--------------------------------------------------------------==

  ! finalize fft (dealloc, ...)
  CALL fft_finalize ( )
  ! finalize ortho
  CALL ortho_finalize ( )

  ! finalize cuda/cublas/cufft
  CALL cp_cuda_finalize ( )

  CALL tihalt(procedureN, isub)

  IF (paral%parent) THEN
     CALL end_swap
  ENDIF

  CALL tipri()

  IF (paral%parent) THEN
     time2 = m_cputime()
     tcpu  = (time2 - time1)
     wclk2 =m_walltime()
     twck  = (wclk2 - wclk1)*1.0e-3_real_8
     CALL ttimp(tcpu,twck)
     CALL m_datum (datx)
     CALL prmem('      CPMD')
     IF (paral%io_parent)&
          WRITE(6,'(/,A,A)') ' PROGRAM CPMD ENDED AT:   ',datx
  ENDIF
  IF (paral%parent) THEN
     ! Print some messages if CPU limit time exceeded of Soft Exit.
     CALL tilimex
     ! IF(EXSOFT) THEN
     ! CALL SOFTEX
     ! ENDIF
  ENDIF

  ! Restore original group settings for the final checks.
  IF (cntl%tqmmm) parai%allgrp=parai%qmmmgrp
  IF (cntl%tpath.OR.tmw) parai%allgrp=supergroup

    call Delete(bicanonicalCpmdConfig)
    if (biCanonicalEnsembleDo) call Delete(bicanonicalCpmdInputConfig)
    
  CALL finalize_cp_grp()
  CALL mp_end()
  ! ==--------------------------------------------------------------==
END SUBROUTINE cpmd
! ==================================================================
