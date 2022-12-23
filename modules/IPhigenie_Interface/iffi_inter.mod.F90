#if 1
MODULE iffi_inter

 USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt, tilimex, tipri, tistart, ttimp
  USE machine, ONLY: m_signal, m_datum, m_walltime, m_cputime
  USE mp_interface, ONLY : mp_end
  USE prng_utils, ONLY : prnginit, repprngu, repprngu_vec, repprngu_vec_cmplx
  USE control_utils, ONLY : control
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
  USE pm_wf_utils, ONLY : pm_wf
  USE pm_gmopts_utils, ONLY : pm_gmopts
  USE pm_mdpt_utils, ONLY : pm_mdpt
  USE prpt_utils, ONLY : prpt
  USE mdpt_utils, ONLY : mdpt
  USE gmopts_utils, ONLY : gmopts
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
  
  USE iffi_types
  USE iffi_elstat_utils, ONLY: CALC_IMPORT,CALC_EXPORT
  USE iffi_comm_utils,   ONLY: IFFI_READ, IFFI_WRITE
  USE iffi_mgmt_utils,   ONLY: ITERATE_WFCT
  USE ions,           ONLY: ions0,ions1
IMPLICIT NONE

PRIVATE

PUBLIC :: IFFI_CPMD_START,IFFI_CPMD_SCF,IFFI_CPMD_END



CONTAINS
!....Interface to the MD simulation program IPHIGENIE
!.....initialization and finalization of DFT calculation
! called by IPHIGENIE:
!..... IFFI_CPMD_START          : manages CPMD startup (called by IPHIGENIE)
!..... IFFI_CPMD_END            : manages CPMD finalize (called by IPHIGENIE)
!..... IFFI_CPMD_SCF            : data exchange and KS iteration (called by IPHIGENIE)

! initialization of iffi/cpmd:
!..... IFFI_INTERPT_INIT        : initialization (INTERPT): first part of INTERPT in interpt_utils.F90
!..... IFFI_INTERFACE_INIT      : initialization : first part of INTERFACE in egointer_utils.F90
!......IFFI_INIT                : initialization of iffi-exclusive objects

! deinitialization of iffi/cpmd:
!..... IFFI_FINALIZE            : deinitialization of iffi-exclusive objects
!..... IFFI_INTERFACE_FINALIZE  : deinitialization : second part of INTERFACE in egointer_utils.F90
!..... IFFI_INTERPT_FINALIZE    : deinitialization (INTERPT): second part of INTERPT in interpt_utils.F90

!***************************************************************************
                

! this is copy/pasted from the routine CPMD in cpmd.F90
! called by IPHIGENIE
! ==================================================================
SUBROUTINE IFFI_CPMD_START(COMM,CSTR)
! ==--------------------------------------------------------------==
 
  USE mp_interface,                    ONLY : mp_comm_world
  USE system,                          ONLY : cnti

  
  IMPLICIT NONE
  CHARACTER(*), PARAMETER                    :: procedureN = 'IFFI_CPMD_START'

  CHARACTER(len=26)                          :: datx
  INTEGER                                    :: isub
  LOGICAL                                    :: tinfo
  REAL(real_8)                               :: tcpu, time1, &!, time2, twck, &
                                                wclk1, wclk2
 !IPHIGENIEIPHIGENIEIPHIGENIEIPHIGENIEIPHIGENIEIPHIGENIE start
      INTEGER IS,IA
      INTEGER NTHREADS
!$    INTEGER   OMP_GET_MAX_THREADS
!$    EXTERNAL  OMP_GET_MAX_THREADS

!     input filename as argument
      INTEGER COMM
      CHARACTER*(*) CSTR

!.....omp check
      NTHREADS=1
!$    NTHREADS=OMP_GET_MAX_THREADS()
      if (NTHREADS.GT.mxno_threads) then
        CALL stopgm(procedureN,'mxno_threads TOO SMALL',0,"")
      endif

!.....set interface type to IPHIGENIE
      cnti%iftype=3
!.....arguments
      cnts%inputfile=CSTR
      mp_comm_world=COMM
      
 !IPHIGENIEIPHIGENIEIPHIGENIEIPHIGENIEIPHIGENIEIPHIGENIE end
 





  !CALL m_signal(24,SOFTEX)
  !CALL m_signal(30,SOFTEX)
  !CALL m_signal(1,SOFTEX)      ! SIGHUP(1), SIGUSR1 (10) and SIGUSR(12) can also
  !CALL m_signal(10,SOFTEX)     ! be caught. VERY useful for queuing systems...
  !CALL m_signal(12,SOFTEX)
  CALL tistart(time1,wclk1)
  !CALL init_fileopen        ! initialize here to avoid Bus Error
  !CALL get_input_name(cnts%inputfile)
  CALL startpa
      
      
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

  ! READ INFORMATION ON XC FUNCTIONAL
  CALL dftin

  ! READ UNIT CELL DEFINITIONS
  CALL sysin

  ! SET UP SUPERCELL
  CALL setsc

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

     ! GENERATE EXCHANGE AND CORRELATION TABULATION VS. DENSITY
     CALL genxc
     ! CALCULATE NUMBER OF PW
     CALL numpw

     IF (cntl%tpath) THEN
        CALL stopgm(procedureN,"Path integral simulations not supported in IFFI/CPMD mode",0,"")
     !   IF (cntl%tpimd) THEN
     !      ! PATH INTEGRAL INPUT
     !      CALL pi_cntl
     !      CALL pi_init
     !   ELSEIF (cntl%tpmin) THEN
     !      ! PATH MINIMISATION INPUT
     !      CALL pm_cntl
     !      CALL pm_init
     !   ENDIF
     ! EXACT FACTORIZATION
     !IF (tshl%txfmqc) THEN
     !   IF (cntl%tddft) THEN
     !      CALL lr_in
     !      CALL tddft_input
     !      cntl%tsymrho=.FALSE.
     !   ENDIF
     ENDIF
     ! MULTIPLE WALKER METADYNAMICS INITIALIZATION
     IF (tmw) CALL stopgm(procedureN,"MULTIPLE WALKER METADYNAMICS not supported in IFFI/CPMD mode",0,"") !CALL mw_init
     IF (cntl%tresponse) THEN
        print *,"WARNING: Linear response not tested in IFFI/CPMD mode. Good luck!"
        CALL nmr_para_p
     ENDIF
     ! INITIALIZE G-VECTORS AND POINTERS
     CALL rinit
     ! FORM FACTORS
     CALL rinforce
     ! INITIALIZE FFT PARAMETERS
     CALL fft_init ( )
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
  ENDIF
  wclk2 =m_walltime()
  tcpu  = (wclk2 - wclk1)*1.0e-3_real_8
  IF (paral%parent) THEN
     IF (paral%io_parent)&
          WRITE(6,'(/,A,T46,F12.2,A8,/)')&
          ' INITIALIZATION TIME:',tcpu,' SECONDS'
  ENDIF



!IPHIGENIEIPHIGENIEIPHIGENIEIPHIGENIEIPHIGENIEIPHIGENIE start
      CALL IFFI_INTERPT_INIT

      CALL IFFI_INTERFACE_INIT !(C0,C2,SC0,PME,GDE,VPP,EIGV)


!.....get total number of qm-atoms
      runinfo%nrqmatoms=0
      DO IS=1,ions1%nsp
!.......for each species IS loop over NA qm-atoms 
        DO IA=1,ions0%na(IS)
          runinfo%nrqmatoms=runinfo%nrqmatoms+1
!         print *,IS,IA,NA(IS),ZV(IS),runinfo%nrqmatoms
        ENDDO
      ENDDO  
      IF (runinfo%nrqmatoms.GT.maxind) THEN
        CALL stopgm(procedureN,"Too many QM atoms! INCREASE maxind iiffi.inc",0,"")
      ENDIF

      CALL IFFI_READ

      CALL IFFI_INIT
      
!.....now perform interface tasks once to have rgyrs
      CALL CALC_IMPORT
      CALL CALC_EXPORT
      CALL IFFI_WRITE


      IF(paral%io_parent) WRITE(6,'(A)') 'INTERFACE| CPMD INIT COMPLETED'
!IPHIGENIEIPHIGENIEIPHIGENIEIPHIGENIEIPHIGENIEIPHIGENIE end 

END SUBROUTINE
!==================================================================



! this is copy/pasted from routine INTERPT in interpt_utils.mod.F90 until INTERFACE is called
!==================================================================
SUBROUTINE IFFI_INTERPT_INIT
!==--------------------------------------------------------------==
  USE ddip,                            ONLY: lenbk
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE fnlalloc_utils,                  ONLY: fnlalloc,&
                                             fnldealloc
  !USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: parai,&
                                             paral
  USE prmem_utils,                     ONLY: prmem
  USE pslo,                            ONLY: pslo_com
  USE rwfopt_utils,                    ONLY: rwfopt
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             ncpw,&
                                             nkpt
  USE utils,                           ONLY: nxxfun
  USE iffi_wfct

  IMPLICIT NONE
  
  ! Variables
    CHARACTER(*), PARAMETER                  :: procedureN = 'interpt'

    INTEGER                                  :: ierr, nc, nc2, nkpntsc0, &
                                                nsc0, nxx
    LOGICAL                                  :: wasdiis
     
     
! ==--------------------------------------------------------------==
! ==  MEMORY ALLOCATION                                           ==
! ==--------------------------------------------------------------==
! ==  C0    : WAVEFUNCTIONS                                       ==
! ==  C2    : GRADIENTS                                           ==
! ==  SC0   : S**(-1/2)*C0 (nonorthogonal orbitals)               ==
! ==  PME   : cntl%diis WF             HNM1 (for CG)                   ==
! ==  GDE   : cntl%diis GRADIENTS                                      ==
! ==  VPP   : WF HESSIAN                                          ==
! ==  EIGV  : ORBITAL ENERGIES                                    ==
! ==--------------------------------------------------------------==

    IF ( cntl%tqmmm ) CALL stopgm("INTERPT","QMMM NOT IMPLEMENTED",& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    nc=crge%n
    nc2=nc
    nsc0=1
    nkpntsc0=1
    IF (pslo_com%tivan.OR.cntl%nonort.OR.(cntl%pcg.AND.cntl%pcgmin)) THEN
       nsc0=crge%n
       nkpntsc0=nkpt%nkpnt
    ENDIF
    IF (cntl%tdipd) THEN
       IF (cntl%tddft) CALL stopgm("INTERPT",&
            "cntl%tddft AND DIPOLE DYNAMICS NOT POSSIBLE",& 
            __LINE__,__FILE__)
       lenbk=nxxfun(crge%n)
       nxx=MAX(lenbk*parai%nproc / nkpt%ngwk ,nc2)
       nsc0=MAX(nxx,nsc0)
       nc2=nxx
    ENDIF

    ALLOCATE(c0(nkpt%ngwk,nc,nkpt%nkpnt),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(c2(nkpt%ngwk,nc2),STAT=ierr)
    !  allocate(c2(nc2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(eigv(crge%n),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(sc0(nkpt%ngwk,nsc0,nkpntsc0),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    IF (cntl%tsde) THEN
       ! TODO check this DUMMY 

       ALLOCATE(vpp(ncpw%ngw),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ELSE IF (cntl%diis) THEN
       ALLOCATE(pme((ncpw%ngw*crge%n+8)*cnti%mdiis/2),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(gde(((ncpw%ngw*crge%n+8)*cnti%mdiis)/2),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(vpp(ncpw%ngw),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ELSE IF (cntl%pcg) THEN
       ALLOCATE(pme(2*ncpw%ngw*crge%n),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)

       ! TODO check this DUMMY 

       ALLOCATE(vpp(ncpw%ngw),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ..FNL and DFNL
    CALL fnlalloc(crge%n,.TRUE.,.FALSE.)
    IF (paral%parent) CALL prmem('   INTERPT')
    ! ==--------------------------------------------------------------==
    IF (cntl%diis.AND.cntl%tpcgfi) THEN
       IF (paral%io_parent)&
            WRITE(6,*) 'INTERPT| SWITCHING TO PCG MINIMIZE'
       wasdiis=.TRUE.
       cntl%diis=.FALSE.
       cntl%pcg=.TRUE.
       cntl%pcgmin=.TRUE.
    ELSE
       wasdiis=.FALSE.
    ENDIF
    CALL rwfopt(c0,c2,sc0,pme,gde,vpp,eigv)
    IF (wasdiis) THEN
       IF (paral%io_parent)&
            WRITE(6,*) 'INTERPT| SWITCHING BACK TO DIIS'
       cntl%diis=.TRUE.
       cntl%pcg=.FALSE.
       cntl%pcgmin=.FALSE.
    ENDIF
    IF (paral%io_parent)&
         CLOSE(5)
    !

END SUBROUTINE
!==================================================================



!this is copy/pasted from routine INTERFACE in egointer_utils.mod.F90 until the md loop starts
!==================================================================
SUBROUTINE IFFI_INTERFACE_INIT
!==--------------------------------------------------------------==

    USE parac,                           ONLY: paral!,&
                                               !parai
    USE coor,                            ONLY: fion,&
                                               taup!,&
                                               !tau0,&
                                               !velp
    USE efld,                            ONLY: extf!,&
                                               !textfld
    USE system,                          ONLY: &
        fpar, maxsys,cnti !, ngw, nhg, nkpt, nnr1!, &
        !parm, spar, parap, cntl, cntr  
    USE store_types,                     ONLY: restart1
    USE ropt,                            ONLY: iteropt
    USE ener,                            ONLY: ener_com   
    USE forcep_utils,                    ONLY: rhoe_psi_size
    USE machine,                         ONLY: m_walltime
    USE zeroing_utils,                   ONLY: zeroing
    USE iffi_types
    USE iffi_wfct
    USE dynit_utils,                     ONLY: dynit
    USE iffi_mgmt_utils,         ONLY: iffi_give_scr_interface
    USE mm_extrap,                       ONLY: cold
    USE system,                          ONLY: nkpt
    USE elct,                            ONLY: crge
 
    IMPLICIT NONE
    
    CHARACTER(*), PARAMETER                  :: procedureN = 'IFFI_INTERFACE_INIT'

    !CHARACTER(len=80)                        :: int_filen, line
    INTEGER                                  ::  ierr, il_psi_1d, &
                                                il_psi_2d, il_rhoe_1d, &
                                                il_rhoe_2d!, irec(100)!i,
    LOGICAL                                  :: eofrun !, fnowf
    REAL(real_8)                             :: ekin1, &
                                                ekin2, ekincp, ekinh1, &
                                                ekinh2,  temp1, &
                                                temp2, time1!, time2 tcpu, !detot, dummy(1), etoto, 
     CHARACTER(len=30)                        :: tag
     INTEGER                                  :: lscr
    
    
    
    CALL rhoe_psi_size(il_rhoe_1d=il_rhoe_1d, il_rhoe_2d=il_rhoe_2d,&
         il_psi_1d=il_psi_1d,il_psi_2d=il_psi_2d)

    ALLOCATE(rhoe(il_rhoe_1d, il_rhoe_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    ALLOCATE(psi(il_psi_1d,il_psi_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent) CLOSE(5)
    time1 =m_walltime()
    ! ..+1 for esp charges

    ALLOCATE(taup(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    ALLOCATE(fion(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    !!ALLOCATE(velp(3,maxsys%nax,maxsys%nsx),STAT=ierr) (doesnt work either in egointer)
    !!IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
    !!     __LINE__,__FILE__)
    CALL zeroing(taup)!,SIZE(taup))
    CALL zeroing(fion)!,SIZE(fion))

    ALLOCATE(extf(fpar%kr1*fpar%kr2s*fpar%kr3s),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

         
    CALL iffi_give_scr_interface(lscr,tag)
    ALLOCATE(scr(lscr),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

         
         
    CALL zeroing(extf)!,kr1*kr2s*kr3s)
    ! ..addition ends
    CALL zeroing(taup)!,3*maxsys%nax*maxsys%nsx)
    !CALL zeroing(velp)!,3*maxsys%nax*maxsys%nsx)
    CALL zeroing(runinfo%boxqm)!,6)
    CALL zeroing(runinfo%boxdum)!,6)
    CALL zeroing(epot1_fix%koloc)!,3*maxind)
    CALL zeroing(epot1_var%locexp)!,10*maxind)
    CALL zeroing(epot1_list%locid)!,maxind)
    CALL zeroing(epot1_list%nincl)!,maxind)
    CALL zeroing(epot1_list%idnear)!,maxnear)
    CALL zeroing(epot1_list%lincl)!,maxind*maxnear)   
    
    !DO i=1,100
    !   epot2%myrag(i)=1.2_real_8
    !ENDDO
    !epot2%myrag(1)=0.85_real_8
    !epot2%myrag(8)=1.23_real_8

    epot1_list%nnear = 0
    epot1_list%nloc  = 0
    
    !FIXME: zeroing woanders

    restart1%restart=.FALSE.
    restart1%rwf    =.TRUE.
    eofrun=.FALSE.
    iteropt%nfi   = 0
    ener_com%ecnstr= 0.0_real_8
    ener_com%erestr= 0.0_real_8
!X    etoto = 0.0_real_8


!   allocating memory for wavefunction extrapolation
    IF (cntl%textrap) THEN
       ALLOCATE(cold(nkpt%ngwk,crge%n,nkpt%nkpnt,cnti%mextra*nkpt%nkpts/nkpt%nkpnt),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF

      
    CALL dynit(ekincp,ekin1,ekin2,temp1,temp2,ekinh1,ekinh2)

    IF (paral%io_parent) WRITE(6,'(A)') ' INTERFACE| ENTERING IPHIGENIE-INTERFACE'
    
 END SUBROUTINE
!==================================================================
 
 


!==================================================================
SUBROUTINE IFFI_INIT
!==--------------------------------------------------------------==
      !USE efld,                            ONLY: extf
      USE system,                          ONLY: fpar
      USE iffi_types
      
      USE zeroing_utils,                   ONLY: zeroing
      
      IMPLICIT NONE      
      CHARACTER(len=*), PARAMETER              :: procedureN = "IFFI_INIT"
      INTEGER ::  ierr
      
               
      ALLOCATE(EXTFMM(fpar%KR1*fpar%KR2S*fpar%KR3S),STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)                
       
      ALLOCATE(VOXEXP(dimlocexp,MYNVOX),STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
           __LINE__,__FILE__)      
      
      ALLOCATE(VOXEXP_MM(dimlocexp,MYNVOX),STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
           __LINE__,__FILE__)

      ALLOCATE(NEARLIST(MYNVOX,maxnearpervo),STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
           __LINE__,__FILE__)

      ALLOCATE(NEARLISTLEN(MYNVOX),STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
           __LINE__,__FILE__)

      ALLOCATE(OLELIST(MYNVOX,maxolepervo),STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
           __LINE__,__FILE__)
      
      ALLOCATE(OLELISTLEN(MYNVOX),STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
           __LINE__,__FILE__)

      ALLOCATE(QMPOTANDFLD(4,maxind),STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
           __LINE__,__FILE__)

      CALL zeroing(QMPOTANDFLD) !FIXME temporary
           
!!      ALLOCATE(ECHRG(NAT+1),STAT=ierr)
!!      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
!!           __LINE__,__FILE__)
      
      ALLOCATE(MMPOTEXP(dimpotexp,maxnear,mxno_threads),STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
           __LINE__,__FILE__)
      
      ALLOCATE(MMPOTEXP_VAR(DIMPOTEXP_VAR,maxnear,mxno_threads),STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
           __LINE__,__FILE__)

      ALLOCATE(VOMULTEXP(dimmultexp,MYNVOX,mxno_threads),STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
           __LINE__,__FILE__)

      ALLOCATE(ATMULTEXP(dimmultexp,maxind,mxno_threads),STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
           __LINE__,__FILE__)

      ALLOCATE(VOXRGYREXP(dimrgyrexp,MYNVOX,mxno_threads),STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
           __LINE__,__FILE__)

      ALLOCATE(ATRGYREXP(dimrgyrexp,maxind,mxno_threads),STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
           __LINE__,__FILE__)
      

      IF (runinfo%meanfieldmode) THEN
        ALLOCATE(EXTFMEAN(fpar%KR1*fpar%KR2S*fpar%KR3S),STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)

        ALLOCATE(EXTFLAST(fpar%KR1*fpar%KR2S*fpar%KR3S),STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
             __LINE__,__FILE__)
      ENDIF
        

      !IF (paral%parent) THEN
        ALLOCATE(MMPOTEXPGLOB(dimpotexp,maxnear),STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)

        ALLOCATE(MMPOTEXPGLOB_VAR(DIMPOTEXP_VAR,maxnear),STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)

        ALLOCATE(ATMULTEXPGLOB(dimmultexp,maxind),STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)


        ALLOCATE(ATRGYREXPGLOB(dimrgyrexp,maxind),STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
      !ENDIF
      
!==--------------------------------------------------------------==
END SUBROUTINE
!==--------------------------------------------------------------==
   


!==================================================================
SUBROUTINE IFFI_FINALIZE
!==--------------------------------------------------------------==
      USE mm_extrap,                       ONLY: cold
      IMPLICIT NONE                  
      
      CHARACTER(*), PARAMETER                  :: procedureN = 'IFFI_FINALIZE'
      INTEGER :: ierr
      
      DEALLOCATE(EXTFMM,STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
           __LINE__,__FILE__)

      DEALLOCATE(VOXEXP,STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
           __LINE__,__FILE__)

      DEALLOCATE(VOXEXP_MM,STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
           __LINE__,__FILE__)

      DEALLOCATE(NEARLIST,STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
           __LINE__,__FILE__)

      DEALLOCATE(NEARLISTLEN,STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
           __LINE__,__FILE__)

      DEALLOCATE(OLELIST,STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
           __LINE__,__FILE__)

      DEALLOCATE(OLELISTLEN,STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
           __LINE__,__FILE__)

      DEALLOCATE(QMPOTANDFLD,STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
           __LINE__,__FILE__)

      !DEALLOCATE(ECHRG,STAT=ierr)
      !IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
      !     __LINE__,__FILE__)

      DEALLOCATE(MMPOTEXP,STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
           __LINE__,__FILE__)

      DEALLOCATE(MMPOTEXP_VAR,STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
           __LINE__,__FILE__)

      DEALLOCATE(VOMULTEXP,STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
           __LINE__,__FILE__)

      DEALLOCATE(ATMULTEXP,STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
           __LINE__,__FILE__)

      DEALLOCATE(VOXRGYREXP,STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
           __LINE__,__FILE__)

      DEALLOCATE(ATRGYREXP,STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
           __LINE__,__FILE__)

      IF (runinfo%meanfieldmode) THEN
        DEALLOCATE(EXTFMEAN,STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)

        DEALLOCATE(EXTFLAST,STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
      ENDIF
      
      !IF (paral%parent) THEN
        DEALLOCATE(MMPOTEXPGLOB,STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)

        DEALLOCATE(MMPOTEXPGLOB_VAR,STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)

        DEALLOCATE(ATMULTEXPGLOB,STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)

        DEALLOCATE(ATRGYREXPGLOB,STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
      !ENDIF
                            
      
      DEALLOCATE(VOXPART,STAT=ierr) !kr1 not kr1s!
      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
           __LINE__,__FILE__)
           
           
      !this was allocated in INIT_VOXEL     
      DEALLOCATE(KOVOX,STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
           __LINE__,__FILE__)
           
      DEALLOCATE(VO2IQ,STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
           __LINE__,__FILE__)
           
      IF (cntl%textrap) THEN     
        DEALLOCATE(cold,STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
             __LINE__,__FILE__)
      ENDIF
           
END SUBROUTINE   
! ==================================================================




!this is copy/pasted from INTERFACE in egointer_utils.mod.F90 from where the md loop ends 
!=================================================================
SUBROUTINE IFFI_INTERFACE_FINALIZE
!==--------------------------------------------------------------==

    USE parac,                           ONLY: paral!,&
                                               !parai
    USE coor,                            ONLY: fion,&
                                               tau0,&
                                               taup,&
                                               velp
    USE efld,                            ONLY: extf!,&
                                               !textfld
    !USE system,                          ONLY: &
     !   cnti, cntl, cntr, kr1, kr2s, kr3s, maxsys, ngw, nhg, nkpt, nnr1!, &
        !parm, spar, parap
    !USE store_types,                     ONLY: restart1
    !USE ropt,                            ONLY: iteropt
    !USE ener,                            ONLY: ener_com   
    USE forcep_utils,                    ONLY: rhoe_psi_size
    USE machine,                         ONLY: m_walltime
    USE iffi_types
    USE iffi_wfct
    USE finalp_utils,                    ONLY: finalp
          
    
    IMPLICIT NONE                  
    CHARACTER(*), PARAMETER                  :: procedureN = 'IFFI_INTERFACE_FINALIZE'
    REAL(real_8)                             :: dummy(1)
    INTEGER                                  :: ierr
    
    !CALL write_irec(irec)
    !CALL zhwwf(2,irec,c0,c2,n,dummy,tau0,velp,taup,iteropt%nfi)
    IF (paral%parent) CALL finalp(tau0,fion,velp,dummy)
    DEALLOCATE(taup,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(fion,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    !DEALLOCATE(velp,STAT=ierr)
    !IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
    !     __LINE__,__FILE__)
    ! ..added: biswas.
    DEALLOCATE(extf,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(psi,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rhoe,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    
END SUBROUTINE
! ==================================================================



!==================================================================
!this is copy/pasted from INTERPT in interpt_utils.mod.F90 from the line on where INTERFACE was called
!==================================================================
SUBROUTINE IFFI_INTERPT_FINALIZE
!==--------------------------------------------------------------==

  USE error_handling,                  ONLY: stopgm
  USE fnlalloc_utils,                  ONLY: fnlalloc,&
                                             fnldealloc
  !USE parac,                           ONLY: parai,&
  !                                           paral
  USE pslo,                            ONLY: pslo_com                                           
  USE system, ONLY : cntl
      
  USE iffi_wfct                                             
        IMPLICIT NONE
        
    CHARACTER(*), PARAMETER                  :: procedureN = 'interpt'

    INTEGER                                  :: ierr!, nc, nc2, nkpntsc0, &
                                                !nsc0, nxx
    !LOGICAL                                  :: wasdiis
        
    DEALLOCATE(c0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(c2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(eigv,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    CALL fnldealloc(.TRUE.,.FALSE.)
    IF (cntl%nonort.OR.pslo_com%tivan.OR.(cntl%pcg.AND.cntl%pcgmin)) THEN
       DEALLOCATE(sc0,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    IF (cntl%diis) THEN
       DEALLOCATE(pme,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(gde,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    DEALLOCATE(vpp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
         
    DEALLOCATE(scr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
!     ==--------------------------------------------------------------==
      RETURN
END SUBROUTINE
!     ==================================================================


! called by IPHIGENIE second part of CPMD in cpmd.F90
! ==================================================================           
SUBROUTINE IFFI_CPMD_END
! ==--------------------------------------------------------------==
      USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
      USE parac, ONLY : paral !,parai
      !USE system, ONLY : cntl
      USE machine, ONLY: m_signal, m_datum, m_walltime, m_cputime

      
      IMPLICIT NONE
      
      CHARACTER(len=26)                          :: datx
      INTEGER                                    :: isub
      REAL(real_8)                               :: tcpu, time1, time2, twck, &
                                                wclk1, wclk2
      
      CALL IFFI_FINALIZE
      
      CALL IFFI_INTERFACE_FINALIZE

      CALL IFFI_INTERPT_FINALIZE
      
!from here again copy/paste code from CPMD in cpmd.F90

  ! finalize fft (dealloc, ...)
  CALL fft_finalize ( )

  !CALL tihalt('CPMD', isub)

  IF (paral%parent) THEN
     CALL end_swap
  ENDIF

  CALL tipri()

  IF (paral%parent) THEN
     TIME1 = 0.0D0
     WCLK1 = 0.0D0
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
  !IF (cntl%tqmmm) parai%allgrp=parai%qmmmgrp
  !IF (cntl%tpath.OR.tmw) parai%allgrp=supergroup

  CALL finalize_cp_grp()
! CALL mp_end() done by iffi

!==--------------------------------------------------------------==
END SUBROUTINE
!==================================================================


! called by IPHIGENIE
!==================================================================
SUBROUTINE IFFI_CPMD_SCF
!==--------------------------------------------------------------==
      USE iffi_types
      USE parac, ONLY : paral
  
      IMPLICIT NONE
      
      !get data from iphigenie
      CALL IFFI_READ

      IF ((.NOT.runinfo%meanfieldmode).OR.epot1_var%updatemeanfield) THEN        
        !calculate external potential on grid
        CALL CALC_IMPORT      

        !iterate KS orbitals
        CALL ITERATE_WFCT 
      ELSE IF (paral%io_parent) THEN
          print *,"SKIP EXTF AND ITERATION"
      ENDIF

      !export electrostatics from rhoe
      CALL CALC_EXPORT

      !put data to iffi
      CALL IFFI_WRITE

!     flush screen output
      !CALL MY_FLUSH(6)

END SUBROUTINE
!     ==================================================================


END MODULE iffi_inter
#endif
