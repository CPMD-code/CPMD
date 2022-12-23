module bicanonicalCpmd 
  USE cnst,                            ONLY: fbohr
  use coor,                            only: fion, tau0, velp
  USE clas,                            ONLY: tclas
  USE adat,                            ONLY: covrad, atwt
  USE dum2_utils,                      ONLY: dum2  
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE fillc_utils,                     ONLY: fillc
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_def,&
                                             fo_app,&
                                             fo_verb
  USE filnmod,                         ONLY: filbod,&
                                             filn
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
!  USE kinds,                           ONLY: real_8

  USE mm_dimmod,                       ONLY: cpat, cpsp
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_send,&
                                             mp_recv,&
                                             mp_group,&
                                             mp_split,&
                                             mp_sync,&
                                             mp_mpi_error_assert
  USE parac,                           ONLY: parai,&
                                             paral
  USE pbc_utils,                       ONLY: pbc                                           
  USE pimd,                            ONLY: grandparent,&
                                               ipcurr,&
                                             nproc_tot,&
                                             parentgroup,&
                                             pc_groups,&
                                             pc_grp,&
                                             pcg_pos,&
                                             supergroup,&
                                             supersource
  USE readsr_utils,                    ONLY: xstring
  USE set_cp_grp_utils,                ONLY: reset_cp_grp
  USE soft,                            ONLY: soft_com
  USE system,                          ONLY: cntr, cnts,&
                                             cntl, maxsys,&
                                             parap, parm, maxsp
  USE testex_utils,                    ONLY: testex_mw
  USE timer,                           ONLY: tihalt,&
                                             tiset
  !transfer cpmd input to other subroutines 
  use bicanonicalConfig, only: bicanonicalConfigType,&
    print
  use bicanonicalCalculationConfig, only: GetNumberOfSpeciesInSecondCanonicalSystem
  use bicanonicalCalculation, only: BicanonicalCalculationType,&
    New, Delete, updateMuXweight, GetChemicalPotential, GetChemicalPotentialAtXWeightOneHalf,&
    GetXWeight, Energies, bicanonicalGradients, Print, PrintDebug
  USE wann,                            ONLY: wan05,&
                                             wannl
  USE zeroing_utils,                   ONLY: zeroing

#ifdef __PARALLEL
  USE mpi
#endif

  
  implicit none
  
  private
  save
  public New, Delete
  public PrintEnergies
  
  PUBLIC bicanonicalEnsembleInitParallel
!  PUBLIC BicanonicalSetFileNameLatest
  public BicanonicalTestExit
  public bicanonicalCpmdConfig
  public bicanonicalCpmdInputConfig
  public CpmdEnergiesGradients
  public biCanonicalEnsembleDo
  public GetNameEnergiesTape, GetNameTrajecTape, GetNameTrajectoryTape,&
    GetNameConstraintTape, GetNameFtrajectoryTape, GetNameMetricTape,&
    GetNameNoseEnergyTape, GetnameNoseTrajecTape, GetNameLatestTape,&
    GetNameEnergiesVarAtTape, GetEnergiesVarAtTape, GetNameGeometryTape,&
    GetNameGeometryXyzTape, GetNameGeometryScaleTape
  public SanityChecks
! public  
  type BicanonicalCpmdType
    integer:: energiesVarAtTape
    character(14) :: nameEnergiesVarAtTape 
    integer:: restartTape
    character(20) :: nameRestartTape 
    integer:: latestTape
    character(13) :: nameLatestTape 
    integer:: energiesTape
    character(13) :: nameEnergiesTape 
    integer:: trajecTape
    character(20) :: nameTrajecTape 
    integer:: trajectoryTape
    character(20) :: nameTrajectoryTape 
    integer:: ftrajectoryTape
    character(20) :: nameFtrajectoryTape 
    integer:: metricTape
    character(20) :: nameMetricTape 
    integer:: constraintTape
    character(20) :: nameConstraintTape 
    integer:: noseEnergyTape
    character(20) :: nameNoseEnergyTape 
    integer:: noseTrajecTape
    character(20) :: nameNoseTrajecTape 
    integer:: geometryTape
    character(20) :: nameGeometryTape 
    integer:: geometryXyzTape
    character(20) :: nameGeometryXyzTape 
    integer:: geometryScaleTape
    character(20) :: nameGeometryScaleTape 
    type (bicanonicalConfigType) :: bicanonicalConfig
    type (bicanonicalCalculationType), pointer :: bicanonicalCalculation=>NULL()
 end type BicanonicalCpmdType
  !
  ! JF TODO this is used as global variable for all the calls in cpmd routines 
  type (BicanonicalCpmdType) :: bicanonicalCpmdConfig
  ! JF TODO this is used as global variable for transfer config from control to cpmd routines 
  type (BicanonicalConfigType) :: bicanonicalCpmdInputConfig
  !private  
  
  !JF We need this type to setup the bicanonionical infrastructure early in the cpmd code when
  ! bicanonicalCpmdConfig was not initialized. 
  type biCanonicalEnsembleType 
    integer :: numberOfCanonicalEnsembles = 2
    integer :: idCanonicalEnsemble
 end type biCanonicalEnsembleType
   type(biCanonicalEnsembleType) :: biCanonicalEnsemble
  logical :: biCanonicalEnsembleDo = .false.
  
  !  public Getter, Setter  
  integer, parameter :: VARIABLE_ATOM_NUMBER_ENERGIES_TAPE = 62
  character(len=14), parameter :: NAME_VARIABLE_ATOM_NUMBER_ENERGIES_TAPE = 'ENERGIES'
  integer, parameter :: OUT_TAPE = 6
  character(len=13), parameter :: NAME_OUT_TAPE = 'OUTPUT_CNF'

  interface New 
    module procedure NewPrivate, NewSplitParallelPrivate
 end interface New

  interface Delete
    module procedure DeletePrivate
 end interface Delete

  interface SanityChecks
    module procedure SanityChecksPrivate
 end interface SanityChecks

  interface PrintEnergies
    module procedure PrintEnergiesPrivate
 end interface PrintEnergies

  interface bicanonicalEnsembleInitParallel
    module procedure bicanonicalEnsembleInitParallelPrivate
 end interface bicanonicalEnsembleInitParallel

!  interface BicanonicalSetFileNameLatest
!    module procedure BicanonicalSetFileNameLatestPrivate
 !  end interface

  interface BicanonicalTestExit
    module procedure BicanonicalTestExitPrivate
 end interface BicanonicalTestExit

  interface CpmdEnergiesGradients 
    module procedure CpmdEnergiesGradientsPrivate
 end interface CpmdEnergiesGradients

  interface getEnergiesVarAtTape
    module procedure getEnergiesVarAtTapePrivate
 end interface getEnergiesVarAtTape

  interface getNameEnergiesVarAtTape
    module procedure getNameEnergiesVarAtTapePrivate
 end interface getNameEnergiesVarAtTape

  interface getNameLatestTape
    module procedure getNameLatestTapePrivate
 end interface getNameLatestTape

  interface getNameEnergiesTape
    module procedure getNameEnergiesTapePrivate
 end interface getNameEnergiesTape

  interface getNameTrajecTape
    module procedure getNameTrajecTapePrivate
 end interface getNameTrajecTape

  interface getNameTrajectoryTape
    module procedure getNameTrajectoryTapePrivate
 end interface getNameTrajectoryTape

  interface getNameConstraintTape
    module procedure getNameConstraintTapePrivate
 end interface getNameConstraintTape

  interface getNameFtrajectoryTape
    module procedure getNameFtrajectoryTapePrivate
 end interface getNameFtrajectoryTape

  interface getNameMetricTape
    module procedure getNameMetricTapePrivate
 end interface getNameMetricTape

  interface getNameNoseEnergyTape
    module procedure getNameNoseEnergyTapePrivate
 end interface getNameNoseEnergyTape

  interface getNameNoseTrajecTape
    module procedure getNameNoseTrajecTapePrivate
 end interface getNameNoseTrajecTape

  interface getNameGeometryTape
    module procedure getNameGeometryTapePrivate
 end interface getNameGeometryTape

  interface getNameGeometryXyzTape
    module procedure getNameGeometryXyzTapePrivate
 end interface getNameGeometryXyzTape

  interface getNameGeometryScaleTape
    module procedure getNameGeometryScaleTapePrivate
 end interface getNameGeometryScaleTape

contains

  subroutine NewSplitParallelPrivate(self)
    type (BicanonicalCpmdType), intent(inout) :: self
    character(*), parameter                  :: procedureN = 'NewSplitParallelPrivate'

    ! Check CPMD input if Bicanonical was requested 
    call RunCheckInputPrivate

    if (.not.biCanonicalEnsembleDo) return
    
    ! Split the procs according to the number of chemical systems 
    call bicanonicalEnsembleInitParallelPrivate
  end subroutine NewSplitParallelPrivate
  
  subroutine NewPrivate(self, config)
    type (BicanonicalConfigType), intent(in) :: config
    type (BicanonicalCpmdType), intent(out) :: self
    !integer :: nspBican2, natBican2, msgid, ierr
    character(*), parameter                  :: procedureN = 'NewPrivate'

    if (.not.biCanonicalEnsembleDo) return
    !check if we have CP MD running 
    if (cntl%tPath .or. cntl%tMdBo .or. tClas .or. cntl%tMdFile .or. &
      cntl%tNabdy .or. cntl%tQmMm .or. cntl%tShop .or. cntl%tPmin)&
      call stopgm(procedureN,&
      'BicanonicalModule implemented for CP MD runtype, only.',& 
       __LINE__,__FILE__)
    if  (.not. (cntl%md .and. .not.cntl%tprcp)) call stopgm(procedureN,&
      'BicanonicalModule implemented for CP MD runtype, only.',& 
       __LINE__,__FILE__)
    
    self%bicanonicalConfig = config
    self%bicanonicalConfig%calculationConfig%numberOfCanonicalEnsembles = biCanonicalEnsemble%numberOfCanonicalEnsembles
    self%bicanonicalConfig%calculationConfig%idCanonicalEnsemble = biCanonicalEnsemble%idCanonicalEnsemble

    !JF Sanity checks and transfer of CPMD stuff 
    ! from large canonical system to small one.  
    call SanityChecksPrivate(self)
    if (grandparent)&
      self%bicanonicalConfig%calculationConfig%atomMassOfExcessSpecies = atwt(ions0%iatyp(maxsys%nsx))
    
    ! initialize the output streams 
    call bicanonicalInitIOPrivate(self)
    !JF TODO DO  add to bicanonicalInitIOPrivate
    CALL bicanonicalAssignRestartFilename(self, GetEnsembleIdPrivate(self))
    !JF TODO However if we add the assignment to bicanonicalInitIOPrivate 
    ! it will not transfer the correct RESTART file name to cpmd routines. 
    ! Debugging welcome
    ! filbod = getNameRestartTapePrivate(self)
    
    self%bicanonicalConfig%calculationConfig%volumeSimulationBox = parm%omega
    
    self%bicanonicalConfig%calculationConfig%ioNode = paral%io_parent
    
    CALL mp_bcast(self%bicanonicalConfig%calculationConfig%atomMassOfExcessSpecies,parai%io_source,supergroup)
    CALL mp_bcast(self%bicanonicalConfig%calculationConfig%numberOfSpeciesInSecondCanonicalSystem,parai%io_source,supergroup)
    CALL mp_bcast(self%bicanonicalConfig%calculationConfig%XWeight,parai%io_source,supergroup)
    CALL mp_bcast(self%bicanonicalConfig%calculationConfig%chemicalPotential,parai%io_source,supergroup)
    CALL mp_bcast(self%bicanonicalConfig%calculationConfig%temperatureGammaFactor,parai%io_source,supergroup)
    
    if (paral%io_parent) call Print(self%bicanonicalConfig, OUT_TAPE)

    allocate(self%bicanonicalCalculation)
    call New(self%bicanonicalCalculation, self%bicanonicalConfig%calculationConfig)

    if (paral%io_parent) call Print(self%bicanonicalCalculation, OUT_TAPE)

  end subroutine NewPrivate

  subroutine DeletePrivate(self)
    type (BicanonicalCpmdType), intent(inout) :: self
    parai%allgrp=supergroup
    if (.not.biCanonicalEnsembleDo)return
    call fileClose(self%energiesVarAtTape)
    !call Delete(self%bicanonicalCalculation)
    !deallocate(self%bicanonicalCalculation)
  end subroutine DeletePrivate

  subroutine RunCheckInputPrivate
    CHARACTER(*), PARAMETER                  :: procedureN = 'RunCheckInputPrivate'
    INTEGER, PARAMETER                       :: IUNIT = 5, LINELEN = 80

    CHARACTER(len=linelen)                   :: line
    INTEGER                                  :: iostat
    INTEGER :: tmp = 0

    biCanonicalEnsembleDo = .false.
    ! ==--------------------------------------------------------------==
    ! 
    ! determine if bicanonicalRunType was requested
    ! 
    IF (parai%cp_me.EQ.0) THEN
       OPEN(unit=iunit,file=cnts%inputfile,status='OLD',iostat=iostat)
       IF (iostat.NE.0)CALL stopgm(procedureN,&
            'Cannot open the input file',& 
            __LINE__,__FILE__)
       outer: DO
          READ(iunit,iostat=iostat,fmt='(A80)') line
          IF (iostat.NE.0) EXIT outer
          IF (INDEX(line,'&CPMD').NE.0) THEN
             inner: DO
                READ(iunit,iostat=iostat,fmt='(A80)') line
                IF (iostat.NE.0) EXIT outer
                IF (INDEX(line,'&END').NE.0) EXIT outer
                IF (INDEX(line,'BICAN').NE.0 .AND.&
                  index(line,'ENS').NE.0) THEN
                     biCanonicalEnsembleDo = .true.
                     tmp = 1 
                ENDIF
             ENDDO inner
          ENDIF
       ENDDO outer
       CLOSE(unit=iunit,iostat=iostat)
       IF (iostat.NE.0)CALL stopgm(procedureN,&
            'Cannot close the input file',& 
            __LINE__,__FILE__)
    ENDIF
    !CALL mp_bcast(cntr%ecut,parai%io_source,parai%cp_grp) ! TODO is it needed?
    call mp_sync(parai%cp_grp)
    call mp_bcast(tmp, parai%io_source, parai%cp_grp)

    !JF TODO This is a temporary work around to keep CPMD running
    !        if BicanonicalEnsemble runtype was requested but with
    !        providing only one single MPI thread. 
    !        In this case we perfom a Canonical ensemble calculation
    !        for the big system, only. 
    ! TODO   The number of MD steps will be set   NFI=1.
    !
    !write (OUT_TAPE,*) "parai%cp_nproc", parai%cp_nproc
    if (parai%cp_nproc /= 1 .and. tmp == 1) then 
      biCanonicalEnsembleDo = .true.
    elseif (parai%cp_nproc == 1 .and. tmp == 1) then
      biCanonicalEnsembleDo = .true.
      biCanonicalEnsemble%numberOfCanonicalEnsembles = 1 
      write (OUT_TAPE, *) 
      write (OUT_TAPE, *) "!!! Warning !!!",&
         "BicanonicalEnsemble runtype was requested ",&
         "providing only one mpi thread, ",&
         "but the implementation requires two or more. ",&
         "We proceed, but with the results being ",& 
         "meaningless and are",&
         " _NOT_ those of the BicanonicalEnsemble."
       !cnti%nomore=1
      write (OUT_TAPE, *) 
   end if
  end subroutine RunCheckInputPrivate


  subroutine SanityChecksPrivate(self)
    type (BicanonicalCpmdType), intent(inout) :: self
    integer,allocatable :: ionsIatypBican2(:),ionsIatomBican2(:)
    integer :: nspBican2, natBican2, msgid, ierr
    real(real_8) :: checkVolumeSimulationBox
    real(real_8), allocatable :: tau0Bican2(:,:,:), velpBican2(:,:,:)
    character(*), parameter                  :: procedureN = 'SanityChecksPrivate'
    
    
    
    !JF TODO This is a temporary work around to keep CPMD running
    !        if BicanonicalEnsemble runtype was requested but with
    !        providing only one single MPI thread. 
    !        In this case we perfom a Canonical ensemble calculation
    !        for the big system, only. 
    ! TODO   The number of MD steps will be set   NFI=1.
    !
    if  (biCanonicalEnsemble%numberOfCanonicalEnsembles == 1 ) return
    
    
    IF (paral%io_parent.AND. biCanonicalEnsemble%idCanonicalEnsemble == 2 ) THEN
      msgid = 16
      call mp_send(ions1%nsp, pc_grp(1), msgid, supergroup )
      msgid = 17
      call mp_send(ions1%nat, pc_grp(1), msgid, supergroup )
      msgid = 18
      call mp_send(ions0%iatyp, maxsys%nsx, pc_grp(1), msgid, supergroup )
      msgid = 19
      call mp_send(ions0%na, maxsys%nsx, pc_grp(1), msgid, supergroup )
      msgid = 20
      call mp_send(self%bicanonicalConfig%calculationConfig%volumeSimulationBox, pc_grp(1), msgid, supergroup )
      msgid = 21
      call mp_send(tau0, 3*maxsys%nax*maxsys%nsx, pc_grp(1), msgid, supergroup )
      msgid = 22
      call mp_send(velp, 3*maxsys%nax*maxsys%nsx, pc_grp(1), msgid, supergroup )
    
    elseif (grandparent) then 
      !Check runtype  

      msgid = 16
      call mp_recv(nspBican2, pc_grp(2), msgid, supergroup )
      !if (nspBican2 /= GetNumberOfSpeciesInSecondCanonicalSystem(self%bicanonicalConfig%calculationConfig) ) &
      !  CALL stopgm(procedureN, 'Input error: Wrong number of species in smaller canonical system speciefied',& 
      !      __LINE__,__FILE__)
      if (nspBican2 + 1  /= ions1%nsp) then
        write (OUT_TAPE,*)  "No. species input _2", nspBican2, "input 1", ions1%nsp
        CALL stopgm(procedureN, 'Wrong number of species in canonical system 2. Must be one less compared to system 1.',& 
          __LINE__,__FILE__)
     end if
      self%bicanonicalConfig%calculationConfig%numberOfSpeciesInSecondCanonicalSystem = nspBican2
    
      msgid = 17
      call mp_recv(natBican2, pc_grp(2), msgid, supergroup )
      if (natBican2 + 1 /= ions1%nat)  then 
        write (OUT_TAPE,*)  "No. of atoms input _2",natBican2,"input 1" ,ions1%nat
        CALL stopgm(procedureN, 'Input error: Number of atoms in canonical system 2 must be one less of system 1.',& 
            __LINE__,__FILE__)
      endif

      
      msgid = 18
      allocate(ionsIatypBican2(maxsp),STAT=ierr)
      IF(ierr/=0) call stopgm(procedureN,'allocation problem', __LINE__,__FILE__)
      call zeroing(ionsIatypBican2)
      call mp_recv(ionsIatypBican2, maxsp, pc_grp(2), msgid, supergroup )
      if ( .not.equalPrivate(ionsIatypBican2(1:nspBican2),ions0%iatyp(1:nspBican2) ) )& 
        CALL stopgm(procedureN, 'Input error: atom types of canonical system 2 differ from that of system 1.',& 
            __LINE__,__FILE__)
      deallocate(ionsIatypBican2, STAT=ierr)
      if(ierr/=0) call stopgm(procedureN,'deallocation problem', __LINE__,__FILE__)

      msgid = 19
      allocate(ionsIatomBican2(maxsp),STAT=ierr)
      IF(ierr/=0) call stopgm(procedureN,'allocation problem', __LINE__,__FILE__)
      call zeroing(ionsIatomBican2)
      call mp_recv(ionsIatomBican2, maxsp, pc_grp(2), msgid, supergroup )
      if ( .not.equalPrivate(ionsIatomBican2(1:nspBican2),ions0%na(1:nspBican2) ) ) &
        CALL stopgm(procedureN, 'Input error: number of atoms in species of canonical system 2 differ from that of system 1.',& 
            __LINE__,__FILE__)
      deallocate(ionsIatomBican2, STAT=ierr)
      if(ierr/=0) call stopgm(procedureN,'deallocation problem', __LINE__,__FILE__)
      
      msgid = 20
      call mp_recv(checkVolumeSimulationBox, pc_grp(2), msgid, supergroup )
      if (3.0_real_8*epsilon(0.0_real_8) < abs (checkVolumeSimulationBox - &
          self%bicanonicalConfig%calculationConfig%volumeSimulationBox )) &
        CALL stopgm(procedureN, 'Input error: Volume of canonical systems different',& 
            __LINE__,__FILE__)

      msgid = 21
      allocate(tau0Bican2(3,maxsys%nax,nspBican2),STAT=ierr)
      IF(ierr/=0) call stopgm(procedureN,'allocation problem', __LINE__,__FILE__)
      call zeroing(tau0Bican2)
      call mp_recv(tau0Bican2,3*maxsys%nax*nspBican2, pc_grp(2), msgid, supergroup )
      if ( .not.globalEqualPrivate(tau0Bican2(1:3,1:maxsys%nax,1:nspBican2),tau0(1:3,1:maxsys%nax,1:nspBican2) ) ) &
        CALL stopgm(procedureN, 'Input error: atom positions of canonical system 2 differ from that of system 1.',& 
            __LINE__,__FILE__)
      deallocate(tau0Bican2, STAT=ierr)
      if(ierr/=0) call stopgm(procedureN,'deallocation problem', __LINE__,__FILE__)

      msgid = 22
      allocate(velpBican2(3,maxsys%nax,nspBican2),STAT=ierr)
      IF(ierr/=0) call stopgm(procedureN,'allocation problem', __LINE__,__FILE__)
      call zeroing(velpBican2)
      call mp_recv(velpBican2,3*maxsys%nax*nspBican2, pc_grp(2), msgid, supergroup )
      if ( .not.globalEqualPrivate(velpBican2(1:3,1:maxsys%nax,1:nspBican2),velp(1:3,1:maxsys%nax,1:nspBican2) ) ) &
        CALL stopgm(procedureN, 'Input error: atom velocities in canonical system 2 differ from that of canocical system 1.',& 
            __LINE__,__FILE__)
      deallocate(velpBican2, STAT=ierr)
      if(ierr/=0) call stopgm(procedureN,'deallocation problem', __LINE__,__FILE__)
      
   end if
    
     !vel must be equal  Test after restart

 end subroutine SanityChecksPrivate

  subroutine Bicanonical_fileName(fileNameIn,fileNameOut,ipl)
    ! ==--------------------------------------------------------------==
    ! == Like some multiple walker metadynamics output files have     ==
    ! == names in the format NAME_#id(.EXTENTION)                     ==
    ! ==--------------------------------------------------------------==
    ! 
    character(len=*), intent(in)             :: fileNameIn
    character(len=*), intent(out)            :: fileNameOut
    integer                                  :: ipl

    character(len=50)                        :: cipNum
    integer                                  :: i1, i2, n1, n2

    fileNameOut=' '
    if (paral%io_parent)&
         write(cipNum,'(i4)') ipl
    call xstring(filenamein,n1,n2)
    call xstring(cipNum,i1,i2)
    fileNameOut=fileNameIn(n1:n2)//cipNum(i1:i2)
    return
  end subroutine bicanonical_filename

  subroutine BicanonicalAssignRestartFilename(self, ipcurr)!,filen,filenbs)
    type (BicanonicalCpmdType), intent(inout) :: self
    ! ==--------------------------------------------------------------==
    ! == Create filenames: RESTART_CNF#ip.  RESTART_CNF#ip 
    ! ==--------------------------------------------------------------==
    ! 
    integer                                  :: ipcurr
    character(len=100)                       :: filRest
    integer                                  :: n1, n2

    if (.not.biCanonicalEnsembleDo)return

    filRest = 'RESTART_CNF'
    call bicanonical_filename(filRest,filn,ipcurr)
    call xstring(filn,n1,n2)
!    filbod=filn(n1:n2)//'.'
    filbod = getNameRestartTapePrivate(self)
    
    return
  end subroutine BicanonicalAssignRestartFilename

  integer function GetEnsembleIdPrivate(self)
    type (BicanonicalCpmdType), intent(inout) :: self
    GetEnsembleIdPrivate = self%bicanonicalConfig%calculationConfig%idCanonicalEnsemble
  end function GetEnsembleIdPrivate

  SUBROUTINE bicanonicalEnsembleInitParallelPrivate
    ! ==--------------------------------------------------------------==
    ! ==    INITIALIZE VARIABLE ATOM NUMBER PARALLEL ENVIRONMENT      ==
    ! ==--------------------------------------------------------------==
    CHARACTER(*), PARAMETER                  :: procedureN = &
      'bicanonicalInitPa'

    CHARACTER(len=12)                        :: fileout
    INTEGER                                  :: color, i, ip, ipp, &
                                                isub, metmp, nhigh, nlow, &
                                                pg_me
    LOGICAL                                  :: ferror
    REAL(real_8)                             :: rproc
    
    
    !if (paral%io_parent) write(*,*) 'Entering ', procedureN
    IF (.NOT.biCanonicalEnsembleDo)RETURN
    !sanity check 
    !if (biCanonicalEnsemble%numberOfCanonicalEnsembles /= 2 ) then
    !   CALL stopgm('bicanonicalInitPara',&
    !     '! Number of canonical Systems not equal two -> code broken!',& 
    !        __LINE__,__FILE__)
    !ENDIF
    CALL tiset(procedureN,isub)


    !Parallel setup
    
    CALL mp_bcast(cntr%ecut,parai%io_source,parai%cp_grp) ! TODO is it needed?

    ! ==--------------------------------------------------------------==
    ! Reassigning PARENTS and SOURCE plus define GRANDPARENT and SUPERSOURCE
    grandparent=paral%io_parent
    supergroup=parai%cp_grp
    supersource=parai%io_source
    nproc_tot=parai%cp_nproc
    pg_me=parai%cp_me
    ! ==--------------------------------------------------------------==
    ! Generate processor groups
    CALL zeroing(parap%pgroup)!,maxcpu)
    CALL zeroing(parap%nlink)!,maxcpu)
    ! Always 1 canoncal system per processor group
    pc_groups = biCanonicalEnsemble%numberOfCanonicalEnsembles
    ! 
    rproc=REAL((parai%cp_nproc/biCanonicalEnsemble%numberOfCanonicalEnsembles),kind=real_8)
    DO i=1,biCanonicalEnsemble%numberOfCanonicalEnsembles
       nlow=NINT((i-1)*rproc)
       pc_grp(i)=nlow
       nhigh=NINT(i*rproc)-1
       IF (parai%cp_me.GE.nlow .AND. parai%cp_me.LE.nhigh) THEN
          parai%cp_nproc=nhigh-nlow+1! number of processors in a group is redifined on all corresponding processors
          parai%mepos=parai%cp_me-nlow! mepos is redifined on all processors of a group
          pcg_pos=i
          biCanonicalEnsemble%idCanonicalEnsemble=i
          DO ip=1,parai%cp_nproc
             parap%pgroup(ip)=ip-1
             ipp=parap%pgroup(ip)
             parap%nlink(ipp)=ip-1
          ENDDO
       ENDIF
    ENDDO

    color=parai%cp_me/parai%cp_nproc
    CALL mp_split(supergroup,color,parai%cp_me,parai%cp_grp,metmp)
    IF (parai%mepos.NE.metmp) THEN
       CALL stopgm('bicanonicalInitPara','! MY_SPLIT RETURNED ERROR!',& 
            __LINE__,__FILE__)
    ENDIF
    parai%cp_me=metmp
    CALL mp_sync(supergroup)
    CALL mp_group(biCanonicalEnsemble%numberOfCanonicalEnsembles,pc_grp,parentgroup,supergroup)
    !
    !     Reset the cp group communicator
    !
    CALL reset_cp_grp()

    ! ==--------------------------------------------------------------==
    IF (paral%io_parent.AND. pg_me+1>parai%cp_nproc) THEN
       !WRITE(fileon,'(I3)') pcg_poS
       !CALL xstring(fileon,i1,i2)
       !fileout='OUTPUT_CNF'//fileon(i1:i2)
       call bicanonical_filename(NAME_OUT_TAPE, fileout, biCanonicalEnsemble%idCanonicalEnsemble)
       CALL fileclose(OUT_TAPE)
       CALL fileOpen(OUT_TAPE, fileout,fo_def,ferror)
#if defined(__Linux) || defined (__ALPHALINUX)
       ! avoid printmemsize messes
       CALL silentstdout
#endif
    ENDIF

    IF (paral%io_parent) THEN
!    IF (grandparent) THEN
       WRITE(OUT_TAPE, '(/,1X,64("*"))')
       WRITE(OUT_TAPE, '(1X,A4,13X,A29,14X,A4)') '**  ','BICANONICAL CALCULATION SPLIT', '**'
       WRITE(OUT_TAPE, '(/,1X,64("*"))')
       WRITE(OUT_TAPE, '(1X,A4,A18,A38,A4)') '**  ','CANONICAL SYSTEM', 'PROCESSORS',&
            '**'
       WRITE(OUT_TAPE, '(1X,64("-"))')
       DO i=1,biCanonicalEnsemble%numberOfCanonicalEnsembles
          nlow=pc_grp(i)
          IF (i.EQ.biCanonicalEnsemble%numberOfCanonicalEnsembles) THEN
             nhigh=nproc_tot-1
          ELSE
             nhigh=pc_grp(i+1)-1
          ENDIF
          WRITE(OUT_TAPE, '(1X,A4,I18,I17,A,I17,A4)') '**  ',i,&
               nlow,' -> ',nhigh, '**'
       ENDDO
       WRITE(OUT_TAPE, '(1X,64("*"),/)')
    ENDIF

    IF (wannl%twann) THEN
       IF (wan05%loc_npgrp>parai%cp_nproc) wan05%loc_npgrp=parai%cp_nproc
    ENDIF

    !JF init file name for second system 
    if (biCanonicalEnsemble%idCanonicalEnsemble == biCanonicalEnsemble%numberOfCanonicalEnsembles) then 
      !filn = cnts%inputfile
      !call xstring(filn,n1,n2)
      !cnts%inputfile = filn(n1:n2)//'_2'
      call setInputFileNameSecondCanoncalSystem(cnts%inputfile,cnts%inputfile)
   end if

    call tihalt(procedureN,isub)
    return
  end subroutine bicanonicalEnsembleInitParallelPrivate

  subroutine setInputFileNameSecondCanoncalSystem(filn,inputfile) 
    character(len=*)  :: filn,inputfile
    integer :: n1,n2

    filn = inputfile
    call xstring(filn,n1,n2)
    inputfile = filn(n1:n2)//'_2'

  end subroutine setInputFileNameSecondCanoncalSystem


  subroutine BicanonicalTestExitPrivate(softExit)
    logical, intent(inout) :: softExit

    IF (.NOT.biCanonicalEnsembleDo)RETURN
    CALL testex_mw(soft_com%exsoft)

  end subroutine bicanonicalTestExitPrivate

  subroutine CpmdEnergiesGradientsPrivate(self, etot)
    type (BicanonicalCpmdType), intent(inout) :: self
    real(real_8), intent(in)  :: eTot 
    CHARACTER(*), PARAMETER                  :: procedureN = 'SetEnergiesGradientsPrivate'
    real(real_8), allocatable :: gradients(:,:,:)
    INTEGER :: ierr, nsp2, msgid
    
    if (.not.biCanonicalEnsembleDo) return
    
    
    !JF TODO This is a temporary work around to keep CPMD running
    !        if BicanonicalEnsemble runtype was requested but with
    !        providing only one single MPI thread. 
    !        In this case we perfom a Canonical ensemble calculation
    !        for the big system, only. 
    ! TODO   The number of MD steps will be set   NFI=1.
    !
    !write (OUT_TAPE,*) "parai%cp_nproc", parai%cp_nproc
    if (parai%cp_nproc == 1 ) then
      write (OUT_TAPE, *) 
      write (OUT_TAPE, *) "!!! Warning !!!",&
         "BicanonicalEnsemble runtype was requested, ",&
         "providing only one mpi thread, ",&
         "but the implementation requires two or more. ",&
         "We proceed, but with the results being ",& 
         "meaningless and are",&
         " _NOT_ those of the BicanonicalEnsemble."
       !cnti%nomore=1
      write (OUT_TAPE, *)  
      return
   end if



    IF (paral%io_parent.AND. biCanonicalEnsemble%idCanonicalEnsemble == 2 ) THEN
      self%bicanonicalCalculation%energyCnf2 = eTot
      msgid = 1
      call mp_send(eTot, pc_grp(1), msgid, supergroup )
      msgid = 2
      call mp_recv(self%bicanonicalCalculation%energyCnf1, pc_grp(1), msgid, supergroup )
    elseif (grandparent) then 
      self%bicanonicalCalculation%energyCnf1 = eTot
      msgid = 1
      call mp_recv(self%bicanonicalCalculation%energyCnf2, pc_grp(2), msgid, supergroup )
      msgid = 2
      call mp_send(eTot, pc_grp(2), msgid, supergroup )
    endif 
  
    if (paral%parent) call Energies(self%bicanonicalCalculation)
   
    nsp2 = GetNumberOfSpeciesInSecondCanonicalSystem(self%bicanonicalCalculation%config)
    
    ALLOCATE(gradients(3, maxsys%nax, maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', __LINE__,__FILE__)
    CALL zeroing(gradients)
    
    if (paral%parent) then
      call bicanonicalGradients(self%bicanonicalCalculation, fion)
   end if
 
    ! JF Apply wall potential to avoid close contacts 
    !
    ! Wall potential is acting only on the extra species in the big system.
    ! The xWeight of the small system is applyed on the wall potential gradients
    ! We assume the big system to be placed on  pg which incudes the grandparent 
    !
    if (grandparent) then
      call ApplyWallPotential(self, gradients)
      gradients(:,:,:) = (1.0_real_8 &
        - GetXWeight(self%bicanonicalCalculation))*gradients(:,:,:)
      fion(:,:,:) = fion(:,:,:) + gradients(:,:,:)
   end if
    
    IF (paral%parent.AND. biCanonicalEnsemble%idCanonicalEnsemble == 2 ) THEN
      msgid = 11
      call mp_send(fion(:,:,1:nsp2),3*maxsys%nax*nsp2, pc_grp(1), msgid, supergroup )
    elseif (grandparent) then 
      msgid = 11
      call mp_recv(gradients(:,:,1:nsp2),3*maxsys%nax*nsp2, pc_grp(2), msgid, supergroup )
      fion(:,:,1:nsp2) =  fion(:,:,1:nsp2) + gradients(:,:,1:nsp2)
   end if
    
    call mp_bcast(fion(:,:,1:nsp2),3*maxsys%nax*nsp2,parai%io_source,supergroup) 
    call mp_sync(parai%cp_grp)

    DEALLOCATE(gradients, STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', __LINE__,__FILE__)
     
    RETURN
  end subroutine CpmdEnergiesGradientsPrivate

  subroutine ApplyWallPotential(self, gradients)
    type (BicanonicalCpmdType), intent(in) :: self
!      
!     Variable atom number set repulsive potential
!     We avoid collisions of atoms in case of no 
!     interaction (zero force in force weighting).
!
!     JF must be called from the large canonical systen (granparent)
!
    character(*), parameter                  :: procedureN = 'ApplyWallPotential'
    real(real_8) :: gradients(3,maxsys%nax,*)
    real(real_8), allocatable :: TSCR(:,:)
  integer ::   ia, ib, isp, iat, is, ierr
  integer,save :: icall = 0 
  integer,save :: nAtomInTSCR_cnf1, nAtomInTSCR_cnf2
  integer,save :: nsp_cnf1, nsp_cnf2
  real(real_8) :: xa(3), xb(3), dab(3), cab(3), eab(3), dd, vwall, dab_(3)
  LOGICAL ::  applywall

 !if (paral%io_parent) write(*,*) 'Entering ', procedureN
  if (.not.biCanonicalEnsembleDo) return
  allocate(tscr(3, maxsys%nax*maxsys%nsx),STAT=ierr)
  if(ierr/=0) call stopgm(procedureN,'allocation problem', __LINE__,__FILE__)
  call zeroing(tscr)
  icall = icall + 1  
  IF ((icall .eq. 1) .and. paral%io_parent) then
    ! we have a mapping of na(is),nsx -> na(is)*nsx  
    ! nsp_cnf1*na(is) and  nsp_cnf2*na(is) ; 
    !
    nsp_cnf1 = maxsys%nsx
    nsp_cnf2 = GetNumberOfSpeciesInSecondCanonicalSystem(self%bicanonicalCalculation%config)
    iat = 0
    do is = 1, nsp_cnf2 
      do ia = 1, ions0%na(is)
        iat = iat + 1
     end do
    end do
    nAtomInTSCR_cnf2 = iat
    do is = nsp_cnf2 + 1, nsp_cnf1
      do ia = 1, ions0%na(is)
        iat = iat + 1
     end do
    end do
    nAtomInTSCR_cnf1 = iat
    !
    if (PrintDebug(self%bicanonicalCalculation)) then
      write (OUT_TAPE,  '(/,2x,a)')' variable atom number wall constraints enabled: will use '
      write (OUT_TAPE, *) "nAtom cnf1: ", nAtomInTSCR_cnf1
      write (OUT_TAPE, *) "nAtom cnf2: ", nAtomInTSCR_cnf2
   end if
  end if
  ! Map coordinates from tau0(1:3,nat,nsp) to tscr(1:3,nat_tot)  
  call dum2(tau0,tscr)
  ! Loop over atoms present in both configurations 
  do ia = 1, nAtomInTSCR_cnf2
  ! Loop over atoms present in the larger configuration exclusively
    do ib = nAtomInTSCR_cnf2 + 1, nAtomInTSCR_cnf1
    !jf dbg write (*,*) ia, ib
    ! set wall parameters :
    ! set wall cutoff distance         
    cab(2) = 1.0_real_8/fbohr
    ! set wall apply switch bond distance         
    cab(1) = cab(2)/10.0_real_8
    ! potential height  
    cab(3) = 1.0e-3_real_8 
    if(ib.eq.ia) call stopgm(procedureN, 'Distance between like atoms specified',& 
      __LINE__,__FILE__)
    ! do coordinate mapping
    call fillc(ia, tscr, xa)
    call fillc(ib, tscr, xb)
    dab(1:3) = xa(1:3) - xb(1:3)
    ! IF ( .NOT. TISOS) 
    dab_(:) = dab(:)
    call pbc(dab_(1), dab_(2), dab_(3), dab(1), dab(2), dab(3), 1, parm%apbc, parm%ibrav)
    dd = dab(1)*dab(1) + dab(2)*dab(2) + dab(3)*dab(3)
    dd = dsqrt(dd)
    !jf dbg write (*,*) "distance ",DD
    if (dd .lt. 1.0e-6_real_8) call stopgm(procedureN, 'Too small distance between atoms',& 
          __LINE__,__FILE__)
    applyWall = .false.
    if (dd .le. cab(1)) applyWall = .true.
    vwall = 0.0_real_8
    xa(1) = dd !just a back up of DD
    if (cab(3) .lt. 0.0_real_8) dd = -1.0_real_8*dd
    if (applyWall) then     !Rij>Rw
      eab(1) = dab(1)/dd    !eij
      eab(2) = dab(2)/dd    !eij
      eab(3) = dab(3)/dd    !eij
      dd = exp(-cab(2)*(dd - cab(1)))
      vwall = cab(3)*(1.0_real_8 - dd)**(2.0_real_8)
      IF (paral%io_parent) &
      write(*,'(a,e18.4,a,2i8,a,F8.4,a,e18.4)')' variable atom number wall potential ',&
        vwall,' (a.u.) is active on atoms ', IA, IB, ' dist:', XA(1),&
        ' force ' , 2.0_real_8*cab(2)*cab(3)*dd*(1.0_real_8 - dd)
        vwall = 2.0_real_8*cab(2)*cab(3)*dd*(1.0_real_8 - dd)
        iat = CPAT(ib)
        isp = CPSP(ib)
        gradients(1:3, iat, isp) = gradients(1:3, iat, isp) + vwall*eab(1:3)
     end if
    end do      
 end do
  deallocate (tscr, STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', __LINE__,__FILE__)
  return
end subroutine ApplyWallPotential

subroutine PrintEnergiesPrivate(self, nfi)
  type (BicanonicalCpmdType), intent(inout) :: self
  CHARACTER(*), PARAMETER                  :: procedureN = 'PrintEnergiePrivate'
  integer, intent(in) :: nfi

    if (.not.biCanonicalEnsembleDo) return
    
    
    !JF TODO This is a temporary work around to keep CPMD running
    !        if BicanonicalEnsemble runtype was requested but with
    !        providing only one single MPI thread. 
    !        In this case we perfom a Canonical ensemble calculation
    !        for the big system, only. 
    ! TODO   The number of MD steps will be set   NFI=1.
    !
    !write (OUT_TAPE,*) "parai%cp_nproc", parai%cp_nproc
    if (parai%cp_nproc == 1 ) then
      write (OUT_TAPE, *) 
      write (OUT_TAPE, *) "!!! Warning !!!",&
         "BicanonicalEnsemble runtype was requested, ",&
         "providing only one mpi thread, ",&
         "but the implementation requires two or more. ",&
         "We proceed, but with the results being ",& 
         "meaningless and are",&
         " _NOT_ those of the BicanonicalEnsemble."
       !cnti%nomore=1
      write (OUT_TAPE, *)  
      return
   end if

    if (grandparent) then
      if ( self%bicanonicalCalculation%config%debug) then
        write (OUT_TAPE, "(4x,'chempot_mu(x=0.5) = ', e20.8)") &
          GetChemicalPotentialAtXWeightOneHalf(self%bicanonicalCalculation)
        write(OUT_TAPE, "(4x,'mu    = ', e20.8)") GetChemicalPotential(self%bicanonicalCalculation)
        write(OUT_TAPE, "(4x,'x     = ', e20.8)") GetXWeight(self%bicanonicalCalculation)
        write(OUT_TAPE, "(4x,'dN = ' ,i4,               4x,'dU = '  , f20.8)")&
          !JF TODO this is hardcoded and may be changed 
        1,&
          self%bicanonicalCalculation%energyCnf1 - self%bicanonicalCalculation%energyCnf2
     end if
      write (getEnergiesVarAtTape(self), '(i10, e20.12, f16.8, f20.10, e20.12)')&
        nfi, GetChemicalPotential(self%bicanonicalCalculation),& 
        GetXWeight(self%bicanonicalCalculation),&
        self%bicanonicalCalculation%energyCnf1 - self%bicanonicalCalculation%energyCnf2,&
        GetChemicalPotentialAtXWeightOneHalf(self%bicanonicalCalculation)
    endif
  end subroutine PrintEnergiesPrivate
  
  subroutine bicanonicalInitIOPrivate(self)
    type (BicanonicalCpmdType) :: self
    integer                                  :: n1, n2
    LOGICAL                                  :: ferror
   
    self%energiesVarAtTape = VARIABLE_ATOM_NUMBER_ENERGIES_TAPE
    self%nameEnergiesVarAtTape = NAME_VARIABLE_ATOM_NUMBER_ENERGIES_TAPE
    if  (grandparent) &
    call fileOpen(self%energiesVarAtTape, self%nameEnergiesVarAtTape, fo_app+fo_verb, fError)

    !JF the tape id are the same as used by other cpmd routines 
    self%energiesTape = 3
    self%trajecTape = 8
    self%trajectoryTape = 4
    self%constraintTape = 31 
    self%ftrajectoryTape = 4 
    self%metricTape = 90
    self%noseTrajecTape = 23
    self%noseEnergyTape = 26
    self%latestTape = 26
    self%geometryTape = 12
    self%geometryXyzTape = 12
    self%geometryScaleTape =  12
    call bicanonical_filename("RESTART_CNF", self%nameRestartTape, GetEnsembleIdPrivate(self))
    call xstring(self%nameRestartTape, n1, n2)
    self%nameRestartTape = self%nameRestartTape(n1:n2)//'.'
    call xstring(self%nameRestartTape, n1, n2)
    call bicanonical_filename("ENERGIES_CNF", self%nameEnergiesTape, GetEnsembleIdPrivate(self))
    call bicanonical_filename("TRAJEC_CNF", self%nameTrajecTape, GetEnsembleIdPrivate(self))
    call xstring(self%nameTrajecTape, n1, n2)
    self%nameTrajecTape = self%nameTrajecTape(n1:n2)//'.xyz'
    call bicanonical_filename("TRAJECTORY_CNF", self%nameTrajectoryTape, GetEnsembleIdPrivate(self))
    call bicanonical_filename("CONSTRAINT_CNF", self%nameConstraintTape, GetEnsembleIdPrivate(self))
    call bicanonical_filename("FTRAJECTORY_CNF", self%nameFtrajectoryTape, GetEnsembleIdPrivate(self))
    call bicanonical_filename("METRIC_CNF", self%nameMetricTape, GetEnsembleIdPrivate(self))
    call bicanonical_filename("NOSE_ENERGY_CNF", self%nameNoseTrajecTape, GetEnsembleIdPrivate(self))
    call bicanonical_filename("NOSE_TRAJEC_CNF", self%nameNoseEnergyTape, GetEnsembleIdPrivate(self))
    call bicanonical_filename("LATEST_CNF", self%nameLatestTape, GetEnsembleIdPrivate(self))
    call bicanonical_filename("GEOMETRY_CNF", self%nameGeometryTape, GetEnsembleIdPrivate(self))
    call bicanonical_filename("GEOMETRY_CNF", self%nameGeometryXyzTape, GetEnsembleIdPrivate(self))
    call xstring(self%nameGeometryXyzTape, n1, n2)
    self%nameGeometryXyzTape = self%nameGeometryXyzTape(n1:n2)//'.xyz'
    call bicanonical_filename("GEOMETRY_CNF", self%nameGeometryScaleTape, GetEnsembleIdPrivate(self))
    call xstring(self%nameGeometryScaleTape, n1, n2)
    self%nameGeometryScaleTape = self%nameGeometryScaleTape(n1:n2)//'.scale'
  end subroutine bicanonicalInitIOPrivate

  integer function getEnergiesVarAtTapePrivate(self) 
    type (BicanonicalCpmdType), intent(in) :: self
    getEnergiesVarAtTapePrivate = VARIABLE_ATOM_NUMBER_ENERGIES_TAPE
  end function getEnergiesVarAtTapePrivate

  function getNameEnergiesVarAtTapePrivate(self) 
    character(20) :: getNameEnergiesVarAtTapePrivate
    type (BicanonicalCpmdType), intent(in) :: self
    getNameEnergiesVarAtTapePrivate = NAME_VARIABLE_ATOM_NUMBER_ENERGIES_TAPE
  end function getNameEnergiesVarAtTapePrivate

  function getNameRestartTapePrivate(self) 
    character(20) :: getNameRestartTapePrivate
    type (BicanonicalCpmdType), intent(in) :: self
    getNameRestartTapePrivate = self%nameRestartTape
  end function getNameRestartTapePrivate

  function getNameLatestTapePrivate(self) 
    character(20) :: getNameLatestTapePrivate
    type (BicanonicalCpmdType), intent(in) :: self
    getNameLatestTapePrivate = self%nameLatestTape
  end function getNameLatestTapePrivate

  function getNameEnergiesTapePrivate(self) 
    character(20) :: getNameEnergiesTapePrivate
    type (BicanonicalCpmdType), intent(in) :: self
    getNameEnergiesTapePrivate = self%nameEnergiesTape
  end function getNameEnergiesTapePrivate

  function getNameTrajecTapePrivate(self) 
    character(20) :: getNameTrajecTapePrivate
    type (BicanonicalCpmdType), intent(in) :: self
    getNameTrajecTapePrivate = self%nameTrajecTape
  end function getNameTrajecTapePrivate

  function getNameTrajectoryTapePrivate(self) 
    character(20) :: getNameTrajectoryTapePrivate
    type (BicanonicalCpmdType), intent(in) :: self
    getNameTrajectoryTapePrivate = self%nameTrajectoryTape
  end function getNameTrajectoryTapePrivate

  function getNameConstraintTapePrivate(self) 
    character(20) :: getNameConstraintTapePrivate
    type (BicanonicalCpmdType), intent(in) :: self
    getNameConstraintTapePrivate = self%nameConstraintTape
  end function getNameConstraintTapePrivate

  function getNameFtrajectoryTapePrivate(self) 
    character(20) :: getNameFtrajectoryTapePrivate
    type (BicanonicalCpmdType), intent(in) :: self
    getNameFtrajectoryTapePrivate = self%nameFtrajectoryTape
  end function getNameFtrajectoryTapePrivate

  function getNameMetricTapePrivate(self) 
    character(20) :: getNameMetricTapePrivate
    type (BicanonicalCpmdType), intent(in) :: self
    getNameMetricTapePrivate = self%nameMetricTape
  end function getNameMetricTapePrivate

  function getNameNoseTrajecTapePrivate(self) 
    character(20) :: getNameNoseTrajecTapePrivate
    type (BicanonicalCpmdType), intent(in) :: self
    getNameNoseTrajecTapePrivate = self%nameNoseTrajecTape
  end function getNameNoseTrajecTapePrivate

  function getNameNoseEnergyTapePrivate(self) 
    character(20) :: getNameNoseEnergyTapePrivate
    type (BicanonicalCpmdType), intent(in) :: self
    getNameNoseEnergyTapePrivate = self%nameNoseEnergyTape
  end function getNameNoseEnergyTapePrivate

  function getNameGeometryXyzTapePrivate(self) 
    character(20) :: getNameGeometryXyzTapePrivate
    type (BicanonicalCpmdType), intent(in) :: self
    getNameGeometryXyzTapePrivate = self%nameGeometryXyzTape
  end function getNameGeometryXyzTapePrivate

  function getNameGeometryTapePrivate(self) 
    character(20) :: getNameGeometryTapePrivate
    type (BicanonicalCpmdType), intent(in) :: self
    getNameGeometryTapePrivate = self%nameGeometryTape
  end function getNameGeometryTapePrivate

  function getNameGeometryScaleTapePrivate(self) 
    character(20) :: getNameGeometryScaleTapePrivate
    type (BicanonicalCpmdType), intent(in) :: self
    getNameGeometryScaleTapePrivate = self%nameGeometryScaleTape
  end function getNameGeometryScaleTapePrivate

  logical function equalPrivate(array1, array2)
    integer, dimension(:), intent(in) :: array1, array2
    integer :: i

    equalPrivate =size(array1) == size(array2)
    if ( equalPrivate ) then
      do i = 1,size(array1)
      equalPrivate = array1(i) == array2(i)
      if ( .not. equalPrivate )exit
      enddo
    endif
  end function equalPrivate

    logical function globalEqualPrivate(ARRAY1 , ARRAY2)
      real(real_8), dimension(:,:,:), intent(in) :: array1, array2
      integer :: j
      globalEqualPrivate = .true.
      do j = 1,size(ARRAY1,3)
        if(all(abs(ARRAY1(:,:,j) - ARRAY2(:,:,j)) < 2*epsilon(0.0_real_8))) cycle
        globalEqualPrivate = .false.
        write (*,*) j, abs(ARRAY1(:,:,j) - ARRAY2(:,:,j)),2*epsilon(0.0_real_8)
        exit
     end do
    end function globalEqualPrivate

  end module bicanonicalCpmd
