module bicanonicalInputReader
      
      ! reads the input of bicanonical ensemble and sets 
      ! it to bicanonicalConfigType

  USE error_handling,                  ONLY: stopgm
  USE kinds, ONLY:  real_8
  use bicanonicalConfig
  implicit none 

  private 
  public bicanonicalInputReaderType
  public New, Delete 
  public ReadInput
  save 

  type bicanonicalInputReaderType
    private 
   ! integer(KINT)                                     :: inputFileUnit
    integer                                           :: inputFileUnit
 end type bicanonicalInputReaderType

  interface New 
    module procedure NewPrivate
 end interface New

  interface Delete
    module procedure DeletePrivate
 end interface Delete

  interface ReadInput
    module procedure ReadInputPrivate
 end interface ReadInput
 
contains 
 
  subroutine NewPrivate(self, inputUnit)
    type (bicanonicalInputReaderType), intent(out) :: self
    !integer(KINT), intent(in)                         :: inputUnit 
    integer, intent(in)                               :: inputUnit 
    self%inputFileUnit = inputUnit
  end subroutine NewPrivate

  subroutine DeletePrivate(self)
    type (bicanonicalInputReaderType), intent(inout) :: self
  end subroutine DeletePrivate
  
  subroutine ReadInputPrivate(self, config)
    type (bicanonicalInputReaderType), intent(inout) :: self
    type (bicanonicalConfigType), intent(inout) :: config
    integer, parameter                          :: linelen = 80
    character(len=linelen)                      :: line

    
      config%calculationConfig%numberOfCanonicalSystems = 2
      backspace(self%inputFileUnit)
      read(self%inputFileUnit, err=1698, end=1699, fmt='(A80)') line
      if ( index(line,'CHEMICALPOT').ne.0 ) then 
        config%calculationConfig%doConstantWeight = .false.
        config%runType = CONSTANT_CHEMICAL_POTENTIAL_BICANONICAL_ENSEMBLE_RUNTYPE
      elseif ( index(line,'XWEIGHT').NE.0 ) then
        config%calculationConfig%doConstantWeight = .true.
        config%runType = CONSTANT_XWEIGHT_BICANONICAL_ENSEMBLE_RUNTYPE
      else
        call stopgm('bicanonicalInputReader', 'ReadInputPrivate: Must specify Runtype',&
 __LINE__,__FILE__)    
     end if
      if ( index(line,'INFO').NE.0) config%calculationConfig%debug = .true.
       
      IF (config%calculationConfig%doConstantWeight) then
        read (self%inputFileUnit, ERR=1698, END=1699, FMT=*) &
          config%calculationConfig%xWeight,&
          config%calculationConfig%temperatureGammaFactor
        if (abs(config%calculationConfig%xWeight - 0.5_real_8 ) >= (0.5_real_8 + 2_real_8*epsilon(0.0_real_8)) ) &
        CALL STOPGM('bicanonicalInputReader', 'Input error: xWeight out of bounds [0.0, 1.00]',&
 __LINE__,__FILE__)     
      else 
        read (self%inputFileUnit, ERR=1698, END=1699, FMT=*) &
          config%calculationConfig%chemicalPotential, &
          config%calculationConfig%temperatureGammaFactor
     end if
      return
1698  CALL STOPGM('bicanonicalInputReader', 'ReadInputPrivate: Error on reading input ', &
__LINE__,__FILE__)     
1699  CALL STOPGM('bicanonicalInputReader', 'ReadInputPrivate: End of file ', &
__LINE__,__FILE__)     
    end subroutine ReadInputPrivate

  end module bicanonicalInputReader
