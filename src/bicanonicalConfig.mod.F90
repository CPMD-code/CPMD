module bicanonicalConfig
  
  USE error_handling,                  ONLY: stopgm
  use bicanonicalCalculationConfig
  implicit none 

  private 
  public bicanonicalConfigType
  public New, Delete, Print 
  public GetNumberOfCanonicalSystems
  public NOT_BICANONICAL_ENSEMBLE_RUNTYPE
  public CONSTANT_CHEMICAL_POTENTIAL_BICANONICAL_ENSEMBLE_RUNTYPE 
  public CONSTANT_XWEIGHT_BICANONICAL_ENSEMBLE_RUNTYPE
  !JF TODO
!  public SetBicanonicalConfigCpmd
!  public GetBicanonicalConfigCpmd
  !save 

  ! run types 
  integer, parameter                                  :: NOT_BICANONICAL_ENSEMBLE_RUNTYPE = 0
  integer, parameter                                  :: CONSTANT_CHEMICAL_POTENTIAL_BICANONICAL_ENSEMBLE_RUNTYPE = 1
  integer, parameter                                  :: CONSTANT_XWEIGHT_BICANONICAL_ENSEMBLE_RUNTYPE = 2
  type bicanonicalConfigType
   ! integer(KINT)                                     :: runType
    integer                                           :: runType
    type (bicanonicalCalculationConfigType)    :: calculationConfig
 end type bicanonicalConfigType

  interface New 
    module procedure NewWithDefaultsPrivate!, NewWithCpmdInputPrivate
 end interface New

  interface Delete
    module procedure DeletePrivate
 end interface Delete

  interface Print
    module procedure PrintPrivate
 end interface Print
 
contains 
 
  subroutine NewWithDefaultsPrivate(self)
    type (bicanonicalConfigType), intent(out)    :: self
    !Defaults
    self%runType = NOT_BICANONICAL_ENSEMBLE_RUNTYPE
    call New(self%calculationConfig)
  end subroutine NewWithDefaultsPrivate

!  subroutine NewWithCpmdInputPrivate(self, cpmdConfig)
!    type (bicanonicalConfigType), intent(out)    :: self
!    type (bicanonicalConfigType), intent(in)     :: cpmdConfig
!    !Defaults
!    self = cpmdConfig
  !  end subroutine

  subroutine DeletePrivate(self)
    type (bicanonicalConfigType), intent(inout)  :: self
    call Delete(self%calculationConfig)
  end subroutine DeletePrivate

  function RunTypeAsString(runType) result (runTypeString)
    !integer(KINT), intent(in)                  :: runType
    integer, intent(in)                        :: runType
    ! character(LCHARS)                         :: runTypeString
    character(len=255)                 :: runTypeString
    select case (runType)
    case (CONSTANT_CHEMICAL_POTENTIAL_BICANONICAL_ENSEMBLE_RUNTYPE)
      runTypeString =  'constant chemical potential'
    case (CONSTANT_XWEIGHT_BICANONICAL_ENSEMBLE_RUNTYPE)
      runTypeString = 'constant bicanonical x weight'
    case (NOT_BICANONICAL_ENSEMBLE_RUNTYPE)
      runTypeString = 'Bicanonical ensemble disabled'
    case default
       call STOPGM('bicanonicalConfig', &
            'Invalid run type in bicanonicalConfig', &
            __LINE__,__FILE__)
    end select
  end function RunTypeAsString

  subroutine PrintPrivate(self, unit)
    type (bicanonicalConfigType), intent(in)   :: self
    !integer(KINT), intent(in)                         :: unit
    integer, intent(in)                               :: unit
    
    if ( self%runType /= NOT_BICANONICAL_ENSEMBLE_RUNTYPE) then
      write(unit,*)
      write(unit, '(";----------------------------------------;")')
      write(unit, '("| General bicanoncal ensemble parameters |")')
      write(unit, '(";----------------------------------------;")')
      write(unit,*)
      write(unit, '(a)')'Sampling bicanonical ensemble at'
      write(unit, '(a,t40,a)')'Run Type:', trim( RunTypeAsString(self%runType) )
      write(unit,*)
      call Print(self%calculationConfig, unit)

   end if
  end subroutine PrintPrivate

  integer function GetNumberOfCanonicalSystems(self)
    type (bicanonicalConfigType), intent(in)   :: self
    GetNumberOfCanonicalSystems = self%calculationConfig%numberOfCanonicalSystems
  end function GetNumberOfCanonicalSystems

end module bicanonicalConfig
