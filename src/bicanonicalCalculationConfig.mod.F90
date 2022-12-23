module bicanonicalCalculationConfig
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  
  implicit none 

  private 
  public bicanonicalCalculationConfigType
  public GetNumberOfSpeciesInSecondCanonicalSystem
  public New, Delete, Print 
  save 

  type bicanonicalCalculationConfigType
    logical                                           :: debug, ioNode 
    logical                                           :: doConstantWeight
    integer                                           :: numberOfCanonicalSystems
    integer                                           :: numberOfSpeciesInSecondCanonicalSystem
    integer                                           :: numberOfCanonicalEnsembles
    integer                                           :: idCanonicalEnsemble
    real(real_8)                                      :: chemicalPotential
    real(real_8)                                      :: xWeight
    real(real_8)                                      :: atomMassOfExcessSpecies
    real(real_8)                                      :: temperatureGammaFactor ! default: 300.0 K
    real(real_8)                                      :: volumeSimulationBox
 end type bicanonicalCalculationConfigType

  interface New 
    module procedure NewPrivate
 end interface New

  interface Delete
    module procedure DeletePrivate
 end interface Delete

  interface Print
    module procedure PrintPrivate
 end interface Print
 
  interface GetNumberOfSpeciesInSecondCanonicalSystem
    module procedure GetNumberOfSpeciesInSecondCanonicalSystemPrivate
 end interface GetNumberOfSpeciesInSecondCanonicalSystem
  
contains 
 
  subroutine NewPrivate(self)
    type (bicanonicalCalculationConfigType), intent(out)  :: self
    !Defaults
    self%debug = .false.
    self%ioNode = .false.
    self%doConstantWeight = .false.
    self%numberOfCanonicalSystems = 2
    self%numberOfSpeciesInSecondCanonicalSystem = 0
    self%chemicalPotential = 0.0_real_8
    self%xWeight = 1.0_real_8 
    self%atomMassOfExcessSpecies=1.0_real_8
    self%temperatureGammaFactor = 300.0_real_8
    self%volumeSimulationBox = 0.0_real_8
  end subroutine NewPrivate

  subroutine DeletePrivate(self)
    type (bicanonicalCalculationConfigType), intent(inout) :: self
  end subroutine DeletePrivate
  
  subroutine PrintPrivate(self, unit)
    type (bicanonicalCalculationConfigType), intent(in) :: self
    !integer(KINT), intent(in)                         :: unit
    integer, intent(in)                               :: unit
    
      write(unit,*)
      write(unit, '(";------------------------------------;")')
      write(unit, '("| bicanonical ensemble md parameters |")')
      write(unit, '(";------------------------------------;")')
      write(unit,*)
      write(unit, '(a,t55,i10)')'no. of canonical systems:', self%numberOfCanonicalSystems
      write(unit, '(a,t55,i10)')'ID  of canonical system :', self%idCanonicalEnsemble
      write(unit, '(a,t55,i10)')'no. of species in sys2:', self%numberOfSpeciesInSecondCanonicalSystem
      write(unit, '(a,t55,f10.4)')'atom mass of excess species in sys1:', self%atomMassOfExcessSpecies
      write(unit, '(a,t55,f10.2)')'temperature used in pre-factors (K):', self%temperatureGammaFactor
      write(unit, '(a,t55,f10.2)')'volume simulation box (bohr^3):', self%volumeSimulationBox
      if (self%doConstantWeight) then 
        write(unit, '(a,t55)')'Constant x (weight) MD run mode'
        write(unit, '(a,t55,f10.5)')'x (weight):', self%xWeight
      else
        write(unit, '(a)')'Constant chemical Potential MD run mode'
        write(unit, '(a,t55,f10.5)')'chemical potential:', self%chemicalPotential
     end if
      write(unit,*)
    end subroutine PrintPrivate

  integer function GetNumberOfSpeciesInSecondCanonicalSystemPrivate(self)
    type (bicanonicalCalculationConfigType), intent(in) :: self
    GetNumberOfSpeciesInSecondCanonicalSystemPrivate = self%numberOfSpeciesInSecondCanonicalSystem
  end function GetNumberOfSpeciesInSecondCanonicalSystemPrivate
  
end module bicanonicalCalculationConfig
