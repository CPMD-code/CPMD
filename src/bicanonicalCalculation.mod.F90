module bicanonicalCalculation
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  use bicanonicalCalculationConfig
!  use multiCanonicalSystem
  implicit none

  private
  save
  public bicanonicalCalculationType
  public New, Delete
  public Energies
  public bicanonicalGradients
  public updateMuXweight
  public GetChemicalPotential
  public GetChemicalPotentialAtXWeightOneHalf
  public GetXWeight
  public Print
  public PrintDebug
!  public PrintResults, PrintConfig
!  public GetEnergy, GetGradient
   
  type bicanonicalCalculationType
    type (bicanonicalCalculationConfigType)      :: config
!    type (MultiCanonicalSystemType), pointer            :: multiCanonicalSystem
    !real(real_8), pointer                               :: gradients(:,:,:)
    real(real_8)                                        :: gammaParameter, beta
    real(real_8)                                        :: chemicalPotential, xWeight
    real(real_8)                                        :: chemicalPotentialAtXWeightOneHalf
    real(real_8)                                        :: energyCnf1, energyCnf2
 end type bicanonicalCalculationType

  interface New 
    module procedure NewPrivate
 end interface New

  interface Delete
    module procedure DeletePrivate
 end interface Delete

  interface Energies
    module procedure EnergiesPrivate
 end interface Energies

  interface bicanonicalGradients
    module procedure bicanonicalGradientsPrivate
 end interface bicanonicalGradients

  interface updateMuXweight
    module procedure updateMuXweightPrivate
 end interface updateMuXweight

  interface GetXWeight
    module procedure GetXWeightPrivate
 end interface GetXWeight

  interface GetChemicalPotential
    module procedure GetChemicalPotentialPrivate
 end interface GetChemicalPotential

  interface GetChemicalPotentialAtXWeightOneHalf
    module procedure GetChemicalPotentialAtXWeightOneHalfPrivate
 end interface GetChemicalPotentialAtXWeightOneHalf
  
  interface Print
    module procedure PrintPrivate
 end interface Print

  interface PrintDebug
    module procedure PrintDebugPrivate
 end interface PrintDebug


 

contains

  subroutine NewPrivate(self, config)
    type (bicanonicalCalculationType), intent(out) :: self
    type (bicanonicalCalculationConfigType), intent(in) :: config
    
    self%config = config 
    self%chemicalPotential = self%config%chemicalPotential
    self%xWeight = self%config%xWeight
    call CalcGammaParameterPrivate(self)
  end subroutine NewPrivate

  subroutine DeletePrivate(self)
    type (bicanonicalCalculationType), intent(inout) :: self
  end subroutine DeletePrivate

  subroutine bicanonicalGradientsPrivate(self, gradients)
    type (bicanonicalCalculationType), intent(inout) :: self 
    real(real_8),intent(inout) :: gradients(:,:,:)
       
    gradients(:,:,:) = GetXWeightPrivate(self)*gradients(:,:,:)
    
  end subroutine bicanonicalGradientsPrivate

  subroutine EnergiesPrivate(self)
    type (bicanonicalCalculationType), intent(inout) :: self
    
    call  updateMuXweight(self)
    
  end subroutine EnergiesPrivate

  subroutine UpdateMuXweightPrivate(self)
    type (bicanonicalCalculationType), intent(inout) :: self 

    if (self%config%doConstantWeight) then
      call CalcMuPrivate(self) 
    else
      call CalcXweightPrivate(self)
   end if

    call CalcMuAtXWeightOneHalfPrivate(self)

  end subroutine UpdateMuXweightPrivate

  subroutine CalcXweightPrivate(self)
    type (bicanonicalCalculationType), intent(inout) :: self 
      
      self%xWeight = 1.0_real_8/(1.0_real_8 + self%gammaParameter*&
        dexp(self%beta*(self%energyCnf1 - self%energyCnf2 - self%chemicalPotential)))
    
    end subroutine CalcXweightPrivate

  subroutine CalcMuPrivate(self)
    type (bicanonicalCalculationType), intent(inout) :: self 

!jf   if CXEQWGT is 0.0 or 1.0, i.e. pure systems the  CHEMPOT_MU
!        formula has poles. 
      if ( self%xWeight<=epsilon(0.0_real_8) )  then
        self%chemicalPotential  =  (self%energyCnf1 - self%energyCnf2) +  huge(1.0_real_8)
      elseif ( self%xWeight>=1.0_real_8 - epsilon(0.0_real_8)) then
        self%chemicalPotential =  (self%energyCnf1 - self%energyCnf2) - huge(1.0_real_8)
      else
        self%chemicalPotential =  (self%energyCnf1 - self%energyCnf2)&
          +  1.0_real_8/self%beta*dlog(self%gammaParameter*self%XWeight/(1.0_real_8 - self%XWeight))
     end if
    
   end subroutine CalcMuPrivate

  subroutine CalcMuAtXWeightOneHalfPrivate(self)
    type (bicanonicalCalculationType), intent(inout) :: self 
    self%chemicalPotentialAtXWeightOneHalf =  (self%energyCnf1 - self%energyCnf2)&
      +  1.0_real_8/self%beta*dlog(self%gammaParameter)
  end subroutine CalcMuAtXWeightOneHalfPrivate

  subroutine CalcGammaParameterPrivate(self)
    type (bicanonicalCalculationType), intent(inout) :: self
    real(real_8), parameter :: hPlanck = 6.6260755e-34_real_8,&
                               kboltz = 0.316681534e-5_real_8,&
                               fbohr = 1._real_8/0.529177210859_real_8
    !                           ry = 13.60569193_real_8,&
    !                           evtokel = 11604.505_real_8
    !                           scmass = 1822.888485_real_8, &
    !                           au_fs= 2.418884326505e-2_real_8,&
    !                           factem = 2._real_8*ry*evtokel
    real(real_8), parameter :: pi = dacos(-1.0_real_8)
    real(real_8) :: lambda 
    
    real(real_8) :: mass
    
    ! see J. Frenzel, B. Meyer, D. Marx (2017)
    !
    mass = self%config%atomMassOfExcessSpecies 
    ! Lambda in SI              
    lambda =1.0e10_real_8*fBohr*hPlanck/dsqrt(2.0_real_8*pi*mass*&
       1.6605405e-27_real_8*1.380658e-23_real_8*self%config%temperatureGammaFactor)
    self%gammaParameter = lambda**3*1.0_real_8/self%config%volumeSimulationBox
    self%beta = 1.0_real_8/(kBoltz*self%config%temperatureGammaFactor)!/factem)
  end subroutine CalcGammaParameterPrivate

  real(real_8) function GetXWeightPrivate(self)
    type (bicanonicalCalculationType), intent(inout) :: self 
    
    if (self%config%idCanonicalEnsemble == 1 ) then 
      GetXWeightPrivate = self%xWeight
    else
     GetXWeightPrivate = 1.0_real_8 - self%xWeight 
    endif
  end function GetXWeightPrivate

  real(real_8) function GetChemicalPotentialPrivate(self)
    type (bicanonicalCalculationType), intent(inout) :: self 
    GetChemicalPotentialPrivate = self%chemicalPotential
  end function GetChemicalPotentialPrivate

  real(real_8) function GetChemicalPotentialAtXWeightOneHalfPrivate(self)
    type (bicanonicalCalculationType), intent(inout) :: self 
    GetChemicalPotentialAtXWeightOneHalfPrivate = self%chemicalPotentialAtXWeightOneHalf
  end function GetChemicalPotentialAtXWeightOneHalfPrivate
  
  logical function IoNodePrivate(self)
    type (bicanonicalCalculationType), intent(inout) :: self 
    IoNodePrivate = self%config%ioNode
  end function IoNodePrivate

  logical function PrintDebugPrivate(self)
    type (bicanonicalCalculationType), intent(in) :: self 
    PrintDebugPrivate = self%config%debug
  end function PrintDebugPrivate

  subroutine PrintPrivate(self, unit)
    type (bicanonicalCalculationType), intent(in) :: self
    !integer(KINT), intent(in)                         :: unit
    integer, intent(in)                               :: unit
    
!      write(unit,*)
!      write(unit, '(";----------------------------------;")')
!      write(unit, '("| VARIABLE_ATOM_NUMBER CALCULATION |")')
!      write(unit, '(";----------------------------------;")')
!      write(unit,*)
      write (unit,'(a,t55,f10.5)') "beta =1/kT (Ha) ", self%beta
      write (unit,'(a          )') "        (T of variable atom number keyword)"
      write (unit,'(a,t55,e10.4)') "gamma = Lambda/PBCVolume",self%gammaParameter
      write (unit,'(a          )') "        (Lambda: free wavelength of exess species in system 1)"
      write (unit,'(a,t55,f10.5)') "-kT ln gamma (Ha)", -1._real_8/self%beta*dlog(self%gammaParameter)
      write(unit,*)
    end subroutine PrintPrivate
end module bicanonicalCalculation
