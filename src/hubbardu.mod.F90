MODULE hubbardu
  USE kinds,                           ONLY: real_8

!  ==================================================================
!  == INCLUDE FILE FOR DFT+U                                       ==
!  ==================================================================
!
!  DFTUP(ORTHO_PAW) -> orthogonal psuedo atomic orbitals
!  DFTUP(NORM_ORTHO_AW) -> normalized orthogonal atomic orbitals (for e.g. inital guess)
!  UATM stores the atom numbers of Habbarad atoms
!  NUATM number of Habbard atoms
!  L_OM is the total length of the array OM
!
!  ==================================================================
!
   IMPLICIT NONE

   INTEGER, PARAMETER :: maxuatm=150,maxluatm=10
   !INTEGER, PARAMETER :: iup_oddrtho=1,iup_norm_ortho=2,iup_norm=3,iup_direct=4
   TYPE :: hubbu_t
      LOGICAL      :: debug
      LOGICAL      :: uverb,tpom
      LOGICAL      :: portho,pnorm
      INTEGER      :: firstcall=0,pfrqom=0
      INTEGER      :: nuatm=0,l_om=0,uatm(maxuatm)=0,nuproj=0
      INTEGER      :: muatm(2,2,maxuatm)=0,nl(maxuatm)=0
      INTEGER      :: l(maxuatm,maxluatm)=0,s(maxuatm,maxluatm)=0
      REAL(real_8) :: ehub = 0.0_real_8
      REAL(real_8) :: u(maxuatm)= 0.0_real_8,a(maxuatm)= 0.0_real_8
      REAL(real_8) :: hs(maxuatm)= 0.0_real_8,bs(maxuatm)= 0.0_real_8
   END TYPE hubbu_t
   TYPE(hubbu_t) :: hubbu

!
!
   REAL(real_8),allocatable, save   :: om(:,:,:,:),fomi(:,:,:,:,:)
!
   COMPLEX(real_8),allocatable,save ::  mycatom(:,:),mycatom0(:,:)
   COMPLEX(real_8),allocatable,save ::  c2u0(:,:)
!
   REAL(real_8),allocatable,save    ::  fion_om(:,:,:)
!  was in atwf.inc in CPMD 3.x implementation
   REAL(real_8),allocatable,save    ::  atrg_u(:,:), atwfr_u(:,:,:),cat_u(:,:,:,:)
END MODULE hubbardu
