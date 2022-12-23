MODULE epot_types
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! .....Declarations
  INTEGER, PARAMETER :: maxind =100 
  INTEGER, PARAMETER :: maxnear=2000 
  INTEGER, PARAMETER :: maxintm=5000 
  INTEGER, PARAMETER :: maxoutr=90000 
  ! .....
  ! .....General data

  ! .....
  ! .....Data for the local expansions

  ! .....Number of near atoms to include




  ! .....List of near atoms to add to each local expansion

  ! 
  ! .....Position of the Local Taylor Expansion 


  ! 
  ! .....Koeffs of the local expansion
  ! .....1   =PHI, 2-4=KX, KY, KZ
  ! .....5-7 =HESSXX,HESSYY,HESSZZ
  ! .....8-10=HESSXY,HESSXZ,HESSYZ


  ! .....Partial Charge, Core Charge, Sigma, CutOffWidth, CutOff



  ! .....FORCES ON MM_ATOMS DUE TO A DENSITY DISTRIBUTION OF THE QM SYSTEM


  ! .....
  ! .....Common blocks             
  ! ==================================================================
  TYPE :: epot1_t
     INTEGER :: nnear
     INTEGER :: nloc
     INTEGER :: locid(maxind)
     INTEGER :: nincl(maxind)
     INTEGER :: idnear(maxnear)
     INTEGER :: lincl(maxind,maxnear)
  END TYPE epot1_t
  TYPE(epot1_t) :: epot1
  TYPE :: epot2_t
     REAL(real_8) :: boxqm(6)
     REAL(real_8) :: boxdum(6)
     REAL(real_8) :: koloc(3,maxind)
     REAL(real_8) :: locexp(10,maxind)
     REAL(real_8) :: datnear(5,maxnear)
     REAL(real_8) :: konear(3,maxnear)
     REAL(real_8) :: myrag(100)
  END TYPE epot2_t
  TYPE(epot2_t) :: epot2
  TYPE :: epot3_t
     INTEGER :: stpno
     INTEGER :: pwlmax
     INTEGER :: intmf
     INTEGER :: outmf
     INTEGER :: nintm
     INTEGER :: noutr
     INTEGER :: idintm(maxintm)
     INTEGER :: idoutr(maxoutr)
  END TYPE epot3_t
  TYPE(epot3_t) :: epot3
  TYPE :: epot4_t
     REAL(real_8) :: zmmnear(maxnear)
     REAL(real_8) :: qintm(5,maxintm)
     REAL(real_8) :: kointm(3,maxintm)
     REAL(real_8) :: qoutr(5,maxoutr)
     REAL(real_8) :: kooutr(3,maxoutr)
  END TYPE epot4_t
  TYPE(epot4_t) :: epot4
  TYPE :: epot5_t
     REAL(real_8) :: gmx_fnear(3,maxnear)
     REAL(real_8) :: gmx_fintm(3,maxintm)
     REAL(real_8) :: gmx_foutr(3,maxoutr)
  END TYPE epot5_t
  TYPE(epot5_t) :: epot5
  ! ..End

END MODULE epot_types
