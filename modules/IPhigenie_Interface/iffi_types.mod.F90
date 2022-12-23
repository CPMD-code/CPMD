
MODULE iffi_types
  USE kinds,        ONLY: real_8
  
  IMPLICIT NONE
  
  INTEGER, PARAMETER :: maxind          = 50          !max number of QM atoms
  INTEGER, PARAMETER :: maxnear         = 5000       !max length of merged nearlist (=C_-3 'NEAR' and C_-2 'OLE' atoms from iffi), MXNO_MERGED_NEAR_ATOMS in iffi  SET ALSO IN epot.inc!!!!
  INTEGER, PARAMETER :: maxnearpervo    = 2000  !max number of near atoms (C_-3) per voxel
  INTEGER, PARAMETER :: maxolepervo     = 5000   !max number of ole atoms (C_-2) per voxel
  INTEGER, PARAMETER :: maxspecies      = 6
  INTEGER, PARAMETER :: mxno_threads    = 1    !max no of omp threads
  
  
  !.....Coefficients of the local expansion
!.....for readable access to LOCEXP entries
   INTEGER, PARAMETER :: KPHI       = 1
   INTEGER, PARAMETER :: KX         = 2
   INTEGER, PARAMETER :: KY         = 3
   INTEGER, PARAMETER :: KZ         = 4
   INTEGER, PARAMETER :: KXX        = 5
   INTEGER, PARAMETER :: KYY        = 6
   INTEGER, PARAMETER :: KXY        = 7
   INTEGER, PARAMETER :: KXZ        = 8
   INTEGER, PARAMETER :: KYZ        = 9
   INTEGER, PARAMETER :: KXXX       = 10
   INTEGER, PARAMETER :: KYYY       = 11
   INTEGER, PARAMETER :: KXYY       = 12
   INTEGER, PARAMETER :: KYXX       = 13
   INTEGER, PARAMETER :: KZXX       = 14
   INTEGER, PARAMETER :: KZYY       = 15
   INTEGER, PARAMETER :: KXYZ       = 16
   INTEGER, PARAMETER :: KXXYY      = 17
   INTEGER, PARAMETER :: KXXZZ      = 18
   INTEGER, PARAMETER :: KYYZZ      = 19
   INTEGER, PARAMETER :: KXXXY      = 20
   INTEGER, PARAMETER :: KXXXZ      = 21
   INTEGER, PARAMETER :: KYYYX      = 22
   INTEGER, PARAMETER :: KYYYZ      = 23
   INTEGER, PARAMETER :: KZZZX      = 24
   INTEGER, PARAMETER :: KZZZY      = 25
   INTEGER, PARAMETER :: dimlocexp  = KZZZY
   

   INTEGER, PARAMETER :: LAD        = 1
   INTEGER, PARAMETER :: PX         = 2
   INTEGER, PARAMETER :: PY         = 3
   INTEGER, PARAMETER :: PZ         = 4
   INTEGER, PARAMETER :: QXX        = 5
   INTEGER, PARAMETER :: QYY        = 6
   INTEGER, PARAMETER :: QXY        = 7
   INTEGER, PARAMETER :: QXZ        = 8
   INTEGER, PARAMETER :: QYZ        = 9
   INTEGER, PARAMETER :: OXYZ       = 10
   INTEGER, PARAMETER :: OXXY       = 11
   INTEGER, PARAMETER :: OXXZ       = 12
   INTEGER, PARAMETER :: OYYX       = 13
   INTEGER, PARAMETER :: OYYZ       = 14
   INTEGER, PARAMETER :: OZZX       = 15
   INTEGER, PARAMETER :: OZZY       = 16
   INTEGER, PARAMETER :: HXXYY      = 17
   INTEGER, PARAMETER :: HXXZZ      = 18
   INTEGER, PARAMETER :: HYYZZ      = 19
   INTEGER, PARAMETER :: HXXXY      = 20
   INTEGER, PARAMETER :: HYYYX      = 21
   INTEGER, PARAMETER :: HXXXZ      = 22
   INTEGER, PARAMETER :: HZZZX      = 23
   INTEGER, PARAMETER :: HYYYZ      = 24
   INTEGER, PARAMETER :: HZZZY      = 25
   INTEGER, PARAMETER :: dimmultexp = HZZZY
   
   INTEGER, PARAMETER :: N_KX       = 1
   INTEGER, PARAMETER :: N_KY       = 2
   INTEGER, PARAMETER :: N_KZ       = 3
   INTEGER, PARAMETER :: PCHARGE    = 4
   INTEGER, PARAMETER :: SIGCHARGE  = 5      
   INTEGER, PARAMETER :: SIGDIPOLE  = 6
   INTEGER, PARAMETER :: CORECHARGE = 7
   INTEGER, PARAMETER :: CO_WIDTH   = 8
   INTEGER, PARAMETER :: CO_VALUE   = 9
   INTEGER, PARAMETER :: sizeneardata = CO_VALUE

   INTEGER, PARAMETER :: N_PX       = 1
   INTEGER, PARAMETER :: N_PY       = 2
   INTEGER, PARAMETER :: N_PZ       = 3
   INTEGER, PARAMETER :: sizeneardata_var = N_PZ

   
!.....for readable access to MMPOTEXP entries
! ....first the fixed quantities..
   INTEGER, PARAMETER :: MPOT   = 1
   INTEGER, PARAMETER :: MEX    = 2
   INTEGER, PARAMETER :: MEY    = 3
   INTEGER, PARAMETER :: MEZ    = 4
   INTEGER, PARAMETER :: MEXX   = 5
   INTEGER, PARAMETER :: MEYY   = 6
   INTEGER, PARAMETER :: MEZZ   = 7
   INTEGER, PARAMETER :: MEXY   = 8
   INTEGER, PARAMETER :: MEXZ   = 9
   INTEGER, PARAMETER :: MEYZ   = 10
   INTEGER, PARAMETER :: dimpotexp  = MEYZ

! ....then the variable quantities..
   INTEGER, PARAMETER :: MEXPFF = 1
   INTEGER, PARAMETER :: MEYPFF = 2
   INTEGER, PARAMETER :: MEZPFF = 3
   INTEGER, PARAMETER :: dimpotexp_var  = MEZPFF
   

   INTEGER, PARAMETER :: RGZERO = 1
   INTEGER, PARAMETER :: RGX    = 2
   INTEGER, PARAMETER :: RGY    = 3
   INTEGER, PARAMETER :: RGZ    = 4
   INTEGER, PARAMETER :: RGXX   = 5
   INTEGER, PARAMETER :: RGYY   = 6
   INTEGER, PARAMETER :: RGZZ   = 7
   INTEGER, PARAMETER :: RGXY   = 8
   INTEGER, PARAMETER :: RGXZ   = 9
   INTEGER, PARAMETER :: RGYZ   = 10
   INTEGER, PARAMETER :: NGAMMA = 11
   INTEGER, PARAMETER :: dimrgyrexp = NGAMMA
  
  
  
!.....fixed quantities (=fixed during pff iteration)............................................................
!     interaction lists (changes each listupdate step)
   TYPE :: epot1_list_t
     INTEGER :: nnear            !number of near atoms
     INTEGER :: nloc             !number of local atomic expansions (=qm atoms)
     INTEGER :: locid(maxind)    !data for the local expansions
     INTEGER :: nincl(maxind)    !number of near atoms
     INTEGER :: idnear(maxnear)  !
     INTEGER :: lincl(maxind,maxnear) !list of near atoms to add to each local expansion
   END TYPE  epot1_list_t
   TYPE(epot1_list_t) epot1_list


   TYPE :: epot1_fix_t
     REAL(real_8) :: koloc(3,maxind)   ! Position of the Local Taylor Expansion (= qm atoms)
     REAL(real_8) :: NEARDATA    (SIZENEARDATA    ,MAXNEAR)  ! data of near atoms (charges,dipoles etc)   
     REAL(real_8) :: boxoffs(3)        ! vector from iffi to cpmd coordinate frame (needed for output of potential/rhoe)
     REAL(real_8) :: boxtrans(3)       ! translation vector when box has been recentered by IPHIGENIE
   END TYPE  epot1_fix_t
   TYPE(epot1_fix_t) epot1_fix
   



   
   
!.....variable quantities (=change during pff iteration)............................................................
   TYPE :: epot1_var_t
     LOGICAL :: TWRITERESTARTFILE
     LOGICAL :: TGRIDPART
     LOGICAL :: TDOEXTRAP
     LOGICAL :: TCALCONLYPOTFROMPFF
     LOGICAL :: TNOPOTANDWFCTUPD
     LOGICAL :: TCALCONLYPFFFIELD
     INTEGER :: CALCESPCHARGES
     LOGICAL :: TUPDATEQMINTLISTS
     LOGICAL :: TUPDATERGYR
     LOGICAL :: UPDATEMEANFIELD
     INTEGER :: NMEAN
   
     REAL(real_8) :: NEARDATA_VAR    (SIZENEARDATA_VAR    ,MAXNEAR) 
     REAL(real_8) :: LOCEXP(DIMLOCEXP,MAXIND)
   END TYPE  epot1_var_t
   TYPE(epot1_var_t) epot1_var
   
   
   
  
  !.....quantities that do not change during the calculation
  TYPE :: runinfo_t
    INTEGER :: nrqmatoms
    INTEGER :: sammtype
    LOGICAL :: pff
    INTEGER :: voxnum(3)
    LOGICAL :: normalmodeanalysis
    LOGICAL :: meanfieldmode
    LOGICAL :: tverbose
    
    REAL(real_8) :: boxqm(6)          !           
    REAL(real_8) :: boxdum(6)         ! box edges 
    REAL(real_8) :: OLERADIUS
    REAL(real_8) :: OLEEPS
    REAL(real_8) :: OLENSIG
    REAL(real_8) :: MYRAG(100)
  END TYPE runinfo_t
  TYPE(runinfo_t) :: runinfo

  
      REAL(real_8), ALLOCATABLE :: EXTFMM(:)
      REAL(real_8), ALLOCATABLE :: VOXEXP(:,:)
      REAL(real_8), ALLOCATABLE :: VOXEXP_MM(:,:)
      INTEGER     , ALLOCATABLE :: NEARLIST(:,:)
      INTEGER     , ALLOCATABLE :: NEARLISTLEN(:)
      INTEGER     , ALLOCATABLE :: OLELIST(:,:)
      INTEGER     , ALLOCATABLE :: OLELISTLEN(:)
      REAL(real_8), ALLOCATABLE :: QMPOTANDFLD(:,:)
      REAL(real_8), ALLOCATABLE :: ECHRG(:)
      REAL(real_8), ALLOCATABLE :: MMPOTEXP(:,:,:)
      REAL(real_8), ALLOCATABLE :: MMPOTEXP_VAR(:,:,:)
      REAL(real_8), ALLOCATABLE :: VOMULTEXP(:,:,:)
      REAL(real_8), ALLOCATABLE :: ATMULTEXP(:,:,:)
      REAL(real_8), ALLOCATABLE :: VOXRGYREXP(:,:,:)
      REAL(real_8), ALLOCATABLE :: ATRGYREXP(:,:,:)
      REAL(real_8), ALLOCATABLE :: EXTFMEAN(:)
      REAL(real_8), ALLOCATABLE :: EXTFLAST(:)
      REAL(real_8), ALLOCATABLE :: MMPOTEXPGLOB(:,:)
      REAL(real_8), ALLOCATABLE :: MMPOTEXPGLOB_VAR(:,:)
      REAL(real_8), ALLOCATABLE :: ATMULTEXPGLOB(:,:)
      REAL(real_8), ALLOCATABLE :: ATRGYREXPGLOB(:,:)
  
  
        
!....variables to handle voxels
      INTEGER :: MYNVOX
!     voxel coordinates r_lambda      
      REAL(real_8), ALLOCATABLE :: KOVOX(:,:)
!     voxel partitioning of grid
      INTEGER     , ALLOCATABLE :: VOXPART(:)
!     mapping voxels->qmatoms
      INTEGER     , ALLOCATABLE :: VO2IQ(:)
  


 
  END MODULE iffi_types
  
#if 0

!.....for passing the filename of the input file from iffi to cpmd
!X      CHARACTER*(200) INTF_FILENAME
!X      COMMON /EPOTFN/ INTF_FILENAME

!.....for printing the external potential (WRITE_EXTF)
      COMPLEX*16 EXTFCOMPLEX,MYVTEMP
      POINTER (IP_EXTFCOMPLEX,EXTFCOMPLEX) 
      POINTER (IP_MYVTEMP,MYVTEMP)
      COMMON /EPOTe/ IP_EXTFCOMPLEX,IP_MYVTEMP

!.....to send iteration steps needed for wavefunction iteration to iffi
      INTEGER ITSTEPS
      COMMON /ITSTEPS0/ ITSTEPS

!X.....General data: number of near atoms, number of local expansions (=qm atoms)
!X      INTEGER NNEAR,NLOC

!X.....Data for the local expansions
!X      INTEGER LOCID(MAXIND)
!X.....Number of near atoms to include
!X      INTEGER NINCL(MAXIND)
!X      INTEGER IDNEAR(MAXNEAR)
!.....List of near atoms to add to each local expansion
!X      INTEGER LINCL(MAXIND,MAXNEAR)

!X     Bounding box (BOXQM) and box (BOXDUM) edges
!X      REAL*8 BOXQM(6),BOXDUM(6)
!X.....Position of the Local Taylor Expansion 
!X      REAL*8 KOLOC(3,MAXIND)
!X.....Translation vector when box has been recentered by EGO
!X      REAL*8 BOXTRANS(3)
!X.....vector from iffi to cpmd coordinate frame (needed for output of potential/rhoe)
!X      REAL*8 BOXOFFS(3)



!.....Coefficients of the local expansion
!.....for readable access to LOCEXP entries
!X      INTEGER KPHI,KX,KY,KZ,KXX,KYY,KXY,KXZ,KYZ
!X      INTEGER KXXX,KYYY,KXYY,KYXX,KZXX,KZYY,KXYZ
!X      INTEGER KXXYY,KXXZZ,KYYZZ,KXXXY,KXXXZ,
!X     &        KYYYX,KYYYZ,KZZZX,KZZZY,DIMLOCEXP
!X      PARAMETER (KPHI=1)
!X      PARAMETER (KX=2,KY=3,KZ=4)
!X      PARAMETER (KXX=5,KYY=6,KXY=7,KXZ=8,KYZ=9)
!X      PARAMETER (KXXX=10,KYYY=11,KXYY=12,KYXX=13,
!X     &           KZXX=14,KZYY=15,KXYZ=16)
!X      PARAMETER (KXXYY=17,KXXZZ=18,KYYZZ=19,KXXXY=20,
!X     &           KXXXZ=21,KYYYX=22,KYYYZ=23,KZZZX=24,
!X     &           KZZZY=25,DIMLOCEXP=25)
!X      REAL*8 LOCEXP(DIMLOCEXP,MAXIND)

!.....Partial Charge, Core Charge, Sigma, CutOffWidth, CutOff
!X      INTEGER N_KX,N_KY,N_KZ,PCHARGE,SIGCHARGE, SIGDIPOLE,
!X     &        CORECHARGE, CO_WIDTH, CO_VALUE,SIZENEARDATA
!X      PARAMETER (N_KX=1,N_KY=2,N_KZ=3,PCHARGE=4,SIGCHARGE=5, 
!X     &           SIGDIPOLE=6, CORECHARGE=7, CO_WIDTH=8, CO_VALUE=9)
!X      PARAMETER(SIZENEARDATA=9)

!X      INTEGER   N_PX,N_PY,N_PZ,SIZENEARDATA_VAR
!X      PARAMETER (N_PX=1,N_PY=2,N_PZ=3)
!X      PARAMETER(SIZENEARDATA_VAR=3)

!X      REAL*8 NEARDATA    (SIZENEARDATA    ,MAXNEAR) 
!X      REAL*8 NEARDATA_VAR(SIZENEARDATA_VAR,MAXNEAR) 

!ZX      REAL*8 MYRAG(100)


!....logicals to steer what should be calculated
!X      LOGICAL TGRIDPART,TDOEXTRAP,TCALCONLYPOTFROMPFF,
!X     &        TNOPOTANDWFCTUPD,TCALCONLYPFFFIELD,TWRITERESTARTFILE,
!X     &        TUPDATEQMINTLISTS, TUPDATERGYR,
!X     &        UPDATEMEANFIELD

!X      INTEGER CALCESPCHARGES



!.....MULTIPOLE MOMENTS OF DFT ATOMS
!X      INTEGER DIMMULTEXP
!X      PARAMETER (DIMMULTEXP=25)
!.....for readable access to VOMULTEXP entries
!X      INTEGER LAD,PX,PY,PZ,QXX,QYY,QXY,QXZ,QYZ
!X      INTEGER OXYZ,OXXY,OXXZ,OYYX,OYYZ,OZZX,OZZY
!X      INTEGER HXXYY,HXXZZ,HYYZZ,HXXXY,HXXXZ,
!X     &        HYYYX,HYYYZ,HZZZX,HZZZY
!X      PARAMETER (LAD=1)
!X      PARAMETER (PX=2,PY=3,PZ=4)
!X      PARAMETER (QXX=5,QYY=6,QXY=7,QXZ=8,QYZ=9)
!X      PARAMETER (OXYZ=10,OXXY=11,OXXZ=12,OYYX=13,
!X     &           OYYZ=14,OZZX=15,OZZY=16)
!X      PARAMETER (HXXYY=17,HXXZZ=18,HYYZZ=19,HXXXY=20,
!X     &           HYYYX=21,HXXXZ=22,HZZZX=23,HYYYZ=24,
!X     &           HZZZY=25)


!.....POTENTIAL,FIELD at positions of QM atoms
!X      REAL*8  QMPOTANDFLD        (4,    MAXIND)
!X      POINTER (IP_QMPOTANDFLD,QMPOTANDFLD)
!X      REAL*8  ECHRG(MAXIND)
!X      POINTER (IP_ECHRG,ECHRG)

!.....POTENTIAL,FIELD,FIELDGRADIENT,PFFFIELD at positions of MM atoms
! ....first the fixed quantities..
!X      INTEGER DIMPOTEXP
!X      PARAMETER (DIMPOTEXP=10)
!.....for readable access to MMPOTEXP entries
!X      INTEGER  MPOT,MEX,MEY,MEZ,MEXX,MEYY,MEZZ,MEXY,MEXZ,MEYZ
!X      PARAMETER (MPOT=1,MEX=2,MEY=3,MEZ=4,MEXX=5,MEYY=6,
!X     &           MEZZ=7,MEXY=8,MEXZ=9,MEYZ=10)

! ....then the variable quantities..
!X      INTEGER DIMPOTEXP_VAR
!X      PARAMETER (DIMPOTEXP_VAR=3)
!.....for readable access to MMPOTEXP entries
!X      INTEGER  MEXPFF,MEYPFF,MEZPFF
!X      PARAMETER (MEXPFF=1,MEYPFF=2,MEZPFF=3)



!    MMPOTEXP_VAR  : QM field at NEAR atoms' dipoles (only quantity needed to update during pmm-scf)
!    MMPOTEXP      : rest of electrostatics at NEAR atoms
!    MMPOTEXPGLOB* : above quantities summed together at master
!X      REAL*8  MMPOTEXP        (DIMPOTEXP,    MAXNEAR,MXNO_THREADS),
!X     &        MMPOTEXP_VAR    (DIMPOTEXP_VAR,MAXNEAR,MXNO_THREADS),
!X     &        MMPOTEXPGLOB    (DIMPOTEXP,    MAXNEAR),
!X     &        MMPOTEXPGLOB_VAR(DIMPOTEXP_VAR,MAXNEAR)

!X      POINTER (IP_MMPOTEXP,MMPOTEXP),
!X     &        (IP_MMPOTEXP_VAR,MMPOTEXP_VAR),
!X     &        (IP_MMPOTEXPGLOB,MMPOTEXPGLOB),
!X     &        (IP_MMPOTEXPGLOB_VAR,MMPOTEXPGLOB_VAR)
     
!    VOMULTEXP     : multipole expansions of VOxels
!    ATMULTEXP     : multipole expansions of QM AToms
!    ATMULTEXPGLOB : above quantity summed together at master
!X      REAL*8   VOMULTEXP,
!X     &         ATMULTEXP(DIMMULTEXP,MAXIND,MXNO_THREADS),
!X     &         ATMULTEXPGLOB(DIMMULTEXP,MAXIND)

!X      POINTER (IP_VOMULTEXP,VOMULTEXP),
!X     &        (IP_ATMULTEXP,ATMULTEXP),
!X     &        (IP_ATMULTEXPGLOB,ATMULTEXPGLOB)
!

!CCCCCCCCCCCCCCCCCCCCCCCCCC
!     variables for gyration moments
!X      INTEGER RGZERO,RGX,RGY,RGZ
!X      INTEGER RGXX,RGYY,RGZZ,RGXY,RGXZ,RGYZ,NGAMMA,DIMRGYREXP
!X      PARAMETER (RGZERO=1,RGX=2,RGY=3,RGZ=4,
!X     &           RGXX=5,RGYY=6,RGZZ=7,RGXY=8,RGXZ=9,RGYZ=10,
!X     &           NGAMMA=11,DIMRGYREXP=11)

!X      REAL*8 ATRGYREXP(DIMRGYREXP,MAXIND,MXNO_THREADS),
!X     &       ATRGYREXPGLOB(DIMRGYREXP,MAXIND)
!X      POINTER (IP_ATRGYREXP,ATRGYREXP),(IP_ATRGYREXPGLOB,ATRGYREXPGLOB)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


!     Integer array to save MM potential on grid during pmm-scf
!X      REAL*8   EXTFMM
!X      POINTER(IP_EXTFMM,EXTFMM)
!X       COMMON /EXTFMM0/ IP_EXTFMM
        
!     for mean field method: 
!X       REAL*8   EXTFMEAN,EXTFLAST
!X       POINTER(IP_EXTFMEAN,EXTFMEAN),(IP_EXTFLAST,EXTFLAST)
!X       COMMON /EXTFMM1/ IP_EXTFMEAN,IP_EXTFLAST
!X      INTEGER NMEAN
        
!     FMM parameters
!X      REAL*8 OLERADIUS,OLEEPS,OLENSIG
!X      INTEGER SAMMTYPE
!X      INTEGER NRQMATOMS
!X      INTEGER VOXNUM(3)

!.....Common blocks             
!     note that in a parallel run all variables in these blocks have to be broadcasted to the
!     nodes by MY_BCAST in INTERFACE_READ (egointer.F).
!     For an integer vector of size N, add N*8/IRAT to the MSGLEN used for broadcasting the block.
!     For a floating point vector of size N, add N*8.

!.....quantities that do not change during the calculation
!     note: ints and logicals in one block!
!X!      COMMON /RUNINFO_int/ NRQMATOMS,SAMMTYPE,PFF,VOXNUM,
!X     &                     NORMALMODEANALYSIS,MEANFIELDMODE,TVERBOSE
!X!      COMMON /RUNINFO_dbl/ BOXQM,BOXDUM,OLERADIUS,OLEEPS,OLENSIG,MYRAG

!X!.....fixed quantities
!X!     interaction lists (changes each listupdate step)
!X      COMMON /EPOT1_fix_int/ NNEAR,NLOC,LOCID,NINCL,IDNEAR,LINCL     

!X     near data (changes each integration step)
!X      COMMON /EPOT2_fix_dbl/ KOLOC,NEARDATA,BOXOFFS,BOXTRANS

!....variable quantities (change each iteration step)
!X      COMMON /EPOT3_var_int/ TWRITERESTARTFILE,TGRIDPART,
!X     &                  TDOEXTRAP,TCALCONLYPOTFROMPFF,TNOPOTANDWFCTUPD,
!X     &                  TCALCONLYPFFFIELD,CALCESPCHARGES,
!X     &                  TUPDATEQMINTLISTS,TUPDATERGYR,
!X     &                  UPDATEMEANFIELD,NMEAN

!X      COMMON /EPOT4_var_dbl/ NEARDATA_VAR,LOCEXP

!X      COMMON /EPOT5/ IP_MMPOTEXP,IP_MMPOTEXP_VAR,IP_VOMULTEXP,
!X     &               IP_MMPOTEXPGLOB,IP_MMPOTEXPGLOB_VAR,
!X     &               IP_ATMULTEXPGLOB,IP_ATMULTEXP,
!X     &               IP_ATRGYREXP,IP_ATRGYREXPGLOB

!X      COMMON /QMPOTANDFLD1/ IP_QMPOTANDFLD,IP_ECHRG
!******************************************************************
!******************************************************************


!....variables to handle voxels
!X      REAL*8 KOVOX !coordinates
!X      POINTER(IP_KOVOX,KOVOX)

!X      REAL*8 VOXEXP,VOXEXP_MM   !local expansions at voxels
!X      POINTER(IP_VOXEXP,VOXEXP),(IP_VOXEXP_MM,VOXEXP_MM)

!X      REAL*8 VOXRGYREXP   !local expansions at voxels
!X      POINTER(IP_VOXRGYREXP,VOXRGYREXP)

!     voxel partitioning of grid
      !XINTEGER VOXPART                           
!X      POINTER(IP_VOXPART,VOXPART)

!X      !XINTEGER MYNVOX         !local voxel list, local number of voxels
!X      INTEGER VO2IQ               !mapping of voxels to qmatoms
!X      POINTER(IP_VO2IQ,VO2IQ)

              
!     near and ole lists of voxels
!X      INTEGER NEARLIST,NEARLISTLEN
!X      POINTER(IP_NEARLIST,NEARLIST),(IP_NEARLISTLEN,NEARLISTLEN)
!X      INTEGER OLELIST,OLELISTLEN
!X      POINTER(IP_OLELIST,OLELIST),(IP_OLELISTLEN,OLELISTLEN)

!X      COMMON /VOXPART1/ IP_KOVOX
!X      COMMON /VOXPART2/ IP_VOXPART,IP_VOXEXP,IP_VOXEXP_MM,
!X     &                  IP_VO2IQ,IP_VOXRGYREXP
!X      COMMON /VOXPART3/ MYNVOX
!X      COMMON /VOXPART4/ IP_NEARLIST,IP_NEARLISTLEN,
!X     &                  IP_OLELIST,IP_OLELISTLEN

#endif
