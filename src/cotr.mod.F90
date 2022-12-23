MODULE cotr
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==================================================================
  ! == CONSTRAINTS OF GEOMETRY FOR NUCLEAR POSITION                 ==
  ! ==--------------------------------------------------------------==
  ! == DMAX   Maximum sum of total possible displacements           ==
  ! == MAXDUM Maximum dummy atoms                                   ==
  ! == MAXLI  Maximum number of atoms given the coordinates of      ==
  ! ==        dummy atom type 2.                                    ==
  ! ==--------------------------------------------------------------==
  REAL(real_8), PARAMETER :: dmax=0.3_real_8 
  INTEGER, PARAMETER :: maxdum=50,maxli=20 
  ! ==--------------------------------------------------------------==
  ! == LFCOM  .TRUE. if fix the center of mass                      ==
  ! == LSHOVE .TRUE. if SHOVE option                                ==
  ! == NODIM  Number of freedom degres                              ==
  ! == NBOAD  Number of bonds which has to be changed               ==
  ! ==        (added or deleted from the empirical hessian          ==
  ! == MCNSTR Number of constraints                                 ==
  ! == BOAD(3,NBOAD) BOAD(1,I) and BOAD(2,I) indicate the bond      ==
  ! ==        between atoms. BOAD(3,I)<0 deleted, >0 added          ==
  ! == LSKCOR(3,NAT) 0 if the atomic coordinate is fixed            ==
  ! == NTCNST(6,MCNSTR) NTCNST(1,:) gives the type                  ==
  ! ==                    1 STRETCH                                 ==
  ! ==                    2 BEND                                    ==
  ! ==                    3 TORSION                                 ==
  ! ==                    4 DISTance           between atoms        ==
  ! ==                      NTCNST(2,:) and NTCNST(3,:)             ==
  ! ==                    5 OUTP                                    ==
  ! ==                    6 COORD                                   ==
  ! ==                    7 DIFFER                                  ==
  ! ==                    8 COORSP                                  ==
  ! ==                    9 COOR_RF                                 ==
  ! ==                   10 BNSWT                                   ==
  ! ==                   11 TOT_COOR                                ==
  ! ==                   12 DISAXIS                                 ==
  ! ==                   13 RESPOS                                  ==
  ! ==                  101 MMFIELD                                 ==
  ! ==                  NTCNST(2:5,:) gives the indexes of atoms    ==
  ! ==                  if needed otherwise = 0                     ==
  ! ==                  NTCNST(6,:) used if LSHOVE                  ==
  ! ==                                                              ==
  ! == LSKPTR(3,NAT)  0 if the atomic coordinate is fixed           ==
  ! ==                  otherwise index between 1 and NODIM         ==
  ! ==--------------------------------------------------------------==
  ! == HESS(NODIM,NODIM) Hessian matrix (second derivatives)        ==
  ! == DTM(NODIM)        Step for steepest descent                  ==
  ! == PMALL(NODIM)      Atomic masses if not fixed otherwise 0.    ==
  ! == KSHOVE(MCNSTR)                                               ==
  ! == GRATE(MCNSTR)     Growth rate for constraint                 ==
  ! == CNSVAL(MCNSTR)    Value of constraint                        ==
  ! == CNPAR(2,MCNSTR)   Parameters for coord. num or bond switch   ==
  ! == ANORM(NODIM,MCNSTR)                                          ==
  ! == ASKEL(NODIM,MCNSTR,12)                                       ==
  ! == mm_ASKEL(MCNSTR,12)                                    ==
  ! == FCSTR(NODIM)      Force ConSTRaint                           ==
  ! == FC(NODIM)                                                    ==
  ! == FV(NODIM)                                                    ==
  ! == CSIGM(MCNSTR) The weight factor for penalty function         ==
  ! == XLAGR(MCNSTR)                                                ==
  ! == YLAGR(MCNSTR)                                                ==
  ! ==--------------------------------------------------------------==

  INTEGER, ALLOCATABLE, SAVE :: boad(:,:)
  INTEGER, ALLOCATABLE, SAVE :: lskcor(:,:)
  INTEGER, ALLOCATABLE, SAVE :: ntcnst(:,:)
  INTEGER, ALLOCATABLE, SAVE :: lskptr(:,:)
  INTEGER, ALLOCATABLE, SAVE :: mm_ASKEL(:,:)

  REAL(real_8), ALLOCATABLE, SAVE :: hess(:,:)
  REAL(real_8), ALLOCATABLE, SAVE :: dtm(:)
  REAL(real_8), ALLOCATABLE, SAVE :: pmall(:)
  REAL(real_8), ALLOCATABLE, SAVE :: kshove(:)
  REAL(real_8), ALLOCATABLE, SAVE :: cnsval(:)
  REAL(real_8), ALLOCATABLE, SAVE :: cnpar(:,:)
  REAL(real_8), ALLOCATABLE, SAVE :: anorm(:,:)
  REAL(real_8), ALLOCATABLE, SAVE :: askel(:,:,:)
  REAL(real_8), ALLOCATABLE, SAVE :: grate(:)
  REAL(real_8), ALLOCATABLE, SAVE :: cnsval_dest(:)
  REAL(real_8), ALLOCATABLE, SAVE :: fcstr(:)
  REAL(real_8), ALLOCATABLE, SAVE :: fc(:)
  REAL(real_8), ALLOCATABLE, SAVE :: fv(:)
  REAL(real_8), ALLOCATABLE, SAVE :: csigm(:)
  REAL(real_8), ALLOCATABLE, SAVE :: xlagr(:)
  REAL(real_8), ALLOCATABLE, SAVE :: ylagr(:)


  ! ==--------------------------------------------------------------==
  ! == RESTRAINTS                                                   ==
  ! ==--------------------------------------------------------------==
  INTEGER, ALLOCATABLE, SAVE :: ntrest(:,:)


  REAL(real_8), ALLOCATABLE, SAVE :: resval(:)
  REAL(real_8), ALLOCATABLE, SAVE :: resfor(:)
  REAL(real_8), ALLOCATABLE, SAVE :: gsrate(:)
  REAL(real_8), ALLOCATABLE, SAVE :: resval_dest(:)
  REAL(real_8), ALLOCATABLE, SAVE :: rskel(:,:,:)
  REAL(real_8), ALLOCATABLE, SAVE :: resc(:)
  REAL(real_8), ALLOCATABLE, SAVE :: resv(:)
  REAL(real_8), ALLOCATABLE, SAVE :: respos(:,:)
  REAL(real_8), ALLOCATABLE, SAVE :: resfdiff(:)
  REAL(real_8), ALLOCATABLE, SAVE :: resm(:,:)
  REAL(real_8), ALLOCATABLE, SAVE :: rsmass(:)
  REAL(real_8), ALLOCATABLE, SAVE :: respar(:,:)

  ! mb   == RESPOS(3,MRESTR) contains the cartesian coordinates of the   ==
  ! mb   ==                  restraint position                          ==

  ! ==--------------------------------------------------------------==
  ! == PATOT    Total atomic masses if not fixed (SUM of PMALL(I))  ==
  ! == The weight factors for the penalty function for:             ==
  ! == DSIGMA    stretchs,                                          ==
  ! == BSIGMA    bends,                                             ==
  ! == TSIGMA    torsions.                                          ==
  ! ==--------------------------------------------------------------==
  REAL(real_8) :: patot,dsigma,bsigma,tsigma
  ! ==--------------------------------------------------------------==
  ! == COORDINATION NUMBER CONSTRAINTS                              ==
  ! == C_KAPPA : width of Fermi function                            ==
  ! == C_RC    : radius of Fermi function                           ==
  ! ==--------------------------------------------------------------==
  REAL(real_8), SAVE :: c_kappa,c_rc
  ! ==--------------------------------------------------------------==
  ! == DUMMY ATOMS : 3 types                                        ==
  ! ==   Type 1 is fixed in place                                   ==
  ! ==   Type 2 lies at the arithmetic mean of the coordinates of   ==
  ! ==          real atoms.                                         ==
  ! ==   Type 3 lies at the center of mass of the coordinates of    ==
  ! ==          real atoms.                                         ==
  ! == DUMMY1(3,MAXDUM) coordinates of type 1 dummy atoms           ==
  ! == DUMMY2(3,MAXDUM) coordinates of type 2 dummy atoms           ==
  ! == DUMMY3(3,MAXDUM) coordinates of type 3 dummy atoms           ==
  ! == LISTDA(MAXDUM,1) Type of the dummy atom                      ==
  ! == LISTDA(MAXDUM,2) Index of the dummy atom                     ==
  ! == LISTD2(MAXLI,MAXDUM) For each dummy atom type 2, give the    ==
  ! ==        indexes of atoms (type 2 = arithmetic mean of coord.) ==
  ! == LISTD3(MAXLI,MAXDUM) For each dummy atom type 3, give the    ==
  ! ==        indexes of atoms (type 3 = center of mass of coord.)  ==
  ! == NDAT   Total number of dummy atoms                           ==
  ! == NDAT1  Number of dummy atoms type 1                          ==
  ! == NDAT2  Number of dummy atoms type 2                          ==
  ! == NDAT3  Number of dummy atoms type 3                          ==
  ! ==--------------------------------------------------------------==


  TYPE :: duat_t
     REAL(real_8) :: dummy1(3,maxdum)
     REAL(real_8) :: dummy2(3,maxdum)
     REAL(real_8) :: dummy3(3,maxdum)
     REAL(real_8) :: dummy4(3,maxdum)
     REAL(real_8) :: weigd4(maxli,maxdum)
     INTEGER :: listda(maxdum,2)
     INTEGER :: listd2(maxli,maxdum)
     INTEGER :: listd3(maxli,maxdum)
     INTEGER :: listd4(maxli,maxdum)
     INTEGER :: ndat
     INTEGER :: ndat1
     INTEGER :: ndat2
     INTEGER :: ndat3
     INTEGER :: ndat4
  END TYPE duat_t
  TYPE(duat_t), SAVE :: duat
  ! ==--------------------------------------------------------------==
  ! == TSDPLINE Minimization along the line of steepest descent     ==
  ! ==          if TSDP=.TRUE.                                      ==
  ! == TINVBFGS Use inverse matrix hessian for BFGS                 ==
  ! == TCGP     .TRUE. Conjugate gradient scheme                    ==
  ! == PXPAR(NODIM) Direction along the minimization                ==
  ! ==--------------------------------------------------------------==

  REAL(real_8), ALLOCATABLE, SAVE :: pxpar(:)

  ! ==================================================================

  TYPE :: cotc0_t
     LOGICAL :: lfcom
     LOGICAL :: lshove
     INTEGER :: nodim
     INTEGER :: nboad
     INTEGER :: mcnstr
     INTEGER :: nctim
  END TYPE cotc0_t
  TYPE(cotc0_t), SAVE :: cotc0
  TYPE :: cotr007_t
     INTEGER :: mrestr
     LOGICAL :: lhyperplane
  END TYPE cotr007_t
  TYPE(cotr007_t), SAVE :: cotr007
  TYPE :: sdpl_t
     LOGICAL :: tsdpline
     LOGICAL :: tcgp
     LOGICAL :: tinvbfgs
  END TYPE sdpl_t
  TYPE(sdpl_t), SAVE :: sdpl

END MODULE cotr
