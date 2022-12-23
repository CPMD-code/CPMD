MODULE symm
  USE kinds,                           ONLY: real_8
  USE system,                          ONLY: maxsp

  IMPLICIT NONE

  ! ==================================================================
  ! == INFORMATION ABOUT SYMMETRY OF THE SYSTEM                     ==
  ! ==================================================================
  ! == TSYMM   SYMMETRIZED COORDINATES                              ==
  ! == TGENC   GENERATE ATOMIC COORDINATES FROM SYMM. UNIQUE ONES   ==
  ! == TMSYM   POINT GROUP MOLECULE                                 ==
  ! == TPGAUTO POINT GROUP AUTO                                     ==
  ! == TORIGIN =.TRUE. IF SYMMORPHIC AND ORIGIN=(0,0,0)(SYM. CENTER)==
  ! == TSYMMORPHIC =.TRUE. IF SYMMORPHIC GROUP (TPGAUTO=.TRUE.)     ==
  ! ==                     OTHERWISE IF ALL FTAU=0._real_8               ==
  ! == TSYMMOR(1:NROT)=.TRUE. IF FTAU(IR)=0._real_8                      ==
  ! == [Broadcasted in SETSYS routine]                              ==
  ! ==--------------------------------------------------------------==

  ! ==--------------------------------------------------------------==
  ! == INDPG   POINT GROUP NUMBER                                   ==
  ! == IHG     POINT GROUP OF THE PRIMITIVE LATTICE, HOLOHEDRAL     ==
  ! ==         GROUP NUMBER:                                        ==
  ! ==          IHG=1 STANDS FOR TRICLINIC    SYSTEM                ==
  ! ==          IHG=2 STANDS FOR MONOCLINIC   SYSTEM                ==
  ! ==          IHG=3 STANDS FOR ORTHORHOMBIC SYSTEM                ==
  ! ==          IHG=4 STANDS FOR TETRAGONAL   SYSTEM                ==
  ! ==          IHG=5 STANDS FOR CUBIC        SYSTEM                ==
  ! ==          IHG=6 STANDS FOR TRIGONAL     SYSTEM                ==
  ! ==          IHG=7 STANDS FOR HEXAGONAL    SYSTEM                ==
  ! == NROT    NUMBER OF ROTATIONS                                  ==
  ! == INVERSION INDEX OF THE INVERSION ROTATION OTHERWISE 0        ==
  ! == NTVEC   NUMBER OF TRANSLATION VECTORS ASSOCIATED WITH 1      ==
  ! == IRTVEC(NAT,NTVEC) ATOM TRANSFORMATION TABLE FOR TRANSLATION  ==
  ! ==         VECTORS                                              ==
  ! == INVE(120) NUMBER OF INVERSE ROTATIONS (GIVEN BY MULTTB)      ==
  ! ==           (USED BY SYMRHO)                                   ==
  ! == MULTAB(120,120) MULTIPLICATION TABLE                         ==
  ! == IRT(120,NAT) ATOM TRANSFORMATION TABLE                       ==
  ! ==         FOR SOLID THE 49TH ROW GIVES INEQUIVALENT ATOMS FOR  ==
  ! ==         TRANSLATIONS (IF NTVEC /= 1)                         ==
  ! == ISYMU(NAT)   NUMBER OF EQUIVALENT ATOMS                      ==
  ! == [Broadcasted in SETSYS routine]                              ==
  ! ==--------------------------------------------------------------==
  INTEGER, ALLOCATABLE :: irt(:,:)
  INTEGER, ALLOCATABLE :: isymu(:)
  INTEGER, ALLOCATABLE :: irtvec(:,:)


  ! ==--------------------------------------------------------------==
  ! == IUN(MAXSP) USED WITH GENERATE COORDINATES OPTION             ==
  ! == STABLE(3,3,48) ROTATION MATRICES (INTEGER)                   ==
  ! == [UnBroadcasted]                                              ==
  ! ==--------------------------------------------------------------==
  INTEGER :: iun(maxsp),stable(3,3,48)
  ! ==--------------------------------------------------------------==
  ! == DELTASYM        REQUIRED ACCURACY (1.e-6_real_8 IS A GOOD VALUE)    ==
  ! ==                 IN DIRECT VECTOR LATTICE BASIS SET           ==
  ! ==                 (SO IN CARTESIAN = ALAT*DELTASYM)            ==
  ! == ORIGIN(3)       IF SYMMORPHIC CENTER OF SYMMETRY             ==
  ! ==                 OTHERWISE (0,0,0)                            ==
  ! == XTABLE(3,3,120) ROTATION MATRICES IN REAL SPACE              ==
  ! ==                 IN CRYSTAL COORDINATES                       ==
  ! == FTABLE(3,3,120) ROTATION MATRICES IN FOURIER SPACE           ==
  ! == FTAU(3,120)     TRANSLATIONS ASSOCIATED WITH XTABLE          ==
  ! == [Broadcasted in SETSYS routine]                              ==
  ! ==--------------------------------------------------------------==
  ! ==--------------------------------------------------------------==
  ! == TVEC(3,NTVEC)   TRANSLATION VECTORS ASSOCIATED WITH 1        ==
  ! == [Broadcasted in SETSYS routine]                              ==
  ! ==--------------------------------------------------------------==
  REAL(real_8), ALLOCATABLE :: tvec(:,:)

  ! ==--------------------------------------------------------------==
  ! == STAG   LABEL OF THE POINT GROUP                              ==
  ! == [UnBroadcasted]                                              ==
  ! ==--------------------------------------------------------------==
  CHARACTER (len=3) :: stag
  ! ==================================================================
  ! == Use in serial/parallel version (SYMRHO)                      ==
  ! == NUMSHEL   Number of shells of equivalent G-vectors           ==
  ! ==           equivalent means symmetry links all G-vectors      ==
  ! == IGEXTINC  Number of G-vectors where RHO(IG)=0 (extinction)   ==
  ! == INDSHEL(NHG) For each G-vector index of the equivalent       ==
  ! ==          shell, i.e. 2 G-vectors related by symmetry have    ==
  ! ==          the same index.                                     ==
  ! ==          If NTVEC/=1 (translation vectors /=0 exist)         ==
  ! ==          0 if RHO(G) has to be nullify                       ==
  ! == ISYSHEL(2,NHG)                                               ==
  ! ==          If TSYMMORPHIC                                      ==
  ! ==          give for IG  if ISYSHEL>0                           ==
  ! ==          or   for -IG if ISYSHEL<0 the multiplicative factor ==
  ! == NINSHEL(NHG) !!!WARNING INTEGER*8                            ==
  ! ==          If TSYMMORPHIC unused                               ==
  ! ==          If .NOT.TSYMMORPHIC                                 ==
  ! ==          For equivalent shell,                               ==
  ! ==          give for each G-vector, the related rotations.      ==
  ! ==          For the rotations IR, the IRth bit sets to 1        ==
  ! == ISISHEL(NHG) !!!WARNING INTEGER*8                            ==
  ! ==          If TSYMMORPHIC unused                               ==
  ! ==          If .NOT.TSYMMORPHIC                                 ==
  ! ==          For equivalent shell, give for each G-vector        ==
  ! ==          if IG or -IG has to take. The IRth bit sets to 1    ==
  ! ==          the component IG is taken for the IRth rotation     ==
  ! ==--------------------------------------------------------------==
  ! == WARNING: you need to use INTEGER*8 in order ot have 64 bits  ==
  ! ==--------------------------------------------------------------==
  INTEGER, ALLOCATABLE :: indshel(:)
  INTEGER, ALLOCATABLE :: isyshel(:,:)
  INTEGER :: numshel,igextinc

  ! ==================================================================
  ! == SCALE OPTION (use in sysin and setsys -- no broadcast)       ==
  ! == SXSCALE Scale factor for x-axis                              ==
  ! == SYSCALE Scale factor for y-axis                              ==
  ! == SZSCALE Scale factor for z-axis                              ==
  ! == TCARTESIAN Coordinates are given in cartesian coordinates    ==
  ! ==            (used only with SCALE option)                     ==
  ! ==--------------------------------------------------------------==
  REAL(real_8) :: sxscale,syscale,szscale
  LOGICAL :: tcartesian
  ! ==================================================================

  TYPE :: symmi_t
     INTEGER :: indpg
     INTEGER :: ihg
     INTEGER :: nrot
     INTEGER :: naxis
     INTEGER :: iunique
     INTEGER :: inversion
     INTEGER :: ntvec
     INTEGER :: inve(120)
     INTEGER :: multab(120,120)
  END TYPE symmi_t
  TYPE(symmi_t) :: symmi
  TYPE :: symmr_t
     REAL(real_8) :: xtable(3,3,120)
     REAL(real_8) :: ftable(3,3,120)
     REAL(real_8) :: ftau(3,120)
     REAL(real_8) :: origin(3)
     REAL(real_8) :: deltasym
  END TYPE symmr_t
  TYPE(symmr_t) :: symmr
  TYPE :: symmt_t
     LOGICAL :: tsymm
     LOGICAL :: tgenc
     LOGICAL :: tmsym
     LOGICAL :: tpgauto
     LOGICAL :: torigin
     LOGICAL :: tsymmorphic
     LOGICAL :: tsymmor(120)
  END TYPE symmt_t
  TYPE(symmt_t) :: symmt

END MODULE symm
