C $Id: pmd.h,v 1.1 2006-12-27 11:29:48 itavern Exp $ -*-fortran-*-
C

C
C CONSTANTS AND VARIABLES TO CONTROL BEHAVIOUR OF PRE MD PROGRAMS

C
C MAX. LENGTH OF TITLE
C MXPTIT SHOULD HAVE THE SAME VALUE AS MAXLNS (toposz.h)
      INTEGER MXPTIT, MDTILE
      PARAMETER (MXPTIT = 10, MDTILE = 80)

      LOGICAL LMDOK
      COMMON /MDDDD/LMDOK

C COMMON DIMENSION
      INTEGER NDIM
      PARAMETER (NDIM = 3)

C LENGTH OF PROCS1 ARRAYS
      INTEGER IDIM, MAMAX
      PARAMETER (IDIM = 2, MAMAX = 26)

C LENGTH OF PROCS2 ARRAYS
      INTEGER MAXREA
      PARAMETER (MAXREA = 5000)

C LENGTH OF PROGCA ARRAYS
      INTEGER MAXNRA, NAT3X
      PARAMETER (MAXNRA = 100, NAT3X = 3)

C LENGTH OF PROGCH ARRAYS
      INTEGER MAXIAT, MHYDIM
      PARAMETER (MAXIAT = 4, MHYDIM = 4)

C LENGTH OF PROSSC ARRAYS
      INTEGER MAXNSR
      PARAMETER (MAXNSR = 70)

C LENGTH OF PROGMT ARRAYS
C SEE INLCUDE 'GMT.H'

C LENGTH OF PROCRY ARRAYS
      INTEGER MAXSYM
      PARAMETER (MAXSYM = 32)

C LENGTH OF PROION ARRAYS
      INTEGER MAXEXS
      PARAMETER (MAXEXS = 20)

C LENGTH OF PROCMT ARRAYS
      INTEGER MAXRES
      PARAMETER (MAXRES = 25)

C----------------------------------------------------
C COMMON BLOCKS
C SEPARATE COMMON BLOCKS FOR INTS/LOGICALS, REALS
C AND CHARS

C CHAR BLOCK 
      CHARACTER*(MDTILE) MDTITL(MXPTIT)
      CHARACTER*(5) APLUS, AMIN, 
     $     AMOL1, AMOL2
      COMMON /PRECAR/ MDTITL, APLUS, AMIN,
     $     AMOL1, AMOL2

C INTEGER BLOCK
      INTEGER
     $     NLRED, 
     $     NTM, NTF, MMAX, NREC,
     $     NREA, NTS, NTX, NTBF, NRE,
     $     NTA, NRA, NR, NTXI, NTXO, NTRA, IG,
     $     IA, JA, KA, LA, MNB,
     $     NTH, NTU, NIAT, IACHT,
     $     NRSC, JRES,
     $     NTXP, NPM, NSMP, NTXS, NSM, NTB,
     $     NRATO,
     $     NSYM,
     $     NTR, NEXSM, JSNSOL, NPLUS, NMIN,
     $     NMOL2,
     $     NSM2, NTLIS, NRESI, IRES, NRDBOX 


      COMMON /PREIN/ 
     $     NLRED,
     $     NTM, NTF, MMAX, NREC(IDIM, MAMAX),
     $     NREA, NTS, NTX, NTBF, NRE(MAXREA),
     $     NTA, NRA, NR, NTXI, NTXO, NTRA, IG,
     $     IA(MAXNRA), JA(MAXNRA), KA(MAXNRA), LA(MAXNRA), MNB(MAXNRA),
     $     NTH, NTU, NIAT, IACHT(MAXIAT, MHYDIM),
     $     NRSC, JRES(MAXNSR),
     $     NTXP, NPM, NSMP, NTXS, NSM, NTB,
     $     NRATO,
     $     NSYM,
     $     NTR, NEXSM, JSNSOL(MAXEXS), NPLUS, NMIN,
     $     NMOL2,
     $     NSM2, NTLIS, NRESI, IRES(MAXRES), NRDBOX 


C REAL BLOCK
      REAL SCALX, SCALI, SCALB, BOX, BETA,
     $     X3, PHI, ALFA, BOND, PHIMIN, PHIMAX,
     $     BETAA, BOXS, DISM, XMIN,
     $     ABC, SYM,
     $     RCUTE, CGPLUS, CGMIN,
     $     XREF, RCUT1, RCUT2

      COMMON /PRERE/ 
     $     SCALX, SCALI, SCALB, BOX(NDIM), BETA(NDIM),
     $     X3(NDIM * NAT3X),
     $     PHI(MAXNRA), ALFA(MAXNRA), BOND(MAXNRA), 
     $     PHIMIN(MAXNRA), PHIMAX(MAXNRA),
     $     BETAA, BOXS, DISM, XMIN(NDIM),
     $     ABC(NDIM), SYM(3, 4, MAXSYM),
     $     RCUTE, CGPLUS, CGMIN,
     $     XREF(NDIM), RCUT1, RCUT2


C-----------END OF COMMON BLOCKS --------------------
C     HERE SOME USEFUL CONSTANTS
C ONE
      INTEGER IONE
      PARAMETER (IONE = 1)

C MASS CODE OF A HYDROGEN (ALSO IN gmt.h)
      INTEGER ISHYD
      PARAMETER (ISHYD = 1)

C A VERY SMALL NUMBER
      REAL EPS
      PARAMETER (EPS = 1.E-6)

C A VERY LARGE NUMBER
      REAL RELBIG
      PARAMETER (RELBIG = 10E24)


C BOLTZMAN CONSTANT FOR KJ/MOL
      REAL BOLKJM
      PARAMETER (BOLKJM = 8.31441E-3)

C BOLTZMAN CONSTANT FOR KCAL/MOL
      REAL BOLKCM
      PARAMETER (BOLKCM = BOLKJM/4.184E0)



