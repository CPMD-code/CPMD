C$Id: topblksz.h,v 1.1 2006-12-27 11:33:07 itavern Exp $  -*-fortran-*-

C define the sizes of molecular topology building blocks
C

C     maximum length of the title
      INTEGER NTITLE
      PARAMETER (NTITLE = 80)

COMMVAR MAXNAT,MAXAT1
C     maximum number of atoms per building block
      INTEGER MAXNAT,MAXAT1
      PARAMETER (MAXNAT = 100,MAXAT1 = MAXNAT+1)
COMMEND

COMMVAR MAXNB
C     maximum number of bonds per building block
      INTEGER MAXNB
      PARAMETER (MAXNB = 100)
COMMEND

COMMVAR MAXNBA
C     maximum number of bond angles per building block
      INTEGER MAXNBA
      PARAMETER (MAXNBA = 100)
COMMEND

COMMVAR MAXIDA
C     maximum number of improper dihedrals per building block
      INTEGER MAXIDA
      PARAMETER (MAXIDA = 50)
COMMEND

COMMVAR MAXNDA
C     maximum number of dihedrals per building block
      INTEGER MAXNDA
      PARAMETER (MAXNDA = 100)
COMMEND

COMMVAR MAXPHD
C     max number of donors per building block
      INTEGER MAXPHD
      PARAMETER (MAXPHD = 6)
COMMEND

COMMVAR MAXPHA
C max number of acceptors per building block
      INTEGER MAXPHA
      PARAMETER (MAXPHA = 6)
COMMEND

COMMVAR MAXNAE
C     max number of exclusions per building block
      INTEGER MAXNAE
      PARAMETER (MAXNAE = 25)
COMMEND

COMMVAR MAXCNS
C     maximum number of constraints in a solvent building blocks
      INTEGER MAXCNS
      PARAMETER (MAXCNS = 30)
COMMEND

