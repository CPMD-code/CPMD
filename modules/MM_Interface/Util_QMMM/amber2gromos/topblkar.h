C$Id: topblkar.h,v 1.1 2006-12-27 11:32:40 itavern Exp $  -*-fortran-*-

C arrays used to hold molecular topology building block data

C--------------------------------------------------------------
C remember to include topblksz.h and toposz.h before this one.
C--------------------------------------------------------------


COMMVAR NRNE
C     this value is read from the linkexclusions block at the
C     beginning of the topology building block file.
      INTEGER NRNE
      COMMON /LNKINT/NRNE
COMMEND


COMMVAR NMRC,NLIN,MB,MCBL,MBAT,MCBAT,MIDAT,MCIAT,MDAT,MCDAT,MPD,MPA,IACM,IMCM
C
C*****solute topology building blocks are read into these arrays
C     by subroutine L<RDSOLU>.
C
      INTEGER NMRC,NLIN
      INTEGER NAT,NB,NBA,NIDA,NDA,NPD,NPA,NTBL,NTBA,NTIA,NTDA
      INTEGER MB(2,MAXNB),MCBL(MAXNB)
      INTEGER MBAT(3,MAXNBA),MCBAT(MAXNBA)
      INTEGER MIDAT(4,MAXIDA),MCIAT(MAXIDA)
      INTEGER MDAT(4,MAXNDA), MCDAT(MAXNDA)
      INTEGER MPD(MAXPHD)
      INTEGER MPA(MAXPHA)

      INTEGER IACM(MAXNAT),IMCM(MAXNAT),ICGM(MAXNAT),MAE(MAXNAT)
      INTEGER MSAE(MAXNAE,MAXNAT)
C LISHYD(I) = 'atom I is a hydrogen'
      LOGICAL LISHYD(MAXNAT)

      COMMON /TOPBLA/NMRC,NAT,NB,NBA,NIDA,NDA,NPD,NPA,NTBL,NTBA,
     $     NTIA,NTDA,NLIN,
     $     MB,MCBL,
     $     MBAT,MCBAT,
     $     MIDAT,MCIAT,
     $     MDAT,MCDAT,
     $     MPD,MPA,
     $     IACM,IMCM,ICGM,MAE,MSAE,
     $     LISHYD

      REAL CGM(MAXNAT)
      COMMON /TOPREL/CGM

C     RNME: the name of the solute building block
      CHARACTER RNME*(MAXRLE)

C     ANM: the atom names of the solute building block
      CHARACTER ANM*(MAXTLE)
      COMMON /TOPCHA/ RNME,ANM(MAXNAT)
COMMEND


COMMVAR IACMS,IMCMS,IJCONM,CGMS,CONM,WASS
C*****solvent build blocks are read into these arrays
C     by subroutine L<RDSOLV>.
C

C     NMATS: the number of atoms in the solvent building block
C     NCONM: the number of constraints in the solvent building block
C     
      INTEGER NMATS,IACMS(MAXNAT),IMCMS(MAXNAT)
      INTEGER NCONM,IJCONM(2,MAXCNS)

      REAL    CGMS(MAXNAT),CONM(MAXCNS),WASS(MAXATT)

      COMMON /TPSLVI/NMATS,IACMS,IMCMS,NCONM,IJCONM
      COMMON /TPSLVR/CGMS,CONM,WASS

C     RNMES: the name of the solvent building block
      CHARACTER RNMES*(MAXRLE)

C     ANNMS: the name of the solvent atoms
      CHARACTER ANMMS*(MAXTLE)
      COMMON /TOPCHA/ RNMES,ANMMS(MAXNAT)
COMMEND
C


