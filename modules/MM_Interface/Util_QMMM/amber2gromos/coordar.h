C $Id: coordar.h,v 1.1 2006-12-27 11:18:21 itavern Exp $ -*-fortran-*-
C GROMOS common blocks for coordinates, restraining and local elevation

C -----------------------------------------------------------
C always include the file coordsz.h before including this one
C -----------------------------------------------------------

C     The variables defining the extent of the
C     arrays used (e.g. NRCON) in these common blocks
C     are initialized to 0 in a BLOCKDATA statement in posio.f
C     These limits are then modified on reading data from
C     file.

COMMVAR X,V,F,XR,SX,XSAVE
C     Coordinate arrays, their use depends on whether
C     Molecular Dynamics, Stochastic Dynamics or Energy Minimisation
C     is performed.
C     X       atom coordinates                 MD,SD,EM
C     F       atom forces                      MD,SD,EM
C     V       atom velocities in MD,SD, scratch array in EM
C     XR      for virial calculation in MD,SD, scratch array in EM
C     SX      stochastic integral array        SD
C     XSAVE   atom coordinates for saving      
C             lowest energy configurations     MD,SD
      REAL X,V,F,XR,SX,XSAVE
      COMMON /COORD/ X(MAXXCO),F(MAXXCO),V(MAXXCO),
     $     XR(MAXXCO),SX(MAXXCO),XSAVE(MAXXCO)
COMMEND

COMMVAR CORTIT
C     CORTIT: a title read from the coordinate file
      CHARACTER *(MXCOTI) CORTIT
      COMMON /CRDTIT/ CORTIT
COMMEND

COMMVAR IAGRP
C     IAGRP is set in L<PROMD> using the values specified in the
C     array L<NRE>.
C     IAGRP is used for a faster lookup of energy contributions in
C     the force calculation routines.
C     In order to save space, however, it is defined as a character
C     array instead of an integer array.
C     The values are written into the array using CHAR() and read
C     using ICHAR(), both of which are standard F77 functions.

      CHARACTER IAGRP
      COMMON /GROUPY/IAGRP(MAXNAT)
COMMEND

COMMVAR JRC,NRCON,NRRST
C*****position re(con)straining
C     We can have EITHER restraining OR constraining.
C     In the event of restraining : NRCON = 0 and NRRST > 0
C     In the event of constraining: NRCON > 0 and NRRST = 0
C
C     In both cases the array JRC holds the indices of re(con)strained
C     atoms. However, the allowed range is different:
C     for constraining JRC(I) 1..L<NRP>
C     for restraining  JRC(I) 1..L<NATTOT>
C
      INTEGER JRC,NRCON,NRRST
      COMMON /POSINT/ JRC(MAXXC),NRCON,NRRST
COMMEND

COMMVAR XC,CXC,XC0,EC0
C*****position re(con)restraining
C     XC  : reference positions
C     CXC : force constants used in the position restraining
C     
      REAL XC,CXC
      COMMON /POSR/ 
     $     XC(MAXXCO),CXC(MAXNAT)
COMMEND


COMMVAR IDR1,JDR1,KDR1,LDR1,ICDR1,IDR2,JDR2,KDR2,LDR2,ICDR2,NDR
C*****atom distance restraining
      INTEGER IDR1,JDR1,KDR1,LDR1,ICDR1
      INTEGER IDR2,JDR2,KDR2,LDR2,ICDR2,NDR
      
      COMMON /DISRIN/
     $     IDR1(MAXNDR),JDR1(MAXNDR),KDR1(MAXNDR),LDR1(MAXNDR),
     $     ICDR1(MAXNDR),
     $     IDR2(MAXNDR),JDR2(MAXNDR),KDR2(MAXNDR),LDR2(MAXNDR),
     $     ICDR2(MAXNDR),
     $     NDR

COMMEND

COMMVAR R0,W0,RIIAVE,DISH,DISC
      REAL R0,W0,RIIAVE,DISH,DISC
      COMMON /DISRFP/ DISH,DISC,R0(MAXNDR),W0(MAXNDR),
     $     RIIAVE(MAXTDR)
C     RIIAVE is used in time averaged distant restraining
COMMEND


COMMVAR IPLR,JPLR,KPLR,LPLR,ICPLR,NDLR,CPLR,PDLR
C*****dihedral restraining
      INTEGER IPLR,JPLR,KPLR,LPLR,ICPLR,NDLR
      COMMON /DIHR/ IPLR(MAXNDL),JPLR(MAXNDL),KPLR(MAXNDL),
     $     LPLR(MAXNDL),ICPLR(MAXNDL),NDLR

      REAL CPLR,PDLR
      COMMON /DIHRFP/ CPLR(MAXNDL),PDLR(MAXNDL)
COMMEND

COMMVAR ICOG,JCOG,NCONG
C     distance constraining
      INTEGER ICOG,JCOG,NCONG
      COMMON /DICO/ICOG(MAXCON),JCOG(MAXCON),NCONG
COMMEND

COMMVAR CONP,FCON
C     FCON is used to store and manipulate the
C     constraint forces
      REAL CONP,FCON
      COMMON /DICOFP/ CONP(MAXCON),FCON(MAXCON)
COMMEND

COMMVAR IPJV,JPJV,KPJV,LPJV,NDJV,CPJV,PJR0,PSJR,AJV,BJV,CJV,COSQAV,COSIAV
C     j-value restraining
      INTEGER NDJV
      INTEGER IPJV, JPJV, KPJV, LPJV
      COMMON /JVINTY/NDJV,
     $     IPJV(MAXJ), JPJV(MAXJ), KPJV(MAXJ), LPJV(MAXJ)

      REAL CPJV,PJR0,PSJR,AJV,BJV,CJV
      REAL COSQAV,COSIAV
      COMMON /JVREL/ CPJV(MAXJ),PJR0(MAXJ),PSJR(MAXJ),
     $     AJV(MAXJ),BJV(MAXJ),CJV(MAXJ),
     $     COSQAV(MAXJ),COSIAV(MAXJ)
COMMEND


COMMVAR NDLE,IPLE,JPLE,KPLE,LPLE,NLECFG,NLEMEM
C*****local elevation

C     The local elevation dihedral angles.
      INTEGER NDLE
      INTEGER IPLE,JPLE,KPLE,LPLE
      INTEGER ILEMEM,NLEVST

C     NLECFG: the number of INTS used to store one configuration
C     NLECFG <= L<MXLECF>
C     NLEMEM: the number of configurations in memory
C     NLEMEM <= L<MLECFG>

      INTEGER NLECFG,NLEMEM

      COMMON /LEINTY/NDLE,NLECFG,NLEMEM,
     $     IPLE(MAXDLE),JPLE(MAXDLE),KPLE(MAXDLE),LPLE(MAXDLE),
     $     ILEMEM(MXLECF,MLECFG),NLEVST(MLECFG)
COMMEND



COMMVAR NSP,NSPM
C     NSP(I) contains the last atom of submolecule I in the solute
C     NSPM: the number of submolecules in the solute
C
C     NSP is defined as MAXNSP+1 in order to allow more efficient code
C     in the nonbonded force calculation routines.
C     (a double IF statement is avoided in the virial calculation)
C     However, only the elements 1..MAXNSP are ever actually filled with
C     useful data.
      INTEGER NSP,NSPM
      COMMON /MNSP/ NSPM,NSP(MAXNSP+1)
COMMEND


COMMVAR C4D
C     Force consts of 4-th dim harmonic osc.
C     If the constants are positive or zero for an atom,
C     it is in 4D.
C     If it is strictly negative, the atom is in 3D.

      REAL C4D
      COMMON /FORCFP/C4D(MAX4DA)
COMMEND

COMMVAR GAM,CC1,CC2,CC3,CC4,CC5,CC6,CC7,CC8,CC9,SWINK,SWINKS
C these used in SD
      REAL GAM,CC1,CC2,CC3,CC4,CC5,CC6,CC7,CC8,CC9
      REAL SWINK,SWINKS
      COMMON /SDFRIC/GAM(MAXNAT),
     $     CC1(MAXNAT),CC2(MAXNAT),CC3(MAXNAT),CC4(MAXNAT),
     $     CC5(MAXNAT),CC6(MAXNAT),CC7(MAXNAT),CC8(MAXNAT),
     $     CC9(MAXNAT),SWINK(MAXNAT),SWINKS(MAXNAT)
COMMEND



