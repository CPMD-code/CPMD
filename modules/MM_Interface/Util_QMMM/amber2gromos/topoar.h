C $Id: topoar.h,v 1.1 2006-12-27 11:33:48 itavern Exp $ -*-fortran-*-
C     topoar.h
C     common block definitions for the GROMOS topology
C
C----------------------------------------------------------
C     Note: include the file toposz.h before including
C     this one.
C----------------------------------------------------------


COMMVAR IAC,IPERT,NRP,MPAC,NRATT,NRATT2
C*****atom types, atom interaction types
C     IAC   : integer atom code
C     IPERT : denotes whether an atom is perturbed or not
C     IPERT is NOT read in by L<RDTOPO>, but is set to L<NOPERT>.
C     The actual reading is performed (optionally) by
C     L<RDPERT> if perturbation is specified in the L<PROMD>
C     control file.
C
C     NRP   : number of atoms per solute molecule
C     MPAC  : interaction matrix
C     NRATT : number of atom types
C     NRATT2: is always NRATT*(NRATT+1)/2
C
      INTEGER IAC,IPERT,NRP,MPAC,NRATT,NRATT2
      COMMON/ATYP/
     $     IAC(MAXNRP),IPERT(MAXNRP),NRP,
     $     MPAC(MAXATT,MAXATT),NRATT,NRATT2
COMMEND

COMMVAR NOPERT
C     a value used to denote that an atom is not perturbed.
C     See array L<IPERT>
      INTEGER NOPERT
      PARAMETER (NOPERT = 0)
COMMEND

COMMVAR FFTYPE,TOPTIT,NTPLNS
C*****title, residue names
C     FFTYPE is an array of length MAXATT of char(MAXTLE)
C     TOPTIT: the title read in from the molecular topology file.
C     NTPLNS: the number of lines in the title
C     FFTYPE: the names of the residues
C
      CHARACTER*(MAXTLE) FFTYPE
      CHARACTER*(MAXTIT) TOPTIT
      COMMON /ATYPCH/
     $     FFTYPE (MAXATT),
     $     TOPTIT(MAXLNS)
      INTEGER NTPLNS
      COMMON /ATPINT/NTPLNS
COMMEND

COMMVAR C12,C6,CS12,CS6
C     C6, C12: normal interaction parameters
C     CS6,CS12: 1-4 interaction parameters
C
      REAL    C12,C6,CS12,CS6
      COMMON /ATYPFP/ 
     $     C12(MXATT2),C6(MXATT2),
     $     CS12(MXATT2),CS6(MXATT2)
COMMEND


COMMVAR PANM
C     PANM: the names of the solute atoms
C     PANM is an array of length MAXNRP of char(MAXNLE)
      CHARACTER*(MAXNLE) PANM
      COMMON /TOP1/
     $    PANM(MAXNRP)
COMMEND



COMMVAR NRAA2, AANM
C*****residue types
C     NRAA2: the number of residues in the topology
C     AANM : the names of the residues
C     AANM is an array of length MAXAA2 of char(MAXRLE)

      INTEGER NRAA2
      COMMON /RTYP/ NRAA2

      CHARACTER*(MAXRLE) AANM
      COMMON /RTYPCH/AANM(MAXAA2)
COMMEND


COMMVAR NBTY, CB,B0,IB,JB,ICB,NBON,IBH,JBH,ICBH,NBONH
C*****bonds
C     NBTY: number of bond types
      INTEGER NBTY
      REAL    CB,B0
      COMMON /BONDFP/
     $     CB(MAXNBT),B0(MAXNBT)

      INTEGER IB,JB,ICB,NBON
      COMMON /BOND/
     $     IB(MAXBON),JB(MAXBON),ICB(MAXBON),NBON,NBTY

      INTEGER IBH,JBH,ICBH,NBONH
      COMMON /BONH/IBH(MAXBNH),JBH(MAXBNH),ICBH(MAXBNH),NBONH
COMMEND


COMMVAR NTTY,CT,T0,ITH,JTH,KTH,ICTH,NTHEH,IT,JT,KT,ICT,NTHE
C*****angles
      INTEGER NTTY
      REAL    CT,T0
      COMMON /ANGLFP/CT(MAXTTY),T0(MAXTTY)

      INTEGER ITH,JTH,KTH,ICTH,NTHEH
      INTEGER IT,JT,KT,ICT,NTHE
      COMMON /ANGL/
     $     ITH(MXQHEH),JTH(MXQHEH),KTH(MXQHEH),ICTH(MXQHEH),NTHEH,
     $     IT(MAXTHE),JT(MAXTHE),KT(MAXTHE),ICT(MAXTHE),NTHE,NTTY
COMMEND


COMMVAR NQTY,CQ,Q0,IQ,JQ,KQ,LQ,ICQ,NQHI,IQH,JQH,KQH,LQH,ICQH,NQHIH
C*****improper dihedrals
      INTEGER NQTY
      REAL    CQ,Q0
      COMMON /IDHAFP/ CQ(MAXQTY),Q0(MAXQTY)

      INTEGER IQ,JQ,KQ,LQ,ICQ,NQHI
      COMMON/IDHA/
     $     IQ(MAXQHI),JQ(MAXQHI),KQ(MAXQHI),LQ(MAXQHI),
     $     ICQ(MAXQHI),NQHI,NQTY

      INTEGER IQH,JQH,KQH,LQH,ICQH,NQHIH
      COMMON/IDHH/
     $     IQH(MAXHIH),JQH(MAXHIH),KQH(MAXHIH),LQH(MAXHIH),
     $     ICQH(MAXHIH),NQHIH
COMMEND


COMMVAR NP,NPTY,CP,PD,IP,JP,KP,LP,ICP,NPHI,IPH,JPH,KPH,LPH,ICPH,NPHIH
C*****dihedrals
      INTEGER NP,NPTY
      REAL    CP,PD
      COMMON /DIHAFP/ CP(MAXPTY),PD(MAXPTY)
      COMMON /DIHTYP/ NP(MAXPTY),NPTY


      INTEGER IP,JP,KP,LP,ICP,NPHI
      COMMON/DIHA/
     $     IP(MAXPHI),JP(MAXPHI),KP(MAXPHI),LP(MAXPHI),
     $     ICP(MAXPHI),NPHI


      INTEGER IPH,JPH,KPH,LPH,ICPH,NPHIH
      COMMON/DIHH/
     $     IPH(MXPHIH),JPH(MXPHIH),KPH(MXPHIH),LPH(MXPHIH),
     $     ICPH(MXPHIH),NPHIH
COMMVAR


COMMVAR MRES,INC,NCAG,INE,KNE,INE14,KNE14,JSNE,NAEX,JSNE14,NAEX14
C*****residue numbers,mass, charge, exclusions, 1-4 interactions,
C     charge group definitions, interaction matrix (mpac)
C---
C     The variables INE,KNE and JSNE are organized as follows:
C     exclusions J of atom I are positioned at
C     JSNE(KNE(I)+1),...JSNE(KNE(I)+INE(I))
C     all J must be > I and in ascending order.
C
C     The variables INE14,KNE14 and JSNE14 are analogous.
C---

      INTEGER MRES,INC,NCAG,INE,KNE,INE14,KNE14
      INTEGER JSNE,NAEX,JSNE14,NAEX14
      COMMON/NON2/
     $     MRES(MAXNRP),
     $     INC(MAXCAG),NCAG,
     $     INE(MAXNRP),KNE(MAXNRP),INE14(MAXNRP),KNE14(MAXNRP),
     $     JSNE(MAXAEX),NAEX,
     $     JSNE14(MXEX14),NAEX14
COMMEND

COMMVAR CG,WINV,WMAS
C     CG  : charge of solute atoms
C     WINV: inverse mass of solute atoms
C     WMAS: mass of solute atoms
      REAL    CG,WINV,WMAS
      COMMON /NON2FP/
     $     CG(MAXNRP),WINV(MAXNRP),WMAS(MAXNRP)
COMMEND


COMMVAR IACS,NRAM,CGS,WINVS,WMASS,ANMS
C*****solvent molecule atom data
C     IACS : integer atom code for solvent atoms
C     NRAM : number of atoms in a solvent molecule
C     CGS  : charge of solvent atoms
C     WINVS: inverse mass of solvent atoms
C     WMASS: mass of solvent atoms
C     ANMS : names of solvent atoms

      INTEGER   IACS,NRAM
      COMMON/SOLV/ IACS(MAXNRS),NRAM

      REAL      CGS,WINVS,WMASS
      COMMON /SOLVFP/ CGS(MAXNRS),WINVS(MAXNRS),WMASS(MAXNRS)

      CHARACTER*(MAXRLE) ANMS
      COMMON /SOLVCH/ ANMS(MAXNRS)
COMMEND


COMMVAR ICONS,JCONS,NCONS,CONS
C*****solvent molecule constraints data
      INTEGER ICONS,JCONS,NCONS
      COMMON/SCON/
     $     ICONS(MXCONS),JCONS(MXCONS),NCONS
      
      REAL CONS
      COMMON /SCONFP/ CONS(MXCONS)
COMMEND

COMMVAR FPEPSI,HBAR
C*****units defined in the topology file
C     EPS0 is the permittivity of vacuum
C     FPEPSI: we store 1/(4* PI* EPS0) 
C     HBAR  : Planck's constant HBAR = H/(2* PI)
      REAL FPEPSI,HBAR
      COMMON /TOPNTS/FPEPSI,HBAR
COMMEND

COMMVAR NPIA,NPIA,NPID,IPIC,NPIT,NPIB,IPIB,JPIB,ICPIB,TPI,CPI,WMCL,BOLPI
C*****path integral topology data
C     NPIA:     number of path integral atoms
C     IPIA:     atom seq. numbers in PI atoms in the classical topology
C     NPID:     number of discretizations per atom
C     IPIC(I):  pseudoparticle number of the 'atom' I: if zero, the
C           'atom' is a classical atom
C     TPI:      temperature as basis for the harmonic 'spring' constants
C     NPIT:     number of harmonic 'spring' constants between
C           pseudoparticles
C     CPI(I):   harmonic 'spring' constant of code I
C     BOLPI:    Boltzmann constant for path integral
C     NPIB:     number of harmonic 'bonds' between pseudoparticles
C     IPIB(I):  'from' atom of 'bond' I
C     JPIB(I):  'to' atom of 'bond' I
C     ICPIB(I): 'bond' code of 'bond' I
 
      INTEGER NPIA, IPIA(MAXPIA), NPID, IPIC(MAXNRP), NPIT, NPIB,
     $     IPIB(MAXPIB), JPIB(MAXPIB), ICPIB(MAXPIB)
      REAL TPI, CPI(MAXPIT), WMCL(MAXPIA), BOLPI
      COMMON /PITINT/ NPIA, IPIA, NPID, IPIC, NPIT, NPIB, IPIB, JPIB,
     $     ICPIB
      COMMON /PITREL/ TPI, CPI, WMCL, BOLPI
COMMEND

