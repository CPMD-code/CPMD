C $Id: pertar.h,v 1.1 2006-12-27 11:28:38 itavern Exp $ -*-fortran-*-
C array definitions for perturbation
C
C--------------------------------------------------------
C REMEMBER TO INCLUDE pertsz.h  BEFORE INCLUDING THIS FILE
C---------------------------------------------------------

COMMVAR NJLA,JLA,IACB,ISCLJ,ISCC,CGB,WMA,WMB
C*****perturbed atoms
C     NJLA:  number of perturbed atoms
C     JLA :  atom sequence number of perturbed atom
C     IACB:  integer atom code for B state
C     ISCLJ: soft core for Lennard-Jones ?
C     ISCC : soft core for Coulomb ?
C     CGB  : charge in B state
C     WMA  : mass in A state
C     WMB  : mas in B state
C
      INTEGER NJLA,JLA,IACB,ISCLJ,ISCC
      REAL CGB,WMA,WMB
      COMMON /PTAINT/NJLA,JLA(MAXPAT),IACB(MAXPAT),
     $     ISCLJ(MAXPAT),ISCC(MAXPAT)
      COMMON /PTAREL/CGB(MAXPAT),WMA(MAXPAT),WMB(MAXPAT)
COMMEND

COMMVAR NEB,IEB,JEB,IETA,IETB
C*****perturbed atom pairs
C     NEB    :number of perturbed atom pairs
C     IEB,JEB: atom sequence numbers of atom pairs
C     IETA, IETB:
C     IETA(1..NEB) = 0 : NO INTERACTION IS CALCULATED FOR THIS PAIR
C                        IN STATE A (1-RLAM)
C                  = 1 : THE 6-12 POTENTIAL PARAMETERS ARE TAKEN FROM
C                        C6 AND C12 IN STATE A
C                  = 2 : THE 6-12 POTENTIAL PARAMETERS ARE TAKEN FROM
C                        CS6 AND CS12 (1-4 INTERACTION) IN STATE A

      INTEGER NEB,IEB,JEB,IETA,IETB
      COMMON /PTPPA/ NEB,IEB(MAXPPA),JEB(MAXPPA),
     $     IETA(MAXPPA),IETB(MAXPPA)

      INTEGER NPTNOA,NPTLJA,NPT14A
      PARAMETER (NPTNOA = 0, NPTLJA = 1, NPT14A = 2)
      INTEGER NPTIMI, NPTIMA
      PARAMETER (NPTIMI = NPTNOA, NPTIMA = NPT14A)
COMMEND


COMMVAR NBONHG,IBHG,JBHG,NCBHG,CBHA,BHA0,CBHB,BHB0
C*****perturbed bonds involving H-atoms
C     NBONHG: number of perturbed bonds involving H-atoms
C     IBHG,JBHG: atom sequence numbers of atoms 
C     NCBHG : bond number to perturb
C     CBHA,CBHB: force constants in states a and B
C     BHA0,BHB0: ideal bond lengths in states A and B
C
      INTEGER NBONHG,IBHG,JBHG,NCBHG
      REAL CBHA,BHA0,CBHB,BHB0
      COMMON /PTBNDH/ NBONHG,
     $     IBHG(MAXPBO),JBHG(MAXPBO),NCBHG(MAXPBO)
      COMMON /PTBREL/ CBHA(MAXPBO),BHA0(MAXPBO),
     $     CBHB(MAXPBO),BHB0(MAXPBO)
COMMEND

COMMVAR NBONG,IBG,JBG,NCBG,CBA,BA0,CBB,BB0
C*****perturbed bonds NOT involving H-atoms
C     NBONG: number of perturbed bonds involving H-atoms
C     IBG,JBG: atom sequence numbers of atoms 
C     NCBG : bond sequence number in mol. topology
C     CBA,CBB: force constants in states a and B
C     BA0,BB0: ideal bond lengths in states A and B
C
      INTEGER NBONG,IBG,JBG,NCBG
      REAL CBA,BA0,CBB,BB0
      COMMON /PTBONG/NBONG,
     $     IBG(MAXPBO),JBG(MAXPBO),NCBG(MAXPBO)
      COMMON /PTBORE/ CBA(MAXPBO),BA0(MAXPBO),
     $     CBB(MAXPBO),BB0(MAXPBO)
COMMEND

COMMVAR NTHEHG,ITHG,JTHG,KTHG,CTHA,THA0,CTHB,THB0
C*****perturbed bond angles involving H-atoms
C     NTHEHG: number of perturbed bond angles involving H-atoms
C     ITHG,JTHG,KTHG: atom sequence numbers
C     THA,THA0 : ideal bond angles in states A and B
C     CTHB,THB0: force constants in states A and B
C
      INTEGER NTHEHG,ITHG,JTHG,KTHG
      REAL CTHA,THA0,CTHB,THB0
      COMMON /PTHEHG/ NTHEHG,
     $     ITHG(MAXPTH),JTHG(MAXPTH),KTHG(MAXPTH)
      COMMON /PHEREL/ CTHA(MAXPTH),THA0(MAXPTH),
     $     CTHB(MAXPTH),THB0(MAXPTH)
COMMEND


COMMVAR NTHEG,ITG,JTG,KTG,CTA,TA0,CTB,TB0
C*****perturbed bond angles NOT involving H-atoms
C     NTHEG: number of perturbed bond angles NOT involving H-atoms
C     ITG,JTG,KTG: atom sequence numbers
C     TA,TA0 : ideal bond angles in states A and B
C     CTB,TB0: force constants in states A and B
C
      INTEGER NTHEG,ITG,JTG,KTG
      REAL CTA,TA0,CTB,TB0
      COMMON /PTTHEG/ NTHEG,
     $     ITG(MAXPTH),JTG(MAXPTH),KTG(MAXPTH)
      COMMON /PTTHRE/ CTA(MAXPTH),TA0(MAXPTH),
     $     CTB(MAXPTH),TB0(MAXPTH)
COMMEND

COMMVAR NQHIHG,IQHG,JQHG,KQHG,LQHG,CQHA,QHA0,CQHB,QHB0
C*****perturbed improper dihedrals involving H-atoms
C     NQHIHG: number of perturbed improper dihedrals involving H-atoms
C     IQHG,JQHG,KQHG,LQHG: atom sequence numbers
C     QHA0,QHB0: ideal angles in states A and B
C     CQHA,CQHB: force constants in states A and B
C
      INTEGER NQHIHG,IQHG,JQHG,KQHG,LQHG
      REAL CQHA,QHA0,CQHB,QHB0
      COMMON /PTQHIH/ NQHIHG,
     $     IQHG(MAXPQH),JQHG(MAXPQH),KQHG(MAXPQH),LQHG(MAXPQH)
      COMMON /PTQHRE/ CQHA(MAXPQH),QHA0(MAXPQH),
     $     CQHB(MAXPQH),QHB0(MAXPQH)
COMMEND


COMMVAR NQHIG,IQG,JQG,KQG,LQG,CQA,QA0,CQB,QB0
C*****perturbed improper dihedrals NOT involving H-atoms
C     NQHIG: number of perturbed improper dihedrals NOT involving H-atoms
C     IQG,JQG,KQG,LQG: atom sequence numbers
C     QA0,QB0: ideal angles in states A and B
C     CQA,CQB: force constants in states A and B
C
      INTEGER NQHIG,IQG,JQG,KQG,LQG
      REAL CQA,QA0,CQB,QB0
      COMMON /PTQHIG/NQHIG,
     $     IQG(MAXPQH),JQG(MAXPQH),KQG(MAXPQH),LQG(MAXPQH)
      COMMON /PTQHRE/
     $     CQA(MAXPQH),QA0(MAXPQH),CQB(MAXPQH),QB0(MAXPQH)
COMMEND

COMMVAR NPHIHG,IPHG,JPHG,KPHG,LPHG,NPHA,NPHB,CPHA,PDHA,CPHB,PDHB
C*****perturbed dihedrals involving H-atoms
C     NPHIHG: number of perturbed dihedrals involving H-atoms
C     IPHG,JPHG,KPHG,LPHG: atom sequence numbers
C     NPHA,NPHB: multiplicity in states A and B
C     PDHA,PDHB: phase shift in staes A anf B
C     CPHA,CPHB: force constants in states A and B
C
      INTEGER NPHIHG,IPHG,JPHG,KPHG,LPHG,NPHA,NPHB
      REAL CPHA,PDHA,CPHB,PDHB
      COMMON /PTPHIH/NPHIHG,
     $     IPHG(MAXPPH),JPHG(MAXPPH),
     $     KPHG(MAXPPH),LPHG(MAXPPH),NPHA(MAXPPH),NPHB(MAXPPH)
      COMMON /PTPTRE/
     $     CPHA(MAXPPH),PDHA(MAXPPH),CPHB(MAXPPH),PDHB(MAXPPH)
COMMEND


COMMVAR NPHIG,IPG,JPG,KPG,LPG,NPA,NPB,CPA,PDA,CPB,PDB
C*****perturbed dihedrals NOT involving H-atoms
C     NPHIG: number of perturbed dihedrals NOT involving H-atoms
C     IPG,JPG,KPG,LPG: atom sequence numbers
C     NPA,NPB: multiplicity in states A and B
C     PDA,PDB: phase shift in staes A anf B
C     CPA,CPB: force constants in states A and B
C
      INTEGER NPHIG,IPG,JPG,KPG,LPG,NPA,NPB
      REAL CPA,PDA,CPB,PDB
      COMMON /PTPHIG/NPHIG,
     $     IPG(MAXPPH),JPG(MAXPPH),
     $     KPG(MAXPPH),LPG(MAXPPH),NPA(MAXPPH),NPB(MAXPPH)
      COMMON /PTPHRE/
     $     CPA(MAXPPH),PDA(MAXPPH),CPB(MAXPPH),PDB(MAXPPH)
COMMEND


COMMVAR NPIBG,IPIBG,JPIBG,ICPIBG
C*****path integral perturbation
C     NPIBG: number of perturbed path integral bonds
C     IPIBG,JPIBG: tom sequence numbers
C     ICPIBG: bond type code
C     CPIA,CPIB: force constants in states A and B
C
      INTEGER NPIBG,IPIBG,JPIBG,ICPIBG
      REAL CPIA, CPIB
      COMMON /PTPIIN/ NPIBG,
     $     IPIBG(MAXPPI),JPIBG(MAXPPI),ICPIBG(MAXPPI)
      COMMON /PTPIRE/ CPIA(MAXPPI),CPIB(MAXPPI)
COMMEND


COMMVAR PTTITL,NUMPTI
C*****title on perturbation input file
      CHARACTER PTTITL*(NPTTIL)
      COMMON /PTITBK/PTTITL(NPTTLN)

C     the actual number of lines in the title
      INTEGER NUMPTI
      COMMON /PTOK/NUMPTI
COMMEND

