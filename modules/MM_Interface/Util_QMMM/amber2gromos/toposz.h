C : toposz.h,v 1.13 1996/11/04 17:37:23 wscott Exp $ -*-fortran-*-
C include file for GROMOS size limits on topology
C

COMMVAR MAXATT,MXATT2
C maximum number of atom types
      INTEGER MAXATT,MXATT2
      PARAMETER (MAXATT = 400)
      PARAMETER (MXATT2 = MAXATT*(MAXATT+1)/2 )
COMMEND

COMMVAR MAXAA2
C maximum number of residues
      INTEGER MAXAA2
      PARAMETER (MAXAA2 = 200000)
COMMEND

COMMVAR MAXNRP,MAXNP2
C maximum number of atoms per solute
      INTEGER MAXNRP,MAXNP2

      PARAMETER (MAXNRP = 400000)
      PARAMETER (MAXNP2 = 2*MAXNRP)
COMMEND

COMMVAR MAXNBT,MAXBNH,MAXBON
C--------BONDS
C maximum number of covalent bond types
      INTEGER MAXNBT
      PARAMETER (MAXNBT = 200)

C maximum number of bonds involving H-atoms in the solute
      INTEGER MAXBNH
      PARAMETER (MAXBNH = 800000)

C maximum number of bonds NOT involving H-atoms in the solute
      INTEGER MAXBON
      PARAMETER (MAXBON = 800000)
COMMEND


COMMVAR MAXTTY,MXQHEH,MAXTHE
C---------BOND ANGLES
C maximum number of bond angle types
      INTEGER MAXTTY
      PARAMETER (MAXTTY = 400)

C maximum number of bond angles involving
C H-atoms in the solute
      INTEGER MXQHEH
      PARAMETER (MXQHEH = 400000)

C maximum number of bond angles NOT involving
C H-atoms in the solute
      INTEGER MAXTHE
      PARAMETER (MAXTHE = 400000)
COMMEND

COMMVAR MAXQTY,MAXHIH,MAXQHI
C---------IMPROPER DIHEDRALS
C maximum number of improper dihedral types
      INTEGER MAXQTY
      PARAMETER (MAXQTY = 200)

C maximum number of improper dihedrals involving
C H-atoms in the solute
      INTEGER MAXHIH
      PARAMETER (MAXHIH = 400000)

C maximum number of improper dihedrals NOT involving 
C H-atoms in the solute
      INTEGER MAXQHI
      PARAMETER (MAXQHI = 400000)
COMMEND


COMMVAR MAXPTY,MXPHIH,MAXPHI
C-----------DIHEDRALS
C maximum number of dihedral types
      INTEGER MAXPTY
      PARAMETER (MAXPTY = 200)

C maximum number of dihedrals involving
C H-atoms in the solute
      INTEGER MXPHIH
      PARAMETER (MXPHIH = 400000)

C maximum number of dihedrals NOT
C involving H-atoms in the solute
      INTEGER MAXPHI
      PARAMETER (MAXPHI = 400000)
COMMEND


COMMVAR MAXCAG,MAXAEX,MXEX14
C maximum number of charge groups in a solute molecule
      INTEGER MAXCAG
      PARAMETER (MAXCAG = 400000)

C maximum total number of exclusions in a solute molecule
      INTEGER MAXAEX
      PARAMETER (MAXAEX = 1000000)

C maximum number of third number atoms in a solute molecule
      INTEGER MXEX14
      PARAMETER (MXEX14 = 500000)
COMMEND


COMMVAR MAXNRS,MXCONS
C maximum number of atoms per solvent molecule
      INTEGER MAXNRS
      PARAMETER (MAXNRS = 3)

C maximum number of solvent constraints
      INTEGER MXCONS
      PARAMETER (MXCONS = 3)
COMMEND

COMMVAR MAXTIT,MAXLNS
C params defining the size of character arrays
C length of title string, and max number of lines allowed
      INTEGER MAXTIT, MAXLNS
      PARAMETER (MAXTIT = 80, MAXLNS = 10)
COMMEND


COMMVAR MAXTLE,MAXRLE,MAXNLE,MXNLE2,MXNLE3,MXNLE4
C----IT IS NOT ADVISED TO CHANGE THE VALUES OF THESE CONSTANTS
C length of atom type names
      INTEGER MAXTLE
      PARAMETER (MAXTLE = 5)

C length of residue type names
      INTEGER MAXRLE
      PARAMETER (MAXRLE = 5)

C length of atom name of solvent and solute atoms
      INTEGER MAXNLE
      PARAMETER (MAXNLE = 5)

C these used for pretty printing...
      INTEGER MXNLE2
      PARAMETER (MXNLE2 = 2*MAXNLE + 1)

      INTEGER MXNLE3
      PARAMETER (MXNLE3 = 3*MAXNLE + 2)

      INTEGER MXNLE4
      PARAMETER (MXNLE4 = 4*MAXNLE + 3)
COMMEND


COMMVAR MAXPIA,MAXPID,MAXPIT,MAXPIB,MAXPIW,MAXWR2
C-----PATH INTEGRAL
C maximum number of discretized atoms
      INTEGER MAXPIA
      PARAMETER (MAXPIA = 100)
 
C maximum number of discretizations
      INTEGER MAXPID
      PARAMETER (MAXPID = 10)
 
C maximum number of path integral 'bond' types, should be MAXPIA
      INTEGER MAXPIT
      PARAMETER (MAXPIT = MAXPIA)
 
C maximum number of path integral 'bonds': MAXPIA*MAXPID
      INTEGER MAXPIB
      PARAMETER (MAXPIB = 1000)
 
C maximum dimension of work arrays
      INTEGER MAXPIW
      PARAMETER (MAXPIW = 1000)
 
C maximum number of atoms forming a bond
      INTEGER MAXWR2
      PARAMETER (MAXWR2 = 4)
COMMEND
