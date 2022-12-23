C $Id: forceszl.h,v 1.1 2006-12-27 11:25:35 itavern Exp $ -*-fortran-*-
C----LARGE VERSION----

C----------------------------------------------------------
C remember to include "coordsz.h" before this file
C----------------------------------------------------------

COMMVAR MAXINB,MAXJNB
C     size of pairlist
C     MAXINB should be at least as large as the total number
C     of charge groups in the system (i.e. NSM + NPM*NCAG)  when
C     using a charge group pairlist.
C
      INTEGER MAXINB,MAXJNB

C     PARAMETER (MAXINB = 1000)
C     PARAMETER (MAXJNB = 6000)

      PARAMETER (MAXINB = 20000)
      PARAMETER (MAXJNB = 2000000)
COMMEND

COMMVAR MXATCG,MAXLAT
C     the maximum number of atoms in a charge group
C     used for local copies of atom parameters in the
C     non-bonded force calculations L<NONBML> and L<NBPML>
C
      INTEGER MXATCG,MAXLAT
      PARAMETER (MXATCG = 30,MAXLAT = MXATCG*MAXDIM)
COMMEND

COMMVAR MXEPST,MXEPFT
C     the length of a title string for writing out Epot values
      INTEGER MXEPST,MXEPFT
      PARAMETER(MXEPST = 16)

C     the length of the format string for writing out Epot values
      PARAMETER (MXEPFT = 48)
COMMEND



COMMVAR IETOT,IKTOT4,IKTOT3,IKSLU4,IKSLU3,IKSUCM,IKSLV4,IKSLV3,IPTOT,MXNRGF
C------------------------------------------------
C     predefined indices for accessing the energy arrays
C     ENER, FREM, FRE3D4
C     NOTE: the order in which the values are stored in the
C     array are significant.
C------------------------------------------------
      
      INTEGER IETOT,IKTOT4,IKTOT3,IKSLU4,IKSLU3,IKSUCM
      INTEGER IKSLV4,IKSLV3,IPTOT,MXNRGF

C     total energy Ekin + EPOT
      PARAMETER (IETOT = 1)

C     total kinetic energy in 4D
      PARAMETER (IKTOT4 = 2)

C     total kinetic energy in 3D
      PARAMETER (IKTOT3 = 3)

C     Ekin of solute in 4D
      PARAMETER (IKSLU4 = 4)

C     Total Ekin of solute in 3D
      PARAMETER (IKSLU3 = 5)

C     Ekin of centre of mass of solute (in 3D)
      PARAMETER (IKSUCM = 6)

C     Ekin of solvent in 4D
      PARAMETER (IKSLV4 = 7)

C     Ekin of solvent in 3D
      PARAMETER (IKSLV3 = 8)

C     total potential energies, i.e.,
C     E(IPTOT) = E(IPCOVT) + E(IPEL) + E(IPLJ) + E(IPRF) + E(IPRC) +E(IPOT4D)
      PARAMETER (IPTOT = 9)

C     this constant defines the number of values of the energy array
C     written to free energy trajectories
      PARAMETER (MXNRGF = IPTOT)
COMMEND

COMMVAR IPBH,IPBN,IPBAH,IPBA,IPQH,IPQE,IPPH,IPPE,IPLJ,IPEL,IPRF,IPRC
C     bonds
      INTEGER IPBH,IPBN
      PARAMETER (IPBH = MXNRGF+1, IPBN = MXNRGF+2)

C     bond angles
      INTEGER IPBAH,IPBA
      PARAMETER (IPBAH = MXNRGF+3, IPBA = MXNRGF+4)

C     improper dihedrals
      INTEGER IPQH,IPQE
      PARAMETER (IPQH = MXNRGF+5, IPQE = MXNRGF+6)

C     dihedrals
      INTEGER IPPH,IPPE
      PARAMETER (IPPH = MXNRGF+7,IPPE = MXNRGF+8)

C     Coulomb, reaction field and Lennard Jones
      INTEGER IPLJ,IPEL,IPRF,IPRC
      PARAMETER (IPLJ = MXNRGF+ 9,IPEL = MXNRGF+10)
      PARAMETER (IPRF = MXNRGF+11,IPRC = MXNRGF+12)
COMMEND


COMMVAR MNCONT, MXCONT
C     these bounds define the values in the array that
C     are used in force to calculate the total potential energy
C
      INTEGER MNCONT, MXCONT
      PARAMETER (MNCONT = IPBH, MXCONT = IPRC)
COMMEND

COMMVAR IPPISP
C     path integral: energy contained in the harmonic 'springs'
C     This energy is NOT part of the potential energy
      INTEGER IPPISP
      PARAMETER (IPPISP = MXCONT+1)
COMMEND

COMMVAR MXEWRT
C     This defines the limit of the energy arrays written to
C     the energy trajectory. Any higher values are only used
C     within the program ( e.g. for printing to screen) but
C     NOT written to file.

      INTEGER MXEWRT
      PARAMETER (MXEWRT = IPPISP)
COMMEND

C     total bond energy
      INTEGER IPBTOT
      PARAMETER (IPBTOT = MXEWRT+1)

C     total bond angle energy
      INTEGER IPBATO
      PARAMETER (IPBATO = MXEWRT+2)

C     total improper dihedral energy
      INTEGER IPQTOT
      PARAMETER (IPQTOT = MXEWRT+3)

C     total dihedral energy
      INTEGER IPPTOT
      PARAMETER (IPPTOT = MXEWRT+4)

C     total epot due to terms with H atoms
      INTEGER IPHTOT
      PARAMETER (IPHTOT = MXEWRT+5)

C     total epot due to terms whith NO H atoms
      INTEGER IPNHTO
      PARAMETER (IPNHTO = MXEWRT+6)

C     total epot due to all covalent terms
C     i.e. E(IPCOVT) = E(IPHTOT) + E(IPNHTO)
      INTEGER IPCOVT
      PARAMETER (IPCOVT = MXEWRT+7)

C     add any new type here.
C     remember to change MAXETBL and possibly
C     MNCONT, MAXCONT and MXEWRT in accordance!!

C     total epot without the "special forces"
      INTEGER IPIPSP
      PARAMETER (IPIPSP = MXEWRT+8)

C     total epot of the special forces
      INTEGER IPSPEC
      PARAMETER (IPSPEC = MXEWRT+9)

C     total nonbonded epot due to charges
C     I.e. EE(IPELEC) = EE(IPEL)+EE(IPRF)+EE(IPRC)

      INTEGER IPELEC
      PARAMETER (IPELEC = MXEWRT+10)

C     total nonbonded epot
C     I.e.: EE(IPNBON) = EE(IPEL)+EE(IPRF)+EE(IPRC)+EE(IPLJ)
      INTEGER IPNBON
      PARAMETER (IPNBON = MXEWRT+11)

C     total kinetic energy in the fourth dimension
C     i.e.: EE(IKSLV4) - EE(IKSLV3) + EE(IKSLU4) - EE(IKSLU3)
      INTEGER IK4THD
      PARAMETER (IK4THD = MXEWRT+12)

C     add any new type here.
C     remember to change MAXETBL and possibly
C     MNCONT, MAXCONT and MXEWRT in accordance!!

      INTEGER MXETBL
      PARAMETER (MXETBL = IK4THD)

C------------------------------------------------
C indices used in array ENERES
C------------------------------------------------
C     position restraint energy
      INTEGER ICPOSR
      PARAMETER (ICPOSR = 1)

C     distance restraint energy
      INTEGER ICDISR
      PARAMETER (ICDISR = 2)

C     dihedral restraint energy
      INTEGER ICDHRE
      PARAMETER (ICDHRE = 3)

C     J value restraining energy
      INTEGER ICJVAL
      PARAMETER (ICJVAL = 4)

C     local elevation energy
      INTEGER ICLOCE
      PARAMETER (ICLOCE = 5)

C     fourth dimension restraining energy
      INTEGER ICRE4D
      PARAMETER (ICRE4D = 6)

C     add any new type here.
C     remember to change MAXCTBL accordingly

      INTEGER MXCTBL
      PARAMETER (MXCTBL = ICRE4D)



C     indices used in array VOLPRT

C     solute 3D scaling factor (internal and rotational)
      INTEGER IVSPIR
      PARAMETER (IVSPIR = 1)

C     solute 3D centre of mass velocity scaling factor
      INTEGER IVSPCM
      PARAMETER (IVSPCM = 2)

C     solvent 3D scaling factor
      INTEGER IVSCLS
      PARAMETER (IVSCLS = 3)

C     total 4th D scaling factor
      INTEGER IVSCL4
      PARAMETER (IVSCL4 = 4)

C     3D box lengths 
      INTEGER IVBOXX,IVBOXY,IVBOXZ
      PARAMETER (IVBOXX = 5,IVBOXY = 6,IVBOXZ = 7)

C     volume of 3D box
      INTEGER IVVOL
      PARAMETER (IVVOL = 8)

C     pressure components
      INTEGER IVPRSX,IVPRSY,IVPRSZ
      PARAMETER (IVPRSX = 9,IVPRSY = 10,IVPRSZ = 11)

C     total pressure
      INTEGER IVPRES
      PARAMETER (IVPRES = 12)

C     kinetic energy components of centre of mass
      INTEGER IVCMX,IVCMY,IVCMZ
      PARAMETER (IVCMX = 13, IVCMY= 14, IVCMZ = 15)

C     total kinetic energy of centre of mass
      INTEGER IVEKCM
      PARAMETER (IVEKCM = 16)

C     virial components
      INTEGER IVVIRX,IVVIRY,IVVIRZ
      PARAMETER (IVVIRX = 17,IVVIRY = 18,IVVIRZ = 19)

C     total virial
      INTEGER IVVIR
      PARAMETER (IVVIR = 20)



C any components larger that this value will NOT be
C written to file in a VOLPRT block
      INTEGER MXVWRT
      PARAMETER (MXVWRT = IVVIR)


      INTEGER IVRLAM,IVRMU
      PARAMETER (IVRLAM = MXVWRT+1, IVRMU = MXVWRT+2)

C     add any new type here.
C     remember to change MXVTBL accordingly

      INTEGER MXVTBL
      PARAMETER (MXVTBL = IVRMU)

