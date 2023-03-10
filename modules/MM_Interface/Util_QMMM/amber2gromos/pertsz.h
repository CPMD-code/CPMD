C $Id: pertsz.h,v 1.1 2006-12-27 11:29:05 itavern Exp $ -*-fortran-*-
C size limits for perturbation files


COMMVAR NPTTIL,NPTTLN
C     maximum length of title string and
C     maximum number of lines of title
      INTEGER NPTTIL,NPTTLN
      PARAMETER (NPTTIL = 80, NPTTLN = 10)
COMMEND

COMMVAR MAXPAT
C     MAXIMUM NUMBER OF PERTURBED ATOMS
      INTEGER MAXPAT
      PARAMETER (MAXPAT = 80)
COMMEND

COMMVAR MAXPPA
C     MAXIMUM NUMBER OF SPECIFIED PERTURBED ATOM PAIRS
      INTEGER MAXPPA
      PARAMETER (MAXPPA = 10)
COMMEND

COMMVAR MAXPBO,MAXPTH,MAXPQH,MAXPPH
C     MAXIMUM NUMBER OF PERTURBED BONDS, BOND ANGLES, IMPROPER
C     DIHEDRALS AND DIHEDRALS

      INTEGER MAXPBO,MAXPTH,MAXPQH,MAXPPH
      PARAMETER (MAXPBO = 5)
      PARAMETER (MAXPTH = 5)
      PARAMETER (MAXPQH = 5)
      PARAMETER (MAXPPH = 5)
COMMEND


COMMVAR MAXPPI
C**** PATH INTEGRAL PERTURBATION
      INTEGER MAXPPI
      PARAMETER (MAXPPI = 10)
COMMEND

