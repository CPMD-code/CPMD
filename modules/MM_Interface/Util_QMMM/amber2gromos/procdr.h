C $Id: procdr.h,v 1.1 2006-12-27 11:30:06 itavern Exp $-*-fortran-*-
C
C     Parameter statements for PROCDR by Andrew Torda.
C     These parameters, where applicable,
C     now have the same values as in toposz.h

C --------------------------------------------------
C remember to include toposz.h before this one
C---------------------------------------------------
C

COMMVAR MAXRES, MAXAT
C     Max residues and atoms that can be read from topology by
C     program L<PROCDR>.
      INTEGER MAXRES, MAXAT
      PARAMETER (MAXRES=MAXAA2, MAXAT=1024)
COMMEND

COMMVAR MAXNOE
C     Max number of NOE's that can be read from input NOE's
C     by porgram L<PROCDR>
      INTEGER MAXNOE
      PARAMETER (MAXNOE=2000)
COMMEND

COMMVAR MAXCHN, MAXTYP
C     Max number of entries in conversion file and max
C     types of residue in porgram L<PORCDR>.
C
      INTEGER MAXCHN, MAXTYP
      PARAMETER (MAXCHN=200, MAXTYP=40)
COMMEND



      INTEGER MAXJNK
      PARAMETER (MAXJNK=MAXBON)

C     Max number of chars that can be on an input line
      INTEGER MAXLIN
      PARAMETER(MAXLIN=80)

C     Max length of a residue name
      INTEGER LENRES
      PARAMETER (LENRES = 8)

C     Max length of an atom name
      INTEGER LENAT
      PARAMETER (LENAT = MAXTLE)

