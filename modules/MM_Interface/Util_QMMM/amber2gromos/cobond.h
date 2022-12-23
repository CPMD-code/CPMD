C $Id: cobond.h,v 1.1 2006-12-27 11:16:29 itavern Exp $ -*-fortran-*-

C predefined values for routines in cobond.f:
C COBOND, ANGLE and DIHANG.


COMMVAR NFDER,NFALL,NFEF,NFEL,NFE,NFNDXL,NFNDXE,NFCMAX,NFCMIN,NFCDMN
C     values for NFCTR used in all subroutines of
C     cobond.f (L<COBOND>, L<ANGLE>, L<DIHANG>)
C     exceptions are indicated

C only DIHANG
      INTEGER NFDER
      PARAMETER (NFDER = -2)

C all
      INTEGER NFALL,NFEF,NFEL,NFE,NFNDXL,NFNDXE
      PARAMETER (NFALL  =-1)
      PARAMETER (NFEF   = 0) 
      PARAMETER (NFEL   = 1)
      PARAMETER (NFE    = 2)
      PARAMETER (NFNDXL = 3)
      PARAMETER (NFNDXE = 4)

      INTEGER NFCMAX, NFCMIN, NFCDMN
      PARAMETER (NFCMAX = NFNDXE, NFCMIN = NFALL, NFCDMN = NFDER)

C NFCTR is   NB means               NM means
C
C  -2        num of d's            num of molecules       dphi/dr in force array,
C                                                         dihedral angle in E
C  -1        num of b/a/d's        num of molecules       E, force, b/a/d E, b/a/d
C   0        num of b's            num of molecules       E, force
C   1        num of b's            num of molecules       E,        b/a/d E, b/a/d 
C   2        num of bonds          num of molecules       E
C   3        seq number            seq number molecule     calc  b/a/d
C   4        seq number            seq number molecule     calc  b/a/d E
C
C b: bond; a: bond angle; d: dihedral angle
COMMEND

COMMVAR IBND4,IBND2,IBNDH,IBNMAX, IBNMIN
C     values for IDIBND in SUBR. L<COBOND>
      INTEGER IBND4,IBND2,IBNDH,IBNMAX, IBNMIN

      PARAMETER (IBND4 = 0)
      PARAMETER (IBND2 = 1)
      PARAMETER (IBNDH = 2)
      PARAMETER (IBNMAX = IBNDH, IBNMIN = IBND4)
COMMEND

COMMVAR IANGC,IANGH,IANMAX,IANMIN
C     values for IDIANG in SUBR. L<ANGLE>
      INTEGER IANGC,IANGH,IANMAX,IANMIN

      PARAMETER (IANGC = 0)
      PARAMETER (IANGH = 1)
      PARAMETER (IANMAX = IANGH, IANMIN = IANGC)
COMMEND

COMMVAR IDICOS,IDINCS
C     values for IDI in SUBR.  L<DIHANG>

      INTEGER IDICOS, IDINCS,IDIMIN, IDIMAX

C this is used for dihedrals
      PARAMETER (IDICOS = 1)
C this is used for improper dihedrals
      PARAMETER (IDINCS = 2)

      PARAMETER (IDIMIN = IDICOS, IDIMAX = IDINCS)
COMMEND


