C $Id: forcear.h,v 1.1 2006-12-27 11:24:34 itavern Exp $ -*-fortran-*-

C-----------------------------------------------------------------
C remember to include "md.h" and "forcesz.h" before this file
C-----------------------------------------------------------------

COMMVAR NRELKP,NUSNRE,NUNRE2
C     These values are initialised in subroutine L<CHKNRE>
C     and used in subroutine L<PRNRG> to print out the
C     energy matrix.
C
C     NRELKP is also used in the nonbonded routines
C     L<NONBML> and L<NBPML> to facilitate
C     generating the energy matrix.

      INTEGER NRELKP,NUSNRE,NUNRE2
      COMMON /EPBLCK/NUSNRE,NUNRE2,NRELKP(MAXNRE,MAXNRE)
COMMEND

COMMVAR EPTSTR,EPTFMT
C     EPTSTR: headings for writing over the energy matrix
C     EPTFMT: format strings used to write out the enery matrix
C     These arrays are initialised in subroutine L<CHKNRE>
C     and used in subroutine L<PRNRG> to print out the
C     energy matrix.
C
      CHARACTER *(MXEPST) EPTSTR
      CHARACTER *(MXEPFT) EPTFMT
      COMMON /EPBLFP/
     $     EPTSTR(MAXNRE),EPTFMT(MAXNRE)
COMMEND


