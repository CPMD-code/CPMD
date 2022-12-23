C $Id: restx.h,v 1.1 2006-12-27 11:31:40 itavern Exp $ -*-fortran-*-
C
C define the control flags for restx.f


COMMVAR NFRALL,NFREF,NFREDT,NFREEE,NFRDNR,NFRENR
C     define the values for the argument NF in subroutine L<RESTX>.
C
C        etot
C          forces
C            rest nrg
C                rest dev
C      NF
C      -1   1   1   1   1
C       0   1   1   0   0
C       1   1   0   1   1
C       2   1
C       3   EC contains dev of nrcth rest
C       4   EC contains nrg of nrcth rest

      INTEGER NFRALL, NFREF,NFREDT,NFREEE,NFRDNR,NFRENR

      PARAMETER (NFRALL =-1)
      PARAMETER (NFREF  = 0)
      PARAMETER (NFREDT = 1)
      PARAMETER (NFREEE = 2)
      PARAMETER (NFRDNR = 3)
      PARAMETER (NFRENR = 4)
COMMEND

