C $Id: ffblock.h,v 1.1 2006-12-27 11:22:44 itavern Exp $ -*-fortran-*-
C block types in a force field file




C define the legal block types here
      INTEGER NFFTIT,NFFEND
      INTEGER NFFIMC,NFFBTY,NFFBAT
      INTEGER NFFDIH,NFFIMD
      INTEGER NFFSAT,NFFAPR
      INTEGER NFFCNT

      INTEGER MINNFF,MAXNFF,NFFERR

      PARAMETER (NFFEND = 1)
      PARAMETER (NFFTIT = 2)
      PARAMETER (NFFIMC = 3)
      PARAMETER (NFFBTY = 4)
      PARAMETER (NFFBAT = 5)
      PARAMETER (NFFDIH = 6)
      PARAMETER (NFFIMD = 7)
      PARAMETER (NFFSAT = 8)
      PARAMETER (NFFAPR = 9)
      PARAMETER (NFFCNT = 10)



C     add other block types here.
C     remember to change MAXNFF in accordance!
      PARAMETER (NFFERR = 0,MINNFF = NFFEND, MAXNFF = NFFCNT)


C     maximum length of block type names
      INTEGER MFFBLE
      PARAMETER (MFFBLE = 16)

      LOGICAL LGOT
      COMMON /FFCOM/LGOT(MAXNFF)

      CHARACTER *(MFFBLE) FFNAME
      COMMON /FFCH/FFNAME(MAXNFF)


