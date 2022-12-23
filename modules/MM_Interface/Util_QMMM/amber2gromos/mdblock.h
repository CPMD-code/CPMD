C $Id: mdblock.h,v 1.1 2006-12-27 11:27:22 itavern Exp $ -*-fortran-*-
C
C     this file defines constants and tables
C     for reading PROMD input files in the GROMOS format.

C     define the legal block types here
      INTEGER NTITLE,NBOUND,NTCOUP,NPCOUP,NNSP,NSTEP
      INTEGER NCOM, NSTART, NSHAKE, NFORCE,NPLIST
      INTEGER NPRINT, NWRITE, NPREST,NDREST
      INTEGER NPERT,NDHRES
      INTEGER NFOURD,NLONGR
      INTEGER NEMBLK,NSDBLK
      INTEGER NPIBLK,NSYST
      INTEGER NJVAL,NLEBLK

      PARAMETER ( NTITLE =  1)
      PARAMETER ( NBOUND =  2)
      PARAMETER ( NTCOUP =  3)
      PARAMETER ( NPCOUP =  4)
      PARAMETER ( NNSP   =  5)
      PARAMETER ( NCOM   =  6)
      PARAMETER ( NSTEP  =  7)
      PARAMETER ( NSHAKE =  8)
      PARAMETER ( NFORCE =  9)
      PARAMETER ( NPLIST = 10)
      PARAMETER ( NPRINT = 11)
      PARAMETER ( NWRITE = 12)
      PARAMETER ( NPREST = 13)
      PARAMETER ( NDREST = 14)
      PARAMETER ( NPERT  = 15)
      PARAMETER ( NDHRES = 16)
      PARAMETER ( NSTART = 17)
      PARAMETER ( NFOURD = 18)
      PARAMETER ( NLONGR = 19)
      PARAMETER ( NEMBLK = 20)
      PARAMETER ( NSDBLK = 21)
      PARAMETER ( NPIBLK = 22)
      PARAMETER ( NSYST  = 23)
      PARAMETER ( NJVAL  = 24)
      PARAMETER ( NLEBLK = 25)

C
C     add other block types here.
C     remember to change MAXBT in accordance!

      INTEGER MINBT, MAXBT,NBTERR
      PARAMETER (NBTERR = 0,MINBT = NTITLE, MAXBT = NLEBLK)


C     maximum length of block type names
      INTEGER MXBTLE
      PARAMETER (MXBTLE = 16)

      LOGICAL LNEDMD,LNEDEM,LGOT
      COMMON /NBT/ LNEDMD(MAXBT),LNEDEM(MAXBT),LGOT(MAXBT)

      CHARACTER* (MXBTLE) BTNAME
      COMMON /NBTCH/ BTNAME(MAXBT)


C     LNEDMD(K) is .TRUE. if a block type is compulsory for MD or
C     SD. LNEDEM(K) is .TRUE. if a block type is necessary for EM.
C     The arrays BTNAME, LNEDMD and LNEDEM are initialized in the
C     BLOCKDATA statement in rdmd.f
C     The array LGOT is initialized by SUBR. RDMD at runtime.




