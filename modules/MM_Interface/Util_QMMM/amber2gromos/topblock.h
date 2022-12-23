C $Id: topblock.h,v 1.1 2006-12-27 11:33:28 itavern Exp $ -*-fortran-*-
C     block types in a topology file
C     this file included by rdtopo.f and wrtopo.f


C define the legal block types here

      INTEGER NPTITL,NATNAM,NRSNAM,NSOLAT
      INTEGER NBNDTY,NBND,NBNDHY
      INTEGER NBATY,NBAN,NBAHY
      INTEGER NIMDTY,NIMDH,NIMD
      INTEGER NDITY,NDIHY,NDIN
      INTEGER NLJBLK,NTPEND,NTPUNT
      INTEGER NSLVBK,NSLVCN
      INTEGER NTPVER,NPITBL


      PARAMETER (NPTITL = 1, NATNAM = 2, NRSNAM = 3, NSOLAT = 4)

      PARAMETER (NBNDTY = 5, NBND = 6, NBNDHY = 7)

      PARAMETER (NBATY = 8, NBAN = 9, NBAHY = 10)

      PARAMETER (NIMDTY = 11, NIMDH = 12, NIMD= 13)

      PARAMETER (NDITY = 14, NDIHY = 15, NDIN= 16)

      PARAMETER (NLJBLK = 17, NTPEND = 18, NTPUNT = 19)

      PARAMETER (NSLVBK = 20, NSLVCN = 21)

      PARAMETER (NTPVER = 22, NPITBL = 23)


C     add other block types here.
C     remember to change MAXTPB in accordance!
      INTEGER NTPERR,MINTPB,MAXTPB
      PARAMETER (NTPERR = 0, MINTPB = NPTITL, MAXTPB = NPITBL)




C     the current topology version which is written to
C     an NTPVER block
      REAL TPVER
      PARAMETER (TPVER = 1.7)


C     maximum length of block type names
      INTEGER MXTBLE
      PARAMETER (MXTBLE = 16)

      LOGICAL LNEED,LGOT
      COMMON /TOPCOM/LNEED(MAXTPB),LGOT(MAXTPB)

      CHARACTER *(MXTBLE) TPNAME
      COMMON /TPCH/TPNAME(MAXTPB)

