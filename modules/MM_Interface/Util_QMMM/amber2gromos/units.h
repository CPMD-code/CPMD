C $Id: units.h,v 1.1 2006-12-27 11:35:49 itavern Exp $ -*-fortran-*-
C      GROMOS UNIT DEFINITIONS
C

C     standard FORTRAN units
      INTEGER ISTDIN,ISTDOT
      PARAMETER (ISTDIN = 5,ISTDOT = 6)

C-------define reserved unit numbers--------------
C     md run control file
      INTEGER IOMDCT
      PARAMETER (IOMDCT = ISTDIN)

C     final coords and velocities
      INTEGER IOXVE
      PARAMETER (IOXVE = 11)

C     coordinates trajectory
      INTEGER IOTRJX
      PARAMETER (IOTRJX = 12)

C     velocities trajectory
      INTEGER IOTRJV
      PARAMETER (IOTRJV = 13)

C reserved for future use
      INTEGER IORSVD
      PARAMETER (IORSVD = 14)

C     energies etc. trajectory
      INTEGER IOTRJE
      PARAMETER (IOTRJE = 15)

C     free energies
      INTEGER IOTRJG
      PARAMETER (IOTRJG = 16)



C input files start from unit 20 upwards
C     formatted GROMOS95 topology files
      INTEGER IOTOPO
      PARAMETER (IOTOPO = 20)

C     initial coordinates and velocities
      INTEGER IOXVI
      PARAMETER (IOXVI = 21)

C     reference coordinates for position re(con)straining
      INTEGER IOREST
      PARAMETER (IOREST = 22)

C     sequence numbers of position re(con)strained atoms
      INTEGER IORSTA
      PARAMETER (IORSTA = 23)

C     distance restrained atom pairs
      INTEGER IORSTP
      PARAMETER (IORSTP = 24)

C     restrained dihedral specifications
      INTEGER IODIHE
      PARAMETER (IODIHE = 25)

C     j-value restraining specifications
      INTEGER IOJVSP
      PARAMETER (IOJVSP = 26)

C     local elevation dihedral specifications
      INTEGER IOLESP
      PARAMETER (IOLESP = 27)

C     sequence numbers indicating atoms in 4D
      INTEGER IO4NDX
      PARAMETER (IO4NDX = 28)

C     friction coefficients for SD
      INTEGER IOGAM
      PARAMETER (IOGAM = 29)

C     data determining perturbation
      INTEGER IOPERT
      PARAMETER (IOPERT = 30)

C     some conversion programs (GROMOS87 --> GROMOS95 file format)
C     read from unit 40
      INTEGER IOCNV
      PARAMETER (IOCNV = 40)


C     numbers 50 to 59 are reserved for developers

      INTEGER IOD00,IOD01,IOD02,IOD03,IOD04
      INTEGER IOD05,IOD06,IOD07,IOD08,IOD09
      PARAMETER (IOD00 = 50)
      PARAMETER (IOD01 = 51)
      PARAMETER (IOD02 = 52)
      PARAMETER (IOD03 = 53)
      PARAMETER (IOD04 = 54)
      PARAMETER (IOD05 = 55)
      PARAMETER (IOD06 = 56)
      PARAMETER (IOD07 = 57)
      PARAMETER (IOD08 = 58)
      PARAMETER (IOD09 = 59)




C------consecutive numbers for table entries
C     The tables are used for mapping a string
C     to a unit number in SUBR. OPNFIL.
C     The table is scanned in a linear fashion which allows
C     the entries to be in any order.

      INTEGER JNMDCT
      INTEGER JNTRJX
      INTEGER JNTRJV
      INTEGER JNTRJE
      INTEGER JNTOPO
      INTEGER JNXVI
      INTEGER JNREST
      INTEGER JNRSTA
      INTEGER JNRSTP
      INTEGER JNPERT
      INTEGER JNDIHE
      INTEGER JNXVE
      INTEGER JND00,JND01,JND02,JND03,JND04
      INTEGER JND05,JND06,JND07,JND08,JND09
      INTEGER JN4COR,JN4NDX
      INTEGER JNCNV
      INTEGER JNGAM
      INTEGER JNJVSP
      INTEGER JNLESP
      INTEGER JNTRJG


      PARAMETER ( JNMDCT =  1)
      PARAMETER ( JNTRJX =  2)
      PARAMETER ( JNTRJV =  3)
      PARAMETER ( JNTRJE =  4)
      PARAMETER ( JNTOPO =  5)
      PARAMETER ( JNXVI  =  6)
      PARAMETER ( JNREST =  7)
      PARAMETER ( JNRSTA =  8)
      PARAMETER ( JNRSTP =  9)
      PARAMETER ( JNPERT = 10)
      PARAMETER ( JNDIHE = 11)
      PARAMETER ( JNXVE  = 12)
      PARAMETER ( JN4COR = 13)
      PARAMETER ( JN4NDX = 14)
      PARAMETER ( JNCNV  = 15)
      PARAMETER ( JNGAM  = 16)
      PARAMETER ( JNJVSP = 17)
      PARAMETER ( JNLESP = 18)
      PARAMETER ( JNTRJG = 19)

C add other types values here.
C remember to change NLSRES accordingly!

      INTEGER NLSRES
      PARAMETER (NLSRES = JNTRJG)

      PARAMETER ( JND00  = NLSRES + 1)
      PARAMETER ( JND01  = NLSRES + 2)
      PARAMETER ( JND02  = NLSRES + 3)
      PARAMETER ( JND03  = NLSRES + 4)
      PARAMETER ( JND04  = NLSRES + 5)
      PARAMETER ( JND05  = NLSRES + 6)
      PARAMETER ( JND06  = NLSRES + 7)
      PARAMETER ( JND07  = NLSRES + 8)
      PARAMETER ( JND08  = NLSRES + 9)
      PARAMETER ( JND09  = NLSRES + 10)

C     the number of reserved units
      INTEGER MAXNTS
      PARAMETER (MAXNTS = JND09)


      
C     The arrays of these common blocks are initialized
C     in a blockdata statement in fileio.f

      INTEGER MAXFNM
      PARAMETER (MAXFNM = 8)

      CHARACTER*(MAXFNM) UNAME
      COMMON /UNCHAR/UNAME(MAXNTS)

      INTEGER IUNUM
      COMMON /UNINT/IUNUM(MAXNTS)
