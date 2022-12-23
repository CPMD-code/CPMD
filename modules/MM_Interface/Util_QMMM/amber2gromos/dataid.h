C $Id: dataid.h,v 1.1 2006-12-27 11:20:08 itavern Exp $ -*-fortran-*-
C GROMOS BLOCK IDENTIFIERS
C

C-----------------------------------------------------
C     Note that these number are never stored to file,
C     instead only strings defined in NAMEID below
C     are ever written. This holds both for
C     the formatted and unformatted case.
C-----------------------------------------------------

C     end of block marker for formatted files
C     not used in unformatted files
      INTEGER IENDID
      PARAMETER (IENDID = 1)

C     format block. Used in unformatted files only.
C     This is used to determine the precision
C     (real4 or real8) of data in the file.
      INTEGER IFRMID
      PARAMETER (IFRMID = 2)

C     title block
      INTEGER  ITITID
      PARAMETER (ITITID = 3)

C     box block
      INTEGER IBOXID
      PARAMETER (IBOXID = 4)

C     positions, only coords X,Y,Z.
C     Formatted case only.
C     If we have NDIM >3 then the W coords are written
C     into a separate block.
      INTEGER IPOSID
      PARAMETER (IPOSID = 5)

      INTEGER IPORID
      PARAMETER (IPORID = 6)

C     positions in the 4th dim
      INTEGER IP4DID
      PARAMETER (IP4DID = 7)

      INTEGER IP4RID
      PARAMETER (IP4RID = 8)

C     position restraining COORDINATES block
      INTEGER IPRCID
      PARAMETER (IPRCID = 9)

C     oblique fractional coordinate block
      INTEGER IPOBID
      PARAMETER (IPOBID = 10)


C     vel block, only coords X,Y,Z.
C     If we have NDIM >3 then the W coords are written
C     into a separate block
      INTEGER IVELID
      PARAMETER (IVELID = 11)

      INTEGER IVRDID
      PARAMETER (IVRDID = 12)


      INTEGER IV4DID
      PARAMETER (IV4DID = 13)

      INTEGER IV4RID
      PARAMETER (IV4RID = 14)


C stochastic dynamics block
      INTEGER ISDIID
      PARAMETER (ISDIID = 15)

      INTEGER ISD4ID
      PARAMETER (ISD4ID = 16)

      INTEGER IPRTID
      PARAMETER (IPRTID = 17)

C     block containing the time averaged distances
C     in distance restraining
      INTEGER IRIIAV
      PARAMETER (IRIIAV = 18)

C     averaged j-value restraint block
      INTEGER ICOSAV
      PARAMETER (ICOSAV = 19)

C     local elevation memory
      INTEGER ILEMID
      PARAMETER (ILEMID = 20)

C     time and step number block
      INTEGER ITIMID
      PARAMETER (ITIMID = 21)

C     isotropic atomic B factors
      INTEGER IABFID
      PARAMETER (IABFID = 22)

C     anisotropic atomic B factors
      INTEGER IANBID
      PARAMETER (IANBID = 23)

C     second moment of coordinate distribution
      INTEGER IP2MID
      PARAMETER (IP2MID = 24)

C     third moment of coordinate distribution
      INTEGER IP3MID
      PARAMETER (IP3MID = 25)

C     fourth moment of coordinate distribution
      INTEGER IP4MID
      PARAMETER (IP4MID = 26)

C     trace of the second moments of coordinate distribution
      INTEGER IP2TID
      PARAMETER (IP2TID = 27)

      INTEGER IQTAID
      PARAMETER (IQTAID = 28)

      INTEGER IQEAID
      PARAMETER (IQEAID = 29)

      INTEGER IQESID
      PARAMETER (IQESID = 30)

      INTEGER IQTSID
      PARAMETER (IQTSID = 31)

      INTEGER IQDSID
      PARAMETER (IQDSID = 32)

      INTEGER IQTCID
      PARAMETER (IQTCID = 33)

      INTEGER IQSDID
      PARAMETER (IQSDID = 34)


      INTEGER IXJVID
      PARAMETER (IXJVID = 35)

      INTEGER IXLEID
      PARAMETER (IXLEID = 36)

      INTEGER IXFRID
      PARAMETER (IXFRID = 37)

      INTEGER IX4DID
      PARAMETER (IX4DID = 38)

C     position restraining ATOM SEQUENCE NUMBER block
      INTEGER IXPRID
      PARAMETER (IXPRID = 39)

C     block containing dihedral restraint data
C     used for dehedral restraining
      INTEGER IDHBLK
      PARAMETER (IDHBLK = 40)

C     block containing distance restraining info
C     ( e.g. definitions of pseudo atoms etc.)
      INTEGER IDRBLK
      PARAMETER (IDRBLK = 41)

C     energies block
      INTEGER INRGID
      PARAMETER (INRGID = 42)

C     volume-pressure block
      INTEGER IDVPRT
      PARAMETER (IDVPRT = 43)

C     free energies block
      INTEGER IDRLAM
      PARAMETER (IDRLAM = 44)

C     3d-4d free energy block
      INTEGER IDRMU
      PARAMETER (IDRMU = 45)


C     solvent statistics block
      INTEGER ISVSID
      PARAMETER (ISVSID = 46)
 
C     difference statistics block
      INTEGER IDFSID
      PARAMETER (IDFSID = 47)
 
C     dipole moment statistics block
      INTEGER IDMSID
      PARAMETER (IDMSID = 48)

C     add other block types here.
C     remember to change MAXIDT in accordance!

C     maximum number of block types
C     IDUKN is used by RDBHDR to denote an unknown block type, i.e.
C     one that is not in the table NAMEID.
C     IDERR is returned on a reading error, such as is
C     encountered at end of file.

      INTEGER IDERR,IDUKN,MINIDT,MAXIDT
      PARAMETER (IDERR = 0,IDUKN = -1,MINIDT = IENDID, MAXIDT = IDMSID)

C     maximum length of block type names
C     It is not recommended to change this constant as doing so
C     would make it difficult to read in unformatted GROMOS
C     files.
      INTEGER MAXIDL
      PARAMETER (MAXIDL = 16)

C NAMEID is the table containing the reserved blocknames
      CHARACTER* (MAXIDL) NAMEID

C IDLN contains the block name of the last block read by
C RDBHDR (see blockio.f)
C The content is only valid if the block type returned
C by RDBHDR is not IDERR.
      CHARACTER IDLN*(MAXIDL)
C
      CHARACTER  FWFAIL*(48)
      CHARACTER  FRFAIL*(48)
      CHARACTER  FMEXP*(48)
      CHARACTER  FMG2*(48)
      CHARACTER  FMIGNO*(28)
C
      CHARACTER  ENDEXP*(20)
      CHARACTER  BONFIL*(24)
      CHARACTER  MAXVAL*(24)
      CHARACTER  EXPFIL*(24)
      CHARACTER  FLGBIN*(24)
      CHARACTER  MTYLIN*(64)
      CHARACTER  STREXP*(24)
C
      COMMON /IDTAB/ NAMEID(MAXIDT),IDLN
      COMMON /IDDTAB/FRFAIL,FWFAIL,FMEXP,FMG2,FMIGNO,
     $     ENDEXP,BONFIL,EXPFIL,FLGBIN,MAXVAL,MTYLIN,STREXP

C     the arrays NAMEID, FWFAIL, FRFAIL FMEXP and FMG2 are initialized
C     in a blockdata statement in blockio.f


C     used to check correct binary formats
C     See WRFMT and RDFMT
      REAL RMAGIC
      PARAMETER (RMAGIC = -17.2E23)

