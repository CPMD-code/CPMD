#if defined(__HASNT_OMP_COLLAPSE)
#define __COLLAPSE2
#else
#define __COLLAPSE2 collapse(2)
#endif
! ==================================================================
      SUBROUTINE RDTOPO
! ==--------------------------------------------------------------==
!     RDTOPO opens the predefined topology file for PROMD (unit 20),
!     calls L<RGTTOP> to read the topology and closes the file.
!     L<PRPTOP> is called in order to prepare the topology
!     for use in force calculation routines.
!    Revision cmb Tokyo/Strasbourg june 2014
! ==--------------------------------------------------------------==
C local vars
      INTEGER IUNIT
C begin
      CALL OPNFIL('TOPO','FORMATTED','OLD',IUNIT)
      IF (IUNIT .LT. 0)THEN
         PRINT *,'RDTOPO: failed to open topology file!'
         PRINT *
         STOP
      ENDIF
      CALL RGTTOP(IUNIT)
      CALL CLSFIL(IUNIT)
      CALL PRPTOP
C end rdtopo
      END

! ==================================================================
      SUBROUTINE GTTOPO(IUNIT)
! ==--------------------------------------------------------------==
C     GTTOPO reads a GROMOS topology file from an open
C     formatted file.
C     If any error occurs on reading, a message is
C     printed and STOP is invoked.
C
C     Typical usage would be:
C     CALL OPNFIL('TOPO','FORMATTED','OLD',IUNIT)
C     IF (IUNIT .LT. 0)THEN
C        PRINT *,PRGSTR,': failed to open topology file!'
C        PRINT *
C        STOP
C     ENDIF
C     CALL GTTOPO(IUNIT)
C     CALL CLSFIL(IUNIT)
C
C     The topology is read into the variables defined
C     in toposz.h and topoar.h .
C     L<PRPTOP> is called in order to prepare the topology
C     for use in force calculation routines.
!    Revision cmb Tokyo/Strasbourg june 2014
! ==--------------------------------------------------------------==
C args
      INTEGER IUNIT
C begin
      CALL RGTTOP(IUNIT)
      CALL PRPTOP
C end gttopo
      END

! ==================================================================
      SUBROUTINE RGTTOP(IUNIT)
! ==--------------------------------------------------------------==
C     SUBROUTINE RGTTOP(IUNIT)
C
C     The name comes from "Raw Get Topology".
C     A topology is read from a formatted open file IUNIT.
C     The topology is read into the variables defined
C     in toposz.h and topoar.h .
C     NO CONVERSION OF VARIABLES IS PERFORMED.
!    Revision cmb Tokyo/Strasbourg june 2014
! ==--------------------------------------------------------------==
C includes
      INCLUDE 'toposz.h'
      INCLUDE 'topoar.h'
      INCLUDE 'fileio.h'
      INCLUDE 'topblock.h'
C args
      INTEGER IUNIT
C local parameters
      INTEGER MAXSTR
      PARAMETER (MAXSTR=20)
C external
      INTEGER NDXTOP
      EXTERNAL NDXTOP
C commons
      LOGICAL LREDTO
      COMMON /TOPOK/LREDTO
C local vars
      INTEGER I,NBLOCK
      LOGICAL LOK,LGTALL
      CHARACTER PRGSTR*(6)
      CHARACTER STR*(MAXSTR)
      REAL TMP
      DATA PRGSTR /'RGTTOP'/
C begin
      CALL INITLN(IUNIT)

!$omp parallel do private(I)
      DO I=1,MAXNRP
         IPERT(I) = NOPERT
      ENDDO

!$omp parallel do private(I)
      DO I=MINTPB,MAXTPB
         LGOT(I) = .FALSE.
      ENDDO

C title
      CALL GETLN(IUNIT,LOK)
      NBLOCK = NDXTOP(FIOLIN)
      IF (NBLOCK .NE. NPTITL) THEN
         CALL FLAGLN(PRGSTR)
         PRINT *,'TITLE block must the first block in the file'
         PRINT *,FIOLIN
         STOP
      ENDIF
      CALL RDTIT(IUNIT,.TRUE.,-100,MAXLNS,TOPTIT,NTPLNS)
      LGOT(NPTITL) = .TRUE.

C units block
      CALL GETLN(IUNIT,LOK)
      NBLOCK = NDXTOP(FIOLIN)
      IF (NBLOCK .NE. NTPUNT) THEN
         CALL FLAGLN(PRGSTR)
         PRINT *,TPNAME(NTPUNT),
     $        ' block must the second block in the file'
         PRINT *,FIOLIN
         STOP
      ENDIF

      CALL GUNT(IUNIT,PRGSTR,FPEPSI,HBAR)
      LGOT(NTPUNT) = .TRUE.

C version block
      CALL GETLN(IUNIT,LOK)
      NBLOCK = NDXTOP(FIOLIN)
      IF (NBLOCK .NE. NTPVER) THEN
         CALL FLAGLN(PRGSTR)
         PRINT *,TPNAME(NTPVER),'block must the third block in the file'
         PRINT *,FIOLIN
         STOP
      ENDIF
      STR = 'TPVER'
      CALL GIMME(IUNIT,PRGSTR,STR)
      CALL CHPREL(TMP,PRGSTR,STR)
      IF (ABS(TMP-TPVER) .GT. 1.0E-03) THEN
         CALL FLAGLN(PRGSTR)
         PRINT *,'wrong version number on file!'
         PRINT '(A,F3.1)','expected : ',TPVER
         PRINT '(A,F3.1)','but read :',TMP
         PRINT *,'retranslate the GROMOS87 topology...'
         STOP
      ENDIF
      CALL CHKEND(IUNIT,PRGSTR)
      LGOT(NTPVER) = .TRUE.


C now read in the other blocks which may be in (almost) any order,
C however, a type block must be read in before any dependent block. E.g.
C bondtypes must be read in before any bonds are specified.

      CALL GETLN(IUNIT,LOK)
      NBLOCK = NDXTOP(FIOLIN)
C while not eof 
 10   IF (LOK .AND. NBLOCK .NE. NTPERR) THEN

         IF (LGOT(NBLOCK)) THEN
            CALL FLAGLN(PRGSTR)
            PRINT 550,TPNAME(NBLOCK)
            STOP
         ENDIF

         IF (NBLOCK .EQ. NATNAM) THEN
            CALL GATYPE(IUNIT,PRGSTR)
         ELSEIF (NBLOCK .EQ. NRSNAM) THEN
            CALL GRES(IUNIT,PRGSTR)
         ELSEIF (NBLOCK .EQ. NSOLAT) THEN
            CALL GSLUAT(IUNIT,PRGSTR)
C bonds
         ELSEIF (NBLOCK .EQ. NBNDTY) THEN
            CALL GBTY(IUNIT,PRGSTR)
         ELSEIF (NBLOCK .EQ. NBNDHY) THEN
            CALL CHKGOT(PRGSTR,NBNDHY,NBNDTY)
            CALL CHKGOT(PRGSTR,NBNDHY,NSOLAT)
            CALL GBND(IUNIT,PRGSTR,NBONH,MAXBNH,NRP,NBTY,IBH,JBH,ICBH)
         ELSEIF (NBLOCK .EQ. NBND) THEN
            CALL CHKGOT(PRGSTR,NBND,NBNDTY)
            CALL CHKGOT(PRGSTR,NBND,NSOLAT)
            CALL GBND(IUNIT,PRGSTR,NBON,MAXBON,NRP,NBTY,IB,JB,ICB)
C bond angles
         ELSEIF (NBLOCK .EQ. NBATY) THEN
            CALL GBANT(IUNIT,PRGSTR)
         ELSEIF (NBLOCK .EQ. NBAHY ) THEN
            CALL CHKGOT(PRGSTR,NBAHY,NBATY)
            CALL CHKGOT(PRGSTR,NBAHY,NSOLAT)
            CALL GBAN(IUNIT,PRGSTR,NTHEH,MXQHEH,
     $           NRP,NTTY,ITH,JTH,KTH,ICTH)
         ELSEIF (NBLOCK .EQ. NBAN ) THEN
            CALL CHKGOT(PRGSTR,NBAN,NBATY)
            CALL CHKGOT(PRGSTR,NBAN,NSOLAT)
            CALL GBAN(IUNIT,PRGSTR,NTHE,MAXTHE,
     $           NRP,NTTY,IT,JT,KT,ICT)
C improper dihedrals
         ELSEIF (NBLOCK .EQ. NIMDTY) THEN
            CALL GIDTY(IUNIT,PRGSTR)
         ELSEIF (NBLOCK .EQ. NIMDH) THEN
            CALL CHKGOT(PRGSTR,NIMDH,NIMDTY)
            CALL CHKGOT(PRGSTR,NIMDH,NSOLAT)
            CALL GID(IUNIT,PRGSTR,
     $           NQHIH,MAXHIH,NRP,NQTY,IQH,JQH,KQH,LQH,ICQH)
         ELSEIF (NBLOCK .EQ. NIMD) THEN
            CALL CHKGOT(PRGSTR,NIMD,NIMDTY)
            CALL CHKGOT(PRGSTR,NIMD,NSOLAT)
            CALL GID(IUNIT,PRGSTR,
     $           NQHI,MAXQHI,NRP,NQTY,IQ,JQ,KQ,LQ,ICQ)
C dihedrals
         ELSEIF (NBLOCK .EQ.  NDITY) THEN
            CALL GDITY(IUNIT,PRGSTR)
         ELSEIF (NBLOCK .EQ.  NDIHY) THEN
            CALL CHKGOT(PRGSTR,NDIHY,NDITY)
            CALL CHKGOT(PRGSTR,NDIHY,NSOLAT)
            CALL GDI(IUNIT,PRGSTR,
     $           NPHIH,MXPHIH,NRP,NPTY,IPH,JPH,KPH,LPH,ICPH)
         ELSEIF (NBLOCK .EQ.  NDIN) THEN
            CALL CHKGOT(PRGSTR,NDIN,NDITY)
            CALL CHKGOT(PRGSTR,NDIN,NSOLAT)
            CALL GDI(IUNIT,PRGSTR,
     $           NPHI,MAXPHI,NRP,NPTY,IP,JP,KP,LP,ICP)
C LJ terms
         ELSEIF (NBLOCK .EQ. NLJBLK) THEN
            CALL CHKGOT(PRGSTR,NLJBLK,NATNAM)
            CALL GLJ(IUNIT,PRGSTR)
C 
         ELSEIF (NBLOCK .EQ.  NSLVBK) THEN
            CALL GSOLV(IUNIT,PRGSTR)
         ELSEIF (NBLOCK .EQ.  NSLVCN) THEN
            CALL CHKGOT(PRGSTR,NSLVCN,NSLVBK)
            CALL GSLVCN(IUNIT,PRGSTR)
C path integral
         ELSEIF (NBLOCK .EQ. NPITBL) THEN
            CALL GPITBL(IUNIT,PRGSTR)
         ELSE
C unknown block
            CALL FLAGLN(PRGSTR)
            PRINT 560,FIOLIN
            STOP
         ENDIF
C 
         LGOT(NBLOCK) = .TRUE.

         CALL GETLN(IUNIT,LOK)
         NBLOCK = NDXTOP(FIOLIN)
         GOTO 10
C end of while loop
      ENDIF

      IF (LOK .AND. NBLOCK .EQ. NTPERR) THEN
C we exited the loop because of an unknown block type, not
C becase we reached the end of file.
         CALL FLAGLN(PRGSTR)
         PRINT 560,FIOLIN
         STOP
      ENDIF

C check we have read all the blocks we need
      LGTALL = .TRUE.
      DO I=MINTPB,MAXTPB
         IF (LNEED(I) .AND. .NOT. LGOT(I)) THEN
            PRINT 570,PRGSTR,TPNAME(I)
            LGTALL = .FALSE.
         ENDIF
      ENDDO

      IF (.NOT. LGTALL) THEN
         STOP
      ENDIF

      LREDTO = .TRUE.

C formats
 550  FORMAT (1X,'blocktype',A,'encountered twice!')
 560  FORMAT (1X,'unknown block type: ',A)
 570  FORMAT (1X,A6,': block ',A,' missing on topology file!')
C end GTTOPO
      END

! ==================================================================
      LOGICAL FUNCTION LTOPOK()
! ==--------------------------------------------------------------==
C     LOGICAL FUNCTION LTOPOK()
C
C     LTOPOK returns .TRUE. if a topology files has been
C     read in successfully. Otherwise .FALSE. is returned.
!    Revision cmb Tokyo/Strasbourg june 2014
! ==--------------------------------------------------------------==
      LOGICAL LREDTO
      COMMON /TOPOK/LREDTO
      LTOPOK = LREDTO
      END

      BLOCKDATA TOPBLK
      INCLUDE 'toposz.h'
      INCLUDE 'topoar.h'
      INCLUDE 'topblock.h'

C     LREDTO is initialized to .FALSE. here
C     and set to .TRUE. after a topology is
C     succesfully read in by RDTOPO.
      LOGICAL LREDTO
      COMMON /TOPOK/LREDTO

      DATA LREDTO /.FALSE./

C init the block names
      DATA TPNAME(NPTITL) / 'TITLE'/
      DATA  LNEED(NPTITL) / .TRUE./

      DATA TPNAME(NATNAM) / 'ATOMTYPENAME'/
      DATA  LNEED(NATNAM) / .FALSE./

      DATA TPNAME(NRSNAM) / 'RESNAME'/
      DATA  LNEED(NRSNAM) / .FALSE./

      DATA TPNAME(NSOLAT) / 'SOLUTEATOM'/
      DATA  LNEED(NSOLAT) / .FALSE./

      DATA TPNAME(NBNDTY) / 'BONDTYPE'/
      DATA  LNEED(NBNDTY) / .FALSE. /

      DATA TPNAME(NBNDHY) / 'BONDH'/
      DATA  LNEED(NBNDHY) / .FALSE./

      DATA TPNAME(NBND)   / 'BOND'/
      DATA  LNEED(NBND)   / .FALSE./

      DATA TPNAME(NBATY)  / 'BONDANGLETYPE'/
      DATA  LNEED(NBATY)  / .FALSE./

      DATA TPNAME(NBAHY)  / 'BONDANGLEH'/
      DATA  LNEED(NBAHY)  / .FALSE./

      DATA TPNAME(NBAN)   / 'BONDANGLE'/
      DATA  LNEED(NBAN)   / .FALSE./

      DATA TPNAME(NIMDTY) / 'IMPDIHEDRALTYPE'/
      DATA  LNEED(NIMDTY) / .FALSE./

      DATA TPNAME(NIMDH)  / 'IMPDIHEDRALH'/
      DATA  LNEED(NIMDH)  / .FALSE./

      DATA TPNAME(NIMD)   / 'IMPDIHEDRAL'/
      DATA  LNEED(NIMD)   / .FALSE./

      DATA TPNAME(NDITY)  / 'DIHEDRALTYPE'/
      DATA  LNEED(NDITY)  / .FALSE./

      DATA TPNAME(NDIHY)  / 'DIHEDRALH'/
      DATA  LNEED(NDIHY)  / .FALSE./

      DATA TPNAME(NDIN)   / 'DIHEDRAL'/
      DATA  LNEED(NDIN)   / .FALSE./

      DATA TPNAME(NLJBLK) / 'LJPARAMETERS'/
      DATA  LNEED(NLJBLK) / .FALSE./

      DATA TPNAME(NTPEND) / 'END'/
      DATA  LNEED(NTPEND) / .FALSE./

      DATA TPNAME(NTPUNT) / 'TOPPHYSCON'/
      DATA  LNEED(NTPUNT) / .TRUE./

      DATA TPNAME(NSLVBK) / 'SOLVENTATOM'/
      DATA  LNEED(NSLVBK) / .FALSE./

      DATA TPNAME(NSLVCN) / 'SOLVENTCONSTR'/
      DATA  LNEED(NSLVCN) / .FALSE./

      DATA TPNAME(NTPVER) / 'TOPVERSION'/
      DATA  LNEED(NTPVER) / .TRUE./

      DATA TPNAME(NPITBL) / 'PATHINTSPEC'/
      DATA  LNEED(NPITBL) / .FALSE./

C init the topology itself
      DATA    NRP /0/
      DATA  NRATT,NRATT2/0,0/
      DATA  NRAA2 /0/
      DATA   NBTY,NBON,NBONH /0,0,0/
      DATA   NTTY,NTHE,NTHEH /0,0,0/
      DATA   NQTY,NQHI,NQHIH /0,0,0/
      DATA   NPTY,NPHI,NPHIH /0,0,0/
      DATA   NCAG,NAEX,NAEX14/0,0,0/
      DATA   NRAM,NCONS      /0,0/

      DATA   FPEPSI,HBAR       /0.0,0.0/
      DATA   NPID /0/

C end blockdata
      END

! ==================================================================
      SUBROUTINE GUNT(IUNIT,PRGSTR,FPEPSI,HBAR)
! ==--------------------------------------------------------------==
C     Subroutine GUNT reads in a TOPPHYSCON block from
C     a molecular topology file.
C     The block contains the variables L<FPEPSI> and L<HBAR>.
C     
C     GUNT is called in subroutine L<RGTTOP>.
!    Revision cmb Tokyo/Strasbourg june 2014
! ==--------------------------------------------------------------==
C args
      INTEGER IUNIT
      CHARACTER PRGSTR*(*)
      REAL FPEPSI,HBAR
C local vars
      CHARACTER STR*(6)
C begin
      STR = 'FPEPSI'
      CALL GIMME(IUNIT,PRGSTR,STR)
      CALL CHPREL(FPEPSI,PRGSTR,STR)
      IF (FPEPSI .LE. 0.0) THEN
         CALL FLAGLN(PRGSTR)
         PRINT *,'invalid value for FPEPSI = ',FPEPSI
         STOP
      ENDIF
      CALL CHKMTY(PRGSTR,STR)

      STR = 'HBAR'
      CALL GIMME(IUNIT,PRGSTR,STR)
      CALL CHPREL(HBAR,PRGSTR,STR)
      IF (HBAR .LE. 0.0) THEN
         CALL FLAGLN(PRGSTR)
         PRINT *,'invalid value for HBAR = ',HBAR
         STOP
      ENDIF
      CALL CHKMTY(PRGSTR,STR)

      CALL CHKEND(IUNIT,PRGSTR)
C end gunt
      END

! ==================================================================
      SUBROUTINE CHKGOT(PRGSTR,NDXRED,NDXCHK)
! ==--------------------------------------------------------------==
C     CHKGOT checks whether a topology block NDXCHK has
C     been read in. If not, a message is written to screen
C     and TOP is called.
C     This routine is used to check that a block NDXCHK
C     has been read in before NDXRED has been read in.
C     It is called by L<RGTTOP> before reading in block
C     NDXRED.
!    Revision cmb Tokyo/Strasbourg june 2014
! ==--------------------------------------------------------------==
      INCLUDE 'topblock.h'
C args
      CHARACTER PRGSTR*(*)
      INTEGER NDXRED,NDXCHK
C begin
      IF (.NOT. LGOT(NDXCHK)) THEN
         CALL FLAGLN(PRGSTR)
         PRINT 500,TPNAME(NDXCHK),TPNAME(NDXRED)
         STOP
      ENDIF
 500  FORMAT (1X,'block ',A8,' must occur before block ',A8)
C end CHKGOT
      END

! ==================================================================
      INTEGER FUNCTION NDXTOP(BTN)
! ==--------------------------------------------------------------==
C args
      CHARACTER *(*)BTN
C includes
      INCLUDE 'topblock.h'
C local vars
      INTEGER I
C begin
      I = MINTPB
C look for BTN in table using while loop
 10   IF (I .LE. MAXTPB .AND. TPNAME(I) .NE. BTN) THEN
         I = I + 1
         GOTO 10
      ENDIF
C     end of while loop
C
      IF (I .GT. MAXTPB) THEN
         NDXTOP = NTPERR
C         PRINT *,'GOT ERROR'
      ELSE
C         PRINT *,'GOT ',TPNAME(I)
         NDXTOP = I
      ENDIF
C end NDXTOP
      END


! ==================================================================
      SUBROUTINE CHKEND(IUNIT,PRGSTR)
! ==--------------------------------------------------------------==
!    Check if the end of the file has been reached
!    (correctly or not)
!    Revision cmb Tokyo/Strasbourg june 2014
! ==--------------------------------------------------------------==
      INCLUDE 'topblock.h'
      INCLUDE 'fileio.h'
C args
      INTEGER IUNIT
      CHARACTER PRGSTR*(*)
C externals
      INTEGER NDXTOP
      EXTERNAL NDXTOP
C local vars
      LOGICAL LOK
      INTEGER NBLOCK
C begin
      CALL GETLN(IUNIT,LOK)
      NBLOCK = NDXTOP(FIOLIN)
      IF (.NOT. LOK .OR. NBLOCK .NE. NTPEND) THEN
         CALL FLAGLN(PRGSTR)
         PRINT '(A20)',FIOLIN
         PRINT *,'END expected'
         STOP
      ENDIF
C end CHKEND
      END

! ==================================================================
      SUBROUTINE GATYPE(IUNIT,PRGSTR)
! ==--------------------------------------------------------------==
!    Check atom type (chemical species) and related name/label
!    Revision cmb Tokyo/Strasbourg june 2014
! ==--------------------------------------------------------------==
      INCLUDE 'toposz.h'
      INCLUDE 'topoar.h'
      INCLUDE 'fileio.h'
      INCLUDE 'formats.h'
C args
      INTEGER IUNIT
      CHARACTER PRGSTR*(*)
C local vars
      INTEGER K
      CHARACTER STR*(14)
C begin
C nratt
      STR = 'NRATT'
      CALL GIMME(IUNIT,PRGSTR,STR)
      CALL CHPINT(NRATT,PRGSTR,STR)
      IF (NRATT .LT. 0 .OR. NRATT .GT. MAXATT) THEN
         CALL FLAGLN(PRGSTR)
         PRINT FMNII,STR,NRATT
         PRINT FMBTII,0,STR,'MAXATT',MAXATT
         STOP
      ENDIF
      CALL CHKMTY(PRGSTR,STR)

C type
      STR = 'ATOM TYPE NAME'
      DO 10 K=1,NRATT
         CALL GIMME(IUNIT,PRGSTR,STR)
         CALL CHPSTR(FFTYPE(K),PRGSTR,STR)
         CALL CHKMTY(PRGSTR,STR)
 10   CONTINUE

      CALL CHKEND(IUNIT,PRGSTR)
C end GATYPE
      END

! ==================================================================
      SUBROUTINE GRES(IUNIT,PRGSTR)
! ==--------------------------------------------------------------==
!    Revision cmb Tokyo/Strasbourg june 2014
! ==--------------------------------------------------------------==
      INCLUDE 'toposz.h'
      INCLUDE 'topoar.h'
      INCLUDE 'fileio.h'
      INCLUDE 'formats.h'
C args
      INTEGER IUNIT
      CHARACTER PRGSTR*(*)
C local vars
      INTEGER K
      CHARACTER STR*(5)
C begin
C nraa2
      STR = 'NRAA2'
      CALL GIMME(IUNIT,PRGSTR,STR)
      CALL CHPINT(NRAA2,PRGSTR,STR)
      IF (NRAA2 .LT. 0 .OR. NRAA2 .GT. MAXAA2) THEN
         CALL FLAGLN(PRGSTR)
         PRINT FMNII,STR,NRAA2
         PRINT FMBTII,0,STR,'MAXAA2',MAXAA2
         STOP
      ENDIF
      CALL CHKMTY(PRGSTR,STR)

C aanm
      STR = 'AANM'
      DO 12 K=1,NRAA2
         CALL GIMME(IUNIT,PRGSTR,STR)
         CALL CHPSTR(AANM(K),PRGSTR,STR)
         CALL CHKMTY(PRGSTR,STR)
 12   CONTINUE

      CALL CHKEND(IUNIT,PRGSTR)
C end GRES
      END

! ==================================================================
      SUBROUTINE GSLUAT(IUNIT,PRGSTR)
! ==--------------------------------------------------------------==
!    Revision cmb Tokyo/Strasbourg june 2014
! ==--------------------------------------------------------------==
      INCLUDE 'toposz.h'
      INCLUDE 'topoar.h'
      INCLUDE 'fileio.h'
      INCLUDE 'formats.h'
C args
      INTEGER IUNIT
      CHARACTER PRGSTR*(*)
C local vars
      INTEGER K,KOLD,MRESOL,II,JJ,ITMP,IATCG,ICHKEXC
      LOGICAL LERR
      CHARACTER STR*(46)
C external
      LOGICAL LISMTY
      EXTERNAL LISMTY
C begin
      LERR = .FALSE.
      NCAG = 0
      NAEX = 0
      NAEX14 = 0
C nrp
      STR = 'NRP'
      CALL GIMME(IUNIT,PRGSTR,STR)
      CALL CHPINT(NRP,PRGSTR,STR)
      IF (NRP .LT. 0 .OR. NRP .GT. MAXNRP) THEN
         CALL FLAGLN(PRGSTR)
         PRINT FMNII,STR,NRP
         PRINT FMBTII,0,STR,'MAXNRP',MAXNRP
         STOP
      ENDIF
      CALL CHKMTY(PRGSTR,STR)

C data consists of two lines:
C     k,mres,panm,iac,mass,cg,tcgc,ine, jsne...
C     ine14, jsne14.....
C
      KOLD = 0
      MRESOL = 0

      DO 14 II=1,NRP
         STR = 'ATNM, MRES, PANM,...'
         CALL GIMME(IUNIT,PRGSTR,STR)
C K
         CALL CHPINT(K,PRGSTR,STR)
C check range
         IF (K .LT. 1 .OR. K .GT. NRP) THEN
            CALL FLAGLN(PRGSTR)
            STR = 'atom sequence number'
            PRINT FMNII,STR,K
            PRINT FMBTII,1,STR,'NRP',NRP
            STOP
         ENDIF
C check order
         IF (KOLD .GT. K) THEN
            CALL FLAGLN(PRGSTR)
            PRINT *,'atoms must be specified in order!'
            PRINT *,'BUT: atom',KOLD,' specified BEFORE atom',K
            STOP
         ENDIF
         IF (KOLD .EQ. K) THEN
            CALL FLAGLN(PRGSTR)
            PRINT *,'atom number',K,' specified twice!'
            STOP
         ENDIF

         CALL CHPINT(MRES(K),PRGSTR,STR)
C check range
         IF (MRES(K) .LT. 1 .OR. MRES(K) .GT. NRAA2) THEN
            CALL FLAGLN(PRGSTR)
            PRINT *,'illegal residue number for atom',K
            PRINT FMBTII,1,'MRES','NRAA2',NRAA2
            PRINT *,'but read MRES = ',MRES(K)
            STOP
         ENDIF

C check order
         IF (MRES(K) .LT. MRESOL) THEN
            CALL FLAGLN(PRGSTR)
            PRINT *,'illegal ordering in residue number before atom',K
            PRINT *,'residue numbers must be a monotonically increasing'
            PRINT *,'function of atom sequence numbers'
            STOP
         ENDIF
         MRESOL = MRES(K)
C panm
         CALL CHPSTR(PANM(K),PRGSTR,STR)
C iac
         CALL CHPINT(IAC(K),PRGSTR,STR)
         IF (IAC(K) .LT. 1 .OR. IAC(K) .GT. NRATT) THEN
            CALL FLAGLN(PRGSTR)
            PRINT *,'illegal atom type in IAC!'
            PRINT *,'IAC(',K,') = ',IAC(K)
            PRINT FMBTII,1,'IAC(I)','NRATT',NRATT
            STOP
         ENDIF

C get mass and invert
         CALL CHPREL(WMAS(K),PRGSTR,STR)
         IF (WMAS(K) .LE. 0.0E+00) THEN
            CALL FLAGLN(PRGSTR)
            PRINT *,'atom',K,' has zero or negative mass!'
            STOP
         ENDIF
         WINV(K) = 1.0E+00/WMAS(K)

C get charge
         CALL CHPREL(CG(K),PRGSTR,STR)

C iatcg --> set NCAG and INC
         CALL CHPINT(IATCG,PRGSTR,STR)
         IF (IATCG .NE. 0 .AND. IATCG .NE. 1)THEN
            CALL FLAGLN(PRGSTR)
            PRINT FMNII,'IATCG',IATCG
            PRINT FMBOO,'IATCG'
            STOP
         ENDIF
         IF (IATCG .EQ. 1) THEN
C this atom terminates a charge group
            NCAG = NCAG+1
            IF (NCAG .GT. MAXCAG) THEN
               CALL FLAGLN(PRGSTR)
               PRINT *,'too many charge groups specified!'
               PRINT *,'maximum allowed: MAXCAG = ',MAXCAG
               STOP
            ENDIF
            INC(NCAG) = K
         ENDIF

C ine and exclusions --> set INE, KNE and JSNE
C ----this here as a reminder
C EXCLUSIONS J OF ATOM I ARE POSITIONED AT    
C                     JSNE(KNE(I)+1),...,JSNE(KNE(I)+INE(I));
C                     IT IS ASSUMED THAT I.LT.J
C ----
         CALL CHPINT(INE(K),PRGSTR,STR)
         IF (INE(K) .LT. 0 .OR. INE(K) .GE. NRP) THEN
            CALL FLAGLN(PRGSTR)
            PRINT 500,K
            PRINT *,'number of exclusions read =',INE(K)
            PRINT FMBTII,1,'INE(I)','NRP',NRP
            PRINT *
            LERR = .TRUE.
         ENDIF

         KNE(K) = NAEX
         DO JJ=1,INE(K)
            CALL CHPINT(ITMP,PRGSTR,STR)
C check range
            IF (ITMP .LT. 1 .OR. ITMP .GT. NRP) THEN
               CALL FLAGLN(PRGSTR)
               PRINT 500,K
               STR = 'atom sequence number J'
               PRINT FMNII,STR,ITMP
               PRINT FMBTII,1,STR,'NRP',NRP
               LERR = .TRUE.
            ENDIF

C make sure ITMP > K
            IF (ITMP .LE. K) THEN
               CALL FLAGLN(PRGSTR)
               PRINT 500,K
               PRINT *,'atom sequence number J of excluded atom'
               PRINT *,'must be larger than I'
               PRINT *,'but read J = ',ITMP
               PRINT *
               LERR = .TRUE.
            ENDIF

C make sure we have exclusions in ascending order
            IF (JJ .GT. 1 .AND. JSNE(NAEX) .GE. ITMP) THEN
               IF (JSNE(NAEX) .EQ. ITMP) THEN
                  CALL FLAGLN(PRGSTR)
                  PRINT 500,K
                  PRINT *,'exclusion specified twice:',ITMP
               ELSE
                  CALL FLAGLN(PRGSTR)
                  PRINT 500,K  
                  PRINT *,'exclusion list must be in ascending order!'
                  PRINT *,'but read ',JSNE(NAEX)
                  PRINT *,'before   ',ITMP
               ENDIF
               LERR = .TRUE.
            ENDIF

C add to list
            NAEX = NAEX + 1
            IF (NAEX .GT. MAXAEX) THEN
               CALL FLAGLN(PRGSTR)
               PRINT *,'too many exclusions specified'
               PRINT *,'maximum allowed: MAXAEX = ',MAXAEX
               STOP
            ENDIF
            JSNE(NAEX) = ITMP
         ENDDO

C make sure the line is empty...
         IF (.NOT. LISMTY()) THEN
            CALL FLAGLN(PRGSTR)
            PRINT 500,K
            PRINT *,'inconsistent number of excluded atoms'
            PRINT *,'number specified :',INE(K)
            PRINT *,'but more data detected on the line'
            LERR = .TRUE.
         ENDIF

C on the next line: 1-4 interactions
         STR = '1-4 interactions: INE14,KNE14(...)'
         CALL GIMME(IUNIT,PRGSTR,STR)

         CALL CHPINT(INE14(K),PRGSTR,STR)
         IF (INE14(K) .LT. 0 .OR. INE14(K) .GE. NRP) THEN
            CALL FLAGLN(PRGSTR)
            PRINT 510,K
            PRINT *,'number of 1-4 interactions read ',INE14(K)
            PRINT FMBTII,1,'INE14(I)','NRP',NRP
            PRINT *
            LERR = .TRUE.
         ENDIF

         KNE14(K) = NAEX14
         DO 30 JJ=1,INE14(K)
            CALL CHPINT(ITMP,PRGSTR,STR)
C check range
            IF (ITMP .LT. 1 .OR. ITMP .GT. NRP) THEN
               CALL FLAGLN(PRGSTR)
               PRINT 510,K
               PRINT *,'read illegal atom sequence number J = ',ITMP
               PRINT FMBTII,1,'J','NRP',NRP
               LERR = .TRUE.
            ENDIF

C make sure ITMP > K
            IF (ITMP .LE. K) THEN
               CALL FLAGLN(PRGSTR)
               PRINT 510,K
               PRINT *,'atom sequence number J must be larger than I'
               PRINT *,'but read J = ',ITMP
               LERR = .TRUE.
            ENDIF

C make sure we have 1-4 interactions in ascending order
            IF (JJ .GT. 1 .AND. JSNE14(NAEX14) .GE. ITMP) THEN
               IF (JSNE14(NAEX14) .EQ. ITMP) THEN
                  CALL FLAGLN(PRGSTR)
                  PRINT 510,K
                  PRINT *,'1-4 interaction specified twice:',ITMP
               ELSE
                  CALL FLAGLN(PRGSTR)
                  PRINT 510,K
                  PRINT *,'1-4 interaction list must be ',
     $                 'in ascending order !'
                  PRINT *,'but read ',JSNE14(NAEX14)
                  PRINT *,'before   ',ITMP
               ENDIF
               LERR = .TRUE.
            ENDIF

C make sure J is not an excluded atom as well
            ICHKEXC = KNE(K)
 117        IF (ICHKEXC .LT. NAEX) THEN
               IF (JSNE(ICHKEXC+1) .EQ. ITMP) THEN
                  CALL FLAGLN(PRGSTR)
                  PRINT 510,K
                  PRINT *,'Atom J is in the exclusion list AND',
     $                 ' in the 1-4 interaction list'
                  PRINT *,'Atom J =',ITMP
                  LERR = .TRUE.
               ENDIF
               ICHKEXC = ICHKEXC + 1
               GOTO 117
            ENDIF

            NAEX14 = NAEX14 + 1
            IF (NAEX14 .GT. MXEX14) THEN
               CALL FLAGLN(PRGSTR)
               PRINT *,'too many 1-4 interactions specified'
               PRINT *,'maximum allowed: MXEX14= ',MXEX14
               STOP
            ENDIF
            JSNE14(NAEX14) = ITMP

 30      CONTINUE

C make sure the line is empty...
         IF (.NOT. LISMTY()) THEN
            CALL FLAGLN(PRGSTR)
            PRINT 510,K
            PRINT *,'inconsistent number of 1-4 atoms'
            PRINT *,'number specified :',INE14(K)
            PRINT *,'but more data detected on the line'
            LERR = .TRUE.
         ENDIF

         KOLD = K
 14   CONTINUE

      IF (LERR) THEN
         STOP
      ENDIF

C check NCAG 
      IF (NCAG .GT. 0) THEN
         IF (INC(NCAG) .NE. NRP) THEN
            CALL FLAGLN(PRGSTR)
            PRINT *,'The last atom in the solute must terminate'
            PRINT *,'a charge group (CGC = 1)'
            STOP
         ENDIF
      ENDIF

      CALL CHKEND(IUNIT,PRGSTR)

 500  FORMAT (1X,'error in exclusion list for atom I =',I5)
 510  FORMAT (1X,'error in 1-4 interaction list for atom I =',I5)
C end GSLUAT
      END

! ==================================================================
      SUBROUTINE GBTY(IUNIT,PRGSTR)
! ==--------------------------------------------------------------==
!    Revision cmb Tokyo/Strasbourg june 2014
! ==--------------------------------------------------------------==
      INCLUDE 'toposz.h'
      INCLUDE 'topoar.h'
      INCLUDE 'formats.h'
C args
      INTEGER IUNIT
      CHARACTER PRGSTR*(*)
C local vars
      INTEGER K
      CHARACTER STR*(5)
C begin
C nbty
      STR = 'NBTY'
      CALL GIMME(IUNIT,PRGSTR,STR)
      CALL CHPINT(NBTY,PRGSTR,STR)
      IF (NBTY .LT. 0 .OR. NBTY .GT. MAXNBT) THEN
         CALL FLAGLN(PRGSTR)
         PRINT FMNII,STR,NBTY
         PRINT FMBTII,0,STR,'MAXNBT',MAXNBT
         STOP
      ENDIF
      CALL CHKMTY(PRGSTR,STR)

C cb, b0
      STR = 'CB,B0'
      DO 18 K=1,NBTY
         CALL GIMME(IUNIT,PRGSTR,STR)
         CALL CHPREL(CB(K),PRGSTR,STR)
         CALL CHPREL(B0(K),PRGSTR,STR)
         CALL CHKMTY(PRGSTR,STR)
         IF (B0(K) .LT. 0.0) THEN
            CALL FLAGLN(PRGSTR)
            PRINT FMNIR,'B0',B0
            PRINT FMGER,'B0'
            STOP
         ENDIF
 18   CONTINUE

      CALL CHKEND(IUNIT,PRGSTR)
C end GBTY
      END

! ==================================================================
      SUBROUTINE GBND(IUNIT,PRGSTR,NBONH,MAXBNH,NRP,NBTY,IBH,JBH,ICBH)
! ==--------------------------------------------------------------==
!    Revision cmb Tokyo/Strasbourg june 2014
! ==--------------------------------------------------------------==
C includes
      INCLUDE 'formats.h'
C args
      INTEGER IUNIT
      CHARACTER PRGSTR*(*)
      INTEGER NBONH,MAXBNH,NRP,NBTY
      INTEGER IBH(MAXBNH),JBH(MAXBNH),ICBH(MAXBNH)
C local vars
      INTEGER I
      LOGICAL LERR
      CHARACTER STR*(40)
C begin
      LERR = .FALSE.
C nbonh
      STR = 'NBON(H)'
      CALL GIMME(IUNIT,PRGSTR,STR)
      CALL CHPINT(NBONH,PRGSTR,STR)
      IF (NBONH .LT. 0 .OR. NBONH .GT. MAXBNH) THEN
         CALL FLAGLN(PRGSTR)
         PRINT FMNII,STR
         PRINT FMBTII,0,STR,'MAX',MAXBNH
         STOP
      ENDIF
      CALL CHKMTY(PRGSTR,STR)

C ibh,jbh,icbh
      STR = 'Atom seq num, Atom seq num, bond type'
      DO 19 I=1,NBONH
         CALL GIMME(IUNIT,PRGSTR,STR)
         CALL CHPINT(IBH(I),PRGSTR,STR)
         CALL CHPINT(JBH(I),PRGSTR,STR)
         CALL CHPINT(ICBH(I),PRGSTR,STR)
         CALL CHKMTY(PRGSTR,STR)

C check ibh,jbh,icbh
         IF (IBH(I) .LT. 1 .OR. IBH(I) .GT. NRP .OR.
     $        JBH(I) .LT. 1 .OR. JBH(I) .GT. NRP) THEN
            LERR = .TRUE.
            CALL FLAGLN(PRGSTR)
            PRINT *,'illegal atom sequence number in bond !'
            PRINT *,'read for atom I:',IBH(I)
            PRINT *,'read for atom J:',JBH(I)
            PRINT FMBTII,1,'I,J','NRP',NRP
            PRINT *
         ENDIF

         IF (ICBH(I) .LT. 1 .OR. ICBH(I) .GT. NBTY) THEN
            LERR = .TRUE.
            CALL FLAGLN(PRGSTR)
            STR = 'bond type'
            PRINT FMNII,STR,ICBH(I)
            PRINT FMBTII,1,STR,'NBTY',NBTY
         ENDIF
 19   CONTINUE

      CALL CHKEND(IUNIT,PRGSTR)

      IF (LERR) THEN
         STOP
      ENDIF
C end GBND
      END

! ==================================================================
      SUBROUTINE GBANT(IUNIT,PRGSTR)
! ==--------------------------------------------------------------==
!    Revision cmb Tokyo/Strasbourg june 2014
! ==--------------------------------------------------------------==
      INCLUDE 'toposz.h'
      INCLUDE 'topoar.h'
      INCLUDE 'formats.h'
C args
      INTEGER IUNIT
      CHARACTER PRGSTR*(*)
C local vars
      INTEGER K
      CHARACTER STR*(5)
C begin
C ntty
      STR = 'NTTY'
      CALL GIMME(IUNIT,PRGSTR,STR)
      CALL CHPINT(NTTY,PRGSTR,STR)
      IF ( NTTY .LT. 0 .OR. NTTY .GT. MAXTTY) THEN
         CALL FLAGLN(PRGSTR)
         PRINT FMNII,STR,NTTY
         PRINT FMBTII,0,STR,'MAXTTY',MAXTTY
         STOP
      ENDIF
      CALL CHKMTY(PRGSTR,STR)

C ct,t0
      STR = 'CT,T0'
      DO K=1,NTTY
         CALL GIMME(IUNIT,PRGSTR,STR)
         CALL CHPREL(CT(K),PRGSTR,STR)
         CALL CHPREL(T0(K),PRGSTR,STR)
         CALL CHKMTY(PRGSTR,STR)
      ENDDO

      CALL CHKEND(IUNIT,PRGSTR)
C end GBANT
      END

! ==================================================================
      SUBROUTINE GBAN(IUNIT,PRGSTR,NTHEH,MXQHEH,NRP
     .               ,NTTY,ITH,JTH,KTH,ICTH)
! ==--------------------------------------------------------------==
!    Revision cmb Tokyo/Strasbourg june 2014
! ==--------------------------------------------------------------==
C includes
      INCLUDE 'formats.h'
C args
      INTEGER IUNIT
      CHARACTER PRGSTR*(*)
C
      INTEGER NTHEH,MXQHEH,NRP,NTTY
      INTEGER ITH(MXQHEH),JTH(MXQHEH)
      INTEGER KTH(MXQHEH),ICTH(MXQHEH)
C local vars
      INTEGER I
      CHARACTER STR*(20)
      LOGICAL LERR
C begin
      LERR = .FALSE.
C ntheh
      STR = 'num of bond angles'
      CALL GIMME(IUNIT,PRGSTR,STR)
      CALL CHPINT(NTHEH,PRGSTR,STR)
      IF (NTHEH .LT. 0 .OR. NTHEH .GT. MXQHEH) THEN
         CALL FLAGLN(PRGSTR)
         PRINT FMNII,STR,NTHEH
         PRINT FMBTII,0,STR,'MAXQHE(H)',MXQHEH
         STOP
      ENDIF
      CALL CHKMTY(PRGSTR,STR)

C ith,jth,kth,icth
      STR = 'ITH,JTH,KTH,ICTH'
      DO I=1,NTHEH
         CALL GIMME(IUNIT,PRGSTR,STR)
         CALL CHPINT(ITH(I),PRGSTR,STR)
         CALL CHPINT(JTH(I),PRGSTR,STR)
         CALL CHPINT(KTH(I),PRGSTR,STR)
         CALL CHPINT(ICTH(I),PRGSTR,STR)
         CALL CHKMTY(PRGSTR,STR)

C check ith,jth,kth,icth
         IF (ITH(I) .LT. 1 .OR. ITH(I) .GT. NRP) THEN
            LERR = .TRUE.
            CALL FLAGLN(PRGSTR)
            PRINT 500,I
            PRINT 510,ITH(I)
            PRINT 520,NRP
         ENDIF

         IF (JTH(I) .LT. 1 .OR. JTH(I) .GT. NRP) THEN
            LERR = .TRUE.
            CALL FLAGLN(PRGSTR)
            PRINT 500,I
            PRINT 510,JTH(I)
            PRINT 520,NRP
         ENDIF

         IF (KTH(I) .LT. 1 .OR. KTH(I) .GT. NRP) THEN
            LERR = .TRUE.
            CALL FLAGLN(PRGSTR)
            PRINT 500
            PRINT 510,KTH(I)
            PRINT 520,NRP
         ENDIF

         IF (ICTH(I) .LT. 1 .OR. ICTH(I) .GT. NTTY) THEN
            LERR = .TRUE.
            CALL FLAGLN(PRGSTR)
            PRINT *,'illegal bond angle type for bond angle',I
            PRINT *,'read BAT:',ICTH(I)
            PRINT FMBTII,1,'BAT','NTTY',NTTY
            PRINT *
         ENDIF

      ENDDO

      CALL CHKEND(IUNIT,PRGSTR)

      IF (LERR) THEN
         STOP
      ENDIF

 500  FORMAT (1X,'atom sequence number out of range for bond angle',I5)
 510  FORMAT (1X,'read ASN:',I5)
 520  FORMAT (1X,'must be 0 < ASN <= NRP = ',I5,//)
C end GBAN
      END

! ==================================================================
      SUBROUTINE GIDTY(IUNIT,PRGSTR)
! ==--------------------------------------------------------------==
!    Revision cmb Tokyo/Strasbourg june 2014
! ==--------------------------------------------------------------==
      INCLUDE 'toposz.h'
      INCLUDE 'topoar.h'
      INCLUDE 'formats.h'
C args
      INTEGER IUNIT
      CHARACTER PRGSTR*(*)
C local vars
      INTEGER K
      CHARACTER STR*(5)
C begin
C nqty
      STR = 'NQTY'
      CALL GIMME(IUNIT,PRGSTR,STR)
      CALL CHPINT(NQTY,PRGSTR,STR)
      IF (NQTY .LT. 0 .OR. NQTY .GT. MAXQTY) THEN
         CALL FLAGLN(PRGSTR)
         PRINT FMNII,STR,NQTY
         PRINT FMBTII,0,STR,'MAXQTY',MAXQTY
         STOP
      ENDIF
      CALL CHKMTY(PRGSTR,STR)

C cq,q0
      STR = 'CQ,Q0'
      DO 52 K=1,NQTY
         CALL GIMME(IUNIT,PRGSTR,STR)
         CALL CHPREL(CQ(K),PRGSTR,STR)
         CALL CHPREL(Q0(K),PRGSTR,STR)
         CALL CHKMTY(PRGSTR,STR)
 52   CONTINUE
      
      CALL CHKEND(IUNIT,PRGSTR)
C end GIDTY
      END

! ==================================================================
      SUBROUTINE GID(IUNIT,PRGSTR,NQHIH,MAXHIH,NRP
     .              ,NQTY,IQH,JQH,KQH,LQH,ICQH)
! ==--------------------------------------------------------------==
!    Revision cmb Tokyo/Strasbourg june 2014
! ==--------------------------------------------------------------==
C includes
      INCLUDE 'formats.h'
C args
      INTEGER IUNIT
      CHARACTER PRGSTR*(*)
C
      INTEGER NQHIH,MAXHIH,NRP,NQTY
      INTEGER IQH(MAXHIH),JQH(MAXHIH),KQH(MAXHIH)
      INTEGER LQH(MAXHIH),ICQH(MAXHIH)
C local vars
      INTEGER I
      CHARACTER STR*(20)
      LOGICAL LERR
C begin
      LERR = .FALSE.
C nqhih
      STR = 'imp. dihedral number'
      CALL GIMME(IUNIT,PRGSTR,STR)
      CALL CHPINT(NQHIH,PRGSTR,STR)
      IF (NQHIH .LT. 0 .OR. NQHIH .GT. MAXHIH) THEN
         CALL FLAGLN(PRGSTR)
         PRINT FMNII,'NQHIH',NQHIH
         PRINT FMBTII,0,'NQHIH','MAXHIH',MAXHIH
         STOP
      ENDIF
      CALL CHKMTY(PRGSTR,STR)

C iqh,jqh,kqh,lqh,icqh
      STR = 'IQH,JQH,KQH,LQH,ICQH'
      DO 54 I=1,NQHIH
         CALL GIMME(IUNIT,PRGSTR,STR)
         CALL CHPINT(IQH(I),PRGSTR,STR)
         CALL CHPINT(JQH(I),PRGSTR,STR)
         CALL CHPINT(KQH(I),PRGSTR,STR)
         CALL CHPINT(LQH(I),PRGSTR,STR)
         CALL CHPINT(ICQH(I),PRGSTR,STR)
         CALL CHKMTY(PRGSTR,STR)

C check iqh,jqh,kqh,lqh,icqh
         IF (IQH(I) .LT. 1 .OR. IQH(I) .GT. NRP) THEN
            LERR = .TRUE.
            CALL FLAGLN(PRGSTR)
            PRINT 500,I
            PRINT 510,IQH(I)
            PRINT 520
         ENDIF

         IF (JQH(I) .LT. 1 .OR. JQH(I) .GT. NRP) THEN
            LERR = .TRUE.
            PRINT 500,I
            PRINT 510,JQH(I)
            PRINT 520
         ENDIF

         IF (KQH(I) .LT. 1 .OR. KQH(I) .GT. NRP) THEN
            LERR = .TRUE.
            CALL FLAGLN(PRGSTR)
            PRINT *,'illegal atom sequence number in KQH!'
            PRINT *,'KQH(',I,') = ',KQH(I)
            PRINT FMBTII,1,'KQH(I)','NRP',NRP
            PRINT *
         ENDIF

         IF (LQH(I) .LT. 1 .OR. LQH(I) .GT. NRP) THEN
            LERR = .TRUE.
            CALL FLAGLN(PRGSTR)
            PRINT 500,I
            PRINT 510,JQH(I)
            PRINT 520
         ENDIF

         IF (ICQH(I) .LT. 1 .OR. ICQH(I) .GT. NQTY) THEN
            LERR = .TRUE.
            CALL FLAGLN(PRGSTR)
            PRINT *,'illegal impr. dihedral type in ',I
            PRINT *,'read IDA:',ICQH(I)
            PRINT FMBTII,1,'IDA','NQTY',NQTY
            PRINT *
         ENDIF

 54   CONTINUE

      CALL CHKEND(IUNIT,PRGSTR)

      IF (LERR) THEN
         STOP
      ENDIF

 500  FORMAT (1X,'atom sequence number out of range for',
     $     'imp. dihedral angle',I5)
 510  FORMAT (1X,'read ASN:',I5)
 520  FORMAT (1X,'must be 0 < ASN <= NRP = ',I5,//)
C end GID
      END

! ==================================================================
      SUBROUTINE GDITY(IUNIT,PRGSTR)
! ==--------------------------------------------------------------==
!    Revision cmb Tokyo/Strasbourg june 2014
! ==--------------------------------------------------------------==
C includes
      INCLUDE 'toposz.h'
      INCLUDE 'topoar.h'
      INCLUDE 'formats.h'
C args
      INTEGER IUNIT
      CHARACTER PRGSTR*(*)
C local vars
      INTEGER K
      CHARACTER STR*(10)
C begin
C npty
      STR = 'NPTY'
      CALL GIMME(IUNIT,PRGSTR,STR)
      CALL CHPINT(NPTY,PRGSTR,STR)
      IF (NPTY .LT. 0 .OR. NPTY .GT. MAXPTY) THEN
         CALL FLAGLN(PRGSTR)
         PRINT FMNII,STR,NPTY
         PRINT FMBTII,0,STR,'MAXPTY',MAXPTY
         STOP
      ENDIF
      CALL CHKMTY(PRGSTR,STR)

C cp,pd,np
      STR = 'CP,PD,NP'
      DO K=1,NPTY
         CALL GIMME(IUNIT,PRGSTR,STR)
         CALL CHPREL(CP(K),PRGSTR,STR)
         CALL CHPREL(PD(K),PRGSTR,STR)
         CALL CHPINT(NP(K),PRGSTR,STR)
         CALL CHKMTY(PRGSTR,STR)

         IF (PD(K) .NE. 1.0 .AND. PD(K) .NE. -1.0) THEN
            CALL FLAGLN(PRGSTR)
            PRINT FMNII,'PD',PD(K)
            PRINT *,'must be -1.0 or +1.0'
            STOP
         ENDIF
      ENDDO

      CALL CHKEND(IUNIT,PRGSTR)
C end GDITY
      END

! ==================================================================
      SUBROUTINE GDI(IUNIT,PRGSTR,NPHIH,MXPHIH,NRP
     .              ,NPTY,IPH,JPH,KPH,LPH,ICPH)
! ==--------------------------------------------------------------==
!    Revision cmb Tokyo/Strasbourg june 2014
! ==--------------------------------------------------------------==
C includes
      INCLUDE 'formats.h'
C args
      INTEGER IUNIT
      CHARACTER PRGSTR*(*)
C local vars
      INTEGER I
      CHARACTER STR*(32)
      LOGICAL LERR
C
      INTEGER NPHIH,MXPHIH,NRP,NPTY
      INTEGER IPH(MXPHIH),JPH(MXPHIH),KPH(MXPHIH)
      INTEGER LPH(MXPHIH),ICPH(MXPHIH)
C begin
C nphih
      LERR = .FALSE.
      STR = 'dihedral angle number'
      CALL GIMME(IUNIT,PRGSTR,STR)
      CALL CHPINT(NPHIH,PRGSTR,STR)
      IF (NPHIH .LT. 0 .OR. NPHIH .GT. MXPHIH) THEN
         CALL FLAGLN(PRGSTR)
         PRINT FMNII,'NPHIH',NPHIH
         PRINT FMBTII,0,'NPHIH','MXPHIH',MXPHIH
         STOP
      ENDIF
      CALL CHKMTY(PRGSTR,STR)

C iph,jph,kph,lph,icph
      STR = 'IPH,JPH,KPH,LPH,ICPH'
      DO I=1,NPHIH
         CALL GIMME(IUNIT,PRGSTR,STR)
         CALL CHPINT(IPH(I),PRGSTR,STR)
         CALL CHPINT(JPH(I),PRGSTR,STR)
         CALL CHPINT(KPH(I),PRGSTR,STR)
         CALL CHPINT(LPH(I),PRGSTR,STR)
         CALL CHPINT(ICPH(I),PRGSTR,STR)
         CALL CHKMTY(PRGSTR,STR)

C check iph,jph,kph,lph,icph
         IF (IPH(I) .LT. 1 .OR. IPH(I) .GT. NRP) THEN
            LERR = .TRUE.
            CALL FLAGLN(PRGSTR)
            PRINT 500,I
            PRINT 510,IPH(I)
            PRINT 520
         ENDIF

         IF (JPH(I) .LT. 1 .OR. JPH(I) .GT. NRP) THEN
            LERR = .TRUE.
            CALL FLAGLN(PRGSTR)
            PRINT 500,I
            PRINT 510,JPH(I)
            PRINT 520
         ENDIF

         IF (KPH(I) .LT. 1 .OR. KPH(I) .GT. NRP) THEN
            LERR = .TRUE.
            CALL FLAGLN(PRGSTR)
            PRINT 500,I
            PRINT 510,KPH(I)
            PRINT 520
         ENDIF

         IF (LPH(I) .LT. 1 .OR. LPH(I) .GT. NRP) THEN
            LERR = .TRUE.
            CALL FLAGLN(PRGSTR)
            PRINT 500,I
            PRINT 510,LPH(I)
            PRINT 520
         ENDIF

         IF (ICPH(I) .LT. 1 .OR. ICPH(I) .GT. NPTY) THEN
            LERR = .TRUE.
            CALL FLAGLN(PRGSTR)
            PRINT *,'illegal dihedral type in ',I
            PRINT *,'read IDA:',IPH(I)
            PRINT FMBTII,1,'IDA','NPTY',NPTY
            PRINT *
         ENDIF

      ENDDO

      CALL CHKEND(IUNIT,PRGSTR)

      IF (LERR) THEN
         STOP
      ENDIF

 500  FORMAT (1X,'atom sequence number out of range for',
     $     'dihedral angle',I5)
 510  FORMAT (1X,'read ASN:',I5)
 520  FORMAT (1X,'must be 0 < ASN <= NRP = ',I5,//)
C end GDI
      END

! ==================================================================
      SUBROUTINE GCHRG(IUNIT,PRGSTR)
! ==--------------------------------------------------------------==
!    Revision cmb Tokyo/Strasbourg june 2014
! ==--------------------------------------------------------------==
C includes
      INCLUDE 'formats.h'
      INCLUDE 'toposz.h'
      INCLUDE 'topoar.h'
C args
      INTEGER IUNIT
      CHARACTER PRGSTR*(*)
C local vars
      INTEGER I
      CHARACTER STR*(5)
      LOGICAL LERR
C begin
      LERR = .FALSE.
C ncag
      STR = 'NCAG'
      CALL GIMME(IUNIT,PRGSTR,STR)
      CALL CHPINT(NCAG,PRGSTR,STR)
      IF (NCAG .LT. 0 .OR. NCAG .GT. MAXCAG) THEN
         CALL FLAGLN(PRGSTR)
         PRINT FMNII,STR,NCAG
         PRINT FMBTII,0,STR,'MAXCAG',MAXCAG
         STOP
      ENDIF
      CALL CHKMTY(PRGSTR,STR)

C inc
      STR = 'INC'
      DO I=1,NCAG
         CALL GIMME(IUNIT,PRGSTR,STR)
         CALL CHPINT(INC(I),PRGSTR,STR)
         CALL CHKMTY(PRGSTR,STR)

C check inc
         IF (INC(I) .LT. 1 .OR. INC(I) .GT. NRP) THEN
            LERR = .TRUE.
            CALL FLAGLN(PRGSTR)
            PRINT *,'illegal atom sequence number in INC!'
            PRINT *,'INC(',I,') = ',INC(I)
            PRINT FMBTII,1,'INC(I)','NRP',NRP
            PRINT *
         ENDIF
      ENDDO

      CALL CHKEND(IUNIT,PRGSTR)

      IF (LERR) THEN
         PRINT *,PRGSTR,': Error in charge group definition!'
         PRINT *
         STOP
      ENDIF

C end GCHRG
      END

! ==================================================================
      SUBROUTINE GLJ(IUNIT,PRGSTR)
! ==--------------------------------------------------------------==
!    Revision cmb Tokyo/Strasbourg june 2014
! ==--------------------------------------------------------------==
      INCLUDE 'toposz.h'
      INCLUDE 'topoar.h'
      INCLUDE 'formats.h'
C args
      INTEGER IUNIT
      CHARACTER PRGSTR*(*)
C local vars
      INTEGER K,ITMP,JTMP,KDEST
      CHARACTER STR*(32)
C begin
C nratt2
      STR = 'NRATT2'
      CALL GIMME(IUNIT,PRGSTR,STR)
      CALL CHPINT(NRATT2,PRGSTR,STR)
      IF (NRATT2 .NE. NRATT*(NRATT+1)/2) THEN
         CALL FLAGLN(PRGSTR)
         PRINT *,'inconsistent value for NRATT2'
         PRINT *,'when reading C12, C6!'
         PRINT *,'read',NRATT2,',expected ',NRATT*(NRATT+1)/2
         STOP
      ENDIF
      CALL CHKMTY(PRGSTR,STR)

C     in this subroutine, we also calculate the MPAC table.
C     first, we set all values of MPAC to 0.
C     while reading, we set the MPAC entries, and check
C     for setting an entry twice.
C     We do not need to check that all entries have been
C     set after reading because if we have read NRATT2 different
C     entries, then we haven't missed any!
!$omp parallel do private(ITMP,JTMP) __COLLAPSE2
      DO ITMP=1,NRATT
         DO JTMP=1,NRATT
            MPAC(ITMP,JTMP) = 0
         ENDDO
      ENDDO

C c12,c6
      STR = 'ITYP,JTYP,C12,C6,CS12,CS6'
      DO K=1,NRATT2
         CALL GIMME(IUNIT,PRGSTR,STR)
         CALL CHPINT(ITMP,PRGSTR,STR)
         IF (ITMP .LT. 1 .OR. ITMP .GT. NRATT) THEN
            CALL FLAGLN(PRGSTR)
            PRINT FMNII,'ITYP',ITMP
            PRINT FMBTII,1,'ITYP','NRATT',NRATT
            STOP
         ENDIF
         CALL CHPINT(JTMP,PRGSTR,STR)
         IF (JTMP .LT. 1 .OR. JTMP .GT. NRATT) THEN
            CALL FLAGLN(PRGSTR)
            PRINT FMNII,'JTYP',JTMP
            PRINT FMBTII,1,'JTYP','NRATT',NRATT
            STOP
         ENDIF

         IF (MPAC(ITMP,JTMP) .NE. 0) THEN
            CALL FLAGLN(PRGSTR)
            PRINT *,'interaction type set twice!'
            PRINT *,'ITYP = ',ITMP
            PRINT *,'JTYP = ',JTMP
            STOP
         ENDIF

         IF (ITMP .LT. JTMP) THEN
            KDEST = ITMP+JTMP*(JTMP-1)/2
         ELSE
            KDEST = JTMP+ITMP*(ITMP-1)/2
         ENDIF
         MPAC(ITMP,JTMP) = KDEST
         MPAC(JTMP,ITMP) = KDEST

         CALL CHPREL(C12(KDEST),PRGSTR,STR)
         CALL CHPREL(C6(KDEST),PRGSTR,STR)
         CALL CHPREL(CS12(KDEST),PRGSTR,STR)
         CALL CHPREL(CS6(KDEST),PRGSTR,STR)
         CALL CHKMTY(PRGSTR,STR)

      ENDDO

      CALL CHKEND(IUNIT,PRGSTR)
C end GLJ
      END

! ==================================================================
      SUBROUTINE GSOLV(IUNIT,PRGSTR)
! ==--------------------------------------------------------------==
!    Revision cmb Tokyo/Strasbourg june 2014
! ==--------------------------------------------------------------==
C includes
      INCLUDE 'formats.h'
      INCLUDE 'toposz.h'
      INCLUDE 'topoar.h'
      INCLUDE 'fileio.h'
C args
      INTEGER IUNIT
      CHARACTER PRGSTR*(*)
C local vars
      INTEGER K,KDEST
      REAL TMPMSS
      CHARACTER STR*(24)
C begin
C nram
      STR = 'NRAM'
      CALL GIMME(IUNIT,PRGSTR,STR)
      CALL CHPINT(NRAM,PRGSTR,STR)
      IF (NRAM .LT. 0 .OR. NRAM .GT. MAXNRS) THEN
         CALL FLAGLN(PRGSTR)
         PRINT FMNII,STR,NRAM
         PRINT FMBTII,0,STR,'MAXNRS',MAXNRS
         STOP
      ENDIF
      CALL CHKMTY(PRGSTR,STR)
   
C winvs,anms,cgs,iacs
      STR = 'K,ANMS,IACS,MASS,CGS'
      DO K=1,NRAM
         CALL GIMME(IUNIT,PRGSTR,STR)
         CALL CHPINT(KDEST,PRGSTR,STR)
         IF (KDEST .LT. 1 .OR. KDEST .GT. NRAM) THEN
            CALL FLAGLN(PRGSTR)
            PRINT *,'inconsistent value for solvent atom number K= ',
     $           KDEST
            PRINT FMBTII,1,'K','NRAM',NRAM
            STOP
         ENDIF
         CALL CHPSTR(ANMS(KDEST),PRGSTR,STR)
         CALL CHPINT(IACS(KDEST),PRGSTR,STR)

         CALL CHPREL(TMPMSS,PRGSTR,STR)
         IF (TMPMSS .LE. 0.0) THEN
            CALL FLAGLN(PRGSTR)
            PRINT *,'mass of solvent atom zero or negative: ',TMPMSS
            STOP
         ENDIF
         WMASS(KDEST)  = TMPMSS
         WINVS(KDEST) = 1.0/TMPMSS

         CALL CHPREL(CGS(KDEST),PRGSTR,STR)
         CALL CHKMTY(PRGSTR,STR)
      ENDDO
      CALL CHKEND(IUNIT,PRGSTR)
C end GSOLV
      END

! ==================================================================
      SUBROUTINE GSLVCN(IUNIT,PRGSTR)
! ==--------------------------------------------------------------==
!    Revision cmb Tokyo/Strasbourg june 2014
! ==--------------------------------------------------------------==
      INCLUDE 'toposz.h'
      INCLUDE 'topoar.h'
      INCLUDE 'formats.h'
C args
      INTEGER IUNIT
      CHARACTER PRGSTR*(*)
C local vars
      INTEGER K
      CHARACTER STR*(24)
C begin
C ncons
      STR = 'NCONS'
      CALL GIMME(IUNIT,PRGSTR,STR)
      CALL CHPINT(NCONS,PRGSTR,STR)
      IF (NCONS .LT. 0 .OR. NCONS .GT. MXCONS)THEN
         CALL FLAGLN(PRGSTR)
         PRINT FMNII,STR,NCONS
         PRINT FMBTII,0,STR,'MXCONS',MXCONS
         STOP
      ENDIF
      CALL CHKMTY(PRGSTR,STR)

C icons,jcons,cons
      STR = 'ICONS,JCONS,CONS'
      DO K=1,NCONS
         CALL GIMME(IUNIT,PRGSTR,STR)
         CALL CHPINT(ICONS(K),PRGSTR,STR)
         CALL CHPINT(JCONS(K),PRGSTR,STR)
C check icons
         IF (ICONS(K) .LT. 1 .OR. ICONS(K) .GT. NRAM) THEN
            CALL FLAGLN(PRGSTR)
            STR = 'ICONS'
            PRINT FMNII,STR,ICONS(K)
            PRINT FMBTII,1,STR,'NRAM',NRAM
            STOP
         ENDIF
C check JCONS
         IF (JCONS(K) .LT. 1 .OR. JCONS(K) .GT. NRAM) THEN
            CALL FLAGLN(PRGSTR)
            STR = 'JCONS'
            PRINT FMNII,STR,ICONS(K)
            PRINT FMBTII,1,STR,'NRAM',NRAM
            STOP
         ENDIF

         CALL CHPREL(CONS(K),PRGSTR,STR)
         IF (CONS(K) .LE. 0.0) THEN
            PRINT FMNIR,'CONS',CONS(K)
            PRINT FMGTR,'CONS'
            STOP
         ENDIF
         CALL CHKMTY(PRGSTR,STR)
      ENDDO

      CALL CHKEND(IUNIT,PRGSTR)
C end GSLVCN
      END

! ==================================================================
      SUBROUTINE GPITBL(IUNIT,PRGSTR)
! ==--------------------------------------------------------------==
!    Revision cmb Tokyo/Strasbourg june 2014
! ==--------------------------------------------------------------==
      INCLUDE 'toposz.h'
      INCLUDE 'topoar.h'
      INCLUDE 'formats.h'
C args
      INTEGER IUNIT
      CHARACTER PRGSTR*(*)
C local vars
      INTEGER K
      CHARACTER STR*(24)
C begin
      STR = 'NPID'
      CALL GIMME(IUNIT,PRGSTR,STR)
      CALL CHPINT(NPID,PRGSTR,STR)
      IF (NPID .LT. 2 .OR. NPID .GT. MAXPID) THEN
         CALL FLAGLN(PRGSTR)
         PRINT FMNII,STR,NPID
         PRINT FMBTII,2,STR,'MAXPID',MAXPID
         STOP
      ENDIF
      CALL CHKMTY(PRGSTR,STR)

      STR = 'IPIC'
      DO K=1,NRP
         CALL GIMME(IUNIT,PRGSTR,STR)
         CALL CHPINT(IPIC(K),PRGSTR,STR)
         CALL CHKMTY(PRGSTR,STR)
      ENDDO

      STR = 'TPI'
      CALL GIMME(IUNIT,PRGSTR,STR)
      CALL CHPREL(TPI,PRGSTR,STR)
      CALL CHKMTY(PRGSTR,STR)

      STR = 'NPIT'
      CALL GIMME(IUNIT,PRGSTR,STR)
      CALL CHPINT(NPIT,PRGSTR,STR)
      CALL CHKMTY(PRGSTR,STR)
      IF (NPIT .GT. MAXPIT .OR. NPIT .LT. 0) THEN
         CALL FLAGLN(PRGSTR)
         PRINT FMNII,STR,NPIT
         PRINT FMBTII,0,STR,'MAXPIT',MAXPIT
         STOP
      ENDIF

      STR = 'CPI'
      DO K=1,NPIT
         CALL GIMME(IUNIT,PRGSTR,STR)
         CALL CHPREL(CPI(K),PRGSTR,STR)
         CALL CHKMTY(PRGSTR,STR)
      ENDDO

      STR = 'NPIB'
      CALL GIMME(IUNIT,PRGSTR,STR)
      CALL CHPINT(NPIB,PRGSTR,STR)
      CALL CHKMTY(PRGSTR,STR)
      IF (NPIB .GT. MAXPIB .OR. NPIB .LT. 0) THEN
         CALL FLAGLN(PRGSTR)
         PRINT FMNII,STR,NPIB
         PRINT FMBTII,0,STR,'MAXPIB',MAXPIB
         STOP
      ENDIF

      STR = 'IPIB, JPIB, ICPIB'
      DO K=1,NPIB
         CALL GIMME(IUNIT,PRGSTR,STR)
         CALL CHPINT(IPIB(K),PRGSTR,STR)
         CALL CHPINT(JPIB(K),PRGSTR,STR)
         CALL CHPINT(ICPIB(K),PRGSTR,STR)
         CALL CHKMTY(PRGSTR,STR)
      ENDDO

      STR = 'NPIA'
      CALL GIMME(IUNIT,PRGSTR,STR)
      CALL CHPINT(NPIA,PRGSTR,STR)
      IF (NPIA .GT. MAXPIA .OR. NPIA .LT. 0) THEN
         CALL FLAGLN(PRGSTR)
         PRINT FMNII,STR,NPIA
         PRINT FMBTII,0,STR,'MAXPIA',MAXPIA
         STOP
      ENDIF
      CALL CHKMTY(PRGSTR,STR)

      STR = 'WMCL'
      DO K=1,NPIA
         CALL GIMME(IUNIT,PRGSTR,STR)
         CALL CHPREL(WMCL(K),PRGSTR,STR)
         CALL CHKMTY(PRGSTR,STR)
      ENDDO

      CALL CHKEND(IUNIT,PRGSTR)
C end GPITBL
      END

! ==================================================================
      SUBROUTINE PRPTOP
! ==--------------------------------------------------------------==
C    Perform the conversion of the topology from data
C    on file to data as it is needed by the force routines.
!    Revision cmb Tokyo/Strasbourg june 2014
! ==--------------------------------------------------------------==
      INCLUDE 'toposz.h'
      INCLUDE 'topoar.h'
C local vars
      REAL CGFAC,PYE,DEG2RD,DGRDSQ
      INTEGER I
C begin
      PYE    = 4.E0*DATAN(1.0D+00)
      DEG2RD = PYE/180.0
      DGRDSQ = DEG2RD*DEG2RD

      CGFAC = SQRT(FPEPSI)

C convert the solute and solvent charges
!$omp parallel do private(I)
      DO I=1,NRP
         CG(I) = CG(I)*CGFAC
      ENDDO
!$omp parallel do private(I)
      DO I=1,NRAM
         CGS(I) = CGS(I)*CGFAC
      ENDDO

C bonds: no conversion needed

C convert bond angles types
!$omp parallel do private(I)
      DO I=1,NTTY
         T0(I) = COS(T0(I)*DEG2RD)
      ENDDO

C convert improper dihedrals
!$omp parallel do private(I)
      DO I=1,NQTY
         CQ(I) = CQ(I)/DGRDSQ
         Q0(I) = Q0(I)*DEG2RD
      ENDDO

C dihedrals: no conversion needed

C solvent: take square of distances
!$omp parallel do private(I)
      DO I=1,NCONS
         CONS(I) = CONS(I)*CONS(I)
      ENDDO
C end prptop
      END

