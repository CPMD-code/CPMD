C $Id: md.h,v 1.1 2006-12-27 11:27:40 itavern Exp $ -*-fortran-*-
C constants and variables to control behaviour of PROMD

COMMVAR MDTITLE, MDTILE,MDTLIN
C length of titles in GROMOS
      INTEGER MDTILE
      PARAMETER (MDTILE = 80)
C number of lines a title block can have
      INTEGER MDTLIN
      PARAMETER (MDTLIN = 10)
COMMEND

COMMVAR MAXNRE,MXNRE2,NRESOL
C maximum length of NRE array for energy groups
      INTEGER MAXNRE
      PARAMETER (MAXNRE = 10)
C
C the length of the arrays needed for summing up
C the contributions
      INTEGER MXNRE2
      PARAMETER (MXNRE2 = MAXNRE*(MAXNRE+1)/2)
C the value used in force calculations
C for solvent contribution
C      INTEGER NRESOL
C      PARAMETER (NRESOL = MAXNRE)
COMMEND


C---------------------------------------------
C for LOGICAL switches (E.G. NTT4), use the following

      INTEGER ITRUE,IFALSE
      PARAMETER (ITRUE =1,IFALSE = 0)

COMMVAR NTX,NTXX,NTXV,NTSX
C vals for NTX (coordinate reading/writing routines)
C initial configuration file formats (tape 21)
C   NTX       X    V   SX
C     1     yes     0  
C     2     yes   yes
C     3     yes   yes  yes

      INTEGER NTXX
      PARAMETER (NTXX = 1)

      INTEGER NTXV
      PARAMETER (NTXV = 2)

      INTEGER NTSX
      PARAMETER (NTSX = 3)

      INTEGER NTXMIN,NTXMAX
      PARAMETER (NTXMIN = NTXX, NTXMAX = NTSX)
COMMEND

COMMVAR NTXO,NTXOUF,NTXOFO
C vals for NTXO
      INTEGER NTXOUF
      PARAMETER (NTXOUF = 0)
      INTEGER NTXOFO
      PARAMETER (NTXOFO = 1)

      INTEGER NTXOMI,NTXOMA
      PARAMETER (NTXOMI = 0, NTXOMA = 1)
COMMEND

COMMVAR NRDBOX,NRDBXN,NRDBXY
C vals for NRDBOX
C a logical for whether to read in the box from tape 5 or tape21
C the actual BOX array and NTB is defined in box.h
C if NRDBOX = NRDBXY then the box is read from tape 21.
      INTEGER NRDBXN
      PARAMETER  (NRDBXN = IFALSE)
      INTEGER NRDBXY
      PARAMETER (NRDBXY = ITRUE)
COMMEND

COMMVAR NTC,NTCNON,NTCDOH,NTCDOB
C vals for NTC shake control)
      INTEGER NTCNON
      PARAMETER (NTCNON = 1)
 
      INTEGER NTCDOH
      PARAMETER (NTCDOH = 2)
 
      INTEGER NTCDOB
      PARAMETER (NTCDOB = 3)
 
      INTEGER NTCMIN,NTCMAX
      PARAMETER (NTCMIN = NTCNON,NTCMAX = NTCDOB)
COMMEND


COMMVAR NTT,NTTOFF,NTTONE,NTTTWO,NTTHRE
C vals for NTT switches -- temperature coupling
      INTEGER NTTOFF, NTTONE, NTTTWO,NTTHRE
C no temperature coupling
      PARAMETER (NTTOFF = 0)
      PARAMETER (NTTONE = 1)
      PARAMETER (NTTTWO = 2)
      PARAMETER (NTTHRE = 3)

      INTEGER NTTMIN, NTTMAX
      PARAMETER (NTTMIN = -NTTHRE, NTTMAX = NTTHRE)
COMMEND


COMMVAR NTP,NTPOFF,NTPISO,NTPANI
C vals for NTP -- pressure coupling
      INTEGER NTPOFF, NTPISO, NTPANI
C no pressure coupling
      PARAMETER (NTPOFF = 0)
C isotropic pos scaling
      PARAMETER (NTPISO = 1)
C anisotropic pos scaling
      PARAMETER (NTPANI = 2)

      INTEGER NTPMIN,NTPMAX
      PARAMETER (NTPMIN = NTPOFF,NTPMAX = NTPANI)
COMMEND


COMMVAR INIT,INSHVX,INSHV,INNOSH,INITCO,INITCO,INITMI,INITMA
C values for INIT -- starting procedures in RUNMD.f
C  INIT
C     shake X
C            shake V
C COM removal if NTCM is 1
C     1   yes    yes   yes
C     2    no    yes   yes
C     3    no     no   yes
C     4    no     no    no

      INTEGER INSHVX
      PARAMETER (INSHVX = 1)

      INTEGER INSHV
      PARAMETER (INSHV = 2)

      INTEGER INNOSH
      PARAMETER (INNOSH=3)

      INTEGER INITCO
      PARAMETER (INITCO = 4)

      INTEGER INITMI, INITMA
      PARAMETER (INITMI = INSHVX, INITMA = INITCO)
COMMEND


COMMVAR NTID, NTIDOF,NTIDNH,NTIDAL,NTIDMI,NTIDMA
C vals for NTID
C     NTID = 0 : IMPROPER DIHEDRAL INTERACTION IS SKIPPED
C          = 1 : IMPROPER DIHEDRAL INTERACTIONS INVOLVING NO H-ATOMS
C                ARE CALCULATED
C          = 2 : IN ADDITION, IMPROPER DIHEDRAL INTERACTIONS INVOLVING
C                H-ATOMS ARE CALCULATED

      INTEGER NTIDOF
      PARAMETER (NTIDOF = 0)

      INTEGER NTIDNH
      PARAMETER (NTIDNH = 1)

      INTEGER NTIDAL
      PARAMETER (NTIDAL = 2)

      INTEGER NTIDMI,NTIDMA
      PARAMETER (NTIDMI = NTIDOF, NTIDMA = NTIDAL)
COMMEND


COMMVAR NTTP, NTPPNO,NTPPYE
C values for NTPP (dihedral angle monitoring, see SUBR. L<FORCE>)
      INTEGER NTPPNO,NTPPYE
      PARAMETER (NTPPNO = IFALSE, NTPPYE = ITRUE)
COMMEND


COMMVAR NTFBNH,NTFBND,NTFANH,NTFANG,NTFIDH,NTFIDE,NTFDHH,NTFIDE
C indices into the L<LTF> array
      INTEGER NTFBNH,NTFBND
      PARAMETER (NTFBNH = 1,NTFBND = 2)

      INTEGER NTFANH, NTFANG
      PARAMETER (NTFANH = 3,NTFANG = 4)

C improper dihedrals
      INTEGER NTFIDH,NTFIDE
      PARAMETER (NTFIDH = 5,NTFIDE = 6)
COMMEND
COMMVAR NTFDHH,NTFDIH,NTFCG, NTFNBN,MAXNTF
C dihedrals
      INTEGER NTFDHH, NTFDIH
      PARAMETER ( NTFDHH = 7, NTFDIH = 8)

      INTEGER NTFCG, NTFNBN
      PARAMETER (NTFCG = 9, NTFNBN = 10)

      INTEGER MAXNTF
      PARAMETER (MAXNTF = NTFNBN)
COMMEND

COMMVAR NTG,NTGOFF, NTGLAM,NTGMU,NTGBOT
C values for NTG
C     0: no perturbation
C     1: calculate dV/dRLAM perturbation
C     2: calculate dV/dRMU perturbation
C     3: calculate both derivatives
C
      INTEGER NTGOFF, NTGLAM,NTGMU,NTGBOT
      PARAMETER (NTGOFF = 0,NTGLAM = 1,NTGMU = 2, NTGBOT = 3)
      INTEGER NTGMIN, NTGMAX
      PARAMETER (NTGMIN = NTGOFF, NTGMAX = NTGBOT)
COMMEND

COMMVAR NRDGL,NRDGLN,NRDGLY
C values for NRDGL
      INTEGER NRDGLN,NRDGLY
      PARAMETER (NRDGLN = IFALSE, NRDGLY = ITRUE)
COMMEND


COMMVAR NTR, NTROFF,NTRCHO,NTRBFA,NTRCON
C     values for NTR (position restraining)
C     NTR
C     0: no position re(con)straining
C     1: position restraining using CHO
C     2: position restraining using CHO/ ATOMIC B-FACTORS
C     3: position constraining

      INTEGER NTROFF,NTRCHO,NTRBFA,NTRCON
      PARAMETER (NTROFF=0, NTRCHO=1,NTRBFA=2,NTRCON=3)
      INTEGER NTRMIN, NTRMAX
      PARAMETER (NTRMIN=NTROFF, NTRMAX = NTRCON)
COMMEND

COMMVAR NRDRX, NRDRXN,NRDRXY
C values for NRDRX
      INTEGER NRDRXN
      PARAMETER (NRDRXN = IFALSE)

      INTEGER NRDRXY
      PARAMETER (NRDRXY = ITRUE)
COMMEND


COMMVAR NTDR,NTDROF,NTDRCD,NTDRW0,NTDRMI, NTDRMA
C values for NTDR (distance restraining)
C     NTDR
C     0 : no distance restraining
C     -1,1 : use CDIS
C     -2,2:  use W0*CDIS
C if NTDR < 0 : use time averageing
C if NTDR > 0 : no time averaging

      INTEGER NTDROF
      PARAMETER (NTDROF = 0)

C remember to check for the absolute values
      INTEGER NTDRCD
      PARAMETER (NTDRCD = 1)

      INTEGER NTDRW0
      PARAMETER (NTDRW0 = 2)

      INTEGER NTDRMI, NTDRMA
      PARAMETER (NTDRMI = -NTDRW0, NTDRMA = NTDRW0)
COMMEND


COMMVAR NTDLR,NTDLRN,NTDLRC,NTDLRF,NTDLR1, NTDLR2
C values for NTDLR (dihedral restraining)
C NTDLR
C     0: off
C     1: perform dihedral restraining with CDLR
C     2: perform dihedral restraining with CDLR* CPLR
C
      INTEGER NTDLRN
      PARAMETER (NTDLRN = 0)

      INTEGER NTDLRC
      PARAMETER (NTDLRC = 1)

      INTEGER NTDLRF
      PARAMETER (NTDLRF = 2)

      INTEGER NTDLR1, NTDLR2
      PARAMETER (NTDLR1 = NTDLRN, NTDLR2 = NTDLRF)
COMMEND


COMMVAR NTPW, NTPWBI,NTPWFO,NTPWMI,NTPWMA
C values for NTPW
      INTEGER NTPWBI
      PARAMETER (NTPWBI = IFALSE)

      INTEGER NTPWFO
      PARAMETER (NTPWFO = ITRUE)

      INTEGER NTPWMI,NTPWMA
      PARAMETER (NTPWMI = NTPWBI, NTPWMA = NTPWFO)
COMMEND


COMMVAR NT4DIM,NT4OFF,NT4CWD,NT4RED,NT4MIN,NT4MAX
C values for NT4DIM
C	0	perform 3D simulation
C	>0	perform 4D simulation
C	1	CW4DA
C	2	CW4DA*C4D
C---
      INTEGER NT4OFF
      PARAMETER (NT4OFF = 0)

      INTEGER NT4CWD
      PARAMETER (NT4CWD = 1)

      INTEGER NT4RED
      PARAMETER (NT4RED = 2)

C     add other values here.
C     remember to update NT4MAX

      INTEGER NT4MIN,NT4MAX
      PARAMETER (NT4MIN = NT4OFF, NT4MAX = NT4RED)
COMMEND


COMMVAR NT4XI,
C values for NT4XI, NT4X0V,NT4XBM,NT4XRM,NT4XRR,NT4XMI,NT4XMA
COMMVERB
C-----------------
C  NT4XI	meaning
C		init coords     init velocities
C	1	set 0	        Maxwell dist.
C	2	Boltzm. dist	Maxwell dist.
C	3	read from 21    Maxwell dist.
C	4	read from 21    read from 21  used for continuation run
C-----------------
COMMEND
      INTEGER NT4X0V
      PARAMETER (NT4X0V = 1)

      INTEGER NT4XBM
      PARAMETER (NT4XBM = 2)

      INTEGER NT4XRM
      PARAMETER (NT4XRM = 3)

      INTEGER NT4XRR
      PARAMETER (NT4XRR = 4)

      INTEGER NT4XMI,NT4XMA
      PARAMETER (NT4XMI = NT4X0V,NT4XMA = NT4XRR)
COMMEND



COMMVAR NTT4,NTT4OF,NTT4ON
C values for NTT4, temperature coupling of 4th dimension
      INTEGER NTT4OF
      INTEGER NTT4ON
      PARAMETER (NTT4OF = IFALSE,NTT4ON = ITRUE)
COMMEND


COMMVAR N4DBON,N4DBAN,N4DDIH,N4DNBD,N4DRST,N4DHAR,MX4DIN
C     ----4D interaction types that can be switched on/off
C         These are indices into the L<NDO4D> array.

C     bonds
      INTEGER N4DBON
      PARAMETER (N4DBON = 1)

C     bond angle
      INTEGER N4DBAN
      PARAMETER (N4DBAN = 2)

C     dihedral angle
      INTEGER N4DDIH
      PARAMETER (N4DDIH = 3)

C     non-bonded
      INTEGER N4DNBD
      PARAMETER (N4DNBD = 4)

C     distance restraining
      INTEGER N4DRST
      PARAMETER (N4DRST = 5)

C     4D harmonic osc.
      INTEGER N4DHAR
      PARAMETER (N4DHAR = 6)

C     add other interaction types here.
C     remember to adjust MX4DIN below!

C     number of 4D interaction types that can be switched on/off
      INTEGER MX4DIN
      PARAMETER (MX4DIN = N4DHAR)
COMMEND


COMMVAR NTCW4D,NTC4NO,NTC4YE
C values for NTCW4D
      INTEGER NTC4NO
      PARAMETER (NTC4NO = IFALSE)

      INTEGER NTC4YE
      PARAMETER (NTC4YE = ITRUE)
COMMEND

COMMVAR NTEM, NTEMSD,NTEMCG,NTEMMI,NTEMMA
C values for NTEM
      INTEGER NTEMSD,NTEMCG
C use steepest descent
      PARAMETER (NTEMSD = 1)
C use conjugate gradients
      PARAMETER (NTEMCG = 2)

      INTEGER NTEMMI,NTEMMA
      PARAMETER (NTEMMI = NTEMSD, NTEMMA = NTEMCG)
COMMEND



COMMVAR NTFR, NTFR0,NTFRFR,NTFRGA,NTFRCA,NTFRMI,NTFRMA
C values for NTFR
COMMVERB
C     NTFR     atomic friction coefficients
C     NTFR0    set to zero
C     NTFRFR   set equal to CFRIC
C     NTFRGA   set equal to CFRIC*GAM (GAM is read from file)
C     NTFRCA   are calculated in subr. FRIC using
C              neighbour atom information
COMMEND
      INTEGER NTFR0
      INTEGER NTFRFR
      INTEGER NTFRGA
      INTEGER NTFRCA
      PARAMETER (NTFR0 = 0)
      PARAMETER (NTFRFR= 1)
      PARAMETER (NTFRGA= 2)
      PARAMETER (NTFRCA= 3)

      INTEGER NTFRMI,NTFRMA
      PARAMETER (NTFRMI = NTFR0,NTFRMA = NTFRCA)
COMMEND


COMMVAR NT4XO,NT4XON,NT4XOY
C values for NT4XO
      INTEGER NT4XON
      INTEGER NT4XOY
      PARAMETER (NT4XON = IFALSE, NT4XOY = ITRUE)
COMMEND


COMMVAR NTPI,NTPION, NTPIOF
C values for NTPI
C NTPION -> path integral on
C NTPIOF -> path integral off

      INTEGER NTPION, NTPIOF
      PARAMETER (NTPION = ITRUE, NTPIOF = IFALSE)
COMMEND


COMMVAR NTJR,NTJROF,NTJCJR,NTJRED,NTJRMI,NTJRMA
C values for NTJR
C     0: off
C     -1,1: use CJR*(acos^^2 + bcos + c)
C     -2,2: use CJR*CPJR*(acos^^2 + bcos + c)
C     NTJR < 0: time averaging
C     NTJR > 0: no time averaging

      INTEGER NTJROF,NTJCJR,NTJRED
      PARAMETER (NTJROF = 0)
      PARAMETER (NTJCJR = 1)
      PARAMETER (NTJRED = 2)

      INTEGER NTJRMI,NTJRMA
      PARAMETER (NTJRMI = -NTJRED, NTJRMA = NTJRED)
COMMEND

COMMVAR NTJRH,NTJRHH, NTJRHF
C values for NTJRH
C    0: use half harmonic potential
C    1: use harmonic potential
      INTEGER NTJRHH, NTJRHF
      PARAMETER (NTJRHH = 0, NTJRHF = 1)
COMMEND

COMMVAR NRDJR,NRDJRY,NRDJRN
C values for NRDJR
C NRDJR
C     0 : don't read in initial averages
C     1 : read in initial JCOUPLINGAV block from unit 21.
      INTEGER NRDJRY,NRDJRN
      PARAMETER (NRDJRY = ITRUE, NRDJRN = IFALSE)
COMMEND


COMMVAR NTLE, NTLEOF,NTLEGS,NTLEQU,NTLEMA,NTLEMI
C values for NTLE
C     0: off
C     1: local elevation with gaussian potential
C     2: local elevation with inverse quadratic
      INTEGER NTLEOF,NTLEGS,NTLEQU
      PARAMETER (NTLEOF = 0)
      PARAMETER (NTLEGS = 1)
      PARAMETER (NTLEQU = 2)

      INTEGER NTLEMA,NTLEMI
      PARAMETER (NTLEMI = NTLEOF,NTLEMA = NTLEQU)
COMMEND

COMMVAR NRDLE,NRDLEY,NRDLEN
C values NRDLE
C     0: don't read LE memory
C     1: read LE memory from unit 21 for continuation run 
      INTEGER NRDLEY,NRDLEN
      PARAMETER (NRDLEY = ITRUE, NRDLEN = IFALSE)
COMMEND


C----------------------------------------------------
C common blocks, for each block in the control file
C separate common blocks for ints/logicals, reals
C and chars

C start block
      INTEGER NPM,NSM,NTX,IG,NTXO
      COMMON /MSTART/ NPM,NSM,NTX,IG,NTXO
      REAL TEMPI,HEAT,BOLTZ
      COMMON /MSRTFP/TEMPI,HEAT,BOLTZ

      INTEGER NRDBOX
      COMMON /MDBOXY/NRDBOX

C tcouple and pcouple block

C we have three ntt switches for 3D
C the 4D switch is separate
      INTEGER NTTPIR,NTTPCM,NTTSLV,NTTNUM
      PARAMETER (NTTPIR = 1,NTTPCM = 2, NTTSLV = 3)
      PARAMETER (NTTNUM = NTTSLV)
      INTEGER NTT,NTP
      COMMON /MBEINT/NTT(NTTNUM),NTP

      REAL TEMP0,TAUT,PRES0,COMP,TAUP
      COMMON /MBEREF/
     $     TEMP0(NTTNUM),TAUT(NTTNUM),PRES0,COMP,TAUP

C nsp block
C      INTEGER NSP
C      COMMON /MNSP/ NSP(MAXNSP)

C com block
      INTEGER NDFMIN,NSCM
      LOGICAL LTCM
      COMMON /MCOM/ NDFMIN,NSCM,LTCM

C step block
      INTEGER NSTLIM,INIT
      COMMON /MSTEP/ NSTLIM,INIT
      REAL T,DT
      COMMON /MSTPFP/ T,DT

C shake block
      INTEGER NTC
      COMMON /MSHAKE/ NTC
      REAL SHKTOL
      COMMON /MSHKFP/ SHKTOL

C force block
      LOGICAL LTF
      INTEGER NRE,NEGR
      COMMON /MFORCE/ LTF(MAXNTF),NEGR,NRE(MAXNRE)

C plist block
      INTEGER NSNB
      LOGICAL LTNB
      COMMON /MPLIST/ NSNB,LTNB
      REAL RCUTP,RSWI2,RCUT2,CELM,RCUTL
      COMMON /MPLSTP/ RCUTP,RSWI2,RCUT2,CELM,RCUTL

C print block
      INTEGER NTPR,NTPP,NTPL
      COMMON /MPRINT/ NTPR,NTPP,NTPL

C write block
      INTEGER NTWX,NTWSE,NTWV,NTWE,NTWG,NTPW
      COMMON /MWRITE/ NTWX,NTWSE,NTWV,NTWE,NTWG,NTPW

C posrest block
      INTEGER NTR,NRDRX
      COMMON /MPOSRE/ NTR,NRDRX
      REAL CHO
      COMMON /MPOSFP/ CHO

C distrest block
      INTEGER NTDR
      LOGICAL LRDDR
      COMMON /MDISRE/NTDR,LRDDR
      REAL CDIS,DR0,TAUDR
      COMMON /MDREFP/ CDIS,DR0,TAUDR

C perturb block
      INTEGER NTG,NRDGL,NLAM,MMU
      COMMON /MPERTU/ NTG,NRDGL,NLAM,MMU
      REAL RLAM,DLAMT,RMU,DMUT
      REAL ALPHLJ,ALPHC
      COMMON /MPRTFP/ RLAM,DLAMT,RMU,DMUT,ALPHLJ,ALPHC

C diherest block
      INTEGER NTDLR
      COMMON /MDIHER/NTDLR
      REAL CDLR
      COMMON /MDIFP/ CDLR

C longrange block
      REAL EPSRF,APPAK,RCRF
      COMMON /MLONGR/EPSRF,APPAK,RCRF

C fourdim block
      INTEGER NT4DIM,NDFMI4,NTT4
      INTEGER NDO4D,NTCW4D,NT4XI,NT4XO

      COMMON /FOURD/ NT4DIM,NDFMI4,NTT4,
     $     NDO4D(MX4DIN),NTCW4D,NT4XI,NT4XO

      REAL CW4DA,TEMP4I,TEMP04,TAUT4,CW4DB,TEMP0B
      COMMON /FOURFP/CW4DA,TEMP4I,TEMP04,TAUT4,CW4DB,TEMP0B


C block title
      CHARACTER MDTITL*(MDTILE)
      COMMON /MTITLE/MDTITL(MDTLIN)
      INTEGER MDTLNS
      COMMON  /MTTITL/MDTLNS

COMMVAR IDOPRO, IDOMD,IDOEM,IDOSD
C values for IDOPRO
C IDOPRO is NOT read in by RDMD, but is set according to
C the existence of certain blocks in the input file.
C If a   MINIMIZE block exists, IDOPRO = IDOEM
C If a STOCHASTIC block exists, IDOPRO = IDOSD
C If neither exists,            IDOPRO = IDOMD
C
C A MINIMIZE and a STOCHASTIC block are mutually exclusive.

C the operations we can do
      INTEGER IDOMD,IDOEM,IDOSD
      PARAMETER (IDOMD = 0,IDOEM = 1,IDOSD = 3)

      INTEGER IDOPRO
C auxiliary for reading in by rdmd and subsequent checking
      COMMON /MDDDD/IDOPRO
COMMEND


C energy min block
      INTEGER NTEM,NCYC
      REAL DELE,DX0,DXM
      COMMON /EMINT/NTEM,NCYC
      COMMON /EMREL/DELE,DX0,DXM

C sd block
      INTEGER NTFR,NSFR,NBREF
      REAL RCUTF,CFRIC,TEMPSD
      COMMON /SDINT/ NTFR,NSFR,NBREF
      COMMON /SDREL/ RCUTF,CFRIC,TEMPSD


C path integral block
      INTEGER NTPI
      COMMON /PIMD_COMMON/NTPI

C J val restraining block
      INTEGER NTJR,NRDJR,NTJRH
      COMMON /JVINT/NTJR,NRDJR,NTJRH

      REAL TAUJR,CJR
      COMMON /JBREL/TAUJR,CJR

C local elevation block
      INTEGER NTLE,NRDLE
      COMMON /LEINT/NTLE,NRDLE

      REAL CWLE
      COMMON /LEREL/CWLE

C-----------end of common blocks --------------------
C here some useful constants
C a very small number
      REAL EPS
      PARAMETER (EPS = 1.E-6)

C a very large number
      REAL RELBIG
      PARAMETER (RELBIG = 10E24)
