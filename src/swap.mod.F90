MODULE swap
  IMPLICIT NONE

  ! ==================================================================
  ! == VARIABLES FOR SWAP FILES                                     ==
  ! == USE FOR K POINTS (BLOCK OPTION)                              ==
  ! ==================================================================
  ! == ISWMAX:  MAXIMAL NUMBER OF DIFFERENT VARIABLES               ==
  ! == ISWFILE: MAXIMAL NUMBER OF SWAP FILES                        ==
  ! ==--------------------------------------------------------------==
  INTEGER, PARAMETER :: iswmax=20,iswfilemax=3 
  ! ==--------------------------------------------------------------==
  ! == IN ORDER TO AVOID RECURSIVE CALL (WITH STOPGM ROUTINE)       ==
  ! ==--------------------------------------------------------------==
  INTEGER, PARAMETER :: i_nostopgm=0,i_stopgm=1 
  ! ==--------------------------------------------------------------==
  ! == ISWFIRST: 0 IF NOT INITIALIZED                               ==
  ! == SMYPID:   STRING USED FOR THE NAME OF SWAP FILES             ==
  ! == NSWBLOCK: NUMBER OF BLOCK FOR EACH ARRAY                     ==
  ! == SWFILE(ISWFILEMAX):  NAME OF SWAP FILES                      ==
  ! == IUNITSW(ISWFILEMAX): UNIT NUMBER                             ==
  ! == ACCESSW(ISWFILEMAX): TYPE OF ACCESS                          ==
  ! ==    S=SEQUENTIAL, REWRITE FORBIDDEN                           ==
  ! ==    D=DIRECT    , REWRITE PERMITTED, BUT DATA SIZE IS FIXED   ==
  ! == ISIDISW(ISWFILEMAX): SIZE OF RECORDS FOR DIRECT ACCESS       ==
  ! == ISWRES(ISWFILEMAX):  RESET NUMBER                            ==
  ! ==--------------------------------------------------------------==
  CHARACTER (len=6):: smypid
  CHARACTER (len=255) :: swfile(iswfilemax)
  CHARACTER (len=1) :: accessw(iswfilemax)
  INTEGER, ALLOCATABLE :: iunitsw(:)
  INTEGER, ALLOCATABLE :: isidisw(:)
  INTEGER, ALLOCATABLE :: iswres(:)

  INTEGER :: iswfirst=0,nswblock

  ! ==--------------------------------------------------------------==
  ! == SWNAME(ISWMAX,ISWFILEMAX):                                   ==
  ! ==    THE POSITION OF NAME OF A VARIABLE GIVES INDEX            ==
  ! ==    IN ANOTHER ARRAYS                                         ==
  ! == ISWMEM(ISWMAX,ISWFILEMAX):                                   ==
  ! ==    GIVES FOR EACH VARIABLE THE INDEX IKPT WHICH IS IN MEMORY ==
  ! == ISWPOS(NBLKP,ISWMAX,ISWFILEMAX):                             ==
  ! ==    FOR AN INDEX GIVES THE RECORD NUMBER FOR EACH BLOCK       ==
  ! ==          OF NKPNT K-POINTS                                   ==
  ! == ISWSIZE(NBLKP*ISWMAX,ISWFILEMAX):                            ==
  ! ==    SIZE OF THE ITH RECORD                                    ==
  ! == ISWACCE(2,NBLKP*ISWMAX,ISWFILEMAX):                          ==
  ! ==    NUMBER OF WRITE(1) AND READ(2) ACCESS FOR ITH RECORD      ==
  ! == ISWCUR(ISWFILEMAX): CURRENT RECORD                           ==
  ! == ISWEND(ISWFILEMAX): TOTAL NUMBER OF VARIABLES                ==
  ! == ISWREC(ISWFILEMAX): TOTAL NUMBER OF RECORDS                  ==
  ! ==--------------------------------------------------------------==
  CHARACTER(len=8) :: swname(iswmax,iswfilemax)
  INTEGER, ALLOCATABLE :: iswmem(:,:)
  INTEGER, ALLOCATABLE :: iswpos(:,:,:)
  INTEGER, ALLOCATABLE :: iswsize(:,:)
  INTEGER, ALLOCATABLE :: iswacce(:,:,:)
  INTEGER, ALLOCATABLE :: iswcur(:)
  INTEGER, ALLOCATABLE :: iswend(:)
  INTEGER, ALLOCATABLE :: iswrec(:)

  ! ==================================================================
  ! == PARAMETER USED WITH TKBCALC AND TKBONE OPTION                ==
  ! ==--------------------------------------------------------------==
  ! == ISW_CALC Index IFILE for ISWACCE (store the number of        ==
  ! ==          calculations for each array)                        ==
  ! ==--------------------------------------------------------------==
  INTEGER, PARAMETER :: isw_calc=2 
  INTEGER, PARAMETER :: isw_hgkp   = 1,isw_hgkm   = 2,isw_maskgw = 3,&
       isw_twnl= 4,isw_emk2   = 5,isw_eigkr  = 6,isw_alm    = 7
  ! ==--------------------------------------------------------------==
  ! == TSWCALC No swap but calculation for arrays (=TKBCALC)        ==
  ! == TSWAPC0 Swap C0 wavefunctions (=.NOT.TKBLONE)                == 
  ! ==--------------------------------------------------------------==
  LOGICAL :: tswcalc,tswapc0
  ! ==================================================================
END MODULE swap
