#include "cpmd_global.h"

MODULE rwswap_utils
  USE utils,                           ONLY: fskip
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

!!$public :: wkpt_swap
!!$public :: rkpt_swap
!!$public :: kpt_noswap
!!$public :: ini_swap
!!$public :: direct_swap
!!$!public :: wri_swap
!!$!public :: rea_swap
!!$public :: go_swap
!!$public :: end_swap
!!$public :: del_swap
!!$public :: inq_swap

!!$contains

END MODULE rwswap_utils

! ==================================================================
! == ROUTINES TO MANAGE SWAP FILES                                ==
! == USED FOR K POINTS (BLOCK OPTIONS)                            ==
! ==--------------------------------------------------------------==
! == IF you need to add an array for k points, you have to change ==
! == three routines:                                              ==
! ==  - WKPT_SWAP                                                 ==
! ==  - RKPT_SWAP                                                 ==
! ==  - KPT_NOSWAP                                                ==
! == and the include file swap.inc.                               ==
! ==--------------------------------------------------------------==
SUBROUTINE wkpt_swap(c0,nstate,ikpt,tag)
  ! ==--------------------------------------------------------------==
  ! == PERFORMS SWAP FOR K-POINTS CALCULATION (WRITE)               ==
  ! == IFILE=1  C0  DIRECT ACCESS                                   ==
  ! == IFILE=2  HGKP HGKM TWNL EMK2 MASKGW SEQUENTIAL ACCESS        ==
  ! ==          INITIALIZED ONCE                                    ==
  ! == IFILE=3  ALM                        SEQUENTIAL ACCESS        ==
  ! ==          REINITIALIZED FOR EACH MOVE OF IONS COORDINATES     ==
  ! ==--------------------------------------------------------------==
  ! == TSWCALC  Calculate array depending on k points               ==
  ! ==          Use IFILE=1 for C0                                  ==
  ! ==          ISWACCE(:,:,ISW_CALC) is used to count the number   ==
  ! ==                 of calculation                               ==
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE system , ONLY:maxsys,ncpw,nkpbl,nkpt
  USE parac, ONLY : paral,parai
  USE elct , ONLY:crge
  USE cppt , ONLY:twnl
  USE nlps , ONLY:imagp,nlm
  USE fint , ONLY:alm,emk2,fint1
  USE swap , ONLY:iswfirst,tswapc0,tswcalc
  USE kpts , ONLY:tkpts
  USE kpnt , ONLY:hgkm,hgkp
  USE sphe , ONLY:maskgw,maskl,tsphere
  USE phfac_utils, ONLY : calc_eigkr
  IMPLICIT NONE
  INTEGER                                    :: nstate
  COMPLEX(real_8) :: c0(nkpt%ngwk,nstate,nkpt%nkpnt)
  INTEGER                                    :: ikpt
  CHARACTER(len=*)                           :: tag

  CHARACTER(len=80)                          :: line
  INTEGER                                    :: i1, nkpoint

  i1=LEN(tag)
  line=tag(1:i1)
  nkpoint=nkpbl(ikpt)
  ! ==--------------------------------------------------------------==
  IF (iswfirst.EQ.0) THEN
     CALL ini_swap(0)
     IF (tswapc0) THEN
        CALL direct_swap(1,2*nkpt%ngwk*crge%n*nkpt%nkpnt)
     ENDIF
  ENDIF
  ! ==--------------------------------------------------------------==
  IF (tkpts%tkall) THEN
     IF (INDEX(line,'HGKP').NE.0)&
          CALL wri_swap(2,hgkp,ncpw%nhg*nkpoint,ikpt,'HGKP')
     IF (INDEX(line,'HGKM').NE.0)&
          CALL wri_swap(2,hgkm,ncpw%nhg*nkpoint,ikpt,'HGKM')
     IF (tsphere.AND.INDEX(line,'MASKGW').NE.0)&
          CALL wri_swap(2,maskgw,2*maskl*nkpoint,ikpt,'MASKGW')
     IF (INDEX(line,'TWNL').NE.0)&
          CALL wri_swap(2,twnl,nkpt%ngwk*maxsys%nhxs*maxsys%nsx*nkpoint,ikpt,'TWNL')
     IF (fint1%ttrot.AND.INDEX(line,'EMK2').NE.0)&
          CALL wri_swap(2,emk2,nkpt%ngwk*nkpoint,ikpt,'EMK2')
     IF (INDEX(line,'EIGKR').NE.0)&
          CALL kpt_noswap(ikpt,1,'EIGKR')
     IF (fint1%ttrot.AND.INDEX(line,'ALM').NE.0)&
          CALL wri_swap(3,alm,imagp*nlm*nlm*nkpoint,ikpt,'ALM')
  ELSEIF (tswcalc) THEN
     ! Calculate each time. Writing is doing nothing.
  ELSE
     CALL kpt_noswap(ikpt,1,line)
  ENDIF
  ! ==--------------------------------------------------------------==
  IF (INDEX(line,'C0').NE.0) THEN
     IF (tswapc0) THEN
        CALL wri_swap(1,c0,2*nkpt%ngwk*nstate*nkpoint,ikpt,'C0')
     ENDIF
  ENDIF
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE wkpt_swap
! ==================================================================
SUBROUTINE rkpt_swap(c0,nstate,ikpt,tag)
  ! ==--------------------------------------------------------------==
  ! == PERFORMS SWAP FOR K-POINTS CALCULATION (READ)                ==
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE rinforce_utils, ONLY : calc_twnl
  USE rkpnt_utils, ONLY : rkpnt, bzmesh, calc_maskgw, calc_hgk
  USE ehpsi_utils, ONLY : calc_emk2
  USE calc_alm_utils, ONLY : calc_k_alm
  USE system , ONLY:maxsys,ncpw,nkpbl,nkpt
  USE parac, ONLY : paral,parai
  USE cppt , ONLY:twnl
  USE nlps , ONLY:imagp,nlm
  USE fint , ONLY:alm,emk2,fint1
  USE kpts , ONLY:tkpts
  USE kpnt , ONLY:hgkm,hgkp
  USE sphe , ONLY:maskgw,maskl,tsphere
  USE swap , ONLY:iswacce,iswmem,isw_alm,isw_calc,isw_eigkr,isw_emk2,isw_hgkm,isw_hgkp,isw_maskgw,isw_twnl,tswapc0,tswcalc
  USE phfac_utils, ONLY : calc_eigkr
  IMPLICIT NONE
  INTEGER                                    :: nstate
  COMPLEX(real_8) :: c0(nkpt%ngwk,nstate,nkpt%nkpnt)
  INTEGER                                    :: ikpt
  CHARACTER(len=*)                           :: tag

  CHARACTER(len=80)                          :: line
  INTEGER                                    :: i1, nkpoint

  i1=LEN(tag)
  line=tag(1:i1)
  nkpoint=nkpbl(ikpt)
  ! ==--------------------------------------------------------------==
  IF (tkpts%tkall) THEN
     IF (INDEX(line,'HGKP').NE.0)&
          CALL rea_swap(2,hgkp,ncpw%nhg*nkpoint,ikpt,'HGKP')
     IF (INDEX(line,'HGKM').NE.0)&
          CALL rea_swap(2,hgkm,ncpw%nhg*nkpoint,ikpt,'HGKM')
     IF (tsphere.AND.INDEX(line,'MASKGW').NE.0)&
          CALL rea_swap(2,maskgw,2*maskl*nkpoint,ikpt,'MASKGW')
     IF (INDEX(line,'TWNL').NE.0)&
          CALL rea_swap(2,twnl,nkpt%ngwk*maxsys%nhxs*maxsys%nsx*nkpoint,ikpt,'TWNL')
     IF (fint1%ttrot.AND.INDEX(line,'EMK2').NE.0)&
          CALL rea_swap(2,emk2,nkpt%ngwk*nkpoint,ikpt,'EMK2')
     IF (INDEX(line,'EIGKR').NE.0)&
          CALL calc_eigkr(ikpt)
     IF (fint1%ttrot.AND.INDEX(line,'ALM').NE.0)&
          CALL rea_swap(3,alm,imagp*nlm*nlm*nkpoint,ikpt,'ALM')
     ! ==------------------------------------------------------------==
  ELSEIF (tswcalc) THEN
1000 CONTINUE
  ! Calculate each time we need.
     IF (INDEX(line,'HGKP').NE.0.OR.INDEX(line,'HGKM').NE.0) THEN
        IF (iswmem(isw_hgkp,isw_calc).NE.ikpt) THEN
           iswacce(2,isw_hgkp,isw_calc)=iswacce(2,isw_hgkp,isw_calc)+1
           iswacce(2,isw_hgkm,isw_calc)=iswacce(2,isw_hgkm,isw_calc)+1
           CALL calc_hgk(ikpt)
           iswmem(isw_hgkp,isw_calc)=ikpt
           iswmem(isw_hgkm,isw_calc)=ikpt
        ENDIF
     ENDIF
     ! Need HGKM and HGKP before.
     IF (tsphere.AND.INDEX(line,'MASKGW').NE.0) THEN
        IF (iswmem(isw_maskgw,isw_calc).NE.ikpt) THEN
           IF (iswmem(isw_hgkp,isw_calc).NE.ikpt) THEN
              ! Add HGKP calculation.
              IF (paral%io_parent)&
                   WRITE(line,'(A,A)') line(1:i1),' HGKP'
              i1=i1+5
              GOTO 1000
           ENDIF
           iswacce(2,isw_maskgw,isw_calc)&
                =iswacce(2,isw_maskgw,isw_calc)+1
           CALL calc_maskgw(ikpt)
           iswmem(isw_maskgw,isw_calc)=ikpt
        ENDIF
     ENDIF
     ! Need HGKM and HGKP and MASKGW before.
     IF (INDEX(line,'TWNL').NE.0) THEN
        IF (iswmem(isw_twnl,isw_calc).NE.ikpt) THEN
           IF (iswmem(isw_hgkp,isw_calc).NE.ikpt) THEN
              ! Add HGKP calculation.
              IF (paral%io_parent)&
                   WRITE(line,'(A,A)') line(1:i1),' HGKP'
              i1=i1+5
              GOTO 1000
           ENDIF
           IF (iswmem(isw_maskgw,isw_calc).NE.ikpt) THEN
              ! Add MASKGW calculation.
              IF (paral%io_parent)&
                   WRITE(line,'(A,A)') line(1:i1),' MASKGW'
              i1=i1+8
              GOTO 1000
           ENDIF
           iswacce(2,isw_twnl,isw_calc)=iswacce(2,isw_twnl,isw_calc)+1
           CALL calc_twnl(ikpt)
           iswmem(isw_twnl,isw_calc)=ikpt
        ENDIF
     ENDIF
     ! Need HGKM and HGKP before.
     IF (fint1%ttrot.AND.INDEX(line,'EMK2').NE.0) THEN
        IF (iswmem(isw_emk2,isw_calc).NE.ikpt) THEN
           IF (iswmem(isw_hgkp,isw_calc).NE.ikpt) THEN
              ! Add HGKP calculation.
              IF (paral%io_parent)&
                   WRITE(line,'(A,A)') line(1:i1),' HGKP'
              i1=i1+5
              GOTO 1000
           ENDIF
           iswacce(2,isw_emk2,isw_calc)=iswacce(2,isw_emk2,isw_calc)+1
           CALL calc_emk2(ikpt)
           iswmem(isw_emk2,isw_calc)=ikpt
        ENDIF
     ENDIF
     IF (INDEX(line,'EIGKR').NE.0) THEN
        IF (iswmem(isw_eigkr,isw_calc).NE.ikpt) THEN
           iswacce(2,isw_eigkr,isw_calc)&
                =iswacce(2,isw_eigkr,isw_calc)+1
           CALL calc_eigkr(ikpt)
           iswmem(isw_eigkr,isw_calc)=ikpt
        ENDIF
     ENDIF
     ! Need TWNL and EIGKR before.
     IF (fint1%ttrot.AND.INDEX(line,'ALM').NE.0) THEN
        IF (iswmem(isw_alm,isw_calc).NE.ikpt) THEN
           IF (iswmem(isw_hgkp,isw_calc).NE.ikpt) THEN
              ! Add HGKP calculation.
              IF (paral%io_parent)&
                   WRITE(line,'(A,A)') line(1:i1),' HGKP'
              i1=i1+5
              GOTO 1000
           ENDIF
           IF (iswmem(isw_eigkr,isw_calc).NE.ikpt) THEN
              ! Add EIGKR calculation.
              IF (paral%io_parent)&
                   WRITE(line,'(A,A)') line(1:i1),' EIGKR'
              i1=i1+6
              GOTO 1000
           ENDIF
           iswacce(2,isw_alm,isw_calc)=iswacce(2,isw_alm,isw_calc)+1
           CALL calc_k_alm(ikpt)
           iswmem(isw_alm,isw_calc)=ikpt
        ENDIF
     ENDIF
     ! ==------------------------------------------------------------==
  ELSE
     CALL kpt_noswap(ikpt,0,line)
  ENDIF
  ! ==--------------------------------------------------------------==
  IF (INDEX(line,'C0').NE.0) THEN
     IF (tswapc0) THEN
        CALL rea_swap(1,c0,2*nkpt%ngwk*nstate*nkpoint,ikpt,'C0')
     ENDIF
  ENDIF
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE rkpt_swap
! ==================================================================
SUBROUTINE kpt_noswap(ikpt,iflag,tag)
  ! ==--------------------------------------------------------------==
  ! == DO NOT USE SWAP FILE, CHANGE ADDRESSES OF POINTER            ==
  ! == IFLAG=1 WRITE                                                ==
  ! ==       0 READ                                                 ==
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE system , ONLY:nkpt
  USE parac, ONLY : paral,parai
  USE fint , ONLY:fint1
  USE kpts , ONLY:tkpts
  USE sphe , ONLY:tsphere
  USE phfac_utils, ONLY : calc_eigkr
  USE zeroing_utils,                   ONLY: zeroing
  IMPLICIT NONE
  INTEGER                                    :: ikpt, iflag
  CHARACTER(len=*)                           :: tag

  CHARACTER(*), PARAMETER                    :: procedureN = 'kpt_noswap'

  CHARACTER(len=80)                          :: line
  INTEGER                                    :: i1, ierr, nblkp1
  INTEGER, ALLOCATABLE                       :: lmaskgw(:)
  INTEGER, ALLOCATABLE, SAVE                 :: lalm(:), leigkr(:), lemk2(:), &
                                                lhgkm(:), lhgkp(:), ltwnl(:)
  INTEGER, SAVE                              :: i_first = 0

  i1=LEN(tag)
  line=tag(1:i1)
  IF (i_first.EQ.0) THEN
     nblkp1=nkpt%nblkp+1
     IF (.NOT.tkpts%tkall) THEN
        IF (fint1%ttrot) THEN
           ALLOCATE(lalm(nblkp1),STAT=ierr)
           IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                __LINE__,__FILE__)
           CALL zeroing(lalm)!,nblkp1)
           lalm(0)=1
        ENDIF
        ALLOCATE(lhgkp(nblkp1),STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
             __LINE__,__FILE__)
        CALL zeroing(lhgkp)!,nblkp1)
        lhgkp(0)=1
        ALLOCATE(lhgkm(nblkp1),STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
             __LINE__,__FILE__)
        CALL zeroing(lhgkm)!,nblkp1)
        lhgkm(0)=1
        ALLOCATE(ltwnl(nblkp1),STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
             __LINE__,__FILE__)
        CALL zeroing(ltwnl)!,nblkp1)
        ltwnl(0)=1
        ALLOCATE(lemk2(nblkp1),STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
             __LINE__,__FILE__)
        CALL zeroing(lemk2)!,nblkp1)
        lemk2(0)=1
        IF (tsphere) THEN
           ALLOCATE(lmaskgw(nblkp1),STAT=ierr)
           IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                __LINE__,__FILE__)
           CALL zeroing(lmaskgw)!,nblkp1)
           lmaskgw(0)=1
        ENDIF
     ENDIF
     ALLOCATE(leigkr(nblkp1),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     CALL zeroing(leigkr)!,nblkp1)
     leigkr(0)=1
     i_first=1
  ENDIF
  IF (.NOT.tkpts%tkall) THEN
     ! TODO what is this?.. 
     IF (tsphere.AND.INDEX(line,'MASKGW').NE.0) THEN
     ENDIF
  ENDIF
  IF (INDEX(line,'EIGKR').NE.0) THEN
     IF (iflag.EQ.0) THEN
        IF (leigkr(0).NE.ikpt) THEN
           CALL calc_eigkr(ikpt)
           leigkr(0)=ikpt
        ENDIF
     ELSE
        leigkr(0)=ikpt
     ENDIF
  ENDIF
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE kpt_noswap
! ==================================================================
SUBROUTINE ini_swap(ifile)
  ! ==================================================================
  ! == INITIALIZES SWAP FILES                                       ==
  ! ==================================================================
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE readsr_utils, ONLY : xstring
  USE system , ONLY:nkpt
  USE parac, ONLY : paral,parai
  USE envj , ONLY:my_pid,tmpdir
  USE kpts , ONLY:tkpts
  USE swap , ONLY:accessw,isidisw,iswacce,iswcur,iswend,iswfilemax,iswfirst,iswmax,&
       iswmem,iswpos,iswrec,iswres,iswsize,iunitsw,nswblock,smypid,swfile,tswapc0,tswcalc
  USE zeroing_utils,                   ONLY: zeroing
#if ! defined(_HASNT_F03_ISO_FORTRAN_ENV)
  USE, INTRINSIC :: ISO_FORTRAN_ENV,  ONLY: FILE_STORAGE_SIZE
#endif

  IMPLICIT NONE

#if defined(_HASNT_F03_ISO_FORTRAN_ENV)
  INTEGER, PARAMETER :: FILE_STORAGE_SIZE = 8
#endif

  INTEGER                                    :: ifile

  CHARACTER(*), PARAMETER                    :: procedureN = 'ini_swap'

  CHARACTER(len=1)                           :: acfile
  INTEGER                                    :: i, i1, i2, ierr, iiisw, iisw, &
                                                irecordsize, isub, lentmp, nsw

  CALL tiset('  INI_SWAP',isub)
  ! ==--------------------------------------------------------------==
  IF (ifile.EQ.0) THEN
     IF (iswfirst.NE.0) THEN
        CALL stopgm(' INI_SWAP',' SWAP ALREADY INITIALIZED',& 
             __LINE__,__FILE__)
     ENDIF
     iswfirst=1
     IF (tkpts%tkblock) THEN
        nswblock=nkpt%nblkp
     ELSE
        nswblock=1
     ENDIF
     ! Initialization of list of data in swap files.
     iiisw=iswfilemax
     nsw  =iswfilemax*iswmax
     iisw =iswfilemax*iswmax*nswblock
     ALLOCATE(iunitsw(iiisw),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     ALLOCATE(isidisw(iiisw),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     ALLOCATE(iswres(iiisw),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     ALLOCATE(iswmem(iswmax,iswfilemax),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     ALLOCATE(iswpos(nswblock,iswmax,iswfilemax),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     ALLOCATE(iswsize(nswblock*iswmax,iswfilemax),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     ALLOCATE(iswacce(2,iswmax*nswblock,iswfilemax),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     ALLOCATE(iswcur(iiisw),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     ALLOCATE(iswend(iiisw),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     ALLOCATE(iswrec(iiisw),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     CALL zeroing(iswmem)!,nsw)
     CALL zeroing(iswpos)!,iisw)
     CALL zeroing(iswsize)!,iisw)
     CALL zeroing(iswacce)!,2*iisw)
     lentmp=INDEX(tmpdir, ' ') - 1
     IF (paral%io_parent)&
          WRITE(smypid,'(I6)') my_pid
     CALL xstring(smypid,i1,i2)
     smypid=smypid(i1:i2)
     DO i=1,iswfilemax
        IF (paral%io_parent)&
             WRITE(swfile(i),'(A,A,I1,A,A)')&
             tmpdir(1:lentmp),'/SWAP_',i,'.',smypid
        iswcur(i)=1
        iswend(i)=0
        iswrec(i)=0
        iunitsw(i)=0
        accessw(i)='S'
        isidisw(i)=0
        iswres(i)=0
     ENDDO
     tswcalc=tkpts%tkbcalc
     tswapc0=.NOT.tkpts%tknoswap
     CALL tihalt('  INI_SWAP',isub)
     RETURN
  ENDIF
  ! ==--------------------------------------------------------------==
  ! Special case if we want to calculate k point per k point(TKNOSWAP)
  IF (ifile.EQ.1.AND.(.NOT.tswapc0)) THEN
     ! Do nothing (no swap file for C0).
     CALL tihalt('  INI_SWAP',isub)
     RETURN
  ENDIF
  ! ==--------------------------------------------------------------==
  ! Open swap files
  acfile=accessw(ifile)
  IF (acfile.EQ.'S') THEN
     iunitsw(ifile)=50+ifile
     IF (paral%io_parent)&
          OPEN(unit=iunitsw(ifile),file=swfile(ifile),access='SEQUENTIAL',&
          form='UNFORMATTED',status='NEW')
  ELSEIF (acfile.EQ.'D') THEN
     IF (isidisw(ifile).LE.0)&
          CALL stopgm('INI_SWAP','SIZE OF RECORD NOT SPECIFIED',& 
          __LINE__,__FILE__)
     iunitsw(ifile)=50+ifile
     ! WARNING: THE RECORDSIZE IS MACHINE/COMPILER-DEPENDENT!!
#if defined(__SUN) || defined(LINUX_IFC) || defined( __IBM) || defined(__HP) || defined(__OSX) || defined(__OSX_IFC) ||defined(__WINNT) || defined(__PRIMEHPC) || defined(__HPC)
     ! these machines specify the length in number of characters.
     irecordsize=8*isidisw(ifile)
#elif defined(__SGI) 
     ! these machines specify the length in number of real words.
     irecordsize=2*isidisw(ifile)
#else
     irecordsize= FILE_STORAGE_SIZE *isidisw(ifile)
#endif        
     IF (paral%io_parent)&
          OPEN(unit=iunitsw(ifile),file=swfile(ifile),&
          access='DIRECT',recl=irecordsize,&
          form='UNFORMATTED',status='NEW')
  ELSE
     CALL stopgm('INI_SWAP','ACCESS SEQUENTIAL OR DIRECT',& 
          __LINE__,__FILE__)
  ENDIF
  CALL tihalt('  INI_SWAP',isub)
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE ini_swap
! ==================================================================
SUBROUTINE direct_swap(ifile,isize)
  ! ==--------------------------------------------------------------==
  ! == SPECIFIED FOR A FILE THE DIRECT ACCESS AND THE SIZE OF RECORD==
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE swap , ONLY:accessw,isidisw,iswfilemax,iswfirst,iunitsw
  USE parac , ONLY:paral
  IMPLICIT NONE
  INTEGER                                    :: ifile, isize

  INTEGER                                    :: iunit

  IF ((iswfirst.EQ.0).OR.(ifile.GT.iswfilemax)) THEN
     IF (paral%io_parent)&
          WRITE(6,*) ' SWAP FILE NUMBER:',ifile
     CALL stopgm(' DIRECT_SWAP','SWAP FILE DOESN''T EXIST! !',& 
          __LINE__,__FILE__)
  ENDIF
  iunit=iunitsw(ifile)
  IF (iunit.EQ.0) THEN
     accessw(ifile)='D'
     IF (isize.LE.0) THEN
        IF (paral%io_parent)&
             WRITE(6,*) ' SIZE OF RECORD=',isize
        CALL stopgm(' DIRECT_SWAP',&
             'SIZE OF RECORD HAS TO BE POSITIVE',& 
             __LINE__,__FILE__)
     ELSE
        isidisw(ifile)=isize
     ENDIF
  ELSE
     IF (paral%io_parent)&
          WRITE(6,*) ' SWAP FILE NUMBER:',ifile,' IS ALREADY INITIALIZED'
     CALL stopgm(' DIRECT_SWAP','DIRECT ACCESS IMPOSSIBLE',& 
          __LINE__,__FILE__)
  ENDIF
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE direct_swap
! ==================================================================
SUBROUTINE go_swap(ifile,iposition)
  ! ==--------------------------------------------------------------==
  ! == GO TO IPOSITION IN THE IFILE IF SEQUENTIAL FILE              ==
  ! == OTHERWISE SET ISWCUR(IFILE) TO IPOSITION                     == 
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE swap , ONLY:accessw,iswcur,iswfilemax,iswfirst,iswrec,iunitsw
  USE parac , ONLY:paral
  USE utils, ONLY : fskip
  IMPLICIT NONE
  INTEGER                                    :: ifile, iposition

  INTEGER                                    :: iunit

  IF (iswfirst.EQ.0.OR.ifile.GT.iswfilemax)&
       CALL stopgm(' GO_SWAP', 'BAD SWAP FILE NUMBER',& 
       __LINE__,__FILE__)
  iunit=iunitsw(ifile)
  IF (accessw(ifile).EQ.'S') THEN
     IF (iunit.EQ.0) THEN
        IF (paral%io_parent)&
             WRITE(6,*) ' SWAP FILE NUMBER:',ifile
        CALL stopgm(' GO_SWAP','SWAP FILE DOESN''T EXIST! !',& 
             __LINE__,__FILE__)
     ENDIF
     IF (iposition.GT.iswrec(ifile)+1) THEN
        CALL stopgm(' GO_SWAP', 'BAD SPECIFIED POSITION',& 
             __LINE__,__FILE__)
     ENDIF
     IF (paral%io_parent)&
          REWIND(iunit)
     CALL fskip(iunit,iposition-1)
  ENDIF
  iswcur(ifile)=iposition
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE go_swap
! ==================================================================
SUBROUTINE end_swap
  ! ==--------------------------------------------------------------==
  ! == DEALLOCATE ALL VARIABLES FOR SWAP AND DELETE FILES           ==
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE swap , ONLY:isidisw,iswacce,iswcur,iswend,iswfilemax,iswfirst,iswmem,&
       iswpos,iswrec,iswres,iswsize,isw_calc,iunitsw,i_stopgm,nswblock,tswcalc
  USE parac , ONLY:paral
  IMPLICIT NONE
  CHARACTER(*), PARAMETER                    :: procedureN = 'end_swap'

  INTEGER                                    :: ierr, ifile, ifirst, iunit

  IF (iswfirst.EQ.0) RETURN
  ifirst=0
  DO ifile=1,iswfilemax
     iunit=iunitsw(ifile)
     IF (iunit.NE.0) THEN
        IF (ifirst.EQ.0) THEN
           IF (paral%io_parent)&
                WRITE(6,'(//,1x,64("="))')
           IF (paral%io_parent)&
                WRITE(6,'(T22,A)') 'SWAP FILES INFORMATION'
           ifirst=1
        ENDIF
        CALL state_swap(ifile)
        CALL del_swap(ifile,i_stopgm)
     ELSE
        IF (ifile.EQ.isw_calc.AND.tswcalc) THEN
           CALL calc_swap
        ENDIF
     ENDIF
  ENDDO
  DEALLOCATE(iunitsw,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
       __LINE__,__FILE__)
  DEALLOCATE(iswmem,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
       __LINE__,__FILE__)
  DEALLOCATE(iswpos,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
       __LINE__,__FILE__)
  DEALLOCATE(iswsize,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
       __LINE__,__FILE__)
  DEALLOCATE(iswacce,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
       __LINE__,__FILE__)
  DEALLOCATE(iswcur,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
       __LINE__,__FILE__)
  DEALLOCATE(iswend,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
       __LINE__,__FILE__)
  DEALLOCATE(iswrec,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
       __LINE__,__FILE__)
  DEALLOCATE(isidisw,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
       __LINE__,__FILE__)
  DEALLOCATE(iswres,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
       __LINE__,__FILE__)
  nswblock=0
  iswfirst=0
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE end_swap
! ==================================================================
SUBROUTINE del_swap(ifile,iflag)
  ! ==--------------------------------------------------------------==
  ! == DELETE SWAP FILE IFILE                                       ==
  ! == IFLAG IS TO AVOID RECURSIVE CALL OF STOPGM                   ==
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE machine, ONLY: m_system
  USE swap , ONLY:accessw,isidisw,iswfilemax,iswfirst,iunitsw,i_nostopgm,swfile
  USE parac , ONLY:paral
  IMPLICIT NONE
  INTEGER                                    :: ifile, iflag

  CHARACTER(len=280)                         :: order
  INTEGER                                    :: iunit

  IF (iflag.EQ.i_nostopgm) THEN
     IF (iswfirst.EQ.0.OR.ifile.GT.iswfilemax)&
          CALL stopgm(' DEL_SWAP', 'BAD SWAP FILE NUMBER',& 
          __LINE__,__FILE__)
  ENDIF
  iunit=iunitsw(ifile)
  IF (iflag.EQ.i_nostopgm) THEN
     IF (iunit.EQ.0) THEN
        IF (paral%io_parent)&
             WRITE(6,*) ' SWAP FILE NUMBER:',ifile
        CALL stopgm(' DEL_SWAP','SWAP FILE DOESN''T EXIST! !',& 
             __LINE__,__FILE__)
     ENDIF
  ENDIF
  IF (paral%io_parent)&
       CLOSE(unit=iunit)
#if defined(__SUN)
  order = '\\rm -f '//swfile(ifile)
#else
  order = 'rm -f '//swfile(ifile)
#endif
  CALL m_system(order)
  CALL res_swap(ifile,iflag)
  iunitsw(ifile)=0
  accessw(ifile)='S'
  isidisw(ifile)=0
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE del_swap
! ==================================================================
SUBROUTINE res_swap(ifile,iflag)
  ! ==--------------------------------------------------------------==
  ! == RESET COUNTERS OF SWAP FILE IFILE                            ==
  ! == NOT CLOSE THE FILE                                           ==
  ! == IFLAG IS TO AVOID RECURSIVE CALL OF STOPGM                   ==
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE swap , ONLY:iswacce,iswcur,iswend,iswfilemax,iswfirst,iswmem,iswpos,&
       iswrec,iswres,iswsize,iunitsw,i_nostopgm
  USE parac , ONLY:paral
  USE zeroing_utils,                   ONLY: zeroing
  IMPLICIT NONE
  INTEGER                                    :: ifile, iflag

  INTEGER                                    :: iunit

  IF (iflag.EQ.i_nostopgm) THEN
     IF (iswfirst.EQ.0.OR.ifile.GT.iswfilemax)&
          CALL stopgm(' RES_SWAP', 'BAD SWAP FILE NUMBER',& 
          __LINE__,__FILE__)
  ENDIF
  iunit=iunitsw(ifile)
  IF (iflag.EQ.i_nostopgm) THEN
     IF (iunit.EQ.0) THEN
        IF (paral%io_parent)&
             WRITE(6,*) ' SWAP FILE NUMBER:',ifile
        CALL stopgm(' RES_SWAP','SWAP FILE DOESN''T EXIST! !',& 
             __LINE__,__FILE__)
     ENDIF
  ENDIF
  CALL zeroing(iswmem(:,ifile))!   ,iswmax)
  CALL zeroing(iswpos(:,:,ifile))! ,  nswblock*iswmax)
  CALL zeroing(iswsize(:,ifile))!  ,  nswblock*iswmax)
  CALL zeroing(iswacce(:,:,ifile))!,2*nswblock*iswmax)
  iswcur(ifile)=1
  iswend(ifile)=0
  iswrec(ifile)=0
  iswres(ifile)=iswres(ifile)+1
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE res_swap
! ==================================================================
SUBROUTINE state_swap(ifile)
  ! ==--------------------------------------------------------------==
  ! == STATE OF SWAP FILE IFILE                                     ==
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE readsr_utils, ONLY : xstring
  USE swap , ONLY:accessw,isidisw,iswacce,iswcur,iswend,iswfilemax,iswfirst,&
       iswmem,iswpos,iswrec,iswres,iswsize,iunitsw,nswblock,swfile,swname
  USE parac , ONLY:paral
  IMPLICIT NONE
  INTEGER                                    :: ifile

  CHARACTER(len=15)                          :: fformat
  INTEGER                                    :: i1, i2, ib, in, iposition, &
                                                iunit

  IF (iswfirst.EQ.0) THEN
     IF (paral%io_parent)&
          WRITE(6,*) ' NO SWAP FILE'
     RETURN
  ENDIF
  IF (ifile.GT.iswfilemax.OR.ifile.LE.0) THEN
     IF (paral%io_parent)&
          WRITE(6,'(A,I6,T51,A15)')&
          ' SWAP FILE NUMBER:',ifile,'DOESN''T EXIST! !'
     RETURN
  ENDIF
  iunit=iunitsw(ifile)
  IF (iunit.EQ.0) THEN
     IF (paral%io_parent)&
          WRITE(6,'(A,I6,T51,A15)')&
          ' SWAP FILE NUMBER:',ifile,'DOESN''T EXIST! !'
     RETURN
  ENDIF
  ! ==--------------------------------------------------------------==
  IF (paral%io_parent)&
       WRITE(6,'(1x,64("="))')
  CALL xstring(swfile(ifile),i1,i2)
  IF (65-(i2-i1).LE.18) THEN
     IF (paral%io_parent)&
          WRITE(fformat,'(A,I2,A)') '(A,/,T',MAX(1,65-(i2-i1)),',A)'
  ELSE
     IF (paral%io_parent)&
          WRITE(fformat,'(A,I2,A)') '(A,T',MAX(18,65-(i2-i1)),',A)'
  ENDIF
  IF (paral%io_parent)&
       WRITE(6,fformat)    ' SWAP FILE NAME: ',swfile(ifile)(i1:i2)
  IF (paral%io_parent)&
       WRITE(6,'(A,T60,I6)') ' SWAP FILE NUMBER:',ifile
  IF (paral%io_parent)&
       WRITE(6,'(A,T60,I6)') ' SWAP FILE UNIT:',iunit
  IF (accessw(ifile).EQ.'S') THEN
     IF (paral%io_parent)&
          WRITE(6,'(" SWAP FILE ACCESS:",T56,"SEQUENTIAL")')
  ELSE
     IF (paral%io_parent)&
          WRITE(6,'(" SWAP FILE ACCESS:",T56,"    DIRECT")')
     IF (paral%io_parent)&
          WRITE(6,'(A,T46,I20)') ' RECORD SIZE:',isidisw(ifile)
  ENDIF
  IF ((iswres(ifile).NE.0).AND.paral%io_parent)&
       WRITE(6,'(A,T60,I6)') ' NUMBER OF RESET:',iswres(ifile)
  IF (paral%io_parent)&
       WRITE(6,'(A,T60,I6)') ' NUMBER OF RECORDS:',iswrec(ifile)
  IF (paral%io_parent)&
       WRITE(6,'(A,T60,I6)') ' CURRENT RECORD:', iswcur(ifile)
  IF (paral%io_parent)&
       WRITE(6,'(A,T60,I6)') ' NUMBER OF VARIABLES:',iswend(ifile)
  IF (paral%io_parent)&
       WRITE(6,'(A,T60,I6)') ' NUMBER OF BLOCK PER VARIABLE:', nswblock
  DO in=1,iswend(ifile)
     CALL xstring(swname(in,ifile),i1,i2)
     IF (paral%io_parent)&
          WRITE(6,'(I6,1X,A,T44,A18,I3,A1)')&
          in,swname(in,ifile)(i1:i2),&
          '(IN MEMORY BLOCK=',iswmem(in,ifile),')'
     DO ib=1,nswblock
        IF (iswpos(ib,in,ifile).NE.0) THEN
           iposition=iswpos(ib,in,ifile)
           IF (paral%io_parent)&
                WRITE(6,'(5X,A,I5,T18,A,I5,T29,A,I8,T43,A,I5,T56,A,I5)')&
                'BLOCK=' ,ib,&
                'REC='   ,iposition,&
                'SIZE='  ,iswsize(iposition,ifile),&
                'WRITE=' ,iswacce(1,iposition,ifile),&
                'READ='  ,iswacce(2,iposition,ifile)
        ENDIF
     ENDDO
  ENDDO
  IF (paral%io_parent)&
       WRITE(6,'(1x,64("="))')
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE state_swap
! ==================================================================
SUBROUTINE calc_swap
  ! ==--------------------------------------------------------------==
  ! == GIVE THE NUMBER OF CALCULATION IF TSWCALC                    ==
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE swap , ONLY:iswacce,isw_alm,isw_calc,isw_eigkr,isw_emk2,isw_hgkm,isw_hgkp,isw_maskgw,isw_twnl,tswcalc
  USE parac , ONLY:paral
  IMPLICIT NONE
! ==--------------------------------------------------------------==

  IF (.NOT.tswcalc) RETURN
  IF (paral%io_parent)&
       WRITE(6,'(1x,64("="))')
  IF (paral%io_parent)&
       WRITE(6,'(1X,A)')&
       'NUMBER OF CALCULATIONS FOR ARRAY DEPENDING ON K POINTS:'
  IF ((iswacce(2,isw_hgkp,isw_calc).NE.0).AND.paral%io_parent)&
       WRITE(6,'(1X,A,T60,I6)')&
       'HGKP  ',iswacce(2,isw_hgkp,isw_calc)
  IF ((iswacce(2,isw_hgkm,  isw_calc).NE.0).AND.paral%io_parent)&
       WRITE(6,'(1X,A,T60,I6)')&
       'HGKM  ',iswacce(2,isw_hgkm,isw_calc)
  IF ((iswacce(2,isw_maskgw,isw_calc).NE.0).AND.paral%io_parent)&
       WRITE(6,'(1X,A,T60,I6)')&
       'MASKGW',iswacce(2,isw_twnl,isw_calc)
  IF ((iswacce(2,isw_twnl,  isw_calc).NE.0).AND.paral%io_parent)&
       WRITE(6,'(1X,A,T60,I6)')&
       'TWNL  ',iswacce(2,isw_twnl,isw_calc)
  IF ((iswacce(2,isw_emk2,  isw_calc).NE.0).AND.paral%io_parent)&
       WRITE(6,'(1X,A,T60,I6)')&
       'EMK2  ',iswacce(2,isw_emk2,isw_calc)
  IF ((iswacce(2,isw_eigkr, isw_calc).NE.0).AND.paral%io_parent)&
       WRITE(6,'(1X,A,T60,I6)')&
       'EIGKR ',iswacce(2,isw_eigkr,isw_calc)
  IF ((iswacce(2,isw_alm,   isw_calc).NE.0).AND.paral%io_parent)&
       WRITE(6,'(1X,A,T60,I6)')&
       'ALM   ',iswacce(2,isw_alm,isw_calc)
  IF (paral%io_parent)&
       WRITE(6,'(1x,64("="))')
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE calc_swap
! ==================================================================
SUBROUTINE inq_swap(kbeg,kend,kinc)
  ! ==--------------------------------------------------------------==
  ! == GIVE THE ORDER TO PERFORM LOOP ON K POINTS BLOCK             ==
  ! == IN ORDER TO REDUCE READ NUMBER                               ==
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE system , ONLY:nkpt
  USE parac, ONLY : paral,parai
  USE kpts , ONLY:tkpts
  USE swap , ONLY:iswfirst,iswmem
  IMPLICIT NONE
  INTEGER                                    :: kbeg, kend, kinc

  IF (.NOT.tkpts%tkblock .OR. iswfirst.EQ.0) THEN
     kbeg=1
     kend=nkpt%nblkp
     kinc=1
     RETURN
  ENDIF
  IF (iswmem(1,1).EQ.nkpt%nblkp) THEN
     kbeg=nkpt%nblkp
     kend=1
     kinc=-1
  ELSE
     kbeg=1
     kend=nkpt%nblkp
     kinc=1
  ENDIF
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE inq_swap
! ==================================================================
SUBROUTINE reskpt_swap
  ! ==--------------------------------------------------------------==
  ! == RESET SWAP FILE IF TKALL                                     ==
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE swap , ONLY:i_nostopgm
  USE kpts , ONLY:tkpts
  IMPLICIT NONE
! ==--------------------------------------------------------------==

  IF (tkpts%tkall) THEN
     CALL res_swap(3,i_nostopgm)
  ENDIF
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE reskpt_swap
! ==================================================================

SUBROUTINE rea_swap(ifile,a,na,ikpt,name)
  ! ==--------------------------------------------------------------==
  ! == READ IN SWAP FILES                                           ==
  ! == IFILE=0 SWAP PERMANENT FILE                                  ==
  ! ==       1 SWAP TEMPORARY FILE                                  ==
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE readsr_utils, ONLY : xstring
  USE swap , ONLY:accessw,iswacce,iswcur,iswend,iswfilemax,iswfirst,&
       iswmem,iswpos,iswsize,iunitsw,nswblock,swname
  USE parac , ONLY:paral
  IMPLICIT NONE
  INTEGER                                    :: ifile, na
  REAL(real_8)                               :: a(na)
  INTEGER                                    :: ikpt
  CHARACTER(len=*)                           :: name

  CHARACTER(len=8)                           :: sname
  INTEGER                                    :: in1, in2, indi, iposition, &
                                                is1, is2, isub, iunit

  CALL tiset('  REA_SWAP',isub)
  ! ==--------------------------------------------------------------==
  in2=LEN(name)
  CALL xstring(name,in1,in2)
  IF (iswfirst.EQ.0.OR.ifile.GT.iswfilemax) THEN
     IF (paral%io_parent)&
          WRITE(6,*) ' SWAP FILE NUMBER:',ifile
     CALL stopgm(' REA_SWAP','SWAP FILE DOESN''T EXIST! !',& 
          __LINE__,__FILE__)
  ENDIF
  iunit=iunitsw(ifile)
  IF (iunit.EQ.0) THEN
     IF (paral%io_parent)&
          WRITE(6,*) ' SWAP FILE NUMBER:',ifile
     CALL stopgm(' REA_SWAP','SWAP FILE DOESN''T EXIST! !',& 
          __LINE__,__FILE__)
  ENDIF
  IF (ikpt.GT.nswblock) THEN
     CALL stopgm(' REA_SWAP','INDEX BLOCK WRONG',& 
          __LINE__,__FILE__)
  ENDIF
  DO indi=1,iswend(ifile)
     sname=swname(indi,ifile)
     CALL xstring(sname,is1,is2)
     IF (sname(is1:is2).EQ.name(in1:in2)) GOTO 30
  ENDDO
  IF ( (iswend(ifile).EQ.0).OR.(indi.GT.iswend(ifile)) ) THEN
     IF (paral%io_parent)&
          WRITE(6,*) ' SWAP FILE NUMBER:',ifile
     IF (paral%io_parent)&
          WRITE(6,*) ' ',name(in1:in2), ' BLOCK=',ikpt
     CALL stopgm(' REA_SWAP', 'READ ERROR',& 
          __LINE__,__FILE__)
  ENDIF
30 CONTINUE
  IF (iswmem(indi,ifile).EQ.ikpt) THEN
     CALL tihalt('  REA_SWAP',isub)
     RETURN
  ENDIF
  iposition=iswpos(ikpt,indi,ifile)
  IF (iposition.EQ.0) THEN
     IF (paral%io_parent)&
          WRITE(6,*) ' SWAP FILE NUMBER:',ifile
     IF (paral%io_parent)&
          WRITE(6,*) ' ',name(in1:in2), ' BLOCK=',ikpt
     CALL stopgm(' REA_SWAP', 'RECORD DOESN''T EXIST! !',& 
          __LINE__,__FILE__)
  ENDIF
  IF (iswsize(iposition,ifile).NE.na) THEN
     IF (paral%io_parent)&
          WRITE(6,*) ' SWAP FILE NUMBER:',ifile
     IF (paral%io_parent)&
          WRITE(6,*) ' ',name(in1:in2), ' BLOCK=',ikpt
     CALL stopgm(' REA_SWAP', 'ARRAY SIZE IS WRONG',& 
          __LINE__,__FILE__)
  ENDIF
  IF (accessw(ifile).EQ.'S') THEN
     CALL go_swap(ifile,iposition)
     IF (paral%io_parent)&
          READ(unit=iunit,err=999,END=998) a
     iswcur(ifile)=iswcur(ifile)+1
  ELSE
     IF (paral%io_parent)&
          READ(unit=iunit,err=999,rec=iposition) a
     iswcur(ifile)=iposition+1
  ENDIF
  iswacce(2,iposition,ifile)=iswacce(2,iposition,ifile)+1
  iswmem(indi,ifile)=ikpt
  CALL tihalt('  REA_SWAP',isub)
  ! ==--------------------------------------------------------------==
  RETURN
998 CONTINUE
  CALL stopgm(' REA_SWAP', 'END OF SWAP FILE',& 
       __LINE__,__FILE__)
999 CONTINUE
  CALL stopgm(' REA_SWAP', 'CAN''T READ IN SWAP FILE',& 
       __LINE__,__FILE__)
END SUBROUTINE rea_swap
! ==================================================================
SUBROUTINE wri_swap(ifile,a,na,ikpt,name)
  ! ==--------------------------------------------------------------==
  ! == WRITE IN SWAP FILES                                          ==
  ! == IFILE=1 SWAP PERMANENT FILE                                  ==
  ! ==       2 SWAP TEMPORARY FILE                                  ==
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE swap , ONLY:accessw,isidisw,iswacce,iswcur,iswend,iswfilemax,iswfirst,&
       iswmax,iswmem,iswpos,iswrec,iswsize,iunitsw,nswblock,swname
  USE parac , ONLY:paral
  USE readsr_utils, ONLY : xstring
  IMPLICIT NONE
  INTEGER                                    :: ifile, na
  REAL(real_8)                               :: a(na)
  INTEGER                                    :: ikpt
  CHARACTER(len=*)                           :: name

  CHARACTER(len=1)                           :: acfile
  CHARACTER(len=8)                           :: sname
  INTEGER                                    :: in1, in2, indi, iposition, &
                                                is1, is2, isub, iunit

  CALL tiset('  WRI_SWAP',isub)
  ! ==--------------------------------------------------------------==
  in2=LEN(name)
  CALL xstring(name,in1,in2)
  IF ((iswfirst.EQ.0).OR.(ifile.GT.iswfilemax)) THEN
     IF (paral%io_parent)&
          WRITE(6,*) ' SWAP FILE NUMBER:',ifile
     CALL stopgm(' WRI_SWAP','SWAP FILE DOESN''T EXIST! !',& 
          __LINE__,__FILE__)
  ENDIF
  iunit=iunitsw(ifile)
  acfile=accessw(ifile)
  IF (iunit.EQ.0) THEN
     CALL ini_swap(ifile)
     iunit=iunitsw(ifile)
  ENDIF
  IF (ikpt.GT.nswblock) THEN
     CALL stopgm(' WRI_SWAP','INDEX IN K POINTS BLOCK WRONG',& 
          __LINE__,__FILE__)
  ENDIF
  DO indi=1,iswend(ifile)
     sname=swname(indi,ifile)
     CALL xstring(sname,is1,is2)
     IF (sname(is1:is2).EQ.name(in1:in2)) GOTO 30
  ENDDO
30 CONTINUE
  IF ( (iswend(ifile).EQ.0).OR.(indi.GT.iswend(ifile)) ) THEN
     iswend(ifile)=iswend(ifile)+1
     IF (iswend(ifile).GT.iswmax) THEN
        CALL stopgm('WRI_SWAP','INCREASE ISWMAX',& 
             __LINE__,__FILE__)
     ENDIF
     swname(iswend(ifile),ifile)=name(in1:in2)//' '
     indi=iswend(ifile)
  ENDIF
  iposition=iswpos(ikpt,indi,ifile)
  IF (iposition.NE.0) THEN
     ! Possible to write if the size is the same.
     IF (acfile.EQ.'S') THEN
        IF (paral%io_parent)&
             WRITE(6,*) ' SWAP FILE NUMBER:',ifile
        IF (paral%io_parent)&
             WRITE(6,*) ' TYPE OF ACCESS: SEQUENTIAL'
        IF (paral%io_parent)&
             WRITE(6,*) ' REWRITE FOR ',name,' (IKPT=',ikpt,')'
        CALL stopgm(' WRI_SWAP', 'REWRITE IS FORBIDDEN! ',& 
             __LINE__,__FILE__)
     ELSE
        IF (iswsize(iposition,ifile).NE.na) THEN
           IF (paral%io_parent)&
                WRITE(6,*) ' SWAP FILE NUMBER:',ifile
           IF (paral%io_parent)&
                WRITE(6,*) ' TYPE OF ACCESS: DIRECT'
           IF (paral%io_parent)&
                WRITE(6,*) ' REWRITE FOR ',name,' (IKPT=',ikpt,')'
           IF (paral%io_parent)&
                WRITE(6,*) ' OLD SIZE=',iswsize(iposition,ifile)
           IF (paral%io_parent)&
                WRITE(6,*) ' NEW SIZE=',na
           CALL stopgm(' WRI_SWAP','WRONG SIZE OF BLOCK',& 
                __LINE__,__FILE__)
        ENDIF
        CALL go_swap(ifile,iposition)
     ENDIF
  ELSE
     IF ( (acfile.EQ.'D').AND.(na.GT.isidisw(ifile)) ) THEN
        IF (paral%io_parent)&
             WRITE(6,*) ' SWAP FILE NUMBER:',ifile
        IF (paral%io_parent)&
             WRITE(6,*) ' TYPE OF ACCESS: DIRECT'
        IF (paral%io_parent)&
             WRITE(6,*) ' REWRITE FOR ',name,' (IKPT=',ikpt,')'
        IF (paral%io_parent)&
             WRITE(6,*) ' RECORD SIZE=',isidisw(ifile)
        IF (paral%io_parent)&
             WRITE(6,*) ' WANTED SIZE=',na
        CALL stopgm('WRI_SWAP','INCREASE RECORD SIE',& 
             __LINE__,__FILE__)
     ENDIF
     CALL go_swap(ifile,iswrec(ifile)+1)
     iswrec(ifile)=iswrec(ifile)+1
     iposition=iswrec(ifile)
  ENDIF
  IF (acfile.EQ.'S') THEN
     IF (paral%io_parent)&
          WRITE(unit=iunit,err=998) a
  ELSE
     ! WARNING: THE RECORDSIZE IS PROCESSOR-DEPENDENT!!
     ! (SEE INI_SWAP)
     IF (paral%io_parent)&
          WRITE(unit=iunit,err=999,rec=iposition) a
  ENDIF
  iposition=iswcur(ifile)
  iswpos(ikpt,indi,ifile)=iposition
  iswsize(iposition,ifile)=na
  iswacce(1,iposition,ifile)=iswacce(1,iposition,ifile)+1
  iswcur(ifile)=iswcur(ifile)+1
  iswmem(indi,ifile)=ikpt
  CALL tihalt('  WRI_SWAP',isub)
  ! ==--------------------------------------------------------------==
  RETURN
998 CONTINUE
  CALL stopgm(' WRI_SWAP', 'CAN''T WRITE IN SEQUENTIAL SWAP FILE',& 
       __LINE__,__FILE__)
999 CONTINUE
  CALL stopgm(' WRI_SWAP', 'CAN''T WRITE IN DIRECT SWAP FILE',& 
       __LINE__,__FILE__)
END SUBROUTINE wri_swap
! ==================================================================

