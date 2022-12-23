MODULE rkpnt_utils
  USE coor,                            ONLY: tau0
  USE cppt,                            ONLY: hg,&
                                             inyh
  USE envj,                            ONLY: hname,&
                                             my_pid,&
                                             tmpdir,&
                                             user
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_def
  USE gvec,                            ONLY: gvec_com
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE k290_utils,                      ONLY: k290
  USE kinds,                           ONLY: real_8
  USE kpnt,                            ONLY: hgkm,&
                                             hgkp,&
                                             rk,&
                                             wk
  USE kpts,                            ONLY: kpts_com,&
                                             tkpts,&
                                             wvk0
  USE machine,                         ONLY: m_datum,&
                                             m_getarg
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai,&
                                             paral
  USE prmem_utils,                     ONLY: prmem
  USE prng_utils,                      ONLY: repprngu
  USE readsr_utils,                    ONLY: xstring
  USE sphe,                            ONLY: gcutwmin,&
                                             maskgw,&
                                             maskl,&
                                             ngw_gamma,&
                                             ngwkps,&
                                             tsphere
  USE symm,                            ONLY: symmr
  USE system,                          ONLY: kpbeg,&
                                             ncpw,&
                                             nkpbl,&
                                             nkpt,&
                                             parm,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rkpnt
  PUBLIC :: bzmesh
  !public :: setmask
  PUBLIC :: calc_maskgw
  PUBLIC :: calc_hgk

CONTAINS

  ! ==================================================================
  SUBROUTINE rkpnt
    ! ==================================================================
    ! == ROUTINE TO SET UP K-POINT MESH AND                           ==
    ! == THE ARRAYS HGKP AND HGKM WHISH HOLD                          ==
    ! == THE QUANTITIES |k+g|^2 and |k-g|^2                           ==
    ! ==================================================================
    ! Variables

    CHARACTER(*), PARAMETER                  :: procedureN = 'rkpnt'

    CHARACTER(len=9)                         :: fformat
    INTEGER                                  :: i1, i2, ierr, ig, ik, ikk, &
                                                ikpt, ir, isub, k
    REAL(real_8)                             :: t, temp(3)

! ==--------------------------------------------------------------==

    CALL tiset('     RKPNT',isub)
    ! ==--------------------------------------------------------------==
    ! == ALLOCATE MEMORY FOR THE PERMANENT ARRAYS                     ==
    ! ==--------------------------------------------------------------==
    ! Set NKPNT = number of kpoints in memory.
    IF (tkpts%tkblock) THEN
       IF (nkpt%nkpnt.GE.nkpt%nkpts) THEN
          nkpt%nkpnt=nkpt%nkpts
          tkpts%tkblock=.FALSE.
          tkpts%tkall=.FALSE.
          tkpts%tknoswap=.FALSE.
          tkpts%tkbcalc=.FALSE.
          nkpt%nblkp=1
       ELSE
          nkpt%nblkp=nkpt%nkpts/nkpt%nkpnt
          IF (MOD(nkpt%nkpts,nkpt%nkpnt).NE.0) nkpt%nblkp=nkpt%nblkp+1
       ENDIF
    ELSE
       nkpt%nkpnt=nkpt%nkpts
       nkpt%nblkp=1
    ENDIF
    ALLOCATE(nkpbl(nkpt%nblkp),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(kpbeg(nkpt%nblkp),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ik=0
    DO ikpt=1,nkpt%nblkp
       kpbeg(ikpt)=ik
       nkpbl(ikpt)=MIN(nkpt%nkpts,nkpt%nkpnt*ikpt)-nkpt%nkpnt*(ikpt-1)
       ik=ik+nkpbl(ikpt)
    ENDDO
    ! ==--------------------------------------------------------------==
    ! ALLOCATE MEMORY FOR THE PERMANENT ARRAYS
    IF (tkpts%tkblock.AND.(tkpts%tkall.OR.tkpts%tkbcalc)) THEN
       kpts_com%nkptall=nkpt%nkpnt
    ELSE
       kpts_com%nkptall=nkpt%nkpts
    ENDIF
    ALLOCATE(hgkp(ncpw%nhg,kpts_com%nkptall),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(hgkm(ncpw%nhg,kpts_com%nkptall),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(hgkp)!,nhg*kpts_com%nkptall)
    CALL zeroing(hgkm)!,nhg*kpts_com%nkptall)
    ! ==--------------------------------------------------------------==
    ! IF USED BLOCK OPTION: MESSAGE
    IF (paral%parent.AND.tkpts%tkblock) THEN
       IF (.NOT.tkpts%tkbcalc.AND..NOT.tkpts%tknoswap) THEN
          CALL stopgm('RKPNT',&
               'BLOCK OPTION NOT AVAILABLE ON PARALLEL MACHINE',& 
               __LINE__,__FILE__)
       ENDIF
       IF (paral%io_parent)&
            WRITE(6,'(/,1x,64("*"))')
       IF (paral%io_parent)&
            WRITE(6,'(A,T60,I6)')&
            ' TOTAL NUMBER OF K POINTS:',nkpt%nkpts
       IF (paral%io_parent)&
            WRITE(6,'(A,T60,I6)')&
            ' CALCULATING K POINTS PER BLOCK OF:',nkpt%nkpnt
       IF (paral%io_parent)&
            WRITE(6,'(A,T60,I6)')&
            ' NUMBER OF BLOCKS:',nkpt%nblkp
       IF (.NOT.tkpts%tknoswap.AND..NOT.tkpts%tkbcalc) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A)')   ' USING SWAP FILES IN DIRECTORY: '
          CALL xstring(tmpdir,i1,i2)
          IF (paral%io_parent)&
               WRITE(fformat,'(A,I2,A)') '(T',MAX(2,65-(i2-i1)),',A)'
          IF (paral%io_parent)&
               WRITE(6,fformat) tmpdir(i1:i2)
       ENDIF
       IF (tkpts%tknoswap) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A)')&
               ' WAVEFUNCTIONS ARE NOT RELEVANT [NOT SAVED].'
       ELSE
          IF (paral%io_parent)&
               WRITE(6,'(A)')&
               ' WAVEFUNCTIONS ARE STORED IN SWAP FILE.'
       ENDIF
       IF (tkpts%tkall) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,A)') ' ALL ARRAYS DEPENDING ON K POINTS ',&
               'ARE STORED IN SWAP FILE.'
       ELSEIF (tkpts%tkbcalc) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,A)')&
               ' ALL ARRAYS DEPENDING ON K POINTS ',&
               'ARE CALCULATED EACH TIME.'
       ENDIF
       IF (paral%io_parent)&
            WRITE(6,'(A)')&
            ' PHASE FACTORS ARE CALCULATED EACH TIME.'
       IF (paral%io_parent)&
            WRITE(6,'(1x,64("*"),/)')
    ENDIF
    ! ==--------------------------------------------------------------==
    ! SCALED OPTION with K-points
    IF (tkpts%tkscale) THEN
       DO ikk=1,nkpt%nkpts
          DO ir=1,3
             temp(ir)=rk(1,ikk)*gvec_com%b1(ir)+rk(2,ikk)*gvec_com%b2(ir)+rk(3,ikk)*gvec_com%b3(ir)
          ENDDO
          rk(1,ikk)=temp(1)
          rk(2,ikk)=temp(2)
          rk(3,ikk)=temp(3)
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    DO ikpt=1,nkpt%nblkp
       CALL calc_hgk(ikpt)
       IF (tkpts%tkblock) CALL wkpt_swap(hgkp,1,ikpt,'HGKP HGKM')
    ENDDO
    ! ==--------------------------------------------------------------==
    ! SET UP ARRAYS FOR SPHERICAL CUTOFF OPTION
    ALLOCATE(ngwkps(nkpt%nkpts),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    IF (tsphere) THEN
       DO ig=1,ncpw%ngw
          IF (hg(ig).GT.gcutwmin) GOTO 1000
       ENDDO
1000   CONTINUE
       maskl=ncpw%ngw-ig+1
       ALLOCATE(maskgw(2,maskl,kpts_com%nkptall),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL setmask(gvec_com%gcutw)
    ELSE
       DO ik=1,nkpt%nkpts
          ngwkps(ik)=spar%ngwks
       ENDDO
    ENDIF
    CALL tihalt('     RKPNT',isub)
    IF (paral%parent) CALL prmem('     RKPNT')
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rkpnt
  ! ==================================================================
  SUBROUTINE bzmesh(nkpoint)
    ! ==--------------------------------------------------------------==
    ! == GENERATE A UNIFORM (MONKHORST-PACK) MESH IN THE KX>0 HALF OF ==
    ! == THE BRILLOUIN ZONE.                                          ==
    ! == REF: HJ MONKHORST AND JD PACK, PRB 13, 5188, (1976).         ==
    ! == USE K290 ROUTINE TO GENERATE SPECIAL K-POINTS                ==
    ! == USE INTERNAL FILES TO COMMUNICATE WITH K290 INSTEAD OF FORT.x==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nkpoint

    CHARACTER(*), PARAMETER                  :: procedureN = 'bzmesh'

    CHARACTER(len=26)                        :: datx
    CHARACTER(len=30)                        :: fformat
    CHARACTER(len=80)                        :: filename
    INTEGER                                  :: i, ia, ierr, ih1, ih2, ikind, &
                                                iout, is, istriz, isum, iu1, &
                                                iu2, l, nat0, nhash, ntvect
    INTEGER, ALLOCATABLE                     :: f0(:,:), includ(:), isc(:), &
                                                list(:), lrot(:,:), lwght(:), &
                                                ty(:)
    LOGICAL                                  :: ferror
    REAL(real_8)                             :: strain(6), sum, weight, &
                                                xrandom(3)
    REAL(real_8), ALLOCATABLE                :: rlist(:,:), rx(:,:), &
                                                tvect(:,:), wvkl(:,:), &
                                                xkapa(:,:)

! Variables
! ==--------------------------------------------------------------==
! ==--------------------------------------------------------------==
! LATTICE STRUCTURE
! SPECIAL POINTS
! Used by SPPT2 and MESH routines.
! ==--------------------------------------------------------------==
! ==--------------------------------------------------------------==
! CALL TISET('    BZMESH',ISUB)
! ==--------------------------------------------------------------==
! == Standard output of K290                                      ==
! ==--------------------------------------------------------------==

    iout=91
    IF (paral%io_parent)&
         CALL fileopen(iout,'KPTS_GENERATION',fo_def,ferror)
    IF (paral%io_parent)&
         WRITE(iout,'(75("*"),/," K290 - SPECIAL K POINTS GENERATION - ")')
    CALL m_datum(datx)
    CALL xstring(hname,ih1,ih2)
    IF (paral%io_parent)&
         WRITE(iout,'(75("*"),/," COMPUTER: ",A," AT ",A26,/)')&
         hname(ih1:ih2),datx
    CALL xstring(user,iu1,iu2)
    IF (paral%io_parent)&
         WRITE(iout,&
         '(" CPMD JOB - (",A,") PID= ",I7," ON ",A,10(" "))')&
         user(iu1:iu2), my_pid, hname(ih1:ih2)
    CALL m_getarg(1,filename)
    CALL xstring(filename,iu1,iu2)
    IF (paral%io_parent)&
         WRITE(fformat,'(A,I2,A)') '(A,T',MAX(21,45-(iu2-iu1)),',A,/)'
    IF (paral%io_parent)&
         WRITE(iout,fformat)   ' FROM THE INPUT FILE : ',filename(iu1:iu2)
    ! ==--------------------------------------------------------------==
    ! == Allocations of arrays                                        ==
    ! ==--------------------------------------------------------------==
    ALLOCATE(xkapa(3,ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(rx(3,ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(tvect(3,ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(ty(ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(isc(ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(f0(49,ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(wvkl(3,nkpoint),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(lwght(nkpoint),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(lrot(48,nkpoint),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(includ(nkpoint),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    nhash = MAX(2000,nkpoint/10)
    ALLOCATE(list((nkpoint+nhash)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(rlist(3,nkpoint),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! == INPUT DATA FOR K290                                          ==
    ! ==--------------------------------------------------------------==
    ! Symmetrization of Monkhorst-Pack mesh
    IF (tkpts%tsymkp) THEN
       istriz = 1
    ELSE
       istriz = -1
    ENDIF
    ! Strain
    DO i = 1,6
       strain(i) = 0._real_8
    ENDDO
    ! If option FULL (symmetry with disorder atoms).
    DO i=1,3
       xrandom(i)=0._real_8
    ENDDO
    nat0=0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          IF (tkpts%tkfull) THEN
             DO i=1,3
                xrandom(i)=0.1_real_8*repprngu()
             ENDDO
          ENDIF
          nat0=nat0+1
          ty(nat0) = is
          DO i=1,3
             xkapa(i,nat0) = tau0(i,ia,is)+xrandom(i)*parm%alat
          ENDDO
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent)&
         WRITE(6,'(/,1X,19("*"),A,16("*"))')&
         ' SPECIAL K-POINTS GENERATION '
    IF (paral%io_parent)&
         WRITE(6,'(" DIMENSIONS ARE:")')
    IF (paral%io_parent)&
         WRITE(6,'(3X," NUMBER OF ATOMS",T58,I8)') ions1%nat
    IF (paral%io_parent)&
         WRITE(6,'(3X," K POINTS MONKHORST-PACK MESH",T42,3I8)')&
         kpts_com%nk1,kpts_com%nk2,kpts_com%nk3
    IF (paral%io_parent)&
         WRITE(6,'(3X," MAXIMAL NUMBER OF K POINTS",T58,I8)') nkpoint
    IF (paral%io_parent)&
         WRITE(6,'(3X," CONSTANT VECTOR SHIFT (MACDONALD)",T42,3F8.3)')&
         wvk0
    IF (tkpts%tsymkp) THEN
       IF (paral%io_parent)&
            WRITE(6,'(3X," SYMMETRIC SPECIAL K POINTS")')
    ELSE
       IF (paral%io_parent)&
            WRITE(6,'(3X," NON SYMMETRIC SPECIAL K POINTS")')
    ENDIF
    IF ((tkpts%tkfull).AND.paral%io_parent)&
         WRITE(6,'(3X,A,T41,A25)')&
         ' FULL MONKHORST-PACK MESH','(ONLY INVERSION SYMMETRY)'
    ! ==--------------------------------------------------------------==
    ! == SUBROUTINE K290                                              ==
    ! ==--------------------------------------------------------------==
    CALL k290(iout,ions1%nat,nkpoint,ions1%nsp,kpts_com%nk1,kpts_com%nk2,kpts_com%nk3,istriz,&
         parm%a1,parm%a2,parm%a3,parm%alat,strain,xkapa,rx,tvect,&
         ty,isc,f0,ntvect,wvk0,wvkl,lwght,lrot,&
         nhash,includ,list,rlist,symmr%deltasym)
    ! ==--------------------------------------------------------------==
    isum=0
    nkpt%nkpts=0
    DO ikind=1,nkpoint
       IF (lwght(ikind).NE.0) THEN
          nkpt%nkpts=ikind
          isum = isum + lwght(nkpt%nkpts)
       ELSE
          GOTO 1000
       ENDIF
    ENDDO
1000 CONTINUE
    sum=REAL(isum,kind=real_8)
    ALLOCATE(wk(nkpt%nkpts),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(rk(3,nkpt%nkpts),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    DO ikind = 1 , nkpt%nkpts
       weight = REAL(lwght(ikind),kind=real_8)/sum
       wk(ikind)=weight
       ! wvkl are in cartesian coordinates.
       DO l = 1 , 3
          rk(l,ikind)=REAL(wvkl(l,ikind),kind=real_8)
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent)&
         WRITE(6,'(A,T15,A51)') ' -->        ',&
         '(SEE THE FILE KPTS_GENERATION FOR MORE INFORMATION)'
    IF (paral%io_parent)&
         WRITE(6,'(1X,A,T58,I8)')&
         'NUMBER OF SPECIAL K POINTS (IN CARTESIAN COORDINATES):',&
         nkpt%nkpts
    IF (paral%io_parent)&
         WRITE(6,'(T5,"NKP",T14,"KX",T29,"KY",T44,"KZ",T58,"WEIGHT")')
    DO ikind = 1 , nkpt%nkpts
       IF (paral%io_parent)&
            WRITE(6,'(2X,I4,3F15.6,F14.6)')&
            ikind,(rk(l,ikind),l=1,3),wk(ikind)
    ENDDO
    IF (paral%io_parent)&
         WRITE(6,'(1X,64("*"),/)')
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent)&
         CALL fileclose(iout)
    DEALLOCATE(xkapa,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rx,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(tvect,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(ty,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(isc,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(f0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(wvkl,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(lwght,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(lrot,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(includ,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(list,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rlist,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! CALL TIHALT('    BZMESH',ISUB)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE bzmesh
  ! ==================================================================
  SUBROUTINE setmask(gcutw)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: gcutw

    INTEGER                                  :: ig, ikpt
    REAL(real_8)                             :: xsum

! Variables
! ==--------------------------------------------------------------==
! MASKGW(1,*,*) mask for  g
! MASKGW(2,*,*) mask for -g

    DO ikpt=1,nkpt%nblkp
       IF (tkpts%tkblock) CALL rkpt_swap(hgkp,1,ikpt,'HGKP HGKM')
       CALL calc_maskgw(ikpt)
       IF (tkpts%tkblock) CALL wkpt_swap(maskgw,1,ikpt,'MASKGW')
    ENDDO
    xsum=0._real_8
    DO ig=1,ncpw%ngw
       IF (hg(ig).LT.gcutw) xsum=xsum+1._real_8
    ENDDO
    CALL mp_sum(xsum,parai%allgrp)
    ngw_gamma=2*NINT(xsum)-1
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE setmask
  ! ==================================================================
  SUBROUTINE calc_maskgw(ikpt)
    ! ==--------------------------------------------------------------==
    ! == CALCULATION OF MASKGW FOR A SET OF NKPNT K POINTS            ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ikpt

    INTEGER                                  :: ig, igg, ik, ikk, nglow
    REAL(real_8)                             :: xsum

! Variables
! ==--------------------------------------------------------------==
! MASKGW(1,*,*) mask for  g
! MASKGW(2,*,*) mask for -g

    nglow=ncpw%ngw-maskl+1
    CALL zeroing(maskgw)!,2*maskl*nkpt%nkpnt)
    DO ik=1,nkpbl(ikpt)
       ikk=kpbeg(ikpt)+ik
       DO ig=nglow,ncpw%ngw
          igg=ig-nglow+1
          IF (hgkp(ig,ik).LT.gvec_com%gcutw) maskgw(1,igg,ik)=1.0_real_8
          IF (hgkm(ig,ik).LT.gvec_com%gcutw) maskgw(2,igg,ik)=1.0_real_8
       ENDDO
       ngwkps(ikk)=2*(ncpw%ngw-maskl)
       DO igg=1,maskl
          IF (maskgw(1,igg,ik).GT.0) ngwkps(ikk)=ngwkps(ikk)+1
          IF (maskgw(2,igg,ik).GT.0) ngwkps(ikk)=ngwkps(ikk)+1
       ENDDO
       xsum=ngwkps(ikk)
       CALL mp_sum(xsum,parai%allgrp)
       ngwkps(ikk)=NINT(xsum)-1
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE calc_maskgw
  ! ==================================================================
  SUBROUTINE calc_hgk(ikpt)
    ! ==--------------------------------------------------------------==
    ! == CALCULATION OF HGKM AND HGKP FOR A SET OF NKPNT K POINTS     ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ikpt

    INTEGER                                  :: i, ig, ik, ikk, ir, j, k, &
                                                nh1, nh2, nh3
    REAL(real_8)                             :: g2m, g2p, t, tm, tp

    nh1=spar%nr1s/2+1
    nh2=spar%nr2s/2+1
    nh3=spar%nr3s/2+1
    DO ik=1,nkpbl(ikpt)
       ikk=kpbeg(ikpt)+ik
       DO ig=1,ncpw%nhg
          i=inyh(1,ig)-nh1
          j=inyh(2,ig)-nh2
          k=inyh(3,ig)-nh3
          g2p=0._real_8
          g2m=0._real_8
          DO ir=1,3
             t=REAL(i,kind=real_8)*gvec_com%b1(ir)+REAL(j,kind=real_8)*gvec_com%b2(ir)+REAL(k,kind=real_8)*gvec_com%b3(ir)
             tp=rk(ir,ikk)+t
             tm=rk(ir,ikk)-t
             g2p=g2p+tp*tp
             g2m=g2m+tm*tm
          ENDDO
          hgkp(ig,ik)=g2p
          hgkm(ig,ik)=g2m
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE calc_hgk
  ! ==================================================================

END MODULE rkpnt_utils
