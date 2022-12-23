MODULE fftnew_utils
  USE cell,                            ONLY: cell_com
  USE cnst,                            ONLY: pi
  USE cp_cuda_types,                   ONLY: cp_cuda_env
  USE cp_cufft_types,                  ONLY: cp_cufft,&
                                             cp_cufft_device_get_ptrs
  USE cppt,                            ONLY: hg,&
                                             indz,&
                                             indzs,&
                                             inyh,&
                                             nzh,&
                                             nzhs
  USE cuda_types,                      ONLY: cuda_memory_t
  USE cuda_utils,                      ONLY: cuda_memcpy_host_to_device
  USE error_handling,                  ONLY: stopgm
  USE fft,                             ONLY: &
       fftpool, fftpoolsize, fpoolv, inzf, inzfp, inzh, inzhp, inzs, inzsp, &
       jgw, jgws, jhg, jhgs, kr2max, kr2min, kr3max, kr3min, lfrm, llr1, &
       lmsq, lmsqmax, lnzf, lnzs, lr1, lr1m, lr1s, lr2, lr2s, lr3, lr3s, &
       lrxpl, lrxpool, lsrm, mfrays, mg, msp, msqf, msqfpool, msqs, msqspool, &
       msrays, mz, ngrm, nhrm, nr1m, nzff, nzffp, nzfs, nzfsp, qr1, qr1s, &
       qr2, qr2max, qr2min, qr2s, qr3, qr3max, qr3min, qr3s, sp5, sp8, sp9, &
       spm
  USE fft_maxfft,                      ONLY: maxfft,&
                                             maxfftn
  USE fftchk_utils,                    ONLY: fftchk
  USE kinds,                           ONLY: real_8
  USE loadpa_utils,                    ONLY: leadim
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum
  USE parac,                           ONLY: parai,&
                                             paral
  USE system,                          ONLY: fpar,&
                                             ncpw,&
                                             parap,&
                                             parm,&
                                             spar
  USE utils,                           ONLY: icopy
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: setfftn
  !public :: rmfftnset
  PUBLIC :: addfftnset
  !public :: setrays

CONTAINS

  ! ==================================================================
  SUBROUTINE setfftn(ipool)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ipool

    INTEGER                                  :: i_device, ip, ipx
    TYPE(cuda_memory_t), POINTER             :: inzs_d, lrxpl_d, msqf_d, &
                                                msqs_d, nzfs_d, sp5_d, sp8_d, &
                                                sp9_d

    IF (ipool.EQ.0) THEN
       msrays   = parai%ngrays
       mfrays   = parai%nhrays
       llr1     = fpar%nnr1
       qr1s     = fpar%kr1s
       qr2s     = fpar%kr2s
       qr3s     = fpar%kr3s
       qr1      = fpar%kr1
       lr1s     = spar%nr1s
       lr2s     = spar%nr2s
       lr3s     = spar%nr3s
       lr1      = parm%nr1
       qr2max   = kr2max
       qr2min   = kr2min
       qr3max   = kr3max
       qr3min   = kr3min
       lsrm     = ngrm
       lfrm     = nhrm
       lr1m     = nr1m
       lmsq     = nhrm
       maxfftn  = maxfft
       jgw      = ncpw%ngw
       jgws     = spar%ngws
       jhg      = ncpw%nhg
       jhgs     = spar%nhgs
       nzff => nzh
       inzf => indz
       nzfs => nzhs
       inzs => indzs
       inzh => inyh

       DO ip=0,parai%nproc-1
          lrxpl(ip,1)    = parap%nrxpl(ip,1)
          lrxpl(ip,2)    = parap%nrxpl(ip,2)
          sp5(ip)        = parap%sparm(5,ip)
          sp8(ip)        = parap%sparm(8,ip)
          sp9(ip)        = parap%sparm(9,ip)
          ipx=lmsq*ip
          IF (nhrm.GT.0) THEN
             CALL icopy(nhrm,msp(1,1,ip+1),1,msqf(ipx+1),1)
             CALL icopy(nhrm,msp(1,2,ip+1),1,msqs(ipx+1),1)
          ENDIF
       ENDDO
    ELSEIF (ipool.GT.0 .AND. ipool.LE.fftpool) THEN
       ! LOAD FROM POOL
       msrays   = fpoolv( 1,ipool)
       mfrays   = fpoolv( 2,ipool)
       llr1     = fpoolv( 3,ipool)
       qr1s     = fpoolv( 4,ipool)
       qr2s     = fpoolv( 5,ipool)
       qr3s     = fpoolv( 6,ipool)
       qr1      = fpoolv( 7,ipool)
       qr2      = fpoolv( 8,ipool)
       qr3      = fpoolv( 9,ipool)
       lr1s     = fpoolv(10,ipool)
       lr2s     = fpoolv(11,ipool)
       lr3s     = fpoolv(12,ipool)
       lr1      = fpoolv(13,ipool)
       lr2      = fpoolv(14,ipool)
       lr3      = fpoolv(15,ipool)
       qr2max   = fpoolv(16,ipool)
       qr2min   = fpoolv(17,ipool)
       qr3max   = fpoolv(18,ipool)
       qr3min   = fpoolv(19,ipool)
       lsrm     = fpoolv(20,ipool)
       lfrm     = fpoolv(21,ipool)
       lr1m     = fpoolv(22,ipool)
       lmsq     = fpoolv(23,ipool)
       maxfftn  = fpoolv(24,ipool)
       jgw      = fpoolv(25,ipool)
       jgws     = fpoolv(26,ipool)
       jhg      = fpoolv(27,ipool)
       jhgs     = fpoolv(28,ipool)
       nzff => nzffp(:, ipool)
       inzf => inzfp(:, ipool)
       nzfs => nzfsp(:, ipool)
       inzs => inzsp(:, ipool)
       inzh => inzhp(:, :, ipool)

       DO ip=0,parai%nproc-1
          lrxpl(ip,1)    = lrxpool(ip,1,ipool)
          lrxpl(ip,2)    = lrxpool(ip,2,ipool)
          sp5(ip)        = spm(5,ip,ipool)
          sp8(ip)        = spm(8,ip,ipool)
          sp9(ip)        = spm(9,ip,ipool)
          ipx=lmsq*ip
          IF (lmsqmax.GT.0) THEN
             CALL icopy(lmsq,msqfpool(1,ip+1,ipool),1,msqf(ipx+1),1)
             CALL icopy(lmsq,msqspool(1,ip+1,ipool),1,msqs(ipx+1),1)
          ENDIF
       ENDDO
    ELSE
       CALL stopgm("SETFFTN","FFTPOOL NOT DEFINED",& 
            __LINE__,__FILE__)
    ENDIF

    !vw copy FFT arrays to GPU memory
    IF( cp_cuda_env%use_fft ) THEN
       DO i_device = 1, cp_cuda_env%fft_n_devices_per_task
          CALL cp_cufft_device_get_ptrs ( cp_cufft, i_device, sp5_d=sp5_d, sp8_d=sp8_d, sp9_d=sp9_d, &
               & msqs_d=msqs_d, msqf_d=msqf_d, lrxpl_d=lrxpl_d, nzfs_d=nzfs_d, inzs_d=inzs_d )
          CALL cuda_memcpy_host_to_device ( sp5, sp5_d )
          CALL cuda_memcpy_host_to_device ( sp8, sp8_d )
          CALL cuda_memcpy_host_to_device ( sp9, sp9_d )
          CALL cuda_memcpy_host_to_device ( lrxpl, lrxpl_d )
          CALL cuda_memcpy_host_to_device ( msqf, msqf_d )
          CALL cuda_memcpy_host_to_device ( msqs, msqs_d )
          CALL cuda_memcpy_host_to_device ( nzfs, nzfs_d )
          CALL cuda_memcpy_host_to_device ( inzs, inzs_d )
       ENDDO
    ENDIF

    ! ==--------------------------------------------------------------==
  END SUBROUTINE setfftn
  ! ==================================================================
  SUBROUTINE rmfftnset(ipool)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ipool

    INTEGER                                  :: i, ip, j

    IF (ipool.GT.0 .AND. ipool.LE.fftpool) THEN
       DO ip=ipool+1,fftpool
          i=ip
          j=ip-1
          CALL icopy(28,fpoolv(1,i),1,fpoolv(1,j),1)
          CALL icopy(lnzf,nzffp(1,i),1,nzffp(1,j),1)
          CALL icopy(lnzf,inzfp(1,i),1,inzfp(1,j),1)
          CALL icopy(lnzs,nzfsp(1,i),1,nzfsp(1,j),1)
          CALL icopy(lnzs,inzsp(1,i),1,inzsp(1,j),1)
          CALL icopy(3*lnzf,inzhp(1,1,i),1,inzhp(1,1,j),1)
          CALL icopy(2*SIZE(lrxpool,1),lrxpool(0,1,i),1,lrxpool(0,1,j),1)
          CALL icopy(9*SIZE(lrxpool,1),spm(1,0,i),1,spm(1,0,j),1)
          IF (lmsqmax.GT.0) THEN
             CALL icopy(lmsqmax*parai%nproc,msqfpool(1,1,i),1,&
                  msqfpool(1,1,j),1)
             CALL icopy(lmsqmax*parai%nproc,msqspool(1,1,i),1,&
                  msqspool(1,1,j),1)
          ENDIF
       ENDDO
       CALL zeroing(fpoolv(:,fftpool))!,28)
       CALL zeroing(nzffp(:,fftpool))!,lnzf)
       CALL zeroing(inzfp(:,fftpool))!,lnzf)
       CALL zeroing(nzfsp(:,fftpool))!,lnzf)
       CALL zeroing(inzsp(:,fftpool))!,lnzf)
       CALL zeroing(inzhp(:,:,fftpool))!,3*lnzf)
       CALL zeroing(lrxpool(:,:,fftpool))!,2*maxcpu+2)
       CALL zeroing(spm(:,:,fftpool))!,9*maxcpu+9)
       IF (lmsqmax.GT.0) THEN
          CALL zeroing(msqfpool(:,:,fftpool))!,lmsqmax*parai%nproc)
          CALL zeroing(msqspool(:,:,fftpool))!,lmsqmax*parai%nproc)
       ENDIF
       fftpool=fftpool-1
    ELSE
       CALL stopgm("RMFFTNSET","FFTPOOL NOT DEFINED",& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
  END SUBROUTINE rmfftnset
  ! ==================================================================
  SUBROUTINE addfftnset(ecutf,ecuts,ipool)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: ecutf, ecuts
    INTEGER                                  :: ipool

    CHARACTER(*), PARAMETER                  :: procedureN = 'addfftnset'

    INTEGER                                  :: i, ierr, ig, j, k, l, lh1, &
                                                lh2, lh3, nh1, nh2, nh3
    INTEGER, SAVE                            :: icount = 0
    REAL(real_8)                             :: aa1, aa2, aa3, rr, xpaim, &
                                                xplanes, xpnow

! ==--------------------------------------------------------------==

    IF (icount.EQ.0) THEN
       lmsqmax=nhrm
       CALL zeroing(lrxpool)!,2*fftpoolsize*(maxcpu+1))
       CALL zeroing(spm)!,9*fftpoolsize*(maxcpu+1))
       l=(lmsqmax*parai%nproc)
       ALLOCATE(msqf(l),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(msqs(l),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       l=(lmsqmax*parai%nproc*fftpoolsize)
       ALLOCATE(msqfpool(lmsqmax,parai%nproc,fftpoolsize),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(msqspool(lmsqmax,parai%nproc,fftpoolsize),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       l=(ncpw%nhg*fftpoolsize)
       ALLOCATE(nzffp(lnzf,l/lnzf),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)

       IF (lnzs > 0) THEN
          ALLOCATE(nzfsp(lnzs,l/lnzs),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ELSE
          ALLOCATE(nzfsp(1,l),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF

       ALLOCATE(inzfp(lnzf,l/lnzf),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)

       IF (lnzs > 0) THEN
          ALLOCATE(inzsp(lnzs,l/lnzs),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ELSE
          ALLOCATE(inzsp(1,l),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF

       l=(3*ncpw%nhg*fftpoolsize)
       ALLOCATE(inzhp(3,lnzf,l/(3*lnzf)),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    icount=1
    IF (ecutf.LT.0._real_8 .OR. ecuts.LT.0._real_8) RETURN
    ! 
    ! IF(TKPNT) CALL STOPGM("ADDFFTNSET","NOT IMPLEMENTED")
    ! 
    fftpool=fftpool+1
    IF (fftpool.GT.fftpoolsize) THEN
       CALL stopgm("ADDFFTNSET","TOO MANY ENTRIES IN POOL",& 
            __LINE__,__FILE__)
    ENDIF
    ipool=fftpool
    ! 
    jhg=ncpw%nhg
    DO ig=2,ncpw%nhg
       IF (parm%tpiba2*hg(ig).GT.ecutf) THEN
          jhg=ig-1
          GOTO 100
       ENDIF
    ENDDO
100 CONTINUE
    rr=REAL(jhg,kind=real_8)
    CALL mp_sum(rr,parai%allgrp)
    jhgs=NINT(rr)
    IF (jhgs.GT.spar%nhgs) CALL stopgm("ADDFFTNSET","JHGS TOO LARGE",& 
         __LINE__,__FILE__)
    jgw=ncpw%nhg
    DO ig=2,ncpw%nhg
       IF (parm%tpiba2*hg(ig).GT.ecuts) THEN
          jgw=ig-1
          GOTO 101
       ENDIF
    ENDDO
101 CONTINUE
    rr=REAL(jgw,kind=real_8)
    CALL mp_sum(rr,parai%allgrp)
    jgws=NINT(rr)
    IF (jgws.GT.spar%ngws) CALL stopgm("ADDFFTNSET","JGWS TOO LARGE",& 
         __LINE__,__FILE__)
    ! 
    aa1=parm%alat
    aa2=parm%alat*cell_com%celldm(2)
    aa3=parm%alat*cell_com%celldm(3)
    lr1s=NINT(aa1/pi*SQRT(ecutf)+0.5_real_8)
    lr2s=NINT(aa2/pi*SQRT(ecutf)+0.5_real_8)
    lr3s=NINT(aa3/pi*SQRT(ecutf)+0.5_real_8)
    ! 
    lr1s=fftchk(lr1s,2)
    lr2s=fftchk(lr2s,2)
    lr3s=fftchk(lr3s,2)
    CALL leadim(lr1s,lr2s,lr3s,qr1s,qr2s,qr3s)
    ! 
    nh1=spar%nr1s/2+1
    nh2=spar%nr2s/2+1
    nh3=spar%nr3s/2+1
    lh1=lr1s/2+1
    lh2=lr2s/2+1
    lh3=lr3s/2+1
    !$omp parallel do private(IG,I,J,K) shared(NH1,NH2,NH3,LH1,LH2,LH3)
    DO ig=1,jhg
       i=inyh(1,ig)-nh1
       j=inyh(2,ig)-nh2
       k=inyh(3,ig)-nh3
       inzhp(1,ig,ipool)=lh1+i
       inzhp(2,ig,ipool)=lh2+j
       inzhp(3,ig,ipool)=lh3+k
    ENDDO
    inzh => inzhp(:,:,ipool)
    ! 
    CALL zeroing(lrxpool(:,:,ipool))!,2*(maxcpu+1))
    xplanes=REAL(lr1s,kind=real_8)
    xpnow=0.0_real_8
    DO i=parai%nproc,1,-1
       xpaim = xpnow + xplanes/parai%nproc
       lrxpool(i-1,1,ipool)=NINT(xpnow)+1
       lrxpool(i-1,2,ipool)=NINT(xpaim)
       IF (NINT(xpaim).GT.lr1s) lrxpool(i-1,2,ipool)=lr1s
       IF (i.EQ.1) lrxpool(i-1,2,ipool)=lr1s
       xpnow = xpaim
    ENDDO
    lr1=lrxpool(parai%mepos,2,ipool)-lrxpool(parai%mepos,1,ipool)+1
    CALL leadim(lr1,lr2s,lr3s,qr1,qr2s,qr3s)
    lr2=lr2s
    lr3=lr3s
    qr2=qr2s
    qr3=qr3s
    llr1=qr1*qr2*qr3
    ! 
    CALL setrays(ipool)
    ! 
    maxfftn = maxfft
    ! 
    fpoolv( 1,ipool) = msrays
    fpoolv( 2,ipool) = mfrays
    fpoolv( 3,ipool) = llr1
    fpoolv( 4,ipool) = qr1s
    fpoolv( 5,ipool) = qr2s
    fpoolv( 6,ipool) = qr3s
    fpoolv( 7,ipool) = qr1
    fpoolv( 8,ipool) = qr2
    fpoolv( 9,ipool) = qr3
    fpoolv(10,ipool) = lr1s
    fpoolv(11,ipool) = lr2s
    fpoolv(12,ipool) = lr3s
    fpoolv(13,ipool) = lr1
    fpoolv(14,ipool) = lr2
    fpoolv(15,ipool) = lr3
    fpoolv(16,ipool) = qr2max
    fpoolv(17,ipool) = qr2min
    fpoolv(18,ipool) = qr3max
    fpoolv(19,ipool) = qr3min
    fpoolv(20,ipool) = lsrm
    fpoolv(21,ipool) = lfrm
    fpoolv(22,ipool) = lr1m
    fpoolv(23,ipool) = lmsq
    fpoolv(24,ipool) = maxfftn
    fpoolv(25,ipool) = jgw
    fpoolv(26,ipool) = jgws
    fpoolv(27,ipool) = jhg
    fpoolv(28,ipool) = jhgs
    ! 
    IF (paral%io_parent) THEN
       WRITE(6,*)
       WRITE(6,'(A,T50,A,I4)') ' ADD NEW FFT SET ',' SET NUMBER ',&
            ipooL
       WRITE(6,'(A,T51,3I5)') ' REAL SPACE GRID ',lr1s,lr2s,lr3S
       WRITE(6,'(A,T20,A,F6.0,T44,A,I10)') ' SPARSE FFT SETUP: ',&
            'CUTOFF [Ry]:',ecuts,'PLANE WAVES:',jgwS
       WRITE(6,'(A,T20,A,F6.0,T44,A,I10)') ' FULL FFT SETUP  : ',&
            'CUTOFF [Ry]:',ecutf,'PLANE WAVES:',jhgS
       WRITE(6,*)
    ENDIF
    ! ==--------------------------------------------------------------==
  END SUBROUTINE addfftnset
  ! ==================================================================
  SUBROUTINE setrays(ipool)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ipool

    CHARACTER(*), PARAMETER                  :: procedureN = 'setrays'

    INTEGER :: i, ierr, ig, ij, img, iny1, iny2, iny3, ip, ipro, ixf, j, &
      jgwl, jhgl, jj, jmg, msglen, mxrp, nh1, nh2, nh3, ny1, ny2, ny3, qr1m
    INTEGER, ALLOCATABLE                     :: mq(:), my(:)

! Variables
! ==--------------------------------------------------------------==
! GATHER ARRAY FOR FFT ALONG X
! SPARSITY FOR FFT ALONG Y

    ALLOCATE(mg(qr2s,qr3s),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(mz((2*qr3s)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(my((2*qr2s)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(mg)!,qr2s*qr3s)
    CALL zeroing(mz)!,2*qr3s)
    CALL zeroing(my)!,2*qr2s)
    nh1=lr1s/2+1
    nh2=lr2s/2+1
    nh3=lr3s/2+1
    DO ig=1,jgw
       ny2=inzh(2,ig)
       ny3=inzh(3,ig)
       iny2=-ny2+2*nh2
       iny3=-ny3+2*nh3
       mg(ny2,ny3)=mg(ny2,ny3)+1
       mg(iny2,iny3)=mg(iny2,iny3)+1
       my(ny2)=my(ny2)+1
       my(iny2)=my(iny2)+1
       mz(ny3)=mz(ny3)+1
       mz(iny3)=mz(iny3)+1
    ENDDO
    CALL mp_sum(my,my(lr2s+1:),lr2s,parai%allgrp)
    CALL icopy(lr2s,my(lr2s+1),1,my,1)
    CALL mp_sum(mz,mz(lr3s+1:),lr3s,parai%allgrp)
    CALL icopy(lr3s,mz(lr3s+1),1,mz,1)
    qr2min=1
    DO i=1,qr2s
       IF (my(i).NE.0) THEN
          qr2min=i
          GOTO 50
       ENDIF
    ENDDO
50  CONTINUE
    qr2max=qr2s
    DO i=qr2s,1,-1
       IF (my(i).NE.0) THEN
          qr2max=i
          GOTO 51
       ENDIF
    ENDDO
51  CONTINUE
    qr3min=1
    DO i=1,qr3
       IF (mz(i).NE.0) THEN
          qr3min=i
          GOTO 52
       ENDIF
    ENDDO
52  CONTINUE
    qr3max=1
    DO i=qr3,1,-1
       IF (mz(i).NE.0) THEN
          qr3max=i
          GOTO 53
       ENDIF
    ENDDO
53  CONTINUE
    ! ==--------------------------------------------------------------==
    img=0
    DO j=1,qr3s
       DO i=1,qr2s
          IF (mg(i,j).NE.0) THEN
             img=img+1
             mg(i,j)=img
          ENDIF
       ENDDO
    ENDDO
    msrays=img
    DO ig=jgw+1,jhg
       ny2=inzh(2,ig)
       ny3=inzh(3,ig)
       iny2=-ny2+2*nh2
       iny3=-ny3+2*nh3
       jmg=mg(ny2,ny3)
       IF (jmg.EQ.0) mg(ny2,ny3)=-1
       jmg=mg(iny2,iny3)
       IF (jmg.EQ.0) mg(iny2,iny3)=-1
    ENDDO
    DO j=1,qr3s
       DO i=1,qr2s
          IF (mg(i,j).LT.0) THEN
             img=img+1
             mg(i,j)=img
          ENDIF
       ENDDO
    ENDDO
    mfrays=img
    ! 
    jhgl=0
    jgwl=0
    CALL zeroing(spm(:,:,ipool))!,9*maxcpu+9)
    spm(1,parai%mepos,ipool)=jhg
    spm(2,parai%mepos,ipool)=jhgl
    spm(3,parai%mepos,ipool)=jgw
    spm(4,parai%mepos,ipool)=jgwl
    spm(5,parai%mepos,ipool)=lr1
    spm(6,parai%mepos,ipool)=lr2
    spm(7,parai%mepos,ipool)=lr3
    spm(8,parai%mepos,ipool)=mfrays
    spm(9,parai%mepos,ipool)=msrays
    ! 
    nzff => nzffp(:,ipool)
    inzf => inzfp(:,ipool)
    nzfs => nzfsp(:,ipool)
    inzs => inzsp(:,ipool)
    DO i=0,parai%nproc-1
       ipro=i
       CALL mp_bcast(spm(:,i,ipool),9,ipro,parai%allgrp)
    ENDDO
    ! MAXIMUM OF LR1, MSRAYS AND MFRAYS FOR MP_INDEX
    lr1m = 0
    lfrm = 0
    lsrm = 0
    DO i=0,parai%nproc-1
       lr1m = MAX(lr1m,spm(5,i,ipool))
       lfrm = MAX(lfrm,spm(8,i,ipool))
       lsrm = MAX(lsrm,spm(9,i,ipool))
    ENDDO
    qr1m=MAX(lr1m+MOD(lr1m+1,2),qr1)
    lmsq=MAX(lfrm,lsrm)
    ! SCATTER ARRAY FOR FFT ALONG X
    ALLOCATE(mq(lmsq*2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(mq)!,2*lmsq)
    DO i=1,qr2s
       DO j=1,qr3s
          ij=mg(i,j)
          IF (ij.GT.0) THEN
             mq(ij)=i
             mq(lmsq+ij)=j
          ENDIF
       ENDDO
    ENDDO
    ! CONCATENATE GATHER/SCATTER ARRAYS
    msglen = lfrm * 8/2
    CALL my_concat(mq(1),msqs,msglen,parai%allgrp)
    CALL my_concat(mq(lmsq+1),msqf,msglen,parai%allgrp)
    ! TRANSLATE I,J TO A SINGLE G/S INDEX
    DO ip=0,parai%nproc-1
       mxrp=parap%sparm(8,ip)
       DO ixf=1,mxrp
          jj=ixf+ip*lmsq
          i=msqs(jj)
          j=msqf(jj)
          msqfpool(ixf,ip+1,ipool)=i+(j-1)*qr2s
          IF (ixf.LE.parap%sparm(9,ip))&
               msqspool(ixf,ip+1,ipool)=i+(j-qr3min)*qr2s
       ENDDO
    ENDDO
    ! REDEFINE NZH AND INDZ FOR COMPRESSED STORAGE
    !$omp parallel do private(IG,NY1,NY2,NY3,INY1,INY2,INY3)
    DO ig=1,jhg
       ny1=inzh(1,ig)
       ny2=inzh(2,ig)
       ny3=inzh(3,ig)
       iny1=-ny1+2*nh1
       iny2=-ny2+2*nh2
       iny3=-ny3+2*nh3
       nzff(ig)=ny1 + (mg(ny2,ny3)-1)*qr1s
       inzf(ig)=iny1 + (mg(iny2,iny3)-1)*qr1s
    ENDDO
    !$omp parallel do private(IG)
    DO ig=1,jgw
       nzfs(ig)=nzff(ig)
       inzs(ig)=inzf(ig)
    ENDDO
    DEALLOCATE(mg,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(my,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(mz,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(mq,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE setrays
  ! ==================================================================

END MODULE fftnew_utils
