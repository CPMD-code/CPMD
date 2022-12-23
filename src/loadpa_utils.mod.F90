#include "cpmd_global.h"

MODULE loadpa_utils
  USE cppt,                            ONLY: hg,&
                                             inyh
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE geq0mod,                         ONLY: geq0
  USE gvec,                            ONLY: epsg,&
                                             epsgx,&
                                             gvec_com
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE isos,                            ONLY: isos1,&
                                             isos3
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE mm_dim_utils,                    ONLY: mm_dim
  USE mm_dimmod,                       ONLY: mm_go_mm,&
                                             mm_revert
  USE mp_interface,                    ONLY: mp_sum,&
                                             mp_sync
  USE parac,                           ONLY: parai,&
                                             paral
  USE prmem_utils,                     ONLY: prmem
  USE sort_utils,                      ONLY: sort2
  USE sphe,                            ONLY: gcutwmax
  USE system,                          ONLY: &
       fpar, iatpe, iatpt, ipept, mapgp, natpe, ncpw, nkpt, norbpe, parap, &
       parm, spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  !$ use omp_lib, only: omp_get_max_threads, omp_get_thread_num

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: loadpa
  PUBLIC :: leadim

CONTAINS

  ! ==================================================================
  SUBROUTINE loadpa
    ! ==--------------------------------------------------------------==
    ! Variables
    CHARACTER(*), PARAMETER                  :: procedureN = 'loadpa'

    INTEGER :: i, i0, ia, iat, icpu, ierr, ig, ihrays, ii, img, in1, in2, &
      in3, iorb, ip, ipp, ir, is, isub, isub2, isub3, isub4, ixrays, izpl, j, &
      j1, j2, jmax, jmin, k, kmax, kmin, mspace, nh1, nh2, nh3, nthreads
    INTEGER, ALLOCATABLE                     :: ihray(:,:), ixray(:,:), &
                                                mgpa(:,:)
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: thread_buff
    LOGICAL                                  :: oldstatus
    REAL(real_8)                             :: g2, sign, t, xpaim, xplanes, &
                                                xpnow, xsaim, xsnow, xstates, &
                                                zpaim, zplanes, zpnow

! ==--------------------------------------------------------------==
! ==  DISTRIBUTION OF PARALLEL WORK                               ==
! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    nthreads = 1
    !$    nthreads = omp_get_max_threads( )
    ALLOCATE(thread_buff(2*nthreads),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    IF (paral%io_parent)WRITE(6,'(/," ",16("PARA"))')
    mspace = (fpar%kr2s*fpar%kr3s)/2 + 1
    ALLOCATE(ixray(fpar%kr2s,fpar%kr3s),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(ihray(fpar%kr2s,fpar%kr3s),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(mgpa(fpar%kr2s,fpar%kr3s),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(iatpe(ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL mp_sync(parai%allgrp)
    ! ==--------------------------------------------------------------==
    ! DISTRIBUTE ATOMS
    ! ==--------------------------------------------------------------==
    xstates=REAL(ions1%nat,kind=real_8)
    xsnow=0.0_real_8
    DO i=parai%nproc,1,-1
       xsaim = xsnow + xstates/parai%nproc
       ipept(1,i-1)=NINT(xsnow)+1
       ipept(2,i-1)=NINT(xsaim)
       IF (NINT(xsaim).GT.ions1%nat) THEN
          ipept(2,i-1)=ions1%nat
       ENDIF
       IF (i.EQ.1) THEN
          ipept(2,i-1)=ions1%nat
       ENDIF
       xsnow = xsaim
    ENDDO

    CALL mm_dim(mm_go_mm,oldstatus)
    ALLOCATE(iatpt(2,ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    iat=0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          iat=iat+1
          iatpt(1,iat)=ia
          iatpt(2,iat)=is
       ENDDO
    ENDDO
    CALL mm_dim(mm_revert,oldstatus)

    DO i=0,parai%nproc-1
       DO j=ipept(1,i),ipept(2,i)
          iatpe(j)=i
       ENDDO
    ENDDO
    natpe=ipept(2,parai%mepos)-ipept(1,parai%mepos)+1
    ! ==--------------------------------------------------------------==
    ! DISTRIBUTE ORBITALS
    ! ==--------------------------------------------------------------==
    xstates=REAL(crge%n,kind=real_8)
    xsnow=0.0_real_8
    DO i=parai%nproc,1,-1
       xsaim = xsnow + xstates/parai%nproc
       parap%nst12(i-1,1)=NINT(xsnow)+1
       parap%nst12(i-1,2)=NINT(xsaim)
       IF (NINT(xsaim).GT.crge%n) THEN
          parap%nst12(i-1,2)=crge%n
       ENDIF
       IF (i.EQ.1) THEN
          parap%nst12(i-1,2)=crge%n
       ENDIF
       xsnow = xsaim
    ENDDO
    norbpe=parap%nst12(parai%mepos,2)-parap%nst12(parai%mepos,1)+1
    ! ==--------------------------------------------------------------==
    ! DISTRIBUTE REAL SPACE YZ-PLANES
    ! ==--------------------------------------------------------------==
    CALL zeroing(parap%nrxpl)!,2*(maxcpu+1))
    xplanes=REAL(spar%nr1s,kind=real_8)
    xpnow=0.0_real_8
    DO i=parai%nproc,1,-1
       xpaim = xpnow + xplanes/parai%nproc
       parap%nrxpl(i-1,1)=NINT(xpnow)+1
       parap%nrxpl(i-1,2)=NINT(xpaim)
       IF (NINT(xpaim).GT.spar%nr1s) THEN
          parap%nrxpl(i-1,2)=spar%nr1s
       ENDIF
       IF (i.EQ.1) THEN
          parap%nrxpl(i-1,2)=spar%nr1s
       ENDIF
       xpnow = xpaim
    ENDDO
    CALL zeroing(parap%nrzpl)!,2*(maxcpu+1))
    IF (isos1%tclust.AND.isos3%ps_type.EQ.1) THEN
       ! DISTRIBUTE REAL SPACE XY-PLANES
       zplanes=REAL(2*spar%nr3s,kind=real_8)
       zpnow=0.0_real_8
       DO i=parai%nproc,1,-1
          zpaim = zpnow + zplanes/parai%nproc
          parap%nrzpl(i-1,1)=NINT(zpnow)+1
          parap%nrzpl(i-1,2)=NINT(zpaim)
          IF (NINT(zpaim).GT.2*spar%nr3s) THEN
             parap%nrzpl(i-1,2)=2*spar%nr3s
          ENDIF
          IF (i.EQ.1) THEN
             parap%nrzpl(i-1,2)=2*spar%nr3s
          ENDIF
          zpnow = zpaim
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    ! DISTRIBUTE G-VECTORS
    ! ==--------------------------------------------------------------==
    CALL tiset(procedureN//'_c',isub4)
    CALL xfft(ixray,gcutwmax)
    CALL xfft(ihray,gvec_com%gcut)
    ixrays = 0
    !$omp parallel do __COLLAPSE2 &
    !$omp             private(I,J) &
    !$omp             shared(fpar) &
    !$omp             reduction(+:IXRAYS)
    DO j=1,fpar%kr3s
       DO i=1,fpar%kr2s
          IF (ixray(i,j).NE.0) THEN
             ixrays=ixrays+1
          ENDIF
       ENDDO
    ENDDO
    !$omp end parallel do
    ipp=0
    DO ii=1,ixrays
       ipp=ipp+1
       ip=MOD(ipp,2*parai%nproc)
       IF (ip.EQ.0) THEN
          ip=2*parai%nproc
       ENDIF
       IF (ip.GT.parai%nproc) THEN
          ip=2*parai%nproc+1-ip
       ENDIF
       CALL iraymax(fpar%kr2s,fpar%kr3s,ixray,i,j,thread_buff,nthreads)
       ixray(i,j)=-ip
    ENDDO
    CALL tihalt(procedureN//'_c',isub4)
    ! ==--------------------------------------------------------------==
    ! MARK ASSOCIATED RAYS   
    ! ==--------------------------------------------------------------==
    CALL tiset(procedureN//'_b',isub3)
    CALL cxfft(ixray,gcutwmax)
    ihrays=0
    !$omp parallel do __COLLAPSE2 &
    !$omp             private(I,J) &
    !$omp             shared(fpar) &
    !$omp             reduction(+:IHRAYS)
    DO j=1,fpar%kr3s
       DO i=1,fpar%kr2s
          IF (ixray(i,j).NE.0) THEN
             ihray(i,j)=ixray(i,j)
          ENDIF
          IF (ihray(i,j).GT.0) THEN
             ihrays=ihrays+1
          ENDIF
       ENDDO
    ENDDO
    !$omp end parallel do
    DO ii=1,ihrays
       ipp=ipp+1
       ip=MOD(ipp,2*parai%nproc)
       IF (ip.EQ.0) THEN
          ip=2*parai%nproc
       ENDIF
       IF (ip.GT.parai%nproc) THEN
          ip=2*parai%nproc+1-ip
       ENDIF
       CALL iraymax(fpar%kr2s,fpar%kr3s,ihray,i,j,thread_buff,nthreads)
       ihray(i,j)=-ip
    ENDDO
    ! MARK ASSOCIATED RAYS   
    CALL cxfft(ihray,gvec_com%gcut)
    CALL tihalt(procedureN//'_b',isub3)
    ! ==--------------------------------------------------------------==
    ! ALLOCATE THE GLOBAL DATA ARRAYS
    ! ==--------------------------------------------------------------==
    ncpw%nhg=INT(REAL(spar%nhgs,kind=real_8)/REAL(parai%nproc,kind=real_8)*1.2_real_8)
    ncpw%nhg=MIN(ncpw%nhg,spar%nhgs)
    ncpw%nhg=MAX(ncpw%nhg,100)
    ncpw%ngw=INT(REAL(spar%ngws,kind=real_8)/REAL(parai%nproc,kind=real_8)*1.2_real_8)
    ncpw%ngw=MIN(ncpw%ngw,spar%ngws)
    ncpw%ngw=MAX(ncpw%ngw,100)
    ALLOCATE(hg(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(inyh(3,ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(mapgp(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(hg)!,nhg)
    CALL zeroing(inyh)!,3*nhg)
    CALL zeroing(mapgp)!,nhg)
    ! ==--------------------------------------------------------------==
    ! SYSTEM PARAMETER PER CPU
    ! ==--------------------------------------------------------------==
    nh1=parm%nr1/2+1
    nh2=parm%nr2/2+1
    nh3=parm%nr3/2+1
    ! ==--------------------------------------------------------------==
    ! TO GET UNIQUE ORDERING OF G-VECTORS, BREAK THE SYMMETRY
    ! ==--------------------------------------------------------------==
    CALL tiset(procedureN//'_a',isub2)
    ! EPSG=EPSGX * ACOS(-1._real_8) *real(NHGS,kind=real_8)
    ! EPSG=GCUT*EPSGX
    epsg=SQRT(REAL(spar%nhgs-1,kind=real_8))*gvec_com%gcut*epsgx
    ig=0
    ncpw%ngw=0
    ncpw%nhg=0
    DO i=0,parm%nr1-1
       jmin=-parm%nr2+1
       jmax=parm%nr2-1
       IF (i.EQ.0) THEN
          jmin=0
       ENDIF
       DO j=jmin,jmax
          kmin=-parm%nr3+1
          kmax=parm%nr3-1
          IF (i.EQ.0.AND.j.EQ.0) THEN
             kmin=0
          ENDIF
          DO k=kmin,kmax
             g2=0._real_8
             DO ir=1,3
                t=REAL(i,kind=real_8)*gvec_com%b1(ir)+REAL(j,kind=real_8)*gvec_com%b2(ir)+REAL(k,kind=real_8)*gvec_com%b3(ir)
                g2=g2+t*t
             ENDDO
             IF (g2.LT.gvec_com%gcut) THEN
                ig=ig+1
                in1=nh1+i
                in2=nh2+j
                in3=nh3+k
                IF (g2.LT.gcutwmax) THEN
                   icpu=ixray(in2,in3)
                   IF (-icpu.EQ.parai%mepos+1) THEN
                      ncpw%ngw=ncpw%ngw+1
                   ENDIF
                   sign=-1._real_8
                ELSE
                   sign=1._real_8
                ENDIF
                icpu=ihray(in2,in3)
                IF (-icpu.EQ.parai%mepos+1) THEN
                   ncpw%nhg=ncpw%nhg+1
                   ! HG(NHG)=G2+SQRT(real(IG-1,kind=real_8))*EPSG*SIGN
                   ! G2*EPSGX*SIGN is the epsilon added (>= epsilon(G2))
                   ! SQRT(FLOAT(IG-1)) is to break the symmetry
                   hg(ncpw%nhg)=g2*(1._real_8+SQRT(REAL(ig-1,kind=real_8))*epsgx*sign)
                   inyh(1,ncpw%nhg)=in1
                   inyh(2,ncpw%nhg)=in2
                   inyh(3,ncpw%nhg)=in3
                ENDIF
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    CALL zeroing(mgpa)!,kr2s*kr3s)
    parai%nhrays=0
    parai%ngrays=0
    DO j1=1,fpar%kr3s
       DO j2=1,fpar%kr2s
          IF (ihray(j2,j1).EQ.-(parai%mepos+1)) parai%nhrays=parai%nhrays+1
          IF (ixray(j2,j1).EQ.-(parai%mepos+1)) THEN
             parai%ngrays=parai%ngrays+1
             mgpa(j2,j1)=parai%ngrays
          ENDIF
       ENDDO
    ENDDO
    img=parai%ngrays
    DO j1=1,fpar%kr3s
       DO j2=1,fpar%kr2s
          IF (ihray(j2,j1).EQ.-(parai%mepos+1).AND.mgpa(j2,j1).EQ.0) THEN
             img=img+1
             mgpa(j2,j1)=img
          ENDIF
       ENDDO
    ENDDO
    ncpw%nhgl=spar%nhgls
    ncpw%ngwl=spar%ngwls
    parm%nr1 =parap%nrxpl(parai%mepos,2)-parap%nrxpl(parai%mepos,1)+1
    parm%nr2 =spar%nr2s
    parm%nr3 =spar%nr3s
    parap%sparm(1,parai%mepos)=ncpw%nhg
    parap%sparm(2,parai%mepos)=ncpw%nhgl
    CALL tihalt(procedureN//'_a',isub2)
    ! ==--------------------------------------------------------------==
    ! IF K POINTS, WE NEED TO DOUBLE DIMENSIONS.
    ! ==--------------------------------------------------------------==
    IF (tkpts%tkpnt) THEN
       nkpt%ngwk=ncpw%ngw*2
       nkpt%nhgk=ncpw%nhg*2
    ELSE
       nkpt%ngwk=ncpw%ngw
       nkpt%nhgk=ncpw%nhg
    ENDIF
    parap%sparm(3,parai%mepos)=ncpw%ngw
    parap%sparm(4,parai%mepos)=ncpw%ngwl
    parap%sparm(5,parai%mepos)=parm%nr1
    parap%sparm(6,parai%mepos)=parm%nr2
    parap%sparm(7,parai%mepos)=parm%nr3
    parap%sparm(8,parai%mepos)=parai%nhrays
    parap%sparm(9,parai%mepos)=parai%ngrays
    CALL my_allgather_i(parap%sparm,SIZE(parap%sparm,1),parai%allgrp)

    IF (paral%io_parent) THEN
       WRITE(6,'(A,A)') '  NCPU     NGW',&
            '     NHG  PLANES  GXRAYS  HXRAYS ORBITALS Z-PLANES'
       DO i=0,parai%nproc-1
          iorb=parap%nst12(i,2)-parap%nst12(i,1)+1
          izpl=parap%nrzpl(i,2)-parap%nrzpl(i,1)+1
          WRITE(6,'(I6,7I8)') i,parap%sparm(3,i),parap%sparm(1,i),parap%sparm(5,i),&
               parap%sparm(9,i),parap%sparm(8,i),iorb,izpL
          ! IF(SPARM(3,I).LE.0) CALL stopgm(procedureN,
          ! *            'NGW .LE. 0')
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    ! DEALLOCATE DATA ARRAYS
    ! ==--------------------------------------------------------------==
    DEALLOCATE(mgpa,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(ixray,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(ihray,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! SORTING OF G-VECTORS
    ! ==--------------------------------------------------------------==
    CALL gsort(hg,inyh,ncpw%nhg)
    DO ig=1,ncpw%nhg
       mapgp(ig)=ig
    ENDDO
    CALL gorder
    ! ==--------------------------------------------------------------==
    geq0=.FALSE.
    i0=0
    DO ig=1,ncpw%ngw
       IF (hg(ig).LT.1.e-5_real_8) THEN
          geq0=.TRUE.
          i0=parai%mepos
       ENDIF
    ENDDO
    CALL mp_sum(i0,parai%igeq0,parai%allgrp)
    IF (paral%io_parent) THEN
       WRITE(6,'(13X,A,I5)') '   G=0 COMPONENT ON PROCESSOR : ',parai%igeq0
       WRITE(6,'(" ",16("PARA"),/)')
       CALL prmem(procedureN)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! LEADING DIMENSIONS OF REAL SPACE ARRAYS
    ! ==--------------------------------------------------------------==
    CALL leadim(parm%nr1,parm%nr2,parm%nr3,fpar%kr1,fpar%kr2,fpar%kr3)
    fpar%nnr1=fpar%kr1*fpar%kr2s*fpar%kr3s
    DEALLOCATE(thread_buff,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE loadpa
  ! ==================================================================
  SUBROUTINE iraymax(n1,n2,ir,i1,i2,thread_buff,nthreads)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: n1, n2, ir(n1,n2), i1, i2, &
                                                nthreads, &
                                                thread_buff(2,0:nthreads-1)

    INTEGER                                  :: im, ithread, j1, j2

! Input
! Output
! ==--------------------------------------------------------------==

    thread_buff(:,:) = 0
    !$omp parallel &
    !$omp          private(J1,J2,IM,ithread) &
    !$omp          shared(N1,N2,thread_buff)
    im=-2**30
    ithread = 0
    !$    ithread = omp_get_thread_num ( )
    !$omp do __COLLAPSE2
    DO j2=1,n2
       DO j1=1,n1
          IF (ir(j1,j2).GE.im) THEN
             im=ir(j1,j2)
             thread_buff(1,ithread) = j1
             thread_buff(2,ithread) = j2
          ENDIF
       ENDDO
    ENDDO
    !$omp enddo
    !$omp end parallel

    ! finally find the max over threads
    i1=thread_buff(1,0)
    i2=thread_buff(2,0)
    im=ir(i1,i2)
    DO ithread=1,nthreads-1
       j1=thread_buff(1,ithread)
       j2=thread_buff(2,ithread)
       ! protect if a thread doesnt have an entry
       IF (j1>0.AND.j2>0) THEN
          IF (ir(j1,j2).GE.im) THEN
             im=ir(j1,j2)
             i1=j1
             i2=j2
          ENDIF
       ENDIF
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE iraymax
  ! ==================================================================
  FUNCTION ixypr(nfull,kr1,kr2,kr3)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nfull, kr1, kr2, kr3, ixypr

    INTEGER                                  :: i2, i3, n2

! ==--------------------------------------------------------------==

    i3=nfull/(kr1*kr2) + 1
    n2=nfull - (i3-1)*kr1*kr2
    i2=n2/kr1 + 1
    ixypr = (i3 - 1)*kr2 + i2
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION ixypr
  ! ==================================================================
  SUBROUTINE leadim(nr1,nr2,nr3,kr1,kr2,kr3)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nr1, nr2, nr3, kr1, kr2, kr3

! ==--------------------------------------------------------------==
! to align things properly for the fft we need that
! KR1=NR1+MOD(NR1,2)
! KR2=NR2+MOD(NR2,2)
! KR3=NR3+MOD(NR3,2)
! instead off that

    kr1=nr1+MOD(nr1+1,2)
    kr2=nr2+MOD(nr2+1,2)
    kr3=nr3+MOD(nr3+1,2)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE leadim
  ! ==================================================================
  SUBROUTINE xfft(iray,gvcut)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: iray(fpar%kr2,fpar%kr3)
    REAL(real_8)                             :: gvcut

    INTEGER                                  :: i, id1, id2, id3, in1, in2, &
                                                in3, ir, ix1, ix2, j, jmax, &
                                                jmin, k, kmax, kmin, nh1, &
                                                nh2, nh3
    REAL(real_8)                             :: g2, t

! ==--------------------------------------------------------------==

    CALL zeroing(iray)!,kr2*kr3)
    nh1=parm%nr1/2+1
    nh2=parm%nr2/2+1
    nh3=parm%nr3/2+1
    DO i=0,parm%nr1-1
       jmin=-parm%nr2+1
       jmax=parm%nr2-1
       IF (i.EQ.0) THEN
          jmin=0
       ENDIF
       DO j=jmin,jmax
          kmin=-parm%nr3+1
          kmax=parm%nr3-1
          IF (i.EQ.0.AND.j.EQ.0) THEN
             kmin=0
          ENDIF
          DO k=kmin,kmax
             g2=0._real_8
             DO ir=1,3
                t=REAL(i,kind=real_8)*gvec_com%b1(ir)+REAL(j,kind=real_8)*gvec_com%b2(ir)+REAL(k,kind=real_8)*gvec_com%b3(ir)
                g2=g2+t*t
             ENDDO
             IF (g2.LT.gvcut) THEN
                in1=nh1+i
                in2=nh2+j
                in3=nh3+k
                id1=2*nh1-in1
                id2=2*nh2-in2
                id3=2*nh3-in3
                ix1=iray(in2,in3)
                ix2=iray(id2,id3)
                IF (ix1.GT.0) THEN
                   iray(in2,in3)=ix1+1
                   IF (ix2.GT.0) THEN
                      IF ((in2.NE.id2).OR.(in3.NE.id3)) THEN
                         IF (paral%io_parent) WRITE(6,*) ' INCONSISTENT MESH?'
                         CALL stopgm('XFFT',' ',& 
                              __LINE__,__FILE__)
                      ENDIF
                   ENDIF
                ELSEIF (ix2.GT.0) THEN
                   iray(id2,id3)=ix2+1
                ELSE
                   iray(in2,in3)=1
                ENDIF
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE xfft
  ! ==================================================================
  SUBROUTINE cxfft(iray,gvcut)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: iray(fpar%kr2,*)
    REAL(real_8)                             :: gvcut

    INTEGER                                  :: i, icpu1, icpu2, id1, id2, &
                                                id3, in1, in2, in3, ir, j, &
                                                jmax, jmin, k, kmax, kmin, &
                                                nh1, nh2, nh3
    REAL(real_8)                             :: g2, t

! ==--------------------------------------------------------------==

    nh1=parm%nr1/2+1
    nh2=parm%nr2/2+1
    nh3=parm%nr3/2+1
    DO i=0,parm%nr1-1
       jmin=-parm%nr2+1
       jmax=parm%nr2-1
       IF (i.EQ.0) THEN
          jmin=0
       ENDIF
       DO j=jmin,jmax
          kmin=-parm%nr3+1
          kmax=parm%nr3-1
          IF (i.EQ.0.AND.j.EQ.0) THEN
             kmin=0
          ENDIF
          DO k=kmin,kmax
             g2=0._real_8
             DO ir=1,3
                t=REAL(i,kind=real_8)*gvec_com%b1(ir)+REAL(j,kind=real_8)*gvec_com%b2(ir)+REAL(k,kind=real_8)*gvec_com%b3(ir)
                g2=g2+t*t
             ENDDO
             IF (g2.LT.gvcut) THEN
                in1=nh1+i
                in2=nh2+j
                in3=nh3+k
                id1=2*nh1-in1
                id2=2*nh2-in2
                id3=2*nh3-in3
                icpu1=iray(in2,in3)
                icpu2=iray(id2,id3)
                IF (icpu2.EQ.0) THEN
                   iray(id2,id3)=icpu1
                ELSEIF (icpu1.EQ.0) THEN
                   iray(in2,in3)=icpu2
                ELSEIF (icpu1.NE.icpu2) THEN
                   IF (paral%io_parent)&
                        WRITE(6,*) ' INCONSISTENT XRAY FIELDS',icpu1,icpu2
                   CALL stopgm('CXFFT',' ',& 
                        __LINE__,__FILE__)
                ENDIF
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE cxfft
  ! ==================================================================
  SUBROUTINE gorder
    ! ==--------------------------------------------------------------==
    ! == CHECK THE G ORDERING USING MAPGP(1:NHG)                      ==
    ! == FOR A MONO-PROCESSOR, SUM_I^NHG = 1/2 * NHGS * (NHGS+1)      ==
    ! == IN PARALLEL VERSION, WE DO THE SAME THINGS:                  ==
    ! == FOR THE IP PROC, MAPGP(IG) = IG                              ==
    ! == FOR EACH JP /= IP, WE ARE LOOKING FOR THE INDEX J:           ==
    ! ==     HG_(in JP)(J) <= HG_(in IP)(IG) < HG_(in JP)(J)          ==
    ! ==     and MAPGP(IG) = MAPGP(IG) + J                            ==
    ! == IF WE SUM ALL MAPGP WE SHOULD HAVE 1/2 * NHGS* (NHGS+1)      ==
    ! ==--------------------------------------------------------------==
    CHARACTER(*), PARAMETER                  :: procedureN = 'gorder'

    INTEGER                                  :: i, ierr, ip, ipp, isub, j, l, &
                                                mho, msg1, nhgmax, nho
    REAL(real_8)                             :: xm, xt
    REAL(real_8), ALLOCATABLE                :: ho(:), hx(:)

! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    nhgmax=0
    DO ip=0,parai%nproc-1
       nhgmax=MAX(nhgmax,parap%sparm(1,ip))
    ENDDO
    ALLOCATE(ho(nhgmax),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(hx(nhgmax),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(ho)!,nhgmax)
    CALL zeroing(hx)!,nhgmax)
    CALL dcopy(ncpw%nhg,hg,1,hx,1)
    msg1 = nhgmax * 8
    CALL mp_sync(parai%allgrp)
    DO ip=1,parai%nproc-1
       ipp=MOD(parai%mepos+ip,parai%nproc)
       CALL my_shift(hx,ho,msg1,parai%mepos,-ip,parai%allgrp)
       nho=parap%sparm(1,ipp)
       l=1
       DO i=1,ncpw%nhg
          IF (hg(i).LT.ho(1)) THEN
             mho=0
          ELSE
             DO j=l,nho-1
                IF (hg(i).GE.ho(j).AND.hg(i).LT.ho(j+1)) THEN
                   l=j
                   mho=j
                   GOTO 100
                ENDIF
             ENDDO
             l=nho
             mho=nho
100          CONTINUE
          ENDIF
          mapgp(i)=mapgp(i)+mho
       ENDDO
    ENDDO
    ! ..check ordering
    xm=0._real_8
    DO i=1,ncpw%ngw
       xm=xm+mapgp(i)
    ENDDO
    CALL mp_sum(xm,parai%allgrp)
    xt=REAL(spar%ngws,kind=real_8)
    xt=0.5_real_8*xt*(xt+1._real_8)-xm
    ! IF(ABS(XT).GT.0.1_real_8.AND.PARENT) THEN
    ! WRITE(6,*) ' GORDER| PROGRAMING ERROR. INFORM THE PROGRAMMER'
    ! CALL STOPGM('GORDER','ERROR IN G-VEC ORDERING (NGW)')
    ! ENDIF
    xm=0._real_8
    DO i=1,ncpw%nhg
       xm=xm+mapgp(i)
    ENDDO
    CALL mp_sum(xm,parai%allgrp)
    xt=REAL(spar%nhgs,kind=real_8)
    xt=0.5_real_8*xt*(xt+1._real_8)-xm
    ! IF(ABS(XT).GT.0.1_real_8.AND.PARENT) THEN
    ! WRITE(6,*) ' GORDER| PROGRAMING ERROR. INFORM THE PROGRAMMER'
    ! CALL STOPGM('GORDER','ERROR IN G-VEC ORDERING (NHG)')
    ! ENDIF
    DEALLOCATE(ho,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(hx,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gorder
  ! ==================================================================

  SUBROUTINE gsort(hg,inyh,nhg)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: inyh(3,*), nhg
    REAL(real_8)                             :: hg(nhg)

    CHARACTER(*), PARAMETER                  :: procedureN = 'gsort'

    INTEGER                                  :: icurr, ierr, ig, it
    INTEGER, ALLOCATABLE                     :: INDEX(:)

! REORDER THE G S IN ORDER OF INCREASING MAGNITUDE

    ALLOCATE(INDEX(nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL sort2(hg,nhg,index)
    DO 25 ig=1,nhg-1
       icurr=ig
30     CONTINUE
       IF (INDEX(icurr).NE.ig) THEN
          it=inyh(1,icurr)
          inyh(1,icurr)=inyh(1,INDEX(icurr))
          inyh(1,INDEX(icurr))=it
          it=inyh(2,icurr)
          inyh(2,icurr)=inyh(2,INDEX(icurr))
          inyh(2,INDEX(icurr))=it
          it=inyh(3,icurr)
          inyh(3,icurr)=inyh(3,INDEX(icurr))
          inyh(3,INDEX(icurr))=it
          it=icurr
          icurr=INDEX(icurr)
          INDEX(it)=it
          IF (INDEX(icurr).EQ.ig) THEN
             INDEX(icurr)=icurr
             GOTO 25
          ENDIF
          GOTO 30
       ENDIF
25  CONTINUE
    DEALLOCATE(index,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gsort

END MODULE loadpa_utils
