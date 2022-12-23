#include "cpmd_global.h"

MODULE jrotation_utils
  USE cnst,                            ONLY: pi
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE machine,                         ONLY: m_flush
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_cart,&
                                             mp_environ,&
                                             mp_max,&
                                             mp_sendrecv,&
                                             mp_sum
  USE parac,                           ONLY: parai,&
                                             paral
  USE prng_utils,                      ONLY: repprngu
  USE sd_wannier_utils,                ONLY: exponentiate
  USE system,                          ONLY: cntl,&
                                             paraw,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: rmatmov,&
                                             rtrans,&
                                             unitmx
  USE wann,                            ONLY: wan05,&
                                             wannc,&
                                             wanni,&
                                             wannr
!!use utils, only : matmov
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: jrotation
  PUBLIC :: set_orbdist
  PUBLIC :: my_set_orbdist
  !public :: jrotation1
  !public :: jrotationp
  !public :: jacobirot_new
  !public :: jacobirot
  !public :: calc_angle
  !public :: xyrot
  !public :: xyrot_row
  !public :: xyrot_col
  !public :: rmrot
  !public :: eberlein
  PUBLIC :: give_scr_jrot
  !public :: xapyb
  !public :: rmupx
  !public :: music

CONTAINS

  ! ==================================================================
  SUBROUTINE jrotation(rotmat,xyzmat,ldx,nstate)
    ! ==--------------------------------------------------------------==
    ! == Jacobi rotation procedure for Wannier functions & centers    ==
    ! ==                                                              ==
    ! ==  Present release:                                            ==
    ! ==             Strasbourg/Tokyo/Zurich/Hyogo, 23 May 2013       ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ldx
    COMPLEX(real_8)                          :: xyzmat(ldx,ldx,*)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: rotmat(nstate,*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'jrotation'

    INTEGER                                  :: i, ierr, imax, isub, j, k, &
                                                loc_nogrp, msglen, my_grp, &
                                                my_nproc
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: loc_nolist, loc_nplist
    LOGICAL                                  :: debug
    LOGICAL, SAVE                            :: is_first = .TRUE.
    REAL(real_8), ALLOCATABLE                :: abc(:,:,:), gmat(:,:)

!nwanopt)
!nstate)
! ==--------------------------------------------------------------==
! ==--------------------------------------------------------------==

    debug=.FALSE.
    CALL tiset(procedureN,isub)
    ! ==--------------------------------------------------------------==
    IF (wanni%w_type.EQ.2 .AND. paral%io_parent) THEN
       WRITE(6,'(A,A)') ' LOCALIZATION| ORBITAL ROTATIONS ',&
            'WITH RESTA FUNCTIONAL NOT POSSIBLE '
       WRITE(6,'(A,A)') ' LOCALIZATION| USE VANDERBILT FUNCTIONAL '
       CALL stopgm(procedureN,'WRONG FUNCTIONAL',& 
            __LINE__,__FILE__)
    ENDIF
    ! ..initialise global rotation matrix
    CALL unitmx(rotmat,nstate)
    IF (wannr%w_ran.GT.0._real_8) THEN
       IF (paral%io_parent) THEN
          ALLOCATE(abc(nstate,nstate,3),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
          ALLOCATE(gmat(nstate,nstate),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
          ! ..randomly distort initial vectors
          DO i=1,nstate
             gmat(i,i)=0._real_8
             DO j=i+1,nstate
                gmat(i,j)=repprngu()*wannr%w_ran
                gmat(j,i)=-gmat(i,j)
             ENDDO
          ENDDO
          CALL exponentiate(gmat,rotmat,nstate)
          CALL xyz_update(rotmat,xyzmat,abc,ldx,nstate)
          DEALLOCATE(gmat,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
               __LINE__,__FILE__)
          DEALLOCATE(abc,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
               __LINE__,__FILE__)
       ENDIF
       CALL mp_bcast(rotmat,nstate**2,parai%io_source,parai%cp_grp)
       DO k=1,wannc%nwanopt
          CALL mp_bcast(xyzmat(:,:,k),ldx*(nstate-1)+nstate,parai%io_source,parai%cp_grp)
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (is_first) THEN
       is_first=.FALSE.

       loc_nogrp = 0
       IF (paral%io_parent.AND.wan05%loc_npgrp.LE.0) CALL stopgm(procedureN,&
            'number of localization processes is less than zero.',& 
            __LINE__,__FILE__)
       IF (paral%io_parent.AND.wan05%loc_npgrp.GT.parai%cp_nproc) THEN
          WRITE(6,'(1X,A)') 'WARNING: the number of'&
               //'localization processes is greater than the '&
               //'available number of processes.'
          WRITE(6,'(1X,A)') 'WARNING: the number of'&
               //'localization processes will be reset! '
          CALL m_flush(6)
          wan05%loc_npgrp=0
       ELSEIF (paral%io_parent.AND.MOD(parai%cp_nproc,parai%cp_nogrp).NE.0) THEN
          WRITE(6,'(1X,A)') 'WARNING: the number of'&
               //'localization processes is not a divisor '&
               //'of the available number of processes.'
          WRITE(6,'(1X,A)') 'WARNING: the number of'&
               //'localization processes will be reset! '
          CALL m_flush(6)
          wan05%loc_npgrp=0
       ELSE
          loc_nogrp = parai%cp_nproc / wan05%loc_npgrp
       ENDIF

       ALLOCATE(loc_nolist(parai%cp_nproc),loc_nplist(parai%cp_nproc),stat=ierr)
       IF (ierr.NE.0) CALL stopgm(procedureN,'Allocation problem',& 
            __LINE__,__FILE__)
       !      CALL mp_cart(mp_comm_world,loc_nogrp,wan05%loc_npgrp,loc_nolist,loc_nplist,&
       CALL mp_cart(parai%cp_grp,loc_nogrp,wan05%loc_npgrp,loc_nolist,loc_nplist,&
            parai%loc_inter_grp,parai%loc_grp)
       DEALLOCATE(loc_nolist,loc_nplist,stat=ierr)
       IF (ierr.NE.0) CALL stopgm(procedureN,'Deallocation problem',& 
            __LINE__,__FILE__)
       CALL mp_environ(parai%loc_grp,parai%loc_nproc,parai%loc_me)
    ENDIF
    my_grp=parai%loc_grp
    my_nproc=parai%loc_nproc


    IF (my_nproc.EQ.1 .OR. cntl%twserial) THEN
       IF (paral%io_parent)&
            CALL jrotation1(rotmat,xyzmat,ldx,nstate)
       CALL mp_bcast(rotmat,nstate**2,parai%io_source,parai%cp_grp)
    ELSE
       CALL set_orbdist(nstate,1,my_nproc,imax)
       CALL jrotationp(rotmat(1,1),xyzmat(1,1,1),ldx,nstate,my_grp)
       ! redist
       CALL mp_bcast(rotmat,nstate**2,parai%io_source,parai%cp_grp)
       DO k=1,wannc%nwanopt
          msglen=MAX(ldx*(nstate-1)+nstate,0) !MAX to avoid negative value
          CALL mp_bcast(xyzmat(:,:,k),msglen,parai%io_source,parai%cp_grp)
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE jrotation
  ! ==================================================================
  SUBROUTINE set_orbdist(nstate,nblock,my_nproc,nbmax)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate, nblock, my_nproc, &
                                                nbmax

    INTEGER                                  :: ip, nx
    REAL(real_8)                             :: xsaim, xsnow, xstates

! ==--------------------------------------------------------------==

    nbmax=0
    xstates=REAL(nblock,kind=real_8)
    IF ((xstates*my_nproc).LT.nstate) THEN
       xstates=REAL(nstate,kind=real_8)/REAL(my_nproc,kind=real_8)
    ENDIF
    xsnow=0.0_real_8
    xsaim=0.0_real_8
    DO ip=1,my_nproc
       xsaim = xsnow + xstates
       paraw%nwa12(ip-1,1)=NINT(xsnow)+1
       paraw%nwa12(ip-1,2)=NINT(xsaim)
       IF (NINT(xsaim).GT.nstate) THEN
          paraw%nwa12(ip-1,2)=nstate
       ENDIF
       IF (NINT(xsnow).GT.nstate) THEN
          paraw%nwa12(ip-1,1)=nstate+1
       ENDIF
       xsnow = xsaim
    ENDDO
    DO ip=0,my_nproc-1
       nx=paraw%nwa12(ip,2)-paraw%nwa12(ip,1)+1
       nbmax=MAX(nbmax,nx)
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE set_orbdist
  ! ==================================================================
  SUBROUTINE my_set_orbdist(nstate,nblock,nbmax)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate, nblock, nbmax

    INTEGER                                  :: ip, n
    REAL(real_8)                             :: xsaim, xsnow, xstates

! ==--------------------------------------------------------------==

    nbmax=0
    xstates=REAL(nstate,kind=real_8)/REAL(parai%nproc,kind=real_8)
    IF (xstates.LT.REAL(nblock,kind=real_8)) xstates=REAL(nblock,kind=real_8)
    xsnow=0.0_real_8
    n=0
    DO ip=1,parai%nproc
       xsaim = xsnow + xstates
       paraw%nwa12(ip-1,1)=NINT(xsnow)+1
       paraw%nwa12(ip-1,2)=NINT(xsaim)
       IF (NINT(xsaim).GT.nstate) paraw%nwa12(ip-1,2)=nstate
       IF (ip.EQ.parai%nproc) paraw%nwa12(ip-1,2)=nstate
       IF (NINT(xsnow).GT.nstate) paraw%nwa12(ip-1,1)=nstate+1
       xsnow = xsaim
       n=n+paraw%nwa12(ip-1,2)-paraw%nwa12(ip-1,1)+1
       nbmax=MAX(nbmax,paraw%nwa12(ip-1,2)-paraw%nwa12(ip-1,1)+1)
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE my_set_orbdist
  ! ==================================================================
  SUBROUTINE jrotation1(rotmat,xyzmat,ldx,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ldx
    COMPLEX(real_8)                          :: xyzmat(ldx,ldx,*)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: rotmat(nstate,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'jrotation1'

    INTEGER                                  :: i, idamax, ierr, imax, iopt, &
                                                isub, j
    INTEGER, ALLOCATABLE                     :: bot(:), top(:)
    LOGICAL                                  :: debug
    REAL(real_8)                             :: gmax
    REAL(real_8), ALLOCATABLE                :: gmat(:)

! ==--------------------------------------------------------------==
! ==--------------------------------------------------------------==

    debug=.FALSE.
    CALL tiset(procedureN,isub)
    ALLOCATE(gmat(nstate*nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    gmax=0._real_8
    ! ==--------------------------------------------------------------==
    ALLOCATE(top(nstate),bot(nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    j=1
    DO i=1,nstate,2
       top(j) = i
       j = j+1
    ENDDO
    j=1
    DO i=2,nstate,2
       bot(j) = i
       j = j+1
    ENDDO
    ! ==--------------------------------------------------------------==
    DO iopt=1,wanni%w_maxs
       ! ..do a Jacobi rotation for each pair of orbitals
       IF (.TRUE.) THEN
          CALL jacobirot(rotmat,xyzmat,ldx,nstate)
       ELSE
          CALL jacobirot_new(rotmat,xyzmat,ldx,nstate,top,bot)
       ENDIF
       ! ..calculate the gradient
       CALL xgradat0(wanni%w_type,gmat,xyzmat,ldx,nstate)
       imax=idamax(nstate*nstate,gmat,1)
       gmax=ABS(gmat(imax))/parm%alat**2
       IF (paral%io_parent.AND.(iopt==1.OR.MOD(iopt,5)==0.OR.gmax.LT.wannr%w_eps))&
            WRITE(6,'(A,I4,A,G12.3)')' WANNIER CODE| IOPT=',iopt,' GMAX=',gmax
       IF (gmax.LT.wannr%w_eps) GOTO 100
    ENDDO
    IF (paral%io_parent)&
         WRITE(6,'(A,G12.3)') ' WANNIER CODE| NO CONVERGENCE (GMAX) ',gmax
100 CONTINUE

    DEALLOCATE(gmat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(top,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(bot,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE jrotation1
  ! ==================================================================
  SUBROUTINE jrotationp(rotmat,xyzmat,ldx,nstate,my_grp)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ldx
    COMPLEX(real_8)                          :: xyzmat(ldx,ldx,*)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: rotmat(nstate,*)
    INTEGER                                  :: my_grp

    CHARACTER(*), PARAMETER                  :: procedureN = 'jrotationp'

    COMPLEX(real_8)                          :: xa
    COMPLEX(real_8), ALLOCATABLE             :: xyzc(:), xyzd(:), xyzm(:)
    INTEGER :: i, i1, i2, ierr, ii, iopt, ip, ipo, is, isub, j, k, kk, kk2, &
      m, mep, mia, mib, mip, mja, mjb, mre, msglen, my_i, my_j, my_me, &
      my_nproc, ngmax, nipo, nn, nomax, nome, npair, nsweep, rcount, scount
    INTEGER, ALLOCATABLE                     :: bot(:), lppl(:,:), top(:)
    LOGICAL                                  :: debug
    REAL(real_8)                             :: dgmax, gcomp, gmax, gmax0
    REAL(real_8), ALLOCATABLE                :: gmat(:), rmat(:)

!nstate)
! ==--------------------------------------------------------------==
! ==--------------------------------------------------------------==

    debug=.FALSE.
    CALL tiset(procedureN,isub)
    gmax=0._real_8
    gmax0=HUGE(0.0_real_8)
    ngmax=0
    CALL mp_environ(my_grp,my_nproc,my_me)
    ALLOCATE(top(nstate),bot(nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    top(:)=HUGE(0)
    bot(:)=HUGE(0)
    ! ==--------------------------------------------------------------==
    nome=paraw%nwa12(my_me,2)-paraw%nwa12(my_me,1)+1
    nomax=0
    DO ip=0,my_nproc-1
       nn=paraw%nwa12(ip,2)-paraw%nwa12(ip,1)+1
       nomax=MAX(nomax,nn)
    ENDDO
    nn=MAX(8*nomax*nomax*wannc%nwanopt,2*nomax*nstate*wannc%nwanopt,&
         2*nomax*my_nproc*wannc%nwanopt)
    ALLOCATE(xyzm(nn),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(xyzc(nn),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(xyzd(nomax*nstate*wannc%nwanopt),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(rmat(4*nomax*nomax),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    nn=MAX(4*nomax*nomax,nstate*nstate)
    ALLOCATE(gmat(nn),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(lppl(2,my_nproc),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL zeroing(gmat)!,nn)
    CALL zeroing(lppl)!,2*my_nproc)
    ! ==--------------------------------------------------------------==
    DO iopt=1,wanni%w_maxs
       nsweep=my_nproc-MOD(my_nproc+1,2)
       npair=(my_nproc+1)/2
       DO is=1,nsweep
          CALL eberlein(is,my_nproc,lppl)
          ! look for my partner in this round
          ipo=-1
          DO ip=1,npair
             IF (lppl(1,ip).EQ.my_me+1) THEN
                ipo=lppl(2,ip)-1
                GOTO 1000
             ELSEIF (lppl(2,ip).EQ.my_me+1) THEN
                ipo=lppl(1,ip)-1
                GOTO 1000
             ENDIF
          ENDDO
1000      CONTINUE
          IF (ipo.GE.0) THEN
             nipo=paraw%nwa12(ipo,2)-paraw%nwa12(ipo,1)+1
          ELSE
             nipo=0
          ENDIF
          IF (nipo*nome.GT.0) THEN
             nn=nipo+nome
             CALL unitmx(rmat,nn)
             ! pack data to be send to partner
             DO k=1,wannc%nwanopt
                kk=(k-1)*nome*nstate
                CALL matmov(nstate,nome,xyzmat(1,paraw%nwa12(my_me,1),k),&
                     ldx,xyzc(kk+1),nstate)
             ENDDO
             ! send/receive data
             scount=wannc%nwanopt*nome*nstate
             rcount=wannc%nwanopt*nipo*nstate
             CALL mp_sendrecv(xyzc,scount,ipo,xyzd,rcount,ipo,my_grp)
             ! build up local problem
             IF (paraw%nwa12(my_me,1).LT.paraw%nwa12(ipo,1)) THEN
                mia=1
                mib=nome+1
                mja=nome*nn+nome+1
                mjb=nome*nn+1
             ELSE
                mia=nipo*nn+nipo+1
                mib=nipo*nn+1
                mja=1
                mjb=nipo+1
             ENDIF
             DO k=1,wannc%nwanopt
                kk=(k-1)*nome*nstate
                kk2=(k-1)*nn*nn
                j=paraw%nwa12(my_me,1)
                CALL matmov(nome,nome,xyzc(j+kk),nstate,&
                     xyzm(kk2+mia),nn)
                j=paraw%nwa12(ipo,1)
                CALL matmov(nipo,nome,xyzc(j+kk),nstate,&
                     xyzm(kk2+mib),nn)
                kk=(k-1)*nipo*nstate
                j=paraw%nwa12(ipo,1)
                CALL matmov(nipo,nipo,xyzd(j+kk),nstate,&
                     xyzm(kk2+mja),nn)
                j=paraw%nwa12(my_me,1)
                CALL matmov(nome,nipo,xyzd(j+kk),nstate,&
                     xyzm(kk2+mjb),nn)
             ENDDO
             ! solve local problem
             IF (.TRUE.) THEN
                CALL jacobirot(rmat,xyzm,nn,nn)
             ELSE
                my_j=1
                DO my_i=1,nn,2
                   top(my_j) = my_i
                   my_j = my_j+1
                ENDDO
                my_j=1
                DO my_i=2,nn,2
                   bot(my_j) = my_i
                   my_j = my_j+1
                ENDDO
                CALL jacobirot_new(rmat,xyzm,nn,nn,top,bot)
                top(1:nn)=HUGE(0)
                bot(1:nn)=HUGE(0)
             ENDIF
             ! send/receive ROTMAT data
             scount=nome*nstate
             rcount=nipo*nstate
             CALL mp_sendrecv(rotmat(:,paraw%nwa12(my_me,1)),scount,ipo,&
                  rotmat(:,paraw%nwa12(ipo,1)),rcount,ipo,my_grp)
             ! update ROTMAT
             CALL dgemm("N","N",nstate,nome,nipo,1._real_8,&
                  rotmat(1,paraw%nwa12(ipo,1)),nstate,rmat(mib),nn,&
                  0._real_8,gmat,nstate)
             CALL dgemm("N","N",nstate,nome,nome,1._real_8,&
                  rotmat(1,paraw%nwa12(my_me,1)),nstate,rmat(mia),nn,&
                  1._real_8,gmat,nstate)
             CALL dcopy(nstate*nome,gmat,1,rotmat(1,paraw%nwa12(my_me,1)),1)
             ! update XYZMAT (first part, multiplication from right)
             DO k=1,wannc%nwanopt
                kk=(k-1)*nome*nstate+1
                kk2=(k-1)*nipo*nstate+1
                CALL xapyb(nstate,xyzd(kk2),rmat(mib),nipo,xyzc(kk),&
                     rmat(mia),nome,xyzm(kk),nn)
             ENDDO
          ELSEIF (nome.GT.0) THEN
             DO k=1,wannc%nwanopt
                kk=(k-1)*nome*nstate
                CALL matmov(nstate,nome,xyzmat(1,paraw%nwa12(my_me,1),k),&
                     ldx,xyzm(kk+1),nstate)
             ENDDO
          ENDIF
          ! distribute rotation matrices and update xyzmat (second part)
          IF (nome.GT.0) THEN
             DO k=1,wannc%nwanopt
                DO i=paraw%nwa12(my_me,1),paraw%nwa12(my_me,2)
                   !CALL zeroing(xyzmat(:,i,k))!,nstate)
                   xyzmat(1:nstate,i,k)=0.0_real_8
                ENDDO
             ENDDO
          ENDIF
          ! transpose rotation matrix and pack for distribution
          IF (nome*nipo.NE.0) THEN
             CALL rtrans(rmat,nn)
             CALL rmatmov(nome,nome,rmat(mia),nn,gmat(1),nome)
             kk=nomax*nomax+1
             CALL rmatmov(nome,nipo,rmat(mjb),nn,gmat(kk),nome)
          ELSE
             CALL zeroing(gmat)!,2*nomax*nomax)
          ENDIF
          DO ip=0,my_nproc-1
             k=2*nomax*nomax+1
             IF (ip.EQ.0) THEN
                CALL dcopy(2*nomax*nomax,gmat,1,gmat(k),1)
                mep = my_me
             ELSE
                msglen=2*nomax*nomax
                mep=MOD(my_me-ip+my_nproc,my_nproc)
                mre=MOD(my_me+ip+my_nproc,my_nproc)
                CALL mp_sendrecv(gmat,msglen,mre,gmat(k:),&
                     msglen,mep,my_grp)
             ENDIF
             mre=paraw%nwa12(mep,2)-paraw%nwa12(mep,1)+1
             IF (nome.GT.0) THEN
                IF (mre.GT.0) THEN
                   ! look for partner of mep
                   mip=-1
                   DO i=1,npair
                      IF (lppl(1,i).EQ.mep+1) THEN
                         mip=lppl(2,i)-1
                         GOTO 1001
                      ELSEIF (lppl(2,i).EQ.mep+1) THEN
                         mip=lppl(1,i)-1
                         GOTO 1001
                      ENDIF
                   ENDDO
1001               CONTINUE
                   IF (mip.GE.0) THEN
                      m=paraw%nwa12(mip,2)-paraw%nwa12(mip,1)+1
                   ELSE
                      m=0
                   ENDIF
                   IF (m.GT.0) THEN
                      j=paraw%nwa12(my_me,1)
                      i1=paraw%nwa12(mep,1)
                      i2=paraw%nwa12(mip,1)
                      k=2*nomax*nomax+1
                      kk=3*nomax*nomax+1
                      DO i=1,wannc%nwanopt
                         ii=(i-1)*nome*nstate
                         CALL rmupx(mre,mre,gmat(k),xyzm(ii+i1),nstate,nome,&
                              xyzmat(i1,j,i),ldx)
                         CALL rmupx(mre,m,gmat(kk),xyzm(ii+i2),nstate,nome,&
                              xyzmat(i1,j,i),ldx)
                      ENDDO
                   ELSE
                      DO k=1,wannc%nwanopt
                         j=paraw%nwa12(my_me,1)
                         i=paraw%nwa12(mep,1)
                         kk=(k-1)*nome*nstate+i
                         CALL matmov(mre,nome,xyzm(kk),nstate,&
                              xyzmat(i,j,k),ldx)
                      ENDDO
                   ENDIF
                ENDIF
             ENDIF
          ENDDO
       ENDDO
       ! ..calculate the max gradient
       ! store the diagonal of xyzmat
       CALL zeroing(xyzc(1:wannc%nwanopt*nomax))!,nwanopt*nomax)
       !$omp parallel do private(K,I,KK,II)
       DO k=1,wannc%nwanopt
          kk=(k-1)*nomax
          DO i=paraw%nwa12(my_me,1),paraw%nwa12(my_me,2)
             ii=i-paraw%nwa12(my_me,1)+1
             xyzc(kk+ii)=xyzmat(i,i,k)
          ENDDO
       ENDDO
       ! collect all the diagonal elements
       msglen=2*wannc%nwanopt*nomax * 8
       CALL my_concat(xyzc,xyzm,msglen,my_grp)
       DO ip=0,my_nproc-1
          j=ip*wannc%nwanopt*nomax
          !$omp parallel do private(K,I,KK,I1,II)
          DO k=1,wannc%nwanopt
             kk=(k-1)*nomax+j
             i1=(k-1)*nstate
             DO i=paraw%nwa12(ip,1),paraw%nwa12(ip,2)
                ii=i-paraw%nwa12(ip,1)+1
                xyzc(i1+i)=xyzm(kk+ii)
             ENDDO
          ENDDO
       ENDDO
       ! calculate gradient
       gmax=0._real_8
       DO j=paraw%nwa12(my_me,1),paraw%nwa12(my_me,2)
          DO i=1,j-1
             gcomp=0._real_8
             DO k=1,wannc%nwanopt
                i1=(k-1)*nstate
                xa=xyzc(i1+j)-xyzc(i1+i)
                gcomp=gcomp+wannc%wwei(k)*&
                     REAL(4._real_8*CONJG(xyzmat(i,j,k))*xa)
             ENDDO
             gmax=MAX(gmax,ABS(gcomp))
          ENDDO
       ENDDO
       ! call my_max_d(GMAX,1,allgrp)
       gmax=ABS(gmax)
       CALL mp_max(gmax,my_grp)
       ! 
       gmax=gmax/parm%alat**2
       dgmax=ABS(gmax-gmax0)
       IF (dgmax.LT.1.0e-2_real_8*wannr%w_eps) THEN
          ngmax=ngmax+1
       ELSE
          ngmax=0
       ENDIF
       gmax0=gmax
       IF (paral%io_parent.AND.(iopt==1.OR.MOD(iopt,5)==0&
            .OR.gmax.LT.wannr%w_eps.OR.ngmax.GE.1)) THEN
          WRITE(6,'(A,I4,A,G12.3)')' WANNIER CODE| IOPT=',iopt,' GMAX=',gmax
       ENDIF

       IF (gmax.LT.wannr%w_eps.OR.ngmax.GE.2) GOTO 100
    ENDDO
    IF (paral%io_parent) WRITE(6,'(A,G12.3)')&
         ' WANNIER CODE| NO CONVERGENCE (GMAX) ',gmax
100 CONTINUE
    ! update the full xyzmat and rotmat
    DO k=1,wannc%nwanopt
       DO i=1,paraw%nwa12(my_me,1)-1
          !CALL zeroing(xyzmat(:,i,k))!,nstate)
          xyzmat(1:nstate,i,k)=0.0_real_8
       ENDDO
       DO i=paraw%nwa12(my_me,2)+1,nstate
          !CALL zeroing(xyzmat(:,i,k))!,nstate)
          xyzmat(1:nstate,i,k)=0.0_real_8
       ENDDO
       i=MAX(((nstate-1)*ldx+nstate),0) ! MAX to avoid negative value
       CALL mp_sum(xyzmat(:,:,k),i,my_grp)
    ENDDO
    DO i=1,paraw%nwa12(my_me,1)-1
       !CALL zeroing(rotmat(:,i))!,nstate)
       rotmat(1:nstate,i)=0._real_8
    ENDDO
    DO i=paraw%nwa12(my_me,2)+1,nstate
       !CALL zeroing(rotmat(:,i))!,nstate)
       rotmat(1:nstate,i)=0._real_8
    ENDDO
    CALL mp_sum(rotmat,nstate*nstate,my_grp)
    DEALLOCATE(top,bot,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(xyzm,xyzc,xyzd,rmat,gmat,lppl,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE jrotationp
  ! ==================================================================
  SUBROUTINE jacobirot_new(rotmat,xyzmat,ldx,nstate,top,bot)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ldx
    COMPLEX(real_8)                          :: xyzmat(ldx,ldx,*)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: rotmat(nstate,nstate)
    INTEGER                                  :: top(nstate), bot(nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'jacobirot_new'
    REAL(real_8), PARAMETER                  :: eps = 1.0e-12_real_8

    INTEGER                                  :: hmn, i, ierr, isub, k, p, q, &
                                                set
    REAL(real_8)                             :: alfa, cg, sg
    REAL(real_8), ALLOCATABLE, DIMENSION(:)  :: alfas

! ==--------------------------------------------------------------==
! orig: 12
! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    ALLOCATE(alfas(nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    ! ..do a Jacobi rotation for each pair of orbitals

    IF (MOD(nstate,2).NE.0) THEN
       hmn = nstate-1
    ELSE
       hmn = nstate
    ENDIF

    DO set=1,hmn-1
       !$omp parallel do private(K,p,q,CG,SG)
       DO k=1,hmn/2
          p = MIN(top(k),bot(k))
          q = MAX(top(k),bot(k))
          ! PRINT *, top(1:2)
          ! PRINT *, bot(1:2)

          ! ..calculate optimal rotation angle
          CALL calc_angle(xyzmat,ldx,p,q,alfas(k))

          IF (ABS(alfas(k)).GT.eps) THEN
             cg=COS(alfas(k))
             sg=SIN(alfas(k))
             CALL rmrot(nstate,rotmat,p,q,cg,sg)
             CALL xyrot_col(nstate,xyzmat,ldx,p,q,cg,sg)
          ENDIF
       ENDDO
       !$omp parallel do private(K,p,q,CG,SG)
       DO k=1,hmn/2
          p = MIN(top(k),bot(k))
          q = MAX(top(k),bot(k))

          IF (ABS(alfas(k)).GT.eps) THEN
             cg=COS(alfas(k))
             sg=SIN(alfas(k))
             CALL xyrot_row(nstate,xyzmat,ldx,p,q,cg,sg)
          ENDIF
       ENDDO

       CALL music(top,bot,hmn)

    ENDDO
    IF (hmn.LT.nstate) THEN
       DO i=1, nstate-1
          CALL calc_angle(xyzmat,ldx,i,nstate,alfa)
          IF (ABS(alfa).GT.eps) THEN
             cg=COS(alfa)
             sg=SIN(alfa)
             CALL rmrot(nstate,rotmat,i,nstate,cg,sg)
             CALL xyrot(nstate,xyzmat,ldx,i,nstate,cg,sg)
          ENDIF
       ENDDO
    ENDIF
    DEALLOCATE(alfas,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE jacobirot_new
  ! ==================================================================
  SUBROUTINE jacobirot(rotmat,xyzmat,ldx,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ldx
    COMPLEX(real_8)                          :: xyzmat(ldx,ldx,*)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: rotmat(nstate,nstate)

    REAL(real_8), PARAMETER                  :: eps = 1.0e-12_real_8

    INTEGER                                  :: i, j
    REAL(real_8)                             :: alfa, cg, sg

! ==--------------------------------------------------------------==
! orig: 12
! ==--------------------------------------------------------------==
! ..do a Jacobi rotation for each pair of orbitals

    DO i=1,nstate
       DO j=i+1,nstate
          ! ..calculate optimal rotation angle
          CALL calc_angle(xyzmat,ldx,i,j,alfa)
          IF (ABS(alfa).GT.eps) THEN
             cg=COS(alfa)
             sg=SIN(alfa)
             CALL rmrot(nstate,rotmat,i,j,cg,sg)
             CALL xyrot(nstate,xyzmat,ldx,i,j,cg,sg)
          ENDIF
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE jacobirot
  ! ==================================================================
  SUBROUTINE calc_angle(xyzmat,ldx,i,j,alfa)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ldx
    COMPLEX(real_8)                          :: xyzmat(ldx,ldx,*)
    INTEGER                                  :: i, j
    REAL(real_8)                             :: alfa

    REAL(real_8), PARAMETER                  :: pi4 = 0.25_real_8*pi 

    COMPLEX(real_8)                          :: yiia, yiib, yija, yijb, yjja, &
                                                yjjb
    INTEGER                                  :: ia
    REAL(real_8)                             :: a12, b12, d2, tan4a

! ==--------------------------------------------------------------==

    a12=0._real_8
    b12=0._real_8
    !CDIR NOVECTOR
    DO ia=1,wannc%nwanopt
       yiia=CONJG(xyzmat(i,i,ia))
       yjja=CONJG(xyzmat(j,j,ia))
       yija=CONJG(xyzmat(i,j,ia))
       yiib=xyzmat(i,i,ia)
       yjjb=xyzmat(j,j,ia)
       yijb=xyzmat(i,j,ia)
       a12=a12 + wannc%wwei(ia)*REAL(yija*(yiib-yjjb))
       b12=b12 + wannc%wwei(ia)*REAL(yija*yijb -&
            0.25_real_8*(yiia-yjja)*(yiib-yjjb))
    ENDDO
    IF (ABS(b12).GT.1.e-10_real_8) THEN
       tan4a=-a12/b12
       alfa=0.25_real_8*ATAN(tan4a)
    ELSEIF (ABS(a12).LT.1.e-10_real_8) THEN
       alfa=0._real_8
       b12=0._real_8
    ELSE
       alfa=pi4
    ENDIF
    d2=a12*SIN(4._real_8*alfa)-b12*COS(4._real_8*alfa)
    IF (d2.LE.0._real_8) alfa=alfa+pi4
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE calc_angle
  ! ==================================================================
  SUBROUTINE xyrot(nstate,xyzmat,ldx,i,j,cg,sg)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate, ldx
    COMPLEX(real_8)                          :: xyzmat(ldx,ldx,*)
    INTEGER                                  :: i, j
    REAL(real_8)                             :: cg, sg

    COMPLEX(real_8)                          :: zi, zj
    INTEGER                                  :: k, l

! ==--------------------------------------------------------------==

    !$omp parallel do private(K,L,ZI,ZJ)
    DO k=1,wannc%nwanopt
       DO l=1,nstate
          zi=cg*xyzmat(l,i,k)+sg*xyzmat(l,j,k)
          zj=-sg*xyzmat(l,i,k)+cg*xyzmat(l,j,k)
          xyzmat(l,i,k)=zi
          xyzmat(l,j,k)=zj
       ENDDO
       DO l=1,nstate
          zi=cg*xyzmat(i,l,k)+sg*xyzmat(j,l,k)
          zj=-sg*xyzmat(i,l,k)+cg*xyzmat(j,l,k)
          xyzmat(i,l,k)=zi
          xyzmat(j,l,k)=zj
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE xyrot
  ! ==================================================================
  SUBROUTINE xyrot_row(nstate,xyzmat,ldx,i,j,cg,sg)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate, ldx
    COMPLEX(real_8)                          :: xyzmat(ldx,ldx,*)
    INTEGER                                  :: i, j
    REAL(real_8)                             :: cg, sg

    COMPLEX(real_8)                          :: zi, zj
    INTEGER                                  :: k, l

! ==--------------------------------------------------------------==

    DO k=1,wannc%nwanopt
       DO l=1,nstate
          zi=cg*xyzmat(i,l,k)+sg*xyzmat(j,l,k)
          zj=-sg*xyzmat(i,l,k)+cg*xyzmat(j,l,k)
          xyzmat(i,l,k)=zi
          xyzmat(j,l,k)=zj
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE xyrot_row
  ! ==================================================================
  SUBROUTINE xyrot_col(nstate,xyzmat,ldx,i,j,cg,sg)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate, ldx
    COMPLEX(real_8)                          :: xyzmat(ldx,ldx,*)
    INTEGER                                  :: i, j
    REAL(real_8)                             :: cg, sg

    COMPLEX(real_8)                          :: zi, zj
    INTEGER                                  :: k, l

! ==--------------------------------------------------------------==

    DO k=1,wannc%nwanopt
       DO l=1,nstate
          zi=cg*xyzmat(l,i,k)+sg*xyzmat(l,j,k)
          zj=-sg*xyzmat(l,i,k)+cg*xyzmat(l,j,k)
          xyzmat(l,i,k)=zi
          xyzmat(l,j,k)=zj
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE xyrot_col
  ! ==================================================================
  SUBROUTINE rmrot(nstate,rotmat,i,j,cg,sg)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    REAL(real_8)                             :: rotmat(nstate,nstate)
    INTEGER                                  :: i, j
    REAL(real_8)                             :: cg, sg

    INTEGER                                  :: l
    REAL(real_8)                             :: zi, zj

! ==--------------------------------------------------------------==

    DO l=1,nstate
       zi=cg*rotmat(l,i)+sg*rotmat(l,j)
       zj=-sg*rotmat(l,i)+cg*rotmat(l,j)
       rotmat(l,i)=zi
       rotmat(l,j)=zj
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rmrot
  ! ==================================================================
  SUBROUTINE eberlein(iter,nlist,lppl)
    INTEGER                                  :: iter, nlist, lppl(2,*)

    INTEGER                                  :: i, ii, jj, ll, npair

! ==--------------------------------------------------------------==

    IF (nlist.LE.1) CALL stopgm('EBERLEIN','NLIST TOO SMALL',& 
         __LINE__,__FILE__)
    npair=(nlist+1)/2
    IF (iter.EQ.1) THEN
       ! ..set up initial ordering
       DO i=1,nlist
          ii=(i+1)/2
          jj=MOD(i+1,2)+1
          lppl(jj,ii)=i
       ENDDO
       IF (MOD(nlist,2).EQ.1) lppl(2,npair)=-1
    ELSEIF (MOD(iter,2).EQ.0) THEN
       ! ..a type shift
       ll=lppl(1,npair)
       DO i=npair,3,-1
          lppl(1,i)=lppl(1,i-1)
       ENDDO
       lppl(1,2)=lppl(2,1)
       lppl(2,1)=ll
    ELSE
       ! ..b type shift
       ll=lppl(2,1)
       DO i=1,npair-1
          lppl(2,i)=lppl(2,i+1)
       ENDDO
       lppl(2,npair)=ll
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE eberlein
  ! ==================================================================
  SUBROUTINE give_scr_jrot(lwann,tag,nstate,serial)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lwann
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate
    LOGICAL                                  :: serial

    INTEGER                                  :: lgmat, nn, nomax

    IF (serial.OR.parai%nproc.EQ.1) THEN
       lgmat=nstate*nstate
       lwann=lgmat+100
    ELSE
       CALL set_orbdist(nstate,1,parai%nproc,nomax)
       nn=MAX(8*nomax*nomax*wannc%nwanopt,2*nomax*nstate*wannc%nwanopt,&
            2*nomax*parai%nproc*wannc%nwanopt)
       lgmat=2*nn+2*nomax*nstate*wannc%nwanopt+4*nomax*nomax
       nn=MAX(4*nomax*nomax,nstate*nstate)
       lgmat=lgmat+nn+2*parai%nproc
       lwann=lgmat+100
    ENDIF
    tag='LWANN'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_jrot
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE xapyb(n,x,ra,na,y,rb,nb,z,lr)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: n, na
    COMPLEX(real_8)                          :: x(n,na)
    INTEGER                                  :: nb
    COMPLEX(real_8)                          :: y(n,nb), z(n,nb)
    INTEGER                                  :: lr
    REAL(real_8)                             :: ra(lr,nb), rb(lr,nb)

    INTEGER                                  :: i, j, k

! ==--------------------------------------------------------------==

    CALL zeroing(z)!,n*nb)
    !$omp parallel do private(J,I,K) __COLLAPSE3
    DO j=1,nb
       DO i=1,n
          DO k=1,na
             z(i,j)=z(i,j)+x(i,k)*ra(k,j)
          ENDDO
       ENDDO
    ENDDO
    !$omp parallel do private(J,I,K) __COLLAPSE3
    DO j=1,nb
       DO i=1,n
          DO k=1,nb
             z(i,j)=z(i,j)+y(i,k)*rb(k,j)
          ENDDO
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE xapyb
  ! ==================================================================
  SUBROUTINE rmupx(n,m,rmat,x,ldx,nx,y,ldy)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: n, m
    REAL(real_8)                             :: rmat(n,m)
    INTEGER                                  :: ldx, nx
    COMPLEX(real_8)                          :: x(ldx,nx)
    INTEGER                                  :: ldy
    COMPLEX(real_8)                          :: y(ldy,nx)

    INTEGER                                  :: i, j, k

! ==--------------------------------------------------------------==

    !$omp parallel do private(I,J,K) __COLLAPSE3
    DO i=1,n
       DO j=1,nx
          DO k=1,m
             y(i,j)=y(i,j)+rmat(i,k)*x(k,j)
          ENDDO
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rmupx
  ! ==================================================================
  SUBROUTINE music(top, bot, n)
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: n
    INTEGER, DIMENSION(1:n)                  :: bot, top

    INTEGER                                  :: k, m
    INTEGER, DIMENSION(n)                    :: tmp_bot, tmp_top

    m = n/2

    DO k=1, m
       IF (k.EQ.1) THEN
          tmp_top(1) = 1
       ELSE IF (k.EQ.2) THEN
          tmp_top(k) = bot(1)
       ELSE IF (k.GT.2) THEN
          tmp_top(k) = top(k-1)
       ENDIF

       IF (k.EQ.m) THEN
          tmp_bot(k) = top(k)
       ELSE
          tmp_bot(k) = bot(k+1)
       ENDIF
    ENDDO

    top(1:m) = tmp_top(1:m)
    bot(1:m) = tmp_bot(1:m)

  END SUBROUTINE music
  ! ==================================================================

END MODULE jrotation_utils
