MODULE wr30wfn_utils
  USE dotp_utils,                      ONLY: dotp
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: int_1,&
                                             int_2,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8
  USE summat_utils,                    ONLY: summat
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: dgive,&
                                             fskip,&
                                             setkwf
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  !!public :: wr30wfn
  !public :: wrwfn
  !!public :: rd30wfn
  !!public :: rdwfns
  !!public :: wrwfns
  !public :: rdwfn
  !!public :: write_n
  !public :: read_n
  !public :: pwtoao
  !public :: aotopw
  PUBLIC :: calcco
  !!public :: putconj

CONTAINS


  ! ==================================================================
  SUBROUTINE calcco(xxmat,xsmat,compl,nattot,nst)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nattot
    REAL(real_8)                             :: xsmat(nattot,nattot)
    INTEGER                                  :: nst
    REAL(real_8)                             :: xxmat(nattot,nst), compl(nst)

    INTEGER                                  :: i, k, l

    DO i=1,nst
       DO k=1,nattot
          DO l=1,nattot
             compl(i)=compl(i)+xxmat(k,i)*xsmat(k,l)*xxmat(l,i)
          ENDDO
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE calcco
  ! ==================================================================


END MODULE wr30wfn_utils



SUBROUTINE pwtoao(nw,ierror,nstate,nkpoint,c0,tau0,icompr,irecord)
  ! ==--------------------------------------------------------------==
  USE dotp_utils, ONLY: dotp
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE setbasis_utils, ONLY : loadc
  USE mp_interface, ONLY: mp_sum
  USE system , ONLY:cntl,maxsys,ncpw,nkpt,spar
  USE parac, ONLY : paral,parai
  USE atwf , ONLY:atwp,catom,loadc_foc_array_size,xsmat,xxmat
  USE ions , ONLY:ions0,ions1
  USE kpts , ONLY:tkpts
  USE phfac_utils, ONLY : phfac
  USE ovlap_utils, ONLY : ovlap
  USE wr30wfn_utils, ONLY : calcco
  USE summat_utils, ONLY : summat
  USE utils, ONLY : setkwf
  USE utils, ONLY : dgive
  USE zeroing_utils,                   ONLY: zeroing
  IMPLICIT NONE
  INTEGER                                    :: nw, ierror, nstate, nkpoint
  COMPLEX(real_8)                            :: c0(nkpt%ngwk,nstate)
  REAL(real_8)                               :: tau0(3,maxsys%nax,maxsys%nsx)
  INTEGER                                    :: icompr, irecord

  CHARACTER(*), PARAMETER                    :: procedureN = 'pwtoao'

  INTEGER                                    :: alloc, i, ia, iaorb, iat, &
                                                ierr, info, is, isca, isub, &
                                                ixx, lscr, natst, nx
  INTEGER(int_8)                             :: fpos_null
  INTEGER, ALLOCATABLE                       :: ipiv(:)
  INTEGER, ALLOCATABLE, DIMENSION(:)         :: icmp
  INTEGER, EXTERNAL                          :: idamax
  LOGICAL                                    :: tlsd2
  REAL(real_8)                               :: foc(loadc_foc_array_size), &
                                                scale, sfc
  REAL(real_8), ALLOCATABLE                  :: compl(:), scr(:)

  CALL tiset('    PWTOAO',isub)
  fpos_null=0
  IF (tkpts%tkpnt)CALL stopgm('WR30FN','K-POINTS NOT IMPLEMENTED',& 
       __LINE__,__FILE__)
  ! Initialize
  CALL phfac(tau0)
  alloc=MAX(atwp%nattot*atwp%nattot,atwp%nattot*nstate)
  ALLOCATE(xxmat(atwp%nattot,alloc/atwp%nattot),STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
       __LINE__,__FILE__)
  ALLOCATE(xsmat(atwp%nattot,alloc/atwp%nattot),STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
       __LINE__,__FILE__)
  CALL zeroing(xsmat)
  ALLOCATE(catom(nkpt%ngwk,atwp%nattot),STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
       __LINE__,__FILE__)
  CALL zeroing(catom)
  lscr=alloc
  ALLOCATE(scr(lscr),STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
       __LINE__,__FILE__)
  CALL zeroing(scr)
  IF (paral%parent) THEN
     alloc=MAX(atwp%nattot,nstate)
     ALLOCATE(ipiv(alloc),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     CALL zeroing(ipiv)
     ALLOCATE(compl(alloc),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     CALL zeroing(compl)
  ENDIF
  ! ==--------------------------------------------------------------==
  ! Calculate Atomic orbitals
  iaorb=1
  iat=0
  DO is=1,ions1%nsp
     DO ia=1,ions0%na(is)
        iat=iat+1
        CALL loadc(catom(1,iaorb),foc,nkpt%ngwk,ncpw%ngw,atwp%nattot-iaorb+1,&
             SIZE(foc),is,iat,natst)
        DO ixx=iaorb,iaorb+natst-1
           sfc=dotp(ncpw%ngw,catom(:,ixx),catom(:,ixx))
           CALL mp_sum(sfc,parai%allgrp)
           sfc=1._real_8/SQRT(sfc)
           CALL dscal(2*ncpw%ngw,sfc,catom(1,ixx),1)
        ENDDO
        iaorb=iaorb+natst
     ENDDO
  ENDDO
  ! Overlap matrices
  tlsd2=cntl%tlsd
  cntl%tlsd=.FALSE.
  ! <AO|AO>
  CALL ovlap(atwp%nattot,xsmat,catom,catom)
  cntl%tlsd=tlsd2
  CALL summat(xsmat,atwp%nattot)
  IF (tkpts%tkpnt) CALL setkwf(ncpw%ngw,atwp%nattot,catom)
  ! <AO|PSI>
  CALL ovlap2(ncpw%ngw,atwp%nattot,nstate,xxmat,catom,c0,.TRUE.)
  CALL mp_sum(xxmat,atwp%nattot*nstate,parai%allgrp)
  IF (paral%io_parent) THEN
     CALL dcopy(atwp%nattot*atwp%nattot,xsmat(1,1),1,scr(1),1)
     ! Calculate expansion coefficients
     CALL dgesv(atwp%nattot,nstate,xsmat,atwp%nattot,ipiv,xxmat,atwp%nattot,info)
     ! Write expansion coefficients to restart file
     irecord=5
     isca=idamax(atwp%nattot*nstate,xxmat,1)
     scale=ABS(dgive(xxmat,isca))
     WRITE(unit=nw,iostat=ierror) irecord,fpos_null
     WRITE(unit=nw,iostat=ierror) nstate,spar%ngwks,icompr,scale
     WRITE(unit=nw,iostat=ierror) atwp%nattot,nstate
     ALLOCATE(icmp(2*atwp%nattot*nstate),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
          __LINE__,__FILE__)
     CALL compress(atwp%nattot*nstate,xxmat,icmp,scale,ABS(icompr),nx)
     WRITE(unit=nw,iostat=ierror) nx
     CALL write_n(nw,ierror,icmp,nx)
     DEALLOCATE(icmp,STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
          __LINE__,__FILE__)
     ! Calculate completness of basis
     CALL calcco(xxmat,scr,compl,atwp%nattot,nstate)
     WRITE(6,'(/,A)') ' COMPLETNESS OF ATOMIC PROJECTION '
     WRITE(6,'(12F5.2)') (compl(i),i=1,nstate)
     DEALLOCATE(ipiv,STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
          __LINE__,__FILE__)
     DEALLOCATE(compl,STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
          __LINE__,__FILE__)
  ENDIF
  ! Deallocate memory
  DEALLOCATE(xxmat,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
       __LINE__,__FILE__)
  DEALLOCATE(xsmat,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
       __LINE__,__FILE__)
  DEALLOCATE(catom,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
       __LINE__,__FILE__)
  DEALLOCATE(scr,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
       __LINE__,__FILE__)
  CALL tihalt('    PWTOAO',isub)
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE pwtoao
! ==================================================================
SUBROUTINE aotopw(nr,c0,nstate,tau0,icompr,scale)
  ! ==--------------------------------------------------------------==
  USE dotp_utils, ONLY: dotp
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE mp_interface, ONLY: mp_sum,mp_bcast
  USE setbasis_utils, ONLY : loadc
  USE prng_utils, ONLY : repprngu_vec, repprngu, prngparaskip
  USE system , ONLY:cntl,maxsys,ncpw,nkpt
  USE parac, ONLY : paral,parai
  USE ions , ONLY:ions0,ions1
  USE atwf , ONLY:atwp,catom,loadc_foc_array_size,xxmat
  USE kpts , ONLY:tkpts
  USE spin , ONLY:spin_mod
  USE gsortho_utils, ONLY : gs_ortho, gs_ortho_c, gsortho
  USE phfac_utils, ONLY : phfac
  IMPLICIT NONE
  INTEGER                                    :: nr, nstate
  COMPLEX(real_8)                            :: c0(nkpt%ngwk,nstate)
  REAL(real_8)                               :: tau0(3,maxsys%nax,maxsys%nsx)
  INTEGER                                    :: icompr
  REAL(real_8)                               :: scale

  CHARACTER(*), PARAMETER                    :: procedureN = 'aotopw'

  INTEGER                                    :: alloc, ia, iaorb, iat, ierr, &
                                                is, isub, ixx, natom, natst, &
                                                nst, nx
  INTEGER, ALLOCATABLE, DIMENSION(:)         :: icmp, isma
  REAL(real_8)                               :: foc(loadc_foc_array_size), sfc
  REAL(real_8), ALLOCATABLE, DIMENSION(:)    :: smat

  CALL tiset('    AOTOPW',isub)
  IF (tkpts%tkpnt)CALL stopgm('WR30FN','K-POINTS NOT IMPLEMENTED',& 
       __LINE__,__FILE__)
  ! Initialize
  CALL phfac(tau0)
  alloc=MAX(atwp%nattot*atwp%nattot,atwp%nattot*nstate)
  ALLOCATE(xxmat(atwp%nattot,alloc/atwp%nattot),STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
       __LINE__,__FILE__)
  ALLOCATE(catom(nkpt%ngwk,atwp%nattot*ncpw%ngw/nkpt%ngwk),STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
       __LINE__,__FILE__)
  ! ==--------------------------------------------------------------==
  ! Calculate Atomic orbitals
  iaorb=1
  iat=0
  DO is=1,ions1%nsp
     DO ia=1,ions0%na(is)
        iat=iat+1
        CALL loadc(c0,foc,ncpw%ngw,ncpw%ngw,nstate,SIZE(foc),is,iat,natst)
        DO ixx=1,natst
           sfc=dotp(ncpw%ngw,c0(:,ixx),c0(:,ixx))
           CALL mp_sum(sfc,parai%allgrp)
           sfc=1._real_8/SQRT(sfc)
           CALL dscal(2*ncpw%ngw,sfc,c0(1,ixx),1)
        ENDDO
        CALL dcopy(2*ncpw%ngw*natst,c0(1,1),1,catom(1,iaorb),1)
        iaorb=iaorb+natst
     ENDDO
  ENDDO
  IF (paral%io_parent) THEN
     ! Read Restart file
     READ(nr) natom,nst
     IF (natom.NE.atwp%nattot) THEN
        WRITE(6,*) ' AOTOPW| ATOMIC ORBITAL BASIS HAS CHANGED'
        CALL stopgm('AOTOPW',' CHANGE OF BASIS',& 
             __LINE__,__FILE__)
     ENDIF
     IF (nstate.NE.nst) THEN
        WRITE(6,*) ' AOTOPW| NBR STATES HAS CHANGED'
        CALL stopgm('AOTOPW',' CHANGE OF BASIS',& 
             __LINE__,__FILE__)
     ENDIF
     CALL repprngu_vec(atwp%nattot*nst,xxmat)
     ALLOCATE(icmp(2*atwp%nattot*nst),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
          __LINE__,__FILE__)
     READ(nr) nx
     CALL read_n(nr,icmp,nx)
     CALL decompr(atwp%nattot*nst,xxmat,icmp,scale,ABS(icompr),nx)
     DEALLOCATE(icmp,STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
          __LINE__,__FILE__)
  ENDIF
  CALL mp_bcast(xxmat,SIZE(xxmat),parai%io_source,parai%cp_grp)
  ! Calculate MO  
  IF (ncpw%ngw.GT.0)&
       CALL dgemm('N','N',2*ncpw%ngw,nstate,atwp%nattot,1.0_real_8,catom(1,1),&
       2*ncpw%ngw,xxmat(1,1),atwp%nattot,0.0_real_8,c0(1,1),2*ncpw%ngw)
  ALLOCATE(smat(nstate**2),isma(nstate),STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
       __LINE__,__FILE__)
  IF (cntl%tlsd) THEN
     CALL gs_ortho(c0,0,c0(:,1),spin_mod%nsup)
     CALL gs_ortho(c0,0,c0(:,spin_mod%nsup+1),spin_mod%nsdown)
  ELSE
     CALL gs_ortho(c0,0,c0,nstate)
  ENDIF
  DEALLOCATE(smat,isma,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
       __LINE__,__FILE__)
  DEALLOCATE(xxmat,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
       __LINE__,__FILE__)
  DEALLOCATE(catom,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
       __LINE__,__FILE__)
  CALL tihalt('    AOTOPW',isub)
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE aotopw
! ==================================================================

! ==================================================================
SUBROUTINE rdwfn(nr,c,ca,cr,icmp,mapw,ngwks0,icompr,scale,tkpnt0,&
     is_read_para)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE mp_interface, ONLY: mp_sync,mp_send, mp_recv
  USE system , ONLY:cntl,mapgp,ncpw,nkpt,parap
  USE parac, ONLY : paral,parai
  USE kpts , ONLY:tkpts
  USE io_utils, ONLY : file_close, file_open, file_seek, file_get_position, file_write_l, &
       file_write_char, file_write_i, file_write_d, file_write_z, file_read_l, &
       file_read_char, file_read_i, file_read_d, file_read_z
  USE zeroing_utils,                   ONLY: zeroing
  USE cp_ieee_interface, ONLY : cp_ieee_is_finite
  IMPLICIT NONE
  INTEGER                                    :: nr
  COMPLEX(real_8)                            :: c(nkpt%ngwk), cr(2*nkpt%ngwk)
  INTEGER                                    :: icmp(*), mapw(*), ngwks0
  COMPLEX(real_8)                            :: ca(ngwks0)
  INTEGER                                    :: icompr
  REAL(real_8)                               :: scale
  LOGICAL                                    :: tkpnt0, is_read_para

  CHARACTER(*), PARAMETER                    :: procedureN = 'rdwfn'

  INTEGER                                    :: i, ig, ip, msgid, msglen, &
                                                ngwip, ngws0, nx, nx_vec(1)
  LOGICAL                                    :: my_io_parent, use_mpi

  IF (is_read_para.OR.paral%cp_inter_io_parent) THEN
     my_io_parent=paral%io_parent
     IF (is_read_para) THEN
        my_io_parent=paral%parent
     ENDIF
     use_mpi=is_read_para.AND.cntl%use_mpi_io
     IF (tkpnt0.AND.tkpts%tkpnt) THEN
        ngws0=ngwks0/2
     ELSE
        ngws0=ngwks0
     ENDIF
     IF (my_io_parent) THEN
        CALL file_read_i(nr,1,nx_vec,use_mpi)
        nx=nx_vec(1)
        CALL file_read_i(nr,2*nx,icmp,use_mpi)
        CALL decompr(2*ngwks0,ca,icmp,scale,icompr,nx)
        CALL zeroing(c)!,nkpt%ngwk)
        CALL zeroing(cr)!,2*nkpt%ngwk)
     ENDIF
     DO ip=0,parai%nproc-1
        IF (my_io_parent) THEN
           IF (parap%pgroup(ip+1).EQ.parai%me) THEN
              !$omp  parallel do schedule(static,1000) &
              !$omp              private(IG) &
              !$omp              shared(NGWS0)
              DO ig=1,parap%sparm(3,ip)
                 IF (mapgp(ig).LE.ngws0) c(ig)=ca(mapgp(ig))
              ENDDO
              !$omp end parallel do
              IF (tkpnt0.AND.tkpts%tkpnt) THEN
                 ngwip=parap%sparm(3,ip)
                 !$omp  parallel do schedule(static,1000) &
                 !$omp              private(IG) &
                 !$omp              shared(NGWS0,NGWIP)
                 DO ig=1,ngwip
                    IF (mapgp(ig).LE.ngws0)&
                         c(ngwip+ig)=ca(ngws0+mapgp(ig))
                 ENDDO
                 !$omp end parallel do
              ENDIF
           ELSE
              msgid=1
              CALL mp_recv(mapw,parap%sparm(3,ip),parap%pgroup(ip+1),msgid,parai%allgrp)
              !$omp  parallel do schedule(static,1000) &
              !$omp              private(IG) &
              !$omp              shared(NGWKS0)
              DO ig=1,parap%sparm(3,ip)
                 IF (mapw(ig).LE.ngwks0) cr(ig)=ca(mapw(ig))
              ENDDO
              !$omp end parallel do
              IF (tkpnt0.AND.tkpts%tkpnt) THEN
                 ngwip=parap%sparm(3,ip)
                 !$omp  parallel do schedule(static,1000) &
                 !$omp              private(IG) &
                 !$omp              shared(NGWS0,NGWIP)
                 DO ig=1,ngwip
                    IF (mapw(ig).LE.ngws0) cr(ngwip+ig)&
                         =ca(ngws0+mapw(ig))
                 ENDDO
                 !$omp end parallel do
              ENDIF
              msgid=2
              IF (tkpnt0.AND.tkpts%tkpnt) THEN
                 msglen = 2 * parap%sparm(3,ip)
              ELSE
                 msglen =     parap%sparm(3,ip)
              ENDIF
              CALL mp_send(cr,msglen,parap%pgroup(ip+1),msgid,parai%allgrp)
           ENDIF
        ELSE
           IF (parap%pgroup(ip+1).EQ.parai%me) THEN
              msgid=1
              CALL mp_send(mapgp,ncpw%ngw,parap%pgroup(1),msgid,parai%allgrp)
              msgid=2
              !msglen = 2 * nkpt%ngwk * 8
              CALL mp_recv(c,nkpt%ngwk,parap%pgroup(1),msgid,parai%allgrp)
           ENDIF
        ENDIF
        CALL mp_sync(parai%allgrp)
     ENDDO
  ENDIF
  ! ==--------------------------------------------------------------==
  DO i = 1, SIZE(c)
     IF( .NOT. cp_ieee_is_finite( c( i ) ) ) THEN
        CALL stopgm(procedureN,'value is not finite',& 
             __LINE__,__FILE__)
     ENDIF
  ENDDO
  RETURN
20 CONTINUE
  CALL stopgm(procedureN,'END OF FILE',& 
       __LINE__,__FILE__)
30 CONTINUE
  CALL stopgm(procedureN,'READ ERROR',& 
       __LINE__,__FILE__)
END SUBROUTINE rdwfn
SUBROUTINE wrwfn(nw,ierror,c,ca,cr,icmp,mapw,icompr,scale,&
     is_read_para)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE mp_interface, ONLY: mp_sync,mp_send, mp_recv
  USE system , ONLY:cntl,mapgp,ncpw,nkpt,spar,parap
  USE parac, ONLY : paral,parai
  USE kpts , ONLY:tkpts
  USE io_utils, ONLY : file_close, file_open, file_seek, file_get_position, file_write_l, &
       file_write_char, file_write_i, file_write_d, file_write_z, file_read_l, &
       file_read_char, file_read_i, file_read_d, file_read_z
  USE gsortho_utils, ONLY : gs_ortho, gs_ortho_c, gs_ortho
  IMPLICIT NONE
  INTEGER                                    :: nw, ierror
  COMPLEX(real_8)                            :: c(nkpt%ngwk), ca(spar%ngwks), &
                                                cr(*)
  INTEGER                                    :: icmp(*), mapw(*), icompr
  REAL(real_8)                               :: scale
  LOGICAL                                    :: is_read_para

  INTEGER                                    :: ig, ip, ipw, msgid, msglen, &
                                                ngwip, nx, nx_vec(1)
  LOGICAL                                    :: my_io_parent, use_mpi

! Variables
! ==--------------------------------------------------------------==
! Need to protect if we cpmd groups

  IF (is_read_para.OR.paral%cp_inter_io_parent) THEN
     my_io_parent=paral%io_parent
     IF (is_read_para) THEN
        my_io_parent=paral%parent
     ENDIF
     use_mpi=is_read_para.AND.cntl%use_mpi_io

     DO ip=0,parai%nproc-1
        IF (my_io_parent) THEN
           IF (parap%pgroup(ip+1).EQ.parai%me) THEN
              !$omp  parallel do schedule(static,1000) &
              !$omp              private(IG)
              DO ig=1,ncpw%ngw
                 cr(ig)=c(ig)
                 mapw(ig)=mapgp(ig)
              ENDDO
              !$omp end parallel do
              IF (tkpts%tkpnt) THEN
                 !$omp  parallel do schedule(static,1000) &
                 !$omp              private(IG)
                 DO ig=ncpw%ngw+1,nkpt%ngwk
                    cr(ig)=c(ig)
                 ENDDO
                 !$omp end parallel do
              ENDIF
           ELSE
              msgid=1
              IF (tkpts%tkpnt) THEN
                 msglen = 2 * parap%sparm(3,ip)
              ELSE
                 msglen =     parap%sparm(3,ip)
              ENDIF
              CALL mp_recv(cr,msglen,parap%pgroup(ip+1),msgid,parai%allgrp)
              msgid=2
              CALL mp_recv(mapw,parap%sparm(3,ip),parap%pgroup(ip+1),msgid,parai%allgrp)
           ENDIF
           !$omp  parallel do schedule(static,1000) &
           !$omp              private(IPW)
           DO ipw=1,parap%sparm(3,ip)
              ca(mapw(ipw))=cr(ipw)
           ENDDO
           !$omp end parallel do
           IF (tkpts%tkpnt) THEN
              ngwip=parap%sparm(3,ip)
              !$omp  parallel do schedule(static,1000) &
              !$omp              private(IPW) &
              !$omp              shared(NGWIP)
              DO ipw=1,ngwip
                 ca(mapw(ipw)+spar%ngws)=cr(ngwip+ipw)
              ENDDO
              !$omp end parallel do
           ENDIF
        ELSE
           IF (parap%pgroup(ip+1).EQ.parai%me) THEN
              msgid=1
              !msglen = 2 * nkpt%ngwk * 8
              CALL mp_send(c,nkpt%ngwk,parap%pgroup(1),msgid,parai%allgrp)
              msgid=2
              CALL mp_send(mapgp,ncpw%ngw,parap%pgroup(1),msgid,parai%allgrp)
           ENDIF
        ENDIF
     ENDDO
     IF (my_io_parent) THEN
        CALL compress(2*spar%ngwks,ca,icmp,scale,icompr,nx)
        ! WRITE(UNIT=NW,IOSTAT=IERROR) NX
        ! CALL WRITE_N(NW,IERROR,ICMP,NX)
        nx_vec(1)=nx
        CALL file_write_i(nw,1,nx_vec,use_mpi)
        CALL file_write_i(nw,2*nx,icmp,use_mpi)
     ENDIF
     CALL mp_sync(parai%allgrp)
  ENDIF
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE wrwfn
! ==================================================================


! ==================================================================
SUBROUTINE wr30wfn(nw,ierror,nstate,c,tau0,tag,irecord)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE mp_interface, ONLY: mp_max
  USE system , ONLY:cnti,maxsys,ncpw,nkpbl,nkpt,spar
  USE parac, ONLY : paral,parai
  USE kpts , ONLY:tkpts
  USE bsym , ONLY:bsfac
  USE utils, ONLY : dgive
  IMPLICIT NONE
  INTEGER                                    :: nw, ierror, nstate
  REAL(real_8) :: c(2*nkpt%ngwk,nstate*nkpt%nkpnt*bsfac), &
      tau0(3,maxsys%nax,maxsys%nsx)
  CHARACTER(len=*)                           :: tag
  INTEGER                                    :: irecord

  CHARACTER(*), PARAMETER                    :: procedureN = 'wr30wfn'

  INTEGER                                    :: i, icompr = 0, idamax, ierr, &
                                                ikpt, isca, isub, kbeg, kend, &
                                                kinc, len, len1, nkpoint
  INTEGER(int_8)                             :: fpos_null
  INTEGER, ALLOCATABLE                       :: icmp(:), mapw(:)
  REAL(real_8)                               :: scale, scale0
  REAL(real_8), ALLOCATABLE                  :: ca(:), cr(:)

  CALL tiset('   WR30WFN',isub)
  fpos_null=0
  IF (cnti%wcompb.LT.0) THEN
     icompr=cnti%wcompb
     CALL pwtoao(nw,ierror,nstate,nkpt%nkpnt,c,tau0,icompr,irecord)
  ELSE
     IF (paral%io_parent) THEN
        len=2*nkpt%ngwk
        len1=ncpw%ngw + 1
        ALLOCATE(cr(4*len),STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
             __LINE__,__FILE__)
        ALLOCATE(mapw(2*len1),STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
             __LINE__,__FILE__)
        ALLOCATE(ca(2*spar%ngwks),STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
             __LINE__,__FILE__)
        ALLOCATE(icmp(2*(2*spar%ngwks+100)),STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
             __LINE__,__FILE__)! allocate as int*8
        icompr=cnti%wcompb
     ELSE
        ! To avoid runtime error regarding 'not allocated' CA            
        ALLOCATE(ca(1),STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
             __LINE__,__FILE__)
        ALLOCATE(cr(1),STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
             __LINE__,__FILE__)
        ALLOCATE(icmp(1),STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
             __LINE__,__FILE__)
        ALLOCATE(mapw(1),STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
             __LINE__,__FILE__)
     ENDIF
     IF (icompr.LT.1) icompr=1
     IF (tag(1:MIN(3,LEN_TRIM(tag))).EQ.'NIL') THEN
        ! NN: For BS two sets of WF welocities 
        irecord=2*nstate*bsfac+1
        isca=idamax(2*nstate*nkpt%ngwk,c,1)
        scale=ABS(dgive(c(1,1),isca))
        CALL mp_max(scale,parai%allgrp)
        IF (paral%io_parent) THEN
           WRITE(unit=nw,iostat=ierror) irecord,fpos_null
           WRITE(unit=nw,iostat=ierror) nstate,spar%ngwks,icompr,scale
        ENDIF
        ! NN: For BS two sets of WF welocities are written 
        DO i=1,nstate*bsfac
           CALL wrwfn(nw,ierror,c(1,i),ca,cr,icmp,mapw,icompr,scale,&
                .FALSE.)
        ENDDO
     ELSE
        irecord=2*nstate*nkpt%nkpts+1
        scale=0._real_8
        CALL inq_swap(kbeg,kend,kinc)
        DO ikpt=kbeg,kend,kinc
           nkpoint=nkpbl(ikpt)
           IF (tkpts%tkblock) CALL rkpt_swap(c,nstate,ikpt,tag)
           isca=idamax(2*nkpt%ngwk*nstate*nkpoint,c,1)
           scale0=ABS(dgive(c(1,1),isca))
           IF (scale0.GT.scale) scale=scale0
        ENDDO
        CALL mp_max(scale,parai%allgrp)

        IF (paral%io_parent) THEN
           WRITE(unit=nw,iostat=ierror) irecord,fpos_null
           WRITE(unit=nw,iostat=ierror) nstate*nkpt%nkpts,spar%ngwks,icompr,scale
        ENDIF
        DO ikpt=1,nkpt%nblkp
           IF (tkpts%tkblock) CALL rkpt_swap(c,nstate,ikpt,tag)
           DO i=1,nstate*nkpbl(ikpt)
              CALL wrwfn(nw,ierror,c(1,i),ca,cr,icmp,mapw,icompr,scale,&
                   .FALSE.)
           ENDDO
        ENDDO
     ENDIF
     IF (paral%io_parent) THEN
        DEALLOCATE(ca,STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
             __LINE__,__FILE__)
        DEALLOCATE(icmp,STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
             __LINE__,__FILE__)
        DEALLOCATE(cr,STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
             __LINE__,__FILE__)
        DEALLOCATE(mapw,STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
             __LINE__,__FILE__)
     ENDIF
  ENDIF
  CALL tihalt('   WR30WFN',isub)
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE wr30wfn

! ==================================================================
SUBROUTINE rd30wfn(nr,c,nstate,tau0,info,&
     tkpnt0,nstate0,nkpts0,tag,fpos_end)
  ! ==--------------------------------------------------------------==
  ! == Read wavefunctions                                           ==
  ! == INPUT:                                                       ==
  ! ==   NR      Logical number of the file                         ==
  ! ==   NSTATE  Number of states                                   ==
  ! ==   TAU0    atomic positions                                   ==
  ! ==   TKPNT0  .TRUE. (has k points)                              ==
  ! ==   NSTATE0 Number of states  for the RESTART file             ==
  ! ==   NKPTS0  Number of kpoints for the RESTART file             ==
  ! ==   TAG     'NIL' wavefunctions                                ==
  ! ==           'C0'  electronic wavefunctions                     ==
  ! == OUTPUT:                                                      ==
  ! ==   C      Wavefunctions                                       ==
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE mp_interface, ONLY: mp_bcast
  USE prng_utils, ONLY : repprngu_vec, repprngu, prngparaskip, repprngu_vec_cmplx
  USE system , ONLY:cnti,cntl,maxsys,ncpw,nkpbl,nkpt,spar
  USE parac, ONLY : paral,parai
  USE kpts , ONLY:tkpts
  USE spin , ONLY:spin_mod
  USE geq0mod , ONLY:geq0
  USE bsym , ONLY:bsclcs
  USE io_utils, ONLY : file_close, file_open, file_seek, file_get_position, file_write_l, &
       file_write_char, file_write_i, file_write_d, file_write_z, file_read_l, &
       file_read_char, file_read_i, file_read_d, file_read_z
  USE kpclean_utils, ONLY : c_clean, r_clean
  USE gsortho_utils, ONLY : gs_ortho, gs_ortho_c
  USE randtowf_utils, ONLY : randtowf
  USE zeroing_utils,                   ONLY: zeroing
  USE utils, ONLY : fskip
  IMPLICIT NONE
  INTEGER                                    :: nr, nstate
  COMPLEX(real_8)                            :: c(nkpt%ngwk,nstate,nkpt%nkpnt)
  REAL(real_8)                               :: tau0(3,maxsys%nax,maxsys%nsx)
  INTEGER                                    :: info
  LOGICAL                                    :: tkpnt0
  INTEGER                                    :: nstate0, nkpts0
  CHARACTER(len=*)                           :: tag
  INTEGER(int_8)                             :: fpos_end

  CHARACTER(*), PARAMETER                    :: procedureN = 'rd30wfn'

  INTEGER :: i, icompr, ierr, ikind, ikpt, ikpt0, isub, lca, len, len1, &
      licmp, n0, ngwa, ngwks0, noortho, northo, nread, nx
  COMPLEX(real_8), ALLOCATABLE               :: smat(:)
  INTEGER, ALLOCATABLE                       :: icmp(:), isma(:), mapw(:)
  REAL(real_8)                               :: scale
  REAL(real_8), ALLOCATABLE                  :: ca(:), cr(:)

!real(real_8), target :: c(2*nkpt%ngwk,nstate,nkpt%nkpnt)
!complex(real_8), pointer :: c_complex(:,:,:)
! Variables
! ==--------------------------------------------------------------==

  CALL tiset(procedureN,isub)
  info=0
  IF (paral%io_parent) READ(nr) n0,ngwks0,icompr,scale
  CALL mp_bcast(n0,parai%io_source,parai%cp_grp)
  CALL mp_bcast(ngwks0,parai%io_source,parai%cp_grp)
  CALL mp_bcast(icompr,parai%io_source,parai%cp_grp)
  cnti%rcompb=icompr
  IF (cnti%rcompb.LT.0) THEN
     CALL aotopw(nr,c,nstate,tau0,icompr,scale)
     GOTO 400
  ENDIF
  ! ==--------------------------------------------------------------==
  IF (paral%io_parent) THEN
     len=2*nkpt%ngwk
     ALLOCATE(cr(4*len),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     len1=2*ncpw%ngw + 1
     ALLOCATE(mapw(len1),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     ngwa=MAX(ngwks0,spar%ngwks)
     lca=MAX(2*ngwa+100,nstate0*MAX(nstate0,nstate))
     ALLOCATE(ca(lca),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     CALL zeroing(ca)!,lca)
     ! Warning: ICMP is allocated as INTEGER*8
     licmp=MAX(2*ngwa+100,(nstate-nstate0)+1)
     ALLOCATE(icmp(2*licmp),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)! ICMP must be allocated as I*8
     ! ICMP is integer*8 (length LICMP).
     CALL zeroing(icmp)!,2*licmp)
  ELSE
     ! To avoid Fortran runtime error 'not allocated'
     ALLOCATE(cr(1),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     ALLOCATE(mapw(1),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     ALLOCATE(ca(1),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     ALLOCATE(icmp(1),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)! ICMP must be allocated as I*8
  ENDIF
  ! ==--------------------------------------------------------------==
  IF (TRIM(tag).EQ.'NIL') THEN
     nread=MIN(nstate,n0)
     ! NN: Modification necessary for two set of WF-Velocities 
     IF (.NOT.cntl%bsymm)THEN
        DO i=1,nread
           CALL rdwfn(nr,c(1,i,1),ca,cr,icmp,mapw,ngwks0,icompr,&
                scale,tkpnt0,.FALSE.)
        ENDDO
     ELSE
        ! NN: TAKING CARE OF HS_STATE_WF_VELOCITIES
        IF (bsclcs.EQ.2)THEN
           ! NN: SKIP BS_WF_VELOCITIES
           IF (paral%io_parent)CALL fskip(nr,2*n0)
           ! NN: READ HS_WF_VELOCITIES
           DO i=1,nread
              CALL rdwfn(nr,c(1,i,1),ca,cr,icmp,mapw,ngwks0,icompr,&
                   scale,tkpnt0,.FALSE.)
           ENDDO
           ! NN: TAKING CARE OF BS_STATE_WF_VELOCITIES 
        ELSE IF (bsclcs.EQ.1)THEN
           DO i=1,nread
              CALL rdwfn(nr,c(1,i,1),ca,cr,icmp,mapw,ngwks0,icompr,&
                   scale,tkpnt0,.FALSE.)
           ENDDO
           ! NN: SKIP HS_WF_VELOCITIES
           IF (paral%io_parent)CALL fskip(nr,2*n0)
        ENDIF
     ENDIF
     nx=nstate-nread
     IF(nx.GT.0) CALL repprngu_vec_cmplx(nkpt%ngwk*nx,c(1,nread+1,1))
     GOTO 200
  ENDIF
  ! ==--------------------------------------------------------------==
  ! == Case for C0 (electronic wavefunctions)                       ==
  ! ==--------------------------------------------------------------==
  ! Check consistency of N0.
  IF (n0.NE.nstate0*nkpts0.AND.paral%io_parent) THEN
     WRITE(6,*) procedureN//'! N0=',n0,' NSTATE0=',nstate0,&
          'NKPST0=',nkpts0
     CALL stopgm(procedureN,&
          'THE NUMBER OF WAVEFUNCTIONS IS NOT CORRECT',& 
          __LINE__,__FILE__)
  ENDIF
  nread=0
  ikpt0=0
  DO ikpt=1,nkpt%nblkp
     DO ikind=1,nkpbl(ikpt)
        ikpt0=ikpt0+1
        IF (ikpt0.LE.nkpts0) THEN
           ! We load a new set of nstates for a new kpoint.
           noortho=0
           DO i=1,nstate
              IF (i.LE.nstate0) THEN
                 ! We load a new state.
                 nread=nread+1
                 CALL rdwfn(nr,c(1,i,ikind),ca,cr,icmp,mapw,ngwks0,&
                      icompr,scale,tkpnt0,.FALSE.)
                 ! We adapt the format of C0.
                 IF (tkpts%tkpnt) THEN
                    IF (.NOT.tkpnt0) THEN
                       ! C(NGW:2*NGW,I,IKIND)=CONJG(C(1:NGW,I,IKIND))
                       CALL dcopy(2*ncpw%ngw,c(1      ,i,ikind),1,&
                            c(1+2*ncpw%ngw,i,ikind),1)
                       CALL dscal(ncpw%ngw,-1._real_8,c(2+2*ncpw%ngw,i,ikind),2)
                       IF (geq0) THEN
                          c(2*ncpw%ngw+1,i,ikind)=0._real_8
                          c(2*ncpw%ngw+2,i,ikind)=0._real_8
                       ENDIF
                    ENDIF
                    !call reshape_inplace(c, (/nkpt%ngwk,nstate,nkpt%nkpnt/), c_complex)
                    CALL c_clean(c(:,i,ikind),1,ikind)
                 ENDIF
              ELSE
                 ! Use a randomized state.
                 !call reshape_inplace(c, (/nkpt%ngwk,nstate,nkpt%nkpnt/), c_complex)
                 CALL randtowf(c(1:,i:,ikind),1,ikind,ikpt0)
              ENDIF
           ENDDO
           ! C has to be real.
           IF (tkpnt0.AND..NOT.tkpts%tkpnt) THEN
              IF (cntl%tlsd) THEN
                 CALL putconj(c,spin_mod%nsup)
                 CALL putconj(c(1,spin_mod%nsup+1,1),spin_mod%nsdown)
              ELSE
                 CALL putconj(c,nstate)
              ENDIF
           ENDIF
           IF (nstate.GT.nstate0) THEN
              noortho=nstate-nstate0
              ! We orthogonalize.
              IF (tkpts%tkpnt) THEN
                 ALLOCATE(smat(nstate*nstate),STAT=ierr)
                 IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                      __LINE__,__FILE__)
                 ALLOCATE(isma(nstate),STAT=ierr)
                 IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                      __LINE__,__FILE__)
                 IF (cntl%tlsd) THEN
                    IF (spin_mod%nsup.GT.nstate0) THEN
                       CALL gs_ortho_c(c(:,:,ikind),nstate0,&
                            c(:,nstate0+1,ikind),spin_mod%nsup,smat)
                       ! We copy UP states  to DOWN states.
                       CALL dcopy(2*nkpt%ngwk*spin_mod%nsdown,c(1,1,ikind),1,&
                            c(1,spin_mod%nsup+1,ikind),1)
                    ELSEIF (spin_mod%nsup.EQ.nstate0) THEN
                       ! We copy UP states  to DOWN states.
                       CALL dcopy(2*nkpt%ngwk*spin_mod%nsdown,c(1,1,ikind),1,&
                            c(1,spin_mod%nsup+1,ikind),1)
                    ELSE
                       ! We orthogonalize only the down states.
                       northo=nstate0-spin_mod%nsup
                       CALL gs_ortho_c(c(:,spin_mod%nsup+1:,ikind),northo,&
                            c(:,spin_mod%nsup+northo+1,ikind),noortho,&
                            smat)
                    ENDIF
                 ELSE
                    CALL gs_ortho_c(c(:,:,ikind),nstate0,&
                         c(:,nstate0+1,ikind),noortho,&
                         SMAT)
                 ENDIF
              ELSE
                 ALLOCATE(smat(nstate*nstate),STAT=ierr)
                 IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                      __LINE__,__FILE__)
                 ALLOCATE(isma(nstate),STAT=ierr)
                 IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                      __LINE__,__FILE__)
                 IF (cntl%tlsd) THEN
                    IF (spin_mod%nsup.GT.nstate0) THEN
                       CALL gs_ortho  (c(:,:,ikind),nstate0,&
                            c(:,nstate0+1,ikind),spin_mod%nsup)
                       ! We copy UP states  to DOWN states.
                       CALL dcopy(2*nkpt%ngwk*spin_mod%nsdown,c(1,1,ikind),1,&
                            c(1,spin_mod%nsup+1,ikind),1)
                    ELSEIF (spin_mod%nsup.EQ.nstate0) THEN
                       ! We copy UP states  to DOWN states.
                       CALL dcopy(2*nkpt%ngwk*spin_mod%nsdown,c(1,1,ikind),1,&
                            c(1,spin_mod%nsup+1,ikind),1)
                    ELSE
                       ! We orthogonalize only the down states.
                       northo=nstate0-spin_mod%nsup
                       CALL gs_ortho  (c(:,spin_mod%nsup+1:,ikind),northo,&
                            c(:,spin_mod%nsup+northo+1,ikind),noortho)
                    ENDIF
                 ELSE
                    CALL gs_ortho  (c(:,:,ikind),nstate0,&
                         c(:,nstate0+1,ikind),noortho)
                 ENDIF
              ENDIF     ! IF(TKPNT)
              DEALLOCATE(smat,STAT=ierr)
              IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
                   __LINE__,__FILE__)
              DEALLOCATE(isma,STAT=ierr)
              IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
                   __LINE__,__FILE__)
           ELSEIF (paral%io_parent) THEN
              ! NSTATE0.GE.NSTATE
              nx=nstate0-nstate
              DO i=1,nx
                 nread=nread+1
                 READ(nr)
                 READ(nr)
              ENDDO
           ENDIF         ! IF(NSTATE.GT.NSTATE0)
        ELSE
           ! We have more k points -> we use the set of nstate
           ! for ikind=1.
           IF (ikind.NE.1) THEN
              CALL dcopy(2*nkpt%ngwk*nstate,c(1,1,1),1,c(1,1,ikind),1)
           ELSE
              ! All the block has to be copied. Do nothing.
              GOTO 100
           ENDIF
        ENDIF             ! IF(IKPT0.LE.KPTS0)
     ENDDO                 ! DO IKIND
100  CONTINUE
     IF (tkpts%tkblock) CALL wkpt_swap(c,nstate,ikpt,tag)
  ENDDO                     ! DO IKPT
  ! ==--------------------------------------------------------------==
200 CONTINUE
  IF (n0.GT.nread.AND.paral%io_parent) THEN
     nx=n0-nread
     IF (cntl%is_in_stream) THEN
        CALL file_seek(nr,.FALSE.,fpos_end)
     ELSE
        CALL fskip(nr,2*nx)
     ENDIF
  ENDIF
  IF (paral%io_parent) THEN
     DEALLOCATE(ca,STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
          __LINE__,__FILE__)
     DEALLOCATE(icmp,STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
          __LINE__,__FILE__)
     DEALLOCATE(cr,STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
          __LINE__,__FILE__)
     DEALLOCATE(mapw,STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
          __LINE__,__FILE__)
  ENDIF
  ! ==--------------------------------------------------------------==
400 CONTINUE
  CALL tihalt(procedureN,isub)
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE rd30wfn

! ==================================================================
SUBROUTINE write_n(nw,ierror,x,n)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  IMPLICIT NONE
  INTEGER                                    :: nw, ierror, n
  REAL(real_8)                               :: x(n)

! ==--------------------------------------------------------------==

  WRITE(unit=nw,iostat=ierror) x
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE write_n
! ==================================================================
SUBROUTINE read_n(nr,x,n)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  IMPLICIT NONE
  INTEGER                                    :: nr, n
  REAL(real_8)                               :: x(n)

! ==--------------------------------------------------------------==

  READ(nr,END=20,err=30) x
  ! ==--------------------------------------------------------------==
  RETURN
20 CONTINUE
  CALL stopgm('READ_N','END OF FILE',& 
       __LINE__,__FILE__)
30 CONTINUE
  CALL stopgm('READ_N','READ ERROR',& 
       __LINE__,__FILE__)
END SUBROUTINE read_n

SUBROUTINE putconj(c0,nstate)
  ! ==--------------------------------------------------------------==
  ! == MULTIPLY BY A COMPLEX TO HAVE C0(1) REAL                     ==
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE mp_interface, ONLY: mp_bcast
  USE system , ONLY:ncpw
  USE parac, ONLY : paral,parai
  USE geq0mod , ONLY:geq0
  USE rgs_utils, ONLY : rgs,rgs_c
  IMPLICIT NONE
  INTEGER                                    :: nstate
  COMPLEX(real_8)                            :: c0(ncpw%ngw,nstate)

  CHARACTER(*), PARAMETER                    :: procedureN = 'putconj'
  REAL(real_8), PARAMETER                    :: eps = 1.e-10_real_8 

  COMPLEX(real_8)                            :: zz
  INTEGER                                    :: i, ierr, nr
  REAL(real_8)                               :: zi
  REAL(real_8), ALLOCATABLE                  :: smat(:)

  nr=0
  DO i=1,nstate
     IF (geq0) THEN
        zz=c0(1,i)
        zi=ABS(AIMAG(c0(1,i)))
     ENDIF
     CALL mp_bcast(zz,parai%igeq0,parai%allgrp)
     CALL mp_bcast(zi,parai%igeq0,parai%allgrp)
     IF (ABS(zz).GT.eps.AND.zi.GT.eps) THEN
        zz=ABS(zz)/zz
        CALL zscal(ncpw%ngw,zz,c0(1,1),1)
        nr=nr+1
     ENDIF
  ENDDO
  IF (nr.GT.0) THEN
     ALLOCATE(smat(nstate*nstate),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     CALL rgs(c0,nstate,smat)
     DEALLOCATE(smat,STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
          __LINE__,__FILE__)
  ENDIF
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE putconj
! ==================================================================

! ==================================================================
SUBROUTINE wrwfns(nw,nstate_to_write,ierror,c,ca,cr,icmp,mapw,&
     icompr,scale)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE mp_interface, ONLY: mp_bcast
  USE system , ONLY:cntl,nkpt
  USE parac, ONLY : paral,parai
  USE part_1d, ONLY: part_1d_get_el_in_blk
  USE part_1d, ONLY: part_1d_nbr_el_in_blk
  USE io_utils, ONLY : file_close, file_open, file_seek, file_get_position, file_write_l, &
       file_write_char, file_write_i, file_write_d, file_write_z, file_read_l, &
       file_read_char, file_read_i, file_read_d, file_read_z
  IMPLICIT NONE
  INTEGER                                    :: nw, nstate_to_write, ierror
  COMPLEX(real_8)                            :: c(nkpt%ngwk,*), ca(*), &
                                                cr(2*nkpt%ngwk)
  INTEGER                                    :: icmp(*), mapw(*), icompr
  REAL(real_8)                               :: scale

  CHARACTER(*), PARAMETER                    :: procedureN = 'wrwfns'

  CHARACTER(len=256)                         :: file_name
  INTEGER                                    :: i, is, isub, my_me, my_nproc, &
                                                nw_para
  INTEGER(int_8)                             :: fpos_cur, fpos_wfn_beg, &
                                                fpos_wfn_end, wfn_size
  LOGICAL                                    :: my_io_parent, &
                                                write_wfn_in_para

  CALL tiset(procedureN,isub)
  ! ==--------------------------------------------------------------==              

  write_wfn_in_para=cntl%is_out_stream.AND.parai%cp_nogrp>1.AND.cntl%use_mpi_io

  fpos_wfn_beg=0;fpos_wfn_end=0;wfn_size=0
  my_me=0;my_nproc=1;my_io_parent=paral%io_parent
  nw_para=nw
  IF (write_wfn_in_para) THEN
     ! get the size of 1 state
     IF (paral%io_parent) CALL file_get_position(nw,.FALSE.,fpos_wfn_beg)

     IF (nstate_to_write>0) THEN
        CALL wrwfn(nw,ierror,c,ca,cr,icmp,mapw,icompr,scale,.FALSE.)
     ENDIF

     IF (paral%io_parent) CALL file_get_position(nw,.FALSE.,&
          wfn_size)
     wfn_size=wfn_size-fpos_wfn_beg
     fpos_wfn_end=fpos_wfn_beg+nstate_to_write*wfn_size

     CALL mp_bcast(fpos_wfn_beg,parai%io_source,parai%cp_grp)
     CALL mp_bcast(fpos_wfn_end,parai%io_source,parai%cp_grp)
     CALL mp_bcast(wfn_size,parai%io_source,parai%cp_grp)
     ! go back at the first state
     IF (paral%io_parent) CALL file_seek(nw,.FALSE.,fpos_wfn_beg)
     ! which group writes what
     my_me=parai%cp_inter_me
     my_nproc=parai%cp_nogrp
     my_io_parent=paral%parent
     ! get file infos
     IF (paral%io_parent) INQUIRE(unit=nw,name=file_name)
     CALL mp_bcast(file_name,parai%io_source,parai%cp_grp)

     ! close file
     IF (paral%io_parent) CALL file_close(nw,.FALSE.)

     ! open para files
     IF (cntl%use_mpi_io) THEN
        ! start at 0 (so -1)
        fpos_wfn_beg=fpos_wfn_beg-1
        IF (my_io_parent) CALL file_open(file_name,'X','X',&
             'READWRITE',.TRUE.,.TRUE.,parai%cp_inter_grp,nw_para)
     ELSE
        IF (.NOT.paral%io_parent.AND.my_io_parent) THEN
           nw_para=666
           CALL file_open(file_name,'UNFORMATTED','OLD',&
                'READWRITE',.TRUE.,.FALSE.,0,nw_para)
        ENDIF
     ENDIF
  ENDIF

  ! ==--------------------------------------------------------------==

  DO is = 1,part_1d_nbr_el_in_blk(nstate_to_write,my_me,my_nproc)
     i = part_1d_get_el_in_blk(is,nstate_to_write,my_me,my_nproc)
     IF (my_io_parent.AND.write_wfn_in_para) THEN
        fpos_cur=(i-1)*wfn_size+fpos_wfn_beg
        CALL file_seek(nw_para,cntl%use_mpi_io,fpos_cur)
     ENDIF

     ! We write a new state.
     CALL wrwfn(nw_para,ierror,c(1,i),ca,cr,icmp,mapw,icompr,scale,&
          write_wfn_in_para)
  ENDDO

  ! ==--------------------------------------------------------------==
  IF (write_wfn_in_para) THEN
     IF (cntl%use_mpi_io) THEN
        IF (my_io_parent) CALL file_close(nw_para,.TRUE.)
     ELSE
        CALL stopgm(procedureN,'we shouldnt get there',& 
             __LINE__,__FILE__)
        IF (.NOT.paral%io_parent.AND.my_io_parent) THEN
           CALL file_close(nw_para,.FALSE.)
        ENDIF
     ENDIF

     ! reopen file and point to the end
     IF (paral%io_parent) THEN
        CALL file_open(file_name,'UNFORMATTED','OLD','READWRITE',&
             cntl%is_out_stream,.FALSE.,0,nw)
        CALL file_seek(nw,.FALSE.,fpos_wfn_end)
     ENDIF


  ENDIF

  ! ==--------------------------------------------------------------==
  CALL tihalt(procedureN,isub)
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE wrwfns

! ==================================================================
SUBROUTINE rdwfns(nr,nstate_to_read,c,ca,cr,icmp,mapw,ngwks0,&
     icompr,scale,tkpnt0,nread)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE mp_interface, ONLY: mp_bcast,mp_sum
  USE system , ONLY:cntl,nkpt
  USE parac, ONLY : paral,parai
  USE part_1d, ONLY: part_1d_get_el_in_blk
  USE part_1d, ONLY: part_1d_nbr_el_in_blk
  USE io_utils, ONLY : file_close, file_open, file_seek, file_get_position, file_write_l, &
       file_write_char, file_write_i, file_write_d, file_write_z, file_read_l, &
       file_read_char, file_read_i, file_read_d, file_read_z
  USE zeroing_utils,                   ONLY: zeroing
  IMPLICIT NONE
  INTEGER                                    :: nr, nstate_to_read
  COMPLEX(real_8)                            :: c(nkpt%ngwk,*), &
                                                cr(2*nkpt%ngwk)
  INTEGER                                    :: icmp(*), mapw(*), ngwks0
  COMPLEX(real_8)                            :: ca(ngwks0)
  INTEGER                                    :: icompr
  REAL(real_8)                               :: scale
  LOGICAL                                    :: tkpnt0
  INTEGER                                    :: nread

  CHARACTER(*), PARAMETER                    :: procedureN = 'rdwfns'

  CHARACTER(len=256)                         :: file_name
  INTEGER                                    :: i, is, isub, my_me, my_nproc, &
                                                nr_para, NREAD_local
  INTEGER(int_8)                             :: fpos_cur, fpos_wfn_beg, &
                                                fpos_wfn_end, wfn_size
  LOGICAL                                    :: my_io_parent, read_wfn_in_para

  CALL tiset(procedureN,isub)
  ! ==--------------------------------------------------------------==              

  NREAD_local=0

  ! ==--------------------------------------------------------------==              
  read_wfn_in_para=cntl%is_in_stream.AND.parai%cp_nogrp>1
  FPOS_WFN_BEG=0;FPOS_WFN_END=0;WFN_SIZE=0
  my_me=0;my_nproc=1;my_io_parent=paral%io_parent
  nr_para=nr
  IF (read_wfn_in_para) THEN
     ! get the size of 1 state
     IF (paral%io_parent) CALL file_get_position(nr,.FALSE.,&
          fpos_wfn_beg)
     IF (nstate_to_read>0) THEN
        CALL rdwfn(nr,c,ca,cr,icmp,mapw,ngwks0,&
             icompr,scale,tkpnt0,.FALSE.)
     ENDIF
     IF (paral%io_parent) CALL file_get_position(nr,.FALSE.,&
          wfn_size)
     wfn_size=wfn_size-fpos_wfn_beg
     fpos_wfn_end=fpos_wfn_beg+nstate_to_read*wfn_size
     CALL mp_bcast(fpos_wfn_beg,parai%io_source,parai%cp_grp)
     CALL mp_bcast(fpos_wfn_end,parai%io_source,parai%cp_grp)
     CALL mp_bcast(wfn_size,parai%io_source,parai%cp_grp)
     ! go back at the first state
     IF (paral%io_parent) CALL file_seek(nr,.FALSE.,fpos_wfn_beg)
     ! which group reads what
     my_me=parai%cp_inter_me
     my_nproc=parai%cp_nogrp
     my_io_parent=paral%parent
     ! open para files
     IF (paral%io_parent) INQUIRE(unit=nr,name=file_name)
     CALL mp_bcast(file_name,parai%io_source,parai%cp_grp)
     IF (cntl%use_mpi_io) THEN
        ! start at 0 (so -1)
        fpos_wfn_beg=fpos_wfn_beg-1
        IF (my_io_parent) CALL file_open(file_name,'X','X',&
             'READ',.TRUE.,.TRUE.,parai%cp_inter_grp,nr_para)
     ELSE
        IF (.NOT.paral%io_parent.AND.my_io_parent) THEN
           nr_para=666
           CALL file_open(file_name,'UNFORMATTED','OLD',&
                'READ',.TRUE.,.FALSE.,0,nr_para)
        ENDIF
     ENDIF
     ! need to zero the state that we can allreduce
     CALL zeroing(c(:,1:nstate_to_read))!,nkpt%ngwk*nstate_to_read)
  ENDIF
  ! ==--------------------------------------------------------------==
  DO is = 1,part_1d_nbr_el_in_blk(nstate_to_read,my_me,my_nproc)
     i = part_1d_get_el_in_blk(is,nstate_to_read,my_me,my_nproc)
     IF (my_io_parent.AND.read_wfn_in_para) THEN
        fpos_cur=(i-1)*wfn_size+fpos_wfn_beg
        CALL file_seek(nr_para,cntl%use_mpi_io,fpos_cur)
     ENDIF
     ! We load a new state.
     NREAD_local=nread_local+1
     CALL rdwfn(nr_para,c(1,i),ca,cr,icmp,mapw,ngwks0,&
          icompr,scale,tkpnt0,read_wfn_in_para)
  ENDDO
  ! ==--------------------------------------------------------------==
  IF (read_wfn_in_para) THEN
     IF (cntl%use_mpi_io) THEN
        IF (my_io_parent) CALL file_close(nr_para,.TRUE.)
     ELSE
        IF (.NOT.paral%io_parent.AND.my_io_parent) THEN
           CALL stopgm(procedureN,'we shouldnt get there',& 
                __LINE__,__FILE__)
           CALL file_close(nr_para,.FALSE.)
        ENDIF
     ENDIF
     IF (paral%io_parent) CALL file_seek(nr,.FALSE.,fpos_wfn_end)
     CALL mp_sum(c,nkpt%ngwk*nstate_to_read,parai%cp_inter_grp)
     CALL mp_sum(NREAD_local,parai%cp_inter_grp)
  ENDIF
  ! ==--------------------------------------------------------------==

  nread=nread+NREAD_local

  ! ==--------------------------------------------------------------==
  CALL tihalt(procedureN,isub)
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE rdwfns
