#include "cpmd_global.h"

MODULE cdft_utils
  USE adat,                            ONLY: elem
  USE atomc_utils,                     ONLY: positx,&
                                             putrho
  USE atwf,                            ONLY: &
       atchg, atrg, atwf_mod, atwfr, atwp, atwr, cdftoc, cdftocl, m1shlx
  USE cdftmod,                         ONLY: &
       atchg2, cdftci, cdftcom, cdfthda, cdfthess, cdftlog, cdftmd, cdftpi, &
       cdftpred, cdftrun, cdftvgrs, chgset, czones, finalchrg, maxsavev, &
       rhol, sccomm, wa, wcfix, wd, wder, wdiff, wg_n, wg_sigma, wgaussl
  USE cnst,                            ONLY: fbohr,&
                                             fpi,&
                                             pi
  USE coor,                            ONLY: tau0
  USE cppt,                            ONLY: nzh
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE fft_maxfft,                      ONLY: maxfftn
  USE fftmain_utils,                   ONLY: fwfftn
  USE fitpack_utils,                   ONLY: curv1,&
                                             curv2,&
                                             curvd
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE isos,                            ONLY: isos1,&
                                             isos3
  USE kinds,                           ONLY: real_8
  USE metr,                            ONLY: metr_com
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_recv,&
                                             mp_send,&
                                             mp_sum
  USE parac,                           ONLY: parai,&
                                             paral
  USE pbc_utils,                       ONLY: pbc
  USE pslo,                            ONLY: pslo_com
  USE qspl,                            ONLY: nsplpo
  USE rhoofr_utils,                    ONLY: rhoofr
  USE rnlsm_utils,                     ONLY: rnlsm
  USE spin,                            ONLY: clsd
  USE system,                          ONLY: &
       cntl, cntr, fpar, iatpt, maxsp, maxsys, ncpw, nkpt, parap, parm, spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing
#include "sizeof.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: cdft_w
  !public :: cdft_chg
  !public :: cdft_int
  PUBLIC :: init_cdft
  PUBLIC :: cdft_adarrays
  PUBLIC :: cdft_reinit
  PUBLIC :: vupdate
  !public :: update_zone
  PUBLIC :: cdft_finalize
  PUBLIC :: write_w
  PUBLIC :: wcdft_restart
  !public :: rcdft_restart
  !public :: cdft_pri
  !public :: deriv_init
  PUBLIC :: cdft_forces
  !public :: cdft_fint
  !public :: atdens_new2
  !public :: ringbuff
  PUBLIC :: v_pred
  !public :: surfe

CONTAINS

  ! TODO: convert all text in write statements to uppercase
  ! ==================================================================
  ! cntl%cdft main source file containing all cntl%md and WF-OPT related
  ! functions H. Oberhofer (ho246@cam.ac.uk) 2009 OMP & optimizations -
  ! MB & TI - Strasbourg/Hyogo 2012
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE cdft_w(rhoe,tau0,chrg)
    ! ==--------------------------------------------------------------==
    ! == Computes weight functions (W) and Charges (CHRG)            ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rhoe(fpar%nnr1,*), &
                                                tau0(:,:,:), chrg(2)

    CHARACTER(*), PARAMETER                  :: procedureN = 'cdft_w'

    INTEGER                                  :: ia, ierr, is, it, mmax
    REAL(real_8)                             :: chgt, chrg1, chrg2, &
                                                datom(nsplpo,2), rcut, &
                                                wn(fpar%nnr1)
    REAL(real_8), ALLOCATABLE                :: arho(:), work(:)

! Variables
! Allocation of local variables

    ALLOCATE(arho(maxsys%mmaxx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(work(maxsys%mmaxx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    IF (.NOT.cdftlog%tspinc)THEN
       CALL dcopy(fpar%nnr1,rhoe,1,rhol,1)
    ELSE
       !$omp parallel do private(is)
       DO is=1,fpar%nnr1
          rhol(is)=rhoe(is,1)-2.0_real_8*rhoe(is,2)
       ENDDO
       !$omp end parallel do
    ENDIF
    ! ..Recalculate weights and charges
    IF (cdftlog%recw) THEN
       IF (paral%io_parent)&
            WRITE(6,*) "RECALCULATING WEIGHTS"
       CALL zeroing(wa)!,nnr1)
       IF (.NOT.wgaussl%thdas)CALL zeroing(wd)!,nnr1)
       CALL zeroing(wn)!,nnr1)
       chrg(1)=0.0_real_8
       chrg(2)=0.0_real_8
       ! ..Get Normalisation WN
       DO is=1,ions1%nsp
          CALL atdens_new2(is,datom,mmax,rcut,arho,work)
          DO ia=1,ions0%na(is)
             CALL putrho(wn,tau0(:,ia,is),&
                  atrg(1,is),datom,mmax,rcut,1._real_8)
          ENDDO
       ENDDO
       ! ..ACCEPTORS
       DO it=1,cdftpi%naccr
          CALL atdens_new2(cdftpi%spa(it),datom,mmax,rcut,arho,work)
          cdftpi%mmx_a=mmax
          CALL cdft_chg(rhol,wa,wn,tau0(:,cdftpi%sna(it),cdftpi%spa(it)),&
               atrg(1,cdftpi%spa(it)),datom,mmax,rcut,chgt)
          chrg(1)=chrg(1)+chgt
       ENDDO
       ! ..DONORS
       DO it=1,cdftpi%ndon
          CALL atdens_new2(cdftpi%spd(it),datom,mmax,rcut,arho,work)
          cdftpi%mmx_d=mmax
          CALL cdft_chg(rhol,wd,wn,tau0(:,cdftpi%snd(it),cdftpi%spd(it)),&
               atrg(1,cdftpi%spd(it)),datom,mmax,rcut,chgt)
          chrg(2)=chrg(2)+chgt
       ENDDO
       IF (.NOT.wgaussl%thdas)THEN
          !$omp parallel do private(is)
          DO is=1,fpar%nnr1
             wdiff(is)=wa(is)-wd(is)
          ENDDO
          !$omp end parallel do
       ELSE
          !$omp parallel do private(is)
          DO is=1,fpar%nnr1
             wdiff(is)=wa(is)
          ENDDO
          !$omp end parallel do
       ENDIF
       cdftlog%recw=.FALSE.
       cdftlog%newgrad=.FALSE.
    ELSE
       ! ..only calculate the charges
       CALL cdft_int(rhol,wa,cdftpi%mmx_a,chrg(1))
       IF (cdftpi%ndon.GT.0)THEN
          CALL cdft_int(rhol,wd,cdftpi%mmx_d,chrg(2))
       ENDIF
    ENDIF
    CALL mp_sum(chrg,2,parai%allgrp)
    chrg(1)=-chrg(1)*parm%omega/(spar%nr1s*spar%nr2s*spar%nr3s)
    IF (cdftpi%ndon.GT.0)THEN
       chrg(2)=-chrg(2)*parm%omega/(spar%nr1s*spar%nr2s*spar%nr3s)
    ELSE
       chrg(2)=0
    ENDIF
    IF (.NOT.cdftlog%tspinc)THEN
       chrg1=chrg(1)
       chrg2=chrg(2)
       !$omp parallel do private(it) reduction(+:CHRG1)
       DO it=1,cdftpi%naccr
          ! !        CHRG(1)=CHRG(1)+ZV(CDFTPI%SPA(IT))
          chrg1=chrg1+ions0%zv(cdftpi%spa(it))
       ENDDO
       !$omp end parallel do
       !$omp parallel do private(it) reduction(+:CHRG2)
       DO it=1,cdftpi%ndon
          ! !        CHRG(2)=CHRG(2)+ZV(CDFTPI%SPD(IT))
          chrg2=chrg2+ions0%zv(cdftpi%spd(it))
       ENDDO
       !$omp end parallel do
       chrg(1)=chrg1
       chrg(2)=chrg2
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Deallocation of local variables
    DEALLOCATE(work,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(arho,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE cdft_w
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE cdft_chg(rhoe,weight,psi,tau,rg,datom,mmax,rcut,chg)
    ! ==--------------------------------------------------------------==
    ! == Computes weight functions (W) and Charges (CHG)              ==
    ! ==--------------------------------------------------------------==
    REAL(real_8) :: rhoe(fpar%kr1,fpar%kr2s,*), weight(fpar%kr1,fpar%kr2s,*), &
      psi(fpar%kr1,fpar%kr2s,*), tau(3), rg(*), datom(nsplpo,2)
    INTEGER                                  :: mmax
    REAL(real_8)                             :: rcut, chg

    INTEGER                                  :: i, i1, igrid, ii, iorgx, &
                                                iorgy, iorgz, ix, j, j1, jj, &
                                                k, k1, kk, mr(3), nptau(3)
    LOGICAL                                  :: tadd
    REAL(real_8)                             :: da1(3,3), r, sc1(3), sc2(3), &
                                                wcur, x, x1, xptau(3), y, y1, &
                                                z, z1

! VARIABLES
! ==--------------------------------------------------------------==
! ..Atom position

    x1=tau(1)
    y1=tau(2)
    z1=tau(3)
    CALL pbc(x1,y1,z1,sc1(1),sc1(2),sc1(3),1,parm%apbc,parm%ibrav)
    CALL dgemv('T',3,3,1._real_8,metr_com%htm1,3,sc1,1,0._real_8,sc2,1)
    IF (sc2(1).LT.0._real_8) sc2(1)=sc2(1)+1.0_real_8
    IF (sc2(2).LT.0._real_8) sc2(2)=sc2(2)+1.0_real_8
    IF (sc2(3).LT.0._real_8) sc2(3)=sc2(3)+1.0_real_8
    nptau(1)=NINT(sc2(1)*spar%nr1s)+1
    nptau(2)=NINT(sc2(2)*spar%nr2s)+1
    nptau(3)=NINT(sc2(3)*spar%nr3s)+1
    sc1(1)=sc2(1)-REAL(nptau(1)-1,kind=real_8)/REAL(spar%nr1s,kind=real_8)
    sc1(2)=sc2(2)-REAL(nptau(2)-1,kind=real_8)/REAL(spar%nr2s,kind=real_8)
    sc1(3)=sc2(3)-REAL(nptau(3)-1,kind=real_8)/REAL(spar%nr3s,kind=real_8)
    IF (nptau(1).GT.spar%nr1s) nptau(1)=1
    IF (nptau(2).GT.spar%nr2s) nptau(2)=1
    IF (nptau(3).GT.spar%nr3s) nptau(3)=1
    CALL dgemv('T',3,3,1._real_8,metr_com%ht,3,sc1,1,0._real_8,xptau,1)
    DO i=1,3
       da1(i,1)=parm%a1(i)/REAL(spar%nr1s,kind=real_8)
       da1(i,2)=parm%a2(i)/REAL(spar%nr2s,kind=real_8)
       da1(i,3)=parm%a3(i)/REAL(spar%nr3s,kind=real_8)
    ENDDO
    ! ..Interaction region
    sc1(1)=rcut
    sc1(2)=rcut
    sc1(3)=rcut
    CALL dgemv('T',3,3,1._real_8,metr_com%htm1,3,sc1,1,0._real_8,sc2,1)
    mr(1)=NINT(sc2(1)*spar%nr1s)+1
    mr(2)=NINT(sc2(2)*spar%nr2s)+1
    mr(3)=NINT(sc2(3)*spar%nr3s)+1
    mr(1)=MIN(mr(1),spar%nr1s/2-1)
    mr(2)=MIN(mr(2),spar%nr2s/2-1)
    mr(3)=MIN(mr(3),spar%nr3s/2-1)
    ! ..Origin of density grid
    iorgx=nptau(1)-mr(1)
    iorgy=nptau(2)-mr(2)
    iorgz=nptau(3)-mr(3)
    IF (iorgx.LT.1) iorgx=iorgx+spar%nr1s
    IF (iorgy.LT.1) iorgy=iorgy+spar%nr2s
    IF (iorgz.LT.1) iorgz=iorgz+spar%nr3s
    ! ..loop over all grid points
    chg=0.0_real_8
    igrid=0
    DO i=0,2*mr(1)
       ix=iorgx+i
       IF (ix.GT.spar%nr1s) ix=ix-spar%nr1s
       IF (ix.GE.parap%nrxpl(parai%mepos,1).AND.ix.LE.parap%nrxpl(parai%mepos,2)) THEN
          ii=ix-parap%nrxpl(parai%mepos,1)+1
          i1=mr(1)-i
          DO j=0,2*mr(2)
             jj=iorgy+j
             IF (jj.GT.spar%nr2s) jj=jj-spar%nr2s
             j1=mr(2)-j
             DO k=0,2*mr(3)
                kk=iorgz+k
                IF (kk.GT.spar%nr3s) kk=kk-spar%nr3s
                k1=mr(3)-k
                CALL positx(x,y,z,i1,j1,k1,da1,xptau(1))
                IF (isos1%tclust.AND..NOT.isos1%ttwod)THEN
                   IF (x1-x.LT.parm%a1(1).AND.y1-y.LT.parm%a2(2).AND.z1-z.LT.parm%a3(3)&
                        .AND.x1-x.GT.0.AND.y1-y.GT.0.AND.z1-z.GT.0)THEN
                      r=SQRT(x*x+y*y+z*z)
                      IF (r.LT.rcut) THEN
                         igrid=igrid+1
                         wcur=curv2(r,mmax,rg,datom(1,1),datom(1,2),0._real_8)/&
                              psi(ii,jj,kk)
                         chg=chg+rhoe(ii,jj,kk)*wcur
                         weight(ii,jj,kk)=weight(ii,jj,kk)+wcur
                      ENDIF
                   ENDIF
                ELSEIF (isos1%tclust.AND.isos1%ttwod)THEN
                   tadd=.FALSE.
                   IF (isos3%snormal.EQ.1) THEN
                      IF (x1-x.LT.parm%a1(1).AND.x1-x.GT.0) tadd=.TRUE.
                   ELSEIF (isos3%snormal.EQ.2) THEN
                      IF (y1-y.LT.parm%a2(2).AND.y1-y.GT.0) tadd=.TRUE.
                   ELSEIF (isos3%snormal.EQ.3)THEN
                      IF (z1-z.LT.parm%a3(3).AND.z1-z.GT.0) tadd=.TRUE.
                   ENDIF
                   IF (tadd) THEN
                      r=SQRT(x*x+y*y+z*z)
                      IF (r.LT.rcut) THEN
                         igrid=igrid+1
                         wcur=curv2(r,mmax,rg,datom(1,1),datom(1,2),0._real_8)/&
                              psi(ii,jj,kk)
                         chg=chg+rhoe(ii,jj,kk)*wcur
                         weight(ii,jj,kk)=weight(ii,jj,kk)+wcur
                      ENDIF
                   ENDIF
                ELSE
                   r=SQRT(x*x+y*y+z*z)
                   IF (r.LT.rcut) THEN
                      igrid=igrid+1
                      wcur=curv2(r,mmax,rg,datom(1,1),datom(1,2),0._real_8)/&
                           psi(ii,jj,kk)
                      chg=chg+rhoe(ii,jj,kk)*wcur
                      weight(ii,jj,kk)=weight(ii,jj,kk)+wcur
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDIF
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE cdft_chg
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE cdft_int(rhoe,weight,mmax,chg)
    ! ==--------------------------------------------------------------==
    ! == Computes only Charges (CHG) by integrating RHO*WEIGHT        ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rhoe(fpar%kr1,fpar%kr2s,*), &
                                                weight(fpar%kr1,fpar%kr2s,*)
    INTEGER                                  :: mmax
    REAL(real_8)                             :: chg

    INTEGER                                  :: igrid, ii, jj, kk

    chg=0.0_real_8
    ! !    IGRID=0
    !$omp parallel do private(II,JJ,KK) reduction(+:CHG) __COLLAPSE3
    DO ii=1,fpar%kr1
       DO jj=1,fpar%kr2s
          DO kk=1,fpar%kr3s
             ! !              IGRID=IGRID+1
             chg=chg+rhoe(ii,jj,kk)*weight(ii,jj,kk)
          ENDDO
       ENDDO
    ENDDO
    !$omp end parallel do
    igrid=fpar%kr1*fpar%kr2s*fpar%kr3s
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE cdft_int
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE init_cdft
    ! ==--------------------------------------------------------------==
    ! == Initialises everything necessary for a cntl%cdft calculation      ==
    ! ==--------------------------------------------------------------==
    CHARACTER(*), PARAMETER                  :: procedureN = 'init_cdft'

    INTEGER                                  :: elc, ierr, is, ish, l, mmax
    LOGICAL, SAVE                            :: cfirst = .TRUE.
    REAL(real_8)                             :: datom(nsplpo,2), rcut
    REAL(real_8), ALLOCATABLE                :: arho(:), work(:)

    COMMON/rbuff/elc

    ! Allocation of local variables
    ALLOCATE(arho(maxsys%mmaxx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(work(maxsys%mmaxx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    cdftvgrs%cdft_v2=cdftcom%cdft_v
    cdfthess%chess(1,1)=1.0_real_8
    cdfthess%chess(1,2)=0.0_real_8
    cdfthess%chess(2,1)=0.0_real_8
    cdfthess%chess(2,2)=1.0_real_8

    IF (cdftlog%tcall)THEN
       IF (cntl%cdft_dekk)THEN
          CALL stopgm('CDFT',&
               "ERROR DUAL CONSTRAINT NOT IMPLEMENTED FOR DEKKER MINIMISATION",& 
               __LINE__,__FILE__)
       ENDIF
       IF (cntl%md)THEN
          CALL stopgm('CDFT',"ERROR DUAL CONSTRAINT cntl%md NOT IMPLEMENTED",& 
               __LINE__,__FILE__)
       ENDIF
       IF (cntl%geopt)THEN
          CALL stopgm('CDFT',"ERROR DUAL CONSTRAINT GEO-OPT NOT IMPLEMENTED",& 
               __LINE__,__FILE__)
       ENDIF
    ENDIF
    IF (cfirst)THEN
       ALLOCATE(rhol(fpar%nnr1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL mp_bcast_byte(cdftpi,size_in_bytes_of(cdftpi),parai%io_source,parai%cp_grp)
       IF (cdftpi%naccr.EQ.0.AND.paral%io_parent) THEN
          CALL stopgm('CDFT',"ERROR ACCEPTOR ATOMS NOT SET",& 
               __LINE__,__FILE__)
       ENDIF
       IF (cntl%md.AND.cdftlog%tspinc) THEN
          CALL stopgm('CDFT',&
               "ERROR SPIN DENSITY CONSTRAINT cntl%md NOT IMPLEMENTED",& 
               __LINE__,__FILE__)
       ENDIF
       IF (cntl%geopt.AND.cdftlog%tspinc) THEN
          CALL stopgm('CDFT',&
               "ERROR SPIN DENSITY CONSTRAINT GEO-OPT NOT IMPLEMENTED",& 
               __LINE__,__FILE__)
       ENDIF
       CALL mp_bcast_byte(wgaussl,size_in_bytes_of(wgaussl),parai%io_source,parai%cp_grp)
       CALL mp_bcast_byte(wcfix%wcut,size_in_bytes_of(wcfix),parai%io_source,parai%cp_grp)
       CALL mp_bcast(cdftvgrs%vgstep,parai%io_source,parai%cp_grp)
       IF (.NOT.cntl%bohr)wcfix%wcut=wcfix%wcut*fbohr
       IF (wgaussl%twgauss)THEN
          IF (wg_n.EQ.1.AND.ions1%nsp.NE.1)THEN
             DO is=2,ions1%nsp
                wg_sigma(is)=wg_sigma(1)
             ENDDO
          ENDIF
          IF (.NOT.cntl%bohr)THEN
             DO is=1,ions1%nsp
                wg_sigma(is)=wg_sigma(is)*fbohr
             ENDDO
          ENDIF
          CALL mp_bcast(wg_sigma,ions1%nsp,parai%io_source,parai%cp_grp)
       ENDIF
       IF (cdftpi%ndon.EQ.0.AND.paral%parent)THEN
          cdftcom%cdft_nc=-1.0_real_8*cdftcom%cdft_nc
          cdftcom%nother=-1.0_real_8*cdftcom%nother
       ENDIF
       atwp%nattot=0
       cdftlog%newgrad=.FALSE.
       cdftlog%dekk=.FALSE.
       DO is=1,ions1%nsp
          atwf_mod%numaor(is)=0
          DO ish=1,atwf_mod%nshell(is)
             l=atwf_mod%lshell(ish,is)
             atwp%nattot=atwp%nattot+ions0%na(is)*(2*l+1)
             atwf_mod%numaor(is)=atwf_mod%numaor(is)+(2*l+1)
          ENDDO
       ENDDO

       CALL cdft_adarrays

       ALLOCATE(wa(fpar%nnr1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(wd(fpar%nnr1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(wdiff(fpar%nnr1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(wdiff)!, nnr1)! avoid uninitialized warning

       IF (cdftlog%tczones)THEN
          IF (cdftcom%cdft_v(1).NE.0.0_real_8)THEN
             cntr%tolog=czones(3,2)
          ELSE
             cntr%tolog=czones(1,2)
          ENDIF
       ENDIF
       CALL mp_bcast_byte(cdftpi,size_in_bytes_of(cdftpi),parai%io_source,parai%cp_grp)
       cfirst=.FALSE.


       IF (cdftlog%thda)cdfthda%hdafirst=.TRUE.

       IF (cdftlog%rcdft)CALL rcdft_restart()
    ENDIF                     ! CFIRST

    cdftlog%recw=.TRUE.
    cdftrun%reconvv=.TRUE.

    IF (cntl%md.OR.cntl%geopt)CALL deriv_init()
    ! IF(.NOT.OCSET)THEN
    ! CALL DCOPY(M1SHLX*MAXSP,OC,1,OCUN,1)
    ! MSGLEN = M1SHLX*(2*MAXSP)*8
    ! CALL MY_BCAST(OCUN,MSGLEN,IO_SOURCE,CP_GRP)
    ! ENDIF

    DO is=1,ions1%nsp
       CALL atdens_new2(is,datom,mmax,rcut,arho,work)
       cdftmd%cdft_rc(is)=rcut
    ENDDO

    IF (paral%io_parent.AND.cntl%md.AND.cdftlog%tpred)THEN
       elc=0
       IF (cdftpred%predord+1.GT.maxsavev)THEN
          WRITE(6,*) 'ERROR IN CDFT PREDICTOR'
          WRITE(6,*) 'PREDICTOR ORDER TOO LARGE, RESETTING ORDER TO',&
               maxsaveV
          cdftpred%predord=maxsavev-1
       ENDIF
    ENDIF

    IF (paral%io_parent)CALL cdft_pri()
    IF (cntl%tsyscomb)CALL mp_bcast(sccomm%n_s0,parai%io_source,parai%cp_grp)

    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE init_cdft
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE cdft_adarrays()
    ! ==--------------------------------------------------------------==
    ! == (Re-)initialise the Donor and Acceptor arrays                ==
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: aindex, is, it, j

    DO it=1,cdftpi%ndon
       aindex=1
       DO is=1,ions1%nsp
          DO j=1,ions0%na(is)
             IF (aindex.EQ.cdftpi%cdft_d(it)) THEN
                cdftpi%spd(it)=is
                cdftpi%snd(it)=j
             ENDIF
             aindex=aindex+1
          ENDDO
       ENDDO
    ENDDO
    DO it=1,cdftpi%naccr
       aindex=1
       DO is=1,ions1%nsp
          DO j=1,ions0%na(is)
             IF (aindex.EQ.cdftpi%cdft_a(it)) THEN
                cdftpi%spa(it)=is
                cdftpi%sna(it)=j
             ENDIF
             aindex=aindex+1
          ENDDO
       ENDDO
    ENDDO

    IF (wgaussl%thdawm)THEN
       DO it=1+cdftpi%ndon,2*cdftpi%ndon
          aindex=1
          DO is=1,ions1%nsp
             DO j=1,ions0%na(is)
                IF (aindex.EQ.cdftpi%cdft_d(it)) THEN
                   cdftpi%spd(it)=is
                   cdftpi%snd(it)=j
                ENDIF
                aindex=aindex+1
             ENDDO
          ENDDO
       ENDDO
       DO it=1+cdftpi%naccr,2*cdftpi%naccr
          aindex=1
          DO is=1,ions1%nsp
             DO j=1,ions0%na(is)
                IF (aindex.EQ.cdftpi%cdft_a(it)) THEN
                   cdftpi%spa(it)=is
                   cdftpi%sna(it)=j
                ENDIF
                aindex=aindex+1
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE cdft_adarrays
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE cdft_reinit()
    ! ==--------------------------------------------------------------==
    ! == Re-initialise the optimisers after an cntl%md oder GeomOpt step   ==
    ! ==--------------------------------------------------------------==

    cdftlog%recw=.TRUE.
    cdftlog%newgrad=.FALSE.
    cdftlog%dekk=.FALSE.
    IF (cdftlog%thda)THEN
       IF (cdftocl%ocset)THEN
          IF (paral%io_parent)&
               WRITE(6,*)"RESETTING INITIAL OCCUPATION NUMBER"
          CALL dcopy(m1shlx*maxsp,cdftoc%oc2,1,atwf_mod%oc,1)
       ENDIF
       IF (chgset)THEN
          IF (paral%io_parent)WRITE(6,*)"Resetting initial charge! "
          CALL dcopy(maxsys%nsx,atchg2,1,atchg,1)
       ENDIF
       IF (cdftlog%rcdft)CALL rcdft_restart()
    ENDIF

    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE cdft_reinit
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE vupdate(c0,psi,rhoe,convtest)
    ! ==--------------------------------------------------------------==
    ! == Calculate the Gradient of W wrt V perform a V optimisation   ==
    ! == Step and check convergence                                   ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8) :: c0(nkpt%ngwk,crge%n,nkpt%nkpnt), psi(maxfftn)
    REAL(real_8)                             :: rhoe(fpar%nnr1,clsd%nlsd)
    LOGICAL                                  :: convtest

    INTEGER                                  :: isub
    REAL(real_8)                             :: aa, chrg_(2), chrs(2), &
                                                hnew(2,2), m, r1, s, sk(2), &
                                                tempv, vga, vsd, yk(2)

! ==--------------------------------------------------------------==

    CALL tiset('   VUPDATE',isub)
    ! ==--------------------------------------------------------------==

    CALL zeroing(chrg_)!,2)
    IF (pslo_com%tivan) CALL rnlsm(c0(:,:,1),crge%n,1,1,.FALSE.)
    CALL rhoofr(c0(:,:,1),rhoe,psi,crge%n)
    CALL mp_sum(ener_com%ekin,parai%allgrp)
    CALL cdft_w(rhoe,tau0,chrg_)
    cdftcom%vgrad(1)=chrg_(2)-chrg_(1)-cdftcom%cdft_nc
    IF (cdftlog%tcall)THEN
       cdftlog%tspinc=.TRUE.
       CALL cdft_w(rhoe,tau0,chrs)
       cdftlog%tspinc=.FALSE.
       cdftcom%vgrad(2)=chrs(2)-chrs(1)-cdftcom%cdft_ns
       ! CDFTLOG%NEWGRAD=.TRUE.
    ELSE
       cdftcom%vgrad(2)=0.0_real_8
    ENDIF
    vga=SQRT(cdftcom%vgrad(1)*cdftcom%vgrad(1)+cdftcom%vgrad(2)*cdftcom%vgrad(2))
    IF (paral%io_parent)THEN
       IF (cdftcom%vcconu.NE.cdftcom%vccon)THEN
          IF (ABS(cdftcom%vgrad(1)).GT.cdftcom%vcconu.OR.ABS(cdftcom%vgrad(1)).GT.cdftcom%vcconu) THEN
             cdftrun%reconvv=.TRUE.
          ELSE IF (.NOT.cdftrun%reconvv) THEN
             convtest=.TRUE.
             finalchrg(1)=chrg_(1)
             finalchrg(2)=chrg_(2)
             IF (.NOT.cdftlog%tcall)THEN
                WRITE(6,*) "VGRAD: ",cdftcom%vgrad(1),"V: ",cdftcom%cdft_v(1)
             ELSE
                WRITE(6,*) "VGRAD: ",cdftcom%vgrad
                WRITE(6,*) "V: ",cdftcom%cdft_v
             ENDIF
          ENDIF
       ELSE
          cdftrun%reconvv=.TRUE.
       ENDIF
       IF (ABS(cdftcom%vgrad(1)).LT.cdftcom%vccon.AND.ABS(cdftcom%vgrad(2)).LT.cdftcom%vccon) THEN
          convtest=.TRUE.
          finalchrg(1)=chrg_(1)
          finalchrg(2)=chrg_(2)
          IF (.NOT.cdftlog%tcall)THEN
             WRITE(6,*) "VGRAD: ",cdftcom%vgrad(1),"V: ",cdftcom%cdft_v(1)
          ELSE
             WRITE(6,*) "VGRAD: ",cdftcom%vgrad
             WRITE(6,*) "V: ",cdftcom%cdft_v
          ENDIF
          cdftrun%reconvv=.FALSE.
       ELSE IF (cdftrun%reconvv) THEN
          ! ..Dekkers Optimisation Method with Newton up to first sign change
          IF (.NOT.cdftlog%tcall)THEN
             IF (cdftlog%newgrad.AND.cntl%cdft_dekk) THEN
                IF (cdftlog%dekk)THEN
                   IF (cdftcom%vgrad(1)*cdftcom%vgrada.LT.0)THEN
                      cdftcom%vgrada=cdftcom%vgrado(1)
                      cdftcom%oldva=cdftcom%oldv(1)
                   ELSE
                      IF (cdftcom%vgrada.LT.vga)THEN
                         tempv=cdftcom%vgrada
                         cdftcom%vgrada=cdftcom%vgrad(1)
                         cdftcom%vgrad=tempv
                         tempv=cdftcom%oldva
                         cdftcom%oldva=cdftcom%cdft_v(1)
                         cdftcom%cdft_v=tempv
                      ENDIF
                   ENDIF
                ELSE
                   IF (cdftcom%vgrad(1)*cdftcom%vgrado(1).LT.0)THEN
                      cdftlog%dekk=.TRUE.
                      cdftcom%vgrada=cdftcom%vgrado(1)
                      cdftcom%oldva=cdftcom%oldv(1)
                   ENDIF
                ENDIF
                vsd=(cdftcom%vgrad(1)-cdftcom%vgrado(1))/(cdftcom%cdft_v(1)-cdftcom%oldv(1))
                cdftcom%oldv=cdftcom%cdft_v
                IF (cdftlog%dekk)THEN
                   m=(cdftcom%oldv(1)+cdftcom%oldva)/2.0_real_8
                   s=cdftcom%oldv(1)-cdftcom%vgrad(1)/vsd
                   IF ((cdftcom%vgrad(1).LT.s.AND.s.LT.m).OR.&
                        (m.LT.s.AND.s.LT.cdftcom%vgrad(1)))THEN
                      cdftcom%cdft_v(1)=s
                   ELSE
                      cdftcom%cdft_v(1)=m
                   ENDIF
                   cdftcom%vgrado=cdftcom%vgrad
                ELSE
                   cdftcom%cdft_v(1)=cdftcom%cdft_v(1)-cdftcom%vgrad(1)/vsd
                   cdftcom%vgrado=cdftcom%vgrad
                ENDIF
             ENDIF
             ! ..Newtons Optimisation Method
             IF (cdftlog%newgrad.AND.cntl%cdft_newt) THEN
                vsd=(cdftcom%vgrad(1)-cdftcom%vgrado(1))/(cdftcom%cdft_v(1)-cdftcom%oldv(1))
                IF (cdftvgrs%maxvmov*ABS(vsd).LT.ABS(cdftcom%vgrad(1))) THEN
                   vsd=SIGN(cdftvgrs%vgstep,vsd)
                ENDIF
                cdftcom%oldv=cdftcom%cdft_v
                cdftcom%cdft_v(1)=cdftcom%cdft_v(1)-cdftcom%vgrad(1)/vsd
                cdftcom%vgrado=cdftcom%vgrad
             ENDIF
          ELSEIF (cdftlog%newgrad.AND.cntl%cdft_newt) THEN
             ! CDFTCOM%OLDV=CDFTCOM%CDFT_V
             ! IF(CDFTCOM%CDFT_V(1).GT.2.0_real_8) THEN
             ! CDFTCOM%CDFT_V(2)=CDFTCOM%CDFT_V(2)+0.1_real_8
             ! CDFTCOM%CDFT_V(1)=CDFTCOM%CDFT_V(2)-0.2_real_8
             ! ELSE
             ! CDFTCOM%CDFT_V(1)=CDFTCOM%CDFT_V(1)+0.1_real_8
             ! ENDIF

             yk=cdftcom%vgrad-cdftcom%vgrado
             sk=cdftcom%cdft_v-cdftcom%oldv
             ! ..update the Hessian
             r1=yk(1)*sk(1)+yk(2)*sk(2)
             aa=cdfthess%chess(1,1)*sk(1)*sk(1)+sk(2)*(2.0_real_8*cdfthess%chess(1,2)*sk(1)&
                  +cdfthess%chess(2,2)*sk(2))
             hnew(1,1)=cdfthess%chess(1,1)+yk(1)*yk(1)/r1-((cdfthess%chess(1,1)*sk(1)&
                  +cdfthess%chess(1,2)*sk(2))**2)/aa
             hnew(1,2)=yk(1)*yk(2)/r1+(cdfthess%chess(1,2)*cdfthess%chess(1,2)&
                  -cdfthess%chess(1,1)*cdfthess%chess(2,2))*sk(1)*sk(2)/aa
             hnew(2,2)=cdfthess%chess(1,1)+yk(2)*yk(2)/r1-((cdfthess%chess(1,2)*sk(1)&
                  +cdfthess%chess(2,2)*sk(2))**2)/aa
             cdfthess%chess(1,1)=hnew(1,1)
             cdfthess%chess(1,2)=hnew(1,2)
             cdfthess%chess(2,1)=hnew(1,2)
             cdfthess%chess(2,2)=hnew(2,2)
             ! ..calculate inverse Hessian
             aa=cdfthess%chess(1,1)*cdfthess%chess(2,2)-cdfthess%chess(1,2)*cdfthess%chess(2,1)
             hnew(1,1)=cdfthess%chess(2,2)/aa
             hnew(1,2)=-cdfthess%chess(1,2)/aa
             hnew(2,1)=-cdfthess%chess(2,1)/aa
             hnew(2,2)=cdfthess%chess(1,1)/aa
             ! ..perform step
             sk(1)=-(hnew(1,1)*cdftcom%vgrad(1)+hnew(1,2)*cdftcom%vgrad(2))
             sk(2)=-(hnew(2,1)*cdftcom%vgrad(1)+hnew(2,2)*cdftcom%vgrad(2))
             cdftcom%oldv=cdftcom%cdft_v
             cdftcom%vgrado=cdftcom%vgrad
             cdftcom%cdft_v=cdftcom%cdft_v+sk/cdftvgrs%vgstep
          ENDIF
          ! ..Gradient Optimisation Method stepsize 0.1 (First step every Method to get a second point for the finite difference)
          IF ((.NOT.cntl%cdft_dekk.AND..NOT.cntl%cdft_newt).OR..NOT.cdftlog%newgrad)THEN
             cdftcom%oldv=cdftcom%cdft_v
             cdftcom%cdft_v=cdftcom%cdft_v+cdftcom%vgrad/cdftvgrs%vgstep
             cdftlog%newgrad=.TRUE.
             cdftcom%vgrado=cdftcom%vgrad
          ENDIF
          convtest=.FALSE.
          IF (cdftlog%tczones)CALL update_zone(cdftcom%vgrad(1)) !vw added array subscript
          IF (paral%io_parent) THEN
             IF (.NOT.cdftlog%tcall)THEN
                WRITE(6,*) "Old V: ",cdftcom%oldv(1),"New V: ",cdftcom%cdft_v(1)
                WRITE(6,*) "VGRAD: ",cdftcom%vgrad(1)
             ELSE
                WRITE(6,*) "Old V: ",cdftcom%oldv
                WRITE(6,*) "New V: ",cdftcom%cdft_v
                WRITE(6,*) "VGRAD: ",cdftcom%vgrad
             ENDIF
             IF (cdftpi%ndon.EQ.1.AND.cdftpi%naccr.EQ.1) THEN
                IF (.NOT.cdftlog%tcall)THEN
                   WRITE(6,'(1X,"Donor (",A,"):",F10.6,T30,'//&
                        '"Acceptor (",A,"):",F10.6)') elem%el(ions0%iatyp(cdftpi%spd(1))),chrg_(2),&
                        elem%el(ions0%iatyp(cdftpi%spa(1))),chrg_(1)
                ELSE
                   WRITE(6,'(1X,"Donor (",A,"):",F8.4,"(",F8.4,T30,")'//&
                        ' Acceptor (",A,"):",F8.4,"(",F8.4,")")') &
                        elem%el(ions0%iatyp(cdftpi%spd(1))),chrg_(2),&
                        chrs(2),elem%el(ions0%iatyp(cdftpi%spa(1))),chrg_(1),chrs(1)
                ENDIF
             ELSE IF (cdftpi%ndon.EQ.0.AND.cdftpi%naccr.EQ.1) THEN
                WRITE(6,*) "ACCEPTOR (",elem%el(ions0%iatyp(cdftpi%spa(1))),"):",chrg_(1)
             ELSE
                WRITE(6,'(1X,"DONOR:",F10.6,T30,"Acceptor:",F10.6)')&
                     chrg_(2),chrg_(1)
             ENDIF
          ENDIF
       ENDIF
    ENDIF
    CALL mp_bcast(cdftcom%cdft_v,SIZE(cdftcom%cdft_v),parai%io_source,parai%cp_grp)
    CALL mp_bcast(convtest,parai%io_source,parai%cp_grp)
    IF (cdftlog%tczones)CALL mp_bcast(cntr%tolog,parai%io_source,parai%cp_grp)
    CALL tihalt('   VUPDATE',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE vupdate
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE update_zone(grad)
    ! ==--------------------------------------------------------------==
    ! == Updates the Tolerance according to the new zone              ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: grad

    REAL(real_8)                             :: grada

    grada=ABS(grad)
    IF (grada<czones(3,1))THEN
       cntr%tolog=czones(3,2)
    ELSEIF (grada<czones(2,1))THEN
       cntr%tolog=czones(2,2)
    ELSE
       cntr%tolog=czones(1,2)
    ENDIF
    IF (paral%io_parent)WRITE(6,*)'TOLOG=',cntr%tolog

    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE update_zone
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE cdft_finalize()
    ! ==--------------------------------------------------------------==
    ! == Write out everything and delete the memory                   ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: it

    IF (paral%io_parent) THEN
       WRITE(6,'(1X,64("-"),/)')
       WRITE(6,'(1X,A)') '*** CDFT FINAL RESULTS ***'
       IF (.NOT.cdftlog%tcall) THEN
          WRITE(6,'(2X,A,T54,D10.4)') "FINAL V: ",cdftcom%oldv(1)
       ELSE
          WRITE(6,'(2X,A,T54,D10.4)') "FINAL V: ",cdftcom%oldv(1)
          WRITE(6,'(2X,T54,D10.4)') cdftcom%oldv(2)
       ENDIF
       WRITE(6,'(2X,A,T54,E10.2E1)') "VGRAD: ",finalchrg(2)-&
            finalchrg(1)-cdftcom%cdft_nc
       WRITE(6,'(2X,"DONOR:",F10.4,T30,"ACCEPTOR:",F10.4)')&
            finalchrg(2),finalchrg(1)
       IF (cdftpi%ndon.GT.0) WRITE(6,'(2X,"DONOR ATOMS")')
       DO it=1,cdftpi%ndon
          WRITE(6,'(3X,I4,"(",A,")")')cdftpi%cdft_d(it),elem%el(ions0%iatyp(cdftpi%spd(it)))
       ENDDO
       WRITE(6,'(2X,"ACCEPTOR ATOMS")')
       DO it=1,cdftpi%naccr
          WRITE(6,'(3X,I4,"(",A,")")')cdftpi%cdft_a(it),elem%el(ions0%iatyp(cdftpi%spa(it)))
       ENDDO
       WRITE(6,'(1X,64("-"),/)')
    ENDIF

    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE cdft_finalize
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE write_w(weight,suff)
    ! ==--------------------------------------------------------------==
    ! == Writes the weight to disk                                    ==
    ! ==--------------------------------------------------------------==
    ! Variables
    REAL(real_8)                             :: weight(fpar%kr1,fpar%kr2,*)
    CHARACTER(len=*)                         :: suff

    CHARACTER(*), PARAMETER                  :: procedureN = 'write_w'

    COMPLEX(real_8), ALLOCATABLE             :: w2(:), wtemp(:)
    INTEGER                                  :: i, ierr, il, ip, ipp, j, k, &
                                                ml, zslice
    REAL(real_8)                             :: da1(3,3), ws
    REAL(real_8), ALLOCATABLE                :: wbuff(:), wfull(:,:,:)

    IF (cntl%cdft_wf)THEN           ! Write out the full weight using the cpmd DENSITY functionality
       ALLOCATE(wbuff(fpar%nnr1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(wtemp(fpar%nnr1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(w2(fpar%nnr1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL dcopy(fpar%nnr1,weight,1,wbuff,1)
       !$omp parallel do private(IL)
       DO il=1,fpar%nnr1
          wtemp(il)=CMPLX(wbuff(il),0._real_8,kind=real_8)
       ENDDO
       !$omp end parallel do
       CALL  fwfftn(wtemp,.FALSE.,parai%allgrp)
       CALL zeroing(w2)!,nhg)
       !$omp parallel do private(IL)
       DO il=1,ncpw%nhg
          w2(il)=wtemp(nzh(il))
       ENDDO
       !$omp end parallel do

       IF (paral%io_parent) THEN
          WRITE(6,'(1X,64("-"),/)')
          WRITE(6,'(2X,"WRITING OUT FULL WEIGHT")')
          WRITE(6,'(2X,"to ",A)')"WEIGHT-"//sufF
       ENDIF
       CALL densto(w2,tau0,"WEIGHT-"//suff)
       IF (paral%io_parent) WRITE(6,'(1X,64("-"),/)')
       ! ==--------------------------------------------------------------==
       ! Deallocation of local variables
       DEALLOCATE(wbuff,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(wtemp,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(w2,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ELSE                      ! Write only a given slice of the weight
       ! ==PARA----------------------------------------------------------==
       IF (paral%parent)THEN
          ALLOCATE(wfull(fpar%kr1s,fpar%kr2s,fpar%kr3s),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(wbuff(3*fpar%nnr1),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          !$omp parallel do private(I,J,K) __COLLAPSE3
          DO i=1,fpar%kr1
             DO j=1,fpar%kr2
                DO k=1,fpar%kr3
                   wfull(parap%nrxpl(parai%me,1)+i-1,j,k)=weight(i,j,k)
                ENDDO
             ENDDO
          ENDDO
          !$omp end parallel do
       ENDIF
       IF (parai%nproc.GT.1)THEN
          IF (paral%parent)THEN
             DO ipp=1,parai%nproc
                ip=parap%pgroup(ipp)
                ml=parap%sparm(5,ip)
                ml=ml+MOD(ml+1,2)
                !msglen=(ml*kr2s*kr3s)*8
                IF (ip.NE.parai%me) THEN
                   CALL mp_recv(wbuff,ml*fpar%kr2s*fpar%kr3s,ip,ip,parai%allgrp)
                   DO i=1,parap%sparm(5,ip)
                      DO j=1,fpar%kr2s
                         DO k=1,fpar%kr3s
                            wfull(parap%nrxpl(ip,1)+i,j,k)&
                                 =wbuff(i+(j-1)*ml+(k-1)*ml*fpar%kr2)
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF
             ENDDO
          ELSE
             !msglen=(nnr1)*8
             CALL mp_send(weight,fpar%nnr1,parap%pgroup(1),parai%me,parai%allgrp)
          ENDIF
       ENDIF
       ! ==PARA----------------------------------------------------------==
       IF (paral%io_parent)THEN
          DO i=1,3
             da1(i,1)=parm%a1(i)*0.529177249_real_8/REAL(spar%nr1s,kind=real_8)
             da1(i,2)=parm%a2(i)*0.529177249_real_8/REAL(spar%nr2s,kind=real_8)
             da1(i,3)=parm%a3(i)*0.529177249_real_8/REAL(spar%nr3s,kind=real_8)
          ENDDO
          IF (cdftcom%wslice.LT.0.0_real_8)THEN
             ws=tau0(3,cdftpi%sna(1),cdftpi%spa(1))/parm%a3(3)
          ELSE
             ws=cdftcom%wslice
          ENDIF
          WRITE(6,'(1X,64("-"),/)')
          WRITE(6,'(2X,"WRITING OUT Weight SLICE AT z=",F8.4)')ws*parm%a3(3)
          WRITE(6,'(2X,"to ",A)')"WEIGHT-"//suff//".dat"
          WRITE(6,'(1X,64("-"),/)')
          zslice=ws*fpar%kr3s
          OPEN(72,file="WEIGHT-"//suff//".dat",status='UNKNOWN')
          DO i=1,fpar%kr1s,cdftci%wstep
             DO j=1,fpar%kr2s,cdftci%wstep
                WRITE(72,*)(i)*da1(1,1),(j)*da1(2,2),wfull(i,j,zslice)
             ENDDO
             WRITE(72,*)
          ENDDO
          CLOSE(72)
       ENDIF
       IF (paral%parent)THEN
          ! ==--------------------------------------------------------------==
          ! Deallocation of local variables
          DEALLOCATE(wbuff,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
          DEALLOCATE(wfull,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
          ! ==--------------------------------------------------------------==
       ENDIF
    ENDIF ! Slice or no slice
    RETURN
  END SUBROUTINE write_w
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE wcdft_restart()
    ! ==--------------------------------------------------------------==
    ! == Writes the CDFT_RESTART file                                 ==
    ! ==--------------------------------------------------------------==
    LOGICAL                                  :: logbuff
    REAL(real_8)                             :: oldv2(2)

    IF (paral%io_parent)&
         OPEN(73,file="CDFT_RESTART",status='UNKNOWN')
    IF (paral%io_parent)&
         REWIND(unit=73)
    IF (cdftlog%thda.AND..NOT.cdfthda%hdafirst)THEN
       IF (paral%io_parent)&
            READ(73,*)logbuff
       IF (cdftlog%tcall)THEN
          IF (paral%io_parent)&
               READ(73,*)oldv2
       ELSE
          IF (paral%io_parent)&
               READ(73,*)oldv2(1)
       ENDIF
       IF (paral%io_parent)&
            REWIND(unit=73)
       IF (paral%io_parent)&
            WRITE(73,*).TRUE. ! There is a second record to be written
       IF (cdftlog%tcall)THEN
          IF (paral%io_parent)&
               WRITE(73,*) oldv2(1),oldv2(2)
          IF (paral%io_parent)&
               WRITE(73,*) cdftcom%cdft_v(1),cdftcom%cdft_v(2)
       ELSE
          IF (paral%io_parent)&
               WRITE(73,*) oldv2(1)
          IF (paral%io_parent)&
               WRITE(73,*) cdftcom%cdft_v(1)
       ENDIF
    ELSE
       IF (paral%io_parent)&
            WRITE(73,*).FALSE.! There is no second record to be read
       IF (cdftlog%tcall)THEN
          IF (paral%io_parent)&
               WRITE(73,*) cdftcom%cdft_v(1),cdftcom%cdft_v(2)
       ELSE
          IF (paral%io_parent)&
               WRITE(73,*) cdftcom%cdft_v(1)
       ENDIF
    ENDIF
    IF (paral%io_parent)&
         CLOSE(73)

    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE wcdft_restart
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE rcdft_restart()
    ! ==--------------------------------------------------------------==
    ! == Writes the CDFT_RESTART file                                 ==
    ! ==--------------------------------------------------------------==
    LOGICAL                                  :: logbuff
    REAL(real_8)                             :: ibuff

    IF (paral%io_parent)THEN
       WRITE(6,*)"READING CDFT_V FROM CDFT_RESTART FILE"
       OPEN(73,file="CDFT_RESTART",status='OLD')
       READ(73,*) logbuff
       IF (.NOT.cdftlog%thda)THEN
          IF (.NOT.cdftlog%tcall)THEN
             READ(73,*) ibuff
             cdftcom%cdft_v(1)=ibuff
             cdftcom%cdft_v(2)=0.0_real_8
          ELSE
             READ(73,*) cdftcom%cdft_v
          ENDIF
       ELSE
          IF (cdfthda%hdafirst)THEN
             cdfthda%hdaresb=logbuff
             IF (.NOT.cdftlog%tcall)THEN
                READ(73,*) ibuff
                cdftcom%cdft_v(1)=ibuff
                cdftcom%cdft_v(2)=0.0_real_8
             ELSE
                READ(73,*) cdftcom%cdft_v
             ENDIF
             IF (logbuff)THEN
                IF (.NOT.cdftlog%tcall)THEN
                   READ(73,*) ibuff
                   cdfthda%vbuff(1)=ibuff
                   cdfthda%vbuff(2)=0.0_real_8
                ELSE
                   READ(73,*) cdfthda%vbuff
                ENDIF
             ELSE
                WRITE(6,*) "NO SECOND STATE INFORMATION FOUND! "
             ENDIF
          ELSE
             IF (.NOT.cdftlog%tvmirr) cdftcom%cdft_v=cdfthda%vbuff
          ENDIF
       ENDIF
       CLOSE(73)
    ENDIF
    CALL mp_bcast(cdftcom%cdft_v,SIZE(cdftcom%cdft_v),parai%io_source,parai%cp_grp)
    CALL mp_bcast_byte(cdfthda,size_in_bytes_of(cdfthda),parai%io_source,parai%cp_grp)

    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rcdft_restart
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE cdft_pri()
    ! ==--------------------------------------------------------------==
    ! == Writes some info to term                                     ==
    ! ==--------------------------------------------------------------==      
    INTEGER                                  :: it

    IF (paral%io_parent)&
         WRITE(6,*)
    IF (paral%io_parent)&
         WRITE(6,'(64("-"))')
    IF (.NOT.cdftlog%tspinc)THEN
       IF (paral%io_parent)&
            WRITE(6,'(A,/)') 'CHARGE CONSTRAINED DFT'
    ELSE
       IF (paral%io_parent)&
            WRITE(6,'(A,/)') 'SPIN DENSITY CONSTRAINED DFT'
    ENDIF
    IF (cdftlog%thda)THEN
       WRITE(6,'(2X,A,/)') '2 STATE MATRIX ELEMENT (HAB) CALCULATION'
    ENDIF
    IF (cntl%cdft_newt)THEN
       IF (paral%io_parent)&
            WRITE(6,'(2X,A,/)')'NEWTON OPTIMISER'
    ELSEIF (cntl%cdft_dekk) THEN
       IF (paral%io_parent)&
            WRITE(6,'(2X,A,/)')'DEKKER OPTIMISER'
    ELSE
       IF (paral%io_parent)&
            WRITE(6,'(2X,A,/)')'GRADIENT OPTIMISER'
    ENDIF
    IF (cntl%tpcgfic) THEN
       IF (paral%io_parent)&
            WRITE(6,'(2X,A,/)')'SWITCH TO CNTL%PCG MINIMIZE AT FIRST CDFT STEP'
    ENDIF
    IF (cdftlog%tczones)THEN
       IF (paral%io_parent)&
            WRITE(6,'(2X,A,/)')'WITH CONVERGENCE ZONES'
       IF (paral%io_parent)&
            WRITE(6,'(3X,E10.2,T20,A1,T28,E10.2,T45,A1,T53,E10.2)')&
            czones(1,2),'|',czones(2,2),'|',czones(3,2)
       IF (paral%io_parent)&
            WRITE(6,'(3X,T14,E10.2E1,T39,E10.2E1)')czones(2,1),czones(3,1)
    ENDIF
    IF (cntl%md)THEN
       IF (cdftlog%ccor)THEN
          IF (paral%io_parent)&
               WRITE(6,'(2X,A)')'MD WITH CUTOFF CORRECTION'
       ELSE
          IF (paral%io_parent)&
               WRITE(6,'(2X,A)')'MD WITHOUT CUTOFF CORRECTION'
       ENDIF
       IF (cdftlog%tpred)THEN
          IF (paral%io_parent)&
               WRITE(6,'(2X,A)')'AND V PREDICTION'
          IF (paral%io_parent)&
               WRITE(6,'(2X,A,T54,I10)')'WITH EXTRAPOLATION ORDER k =',&
               cdftpred%predord
       ELSE
          IF (paral%io_parent)&
               WRITE(6,'(/)')
       ENDIF
    ENDIF
    IF (cntl%geopt)THEN
       IF (cdftlog%ccor)THEN
          IF (paral%io_parent)&
               WRITE(6,'(2X,A,/)')'GEOPT WITH CUTOFF CORRECTION'
       ELSE
          IF (paral%io_parent)&
               WRITE(6,'(2X,A,/)')'GEOPT WITHOUT CUTOFF CORRECTION'
       ENDIF
    ENDIF
    IF (paral%io_parent)&
         WRITE(6,'(1X,A,/)') 'PARAMETERS:'
    IF (cdftlog%thda.AND.wgaussl%thdawm)THEN
       IF (paral%io_parent)&
            WRITE(6,'(2X,A)')'STATE A:'
    ENDIF
    IF (cdftpi%ndon.GT.0.AND.paral%io_parent) THEN
       WRITE(6,'(3X,A,T54,I10)')'DONOR ATOM(S):',cdftpi%ndon
    ENDIF
    DO it=1,cdftpi%ndon
       IF (paral%io_parent)&
            WRITE(6,'(4X,T50,I4,T60,A)')cdftpi%cdft_d(it),elem%el(ions0%iatyp(cdftpi%spd(it)))
    ENDDO
    IF (paral%io_parent)&
         WRITE(6,'(3X,A,T54,I10)')'ACCEPTOR ATOM(S):',cdftpi%naccr
    DO it=1,cdftpi%naccr
       IF (paral%io_parent)&
            WRITE(6,'(4X,T50,I4,T60,A)')cdftpi%cdft_a(it),elem%el(ions0%iatyp(cdftpi%spa(it)))
    ENDDO
    IF (cdftlog%thda.AND..NOT.wgaussl%thdawm)THEN
       WRITE(6,'(2X,A)')'STATE A:'
    ENDIF
    IF (cdftpi%ndon.GT.0)THEN
       IF (paral%io_parent)&
            WRITE(6,'(3X,A,T54,F10.3,/)')'CHARGE DIFFERENCE:',cdftcom%cdft_nc
    ELSE
       IF (paral%io_parent)&
            WRITE(6,'(3X,A,T54,F10.3,/)')'CHARGE OF ACCEPTOR:',-cdftcom%cdft_nc
    ENDIF
    IF (cdftlog%thda.AND.wgaussl%thdawm.AND.paral%io_parent)THEN
       WRITE(6,'(2X,A,/)')'STATE B:'
       IF (cdftpi%ndon.GT.0)WRITE(6,'(3X,A,T54,I10)')'DONOR ATOM(S):',cdftpi%ndon
       DO it=1,cdftpi%ndon
          WRITE(6,'(4X,T50,I4,T60,A)')cdftpi%cdft_d(it+cdftpi%ndon),&
               elem%el(ions0%iatyp(cdftpi%spd(it+cdftpi%ndon)))
       ENDDO
       WRITE(6,'(3X,A,T54,I10)')'ACCEPTOR ATOM(S):',cdftpi%naccr
       DO it=1,cdftpi%naccr
          WRITE(6,'(4X,T50,I4,T60,A)')cdftpi%cdft_a(it+cdftpi%naccr),&
               elem%el(ions0%iatyp(cdftpi%spa(it+cdftpi%naccr)))
       ENDDO
    ENDIF
    IF (cdftlog%thda.AND..NOT.wgaussl%thdawm.AND.paral%io_parent)THEN
       WRITE(6,'(2X,A)')'STATE B:'
    ENDIF
    IF (cdftpi%ndon.GT.0)THEN
       IF (paral%io_parent)&
            WRITE(6,'(3X,A,T54,F10.3,/)')'CHARGE DIFFERENCE:',cdftcom%nother
    ELSE
       IF (paral%io_parent)&
            WRITE(6,'(3X,A,T54,F10.3,/)')'CHARGE OF ACCEPTOR:',-cdftcom%nother
    ENDIF

    IF (paral%io_parent)&
         WRITE(6,'(2X,A,T54,I10)')'MAXIMUM NUMBER OF OPTIMISATION STEPS:'&
         ,cdftci%cdft_end
    IF (paral%io_parent)&
         WRITE(6,'(2X,A,T54,F10.3)')'INITIAL LAGRANGE MULTIPLIER(S):',&
         cdftcom%cdft_v(1)
    IF (cdftlog%tcall.AND.paral%io_parent)WRITE(6,'(2X,T54,F10.3)')cdftcom%cdft_v(2)
    IF (paral%io_parent)&
         WRITE(6,'(2X,A,T44,E10.2E1,E10.2E1)')'CONSTRAINT CONVERGENCE:',&
         cdftcom%vccon,cdftcom%vcconu
    IF (wgaussl%twgauss)THEN
       IF (paral%io_parent)&
            WRITE(6,*)
       IF (paral%io_parent)&
            WRITE(6,'(2X,A)')"GAUSSIAN WEIGHT"
       IF (paral%io_parent)&
            WRITE(6,'(2X,A)')"SIGMA:"
       DO it=1,ions1%nsp
          IF (cntl%bohr)THEN
             IF (paral%io_parent)&
                  WRITE(6,'(3X,T35,I4,T45,A,T54,F10.3)')it,elem%el(ions0%iatyp(it)),&
                  wg_sigma(it)
          ELSE
             IF (paral%io_parent)&
                  WRITE(6,'(3X,T35,I4,T45,A,T54,F10.3)')it,elem%el(ions0%iatyp(it)),&
                  wg_sigma(it)/fbohr
          ENDIF
       ENDDO
    ENDIF
    IF (paral%io_parent)&
         WRITE(6,*)
    IF (wcfix%wcut.GT.0.0_real_8)THEN
       IF (cntl%bohr)THEN
          IF (paral%io_parent)&
               WRITE(6,'(2X,A,T54,F10.3)')"WEIGHT CUTOFF",wcfix%wcut
       ELSE
          IF (paral%io_parent)&
               WRITE(6,'(2X,A,T54,F10.3)')"WEIGHT CUTOFF",wcfix%wcut/fbohr
       ENDIF
    ELSE
       IF (paral%io_parent)&
            WRITE(6,'(2X,A)')"VARIABLE CUTOFF"
       DO it=1,ions1%nsp
          IF (cntl%bohr)THEN
             IF (paral%io_parent)&
                  WRITE(6,'(3X,T35,I4,T45,A,T54,F10.3)')it,elem%el(ions0%iatyp(it)),&
                  cdftmd%cdft_rc(it)
          ELSE
             IF (paral%io_parent)&
                  WRITE(6,'(3X,T35,I4,T45,A,T54,F10.3)')it,elem%el(ions0%iatyp(it)),&
                  cdftmd%cdft_rc(it)/fbohr
          ENDIF
       ENDDO

    ENDIF
    IF (paral%io_parent)&
         WRITE(6,'(64("-"))')
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE cdft_pri
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE deriv_init()
    ! ==--------------------------------------------------------------==
    ! == Setup the weight derivatives                                 ==
    ! ==--------------------------------------------------------------==

    CHARACTER(*), PARAMETER                  :: procedureN = 'deriv_init'

    INTEGER                                  :: ierr, ir, is, ish, mmax
    REAL(real_8)                             :: arho(maxsys%mmaxx), &
                                                datom(nsplpo,2), dx, dy, dz, &
                                                fpion, work(maxsys%mmaxx)

    ALLOCATE(wder(maxsys%mmaxx,maxsp,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    dx=parm%a1(1)/REAL(spar%nr1s,kind=real_8)
    dy=parm%a2(2)/REAL(spar%nr2s,kind=real_8)
    dz=parm%a3(3)/REAL(spar%nr3s,kind=real_8)
    cdftmd%cdft_shell=SQRT(dx*dx+dy*dy+dz*dz)

    DO is=1,ions1%nsp
       ! Set everything up
       CALL zeroing(arho)!,maxsys%mmaxx)
       CALL zeroing(datom)!,nsplpo*2)
       CALL zeroing(work)!,maxsys%mmaxx)
       mmax=atwr%meshat(is)
       DO ish=1,atwf_mod%nshell(is)
          !$omp parallel do private(IR) shared(IS,ISH)
          DO ir=1,mmax
             arho(ir)=arho(ir)+&
                  cdftoc%ocun(ish,is)*(atwfr(ir,ish,is)/atrg(ir,is))**2
          ENDDO
          !$omp end parallel do
       ENDDO
       mmax=MIN(mmax,nsplpo)
       cdftmd%cdft_mmax(is)=mmax
       fpion=1._real_8/fpi
       !$omp parallel do private(IR)
       DO ir=1,mmax
          datom(ir,1)=arho(ir)*fpion
       ENDDO
       !$omp end parallel do
       CALL curv1(mmax,atrg(1,is),datom(1,1),0._real_8,0._real_8,3,datom(1,2),&
            work,0._real_8,ierr)
       ! Construct the first derivative for the cutoff correction
       !$omp parallel do private(IR) shared(IS)
       DO ir=1,mmax
          wder(ir,is,1)=atrg(ir,is)
          wder(ir,is,2)=curvd(atrg(ir,is),mmax,atrg(1,is),datom(1,1),&
               datom(1,2),0._real_8)
       ENDDO
       !$omp end parallel do
       CALL zeroing(work)!,maxsys%mmaxx)
       CALL curv1(mmax,wder(1,is,1),wder(1,is,2),0._real_8,0._real_8,3,&
            wder(1,is,3),work,0._real_8,ierr)
    ENDDO

    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE deriv_init
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE cdft_forces(fion,tau0,rho)
    ! ==--------------------------------------------------------------==
    ! == Calculate the Bias Forces                                    ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: fion(:,:,:), tau0(:,:,:), &
                                                rho(:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'cdft_forces'

    INTEGER                                  :: donac, ia, ierr, is, isa, &
                                                isub, it, mmax
    REAL(real_8)                             :: datom(nsplpo,2), fat(3), &
                                                fb(3), mult, rcut, &
                                                wn(fpar%nnr1)
    REAL(real_8), ALLOCATABLE                :: arho(:), work(:)

    CALL tiset('CDFT-FORCE',isub)
    ALLOCATE(arho(maxsys%mmaxx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(work(maxsys%mmaxx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    ! Get Normalisation
    CALL zeroing(wn)!,nnr1)
    DO is=1,ions1%nsp
       CALL atdens_new2(is,datom,mmax,rcut,arho,work)
       DO ia=1,ions0%na(is)
          CALL putrho(wn,tau0(:,ia,is),atrg(1,is),datom,mmax,rcut,1._real_8)
       ENDDO
    ENDDO
    ! Loop over all Atoms
    mult=-parm%omega*cdftcom%cdft_v(1)/(spar%nr1s*spar%nr2s*spar%nr3s)
    DO isa=1,ions1%nat
       ia=iatpt(1,isa)
       is=iatpt(2,isa)
       ! Find out if its in Donor or Acceptor
       donac=0
       DO it=1,cdftpi%naccr
          IF (isa.EQ.cdftpi%cdft_a(it))THEN
             donac=1
             GOTO 101
          ENDIF
       ENDDO
       DO it=1,cdftpi%ndon
          IF (isa.EQ.cdftpi%cdft_d(it))THEN
             donac=2
             GOTO 101
          ENDIF
       ENDDO
101    CONTINUE
       CALL cdft_fint(rho,wn,wdiff,tau0(:,ia,is),cdftmd%cdft_mmax(is),fat,&
            fb,donac,is,cdftmd%cdft_rc(is))
       fion(1,ia,is)=fion(1,ia,is)-mult*fat(1)-cdftcom%cdft_v(1)*fb(1)
       fion(2,ia,is)=fion(2,ia,is)-mult*fat(2)-cdftcom%cdft_v(1)*fb(2)
       fion(3,ia,is)=fion(3,ia,is)-mult*fat(3)-cdftcom%cdft_v(1)*fb(3)
    ENDDO
    DEALLOCATE(arho,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(work,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    CALL tihalt('CDFT-FORCE',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE cdft_forces
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE cdft_fint(rhoe_,psi_,weight_,tau,mmax,fat,fb,donac,sp,&
       rcut)
    ! ==--------------------------------------------------------------==
    ! == Computes the integral for the Bias Forces                    ==
    ! ==--------------------------------------------------------------==
    REAL(real_8), TARGET                     :: rhoe_(:), psi_(:), weight_(:)
    REAL(real_8)                             :: tau(:)
    INTEGER                                  :: mmax
    REAL(real_8)                             :: fat(:), fb(:)
    INTEGER                                  :: donac, sp
    REAL(real_8)                             :: rcut

    INTEGER                                  :: i, i1, ii, iorgx, iorgy, &
                                                iorgz, ix, j, j1, jj, k, k1, &
                                                kk, mr(3), nptau(3)
    REAL(real_8)                             :: da1(3,3), r, rcs, sc1(3), &
                                                sc2(3), wbuf, wbuf2, wn, x, &
                                                x1, xptau(3), y, y1, z, z1
    REAL(real_8), POINTER                    :: psi(:,:,:), rhoe(:,:,:), &
                                                weight(:,:,:)

! Variables
! ==--------------------------------------------------------------==
! ..Atom position

    rhoe(1:fpar%kr1,1:fpar%kr2s,1:SIZE(rhoe_)/(fpar%kr1*fpar%kr2s)) => rhoe_
    psi(1:fpar%kr1,1:fpar%kr2s,1:SIZE(psi_)/(fpar%kr1*fpar%kr2s)) => psi_
    weight(1:fpar%kr1,1:fpar%kr2s,1:SIZE(weight_)/(fpar%kr1*fpar%kr2s)) => weight_

    x1=tau(1)
    y1=tau(2)
    z1=tau(3)
    CALL pbc(x1,y1,z1,sc1(1),sc1(2),sc1(3),1,parm%apbc,parm%ibrav)
    CALL dgemv('T',3,3,1._real_8,metr_com%htm1,3,sc1,1,0._real_8,sc2,1)
    IF (sc2(1).LT.0._real_8) sc2(1)=sc2(1)+1.0_real_8
    IF (sc2(2).LT.0._real_8) sc2(2)=sc2(2)+1.0_real_8
    IF (sc2(3).LT.0._real_8) sc2(3)=sc2(3)+1.0_real_8
    nptau(1)=NINT(sc2(1)*spar%nr1s)+1
    nptau(2)=NINT(sc2(2)*spar%nr2s)+1
    nptau(3)=NINT(sc2(3)*spar%nr3s)+1
    sc1(1)=sc2(1)-REAL(nptau(1)-1,kind=real_8)/REAL(spar%nr1s,kind=real_8)
    sc1(2)=sc2(2)-REAL(nptau(2)-1,kind=real_8)/REAL(spar%nr2s,kind=real_8)
    sc1(3)=sc2(3)-REAL(nptau(3)-1,kind=real_8)/REAL(spar%nr3s,kind=real_8)
    IF (nptau(1).GT.spar%nr1s) nptau(1)=1
    IF (nptau(2).GT.spar%nr2s) nptau(2)=1
    IF (nptau(3).GT.spar%nr3s) nptau(3)=1
    CALL dgemv('T',3,3,1._real_8,metr_com%ht,3,sc1,1,0._real_8,xptau,1)
    DO i=1,3
       da1(i,1)=parm%a1(i)/REAL(spar%nr1s,kind=real_8)
       da1(i,2)=parm%a2(i)/REAL(spar%nr2s,kind=real_8)
       da1(i,3)=parm%a3(i)/REAL(spar%nr3s,kind=real_8)
    ENDDO
    ! ..Interaction region
    sc1(1)=rcut
    sc1(2)=rcut
    sc1(3)=rcut
    CALL dgemv('T',3,3,1._real_8,metr_com%htm1,3,sc1,1,0._real_8,sc2,1)
    mr(1)=NINT(sc2(1)*spar%nr1s)+1
    mr(2)=NINT(sc2(2)*spar%nr2s)+1
    mr(3)=NINT(sc2(3)*spar%nr3s)+1
    mr(1)=MIN(mr(1),spar%nr1s/2-1)
    mr(2)=MIN(mr(2),spar%nr2s/2-1)
    mr(3)=MIN(mr(3),spar%nr3s/2-1)
    ! ..Origin of density grid
    iorgx=nptau(1)-mr(1)
    iorgy=nptau(2)-mr(2)
    iorgz=nptau(3)-mr(3)
    IF (iorgx.LT.1) iorgx=iorgx+spar%nr1s
    IF (iorgy.LT.1) iorgy=iorgy+spar%nr2s
    IF (iorgz.LT.1) iorgz=iorgz+spar%nr3s
    ! ..loop over all grid points
    DO i=1,3
       fat(i)=0.0_real_8
       fb(i)=0.0_real_8
    ENDDO
    rcs=rcut-cdftmd%cdft_shell
    DO i=0,2*mr(1)
       ix=iorgx+i
       IF (ix.GT.spar%nr1s) ix=ix-spar%nr1s
       IF (ix.GE.parap%nrxpl(parai%mepos,1).AND.ix.LE.parap%nrxpl(parai%mepos,2)) THEN
          ii=ix-parap%nrxpl(parai%mepos,1)+1
          i1=mr(1)-i
          DO j=0,2*mr(2)
             jj=iorgy+j
             IF (jj.GT.spar%nr2s) jj=jj-spar%nr2s
             j1=mr(2)-j
             DO k=0,2*mr(3)
                kk=iorgz+k
                IF (kk.GT.spar%nr3s) kk=kk-spar%nr3s
                k1=mr(3)-k
                CALL positx(x,y,z,i1,j1,k1,da1,xptau(1))
                IF (isos1%tclust)THEN
                   IF (x1-x.GE.parm%a1(1).OR.y1-y.GE.parm%a2(2).OR.z1-z.GE.parm%a3(3)&
                        .OR.x1-x.LE.0.OR.y1-y.LE.0.OR.z1-z.LE.0)THEN
                      GOTO 201
                   ENDIF
                ENDIF
                r=SQRT(x*x+y*y+z*z)
                wn=psi(ii,jj,kk)
                IF (r.LT.rcut.AND.r.GT.1.0e-10_real_8)THEN
                   IF (donac.EQ.0)THEN
                      wbuf=weight(ii,jj,kk)
                   ELSEIF (donac.EQ.1)THEN
                      wbuf=weight(ii,jj,kk)-1! sign change due to optimisation ...
                   ELSEIF (donac.EQ.2)THEN
                      wbuf=1+weight(ii,jj,kk)! sign change due to optimisation ...
                   ENDIF
                   wbuf2=curv2(r,mmax,wder(1,sp,1),wder(1,sp,2),&
                        wder(1,sp,3),0._real_8)*wbuf*rhoe(ii,jj,kk)/wn
                   fat(1)=fat(1)+wbuf2*x/r
                   fat(2)=fat(2)+wbuf2*y/r
                   fat(3)=fat(3)+wbuf2*z/r
                   IF (cdftlog%ccor.AND.r.GE.rcs) THEN
                      wbuf2=rhoe(ii,jj,kk)*wbuf*surfe(x,y,z,rcut,da1)/wn
                      fb(1)=fb(1)+wbuf2*x
                      fb(2)=fb(2)+wbuf2*y
                      fb(3)=fb(3)+wbuf2*z
                   ENDIF
                ENDIF
201             CONTINUE
             ENDDO
          ENDDO
       ENDIF
    ENDDO
    ! !    DO I=1,3
    ! !      FB(I)=CDFTMD%RHOCUT(SP)*RCUT*FB(I)
    ! !    ENDDO
    fb(1)=cdftmd%rhocut(sp)*rcut*fb(1)
    fb(2)=cdftmd%rhocut(sp)*rcut*fb(2)
    fb(3)=cdftmd%rhocut(sp)*rcut*fb(3)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE cdft_fint
  ! ==================================================================

  ! ==================================================================
  FUNCTION surfe(x,y,z,r,dx)
    ! ==--------------------------------------------------------------==
    ! == Calculate the Surface Element                                ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: x, y, z, r, dx(3,3), surfe

    surfe=SIGN(1.0_real_8,z)*&
         (y*dx(1,1)*dx(3,3)-x*dx(2,2)*dx(3,3))/(r**3-r*z**2)
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION surfe
  ! ==================================================================
  SUBROUTINE atdens_new2(is,datom,mmax,rcut,arho,work)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: is
    REAL(real_8)                             :: datom(nsplpo,2)
    INTEGER                                  :: mmax
    REAL(real_8)                             :: rcut, arho(maxsys%mmaxx), &
                                                work(maxsys%mmaxx)

    INTEGER                                  :: ierr, ir, ish, ocspec
    REAL(real_8)                             :: fpion, isig, norm

    CALL zeroing(arho)!,maxsys%mmaxx)
    mmax=atwr%meshat(is)
    ! ..Sum density
    isig=wg_sigma(is)*SQRT(2.0_real_8)

    IF (wgaussl%twgauss)THEN
       norm=1.0_real_8/(isig*SQRT(pi))
       ocspec=0
       !$omp  parallel do private(ISH) &
       !$omp  reduction(+:OCSPEC)
       DO ish=1,atwf_mod%nshell(is)
          ocspec=ocspec+cdftoc%ocun(ish,is)
       ENDDO
       !$omp  end parallel do
       !$omp  parallel do private(IR) 
       DO ir=1,mmax
          arho(ir)=ocspec*norm*EXP(-(atrg(ir,is)/isig)**2)
       ENDDO
       !$omp  end parallel do
    ELSE
       DO ish=1,atwf_mod%nshell(is)
          !$omp  parallel do private(IR) shared(IS,ISH)
          DO ir=1,mmax
             arho(ir)=arho(ir)+&
                  cdftoc%ocun(ish,is)*(atwfr(ir,ish,is)/atrg(ir,is))**2
          ENDDO
          !$omp  end parallel do
       ENDDO
    ENDIF
    mmax=MIN(mmax,nsplpo)
    fpion=1._real_8/fpi
    !$omp  parallel do private(IR)
    DO ir=1,mmax
       datom(ir,1)=arho(ir)*fpion
    ENDDO
    !$omp  end parallel do
    CALL curv1(mmax,atrg(1,is),datom(1,1),0._real_8,0._real_8,3,datom(1,2),work,&
         0._real_8,ierr)
    IF (wcfix%wcut.GT.0)THEN
       rcut=wcfix%wcut
       cdftmd%rhocut(is)=curv2(rcut,mmax,atrg(1,is),datom(1,1),datom(1,2),&
            0._real_8)
    ELSE
       DO ir=1,mmax
          IF (datom(ir,1).LT.1.e-6_real_8.AND.atrg(ir,is).GT.1._real_8) THEN
             rcut=atrg(ir,is)
             cdftmd%rhocut(is)=datom(ir,1)
             GOTO 100
          ENDIF
       ENDDO
       rcut=atrg(mmax,is)
       IF (paral%io_parent)&
            WRITE(6,*) &
            "REFERENCE DENSITY DECAYS TOO SLOW, Rc SET TO MAXIMUM! "
100    CONTINUE
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE atdens_new2
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE ringbuff(vin,full)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: vin
    LOGICAL                                  :: full

    INTEGER                                  :: elc, i

! INPUT V VALUE
! returns true if the number of elements in the buffer=CDFTPRED%PREDORD
! Variables
! number of elements in the buffer

    COMMON/rbuff/elc

    full=.FALSE.

    IF (elc.LT.cdftpred%predord+1)THEN
       cdftpred%vpredbuf(elc+1)=vin
       elc=elc+1
    ELSE
       DO i=1,cdftpred%predord
          cdftpred%vpredbuf(i)=cdftpred%vpredbuf(i+1)
       ENDDO
       cdftpred%vpredbuf(cdftpred%predord+1)=vin
    ENDIF

    IF (elc.EQ.cdftpred%predord+1)full=.TRUE.
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE ringbuff
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE v_pred()
    ! ==--------------------------------------------------------------==
    ! Variables
    INTEGER                                  :: i, j
    LOGICAL                                  :: full
    REAL(real_8)                             :: l, v

! V VALUE

    v=0.0_real_8

    IF (paral%io_parent)THEN
       CALL ringbuff(cdftcom%cdft_v(1),full)
       IF (full) THEN
          DO j=0,cdftpred%predord
             l=1.0_real_8
             DO i=0,cdftpred%predord
                IF (i.NE.j)THEN
                   l=l*REAL(cdftpred%predord+1-i,kind=real_8)/REAL(j-i,kind=real_8)
                ENDIF
             ENDDO
             v=v+l*cdftpred%vpredbuf(j+1)
          ENDDO
          cdftcom%cdft_v(1)=v
          WRITE(6,*)"PREDICTED V TO BE: ",V
       ELSE
          WRITE(6,*)"TOO FEW VALUES FOR PREDICTION,"
          WRITE(6,*)"USING OLD V! "
       ENDIF
    ENDIF
    CALL mp_bcast(cdftcom%cdft_v,SIZE(cdftcom%cdft_v),parai%io_source,parai%cp_grp)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE v_pred
  ! ==================================================================

END MODULE cdft_utils
