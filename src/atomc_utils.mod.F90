MODULE atomc_utils
  USE atwf,                            ONLY: atrg,&
                                             atwf_mod,&
                                             atwfr,&
                                             atwp,&
                                             atwr
  USE cdftmod,                         ONLY: cdftmd
  USE cnst,                            ONLY: fpi
  USE error_handling,                  ONLY: stopgm
  USE fitpack_utils,                   ONLY: curv1,&
                                             curv2
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE isos,                            ONLY: isos1,&
                                             isos3
  USE kinds,                           ONLY: real_8
  USE metr,                            ONLY: metr_com
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai
  USE pbc_utils,                       ONLY: pbc
  USE qspl,                            ONLY: nsplpo
  USE system,                          ONLY: fpar,&
                                             maxsys,&
                                             parap,&
                                             parm,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: atomc
  !public :: atdens_new
  PUBLIC :: putrho
  PUBLIC :: positx
  !public :: getchg

CONTAINS

  ! ==================================================================
  SUBROUTINE atomc(rhoe,psi,tau0,achrg)
    ! ==--------------------------------------------------------------==
    ! == COMPUTES ATOMIC CHARGES BY INTEGRATION WITH WEIGHT FUNCTIONS ==
    ! ==--------------------------------------------------------------==

    REAL(real_8)                             :: rhoe(:), psi(:), tau0(:,:,:), &
                                                achrg(:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'atomc'

    INTEGER                                  :: ia, iat, ierr, is, ish, isub, &
                                                l, mmax
    INTEGER, SAVE                            :: ifirst = 0
    REAL(real_8)                             :: chrg, rcut
    REAL(real_8), ALLOCATABLE                :: arho(:), work(:)
    REAL(real_8), ALLOCATABLE, SAVE          :: datom(:,:)

!(3,maxsys%nax,maxsys%nsx), &
!ions1%nat)
! ==--------------------------------------------------------------==

    CALL tiset('     ATOMC',isub)
    ! ..initialization
    IF (ifirst.EQ.0) THEN
       atwp%nattot=0
       DO is=1,ions1%nsp
          atwf_mod%numaor(is)=0
          DO ish=1,atwf_mod%nshell(is)
             l=atwf_mod%lshell(ish,is)
             atwp%nattot=atwp%nattot+ions0%na(is)*(2*l+1)
             atwf_mod%numaor(is)=atwf_mod%numaor(is)+(2*l+1)
          ENDDO
       ENDDO
       ALLOCATE(datom(nsplpo,2),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ifirst=1
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Allocation of local variables
    ALLOCATE(arho(maxsys%mmaxx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(work(maxsys%mmaxx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL zeroing(psi)!,nnr1)
    ! ..sum all partition functions in psi
    DO is=1,ions1%nsp
       CALL atdens_new(is,datom,mmax,rcut,arho,work)
       DO ia=1,ions0%na(is)
          CALL putrho(psi,tau0(:,ia,is),atrg(:,is),datom,mmax,rcut,1._real_8)
       ENDDO
    ENDDO
    ! ..charge integration
    iat=0
    DO is=1,ions1%nsp
       CALL atdens_new(is,datom,mmax,rcut,arho,work)
       DO ia=1,ions0%na(is)
          iat=iat+1
          CALL getchg(rhoe,psi,tau0(:,ia,is),atrg(:,is),datom,&
               mmax,rcut,chrg)
          achrg(iat)=chrg
       ENDDO
    ENDDO
    CALL mp_sum(achrg,ions1%nat,parai%allgrp)
    iat=0
    chrg=parm%omega/REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          iat=iat+1
          achrg(iat)=ions0%zv(is)-achrg(iat)*chrg
       ENDDO
    ENDDO
    CALL tihalt('     ATOMC',isub)
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
  END SUBROUTINE atomc
  ! ==================================================================
  SUBROUTINE atdens_new(is,datom,mmax,rcut,arho,work)
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: is
    REAL(real_8)                             :: datom(nsplpo,2)
    INTEGER                                  :: mmax
    REAL(real_8)                             :: rcut, arho(maxsys%mmaxx), &
                                                work(maxsys%mmaxx)

    INTEGER                                  :: ierr, ir, ish

    CALL zeroing(arho)!,maxsys%mmaxx)
    mmax=atwr%meshat(is)
    ! ..Sum density
    DO ish=1,atwf_mod%nshell(is)
       DO ir=1,mmax
          arho(ir)=arho(ir)+&
               atwf_mod%oc(ish,is) * (atwfr(ir,ish,is)/atrg(ir,is))**2
       ENDDO
    ENDDO
    mmax=MIN(mmax,nsplpo)
    !$omp parallel do private(IR)
    DO ir=1,mmax
       datom(ir,1)=arho(ir)/fpi
    ENDDO
    CALL curv1(mmax,atrg(1,is),datom(1,1),0._real_8,0._real_8,3,datom(1,2),&
         work,0._real_8,ierr)
    DO ir=1,mmax
       IF (datom(ir,1).LT.1.e-6_real_8.AND.atrg(ir,is).GT.1._real_8) THEN
          ! IF(DATOM(IR,1).LT.1.e-14_real_8.AND.ATRG(IR,IS).GT.1._real_8) THEN
          rcut=atrg(ir,is)
          cdftmd%rhocut(is)=datom(ir,1)
          GOTO 100
       ENDIF
    ENDDO
    rcut=atrg(mmax,is)
100 CONTINUE
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE atdens_new
  ! ==================================================================
  SUBROUTINE putrho(rhoe,tau,rg,datom,mmax,rcut,cxscal)
    ! ==--------------------------------------------------------------==

    REAL(real_8)                             :: rhoe(fpar%kr1,fpar%kr2,*), &
                                                tau(*), rg(*), datom(nsplpo,2)
    INTEGER                                  :: mmax
    REAL(real_8)                             :: rcut, cxscal

    INTEGER                                  :: i, i1, igrid, ii, iorgx, &
                                                iorgy, iorgz, ix, j, j1, jj, &
                                                k, k1, kk, mr(3), nptau(3)
    LOGICAL                                  :: tadd
    REAL(real_8)                             :: da1(3,3), r, sc1(3), sc2(3), &
                                                x, x1, xptau(3), y, y1, z, z1

! Variables
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
                         ! ..Add density to grid
                         igrid=igrid+1
                         rhoe(ii,jj,kk)=rhoe(ii,jj,kk)+&
                              cxscal*curv2(r,mmax,rg,datom(1,1),datom(1,2),0._real_8)
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
                         ! ..Add density to grid
                         igrid=igrid+1
                         rhoe(ii,jj,kk)=rhoe(ii,jj,kk)+&
                              cxscal*curv2(r,mmax,rg,datom(1,1),datom(1,2),0._real_8)
                      ENDIF
                   ENDIF
                ELSE
                   r=SQRT(x*x+y*y+z*z)
                   IF (r.LT.rcut) THEN
                      ! ..Add density to grid
                      igrid=igrid+1
                      rhoe(ii,jj,kk)=rhoe(ii,jj,kk)+&
                           cxscal*curv2(r,mmax,rg,datom(1,1),datom(1,2),0._real_8)
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDIF
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE putrho
  ! ==================================================================
  SUBROUTINE positx(x,y,z,i,j,k,da,xp)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: x, y, z
    INTEGER                                  :: i, j, k
    REAL(real_8)                             :: da(3,3), xp(3)

! ==--------------------------------------------------------------==

    x=i*da(1,1)+j*da(1,2)+k*da(1,3)+xp(1)
    y=i*da(2,1)+j*da(2,2)+k*da(2,3)+xp(2)
    z=i*da(3,1)+j*da(3,2)+k*da(3,3)+xp(3)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE positx
  ! ==================================================================
  SUBROUTINE getchg(rhoe,psi,tau,rg,datom,mmax,rcut,chg)
    ! ==--------------------------------------------------------------==

    REAL(real_8)                             :: rhoe(fpar%kr1,fpar%kr2,*), &
                                                psi(fpar%kr1,fpar%kr2,*), &
                                                tau(3), rg(*), datom(nsplpo,2)
    INTEGER                                  :: mmax
    REAL(real_8)                             :: rcut, chg

    INTEGER                                  :: i, i1, igrid, ii, iorgx, &
                                                iorgy, iorgz, ix, j, j1, jj, &
                                                k, k1, kk, mr(3), nptau(3)
    REAL(real_8)                             :: da1(3,3), r, sc1(3), sc2(3), &
                                                x, x1, xptau(3), y, y1, z, z1

! Variables
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
                r=SQRT(x*x+y*y+z*z)
                IF (r.LT.rcut.AND.psi(ii,jj,kk).NE.0.0_real_8) THEN
                   igrid=igrid+1
                   chg=chg+rhoe(ii,jj,kk)/psi(ii,jj,kk)*&
                        curv2(r,mmax,rg,datom(1,1),datom(1,2),0._real_8)
                ENDIF
             ENDDO
          ENDDO
       ENDIF
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE getchg
  ! ==================================================================

END MODULE atomc_utils
