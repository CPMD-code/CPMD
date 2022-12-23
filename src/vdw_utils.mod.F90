#include "cpmd_global.h"

MODULE vdw_utils
  USE adat,                            ONLY: covrad
  USE cnst,                            ONLY: fbohr
  USE dftd3_driver,                    ONLY: dftd3
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_def,&
                                             fo_old
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE meta_multiple_walkers_utils,     ONLY: mw_filename
  USE metr,                            ONLY: metr_com
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum
  USE mw,                              ONLY: mwi,&
                                             tmw
  USE parac,                           ONLY: parai,&
                                             paral
  USE pbc_utils,                       ONLY: pbc, dist_pbc
  USE pimd,                            ONLY: ipcurr,&
                                             np_local,&
                                             np_low
  USE sort_utils,                      ONLY: sort2
  USE strs,                            ONLY: alpha,&
                                             beta
  USE system,                          ONLY: cntl,&
                                             iatpt,&
                                             maxsys,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE vdwcmod,                         ONLY: &
       boadwf, icontfragw, icontfragwi, ifragdata, ifragw, iwfcref, natwfcx, &
       nfrags, nfragx, npt12, nwfcx, radfrag, rwann, rwfc, spr, swann, &
       taufrag, tauref, twannupx, vdwi, vdwr, vdwwfi, vdwwfl, vdwwfr, wwfcref, &
       tdftd3
  USE wrgeo_utils,                     ONLY: wrgeof
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: vdw
  PUBLIC :: vdw_wf
  PUBLIC :: wwannier
  PUBLIC :: rwannier

CONTAINS

  SUBROUTINE vdw(tau0,nvdw,idvdw,ivdw,jvdw,vdwst,vdwrm,vdwbe,&
       VDWEPS,S6GRIM,NXVDW,NYVDW,NZVDW,EVDW,FION,DEVDW)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! == THE EMPIRICAL VAN DER WAALS CORRECTION TO THE TOTAL ENERGY,  == 
    ! == IONIC FORCE, AND STRESS TENSOR ACCORDING TO THE FORMUALTION  ==
    ! == OF R. LE SAR, J. PHYS. CHEM. 88, p. 4272 (1984)              ==
    ! ==               J. CHEM. PHYS. 86, p. 1485 (1987)              ==
    ! == Modified as in:                                              == 
    ! ==   M. Elstner et al. J. Chem. Phys. 114, p. 5149 (2001)       ==
    ! ==   or Grimme, J. Comput. Chem. 27, 1787, 2006 (IDVDW=3)       ==
    ! == Present release: Munich/Philadelphia/Tsukuba, 11 May 2010    ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:)
    INTEGER                                  :: nvdw
    INTEGER, ALLOCATABLE                     :: idvdw(:), ivdw(:), jvdw(:)
    REAL(real_8), ALLOCATABLE                :: vdwst(:), vdwrm(:), vdwbe(:)
    REAL(real_8)                             :: VDWEPS, S6GRIM
    INTEGER                                  :: nxvdw, nyvdw, nzvdw
    REAL(real_8)                             :: EVDW, fion(:,:,:), devdw(:)

    INTEGER                                  :: iv, jv, l, l_upper, m, &
                                                m_upper, n, nx, ny, nz
    REAL(real_8)                             :: aexp, ffac, &
                                                reps = 1.0e-2_real_8, rexp, &
                                                rfac, rlm1, rlm2, rlm6, &
                                                rlp(3), rphi, rpow, xlm, ylm, &
                                                zlm, xlm_, ylm_, zlm_

! ==--------------------------------------------------------------==

    evdw=0.0_real_8
    CALL zeroing(devdw)!,6)
    IF (paral%parent) THEN
       IF (tdftd3) THEN
          ! dftd3
          CALL dftd3(tau0,nxvdw,nyvdw,nzvdw,evdw,fion,devdw)
       ELSE
          ! dftd2
          DO n=1,nvdw
             iv=ivdw(n)
             jv=jvdw(n)
             l_upper=ions0%na(iv)
             m_upper=ions0%na(jv)
             DO l=1,l_upper
                IF (iv.EQ.jv) THEN
                   m_upper=l
                ENDIF
                DO m=1,m_upper
                   xlm_=tau0(1,l,iv)-tau0(1,m,jv)
                   ylm_=tau0(2,l,iv)-tau0(2,m,jv)
                   zlm_=tau0(3,l,iv)-tau0(3,m,jv)
                   CALL pbc(xlm_,ylm_,zlm_,xlm,ylm,zlm,1,parm%apbc,parm%ibrav)
                   DO nx=-nxvdw,nxvdw
                      DO ny=-nyvdw,nyvdw
                         DO nz=-nzvdw,nzvdw
                            rlp(1)=xlm+nx*metr_com%ht(1,1)+ny*metr_com%ht(2,1)+nz*metr_com%ht(3,1)
                            rlp(2)=ylm+nx*metr_com%ht(1,2)+ny*metr_com%ht(2,2)+nz*metr_com%ht(3,2)
                            rlp(3)=zlm+nx*metr_com%ht(1,3)+ny*metr_com%ht(2,3)+nz*metr_com%ht(3,3)
                            rlm2=rlp(1)*rlp(1)+rlp(2)*rlp(2)+rlp(3)*rlp(3)
                            rlm1=SQRT(rlm2)
                            IF (rlm1.LT.reps) THEN
                               rphi=0.0_real_8
                               ffac=0.0_real_8
                            ELSE IF (idvdw(n).EQ.1) THEN
                               rlm6=rlm2*rlm2*rlm2
                               rpow=vdwbe(n)*(rlm1/vdwrm(n))**7
                               IF (ABS(rpow).LT.1.0_real_8/vdweps) THEN
                                  rexp=EXP(-rpow)
                               ELSE
                                  rexp=0.0_real_8
                               ENDIF
                               aexp=1.0_real_8-rexp
                               rphi=vdwst(n)/rlm6*aexp**4
                               rfac=-6.0_real_8/rlm1+28.0_real_8*rpow*rexp/(rlm1*aexp)
                               ffac=rphi*rfac/rlm1
                            ELSE IF (idvdw(n).EQ.2) THEN
                               rlm6=rlm2*rlm2*rlm2
                               rpow=vdwbe(n)*(rlm1-vdwrm(n))
                               aexp=vdwst(n)*0.5_real_8
                               IF (ABS(rpow).LT.1.0_real_8/vdweps) THEN
                                  rphi=aexp*(1.0_real_8+TANH(rpow))/rlm6
                                  rfac=-6.0_real_8/rlm1+vdwbe(n)/(COSH(rpow)**2*(1.0_real_8+TANH(rpow)))
                                  ffac=rphi*rfac/rlm1
                               ELSE
                                  rphi=0.0_real_8
                                  ffac=0.0_real_8
                               ENDIF
                               ! -----------------------------------------------------------------------
                               ! Dispersion energy and forces according to Grimme
                            ELSE IF (idvdw(n).EQ.3) THEN
                               rlm6=rlm2*rlm2*rlm2
                               rpow=vdwbe(n)*(rlm1/vdwrm(n)-1.0_real_8)
                               IF (ABS(rpow).LT.1.0_real_8/vdweps) THEN
                                  rexp=1.0_real_8+EXP(-rpow)
                                  rphi=-s6grim*vdwst(n)/(rlm6*rexp)
                                  ffac=s6grim*vdwst(n)/(rlm6*rexp*rlm1)&
                                       *(6.0_real_8/rlm1-vdwbe(n)*(rexp-1.0_real_8)&
                                       /(vdwrm(n)*rexp))
                               ELSE
                                  rphi=0.0_real_8
                                  ffac=0.0_real_8
                               ENDIF
                               ! ------------------------------------------------------------------------
                            ENDIF
                            ! Energy
                            evdw=evdw+rphi
                            ! Ionic forces 
                            fion(1,l,iv)=fion(1,l,iv)-ffac*rlp(1)
                            fion(2,l,iv)=fion(2,l,iv)-ffac*rlp(2)
                            fion(3,l,iv)=fion(3,l,iv)-ffac*rlp(3)
                            fion(1,m,jv)=fion(1,m,jv)+ffac*rlp(1)
                            fion(2,m,jv)=fion(2,m,jv)+ffac*rlp(2)
                            fion(3,m,jv)=fion(3,m,jv)+ffac*rlp(3)
                            ! Stress tensor 
                            IF (.NOT.(ANY(alpha==0).OR.ANY(beta==0))) THEN
                               devdw(1)=devdw(1)&
                                    +FFAC*RLP(ALPHA(1))*RLP(BETA(1))
                               devdw(2)=devdw(2)&
                                    +FFAC*RLP(ALPHA(2))*RLP(BETA(2))
                               devdw(3)=devdw(3)&
                                    +FFAC*RLP(ALPHA(3))*RLP(BETA(3))
                               devdw(4)=devdw(4)&
                                    +FFAC*RLP(ALPHA(4))*RLP(BETA(4))
                               devdw(5)=devdw(5)&
                                    +FFAC*RLP(ALPHA(5))*RLP(BETA(5))
                               devdw(6)=devdw(6)&
                                    +FFAC*RLP(ALPHA(6))*RLP(BETA(6))
                            ENDIF
                         ENDDO
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDIF
    ENDIF
    ! ==================================================================
    RETURN
  END SUBROUTINE vdw
  ! ==================================================================
  ! ==     Next part = Wannier Functions/Centers Stuff for vdW      ==
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE vdw_wf(tau0,fion,nstate)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! == VAN DER WAALS CORRECTION BASED ON THE USE OF MAXIMALLY       == 
    ! == LOCALIZED WANNIER FUNCTIONS AS PROPOSED IN                   ==
    ! == P. L. SILVESTRELLI, PHYS. REV. LETT. 100, 053002 (2008)      ==
    ! ==                     J. PHYS. CHEM. A 113, p.5224 (2009)      ==
    ! == A SECOND APPROACH IS ALSO IMPLEMENTED ACCORDING TO           ==
    ! ==                     PHYS. REV. B 85, 073101 (2012)           == 
    ! == Present Release: Strasbourg/Padova/Tokyo   15 May 2012       ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:), fion(:,:,:)
    INTEGER                                  :: nstate

    CHARACTER(*), PARAMETER                  :: procedureN = 'vdw_wf'

    INTEGER                                  :: i, i1, ia, ib, icfr, ierr, &
                                                in, ip, ipx, is, isub, jcfr, &
                                                l, m, numx, nx, ny, nz
    INTEGER, SAVE                            :: ifirst = 0
    REAL(real_8) :: c6, c6fac, c6vdw, datauf, distn, distn6, dr1, dr2, dr3, &
      dxmax, dxnorm, efac, ffac, reps = 1.0_real_8, rlp(3), seffi, seffi3, &
      seffi32, seffj, seffj3, seffj32, si3, si3d2, sj3, temp1, temp2, temp3, &
      veff1, vol1, voleffwan, vtot1, xlm, ylm, zlm, xlm_, ylm_, zlm_
    REAL(real_8), ALLOCATABLE                :: fionwf(:,:,:), fove(:)
    REAL(real_8), ALLOCATABLE, SAVE          :: c6all(:,:,:)

! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    vdwr%evdw=0.0_real_8
    CALL zeroing(vdwr%devdw)!,6)
    !   IF (.NOT.paral%io_parent) GOTO 9999
    IF (cntl%tpath) THEN
       numx=np_local; ipx=ipcurr-np_low+1
    ELSE
       numx=1; ipx=1
    ENDIF
    IF (cntl%tpath) vdwwfl%twannup=twannupx(ipx)
    IF (vdwwfl%twannup) THEN
       IF (ifirst.EQ.0) THEN
          ALLOCATE(c6all(nstate,nstate,numx),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ifirst=1
       ENDIF
       ! Construct fragments acoording to the selected criterion 
       IF (vdwwfi%icriteri.EQ.0) THEN    ! the z-level criterion
          CALL fragment_zlevel(nstate,rwann(:,:,ipx),swann(:,ipx),&
               tauref(:,:,:,ipx),icontfragw(:,ipx),icontfragwi(:,:,ipx),&
               iwfcref(:,:,:,ipx),ifragw(:,:,ipx),wwfcref(:,:,:,ipx),&
               rwfc(:,:,:,ipx),spr(:,:,ipx),taufrag(:,:,:,ipx))
       ELSE IF (vdwwfi%icriteri.EQ.1) THEN! the radius criterion
          CALL fragment_radius(nstate,rwann(:,:,ipx),swann(:,ipx),&
               tauref(:,:,:,ipx),icontfragw(:,ipx),icontfragwi(:,:,ipx),&
               iwfcref(:,:,:,ipx),ifragw(:,:,ipx),wwfcref(:,:,:,ipx),&
               rwfc(:,:,:,ipx),spr(:,:,ipx),taufrag(:,:,:,ipx))
       ELSE IF (vdwwfi%icriteri.EQ.2) THEN! the bond criterion
          CALL fragment_bond(nstate,rwann(:,:,ipx),swann(:,ipx),&
               tauref(:,:,:,ipx),icontfragw(:,ipx),icontfragwi(:,:,ipx),&
               iwfcref(:,:,:,ipx),ifragw(:,:,ipx),wwfcref(:,:,:,ipx),&
               rwfc(:,:,:,ipx),spr(:,:,ipx),taufrag(:,:,:,ipx))
       ELSE
          CALL stopgm(procedureN,'WRONG ICRITERI! ',& 
               __LINE__,__FILE__)
       ENDIF
       nfrags(ipx)=vdwwfi%multifrag
    ENDIF
    ! Compute reference coordinates for each WF and move the positions of
    ! the WFCs according to the displacement of their own reference points 
    !$omp parallel do private(ICFR,JCFR,I,DR1,DR2,DR3,TEMP1,TEMP2,TEMP3, &
    !$omp  IP,IN,IA,IS)
    DO icfr=1,nfrags(ipx)
       DO jcfr=1,icontfragw(icfr,ipx)
          dr1=taufrag(1,jcfr,icfr,ipx)
          dr2=taufrag(2,jcfr,icfr,ipx)
          dr3=taufrag(3,jcfr,icfr,ipx)
          temp1=0.0_real_8
          temp2=0.0_real_8
          temp3=0.0_real_8
          ip=icontfragwi(jcfr,icfr,ipx)
          DO i=1,ip
             in=iwfcref(i,jcfr,icfr,ipx)
             ia=iatpt(1,in)
             is=iatpt(2,in)
             ! temp1=temp1+tau0(1,ia,is)
             ! temp2=temp2+tau0(2,ia,is)
             ! temp3=temp3+tau0(3,ia,is)
             temp1=temp1+tau0(1,ia,is)*wwfcref(i,jcfr,icfr,ipx)
             temp2=temp2+tau0(2,ia,is)*wwfcref(i,jcfr,icfr,ipx)
             temp3=temp3+tau0(3,ia,is)*wwfcref(i,jcfr,icfr,ipx)
          ENDDO
          ! taufrag(1,jcfr,icfr)=temp1/REAL(ip,kind=real_8)
          ! taufrag(2,jcfr,icfr)=temp2/REAL(ip,kind=real_8)
          ! taufrag(3,jcfr,icfr)=temp3/REAL(ip,kind=real_8)
          taufrag(1,jcfr,icfr,ipx)=temp1
          taufrag(2,jcfr,icfr,ipx)=temp2
          taufrag(3,jcfr,icfr,ipx)=temp3
          dr1=taufrag(1,jcfr,icfr,ipx)-dr1
          dr2=taufrag(2,jcfr,icfr,ipx)-dr2
          dr3=taufrag(3,jcfr,icfr,ipx)-dr3
          rwfc(1,jcfr,icfr,ipx)=rwfc(1,jcfr,icfr,ipx)+dr1
          rwfc(2,jcfr,icfr,ipx)=rwfc(2,jcfr,icfr,ipx)+dr2
          rwfc(3,jcfr,icfr,ipx)=rwfc(3,jcfr,icfr,ipx)+dr3
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    ! Compute C6 coefficients
    IF (vdwwfi%iswitchvdw.EQ.2) THEN
       ! vdw-WF2: enforce exact H polarizability
       ! Compute intrafragment overlap
       IF (vdwwfl%twannup) THEN
          ALLOCATE(fove(nfragx),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          DO icfr=1,nfrags(ipx)      ! scan all fragments
             i1=icontfragw(icfr,ipx) ! fragment centers
             veff1=0.0_real_8
             vtot1=0.0_real_8
             DO ib=1,i1
                CALL volwan(rwfc(1,1,icfr,ipx),spr(1,icfr,ipx),i1,ib,&
                     voleffwan,vol1)
                veff1=veff1+voleffwan
                vtot1=vtot1+vol1
             ENDDO
             fove(icfr)=(veff1/vtot1)**(1._real_8/3._real_8)
          ENDDO
          ! Compute C6
          c6fac=1.2089117_real_8*SQRT(REAL(vdwwfi%nelpwf,kind=real_8))! bugfix
          DO icfr=1,nfrags(ipx)
             DO jcfr=icfr,nfrags(ipx)
                DO l=1,icontfragw(icfr,ipx)
                   DO m=1,icontfragw(jcfr,ipx)
                      seffi=spr(l,icfr,ipx)*fove(icfr)
                      seffj=spr(m,jcfr,ipx)*fove(jcfr)
                      seffi3=seffi*seffi*seffi
                      seffj3=seffj*seffj*seffj
                      seffi32=SQRT(seffi3)
                      seffj32=SQRT(seffj3)
                      c6=c6fac*seffi3*seffj3/(seffj32+seffi32)
                      c6all(ifragw(m,jcfr,ipx),ifragw(l,icfr,ipx),ipx)=c6
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
          DEALLOCATE(fove,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
       ELSE
          ! Keep using already computed C6ALL
       ENDIF
    ELSE IF (vdwwfi%iswitchvdw.EQ.1) THEN
       ! vdw-WF
       IF (vdwwfl%twannup) THEN
          c6fac=0.5_real_8*SQRT(REAL(vdwwfi%nelpwf,kind=real_8))/(3._real_8**1.25_real_8)
          DO icfr=1,nfrags(ipx)
             DO jcfr=icfr,nfrags(ipx)
                DO l=1,icontfragw(icfr,ipx)
                   DO m=1,icontfragw(jcfr,ipx)
                      si3=spr(l,icfr,ipx)*spr(l,icfr,ipx)*spr(l,icfr,ipx)
                      sj3=spr(m,jcfr,ipx)*spr(m,jcfr,ipx)*spr(m,jcfr,ipx)
                      si3d2=SQRT(si3)
                      c6=c6fac*si3d2*sj3*doubles(spr(l,icfr,ipx),spr(m,jcfr,ipx))
                      c6all(ifragw(m,jcfr,ipx),ifragw(l,icfr,ipx),ipx)=c6
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ELSE
          ! Keep using already computed C6ALL
       ENDIF
    ELSE
       CALL stopgm(procedureN,'WRONG ISWITCHVDW !',& 
            __LINE__,__FILE__)
    ENDIF
    ! Print computed C6 coefficients 
    IF (paral%io_parent.AND.vdwwfl%tpc6) THEN
       WRITE(6,'(/,4(1X,A),4X,A,3X,4(1X,A),4X,A,3X)')&
            'IFRAG','IWANN','IFRAG','IWANN','C6(A.U.)',&
            'IFRAG','IWANN','IFRAG','IWANN','C6(A.U.)'
       DO icfr=1,nfrags(ipx)
          DO jcfr=icfr,nfrags(ipx)
             DO l=1,icontfragw(icfr,ipx)
                DO m=1,icontfragw(jcfr,ipx),2
                   IF (m.EQ.icontfragw(jcfr,ipx)) THEN
                      WRITE(6,'(4I6,F12.5)')&
                           icfr,l,jcfr,m,c6all(ifragw(m,jcfr,ipx),ifragw(l,icfr,ipx),ipx)
                   ELSE
                      WRITE(6,'(4I6,F12.5,3X,4I6,F12.5)')&
                           icfr,l,jcfr,m,c6all(ifragw(m,jcfr,ipx),ifragw(l,icfr,ipx),ipx),&
                           icfr,l,jcfr,m+1,c6all(ifragw(m+1,jcfr,ipx),ifragw(l,icfr,ipx),ipx)
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Compute the van der Waals corrections
    IF (paral%io_parent.AND.vdwwfl%tpforce) THEN
       ALLOCATE(fionwf(3,maxsys%nax,maxsys%nsx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(fionwf)!,3*maxsys%nax*maxsys%nsx)
       CALL dcopy(3*maxsys%nax*maxsys%nsx,fion,1,fionwf,1)
    ENDIF
    IF (paral%io_parent) THEN
       DO icfr=1,nfrags(ipx)
          DO jcfr=icfr,nfrags(ipx)
             DO l=1,icontfragw(icfr,ipx)
                DO m=1,icontfragw(jcfr,ipx)
                   xlm_=rwfc(1,l,icfr,ipx)-rwfc(1,m,jcfr,ipx)
                   ylm_=rwfc(2,l,icfr,ipx)-rwfc(2,m,jcfr,ipx)
                   zlm_=rwfc(3,l,icfr,ipx)-rwfc(3,m,jcfr,ipx)
                   CALL pbc(xlm_,ylm_,zlm_,xlm,ylm,zlm,1,parm%apbc,parm%ibrav)
                   DO nx=-vdwi%nxvdw,vdwi%nxvdw
                      DO ny=-vdwi%nyvdw,vdwi%nyvdw
                         DO nz=-vdwi%nzvdw,vdwi%nzvdw
                            ! Exclude self-interactions: better here than icfr=1,vdwwfi%multifrag-1
                            IF (icfr==jcfr.AND.nx*nx+ny*ny+nz*nz==0) CYCLE
                            rlp(1)=xlm+nx*metr_com%ht(1,1)+ny*metr_com%ht(2,1)+nz*metr_com%ht(3,1)
                            rlp(2)=ylm+nx*metr_com%ht(1,2)+ny*metr_com%ht(2,2)+nz*metr_com%ht(3,2)
                            rlp(3)=zlm+nx*metr_com%ht(1,3)+ny*metr_com%ht(2,3)+nz*metr_com%ht(3,3)
                            distn=SQRT(rlp(1)*rlp(1)+rlp(2)*rlp(2)+rlp(3)*rlp(3))
                            distn6=distn*distn*distn*distn*distn*distn
                            c6vdw=-c6all(ifragw(m,jcfr,ipx),ifragw(l,icfr,ipx),ipx)
                            IF (distn.LT.reps) THEN
                               efac=0.0_real_8
                               ffac=0.0_real_8
                            ELSE
                               IF (vdwwfl%tdampda) THEN
                                  efac=c6vdw/distn6&
                                       +potda(distn,vdwwfi%nelpwf,spr(l,icfr,ipx),spr(m,jcfr,ipx))
                                  ffac=(-6._real_8*c6vdw/(distn6*distn)&
                                       +fpotda(distn,vdwwfi%nelpwf,spr(l,icfr,ipx),spr(m,jcfr,ipx)))/distn
                               ELSE
                                  efac=c6vdw&
                                       *POT59(DISTN,vdwwfr%a6,SPR(L,ICFR,ipx),SPR(M,JCFR,ipx),&
                                       vdwwfi%iswitchvdw)/DISTN6
                                  ffac=c6vdw&
                                       *FPOT59(DISTN,vdwwfr%a6,SPR(L,ICFR,ipx),SPR(M,JCFR,ipx),&
                                       vdwwfi%iswitchvdw)/DISTN
                               ENDIF
                            ENDIF
                            ! Energy
                            vdwr%evdw=vdwr%evdw+efac
                            ! Ionic forces
                            ip=icontfragwi(l,icfr,ipx)
                            DO i=1,ip
                               in=iwfcref(i,l,icfr,ipx)
                               ia=iatpt(1,in)
                               is=iatpt(2,in)
                               ! fion(1,ia,is)=fion(1,ia,is)-ffac*rlp(1)/REAL(ip,kind=real_8)
                               ! fion(2,ia,is)=fion(2,ia,is)-ffac*rlp(2)/REAL(ip,kind=real_8)
                               ! fion(3,ia,is)=fion(3,ia,is)-ffac*rlp(3)/REAL(ip,kind=real_8)
                               fion(1,ia,is)=fion(1,ia,is)-ffac*rlp(1)*wwfcref(i,l,icfr,ipx)
                               fion(2,ia,is)=fion(2,ia,is)-ffac*rlp(2)*wwfcref(i,l,icfr,ipx)
                               fion(3,ia,is)=fion(3,ia,is)-ffac*rlp(3)*wwfcref(i,l,icfr,ipx)
                            ENDDO
                            ip=icontfragwi(m,jcfr,ipx)
                            DO i=1,ip
                               in=iwfcref(i,m,jcfr,ipx)
                               ia=iatpt(1,in)
                               is=iatpt(2,in)
                               ! fion(1,ia,is)=fion(1,ia,is)+ffac*rlp(1)/REAL(ip,kind=real_8)
                               ! fion(2,ia,is)=fion(2,ia,is)+ffac*rlp(2)/REAL(ip,kind=real_8)
                               ! fion(3,ia,is)=fion(3,ia,is)+ffac*rlp(3)/REAL(ip,kind=real_8)
                               fion(1,ia,is)=fion(1,ia,is)+ffac*rlp(1)*wwfcref(i,m,jcfr,ipx)
                               fion(2,ia,is)=fion(2,ia,is)+ffac*rlp(2)*wwfcref(i,m,jcfr,ipx)
                               fion(3,ia,is)=fion(3,ia,is)+ffac*rlp(3)*wwfcref(i,m,jcfr,ipx)
                            ENDDO
                            ! Stress tensor
                            IF (.NOT.(ANY(alpha==0).OR.ANY(beta==0))) THEN
                               vdwr%devdw(1)=vdwr%devdw(1)&
                                    +FFAC*RLP(ALPHA(1))*RLP(BETA(1))
                               vdwr%devdw(2)=vdwr%devdw(2)&
                                    +FFAC*RLP(ALPHA(2))*RLP(BETA(2))
                               vdwr%devdw(3)=vdwr%devdw(3)&
                                    +FFAC*RLP(ALPHA(3))*RLP(BETA(3))
                               vdwr%devdw(4)=vdwr%devdw(4)&
                                    +FFAC*RLP(ALPHA(4))*RLP(BETA(4))
                               vdwr%devdw(5)=vdwr%devdw(5)&
                                    +FFAC*RLP(ALPHA(5))*RLP(BETA(5))
                               vdwr%devdw(6)=vdwr%devdw(6)&
                                    +FFAC*RLP(ALPHA(6))*RLP(BETA(6))
                            ENDIF
                         ENDDO
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    IF (paral%io_parent.AND.vdwwfl%tpforce) THEN
       CALL daxpy(3*maxsys%nax*maxsys%nsx,-1.0_real_8,fion,1,fionwf,1)
       CALL dscal(3*maxsys%nax*maxsys%nsx,-1.0_real_8,fionwf,1)
       WRITE(6,'(/,1X,A)') 'VAN DER WAALS CORRECTIONS TO IONIC FORCES'
       CALL wrgeof(tau0,fionwf)
       DEALLOCATE(fionwf,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Compute the displacement of ions w.r.t. reference coordinates
    dxmax=0.0_real_8
    dxnorm=0.0_real_8
    !$omp parallel do private(I,IA,IS,DR1,DR2,DR3,DATAUF) &
    !$omp  reduction(+:DXNORM) reduction(MAX:DXMAX)
    DO i=1,ions1%nat
       ia=iatpt(1,i)
       is=iatpt(2,i)
       dr1=tau0(1,ia,is)-tauref(1,ia,is,ipx)
       dr2=tau0(2,ia,is)-tauref(2,ia,is,ipx)
       dr3=tau0(3,ia,is)-tauref(3,ia,is,ipx)
       datauf=SQRT(dr1*dr1+dr2*dr2+dr3*dr3)
       dxnorm=dxnorm+datauf
       IF (datauf.GT.dxmax) dxmax=datauf
    ENDDO
    dxnorm=dxnorm/REAL(ions1%nat,kind=real_8)
    ! Set TWANNUP=.TRUE. if displacement of ions exceeds the tolerance
    vdwwfl%twannup=.FALSE.
    IF (dxmax.GT.vdwwfr%tolwann) vdwwfl%twannup=.TRUE.
    IF (paral%io_parent.AND.vdwwfl%tpinfo) THEN
       IF (vdwwfl%twannup) THEN
          WRITE(6,'(/,T3,A,E9.3,A,E9.3,A,/)')&
               'DXMAX= ',dxmax,' DXNORM= ',dxnorm,&
               '; RECOMPUTE WANNIER FUNCTIONS'
       ELSE
          WRITE(6,'(/,T3,A,E9.3,A,E9.3,A,/)')&
               'DXMAX= ',dxmax,' DXNORM= ',dxnorm,&
               '; EXTRAPOLATE WANNIER CENTERS'
       ENDIF
    ENDIF
9999 CONTINUE
    CALL mp_bcast(vdwwfl%twannup,parai%io_source,parai%cp_grp)
    twannupx(ipx)=vdwwfl%twannup

    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE vdw_wf
  ! ==================================================================
  REAL(real_8) FUNCTION fpotda(r,nel,s1,s2)
    ! ==--------------------------------------------------------------==
    IMPLICIT NONE
    ! Arguments
    INTEGER :: nel
    REAL(real_8) :: r,s1,s2
    ! Local variables
    REAL(real_8) :: q,ss,f1,f2,es
    ! ==--------------------------------------------------------------==
    q=REAL(nel,KIND=real_8)
    ss=s1*s1+s2*s2
    f1=s1*s2/ss
    f2=r*r/ss
    es=EXP(-1.5_real_8*f2)
    fpotda=-4._real_8*q**2*f1**3*(1._real_8+3._real_8*f2)*es/(r*r)
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION fpotda
  ! ==================================================================
  REAL(real_8) FUNCTION potda(r,nel,s1,s2)
    ! ==--------------------------------------------------------------==
    IMPLICIT NONE
    ! Arguments
    INTEGER :: nel
    REAL(real_8) :: r,s1,s2
    ! Local variables
    REAL(real_8) :: q,ss,f1,f2,es
    ! ==--------------------------------------------------------------==
    q=REAL(nel,KIND=real_8)
    ss=s1*s1+s2*s2
    f1=s1*s2/ss
    f2=r*r/ss
    es=EXP(-1.5_real_8*f2)
    potda=4._real_8*q**2*f1**3*es/r
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION potda
  ! ==================================================================
  REAL(real_8) FUNCTION fpot59(r,a,s1,s2,iswitchvdw)
    ! ==--------------------------------------------------------------==
    IMPLICIT NONE
    ! Arguments
    INTEGER :: iswitchvdw
    REAL(real_8) :: a,r,s1,s2
    ! Local variables
    REAL(real_8) :: ai,aj,ar,beta,bt,cobt,sibt
    REAL(real_8) :: dfonds,dlogso,dsondr,dsondrpar,dsondbeta
    REAL(real_8) :: ebt,embt,es,facpar,fdamp,rS,soij,tau
    ! ==--------------------------------------------------------------==
    IF (a.GT.1.e-5_real_8) THEN
       ! 
       IF (iswitchvdw.EQ.1) THEN
          rS=SQRT(3._real_8)*((0.769_real_8+0.5_real_8*LOG(s1))*s1&
               +(0.769_real_8+0.5_real_8*LOG(S2))*S2)
       ELSEIF (iswitchvdw.EQ.2) THEN
          rS=1.309_real_8*(s1+s2)
       ENDIF
       es=EXP(-a*(r/rS-1._real_8))
       fdamp=1._real_8/(1._real_8+es)
       facpar=-6._real_8+a*r/rS*fdamp*es
       fpot59=fdamp/(r**7)*facpar
    ELSE
       ! 
       ! Damping function after Mol. Phys. 107, 999 (2009):
       ! 
       IF (iswitchvdw.EQ.2) THEN
          rS=1.309_real_8*(s1+s2)
          es=EXP(-a*(r/rS-1._real_8))
          fdamp=1._real_8/(1._real_8+es)
          facpar=-6._real_8+a*r/rS*fdamp*es
          fpot59=fdamp/(r**7)*facpar
       ELSEIF (iswitchvdw.EQ.1) THEN
          ai=SQRT(3._real_8)/s1
          aj=SQRT(3._real_8)/s2
          beta=0.5_real_8*(ai+aj)*r
          tau=(ai-aj)/(ai+aj)
          IF (ABS(tau).LT.1.e-05_real_8) THEN
             ar=ai*r
             soij=EXP(-ar)*(1._real_8+ar+ar*ar/3._real_8)
             dsondr=-ai*ai/3._real_8*EXP(-ar)*r*(1._real_8+ar)
          ELSE
             bt=beta*tau
             cobt=COSH(bt)
             sibt=SINH(bt)
             soij=cobt*beta/tau+(1._real_8+beta-1._real_8/tau/tau)*sibt
             soij=soij*(ai*aj)**1.5_real_8*r**3/beta**4/tau*EXP(-beta)
             ebt=EXP(bt)
             embt=EXP(-bt)
             dsondrpar=3._real_8*soij/r
             dsondbeta=ebt*(beta*(1._real_8+tau)+1._real_8+tau)+&
                  EMBT*(BETA*(TAU-1._real_8)-1._real_8+TAU)
             dsondbeta=dsondbeta*(ai*aj)**1.5_real_8*r**3/beta**4&
                  /TAU*EXP(-BETA)
             dsondbeta=-(beta+4._real_8/beta)*soij+0.5_real_8*dsondbeta
             dsondr=dsondrpar+dsondbeta*(ai+aj)*0.5_real_8
          ENDIF
          dlogso=LOG(ABS(soij))
          fdamp=1._real_8-soij*soij*(1._real_8-2._real_8*dlogso*(1._real_8-dlogso))
          dfonds=-4._real_8*soij*dlogso*dlogso
          fpot59=(-6._real_8*fdamp+r*dfonds*dsondr)/(r**7)
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION fpot59
  ! ==================================================================
  REAL(real_8) FUNCTION pot59(r,a,s1,s2,iswitchvdw)
    ! ==--------------------------------------------------------------==
    ! ***-1: f(R)=1/[1+exp(-a*(R/Rs-1))]     
    IMPLICIT NONE
    ! Arguments
    INTEGER :: iswitchvdw
    REAL(real_8) :: a,r,s1,s2
    ! Local variables
    REAL(real_8) :: ai,aj,ar,beta,cobt,sibt,es,rS,soij,tau
    REAL(real_8) :: dlogso,root3
    ! ==--------------------------------------------------------------==
    root3=SQRT(3._real_8)
    ! 
    IF (a.GT.1.e-5_real_8) THEN
       IF (iswitchvdw.EQ.1) THEN
          rS=root3*((0.769_real_8+0.5_real_8*LOG(s1))*s1&
               +(0.769_real_8+0.5_real_8*LOG(S2))*S2)
       ELSE IF (iswitchvdw.EQ.2) THEN
          rS=1.309_real_8*(s1+s2)
       ENDIF
       es=EXP(-a*(r/rS-1._real_8))
       pot59=1._real_8/(1._real_8+es)
    ELSE
       ! 
       ! Damping function after Mol. Phys. 107, 999 (2009):
       ! 
       IF (iswitchvdw.EQ.2) THEN
          rS=1.309_real_8*(s1+s2)
          es=EXP(-a*(r/rS-1._real_8))
          pot59=1._real_8/(1._real_8+es)
       ELSE IF (iswitchvdw.EQ.1) THEN
          ai=root3/s1
          aj=root3/s2
          beta=0.5_real_8*(ai+aj)*r
          tau=(ai-aj)/(ai+aj)
          IF (ABS(tau).LT.1.e-05_real_8) THEN
             ar=ai*r
             soij=EXP(-ar)*(1._real_8+ar+ar*ar/3.0_real_8)
          ELSE
             cobt=COSH(beta*tau)
             sibt=SINH(beta*tau)
             soij=cobt*beta/tau+(1._real_8+beta-1._real_8/tau/tau)*sibt
             soij=soij*(ai*aj)**1.5_real_8*r**3/beta**4/tau*EXP(-beta)
          ENDIF
          dlogso=LOG(ABS(soij))
          pot59=1._real_8-soij*soij*(1._real_8-2._real_8*dlogso*(1._real_8-dlogso))
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION pot59
  ! ==================================================================
  FUNCTION doubles(si,sj)
    ! ==--------------------------------------------------------------==
    ! == This routine gives an estimate of the C6 Van der Waals
    ! == coefficient using the information relative to the
    ! == Wannier functions and solving the original 6-D
    ! == integral :
    ! == const.* int int dr dr' (w(r)*w(r'))/(w(r)+w(r'))/|r-r'|**6
    ! == by using spherical coordinates and thus reducing to a
    ! == 2-D integral, which is evaluated by Simpson method.
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: si, sj, doubles

    INTEGER, PARAMETER                       :: npoints = 51

    INTEGER                                  :: i, imax, j
    REAL(real_8)                             :: a, beta, frac23, frac43, hi, &
                                                hj, sum, wsum, xc, yc
    REAL(real_8), DIMENSION(npoints)         :: wi, wj, x, y

! ==--------------------------------------------------------------==
! 
! define cutoffs :
! 

    xc=2.307978_real_8+1.5_real_8*LOG(si)
    yc=2.307978_real_8+1.5_real_8*LOG(sj)
    ! 
    frac43=4._real_8/3._real_8
    frac23=2._real_8/3._real_8
    beta=(si/sj)**1.5_real_8
    a=1._real_8/beta
    hi=xc/REAL(npoints-1,kind=real_8)
    hj=yc/REAL(npoints-1,kind=real_8)
    x(1)=0._real_8
    y(1)=0._real_8
    wi(1)=hi/3._real_8
    wj(1)=hj/3._real_8
    wsum=0.0_real_8
    !$omp parallel do private(i) reduction(+:wsum)
    DO i=2,npoints-1
       x(i)=REAL(i-1,kind=real_8)*hi
       y(i)=REAL(i-1,kind=real_8)*hj
       IF (MOD(i,2).EQ.0) THEN
          wi(i)=hi*frac43
          wj(i)=hj*frac43
       ELSE
          wi(i)=hi*frac23
          wj(i)=hj*frac23
       ENDIF
       wsum=wsum+wi(i)+wj(i)
    ENDDO
    x(npoints)=REAL(npoints-1,kind=real_8)*hi
    y(npoints)=REAL(npoints-1,kind=real_8)*hj
    wi(npoints)=hi/3._real_8
    wj(npoints)=hj/3._real_8
    wsum=wsum+wi(1)+wj(1)+wi(npoints)+wj(npoints)
    ! dbg   WRITE(6,*) ' Weights sum check = ',wsum
    CALL set_ptdist(npoints,1,parai%nproc,imax)
    sum=0.0_real_8
    !$omp parallel do private(i,j) reduction(+:sum) __COLLAPSE2
    !   DO i=1,npoints
    DO i=npt12(parai%mepos,1),npt12(parai%mepos,2)
       DO j=1,npoints
          sum=sum+wi(i)*wj(j)*effe(x(i),y(j),a)
       ENDDO
    ENDDO
    CALL mp_sum(sum,parai%allgrp)
    doubles=sum
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION doubles
  ! ==================================================================
  FUNCTION effe(x,y,a)
    ! ==--------------------------------------------------------------==
    ! 
    ! decaying exponentials:
    REAL(real_8)                             :: x, y, a, effe

    REAL(real_8)                             :: ex, ey

! ==--------------------------------------------------------------==

    ex=EXP(-x)
    ey=EXP(-y)
    effe=x*x*y*y*ex*ey/(a*ex+ey)
    ! 
    ! gaussians:
    ! REAL*8 x,y,a,EFFE,x2,y2,ex2,ey2
    ! x2=x*x
    ! y2=y*y
    ! ex2=exp(-x2)
    ! ey2=exp(-y2)
    ! EFFE=x2*y2*ex2*ey2/(ex2+ey2)
    ! 
    ! spheres:
    ! REAL*8 x,y,a,EFFE
    ! REAL*8 pi,s,beta,s3d2,w,wx,wy,wsum
    ! pi=3.1415926535_real_8
    ! s=SQRT(3._real_8)
    ! beta=100.
    ! s3d2=s*SQRT(s)
    ! w=0.5*SQRT(3./pi)/s3d2
    ! wx=0.
    ! wy=0.
    ! IF(x.le.s) wx=w
    ! IF(y.le.s) wy=w
    ! c      wx=w/(exp(beta*(x-s))+1.)
    ! c      wy=w/(exp(beta*(y-s))+1.)
    ! wsum=wx+wy
    ! IF(wsum.eq.0.) THEN
    ! EFFE=0.
    ! ELSE
    ! EFFE=3.*SQRT(pi)*x*x*y*y*wx*wy/wsum
    ! ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION effe
  ! ==================================================================
  SUBROUTINE volwan(rw,sw,i1,iin,volume,vfree)
    ! ==--------------------------------------------------------------==
    ! == I1  = number of centers in the fragment with rw and sw       ==
    ! == IIN = center number for overlap calculation                  ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rw(3,*), sw(*)
    INTEGER                                  :: i1, iin
    REAL(real_8)                             :: volume, vfree

    INTEGER, PARAMETER                       :: npoints = 51

    INTEGER                                  :: imax, inters, la, nx, ny, nz
    REAL(real_8)                             :: dens, difac, distn, distp, &
                                                dmax, dmez, dv, dx, dxmez, &
                                                ff, r(3), rwr(3,i1), s12, &
                                                swr(i1), wt1, wt2

! ==--------------------------------------------------------------==
! Inters contains information on intersections

    inters=0
    ! Scan fragments to find overlap with iin and copy in reduced vectors
    DO la=1,i1
       IF (iin.NE.la) THEN
          s12=sw(iin)+sw(la)     ! sum the spreads
          distn=dist_pbc(rw(:,iin),rw(:,la),.TRUE.)! iin-la distance
          IF (distn.LT.s12) THEN  ! intersection -> save neighbors
             inters=inters+1    ! intersection number
             rwr(1,inters)=rw(1,la)! spread & coord of intersecting guys
             rwr(2,inters)=rw(2,la)
             rwr(3,inters)=rw(3,la)
             swr(inters)=sw(la)
          ENDIF
       ENDIF
    ENDDO
    ! Integral
    dmez=sw(iin)      ! cube having side=2*spread centered in rwfc1
    dmax=dmez*2.0_real_8
    dx=dmax/REAL(npoints,kind=real_8)
    dxmez=dx*0.5_real_8    ! half lattice size
    dv=dx*dx*dx
    ! 
    volume=0.0_real_8
    vfree=0.0_real_8
    CALL set_ptdist(npoints,1,parai%nproc,imax)
    !$omp parallel do private(NX,NY,NZ,LA,R,DISTN,DENS,DIFAC,DISTP &
    !$omp  ,WT1,WT2,FF) reduction(+:VOLUME,VFREE)
    !   DO nx=1,npoints
    DO nx=npt12(parai%mepos,1),npt12(parai%mepos,2)
       r(1)=REAL(nx-1,kind=real_8)*dx-dmez+dxmez+rw(1,iin)
       DO ny=1,npoints
          r(2)=REAL(ny-1,kind=real_8)*dx-dmez+dxmez+rw(2,iin)
          DO nz=1,npoints
             r(3)=REAL(nz-1,kind=real_8)*dx-dmez+dxmez+rw(3,iin)
             distn=dist_pbc(rw(:,iin),r(:),.TRUE.)! distance iin-r
             dens=EXP(-2.0_real_8*SQRT(3._real_8)*distn/dmez)
             IF (distn.LE.vdwwfr%tolref*dmez) THEN! inside the volume wan-iin
                difac=1.0_real_8       ! weight=1 if no overlap
                ! Account only of centers having intersection with the center under consideration
                DO la=1,inters    ! overlap re-weighting
                   distp=dist_pbc(rwr(:,la),r(:),.TRUE.)! distance from r
                   IF (distp.LE.vdwwfr%tolref*swr(la)) difac=difac+1.0_real_8
                ENDDO
                wt1=1.0_real_8/REAL(difac,kind=real_8)! weight the volume divided among fragments
                wt2=wt1*wt1
                ! c              FF=DENS*DISTN*DISTN*DISTN !additional weight tkatchenko-overlap
                ff=1.0_real_8
                volume=volume+wt2*dv*ff
                vfree=vfree+wt1*dv*ff
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    CALL mp_sum(volume,parai%allgrp)
    CALL mp_sum(vfree,parai%allgrp)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE volwan
  ! ==================================================================
  SUBROUTINE fragment_zlevel(nstate,rwann,swann,tauref,icontfragw,&
       icontfragwi,iwfcref,ifragw,wwfcref,rwfc,spr,taufrag)
    ! ==--------------------------------------------------------------==
    ! ==      CONSTRUCTS FRAGMENTS BY ADOPTING Z-LEVEL CRITERION      ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    REAL(real_8)                             :: rwann(:,:), swann(:), &
                                                tauref(:,:,:)
    INTEGER                                  :: icontfragw(:), &
                                                icontfragwi(:,:), &
                                                iwfcref(:,:,:), ifragw(:,:)
    REAL(real_8)                             :: wwfcref(:,:,:), rwfc(:,:,:), &
                                                spr(:,:), taufrag(:,:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'fragment_zlevel'

    INTEGER                                  :: i, ia, icontfragi(2), ierr, &
                                                ipar, is, jcfr, jpar
    INTEGER, ALLOCATABLE                     :: ind(:)
    REAL(real_8)                             :: zlm
    REAL(real_8), ALLOCATABLE                :: dcwfc(:)

! ==--------------------------------------------------------------==

    ALLOCATE(dcwfc(ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(ind(ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    icontfragi(1)=0
    icontfragi(2)=0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          zlm=tauref(3,ia,is)
          IF (zlm.LT.vdwwfr%zlevel) THEN
             icontfragi(1)=icontfragi(1)+1
          ELSE
             icontfragi(2)=icontfragi(2)+1
          ENDIF
       ENDDO
    ENDDO
    ! Identify WFs included in each fragment
    icontfragw(1)=0
    icontfragw(2)=0
    CALL zeroing(icontfragwi(:,1:2))!,nwfcx*2)
    DO jcfr=1,nstate
       zlm=rwann(3,jcfr)
       IF (zlm.LT.vdwwfr%zlevel) THEN
          icontfragw(1)=icontfragw(1)+1
          ipar=icontfragw(1)
          rwfc(1,ipar,1)=rwann(1,jcfr)
          rwfc(2,ipar,1)=rwann(2,jcfr)
          rwfc(3,ipar,1)=rwann(3,jcfr)
          spr(ipar,1)=swann(jcfr)
          ifragw(ipar,1)=jcfr
          DO i=1,ions1%nat
             ia=iatpt(1,i)
             is=iatpt(2,i)
             dcwfc(i)=dist_pbc(rwann(:,jcfr),tauref(:,ia,is),.TRUE.)
          ENDDO
          CALL sort2(dcwfc,ions1%nat,ind)
          icontfragwi(ipar,1)=icontfragwi(ipar,1)+1
          jpar=icontfragwi(ipar,1)
          iwfcref(jpar,ipar,1)=ind(1)
          DO i=2,ions1%nat
             IF (dcwfc(i).LE.swann(jcfr)*vdwwfr%tolref) THEN
                icontfragwi(ipar,1)=icontfragwi(ipar,1)+1
                jpar=icontfragwi(ipar,1)
                IF (jpar.GT.natwfcx)&
                     CALL stopgm(procedureN,'NATWFCX TOO SMALL ! ',& 
                     __LINE__,__FILE__)
                iwfcref(jpar,ipar,1)=ind(i)
             ENDIF
          ENDDO
       ELSE
          icontfragw(2)=icontfragw(2)+1
          ipar=icontfragw(2)
          rwfc(1,ipar,2)=rwann(1,jcfr)
          rwfc(2,ipar,2)=rwann(2,jcfr)
          rwfc(3,ipar,2)=rwann(3,jcfr)
          spr(ipar,2)=swann(jcfr)
          ifragw(ipar,2)=jcfr
          DO i=1,ions1%nat
             ia=iatpt(1,i)
             is=iatpt(2,i)
             dcwfc(i)=dist_pbc(rwann(:,jcfr),tauref(:,ia,is),.TRUE.)
          ENDDO
          CALL sort2(dcwfc,ions1%nat,ind)
          icontfragwi(ipar,2)=icontfragwi(ipar,2)+1
          jpar=icontfragwi(ipar,2)
          iwfcref(jpar,ipar,2)=ind(1)
          DO i=2,ions1%nat
             IF (dcwfc(i).LE.swann(jcfr)*vdwwfr%tolref) THEN
                icontfragwi(ipar,2)=icontfragwi(ipar,2)+1
                jpar=icontfragwi(ipar,2)
                IF (jpar.GT.natwfcx)&
                     CALL stopgm(procedureN,'NATWFCX TOO SMALL ! ',& 
                     __LINE__,__FILE__)
                iwfcref(jpar,ipar,2)=ind(i)
             ENDIF
          ENDDO
       ENDIF
    ENDDO
    ! Compute reference coordinates of each WFC
    ! DO icfr=1,2
    !    DO jcfr=1,icontfragw(icfr)
    !       temp1=0.0_real_8
    !       temp2=0.0_real_8
    !       temp3=0.0_real_8
    !       ipar=icontfragwi(jcfr,icfr)
    !       DO i=1,ipar
    !          in=iwfcref(i,jcfr,icfr)
    !          ia=iatpt(1,in)
    !          is=iatpt(2,in)
    !          temp1=temp1+tauref(1,ia,is)
    !          temp2=temp2+tauref(2,ia,is)
    !          temp3=temp3+tauref(3,ia,is)
    !       ENDDO
    !       taufrag(1,jcfr,icfr)=temp1/REAL(ipar,kind=real_8)
    !       taufrag(2,jcfr,icfr)=temp2/REAL(ipar,kind=real_8)
    !       taufrag(3,jcfr,icfr)=temp3/REAL(ipar,kind=real_8)
    !    ENDDO
    ! ENDDO
    CALL set_taufrag(vdwwfi%multifrag,icontfragw,icontfragwi,iwfcref,&
         tauref,rwfc,wwfcref,taufrag)
    ! Print information of fragments
    CALL print_assignment(vdwwfi%multifrag,nwfcx,natwfcx,&
         icontfragi,icontfragw,icontfragwi,&
         iwfcref,rwfc,spr,taufrag,&
         vdwwfl%tpinfo,vdwwfl%tpfrag)
    ! ==--------------------------------------------------------------==
    DEALLOCATE(dcwfc,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(ind,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE fragment_zlevel
  ! ==================================================================
  SUBROUTINE fragment_radius(nstate,rwann,swann,tauref,icontfragw,&
       icontfragwi,iwfcref,ifragw,wwfcref,rwfc,spr,taufrag)
    ! ==--------------------------------------------------------------==
    ! ==      CONSTRUCTS FRAGMENTS BY ADOPTING RADIUS CRITERION       ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    REAL(real_8)                             :: rwann(:,:), swann(:), &
                                                tauref(:,:,:)
    INTEGER                                  :: icontfragw(:), &
                                                icontfragwi(:,:), &
                                                iwfcref(:,:,:), ifragw(:,:)
    REAL(real_8)                             :: wwfcref(:,:,:), rwfc(:,:,:), &
                                                spr(:,:), taufrag(:,:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'fragment_radius'

    INTEGER                                  :: i, ia, icfr, &
                                                icontfragi(vdwwfi%multifrag), &
                                                ierr, in, ipar, is, j, jcfr, &
                                                jpar, m
    INTEGER, ALLOCATABLE                     :: ind(:)
    REAL(real_8)                             :: distfr
    REAL(real_8), ALLOCATABLE                :: dcwfc(:)

! ==--------------------------------------------------------------==

    ALLOCATE(dcwfc(ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(ind(ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(icontfragwi)!,nwfcx*vdwwfi%multifrag)
    DO icfr=1,vdwwfi%multifrag
       icontfragi(icfr)=0
       icontfragw(icfr)=0
       in=ifragdata(icfr)
       m=iatpt(1,in)
       j=iatpt(2,in)
       ! Count the number of ions included in each fragment
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             distfr=dist_pbc(tauref(:,ia,is),tauref(:,m,j),.TRUE.)
             IF (distfr.LE.radfrag(icfr)) THEN
                icontfragi(icfr)=icontfragi(icfr)+1
             ENDIF
          ENDDO
       ENDDO
       ! Identify WFs included in each fragment
       DO jcfr=1,nstate
          distfr=dist_pbc(rwann(:,jcfr),tauref(:,m,j),.TRUE.)
          IF (distfr.LE.radfrag(icfr)) THEN
             icontfragw(icfr)=icontfragw(icfr)+1
             ipar=icontfragw(icfr)
             rwfc(1,ipar,icfr)=rwann(1,jcfr)
             rwfc(2,ipar,icfr)=rwann(2,jcfr)
             rwfc(3,ipar,icfr)=rwann(3,jcfr)
             spr(ipar,icfr)=swann(jcfr)
             ifragw(ipar,icfr)=jcfr
             DO i=1,ions1%nat
                ia=iatpt(1,i)
                is=iatpt(2,i)
                dcwfc(i)=dist_pbc(rwann(:,jcfr),tauref(:,ia,is),.TRUE.)
             ENDDO
             CALL sort2(dcwfc,ions1%nat,ind)
             icontfragwi(ipar,icfr)=icontfragwi(ipar,icfr)+1
             jpar=icontfragwi(ipar,icfr)
             iwfcref(jpar,ipar,icfr)=ind(1)
             DO i=2,ions1%nat
                IF (dcwfc(i).LE.swann(jcfr)*vdwwfr%tolref) THEN
                   icontfragwi(ipar,icfr)=icontfragwi(ipar,icfr)+1
                   jpar=icontfragwi(ipar,icfr)
                   IF (jpar.GT.natwfcx)&
                        CALL stopgm(procedureN,'NATWFCX TOO SMALL ! ',& 
                        __LINE__,__FILE__)
                   iwfcref(jpar,ipar,icfr)=ind(i)
                ENDIF
             ENDDO
          ENDIF
       ENDDO
    ENDDO
    ! Compute reference coordinates of each WFC
    ! DO icfr=1,vdwwfi%multifrag
    !    DO jcfr=1,icontfragw(icfr)
    !       temp1=0.0_real_8
    !       temp2=0.0_real_8
    !       temp3=0.0_real_8
    !       ipar=icontfragwi(jcfr,icfr)
    !       DO i=1,ipar
    !          in=iwfcref(i,jcfr,icfr)
    !          ia=iatpt(1,in)
    !          is=iatpt(2,in)
    !          temp1=temp1+tauref(1,ia,is)
    !          temp2=temp2+tauref(2,ia,is)
    !          temp3=temp3+tauref(3,ia,is)
    !       ENDDO
    !       taufrag(1,jcfr,icfr)=temp1/REAL(ipar,kind=real_8)
    !       taufrag(2,jcfr,icfr)=temp2/REAL(ipar,kind=real_8)
    !       taufrag(3,jcfr,icfr)=temp3/REAL(ipar,kind=real_8)
    !    ENDDO
    ! ENDDO
    CALL set_taufrag(vdwwfi%multifrag,icontfragw,icontfragwi,iwfcref,&
         tauref,rwfc,wwfcref,taufrag)
    ! Print information of fragments
    CALL print_assignment(vdwwfi%multifrag,nwfcx,natwfcx,&
         icontfragi,icontfragw,icontfragwi,&
         iwfcref,rwfc,spr,taufrag,&
         vdwwfl%tpinfo,vdwwfl%tpfrag)
    ! ==--------------------------------------------------------------==
    DEALLOCATE(dcwfc,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(ind,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE fragment_radius
  ! ==================================================================
  SUBROUTINE fragment_bond(nstate,rwann,swann,tauref,icontfragw,&
       icontfragwi,iwfcref,ifragw,wwfcref,rwfc,spr,taufrag)
    ! ==--------------------------------------------------------------==
    ! ==      CONSTRUCTS FRAGMENTS BY ADOPTING BOND CRITERION         ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    REAL(real_8)                             :: rwann(:,:), swann(:), &
                                                tauref(:,:,:)
    INTEGER                                  :: icontfragw(:), &
                                                icontfragwi(:,:), &
                                                iwfcref(:,:,:), ifragw(:,:)
    REAL(real_8)                             :: wwfcref(:,:,:), rwfc(:,:,:), &
                                                spr(:,:), taufrag(:,:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'fragment_bond'
    INTEGER, PARAMETER                       :: maxl = 2

    INTEGER                                  :: i, ia, icfr, ierr, ii, ipar, &
                                                is, j, jcfr, k, l, nmpa
    INTEGER, ALLOCATABLE                     :: ibind(:), icontfragi(:), &
                                                imol(:), ind(:), iz(:)
    REAL(real_8)                             :: cdij, dij
    REAL(real_8), ALLOCATABLE                :: c(:,:), dcwfc(:)

! ==--------------------------------------------------------------==

    nmpa=ions1%nat*(ions1%nat+1)/2
    ALLOCATE(iz(ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(c(3,ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(dcwfc(ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(ibind(nmpa),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(imol(ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(ind(ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(icontfragi(ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(ibind)!,nmpa)
    CALL zeroing(imol)!,ions1%nat)
    DO i=1,ions1%nat
       ia=iatpt(1,i)
       is=iatpt(2,i)
       iz(i)=ions0%iatyp(is)
       c(1,i)=tauref(1,ia,is)
       c(2,i)=tauref(2,ia,is)
       c(3,i)=tauref(3,ia,is)
    ENDDO
    ! Normal bonds
    DO i=1,ions1%nat
       DO j=i+1,ions1%nat
          ii=i+j*(j-1)/2
          dij=dist_pbc(c(:,i),c(:,j),.TRUE.)
          cdij=vdwwfr%xmfacwf*(covrad(iz(i))+covrad(iz(j)))*fbohr
          IF (cdij.GE.dij) THEN
             ibind(ii)=1
          ELSE
             ibind(ii)=0
          ENDIF
       ENDDO
    ENDDO
    ! Add or delete bonds
    IF (vdwwfi%nboadwf.GT.0) THEN
       DO k=1,vdwwfi%nboadwf
          i=boadwf(1,k)
          j=boadwf(2,k)
          l=boadwf(3,k)
          IF (i.GT.j) THEN
             ii=j+i*(i-1)/2
          ELSE
             ii=i+j*(j-1)/2
          ENDIF
          IF (l.LT.0) ibind(ii)=0
          IF (l.GT.0) ibind(ii)=1
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Search for fragments
    vdwwfi%multifrag=0
10  CONTINUE
    DO i=1,ions1%nat
       IF (imol(i).EQ.0) THEN
          vdwwfi%multifrag=vdwwfi%multifrag+1
          imol(i)=vdwwfi%multifrag
          GOTO 20
       ENDIF
    ENDDO
    GOTO 30
20  CONTINUE
    DO l=1,maxl
       ! Scan in ascending order
       DO i=1,ions1%nat
          IF (imol(i).EQ.vdwwfi%multifrag) THEN
             DO j=1,i-1
                ii=j+i*(i-1)/2
                IF (ibind(ii).EQ.1) imol(j)=vdwwfi%multifrag
             ENDDO
             DO j=i+1,ions1%nat
                ii=i+j*(j-1)/2
                IF (ibind(ii).EQ.1) imol(j)=vdwwfi%multifrag
             ENDDO
          ENDIF
       ENDDO
       ! Rescan in descending order
       DO i=ions1%nat,1,-1
          IF (imol(i).EQ.vdwwfi%multifrag) THEN
             DO j=1,i-1
                ii=j+i*(i-1)/2
                IF (ibind(ii).EQ.1) imol(j)=vdwwfi%multifrag
             ENDDO
             DO j=i+1,ions1%nat
                ii=i+j*(j-1)/2
                IF (ibind(ii).EQ.1) imol(j)=vdwwfi%multifrag
             ENDDO
          ENDIF
       ENDDO
    ENDDO
    GOTO 10
30  CONTINUE
    CALL zeroing(icontfragi(1:vdwwfi%multifrag))!,vdwwfi%multifrag)
    DO i=1,ions1%nat
       icfr=imol(i)
       icontfragi(icfr)=icontfragi(icfr)+1
    ENDDO
    ! ==--------------------------------------------------------------==
    ! Assign each wfc to one of fragments
    CALL zeroing(icontfragw(1:vdwwfi%multifrag))!,vdwwfi%multifrag)
    CALL zeroing(icontfragwi)!,nwfcx*vdwwfi%multifrag)
    DO j=1,nstate
       DO i=1,ions1%nat
          dcwfc(i)=dist_pbc(rwann(:,j),c(:,i),.TRUE.)
       ENDDO
       CALL sort2(dcwfc,ions1%nat,ind)
       ipar=ind(1)
       icfr=imol(ipar)
       icontfragw(icfr)=icontfragw(icfr)+1
       jcfr=icontfragw(icfr)
       rwfc(1,jcfr,icfr)=rwann(1,j)
       rwfc(2,jcfr,icfr)=rwann(2,j)
       rwfc(3,jcfr,icfr)=rwann(3,j)
       spr(jcfr,icfr)=swann(j)
       ifragw(jcfr,icfr)=j
       icontfragwi(jcfr,icfr)=icontfragwi(jcfr,icfr)+1
       ipar=icontfragwi(jcfr,icfr)
       iwfcref(ipar,jcfr,icfr)=ind(1)
       DO i=2,ions1%nat
          IF (dcwfc(i).LE.swann(j)*vdwwfr%tolref) THEN
             icontfragwi(jcfr,icfr)=icontfragwi(jcfr,icfr)+1
             ipar=icontfragwi(jcfr,icfr)
             IF (ipar.GT.natwfcx)&
                  CALL stopgm(procedureN,'NATWFCX TOO SMALL ! ',& 
                  __LINE__,__FILE__)
             iwfcref(ipar,jcfr,icfr)=ind(i)
          ENDIF
       ENDDO
    ENDDO
    ! Compute reference coordinates of each WFC
    ! DO icfr=1,vdwwfi%multifrag
    !    DO jcfr=1,icontfragw(icfr)
    !       temp1=0.0_real_8
    !       temp2=0.0_real_8
    !       temp3=0.0_real_8
    !       ipar=icontfragwi(jcfr,icfr)
    !       DO i=1,ipar
    !          in=iwfcref(i,jcfr,icfr)
    !          ia=iatpt(1,in)
    !          is=iatpt(2,in)
    !          temp1=temp1+tauref(1,ia,is)
    !          temp2=temp2+tauref(2,ia,is)
    !          temp3=temp3+tauref(3,ia,is)
    !       ENDDO
    !       taufrag(1,jcfr,icfr)=temp1/REAL(ipar,kind=real_8)
    !       taufrag(2,jcfr,icfr)=temp2/REAL(ipar,kind=real_8)
    !       taufrag(3,jcfr,icfr)=temp3/REAL(ipar,kind=real_8)
    !    ENDDO
    ! ENDDO
    CALL set_taufrag(vdwwfi%multifrag,icontfragw,icontfragwi,iwfcref,&
         tauref,rwfc,wwfcref,taufrag)
    ! Print out information of fragments
    CALL print_assignment(vdwwfi%multifrag,nwfcx,natwfcx,&
         icontfragi,icontfragw,icontfragwi,&
         iwfcref,rwfc,spr,taufrag,&
         vdwwfl%tpinfo,vdwwfl%tpfrag)
    ! ==--------------------------------------------------------------==
    DEALLOCATE(iz,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(c,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(dcwfc,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(ibind,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(imol,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(ind,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(icontfragi,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE fragment_bond
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE print_assignment(mfrag,nwfcx,natx,ifragi,ifragw,&
       ifragwi,iref,rwfc,spr,tauf,&
       tpinfo,tpfrag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: mfrag, nwfcx, natx, &
                                                ifragi(*), ifragw(*), &
                                                ifragwi(nwfcx,*), &
                                                iref(natx,nwfcx,*)
    REAL(real_8)                             :: rwfc(3,nwfcx,*), &
                                                spr(nwfcx,*), tauf(3,nwfcx,*)
    LOGICAL                                  :: tpinfo, tpfrag

    INTEGER                                  :: i, j, k

! ==--------------------------------------------------------------==

    IF (paral%io_parent.AND.tpinfo) THEN
       WRITE(6,'(/,1X,23("*"),A,23("*"))')&
            ' VDW-WANNIER INITIALIZATION '
       WRITE(6,'(1X,A,T70,I6)') 'TOTAL NUMBER OF FRAGMENTS:',mfrag
       DO i=1,mfrag
          WRITE(6,'(1X,A,1X,I6,T21,A,T40,I6,T47,A,T53,A,T59,I6,T66,A)')&
               'FRAGMENT NR.',i,':',ifragi(i),'IONS','AND',ifragw(i),&
               'WF CENTERS'
          IF (tpfrag) THEN
             WRITE(6,'(1X,A,T16,A,T35,A,T48,A,T56,A)')&
                  'IWANN','CENTER(BOHR)','SPREAD(BOHR)','R',&
                  'REFERENCE(BOHR)'
             DO j=1,ifragw(i)
                WRITE(6,'(I6,4F9.4,I6,3F9.4)')&
                     j,(rwfc(k,j,i),k=1,3),spr(j,i),ifragwi(j,i),&
                     (tauf(k,j,i),k=1,3)
             ENDDO
          ENDIF
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE print_assignment
  ! ==================================================================
  SUBROUTINE wwannier(nstate,tau0,rwann,swann)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    REAL(real_8)                             :: tau0(:,:,:), rwann(:,:), &
                                                swann(:)

    CHARACTER(len=100)                       :: filen
    INTEGER                                  :: i, ia, ipx, is, j
    LOGICAL                                  :: erread

! ==--------------------------------------------------------------==

    filen='WANNIER_RESTART'
    IF (cntl%tpath.OR.tmw) THEN
       IF (cntl%tpath) THEN
          ipx=ipcurr
       ELSE IF (tmw) THEN
          ipx=mwi%walker_id
       ENDIF
       CALL mw_filename('WANNIER_RESTART_',filen,ipx)
    ENDIF
    CALL fileopen(83,filen,fo_def,erread)
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          WRITE(83,*) (tau0(i,ia,is),i=1,3)
       ENDDO
    ENDDO
    WRITE(83,*) nstate
    DO i=1,nstate
       WRITE(83,*) (rwann(j,i),j=1,3),swann(i)
    ENDDO
    CALL fileclose(83)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE wwannier
  ! ==================================================================
  SUBROUTINE rwannier(nstate,tau0,rwann,swann,trwannc)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    REAL(real_8)                             :: tau0(:,:,:), rwann(:,:), &
                                                swann(:)
    LOGICAL                                  :: trwannc

    CHARACTER(*), PARAMETER                  :: procedureN = 'rwannier'

    CHARACTER(len=100)                       :: filen
    INTEGER                                  :: i, ia, ipx, is, j, nwfc
    LOGICAL                                  :: erread

! ==--------------------------------------------------------------==

    filen='WANNIER_RESTART'
    IF (cntl%tpath.OR.tmw) THEN
       IF (cntl%tpath) THEN
          ipx=ipcurr
       ELSE IF (tmw) THEN
          ipx=mwi%walker_id
       ENDIF
       CALL mw_filename('WANNIER_RESTART_',filen,ipx)
    ENDIF
    IF (trwannc) THEN
       CALL fileopen(83,filen,fo_old,erread)
       IF (.NOT.erread) THEN
          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                READ(83,err=99,END=99,fmt=*) (tau0(i,ia,is),i=1,3)
             ENDDO
          ENDDO
          READ(83,err=99,END=99,fmt=*) nwfc
          IF (nwfc.NE.nstate)&
               CALL stopgm(procedureN,'WRONG NUMBER OF WANNIER CENTERS ! ',& 
               __LINE__,__FILE__)
          DO i=1,nwfc
             READ(83,err=99,END=99,fmt=*) (rwann(j,i),j=1,3),swann(i)
          ENDDO
          CALL fileclose(83)
          trwannc=.TRUE.
       ELSE
          trwannc=.FALSE.
       ENDIF
       GOTO 100
    ENDIF
99  CONTINUE
    trwannc=.FALSE.
100 CONTINUE
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rwannier
  ! ==================================================================
  SUBROUTINE set_ptdist(nstate,nblock,my_nproc,nbmax)
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
       npt12(ip-1,1)=NINT(xsnow)+1
       npt12(ip-1,2)=NINT(xsaim)
       IF (NINT(xsaim).GT.nstate) THEN
          npt12(ip-1,2)=nstate
       ENDIF
       IF (NINT(xsnow).GT.nstate) THEN
          npt12(ip-1,1)=nstate+1
       ENDIF
       xsnow = xsaim
    ENDDO
    DO ip=0,my_nproc-1
       nx=npt12(ip,2)-npt12(ip,1)+1
       nbmax=MAX(nbmax,nx)
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE set_ptdist
  ! ==================================================================
  SUBROUTINE set_taufrag(mfrag,icontfragw,icontfragwi,iwfcref,&
       tauref,rwfc,wwfcref,taufrag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: mfrag, icontfragw(:), &
                                                icontfragwi(:,:), &
                                                iwfcref(:,:,:)
    REAL(real_8)                             :: tauref(:,:,:), rwfc(:,:,:), &
                                                wwfcref(:,:,:), taufrag(:,:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'set_taufrag'
    REAL(real_8), PARAMETER                  :: tols = 1.e-12_real_8

    INTEGER                                  :: i, ia, icfr, ierr, in, info, &
                                                ipar, is, jcfr, lwork, &
                                                maxdim, maxipar, mindim, nrank
    REAL(real_8)                             :: temp1, temp2, temp3, w, a_(3)
    REAL(real_8), ALLOCATABLE                :: a(:,:), b(:), s(:), work(:)

! ==--------------------------------------------------------------==

    IF (vdwwfl%treffit) THEN
       ! LLS fit to WF centers 
       maxipar=0
       DO icfr=1,mfrag
          DO jcfr=1,icontfragw(icfr)
             ipar=icontfragwi(jcfr,icfr)
             maxipar=MAX(maxipar,ipar)
          ENDDO
       ENDDO
       ALLOCATE(a(4,maxipar+1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       maxdim=MAX(4,maxipar+1)
       ALLOCATE(b(maxdim),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       mindim=MIN(4,maxipar+1)
       ALLOCATE(s(mindim),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       lwork=3*mindim+MAX(2*mindim,maxdim)
       ALLOCATE(work(lwork),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       DO icfr=1,mfrag
          DO jcfr=1,icontfragw(icfr)
             ipar=icontfragwi(jcfr,icfr)
             IF (ipar==1) THEN
                in=iwfcref(1,jcfr,icfr)
                ia=iatpt(1,in)
                is=iatpt(2,in)
                taufrag(:,jcfr,icfr)=tauref(:,ia,is)
                wwfcref(1,jcfr,icfr)=1.0_real_8
             ELSE
                DO i=1,ipar
                   in=iwfcref(i,jcfr,icfr)
                   ia=iatpt(1,in)
                   is=iatpt(2,in)
                   a_(1:3)=tauref(1:3,ia,is)-rwfc(1:3,jcfr,icfr)
                   CALL pbc(a_(1),a_(2),a_(3),a(1,i),a(2,i),a(3,i),1,&
                        parm%apbc,parm%ibrav)
                   a(4,i)=1.0_real_8
                ENDDO
                b(1:3)=0.0_real_8
                b(4)=1.0_real_8
                CALL dgelss(4,ipar,1,a,4,b,maxdim,s,tols,nrank,work,lwork,info)
                IF (info/=0) CALL stopgm(procedureN,'DGELSS ENDED WITH INFO.NE.0',&
                     __LINE__,__FILE__)
                temp1=0.0_real_8
                temp2=0.0_real_8
                temp3=0.0_real_8
                DO i=1,ipar
                   in=iwfcref(i,jcfr,icfr)
                   ia=iatpt(1,in)
                   is=iatpt(2,in)
                   temp1=temp1+tauref(1,ia,is)*b(i)
                   temp2=temp2+tauref(2,ia,is)*b(i)
                   temp3=temp3+tauref(3,ia,is)*b(i)
                   wwfcref(i,jcfr,icfr)=b(i)
                ENDDO
                taufrag(1,jcfr,icfr)=temp1
                taufrag(2,jcfr,icfr)=temp2
                taufrag(3,jcfr,icfr)=temp3
             ENDIF
          ENDDO
       ENDDO
       DEALLOCATE(a,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(b,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(s,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(work,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ELSE
       ! atom/bond centers
       DO icfr=1,mfrag
          DO jcfr=1,icontfragw(icfr)
             temp1=0.0_real_8
             temp2=0.0_real_8
             temp3=0.0_real_8
             ipar=icontfragwi(jcfr,icfr)
             w=1.0_real_8/REAL(ipar,kind=real_8)
             DO i=1,ipar
                in=iwfcref(i,jcfr,icfr)
                ia=iatpt(1,in)
                is=iatpt(2,in)
                temp1=temp1+tauref(1,ia,is)
                temp2=temp2+tauref(2,ia,is)
                temp3=temp3+tauref(3,ia,is)
                wwfcref(i,jcfr,icfr)=w
             ENDDO
             taufrag(1,jcfr,icfr)=temp1*w
             taufrag(2,jcfr,icfr)=temp2*w
             taufrag(3,jcfr,icfr)=temp3*w
          ENDDO
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE set_taufrag
  ! ==================================================================

END MODULE vdw_utils
