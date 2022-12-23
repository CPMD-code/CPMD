#ifdef __SR11000
!option OPT(O(3))
#endif

MODULE setbasis_utils
  USE adat,                            ONLY: elem,&
                                             nelcon
  USE atom,                            ONLY: ecpfiles
  USE atwf,                            ONLY: &
       atrg, atrg_epr, atwf_mod, atwfr, atwfr_epr, atwp, atwr, cat, catord, &
       cdftoc, cdftocl, m1shlx, natsave
  USE cdftmod,                         ONLY: cdftlog
  USE cnst,                            ONLY: fpi
  USE cppt,                            ONLY: gk,&
                                             hg,&
                                             inyh
  USE error_handling,                  ONLY: stopgm
  USE fitpack_utils,                   ONLY: curv1,&
                                             curv2
  USE gvec,                            ONLY: gvec_com
  USE inscan_utils,                    ONLY: inscan
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE lsfbtr_utils,                    ONLY: lsfbtr
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sync
  USE mw,                              ONLY: tmw
  USE parac,                           ONLY: parai,&
                                             paral
  USE pimd,                            ONLY: grandparent,&
                                             supergroup,&
                                             supersource
  USE pslo,                            ONLY: pslo_com
  USE qspl,                            ONLY: ggng,&
                                             ggnh,&
                                             nsplpo
  USE readsr_utils,                    ONLY: readsi,&
                                             xstring
  USE recpnew_utils,                   ONLY: ckgrid,&
                                             get_pplib,&
                                             tgrid
  USE response_pmod,                   ONLY: epr_options,&
                                             nmr_para,&
                                             response1
  USE sfac,                            ONLY: ei1,&
                                             ei2,&
                                             ei3,&
                                             eigr,&
                                             eigrb
  USE sphe,                            ONLY: gcutka
  USE system,                          ONLY: cntl,&
                                             maxsys,&
                                             ncpw,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: icopy
  USE vdbp,                            ONLY: maxsp,&
                                             ncpr1,&
                                             vdb_pawf,&
                                             vdb_r
  USE ylmr_utils,                      ONLY: ylmr
  USE zeroing_utils,                   ONLY: zeroing
#include "sizeof.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: setbasis
  PUBLIC :: loadc
  PUBLIC :: rtog
  PUBLIC :: storps

CONTAINS

  ! ==================================================================
  SUBROUTINE setbasis
    ! ==--------------------------------------------------------------==
    ! == SET ATOMIC BASIS         (call once)                         ==
    ! ==--------------------------------------------------------------==
    ! Variables
    CHARACTER(*), PARAMETER                  :: procedureN = 'setbasis'
    CHARACTER(len=1), DIMENSION(0:3), &
      PARAMETER                              :: lq = (/'S','P','D','F'/)
    INTEGER, PARAMETER                       :: ifnum = 21 

    INTEGER                                  :: i, iatom, ibtype(maxsp), &
                                                ierr, iogroup, iosource, ir, &
                                                is, isub, l, m, ms, n, ne(4,7)
    INTEGER, SAVE                            :: ifirst = 0
    LOGICAL                                  :: ionode
    REAL(real_8)                             :: dclog, valch, zc

! ==--------------------------------------------------------------==

    IF (ifirst.NE.0) RETURN
    ! ==--------------------------------------------------------------==
    CALL tiset(procedureN,isub)
    IF (cntl%tpath.OR.tmw) THEN
       ionode=grandparent
       iosource=supersource
       iogroup=supergroup
       ! dsebasti
    ELSEIF (response1%tnmr) THEN
       ionode=nmr_para%nmr_superparent
       iosource=nmr_para%nmr_supersource
       iogroup=nmr_para%nmr_supergroup
       ! dsebasti
    ELSE
       ionode=paral%io_parent
       iosource=parai%io_source
       iogroup=parai%cp_grp
    ENDIF
    ALLOCATE(atwfr(maxsys%mmaxx,m1shlx,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(atrg(maxsys%mmaxx,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    IF (epr_options%tepr_hyp) THEN
       ALLOCATE(atwfr_epr(maxsys%mmaxx,m1shlx,maxsys%nsx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(atrg_epr(maxsys%mmaxx,maxsys%nsx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ifirst=1
    IF (ionode) THEN
       WRITE(6,'(/,A)') ' GENERATE ATOMIC BASIS SET'        
       ! ==--------------------------------------------------------------==
       ! SLATER MINIMAL BASIS
       DO is=1,ions1%nsp
          ibtype(is)=1
          iatom=ions0%iatyp(is)
          valch=ions0%zv(is)
          zc=REAL(iatom,kind=real_8)-valch
          CALL icopy(28,nelcon(1,1,iatom),1,ne,1)
          ms=0
          m=0
          DO n=1,7
             DO l=1,4
                m=m+ne(l,n)
                IF (m.GT.NINT(zc).AND.ne(l,n).NE.0) THEN
                   ms=ms+1
                   atwf_mod%nqsto(ms,is)=n
                   atwf_mod%lshell(ms,is)=l-1
                   atwf_mod%stoexp(ms,is)=srules(iatom,ne,n,l-1)
                   atwf_mod%oc(ms,is)=ne(l,n)
                ENDIF
             ENDDO
          ENDDO
          atwf_mod%nshell(is)=ms
          atwr%meshat(is)=256
          dclog=1.049999881_real_8
          atwr%clogat(is)=LOG(dclog)
          atrg(1,is)=0.7142857e-03_real_8
          DO ir=2,atwr%meshat(is)
             atrg(ir,is)=dclog*atrg(ir-1,is)
          ENDDO
          DO i=1,atwf_mod%nshell(is)
             l=atwf_mod%lshell(i,is)
             CALL storps(atrg(1,is),atwfr(1,i,is),atwr%meshat(is),&
                  atwf_mod%nqsto(i,is),l,atwf_mod%stoexp(i,is))
          ENDDO
       ENDDO
       IF (cntl%cdft)THEN
          DO is=1,ions1%nsp
             CALL dcopy(atwf_mod%nshell(is),atwf_mod%oc(1,is),1,cdftoc%ocun(1,is),1)
          ENDDO
       ENDIF
       ! LOOK FOR SECTION BASIS ON INPUT FILE
       CALL read_basis(ibtype)
       ! Print Basis Set Info
       DO is=1,ions1%nsp
          iatom=ions0%iatyp(is)
          IF (ibtype(is).EQ.1) THEN
             WRITE(6,'(5X,A,8X,A)') elem%el(iatom),'SLATER ORBITALS'
             DO i=1,atwf_mod%nshell(is)
                WRITE(6,&
                     '(8X,I1,A1,8X,"ALPHA=",F9.4,T40,"OCCUPATION=",F5.2)')&
                     atwf_mod%nqsto(i,is),lq(atwf_mod%lshell(i,is)),atwf_mod%stoexp(i,is),atwf_mod%oc(i,is)
             ENDDO
          ELSEIF (ibtype(is).EQ.2) THEN
             WRITE(6,'(5X,A,8X,A)') elem%el(iatom),'NUMERICAL ORBITALS'
             DO i=1,atwf_mod%nshell(is)
                WRITE(6,&
                     '(8X,"L VALUE=",A,6X,T40,"OCCUPATION=",F5.2)')&
                     lq(atwf_mod%lshell(i,is)),atwf_mod%oc(i,is)
             ENDDO
          ELSEIF (ibtype(is).EQ.3) THEN
             WRITE(6,'(5X,A,8X,A)') elem%el(iatom),'PSEUDO ATOMIC ORBITALS'
             DO i=1,atwf_mod%nshell(is)
                WRITE(6,&
                     '(8X,"L VALUE=",A,6X,T40,"OCCUPATION=",F5.2)')&
                     lq(atwf_mod%lshell(i,is)),atwf_mod%oc(i,is)
             ENDDO
          ENDIF
       ENDDO
       WRITE(6,*)
    ENDIF
    ! Parallel
    CALL mp_sync(iogroup)
    CALL mp_bcast(atwf_mod%nshell,SIZE(atwf_mod%nshell),iosource,iogroup)
    CALL mp_bcast_byte(atwf_mod, size_in_bytes_of(atwf_mod),iosource,iogroup)
    CALL mp_bcast_byte(atwr, size_in_bytes_of(atwr),iosource,iogroup)
    CALL mp_bcast(atwr%clogat,SIZE(atwr%clogat),iosource,iogroup)
    CALL mp_bcast(atwfr,SIZE(atwfr),iosource,iogroup)
    CALL mp_bcast(atrg,SIZE(atrg),iosource,iogroup)
    IF (epr_options%tepr_hyp) THEN
       CALL mp_bcast(atwfr_epr,SIZE(atwfr_epr),iosource,iogroup)
       CALL mp_bcast(atrg_epr,SIZE(atrg_epr),iosource,iogroup)
       CALL mp_bcast(atwr%mesh_epr,SIZE(atwr%mesh_epr),iosource,iogroup)
    ENDIF
    IF (cntl%cdft)THEN
       CALL mp_bcast_byte(cdftoc, size_in_bytes_of(cdftoc),iosource,iogroup)
       CALL mp_bcast_byte(cdftocl, size_in_bytes_of(cdftocl),iosource,iogroup)
    ENDIF
    ! Parallel
    ! Initalization of NATTOT and M1SHL
    atwp%m1shl=0
    atwp%nattot=0
    atwp%numaormax=0
    DO is=1,ions1%nsp
       IF (atwp%m1shl.LT.atwf_mod%nshell(is)) atwp%m1shl=atwf_mod%nshell(is)
       atwf_mod%numaor(is)=0
       DO i=1,atwf_mod%nshell(is)
          l=atwf_mod%lshell(i,is)
          atwp%nattot=atwp%nattot+ions0%na(is)*(2*l+1)
          atwf_mod%numaor(is)=atwf_mod%numaor(is)+(2*l+1)
       ENDDO
       atwp%numaormax=MAX(atwp%numaormax,atwf_mod%numaor(is))
    ENDDO
    ALLOCATE(cat(nsplpo,2,atwp%m1shl,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    IF (cntl%tipao) THEN
       ALLOCATE(catord(m1shlx,maxsys%nax*maxsys%nsx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)! TODO check dimensions
    ENDIF
    ! Transform to G-Space
    CALL rtog(0)
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    RETURN
  END SUBROUTINE setbasis
  ! ==================================================================
  SUBROUTINE rtog(ity)
    ! ==--------------------------------------------------------------==
    ! == Transform atomic orbital in R-space to G-space               ==
    ! == Calculate CAT(NSPLPO,2,M1SHL,maxsys%nsx) for each species           ==
    ! ==--------------------------------------------------------------==
    ! 
    INTEGER                                  :: ity

    CHARACTER(*), PARAMETER                  :: procedureN = 'rtog'

    COMPLEX(real_8), ALLOCATABLE             :: work(:,:)
    INTEGER                                  :: ierr, il, ir, is, ishell, l, &
                                                mmax, n2, n22, nwork
    LOGICAL                                  :: ionode, saved
    REAL(real_8)                             :: disc, gmin, rmin, xmax
    REAL(real_8), ALLOCATABLE                :: bsint(:), fint(:), gg(:), &
                                                temp(:,:)

    atwf_mod%nbcut=ity
    IF (cntl%tpath) THEN
       ionode=grandparent
       ! dsebasti
    ELSEIF (response1%tnmr) THEN
       ionode=nmr_para%nmr_superparent
       ! dsebasti
    ELSE
       ionode=paral%io_parent
    ENDIF
    nwork=2*maxsys%mmaxx
    ! Allocation of local arrays
    ALLOCATE(work(nwork,5),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(fint(nwork),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(bsint(nwork),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(gg(nwork),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(temp(maxsys%mmaxx,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL zeroing(cat)!,2*nsplpo*atwp%m1shl*maxsys%nsx)
    DO is=1,ions1%nsp
       mmax=atwr%meshat(is)
       xmax=mmax
       n2=NINT(LOG(xmax)/LOG(2._real_8)+0.499999_real_8)
       n22=2**n2
       rmin=LOG(atrg(1,is))
       IF (ity.EQ.1) THEN
          gmin=LOG(SQRT(gvec_com%gcut)*parm%tpiba)-(mmax-1)*atwr%clogat(is)
       ELSE
          gmin=LOG(SQRT(gvec_com%gcutw+gcutka)*parm%tpiba)-(mmax-1)*atwr%clogat(is)
       ENDIF
       DO il=1,mmax
          gg(il)=(EXP(gmin+(il-1)*atwr%clogat(is))/parm%tpiba)**2
       ENDDO
       DO ishell=1,atwf_mod%nshell(is)
          saved=.FALSE.
          l=atwf_mod%lshell(ishell,is)
          CALL zeroing(fint(1:n22))!,n22)
          DO ir=1,mmax
             fint(ir)=fpi*atwfr(ir,ishell,is)/atrg(ir,is)
          ENDDO
          ! Fourier transformation
          CALL lsfbtr(fint,bsint,l,rmin,gmin,atwr%clogat(is),n2,saved,&
               work(1,1),work(1,2),work(1,3),work(1,5),nwork,disc)
          IF (disc.GT.1.e-5_real_8.AND.ionode) THEN
             IF (paral%io_parent)&
                  WRITE(6,*) ' WARNING! LSFBTR ACCURACY ONLY ',disC
          ENDIF
          IF (ity.EQ.1) THEN
             CALL tgrid(gg,mmax,ggnh,nsplpo,bsint,&
                  maxsys%mmaxx,temp(1,1),temp(1,2),temp(1,3))
             CALL dcopy(nsplpo,bsint(1),1,cat(1,1,ishell,is),1)
             IF (l.GT.0.AND.ggnh(1).LT.1.e-12_real_8) cat(1,1,ishell,is)=0.0_real_8
             CALL curv1(nsplpo,ggnh,cat(1,1,ishell,is),0.0_real_8,0.0_real_8,3,&
                  cat(1,2,ishell,is),temp,0.0_real_8,ierr)
          ELSE
             CALL tgrid(gg,mmax,ggng,nsplpo,bsint,&
                  maxsys%mmaxx,temp(1,1),temp(1,2),temp(1,3))
             CALL dcopy(nsplpo,bsint(1),1,cat(1,1,ishell,is),1)
             IF (l.GT.0.AND.ggng(1).LT.1.e-12_real_8) cat(1,1,ishell,is)=0.0_real_8
             CALL curv1(nsplpo,ggng,cat(1,1,ishell,is),0.0_real_8,0.0_real_8,3,&
                  cat(1,2,ishell,is),temp,0.0_real_8,ierr)
          ENDIF
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    ! Deallocation of local arrays
    DEALLOCATE(work,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(fint,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(bsint,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(gg,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(temp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rtog
  ! ==================================================================
  SUBROUTINE storps(r,rwfn,mesh,n,l,alpha)
    ! ==--------------------------------------------------------------==
    ! == CALCULATE R*STO ON A RADIAL GRID                             ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: mesh
    REAL(real_8)                             :: rwfn(mesh), r(mesh)
    INTEGER                                  :: n, l
    REAL(real_8)                             :: alpha

    REAL(real_8), DIMENSION(0:14), PARAMETER :: fac = (/ 1._real_8, 1._real_8,&
      2._real_8, 6._real_8, 24._real_8, 120._real_8, 720._real_8,5040._real_8,&
      40320._real_8, 362880._real_8, 36288.e2_real_8, 399168.e2_real_8,&
      4.790016e8_real_8, 6.2270208e9_real_8, 8.7178291e10_real_8/)
    REAL(real_8), PARAMETER                  :: xrmin = -400._real_8 

    INTEGER                                  :: ir
    REAL(real_8)                             :: an, xr

! ==--------------------------------------------------------------==

    an=SQRT((2._real_8*alpha)**(2*n+1)/fac(2*n))
    DO ir=1,mesh
       xr=-alpha*r(ir)
       IF (xr.LT.xrmin) THEN
          rwfn(ir) = 0._real_8
       ELSE
          rwfn(ir) = an*r(ir)**n * EXP(-alpha*r(ir))
       ENDIF
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE storps
  ! ==================================================================
  SUBROUTINE gtorps(r,rwfn,mesh,l,gexp,gcoe,nalp)
    ! ==--------------------------------------------------------------==
    ! == CALCULATE R*SUM(I)GTO_I ON A RADIAL GRID                     ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: mesh
    REAL(real_8)                             :: rwfn(mesh), r(mesh)
    INTEGER                                  :: l, nalp
    REAL(real_8)                             :: gcoe(nalp), gexp(nalp)

    REAL(real_8), DIMENSION(4), PARAMETER :: &
      dfac = (/ 1._real_8,3._real_8,15._real_8,105._real_8/)

    INTEGER                                  :: ic, ir
    REAL(real_8)                             :: an, pi

! ==--------------------------------------------------------------==

    pi=ACOS(-1._real_8)
    CALL zeroing(rwfn)!,mesh)
    DO ic=1,nalp
       an=2._real_8**(l+2)*dfac(l+1)**(-0.5_real_8)*(2._real_8*pi)**(-0.25_real_8)*&
            gexp(ic)**((2._real_8*l+3._real_8)/4._real_8)
       an=an*gcoe(ic)
       DO ir=1,mesh
          rwfn(ir) = rwfn(ir) +&
               an*r(ir)**(l+1) * EXP(-gexp(ic)*r(ir)*r(ir))
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gtorps
  ! ==================================================================
  FUNCTION srules(natom,ne,n,l)
    ! ==--------------------------------------------------------------==
    ! == Get Slater exponents from the Clementi-Raimondi table or     ==
    ! == Calculate Slater exponents by Slaters rule                   ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: natom, ne(4,7), n, l
    REAL(real_8)                             :: srules

    REAL(real_8), DIMENSION(2, 3, 18), PARAMETER :: cr = RESHAPE((/1.0_real_8,&
      0.0_real_8,0.0_real_8,0.0_real_8,0.0_real_8,0.0_real_8, 1.6875_real_8,&
      0.0_real_8,0.0_real_8,0.0_real_8,0.0_real_8,0.0_real_8,2.6906_real_8,&
      0.0_real_8,0.6396_real_8,0.0_real_8,0.0_real_8,0.0_real_8,3.6848_real_8,&
      0.0_real_8,0.9560_real_8,0.0_real_8,0.0_real_8,0.0_real_8,4.6795_real_8,&
      0.0_real_8,1.2881_real_8,1.2107_real_8,0.0_real_8,0.0_real_8,&
      5.6727_real_8,0.0_real_8,1.6083_real_8,1.5679_real_8,0.0_real_8,&
      0.0_real_8,6.6651_real_8,0.0_real_8,1.9237_real_8,1.9170_real_8,&
      0.0_real_8,0.0_real_8,7.6579_real_8,0.0_real_8,2.2458_real_8,&
      2.2266_real_8,0.0_real_8,0.0_real_8,8.6501_real_8,0.0_real_8,&
      2.5638_real_8,2.5500_real_8,0.0_real_8,0.0_real_8,9.6421_real_8,&
      0.0_real_8,2.8792_real_8,2.8792_real_8,0.0_real_8,0.0_real_8,&
      10.6259_real_8,0._real_8,3.2857_real_8,3.4009_real_8,0.8359_real_8,&
      0._real_8,11.6089_real_8,0.0_real_8,3.6960_real_8,3.9129_real_8,&
      1.1025_real_8,0.0_real_8,12.5910_real_8,0.0_real_8,4.1068_real_8,&
      4.4817_real_8,1.3724_real_8,1.3552_real_8,13.5724_real_8,0.0_real_8,&
      5.5100_real_8,4.9725_real_8,1.6344_real_8,1.4284_real_8,14.5578_real_8,&
      0.0_real_8,4.9125_real_8,5.4806_real_8,1.8806_real_8,1.6288_real_8,&
      15.5409_real_8,0.0_real_8,5.3144_real_8,5.9885_real_8,2.1223_real_8,&
      1.8273_real_8,16.5239_real_8,0.0_real_8,5.7152_real_8,6.4966_real_8,&
      2.3561_real_8,2.0387_real_8,17.5075_real_8,0.0_real_8,6.1152_real_8,&
      7.0041_real_8,2.5856_real_8,2.2547_real_8/), (/2,3,18/))
    REAL(real_8), DIMENSION(7), PARAMETER :: xns = (/1.0_real_8,2.0_real_8,&
      3.0_real_8,3.7_real_8,4.0_real_8,4.2_real_8,4.4_real_8/)

    INTEGER                                  :: i, l1, l2, m, m1, m2
    REAL(real_8)                             :: s

! ==--------------------------------------------------------------==
! Try the Clementi-Raimondi table

    IF (natom.LE.18.AND.l.LE.1) THEN
       srules=cr(l+1,n,natom)
       RETURN
    ENDIF
    ! Calculate the shielding
    s=0.0_real_8
    ! The complete shell
    l1=l+1
    IF (l1.EQ.1) l2=2
    IF (l1.EQ.2) l2=1
    IF (l1.EQ.3) l2=4
    IF (l1.EQ.4) l2=3
    ! Rule a) no contribution from shells further out
    ! Rule b) 0.35 (1s 0.3) from each other electron in the same shell
    IF (n.EQ.1) THEN
       m=ne(1,1)
       s=s+0.3_real_8*(m-1)
    ELSE
       m=ne(l1,n)+ne(l2,n)
       s=s+0.35_real_8*(m-1)
    ENDIF
    ! Rule c) if (s,p) shell 0.85 from each electron with n-1, and 1.0
    ! from all electrons further in
    IF (l1+l2.EQ.3) THEN
       IF (n.GT.1) THEN
          m1=ne(1,n-1)+ne(2,n-1)+ne(3,n-1)+ne(4,n-1)
          m2=0
          DO i=1,n-2
             m2=m2+ne(1,i)+ne(2,i)+ne(3,i)+ne(4,i)
          ENDDO
          s=s+0.85_real_8*m1+1.0_real_8*m2
       ENDIF
    ELSE
       ! Rule d) if (d,f) shell 1.0 from each electron inside
       m=0
       DO i=1,n-1
          m=m+ne(1,i)+ne(2,i)+ne(3,i)+ne(4,i)
       ENDDO
       s=s+1.0_real_8*m
    ENDIF
    ! Slater exponent is (Z-S)/NS
    srules = (REAL(natom,kind=real_8) - s)/xns(n)
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION srules
  ! ==================================================================
  SUBROUTINE read_basis(ibtype)
    ! ==--------------------------------------------------------------==
    ! ==  THIS ROUTINE READS THE SECTION &BASIS &END ON UNIT IUNIT    ==
    ! ==--------------------------------------------------------------==
    ! ==  THE INPUT HAS TO BE IN THE FOLLOWING FORMAT                 ==
    ! ==     &BASIS                                                   ==
    ! ==        BASIS ATOM TYPE 1                                     ==
    ! ==        ...                                                   ==
    ! ==        BASIS ATOM TYPE N                                     ==
    ! ==     &END                                                     ==
    ! ==--------------------------------------------------------------==
    ! ==  DEFINITION OF BASIS                                         ==
    ! ==                                                              ==
    ! ==  STO TYPE:                                                   ==
    ! ==  ---------                                                   ==
    ! ==      SLATER     nshell [OCCUPATION]                          ==
    ! ==       n   l     exp                                          ==
    ! ==       ...                                                    ==
    ! ==       n   l     exp                                          ==
    ! ==       [f1 f2 ... fn]                                         ==
    ! ==                                                              ==
    ! ==  NUMERICAL FUNCTIONS:                                        ==
    ! ==  --------------------                                        ==
    ! ==      *filename  nshell  [FORMAT=n]  [OCCUPATION]             ==
    ! ==      l1 ... ln                                               ==
    ! ==       [f1 f2 ... fn]                                         ==
    ! ==                                                              ==
    ! ==  GAUSSIAN BASIS FUNCTIONS:                                   ==
    ! ==  -------------------------                                   ==
    ! ==      *filename  nshell  GAUSSIAN [OCCUPATION]                ==
    ! ==      l1 ... ln                                               ==
    ! ==       [f1 f2 ... fn]                                         ==
    ! ==                                                              ==
    ! ==  PSEUDO ATOMIC ORBITALS:                                     ==
    ! ==  -----------------------                                     ==
    ! ==      PSEUDO AO  nshell  [OCCUPATION]                         ==
    ! ==      l1 ... ln     lx=-1 --> skip                            ==
    ! ==       [f1 f2 ... fn]                                         ==
    ! ==                                                              ==
    ! ==  SKIP ATOM TYPE:                                             ==
    ! ==  ---------------                                             ==
    ! ==      SKIP ATOM                                               ==
    ! ==                                                              ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ibtype(*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'read_basis'
    CHARACTER(len=1), DIMENSION(0:4), &
      PARAMETER                              :: shln = (/'S','P','D','F','G'/)
    INTEGER, PARAMETER                       :: ifnum = 21 

    CHARACTER(len=120)                       :: ecplib
    CHARACTER(len=200)                       :: filen, fnames
    CHARACTER(len=80)                        :: line
    INTEGER :: i, ia, iaa, ias, iatom, ie, ieold, ierr, iform, il, ins, iout, &
      ir, is, iunit, j, l, lenecp, lgmax, lsold(m1shlx), meshat2, nalp(4), &
      nfmax(4), nold, npoint, ntt(4)
    LOGICAL                                  :: erread, exists, toccup
    REAL(real_8)                             :: clogat2, dclog, gcoe(30,5,4), &
                                                gexp(30,4), ocold(m1shlx), &
                                                rp1, rp2
    REAL(real_8), ALLOCATABLE                :: atrg2(:), rp(:), temp(:,:)

! ==--------------------------------------------------------------==

    iunit=5
    ierr=inscan(iunit,'&BASIS')
    IF (ierr.NE.0) RETURN
    is=0
    ! ==--------------------------------------------------------------==
    ! Allocation of local variables
    ALLOCATE(atrg2(maxsys%mmaxx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(rp(maxsys%mmaxx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(temp(maxsys%mmaxx,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL get_pplib(ecplib,lenecp)
10  CONTINUE
    IF (paral%io_parent)&
         READ(iunit,err=99,END=99,fmt='(A80)') line
    IF (INDEX(line,'&END').NE.0) GOTO 100
    IF (INDEX(line,'SLATER').NE.0) THEN
       ! STO
       is=is+1
       nold=atwf_mod%nshell(is)
       CALL icopy(nold,atwf_mod%lshell(1,is),1,lsold,1)
       CALL dcopy(nold,atwf_mod%oc(1,is),1,ocold(1),1)
       IF (is.GT.ions1%nsp) CALL stopgm('READ_BASIS','NSP',& 
            __LINE__,__FILE__)
       toccup=.FALSE.
       IF (INDEX(line,'OCCUPATION').NE.0) toccup=.TRUE.
       ! Read STO
       CALL xstring(line,ia,ie)
       ! Read NSHELL(IS)
       CALL readsi(line,ie+1,iout,atwf_mod%nshell(is),erread)
       IF (erread.OR.atwf_mod%nshell(is).LE.0) THEN
          IF (paral%io_parent)&
               WRITE(6,*) 'READ_BASIS! STO FOR IS=',iS
          CALL stopgm('READ_BASIS','BAD VALUE FOR NSHELL',& 
               __LINE__,__FILE__)
       ENDIF
       DO i=1,atwf_mod%nshell(is)
          IF (paral%io_parent)&
               READ(iunit,END=99,err=99,fmt=*)&
               atwf_mod%nqsto(i,is),atwf_mod%lshell(i,is),atwf_mod%stoexp(i,is)
       ENDDO
       atwr%meshat(is)=256
       dclog=1.050_real_8
       atwr%clogat(is)=LOG(dclog)
       atrg(1,is)=0.7142857e-03_real_8
       DO ir=2,atwr%meshat(is)
          atrg(ir,is)=dclog*atrg(ir-1,is)
       ENDDO
       DO i=1,atwf_mod%nshell(is)
          l=atwf_mod%lshell(i,is)
          CALL storps(atrg(1,is),atwfr(1,i,is),atwr%meshat(is),&
               atwf_mod%nqsto(i,is),l,atwf_mod%stoexp(i,is))
       ENDDO
       IF (toccup) THEN
          IF (paral%io_parent)&
               READ(iunit,END=99,err=99,fmt=*) (atwf_mod%oc(i,is),i=1,atwf_mod%nshell(is))
          IF (cdftlog%thda)THEN
             cdftocl%ocset=.TRUE.
             IF (paral%io_parent)&
                  READ(iunit,END=99,err=99,fmt=*)&
                  (cdftoc%oc2(i,is),i=1,atwf_mod%nshell(is))
          ENDIF
          IF (cntl%cdft)THEN
             CALL filloc(nold,lsold,ocold,atwf_mod%nshell(is),&
                  atwf_mod%lshell(1,is),cdftoc%ocun(1,is))
          ENDIF
       ELSE
          CALL filloc(nold,lsold,ocold,atwf_mod%nshell(is),&
               atwf_mod%lshell(1,is),atwf_mod%oc(1,is))
       ENDIF
       GOTO 10
    ELSE IF (INDEX(line,'*').NE.0) THEN
       ! Numerical AO or Gaussian AO
       ias=INDEX(line,'*')
       line(ias:ias)=' '
       is=is+1
       ibtype(is)=2
       nold=atwf_mod%nshell(is)
       CALL icopy(nold,atwf_mod%lshell(1,is),1,lsold,1)
       CALL dcopy(nold,atwf_mod%oc(1,is),1,ocold(1),1)
       IF (is.GT.ions1%nsp) CALL stopgm('READ_BASIS','NSP',& 
            __LINE__,__FILE__)
       toccup=.FALSE.
       IF (INDEX(line,'OCCUPATION').NE.0) toccup=.TRUE.
       ! Read *filename
       CALL xstring(line,ia,ie)
       ! Read NSHELL(IS)
       CALL readsi(line,ie+1,iout,atwf_mod%nshell(is),erread)
       IF (erread.OR.atwf_mod%nshell(is).LE.0) THEN
          IF (paral%io_parent)&
               WRITE(6,*) 'READ_BASIS! ATOMIC ORBITAL FOR IS=',iS
          CALL stopgm('READ_BASIS','BAD VALUE FOR NSHELL',& 
               __LINE__,__FILE__)
       ENDIF
       IF (paral%io_parent)&
            READ(iunit,*) (atwf_mod%lshell(i,is),i=1,atwf_mod%nshell(is))
       filen=ecplib(1:lenecp)//line(ia:ie)
       IF (paral%io_parent)&
            INQUIRE(file=filen,exist=exists)
       IF (.NOT.exists) THEN
          IF (paral%io_parent)&
               WRITE(6,*) ' READ_BASIS: FILE ',filen,' NOT FOUND'
          CALL stopgm('READ_BASIS',' ',& 
               __LINE__,__FILE__)
       ENDIF
       IF (paral%io_parent)&
            OPEN(unit=ifnum,file=filen,status='UNKNOWN')
       IF (paral%io_parent)&
            REWIND(ifnum)
       IF (INDEX(line,'GAUSS').NE.0) THEN
          ! Gaussians
          IF (paral%io_parent)&
               READ(ifnum,*)
          IF (paral%io_parent)&
               READ(ifnum,*) lgmax
          DO il=1,lgmax
             IF (paral%io_parent)&
                  READ(ifnum,*)
             IF (paral%io_parent)&
                  READ(ifnum,*) nfmax(il),nalp(il)
             IF (paral%io_parent)&
                  READ(ifnum,*)
             IF (paral%io_parent)&
                  READ(ifnum,*) (gexp(i,il),i=1,nalp(il))
             DO j=1,nfmax(il)
                IF (paral%io_parent)&
                     READ(ifnum,*) (gcoe(i,j,il),i=1,nalp(il))
             ENDDO
          ENDDO
          atwr%meshat(is)=256
          dclog=1.050_real_8
          atwr%clogat(is)=LOG(dclog)
          atrg(1,is)=0.7142857e-03_real_8
          DO ir=2,atwr%meshat(is)
             atrg(ir,is)=dclog*atrg(ir-1,is)
          ENDDO
          DO i=1,4
             ntt(i)=0
          ENDDO
          DO i=1,atwf_mod%nshell(is)
             l=atwf_mod%lshell(i,is)+1
             ntt(l)=ntt(l)+1
             CALL gtorps(atrg(1,is),atwfr(1,i,is),atwr%meshat(is),&
                  l-1,gexp(1,l),gcoe(1,ntt(l),l),nalp(l))
          ENDDO
       ELSE
          ! Numerical
          iform=1
          IF (INDEX(line,'FORMAT=').NE.0) THEN
             iaa=INDEX(line,'FORMAT=')
             CALL readsi(line,iaa+7,iout,iform,erread)
             IF (erread) THEN
                CALL stopgm('READ_BASIS','ERROR WHEN READING FORMAT',& 
                     __LINE__,__FILE__)
             ENDIF
          ENDIF
          IF (iform.EQ.1) THEN
             IF (paral%io_parent)&
                  READ(ifnum,*) atwr%meshat(is),atwr%clogat(is)
             IF (atwr%clogat(is).GT.1.0_real_8) atwr%clogat(is)=LOG(atwr%clogat(is))
             DO ir=1,atwr%meshat(is)
                IF (paral%io_parent)&
                     READ(ifnum,*) atrg(ir,is),&
                     (atwfr(ir,il,is),il=1,atwf_mod%nshell(is))
             ENDDO
          ELSEIF (iform.EQ.2) THEN
             DO il=1,atwf_mod%nshell(is)
                IF (paral%io_parent)&
                     READ(ifnum) atwr%meshat(is),atwr%clogat(is),&
                     (atrg(ir,is),atwfr(ir,il,is),ir=1,atwr%meshat(is))
             ENDDO
             IF (atwr%clogat(is).GT.1.0_real_8) atwr%clogat(is)=LOG(atwr%clogat(is))
          ELSEIF (iform.EQ.3) THEN
             IF (paral%io_parent)&
                  READ(ifnum,*)
             DO il=1,atwf_mod%nshell(is)
                IF (paral%io_parent)&
                     READ(ifnum,*) atwr%meshat(is),atwr%clogat(is)
                DO ir=1,atwr%meshat(is)
                   IF (paral%io_parent)&
                        READ(ifnum,*) atrg(ir,is),atwfr(ir,il,is)
                ENDDO
             ENDDO
             IF (atwr%clogat(is).GT.1.0_real_8) atwr%clogat(is)=LOG(atwr%clogat(is))
          ELSE
             IF (paral%io_parent)&
                  WRITE(6,*) ' READ_BASIS: THIS FORMAT IS NOT PROGRAMMED'
             CALL stopgm('READ_BASIS',' ',& 
                  __LINE__,__FILE__)
          ENDIF
          IF (paral%io_parent)&
               CLOSE(ifnum)
          meshat2=256
          rp1=atrg(1,is)
          rp2=atrg(atwr%meshat(is),is)
          CALL ckgrid(rp1,rp2,atrg2,meshat2,clogat2)
          DO il=1,atwf_mod%nshell(is)
             CALL tgrid(atrg(1,is),atwr%meshat(is),atrg2,&
                  meshat2,atwfr(1,il,is),&
                  maxsys%mmaxx,temp(1,1),temp(1,2),temp(1,3))
          ENDDO
          CALL dcopy(meshat2,atrg2(1),1,atrg(1,is),1)
          atwr%clogat(is)=clogat2
          atwr%meshat(is)=meshat2
       ENDIF
       IF (toccup) THEN
          IF (paral%io_parent)&
               READ(iunit,END=99,err=99,fmt=*) (atwf_mod%oc(i,is),i=1,atwf_mod%nshell(is))
          IF (cdftlog%thda)THEN
             cdftocl%ocset=.TRUE.
             IF (paral%io_parent)&
                  READ(iunit,END=99,err=99,fmt=*)&
                  (cdftoc%oc2(i,is),i=1,atwf_mod%nshell(is))
          ENDIF
          IF (cntl%cdft)THEN
             CALL filloc(nold,lsold,ocold,atwf_mod%nshell(is),&
                  atwf_mod%lshell(1,is),cdftoc%ocun(1,is))
          ENDIF
       ELSE
          CALL filloc(nold,lsold,ocold,atwf_mod%nshell(is),&
               atwf_mod%lshell(1,is),atwf_mod%oc(1,is))
       ENDIF
       GOTO 10
    ELSE IF( (INDEX(line,'PSEUDO').NE.0).AND.&
         (INDEX(line,'AO').NE.0) ) THEN
       ! Pseudo AO
       is=is+1
       ibtype(is)=3
       nold=atwf_mod%nshell(is)
       CALL icopy(nold,atwf_mod%lshell(1,is),1,lsold,1)
       CALL dcopy(nold,atwf_mod%oc(1,is),1,ocold(1),1)
       IF (is.GT.ions1%nsp) CALL stopgm('READ_BASIS','NSP',& 
            __LINE__,__FILE__)
       toccup=.FALSE.
       IF (INDEX(line,'OCCUPATION').NE.0) toccup=.TRUE.
       ! Read PSEUDO
       CALL xstring(line,ia,ie)
       ieold=ie
       ! Read AO
       CALL xstring(line(ie+1:80),ia,ie)
       ie=ieold+ie+1-1
       ! Read NSHELL
       CALL readsi(line,ie+1,iout,atwf_mod%nshell(is),erread)
       IF (erread.OR.atwf_mod%nshell(is).LE.0) THEN
          IF (paral%io_parent)&
               WRITE(6,'(1X,A,I2)') 'READ_BASIS! PSEUDO AO FOR IS=',iS
          CALL stopgm('READ_BASIS','BAD VALUE FOR NSHELL',& 
               __LINE__,__FILE__)
       ENDIF
       IF (paral%io_parent)&
            READ(iunit,*) (atwf_mod%lshell(i,is),i=1,atwf_mod%nshell(is))
       iatom=ions0%iatyp(is)
       fnames=ecpfiles(is)
       CALL xstring(fnames,ia,ie)
       ! Read also for USPP
       IF (pslo_com%tvan(is))THEN
          npoint=ncpr1%meshva(is)-1
          DO ir=1,npoint
             atrg(ir,is)=vdb_r(ir+1,is)
             DO il=1,atwf_mod%nshell(is)
                atwfr(ir,il,is)=vdb_pawf(is,ir+1,il)
             ENDDO
          ENDDO
       ELSE
          IF (paral%io_parent)&
               OPEN(unit=ifnum,file=fnames(ia:ie),status='UNKNOWN')
          ierr=inscan(ifnum,'&WAVEFUNCTION')
          IF (ierr.NE.0)&
               CALL stopgm('READ_BASIS','SECTION WAVEFUNCTION',& 
               __LINE__,__FILE__)
          IF (paral%io_parent)&
               READ(ifnum,*) npoint
          DO ir=1,npoint
             IF (paral%io_parent)&
                  READ(ifnum,err=99,fmt=*) atrg(ir,is),&
                  (atwfr(ir,il,is),il=1,atwf_mod%nshell(is))
          ENDDO
       ENDIF
       IF (epr_options%tepr_hyp) THEN
          ! RDeclerck
          atwr%mesh_epr(is)=npoint
          atrg_epr(1,is)=1
          CALL dcopy(maxsys%mmaxx,atrg(1,is),1,&
               atrg_epr(1,is),1)
          CALL dcopy(maxsys%mmaxx*atwf_mod%nshell(is),atwfr(1,1,is),1,&
               atwfr_epr(1,1,is),1)
          ! RDeclerck
       ENDIF
       atwr%meshat(is)=256
       CALL dcopy(npoint,atrg(1,is),1,rp(1),1)
       rp1=rp(1)
       rp2=rp(npoint)
       CALL ckgrid(rp1,rp2,atrg(1,is),atwr%meshat(is),atwr%clogat(is))
       DO il=1,atwf_mod%nshell(is)
          CALL tgrid(rp,npoint,atrg(1,is),atwr%meshat(is),&
               atwfr(1,il,is),&
               maxsys%mmaxx,temp(1,1),temp(1,2),temp(1,3))
       ENDDO
       ins=0
       DO il=1,atwf_mod%nshell(is)
          IF (atwf_mod%lshell(il,is).GE.0) THEN
             ins=ins+1
             IF (ins.NE.il) THEN
                atwf_mod%lshell(ins,is)=atwf_mod%lshell(il,is)
                CALL dcopy(atwr%meshat(is),atwfr(1,il,is),1,&
                     atwfr(1,ins,is),1)
             ENDIF
          ENDIF
       ENDDO
       atwf_mod%nshell(is)=ins
       IF (paral%io_parent)&
            CLOSE(ifnum)
       IF (toccup) THEN
          IF (paral%io_parent)&
               READ(iunit,END=99,err=99,fmt=*) (atwf_mod%oc(i,is),i=1,atwf_mod%nshell(is))
          IF (cdftlog%thda)THEN
             cdftocl%ocset=.TRUE.
             IF (paral%io_parent)&
                  READ(iunit,END=99,err=99,fmt=*)&
                  (cdftoc%oc2(i,is),i=1,atwf_mod%nshell(is))
          ENDIF
          IF (cntl%cdft)THEN
             CALL filloc(nold,lsold,ocold,atwf_mod%nshell(is),&
                  atwf_mod%lshell(1,is),cdftoc%ocun(1,is))
          ENDIF
       ELSE
          CALL filloc(nold,lsold,ocold,atwf_mod%nshell(is),&
               atwf_mod%lshell(1,is),atwf_mod%oc(1,is))
       ENDIF
       GOTO 10
    ELSE IF (INDEX(line,'SKIP').NE.0) THEN
       ! Skip this atom type
       is=is+1
       IF (is.GT.ions1%nsp) CALL stopgm('READ_BASIS','NSP',& 
            __LINE__,__FILE__)
       GOTO 10
    ELSE
       ! Dummy line
       GOTO 10
    ENDIF
    ! ==--------------------------------------------------------------==
99  CONTINUE
    IF (paral%io_parent)&
         WRITE(6,*) ' READ_BASIS  : ERROR IN READING INPUT FILE'
    CALL stopgm('READ_BASIS',' ',& 
         __LINE__,__FILE__)
100 CONTINUE
    ! ==--------------------------------------------------------------==
    ! Deallocation of local variables
    DEALLOCATE(atrg2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(temp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE read_basis
  ! ==================================================================
  SUBROUTINE filloc(nold,lsold,ocold,nshell,lshell,oc)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nold, lsold(nold)
    REAL(real_8)                             :: ocold(nold)
    INTEGER                                  :: nshell, lshell(nshell)
    REAL(real_8)                             :: oc(nshell)

    INTEGER                                  :: i, is, j, l
    REAL(real_8)                             :: f

    DO is=1,nshell
       oc(is)=0.0_real_8
    ENDDO
    DO i=1,nold
       l=lsold(i)
       f=ocold(i)
       DO j=1,nshell
          IF (l.EQ.lshell(j).AND.oc(j).EQ.0.0_real_8) THEN
             oc(j)=f
             GOTO 10
          ENDIF
       ENDDO
10     CONTINUE
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE filloc
  ! ==================================================================
  SUBROUTINE loadc(c0,foc,ldc,nngw,ndim,ld_foc,is,iat,natst)
    ! ==--------------------------------------------------------------==
    ! == Calculate atomic orbital in plane wave basis.                ==
    ! == 
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ldc, nngw, ndim
    COMPLEX(real_8)                          :: c0(ldc,ndim)
    INTEGER                                  :: ld_foc
    REAL(real_8)                             :: foc(ld_foc)
    INTEGER                                  :: is, iat, natst

    INTEGER, DIMENSION(5), PARAMETER         :: lpp = (/0,1,4,9,16/)

    COMPLEX(real_8)                          :: ci, ei123
    INTEGER                                  :: ig, ish, isub, iv, l, ll, &
                                                lxx, ly
    REAL(real_8)                             :: cc, vol

! Primitive method
! ==--------------------------------------------------------------==

    CALL tiset('     LOADC',isub)
    natsave = 0
    vol=1._real_8/SQRT(parm%omega)
    lxx=0
    DO ish=1,atwf_mod%nshell(is)
       l=atwf_mod%lshell(ish,is)
       ci=(0.0_real_8,1.0_real_8)**l
       ll=2*l+1
       IF (cntl%tipao)THEN
          catord(ish,iat)=lxx+1+natsave
       ENDIF
       DO iv=1,ll
          lxx=lxx+1
          CALL zeroing(c0(:,lxx))!,ldc)
          ly=lpp(l+1)+iv
          IF (atwf_mod%nbcut.EQ.1) THEN
             IF (cntl%bigmem) THEN
                !$omp     parallel do private(IG,CC) shared(LY,NSPLPO)
#ifdef __SR8000
                !poption parallel
#endif
                DO ig=1,nngw
                   cc=curv2(hg(ig),nsplpo,ggnh(1),cat(1,1,ish,is),&
                        cat(1,2,ish,is),0.0_real_8)*vol
                   c0(ig,lxx)=ci*ylmr(ly,ig,gk(1,1))*cc*eigrb(ig,iat)
                ENDDO
             ELSE
                !$omp     parallel do private(IG,CC,EI123) shared(LY,NSPLPO)
#ifdef __SR8000
                !poption parallel
#endif
                DO ig=1,nngw
                   cc=curv2(hg(ig),nsplpo,ggnh(1),cat(1,1,ish,is),&
                        cat(1,2,ish,is),0.0_real_8)*vol
                   ei123=ei1(iat,inyh(1,ig))*ei2(iat,inyh(2,ig))*&
                        ei3(iat,inyh(3,ig))
                   c0(ig,lxx)=ci*ylmr(ly,ig,gk(1,1))*cc*ei123
                ENDDO
             ENDIF
          ELSEIF (atwf_mod%nbcut.EQ.0) THEN
             IF (nngw.GT.ncpw%ngw) CALL stopgm('LOADC','NGW',& 
                  __LINE__,__FILE__)
             !$omp     parallel do private(IG,CC) shared(LY,NSPLPO)
#ifdef __SR8000
             !poption parallel
#endif
             DO ig=1,nngw
                cc=curv2(hg(ig),nsplpo,ggng(1),cat(1,1,ish,is),&
                     cat(1,2,ish,is),0.0_real_8)*vol
                c0(ig,lxx)=ci*ylmr(ly,ig,gk(1,1))*cc*eigr(ig,iat,1)
             ENDDO
          ELSE
             CALL stopgm('LOADC','NBCUT',& 
                  __LINE__,__FILE__)
          ENDIF
          foc(lxx)=atwf_mod%oc(ish,is)/REAL(ll,kind=real_8)
       ENDDO
    ENDDO
    natst=lxx
    natsave=natsave+natst
    CALL tihalt('     LOADC',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE loadc
  ! ==================================================================

END MODULE setbasis_utils
