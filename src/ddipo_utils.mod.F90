MODULE ddipo_utils
  USE atwf,                            ONLY: atwp,&
                                             loadc_foc_array_size
  USE ddip,                            ONLY: ngg1,&
                                             ngg2,&
                                             ngg3,&
                                             ngwmax,&
                                             pdipole,&
                                             pdipolt
  USE dipomod,                         ONLY: moment
  USE dotp_utils,                      ONLY: dotp
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE forcedr_utils,                   ONLY: give_scr_forcedr
  USE gvec,                            ONLY: epsg,&
                                             gvec_com
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE jrotation_utils,                 ONLY: give_scr_jrot,&
                                             jrotation
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE latgen_utils,                    ONLY: latgen
  USE metr,                            ONLY: metr_com
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum
  USE numpw_utils,                     ONLY: gkpwf
  USE opeigr_utils,                    ONLY: give_scr_opeigr,&
                                             opeigr
  USE parac,                           ONLY: parai,&
                                             paral
  USE prcp,                            ONLY: prcp_com
  USE rggen_utils,                     ONLY: recips
  USE rotate_utils,                    ONLY: rotate
  USE sd_wannier_utils,                ONLY: give_scr_sdwann
  USE setbasis_utils,                  ONLY: loadc
  USE sort_utils,                      ONLY: sort2
  USE sphe,                            ONLY: gcutka,&
                                             gcutwmax,&
                                             gcutwmin,&
                                             tsphere
  USE spin,                            ONLY: lspin2,&
                                             spin_mod
  USE system,                          ONLY: cntl,&
                                             mapgp,&
                                             ncpw,&
                                             nkpt,&
                                             parap,&
                                             parm,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: icopy,&
                                             rmatmov
  USE wann,                            ONLY: wan05,&
                                             wannc,&
                                             wanni,&
                                             wannl,&
                                             wannr
  USE wannier_center_utils,            ONLY: wannier_center
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ddipo
  PUBLIC :: give_scr_ddipo
  PUBLIC :: setdip
  PUBLIC :: set_operator

CONTAINS

  ! ==================================================================
  SUBROUTINE ddipo(tau0,c0,cm,c2,sc0,nstate,center)
    ! ==--------------------------------------------------------------==
    ! == CALCULATES THE DIPOLE MOMENT                                 ==
    ! == Modified (CPLS) to compute Wannier functions (August 1997)   ==
    ! ==                                                              ==
    ! == N.B. This version introduces the 'WANNIER_SD' subroutine for ==
    ! ==      the SD procedure and allows to work with the Wannier    ==
    ! ==      functions 'CWANN' using a PARALLEL MACHINE.             ==
    ! == N.B. This version works for NON-CUBIC CELLS too !            ==
    ! ==      Hence x,y,z must be interpreted as 1,2,3                ==
    ! ==      IT GIVES NON-DISTORTED WANNIER CENTERS IN ALL THE       ==
    ! ==      CASES (see 'CNC') :                                     ==
    ! ==                          SIMPLE CUBIC (IBRAV=1)              ==
    ! ==                          TETRAGONAL   (IBRAV=6)              ==
    ! ==                          ORTHOROMBIC  (IBRAV=8,0)            ==
    ! ==                          BCC          (IBRAV=3)              ==
    ! ==                          FCC          (IBRAV=2)              ==
    ! ==                          HEXAGONAL    (IBRAV=4)              ==
    ! ==                          MONOCLINIC   (IBRAV=12)             ==
    ! ==                          TRICLINIC    (IBRAV=14)             ==
    ! ==                                                              ==
    ! == N.B. In this version <i| exp(-Gx*x)*exp(-Gy*y) |i> is        ==
    ! ==      computed BEFORE the unitary rotation (it is useful      ==
    ! ==      for use with GENERAL NON-CUBIC CELLS) and then it is    ==
    ! ==      rotated at the end of the calculation (see CNEW).       ==
    ! ==                                                              ==
    ! == 20/08/1998 Modified by Simone Raugei                         ==
    ! ==            3d ---> 3f,g...                                   ==
    ! ==                                                              ==
    ! == 01/09/1998 Modified by Simone Raugei                         ==
    ! ==            In cntl%md, the WF are updated starting from the last  ==
    ! ==            rotations sequence                                ==
    ! ==                                                              ==
    ! == 03/09/1998 Modified by Simone Raugei                         ==
    ! ==            The subroutines were rearranged                   ==
    ! ==                                                              ==
    ! == 05/09/1998 Modified by Simone Raugei                         ==
    ! ==            A new interface was  added                        ==
    ! ==                                                              ==
    ! == 10/09/1998 Modified by Simone Raugei                         ==
    ! ==            Boero s modifications fo VPP were added           ==
    ! ==                                                              ==
    ! == 14/09/1998 Modified by Simone Raugei                         ==
    ! ==            A different functional involving                  ==
    ! ==            log(|<...>|^2) terms (Resta, Sorella) is used     ==
    ! ==            (see wanlog) (P. L. Silvestrelli)                 ==
    ! ==                                                              ==
    ! == 26/09/1998 Modified by Simone Raugei                         ==
    ! ==            In LSD calculations  the omega maximization is    ==
    ! ==            performed in separately in alpha and beta         ==
    ! ==            subspaces.                                        ==
    ! ==                                                              ==
    ! == 30/09/1998 New modifications for VPP-Fujitsu added           ==
    ! ==            (M. Boero)                                        ==
    ! ==                                                              ==
    ! == 04/12/1998 Rewritten for Version 3.3 (J. Hutter)             ==
    ! ==                                                              ==
    ! == 01/07/1999 Arbitrary symmetry code for Wannier functions     ==
    ! ==            using metric tensor                               ==
    ! ==                                                              ==
    ! == 04/10/1999 Fixed problem in G-vectors for the variable       ==
    ! ==            cell case                                         ==
    ! == 25/11/1999 Fixed dipole moment calculation                   ==
    ! ==--------------------------------------------------------------==
    ! input
    REAL(real_8)                             :: tau0(:,:,:)
    COMPLEX(real_8)                          :: c0(:,:), cm(:,:), c2(:,:), &
                                                sc0(*)
    INTEGER                                  :: nstate
    REAL(real_8), INTENT(out)                :: center(:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'ddipo'

    COMPLEX(real_8)                          :: dd1, dd2, dd3, ddx, det, &
                                                ephase1, ephase2, ephase3
    COMPLEX(real_8), ALLOCATABLE             :: ddmat(:,:)
    COMPLEX(real_8), ALLOCATABLE, &
      DIMENSION(:, :)                        :: lbasis
    COMPLEX(real_8), ALLOCATABLE, SAVE       :: xyzmat(:,:,:)
    INTEGER                                  :: i, ia, iaorb, iat, ierr, &
                                                info, is, isub, ixx, n1, &
                                                natst, nmcol, nxmax, w_opt_tmp
    INTEGER, ALLOCATABLE, SAVE               :: mapcol(:), mapful(:)
    INTEGER, SAVE                            :: i_comp = 0, ifirst = 0
    LOGICAL                                  :: wannier_recomput
    REAL(real_8) :: cost, d1, d2, d3, fac, foc(loadc_foc_array_size), omegon, &
      p(3), phase1, phase2, phase3, sfc, tpi, tpzz, x0(3), zvtot
    REAL(real_8), ALLOCATABLE, DIMENSION(:)  :: ovl, s, work
    REAL(real_8), ALLOCATABLE, &
      DIMENSION(:, :)                        :: u, vt
    REAL(real_8), ALLOCATABLE, SAVE          :: rotlsd(:), rotmat(:,:)

! ==--------------------------------------------------------------==
! ..Wannier functions
! ==--------------------------------------------------------------==
! ==--------------------------------------------------------------==
! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)

    w_opt_tmp=wanni%w_opt
    IF (MOD(i_comp,wan05%loc_relocalize_every)==0 .AND. wan05%loc_relocalize) wanni%w_opt=2

    wannier_recomput=MOD(i_comp,wan05%loc_recompute_dipole_matrices_every)==0
    i_comp=i_comp+1
    IF (paral%io_parent.AND.wannier_recomput) &
         WRITE(6,*) procedureN//': RECOMPUTE DIPOLE MATRICES'

    ! ..initialization on first call to subroutine
    IF (ifirst.EQ.0) THEN
       n1=0
       DO i=0,parai%nproc-1
          n1=MAX(n1,parap%sparm(3,i))
       ENDDO
       ngwmax=n1
       ALLOCATE(mapful(2*(spar%ngws)),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       nmcol=(parai%nproc*ngwmax)
       ALLOCATE(mapcol(nmcol),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL setdip(mapful,mapcol)
       ! ..initialization for Wannier function part
       IF (wannl%twann) THEN
          CALL set_operator(.TRUE.)
          ALLOCATE(xyzmat(nstate,nstate,wannc%nwanopt),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          CALL zeroing(xyzmat)!,nstate*nstate*wannc%nwanopt)
          ALLOCATE(rotmat(nstate,nstate),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          CALL zeroing(rotmat)
          IF (cntl%tlsd) THEN
             nxmax=MAX(spin_mod%nsup,spin_mod%nsdown)
             ALLOCATE(rotlsd(nxmax*nxmax),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             CALL zeroing(rotlsd)
          ELSEIF (lspin2%tlse) THEN
             nxmax=nstate-2
             ALLOCATE(rotlsd(nxmax*nxmax),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             CALL zeroing(rotlsd)
          ENDIF
       ENDIF
       ifirst=1
    ENDIF
    IF (cntl%tlsd) THEN
       fac=1._real_8
    ELSEIF (lspin2%tlse) THEN
       fac=1._real_8
    ELSE
       fac=2._real_8
    ENDIF
    omegon=1.0_real_8/parm%omega
    ! ..ionic contribution
    ! ..center of charge
    p(1)=0._real_8
    p(2)=0._real_8
    p(3)=0._real_8
    x0(1)=0._real_8
    x0(2)=0._real_8
    x0(3)=0._real_8
    zvtot=0._real_8
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          x0(1)=x0(1)+ions0%zv(is)*tau0(1,ia,is)
          x0(2)=x0(2)+ions0%zv(is)*tau0(2,ia,is)
          x0(3)=x0(3)+ions0%zv(is)*tau0(3,ia,is)
       ENDDO
       zvtot=zvtot+ions0%na(is)*ions0%zv(is)
    ENDDO
    x0(1)=x0(1)/zvtot
    x0(2)=x0(2)/zvtot
    x0(3)=x0(3)/zvtot
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          p(1)=p(1)+ions0%zv(is)*(tau0(1,ia,is)-x0(1))
          p(2)=p(2)+ions0%zv(is)*(tau0(2,ia,is)-x0(2))
          p(3)=p(3)+ions0%zv(is)*(tau0(3,ia,is)-x0(3))
       ENDDO
    ENDDO
    p(1)=p(1)*omegon
    p(2)=p(2)*omegon
    p(3)=p(3)*omegon
    IF (ABS(p(1)+p(2)+p(3)) .GT. 1.e-10_real_8)&
         CALL stopgm(procedureN//"|","REFERENCE POINT",& 
         __LINE__,__FILE__)
    ! ..phase factors due to centor of charge motion
    tpzz=parm%tpiba*nstate
    phase1=(gvec_com%b1(1)*x0(1)+gvec_com%b1(2)*x0(2)+gvec_com%b1(3)*x0(3))*tpzz
    phase2=(gvec_com%b2(1)*x0(1)+gvec_com%b2(2)*x0(2)+gvec_com%b2(3)*x0(3))*tpzz
    phase3=(gvec_com%b3(1)*x0(1)+gvec_com%b3(2)*x0(2)+gvec_com%b3(3)*x0(3))*tpzz
    ephase1=CMPLX(COS(phase1),SIN(phase1),kind=real_8)
    ephase2=CMPLX(COS(phase2),SIN(phase2),kind=real_8)
    ephase3=CMPLX(COS(phase3),SIN(phase3),kind=real_8)
    ! ..electronic contribution
    ALLOCATE(ddmat(nstate, nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    IF (wannier_recomput) THEN
       CALL opeigr(c0,c2,sc0,nstate,mapful,mapcol,ddmat,&
            1,0,dd1)
       IF (wannl%twann) CALL dcopy(2*nstate*nstate,ddmat,1,xyzmat(1,1,1),1)
       det=dd1*ephase1
       d1=ATAN(AIMAG(det)/REAL(det))
       CALL opeigr(c0,c2,sc0,nstate,mapful,mapcol,ddmat,&
            2,0,dd2)
       IF (wannl%twann) CALL dcopy(2*nstate*nstate,ddmat,1,xyzmat(1,1,2),1)
       det=dd2*ephase2
       d2=ATAN(AIMAG(det)/REAL(det))
       CALL opeigr(c0,c2,sc0,nstate,mapful,mapcol,ddmat,&
            3,0,dd3)
       IF (wannl%twann) CALL dcopy(2*nstate*nstate,ddmat,1,xyzmat(1,1,3),1)
       det=dd3*ephase3
       d3=ATAN(AIMAG(det)/REAL(det))
       ! ..Sum up dipoles
       tpi=2._real_8*ACOS(-1._real_8)
       cost=fac*omegon/tpi
       DO i=1,3
          pdipole(i)=cost*(d1*parm%a1(i)+d2*parm%a2(i)+d3*parm%a3(i))
          pdipolt(i)=pdipole(i) + p(i)
          moment%dmom(i)=parm%omega*pdipolt(i)! cmb
       ENDDO
    ENDIF
    ! Wannier function calculation
    IF (wannl%twann) THEN
       ! Additional operators
       IF (wannier_recomput) THEN
          DO i=4,wannc%nwanopt
             CALL opeigr(c0,c2,sc0,nstate,mapful,mapcol,ddmat,&
                  wannc%iow(1,i-3),wannc%iow(2,i-3),ddx)
             CALL dcopy(2*nstate*nstate,ddmat,1,xyzmat(1,1,i),1)
          ENDDO
       ENDIF
       ! Optimisers
       IF (wanni%w_opt.EQ.1) THEN
          IF (paral%parent) THEN
             ! steepest descent
             IF (cntl%tlsd) THEN
                CALL zeroing(rotmat)!,nstate*nstate)
                nxmax=MAX(spin_mod%nsup,spin_mod%nsdown)
                CALL sd_wannier(rotlsd,xyzmat,nstate,spin_mod%nsup)
                CALL rmatmov(spin_mod%nsup,spin_mod%nsup,rotlsd,spin_mod%nsup,rotmat,nstate)
                CALL sd_wannier(rotlsd,xyzmat(spin_mod%nsup+1,spin_mod%nsup+1,1),&
                     nstate,spin_mod%nsdown)
                CALL rmatmov(spin_mod%nsdown,spin_mod%nsdown,rotlsd,spin_mod%nsdown,&
                     rotmat(spin_mod%nsup+1,spin_mod%nsup+1),nstate)
             ELSEIF (lspin2%tlse) THEN
                nxmax=nstate-2
                CALL sd_wannier(rotlsd,xyzmat,nstate,nxmax)
                CALL zeroing(rotmat)!,nstate*nstate)
                CALL rmatmov(nxmax,nxmax,rotlsd,nxmax,rotmat,nstate)
                rotmat(nxmax+1,nxmax+1)=1._real_8
                rotmat(nstate,nstate)=1._real_8
             ELSE
                CALL sd_wannier(rotmat,xyzmat,nstate,nstate)
             ENDIF
          ENDIF
          CALL mp_bcast(rotmat,SIZE(rotmat),parai%source,parai%allgrp)
       ELSEIF (wanni%w_opt.EQ.2) THEN
          ! 
          ! Jacobi rotations
          IF (cntl%tlsd) THEN
             CALL zeroing(rotmat)!,nstate*nstate)
             nxmax=MAX(spin_mod%nsup,spin_mod%nsdown)
             CALL jrotation(rotlsd,xyzmat(1,1,1),nstate,spin_mod%nsup)
             CALL rmatmov(spin_mod%nsup,spin_mod%nsup,rotlsd,spin_mod%nsup,rotmat,nstate)
             CALL jrotation(rotlsd,xyzmat(spin_mod%nsup+1,spin_mod%nsup+1,1),&
                  nstate,spin_mod%nsdown)
             CALL rmatmov(spin_mod%nsdown,spin_mod%nsdown,rotlsd,spin_mod%nsdown,&
                  rotmat(spin_mod%nsup+1,spin_mod%nsup+1),nstate)
          ELSEIF (lspin2%tlse) THEN
             nxmax=nstate-2
             CALL jrotation(rotlsd,xyzmat,nstate,nxmax)
             CALL zeroing(rotmat)!,nstate*nstate)
             CALL rmatmov(nxmax,nxmax,rotlsd,nxmax,rotmat,nstate)
             rotmat(nxmax+1,nxmax+1)=1._real_8
             rotmat(nstate,nstate)=1._real_8
          ELSE
             CALL jrotation(rotmat,xyzmat,nstate,nstate)
          ENDIF
       ELSEIF (wanni%w_opt.EQ.3) THEN
          ! SVD


          IF (cntl%tlsd) CALL stopgm(procedureN,'NYI',& 
               __LINE__,__FILE__)


          ! mem alloc
          ALLOCATE(ovl(atwp%nattot*nstate), lbasis(nkpt%ngwk,atwp%nattot),&
               s(nstate), u(nstate,nstate), vt(nstate,nstate),&
               work(50*nstate), stat=ierr)
          IF (ierr/=0) CALL stopgm(procedureN,'Allocation problem',& 
               __LINE__,__FILE__)


          ! LOAD ATOMIC GUESS TO PW BASIS
          iaorb=1
          iat=0
          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                iat=iat+1
                CALL loadc(lbasis(1,iaorb),foc,ncpw%ngw,ncpw%ngw,atwp%nattot-iaorb+1,&
                     SIZE(foc),is,iat,natst)
                DO ixx=iaorb,iaorb+natst-1
                   sfc=dotp(ncpw%ngw,lbasis(:,ixx),lbasis(:,ixx))
                   CALL mp_sum(sfc,parai%allgrp)
                   IF (sfc.EQ.0._real_8) THEN
                      CALL stopgm(procedureN,'WRONG ATOMIC ORBITAL',& 
                           __LINE__,__FILE__)
                   ELSE
                      sfc=1._real_8/SQRT(sfc)
                   ENDIF
                   CALL dscal(2*ncpw%ngw,sfc,lbasis(1,ixx),1)
                ENDDO
                iaorb=iaorb+natst
             ENDDO
          ENDDO


          ! COMPUTE OVERLAP WITH *LOCALIZED* WAVEFUNCTIONS: OVL = CA' * C0
          ovl(:)=0.0_real_8
          CALL ovlap2(ncpw%ngw,atwp%nattot,nstate,ovl,lbasis,c0,.TRUE.)
          CALL mp_sum(ovl,atwp%nattot*nstate,parai%allgrp)


          ! PROJECT: CA * (CA'*C0)
          CALL dgemm('N','N',2*ncpw%ngw,nstate,atwp%nattot,1.0_real_8,lbasis,2*ncpw%ngw,&
               ovl,atwp%nattot,0.0_real_8,c2,2*ncpw%ngw)


          ! COMPUTE OVERLAP WITH LOW RANK LOCALIZED WAVEFUNCTIONS: OVL = C0 * ( CA * (CA'*C0) )
          ovl(:)=0.0_real_8
          CALL ovlap2(ncpw%ngw,nstate,nstate,ovl,c0,c2,.TRUE.)
          CALL mp_sum(ovl,nstate**2,parai%allgrp)

          ! COMPUTE SVD OF OVERLAP MATRIX: SVD( OVL )
          IF (paral%io_parent) WRITE(6,*) "COMPUTE SVD"
          CALL dgesvd('A','A', nstate,nstate,ovl,nstate,s,u,&
               nstate,vt,nstate,work,SIZE(work),info)


          ! COMPUTE ROTATION MATRIX
          IF (paral%io_parent) WRITE(6,*) "COMPUTE ROTATION MATRIX"
          CALL dgemm("N","N",nstate,nstate,nstate,1._real_8,u,nstate,vt,&
               nstate,0._real_8,rotmat,nstate)


          ! BCAST TO ALL PROCS
          CALL mp_bcast(rotmat,SIZE(rotmat),parai%io_source,parai%cp_grp)


          ! mem dealloc
          DEALLOCATE(ovl,lbasis,s,u,vt,work, stat=ierr)
          IF (ierr/=0) CALL stopgm(procedureN,'Deallocation problem',& 
               __LINE__,__FILE__)


       ELSE
          CALL stopgm(procedureN,'not implemented',& 
               __LINE__,__FILE__)
       ENDIF
       ! reset randomization variable (only randomize once)
       wannr%w_ran=0._real_8
       ! rotate the wavefunctions and velocities 
       ! to the Wannier representation
       CALL rotate(1._real_8,c0,0._real_8,c2,rotmat,nstate,&
            2*ncpw%ngw,cntl%tlsd,spin_mod%nsup,spin_mod%nsdown)
       CALL dcopy(2*ncpw%ngw*nstate,c2,1,c0,1)
       CALL rotate(1._real_8,cm,0._real_8,c2,rotmat,nstate,&
            2*ncpw%ngw,cntl%tlsd,spin_mod%nsup,spin_mod%nsdown)
       CALL dcopy(2*ncpw%ngw*nstate,c2,1,cm,1)
       ! calculate the centers and spread of the wannier functions
       CALL wannier_center(xyzmat,nstate,nstate,center,tau0)
    ENDIF
    DEALLOCATE(ddmat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE ddipo
  ! ==================================================================
  SUBROUTINE give_scr_ddipo(lddipo,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lddipo
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: lforcedr, lopeigr, lself, &
                                                lsetdip, lwanc, lwann, nmax, &
                                                nstate

    nstate=crge%n
    nmax=nstate
    IF (cntl%tlsd) nmax=MAX(spin_mod%nsup,spin_mod%nsdown)
    lself=2*nstate*nstate
    lwanc=4*nstate
    CALL give_scr_opeigr(lopeigr,tag,nstate)
    lwann=0
    lforcedr=0
    IF (wannl%twann)THEN
       IF (wannc%nwanopt.EQ.0) CALL set_operator(.FALSE.)
       IF (wanni%w_opt.EQ.1) THEN
          CALL give_scr_sdwann(lwann,tag,nmax)
       ELSEIF (wanni%w_opt.EQ.2) THEN
          CALL give_scr_jrot(lwann,tag,nmax,.FALSE.)
       ENDIF
       IF (wannl%twmol.OR.wannl%twdos) THEN
          CALL give_scr_forcedr(lforcedr,tag,crge%n,.FALSE.,.TRUE.)
       ENDIF
       lwann=MAX(lwann,lforcedr)
    ENDIF
    lsetdip=spar%ngws/2+1
    lddipo=MAX(lsetdip,lself+lopeigr,lself+lwann,lself+lwanc)+100
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_ddipo
  ! ==================================================================
  SUBROUTINE setdip(mapful,mapcol)
    ! NOTE: MAPCOL is unused 
    ! ==--------------------------------------------------------------==
    ! input
    INTEGER                                  :: mapful(2,*), mapcol(ngwmax,*)

    INTEGER                                  :: blklen, i, i1, i2, i3, icurr, &
                                                ierr, ig, isub, it, j, j1, &
                                                j2, j3, k, kmax, kmin, ngg1h, &
                                                ngg2h, ngg3h
    INTEGER, ALLOCATABLE                     :: iny(:,:)
    REAL(real_8)                             :: aux1(3), aux2(3), aux3(3), &
                                                bax1(3), bax2(3), bax3(3), &
                                                g2, omegan, t1, t2, t3
    REAL(real_8), ALLOCATABLE                :: hgt(:)

#ifdef __VECTOR 
    INTEGER :: ijk,iesri,iesrj,iesrk,IALL
#endif
    INTEGER, ALLOCATABLE :: iscr(:)
    CHARACTER(*),PARAMETER::procedureN='SETDIP'
    ! ==--------------------------------------------------------------==
    CALL tiset('    SETDIP',isub)
    ALLOCATE(iny(3,spar%ngws),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(hgt(spar%ngws),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(iscr(2*spar%ngws),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL icopy(ncpw%ngw,mapgp,1,iscr,1)
    CALL zeroing(iscr(ncpw%ngw+1:ncpw%ngw+ngwmax))!,ngwmax-ngw)
    blklen=ngwmax*(8/2)
    CALL my_concat(iscr,mapcol,blklen,parai%allgrp)
    ! ..construction of g-vectors
    IF (cntl%tprcp.OR.cntl%tpres) THEN   ! variable cell case or stress tensor
       CALL latgen(parm%ibrav,prcp_com%cellrf,aux1,aux2,aux3,omegan)
       CALL recips(parm%alat,aux1,aux2,aux3,bax1,bax2,bax3)
    ELSE                      ! fixed cell case
       bax1(1)=gvec_com%b1(1)
       bax1(2)=gvec_com%b1(2)
       bax1(3)=gvec_com%b1(3)
       bax2(1)=gvec_com%b2(1)
       bax2(2)=gvec_com%b2(2)
       bax2(3)=gvec_com%b2(3)
       bax3(1)=gvec_com%b3(1)
       bax3(2)=gvec_com%b3(2)
       bax3(3)=gvec_com%b3(3)
    ENDIF
#ifdef __VECTOR
    iesri=(2*spar%nr3s-1)
    iesrj=(2*spar%nr2s-1)*iesri
    iesrk=(spar%nr1s-1)*iesrj
    IALL=iesrk
#endif
    ! Maximum K point amplitude (Kmax) in Brillouin zone.
    gcutka=gkpwf(bax1,bax2,bax3,tkpts%tkpnt.AND.tsphere)
    IF (gcutka.EQ.0._real_8) THEN
       gcutwmin=gvec_com%gcutw
       gcutwmax=gvec_com%gcutw
    ELSE
       ! For each k and each G, |k+G|^2 < GCUTW
       ! so we define GCUTWMIN and GCUTWMAX as:
       ! GCUTWMIN < |G|^2 < GCUTWMAX we have to apply a mask.
       gcutwmin = MAX(0._real_8,SQRT(gvec_com%gcutw)-SQRT(gcutka))**2
       gcutwmax = (SQRT(gvec_com%gcutw)+SQRT(gcutka))**2
    ENDIF
    ig=0
    ! mb...I=0 
    DO j=0,spar%nr2s-1
       kmin=-spar%nr3s+1
       kmax=spar%nr3s-1
       IF (j.EQ.0) kmin=0
       DO k=kmin,kmax
          t1=REAL(j,kind=real_8)*bax2(1)+REAL(k,kind=real_8)*bax3(1)
          t2=REAL(j,kind=real_8)*bax2(2)+REAL(k,kind=real_8)*bax3(2)
          t3=REAL(j,kind=real_8)*bax2(3)+REAL(k,kind=real_8)*bax3(3)
          g2=t1*t1+t2*t2+t3*t3
          IF (g2.LT.gcutwmax) THEN
             ig=ig+1
             hgt(ig)=g2-SQRT(REAL(ig-1,kind=real_8))*epsg
             iny(1,ig)=0
             iny(2,ig)=j
             iny(3,ig)=k
          ENDIF
       ENDDO
    ENDDO
    ! mb...I<>0
#ifdef __VECTOR 
    DO ijk=1,IALL
       i=1+INT((ijk-1)/iesrj)
       j=1-spar%nr2s+INT(MOD(ijk-1,iesrj)/iesri)
       k=1-spar%nr3s+MOD(ijk-1,iesri)
       t1=REAL(i,kind=real_8)*bax1(1)+REAL(j,kind=real_8)*bax2(1)+REAL(k,kind=real_8)*bax3(1)
       t2=REAL(i,kind=real_8)*bax1(2)+REAL(j,kind=real_8)*bax2(2)+REAL(k,kind=real_8)*bax3(2)
       t3=REAL(i,kind=real_8)*bax1(3)+REAL(j,kind=real_8)*bax2(3)+REAL(k,kind=real_8)*bax3(3)
       g2=t1*t1+t2*t2+t3*t3
       IF (g2.LT.gcutwmax) THEN
          ig=ig+1
          hgt(ig)=g2-SQRT(REAL(ig-1,kind=real_8))*epsg
          iny(1,ig)=i
          iny(2,ig)=j
          iny(3,ig)=k
       ENDIF
    ENDDO
#else 
    DO i=1,spar%nr1s-1
       DO j=-spar%nr2s+1,spar%nr2s-1
          DO k=-spar%nr3s+1,spar%nr3s-1
             t1=REAL(i,kind=real_8)*bax1(1)+REAL(j,kind=real_8)*bax2(1)+REAL(k,kind=real_8)*bax3(1)
             t2=REAL(i,kind=real_8)*bax1(2)+REAL(j,kind=real_8)*bax2(2)+REAL(k,kind=real_8)*bax3(2)
             t3=REAL(i,kind=real_8)*bax1(3)+REAL(j,kind=real_8)*bax2(3)+REAL(k,kind=real_8)*bax3(3)
             g2=t1*t1+t2*t2+t3*t3
             IF (g2.LT.gcutwmax) THEN
                ig=ig+1
                hgt(ig)=g2-SQRT(REAL(ig-1,kind=real_8))*epsg
                iny(1,ig)=i
                iny(2,ig)=j
                iny(3,ig)=k
             ENDIF
          ENDDO
       ENDDO
    ENDDO
#endif
    IF (ig.NE.spar%ngws) CALL stopgm('SETDIP',&
         'INCONSISTENT NUMBER OF PLANE WAVES',& 
         __LINE__,__FILE__)
    CALL sort2(hgt,spar%ngws,iscr)
    DO ig=1,spar%ngws-1
       icurr=ig
30     IF (iscr(icurr).NE.ig) THEN
          it=iny(1,icurr)
          iny(1,icurr)=iny(1,iscr(icurr))
          iny(1,iscr(icurr))=it
          it=iny(2,icurr)
          iny(2,icurr)=iny(2,iscr(icurr))
          iny(2,iscr(icurr))=it
          it=iny(3,icurr)
          iny(3,icurr)=iny(3,iscr(icurr))
          iny(3,iscr(icurr))=it
          it=icurr
          icurr=iscr(icurr)
          iscr(it)=it
          IF (iscr(icurr).EQ.ig) THEN
             iscr(icurr)=icurr
             GOTO 25
          ENDIF
          GOTO 30
       ENDIF
25     CONTINUE
    ENDDO
    ngg1=0
    ngg2=0
    ngg3=0
#ifdef __SR8000
    !poption parallel, tlocal(IG), max(NGG1,NGG2,NGG3)
#endif
    !$omp parallel do private(IG) &
    !$omp             reduction(max:NGG1,NGG2,NGG3)
    DO ig=1,spar%ngws
       ngg1=MAX(ABS(iny(1,ig)),ngg1)
       ngg2=MAX(ABS(iny(2,ig)),ngg2)
       ngg3=MAX(ABS(iny(3,ig)),ngg3)
    ENDDO
    ngg1h=ngg1+2
    ngg2h=ngg2+2
    ngg3h=ngg3+2
    ngg1=2*(ngg1+1)+1
    ngg2=2*(ngg2+1)+1
    ngg3=2*(ngg3+1)+1
#ifdef __SR8000
    !poption parallel, tlocal(IG,I1,I2,I3,J1,J2,J3)
#endif
    !$omp parallel do private(IG,I1,I2,I3,J1,J2,J3)
    DO ig=1,spar%ngws
       i1=iny(1,ig)+ngg1h
       i2=iny(2,ig)+ngg2h
       i3=iny(3,ig)+ngg3h
       j1=-iny(1,ig)+ngg1h
       j2=-iny(2,ig)+ngg2h
       j3=-iny(3,ig)+ngg3h
       mapful(1,ig)=i1+i2*ngg1+i3*(ngg1*ngg2)
       mapful(2,ig)=j1+j2*ngg1+j3*(ngg1*ngg2)
    ENDDO
    DEALLOCATE(hgt,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(iny,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(iscr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    CALL tihalt('    SETDIP',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE setdip
  ! ==================================================================
  SUBROUTINE set_operator(tpri)
    ! ==--------------------------------------------------------------==
    LOGICAL                                  :: tpri

    CHARACTER(len=2), DIMENSION(1:3)         :: xyz = (/' X',' Y',' Z'/)
    INTEGER                                  :: i, k, l
    REAL(real_8)                             :: metric(3,3)

    CALL dgemm('N','T',3,3,3,1._real_8,metr_com%ht,3,metr_com%ht,3,0._real_8,metric,3)
    CALL zeroing(wannc%iow)!,6)
    CALL zeroing(wannc%wwei)!,6)
    wannc%nwanopt=6

    wannc%wwei(1)=metric(1,1)-metric(1,2)-metric(1,3)
    wannc%wwei(2)=metric(2,2)-metric(1,2)-metric(2,3)
    wannc%wwei(3)=metric(3,3)-metric(1,3)-metric(2,3)
    wannc%wwei(4)=metric(1,2)
    wannc%wwei(5)=metric(1,3)
    wannc%wwei(6)=metric(2,3)
    wannc%iow(1,1)=1
    wannc%iow(2,1)=2
    wannc%iow(1,2)=1
    wannc%iow(2,2)=3
    wannc%iow(1,3)=2
    wannc%iow(2,3)=3

    ! RV----quadrupole------------
    ! xx
    DO k=1,3
       DO l=1,3
          metric(k,l)=metr_com%ht(k,1)*metr_com%ht(l,1)
       ENDDO
    ENDDO
    wannc%wquadi(1,1)=metric(1,1)-metric(1,2)-metric(1,3)
    wannc%wquadi(1,2)=metric(2,2)-metric(1,2)-metric(2,3)
    wannc%wquadi(1,3)=metric(3,3)-metric(1,3)-metric(2,3)
    wannc%wquadi(1,4)=metric(1,2)
    wannc%wquadi(1,5)=metric(1,3)
    wannc%wquadi(1,6)=metric(2,3)

    ! xy 

    DO k=1,3
       DO l=1,3
          metric(k,l)=0.5_real_8*(metr_com%ht(k,1)*metr_com%ht(l,2)+metr_com%ht(k,2)*metr_com%ht(l,1))
       ENDDO
    ENDDO
    wannc%wquadi(2,1)=metric(1,1)-metric(1,2)-metric(1,3)
    wannc%wquadi(2,2)=metric(2,2)-metric(1,2)-metric(2,3)
    wannc%wquadi(2,3)=metric(3,3)-metric(1,3)-metric(2,3)
    wannc%wquadi(2,4)=metric(1,2)
    wannc%wquadi(2,5)=metric(1,3)
    wannc%wquadi(2,6)=metric(2,3)

    ! xz

    DO k=1,3
       DO l=1,3
          metric(k,l)=0.5_real_8*(metr_com%ht(k,1)*metr_com%ht(l,3)+metr_com%ht(k,3)*metr_com%ht(l,1))
       ENDDO
    ENDDO
    wannc%wquadi(3,1)=metric(1,1)-metric(1,2)-metric(1,3)
    wannc%wquadi(3,2)=metric(2,2)-metric(1,2)-metric(2,3)
    wannc%wquadi(3,3)=metric(3,3)-metric(1,3)-metric(2,3)
    wannc%wquadi(3,4)=metric(1,2)
    wannc%wquadi(3,5)=metric(1,3)
    wannc%wquadi(3,6)=metric(2,3)

    ! yy
    DO k=1,3
       DO l=1,3
          metric(k,l)=metr_com%ht(k,2)*metr_com%ht(l,2)
       ENDDO
    ENDDO
    wannc%wquadi(4,1)=metric(1,1)-metric(1,2)-metric(1,3)
    wannc%wquadi(4,2)=metric(2,2)-metric(1,2)-metric(2,3)
    wannc%wquadi(4,3)=metric(3,3)-metric(1,3)-metric(2,3)
    wannc%wquadi(4,4)=metric(1,2)
    wannc%wquadi(4,5)=metric(1,3)
    wannc%wquadi(4,6)=metric(2,3)

    ! yz

    DO k=1,3
       DO l=1,3
          metric(k,l)=0.5_real_8*(metr_com%ht(k,2)*metr_com%ht(l,3)+metr_com%ht(k,3)*metr_com%ht(l,2))
       ENDDO
    ENDDO
    wannc%wquadi(5,1)=metric(1,1)-metric(1,2)-metric(1,3)
    wannc%wquadi(5,2)=metric(2,2)-metric(1,2)-metric(2,3)
    wannc%wquadi(5,3)=metric(3,3)-metric(1,3)-metric(2,3)
    wannc%wquadi(5,4)=metric(1,2)
    wannc%wquadi(5,5)=metric(1,3)
    wannc%wquadi(5,6)=metric(2,3)

    ! zz
    DO k=1,3
       DO l=1,3
          metric(k,l)=metr_com%ht(k,3)*metr_com%ht(l,3)
       ENDDO
    ENDDO
    wannc%wquadi(6,1)=metric(1,1)-metric(1,2)-metric(1,3)
    wannc%wquadi(6,2)=metric(2,2)-metric(1,2)-metric(2,3)
    wannc%wquadi(6,3)=metric(3,3)-metric(1,3)-metric(2,3)
    wannc%wquadi(6,4)=metric(1,2)
    wannc%wquadi(6,5)=metric(1,3)
    wannc%wquadi(6,6)=metric(2,3)
    ! RV----quadrupole------------

    IF (paral%parent.AND.tpri) THEN
       DO i=1,3
          IF (paral%io_parent)&
               WRITE(6,'(T9,A,T19,I3,T30,A,T40,A,T56,F10.4)')&
               ' OPERATOR:',i,xyz(i),' WEIGHT/L^2=',wannc%wwei(i)/parm%alat**2
       ENDDO
       DO i=4,6
          IF (paral%io_parent)&
               WRITE(6,'(T9,A,T19,I3,T25,A,T30,A,T40,A,T56,F10.4)')&
               ' OPERATOR:',i,xyz(wannc%iow(1,i-3)),xyz(wannc%iow(2,i-3)),&
               ' WEIGHT/L^2=',wannc%wwei(i)/parm%alat**2
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE set_operator
  ! ==================================================================

END MODULE ddipo_utils
