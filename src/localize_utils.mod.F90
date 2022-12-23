MODULE localize_utils
  USE atwf,                            ONLY: atwp,&
                                             loadc_foc_array_size
  USE ddip,                            ONLY: ngwmax
  USE ddipo_utils,                     ONLY: set_operator,&
                                             setdip
  USE dotp_utils,                      ONLY: dotp
  USE error_handling,                  ONLY: stopgm
  USE forcedr_driver,                  ONLY: forcedr
  USE forcep_utils,                    ONLY: rhoe_psi_size
  USE gvec,                            ONLY: gvec_com
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE jrotation_utils,                 ONLY: jrotation
  USE kinds,                           ONLY: real_8
  USE linres,                          ONLY: td03,&
                                             wcent
  USE molorb_utils,                    ONLY: molorb
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum
  USE nmr_position_p_utils,            ONLY: bound
  USE opeigr_utils,                    ONLY: opeigr
  USE parac,                           ONLY: parai,&
                                             paral
  USE response_pmod,                   ONLY: lower_left,&
                                             lower_left_value,&
                                             response1,&
                                             voa_data,&
                                             wanniercenters
  USE rmas,                            ONLY: rmass
  USE rotate_utils,                    ONLY: rotate
  USE setbasis_utils,                  ONLY: loadc
  USE spin,                            ONLY: spin_mod
  USE system,                          ONLY: cntl,&
                                             maxsys,&
                                             ncpw,&
                                             nkpt,&
                                             parap,&
                                             parm,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: rmatmov
  USE wann,                            ONLY: wan05,&
                                             wannc,&
                                             wanni,&
                                             wannl,&
                                             wannr
  USE wannier_center_utils,            ONLY: wannier_center
  USE wannier_print_utils,             ONLY: wannier_print
  USE wc_dos_utils,                    ONLY: wc_dos
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: localize
  PUBLIC :: localize2
  PUBLIC :: wc_print
  !public :: wannier_save

CONTAINS

  ! ==================================================================
  SUBROUTINE localize(tau0,c0,c2,sc0,nstate)
    ! ==--------------------------------------------------------------==
    ! == CALCULATES LOCALIZED FUNCTIONS                               ==
    ! ==--------------------------------------------------------------==
    ! input
    REAL(real_8)                             :: tau0(:,:,:)
    COMPLEX(real_8)                          :: sc0(*)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c2(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'localize'

    COMPLEX(real_8)                          :: dd1, dd2, dd3
    COMPLEX(real_8), ALLOCATABLE             :: ddmat(:,:), psi(:,:)
    COMPLEX(real_8), ALLOCATABLE, &
      DIMENSION(:, :)                        :: lbasis
    COMPLEX(real_8), ALLOCATABLE, SAVE       :: xyzmat(:,:,:)
    INTEGER :: i, ia, iaorb, iat, ierr, il_psi_1d, il_psi_2d, il_rhoe_1d, &
      il_rhoe_2d, info, is, isub, ix1, ix2, ixl, ixx, n1, natst, nmcol, nxmax
    INTEGER, ALLOCATABLE, SAVE               :: mapcol(:), mapful(:)
    REAL(real_8)                             :: foc(loadc_foc_array_size), sfc
    REAL(real_8), ALLOCATABLE                :: center(:,:), eigv(:), &
                                                fion(:,:,:), rhoe(:,:)
    REAL(real_8), ALLOCATABLE, DIMENSION(:)  :: ovl, s, work
    REAL(real_8), ALLOCATABLE, &
      DIMENSION(:, :)                        :: u, vt
    REAL(real_8), ALLOCATABLE, SAVE          :: rotlsd(:), rotmat(:,:)

! 
! ==--------------------------------------------------------------==
! ..Wannier functions
! ==--------------------------------------------------------------==
! ==--------------------------------------------------------------==
! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    IF (paral%io_parent) THEN
       WRITE(6,'(/,1X,64("*"))')
       WRITE(6,'(1X,A)') 'LOCALIZATION OF WAVEFUNCTION'
       IF (wanni%w_type.EQ.1)&
            WRITE(6,'(A)')   '         FUNCTIONAL: VANDERBILT |M|^2'
       IF (wanni%w_type.EQ.2)&
            WRITE(6,'(A)')   '         FUNCTIONAL: RESTA Log|M|^2'
       IF (wanni%w_opt.EQ.1)&
            WRITE(6,'(A)') '         OPTIMISATION: STEEPEST DESCENT'
       IF (wanni%w_opt.EQ.2)&
            WRITE(6,'(A)') '         OPTIMISATION: JACOBI ROTATION'
       WRITE(6,'(A,T56,1PE10.4)')&
            '         CONVERGENCE CRITERIA:',wannr%w_eps
       WRITE(6,'(A,T56,I10)')&
            '         MAXIMUM # OF STEPS:',wanni%w_maxs
       WRITE(6,'(A,T56,1PE10.4)')&
            '         RANDOMIZATION AMPLITUDE:',wannr%w_ran
       WRITE(6,'(A,T56,1PE10.4)')&
            '         STEP SIZE:',wannr%w_step
    ENDIF
    ! ..initialization
    n1=0
    DO i=0,parai%nproc-1
       n1=MAX(n1,parap%sparm(3,i))
    ENDDO
    ngwmax=n1
    ALLOCATE(mapful(2*spar%ngws),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    nmcol=parai%nproc*ngwmax
    ALLOCATE(mapcol(nmcol),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL setdip(mapful,mapcol)
    ! ..initialization for Wannier function part
    CALL set_operator(.TRUE.)
    ALLOCATE(xyzmat(nstate,nstate,wannc%nwanopt),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(xyzmat)!,nstate*nstate*nwanopt)
    ALLOCATE(rotmat(nstate,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    IF (cntl%tlsd) THEN
       nxmax=MAX(spin_mod%nsup,spin_mod%nsdown)
       ALLOCATE(rotlsd(nxmax*nxmax),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! ..electronic contribution
    ALLOCATE(ddmat(nstate, nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL opeigr(c0,c2,sc0,nstate,mapful,mapcol,ddmat,&
         1,0,dd1)
    CALL dcopy(2*nstate*nstate,ddmat,1,xyzmat(1,1,1),1)
    CALL opeigr(c0,c2,sc0,nstate,mapful,mapcol,ddmat,&
         2,0,dd2)
    CALL dcopy(2*nstate*nstate,ddmat,1,xyzmat(1,1,2),1)
    CALL opeigr(c0,c2,sc0,nstate,mapful,mapcol,ddmat,&
         3,0,dd3)
    CALL dcopy(2*nstate*nstate,ddmat,1,xyzmat(1,1,3),1)
    DO i=4,wannc%nwanopt
       CALL opeigr(c0,c2,sc0,nstate,mapful,mapcol,ddmat,&
            wannc%iow(1,i-3),wannc%iow(2,i-3),dd1)
       CALL dcopy(2*nstate*nstate,ddmat,1,xyzmat(1,1,i),1)
    ENDDO
    DEALLOCATE(ddmat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)

    ! Calculation of the spread BEFORE the actual localization
    IF (paral%io_parent) THEN
       WRITE(6,'(/,1X,64("*"))')
       WRITE(6,'(" ****",6X,A,6X,"****",/)')&
            'CENTERS AND SPREAD BEFORE THE OPTIMIZATION'
    ENDIF
    ALLOCATE(center(4, nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    CALL wannier_center(xyzmat,nstate,nstate,center,tau0)
    IF (paral%io_parent) CALL wc_print(nstate,center)

    ! Wannier function calculation
    ! Optimisers
    IF (wanni%w_opt.EQ.1) THEN
       IF (paral%io_parent) THEN
          ! steepest descent
          IF (cntl%tlsd) THEN
             CALL zeroing(rotmat)!,nstate*nstate)
             nxmax=MAX(spin_mod%nsup,spin_mod%nsdown)
             ix1=1
             ix2=ix1+nxmax*nxmax
             ixl=ix2+nxmax*nxmax
             CALL sd_wannier(rotlsd,xyzmat(1,1,1),nstate,spin_mod%nsup)
             CALL rmatmov(spin_mod%nsup,spin_mod%nsup,rotlsd,spin_mod%nsup,rotmat,nstate)
             CALL sd_wannier(rotlsd,xyzmat(spin_mod%nsup+1,spin_mod%nsup+1,1),nstate,&
                  spin_mod%nsdown)
             CALL rmatmov(spin_mod%nsdown,spin_mod%nsdown,rotlsd,spin_mod%nsdown,&
                  rotmat(spin_mod%nsup+1,spin_mod%nsup+1),nstate)
          ELSE
             ix1=1
             ix2=ix1+nstate*nstate
             ixl=ix2+nstate*nstate
             CALL sd_wannier(rotmat,xyzmat,nstate,nstate)
          ENDIF
       ENDIF
       CALL mp_bcast(rotmat,SIZE(rotmat),parai%io_source,parai%cp_grp)
    ELSEIF (wanni%w_opt.EQ.2) THEN
       ! Jacobi rotations
       IF (cntl%tlsd) THEN
          CALL zeroing(rotmat)!,nstate*nstate)
          nxmax=MAX(spin_mod%nsup,spin_mod%nsdown)
          IF(paral%io_parent) WRITE(6,*) procedureN//': localize spin up'
          CALL jrotation(rotlsd(1),xyzmat(1,1,1),nstate,spin_mod%nsup)
          CALL rmatmov(spin_mod%nsup,spin_mod%nsup,rotlsd,spin_mod%nsup,rotmat,nstate)
          IF(paral%io_parent) WRITE(6,*) procedureN//': localize spin down'
          CALL jrotation(rotlsd(1),xyzmat(spin_mod%nsup+1,spin_mod%nsup+1,1),nstate,spin_mod%nsdown)
          CALL rmatmov(spin_mod%nsdown,spin_mod%nsdown,rotlsd,spin_mod%nsdown,&
               rotmat(spin_mod%nsup+1,spin_mod%nsup+1),nstate)
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
    ! rotate the wavefunctions to the Wannier representation
    CALL rotate(1._real_8,c0,0._real_8,c2,rotmat,nstate,2*ncpw%ngw,cntl%tlsd,spin_mod%nsup,spin_mod%nsdown)
    CALL dcopy(2*ncpw%ngw*nstate,c2,1,c0,1)
    CALL rhoe_psi_size(il_psi_1d=il_psi_1d,il_psi_2d=il_psi_2d,&
         il_rhoe_1d=il_rhoe_1d,il_rhoe_2d=il_rhoe_2d)
    ALLOCATE(psi(il_psi_1d,il_psi_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(rhoe(il_rhoe_1d,il_rhoe_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(fion(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(fion)!,3*maxsys%nax*maxsys%nsx)
    IF (wannl%twmol.OR.wannl%twdos) THEN
       ALLOCATE(eigv(nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL forcedr(c0,c2,sc0,rhoe,psi,tau0,fion,eigv,nstate,&
            1,.FALSE.,.TRUE.)
    ENDIF
    ! calculate the centers and spread of the wannier functions
    CALL wannier_center(xyzmat,nstate,nstate,center,tau0)
    IF (paral%io_parent) CALL wc_print(nstate,center)
    IF (cntl%tresponse) CALL wannier_save(nstate,center)
    !CSOC[
    IF (cntl%tddft.OR.cntl%tspec.OR.cntl%tsoc) THEN
       !CSOC]
       IF (td03%treorder.OR.td03%molstat) THEN
          CALL dcopy(nstate,center(1,1),4,wcent(1,1),3)
          CALL dcopy(nstate,center(2,1),4,wcent(2,1),3)
          CALL dcopy(nstate,center(3,1),4,wcent(3,1),3)
       ENDIF
    ENDIF
    IF (wannl%twmol) THEN
       CALL molorb(c0,c2,tau0,nstate,center)
       CALL forcedr(c0,c2,sc0,rhoe,psi,tau0,fion,eigv,nstate,&
            1,.FALSE.,.TRUE.)
    ENDIF
    IF (wannl%twdos) CALL wc_dos(c0,c2,nstate,center)

    ! print selected Wannier functions
    CALL wannier_print(1,c0,tau0,nstate,psi(:,1),center)

    ! save rotmat for voa
    IF(response1%tvoa) THEN
       voa_data%rotmat=rotmat
    ENDIF

    DEALLOCATE(center,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)

    DEALLOCATE(psi,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(mapful,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(mapcol,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(xyzmat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rotmat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (cntl%tlsd) DEALLOCATE(rotlsd,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE localize
  ! ==================================================================
  SUBROUTINE localize2(tau0,c0,c2,sc0,nstate)
    ! ==--------------------------------------------------------------==
    ! == CALCULATES LOCALIZED FUNCTIONS                               ==
    ! ==--------------------------------------------------------------==
    ! input
    REAL(real_8)                             :: tau0(:,:,:)
    COMPLEX(real_8)                          :: sc0(*)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c2(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'localize2'

    COMPLEX(real_8)                          :: dd1, dd2, dd3
    COMPLEX(real_8), ALLOCATABLE             :: ddmat(:,:)
    COMPLEX(real_8), ALLOCATABLE, &
      DIMENSION(:, :)                        :: lbasis
    COMPLEX(real_8), ALLOCATABLE, SAVE       :: xyzmat(:,:,:)
    INTEGER                                  :: i, ia, iaorb, iat, ierr, &
                                                info, is, isub, ix1, ix2, &
                                                ixl, ixx, n1, natst, nmcol, &
                                                nxmax, w_opt_tmp
    INTEGER, ALLOCATABLE, SAVE               :: mapcol(:), mapful(:)
    INTEGER, SAVE                            :: i_comp = 0, ifirst = 0
    LOGICAL                                  :: wannier_recomput
    REAL(real_8)                             :: fac, &
                                                foc(loadc_foc_array_size), &
                                                omegon, pmaton, sfc
    REAL(real_8), ALLOCATABLE                :: center(:,:)
    REAL(real_8), ALLOCATABLE, DIMENSION(:)  :: ovl, s, work
    REAL(real_8), ALLOCATABLE, &
      DIMENSION(:, :)                        :: u, vt
    REAL(real_8), ALLOCATABLE, SAVE          :: rotlsd(:), rotmat(:,:)

    CALL tiset(procedureN,isub)
    w_opt_tmp=wanni%w_opt
    IF(MOD(i_comp,wan05%loc_relocalize_every)==0.AND.wan05%loc_relocalize) wanni%w_opt=2

    wannier_recomput=MOD(i_comp,wan05%loc_recompute_dipole_matrices_every)==0
    i_comp=i_comp+1
    IF (paral%io_parent.AND.wannier_recomput) &
         WRITE(6,*) procedureN//': RECOMPUTE DIPOLE MATRICES'
    IF (ifirst.EQ.0) THEN
       ifirst=1
       ! ..initialization
       n1=0
       DO i=0,parai%nproc-1
          n1=MAX(n1,parap%sparm(3,i))
       ENDDO
       ngwmax=n1
       ALLOCATE(mapful(2*spar%ngws),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       nmcol=parai%nproc*ngwmax
       ALLOCATE(mapcol(nmcol),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL setdip(mapful,mapcol)
       ! ..initialization for Wannier function part
       CALL set_operator(.TRUE.)
       ALLOCATE(xyzmat(nstate,nstate,wannc%nwanopt),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(rotmat(nstate,nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(xyzmat)!,nstate*nstate*nwanopt)
       IF (cntl%tlsd) THEN
          nxmax=MAX(spin_mod%nsup,spin_mod%nsdown)
          ALLOCATE(rotlsd(nxmax*nxmax),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
    ENDIF
    IF (cntl%tlsd) THEN
       fac=1._real_8
    ELSE
       fac=2._real_8
    ENDIF
    pmaton=1.0_real_8/rmass%pmat0
    omegon=1.0_real_8/parm%omega
    ! ..electronic contribution
    IF (wannier_recomput) THEN
       ALLOCATE(ddmat(nstate,nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       CALL opeigr(c0,c2,sc0,nstate,mapful,mapcol,ddmat,1,0,dd1)
       CALL dcopy(2*nstate*nstate,ddmat,1,xyzmat(1,1,1),1)
       CALL opeigr(c0,c2,sc0,nstate,mapful,mapcol,ddmat,2,0,dd2)
       CALL dcopy(2*nstate*nstate,ddmat,1,xyzmat(1,1,2),1)
       CALL opeigr(c0,c2,sc0,nstate,mapful,mapcol,ddmat,3,0,dd3)
       CALL dcopy(2*nstate*nstate,ddmat,1,xyzmat(1,1,3),1)
       DO i=4,wannc%nwanopt
          CALL opeigr(c0,c2,sc0,nstate,mapful,mapcol,ddmat,&
               wannc%iow(1,i-3),wannc%iow(2,i-3),dd1)
          CALL dcopy(2*nstate*nstate,ddmat,1,xyzmat(1,1,i),1)
       ENDDO
       DEALLOCATE(ddmat,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    ENDIF

    ! Wannier function calculation
    ! Optimisers
    IF (wanni%w_opt.EQ.1) THEN
       IF (paral%io_parent) THEN
          ! steepest descent
          IF (cntl%tlsd) THEN
             CALL zeroing(rotmat)!,nstate*nstate)
             nxmax=MAX(spin_mod%nsup,spin_mod%nsdown)
             ix1=1
             ix2=ix1+nxmax*nxmax
             ixl=ix2+nxmax*nxmax
             CALL sd_wannier(rotlsd,xyzmat(1,1,1),nstate,spin_mod%nsup)
             CALL rmatmov(spin_mod%nsup,spin_mod%nsup,rotlsd,spin_mod%nsup,rotmat,nstate)
             CALL sd_wannier(rotlsd,xyzmat(spin_mod%nsup+1,spin_mod%nsup+1,1),nstate,&
                  spin_mod%nsdown)
             CALL rmatmov(spin_mod%nsdown,spin_mod%nsdown,rotlsd,spin_mod%nsdown,&
                  rotmat(spin_mod%nsup+1,spin_mod%nsup+1),nstate)
          ELSE
             ix1=1
             ix2=ix1+nstate*nstate
             ixl=ix2+nstate*nstate
             CALL sd_wannier(rotmat,xyzmat,nstate,nstate)
          ENDIF
       ENDIF
       CALL mp_bcast(rotmat,SIZE(rotmat),parai%io_source,parai%cp_grp)
    ELSEIF (wanni%w_opt.EQ.2) THEN
       ! Jacobi rotations
       IF (cntl%tlsd) THEN
          CALL zeroing(rotmat)!,nstate*nstate)
          nxmax=MAX(spin_mod%nsup,spin_mod%nsdown)
          IF(paral%io_parent) WRITE(6,*) procedureN//': localize spin up'
          CALL jrotation(rotlsd,xyzmat(1,1,1),nstate,spin_mod%nsup)
          CALL rmatmov(spin_mod%nsup,spin_mod%nsup,rotlsd,spin_mod%nsup,rotmat,nstate)
          IF(paral%io_parent) WRITE(6,*) procedureN//': localize spin down'
          CALL jrotation(rotlsd,xyzmat(spin_mod%nsup+1,spin_mod%nsup+1,1),nstate,spin_mod%nsdown)
          CALL rmatmov(spin_mod%nsdown,spin_mod%nsdown,rotlsd,spin_mod%nsdown,&
               rotmat(spin_mod%nsup+1,spin_mod%nsup+1) ,nstate)
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
    ! rotate the wavefunctions to the Wannier representation
    ! IF(cp_nogrp.GT.1) THEN
    ! !quickfix: we need to make sure that each groups see the same rotation matrix
    ! !need to modify jrotation to use cp_grp
    ! MSGLEN=NSTATE*NSTATE * 8
    ! CALL MY_BCAST(ROTMAT,MSGLEN,0,CP_INTER_GRP)
    ! ENDIF
    CALL rotate(1._real_8,c0,0._real_8,c2,rotmat,nstate,&
         2*ncpw%ngw,cntl%tlsd,spin_mod%nsup,spin_mod%nsdown)
    CALL dcopy(2*ncpw%ngw*nstate,c2,1,c0,1)
    ! calculate the centers and spread of the wannier functions
    ALLOCATE(center(4,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL wannier_center(xyzmat,nstate,nstate,center,tau0)
    DEALLOCATE(center,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    wanni%w_opt=w_opt_tmp
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE localize2
  ! ==================================================================
  SUBROUTINE wc_print(nstate,center)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    REAL(real_8)                             :: center(:,:)

    INTEGER                                  :: i
    REAL(real_8)                             :: sum

!(4,nstate)
! ==--------------------------------------------------------------==

    WRITE(6,'(/,1X,64("*"))')
    WRITE(6,'(" *",12X,A,18X,A,1X,"*")')&
         ' WANNIER CENTERS  ','<X^2> - <X>^2'
    WRITE(6,'(1X,64("*"))')
    sum = 0.0_real_8
    DO i=1,nstate
       WRITE(6,'(3F12.4,17X,F12.4)')&
            center(1,i),center(2,i),center(3,i),center(4,i)
       sum = sum+center(4,i)
    ENDDO
    WRITE(6,'(/,1X,"  TOTAL SPREAD OF THE SYSTEM ",F12.4)') suM
    WRITE(6,'(1X,64("*"),/)')
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE wc_print
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE wannier_save(nstate,center)
    ! saves the wannier centers for use in the NMR/EPR routines. ONLY
    ! possible when called from NMR_P/EPR_P!!! 
    ! (Needs memory to be allocated)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    REAL(real_8)                             :: center(:,:)

    INTEGER                                  :: center_coord, is
    REAL(real_8)                             :: projection

!(4,nstate)
! ==--------------------------------------------------------------==

    IF (.NOT. (response1%tnmr.OR.response1%tepr)) RETURN

    IF (paral%parent) THEN
       IF (response1%tnmr .OR. response1%tepr) THEN
          CALL dcopy(4*nstate,center,1,wanniercenters,1)
          DO is=1,nstate
             ! FIRST DIMENSION (nr1s)
             projection = center(1,is) * gvec_com%b1(1)  +&
                  center(2,is) * gvec_com%b1(2)+&
                  center(3,is) * gvec_com%b1(3)
             projection = projection / parm%alat * REAL(spar%nr1s,kind=real_8) + 1
             center_coord = bound( NINT(projection) ,spar%nr1s)
             lower_left(1,is) = bound(center_coord - spar%nr1s/2,spar%nr1s)
             lower_left_value(1,is) = - spar%nr1s/2
             ! SECOND DIMENSION (nr2s)
             projection = center(1,is) * gvec_com%b2(1)  +&
                  center(2,is) * gvec_com%b2(2)+&
                  center(3,is) * gvec_com%b2(3)
             projection = projection / parm%alat * REAL(spar%nr2s,kind=real_8) + 1
             center_coord = bound( NINT(projection) ,spar%nr2s)
             lower_left(2,is) = bound(center_coord - spar%nr2s/2,spar%nr2s)
             lower_left_value(2,is) = - spar%nr2s/2
             ! THIRD DIMENSION (nr3s)
             projection = center(1,is) * gvec_com%b3(1)  +&
                  center(2,is) * gvec_com%b3(2)+&
                  center(3,is) * gvec_com%b3(3)
             projection = projection / parm%alat * REAL(spar%nr3s,kind=real_8) + 1
             center_coord = bound( NINT(projection) ,spar%nr3s)
             lower_left(3,is) = bound(center_coord - spar%nr3s/2,spar%nr3s)
             lower_left_value(3,is) = - spar%nr3s/2
          ENDDO
       ENDIF
    ENDIF

    CALL mp_bcast(lower_left_value,SIZE(lower_left_value),parai%source,parai%allgrp)
    CALL mp_bcast(lower_left,SIZE(lower_left),parai%source,parai%allgrp)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE wannier_save
  ! ==================================================================


END MODULE localize_utils
