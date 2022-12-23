MODULE printp_utils
  USE adat,                            ONLY: elem
  USE cell,                            ONLY: cell_com
  USE clas,                            ONLY: clas3,&
                                             clas4,&
                                             clas8,&
                                             claspres,&
                                             clasv
  USE cnst,                            ONLY: au_kb,&
                                             fbohr,&
                                             scmass
  USE coninp_utils,                    ONLY: raddeg
  USE cotr,                            ONLY: &
       cotc0, cotr007, duat, fv, ntcnst, ntrest, resfdiff, resm, resv, &
       resval, xlagr
  USE ddip,                            ONLY: pdipole,&
                                             pdipolt
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_app,&
                                             fo_def,&
                                             fo_ufo,&
                                             fo_verb,&
                                             fo_vmark
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_4,&
                                             real_8
  USE meta_multiple_walkers_utils,     ONLY: mw_filename
  USE metr,                            ONLY: metr_com
  USE mfep,                            ONLY: mfepi
  USE mm_dim_utils,                    ONLY: mm_dim
  USE mm_dimmod,                       ONLY: clsaabox,&
                                             mm_go_mm,&
                                             mm_revert,&
                                             nat_grm
  USE movi,                            ONLY: imtyp
  USE mw,                              ONLY: mwi,&
                                             tmw
  USE parac,                           ONLY: paral
  USE pimd,                            ONLY: ipcurr
  USE printpmod,                       ONLY: if1,&
                                             if2,&
                                             if3
  USE readsr_utils,                    ONLY: xstring
  USE rmas,                            ONLY: rmass
  USE ropt,                            ONLY: iteropt,&
                                             ropt_mod
  USE store_types,                     ONLY: cprint,&
                                             rout1,&
                                             trajsmall,&
                                             trajsmalln
  USE strs,                            ONLY: paiu
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             cntr,&
                                             iatpt,&
                                             parm
  use bicanonicalCpmd, only: bicanonicalCpmdConfig,& 
  biCanonicalEnsembleDo, getNameConstraintTape, getNameFtrajectoryTape,&
  getNamemetricTape, getNameTrajecTape, getNameTrajectoryTape

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: printp
  PUBLIC :: printp2
  PUBLIC :: print_mts_forces

CONTAINS

  ! ==================================================================
  SUBROUTINE printp(taur,taup,velp)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: taur(:,:,:), taup(:,:,:), &
                                                velp(:,:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'printp'
    CHARACTER(len=4), PARAMETER              :: dcdtype = 'CORD'
    CHARACTER(len=80), PARAMETER :: &
      dcdtitle = 'CHARMM format DCD trajectory written by CPMD. '

    CHARACTER(len=20) :: f1 = 'MOVIE     ', f10 = 'METRIC    ', &
      f2 = 'TRAJECTORY', f3 = 'CONSTRAINT', f4 = 'STRESS    ', &
      f5 = 'CELL      ', f6 = 'DIPOLE    ', f7 = 'STRECL    ', &
      f8 = 'TRAJEC.xyz', f9 = 'TRAJEC.dcd'
    CHARACTER(len=4)                         :: cstring
    CHARACTER(len=80)                        :: ftmp
    INTEGER                                  :: i, i1, i2, ia, ic, id, ierr, &
                                                is, it0, ityp, iz0, j, k, l, &
                                                naa
    INTEGER, ALLOCATABLE                     :: gr_iat(:)
    INTEGER, SAVE                            :: ifirst = 1
    LOGICAL                                  :: ferror, status
    REAL(real_4)                             :: dcddelta
    REAL(real_4), ALLOCATABLE                :: dcdx(:), dcdy(:), dcdz(:)
    REAL(real_8)                             :: const, cval, dcdcell(6), &
                                                fact, fval, out(3,3), x(3)
    REAL(real_8), ALLOCATABLE                :: gr_tau(:,:), gr_vel(:,:)

    IF (ifirst.EQ.1) THEN
       if1=fo_verb
       if2=fo_vmark
       if3=fo_app+fo_vmark
       ifirst=0
    ENDIF
    ! Setup file-open mode and file-names for path minimization
    IF (cntl%tpmin) THEN
       IF (MOD(iteropt%nfi,cnti%nomore).EQ.1) THEN
          if3=fo_def
       ELSE
          if3=fo_app
       ENDIF
       ftmp='TRAJECTORY_'
       CALL mw_filename(ftmp,f2,ipcurr)
       CALL xstring(f2,i1,i2)
       ftmp=f2(i1:i2)//'.'
       CALL mw_filename(ftmp,f2,mfepi%istring)
       ftmp='TRAJEC_'
       CALL mw_filename(ftmp,f8,ipcurr)
       CALL xstring(f8,i1,i2)
       ftmp=f8(i1:i2)//'.'
       CALL mw_filename(ftmp,f8,mfepi%istring)
       CALL xstring(f8,i1,i2)
       ftmp=f8(i1:i2)//".xyz"
       f8=ftmp
       ftmp='CONSTRAINT_'
       CALL mw_filename(ftmp,f3,ipcurr)
       CALL xstring(f3,i1,i2)
       ftmp=f3(i1:i2)//'.'
       CALL mw_filename(ftmp,f3,mfepi%istring)
       ftmp='METRIC_'
       CALL mw_filename(ftmp,f10,ipcurr)
       CALL xstring(f10,i1,i2)
       ftmp=f10(i1:i2)//'.'
       CALL mw_filename(ftmp,f10,mfepi%istring)
    ENDIF

    CALL mm_dim(mm_go_mm,status)
    ALLOCATE(gr_tau(3,ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(gr_vel(3,ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(gr_iat(ions1%nat+1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! If multiple walkers TRAJECTORY to TRAJECTORY_IP and TRAJEC.xyz to TRAJEC_IP.xyz
    IF (tmw)THEN
       ftmp='TRAJECTORY_'
       CALL mw_filename(ftmp,f2,mwi%walker_id)
       ftmp='TRAJEC_'
       CALL mw_filename(ftmp,f8,mwi%walker_id)
       CALL xstring(f8,i1,i2)
       ftmp=f8(i1:i2)//".xyz"
       f8=ftmp
       ftmp='CONSTRAINT_'
       CALL mw_filename(ftmp,f3,mwi%walker_id)
       ftmp='METRIC_'
       CALL mw_filename(ftmp,f10,mwi%walker_id)
    ENDIF
    if (biCanonicalEnsembleDo) then 
       f2 = getNameTrajectoryTape(bicanonicalCpmdConfig)
       f8 = getNameTrajecTape(bicanonicalCpmdConfig)
       f3 = getNameConstraintTape(bicanonicalCpmdConfig)
       f10 = getNameMetricTape(bicanonicalCpmdConfig)
    ENDIF

    ! ==--------------------------------------------------------------==
    ! ==  OPEN THE TRAJECTORY AND MOVIE FILES                         ==
    ! ==--------------------------------------------------------------==
    IF (ropt_mod%rprint) THEN
       IF (paral%io_parent)&
            CALL fileopen(4,f2,fo_app+if2,ferror)
    ENDIF
    IF (ropt_mod%rprint.OR.(ropt_mod%txyz.AND.rout1%xtout).OR.(ropt_mod%tdcd.AND.rout1%dcout)) THEN
       IF (cotc0%mcnstr.GT.0.OR.cotr007%mrestr.GT.0) THEN
          IF (paral%io_parent)&
               CALL fileopen(31,f3,if3,ferror)
          IF ((cotr007%mrestr.GT.0).AND.paral%io_parent)&
               CALL fileopen(90,f10,if3,ferror)
          if3=fo_app
       ENDIF
    ENDIF
    ! CHECK IF WE NEED TO WRITE A DCD HEADER
    ! FIXME: ideally we would also check if an existing dcd file
    ! matches the current system and rename/stop if it does not.
    ! FIXME: also we don't update the number of steps in the dcd file.
    ! AK 2005/08/21
    IF (rout1%dcout) THEN
       ferror=.TRUE.
       IF (paral%io_parent)&
            INQUIRE(file=f9,exist=ferror)
       IF (.NOT.ferror) THEN
          IF (paral%io_parent)&
               CALL fileopen(35,f9,fo_app+fo_ufo,ferror)
          dcddelta=cntr%delt_ions
          IF (paral%io_parent)&
               WRITE(35) dcdtype,0,iteropt%nfi,cnti%ntraj,0,0,0,0,0,0,dcddelta,1,&
               0,0,0,0,0,0,0,0,24
          IF (paral%io_parent)&
               WRITE(35) 2,dcdtitle,dcdtitle
          IF (paral%io_parent)&
               WRITE(35) ions1%nat
          IF (paral%io_parent)&
               CALL fileclose(35)
       ENDIF
    ENDIF
    IF ((ropt_mod%txyz.AND.rout1%xtout).AND.paral%io_parent)&
         CALL fileopen( 8,f8,fo_app+if1,ferror)
    IF ((ropt_mod%tdcd.AND.rout1%dcout).AND.paral%io_parent)&
         CALL fileopen(35,f9,fo_app+fo_ufo+if1,ferror)
    IF ((ropt_mod%movie).AND.paral%io_parent)&
         CALL fileopen(11,f1,fo_app+if2,ferror)
    IF ((ropt_mod%calstc).AND.paral%io_parent)&
         CALL fileopen(55,f7,fo_app+if2,ferror)
    IF ((cntl%caldip).AND.paral%io_parent)&
         CALL fileopen(34,f6,fo_app+if2,ferror)
    ! ..Verbose+Mark only on first encounter.
    if1=0
    if2=0
    ! get coordinates and velocities into gromos ordering
    i=0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          i=i+1
          gr_iat(nat_grm(i))=ions0%iatyp(is)
          DO k=1,3
             gr_tau(k,nat_grm(i))=taup(k,ia,is)-clsaabox%mm_c_trans(k)
             gr_vel(k,nat_grm(i))=velp(k,ia,is)
          ENDDO
       ENDDO
    ENDDO
    ! ..Store ionic coordinates and velocities for statistics
    IF (ropt_mod%rprint) THEN
       IF (cntl%tsampl) THEN
          DO i=cprint%minwriteatom,cprint%maxwriteatom
             IF (.NOT.trajsmall .OR. i.LE.trajsmalln) THEN
                IF (cprint%twritebintrajectory) THEN
                   IF (paral%io_parent)&
                        WRITE(4) iteropt%nfi,(gr_tau(k,i),k=1,3),(gr_vel(k,i),k=1,3)
                ELSE
                   IF (paral%io_parent)&
                        WRITE(4,'(I7,6(2X,F22.14))')&
                        iteropt%nfi,(gr_tau(k,i),k=1,3),(gr_vel(k,i),k=1,3)
                ENDIF
             ENDIF
          ENDDO
       ELSE
          DO i=cprint%minwriteatom,cprint%maxwriteatom
             IF (.NOT.trajsmall .OR. i.LE.trajsmalln) THEN
                IF (cprint%twritebintrajectory) THEN
                   IF (paral%io_parent)&
                        WRITE(4) iteropt%nfi,(gr_tau(k,i),k=1,3),(gr_vel(k,i),k=1,3)
                ELSE
                   IF (paral%io_parent)&
                        WRITE(4,'(I7,6(2X,F22.14))')&
                        iteropt%nfi,(gr_tau(k,i),k=1,3),(gr_vel(k,i),k=1,3)
                ENDIF
             ENDIF
          ENDDO
       ENDIF
       IF (paral%io_parent)&
            CALL fileclose(4)
    ENDIF
    IF (ropt_mod%rprint.OR.(ropt_mod%txyz.AND.rout1%xtout).OR.(ropt_mod%tdcd.AND.rout1%dcout)) THEN
       IF (cntl%tprcp.AND..NOT.(cntl%tpath.AND.cntl%tpimd)) THEN
          IF (paral%io_parent)&
               CALL fileopen(32,f5,fo_app+if2,ferror)
          IF (paral%io_parent)&
               WRITE(32,*) '  CELL PARAMETERS at Step:', iteropt%nfi
          DO i=1,3
             IF (paral%io_parent)&
                  WRITE(32,'(3(1X,F14.6),8X,3(1X,F12.6))')&
                  (metr_com%ht(i,j),j=1,3),(metr_com%htvel(i,j),j=1,3)
          ENDDO
          IF (paral%io_parent)&
               CALL fileclose(32)
       ENDIF
       IF (cotc0%mcnstr.GT.0.OR.cotr007%mrestr.GT.0) THEN
          DO j=1,cotc0%mcnstr
             fval=fv(j)
             ityp=ntcnst(1,j)
             IF (ityp.EQ.2.OR.ityp.EQ.3.OR.ityp.EQ.5) CALL raddeg(fval,1)
             IF (paral%io_parent)&
                  WRITE(31,'(I7,2X,I4,5X,2(1PE20.10))') iteropt%nfi,j,xlagr(j),fval
          ENDDO
          DO j=1,cotr007%mrestr
             fval=resv(j)
             cval=resval(j)
             ityp=ntrest(1,j)
             IF (ityp.EQ.2.OR.ityp.EQ.3.OR.ityp.EQ.5) CALL raddeg(fval,1)
             IF (ityp.EQ.2.OR.ityp.EQ.3.OR.ityp.EQ.5) CALL raddeg(cval,1)
             fval=fval-cval
             IF (paral%io_parent)&
                  WRITE(31,'(I7,2X,I4,5X,3(1PE20.10)," R")') iteropt%nfi,j,cval,fval,&
                  resfdiff(j)
          ENDDO
          IF (paral%io_parent) WRITE(cstring,'(I4)') cotr007%mrestr
          CALL xstring(cstring,i1,i2)
          DO i=1,cotr007%mrestr
             IF (paral%io_parent)&
                  WRITE(90,'('//cstring(i1:i2)//'(1PE20.10))')&
                  (resm(i,j),j=1,cotr007%mrestr)
          ENDDO
          IF (paral%io_parent)&
               CALL fileclose(31)
          IF ((cotr007%mrestr.GT.0).AND.paral%io_parent)&
               CALL fileclose(90)
       ENDIF
    ENDIF
    ! ..Write xyz file output
    IF (ropt_mod%txyz) THEN
       IF (paral%io_parent)&
            WRITE(8,*) ions1%nat+duat%ndat
       IF (paral%io_parent)&
            WRITE(8,*) 'STEP:', iteropt%nfi
       DO i=1,ions1%nat
          IF (paral%io_parent)&
               WRITE(8,'(A2,3(2X,F12.6))')&
               elem%el(gr_iat(i)),(gr_tau(k,i)/fbohr,k=1,3)
       ENDDO
       ! Dummy atoms (type 1-4)
       DO i=ions1%nat+1,ions1%nat+duat%ndat
          naa=i-ions1%nat
          ityp=duat%listda(naa,1)
          id=duat%listda(naa,2)
          IF (ityp.EQ.4) THEN
             x(1)=duat%dummy4(1,id)
             x(2)=duat%dummy4(2,id)
             x(3)=duat%dummy4(3,id)
          ELSEIF (ityp.EQ.3) THEN
             x(1)=duat%dummy3(1,id)
             x(2)=duat%dummy3(2,id)
             x(3)=duat%dummy3(3,id)
          ELSEIF (ityp.EQ.2) THEN
             x(1)=duat%dummy2(1,id)
             x(2)=duat%dummy2(2,id)
             x(3)=duat%dummy2(3,id)
          ELSEIF (ityp.EQ.1) THEN
             x(1)=duat%dummy1(1,id)
             x(2)=duat%dummy1(2,id)
             x(3)=duat%dummy1(3,id)
          ENDIF
          IF (paral%io_parent)&
               WRITE(8,'(A2,3(2X,F12.6))')&
               'X',(x(k)/fbohr,k=1,3)
       ENDDO
       IF (paral%io_parent)&
            CALL fileclose(8)
    ENDIF
    ! ..Write dcd file output
    ! ..FIXME: we don't update the number of frames field in the header (yet).
    IF (ropt_mod%tdcd) THEN
       ALLOCATE(dcdx(ions1%nat),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(dcdy(ions1%nat),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(dcdz(ions1%nat),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       DO i=1,ions1%nat
          dcdx(i)=gr_tau(1,i)/fbohr
          dcdy(i)=gr_tau(2,i)/fbohr
          dcdz(i)=gr_tau(3,i)/fbohr
       ENDDO
       ! store unitcell. use classical box for QM/MM
       IF (cntl%tqmmm) THEN
          dcdcell(1)=clsaabox%box_au(1)/fbohr
          dcdcell(3)=clsaabox%box_au(2)/fbohr
          dcdcell(6)=clsaabox%box_au(3)/fbohr
          dcdcell(2)=0.0_real_8
          dcdcell(4)=0.0_real_8
          dcdcell(5)=0.0_real_8
       ELSE
          ! FIXME: this does not work for variable cell calculations.
          ! we overried it later.
          dcdcell(1)=cell_com%celldm(1)/fbohr
          dcdcell(3)=cell_com%celldm(1)*cell_com%celldm(2)/fbohr
          dcdcell(6)=cell_com%celldm(1)*cell_com%celldm(3)/fbohr
          dcdcell(2)=cell_com%celldm(4)
          dcdcell(4)=cell_com%celldm(5)
          dcdcell(5)=cell_com%celldm(6)

          IF (cntl%tprcp) THEN
             dcdcell(1)=SQRT(metr_com%ht(1,1)**2 + metr_com%ht(1,2)**2 + metr_com%ht(1,3)**2)/fbohr
             dcdcell(3)=SQRT(metr_com%ht(2,1)**2 + metr_com%ht(2,2)**2 + metr_com%ht(2,3)**2)/fbohr
             dcdcell(6)=SQRT(metr_com%ht(3,1)**2 + metr_com%ht(3,2)**2 + metr_com%ht(3,3)**2)/fbohr
             dcdcell(2)=(metr_com%ht(1,1)*metr_com%ht(2,1)+metr_com%ht(1,2)*metr_com%ht(2,2)+metr_com%ht(1,3)*metr_com%ht(2,3))&
                  /(dcdcell(1)*dcdcell(3))
             dcdcell(4)=(metr_com%ht(1,1)*metr_com%ht(3,1)+metr_com%ht(1,2)*metr_com%ht(3,2)+metr_com%ht(1,3)*metr_com%ht(3,3))&
                  /(dcdcell(1)*dcdcell(6))
             dcdcell(5)=(metr_com%ht(2,1)*metr_com%ht(3,1)+metr_com%ht(2,2)*metr_com%ht(3,2)+metr_com%ht(2,3)*metr_com%ht(3,3))&
                  /(dcdcell(3)*dcdcell(6))
          ENDIF
       ENDIF
       IF (paral%io_parent)&
            WRITE(35) (dcdcell(i),i=1,6)
       IF (paral%io_parent)&
            WRITE(35) (dcdx(i),i=1,ions1%nat)
       IF (paral%io_parent)&
            WRITE(35) (dcdy(i),i=1,ions1%nat)
       IF (paral%io_parent)&
            WRITE(35) (dcdz(i),i=1,ions1%nat)
       DEALLOCATE(dcdx,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(dcdy,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(dcdz,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       IF (paral%io_parent)&
            CALL fileclose(35)
    ENDIF
    ! ..Write Movie File
    ! FIXME: AK 2005/05/25 not yet converted to GROMOS ordering.
    IF (ropt_mod%movie)THEN
       DO is=1,ions1%nsp
          iz0=imtyp(is)
          DO ia=1,ions0%na(is)
             IF (ions0%iatyp(is).EQ.4) THEN
                ! ..Change Be -> Na ; problem with movie
                it0=11
                IF (paral%io_parent)&
                     WRITE(11,'(3(2X,F12.4),2I4)')&
                     (taup(k,ia,is)/fbohr,k=1,3),it0,iz0
             ELSE
                IF (paral%io_parent)&
                     WRITE(11,'(3(2X,F12.4),2I4)')&
                     (taup(k,ia,is)/fbohr,k=1,3),ions0%iatyp(is),iz0
             ENDIF
          ENDDO
       ENDDO
       IF (paral%io_parent)&
            CALL fileclose(11)
    ENDIF
    ! ..Write Stress tensor
    IF (.NOT.(cntl%tpath.AND.cntl%tpimd).AND.ropt_mod%calste.AND.(MOD(iteropt%nfi-1,cnti%npres).EQ.0)) THEN
       CALL dcopy(9,paiu,1,out,1)
       !$omp parallel do private(I,IA,IS,FACT) reduction(+:OUT)
       DO i=1,ions1%nat
          ia=iatpt(1,i)
          is=iatpt(2,i)
          fact=rmass%pma(is)
          out(1,1)=out(1,1)+fact*velp(1,ia,is)*velp(1,ia,is)
          out(2,1)=out(2,1)+fact*velp(2,ia,is)*velp(1,ia,is)
          out(3,1)=out(3,1)+fact*velp(3,ia,is)*velp(1,ia,is)
          out(1,2)=out(1,2)+fact*velp(1,ia,is)*velp(2,ia,is)
          out(2,2)=out(2,2)+fact*velp(2,ia,is)*velp(2,ia,is)
          out(3,2)=out(3,2)+fact*velp(3,ia,is)*velp(2,ia,is)
          out(1,3)=out(1,3)+fact*velp(1,ia,is)*velp(3,ia,is)
          out(2,3)=out(2,3)+fact*velp(2,ia,is)*velp(3,ia,is)
          out(3,3)=out(3,3)+fact*velp(3,ia,is)*velp(3,ia,is)
       ENDDO
       ! We give th true total stress tensor. (T.D.)
       ! DO I=1,3
       ! OUT(I,I)=OUT(I,I)+DRUCK*OMEGA
       ! ENDDO
       IF (paral%io_parent)&
            CALL fileopen(33,f4,fo_app+if2,ferror)
       IF (paral%io_parent)&
            WRITE(33,*) '  TOTAL STRESS TENSOR (kB): Step:', iteropt%nfi
       DO i=1,3
          IF (paral%io_parent)&
               WRITE(33,'(5X,3(F20.8))') ((out(i,j)/parm%omega)*au_kb,j=1,3)
       ENDDO
       IF (paral%io_parent)&
            CALL fileclose(33)
    ENDIF
    ! WRITE CLASSICAL STRESS TENSOR
    ! we add the ideal gas contribution of the classical particles
    IF (ropt_mod%calstc) THEN
       DO ic=1,clas3%ncltyp
          IF (clas3%is_qm(ic).EQ.0) THEN
             ! classical particle
             const=clas4%clmas(ic)*scmass
             DO ia=clas3%ncrang(1,ic),clas3%ncrang(2,ic)
                DO k=1,3
                   DO l=1,3
                      claspres(k,l)=claspres(k,l)+const*clasv(k,ia)*clasv(l,ia)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDDO
       IF (paral%io_parent)&
            WRITE(55,*) 'CLASSICAL STRESS TENSOR (kB): Step:', iteropt%nfi
       ! WRITE(57,*) 'CLASSICAL STRESS TENSOR (kB): Step:', NFI
       DO i=1,3
          IF (paral%io_parent)&
               WRITE(55,'(5X,3(F20.8))') ((claspres(i,j)/clas8%clomega)*au_kb,j=1,3)
          ! WRITE(57,'(5X,3(F20.8))') ((claspres(I,J)/CLOMEGA)*AU_KB,J=1,3)
       ENDDO
       IF (paral%io_parent)&
            CALL fileclose(55)
    ENDIF
    ! ..Write Dipole Moments
    IF (cntl%caldip) THEN
       IF (paral%io_parent)&
            WRITE(34,'(I7,6(2X,1PE16.6E2))')iteropt%nfi,(pdipole(i),i=1,3),&
            (pdipolt(i),i=1,3)
       IF (paral%io_parent)&
            CALL fileclose(34)
    ENDIF
    CALL mm_dim(mm_revert,status)
    DEALLOCATE(gr_tau,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(gr_vel,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(gr_iat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE printp
  ! ==================================================================


  ! ==================================================================
  SUBROUTINE printp2(taur,taup,velp,fion)
    ! ==--------------------------------------------------------------==

    REAL(real_8)                             :: taur(:,:,:), taup(:,:,:), &
                                                velp(:,:,:), fion(:,:,:)

    CHARACTER(len=11), PARAMETER             :: f2default = 'FTRAJECTORY'
    CHARACTER(len=20) :: f2

    INTEGER                                  :: j, k, l
    INTEGER, SAVE                            :: if2 = fo_vmark
    LOGICAL                                  :: ferror

! ==--------------------------------------------------------------==
! ==  OPEN THE TRAJECTORY AND MOVIE FILES                         ==
! ==--------------------------------------------------------------==
    
    if (biCanonicalEnsembleDo) then 
       f2 = getNameFtrajectoryTape(bicanonicalCpmdConfig)
     else
       f2 = f2default
    ENDIF
    IF (ropt_mod%rprint) THEN
       IF (paral%io_parent)&
            CALL fileopen(4,f2,fo_app+if2,ferror)
    ENDIF
    if2=0
    ! ..Store ionic coordinates and velocities for statistics
    ! FIXME: AK 2005/05/25 not yet converted to GROMOS ordering.
    IF (ropt_mod%rprint) THEN
       IF (cntl%tsampl) THEN
          l=0
          DO k=1,ions1%nsp
             DO j=1,ions0%na(k)
                l=l+1
                IF (l.GE.cprint%minwriteatom.AND.l.LE.cprint%maxwriteatom) THEN
                   IF (.NOT.trajsmall .OR. l.LE.trajsmalln) THEN
                      IF (cprint%twritebintrajectory) THEN
                         IF (paral%io_parent)&
                              WRITE(4)&
                              iteropt%nfi,&
                              taup(1,j,k),taup(2,j,k),taup(3,j,k),&
                              velp(1,j,k),velp(2,j,k),velp(3,j,k),&
                              fion(1,j,k),fion(2,j,k),fion(3,j,k)
                      ELSE
                         IF (paral%io_parent)&
                              WRITE(4,'(I7,9(2X,F22.14))')&
                              iteropt%nfi,&
                              taup(1,j,k),taup(2,j,k),taup(3,j,k),&
                              velp(1,j,k),velp(2,j,k),velp(3,j,k),&
                              fion(1,j,k),fion(2,j,k),fion(3,j,k)
                      ENDIF
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ELSE
          l=0
          DO k=1,ions1%nsp
             DO j=1,ions0%na(k)
                l=l+1
                IF (l.GE.cprint%minwriteatom.AND.l.LE.cprint%maxwriteatom) THEN
                   IF (.NOT.trajsmall .OR. l.LE.trajsmalln) THEN
                      IF (cprint%twritebintrajectory) THEN
                         IF (paral%io_parent)&
                              WRITE(4)&
                              iteropt%nfi,&
                              taup(1,j,k),taup(2,j,k),taup(3,j,k),&
                              velp(1,j,k),velp(2,j,k),velp(3,j,k),&
                              fion(1,j,k),fion(2,j,k),fion(3,j,k)
                      ELSE
                         IF (paral%io_parent)&
                              WRITE(4,'(I7,9(2X,F22.14))')&
                              iteropt%nfi,&
                              taup(1,j,k),taup(2,j,k),taup(3,j,k),&
                              velp(1,j,k),velp(2,j,k),velp(3,j,k),&
                              fion(1,j,k),fion(2,j,k),fion(3,j,k)
                      ENDIF
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDIF
       IF (paral%io_parent)&
            CALL fileclose(4)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE printp2
  ! ==================================================================


  ! Purpose: print MTS low or high level forces depending on job input
  !      The effective forces which can be a combination of low and high levels
  !      are already printed to the FTRAJECTORY file.
  !      The file is printed in xyz format. (N atoms, title line, xyz data)
  !
  ! Author: Pablo Baudin
  ! Date: May 2018
  subroutine print_mts_forces(forces, title, job)

     implicit none
     !> forces to be printed
     real(real_8), intent(in) :: forces(:,:,:)
     !> title string to be printed
     character(len=100), intent(in) :: title
     !> Job type 'HIGH' or 'LOW' level
     character(len=*), intent(in) :: job

     real(real_8), allocatable :: gr_forces(:,:)
     integer, allocatable :: gr_iat(:)
     integer :: i, ia, is, k, funit, ierr
     logical :: ferror
     character(len=20) :: fname
     character(*), parameter :: procedureN = 'print_mts_forces'

     if (paral%io_parent) then
        ! get file name
        select case(job)
        case('LOW')
           fname = 'MTS_LOW_FORCES.xyz '
        case('HIGH')
           fname = 'MTS_HIGH_FORCES.xyz'
        end select


        ! allocate gromos arrays
        ALLOCATE(gr_iat(ions1%nat+1),STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'allocation problem: gr_iat',&
           __LINE__,__FILE__)
        ALLOCATE(gr_forces(3,ions1%nat),STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'allocation problem: gr_forces',&
           __LINE__,__FILE__)

        ! get forces into gromos ordering
        i=0
        do is=1,ions1%nsp
           do ia=1,ions0%na(is)
              i=i+1
              gr_iat(nat_grm(i))=ions0%iatyp(is)
              do k=1,3
                 gr_forces(k,nat_grm(i))=forces(k,ia,is)
              end do
           end do
        end do


        ! WRITE FILE:
        ! -----------
        !
        funit = 4
        call fileopen(funit,fname,fo_app,ferror)
        !
        ! write number of atoms
        write(funit,*) ions1%nat

        ! write title line
        write(funit,*) title

        ! write ionic forces
        do i=1,ions1%nat
           write(funit,'(a2,3(2x,f22.14))') elem%el(gr_iat(i)),(gr_forces(k,i),k=1,3)
        enddo
        call fileclose(funit)


        ! deallocate gromos arrays
        DEALLOCATE(gr_forces,STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem: gr_forces',&
           __LINE__,__FILE__)
        DEALLOCATE(gr_iat,STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem: gr_iat',&
           __LINE__,__FILE__)
     end if

  end subroutine print_mts_forces

END MODULE printp_utils
