MODULE recpnew_utils
  USE adat,                            ONLY: atwt,&
                                             defrag,&
                                             elem
  USE atom,                            ONLY: atom_common,&
                                             ecpfiles,&
                                             rps,&
                                             rv,&
                                             rw,&
                                             vr
  USE cp_xc_utils,                     ONLY: cp_xc_functional_env
  USE dpot,                            ONLY: dpot_mod
  USE error_handling,                  ONLY: stopgm
  USE fitpack_utils,                   ONLY: curv1,&
                                             curv2
  USE func,                            ONLY: func1,&
                                             func2
  USE inscan_utils,                    ONLY: inscan
  USE ions,                            ONLY: al,&
                                             bl,&
                                             ions0,&
                                             ions1,&
                                             rcl
  USE kinds,                           ONLY: real_8
  USE machine,                         ONLY: m_getarg,&
                                             m_getenv,&
                                             m_iargc
  USE nlcc,                            ONLY: corecg,&
                                             corei,&
                                             corel,&
                                             corer,&
                                             rcgrid
  USE nlps,                            ONLY: nghcom,&
                                             nghtol,&
                                             nlps_com
  USE parac,                           ONLY: paral
  USE pslo,                            ONLY: pslo_com
  USE ragg,                            ONLY: raggio
  USE readsr_utils,                    ONLY: readsi,&
                                             readsr,&
                                             xstring
  USE readvan_utils,                   ONLY: readvan
  USE rmas,                            ONLY: rmass
  USE sgpp,                            ONLY: mcfun,&
                                             mpro,&
                                             sgpp1,&
                                             sgpp2
  USE special_functions,               ONLY: cp_erf
  USE system,                          ONLY: cntl,&
                                             lx,&
                                             maxsys,&
                                             nbrx
  USE vdbp,                            ONLY: &
       betar, dion, ncpr1, qfunc, qqq, qrl, r, rab, rsatom, rscore, ru, &
       rucore, vdb_pawf, vdb_r
  USE vdbt,                            ONLY: itmax,&
                                             vdbti

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: recpnew
  PUBLIC :: allelec
  PUBLIC :: recpeam
  PUBLIC :: tgrid
  PUBLIC :: ckgrid
  PUBLIC :: get_pplib

CONTAINS

  ! ==================================================================
  SUBROUTINE recpnew(isp,ecpnam)
    ! ==--------------------------------------------------------------==
    ! ==  Reads Pseudopotential Input (New Format)                    ==
    ! ==--------------------------------------------------------------==
    ! == INPUT:                                                       ==
    ! == ISP     Species index                                        ==
    ! == ECPNAM  Filename of pseudopotential file                     ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: isp
    CHARACTER(len=*)                         :: ecpnam

    CHARACTER(*), PARAMETER                  :: procedureN = 'recpnew'
    INTEGER, PARAMETER                       :: iunit = 21 

    CHARACTER(len=120)                       :: ecplib
    CHARACTER(len=200)                       :: fnames
    CHARACTER(len=80)                        :: line
    INTEGER                                  :: i, ia, ie, ierr, ii, il, &
                                                inghcom, iostat, iout, ir, &
                                                it, iv, j, k, l, lenecp, lp, &
                                                lread, m, meshv
    INTEGER, SAVE                            :: ivan1 = 0
    LOGICAL                                  :: erread, exists, tsca
    REAL(real_8)                             :: dal, dxc, excf, &
                                                hmat(mpro*mpro), s1, xctest
    REAL(real_8), ALLOCATABLE                :: c(:), fc(:), rr(:), temp(:)

! ==--------------------------------------------------------------==
! Allocation of local variables

    ALLOCATE(rr(maxsys%mmaxx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(c(maxsys%mmaxx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(fc(maxsys%mmaxx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(temp(maxsys%mmaxx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    sgpp1%tsgp(isp)=.FALSE.
    pslo_com%tnum(isp)=.FALSE.
    pslo_com%tlog(isp)=.FALSE.
    tsca=.FALSE.
    CALL get_pplib(ecplib,lenecp)
    CALL xstring(ecpnam,ia,ie)
    fnames=ecplib(1:lenecp)//ecpnam(ia:ie)
    IF (paral%io_parent)&
         INQUIRE(file=fnames,exist=exists)
    IF (.NOT.exists) THEN
       fnames=ecpnam(ia:ie)
       IF (paral%io_parent)&
            INQUIRE(file=fnames,exist=exists)
       IF (.NOT.exists) THEN
          IF (paral%io_parent)&
               WRITE(6,*) ' RECPNEW| ECPFILE NOT FOUND ',fnames
          CALL stopgm('RECPNEW',' ',& 
               __LINE__,__FILE__)
       ENDIF
    ENDIF
    ecpfiles(isp)=fnames
    IF (pslo_com%tvan(isp)) THEN
       ! ..Vanderbilt PP in binary or formatted form
       IF (ivan1.EQ.0) THEN
          ivan1=1
          ALLOCATE(rscore(maxsys%mmaxx,maxsys%nsx),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(dion(nbrx,nbrx,maxsys%nsx),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(betar(maxsys%mmaxx,nbrx,maxsys%nsx),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(qqq(nbrx,nbrx,maxsys%nsx),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(qfunc(maxsys%mmaxx,nbrx,nbrx,maxsys%nsx),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(qrl(maxsys%mmaxx,nbrx,nbrx,lx,maxsys%nsx),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(r(maxsys%mmaxx,maxsys%nsx),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(rucore(maxsys%mmaxx,nbrx,maxsys%nsx),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(ru(maxsys%mmaxx,maxsys%nsx),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(rab(maxsys%mmaxx,maxsys%nsx),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(rsatom(maxsys%mmaxx,maxsys%nsx),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(vdb_pawf(maxsys%nsx,maxsys%mmaxx,nbrx),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(vdb_r(maxsys%mmaxx,maxsys%nsx),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
       pslo_com%tnum(isp)=.TRUE.
       dpot_mod%tkb(isp)=.FALSE.
       ions0%igau(isp)=0
       ! 
       CALL readvan(isp,fnames)
       ! 
       IF (ncpr1%ifpcor(isp).EQ.1) corel%tnlcc(isp) = .TRUE.
       IF (rmass%pma0(isp).EQ.0.0_real_8)rmass%pma0(isp)=atwt(ions0%iatyp(isp))
       raggio(isp)=defrag(ions0%iatyp(isp))
       RETURN
    ENDIF
    IF (paral%io_parent)&
         OPEN(unit=iunit,file=fnames,status='OLD')
    IF (paral%io_parent)&
         REWIND(iunit)
    ! ==--------------------------------------------------------------==
    ! ..General info about PP and atom
    ierr=inscan(iunit,'&ATOM')
    IF (ierr.NE.0) THEN
       IF (paral%io_parent)&
            WRITE(6,*) ' RECPNEW| &ATOM SECTION NOT FOUND '
       CALL stopgm('RECPNEW',' ',& 
            __LINE__,__FILE__)
    ENDIF
    IF (paral%io_parent)&
         READ(iunit,END=20,err=20,fmt='(A)') line
    ia=INDEX(line,'=')+1
    CALL readsi(line,ia,iout,ions0%iatyp(isp),erread)
    IF (paral%io_parent)&
         READ(iunit,END=20,err=20,fmt='(A)') line
    ia=INDEX(line,'=')+1
    CALL readsr(line,ia,iout,ions0%zv(isp),erread)
    IF (paral%io_parent)&
         READ(iunit,END=20,err=20,fmt='(A)') line
    ia=INDEX(line,'=')+1
    CALL readsr(line,ia,iout,excf,erread)
    ia=iout+1
    CALL readsr(line,ia,iout,s1,erread)
    dal=ABS(s1-func2%salpha)
    IF (paral%io_parent)&
         READ(iunit,END=20,err=20,fmt='(A)') line
    IF (rmass%pma0(isp).EQ.0.0_real_8)rmass%pma0(isp)=atwt(ions0%iatyp(isp))
    raggio(isp)=defrag(ions0%iatyp(isp))
    excf=excf-100._real_8*(INT(excf)/100)
    IF (cntl%use_xc_driver) THEN
        xctest = cp_xc_functional_env%get_pseudo_code()
    ELSE
        xctest = 10._real_8*func1%mgcx+func1%mgcc
    END IF
    dxc=ABS(xctest-excf)
    IF (dxc.GT.0.1_real_8.OR.dal.GT.1.e-4_real_8) THEN
       CALL xstring(ecpnam,ia,ie)
       IF (paral%io_parent)&
            WRITE(6,'(/,1X,64("!"))')
       IF (paral%io_parent)&
            WRITE(6,'(A,A)') ' WARNING! XC FUNCTIONALS INCONSISTENT FOR ',&
            ecpnam(ia:ie)
       IF (paral%io_parent)&
            WRITE(6,'(1X,64("!"))')
    ENDIF
    ! ..The info section
    ierr=inscan(iunit,'&INFO')
    it=0
    DO i=1,60
       IF (paral%io_parent)&
            READ(iunit,END=20,err=20,fmt='(A)') line
       IF (INDEX(line,'&END').NE.0) GOTO 10
       it=it+1
       vdbti(it,isp)=line(1:66)
    ENDDO
10  CONTINUE
    itmax(isp)=it
    ! ..Potential section
    ierr=inscan(iunit,'&POTENTIAL')
    IF (ierr.NE.0) THEN
       IF (paral%io_parent)&
            WRITE(6,*) ' RECPNEW| &POTENTIAL SECTION NOT FOUND '
       CALL stopgm('RECPNEW',' ',& 
            __LINE__,__FILE__)
    ENDIF
    ! ..Normconserving PP
    IF (paral%io_parent)&
         READ(iunit,END=20,err=20,fmt='(A)') line
    IF (INDEX(line,'BHS').NE.0) THEN
       ! ..Bachelet type
       ions0%igau(isp)=3
       pslo_com%tnum(isp)=.FALSE.
       IF (paral%io_parent)&
            READ(iunit,*) ions0%wrc(isp,1),ions0%rc(isp,1),ions0%wrc(isp,2),ions0%rc(isp,2)
       DO i=1,dpot_mod%lmax(isp)
          DO j=1,ions0%igau(isp)
             IF (paral%io_parent)&
                  READ(iunit,*) rcl(j,isp,i),al(j,isp,i),bl(j,isp,i)
          ENDDO
       ENDDO
    ELSEIF (INDEX(line,'CAR').NE.0) THEN
       ! ..VonBarth-Car type
       ions0%igau(isp)=1
       pslo_com%tnum(isp)=.FALSE.
       IF (paral%io_parent)&
            READ(iunit,*) ions0%rc(isp,1)
       DO i=1,dpot_mod%lmax(isp)
          IF (paral%io_parent)&
               READ(iunit,*) rcl(1,isp,i),al(1,isp,i),bl(1,isp,i)
       ENDDO
    ELSEIF (INDEX(line,'GOEDECK').NE.0) THEN
       ! ..Stefan Goedecker PP
       sgpp1%tsgp(isp)=.TRUE.
       pslo_com%tnum(isp)=.FALSE.
       dpot_mod%tkb(isp)=.FALSE.
       ions0%igau(isp)=0
       IF (paral%io_parent)&
            READ(iunit,*) dpot_mod%lmax(isp)
       dpot_mod%lmax(isp)=dpot_mod%lmax(isp)+1
       dpot_mod%lloc(isp)=dpot_mod%lmax(isp)
       dpot_mod%lskip(isp)=dpot_mod%lmax(isp)+5
       IF (paral%io_parent)&
            READ(iunit,*) sgpp2%rcsg(isp)
       IF (paral%io_parent)&
            READ(iunit,*) sgpp1%nclsg(isp),(sgpp2%clsg(i,isp),i=1,sgpp1%nclsg(isp))
       IF (sgpp1%nclsg(isp).GT.mcfun) CALL stopgm('RECPNEW','C_FUNCTION',& 
            __LINE__,__FILE__)
       iv=0
       lp=0
       DO l=1,dpot_mod%lmax(isp)-1
          IF (paral%io_parent)&
               READ(iunit,*) sgpp2%rcnl(l,isp),sgpp2%npro(l,isp),(hmat(i),i=1,&
               (sgpp2%npro(l,isp)*(sgpp2%npro(l,isp)+1))/2)
          IF (sgpp2%npro(l,isp).GT.mpro) CALL stopgm('RECPNEW','HMAT',& 
               __LINE__,__FILE__)
          ii=0
          DO i=1,sgpp2%npro(l,isp)
             DO j=i,sgpp2%npro(l,isp)
                ii=ii+1
                sgpp2%hlsg(i,j,l,isp)=hmat(ii)
                sgpp2%hlsg(j,i,l,isp)=hmat(ii)
             ENDDO
          ENDDO
          DO m=1,2*l-1
             lp=lp+1
             DO k=1,sgpp2%npro(l,isp)
                iv=iv+1
                nghtol(iv,isp)=l-1
                sgpp2%lpval(iv,isp)=lp
                sgpp2%lfval(iv,isp)=k
             ENDDO
          ENDDO
       ENDDO
       nlps_com%ngh(isp)=iv
    ELSE
       ! ..Numeric
       ! Logscale option.
       IF (INDEX(line,'LOGSCALE').NE.0.OR.&
            (INDEX(line,'ECHELLE').NE.0&
            .AND.INDEX(line,'LOG').NE.0)) THEN
          ! Use logscale defined in the pp file (do not interpolate)
          tsca=.TRUE.
          IF (paral%io_parent)&
               READ(iunit,END=20,err=20,fmt='(A)') line
       ENDIF
       pslo_com%tnum(isp)=.TRUE.
       ions0%igau(isp)=0
       ! Read number of points
       CALL readsi(line,1,iout,meshv,erread)
       IF (meshv.GT.maxsys%mmaxx) THEN
          IF (paral%io_parent)&
               WRITE(6,*) ' RECPNEW! Potential SECTION'
          CALL xstring(ecpnam,ia,ie)
          IF (paral%io_parent)&
               WRITE(6,*) ' RECPNEW! MESH FOR ',ecpnam(ia:ie),' IS ',&
               meshv
          IF (paral%io_parent)&
               WRITE(6,*) ' RECPNEW! MAX NUMBER OF SPLINE POINTS:',maxsys%mmaxx
          IF (paral%io_parent)&
               WRITE(6,*) ' RECPNEW! INCREASE SPLINE POINTS NUMBER'
          CALL stopgm('RECPNEW','MESH TOO BIG',& 
               __LINE__,__FILE__)
       ENDIF
       IF (paral%io_parent) THEN
          DO j=1,meshv
             READ(iunit,*,iostat=iostat) rr(j),(vr(j,isp,il),il=1,dpot_mod%lmax(isp))
             IF(iostat/=0) CALL stopgm(procedureN,'wrong PP file format or wrong PP input parameters',& 
                  __LINE__,__FILE__)
          ENDDO
       ENDIF
       CALL ckgrid(rr(1),rr(meshv),rv(1,isp),meshv,atom_common%clogvp(isp))
       IF (.NOT.tsca) THEN
          ! Interpolate data in another grid
          DO il=1,dpot_mod%lmax(isp)
             CALL tgrid(rr,meshv,rv(1,isp),meshv,vr(1,isp,il),&
                  maxsys%mmaxx,c,fc,temp)
          ENDDO
       ELSE
          ! Do not interpolate data. Use the same grid.
          DO j=1, meshv
             rv(j,isp)=rr(j)
          ENDDO
       ENDIF
       atom_common%meshvp(isp)=meshv
    ENDIF
    ! ..Put analytic potentials on a logarithmic mesh
    IF (.NOT.(pslo_com%tnum(isp).OR.sgpp1%tsgp(isp)))&
         CALL agrid(isp,rr)
    ! ..Read Atomic Wavefunction (if Kleinman-Bylander Integration)
    IF (dpot_mod%tkb(isp)) THEN
       lread=dpot_mod%lmax(isp)
       IF (dpot_mod%lloc(isp).EQ.dpot_mod%lmax(isp)) lread=lread-1
       IF (lread.GT.0) THEN
          ierr=inscan(iunit,'&WAVEFUNCTION')
          IF (ierr.NE.0) THEN
             IF (paral%io_parent)&
                  WRITE(6,*) ' RECPNEW| &WAVEFUNCTION SECTION NOT FOUND '
             CALL stopgm('RECPNEW',' ',& 
                  __LINE__,__FILE__)
          ENDIF
          IF (paral%io_parent)&
               READ(iunit,*) meshv
          IF (meshv.GT.maxsys%mmaxx) THEN
             IF (paral%io_parent)&
                  WRITE(6,*) ' RECPNEW! Wavefunction SECTION'
             CALL xstring(ecpnam,ia,ie)
             IF (paral%io_parent)&
                  WRITE(6,*) ' RECPNEW! MESH FOR ',ecpnam(ia:ie),' IS ',&
                  meshv
             IF (paral%io_parent)&
                  WRITE(6,*) ' RECPNEW! MAX NUMBER OF SPLINE POINTS:',maxsys%mmaxx
             IF (paral%io_parent)&
                  WRITE(6,*) ' RECPNEW! INCREASE SPLINE POINTS NUMBER'
             CALL stopgm('RECPNEW','MESH TOO BIG',& 
                  __LINE__,__FILE__)
          ENDIF
          DO j=1,meshv
             IF (paral%io_parent)&
                  READ(iunit,*) rr(j),(rps(j,isp,il),il=1,lread)
          ENDDO
          ! ..We need same grid for wavefunctions and potentials
          IF (tsca) THEN
             ! Check if the mesh is the same for potential and wavefunctions
             IF (meshv.NE.atom_common%meshvp(isp)) THEN
                CALL xstring(fnames,ia,ie)
                IF (paral%io_parent)&
                     WRITE(6,*) ' RECPNEW! YOU USE LOGSCALE OPTION IN THE',&
                     ' PSEUDOPOTENTIAL FILE:',fnames(ia:ie)
                IF (paral%io_parent)&
                     WRITE(6,*) ' RECPNEW! Potential MESH=',atom_common%meshvp(isp)
                IF (paral%io_parent)&
                     WRITE(6,*) ' RECPNEW! Wavefunction MESH=',atom_common%meshw(isp)
                CALL stopgm('RECPNEW','THE LOGSCALES ARE DIFFERENT'//&
                     'FOR Potential AND WAVEFUNCTIONS',& 
                     __LINE__,__FILE__)
             ENDIF
          ENDIF
          atom_common%clogw(isp)=atom_common%clogvp(isp)
          atom_common%meshw(isp)=atom_common%meshvp(isp)
          CALL dcopy(atom_common%meshw(isp),rv(1,isp),1,rw(1,isp),1)
          IF (.NOT.tsca) THEN
             DO il=1,lread
                CALL tgrid(rr,meshv,rw(1,isp),atom_common%meshw(isp),rps(1,isp,il),&
                     maxsys%mmaxx,c,fc,temp)
             ENDDO
          ENDIF
       ELSE
          atom_common%clogw(isp)=atom_common%clogvp(isp)
          atom_common%meshw(isp)=atom_common%meshvp(isp)
       ENDIF
       ! ..aa   To decrease the number of non-local projectors
       iv=0
       DO il=1,dpot_mod%lmax(ions1%nsp)
          IF (il.NE.dpot_mod%lloc(ions1%nsp).AND.il.NE.dpot_mod%lskip(ions1%nsp)) THEN
             l=il-1
             nlps_com%ngh(ions1%nsp)=nlps_com%ngh(ions1%nsp)+(2*l+1)
             IF (l.EQ.0) THEN
                inghcom=0
             ELSE
                inghcom=l*l
             ENDIF
             DO j=1,2*l+1
                iv=iv+1
                nghtol(iv,ions1%nsp)=l
                nghcom(iv,ions1%nsp)=inghcom+j
             ENDDO
          ENDIF
       ENDDO
       ! ..aa 
    ENDIF
    ! ..Nonlinear core correction
    IF (corel%tnlcc(isp)) THEN
       ierr=inscan(iunit,'&NLCC')
       IF (ierr.NE.0) THEN
          WRITE(6,*) ' RECPNEW: &NLCC SECTION NOT FOUND'
          CALL stopgm('RECPNEW',' ',& 
               __LINE__,__FILE__)
       ENDIF
       READ(iunit,'(A)') linE
       IF (INDEX(line,'ANAL').NE.0) THEN
          ! ..Analytic definition of core charges
          corei%nlcct(isp)=1
          READ(iunit,*) corer%anlcc(isp),corer%bnlcc(isp),corer%enlcc(isp)
       ELSEIF (INDEX(line,'NUME').NE.0) THEN
          ! ..Numerical definition of core charges
          corei%nlcct(isp)=2
          READ(iunit,*) meshv
          IF (meshv.GT.maxsys%mmaxx) THEN
             WRITE(6,*) ' RECPNEW! NONLINEAR CORE CORRECTION SECTION'
             CALL xstring(ecpnam,ia,ie)
             WRITE(6,*) ' RECPNEW! MESH FOR ',ECPNAM(IA:IE),' IS ',&
                  meshV
             WRITE(6,*) ' RECPNEW! MAX NUMBER OF SPLINE POINTS:',maxsys%mmaxx
             WRITE(6,*) ' RECPNEW! INCREASE SPLINE POINTS NUMBER'
             CALL stopgm('RECPNEW','MESH TOO BIG',& 
                  __LINE__,__FILE__)
          ENDIF
          DO ir=1,meshv
             READ(iunit,*) rr(ir),corecg(ir,isp)
          ENDDO
          CALL ckgrid(rr(1),rr(meshv),rcgrid(1,isp),meshv,corer%clogcc(isp))
          CALL tgrid(rr,meshv,rcgrid(1,isp),meshv,corecg(1,isp),&
               maxsys%mmaxx,c,fc,temp)
          corei%meshcc(isp)=meshv
       ELSE
          IF (paral%io_parent)&
               WRITE(6,*) ' TYPE OF NLCC DEFINITION NOT CLEAR'
          IF (paral%io_parent)&
               WRITE(6,*) line
          CALL stopgm('RNLCC',' ',& 
               __LINE__,__FILE__)
       ENDIF
    ELSE
       ierr=inscan(iunit,'&NLCC')
       IF (ierr.EQ.0) THEN
          IF (paral%io_parent)&
               WRITE(6,'(/,1X,64("! "))')
          IF (paral%io_parent)&
               WRITE(6,'(1X,A)')&
               'WARNING! PSEUDOPOTENTIAL GENERATED WITH NLCC '
          CALL xstring(ecpnam,ia,ie)
          IF (paral%io_parent)&
               WRITE(6,'(1X,A,A)') 'WARNING! FOR ',ecpnam(ia:ie)
          IF (paral%io_parent)&
               WRITE(6,'(1X,64("! "))')
       ENDIF
    ENDIF
    IF (paral%io_parent)&
         CLOSE(iunit)
    ! ==--------------------------------------------------------------==
    ! Deallocation of local variables
    DEALLOCATE(rr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(c,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(fc,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(temp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
    ! ==--------------------------------------------------------------==
20  CONTINUE
    IF (paral%io_parent)&
         WRITE(6,*) ' ERROR WHILE READING PP DEFINITIONS'
    CALL stopgm('RECPNEW',' ',& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE recpnew
  ! ==================================================================
  SUBROUTINE ckgrid(rp1,rpl,rw,mesh,clog)
    ! ==--------------------------------------------------------------==
    ! == Construct logarithmic grid from RP1 to RPL with MESH points  ==
    ! == in RW(1:MESH).                                               ==
    ! == CLOG is the step                                             ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rp1, rpl
    INTEGER                                  :: mesh
    REAL(real_8)                             :: rw(mesh), clog

    INTEGER                                  :: ir
    REAL(real_8)                             :: cf

    clog=LOG(rpl/rp1)/REAL(mesh-1,kind=real_8)
    cf=EXP(clog)
    rw(1)=rp1
    DO ir=2,mesh
       rw(ir)=rw(ir-1)*cf
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE ckgrid
  ! ==================================================================
  SUBROUTINE tgrid(rp,mp,rw,mw,f,mmaxx,c,fc,temp)
    ! ==--------------------------------------------------------------==
    ! == Interpolate F(1:MP) from RP(1:MP) grid to RW(1:MW)           ==
    ! == Output: F(1:MW)                                              ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: mp
    REAL(real_8)                             :: rp(mp)
    INTEGER                                  :: mw
    REAL(real_8)                             :: rw(mw), f(*)
    INTEGER                                  :: mmaxx
    REAL(real_8)                             :: c(mmaxx), fc(mmaxx), &
                                                temp(mmaxx)

    INTEGER                                  :: i, ierr

    CALL dcopy(mp,f(1),1,fc(1),1)
    CALL curv1(mp,rp,fc,0.0_real_8,0.0_real_8,3,c,temp,0.0_real_8,ierr)
    DO i=1,mw
       f(i)=curv2(rw(i),mp,rp,fc,c,0.0_real_8)
    ENDDO
    RETURN
  END SUBROUTINE tgrid
  ! ==================================================================
  SUBROUTINE agrid(isp,rr)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: isp
    REAL(real_8)                             :: rr(*)

    INTEGER                                  :: ib, il, ir, meshv
    REAL(real_8)                             :: rr1, rr2, rx

    IF (ions0%igau(isp).GT.1) THEN
       ions0%rc(isp,1)=1._real_8/ions0%rc(isp,1)**0.5_real_8
       ions0%rc(isp,2)=1._real_8/ions0%rc(isp,2)**0.5_real_8
       DO il=1,dpot_mod%lmax(isp)
          DO ib=1,ions0%igau(isp)
             rcl(ib,isp,il) = 1._real_8/rcl(ib,isp,il)**0.5_real_8
          ENDDO
       ENDDO
    ENDIF
    atom_common%clogvp(isp)=LOG(1.025_real_8)
    meshv=INT(LOG(7200.0_real_8*REAL(ions0%iatyp(isp),kind=real_8))/atom_common%clogvp(isp))
    IF (meshv.GT.maxsys%mmaxx) THEN
       atom_common%clogvp(isp)=REAL(meshv,kind=real_8)/REAL(maxsys%mmaxx,kind=real_8)*atom_common%clogvp(isp)
       meshv=maxsys%mmaxx
    ENDIF
    rr1=0.00625_real_8/REAL(ions0%iatyp(isp),kind=real_8)
    rx=LOG(rr1)+meshv*atom_common%clogvp(isp)
    rr2=EXP(rx)
    CALL ckgrid(rr1,rr2,rv(1,isp),meshv,atom_common%clogvp(isp))
    IF (ions0%igau(isp).EQ.1) THEN
       DO il=1,dpot_mod%lmax(isp)
          DO ir=1,meshv
             rx=rv(ir,isp)
             vr(ir,isp,il)=-ions0%zv(isp)/rx*cp_erf(rx/ions0%rc(isp,1))&
                  +(al(1,isp,il)+bl(1,isp,il)*rx*rx)&
                  *exponential(-(rx/rcl(1,isp,il))**2)
          ENDDO
       ENDDO
    ELSEIF (ions0%igau(isp).EQ.3) THEN
       DO il=1,dpot_mod%lmax(isp)
          DO ir=1,meshv
             rx=rv(ir,isp)
             vr(ir,isp,il)=-ions0%zv(isp)/rx*&
                  (ions0%wrc(isp,1)*cp_erf(rx/ions0%rc(isp,1))+&
                  ions0%wrc(isp,2)*cp_erf(rx/ions0%rc(isp,2)))&
                  +(al(1,isp,il)+bl(1,isp,il)*rx*rx)&
                  *exponential(-(rx/rcl(1,isp,il))**2)&
                  +(al(2,isp,il)+bl(2,isp,il)*rx*rx)&
                  *exponential(-(rx/rcl(2,isp,il))**2)&
                  +(al(3,isp,il)+bl(3,isp,il)*rx*rx)&
                  *exponential(-(rx/rcl(3,isp,il))**2)
          ENDDO
       ENDDO
    ELSE
       DO ir=1,meshv
          rx=rv(ir,isp)
          vr(ir,isp,1)=-ions0%zv(isp)/rx
       ENDDO
    ENDIF
    atom_common%meshvp(isp)=meshv
    pslo_com%tnum(isp)=.TRUE.
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE agrid
  ! ==================================================================
  FUNCTION exponential(x)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: x, exponential

    REAL(real_8), PARAMETER                  :: xmin = -400._real_8 

! ==--------------------------------------------------------------==

    IF (x.LT.xmin) THEN
       exponential = 0._real_8
    ELSE
       exponential = EXP(x)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION exponential
  ! ==================================================================
  SUBROUTINE allelec(isp,ecpnam)
    ! ==--------------------------------------------------------------==
    ! ==  GENERATES 1/R Potential                                     ==
    ! ==--------------------------------------------------------------==
    ! == INPUT:                                                       ==
    ! == ISP     Species index                                        ==
    ! == ECPNAM  Atom name                                            ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: isp
    CHARACTER(len=40)                        :: ecpnam

    CHARACTER(*), PARAMETER                  :: procedureN = 'allelec'

    CHARACTER(len=2)                         :: en
    INTEGER                                  :: i, ia, ie, ierr, ifound
    REAL(real_8), ALLOCATABLE                :: rr(:)

    ALLOCATE(rr(maxsys%mmaxx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    sgpp1%tsgp(isp)=.FALSE.
    pslo_com%tvan(isp)=.FALSE.
    pslo_com%tnum(isp)=.TRUE.
    pslo_com%tlog(isp)=.FALSE.
    dpot_mod%tkb(isp)=.FALSE.
    ecpfiles(isp)="NONE"
    dpot_mod%lmax(isp)=1
    ions0%igau(isp)=0
    CALL xstring(ecpnam,ia,ie)
    IF (ia.EQ.ie) THEN
       en=" "//ecpnam(ia:ie)
    ELSEIF (ie.EQ.ia+1) THEN
       en=ecpnam(ia:ie)
    ELSE
       CALL stopgm("ALLELEC","ATOM NAME "//ecpnam(ia:ie),& 
            __LINE__,__FILE__)
    ENDIF
    ifound=0
    DO i=1,99
       IF (en .EQ. elem%el(i)) ifound=i
    ENDDO
    IF (ifound.EQ.0) CALL stopgm("ALLELEC","ATOM NAME "//ecpnam(ia:ie)&
         ,& 
         __LINE__,__FILE__)
    ions0%iatyp(isp)=ifound
    IF (rmass%pma0(isp).EQ.0.0_real_8)rmass%pma0(isp)=atwt(ions0%iatyp(isp))
    raggio(isp)=defrag(ions0%iatyp(isp))
    ions0%zv(isp)=REAL(ifound,kind=real_8)
    vdbti(1,isp)="================================="//&
         "================================="
    vdbti(2,isp)=" "
    vdbti(3,isp)="                      FULL Potential ("//en//")"
    vdbti(4,isp)=" "
    vdbti(5,isp)="================================="//&
         "================================="
    itmax(isp)=5
    CALL agrid(isp,rr)
    DEALLOCATE(rr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE allelec
  ! ==================================================================
  SUBROUTINE recpeam(isp,ecpnam)
    ! ==--------------------------------------------------------------==
    ! ==  Reads Pseudopotential For EAM Atoms                         ==
    ! ==--------------------------------------------------------------==
    ! == INPUT:                                                       ==
    ! == ISP     Species index                                        ==
    ! == ECPNAM  Filename of pseudopotential file                     ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: isp
    CHARACTER(len=40)                        :: ecpnam

    CHARACTER(*), PARAMETER                  :: procedureN = 'recpeam'
    INTEGER, PARAMETER                       :: iunit = 21 

    CHARACTER(len=120)                       :: ecplib
    CHARACTER(len=200)                       :: fnames
    CHARACTER(len=80)                        :: line
    INTEGER                                  :: i, ia, ie, ierr, iout, it, &
                                                lenecp
    LOGICAL                                  :: erread, exists
    REAL(real_8), ALLOCATABLE                :: rr(:)

! ==--------------------------------------------------------------==

    ALLOCATE(rr(maxsys%mmaxx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    sgpp1%tsgp(isp)=.FALSE.
    pslo_com%tvan(isp)=.FALSE.
    pslo_com%tnum(isp)=.FALSE.
    pslo_com%tlog(isp)=.FALSE.
    CALL get_pplib(ecplib,lenecp)
    CALL xstring(ecpnam,ia,ie)
    fnames=ecplib(1:lenecp)//ecpnam(ia:ie)
    IF (paral%io_parent)&
         INQUIRE(file=fnames,exist=exists)
    IF (.NOT.exists) THEN
       fnames=ecpnam(ia:ie)
       IF (paral%io_parent)&
            INQUIRE(file=fnames,exist=exists)
       IF (.NOT.exists) THEN
          IF (paral%io_parent)&
               WRITE(6,*) ' RECPEAM| ECPFILE NOT FOUND ',fnames
          CALL stopgm('RECPEAM',' ',& 
               __LINE__,__FILE__)
       ENDIF
    ENDIF
    ecpfiles(isp)=fnames
    IF (paral%io_parent)&
         OPEN(unit=iunit,file=fnames,status='OLD')
    IF (paral%io_parent)&
         REWIND(iunit)
    ! ==--------------------------------------------------------------==
    ! ..General info about PP and atom
    ierr=inscan(iunit,'&ATOM')
    IF (ierr.NE.0) THEN
       IF (paral%io_parent)&
            WRITE(6,*) ' RECPEAM| &ATOM SECTION NOT FOUND '
       CALL stopgm('RECPEAM',' ',& 
            __LINE__,__FILE__)
    ENDIF
    IF (paral%io_parent)&
         READ(iunit,END=20,err=20,fmt='(A)') line
    ia=INDEX(line,'=')+1
    CALL readsi(line,ia,iout,ions0%iatyp(isp),erread)
    IF (paral%io_parent)&
         READ(iunit,END=20,err=20,fmt='(A)') line
    ia=INDEX(line,'=')+1
    CALL readsr(line,ia,iout,ions0%zv(isp),erread)
    IF (rmass%pma0(isp).EQ.0.0_real_8)rmass%pma0(isp)=atwt(ions0%iatyp(isp))
    raggio(isp)=defrag(ions0%iatyp(isp))
    ! ..The info section
    ierr=inscan(iunit,'&INFO')
    it=0
    DO i=1,60
       IF (paral%io_parent)&
            READ(iunit,END=20,err=20,fmt='(A)') line
       IF (INDEX(line,'&END').NE.0) GOTO 10
       it=it+1
       vdbti(it,isp)=line(1:66)
    ENDDO
10  CONTINUE
    itmax(isp)=it
    ! ..Potential section
    ierr=inscan(iunit,'&Potential')
    IF (ierr.NE.0) THEN
       IF (paral%io_parent)&
            WRITE(6,*) ' RECPEAM| &Potential SECTION NOT FOUND '
       CALL stopgm('RECPEAM',' ',& 
            __LINE__,__FILE__)
    ENDIF
    ! ..VonBarth-Car type
    ions0%igau(isp)=1
    pslo_com%tnum(isp)=.FALSE.
    IF (paral%io_parent)&
         READ(iunit,*) ions0%rc(isp,1)
    DO i=1,dpot_mod%lmax(isp)
       IF (paral%io_parent)&
            READ(iunit,*) rcl(1,isp,i),al(1,isp,i),bl(1,isp,i)
    ENDDO
    CALL agrid(isp,rr)
    IF (ions0%rc(isp,1).LT.1.e-10_real_8) ions0%igau(isp)=-1
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent)&
         CLOSE(iunit)
    DEALLOCATE(rr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    RETURN
    ! ==--------------------------------------------------------------==
20  CONTINUE
    IF (paral%io_parent)&
         WRITE(6,*) ' ERROR WHILE READING PP DEFINITIONS'
    CALL stopgm('RECPEAM',' ',& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE recpeam
  ! ==================================================================
  SUBROUTINE get_pplib(ecplib,lenecp)
    CHARACTER(len=*)                         :: ecplib
    INTEGER                                  :: lenecp

    INTEGER                                  :: icarg

! 
! 

    icarg=m_iargc()
    IF (icarg.GT.1) THEN
       CALL m_getarg(2,ecplib)
       lenecp=INDEX(ecplib,' ')-1
       ecplib=ecplib(1:lenecp)//'/'
       lenecp=lenecp+1
    ELSE
       CALL m_getenv('CPMD_PP_LIBRARY_PATH',ecplib)
       lenecp=INDEX(ecplib,' ')-1
       IF (lenecp.EQ.0) THEN
          CALL m_getenv('PP_LIBRARY_PATH',ecplib)
          lenecp=INDEX(ecplib,' ')-1
          IF (lenecp.EQ.0) THEN
             ecplib='.'
             lenecp=1
          ENDIF
          ecplib=ecplib(1:lenecp)//'/'
          lenecp=lenecp+1
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE get_pplib
  ! ==================================================================
END MODULE recpnew_utils
