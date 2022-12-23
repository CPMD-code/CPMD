MODULE pm_mdpt_utils
  USE atwf,                            ONLY: atwp
  USE bsym,                            ONLY: bsfac
  USE cl_init_utils,                   ONLY: cl_init
  USE clas,                            ONLY: tclas
  USE cnst,                            ONLY: fbohr
  USE coninp_utils,                    ONLY: raddeg
  USE coor,                            ONLY: tau0
  USE cotr,                            ONLY: cotc0,&
                                             cotr007,&
                                             ntrest,&
                                             resval
  USE ddip,                            ONLY: lenbk
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_def,&
                                             fo_nofp,&
                                             fo_old,&
                                             fo_verb,&
                                             fo_vmark
  USE filnmod,                         ONLY: filbod,&
                                             filn
  USE fnlalloc_utils,                  ONLY: fnlalloc,&
                                             fnldealloc
  USE fusion_utils,                    ONLY: fusion,&
                                             separate
  USE jrotation_utils,                 ONLY: set_orbdist
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE linres,                          ONLY: &
       clrwf, nlinr, nlinw, nolr, nua, nub, td01, td03, urot
  USE lr_in_utils,                     ONLY: lr_in,&
                                             tddft_input
  USE md_driver,                       ONLY: mddiag
  USE mdclas_utils,                    ONLY: mdclas
  USE mdfile_utils,                    ONLY: mdfile
  USE mdmain_utils,                    ONLY: mdmain
  USE mdshop_bo_utils,                 ONLY: mdshop_bo
  USE mdshop_cp_utils,                 ONLY: mdshop_cp
  USE meta_multiple_walkers_utils,     ONLY: mw_filename
  USE mfep,                            ONLY: alpha_pmin,&
                                             factor,&
                                             irest,&
                                             lcvsp,&
                                             mfepi,&
                                             mfepl,&
                                             tolds
  USE mm_dim_utils,                    ONLY: mm_dim
  USE mm_dimmod,                       ONLY: mm_go_qm,&
                                             mm_revert
  USE mm_mddiag_utils,                 ONLY: mm_mddiag
  USE mm_mdmain_utils,                 ONLY: mm_mdmain
  USE mm_mdshop_bo_utils,              ONLY: mm_mdshop_bo
  USE mm_mdshop_cp_utils,              ONLY: mm_mdshop_cp
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum,&
                                             mp_sync
  USE parac,                           ONLY: parai,&
                                             paral
  USE pimd,                            ONLY: &
       grandparent, ipcurr, np_high, np_low, pimd1, pimd3, supergroup, &
       supersource, trep
  USE printpmod,                       ONLY: if1,&
                                             if2,&
                                             if3
  USE prmem_utils,                     ONLY: prmem
  USE pslo,                            ONLY: pslo_com
  USE readsr_utils,                    ONLY: xstring
  USE rreadf_utils,                    ONLY: rreadf
  USE soft,                            ONLY: soft_com
  USE store_types,                     ONLY: restart1
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             maxsys,&
                                             ncpw,&
                                             nkpt
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: nxxfun
  USE vdwcmod,                         ONLY: vdwl
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: pm_mdpt

CONTAINS

  ! ==================================================================
  SUBROUTINE pm_mdpt
    ! ==--------------------------------------------------------------==

    ! Variables
    CHARACTER(*), PARAMETER                  :: procedureN = 'pm_mdpt'
    CHARACTER(len=20), PARAMETER :: f1 = 'LATEST_STRING', f2 = 'STRING.', &
      f3 = 'CONSTRAINT_', f4 = 'METRIC_', f5 = 'PATHMIN'
    INTEGER, PARAMETER                       :: iif1 = fo_verb, &
                                                iif2 = fo_vmark, iif3 = fo_def

    CHARACTER(len=1)                         :: r
    CHARACTER(len=20)                        :: filename
    CHARACTER(len=80)                        :: cflbod, cipnum, line
    COMPLEX(real_8), ALLOCATABLE             :: c0(:,:,:), c2(:,:,:), cm(:), &
                                                sc0(:)
    COMPLEX(real_8), ALLOCATABLE, TARGET     :: c1(:,:,:,:)
    INTEGER :: i, i1, i2, idummy, ierr, ilmax, ilmin, iloop, ip, istep, isub, &
      iunit, j, jcv, k, maxr, minr, mm, ms, n1, n2, nact, nc0, nc0_1d, &
      nc0_2d, nc0_3d, nc1, nc2, nc2_1d, nc2_2d, nc2_3d, ncm, ncol, nmin, &
      nsamples, nsc0, nstate, nstx, nus, nxx, s1, s2, sigma
    LOGICAL                                  :: convge, ferror, statusdummy, &
                                                tmdebug
    REAL(real_8)                             :: cvalf, de, demax, denorm, &
                                                dist, dx, dxmax, dxnorm, &
                                                esum, fact, fvalf, length, &
                                                norm_
    REAL(real_8), ALLOCATABLE :: fbym(:,:), fcol(:), fcolav(:,:), fed(:), &
      fedp(:), gamx(:), gamy(:), l(:), metric(:,:), metricav(:,:,:), no(:,:), &
      pij(:,:,:), resvip(:,:), resvipd(:), resvipo(:,:), resvips(:,:), &
      resvout(:,:), s(:), vpp(:)

    CALL tiset(procedureN,isub)

    tmdebug=.FALSE.
    IF (cntl%bohr) factor=factor*fbohr ! for compatibility with previous version
    ! ==--------------------------------------------------------------==
    ! ==  MEMORY ALLOCATION                                           ==
    ! ==--------------------------------------------------------------==
    ! ==  C0         WF POSITIONS                                     ==
    ! ==  CM         WF VELOCITIES                                    ==
    ! ==  C2         WF FORCES                                        ==
    ! ==  SC0        S*C0 (only Vanderbilt)                           ==
    ! ==  GAMX       CONSTRAINT MATRIX POSITIONS                      ==
    ! ==  GAMY       CONSTRAINT MATRIX VELOCITIES                     ==
    ! ==--------------------------------------------------------------==
    IF ( paral%qmnode ) THEN
       IF (cntl%tddft) THEN
          CALL lr_in
          CALL tddft_input
          cntl%tsymrho=.FALSE.
       ENDIF
       nstate=crge%n
       ! McB
       ! rv        IF ( cntl%tshop ) NSTATE=NST_S0+NST_S1
       ! McB
       nc1=1
       nua=1
       nub=1
       nlinw=0
       nlinr=0
       nus=1
       IF (cntl%tdiag) THEN
          IF (cntl%tddft) CALL stopgm('PM_MDPT',&
               'CNTL%TDDFT.AND.TDIAG NOT POSSIBLE',& 
               __LINE__,__FILE__)
          IF (cntl%tdavi) THEN
             IF (cnti%ndavv.EQ.0) THEN
                cnti%ndavv=2*nstate
             ELSE
                cnti%ndavv=MIN(cnti%ndavv+nstate,2*nstate)
             ENDIF
             nstate=cnti%ndavv
          ENDIF
          ! nc0=2*nkpt%ngwk*nstate*nkpt%nkpnt+8
          nc0_1d = 2 * nkpt%ngwk
          nc0_2d = nstate
          nc0_3d = nkpt%nkpnt
          nc0=nc0_1d*nc0_2d*nc0_3d
          IF (cntl%tfrsblk)THEN
             ncm=2*nkpt%ngwk*MAX(nstate,(cnti%nkry_max+1)*cnti%nkry_block)+8
          ELSEIF (cntl%diis) THEN
             ncm=(nkpt%ngwk*nstate+8)*cnti%mdiis*nkpt%nkpnt+&
                  ((nkpt%ngwk*nstate+8)*cnti%mdiis*nkpt%nkpnt)/4
          ELSE
             ncm=2*nkpt%ngwk*MAX(nstate,cnti%nkry_max*cnti%nkry_block)+8
          ENDIF
          nc2_1d = nkpt%ngwk
          nc2_2d = MAX(atwp%numaormax,nstate)
          nc2_3d = 1
          nc2 = nc2_1d * nc2_2d * nc2_3d
          !nc2=2*nkpt%ngwk*MAX(atwp%numaormax,nstate)+8
          IF (cntl%tdavi)ncm=ncm*2
          ALLOCATE(vpp(nkpt%ngwk),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ELSEIF (cntl%tmdbo.OR.cntl%tmdfile) THEN
          IF (cntl%tsde) THEN
             ncm=8
          ELSE IF (cntl%diis) THEN
             ncm=(nkpt%ngwk*nstate+8)*cnti%mdiis*nkpt%nkpnt+&
                  ((nkpt%ngwk*nstate+8)*cnti%mdiis*nkpt%nkpnt)/4+100
          ELSE IF (cntl%pcg) THEN
             ncm=2*nkpt%ngwk*nstate*nkpt%nkpnt+8
          ENDIF
          ! CB: Need memory for two wf in BS case          
          ! nc0=(2*nkpt%ngwk*nstate*nkpt%nkpnt+8)*bsfac
          nc0_1d = nkpt%ngwk
          nc0_2d = nstate * bsfac
          nc0_3d = nkpt%nkpnt
          nc0=nc0_1d*nc0_2d*nc0_3d
          nc2_1d = nkpt%ngwk
          nc2_2d = nstate
          nc2_3d = bsfac
          nc2 = nc2_1d * nc2_2d* nc2_3d
          !nc2=(2*nkpt%ngwk*nstate+8)*bsfac
          ALLOCATE(vpp(nkpt%ngwk),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          IF (cntl%tddft) THEN
             DO i=1,nstate
                IF (crge%f(i,1).LT.1.e-6_real_8) CALL stopgm('TDDFT',&
                     'ALL STATES HAVE TO BE OCCUPIED',& 
                     __LINE__,__FILE__)
             ENDDO
             ms=MAX(td01%ns_mix,td01%ns_sin,td01%ns_tri)
             IF (td03%tda.AND.td01%ldiag.EQ.1) THEN
                mm=1
             ELSE
                mm=2
             ENDIF
             IF (td03%tda.AND.td01%msubta.GT.0) THEN
                nact=td01%msubta
             ELSE
                nact=nstate
             ENDIF
             nolr=nact
             nc1=2*ncpw%ngw*nact*ms*mm
             nua=nstate
             nub=nact
             nlinw=ms
             nlinr=0
             nus=2*nua*nub
          ENDIF
       ELSE
          IF (cntl%tddft) CALL stopgm('PM_MDPT',&
               'TDDFT.AND.CP NOT POSSIBLE',& 
               __LINE__,__FILE__)
          ! NN  For the BS-CPMD, HS & BS wavefunctions are propagated at the 
          ! same time, thus the array sizes of C0,CM and C2 are multiplied by 2 
          ! nc0=(2*nkpt%ngwk*nstate*nkpt%nkpnt+8)*bsfac
          nc0_1d = nkpt%ngwk
          nc0_2d = nstate * bsfac
          nc0_3d = nkpt%nkpnt
          nc0=nc0_1d*nc0_2d*nc0_3d
          ncm=(2*nkpt%ngwk*nstate*nkpt%nkpnt+8)*bsfac
          nc2_1d = nkpt%ngwk
          nc2_2d = nstate
          nc2_3d = bsfac
          nc2 = nc2_1d * nc2_2d* nc2_3d
          !nc2=(2*nkpt%ngwk*nstate+8)*bsfac
       ENDIF
       nsc0=1
       IF (pslo_com%tivan.OR.cntl%nonort.OR.(cntl%quenchb.AND.cntl%pcg.AND.cntl%pcgmin)) THEN
          nsc0=2*nkpt%ngwk*nstate*nkpt%nkpnt*bsfac
       ENDIF
       IF (cntl%tdipd.OR.vdwl%vdwd) THEN
          IF (cntl%tddft) CALL stopgm('PM_MDPT',&
               'TDDFT AND DIPOLE DYNAMICS NOT POSSIBLE',& 
               __LINE__,__FILE__)
          lenbk=nxxfun(nstate)
          nxx=MAX(2*lenbk*parai%nproc,nc2)
          nsc0=MAX(nxx,nsc0)
          nc2_3d = FLOOR(REAL(nxx,real_8)/REAL(nc2_1d*nc2_2d,real_8))
          nc2 = nc2_1d * nc2_2d* nc2_3d
          !nc2=nxx
       ENDIF
       IF (cntl%tshop) THEN
          nc0=2*nc0
          nc0_2d = nc0_2d * 2
          ncm=2*ncm
          nc2_3d = nc2_3d * 2
          nc2 = nc2_1d * nc2_2d* nc2_3d
          !nc2=2*nc2
       ENDIF
       ALLOCATE(c0(nc0_1d,nc0_2d,nc0_3d),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(c0)
       ALLOCATE(cm(ncm),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(cm)

       IF (nc1 > 1) THEN
          ALLOCATE(c1(ncpw%ngw,nact,ms,mm),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)! (NGW,NOLR,NLINW,2)
       ELSE
          ALLOCATE(c1(1,1,1,1),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
       CALL zeroing(c1)
       ALLOCATE(c2(nc2_1d,nc2_2d,nc2_3d),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(c2)
       ALLOCATE(sc0(nsc0),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(urot(nua,nub,2),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)! TODO check dimensions
       CALL zeroing(urot)!,nus)
       clrwf => c1
       IF (cntl%tdmal) THEN
          CALL set_orbdist(nstate,cnti%nstblk,parai%nproc,nstx)
          ALLOCATE(gamx(nstate*nstx*bsfac),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(gamy(nstate*nstx*bsfac),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ELSE
          ALLOCATE(gamx(nstate*nstate*bsfac),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(gamy(nstate*nstate*bsfac),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
       ! ==--------------------------------------------------------------==
       ! FNL and DFNL
       CALL mm_dim(mm_go_qm,statusdummy)
       CALL fnlalloc(nstate,.TRUE.,cntl%tpres)
       CALL mm_dim(mm_revert,statusdummy)
    ELSE
       nmin = 10
       ALLOCATE(c0(nmin,1,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(c0)
       ALLOCATE(cm(nmin),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(cm)
       ALLOCATE(c1(nmin,1,1,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(c1)
       ALLOCATE(c2(nmin,1,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(c2)
       ALLOCATE(sc0(nmin),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(urot(nmin,1,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)! TODO check dimensions
       ALLOCATE(gamx(nmin),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(gamy(nmin),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(vpp(nmin),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent) CALL prmem('      MDPT')
    ! ==--------------------------------------------------------------==

    IF ( cntl%tshop.AND.(cntl%tfusi.OR.cntl%tsep) ) THEN
       CALL stopgm('MDPT','CNTL%TSHOP.AND.(CNTL%TFUSI.OR.TSEP) NOT POSSIBLE',& 
            __LINE__,__FILE__)
       IF (cntl%tfusi) THEN
          CALL fusion(c0,cm)
       ELSEIF (cntl%tsep) THEN
          CALL separate(c0,cm)
       ENDIF
       go to 500
    ENDIF

    IF ( cntl%tddft.AND.cntl%tshop ) THEN
       CALL stopgm('MDPT','CNTL%TDDFT.AND.TSHOP NOT POSSIBLE',& 
            __LINE__,__FILE__)
    ENDIF


    ! ==--------------------------------------------------------------==
    ! Initialisation of the PATH
    IF (pimd1%tinit) THEN
       ! Generate replica coordinates
       ALLOCATE(trep(3,maxsys%nax,maxsys%nsx,pimd3%np_total),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       IF (pimd1%rrep) THEN
          CALL stopgm('PM_WF','GENERATE REPLICA NOT IMPLEMENTED',& 
               __LINE__,__FILE__)
          ! Classical test 
          IF (pimd1%testpi) THEN
             CALL dcopy(3*maxsys%nax*maxsys%nsx,tau0,1,trep(1,1,1,1),1)
          ENDIF
       ELSEIF (pimd1%repread) THEN
          CALL rreadf
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    ! fpathsave=fo_info%fpath
    ! iapathsave=fo_info%iapath
    ! iepathsave=fo_info%iepath
    ! ==--------------------------------------------------------------==
    ALLOCATE(resvipo(cotr007%mrestr,pimd3%np_total),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(resvout(cotr007%mrestr,pimd3%np_total),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    IF (grandparent) THEN
       DO i=1,cotr007%mrestr
          DO j=1,pimd3%np_total
             resvipo(i,j)=resval(i)
          ENDDO
       ENDDO
    ENDIF
    ALLOCATE(fcol(cotr007%mrestr),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(fcolav(cotr007%mrestr,pimd3%np_total),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(fbym(cotr007%mrestr,pimd3%np_total),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(metric(cotr007%mrestr,cotr007%mrestr),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(metricav(cotr007%mrestr,cotr007%mrestr,pimd3%np_total),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    ALLOCATE(resvip(cotr007%mrestr,pimd3%np_total),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(resvips(cotr007%mrestr,pimd3%np_total),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(resvipd(cotr007%mrestr),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    IF (mfepl%projout) THEN
       ALLOCATE(no(cotr007%mrestr,pimd3%np_total),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(pij(cotr007%mrestr,cotr007%mrestr,pimd3%np_total),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF

    ALLOCATE(l(pimd3%np_total),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(s(pimd3%np_total),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    ALLOCATE(fed(pimd3%np_total),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(fedp(pimd3%np_total),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    IF (grandparent) THEN
       IF (mfepl%rlate_string) THEN
          ! Read LATEST_STRING file
          CALL xstring(f1,s1,s2)
          CALL fileopen(97,f1,fo_old+fo_nofp,ferror)
          IF (ferror) GOTO 100
          READ(unit=97,END=100,err=100,fmt=*) mfepi%nloopold
          CALL fileclose(97)
       ENDIF
       ! Read latest STRING.num file
       CALL mw_filename(f2,filename,mfepi%nloopold)
       CALL xstring(filename,s1,s2)
       CALL fileopen(95,filename,fo_old+fo_nofp,ferror)
       IF (ferror) GOTO 100
       DO i=1,pimd3%np_total
          READ(95,*) idummy,(resvipo(irest(j),i),j=1,mfepi%ncvsp)
          DO j=1,mfepi%ncvsp
             jcv=irest(j)
             IF (ntrest(1,jcv).EQ.1 .OR.&
                  ntrest(1,jcv).EQ.4 .OR.&
                  ntrest(1,jcv).EQ.7 .OR.&
                  ntrest(1,jcv).EQ.12) THEN
                IF (.NOT.cntl%bohr) resvipo(jcv,i)=resvipo(jcv,i)*fbohr
             ELSEIF(ntrest(1,jcv).EQ.2 .OR.&
                  ntrest(1,jcv).EQ.3 .OR.&
                  ntrest(1,jcv).EQ.5) THEN
                CALL raddeg(resvipo(jcv,i),-1)
             ENDIF
          ENDDO
       ENDDO
       CALL fileclose(95)
    ENDIF
    CALL mp_bcast(mfepi%nloopold,supersource,supergroup)
    ! ==--------------------------------------------------------------==
    ilmin=mfepi%nloopold+1
    ilmax=mfepi%nloopold+mfepi%nloop
    DO iloop=ilmin,ilmax
       IF (paral%io_parent) THEN
          WRITE(6,'(/,1X,64("*"))')
          WRITE(6,'(" *",62X,"*")')
          WRITE(6,'(" *",23X,A,I5,24X,"*")') ' LOOP    :',iloop
          WRITE(6,'(" *",62X,"*")')
          WRITE(6,'(1X,64("*"),/)')
       ENDIF
       mfepi%istring=iloop
       convge=.FALSE.
       !IF (paral%parent) THEN
       !CALL mp_bcast(resvipo,SIZE(resvipo),supersource,parentgroup)
       CALL mp_bcast(resvipo,SIZE(resvipo),supersource,supergroup)
       IF (paral%io_parent) THEN
          WRITE(6,'(T7,A,T12,A,T35,A)') 'IP','IREST','RESV'
          DO i=1,pimd3%np_total
             DO j=1,mfepi%ncvsp
                jcv=irest(j)
                WRITE(6,'(2(1X,I7),F22.14)') i,jcv,resvipo(jcv,i)
             ENDDO
          ENDDO
       ENDIF
       !ENDIF

       CALL zeroing(metricav)!,cotr007%mrestr*cotr007%mrestr*pimd3%np_total)
       CALL zeroing(fcolav)!,cotr007%mrestr*pimd3%np_total)
       ! ==--------------------------------------------------------------==
       DO ip=np_low,np_high
          ! IPCUUR is in the include file pimd.inc.
          ipcurr=ip
          IF (pimd1%tinit) THEN
             CALL dcopy(3*maxsys%nax*maxsys%nsx,trep(1,1,1,ip),1,tau0,1)
          ENDIF
          IF (paral%parent) CALL dcopy(cotr007%mrestr,resvipo(1,ip),1,resval(1),1)
          ! Construct filenames
          cflbod='RESTART_'
          WRITE(cipnum,'(I4)') ip
          CALL xstring(cflbod,n1,n2)
          CALL xstring(cipnum,i1,i2)
          filbod=cflbod(n1:n2)//cipnum(i1:i2)//'.1'
          IF (restart1%rlate.OR.iloop.GT.1) THEN
             filn=filbod
          ELSE
             filn=cflbod(n1:n2)//cipnum(i1:i2)
          ENDIF
          ! Construct filepath
          !       cflbod=fpathsave(iapathsave:iepathsave-1)&
          !            //'_'//cipnum(i1:i2)//'/'
          !       write(fo_info%fpath,'(A)') cflbod
          !       call xstring(fo_info%fpath,fo_info%iapath,fo_info%iepath)

          IF (paral%io_parent) THEN
             WRITE(6,'(/,1X,64("*"))')
             WRITE(6,'(" *",62X,"*")')
             WRITE(6,'(" *",23X,A,I5,24X,"*")') ' REPLICA :',iP
             WRITE(6,'(" *",62X,"*")')
             WRITE(6,'(1X,64("*"),/)')
          ENDIF

          if1=iif1
          if2=iif2
          if3=iif3

          IF (tclas) THEN
             IF (tkpts%tkpnt)CALL stopgm('PM_MDPT','K-POINTS NOT IMPLEMENTED',& 
                  __LINE__,__FILE__)
             IF (cntl%tddft)CALL stopgm('PM_MDPT','TDDFT NOT POSSIBLE',& 
                  __LINE__,__FILE__)
             iunit=5
             CALL cl_init(iunit)
             CALL mdclas(c0,cm,c2,sc0,gamx,gamy)
          ELSEIF (cntl%tmdbo) THEN
             IF ( cntl%tqmmm ) THEN
                IF ( cntl%tshop ) THEN
                   CALL stopgm('PM_MDPT',&
                        'CNTL%TMDBO.AND.CNTL%TQMMM.AND.TSHOP NOT POSSIBLE',& 
                        __LINE__,__FILE__)
                   CALL mm_mdshop_bo(c0,cm,c1,c2,sc0,vpp,gamx,gamy)
                ELSE
                   CALL mm_mddiag(c0,cm,c1,c2,sc0,vpp,gamx,gamy)
                ENDIF
             ELSE
                IF ( cntl%tshop ) THEN
                   CALL stopgm('PM_MDPT','CNTL%TMDBO.AND.TSHOP NOT POSSIBLE',& 
                        __LINE__,__FILE__)
                   CALL mdshop_bo(c0,cm,c1,c2,sc0,vpp,gamx,gamy)
                ELSE
                   CALL mddiag(c0,cm,c1,c2,sc0,vpp,gamx,gamy)
                ENDIF
             ENDIF
          ELSEIF (cntl%tmdfile) THEN
             CALL mdfile(c0,cm,c1,c2,sc0,vpp,gamx,gamy)
          ELSE! CP dynamics
             IF (tkpts%tkpnt)CALL stopgm('PM_MDPT','K-POINTS NOT IMPLEMENTED',& 
                  __LINE__,__FILE__)
             IF (cntl%tddft)CALL stopgm('PM_MDPT','TDDFT.AND.CP NOT POSSIBLE',& 
                  __LINE__,__FILE__)
             IF ( cntl%tqmmm ) THEN
                IF ( cntl%tshop ) THEN
                   CALL stopgm('PM_MDPT','CNTL%TQMMM.AND.TSHOP NOT POSSIBLE',& 
                        __LINE__,__FILE__)
                   CALL mm_mdshop_cp(c0,cm,c2,sc0,gamx,gamy)
                ELSE
                   CALL mm_mdmain(c0,cm,c2,sc0,gamx,gamy)
                ENDIF
             ELSE
                IF ( cntl%tshop ) THEN
                   CALL stopgm('PM_MDPT','SURFACE HOPPING NOT POSSIBLE',& 
                        __LINE__,__FILE__)
                   CALL mdshop_cp(c0,cm,c2,sc0,gamx,gamy)
                ELSE
                   CALL mdmain(c0,cm,c2,sc0,gamx,gamy)
                ENDIF
             ENDIF
          ENDIF

          IF (soft_com%exnomore) THEN
             soft_com%exsoft=.FALSE.
             soft_com%exnomore=.FALSE.
          ELSE IF (soft_com%exsoft) THEN
             GOTO 500
          ENDIF

          IF (paral%io_parent) THEN
             ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! averages in each directory
             ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! Read CONSTRAINT file
             CALL mw_filename(f3,cflbod,ipcurr)
             CALL xstring(cflbod,s1,s2)
             cflbod=cflbod(s1:s2)//'.'
             CALL mw_filename(cflbod,filename,iloop)
             CALL fileopen(31,filename,fo_old,ferror)
             IF (ferror) GOTO 100
             DO i=1,mfepi%nequi
                DO j=1,cotc0%mcnstr
                   READ(unit=31,END=200,err=200,fmt=*) line
                ENDDO
                DO j=1,cotr007%mrestr
                   READ(unit=31,END=200,err=200,fmt=*) line
                ENDDO
             ENDDO
             nsamples=0
             DO WHILE(.TRUE.)
                DO j=1,cotc0%mcnstr
                   READ(unit=31,END=20,err=200,fmt=*)&
                        istep,idummy,cvalf,fvalf
                ENDDO
                READ(unit=31,END=20,err=200,fmt=*)&
                     istep,idummy,cvalf,fvalf,fcol(1),r
                IF (lcvsp(1)) fcolav(1,ip)=fcolav(1,ip)+fcol(1)
                DO j=2,cotr007%mrestr
                   READ(unit=31,END=200,err=200,fmt=*)&
                        istep,idummy,cvalf,fvalf,fcol(j),r
                   IF (lcvsp(j)) fcolav(j,ip)=fcolav(j,ip)+fcol(j)
                ENDDO
                nsamples=nsamples+1
             ENDDO
20           CONTINUE
             CALL fileclose(31)
             ! Read METRIC file
             CALL mw_filename(f4,cflbod,ipcurr)
             CALL xstring(cflbod,s1,s2)
             cflbod=cflbod(s1:s2)//'.'
             CALL mw_filename(cflbod,filename,iloop)
             CALL fileopen(90,filename,fo_old,ferror)
             IF (ferror) GOTO 100
             DO i=1,mfepi%nequi
                DO j=1,cotr007%mrestr
                   READ(unit=90,END=200,err=200,fmt=*) line
                ENDDO
             ENDDO
             DO i=1,nsamples
                DO j=1,cotr007%mrestr
                   READ(unit=90,END=200,err=200,fmt=*)&
                        (metric(j,k),k=1,cotr007%mrestr)
                   DO k=1,cotr007%mrestr
                      IF (lcvsp(j).AND.lcvsp(k))&
                           metricav(j,k,ip)=metricav(j,k,ip)+metric(j,k)
                   ENDDO
                ENDDO
             ENDDO
             CALL fileclose(90)

             IF (nsamples.EQ.0) CALL stopgm('PM_MDPT',&
                  'FAIL IN SAMPLING MEAN FORCE',& 
                  __LINE__,__FILE__)
             fact=1.0_real_8/REAL(nsamples,kind=real_8)
             CALL dscal(cotr007%mrestr,-1.0_real_8*fact,fcolav(1,ip),1)
             CALL dscal(cotr007%mrestr*cotr007%mrestr,fact,metricav(1,1,ip),1)

             WRITE(6,'(20X,A)') ' <<<<< AVERAGED GRADIENTS >>>>>'
             WRITE(6,'(3(1PE20.10))')&
                  (fcolav(irest(j),ip),j=1,mfepi%ncvsp)
             WRITE(6,'(/,22X,A)') ' <<<<< AVERAGED METRIC >>>>>'
             DO k=1,mfepi%ncvsp
                WRITE(6,'(3(1PE20.10))')&
                     (metricav(irest(j),irest(k),ip),j=1,mfepi%ncvsp)
             ENDDO
          ENDIF

          IF (soft_com%exsoft) GOTO 500

       ENDDO
       ! ==--------------------------------------------------------------==

       CALL mp_sync(supergroup)

       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! collect averages
       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       CALL mp_sum(fcolav,cotr007%mrestr*pimd3%np_total,supergroup)
       CALL mp_sum(metricav,cotr007%mrestr*cotr007%mrestr*pimd3%np_total,supergroup)

       IF (grandparent) THEN
          ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! grandparent updates string
          ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          ! multiply force and metric
          WRITE(line,'(I4)') cotr007%mrestr
          CALL xstring(line,i1,i2)
          DO i=1,pimd3%np_total
             CALL dgemm('N','N',cotr007%mrestr,1,cotr007%mrestr,1.0_real_8,&
                  metricav(1,1,i),cotr007%mrestr,&
                  fcolav(1,i),cotr007%mrestr,0.0_real_8,fbym(1,i),cotr007%mrestr)
             IF (tmdebug) THEN
                WRITE(6,'(I7,'//line(i1:i2)//'(2X,F22.14))')&
                     i,(fcolav(j,i),j=1,cotr007%mrestr)
             ENDIF
          ENDDO

          ! estimate free energy difference along the path
          CALL fileopen(98,f5,fo_def+fo_nofp,ferror)
          WRITE(6,'(/,18X,A)') ' <<<<< FREE ENERGY DIFFERENCE >>>>>'
          WRITE(6,'(T7,A,T21,A,T32,A,T43,A,T56,A,T67,A)')&
               'IP','FED','RESV','GRADV','RESV','GRADV'
          fed(1)=0.0_real_8
          ncol=NINT(REAL(mfepi%ncvsp,kind=real_8)*0.5_real_8)
          minr=1
          maxr=MIN(minr+1,mfepi%ncvsp)
          WRITE(6,'(1X,I7,F15.9,T24,4F12.6)')&
               1,fed(1),(resvipo(irest(j),1),fcolav(irest(j),1),j=minr,maxr)
          DO k=2,ncol
             minr=maxr+1
             maxr=MIN(minr+1,mfepi%ncvsp)
             WRITE(6,'(T24,4F12.6)')&
                  (resvipo(irest(j),1),fcolav(irest(j),1),j=minr,maxr)
          ENDDO
          WRITE(line,'(I4)') 2*mfepi%ncvsp+1
          CALL xstring(line,i1,i2)
          WRITE(98,'(I7,'//line(i1:i2)//'F15.9)')&
               1,fed(1),(resvipo(irest(j),1),j=1,mfepi%ncvsp),&
               (fcolav(irest(j),1),j=1,mfepi%ncvsp)
          DO i=2,pimd3%np_total
             esum=0.0_real_8
             DO j=1,mfepi%ncvsp
                jcv=irest(j)
                esum=esum+0.5_real_8*(resvipo(jcv,i)-resvipo(jcv,i-1))&
                     *(fcolav(jcv,i)+fcolav(jcv,i-1))
             ENDDO
             fed(i)=fed(i-1)+esum
             minr=1
             maxr=MIN(minr+1,mfepi%ncvsp)
             WRITE(6,'(1X,I7,F15.9,T24,4F12.6)')&
                  i,fed(i),(resvipo(irest(j),i),&
                  fcolav(irest(j),i),j=minr,maxr)
             DO k=2,ncol
                minr=maxr+1
                maxr=MIN(minr+1,mfepi%ncvsp)
                WRITE(6,'(T24,4F12.6)')&
                     (resvipo(irest(j),i),fcolav(irest(j),i),j=minr,maxr)
             ENDDO
             WRITE(98,'(I7,'//line(i1:i2)//'F15.9)')&
                  i,fed(i),(resvipo(irest(j),i),j=1,mfepi%ncvsp),&
                  (fcolav(irest(j),i),j=1,mfepi%ncvsp)
          ENDDO
          CALL fileclose(98)

          IF (iloop.GT.ilmin) THEN
             demax=0.0_real_8
             denorm=0.0_real_8
             DO i=1,pimd3%np_total
                de=ABS(fed(i)-fedp(i))
                denorm=denorm+de*de
                demax=MAX(de,demax)
             ENDDO
             denorm=denorm/REAL(pimd3%np_total,kind=real_8)
             WRITE(6,'(/,1X,A,/,3X,A,1PE14.3,2X,A,1PE14.3/)')&
                  ' DEVIATION OF FREE ENERGY DIFFERENCES: ',&
                  ' DEMAX = ',demax,' DENORM = ',denorM
          ENDIF
          CALL dcopy(pimd3%np_total,fed,1,fedp,1)

          ! evolution of collective variables

          IF (mfepl%projout) THEN
             ! with projecting out of forces
             DO j=1,cotr007%mrestr
                no(j,1)=resvipo(j,2)-resvipo(j,1)
                no(j,pimd3%np_total)=resvipo(j,pimd3%np_total-1)-resvipo(j,pimd3%np_total)
             ENDDO
             DO i=2,pimd3%np_total-1
                DO j=1,cotr007%mrestr
                   no(j,i)=resvipo(j,i+1)-resvipo(j,i)
                ENDDO
                CALL iscal(no(1,i),fbym(1,i),cotr007%mrestr,sigma)
                DO j=1,cotr007%mrestr
                   no(j,i)=resvipo(j,i+sigma)-resvipo(j,i)
                ENDDO
             ENDDO
             DO i=1,pimd3%np_total
                CALL norme(no(1,i),cotr007%mrestr,norm_)
                CALL dscal(cotr007%mrestr,1.0_real_8/norm_,no(1,i),1)
             ENDDO
             DO i=1,pimd3%np_total
                DO j=1,cotr007%mrestr
                   DO k=1,cotr007%mrestr
                      pij(j,k,i)=-1.0_real_8*no(j,i)*no(k,i)
                   ENDDO
                   pij(j,j,i)=1.0_real_8+pij(j,j,i)
                ENDDO
                CALL dgemm('N','N',cotr007%mrestr,1,cotr007%mrestr,1.0_real_8,&
                     pij(1,1,i),cotr007%mrestr,&
                     fbym(1,i),cotr007%mrestr,0.0_real_8,fcolav(1,i),cotr007%mrestr)
             ENDDO

             DO i=1,cotr007%mrestr
                DO j=2,pimd3%np_total-1
                   resvip(i,j)=resvipo(i,j)-factor*fcolav(i,j)
                ENDDO
                resvip(i,1)=resvipo(i,1)
                resvip(i,pimd3%np_total)=resvipo(i,pimd3%np_total)
             ENDDO

          ELSE
             ! without projecting out forces
             DO i=1,cotr007%mrestr
                DO j=2,pimd3%np_total-1
                   resvip(i,j)=resvipo(i,j)-factor*fbym(i,j)
                ENDDO
                resvip(i,1)=resvipo(i,1)
                resvip(i,pimd3%np_total)=resvipo(i,pimd3%np_total)
             ENDDO
          ENDIF

          ! smoothing
          DO i=1,cotr007%mrestr
             DO j=2,pimd3%np_total-1
                resvips(i,j)=alpha_pmin*resvip(i,j-1)&
                     +(1.0_real_8-2.0_real_8*alpha_pmin)*resvip(i,j)&
                     +alpha_pmin*resvip(i,j+1)
             ENDDO
             DO j=2,pimd3%np_total-1
                resvip(i,j)=resvips(i,j)
             ENDDO
          ENDDO

          ! reparametrization
          CALL zeroing(l)!,pimd3%np_total)
          DO i=2,pimd3%np_total
             CALL distcurv(resvip(1,i),resvip(1,i-1),cotr007%mrestr,dist)
             l(i)=l(i-1)+dist
          ENDDO
          length=l(pimd3%np_total)

          DO i=1,pimd3%np_total
             s(i)=(REAL(i,kind=real_8)-1.0_real_8)*length/REAL(pimd3%np_total-1,kind=real_8)
          ENDDO

          DO j=1,cotr007%mrestr     ! first point and last point
             resvips(j,1)=resvip(j,1)
             resvips(j,pimd3%np_total)=resvip(j,pimd3%np_total)
          ENDDO

          DO i=2,pimd3%np_total-1
             DO j=2,pimd3%np_total
                IF (l(j).GE.s(i)) THEN
                   CALL distcurv(resvip(1,j-1),resvip(1,j),&
                        cotr007%mrestr,dist)
                   DO k=1,cotr007%mrestr
                      resvipd(k)=resvip(k,j)-resvip(k,j-1)
                      resvips(k,i)=(s(i)-l(j-1))*resvipd(k)/dist&
                           +resvip(k,j-1)
                   ENDDO
                   GOTO 10
                ENDIF
             ENDDO
10           CONTINUE
          ENDDO

          ! check convergence on string positions
          dxmax=0.0_real_8
          dxnorm=0.0_real_8
          DO i=1,pimd3%np_total
             DO j=1,mfepi%ncvsp
                jcv=irest(j)
                dx=ABS(resvips(jcv,i)-resvipo(jcv,i))
                dxnorm=dxnorm+dx*dx
                dxmax=MAX(dx,dxmax)
             ENDDO
          ENDDO
          dxnorm=dxnorm/REAL(pimd3%np_total*mfepi%ncvsp,kind=real_8)
          WRITE(6,'(/,1X,A,/,3X,A,1PE14.3,2X,A,1PE14.3/)')&
               ' DEVIATION OF STRING POSITIONS: ',&
               ' DXMAX = ',dxmax,' DXNORM = ',dxnorM
          IF (dxmax.LT.tolds) THEN
             convge=.TRUE.
             GOTO 30
          ENDIF

          CALL dcopy(cotr007%mrestr*pimd3%np_total,resvips,1,resvipo,1)
          CALL dcopy(cotr007%mrestr*pimd3%np_total,resvips,1,resvout,1)

          ! write in STRING.num and LATEST_STRING for new run and trajectory
          CALL mw_filename(f2,filename,iloop)
          CALL xstring(filename,s1,s2)
          CALL fileopen(96,filename,fo_def+fo_nofp,ferror)
          WRITE(line,'(I4)') mfepi%ncvsp
          CALL xstring(line,i1,i2)
          DO i=1,pimd3%np_total
             DO j=1,mfepi%ncvsp
                jcv=irest(j)
                IF (ntrest(1,jcv).EQ.1 .OR.&
                     ntrest(1,jcv).EQ.4 .OR.&
                     ntrest(1,jcv).EQ.7 .OR.&
                     ntrest(1,jcv).EQ.12) THEN
                   IF (.NOT.cntl%bohr) resvout(jcv,i)=resvout(jcv,i)/fbohr
                ELSEIF (ntrest(1,jcv).EQ.2 .OR.&
                     ntrest(1,jcv).EQ.3 .OR.&
                     ntrest(1,jcv).EQ.5) THEN
                   CALL raddeg(resvout(jcv,i),1)
                ENDIF
             ENDDO
             WRITE(96,'(I7,'//line(i1:i2)//'(2X,F22.14))')&
                  i,(resvout(irest(j),i),j=1,mfepi%ncvsp)
          ENDDO
          CALL fileclose(96)
          CALL fileopen(97,f1,fo_def+fo_nofp,ferror)
          WRITE(97,*) iloop
          CALL fileclose(97)
       ENDIF
30     CONTINUE
       CALL mp_sync(supergroup)
       CALL mp_bcast(convge,supersource,supergroup)
       IF (convge) GOTO 500
    ENDDO
    ! ==--------------------------------------------------------------==

500 CONTINUE
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    DEALLOCATE(fedp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(fed,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(s,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(l,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (mfepl%projout) THEN
       DEALLOCATE(pij,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(no,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    DEALLOCATE(resvipd,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(resvips,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(resvip,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(metricav,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(metric,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(fbym,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(fcolav,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(fcol,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(resvipo,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF ( paral%qmnode ) THEN
       CALL fnldealloc(.TRUE.,cntl%tpres)
    ELSE
       DEALLOCATE(vpp,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    DEALLOCATE(gamy,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(gamx,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(urot,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(sc0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(c2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(c1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(cm,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(c0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    GOTO 9999
    ! ==--------------------------------------------------------------==
100 CONTINUE
    CALL stopgm('PM_MDPT','REQUIRED FILE '//filename(s1:s2)//&
         ' NOT PRESENT OR NOT READABLE',& 
         __LINE__,__FILE__)
200 CONTINUE
    CALL stopgm('PM_MDPT','ERROR WHILE READING '//filename(s1:s2)//&
         ' FILE',& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
9999 CONTINUE
    CALL tihalt(procedureN,isub)
    RETURN
  END SUBROUTINE pm_mdpt
  ! ==================================================================
  SUBROUTINE distcurv(c1,c2,ncolvar,dist)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ncolvar
    REAL(real_8)                             :: c2(ncolvar), c1(ncolvar), dist

    INTEGER                                  :: i
    REAL(real_8)                             :: c12

    dist=0.0_real_8
    DO i=1,ncolvar
       c12=c1(i)-c2(i)! access only once c1(i) and c2(i)
       dist=dist+c12*c12
    ENDDO
    dist=SQRT(dist)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE distcurv
  ! ==================================================================
  SUBROUTINE iscal(c1,c2,ncolvar,sigma)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ncolvar
    REAL(real_8)                             :: c2(ncolvar), c1(ncolvar)
    INTEGER                                  :: sigma

    INTEGER                                  :: i
    REAL(real_8)                             :: scal

    scal=0.0_real_8
    DO i=1,ncolvar
       scal=scal+c1(i)*c2(i)
    ENDDO
    IF (scal.GT.0) THEN
       sigma=1
    ELSE
       sigma=-1
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE iscal
  ! ==================================================================
  SUBROUTINE norme(c1,ncolvar,norm)
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ncolvar
    REAL(real_8)                             :: c1(ncolvar), norm

    INTEGER                                  :: i

    norm=0.0_real_8
    DO i=1,ncolvar
       norm=norm+c1(i)*c1(i)
    ENDDO
    norm=SQRT(norm)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE norme
END MODULE pm_mdpt_utils
