MODULE pm_gmopts_utils
  USE atwf,                            ONLY: atwp
  USE coor,                            ONLY: tau0,&
                                             velp
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE filnmod,                         ONLY: filbod,&
                                             filn
  USE fnlalloc_utils,                  ONLY: fnlalloc,&
                                             fnldealloc
  USE kinds,                           ONLY: real_8
  USE linres,                          ONLY: &
       clrwf, nlinr, nlinw, nolr, nua, nub, td01, td03, urot
  USE lr_in_utils,                     ONLY: lr_in,&
                                             tddft_input
  USE parac,                           ONLY: paral
  USE pimd,                            ONLY: ipcurr,&
                                             np_high,&
                                             np_low,&
                                             pimd1,&
                                             pimd3,&
                                             trep
  USE prmem_utils,                     ONLY: prmem
  USE pslo,                            ONLY: pslo_com
  USE readsr_utils,                    ONLY: xstring
  USE rreadf_utils,                    ONLY: rreadf
  USE store_types,                     ONLY: restart1
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             maxsys,&
                                             ncpw,&
                                             nkpt
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: pm_gmopts

CONTAINS

  ! ==================================================================
  SUBROUTINE pm_gmopts
    ! ==--------------------------------------------------------------==
    CHARACTER(*), PARAMETER                  :: procedureN = 'pm_gmopts'

    CHARACTER(len=12)                        :: cflbod, cipnum
    COMPLEX(real_8), ALLOCATABLE             :: c0(:), c2(:)
    COMPLEX(real_8), ALLOCATABLE, TARGET     :: c1(:,:,:,:)
    INTEGER                                  :: i, i1, i2, ierr, ip, isub, &
                                                mm, ms, n1, n2, nact, nc0, &
                                                nc1, nc2, ncm, neigv, ngde, &
                                                npme, nsc0, nstate, nus, nx
    REAL(real_8), ALLOCATABLE                :: cm(:), eigv(:), gde(:), &
                                                pme(:), sc0(:), vpp(:)

! TODO refactor this 
! ==--------------------------------------------------------------==
! ==  MEMORY ALLOCATION                                           ==
! ==--------------------------------------------------------------==
! ==                      WAVEFUNCTIONS                           ==
! ==  C0    : WAVEFUNCTIONS                                       ==
! ==  C2    : GRADIENTS                                           ==
! ==  CM    : VELOCITIES  OR BACKUP WAVEFUNCTION                  ==
! ==  SC0   : S**(-1/2)*C0 (nonorthogonal orbitals)               ==
! ==  PME   : cntl%diis WF             HNM1 (for CG)              ==
! ==  GDE   : cntl%diis GRADIENTS                                 ==
! ==  VPP   : WF HESSIAN                                          ==
! ==  EIGV  : ORBITAL ENERGIES                                    ==
! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    IF (cntl%tddft) THEN
       CALL lr_in
       CALL tddft_input
       cntl%tsymrho=.FALSE.
    ENDIF
    nstate=crge%n
    IF (cntl%tdiag.AND.cntl%tdavi) THEN
       IF (cnti%ndavv.EQ.0) THEN
          cnti%ndavv=2*nstate
       ELSE
          cnti%ndavv=MIN(cnti%ndavv+nstate,2*nstate)
       ENDIF
       nstate=cnti%ndavv
    ENDIF
    nc0=2*nkpt%ngwk*nstate*nkpt%nkpnt+8
    neigv=nstate*nkpt%nkpts
    nc2=2*nkpt%ngwk*MAX(atwp%numaormax,nstate)+8
    IF (cntl%tdiag) THEN
       IF (cntl%tfrsblk) THEN
          ncm=2*nkpt%ngwk*MAX(nstate,(cnti%nkry_max+1)*cnti%nkry_block)+8
       ELSEIF (cntl%diis) THEN
          ncm=(nkpt%ngwk*nstate+8)*cnti%mdiis*nkpt%nkpnt+&
               ((nkpt%ngwk*nstate+8)*cnti%mdiis*nkpt%nkpnt)/4
          ncm = ncm*2
       ELSE
          ncm=2*nkpt%ngwk*MAX(nstate,cnti%nkry_max*cnti%nkry_block)+8
       ENDIF
    ELSE
       ncm = 1
       IF (cntl%simul) THEN
          ncm=nc0
       ENDIF
    ENDIF
    IF (pslo_com%tivan.OR.cntl%nonort.OR.(cntl%pcg.AND.cntl%pcgmin)) THEN
       nsc0=MAX(2*nkpt%ngwk*nstate*nkpt%nkpnt,nc2)
    ELSEIF (cntl%tdiag .AND.cntl%diis) THEN
       nsc0=nc2
    ELSE
       nsc0=1
    ENDIF
    IF (cntl%tdebfor) THEN
       ncm = nc0
    ENDIF
    IF (cntl%tddft) THEN
       DO i=1,nstate
          IF (crge%f(i,1).LT.1.e-6_real_8) THEN
             CALL stopgm("cntl%tddft","ALL STATES HAVE TO BE OCCUPIED",& 
                  __LINE__,__FILE__)
          ENDIF
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
    ELSE
       nc1=1
       nua=1
       nub=1
       nlinw=0
       nlinr=0
       nus=1
    ENDIF
    ALLOCATE(c0(nc0),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    IF (nc1 > 1) THEN
       ALLOCATE(c1(ncpw%ngw,nact,ms,mm),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)! (NGW,NOLR,NLINW,2)
    ELSE
       ALLOCATE(c1(1,1,1,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    clrwf => c1

    ALLOCATE(c2(nc2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(cm(ncm),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(eigv(neigv),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(sc0(nsc0),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(vpp(nkpt%ngwk),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(urot(nua,nub,2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__) ! TODO check dimensions
    CALL zeroing(urot)!,nus)
    IF (cntl%tsde) THEN
       npme = 1
       ngde = 1
    ELSE IF (cntl%diis .AND. .NOT. cntl%tdiag) THEN
       npme = (nkpt%ngwk*nstate+8)*cnti%mdiis*nkpt%nkpnt
       ngde = ((nkpt%ngwk*nstate+8)*cnti%mdiis*nkpt%nkpnt)/4
    ELSE IF (cntl%pcg) THEN
       npme = 2*nkpt%ngwk*nstate*nkpt%nkpnt
       ngde = 1
    ELSE IF (cntl%tdiag) THEN
       IF (cntl%tdavi) THEN
          npme = nc0
          ngde = nc0
       ELSE
          npme = 1
          ngde = 1
       ENDIF
    ELSE
       IF (paral%io_parent) THEN
          WRITE(6,*) ' WRONG OPTION FOR WAVEFUNCTION OPTIMIZATION'
       ENDIF
       CALL stopgm('PM_GMOPTS',' ',& 
            __LINE__,__FILE__)
    ENDIF
    ALLOCATE(pme(npme),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(gde(ngde),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    nx=nstate
    IF (cntl%tdavi) THEN
       nx=cnti%ndavv
    ENDIF
    ! FNL and DFNL
    CALL fnlalloc(nstate,.TRUE.,cntl%tpres)
    IF (paral%io_parent) THEN
       CALL prmem(' PM_GMOPTS')
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
    DO ip=np_low,np_high
       ! IPCUUR is in the include file pimd.inc.
       ipcurr=ip
       IF (pimd1%tinit) THEN
          CALL dcopy(3*maxsys%nax*maxsys%nsx,trep(1,1,1,ipcurr),1,tau0,1)
       ENDIF
       ! Construct filenames
       cflbod='RESTART_'
       IF (paral%io_parent)&
            WRITE(cipnum,'(I4)') ipcurr
       CALL xstring(cflbod,n1,n2)
       CALL xstring(cipnum,i1,i2)
       filbod=cflbod(n1:n2)//cipnum(i1:i2)//'.1'
       IF (restart1%rlate) THEN
          filn=filbod
       ELSE
          filn=cflbod(n1:n2)//cipnum(i1:i2)
       ENDIF
       IF (paral%io_parent) THEN
          WRITE(6,'(/,1X,64("*"))')
          WRITE(6,'(" *",62X,"*")')
          WRITE(6,'(" *",23X,A,I5,24X,"*")') ' REPLICA :',ipcurr
          WRITE(6,'(" *",62X,"*")')
          WRITE(6,'(1X,64("*"),/)')
       ENDIF
       CALL rgmopt(c0,c1,c2,cm,sc0,pme,gde,vpp,eigv)
    ENDDO
    ! ==--------------------------------------------------------------==
    ! Deallocations
    IF (cntl%tmemchk) THEN
    ENDIF
    DEALLOCATE(c0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(c1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(c2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(cm,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(eigv,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(urot,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    CALL fnldealloc(.TRUE.,cntl%tpres)
    DEALLOCATE(sc0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(vpp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(pme,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(gde,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(velp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    RETURN
  END SUBROUTINE pm_gmopts
  ! ==================================================================

END MODULE pm_gmopts_utils
