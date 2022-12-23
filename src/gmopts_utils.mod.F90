MODULE gmopts_utils
  USE atwf,                            ONLY: atwp
  USE bsym,                            ONLY: bsfac
  USE coor,                            ONLY: velp
  USE ddip,                            ONLY: lenbk
  USE debfor_utils,                    ONLY: debfor
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE fnlalloc_utils,                  ONLY: fnlalloc,&
                                             fnldealloc
  USE kinds,                           ONLY: real_8
  USE linres,                          ONLY: &
       clrwf, nlinr, nlinw, nolr, nua, nub, td01, td03, urot
  USE lr_in_utils,                     ONLY: lr_in,&
                                             tddft_input
  USE parac,                           ONLY: parai,&
                                             paral
  USE prmem_utils,                     ONLY: prmem
  USE pslo,                            ONLY: pslo_com
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             ncpw,&
                                             nkpt
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: nxxfun
  USE vdwcmod,                         ONLY: vdwl
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: gmopts

CONTAINS

  ! ==================================================================
  SUBROUTINE gmopts
    ! ==--------------------------------------------------------------==
    ! Variables

    CHARACTER(*), PARAMETER                  :: procedureN = 'gmopts'

    COMPLEX(real_8), ALLOCATABLE             :: c0(:), c2(:), cm(:), gde(:), &
                                                pme(:), sc0(:)
    COMPLEX(real_8), ALLOCATABLE, TARGET     :: c1(:,:,:,:)
    INTEGER                                  :: i, ierr, isub, mm, ms, nact, &
                                                nc0, nc1, nc2, ncm, neigv, &
                                                ngde, npme, nsc0, nstate, &
                                                nus, nxx
    REAL(real_8), ALLOCATABLE                :: eigv(:), vpp(:)

! TODO refactor CLRWF
! ==--------------------------------------------------------------==
! ==  MEMORY ALLOCATION                                           ==
! ==--------------------------------------------------------------==
! ==                      WAVEFUNCTIONS                           ==
! ==  C0    : WAVEFUNCTIONS                                       ==
! ==  C2    : GRADIENTS                                           ==
! ==  CM    : VELOCITIES  OR BACKUP WAVEFUNCTION                  ==
! ==  SC0   : S**(-1/2)*C0 (nonorthogonal orbitals)               ==
! ==  PME   : DIIS WF             HNM1 (for CG)                   ==
! ==  GDE   : DIIS GRADIENTS                                      ==
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
    nc0=(2*nkpt%ngwk*nstate*nkpt%nkpnt+8)*bsfac
    neigv=nstate*nkpt%nkpts*bsfac
    nc2=(2*nkpt%ngwk*MAX(atwp%numaormax,nstate)+8)*bsfac
    IF (cntl%tdiag) THEN
       IF (cntl%tfrsblk) THEN
          ncm=nkpt%ngwk*MAX(nstate,(cnti%nkry_max+1)*cnti%nkry_block)+8
       ELSEIF (cntl%diis) THEN
          ncm=(nkpt%ngwk*nstate+8)*cnti%mdiis*nkpt%nkpnt+&
               ((nkpt%ngwk*nstate+8)*cnti%mdiis*nkpt%nkpnt)/4
       ELSE
          ncm=nkpt%ngwk*MAX(nstate,cnti%nkry_max*cnti%nkry_block)+8
       ENDIF
    ELSE
       ncm = 1
       IF (cntl%simul) THEN
          ncm=nc0
       ENDIF
    ENDIF
    IF (pslo_com%tivan.OR.cntl%nonort.OR.(cntl%pcg.AND.cntl%pcgmin)) THEN
       nsc0=MAX(nkpt%ngwk*nstate*nkpt%nkpnt,nc2)
    ELSEIF (cntl%tdiag .AND.cntl%diis) THEN
       nsc0=nc2
    ELSE
       nsc0=1
    ENDIF
    IF (cntl%tdipd.OR.vdwl%vdwd) THEN
       lenbk=nxxfun(nstate)
       nxx=MAX(lenbk*parai%nproc,nc2)
       nsc0=MAX(nxx,nsc0)
       nc2=MAX(nxx,nc2)
    ENDIF
    IF (cntl%tdebfor) THEN
       ncm = nc0
    ENDIF
    IF (cntl%tddft) THEN
       DO i=1,nstate
          IF (crge%f(i,1).LT.1.e-6_real_8) THEN
             ! TODO make sure F(I) == F(I,1) here
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

    IF (nc1 /= 1) THEN
       ALLOCATE(c1(ncpw%ngw,nact,ms,mm),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ELSE
       ALLOCATE(c1(1,1,1,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF

    ! TODO refactor CLRWF 
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
    ALLOCATE(sc0(nsc0*bsfac),STAT=ierr)
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
       npme = (nkpt%ngwk*nstate+8)*cnti%mdiis*nkpt%nkpnt /2
       ngde = ((nkpt%ngwk*nstate+8)*cnti%mdiis*nkpt%nkpnt)/8
    ELSE IF (cntl%pcg) THEN
       npme = nkpt%ngwk*nstate*nkpt%nkpnt
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
       IF (paral%parent) THEN
          IF (paral%io_parent)&
               WRITE(6,*) ' WRONG OPTION FOR WAVEFUNCTION OPTIMIZATION'
       ENDIF
       CALL stopgm('GMOPTS',' ',& 
            __LINE__,__FILE__)
    ENDIF
    ALLOCATE(pme(npme),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(gde(ngde),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! FNL and DFNL
    CALL fnlalloc(nstate,.TRUE.,cntl%tpres)
    IF (paral%parent) THEN
       CALL prmem('    GMOPTS')
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (cntl%tdebfor) THEN
       cntl%tsymrho=.FALSE.
       CALL debfor(c0,c1,c2,cm,sc0,pme,gde,vpp,eigv)
    ELSE
       CALL rgmopt(c0,c1,c2,cm,sc0,pme,gde,vpp,eigv)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Deallocations
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
  END SUBROUTINE gmopts
  ! ==================================================================

END MODULE gmopts_utils
