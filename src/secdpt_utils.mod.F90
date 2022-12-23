MODULE secdpt_utils
  USE atwf,                            ONLY: atwp
  USE bsym,                            ONLY: bsfac
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE fnlalloc_utils,                  ONLY: fnlalloc,&
                                             fnldealloc
  USE ions,                            ONLY: ions1
  USE kinds,                           ONLY: real_8
  USE linres,                          ONLY: &
       clrwf, nlinr, nlinw, nolr, nua, nub, td01, td03, urot
  USE lr_in_utils,                     ONLY: lr_in,&
                                             tddft_input
  USE parac,                           ONLY: parai,&
                                             paral
  USE prmem_utils,                     ONLY: prmem
  USE prop,                            ONLY: prop1
  USE sdlinres_utils,                  ONLY: sdlinres
  USE secder_driver,                   ONLY: secder
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             ncpw,&
                                             nkpt
  USE utils,                           ONLY: nxxfun
  USE vibana_utils,                    ONLY: hesfile,&
                                             vibana
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: secdpt

CONTAINS

  ! ==================================================================
  SUBROUTINE secdpt
    ! ==--------------------------------------------------------------==
    ! ==  VIBRATIONAL ANALYSIS                                        ==
    ! ==--------------------------------------------------------------==
    ! Variables
    CHARACTER(*), PARAMETER                  :: procedureN = 'secdpt'

    COMPLEX(real_8), ALLOCATABLE             :: c0(:,:,:), c2(:,:,:), cg0(:), &
                                                gde(:), pme(:), sc0(:)
    COMPLEX(real_8), ALLOCATABLE, TARGET     :: c1(:,:,:,:)
    INTEGER :: i, ierr, lenbk, mm, ms, nact, nc, nc0, nc0_1d, nc0_2d, nc0_3d, &
      nc1, nc2, nc2_1d, nc2_2d, nc2_3d, neig, ngde, npme, nsc0, nstate, nus, &
      nvpp, nxx
    REAL(real_8), ALLOCATABLE                :: eigv(:), sder(:,:), vpp(:)

! TODO refactor this 
! ==--------------------------------------------------------------==
! ==  MEMORY ALLOCATION                                           ==
! ==--------------------------------------------------------------==
! ==                      WAVEFUNCTIONS                           ==
! ==  CG0   : WAVEFUNCTIONS AT EQUILIBRIUM GEOMETRY               ==
! ==  C0    : WAVEFUNCTIONS                                       ==
! ==  C2    : GRADIENTS                                           ==
! ==  SC0   : S**(-1/2)*C0 (nonorthogonal orbitals)               ==
! ==  EIGV  : ORBITAL ENERGIES                                    ==
! ==--------------------------------------------------------------==

    IF ( cntl%tqmmm ) CALL stopgm("SECDPT","QMMM NOT IMPLEMENTED",& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    IF (cntl%tddft) THEN
       CALL lr_in
       CALL tddft_input
    ENDIF
    nc=2*ncpw%ngw*crge%n+8
    nstate=crge%n
    nsc0=1
    nc0_1d = 1
    nc0_2d = 1
    nc0_3d = 1
    IF (cntl%tdiag) THEN
       IF (cntl%tdavi) CALL stopgm('SECDPT',' TDAVI NOT IMPLEMENTED',& 
            __LINE__,__FILE__)
       nc0=(2*nkpt%ngwk*nstate*nkpt%nkpnt+8)*bsfac
       nc0_1d = nkpt%ngwk
       nc0_2d = nstate
       nc0_3d = nkpt%nkpnt*bsfac !vw factor of 2 needed?
       nc2=(2*nkpt%ngwk*MAX(atwp%numaormax,nstate)+8)*bsfac
       nc2_1d = nkpt%ngwk
       nc2_2d = MAX(atwp%numaormax,nstate)
       nc2_3d = bsfac
       neig=nstate*nkpt%nkpts*bsfac
    ELSE
       nc0=(2*nkpt%ngwk*nstate*nkpt%nkpnt+8)*bsfac
       nc0_1d = nkpt%ngwk
       nc0_2d = nstate
       nc0_3d = nkpt%nkpnt*bsfac !vw factor of 2 needed?
       nc2=(2*nkpt%ngwk*MAX(atwp%numaormax,nstate*nkpt%nkpnt)+8)*bsfac
       nc2_1d = nkpt%ngwk
       nc2_2d = MAX(atwp%numaormax,nstate)
       nc2_3d = bsfac*nkpt%nkpnt
       nsc0=MAX(2*nkpt%ngwk*crge%n*nkpt%nkpnt*bsfac,nc2)
       neig=nstate*nkpt%nkpts*bsfac
    ENDIF
    IF (prop1%dberry) THEN
       lenbk=nxxfun(nstate)
       nxx=MAX(2*lenbk*parai%nproc,nc2)
       nsc0=MAX(nxx,nsc0)
       nc2=MAX(nxx,nc2)
    ENDIF
    IF (cntl%tddft) THEN
       DO i=1,nstate
          IF (crge%f(i,1).LT.1.e-6_real_8) CALL stopgm("cntl%tddft",&
               "ALL STATES HAVE TO BEOCCUPIED",& 
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
    ELSE
       nc1=1
       nua=1
       nub=1
       nlinw=0
       nlinr=0
       nus=1
    ENDIF
    ALLOCATE(cg0(nc0),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(c0(nc0_1d,nc0_2d,nc0_3d),STAT=ierr)
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

    ALLOCATE(c2(nc2_1d,nc2_2d,nc2_3d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(sc0(nsc0),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(eigv(neig),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    ALLOCATE(urot(nua,nub,2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    CALL zeroing(urot)!,nus)
    CALL zeroing(c2)

    ! TODO refactor CLRWF
    clrwf => c1
    ! ==--------------------------------------------------------------==
    ! ==  PME   : cntl%diis WF             HNM1 (for CG)                   ==
    ! ==  GDE   : cntl%diis GRADIENTS                                      ==
    ! ==  VPP   : WF HESSIAN                                          ==
    ! ==--------------------------------------------------------------==
    IF (.NOT.cntl%tsdan) THEN
       IF (cntl%tsde) THEN
          npme = 1
          ngde = 1
          nvpp = 1
       ELSE IF (cntl%diis) THEN
          npme = (nkpt%ngwk*nstate+8)*cnti%mdiis*nkpt%nkpnt
          ngde = ((nkpt%ngwk*nstate+8)*cnti%mdiis*nkpt%nkpnt)/4
          nvpp = nkpt%ngwk
       ELSE IF (cntl%pcg) THEN
          npme = 2*nkpt%ngwk*nstate*nkpt%nkpnt
          ngde = 1
          nvpp = nkpt%ngwk
       ELSE IF (cntl%tdiag) THEN
          IF (cntl%tfrsblk) THEN
             npme=2*nkpt%ngwk*MAX(nstate,(cnti%nkry_max+1)*cnti%nkry_block)
          ELSE
             npme=2*nkpt%ngwk*MAX(nstate,cnti%nkry_max*cnti%nkry_block)
          ENDIF
          ngde = 1
          nvpp = nkpt%ngwk
       ELSE
          IF (paral%io_parent)&
               WRITE(6,*) ' WRONG OPTION FOR WAVEFUNCTION OPTIMIZATION'
          CALL stopgm('SECDPT',' ',& 
               __LINE__,__FILE__)
       ENDIF
       ALLOCATE(pme(npme),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(gde(ngde),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(vpp(nvpp),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! FNL..DFNL
    CALL fnlalloc(crge%n,.TRUE.,.FALSE.)
    IF (paral%parent) CALL prmem('    SECDPT')
    ! ==--------------------------------------------------------------==
    ! ..from a stored Hessian
    IF (cntl%tsdin) THEN
       IF (cntl%bsymm)CALL stopgm('SECDPT',&
            'BSYMM AND STORED_HESSIAN NOT IMPLEMENTED',& 
            __LINE__,__FILE__)

       IF (paral%parent)  THEN
          ALLOCATE(sder(3*ions1%nat,3*ions1%nat),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ELSE
          ! Avoid 'not allocated' runtime error
          ALLOCATE(sder(1,1),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
       ENDIF
       CALL hesfile(sder,eigv,c0,c2,sc0)
       CALL vibana(sder)
       IF (paral%parent) DEALLOCATE(sder,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ELSE IF (cntl%tsdan) THEN
       ! ..linear response
       IF (cntl%tddft) CALL stopgm("SECDPT","cntl%tddft NOT IMPLEMENTED",& 
            __LINE__,__FILE__)
       IF (cntl%bsymm)CALL stopgm('SECDPT',&
            'BSYMM AND LINEAR RESPONSE NOT IMPLEMENTED',& 
            __LINE__,__FILE__)
       CALL sdlinres(c0,c2,sc0,cg0,eigv)
    ELSE
       ! ..finite differences
       IF (paral%io_parent)&
            WRITE(6,'(A)') ' START FINITE DIFFERENCES'
       CALL secder(c0,c1(:,:,:,1),c2,sc0,cg0,pme,gde,vpp,eigv)
    ENDIF
    ! ==--------------------------------------------------------------==
    DEALLOCATE(cg0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(c0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(c1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(c2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(eigv,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(sc0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(urot,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    CALL fnldealloc(.TRUE.,.FALSE.)
    IF (.NOT.cntl%tsdan) THEN
       DEALLOCATE(pme,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(gde,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(vpp,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE secdpt
  ! ==================================================================

END MODULE secdpt_utils
