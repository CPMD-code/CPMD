MODULE interpt_utils
  USE ddip,                            ONLY: lenbk
  USE egointer_utils,                  ONLY: INTERFACE
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE fnlalloc_utils,                  ONLY: fnlalloc,&
                                             fnldealloc
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: parai,&
                                             paral
  USE prmem_utils,                     ONLY: prmem
  USE pslo,                            ONLY: pslo_com
  USE rwfopt_utils,                    ONLY: rwfopt
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             ncpw,&
                                             nkpt
  USE utils,                           ONLY: nxxfun

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: interpt

CONTAINS

  ! ==================================================================
  SUBROUTINE interpt
    ! ==--------------------------------------------------------------==
    ! Variables
    CHARACTER(*), PARAMETER                  :: procedureN = 'interpt'

    COMPLEX(real_8), ALLOCATABLE             :: c0(:,:,:), c2(:,:), gde(:), &
                                                pme(:), sc0(:,:,:)
    INTEGER                                  :: ierr, nc, nc2, nkpntsc0, &
                                                nsc0, nxx
    LOGICAL                                  :: wasdiis
    REAL(real_8), ALLOCATABLE                :: eigv(:), vpp(:)

! ==--------------------------------------------------------------==
! ==  MEMORY ALLOCATION                                           ==
! ==--------------------------------------------------------------==
! ==  C0    : WAVEFUNCTIONS                                       ==
! ==  C2    : GRADIENTS                                           ==
! ==  SC0   : S**(-1/2)*C0 (nonorthogonal orbitals)               ==
! ==  PME   : cntl%diis WF             HNM1 (for CG)                   ==
! ==  GDE   : cntl%diis GRADIENTS                                      ==
! ==  VPP   : WF HESSIAN                                          ==
! ==  EIGV  : ORBITAL ENERGIES                                    ==
! ==--------------------------------------------------------------==

    IF ( cntl%tqmmm ) CALL stopgm("INTERPT","QMMM NOT IMPLEMENTED",& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    nc=crge%n
    nc2=nc
    nsc0=1
    nkpntsc0=1
    IF (pslo_com%tivan.OR.cntl%nonort.OR.(cntl%pcg.AND.cntl%pcgmin)) THEN
       nsc0=crge%n
       nkpntsc0=nkpt%nkpnt
    ENDIF
    IF (cntl%tdipd) THEN
       IF (cntl%tddft) CALL stopgm("INTERPT",&
            "cntl%tddft AND DIPOLE DYNAMICS NOT POSSIBLE",& 
            __LINE__,__FILE__)
       lenbk=nxxfun(crge%n)
       nxx=MAX(lenbk*parai%nproc / nkpt%ngwk ,nc2)
       nsc0=MAX(nxx,nsc0)
       nc2=nxx
    ENDIF

    ALLOCATE(c0(nkpt%ngwk,nc,nkpt%nkpnt),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(c2(nkpt%ngwk,nc2),STAT=ierr)
    !  allocate(c2(nc2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(eigv(crge%n),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(sc0(nkpt%ngwk,nsc0,nkpntsc0),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    IF (cntl%tsde) THEN
       ! TODO check this DUMMY 

       ALLOCATE(vpp(ncpw%ngw),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ELSE IF (cntl%diis) THEN
       ALLOCATE(pme((ncpw%ngw*crge%n+8)*cnti%mdiis/2),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(gde(((ncpw%ngw*crge%n+8)*cnti%mdiis)/2),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(vpp(ncpw%ngw),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ELSE IF (cntl%pcg) THEN
       ALLOCATE(pme(2*ncpw%ngw*crge%n),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)

       ! TODO check this DUMMY 

       ALLOCATE(vpp(ncpw%ngw),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ..FNL and DFNL
    CALL fnlalloc(crge%n,.TRUE.,.FALSE.)
    IF (paral%parent) CALL prmem('   INTERPT')
    ! ==--------------------------------------------------------------==
    IF (cntl%diis.AND.cntl%tpcgfi) THEN
       IF (paral%io_parent)&
            WRITE(6,*) 'INTERPT| SWITCHING TO PCG MINIMIZE'
       wasdiis=.TRUE.
       cntl%diis=.FALSE.
       cntl%pcg=.TRUE.
       cntl%pcgmin=.TRUE.
    ELSE
       wasdiis=.FALSE.
    ENDIF
    CALL rwfopt(c0,c2,sc0,pme,gde,vpp,eigv)
    IF (wasdiis) THEN
       IF (paral%io_parent)&
            WRITE(6,*) 'INTERPT| SWITCHING BACK TO DIIS'
       cntl%diis=.TRUE.
       cntl%pcg=.FALSE.
       cntl%pcgmin=.FALSE.
    ENDIF
    IF (paral%io_parent)&
         CLOSE(5)
    ! ==--------------------------------------------------------------==
    CALL INTERFACE(c0,c2,sc0,pme,gde,vpp,eigv)
    ! ==--------------------------------------------------------------==
    DEALLOCATE(c0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(c2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(eigv,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    CALL fnldealloc(.TRUE.,.FALSE.)
    IF (cntl%nonort.OR.pslo_com%tivan.OR.(cntl%pcg.AND.cntl%pcgmin)) THEN
       DEALLOCATE(sc0,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    IF (cntl%diis) THEN
       DEALLOCATE(pme,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(gde,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    DEALLOCATE(vpp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE interpt
  ! ==================================================================

END MODULE interpt_utils
