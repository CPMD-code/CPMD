MODULE wfopts_utils
  USE atwf,                            ONLY: atwp
  USE bswfo_utils,                     ONLY: bs_wfo
  USE cl_init_utils,                   ONLY: cl_init
  USE clas,                            ONLY: tclas
  USE coor,                            ONLY: tau0
  USE cplngsmod,                       ONLY: tcpl,&
                                             tcpllr
  USE ddip,                            ONLY: lenbk
  USE do_perturbation_p_utils,         ONLY: do_perturbation
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE fnlalloc_utils,                  ONLY: fnlalloc,&
                                             fnldealloc
  USE ions,                            ONLY: ions1
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE lr_in_utils,                     ONLY: lr_in
  USE nlcc,                            ONLY: corel,&
                                             vnlcc,&
                                             vnlt
  USE parac,                           ONLY: parai,&
                                             paral
  USE phfac_utils,                     ONLY: phfac
  USE prmem_utils,                     ONLY: prmem
  USE pslo,                            ONLY: pslo_com
  USE rwfopt_utils,                    ONLY: rwfopt
  USE sfac,                            ONLY: eigr,&
                                             eigrb
  USE soft,                            ONLY: soft_com
  USE spin,                            ONLY: lspin2
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

  PUBLIC :: wfopts

CONTAINS

  ! ==================================================================
  SUBROUTINE wfopts
    ! ==--------------------------------------------------------------==
    ! Variables
    CHARACTER(*), PARAMETER                  :: procedureN = 'wfopts'

    COMPLEX(real_8), ALLOCATABLE             :: c0(:,:,:), c2(:,:), gde(:), &
                                                pme(:), sc0(:,:,:)
    INTEGER                                  :: ierr, isub, iunit, nc0, nc2, &
                                                nce, neig, ngde, nkptc0, &
                                                nkptsc0, npme, nsc0, nstate, &
                                                nvpp, nxx
    REAL(real_8), ALLOCATABLE                :: eigv(:), vpp(:)

! ==--------------------------------------------------------------==
! ==  GDE   : cntl%diis GRADIENTS                 CR   (for DIAG)      ==
! ==  VPP   : WF HESSIAN                                          ==
! ==  EIGV  : ORBITAL ENERGIES                                    ==
! ==--------------------------------------------------------------==
! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)


    nstate=crge%n
    ! Minimum requirements for all methods
    nc0  = nstate
    nkptc0 = nkpt%nkpts
    nc2  = MAX(atwp%numaormax,nstate)
    nsc0 = 1
    nkptsc0 = 1
    nvpp = 1
    npme = 1
    ngde = 1
    neig = nstate
    ! Special requirements for diagonalization of wavefunction
    IF (cntl%tdiag.OR.cntl%tdiagopt) THEN
       IF (cntl%tdavi) THEN
          IF (cnti%ndavv.EQ.0) THEN
             cnti%ndavv=2*nstate
          ELSE
             cnti%ndavv=MIN(cnti%ndavv+nstate,2*nstate)
          ENDIF
          nc0  = MAX (nc0,cnti%ndavv)
          nc2  = MAX (nc2,cnti%ndavv)
          nvpp = MAX (nvpp, nkpt%ngwk)
          npme = MAX (npme, nkpt%ngwk*cnti%ndavv)
          ngde = MAX (ngde, nkpt%ngwk*cnti%ndavv)
          neig = MAX (neig, cnti%ndavv)
       ELSE IF (cntl%tfrho_upw) THEN
          neig = MAX (neig, nstate*nkpt%nkpts)
          IF (cntl%diis) THEN
             npme = MAX (npme, (nkpt%ngwk*nstate+8)*cnti%mdiis*nkpt%nkpnt)
             ngde = MAX (ngde, ((nkpt%ngwk*nstate+8)*cnti%mdiis*nkpt%nkpnt)/4)
             nvpp = MAX (nvpp, nkpt%ngwk)
          ELSE IF (cntl%pcg) THEN
             npme = MAX (npme, nkpt%ngwk*nstate*nkpt%nkpnt)
             nvpp = MAX (nvpp, nkpt%ngwk)
          ENDIF
       ELSE
          IF (cntl%tfrsblk) THEN
             nce = 2*nkpt%ngwk*MAX(nstate,(cnti%nkry_max+1)*cnti%nkry_block)
          ELSE
             nce = 2*nkpt%ngwk*MAX(nstate,cnti%nkry_max*cnti%nkry_block)
          ENDIF
          npme = MAX (npme, nce)
          nvpp = MAX (nvpp, ncpw%ngw)
          neig = MAX (neig, nstate*nkpt%nkpts)
       ENDIF
    ENDIF
    ! Minimum scratch space if used at all
    IF (pslo_com%tivan.OR.cntl%nonort.OR.(cntl%pcg.AND.cntl%pcgmin).OR.cntl%diis.OR.cntl%tfrho_upw) THEN
       nsc0 = MAX (nsc0,nstate)
       IF (cntl%tdavi)nsc0=MAX(nsc0,cnti%ndavv)
    ENDIF
    IF (cntl%tdipd.OR.vdwl%vdwd) THEN
       lenbk=nxxfun(nstate)
       nxx=MAX(2*lenbk*parai%nproc/nkpt%ngwk,nc2 ) !vw nkpt%ngwk added while converting (*) -> (:)
       nsc0=MAX(nxx,nsc0)
       nc2=MAX(nxx,nc2)
    ENDIF
    ! Special requirements for optimization of wavefunction
    IF (.NOT.cntl%tdiag) THEN
       neig = MAX (neig, nstate*nkpt%nkpts)
       IF (cntl%nonort) THEN
          nsc0 = MAX (nsc0, nstate)
          nkptsc0 = nkpt%nkpnt
       ELSEIF (tkpts%tkpnt) THEN
          nsc0 = MAX (nsc0, nc2)
       ENDIF
       IF (cntl%tsde) THEN
          nvpp = MAX (nvpp, nkpt%ngwk)
       ELSE IF (cntl%diis) THEN
          npme = MAX (npme, (nkpt%ngwk*nstate+8)*cnti%mdiis*nkpt%nkpnt)
          ngde = MAX (ngde, ((nkpt%ngwk*nstate+8)*cnti%mdiis*nkpt%nkpnt)/4)
          nvpp = MAX (nvpp, nkpt%ngwk)
       ELSE IF (cntl%pcg) THEN
          npme = MAX (npme, nkpt%ngwk*nstate*nkpt%nkpnt)
          nvpp = MAX (nvpp, nkpt%ngwk)
       ELSE
          IF (paral%io_parent)&
               WRITE(6,*) ' WRONG OPTION FOR WAVEFUNCTION OPTIMIZATION'
          CALL stopgm('WFOPTS',' OPTION      ',& 
               __LINE__,__FILE__)
       ENDIF
    ENDIF
    ! Broken symmetry state
    IF (cntl%bsymm) THEN
       nkptc0 = 2
       neig = 2*neig
    ENDIF
    ! Linear response needs KS matrix for LSE (no canonical states)
    IF (lspin2%tlse.AND.tcpl.AND.tcpllr) neig = nstate*nstate*nkpt%nkpts
    ! Allocate the memory
    ALLOCATE(c0(nkpt%ngwk,nc0,nkptc0),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(c0)
    ALLOCATE(c2(nkpt%ngwk,nc2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(c2)
    ALLOCATE(sc0(nkpt%ngwk,nsc0,nkptsc0),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(sc0)
    ALLOCATE(pme(npme),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(gde(ngde),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(vpp(nvpp),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(eigv(neig),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(eigv)!,neig)
    ! ==--------------------------------------------------------------==
    IF (tclas) THEN
       iunit=5
       CALL cl_init(iunit)
    ENDIF
    IF (cntl%tdavi) nstate=cnti%ndavv
    ! FNL and DFNL
    CALL fnlalloc(nstate,.TRUE.,cntl%tpres)
    IF (paral%parent) CALL prmem('    WFOPTS')
    ! ==--------------------------------------------------------------==
    IF (tcpl) THEN
       IF (tcpllr) CALL lr_in
    ENDIF
    ! ==--------------------------------------------------------------==
    ! CB
    IF (.NOT.cntl%bsymm) THEN
       CALL rwfopt(c0,c2,sc0,pme,gde,vpp,eigv)
    ELSE
       CALL bs_wfo(c0,c2,sc0,pme,gde,vpp,eigv)
    ENDIF
    ! ==--------------------------------------------------------------==

    DEALLOCATE(eigv,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(sc0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(pme,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(gde,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(vpp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    ! The following is due to a bad placing of the EIGR
    ! arrays in memory: By deallocating and re-allocating,
    ! they are placed at the first available memory location.
    IF (cntl%tresponse .AND. .NOT. soft_com%exsoft) THEN
       DEALLOCATE(eigr,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(eigr(ncpw%ngw,ions1%nat,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       IF (cntl%bigmem) THEN
          DEALLOCATE(eigrb,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(eigrb(ncpw%nhg,ions1%nat),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
       CALL phfac(tau0)

       CALL do_perturbation(c0,c2,nstate)
    ENDIF

    ! Deallocate vnlt and vnlcc allocated in initrun routine
    IF (corel%tinlc) THEN
       DEALLOCATE(vnlt,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(vnlcc,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    CALL fnldealloc(.TRUE.,cntl%tpres)
    DEALLOCATE(c0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(c2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
  END SUBROUTINE wfopts
  ! ==================================================================
END MODULE wfopts_utils
