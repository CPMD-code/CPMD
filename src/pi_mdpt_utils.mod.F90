MODULE pi_mdpt_utils
  USE atwf,                            ONLY: atwp
  USE ddip,                            ONLY: lenbk
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE fnlalloc_utils,                  ONLY: fnlalloc,&
                                             fnldealloc
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE mp_interface,                    ONLY: mp_sync
  USE parac,                           ONLY: parai,&
                                             paral
  USE pi_diag_utils,                   ONLY: pi_diag
  USE pi_md_utils,                     ONLY: pi_md
  USE pimd,                            ONLY: np_local,&
                                             pc_groups,&
                                             supergroup
  USE prmem_utils,                     ONLY: prmem
  USE pslo,                            ONLY: pslo_com
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             ncpw,&
                                             nkpt
  USE testex_utils,                    ONLY: testex
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: nxxfun
  USE vdwcmod,                         ONLY: vdwl
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: pi_mdpt

CONTAINS

  ! ==================================================================
  SUBROUTINE pi_mdpt
    ! ==--------------------------------------------------------------==
    ! Variables

    CHARACTER(*), PARAMETER                  :: procedureN = 'pi_mdpt'

    COMPLEX(real_8), ALLOCATABLE             :: c0(:,:,:), c2(:), cm(:), &
                                                sc0(:)
    INTEGER                                  :: ierr, isub, nc0, nc0_1d, &
                                                nc0_2d, nc0_3d, nc2, ncm, &
                                                nsc0, nstate, nxx
    LOGICAL                                  :: exsoft
    REAL(real_8), ALLOCATABLE                :: gamx(:), gamy(:), vpp(:)

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

    CALL tiset(procedureN,isub)
    IF ( cntl%tqmmm ) CALL stopgm("PI_MDPT","QMMM NOT IMPLEMENTED",& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    nstate=crge%n
    IF (cntl%tdiag) THEN
       IF (cntl%tdavi) THEN
          IF (cnti%ndavv.EQ.0) THEN
             cnti%ndavv=2*nstate
          ELSE
             cnti%ndavv=MIN(cnti%ndavv+nstate,2*nstate)
          ENDIF
          nstate=cnti%ndavv
       ENDIF
       nc0_1d = 2 * nkpt%ngwk ! Factor 2 for consistentcy with mdpt
       nc0_2d = nstate
       nc0_3d = np_local*nkpt%nkpnt
       nc0 = nc0_1d*nc0_2d*nc0_3d
       IF (cntl%tfrsblk)THEN
          ncm=nkpt%ngwk*MAX(nstate,(cnti%nkry_max+1)*cnti%nkry_block)+8
       ELSE
          ncm=nkpt%ngwk*MAX(nstate,cnti%nkry_max*cnti%nkry_block)+8
       ENDIF
       nc2=nkpt%ngwk*MAX(atwp%numaormax,nstate)+8
       IF (cntl%tdavi) ncm=ncm*2
       IF (cntl%tlanc) ncm=MAX(ncm,ncpw%ngw*8)
       ALLOCATE(vpp(nkpt%ngwk),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ELSEIF (cntl%tmdbo) THEN
       IF (cntl%tsde) THEN
          ncm=8
       ELSE IF (cntl%diis) THEN
          ncm=(nkpt%ngwk*nstate+8)*cnti%mdiis+((nkpt%ngwk*nstate+8)*cnti%mdiis)/4+100
       ELSE IF (cntl%pcg) THEN
          ncm=nkpt%ngwk*nstate
       ENDIF
       nc0_1d = nkpt%ngwk
       nc0_2d = nstate
       nc0_3d = np_local*nkpt%nkpnt
       nc0 = nc0_1d*nc0_2d*nc0_3d
       !!nc0=np_local*(nkpt%ngwk*nstate*nkpt%nkpnt+8)
       nc2=nkpt%ngwk*nstate+8
       ALLOCATE(vpp(nkpt%ngwk),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ELSE
       nc0_1d = nkpt%ngwk
       nc0_2d = nstate
       nc0_3d = np_local*nkpt%nkpnt
       nc0 = nc0_1d*nc0_2d*nc0_3d
       !!nc0=np_local*(nkpt%ngwk*nstate*nkpt%nkpnt+8)
       ncm=np_local*(nkpt%ngwk*nstate*nkpt%nkpnt+8)
       nc2=np_local*(nkpt%ngwk*nstate+8)
    ENDIF
    nsc0=1
    IF (pslo_com%tivan.OR.cntl%nonort.OR.(cntl%quenchb.AND.cntl%pcg.AND.cntl%pcgmin)) THEN
       nsc0=nkpt%ngwk*nstate*nkpt%nkpnt
    ENDIF
    IF (cntl%tdipd.OR.vdwl%vdwd) THEN
       IF (cntl%tddft) CALL stopgm("PI_MDPT","cntl%tddft AND DIPOLE DYNAMICS NOT POSSIBLE",& 
            __LINE__,__FILE__)
       lenbk=nxxfun(nstate)
       nxx=MAX(2*lenbk*parai%nproc,nc2)
       nsc0=MAX(nxx,nsc0)
       nc2=nxx
    ENDIF
    ALLOCATE(c0(nc0_1d,nc0_2d,nc0_3d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(cm(ncm),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(c2(nc2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ALLOCATE(sc0(ncpw%ngw*nstate),STAT=ierr)
    ALLOCATE(sc0(nsc0),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(gamx(nstate*nstate*np_local),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(gamy(nstate*nstate*np_local),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(c0)
    CALL zeroing(cm)
    CALL zeroing(c2)
    CALL zeroing(sc0)
    CALL zeroing(gamx)
    CALL zeroing(gamy)
    ! ==--------------------------------------------------------------==
    ! ..FNL
    CALL fnlalloc(nstate,.TRUE.,.FALSE.)
    IF (paral%parent) CALL prmem('   PI_MDPT')
    ! ==--------------------------------------------------------------==
    IF (cntl%tmdbo) THEN
       CALL pi_diag(c0,cm,c2,sc0,vpp,gamx,gamy)
    ELSE
       IF (tkpts%tkpnt)CALL stopgm('PI_MDPT','K-POINTS NOT IMPLEMENTED',& 
            __LINE__,__FILE__)
       CALL pi_md(c0,cm,c2,sc0,gamx,gamy)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! clean up EXIT file for processor groups
    CALL mp_sync(supergroup)
    pc_groups=1
    CALL testex(exsoft)
    ! ==--------------------------------------------------------------==
    DEALLOCATE(c0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(cm,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(c2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(sc0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (cntl%tdiag.OR.cntl%tmdbo) DEALLOCATE(vpp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(gamx,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(gamy,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    CALL fnldealloc(.TRUE.,.FALSE.)
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    RETURN
  END SUBROUTINE pi_mdpt
  ! ==================================================================

END MODULE pi_mdpt_utils
