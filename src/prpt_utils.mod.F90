MODULE prpt_utils
  USE atwf,                            ONLY: atwp
  USE ddip,                            ONLY: lenbk
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE fnlalloc_utils,                  ONLY: fnlalloc,&
                                             fnldealloc
  USE jrotation_utils,                 ONLY: set_orbdist
  USE kinds,                           ONLY: real_8
  USE npt_md_utils,                    ONLY: npt_bomd,&
                                             npt_cpmd
  USE parac,                           ONLY: parai,&
                                             paral
  USE prbomd_utils,                    ONLY: prbomd
  USE prcpmd_utils,                    ONLY: prcpmd
  USE prmdfile_utils,                  ONLY: prmdfile
  USE prmem_utils,                     ONLY: prmem
  USE pslo,                            ONLY: pslo_com
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             ncpw,&
                                             nkpt
  USE utils,                           ONLY: nxxfun
  USE vdwcmod,                         ONLY: vdwl
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: prpt

CONTAINS

  ! ==================================================================
  SUBROUTINE prpt
    ! ==--------------------------------------------------------------==
    ! Variables
    CHARACTER(*), PARAMETER                  :: procedureN = 'prpt'

    COMPLEX(real_8), ALLOCATABLE             :: c0(:), c1(:), c2(:,:), cm(:), &
                                                sc0(:)
    INTEGER                                  :: ierr, nc, nc1, nc2, nc2_1d, &
                                                nc2_2d, ncm, ncs, nstate, &
                                                nstx, nvpp
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

    IF (cntl%tqmmm) CALL stopgm('PRPT','QM/MM NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    IF (cntl%bsymm) CALL stopgm('PRPT','BROKEN SYMMETRY NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    nstate = crge%n
    ! ==--------------------------------------------------------------==
    nc=ncpw%ngw*nstate+8
    ncm=nc
    nc1=1
    nvpp=8
    lenbk=0
    IF (cntl%tdipd.OR.vdwl%vdwd) THEN
       lenbk=nxxfun(nstate)
       nc2=MAX(2*lenbk*parai%nproc,nc)
       nc2_1d = ncpw%ngw
       nc2_2d = MAX(nstate, lenbk*parai%nproc/ncpw%ngw )
       ncs=nc2
    ELSE
       nc=ncpw%ngw*crge%n+8
       nc2=ncpw%ngw*crge%n+8
       nc2_1d = ncpw%ngw
       nc2_2d = crge%n
       ncs=1
       IF (pslo_com%tivan.OR.cntl%nonort) ncs=nc
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (cntl%tmdbo.OR.cntl%tmdfile) THEN
       nvpp=2*nkpt%ngwk
       IF (cntl%tdiag) THEN
          IF (cntl%tdavi) THEN
             IF (cnti%ndavv.EQ.0) THEN
                cnti%ndavv=2*nstate
             ELSE
                cnti%ndavv=MIN(cnti%ndavv+nstate,2*nstate)
             ENDIF
             nstate=cnti%ndavv
          ENDIF
          IF (cntl%tfrsblk)THEN
             ncm=nkpt%ngwk*MAX(nstate,(cnti%nkry_max+1)*cnti%nkry_block)+8
          ELSEIF (cntl%diis) THEN
             ncm=(nkpt%ngwk*nstate+8)*cnti%mdiis*nkpt%nkpnt+&
                  ((2*nkpt%ngwk*nstate+8)*cnti%mdiis*nkpt%nkpnt)/4
          ELSE
             ncm=nkpt%ngwk*MAX(nstate,cnti%nkry_max*cnti%nkry_block)+8
          ENDIF
          nc2=nkpt%ngwk*MAX(atwp%numaormax,nstate)+8
          nc2_1d = nkpt%ngwk
          nc2_2d = MAX(atwp%numaormax,nstate)
          IF (cntl%tdavi)ncm=ncm*2
       ELSE 
          IF (cntl%diis) THEN
             ncm=(nkpt%ngwk*nstate+8)*cnti%mdiis+((nkpt%ngwk*nstate+8)*cnti%mdiis)/4+100
          ELSEIF (cntl%pcg) THEN
             ncm=2*nkpt%ngwk*nstate*nkpt%nkpnt+8
          ENDIF
          nc=2*nkpt%ngwk*nstate*nkpt%nkpnt+8
          nc2=MAX(2*nkpt%ngwk*nstate+8,nc2)
          nc2_1d = nkpt%ngwk
          nc2_2d = nstate
          ncs=nc2
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    ALLOCATE(c0(nc),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(c0)
    ALLOCATE(c1(nc1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(c1)
    ALLOCATE(c2(nc2_1d,nc2_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(c2)
    ALLOCATE(cm(ncm),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(cm)
    ALLOCATE(sc0(ncs),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(sc0)
    ALLOCATE(vpp(nvpp),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(vpp)
    IF (cntl%tdmal) THEN
       CALL set_orbdist(nstate,cnti%nstblk,parai%nproc,nstx)
       ALLOCATE(gamx(nstate*nstx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(gamy(nstate*nstx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ELSE
       ALLOCATE(gamx(nstate*nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(gamy(nstate*nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! FNL
    CALL fnlalloc(crge%n,.TRUE.,.TRUE.)
    IF (paral%parent) CALL prmem('      MDPT')
    ! ==--------------------------------------------------------------==
    IF (cntl%tmdbo) THEN
       IF (cntl%tpenn.OR.cntl%tshock) THEN
          CALL npt_bomd(c0,cm,c2,sc0,vpp,gamx,gamy)
       ELSE
          CALL prbomd(c0,cm,c2,sc0,vpp,gamx,gamy)
       ENDIF
    ELSEIF(cntl%tmdfile) THEN
       CALL prmdfile(c0,cm,c1,c2,sc0,vpp,gamx,gamy)
    ELSE
       IF (cntl%tpenn.OR.cntl%tshock) THEN
          CALL npt_cpmd(c0,cm,c2,sc0,gamx,gamy)
       ELSE
          CALL prcpmd(c0,cm,c2,sc0,gamx,gamy)
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    DEALLOCATE(c0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(cm,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(c1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(c2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(sc0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(vpp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(gamx,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(gamy,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    CALL fnldealloc(.TRUE.,.TRUE.)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE prpt
  ! ==================================================================
END MODULE prpt_utils
