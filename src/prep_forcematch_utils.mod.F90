MODULE prep_forcematch_utils
  USE bsym,                            ONLY: bsfac
  USE ddip,                            ONLY: lenbk
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE fnlalloc_utils,                  ONLY: fnlalloc,&
                                             fnldealloc
  USE jrotation_utils,                 ONLY: set_orbdist
  USE kinds,                           ONLY: real_8
  USE linres,                          ONLY: clrwf,&
                                             nlinr,&
                                             nlinw,&
                                             nua,&
                                             nub,&
                                             urot
  USE mm_dim_utils,                    ONLY: mm_dim
  USE mm_dimmod,                       ONLY: mm_go_qm,&
                                             mm_revert
  USE mm_forcematch_utils,             ONLY: mm_forcematch
  USE parac,                           ONLY: parai,&
                                             paral
  USE prmem_utils,                     ONLY: prmem
  USE pslo,                            ONLY: pslo_com
  USE shop_rest_2,                     ONLY: c0old,&
                                             c0old2
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             nkpt
  USE utils,                           ONLY: nxxfun
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: prep_forcematch

CONTAINS

  ! ==================================================================
  SUBROUTINE prep_forcematch
    ! ==--------------------------------------------------------------==
    ! Variables
    CHARACTER(*), PARAMETER                  :: procedureN = 'prep_forcematch'

    COMPLEX(real_8), ALLOCATABLE             :: c0(:), c2(:), cm(:), sc0(:)
    COMPLEX(real_8), ALLOCATABLE, TARGET     :: c1(:,:,:,:)
    INTEGER                                  :: ierr, nc0, nc1, nc2, ncm, &
                                                nmin, nsc0, nstate, nstx, &
                                                nus, nxx
    LOGICAL                                  :: statusdummy
    REAL(real_8), ALLOCATABLE                :: gamx(:), gamy(:), vpp(:)

! defined 4D for
! CLRWF pointer
! assignment
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
       nstate=crge%n
       nc1=1
       nua=1
       nub=1
       nlinw=0
       nlinr=0
       nus=1

       IF (cntl%tsde) THEN
          ncm=8
       ELSE IF (cntl%diis) THEN
          ncm=(nkpt%ngwk*nstate+8)*cnti%mdiis*nkpt%nkpnt+&
               ((nkpt%ngwk*nstate+8)*cnti%mdiis*nkpt%nkpnt)/4+100
       ELSE IF (cntl%pcg) THEN
          ncm=2*nkpt%ngwk*nstate*nkpt%nkpnt+8
       ENDIF
       ! CB: Need memory for two wf in BS case          
       nc0=(2*nkpt%ngwk*nstate*nkpt%nkpnt+8)*bsfac
       nc2=(2*nkpt%ngwk*nstate+8)*bsfac
       ALLOCATE(vpp(nkpt%ngwk),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       nsc0=1
       IF (pslo_com%tivan.OR.cntl%nonort.OR.(cntl%quenchb.AND.cntl%pcg.AND.cntl%pcgmin)) THEN
          nsc0=2*nkpt%ngwk*nstate*nkpt%nkpnt*bsfac
       ENDIF
       IF (cntl%tdipd) THEN
          lenbk=nxxfun(nstate)
          nxx=MAX(2*lenbk*parai%nproc,nc2)
          nsc0=MAX(nxx,nsc0)
          nc2=nxx
       ENDIF
       nmin=10
       ! CALL MEMORY(IP_C0OLD,NMIN,'C0OLD')
       ! CALL MEMORY(IP_C0OLD2,NMIN,'C0OLD2')
       ! CALL MEMORY(IP_C0,NC0,'C0')
       ! CALL MEMORY(IP_CM,NCM,'CM')
       ! CALL MEMORY(IP_C1,NC1,'C1')
       ! CALL MEMORY(IP_C2,NC2,'C2')
       ! CALL MEMORY(IP_SC0,NSC0,'SC0')
       ! CALL MEMORY(IP_UROT,NUS,'UROT')
       ALLOCATE(c0old(nmin,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(c0old2(nmin,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(c0(nc0),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(cm(ncm),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(c2(nc2),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(sc0(nsc0),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(urot(nus,1,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)! UROT should be 3D,
       ! but here we know only
       ! total size required
       CALL zeroing(urot)!,nus)
       ALLOCATE(c1(nc1,1,1,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       clrwf => c1
       IF (cntl%tdmal) THEN
          CALL set_orbdist(nstate,cnti%nstblk,parai%nproc,nstx)
          ! CALL MEMORY(IP_GAMX,NSTATE*NSTX*BSFAC,'GAMX') 
          ! CALL MEMORY(IP_GAMY,NSTATE*NSTX*BSFAC,'GAMY') 
          ALLOCATE(gamx(nstate*nstx*bsfac),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(gamy(nstate*nstx*bsfac),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ELSE
          ! CALL MEMORY(IP_GAMX,NSTATE*NSTATE*BSFAC,'GAMX') 
          ! CALL MEMORY(IP_GAMY,NSTATE*NSTATE*BSFAC,'GAMY') 
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
       ! 
       CALL fnlalloc(nstate,.TRUE.,cntl%tpres)
       CALL mm_dim(mm_revert,statusdummy)
    ELSE
       nmin = 10
       ! CALL MEMORY(IP_C0,NMIN,'C0')
       ! CALL MEMORY(IP_CM,NMIN,'CM')
       ! CALL MEMORY(IP_C1,NMIN,'C1')
       ! CALL MEMORY(IP_C2,NMIN,'C2')
       ! CALL MEMORY(IP_SC0,NMIN,'SC0')
       ! CALL MEMORY(IP_UROT,NMIN,'UROT')
       ! CALL MEMORY(IP_GAMX,NMIN,'GAMX')
       ! CALL MEMORY(IP_GAMY,NMIN,'GAMY')
       ! CALL MEMORY(IP_VPP,NMIN,'VPP')
       ! 
       ! CALL MEMORY(IP_C0OLD,NMIN,'C0OLD')
       ! CALL MEMORY(IP_C0OLD2,NMIN,'C0OLD2')
       ALLOCATE(c0(nmin),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(cm(nmin),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(c1(nmin,1,1,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)! C1 is defined 4D for
       ! CLRWF pointer
       ! assignemnt
       ALLOCATE(c2(nmin),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(sc0(nmin),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(urot(nmin,1,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(gamx(nmin),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(gamy(nmin),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(vpp(nmin),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)

       ALLOCATE(c0old(nmin,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(c0old2(nmin,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ! 
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (paral%parent) CALL prmem('      MDPT')
    ! ==--------------------------------------------------------------==
    ! 
    IF (cntl%tqmmm) THEN
       CALL mm_forcematch(c0,cm,c1,c2,sc0,vpp,gamx,gamy)
    ELSE
       ! CALL FORCEMATCH(C0,CM,C1,C2,SC0,VPP,GAMX,GAMY)
       ! non-QMMM force matching not implemented yet
       CALL stopgm('prep_forcematch', 'force matching only with QMMM! ',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    ! CALL FREEM(IP_C0)
    ! CALL FREEM(IP_CM)
    ! CALL FREEM(IP_C1)
    ! CALL FREEM(IP_C2)
    ! CALL FREEM(IP_SC0)
    ! CALL FREEM(IP_UROT)
    ! CALL FREEM(IP_GAMX)
    ! CALL FREEM(IP_GAMY)
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
    DEALLOCATE(urot,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(gamx,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(gamy,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF ( paral%qmnode ) THEN
       CALL fnldealloc(.TRUE.,cntl%tpres)
    ELSE
       ! CALL FREEM(IP_VPP)
       DEALLOCATE(vpp,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE prep_forcematch
  ! ==================================================================
END MODULE prep_forcematch_utils
