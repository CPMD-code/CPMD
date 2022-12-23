MODULE mixing_r_utils
  USE anderson_utils,                  ONLY: anderson,&
                                             change_xmix
  USE andr,                            ONLY: andr2,&
                                             andr3
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE rhodiis_utils,                   ONLY: change_nrdiis,&
                                             rhodiis
  USE spin,                            ONLY: clsd
  USE system,                          ONLY: cntl,&
                                             fpar
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mixing_r

CONTAINS

  ! ==================================================================
  SUBROUTINE mixing_r(nfr,gemax,rin0,rout0,rmix,thl)
    ! ==--------------------------------------------------------------==
    ! == Routine to perform density mixing.                           ==
    ! == RIN0 is input density in real space                          ==
    ! == ROUT0 is output density in real space                        ==
    ! == RMIX is mixed density  in real space                         ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nfr
    REAL(real_8) :: gemax, rin0(fpar%nnr1,clsd%nlsd), &
      rout0(fpar%nnr1,clsd%nlsd), rmix(fpar%nnr1,clsd%nlsd), thl(2)

    CHARACTER(*), PARAMETER                  :: procedureN = 'mixing_r'

    INTEGER                                  :: ierr, nnx
    INTEGER, SAVE                            :: ifirst = 0
    REAL(real_8), ALLOCATABLE                :: dm(:,:), vb(:)
    REAL(real_8), ALLOCATABLE, SAVE          :: dmat(:,:,:)
    REAL(real_8), ALLOCATABLE, SAVE, TARGET  :: driter(:,:), riter(:,:)
    REAL(real_8), POINTER, SAVE              :: rinm1(:,:), routm1(:,:)

! ==--------------------------------------------------------------==
! Test SCR

    ALLOCATE(vb(andr3%nrdiismax+1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(dm(andr3%nrdiismax+1, andr3%nrdiismax+1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! Allocation of SCR
    ! ==--------------------------------------------------------------==
    nnx=fpar%nnr1*clsd%nlsd
    IF (ifirst.EQ.0) THEN
       ifirst=1
       IF (andr3%nrdiismax.GT.0) THEN
          ALLOCATE(riter(andr3%nrdiismax*fpar%nnr1,clsd%nlsd),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(driter(andr3%nrdiismax*fpar%nnr1,clsd%nlsd),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(dmat(andr3%nrdiismax,andr3%nrdiismax,clsd%nlsd),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ! We can use Anderson mixing and after cntl%diis mixing.
          rinm1 => riter
          routm1 => driter
       ELSE
          ALLOCATE(rinm1(fpar%nnr1,clsd%nlsd),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)

          ALLOCATE(routm1(fpar%nnr1,clsd%nlsd),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
       CALL zeroing(rinm1)!,nnx)
       CALL zeroing(routm1)!,nnx)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Initialization for density mixing (needed in cntl%md).
    IF (nfr.EQ.1) THEN
       IF (andr3%nrdiismax.GT.0) THEN
          andr3%inrdiis=1
          CALL rhodiis(0,rmix(1,1),rin0(1,1),rout0(1,1),riter(1,1),&
               driter(1,1),andr2%andrmix,dmat(1,1,1),dm,vb,andr3%nrdiismax,0,1,thl(1))
          IF (cntl%tlsd) THEN
             CALL rhodiis(0,rmix(1,2),rin0(1,2),rout0(1,2),riter(1,2),&
                  driter(1,2),andr2%andrmix,dmat(1,1,2),dm,vb,andr3%nrdiismax,0,2,thl(1))
          ENDIF
       ENDIF
    ENDIF
    ! Change parameter if needed.
    CALL change_xmix(andr2%andrmix,gemax)
    IF (andr3%nrdiismax.GT.0) CALL change_nrdiis(gemax)
    IF (cntl%tlsd) THEN
       IF (andr3%nrdiis.GT.0) THEN
          CALL rhodiis(fpar%nnr1,rmix(1,1),rin0(1,1),rout0(1,1),riter(1,1),&
               driter(1,1),andr2%andrmix,dmat(1,1,1),dm,vb,andr3%nrdiismax,andr3%nrdiis,1,thl(&
               1))
          CALL rhodiis(fpar%nnr1,rmix(1,2),rin0(1,2),rout0(1,2),riter(1,2),&
               driter(1,2),andr2%andrmix,dmat(1,1,2),dm,vb,andr3%nrdiismax,andr3%nrdiis,2,thl(&
               2))
       ELSE
          CALL anderson(andr2%andrmix,nfr,rinm1(:,1),routm1(:,1),rin0(:,1),&
               rout0(:,1),rmix(:,1),fpar%nnr1,thl(1))
          CALL anderson(andr2%andrmix,nfr,rinm1(:,2),routm1(:,2),rin0(:,2),&
               rout0(:,2),rmix(:,2),fpar%nnr1,thl(2))
          CALL dcopy(nnx,rin0(1,1),1,rinm1(1,1),1)
          CALL dcopy(nnx,rout0(1,1),1,routm1(1,1),1)
       ENDIF
    ELSE
       IF (andr3%nrdiis.GT.0) THEN
          CALL rhodiis(fpar%nnr1,rmix,rin0,rout0,riter,driter,andr2%andrmix,dmat,&
               dm,vb,andr3%nrdiismax,andr3%nrdiis,1,thl(1))
          ! THL(1)=0._real_8
       ELSE
          CALL anderson(andr2%andrmix,nfr,rinm1,routm1,rin0,rout0,rmix,fpar%nnr1,&
               thl(1))
          CALL dcopy(nnx,rin0(1,1),1,rinm1(1,1),1)
          CALL dcopy(nnx,rout0(1,1),1,routm1(1,1),1)
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    DEALLOCATE(vb,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(dm,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE mixing_r
  ! ==================================================================


END MODULE mixing_r_utils
