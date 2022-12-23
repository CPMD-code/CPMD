MODULE fnlalloc_utils
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions1
  USE kinds,                           ONLY: real_8
  USE nlps,                            ONLY: imagp,&
                                             ndfnl
  USE parac,                           ONLY: parai
  USE sfac,                            ONLY: dfnl,&
                                             fnl,&
                                             fnl2,&
                                             ldf1,&
                                             ldf2,&
                                             tfnl2
  USE system,                          ONLY: cntl,&
                                             maxsys,&
                                             natpe,&
                                             nkpt,&
                                             norbpe,&
                                             parap
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: fnlalloc
  PUBLIC :: fnldealloc
  PUBLIC :: fnl_set

CONTAINS

  ! ==================================================================
  SUBROUTINE fnlalloc(nstate,tfor,tstress)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    LOGICAL                                  :: tfor, tstress

    CHARACTER(*), PARAMETER                  :: procedureN = 'fnlalloc'

    INTEGER                                  :: ierr, ldfnl, lfnl

    IF (cntl%tfdist) THEN
       tfnl2=.TRUE.

       lfnl=MAX(1,imagp*nstate*natpe*maxsys%nhxs*nkpt%nkpnt)
       IF (lfnl == 1) THEN
          ALLOCATE(fnl(1,1,1,1,1),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ELSE
          ALLOCATE(fnl(imagp,natpe,maxsys%nhxs,nstate,nkpt%nkpnt),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF

       CALL zeroing(fnl)!,lfnl)

       lfnl=MAX(1,imagp*norbpe*ions1%nat*maxsys%nhxs*nkpt%nkpnt)
       IF (lfnl == 1) THEN
          ALLOCATE(fnl2(1,1,1,1,1),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ELSE
          ALLOCATE(fnl2(imagp,ions1%nat,maxsys%nhxs,norbpe,nkpt%nkpnt),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF

       CALL zeroing(fnl2)!,lfnl)
       ldf1=natpe
       ldf2=norbpe

    ELSE

       tfnl2=.FALSE.
       lfnl=MAX(1,imagp*nstate*ions1%nat*maxsys%nhxs*nkpt%nkpnt)
       IF (lfnl == 1) THEN
          ALLOCATE(fnl(1,1,1,1,1),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ELSE
          ALLOCATE(fnl(imagp,ions1%nat,maxsys%nhxs,nstate,nkpt%nkpnt),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF

       CALL zeroing(fnl)!,lfnl)

       fnl2 => fnl

       ldf1=ions1%nat
       ldf2=nstate

    ENDIF
    ! DFNL
    IF ( cntl%tshop ) THEN
       ndfnl=INT(nstate/parai%nproc)+1
    ELSE
       ndfnl=parap%nst12(parai%mepos,2)-parap%nst12(parai%mepos,1)+1
    ENDIF
    ldfnl=imagp*3*ions1%nat*maxsys%nhxs*ndfnl*nkpt%nkpnt
    IF (ldfnl.LE.0) ldfnl=1
    IF (ldfnl == 1) THEN
       ALLOCATE(dfnl(1,1,1,1,1,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ELSE
       ALLOCATE(dfnl(imagp,ions1%nat,maxsys%nhxs,3,ndfnl,nkpt%nkpnt),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE fnlalloc
  ! ==================================================================
  SUBROUTINE fnldealloc(tfor,tstress)
    ! ==--------------------------------------------------------------==
    LOGICAL                                  :: tfor, tstress

    CHARACTER(*), PARAMETER                  :: procedureN = 'fnldealloc'

    INTEGER                                  :: ierr

! ==--------------------------------------------------------------==

    DEALLOCATE(fnl,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (tfnl2) DEALLOCATE(fnl2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(dfnl,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE fnldealloc
  ! ==================================================================
  SUBROUTINE fnl_set(tag)
    ! ==--------------------------------------------------------------==
    ! C.BEKAS

    CHARACTER(len=*)                         :: tag

    INTEGER                                  :: ldfd
    INTEGER, SAVE                            :: ldf1b, ldf2b, ndfnlb
    LOGICAL, SAVE                            :: tfnl2b, tfnl2d
    REAL(real_8), POINTER                    :: xdum(:,:,:,:,:), &
                                                xdum2(:,:,:,:,:,:)
    REAL(real_8), POINTER, SAVE              :: xdfnl(:,:,:,:,:,:), &
                                                xfnl(:,:,:,:,:), &
                                                xfnl2(:,:,:,:,:)

    IF (INDEX(tag,'SAVE').NE.0) THEN

       xfnl => fnl
       xfnl2 => fnl2
       xdfnl => dfnl

       ldf1b=ldf1
       ldf2b=ldf2
       ndfnlb=ndfnl
       tfnl2b=tfnl2
    ELSEIF (INDEX(tag,'RECV').NE.0) THEN
       fnl => xfnl
       fnl2 => xfnl2
       dfnl => xdfnl

       ldf1=ldf1b
       ldf2=ldf2b
       ndfnl=ndfnlb
       tfnl2=tfnl2b
    ELSEIF (INDEX(tag,'SWITCH').NE.0) THEN
       xdum => fnl
       fnl => xfnl
       xfnl => xdum

       xdum => fnl2
       fnl2 => xfnl2
       xfnl2 => xdum

       xdum2 => dfnl
       dfnl => xdfnl
       xdfnl => xdum2

       ldfd=ldf1
       ldf1=ldf1b
       ldf1b=ldfd
       ldfd=ldf2
       ldf2=ldf2b
       ldf2b=ldfd
       ldfd=ndfnl
       ndfnl=ndfnlb
       ndfnlb=ldfd
       tfnl2d=tfnl2
       tfnl2=tfnl2b
       tfnl2b=tfnl2d
    ELSEIF (INDEX(tag,'MIX').NE.0) THEN
       xdum2 => dfnl
       dfnl => xdfnl
       xdfnl => xdum2

       ldfd=ndfnl
       ndfnl=ndfnlb
       ndfnlb=ldfd
    ELSE
       CALL stopgm('FNL_SET','INVALID TAG',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE fnl_set
  ! ==================================================================

END MODULE fnlalloc_utils
