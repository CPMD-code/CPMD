MODULE pm_init_utils
  USE cotr,                            ONLY: cotr007
  USE envj,                            ONLY: tjlimit
  USE error_handling,                  ONLY: stopgm
  USE mfep,                            ONLY: irest,&
                                             lcvsp,&
                                             mfepi
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sync
  USE mp_multiple_comm_init,           ONLY: multiple_comm_init
  USE parac,                           ONLY: parai,&
                                             paral
  USE pimd,                            ONLY: grandparent,&
                                             parentgroup,&
                                             supergroup,&
                                             supersource
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE wann,                            ONLY: wan05,&
                                             wannl

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: pm_init

CONTAINS

  ! ==================================================================
  SUBROUTINE pm_init
    ! ==--------------------------------------------------------------==
    ! ==    INITIALIZE PATH MINIMISATION RUNS                         ==
    ! ==--------------------------------------------------------------==

    ! Variables
    CHARACTER(*), PARAMETER                  :: procedureN = 'pm_init'

    INTEGER                                  :: i, ierr, isub

    CALL tiset(procedureN,isub)
    CALL multiple_comm_init("REPLICA")

    ! Broadcast maximum time for the job
    IF (paral%io_parent) CALL mp_bcast(tjlimit,supersource,parentgroup)

    ! List and flags for collective variables 
    !IF (paral%io_parent) THEN
    IF (mfepi%ncvsp.GE.1.AND.mfepi%ncvsp.LE.cotr007%mrestr) THEN
       IF (.NOT.grandparent) THEN
          ALLOCATE(irest(mfepi%ncvsp),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
       !CALL mp_bcast(irest,SIZE(irest),supersource,parentgroup)
       CALL mp_bcast(irest,SIZE(irest),supersource,supergroup)
    ELSE
       CALL stopgm('PM_INIT',' WRONG NUMBER OF NCVSP ',& 
            __LINE__,__FILE__)
    ENDIF
    ALLOCATE(lcvsp(cotr007%mrestr),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    DO i=1,cotr007%mrestr
       lcvsp(i)=.FALSE.
    ENDDO
    DO i=1,mfepi%ncvsp
       lcvsp(irest(i))=.TRUE.
    ENDDO

    IF (wannl%twann) THEN
       IF (wan05%loc_npgrp>parai%cp_nproc) wan05%loc_npgrp=parai%cp_nproc
    ENDIF

    CALL mp_sync(supergroup)
    CALL tihalt(procedureN,isub)
    RETURN
  END SUBROUTINE pm_init
  ! ==================================================================

END MODULE pm_init_utils
