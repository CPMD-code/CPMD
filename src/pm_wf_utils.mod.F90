MODULE pm_wf_utils
  USE coor,                            ONLY: tau0
  USE error_handling,                  ONLY: stopgm
  USE filnmod,                         ONLY: filbod,&
                                             filn
  USE parac,                           ONLY: paral
  USE pimd,                            ONLY: ipcurr,&
                                             np_high,&
                                             np_low,&
                                             pimd1,&
                                             pimd3,&
                                             trep,&
                                             trepnm
  USE readsr_utils,                    ONLY: xstring
  USE rreadf_utils,                    ONLY: rreadf
  USE store_types,                     ONLY: restart1
  USE system,                          ONLY: cntl,&
                                             maxsys
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE wfopts_utils,                    ONLY: wfopts

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: pm_wf

CONTAINS

  ! ==================================================================
  SUBROUTINE pm_wf
    ! ==--------------------------------------------------------------==
    ! ==    CONTROL WAVEFUNCTION OPTIMIZATION FOR PATH MINIMISATION   ==
    ! ==--------------------------------------------------------------==
    ! Variables
    CHARACTER(*), PARAMETER                  :: procedureN = 'pm_wf'

    CHARACTER(len=12)                        :: cflbod, cipnum
    INTEGER                                  :: i1, i2, ierr, ip, isub, n1, n2

    CALL tiset(procedureN,isub)

    ! ==--------------------------------------------------------------==
    IF (cntl%tqmmm) THEN
       CALL stopgm('PM_WF','QMMM NOT IMPLEMENTED',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (pimd1%tinit) THEN
       ! Generate Replica coordinates
       ALLOCATE(trep(3,maxsys%nax,maxsys%nsx,pimd3%np_total),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       IF (pimd1%tread_cent.OR.pimd1%tpinm) THEN
          ALLOCATE(trepnm(3,maxsys%nax,maxsys%nsx,pimd3%np_total),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
       IF (pimd1%rrep) THEN
          CALL stopgm('PM_WF','GENERATE REPLICA NOT IMPLEMENTED',& 
               __LINE__,__FILE__)
          ! Classical test 
          IF (pimd1%testpi) THEN
             CALL dcopy(3*maxsys%nax*maxsys%nsx,tau0,1,trep(1,1,1,1),1)
          ENDIF
       ELSEIF (pimd1%repread) THEN
          CALL rreadf
       ENDIF
    ENDIF
    ! Optimize wavefunction for replica index IPCURR.
    DO ip=np_low,np_high
       ! IPCUUR is in the include file pimd.inc.
       ipcurr=ip
       IF (pimd1%tinit) THEN
          CALL dcopy(3*maxsys%nax*maxsys%nsx,trep(1,1,1,ipcurr),1,tau0,1)
       ENDIF
       ! Construct filenames
       cflbod='RESTART_'
       IF (paral%io_parent)&
            WRITE(cipnum,'(I4)') ipcurr
       CALL xstring(cflbod,n1,n2)
       CALL xstring(cipnum,i1,i2)
       filbod=cflbod(n1:n2)//cipnum(i1:i2)//'.1'
       IF (restart1%rlate) THEN
          filn=filbod
       ELSE
          filn=cflbod(n1:n2)//cipnum(i1:i2)
       ENDIF
       IF (paral%io_parent) THEN
          WRITE(6,'(/,1X,64("*"))')
          WRITE(6,'(" *",62X,"*")')
          WRITE(6,'(" *",23X,A,I5,24X,"*")') ' REPLICA :',ipcurr
          WRITE(6,'(" *",62X,"*")')
          WRITE(6,'(1X,64("*"),/)')
       ENDIF
       CALL wfopts
    ENDDO
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    RETURN
  END SUBROUTINE pm_wf
  ! ==================================================================

END MODULE pm_wf_utils
