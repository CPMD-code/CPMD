MODULE pi_wf_utils
  USE coor,                            ONLY: tau0
  USE error_handling,                  ONLY: stopgm
  USE filnmod,                         ONLY: filbod,&
                                             filn
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_bcast
  USE parac,                           ONLY: paral,&
                                             parai
  USE pimd,                            ONLY: ipcurr,&
                                             np_high,&
                                             np_low,&
                                             pimd1,&
                                             pimd3,&
                                             trep,&
                                             trepnm
  USE readsr_utils,                    ONLY: xstring
  USE repgen_utils,                    ONLY: repgen
  USE rreadf_utils,                    ONLY: rreadf
  USE store_types,                     ONLY: restart1
  USE system,                          ONLY: cntl,&
                                             maxsys
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE wfopts_utils,                    ONLY: wfopts

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: pi_wf
  PUBLIC :: initrep

CONTAINS

  ! ==================================================================
  SUBROUTINE pi_wf
    ! ==--------------------------------------------------------------==
    ! ==    CONTROL WAVEFUNCTION OPTIMIZATION FOR PATH INTEGRAL RUNS  ==
    ! ==--------------------------------------------------------------==
    ! Variables
    CHARACTER(*), PARAMETER                  :: procedureN = 'pi_wf'

    CHARACTER(len=12)                        :: cflbod, cipnum
    INTEGER                                  :: i1, i2, ierr, ip, isub, n1, n2

! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    IF ( cntl%tqmmm ) CALL stopgm("PI_WF","QMMM NOT IMPLEMENTED",& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    IF (pimd1%tinit) THEN
       ! ..Generate Replica coordinates
       ALLOCATE(trep(3,maxsys%nax,maxsys%nsx,pimd3%np_total),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       IF (pimd1%tread_cent.OR.pimd1%tpinm) THEN
          ALLOCATE(trepnm(3,maxsys%nax,maxsys%nsx,pimd3%np_total),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
       IF (pimd1%rrep) THEN
          CALL repgen(tau0)
          ! DM1
          ! ..Classical test 
          IF (pimd1%testpi) THEN
             CALL dcopy(3*maxsys%nax*maxsys%nsx,tau0,1,trep(1,1,1,1),1)
          ENDIF
          ! DM2
       ELSEIF (pimd1%repread) THEN
          CALL rreadf
       ENDIF
    ENDIF
    ! ..optimize wavefunction for Trotter index IPCURR
    DO ip=np_low,np_high
       ! IPCUUR is in the include file pimd.inc.
       ipcurr=ip
       IF (pimd1%tinit) THEN
          CALL dcopy(3*maxsys%nax*maxsys%nsx,trep(1,1,1,ipcurr),1,tau0,1)
       ENDIF
       ! ..Construct filenames
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
       IF (paral%parent) THEN
          IF (paral%io_parent)&
               WRITE(6,'(/,1X,64("*"))')
          IF (paral%io_parent)&
               WRITE(6,'(" *",62X,"*")')
          IF (paral%io_parent)&
               WRITE(6,'(" *",23X,A,I5,24X,"*")') ' REPLICA :',ipcurr
          IF (paral%io_parent)&
               WRITE(6,'(" *",62X,"*")')
          IF (paral%io_parent)&
               WRITE(6,'(1X,64("*"),/)')
       ENDIF
       CALL wfopts
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
    CALL tihalt(procedureN,isub)
  END SUBROUTINE pi_wf
  ! ==================================================================
  SUBROUTINE initrep(pitaup)
    ! ==--------------------------------------------------------------==
    ! == initialize replica coordinates in the initialization process ==
    ! == of path integral MD runs                                     ==
    ! ==--------------------------------------------------------------==
    ! Variables
    REAL(real_8)                             :: pitaup(:,:,:,:)
    CHARACTER(*), PARAMETER                  :: procedureN = 'initrep'
    INTEGER                                  :: ierr
    ! ==--------------------------------------------------------------==
    ! ..Generate Replica coordinates
    ALLOCATE(trep(3,maxsys%nax,maxsys%nsx,pimd3%np_total),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    IF (pimd1%tread_cent.OR.pimd1%tpinm) THEN
       ALLOCATE(trepnm(3,maxsys%nax,maxsys%nsx,pimd3%np_total),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    IF (pimd1%rrep) THEN
       CALL repgen(tau0)
       ! DM1
       ! ..Classical test
       IF (pimd1%testpi) THEN
          CALL dcopy(3*maxsys%nax*maxsys%nsx,tau0,1,trep(1,1,1,1),1)
       ENDIF
       ! DM2
    ELSEIF (pimd1%repread) THEN
       CALL rreadf
    ENDIF
    ! ..Copy trep to pitaup and broadcast in cp_grp 
    DO ipcurr=np_low,np_high
       IF (paral%io_parent) THEN
          CALL dcopy(3*maxsys%nax*maxsys%nsx,trep(1,1,1,ipcurr),1,&
             pitaup(1,1,1,ipcurr),1)
       ENDIF
       CALL mp_bcast(pitaup(:,:,:,ipcurr),3*maxsys%nax*maxsys%nsx,&
          parai%io_source,parai%cp_grp)
    ENDDO
    DEALLOCATE(trep,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (pimd1%tread_cent.OR.pimd1%tpinm) THEN
       DEALLOCATE(trepnm,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE initrep
  ! ==================================================================

END MODULE pi_wf_utils
