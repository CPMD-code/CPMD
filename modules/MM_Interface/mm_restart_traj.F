C     ==================================================================
      SUBROUTINE mm_restart_traj(TAU,VEL,IREC)
      use kinds, only: real_4, real_8, int_1, int_2, int_4, int_8
      use error_handling, only: stopgm
C     ==--------------------------------------------------------------==
cmb - Revised on 23 Feb. 2005
      use system
      use mp_interface, only: mp_bcast
      use parac
      use ions
      use store_types
      use mm_input
      use mm_dim_utils, only : mm_dim
      use mm_dimmod
      use fileopenmod
      use fileopen_utils, only : fileclose,fileopen

      IMPLICIT NONE
      real*8 TAU(3,maxsys%nax,*)
      real*8 VEL(3,maxsys%nax,*)
      REAL*8, ALLOCATABLE :: GR_TAU(:,:),GR_VEL(:,:)
      INTEGER, ALLOCATABLE :: GR_IAT(:)
      integer IREC(100)
!
      integer i,j,iframe,ix
      integer MSGLEN,ierr
      LOGICAL FERROR,OLDSTATUS
      character(*),parameter :: procedureN='mm_restart_traj'
!
      CALL mm_dim(mm_go_mm,OLDSTATUS)
      allocate(GR_TAU(3,ions1%nat),STAT=ierr)
      if(ierr/=0) call stopgm(procedureN,'allocation problem',
     & __LINE__,'__file__')
      allocate(GR_VEL(3,ions1%nat),STAT=ierr)
      if(ierr/=0) call stopgm(procedureN,'allocation problem',
     & __LINE__,'__file__')
      allocate(GR_IAT((ions1%nat+1)/2),STAT=ierr)
      if(ierr/=0) call stopgm(procedureN,'allocation problem',
     & __LINE__,'__file__')

      IF(paral%io_parent) THEN
        write(6,*) 'RESTARTING FROM TRAJECTORY FILE ',
     &        rtr_c%file_traj_restart
        write(6,*) 'FRAME',rtr_i%iframe_restart
        CALL FILEOPEN(15,rtr_c%file_traj_restart,FO_OLD,FERROR)
        IF(FERROR) call stopgm('MM_RESTART_TRAJ','TRAJECTORY NOT FOUND',
     & __LINE__,"__file__")

        DO i=1,ions1%nat
C     search to requested frame           
           IF(i.EQ.1) THEN
              DO j=1,2147483647-1 ! max value for integer*4 -1 for the last iteration
                 READ(15,*,end=145,err=146) iframe,
     &                (GR_TAU(ix,i),ix=1,3),(GR_VEL(ix,i),ix=1,3)
                 IF(iframe.eq.rtr_i%iframe_restart) GOTO 140
              END DO
           ELSE
              READ(15,*,end=145,err=146) iframe,
     &             (GR_TAU(ix,i),ix=1,3),(GR_VEL(ix,i),ix=1,3)
           ENDIF
 140       CONTINUE
C     convert from gromos to internal coordinate ordering
           DO j=1,3
              TAU(j,cpat(i),cpsp(i))=GR_TAU(j,i)
              VEL(j,cpat(i),cpsp(i))=rtr_r%rev_vel*GR_VEL(j,i)
           END DO
           IF(iframe.ne.rtr_i%iframe_restart) THEN
              WRITE(6,*) 'CORRUPTED OR WRONG TRAJECTORY FILE !'
              GOTO 146
           ENDIF
        END DO
        CALL FILECLOSE(15)
        CALL FILEOPEN(15,'MM_CELL_TRANS',FO_OLD,FERROR)
        IF(FERROR) THEN
           WRITE(6,*) 'MM_RESTART_TRAJ| MM_CELL_TRANS FILE NOT FOUND'
           clsaabox%MM_C_TRANS(1)=0.0d0
           clsaabox%MM_C_TRANS(2)=0.0d0
           clsaabox%MM_C_TRANS(3)=0.0d0
        ELSE
           DO j=1,2147483647-1 ! max value for integer*4 -1 for the last iteration
              READ(15,*,end=145,err=146) iframe,
     &             (clsaabox%MM_C_TRANS(ix),ix=1,3)
              IF(iframe.eq.rtr_i%iframe_restart) GOTO 150
           END DO
 150       CONTINUE
           CALL FILECLOSE(15)
        ENDIF
      ENDIF                     ! parent
      
      IREC(IREC_VEL)=1
      CALL Mp_BCAST(TAU,3*maxsys%nax*maxsys%nsx,parai%io_source,
     &     parai%cp_grp)
      CALL Mp_BCAST(VEL,3*maxsys%nax*maxsys%nsx,parai%io_source,
     &     parai%cp_grp)
      CALL mm_center(TAU,clsaabox%mm_c_trans,.TRUE.)

      CALL mm_dim(mm_revert,OLDSTATUS)
      deallocate(GR_TAU,STAT=ierr)
      if(ierr/=0) call stopgm(procedureN,'deallocation problem',
     & __LINE__,'__file__')
      deallocate(GR_VEL,STAT=ierr)
      if(ierr/=0) call stopgm(procedureN,'deallocation problem',
     & __LINE__,'__file__')
      deallocate(GR_IAT,STAT=ierr)
      if(ierr/=0) call stopgm(procedureN,'deallocation problem',
     & __LINE__,'__file__')

      RETURN 
      
 145  call stopgm('MM_RESTART_TRAJ','FRAME NOT FOUND IN TRAJECTORY',
     & __LINE__,"__file__")
 146  call stopgm('MM_RESTART_TRAJ','ERROR READING TRAJECTORY',
     & __LINE__,"__file__")
      END
