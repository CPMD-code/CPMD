MODULE mp_multiple_comm_init
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_def,&
                                             fo_scratch
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_environ,&
                                             mp_group,&
                                             mp_max,&
                                             mp_split
  USE parac,                           ONLY: parai,&
                                             paral
  USE pimd,                            ONLY: &
       grandparent, ipm1, ipp1, np_high, np_local, np_low, nproc_tot, num_np, &
       parentgroup, pc_groups, pc_grp, pcg_pos, pimd3, supergroup, supersource
  USE readsr_utils,                    ONLY: xstring
  USE set_cp_grp_utils,                ONLY: reset_cp_grp
  USE system,                          ONLY: cntr,&
                                             parap
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: multiple_comm_init

CONTAINS

  ! ==================================================================
  SUBROUTINE multiple_comm_init(label)
    ! ==--------------------------------------------------------------==
    ! ==    INITIALIZE MULTIPLE COMMUNICATORS
    ! ==--------------------------------------------------------------==

    CHARACTER(LEN=*), INTENT(IN)             :: label

    CHARACTER(*), PARAMETER :: procedureN = 'multiple_comm_init'

    CHARACTER(len=12)                        :: fileon, fileout
    INTEGER :: COLOR, i, i1, i2, ip, ipp, isub, j, me_tmp, mhigh, mlow, &
      my_sg_me, my_sg_nproc, nhigh, nlow, nplocal, npp, pg_me, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: my_pc_grp
    LOGICAL                                  :: ferror
    REAL(real_8)                             :: rproc

    CALL tiset(procedureN,isub)
    CALL mp_bcast(cntr%ecut,parai%io_source,parai%cp_grp) ! TODO is it needed?

    ! CP_GRP becomes SUPERGROUP
    grandparent=paral%io_parent
    supergroup=parai%cp_grp
    supersource=parai%io_source
    nproc_tot=parai%cp_nproc
    pg_me=parai%cp_me
    ! ..Number of processor groups
    IF (paral%io_parent.AND.pc_groups.EQ.0) THEN
       DO i=1,parai%cp_nproc
          IF (MOD(parai%cp_nproc,i).EQ.0) THEN
             npp=parai%cp_nproc/i
             IF (MOD(pimd3%np_total,npp).EQ.0) THEN
                pc_groups=npp
                GOTO 100
             ENDIF
          ENDIF
       ENDDO
       pc_groups=1
100    CONTINUE
    ENDIF
    CALL mp_bcast(pc_groups,parai%io_source,parai%cp_grp)

    pc_grp(1)=0
    pcg_pos=1
    ! ..generate processor groups
    CALL zeroing(parap%pgroup)!,maxcpu)
    CALL zeroing(parap%nlink)!,maxcpu)
    rproc=REAL(parai%cp_nproc,kind=real_8)/pc_groups
    DO i=1,pc_groups
       nlow=NINT((i-1)*rproc)
       pc_grp(i)=nlow
       nhigh=NINT(i*rproc)-1
       IF (parai%cp_me.GE.nlow .AND. parai%cp_me.LE.nhigh) THEN
          parai%cp_nproc=nhigh-nlow+1
          parai%mepos=parai%cp_me-nlow
          pcg_pos=i
          DO ip=1,parai%cp_nproc
             parap%pgroup(ip)=ip-1
             ipp=parap%pgroup(ip)
             parap%nlink(ipp)=ip-1
          ENDDO
       ENDIF
    ENDDO


    color=parai%cp_me/parai%cp_nproc
    CALL mp_split(supergroup,color,parai%cp_me,parai%cp_grp,me_tmp)
    IF(parai%mepos.NE.me_tmp) THEN
       CALL stopgm(procedureN,'!MY_SPLIT RETURNED ERROR!',& 
            __LINE__,__FILE__)
    END IF
    parai%cp_me=me_tmp


    !     Reset the cp group communicator
    !
    CALL reset_cp_grp()


    !vw>>>
    !vw !!! need to set the parentgroup AFTER the cp_grp are set !!!
    !vw !!! this is due to the potential reordering of the tasks by the underlying mpi library !!!
    !vw not needed CALL mp_sync(supergroup)
    CALL mp_environ(supergroup ,my_sg_nproc,my_sg_me)

    ALLOCATE(my_pc_grp(0:my_sg_nproc-1),stat=stat)
    IF (stat.NE.0) CALL stopgm(procedureN,'Allocation problem',&
         __LINE__,__FILE__)

    my_pc_grp(:) = -1
    IF( paral%io_parent ) my_pc_grp(my_sg_me) = my_sg_me
    CALL mp_max(my_pc_grp,SIZE(my_pc_grp),supergroup)
    j = 0
    DO i = 0, my_sg_nproc-1
       IF( my_pc_grp(i) /= -1 ) THEN
          j = j + 1; pc_grp(j) = my_pc_grp(i)
       ENDIF
    ENDDO
    IF( j /= pc_groups ) CALL stopgm(procedureN,'Something wrong here',&
         __LINE__,__FILE__)

    CALL mp_group(pc_groups,pc_grp,parentgroup,supergroup)

    DEALLOCATE(my_pc_grp,stat=stat)
    IF (stat.NE.0) CALL stopgm(procedureN,'Deallocation problem',&
         __LINE__,__FILE__)
    !vw <<<


    IF (pimd3%loutfn.EQ.0) THEN
       IF (.NOT.grandparent) THEN
#if defined(__WINNT)
          WRITE(fileon,'(I3)') pcg_pos
          CALL xstring(fileon,i1,i2)
          fileout='OUTPUT_'//fileon(i1:i2)
          CALL fileclose(6)
          CALL fileopen(6,fileout,FO_DEF,FERROR)
#else
          CALL fileclose(6)
          CALL fileopen(6,'',FO_SCRATCH,FERROR)
#if defined(__Linux) || defined (__ALPHALINUX)
          ! AK: we need to do the same with the c-style stdout for printmemsize
          CALL silentstdout
#endif
#endif
       ENDIF
    ELSEIF (pimd3%loutfn.EQ.1) THEN
       IF (paral%io_parent.AND..NOT.grandparent) THEN
          WRITE(fileon,'(I3)') pcg_pos
          CALL xstring(fileon,i1,i2)
          fileout='OUTPUT_'//fileon(i1:i2)
          CALL fileclose(6)
          CALL fileopen(6,fileout,FO_DEF,FERROR)
#if defined(__Linux) || defined (__ALPHALINUX)
          ! avoid printmemsize messes
          CALL silentstdout
#endif
       ELSEIF (.NOT.grandparent) THEN
#if defined(__WINNT)
          WRITE(fileon,'(I3)') pcg_pos
          CALL xstring(fileon,i1,i2)
          fileout='OUTPUT_'//fileon(i1:i2)
          CALL fileclose(6)
          CALL fileopen(6,fileout,FO_DEF,FERROR)
#else
          CALL fileclose(6)
          CALL fileopen(6,'',FO_SCRATCH,FERROR)
#if defined(__Linux) || defined (__ALPHALINUX)
          ! AK: we need to do the same with the c-style stdout for printmemsize
          CALL silentstdout
#endif
#endif
       ENDIF
    ELSEIF (pimd3%loutfn.EQ.2) THEN
       IF (.NOT.grandparent) THEN
          WRITE(fileon,'(I3)') pg_me
          CALL xstring(fileon,i1,i2)
          fileout='OUTPUT_'//fileon(i1:i2)
          CALL fileclose(6)
          CALL fileopen(6,fileout,FO_DEF,FERROR)
#if defined(__Linux) || defined (__ALPHALINUX)
          ! avoid printmemsize messes
          CALL silentstdout
#endif
       ENDIF
    ENDIF
    ! ..establish connection Replica/Trotter index with pc-group
    rproc=REAL(pimd3%np_total,kind=real_8)/pc_groups
    DO i=1,pc_groups
       nlow=NINT((i-1)*rproc) + 1
       nhigh=NINT(i*rproc)
       nplocal=nhigh-nlow+1
       DO j=nlow,nhigh
          num_np(j)=i
       ENDDO
       IF (pcg_pos.EQ.i) THEN
          np_low=nlow
          np_high=nhigh
          np_local=nplocal
       ENDIF
    ENDDO
    ! ..Connectivity of replica
    DO ip=1,pimd3%np_total
       ipp1(ip)=ip+1
       ipm1(ip)=ip-1
    ENDDO
    ipp1(pimd3%np_total)=1
    ipm1(1)=pimd3%np_total
    ! ..Print some info
    IF (grandparent) THEN
       IF (paral%io_parent)&
            WRITE(6,'(8X,A,A)') 'PC_GROUP',&
            '        PROCESSORS        '//TRIM(label)//' INDICES'
       rproc=REAL(pimd3%np_total,kind=real_8)/pc_groups
       DO i=1,pc_groups
          nlow=NINT((i-1)*rproc) + 1
          nhigh=NINT(i*rproc)
          mlow=pc_grp(i)
          IF (i.EQ.pc_groups) THEN
             mhigh=nproc_tot-1
          ELSE
             mhigh=pc_grp(i+1)-1
          ENDIF
          IF (paral%io_parent)&
               WRITE(6,'(9X,I4,10X,I4,A,I4,12X,I3,A,I3)') i,&
               mlow,' - ',mhigh,nlow,' - ',nhigh
       ENDDO
       IF (paral%io_parent)&
            WRITE(6,'(A)') 
    ENDIF

    CALL tihalt(procedureN,isub)

  END SUBROUTINE multiple_comm_init
  ! ==================================================================

END MODULE mp_multiple_comm_init
