#define PRINT_GROUP_INFOS .FALSE.

MODULE set_cp_grp_utils
  USE cp_grp_utils,                    ONLY: cp_grp_get_cp_rank
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: long_string_length
  USE machine,                         ONLY: m_compiler_options,&
                                             m_compiler_version,&
                                             m_flush
  USE mp_interface,                    ONLY: mp_cart,&
                                             mp_comm_world,&
                                             mp_environ,&
                                             mp_get_library_version,&
                                             mp_get_version,&
                                             mp_max,&
                                             mp_sum,&
                                             mp_sync
  USE parac,                           ONLY: parai,&
                                             paral
  USE pimd,                            ONLY: grandparent

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: set_cp_grp
  PUBLIC :: finalize_cp_grp
  PUBLIC :: reset_cp_grp

CONTAINS

  ! ==================================================================
  SUBROUTINE set_cp_grp()
    ! ==--------------------------------------------------------------==
    CHARACTER(*), PARAMETER                  :: procedureN = 'set_cp_grp'

    CHARACTER(2*long_string_length)          :: compiler_options, &
                                                compiler_version
    CHARACTER(long_string_length)            :: libversion
    INTEGER                                  :: cp_inter_nproc, cp_npgrp, &
                                                ierr, iprc, subversion, &
                                                version
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: cp_nolist, cp_nplist

! ==--------------------------------------------------------------==
! 
! check that there is enough processes for the number of groups
! 

    cp_npgrp = 0
    IF (parai%cp_me.EQ.0.AND.parai%cp_nogrp.LE.0) CALL stopgm(proceduren,&
         'number of cp_groups is less than zero.',& 
         __LINE__,__FILE__)
    IF (parai%cp_nogrp.GT.parai%cp_nproc) THEN
       IF (parai%cp_me.EQ.0) WRITE(6,'(1x,a)') 'warning: the number of '&
            //'cp_groups is greater than the available number '&
            //'of processes.'
       IF (parai%cp_me.EQ.0) WRITE(6,'(1x,a)') 'warning: the number of '&
            //'cp_groups will be reset! '
       CALL m_flush(6)
       parai%cp_nogrp = 0
    ELSEIF (MOD(parai%cp_nproc,parai%cp_nogrp).NE.0) THEN
       IF (parai%cp_me.EQ.0) WRITE(6,'(1x,a)') 'warning: the number of '&
            //'cp_groups is not a divisor of the available '&
            //'number of processes.'
       IF (parai%cp_me.EQ.0) WRITE(6,'(1x,a)') 'warning: the number of '&
            //'cp_groups will be reset! '
       CALL m_flush(6)
       parai%cp_nogrp = 0
    ELSE
       cp_npgrp = parai%cp_nproc / parai%cp_nogrp
    ENDIF

    ! 
    ! set the cp group communicator
    ! 
    ALLOCATE(cp_nolist(parai%cp_nproc),cp_nplist(parai%cp_nproc),stat=ierr)
    IF (ierr.NE.0) CALL stopgm(proceduren,'allocation problem',& 
         __LINE__,__FILE__)
    CALL mp_cart(mp_comm_world,parai%cp_nogrp,cp_npgrp,cp_nolist,cp_nplist,&
         parai%cp_inter_grp,parai%allgrp)
    DEALLOCATE(cp_nolist,cp_nplist,stat=ierr)
    IF (ierr.NE.0) CALL stopgm(proceduren,'deallocation problem',& 
         __LINE__,__FILE__)
    CALL mp_environ(parai%cp_inter_grp,cp_inter_nproc,parai%cp_inter_me)
    CALL mp_environ(parai%allgrp,parai%nproc,parai%me)
    parai%mepos = parai%me ! set this usefull variable as well

    ! 
    ! set up the source process within the group (allgrp communicator) 
    ! 
    parai%source = 0                ! range from 0 to nproc - 1
    paral%parent = parai%me.EQ.parai%source

    ! 
    ! set up the io process (cp_grp communicator)
    ! 
    ! io_source = 0            ! should range from 0 to cp_nproc - 1
    ! ! and independant of source
    ! ! but more work woould be needed !
    ! >>> so we hack
    parai%io_source = 0
    IF (parai%me.EQ.parai%source.AND.parai%cp_inter_me.EQ.0) parai%io_source = parai%cp_me
    CALL mp_max(parai%io_source,parai%cp_grp)
    ! <<<
    paral%io_parent = parai%cp_me.EQ.parai%io_source

    ! 
    ! set the source process for the intergroup communicator
    ! 
    parai%cp_inter_io_source = 0
    IF (parai%cp_me.EQ.parai%io_source) parai%cp_inter_io_source = parai%cp_inter_me
    CALL mp_max(parai%cp_inter_io_source,parai%cp_grp)
    paral%cp_inter_io_parent = parai%cp_inter_me.EQ.parai%cp_inter_io_source

    ! 
    ! set mapping (me,cp_inter_me) to cp_me
    ! 
    IF (ALLOCATED(cp_grp_get_cp_rank)) CALL stopgm(proceduren,&
         'this should be already allocated',& 
         __LINE__,__FILE__)
    ALLOCATE(cp_grp_get_cp_rank(0:parai%nproc-1,0:cp_inter_nproc),&
         stat=ierr)
    IF (ierr.NE.0) CALL stopgm(proceduren,'allocation problem',& 
         __LINE__,__FILE__)
    cp_grp_get_cp_rank = 0
    cp_grp_get_cp_rank(parai%me,parai%cp_inter_me) = parai%cp_me
    CALL mp_sum(cp_grp_get_cp_rank,parai%cp_nproc,parai%cp_grp)

    ! 
    ! some tests
    ! 
    IF (cp_inter_nproc.NE.parai%cp_nogrp) CALL stopgm(proceduren,&
         'something wrong here',& 
         __LINE__,__FILE__)
    IF (cp_npgrp.NE.parai%nproc) CALL stopgm(proceduren,&
         'something wrong here',& 
         __LINE__,__FILE__)

    IF (paral%io_parent) THEN
       CALL m_compiler_version( compiler_version )
       CALL m_compiler_options( compiler_options )
       CALL mp_get_version( version, subversion )
       CALL mp_get_library_version( libversion )
       WRITE(6,'(1x,a,a)') 'compiler version: ',TRIM(compiler_version)
       WRITE(6,'(1x,a,a)') 'compiler options: ',TRIM(compiler_options)
       WRITE(6,'(1x,a,a)'      ) 'MPI version of the library: ',TRIM(libversion)
       WRITE(6,'(1x,a,i0,a,i0)') 'MPI version of the standard: ',version, '.', subversion
       WRITE(6,'(1x,2(a,i0),a)') 'cp_groups: we are using a ',&
            parai%cp_nogrp,' x ',cp_npgrp,' grid (groups x nprocs).'
    ENDIF

    IF (PRINT_GROUP_INFOS) THEN
       DO iprc=0,parai%cp_nproc-1
          IF (iprc.EQ.parai%cp_me) THEN
             WRITE(6,'(3(a,2(1x,i0)),3(a,l2))')&
                  ' cp_grp',parai%cp_nproc,parai%cp_me,&
                  ' allgrp',parai%nproc,parai%me,&
                  ' cp_inter_grp',cp_inter_nproc,parai%cp_inter_me,&
                  ' parent',paral%parent,' cp_inter_io_parent',paral%cp_inter_io_parent,&
                  ' io_parent',paral%io_parent
          ENDIF
          CALL m_flush(6)
          CALL mp_sync(parai%cp_grp)
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
  END SUBROUTINE set_cp_grp


  SUBROUTINE reset_cp_grp()
    !     ==--------------------------------------------------------------==
    CHARACTER(*), PARAMETER                  :: procedureN = 'reset_cp_grp'

    INTEGER                                  :: cp_inter_nproc, cp_npgrp, &
                                                ierr, iprc
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: cp_nolist, cp_nplist

!     ==--------------------------------------------------------------==
!
!     check if there are enough processes for the number of groups
!

    cp_npgrp = 0
    IF(parai%cp_me.EQ.0.AND.parai%cp_nogrp.LE.0) CALL stopgm(proceduren, &
         &     'number of cp_groups is less than zero.',&
         & __LINE__,__FILE__)
    IF(parai%cp_nogrp.GT.parai%cp_nproc) THEN
       IF(parai%cp_me.EQ.0) WRITE(6,'(1x,a)') 'warning: the number of ' &
            &        //'cp_groups is greater than the available number ' &
            &        //'of processes.'
       IF(parai%cp_me.EQ.0) WRITE(6,'(1x,a)') 'warning: the number of '&
            &        //'cp_groups will be reset!'
       CALL m_flush(6)
       parai%cp_nogrp = 0
    ELSEIF(MOD(parai%cp_nproc,parai%cp_nogrp).NE.0) THEN
       IF(parai%cp_me.EQ.0) WRITE(6,'(1x,a)') 'warning: the number of '&
            &        //'cp_groups is not a divisor of the available '&
            &        //'number of processes.'
       IF(parai%cp_me.EQ.0) WRITE(6,'(1x,a)') 'warning: the number of '&
            &        //'cp_groups will be reset!'
       CALL m_flush(6)
       parai%cp_nogrp = 0
    ELSE
       cp_npgrp = parai%cp_nproc / parai%cp_nogrp
    ENDIF
    !
    !     reset the cp group communicator
    !
    ALLOCATE(cp_nolist(parai%cp_nproc),cp_nplist(parai%cp_nproc),stat=ierr)
    IF(ierr.NE.0) CALL stopgm(proceduren,'allocation problem',&
         & __LINE__,__FILE__)
    CALL mp_cart(parai%cp_grp,parai%cp_nogrp,cp_npgrp,cp_nolist,cp_nplist,&
         &     parai%cp_inter_grp,parai%allgrp)
    DEALLOCATE(cp_nolist,cp_nplist,stat=ierr)
    IF(ierr.NE.0) CALL stopgm(proceduren,'deallocation problem',&
         & __LINE__,__FILE__)
    CALL mp_environ(parai%cp_inter_grp,cp_inter_nproc,parai%cp_inter_me)
    CALL mp_environ(parai%allgrp,parai%nproc,parai%me)
    parai%mepos = parai%me !set this useful variable as well
    !
    !     set up the source process within the group (allgrp communicator) 
    !
    parai%source = 0                ! range from 0 to nproc - 1
    paral%parent = parai%me.EQ.parai%source
    !
    !     set up the io process (cp_grp communicator)
    !
    !      io_source = 0            ! should range from 0 to cp_nproc - 1
    !                               ! and independant of source...
    !                               ! but additional work is likely be needed !
    !     >>> so we hack
    parai%io_source = 0
    IF(parai%me.EQ.parai%source.AND.parai%cp_inter_me.EQ.0) parai%io_source = parai%cp_me
    CALL mp_max(parai%io_source,parai%cp_grp)
    !     <<<
    paral%io_parent = parai%cp_me.EQ.parai%io_source
    !
    !     set the source process for the intergroup communicator
    !
    parai%cp_inter_io_source = 0
    IF(parai%cp_me.EQ.parai%io_source) parai%cp_inter_io_source = parai%cp_inter_me
    CALL mp_max(parai%cp_inter_io_source,parai%cp_grp)
    paral%cp_inter_io_parent = parai%cp_inter_me.EQ.parai%cp_inter_io_source
    !
    !     set mapping (me,cp_inter_me) to cp_me
    !
    CALL finalize_cp_grp()
    ALLOCATE(cp_grp_get_cp_rank(0:parai%nproc-1,0:cp_inter_nproc),&
         &     stat=ierr)
    IF(ierr.NE.0) CALL stopgm(proceduren,'allocation problem',&
         & __LINE__,__FILE__)
    cp_grp_get_cp_rank = 0
    cp_grp_get_cp_rank(parai%me,parai%cp_inter_me) = parai%cp_me
    CALL mp_sum(cp_grp_get_cp_rank,parai%cp_nproc,parai%cp_grp)
    !
    !     some tests on number of groups/processors deemed essential 
    !
    IF(cp_inter_nproc.NE.parai%cp_nogrp) CALL stopgm(proceduren,&
         &     'inconsistent cp_inter_nproc/cp_nogrp',&
         & __LINE__,__FILE__)
    IF(cp_npgrp.NE.parai%nproc) CALL stopgm(proceduren,&
         &     'inconsistent cp_npgrp/nproc',&
         & __LINE__,__FILE__)

    IF(grandparent) THEN
       WRITE(6,'(1x,2(a,i0),a)') 'CP_GROUPS: USING A ',&
            &   parai%cp_nogrp,' x ',cp_npgrp,' GRID (GROUPS X NPROCS) PER PC_GROUP.'
       CALL m_flush(6)
    ENDIF

    IF(PRINT_GROUP_INFOS) THEN
       DO iprc=0,parai%cp_nproc-1
          IF(iprc.EQ.parai%cp_me) THEN
             WRITE(6,'(3(a,2(1x,i0)),3(a,l2))')&
                  &        ' cp_grp',parai%cp_nproc,parai%cp_me,&
                  &        ' allgrp',parai%nproc,parai%me,&
                  &        ' cp_inter_grp',cp_inter_nproc,parai%cp_inter_me,&
                  &        ' parent',paral%parent,' cp_inter_io_parent',paral%cp_inter_io_parent,&
                  &        ' io_parent',paral%io_parent
          ENDIF
          CALL m_flush(6)
          CALL mp_sync(parai%cp_grp)
       ENDDO
    ENDIF
    !     ==--------------------------------------------------------------==
  END SUBROUTINE reset_cp_grp

  ! ==================================================================
  SUBROUTINE finalize_cp_grp()
    ! ==--------------------------------------------------------------==
    CHARACTER(*), PARAMETER                  :: procedureN = 'finalize_cp_grp'

    INTEGER                                  :: ierr

! ==--------------------------------------------------------------==

    IF (ALLOCATED(cp_grp_get_cp_rank)) THEN
       DEALLOCATE(cp_grp_get_cp_rank,stat=ierr)
       IF (ierr.NE.0) CALL stopgm(proceduren,'deallocation problem',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
  END SUBROUTINE finalize_cp_grp
  ! ==================================================================

END MODULE set_cp_grp_utils
