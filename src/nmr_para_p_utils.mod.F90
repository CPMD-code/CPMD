MODULE nmr_para_p_utils
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_group,&
                                             mp_sync
  USE parac,                           ONLY: parai,&
                                             paral
  USE response_pmod,                   ONLY: nmr_para
  USE system,                          ONLY: parap

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: nmr_para_p
  PUBLIC :: rewrite_output

CONTAINS

  ! ==================================================================
  SUBROUTINE nmr_para_p
    ! ==--------------------------------------------------------------==

    CHARACTER(len=1), DIMENSION(0:5), PARAMETER :: &
      zahlenchars = (/'0','1','2','3','4','5'/)

    INTEGER                                  :: i, myfirst, nproc_L, nproc_p

! SAVE OLD VALUES:

    nmr_para%nmr_superparent   = paral%parent
    nmr_para%nmr_supergroup    = parai%allgrp
    nmr_para%nmr_total_nproc   = parai%nproc
    nmr_para%nmr_supersource   = parai%source


    IF (parai%nproc/REAL(nmr_para%nmr_threads,kind=real_8) .LT. 2.99_real_8)&
         nmr_para%nmr_threads = 1


    IF (paral%io_parent)&
         WRITE (6,&
         '("RESPONSE SUPERPARALLEL: Initializing ",I1,'//&
         '" parallel threads.")') nmr_para%nmr_threads

    ! It is problematic if a groundstate-wf-optimization AND a NMR 
    ! perturbation run is desired. Then only one thread is possible.
    ! The thread business only works well if nothing (except a response) is
    ! calculated.
    ! ==--------------------------------------------------------------==
    IF (nmr_para%nmr_threads.EQ.6) THEN
       ! initialization
       DO i=0,nmr_para%nmr_total_nproc
          parap%pgroup(i) = -1
          parap%nlink(i)  = -1
       ENDDO
       DO i=0,6
          nmr_para%parents(i)=-1
       ENDDO

       ! repartitioning of processors. The p-tasks get 20%-30% more
       ! cpus than the L-tasks.
       nproc_L = INT( nmr_para%nmr_total_nproc / 6.6_real_8 )
       nproc_p = INT((nmr_para%nmr_total_nproc - 3 * nproc_L)/3.0_real_8)

       nmr_para%parents(1) = 0
       nmr_para%parents(2) =   nproc_p
       nmr_para%parents(3) = 2*nproc_p
       nmr_para%parents(4) = 3*nproc_p
       nmr_para%parents(5) = 3*nproc_p +   nproc_L
       nmr_para%parents(6) = 3*nproc_p + 2*nproc_L
       nmr_para%parents(7) = nmr_para%nmr_total_nproc
       parai%nproc = 0
       DO i=1,nmr_para%nmr_threads
          IF (nmr_para%parents(i).EQ.parai%me) paral%parent = .TRUE.
          IF (nmr_para%parents(i).LE.parai%me .AND. nmr_para%parents(i+1).GT.parai%me) THEN
             nmr_para%nmr_mygroup = i-1! groups go from 0 to 5
             myfirst     = nmr_para%parents(i)
             IF (parai%nproc .NE. 0) CALL stopgm('NMR_PARA',&
                  'FATAL: Belonging to two threads.',& 
                  __LINE__,__FILE__)
             parai%nproc       = nmr_para%parents(i+1)-nmr_para%parents(i)
          ENDIF
       ENDDO
       IF (parai%nproc .EQ. 0) CALL stopgm('NMR_PARA',&
            'FATAL: Belonging to NO thread at all.',& 
            __LINE__,__FILE__)
       DO i=1,parai%nproc
          parap%pgroup(i) = myfirst + i-1
          parap%nlink(myfirst + i-1) = i-1
       ENDDO
       parai%mepos = parai%me - myfirst
       parai%source = 0
       CALL mp_group(parai%nproc,parap%pgroup,parai%allgrp,nmr_para%nmr_supergroup)
       ! ==--------------------------------------------------------------==
    ELSEIF (nmr_para%nmr_threads.EQ.3) THEN
       ! initialization
       DO i=0,nmr_para%nmr_total_nproc
          parap%pgroup(i) = -1
          parap%nlink(i)  = -1
       ENDDO
       DO i=0,6
          nmr_para%parents(i)=-1
       ENDDO

       nproc_L = INT( nmr_para%nmr_total_nproc / 3.0_real_8 )

       nmr_para%parents(1) = 0
       nmr_para%parents(2) =   nproc_L
       nmr_para%parents(3) = 2*nproc_L
       nmr_para%parents(4) = nmr_para%nmr_total_nproc

       parai%nproc = 0
       DO i=1,nmr_para%nmr_threads
          IF (nmr_para%parents(i).EQ.parai%me) paral%parent = .TRUE.
          IF (nmr_para%parents(i).LE.parai%me .AND. nmr_para%parents(i+1).GT.parai%me) THEN
             nmr_para%nmr_mygroup = i-1! groups go from 0 to 2
             myfirst     = nmr_para%parents(i)
             IF (parai%nproc .NE. 0) CALL stopgm('NMR_PARA',&
                  'FATAL: Belonging to two threads.',& 
                  __LINE__,__FILE__)
             parai%nproc       = nmr_para%parents(i+1)-nmr_para%parents(i)
          ENDIF
       ENDDO
       IF (parai%nproc .EQ. 0) CALL stopgm('NMR_PARA',&
            'FATAL: Belonging to NO thread at all.',& 
            __LINE__,__FILE__)
       DO i=1,parai%nproc
          parap%pgroup(i) = myfirst + i-1
          parap%nlink(myfirst + i-1) = i-1
       ENDDO
       parai%mepos = parai%me - myfirst
       parai%source = 0
       CALL mp_group(parai%nproc,parap%pgroup,parai%allgrp,nmr_para%nmr_supergroup)
       ! ==--------------------------------------------------------------==
    ELSEIF (nmr_para%nmr_threads.EQ.1) THEN
       ! nothing to do. Except:
       nmr_para%nmr_mygroup = 0
       nmr_para%parents(1)  = 0
    ELSE
       CALL stopgm('NMR_PARA','INVALID NUMBER OF THREADS.',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==

    ! Informative output:
    DO i=1,nmr_para%nmr_threads
       IF (parai%me.EQ.nmr_para%parents(i)) THEN
          IF (nmr_para%nmr_superparent) THEN
             IF (paral%io_parent)&
                  WRITE (6,'(A,I3,A,I1,A,I3,A)')&
                  "PROCESSOR ",parai%me," (GROUP ",nmr_para%nmr_mygroup,&
                  " with ",parai%nproc," procs): I am the superparent. "
          ELSE
             IF (paral%io_parent)&
                  WRITE (6,'(A,I3,A,I1,A,I3,A,A5,A)')&
                  "PROCESSOR ",parai%me," (GROUP ",nmr_para%nmr_mygroup,&
                  " with ",parai%nproc," procs): Output --> out_",&
                  zahlenchars(nmr_para%nmr_mygroup),"."

             IF (paral%io_parent)&
                  CLOSE(6)
             IF (paral%io_parent)&
                  OPEN(unit=6,file='out_'//zahlenchars(nmr_para%nmr_mygroup))
          ENDIF
       ENDIF
       CALL mp_sync(nmr_para%nmr_supergroup)
    ENDDO


    RETURN
  END SUBROUTINE nmr_para_p
  ! ==================================================================
  SUBROUTINE rewrite_output
    ! ==--------------------------------------------------------------==
    CHARACTER(len=1), DIMENSION(0:5), PARAMETER :: &
      zahlenchars = (/'0','1','2','3','4','5'/)

    CHARACTER(len=128)                       :: transfer
    INTEGER                                  :: i

97  FORMAT (10("="),10(" "),a15,i1,10(" "),20("="))
98  FORMAT (66("="))

    IF (paral%io_parent .AND. .NOT. nmr_para%nmr_superparent)&
         CLOSE(6)
    CALL mp_sync(nmr_para%nmr_supergroup)

    IF (nmr_para%nmr_superparent) THEN
       DO i=1,nmr_para%nmr_threads-1
          IF (paral%io_parent)&
               WRITE (6,98)
          IF (paral%io_parent)&
               WRITE (6,97) 'OUTPUT GROUP ',i
          IF (paral%io_parent)&
               WRITE (6,98)
          IF (paral%io_parent)&
               OPEN(unit=9,file='out_'//zahlenchars(i))

20        READ (unit=9,fmt='(A128)',END=29) transfer
          IF (paral%io_parent)&
               WRITE (6,'(A)') TRANSFER(1:66)
          GOTO 20


29        IF (paral%io_parent)  CLOSE(unit=9)

          ! Delete files:
          IF (paral%io_parent)&
               OPEN(unit=9,file='out_'//zahlenchars(i))
          IF (paral%io_parent)&
               CLOSE(unit=9,status='DELETE')

       ENDDO
    ENDIF

    RETURN
  END SUBROUTINE rewrite_output
  ! ==================================================================

END MODULE nmr_para_p_utils
