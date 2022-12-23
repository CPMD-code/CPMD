MODULE startpa_utils
  USE error_handling,                  ONLY: stopgm
  USE fft_utils,                       ONLY: alloc_fft
  USE geq0mod,                         ONLY: geq0
  USE machine,                         ONLY: m_getarg,&
                                             m_iargc
  USE mp_interface,                    ONLY: mp_environ,&
                                             mp_start,&
                                             mp_sync,&
                                             mp_task_query
  USE parac,                           ONLY: parai
  USE pimd_utils,                      ONLY: alloc_pimd
  USE readsr_utils,                    ONLY: readsi
  USE set_cp_grp_utils,                ONLY: set_cp_grp
  USE system,                          ONLY: cnts,&
                                             parap
  USE system_utils,                    ONLY: alloc_system
  USE utils,                           ONLY: numcpus
  USE vdwcmod_utils,                   ONLY: alloc_vdwcmod
  USE zeroing_utils,                   ONLY: zeroing
#include "sizeof.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: startpa

CONTAINS

  ! ==================================================================
  SUBROUTINE startpa
    ! ==--------------------------------------------------------------==
    CHARACTER(*), PARAMETER                  :: procedureN = 'startpa'
    INTEGER, PARAMETER                       :: IUNIT = 5, LINELEN = 80, &
                                                PATHLEN = 255 

    CHARACTER(len=linelen)                   :: line
    INTEGER                                  :: ie, iostat, ip, ipp
    LOGICAL                                  :: erread

#ifdef __ES
    INTEGER :: icarg
    INTEGER, PARAMETER :: linelen=130 
    CHARACTER (len=linelen) :: filename
    INTEGER :: lef,iend,i
    CHARACTER(len=130) :: exit_fpath ! PATH TO EXIT FILE ON ES
    COMMON/fpath/ exit_fpath
#endif
    ! ==--------------------------------------------------------------==
    ! PGROUP(I) : Pointer position --> PE number
    ! NLINK(I)  : Pointer PE number --> position - 1
    CALL mp_start
    CALL mp_task_query(parai%cp_grp)
    CALL mp_environ(parai%cp_grp,parai%cp_nproc,parai%cp_me)
    ! 
    parai%allgrp=-123456789;parai%nproc=-1234567;parai%me=-123456;parai%mepos=-12345;

    ! ==--------------------------------------------------------------==
    ! 
    ! get the number of cp groups from input
    ! 
    parai%cp_nogrp = 1
    IF (parai%cp_me.EQ.0) THEN
       OPEN(unit=iunit,file=cnts%inputfile,status='OLD',iostat=iostat)
       IF (iostat.NE.0)CALL stopgm(procedureN,&
            'Cannot open the input file',& 
            __LINE__,__FILE__)
       outer: DO
          READ(iunit,iostat=iostat,fmt='(A80)') line
          IF (iostat.NE.0) EXIT outer
          IF (INDEX(line,'&CPMD').NE.0) THEN
             inner: DO
                READ(iunit,iostat=iostat,fmt='(A80)') linE
                IF (iostat.NE.0) EXIT outer
                IF (INDEX(line,'&END').NE.0) EXIT outer
                IF (INDEX(line,'CP_GROUPS').NE.0) THEN
                   READ(iunit,iostat=iostat,fmt='(A80)') linE
                   IF (iostat.NE.0)CALL stopgm(procedureN,&
                        'error while reading a line',& 
                        __LINE__,__FILE__)
                   CALL readsi(line,1,ie,parai%cp_nogrp,erread)
                   IF (erread)CALL stopgm(procedureN,&
                        'error while reading an integer',& 
                        __LINE__,__FILE__)
                ENDIF
             ENDDO inner
          ENDIF
       ENDDO outer
       CLOSE(unit=iunit,iostat=iostat)
       IF (iostat.NE.0)CALL stopgm(procedureN,&
            'Cannot close the input file',& 
            __LINE__,__FILE__)
    ENDIF
    CALL mp_bcast_byte(parai%cp_nogrp,size_in_bytes_of(parai%cp_nogrp),0,parai%cp_grp)
    ! ==--------------------------------------------------------------==
    ! 
    ! set the cp groups, parent, io_parent, source and io_source
    ! 
    CALL set_cp_grp()

    ! ==--------------------------------------------------------------==
    ! alloc some global arrays in system (alloc with max number of procs ie parai%cp_nproc)
    CALL alloc_system( parai%cp_nproc )
    CALL alloc_pimd( parai%cp_nproc )
    CALL alloc_fft( parai%cp_nproc )
    CALL alloc_vdwcmod( parai%cp_nproc )

    ! ==--------------------------------------------------------------==
    ! 
    ! some setup
    ! 
    geq0=.TRUE.
    CALL zeroing(parap%pgroup)
    DO ip=1,parai%cp_nproc
       parap%pgroup(ip) = ip-1
    ENDDO
    CALL mp_sync(parai%cp_grp)
    CALL zeroing(parap%nlink)
    DO ip=1,parai%cp_nproc
       ipp=parap%pgroup(ip)
       parap%nlink(ipp)=ip-1
    ENDDO


    CALL numcpus(parai%ncpus)
    ! 
#ifdef __ES
    icarg=m_iargc()
    IF (icarg.GT.2) THEN
       CALL m_getarg(3,filename)
    ELSE
       filename='./cpmd.log'
    ENDIF
    lef=LEN(TRIM(filename))
    iend=lef
    DO i=lef,1,-1
       IF (filename(i:i).EQ.'/') THEN
          iend=i
          EXIT
       ENDIF
    ENDDO
    IF (iend.EQ.lef) THEN
       CALL stopgm("STARTPA","EXIT_FPATH UNDEFINED",& 
            __LINE__,__FILE__)
    ENDIF
    exit_fpath=filename(1:iend)
#endif
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE startpa
  ! ==================================================================

END MODULE startpa_utils
