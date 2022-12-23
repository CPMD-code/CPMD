SUBROUTINE mp_bcast_byte(DATA,n,root,comm)
  ! ==--------------------------------------------------------------==
  ! == Wrapper to mpi_bcast                                         ==
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE mp_interface, ONLY: mp_mpi_error_assert
#ifdef __PARALLEL
  USE mpi
#endif
  IMPLICIT NONE
  ! Arguments
  INTEGER :: DATA(*)
  INTEGER :: n,root,comm
  ! Variables
  INTEGER :: ierr
  CHARACTER(*),PARAMETER :: procedureN='mp_bcast_byte'
#ifdef __PARALLEL
  ! ==--------------------------------------------------------------==
  CALL mpi_bcast(DATA,n,mpi_byte,root,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE mp_bcast_byte


SUBROUTINE my_stopall(code)
  ! ==--------------------------------------------------------------==
  ! == STOP ALL PROCESSORS                                          ==
  ! ==--------------------------------------------------------------==
  USE mp_interface, ONLY: mp_comm_world
#ifdef __PARALLEL
  USE mpi
#endif
  IMPLICIT NONE
  ! Arguments
  INTEGER :: code
#ifdef __PARALLEL
  ! Variables
  INTEGER :: ierr
  LOGICAL :: is_init
  ! ==--------------------------------------------------------------==
  CALL mpi_initialized(is_init,ierr)
  if( is_init ) then
     CALL mpi_abort(mp_comm_world,code,ierr)
  else
     stop 1
  endif
#endif
  STOP 1
  ! ==--------------------------------------------------------------==
END SUBROUTINE my_stopall

! ==================================================================
SUBROUTINE my_concat(outmsg,inmsg,blklen,gid)
  ! ==--------------------------------------------------------------==
  ! == Concat all OUTMSG from GID group_ processor into INMSG        ==
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE mp_interface, ONLY: mp_mpi_error_assert
  USE utils, ONLY : icopy
#ifdef __PARALLEL
  USE mpi
#endif
  IMPLICIT NONE
  ! Arguments
  INTEGER :: outmsg(*),inmsg(*),blklen,gid
  CHARACTER(*),PARAMETER::procedureN='my_concat'
#ifdef __PARALLEL
  ! Variables
  INTEGER :: ierr
  ! ==--------------------------------------------------------------==
  CALL mpi_allgather(outmsg,blklen,mpi_byte,inmsg,blklen,&
       mpi_byte,gid,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
  ! Variables
  INTEGER :: len
  ! ==--------------------------------------------------------------==
  len = blklen/(8/2)
  CALL icopy(len,outmsg,1,inmsg,1)
#endif
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE my_concat
! ==================================================================
SUBROUTINE my_concatv(outmsg,inmsg,blklen,recvcnt,recvdispl,gid)
  ! ==--------------------------------------------------------------==
  ! == Concat all OUTMSG of different lengths from GID group_        ==
  ! == processor into INMSG. NOTE DOUBLE PRECISION ARRAYS USED      ==
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE mp_interface, ONLY: mp_mpi_error_assert
#ifdef __PARALLEL
  USE mpi
#endif
  IMPLICIT NONE
  ! Arguments
  REAL(real_8) :: outmsg(*),inmsg(*)
  INTEGER :: blklen,gid,recvcnt(*),recvdispl(*)
  CHARACTER(*),PARAMETER::procedureN='my_concatv'
#ifdef __PARALLEL
  ! Variables
  INTEGER :: ierr
  ! ==--------------------------------------------------------------==
  CALL mpi_allgatherv(outmsg,blklen,mpi_double_precision,inmsg,&
       recvcnt,recvdispl,mpi_double_precision&
       ,gid,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
  ! ==--------------------------------------------------------------==
  CALL dcopy(blklen,outmsg,1,inmsg(1+recvdispl(1)),1)
#endif
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE my_concatv
! ==================================================================
! ==================================================================
SUBROUTINE my_source_concatv(outmsg,inmsg,blklen,recvcnt,&
     recvdispl,root,gid)
  ! ==--------------------------------------------------------------==
  ! == Concat all OUTMSG of different lengths from GID group_        ==
  ! == processor into INMSG. NOTE DOUBLE PRECISION ARRAYS USED      ==
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE mp_interface, ONLY: mp_mpi_error_assert
#ifdef __PARALLEL
  USE mpi
#endif
  IMPLICIT NONE
  ! Arguments
  REAL(real_8) :: outmsg(*),inmsg(*)
  INTEGER :: blklen,gid,recvcnt(*),recvdispl(*),root
  CHARACTER(*),PARAMETER::procedureN='my_source_concatv'
#ifdef __PARALLEL
  ! Variables
  INTEGER :: ierr
  INTEGER :: err_msg_len, ierr_err
  ! ==--------------------------------------------------------------==
  CALL mpi_gatherv(outmsg,blklen,mpi_double_precision,inmsg,&
       recvcnt,recvdispl,mpi_double_precision&
       ,root,gid,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
  ! ==--------------------------------------------------------------==
  CALL dcopy(blklen,outmsg,1,inmsg(1+recvdispl(1)),1)
#endif
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE my_source_concatv
! ==================================================================
SUBROUTINE my_trans(outmsg,inmsg,blklen,group_)
  ! ==--------------------------------------------------------------==
  ! == Matrix Transpose: group_ 1 : all processors                   ==
  ! ==                   group_ 2 : orbital split                    ==
  ! ==                   group_ 3 : plane wave split                 ==
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE machine, ONLY: m_walltime
  USE mp_interface, ONLY: mp_mpi_error_assert
  USE system , ONLY:group
  USE parac, ONLY : paral,parai
  USE pstat , ONLY:cmcal,cmlen,cmtim,ipar_aall
  USE utils, ONLY : icopy
#ifdef __PARALLEL
  USE mpi
#endif
  IMPLICIT NONE
  ! Arguments
  INTEGER :: outmsg(*),inmsg(*),blklen,gid,group_
  ! Variables
  INTEGER :: nn
  REAL(real_8) :: tim1,tim2
  CHARACTER(*),PARAMETER::procedureN='my_trans'
#ifdef __PARALLEL
  INTEGER :: ierr
#else
  INTEGER :: len
#endif
  ! ==--------------------------------------------------------------==
  cmcal(ipar_aall)=cmcal(ipar_aall)+1.0_real_8
  tim1=m_walltime()
  IF (group_.EQ.1) THEN
     gid=parai%allgrp
     nn=parai%nproc
  ELSEIF (group_.EQ.2) THEN
     gid=group%meogrp
     nn=group%nogrp
  ELSEIF (group_.EQ.3) THEN
     gid=group%mepgrp
     nn=group%npgrp
  ENDIF
#ifdef __PARALLEL
  CALL mpi_alltoall(outmsg,blklen,mpi_byte,inmsg,blklen,&
       mpi_byte,gid,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
  len = blklen/(8/2)
  CALL icopy(len,outmsg(1),1,inmsg(1),1)
#endif
  tim2=m_walltime()
  cmtim(ipar_aall)=cmtim(ipar_aall)+tim2-tim1
  cmlen(ipar_aall)=cmlen(ipar_aall)+blklen*nn
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE my_trans
! ==================================================================
SUBROUTINE my_allgather_i(DATA,n,comm)
  ! ==--------------------------------------------------------------==
  ! == Wrapper to mpi_allgather                                     ==
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE mp_interface, ONLY: mp_mpi_error_assert
#ifdef __PARALLEL
  USE mpi
#endif
  IMPLICIT NONE
  ! Arguments
  INTEGER :: DATA(*)
  INTEGER :: n,comm
  ! Variables
  CHARACTER(*),PARAMETER :: procedureN='MY_ALLGATHER_I'
#ifdef __PARALLEL
  INTEGER :: ierr
  ! ==--------------------------------------------------------------==
  CALL mpi_allgather(mpi_in_place,n,mpi_integer,DATA,n,&
       mpi_integer,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE my_allgather_i
! ==================================================================
SUBROUTINE my_shift(msend,mrecv,msglen,mep,ip,gid)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE mp_interface, ONLY: mp_mpi_error_assert
  USE parac, ONLY : paral,parai
  USE pstat , ONLY:cmcal,cmlen,cmtim,ipar_aall, ipar_send
#ifdef __PARALLEL
  USE mpi
#endif
  IMPLICIT NONE
  ! Arguments
  INTEGER :: msglen,mep,ip,gid
  REAL(real_8) :: msend(*),mrecv(*)
  CHARACTER(*),PARAMETER::procedureN='my_shift'
#ifdef __PARALLEL
  ! Variables
  INTEGER :: status(mpi_status_size),ipsend,iprecv,&
       whoami,howmany,&
       itype,irequest,ierr
  ! ==--------------------------------------------------------------==
  CALL mpi_comm_rank(gid,whoami,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  CALL mpi_comm_size(gid,howmany,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  ipsend=MOD(whoami+ip+howmany,howmany)
  iprecv=MOD(whoami-ip+howmany,howmany)

  cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  cmlen(ipar_send)=cmlen(ipar_send)+msglen
  itype=1
  irequest=1
  CALL mpi_isend(msend,msglen,mpi_byte,ipsend,itype,&
       gid,irequest,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  CALL mpi_recv(mrecv,msglen,mpi_byte,iprecv,itype,&
       gid,status,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  CALL mpi_wait(irequest,status,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
  CALL stopgm('MY_SHIFT','NOT INSTALLED',& 
       __LINE__,__FILE__)
#endif
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE my_shift
