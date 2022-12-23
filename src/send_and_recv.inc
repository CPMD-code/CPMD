SUBROUTINE mp_send_char_r0(DATA, dest,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments 
  INTEGER, INTENT(in) ::  dest,tag,comm
  CHARACTER(*), INTENT(in)  :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  CHARACTER(*),PARAMETER::procedureN='mp_isend_char_r0'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+len(data)*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_send(DATA,LEN(DATA),MPI_CHARACTER,dest,tag,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_send_char_r0

SUBROUTINE mp_recv_char_r0(DATA, src,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments
  INTEGER, INTENT(in) ::  src,tag,comm
  CHARACTER(*), INTENT(out)  :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_recv_char_r0'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+len(data)*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_recv(DATA,LEN(DATA),MPI_CHARACTER,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_recv_char_r0
SUBROUTINE mp_send_char_r1(DATA,n, dest,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments 
  INTEGER, INTENT(in) :: n, dest,tag,comm
  CHARACTER(*), INTENT(in) , DIMENSION(*) :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_isend_char_r1'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+len(data)*n*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_send(DATA,LEN(DATA)*n,MPI_CHARACTER,dest,tag,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_send_char_r1

SUBROUTINE mp_recv_char_r1(DATA,n, src,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments
  INTEGER, INTENT(in) :: n, src,tag,comm
  CHARACTER(*), INTENT(out) , DIMENSION(*) :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_recv_char_r1'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+len(data)*n*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_recv(DATA,LEN(DATA)*n,MPI_CHARACTER,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_recv_char_r1
SUBROUTINE mp_send_int4_r0(DATA, dest,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments 
  INTEGER, INTENT(in) ::  dest,tag,comm
  INTEGER(int_4), INTENT(in)  :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_isend_int4_r0'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+1*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_send(DATA,1,MPI_INTEGER,dest,tag,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_send_int4_r0

SUBROUTINE mp_recv_int4_r0(DATA, src,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments
  INTEGER, INTENT(in) ::  src,tag,comm
  INTEGER(int_4), INTENT(out)  :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_recv_int4_r0'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+1*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_recv(DATA,1,MPI_INTEGER,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_recv_int4_r0
SUBROUTINE mp_send_int4_r1(DATA,n, dest,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments 
  INTEGER, INTENT(in) :: n, dest,tag,comm
  INTEGER(int_4), INTENT(in) , DIMENSION(*) :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_isend_int4_r1'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+n*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_send(DATA,n,MPI_INTEGER,dest,tag,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_send_int4_r1

SUBROUTINE mp_recv_int4_r1(DATA,n, src,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments
  INTEGER, INTENT(in) :: n, src,tag,comm
  INTEGER(int_4), INTENT(out) , DIMENSION(*) :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_recv_int4_r1'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+n*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_recv(DATA,n,MPI_INTEGER,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_recv_int4_r1
SUBROUTINE mp_send_int4_r2(DATA,n, dest,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments 
  INTEGER, INTENT(in) :: n, dest,tag,comm
  INTEGER(int_4), INTENT(in) , DIMENSION(1,*) :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_isend_int4_r2'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+n*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_send(DATA,n,MPI_INTEGER,dest,tag,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_send_int4_r2

SUBROUTINE mp_recv_int4_r2(DATA,n, src,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments
  INTEGER, INTENT(in) :: n, src,tag,comm
  INTEGER(int_4), INTENT(out) , DIMENSION(1,*) :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_recv_int4_r2'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+n*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_recv(DATA,n,MPI_INTEGER,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_recv_int4_r2
SUBROUTINE mp_send_int8_r0(DATA, dest,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments 
  INTEGER, INTENT(in) ::  dest,tag,comm
  INTEGER(int_8), INTENT(in)  :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_isend_int8_r0'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+1*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_send(DATA,1,MPI_INTEGER8,dest,tag,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_send_int8_r0

SUBROUTINE mp_recv_int8_r0(DATA, src,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments
  INTEGER, INTENT(in) ::  src,tag,comm
  INTEGER(int_8), INTENT(out)  :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_recv_int8_r0'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+1*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_recv(DATA,1,MPI_INTEGER8,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_recv_int8_r0
SUBROUTINE mp_send_int8_r1(DATA,n, dest,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments 
  INTEGER, INTENT(in) :: n, dest,tag,comm
  INTEGER(int_8), INTENT(in) , DIMENSION(*) :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_isend_int8_r1'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+n*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_send(DATA,n,MPI_INTEGER8,dest,tag,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_send_int8_r1

SUBROUTINE mp_recv_int8_r1(DATA,n, src,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments
  INTEGER, INTENT(in) :: n, src,tag,comm
  INTEGER(int_8), INTENT(out) , DIMENSION(*) :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_recv_int8_r1'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+n*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_recv(DATA,n,MPI_INTEGER8,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_recv_int8_r1
SUBROUTINE mp_send_real8_r0(DATA, dest,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments 
  INTEGER, INTENT(in) ::  dest,tag,comm
  REAL(real_8), INTENT(in)  :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_isend_real8_r0'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+1*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_send(DATA,1,MPI_DOUBLE_PRECISION,dest,tag,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_send_real8_r0

SUBROUTINE mp_recv_real8_r0(DATA, src,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments
  INTEGER, INTENT(in) ::  src,tag,comm
  REAL(real_8), INTENT(out)  :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_recv_real8_r0'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+1*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_recv(DATA,1,MPI_DOUBLE_PRECISION,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_recv_real8_r0
SUBROUTINE mp_send_real8_r1(DATA,n, dest,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments 
  INTEGER, INTENT(in) :: n, dest,tag,comm
  REAL(real_8), INTENT(in) , DIMENSION(*) :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_isend_real8_r1'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+n*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_send(DATA,n,MPI_DOUBLE_PRECISION,dest,tag,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_send_real8_r1

SUBROUTINE mp_recv_real8_r1(DATA,n, src,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments
  INTEGER, INTENT(in) :: n, src,tag,comm
  REAL(real_8), INTENT(out) , DIMENSION(*) :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_recv_real8_r1'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+n*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_recv(DATA,n,MPI_DOUBLE_PRECISION,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_recv_real8_r1
SUBROUTINE mp_send_real8_r2(DATA,n, dest,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments 
  INTEGER, INTENT(in) :: n, dest,tag,comm
  REAL(real_8), INTENT(in) , DIMENSION(1,*) :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_isend_real8_r2'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+n*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_send(DATA,n,MPI_DOUBLE_PRECISION,dest,tag,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_send_real8_r2

SUBROUTINE mp_recv_real8_r2(DATA,n, src,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments
  INTEGER, INTENT(in) :: n, src,tag,comm
  REAL(real_8), INTENT(out) , DIMENSION(1,*) :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_recv_real8_r2'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+n*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_recv(DATA,n,MPI_DOUBLE_PRECISION,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_recv_real8_r2
SUBROUTINE mp_send_real8_r3(DATA,n, dest,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments 
  INTEGER, INTENT(in) :: n, dest,tag,comm
  REAL(real_8), INTENT(in) , DIMENSION(1,1,*) :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_isend_real8_r3'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+n*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_send(DATA,n,MPI_DOUBLE_PRECISION,dest,tag,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_send_real8_r3

SUBROUTINE mp_recv_real8_r3(DATA,n, src,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments
  INTEGER, INTENT(in) :: n, src,tag,comm
  REAL(real_8), INTENT(out) , DIMENSION(1,1,*) :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_recv_real8_r3'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+n*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_recv(DATA,n,MPI_DOUBLE_PRECISION,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_recv_real8_r3
SUBROUTINE mp_send_real8_r4(DATA,n, dest,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments 
  INTEGER, INTENT(in) :: n, dest,tag,comm
  REAL(real_8), INTENT(in) , DIMENSION(1,1,1,*) :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_isend_real8_r4'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+n*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_send(DATA,n,MPI_DOUBLE_PRECISION,dest,tag,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_send_real8_r4

SUBROUTINE mp_recv_real8_r4(DATA,n, src,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments
  INTEGER, INTENT(in) :: n, src,tag,comm
  REAL(real_8), INTENT(out) , DIMENSION(1,1,1,*) :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_recv_real8_r4'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+n*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_recv(DATA,n,MPI_DOUBLE_PRECISION,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_recv_real8_r4
SUBROUTINE mp_send_real8_r5(DATA,n, dest,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments 
  INTEGER, INTENT(in) :: n, dest,tag,comm
  REAL(real_8), INTENT(in) , DIMENSION(1,1,1,1,*) :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_isend_real8_r5'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+n*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_send(DATA,n,MPI_DOUBLE_PRECISION,dest,tag,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_send_real8_r5

SUBROUTINE mp_recv_real8_r5(DATA,n, src,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments
  INTEGER, INTENT(in) :: n, src,tag,comm
  REAL(real_8), INTENT(out) , DIMENSION(1,1,1,1,*) :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_recv_real8_r5'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+n*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_recv(DATA,n,MPI_DOUBLE_PRECISION,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_recv_real8_r5
SUBROUTINE mp_send_complex4_r0(DATA, dest,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments 
  INTEGER, INTENT(in) ::  dest,tag,comm
  COMPLEX(real_4), INTENT(in)  :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_isend_complex4_r0'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+1*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_send(DATA,1,MPI_COMPLEX,dest,tag,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_send_complex4_r0

SUBROUTINE mp_recv_complex4_r0(DATA, src,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments
  INTEGER, INTENT(in) ::  src,tag,comm
  COMPLEX(real_4), INTENT(out)  :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_recv_complex4_r0'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+1*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_recv(DATA,1,MPI_COMPLEX,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_recv_complex4_r0
SUBROUTINE mp_send_complex4_r1(DATA,n, dest,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments 
  INTEGER, INTENT(in) :: n, dest,tag,comm
  COMPLEX(real_4), INTENT(in) , DIMENSION(*) :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_isend_complex4_r1'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+n*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_send(DATA,n,MPI_COMPLEX,dest,tag,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_send_complex4_r1

SUBROUTINE mp_recv_complex4_r1(DATA,n, src,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments
  INTEGER, INTENT(in) :: n, src,tag,comm
  COMPLEX(real_4), INTENT(out) , DIMENSION(*) :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_recv_complex4_r1'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+n*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_recv(DATA,n,MPI_COMPLEX,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_recv_complex4_r1
SUBROUTINE mp_send_complex4_r2(DATA,n, dest,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments 
  INTEGER, INTENT(in) :: n, dest,tag,comm
  COMPLEX(real_4), INTENT(in) , DIMENSION(1,*) :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_isend_complex4_r2'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+n*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_send(DATA,n,MPI_COMPLEX,dest,tag,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_send_complex4_r2

SUBROUTINE mp_recv_complex4_r2(DATA,n, src,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments
  INTEGER, INTENT(in) :: n, src,tag,comm
  COMPLEX(real_4), INTENT(out) , DIMENSION(1,*) :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_recv_complex4_r2'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+n*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_recv(DATA,n,MPI_COMPLEX,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_recv_complex4_r2
SUBROUTINE mp_send_complex4_r3(DATA,n, dest,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments 
  INTEGER, INTENT(in) :: n, dest,tag,comm
  COMPLEX(real_4), INTENT(in) , DIMENSION(1,1,*) :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_isend_complex4_r3'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+n*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_send(DATA,n,MPI_COMPLEX,dest,tag,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_send_complex4_r3

SUBROUTINE mp_recv_complex4_r3(DATA,n, src,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments
  INTEGER, INTENT(in) :: n, src,tag,comm
  COMPLEX(real_4), INTENT(out) , DIMENSION(1,1,*) :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_recv_complex4_r3'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+n*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_recv(DATA,n,MPI_COMPLEX,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_recv_complex4_r3
SUBROUTINE mp_send_complex4_r4(DATA,n, dest,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments 
  INTEGER, INTENT(in) :: n, dest,tag,comm
  COMPLEX(real_4), INTENT(in) , DIMENSION(1,1,1,*) :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_isend_complex4_r4'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+n*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_send(DATA,n,MPI_COMPLEX,dest,tag,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_send_complex4_r4

SUBROUTINE mp_recv_complex4_r4(DATA,n, src,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments
  INTEGER, INTENT(in) :: n, src,tag,comm
  COMPLEX(real_4), INTENT(out) , DIMENSION(1,1,1,*) :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_recv_complex4_r4'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+n*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_recv(DATA,n,MPI_COMPLEX,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_recv_complex4_r4
SUBROUTINE mp_send_complex8_r0(DATA, dest,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments 
  INTEGER, INTENT(in) ::  dest,tag,comm
  COMPLEX(real_8), INTENT(in)  :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_isend_complex8_r0'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+1*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_send(DATA,1,MPI_DOUBLE_COMPLEX,dest,tag,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_send_complex8_r0

SUBROUTINE mp_recv_complex8_r0(DATA, src,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments
  INTEGER, INTENT(in) ::  src,tag,comm
  COMPLEX(real_8), INTENT(out)  :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_recv_complex8_r0'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+1*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_recv(DATA,1,MPI_DOUBLE_COMPLEX,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_recv_complex8_r0
SUBROUTINE mp_send_complex8_r1(DATA,n, dest,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments 
  INTEGER, INTENT(in) :: n, dest,tag,comm
  COMPLEX(real_8), INTENT(in) , DIMENSION(*) :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_isend_complex8_r1'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+n*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_send(DATA,n,MPI_DOUBLE_COMPLEX,dest,tag,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_send_complex8_r1

SUBROUTINE mp_recv_complex8_r1(DATA,n, src,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments
  INTEGER, INTENT(in) :: n, src,tag,comm
  COMPLEX(real_8), INTENT(out) , DIMENSION(*) :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_recv_complex8_r1'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+n*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_recv(DATA,n,MPI_DOUBLE_COMPLEX,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_recv_complex8_r1
SUBROUTINE mp_send_complex8_r2(DATA,n, dest,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments 
  INTEGER, INTENT(in) :: n, dest,tag,comm
  COMPLEX(real_8), INTENT(in) , DIMENSION(1,*) :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_isend_complex8_r2'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+n*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_send(DATA,n,MPI_DOUBLE_COMPLEX,dest,tag,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_send_complex8_r2

SUBROUTINE mp_recv_complex8_r2(DATA,n, src,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments
  INTEGER, INTENT(in) :: n, src,tag,comm
  COMPLEX(real_8), INTENT(out) , DIMENSION(1,*) :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_recv_complex8_r2'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+n*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_recv(DATA,n,MPI_DOUBLE_COMPLEX,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_recv_complex8_r2
SUBROUTINE mp_send_complex8_r3(DATA,n, dest,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments 
  INTEGER, INTENT(in) :: n, dest,tag,comm
  COMPLEX(real_8), INTENT(in) , DIMENSION(1,1,*) :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_isend_complex8_r3'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+n*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_send(DATA,n,MPI_DOUBLE_COMPLEX,dest,tag,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_send_complex8_r3

SUBROUTINE mp_recv_complex8_r3(DATA,n, src,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments
  INTEGER, INTENT(in) :: n, src,tag,comm
  COMPLEX(real_8), INTENT(out) , DIMENSION(1,1,*) :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_recv_complex8_r3'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+n*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_recv(DATA,n,MPI_DOUBLE_COMPLEX,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_recv_complex8_r3
SUBROUTINE mp_send_complex8_r4(DATA,n, dest,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments 
  INTEGER, INTENT(in) :: n, dest,tag,comm
  COMPLEX(real_8), INTENT(in) , DIMENSION(1,1,1,*) :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_isend_complex8_r4'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+n*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_send(DATA,n,MPI_DOUBLE_COMPLEX,dest,tag,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_send_complex8_r4

SUBROUTINE mp_recv_complex8_r4(DATA,n, src,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments
  INTEGER, INTENT(in) :: n, src,tag,comm
  COMPLEX(real_8), INTENT(out) , DIMENSION(1,1,1,*) :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_recv_complex8_r4'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+n*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_recv(DATA,n,MPI_DOUBLE_COMPLEX,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_recv_complex8_r4
SUBROUTINE mp_send_logical_r0(DATA, dest,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments 
  INTEGER, INTENT(in) ::  dest,tag,comm
  LOGICAL, INTENT(in)  :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_isend_logical_r0'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+1*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_send(DATA,1,MPI_LOGICAL,dest,tag,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_send_logical_r0

SUBROUTINE mp_recv_logical_r0(DATA, src,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments
  INTEGER, INTENT(in) ::  src,tag,comm
  LOGICAL, INTENT(out)  :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_recv_logical_r0'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+1*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_recv(DATA,1,MPI_LOGICAL,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_recv_logical_r0
SUBROUTINE mp_send_logical_r1(DATA,n, dest,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments 
  INTEGER, INTENT(in) :: n, dest,tag,comm
  LOGICAL, INTENT(in) , DIMENSION(*) :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_isend_logical_r1'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+n*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_send(DATA,n,MPI_LOGICAL,dest,tag,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_send_logical_r1

SUBROUTINE mp_recv_logical_r1(DATA,n, src,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments
  INTEGER, INTENT(in) :: n, src,tag,comm
  LOGICAL, INTENT(out) , DIMENSION(*) :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_recv_logical_r1'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+n*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_recv(DATA,n,MPI_LOGICAL,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_recv_logical_r1
SUBROUTINE mp_send_logical_r2(DATA,n, dest,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments 
  INTEGER, INTENT(in) :: n, dest,tag,comm
  LOGICAL, INTENT(in) , DIMENSION(1,*) :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_isend_logical_r2'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+n*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_send(DATA,n,MPI_LOGICAL,dest,tag,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_send_logical_r2

SUBROUTINE mp_recv_logical_r2(DATA,n, src,tag,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments
  INTEGER, INTENT(in) :: n, src,tag,comm
  LOGICAL, INTENT(out) , DIMENSION(1,*) :: DATA
  ! Variables
  INTEGER :: ierr
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_recv_logical_r2'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+n*{TYPE_SIZE_IN_BYTES}
  !tim1=timef()
#ifdef __PARALLEL
  CALL mpi_recv(DATA,n,MPI_LOGICAL,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_recv_logical_r2
