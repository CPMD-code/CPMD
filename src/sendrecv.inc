SUBROUTINE mp_sendrecv_int4_r1(data_in,n_in,dest,data_out,n_out,src,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments
  INTEGER, INTENT(in) :: n_in,n_out,src,dest,comm
  INTEGER(int_4), INTENT(in) , DIMENSION(*) :: data_in
  INTEGER(int_4), INTENT(out) , DIMENSION(*) :: data_out
  ! Variables
  INTEGER :: ierr,tag_in,tag_out
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_sendrecv_int4_r1'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+scount
  !tim1=timef()
#ifdef __PARALLEL
  tag_in=1;tag_out=1
  CALL mpi_sendrecv(data_in,n_in,MPI_INTEGER,dest,tag_in,&
       data_out,n_out,MPI_INTEGER,src,tag_out,comm,&
       MPI_STATUS_IGNORE,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
  CALL mp_mpi_error_assert(1,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_sendrecv_int4_r1
SUBROUTINE mp_sendrecv_int4_r2(data_in,n_in,dest,data_out,n_out,src,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments
  INTEGER, INTENT(in) :: n_in,n_out,src,dest,comm
  INTEGER(int_4), INTENT(in) , DIMENSION(1,*) :: data_in
  INTEGER(int_4), INTENT(out) , DIMENSION(1,*) :: data_out
  ! Variables
  INTEGER :: ierr,tag_in,tag_out
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_sendrecv_int4_r2'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+scount
  !tim1=timef()
#ifdef __PARALLEL
  tag_in=1;tag_out=1
  CALL mpi_sendrecv(data_in,n_in,MPI_INTEGER,dest,tag_in,&
       data_out,n_out,MPI_INTEGER,src,tag_out,comm,&
       MPI_STATUS_IGNORE,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
  CALL mp_mpi_error_assert(1,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_sendrecv_int4_r2
SUBROUTINE mp_sendrecv_int8_r1(data_in,n_in,dest,data_out,n_out,src,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments
  INTEGER, INTENT(in) :: n_in,n_out,src,dest,comm
  INTEGER(int_8), INTENT(in) , DIMENSION(*) :: data_in
  INTEGER(int_8), INTENT(out) , DIMENSION(*) :: data_out
  ! Variables
  INTEGER :: ierr,tag_in,tag_out
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_sendrecv_int8_r1'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+scount
  !tim1=timef()
#ifdef __PARALLEL
  tag_in=1;tag_out=1
  CALL mpi_sendrecv(data_in,n_in,MPI_INTEGER8,dest,tag_in,&
       data_out,n_out,MPI_INTEGER8,src,tag_out,comm,&
       MPI_STATUS_IGNORE,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
  CALL mp_mpi_error_assert(1,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_sendrecv_int8_r1
SUBROUTINE mp_sendrecv_real8_r1(data_in,n_in,dest,data_out,n_out,src,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments
  INTEGER, INTENT(in) :: n_in,n_out,src,dest,comm
  REAL(real_8), INTENT(in) , DIMENSION(*) :: data_in
  REAL(real_8), INTENT(out) , DIMENSION(*) :: data_out
  ! Variables
  INTEGER :: ierr,tag_in,tag_out
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_sendrecv_real8_r1'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+scount
  !tim1=timef()
#ifdef __PARALLEL
  tag_in=1;tag_out=1
  CALL mpi_sendrecv(data_in,n_in,MPI_DOUBLE_PRECISION,dest,tag_in,&
       data_out,n_out,MPI_DOUBLE_PRECISION,src,tag_out,comm,&
       MPI_STATUS_IGNORE,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
  CALL mp_mpi_error_assert(1,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_sendrecv_real8_r1
SUBROUTINE mp_sendrecv_real8_r2(data_in,n_in,dest,data_out,n_out,src,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments
  INTEGER, INTENT(in) :: n_in,n_out,src,dest,comm
  REAL(real_8), INTENT(in) , DIMENSION(1,*) :: data_in
  REAL(real_8), INTENT(out) , DIMENSION(1,*) :: data_out
  ! Variables
  INTEGER :: ierr,tag_in,tag_out
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_sendrecv_real8_r2'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+scount
  !tim1=timef()
#ifdef __PARALLEL
  tag_in=1;tag_out=1
  CALL mpi_sendrecv(data_in,n_in,MPI_DOUBLE_PRECISION,dest,tag_in,&
       data_out,n_out,MPI_DOUBLE_PRECISION,src,tag_out,comm,&
       MPI_STATUS_IGNORE,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
  CALL mp_mpi_error_assert(1,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_sendrecv_real8_r2
SUBROUTINE mp_sendrecv_real8_r3(data_in,n_in,dest,data_out,n_out,src,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments
  INTEGER, INTENT(in) :: n_in,n_out,src,dest,comm
  REAL(real_8), INTENT(in) , DIMENSION(1,1,*) :: data_in
  REAL(real_8), INTENT(out) , DIMENSION(1,1,*) :: data_out
  ! Variables
  INTEGER :: ierr,tag_in,tag_out
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_sendrecv_real8_r3'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+scount
  !tim1=timef()
#ifdef __PARALLEL
  tag_in=1;tag_out=1
  CALL mpi_sendrecv(data_in,n_in,MPI_DOUBLE_PRECISION,dest,tag_in,&
       data_out,n_out,MPI_DOUBLE_PRECISION,src,tag_out,comm,&
       MPI_STATUS_IGNORE,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
  CALL mp_mpi_error_assert(1,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_sendrecv_real8_r3
SUBROUTINE mp_sendrecv_real8_r4(data_in,n_in,dest,data_out,n_out,src,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments
  INTEGER, INTENT(in) :: n_in,n_out,src,dest,comm
  REAL(real_8), INTENT(in) , DIMENSION(1,1,1,*) :: data_in
  REAL(real_8), INTENT(out) , DIMENSION(1,1,1,*) :: data_out
  ! Variables
  INTEGER :: ierr,tag_in,tag_out
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_sendrecv_real8_r4'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+scount
  !tim1=timef()
#ifdef __PARALLEL
  tag_in=1;tag_out=1
  CALL mpi_sendrecv(data_in,n_in,MPI_DOUBLE_PRECISION,dest,tag_in,&
       data_out,n_out,MPI_DOUBLE_PRECISION,src,tag_out,comm,&
       MPI_STATUS_IGNORE,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
  CALL mp_mpi_error_assert(1,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_sendrecv_real8_r4
SUBROUTINE mp_sendrecv_complex8_r1(data_in,n_in,dest,data_out,n_out,src,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments
  INTEGER, INTENT(in) :: n_in,n_out,src,dest,comm
  COMPLEX(real_8), INTENT(in) , DIMENSION(*) :: data_in
  COMPLEX(real_8), INTENT(out) , DIMENSION(*) :: data_out
  ! Variables
  INTEGER :: ierr,tag_in,tag_out
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_sendrecv_complex8_r1'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+scount
  !tim1=timef()
#ifdef __PARALLEL
  tag_in=1;tag_out=1
  CALL mpi_sendrecv(data_in,n_in,MPI_DOUBLE_COMPLEX,dest,tag_in,&
       data_out,n_out,MPI_DOUBLE_COMPLEX,src,tag_out,comm,&
       MPI_STATUS_IGNORE,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
  CALL mp_mpi_error_assert(1,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_sendrecv_complex8_r1
SUBROUTINE mp_sendrecv_complex8_r2(data_in,n_in,dest,data_out,n_out,src,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments
  INTEGER, INTENT(in) :: n_in,n_out,src,dest,comm
  COMPLEX(real_8), INTENT(in) , DIMENSION(1,*) :: data_in
  COMPLEX(real_8), INTENT(out) , DIMENSION(1,*) :: data_out
  ! Variables
  INTEGER :: ierr,tag_in,tag_out
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_sendrecv_complex8_r2'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+scount
  !tim1=timef()
#ifdef __PARALLEL
  tag_in=1;tag_out=1
  CALL mpi_sendrecv(data_in,n_in,MPI_DOUBLE_COMPLEX,dest,tag_in,&
       data_out,n_out,MPI_DOUBLE_COMPLEX,src,tag_out,comm,&
       MPI_STATUS_IGNORE,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
  CALL mp_mpi_error_assert(1,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_sendrecv_complex8_r2
SUBROUTINE mp_sendrecv_complex8_r3(data_in,n_in,dest,data_out,n_out,src,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments
  INTEGER, INTENT(in) :: n_in,n_out,src,dest,comm
  COMPLEX(real_8), INTENT(in) , DIMENSION(1,1,*) :: data_in
  COMPLEX(real_8), INTENT(out) , DIMENSION(1,1,*) :: data_out
  ! Variables
  INTEGER :: ierr,tag_in,tag_out
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_sendrecv_complex8_r3'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+scount
  !tim1=timef()
#ifdef __PARALLEL
  tag_in=1;tag_out=1
  CALL mpi_sendrecv(data_in,n_in,MPI_DOUBLE_COMPLEX,dest,tag_in,&
       data_out,n_out,MPI_DOUBLE_COMPLEX,src,tag_out,comm,&
       MPI_STATUS_IGNORE,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
  CALL mp_mpi_error_assert(1,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_sendrecv_complex8_r3
SUBROUTINE mp_sendrecv_complex8_r4(data_in,n_in,dest,data_out,n_out,src,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE pstat
  ! Arguments
  INTEGER, INTENT(in) :: n_in,n_out,src,dest,comm
  COMPLEX(real_8), INTENT(in) , DIMENSION(1,1,1,*) :: data_in
  COMPLEX(real_8), INTENT(out) , DIMENSION(1,1,1,*) :: data_out
  ! Variables
  INTEGER :: ierr,tag_in,tag_out
  REAL(real_8) :: tim1,tim2
  REAL(real_8), EXTERNAL :: timef
  CHARACTER(*),PARAMETER::procedureN='mp_sendrecv_complex8_r4'
  !cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  !cmlen(ipar_send)=cmlen(ipar_send)+scount
  !tim1=timef()
#ifdef __PARALLEL
  tag_in=1;tag_out=1
  CALL mpi_sendrecv(data_in,n_in,MPI_DOUBLE_COMPLEX,dest,tag_in,&
       data_out,n_out,MPI_DOUBLE_COMPLEX,src,tag_out,comm,&
       MPI_STATUS_IGNORE,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
  CALL mp_mpi_error_assert(1,procedureN,__LINE__,__FILE__)
#endif
  !tim2=timef()
  !cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
END SUBROUTINE mp_sendrecv_complex8_r4
