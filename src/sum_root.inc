SUBROUTINE mp_sum_root_int4_r0(data_in,data_out, root, comm)
#ifdef __PARALLEL
  USE mpi
#endif
  ! Arguments
  INTEGER, INTENT(in) ::  root, comm
  INTEGER(int_4), INTENT(in)  :: data_in
  INTEGER(int_4), INTENT(inout)  :: data_out
  ! Variables
  CHARACTER(*),PARAMETER :: procedureN='mp_sum_root_int4_r0'
#ifdef __PARALLEL
  INTEGER :: ierr
  CALL mpi_reduce(data_in,data_out,1,MPI_INTEGER,&
       MPI_SUM,root,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
  data_out =data_in 
#endif
END SUBROUTINE mp_sum_root_int4_r0
SUBROUTINE mp_sum_root_int4_r1(data_in,data_out,n, root, comm)
#ifdef __PARALLEL
  USE mpi
#endif
  ! Arguments
  INTEGER, INTENT(in) :: n, root, comm
  INTEGER(int_4), INTENT(in) , DIMENSION(*) :: data_in
  INTEGER(int_4), INTENT(inout) , DIMENSION(*) :: data_out
  ! Variables
  CHARACTER(*),PARAMETER :: procedureN='mp_sum_root_int4_r1'
#ifdef __PARALLEL
  INTEGER :: ierr
  CALL mpi_reduce(data_in,data_out,n,MPI_INTEGER,&
       MPI_SUM,root,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
  data_out (1:n)=data_in (1:n)
#endif
END SUBROUTINE mp_sum_root_int4_r1
SUBROUTINE mp_sum_root_int4_r2(data_in,data_out,n, root, comm)
#ifdef __PARALLEL
  USE mpi
#endif
  ! Arguments
  INTEGER, INTENT(in) :: n, root, comm
  INTEGER(int_4), INTENT(in) , DIMENSION(1,*) :: data_in
  INTEGER(int_4), INTENT(inout) , DIMENSION(1,*) :: data_out
  ! Variables
  CHARACTER(*),PARAMETER :: procedureN='mp_sum_root_int4_r2'
#ifdef __PARALLEL
  INTEGER :: ierr
  CALL mpi_reduce(data_in,data_out,n,MPI_INTEGER,&
       MPI_SUM,root,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
  data_out (1,1:n)=data_in (1,1:n)
#endif
END SUBROUTINE mp_sum_root_int4_r2
SUBROUTINE mp_sum_root_real8_r0(data_in,data_out, root, comm)
#ifdef __PARALLEL
  USE mpi
#endif
  ! Arguments
  INTEGER, INTENT(in) ::  root, comm
  REAL(real_8), INTENT(in)  :: data_in
  REAL(real_8), INTENT(inout)  :: data_out
  ! Variables
  CHARACTER(*),PARAMETER :: procedureN='mp_sum_root_real8_r0'
#ifdef __PARALLEL
  INTEGER :: ierr
  CALL mpi_reduce(data_in,data_out,1,MPI_DOUBLE_PRECISION,&
       MPI_SUM,root,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
  data_out =data_in 
#endif
END SUBROUTINE mp_sum_root_real8_r0
SUBROUTINE mp_sum_root_real8_r1(data_in,data_out,n, root, comm)
#ifdef __PARALLEL
  USE mpi
#endif
  ! Arguments
  INTEGER, INTENT(in) :: n, root, comm
  REAL(real_8), INTENT(in) , DIMENSION(*) :: data_in
  REAL(real_8), INTENT(inout) , DIMENSION(*) :: data_out
  ! Variables
  CHARACTER(*),PARAMETER :: procedureN='mp_sum_root_real8_r1'
#ifdef __PARALLEL
  INTEGER :: ierr
  CALL mpi_reduce(data_in,data_out,n,MPI_DOUBLE_PRECISION,&
       MPI_SUM,root,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
  data_out (1:n)=data_in (1:n)
#endif
END SUBROUTINE mp_sum_root_real8_r1
SUBROUTINE mp_sum_root_real8_r2(data_in,data_out,n, root, comm)
#ifdef __PARALLEL
  USE mpi
#endif
  ! Arguments
  INTEGER, INTENT(in) :: n, root, comm
  REAL(real_8), INTENT(in) , DIMENSION(1,*) :: data_in
  REAL(real_8), INTENT(inout) , DIMENSION(1,*) :: data_out
  ! Variables
  CHARACTER(*),PARAMETER :: procedureN='mp_sum_root_real8_r2'
#ifdef __PARALLEL
  INTEGER :: ierr
  CALL mpi_reduce(data_in,data_out,n,MPI_DOUBLE_PRECISION,&
       MPI_SUM,root,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
  data_out (1,1:n)=data_in (1,1:n)
#endif
END SUBROUTINE mp_sum_root_real8_r2
SUBROUTINE mp_sum_root_complex8_r0(data_in,data_out, root, comm)
#ifdef __PARALLEL
  USE mpi
#endif
  ! Arguments
  INTEGER, INTENT(in) ::  root, comm
  COMPLEX(real_8), INTENT(in)  :: data_in
  COMPLEX(real_8), INTENT(inout)  :: data_out
  ! Variables
  CHARACTER(*),PARAMETER :: procedureN='mp_sum_root_complex8_r0'
#ifdef __PARALLEL
  INTEGER :: ierr
  CALL mpi_reduce(data_in,data_out,1,MPI_DOUBLE_COMPLEX,&
       MPI_SUM,root,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
  data_out =data_in 
#endif
END SUBROUTINE mp_sum_root_complex8_r0
SUBROUTINE mp_sum_root_complex8_r1(data_in,data_out,n, root, comm)
#ifdef __PARALLEL
  USE mpi
#endif
  ! Arguments
  INTEGER, INTENT(in) :: n, root, comm
  COMPLEX(real_8), INTENT(in) , DIMENSION(*) :: data_in
  COMPLEX(real_8), INTENT(inout) , DIMENSION(*) :: data_out
  ! Variables
  CHARACTER(*),PARAMETER :: procedureN='mp_sum_root_complex8_r1'
#ifdef __PARALLEL
  INTEGER :: ierr
  CALL mpi_reduce(data_in,data_out,n,MPI_DOUBLE_COMPLEX,&
       MPI_SUM,root,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
  data_out (1:n)=data_in (1:n)
#endif
END SUBROUTINE mp_sum_root_complex8_r1
SUBROUTINE mp_sum_root_complex8_r2(data_in,data_out,n, root, comm)
#ifdef __PARALLEL
  USE mpi
#endif
  ! Arguments
  INTEGER, INTENT(in) :: n, root, comm
  COMPLEX(real_8), INTENT(in) , DIMENSION(1,*) :: data_in
  COMPLEX(real_8), INTENT(inout) , DIMENSION(1,*) :: data_out
  ! Variables
  CHARACTER(*),PARAMETER :: procedureN='mp_sum_root_complex8_r2'
#ifdef __PARALLEL
  INTEGER :: ierr
  CALL mpi_reduce(data_in,data_out,n,MPI_DOUBLE_COMPLEX,&
       MPI_SUM,root,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
  data_out (1,1:n)=data_in (1,1:n)
#endif
END SUBROUTINE mp_sum_root_complex8_r2
