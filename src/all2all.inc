subroutine mp_all2all_int4_r1(data_in, data_out, n, comm )
#ifdef __PARALLEL
  use mpi
#endif
  ! Arguments
  integer, intent(in) :: n, comm
  integer(int_4), intent(inout) , dimension(*) :: data_out, data_in
  ! Variables
  integer :: ierr
  character(*),parameter :: procedureN='mp_all2all_int4_r1'
#ifdef __PARALLEL

  CALL mpi_alltoall(data_in,n,MPI_INTEGER,data_out,n,MPI_INTEGER,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
  data_out(1:n)=data_in(1:n)
#endif
end subroutine mp_all2all_int4_r1
subroutine mp_all2all_int8_r1(data_in, data_out, n, comm )
#ifdef __PARALLEL
  use mpi
#endif
  ! Arguments
  integer, intent(in) :: n, comm
  integer(int_8), intent(inout) , dimension(*) :: data_out, data_in
  ! Variables
  integer :: ierr
  character(*),parameter :: procedureN='mp_all2all_int8_r1'
#ifdef __PARALLEL

  CALL mpi_alltoall(data_in,n,MPI_INTEGER8,data_out,n,MPI_INTEGER8,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
  data_out(1:n)=data_in(1:n)
#endif
end subroutine mp_all2all_int8_r1
subroutine mp_all2all_real8_r1(data_in, data_out, n, comm )
#ifdef __PARALLEL
  use mpi
#endif
  ! Arguments
  integer, intent(in) :: n, comm
  real(real_8), intent(inout) , dimension(*) :: data_out, data_in
  ! Variables
  integer :: ierr
  character(*),parameter :: procedureN='mp_all2all_real8_r1'
#ifdef __PARALLEL

  CALL mpi_alltoall(data_in,n,MPI_DOUBLE_PRECISION,data_out,n,MPI_DOUBLE_PRECISION,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
  data_out(1:n)=data_in(1:n)
#endif
end subroutine mp_all2all_real8_r1
subroutine mp_all2all_complex4_r1(data_in, data_out, n, comm )
#ifdef __PARALLEL
  use mpi
#endif
  ! Arguments
  integer, intent(in) :: n, comm
  complex(real_4), intent(inout) , dimension(*) :: data_out, data_in
  ! Variables
  integer :: ierr
  character(*),parameter :: procedureN='mp_all2all_complex4_r1'
#ifdef __PARALLEL

  CALL mpi_alltoall(data_in,n,MPI_COMPLEX,data_out,n,MPI_COMPLEX,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
  data_out(1:n)=data_in(1:n)
#endif
end subroutine mp_all2all_complex4_r1
subroutine mp_all2all_complex8_r1(data_in, data_out, n, comm )
#ifdef __PARALLEL
  use mpi
#endif
  ! Arguments
  integer, intent(in) :: n, comm
  complex(real_8), intent(inout) , dimension(*) :: data_out, data_in
  ! Variables
  integer :: ierr
  character(*),parameter :: procedureN='mp_all2all_complex8_r1'
#ifdef __PARALLEL

  CALL mpi_alltoall(data_in,n,MPI_DOUBLE_COMPLEX,data_out,n,MPI_DOUBLE_COMPLEX,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
  data_out(1:n)=data_in(1:n)
#endif
end subroutine mp_all2all_complex8_r1
