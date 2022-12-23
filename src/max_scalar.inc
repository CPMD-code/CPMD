SUBROUTINE mp_max_int4_r0(DATA,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  IMPLICIT NONE
  ! Arguments
  INTEGER, INTENT(in) :: comm
  INTEGER(int_4), INTENT(inout) :: DATA
  ! Variables
  CHARACTER(*),PARAMETER :: procedureN='mp_max_int4_r0'
#ifdef __PARALLEL
  INTEGER :: ierr
  INTEGER(int_4) :: stack_buff
  CALL mpi_allreduce(DATA,stack_buff,1,MPI_INTEGER,&
         MPI_MAX,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  DATA=stack_buff
#endif
END SUBROUTINE mp_max_int4_r0
SUBROUTINE mp_max_int8_r0(DATA,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  IMPLICIT NONE
  ! Arguments
  INTEGER, INTENT(in) :: comm
  INTEGER(int_8), INTENT(inout) :: DATA
  ! Variables
  CHARACTER(*),PARAMETER :: procedureN='mp_max_int8_r0'
#ifdef __PARALLEL
  INTEGER :: ierr
  INTEGER(int_8) :: stack_buff
  CALL mpi_allreduce(DATA,stack_buff,1,MPI_INTEGER8,&
         MPI_MAX,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  DATA=stack_buff
#endif
END SUBROUTINE mp_max_int8_r0
SUBROUTINE mp_max_real8_r0(DATA,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  IMPLICIT NONE
  ! Arguments
  INTEGER, INTENT(in) :: comm
  REAL(real_8), INTENT(inout) :: DATA
  ! Variables
  CHARACTER(*),PARAMETER :: procedureN='mp_max_real8_r0'
#ifdef __PARALLEL
  INTEGER :: ierr
  REAL(real_8) :: stack_buff
  CALL mpi_allreduce(DATA,stack_buff,1,MPI_DOUBLE_PRECISION,&
         MPI_MAX,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  DATA=stack_buff
#endif
END SUBROUTINE mp_max_real8_r0
SUBROUTINE mp_max_complex8_r0(DATA,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  IMPLICIT NONE
  ! Arguments
  INTEGER, INTENT(in) :: comm
  COMPLEX(real_8), INTENT(inout) :: DATA
  ! Variables
  CHARACTER(*),PARAMETER :: procedureN='mp_max_complex8_r0'
#ifdef __PARALLEL
  INTEGER :: ierr
  COMPLEX(real_8) :: stack_buff
  CALL mpi_allreduce(DATA,stack_buff,1,MPI_DOUBLE_COMPLEX,&
         MPI_MAX,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  DATA=stack_buff
#endif
END SUBROUTINE mp_max_complex8_r0
