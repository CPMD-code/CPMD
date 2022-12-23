SUBROUTINE mp_bcast_char_r0(DATA,  root,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  ! Arguments
  INTEGER, INTENT(in) ::  root,comm
  CHARACTER(*), INTENT(inout)  :: DATA
  ! Variables
  INTEGER :: ierr
  CHARACTER(*),PARAMETER :: procedureN='mp_bcast_char_r0'
#ifdef __PARALLEL
  CALL mpi_bcast(DATA,LEN(DATA),MPI_CHARACTER,root,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
END SUBROUTINE mp_bcast_char_r0
SUBROUTINE mp_bcast_char_r1(DATA, n, root,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  ! Arguments
  INTEGER, INTENT(in) :: n, root,comm
  CHARACTER(*), INTENT(inout) , DIMENSION(*) :: DATA
  ! Variables
  INTEGER :: ierr
  CHARACTER(*),PARAMETER :: procedureN='mp_bcast_char_r1'
#ifdef __PARALLEL
  CALL mpi_bcast(DATA,LEN(DATA)*n,MPI_CHARACTER,root,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
END SUBROUTINE mp_bcast_char_r1
SUBROUTINE mp_bcast_int4_r0(DATA,  root,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  ! Arguments
  INTEGER, INTENT(in) ::  root,comm
  INTEGER(int_4), INTENT(inout)  :: DATA
  ! Variables
  INTEGER :: ierr
  CHARACTER(*),PARAMETER :: procedureN='mp_bcast_int4_r0'
#ifdef __PARALLEL
  CALL mpi_bcast(DATA,1,MPI_INTEGER,root,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
END SUBROUTINE mp_bcast_int4_r0
SUBROUTINE mp_bcast_int4_r1(DATA, n, root,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  ! Arguments
  INTEGER, INTENT(in) :: n, root,comm
  INTEGER(int_4), INTENT(inout) , DIMENSION(*) :: DATA
  ! Variables
  INTEGER :: ierr
  CHARACTER(*),PARAMETER :: procedureN='mp_bcast_int4_r1'
#ifdef __PARALLEL
  CALL mpi_bcast(DATA,n,MPI_INTEGER,root,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
END SUBROUTINE mp_bcast_int4_r1
SUBROUTINE mp_bcast_int4_r2(DATA, n, root,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  ! Arguments
  INTEGER, INTENT(in) :: n, root,comm
  INTEGER(int_4), INTENT(inout) , DIMENSION(1,*) :: DATA
  ! Variables
  INTEGER :: ierr
  CHARACTER(*),PARAMETER :: procedureN='mp_bcast_int4_r2'
#ifdef __PARALLEL
  CALL mpi_bcast(DATA,n,MPI_INTEGER,root,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
END SUBROUTINE mp_bcast_int4_r2
SUBROUTINE mp_bcast_int8_r0(DATA,  root,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  ! Arguments
  INTEGER, INTENT(in) ::  root,comm
  INTEGER(int_8), INTENT(inout)  :: DATA
  ! Variables
  INTEGER :: ierr
  CHARACTER(*),PARAMETER :: procedureN='mp_bcast_int8_r0'
#ifdef __PARALLEL
  CALL mpi_bcast(DATA,1,MPI_INTEGER8,root,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
END SUBROUTINE mp_bcast_int8_r0
SUBROUTINE mp_bcast_int8_r1(DATA, n, root,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  ! Arguments
  INTEGER, INTENT(in) :: n, root,comm
  INTEGER(int_8), INTENT(inout) , DIMENSION(*) :: DATA
  ! Variables
  INTEGER :: ierr
  CHARACTER(*),PARAMETER :: procedureN='mp_bcast_int8_r1'
#ifdef __PARALLEL
  CALL mpi_bcast(DATA,n,MPI_INTEGER8,root,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
END SUBROUTINE mp_bcast_int8_r1
SUBROUTINE mp_bcast_real8_r0(DATA,  root,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  ! Arguments
  INTEGER, INTENT(in) ::  root,comm
  REAL(real_8), INTENT(inout)  :: DATA
  ! Variables
  INTEGER :: ierr
  CHARACTER(*),PARAMETER :: procedureN='mp_bcast_real8_r0'
#ifdef __PARALLEL
  CALL mpi_bcast(DATA,1,MPI_DOUBLE_PRECISION,root,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
END SUBROUTINE mp_bcast_real8_r0
SUBROUTINE mp_bcast_real8_r1(DATA, n, root,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  ! Arguments
  INTEGER, INTENT(in) :: n, root,comm
  REAL(real_8), INTENT(inout) , DIMENSION(*) :: DATA
  ! Variables
  INTEGER :: ierr
  CHARACTER(*),PARAMETER :: procedureN='mp_bcast_real8_r1'
#ifdef __PARALLEL
  CALL mpi_bcast(DATA,n,MPI_DOUBLE_PRECISION,root,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
END SUBROUTINE mp_bcast_real8_r1
SUBROUTINE mp_bcast_real8_r2(DATA, n, root,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  ! Arguments
  INTEGER, INTENT(in) :: n, root,comm
  REAL(real_8), INTENT(inout) , DIMENSION(1,*) :: DATA
  ! Variables
  INTEGER :: ierr
  CHARACTER(*),PARAMETER :: procedureN='mp_bcast_real8_r2'
#ifdef __PARALLEL
  CALL mpi_bcast(DATA,n,MPI_DOUBLE_PRECISION,root,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
END SUBROUTINE mp_bcast_real8_r2
SUBROUTINE mp_bcast_real8_r3(DATA, n, root,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  ! Arguments
  INTEGER, INTENT(in) :: n, root,comm
  REAL(real_8), INTENT(inout) , DIMENSION(1,1,*) :: DATA
  ! Variables
  INTEGER :: ierr
  CHARACTER(*),PARAMETER :: procedureN='mp_bcast_real8_r3'
#ifdef __PARALLEL
  CALL mpi_bcast(DATA,n,MPI_DOUBLE_PRECISION,root,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
END SUBROUTINE mp_bcast_real8_r3
SUBROUTINE mp_bcast_real8_r4(DATA, n, root,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  ! Arguments
  INTEGER, INTENT(in) :: n, root,comm
  REAL(real_8), INTENT(inout) , DIMENSION(1,1,1,*) :: DATA
  ! Variables
  INTEGER :: ierr
  CHARACTER(*),PARAMETER :: procedureN='mp_bcast_real8_r4'
#ifdef __PARALLEL
  CALL mpi_bcast(DATA,n,MPI_DOUBLE_PRECISION,root,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
END SUBROUTINE mp_bcast_real8_r4
SUBROUTINE mp_bcast_real8_r5(DATA, n, root,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  ! Arguments
  INTEGER, INTENT(in) :: n, root,comm
  REAL(real_8), INTENT(inout) , DIMENSION(1,1,1,1,*) :: DATA
  ! Variables
  INTEGER :: ierr
  CHARACTER(*),PARAMETER :: procedureN='mp_bcast_real8_r5'
#ifdef __PARALLEL
  CALL mpi_bcast(DATA,n,MPI_DOUBLE_PRECISION,root,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
END SUBROUTINE mp_bcast_real8_r5
SUBROUTINE mp_bcast_complex8_r0(DATA,  root,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  ! Arguments
  INTEGER, INTENT(in) ::  root,comm
  COMPLEX(real_8), INTENT(inout)  :: DATA
  ! Variables
  INTEGER :: ierr
  CHARACTER(*),PARAMETER :: procedureN='mp_bcast_complex8_r0'
#ifdef __PARALLEL
  CALL mpi_bcast(DATA,1,MPI_DOUBLE_COMPLEX,root,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
END SUBROUTINE mp_bcast_complex8_r0
SUBROUTINE mp_bcast_complex8_r1(DATA, n, root,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  ! Arguments
  INTEGER, INTENT(in) :: n, root,comm
  COMPLEX(real_8), INTENT(inout) , DIMENSION(*) :: DATA
  ! Variables
  INTEGER :: ierr
  CHARACTER(*),PARAMETER :: procedureN='mp_bcast_complex8_r1'
#ifdef __PARALLEL
  CALL mpi_bcast(DATA,n,MPI_DOUBLE_COMPLEX,root,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
END SUBROUTINE mp_bcast_complex8_r1
SUBROUTINE mp_bcast_complex8_r2(DATA, n, root,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  ! Arguments
  INTEGER, INTENT(in) :: n, root,comm
  COMPLEX(real_8), INTENT(inout) , DIMENSION(1,*) :: DATA
  ! Variables
  INTEGER :: ierr
  CHARACTER(*),PARAMETER :: procedureN='mp_bcast_complex8_r2'
#ifdef __PARALLEL
  CALL mpi_bcast(DATA,n,MPI_DOUBLE_COMPLEX,root,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
END SUBROUTINE mp_bcast_complex8_r2
SUBROUTINE mp_bcast_complex8_r3(DATA, n, root,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  ! Arguments
  INTEGER, INTENT(in) :: n, root,comm
  COMPLEX(real_8), INTENT(inout) , DIMENSION(1,1,*) :: DATA
  ! Variables
  INTEGER :: ierr
  CHARACTER(*),PARAMETER :: procedureN='mp_bcast_complex8_r3'
#ifdef __PARALLEL
  CALL mpi_bcast(DATA,n,MPI_DOUBLE_COMPLEX,root,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
END SUBROUTINE mp_bcast_complex8_r3
SUBROUTINE mp_bcast_complex8_r4(DATA, n, root,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  ! Arguments
  INTEGER, INTENT(in) :: n, root,comm
  COMPLEX(real_8), INTENT(inout) , DIMENSION(1,1,1,*) :: DATA
  ! Variables
  INTEGER :: ierr
  CHARACTER(*),PARAMETER :: procedureN='mp_bcast_complex8_r4'
#ifdef __PARALLEL
  CALL mpi_bcast(DATA,n,MPI_DOUBLE_COMPLEX,root,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
END SUBROUTINE mp_bcast_complex8_r4
SUBROUTINE mp_bcast_logical_r0(DATA,  root,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  ! Arguments
  INTEGER, INTENT(in) ::  root,comm
  LOGICAL, INTENT(inout)  :: DATA
  ! Variables
  INTEGER :: ierr
  CHARACTER(*),PARAMETER :: procedureN='mp_bcast_logical_r0'
#ifdef __PARALLEL
  CALL mpi_bcast(DATA,1,MPI_LOGICAL,root,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
END SUBROUTINE mp_bcast_logical_r0
SUBROUTINE mp_bcast_logical_r1(DATA, n, root,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  ! Arguments
  INTEGER, INTENT(in) :: n, root,comm
  LOGICAL, INTENT(inout) , DIMENSION(*) :: DATA
  ! Variables
  INTEGER :: ierr
  CHARACTER(*),PARAMETER :: procedureN='mp_bcast_logical_r1'
#ifdef __PARALLEL
  CALL mpi_bcast(DATA,n,MPI_LOGICAL,root,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
END SUBROUTINE mp_bcast_logical_r1
SUBROUTINE mp_bcast_logical_r2(DATA, n, root,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  ! Arguments
  INTEGER, INTENT(in) :: n, root,comm
  LOGICAL, INTENT(inout) , DIMENSION(1,*) :: DATA
  ! Variables
  INTEGER :: ierr
  CHARACTER(*),PARAMETER :: procedureN='mp_bcast_logical_r2'
#ifdef __PARALLEL
  CALL mpi_bcast(DATA,n,MPI_LOGICAL,root,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
END SUBROUTINE mp_bcast_logical_r2
