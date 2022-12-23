SUBROUTINE mp_max_int4_r1(DATA,n,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE para_global, ONLY : para_use_mpi_in_place
  USE para_global, ONLY : para_buff_size
  USE para_global, ONLY : para_stack_buff_size
  IMPLICIT NONE
  ! Arguments
  INTEGER, INTENT(in) :: n,comm
  INTEGER(int_4), INTENT(inout) , DIMENSION(*) :: DATA
  ! Variables
  CHARACTER(*),PARAMETER :: procedureN='mp_max_int4_r1'
#ifdef __PARALLEL
  INTEGER :: ierr,m,i,stat,buff_size
  INTEGER(int_4),DIMENSION(:),ALLOCATABLE :: buff
  INTEGER(int_4),DIMENSION(para_stack_buff_size) :: stack_buff
  IF (n<=para_stack_buff_size) THEN
    CALL mpi_allreduce(DATA,stack_buff,n,MPI_INTEGER,&
         MPI_MAX,comm,ierr)
    CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
    DATA(1:n)=stack_buff(1:n)
  ELSE
    IF (para_use_mpi_in_place) THEN
      CALL mpi_allreduce(MPI_IN_PLACE,DATA,n,MPI_INTEGER,&
           MPI_MAX,comm,ierr)
      CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
    ELSE
      buff_size=MIN(n,para_buff_size)
      ALLOCATE(buff(buff_size),stat=stat)
      IF (stat.ne.0) CALL stopgm(procedureN,'allocation problem',& 
           __LINE__,__FILE__)
      DO i=1,n,buff_size
        m=MIN(buff_size,n-i+1)
        CALL mpi_allreduce(DATA(i),buff,m,MPI_INTEGER,&
             MPI_MAX,comm,ierr)
        CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
        DATA(i:i+m-1)=buff(1:m)
      ENDDO
      DEALLOCATE(buff,stat=stat)
      IF (stat.ne.0) CALL stopgm(procedureN,'deallocation problem',&
           __LINE__,__FILE__)
    ENDIF
  ENDIF
#endif
END SUBROUTINE mp_max_int4_r1
SUBROUTINE mp_max_real8_r1(DATA,n,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE para_global, ONLY : para_use_mpi_in_place
  USE para_global, ONLY : para_buff_size
  USE para_global, ONLY : para_stack_buff_size
  IMPLICIT NONE
  ! Arguments
  INTEGER, INTENT(in) :: n,comm
  REAL(real_8), INTENT(inout) , DIMENSION(*) :: DATA
  ! Variables
  CHARACTER(*),PARAMETER :: procedureN='mp_max_real8_r1'
#ifdef __PARALLEL
  INTEGER :: ierr,m,i,stat,buff_size
  REAL(real_8),DIMENSION(:),ALLOCATABLE :: buff
  REAL(real_8),DIMENSION(para_stack_buff_size) :: stack_buff
  IF (n<=para_stack_buff_size) THEN
    CALL mpi_allreduce(DATA,stack_buff,n,MPI_DOUBLE_PRECISION,&
         MPI_MAX,comm,ierr)
    CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
    DATA(1:n)=stack_buff(1:n)
  ELSE
    IF (para_use_mpi_in_place) THEN
      CALL mpi_allreduce(MPI_IN_PLACE,DATA,n,MPI_DOUBLE_PRECISION,&
           MPI_MAX,comm,ierr)
      CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
    ELSE
      buff_size=MIN(n,para_buff_size)
      ALLOCATE(buff(buff_size),stat=stat)
      IF (stat.ne.0) CALL stopgm(procedureN,'allocation problem',& 
           __LINE__,__FILE__)
      DO i=1,n,buff_size
        m=MIN(buff_size,n-i+1)
        CALL mpi_allreduce(DATA(i),buff,m,MPI_DOUBLE_PRECISION,&
             MPI_MAX,comm,ierr)
        CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
        DATA(i:i+m-1)=buff(1:m)
      ENDDO
      DEALLOCATE(buff,stat=stat)
      IF (stat.ne.0) CALL stopgm(procedureN,'deallocation problem',&
           __LINE__,__FILE__)
    ENDIF
  ENDIF
#endif
END SUBROUTINE mp_max_real8_r1
SUBROUTINE mp_max_complex8_r1(DATA,n,comm)
#ifdef __PARALLEL
  USE mpi
#endif
  USE para_global, ONLY : para_use_mpi_in_place
  USE para_global, ONLY : para_buff_size
  USE para_global, ONLY : para_stack_buff_size
  IMPLICIT NONE
  ! Arguments
  INTEGER, INTENT(in) :: n,comm
  COMPLEX(real_8), INTENT(inout) , DIMENSION(*) :: DATA
  ! Variables
  CHARACTER(*),PARAMETER :: procedureN='mp_max_complex8_r1'
#ifdef __PARALLEL
  INTEGER :: ierr,m,i,stat,buff_size
  COMPLEX(real_8),DIMENSION(:),ALLOCATABLE :: buff
  COMPLEX(real_8),DIMENSION(para_stack_buff_size) :: stack_buff
  IF (n<=para_stack_buff_size) THEN
    CALL mpi_allreduce(DATA,stack_buff,n,MPI_DOUBLE_COMPLEX,&
         MPI_MAX,comm,ierr)
    CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
    DATA(1:n)=stack_buff(1:n)
  ELSE
    IF (para_use_mpi_in_place) THEN
      CALL mpi_allreduce(MPI_IN_PLACE,DATA,n,MPI_DOUBLE_COMPLEX,&
           MPI_MAX,comm,ierr)
      CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
    ELSE
      buff_size=MIN(n,para_buff_size)
      ALLOCATE(buff(buff_size),stat=stat)
      IF (stat.ne.0) CALL stopgm(procedureN,'allocation problem',& 
           __LINE__,__FILE__)
      DO i=1,n,buff_size
        m=MIN(buff_size,n-i+1)
        CALL mpi_allreduce(DATA(i),buff,m,MPI_DOUBLE_COMPLEX,&
             MPI_MAX,comm,ierr)
        CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
        DATA(i:i+m-1)=buff(1:m)
      ENDDO
      DEALLOCATE(buff,stat=stat)
      IF (stat.ne.0) CALL stopgm(procedureN,'deallocation problem',&
           __LINE__,__FILE__)
    ENDIF
  ENDIF
#endif
END SUBROUTINE mp_max_complex8_r1
