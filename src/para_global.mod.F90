MODULE para_global
  IMPLICIT NONE

  PRIVATE

  ! allows MPI to allocate memory for reduction operations 
  LOGICAL, PUBLIC :: para_use_mpi_in_place=.FALSE.

  ! maximum buffer size for reduction operations (require an allocation)
  INTEGER, PUBLIC :: para_buff_size=2**16

  ! stack buffer size for reduction operations
  INTEGER, PUBLIC :: para_stack_buff_size=2**8

END MODULE para_global
