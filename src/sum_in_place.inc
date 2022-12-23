subroutine mp_sum_in_place_int4_r1(data,n,comm)
#ifdef __PARALLEL
  use mpi
#endif
  use para_global, only : para_use_mpi_in_place
  use para_global, only : para_buff_size
  use para_global, only : para_stack_buff_size
  implicit none
  ! Arguments
  integer, intent(in) :: n,comm
  integer(int_4), intent(inout) , dimension(*) :: data
  ! Variables
  character(*),parameter :: procedureN='mp_sum_in_place_int4_r1'
#ifdef __PARALLEL
  integer :: ierr,m,i,stat,buff_size
  integer(int_4),dimension(:),allocatable :: buff
  integer(int_4),dimension(para_stack_buff_size) :: stack_buff
  if (n<=para_stack_buff_size) then
    call mpi_allreduce(data,stack_buff,n,MPI_INTEGER,&
         MPI_SUM,comm,ierr)
    call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
    data(1:n)=stack_buff(1:n)
  else
    if (para_use_mpi_in_place) then
      call mpi_allreduce(MPI_IN_PLACE,data,n,MPI_INTEGER,&
           MPI_SUM,comm,ierr)
      call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
    else
      buff_size=min(n,para_buff_size)
      allocate(buff(buff_size),stat=stat)
      if (stat.ne.0) call stopgm(procedureN,'allocation problem',& 
           __LINE__,__FILE__)
      do i=1,n,buff_size
        m=min(buff_size,n-i+1)
        call mpi_allreduce(data(i),buff,m,MPI_INTEGER,&
             MPI_SUM,comm,ierr)
        call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
        data(i:i+m-1)=buff(1:m)
      enddo
      deallocate(buff,stat=stat)
      if (stat.ne.0) call stopgm(procedureN,'deallocation problem',&
           __LINE__,__FILE__)
    endif
  endif
#endif
end subroutine mp_sum_in_place_int4_r1
subroutine mp_sum_in_place_int4_r2(data,n,comm)
#ifdef __PARALLEL
  use mpi
#endif
  use para_global, only : para_use_mpi_in_place
  use para_global, only : para_buff_size
  use para_global, only : para_stack_buff_size
  implicit none
  ! Arguments
  integer, intent(in) :: n,comm
  integer(int_4), intent(inout) , dimension(1,*) :: data
  ! Variables
  character(*),parameter :: procedureN='mp_sum_in_place_int4_r2'
#ifdef __PARALLEL
  integer :: ierr,m,i,stat,buff_size
  integer(int_4),dimension(:),allocatable :: buff
  integer(int_4),dimension(para_stack_buff_size) :: stack_buff
  if (n<=para_stack_buff_size) then
    call mpi_allreduce(data,stack_buff,n,MPI_INTEGER,&
         MPI_SUM,comm,ierr)
    call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
    data(1,1:n)=stack_buff(1:n)
  else
    if (para_use_mpi_in_place) then
      call mpi_allreduce(MPI_IN_PLACE,data,n,MPI_INTEGER,&
           MPI_SUM,comm,ierr)
      call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
    else
      buff_size=min(n,para_buff_size)
      allocate(buff(buff_size),stat=stat)
      if (stat.ne.0) call stopgm(procedureN,'allocation problem',& 
           __LINE__,__FILE__)
      do i=1,n,buff_size
        m=min(buff_size,n-i+1)
        call mpi_allreduce(data(1,i),buff,m,MPI_INTEGER,&
             MPI_SUM,comm,ierr)
        call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
        data(1,i:i+m-1)=buff(1:m)
      enddo
      deallocate(buff,stat=stat)
      if (stat.ne.0) call stopgm(procedureN,'deallocation problem',&
           __LINE__,__FILE__)
    endif
  endif
#endif
end subroutine mp_sum_in_place_int4_r2
subroutine mp_sum_in_place_int8_r1(data,n,comm)
#ifdef __PARALLEL
  use mpi
#endif
  use para_global, only : para_use_mpi_in_place
  use para_global, only : para_buff_size
  use para_global, only : para_stack_buff_size
  implicit none
  ! Arguments
  integer, intent(in) :: n,comm
  integer(int_8), intent(inout) , dimension(*) :: data
  ! Variables
  character(*),parameter :: procedureN='mp_sum_in_place_int8_r1'
#ifdef __PARALLEL
  integer :: ierr,m,i,stat,buff_size
  integer(int_8),dimension(:),allocatable :: buff
  integer(int_8),dimension(para_stack_buff_size) :: stack_buff
  if (n<=para_stack_buff_size) then
    call mpi_allreduce(data,stack_buff,n,MPI_INTEGER8,&
         MPI_SUM,comm,ierr)
    call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
    data(1:n)=stack_buff(1:n)
  else
    if (para_use_mpi_in_place) then
      call mpi_allreduce(MPI_IN_PLACE,data,n,MPI_INTEGER8,&
           MPI_SUM,comm,ierr)
      call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
    else
      buff_size=min(n,para_buff_size)
      allocate(buff(buff_size),stat=stat)
      if (stat.ne.0) call stopgm(procedureN,'allocation problem',& 
           __LINE__,__FILE__)
      do i=1,n,buff_size
        m=min(buff_size,n-i+1)
        call mpi_allreduce(data(i),buff,m,MPI_INTEGER8,&
             MPI_SUM,comm,ierr)
        call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
        data(i:i+m-1)=buff(1:m)
      enddo
      deallocate(buff,stat=stat)
      if (stat.ne.0) call stopgm(procedureN,'deallocation problem',&
           __LINE__,__FILE__)
    endif
  endif
#endif
end subroutine mp_sum_in_place_int8_r1
subroutine mp_sum_in_place_real8_r1(data,n,comm)
#ifdef __PARALLEL
  use mpi
#endif
  use para_global, only : para_use_mpi_in_place
  use para_global, only : para_buff_size
  use para_global, only : para_stack_buff_size
  implicit none
  ! Arguments
  integer, intent(in) :: n,comm
  real(real_8), intent(inout) , dimension(*) :: data
  ! Variables
  character(*),parameter :: procedureN='mp_sum_in_place_real8_r1'
#ifdef __PARALLEL
  integer :: ierr,m,i,stat,buff_size
  real(real_8),dimension(:),allocatable :: buff
  real(real_8),dimension(para_stack_buff_size) :: stack_buff
  if (n<=para_stack_buff_size) then
    call mpi_allreduce(data,stack_buff,n,MPI_DOUBLE_PRECISION,&
         MPI_SUM,comm,ierr)
    call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
    data(1:n)=stack_buff(1:n)
  else
    if (para_use_mpi_in_place) then
      call mpi_allreduce(MPI_IN_PLACE,data,n,MPI_DOUBLE_PRECISION,&
           MPI_SUM,comm,ierr)
      call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
    else
      buff_size=min(n,para_buff_size)
      allocate(buff(buff_size),stat=stat)
      if (stat.ne.0) call stopgm(procedureN,'allocation problem',& 
           __LINE__,__FILE__)
      do i=1,n,buff_size
        m=min(buff_size,n-i+1)
        call mpi_allreduce(data(i),buff,m,MPI_DOUBLE_PRECISION,&
             MPI_SUM,comm,ierr)
        call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
        data(i:i+m-1)=buff(1:m)
      enddo
      deallocate(buff,stat=stat)
      if (stat.ne.0) call stopgm(procedureN,'deallocation problem',&
           __LINE__,__FILE__)
    endif
  endif
#endif
end subroutine mp_sum_in_place_real8_r1
subroutine mp_sum_in_place_real8_r2(data,n,comm)
#ifdef __PARALLEL
  use mpi
#endif
  use para_global, only : para_use_mpi_in_place
  use para_global, only : para_buff_size
  use para_global, only : para_stack_buff_size
  implicit none
  ! Arguments
  integer, intent(in) :: n,comm
  real(real_8), intent(inout) , dimension(1,*) :: data
  ! Variables
  character(*),parameter :: procedureN='mp_sum_in_place_real8_r2'
#ifdef __PARALLEL
  integer :: ierr,m,i,stat,buff_size
  real(real_8),dimension(:),allocatable :: buff
  real(real_8),dimension(para_stack_buff_size) :: stack_buff
  if (n<=para_stack_buff_size) then
    call mpi_allreduce(data,stack_buff,n,MPI_DOUBLE_PRECISION,&
         MPI_SUM,comm,ierr)
    call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
    data(1,1:n)=stack_buff(1:n)
  else
    if (para_use_mpi_in_place) then
      call mpi_allreduce(MPI_IN_PLACE,data,n,MPI_DOUBLE_PRECISION,&
           MPI_SUM,comm,ierr)
      call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
    else
      buff_size=min(n,para_buff_size)
      allocate(buff(buff_size),stat=stat)
      if (stat.ne.0) call stopgm(procedureN,'allocation problem',& 
           __LINE__,__FILE__)
      do i=1,n,buff_size
        m=min(buff_size,n-i+1)
        call mpi_allreduce(data(1,i),buff,m,MPI_DOUBLE_PRECISION,&
             MPI_SUM,comm,ierr)
        call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
        data(1,i:i+m-1)=buff(1:m)
      enddo
      deallocate(buff,stat=stat)
      if (stat.ne.0) call stopgm(procedureN,'deallocation problem',&
           __LINE__,__FILE__)
    endif
  endif
#endif
end subroutine mp_sum_in_place_real8_r2
subroutine mp_sum_in_place_real8_r3(data,n,comm)
#ifdef __PARALLEL
  use mpi
#endif
  use para_global, only : para_use_mpi_in_place
  use para_global, only : para_buff_size
  use para_global, only : para_stack_buff_size
  implicit none
  ! Arguments
  integer, intent(in) :: n,comm
  real(real_8), intent(inout) , dimension(1,1,*) :: data
  ! Variables
  character(*),parameter :: procedureN='mp_sum_in_place_real8_r3'
#ifdef __PARALLEL
  integer :: ierr,m,i,stat,buff_size
  real(real_8),dimension(:),allocatable :: buff
  real(real_8),dimension(para_stack_buff_size) :: stack_buff
  if (n<=para_stack_buff_size) then
    call mpi_allreduce(data,stack_buff,n,MPI_DOUBLE_PRECISION,&
         MPI_SUM,comm,ierr)
    call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
    data(1,1,1:n)=stack_buff(1:n)
  else
    if (para_use_mpi_in_place) then
      call mpi_allreduce(MPI_IN_PLACE,data,n,MPI_DOUBLE_PRECISION,&
           MPI_SUM,comm,ierr)
      call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
    else
      buff_size=min(n,para_buff_size)
      allocate(buff(buff_size),stat=stat)
      if (stat.ne.0) call stopgm(procedureN,'allocation problem',& 
           __LINE__,__FILE__)
      do i=1,n,buff_size
        m=min(buff_size,n-i+1)
        call mpi_allreduce(data(1,1,i),buff,m,MPI_DOUBLE_PRECISION,&
             MPI_SUM,comm,ierr)
        call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
        data(1,1,i:i+m-1)=buff(1:m)
      enddo
      deallocate(buff,stat=stat)
      if (stat.ne.0) call stopgm(procedureN,'deallocation problem',&
           __LINE__,__FILE__)
    endif
  endif
#endif
end subroutine mp_sum_in_place_real8_r3
subroutine mp_sum_in_place_real8_r4(data,n,comm)
#ifdef __PARALLEL
  use mpi
#endif
  use para_global, only : para_use_mpi_in_place
  use para_global, only : para_buff_size
  use para_global, only : para_stack_buff_size
  implicit none
  ! Arguments
  integer, intent(in) :: n,comm
  real(real_8), intent(inout) , dimension(1,1,1,*) :: data
  ! Variables
  character(*),parameter :: procedureN='mp_sum_in_place_real8_r4'
#ifdef __PARALLEL
  integer :: ierr,m,i,stat,buff_size
  real(real_8),dimension(:),allocatable :: buff
  real(real_8),dimension(para_stack_buff_size) :: stack_buff
  if (n<=para_stack_buff_size) then
    call mpi_allreduce(data,stack_buff,n,MPI_DOUBLE_PRECISION,&
         MPI_SUM,comm,ierr)
    call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
    data(1,1,1,1:n)=stack_buff(1:n)
  else
    if (para_use_mpi_in_place) then
      call mpi_allreduce(MPI_IN_PLACE,data,n,MPI_DOUBLE_PRECISION,&
           MPI_SUM,comm,ierr)
      call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
    else
      buff_size=min(n,para_buff_size)
      allocate(buff(buff_size),stat=stat)
      if (stat.ne.0) call stopgm(procedureN,'allocation problem',& 
           __LINE__,__FILE__)
      do i=1,n,buff_size
        m=min(buff_size,n-i+1)
        call mpi_allreduce(data(1,1,1,i),buff,m,MPI_DOUBLE_PRECISION,&
             MPI_SUM,comm,ierr)
        call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
        data(1,1,1,i:i+m-1)=buff(1:m)
      enddo
      deallocate(buff,stat=stat)
      if (stat.ne.0) call stopgm(procedureN,'deallocation problem',&
           __LINE__,__FILE__)
    endif
  endif
#endif
end subroutine mp_sum_in_place_real8_r4
subroutine mp_sum_in_place_real8_r5(data,n,comm)
#ifdef __PARALLEL
  use mpi
#endif
  use para_global, only : para_use_mpi_in_place
  use para_global, only : para_buff_size
  use para_global, only : para_stack_buff_size
  implicit none
  ! Arguments
  integer, intent(in) :: n,comm
  real(real_8), intent(inout) , dimension(1,1,1,1,*) :: data
  ! Variables
  character(*),parameter :: procedureN='mp_sum_in_place_real8_r5'
#ifdef __PARALLEL
  integer :: ierr,m,i,stat,buff_size
  real(real_8),dimension(:),allocatable :: buff
  real(real_8),dimension(para_stack_buff_size) :: stack_buff
  if (n<=para_stack_buff_size) then
    call mpi_allreduce(data,stack_buff,n,MPI_DOUBLE_PRECISION,&
         MPI_SUM,comm,ierr)
    call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
    data(1,1,1,1,1:n)=stack_buff(1:n)
  else
    if (para_use_mpi_in_place) then
      call mpi_allreduce(MPI_IN_PLACE,data,n,MPI_DOUBLE_PRECISION,&
           MPI_SUM,comm,ierr)
      call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
    else
      buff_size=min(n,para_buff_size)
      allocate(buff(buff_size),stat=stat)
      if (stat.ne.0) call stopgm(procedureN,'allocation problem',& 
           __LINE__,__FILE__)
      do i=1,n,buff_size
        m=min(buff_size,n-i+1)
        call mpi_allreduce(data(1,1,1,1,i),buff,m,MPI_DOUBLE_PRECISION,&
             MPI_SUM,comm,ierr)
        call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
        data(1,1,1,1,i:i+m-1)=buff(1:m)
      enddo
      deallocate(buff,stat=stat)
      if (stat.ne.0) call stopgm(procedureN,'deallocation problem',&
           __LINE__,__FILE__)
    endif
  endif
#endif
end subroutine mp_sum_in_place_real8_r5
subroutine mp_sum_in_place_complex8_r1(data,n,comm)
#ifdef __PARALLEL
  use mpi
#endif
  use para_global, only : para_use_mpi_in_place
  use para_global, only : para_buff_size
  use para_global, only : para_stack_buff_size
  implicit none
  ! Arguments
  integer, intent(in) :: n,comm
  complex(real_8), intent(inout) , dimension(*) :: data
  ! Variables
  character(*),parameter :: procedureN='mp_sum_in_place_complex8_r1'
#ifdef __PARALLEL
  integer :: ierr,m,i,stat,buff_size
  complex(real_8),dimension(:),allocatable :: buff
  complex(real_8),dimension(para_stack_buff_size) :: stack_buff
  if (n<=para_stack_buff_size) then
    call mpi_allreduce(data,stack_buff,n,MPI_DOUBLE_COMPLEX,&
         MPI_SUM,comm,ierr)
    call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
    data(1:n)=stack_buff(1:n)
  else
    if (para_use_mpi_in_place) then
      call mpi_allreduce(MPI_IN_PLACE,data,n,MPI_DOUBLE_COMPLEX,&
           MPI_SUM,comm,ierr)
      call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
    else
      buff_size=min(n,para_buff_size)
      allocate(buff(buff_size),stat=stat)
      if (stat.ne.0) call stopgm(procedureN,'allocation problem',& 
           __LINE__,__FILE__)
      do i=1,n,buff_size
        m=min(buff_size,n-i+1)
        call mpi_allreduce(data(i),buff,m,MPI_DOUBLE_COMPLEX,&
             MPI_SUM,comm,ierr)
        call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
        data(i:i+m-1)=buff(1:m)
      enddo
      deallocate(buff,stat=stat)
      if (stat.ne.0) call stopgm(procedureN,'deallocation problem',&
           __LINE__,__FILE__)
    endif
  endif
#endif
end subroutine mp_sum_in_place_complex8_r1
subroutine mp_sum_in_place_complex8_r2(data,n,comm)
#ifdef __PARALLEL
  use mpi
#endif
  use para_global, only : para_use_mpi_in_place
  use para_global, only : para_buff_size
  use para_global, only : para_stack_buff_size
  implicit none
  ! Arguments
  integer, intent(in) :: n,comm
  complex(real_8), intent(inout) , dimension(1,*) :: data
  ! Variables
  character(*),parameter :: procedureN='mp_sum_in_place_complex8_r2'
#ifdef __PARALLEL
  integer :: ierr,m,i,stat,buff_size
  complex(real_8),dimension(:),allocatable :: buff
  complex(real_8),dimension(para_stack_buff_size) :: stack_buff
  if (n<=para_stack_buff_size) then
    call mpi_allreduce(data,stack_buff,n,MPI_DOUBLE_COMPLEX,&
         MPI_SUM,comm,ierr)
    call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
    data(1,1:n)=stack_buff(1:n)
  else
    if (para_use_mpi_in_place) then
      call mpi_allreduce(MPI_IN_PLACE,data,n,MPI_DOUBLE_COMPLEX,&
           MPI_SUM,comm,ierr)
      call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
    else
      buff_size=min(n,para_buff_size)
      allocate(buff(buff_size),stat=stat)
      if (stat.ne.0) call stopgm(procedureN,'allocation problem',& 
           __LINE__,__FILE__)
      do i=1,n,buff_size
        m=min(buff_size,n-i+1)
        call mpi_allreduce(data(1,i),buff,m,MPI_DOUBLE_COMPLEX,&
             MPI_SUM,comm,ierr)
        call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
        data(1,i:i+m-1)=buff(1:m)
      enddo
      deallocate(buff,stat=stat)
      if (stat.ne.0) call stopgm(procedureN,'deallocation problem',&
           __LINE__,__FILE__)
    endif
  endif
#endif
end subroutine mp_sum_in_place_complex8_r2
subroutine mp_sum_in_place_complex8_r3(data,n,comm)
#ifdef __PARALLEL
  use mpi
#endif
  use para_global, only : para_use_mpi_in_place
  use para_global, only : para_buff_size
  use para_global, only : para_stack_buff_size
  implicit none
  ! Arguments
  integer, intent(in) :: n,comm
  complex(real_8), intent(inout) , dimension(1,1,*) :: data
  ! Variables
  character(*),parameter :: procedureN='mp_sum_in_place_complex8_r3'
#ifdef __PARALLEL
  integer :: ierr,m,i,stat,buff_size
  complex(real_8),dimension(:),allocatable :: buff
  complex(real_8),dimension(para_stack_buff_size) :: stack_buff
  if (n<=para_stack_buff_size) then
    call mpi_allreduce(data,stack_buff,n,MPI_DOUBLE_COMPLEX,&
         MPI_SUM,comm,ierr)
    call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
    data(1,1,1:n)=stack_buff(1:n)
  else
    if (para_use_mpi_in_place) then
      call mpi_allreduce(MPI_IN_PLACE,data,n,MPI_DOUBLE_COMPLEX,&
           MPI_SUM,comm,ierr)
      call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
    else
      buff_size=min(n,para_buff_size)
      allocate(buff(buff_size),stat=stat)
      if (stat.ne.0) call stopgm(procedureN,'allocation problem',& 
           __LINE__,__FILE__)
      do i=1,n,buff_size
        m=min(buff_size,n-i+1)
        call mpi_allreduce(data(1,1,i),buff,m,MPI_DOUBLE_COMPLEX,&
             MPI_SUM,comm,ierr)
        call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
        data(1,1,i:i+m-1)=buff(1:m)
      enddo
      deallocate(buff,stat=stat)
      if (stat.ne.0) call stopgm(procedureN,'deallocation problem',&
           __LINE__,__FILE__)
    endif
  endif
#endif
end subroutine mp_sum_in_place_complex8_r3
subroutine mp_sum_in_place_complex8_r4(data,n,comm)
#ifdef __PARALLEL
  use mpi
#endif
  use para_global, only : para_use_mpi_in_place
  use para_global, only : para_buff_size
  use para_global, only : para_stack_buff_size
  implicit none
  ! Arguments
  integer, intent(in) :: n,comm
  complex(real_8), intent(inout) , dimension(1,1,1,*) :: data
  ! Variables
  character(*),parameter :: procedureN='mp_sum_in_place_complex8_r4'
#ifdef __PARALLEL
  integer :: ierr,m,i,stat,buff_size
  complex(real_8),dimension(:),allocatable :: buff
  complex(real_8),dimension(para_stack_buff_size) :: stack_buff
  if (n<=para_stack_buff_size) then
    call mpi_allreduce(data,stack_buff,n,MPI_DOUBLE_COMPLEX,&
         MPI_SUM,comm,ierr)
    call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
    data(1,1,1,1:n)=stack_buff(1:n)
  else
    if (para_use_mpi_in_place) then
      call mpi_allreduce(MPI_IN_PLACE,data,n,MPI_DOUBLE_COMPLEX,&
           MPI_SUM,comm,ierr)
      call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
    else
      buff_size=min(n,para_buff_size)
      allocate(buff(buff_size),stat=stat)
      if (stat.ne.0) call stopgm(procedureN,'allocation problem',& 
           __LINE__,__FILE__)
      do i=1,n,buff_size
        m=min(buff_size,n-i+1)
        call mpi_allreduce(data(1,1,1,i),buff,m,MPI_DOUBLE_COMPLEX,&
             MPI_SUM,comm,ierr)
        call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
        data(1,1,1,i:i+m-1)=buff(1:m)
      enddo
      deallocate(buff,stat=stat)
      if (stat.ne.0) call stopgm(procedureN,'deallocation problem',&
           __LINE__,__FILE__)
    endif
  endif
#endif
end subroutine mp_sum_in_place_complex8_r4
