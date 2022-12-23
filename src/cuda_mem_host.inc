SUBROUTINE cuda_alloc_host_int4_r1 ( mem, sizes, lbs )
  integer(int_4) , dimension(:) , POINTER __CONTIGUOUS    :: mem
  INTEGER, INTENT(IN)          :: sizes
  INTEGER, INTENT(IN), OPTIONAL  :: lbs

  CHARACTER(*), PARAMETER :: procedureN = 'cuda_alloc_host_int4_r1'
  integer(int_4) , dimension(:) , POINTER __CONTIGUOUS    :: mem_tmp
  INTEGER(KIND(cudaSuccess))           :: c_status
  INTEGER(C_SIZE_T)                    :: c_n_bytes
  TYPE( C_PTR )                        :: cptr

#if defined(_HAS_CUDA)

  c_n_bytes = INT( sizes, C_SIZE_T ) * INT( sizeof_int_4, C_SIZE_T )

  c_status = cudaMallocHost( cptr, c_n_bytes )
  IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))//&
         ' trying to allocate '//TRIM(int2str( INT(c_n_bytes,int_8) ))//' bytes'&
         //". "//trim(cuda_get_error_string(c_status)),&
       __LINE__,__FILE__)

  CALL C_F_POINTER( cptr, mem_tmp, [ sizes ] )
  IF(PRESENT( lbs ) ) THEN
     mem (lbs:) => mem_tmp
  ELSE
     mem => mem_tmp
  ENDIF

#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

END SUBROUTINE cuda_alloc_host_int4_r1


SUBROUTINE cuda_alloc_host_portable_int4_r1 ( mem, sizes, lbs )
  integer(int_4) , dimension(:) , POINTER __CONTIGUOUS    :: mem
  INTEGER, INTENT(IN)          :: sizes
  INTEGER, INTENT(IN), OPTIONAL  :: lbs

  CHARACTER(*), PARAMETER :: procedureN = 'cuda_alloc_host_portable_int4_r1'
  integer(int_4) , dimension(:) , POINTER __CONTIGUOUS    :: mem_tmp
  INTEGER(KIND(cudaSuccess))           :: c_status
  INTEGER(C_SIZE_T)                    :: c_n_bytes
  TYPE( C_PTR )                        :: cptr

#if defined(_HAS_CUDA)

  c_n_bytes = INT( sizes, C_SIZE_T ) * INT( sizeof_int_4, C_SIZE_T )

  c_status = cudaHostAlloc( cptr, c_n_bytes, cudaHostAllocPortable )
  IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))//&
         ' trying to allocate '//TRIM(int2str( INT(c_n_bytes,int_8) ))//' bytes'&
         //". "//trim(cuda_get_error_string(c_status)),&
       __LINE__,__FILE__)

  CALL C_F_POINTER( cptr, mem_tmp, [ sizes ] )
  IF(PRESENT( lbs ) ) THEN
     mem (lbs:) => mem_tmp
  ELSE
     mem => mem_tmp
  ENDIF

#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

END SUBROUTINE cuda_alloc_host_portable_int4_r1


SUBROUTINE cuda_dealloc_host_int4_r1 ( mem )
  integer(int_4) , dimension(:), POINTER __CONTIGUOUS  :: mem

  CHARACTER(*), PARAMETER               :: procedureN = 'cuda_dealloc_host_int4_r1'
  INTEGER(KIND(cudaSuccess))            :: c_status

#if defined(_HAS_CUDA)

  c_status = cudaFreeHost( C_LOC( mem ) )
  IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))&
         //". "//trim(cuda_get_error_string(c_status)),&
       __LINE__,__FILE__)

  NULLIFY( mem )

#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

END SUBROUTINE cuda_dealloc_host_int4_r1


SUBROUTINE cuda_host_register_int4_r1 ( mem )
  integer(int_4) , dimension(:), INTENT(IN), TARGET :: mem

  CHARACTER(*), PARAMETER :: procedureN = 'cuda_host_register_int4_r1'

  INTEGER(C_SIZE_T)                        :: c_n_bytes
  INTEGER(KIND(cudaSuccess))               :: c_status

#if defined(_HAS_CUDA)

  c_n_bytes = INT( SIZE( mem ), C_SIZE_T ) * INT( sizeof_complex_8, C_SIZE_T )

  c_status = cudaHostRegister ( C_LOC( mem ), c_n_bytes, cudaHostRegisterPortable )
  IF( c_status /= cudaErrorHostMemoryAlreadyRegistered ) THEN
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)
  ENDIF

#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

END SUBROUTINE cuda_host_register_int4_r1

SUBROUTINE cuda_host_unregister_int4_r1 ( mem )
  integer(int_4) , dimension(:), INTENT(IN), TARGET :: mem

  CHARACTER(*), PARAMETER :: procedureN = 'cuda_host_unregister_int4_r1'

  INTEGER(KIND(cudaSuccess))               :: c_status

#if defined(_HAS_CUDA)

  c_status = cudaHostUnregister ( C_LOC( mem ) )
  IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))&
         //". "//trim(cuda_get_error_string(c_status)),&
       __LINE__,__FILE__)

#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

END SUBROUTINE cuda_host_unregister_int4_r1
SUBROUTINE cuda_alloc_host_int4_r2 ( mem, sizes, lbs )
  integer(int_4) , dimension(:,:) , POINTER __CONTIGUOUS    :: mem
  INTEGER, INTENT(IN) , DIMENSION(:)         :: sizes
  INTEGER, INTENT(IN), OPTIONAL , DIMENSION(:) :: lbs

  CHARACTER(*), PARAMETER :: procedureN = 'cuda_alloc_host_int4_r2'
  integer(int_4) , dimension(:,:) , POINTER __CONTIGUOUS    :: mem_tmp
  INTEGER(KIND(cudaSuccess))           :: c_status
  INTEGER(C_SIZE_T)                    :: c_n_bytes
  TYPE( C_PTR )                        :: cptr

#if defined(_HAS_CUDA)

  c_n_bytes = INT( PRODUCT(sizes), C_SIZE_T ) * INT( sizeof_int_4, C_SIZE_T )

  c_status = cudaMallocHost( cptr, c_n_bytes )
  IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))//&
         ' trying to allocate '//TRIM(int2str( INT(c_n_bytes,int_8) ))//' bytes'&
         //". "//trim(cuda_get_error_string(c_status)),&
       __LINE__,__FILE__)

  CALL C_F_POINTER( cptr, mem_tmp, [ sizes ] )
  IF(PRESENT( lbs ) ) THEN
     mem (lbs(1):,lbs(2):) => mem_tmp
  ELSE
     mem => mem_tmp
  ENDIF

#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

END SUBROUTINE cuda_alloc_host_int4_r2


SUBROUTINE cuda_alloc_host_portable_int4_r2 ( mem, sizes, lbs )
  integer(int_4) , dimension(:,:) , POINTER __CONTIGUOUS    :: mem
  INTEGER, INTENT(IN) , DIMENSION(:)         :: sizes
  INTEGER, INTENT(IN), OPTIONAL , DIMENSION(:) :: lbs

  CHARACTER(*), PARAMETER :: procedureN = 'cuda_alloc_host_portable_int4_r2'
  integer(int_4) , dimension(:,:) , POINTER __CONTIGUOUS    :: mem_tmp
  INTEGER(KIND(cudaSuccess))           :: c_status
  INTEGER(C_SIZE_T)                    :: c_n_bytes
  TYPE( C_PTR )                        :: cptr

#if defined(_HAS_CUDA)

  c_n_bytes = INT( PRODUCT(sizes), C_SIZE_T ) * INT( sizeof_int_4, C_SIZE_T )

  c_status = cudaHostAlloc( cptr, c_n_bytes, cudaHostAllocPortable )
  IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))//&
         ' trying to allocate '//TRIM(int2str( INT(c_n_bytes,int_8) ))//' bytes'&
         //". "//trim(cuda_get_error_string(c_status)),&
       __LINE__,__FILE__)

  CALL C_F_POINTER( cptr, mem_tmp, [ sizes ] )
  IF(PRESENT( lbs ) ) THEN
     mem (lbs(1):,lbs(2):) => mem_tmp
  ELSE
     mem => mem_tmp
  ENDIF

#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

END SUBROUTINE cuda_alloc_host_portable_int4_r2


SUBROUTINE cuda_dealloc_host_int4_r2 ( mem )
  integer(int_4) , dimension(:,:), POINTER __CONTIGUOUS  :: mem

  CHARACTER(*), PARAMETER               :: procedureN = 'cuda_dealloc_host_int4_r2'
  INTEGER(KIND(cudaSuccess))            :: c_status

#if defined(_HAS_CUDA)

  c_status = cudaFreeHost( C_LOC( mem ) )
  IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))&
         //". "//trim(cuda_get_error_string(c_status)),&
       __LINE__,__FILE__)

  NULLIFY( mem )

#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

END SUBROUTINE cuda_dealloc_host_int4_r2


SUBROUTINE cuda_host_register_int4_r2 ( mem )
  integer(int_4) , dimension(:,:), INTENT(IN), TARGET :: mem

  CHARACTER(*), PARAMETER :: procedureN = 'cuda_host_register_int4_r2'

  INTEGER(C_SIZE_T)                        :: c_n_bytes
  INTEGER(KIND(cudaSuccess))               :: c_status

#if defined(_HAS_CUDA)

  c_n_bytes = INT( SIZE( mem ), C_SIZE_T ) * INT( sizeof_complex_8, C_SIZE_T )

  c_status = cudaHostRegister ( C_LOC( mem ), c_n_bytes, cudaHostRegisterPortable )
  IF( c_status /= cudaErrorHostMemoryAlreadyRegistered ) THEN
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)
  ENDIF

#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

END SUBROUTINE cuda_host_register_int4_r2

SUBROUTINE cuda_host_unregister_int4_r2 ( mem )
  integer(int_4) , dimension(:,:), INTENT(IN), TARGET :: mem

  CHARACTER(*), PARAMETER :: procedureN = 'cuda_host_unregister_int4_r2'

  INTEGER(KIND(cudaSuccess))               :: c_status

#if defined(_HAS_CUDA)

  c_status = cudaHostUnregister ( C_LOC( mem ) )
  IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))&
         //". "//trim(cuda_get_error_string(c_status)),&
       __LINE__,__FILE__)

#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

END SUBROUTINE cuda_host_unregister_int4_r2
SUBROUTINE cuda_alloc_host_real8_r1 ( mem, sizes, lbs )
  real(real_8) , dimension(:) , POINTER __CONTIGUOUS    :: mem
  INTEGER, INTENT(IN)          :: sizes
  INTEGER, INTENT(IN), OPTIONAL  :: lbs

  CHARACTER(*), PARAMETER :: procedureN = 'cuda_alloc_host_real8_r1'
  real(real_8) , dimension(:) , POINTER __CONTIGUOUS    :: mem_tmp
  INTEGER(KIND(cudaSuccess))           :: c_status
  INTEGER(C_SIZE_T)                    :: c_n_bytes
  TYPE( C_PTR )                        :: cptr

#if defined(_HAS_CUDA)

  c_n_bytes = INT( sizes, C_SIZE_T ) * INT( sizeof_real_8, C_SIZE_T )

  c_status = cudaMallocHost( cptr, c_n_bytes )
  IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))//&
         ' trying to allocate '//TRIM(int2str( INT(c_n_bytes,int_8) ))//' bytes'&
         //". "//trim(cuda_get_error_string(c_status)),&
       __LINE__,__FILE__)

  CALL C_F_POINTER( cptr, mem_tmp, [ sizes ] )
  IF(PRESENT( lbs ) ) THEN
     mem (lbs:) => mem_tmp
  ELSE
     mem => mem_tmp
  ENDIF

#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

END SUBROUTINE cuda_alloc_host_real8_r1


SUBROUTINE cuda_alloc_host_portable_real8_r1 ( mem, sizes, lbs )
  real(real_8) , dimension(:) , POINTER __CONTIGUOUS    :: mem
  INTEGER, INTENT(IN)          :: sizes
  INTEGER, INTENT(IN), OPTIONAL  :: lbs

  CHARACTER(*), PARAMETER :: procedureN = 'cuda_alloc_host_portable_real8_r1'
  real(real_8) , dimension(:) , POINTER __CONTIGUOUS    :: mem_tmp
  INTEGER(KIND(cudaSuccess))           :: c_status
  INTEGER(C_SIZE_T)                    :: c_n_bytes
  TYPE( C_PTR )                        :: cptr

#if defined(_HAS_CUDA)

  c_n_bytes = INT( sizes, C_SIZE_T ) * INT( sizeof_real_8, C_SIZE_T )

  c_status = cudaHostAlloc( cptr, c_n_bytes, cudaHostAllocPortable )
  IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))//&
         ' trying to allocate '//TRIM(int2str( INT(c_n_bytes,int_8) ))//' bytes'&
         //". "//trim(cuda_get_error_string(c_status)),&
       __LINE__,__FILE__)

  CALL C_F_POINTER( cptr, mem_tmp, [ sizes ] )
  IF(PRESENT( lbs ) ) THEN
     mem (lbs:) => mem_tmp
  ELSE
     mem => mem_tmp
  ENDIF

#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

END SUBROUTINE cuda_alloc_host_portable_real8_r1


SUBROUTINE cuda_dealloc_host_real8_r1 ( mem )
  real(real_8) , dimension(:), POINTER __CONTIGUOUS  :: mem

  CHARACTER(*), PARAMETER               :: procedureN = 'cuda_dealloc_host_real8_r1'
  INTEGER(KIND(cudaSuccess))            :: c_status

#if defined(_HAS_CUDA)

  c_status = cudaFreeHost( C_LOC( mem ) )
  IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))&
         //". "//trim(cuda_get_error_string(c_status)),&
       __LINE__,__FILE__)

  NULLIFY( mem )

#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

END SUBROUTINE cuda_dealloc_host_real8_r1


SUBROUTINE cuda_host_register_real8_r1 ( mem )
  real(real_8) , dimension(:), INTENT(IN), TARGET :: mem

  CHARACTER(*), PARAMETER :: procedureN = 'cuda_host_register_real8_r1'

  INTEGER(C_SIZE_T)                        :: c_n_bytes
  INTEGER(KIND(cudaSuccess))               :: c_status

#if defined(_HAS_CUDA)

  c_n_bytes = INT( SIZE( mem ), C_SIZE_T ) * INT( sizeof_complex_8, C_SIZE_T )

  c_status = cudaHostRegister ( C_LOC( mem ), c_n_bytes, cudaHostRegisterPortable )
  IF( c_status /= cudaErrorHostMemoryAlreadyRegistered ) THEN
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)
  ENDIF

#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

END SUBROUTINE cuda_host_register_real8_r1

SUBROUTINE cuda_host_unregister_real8_r1 ( mem )
  real(real_8) , dimension(:), INTENT(IN), TARGET :: mem

  CHARACTER(*), PARAMETER :: procedureN = 'cuda_host_unregister_real8_r1'

  INTEGER(KIND(cudaSuccess))               :: c_status

#if defined(_HAS_CUDA)

  c_status = cudaHostUnregister ( C_LOC( mem ) )
  IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))&
         //". "//trim(cuda_get_error_string(c_status)),&
       __LINE__,__FILE__)

#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

END SUBROUTINE cuda_host_unregister_real8_r1
SUBROUTINE cuda_alloc_host_real8_r2 ( mem, sizes, lbs )
  real(real_8) , dimension(:,:) , POINTER __CONTIGUOUS    :: mem
  INTEGER, INTENT(IN) , DIMENSION(:)         :: sizes
  INTEGER, INTENT(IN), OPTIONAL , DIMENSION(:) :: lbs

  CHARACTER(*), PARAMETER :: procedureN = 'cuda_alloc_host_real8_r2'
  real(real_8) , dimension(:,:) , POINTER __CONTIGUOUS    :: mem_tmp
  INTEGER(KIND(cudaSuccess))           :: c_status
  INTEGER(C_SIZE_T)                    :: c_n_bytes
  TYPE( C_PTR )                        :: cptr

#if defined(_HAS_CUDA)

  c_n_bytes = INT( PRODUCT(sizes), C_SIZE_T ) * INT( sizeof_real_8, C_SIZE_T )

  c_status = cudaMallocHost( cptr, c_n_bytes )
  IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))//&
         ' trying to allocate '//TRIM(int2str( INT(c_n_bytes,int_8) ))//' bytes'&
         //". "//trim(cuda_get_error_string(c_status)),&
       __LINE__,__FILE__)

  CALL C_F_POINTER( cptr, mem_tmp, [ sizes ] )
  IF(PRESENT( lbs ) ) THEN
     mem (lbs(1):,lbs(2):) => mem_tmp
  ELSE
     mem => mem_tmp
  ENDIF

#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

END SUBROUTINE cuda_alloc_host_real8_r2


SUBROUTINE cuda_alloc_host_portable_real8_r2 ( mem, sizes, lbs )
  real(real_8) , dimension(:,:) , POINTER __CONTIGUOUS    :: mem
  INTEGER, INTENT(IN) , DIMENSION(:)         :: sizes
  INTEGER, INTENT(IN), OPTIONAL , DIMENSION(:) :: lbs

  CHARACTER(*), PARAMETER :: procedureN = 'cuda_alloc_host_portable_real8_r2'
  real(real_8) , dimension(:,:) , POINTER __CONTIGUOUS    :: mem_tmp
  INTEGER(KIND(cudaSuccess))           :: c_status
  INTEGER(C_SIZE_T)                    :: c_n_bytes
  TYPE( C_PTR )                        :: cptr

#if defined(_HAS_CUDA)

  c_n_bytes = INT( PRODUCT(sizes), C_SIZE_T ) * INT( sizeof_real_8, C_SIZE_T )

  c_status = cudaHostAlloc( cptr, c_n_bytes, cudaHostAllocPortable )
  IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))//&
         ' trying to allocate '//TRIM(int2str( INT(c_n_bytes,int_8) ))//' bytes'&
         //". "//trim(cuda_get_error_string(c_status)),&
       __LINE__,__FILE__)

  CALL C_F_POINTER( cptr, mem_tmp, [ sizes ] )
  IF(PRESENT( lbs ) ) THEN
     mem (lbs(1):,lbs(2):) => mem_tmp
  ELSE
     mem => mem_tmp
  ENDIF

#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

END SUBROUTINE cuda_alloc_host_portable_real8_r2


SUBROUTINE cuda_dealloc_host_real8_r2 ( mem )
  real(real_8) , dimension(:,:), POINTER __CONTIGUOUS  :: mem

  CHARACTER(*), PARAMETER               :: procedureN = 'cuda_dealloc_host_real8_r2'
  INTEGER(KIND(cudaSuccess))            :: c_status

#if defined(_HAS_CUDA)

  c_status = cudaFreeHost( C_LOC( mem ) )
  IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))&
         //". "//trim(cuda_get_error_string(c_status)),&
       __LINE__,__FILE__)

  NULLIFY( mem )

#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

END SUBROUTINE cuda_dealloc_host_real8_r2


SUBROUTINE cuda_host_register_real8_r2 ( mem )
  real(real_8) , dimension(:,:), INTENT(IN), TARGET :: mem

  CHARACTER(*), PARAMETER :: procedureN = 'cuda_host_register_real8_r2'

  INTEGER(C_SIZE_T)                        :: c_n_bytes
  INTEGER(KIND(cudaSuccess))               :: c_status

#if defined(_HAS_CUDA)

  c_n_bytes = INT( SIZE( mem ), C_SIZE_T ) * INT( sizeof_complex_8, C_SIZE_T )

  c_status = cudaHostRegister ( C_LOC( mem ), c_n_bytes, cudaHostRegisterPortable )
  IF( c_status /= cudaErrorHostMemoryAlreadyRegistered ) THEN
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)
  ENDIF

#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

END SUBROUTINE cuda_host_register_real8_r2

SUBROUTINE cuda_host_unregister_real8_r2 ( mem )
  real(real_8) , dimension(:,:), INTENT(IN), TARGET :: mem

  CHARACTER(*), PARAMETER :: procedureN = 'cuda_host_unregister_real8_r2'

  INTEGER(KIND(cudaSuccess))               :: c_status

#if defined(_HAS_CUDA)

  c_status = cudaHostUnregister ( C_LOC( mem ) )
  IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))&
         //". "//trim(cuda_get_error_string(c_status)),&
       __LINE__,__FILE__)

#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

END SUBROUTINE cuda_host_unregister_real8_r2
SUBROUTINE cuda_alloc_host_real8_r3 ( mem, sizes, lbs )
  real(real_8) , dimension(:,:,:) , POINTER __CONTIGUOUS    :: mem
  INTEGER, INTENT(IN) , DIMENSION(:)         :: sizes
  INTEGER, INTENT(IN), OPTIONAL , DIMENSION(:) :: lbs

  CHARACTER(*), PARAMETER :: procedureN = 'cuda_alloc_host_real8_r3'
  real(real_8) , dimension(:,:,:) , POINTER __CONTIGUOUS    :: mem_tmp
  INTEGER(KIND(cudaSuccess))           :: c_status
  INTEGER(C_SIZE_T)                    :: c_n_bytes
  TYPE( C_PTR )                        :: cptr

#if defined(_HAS_CUDA)

  c_n_bytes = INT( PRODUCT(sizes), C_SIZE_T ) * INT( sizeof_real_8, C_SIZE_T )

  c_status = cudaMallocHost( cptr, c_n_bytes )
  IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))//&
         ' trying to allocate '//TRIM(int2str( INT(c_n_bytes,int_8) ))//' bytes'&
         //". "//trim(cuda_get_error_string(c_status)),&
       __LINE__,__FILE__)

  CALL C_F_POINTER( cptr, mem_tmp, [ sizes ] )
  IF(PRESENT( lbs ) ) THEN
     mem (lbs(1):,lbs(2):,lbs(3):) => mem_tmp
  ELSE
     mem => mem_tmp
  ENDIF

#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

END SUBROUTINE cuda_alloc_host_real8_r3


SUBROUTINE cuda_alloc_host_portable_real8_r3 ( mem, sizes, lbs )
  real(real_8) , dimension(:,:,:) , POINTER __CONTIGUOUS    :: mem
  INTEGER, INTENT(IN) , DIMENSION(:)         :: sizes
  INTEGER, INTENT(IN), OPTIONAL , DIMENSION(:) :: lbs

  CHARACTER(*), PARAMETER :: procedureN = 'cuda_alloc_host_portable_real8_r3'
  real(real_8) , dimension(:,:,:) , POINTER __CONTIGUOUS    :: mem_tmp
  INTEGER(KIND(cudaSuccess))           :: c_status
  INTEGER(C_SIZE_T)                    :: c_n_bytes
  TYPE( C_PTR )                        :: cptr

#if defined(_HAS_CUDA)

  c_n_bytes = INT( PRODUCT(sizes), C_SIZE_T ) * INT( sizeof_real_8, C_SIZE_T )

  c_status = cudaHostAlloc( cptr, c_n_bytes, cudaHostAllocPortable )
  IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))//&
         ' trying to allocate '//TRIM(int2str( INT(c_n_bytes,int_8) ))//' bytes'&
         //". "//trim(cuda_get_error_string(c_status)),&
       __LINE__,__FILE__)

  CALL C_F_POINTER( cptr, mem_tmp, [ sizes ] )
  IF(PRESENT( lbs ) ) THEN
     mem (lbs(1):,lbs(2):,lbs(3):) => mem_tmp
  ELSE
     mem => mem_tmp
  ENDIF

#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

END SUBROUTINE cuda_alloc_host_portable_real8_r3


SUBROUTINE cuda_dealloc_host_real8_r3 ( mem )
  real(real_8) , dimension(:,:,:), POINTER __CONTIGUOUS  :: mem

  CHARACTER(*), PARAMETER               :: procedureN = 'cuda_dealloc_host_real8_r3'
  INTEGER(KIND(cudaSuccess))            :: c_status

#if defined(_HAS_CUDA)

  c_status = cudaFreeHost( C_LOC( mem ) )
  IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))&
         //". "//trim(cuda_get_error_string(c_status)),&
       __LINE__,__FILE__)

  NULLIFY( mem )

#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

END SUBROUTINE cuda_dealloc_host_real8_r3


SUBROUTINE cuda_host_register_real8_r3 ( mem )
  real(real_8) , dimension(:,:,:), INTENT(IN), TARGET :: mem

  CHARACTER(*), PARAMETER :: procedureN = 'cuda_host_register_real8_r3'

  INTEGER(C_SIZE_T)                        :: c_n_bytes
  INTEGER(KIND(cudaSuccess))               :: c_status

#if defined(_HAS_CUDA)

  c_n_bytes = INT( SIZE( mem ), C_SIZE_T ) * INT( sizeof_complex_8, C_SIZE_T )

  c_status = cudaHostRegister ( C_LOC( mem ), c_n_bytes, cudaHostRegisterPortable )
  IF( c_status /= cudaErrorHostMemoryAlreadyRegistered ) THEN
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)
  ENDIF

#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

END SUBROUTINE cuda_host_register_real8_r3

SUBROUTINE cuda_host_unregister_real8_r3 ( mem )
  real(real_8) , dimension(:,:,:), INTENT(IN), TARGET :: mem

  CHARACTER(*), PARAMETER :: procedureN = 'cuda_host_unregister_real8_r3'

  INTEGER(KIND(cudaSuccess))               :: c_status

#if defined(_HAS_CUDA)

  c_status = cudaHostUnregister ( C_LOC( mem ) )
  IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))&
         //". "//trim(cuda_get_error_string(c_status)),&
       __LINE__,__FILE__)

#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

END SUBROUTINE cuda_host_unregister_real8_r3
SUBROUTINE cuda_alloc_host_complex8_r1 ( mem, sizes, lbs )
  complex(real_8) , dimension(:) , POINTER __CONTIGUOUS    :: mem
  INTEGER, INTENT(IN)          :: sizes
  INTEGER, INTENT(IN), OPTIONAL  :: lbs

  CHARACTER(*), PARAMETER :: procedureN = 'cuda_alloc_host_complex8_r1'
  complex(real_8) , dimension(:) , POINTER __CONTIGUOUS    :: mem_tmp
  INTEGER(KIND(cudaSuccess))           :: c_status
  INTEGER(C_SIZE_T)                    :: c_n_bytes
  TYPE( C_PTR )                        :: cptr

#if defined(_HAS_CUDA)

  c_n_bytes = INT( sizes, C_SIZE_T ) * INT( sizeof_complex_8, C_SIZE_T )

  c_status = cudaMallocHost( cptr, c_n_bytes )
  IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))//&
         ' trying to allocate '//TRIM(int2str( INT(c_n_bytes,int_8) ))//' bytes'&
         //". "//trim(cuda_get_error_string(c_status)),&
       __LINE__,__FILE__)

  CALL C_F_POINTER( cptr, mem_tmp, [ sizes ] )
  IF(PRESENT( lbs ) ) THEN
     mem (lbs:) => mem_tmp
  ELSE
     mem => mem_tmp
  ENDIF

#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

END SUBROUTINE cuda_alloc_host_complex8_r1


SUBROUTINE cuda_alloc_host_portable_complex8_r1 ( mem, sizes, lbs )
  complex(real_8) , dimension(:) , POINTER __CONTIGUOUS    :: mem
  INTEGER, INTENT(IN)          :: sizes
  INTEGER, INTENT(IN), OPTIONAL  :: lbs

  CHARACTER(*), PARAMETER :: procedureN = 'cuda_alloc_host_portable_complex8_r1'
  complex(real_8) , dimension(:) , POINTER __CONTIGUOUS    :: mem_tmp
  INTEGER(KIND(cudaSuccess))           :: c_status
  INTEGER(C_SIZE_T)                    :: c_n_bytes
  TYPE( C_PTR )                        :: cptr

#if defined(_HAS_CUDA)

  c_n_bytes = INT( sizes, C_SIZE_T ) * INT( sizeof_complex_8, C_SIZE_T )

  c_status = cudaHostAlloc( cptr, c_n_bytes, cudaHostAllocPortable )
  IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))//&
         ' trying to allocate '//TRIM(int2str( INT(c_n_bytes,int_8) ))//' bytes'&
         //". "//trim(cuda_get_error_string(c_status)),&
       __LINE__,__FILE__)

  CALL C_F_POINTER( cptr, mem_tmp, [ sizes ] )
  IF(PRESENT( lbs ) ) THEN
     mem (lbs:) => mem_tmp
  ELSE
     mem => mem_tmp
  ENDIF

#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

END SUBROUTINE cuda_alloc_host_portable_complex8_r1


SUBROUTINE cuda_dealloc_host_complex8_r1 ( mem )
  complex(real_8) , dimension(:), POINTER __CONTIGUOUS  :: mem

  CHARACTER(*), PARAMETER               :: procedureN = 'cuda_dealloc_host_complex8_r1'
  INTEGER(KIND(cudaSuccess))            :: c_status

#if defined(_HAS_CUDA)

  c_status = cudaFreeHost( C_LOC( mem ) )
  IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))&
         //". "//trim(cuda_get_error_string(c_status)),&
       __LINE__,__FILE__)

  NULLIFY( mem )

#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

END SUBROUTINE cuda_dealloc_host_complex8_r1


SUBROUTINE cuda_host_register_complex8_r1 ( mem )
  complex(real_8) , dimension(:), INTENT(IN), TARGET :: mem

  CHARACTER(*), PARAMETER :: procedureN = 'cuda_host_register_complex8_r1'

  INTEGER(C_SIZE_T)                        :: c_n_bytes
  INTEGER(KIND(cudaSuccess))               :: c_status

#if defined(_HAS_CUDA)

  c_n_bytes = INT( SIZE( mem ), C_SIZE_T ) * INT( sizeof_complex_8, C_SIZE_T )

  c_status = cudaHostRegister ( C_LOC( mem ), c_n_bytes, cudaHostRegisterPortable )
  IF( c_status /= cudaErrorHostMemoryAlreadyRegistered ) THEN
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)
  ENDIF

#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

END SUBROUTINE cuda_host_register_complex8_r1

SUBROUTINE cuda_host_unregister_complex8_r1 ( mem )
  complex(real_8) , dimension(:), INTENT(IN), TARGET :: mem

  CHARACTER(*), PARAMETER :: procedureN = 'cuda_host_unregister_complex8_r1'

  INTEGER(KIND(cudaSuccess))               :: c_status

#if defined(_HAS_CUDA)

  c_status = cudaHostUnregister ( C_LOC( mem ) )
  IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))&
         //". "//trim(cuda_get_error_string(c_status)),&
       __LINE__,__FILE__)

#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

END SUBROUTINE cuda_host_unregister_complex8_r1
SUBROUTINE cuda_alloc_host_complex8_r2 ( mem, sizes, lbs )
  complex(real_8) , dimension(:,:) , POINTER __CONTIGUOUS    :: mem
  INTEGER, INTENT(IN) , DIMENSION(:)         :: sizes
  INTEGER, INTENT(IN), OPTIONAL , DIMENSION(:) :: lbs

  CHARACTER(*), PARAMETER :: procedureN = 'cuda_alloc_host_complex8_r2'
  complex(real_8) , dimension(:,:) , POINTER __CONTIGUOUS    :: mem_tmp
  INTEGER(KIND(cudaSuccess))           :: c_status
  INTEGER(C_SIZE_T)                    :: c_n_bytes
  TYPE( C_PTR )                        :: cptr

#if defined(_HAS_CUDA)

  c_n_bytes = INT( PRODUCT(sizes), C_SIZE_T ) * INT( sizeof_complex_8, C_SIZE_T )

  c_status = cudaMallocHost( cptr, c_n_bytes )
  IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))//&
         ' trying to allocate '//TRIM(int2str( INT(c_n_bytes,int_8) ))//' bytes'&
         //". "//trim(cuda_get_error_string(c_status)),&
       __LINE__,__FILE__)

  CALL C_F_POINTER( cptr, mem_tmp, [ sizes ] )
  IF(PRESENT( lbs ) ) THEN
     mem (lbs(1):,lbs(2):) => mem_tmp
  ELSE
     mem => mem_tmp
  ENDIF

#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

END SUBROUTINE cuda_alloc_host_complex8_r2


SUBROUTINE cuda_alloc_host_portable_complex8_r2 ( mem, sizes, lbs )
  complex(real_8) , dimension(:,:) , POINTER __CONTIGUOUS    :: mem
  INTEGER, INTENT(IN) , DIMENSION(:)         :: sizes
  INTEGER, INTENT(IN), OPTIONAL , DIMENSION(:) :: lbs

  CHARACTER(*), PARAMETER :: procedureN = 'cuda_alloc_host_portable_complex8_r2'
  complex(real_8) , dimension(:,:) , POINTER __CONTIGUOUS    :: mem_tmp
  INTEGER(KIND(cudaSuccess))           :: c_status
  INTEGER(C_SIZE_T)                    :: c_n_bytes
  TYPE( C_PTR )                        :: cptr

#if defined(_HAS_CUDA)

  c_n_bytes = INT( PRODUCT(sizes), C_SIZE_T ) * INT( sizeof_complex_8, C_SIZE_T )

  c_status = cudaHostAlloc( cptr, c_n_bytes, cudaHostAllocPortable )
  IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))//&
         ' trying to allocate '//TRIM(int2str( INT(c_n_bytes,int_8) ))//' bytes'&
         //". "//trim(cuda_get_error_string(c_status)),&
       __LINE__,__FILE__)

  CALL C_F_POINTER( cptr, mem_tmp, [ sizes ] )
  IF(PRESENT( lbs ) ) THEN
     mem (lbs(1):,lbs(2):) => mem_tmp
  ELSE
     mem => mem_tmp
  ENDIF

#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

END SUBROUTINE cuda_alloc_host_portable_complex8_r2


SUBROUTINE cuda_dealloc_host_complex8_r2 ( mem )
  complex(real_8) , dimension(:,:), POINTER __CONTIGUOUS  :: mem

  CHARACTER(*), PARAMETER               :: procedureN = 'cuda_dealloc_host_complex8_r2'
  INTEGER(KIND(cudaSuccess))            :: c_status

#if defined(_HAS_CUDA)

  c_status = cudaFreeHost( C_LOC( mem ) )
  IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))&
         //". "//trim(cuda_get_error_string(c_status)),&
       __LINE__,__FILE__)

  NULLIFY( mem )

#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

END SUBROUTINE cuda_dealloc_host_complex8_r2


SUBROUTINE cuda_host_register_complex8_r2 ( mem )
  complex(real_8) , dimension(:,:), INTENT(IN), TARGET :: mem

  CHARACTER(*), PARAMETER :: procedureN = 'cuda_host_register_complex8_r2'

  INTEGER(C_SIZE_T)                        :: c_n_bytes
  INTEGER(KIND(cudaSuccess))               :: c_status

#if defined(_HAS_CUDA)

  c_n_bytes = INT( SIZE( mem ), C_SIZE_T ) * INT( sizeof_complex_8, C_SIZE_T )

  c_status = cudaHostRegister ( C_LOC( mem ), c_n_bytes, cudaHostRegisterPortable )
  IF( c_status /= cudaErrorHostMemoryAlreadyRegistered ) THEN
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)
  ENDIF

#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

END SUBROUTINE cuda_host_register_complex8_r2

SUBROUTINE cuda_host_unregister_complex8_r2 ( mem )
  complex(real_8) , dimension(:,:), INTENT(IN), TARGET :: mem

  CHARACTER(*), PARAMETER :: procedureN = 'cuda_host_unregister_complex8_r2'

  INTEGER(KIND(cudaSuccess))               :: c_status

#if defined(_HAS_CUDA)

  c_status = cudaHostUnregister ( C_LOC( mem ) )
  IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//TRIM(int2str( INT( c_status ) ))&
         //". "//trim(cuda_get_error_string(c_status)),&
       __LINE__,__FILE__)

#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

END SUBROUTINE cuda_host_unregister_complex8_r2
