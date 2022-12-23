subroutine cuda_memcpy_host_to_device_int4_r1 ( a_host, a_devc )
  integer(int_4) , dimension(:), intent(in), target :: a_host
  TYPE( cuda_memory_t ) :: a_devc

  INTEGER :: n
  INTEGER(KIND(cudaSuccess)) :: c_status
  INTEGER( C_SIZE_T ) :: c_n_bytes
  CHARACTER(*), PARAMETER :: procedureN = 'cuda_memcpy_host_to_device_int4_r1'

#if defined(_HAS_CUDA)

  IF( .NOT.a_devc%init ) CALL stopgm(procedureN,'memory not allocated',&
       __LINE__,__FILE__)

  CALL cuda_check_device ( a_devc%device, procedureN )

  n = SIZE( a_host )

  IF( n > 0 ) THEN

     c_n_bytes = INT( n, C_SIZE_T )* INT( sizeof_int_4, C_SIZE_T )

     c_status = cudaMemcpy( a_devc%ptr, C_LOC( a_host ), c_n_bytes, cudaMemcpyHostToDevice )
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//int2str( INT( c_status ) )&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)

  ENDIF
#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

end subroutine cuda_memcpy_host_to_device_int4_r1


subroutine cuda_memcpy_device_to_host_int4_r1 ( b_devc, b_host )
  integer(int_4) , dimension(:), intent(inout), target :: b_host
  TYPE( cuda_memory_t ) :: b_devc

  INTEGER :: n
  INTEGER(KIND(cudaSuccess)) :: c_status
  INTEGER( C_SIZE_T ) :: c_n_bytes
  CHARACTER(*), PARAMETER :: procedureN = 'cuda_memcpy_device_to_host_int4_r1'

#if defined(_HAS_CUDA)

  IF( .NOT.b_devc%init ) CALL stopgm(procedureN,'memory not allocated',&
       __LINE__,__FILE__)

  CALL cuda_check_device ( b_devc%device, procedureN )

  n = SIZE( b_host )

  IF( n > 0 ) THEN

     c_n_bytes = INT( n, C_SIZE_T )* INT( sizeof_int_4, C_SIZE_T )

     c_status = cudaMemcpy( C_LOC( b_host ), b_devc%ptr, c_n_bytes, cudaMemcpyDeviceToHost )
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//int2str( INT( c_status ) )&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)

  ENDIF
#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

end subroutine cuda_memcpy_device_to_host_int4_r1


subroutine cuda_memcpy_async_host_to_device_int4_r1 ( a_host, a_devc, stream )
  integer(int_4) , dimension(:), intent(in), target :: a_host
  TYPE( cuda_memory_t ) :: a_devc
  TYPE( cuda_stream_t ) :: stream

  INTEGER :: n
  INTEGER(KIND(cudaSuccess)) :: c_status
  INTEGER( C_SIZE_T ) :: c_n_bytes
  CHARACTER(*), PARAMETER :: procedureN = 'cuda_memcpy_async_host_to_device_int4_r1'

#if defined(_HAS_CUDA)

  IF( .NOT.a_devc%init ) CALL stopgm(procedureN,'memory not allocated',&
       __LINE__,__FILE__)
  IF( .NOT. stream%init ) CALL stopgm(procedureN,'stream not initialized',&
       __LINE__,__FILE__)
  IF( a_devc%device /= stream%device ) CALL stopgm(procedureN,'devices not consistent',&
       __LINE__,__FILE__)
  CALL cuda_check_device ( a_devc%device, procedureN )

  n = SIZE( a_host )

  IF( n > 0 ) THEN

     c_n_bytes = INT( n, C_SIZE_T )* INT( sizeof_int_4, C_SIZE_T )

     c_status = cudaMemcpyAsync( a_devc%ptr, C_LOC( a_host ), c_n_bytes, cudaMemcpyHostToDevice, stream%s )
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//int2str( INT( c_status ) )&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)

  ENDIF
#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

end subroutine cuda_memcpy_async_host_to_device_int4_r1


subroutine cuda_memcpy_async_device_to_host_int4_r1 ( b_devc, b_host, stream )
  TYPE( cuda_memory_t ) :: b_devc
  integer(int_4) , dimension(:), intent(inout), target :: b_host
  TYPE( cuda_stream_t ) :: stream

  INTEGER :: n
  INTEGER(KIND(cudaSuccess)) :: c_status
  INTEGER( C_SIZE_T ) :: c_n_bytes
  CHARACTER(*), PARAMETER :: procedureN = 'cuda_memcpy_async_device_to_host_int4_r1'

#if defined(_HAS_CUDA)

  IF( .NOT.b_devc%init ) CALL stopgm(procedureN,'memory not allocated',&
       __LINE__,__FILE__)
  IF( .NOT. stream%init ) CALL stopgm(procedureN,'stream not initialized',&
       __LINE__,__FILE__)
  IF( b_devc%device /= stream%device ) CALL stopgm(procedureN,'devices not consistent',&
       __LINE__,__FILE__)
  CALL cuda_check_device ( b_devc%device, procedureN )

  n = SIZE( b_host )

  IF( n > 0 ) THEN

     c_n_bytes = INT( n, C_SIZE_T )* INT( sizeof_int_4, C_SIZE_T )

     c_status = cudaMemcpyAsync( C_LOC( b_host ), b_devc%ptr, c_n_bytes, cudaMemcpyDeviceToHost, stream%s )
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//int2str( INT( c_status ) )&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)

  ENDIF
#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

end subroutine cuda_memcpy_async_device_to_host_int4_r1


subroutine cuda_memcpy_host_to_device_int4_r2 ( a_host, a_devc )
  integer(int_4) , dimension(:,:), intent(in), target :: a_host
  TYPE( cuda_memory_t ) :: a_devc

  INTEGER :: n
  INTEGER(KIND(cudaSuccess)) :: c_status
  INTEGER( C_SIZE_T ) :: c_n_bytes
  CHARACTER(*), PARAMETER :: procedureN = 'cuda_memcpy_host_to_device_int4_r2'

#if defined(_HAS_CUDA)

  IF( .NOT.a_devc%init ) CALL stopgm(procedureN,'memory not allocated',&
       __LINE__,__FILE__)

  CALL cuda_check_device ( a_devc%device, procedureN )

  n = SIZE( a_host )

  IF( n > 0 ) THEN

     c_n_bytes = INT( n, C_SIZE_T )* INT( sizeof_int_4, C_SIZE_T )

     c_status = cudaMemcpy( a_devc%ptr, C_LOC( a_host ), c_n_bytes, cudaMemcpyHostToDevice )
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//int2str( INT( c_status ) )&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)

  ENDIF
#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

end subroutine cuda_memcpy_host_to_device_int4_r2


subroutine cuda_memcpy_device_to_host_int4_r2 ( b_devc, b_host )
  integer(int_4) , dimension(:,:), intent(inout), target :: b_host
  TYPE( cuda_memory_t ) :: b_devc

  INTEGER :: n
  INTEGER(KIND(cudaSuccess)) :: c_status
  INTEGER( C_SIZE_T ) :: c_n_bytes
  CHARACTER(*), PARAMETER :: procedureN = 'cuda_memcpy_device_to_host_int4_r2'

#if defined(_HAS_CUDA)

  IF( .NOT.b_devc%init ) CALL stopgm(procedureN,'memory not allocated',&
       __LINE__,__FILE__)

  CALL cuda_check_device ( b_devc%device, procedureN )

  n = SIZE( b_host )

  IF( n > 0 ) THEN

     c_n_bytes = INT( n, C_SIZE_T )* INT( sizeof_int_4, C_SIZE_T )

     c_status = cudaMemcpy( C_LOC( b_host ), b_devc%ptr, c_n_bytes, cudaMemcpyDeviceToHost )
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//int2str( INT( c_status ) )&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)

  ENDIF
#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

end subroutine cuda_memcpy_device_to_host_int4_r2


subroutine cuda_memcpy_async_host_to_device_int4_r2 ( a_host, a_devc, stream )
  integer(int_4) , dimension(:,:), intent(in), target :: a_host
  TYPE( cuda_memory_t ) :: a_devc
  TYPE( cuda_stream_t ) :: stream

  INTEGER :: n
  INTEGER(KIND(cudaSuccess)) :: c_status
  INTEGER( C_SIZE_T ) :: c_n_bytes
  CHARACTER(*), PARAMETER :: procedureN = 'cuda_memcpy_async_host_to_device_int4_r2'

#if defined(_HAS_CUDA)

  IF( .NOT.a_devc%init ) CALL stopgm(procedureN,'memory not allocated',&
       __LINE__,__FILE__)
  IF( .NOT. stream%init ) CALL stopgm(procedureN,'stream not initialized',&
       __LINE__,__FILE__)
  IF( a_devc%device /= stream%device ) CALL stopgm(procedureN,'devices not consistent',&
       __LINE__,__FILE__)
  CALL cuda_check_device ( a_devc%device, procedureN )

  n = SIZE( a_host )

  IF( n > 0 ) THEN

     c_n_bytes = INT( n, C_SIZE_T )* INT( sizeof_int_4, C_SIZE_T )

     c_status = cudaMemcpyAsync( a_devc%ptr, C_LOC( a_host ), c_n_bytes, cudaMemcpyHostToDevice, stream%s )
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//int2str( INT( c_status ) )&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)

  ENDIF
#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

end subroutine cuda_memcpy_async_host_to_device_int4_r2


subroutine cuda_memcpy_async_device_to_host_int4_r2 ( b_devc, b_host, stream )
  TYPE( cuda_memory_t ) :: b_devc
  integer(int_4) , dimension(:,:), intent(inout), target :: b_host
  TYPE( cuda_stream_t ) :: stream

  INTEGER :: n
  INTEGER(KIND(cudaSuccess)) :: c_status
  INTEGER( C_SIZE_T ) :: c_n_bytes
  CHARACTER(*), PARAMETER :: procedureN = 'cuda_memcpy_async_device_to_host_int4_r2'

#if defined(_HAS_CUDA)

  IF( .NOT.b_devc%init ) CALL stopgm(procedureN,'memory not allocated',&
       __LINE__,__FILE__)
  IF( .NOT. stream%init ) CALL stopgm(procedureN,'stream not initialized',&
       __LINE__,__FILE__)
  IF( b_devc%device /= stream%device ) CALL stopgm(procedureN,'devices not consistent',&
       __LINE__,__FILE__)
  CALL cuda_check_device ( b_devc%device, procedureN )

  n = SIZE( b_host )

  IF( n > 0 ) THEN

     c_n_bytes = INT( n, C_SIZE_T )* INT( sizeof_int_4, C_SIZE_T )

     c_status = cudaMemcpyAsync( C_LOC( b_host ), b_devc%ptr, c_n_bytes, cudaMemcpyDeviceToHost, stream%s )
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//int2str( INT( c_status ) )&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)

  ENDIF
#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

end subroutine cuda_memcpy_async_device_to_host_int4_r2


subroutine cuda_memcpy_host_to_device_real8_r1 ( a_host, a_devc )
  real(real_8) , dimension(:), intent(in), target :: a_host
  TYPE( cuda_memory_t ) :: a_devc

  INTEGER :: n
  INTEGER(KIND(cudaSuccess)) :: c_status
  INTEGER( C_SIZE_T ) :: c_n_bytes
  CHARACTER(*), PARAMETER :: procedureN = 'cuda_memcpy_host_to_device_real8_r1'

#if defined(_HAS_CUDA)

  IF( .NOT.a_devc%init ) CALL stopgm(procedureN,'memory not allocated',&
       __LINE__,__FILE__)

  CALL cuda_check_device ( a_devc%device, procedureN )

  n = SIZE( a_host )

  IF( n > 0 ) THEN

     c_n_bytes = INT( n, C_SIZE_T )* INT( sizeof_real_8, C_SIZE_T )

     c_status = cudaMemcpy( a_devc%ptr, C_LOC( a_host ), c_n_bytes, cudaMemcpyHostToDevice )
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//int2str( INT( c_status ) )&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)

  ENDIF
#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

end subroutine cuda_memcpy_host_to_device_real8_r1


subroutine cuda_memcpy_device_to_host_real8_r1 ( b_devc, b_host )
  real(real_8) , dimension(:), intent(inout), target :: b_host
  TYPE( cuda_memory_t ) :: b_devc

  INTEGER :: n
  INTEGER(KIND(cudaSuccess)) :: c_status
  INTEGER( C_SIZE_T ) :: c_n_bytes
  CHARACTER(*), PARAMETER :: procedureN = 'cuda_memcpy_device_to_host_real8_r1'

#if defined(_HAS_CUDA)

  IF( .NOT.b_devc%init ) CALL stopgm(procedureN,'memory not allocated',&
       __LINE__,__FILE__)

  CALL cuda_check_device ( b_devc%device, procedureN )

  n = SIZE( b_host )

  IF( n > 0 ) THEN

     c_n_bytes = INT( n, C_SIZE_T )* INT( sizeof_real_8, C_SIZE_T )

     c_status = cudaMemcpy( C_LOC( b_host ), b_devc%ptr, c_n_bytes, cudaMemcpyDeviceToHost )
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//int2str( INT( c_status ) )&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)

  ENDIF
#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

end subroutine cuda_memcpy_device_to_host_real8_r1


subroutine cuda_memcpy_async_host_to_device_real8_r1 ( a_host, a_devc, stream )
  real(real_8) , dimension(:), intent(in), target :: a_host
  TYPE( cuda_memory_t ) :: a_devc
  TYPE( cuda_stream_t ) :: stream

  INTEGER :: n
  INTEGER(KIND(cudaSuccess)) :: c_status
  INTEGER( C_SIZE_T ) :: c_n_bytes
  CHARACTER(*), PARAMETER :: procedureN = 'cuda_memcpy_async_host_to_device_real8_r1'

#if defined(_HAS_CUDA)

  IF( .NOT.a_devc%init ) CALL stopgm(procedureN,'memory not allocated',&
       __LINE__,__FILE__)
  IF( .NOT. stream%init ) CALL stopgm(procedureN,'stream not initialized',&
       __LINE__,__FILE__)
  IF( a_devc%device /= stream%device ) CALL stopgm(procedureN,'devices not consistent',&
       __LINE__,__FILE__)
  CALL cuda_check_device ( a_devc%device, procedureN )

  n = SIZE( a_host )

  IF( n > 0 ) THEN

     c_n_bytes = INT( n, C_SIZE_T )* INT( sizeof_real_8, C_SIZE_T )

     c_status = cudaMemcpyAsync( a_devc%ptr, C_LOC( a_host ), c_n_bytes, cudaMemcpyHostToDevice, stream%s )
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//int2str( INT( c_status ) )&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)

  ENDIF
#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

end subroutine cuda_memcpy_async_host_to_device_real8_r1


subroutine cuda_memcpy_async_device_to_host_real8_r1 ( b_devc, b_host, stream )
  TYPE( cuda_memory_t ) :: b_devc
  real(real_8) , dimension(:), intent(inout), target :: b_host
  TYPE( cuda_stream_t ) :: stream

  INTEGER :: n
  INTEGER(KIND(cudaSuccess)) :: c_status
  INTEGER( C_SIZE_T ) :: c_n_bytes
  CHARACTER(*), PARAMETER :: procedureN = 'cuda_memcpy_async_device_to_host_real8_r1'

#if defined(_HAS_CUDA)

  IF( .NOT.b_devc%init ) CALL stopgm(procedureN,'memory not allocated',&
       __LINE__,__FILE__)
  IF( .NOT. stream%init ) CALL stopgm(procedureN,'stream not initialized',&
       __LINE__,__FILE__)
  IF( b_devc%device /= stream%device ) CALL stopgm(procedureN,'devices not consistent',&
       __LINE__,__FILE__)
  CALL cuda_check_device ( b_devc%device, procedureN )

  n = SIZE( b_host )

  IF( n > 0 ) THEN

     c_n_bytes = INT( n, C_SIZE_T )* INT( sizeof_real_8, C_SIZE_T )

     c_status = cudaMemcpyAsync( C_LOC( b_host ), b_devc%ptr, c_n_bytes, cudaMemcpyDeviceToHost, stream%s )
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//int2str( INT( c_status ) )&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)

  ENDIF
#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

end subroutine cuda_memcpy_async_device_to_host_real8_r1


subroutine cuda_memcpy_host_to_device_real8_r2 ( a_host, a_devc )
  real(real_8) , dimension(:,:), intent(in), target :: a_host
  TYPE( cuda_memory_t ) :: a_devc

  INTEGER :: n
  INTEGER(KIND(cudaSuccess)) :: c_status
  INTEGER( C_SIZE_T ) :: c_n_bytes
  CHARACTER(*), PARAMETER :: procedureN = 'cuda_memcpy_host_to_device_real8_r2'

#if defined(_HAS_CUDA)

  IF( .NOT.a_devc%init ) CALL stopgm(procedureN,'memory not allocated',&
       __LINE__,__FILE__)

  CALL cuda_check_device ( a_devc%device, procedureN )

  n = SIZE( a_host )

  IF( n > 0 ) THEN

     c_n_bytes = INT( n, C_SIZE_T )* INT( sizeof_real_8, C_SIZE_T )

     c_status = cudaMemcpy( a_devc%ptr, C_LOC( a_host ), c_n_bytes, cudaMemcpyHostToDevice )
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//int2str( INT( c_status ) )&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)

  ENDIF
#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

end subroutine cuda_memcpy_host_to_device_real8_r2


subroutine cuda_memcpy_device_to_host_real8_r2 ( b_devc, b_host )
  real(real_8) , dimension(:,:), intent(inout), target :: b_host
  TYPE( cuda_memory_t ) :: b_devc

  INTEGER :: n
  INTEGER(KIND(cudaSuccess)) :: c_status
  INTEGER( C_SIZE_T ) :: c_n_bytes
  CHARACTER(*), PARAMETER :: procedureN = 'cuda_memcpy_device_to_host_real8_r2'

#if defined(_HAS_CUDA)

  IF( .NOT.b_devc%init ) CALL stopgm(procedureN,'memory not allocated',&
       __LINE__,__FILE__)

  CALL cuda_check_device ( b_devc%device, procedureN )

  n = SIZE( b_host )

  IF( n > 0 ) THEN

     c_n_bytes = INT( n, C_SIZE_T )* INT( sizeof_real_8, C_SIZE_T )

     c_status = cudaMemcpy( C_LOC( b_host ), b_devc%ptr, c_n_bytes, cudaMemcpyDeviceToHost )
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//int2str( INT( c_status ) )&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)

  ENDIF
#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

end subroutine cuda_memcpy_device_to_host_real8_r2


subroutine cuda_memcpy_async_host_to_device_real8_r2 ( a_host, a_devc, stream )
  real(real_8) , dimension(:,:), intent(in), target :: a_host
  TYPE( cuda_memory_t ) :: a_devc
  TYPE( cuda_stream_t ) :: stream

  INTEGER :: n
  INTEGER(KIND(cudaSuccess)) :: c_status
  INTEGER( C_SIZE_T ) :: c_n_bytes
  CHARACTER(*), PARAMETER :: procedureN = 'cuda_memcpy_async_host_to_device_real8_r2'

#if defined(_HAS_CUDA)

  IF( .NOT.a_devc%init ) CALL stopgm(procedureN,'memory not allocated',&
       __LINE__,__FILE__)
  IF( .NOT. stream%init ) CALL stopgm(procedureN,'stream not initialized',&
       __LINE__,__FILE__)
  IF( a_devc%device /= stream%device ) CALL stopgm(procedureN,'devices not consistent',&
       __LINE__,__FILE__)
  CALL cuda_check_device ( a_devc%device, procedureN )

  n = SIZE( a_host )

  IF( n > 0 ) THEN

     c_n_bytes = INT( n, C_SIZE_T )* INT( sizeof_real_8, C_SIZE_T )

     c_status = cudaMemcpyAsync( a_devc%ptr, C_LOC( a_host ), c_n_bytes, cudaMemcpyHostToDevice, stream%s )
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//int2str( INT( c_status ) )&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)

  ENDIF
#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

end subroutine cuda_memcpy_async_host_to_device_real8_r2


subroutine cuda_memcpy_async_device_to_host_real8_r2 ( b_devc, b_host, stream )
  TYPE( cuda_memory_t ) :: b_devc
  real(real_8) , dimension(:,:), intent(inout), target :: b_host
  TYPE( cuda_stream_t ) :: stream

  INTEGER :: n
  INTEGER(KIND(cudaSuccess)) :: c_status
  INTEGER( C_SIZE_T ) :: c_n_bytes
  CHARACTER(*), PARAMETER :: procedureN = 'cuda_memcpy_async_device_to_host_real8_r2'

#if defined(_HAS_CUDA)

  IF( .NOT.b_devc%init ) CALL stopgm(procedureN,'memory not allocated',&
       __LINE__,__FILE__)
  IF( .NOT. stream%init ) CALL stopgm(procedureN,'stream not initialized',&
       __LINE__,__FILE__)
  IF( b_devc%device /= stream%device ) CALL stopgm(procedureN,'devices not consistent',&
       __LINE__,__FILE__)
  CALL cuda_check_device ( b_devc%device, procedureN )

  n = SIZE( b_host )

  IF( n > 0 ) THEN

     c_n_bytes = INT( n, C_SIZE_T )* INT( sizeof_real_8, C_SIZE_T )

     c_status = cudaMemcpyAsync( C_LOC( b_host ), b_devc%ptr, c_n_bytes, cudaMemcpyDeviceToHost, stream%s )
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//int2str( INT( c_status ) )&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)

  ENDIF
#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

end subroutine cuda_memcpy_async_device_to_host_real8_r2


subroutine cuda_memcpy_host_to_device_complex8_r1 ( a_host, a_devc )
  complex(real_8) , dimension(:), intent(in), target :: a_host
  TYPE( cuda_memory_t ) :: a_devc

  INTEGER :: n
  INTEGER(KIND(cudaSuccess)) :: c_status
  INTEGER( C_SIZE_T ) :: c_n_bytes
  CHARACTER(*), PARAMETER :: procedureN = 'cuda_memcpy_host_to_device_complex8_r1'

#if defined(_HAS_CUDA)

  IF( .NOT.a_devc%init ) CALL stopgm(procedureN,'memory not allocated',&
       __LINE__,__FILE__)

  CALL cuda_check_device ( a_devc%device, procedureN )

  n = SIZE( a_host )

  IF( n > 0 ) THEN

     c_n_bytes = INT( n, C_SIZE_T )* INT( sizeof_complex_8, C_SIZE_T )

     c_status = cudaMemcpy( a_devc%ptr, C_LOC( a_host ), c_n_bytes, cudaMemcpyHostToDevice )
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//int2str( INT( c_status ) )&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)

  ENDIF
#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

end subroutine cuda_memcpy_host_to_device_complex8_r1


subroutine cuda_memcpy_device_to_host_complex8_r1 ( b_devc, b_host )
  complex(real_8) , dimension(:), intent(inout), target :: b_host
  TYPE( cuda_memory_t ) :: b_devc

  INTEGER :: n
  INTEGER(KIND(cudaSuccess)) :: c_status
  INTEGER( C_SIZE_T ) :: c_n_bytes
  CHARACTER(*), PARAMETER :: procedureN = 'cuda_memcpy_device_to_host_complex8_r1'

#if defined(_HAS_CUDA)

  IF( .NOT.b_devc%init ) CALL stopgm(procedureN,'memory not allocated',&
       __LINE__,__FILE__)

  CALL cuda_check_device ( b_devc%device, procedureN )

  n = SIZE( b_host )

  IF( n > 0 ) THEN

     c_n_bytes = INT( n, C_SIZE_T )* INT( sizeof_complex_8, C_SIZE_T )

     c_status = cudaMemcpy( C_LOC( b_host ), b_devc%ptr, c_n_bytes, cudaMemcpyDeviceToHost )
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//int2str( INT( c_status ) )&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)

  ENDIF
#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

end subroutine cuda_memcpy_device_to_host_complex8_r1


subroutine cuda_memcpy_async_host_to_device_complex8_r1 ( a_host, a_devc, stream )
  complex(real_8) , dimension(:), intent(in), target :: a_host
  TYPE( cuda_memory_t ) :: a_devc
  TYPE( cuda_stream_t ) :: stream

  INTEGER :: n
  INTEGER(KIND(cudaSuccess)) :: c_status
  INTEGER( C_SIZE_T ) :: c_n_bytes
  CHARACTER(*), PARAMETER :: procedureN = 'cuda_memcpy_async_host_to_device_complex8_r1'

#if defined(_HAS_CUDA)

  IF( .NOT.a_devc%init ) CALL stopgm(procedureN,'memory not allocated',&
       __LINE__,__FILE__)
  IF( .NOT. stream%init ) CALL stopgm(procedureN,'stream not initialized',&
       __LINE__,__FILE__)
  IF( a_devc%device /= stream%device ) CALL stopgm(procedureN,'devices not consistent',&
       __LINE__,__FILE__)
  CALL cuda_check_device ( a_devc%device, procedureN )

  n = SIZE( a_host )

  IF( n > 0 ) THEN

     c_n_bytes = INT( n, C_SIZE_T )* INT( sizeof_complex_8, C_SIZE_T )

     c_status = cudaMemcpyAsync( a_devc%ptr, C_LOC( a_host ), c_n_bytes, cudaMemcpyHostToDevice, stream%s )
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//int2str( INT( c_status ) )&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)

  ENDIF
#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

end subroutine cuda_memcpy_async_host_to_device_complex8_r1


subroutine cuda_memcpy_async_device_to_host_complex8_r1 ( b_devc, b_host, stream )
  TYPE( cuda_memory_t ) :: b_devc
  complex(real_8) , dimension(:), intent(inout), target :: b_host
  TYPE( cuda_stream_t ) :: stream

  INTEGER :: n
  INTEGER(KIND(cudaSuccess)) :: c_status
  INTEGER( C_SIZE_T ) :: c_n_bytes
  CHARACTER(*), PARAMETER :: procedureN = 'cuda_memcpy_async_device_to_host_complex8_r1'

#if defined(_HAS_CUDA)

  IF( .NOT.b_devc%init ) CALL stopgm(procedureN,'memory not allocated',&
       __LINE__,__FILE__)
  IF( .NOT. stream%init ) CALL stopgm(procedureN,'stream not initialized',&
       __LINE__,__FILE__)
  IF( b_devc%device /= stream%device ) CALL stopgm(procedureN,'devices not consistent',&
       __LINE__,__FILE__)
  CALL cuda_check_device ( b_devc%device, procedureN )

  n = SIZE( b_host )

  IF( n > 0 ) THEN

     c_n_bytes = INT( n, C_SIZE_T )* INT( sizeof_complex_8, C_SIZE_T )

     c_status = cudaMemcpyAsync( C_LOC( b_host ), b_devc%ptr, c_n_bytes, cudaMemcpyDeviceToHost, stream%s )
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//int2str( INT( c_status ) )&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)

  ENDIF
#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

end subroutine cuda_memcpy_async_device_to_host_complex8_r1


subroutine cuda_memcpy_host_to_device_complex8_r2 ( a_host, a_devc )
  complex(real_8) , dimension(:,:), intent(in), target :: a_host
  TYPE( cuda_memory_t ) :: a_devc

  INTEGER :: n
  INTEGER(KIND(cudaSuccess)) :: c_status
  INTEGER( C_SIZE_T ) :: c_n_bytes
  CHARACTER(*), PARAMETER :: procedureN = 'cuda_memcpy_host_to_device_complex8_r2'

#if defined(_HAS_CUDA)

  IF( .NOT.a_devc%init ) CALL stopgm(procedureN,'memory not allocated',&
       __LINE__,__FILE__)

  CALL cuda_check_device ( a_devc%device, procedureN )

  n = SIZE( a_host )

  IF( n > 0 ) THEN

     c_n_bytes = INT( n, C_SIZE_T )* INT( sizeof_complex_8, C_SIZE_T )

     c_status = cudaMemcpy( a_devc%ptr, C_LOC( a_host ), c_n_bytes, cudaMemcpyHostToDevice )
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//int2str( INT( c_status ) )&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)

  ENDIF
#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

end subroutine cuda_memcpy_host_to_device_complex8_r2


subroutine cuda_memcpy_device_to_host_complex8_r2 ( b_devc, b_host )
  complex(real_8) , dimension(:,:), intent(inout), target :: b_host
  TYPE( cuda_memory_t ) :: b_devc

  INTEGER :: n
  INTEGER(KIND(cudaSuccess)) :: c_status
  INTEGER( C_SIZE_T ) :: c_n_bytes
  CHARACTER(*), PARAMETER :: procedureN = 'cuda_memcpy_device_to_host_complex8_r2'

#if defined(_HAS_CUDA)

  IF( .NOT.b_devc%init ) CALL stopgm(procedureN,'memory not allocated',&
       __LINE__,__FILE__)

  CALL cuda_check_device ( b_devc%device, procedureN )

  n = SIZE( b_host )

  IF( n > 0 ) THEN

     c_n_bytes = INT( n, C_SIZE_T )* INT( sizeof_complex_8, C_SIZE_T )

     c_status = cudaMemcpy( C_LOC( b_host ), b_devc%ptr, c_n_bytes, cudaMemcpyDeviceToHost )
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//int2str( INT( c_status ) )&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)

  ENDIF
#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

end subroutine cuda_memcpy_device_to_host_complex8_r2


subroutine cuda_memcpy_async_host_to_device_complex8_r2 ( a_host, a_devc, stream )
  complex(real_8) , dimension(:,:), intent(in), target :: a_host
  TYPE( cuda_memory_t ) :: a_devc
  TYPE( cuda_stream_t ) :: stream

  INTEGER :: n
  INTEGER(KIND(cudaSuccess)) :: c_status
  INTEGER( C_SIZE_T ) :: c_n_bytes
  CHARACTER(*), PARAMETER :: procedureN = 'cuda_memcpy_async_host_to_device_complex8_r2'

#if defined(_HAS_CUDA)

  IF( .NOT.a_devc%init ) CALL stopgm(procedureN,'memory not allocated',&
       __LINE__,__FILE__)
  IF( .NOT. stream%init ) CALL stopgm(procedureN,'stream not initialized',&
       __LINE__,__FILE__)
  IF( a_devc%device /= stream%device ) CALL stopgm(procedureN,'devices not consistent',&
       __LINE__,__FILE__)
  CALL cuda_check_device ( a_devc%device, procedureN )

  n = SIZE( a_host )

  IF( n > 0 ) THEN

     c_n_bytes = INT( n, C_SIZE_T )* INT( sizeof_complex_8, C_SIZE_T )

     c_status = cudaMemcpyAsync( a_devc%ptr, C_LOC( a_host ), c_n_bytes, cudaMemcpyHostToDevice, stream%s )
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//int2str( INT( c_status ) )&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)

  ENDIF
#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

end subroutine cuda_memcpy_async_host_to_device_complex8_r2


subroutine cuda_memcpy_async_device_to_host_complex8_r2 ( b_devc, b_host, stream )
  TYPE( cuda_memory_t ) :: b_devc
  complex(real_8) , dimension(:,:), intent(inout), target :: b_host
  TYPE( cuda_stream_t ) :: stream

  INTEGER :: n
  INTEGER(KIND(cudaSuccess)) :: c_status
  INTEGER( C_SIZE_T ) :: c_n_bytes
  CHARACTER(*), PARAMETER :: procedureN = 'cuda_memcpy_async_device_to_host_complex8_r2'

#if defined(_HAS_CUDA)

  IF( .NOT.b_devc%init ) CALL stopgm(procedureN,'memory not allocated',&
       __LINE__,__FILE__)
  IF( .NOT. stream%init ) CALL stopgm(procedureN,'stream not initialized',&
       __LINE__,__FILE__)
  IF( b_devc%device /= stream%device ) CALL stopgm(procedureN,'devices not consistent',&
       __LINE__,__FILE__)
  CALL cuda_check_device ( b_devc%device, procedureN )

  n = SIZE( b_host )

  IF( n > 0 ) THEN

     c_n_bytes = INT( n, C_SIZE_T )* INT( sizeof_complex_8, C_SIZE_T )

     c_status = cudaMemcpyAsync( C_LOC( b_host ), b_devc%ptr, c_n_bytes, cudaMemcpyDeviceToHost, stream%s )
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//int2str( INT( c_status ) )&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)

  ENDIF
#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

end subroutine cuda_memcpy_async_device_to_host_complex8_r2


