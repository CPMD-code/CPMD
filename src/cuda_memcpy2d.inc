subroutine cuda_memcpy2d_host_to_device_int4_r1 ( host, lda_host, devc, lda_devc, m, n )
  integer(int_4) , dimension(:), intent(in), target __CONTIGUOUS :: host
  TYPE( cuda_memory_t ) :: devc
  INTEGER, INTENT(IN) :: lda_host, lda_devc, m, n

  INTEGER(KIND(cudaSuccess)) :: c_status
  INTEGER( C_SIZE_T ) :: c_n, c_m_bytes, c_lda_host_bytes, c_lda_devc_bytes
  CHARACTER(*), PARAMETER :: procedureN = 'cuda_memcpy2d_host_to_device_int4_r1'

#if defined(_HAS_CUDA)

  IF( .NOT. devc%init ) CALL stopgm(procedureN,'memory not allocated',&
       __LINE__,__FILE__)
  CALL cuda_check_device ( devc%device, procedureN )

  IF( m * n > 0 ) THEN

     c_m_bytes = INT( m, C_SIZE_T )* INT( sizeof_int_4, C_SIZE_T )
     c_n = INT( n, C_SIZE_T )
     c_lda_host_bytes = INT( lda_host, C_SIZE_T )* INT( sizeof_int_4, C_SIZE_T )
     c_lda_devc_bytes = INT( lda_devc, C_SIZE_T )* INT( sizeof_int_4, C_SIZE_T )

     c_status = cudaMemcpy2D( devc%ptr, c_lda_devc_bytes, C_LOC( host ), c_lda_host_bytes, &
          &                   c_m_bytes, c_n, cudaMemcpyHostToDevice )
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//int2str( INT( c_status ) )&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)

  ENDIF
#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

end subroutine cuda_memcpy2d_host_to_device_int4_r1


subroutine cuda_memcpy2d_device_to_host_int4_r1 ( devc, lda_devc, host, lda_host, m, n )
  integer(int_4) , dimension(:), intent(inout), target __CONTIGUOUS :: host
  TYPE( cuda_memory_t ) :: devc
  INTEGER, INTENT(IN) :: lda_host, lda_devc, m, n

  INTEGER(KIND(cudaSuccess)) :: c_status
  INTEGER( C_SIZE_T ) :: c_n, c_m_bytes, c_lda_host_bytes, c_lda_devc_bytes
  CHARACTER(*), PARAMETER :: procedureN = 'cuda_memcpy2d_device_to_host_int4_r1'

#if defined(_HAS_CUDA)

  IF( .NOT. devc%init ) CALL stopgm(procedureN,'memory not allocated',&
       __LINE__,__FILE__)
  CALL cuda_check_device ( devc%device, procedureN )

  IF( m * n > 0 ) THEN

     c_m_bytes = INT( m, C_SIZE_T )* INT( sizeof_int_4, C_SIZE_T )
     c_n = INT( n, C_SIZE_T )
     c_lda_host_bytes = INT( lda_host, C_SIZE_T )* INT( sizeof_int_4, C_SIZE_T )
     c_lda_devc_bytes = INT( lda_devc, C_SIZE_T )* INT( sizeof_int_4, C_SIZE_T )

     c_status = cudaMemcpy2D( C_LOC( host ), c_lda_host_bytes, devc%ptr, c_lda_devc_bytes, &
          &                   c_m_bytes, c_n, cudaMemcpyDeviceToHost )
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//int2str( INT( c_status ) )&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)

  ENDIF
#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

end subroutine cuda_memcpy2d_device_to_host_int4_r1


subroutine cuda_memcpy2d_async_host_to_device_int4_r1 ( host, lda_host, devc, lda_devc, m, n, stream )
  integer(int_4) , dimension(:), intent(in), target __CONTIGUOUS :: host
  TYPE( cuda_memory_t ) :: devc
  INTEGER, INTENT(IN) :: lda_host, lda_devc, m, n
  TYPE( cuda_stream_t ) :: stream

  INTEGER(KIND(cudaSuccess)) :: c_status
  INTEGER( C_SIZE_T ) :: c_n, c_m_bytes, c_lda_host_bytes, c_lda_devc_bytes
  CHARACTER(*), PARAMETER :: procedureN = 'cuda_memcpy2d_async_host_to_device_int4_r1'

#if defined(_HAS_CUDA)

  IF( .NOT. devc%init ) CALL stopgm(procedureN,'memory not allocated',&
       __LINE__,__FILE__)
  IF( .NOT. stream%init ) CALL stopgm(procedureN,'stream not initialized',&
       __LINE__,__FILE__)
  IF( devc%device /= stream%device ) CALL stopgm(procedureN,'devices not consistent',&
       __LINE__,__FILE__)
  CALL cuda_check_device ( devc%device, procedureN )

  IF( m * n > 0 ) THEN

     c_m_bytes = INT( m, C_SIZE_T )* INT( sizeof_int_4, C_SIZE_T )
     c_n = INT( n, C_SIZE_T )
     c_lda_host_bytes = INT( lda_host, C_SIZE_T )* INT( sizeof_int_4, C_SIZE_T )
     c_lda_devc_bytes = INT( lda_devc, C_SIZE_T )* INT( sizeof_int_4, C_SIZE_T )

     c_status = cudaMemcpy2DAsync( devc%ptr, c_lda_devc_bytes, C_LOC( host ), c_lda_host_bytes, &
          &                   c_m_bytes, c_n, cudaMemcpyHostToDevice, stream%s )
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//int2str( INT( c_status ) )&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)

  ENDIF
#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

end subroutine cuda_memcpy2d_async_host_to_device_int4_r1


subroutine cuda_memcpy2d_async_device_to_host_int4_r1 ( devc, lda_devc, host, lda_host, m, n, stream )
  integer(int_4) , dimension(:), intent(inout), target __CONTIGUOUS :: host
  TYPE( cuda_memory_t ) :: devc
  INTEGER, INTENT(IN) :: lda_host, lda_devc, m, n
  TYPE( cuda_stream_t ) :: stream

  INTEGER(KIND(cudaSuccess)) :: c_status
  INTEGER( C_SIZE_T ) :: c_n, c_m_bytes, c_lda_host_bytes, c_lda_devc_bytes
  CHARACTER(*), PARAMETER :: procedureN = 'cuda_memcpy2d_async_device_to_host_int4_r1'

#if defined(_HAS_CUDA)

  IF( .NOT. devc%init ) CALL stopgm(procedureN,'memory not allocated',&
       __LINE__,__FILE__)
  IF( .NOT. stream%init ) CALL stopgm(procedureN,'stream not initialized',&
       __LINE__,__FILE__)
  IF( devc%device /= stream%device ) CALL stopgm(procedureN,'devices not consistent',&
       __LINE__,__FILE__)
  CALL cuda_check_device ( devc%device, procedureN )
  
  IF( m * n > 0 ) THEN

     c_m_bytes = INT( m, C_SIZE_T )* INT( sizeof_int_4, C_SIZE_T )
     c_n = INT( n, C_SIZE_T )
     c_lda_host_bytes = INT( lda_host, C_SIZE_T )* INT( sizeof_int_4, C_SIZE_T )
     c_lda_devc_bytes = INT( lda_devc, C_SIZE_T )* INT( sizeof_int_4, C_SIZE_T )

     c_status = cudaMemcpy2DAsync( C_LOC( host ), c_lda_host_bytes, devc%ptr, c_lda_devc_bytes, &
          &                   c_m_bytes, c_n, cudaMemcpyDeviceToHost, stream%s )
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//int2str( INT( c_status ) )&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)

  ENDIF
#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

end subroutine cuda_memcpy2d_async_device_to_host_int4_r1
subroutine cuda_memcpy2d_host_to_device_int4_r2 ( host, lda_host, devc, lda_devc, m, n )
  integer(int_4) , dimension(:,:), intent(in), target __CONTIGUOUS :: host
  TYPE( cuda_memory_t ) :: devc
  INTEGER, INTENT(IN) :: lda_host, lda_devc, m, n

  INTEGER(KIND(cudaSuccess)) :: c_status
  INTEGER( C_SIZE_T ) :: c_n, c_m_bytes, c_lda_host_bytes, c_lda_devc_bytes
  CHARACTER(*), PARAMETER :: procedureN = 'cuda_memcpy2d_host_to_device_int4_r2'

#if defined(_HAS_CUDA)

  IF( .NOT. devc%init ) CALL stopgm(procedureN,'memory not allocated',&
       __LINE__,__FILE__)
  CALL cuda_check_device ( devc%device, procedureN )

  IF( m * n > 0 ) THEN

     c_m_bytes = INT( m, C_SIZE_T )* INT( sizeof_int_4, C_SIZE_T )
     c_n = INT( n, C_SIZE_T )
     c_lda_host_bytes = INT( lda_host, C_SIZE_T )* INT( sizeof_int_4, C_SIZE_T )
     c_lda_devc_bytes = INT( lda_devc, C_SIZE_T )* INT( sizeof_int_4, C_SIZE_T )

     c_status = cudaMemcpy2D( devc%ptr, c_lda_devc_bytes, C_LOC( host ), c_lda_host_bytes, &
          &                   c_m_bytes, c_n, cudaMemcpyHostToDevice )
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//int2str( INT( c_status ) )&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)

  ENDIF
#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

end subroutine cuda_memcpy2d_host_to_device_int4_r2


subroutine cuda_memcpy2d_device_to_host_int4_r2 ( devc, lda_devc, host, lda_host, m, n )
  integer(int_4) , dimension(:,:), intent(inout), target __CONTIGUOUS :: host
  TYPE( cuda_memory_t ) :: devc
  INTEGER, INTENT(IN) :: lda_host, lda_devc, m, n

  INTEGER(KIND(cudaSuccess)) :: c_status
  INTEGER( C_SIZE_T ) :: c_n, c_m_bytes, c_lda_host_bytes, c_lda_devc_bytes
  CHARACTER(*), PARAMETER :: procedureN = 'cuda_memcpy2d_device_to_host_int4_r2'

#if defined(_HAS_CUDA)

  IF( .NOT. devc%init ) CALL stopgm(procedureN,'memory not allocated',&
       __LINE__,__FILE__)
  CALL cuda_check_device ( devc%device, procedureN )

  IF( m * n > 0 ) THEN

     c_m_bytes = INT( m, C_SIZE_T )* INT( sizeof_int_4, C_SIZE_T )
     c_n = INT( n, C_SIZE_T )
     c_lda_host_bytes = INT( lda_host, C_SIZE_T )* INT( sizeof_int_4, C_SIZE_T )
     c_lda_devc_bytes = INT( lda_devc, C_SIZE_T )* INT( sizeof_int_4, C_SIZE_T )

     c_status = cudaMemcpy2D( C_LOC( host ), c_lda_host_bytes, devc%ptr, c_lda_devc_bytes, &
          &                   c_m_bytes, c_n, cudaMemcpyDeviceToHost )
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//int2str( INT( c_status ) )&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)

  ENDIF
#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

end subroutine cuda_memcpy2d_device_to_host_int4_r2


subroutine cuda_memcpy2d_async_host_to_device_int4_r2 ( host, lda_host, devc, lda_devc, m, n, stream )
  integer(int_4) , dimension(:,:), intent(in), target __CONTIGUOUS :: host
  TYPE( cuda_memory_t ) :: devc
  INTEGER, INTENT(IN) :: lda_host, lda_devc, m, n
  TYPE( cuda_stream_t ) :: stream

  INTEGER(KIND(cudaSuccess)) :: c_status
  INTEGER( C_SIZE_T ) :: c_n, c_m_bytes, c_lda_host_bytes, c_lda_devc_bytes
  CHARACTER(*), PARAMETER :: procedureN = 'cuda_memcpy2d_async_host_to_device_int4_r2'

#if defined(_HAS_CUDA)

  IF( .NOT. devc%init ) CALL stopgm(procedureN,'memory not allocated',&
       __LINE__,__FILE__)
  IF( .NOT. stream%init ) CALL stopgm(procedureN,'stream not initialized',&
       __LINE__,__FILE__)
  IF( devc%device /= stream%device ) CALL stopgm(procedureN,'devices not consistent',&
       __LINE__,__FILE__)
  CALL cuda_check_device ( devc%device, procedureN )

  IF( m * n > 0 ) THEN

     c_m_bytes = INT( m, C_SIZE_T )* INT( sizeof_int_4, C_SIZE_T )
     c_n = INT( n, C_SIZE_T )
     c_lda_host_bytes = INT( lda_host, C_SIZE_T )* INT( sizeof_int_4, C_SIZE_T )
     c_lda_devc_bytes = INT( lda_devc, C_SIZE_T )* INT( sizeof_int_4, C_SIZE_T )

     c_status = cudaMemcpy2DAsync( devc%ptr, c_lda_devc_bytes, C_LOC( host ), c_lda_host_bytes, &
          &                   c_m_bytes, c_n, cudaMemcpyHostToDevice, stream%s )
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//int2str( INT( c_status ) )&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)

  ENDIF
#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

end subroutine cuda_memcpy2d_async_host_to_device_int4_r2


subroutine cuda_memcpy2d_async_device_to_host_int4_r2 ( devc, lda_devc, host, lda_host, m, n, stream )
  integer(int_4) , dimension(:,:), intent(inout), target __CONTIGUOUS :: host
  TYPE( cuda_memory_t ) :: devc
  INTEGER, INTENT(IN) :: lda_host, lda_devc, m, n
  TYPE( cuda_stream_t ) :: stream

  INTEGER(KIND(cudaSuccess)) :: c_status
  INTEGER( C_SIZE_T ) :: c_n, c_m_bytes, c_lda_host_bytes, c_lda_devc_bytes
  CHARACTER(*), PARAMETER :: procedureN = 'cuda_memcpy2d_async_device_to_host_int4_r2'

#if defined(_HAS_CUDA)

  IF( .NOT. devc%init ) CALL stopgm(procedureN,'memory not allocated',&
       __LINE__,__FILE__)
  IF( .NOT. stream%init ) CALL stopgm(procedureN,'stream not initialized',&
       __LINE__,__FILE__)
  IF( devc%device /= stream%device ) CALL stopgm(procedureN,'devices not consistent',&
       __LINE__,__FILE__)
  CALL cuda_check_device ( devc%device, procedureN )
  
  IF( m * n > 0 ) THEN

     c_m_bytes = INT( m, C_SIZE_T )* INT( sizeof_int_4, C_SIZE_T )
     c_n = INT( n, C_SIZE_T )
     c_lda_host_bytes = INT( lda_host, C_SIZE_T )* INT( sizeof_int_4, C_SIZE_T )
     c_lda_devc_bytes = INT( lda_devc, C_SIZE_T )* INT( sizeof_int_4, C_SIZE_T )

     c_status = cudaMemcpy2DAsync( C_LOC( host ), c_lda_host_bytes, devc%ptr, c_lda_devc_bytes, &
          &                   c_m_bytes, c_n, cudaMemcpyDeviceToHost, stream%s )
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//int2str( INT( c_status ) )&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)

  ENDIF
#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

end subroutine cuda_memcpy2d_async_device_to_host_int4_r2
subroutine cuda_memcpy2d_host_to_device_real8_r1 ( host, lda_host, devc, lda_devc, m, n )
  real(real_8) , dimension(:), intent(in), target __CONTIGUOUS :: host
  TYPE( cuda_memory_t ) :: devc
  INTEGER, INTENT(IN) :: lda_host, lda_devc, m, n

  INTEGER(KIND(cudaSuccess)) :: c_status
  INTEGER( C_SIZE_T ) :: c_n, c_m_bytes, c_lda_host_bytes, c_lda_devc_bytes
  CHARACTER(*), PARAMETER :: procedureN = 'cuda_memcpy2d_host_to_device_real8_r1'

#if defined(_HAS_CUDA)

  IF( .NOT. devc%init ) CALL stopgm(procedureN,'memory not allocated',&
       __LINE__,__FILE__)
  CALL cuda_check_device ( devc%device, procedureN )

  IF( m * n > 0 ) THEN

     c_m_bytes = INT( m, C_SIZE_T )* INT( sizeof_real_8, C_SIZE_T )
     c_n = INT( n, C_SIZE_T )
     c_lda_host_bytes = INT( lda_host, C_SIZE_T )* INT( sizeof_real_8, C_SIZE_T )
     c_lda_devc_bytes = INT( lda_devc, C_SIZE_T )* INT( sizeof_real_8, C_SIZE_T )

     c_status = cudaMemcpy2D( devc%ptr, c_lda_devc_bytes, C_LOC( host ), c_lda_host_bytes, &
          &                   c_m_bytes, c_n, cudaMemcpyHostToDevice )
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//int2str( INT( c_status ) )&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)

  ENDIF
#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

end subroutine cuda_memcpy2d_host_to_device_real8_r1


subroutine cuda_memcpy2d_device_to_host_real8_r1 ( devc, lda_devc, host, lda_host, m, n )
  real(real_8) , dimension(:), intent(inout), target __CONTIGUOUS :: host
  TYPE( cuda_memory_t ) :: devc
  INTEGER, INTENT(IN) :: lda_host, lda_devc, m, n

  INTEGER(KIND(cudaSuccess)) :: c_status
  INTEGER( C_SIZE_T ) :: c_n, c_m_bytes, c_lda_host_bytes, c_lda_devc_bytes
  CHARACTER(*), PARAMETER :: procedureN = 'cuda_memcpy2d_device_to_host_real8_r1'

#if defined(_HAS_CUDA)

  IF( .NOT. devc%init ) CALL stopgm(procedureN,'memory not allocated',&
       __LINE__,__FILE__)
  CALL cuda_check_device ( devc%device, procedureN )

  IF( m * n > 0 ) THEN

     c_m_bytes = INT( m, C_SIZE_T )* INT( sizeof_real_8, C_SIZE_T )
     c_n = INT( n, C_SIZE_T )
     c_lda_host_bytes = INT( lda_host, C_SIZE_T )* INT( sizeof_real_8, C_SIZE_T )
     c_lda_devc_bytes = INT( lda_devc, C_SIZE_T )* INT( sizeof_real_8, C_SIZE_T )

     c_status = cudaMemcpy2D( C_LOC( host ), c_lda_host_bytes, devc%ptr, c_lda_devc_bytes, &
          &                   c_m_bytes, c_n, cudaMemcpyDeviceToHost )
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//int2str( INT( c_status ) )&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)

  ENDIF
#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

end subroutine cuda_memcpy2d_device_to_host_real8_r1


subroutine cuda_memcpy2d_async_host_to_device_real8_r1 ( host, lda_host, devc, lda_devc, m, n, stream )
  real(real_8) , dimension(:), intent(in), target __CONTIGUOUS :: host
  TYPE( cuda_memory_t ) :: devc
  INTEGER, INTENT(IN) :: lda_host, lda_devc, m, n
  TYPE( cuda_stream_t ) :: stream

  INTEGER(KIND(cudaSuccess)) :: c_status
  INTEGER( C_SIZE_T ) :: c_n, c_m_bytes, c_lda_host_bytes, c_lda_devc_bytes
  CHARACTER(*), PARAMETER :: procedureN = 'cuda_memcpy2d_async_host_to_device_real8_r1'

#if defined(_HAS_CUDA)

  IF( .NOT. devc%init ) CALL stopgm(procedureN,'memory not allocated',&
       __LINE__,__FILE__)
  IF( .NOT. stream%init ) CALL stopgm(procedureN,'stream not initialized',&
       __LINE__,__FILE__)
  IF( devc%device /= stream%device ) CALL stopgm(procedureN,'devices not consistent',&
       __LINE__,__FILE__)
  CALL cuda_check_device ( devc%device, procedureN )

  IF( m * n > 0 ) THEN

     c_m_bytes = INT( m, C_SIZE_T )* INT( sizeof_real_8, C_SIZE_T )
     c_n = INT( n, C_SIZE_T )
     c_lda_host_bytes = INT( lda_host, C_SIZE_T )* INT( sizeof_real_8, C_SIZE_T )
     c_lda_devc_bytes = INT( lda_devc, C_SIZE_T )* INT( sizeof_real_8, C_SIZE_T )

     c_status = cudaMemcpy2DAsync( devc%ptr, c_lda_devc_bytes, C_LOC( host ), c_lda_host_bytes, &
          &                   c_m_bytes, c_n, cudaMemcpyHostToDevice, stream%s )
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//int2str( INT( c_status ) )&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)

  ENDIF
#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

end subroutine cuda_memcpy2d_async_host_to_device_real8_r1


subroutine cuda_memcpy2d_async_device_to_host_real8_r1 ( devc, lda_devc, host, lda_host, m, n, stream )
  real(real_8) , dimension(:), intent(inout), target __CONTIGUOUS :: host
  TYPE( cuda_memory_t ) :: devc
  INTEGER, INTENT(IN) :: lda_host, lda_devc, m, n
  TYPE( cuda_stream_t ) :: stream

  INTEGER(KIND(cudaSuccess)) :: c_status
  INTEGER( C_SIZE_T ) :: c_n, c_m_bytes, c_lda_host_bytes, c_lda_devc_bytes
  CHARACTER(*), PARAMETER :: procedureN = 'cuda_memcpy2d_async_device_to_host_real8_r1'

#if defined(_HAS_CUDA)

  IF( .NOT. devc%init ) CALL stopgm(procedureN,'memory not allocated',&
       __LINE__,__FILE__)
  IF( .NOT. stream%init ) CALL stopgm(procedureN,'stream not initialized',&
       __LINE__,__FILE__)
  IF( devc%device /= stream%device ) CALL stopgm(procedureN,'devices not consistent',&
       __LINE__,__FILE__)
  CALL cuda_check_device ( devc%device, procedureN )
  
  IF( m * n > 0 ) THEN

     c_m_bytes = INT( m, C_SIZE_T )* INT( sizeof_real_8, C_SIZE_T )
     c_n = INT( n, C_SIZE_T )
     c_lda_host_bytes = INT( lda_host, C_SIZE_T )* INT( sizeof_real_8, C_SIZE_T )
     c_lda_devc_bytes = INT( lda_devc, C_SIZE_T )* INT( sizeof_real_8, C_SIZE_T )

     c_status = cudaMemcpy2DAsync( C_LOC( host ), c_lda_host_bytes, devc%ptr, c_lda_devc_bytes, &
          &                   c_m_bytes, c_n, cudaMemcpyDeviceToHost, stream%s )
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//int2str( INT( c_status ) )&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)

  ENDIF
#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

end subroutine cuda_memcpy2d_async_device_to_host_real8_r1
subroutine cuda_memcpy2d_host_to_device_real8_r2 ( host, lda_host, devc, lda_devc, m, n )
  real(real_8) , dimension(:,:), intent(in), target __CONTIGUOUS :: host
  TYPE( cuda_memory_t ) :: devc
  INTEGER, INTENT(IN) :: lda_host, lda_devc, m, n

  INTEGER(KIND(cudaSuccess)) :: c_status
  INTEGER( C_SIZE_T ) :: c_n, c_m_bytes, c_lda_host_bytes, c_lda_devc_bytes
  CHARACTER(*), PARAMETER :: procedureN = 'cuda_memcpy2d_host_to_device_real8_r2'

#if defined(_HAS_CUDA)

  IF( .NOT. devc%init ) CALL stopgm(procedureN,'memory not allocated',&
       __LINE__,__FILE__)
  CALL cuda_check_device ( devc%device, procedureN )

  IF( m * n > 0 ) THEN

     c_m_bytes = INT( m, C_SIZE_T )* INT( sizeof_real_8, C_SIZE_T )
     c_n = INT( n, C_SIZE_T )
     c_lda_host_bytes = INT( lda_host, C_SIZE_T )* INT( sizeof_real_8, C_SIZE_T )
     c_lda_devc_bytes = INT( lda_devc, C_SIZE_T )* INT( sizeof_real_8, C_SIZE_T )

     c_status = cudaMemcpy2D( devc%ptr, c_lda_devc_bytes, C_LOC( host ), c_lda_host_bytes, &
          &                   c_m_bytes, c_n, cudaMemcpyHostToDevice )
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//int2str( INT( c_status ) )&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)

  ENDIF
#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

end subroutine cuda_memcpy2d_host_to_device_real8_r2


subroutine cuda_memcpy2d_device_to_host_real8_r2 ( devc, lda_devc, host, lda_host, m, n )
  real(real_8) , dimension(:,:), intent(inout), target __CONTIGUOUS :: host
  TYPE( cuda_memory_t ) :: devc
  INTEGER, INTENT(IN) :: lda_host, lda_devc, m, n

  INTEGER(KIND(cudaSuccess)) :: c_status
  INTEGER( C_SIZE_T ) :: c_n, c_m_bytes, c_lda_host_bytes, c_lda_devc_bytes
  CHARACTER(*), PARAMETER :: procedureN = 'cuda_memcpy2d_device_to_host_real8_r2'

#if defined(_HAS_CUDA)

  IF( .NOT. devc%init ) CALL stopgm(procedureN,'memory not allocated',&
       __LINE__,__FILE__)
  CALL cuda_check_device ( devc%device, procedureN )

  IF( m * n > 0 ) THEN

     c_m_bytes = INT( m, C_SIZE_T )* INT( sizeof_real_8, C_SIZE_T )
     c_n = INT( n, C_SIZE_T )
     c_lda_host_bytes = INT( lda_host, C_SIZE_T )* INT( sizeof_real_8, C_SIZE_T )
     c_lda_devc_bytes = INT( lda_devc, C_SIZE_T )* INT( sizeof_real_8, C_SIZE_T )

     c_status = cudaMemcpy2D( C_LOC( host ), c_lda_host_bytes, devc%ptr, c_lda_devc_bytes, &
          &                   c_m_bytes, c_n, cudaMemcpyDeviceToHost )
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//int2str( INT( c_status ) )&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)

  ENDIF
#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

end subroutine cuda_memcpy2d_device_to_host_real8_r2


subroutine cuda_memcpy2d_async_host_to_device_real8_r2 ( host, lda_host, devc, lda_devc, m, n, stream )
  real(real_8) , dimension(:,:), intent(in), target __CONTIGUOUS :: host
  TYPE( cuda_memory_t ) :: devc
  INTEGER, INTENT(IN) :: lda_host, lda_devc, m, n
  TYPE( cuda_stream_t ) :: stream

  INTEGER(KIND(cudaSuccess)) :: c_status
  INTEGER( C_SIZE_T ) :: c_n, c_m_bytes, c_lda_host_bytes, c_lda_devc_bytes
  CHARACTER(*), PARAMETER :: procedureN = 'cuda_memcpy2d_async_host_to_device_real8_r2'

#if defined(_HAS_CUDA)

  IF( .NOT. devc%init ) CALL stopgm(procedureN,'memory not allocated',&
       __LINE__,__FILE__)
  IF( .NOT. stream%init ) CALL stopgm(procedureN,'stream not initialized',&
       __LINE__,__FILE__)
  IF( devc%device /= stream%device ) CALL stopgm(procedureN,'devices not consistent',&
       __LINE__,__FILE__)
  CALL cuda_check_device ( devc%device, procedureN )

  IF( m * n > 0 ) THEN

     c_m_bytes = INT( m, C_SIZE_T )* INT( sizeof_real_8, C_SIZE_T )
     c_n = INT( n, C_SIZE_T )
     c_lda_host_bytes = INT( lda_host, C_SIZE_T )* INT( sizeof_real_8, C_SIZE_T )
     c_lda_devc_bytes = INT( lda_devc, C_SIZE_T )* INT( sizeof_real_8, C_SIZE_T )

     c_status = cudaMemcpy2DAsync( devc%ptr, c_lda_devc_bytes, C_LOC( host ), c_lda_host_bytes, &
          &                   c_m_bytes, c_n, cudaMemcpyHostToDevice, stream%s )
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//int2str( INT( c_status ) )&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)

  ENDIF
#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

end subroutine cuda_memcpy2d_async_host_to_device_real8_r2


subroutine cuda_memcpy2d_async_device_to_host_real8_r2 ( devc, lda_devc, host, lda_host, m, n, stream )
  real(real_8) , dimension(:,:), intent(inout), target __CONTIGUOUS :: host
  TYPE( cuda_memory_t ) :: devc
  INTEGER, INTENT(IN) :: lda_host, lda_devc, m, n
  TYPE( cuda_stream_t ) :: stream

  INTEGER(KIND(cudaSuccess)) :: c_status
  INTEGER( C_SIZE_T ) :: c_n, c_m_bytes, c_lda_host_bytes, c_lda_devc_bytes
  CHARACTER(*), PARAMETER :: procedureN = 'cuda_memcpy2d_async_device_to_host_real8_r2'

#if defined(_HAS_CUDA)

  IF( .NOT. devc%init ) CALL stopgm(procedureN,'memory not allocated',&
       __LINE__,__FILE__)
  IF( .NOT. stream%init ) CALL stopgm(procedureN,'stream not initialized',&
       __LINE__,__FILE__)
  IF( devc%device /= stream%device ) CALL stopgm(procedureN,'devices not consistent',&
       __LINE__,__FILE__)
  CALL cuda_check_device ( devc%device, procedureN )
  
  IF( m * n > 0 ) THEN

     c_m_bytes = INT( m, C_SIZE_T )* INT( sizeof_real_8, C_SIZE_T )
     c_n = INT( n, C_SIZE_T )
     c_lda_host_bytes = INT( lda_host, C_SIZE_T )* INT( sizeof_real_8, C_SIZE_T )
     c_lda_devc_bytes = INT( lda_devc, C_SIZE_T )* INT( sizeof_real_8, C_SIZE_T )

     c_status = cudaMemcpy2DAsync( C_LOC( host ), c_lda_host_bytes, devc%ptr, c_lda_devc_bytes, &
          &                   c_m_bytes, c_n, cudaMemcpyDeviceToHost, stream%s )
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//int2str( INT( c_status ) )&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)

  ENDIF
#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

end subroutine cuda_memcpy2d_async_device_to_host_real8_r2
subroutine cuda_memcpy2d_host_to_device_complex8_r1 ( host, lda_host, devc, lda_devc, m, n )
  complex(real_8) , dimension(:), intent(in), target __CONTIGUOUS :: host
  TYPE( cuda_memory_t ) :: devc
  INTEGER, INTENT(IN) :: lda_host, lda_devc, m, n

  INTEGER(KIND(cudaSuccess)) :: c_status
  INTEGER( C_SIZE_T ) :: c_n, c_m_bytes, c_lda_host_bytes, c_lda_devc_bytes
  CHARACTER(*), PARAMETER :: procedureN = 'cuda_memcpy2d_host_to_device_complex8_r1'

#if defined(_HAS_CUDA)

  IF( .NOT. devc%init ) CALL stopgm(procedureN,'memory not allocated',&
       __LINE__,__FILE__)
  CALL cuda_check_device ( devc%device, procedureN )

  IF( m * n > 0 ) THEN

     c_m_bytes = INT( m, C_SIZE_T )* INT( sizeof_complex_8, C_SIZE_T )
     c_n = INT( n, C_SIZE_T )
     c_lda_host_bytes = INT( lda_host, C_SIZE_T )* INT( sizeof_complex_8, C_SIZE_T )
     c_lda_devc_bytes = INT( lda_devc, C_SIZE_T )* INT( sizeof_complex_8, C_SIZE_T )

     c_status = cudaMemcpy2D( devc%ptr, c_lda_devc_bytes, C_LOC( host ), c_lda_host_bytes, &
          &                   c_m_bytes, c_n, cudaMemcpyHostToDevice )
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//int2str( INT( c_status ) )&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)

  ENDIF
#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

end subroutine cuda_memcpy2d_host_to_device_complex8_r1


subroutine cuda_memcpy2d_device_to_host_complex8_r1 ( devc, lda_devc, host, lda_host, m, n )
  complex(real_8) , dimension(:), intent(inout), target __CONTIGUOUS :: host
  TYPE( cuda_memory_t ) :: devc
  INTEGER, INTENT(IN) :: lda_host, lda_devc, m, n

  INTEGER(KIND(cudaSuccess)) :: c_status
  INTEGER( C_SIZE_T ) :: c_n, c_m_bytes, c_lda_host_bytes, c_lda_devc_bytes
  CHARACTER(*), PARAMETER :: procedureN = 'cuda_memcpy2d_device_to_host_complex8_r1'

#if defined(_HAS_CUDA)

  IF( .NOT. devc%init ) CALL stopgm(procedureN,'memory not allocated',&
       __LINE__,__FILE__)
  CALL cuda_check_device ( devc%device, procedureN )

  IF( m * n > 0 ) THEN

     c_m_bytes = INT( m, C_SIZE_T )* INT( sizeof_complex_8, C_SIZE_T )
     c_n = INT( n, C_SIZE_T )
     c_lda_host_bytes = INT( lda_host, C_SIZE_T )* INT( sizeof_complex_8, C_SIZE_T )
     c_lda_devc_bytes = INT( lda_devc, C_SIZE_T )* INT( sizeof_complex_8, C_SIZE_T )

     c_status = cudaMemcpy2D( C_LOC( host ), c_lda_host_bytes, devc%ptr, c_lda_devc_bytes, &
          &                   c_m_bytes, c_n, cudaMemcpyDeviceToHost )
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//int2str( INT( c_status ) )&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)

  ENDIF
#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

end subroutine cuda_memcpy2d_device_to_host_complex8_r1


subroutine cuda_memcpy2d_async_host_to_device_complex8_r1 ( host, lda_host, devc, lda_devc, m, n, stream )
  complex(real_8) , dimension(:), intent(in), target __CONTIGUOUS :: host
  TYPE( cuda_memory_t ) :: devc
  INTEGER, INTENT(IN) :: lda_host, lda_devc, m, n
  TYPE( cuda_stream_t ) :: stream

  INTEGER(KIND(cudaSuccess)) :: c_status
  INTEGER( C_SIZE_T ) :: c_n, c_m_bytes, c_lda_host_bytes, c_lda_devc_bytes
  CHARACTER(*), PARAMETER :: procedureN = 'cuda_memcpy2d_async_host_to_device_complex8_r1'

#if defined(_HAS_CUDA)

  IF( .NOT. devc%init ) CALL stopgm(procedureN,'memory not allocated',&
       __LINE__,__FILE__)
  IF( .NOT. stream%init ) CALL stopgm(procedureN,'stream not initialized',&
       __LINE__,__FILE__)
  IF( devc%device /= stream%device ) CALL stopgm(procedureN,'devices not consistent',&
       __LINE__,__FILE__)
  CALL cuda_check_device ( devc%device, procedureN )

  IF( m * n > 0 ) THEN

     c_m_bytes = INT( m, C_SIZE_T )* INT( sizeof_complex_8, C_SIZE_T )
     c_n = INT( n, C_SIZE_T )
     c_lda_host_bytes = INT( lda_host, C_SIZE_T )* INT( sizeof_complex_8, C_SIZE_T )
     c_lda_devc_bytes = INT( lda_devc, C_SIZE_T )* INT( sizeof_complex_8, C_SIZE_T )

     c_status = cudaMemcpy2DAsync( devc%ptr, c_lda_devc_bytes, C_LOC( host ), c_lda_host_bytes, &
          &                   c_m_bytes, c_n, cudaMemcpyHostToDevice, stream%s )
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//int2str( INT( c_status ) )&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)

  ENDIF
#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

end subroutine cuda_memcpy2d_async_host_to_device_complex8_r1


subroutine cuda_memcpy2d_async_device_to_host_complex8_r1 ( devc, lda_devc, host, lda_host, m, n, stream )
  complex(real_8) , dimension(:), intent(inout), target __CONTIGUOUS :: host
  TYPE( cuda_memory_t ) :: devc
  INTEGER, INTENT(IN) :: lda_host, lda_devc, m, n
  TYPE( cuda_stream_t ) :: stream

  INTEGER(KIND(cudaSuccess)) :: c_status
  INTEGER( C_SIZE_T ) :: c_n, c_m_bytes, c_lda_host_bytes, c_lda_devc_bytes
  CHARACTER(*), PARAMETER :: procedureN = 'cuda_memcpy2d_async_device_to_host_complex8_r1'

#if defined(_HAS_CUDA)

  IF( .NOT. devc%init ) CALL stopgm(procedureN,'memory not allocated',&
       __LINE__,__FILE__)
  IF( .NOT. stream%init ) CALL stopgm(procedureN,'stream not initialized',&
       __LINE__,__FILE__)
  IF( devc%device /= stream%device ) CALL stopgm(procedureN,'devices not consistent',&
       __LINE__,__FILE__)
  CALL cuda_check_device ( devc%device, procedureN )
  
  IF( m * n > 0 ) THEN

     c_m_bytes = INT( m, C_SIZE_T )* INT( sizeof_complex_8, C_SIZE_T )
     c_n = INT( n, C_SIZE_T )
     c_lda_host_bytes = INT( lda_host, C_SIZE_T )* INT( sizeof_complex_8, C_SIZE_T )
     c_lda_devc_bytes = INT( lda_devc, C_SIZE_T )* INT( sizeof_complex_8, C_SIZE_T )

     c_status = cudaMemcpy2DAsync( C_LOC( host ), c_lda_host_bytes, devc%ptr, c_lda_devc_bytes, &
          &                   c_m_bytes, c_n, cudaMemcpyDeviceToHost, stream%s )
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//int2str( INT( c_status ) )&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)

  ENDIF
#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

end subroutine cuda_memcpy2d_async_device_to_host_complex8_r1
subroutine cuda_memcpy2d_host_to_device_complex8_r2 ( host, lda_host, devc, lda_devc, m, n )
  complex(real_8) , dimension(:,:), intent(in), target __CONTIGUOUS :: host
  TYPE( cuda_memory_t ) :: devc
  INTEGER, INTENT(IN) :: lda_host, lda_devc, m, n

  INTEGER(KIND(cudaSuccess)) :: c_status
  INTEGER( C_SIZE_T ) :: c_n, c_m_bytes, c_lda_host_bytes, c_lda_devc_bytes
  CHARACTER(*), PARAMETER :: procedureN = 'cuda_memcpy2d_host_to_device_complex8_r2'

#if defined(_HAS_CUDA)

  IF( .NOT. devc%init ) CALL stopgm(procedureN,'memory not allocated',&
       __LINE__,__FILE__)
  CALL cuda_check_device ( devc%device, procedureN )

  IF( m * n > 0 ) THEN

     c_m_bytes = INT( m, C_SIZE_T )* INT( sizeof_complex_8, C_SIZE_T )
     c_n = INT( n, C_SIZE_T )
     c_lda_host_bytes = INT( lda_host, C_SIZE_T )* INT( sizeof_complex_8, C_SIZE_T )
     c_lda_devc_bytes = INT( lda_devc, C_SIZE_T )* INT( sizeof_complex_8, C_SIZE_T )

     c_status = cudaMemcpy2D( devc%ptr, c_lda_devc_bytes, C_LOC( host ), c_lda_host_bytes, &
          &                   c_m_bytes, c_n, cudaMemcpyHostToDevice )
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//int2str( INT( c_status ) )&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)

  ENDIF
#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

end subroutine cuda_memcpy2d_host_to_device_complex8_r2


subroutine cuda_memcpy2d_device_to_host_complex8_r2 ( devc, lda_devc, host, lda_host, m, n )
  complex(real_8) , dimension(:,:), intent(inout), target __CONTIGUOUS :: host
  TYPE( cuda_memory_t ) :: devc
  INTEGER, INTENT(IN) :: lda_host, lda_devc, m, n

  INTEGER(KIND(cudaSuccess)) :: c_status
  INTEGER( C_SIZE_T ) :: c_n, c_m_bytes, c_lda_host_bytes, c_lda_devc_bytes
  CHARACTER(*), PARAMETER :: procedureN = 'cuda_memcpy2d_device_to_host_complex8_r2'

#if defined(_HAS_CUDA)

  IF( .NOT. devc%init ) CALL stopgm(procedureN,'memory not allocated',&
       __LINE__,__FILE__)
  CALL cuda_check_device ( devc%device, procedureN )

  IF( m * n > 0 ) THEN

     c_m_bytes = INT( m, C_SIZE_T )* INT( sizeof_complex_8, C_SIZE_T )
     c_n = INT( n, C_SIZE_T )
     c_lda_host_bytes = INT( lda_host, C_SIZE_T )* INT( sizeof_complex_8, C_SIZE_T )
     c_lda_devc_bytes = INT( lda_devc, C_SIZE_T )* INT( sizeof_complex_8, C_SIZE_T )

     c_status = cudaMemcpy2D( C_LOC( host ), c_lda_host_bytes, devc%ptr, c_lda_devc_bytes, &
          &                   c_m_bytes, c_n, cudaMemcpyDeviceToHost )
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//int2str( INT( c_status ) )&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)

  ENDIF
#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

end subroutine cuda_memcpy2d_device_to_host_complex8_r2


subroutine cuda_memcpy2d_async_host_to_device_complex8_r2 ( host, lda_host, devc, lda_devc, m, n, stream )
  complex(real_8) , dimension(:,:), intent(in), target __CONTIGUOUS :: host
  TYPE( cuda_memory_t ) :: devc
  INTEGER, INTENT(IN) :: lda_host, lda_devc, m, n
  TYPE( cuda_stream_t ) :: stream

  INTEGER(KIND(cudaSuccess)) :: c_status
  INTEGER( C_SIZE_T ) :: c_n, c_m_bytes, c_lda_host_bytes, c_lda_devc_bytes
  CHARACTER(*), PARAMETER :: procedureN = 'cuda_memcpy2d_async_host_to_device_complex8_r2'

#if defined(_HAS_CUDA)

  IF( .NOT. devc%init ) CALL stopgm(procedureN,'memory not allocated',&
       __LINE__,__FILE__)
  IF( .NOT. stream%init ) CALL stopgm(procedureN,'stream not initialized',&
       __LINE__,__FILE__)
  IF( devc%device /= stream%device ) CALL stopgm(procedureN,'devices not consistent',&
       __LINE__,__FILE__)
  CALL cuda_check_device ( devc%device, procedureN )

  IF( m * n > 0 ) THEN

     c_m_bytes = INT( m, C_SIZE_T )* INT( sizeof_complex_8, C_SIZE_T )
     c_n = INT( n, C_SIZE_T )
     c_lda_host_bytes = INT( lda_host, C_SIZE_T )* INT( sizeof_complex_8, C_SIZE_T )
     c_lda_devc_bytes = INT( lda_devc, C_SIZE_T )* INT( sizeof_complex_8, C_SIZE_T )

     c_status = cudaMemcpy2DAsync( devc%ptr, c_lda_devc_bytes, C_LOC( host ), c_lda_host_bytes, &
          &                   c_m_bytes, c_n, cudaMemcpyHostToDevice, stream%s )
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//int2str( INT( c_status ) )&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)

  ENDIF
#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

end subroutine cuda_memcpy2d_async_host_to_device_complex8_r2


subroutine cuda_memcpy2d_async_device_to_host_complex8_r2 ( devc, lda_devc, host, lda_host, m, n, stream )
  complex(real_8) , dimension(:,:), intent(inout), target __CONTIGUOUS :: host
  TYPE( cuda_memory_t ) :: devc
  INTEGER, INTENT(IN) :: lda_host, lda_devc, m, n
  TYPE( cuda_stream_t ) :: stream

  INTEGER(KIND(cudaSuccess)) :: c_status
  INTEGER( C_SIZE_T ) :: c_n, c_m_bytes, c_lda_host_bytes, c_lda_devc_bytes
  CHARACTER(*), PARAMETER :: procedureN = 'cuda_memcpy2d_async_device_to_host_complex8_r2'

#if defined(_HAS_CUDA)

  IF( .NOT. devc%init ) CALL stopgm(procedureN,'memory not allocated',&
       __LINE__,__FILE__)
  IF( .NOT. stream%init ) CALL stopgm(procedureN,'stream not initialized',&
       __LINE__,__FILE__)
  IF( devc%device /= stream%device ) CALL stopgm(procedureN,'devices not consistent',&
       __LINE__,__FILE__)
  CALL cuda_check_device ( devc%device, procedureN )
  
  IF( m * n > 0 ) THEN

     c_m_bytes = INT( m, C_SIZE_T )* INT( sizeof_complex_8, C_SIZE_T )
     c_n = INT( n, C_SIZE_T )
     c_lda_host_bytes = INT( lda_host, C_SIZE_T )* INT( sizeof_complex_8, C_SIZE_T )
     c_lda_devc_bytes = INT( lda_devc, C_SIZE_T )* INT( sizeof_complex_8, C_SIZE_T )

     c_status = cudaMemcpy2DAsync( C_LOC( host ), c_lda_host_bytes, devc%ptr, c_lda_devc_bytes, &
          &                   c_m_bytes, c_n, cudaMemcpyDeviceToHost, stream%s )
     IF( c_status /= cudaSuccess ) CALL stopgm(procedureN,"cuda error: "//int2str( INT( c_status ) )&
         //". "//trim(cuda_get_error_string(c_status)),&
          __LINE__,__FILE__)

  ENDIF
#else

  CALL stopgm(procedureN,"no cuda available",&
       __LINE__,__FILE__)

#endif

end subroutine cuda_memcpy2d_async_device_to_host_complex8_r2
