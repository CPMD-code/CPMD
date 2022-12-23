MODULE cp_cudensity_utils

  USE cp_cufft_types,                  ONLY: cp_cufft_stream_get_ptrs,&
                                             cp_cufft_t
  USE cp_curho_types,                  ONLY: cp_curho_stream_get_ptrs,&
                                             cp_curho_t
  USE cuda_types,                      ONLY: cuda_memory_t,&
                                             cuda_stream_t
  USE cuda_utils,                      ONLY: cuda_memcpy_device_to_host
  USE cuuser_utils,                    ONLY: CuUser_Density_Sum
  USE kinds,                           ONLY: real_8
  USE thread_view_types,               ONLY: thread_view_t

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: cp_cubuild_density_sum
  PUBLIC :: cp_cubuild_density_copy_to_host


CONTAINS


  SUBROUTINE cp_cubuild_density_sum ( alpha_real, alpha_imag, n, cp_cufft, cp_curho, thread_view )
    REAL(real_8), INTENT(IN)                 :: alpha_real, alpha_imag
    INTEGER, INTENT(IN)                      :: n
    TYPE(cp_cufft_t), INTENT(IN), TARGET     :: cp_cufft
    TYPE(cp_curho_t), INTENT(INOUT), TARGET  :: cp_curho
    TYPE(thread_view_t), INTENT(IN)          :: thread_view

    INTEGER                                  :: device_idx, stream_idx
    TYPE(cuda_memory_t), POINTER             :: psi_d, rho_d
    TYPE(cuda_stream_t), POINTER             :: stream_p

    device_idx = thread_view%device_idx
    stream_idx = thread_view%stream_idx - 1 !vw FIX that -1

    !vw here t1 shall contain the data !
    CALL cp_cufft_stream_get_ptrs ( cp_cufft, device_idx, stream_idx, stream=stream_p, t1_d=psi_d )
    CALL cp_curho_stream_get_ptrs ( cp_curho, device_idx, stream_idx, rho_d=rho_d )
    CALL CuUser_Density_Sum ( alpha_real, alpha_imag, psi_d, rho_d, n, stream_p )

  END SUBROUTINE cp_cubuild_density_sum


  SUBROUTINE cp_cubuild_density_copy_to_host ( cp_curho, rhoes, thread_views )
    TYPE(cp_curho_t), INTENT(IN), TARGET     :: cp_curho
    REAL(real_8), DIMENSION(:, :, :), &
      INTENT(INOUT)                          :: rhoes
    TYPE(thread_view_t), DIMENSION(:), &
      INTENT(IN)                             :: thread_views

    INTEGER                                  :: device_idx, i_stream, &
                                                stream_idx
    TYPE(cuda_memory_t), POINTER             :: rho_d
    TYPE(thread_view_t)                      :: thread_view

    DO i_stream = 1, SIZE( thread_views )
       thread_view = thread_views( i_stream )
       device_idx = thread_view%device_idx
       stream_idx = thread_view%stream_idx - 1 !vw FIX that -1
       CALL cp_curho_stream_get_ptrs ( cp_curho, device_idx, stream_idx, rho_d=rho_d )
       CALL cuda_memcpy_device_to_host ( rho_d, rhoes( :, 1, i_stream ) )
    ENDDO

  END SUBROUTINE cp_cubuild_density_copy_to_host


END MODULE cp_cudensity_utils
