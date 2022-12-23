MODULE density_utils
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: build_density_real
  PUBLIC :: build_density_imag
  PUBLIC :: build_density_sum

CONTAINS

  ! ==================================================================
  SUBROUTINE build_density_real(alpha,psi,rho,n)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: alpha
    COMPLEX(real_8)                          :: psi(:)
    REAL(real_8)                             :: rho(:)
    INTEGER                                  :: n

    INTEGER                                  :: l

    !$omp parallel do private(L)
#ifdef __SR8000
    !poption parallel
#endif
#ifdef _vpp_
    !OCL NOVREC
#endif
    DO l=1,n
       rho(l)=rho(l)+alpha*REAL(psi(l))**2
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE build_density_real
  ! ==================================================================
  SUBROUTINE build_density_imag(alpha,psi,rho,n)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: alpha
    COMPLEX(real_8)                          :: psi(:)
    REAL(real_8)                             :: rho(:)
    INTEGER                                  :: n

    INTEGER                                  :: l

    !$omp parallel do private(L)
#ifdef __SR8000
    !poption parallel
#endif
#ifdef _vpp_
    !OCL NOVREC
#endif
    DO l=1,n
       rho(l)=rho(l)+alpha*AIMAG(psi(l))**2
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE build_density_imag
  ! ==================================================================
  SUBROUTINE build_density_sum(alpha_real,alpha_imag,psi,rho,n)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: alpha_real, alpha_imag
    COMPLEX(real_8)                          :: psi(:)
    REAL(real_8)                             :: rho(:)
    INTEGER                                  :: n

    INTEGER                                  :: l

    !$omp parallel do private(L)
#ifdef __SR8000
    !poption parallel
#endif
#ifdef _vpp_
    !OCL NOVREC
#endif
    DO l=1,n
       rho(l)=rho(l)+alpha_real*REAL(psi(l))**2&
            +alpha_imag*AIMAG(psi(l))**2
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE build_density_sum

END MODULE density_utils
