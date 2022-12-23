#include "cpmd_global.h"

MODULE lxc_utils

  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8

  USE, INTRINSIC :: iso_c_binding,     ONLY: C_PTR,C_NULL_PTR

#if defined(_HAS_LIBXC)

  USE xc_f03_lib_m, ONLY: xc_f03_func_t, xc_f03_func_info_t, &
       xc_f03_lda_exc_vxc, xc_f03_gga_exc_vxc, xc_f03_func_info_get_family, xc_f03_func_get_info, xc_f03_func_init, &
       xc_f03_functional_get_number, &
       XC_FAMILY_LDA, XC_FAMILY_GGA, XC_FAMILY_HYB_GGA, XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA,  &
       XC_POLARIZED, XC_UNPOLARIZED, XC_GGA_X_PBE, XC_GGA_C_PBE

#endif

  IMPLICIT NONE

#if defined(_HAS_LIBXC)

#else

  TYPE :: xc_f03_func_info_t
     PRIVATE
     TYPE(C_PTR) :: ptr = C_NULL_PTR
  END TYPE xc_f03_func_info_t

  TYPE :: xc_f03_func_t
     PRIVATE
     TYPE(C_PTR) :: ptr = C_NULL_PTR
  END TYPE xc_f03_func_t

  INTEGER, PARAMETER :: XC_FAMILY_LDA = 0
  INTEGER, PARAMETER :: XC_FAMILY_GGA = 1
  INTEGER, PARAMETER :: XC_FAMILY_HYB_GGA = 2
  INTEGER, PARAMETER :: XC_FAMILY_MGGA = 3
  INTEGER, PARAMETER :: XC_FAMILY_HYB_MGGA = 4

  INTEGER, PARAMETER :: XC_UNPOLARIZED = 1
  INTEGER, PARAMETER :: XC_POLARIZED = 2

  INTEGER, PARAMETER :: XC_GGA_X_PBE = 1
  INTEGER, PARAMETER :: XC_GGA_C_PBE = 2

#endif

  PUBLIC :: xc_f03_func_t
  PUBLIC :: xc_f03_func_info_t

  PUBLIC :: XC_FAMILY_LDA
  PUBLIC :: XC_FAMILY_GGA
  PUBLIC :: XC_FAMILY_HYB_GGA
  PUBLIC :: XC_UNPOLARIZED
  PUBLIC :: XC_GGA_X_PBE, XC_GGA_C_PBE

  PUBLIC :: lxc_func_get_info
  PUBLIC :: lxc_func_init
  PUBLIC :: lxc_gga_exc_vxc
  PUBLIC :: lxc_lda_exc_vxc
  PUBLIC :: lxc_func_info_get_family
  PUBLIC :: lxc_functional_get_number


CONTAINS

  FUNCTION lxc_functional_get_number( func_string ) RESULT(number)
    CHARACTER(len=*), INTENT(in)             :: func_string
    INTEGER                                  :: number

    CHARACTER(*), PARAMETER :: procedureN = 'lxc_functional_get_number'

#if defined(_HAS_LIBXC)

    number = xc_f03_functional_get_number ( func_string )

#else

    CALL stopgm(procedureN,"no libxc available",&
         __LINE__,__FILE__)

#endif

  END FUNCTION lxc_functional_get_number


  FUNCTION lxc_func_get_info(p) RESULT(info)
    TYPE(xc_f03_func_t), INTENT(in)          :: p
    TYPE(xc_f03_func_info_t)                 :: info

    CHARACTER(*), PARAMETER :: procedureN = 'lxc_func_get_info'

#if defined(_HAS_LIBXC)

    info = xc_f03_func_get_info(p)

#else

    CALL stopgm(procedureN,"no libxc available",&
         __LINE__,__FILE__)

#endif

  END FUNCTION lxc_func_get_info

  SUBROUTINE lxc_func_init(p,functional,nspin)
    TYPE(xc_f03_func_t), INTENT(inout)       :: p
    INTEGER, INTENT(in)                      :: functional, nspin

    CHARACTER(*), PARAMETER                  :: procedureN = 'lxc_func_init'

#if defined(_HAS_LIBXC)

    CALL xc_f03_func_init(p,functional,nspin)

#else

    CALL stopgm(procedureN,"no libxc available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE lxc_func_init

  SUBROUTINE lxc_lda_exc_vxc(p, np, rho, zk, vrho)
    TYPE(xc_f03_func_t), INTENT(in)          :: p
    INTEGER, INTENT(in)                      :: np
    REAL(real_8), INTENT(in)                 :: rho(*)
    REAL(real_8), INTENT(out)                :: zk(*), vrho(*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'lxc_lda_exc_vxc'

#if defined(_HAS_LIBXC)

    CALL xc_f03_lda_exc_vxc( p, np, rho, zk, vrho )

#else

    CALL stopgm(procedureN,"no libxc available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE lxc_lda_exc_vxc


  SUBROUTINE lxc_gga_exc_vxc(p, np, rho, sigma, zk, vrho, vsigma)
    TYPE(xc_f03_func_t), INTENT(in)          :: p
    INTEGER, INTENT(in)                      :: np
    REAL(real_8), INTENT(in)                 :: rho(*), sigma(*)
    REAL(real_8), INTENT(out)                :: zk(*), vrho(*), vsigma(*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'lxc_gga_exc_vxc'

#if defined(_HAS_LIBXC)

    CALL xc_f03_gga_exc_vxc( p, np, rho, sigma, zk, vrho, vsigma )

#else

    CALL stopgm(procedureN,"no libxc available",&
         __LINE__,__FILE__)

#endif

  END SUBROUTINE lxc_gga_exc_vxc


  FUNCTION lxc_func_info_get_family(info) RESULT(family)
    TYPE(xc_f03_func_info_t), INTENT(in)     :: info
    INTEGER                                  :: family

    CHARACTER(*), PARAMETER :: procedureN = 'lxc_func_info_get_family'

    family = 0
#if defined(_HAS_LIBXC)

    family = xc_f03_func_info_get_family( info )

#else

    CALL stopgm(procedureN,"no libxc available",&
         __LINE__,__FILE__)

#endif

  END FUNCTION lxc_func_info_get_family


END MODULE lxc_utils
