#if defined(__CHECKS_IEEE_ARITHMETIC)
#define CHECKS_IEEE_ARITHMETIC
#endif

MODULE cp_ieee_interface

#ifdef CHECKS_IEEE_ARITHMETIC
  USE, INTRINSIC :: ieee_arithmetic, ONLY : ieee_is_finite
#endif
  USE kinds, ONLY : real_8

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: cp_ieee_is_finite

  INTERFACE cp_ieee_is_finite
     MODULE PROCEDURE cp_ieee_is_finite_real8_r0
     MODULE PROCEDURE cp_ieee_is_finite_complex8_r0
  END INTERFACE cp_ieee_is_finite

CONTAINS

  FUNCTION cp_ieee_is_finite_real8_r0(x) RESULT(reslt)
    REAL(real_8), INTENT(in)                 :: x
    LOGICAL                                  :: reslt

#ifdef CHECKS_IEEE_ARITHMETIC
    reslt = ieee_is_finite( x )
#else
    reslt = .TRUE.
#endif

  END FUNCTION cp_ieee_is_finite_real8_r0

  FUNCTION cp_ieee_is_finite_complex8_r0(z) RESULT(reslt)
    COMPLEX(real_8), INTENT(in)              :: z
    LOGICAL                                  :: reslt

    REAL(real_8)                             :: im, re

    re = REAL(z,kind=real_8)
    im = AIMAG(z)
#ifdef CHECKS_IEEE_ARITHMETIC
    reslt = ieee_is_finite( re ) .AND. ieee_is_finite( im )
#else
    reslt = .TRUE.
#endif

  END FUNCTION cp_ieee_is_finite_complex8_r0

END MODULE cp_ieee_interface
