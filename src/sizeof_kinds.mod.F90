MODULE sizeof_kinds
  USE kinds,                           ONLY: int_2,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8
#include "sizeof.h"

  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER, PUBLIC :: sizeof_int = size_in_bytes_of( 0 )

  INTEGER, PARAMETER, PUBLIC :: sizeof_int_1 = size_in_bytes_of( 0_int_1 )
  INTEGER, PARAMETER, PUBLIC :: sizeof_int_2 = size_in_bytes_of( 0_int_2 )
  INTEGER, PARAMETER, PUBLIC :: sizeof_int_4 = size_in_bytes_of( 0_int_4 )
  INTEGER, PARAMETER, PUBLIC :: sizeof_int_8 = size_in_bytes_of( 0_int_8 )

  INTEGER, PARAMETER, PUBLIC :: sizeof_real_4 = size_in_bytes_of( 0.0_real_4 )
  INTEGER, PARAMETER, PUBLIC :: sizeof_real_8 = size_in_bytes_of( 0.0_real_8 )

  INTEGER, PARAMETER, PUBLIC :: sizeof_complex_4 = size_in_bytes_of( (0.0_real_4,0.0_real_4) )
  INTEGER, PARAMETER, PUBLIC :: sizeof_complex_8 = size_in_bytes_of( (0.0_real_8,0.0_real_8) )

END MODULE sizeof_kinds
