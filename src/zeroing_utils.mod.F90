MODULE zeroing_utils
  USE azzero_utils,                    ONLY: azzero_dont_use_me_anymore,&
                                             i8azzero_dont_use_me_anymore,&
                                             iazzero_dont_use_me_anymore,&
                                             zazzero_dont_use_me_anymore
  USE kinds,                           ONLY: int_4,&
                                             int_8,&
                                             real_8

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: zeroing

  !
  ! interfaces
  !
  INTERFACE zeroing
     MODULE PROCEDURE zeroing_int4_r1
     MODULE PROCEDURE zeroing_int4_r2
     MODULE PROCEDURE zeroing_int4_r3
     MODULE PROCEDURE zeroing_int4_r4
     MODULE PROCEDURE zeroing_int4_r5
     MODULE PROCEDURE zeroing_int4_r6
     MODULE PROCEDURE zeroing_int4_r7

     MODULE PROCEDURE zeroing_int8_r1
     MODULE PROCEDURE zeroing_int8_r2
     MODULE PROCEDURE zeroing_int8_r3
     MODULE PROCEDURE zeroing_int8_r4
     MODULE PROCEDURE zeroing_int8_r5
     MODULE PROCEDURE zeroing_int8_r6
     MODULE PROCEDURE zeroing_int8_r7

     MODULE PROCEDURE zeroing_real8_r1
     MODULE PROCEDURE zeroing_real8_r2
     MODULE PROCEDURE zeroing_real8_r3
     MODULE PROCEDURE zeroing_real8_r4
     MODULE PROCEDURE zeroing_real8_r5
     MODULE PROCEDURE zeroing_real8_r6
     MODULE PROCEDURE zeroing_real8_r7

     MODULE PROCEDURE zeroing_complex8_r1
     MODULE PROCEDURE zeroing_complex8_r2
     MODULE PROCEDURE zeroing_complex8_r3
     MODULE PROCEDURE zeroing_complex8_r4
     MODULE PROCEDURE zeroing_complex8_r5
     MODULE PROCEDURE zeroing_complex8_r6
     MODULE PROCEDURE zeroing_complex8_r7
  END INTERFACE zeroing


CONTAINS


  ! include file for the interfaces

#include "zeroing.inc"


END MODULE zeroing_utils
