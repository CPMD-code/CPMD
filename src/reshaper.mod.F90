MODULE reshaper
  USE kinds,                           ONLY: int_1,&
                                             int_2,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8

  IMPLICIT NONE

  INTERFACE reshape_inplace
     MODULE PROCEDURE reshape_inplace_r3_c3

     MODULE PROCEDURE reshape_inplace_c2_c2
     MODULE PROCEDURE reshape_inplace_c2_r3
     MODULE PROCEDURE reshape_inplace_c2_c1
     MODULE PROCEDURE reshape_inplace_c3_c1
     MODULE PROCEDURE reshape_inplace_c1_c3
     MODULE PROCEDURE reshape_inplace_c2_c3
     MODULE PROCEDURE reshape_inplace_c3_c2

     MODULE PROCEDURE reshape_inplace_c1_r1
     MODULE PROCEDURE reshape_inplace_c2_r2
     MODULE PROCEDURE reshape_inplace_c2_r1
     MODULE PROCEDURE reshape_inplace_c1_r2

     MODULE PROCEDURE reshape_inplace_r1_c1
     MODULE PROCEDURE reshape_inplace_r3_r2
     MODULE PROCEDURE reshape_inplace_r2_r3
     MODULE PROCEDURE reshape_inplace_r3_r4
     MODULE PROCEDURE reshape_inplace_r2_r1
     MODULE PROCEDURE reshape_inplace_r2_c1
     MODULE PROCEDURE reshape_inplace_r2_c2
     MODULE PROCEDURE reshape_inplace_r1_r2
     MODULE PROCEDURE reshape_inplace_r1_r3

     MODULE PROCEDURE reshape_inplace_i1_r1
     MODULE PROCEDURE reshape_inplace_i2_i1

     MODULE PROCEDURE reshape_inplace_c41_c82
     MODULE PROCEDURE reshape_inplace_c41_c81

     ! TODO add procedures for other ranks and types
  END INTERFACE reshape_inplace

  INTERFACE type_cast
     MODULE PROCEDURE type_cast_i8_i4
     MODULE PROCEDURE type_cast_i8_i2
     MODULE PROCEDURE type_cast_i8_i1
     MODULE PROCEDURE type_cast_r1_c1
     MODULE PROCEDURE type_cast_c1_r1
     MODULE PROCEDURE type_cast_c8_c4
  END INTERFACE type_cast

  PUBLIC :: reshape_inplace, type_cast

CONTAINS

  SUBROUTINE type_cast_i8_i4(src, size_src, dst)
    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE
    INTEGER(int_8), TARGET, INTENT(in) :: src(*)
    INTEGER, INTENT(in) :: size_src
    INTEGER(int_4), POINTER, INTENT(out) :: dst(:)

    TYPE(c_ptr) :: loc_src

    loc_src = C_LOC(src)
    CALL C_F_POINTER(loc_src, dst, (/size_src/))
  END SUBROUTINE type_cast_i8_i4

  SUBROUTINE type_cast_i8_i2(src, size_src, dst)
    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE
    INTEGER(int_8), TARGET, INTENT(in) :: src(*)
    INTEGER, INTENT(in) :: size_src
    INTEGER(int_2), POINTER, INTENT(out) :: dst(:)

    TYPE(c_ptr) :: loc_src

    loc_src = C_LOC(src)
    CALL C_F_POINTER(loc_src, dst, (/size_src/))
  END SUBROUTINE type_cast_i8_i2

  SUBROUTINE type_cast_i8_i1(src, size_src, dst)
    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE
    INTEGER(int_8), TARGET, INTENT(in) :: src(*)
    INTEGER, INTENT(in) :: size_src
    INTEGER(int_1), POINTER, INTENT(out) :: dst(:)

    TYPE(c_ptr) :: loc_src

    loc_src = C_LOC(src)
    CALL C_F_POINTER(loc_src, dst, (/size_src/))
  END SUBROUTINE type_cast_i8_i1

  SUBROUTINE type_cast_r1_c1(src, size_src, dst)
    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE
    REAL(real_8), TARGET, INTENT(in) :: src(*)
    INTEGER, INTENT(in) :: size_src
    COMPLEX(real_8), POINTER, INTENT(out) :: dst(:)

    TYPE(c_ptr) :: loc_src

    loc_src = C_LOC(src)
    CALL C_F_POINTER(loc_src, dst, (/size_src/2/))
  END SUBROUTINE type_cast_r1_c1

  SUBROUTINE type_cast_c1_r1(src, size_src, dst)
    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE
    COMPLEX(real_8), TARGET, INTENT(in) :: src(*)
    INTEGER, INTENT(in) :: size_src
    REAL(real_8), POINTER, INTENT(out) :: dst(:)

    TYPE(c_ptr) :: loc_src

    loc_src = C_LOC(src)
    CALL C_F_POINTER(loc_src, dst, (/size_src*2/))
  END SUBROUTINE type_cast_c1_r1

  SUBROUTINE type_cast_c8_c4(src, size_src, dst)
    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE
    COMPLEX(real_8), TARGET, INTENT(in) :: src(*)
    INTEGER, INTENT(in) :: size_src
    COMPLEX(real_4), POINTER, INTENT(out) :: dst(:)

    TYPE(c_ptr) :: loc_src

    loc_src = C_LOC(src)
    CALL C_F_POINTER(loc_src, dst, (/size_src*2/))
  END SUBROUTINE type_cast_c8_c4

  SUBROUTINE reshape_inplace_c2_c2(tgt_x, new_shape, ptr_x)
    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE
    COMPLEX(real_8), TARGET, INTENT(in) :: tgt_x(1,*) ! 1 is dummy dimension
    INTEGER, INTENT(in) :: new_shape(:)
    COMPLEX(real_8), POINTER, INTENT(out) :: ptr_x(:,:)

    TYPE(c_ptr) :: loc_x

    loc_x = C_LOC(tgt_x)
    CALL C_F_POINTER(loc_x, ptr_x, new_shape)
  END SUBROUTINE reshape_inplace_c2_c2

  SUBROUTINE reshape_inplace_c2_c1(tgt_x, new_shape, ptr_x)
    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE
    COMPLEX(real_8), TARGET, INTENT(in) :: tgt_x(1,*) ! 1 is dummy dimension
    INTEGER, INTENT(in) :: new_shape(:)
    COMPLEX(real_8), POINTER, INTENT(out) :: ptr_x(:)

    TYPE(c_ptr) :: loc_x

    loc_x = C_LOC(tgt_x)
    CALL C_F_POINTER(loc_x, ptr_x, new_shape)
  END SUBROUTINE reshape_inplace_c2_c1

  SUBROUTINE reshape_inplace_c41_c82(tgt_x, new_shape, ptr_x)
    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE
    COMPLEX(real_8), TARGET, INTENT(in) :: tgt_x(1,*) ! 1 is dummy dimension
    INTEGER, INTENT(in) :: new_shape(:)
    COMPLEX(real_4), POINTER, INTENT(out) :: ptr_x(:)

    TYPE(c_ptr) :: loc_x

    loc_x = C_LOC(tgt_x)
    CALL C_F_POINTER(loc_x, ptr_x, new_shape)
  END SUBROUTINE reshape_inplace_c41_c82

  SUBROUTINE reshape_inplace_c41_c81(tgt_x, new_shape, ptr_x)
    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE
    COMPLEX(real_8), TARGET, INTENT(in) :: tgt_x(*)
    INTEGER, INTENT(in) :: new_shape(:)
    COMPLEX(real_4), POINTER, INTENT(out) :: ptr_x(:)

    TYPE(c_ptr) :: loc_x

    loc_x = C_LOC(tgt_x)
    CALL C_F_POINTER(loc_x, ptr_x, new_shape)
  END SUBROUTINE reshape_inplace_c41_c81

  SUBROUTINE reshape_inplace_c3_c1(tgt_x, new_shape, ptr_x)
    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE
    COMPLEX(real_8), TARGET, INTENT(in) :: tgt_x(1,1,*) ! 1 is dummy dimension
    INTEGER, INTENT(in) :: new_shape(:)
    COMPLEX(real_8), POINTER, INTENT(out) :: ptr_x(:)

    TYPE(c_ptr) :: loc_x

    loc_x = C_LOC(tgt_x)
    CALL C_F_POINTER(loc_x, ptr_x, new_shape)
  END SUBROUTINE reshape_inplace_c3_c1

  SUBROUTINE reshape_inplace_r3_c3(tgt_x, new_shape, ptr_x)
    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE
    REAL(real_8), TARGET, INTENT(in) :: tgt_x(1,1,*)
    INTEGER, INTENT(in) :: new_shape(:)
    COMPLEX(real_8), POINTER, INTENT(out) :: ptr_x(:,:,:)

    TYPE(c_ptr) :: loc_x

    loc_x = C_LOC(tgt_x)
    CALL C_F_POINTER(loc_x, ptr_x, new_shape)
  END SUBROUTINE reshape_inplace_r3_c3

  SUBROUTINE reshape_inplace_c2_r3(tgt_x, new_shape, ptr_x)
    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE
    COMPLEX(real_8), TARGET, INTENT(in) :: tgt_x(1,*)
    INTEGER, INTENT(in) :: new_shape(:)
    REAL(real_8), POINTER, INTENT(out) :: ptr_x(:,:,:)

    TYPE(c_ptr) :: loc_x

    loc_x = C_LOC(tgt_x)
    CALL C_F_POINTER(loc_x, ptr_x, new_shape)
  END SUBROUTINE reshape_inplace_c2_r3

  SUBROUTINE reshape_inplace_c2_c3(tgt_x, new_shape, ptr_x)
    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE
    COMPLEX(real_8), TARGET, INTENT(in) :: tgt_x(1,*)
    INTEGER, INTENT(in) :: new_shape(:)
    COMPLEX(real_8), POINTER, INTENT(out) :: ptr_x(:,:,:)

    TYPE(c_ptr) :: loc_x

    loc_x = C_LOC(tgt_x)
    CALL C_F_POINTER(loc_x, ptr_x, new_shape)
  END SUBROUTINE reshape_inplace_c2_c3

  SUBROUTINE reshape_inplace_c3_c2(tgt_x, new_shape, ptr_x)
    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE
    COMPLEX(real_8), TARGET, INTENT(in) :: tgt_x(1,1,*)
    INTEGER, INTENT(in) :: new_shape(:)
    COMPLEX(real_8), POINTER, INTENT(out) :: ptr_x(:,:)

    TYPE(c_ptr) :: loc_x

    loc_x = C_LOC(tgt_x)
    CALL C_F_POINTER(loc_x, ptr_x, new_shape)
  END SUBROUTINE reshape_inplace_c3_c2

  SUBROUTINE reshape_inplace_c1_r1(tgt_x, new_shape, ptr_x)
    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE
    COMPLEX(real_8), TARGET, INTENT(in) :: tgt_x(*)
    INTEGER, INTENT(in) :: new_shape(:)
    REAL(real_8), POINTER, INTENT(out) :: ptr_x(:)

    TYPE(c_ptr) :: loc_x

    loc_x = C_LOC(tgt_x)
    CALL C_F_POINTER(loc_x, ptr_x, new_shape)
  END SUBROUTINE reshape_inplace_c1_r1

  SUBROUTINE reshape_inplace_i1_r1(tgt_x, new_shape, ptr_x)
    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE
    INTEGER, TARGET, INTENT(in) :: tgt_x(*)
    INTEGER, INTENT(in) :: new_shape(:)
    REAL(real_8), POINTER, INTENT(out) :: ptr_x(:)

    TYPE(c_ptr) :: loc_x

    loc_x = C_LOC(tgt_x)
    CALL C_F_POINTER(loc_x, ptr_x, new_shape)
  END SUBROUTINE reshape_inplace_i1_r1

  SUBROUTINE reshape_inplace_i2_i1(tgt_x, new_shape, ptr_x)
    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE
    INTEGER, TARGET, INTENT(in) :: tgt_x(1,*)
    INTEGER, INTENT(in) :: new_shape(:)
    INTEGER, POINTER, INTENT(out) :: ptr_x(:)

    TYPE(c_ptr) :: loc_x

    loc_x = C_LOC(tgt_x)
    CALL C_F_POINTER(loc_x, ptr_x, new_shape)
  END SUBROUTINE reshape_inplace_i2_i1

  SUBROUTINE reshape_inplace_c1_r2(tgt_x, new_shape, ptr_x)
    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE
    COMPLEX(real_8), TARGET, INTENT(in) :: tgt_x(*)
    INTEGER, INTENT(in) :: new_shape(:)
    REAL(real_8), POINTER, INTENT(out) :: ptr_x(:,:)

    TYPE(c_ptr) :: loc_x

    loc_x = C_LOC(tgt_x)
    CALL C_F_POINTER(loc_x, ptr_x, new_shape)
  END SUBROUTINE reshape_inplace_c1_r2

  SUBROUTINE reshape_inplace_c1_c3(tgt_x, new_shape, ptr_x)
    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE
    COMPLEX(real_8), TARGET, INTENT(in) :: tgt_x(*)
    INTEGER, INTENT(in) :: new_shape(:)
    COMPLEX(real_8), POINTER, INTENT(out) :: ptr_x(:,:,:)

    TYPE(c_ptr) :: loc_x

    loc_x = C_LOC(tgt_x)
    CALL C_F_POINTER(loc_x, ptr_x, new_shape)
  END SUBROUTINE reshape_inplace_c1_c3

  SUBROUTINE reshape_inplace_r1_c1(tgt_x, new_shape, ptr_x)
    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE
    REAL(real_8), TARGET, INTENT(in) :: tgt_x(*)
    INTEGER, INTENT(in) :: new_shape(:)
    COMPLEX(real_8), POINTER, INTENT(out) :: ptr_x(:)

    TYPE(c_ptr) :: loc_x

    loc_x = C_LOC(tgt_x)
    CALL C_F_POINTER(loc_x, ptr_x, new_shape)
  END SUBROUTINE reshape_inplace_r1_c1

  SUBROUTINE reshape_inplace_c2_r2(tgt_x, new_shape, ptr_x)
    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE
    COMPLEX(real_8), TARGET, INTENT(in) :: tgt_x(1,*)
    INTEGER, INTENT(in) :: new_shape(:)
    REAL(real_8), POINTER, INTENT(out) :: ptr_x(:,:)

    TYPE(c_ptr) :: loc_x

    loc_x = C_LOC(tgt_x)
    CALL C_F_POINTER(loc_x, ptr_x, new_shape)
  END SUBROUTINE reshape_inplace_c2_r2

  SUBROUTINE reshape_inplace_c2_r1(tgt_x, new_shape, ptr_x)
    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE
    COMPLEX(real_8), TARGET, INTENT(in) :: tgt_x(1,*)
    INTEGER, INTENT(in) :: new_shape(:)
    REAL(real_8), POINTER, INTENT(out) :: ptr_x(:)

    TYPE(c_ptr) :: loc_x

    loc_x = C_LOC(tgt_x)
    CALL C_F_POINTER(loc_x, ptr_x, new_shape)
  END SUBROUTINE reshape_inplace_c2_r1

  SUBROUTINE reshape_inplace_r2_r3(tgt_x, new_shape, ptr_x)
    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE
    REAL(real_8), TARGET, INTENT(in) :: tgt_x(1,*)
    INTEGER, INTENT(in) :: new_shape(:)
    REAL(real_8), POINTER, INTENT(out) :: ptr_x(:,:,:)

    TYPE(c_ptr) :: loc_x

    loc_x = C_LOC(tgt_x)
    CALL C_F_POINTER(loc_x, ptr_x, new_shape)
  END SUBROUTINE reshape_inplace_r2_r3

  SUBROUTINE reshape_inplace_r3_r4(tgt_x, new_shape, ptr_x)
    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE
    REAL(real_8), TARGET, INTENT(in) :: tgt_x(1,1,*)
    INTEGER, INTENT(in) :: new_shape(:)
    REAL(real_8), POINTER, INTENT(out) :: ptr_x(:,:,:,:)

    TYPE(c_ptr) :: loc_x

    loc_x = C_LOC(tgt_x)
    CALL C_F_POINTER(loc_x, ptr_x, new_shape)
  END SUBROUTINE reshape_inplace_r3_r4

  SUBROUTINE reshape_inplace_r3_r2(tgt_x, new_shape, ptr_x)
    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE
    REAL(real_8), TARGET, INTENT(in) :: tgt_x(1,1,*)
    INTEGER, INTENT(in) :: new_shape(:)
    REAL(real_8), POINTER, INTENT(out) :: ptr_x(:,:)

    TYPE(c_ptr) :: loc_x

    loc_x = C_LOC(tgt_x)
    CALL C_F_POINTER(loc_x, ptr_x, new_shape)
  END SUBROUTINE reshape_inplace_r3_r2

  SUBROUTINE reshape_inplace_r2_r1(tgt_x, new_shape, ptr_x)
    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE
    REAL(real_8), TARGET, INTENT(in) :: tgt_x(1,*)
    INTEGER, INTENT(in) :: new_shape(:)
    REAL(real_8), POINTER, INTENT(out) :: ptr_x(:)

    TYPE(c_ptr) :: loc_x

    loc_x = C_LOC(tgt_x)
    CALL C_F_POINTER(loc_x, ptr_x, new_shape)
  END SUBROUTINE reshape_inplace_r2_r1

  SUBROUTINE reshape_inplace_r2_c1(tgt_x, new_shape, ptr_x)
    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE
    REAL(real_8), TARGET, INTENT(in) :: tgt_x(1,*)
    INTEGER, INTENT(in) :: new_shape(:)
    COMPLEX(real_8), POINTER, INTENT(out) :: ptr_x(:)

    TYPE(c_ptr) :: loc_x

    loc_x = C_LOC(tgt_x)
    CALL C_F_POINTER(loc_x, ptr_x, new_shape)
  END SUBROUTINE reshape_inplace_r2_c1

  SUBROUTINE reshape_inplace_r2_c2(tgt_x, new_shape, ptr_x)
    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE
    REAL(real_8), TARGET, INTENT(in) :: tgt_x(1,*)
    INTEGER, INTENT(in) :: new_shape(:)
    COMPLEX(real_8), POINTER, INTENT(out) :: ptr_x(:,:)

    TYPE(c_ptr) :: loc_x

    loc_x = C_LOC(tgt_x)
    CALL C_F_POINTER(loc_x, ptr_x, new_shape)
  END SUBROUTINE reshape_inplace_r2_c2

  SUBROUTINE reshape_inplace_r1_r2(tgt_x, new_shape, ptr_x)
    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE
    REAL(real_8), TARGET, INTENT(in) :: tgt_x(*)
    INTEGER, INTENT(in) :: new_shape(:)
    REAL(real_8), POINTER, INTENT(out) :: ptr_x(:,:)

    TYPE(c_ptr) :: loc_x

    loc_x = C_LOC(tgt_x)
    CALL C_F_POINTER(loc_x, ptr_x, new_shape)
  END SUBROUTINE reshape_inplace_r1_r2

  SUBROUTINE reshape_inplace_r1_r3(tgt_x, new_shape, ptr_x)
    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE
    REAL(real_8), TARGET, INTENT(in) :: tgt_x(*)
    INTEGER, INTENT(in) :: new_shape(:)
    REAL(real_8), POINTER, INTENT(out) :: ptr_x(:,:,:)

    TYPE(c_ptr) :: loc_x

    loc_x = C_LOC(tgt_x)
    CALL C_F_POINTER(loc_x, ptr_x, new_shape)
  END SUBROUTINE reshape_inplace_r1_r3

END MODULE reshaper
