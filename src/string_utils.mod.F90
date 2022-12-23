MODULE string_utils

  USE kinds,                           ONLY: int_8
  USE, INTRINSIC :: ISO_C_BINDING,  ONLY: c_ptr, c_size_t, c_char, c_null_char

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: int2str
  PUBLIC :: str2int
  public :: c_strlen
  PUBLIC :: remove_c_null_char

  INTEGER, PARAMETER, PRIVATE :: max_string_length=12


  INTERFACE int2str
     MODULE PROCEDURE int2str_int4
     MODULE PROCEDURE int2str_int8
  END INTERFACE int2str

  interface 
     function c_strlen( s ) result(result) bind(C,name="strlen")
       import c_ptr, c_size_t
       integer(c_size_t) :: result
       type(c_ptr), value, intent(in) :: s  !character(len=*), intent(in)
     end function c_strlen
  end interface


CONTAINS


  FUNCTION remove_c_null_char ( c_str ) RESULT( reslt )
    CHARACTER(c_char), DIMENSION(:), &
         INTENT(IN)                             :: c_str
    CHARACTER(len=SIZE(c_str))               :: reslt
    
    INTEGER                                  :: i

    reslt = ''
    DO i = 1,SIZE( c_str )
       IF( c_str(i) == C_NULL_CHAR ) EXIT
       reslt(i:i) = c_str(i)
    ENDDO
  END FUNCTION remove_c_null_char


  ! ==================================================================
  FUNCTION int2str_int4(i) RESULT(s)
    INTEGER, INTENT(in)                      :: i
    CHARACTER(len=max_string_length)         :: s

! ==--------------------------------------------------------------==

    WRITE (unit=s,fmt='(I0)') i

    ! ==--------------------------------------------------------------==
  END FUNCTION int2str_int4

  ! ==================================================================
  FUNCTION int2str_int8(i) RESULT(s)
    INTEGER(int_8), INTENT(in)               :: i
    CHARACTER(len=max_string_length)         :: s

! ==--------------------------------------------------------------==

    WRITE (unit=s,fmt='(I0)') i

    ! ==--------------------------------------------------------------==
  END FUNCTION int2str_int8

  ! ==================================================================
  FUNCTION str2int(s) RESULT(i)
    CHARACTER(len=*), INTENT(in)             :: s
    INTEGER                                  :: i

! ==--------------------------------------------------------------==

    READ (unit=s,fmt='(I10)') i

    ! ==--------------------------------------------------------------==
  END FUNCTION str2int

  ! ==================================================================

END MODULE string_utils
