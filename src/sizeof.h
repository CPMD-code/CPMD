
#if defined (__HAS_SIZEOF)

#define size_in_bytes_of(a) (INT(sizeof(a)))

#else

#define size_in_bytes_of(a) (storage_size(a)/storage_size(1_int_1))

#endif

USE kinds, ONLY : int_1
