  SUBROUTINE i4_array_alloc_r1 ( this, n, lb )
    INTEGER(int_4), DIMENSION( : ), ALLOCATABLE, INTENT( INOUT ) :: this
    INTEGER, INTENT( IN ) :: n 
    INTEGER, INTENT( IN ), OPTIONAL :: lb 
    INTEGER :: istat, my_lb 
    CHARACTER(*), PARAMETER :: procedureN = 'i4_array_alloc_r1'
    IF(ALLOCATED(this)) CALL stopgm( procedureN, 'array already allocated',&
         __LINE__,__FILE__)
    my_lb = 1
    IF( PRESENT( lb ) ) my_lb = lb
    ALLOCATE( this( my_lb:my_lb+n-1 ), &
         STAT=istat )
    IF(istat/=0) CALL stopgm( procedureN,'allocation problem',&
         __LINE__,__FILE__)
  END SUBROUTINE i4_array_alloc_r1

  SUBROUTINE i4_array_dealloc_r1 ( this )
    INTEGER(int_4), DIMENSION( : ), ALLOCATABLE, INTENT( INOUT ) :: this
    INTEGER :: istat
    CHARACTER(*), PARAMETER :: procedureN = 'i4_array_dealloc_r1'
    IF( ALLOCATED( this ) ) THEN
       DEALLOCATE( this, STAT=istat )
       IF(istat/=0) CALL stopgm( procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
  END SUBROUTINE i4_array_dealloc_r1

  SUBROUTINE i4_array_realloc_r1 ( this, n, keep_data )
    INTEGER(int_4), DIMENSION( : ), ALLOCATABLE, INTENT( INOUT ) :: this
    INTEGER, INTENT( IN ) :: n 
    LOGICAL, INTENT( IN ), OPTIONAL :: keep_data
    INTEGER :: ub , nb 
    INTEGER(int_4), DIMENSION( : ), ALLOCATABLE :: tmp
    LOGICAL :: my_keep_data
    CHARACTER(*), PARAMETER :: procedureN = 'i4_array_realloc_r1'
    my_keep_data = .FALSE.
    IF( PRESENT( keep_data ) ) my_keep_data = keep_data
    IF (ALLOCATED(this)) THEN
       IF( my_keep_data ) THEN
          ub = MIN( n, UBOUND( this, 1 ) )
          nb=MAX( n, UBOUND( this, 1 ) )
          CALL array_alloc ( tmp, nb )
          tmp( 1:ub ) = this( 1:ub )
          CALL MOVE_ALLOC( tmp, this )
       ELSE
          CALL array_dealloc ( this )
          CALL array_alloc ( this, n )
       ENDIF
    ELSE
        CALL array_alloc ( this, n )
    ENDIF
  END SUBROUTINE i4_array_realloc_r1

  SUBROUTINE i4_array_ensure_size_r1 ( this, n, keep_data, growing_factor )
    INTEGER(int_4), DIMENSION( : ), ALLOCATABLE, INTENT( INOUT ) :: this
    INTEGER, INTENT( IN ) :: n 
    LOGICAL, INTENT( IN ), OPTIONAL :: keep_data
    REAL( real_8 ), INTENT( IN ), OPTIONAL :: growing_factor 
    INTEGER :: new_n 
    CHARACTER(*), PARAMETER :: procedureN = 'i4_array_ensure_size_r1'
    IF( ALLOCATED( this ) ) THEN
       IF( SIZE( this ) < n ) THEN
          IF( PRESENT( growing_factor ) ) THEN
             new_n = INT( n * growing_factor )
          ELSE
             new_n = INT( n * array_default_growing_factor )
          ENDIF
          CALL array_realloc ( this, new_n, keep_data=keep_data )
       ENDIF
    ELSE
       CALL array_alloc ( this, n )
    ENDIF
  END SUBROUTINE i4_array_ensure_size_r1
  SUBROUTINE i4_array_alloc_r2 ( this, n, lb )
    INTEGER(int_4), DIMENSION( :,: ), ALLOCATABLE, INTENT( INOUT ) :: this
    INTEGER, INTENT( IN ) :: n (2)
    INTEGER, INTENT( IN ), OPTIONAL :: lb (2)
    INTEGER :: istat, my_lb (2)
    CHARACTER(*), PARAMETER :: procedureN = 'i4_array_alloc_r2'
    IF(ALLOCATED(this)) CALL stopgm( procedureN, 'array already allocated',&
         __LINE__,__FILE__)
    my_lb = 1
    IF( PRESENT( lb ) ) my_lb = lb
    ALLOCATE( this( my_lb(1):my_lb(1)+n(1)-1, my_lb(2):my_lb(2)+n(2)-1 ), &
         STAT=istat )
    IF(istat/=0) CALL stopgm( procedureN,'allocation problem',&
         __LINE__,__FILE__)
  END SUBROUTINE i4_array_alloc_r2

  SUBROUTINE i4_array_dealloc_r2 ( this )
    INTEGER(int_4), DIMENSION( :,: ), ALLOCATABLE, INTENT( INOUT ) :: this
    INTEGER :: istat
    CHARACTER(*), PARAMETER :: procedureN = 'i4_array_dealloc_r2'
    IF( ALLOCATED( this ) ) THEN
       DEALLOCATE( this, STAT=istat )
       IF(istat/=0) CALL stopgm( procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
  END SUBROUTINE i4_array_dealloc_r2

  SUBROUTINE i4_array_realloc_r2 ( this, n, keep_data )
    INTEGER(int_4), DIMENSION( :,: ), ALLOCATABLE, INTENT( INOUT ) :: this
    INTEGER, INTENT( IN ) :: n (2)
    LOGICAL, INTENT( IN ), OPTIONAL :: keep_data
    INTEGER :: ub (2), nb (2)
    INTEGER(int_4), DIMENSION( :,: ), ALLOCATABLE :: tmp
    LOGICAL :: my_keep_data
    CHARACTER(*), PARAMETER :: procedureN = 'i4_array_realloc_r2'
    my_keep_data = .FALSE.
    IF( PRESENT( keep_data ) ) my_keep_data = keep_data
    IF (ALLOCATED(this)) THEN
       IF( my_keep_data ) THEN
          ub(:) = MIN( n(:), UBOUND( this ) )
          nb(:)=MAX( n(:), UBOUND( this ) )
          CALL array_alloc ( tmp, nb )
          tmp( 1:ub(1), 1:ub(2) ) = this( 1:ub(1), 1:ub(2) )
          CALL MOVE_ALLOC( tmp, this )
       ELSE
          CALL array_dealloc ( this )
          CALL array_alloc ( this, n )
       ENDIF
    ELSE
        CALL array_alloc ( this, n )
    ENDIF
  END SUBROUTINE i4_array_realloc_r2

  SUBROUTINE i4_array_ensure_size_r2 ( this, n, keep_data, growing_factor )
    INTEGER(int_4), DIMENSION( :,: ), ALLOCATABLE, INTENT( INOUT ) :: this
    INTEGER, INTENT( IN ) :: n (2)
    LOGICAL, INTENT( IN ), OPTIONAL :: keep_data
    REAL( real_8 ), INTENT( IN ), OPTIONAL :: growing_factor (2)
    INTEGER :: new_n (2)
    CHARACTER(*), PARAMETER :: procedureN = 'i4_array_ensure_size_r2'
    IF( ALLOCATED( this ) ) THEN
       IF( ANY( SHAPE( this ) < n ) ) THEN
          IF( PRESENT( growing_factor ) ) THEN
             new_n = INT( n * growing_factor )
          ELSE
             new_n = INT( n * array_default_growing_factor )
          ENDIF
          CALL array_realloc ( this, new_n, keep_data=keep_data )
       ENDIF
    ELSE
       CALL array_alloc ( this, n )
    ENDIF
  END SUBROUTINE i4_array_ensure_size_r2
  SUBROUTINE l_array_alloc_r1 ( this, n, lb )
    LOGICAL, DIMENSION( : ), ALLOCATABLE, INTENT( INOUT ) :: this
    INTEGER, INTENT( IN ) :: n 
    INTEGER, INTENT( IN ), OPTIONAL :: lb 
    INTEGER :: istat, my_lb 
    CHARACTER(*), PARAMETER :: procedureN = 'l_array_alloc_r1'
    IF(ALLOCATED(this)) CALL stopgm( procedureN, 'array already allocated',&
         __LINE__,__FILE__)
    my_lb = 1
    IF( PRESENT( lb ) ) my_lb = lb
    ALLOCATE( this( my_lb:my_lb+n-1 ), &
         STAT=istat )
    IF(istat/=0) CALL stopgm( procedureN,'allocation problem',&
         __LINE__,__FILE__)
  END SUBROUTINE l_array_alloc_r1

  SUBROUTINE l_array_dealloc_r1 ( this )
    LOGICAL, DIMENSION( : ), ALLOCATABLE, INTENT( INOUT ) :: this
    INTEGER :: istat
    CHARACTER(*), PARAMETER :: procedureN = 'l_array_dealloc_r1'
    IF( ALLOCATED( this ) ) THEN
       DEALLOCATE( this, STAT=istat )
       IF(istat/=0) CALL stopgm( procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
  END SUBROUTINE l_array_dealloc_r1

  SUBROUTINE l_array_realloc_r1 ( this, n, keep_data )
    LOGICAL, DIMENSION( : ), ALLOCATABLE, INTENT( INOUT ) :: this
    INTEGER, INTENT( IN ) :: n 
    LOGICAL, INTENT( IN ), OPTIONAL :: keep_data
    INTEGER :: ub , nb 
    LOGICAL, DIMENSION( : ), ALLOCATABLE :: tmp
    LOGICAL :: my_keep_data
    CHARACTER(*), PARAMETER :: procedureN = 'l_array_realloc_r1'
    my_keep_data = .FALSE.
    IF( PRESENT( keep_data ) ) my_keep_data = keep_data
    IF (ALLOCATED(this)) THEN
       IF( my_keep_data ) THEN
          ub = MIN( n, UBOUND( this, 1 ) )
          nb=MAX( n, UBOUND( this, 1 ) )
          CALL array_alloc ( tmp, nb )
          tmp( 1:ub ) = this( 1:ub )
          CALL MOVE_ALLOC( tmp, this )
       ELSE
          CALL array_dealloc ( this )
          CALL array_alloc ( this, n )
       ENDIF
    ELSE
        CALL array_alloc ( this, n )
    ENDIF
  END SUBROUTINE l_array_realloc_r1

  SUBROUTINE l_array_ensure_size_r1 ( this, n, keep_data, growing_factor )
    LOGICAL, DIMENSION( : ), ALLOCATABLE, INTENT( INOUT ) :: this
    INTEGER, INTENT( IN ) :: n 
    LOGICAL, INTENT( IN ), OPTIONAL :: keep_data
    REAL( real_8 ), INTENT( IN ), OPTIONAL :: growing_factor 
    INTEGER :: new_n 
    CHARACTER(*), PARAMETER :: procedureN = 'l_array_ensure_size_r1'
    IF( ALLOCATED( this ) ) THEN
       IF( SIZE( this ) < n ) THEN
          IF( PRESENT( growing_factor ) ) THEN
             new_n = INT( n * growing_factor )
          ELSE
             new_n = INT( n * array_default_growing_factor )
          ENDIF
          CALL array_realloc ( this, new_n, keep_data=keep_data )
       ENDIF
    ELSE
       CALL array_alloc ( this, n )
    ENDIF
  END SUBROUTINE l_array_ensure_size_r1
  SUBROUTINE l_array_alloc_r2 ( this, n, lb )
    LOGICAL, DIMENSION( :,: ), ALLOCATABLE, INTENT( INOUT ) :: this
    INTEGER, INTENT( IN ) :: n (2)
    INTEGER, INTENT( IN ), OPTIONAL :: lb (2)
    INTEGER :: istat, my_lb (2)
    CHARACTER(*), PARAMETER :: procedureN = 'l_array_alloc_r2'
    IF(ALLOCATED(this)) CALL stopgm( procedureN, 'array already allocated',&
         __LINE__,__FILE__)
    my_lb = 1
    IF( PRESENT( lb ) ) my_lb = lb
    ALLOCATE( this( my_lb(1):my_lb(1)+n(1)-1, my_lb(2):my_lb(2)+n(2)-1 ), &
         STAT=istat )
    IF(istat/=0) CALL stopgm( procedureN,'allocation problem',&
         __LINE__,__FILE__)
  END SUBROUTINE l_array_alloc_r2

  SUBROUTINE l_array_dealloc_r2 ( this )
    LOGICAL, DIMENSION( :,: ), ALLOCATABLE, INTENT( INOUT ) :: this
    INTEGER :: istat
    CHARACTER(*), PARAMETER :: procedureN = 'l_array_dealloc_r2'
    IF( ALLOCATED( this ) ) THEN
       DEALLOCATE( this, STAT=istat )
       IF(istat/=0) CALL stopgm( procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
  END SUBROUTINE l_array_dealloc_r2

  SUBROUTINE l_array_realloc_r2 ( this, n, keep_data )
    LOGICAL, DIMENSION( :,: ), ALLOCATABLE, INTENT( INOUT ) :: this
    INTEGER, INTENT( IN ) :: n (2)
    LOGICAL, INTENT( IN ), OPTIONAL :: keep_data
    INTEGER :: ub (2), nb (2)
    LOGICAL, DIMENSION( :,: ), ALLOCATABLE :: tmp
    LOGICAL :: my_keep_data
    CHARACTER(*), PARAMETER :: procedureN = 'l_array_realloc_r2'
    my_keep_data = .FALSE.
    IF( PRESENT( keep_data ) ) my_keep_data = keep_data
    IF (ALLOCATED(this)) THEN
       IF( my_keep_data ) THEN
          ub(:) = MIN( n(:), UBOUND( this ) )
          nb(:)=MAX( n(:), UBOUND( this ) )
          CALL array_alloc ( tmp, nb )
          tmp( 1:ub(1), 1:ub(2) ) = this( 1:ub(1), 1:ub(2) )
          CALL MOVE_ALLOC( tmp, this )
       ELSE
          CALL array_dealloc ( this )
          CALL array_alloc ( this, n )
       ENDIF
    ELSE
        CALL array_alloc ( this, n )
    ENDIF
  END SUBROUTINE l_array_realloc_r2

  SUBROUTINE l_array_ensure_size_r2 ( this, n, keep_data, growing_factor )
    LOGICAL, DIMENSION( :,: ), ALLOCATABLE, INTENT( INOUT ) :: this
    INTEGER, INTENT( IN ) :: n (2)
    LOGICAL, INTENT( IN ), OPTIONAL :: keep_data
    REAL( real_8 ), INTENT( IN ), OPTIONAL :: growing_factor (2)
    INTEGER :: new_n (2)
    CHARACTER(*), PARAMETER :: procedureN = 'l_array_ensure_size_r2'
    IF( ALLOCATED( this ) ) THEN
       IF( ANY( SHAPE( this ) < n ) ) THEN
          IF( PRESENT( growing_factor ) ) THEN
             new_n = INT( n * growing_factor )
          ELSE
             new_n = INT( n * array_default_growing_factor )
          ENDIF
          CALL array_realloc ( this, new_n, keep_data=keep_data )
       ENDIF
    ELSE
       CALL array_alloc ( this, n )
    ENDIF
  END SUBROUTINE l_array_ensure_size_r2
