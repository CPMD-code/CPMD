MODULE min_heap
  USE error_handling,                  ONLY: stopgm

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: heap_new
  PUBLIC :: heap_release
  PUBLIC :: heap_fill
  PUBLIC :: heap_get_first
  PUBLIC :: heap_get
  PUBLIC :: heap_reset_first
  PUBLIC :: heap_reset

  TYPE, PUBLIC :: heap_t
     INTEGER :: n = 0
     INTEGER,DIMENSION(:),ALLOCATABLE :: vals
     INTEGER,DIMENSION(:),ALLOCATABLE :: keys
  END TYPE heap_t

CONTAINS

  SUBROUTINE heap_new(heap,n)
    TYPE(heap_t), INTENT(inout)              :: heap
    INTEGER, INTENT(in)                      :: n

    CHARACTER(*), PARAMETER                  :: procedureN = 'heap_new'

    INTEGER                                  :: ierr

    heap%n = n
    ALLOCATE(heap%vals(n),heap%keys(n),stat=ierr)
    IF (ierr.NE.0) CALL stopgm(procedureN,'Allocation problem',& 
         __LINE__,__FILE__)
  END SUBROUTINE heap_new

  SUBROUTINE heap_release(heap)
    TYPE(heap_t), INTENT(inout)              :: heap

    CHARACTER(*), PARAMETER                  :: procedureN = 'heap_release'

    INTEGER                                  :: ierr

    heap%n = 0
    DEALLOCATE(heap%vals,heap%keys,stat=ierr)
    IF (ierr.NE.0) CALL stopgm(procedureN,'Deallocation problem',& 
         __LINE__,__FILE__)
  END SUBROUTINE heap_release

  SUBROUTINE heap_fill(heap,keys,vals)
    TYPE(heap_t), INTENT(inout)              :: heap
    INTEGER, DIMENSION(:), INTENT(in)        :: keys, vals

    CHARACTER(*), PARAMETER                  :: procedureN = 'heap_fill'

    INTEGER                                  :: heap_size, i

    IF (heap%n.NE.SIZE(vals).OR.heap%n.NE.SIZE(keys))&
         CALL stopgm(procedureN,'size not concistent',& 
         __LINE__,__FILE__)
    heap%vals(:) = vals(:)
    heap%keys(:) = keys(:)
    heap_size = heap%n
    i = heap_size/2
    DO WHILE(i.GE.1)
       CALL heap_adjust_min(heap%keys,heap%vals,i,heap_size)
       i = i-1
    ENDDO
  END SUBROUTINE heap_fill

  SUBROUTINE heap_get_first(heap,key,val)
    TYPE(heap_t), INTENT(inout)              :: heap
    INTEGER, INTENT(out)                     :: key, val

    CHARACTER(*), PARAMETER                  :: procedureN = 'heap_get_first'

    CALL heap_get(heap,1,key,val)
  END SUBROUTINE heap_get_first

  SUBROUTINE heap_get(heap,i,key,val)
    TYPE(heap_t), INTENT(inout)              :: heap
    INTEGER, INTENT(in)                      :: i
    INTEGER, INTENT(out)                     :: key, val

    CHARACTER(*), PARAMETER                  :: procedureN = 'heap_get'

    IF (heap%n.LT.1) CALL stopgm(procedureN,'size not concistent',& 
         __LINE__,__FILE__)
    IF (i.GT.heap%n) CALL stopgm(procedureN,'index not concistent',& 
         __LINE__,__FILE__)
    IF (i.LT.1) CALL stopgm(procedureN,'index not concistent',& 
         __LINE__,__FILE__)
    key = heap%keys(i)
    val = heap%vals(key)
  END SUBROUTINE heap_get

  SUBROUTINE heap_reset_first(heap,val)
    TYPE(heap_t), INTENT(inout)              :: heap
    INTEGER, INTENT(in)                      :: val

    CHARACTER(*), PARAMETER :: procedureN = 'heap_reset_first'

    CALL heap_reset(heap,1,val)
  END SUBROUTINE heap_reset_first

  SUBROUTINE heap_reset(heap,i,val)
    TYPE(heap_t), INTENT(inout)              :: heap
    INTEGER, INTENT(in)                      :: i, val

    CHARACTER(*), PARAMETER                  :: procedureN = 'heap_reset'

    IF (heap%n.LT.1) CALL stopgm(procedureN,'size not concistent',& 
         __LINE__,__FILE__)
    heap%vals(heap%keys(i)) = val
    CALL heap_adjust_min(heap%keys,heap%vals,i,heap%n)
  END SUBROUTINE heap_reset

  SUBROUTINE heap_adjust_min(heap,VALUE,root,size)
    ! Adjusts a heap. Operation also known as siftdown.
    ! The assumption is that the sons of ROOT are already
    ! roots of subtrees that are in heap order. Then it adjusts
    ! HEAP so that ROOT also gets in heap order.
    ! Please note that the HEAP array must have values
    ! in it, i.e. pointers to VALUE.
    INTEGER, DIMENSION(:), INTENT(inout)     :: heap
    INTEGER, DIMENSION(:), INTENT(in)        :: VALUE
    INTEGER, INTENT(in)                      :: root, size

    INTEGER                                  :: e, j, k

    e=heap(root); k=VALUE(heap(root))
    j=2*root
    DO WHILE(j.LE.size)
       IF (j.LT.size) THEN
          IF (VALUE(heap(j)).GT.VALUE(heap(j+1)))j=j+1
       ENDIF
       IF (k.LE.VALUE(heap(j)))EXIT
       heap(j/2) = heap(j)
       j=j*2
    ENDDO
    heap(j/2)=e
  END SUBROUTINE heap_adjust_min

END MODULE min_heap
