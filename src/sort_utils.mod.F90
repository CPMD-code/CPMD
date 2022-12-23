MODULE sort_utils
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: sort2
  PUBLIC :: sort2i

CONTAINS

  SUBROUTINE sort2 ( arr, n, index )


    INTEGER                                  :: n
    REAL(real_8)                             :: arr(1:n)
    INTEGER                                  :: INDEX(1:n)

    INTEGER, PARAMETER                       :: m = 7, nstack = 50  

    INTEGER                                  :: i, ib, ir, istack(1:nstack), &
                                                itemp, j, jstack, k, l
    REAL(real_8)                             :: a, temp

! -----------------------------------------------------------------------

    DO i = 1, n
       INDEX(i) = i
    ENDDO
    jstack = 0
    l = 1
    ir = n
1   IF (ir-l.LT.m) THEN
       DO j = l + 1, ir
          a = arr(j)
          ib = INDEX(j)
          DO i = j - 1, 1, -1
             IF (arr(i).LE.a) go to 2
             arr(i+1) = arr(i)
             INDEX(i+1) = INDEX(i)
          ENDDO
          i = 0
2         arr(i+1) = a
          INDEX(i+1) = ib
       ENDDO
       IF (jstack.EQ.0) RETURN
       ir = istack(jstack)
       l = istack(jstack-1)
       jstack = jstack - 2
    ELSE
       k = (l+ir)/2
       temp = arr(k)
       arr(k) = arr(l+1)
       arr(l+1) = temp
       itemp = INDEX(k)
       INDEX(k) = INDEX(l+1)
       INDEX(l+1) = itemp
       IF (arr(l+1).GT.arr(ir)) THEN
          temp = arr(l+1)
          arr(l+1) = arr(ir)
          arr(ir) = temp
          itemp = INDEX(l+1)
          INDEX(l+1) = INDEX(ir)
          INDEX(ir) = itemp
       ENDIF
       IF (arr(l).GT.arr(ir)) THEN
          temp = arr(l)
          arr(l) = arr(ir)
          arr(ir) = temp
          itemp = INDEX(l)
          INDEX(l) = INDEX(ir)
          INDEX(ir) = itemp
       ENDIF
       IF (arr(l+1).GT.arr(l)) THEN
          temp = arr(l+1)
          arr(l+1) = arr(l)
          arr(l) = temp
          itemp = INDEX(l+1)
          INDEX(l+1) = INDEX(l)
          INDEX(l) = itemp
       ENDIF
       i = l + 1
       j = ir
       a = arr(l)
       ib = INDEX(l)
3      CONTINUE
       i = i + 1
       IF (arr(i).LT.a) go to 3
4      CONTINUE
       j = j - 1
       IF (arr(j).GT.a) go to 4
       IF (j.LT.i) go to 5
       temp = arr(i)
       arr(i) = arr(j)
       arr(j) = temp
       itemp = INDEX(i)
       INDEX(i) = INDEX(j)
       INDEX(j) = itemp
       go to 3
5      arr(l) = arr(j)
       arr(j) = a
       INDEX(l) = INDEX(j)
       INDEX(j) = ib
       jstack = jstack + 2
       IF (jstack.GT.nstack) CALL stopgm('SORT2','Nstack too small',& 
            __LINE__,__FILE__)
       IF (ir-i+1.GE.j-l) THEN
          istack(jstack) = ir
          istack(jstack-1) = i
          ir = j - 1
       ELSE
          istack(jstack) = j - 1
          istack(jstack-1) = l
          l = i
       ENDIF
    ENDIF

    go to 1

  END SUBROUTINE sort2

  ! ******************************************************************************

  SUBROUTINE sort2i ( iarr, n, index )


    INTEGER                                  :: n, iarr(1:n), INDEX(1:n)

    INTEGER, PARAMETER                       :: m = 7, nstack = 50  

    INTEGER                                  :: a, i, ib, ir, &
                                                istack(1:nstack), itemp, j, &
                                                jstack, k, l, temp

! ------------------------------------------------------------------------------

    DO i = 1, n
       INDEX(i) = i
    ENDDO
    jstack = 0
    l = 1
    ir = n
1   IF (ir-l.LT.m) THEN
       DO j = l + 1, ir
          a = iarr(j)
          ib = INDEX(j)
          DO i = j - 1, 1, -1
             IF (iarr(i).LE.a) go to 2
             iarr(i+1) = iarr(i)
             INDEX(i+1) = INDEX(i)
          ENDDO
          i = 0
2         iarr(i+1) = a
          INDEX(i+1) = ib
       ENDDO
       IF (jstack.EQ.0) RETURN
       ir = istack(jstack)
       l = istack(jstack-1)
       jstack = jstack - 2
    ELSE
       k = (l+ir)/2
       temp = iarr(k)
       iarr(k) = iarr(l+1)
       iarr(l+1) = temp
       itemp = INDEX(k)
       INDEX(k) = INDEX(l+1)
       INDEX(l+1) = itemp
       IF (iarr(l+1).GT.iarr(ir)) THEN
          temp = iarr(l+1)
          iarr(l+1) = iarr(ir)
          iarr(ir) = temp
          itemp = INDEX(l+1)
          INDEX(l+1) = INDEX(ir)
          INDEX(ir) = itemp
       ENDIF
       IF (iarr(l).GT.iarr(ir)) THEN
          temp = iarr(l)
          iarr(l) = iarr(ir)
          iarr(ir) = temp
          itemp = INDEX(l)
          INDEX(l) = INDEX(ir)
          INDEX(ir) = itemp
       ENDIF
       IF (iarr(l+1).GT.iarr(l)) THEN
          temp = iarr(l+1)
          iarr(l+1) = iarr(l)
          iarr(l) = temp
          itemp = INDEX(l+1)
          INDEX(l+1) = INDEX(l)
          INDEX(l) = itemp
       ENDIF
       i = l + 1
       j = ir
       a = iarr(l)
       ib = INDEX(l)
3      CONTINUE
       i = i + 1
       IF (iarr(i).LT.a) go to 3
4      CONTINUE
       j = j - 1
       IF (iarr(j).GT.a) go to 4
       IF (j.LT.i) go to 5
       temp = iarr(i)
       iarr(i) = iarr(j)
       iarr(j) = temp
       itemp = INDEX(i)
       INDEX(i) = INDEX(j)
       INDEX(j) = itemp
       go to 3
5      iarr(l) = iarr(j)
       iarr(j) = a
       INDEX(l) = INDEX(j)
       INDEX(j) = ib
       jstack = jstack + 2
       IF (jstack.GT.nstack) CALL stopgm('SORT2I','Nstack too small'&
            ,& 
            __LINE__,__FILE__)
       IF (ir-i+1.GE.j-l) THEN
          istack(jstack) = ir
          istack(jstack-1) = i
          ir = j - 1
       ELSE
          istack(jstack) = j - 1
          istack(jstack-1) = l
          l = i
       ENDIF
    ENDIF

    go to 1

  END SUBROUTINE sort2i

END MODULE sort_utils
