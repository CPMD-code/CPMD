#undef __UNROLL
! #define __SAFE
#undef __SAFE
! ----------------------------------------------------------------------
SUBROUTINE  dscal(n,da,dx,incx)
  ! ----------------------------------------------------------------------
  ! 
  ! scales a vector by a constant.
  ! uses unrolled loops for increment equal to one.
  ! jack dongarra, linpack, 3/11/78.
  ! modified 3/93 to return if incx .le. 0.
  ! modified 12/3/93, array(1) declarations changed to array(*)
  ! 
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  DOUBLE PRECISION da,dx(*)
  INTEGER :: i,incx,m,mp1,n,nincx
  ! 
#ifdef __SAFE
  IF ( n.GT.0 .AND. incx.GT.0 ) THEN
#endif
     IF (incx.EQ.1) THEN
        ! 
        ! code for increment equal to 1
        ! 
        ! 
        ! clean-up loop
        ! 
#ifdef __UNROLL
        m = MOD(n,5)
        IF ( m .NE. 0 ) THEN
           DO i = 1,m
              dx(i) = da*dx(i)
           ENDDO
        ENDIF
        IF ( n .GE. 5 ) THEN
           mp1 = m + 1
           DO i = mp1,n,5
              dx(i) = da*dx(i)
              dx(i + 1) = da*dx(i + 1)
              dx(i + 2) = da*dx(i + 2)
              dx(i + 3) = da*dx(i + 3)
              dx(i + 4) = da*dx(i + 4)
           ENDDO
        ENDIF
#else
        DO i = 1,n
           dx(i) = da*dx(i)
        ENDDO
#endif
     ELSE
        ! 
        ! code for increment not equal to 1
        ! 
        nincx = n*incx
        DO i = 1,nincx,incx
           dx(i) = da*dx(i)
        ENDDO
     ENDIF
#ifdef __SAFE
  ENDIF
#endif
  RETURN
END SUBROUTINE dscal

! ----------------------------------------------------------------------
SUBROUTINE  zdscal(n,da,zx,incx)
  ! ----------------------------------------------------------------------
  ! 
  ! scales a vector by a constant.
  ! jack dongarra, 3/11/78.
  ! modified 3/93 to return if incx .le. 0.
  ! modified 12/3/93, array(1) declarations changed to array(*)
  ! 
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  DOUBLE COMPLEX zx(*)
  DOUBLE PRECISION da
  INTEGER :: i,incx,ix,n
  ! 
#ifdef __SAFE
  IF ( n.LE.0 .OR. incx.LE.0 )RETURN
#endif
  IF (incx.EQ.1)go to 20
  ! 
  ! code for increment not equal to 1
  ! 
  ix = 1
  DO 10 i = 1,n
     zx(ix) = CMPLX(da,0.0_real_8,kind=real_8)*zx(ix)
     ix = ix + incx
10   CONTINUE
     RETURN
     ! 
     ! code for increment equal to 1
     ! 
20   DO 30 i = 1,n
        zx(i) = CMPLX(da,0.0_real_8,kind=real_8)*zx(i)
30      CONTINUE
        RETURN
     ENDDO

     ! ----------------------------------------------------------------------
     SUBROUTINE  dcopy(n,dx,incx,dy,incy)
       ! ----------------------------------------------------------------------
       ! 
       ! copies a vector, x, to a vector, y.
       ! uses unrolled loops for increments equal to one.
       ! jack dongarra, linpack, 3/11/78.
       ! modified 12/3/93, array(1) declarations changed to array(*)
       ! 
       USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
       USE error_handling, ONLY: stopgm
       USE timer, ONLY: tiset, tihalt
       DOUBLE PRECISION dx(*),dy(*)
       INTEGER :: i,incx,incy,ix,iy,m,mp1,n
       ! 
#ifdef __SAFE
       IF (n.GT.0)THEN
#endif
          IF (incx.NE.1.OR.incy.NE.1) THEN
             ! 
             ! code for unequal increments or equal increments
             ! not equal to 1
             ! 
             ix = 1
             iy = 1
             IF (incx.LT.0)ix = (-n+1)*incx + 1
             IF (incy.LT.0)iy = (-n+1)*incy + 1
             DO i = 1,n
                dy(iy) = dx(ix)
                ix = ix + incx
                iy = iy + incy
             ENDDO
          ELSE
             ! 
             ! code for both increments equal to 1
             ! 
             ! 
             ! clean-up loop
             ! 
#ifdef __UNROLL
             m = MOD(n,7)
             IF ( m .NE. 0 ) THEN
                DO i = 1,m
                   dy(i) = dx(i)
                ENDDO
             ENDIF
             IF ( n .GE. 7 ) THEN
                mp1 = m + 1
                DO i = mp1,n,7
                   dy(i) = dx(i)
                   dy(i + 1) = dx(i + 1)
                   dy(i + 2) = dx(i + 2)
                   dy(i + 3) = dx(i + 3)
                   dy(i + 4) = dx(i + 4)
                   dy(i + 5) = dx(i + 5)
                   dy(i + 6) = dx(i + 6)
                ENDDO
             ENDIF
#else
             DO i = 1,n
                dy(i) = dx(i)
             ENDDO
#endif
          ENDIF
#ifdef __SAFE
       ENDIF
#endif
       RETURN
     END SUBROUTINE dcopy

     ! ----------------------------------------------------------------------
     SUBROUTINE  zcopy(n,dx,incx,dy,incy)
       ! ----------------------------------------------------------------------
       ! 
       ! copies a vector, x, to a vector, y.
       ! uses unrolled loops for increments equal to one.
       ! jack dongarra, linpack, 3/11/78.
       ! modified 12/3/93, array(1) declarations changed to array(*)
       ! 
       USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
       USE error_handling, ONLY: stopgm
       USE timer, ONLY: tiset, tihalt
       COMPLEX(real_8) :: dx(*),dy(*)
       INTEGER :: i,incx,incy,ix,iy,m,mp1,n
       ! 
#ifdef __SAFE
       IF (n.GT.0)THEN
#endif
          IF (incx.NE.1.OR.incy.NE.1) THEN
             ! 
             ! code for unequal increments or equal increments
             ! not equal to 1
             ! 
             ix = 1
             iy = 1
             IF (incx.LT.0)ix = (-n+1)*incx + 1
             IF (incy.LT.0)iy = (-n+1)*incy + 1
             DO i = 1,n
                dy(iy) = dx(ix)
                ix = ix + incx
                iy = iy + incy
             ENDDO
          ELSE
             ! 
             ! code for both increments equal to 1
             ! 
             ! 
             ! clean-up loop
             ! 
#ifdef __UNROLL
             m = MOD(n,7)
             IF ( m .NE. 0 ) THEN
                DO i = 1,m
                   dy(i) = dx(i)
                ENDDO
             ENDIF
             IF ( n .GE. 7 ) THEN
                mp1 = m + 1
                DO i = mp1,n,7
                   dy(i) = dx(i)
                   dy(i + 1) = dx(i + 1)
                   dy(i + 2) = dx(i + 2)
                   dy(i + 3) = dx(i + 3)
                   dy(i + 4) = dx(i + 4)
                   dy(i + 5) = dx(i + 5)
                   dy(i + 6) = dx(i + 6)
                ENDDO
             ENDIF
#else
             DO i = 1,n
                dy(i) = dx(i)
             ENDDO
#endif
          ENDIF
#ifdef __SAFE
       ENDIF
#endif
       RETURN
     END SUBROUTINE zcopy

     ! ----------------------------------------------------------------------
     DOUBLE PRECISION FUNCTION ddot(n,dx,incx,dy,incy)
       ! ----------------------------------------------------------------------
       ! 
       ! forms the dot product of two vectors.
       ! uses unrolled loops for increments equal to one.
       ! jack dongarra, linpack, 3/11/78.
       ! modified 12/3/93, array(1) declarations changed to array(*)
       ! 
       DOUBLE PRECISION dx(*),dy(*),dtemp
       INTEGER :: i,incx,incy,ix,iy,m,mp1,n
       ! 
       ddot = 0.0_real_8
       dtemp = 0.0_real_8
#ifdef __SAFE
       IF (n.LE.0)RETURN
#endif
       IF (incx.EQ.1.AND.incy.EQ.1)go to 20
       ! 
       ! code for unequal increments or equal increments
       ! not equal to 1
       ! 
       ix = 1
       iy = 1
       IF (incx.LT.0)ix = (-n+1)*incx + 1
       IF (incy.LT.0)iy = (-n+1)*incy + 1
       DO 10 i = 1,n
          dtemp = dtemp + dx(ix)*dy(iy)
          ix = ix + incx
          iy = iy + incy
10        CONTINUE
          ddot = dtemp
          RETURN
          ! 
          ! code for both increments equal to 1
          ! 
          ! 
          ! clean-up loop
          ! 
20        m = MOD(n,5)
          IF ( m .EQ. 0 ) go to 40
          DO 30 i = 1,m
             dtemp = dtemp + dx(i)*dy(i)
30           CONTINUE
             IF ( n .LT. 5 ) go to 60
40           mp1 = m + 1
             DO 50 i = mp1,n,5
                dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +&
                     dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
50              CONTINUE
60              ddot = dtemp
                RETURN
             ENDDO

             ! ----------------------------------------------------------------------
             SUBROUTINE  drot (n,dx,incx,dy,incy,c,s)
               ! ----------------------------------------------------------------------
               ! 
               ! applies a plane rotation.
               ! jack dongarra, linpack, 3/11/78.
               ! modified 12/3/93, array(1) declarations changed to array(*)
               ! 
               USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
               USE error_handling, ONLY: stopgm
               USE timer, ONLY: tiset, tihalt
               DOUBLE PRECISION dx(*),dy(*),dtemp,c,s
               INTEGER :: i,incx,incy,ix,iy,n
               ! 
#ifdef __SAFE
               IF (n.LE.0)RETURN
#endif
               IF (incx.EQ.1.AND.incy.EQ.1)go to 20
               ! 
               ! code for unequal increments or equal increments not equal
               ! to 1
               ! 
               ix = 1
               iy = 1
               IF (incx.LT.0)ix = (-n+1)*incx + 1
               IF (incy.LT.0)iy = (-n+1)*incy + 1
               DO 10 i = 1,n
                  dtemp = c*dx(ix) + s*dy(iy)
                  dy(iy) = c*dy(iy) - s*dx(ix)
                  dx(ix) = dtemp
                  ix = ix + incx
                  iy = iy + incy
10                CONTINUE
                  RETURN
                  ! 
                  ! code for both increments equal to 1
                  ! 
20                DO 30 i = 1,n
                     dtemp = c*dx(i) + s*dy(i)
                     dy(i) = c*dy(i) - s*dx(i)
                     dx(i) = dtemp
30                   CONTINUE
                     RETURN
                  ENDDO


