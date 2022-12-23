MODULE symmetry_utils
  USE error_handling,                  ONLY: stopgm
  USE parac,                           ONLY: paral

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: pgsym
  PUBLIC :: fccsym
  PUBLIC :: bccsym
  PUBLIC :: trgsym
  PUBLIC :: bctsym
  PUBLIC :: alfasym

CONTAINS

  ! ==================================================================
  SUBROUTINE pgsym(indpg,s,ntot)
    ! ==--------------------------------------------------------------==
    ! sets up the point-group rotation matrices for all 32 point groups
    ! for the simple (p) bravais lattices.
    ! These matrices all have integer elements, and express 
    ! the relations between the primitive basis vectors 
    ! due to each point-group rotation.  
    ! The special axis is always the z-axis, one basal-plane
    ! vector is always along x, and the other basal-plane vector is at
    ! angle beta for monoclinic (beta is not actually used), at 120
    ! degrees for trigonal and hexagonal(p) groups, and at 90 degrees
    ! for remaining groups.
    ! subroutine trgsym takes care of the trigonal(r) groups, subroutine
    ! fccsym takes care of the cubic (f) groups, subroutine bccsym takes
    ! care of the cubic (i) groups, and subroutine bctsym takes care of
    ! the tetragonal (i) groups.  for these latter
    ! groups, the crystallographic vectors are chosen
    ! differently.  for the triclinic groups, the matrices are
    ! simply identity and inversion.
    ! ==--------------------------------------------------------------==
    ! indpg is the point-group index, defined as follows:
    ! ==--------------------------------------------------------------==
    ! indpg   group    indpg   group     indpg   group     indpg   group
    ! 1     1 (c1)      9     3m (c3v)    17   4/mmm(d4h)   25   222(d2)
    ! 2    <1>(ci)     10     <3>m(d3d)   18   6 (c6)       26  mm2(c2v)
    ! 3     2 (c2)     11     4 (c4)      19   <6>(c3h)     27  mmm(d2h)
    ! 4     m (c1h)    12     <4>(s4)     20   6/m(c6h)     28  23 (t)
    ! 5     2/m(c2h)   13     4/m(c4h)    21   622(d6)      29  m3 (th)
    ! 6     3 (c3)     14     422(d4)     22   6mm(c6v)     30  432 (o)
    ! 7     <3>(c3i)   15     4mm(c4v)    23   <6>m2(d3h)   31 <4>3m(td)
    ! 8     32 (d3)    16    <4>2m(d2d)   24    6/mmm(d6h)  32  m3m(oh)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: indpg, s(3,3,48), ntot

    INTEGER                                  :: i, j, jind, k, kind, sum

    IF (indpg.LT.1) go to 100
    DO k=1,48
       DO j=1,3
          DO i=1,3
             s(i,j,k)=0
          ENDDO
       ENDDO
    ENDDO
    DO i=1,3
       s(i,i,1)=1
    ENDDO
    ! ==--------------------------------------------------------------==
    ! == Triclinic groups.                                            ==
    ! ==--------------------------------------------------------------==
    IF (indpg.EQ.1) THEN
       ! point group 1. [1]
       ntot=1
       RETURN
    ELSEIF (indpg.EQ.2) THEN
       ! add inversion for point group <1>. [2]
       ntot=2
       DO i=1,3
          s(i,i,2)=-1
       ENDDO
       RETURN
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == Monoclinic groups.                                           ==
    ! ==--------------------------------------------------------------==
    IF (indpg.LE.5) THEN
       ntot=2
       DO i=1,3
          s(i,i,2)=1
          s(i,i,3)=1
          s(i,i,4)=1
       ENDDO
       IF (indpg.EQ.3) THEN
          ! point group 2. [3]
          s(1,1,2)=-1
          s(2,2,2)=-1
          RETURN
       ELSEIF (indpg.EQ.4) THEN
          ! point group m. [4]
          s(3,3,2)=-1
          RETURN
       ELSEIF (indpg.EQ.5) THEN
          ntot=4
          ! point group 2/m. [5]
          s(1,1,2)=-1
          s(2,2,2)=-1
          s(3,3,3)=-1
          DO i=1,3
             s(i,i,4)=-1
          ENDDO
          RETURN
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == Trigonal groups.                                             ==
    ! ==--------------------------------------------------------------==
    IF (indpg.LE.10) THEN
       ! point group 3. [6]
       ntot=3
       s(1,2,2)=1
       s(2,1,2)=-1
       s(2,2,2)=-1
       s(3,3,2)=1
       s(1,1,3)=-1
       s(1,2,3)=-1
       s(2,1,3)=1
       s(3,3,3)=1
       IF (indpg.EQ.6) THEN
          RETURN
       ELSEIF (indpg.EQ.8.OR.indpg.EQ.9) THEN
          ntot=6
          ! set up c2x. 
          s(1,1,4)=1
          s(3,3,4)=-1
          s(2,1,4)=-1
          s(2,2,4)=-1
          IF (indpg.EQ.9) THEN
             ! set up mx.
             DO j=1,3
                DO i=1,3
                   s(i,j,4)=-s(i,j,4)
                ENDDO
             ENDDO
             ! Unnecessary (Thierry Deutsch 13/1/98)
             s(3,3,4)=1
          ENDIF
          ! point group 32 =3(x)c2x [8] or 3m = 3(x)mx [9].
          DO kind=2,3
             DO i=1,3
                DO j=1,3
                   sum=0
                   DO k=1,3
                      sum=sum+s(i,k,4)*s(k,j,kind)
                   ENDDO
                   s(i,j,kind+3)=sum
                ENDDO
             ENDDO
          ENDDO
          RETURN
       ELSEIF (indpg.EQ.7.OR.indpg.EQ.10) THEN
          ntot=6
          ! point group <3> = 3(x)i [7]
          DO kind=1,3
             DO j=1,3
                DO i=1,3
                   s(i,j,kind+3)=-s(i,j,kind)
                ENDDO
             ENDDO
          ENDDO
          IF (indpg.EQ.7) THEN
             RETURN
          ELSE
             ! point group <3>m = <3>(x)mx [10]
             ntot=12
             s(1,1,7)=-1
             s(2,1,7)=1
             s(2,2,7)=1
             s(3,3,7)=1
             DO kind=2,6
                DO j=1,3
                   DO i=1,3
                      sum=0
                      DO k=1,3
                         sum=sum+s(i,k,7)*s(k,j,kind)
                      ENDDO
                      s(i,j,kind+6)=sum
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDIF
       RETURN
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == Tetragonal groups.                                           ==
    ! ==--------------------------------------------------------------==
    IF (indpg.LE.17) THEN
       ! point group 4. [11]
       ntot=4
       s(1,2,2)=1
       s(2,1,2)=-1
       s(3,3,2)=1
       s(1,1,3)=-1
       s(2,2,3)=-1
       s(3,3,3)=1
       s(1,2,4)=-1
       s(2,1,4)=1
       s(3,3,4)=1
       IF (indpg.EQ.11) THEN
          RETURN
       ELSEIF (indpg.EQ.12.OR.indpg.EQ.16) THEN
          ! point group <4>. [12]
          DO kind=2,4,2
             DO j=1,3
                DO i=1,3
                   s(i,j,kind)=-s(i,j,kind)
                ENDDO
             ENDDO
          ENDDO
          IF (indpg.EQ.12) THEN
             RETURN
          ELSEIF (indpg.EQ.16) THEN
             ! point group <4>2m = <4>(x)c2x. [16]
             ntot=8
             ! the following two cards are wrong (original janak s program)
             ! s(1,1,5)= 1
             ! s(2,2,5)=-1
             ! the following two are the correct cards 
             ! (stefano de gironcoli dixit ad 1989, 26 may)
             s(1,2,5)= 1
             s(2,1,5)= 1
             s(3,3,5)=-1
             DO kind=2,4
                DO j=1,3
                   DO i=1,3
                      sum=0
                      DO k=1,3
                         sum=sum+s(i,k,5)*s(k,j,kind)
                      ENDDO
                      s(i,j,kind+4)=sum
                   ENDDO
                ENDDO
             ENDDO
             RETURN
          ENDIF
       ELSEIF(indpg.EQ.13.OR.indpg.EQ.14.OR.&
            indpg.EQ.15.OR.indpg.EQ.17) THEN
          ntot=8
          IF (indpg.EQ.13) THEN
             ! point group 4/m = mz(x)4 [13]
             s(1,1,5)=1
             s(2,2,5)=1
             s(3,3,5)=-1
          ELSEIF (indpg.EQ.15) THEN
             ! point group 4mm = mx(x)4 [15]
             s(1,1,5)=-1
             s(2,2,5)=1
             s(3,3,5)=1
          ELSEIF (indpg.EQ.14.OR.indpg.EQ.17) THEN
             ! point groups 422 = c2x(x)4 or 4/mmm = i(x)422. [14]
             s(1,1,5)=1
             s(2,2,5)=-1
             s(3,3,5)=-1
          ENDIF
          ! Add 3 rotations products of 2, 3, 4 with 5.
          DO kind=2,4
             DO j=1,3
                DO i=1,3
                   sum=0
                   DO k=1,3
                      sum=sum+s(i,k,5)*s(k,j,kind)
                   ENDDO
                   s(i,j,kind+4)=sum
                ENDDO
             ENDDO
          ENDDO
          IF (indpg.EQ.14) THEN
             RETURN
          ELSEIF (indpg.EQ.17) THEN
             ! Add inversion and products. [17]
             ntot=16
             DO kind=1,8
                DO j=1,3
                   DO i=1,3
                      s(i,j,kind+8)=-s(i,j,kind)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDIF
       RETURN
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == Hexagonal point groups.                                      ==
    ! ==--------------------------------------------------------------==
    IF (indpg.LE.24) THEN
       ntot=6
       ! point group 6. [18]
       s(1,1,2)=1
       s(1,2,2)=1
       s(2,1,2)=-1
       s(1,2,3)=1
       s(2,1,3)=-1
       s(2,2,3)=-1
       s(1,1,4)=-1
       s(2,2,4)=-1
       s(1,1,5)=-1
       s(1,2,5)=-1
       s(2,1,5)=1
       s(1,2,6)=-1
       s(2,1,6)=1
       s(2,2,6)=1
       DO i=2,6
          s(3,3,i)=1
       ENDDO
       IF (indpg.EQ.18) THEN
          RETURN
       ELSEIF (indpg.EQ.19.OR.indpg.EQ.23) THEN
          ! point group <6>. [19]
          DO kind=2,6,2
             DO j=1,3
                DO i=1,3
                   s(i,j,kind)=-s(i,j,kind)
                ENDDO
             ENDDO
          ENDDO
          IF (indpg.EQ.19) RETURN
          ! point group <6>m2 = c2y(x)<6>. [23]
          ntot=12
          s(1,1,7)=-1
          ! There was a missing card (Thierry Deutsch 13/1/98):
          s(1,2,7)=-1
          s(2,2,7)=1
          s(3,3,7)=-1
          DO kind=2,6
             DO j=1,3
                DO i=1,3
                   sum=0
                   DO k=1,3
                      sum=sum+s(i,k,7)*s(k,j,kind)
                   ENDDO
                   s(i,j,kind+6)=sum
                ENDDO
             ENDDO
          ENDDO
          RETURN
       ENDIF
       ntot=12
       IF (indpg.EQ.20) THEN
          ! point group 6/m = i(x)6 [20]
          DO kind=1,6
             DO j=1,3
                DO i=1,3
                   s(i,j,kind+6)=-s(i,j,kind)
                ENDDO
             ENDDO
          ENDDO
          RETURN
       ELSEIF (indpg.EQ.22) THEN
          ! point group 6mm = mx(x)6 [22]
          ! the following two cards are wrong (original janak s program)
          ! s(1,1,7)=-1
          ! s(2,2,7)=1
          ! the following two are the correct cards
          ! (lausanne, 10 june 1988)
          s(1,2,7)=1
          s(2,1,7)=1
          s(3,3,7)=1
       ELSE
          ! point groups 622 = c2x(x)6 [21] or 6/mmm = i(x)622. [24]  
          ! the following two cards are wrong (original janak s program)
          ! 48 s(1,1,7)=1
          ! s(2,2,7)=-1
          ! the following two are the correct cards as sgamed by s. massidda
          s(1,2,7)=1
          s(2,1,7)=1
          s(3,3,7)=-1
       ENDIF
       ! Add products 2, 3, 4, 5, 6 with 7
       DO kind=2,6
          DO j=1,3
             DO i=1,3
                sum=0
                DO k=1,3
                   sum=sum+s(i,k,7)*s(k,j,kind)
                ENDDO
                s(i,j,kind+6)=sum
             ENDDO
          ENDDO
       ENDDO
       IF (indpg.EQ.24) THEN
          RETURN
       ELSE                  ! [21]
          ! Add inversion and products
          ntot=24
          DO kind=1,12
             DO j=1,3
                DO i=1,3
                   s(i,j,kind+12)=-s(i,j,kind)
                ENDDO
             ENDDO
          ENDDO
       ENDIF
       RETURN
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == Orthorhombic point groups.                                   ==
    ! ==--------------------------------------------------------------==
    IF (indpg.LE.27) THEN
       ntot=4
       IF (indpg.EQ.26) THEN
          ! point group mm2. [26]
          DO kind=2,4
             DO i=1,3
                s(i,i,kind)=1
             ENDDO
          ENDDO
          s(1,1,2)=-1
          s(2,2,3)=-1
          s(1,1,4)=-1
          s(2,2,4)=-1
          RETURN
       ENDIF
       ! point group 222. [25]
       DO kind=2,4
          DO i=1,3
             s(i,i,kind)=-1
          ENDDO
       ENDDO
       s(1,1,2)=1
       s(2,2,3)=1
       s(3,3,4)=1
       IF (indpg.EQ.25) THEN
          RETURN
       ELSEIF (indpg.EQ.27) THEN
          ! point group mmm = i(x)222 [27]
          ntot=8
          DO kind=1,4
             DO j=1,3
                DO i=1,3
                   s(i,j,kind+4)=-s(i,j,kind)
                ENDDO
             ENDDO
          ENDDO
       ENDIF
       RETURN
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == Cubic point groups.                                          ==
    ! ==--------------------------------------------------------------==
    IF (indpg.LE.32) THEN
       ntot=12
       ! point group 23. [28]
       DO kind=2,4
          DO i=1,3
             s(i,i,kind)=-1
          ENDDO
       ENDDO
       s(1,1,2)=1
       s(2,2,3)=1
       s(3,3,4)=1
       s(1,2,5)=1
       s(2,3,5)=1
       s(3,1,5)=1
       s(1,3,6)=1
       s(2,1,6)=1
       s(3,2,6)=1
       DO kind=5,6
          DO jind=2,4
             DO j=1,3
                DO i=1,3
                   sum=0
                   DO k=1,3
                      sum=sum+s(i,k,jind)*s(k,j,kind)
                   ENDDO
                   s(i,j,kind+2*jind-2)=sum
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       IF (indpg.EQ.28) RETURN
       ntot=24
       IF (indpg.EQ.29) THEN
          ! point group m3 = i(x)23 [29]
          DO kind=1,12
             DO j=1,3
                DO i=1,3
                   s(i,j,kind+12)=-s(i,j,kind)
                ENDDO
             ENDDO
          ENDDO
          RETURN
       ENDIF
       IF (indpg.EQ.30.OR.indpg.EQ.32) THEN
          ! point group 432 = c4x(x)23 [30] [32]
          s(1,1,13)=1
          s(2,3,13)=1
          s(3,2,13)=-1
          DO kind=2,12
             DO j=1,3
                DO i=1,3
                   sum=0
                   DO k=1,3
                      sum=sum+s(i,k,13)*s(k,j,kind)
                   ENDDO
                   s(i,j,kind+12)=sum
                ENDDO
             ENDDO
          ENDDO
          IF (indpg.EQ.30) RETURN
          ! point group m3m = i(x)432 [32]
          ntot=48
          DO kind=1,24
             DO j=1,3
                DO i=1,3
                   s(i,j,kind+24)=-s(i,j,kind)
                ENDDO
             ENDDO
          ENDDO
          RETURN
       ENDIF
       ! point <4>3m = <c4x>(x)23 [31]
       s(1,1,13)=-1
       s(2,3,13)=-1
       s(3,2,13)=1
       DO kind=2,12
          DO j=1,3
             DO i=1,3
                sum=0
                DO k=1,3
                   sum=sum+s(i,k,13)*s(k,j,kind)
                ENDDO
                s(i,j,kind+12)=sum
             ENDDO
          ENDDO
       ENDDO
       RETURN
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == If indpg exceeds 32 or is less than 1.                       ==
    ! ==--------------------------------------------------------------==
100 CONTINUE
    IF (paral%io_parent)&
         WRITE(6,110) indpg
110 FORMAT(' point-group index',i4,' is not between 1 and 32.',&
         ' stopping')
    CALL stopgm('pgsym',' ',& 
         __LINE__,__FILE__)
  END SUBROUTINE pgsym
  ! ==================================================================
  SUBROUTINE fccsym(indpg,s,ntot)
    ! ==--------------------------------------------------------------==
    ! point group rotation matrices for the fcc bravais lattice.  it is
    ! assumed that a1=(a/2)(-1,0,1), a2=(a/2)(0,1,1), a3=(a/2)(-1,1,0).
    ! the rotation matrices express how these basis vectors are trans-
    ! formed into one another by point-group rotations.
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: indpg, s(3,3,48), ntot

    INTEGER                                  :: i, j, jind, k, kind, sum

    IF (indpg.LT.28.OR.indpg.GT.32) go to 100
    DO k=1,48
       DO j=1,3
          DO i=1,3
             s(i,j,k)=0
          ENDDO
       ENDDO
    ENDDO
    DO i=1,3
       s(i,i,1)=1
    ENDDO
    ntot=12
    ! point group 23.
    s(1,2,2)=-1
    s(1,3,2)=1
    s(2,2,2)=-1
    s(3,1,2)=1
    s(3,2,2)=-1
    s(1,1,3)=-1
    s(2,1,3)=-1
    s(2,3,3)=1
    s(3,1,3)=-1
    s(3,2,3)=1
    s(1,2,4)=1
    s(1,3,4)=-1
    s(2,1,4)=1
    s(2,3,4)=-1
    s(3,3,4)=-1
    s(1,3,5)=-1
    s(2,2,5)=1
    s(2,3,5)=-1
    s(3,1,5)=1
    s(3,3,5)=-1
    s(1,1,6)=-1
    s(1,3,6)=1
    s(2,1,6)=-1
    s(2,2,6)=1
    s(3,1,6)=-1
    DO kind=5,6
       DO jind=2,4
          DO j=1,3
             DO i=1,3
                sum=0
                DO k=1,3
                   sum=sum+s(i,k,jind)*s(k,j,kind)
                ENDDO
                s(i,j,kind+2*jind-2)=sum
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    IF (indpg.EQ.28) RETURN
    ntot=24
    IF (indpg.NE.29) go to 76
    ! point group m3 = i(x)23
    DO kind=1,12
       DO j=1,3
          DO i=1,3
             s(i,j,kind+12)=-s(i,j,kind)
          ENDDO
       ENDDO
    ENDDO
    RETURN
76  CONTINUE
    IF (indpg.EQ.31) go to 80
    ! point group 432 = c4x(x)23
    s(1,3,13)=1
    s(2,1,13)=-1
    s(2,3,13)=1
    s(3,2,13)=-1
    s(3,3,13)=1
    DO kind=2,12
       DO j=1,3
          DO i=1,3
             sum=0
             DO k=1,3
                sum=sum+s(i,k,13)*s(k,j,kind)
             ENDDO
             s(i,j,kind+12)=sum
          ENDDO
       ENDDO
    ENDDO
    IF (indpg.EQ.30) RETURN
    ! point group m3m = i(x)432
    ntot=48
    DO kind=1,24
       DO j=1,3
          DO i=1,3
             s(i,j,kind+24)=-s(i,j,kind)
          ENDDO
       ENDDO
    ENDDO
    RETURN
    ! point <4>3m = <c4x>(x)23
80  CONTINUE
    s(1,3,13)=-1
    s(2,1,13)=1
    s(2,3,13)=-1
    s(3,2,13)=1
    s(3,3,13)=-1
    DO kind=2,12
       DO j=1,3
          DO i=1,3
             sum=0
             DO k=1,3
                sum=sum+s(i,k,13)*s(k,j,kind)
             ENDDO
             s(i,j,kind+12)=sum
          ENDDO
       ENDDO
    ENDDO
    RETURN
    ! here if indpg exceeds 32 or is less than 28.
100 CONTINUE
    IF (paral%io_parent)&
         WRITE(6,110) indpg
110 FORMAT(' point-group index',i4,&
         ' is not between 28 and 32 for a cubic group. stopping')
    CALL stopgm('fccsym',' ',& 
         __LINE__,__FILE__)
  END SUBROUTINE fccsym
  ! ==================================================================
  SUBROUTINE bccsym(indpg,s,ntot)
    ! ==--------------------------------------------------------------==
    ! point group rotation matrices for cubic (bcc) bravais lattice.
    ! it is assumed that a1=(a/2)(1,1,1), a2=(a/2)(-1,1,1),
    ! a3=(a/2)(-1,-1,1).  the rotation matrices express how these basis
    ! vectors are transformed into one another by point-group rotations.
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: indpg, s(3,3,48), ntot

    INTEGER                                  :: i, j, jind, k, kind, sum

    DO k=1,48
       DO j=1,3
          DO i=1,3
             s(i,j,k)=0
          ENDDO
       ENDDO
    ENDDO
    DO i=1,3
       s(i,i,1)=1
    ENDDO
    IF (indpg.LT.28.OR.indpg.GT.32) go to 100
    ntot=12
    ! point group 23.
    s(1,2,2)=-1
    s(2,1,2)=-1
    s(3,1,2)=-1
    s(3,2,2)=1
    s(3,3,2)=-1
    s(1,1,3)=-1
    s(1,2,3)=1
    s(1,3,3)=-1
    s(2,3,3)=-1
    s(3,2,3)=-1
    s(1,3,4)=1
    s(2,1,4)=1
    s(2,2,4)=-1
    s(2,3,4)=1
    s(3,1,4)=1
    s(1,1,5)=1
    s(2,1,5)=1
    s(2,2,5)=-1
    s(2,3,5)=1
    s(3,2,5)=-1
    s(1,1,6)=1
    s(2,3,6)=-1
    s(3,1,6)=-1
    s(3,2,6)=1
    s(3,3,6)=-1
    DO kind=5,6
       DO jind=2,4
          DO j=1,3
             DO i=1,3
                sum=0
                DO k=1,3
                   sum=sum+s(i,k,jind)*s(k,j,kind)
                ENDDO
                s(i,j,kind+2*jind-2)=sum
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    IF (indpg.EQ.28) RETURN
    ntot=24
    IF (indpg.NE.29) go to 76
    ! point group m3 = i(x)23
    DO kind=1,12
       DO j=1,3
          DO i=1,3
             s(i,j,kind+12)=-s(i,j,kind)
          ENDDO
       ENDDO
    ENDDO
    RETURN
76  CONTINUE
    IF (indpg.EQ.31) go to 80
    ! point group 432 = c4x(x)23
    s(1,3,13)=-1
    s(2,1,13)=-1
    s(2,2,13)=1
    s(2,3,13)=-1
    s(3,2,13)=1
    DO kind=2,12
       DO j=1,3
          DO i=1,3
             sum=0
             DO k=1,3
                sum=sum+s(i,k,13)*s(k,j,kind)
             ENDDO
             s(i,j,kind+12)=sum
          ENDDO
       ENDDO
    ENDDO
    IF (indpg.EQ.30) RETURN
    ! point group m3m = i(x)432
    ntot=48
    DO kind=1,24
       DO j=1,3
          DO i=1,3
             s(i,j,kind+24)=-s(i,j,kind)
          ENDDO
       ENDDO
    ENDDO
    RETURN
    ! point <4>3m = <c4x>(x)23
80  CONTINUE
    s(1,3,13)=1
    s(2,1,13)=1
    s(2,2,13)=-1
    s(2,3,13)=1
    s(3,2,13)=-1
    DO kind=2,12
       DO j=1,3
          DO i=1,3
             sum=0
             DO k=1,3
                sum=sum+s(i,k,13)*s(k,j,kind)
             ENDDO
             s(i,j,kind+12)=sum
          ENDDO
       ENDDO
    ENDDO
    RETURN
    ! here if indpg exceeds 32 or is less than 28.
100 CONTINUE
    IF (paral%io_parent)&
         WRITE(6,110) indpg
110 FORMAT(' point-group index',i4,&
         ' is not between 28 and 32 for a cubic group. stopping')
    CALL stopgm('bccsym',' ',& 
         __LINE__,__FILE__)
  END SUBROUTINE bccsym
  ! ==================================================================
  SUBROUTINE trgsym(indpg,s,ntot)
    ! ==--------------------------------------------------------------==
    ! sets up the rotation matrices for trigonal(r) groups.  for these
    ! groups, the z-axis is chosen as the 3-fold axis, but the
    ! crystallographic vectors form a three-fold star around the z-axis,
    ! and the primitive cell is a simple rhombohedron.
    ! if c is the cosine of the angle between any pair of crystallo-
    ! graphic vectors, and if tx=sqrt((1-c)/2), ty=sqrt((1-c)/6),
    ! tz=sqrt((1+2c)/3), the crystallographic vectors are
    ! a1=a(0,2ty,tz),  a2=a(tx,-ty,tz),  a3=a(-tx,-ty,tz).
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: indpg, s(3,3,48), ntot

    INTEGER                                  :: i, j, k, kind, sum

    IF (indpg.LT.6.OR.indpg.GT.10) go to 100
    DO k=1,48
       DO j=1,3
          DO i=1,3
             s(i,j,k)=0
          ENDDO
       ENDDO
    ENDDO
    DO i=1,3
       s(i,i,1)=1
    ENDDO
    ! point group 3.
    ntot=3
    s(1,3,2)=1
    s(2,1,2)=1
    s(3,2,2)=1
    s(1,2,3)=1
    s(2,3,3)=1
    s(3,1,3)=1
    IF (indpg.EQ.6) RETURN
    IF (indpg.EQ.7.OR.indpg.EQ.10) go to 15
    ntot=6
    ! set up c2x.
    s(1,1,4)=-1
    s(3,2,4)=-1
    s(2,3,4)=-1
    IF (indpg.EQ.8) go to 12
    ! set up mx.
    DO j=1,3
       DO i=1,3
          s(i,j,4)=-s(i,j,4)
       ENDDO
    ENDDO
    ! point group 32 =3(x)c2x or 3m = 3(x)mx.
12  CONTINUE
    DO kind=2,3
       DO i=1,3
          DO j=1,3
             sum=0
             DO k=1,3
                sum=sum+s(i,k,4)*s(k,j,kind)
             ENDDO
             s(i,j,kind+3)=sum
          ENDDO
       ENDDO
    ENDDO
    RETURN
15  CONTINUE
    ntot=6
    ! point group <3> = 3(x)i
    DO kind=1,3
       DO j=1,3
          DO i=1,3
             s(i,j,kind+3)=-s(i,j,kind)
          ENDDO
       ENDDO
    ENDDO
    IF (indpg.EQ.7) RETURN
    ! point group <3>m = <3>(x)mx
    ntot=12
    s(1,1,7)=1
    s(2,3,7)=1
    s(3,2,7)=1
    DO kind=2,6
       DO j=1,3
          DO i=1,3
             sum=0
             DO k=1,3
                sum=sum+s(i,k,7)*s(k,j,kind)
             ENDDO
             s(i,j,kind+6)=sum
          ENDDO
       ENDDO
    ENDDO
    RETURN
100 CONTINUE
    IF (paral%io_parent)&
         WRITE(6,110) indpg
110 FORMAT(' point group',i3,' is not trigonal. stopping')
    CALL stopgm('trgsym',' ',& 
         __LINE__,__FILE__)
  END SUBROUTINE trgsym
  ! ==================================================================
  SUBROUTINE bctsym(indpg,s,ntot)
    ! ==--------------------------------------------------------------==
    ! point group rotation matrices for tetragonal (i) bravais lattices.
    ! it is assumed that a1=(a/2,a/2,c/2), a2=(a/2,-a/2,c/2), and
    ! a3=(-a/2,-a/2,c/2).  the rotation matrices express how these
    ! basis vectors are transformed into one another by point group
    ! rotations.
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: indpg, s(3,3,48), ntot

    INTEGER                                  :: i, j, k, kind, sum

    DO k=1,48
       DO j=1,3
          DO i=1,3
             s(i,j,k)=0
          ENDDO
       ENDDO
    ENDDO
    DO i=1,3
       s(i,i,1)=1
    ENDDO
    IF (indpg.LT.11.OR.indpg.GT.17) go to 40
    ! point group 4.
    ntot=4
    s(1,2,2)=1
    s(2,3,2)=1
    s(3,1,2)=1
    s(3,2,2)=-1
    s(3,3,2)=1
    s(1,3,3)=1
    s(2,1,3)=1
    s(2,2,3)=-1
    s(2,3,3)=1
    s(3,1,3)=1
    s(1,1,4)=1
    s(1,2,4)=-1
    s(1,3,4)=1
    s(2,1,4)=1
    s(3,2,4)=1
    IF (indpg.EQ.11) RETURN
    IF (indpg.NE.12.AND.indpg.NE.16) go to 25
    ! point group <4>.
    DO kind=2,4,2
       DO j=1,3
          DO i=1,3
             s(i,j,kind)=-s(i,j,kind)
          ENDDO
       ENDDO
    ENDDO
    IF (indpg.EQ.12) RETURN
    ! point group <4>2m = <4>(x)c2x.
    ntot=8
    s(1,1,5)=-1
    s(1,2,5)=1
    s(1,3,5)=-1
    s(2,3,5)=-1
    s(3,2,5)=-1
    DO kind=2,4
       DO j=1,3
          DO i=1,3
             sum=0
             DO k=1,3
                sum=sum+s(i,k,5)*s(k,j,kind)
             ENDDO
             s(i,j,kind+4)=sum
          ENDDO
       ENDDO
    ENDDO
    RETURN
25  CONTINUE
    ntot=8
    IF (indpg.NE.13) go to 26
    ! point group 4/m = mz(x)4
    s(1,3,5)=-1
    s(2,1,5)=-1
    s(2,2,5)=1
    s(2,3,5)=-1
    s(3,1,5)=-1
    go to 30
26  CONTINUE
    IF (indpg.NE.15) go to 27
    ! point group 4mm = mx(x)4.
    s(1,1,5)=1
    s(1,2,5)=-1
    s(1,3,5)=1
    s(2,3,5)=1
    s(3,2,5)=1
    go to 30
    ! point groups 422 = c2x(x)4 and 4/mmm = i(x)422.
27  CONTINUE
    s(1,1,5)=-1
    s(1,2,5)=1
    s(1,3,5)=-1
    s(2,3,5)=-1
    s(3,2,5)=-1
30  CONTINUE
    DO kind=2,4
       DO j=1,3
          DO i=1,3
             sum=0
             DO k=1,3
                sum=sum+s(i,k,5)*s(k,j,kind)
             ENDDO
             s(i,j,kind+4)=sum
          ENDDO
       ENDDO
    ENDDO
    IF (indpg.NE.17) RETURN
    ntot=16
    DO kind=1,8
       DO j=1,3
          DO i=1,3
             s(i,j,kind+8)=-s(i,j,kind)
          ENDDO
       ENDDO
    ENDDO
    RETURN
40  CONTINUE
    IF (paral%io_parent)&
         WRITE(6,45) indpg
45  FORMAT(' point-group index',i3,&
         ' is not between 11 and 17 for a tetragonal group. stopping')
    CALL stopgm('bctsym',' ',& 
         __LINE__,__FILE__)
  END SUBROUTINE bctsym
  ! ==================================================================
  SUBROUTINE alfasym(indpg,s,ntot)
    ! ==--------------------------------------------------------------==
    ! SUBROUTINE PER LE SIMMETRIE DI ALFA-GA.
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: indpg, s(3,3,48), ntot

    INTEGER                                  :: i, j, k

    IF ((indpg.NE.1).AND.(indpg.NE.25)) go to 100
    DO k=1,8
       DO i=1,3
          DO j=1,3
             s(i,j,k)=0
          ENDDO
       ENDDO
    ENDDO
    DO i=1,3
       s(i,i,1)=1
    ENDDO
    IF (indpg.EQ.25) go to 20
    ntot=1
    RETURN
20  CONTINUE
    ntot=4
    ! R2
    s(1,1,2)=-1
    s(2,3,2)= 1
    s(3,2,2)= 1
    ! R3
    s(1,1,3)=-1
    s(2,3,3)=-1
    s(3,2,3)=-1
    ! R4
    s(1,1,4)=1
    s(2,2,4)=-1
    s(3,3,4)=-1
    RETURN
100 CONTINUE
    IF (paral%io_parent)&
         WRITE(6,110)indpg
110 FORMAT(' POINT-GROUP INDEX ',i4,&
         ' DOES NOT CORRESPOND TO ALFA-GA, STOPPING')
    CALL stopgm('ALFASYM',' ',& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE alfasym
  ! ==================================================================

END MODULE symmetry_utils
