! ==================================================================
SUBROUTINE dgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
  ! ==--------------------------------------------------------------==
  ! == Purpose                                                      ==
  ! == =======                                                      ==
  ! ==                                                              ==
  ! == DGEMM performs one of the matrix-matrix operations           ==
  ! ==                                                              ==
  ! == C := alpha*op(A)*op(B)+beta*C,                               ==
  ! ==                                                              ==
  ! ==   where op(X) is one of                                      ==
  ! ==                                                              ==
  ! ==   op(X) = X or op(X) = X',                                   ==
  ! ==                                                              ==
  ! == alpha and beta are scalars, and A, B and C are matrices,     ==
  ! == with op(A) an m by k matrix, op(B) a k by n matrix           ==
  ! == and C an m by n matrix.                                      ==
  ! ==                                                              ==
  ! == Parameters                                                   ==
  ! == ==========                                                   ==
  ! ==                                                              ==
  ! == TRANSA - CHARACTER*1.                                        ==
  ! ==   On entry, TRANSA specifies the form of op(A) to be used    ==
  ! ==   in the matrix multiplication as follows:                   ==
  ! ==                                                              ==
  ! ==     TRANSA = 'N' or 'n', op(A) = A.                          ==
  ! ==                                                              ==
  ! ==     TRANSA = 'T' or 't', op(A) = A'.                         ==
  ! ==                                                              ==
  ! ==     TRANSA = 'C' or 'c', op(A) = A'.                         ==
  ! ==                                                              ==
  ! ==   Unchanged on exit.                                         ==
  ! ==                                                              ==
  ! == TRANSB - CHARACTER*1.                                        ==
  ! ==   On entry, TRANSB specifies the form of op(B) to be used    ==
  ! ==   in the matrix multiplication as follows:                   ==
  ! ==                                                              ==
  ! ==     TRANSB = 'N' or 'n', op(B) = B.                          ==
  ! ==                                                              ==
  ! ==     TRANSB = 'T' or 't', op(B) = B'.                         ==
  ! ==                                                              ==
  ! ==     TRANSB = 'C' or 'c', op(B) = B'.                         ==
  ! ==                                                              ==
  ! ==   Unchanged on exit.                                         ==
  ! ==                                                              ==
  ! == M - INTEGER.                                                 ==
  ! ==   On entry, M specifies the number of rows of the matrix     ==
  ! ==   op(A) and of the matrix C. M must be at least zero.        ==
  ! ==   Unchanged on exit.                                         ==
  ! ==                                                              ==
  ! == N - INTEGER.                                                 ==
  ! ==   On entry, N specifies the number of columns of the matrix  ==
  ! ==   op(B) and the number of columns of the matrix C. N must    ==
  ! ==   be at least zero.                                          ==
  ! ==   Unchanged on exit.                                         ==
  ! ==                                                              ==
  ! == K - INTEGER.                                                 ==
  ! ==   On entry, K specifies the number of columns of the matrix  ==
  ! ==   op(A) and the number of rows of the matrix op(B). K must   ==
  ! ==   be at least  zero.                                         ==
  ! ==   Unchanged on exit.                                         ==
  ! ==                                                              ==
  ! == ALPHA - DOUBLE PRECISION.                                    ==
  ! ==   On entry, ALPHA specifies the scalar alpha.                ==
  ! ==   Unchanged on exit.                                         ==
  ! ==                                                              ==
  ! == A - DOUBLE PRECISION array of DIMENSION (LDA,ka), where ka   ==
  ! ==   is k when TRANSA = 'N' or 'n', and is m otherwise.         ==
  ! ==   Before entry with TRANSA = 'N' or 'n', the leading m by k  ==
  ! ==   part of the array A must contain the matrix A, otherwise   ==
  ! ==   the leading k by m part of the array A must contain the    ==
  ! ==   matrix A.                                                  ==
  ! ==   Unchanged on exit.                                         ==
  ! ==                                                              ==
  ! == LDA - INTEGER.                                               ==
  ! ==   On entry, LDA specifies the first dimension of A as        ==
  ! ==   declared in the calling (sub) program. When TRANSA = 'N'   ==
  ! ==   or 'n' then LDA must be at least max(1,m), otherwise LDA   ==
  ! ==   must be at least max(1,k).                                 ==
  ! ==   Unchanged on exit.                                         ==
  ! ==                                                              ==
  ! == B - DOUBLE PRECISION array of DIMENSION (LDB,kb), where kb   ==
  ! ==   is n when TRANSB = 'N' or 'n', and is k otherwise.         ==
  ! ==   Before entry with TRANSB = 'N' or 'n', the leading k by n  ==
  ! ==   part of the array B must contain the matrix B, otherwise   ==
  ! ==   the leading n by k part of the array B must contain the    ==
  ! ==   matrix B.                                                  ==
  ! ==   Unchanged on exit.                                         ==
  ! ==                                                              ==
  ! == LDB - INTEGER.                                               ==
  ! ==   On entry, LDB specifies the first dimension of B as        ==
  ! ==   declared in the calling (sub) program. When TRANSB = 'N'   ==
  ! ==   or 'n' then LDB must be at least max(1,k), otherwise LDB   ==
  ! ==   must be at least max(1,n).                                 ==
  ! ==   Unchanged on exit.                                         ==
  ! ==                                                              ==
  ! == BETA - DOUBLE PRECISION.                                     ==
  ! ==   On entry, BETA specifies the scalar beta. When BETA is     ==
  ! ==   supplied as zero then C need not be set on input.          ==
  ! ==   Unchanged on exit.                                         ==
  ! ==                                                              ==
  ! == C - DOUBLE PRECISION array of DIMENSION (LDC,n). Before      ==
  ! ==   entry, the leading m by n part of the array C must contain ==
  ! ==   the matrix C, except when beta is zero, in which case C    ==
  ! ==   need not be set on entry. On exit, the array C is          ==
  ! ==   overwritten by the m by n matrix                           ==
  ! ==   (alpha*op(A)*op(B)+beta*C).                                ==
  ! ==                                                              ==
  ! == LDC - INTEGER.                                               ==
  ! ==   On entry, LDC specifies the first dimension of C as        ==
  ! ==   declared in the calling (sub) program. LDC must be         ==
  ! ==   at least max(1,m).                                         ==
  ! ==   Unchanged on exit.                                         ==
  ! ==                                                              ==
  ! ==                                                              ==
  ! == Level 3 Blas routine.                                        ==
  ! ==                                                              ==
  ! ==  -- Written on 8-February-1989.                              ==
  ! ==     Jack Dongarra, Argonne National Laboratory.              ==
  ! ==     Iain Duff, AERE Harwell.                                 ==
  ! ==     Jeremy Du Croz, Numerical Algorithms Group Ltd.          ==
  ! ==     Sven Hammarling, Numerical Algorithms Group Ltd.         ==
  ! ==--------------------------------------------------------------==
  ! Scalar Arguments
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  CHARACTER(len=1) :: transa,transb
  INTEGER :: m,n,k,lda,ldb,ldc
  DOUBLE PRECISION alpha,beta
  ! Array Arguments
  DOUBLE PRECISION a(lda,*),b(ldb,*),c(ldc,*)
  ! External Functions
  LOGICAL :: lsame
  EXTERNAL         lsame
  ! External Subroutines
  EXTERNAL         xerbla
  ! Intrinsic Functions
  INTRINSIC        max
  ! Local Arrays 
  DOUBLE PRECISION a1(m,k),b1(k,n),c0(ldc,n),c1(n,m)
  ! Local Scalars
  LOGICAL :: nota,notb
  INTEGER :: i,info,j,l,ncola,nrowa,nrowb,ier
  DOUBLE PRECISION t0
  ! Parameters
  DOUBLE PRECISION one,zero
  PARAMETER       (one=1.0e+0_real_8,zero=0.0e+0_real_8)
  ! ==--------------------------------------------------------------==
  ! 
  ! Set NOTA and NOTB as true if A and B respectively are not
  ! transposed and set NROWA, NCOLA and NROWB as the number of 
  ! rows and columns of A and the number of rows of B respectively.
  ! 
  nota=lsame(transa,'N')
  notb=lsame(transb,'N')
  IF (nota) THEN
     nrowa=m
     ncola=k
  ELSE
     nrowa=k
     ncola=m
  ENDIF
  IF (notb) THEN
     nrowb=k
  ELSE
     nrowb=n
  ENDIF
  ! 
  ! Test the input parameters
  ! 
  info=0
  IF ((.NOT.nota).AND.&
       (.NOT.LSAME(TRANSA,'C')).AND.&
       (.NOT.LSAME(TRANSA,'T'))) THEN
     info=1
  ELSE IF((.NOT.notb).AND.&
       (.NOT.LSAME(TRANSB,'C')).AND.&
       (.NOT.LSAME(TRANSB,'T'))) THEN
     info=2
  ELSE IF (m.LT.0) THEN
     info=3
  ELSE IF (n.LT.0) THEN
     info=4
  ELSE IF (k.LT.0)THEN
     info=5
  ELSE IF (lda.LT.MAX(1,nrowa)) THEN
     info=8
  ELSE IF (ldb.LT.MAX(1,nrowb)) THEN
     info=10
  ELSE IF (ldc.LT.MAX(1,m)) THEN
     info=13
  ENDIF
  IF (info.NE.0) THEN
     CALL xerbla('DGEMM ',info)
     RETURN
  ENDIF
  ! 
  ! Quick return if possible
  ! 
  IF ((m.EQ.0).OR.(n.EQ.0).OR.&
       (((ALPHA.EQ.ZERO).OR.(K.EQ.0)).AND.(BETA.EQ.ONE)))&
       RETURN
  ! 
  ! And if alpha.eq.zero
  ! 
  IF (alpha.EQ.zero) THEN
     IF (beta.EQ.zero) THEN
        !$omp parallel do private(J,I)
        DO j=1,n
           DO i=1,m
              c(i,j)=zero
           ENDDO
        ENDDO
     ELSE
        !$omp parallel do private(J,I)
        DO j=1,n
           DO i=1,m
              c(i,j)=beta*c(i,j)
           ENDDO
        ENDDO
     ENDIF
     RETURN
  ENDIF
  ! 
  ! Start the operations
  ! 
  IF (notb) THEN
     IF (nota) THEN
        IF (beta.EQ.zero) THEN
           ! ==--------------------------------------------------------------==
           ! == Form C := alpha*A*B                                          ==
           ! ==--------------------------------------------------------------==
           IF (k.LT.2.OR.m.LT.2.OR.n.LT.2) THEN
              !$omp parallel do private(j,i,l,t0)
              DO j=1,n
                 DO i=1,m
                    c(i,j)=zero
                 ENDDO
                 DO l=1,k
                    t0=alpha*b(l,j)
                    DO i=1,m
                       c(i,j)=c(i,j)+t0*a(i,l)
                    ENDDO
                 ENDDO
              ENDDO
           ELSE
              CALL hdmffm(a,m,k,lda,b,n,ldb,c,ldc,ier)
              !$omp parallel do private(j,i)
              DO j=1,n
                 DO i=1,m
                    c(i,j)=alpha*c(i,j)
                 ENDDO
              ENDDO
           ENDIF
        ELSE
           ! ==--------------------------------------------------------------==
           ! == Form C := alpha*A*B+beta*C                                   ==
           ! ==--------------------------------------------------------------==
           IF (k.LT.2.OR.m.LT.2.OR.n.LT.2) THEN
              !$omp parallel do private(j,i,l,t0)
              DO j=1,n
                 DO i=1,m
                    c(i,j)=beta*c(i,j)
                 ENDDO
                 DO l=1,k
                    t0=alpha*b(l,j)
                    DO i=1,m
                       c(i,j)=c(i,j)+t0*a(i,l)
                    ENDDO
                 ENDDO
              ENDDO
           ELSE
              CALL hdmffm(a,m,k,lda,b,n,ldb,c0,ldc,ier)
              !$omp parallel do private(j,i)
              DO j=1,n
                 DO i=1,m
                    c(i,j)=alpha*c0(i,j)+beta*c(i,j)
                 ENDDO
              ENDDO
           ENDIF
        ENDIF
     ELSE
        IF (beta.EQ.zero) THEN
           ! ==--------------------------------------------------------------==
           ! == Form C := alpha*A'*B                                         ==
           ! ==--------------------------------------------------------------==
           IF (l.LT.2.OR.m.LT.2.OR.n.LT.2) THEN
              !$omp parallel do private(j,i,l,t0)
              DO j=1,n
                 DO i=1,m
                    t0=zero
                    DO l=1,k
                       t0=t0+a(l,i)*b(l,j)
                    ENDDO
                    c(i,j)=alpha*t0
                 ENDDO
              ENDDO
           ELSE
              !$omp parallel do private(i,l)
              DO i=1,m
                 DO l=1,k
                    a1(i,l)=a(l,i)
                 ENDDO
              ENDDO
              CALL hdmffm(a1,m,k,m,b,n,ldb,c,ldc,ier)
              !$omp parallel do private(j,i)
              DO j=1,n
                 DO i=1,m
                    c(i,j)=alpha*c(i,j)
                 ENDDO
              ENDDO
           ENDIF
        ELSE
           ! ==--------------------------------------------------------------==
           ! == Form C := alpha*A'*B+beta*C                                  ==
           ! ==--------------------------------------------------------------==
           IF (k.LT.2.OR.m.LT.2.OR.n.LT.2) THEN
              !$omp parallel do private(j,i,l,t0)
              DO j=1,n
                 DO i=1,m
                    t0=zero
                    DO l=1,k
                       t0=t0+a(l,i)*b(l,j)
                    ENDDO
                    c(i,j)=alpha*t0+beta*c(i,j)
                 ENDDO
              ENDDO
           ELSE
              !$omp parallel do private(i,l)
              DO i=1,m
                 DO l=1,k
                    a1(i,l)=a(l,i)
                 ENDDO
              ENDDO
              CALL hdmffm(a1,m,k,m,b,n,ldb,c0,ldc,ier)
              !$omp parallel do private(j,i)
              DO j=1,n
                 DO i=1,m
                    c(i,j)=alpha*c0(i,j)+beta*c(i,j)
                 ENDDO
              ENDDO
           ENDIF
        ENDIF
     ENDIF
  ELSE
     IF (nota)THEN
        IF (beta.EQ.zero) THEN
           ! ==--------------------------------------------------------------==
           ! == Form C := alpha*A*B'                                         ==
           ! ==--------------------------------------------------------------==
           IF (k.LT.2.OR.m.LT.2.OR.n.LT.2) THEN
              !$omp parallel do private(j,i,l,t0)
              DO j=1,n
                 DO i=1,m
                    c(i,j)=zero
                 ENDDO
                 DO l=1,k
                    t0=alpha*b(j,l)
                    DO i=1,m
                       c(i,j)=c(i,j)+t0*a(i,l)
                    ENDDO
                 ENDDO
              ENDDO
           ELSE
              !$omp parallel do private(j,l)
              DO j=1,n
                 DO l=1,k
                    b1(l,j)=b(j,l)
                 ENDDO
              ENDDO
              CALL hdmffm(a,m,k,lda,b1,n,k,c,ldc,ier)
              !$omp parallel do private(j,i)
              DO j=1,n
                 DO i=1,m
                    c(i,j)=alpha*c(i,j)
                 ENDDO
              ENDDO
           ENDIF
        ELSE
           ! ==--------------------------------------------------------------==
           ! == Form C := alpha*A*B'+beta*C                                  ==
           ! ==--------------------------------------------------------------==
           IF (k.LT.2.OR.m.LT.2.OR.n.LT.2) THEN
              !$omp parallel do private(j,i,l,t0)
              DO j=1,n
                 DO i=1,m
                    c(i,j)=beta*c(i,j)
                 ENDDO
                 DO l=1,k
                    t0=alpha*b(j,l)
                    DO i=1,m
                       c(i,j)=c(i,j)+t0*a(i,l)
                    ENDDO
                 ENDDO
              ENDDO
           ELSE
              !$omp parallel do private(j,l)
              DO j=1,n
                 DO l=1,k
                    b1(l,j)=b(j,l)
                 ENDDO
              ENDDO
              CALL hdmffm(a,m,k,lda,b1,n,k,c0,ldc,ier)
              !$omp parallel do private(j,i)
              DO j=1,n
                 DO i=1,m
                    c(i,j)=alpha*c0(i,j)+beta*c(i,j)
                 ENDDO
              ENDDO
           ENDIF
        ENDIF
     ELSE
        IF (beta.EQ.zero) THEN
           ! ==--------------------------------------------------------------==
           ! == Form  C := alpha*A'*B'                                       ==
           ! ==--------------------------------------------------------------==
           IF (k.LT.2.OR.m.LT.2.OR.n.LT.2) THEN
              !$omp parallel do private(j,i,l,t0)
              DO j=1,n
                 DO i=1,m
                    t0=zero
                    DO l=1,k
                       t0=t0+a(l,i)*b(j,l)
                    ENDDO
                    c(i,j)=alpha*t0
                 ENDDO
              ENDDO
           ELSE
              CALL hdmffm(b,n,k,ldb,a,m,lda,c1,n,ier)
              !$omp parallel do private(j,i)
              DO j=1,n
                 DO i=1,m
                    c(i,j)=alpha*c1(j,i)
                 ENDDO
              ENDDO
           ENDIF
        ELSE
           ! ==--------------------------------------------------------------==
           ! == Form C := alpha*A'*B'+beta*C                                 ==
           ! ==--------------------------------------------------------------==
           IF (k.LT.2.OR.m.LT.2.OR.n.LT.2) THEN
              !$omp parallel do private(j,i,l,t0)
              DO j=1,n
                 DO i=1,m
                    t0=zero
                    DO l=1,k
                       t0=t0+a(l,i)*b(j,l)
                    ENDDO
                    c(i,j)=alpha*t0+beta*c(i,j)
                 ENDDO
              ENDDO
           ELSE
              CALL hdmffm(b,n,k,ldb,a,m,lda,c1,n,ier)
              !$omp parallel do private(j,i)
              DO j=1,n
                 DO i=1,m
                    c(i,j)=alpha*c1(j,i)+beta*c(i,j)
                 ENDDO
              ENDDO
           ENDIF
        ENDIF
     ENDIF
  ENDIF
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE dgemm
! ==================================================================

! ==================================================================
SUBROUTINE dgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
  ! ==--------------------------------------------------------------==
  ! == Purpose                                                      ==
  ! == =======                                                      ==
  ! ==                                                              ==
  ! == DGEMV performs one of the matrix-vector operations           ==
  ! ==                                                              ==
  ! == y := alpha*A*x + beta*y, or y := alpha*A'*x + beta*y,        ==
  ! ==                                                              ==
  ! == where alpha and beta are scalars, x and y are vectors        ==
  ! == and A is an m by n matrix.                                   ==
  ! ==                                                              ==
  ! == Parameters                                                   ==
  ! == ==========                                                   ==
  ! ==                                                              ==
  ! == TRANS  - CHARACTER*1.                                        ==
  ! ==    On entry, TRANS specifies the operation to be performed   ==
  ! ==    as follows:                                               ==
  ! ==                                                              ==
  ! ==       TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.          ==
  ! ==                                                              ==
  ! ==       TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.         ==
  ! ==                                                              ==
  ! ==       TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.         ==
  ! ==                                                              ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == M      - INTEGER.                                            ==
  ! ==    On entry, M specifies the number of rows of the matrix A. ==
  ! ==    M must be at least zero.                                  ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == N      - INTEGER.                                            ==
  ! ==    On entry, N specifies the number of columns of            ==
  ! ==    the matrix A. N must be at least zero.                    ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == ALPHA  - DOUBLE PRECISION.                                   ==
  ! ==    On entry, ALPHA specifies the scalar alpha.               ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == A      - DOUBLE PRECISION array of DIMENSION (LDA,n).        ==
  ! ==    Before entry, the leading m by n part of the array A      ==
  ! ==    must contain the matrix of coefficients.                  ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == LDA    - INTEGER.                                            ==
  ! ==    On entry, LDA specifies the first dimension of A          ==
  ! ==    as declared in the calling (sub) program.                 ==
  ! ==    LDA must be at least max(1,m).                            ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == X      - DOUBLE PRECISION array of DIMENSION at least        ==
  ! ==    (1+(n-1)*abs(INCX)) when TRANS = 'N' or 'n' and           ==
  ! ==    at least (1+(m-1)*abs(INCX)) otherwise. Before entry,     ==
  ! ==    the incremented array X must contain the vector x.        ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == INCX   - INTEGER.                                            ==
  ! ==    On entry, INCX specifies the increment for the elements   ==
  ! ==    of X. INCX must not be zero.                              ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == BETA   - DOUBLE PRECISION.                                   ==
  ! ==    On entry, BETA specifies the scalar beta. When BETA is    ==
  ! ==    supplied as zero then Y need not be set on input.         ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == Y      - DOUBLE PRECISION array of DIMENSION at least        ==
  ! ==    (1+(m-1)*abs(INCY)) when TRANS = 'N' or 'n' and at least  ==
  ! ==    (1+(n-1)*abs(INCY)) otherwise. Before entry with BETA     ==
  ! ==    non-zero, the incremented array Y must contain the vector ==
  ! ==    y. On exit, Y is overwritten by the updated vector y.     ==
  ! ==                                                              ==
  ! == INCY   - INTEGER.                                            ==
  ! ==    On entry, INCY specifies the increment for the elements   ==
  ! ==    of Y. INCY must not be zero.                              ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! ==                                                              ==
  ! == Level 2 Blas routine.                                        ==
  ! ==                                                              ==
  ! == -- Written on 22-October-1986.                               ==
  ! ==    Jack Dongarra, Argonne National Lab.                      ==
  ! ==    Jeremy Du Croz, Nag Central Office.                       ==
  ! ==    Sven Hammarling, Nag Central Office.                      ==
  ! ==    Richard Hanson, Sandia National Labs                      ==
  ! ==--------------------------------------------------------------==
  ! Scalar Arguments
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  DOUBLE PRECISION alpha,beta
  INTEGER :: incx,incy,lda,m,n
  CHARACTER(len=1) :: trans
  ! Array Arguments
  DOUBLE PRECISION a(lda,*),x(*),y(*)
  ! Parameters
  DOUBLE PRECISION one,zero
  PARAMETER       (one=1.0e+0_real_8,zero=0.0e+0_real_8)
  ! Local Scalars
  DOUBLE PRECISION temp
  INTEGER :: i,info,ix,iy,j,jx,jy,kx,ky,lenx,leny
  ! External Functions
  LOGICAL :: lsame
  EXTERNAL         lsame
  ! External Subroutines
  EXTERNAL         xerbla
  ! Intrinsic Functions
  INTRINSIC        max
  ! Auxiliary variables
  DOUBLE PRECISION t0,t1,t2,t3,t4,t5,t6,t7
  INTEGER :: modj,jx0,jx1,jx2,jx3,jx4,jx5,jx6,jx7
  INTEGER :: ijy(n),iix(m),ijx(n),iiy(m)
  ! ==---------------------------------------------------------------=
  ! 
  ! Test the input parameters
  ! 
  info=0
  IF (.NOT.lsame(trans,'N').AND.&
       .NOT.LSAME(TRANS,'T').AND.&
       .NOT.LSAME(TRANS,'C')) THEN
     info=1
  ELSE IF (m.LT.0) THEN
     info=2
  ELSE IF (n.LT.0) THEN
     info=3
  ELSE IF (lda.LT.MAX(1,m)) THEN
     info=6
  ELSE IF (incx.EQ.0) THEN
     info=8
  ELSE IF (incy.EQ.0) THEN
     info=11
  ENDIF
  IF (info.NE.0) THEN
     CALL xerbla('DGEMV ',info)
     RETURN
  ENDIF
  ! 
  ! Quick return if possible 
  ! 
  IF ((m.EQ.0).OR.(n.EQ.0).OR.&
       ((ALPHA.EQ.ZERO).AND.(BETA.EQ.ONE)))  RETURN
  ! 
  ! Set LENX and LENY, the lengths of the vectors x and y, 
  ! and set up the start points in X and Y.
  ! 
  IF (lsame(trans,'N')) THEN
     lenx=n
     leny=m
  ELSE
     lenx=m
     leny=n
  ENDIF
  IF (incx.GT.0) THEN
     kx=1
  ELSE
     kx=1-(lenx-1)*incx
  ENDIF
  IF (incy.GT.0) THEN
     ky=1
  ELSE
     ky=1-(leny-1)*incy
  ENDIF
  ! 
  ! Start the operations. In this version the elements of A 
  ! are accessed sequentially with one pass through A.
  ! 
  ! First form y := beta*y.
  ! 
  IF (beta.NE.one) THEN
     IF (incy.EQ.1) THEN
        IF (beta.EQ.zero) THEN
           !$omp parallel do private(I)
           DO i=1,leny
              y(i)=zero
           ENDDO
        ELSE
           !$omp parallel do private(I)
           DO i=1,leny
              y(i)=beta*y(i)
           ENDDO
        ENDIF
     ELSE
        IF (beta.EQ.zero) THEN
           !$omp parallel do private(I,IY)
           DO i=1,leny
              iy=ky+(i-1)*incy
              y(iy)=zero
           ENDDO
        ELSE
           !$omp parallel do private(I,IY)
           !OCL NOVREC(Y)
           !CDIR NODEP(Y)
           DO i=1,leny
              iy=ky+(i-1)*incy
              y(iy)=beta*y(iy)
           ENDDO
        ENDIF
     ENDIF
  ENDIF
  IF (alpha.EQ.zero) RETURN
  IF (lsame(trans,'N')) THEN
     ! ==--------------------------------------------------------------==
     ! == Form y := alpha*A*x + y                                      ==
     ! ==--------------------------------------------------------------==
     IF (incy.EQ.1) THEN
        IF (incx.EQ.1) THEN
           modj=MOD(n,8)
           DO j=1,modj
              temp=alpha*x(j)
              !$omp parallel do private(I)
              DO i=1,m
                 y(i)=y(i)+temp*a(i,j)
              ENDDO
           ENDDO
           DO j=modj+1,n,8
              t0=alpha*x(j  )
              t1=alpha*x(j+1)
              t2=alpha*x(j+2)
              t3=alpha*x(j+3)
              t4=alpha*x(j+4)
              t5=alpha*x(j+5)
              t6=alpha*x(j+6)
              t7=alpha*x(j+7)
              !$omp parallel do private(I)
              DO i=1,m
                 y(i)=y(i)&
                      +t0*a(i,j  )+t1*a(i,j+1)+t2*a(i,j+2)+t3*a(i,j+3)&
                      +t4*a(i,j+4)+t5*a(i,j+5)+t6*a(i,j+6)+t7*a(i,j+7)
              ENDDO
           ENDDO
        ELSE
           ! ..Fill list vector IJX
           ! ..IJX(1:N) -> {KX,KX+INCX,KX+2*INCX,...,KX+(N-1)*INCX}
           DO j=1,n
              ijx(j)=kx+(j-1)*incx
           ENDDO
           modj=MOD(n,8)
           DO j=1,modj
              jx=ijx(j)
              temp=alpha*x(jx)
              !$omp parallel do private(I)
              DO i=1,m
                 y(i)=y(i)+temp*a(i,j)
              ENDDO
           ENDDO
           DO j=modj+1,n,8
              jx0=ijx(j  )
              jx1=ijx(j+1)
              jx2=ijx(j+2)
              jx3=ijx(j+3)
              jx4=ijx(j+4)
              jx5=ijx(j+5)
              jx6=ijx(j+6)
              jx7=ijx(j+7)
              t0=alpha*x(jx0)
              t1=alpha*x(jx1)
              t2=alpha*x(jx2)
              t3=alpha*x(jx3)
              t4=alpha*x(jx4)
              t5=alpha*x(jx5)
              t6=alpha*x(jx6)
              t7=alpha*x(jx7)
              !$omp parallel do private(I)
              DO i=1,m
                 y(i)=y(i)&
                      +t0*a(i,j  )+t1*a(i,j+1)+t2*a(i,j+2)+t3*a(i,j+3)&
                      +t4*a(i,j+4)+t5*a(i,j+5)+t6*a(i,j+6)+t7*a(i,j+7)
              ENDDO
           ENDDO
        ENDIF
     ELSE
        ! ..Fill list vector IIY
        ! ..IIY(1:M) -> {KY,KY+INCY,KY+2*INCY,...,KY+(M-1)*INCY}
        DO j=1,m
           iiy(j)=ky+(j-1)*incy
        ENDDO
        IF (incx.EQ.1) THEN
           !OCL NOVREC(Y)
           DO j=1,n
              temp=alpha*x(j)
              !$omp parallel do private(I,IY)
              !CDIR NODEP(Y)
              DO i=1,m
                 iy=iiy(i)
                 y(iy)=y(iy)+temp*a(i,j)
              ENDDO
           ENDDO
        ELSE
           ! ..Fill list vector IJX
           ! ..IJX(1:N) -> {KX,KX+INCX,KX+2*INCX,...,KX+(N-1)*INCX}
           DO j=1,n
              ijx(j)=kx+(j-1)*incx
           ENDDO
           !OCL NOVREC(Y)
           DO j=1,n
              jx=ijx(j)
              temp=alpha*x(jx)
              !$omp parallel do private(I,IY)
              !CDIR NODEP(Y)
              DO i=1,m
                 iy=iiy(i)
                 y(iy)=y(iy)+temp*a(i,j)
              ENDDO
           ENDDO
        ENDIF
     ENDIF
  ELSE
     ! ==--------------------------------------------------------------==
     ! == Form y := alpha*A'*x + y                                     ==
     ! ==--------------------------------------------------------------==
     ! ..Fill list vector IJY
     ! ..IJY(1:N) -> {KY,KY+INCY,KY+2*INCY,...,KY+(N-1)*INCY}
     DO j=1,n
        ijy(j)=ky+(j-1)*incy
     ENDDO
     IF (incx.EQ.1) THEN
        !$omp parallel do private(J,I,JY,TEMP) 
        !OCL NOVREC(Y)
        !CDIR NODEP(Y)
        DO j=1,n
           temp=zero
           DO i=1,m
              temp=temp+a(i,j)*x(i)
           ENDDO
           jy=ijy(j)
           y(jy)=y(jy)+alpha*temp
        ENDDO
     ELSE
        ! ..Fill list vector IIX
        ! ..IIX(1:M) -> {KX,KX+INCX,KX+2*INCX,...,KX+(M-1)*INCX}
        DO j=1,m
           iix(j)=kx+(j-1)*incx
        ENDDO
        !$omp parallel do private(J,I,IX,JY,TEMP) 
        !OCL NOVREC(Y)
        !CDIR NODEP(Y)
        DO j=1,n
           temp=zero
           DO i=1,m
              ix=iix(i)
              temp=temp+a(i,j)*x(ix)
           ENDDO
           jy=ijy(j)
           y(jy)=y(jy)+alpha*temp
        ENDDO
     ENDIF
  ENDIF
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE dgemv
! ==================================================================

! ==================================================================
SUBROUTINE dger(m,n,alpha,x,incx,y,incy,a,lda)
  ! ==--------------------------------------------------------------==
  ! == Purpose                                                      ==
  ! == =======                                                      ==
  ! ==                                                              ==
  ! == DGER   performs the rank 1 operation                         ==
  ! ==                                                              ==
  ! == A := alpha*x*y' + A,                                         ==
  ! ==                                                              ==
  ! == where alpha is a scalar, x is an m element vector,           ==
  ! == y is an n element vector and A is an m by n matrix.          ==
  ! ==                                                              ==
  ! == Parameters                                                   ==
  ! == ==========                                                   ==
  ! ==                                                              ==
  ! == M      - INTEGER.                                            ==
  ! ==    On entry, M specifies the number of rows of the matrix A. ==
  ! ==    M must be at least zero.                                  ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == N      - INTEGER.                                            ==
  ! ==    On entry, N specifies the number of columns               ==
  ! ==    of the matrix A. N must be at least zero.                 ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == ALPHA  - DOUBLE PRECISION.                                   ==
  ! ==    On entry, ALPHA specifies the scalar alpha.               ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == X      - DOUBLE PRECISION array of dimension at least        ==
  ! ==    (1+(m-1)*abs(INCX)). Before entry, the incremented        ==
  ! ==    array X must contain the m element vector x.              ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == INCX   - INTEGER.                                            ==
  ! ==    On entry, INCX specifies the increment for the elements   ==
  ! ==    of X. INCX must not be zero.                              ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == Y      - DOUBLE PRECISION array of dimension at least        ==
  ! ==    (1+(n-1)*abs(INCY)). Before entry, the incremented        ==
  ! ==    array Y must contain the n element vector y.              ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == INCY   - INTEGER.                                            ==
  ! ==    On entry, INCY specifies the increment for the elements   ==
  ! ==    of Y. INCY must not be zero.                              ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == A      - DOUBLE PRECISION array of DIMENSION (LDA,n).        ==
  ! ==    Before entry, the leading m by n part of the array A      ==
  ! ==    must contain the matrix of coefficients. On exit, A       ==
  ! ==    is overwritten by the updated matrix.                     ==
  ! ==                                                              ==
  ! == LDA    - INTEGER.                                            ==
  ! ==    On entry, LDA specifies the first dimension of A as       ==
  ! ==    declared in the calling (sub) program. LDA must be        ==
  ! ==    at least max(1,m).                                        ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! ==                                                              ==
  ! == Level 2 Blas routine.                                        ==
  ! ==                                                              ==
  ! == -- Written on 22-October-1986.                               ==
  ! ==    Jack Dongarra, Argonne National Lab.                      ==
  ! ==    Jeremy Du Croz, Nag Central Office.                       ==
  ! ==    Sven Hammarling, Nag Central Office.                      ==
  ! ==    Richard Hanson, Sandia National Labs.                     ==
  ! ==--------------------------------------------------------------==
  ! Scalar Arguments
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  DOUBLE PRECISION alpha
  INTEGER :: incx,incy,lda,m,n
  ! Array Arguments
  DOUBLE PRECISION a(lda,*),x(*),y(*)
  ! Parameters
  DOUBLE PRECISION zero
  PARAMETER       (zero=0.0e+0_real_8)
  ! Local Scalars
  DOUBLE PRECISION temp
  INTEGER :: i,info,ix,j,jy,kx
  ! External Subroutines
  EXTERNAL         xerbla
  ! Intrinsic Functions
  INTRINSIC        max
  ! ==--------------------------------------------------------------==
  ! 
  ! Test the input parameters
  ! 
  info=0
  IF (m.LT.0) THEN
     info=1
  ELSE IF (n.LT.0) THEN
     info=2
  ELSE IF (incx.EQ.0) THEN
     info=5
  ELSE IF (incy.EQ.0) THEN
     info=7
  ELSE IF (lda.LT.MAX(1,m)) THEN
     info=9
  ENDIF
  IF (info.NE.0) THEN
     CALL xerbla('DGER  ',info)
     RETURN
  ENDIF
  ! 
  ! Quick return if possible
  ! 
  IF ((m.EQ.0).OR.(n.EQ.0).OR.(alpha.EQ.zero)) RETURN
  ! 
  ! Start the operations. In this version the elements of A 
  ! are accessed sequentially with one pass through A.
  ! 
  IF (incy.GT.0) THEN
     jy=1
  ELSE
     jy=1-(n-1)*incy
  ENDIF
  IF (incx.EQ.1) THEN
     !$omp parallel do private(J,I,TEMP)
     DO j=1,n
        temp=alpha*y(jy+(j-1)*incy)
        DO i=1,m
           a(i,j)=a(i,j)+x(i)*temp
        ENDDO
     ENDDO
  ELSE
     IF (incx.GT.0) THEN
        kx=1
     ELSE
        kx=1-(m-1)*incx
     ENDIF
     !$omp parallel do private(J,I,TEMP)
     DO j=1,n
        temp=alpha*y(jy+(j-1)*incy)
        DO i=1,m
           a(i,j)=a(i,j)+x(kx+(i-1)*incx)*temp
        ENDDO
     ENDDO
  ENDIF
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE dger
! ==================================================================

! ==================================================================
SUBROUTINE dsyrk(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
  ! ==--------------------------------------------------------------==
  ! == Purpose                                                      ==
  ! == =======                                                      ==
  ! ==                                                              ==
  ! == DSYRK  performs one of the symmetric rank k operations       ==
  ! ==                                                              ==
  ! == C := alpha*A*A' + beta*C,                                    ==
  ! ==                                                              ==
  ! == or                                                           ==
  ! ==                                                              ==
  ! == C := alpha*A'*A + beta*C,                                    ==
  ! ==                                                              ==
  ! == where alpha and beta  are scalars, C is an n by n symmetric  ==
  ! == matrix and A is an n by k matrix in the first case and       ==
  ! == a k by n matrix in the second case.                          ==
  ! ==                                                              ==
  ! == Parameters                                                   ==
  ! == ==========                                                   ==
  ! ==                                                              ==
  ! == UPLO   - CHARACTER*1.                                        ==
  ! ==    On entry, UPLO specifies whether the upper or lower       ==
  ! ==    triangular part of the array C is to be referenced        ==
  ! ==    as follows:                                               ==
  ! ==                                                              ==
  ! ==       UPLO = 'U' or 'u' Only the upper triangular part of C  ==
  ! ==                         is to be referenced.                 ==
  ! ==                                                              ==
  ! ==       UPLO = 'L' or 'l' Only the lower triangular part of C  ==
  ! ==                         is to be referenced.                 ==
  ! ==                                                              ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == TRANS  - CHARACTER*1.                                        ==
  ! ==    On entry, TRANS specifies the operation to be performed   ==
  ! ==    as follows:                                               ==
  ! ==                                                              ==
  ! ==       TRANS = 'N' or 'n'   C := alpha*A*A' + beta*C.         ==
  ! ==                                                              ==
  ! ==       TRANS = 'T' or 't'   C := alpha*A'*A + beta*C.         ==
  ! ==                                                              ==
  ! ==       TRANS = 'C' or 'c'   C := alpha*A'*A + beta*C.         ==
  ! ==                                                              ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == N      - INTEGER.                                            ==
  ! ==    On entry, N specifies the order of the matrix C.          ==
  ! ==    N must be at least zero.                                  ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == K      - INTEGER.                                            ==
  ! ==    On entry with TRANS = 'N' or 'n', K specifies the number  ==
  ! ==    of columns of the matrix A, and on entry with             ==
  ! ==    TRANS = 'T' or 't' or 'C' or 'c', K specifies the number  ==
  ! ==    of rows of the matrix  A.  K must be at least zero.       ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == ALPHA  - DOUBLE PRECISION.                                   ==
  ! ==    On entry, ALPHA specifies the scalar alpha.               ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == A      - DOUBLE PRECISION array of DIMENSION (LDA,ka),       ==
  ! ==    where ka is k when TRANS = 'N' or 'n', and is n otherwise.==
  ! ==    Before entry with TRANS = 'N' or 'n', the leading n by k  ==
  ! ==    part of the array A must contain the matrix A, otherwise  ==
  ! ==    the leading k by n part of the array A must contain the   ==
  ! ==    matrix A.                                                 ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == LDA    - INTEGER.                                            ==
  ! ==    On entry, LDA specifies the first dimension of A as       ==
  ! ==    declared in the calling (sub) program. When TRANS = 'N'   ==
  ! ==    or 'n' then LDA must be at least max(1,n), otherwise      ==
  ! ==    LDA must be at least max(1 k).                            ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == BETA   - DOUBLE PRECISION.                                   ==
  ! ==    On entry, BETA specifies the scalar beta.                 ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == C      - DOUBLE PRECISION array of DIMENSION (LDC,n).        ==
  ! ==    Before entry with UPLO = 'U' or 'u', the leading n by n   ==
  ! ==    upper triangular part of the array C must contain the     ==
  ! ==    upper triangular part of the symmetric matrix and the     ==
  ! ==    strictly lower triangular part of C is not referenced.    ==
  ! ==    On exit, the upper triangular part of the array C is      ==
  ! ==    overwritten by the upper triangular part of the updated   ==
  ! ==    matrix. Before entry with UPLO = 'L' or 'l', the leading  ==
  ! ==    n by n lower triangular part of the array C must contain  ==
  ! ==    the lower triangular part of the symmetric matrix and the ==
  ! ==    strictly upper triangular part of C is not referenced.    ==
  ! ==    On exit, the lower triangular part of the array C is      ==
  ! ==    overwritten by the lower triangular part of the updated   ==
  ! ==    matrix.                                                   ==
  ! ==                                                              ==
  ! == LDC    - INTEGER.                                            ==
  ! ==    On entry, LDC specifies the first dimension of C as       ==
  ! ==    declared in the calling (sub) program. LDC must be at     ==
  ! ==    least max(1,n).                                           ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! ==                                                              ==
  ! == Level 3 Blas routine.                                        ==
  ! ==                                                              ==
  ! == -- Written on 8-February-1989.                               ==
  ! ==    Jack Dongarra, Argonne National Laboratory.               ==
  ! ==    Iain Duff, AERE Harwell.                                  ==
  ! ==    Jeremy Du Croz, Numerical Algorithms Group Ltd.           ==
  ! ==    Sven Hammarling, Numerical Algorithms Group Ltd.          ==
  ! ==--------------------------------------------------------------==
  ! Scalar Arguments
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  CHARACTER(len=1) :: uplo,trans
  INTEGER :: n,k,lda,ldc
  DOUBLE PRECISION alpha,beta
  ! Array Arguments
  DOUBLE PRECISION a(lda,*),c(ldc,*)
  ! External Functions
  LOGICAL :: lsame
  EXTERNAL         lsame
  ! External Subroutines
  EXTERNAL         xerbla
  ! Intrinsic Functions
  INTRINSIC        max
  ! Local Scalars
  LOGICAL :: upper
  INTEGER :: i,info,j,l,nrowa
  DOUBLE PRECISION temp
  ! Parameters
  DOUBLE PRECISION one,zero
  PARAMETER       (one=1.0e+0_real_8,zero=0.0e+0_real_8)
  ! ==--------------------------------------------------------------==
  ! 
  ! Test the input parameters
  ! 
  IF (lsame(trans,'N')) THEN
     nrowa=n
  ELSE
     nrowa=k
  ENDIF
  upper=lsame(uplo,'U')
  ! 
  info=0
  IF ((.NOT.upper).AND.&
       (.NOT.LSAME(UPLO,'L'))) THEN
     info=1
  ELSE IF((.NOT.lsame(trans,'N')).AND.&
       (.NOT.LSAME(TRANS,'T')).AND.&
       (.NOT.LSAME(TRANS,'C'))) THEN
     info=2
  ELSE IF (n.LT.0) THEN
     info=3
  ELSE IF (k.LT.0) THEN
     info=4
  ELSE IF (lda.LT.MAX(1,nrowa)) THEN
     info=7
  ELSE IF (ldc.LT.MAX(1,n)) THEN
     info=10
  ENDIF
  IF (info.NE.0) THEN
     CALL xerbla('DSYRK ',info)
     RETURN
  ENDIF
  ! 
  ! Quick return if possible
  ! 
  IF ((n.EQ.0).OR.&
       (((ALPHA.EQ.ZERO).OR.(K.EQ.0)).AND.(BETA.EQ.ONE)))&
       RETURN
  ! 
  ! And when alpha.eq.zero
  ! 
  IF (alpha.EQ.zero) THEN
     IF (upper) THEN
        IF (beta.EQ.zero) THEN
           !$omp parallel do private(J,I) schedule(static,1)
           DO j=1,n
              DO i=1,j
                 c(i,j)=zero
              ENDDO
           ENDDO
        ELSE
           !$omp parallel do private(J,I) schedule(static,1)
           DO j=1,n
              DO i=1,j
                 c(i,j)=beta*c(i,j)
              ENDDO
           ENDDO
        ENDIF
     ELSE
        IF (beta.EQ.zero) THEN
           !$omp parallel do private(J,I) schedule(static,1)
           DO j=1,n
              DO i=j,n
                 c(i,j)=zero
              ENDDO
           ENDDO
        ELSE
           !$omp parallel do private(J,I) schedule(static,1)
           DO j=1,n
              DO i=j,n
                 c(i,j)=beta*c(i,j)
              ENDDO
           ENDDO
        ENDIF
     ENDIF
     RETURN
  ENDIF
  ! 
  ! Start the operations
  ! 
  IF (lsame(trans,'N')) THEN
     IF (upper) THEN
        IF (beta.EQ.zero) THEN
           ! ==--------------------------------------------------------------==
           ! ==  Form C := alpha*A*A' (upper triangle referenced)            ==
           ! ==--------------------------------------------------------------==
           !$omp parallel do private(J,I,L,TEMP) schedule(static,1)
           DO j=1,n
              DO i=1,j
                 c(i,j)=zero
              ENDDO
              DO l=1,k
                 temp=alpha*a(j,l)
                 DO i=1,j
                    c(i,j)=c(i,j)+temp*a(i,l)
                 ENDDO
              ENDDO
           ENDDO
        ELSE IF (beta.EQ.one) THEN
           ! ==--------------------------------------------------------------==
           ! ==  Form C := alpha*A*A' + C (upper triangle referenced)        ==
           ! ==--------------------------------------------------------------==
           !$omp parallel do private(J,I,L,TEMP) schedule(static,1)
           DO j=1,n
              DO l=1,k
                 temp=alpha*a(j,l)
                 DO i=1,j
                    c(i,j)=c(i,j)+temp*a(i,l)
                 ENDDO
              ENDDO
           ENDDO
        ELSE
           ! ==--------------------------------------------------------------==
           ! ==  Form C := alpha*A*A' + beta*C (upper triangle referenced)   ==
           ! ==--------------------------------------------------------------==
           !$omp parallel do private(J,I,L,TEMP) schedule(static,1)
           DO j=1,n
              DO i=1,j
                 c(i,j)=beta*c(i,j)
              ENDDO
              DO l=1,k
                 temp=alpha*a(j,l)
                 DO i=1,j
                    c(i,j)=c(i,j)+temp*a(i,l)
                 ENDDO
              ENDDO
           ENDDO
        ENDIF
     ELSE
        IF (beta.EQ.zero) THEN
           ! ==--------------------------------------------------------------==
           ! ==  Form C := alpha*A*A' (lower triangle referenced)            ==
           ! ==--------------------------------------------------------------==
           !$omp parallel do private(J,I,L,TEMP) schedule(static,1)
           DO j=1,n
              DO i=j,n
                 c(i,j)=zero
              ENDDO
              DO l=1,k
                 temp=alpha*a(j,l)
                 DO i=j,n
                    c(i,j)=c(i,j)+temp*a(i,l)
                 ENDDO
              ENDDO
           ENDDO
        ELSE IF (beta.EQ.one) THEN
           ! ==--------------------------------------------------------------==
           ! ==  Form C := alpha*A*A' + C (lower triangle referenced)        ==
           ! ==--------------------------------------------------------------==
           !$omp parallel do private(J,I,L,TEMP) schedule(static,1)
           DO j=1,n
              DO l=1,k
                 temp=alpha*a(j,l)
                 DO i=j,n
                    c(i,j)=c(i,j)+temp*a(i,l)
                 ENDDO
              ENDDO
           ENDDO
        ELSE
           ! ==--------------------------------------------------------------==
           ! ==  Form C := alpha*A*A' + beta*C (lower triangle referenced)   ==
           ! ==--------------------------------------------------------------==
           !$omp parallel do private(J,I,L,TEMP) schedule(static,1)
           DO j=1,n
              DO i=j,n
                 c(i,j)=beta*c(i,j)
              ENDDO
              DO l=1,k
                 temp=alpha*a(j,l)
                 DO i=j,n
                    c(i,j)=c(i,j)+temp*a(i,l)
                 ENDDO
              ENDDO
           ENDDO
        ENDIF
     ENDIF
  ELSE
     IF (upper) THEN
        IF (beta.EQ.zero) THEN
           ! ==--------------------------------------------------------------==
           ! ==  Form C := alpha*A'*A (upper triangle referenced)            ==
           ! ==--------------------------------------------------------------==
           !$omp parallel do private(J,I,L,TEMP) schedule(static,1)
           DO j=1,n
              DO i=1,j
                 temp=zero
                 DO l=1,k
                    temp=temp+a(l,i)*a(l,j)
                 ENDDO
                 c(i,j)=alpha*temp
              ENDDO
           ENDDO
        ELSE IF (beta.EQ.one) THEN
           ! ==--------------------------------------------------------------==
           ! ==  Form C := alpha*A'*A + C (upper triangle referenced)        ==
           ! ==--------------------------------------------------------------==
           !$omp parallel do private(J,I,L,TEMP) schedule(static,1)
           DO j=1,n
              DO i=1,j
                 temp=zero
                 DO l=1,k
                    temp=temp+a(l,i)*a(l,j)
                 ENDDO
                 c(i,j)=alpha*temp+c(i,j)
              ENDDO
           ENDDO
        ELSE
           ! ==--------------------------------------------------------------==
           ! ==  Form C := alpha*A'*A + beta*C (upper triangle referenced)   ==
           ! ==--------------------------------------------------------------==
           !$omp parallel do private(J,I,L,TEMP) schedule(static,1)
           DO j=1,n
              DO i=1,j
                 temp=zero
                 DO l=1,k
                    temp=temp+a(l,i)*a(l,j)
                 ENDDO
                 c(i,j)=alpha*temp+beta*c(i,j)
              ENDDO
           ENDDO
        ENDIF
     ELSE
        IF (beta.EQ.zero) THEN
           ! ==--------------------------------------------------------------==
           ! ==  Form C := alpha*A'*A (lower triangle referenced)            ==
           ! ==--------------------------------------------------------------==
           !$omp parallel do private(J,I,L,TEMP) schedule(static,1)
           DO j=1,n
              DO i=j,n
                 temp=zero
                 DO l=1,k
                    temp=temp+a(l,i)*a(l,j)
                 ENDDO
                 c(i,j)=alpha*temp
              ENDDO
           ENDDO
        ELSE IF (beta.EQ.one) THEN
           ! ==--------------------------------------------------------------==
           ! ==  Form C := alpha*A'*A + C (lower triangle referenced)        ==
           ! ==--------------------------------------------------------------==
           !$omp parallel do private(J,I,L,TEMP) schedule(static,1)
           DO j=1,n
              DO i=j,n
                 temp=zero
                 DO l=1,k
                    temp=temp+a(l,i)*a(l,j)
                 ENDDO
                 c(i,j)=alpha*temp+c(i,j)
              ENDDO
           ENDDO
        ELSE
           ! ==--------------------------------------------------------------==
           ! ==  Form C := alpha*A'*A + beta*C (lower triangle referenced)   ==
           ! ==--------------------------------------------------------------==
           !$omp parallel do private(J,I,L,TEMP) schedule(static,1)
           DO j=1,n
              DO i=j,n
                 temp=zero
                 DO l=1,k
                    temp=temp+a(l,i)*a(l,j)
                 ENDDO
                 c(i,j)=alpha*temp+beta*c(i,j)
              ENDDO
           ENDDO
        ENDIF
     ENDIF
  ENDIF
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE dsyrk
! ==================================================================

! ==================================================================
SUBROUTINE dtrmm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
  ! ==--------------------------------------------------------------==
  ! == Purpose                                                      ==
  ! == =======                                                      ==
  ! ==                                                              ==
  ! == DTRMM performs one of the matrix-matrix operations           ==
  ! ==                                                              ==
  ! == B := alpha*op(A)*B,   or   B := alpha*B*op(A),               ==
  ! ==                                                              ==
  ! == where alpha is a scalar, B is an m by n matrix, A is a unit, ==
  ! == or non-unit, upper or lower triangular matrix and op(A) is   ==
  ! == one of                                                       ==
  ! ==                                                              ==
  ! == op(A) = A   or   op(A) = A'.                                 ==
  ! ==                                                              ==
  ! == Parameters                                                   ==
  ! == ==========                                                   ==
  ! ==                                                              ==
  ! == SIDE   - CHARACTER*1.                                        ==
  ! ==    On entry, SIDE specifies whether op(A) multiplies B       ==
  ! ==    form the left or right as follows:                        ==
  ! ==                                                              ==
  ! ==       SIDE = 'L' or 'l'   B := alpha*op(A)*B.                ==
  ! ==                                                              ==
  ! ==       SIDE = 'R' or 'r'   B := alpha*B*op(A).                ==
  ! ==                                                              ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == UPLO   - CHARACTER*1.                                        ==
  ! ==    On entry, UPLO specifies whether the matrix A is an upper ==
  ! ==    or lower triangular matrix as follows:                    ==
  ! ==                                                              ==
  ! ==       UPLO = 'U' or 'u'  A is an upper triangular matrix.    ==
  ! ==                                                              ==
  ! ==       UPLO = 'L' or 'l'  A is a lower triangular matrix.     ==
  ! ==                                                              ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == TRANSA - CHARACTER*1.                                        ==
  ! ==    On entry, TRANSA specifies the form of op(A) to be        ==
  ! ==    used in the matrix multiplication as follows:             ==
  ! ==                                                              ==
  ! ==       TRANSA = 'N' or 'n'   op(A) = A.                       ==
  ! ==                                                              ==
  ! ==       TRANSA = 'T' or 't'   op(A) = A'.                      ==
  ! ==                                                              ==
  ! ==       TRANSA = 'C' or 'c'   op(A) = A'.                      ==
  ! ==                                                              ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == DIAG   - CHARACTER*1.                                        ==
  ! ==    On entry, DIAG specifies whether or not A is unit         ==
  ! ==    triangular as follows:                                    ==
  ! ==                                                              ==
  ! ==       DIAG = 'U' or 'u' A is assumed to be unit triangular.  ==
  ! ==                                                              ==
  ! ==       DIAG = 'N' or 'n' A is not assumed to be unit          ==
  ! ==                         triangular.                          ==
  ! ==                                                              ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == M      - INTEGER.                                            ==
  ! ==    On entry, M specifies the number of rows of B. M must     ==
  ! ==    be at least zero.                                         ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == N      - INTEGER.                                            ==
  ! ==    On entry, N specifies the number of columns of B. N must  ==
  ! ==    be at least zero.                                         ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == ALPHA  - DOUBLE PRECISION.                                   ==
  ! ==    On entry, ALPHA specifies the scalar alpha. When alpha    ==
  ! ==    is zero then A is not referenced and B need not be set    ==
  ! ==    before entry.                                             ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == A      - DOUBLE PRECISION array of DIMENSION (LDA,k),        ==
  ! ==    where k is m when SIDE = 'L' or 'l' and is n when         ==
  ! ==    SIDE = 'R' or 'r'. Before entry with UPLO = 'U' or 'u',   ==
  ! ==    the leading k by k upper triangular part of the array A   ==
  ! ==    must contain the upper triangular matrix and the strictly ==
  ! ==    lower triangular part of A is not referenced.             ==
  ! ==    Before entry with UPLO = 'L' or 'l', the leading k by k   ==
  ! ==    lower triangular part of the array A must contain the     ==
  ! ==    lower triangular matrix and the strictly upper triangular ==
  ! ==    part of A is not referenced.                              ==
  ! ==    Note that when DIAG = 'U' or 'u', the diagonal elements   ==
  ! ==    of A are not referenced either, but are assumed to be     ==
  ! ==    unity.                                                    ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == LDA    - INTEGER.                                            ==
  ! ==    On entry, LDA specifies the first dimension of A as       ==
  ! ==    declared in the calling (sub) program. When SIDE = 'L'    ==
  ! ==    or 'l'  then LDA  must be at least max(1,m), when         ==
  ! ==    SIDE = 'R' or 'r' then LDA must be at least max(1,n).     ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == B      - DOUBLE PRECISION array of DIMENSION (LDB,n).        ==
  ! ==    Before entry, the leading m by n part of the array B must ==
  ! ==    contain the matrix B, and on exit is overwritten by the   ==
  ! ==    transformed matrix.                                       ==
  ! ==                                                              ==
  ! == LDB    - INTEGER.                                            ==
  ! ==    On entry, LDB specifies the first dimension of B as       ==
  ! ==    declared in the calling (sub) program. LDB must be at     ==
  ! ==    least max(1,m).                                           ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! ==                                                              ==
  ! == Level 3 Blas routine.                                        ==
  ! ==                                                              ==
  ! == -- Written on 8-February-1989.                               ==
  ! ==    Jack Dongarra, Argonne National Laboratory.               ==
  ! ==    Iain Duff, AERE Harwell.                                  ==
  ! ==    Jeremy Du Croz, Numerical Algorithms Group Ltd.           ==
  ! ==    Sven Hammarling, Numerical Algorithms Group Ltd.          ==
  ! ==--------------------------------------------------------------==
  ! Scalar Arguments
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  CHARACTER(len=1) :: side,uplo,transa,diag
  INTEGER :: m,n,lda,ldb
  DOUBLE PRECISION alpha
  ! Array Arguments
  DOUBLE PRECISION a(lda,*),b(ldb,*)
  ! External Functions
  LOGICAL :: lsame
  EXTERNAL         lsame
  ! External Subroutines
  EXTERNAL         xerbla
  ! Intrinsic Functions
  INTRINSIC        max
  ! Local Scalars 
  LOGICAL :: lside,nounit,upper
  INTEGER :: i,info,j,k,nrowa
  DOUBLE PRECISION temp
  ! Parameters
  DOUBLE PRECISION one,zero
  PARAMETER       (one=1.0e+0_real_8,zero=0.0e+0_real_8)
  ! ==--------------------------------------------------------------==
  ! 
  ! Test the input parameters
  ! 
  lside=lsame(side,'L')
  IF (lside) THEN
     nrowa=m
  ELSE
     nrowa=n
  ENDIF
  nounit=lsame(diag,'N')
  upper =lsame(uplo,'U')
  ! 
  info=0
  IF ((.NOT.lside).AND.&
       (.NOT.LSAME(SIDE,'R'))) THEN
     info=1
  ELSE IF((.NOT.upper).AND.&
       (.NOT.LSAME(UPLO,'L'))) THEN
     info=2
  ELSE IF((.NOT.lsame(transa,'N')).AND.&
       (.NOT.LSAME(TRANSA,'T')).AND.&
       (.NOT.LSAME(TRANSA,'C'))) THEN
     info=3
  ELSE IF((.NOT.lsame(diag,'U')).AND.&
       (.NOT.LSAME(DIAG,'N'))) THEN
     info=4
  ELSE IF (m.LT.0) THEN
     info=5
  ELSE IF (n.LT.0) THEN
     info=6
  ELSE IF (lda.LT.MAX(1,nrowa)) THEN
     info=9
  ELSE IF (ldb.LT.MAX(1,m)) THEN
     info=11
  ENDIF
  IF (info.NE.0) THEN
     CALL xerbla('DTRMM ',info)
     RETURN
  ENDIF
  ! 
  ! Quick return if possible
  ! 
  IF (n.EQ.0) RETURN
  ! 
  ! And when alpha.eq.zero
  ! 
  IF (alpha.EQ.zero) THEN
     !$omp parallel do private(J,I)
     DO j=1,n
        DO i=1,m
           b(i,j)=zero
        ENDDO
     ENDDO
     RETURN
  ENDIF
  ! 
  ! Start the operations
  ! 
  IF (lside) THEN
     IF (lsame(transa,'N')) THEN
        ! ==--------------------------------------------------------------==
        ! ==  Form B := alpha*A*B                                         ==
        ! ==--------------------------------------------------------------==
        IF (upper) THEN
           !$omp parallel do private(J,K,I,TEMP) 
           DO j=1,n
              DO k=1,m
                 ! !              IF(B(K,J).NE.ZERO) THEN
                 temp=alpha*b(k,j)
                 DO i=1,k-1
                    b(i,j)=b(i,j)+temp*a(i,k)
                 ENDDO
                 IF (nounit) temp=temp*a(k,k)
                 b(k,j)=temp
                 ! !              ENDIF
              ENDDO
           ENDDO
        ELSE
           !$omp parallel do private(J,K,I,TEMP) 
           DO j=1,n
              DO k=m,1,-1
                 ! !              IF(B(K,J).NE.ZERO) THEN
                 temp=alpha*b(k,j)
                 b(k,j)=temp
                 IF (nounit) b(k,j)=b(k,j)*a(k,k)
                 DO i=k+1,m
                    b(i,j)=b(i,j)+temp*a(i,k)
                 ENDDO
                 ! !              ENDIF
              ENDDO
           ENDDO
        ENDIF
     ELSE
        ! ==--------------------------------------------------------------==
        ! ==  Form B := alpha*A'*B                                        ==
        ! ==--------------------------------------------------------------==
        IF (upper) THEN
           !$omp parallel do private(J,I,K,TEMP) 
           DO j=1,n
              DO i=m,1,-1
                 temp=b(i,j)
                 IF (nounit) temp=temp*a(i,i)
                 DO k=1,i-1
                    temp=temp+a(k,i)*b(k,j)
                 ENDDO
                 b(i,j)=alpha*temp
              ENDDO
           ENDDO
        ELSE
           !$omp parallel do private(J,I,K,TEMP) 
           DO j=1,n
              DO i=1,m
                 temp=b(i,j)
                 IF (nounit) temp=temp*a(i,i)
                 DO k=i+1,m
                    temp=temp+a(k,i)*b(k,j)
                 ENDDO
                 b(i,j)=alpha*temp
              ENDDO
           ENDDO
        ENDIF
     ENDIF
  ELSE
     IF (lsame(transa,'N')) THEN
        ! ==--------------------------------------------------------------==
        ! ==  Form B := alpha*B*A                                         ==
        ! ==--------------------------------------------------------------==
        IF (upper) THEN
           !$omp parallel do private(J,I,K,TEMP) schedule(static,1)
           DO j=n,1,-1
              temp=alpha
              IF (nounit) temp=temp*a(j,j)
              DO i=1,m
                 b(i,j)=temp*b(i,j)
              ENDDO
              DO k=1,j-1
                 temp=alpha*a(k,j)
                 DO i=1,m
                    b(i,j)=b(i,j)+temp*b(i,k)
                 ENDDO
              ENDDO
           ENDDO
        ELSE
           !$omp parallel do private(J,I,K,TEMP) schedule(static,1)
           DO j=1,n
              temp=alpha
              IF (nounit) temp=temp*a(j,j)
              DO i=1,m
                 b(i,j)=temp*b(i,j)
              ENDDO
              DO k=j+1,n
                 temp=alpha*a(k,j)
                 DO i=1,m
                    b(i,j)=b(i,j)+temp*b(i,k)
                 ENDDO
              ENDDO
           ENDDO
        ENDIF
     ELSE
        ! ==--------------------------------------------------------------==
        ! ==  Form B := alpha*B*A'                                        ==
        ! ==--------------------------------------------------------------==
        IF (upper) THEN
           !$omp parallel do private(J,I,K,TEMP) schedule(static,1)
           DO k=1,n
              DO j=1,k-1
                 temp=alpha*a(j,k)
                 DO i=1,m
                    b(i,j)=b(i,j)+temp*b(i,k)
                 ENDDO
              ENDDO
              temp=alpha
              IF (nounit) temp=temp*a(k,k)
              IF (temp.NE.one) THEN
                 DO i=1,m
                    b(i,k)=temp*b(i,k)
                 ENDDO
              ENDIF
           ENDDO
        ELSE
           !$omp parallel do private(J,I,K,TEMP) schedule(static,1)
           DO k=n,1,-1
              DO j=k+1,n
                 temp=alpha*a(j,k)
                 DO i=1,m
                    b(i,j)=b(i,j)+temp*b(i,k)
                 ENDDO
              ENDDO
              temp=alpha
              IF (nounit) temp=temp*a(k,k)
              IF (temp.NE.one) THEN
                 DO i=1,m
                    b(i,k)=temp*b(i,k)
                 ENDDO
              ENDIF
           ENDDO
        ENDIF
     ENDIF
  ENDIF
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE dtrmm
! ==================================================================

! ==================================================================
SUBROUTINE zgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
  ! ==--------------------------------------------------------------==
  ! == Purpose                                                      ==
  ! == =======                                                      ==
  ! ==                                                              ==
  ! == ZGEMM performs one of the matrix-matrix operations           ==
  ! ==                                                              ==
  ! == C := alpha*op(A)*op(B) + beta*C,                             ==
  ! ==                                                              ==
  ! == where op(X) is one of                                        ==
  ! ==                                                              ==
  ! == op(X) = X  or  op(X) = X'  or  op(X) = conjg(X'),            ==
  ! ==                                                              ==
  ! == alpha and beta are scalars, and A, B and C are matrices,     ==
  ! == with op(A) an m by k matrix, op(B) a k by n matrix and C     ==
  ! == an m by n matrix.                                            ==
  ! ==                                                              ==
  ! == Parameters                                                   ==
  ! == ==========                                                   ==
  ! ==                                                              ==
  ! == TRANSA - CHARACTER*1.                                        ==
  ! ==    On entry,TRANSA specifies the form of op(A) to be         ==
  ! ==    used in the matrix multiplication as follows:             ==
  ! ==                                                              ==
  ! ==       TRANSA = 'N' or 'n',  op(A) = A.                       ==
  ! ==                                                              ==
  ! ==       TRANSA = 'T' or 't',  op(A) = A'.                      ==
  ! ==                                                              ==
  ! ==       TRANSA = 'C' or 'c',  op(A) = conjg(A').               ==
  ! ==                                                              ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == TRANSB - CHARACTER*1.                                        ==
  ! ==    On entry, TRANSB specifies the form of op(B) to be        == 
  ! ==    used in the matrix multiplication as follows:             ==
  ! ==                                                              ==
  ! ==       TRANSB = 'N' or 'n',  op(B) = B.                       ==
  ! ==                                                              ==
  ! ==       TRANSB = 'T' or 't',  op(B) = B'.                      ==
  ! ==                                                              ==
  ! ==       TRANSB = 'C' or 'c',  op(B) = conjg(B').               ==
  ! ==                                                              ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == M      - INTEGER.                                            ==
  ! ==    On entry, M specifies the number of rows of the matrix    ==
  ! ==    op(A) and of the matrix C. M must be at least zero.       ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == N      - INTEGER.                                            ==
  ! ==    On entry, N specifies the number of columns of the matrix ==
  ! ==    op(B) and the number of columns of the matrix C. N must   ==
  ! ==    be at least zero.                                         ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == K      - INTEGER.                                            ==
  ! ==    On entry, K specifies the number of columns of the matrix ==
  ! ==    op(A) and the number of rows of the matrix op(B).         ==
  ! ==    K must be at least zero.                                  ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == ALPHA  - complex(8) :: .                                         ==
  ! ==    On entry, ALPHA specifies the scalar alpha.               ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == A      - complex(8) :: array of DIMENSION (LDA,ka), where ka    ==      
  ! ==    is k when TRANSA = 'N' or 'n', and is m otherwise.        ==
  ! ==    Before entry with TRANSA = 'N' or 'n', the leading m by k ==
  ! ==    part of the array A must contain the matrix A, otherwise  ==
  ! ==    the leading k by m part of the array A must contain the   ==
  ! ==    matrix A.                                                 ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == LDA    - INTEGER.                                            ==
  ! ==    On entry, LDA specifies the first dimension of A as       ==
  ! ==    declared in the calling (sub) program. When               ==
  ! ==    TRANSA = 'N' or 'n' then LDA must be at least             ==
  ! ==    max(1,m), otherwise LDA must be at least max(1,k).        ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == B      - complex(8) :: array of DIMENSION (LDB,kb), where kb    ==
  ! ==    is n when TRANSB = 'N' or 'n', and is k otherwise.        ==
  ! ==    Before entry with TRANSB = 'N' or 'n', the leading k by n ==
  ! ==    part of the array B must contain the matrix B, otherwise  ==
  ! ==    the leading n by k part of the array B must contain the   ==
  ! ==    matrix B.                                                 ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == LDB    - INTEGER.                                            ==
  ! ==    On entry, LDB specifies the first dimension of B as       ==
  ! ==    declared in the calling (sub) program. When TRANSB = 'N'  ==
  ! ==    or 'n' then LDB must be at least max(1,k), otherwise      ==
  ! ==    LDB must be at least max(1,n).                            ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == BETA   - complex(8) :: .                                         ==
  ! ==    On entry, BETA specifies the scalar beta. When BETA is    ==
  ! ==    supplied as zero then C need not be set on input.         ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == C      - complex(8) :: array of DIMENSION (LDC,n).              ==
  ! ==    Before entry, the leading m by n part of the array C must ==
  ! ==    contain the matrix C, except when beta is zero, in which  ==
  ! ==    case C need not be set on entry. On exit, the array C is  ==
  ! ==    overwritten by the m by n matrix                          ==
  ! ==    (alpha*op(A)*op(B) + beta*C).                             ==
  ! ==                                                              ==
  ! == LDC    - INTEGER.                                            ==
  ! ==    On entry, LDC specifies the first dimension of C as       ==
  ! ==    declared in the calling (sub) program. LDC must be at     ==
  ! ==    least max(1,m).                                           ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! ==                                                              ==
  ! == Level 3 Blas routine.                                        ==
  ! ==                                                              ==
  ! == -- Written on 8-February-1989.                               ==
  ! ==    Jack Dongarra, Argonne National Laboratory.               ==
  ! ==    Iain Duff, AERE Harwell.                                  ==
  ! ==    Jeremy Du Croz, Numerical Algorithms Group Ltd.           ==
  ! ==    Sven Hammarling, Numerical Algorithms Group Ltd.          ==
  ! ==--------------------------------------------------------------==
  ! Scalar Arguments
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  CHARACTER(len=1) :: transa,transb
  INTEGER :: m,n,k,lda,ldb,ldc
  COMPLEX(real_8) :: alpha,beta
  ! Array Arguments
  COMPLEX(real_8) :: a(lda,*),b(ldb,*),c(ldc,*)
  ! External Functions
  LOGICAL :: lsame
  EXTERNAL    lsame
  ! External Subroutines
  EXTERNAL    xerbla
  ! Intrinsic Functions
  INTRINSIC   conjg,max
  ! Local Scalars
  LOGICAL :: conja,conjb,nota,notb
  INTEGER :: i,info,j,l,ncola,nrowa,nrowb
  COMPLEX(real_8) :: temp
  ! Parameters
  COMPLEX(real_8) :: one
  PARAMETER  (one=(1.0e+0_real_8,0.0e+0_real_8))
  COMPLEX(real_8) :: zero
  PARAMETER  (zero=(0.0e+0_real_8,0.0e+0_real_8))
  ! ==--------------------------------------------------------------==
  ! 
  ! Set NOTA and NOTB as true if A and B respectively are not
  ! conjugated or transposed, set CONJA and CONJB as true if A and
  ! B respectively are to be transposed but not conjugated and set
  ! NROWA, NCOLA and NROWB as the number of rows and columns of A
  ! and the number of rows of B respectively.
  ! 
  nota =lsame(transa,'N')
  notb =lsame(transb,'N')
  conja=lsame(transa,'C')
  conjb=lsame(transb,'C')
  IF (nota) THEN
     nrowa=m
     ncola=k
  ELSE
     nrowa=k
     ncola=m
  ENDIF
  IF (notb) THEN
     nrowb=k
  ELSE
     nrowb=n
  ENDIF
  ! 
  ! Test the input parameters
  ! 
  info=0
  IF ((.NOT.nota).AND.&
       (.NOT.CONJA).AND.&
       (.NOT.LSAME(TRANSA,'T'))) THEN
     info=1
  ELSE IF((.NOT.notb).AND.&
       (.NOT.CONJB).AND.&
       (.NOT.LSAME(TRANSB,'T'))) THEN
     info=2
  ELSE IF (m.LT.0) THEN
     info=3
  ELSE IF (n.LT.0) THEN
     info=4
  ELSE IF (k.LT.0) THEN
     info=5
  ELSE IF (lda.LT.MAX(1,nrowa)) THEN
     info=8
  ELSE IF (ldb.LT.MAX(1,nrowb)) THEN
     info=10
  ELSE IF (ldc.LT.MAX(1,m)) THEN
     info=13
  ENDIF
  IF (info.NE.0) THEN
     CALL xerbla('ZGEMM ',info)
     RETURN
  ENDIF
  ! 
  ! Quick return if possible
  ! 
  IF ((m.EQ.0).OR.(n.EQ.0).OR.&
       (((ALPHA.EQ.ZERO).OR.(K.EQ.0)).AND.(BETA.EQ.ONE)))&
       RETURN
  ! 
  ! And when alpha.eq.zero
  ! 
  IF (alpha.EQ.zero) THEN
     IF (beta.EQ.zero) THEN
        !$omp parallel do private(J,I) 
        DO j=1,n
           DO i=1,m
              c(i,j)=zero
           ENDDO
        ENDDO
     ELSE
        !$omp parallel do private(J,I) 
        DO j=1,n
           DO i=1,m
              c(i,j)=beta*c(i,j)
           ENDDO
        ENDDO
     ENDIF
     RETURN
  ENDIF
  ! 
  ! Start the operations
  ! 
  IF (notb) THEN
     IF (nota) THEN
        ! ==--------------------------------------------------------------==
        ! ==  Form C := alpha*A*B + beta*C                                ==
        ! ==--------------------------------------------------------------==
        IF (beta.EQ.zero) THEN
           !$omp parallel do private(J,I,L,TEMP)
           DO j=1,n
              DO i=1,m
                 c(i,j)=zero
              ENDDO
              DO l=1,k
                 temp=alpha*b(l,j)
                 DO i=1,m
                    c(i,j)=c(i,j)+temp*a(i,l)
                 ENDDO
              ENDDO
           ENDDO
        ELSE IF (beta.EQ.one) THEN
           !$omp parallel do private(J,I,L,TEMP)
           DO j=1,n
              DO l=1,k
                 temp=alpha*b(l,j)
                 DO i=1,m
                    c(i,j)=c(i,j)+temp*a(i,l)
                 ENDDO
              ENDDO
           ENDDO
        ELSE
           !$omp parallel do private(J,I,L,TEMP)
           DO j=1,n
              DO i=1,m
                 c(i,j)=beta*c(i,j)
              ENDDO
              DO l=1,k
                 temp=alpha*b(l,j)
                 DO i=1,m
                    c(i,j)=c(i,j)+temp*a(i,l)
                 ENDDO
              ENDDO
           ENDDO
        ENDIF
     ELSE IF (conja) THEN
        ! ==--------------------------------------------------------------==
        ! ==  Form C := alpha*conjg(A')*B + beta*C                        ==
        ! ==--------------------------------------------------------------==
        IF (beta.EQ.zero) THEN
           !$omp parallel do private(J,I,L,TEMP)
           DO j=1,n
              DO i=1,m
                 temp=zero
                 DO l=1,k
                    temp=temp+CONJG(a(l,i))*b(l,j)
                 ENDDO
                 c(i,j)=alpha*temp
              ENDDO
           ENDDO
        ELSE IF (beta.EQ.one) THEN
           !$omp parallel do private(J,I,L,TEMP)
           DO j=1,n
              DO i=1,m
                 temp=zero
                 DO l=1,k
                    temp=temp+CONJG(a(l,i))*b(l,j)
                 ENDDO
                 c(i,j)=alpha*temp+c(i,j)
              ENDDO
           ENDDO
        ELSE
           !$omp parallel do private(J,I,L,TEMP)
           DO j=1,n
              DO i=1,m
                 temp=zero
                 DO l=1,k
                    temp=temp+CONJG(a(l,i))*b(l,j)
                 ENDDO
                 c(i,j)=alpha*temp+beta*c(i,j)
              ENDDO
           ENDDO
        ENDIF
     ELSE
        ! ==--------------------------------------------------------------==
        ! ==  Form C := alpha*A'*B + beta*C                               ==
        ! ==--------------------------------------------------------------==
        IF (beta.EQ.zero) THEN
           !$omp parallel do private(J,I,L,TEMP)
           DO j=1,n
              DO i=1,m
                 temp=zero
                 DO l=1,k
                    temp=temp+a(l,i)*b(l,j)
                 ENDDO
                 c(i,j)=alpha*temp
              ENDDO
           ENDDO
        ELSE IF (beta.EQ.one) THEN
           !$omp parallel do private(J,I,L,TEMP)
           DO j=1,n
              DO i=1,m
                 temp=zero
                 DO l=1,k
                    temp=temp+a(l,i)*b(l,j)
                 ENDDO
                 c(i,j)=alpha*temp+c(i,j)
              ENDDO
           ENDDO
        ELSE
           !$omp parallel do private(J,I,L,TEMP)
           DO j=1,n
              DO i=1,m
                 temp=zero
                 DO l=1,k
                    temp=temp+a(l,i)*b(l,j)
                 ENDDO
                 c(i,j)=alpha*temp+beta*c(i,j)
              ENDDO
           ENDDO
        ENDIF
     ENDIF
  ELSE IF (nota) THEN
     IF (conjb) THEN
        ! ==--------------------------------------------------------------==
        ! ==  Form C := alpha*A*conjg(B') + beta*C                        ==
        ! ==--------------------------------------------------------------==
        IF (beta.EQ.zero) THEN
           !$omp parallel do private(J,I,L,TEMP)
           DO j=1,n
              DO i=1,m
                 c(i,j)=zero
              ENDDO
              DO l=1,k
                 temp=alpha*CONJG(b(j,l))
                 DO i=1,m
                    c(i,j)=c(i,j)+temp*a(i,l)
                 ENDDO
              ENDDO
           ENDDO
        ELSE IF (beta.EQ.one) THEN
           !$omp parallel do private(J,I,L,TEMP)
           DO j=1,n
              DO l=1,k
                 temp=alpha*CONJG(b(j,l))
                 DO i=1,m
                    c(i,j)=c(i,j)+temp*a(i,l)
                 ENDDO
              ENDDO
           ENDDO
        ELSE
           !$omp parallel do private(J,I,L,TEMP)
           DO j=1,n
              DO i=1,m
                 c(i,j)=beta*c(i,j)
              ENDDO
              DO l=1,k
                 temp=alpha*CONJG(b(j,l))
                 DO i=1,m
                    c(i,j)=c(i,j)+temp*a(i,l)
                 ENDDO
              ENDDO
           ENDDO
        ENDIF
     ELSE
        ! ==--------------------------------------------------------------==
        ! ==  Form C := alpha*A*B' + beta*C                               ==
        ! ==--------------------------------------------------------------==
        IF (beta.EQ.zero) THEN
           !$omp parallel do private(J,I,L,TEMP)
           DO j=1,n
              DO i=1,m
                 c(i,j)=zero
              ENDDO
              DO l=1,k
                 temp=alpha*b(j,l)
                 DO i=1,m
                    c(i,j)=c(i,j)+temp*a(i,l)
                 ENDDO
              ENDDO
           ENDDO
        ELSE IF (beta.EQ.one) THEN
           !$omp parallel do private(J,I,L,TEMP)
           DO j=1,n
              DO l=1,k
                 temp=alpha*b(j,l)
                 DO i=1,m
                    c(i,j)=c(i,j)+temp*a(i,l)
                 ENDDO
              ENDDO
           ENDDO
        ELSE
           !$omp parallel do private(J,I,L,TEMP)
           DO j=1,n
              DO i=1,m
                 c(i,j)=beta*c(i,j)
              ENDDO
              DO l=1,k
                 temp=alpha*b(j,l)
                 DO i=1,m
                    c(i,j)=c(i,j)+temp*a(i,l)
                 ENDDO
              ENDDO
           ENDDO
        ENDIF
     ENDIF
  ELSE IF (conja) THEN
     IF (conjb) THEN
        ! ==--------------------------------------------------------------==
        ! ==  Form C := alpha*conjg(A')*conjg(B') + beta*C                ==
        ! ==--------------------------------------------------------------==
        IF (beta.EQ.zero) THEN
           !$omp parallel do private(J,I,L,TEMP)
           DO j=1,n
              DO i=1,m
                 temp=zero
                 DO l=1,k
                    temp=temp+CONJG(a(l,i))*CONJG(b(j,l))
                 ENDDO
                 c(i,j)=alpha*temp
              ENDDO
           ENDDO
        ELSE IF (beta.EQ.one) THEN
           !$omp parallel do private(J,I,L,TEMP)
           DO j=1,n
              DO i=1,m
                 temp=zero
                 DO l=1,k
                    temp=temp+CONJG(a(l,i))*CONJG(b(j,l))
                 ENDDO
                 c(i,j)=alpha*temp+c(i,j)
              ENDDO
           ENDDO
        ELSE
           !$omp parallel do private(J,I,L,TEMP)
           DO j=1,n
              DO i=1,m
                 temp=zero
                 DO l=1,k
                    temp=temp+CONJG(a(l,i))*CONJG(b(j,l))
                 ENDDO
                 c(i,j)=alpha*temp+beta*c(i,j)
              ENDDO
           ENDDO
        ENDIF
     ELSE
        ! ==--------------------------------------------------------------==
        ! ==  Form C := alpha*conjg(A')*B' + beta*C                       ==
        ! ==--------------------------------------------------------------==
        IF (beta.EQ.zero) THEN
           !$omp parallel do private(J,I,L,TEMP)
           DO j=1,n
              DO i=1,m
                 temp=zero
                 DO l=1,k
                    temp=temp+CONJG(a(l,i))*b(j,l)
                 ENDDO
                 c(i,j)=alpha*temp
              ENDDO
           ENDDO
        ELSE IF (beta.EQ.one) THEN
           !$omp parallel do private(J,I,L,TEMP)
           DO j=1,n
              DO i=1,m
                 temp=zero
                 DO l=1,k
                    temp=temp+CONJG(a(l,i))*b(j,l)
                 ENDDO
                 c(i,j)=alpha*temp+c(i,j)
              ENDDO
           ENDDO
        ELSE
           !$omp parallel do private(J,I,L,TEMP)
           DO j=1,n
              DO i=1,m
                 temp=zero
                 DO l=1,k
                    temp=temp+CONJG(a(l,i))*b(j,l)
                 ENDDO
                 c(i,j)=alpha*temp+beta*c(i,j)
              ENDDO
           ENDDO
        ENDIF
     ENDIF
  ELSE
     IF (conjb) THEN
        ! ==--------------------------------------------------------------==
        ! ==  Form C := alpha*A'*conjg(B') + beta*C                       ==
        ! ==--------------------------------------------------------------==
        IF (beta.EQ.zero) THEN
           !$omp parallel do private(J,I,L,TEMP)
           DO j=1,n
              DO i=1,m
                 temp=zero
                 DO l=1,k
                    temp=temp+a(l,i)*CONJG(b(j,l))
                 ENDDO
                 c(i,j)=alpha*temp
              ENDDO
           ENDDO
        ELSE IF (beta.EQ.one) THEN
           !$omp parallel do private(J,I,L,TEMP)
           DO j=1,n
              DO i=1,m
                 temp=zero
                 DO l=1,k
                    temp=temp+a(l,i)*CONJG(b(j,l))
                 ENDDO
                 c(i,j)=alpha*temp+c(i,j)
              ENDDO
           ENDDO
        ELSE
           !$omp parallel do private(J,I,L,TEMP)
           DO j=1,n
              DO i=1,m
                 temp=zero
                 DO l=1,k
                    temp=temp+a(l,i)*CONJG(b(j,l))
                 ENDDO
                 c(i,j)=alpha*temp+beta*c(i,j)
              ENDDO
           ENDDO
        ENDIF
     ELSE
        ! ==--------------------------------------------------------------==
        ! ==  Form C := alpha*A'*B' + beta*C                              ==
        ! ==--------------------------------------------------------------==
        IF (beta.EQ.zero) THEN
           !$omp parallel do private(J,I,L,TEMP)
           DO j=1,n
              DO i=1,m
                 temp=zero
                 DO l=1,k
                    temp=temp+a(l,i)*b(j,l)
                 ENDDO
                 c(i,j)=alpha*temp
              ENDDO
           ENDDO
        ELSE IF (beta.EQ.one) THEN
           !$omp parallel do private(J,I,L,TEMP)
           DO j=1,n
              DO i=1,m
                 temp=zero
                 DO l=1,k
                    temp=temp+a(l,i)*b(j,l)
                 ENDDO
                 c(i,j)=alpha*temp+c(i,j)
              ENDDO
           ENDDO
        ELSE
           !$omp parallel do private(J,I,L,TEMP)
           DO j=1,n
              DO i=1,m
                 temp=zero
                 DO l=1,k
                    temp=temp+a(l,i)*b(j,l)
                 ENDDO
                 c(i,j)=alpha*temp+c(i,j)
              ENDDO
           ENDDO
        ENDIF
     ENDIF
  ENDIF
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE zgemm
! ==================================================================

! ==================================================================
SUBROUTINE zgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
  ! ==--------------------------------------------------------------==
  ! == Purpose                                                      ==
  ! == =======                                                      ==
  ! ==                                                              ==
  ! == ZGEMV performs one of the matrix-vector operations           ==
  ! ==                                                              ==
  ! == y := alpha*A*x + beta*y, or y := alpha*A'*x + beta*y, or     ==
  ! ==                                                              ==
  ! == y := alpha*conjg(A')*x + beta*y,                             ==
  ! ==                                                              ==
  ! == where alpha and beta are scalars, x and y are vectors        ==
  ! == and A is an m by n matrix.                                   ==
  ! ==                                                              ==
  ! == Parameters                                                   ==
  ! == ==========                                                   ==
  ! ==                                                              ==
  ! == TRANS  - CHARACTER*1.                                        ==
  ! ==    On entry, TRANS specifies the operation to be performed   ==
  ! ==    as follows:                                               ==
  ! ==                                                              ==
  ! ==       TRANS = 'N' or 'n'  y := alpha*A*x + beta*y.           ==
  ! ==                                                              ==
  ! ==       TRANS = 'T' or 't'  y := alpha*A'*x + beta*y.          ==
  ! ==                                                              ==
  ! ==       TRANS = 'C' or 'c'  y := alpha*conjg(A')*x + beta*y.   ==
  ! ==                                                              ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == M      - INTEGER.                                            ==
  ! ==    On entry, M specifies the number of rows of the matrix A. ==
  ! ==    M must be at least zero.                                  ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == N      - INTEGER.                                            ==
  ! ==    On entry, N specifies the number of columns of the matrix ==
  ! ==    A. N must be at least zero.                               ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == ALPHA  - complex(8) :: .                                         ==
  ! ==    On entry, ALPHA specifies the scalar alpha.               ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == A      - complex(8) :: array of DIMENSION (LDA,n).              ==
  ! ==    Before entry, the leading m by n part of the array A      ==
  ! ==    must contain the matrix of coefficients.                  ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == LDA    - INTEGER.                                            ==
  ! ==    On entry, LDA specifies the first dimension of A as       ==
  ! ==    declared in the calling (sub) program. LDA must be at     ==
  ! ==    least max(1,m).                                           ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == X      - complex(8) :: array of DIMENSION at least              ==
  ! ==    (1+(n-1)*abs(INCX)) when TRANS = 'N' or 'n' and at least  ==
  ! ==    (1+(m-1)*abs(INCX)) otherwise. Before entry, the          ==
  ! ==    incremented array X must contain the vector x.            ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == INCX   - INTEGER.                                            ==
  ! ==    On entry, INCX specifies the increment for the elements   ==
  ! ==    of X. INCX must not be zero.                              ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == BETA   - complex(8) :: .                                         ==
  ! ==    On entry, BETA specifies the scalar beta. When BETA is    ==
  ! ==    supplied as zero then Y need not be set on input.         ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! == Y      - complex(8) :: array of DIMENSION at least              ==
  ! ==    (1+(m-1)*abs(INCY)) when TRANS = 'N' or 'n' and at least  ==
  ! ==    (1+(n-1)*abs(INCY)) otherwise. Before entry with BETA     ==
  ! ==    non-zero, the incremented array Y must contain the vector ==
  ! ==    y. On exit, Y is overwritten by the updated vector y.     ==
  ! ==                                                              ==
  ! == INCY   - INTEGER.                                            ==
  ! ==    On entry, INCY specifies the increment for the elements   ==
  ! ==    of Y. INCY must not be zero.                              ==
  ! ==    Unchanged on exit.                                        ==
  ! ==                                                              ==
  ! ==                                                              ==
  ! == Level 2 Blas routine.                                        ==
  ! ==                                                              ==
  ! == -- Written on 22-October-1986.                               ==
  ! ==    Jack Dongarra, Argonne National Lab.                      ==
  ! ==    Jeremy Du Croz, Nag Central Office.                       ==
  ! ==    Sven Hammarling, Nag Central Office.                      ==
  ! ==    Richard Hanson, Sandia National Labs.                     ==
  ! ==--------------------------------------------------------------==
  ! Scalar Arguments
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  COMPLEX(real_8) :: alpha,beta
  INTEGER :: incx,incy,lda,m,n
  CHARACTER(len=1) :: trans
  ! Array Arguments
  COMPLEX(real_8) :: a(lda,*),x(*),y(*)
  ! Parameters
  COMPLEX(real_8) :: one
  PARAMETER  (one =(1.0e+0_real_8,0.0e+0_real_8))
  COMPLEX(real_8) :: zero
  PARAMETER  (zero=(0.0e+0_real_8,0.0e+0_real_8))
  ! Local Scalars
  COMPLEX(real_8) :: temp
  INTEGER :: i,info,ix,iy,j,jx,jy,kx,ky,lenx,leny
  LOGICAL :: noconj
  ! External Functions
  LOGICAL :: lsame
  EXTERNAL    lsame
  ! External Subroutines
  EXTERNAL    xerbla
  ! Intrinsic Functions
  INTRINSIC   conjg,max
  ! Auxiliary variables
  INTEGER :: ijy(n),iix(m),ijx(n),iiy(m)
  ! ==--------------------------------------------------------------==
  ! 
  ! Test the input parameters
  ! 
  info=0
  IF (.NOT.lsame(trans,'N').AND.&
       .NOT.LSAME(TRANS,'T').AND.&
       .NOT.LSAME(TRANS,'C')) THEN
     info=1
  ELSE IF (m.LT.0) THEN
     info=2
  ELSE IF (n.LT.0) THEN
     info=3
  ELSE IF (lda.LT.MAX(1,m)) THEN
     info=6
  ELSE IF (incx.EQ.0) THEN
     info=8
  ELSE IF (incy.EQ.0) THEN
     info=11
  ENDIF
  IF (info.NE.0) THEN
     CALL xerbla('ZGEMV ',info)
     RETURN
  ENDIF
  ! 
  ! Quick return if possible
  ! 
  IF ((m.EQ.0).OR.(n.EQ.0).OR.&
       ((ALPHA.EQ.ZERO).AND.(BETA.EQ.ONE))) RETURN
  ! 
  noconj=lsame(trans,'T')
  ! 
  ! Set LENX and LENY, the lengths of the vectors x and y, 
  ! and set up the start points in X and Y.
  ! 
  IF (lsame(trans,'N')) THEN
     lenx=n
     leny=m
  ELSE
     lenx=m
     leny=n
  ENDIF
  IF (incx.GT.0) THEN
     kx=1
  ELSE
     kx=1-(lenx-1)*incx
  ENDIF
  IF (incy.GT.0) THEN
     ky=1
  ELSE
     ky=1-(leny-1)*incy
  ENDIF
  ! 
  ! Start the operations. In this version the elements of A
  ! are accessed sequentially with one pass through A.
  ! 
  ! First form  y := beta*y.
  ! 
  IF (beta.NE.one) THEN
     IF (incy.EQ.1) THEN
        IF (beta.EQ.zero) THEN
           !$omp parallel do private(I)
           DO i=1,leny
              y(i)=zero
           ENDDO
        ELSE
           !$omp parallel do private(I)
           DO i=1,leny
              y(i)=beta*y(i)
           ENDDO
        ENDIF
     ELSE
        IF (beta.EQ.zero) THEN
           !$omp parallel do private(I,IY)
           DO i=1,leny
              iy=ky+(i-1)*incy
              y(iy)=zero
           ENDDO
        ELSE
           !$omp parallel do private(I,IY)
           !OCL NOVREC(Y)
           !CDIR NODEP(Y)
           DO i=1,leny
              iy=ky+(i-1)*incy
              y(iy)=beta*y(iy)
           ENDDO
        ENDIF
     ENDIF
  ENDIF
  IF (alpha.EQ.zero) RETURN
  IF (lsame(trans,'N')) THEN
     ! ==--------------------------------------------------------------==
     ! ==  Form y := alpha*A*x + y                                     ==
     ! ==--------------------------------------------------------------==
     ! ..Fill list vector IJX
     ! ..IJX(1:N) -> {KX,KX+INCX,KX+2*INCX,...,KX+(N-1)*INCX}
     DO j=1,n
        ijx(j)=kx+(j-1)*incx
     ENDDO
     IF (incy.EQ.1) THEN
        DO j=1,n
           jx=ijx(j)
           temp=alpha*x(jx)
           !$omp parallel do private(I)
           DO i=1,m
              y(i)=y(i)+temp*a(i,j)
           ENDDO
        ENDDO
     ELSE
        ! ..Fill list vector IIY
        ! ..IIY(1:M) -> {KY,KY+INCY,KY+2*INCY,...,KY+(M-1)*INCY}
        DO j=1,m
           iiy(j)=ky+(j-1)*incy
        ENDDO
        !OCL NOVREC(Y)
        DO j=1,n
           jx=ijx(j)
           temp=alpha*x(jx)
           !$omp parallel do private(I,IY)
           !CDIR NODEP(Y)
           DO i=1,m
              iy=iiy(i)
              y(iy)=y(iy)+temp*a(i,j)
           ENDDO
        ENDDO
     ENDIF
  ELSE
     ! ==--------------------------------------------------------------==
     ! ==  Form y := alpha*A'*x + y or y := alpha*conjg( A' )*x + y    ==
     ! ==--------------------------------------------------------------==
     ! ..Fill list vector IJY
     ! ..IJY(1:N) -> {KY,KY+INCY,KY+2*INCY,...,KY+(N-1)*INCY}
     DO j=1,n
        ijy(j)=ky+(j-1)*incy
     ENDDO
     IF (incx.EQ.1) THEN
        IF (noconj) THEN
           !$omp parallel do private(J,I,TEMP,JY)
           !OCL NOVREC(Y)
           !CDIR NODEP(Y)
           DO j=1,n
              temp=zero
              DO i=1,m
                 temp=temp+a(i,j)*x(i)
              ENDDO
              jy=ijy(j)
              y(jy)=y(jy)+alpha*temp
           ENDDO
        ELSE
           !$omp parallel do private(J,I,TEMP,JY)
           !OCL NOVREC(Y)
           !CDIR NODEP(Y)
           DO j=1,n
              temp=zero
              DO i=1,m
                 temp=temp+CONJG(a(i,j))*x(i)
              ENDDO
              jy=ijy(j)
              y(jy)=y(jy)+alpha*temp
           ENDDO
        ENDIF
     ELSE
        ! ..Fill list vector IIX
        ! ..IIX(1:M) -> {KX,KX+INCX,KX+2*INCX,...,KX+(M-1)*INCX}
        DO j=1,m
           iix(j)=kx+(j-1)*incx
        ENDDO
        IF (noconj) THEN
           !$omp parallel do private(J,I,IX,JY,TEMP)
           !OCL NOVREC(Y)
           !CDIR NODEP(Y)
           DO j=1,n
              temp=zero
              DO i=1,m
                 ix=iix(i)
                 temp=temp+a(i,j)*x(ix)
              ENDDO
              jy=ijy(j)
              y(jy)=y(jy)+alpha*temp
           ENDDO
        ELSE
           !$omp parallel do private(J,I,IX,JY,TEMP)
           !OCL NOVREC(Y)
           !CDIR NODEP(Y)
           DO j=1,n
              temp=zero
              DO i=1,m
                 ix=iix(i)
                 temp=temp+CONJG(a(i,j))*x(ix)
              ENDDO
              jy=ijy(j)
              y(jy)=y(jy)+alpha*temp
           ENDDO
        ENDIF
     ENDIF
  ENDIF
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE zgemv
! ==================================================================

! ==================================================================
SUBROUTINE dsymm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
  ! ==--------------------------------------------------------------==
  ! == Purpose                                                      ==
  ! == =======                                                      ==
  ! ==                                                              ==
  ! == DSYMM  performs one of the matrix-matrix operations          ==
  ! ==                                                              ==
  ! == C := alpha*A*B + beta*C,                                     ==
  ! ==                                                              ==
  ! ==   or                                                         ==
  ! ==                                                              ==
  ! == C := alpha*B*A + beta*C,                                     ==
  ! ==                                                              ==
  ! == where alpha and beta are scalars, A is a symmetric matrix    ==
  ! == and  B and C are  m by n matrices.                           ==
  ! ==                                                              ==
  ! == Parameters                                                   ==
  ! == ==========                                                   ==
  ! ==                                                              ==
  ! == SIDE   - CHARACTER*1.                                        ==
  ! ==   On entry, SIDE specifies whether the symmetric matrix A    ==
  ! ==   appears on the left or right in the operation as follows:  ==
  ! ==                                                              ==
  ! ==   SIDE = 'L' or 'l'   C := alpha*A*B + beta*C,               ==
  ! ==                                                              ==
  ! ==   SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,               ==
  ! ==                                                              ==
  ! ==   Unchanged on exit.                                         ==
  ! ==                                                              ==
  ! == UPLO   - CHARACTER*1.                                        ==
  ! ==   On  entry, UPLO specifies whether the upper or lower       ==
  ! ==   triangular part of the symmetric matrix A is to be         ==
  ! ==   referenced as follows:                                     ==
  ! ==                                                              ==
  ! ==   UPLO = 'U' or 'u' Only the upper triangular part of the    ==
  ! ==                     symmetric matrix is to be referenced.    ==
  ! ==                                                              ==
  ! ==   UPLO = 'L' or 'l'   Only the lower triangular part of the  ==
  ! ==                       symmetric matrix is to be referenced.  ==
  ! ==                                                              ==
  ! ==   Unchanged on exit.                                         ==
  ! ==                                                              ==
  ! == M      - INTEGER.                                            ==
  ! ==   On entry, M specifies the number of rows of the matrix C.  ==
  ! ==   M  must be at least zero.                                  ==
  ! ==   Unchanged on exit.                                         ==
  ! ==                                                              ==
  ! == N      - INTEGER.                                            ==
  ! ==   On entry, N specifies the number of columns of             ==
  ! ==   the matrix C. N must be at least zero.                     ==
  ! ==   Unchanged on exit.                                         ==
  ! ==                                                              ==
  ! == ALPHA  - DOUBLE PRECISION.                                   ==
  ! ==   On entry, ALPHA specifies the scalar alpha.                ==
  ! ==   Unchanged on exit.                                         ==
  ! ==                                                              ==
  ! == A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ),    ==
  ! ==   where ka is m when SIDE = 'L' or 'l' and is n otherwise.   ==
  ! ==   Before entry with SIDE = 'L' or 'l', the m by m part of    ==
  ! ==   the array A must contain the symmetric matrix, such that   ==
  ! ==   when UPLO = 'U' or 'u', the leading m by m upper           ==
  ! ==   triangular part of the array A must contain the upper      ==
  ! ==   triangular part of the symmetric matrix and the strictly   ==  
  ! ==   lower triangular part of  A  is not referenced, and        ==
  ! ==   when UPLO = 'L' or 'l', the leading m by m lower           ==
  ! ==   triangular part of the array A must contain the lower      ==
  ! ==   triangular part of the symmetric matrix and the strictly   ==
  ! ==   upper triangular part of A is not referenced.              ==
  ! ==   Before entry with SIDE = 'R' or 'r', the n by n part of    ==
  ! ==   the array A must contain the symmetric matrix, such that   ==
  ! ==   when UPLO = 'U' or 'u', the leading n by n upper           ==
  ! ==   triangular part of the array A must contain the upper      ==
  ! ==   triangular part of the symmetric matrix and the strictly   ==
  ! ==   lower triangular part of A is not referenced, and when     ==
  ! ==   UPLO = 'L' or 'l', the leading n by n lower triangular     ==
  ! ==   part of the array A must contain the lower triangular part ==
  ! ==   of the symmetric matrix and the strictly upper triangular  ==
  ! ==   part of A is not referenced.                               ==
  ! ==   Unchanged on exit.                                         ==
  ! ==                                                              ==
  ! == LDA    - INTEGER.                                            ==
  ! ==   On entry, LDA specifies the first dimension of A as        ==
  ! ==   declared in the calling (sub) program.                     == 
  ! ==   When SIDE = 'L' or 'l' then                                ==
  ! ==   LDA must be at least max( 1, m ), otherwise LDA must be at ==
  ! ==   least max( 1, n ).                                         ==
  ! ==   Unchanged on exit.                                         ==
  ! ==                                                              ==
  ! == B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).     ==
  ! ==   Before entry, the leading m by n part of the array B must  ==
  ! ==   contain the matrix B.                                      ==
  ! ==   Unchanged on exit.                                         ==
  ! ==                                                              ==
  ! == LDB    - INTEGER.                                            ==
  ! ==   On entry, LDB specifies the first dimension of B as        ==
  ! ==   declared in the calling (sub) program. LDB must be         ==
  ! ==   at least max( 1, m ).                                      ==
  ! ==   Unchanged on exit.                                         ==
  ! ==                                                              ==
  ! == BETA   - DOUBLE PRECISION.                                   ==
  ! ==   On entry, BETA specifies the scalar beta. When BETA is     ==
  ! ==   supplied as zero then C need not be set on input.          ==
  ! ==   Unchanged on exit.                                         ==
  ! ==                                                              ==
  ! == C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).     ==
  ! ==   Before entry, the leading m by n part of the array C must  ==
  ! ==   contain the matrix C, except when beta is zero, in which   ==
  ! ==   case C need not be set on entry.                           ==
  ! ==   On exit, the array C is overwritten by the m by n updated  ==
  ! ==   matrix.                                                    ==
  ! ==                                                              ==
  ! == LDC    - INTEGER.                                            ==
  ! ==   On entry, LDC specifies the first dimension of C as        ==
  ! ==   declared in the calling (sub) program. LDC must be         ==
  ! ==   at least max( 1, m ).                                      ==
  ! ==   Unchanged on exit.                                         ==
  ! ==                                                              ==
  ! ==                                                              ==
  ! == Level 3 Blas routine.                                        ==
  ! ==                                                              ==
  ! == -- Written on 8-February-1989.                               ==
  ! ==    Jack Dongarra, Argonne National Laboratory.               ==
  ! ==    Iain Duff, AERE Harwell.                                  ==
  ! ==    Jeremy Du Croz, Numerical Algorithms Group Ltd.           ==
  ! ==    Sven Hammarling, Numerical Algorithms Group Ltd.          ==
  ! ==--------------------------------------------------------------==
  ! Scalar Arguments
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  CHARACTER(len=1) :: side,uplo
  INTEGER :: m,n,lda,ldb,ldc
  DOUBLE PRECISION alpha,beta
  ! Array Arguments
  DOUBLE PRECISION a(lda,*),b(ldb,*),c(ldc,*)
  ! External Functions
  LOGICAL :: lsame
  EXTERNAL         lsame
  ! External Subroutines
  EXTERNAL         xerbla
  ! Intrinsic Functions
  INTRINSIC        max
  ! Local Scalars
  LOGICAL :: upper
  INTEGER :: i,info,j,k,nrowa
  DOUBLE PRECISION temp1,temp2
  ! Parameters
  DOUBLE PRECISION one,zero
  PARAMETER       (one=1.0e+0_real_8,zero=0.0e+0_real_8)
  ! ==--------------------------------------------------------------==
  ! 
  ! Set NROWA as the number of rows of A.
  ! 
  IF (lsame(side,'L')) THEN
     nrowa = m
  ELSE
     nrowa = n
  ENDIF
  upper = lsame(uplo,'U')
  ! 
  ! Test the input parameters.
  ! 
  info = 0
  IF ((.NOT.lsame(side,'L')).AND.&
       (.NOT.LSAME(SIDE,'R'))) THEN
     info = 1
  ELSE IF((.NOT.upper).AND.&
       (.NOT.LSAME(UPLO,'L'))) THEN
     info = 2
  ELSE IF (m.LT.0) THEN
     info = 3
  ELSE IF (n.LT.0) THEN
     info = 4
  ELSE IF (lda.LT.MAX(1,nrowa)) THEN
     info = 7
  ELSE IF (ldb.LT.MAX(1,m)) THEN
     info = 9
  ELSE IF (ldc.LT.MAX(1,m)) THEN
     info = 12
  ENDIF
  IF (info.NE.0) THEN
     CALL xerbla('DSYMM ',info)
     RETURN
  ENDIF
  ! 
  ! Quick return if possible.
  ! 
  IF ((m.EQ.0).OR.(n.EQ.0).OR.&
       ((ALPHA.EQ.ZERO).AND.(BETA.EQ.ONE)))&
       RETURN
  ! 
  ! And when alpha.eq.zero.
  ! 
  IF (alpha.EQ.zero) THEN
     IF (beta.EQ.zero) THEN
        !$omp parallel do private(J,I)
        DO j=1,n
           DO i=1,m
              c(i,j)=zero
           ENDDO
        ENDDO
     ELSE
        !$omp parallel do private(J,I)
        DO j=1,n
           DO i=1,m
              c(i,j)=beta*c(i,j)
           ENDDO
        ENDDO
     ENDIF
     RETURN
  ENDIF
  ! 
  ! Start the operations.
  ! 
  IF (lsame(side,'L')) THEN
     ! 
     ! Form  C := alpha*A*B + beta*C.
     ! 
     IF (upper) THEN
        !$omp parallel do private(J,I,K,TEMP1,TEMP2)
        DO j=1,n
           DO i=1,m
              temp1=alpha*b(i,j)
              temp2=zero
              DO k=1,i-1
                 c(k,j)=c(k,j)+temp1*a(k,i)
                 temp2=temp2+b(k,j)*a(k,i)
              ENDDO
              IF (beta.EQ.zero) THEN
                 c(i,j)=temp1*a(i,i)+alpha*temp2
              ELSE
                 c(i,j)=beta*c(i,j)+temp1*a(i,i)+alpha*temp2
              ENDIF
           ENDDO
        ENDDO
     ELSE
        !$omp parallel do private(J,I,K,TEMP1,TEMP2)
        DO j=1,n
           DO i=m,1,-1
              temp1=alpha*b(i,j)
              temp2=zero
              DO k=i+1,m
                 c(k,j)=c(k,j)+temp1*a(k,i)
                 temp2=temp2+b(k,j)*a(k,i)
              ENDDO
              IF (beta.EQ.zero) THEN
                 c(i,j)=temp1*a(i,i)+alpha*temp2
              ELSE
                 c(i,j)=beta*c(i,j)+temp1*a(i,i)+alpha*temp2
              ENDIF
           ENDDO
        ENDDO
     ENDIF
  ELSE
     ! 
     ! Form  C := alpha*B*A + beta*C.
     ! 
     !$omp parallel do private(J,I,K,TEMP1,TEMP2)
     DO j=1,n
        temp1=alpha*a(j,j)
        IF (beta.EQ.zero) THEN
           DO i=1,m
              c(i,j)=temp1*b(i,j)
           ENDDO
        ELSE
           DO i=1,m
              c(i,j)=beta*c(i,j)+temp1*b(i,j)
           ENDDO
        ENDIF
        DO k=1,j-1
           IF (upper) THEN
              temp1=alpha*a(k,j)
           ELSE
              temp1=alpha*a(j,k)
           ENDIF
           DO i=1,m
              c(i,j)=c(i,j)+temp1*b(i,k)
           ENDDO
        ENDDO
        DO k=j+1,n
           IF (upper) THEN
              temp1=alpha*a(j,k)
           ELSE
              temp1=alpha*a(k,j)
           ENDIF
           DO i=1,m
              c(i,j)=c(i,j)+temp1*b(i,k)
           ENDDO
        ENDDO
     ENDDO
  ENDIF
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE dsymm
! ==================================================================

! ==================================================================
SUBROUTINE daxpy(n,da,dx,incx,dy,incy)
  ! ==--------------------------------------------------------------==
  ! constant times a vector plus a vector.
  ! uses unrolled loops for increments equal to one.
  ! jack dongarra, linpack, 3/11/78.
  ! modified 12/3/93, array(1) declarations changed to array(*)
  ! 
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  DOUBLE PRECISION dx(*),dy(*),da
  INTEGER :: i,incx,incy,ix,iy,n
  ! 
  IF (n.LE.0) RETURN
  IF (da .EQ. 0.0_real_8) RETURN
  IF (incx.EQ.1.AND.incy.EQ.1) THEN
     ! 
     ! code for both increments equal to 1
     ! 
     !$omp parallel do private(i)
     DO i=1,n
        dy(i) = dy(i) + da*dx(i)
     ENDDO
     RETURN
  ELSE
     ! 
     ! code for unequal increments or equal increments
     ! not equal to 1
     ! 
     ix = 1
     iy = 1
     IF (incx.LT.0) ix = (-n+1)*incx + 1
     IF (incy.LT.0) iy = (-n+1)*incy + 1
     !$omp parallel do private(i)
     DO i = 1,n
        dy(iy+(i-1)*incy) = dy(iy+(i-1)*incy) + da*dx(ix+(i-1)*incx)
     ENDDO
     RETURN
  ENDIF
  ! ==--------------------------------------------------------------==
END SUBROUTINE daxpy
! ==================================================================

! ==================================================================
SUBROUTINE  dcopy(n,dx,incx,dy,incy)
  ! ==--------------------------------------------------------------==
  ! copies a vector, x, to a vector, y.
  ! uses unrolled loops for increments equal to one.
  ! jack dongarra, linpack, 3/11/78.
  ! modified 12/3/93, array(1) declarations changed to array(*)
  ! 
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  DOUBLE PRECISION dx(*),dy(*)
  INTEGER :: i,incx,incy,ix,iy,n
  ! 
  IF (n.LE.0) RETURN
  IF (incx.EQ.1.AND.incy.EQ.1) THEN
     ! 
     ! code for both increments equal to 1
     ! 
     !$omp parallel do private(i)
     DO i=1,n
        dy(i) = dx(i)
     ENDDO
     RETURN
  ELSE
     ! 
     ! code for unequal increments or equal increments
     ! not equal to 1
     ! 
     ix = 1
     iy = 1
     IF (incx.LT.0) ix = (-n+1)*incx + 1
     IF (incy.LT.0) iy = (-n+1)*incy + 1
     !$omp parallel do private(i)
     DO i = 1,n
        dy(iy+(i-1)*incy) = dx(ix+(i-1)*incx)
     ENDDO
     RETURN
  ENDIF
  ! ==--------------------------------------------------------------==
END SUBROUTINE dcopy
! ==================================================================

! ==================================================================
DOUBLE PRECISION FUNCTION ddot(n,dx,incx,dy,incy)
  ! ==--------------------------------------------------------------==
  ! forms the dot product of two vectors.
  ! uses unrolled loops for increments equal to one.
  ! jack dongarra, linpack, 3/11/78.
  ! modified 12/3/93, array(1) declarations changed to array(*)
  ! 
  DOUBLE PRECISION dx(*),dy(*),dtemp
  INTEGER :: i,incx,incy,ix,iy,n
  ! 
  ddot = 0.0_real_8
  dtemp = 0.0_real_8
  IF (n.LE.0) RETURN
  IF (incx.EQ.1.AND.incy.EQ.1) THEN
     ! 
     ! code for both increments equal to 1
     ! 
     !$omp parallel do private(i) reduction(+:dtemp)
     DO i=1,n
        dtemp = dtemp + dx(i)*dy(i)
     ENDDO
     ddot = dtemp
     RETURN
  ELSE
     ! 
     ! code for unequal increments or equal increments
     ! not equal to 1
     ! 
     ix = 1
     iy = 1
     IF (incx.LT.0) ix = (-n+1)*incx + 1
     IF (incy.LT.0) iy = (-n+1)*incy + 1
     !$omp parallel do private(i) reduction(+:dtemp)
     DO i = 1,n
        dtemp = dtemp + dx(ix+(i-1)*incx)*dy(iy+(i-1)*incy)
     ENDDO
     ddot = dtemp
     RETURN
  ENDIF
  ! ==--------------------------------------------------------------==
END FUNCTION ddot
! ==================================================================

! ==================================================================
SUBROUTINE  dscal(n,da,dx,incx)
  ! ==--------------------------------------------------------------==
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
  IF ( n.LE.0 .OR. incx.LE.0 ) RETURN
  IF (incx.EQ.1) THEN
     ! 
     ! code for increment equal to 1
     ! 
     !$omp parallel do private(i)
     DO i=1,n
        dx(i) = da*dx(i)
     ENDDO
     RETURN
  ELSE
     ! 
     ! code for increment not equal to 1
     ! 
     nincx = n*incx
     !$omp parallel do private(i)
     DO i = 1,nincx,incx
        dx(i) = da*dx(i)
     ENDDO
     RETURN
  ENDIF
  ! ==--------------------------------------------------------------==
END SUBROUTINE dscal
! ==================================================================
