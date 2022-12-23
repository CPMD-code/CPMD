C     ==================================================================
      SUBROUTINE mm_invert(V,LV,N,W)
      use kinds, only: real_4, real_8, int_1, int_2, int_4, int_8
      use error_handling, only: stopgm
      use timer, only: tiset, tihalt
      IMPLICIT NONE
C     ==--------------------------------------------------------------==
C      ________________________________________________________
C     |                                                        |
C     |     INVERT A GENERAL MATRIX WITH COMPLETE PIVOTING     |
C     |                                                        |
C     |    INPUT:                                              |
C     |                                                        |
C     |         V     --ARRAY CONTAINING MATRIX                |
C     |                                                        |
C     |         LV    --LEADING (ROW) DIMENSION OF ARRAY V     |
C     |                                                        |
C     |         N     --DIMENSION OF MATRIX STORED IN ARRAY V  |
C     |                                                        |
C     |         W     --WORK ARRAY WITH AT LEAST 2N ELEMENTS   |
C     |                                                        |
C     |    OUTPUT:                                             |
C     |                                                        |
C     |         V     --INVERSE                                |
C     |                                                        |
C     |    BUILTIN FUNCTIONS: ABS                              |
C     |________________________________________________________|
C
      INTEGER W(*),LV
      REAL*8 V(LV,*),S,T
      INTEGER H,I,J,K,L,M,N,O,P,Q,ISUB
      CALL TISET(' MM_INVERT',ISUB)
      IF ( N .EQ. 1 ) GOTO 120
      O = N + 1
      L = 0
      M = 1
10    IF ( L .EQ. N ) GOTO 90
      K = L
      L = M
      M = M + 1
C     ---------------------------------------
C     |*** FIND PIVOT AND START ROW SWAP ***|
C     ---------------------------------------
      P = L
      Q = L
      S = DABS(V(L,L))
      DO  H = L,N
        DO  I = L,N
          T = DABS(V(I,H))
          IF ( T .LE. S ) GOTO 20
          P = I
          Q = H
          S = T
        ENDDO
      ENDDO
20    CONTINUE
      W(N+L) = P
      W(O-L) = Q
      DO I = 1,N
        T = V(I,L)
        V(I,L) = V(I,Q)
        V(I,Q) = T
      ENDDO
      S = V(P,L)
      V(P,L) = V(L,L)
      IF ( S .EQ. 0.D0 ) GOTO 130
C     -----------------------------
C     |*** COMPUTE MULTIPLIERS ***|
C     -----------------------------
      V(L,L) = -1.D0
      S = 1.D0/S
!$OMP parallel do private(I) shared(S)
#ifdef __SR11000
*poption parallel, tlocal(I)
#endif
      DO I = 1,N
        V(I,L) = -S*V(I,L)
      ENDDO
      J = L
50    J = J + 1
      IF ( J .GT. N ) J = 1
      IF ( J .EQ. L ) GOTO 10
      T = V(P,J)
      V(P,J) = V(L,J)
      V(L,J) = T
      IF ( T .EQ. 0.D0 ) GOTO 50
C     ------------------------------
C     |*** ELIMINATE BY COLUMNS ***|
C     ------------------------------
      IF ( K .EQ. 0 ) GOTO 70
      DO I = 1,K
        V(I,J) = V(I,J) + T*V(I,L)
      ENDDO
70    V(L,J) = S*T
      IF ( M .GT. N ) GOTO 50
      DO I = M,N
        V(I,J) = V(I,J) + T*V(I,L)
      ENDDO
      GOTO 50
C     -----------------------
C     |*** PIVOT COLUMNS ***|
C     -----------------------
90    L = W(K+N)
      DO I = 1,N
        T = V(I,L)
        V(I,L) = V(I,K)
        V(I,K) = T
      ENDDO
      K = K - 1
      IF ( K .GT. 0 ) GOTO 90
C     --------------------
C     |*** PIVOT ROWS ***|
C     --------------------
      DO J = 1,N
        DO I = 2,N
           P = W(I)
           H = O - I
           T = V(P,J)
           V(P,J) = V(H,J)
           V(H,J) = T
        ENDDO
      ENDDO
      CALL TIHALT(' MM_INVERT',ISUB)
      RETURN
120   IF ( V(1,1) .EQ. 0.D0 ) GOTO 130
      V(1,1) = 1.d0/V(1,1)
      CALL TIHALT(' MM_INVERT',ISUB)
      RETURN
130   WRITE(6,*) 'MATRIX HAS NO INVERSE',P,L
      call stopgm('MM_INVERT','DEGENERACY PROBLEM',
     & __LINE__,"__file__")
      END
