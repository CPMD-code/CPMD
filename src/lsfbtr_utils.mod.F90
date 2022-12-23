MODULE lsfbtr_utils
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: lsfbtr

CONTAINS

#if defined(__PRIMERGY)
  !OCL NOPACKED
#endif
  ! ==================================================================
  SUBROUTINE lsfbtr(f,g,l,rhomin,kapmin,h,n,saved,xa,ta,&
       TRBE,WW,NDIM,DISC)
    ! ==--------------------------------------------------------------==
    ! THIS PROGRAM CALCULATES THE SPHERICAL BESSEL TRANSFORM OF ORDER      
    ! L FOR A FUNCTION WHOSE VALUES ARE GIVEN IN THE ARRAY F FOR 2**N      
    ! VALUES OF R.  THE VALUES OF R ARE GIVEN BY R(I)=EXP(RHOMIN+(I-1)*H)  
    ! WHERE I=1,2,...,NT AND NT=2**N.  THE VALUES OF THE TRANSFORM         
    ! ARE RETURNED IN THE ARRAY G FOR NT K VALUES.  THE VALUES OF K        
    ! ARE EXP(KAPMIN+(I-1)*H), I=1,2,...,NT.                               
    ! 
    ! TWO DIFFERENT SETS OF RESULTS ARE CALCULATED.  ONE SET IS ACCU-      
    ! RATE FOR LARGE VALUES OF K AND THE OTHER IS ACCURATE FOR K CLOSE     
    ! TO ZERO.  THE TWO SETS ARE JOINED AT THE K VALUE FOR WHICH THE       
    ! DIFFERENCE BETWEEN THEM IS LEAST.  FOR L=0 OR L=1 THE SMALL K        
    ! RESULTS ARE CALCULATED USING SIEGMANS METHOD.                        
    ! 
    ! IF THE LOGICAL VARIABLE SAVED IS .FALSE. AUXILIARY DATA IN THE       
    ! ARRAYS TA, TRBE, AND WW ARE CALCULATED ANEW.  THIS MUST BE DONE      
    ! IF THE MESH PARAMETERS N, H, RHOMIN OR KAPMIN HAVE CHANGED FROM      
    ! THE PREVIOUS CALL TO THE SUBROUTINE AND ON THE FIRST CALL TO         
    ! THE SUBROUTINE.  OTHERWISE A SUBSTANTIAL IMPROVEMENT IN SPEED        
    ! CAN BE OBTAINED BY TAKING SAVED=.TRUE..                              
    ! 
    ! THE VARIABLE DISC IS THE DIFFERENCE BETWEEN THE SMALL K AND          
    ! LARGE K VALUES FOR THE RESULT AT THE MATCHING POINT OF MINIMUM       
    ! DISCREPANCY.  IT PROVIDES A ROUGH GUIDE TO THE ABSOLUTE ACC-         
    ! URACY OF THE RESULTS AT THE MAXIMUM VALUE OF THE TRANSFORM.          
    ! IT IS, HOWEVER, A GREAT OVERESTIMATE OF THE ABSOLUTE ERROR AT        
    ! LARGE K VALUES, AND SMALL K VALUES IN THE CASE OF NON-ZERO L.        
    ! 
    ! XA, TA, TRBE AND WW ARE COMPLEX WORKING AREA ARRAYS DIMENSIONED      
    ! AS XA(NDIM), TA(NDIM), TRBE(NDIM,2) AND WW(NDIM) WHERE NDIM          
    ! MUST BE AT LEAST AS LARGE AS 2**N.  AS NOTED, IT MAY BE USEFUL       
    ! TO PRESERVE THE CONTENTS OF TA, TRBE AND WW BUT THE CONTENTS OF      
    ! XA ARE NOT OF VALUE.                                                 
    ! 
    INTEGER                                  :: l
    REAL(real_8)                             :: rhomin, kapmin, h
    INTEGER                                  :: n
    LOGICAL                                  :: saved
    INTEGER                                  :: ndim
    COMPLEX(real_8)                          :: ww(ndim), trbe(ndim,2), &
                                                ta(ndim), xa(ndim)
    REAL(real_8)                             :: g(ndim), f(ndim), disc

    COMPLEX(real_8)                          :: W1, W2, Y1, Y2, ZZM
    INTEGER                                  :: i, ic, ii, ij, ijk, jj, lp, &
                                                na, nh, nhp, nt
    REAL(real_8)                             :: aa, an, bb, cc, cd, cf, cl, &
                                                cm, d1, d2, dt, phi, pi, rr, &
                                                s, t, te, xx, yy

    nt=2**n
    IF (nt.GT.ndim) go to 131
    nh=nt/2
    nhp=nh+1
    pi=4._real_8*ATAN(1._real_8)
    an=nt
    dt=2._real_8*pi/(h*an)
    cm=SQRT(2._real_8*pi)/an
    cl=0.25_real_8*h/an
    IF (.NOT.saved) go to 100
    ! 
    ! CALCULATE RESULTS ACCURATE AT SMALL K VALUES                         
    ! 
1   IF (l.LT.2) go to 11
    aa=EXP((l+1.5_real_8)*rhomin)
    bb=EXP((l+1.5_real_8)*h)
    DO i=1,nt
       xa(i)=f(i)*aa
       aa=aa*bb
    ENDDO
    CALL nlogn(n,xa,ww,ndim)
    DO i=1,nh
       xa(i)=ta(i)*xa(i)
    ENDDO
    DO 4 jj=1,l
       aa=2*jj-0.5_real_8
       t=0._real_8
       zzm=1._real_8/CMPLX(aa,t,kind=real_8)
       DO i=1,nh
          xa(i)=zzm*xa(i)
          t=t+dt
       ENDDO
4   CONTINUE
    DO i=nhp,nt
       xa(i)=0._real_8
    ENDDO
    CALL nlogn(n,xa,ww,ndim)
    aa=EXP((l-1.5_real_8)*kapmin)*cm
    bb=EXP((l-1.5_real_8)*h)
    DO i=1,nt
       g(i)=REAL(aa*xa(i))
       aa=aa*bb
    ENDDO
    go to 21
    ! 
    ! FOR L=0 OR 1 CALCULATE RESULTS ACCURATE AT SMALL K VALUES USING      
    ! SIEGMANS METHOD                                                      
    ! 
11  lp=l+1
    DO i=nh,nt
       xa(i)=0._real_8
    ENDDO
    aa=EXP(3._real_8*rhomin)
    bb=EXP(3._real_8*h)
    DO i=1,nh
       xx=aa*f(2*i-1)
       aa=aa*bb
       yy=aa*f(2*i)
       xa(i)=CMPLX(xx,yy,kind=real_8)
       aa=aa*bb
    ENDDO
    CALL nlogn(n,xa,ww,ndim)
    y1=trbe(1,lp)*xa(1)
    y2=trbe(1,lp)*CONJG(xa(1))
    xa(1)=2._real_8*(y1+y2+CONJG(y2-y1))
    xa(nhp)=4._real_8*CONJG(xa(nhp)*trbe(nhp,lp))
    DO i=2,nh
       ic=nt-i+2
       y1=xa(i)
       y2=CONJG(xa(ic))
       w1=y1+y2
       w2=ww(i+nh)*(y1-y2)
       y1=(w1-w2)*trbe(i,lp)
       y2=(w1+w2)*CONJG(trbe(ic,lp))
       w1=y1+y2
       w2=ww(i+nh)*(y1-y2)
       xa(i)=w1+w2
       xa(ic)=CONJG(w1-w2)
    ENDDO
    CALL nlogn(n,xa,ww,ndim)
    DO i=1,nh
       g(2*i-1)=cl*REAL(xa(i))
       g(2*i)=cl*AIMAG(xa(i))
    ENDDO
    ! 
    ! CALCULATE RESULTS ACCURATE AT LARGE K VALUES                         
    ! 
21  aa=EXP(1.5_real_8*rhomin)
    bb=EXP(1.5_real_8*h)
    DO i=1,nt
       xa(i)=aa*f(i)
       aa=aa*bb
    ENDDO
    CALL nlogn(n,xa,ww,ndim)
    ij=MOD(l,2)
    ijk=MOD(l,4)
    DO i=1,nh
       y1=xa(i)*ta(i+ij*nh)
       IF (ijk.GT.1) y1=-y1
       xa(i)=y1
    ENDDO
    IF (l.EQ.0) go to 24
    DO 25 jj=1,l
       aa=2*jj-l-0.5_real_8
       bb=jj-0.5_real_8
       t=0._real_8
       DO i=1,nh
          xa(i)=xa(i)*CMPLX(bb,-t,kind=real_8)/CMPLX(aa,t,kind=real_8)
          t=t+dt
       ENDDO
25  CONTINUE
24  DO i=nhp,nt
       xa(i)=0._real_8
    ENDDO
    CALL nlogn(n,xa,ww,ndim)
    aa=EXP(-1.5_real_8*kapmin)*cm
    bb=EXP(-1.5_real_8*h)
    DO i=1,nt
       xa(i)=aa*xa(i)
       aa=aa*bb
    ENDDO
    ! 
    ! FIND POINT OF MINIMUM DISCREPANCY BETWEEN SMALL K AND LARGE K VALUES 
    ! OF TRANSFORM AND CONSTRUCT BEST APPROXIMATE RESULT.                  
    ! 
    d1=ABS(g(1)-REAL(xa(1)))
    d2=ABS(g(2)-REAL(xa(2)))
    te=d1+d2
    ii=2
    d1=d2
    DO 31 i=3,nt
       d2=ABS(g(i)-REAL(xa(i)))
       aa=d1+d2
       IF (aa.GE.te) go to 32
       ii=i
       te=aa
32     d1=d2
31  CONTINUE
    DO i=ii,nt
       g(i)=REAL(xa(i))
    ENDDO
    disc=te
    RETURN
    ! 
    ! INITIALIZE AUXILIARY DATA FOR TRANSFORM                              
    ! 
    ! THE NEXT INSTRUCTIONS INITIALIZE THE ARRAY TA USED IN CALCULATING    
    ! THE TRANSFORM AT LARGE K VALUES AND SMALL K VALUES FOR L GREATER     
    ! THAN 1.                                                              
    ! 
100 CONTINUE
    y1=1._real_8
    aa=(rhomin+kapmin)*dt
    y2=CMPLX(COS(aa),SIN(aa),kind=real_8)
    DO i=1,nh
       t=(i-1)*dt
       s=0._real_8
       rr=SQRT(110.25_real_8+t*t)
       phi=ATAN(t/10.5_real_8)
       DO na=1,10
          s=s+ATAN(t/(na-0.5_real_8))
       ENDDO
       s=s-t*LOG(rr)+t-10._real_8*phi+SIN(phi)/(12._real_8*rr)
       s=s-SIN(3._real_8*phi)/(360._real_8*rr**3)&
            +SIN(5._real_8*PHI)/(1260._real_8*RR**5)&
            -SIN(7._real_8*PHI)/(1680._real_8*RR**7)
       xx=EXP(pi*t)
       xx=ATAN((xx-1._real_8)/(xx+1._real_8))
       cc=s-xx
       ta(i)=y1*CMPLX(COS(cc),SIN(cc),kind=real_8)
       cc=s+xx
       ta(i+nh)=y1*CMPLX(COS(cc),SIN(cc),kind=real_8)
       y1=y1*y2
    ENDDO
    ta(1)=ta(1)/2._real_8
    ta(1+nh)=ta(1+nh)/2._real_8
    ta(nh)=ta(nh)/2._real_8
    ta(nt)=ta(nt)/2._real_8
    ! 
    ! THE NEXT INSTRUCTIONS INITIALIZE THE ARRAY WW USED BY                
    ! THE NLOGN SUBROUTINE.  THE ELEMENTS IN THE SECOND HALF OF THE        
    ! ARRAY WW ARE USED IN THE IMPLEMENTATION OF SIEGMANS METHOD.          
    ! 
    DO i=1,nh
       xx=(i-1)*pi/an
       ww(i+nh)=CMPLX(-SIN(xx),COS(xx),kind=real_8)
       xx=2._real_8*xx
       ww(i)=CMPLX(COS(xx),SIN(xx),kind=real_8)
    ENDDO
    ! 
    ! THE NEXT INSTRUCTIONS INITIALIZE THE ARRAY TRBE USED IN THE          
    ! IMPLEMENTATION OF SIEGMANS METHOD.                                   
    ! 
    DO i=1,nt
       xa(i)=0._real_8
    ENDDO
    cf=EXP(h)
    xx=EXP(rhomin+kapmin)
    DO 113 i=1,nt
       aa=SIN(xx)/xx
       xx=cf*xx
       bb=SIN(xx)/xx
       xa(i)=CMPLX(aa,bb,kind=real_8)
       IF (xx.GT.1.e8_real_8) go to 114
       xx=xx*cf
113 CONTINUE
114 CALL nlogn(n,xa,ww,ndim)
    trbe(1,1)=xa(1)
    trbe(nhp,1)=CONJG(xa(nhp))
    DO i=2,nh
       ic=nt-i+2
       y1=xa(i)
       y2=CONJG(xa(ic))
       w1=y1+y2
       w2=ww(nh+i)*(y1-y2)
       y1=w1-w2
       y2=CONJG(w1+w2)
       trbe(i,1) =0.5_real_8*CONJG(y1)
       trbe(ic,1)=0.5_real_8*CONJG(y2)
    ENDDO
    DO i=1,nt
       xa(i)=0._real_8
    ENDDO
    xx=EXP(rhomin+kapmin)
    DO 121 i=1,nt
       IF (xx.LT.0.1_real_8) go to 122
       aa=(SIN(xx)/xx-COS(xx))/xx
       xx=cf*xx
       bb=(SIN(xx)/xx-COS(xx))/xx
       xa(i)=CMPLX(aa,bb,kind=real_8)
       IF (xx.GT.1.e8_real_8) go to 123
       go to 121
122    cc=xx*xx/2._real_8
       cd=1._real_8-cc/5._real_8+cc*cc/70._real_8-cc*cc*cc/1890._real_8+cc**4/83160._real_8
       aa=xx*cd/3._real_8
       xx=cf*xx
       cc=xx*xx/2._real_8
       cd=1._real_8-cc/5._real_8+cc*cc/70._real_8-cc*cc*cc/1890._real_8+cc**4/83160._real_8
       bb=xx*cd/3._real_8
       xa(i)=CMPLX(aa,bb,kind=real_8)
121    xx=xx*cf
    CONTINUE
123 CALL nlogn(n,xa,ww,ndim)
    trbe(1,2)=xa(1)
    trbe(nhp,2)=CONJG(xa(nhp))
    DO i=2,nh
       ic=nt-i+2
       y1=xa(i)
       y2=CONJG(xa(ic))
       w1=y1+y2
       w2=ww(nh+i)*(y1-y2)
       y1=w1-w2
       y2=CONJG(w1+w2)
       trbe(i,2) =0.5_real_8*CONJG(y1)
       trbe(ic,2)=0.5_real_8*CONJG(y2)
    ENDDO
    go to 1
131 WRITE (6,132) ndim,nt
132 FORMAT(1x,'SUBROUTINE LSFBTR NOT EXECUTED. DIMENSION',&
         ' NDIM TOO SMALL: NDIM =',i10,'  NT = 2**N =',i10)
    RETURN
  END SUBROUTINE lsfbtr
#if defined(__PRIMERGY)
  !OCL NOPACKED
#endif
  ! ==================================================================
  SUBROUTINE nlogn(n,x,ww,ndim)
    ! ==--------------------------------------------------------------==
    ! FAST FOURIER TRANSFORM ROUTINE                                      
    ! 
    INTEGER                                  :: n, ndim
    COMPLEX(real_8)                          :: ww(ndim), x(ndim)

    COMPLEX(real_8)                          :: hold, q, wk
    INTEGER                                  :: i, iblock, ii, istart, j, jh, &
                                                k, l, lbhalf, lblock, lx, &
                                                mm(15), nblock

    DO i=1,n
       mm(i)=2**(n-i)
    ENDDO
    lx=2*mm(1)
    DO 4 l=1,n
       nblock=2**(l-1)
       lblock=lx/nblock
       lbhalf=lblock/2
       k=0
       DO 4 iblock=1,nblock
          wk=ww(k+1)
          istart=lblock*(iblock-1)
          DO 2 i=1,lbhalf
             j=istart+i
             jh=j+lbhalf
             q=x(jh)*wk
             x(jh)=x(j)-q
             x(j)=x(j)+q
2         CONTINUE
          DO 3 i=2,n
             ii=i
             IF (k.LT.mm(i)) go to 4
             k=k-mm(i)
3         CONTINUE
4         k=k+mm(ii)
       CONTINUE
    CONTINUE
    k=0
    DO 7 j=1,lx
       IF (k.LT.j) go to 5
       hold= x(j)
       x(j)=x(k+1)
       x(k+1)= hold
5      DO 6 i=1,n
          ii=i
          IF (k.LT.mm(i)) go to 7
          k=k-mm(i)
6      CONTINUE
7      k=k+mm(ii)
    RETURN
  END SUBROUTINE nlogn

END MODULE lsfbtr_utils
