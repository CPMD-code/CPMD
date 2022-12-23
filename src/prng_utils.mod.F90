

! a minimal implementation of L'Ecuyer MRG32k3a
! with wrappers to simplify calling from around CPMD
MODULE prng_utils
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: parai
  USE prng,                            ONLY: prng_com
  USE system,                          ONLY: cnti

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: prnginit
  PUBLIC :: repprngu
  PUBLIC :: repprngg
  PUBLIC :: repprngu_vec
  PUBLIC :: repprngu_vec_cmplx, repprngu_mat_cmplx
  PUBLIC :: prngskip
  PUBLIC :: prngparaskip

CONTAINS


  ! /* Does a matrix-vector multiply modulo m. used to skip forward */
  SUBROUTINE mvmultm( a,s,v,m)
    REAL(real_8)                             :: a(3,3), s(3), v(3), m

    REAL(real_8), PARAMETER :: two17 = 131072.0_real_8                 , &
      two53 = 9007199254740992.0_real_8       

    INTEGER                                  :: i, j
    REAL(real_8)                             :: a2, ad217, nx, x(3)

    DO i=1,3
       x(i)=0
       DO j=1,3
          a2=a(i,j)
          nx=x(i)+a2*s(j)
          IF ((nx .GE. two53) .OR. (nx .LE. -two53)) THEN
             ad217 = INT(a2 / two17)
             a2 = a2 - ad217 * two17
             nx = ad217 * s(j)
             nx = (nx - INT(nx / m) * m) * two17 + a2 * s(j) + x(i)
          ENDIF
          nx = nx - INT(nx / m)*m
          IF (nx.LT.0.0_real_8) THEN
             x(i)=nx+m
          ELSE
             x(i)=nx
          ENDIF
       ENDDO
    ENDDO
    v=x
  END SUBROUTINE mvmultm

  ! extracts a uniform prn within (0,1]
  FUNCTION prngunif(sprng) RESULT(prngunif_r)
    REAL(real_8)                             :: sprng(2,3), prngunif_r

    LOGICAL, PARAMETER                       :: anti = .FALSE.  
    REAL(real_8), PARAMETER :: a12 = 1403580.0_real_8                , &
      a13n = 810728.0_real_8                , &
      a21 = 527612.0_real_8                , &
      a23n = 1370589.0_real_8                , &
      m1 = 4294967087.0_real_8               , &
      m2 = 4294944443.0_real_8               

    REAL(real_8)                             :: p1, p2, u

! advances stream by one

    p1=a12*sprng(1,2)-a13n*sprng(1,1)
    p1=p1-INT(p1/m1)*m1
    IF (p1.LT.0.0_real_8) p1=p1+m1

    sprng(1,1)=sprng(1,2)
    sprng(1,2)=sprng(1,3)
    sprng(1,3)=p1

    p2=a21*sprng(2,3)-a23n*sprng(2,1)
    p2=p2-INT(p2/m2)*m2
    IF (p2.LT.0.0_real_8) p2=p2+m2

    sprng(2,1)=sprng(2,2)
    sprng(2,2)=sprng(2,3)
    sprng(2,3)=p2

    ! computes a PRN out of the status
    IF (p1.GT.p2) THEN
       u=(p1-p2)*2.328306549295727688e-10_real_8
    ELSE
       u=(p1-p2+m1)*2.328306549295727688e-10_real_8
    ENDIF
    IF (anti) THEN
       prngunif_r=1-u
    ELSE
       prngunif_r=u
    ENDIF

  END FUNCTION prngunif

  ! Initializes the status of the PRNG using a simple linear congruential
  ! generator
  SUBROUTINE prngseed(seed,sprng)
    INTEGER                                  :: seed
    REAL(real_8)                             :: sprng(2,3)

    REAL(real_8), PARAMETER :: m1 = 4294967087.0_real_8               , &
      m2 = 4294944443.0_real_8               

    INTEGER                                  :: i, last

    last = seed
    DO i=1,3
       last= IEOR(last,ISHFT(last,-30))
       last=last*1812433253 + i
       sprng(1,i)=last
       IF (sprng(1,i).LT.0) sprng(1,i)=-sprng(1,i)
       DO WHILE (sprng(1,i).GT.m1)
          sprng(1,i)=sprng(1,i)-m1
       ENDDO
    ENDDO
    DO i=1,3
       last= IEOR(last,ISHFT(last,-30))
       last=last*1812433253 + i
       sprng(2,i)=last
       IF (sprng(2,i).LT.0) sprng(2,i)=-sprng(2,i)
       DO WHILE (sprng(2,i).GT.m2)
          sprng(2,i)=sprng(2,i)-m2
       ENDDO
    ENDDO
  END SUBROUTINE prngseed

  ! Fast forward across the stream, to get an 'independent' sequence
  SUBROUTINE prngskip(sprng,nsprng)
    REAL(real_8)                             :: sprng(2,3), nsprng(2,3)

    REAL(real_8), DIMENSION(3, 3), PARAMETER :: shift1 = RESHAPE((/&
      82758667.0_real_8, 3672831523.0_real_8, 3672091415.0_real_8,&
      1871391091.0_real_8,   69195019.0_real_8, 3528743235.0_real_8,&
      4127413238.0_real_8, 1871391091.0_real_8,   69195019.0_real_8/),(/3,3/))&
      , shift2 = RESHAPE((/1511326704.0_real_8, 4292754251.0_real_8, &
      3859662829.0_real_8,3759209742.0_real_8, 1511326704.0_real_8, &
      4292754251.0_real_8,1610795712.0_real_8, 3889917532.0_real_8, &
      3708466080.0_real_8/),(/3,3/))
    REAL(real_8), PARAMETER :: m1 = 4294967087.0_real_8, &
      m2 = 4294944443.0_real_8               

    REAL(real_8)                             :: tmp(2,3)

    tmp=sprng
    CALL mvmultm(shift1,tmp(1,:),tmp(1,:),m1)
    CALL mvmultm(shift2,tmp(2,:),tmp(2,:),m2)
    nsprng=tmp
  END SUBROUTINE prngskip

  ! here we use a ratio-of-uniforms rather than box-muller, because 
  ! it is as fast (or faster) and we avoid storing buffers (which makes
  ! restart more compact)
  REAL(real_8) FUNCTION prnggauss(sprng)
    REAL(real_8) :: sprng(2,3)

    REAL(real_8) :: u,v,x,y,q
    DO
       u = prngunif(sprng)
       v = prngunif(sprng)
       v = 1.7156_real_8 * (v-0.5_real_8)

       x = u - 0.449871_real_8
       y = ABS(v)+0.386595_real_8
       q = x*x + y*(0.19600_real_8*y-0.25472_real_8*x)
       IF (q .LT. 0.27597_real_8)  EXIT
       IF (q .GT. 0.27846_real_8) CYCLE
       IF (v*v.LT.-4.0_real_8*LOG(u)*u*u) EXIT
    ENDDO
    prnggauss=v/u
  END FUNCTION prnggauss

  SUBROUTINE prnginit()
    ! TEST OF THE IMPLEMENTATION : the sum must be ~ 5001090.95 (see
    ! L'Ecuyer, Op.Res. 47, 159 (1999))
    ! REAL(8) ss
    ! REPSEED=12345
    ! ss=0.0_real_8
    ! do i=1,10000000
    ! ss=ss+REPPRNGU()
    ! enddo
    ! write(6,*) "PRNG SUM IS ", ss

    CALL prngseed(cnti%iprng, prng_com%repseed)
    CALL prngskip(prng_com%repseed,prng_com%paraseed)
    CALL prngparaskip(prng_com%paraseed,prng_com%paraseed)
  END SUBROUTINE prnginit

  SUBROUTINE prngparaskip(sprng,nsprng)
    REAL(real_8)                             :: sprng(2,3), nsprng(2,3)

    INTEGER                                  :: i
    REAL(real_8)                             :: tmp(2,3)

    tmp=sprng
    DO i=0,parai%me-1
       CALL prngskip(tmp,tmp)
    ENDDO
    nsprng=tmp
  END SUBROUTINE prngparaskip

  ! this set of functions generate uniform (*U) and Gaussian (*G) 
  ! pseudo random numbers. 
  ! REP* versions are initialized so as to generate the same stream on all
  ! nodes
  ! PARA* versions generate independent streams on all nodes.
  ! both functions should be called always by all instances of the code: 
  ! in this way, restart should work properly.
  ! NB: RUNS WITH DIFFERENT NUMBER OF PROCESSORS WILL GIVE DIFFERENT 
  ! (but statistically equivalent) RESULTS
  REAL(real_8) FUNCTION repprngu()
    repprngu=prngunif(prng_com%repseed)
  END FUNCTION repprngu
  SUBROUTINE repprngu_vec(n,a)
    INTEGER                                  :: n
    REAL(real_8), DIMENSION(n)               :: a

    INTEGER                                  :: i

    DO i=1,n
       a(i)=repprngu()
    ENDDO
    RETURN
  END SUBROUTINE repprngu_vec

  SUBROUTINE repprngu_vec_cmplx(n,a)
    INTEGER                                  :: n
    COMPLEX(real_8), DIMENSION(n)            :: a

    INTEGER                                  :: i
    REAL(real_8)                             :: re(2)

    DO i=1,n
       CALL repprngu_vec(2,re)
       a(i)=CMPLX(re(1),re(2),kind=real_8)
    ENDDO
  END SUBROUTINE repprngu_vec_cmplx

  SUBROUTINE repprngu_mat_cmplx(m,n,a)
    INTEGER, INTENT(IN)                      :: m, n
    COMPLEX(real_8), DIMENSION(:, :)         :: a

    INTEGER                                  :: i, j
    REAL(real_8)                             :: re(2)

    DO j=1,n!size(a,2)
       DO i=1,m!size(a,1)
          CALL repprngu_vec(2,re)
          a(i,j)=CMPLX(re(1),re(2),kind=real_8)
       ENDDO
    ENDDO
  END SUBROUTINE repprngu_mat_cmplx

  REAL(real_8) FUNCTION repprngg()
    repprngg=prnggauss(prng_com%repseed)
  END FUNCTION repprngg

  REAL(real_8) FUNCTION paraprngu()
    paraprngu=prngunif(prng_com%paraseed)
  END FUNCTION paraprngu

  REAL(real_8) FUNCTION paraprngg()
    paraprngg=prnggauss(prng_com%paraseed)
  END FUNCTION paraprngg

END MODULE prng_utils





