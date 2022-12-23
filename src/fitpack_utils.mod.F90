#ifdef __SR8000
!option MP(P(0)), LANGLVL(SAVE(0))
#endif


MODULE fitpack_utils
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: curv1
  PUBLIC :: curv2
  PUBLIC :: curvd
  PUBLIC :: curvi

CONTAINS

  ! From inet!cs.utexas.edu!cline Tue Oct 3117:10:31 CST 1989
  ! Received: from mojave.cs.utexas.edu by cs.utexas.edu (5.59/1.44)
  ! id AA29509; Tue, 31 Oct 8917:11:51 CST
  ! Posted-Date: Tue, 31 Oct 8917:10:31 CST
  ! Message-Id: <8910312310.AA04442@mojave.cs.utexas.edu>
  ! Received: by mojave.cs.utexas.edu (14.5/1.4-Client)
  ! id AA04442; Tue, 31 Oct 8917:10:34 cst
  ! Date: Tue, 31 Oct 8917:10:31 CST
  ! From: cline@cs.utexas.edu (Alan Cline)
  ! To: ehg@research.att.com
  ! Subject: New FITPACK Subset for netlib
  ! 
  ! 
  ! This new version of FITPACK distributed by netlib is about 20% of 
  ! the total package in terms of characters, lines of code, and num-
  ! ber of subprograms. However, these 25 subprograms represent about
  ! 95% of usages of the package.  What has been omitted are such ca-
  ! pabilities as:
  ! 1. Automatic tension determination,
  ! 2. Derivatives, arclengths, and enclosed areas for planar 
  ! curves,
  ! 3. Three dimensional curves,
  ! 4. Special surface fitting using equispacing assumptions,
  ! 5. Surface fitting in annular, wedge, polar, toroidal, lunar,
  ! and spherical geometries,
  ! 6. B-splines in tension generation and usage,
  ! 7. General surface fitting in three dimensional space.
  ! 
  ! (The code previously circulated in netlib is less than 10% of the
  ! total  package  and is more than a decade old.  Its usage is dis-
  ! couraged.)
  ! 
  ! Please note:  Two versions of the subroutine snhcsh are included.
  ! Both serve the same purpose:  obtaining approximations to certain
  ! hyperbolic trigonometric-like functions.  The first is less accu-
  ! rate (but more efficient) than the second.  Installers should se- 
  ! lect the one with the precision they desire.
  ! 
  ! Interested parties can obtain the entire package on disk or  tape
  ! from Pleasant  Valley Software, 8603 Altus Cove, Austin TX (USA),
  ! 78759 at a cost of $495 US. A 340 page manual  is  available  for
  ! 30  US  per  copy.  The  package  includes  examples and machine
  ! readable documentation.
  ! 
  SUBROUTINE curv1 (n,x,y,slp1,slpn,islpsw,yp,temp,&
       sigma,ierr)
    INTEGER                                  :: n
    REAL(real_8)                             :: x(n), y(n), slp1, slpn
    INTEGER                                  :: islpsw
    REAL(real_8)                             :: yp(n), temp(n), sigma
    INTEGER                                  :: ierr

    INTEGER                                  :: i, ibak, nm1, np1
    REAL(real_8)                             :: c1, c2, c3, delx1, delx2, &
                                                delxn, delxnm, diag, diag1, &
                                                diag2, dx1, dx2, sdiag1, &
                                                sdiag2, sigmap, slpp1, slppn

! 
! coded by alan kaylor cline
! from fitpack -- january 26, 1987
! a curve and surface fitting package
! a product of pleasant valley software
! 8603 altus cove, austin, texas 78759, usa
! 
! this subroutine determines the parameters necessary to
! compute an interpolatory spline under tension through
! a sequence of functional values. the slopes at the two
! ends of the curve may be specified or omitted.  for actual
! computation of points on the curve it is necessary to call
! the function curv2.
! 
! on input--
! 
! n is the number of values to be interpolated (n.ge.2).
! 
! x is an array of the n increasing abscissae of the
! functional values.
! 
! y is an array of the n ordinates of the values, (i. e.
! y(k) is the functional value corresponding to x(k) ).
! 
! slp1 and slpn contain the desired values for the first
! derivative of the curve at x(1) and x(n), respectively.
! the user may omit values for either or both of these
! parameters and signal this with islpsw.
! 
! islpsw contains a switch indicating which slope data
! should be used and which should be estimated by this
! subroutine,
! = 0 if slp1 and slpn are to be used,
! = 1 if slp1 is to be used but not slpn,
! = 2 if slpn is to be used but not slp1,
! = 3 if both slp1 and slpn are to be estimated
! internally.
! 
! yp is an array of length at least n.
! 
! temp is an array of length at least n which is used for
! scratch storage.
! 
! and
! 
! sigma contains the tension factor. this value indicates
! the curviness desired. if abs(sigma) is nearly zero
! (e.g. .001) the resulting curve is approximately a
! cubic spline. if abs(sigma) is large (e.g. 50.) the
! resulting curve is nearly a polygonal line. if sigma
! equals zero a cubic spline results.  a standard value
! for sigma is approximately 1. in absolute value.
! 
! on output--
! 
! yp contains the values of the second derivative of the
! curve at the given nodes.
! 
! ierr contains an error flag,
! = 0 for normal return,
! = 1 if n is less than 2,
! = 2 if x-values are not strictly increasing.
! 
! and
! 
! n, x, y, slp1, slpn, islpsw and sigma are unaltered.
! 
! this subroutine references package modules ceez, terms,
! and snhcsh.
! 
! -----------------------------------------------------------
! 

    nm1 = n-1
    np1 = n+1
    ierr = 0
    IF (n .LE. 1) go to 8
    IF (x(n) .LE. x(1)) go to 9
    ! 
    ! denormalize tension factor
    ! 
    sigmap = ABS(sigma)*REAL(n-1,kind=real_8)/(x(n)-x(1))
    ! 
    ! approximate end slopes
    ! 
    IF (islpsw .GE. 2) go to 1
    slpp1 = slp1
    go to 2
1   delx1 = x(2)-x(1)
    delx2 = delx1+delx1
    IF (n .GT. 2) delx2 = x(3)-x(1)
    IF (delx1 .LE. 0._real_8 .OR. delx2 .LE. delx1) go to 9
    CALL ceez (delx1,delx2,sigmap,c1,c2,c3,n)
    slpp1 = c1*y(1)+c2*y(2)
    IF (n .GT. 2) slpp1 = slpp1+c3*y(3)
2   IF (islpsw .EQ. 1 .OR. islpsw .EQ. 3) go to 3
    slppn = slpn
    go to 4
3   delxn = x(n)-x(nm1)
    delxnm = delxn+delxn
    IF (n .GT. 2) delxnm = x(n)-x(n-2)
    IF (delxn .LE. 0._real_8 .OR. delxnm .LE. delxn) go to 9
    CALL ceez (-delxn,-delxnm,sigmap,c1,c2,c3,n)
    slppn = c1*y(n)+c2*y(nm1)
    IF (n .GT. 2) slppn = slppn+c3*y(n-2)
    ! 
    ! set up right hand side and tridiagonal system for yp and
    ! perform forward elimination
    ! 
4   delx1 = x(2)-x(1)
    IF (delx1 .LE. 0._real_8) go to 9
    dx1 = (y(2)-y(1))/delx1
    CALL terms (diag1,sdiag1,sigmap,delx1)
    yp(1) = (dx1-slpp1)/diag1
    temp(1) = sdiag1/diag1
    IF (n .EQ. 2) go to 6
    DO 5 i = 2,nm1
       delx2 = x(i+1)-x(i)
       IF (delx2 .LE. 0._real_8) go to 9
       dx2 = (y(i+1)-y(i))/delx2
       CALL terms (diag2,sdiag2,sigmap,delx2)
       diag = diag1+diag2-sdiag1*temp(i-1)
       yp(i) = (dx2-dx1-sdiag1*yp(i-1))/diag
       temp(i) = sdiag2/diag
       dx1 = dx2
       diag1 = diag2
5   sdiag1 = sdiag2
6   diag = diag1-sdiag1*temp(nm1)
    yp(n) = (slppn-dx1-sdiag1*yp(nm1))/diag
    ! 
    ! perform back substitution
    ! 
    DO 7 i = 2,n
       ibak = np1-i
7   yp(ibak) = yp(ibak)-temp(ibak)*yp(ibak+1)
    RETURN
    ! 
    ! too few points
    ! 
8   ierr = 1
    RETURN
    ! 
    ! x-values not strictly increasing
    ! 
9   ierr = 2
    RETURN
  END SUBROUTINE curv1


  SUBROUTINE curvs (n,x,y,d,isw,s,eps,ys,ysp,sigma,temp,&
       ierr)
    INTEGER                                  :: n
    REAL(real_8)                             :: x(n), y(n), d(n)
    INTEGER                                  :: isw
    REAL(real_8)                             :: s, eps, ys(n), ysp(n), sigma, &
                                                temp(n,9)
    INTEGER                                  :: ierr

! 
! coded by alan kaylor cline
! from fitpack -- january 26, 1987
! a curve and surface fitting package
! a product of pleasant valley software
! 8603 altus cove, austin, texas 78759, usa
! 
! this subroutine determines the parameters necessary to
! compute a smoothing spline under tension. for a given
! increasing sequence of abscissae (x(i)), i = 1,..., n and
! associated ordinates (y(i)), i = 1,..., n, the function
! determined minimizes the summation from i = 1 to n-1 of
! the square of the second derivative of f plus sigma
! squared times the difference of the first derivative of f
! and (f(x(i+1))-f(x(i)))/(x(i+1)-x(i)) squared, over all
! functions f with two continuous derivatives such that the
! summation of the square of (f(x(i))-y(i))/d(i) is less
! than or equal to a given constant s, where (d(i)), i = 1,
! ..., n are a given set of observation weights. the
! function determined is a spline under tension with third
! derivative discontinuities at (x(i)), i = 2,..., n-1. for
! actual computation of points on the curve it is necessary
! to call the function curv2. the determination of the curve
! is performed by subroutine curvss, the subroutine curvs
! only decomposes the workspace for curvss.
! 
! on input--
! 
! n is the number of values to be smoothed (n.ge.2).
! 
! x is an array of the n increasing abscissae of the
! values to be smoothed.
! 
! y is an array of the n ordinates of the values to be
! smoothed, (i. e. y(k) is the functional value
! corresponding to x(k) ).
! 
! d is a parameter containing the observation weights.
! this may either be an array of length n or a scalar
! (interpreted as a constant). the value of d
! corresponding to the observation (x(k),y(k)) should
! be an approximation to the standard deviation of error.
! 
! isw contains a switch indicating whether the parameter
! d is to be considered a vector or a scalar,
! = 0 if d is an array of length n,
! = 1 if d is a scalar.
! 
! s contains the value controlling the smoothing. this
! must be non-negative. for s equal to zero, the
! subroutine does interpolation, larger values lead to
! smoother funtions. if parameter d contains standard
! deviation estimates, a reasonable value for s is
! float(n).
! 
! eps contains a tolerance on the relative precision to
! which s is to be interpreted. this must be greater than
! or equal to zero and less than or equal to one. a
! reasonable value for eps is sqrt(2./float(n)).
! 
! ys is an array of length at least n.
! 
! ysp is an array of length at least n.
! 
! sigma contains the tension factor. this value indicates
! the degree to which the first derivative part of the
! smoothing functional is emphasized. if sigma is nearly
! zero (e. g. .001) the resulting curve is approximately a
! cubic spline. if sigma is large (e. g. 50.) the
! resulting curve is nearly a polygonal line. if sigma
! equals zero a cubic spline results. a standard value for
! sigma is approximately 1.
! 
! and
! 
! temp is an array of length at least 9*n which is used
! for scratch storage.
! 
! on output--
! 
! ys contains the smoothed ordinate values.
! 
! ysp contains the values of the second derivative of the
! smoothed curve at the given nodes.
! 
! ierr contains an error flag,
! = 0 for normal return,
! = 1 if n is less than 2,
! = 2 if s is negative,
! = 3 if eps is negative or greater than one,
! = 4 if x-values are not strictly increasing,
! = 5 if a d-value is non-positive.
! 
! and
! 
! n, x, y, d, isw, s, eps, and sigma are unaltered.
! 
! this subroutine references package modules curvss, terms,
! and snhcsh.
! 
! -----------------------------------------------------------
! 
! decompose temp into nine arrays and call curvss
! 

    CALL curvss (n,x,y,d,isw,s,eps,ys,ysp,sigma,temp(1,1),&
         temp(1,2),temp(1,3),temp(1,4),temp(1,5),&
         temp(1,6),temp(1,7),temp(1,8),temp(1,9),&
         ierr)
    RETURN
  END SUBROUTINE curvs

  REAL(real_8) FUNCTION curv2 (t,n,x,y,yp,sigma)
    INTEGER :: n
    REAL(real_8) :: t,x(n),y(n),yp(n),sigma
    REAL(real_8) :: del1, del2, dels, dummy, s1, s2, sigdel, sigmap, ss, sum
    INTEGER :: i, im1 
    ! 
    ! coded by alan kaylor cline
    ! from fitpack -- january 26, 1987
    ! a curve and surface fitting package
    ! a product of pleasant valley software
    ! 8603 altus cove, austin, texas 78759, usa
    ! 
    ! this function interpolates a curve at a given point
    ! using a spline under tension. the subroutine curv1 should
    ! be called earlier to determine certain necessary
    ! parameters.
    ! 
    ! on input--
    ! 
    ! t contains a real value to be mapped onto the interpo-
    ! lating curve.
    ! 
    ! n contains the number of points which were specified to
    ! determine the curve.
    ! 
    ! x and y are arrays containing the abscissae and
    ! ordinates, respectively, of the specified points.
    ! 
    ! yp is an array of second derivative values of the curve
    ! at the nodes.
    ! 
    ! and
    ! 
    ! sigma contains the tension factor (its sign is ignored).
    ! 
    ! the parameters n, x, y, yp, and sigma should be input
    ! unaltered from the output of curv1.
    ! 
    ! on output--
    ! 
    ! curv2 contains the interpolated value.
    ! 
    ! none of the input parameters are altered.
    ! 
    ! this function references package modules intrvl and
    ! snhcsh.
    ! 
    ! -----------------------------------------------------------
    ! 
    ! determine interval
    ! 
    !$omp critical
    im1 = intrvl(t,x,n)
    !$omp end critical
    i = im1+1
    ! 
    ! denormalize tension factor
    ! 
    sigmap = ABS(sigma)*REAL(n-1,kind=real_8)/(x(n)-x(1))
    ! 
    ! set up and perform interpolation
    ! 
    del1 = t-x(im1)
    del2 = x(i)-t
    dels = x(i)-x(im1)
    sum = (y(i)*del1+y(im1)*del2)/dels
    IF (sigmap .NE. 0._real_8) go to 1
    curv2 = sum-del1*del2*(yp(i)*(del1+dels)+yp(im1)*&
         (del2+dels))/(6._real_8*dels)
    RETURN
1   sigdel = sigmap*dels
    CALL snhcsh (ss,dummy,sigdel,-1)
    CALL snhcsh (s1,dummy,sigmap*del1,-1)
    CALL snhcsh (s2,dummy,sigmap*del2,-1)
    curv2 = sum+(yp(i)*del1*(s1-ss)+yp(im1)*del2*(s2-ss))/&
         (sigdel*sigmap*(1._real_8+ss))
    RETURN
  END FUNCTION curv2

  REAL(real_8) FUNCTION curvd (t,n,x,y,yp,sigma)
    INTEGER :: n
    REAL(real_8) :: t,x(n),y(n),yp(n),sigma
    REAL(real_8) :: c1, c2, del1, del2, dels, dummy, sigdel, sigmap, ss, sum
    INTEGER :: i, im1
    ! 
    ! coded by alan kaylor cline
    ! from fitpack -- january 26, 1987
    ! a curve and surface fitting package
    ! a product of pleasant valley software
    ! 8603 altus cove, austin, texas 78759, usa
    ! 
    ! this function differentiates a curve at a given point
    ! using a spline under tension. the subroutine curv1 should
    ! be called earlier to determine certain necessary
    ! parameters.
    ! 
    ! on input--
    ! 
    ! t contains a real value at which the derivative is to be
    ! determined.
    ! 
    ! n contains the number of points which were specified to
    ! determine the curve.
    ! 
    ! x and y are arrays containing the abscissae and
    ! ordinates, respectively, of the specified points.
    ! 
    ! yp is an array of second derivative values of the curve
    ! at the nodes.
    ! 
    ! and
    ! 
    ! sigma contains the tension factor (its sign is ignored).
    ! 
    ! the parameters n, x, y, yp, and sigma should be input
    ! unaltered from the output of curv1.
    ! 
    ! on output--
    ! 
    ! curvd contains the derivative value.
    ! 
    ! none of the input parameters are altered.
    ! 
    ! this function references package modules intrvl and
    ! snhcsh.
    ! 
    ! -----------------------------------------------------------
    ! 
    ! determine interval
    ! 
    !$omp critical
    im1 = intrvl(t,x,n)
    !$omp end critical
    i = im1+1
    ! 
    ! denormalize tension factor
    ! 
    sigmap = ABS(sigma)*REAL(n-1,kind=real_8)/(x(n)-x(1))
    ! 
    ! set up and perform differentiation
    ! 
    del1 = t-x(im1)
    del2 = x(i)-t
    dels = x(i)-x(im1)
    sum = (y(i)-y(im1))/dels
    IF (sigmap .NE. 0._real_8) go to 1
    curvd = sum+(yp(i)*(2._real_8*del1*del1-del2*(del1+dels))-&
         yp(im1)*(2._real_8*del2*del2-del1*(del2+dels)))&
         /(6._real_8*dels)
    RETURN
1   sigdel = sigmap*dels
    CALL snhcsh (ss,dummy,sigdel,-1)
    CALL snhcsh (dummy,c1,sigmap*del1,1)
    CALL snhcsh (dummy,c2,sigmap*del2,1)
    curvd = sum+(yp(i)*(c1-ss)-yp(im1)*(c2-ss))/&
         (sigdel*sigmap*(1._real_8+ss))
    RETURN
  END FUNCTION curvd

  REAL(real_8) FUNCTION curvi (xl,xu,n,x,y,yp,sigma)
    INTEGER :: n
    REAL(real_8) :: xl,xu,x(n),y(n),yp(n),sigma
    REAL(real_8) :: c1,c2,cl1,cl2,cs,cu1,cu2,del1,del2,deli,dell1,dell2,dels,&
         delu1,delu2,dummy,sigmap,ss,ssign,sum,t1,t2,xxl,xxu
    INTEGER :: i,il,ilm1,ilp1,iu,ium1

    ! 
    ! coded by alan kaylor cline
    ! from fitpack -- january 26, 1987
    ! a curve and surface fitting package
    ! a product of pleasant valley software
    ! 8603 altus cove, austin, texas 78759, usa
    ! 
    ! this function integrates a curve specified by a spline
    ! under tension between two given limits. the subroutine
    ! curv1 should be called earlier to determine necessary
    ! parameters.
    ! 
    ! on input--
    ! 
    ! xl and xu contain the upper and lower limits of inte-
    ! gration, respectively. (sl need not be less than or
    ! equal to xu, curvi (xl,xu,...) .eq. -curvi (xu,xl,...) ).
    ! 
    ! n contains the number of points which were specified to
    ! determine the curve.
    ! 
    ! x and y are arrays containing the abscissae and
    ! ordinates, respectively, of the specified points.
    ! 
    ! yp is an array from subroutine curv1 containing
    ! the values of the second derivatives at the nodes.
    ! 
    ! and
    ! 
    ! sigma contains the tension factor (its sign is ignored).
    ! 
    ! the parameters n, x, y, yp, and sigma should be input
    ! unaltered from the output of curv1.
    ! 
    ! on output--
    ! 
    ! curvi contains the integral value.
    ! 
    ! none of the input parameters are altered.
    ! 
    ! this function references package modules intrvl and
    ! snhcsh.
    ! 
    ! -----------------------------------------------------------
    ! 
    ! denormalize tension factor
    ! 
    sigmap = ABS(sigma)*REAL(n-1,kind=real_8)/(x(n)-x(1))
    ! 
    ! determine actual upper and lower bounds
    ! 
    xxl = xl
    xxu = xu
    ssign = 1._real_8
    IF (xl .LT. xu) go to 1
    xxl = xu
    xxu = xl
    ssign = -1._real_8
    IF (xl .GT. xu) go to 1
    ! 
    ! return zero if xl .eq. xu
    ! 
    curvi = 0._real_8
    RETURN
    ! 
    ! search for cntl%proper intervals
    ! 
1   CONTINUE
    !$omp critical
    ilm1 = intrvl (xxl,x,n)
    !$omp end critical
    il = ilm1+1
    !$omp critical
    ium1 = intrvl (xxu,x,n)
    !$omp end critical
    iu = ium1+1
    IF (il .EQ. iu) go to 8
    ! 
    ! integrate from xxl to x(il)
    ! 
    sum = 0._real_8
    IF (xxl .EQ. x(il)) go to 3
    del1 = xxl-x(ilm1)
    del2 = x(il)-xxl
    dels = x(il)-x(ilm1)
    t1 = (del1+dels)*del2/(2._real_8*dels)
    t2 = del2*del2/(2._real_8*dels)
    sum = t1*y(il)+t2*y(ilm1)
    IF (sigma .EQ. 0._real_8) go to 2
    CALL snhcsh (dummy,c1,sigmap*del1,2)
    CALL snhcsh (dummy,c2,sigmap*del2,2)
    CALL snhcsh (ss,cs,sigmap*dels,3)
    sum = sum+((dels*dels*(cs-ss/2._real_8)-del1*del1*(c1-ss/2._real_8))&
         *yp(il)+del2*del2*(c2-ss/2._real_8)*yp(ilm1))/&
         (sigmap*sigmap*dels*(1._real_8+ss))
    go to 3
2   sum = sum-t1*t1*dels*yp(il)/6._real_8&
         -t2*(del1*(del2+dels)+dels*dels)*yp(ilm1)/12._real_8
    ! 
    ! integrate over interior intervals
    ! 
3   IF (iu-il .EQ. 1) go to 6
    ilp1 = il+1
    DO 5 i = ilp1,ium1
       dels = x(i)-x(i-1)
       sum = sum+(y(i)+y(i-1))*dels/2._real_8
       IF (sigma .EQ. 0._real_8) go to 4
       CALL snhcsh (ss,cs,sigmap*dels,3)
       sum = sum+(yp(i)+yp(i-1))*dels*(cs-ss/2._real_8)/&
            (sigmap*sigmap*(1._real_8+ss))
       go to 5
4      sum = sum-(yp(i)+yp(i-1))*dels*dels*dels/24._real_8
5   CONTINUE
    ! 
    ! integrate from x(iu-1) to xxu
    ! 
6   IF (xxu .EQ. x(ium1)) go to 10
    del1 = xxu-x(ium1)
    del2 = x(iu)-xxu
    dels = x(iu)-x(ium1)
    t1 = del1*del1/(2._real_8*dels)
    t2 = (del2+dels)*del1/(2._real_8*dels)
    sum = sum+t1*y(iu)+t2*y(ium1)
    IF (sigma .EQ. 0._real_8) go to 7
    CALL snhcsh (dummy,c1,sigmap*del1,2)
    CALL snhcsh (dummy,c2,sigmap*del2,2)
    CALL snhcsh (ss,cs,sigmap*dels,3)
    sum = sum+(yp(iu)*del1*del1*(c1-ss/2._real_8)+yp(ium1)*&
         (dels*dels*(cs-ss/2._real_8)-del2*del2*(c2-ss/2._real_8)))&
         /(sigmap*sigmap*dels*(1._real_8+ss))
    go to 10
7   sum = sum-t1*(del2*(del1+dels)+dels*dels)*yp(iu)/12._real_8&
         -t2*t2*dels*yp(ium1)/6._real_8
    go to 10
    ! 
    ! integrate from xxl to xxu
    ! 
8   delu1 = xxu-x(ium1)
    delu2 = x(iu)-xxu
    dell1 = xxl-x(ium1)
    dell2 = x(iu)-xxl
    dels = x(iu)-x(ium1)
    deli = xxu-xxl
    t1 = (delu1+dell1)*deli/(2._real_8*dels)
    t2 = (delu2+dell2)*deli/(2._real_8*dels)
    sum = t1*y(iu)+t2*y(ium1)
    IF (sigma .EQ. 0._real_8) go to 9
    CALL snhcsh (dummy,cu1,sigmap*delu1,2)
    CALL snhcsh (dummy,cu2,sigmap*delu2,2)
    CALL snhcsh (dummy,cl1,sigmap*dell1,2)
    CALL snhcsh (dummy,cl2,sigmap*dell2,2)
    CALL snhcsh (ss,dummy,sigmap*dels,-1)
    sum = sum+(yp(iu)*(delu1*delu1*(cu1-ss/2._real_8)&
         -dell1*dell1*(cl1-ss/2._real_8))&
         +yp(ium1)*(dell2*dell2*(cl2-ss/2._real_8)&
         -delu2*delu2*(cu2-ss/2._real_8)))/&
         (sigmap*sigmap*dels*(1._real_8+ss))
    go to 10
9   sum = sum-t1*(delu2*(dels+delu1)+dell2*(dels+dell1))*&
         yp(iu)/12._real_8&
         -t2*(dell1*(dels+dell2)+delu1*(dels+delu2))*&
         yp(ium1)/12._real_8
    ! 
    ! correct sign and return
    ! 
10  curvi = ssign*sum
    RETURN
  END FUNCTION curvi

  SUBROUTINE curvp1 (n,x,y,p,yp,temp,sigma,ierr)
    INTEGER                                  :: n
    REAL(real_8)                             :: x(n), y(n), p, yp(n), &
                                                temp(1), sigma
    INTEGER                                  :: ierr

    INTEGER                                  :: i, ibak, nm1, np1, npi, npibak
    REAL(real_8)                             :: delx1, delx2, diag, diag1, &
                                                diag2, dx1, dx2, sdiag1, &
                                                sdiag2, sigmap, ypn

! 
! coded by alan kaylor cline
! from fitpack -- january 26, 1987
! a curve and surface fitting package
! a product of pleasant valley software
! 8603 altus cove, austin, texas 78759, usa
! 
! this subroutine determines the parameters necessary to
! compute a periodic interpolatory spline under tension
! through a sequence of functional values. for actual ends
! of the curve may be specified or omitted.  for actual
! computation of points on the curve it is necessary to call
! the function curvp2.
! 
! on input--
! 
! n is the number of values to be interpolated (n.ge.2).
! 
! x is an array of the n increasing abscissae of the
! functional values.
! 
! y is an array of the n ordinates of the values, (i. e.
! y(k) is the functional value corresponding to x(k) ).
! 
! p is the period (p .gt. x(n)-x(1)).
! 
! yp is an array of length at least n.
! 
! temp is an array of length at least 2*n which is used
! for scratch storage.
! 
! and
! 
! sigma contains the tension factor.  this value indicates
! the curviness desired. if abs(sigma) is nearly zero
! (e.g. .001) the resulting curve is approximately a
! cubic spline. if abs(sigma) is large (e.g. 50.) the
! resulting curve is nearly a polygonal line. if sigma
! equals zero a cubic spline results.  a standard value
! for sigma is approximately 1. in absolute value.
! 
! on output--
! 
! yp contains the values of the second derivative of the
! curve at the given nodes.
! 
! ierr contains an error flag,
! = 0 for normal return,
! = 1 if n is less than 2,
! = 2 if p is less than or equal to x(n)-x(1),
! = 3 if x-values are not strictly increasing.
! 
! and
! 
! n, x, y, and sigma are unaltered.
! 
! this subroutine references package modules terms and
! snhcsh.
! 
! -----------------------------------------------------------
! 

    nm1 = n-1
    np1 = n+1
    ierr = 0
    IF (n .LE. 1) go to 6
    IF (p .LE. x(n)-x(1) .OR. p .LE. 0._real_8) go to 7
    ! 
    ! denormalize tension factor
    ! 
    sigmap = ABS(sigma)*REAL(n,kind=real_8)/p
    ! 
    ! set up right hand side and tridiagonal system for yp and
    ! perform forward elimination
    ! 
    delx1 = p-(x(n)-x(1))
    dx1 = (y(1)-y(n))/delx1
    CALL terms (diag1,sdiag1,sigmap,delx1)
    delx2 = x(2)-x(1)
    IF (delx2 .LE. 0._real_8) go to 8
    dx2 = (y(2)-y(1))/delx2
    CALL terms (diag2,sdiag2,sigmap,delx2)
    diag = diag1+diag2
    yp(1) = (dx2-dx1)/diag
    temp(np1) = -sdiag1/diag
    temp(1) = sdiag2/diag
    dx1 = dx2
    diag1 = diag2
    sdiag1 = sdiag2
    IF (n .EQ. 2) go to 2
    DO 1 i = 2,nm1
       npi = n+i
       delx2 = x(i+1)-x(i)
       IF (delx2 .LE. 0._real_8) go to 8
       dx2 = (y(i+1)-y(i))/delx2
       CALL terms (diag2,sdiag2,sigmap,delx2)
       diag = diag1+diag2-sdiag1*temp(i-1)
       yp(i) = (dx2-dx1-sdiag1*yp(i-1))/diag
       temp(npi) = -temp(npi-1)*sdiag1/diag
       temp(i) = sdiag2/diag
       dx1 = dx2
       diag1 = diag2
1   sdiag1 = sdiag2
2   delx2 = p-(x(n)-x(1))
    dx2 = (y(1)-y(n))/delx2
    CALL terms (diag2,sdiag2,sigmap,delx2)
    yp(n) = dx2-dx1
    temp(nm1) = temp(2*n-1)-temp(nm1)
    IF (n .EQ. 2) go to 4
    ! 
    ! perform first step of back substitution
    ! 
    DO 3 i = 3,n
       ibak = np1-i
       npibak =n+ibak
       yp(ibak) = yp(ibak)-temp(ibak)*yp(ibak+1)
3   temp(ibak) =temp(npibak)-temp(ibak)*temp(ibak+1)
4   yp(n) = (yp(n)-sdiag2*yp(1)-sdiag1*yp(nm1))/&
         (diag1+diag2+sdiag2*temp(1)+sdiag1*temp(nm1))
    ! 
    ! perform second step of back substitution
    ! 
    ypn =   yp(n)
    DO 5 i = 1,nm1
5   yp(i) = yp(i)+temp(i)*ypn
    RETURN
    ! 
    ! too few points
    ! 
6   ierr = 1
    RETURN
    ! 
    ! period too small
    ! 
7   ierr = 2
    RETURN
    ! 
    ! x-values not strictly increasing
    ! 
8   ierr = 3
    RETURN
  END SUBROUTINE curvp1

  SUBROUTINE curvps (n,x,y,p,d,isw,s,eps,ys,ysp,sigma,&
       temp,ierr)
    ! 
    INTEGER                                  :: n
    REAL(real_8)                             :: x(n), y(n), p, d(n)
    INTEGER                                  :: isw
    REAL(real_8)                             :: s, eps, ys(n), ysp(n), sigma, &
                                                temp(n,11)
    INTEGER                                  :: ierr

! 
! coded by alan kaylor cline
! from fitpack -- january 26, 1987
! a curve and surface fitting package
! a product of pleasant valley software
! 8603 altus cove, austin, texas 78759, usa
! 
! this subroutine determines the parameters necessary to
! compute a periodic smoothing spline under tension. for a
! given increasing sequence of abscissae (x(i)), i = 1,...,n
! and associated ordinates (y(i)), i = 1,...,n, letting p be
! the period, x(n+1) = x(1)+p, and y(n+1) = y(1), the
! function determined minimizes the summation from i = 1 to
! n of the square of the second derivative of f plus sigma
! squared times the difference of the first derivative of f
! and (f(x(i+1))-f(x(i)))/(x(i+1)-x(i)) squared, over all
! functions f with period p and two continuous derivatives
! such that the summation of the square of
! (f(x(i))-y(i))/d(i) is less than or equal to a given
! constant s, where (d(i)), i = 1,...,n are a given set of
! observation weights. the function determined is a periodic
! spline under tension with third derivative discontinuities
! at (x(i)) i = 1,...,n (and all periodic translations of
! these values). for actual computation of points on the
! curve it is necessary to call the function curvp2. the
! determination of the curve is performed by subroutine
! curvpp, the subroutin curvps only decomposes the workspace
! for curvpp.
! 
! on input--
! 
! n is the number of values to be smoothed (n.ge.2).
! 
! x is an array of the n increasing abscissae of the
! values to be smoothed.
! 
! y is an array of the n ordinates of the values to be
! smoothed, (i. e. y(k) is the functional value
! corresponding to x(k) ).
! 
! p is the period (p .gt. x(n)-x(1)).
! 
! d is a parameter containing the observation weights.
! this may either be an array of length n or a scalar
! (interpreted as a constant). the value of d
! corresponding to the observation (x(k),y(k)) should
! be an approximation to the standard deviation of error.
! 
! isw contains a switch indicating whether the parameter
! d is to be considered a vector or a scalar,
! = 0 if d is an array of length n,
! = 1 if d is a scalar.
! 
! s contains the value controlling the smoothing. this
! must be non-negative. for s equal to zero, the
! subroutine does interpolation, larger values lead to
! smoother funtions. if parameter d contains standard
! deviation estimates, a reasonable value for s is
! float(n).
! 
! eps contains a tolerance on the relative precision to
! which s is to be interpreted. this must be greater than
! or equal to zero and less than or equal to one. a
! reasonable value for eps is sqrt(2./float(n)).
! 
! ys is an array of length at least n.
! 
! ysp is an array of length at least n.
! 
! sigma contains the tension factor. this value indicates
! the degree to which the first derivative part of the
! smoothing functional is emphasized. if sigma is nearly
! zero (e. g. .001) the resulting curve is approximately a
! cubic spline. if sigma is large (e. g. 50.) the
! resulting curve is nearly a polygonal line. if sigma
! equals zero a cubic spline results. a standard value for
! sigma is approximately 1.
! 
! and
! 
! temp is an array of length at least 11*n which is used
! for scratch storage.
! 
! on output--
! 
! ys contains the smoothed ordinate values.
! 
! ysp contains the values of the second derivative of the
! smoothed curve at the given nodes.
! 
! ierr contains an error flag,
! = 0 for normal return,
! = 1 if n is less than 2,
! = 2 if s is negative,
! = 3 if eps is negative or greater than one,
! = 4 if x-values are not strictly increasing,
! = 5 if a d-value is non-positive,
! = 6 if p is less than or equal to x(n)-x(1).
! 
! and
! 
! n, x, y, p, d, isw, s, eps, and sigma are unaltered.
! 
! this subroutine references package modules curvpp, terms,
! and snhcsh.
! 
! -----------------------------------------------------------
! 
! decompose temp into eleven arrays and call curvpp
! 

    CALL curvpp (n,x,y,p,d,isw,s,eps,ys,ysp,sigma,&
         temp(1,1),temp(1,2),temp(1,3),temp(1,4),&
         temp(1,5),temp(1,6),temp(1,7),temp(1,8),&
         temp(1,9),temp(1,10),temp(1,11),ierr)
    RETURN
  END SUBROUTINE curvps


  REAL(real_8) FUNCTION curvp2 (t,n,x,y,p,yp,sigma)
    INTEGER :: n
    REAL(real_8) :: t,x(n),y(n),p,yp(n),sigma
    REAL(real_8) :: del1,del2,dels,dummy,s1,s2,sigdel,sigmap,ss,sum,tp
    INTEGER :: i,im1

    ! 
    ! coded by alan kaylor cline
    ! from fitpack -- january 26, 1987
    ! a curve and surface fitting package
    ! a product of pleasant valley software
    ! 8603 altus cove, austin, texas 78759, usa
    ! 
    ! this function interpolates a curve at a given point using
    ! a periodic spline under tension. the subroutine curvp1
    ! should be called earlier to determine certain necessary
    ! parameters.
    ! 
    ! on input--
    ! 
    ! t contains a real value to be mapped onto the interpo-
    ! lating curve.
    ! 
    ! n contains the number of points which were specified to
    ! determine the curve.
    ! 
    ! x and y are arrays containing the abscissae and
    ! ordinates, respectively, of the specified points.
    ! 
    ! p contains the period.
    ! 
    ! yp is an array of second derivative values of the curve
    ! at the nodes.
    ! 
    ! and
    ! 
    ! sigma contains the tension factor (its sign is ignored).
    ! 
    ! the parameters n, x, y, p, yp, and sigma should be input
    ! unaltered from the output of curvp1.
    ! 
    ! on output--
    ! 
    ! curvp2 contains the interpolated value.
    ! 
    ! none of the input parameters are altered.
    ! 
    ! this function references package modules intrvp and
    ! snhcsh.
    ! 
    ! -----------------------------------------------------------
    ! 
    ! determine interval
    ! 
    !$omp critical
    im1 = intrvp (t,x,n,p,tp)
    !$omp end critical
    i = im1+1
    ! 
    ! denormalize tension factor
    ! 
    sigmap = ABS(sigma)*REAL(n,kind=real_8)/p
    ! 
    ! set up and perform interpolation
    ! 
    del1 = tp-x(im1)
    IF (im1 .EQ. n) go to 1
    del2 = x(i)-tp
    dels = x(i)-x(im1)
    go to 2
1   i = 1
    del2 = x(1)+p-tp
    dels = p-(x(n)-x(1))
2   sum = (y(i)*del1+y(im1)*del2)/dels
    IF (sigmap .NE. 0._real_8) go to 3
    curvp2 = sum-del1*del2*(yp(i)*(del1+dels)+yp(im1)*&
         (del2+dels))/(6._real_8*dels)
    RETURN
3   sigdel = sigmap*dels
    CALL snhcsh (ss,dummy,sigdel,-1)
    CALL snhcsh (s1,dummy,sigmap*del1,-1)
    CALL snhcsh (s2,dummy,sigmap*del2,-1)
    curvp2 = sum+(yp(i)*del1*(s1-ss)+yp(im1)*del2*(s2-ss))/&
         (sigdel*sigmap*(1._real_8+ss))
    RETURN
  END FUNCTION curvp2

  REAL(real_8)  FUNCTION curvpi (xl,xu,n,x,y,p,yp,sigma)
    INTEGER :: n
    REAL(real_8) :: xl,xu,x(n),y(n),p,yp(n),sigma
    REAL(real_8) :: c1,c2,cl1,cl2,cs,cu1,cu2,del1,del2,deli,dell1,dell2,dels,delu1,&
         delu2,dummy,s1,s2,s3,s4,s5,s6,s7,s8,si,sigmap,so,ss,t1,t2,x1pp,xil,xiu,xsave,xxl,xxu
    INTEGER :: i,ideltp,ii,il,ilm1,ilp1,im1,isave,isign,iu,ium1,iup1,lper,np1

    ! 
    ! coded by alan kaylor cline
    ! from fitpack -- january 26, 1987
    ! a curve and surface fitting package
    ! a product of pleasant valley software
    ! 8603 altus cove, austin, texas 78759, usa
    ! 
    ! this function integrates a curve specified by a periodic
    ! spline under tension between two given limits. the
    ! subroutine curvp1 should be called earlier to determine
    ! necessary parameters.
    ! 
    ! on input--
    ! 
    ! xl and xu contain the upper and lower limits of inte-
    ! gration, respectively. (sl need not be less than or
    ! equal to xu, curvpi (xl,xu,...) .eq. -curvpi (xu,xl,...) ).
    ! 
    ! n contains the number of points which were specified to
    ! determine the curve.
    ! 
    ! x and y are arrays containing the abscissae and
    ! ordinates, respectively, of the specified points.
    ! 
    ! p contains the period.
    ! 
    ! yp is an array from subroutine curvp1 containing
    ! the values of the second derivatives at the nodes.
    ! 
    ! and
    ! 
    ! sigma contains the tension factor (its sign is ignored).
    ! 
    ! the parameters n, x, y, p, yp, and sigma should be input
    ! unaltered from the output of curvp1.
    ! 
    ! on output--
    ! 
    ! 
    ! curvpi contains the integral value.
    ! 
    ! none of the input parameters are altered.
    ! 
    ! this function references package modules intrvp and
    ! snhcsh.
    ! 
    ! --------------------------------------------------------------
    ! 
    INTEGER :: uper
    LOGICAL :: bdy
    ! 
    ! denormalize tension factor
    ! 
    sigmap = ABS(sigma)*REAL(n,kind=real_8)/p
    ! 
    ! determine actual upper and lower bounds
    ! 
    x1pp = x(1)+p
    isign = 1
    !$omp critical
    ilm1 = intrvp (xl,x,n,p,xxl)
    !$omp end critical
    lper = INT((xl-x(1))/p)
    IF (xl .LT. x(1)) lper = lper-1
    !$omp critical
    ium1 = intrvp (xu,x,n,p,xxu)
    !$omp end critical
    uper = INT((xu-x(1))/p)
    IF (xu .LT. x(1)) uper = uper-1
    ideltp = uper-lper
    bdy = REAL(ideltp,kind=real_8)*(xxu-xxl) .LT. 0._real_8
    IF ((ideltp .EQ. 0 .AND. xxu .LT. xxl) .OR.&
         ideltp .LT. 0) isign = -1
    IF (bdy) ideltp = ideltp-isign
    IF (xxu .GE. xxl) go to 1
    xsave = xxl
    xxl = xxu
    xxu = xsave
    isave = ilm1
    ilm1 = ium1
    ium1 = isave
1   il = ilm1+1
    IF (ilm1 .EQ. n) il = 1
    xil = x(il)
    IF (ilm1 .EQ. n) xil = x1pp
    iu = ium1+1
    IF (ium1 .EQ. n) iu = 1
    xiu = x(iu)
    IF (ium1 .EQ. n) xiu = x1pp
    s1 = 0._real_8
    IF (ilm1 .EQ. 1 .OR. (ideltp .EQ. 0 .AND.&
         .NOT. bdy)) go to 4
    ! 
    ! integrate from x(1) to x(ilm1), store in s1
    ! 
    DO 3 i = 2,ilm1
       dels = x(i)-x(i-1)
       s1 = s1+(y(i)+y(i-1))*dels/2._real_8
       IF (sigma .EQ. 0._real_8) go to 2
       CALL snhcsh (ss,cs,sigmap*dels,3)
       s1 = s1+(yp(i)+yp(i-1))*dels*(cs-ss/2._real_8)/&
            (sigmap*sigmap*(1._real_8+ss))
       go to 3
2      s1 = s1-(yp(i)+yp(i-1))*dels*dels*dels/24._real_8
3   CONTINUE
4   s2 = 0._real_8
    IF (x(ilm1) .GE. xxl .OR. (ideltp .EQ. 0&
         .AND. .NOT. bdy)) go to 6
    ! 
    ! integrate from x(ilm1) to xxl, store in s2
    ! 
    del1 = xxl-x(ilm1)
    del2 = xil-xxl
    dels = xil-x(ilm1)
    t1 = del1*del1/(2._real_8*dels)
    t2 = (del2+dels)*del1/(2._real_8*dels)
    s2 = t1*y(il)+t2*y(ilm1)
    IF (sigma .EQ. 0._real_8) go to 5
    CALL snhcsh (dummy,c1,sigmap*del1,2)
    CALL snhcsh (dummy,c2,sigmap*del2,2)
    CALL snhcsh (ss,cs,sigmap*dels,3)
    s2 = s2+(yp(il)*del1*del1*(c1-ss/2._real_8)+yp(ilm1)*&
         (dels*dels*(cs-ss/2._real_8)-del2*del2*(c2-ss/2._real_8)))&
         /(sigmap*sigmap*dels*(1._real_8+ss))
    go to 6
5   s2 = s2-t1*(del2*(del1+dels)&
         +dels*dels)*yp(il)/12._real_8&
         -t2*t2*dels*yp(ilm1)/6._real_8
6   s3 = 0._real_8
    IF (xxl .GE. xil .OR. (ideltp .EQ. 0 .AND. bdy)&
         .OR. ilm1 .EQ. ium1) go to 8
    ! 
    ! integrate from xxl to xil, store in s3
    ! 
    del1 = xxl-x(ilm1)
    del2 = xil-xxl
    dels = xil-x(ilm1)
    t1 = (del1+dels)*del2/(2._real_8*dels)
    t2 = del2*del2/(2._real_8*dels)
    s3 = t1*y(il)+t2*y(ilm1)
    IF (sigma .EQ. 0._real_8) go to 7
    CALL snhcsh (dummy,c1,sigmap*del1,2)
    CALL snhcsh (dummy,c2,sigmap*del2,2)
    CALL snhcsh (ss,cs,sigmap*dels,3)
    s3 = s3+((dels*dels*(cs-ss/2._real_8)-del1*del1*(c1-ss/2._real_8))&
         *yp(il)+del2*del2*(c2-ss/2._real_8)*yp(ilm1))/&
         (sigmap*sigmap*dels*(1._real_8+ss))
    go to 8
7   s3 = s3-t1*t1*dels*yp(il)/6._real_8&
         -t2*(del1*(del2+dels)+dels*dels)*&
         yp(ilm1)/12._real_8
8   s4 = 0._real_8
    IF (ilm1 .GE. ium1-1 .OR. (ideltp .EQ. 0 .AND. bdy))&
         go to 11
    ! 
    ! integrate from xil to x(ium1), store in s4
    ! 
    ilp1 = il+1
    DO 10 i = ilp1,ium1
       dels = x(i)-x(i-1)
       s4 = s4+(y(i)+y(i-1))*dels/2._real_8
       IF (sigma .EQ. 0._real_8) go to 9
       CALL snhcsh (ss,cs,sigmap*dels,3)
       s4 = s4+(yp(i)+yp(i-1))*dels*(cs-ss/2._real_8)/&
            (sigmap*sigmap*(1._real_8+ss))
       go to 10
9      s4 = s4-(yp(i)+yp(i-1))*dels*dels*dels/24._real_8
10  CONTINUE
11  s5 = 0._real_8
    IF (x(ium1) .GE. xxu .OR. (ideltp .EQ. 0 .AND. bdy)&
         .OR. ilm1 .EQ. ium1) go to 13
    ! 
    ! integrate from x(ium1) to xxu, store in s5
    ! 
    del1 = xxu-x(ium1)
    del2 = xiu-xxu
    dels = xiu-x(ium1)
    t1 = del1*del1/(2._real_8*dels)
    t2 = (del2+dels)*del1/(2._real_8*dels)
    s5 = t1*y(iu)+t2*y(ium1)
    IF (sigma .EQ. 0._real_8) go to 12
    CALL snhcsh (dummy,c1,sigmap*del1,2)
    CALL snhcsh (dummy,c2,sigmap*del2,2)
    CALL snhcsh (ss,cs,sigmap*dels,3)
    s5 = s5+(yp(iu)*del1*del1*(c1-ss/2._real_8)+yp(ium1)*&
         (dels*dels*(cs-ss/2._real_8)-del2*del2*(c2-ss/2._real_8)))&
         /(sigmap*sigmap*dels*(1._real_8+ss))
    go to 13
12  s5 = s5-t1*(del2*(del1+dels)&
         +dels*dels)*yp(iu)/12._real_8&
         -t2*t2*dels*yp(ium1)/6._real_8
13  s6 = 0._real_8
    IF (xxu .GE. xiu .OR. (ideltp .EQ. 0 .AND.&
         .NOT. bdy)) go to 15
    ! 
    ! integrate from xxu to xiu, store in s6
    ! 
    del1 = xxu-x(ium1)
    del2 = xiu-xxu
    dels = xiu-x(ium1)
    t1 = (del1+dels)*del2/(2._real_8*dels)
    t2 = del2*del2/(2._real_8*dels)
    s6 = t1*y(iu)+t2*y(ium1)
    IF (sigma .EQ. 0._real_8) go to 14
    CALL snhcsh (dummy,c1,sigmap*del1,2)
    CALL snhcsh (dummy,c2,sigmap*del2,2)
    CALL snhcsh (ss,cs,sigmap*dels,3)
    s6 = s6+((dels*dels*(cs-ss/2._real_8)-del1*del1*(c1-ss/2._real_8))&
         *yp(iu)+del2*del2*(c2-ss/2._real_8)*yp(ium1))/&
         (sigmap*sigmap*dels*(1._real_8+ss))
    go to 15
14  s6 = s6-t1*t1*dels*yp(iu)/6._real_8&
         -t2*(del1*(del2+dels)+dels*dels)*&
         yp(ium1)/12._real_8
15  s7 = 0._real_8
    IF (iu .EQ. 1 .OR. (ideltp .EQ. 0 .AND. .NOT. bdy))&
         go to 18
    ! 
    ! integrate from xiu to x1pp, store in s7
    ! 
    np1 = n+1
    iup1 = iu+1
    DO 17 ii = iup1,np1
       im1 = ii-1
       i = ii
       IF (i .EQ. np1) i=1
       dels = x(i)-x(im1)
       IF (dels .LE. 0._real_8) dels=dels+p
       s7 = s7+(y(i)+y(im1))*dels/2._real_8
       IF (sigma .EQ. 0._real_8) go to 16
       CALL snhcsh (ss,cs,sigmap*dels,3)
       s7 = s7+(yp(i)+yp(im1))*dels*(cs-ss/2._real_8)/&
            (sigmap*sigmap*(1._real_8+ss))
       go to 17
16     s7 = s7-(yp(i)+yp(im1))*dels*dels*dels/24._real_8
17  CONTINUE
18  s8 = 0._real_8
    IF (ilm1 .LT. ium1 .OR. (ideltp .EQ. 0 .AND. bdy))&
         go to 20
    ! 
    ! integrate from xxl to xxu, store in s8
    ! 
    delu1 = xxu-x(ium1)
    delu2 = xiu-xxu
    dell1 = xxl-x(ium1)
    dell2 = xiu-xxl
    dels = xiu-x(ium1)
    deli = xxu-xxl
    t1 = (delu1+dell1)*deli/(2._real_8*dels)
    t2 = (delu2+dell2)*deli/(2._real_8*dels)
    s8 = t1*y(iu)+t2*y(ium1)
    IF (sigma .EQ. 0._real_8) go to 19
    CALL snhcsh (dummy,cu1,sigmap*delu1,2)
    CALL snhcsh (dummy,cu2,sigmap*delu2,2)
    CALL snhcsh (dummy,cl1,sigmap*dell1,2)
    CALL snhcsh (dummy,cl2,sigmap*dell2,2)
    CALL snhcsh (ss,dummy,sigmap*dels,-1)
    s8 = s8+(yp(iu)*(delu1*delu1*(cu1-ss/2._real_8)&
         -dell1*dell1*(cl1-ss/2._real_8))&
         +yp(ium1)*(dell2*dell2*(cl2-ss/2._real_8)&
         -delu2*delu2*(cu2-ss/2._real_8)))/&
         (sigmap*sigmap*dels*(1._real_8+ss))
    go to 20
19  s8 = s8-t1*(delu2*(dels+delu1)&
         +dell2*(dels+dell1))*yp(iu)/12._real_8&
         -t2*(dell1*(dels+dell2)&
         +delu1*(dels+delu2))*yp(ium1)/12._real_8
20  so = s1+s2+s6+s7
    si = s3+s4+s5+s8
    IF (bdy) go to 21
    curvpi = REAL(ideltp,kind=real_8)*(so+si)+REAL(isign,kind=real_8)*si
    RETURN
21  curvpi = REAL(ideltp,kind=real_8)*(so+si)+REAL(isign,kind=real_8)*so
    RETURN
  END FUNCTION curvpi

  SUBROUTINE ceez (del1,del2,sigma,c1,c2,c3,n)
    REAL(real_8)                             :: del1, del2, sigma, c1, c2, c3
    INTEGER                                  :: n

    REAL(real_8)                             :: coshm1, coshm2, del, delm, &
                                                delp, denom, dummy, sinhmm, &
                                                sinhmp

! 
! coded by alan kaylor cline
! from fitpack -- january 26, 1987
! a curve and surface fitting package
! a product of pleasant valley software
! 8603 altus cove, austin, texas 78759, usa
! 
! this subroutine determines the coefficients c1, c2, and c3
! used to determine endpoint slopes. specifically, if
! function values y1, y2, and y3 are given at points x1, x2,
! and x3, respectively, the quantity c1*y1 + c2*y2 + c3*y3
! is the value of the derivative at x1 of a spline under
! tension (with tension factor sigma) passing through the
! three points and having third derivative equal to zero at
! x1. optionally, only two values, c1 and c2 are determined.
! 
! on input--
! 
! del1 is x2-x1 (.gt. 0._real_8).
! 
! del2 is x3-x1 (.gt. 0._real_8). if n .eq. 2, this parameter is
! ignored.
! 
! sigma is the tension factor.
! 
! and
! 
! n is a switch indicating the number of coefficients to
! be returned. if n .eq. 2 only two coefficients are
! returned. otherwise all three are returned.
! 
! on output--
! 
! c1, c2, and c3 contain the coefficients.
! 
! none of the input parameters are altered.
! 
! this subroutine references package module snhcsh.
! 
! -----------------------------------------------------------
! 

    IF (n .EQ. 2) go to 2
    IF (sigma .NE. 0._real_8) go to 1
    del = del2-del1
    ! 
    ! tension .eq. 0.
    ! 
    c1 = -(del1+del2)/(del1*del2)
    c2 = del2/(del1*del)
    c3 = -del1/(del2*del)
    RETURN
    ! 
    ! tension .ne. 0.
    ! 
1   CALL snhcsh (dummy,coshm1,sigma*del1,1)
    CALL snhcsh (dummy,coshm2,sigma*del2,1)
    delp = sigma*(del2+del1)/2._real_8
    delm = sigma*(del2-del1)/2._real_8
    CALL snhcsh (sinhmp,dummy,delp,-1)
    CALL snhcsh (sinhmm,dummy,delm,-1)
    denom = coshm1*(del2-del1)-2._real_8*del1*delp*delm*&
         (1._real_8+sinhmp)*(1._real_8+sinhmm)
    c1 = 2._real_8*delp*delm*(1._real_8+sinhmp)*(1._real_8+sinhmm)/denom
    c2 = -coshm2/denom
    c3 = coshm1/denom
    RETURN
    ! 
    ! two coefficients
    ! 
2   c1 = -1._real_8/del1
    c2 = -c1
    RETURN
  END SUBROUTINE ceez

  SUBROUTINE curvpp (n,x,y,p,d,isw,s,eps,ys,ysp,sigma,&
       td,tsd1,hd,hsd1,hsd2,rd,rsd1,rsd2,&
       rnm1,rn,v,ierr)
    INTEGER                                  :: n
    REAL(real_8)                             :: x(n), y(n), p, d(n)
    INTEGER                                  :: isw
    REAL(real_8) :: s, eps, ys(n), ysp(n), sigma, td(n), tsd1(n), hd(n), &
      hsd1(n), hsd2(n), rd(n), rsd1(n), rsd2(n), rnm1(n), rn(n), v(n)
    INTEGER                                  :: ierr

    INTEGER                                  :: i, ibak, im1, ip1, nm1, nm2, &
                                                nm3
    REAL(real_8) :: alpha, alphap, beta, betap, betapp, con, delxi, delxi1, &
      delyi, delyi1, di, dim1, disq, f, g, h, hdi, hdim1, hsd11, hsd1ip, &
      hsd1p, q, rdn, rnm1sm, rnm1t, rnsm, rnt, rsd1i, rsd2i, sigmap, sl, &
      step, su, sum, sum2, sumd, sumn, sumnm1, sumy, tui, wi, wim1, wim2, &
      yspn, yspnm1

! 
! coded by alan kaylor cline
! from fitpack -- january 26, 1987
! a curve and surface fitting package
! a product of pleasant valley software
! 8603 altus cove, austin, texas 78759, usa
! 
! this subroutine determines the parameters necessary to
! compute a periodic smoothing spline under tension. for a
! given increasing sequence of abscissae (x(i)), i = 1,...,n
! and associated ordinates (y(i)), i = 1,...,n, letting p be
! the period, x(n+1) = x(1)+p, and y(n+1) = y(1), the
! function determined minimizes the summation from i = 1 to
! n of the square of the second derivative of f plus sigma
! squared times the difference of the first derivative of f
! and (f(x(i+1))-f(x(i)))/(x(i+1)-x(i)) squared, over all
! functions f with period p and two continuous derivatives
! such that the summation of the square of
! (f(x(i))-y(i))/d(i) is less than or equal to a given
! constant s, where (d(i)), i = 1,...,n are a given set of
! observation weights. the function determined is a periodic
! spline under tension with third derivative discontinuities
! at (x(i)) i = 1,...,n (and all periodic translations of
! these values). for actual computation of points on the
! curve it is necessary to call the function curvp2.
! 
! on input--
! 
! n is the number of values to be smoothed (n.ge.2).
! 
! x is an array of the n increasing abscissae of the
! values to be smoothed.
! 
! y is an array of the n ordinates of the values to be
! smoothed, (i. e. y(k) is the functional value
! corresponding to x(k) ).
! 
! p is the period (p .gt. x(n)-x(1)).
! 
! d is a parameter containing the observation weights.
! this may either be an array of length n or a scalar
! (interpreted as a constant). the value of d
! corresponding to the observation (x(k),y(k)) should
! be an approximation to the standard deviation of error.
! 
! isw contains a switch indicating whether the parameter
! d is to be considered a vector or a scalar,
! = 0 if d is an array of length n,
! = 1 if d is a scalar.
! 
! s contains the value controlling the smoothing. this
! must be non-negative. for s equal to zero, the
! subroutine does interpolation, larger values lead to
! smoother funtions. if parameter d contains standard
! deviation estimates, a reasonable value for s is
! float(n).
! 
! eps contains a tolerance on the relative precision to
! which s is to be interpreted. this must be greater than
! or equal to zero and less than equal or equal to one. a
! reasonable value for eps is sqrt(2./float(n)).
! 
! ys is an array of length at least n.
! 
! ysp is an array of length at least n.
! 
! sigma contains the tension factor. this value indicates
! the degree to which the first derivative part of the
! smoothing functional is emphasized. if sigma is nearly
! zero (e. g. .001) the resulting curve is approximately a
! cubic spline. if sigma is large (e. g. 50.) the
! resulting curve is nearly a polygonal line. if sigma
! equals zero a cubic spline results. a standard value for
! sigma is approximately 1.
! 
! and
! 
! td, tsd1, hd, hsd1, hsd2, rd, rsd1, rsd2, rnm1, rn, and
! v are arrays of length at least n which are used for
! scratch storage.
! 
! on output--
! 
! ys contains the smoothed ordinate values.
! 
! ysp contains the values of the second derivative of the
! smoothed curve at the given nodes.
! 
! ierr contains an error flag,
! = 0 for normal return,
! = 1 if n is less than 2,
! = 2 if s is negative,
! = 3 if eps is negative or greater than one,
! = 4 if x-values are not strictly increasing,
! = 5 if a d-value is non-positive,
! = 6 if p is less than or equal to x(n)-x(1).
! 
! and
! 
! n, x, y, d, isw, s, eps, and sigma are unaltered.
! 
! this subroutine references package modules terms and
! snhcsh.
! 
! -----------------------------------------------------------
! 

    IF (n .LT. 2) go to 25
    IF (s .LT. 0._real_8) go to 26
    IF (eps .LT. 0._real_8 .OR. eps .GT. 1._real_8) go to 27
    IF (p .LE. x(n)-x(1)) go to 30
    ierr = 0
    q = 0._real_8
    rsd1(1) = 0._real_8
    rsd2(1) = 0._real_8
    rsd2(2) = 0._real_8
    rsd1(n-1) = 0._real_8
    rsd2(n-1) = 0._real_8
    rsd2(n) = 0._real_8
    ! 
    ! denormalize tension factor
    ! 
    sigmap = ABS(sigma)*REAL(n,kind=real_8)/p
    ! 
    ! form t matrix and second differences of y into ys
    ! 
    nm1 = n-1
    nm2 = n-2
    nm3 = n-3
    delxi1 = x(1)+p-x(n)
    delyi1 = (y(1)-y(n))/delxi1
    CALL terms (dim1,tsd1(1),sigmap,delxi1)
    hsd1(1) = 1._real_8/delxi1
    DO 1 i = 1,n
       ip1 = i+1
       IF (i .EQ. n) ip1 = 1
       delxi = x(ip1)-x(i)
       IF (i .EQ. n) delxi = x(1)+p-x(n)
       IF (delxi .LE. 0._real_8) go to 28
       delyi = (y(ip1)-y(i))/delxi
       ys(i) = delyi-delyi1
       CALL terms (di,tsd1(ip1),sigmap,delxi)
       td(i) = di+dim1
       hd(i) = -(1._real_8/delxi+1._real_8/delxi1)
       hsd1(ip1) = 1._real_8/delxi
       delxi1 = delxi
       delyi1 = delyi
1   dim1 = di
    hsd11 = hsd1(1)
    IF (n .GE. 3) go to 2
    tsd1(2) = tsd1(1)+tsd1(2)
    tsd1(1) = 0._real_8
    hsd1(2) = hsd1(1)+hsd1(2)
    hsd1(1) = 0._real_8
    ! 
    ! calculate lower and upper tolerances
    ! 
2   sl = s*(1._real_8-eps)
    su = s*(1._real_8+eps)
    IF (d(1) .LE. 0._real_8) go to 29
    IF (isw .EQ. 1) go to 5
    ! 
    ! form h matrix - d array
    ! 
    betapp = hsd1(n)*d(n)*d(n)
    betap = hsd1(1)*d(1)*d(1)
    alphap = hd(n)*d(n)*d(n)
    im1 = n
    sumd = 0._real_8
    sumy = 0._real_8
    DO 3 i = 1,n
       disq = d(i)*d(i)
       sumd = sumd+1._real_8/disq
       sumy = sumy+y(i)/disq
       ip1 = i+1
       IF (i .EQ. n) ip1 = 1
       alpha = hd(i)*disq
       IF (d(ip1) .LE. 0._real_8) go to 29
       hsd1ip = hsd1(ip1)
       IF (i .EQ. n) hsd1ip = hsd11
       beta = hsd1ip*d(ip1)*d(ip1)
       hd(i) = (hsd1(i)*d(im1))**2+alpha*hd(i)&
            +beta*hsd1ip
       hsd2(i) = hsd1(i)*betapp
       hsd1(i) = hsd1(i)*(alpha+alphap)
       im1 = i
       alphap = alpha
       betapp = betap
3   betap = beta
    IF (n .EQ. 3) hsd1(3) = hsd1(3)+hsd2(2)
    ! 
    ! test for straight line fit
    ! 
    con = sumy/sumd
    sum = 0._real_8
    DO 4 i = 1,n
4   sum = sum+((y(i)-con)/d(i))**2
    IF (sum .LE. su) go to 23
    go to 8
    ! 
    ! form h matrix - d constant
    ! 
5   sl = d(1)*d(1)*sl
    su = d(1)*d(1)*su
    hsd1p = hsd1(n)
    hdim1 = hd(n)
    sumy = 0._real_8
    DO 6 i = 1,n
       sumy = sumy+y(i)
       hsd1ip = hsd11
       IF (i .LT. n) hsd1ip = hsd1(i+1)
       hdi = hd(i)
       hd(i) = hsd1(i)*hsd1(i)+hdi*hdi+hsd1ip*hsd1ip
       hsd2(i) = hsd1(i)*hsd1p
       hsd1p = hsd1(i)
       hsd1(i) = hsd1p*(hdi+hdim1)
6   hdim1 = hdi
    IF (n .EQ. 3) hsd1(3) = hsd1(3)+hsd2(2)
    ! 
    ! test for straight line fit
    ! 
    con = sumy/REAL(n,kind=real_8)
    sum = 0._real_8
    DO 7 i = 1,n
7   sum = sum+(y(i)-con)**2
    IF (sum .LE. su) go to 23
    ! 
    ! top of iteration
    ! cholesky factorization of q*t+h into r
    ! 
    ! 
    ! i = 1
    ! 
8   rd(1) = 1._real_8/(q*td(1)+hd(1))
    rnm1(1) = hsd2(1)
    yspnm1 = ys(nm1)
    rn(1) = q*tsd1(1)+hsd1(1)
    yspn = ys(n)
    ysp(1) = ys(1)
    rsd1i = q*tsd1(2)+hsd1(2)
    rsd1(2) = rsd1i*rd(1)
    sumnm1 = 0._real_8
    sum2 = 0._real_8
    sumn = 0._real_8
    IF (n .EQ. 3) go to 11
    IF (n .EQ. 2) go to 12
    ! 
    ! i = 2
    ! 
    rd(2) = 1._real_8/(q*td(2)+hd(2)-rsd1i*rsd1(2))
    rnm1(2) = -rnm1(1)*rsd1(2)
    rn(2) = hsd2(2)-rn(1)*rsd1(2)
    ysp(2) = ys(2)-rsd1(2)*ysp(1)
    IF (n .EQ. 4) go to 10
    DO 9 i = 3,nm2
       rsd2i = hsd2(i)
       rsd1i = q*tsd1(i)+hsd1(i)-rsd2i*rsd1(i-1)
       rsd2(i) = rsd2i*rd(i-2)
       rsd1(i) = rsd1i*rd(i-1)
       rd(i) = 1._real_8/(q*td(i)+hd(i)-rsd1i*rsd1(i)&
            -rsd2i*rsd2(i))
       rnm1(i) = -rnm1(i-2)*rsd2(i)-rnm1(i-1)*rsd1(i)
       rnm1t = rnm1(i-2)*rd(i-2)
       sumnm1 = sumnm1+rnm1t*rnm1(i-2)
       rnm1(i-2) = rnm1t
       sum2 = sum2+rnm1t*rn(i-2)
       yspnm1 = yspnm1-rnm1t*ysp(i-2)
       rn(i) = -rn(i-2)*rsd2(i)-rn(i-1)*rsd1(i)
       rnt = rn(i-2)*rd(i-2)
       sumn = sumn+rnt*rn(i-2)
       rn(i-2) = rnt
       yspn = yspn-rnt*ysp(i-2)
9   ysp(i) = ys(i)-rsd1(i)*ysp(i-1)-rsd2(i)*ysp(i-2)
    ! 
    ! i = n-3
    ! 
10  rnm1(nm3) = hsd2(nm1)+rnm1(nm3)
    rnm1(nm2) = rnm1(nm2)-hsd2(nm1)*rsd1(nm2)
    rnm1t = rnm1(nm3)*rd(nm3)
    sumnm1 = sumnm1+rnm1t*rnm1(nm3)
    rnm1(nm3) = rnm1t
    sum2 = sum2+rnm1t*rn(nm3)
    yspnm1 = yspnm1-rnm1t*ysp(nm3)
    rnt = rn(nm3)*rd(nm3)
    sumn = sumn+rnt*rn(nm3)
    rn(nm3) = rnt
    yspn = yspn-rnt*ysp(nm3)
    ! 
    ! i = n-2
    ! 
11  rnm1(nm2) = q*tsd1(nm1)+hsd1(nm1)+rnm1(nm2)
    rnm1t = rnm1(nm2)*rd(nm2)
    sumnm1 = sumnm1+rnm1t*rnm1(nm2)
    rnm1(nm2) = rnm1t
    rn(nm2) = hsd2(n)+rn(nm2)
    sum2 = sum2+rnm1t*rn(nm2)
    yspnm1 = yspnm1-rnm1t*ysp(nm2)
    rnt = rn(nm2)*rd(nm2)
    sumn = sumn+rnt*rn(nm2)
    rn(nm2) = rnt
    yspn = yspn-rnt*ysp(nm2)
    ! 
    ! i = n-1
    ! 
12  rd(nm1) = 1._real_8/(q*td(nm1)+hd(nm1)-sumnm1)
    ysp(nm1) = yspnm1
    rn(nm1) = q*tsd1(n)+hsd1(n)-sum2
    rnt = rn(nm1)*rd(nm1)
    sumn = sumn+rnt*rn(nm1)
    rn(nm1) = rnt
    yspn = yspn-rnt*ysp(nm1)
    ! 
    ! i = n
    ! 
    rdn = q*td(n)+hd(n)-sumn
    rd(n) = 0._real_8
    IF (rdn .GT. 0._real_8) rd(n) = 1._real_8/rdn
    ysp(n) = yspn
    ! 
    ! back solve of r(transpose)* r * ysp = ys
    ! 
    ysp(n) = rd(n)*ysp(n)
    ysp(nm1) = rd(nm1)*ysp(nm1)-rn(nm1)*ysp(n)
    IF (n .EQ. 2) go to 14
    yspn = ysp(n)
    yspnm1 = ysp(nm1)
    DO 13 ibak = 1,nm2
       i = nm1-ibak
13  ysp(i) = rd(i)*ysp(i)-rsd1(i+1)*ysp(i+1)&
         -rsd2(i+2)*ysp(i+2)-rnm1(i)*yspnm1&
         -rn(i)*yspn
14  sum = 0._real_8
    delyi1 = (ysp(1)-ysp(n))/(x(1)+p-x(n))
    IF (isw .EQ. 1) go to 16
    ! 
    ! calculation of residual norm
    ! - d array
    ! 
    DO 15 i = 1,nm1
       delyi = (ysp(i+1)-ysp(i))/(x(i+1)-x(i))
       v(i) = (delyi-delyi1)*d(i)*d(i)
       sum = sum+v(i)*(delyi-delyi1)
15  delyi1 = delyi
    delyi = (ysp(1)-ysp(n))/(x(1)+p-x(n))
    v(n) = (delyi-delyi1)*d(n)*d(n)
    go to 18
    ! 
    ! calculation of residual norm
    ! - d constant
    ! 
16  DO 17 i = 1,nm1
       delyi = (ysp(i+1)-ysp(i))/(x(i+1)-x(i))
       v(i) = delyi-delyi1
       sum = sum+v(i)*(delyi-delyi1)
17  delyi1 = delyi
    delyi = (ysp(1)-ysp(n))/(x(1)+p-x(n))
    v(n) = delyi-delyi1
18  sum = sum+v(n)*(delyi-delyi1)
    ! 
    ! test for convergence
    ! 
    IF (sum .LE. su .AND. sum .GE. sl .AND. q .GT. 0._real_8)&
         go to 21
    ! 
    ! calculation of newton correction
    ! 
    f = 0._real_8
    g = 0._real_8
    rnm1sm = 0._real_8
    rnsm = 0._real_8
    im1 = n
    IF (n .EQ. 2) go to 20
    wim2 = 0._real_8
    wim1 = 0._real_8
    DO 19 i = 1,nm2
       tui = tsd1(i)*ysp(im1)+td(i)*ysp(i)&
            +tsd1(i+1)*ysp(i+1)
       wi = tui-rsd1(i)*wim1-rsd2(i)*wim2
       rnm1sm = rnm1sm-rnm1(i)*wi
       rnsm = rnsm-rn(i)*wi
       f = f+tui*ysp(i)
       g = g+wi*wi*rd(i)
       im1 = i
       wim2 = wim1
19  wim1 = wi
20  tui = tsd1(nm1)*ysp(im1)+td(nm1)*ysp(nm1)&
         +tsd1(n)*ysp(n)
    wi = tui+rnm1sm
    f = f+tui*ysp(nm1)
    g = g+wi*wi*rd(nm1)
    tui = tsd1(n)*ysp(nm1)+td(n)*ysp(n)&
         +tsd1(1)*ysp(1)
    wi = tui+rnsm-rn(nm1)*wi
    f = f+tui*ysp(n)
    g = g+wi*wi*rd(n)
    h = f-q*g
    IF (h .LE. 0._real_8 .AND. q .GT. 0._real_8) go to 21
    ! 
    ! update q - newton step
    ! 
    step = (sum-SQRT(sum*sl))/h
    IF (sl .NE. 0._real_8) step = step*SQRT(sum/sl)
    q = q+step
    go to 8
    ! 
    ! store smoothed y-values and second derivatives
    ! 
21  DO 22 i = 1,n
       ys(i) = y(i)-v(i)
22  ysp(i) = q*ysp(i)
    RETURN
    ! 
    ! store constant ys and zero ysp
    ! 
23  DO 24 i = 1,n
       ys(i) = con
24  ysp(i) = 0._real_8
    RETURN
    ! 
    ! n less than 2
    ! 
25  ierr = 1
    RETURN
    ! 
    ! s negative
    ! 
26  ierr = 2
    RETURN
    ! 
    ! eps negative or greater than 1
    ! 
27  ierr = 3
    RETURN
    ! 
    ! x-values not strictly increasing
    ! 
28  ierr = 4
    RETURN
    ! 
    ! weight non-positive
    ! 
29  ierr = 5
    RETURN
    ! 
    ! incorrect period
    ! 
30  ierr = 6
    RETURN
  END SUBROUTINE curvpp

  SUBROUTINE curvss (n,x,y,d,isw,s,eps,ys,ysp,sigma,td,&
       tsd1,hd,hsd1,hsd2,rd,rsd1,rsd2,v,&
       ierr)
    INTEGER                                  :: n
    REAL(real_8)                             :: x(n), y(n), d(n)
    INTEGER                                  :: isw
    REAL(real_8)                             :: s, eps, ys(n), ysp(n), sigma, &
                                                td(n), tsd1(n), hd(n), &
                                                hsd1(n), hsd2(n), rd(n), &
                                                rsd1(n), rsd2(n), v(n)
    INTEGER                                  :: ierr

    INTEGER                                  :: i, ibak, nm1, nm3
    REAL(real_8) :: alpha, alphap, beta, betap, betapp, delxi, delxi1, delyi, &
      delyi1, di, dim1, f, g, h, hdi, hdim1, hsd1p, p, rdim1, rsd1i, rsd2i, &
      sigmap, sl, step, su, sum, tui, wi, wim1, wim2, yspim2

! 
! coded by alan kaylor cline
! from fitpack -- january 26, 1987
! a curve and surface fitting package
! a product of pleasant valley software
! 8603 altus cove, austin, texas 78759, usa
! 
! this subroutine determines the parameters necessary to
! compute a smoothing spline under tension. for a given
! increasing sequence of abscissae (x(i)), i = 1,..., n and
! associated ordinates (y(i)), i = 1,..., n, the function
! determined minimizes the summation from i = 1 to n-1 of
! the square of the second derivative of f plus sigma
! squared times the difference of the first derivative of f
! and (f(x(i+1))-f(x(i)))/(x(i+1)-x(i)) squared, over all
! functions f with two continuous derivatives such that the
! summation of the square of (f(x(i))-y(i))/d(i) is less
! than or equal to a given constant s, where (d(i)), i = 1,
! ..., n are a given set of observation weights. the
! function determined is a spline under tension with third
! derivative discontinuities at (x(i)), i = 2,..., n-1. for
! actual computation of points on the curve it is necessary
! to call the function curv2.
! 
! on input--
! 
! n is the number of values to be smoothed (n.ge.2).
! 
! x is an array of the n increasing abscissae of the
! values to be smoothed.
! 
! y is an array of the n ordinates of the values to be
! smoothed, (i. e. y(k) is the functional value
! corresponding to x(k) ).
! 
! d is a parameter containing the observation weights.
! this may either be an array of length n or a scalar
! (interpreted as a constant). the value of d
! corresponding to the observation (x(k),y(k)) should
! be an approximation to the standard deviation of error.
! 
! isw contains a switch indicating whether the parameter
! d is to be considered a vector or a scalar,
! = 0 if d is an array of length n,
! = 1 if d is a scalar.
! 
! s contains the value controlling the smoothing. this
! must be non-negative. for s equal to zero, the
! subroutine does interpolation, larger values lead to
! smoother funtions. if parameter d contains standard
! deviation estimates, a reasonable value for s is
! float(n).
! 
! eps contains a tolerance on the relative precision to
! which s is to be interpreted. this must be greater than
! or equal to zero and less than equal or equal to one. a
! reasonable value for eps is sqrt(2./float(n)).
! 
! ys is an array of length at least n.
! 
! ysp is an array of length at least n.
! 
! sigma contains the tension factor. this value indicates
! the degree to which the first derivative part of the
! smoothing functional is emphasized. if sigma is nearly
! zero (e. g. .001) the resulting curve is approximately a
! cubic spline. if sigma is large (e. g. 50.) the
! resulting curve is nearly a polygonal line. if sigma
! equals zero a cubic spline results. a standard value for
! sigma is approximately 1.
! 
! and
! 
! td, tsd1, hd, hsd1, hsd2, rd, rsd1, rsd2, and v are
! arrays of length at least n which are used for scratch
! storage.
! 
! on output--
! 
! ys contains the smoothed ordinate values.
! 
! ysp contains the values of the second derivative of the
! smoothed curve at the given nodes.
! 
! ierr contains an error flag,
! = 0 for normal return,
! = 1 if n is less than 2,
! = 2 if s is negative,
! = 3 if eps is negative or greater than one,
! = 4 if x-values are not strictly increasing,
! = 5 if a d-value is non-positive.
! 
! and
! 
! n, x, y, d, isw, s, eps, and sigma are unaltered.
! 
! this subroutine references package modules terms and
! snhcsh.
! 
! -----------------------------------------------------------
! 

    IF (n .LT. 2) go to 16
    IF (s .LT. 0._real_8) go to 17
    IF (eps .LT. 0._real_8 .OR. eps .GT. 1._real_8) go to 18
    ierr = 0
    p = 0._real_8
    v(1) = 0._real_8
    v(n) = 0._real_8
    ysp(1) = 0._real_8
    ysp(n) = 0._real_8
    IF (n .EQ. 2) go to 14
    rsd1(1) = 0._real_8
    rd(1) = 0._real_8
    rsd2(n) = 0._real_8
    rdim1 = 0._real_8
    yspim2 = 0._real_8
    ! 
    ! denormalize tension factor
    ! 
    sigmap = ABS(sigma)*REAL(n-1,kind=real_8)/(x(n)-x(1))
    ! 
    ! form t matrix and second differences of y into ys
    ! 
    nm1 = n-1
    nm3 = n-3
    delxi1 = 1._real_8
    delyi1 = 0._real_8
    dim1 = 0._real_8
    DO 1 i = 1,nm1
       delxi = x(i+1)-x(i)
       IF (delxi .LE. 0._real_8) go to 19
       delyi = (y(i+1)-y(i))/delxi
       ys(i) = delyi-delyi1
       CALL terms (di,tsd1(i+1),sigmap,delxi)
       td(i) = di+dim1
       hd(i) = -(1._real_8/delxi+1._real_8/delxi1)
       hsd1(i+1) = 1._real_8/delxi
       delxi1 = delxi
       delyi1 = delyi
1   dim1 = di
    ! 
    ! calculate lower and upper tolerances
    ! 
    sl = s*(1._real_8-eps)
    su = s*(1._real_8+eps)
    IF (isw .EQ. 1) go to 3
    ! 
    ! form h matrix - d array
    ! 
    IF (d(1) .LE. 0._real_8 .OR. d(2) .LE. 0._real_8) go to 20
    betapp = 0._real_8
    betap = 0._real_8
    alphap = 0._real_8
    DO 2 i = 2,nm1
       alpha = hd(i)*d(i)*d(i)
       IF (d(i+1) .LE. 0._real_8) go to 20
       beta = hsd1(i+1)*d(i+1)*d(i+1)
       hd(i) = (hsd1(i)*d(i-1))**2+alpha*hd(i)&
            +beta*hsd1(i+1)
       hsd2(i) = hsd1(i)*betapp
       hsd1(i) = hsd1(i)*(alpha+alphap)
       alphap = alpha
       betapp = betap
2   betap = beta
    go to 5
    ! 
    ! form h matrix - d constant
    ! 
3   IF (d(1) .LE. 0._real_8) go to 20
    sl = d(1)*d(1)*sl
    su = d(1)*d(1)*su
    hsd1p = 0._real_8
    hdim1 = 0._real_8
    DO 4 i = 2,nm1
       hdi = hd(i)
       hd(i) = hsd1(i)*hsd1(i)+hdi*hdi+hsd1(i+1)*hsd1(i+1)
       hsd2(i) = hsd1(i)*hsd1p
       hsd1p = hsd1(i)
       hsd1(i) = hsd1p*(hdi+hdim1)
4   hdim1 = hdi
    ! 
    ! top of iteration
    ! cholesky factorization of p*t+h into r
    ! 
5   DO 6 i = 2,nm1
       rsd2i = hsd2(i)
       rsd1i = p*tsd1(i)+hsd1(i)-rsd2i*rsd1(i-1)
       rsd2(i) = rsd2i*rdim1
       rdim1 = rd(i-1)
       rsd1(i) = rsd1i*rdim1
       rd(i) = 1._real_8/(p*td(i)+hd(i)-rsd1i*rsd1(i)&
            -rsd2i*rsd2(i))
       ysp(i) = ys(i)-rsd1(i)*ysp(i-1)-rsd2(i)*yspim2
6   yspim2 = ysp(i-1)
    ! 
    ! back solve of r(transpose)* r * ysp = ys
    ! 
    ysp(nm1) = rd(nm1)*ysp(nm1)
    IF (n .EQ. 3) go to 8
    DO 7 ibak = 1,nm3
       i = nm1-ibak
7   ysp(i) = rd(i)*ysp(i)-rsd1(i+1)*ysp(i+1)&
         -rsd2(i+2)*ysp(i+2)
8   sum = 0._real_8
    delyi1 = 0._real_8
    IF (isw .EQ. 1) go to 10
    ! 
    ! calculation of residual norm
    ! - d array
    ! 
    DO 9 i = 1,nm1
       delyi = (ysp(i+1)-ysp(i))/(x(i+1)-x(i))
       v(i) = (delyi-delyi1)*d(i)*d(i)
       sum = sum+v(i)*(delyi-delyi1)
9   delyi1 = delyi
    v(n) = -delyi1*d(n)*d(n)
    go to 12
    ! 
    ! calculation of residual norm
    ! - d constant
    ! 
10  DO 11 i = 1,nm1
       delyi = (ysp(i+1)-ysp(i))/(x(i+1)-x(i))
       v(i) = delyi-delyi1
       sum = sum+v(i)*(delyi-delyi1)
11  delyi1 = delyi
    v(n) = -delyi1
12  sum = sum-v(n)*delyi1
    ! 
    ! test for convergence
    ! 
    IF (sum .LE. su) go to 14
    ! 
    ! calculation of newton correction
    ! 
    f = 0._real_8
    g = 0._real_8
    wim2 = 0._real_8
    wim1 = 0._real_8
    DO 13 i = 2,nm1
       tui = tsd1(i)*ysp(i-1)+td(i)*ysp(i)&
            +tsd1(i+1)*ysp(i+1)
       wi = tui-rsd1(i)*wim1-rsd2(i)*wim2
       f = f+tui*ysp(i)
       g = g+wi*wi*rd(i)
       wim2 = wim1
13  wim1 = wi
    h = f-p*g
    IF (h .LE. 0._real_8) go to 14
    ! 
    ! update p - newton step
    ! 
    step = (sum-SQRT(sum*sl))/h
    IF (sl .NE. 0._real_8) step = step*SQRT(sum/sl)
    p = p+step
    go to 5
    ! 
    ! store smoothed y-values and second derivatives
    ! 
14  DO 15 i = 1,n
       ys(i) = y(i)-v(i)
15  ysp(i) = p*ysp(i)
    RETURN
    ! 
    ! n less than 2
    ! 
16  ierr = 1
    RETURN
    ! 
    ! s negative
    ! 
17  ierr = 2
    RETURN
    ! 
    ! eps negative or greater than 1
    ! 
18  ierr = 3
    RETURN
    ! 
    ! x-values not strictly increasing
    ! 
19  ierr = 4
    RETURN
    ! 
    ! weight non-positive
    ! 
20  ierr = 5
    RETURN
  END SUBROUTINE curvss

  REAL(real_8)  FUNCTION intrvl (t,x,n)
    INTEGER :: n
    REAL(real_8) :: t,x(n)
    REAL(real_8) :: tt
    INTEGER :: i=1,ih,il

    ! 
    ! coded by alan kaylor cline
    ! from fitpack -- january 26, 1987
    ! a curve and surface fitting package
    ! a product of pleasant valley software
    ! 8603 altus cove, austin, texas 78759, usa
    ! 
    ! this function determines the index of the interval
    ! (determined by a given increasing sequence) in which
    ! a given value lies.
    ! 
    ! on input--
    ! 
    ! t is the given value.
    ! 
    ! x is a vector of strictly increasing values.
    ! 
    ! and
    ! 
    ! n is the length of x (n .ge. 2).
    ! 
    ! on output--
    ! 
    ! intrvl returns an integer i such that
    ! 
    ! i =  1       if         e   t .lt. x(2)  ,
    ! i =  n-1     if x(n-1) .le. t            ,
    ! otherwise       x(i)  .le. t .le. x(i+1),
    ! 
    ! none of the input parameters are altered.
    ! 
    ! -----------------------------------------------------------
    ! 
!!!   SAVE i
    ! 
    tt = t
    ! 
    ! check for illegal i
    ! 
    IF (i .GE. n) i = n/2
    ! 
    ! check old interval and extremes
    ! 
    IF (tt .LT. x(i)) THEN
       IF (tt .LE. x(2)) THEN
          i = 1
          intrvl = 1
          RETURN
       ELSE
          il = 2
          ih = i
       ENDIF
    ELSE IF (tt .LE. x(i+1)) THEN
       intrvl = i
       RETURN
    ELSE IF (tt .GE. x(n-1)) THEN
       i = n-1
       intrvl = n-1
       RETURN
    ELSE
       il = i+1
       ih = n-1
    ENDIF
    ! 
    ! binary search loop
    ! 
1   i = (il+ih)/2
    IF (tt .LT. x(i)) THEN
       ih = i
    ELSE IF (tt .GT. x(i+1)) THEN
       il = i+1
    ELSE
       intrvl = i
       RETURN
    ENDIF
    go to 1
  END FUNCTION intrvl

  REAL(real_8)  FUNCTION intrvp (t,x,n,p,tp)
    INTEGER :: n
    REAL(real_8) :: t,x(n),p,tp
    REAL(real_8) :: tt
    INTEGER :: ih,il,nper
    ! 
    ! coded by alan kaylor cline
    ! from fitpack -- january 26, 1987
    ! a curve and surface fitting package
    ! a product of pleasant valley software
    ! 8603 altus cove, austin, texas 78759, usa
    ! 
    ! this function determines the index of the interval
    ! (determined by a given increasing sequence) in which a
    ! given value lies, after translating the value to within
    ! the correct period.  it also returns this translated value.
    ! 
    ! on input--
    ! 
    ! t is the given value.
    ! 
    ! x is a vector of strictly increasing values.
    ! 
    ! n is the length of x (n .ge. 2).
    ! 
    ! and
    ! 
    ! p contains the period.
    ! 
    ! on output--
    ! 
    ! tp contains a translated value of t (i. e. x(1) .le. tp,
    ! tp .lt. x(1)+p, and tp = t + k*p for some integer k).
    ! 
    ! intrvl returns an integer i such that
    ! 
    ! i = 1       if             tp .lt. x(2)  ,
    ! i = n       if   x(n) .le. tp            ,
    ! otherwise       x(i)  .le. tp .lt. x(i+1),
    ! 
    ! none of the input parameters are altered.
    ! 
    ! -----------------------------------------------------------
    ! 
#if (! defined __NEC)
    INTEGER :: i
    COMMON /fitpack/ i
#else
!!!   SAVE i
    INTEGER :: i=1
#endif
    ! 
    nper = (t-x(1))/p
    tp = t-REAL(nper,kind=real_8)*p
    IF (tp .LT. x(1)) tp = tp+p
    tt = tp
    ! 
    ! check for illegal i
    ! 
    IF (i .GE. n) i = n/2
    ! 
    ! check old interval and extremes
    ! 
    IF (tt .LT. x(i)) THEN
       IF (tt .LE. x(2)) THEN
          i = 1
          intrvp = 1
          RETURN
       ELSE
          il = 2
          ih = i
       ENDIF
    ELSE IF (tt .LE. x(i+1)) THEN
       intrvp = i
       RETURN
    ELSE IF (tt .GE. x(n)) THEN
       i = n
       intrvp = n
       RETURN
    ELSE
       il = i+1
       ih = n
    ENDIF
    ! 
    ! binary search loop
    ! 
1   i = (il+ih)/2
    IF (tt .LT. x(i)) THEN
       ih = i
    ELSE IF (tt .GT. x(i+1)) THEN
       il = i+1
    ELSE
       intrvp = i
       RETURN
    ENDIF
    go to 1
  END FUNCTION intrvp

  SUBROUTINE x_snhcsh (sinhm,coshm,x,isw)
    REAL(real_8)                             :: sinhm, coshm, x
    INTEGER                                  :: isw

    REAL(real_8), PARAMETER :: cp0 = 0.5000000_real_8 , &
      cp1 = 0.4166665e-1_real_8 , cp2 = 0.1388967e-2_real_8 , &
      cp3 = 0.2472673e-4_real_8 , cp4 = 0.2982628e-6_real_8 , &
      sp10 = 0.1666665_real_8 , sp11 = 0.8334261e-2_real_8 , &
      sp12 = 0.1975135e-3_real_8 , sp13 = 0.3029390e-5_real_8 , &
      sp20 = 0.1667035_real_8 , sp21 = 0.8315072e-2_real_8 , &
      sp22 = 0.2018107e-3_real_8 , sp23 = 0.2459974e-5_real_8 , &
      sp24 = 0.3693467e-7_real_8 , sp31 = 0.4001477e-1_real_8 , &
      sp32 = 0.6646307e-3_real_8 , sp33 = 0.6666558e-5_real_8 , &
      sp41 = 0.9868757e-1_real_8 , sp42 = 0.2729702e-3_real_8 , &
      sp43 = 0.2311816e-4_real_8 
    REAL(real_8), PARAMETER :: sq30 = 0.6017497e1_real_8 , &
      sq31 = -0.6372739e-1_real_8 , sq32 = 0.2037930e-3_real_8 , &
      sq40 = 0.9110034e1_real_8 , sq41 = -0.7549779e-1_real_8 , &
      sq42 = 0.1776637e-3_real_8 

    REAL(real_8)                             :: ax, expx, xs

! 
! coded by alan kaylor cline
! from fitpack -- january 26, 1987
! a curve and surface fitting package
! a product of pleasant valley software
! 8603 altus cove, austin, texas 78759, usa
! 
! this subroutine returns approximations to
! sinhm(x) = sinh(x)/x-1
! coshm(x) = cosh(x)-1
! and
! coshmm(x) = (cosh(x)-1-x*x/2)/(x*x)
! with relative error less than 1.0e-6
! 
! on input--
! 
! x contains the value of the independent variable.
! 
! isw indicates the function desired
! = -1 if only sinhm is desired,
! =  0 if both sinhm and coshm are desired,
! =  1 if only coshm is desired,
! =  2 if only coshmm is desired,
! =  3 if both sinhm and coshmm are desired.
! 
! on output--
! 
! sinhm contains the value of sinhm(x) if isw .le. 0 or
! isw .eq. 3 (sinhm is unaltered if isw .eq.1 or isw .eq.
! 2).
! 
! coshm contains the value of coshm(x) if isw .eq. 0 or
! isw .eq. 1 and contains the value of coshmm(x) if isw
! .ge. 2 (coshm is unaltered if isw .eq. -1).
! 
! and
! 
! x and isw are unaltered.
! 
! -----------------------------------------------------------
! 
! 

    ax = ABS(x)
    IF (isw .GE. 0) go to 5
    ! 
    ! sinhm approximation
    ! 
    IF (ax .GT. 4.45_real_8) go to 2
    xs = ax*ax
    IF (ax .GT. 2.3_real_8) go to 1
    ! 
    ! sinhm approximation on (0.,2.3)
    ! 
    sinhm = xs*(((sp13*xs+sp12)*xs+sp11)*xs+sp10)
    RETURN
    ! 
    ! sinhm approximation on (2.3,4.45)
    ! 
1   sinhm = xs*((((sp24*xs+sp23)*xs+sp22)*xs+sp21)&
         *xs+sp20)
    RETURN
2   IF (ax .GT. 7.65_real_8) go to 3
    ! 
    ! sinhm approximation on (4.45,7.65)
    ! 
    xs = ax*ax
    sinhm = xs*(((sp33*xs+sp32)*xs+sp31)*xs+1._real_8)/&
         ((sq32*xs+sq31)*xs+sq30)
    RETURN
3   IF (ax .GT. 10.1_real_8) go to 4
    ! 
    ! sinhm approximation on (7.65,10.1)
    ! 
    xs = ax*ax
    sinhm = xs*(((sp43*xs+sp42)*xs+sp41)*xs+1._real_8)/&
         ((sq42*xs+sq41)*xs+sq40)
    RETURN
    ! 
    ! sinhm approximation above 10.1
    ! 
4   sinhm = EXP(ax)/(ax+ax)-1._real_8
    RETURN
    ! 
    ! coshm and (possibly) sinhm approximation
    ! 
5   IF (isw .GE. 2) go to 7
    IF (ax .GT. 2.3_real_8) go to 6
    xs = ax*ax
    coshm = xs*((((cp4*xs+cp3)*xs+cp2)*xs+cp1)*xs+cp0)
    IF (isw .EQ. 0) sinhm = xs*(((sp13*xs+sp12)*xs+sp11)&
         *xs+sp10)
    RETURN
6   expx = EXP(ax)
    coshm = (expx+1._real_8/expx)/2._real_8-1._real_8
    IF (isw .EQ. 0) sinhm = (expx-1._real_8/expx)/(ax+ax)-1._real_8
    RETURN
    ! 
    ! coshmm and (possibly) sinhm approximation
    ! 
7   xs = ax*ax
    IF (ax .GT. 2.3_real_8) go to 8
    coshm = xs*(((cp4*xs+cp3)*xs+cp2)*xs+cp1)
    IF (isw .EQ. 3) sinhm = xs*(((sp13*xs+sp12)*xs+sp11)&
         *xs+sp10)
    RETURN
8   expx = EXP(ax)
    coshm = ((expx+1._real_8/expx-xs)/2._real_8-1._real_8)/xs
    IF (isw .EQ. 3) sinhm = (expx-1._real_8/expx)/(ax+ax)-1._real_8
    RETURN
  END SUBROUTINE x_snhcsh

  SUBROUTINE snhcsh (sinhm,coshm,x,isw)
    REAL(real_8)                             :: sinhm, coshm, x
    INTEGER                                  :: isw

    REAL(real_8), PARAMETER :: cp1 = 0.744437205569040e-1_real_8 , &
      cp2 = 0.206270719503934e-2_real_8 , cp3 = 0.270540125846525e-4_real_8 , &
      cp4 = 0.181666923620944e-6_real_8 , cp5 = 0.552200614584744e-9_real_8 , &
      cq0 = 0.200000000000000e+1_real_8 , &
      cq1 = -0.177792255528382e-1_real_8 , &
      cq2 = 0.514609638642689e-4_real_8 , &
      sp11 = 0.398088289992973e-1_real_8 , &
      sp12 = 0.715314759211209e-3_real_8 , &
      sp13 = 0.612189863171694e-5_real_8 , &
      sp14 = 0.227581660976348e-7_real_8 , &
      sp21 = 0.425024142813226e-1_real_8 , &
      sp22 = 0.833264803327242e-3_real_8 , sp23 = 0.849213455598455e-5_real_8 
    REAL(real_8), PARAMETER :: sp24 = 0.473731823101666e-7_real_8 , &
      sp25 = 0.129094158037272e-9_real_8 , &
      sp31 = 0.428888148791777e-1_real_8 , &
      sp32 = 0.850447617691392e-3_real_8 , &
      sp33 = 0.884775635776784e-5_real_8 , &
      sp34 = 0.511529451668737e-7_real_8 , &
      sp35 = 0.155193945864942e-9_real_8 , &
      sp41 = 0.432535234960858e-1_real_8 , &
      sp42 = 0.866559391672985e-3_real_8 , &
      sp43 = 0.920119535795222e-5_real_8 , &
      sp44 = 0.545792817714192e-7_real_8 , &
      sp45 = 0.188070632058331e-9_real_8 , &
      sq10 = 0.599999999999986e+1_real_8 , &
      sq11 = -0.611470260009508e-1_real_8 , &
      sq12 = 0.206382701413725e-3_real_8 
    REAL(real_8), PARAMETER :: sq20 = 0.600000000268619e+1_real_8 , &
      sq21 = -0.449855169512505e-1_real_8 , &
      sq22 = 0.106008515744821e-3_real_8 , &
      sq30 = 0.600000145086489e+1_real_8 , &
      sq31 = 0.426677570538507e-1_real_8 , &
      sq32 = 0.933128831061610e-4_real_8 , &
      sq40 = 0.600005006283834e+1_real_8 , &
      sq41 = -0.404938841672262e-1_real_8 , &
      sq42 = 0.824891748820670e-4_real_8 , &
      zp1 = 0.244515150174258e-1_real_8 , zp2 = 0.324851059327161e-3_real_8 , &
      zp3 = 0.218274535686385e-5_real_8 , zp4 = 0.664418805876835e-8_real_8 , &
      zq0 = 0.240000000000000e+2_real_8 , zq1 = -0.213163639579425_real_8 
    REAL(real_8), PARAMETER :: zq2 = 0.616165782306621e-3_real_8 

    REAL(real_8)                             :: ax, expx, xs

! 
! coded by alan kaylor cline
! from fitpack -- january 26, 1987
! a curve and surface fitting package
! a product of pleasant valley software
! 8603 altus cove, austin, texas 78759, usa
! 
! this subroutine returns approximations to
! sinhm(x) = sinh(x)/x-1
! coshm(x) = cosh(x)-1
! and
! coshmm(x) = (cosh(x)-1-x*x/2)/(x*x)
! with relative error less than 4.0e-14.
! 
! on input--
! 
! x contains the value of the independent variable.
! 
! isw indicates the function desired
! = -1 if only sinhm is desired,
! =  0 if both sinhm and coshm are desired,
! =  1 if only coshm is desired,
! =  2 if only coshmm is desired,
! =  3 if both sinhm and coshmm are desired.
! 
! on output--
! 
! sinhm contains the value of sinhm(x) if isw .le. 0 or
! isw .eq. 3 (sinhm is unaltered if isw .eq.1 or isw .eq.
! 2).
! 
! coshm contains the value of coshm(x) if isw .eq. 0 or
! isw .eq. 1 and contains the value of coshmm(x) if isw
! .ge. 2 (coshm is unaltered if isw .eq. -1).
! 
! and
! 
! x and isw are unaltered.
! 
! -----------------------------------------------------------

    ax = ABS(x)
    IF (isw .GE. 0) go to 5
    ! 
    ! sinhm approximation
    ! 
    IF (ax .GT. 3.9_real_8) go to 2
    xs = ax*ax
    IF (ax .GT. 2.2_real_8) go to 1
    ! 
    ! sinhm approximation on (0.,2.2)
    ! 
    sinhm = xs*((((sp14*xs+sp13)*xs+sp12)*xs+sp11)*xs+1._real_8)/&
         ((sq12*xs+sq11)*xs+sq10)
    RETURN
    ! 
    ! sinhm approximation on (2.2,3.9)
    ! 
1   sinhm = xs*(((((sp25*xs+sp24)*xs+sp23)*xs+sp22)*xs+sp21)&
         *xs+1._real_8)/((sq22*xs+sq21)*xs+sq20)
    RETURN
2   IF (ax .GT. 5.1_real_8) go to 3
    ! 
    ! sinhm approximation on (3.9,5.1)
    ! 
    xs = ax*ax
    sinhm = xs*(((((sp35*xs+sp34)*xs+sp33)*xs+sp32)*xs+sp31)&
         *xs+1._real_8)/((sq32*xs+sq31)*xs+sq30)
    RETURN
3   IF (ax .GT. 6.1_real_8) go to 4
    ! 
    ! sinhm approximation on (5.1,6.1)
    ! 
    xs = ax*ax
    sinhm = xs*(((((sp45*xs+sp44)*xs+sp43)*xs+sp42)*xs+sp41)&
         *xs+1._real_8)/((sq42*xs+sq41)*xs+sq40)
    RETURN
    ! 
    ! sinhm approximation above 6.1
    ! 
4   expx = EXP(ax)
    sinhm = (expx-1._real_8/expx)/(ax+ax)-1._real_8
    RETURN
    ! 
    ! coshm and (possibly) sinhm approximation
    ! 
5   IF (isw .GE. 2) go to 7
    IF (ax .GT. 2.2_real_8) go to 6
    xs = ax*ax
    coshm = xs*(((((cp5*xs+cp4)*xs+cp3)*xs+cp2)*xs+cp1)&
         *xs+1._real_8)/((cq2*xs+cq1)*xs+cq0)
    IF (isw .EQ. 0) sinhm = xs*((((sp14*xs+sp13)*xs+sp12)&
         *xs+sp11)*xs+1._real_8)/((sq12*xs+sq11)*xs+sq10)
    RETURN
6   expx = EXP(ax)
    coshm = (expx+1._real_8/expx)/2._real_8-1._real_8
    IF (isw .EQ. 0) sinhm = (expx-1._real_8/expx)/(ax+ax)-1._real_8
    RETURN
    ! 
    ! coshmm and (possibly) sinhm approximation
    ! 
7   xs = ax*ax
    IF (ax .GT. 2.2_real_8) go to 8
    coshm = xs*((((zp4*xs+zp3)*xs+zp2)*xs+zp1)*xs+1._real_8)/&
         ((zq2*xs+zq1)*xs+zq0)
    IF (isw .EQ. 3) sinhm = xs*((((sp14*xs+sp13)*xs+sp12)&
         *xs+sp11)*xs+1._real_8)/((sq12*xs+sq11)*xs+sq10)
    RETURN
8   expx = EXP(ax)
    coshm = ((expx+1._real_8/expx-xs)/2._real_8-1._real_8)/xs
    IF (isw .EQ. 3) sinhm = (expx-1._real_8/expx)/(ax+ax)-1._real_8
    RETURN
  END SUBROUTINE snhcsh

  SUBROUTINE terms (diag,sdiag,sigma,del)
    REAL(real_8)                             :: diag, sdiag, sigma, del

    REAL(real_8)                             :: coshm, denom, sigdel, sinhm

! 
! coded by alan kaylor cline
! from fitpack -- january 26, 1987
! a curve and surface fitting package
! a product of pleasant valley software
! 8603 altus cove, austin, texas 78759, usa
! 
! this subroutine computes the diagonal and superdiagonal
! terms of the tridiagonal linear system associated with
! spline under tension interpolation.
! 
! on input--
! 
! sigma contains the tension factor.
! 
! and
! 
! del contains the step size.
! 
! on output--
! 
! sigma*del*cosh(sigma*del) - sinh(sigma*del)
! diag = del*--------------------------------------------.
! (sigma*del)**2 * sinh(sigma*del)
! 
! sinh(sigma*del) - sigma*del
! sdiag = del*----------------------------------.
! (sigma*del)**2 * sinh(sigma*del)
! 
! and
! 
! sigma and del are unaltered.
! 
! this subroutine references package module snhcsh.
! 
! -----------------------------------------------------------
! 

    IF (sigma .NE. 0._real_8) go to 1
    diag = del/3._real_8
    sdiag = del/6._real_8
    RETURN
1   sigdel = sigma*del
    CALL snhcsh (sinhm,coshm,sigdel,0)
    denom = sigma*sigdel*(1._real_8+sinhm)
    diag = (coshm-sinhm)/denom
    sdiag = sinhm/denom
    RETURN
  END SUBROUTINE terms


END MODULE fitpack_utils
