#ifdef __SR8000
!option MP(P(0)), LANGLVL(SAVE(0))
#endif
MODULE ylmr_utils
  USE bessm_utils,                     ONLY: bessl
  USE cnst,                            ONLY: fpi
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE ylmr2_utils,                     ONLY: ylmr2

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: bess
  PUBLIC :: ylmr
  PUBLIC :: ylmkr

CONTAINS

  ! ==================================================================
  FUNCTION ylmr(l,ig,gk)
    ! ==----------------------------------------------------------------==
    ! == REAL SPHERICAL HARMONICS,  L IS COMBINED INDEX                 ==
    ! == FOR lm  (L=1,2...9) ORDER:                                     ==
    ! == s, p_x, p_z, p_y, d_xy, d_xz, d_z^2, d_yz, d_x^2-y^2           ==
    ! ==----------------------------------------------------------------==
    INTEGER                                  :: l, ig
    REAL(real_8)                             :: gk(3,*), ylmr

    REAL(real_8)                             :: hg(1), r, x, y, ylm(1), z

    x = gk(1,ig)
    y = gk(2,ig)
    z = gk(3,ig)
    r = MAX(SQRT(x*x+y*y+z*z),1.0e-6_real_8)
    x = x/r
    y = y/r
    z = z/r
    IF (l.EQ.1) THEN
       ylmr = SQRT(1.0_real_8/fpi)
    ELSE IF (l.EQ.2) THEN
       ylmr = SQRT(3.0_real_8/fpi)*x
    ELSE IF (l.EQ.3) THEN
       ylmr = SQRT(3.0_real_8/fpi)*z
    ELSE IF (l.EQ.4) THEN
       ylmr = SQRT(3.0_real_8/fpi)*y
    ELSE IF (l.EQ.5) THEN
       ylmr = SQRT(15.0_real_8/fpi)*x*y
    ELSE IF (l.EQ.6) THEN
       ylmr = SQRT(15.0_real_8/fpi)*x*z
    ELSE IF (l.EQ.7) THEN
       ylmr = SQRT(5.0_real_8/fpi/4._real_8)*(2._real_8*z*z-x*x-y*y)
    ELSE IF (l.EQ.8) THEN
       ylmr = SQRT(15.0_real_8/fpi)*y*z
    ELSE IF (l.EQ.9) THEN
       ylmr = SQRT(15.0_real_8/fpi/4._real_8)*(x*x-y*y)
    ELSE IF (l.GE.10) THEN
       hg=r*r
       CALL ylmr2(l,1,hg,gk(:,ig),ylm)
       ylmr=ylm(1)
    ENDIF
    ! ==----------------------------------------------------------------==
    RETURN
  END FUNCTION ylmr
  ! ==================================================================
  SUBROUTINE bess(xg,l,mmax,r,jl)
    ! ==----------------------------------------------------------------==
    ! == CALCULATES SPHERICAL BESSEL FUNCTIONS                          ==
    ! ==----------------------------------------------------------------==
    REAL(real_8)                             :: xg
    INTEGER                                  :: l, mmax
    REAL(real_8)                             :: r(mmax), jl(mmax)

    REAL(real_8), PARAMETER                  :: eps = 1.e-8_real_8 

    INTEGER                                  :: ir
    REAL(real_8)                             :: xrg

! ==----------------------------------------------------------------==

    !$omp parallel do private(IR,XRG) 
    DO ir=1,mmax
       xrg=r(ir)*xg
       jl(ir)=bessl(l-1,xrg)
    ENDDO
    ! ==----------------------------------------------------------------==
    RETURN
  END SUBROUTINE bess
  ! ====================================================================
  FUNCTION ylmkr(l,ig,gk,rk,iplus)
    ! ==----------------------------------------------------------------==
    ! == REAL SPHERICAL HARMONICS,  L IS COMBINED INDEX FOR lm          ==
    ! == (L=1,2...9) ORDER:                                             ==
    ! == s, p_x, p_z, p_y, d_xy, d_xz, d_z^2, d_yz, d_x^2-y^2           ==
    ! ==----------------------------------------------------------------==
    INTEGER                                  :: l, ig
    REAL(real_8)                             :: gk(3,*), rk(3)
    INTEGER                                  :: iplus
    REAL(real_8)                             :: ylmkr

    REAL(real_8)                             :: gkk(3), hg(1), r, x, y, &
                                                ylm(1), z

    IF (iplus.EQ.1) THEN
       x = rk(1)+gk(1,ig)
       y = rk(2)+gk(2,ig)
       z = rk(3)+gk(3,ig)
    ELSEIF (iplus.EQ.-1) THEN
       x = rk(1)-gk(1,ig)
       y = rk(2)-gk(2,ig)
       z = rk(3)-gk(3,ig)
    ELSE
       CALL stopgm(' YLMRK','IPLUS WRONG',& 
            __LINE__,__FILE__)
    ENDIF
    gkk(1)=x
    gkk(2)=y
    gkk(3)=z
    r = MAX(SQRT(x*x+y*y+z*z),1.0e-6_real_8)
    x = x/r
    y = y/r
    z = z/r
    IF (l.EQ.1) THEN
       ylmkr = SQRT(1._real_8/fpi)
    ELSE IF (l.EQ.2) THEN
       ylmkr = SQRT(3._real_8/fpi)*x
    ELSE IF (l.EQ.3) THEN
       ylmkr = SQRT(3._real_8/fpi)*z
    ELSE IF (l.EQ.4) THEN
       ylmkr = SQRT(3._real_8/fpi)*y
    ELSE IF (l.EQ.5) THEN
       ylmkr = SQRT(15._real_8/fpi)*x*y
    ELSE IF (l.EQ.6) THEN
       ylmkr = SQRT(15._real_8/fpi)*x*z
    ELSE IF (l.EQ.7) THEN
       ylmkr = SQRT(5._real_8/fpi/4._real_8)*(2._real_8*z*z-x*x-y*y)
    ELSE IF (l.EQ.8) THEN
       ylmkr = SQRT(15._real_8/fpi)*y*z
    ELSE IF (l.EQ.9) THEN
       ylmkr = SQRT(15._real_8/fpi/4._real_8)*(x*x-y*y)
    ELSE IF (l.GE.10) THEN
       hg=r*r
       CALL ylmr2(l,1,hg,gkk(:),ylm)
       ylmkr=ylm(1)
    ENDIF
    ! ==----------------------------------------------------------------==
    RETURN
  END FUNCTION ylmkr
  ! ====================================================================

END MODULE ylmr_utils
