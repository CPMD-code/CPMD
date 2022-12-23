! ==================================================================
MODULE x_hjs
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: x_r_hjs_lc
  PUBLIC :: x_u_hjs_lc

  ! PBE parametrization
  REAL(real_8),PARAMETER,PRIVATE :: a=0.757211_real_8
  REAL(real_8),PARAMETER,PRIVATE :: b=-0.106364_real_8
  REAL(real_8),PARAMETER,PRIVATE :: c=-0.118649_real_8
  REAL(real_8),PARAMETER,PRIVATE :: d=0.609650_real_8
  REAL(real_8),PARAMETER,PRIVATE :: e=-0.0477963_real_8
  REAL(real_8),PARAMETER,PRIVATE :: a2=0.0159941_real_8
  REAL(real_8),PARAMETER,PRIVATE :: a3=0.0852995_real_8
  REAL(real_8),PARAMETER,PRIVATE :: a4=-0.160368_real_8
  REAL(real_8),PARAMETER,PRIVATE :: a5=0.152645_real_8
  REAL(real_8),PARAMETER,PRIVATE :: a6=-0.0971263_real_8
  REAL(real_8),PARAMETER,PRIVATE :: a7=0.0422061_real_8
  REAL(real_8),PARAMETER,PRIVATE :: b1=5.33319_real_8
  REAL(real_8),PARAMETER,PRIVATE :: b2=-12.4780_real_8
  REAL(real_8),PARAMETER,PRIVATE :: b3=11.0988_real_8
  REAL(real_8),PARAMETER,PRIVATE :: b4=-5.11013_real_8
  REAL(real_8),PARAMETER,PRIVATE :: b5=1.71468_real_8
  REAL(real_8),PARAMETER,PRIVATE :: b6=-0.610380_real_8
  REAL(real_8),PARAMETER,PRIVATE :: b7=0.307555_real_8
  REAL(real_8),PARAMETER,PRIVATE :: b8=-0.0770547_real_8
  REAL(real_8),PARAMETER,PRIVATE :: b9=0.0334840_real_8
  REAL(real_8),PARAMETER,PRIVATE :: s0=2.0_real_8

CONTAINS
  ! ==================================================================
  SUBROUTINE x_r_hjs_lc(r0,rho,grho,sx,v1x,v2x)
    ! ==--------------------------------------------------------------==
    REAL(real_8), INTENT(in)                 :: r0, rho, grho
    REAL(real_8), INTENT(out)                :: sx, v1x, v2x

    CHARACTER(len=*), PARAMETER              :: procedureN = 'x_r_hjs_lc'
    REAL(real_8), PARAMETER                  :: tol = 1.0e-20_real_8 

    REAL(real_8) :: Ax, dEGdn, dEGdr, dEnerdn, dEnerdr, detadn, detadr, dFdn, &
      dFdr, dFxdn, dFxdr, dkfdr, dlambadadn, dlambadadr, dsdn, dsdr, dy0dr, &
      dzetadn, dzetadr, eg, Ener, eta, f, Fx, kf, lambada, s, y0, zeta
    REAL(real_8), EXTERNAL                   :: dei

    IF (r0.LT.1.0e-3_real_8) CALL stopgm(procedureN,'RO too small',& 
         __LINE__,__FILE__)
    IF (rho > tol .AND. grho > tol) THEN
       INCLUDE 'rhjsx.inc'
       sx = Ener
       v1x = dEnerdr
       v2x = 2.0_real_8*dEnerdn ! multiply by 2, cause we take derivative vs
       ! (nabla rho) and not (nabla rho)^2
    ELSE
       sx = 0.0_real_8
       v1x = 0.0_real_8
       v2x = 0.0_real_8
    ENDIF
    ! ==--------------------------------------------------------------==
  END SUBROUTINE x_r_hjs_lc
  ! ==================================================================
  SUBROUTINE x_u_hjs_lc(r0,rhoa,rhob,grhoaa,grhobb,&
       sx,v1xa,v1xb,v2xa,v2xb,v2xab)
    ! ==--------------------------------------------------------------==
    REAL(real_8), INTENT(in)                 :: r0, rhoa, rhob, grhoaa, grhobb
    REAL(real_8), INTENT(out)                :: sx, v1xa, v1xb, v2xa, v2xb, &
                                                v2xab

    CHARACTER(len=*), PARAMETER              :: procedureN = 'x_u_hjs_lc'
    REAL(real_8), PARAMETER                  :: tol = 1.0e-20_real_8 

    REAL(real_8) :: Ax, dEGadnaa, dEGadra, dEGbdnbb, dEGbdrb, dEnerdnaa, &
      dEnerdnbb, dEnerdra, dEnerdrb, detaadnaa, detaadra, detabdnbb, &
      detabdrb, dFadnaa, dFadra, dFbdnbb, dFbdrb, dFxadnaa, dFxadra, &
      dFxbdnbb, dFxbdrb, dkfadra, dkfbdrb, dlambadaadnaa, dlambadaadra, &
      dlambadabdnbb, dlambadabdrb, dsadnaa, dsadra, dsbdnbb, dsbdrb, dy0adra, &
      dy0bdrb, dzetaadnaa, dzetaadra, dzetabdnbb, dzetabdrb, EGa, EGb, Ener, &
      etaa, etab, Fa, Fb, Fxa, Fxb, kfa, kfb, lambadaa, lambadab, sa, sb, &
      y0a, y0b, zetaa, zetab
    REAL(real_8), EXTERNAL                   :: dei

    IF (r0.LT.1.0e-3_real_8) CALL stopgm(procedureN,'RO too small',& 
         __LINE__,__FILE__)
    INCLUDE 'uhjsx.inc'
    sx = Ener
    v1xa = dEnerdra
    v1xb = dEnerdrb
    v2xa = 2.0_real_8*dEnerdnaa ! multiply by 2, cause we take derivative vs
    ! (nabla rho) and not (nabla rho)**2
    v2xb = 2.0_real_8*dEnerdnbb ! multiply by 2, cause we take derivative vs
    ! (nabla rho) and not (nabla rho)**2
    v2xab = 0.0_real_8
    ! ==--------------------------------------------------------------==
  END SUBROUTINE x_u_hjs_lc
  ! ==================================================================
END MODULE x_hjs
! ==================================================================
