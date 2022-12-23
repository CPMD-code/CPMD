MODULE ylmr2_utils
  USE cnst,                            ONLY: fpi
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: paral

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ylmr2

CONTAINS

  SUBROUTINE ylmr2(l,nhg,hg,gk,ylm)
    ! ==--------------------------------------------------------------==
    ! REAL SPHERICAL HARMONICS,  L IS COMBINED INDEX FOR LM (L=1,2...25)
    ! ORDER:  S, P_X, P_Z, P_Y, D_XY, D_XZ, D_Z^2, D_YZ, D_X^2-Y^2  ....
    ! THE REAL SPHERICAL HARMONICS USED HERE FORM BASES FOR THE
    ! IRREDUCIBLE REPRESENTATIONS OF THE GROUP O
    ! 
    ! SEE WIESSBLUTH 'ATOMS AND MOLECULES' PAGES 128-130
    ! ERRORS IN WEISSBLUTH HAVE BEEN CORRECTED:
    ! 1.) ELIMINATION OF THE 7 FROM L=20
    ! 2.) ADDITION OF THE FACTOR 1./SQRT(12.) TO L=25
    ! 
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: l, nhg
    REAL(real_8)                             :: hg(*), gk(3,*), ylm(*)

    REAL(real_8), PARAMETER                  :: thr = 1.0e-6_real_8

    INTEGER                                  :: ig
    REAL(real_8)                             :: c

! ==--------------------------------------------------------------==
! NOTE :   Y_LM (G=0) = SQRT(FPI)  WHEN L=0  AND  = 0  WHEN L>0
! ==--------------------------------------------------------------==

    IF (l.EQ.1) THEN
       c=SQRT(1._real_8/fpi)
       DO ig=1,nhg
          ylm(ig) = c
       ENDDO
    ELSE IF (l.EQ.2) THEN
       c=SQRT(3._real_8/fpi)
       !$omp parallel do private(IG) shared(YLM)
       DO ig=1,nhg
          IF (ABS(hg(ig)).GE.thr) THEN
             ylm(ig) = c*gk(1,ig)/SQRT(hg(ig))   !   X
          ELSE
             ylm(ig) = 0.0_real_8
             hg(ig)  = 0.0_real_8
          ENDIF
       ENDDO
    ELSE IF (l.EQ.3) THEN
       c=SQRT(3._real_8/fpi)
       !$omp parallel do private(IG) shared(YLM)
       DO ig=1,nhg
          IF (ABS(hg(ig)).GE.thr) THEN
             ylm(ig) = c*gk(3,ig)/SQRT(hg(ig))   !   Z
          ELSE
             ylm(ig) = 0.0_real_8
             hg(ig)  = 0.0_real_8
          ENDIF
       ENDDO
    ELSE IF (l.EQ.4) THEN
       c=SQRT(3._real_8/fpi)
       !$omp parallel do private(IG) shared(YLM)
       DO ig=1,nhg
          IF (ABS(hg(ig)).GE.thr) THEN
             ylm(ig) = c*gk(2,ig)/SQRT(hg(ig))! Y
          ELSE
             ylm(ig) = 0.0_real_8
             hg(ig)  = 0.0_real_8
          ENDIF
       ENDDO
    ELSE IF (l.EQ.5) THEN
       c=SQRT(15._real_8/fpi)
       !$omp parallel do private(IG) shared(YLM)
       DO ig=1,nhg
          IF (ABS(hg(ig)).GE.thr) THEN
             ylm(ig) = c*gk(1,ig)*gk(2,ig)/hg(ig)! X*Y
          ELSE
             ylm(ig) = 0.0_real_8
             hg(ig)  = 0.0_real_8
          ENDIF
       ENDDO
    ELSE IF (l.EQ.6) THEN
       c=SQRT(15._real_8/fpi)
       !$omp parallel do private(IG) shared(YLM)
       DO ig=1,nhg
          IF (ABS(hg(ig)).GE.thr) THEN
             ylm(ig) = c*gk(1,ig)*gk(3,ig)/hg(ig)! X*Z
          ELSE
             ylm(ig) = 0.0_real_8
             hg(ig)  = 0.0_real_8
          ENDIF
       ENDDO
    ELSE IF (l.EQ.7) THEN
       c=SQRT(5._real_8/fpi/4._real_8)
       !$omp parallel do private(IG) shared(YLM)
       DO ig=1,nhg
          IF (ABS(hg(ig)).GE.thr) THEN
             ylm(ig) = c*(3._real_8*gk(3,ig)**2/hg(ig)-1._real_8)! (3.*Z*Z-1.0)
          ELSE
             ylm(ig) = 0.0_real_8
             hg(ig)  = 0.0_real_8
          ENDIF
       ENDDO
    ELSE IF (l.EQ.8) THEN
       c=SQRT(15._real_8/fpi)
       !$omp parallel do private(IG) shared(YLM)
       DO ig=1,nhg
          IF (ABS(hg(ig)).GE.thr) THEN
             ylm(ig) = c*gk(2,ig)*gk(3,ig)/hg(ig)! Y*Z
          ELSE
             ylm(ig) = 0.0_real_8
             hg(ig)  = 0.0_real_8
          ENDIF
       ENDDO
    ELSE IF (l.EQ.9) THEN
       c=SQRT(15._real_8/fpi/4._real_8)
       !$omp parallel do private(IG) shared(YLM)
       DO ig=1,nhg
          IF (ABS(hg(ig)).GE.thr) THEN
             ylm(ig) = c*(gk(1,ig)**2-gk(2,ig)**2)/hg(ig)! X*X-Y*Y
          ELSE
             ylm(ig) = 0.0_real_8
             hg(ig)  = 0.0_real_8
          ENDIF
       ENDDO
    ELSE IF (l.EQ.10) THEN
       c=SQRT(7._real_8/fpi)*5._real_8/2._real_8
       !$omp parallel do private(IG) shared(YLM)
       DO ig=1,nhg
          IF (ABS(hg(ig)).GE.thr) THEN
             ylm(ig) = c*gk(1,ig)*(gk(1,ig)**2-0.6_real_8*hg(ig))/&
                  (HG(IG)*SQRT(HG(IG)))                ! X(X^2-3R^2/5)
          ELSE
             ylm(ig) = 0.0_real_8
             hg(ig)  = 0.0_real_8
          ENDIF
       ENDDO
    ELSE IF (l.EQ.11) THEN
       c=SQRT(7._real_8/fpi)*5._real_8/2._real_8
       !$omp parallel do private(IG) shared(YLM)
       DO ig=1,nhg                              ! Y(Y^2-3R^2/5)
          IF (ABS(hg(ig)).GE.thr) THEN
             ylm(ig) = c*gk(2,ig)*(gk(2,ig)**2-0.6_real_8*hg(ig))/&
                  (HG(IG)*SQRT(HG(IG)))
          ELSE
             ylm(ig) = 0.0_real_8
             hg(ig)  = 0.0_real_8
          ENDIF
       ENDDO
    ELSE IF (l.EQ.12) THEN
       c=SQRT(7._real_8*15._real_8/fpi)
       !$omp parallel do private(IG) shared(YLM)
       DO ig=1,nhg                              ! XYZ
          IF (ABS(hg(ig)).GE.thr) THEN
             ylm(ig) = c*gk(1,ig)*gk(2,ig)*gk(3,ig)/(hg(ig)*SQRT(hg(ig)))
          ELSE
             ylm(ig) = 0.0_real_8
             hg(ig)  = 0.0_real_8
          ENDIF
       ENDDO
    ELSE IF (l.EQ.13) THEN
       c=SQRT(7._real_8/fpi)*5._real_8/2._real_8
       !$omp parallel do private(IG) shared(YLM)
       DO ig=1,nhg                              ! Z(Z^2-.6R^2)
          IF (ABS(hg(ig)).GE.thr) THEN
             ylm(ig) = c*gk(3,ig)*(gk(3,ig)**2-0.6_real_8*hg(ig))/&
                  (HG(IG)*SQRT(HG(IG)))
          ELSE
             ylm(ig) = 0.0_real_8
             hg(ig)  = 0.0_real_8
          ENDIF
       ENDDO
    ELSE IF (l.EQ.14) THEN
       c=SQRT(7._real_8*15._real_8/fpi)/2._real_8
       !$omp parallel do private(IG) shared(YLM)
       DO ig=1,nhg                              ! Z(X^2-Y^2)
          IF (ABS(hg(ig)).GE.thr) THEN
             ylm(ig) = c*gk(3,ig)*(gk(1,ig)**2-gk(2,ig)**2)/&
                  (HG(IG)*SQRT(HG(IG)))
          ELSE
             ylm(ig) = 0.0_real_8
             hg(ig)  = 0.0_real_8
          ENDIF
       ENDDO
    ELSE IF (l.EQ.15) THEN
       c=SQRT(7._real_8*15._real_8/fpi)/2._real_8
       !$omp parallel do private(IG) shared(YLM)
       DO ig=1,nhg                              ! Y(Z^2-X^2)
          IF (ABS(hg(ig)).GE.thr) THEN
             ylm(ig) = c*gk(2,ig)*(gk(3,ig)**2-gk(1,ig)**2)/&
                  (HG(IG)*SQRT(HG(IG)))
          ELSE
             ylm(ig) = 0.0_real_8
             hg(ig)  = 0.0_real_8
          ENDIF
       ENDDO
    ELSE IF (l.EQ.16) THEN
       c=SQRT(7._real_8*15._real_8/fpi)/2._real_8
       !$omp parallel do private(IG) shared(YLM)
       DO ig=1,nhg                              ! X(Y^2-Z^2)
          IF (ABS(hg(ig)).GE.thr) THEN
             ylm(ig) = c*gk(1,ig)*(gk(2,ig)**2-gk(3,ig)**2)/&
                  (HG(IG)*SQRT(HG(IG)))
          ELSE
             ylm(ig) = 0.0_real_8
             hg(ig)  = 0.0_real_8
          ENDIF
       ENDDO
    ELSE IF (l.EQ.17) THEN
       c=SQRT(3._real_8*7._real_8/fpi)*5._real_8/4._real_8
       !$omp parallel do private(IG) shared(YLM)
       DO ig=1,nhg                              ! A1
          IF (ABS(hg(ig)).GE.thr) THEN
             ylm(ig) = c*((gk(1,ig)**4+gk(2,ig)**4+gk(3,ig)**4)/&
                  (HG(IG)*HG(IG))-0.6_real_8)
          ELSE
             ylm(ig) = 0.0_real_8
             hg(ig)  = 0.0_real_8
          ENDIF
       ENDDO
    ELSE IF (l.EQ.18) THEN
       c=SQRT(9._real_8*35._real_8/fpi)/2._real_8
       !$omp parallel do private(IG) shared(YLM)
       DO ig=1,nhg                              ! YZ(Y^2-Z^2)
          IF (ABS(hg(ig)).GE.thr) THEN
             ylm(ig) = c*gk(2,ig)*gk(3,ig)*(gk(2,ig)**2-gk(3,ig)**2)/&
                  (HG(IG)*HG(IG))
          ELSE
             ylm(ig) = 0.0_real_8
             hg(ig)  = 0.0_real_8
          ENDIF
       ENDDO
    ELSE IF (l.EQ.19) THEN
       c=SQRT(9._real_8*35._real_8/fpi)/2._real_8
       !$omp parallel do private(IG) shared(YLM)
       DO ig=1,nhg                              ! ZX(Z^2-X^2)
          IF (ABS(hg(ig)).GE.thr) THEN
             ylm(ig) = c*gk(1,ig)*gk(3,ig)*(gk(3,ig)**2-gk(1,ig)**2)/&
                  (HG(IG)*HG(IG))
          ELSE
             ylm(ig) = 0.0_real_8
             hg(ig)  = 0.0_real_8
          ENDIF
       ENDDO
    ELSE IF (l.EQ.20) THEN
       c=SQRT(9._real_8*5._real_8/fpi)/4._real_8
       !$omp parallel do private(IG) shared(YLM)
       DO ig=1,nhg                              ! E\EPSILON
          IF (ABS(hg(ig)).GE.thr) THEN
             ylm(ig) = c*((gk(1,ig)**4-gk(2,ig)**4)-&
                  6._real_8*GK(3,IG)**2*(GK(1,IG)**2-GK(2,IG)**2))/(HG(IG)*HG(IG))
          ELSE
             ylm(ig) = 0.0_real_8
             hg(ig)  = 0.0_real_8
          ENDIF
       ENDDO
    ELSE IF (l.EQ.21) THEN
       c=SQRT(9._real_8*35._real_8/fpi)/2._real_8
       !$omp parallel do private(IG) shared(YLM)
       DO ig=1,nhg                              ! XY(X^2-Y^2)
          IF (ABS(hg(ig)).GE.thr) THEN
             ylm(ig) = c*gk(1,ig)*gk(2,ig)*(gk(1,ig)**2-gk(2,ig)**2)/&
                  (HG(IG)*HG(IG))
          ELSE
             ylm(ig) = 0.0_real_8
             hg(ig)  = 0.0_real_8
          ENDIF
       ENDDO
    ELSE IF (l.EQ.22) THEN
       c=SQRT(9._real_8*5._real_8/fpi)*7._real_8/2._real_8
       !$omp parallel do private(IG) shared(YLM)
       DO ig=1,nhg                              ! XY(Z^2-1/7*R^2)
          IF (ABS(hg(ig)).GE.thr) THEN
             ylm(ig) = c*gk(1,ig)*gk(2,ig)*(gk(3,ig)**2-hg(ig)/7._real_8)/&
                  (HG(IG)*HG(IG))
          ELSE
             ylm(ig) = 0.0_real_8
             hg(ig)  = 0.0_real_8
          ENDIF
       ENDDO
    ELSE IF (l.EQ.23) THEN
       c=SQRT(9._real_8*5._real_8/fpi)*7._real_8/2._real_8
       !$omp parallel do private(IG) shared(YLM)
       DO ig=1,nhg                              ! ZX(Y^2-1/7*R^2)
          IF (ABS(hg(ig)).GE.thr) THEN
             ylm(ig) = c*gk(1,ig)*gk(3,ig)*(gk(2,ig)**2-hg(ig)/7._real_8)/&
                  (HG(IG)*HG(IG))
          ELSE
             ylm(ig) = 0.0_real_8
             hg(ig)  = 0.0_real_8
          ENDIF
       ENDDO
    ELSE IF (l.EQ.24) THEN
       c=SQRT(9._real_8*5._real_8/fpi)*7._real_8/2._real_8
       !$omp parallel do private(IG) shared(YLM)
       DO ig=1,nhg                              ! YZ(X^2-1/7*R^2)
          IF (ABS(hg(ig)).GE.thr) THEN
             ylm(ig) = c*gk(2,ig)*gk(3,ig)*(gk(1,ig)**2-hg(ig)/7._real_8)/&
                  (HG(IG)*HG(IG))
          ELSE
             ylm(ig) = 0.0_real_8
             hg(ig)  = 0.0_real_8
          ENDIF
       ENDDO
    ELSE IF (l.EQ.25) THEN
       c=SQRT(9._real_8*5._real_8/fpi/3._real_8)*7._real_8/2._real_8
       !$omp parallel do private(IG) shared(YLM)
       DO ig=1,nhg                              ! E\THETA
          IF (ABS(hg(ig)).GE.thr) THEN
             ylm(ig) = c*( gk(3,ig)**4-0.5_real_8*(gk(1,ig)**4+gk(2,ig)**4)-&
                  6._real_8/7._real_8*HG(IG)&
                  *(gk(3,ig)**2-0.5_real_8*(gk(1,ig)**2+gk(2,ig)**2) ))&
                  /( HG(IG)*HG(IG) )
          ELSE
             ylm(ig) = 0.0_real_8
             hg(ig)  = 0.0_real_8
          ENDIF
       ENDDO
    ELSE IF (l.GE.26) THEN
       IF (paral%io_parent)&
            WRITE(6,*) 'GIVEN L=',l,' [MAX L = 26]'
       CALL stopgm(' YLMR2',' HIGHER L NOT PROGRAMMED ',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE ylmr2
  ! ==================================================================

END MODULE ylmr2_utils
