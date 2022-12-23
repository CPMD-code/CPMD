MODULE dylmr_utils
  USE cnst,                            ONLY: fpi
  USE error_handling,                  ONLY: stopgm
  USE geq0mod,                         ONLY: geq0
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: paral
  USE str2,                            ONLY: gagk
  USE strs,                            ONLY: alpha,&
                                             beta,&
                                             delta
  USE system,                          ONLY: ncpw,&
                                             parm
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: dylmr
  PUBLIC :: dbess
  PUBLIC :: dylmrx
  PUBLIC :: dylmrkx

CONTAINS

#ifdef __SR8000
  !option MP(P(0)), LANGLVL(SAVE(0))
#endif
  ! ==================================================================
  SUBROUTINE dylmr(l,gk,dylm,kk)
    ! ==--------------------------------------------------------------==
    ! REAL SPHERICAL HARMONICS,  L IS COMBINED INDEX FOR LM  (L=1,2...9)
    ! ORDER:  S, P_X, P_Z, P_Y, D_XY, D_XZ, D_Z^2, D_YZ, D_X^2-Y^2
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: l
    REAL(real_8)                             :: gk(3,ncpw%nhg), dylm(ncpw%nhg)
    INTEGER                                  :: kk

    REAL(real_8), PARAMETER                  :: divcut = 1.0e-6_real_8 

    INTEGER                                  :: ig, is1
    LOGICAL, SAVE                            :: twarnc = .TRUE.
    REAL(real_8)                             :: hgg, r, r2, x(3), c

    is1=1
    IF (geq0.AND.l.GT.1) THEN
       hgg = SQRT(gk(1,1)*gk(1,1)+gk(2,1)*gk(2,1)+gk(3,1)*gk(3,1))
       IF (hgg.LT.divcut) THEN
          is1=2
          dylm(1)=0.0_real_8
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (l.EQ.1) THEN
       DO ig=1,ncpw%nhg
          dylm(ig) = 0._real_8
       ENDDO
    ELSE IF (l.EQ.2) THEN
       !$omp parallel do private(IG,X,R,R2) shared(DYLM)
       DO ig=is1,ncpw%nhg
          x(1) = gk(1,ig)
          x(2) = gk(2,ig)
          x(3) = gk(3,ig)
          r = MAX(SQRT(x(1)*x(1)+x(2)*x(2)+x(3)*x(3)),divcut)
          r2 = r*r*parm%tpiba2
          x(1) = x(1)/r
          x(2) = x(2)/r
          x(3) = x(3)/r
          dylm(ig) = SQRT(3._real_8/fpi)*( x(1)*gagk(ig,kk)/r2 -&
               x(alpha(kk))*delta(1,beta(kk)) )
       ENDDO
    ELSE IF (l.EQ.3) THEN
       !$omp parallel do private(IG,X,R,R2) shared(DYLM)
       DO ig=is1,ncpw%nhg
          x(1) = gk(1,ig)
          x(2) = gk(2,ig)
          x(3) = gk(3,ig)
          r = MAX(SQRT(x(1)*x(1)+x(2)*x(2)+x(3)*x(3)),divcut)
          r2 = r*r*parm%tpiba2
          x(1) = x(1)/r
          x(2) = x(2)/r
          x(3) = x(3)/r
          dylm(ig) = SQRT(3._real_8/fpi)*( x(3)*gagk(ig,kk)/r2 -&
               x(alpha(kk))*delta(3,beta(kk)) )
       ENDDO
    ELSE IF (l.EQ.4) THEN
       !$omp parallel do private(IG,X,R,R2) shared(DYLM) 
       DO ig=is1,ncpw%nhg
          x(1) = gk(1,ig)
          x(2) = gk(2,ig)
          x(3) = gk(3,ig)
          r = MAX(SQRT(x(1)*x(1)+x(2)*x(2)+x(3)*x(3)),divcut)
          r2 = r*r*parm%tpiba2
          x(1) = x(1)/r
          x(2) = x(2)/r
          x(3) = x(3)/r
          dylm(ig) = SQRT(3._real_8/fpi)*( x(2)*gagk(ig,kk)/r2 -&
               x(alpha(kk))*delta(2,beta(kk)) )
       ENDDO
    ELSE IF (l.EQ.5) THEN
       !$omp parallel do private(IG,X,R,R2) shared(DYLM)
       DO ig=is1,ncpw%nhg
          x(1) = gk(1,ig)
          x(2) = gk(2,ig)
          x(3) = gk(3,ig)
          r = MAX(SQRT(x(1)*x(1)+x(2)*x(2)+x(3)*x(3)),divcut)
          r2 = r*r*parm%tpiba2
          x(1) = x(1)/r
          x(2) = x(2)/r
          x(3) = x(3)/r
          dylm(ig) = SQRT(15._real_8/fpi)*( 2._real_8*x(1)*x(2)*gagk(ig,kk)/r2 -&
               x(alpha(kk))*x(2)*delta(1,beta(kk)) -&
               x(alpha(kk))*x(1)*delta(2,beta(kk)) )
       ENDDO
    ELSE IF (l.EQ.6) THEN
       !$omp parallel do private(IG,X,R,R2) shared(DYLM)
       DO ig=is1,ncpw%nhg
          x(1) = gk(1,ig)
          x(2) = gk(2,ig)
          x(3) = gk(3,ig)
          r = MAX(SQRT(x(1)*x(1)+x(2)*x(2)+x(3)*x(3)),divcut)
          r2 = r*r*parm%tpiba2
          x(1) = x(1)/r
          x(2) = x(2)/r
          x(3) = x(3)/r
          dylm(ig) = SQRT(15._real_8/fpi)*( 2._real_8*x(1)*x(3)*gagk(ig,kk)/r2 -&
               x(alpha(kk))*x(3)*delta(1,beta(kk)) -&
               x(alpha(kk))*x(1)*delta(3,beta(kk)) )
       ENDDO
    ELSE IF (l.EQ.7) THEN
       !$omp parallel do private(IG,X,R,R2) shared(DYLM)
       DO ig=is1,ncpw%nhg
          x(1) = gk(1,ig)
          x(2) = gk(2,ig)
          x(3) = gk(3,ig)
          r = MAX(SQRT(x(1)*x(1)+x(2)*x(2)+x(3)*x(3)),divcut)
          r2 = r*r*parm%tpiba2
          x(1) = x(1)/r
          x(2) = x(2)/r
          x(3) = x(3)/r
          dylm(ig) = 6._real_8*SQRT(5._real_8/fpi/4._real_8)*( x(3)*x(3)*&
               gagk(ig,kk)/r2 -&
               x(alpha(kk))*x(3)*delta(3,beta(kk)) )
       ENDDO
    ELSE IF (l.EQ.8) THEN
       !$omp parallel do private(IG,X,R,R2) shared(DYLM)
       DO ig=is1,ncpw%nhg
          x(1) = gk(1,ig)
          x(2) = gk(2,ig)
          x(3) = gk(3,ig)
          r = MAX(SQRT(x(1)*x(1)+x(2)*x(2)+x(3)*x(3)),divcut)
          r2 = r*r*parm%tpiba2
          x(1) = x(1)/r
          x(2) = x(2)/r
          x(3) = x(3)/r
          dylm(ig) = SQRT(15._real_8/fpi)*(2._real_8*x(2)*x(3)*gagk(ig,kk)/r2 -&
               x(alpha(kk))*x(3)*delta(2,beta(kk)) -&
               x(alpha(kk))*x(2)*delta(3,beta(kk)) )
       ENDDO
    ELSE IF (l.EQ.9) THEN
       !$omp parallel do private(IG,X,R,R2) shared(DYLM)
       DO ig=is1,ncpw%nhg
          x(1) = gk(1,ig)
          x(2) = gk(2,ig)
          x(3) = gk(3,ig)
          r = MAX(SQRT(x(1)*x(1)+x(2)*x(2)+x(3)*x(3)),divcut)
          r2 = r*r*parm%tpiba2
          x(1) = x(1)/r
          x(2) = x(2)/r
          x(3) = x(3)/r
          dylm(ig) = 2._real_8*SQRT(15._real_8/fpi/4._real_8)*( (x(1)*x(1)-x(2)*x(2))*&
               gagk(ig,kk)/r2 - x(alpha(kk))*x(1)*delta(1,beta(kk)) +&
               x(alpha(kk))*x(2)*delta(2,beta(kk)) )
       ENDDO
       ! 
       ! FIXME:
       ! higher L values are not (yet) programmed. we try to work around it
       ! by setting DYLM in those cases to zero. since this subroutine is 
       ! only used for calculating the vanderbilt augmentation charge 
       ! contribution to the stress tensor, it might be good enough to do
       ! an approximate cell relaxation (the bulk of the contributions
       ! will be from the lower L cases anyway).
       ! 
       ! <axel.kohlmeyer@theochem.ruhr-uni-bochum.de>, 12/2003
       ! 
!   ELSE
!      IF (paral%parent.AND.twarnc) THEN
!         twarnc = .FALSE.
!         IF (paral%io_parent) THEN
!            WRITE(6,'(A,I2,A)') ' DYLMR| BIG FAT WARNING: L=',l,&
!                 ' IS NOT SUPPORTED (MAX 9)'
!            WRITE(6,'(A,A)') ' DYLMR| TRYING TO CONTINUE BY',&
!                 ' SETTING DYLM() TO 0.0 FOR L>9'
!            WRITE(6,'(A)')' DYLMR| STRESS TENSOR WILL NOT BE ACCURATE'
!         ENDIF
!      ENDIF
!      DO ig=is1,ncpw%nhg
!         dylm(ig) = 0.0_real_8
!      ENDDO
!   ENDIF
    ELSE IF (l.EQ.10) THEN
       c=SQRT(7._real_8/fpi)*5._real_8/2._real_8
       !$OMP parallel do private(IG,X,R,R2) shared(DYLM)
       DO ig=is1,ncpw%nhg
          x(1) = gk(1,ig)
          x(2) = gk(2,ig)
          x(3) = gk(3,ig)
          r = MAX(SQRT(x(1)*x(1)+x(2)*x(2)+x(3)*x(3)),divcut)
          r2 = r*r*parm%tpiba2
          x(1) = x(1)/r
          x(2) = x(2)/r
          x(3) = x(3)/r
          dylm(ig) = c*( &
               3.0_real_8*x(1)*(x(1)*x(1)-0.6_real_8)*gagk(ig,kk)/r2 &
              -0.6_real_8*x(alpha(kk))*( &
              (3.0_real_8*x(1)*x(1)-1.0_real_8)*delta(1,beta(kk)) &
              -2.0_real_8*x(1)*x(2)*delta(2,beta(kk)) &
              -2.0_real_8*x(1)*x(3)*delta(3,beta(kk))))
       ENDDO
    ELSE IF (l.EQ.11) THEN
       c=SQRT(7._real_8/fpi)*5._real_8/2._real_8
       !$OMP parallel do private(IG,X,R,R2) shared(DYLM)
       DO ig=is1,ncpw%nhg
          x(1) = gk(1,ig)
          x(2) = gk(2,ig)
          x(3) = gk(3,ig)
          r = MAX(SQRT(x(1)*x(1)+x(2)*x(2)+x(3)*x(3)),divcut)
          r2 = r*r*parm%tpiba2
          x(1) = x(1)/r
          x(2) = x(2)/r
          x(3) = x(3)/r
          dylm(ig) = c*( &
               3.0_real_8*x(2)*(x(2)*x(2)-0.6_real_8)*gagk(ig,kk)/r2 &
              -0.6_real_8*x(alpha(kk))*( &
              -2.0_real_8*x(2)*x(1)*delta(1,beta(kk)) &
             +(3.0_real_8*x(2)*x(2)-1.0_real_8)*delta(2,beta(kk)) &
              -2.0_real_8*x(2)*x(3)*delta(3,beta(kk))))
       ENDDO
    ELSE IF (l.EQ.12) THEN
       c=SQRT(7._real_8*15._real_8/fpi)
       !$OMP parallel do private(IG,X,R,R2) shared(DYLM)
       DO ig=is1,ncpw%nhg
          x(1) = gk(1,ig)
          x(2) = gk(2,ig)
          x(3) = gk(3,ig)
          r = MAX(SQRT(x(1)*x(1)+x(2)*x(2)+x(3)*x(3)),divcut)
          r2 = r*r*parm%tpiba2
          x(1) = x(1)/r
          x(2) = x(2)/r
          x(3) = x(3)/r
          dylm(ig) = c*( &
               3.0_real_8*x(1)*x(2)*x(3)*gagk(ig,kk)/r2 &
              -x(alpha(kk))*( &
               x(2)*x(3)*delta(1,beta(kk)) &
              +x(1)*x(3)*delta(2,beta(kk)) &
              +x(1)*x(2)*delta(3,beta(kk))))
       ENDDO
    ELSE IF (l.EQ.13) THEN
       c=SQRT(7._real_8/fpi)*5._real_8/2._real_8
       !$OMP parallel do private(IG,X,R,R2) shared(DYLM)
       DO ig=is1,ncpw%nhg
          x(1) = gk(1,ig)
          x(2) = gk(2,ig)
          x(3) = gk(3,ig)
          r = MAX(SQRT(x(1)*x(1)+x(2)*x(2)+x(3)*x(3)),divcut)
          r2 = r*r*parm%tpiba2
          x(1) = x(1)/r
          x(2) = x(2)/r
          x(3) = x(3)/r
          dylm(ig) = c*( &
               3.0_real_8*x(3)*(x(3)*x(3)-0.6_real_8)*gagk(ig,kk)/r2 &
              -0.6_real_8*x(alpha(kk))*( &
              -2.0_real_8*x(3)*x(1)*delta(1,beta(kk)) &
              -2.0_real_8*x(3)*x(2)*delta(2,beta(kk)) &
             +(3.0_real_8*x(3)*x(3)-1.0_real_8)*delta(3,beta(kk))))
       ENDDO
    ELSE IF (l.EQ.14) THEN
       c=SQRT(7._real_8*15._real_8/fpi)/2._real_8
       !$OMP parallel do private(IG,X,R,R2) shared(DYLM)
       DO ig=is1,ncpw%nhg
          x(1) = gk(1,ig)
          x(2) = gk(2,ig)
          x(3) = gk(3,ig)
          r = MAX(SQRT(x(1)*x(1)+x(2)*x(2)+x(3)*x(3)),divcut)
          r2 = r*r*parm%tpiba2
          x(1) = x(1)/r
          x(2) = x(2)/r
          x(3) = x(3)/r
          dylm(ig) = c*( &
               3.0_real_8*x(3)*(x(1)*x(1)-x(2)*x(2))*gagk(ig,kk)/r2 &
              -x(alpha(kk))*( &
               2.0_real_8*x(3)*x(1)*delta(1,beta(kk)) &
              -2.0_real_8*x(3)*x(2)*delta(2,beta(kk)) &
              +(x(1)*x(1)-x(2)*x(2))*delta(3,beta(kk))))
       ENDDO
    ELSE IF (l.EQ.15) THEN
       c=SQRT(7._real_8*15._real_8/fpi)/2._real_8
       !$OMP parallel do private(IG,X,R,R2) shared(DYLM)
       DO ig=is1,ncpw%nhg
          x(1) = gk(1,ig)
          x(2) = gk(2,ig)
          x(3) = gk(3,ig)
          r = MAX(SQRT(x(1)*x(1)+x(2)*x(2)+x(3)*x(3)),divcut)
          r2 = r*r*parm%tpiba2
          x(1) = x(1)/r
          x(2) = x(2)/r
          x(3) = x(3)/r
          dylm(ig) = c*( &
               3.0_real_8*x(2)*(x(3)*x(3)-x(1)*x(1))*gagk(ig,kk)/r2 &
              -x(alpha(kk))*( &
              -2.0_real_8*x(2)*x(1)*delta(1,beta(kk)) &
              +(x(3)*x(3)-x(1)*x(1))*delta(2,beta(kk)) &
              +2.0_real_8*x(2)*x(3)*delta(3,beta(kk))))
       ENDDO
    ELSE IF (l.EQ.16) THEN
       c=SQRT(7._real_8*15._real_8/fpi)/2._real_8
       !$OMP parallel do private(IG,X,R,R2) shared(DYLM)
       DO ig=is1,ncpw%nhg
          x(1) = gk(1,ig)
          x(2) = gk(2,ig)
          x(3) = gk(3,ig)
          r = MAX(SQRT(x(1)*x(1)+x(2)*x(2)+x(3)*x(3)),divcut)
          r2 = r*r*parm%tpiba2
          x(1) = x(1)/r
          x(2) = x(2)/r
          x(3) = x(3)/r
          dylm(ig) = c*( &
               3.0_real_8*x(1)*(x(2)*x(2)-x(3)*x(3))*gagk(ig,kk)/r2 &
              -x(alpha(kk))*( &
               (x(2)*x(2)-x(3)*x(3))*delta(1,beta(kk)) &
              +2.0_real_8*x(1)*x(2)*delta(2,beta(kk)) &
              -2.0_real_8*x(1)*x(3)*delta(3,beta(kk))))
       ENDDO
    ELSE IF (l.EQ.17) THEN
       c=SQRT(3._real_8*7._real_8/fpi)*5._real_8/4._real_8
       !$OMP parallel do private(IG,X,R,R2) shared(DYLM)
       DO ig=is1,ncpw%nhg
          x(1) = gk(1,ig)
          x(2) = gk(2,ig)
          x(3) = gk(3,ig)
          r = MAX(SQRT(x(1)*x(1)+x(2)*x(2)+x(3)*x(3)),divcut)
          r2 = r*r*parm%tpiba2
          x(1) = x(1)/r
          x(2) = x(2)/r
          x(3) = x(3)/r
          dylm(ig) = c*( &
               4.0_real_8*(x(1)**4+x(2)**4+x(3)**4)*gagk(ig,kk)/r2 &
              -4.0_real_8*x(alpha(kk))*( &
               x(1)**3*delta(1,beta(kk)) &
              +x(2)**3*delta(2,beta(kk)) &
              +x(3)**3*delta(3,beta(kk))))
       ENDDO
    ELSE IF (l.EQ.18) THEN
       c=SQRT(9._real_8*35._real_8/fpi)/2._real_8
       !$OMP parallel do private(IG,X,R,R2) shared(DYLM)
       DO ig=is1,ncpw%nhg
          x(1) = gk(1,ig)
          x(2) = gk(2,ig)
          x(3) = gk(3,ig)
          r = MAX(SQRT(x(1)*x(1)+x(2)*x(2)+x(3)*x(3)),divcut)
          r2 = r*r*parm%tpiba2
          x(1) = x(1)/r
          x(2) = x(2)/r
          x(3) = x(3)/r
          dylm(ig) = c*( &
               4.0_real_8*x(2)*x(3)*(x(2)*x(2)-x(3)*x(3))*gagk(ig,kk)/r2 &
              -x(alpha(kk))*( &
               (3.0_real_8*x(2)*x(2)-x(3)*x(3))*x(3)*delta(2,beta(kk)) &
              -(3.0_real_8*x(3)*x(3)-x(2)*x(2))*x(2)*delta(3,beta(kk))))
       ENDDO
    ELSE IF (l.EQ.19) THEN
       c=SQRT(9._real_8*35._real_8/fpi)/2._real_8
       !$OMP parallel do private(IG,X,R,R2) shared(DYLM)
       DO ig=is1,ncpw%nhg
          x(1) = gk(1,ig)
          x(2) = gk(2,ig)
          x(3) = gk(3,ig)
          r = MAX(SQRT(x(1)*x(1)+x(2)*x(2)+x(3)*x(3)),divcut)
          r2 = r*r*parm%tpiba2
          x(1) = x(1)/r
          x(2) = x(2)/r
          x(3) = x(3)/r
          dylm(ig) = c*( &
               4.0_real_8*x(1)*x(3)*(x(3)*x(3)-x(1)*x(1))*gagk(ig,kk)/r2 &
              -x(alpha(kk))*( &
              -(3.0_real_8*x(1)*x(1)-x(3)*x(3))*x(3)*delta(1,beta(kk)) &
              +(3.0_real_8*x(3)*x(3)-x(1)*x(1))*x(1)*delta(3,beta(kk))))
       ENDDO
    ELSE IF(l.EQ.20) THEN
       c=SQRT(9._real_8*5._real_8/fpi)/4._real_8
       !$OMP parallel do private(IG,X,R,R2) shared(DYLM)
       DO ig=is1,ncpw%nhg
          x(1) = gk(1,ig)
          x(2) = gk(2,ig)
          x(3) = gk(3,ig)
          r = MAX(SQRT(x(1)*x(1)+x(2)*x(2)+x(3)*x(3)),divcut)
          r2 = r*r*parm%tpiba2
          x(1) = x(1)/r
          x(2) = x(2)/r
          x(3) = x(3)/r
          dylm(ig) = c*( &
               4.0_real_8*((x(1)**4-x(2)**4)-6.0_real_8*x(3)**2*(x(1)**2-x(2)**2)) &
               *gagk(ig,kk)/r2 -x(alpha(kk))*( &
               (4.0_real_8*x(1)**3-12.0_real_8*x(3)**2*x(1))*delta(1,beta(kk)) &
              -(4.0_real_8*x(2)**3-12.0_real_8*x(3)**2*x(2))*delta(2,beta(kk)) &
              -12.0_real_8*x(3)*(x(1)**2-x(2)**2)*delta(3,beta(kk))))
       ENDDO
    ELSE IF (l.EQ.21) THEN
       c=SQRT(9._real_8*35._real_8/fpi)/2._real_8
       !$OMP parallel do private(IG,X,R,R2) shared(DYLM)
       DO ig=is1,ncpw%nhg
          x(1) = gk(1,ig)
          x(2) = gk(2,ig)
          x(3) = gk(3,ig)
          r = MAX(SQRT(x(1)*x(1)+x(2)*x(2)+x(3)*x(3)),divcut)
          r2 = r*r*parm%tpiba2
          x(1) = x(1)/r
          x(2) = x(2)/r
          x(3) = x(3)/r
          dylm(ig) = c*( &
               4.0_real_8*x(1)*x(2)*(x(1)*x(1)-x(2)*x(2))*gagk(ig,kk)/r2 &
              -x(alpha(kk))*( &
              +(3.0_real_8*x(1)*x(1)-x(2)*x(2))*x(2)*delta(1,beta(kk)) &
              -(3.0_real_8*x(2)*x(2)-x(1)*x(1))*x(1)*delta(2,beta(kk))))
       ENDDO
   ELSE IF (l.EQ.22) THEN
       c=SQRT(9._real_8*5._real_8/fpi)*7._real_8/2._real_8
       !$OMP parallel do private(IG,X,R,R2) shared(DYLM)
       DO ig=is1,ncpw%nhg
          x(1) = gk(1,ig)
          x(2) = gk(2,ig)
          x(3) = gk(3,ig)
          r = MAX(SQRT(x(1)*x(1)+x(2)*x(2)+x(3)*x(3)),divcut)
          r2 = r*r*parm%tpiba2
          x(1) = x(1)/r
          x(2) = x(2)/r
          x(3) = x(3)/r
          dylm(ig) = c/7.0_real_8*( &
               4.0_real_8*x(1)*x(2)*(7.0_real_8*x(3)**2-1.0_real_8)*gagk(ig,kk)/r2 &
              -x(alpha(kk))*( &
               x(2)*(7.0_real_8*x(3)**2-2.0_real_8*x(1)**2-1.0_real_8)*delta(1,beta(kk)) &
              +x(1)*(7.0_real_8*x(3)**2-2.0_real_8*x(2)**2-1.0_real_8)*delta(2,beta(kk)) &
              +12.0_real_8*x(1)*x(2)*x(3)*delta(3,beta(kk))))
       ENDDO
    ELSE IF (l.EQ.23) THEN
       c=SQRT(9._real_8*5._real_8/fpi)*7._real_8/2._real_8
       !$OMP parallel do private(IG,X,R,R2) shared(DYLM)
       DO ig=is1,ncpw%nhg
          x(1) = gk(1,ig)
          x(2) = gk(2,ig)
          x(3) = gk(3,ig)
          r = MAX(SQRT(x(1)*x(1)+x(2)*x(2)+x(3)*x(3)),divcut)
          r2 = r*r*parm%tpiba2
          x(1) = x(1)/r
          x(2) = x(2)/r
          x(3) = x(3)/r
          dylm(ig) = c/7.0_real_8*( &
               4.0_real_8*x(3)*x(1)*(7.0_real_8*x(2)**2-1.0_real_8)*gagk(ig,kk)/r2 &
              -x(alpha(kk))*( &
               x(3)*(7.0_real_8*x(2)**2-2.0_real_8*x(1)**2-1.0_real_8)*delta(1,beta(kk)) &
              +12.0_real_8*x(1)*x(2)*x(3)*delta(2,beta(kk)) &
              +x(1)*(7.0_real_8*x(2)**2-2.0_real_8*x(3)**2-1.0_real_8)*delta(3,beta(kk))))
       ENDDO
    ELSE IF (l.EQ.24) THEN
       c=SQRT(9._real_8*5._real_8/fpi)*7._real_8/2._real_8
       !$OMP parallel do private(IG,X,R,R2) shared(DYLM)
       DO ig=is1,ncpw%nhg
          x(1) = gk(1,ig)
          x(2) = gk(2,ig)
          x(3) = gk(3,ig)
          r = MAX(SQRT(x(1)*x(1)+x(2)*x(2)+x(3)*x(3)),divcut)
          r2 = r*r*parm%tpiba2
          x(1) = x(1)/r
          x(2) = x(2)/r
          x(3) = x(3)/r
          dylm(ig) = c/7.0_real_8*( &
               4.0_real_8*x(2)*x(3)*(7.0_real_8*x(1)**2-1.0_real_8)*gagk(ig,kk)/r2 &
              -x(alpha(kk))*( &
              +12.0_real_8*x(1)*x(2)*x(3)*delta(1,beta(kk)) &
              +x(3)*(7.0_real_8*x(1)**2-2.0_real_8*x(2)**2-1.0_real_8)*delta(2,beta(kk)) &
              +x(2)*(7.0_real_8*x(1)**2-2.0_real_8*x(3)**2-1.0_real_8)*delta(3,beta(kk))))
       ENDDO
    ELSE IF (l.EQ.25) THEN
       c=SQRT(9._real_8*5._real_8/fpi/3._real_8)*7._real_8/2._real_8
       !$OMP parallel do private(IG,X,R,R2) shared(DYLM)
       DO ig=is1,ncpw%nhg
          x(1) = gk(1,ig)
          x(2) = gk(2,ig)
          x(3) = gk(3,ig)
          r = MAX(SQRT(x(1)*x(1)+x(2)*x(2)+x(3)*x(3)),divcut)
          r2 = r*r*parm%tpiba2
          x(1) = x(1)/r
          x(2) = x(2)/r
          x(3) = x(3)/r
          dylm(ig) = c/7.0_real_8*( &
               4.0_real_8*(7.0_real_8*x(3)**4-3.5_real_8*(x(1)**4+x(2)**4) &
              -6.0_real_8*(x(3)**2-0.5_real_8*(x(1)**2+x(2)**2)))*gagk(ig,kk)/r2 &
              -x(alpha(kk))*( &
               x(1)*(-14._real_8*x(1)**2-18._real_8*x(3)**2+12.0_real_8)*delta(1,beta(kk)) &
              +x(2)*(-14._real_8*x(2)**2-18._real_8*x(3)**2+12.0_real_8)*delta(2,beta(kk)) &
              +x(3)*(28.0_real_8*x(3)**2-18._real_8*x(3)**2-6.0_real_8)*delta(3,beta(kk))))
       ENDDO
    ELSE IF (l.GE.26) THEN
       CALL stopgm('DYLMR',' HIGHER L NOT PROGRAMMED  ',&
            __LINE__,__FILE__)
    ENDIF
  
    RETURN
  END SUBROUTINE dylmr
  ! ==================================================================
  SUBROUTINE dbess(xg,l,mmax,r,djl)
    ! ==--------------------------------------------------------------==
    ! CALCULATES DERIVATIVES OF SPHERICAL BESSEL FUNCTIONS  j_l(Gr)
    ! WITH RESPECT TO h_alpha,beta (WITHOUT THE FACTOR GAGK(KK,IG)*HTM1)
    ! NAMELY, WE COMPUTE djl(x) = -x * d(j_l(x)/dx)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: xg
    INTEGER                                  :: l, mmax
    REAL(real_8)                             :: r(*), djl(mmax)

    REAL(real_8), PARAMETER                  :: eps = 1.e-8_real_8 

    INTEGER                                  :: ir
    REAL(real_8)                             :: xrg

! ==--------------------------------------------------------------==

    IF (l.EQ.1) THEN                      ! S  PART
       IF (ABS(xg).LT.eps) THEN
          CALL zeroing(djl)!,mmax)
       ELSE
          djl(1) = 0._real_8
          !$omp parallel do private(IR,XRG) shared(DJL)
          DO ir=2,mmax
             xrg=r(ir)*xg
             djl(ir) = SIN(xrg)/xrg-COS(xrg)
          ENDDO
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (l.EQ.2) THEN                      ! P  PART
       IF (ABS(xg).LT.eps) THEN
          CALL zeroing(djl)!,mmax)
       ELSE
          djl(1) = 0._real_8
          !$omp parallel do private(IR,XRG) shared(DJL)
          DO ir=2,mmax
             xrg=r(ir)*xg
             djl(ir) = 2._real_8*(SIN(xrg)/xrg-COS(xrg))/xrg - SIN(xrg)
          ENDDO
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (l.EQ.3) THEN                      ! D  PART
       IF (ABS(xg).LT.eps) THEN
          CALL zeroing(djl)!,mmax)
       ELSE
          djl(1) = 0._real_8
          !$omp parallel do private(IR,XRG) shared(DJL)
          DO ir=2,mmax
             xrg=r(ir)*xg
             djl(ir) = ( SIN(xrg)*(9._real_8/(xrg*xrg)-4._real_8)  &
                     - 9._real_8*COS(xrg)/xrg ) /xrg + COS(xrg)
          ENDDO
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (l.EQ.4) THEN                      ! F  PART
       IF (ABS(xg).LT.eps) THEN
          CALL zeroing(djl)!,mmax)
       ELSE
          djl(1) = 0._real_8
          !$omp parallel do private(IR,XRG) shared(DJL)
          DO ir=2,mmax
            xrg=r(ir)*xg
            djl(ir) = SIN(xrg)*(60._real_8/(xrg*xrg)-27._real_8)/(xrg*xrg) &
                    - COS(xrg)*(60._real_8/(xrg*xrg)-7._real_8)/xrg        &
                    + SIN(xrg)        
          ENDDO
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (l.EQ.5) THEN                      ! G  PART
       IF (ABS(xg).LT.eps) THEN
          CALL zeroing(djl)!,mmax)
       ELSE
          djl(1) = 0._real_8
          !$omp parallel do private(IR,XRG) shared(DJL)
          DO ir=2,mmax
            xrg=r(ir)*xg
            djl(ir) = SIN(xrg)*(525._real_8/(xrg*xrg*xrg*xrg)-240._real_8/(xrg*xrg)+11._real_8)/xrg &
                    - COS(xrg)*(525._real_8/(xrg*xrg)-65._real_8)/(xrg*xrg) &
                    - COS(xrg)        
          ENDDO
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (l.LE.0 .OR. l.GE.6) THEN
       CALL stopgm('DBESS',' L NOT PROGRAMMED  ',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE dbess
  ! ==================================================================
  FUNCTION dylmrx(l,ig,gk,kk)
    ! ==--------------------------------------------------------------==
    ! REAL SPHERICAL HARMONICS,  L IS COMBINED INDEX FOR LM  (L=1,2...9)
    ! ORDER:  S, P_X, P_Z, P_Y, D_XY, D_XZ, D_Z^2, D_YZ, D_X^2-Y^2
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: l, ig
    REAL(real_8)                             :: gk(3,ncpw%nhg)
    INTEGER                                  :: kk
    REAL(real_8)                             :: dylmrx

    REAL(real_8)                             :: c, r, r2, x(3)

    IF (l.GE.26) THEN
       CALL stopgm('DYLMR',' L NOT PROGRAMMED  ',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    x(1) = gk(1,ig)
    x(2) = gk(2,ig)
    x(3) = gk(3,ig)
    r = SQRT(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
    IF (r.LT.1.e-6_real_8) THEN
       dylmrx=0.0_real_8
       RETURN
    ENDIF
    r2 = r*r*parm%tpiba2
    x(1) = x(1)/r
    x(2) = x(2)/r
    x(3) = x(3)/r
    ! ==--------------------------------------------------------------==
    IF (l.EQ.1) THEN
       dylmrx = 0._real_8
    ELSE IF (l.EQ.2) THEN
       dylmrx = SQRT(3._real_8/fpi)*( x(1)*gagk(ig,kk)/r2 -&
            x(alpha(kk))*delta(1,beta(kk)) )
    ELSE IF (l.EQ.3) THEN
       dylmrx = SQRT(3._real_8/fpi)*( x(3)*gagk(ig,kk)/r2 -&
            x(alpha(kk))*delta(3,beta(kk)) )
    ELSE IF (l.EQ.4) THEN
       dylmrx = SQRT(3._real_8/fpi)*( x(2)*gagk(ig,kk)/r2 -&
            x(alpha(kk))*delta(2,beta(kk)) )
    ELSE IF (l.EQ.5) THEN
       dylmrx = SQRT(15._real_8/fpi)*( 2._real_8*x(1)*x(2)*gagk(ig,kk)/r2 -&
            x(alpha(kk))*x(2)*delta(1,beta(kk)) -&
            x(alpha(kk))*x(1)*delta(2,beta(kk)) )
    ELSE IF (l.EQ.6) THEN
       dylmrx = SQRT(15._real_8/fpi)*( 2._real_8*x(1)*x(3)*gagk(ig,kk)/r2 -&
            x(alpha(kk))*x(3)*delta(1,beta(kk)) -&
            x(alpha(kk))*x(1)*delta(3,beta(kk)) )
    ELSE IF (l.EQ.7) THEN
       dylmrx = 6._real_8*SQRT(5._real_8/fpi/4._real_8)*( x(3)*x(3)*gagk(ig,kk)/r2 -&
            x(alpha(kk))*x(3)*delta(3,beta(kk)) )
    ELSE IF (l.EQ.8) THEN
       dylmrx = SQRT(15._real_8/fpi)*(2._real_8*x(2)*x(3)*gagk(ig,kk)/r2 -&
            x(alpha(kk))*x(3)*delta(2,beta(kk)) -&
            x(alpha(kk))*x(2)*delta(3,beta(kk)) )
    ELSE IF (l.EQ.9) THEN
       dylmrx = 2._real_8*SQRT(15._real_8/fpi/4._real_8)*( (x(1)*x(1)-x(2)*x(2))*&
            gagk(ig,kk)/r2 - x(alpha(kk))*x(1)*delta(1,beta(kk)) +&
            x(alpha(kk))*x(2)*delta(2,beta(kk)) )
    ELSE IF (l.EQ.10) THEN
       c=SQRT(7._real_8/fpi)*5._real_8/2._real_8
       dylmrx = c*( &
            3.0_real_8*x(1)*(x(1)*x(1)-0.6_real_8)*gagk(ig,kk)/r2 &
            -0.6_real_8*x(alpha(kk))*( &
            (3.0_real_8*x(1)*x(1)-1.0_real_8)*delta(1,beta(kk)) &
            -2.0_real_8*x(1)*x(2)*delta(2,beta(kk)) &
            -2.0_real_8*x(1)*x(3)*delta(3,beta(kk))))
    ELSE IF (l.EQ.11) THEN
       c=SQRT(7._real_8/fpi)*5._real_8/2._real_8
       dylmrx = c*( &
            3.0_real_8*x(2)*(x(2)*x(2)-0.6_real_8)*gagk(ig,kk)/r2 &
            -0.6_real_8*x(alpha(kk))*( &
            -2.0_real_8*x(2)*x(1)*delta(1,beta(kk)) &
            +(3.0_real_8*x(2)*x(2)-1.0_real_8)*delta(2,beta(kk)) &
            -2.0_real_8*x(2)*x(3)*delta(3,beta(kk))))
    ELSE IF (l.EQ.12) THEN
       c=SQRT(7._real_8*15._real_8/fpi)
       dylmrx = c*( &
            3.0_real_8*x(1)*x(2)*x(3)*gagk(ig,kk)/r2 &
            -x(alpha(kk))*( &
            x(2)*x(3)*delta(1,beta(kk)) &
            +x(1)*x(3)*delta(2,beta(kk)) &
            +x(1)*x(2)*delta(3,beta(kk))))
    ELSE IF (l.EQ.13) THEN
       c=SQRT(7._real_8/fpi)*5._real_8/2._real_8
       dylmrx = c*( &
            3.0_real_8*x(3)*(x(3)*x(3)-0.6_real_8)*gagk(ig,kk)/r2 &
            -0.6_real_8*x(alpha(kk))*( &
            -2.0_real_8*x(3)*x(1)*delta(1,beta(kk)) &
            -2.0_real_8*x(3)*x(2)*delta(2,beta(kk)) &
            +(3.0_real_8*x(3)*x(3)-1.0_real_8)*delta(3,beta(kk))))
    ELSE IF (l.EQ.14) THEN
       c=SQRT(7._real_8*15._real_8/fpi)/2._real_8
       dylmrx = c*( &
            3.0_real_8*x(3)*(x(1)*x(1)-x(2)*x(2))*gagk(ig,kk)/r2 &
            -x(alpha(kk))*( &
            2.0_real_8*x(3)*x(1)*delta(1,beta(kk)) &
            -2.0_real_8*x(3)*x(2)*delta(2,beta(kk)) &
            +(x(1)*x(1)-x(2)*x(2))*delta(3,beta(kk))))
    ELSE IF (l.EQ.15) THEN
       c=SQRT(7._real_8*15._real_8/fpi)/2._real_8
       dylmrx = c*( &
            3.0_real_8*x(2)*(x(3)*x(3)-x(1)*x(1))*gagk(ig,kk)/r2 &
            -x(alpha(kk))*( &
            -2.0_real_8*x(2)*x(1)*delta(1,beta(kk)) &
            +(x(3)*x(3)-x(1)*x(1))*delta(2,beta(kk)) &
            +2.0_real_8*x(2)*x(3)*delta(3,beta(kk))))
    ELSE IF (l.EQ.16) THEN
       c=SQRT(7._real_8*15._real_8/fpi)/2._real_8
       dylmrx = c*( &
            3.0_real_8*x(1)*(x(2)*x(2)-x(3)*x(3))*gagk(ig,kk)/r2 &
            -x(alpha(kk))*( &
            (x(2)*x(2)-x(3)*x(3))*delta(1,beta(kk)) &
            +2.0_real_8*x(1)*x(2)*delta(2,beta(kk)) &
            -2.0_real_8*x(1)*x(3)*delta(3,beta(kk))))
    ELSE IF (l.EQ.17) THEN
       c=SQRT(3._real_8*7._real_8/fpi)*5._real_8/4._real_8
       dylmrx = c*( &
            4.0_real_8*(x(1)**4+x(2)**4+x(3)**4)*gagk(ig,kk)/r2 &
            -4.0_real_8*x(alpha(kk))*( &
            x(1)**3*delta(1,beta(kk)) &
            +x(2)**3*delta(2,beta(kk)) &
            +x(3)**3*delta(3,beta(kk))))
    ELSE IF (l.EQ.18) THEN
       c=SQRT(9._real_8*35._real_8/fpi)/2._real_8
       dylmrx = c*( &
            4.0_real_8*x(2)*x(3)*(x(2)*x(2)-x(3)*x(3))*gagk(ig,kk)/r2 &
            -x(alpha(kk))*( &
            (3.0_real_8*x(2)*x(2)-x(3)*x(3))*x(3)*delta(2,beta(kk)) &
            -(3.0_real_8*x(3)*x(3)-x(2)*x(2))*x(2)*delta(3,beta(kk))))
    ELSE IF (l.EQ.19) THEN
       c=SQRT(9._real_8*35._real_8/fpi)/2._real_8
       dylmrx = c*( &
            4.0_real_8*x(1)*x(3)*(x(3)*x(3)-x(1)*x(1))*gagk(ig,kk)/r2 &
            -x(alpha(kk))*( &
            -(3.0_real_8*x(1)*x(1)-x(3)*x(3))*x(3)*delta(1,beta(kk)) &
            +(3.0_real_8*x(3)*x(3)-x(1)*x(1))*x(1)*delta(3,beta(kk))))
    ELSE IF (l.EQ.20) THEN
       c=SQRT(9._real_8*5._real_8/fpi)/4._real_8
       dylmrx = c*( &
            4.0_real_8*((x(1)**4-x(2)**4)-6.0_real_8*x(3)**2*(x(1)**2-x(2)**2)) &
            *gagk(ig,kk)/r2 -x(alpha(kk))*( &
            (4.0_real_8*x(1)**3-12.0_real_8*x(3)**2*x(1))*delta(1,beta(kk)) &
            -(4.0_real_8*x(2)**3-12.0_real_8*x(3)**2*x(2))*delta(2,beta(kk)) &
            -12.0_real_8*x(3)*(x(1)**2-x(2)**2)*delta(3,beta(kk))))
    ELSE IF (l.EQ.21) THEN
       c=SQRT(9._real_8*35._real_8/fpi)/2._real_8
       dylmrx = c*( &
            4.0_real_8*x(1)*x(2)*(x(1)*x(1)-x(2)*x(2))*gagk(ig,kk)/r2 &
            -x(alpha(kk))*( &
            +(3.0_real_8*x(1)*x(1)-x(2)*x(2))*x(2)*delta(1,beta(kk)) &
            -(3.0_real_8*x(2)*x(2)-x(1)*x(1))*x(1)*delta(2,beta(kk))))
    ELSE IF (l.EQ.22) THEN
       c=SQRT(9._real_8*5._real_8/fpi)*7._real_8/2._real_8
       dylmrx = c/7.0_real_8*( &
            4.0_real_8*x(1)*x(2)*(7.0_real_8*x(3)**2-1.0_real_8)*gagk(ig,kk)/r2 &
            -x(alpha(kk))*( &
            x(2)*(7.0_real_8*x(3)**2-2.0_real_8*x(1)**2-1.0_real_8)*delta(1,beta(kk)) &
            +x(1)*(7.0_real_8*x(3)**2-2.0_real_8*x(2)**2-1.0_real_8)*delta(2,beta(kk)) &
            +12.0_real_8*x(1)*x(2)*x(3)*delta(3,beta(kk))))
    ELSE IF (l.EQ.23) THEN
       c=SQRT(9._real_8*5._real_8/fpi)*7._real_8/2._real_8
       dylmrx = c/7.0_real_8*( &
            4.0_real_8*x(3)*x(1)*(7.0_real_8*x(2)**2-1.0_real_8)*gagk(ig,kk)/r2 &
            -x(alpha(kk))*( &
            x(3)*(7.0_real_8*x(2)**2-2.0_real_8*x(1)**2-1.0_real_8)*delta(1,beta(kk)) &
            +12.0_real_8*x(1)*x(2)*x(3)*delta(2,beta(kk)) &
            +x(1)*(7.0_real_8*x(2)**2-2.0_real_8*x(3)**2-1.0_real_8)*delta(3,beta(kk))))
    ELSE IF (l.EQ.24) THEN
       c=SQRT(9._real_8*5._real_8/fpi)*7._real_8/2._real_8
       dylmrx = c/7.0_real_8*( &
            4.0_real_8*x(2)*x(3)*(7.0_real_8*x(1)**2-1.0_real_8)*gagk(ig,kk)/r2 &
            -x(alpha(kk))*( &
            +12.0_real_8*x(1)*x(2)*x(3)*delta(1,beta(kk)) &
            +x(3)*(7.0_real_8*x(1)**2-2.0_real_8*x(2)**2-1.0_real_8)*delta(2,beta(kk)) &
            +x(2)*(7.0_real_8*x(1)**2-2.0_real_8*x(3)**2-1.0_real_8)*delta(3,beta(kk))))
    ELSE IF (l.EQ.25) THEN
       c=SQRT(9._real_8*5._real_8/fpi/3._real_8)*7._real_8/2._real_8
       dylmrx = c/7.0_real_8*( &
            4.0_real_8*(7.0_real_8*x(3)**4-3.5_real_8*(x(1)**4+x(2)**4) &
            -6.0_real_8*(x(3)**2-0.5_real_8*(x(1)**2+x(2)**2)))*gagk(ig,kk)/r2 &
            -x(alpha(kk))*( &
            x(1)*(-14.0_real_8*x(1)**2-18._real_8*x(3)**2+12.0_real_8)*delta(1,beta(kk)) &
            +x(2)*(-14.0_real_8*x(2)**2-18._real_8*x(3)**2+12.0_real_8)*delta(2,beta(kk)) &
            +x(3)*( 28.0_real_8*x(3)**2-18._real_8*x(3)**2- 6.0_real_8)*delta(3,beta(kk))))
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION dylmrx
  ! ==================================================================
  FUNCTION dylmrkx(l,ig,gk,rk,gagkk,kk,iplus)
    ! ==--------------------------------------------------------------==
    ! == REAL SPHERICAL HARMONICS,  L IS COMBINED INDEX FOR LM        ==
    ! == (L=1,2...9) ORDER:                                           ==
    ! == S, P_X, P_Z, P_Y, D_XY, D_XZ, D_Z^2, D_YZ, D_X^2-Y^2         ==
    ! ==--------------------------------------------------------------==
    ! == IG: components number in GK and GAGK array                   ==
    ! == RK: k points components                                      ==
    ! == KK: from 1 to 6                                              ==
    ! == IPLUS: -1 or 1 (RK+GK) or (RK-GK)                            ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: l, ig
    REAL(real_8)                             :: gk(3,ncpw%ngw), rk(3), &
                                                gagkk(ncpw%ngw,6)
    INTEGER                                  :: kk, iplus
    REAL(real_8)                             :: dylmrkx

    REAL(real_8)                             :: c, r, r2, x(3)

    IF (l.GE.26) THEN
       CALL stopgm('DYLMKRX',' L NOT PROGRAMMED  ',& 
            __LINE__,__FILE__)
    ENDIF
    IF (iplus.EQ.1) THEN
       x(1) = rk(1)+gk(1,ig)
       x(2) = rk(2)+gk(2,ig)
       x(3) = rk(3)+gk(3,ig)
    ELSEIF (iplus.EQ.-1) THEN
       x(1) = rk(1)-gk(1,ig)
       x(2) = rk(2)-gk(2,ig)
       x(3) = rk(3)-gk(3,ig)
    ELSE
       CALL stopgm(' DYLMRKX','IPLUS WRONG',& 
            __LINE__,__FILE__)
    ENDIF
    r = SQRT(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
    IF (r.LT.1.e-6_real_8) THEN
       dylmrkx=0.0_real_8
       RETURN
    ENDIF
    r2 = r*r*parm%tpiba2
    x(1) = x(1)/r
    x(2) = x(2)/r
    x(3) = x(3)/r

    IF (l.EQ.1) THEN
       dylmrkx = 0._real_8
    ELSE IF (l.EQ.2) THEN
       dylmrkx = SQRT(3._real_8/fpi)*( x(1)*gagkk(ig,kk)/r2 -&
            x(alpha(kk))*delta(1,beta(kk)) )
    ELSE IF (l.EQ.3) THEN
       dylmrkx = SQRT(3._real_8/fpi)*( x(3)*gagkk(ig,kk)/r2 -&
            x(alpha(kk))*delta(3,beta(kk)) )
    ELSE IF (l.EQ.4) THEN
       dylmrkx = SQRT(3._real_8/fpi)*( x(2)*gagkk(ig,kk)/r2 -&
            x(alpha(kk))*delta(2,beta(kk)) )
    ELSE IF (l.EQ.5) THEN
       dylmrkx = SQRT(15._real_8/fpi)*( 2._real_8*x(1)*x(2)*gagkk(ig,kk)/r2 -&
            x(alpha(kk))*x(2)*delta(1,beta(kk)) -&
            x(alpha(kk))*x(1)*delta(2,beta(kk)) )
    ELSE IF (l.EQ.6) THEN
       dylmrkx = SQRT(15._real_8/fpi)*( 2._real_8*x(1)*x(3)*gagkk(ig,kk)/r2 -&
            x(alpha(kk))*x(3)*delta(1,beta(kk)) -&
            x(alpha(kk))*x(1)*delta(3,beta(kk)) )
    ELSE IF (l.EQ.7) THEN
       dylmrkx = 6._real_8*SQRT(5._real_8/fpi/4._real_8)*( x(3)*x(3)*gagkk(ig,kk)/r2&
            - X(ALPHA(KK))*X(3)*DELTA(3,BETA(KK)) )
    ELSE IF (l.EQ.8) THEN
       dylmrkx = SQRT(15._real_8/fpi)*(2._real_8*x(2)*x(3)*gagkk(ig,kk)/r2 -&
            x(alpha(kk))*x(3)*delta(2,beta(kk)) -&
            x(alpha(kk))*x(2)*delta(3,beta(kk)) )
    ELSE IF (l.EQ.9) THEN
       dylmrkx = 2._real_8*SQRT(15._real_8/fpi/4._real_8)*( (x(1)*x(1)-x(2)*x(2))*&
            gagkk(ig,kk)/r2 - x(alpha(kk))*x(1)*delta(1,beta(kk)) +&
            x(alpha(kk))*x(2)*delta(2,beta(kk)) )
    ELSE IF (l.EQ.10) THEN
       c=SQRT(7._real_8/fpi)*5._real_8/2._real_8
       dylmrkx = c*( &
            3.0_real_8*x(1)*(x(1)*x(1)-0.6_real_8)*gagkk(ig,kk)/r2 &
            -0.6_real_8*x(alpha(kk))*( &
            (3.0_real_8*x(1)*x(1)-1.0_real_8)*delta(1,beta(kk)) &
            -2.0_real_8*x(1)*x(2)*delta(2,beta(kk)) &
            -2.0_real_8*x(1)*x(3)*delta(3,beta(kk))))
    ELSE IF (l.EQ.11) THEN
       c=SQRT(7._real_8/fpi)*5._real_8/2._real_8
       dylmrkx = c*( &
            3.0_real_8*x(2)*(x(2)*x(2)-0.6_real_8)*gagkk(ig,kk)/r2 &
            -0.6_real_8*x(alpha(kk))*( &
            -2.0_real_8*x(2)*x(1)*delta(1,beta(kk)) &
            +(3.0_real_8*x(2)*x(2)-1.0_real_8)*delta(2,beta(kk)) &
            -2.0_real_8*x(2)*x(3)*delta(3,beta(kk))))
    ELSE IF (l.EQ.12) THEN
       c=SQRT(7._real_8*15._real_8/fpi)
       dylmrkx = c*( &
            3.0_real_8*x(1)*x(2)*x(3)*gagkk(ig,kk)/r2 &
            -x(alpha(kk))*( &
            x(2)*x(3)*delta(1,beta(kk)) &
            +x(1)*x(3)*delta(2,beta(kk)) &
            +x(1)*x(2)*delta(3,beta(kk))))
    ELSE IF (l.EQ.13) THEN
       c=SQRT(7._real_8/fpi)*5._real_8/2._real_8
       dylmrkx = c*( &
            3.0_real_8*x(3)*(x(3)*x(3)-0.6_real_8)*gagkk(ig,kk)/r2 &
            -0.6_real_8*x(alpha(kk))*( &
            -2.0_real_8*x(3)*x(1)*delta(1,beta(kk)) &
            -2.0_real_8*x(3)*x(2)*delta(2,beta(kk)) &
            +(3.0_real_8*x(3)*x(3)-1.0_real_8)*delta(3,beta(kk))))
    ELSE IF (l.EQ.14) THEN
       c=SQRT(7._real_8*15._real_8/fpi)/2._real_8
       dylmrkx = c*( &
            3.0_real_8*x(3)*(x(1)*x(1)-x(2)*x(2))*gagkk(ig,kk)/r2 &
            -x(alpha(kk))*( &
            2.0_real_8*x(3)*x(1)*delta(1,beta(kk)) &
            -2.0_real_8*x(3)*x(2)*delta(2,beta(kk)) &
            +(x(1)*x(1)-x(2)*x(2))*delta(3,beta(kk))))
    ELSE IF (l.EQ.15) THEN
       c=SQRT(7._real_8*15._real_8/fpi)/2._real_8
       dylmrkx = c*( &
            3.0_real_8*x(2)*(x(3)*x(3)-x(1)*x(1))*gagkk(ig,kk)/r2 &
            -x(alpha(kk))*( &
            -2.0_real_8*x(2)*x(1)*delta(1,beta(kk)) &
            +(x(3)*x(3)-x(1)*x(1))*delta(2,beta(kk)) &
            +2.0_real_8*x(2)*x(3)*delta(3,beta(kk))))
    ELSE IF (l.EQ.16) THEN
       c=SQRT(7._real_8*15._real_8/fpi)/2._real_8
       dylmrkx = c*( &
            3.0_real_8*x(1)*(x(2)*x(2)-x(3)*x(3))*gagkk(ig,kk)/r2 &
            -x(alpha(kk))*( &
            (x(2)*x(2)-x(3)*x(3))*delta(1,beta(kk)) &
            +2.0_real_8*x(1)*x(2)*delta(2,beta(kk)) &
            -2.0_real_8*x(1)*x(3)*delta(3,beta(kk))))
    ELSE IF (l.EQ.17) THEN
       c=SQRT(3._real_8*7._real_8/fpi)*5._real_8/4._real_8
       dylmrkx = c*( &
            4.0_real_8*(x(1)**4+x(2)**4+x(3)**4)*gagkk(ig,kk)/r2 &
            -4.0_real_8*x(alpha(kk))*( &
            x(1)**3*delta(1,beta(kk)) &
            +x(2)**3*delta(2,beta(kk)) &
            +x(3)**3*delta(3,beta(kk))))
    ELSE IF (l.EQ.18) THEN
       c=SQRT(9._real_8*35._real_8/fpi)/2._real_8
       dylmrkx = c*( &
            4.0_real_8*x(2)*x(3)*(x(2)*x(2)-x(3)*x(3))*gagkk(ig,kk)/r2 &
            -x(alpha(kk))*( &
            (3.0_real_8*x(2)*x(2)-x(3)*x(3))*x(3)*delta(2,beta(kk)) &
            -(3.0_real_8*x(3)*x(3)-x(2)*x(2))*x(2)*delta(3,beta(kk))))
    ELSE IF (l.EQ.19) THEN
       c=SQRT(9._real_8*35._real_8/fpi)/2._real_8
       dylmrkx = c*( &
            4.0_real_8*x(1)*x(3)*(x(3)*x(3)-x(1)*x(1))*gagkk(ig,kk)/r2 &
            -x(alpha(kk))*( &
            -(3.0_real_8*x(1)*x(1)-x(3)*x(3))*x(3)*delta(1,beta(kk)) &
            +(3.0_real_8*x(3)*x(3)-x(1)*x(1))*x(1)*delta(3,beta(kk))))
    ELSE IF (l.EQ.20) THEN
       c=SQRT(9._real_8*5._real_8/fpi)/4._real_8
       dylmrkx = c*( &
            4.0_real_8*((x(1)**4-x(2)**4)-6.0_real_8*x(3)**2*(x(1)**2-x(2)**2)) &
            *gagkk(ig,kk)/r2 -x(alpha(kk))*( &
            (4.0_real_8*x(1)**3-12.0_real_8*x(3)**2*x(1))*delta(1,beta(kk)) &
            -(4.0_real_8*x(2)**3-12.0_real_8*x(3)**2*x(2))*delta(2,beta(kk)) &
            -12.0_real_8*x(3)*(x(1)**2-x(2)**2)*delta(3,beta(kk))))
    ELSE IF (l.EQ.21) THEN
       c=SQRT(9._real_8*35._real_8/fpi)/2._real_8
       dylmrkx = c*( &
            4.0_real_8*x(1)*x(2)*(x(1)*x(1)-x(2)*x(2))*gagkk(ig,kk)/r2 &
            -x(alpha(kk))*( &
            +(3.0_real_8*x(1)*x(1)-x(2)*x(2))*x(2)*delta(1,beta(kk)) &
            -(3.0_real_8*x(2)*x(2)-x(1)*x(1))*x(1)*delta(2,beta(kk))))
    ELSE IF (l.EQ.22) THEN
       c=SQRT(9._real_8*5._real_8/fpi)*7._real_8/2._real_8
       dylmrkx = c/7.0_real_8*( &
            4.0_real_8*x(1)*x(2)*(7.0_real_8*x(3)**2-1.0_real_8)*gagkk(ig,kk)/r2 &
            -x(alpha(kk))*( &
            x(2)*(7.0_real_8*x(3)**2-2.0_real_8*x(1)**2-1.0_real_8)*delta(1,beta(kk)) &
            +x(1)*(7.0_real_8*x(3)**2-2.0_real_8*x(2)**2-1.0_real_8)*delta(2,beta(kk)) &
            +12.0_real_8*x(1)*x(2)*x(3)*delta(3,beta(kk))))
    ELSE IF (l.EQ.23) THEN
       c=SQRT(9._real_8*5._real_8/fpi)*7._real_8/2._real_8
       dylmrkx = c/7.0_real_8*( &
            4.0_real_8*x(3)*x(1)*(7.0_real_8*x(2)**2-1.0_real_8)*gagkk(ig,kk)/r2 &
            -x(alpha(kk))*( &
            x(3)*(7.0_real_8*x(2)**2-2.0_real_8*x(1)**2-1.0_real_8)*delta(1,beta(kk)) &
            +12.0_real_8*x(1)*x(2)*x(3)*delta(2,beta(kk)) &
            +x(1)*(7.0_real_8*x(2)**2-2.0_real_8*x(3)**2-1.0_real_8)*delta(3,beta(kk))))
    ELSE IF (l.EQ.24) THEN
       c=SQRT(9._real_8*5._real_8/fpi)*7._real_8/2._real_8
       dylmrkx = c/7.0_real_8*( &
            4.0_real_8*x(2)*x(3)*(7.0_real_8*x(1)**2-1.0_real_8)*gagkk(ig,kk)/r2 &
            -x(alpha(kk))*( &
            +12.0_real_8*x(1)*x(2)*x(3)*delta(1,beta(kk)) &
            +x(3)*(7.0_real_8*x(1)**2-2.0_real_8*x(2)**2-1.0_real_8)*delta(2,beta(kk)) &
            +x(2)*(7.0_real_8*x(1)**2-2.0_real_8*x(3)**2-1.0_real_8)*delta(3,beta(kk))))
    ELSE IF (l.EQ.25) THEN
       c=SQRT(9._real_8*5._real_8/fpi/3._real_8)*7._real_8/2._real_8
       dylmrkx = c/7.0_real_8*( &
            4.0_real_8*(7.0_real_8*x(3)**4-3.5_real_8*(x(1)**4+x(2)**4) &
            -6.0_real_8*(x(3)**2-0.5_real_8*(x(1)**2+x(2)**2)))*gagkk(ig,kk)/r2 &
            -x(alpha(kk))*( &
            x(1)*(-14.0_real_8*x(1)**2-18._real_8*x(3)**2+12.0_real_8)*delta(1,beta(kk)) &
            +x(2)*(-14.0_real_8*x(2)**2-18._real_8*x(3)**2+12.0_real_8)*delta(2,beta(kk)) &
            +x(3)*( 28.0_real_8*x(3)**2-18._real_8*x(3)**2- 6.0_real_8)*delta(3,beta(kk))))
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION dylmrkx
  ! ==--------------------------------------------------------------==


END MODULE dylmr_utils
