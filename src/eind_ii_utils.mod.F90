MODULE eind_ii_utils
  USE cnst,                            ONLY: pi
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE metr,                            ONLY: metr_com
  USE parac,                           ONLY: paral
  USE pbc_utils,                       ONLY: pbc
  USE ragg,                            ONLY: raggio
  USE special_functions,               ONLY: cp_erfc
  USE system,                          ONLY: parm
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: eind_ii

CONTAINS

  ! ==================================================================
  SUBROUTINE eind_ii(tau0,is,ia,ka,eind,iesr)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! == CONTRIBUTION OF ESR TO THE SECOND ORDER ENERGY               ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:)
    INTEGER                                  :: is, ia, ka
    REAL(real_8)                             :: eind
    INTEGER                                  :: iesr

    REAL(real_8), PARAMETER                  :: argmax = 20._real_8 

    INTEGER                                  :: ishft, isub, ix, iy, iz, j, &
                                                k, l, m
    LOGICAL                                  :: tzero
    REAL(real_8)                             :: aarg, arg, dfd, dsd, dspi, &
                                                earg, erre2, esrtzero, rckj, &
                                                rlm, rxlm(3), xlm, ylm, zlm, &
                                                zv2, xlm_, ylm_, zlm_

! ==--------------------------------------------------------------==

    CALL tiset('   EIND_II',isub)
    dspi=SQRT(pi)
    IF (paral%parent) THEN
       k=is
       DO j=1,ions1%nsp
          zv2=ions0%zv(k)*ions0%zv(j)
          rckj=SQRT(raggio(k)*raggio(k)+raggio(j)*raggio(j))
          l=ia
          DO m=1,ions0%na(j)
             IF (l.EQ.m.AND.k.EQ.j) THEN
                xlm=0._real_8
                ylm=0._real_8
                zlm=0._real_8
                tzero=.TRUE.
             ELSE
                tzero=.FALSE.
                xlm_=tau0(1,l,k)-tau0(1,m,j)
                ylm_=tau0(2,l,k)-tau0(2,m,j)
                zlm_=tau0(3,l,k)-tau0(3,m,j)
                CALL pbc(xlm_,ylm_,zlm_,xlm,ylm,zlm,1,parm%apbc,parm%ibrav)
             ENDIF
             DO ix=-iesr,iesr
                DO iy=-iesr,iesr
                   DO iz=-iesr,iesr
                      ishft=ix*ix+iy*iy+iz*iz
                      IF (.NOT.(tzero.AND.ishft.EQ.0)) THEN
                         rxlm(1)=xlm+ix*metr_com%ht(1,1)+iy*metr_com%ht(2,1)+iz*metr_com%ht(3,1)
                         rxlm(2)=ylm+ix*metr_com%ht(1,2)+iy*metr_com%ht(2,2)+iz*metr_com%ht(3,2)
                         rxlm(3)=zlm+ix*metr_com%ht(1,3)+iy*metr_com%ht(2,3)+iz*metr_com%ht(3,3)
                         erre2=rxlm(1)**2+rxlm(2)**2+rxlm(3)**2
                         rlm=SQRT(erre2)
                         arg=rlm/rckj
                         IF (arg.LE.argmax) THEN! ADDESR,ADDPRE /= 0
                            IF (tzero) THEN
                               esrtzero=0.5_real_8*2._real_8*zv2/erre2
                            ELSE
                               esrtzero=1._real_8*2._real_8*zv2/erre2
                            ENDIF
                            earg=cp_erfc(arg)/(2._real_8*rlm)
                            aarg=EXP(-arg*arg)/(dspi*rckj)
                            dfd=-esrtzero*(earg + aarg)
                            dsd=esrtzero*2._real_8*(3._real_8*earg/(2._real_8*erre2)+&
                                 aarg*(1.5_real_8/erre2+1._real_8/(rckj*rckj)))
                            eind=eind+dfd+dsd*rxlm(ka)*rxlm(ka)
                         ENDIF
                      ENDIF
                   ENDDO! ixc
                ENDDO! iyc
             ENDDO ! izc
          ENDDO
       ENDDO
    ENDIF
    CALL tihalt('   EIND_II',isub)
    ! ==================================================================
    RETURN
  END SUBROUTINE eind_ii
  ! ==================================================================

END MODULE eind_ii_utils
