MODULE rpiiint_utils
  USE cnst,                            ONLY: pi
  USE eam,                             ONLY: tieam
  USE eam_pot_utils,                   ONLY: eam_pot
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE metr,                            ONLY: metr_com
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai,&
                                             paral
  USE pbc_utils,                       ONLY: pbc
  USE ragg,                            ONLY: raggio
  USE special_functions,               ONLY: cp_erfc
  USE system,                          ONLY: iatpe,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rpiiint

CONTAINS

  ! ==================================================================
  SUBROUTINE rpiiint(esr,tau0,fion,iesr,tfor)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! == THE RESIDUAL PART OF THE ION-ION INTERACTION (ESR)           ==
    ! == DUE TO THE OVERLAP OF THE SMEARED IONIC CHARGE DENSITIES     ==
    ! == (ionic point charges replaced by gaussian charge distrib.)   ==
    ! == CORRESPONDING TO DIFFERENT ATOMS.                            ==
    ! == ESR depends only on TAU0 (RAGGIO and VALENCE CHARGES)        ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: esr, tau0(:,:,:), fion(:,:,:)
    INTEGER                                  :: iesr
    LOGICAL                                  :: tfor

    REAL(real_8), PARAMETER                  :: argmax = 20._real_8 

    INTEGER                                  :: iat, inf, ishft, isub, ix, &
                                                iy, iz, j, k, l, lax, m
    INTEGER, SAVE                            :: iflag = 0
    LOGICAL                                  :: tzero
    REAL(real_8) :: addesr, addpre, arg, erre2, esrtzero, rckj, repand, rlm, &
      rxlm1, rxlm2, rxlm3, tfion1, tfion2, tfion3, tfion4, tfion5, tfion6, &
      xlm, ylm, zlm, zv2, xlm_, ylm_, zlm_

#ifdef __VECTOR 
    INTEGER :: ixyz,IALL,iesri,iesrj,iesrk
#endif 
    ! ==--------------------------------------------------------------==
    IF (iflag.EQ.0.AND.paral%parent)  THEN
       IF (paral%io_parent)&
            WRITE(6,'(A,T52,I2,A,I2,A,I2,A)')&
            ' EWALD| SUM IN REAL SPACE OVER ',&
            2*iesr+1,'*',2*iesr+1,'*',2*iesr+1,' CELLS'
       iflag=1
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL tiset('   RPIIINT',isub)
#ifdef __VECTOR 
    iesri=(2*iesr+1)
    iesrj=(2*iesr+1)*iesri
    iesrk=(2*iesr+1)*iesrj
    IALL=iesrk
#endif 
    esr=0._real_8
    iat=0
    DO k=1,ions1%nsp
       DO j=k,ions1%nsp
          zv2=ions0%zv(k)*ions0%zv(j)
          IF (ABS(zv2).LT.1.e-10_real_8) GOTO 2000
          rckj=SQRT(raggio(k)*raggio(k)+raggio(j)*raggio(j))
          lax=ions0%na(k)
          DO l=1,lax
             IF (iatpe(iat+l).NE.parai%mepos) GOTO 1000
             inf=1
             IF (k.EQ.j)inf=l
             DO m=inf,ions0%na(j)
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
                IF (tfor) THEN
                   tfion1 = fion(1,l,k)
                   tfion2 = fion(2,l,k)
                   tfion3 = fion(3,l,k)
                   tfion4 = fion(1,m,j)
                   tfion5 = fion(2,m,j)
                   tfion6 = fion(3,m,j)
                ENDIF
#ifdef __VECTOR 
                DO ixyz=1,IALL
                   ix=-iesr+INT(MOD(ixyz-1,iesrk)/iesrj)
                   iy=-iesr+INT(MOD(ixyz-1,iesrj)/iesri)
                   iz=-iesr+MOD(ixyz-1,iesri)
#else 
                   DO ix=-iesr,iesr
                      DO iy=-iesr,iesr
                         DO iz=-iesr,iesr
#endif 
                            ishft=ix*ix+iy*iy+iz*iz
                            IF (.NOT.(tzero.AND.ishft.EQ.0)) THEN
                               rxlm1=xlm+ix*metr_com%ht(1,1)+iy*metr_com%ht(2,1)+iz*metr_com%ht(3,1)
                               rxlm2=ylm+ix*metr_com%ht(1,2)+iy*metr_com%ht(2,2)+iz*metr_com%ht(3,2)
                               rxlm3=zlm+ix*metr_com%ht(1,3)+iy*metr_com%ht(2,3)+iz*metr_com%ht(3,3)
                               erre2=rxlm1*rxlm1+rxlm2*rxlm2+rxlm3*rxlm3
                               rlm=SQRT(erre2)
                               arg=rlm/rckj
                               IF (arg.LE.argmax) THEN! ADDESR,ADDPRE /= 0
                                  IF (tzero) THEN
                                     esrtzero=0.5_real_8
                                  ELSE
                                     esrtzero=1._real_8
                                  ENDIF
                                  addesr=zv2*cp_erfc(arg)/rlm
                                  esr=esr+addesr*esrtzero
                                  IF (tfor) THEN
                                     addpre=(2._real_8*zv2/SQRT(pi))*&
                                          EXP(-arg*arg)/rckj
                                     repand=esrtzero*(addesr+addpre)/erre2
                                     tfion1=tfion1+repand*rxlm1
                                     tfion2=tfion2+repand*rxlm2
                                     tfion3=tfion3+repand*rxlm3
                                     tfion4=tfion4-repand*rxlm1
                                     tfion5=tfion5-repand*rxlm2
                                     tfion6=tfion6-repand*rxlm3
                                  ENDIF
                               ENDIF
                            ENDIF
#ifdef __VECTOR 
                         ENDDO! ixyz
#else 
                      ENDDO! ixc
                   ENDDO! iyc
                ENDDO! izc
#endif 
                IF (tfor) THEN
                   fion(1,l,k) = tfion1
                   fion(2,l,k) = tfion2
                   fion(3,l,k) = tfion3
                   fion(1,m,j) = tfion4
                   fion(2,m,j) = tfion5
                   fion(3,m,j) = tfion6
                ENDIF
             ENDDO
1000         CONTINUE
          ENDDO
2000      CONTINUE
       ENDDO
       iat=iat+ions0%na(k)
    ENDDO
    ! 
    ! Embedded Atom Model
    ! 
    IF (tieam) THEN
       CALL eam_pot(esr,tau0,iesr,fion,tfor)
    ENDIF
    ! 
    CALL mp_sum(esr,parai%allgrp)
    IF (.NOT.paral%parent) esr=0._real_8
    CALL tihalt('   RPIIINT',isub)
    ! ==================================================================
    RETURN
  END SUBROUTINE rpiiint
  ! ==================================================================

END MODULE rpiiint_utils
