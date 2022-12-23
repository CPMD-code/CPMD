MODULE vlocst_utils
  USE cppt,                            ONLY: nzh
  USE ffsum_utils,                     ONLY: ffsum
  USE fft_maxfft,                      ONLY: maxfft
  USE geq0mod,                         ONLY: geq0
  USE kinds,                           ONLY: real_8
  USE pslo,                            ONLY: pslo_com
  USE str2,                            ONLY: dqg,&
                                             drhovg,&
                                             gagk
  USE strs,                            ONLY: alpha,&
                                             beta,&
                                             delta,&
                                             deps,&
                                             depst
  USE system,                          ONLY: ncpw,&
                                             parm,&
                                             cntl
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: vlocst

CONTAINS

  ! ==================================================================
  SUBROUTINE vlocst(epseu,v,eivps)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES DEPS(1:6)                    ==
    ! == HARTREE ENERGY CONTRIBUTION TO THE STRESS TENSOR             ==
    ! == Nielsen and Martin, PRB 32,p 3793 (1985)                     ==
    ! ==--------------------------------------------------------------==
    ! == V(1:MAXFFT) [in] electronic density in reciprocal space      ==
    ! == EIVPS(1:NHG): [in] phase factor times local pseudopotential  ==
    ! ==             used for Vanderbilt pp                           ==
    ! == EPSEU: Contribution of local pseudopotential part            ==
    ! == Use:                                                         ==
    ! ==      GAGK(1:NHG,1:6) GK_alpha * GK_beta                      ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: epseu
    COMPLEX(real_8)                          :: v(maxfft), eivps(ncpw%nhg)

    INTEGER                                  :: ig, ig1, isub, kk

#ifdef __VECTOR 
    COMPLEX(real_8)                          :: temp,t1,t2,t3,t4,t5,t6
#endif 
    ! ==--------------------------------------------------------------==
    CALL tiset('    VLOCST',isub)
    CALL zeroing(dqg)!,nhg)
    ! Calculate derivative local pp part*structure factor
    CALL ffsum(dqg)
    ! Sum excluding G=0
    IF (geq0) dqg(1,1)=(0.0_real_8,0.0_real_8)
#if defined(__VECTOR)
    t1=(0.0_real_8,0.0_real_8)
    t2=(0.0_real_8,0.0_real_8)
    t3=(0.0_real_8,0.0_real_8)
    t4=(0.0_real_8,0.0_real_8)
    t5=(0.0_real_8,0.0_real_8)
    t6=(0.0_real_8,0.0_real_8)
#ifdef __SR8000
    !poption parallel, tlocal(IG) 
    !poption psum(T1,T2,T3,T4,T5,T6)
#endif 
    !$omp parallel do private(IG) &
    !$omp  reduction(+:T1,T2,T3,T4,T5,T6)
    DO ig=1,ncpw%nhg
       t1=t1+CONJG(v(nzh(ig)))*dqg(ig,1)*gagk(ig,1)
       t2=t2+CONJG(v(nzh(ig)))*dqg(ig,1)*gagk(ig,2)
       t3=t3+CONJG(v(nzh(ig)))*dqg(ig,1)*gagk(ig,3)
       t4=t4+CONJG(v(nzh(ig)))*dqg(ig,1)*gagk(ig,4)
       t5=t5+CONJG(v(nzh(ig)))*dqg(ig,1)*gagk(ig,5)
       t6=t6+CONJG(v(nzh(ig)))*dqg(ig,1)*gagk(ig,6)
    ENDDO
    depst(1)=depst(1)+t1
    depst(2)=depst(2)+t2
    depst(3)=depst(3)+t3
    depst(4)=depst(4)+t4
    depst(5)=depst(5)+t5
    depst(6)=depst(6)+t6
#else 
    !$omp parallel do private(KK,IG) reduction(+:depst)
    DO ig=1,ncpw%nhg
       DO kk=1,6
          depst(kk)=depst(kk)+CONJG(v(nzh(ig)))*dqg(ig,1)*gagk(ig,kk)
       ENDDO
    ENDDO
#endif 
    IF (pslo_com%tivan) THEN
       DO kk=1,6
          ig1=1
          IF (geq0) ig1=2
          IF (delta(alpha(kk),beta(kk)).GT.0.5_real_8) THEN
#ifdef __VECTOR 
             temp=(0.0_real_8,0.0_real_8)
#ifdef __SR8000
             !poption parallel, tlocal(IG), psum(TEMP)
#endif 
             !$omp parallel do private(IG) reduction(+:TEMP)
             DO ig=ig1,ncpw%nhg
                temp=temp+CONJG(drhovg(ig,kk))*eivps(ig)
                IF (cntl%tlsd) temp=temp+CONJG(drhovg(ig,6+kk))*eivps(ig)
             ENDDO
             depst(kk)=depst(kk)+temp
#else 
             !$omp parallel do private(ig) reduction(+:depst)
             DO ig=ig1,ncpw%nhg
                depst(kk)=depst(kk)+CONJG(drhovg(ig,kk))*eivps(ig)
                IF (cntl%tlsd) depst(kk)=depst(kk)+CONJG(drhovg(ig,6+kk))*eivps(ig)
             ENDDO
#endif 
          ENDIF
       ENDDO
    ENDIF
    ! Factor 2._real_8 with GAGK
#ifdef __VECTOR 
    depst(1)=2.0_real_8*depst(1)
    depst(2)=2.0_real_8*depst(2)
    depst(3)=2.0_real_8*depst(3)
    depst(4)=2.0_real_8*depst(4)
    depst(5)=2.0_real_8*depst(5)
    depst(6)=2.0_real_8*depst(6)
#else 
    DO kk=1,6
       depst(kk)=2.0_real_8*depst(kk)
    ENDDO
#endif 
    ! Add contribution of local part (average non-Coulomb part) 
#ifdef __VECTOR 
    deps(1)=-epseu*delta(alpha(1),beta(1))+parm%omega*REAL(depst(1))
    deps(2)=-epseu*delta(alpha(2),beta(2))+parm%omega*REAL(depst(2))
    deps(3)=-epseu*delta(alpha(3),beta(3))+parm%omega*REAL(depst(3))
    deps(4)=-epseu*delta(alpha(4),beta(4))+parm%omega*REAL(depst(4))
    deps(5)=-epseu*delta(alpha(5),beta(5))+parm%omega*REAL(depst(5))
    deps(6)=-epseu*delta(alpha(6),beta(6))+parm%omega*REAL(depst(6))
#else 
    DO kk=1,6
       deps(kk) = -epseu*delta(alpha(kk),beta(kk)) +&
            parm%omega*REAL(depst(kk))
    ENDDO
#endif 
    CALL tihalt('    VLOCST',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE vlocst
  ! ==================================================================

END MODULE vlocst_utils
