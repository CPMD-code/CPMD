#include "cpmd_global.h"

MODULE nlsm1_s_utils
  USE geq0mod,                         ONLY: geq0
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE kpnt,                            ONLY: eigkr
  USE kpts,                            ONLY: tkpts
  USE mp_interface,                    ONLY: mp_sum
  USE nlps,                            ONLY: imagp,&
                                             nghtol,&
                                             nlps_com
  USE parac,                           ONLY: parai
  USE str2,                            ONLY: becs,&
                                             betaa
  USE system,                          ONLY: maxsys,&
                                             ncpw,&
                                             nkpt,&
                                             parap
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: nlsm1_s

CONTAINS
  ! ==================================================================
  SUBROUTINE nlsm1_s(c0,scr,dai,nstate,ikind)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==  THE ARRAY BECS (THE DERIVATIVE OF FNL WITH RESPECT TO H)    ==
    ! ==  USE BETAA (DERIVATIVE OF TWNL) WHICH IS COMPUTED            ==
    ! ==            IN THE SUBROUTINE PUTBET                          ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: scr(nkpt%ngwk,ions1%nat)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: dai(imagp,maxsys%nax,nstate)
    COMPLEX(real_8)                          :: c0(nkpt%ngwk,nstate)
    INTEGER                                  :: ikind

    COMPLEX(real_8), PARAMETER               :: zone = (1._real_8,0._real_8) ,&
                                                zzero = (0._real_8,0._real_8)

    COMPLEX(real_8)                          :: ci
    INTEGER                                  :: i, ia, ig, ii, is, isa0, &
                                                isub, iv, kk
    REAL(real_8)                             :: cii, cir, ei, er, t

! ==--------------------------------------------------------------==

    CALL tiset('   NLSM1_S',isub)
    isa0=0
    DO is=1,ions1%nsp
       DO iv=1,nlps_com%ngh(is)
          ci=(0.0_real_8,-1.0_real_8)**nghtol(iv,is)
          cir=REAL(ci)
          cii=AIMAG(ci)
          DO kk=1,6
             ! MAKE USE OF THE SPECIAL STRUCTURE OF CI
             IF (ABS(cir).GT.0.5_real_8) THEN
                ! CI IS REAL
#ifdef __SR8000
                !poption parallel, tlocal(IA,IG,ER,EI,T)
#endif 
                !$OMP parallel do private(IA,IG,ER,EI,T) __COLLAPSE2
                DO ia=1,ions0%na(is)
                   !mb              isa=isa0+ia
                   DO ig=1,nkpt%ngwk
                      !mb                er=dreal(eigkr(ig,isa,ikind))
                      !mb                ei=dimag(eigkr(ig,isa,ikind))
                      er=REAL(eigkr(ig,isa0+ia,ikind),kind=real_8)
                      ei=AIMAG(eigkr(ig,isa0+ia,ikind))
                      t=betaa(ig,iv,kk,is,ikind)*cir
                      scr(ig,ia)=CMPLX(t*er,t*ei,kind=real_8)
                   ENDDO
                ENDDO
             ELSE
                ! CI IS IMAGINARY
#ifdef __SR8000
                !poption parallel, tlocal(IA,IG,ER,EI,T)
#endif 
                !$omp parallel do private(IA,IG,ER,EI,T) __COLLAPSE2
                DO ia=1,ions0%na(is)
                   !mb              isa=isa0+ia
                   DO ig=1,nkpt%ngwk
                      !mb                er=dreal(eigkr(ig,isa,ikind))
                      !mb                ei=dimag(eigkr(ig,isa,ikind))
                      er=REAL(eigkr(ig,isa0+ia,ikind),kind=real_8)
                      ei=AIMAG(eigkr(ig,isa0+ia,ikind))
                      t=betaa(ig,iv,kk,is,ikind)*cii
                      scr(ig,ia)=CMPLX(-t*ei,t*er,kind=real_8)
                   ENDDO
                ENDDO
             ENDIF
             IF (geq0) THEN
                IF (tkpts%tkpnt) THEN
                   !$omp parallel do private(IA)
                   DO ia=1,ions0%na(is)
                      scr(ncpw%ngw+1,ia)=CMPLX(0._real_8,0._real_8,kind=real_8)
                   ENDDO
                ELSE
                   !$omp parallel do private(IA)
                   DO ia=1,ions0%na(is)
                      scr(1,ia)=0.5_real_8*scr(1,ia)
                   ENDDO
                ENDIF
             ENDIF
             IF (tkpts%tkpnt) THEN
                CALL zgemm('C','N',ions0%na(is),nstate,nkpt%ngwk,zone,scr(1,1),&
                     nkpt%ngwk,c0(1,1),nkpt%ngwk,zzero,&
                     dai(1,1,1),maxsys%nax)
             ELSE
                CALL dgemm('T','N',ions0%na(is),nstate,2*ncpw%ngw,2.0_real_8,scr(1,1),&
                     2*ncpw%ngw,c0(1,1),2*ncpw%ngw,0.0_real_8,dai(1,1,1),maxsys%nax)
             ENDIF
             CALL mp_sum(dai,imagp*maxsys%nax*nstate,parai%allgrp)
#ifdef __VECTOR 
             IF (imagp.EQ.1) THEN
#ifdef __SR8000
                !poption parallel, tlocal(I,IA,II)
#endif 
                !$omp parallel do private(I,IA,II)
                DO i=parap%nst12(parai%mepos,1),parap%nst12(parai%mepos,2)
                   ii=i-parap%nst12(parai%mepos,1)+1
                   DO ia=1,ions0%na(is)
                      becs(1,isa0+ia,iv,kk,ii,ikind)=dai(1,ia,i)
                   ENDDO
                ENDDO
             ELSE
#ifdef __SR8000
                !poption parallel, tlocal(I,IA,II)
#endif 
                !$omp parallel do private(I,IA,II)
                DO i=parap%nst12(parai%mepos,1),parap%nst12(parai%mepos,2)
                   ii=i-parap%nst12(parai%mepos,1)+1
                   DO ia=1,ions0%na(is)
                      becs(1,isa0+ia,iv,kk,ii,ikind)=dai(1,ia,i)
                      becs(2,isa0+ia,iv,kk,ii,ikind)=dai(2,ia,i)
                   ENDDO
                ENDDO
             ENDIF
#else 
             DO i=parap%nst12(parai%mepos,1),parap%nst12(parai%mepos,2)
                ii=i-parap%nst12(parai%mepos,1)+1
                CALL dcopy(imagp*ions0%na(is),dai(1,1,i),1,&
                     becs(1,isa0+1,iv,kk,ii,ikind),1)
             ENDDO
#endif 
          ENDDO
       ENDDO
       isa0=isa0+ions0%na(is)
    ENDDO
    CALL tihalt('   NLSM1_S',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE nlsm1_s
  ! ==================================================================

END MODULE nlsm1_s_utils
