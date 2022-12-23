MODULE hesele_utils
  USE cppt,                            ONLY: hg,&
                                             twnl
  USE cvan,                            ONLY: deeq,&
                                             dvan
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_bcast
  USE nlps,                            ONLY: nghtol,&
                                             nlps_com,&
                                             wsg
  USE parac,                           ONLY: parai
  USE pslo,                            ONLY: pslo_com
  USE sgpp,                            ONLY: sgpp1,&
                                             sgpp2
  USE simulmod,                        ONLY: vploc
  USE system,                          ONLY: cntl,&
                                             cntr,&
                                             ncpw,&
                                             parm

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: hesele

CONTAINS

  ! ==================================================================
  SUBROUTINE hesele (svar2,vpp)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: svar2, vpp(ncpw%ngw)

    INTEGER                                  :: ia, ig, is, isa0, iv, j, jv, &
                                                ki, kj, l, l2, li, lj
    REAL(real_8)                             :: dd, emp, ftpi, vp

! Variables
! ==--------------------------------------------------------------==
! ==  Calculate   H0 : The diagonal approximation to the          ==
! ==                   2nd derivative matrix                      ==
! ==--------------------------------------------------------------==

    IF (cntl%prec) THEN
       ! ==--------------------------------------------------------------==
       ! ==  Local Potential                                             ==
       ! ==--------------------------------------------------------------==
       CALL mp_bcast(vploc,parai%igeq0,parai%allgrp)
       vp=vploc
       ! ==--------------------------------------------------------------==
       ! ==  Nonlocal Potential                                          ==
       ! ==--------------------------------------------------------------==
#if ! defined(__VECTOR)
       !$omp parallel do private(IG,IS,IV,JV,IA,ISA0,DD,L,KI,LI,L2,LJ,KJ) 
#else
       !$omp parallel do private(IG,IS,IV,JV,IA,ISA0,DD,L,KI,LI,L2,LJ,KJ) &
       !$omp  schedule(static)
#endif
       DO ig=1,ncpw%ngw
          vpp(ig)=vp
          isa0=0
          DO is=1,ions1%nsp
             IF (pslo_com%tvan(is)) THEN
                DO iv=1,nlps_com%ngh(is)
                   DO jv=1,nlps_com%ngh(is)
                      dd=dvan(iv,jv,is)
                      DO ia=1,ions0%na(is)
                         dd=dd+deeq(isa0+ia,iv,jv,1)
                      ENDDO
                      IF (cntl%tlsd) THEN
                         DO ia=1,ions0%na(is)
                            dd=dd+deeq(isa0+ia,iv,jv,2)
                         ENDDO
                         dd=0.5_real_8*dd
                      ENDIF
                      vpp(ig)=vpp(ig)+dd*twnl(ig,iv,is,1)*twnl(ig,jv,is,1)
                   ENDDO
                ENDDO
             ELSEIF (sgpp1%tsgp(is)) THEN
                DO iv=1,nlps_com%ngh(is)
                   l=nghtol(iv,is)+1
                   ki=sgpp2%lfval(iv,is)
                   li=sgpp2%lpval(iv,is)
                   DO jv=1,nlps_com%ngh(is)
                      l2=nghtol(jv,is)+1
                      lj=sgpp2%lpval(jv,is)
                      IF (l.EQ.l2.AND.li.EQ.lj) THEN
                         kj=sgpp2%lfval(jv,is)
                         dd=sgpp2%hlsg(ki,kj,l,is)*ions0%na(is)
                         vpp(ig)=vpp(ig)+dd*twnl(ig,iv,is,1)*twnl(ig,jv,is,1)
                      ENDIF
                   ENDDO
                ENDDO
             ELSE
                DO iv=1,nlps_com%ngh(is)
                   vpp(ig)=vpp(ig)+ions0%na(is)*wsg(is,iv)&
                        *twnl(ig,iv,is,1)*twnl(ig,iv,is,1)
                ENDDO
             ENDIF
             isa0=isa0+ions0%na(is)
          ENDDO
       ENDDO
       ! ==--------------------------------------------------------------==
       ! ==  Kinetic Energy                                              ==
       ! ==--------------------------------------------------------------==
       ftpi=0.5_real_8*parm%tpiba2
#if ! defined(__VECTOR)
       !$omp parallel do private(J,EMP) 
#else
       !$omp parallel do private(J,EMP) &
       !$omp  schedule(static)
#endif
       DO j=1,ncpw%ngw
          vpp(j)=ftpi*hg(j)+vpp(j)
          emp=cntr%hthrs
          IF (vpp(j).LT.emp) vpp(j)=emp
       ENDDO
#if ! defined(__VECTOR)
       !$omp parallel do private(J) 
#else
       !$omp parallel do private(J) &
       !$omp  schedule(static)
#endif
       DO j=1,ncpw%ngw
          vpp(j)=1.0_real_8/vpp(j)
       ENDDO
    ELSE
#if ! defined(__VECTOR)
       !$omp parallel do private(J) 
#else
       !$omp parallel do private(J) &
       !$omp  schedule(static)
#endif
       DO j=1,ncpw%ngw
          vpp(j)=svar2
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE hesele
  ! ==================================================================

END MODULE hesele_utils
