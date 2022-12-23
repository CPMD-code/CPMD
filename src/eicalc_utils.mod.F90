#include "cpmd_global.h"

MODULE eicalc_utils
  USE cppt,                            ONLY: gk,&
                                             inyh,&
                                             rhops,&
                                             vps
  USE ions,                            ONLY: ions1
  USE kinds,                           ONLY: real_8
  USE nvtx_utils
  USE sfac,                            ONLY: ei1,&
                                             ei2,&
                                             ei3,&
                                             eigrb
  USE system,                          ONLY: cntl,&
                                             iatpt,&
                                             ncpw,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: eicalc
  PUBLIC :: eicalc1

CONTAINS

  ! ==================================================================
  SUBROUTINE eicalc(eivps,eirop)
    ! ==--------------------------------------------------------------==
    ! == EIVPS : phase factor times local pseudopotential  (VPS)      ==
    ! == EIROP : phase factor times Gaussian charge distributions     ==
    ! ==         which replaced ionic point charges (RHOPS)           ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: eivps(:), eirop(:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'eicalc'

    COMPLEX(real_8)                          :: ei123
    INTEGER                                  :: ia, ig, is, isa, isub
    REAL(real_8)                             :: ei, er

    CALL tiset(procedureN,isub)
    __NVTX_TIMER_START ( procedureN )

#if defined(__VECTOR) || defined(__SR11KIBM) || defined(__PRIMEHPC)
#if defined __NEC || defined __PRIMEHPC
    CALL zeroing(eivps)!,nhg)
    CALL zeroing(eirop)!,nhg)
    IF (cntl%bigmem) THEN
       !$omp parallel do private(IG,ISA,IS)
       DO ig=1,ncpw%nhg
          DO isa = 1, ions1%nat
             is=iatpt(2,isa)
             eivps(ig)=eivps(ig)+vps(is,ig)*eigrb(ig,isa)
             eirop(ig)=eirop(ig)+rhops(is,ig)*eigrb(ig,isa)
          ENDDO
       ENDDO
    ELSE
       !$omp parallel do private(IG,ISA,IS,EI123,ER,EI)
       DO ig=1,ncpw%nhg
          DO isa = 1, ions1%nat
             is=iatpt(2,isa)
             ei123=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
                  ei3(isa,inyh(3,ig))
             er=REAL(ei123)
             ei=AIMAG(ei123)
             eivps(ig)=eivps(ig)+CMPLX(er*vps(is,ig),ei*vps(is,ig),kind=real_8)
             eirop(ig)=eirop(ig)+&
                  CMPLX(er*rhops(is,ig),ei*rhops(is,ig),kind=real_8)
          ENDDO
       ENDDO
    ENDIF
#elif defined(_vpp_)
    CALL zeroing(eivps)!,nhg)
    CALL zeroing(eirop)!,nhg)
    IF (cntl%bigmem) THEN
       DO ig=1,ncpw%nhg
          DO isa=1,ions1%nat
             is=iatpt(2,isa)
             eivps(ig)=eivps(ig)+vps(is,ig)*eigrb(ig,isa)
             eirop(ig)=eirop(ig)+rhops(is,ig)*eigrb(ig,isa)
          ENDDO
       ENDDO
    ELSE
       DO ig=1,ncpw%nhg
          DO isa=1,ions1%nat
             is=iatpt(2,isa)
             ei123=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
                  ei3(isa,inyh(3,ig))
             eivps(ig)=eivps(ig)+vps(is,ig)*ei123
             eirop(ig)=eirop(ig)+rhops(is,ig)*ei123
          ENDDO
       ENDDO
    ENDIF
#else
    CALL zeroing(eivps)!,nhg)
    CALL zeroing(eirop)!,nhg)
    IF (cntl%bigmem) THEN
       !$omp parallel do private(IG,IS,IA,ISA)
#ifdef __SR8000
       !poption parallel
       !voption indep(RHOPS,EIGRB,VPS,EIROP,EIVPS)
#endif
       DO ig=1,ncpw%nhg
          DO isa=1,ions1%nat
             ia=iatpt(1,isa)
             is=iatpt(2,isa)
             eivps(ig)=eivps(ig)+vps(is,ig)*eigrb(ig,isa)
             eirop(ig)=eirop(ig)+rhops(is,ig)*eigrb(ig,isa)
          ENDDO
       ENDDO
    ELSE
       !$omp parallel do private(IG,IS,IA,ISA,EI123,ER,EI)
#ifdef __SR8000
       !poption parallel
       !voption indep(EI3,EI2,EI1,INYH,RHOPS,VPS,EIROP,EIVPS)
#endif
       DO ig=1,ncpw%nhg
          DO isa=1,ions1%nat
             ia=iatpt(1,isa)
             is=iatpt(2,isa)
             ei123=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
                  ei3(isa,inyh(3,ig))
             er=REAL(ei123)
             ei=AIMAG(ei123)
             eivps(ig)=eivps(ig)+CMPLX(er*vps(is,ig),ei*vps(is,ig),kind=real_8)
             eirop(ig)=eirop(ig)+&
                  CMPLX(er*rhops(is,ig),ei*rhops(is,ig),kind=real_8)
          ENDDO
       ENDDO
    ENDIF
#endif
#else
    CALL zeroing(eivps)!,nhg)
    CALL zeroing(eirop)!,nhg)
    IF (cntl%bigmem) THEN
       !$omp parallel do private(IG,IS,IA,ISA) shared(EIVPS,EIROP)
       DO ig=1,ncpw%nhg
          DO isa=1,ions1%nat
             ia=iatpt(1,isa)
             is=iatpt(2,isa)
             eivps(ig)=eivps(ig)+vps(is,ig)*eigrb(ig,isa)
             eirop(ig)=eirop(ig)+rhops(is,ig)*eigrb(ig,isa)
          ENDDO
       ENDDO
    ELSE
       !$omp parallel do private(IG,ISA,IA,IS,ER,EI,EI123) shared(EIVPS,EIROP)
       DO ig=1,ncpw%nhg
          DO isa=1,ions1%nat
             ia=iatpt(1,isa) ! acm: this is never used...
             is=iatpt(2,isa)
             ei123=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
                  ei3(isa,inyh(3,ig))
             er=REAL(ei123)
             ei=AIMAG(ei123)
             eivps(ig)=eivps(ig)+CMPLX(er*vps(is,ig),ei*vps(is,ig),kind=real_8)
             eirop(ig)=eirop(ig)+&
                  CMPLX(er*rhops(is,ig),ei*rhops(is,ig),kind=real_8)
          ENDDO
       ENDDO
    ENDIF
#endif

    __NVTX_TIMER_STOP
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE eicalc
  ! ==================================================================
  SUBROUTINE eicalc1(k,is,isa,eivps1,eirop1)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: k, is, isa
    COMPLEX(real_8)                          :: eivps1(ncpw%nhg), &
                                                eirop1(ncpw%nhg)

    COMPLEX(real_8)                          :: ei123, g
    INTEGER                                  :: ig, isub
    REAL(real_8)                             :: ei, er

    CALL tiset('   EICALC1',isub)
    IF (cntl%bigmem) THEN
       !$omp parallel do private(IG,G)
       DO ig=1,ncpw%nhg
          g=CMPLX(0._real_8,-gk(k,ig),kind=real_8)*parm%tpiba*eigrb(ig,isa)
          eivps1(ig)=vps(is,ig)*g
          eirop1(ig)=rhops(is,ig)*g
       ENDDO
    ELSE
       !$omp parallel do private(IG,G,EI123,ER,EI)
       DO ig=1,ncpw%nhg
          g=CMPLX(0._real_8,-gk(k,ig),kind=real_8)*parm%tpiba
          ei123=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
               ei3(isa,inyh(3,ig))
          er=REAL(ei123)
          ei=AIMAG(ei123)
          eivps1(ig)=CMPLX(er*vps(is,ig),ei*vps(is,ig),kind=real_8)*g
          eirop1(ig)=CMPLX(er*rhops(is,ig),ei*rhops(is,ig),kind=real_8)*g
       ENDDO
    ENDIF
    CALL tihalt('   EICALC1',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE eicalc1
  ! ==================================================================

END MODULE eicalc_utils
