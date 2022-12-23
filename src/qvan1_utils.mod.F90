#ifdef __SR8000
!option MP(P(0)), LANGLVL(SAVE(0))
#endif

MODULE qvan1_utils
  USE aavan,                           ONLY: ap,&
                                             indv,&
                                             lpl,&
                                             lpx
  USE cnst,                            ONLY: fpi
  USE cvan,                            ONLY: nelev
  USE kinds,                           ONLY: real_8
  USE nlps,                            ONLY: nlps_com
  USE system,                          ONLY: nbrx

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: qvan1

CONTAINS

  ! =====================================================================
  SUBROUTINE qvan1(iv,jv,is,qrd0,qg0)
    ! ==-----------------------------------------------------------------==
    ! ==  Q(G,L,K) = SUM_LM (-I)^L AP(LM,L,K) YR_LM(G^) QRAD(G,L,L,K)    ==
    ! ==-----------------------------------------------------------------==
    INTEGER                                  :: iv, jv, is
    REAL(real_8)                             :: qrd0(nbrx,nbrx), qg0

    INTEGER                                  :: i, istep, IVl, ivs, jvl, jvs, &
                                                l, lp
    REAL(real_8)                             :: sigr, ylm0

    ivs=indv(iv,is)
    jvs=indv(jv,is)
    istep=nlps_com%ngh(is)/nelev(is)
    IVl=1+MOD(iv-1,istep)
    jvl=1+MOD(jv-1,istep)
    qg0=0.0_real_8
    ylm0=SQRT(1._real_8/fpi)
    DO i=1,lpx(IVl,jvl)
       lp=lpl(IVl,jvl,i)
       IF (lp.EQ.1) THEN
          l=1
          sigr=ap(lp,IVl,jvl)
          qg0=qg0+sigr*ylm0*qrd0(ivs,jvs)
       ENDIF
    ENDDO
    ! ==-----------------------------------------------------------------==
    RETURN
  END SUBROUTINE qvan1
  ! =====================================================================

END MODULE qvan1_utils
