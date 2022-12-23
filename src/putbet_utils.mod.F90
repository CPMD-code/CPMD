MODULE putbet_utils
  USE cppt,                            ONLY: gk,&
                                             hg
  USE dpot,                            ONLY: dpot_mod
  USE dylmr_utils,                     ONLY: dylmrkx,&
                                             dylmrx
  USE fitpack_utils,                   ONLY: curv2,&
                                             curvd
  USE geq0mod,                         ONLY: geq0
  USE ions,                            ONLY: ions1
  USE kinds,                           ONLY: real_8
  USE kpnt,                            ONLY: hgkm,&
                                             hgkp,&
                                             rk
  USE kpts,                            ONLY: tkpts
  USE nlps,                            ONLY: nghcom,&
                                             nlps_com
  USE pslo,                            ONLY: pslo_com
  USE qspl,                            ONLY: ggng,&
                                             nsplpo,&
                                             twns
  USE sgpp,                            ONLY: sgpp1,&
                                             sgpp2
  USE str2,                            ONLY: betaa,&
                                             gagk,&
                                             gagkm,&
                                             gagkp
  USE system,                          ONLY: kpbeg,&
                                             ncpw,&
                                             nkpbl,&
                                             parm
  USE vdbp,                            ONLY: ncpr1
  USE ylmr_utils,                      ONLY: ylmkr,&
                                             ylmr

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: putbet

CONTAINS

  ! ==================================================================
  SUBROUTINE putbet(ikpt)
    ! ==--------------------------------------------------------------==
    ! == COMPUTES THE ARRAY BETAA WHICH IS USED IN NLSM1_S ROUTINE    ==
    ! == BETAA IS THE DERIVATIVE OF TWNL VERSUS STRAIN                ==
    ! == R -> (1 + e) R                                               ==
    ! == G -> (1 - e) G                                               ==
    ! ==--------------------------------------------------------------==
    ! == IKPT is for block option the block index                     ==
    ! ==         otherwise is not useful                              ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ikpt

    INTEGER                                  :: ig, ig1, ikind, ikk, is, &
                                                istep, iv, kk, lp, nkpoint
    REAL(real_8)                             :: aux1, aux2, dw, tw, vol, &
                                                voltpi

    nkpoint=nkpbl(ikpt)
    DO ikind=1,nkpoint
       ikk=kpbeg(ikpt)+ikind
       vol=1._real_8/SQRT(parm%omega)
       voltpi=vol/parm%tpiba2
       DO is=1,ions1%nsp
          DO iv=1,nlps_com%ngh(is)
             IF (pslo_com%tvan(is)) THEN
                istep=ncpr1%nvales(is)*ncpr1%nvales(is)
                lp=1+MOD(iv-1,istep)
             ELSEIF (dpot_mod%tkb(is)) THEN
                lp=nghcom(iv,is)
             ELSEIF (sgpp1%tsgp(is)) THEN
                lp=sgpp2%lpval(iv,is)
             ELSE
                lp = (iv-1)/NINT(nlps_com%rmaxn(is)) + 1
             ENDIF
             ! 0-component
             IF (geq0) THEN
                IF (tkpts%tkpnt) THEN
                   IF (  (rk(1,ikind).EQ.0._real_8).AND.&
                        (rk(2,ikind).EQ.0._real_8).AND.&
                        (rk(3,ikind).EQ.0._real_8) ) THEN
                      ig1=2
                      tw=twns(1,1,iv,is)
                      DO kk=1,6
                         betaa(    1,iv,kk,is,ikind)=dylmrkx(lp,1,gk(1,1),&
                              rk(1,ikk),gagkp(1,1,ikind),kk, 1)*tw*vol
                         betaa(ncpw%ngw+1,iv,kk,is,ikind)=dylmrkx(lp,1,gk(1,1),&
                              rk(1,ikk),gagkm(1,1,ikind),kk,-1)*tw*vol
                      ENDDO
                   ELSE
                      ig1=1
                   ENDIF
                ELSE
                   ig1=2
                   tw=twns(1,1,iv,is)
                   DO kk=1,6
                      betaa(1,iv,kk,is,ikind)=dylmrx(lp,1,gk(1,1),kk)*tw*vol
                   ENDDO
                ENDIF
             ELSE
                ig1=1
             ENDIF
             ! K points version
             IF (tkpts%tkpnt) THEN
                DO ig=ig1,ncpw%ngw
                   ! RK+GK (GAGKP)
                   tw=curv2(hgkp(ig,ikind),nsplpo,ggng(1),&
                        twns(1,1,iv,is),twns(1,2,iv,is),0.0_real_8)*vol
                   dw=-2._real_8*curvd(hgkp(ig,ikind),nsplpo,ggng(1),&
                        twns(1,1,iv,is),twns(1,2,iv,is),0.0_real_8)&
                        *voltpi
                   DO kk=1,6
                      aux1=dylmrkx(lp,ig,gk(1,1),rk(1,ikk),gagkp(1,1,ikind),&
                           kk, 1)*tw
                      aux2=ylmkr(lp,ig,gk(1,1),rk(1,ikk), 1)*&
                           dw*gagkp(ig,kk,ikind)
                      betaa(ig,iv,kk,is,ikind)=aux1+aux2
                   ENDDO
                   ! RK-GK (GAGKM)
                   tw=curv2(hgkm(ig,ikind),nsplpo,ggng(1),&
                        twns(1,1,iv,is),twns(1,2,iv,is),0.0_real_8)*vol
                   dw=-2._real_8*curvd(hgkm(ig,ikind),nsplpo,ggng(1),&
                        twns(1,1,iv,is),twns(1,2,iv,is),0.0_real_8)&
                        *voltpi
                   DO kk=1,6
                      aux1=dylmrkx(lp,ig,gk(1,1),rk(1,ikk),gagkm(1,1,ikind),&
                           kk,-1)*tw
                      aux2=ylmkr(lp,ig,gk(1,1),rk(1,ikk),-1)*&
                           dw*gagkm(ig,kk,ikind)
                      betaa(ncpw%ngw+ig,iv,kk,is,ikind)=aux1+aux2
                   ENDDO
                ENDDO
             ELSE
                ! Gamma point
#ifdef __SR8000
                !poption parallel, tlocal(IG,KK,TW,DW,AUX1,AUX2)
#endif 
                !$omp parallel do private(IG,KK,TW,DW,AUX1,AUX2)
                DO ig=ig1,ncpw%ngw
                   tw=curv2(hg(ig),nsplpo,ggng(1),twns(1,1,iv,is),&
                        twns(1,2,iv,is),0.0_real_8)*vol
                   dw=-2._real_8*curvd(hg(ig),nsplpo,ggng(1),twns(1,1,iv,is),&
                        twns(1,2,iv,is),0.0_real_8)&
                        *voltpi
                   DO kk=1,6
                      aux1=dylmrx(lp,ig,gk(1,1),kk)*tw
                      aux2=ylmr(lp,ig,gk(1,1))*dw*gagk(ig,kk)
                      betaa(ig,iv,kk,is,ikind)=aux1+aux2
                   ENDDO
                ENDDO
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE putbet
  ! ==================================================================

END MODULE putbet_utils
