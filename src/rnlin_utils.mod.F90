MODULE rnlin_utils
  USE atom,                            ONLY: atom_common,&
                                             gnl,&
                                             rps,&
                                             rv,&
                                             rw
  USE bessm_utils,                     ONLY: bessl
  USE cnst,                            ONLY: fpi,&
                                             pi
  USE dcacp_utils,                     ONLY: dcacp
  USE dpot,                            ONLY: dpot_mod
  USE error_handling,                  ONLY: stopgm
  USE fitpack_utils,                   ONLY: curv1
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE nlps,                            ONLY: nghtol,&
                                             nlps_com,&
                                             rgh
  USE parac,                           ONLY: parai
  USE pslo,                            ONLY: pslo_com
  USE qspl,                            ONLY: ggng,&
                                             nsplpa,&
                                             nsplpe,&
                                             nsplpo,&
                                             twns
  USE sgpp,                            ONLY: sgpp1,&
                                             sgpp2
  USE system,                          ONLY: parm
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: simpsn
  USE ylmr_utils,                      ONLY: bess
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: nlin

CONTAINS

  ! ==================================================================
  SUBROUTINE nlin(is,rs1,rs2)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: is
    REAL(real_8)                             :: rs1(nsplpo), rs2(nsplpo)

    INTEGER                                  :: isub

    CALL tiset('      NLIN',isub)
    IF (sgpp1%tsgp(is)) THEN
       IF (pslo_com%tnum(is)) THEN
          CALL unlin(is,rs1,rs2)
       ELSE
          CALL sgnlin(is,rs1,rs2)
       ENDIF
    ELSE
       CALL gnlin(is,rs1,rs2)
    ENDIF
    CALL tihalt('      NLIN',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE nlin
  ! ==================================================================
  SUBROUTINE sgnlin(is,rs1,rs2)
    ! ==--------------------------------------------------------------==
    ! ==   COMPUTES TWNS FOR NORMCONSERVING PSEUDO-POTENTIALS         ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: is
    REAL(real_8)                             :: rs1(nsplpo), rs2(nsplpo)

    INTEGER                                  :: ierr, il, iv, k, l
    REAL(real_8)                             :: gg, ggl, ggr, rc2

    DO iv=1,nlps_com%ngh(is)
       CALL zeroing(twns(:,1,iv,is))!,nsplpo)
       l=nghtol(iv,is)+1
       rc2=sgpp2%rcnl(l,is)**2
       !$omp parallel do private(IL,GG,GGL)
       DO il=nsplpa,nsplpe
          ggl=ggng(il)
          gg=ggl*parm%tpiba2
          rs1(il)=EXP(-0.5_real_8*gg*rc2)*pi**1.25_real_8
       ENDDO
       k=sgpp2%lfval(iv,is)
       IF (l.EQ.1) THEN
          IF (k.EQ.1) THEN
             !$omp parallel do private(IL,GG)
             DO il=nsplpa,nsplpe
                gg=ggng(il)*parm%tpiba2
                twns(il,1,iv,is)=4._real_8*sgpp2%rcnl(1,is)&
                     *SQRT(2._real_8*sgpp2%rcnl(1,is))*rs1(il)
             ENDDO
          ELSEIF (k.EQ.2) THEN
             !$omp parallel do private(IL,GG)
             DO il=nsplpa,nsplpe
                gg=ggng(il)*parm%tpiba2
                twns(il,1,iv,is)=8._real_8*sgpp2%rcnl(1,is)&
                     *SQRT(2._real_8/15._real_8*sgpp2%rcnl(1,is))*rs1(il)*(3._real_8-gg*rc2)
             ENDDO
          ELSEIF (k.EQ.3) THEN
             !$omp parallel do private(IL,GG,GGR)
             DO il=nsplpa,nsplpe
                gg=ggng(il)*parm%tpiba2
                ggr=gg*rc2*0.5_real_8
                twns(il,1,iv,is)=64._real_8*sgpp2%rcnl(l,is)&
                     *SQRT(2._real_8/945._real_8*sgpp2%rcnl(l,is))*&
                     rs1(il)*(3.75_real_8-5._real_8*ggr+ggr*ggr)
             ENDDO
          ELSE
             CALL stopgm('SGNLIN','K NOT PROGRAMMED',& 
                  __LINE__,__FILE__)
          ENDIF
       ELSEIF (l.EQ.2) THEN
          IF (k.EQ.1) THEN
             !$omp parallel do private(IL,GG)
             DO il=nsplpa,nsplpe
                gg=ggng(il)*parm%tpiba2
                twns(il,1,iv,is)=8._real_8*rc2*SQRT(sgpp2%rcnl(2,is)/3._real_8)*&
                     SQRT(gg)*rs1(il)
             ENDDO
          ELSEIF (k.EQ.2) THEN
             !$omp parallel do private(IL,GG)
             DO il=nsplpa,nsplpe
                gg=ggng(il)*parm%tpiba2
                twns(il,1,iv,is)=16._real_8*rc2*SQRT(sgpp2%rcnl(2,is)/105._real_8)*&
                     rs1(il)*SQRT(gg)*(5._real_8-gg*rc2)
             ENDDO
          ELSEIF (k.EQ.3) THEN
             !$omp parallel do private(IL,GG,GGR)
             DO il=nsplpa,nsplpe
                gg=ggng(il)*parm%tpiba2
                ggr=gg*rc2*0.5_real_8
                twns(il,1,iv,is)=128._real_8*rc2*SQRT(sgpp2%rcnl(l,is)/10395._real_8)*&
                     rs1(il)*SQRT(gg)*(8.75_real_8-7._real_8*ggr+ggr*ggr)
             ENDDO
          ELSE
             CALL stopgm('SGNLIN','K NOT PROGRAMMED',& 
                  __LINE__,__FILE__)
          ENDIF
       ELSEIF (l.EQ.3) THEN
          IF (k.EQ.1) THEN
             !$omp parallel do private(IL,GG)
             DO il=nsplpa,nsplpe
                gg=ggng(il)*parm%tpiba2
                twns(il,1,iv,is)=8._real_8*rc2*sgpp2%rcnl(l,is)*&
                     SQRT(2._real_8*sgpp2%rcnl(l,is)/15._real_8)*rs1(il)*gg
             ENDDO
          ELSEIF (k.EQ.2) THEN
             !$omp parallel do private(IL,GG,GGR)
             DO il=nsplpa,nsplpe
                gg=ggng(il)*parm%tpiba2
                ggr=gg*rc2
                twns(il,1,iv,is)=16._real_8*rc2*sgpp2%rcnl(l,is)*&
                     SQRT(2._real_8*sgpp2%rcnl(l,is)/945._real_8)*rs1(il)*gg*(7._real_8-ggr)
             ENDDO
          ELSEIF (k.EQ.3) THEN
             !$omp parallel do private(IL,GG,GGR)
             DO il=nsplpa,nsplpe
                gg=ggng(il)*parm%tpiba2
                ggr=gg*rc2*0.5_real_8
                twns(il,1,iv,is)=128._real_8*rc2&
                     *SQRT(2._real_8*sgpp2%rcnl(l,is)/135135._real_8)*&
                     sgpp2%rcnl(l,is)*rs1(il)*gg*(15.75_real_8-9._real_8*ggr+ggr*ggr)
             ENDDO
          ELSE
             CALL stopgm('SGNLIN','K NOT PROGRAMMED',& 
                  __LINE__,__FILE__)
          ENDIF
       ELSEIF (l.EQ.4) THEN
          IF (k.EQ.1) THEN
             !$omp parallel do private(IL,GG)
             DO il=nsplpa,nsplpe
                gg=ggng(il)*parm%tpiba2
                twns(il,1,iv,is)=16._real_8*rc2*rc2*gg*SQRT(gg)*rs1(il)&
                     *SQRT(sgpp2%rcnl(l,is)/105._real_8)
             ENDDO
          ELSEIF (k.EQ.2) THEN
             !$omp parallel do private(IL,GG,GGR)
             DO il=nsplpa,nsplpe
                gg=ggng(il)*parm%tpiba2
                ggr=gg*rc2
                twns(il,1,iv,is)=32._real_8*rc2*rc2*gg*SQRT(gg)*rs1(il)&
                     *SQRT(sgpp2%rcnl(l,is)/10395._real_8)*(9._real_8-ggr)
             ENDDO
          ELSEIF (k.EQ.3) THEN
             !$omp parallel do private(IL,GG,GGR)
             DO il=nsplpa,nsplpe
                gg=ggng(il)*parm%tpiba2
                ggr=gg*rc2*0.5_real_8
                twns(il,1,iv,is)=256._real_8*rc2*rc2*gg*SQRT(gg)*rs1(il)&
                     *SQRT(sgpp2%rcnl(l,is)/2027025._real_8)&
                     *(24.75_real_8-11._real_8*ggr+ggr*ggr)
             ENDDO
          ELSE
             CALL stopgm('SGNLIN','K NOT PROGRAMMED',& 
                  __LINE__,__FILE__)
          ENDIF
       ENDIF
       CALL mp_sum(twns(:,:,iv,is),nsplpo,parai%allgrp)
       CALL curv1(nsplpo,ggng,twns(1,1,iv,is),&
            0.0_real_8,0.0_real_8,3,twns(1,2,iv,is),rs2,0.0_real_8,ierr)
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE sgnlin
  ! ==================================================================
  SUBROUTINE gnlin(is,jl,fint)
    ! ==--------------------------------------------------------------==
    ! ==   COMPUTES TWNS FOR NORMCONSERVING PSEUDO-POTENTIALS         ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: is
    REAL(real_8)                             :: jl(atom_common%meshw(is)), &
                                                fint(atom_common%meshw(is))

    INTEGER                                  :: ierr, il, ir, iv, l, lold, m, &
                                                mmax, nghl
    REAL(real_8)                             :: c, ggl, tmp, xg, xgr

    IF (dpot_mod%tkb(is)) THEN
       ! Kleinman-Bylander scheme
       lold=-1
       mmax=atom_common%meshw(is)
       c=atom_common%clogw(is)*fpi
       DO iv=1,nlps_com%ngh(is)
          l=nghtol(iv,is)+1
          IF (l.EQ.lold) THEN
             CALL dcopy(2*nsplpo,twns(1,1,iv-1,is),1,twns(1,1,iv,is),1)
          ELSE
             CALL zeroing(twns(:,1,iv,is))!,nsplpo)
             IF ( l /= dpot_mod%lskip(is) .and. .not. dcacp%skip(l,is) ) THEN
                DO il=nsplpa,nsplpe
                   ggl=ggng(il)
                   xg=SQRT(ggl)*parm%tpiba
                   CALL bess(xg,l,mmax,rw(1,is),jl)
                   DO ir=1,mmax
                      fint(ir)=rw(ir,is)**2*rps(ir,is,l)*gnl(ir,is,l)*jl(ir)
                   ENDDO
                   CALL simpsn(mmax,fint,tmp)
                   twns(il,1,iv,is)=c*tmp
                ENDDO
                CALL mp_sum(twns(:,:,iv,is),nsplpo,parai%allgrp)
                CALL curv1(nsplpo,ggng,twns(1,1,iv,is),&
                     0.0_real_8,0.0_real_8,3,twns(1,2,iv,is),fint,0.0_real_8,ierr)
             ENDIF
          ENDIF
          lold=l
       ENDDO
    ELSE
       ! Gauss-Hermit integration;
       ! Initialize spherical Bessel functions
       nghl=NINT(nlps_com%rmaxn(is))
       c=fpi
       DO iv=1,nlps_com%ngh(is)
          CALL zeroing(twns(:,1,iv,is))!,nsplpo)
          m=MOD(iv-1,nghl)+1
          l=nghtol(iv,is)+1
          IF ( l /= dpot_mod%lskip(is) ) THEN
             !$omp parallel do private(IL,GGL,XGR)
             DO il=nsplpa,nsplpe
                ggl=ggng(il)
                xgr=SQRT(ggl)*parm%tpiba*rgh(m,is)
                twns(il,1,iv,is)=c*bessl(l-1,xgr)
             ENDDO
             CALL mp_sum(twns(:,:,iv,is),nsplpo,parai%allgrp)
             CALL curv1(nsplpo,ggng,twns(1,1,iv,is),&
                  0.0_real_8,0.0_real_8,3,twns(1,2,iv,is),fint,0.0_real_8,ierr)
          ENDIF
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gnlin
  ! ==================================================================
  SUBROUTINE unlin(is,jl,fint)
    ! ==--------------------------------------------------------------==
    ! ==   COMPUTES TWNS FOR NORMCONSERVING PSEUDO-POTENTIALS         ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: is
    REAL(real_8)                             :: jl(atom_common%meshvp(is)), &
                                                fint(atom_common%meshvp(is))

    INTEGER                                  :: ierr, il, ir, iv, l, lold, &
                                                mmax
    REAL(real_8)                             :: c, ggl, tmp, xg

    lold=-1
    mmax=atom_common%meshvp(is)
    c=fpi
    DO iv=1,nlps_com%ngh(is)
       l=nghtol(iv,is)+1
       IF (l.EQ.lold) THEN
          CALL dcopy(2*nsplpo,twns(1,1,iv-1,is),1,twns(1,1,iv,is),1)
       ELSE
          CALL zeroing(twns(:,1,iv,is))!,nsplpo)
          DO il=nsplpa,nsplpe
             ggl=ggng(il)
             xg=SQRT(ggl)*parm%tpiba
             CALL bess(xg,l,mmax,rv(1,is),jl)
             DO ir=1,mmax
                fint(ir)=rw(ir,is)*rv(ir,is)*gnl(ir,is,l)*jl(ir)
             ENDDO
             CALL simpsn(mmax,fint,tmp)
             twns(il,1,iv,is)=c*tmp
          ENDDO
          CALL mp_sum(twns(:,:,iv,is),nsplpo,parai%allgrp)
          CALL curv1(nsplpo,ggng,twns(1,1,iv,is),&
               0.0_real_8,0.0_real_8,3,twns(1,2,iv,is),fint,0.0_real_8,ierr)
       ENDIF
       lold=l
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE unlin
  ! ==================================================================

END MODULE rnlin_utils
