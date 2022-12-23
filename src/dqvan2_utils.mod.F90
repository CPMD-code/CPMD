MODULE dqvan2_utils
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

!!$private

!!$public :: dqvan2

!!$contains

END MODULE dqvan2_utils

! =================================================================
SUBROUTINE dqvan2(iv,jv,is,qrada,ylm,gagk,qg,dqg,spline)
  ! ==-------------------------------------------------------------==
  ! ==  DERIVATIVE VERSION OF QVAN2                                ==
  ! ==-------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE system , ONLY:lx,nbrx,ncpw,parm
  USE parac, ONLY : paral,parai
  USE nlps , ONLY:nlps_com
  USE cvan , ONLY:nelev
  USE geq0mod , ONLY:geq0
  USE qspl , ONLY:ggnh,nqdim,nsplpo,qspl1
  USE cppt , ONLY:gk,hg,igl,qrad,ylmb
  USE aavan , ONLY:ap,indv,lpl,lpx,nlx
  USE dylmr_utils, ONLY : dylmr, dbess, dylmrx, dylmrkx
  USE ylmr2_utils, ONLY : ylmr2
  USE fitpack_utils, ONLY : curv2, curvd
  USE zeroing_utils,                   ONLY: zeroing
  IMPLICIT NONE
  INTEGER                                    :: iv, jv, is
  REAL(real_8) :: qrada(nqdim,nbrx,nbrx,lx,*), ylm(*), gagk(ncpw%nhg,*), &
      qg(2,ncpw%nhg), dqg(2,ncpw%nhg,6), spline(*)

  COMPLEX(real_8)                            :: sig
  INTEGER                                    :: i, ig, is1, istep, isub, ivl, &
                                                ivs, jvl, jvs, kk, l, lp
  REAL(real_8)                               :: sigi, sigr, sigx, vol

! indeed complex and 2D
! Variables
! ==-------------------------------------------------------------==
! ==  IV  = 1..8    ! s_1 p_x1 p_z1 p_y1 s_2 p_x2 p_z2 p_y2      ==
! ==  IVS = 1..4    ! s_1 s_2 p_1 p_2                            ==
! ==  IVL = 1..4    ! s p_x p_z p_y                              ==
! ==                                                             ==
! == NOTE :   IV  = 1..8 (sppp sppp)                             ==
! ==    IVS = 1..4 (sspp) OR 1..2 (sp)                           ==
! ==    IVL = 1..4 (sppp)                                        ==
! ==-------------------------------------------------------------==

  CALL tiset('    DQVAN2',isub)
  ! ==-------------------------------------------------------------==
  vol=1._real_8/parm%omega
  ivs=indv(iv,is)
  jvs=indv(jv,is)
  istep=nlps_com%ngh(is)/nelev(is)
  ivl=1+MOD(iv-1,istep)
  jvl=1+MOD(jv-1,istep)
  IF (ivl.GT.nlx)  CALL stopgm(' DQVAN2 ',' IVL.GT.NLX  ',& 
       __LINE__,__FILE__)
  IF (jvl.GT.nlx)  CALL stopgm(' DQVAN2 ',' JVL.GT.NLX  ',& 
       __LINE__,__FILE__)
  IF (ivs.GT.nbrx) CALL stopgm(' DQVAN2 ',' IVS.GT.NBRX ',& 
       __LINE__,__FILE__)
  IF (jvs.GT.nbrx) CALL stopgm(' DQVAN2 ',' JVS.GT.NBRX ',& 
       __LINE__,__FILE__)
  CALL zeroing(qg)!,2*nhg)
  CALL zeroing(dqg)!,2*nhg*6)
  DO i=1,lpx(ivl,jvl)
     lp=lpl(ivl,jvl,i)
     ! EXTRACTION OF ANGULAR MOMENTUM L FROM LP:
     IF (lp.GE.26) CALL stopgm(' DQVAN2 ',' LP.GE.26 ',& 
          __LINE__,__FILE__)
     l=NINT(SQRT(REAL(lp,kind=real_8))+0.499999_real_8)
     ! 
     sig=(0._real_8,-1._real_8)**(l-1)
     CALL ylmr2(lp,ncpw%nhg,hg,gk,ylm)
     sigr=REAL(sig)
     sigi=AIMAG(sig)
     IF (qspl1%qspline) THEN
        !$omp parallel do private(IG) shared(SPLINE)
        DO ig=1,ncpw%nhg
           spline(ig)=curv2(hg(ig),nsplpo,ggnh(1),&
                qrad(1,ivs,jvs,l,is),&
                qrad(nsplpo+1,ivs,jvs,l,is),0.0_real_8)
        ENDDO
        IF (ABS(sigr).GT.0.5_real_8) THEN
           sigx=vol*sigr*ap(lp,ivl,jvl)
           !$omp parallel do private(IG) shared(QG)
           DO ig=1,ncpw%nhg
              qg(1,ig)=qg(1,ig)+sigx*ylmb(ig,lp,1)*spline(ig)
           ENDDO
           !$omp parallel do private(IG) shared(SPLINE)
           DO ig=1,ncpw%nhg
              spline(ig)=-2._real_8*hg(ig)*curvd(hg(ig),nsplpo,ggnh(1),&
                   qrad(1,ivs,jvs,l,is),&
                   qrad(nsplpo+1,ivs,jvs,l,is),0.0_real_8)
           ENDDO
           is1=1
           IF (geq0) is1=2
           !$omp parallel do private(KK,IG) shared(DQG)
           DO kk=1,6
              DO ig=is1,ncpw%nhg
                 dqg(1,ig,kk)=dqg(1,ig,kk)+sigx*ylmb(ig,lp,1)*&
                      spline(ig)*gagk(ig,kk)/(hg(ig)*parm%tpiba2)
              ENDDO
           ENDDO
        ELSE
           sigx=vol*sigi*ap(lp,ivl,jvl)
           !$omp parallel do private(IG) shared(QG)
           DO ig=1,ncpw%nhg
              qg(2,ig)=qg(2,ig)+sigx*ylmb(ig,lp,1)*spline(ig)
           ENDDO
           !$omp parallel do private(IG) shared(SPLINE)
           DO ig=1,ncpw%nhg
              spline(ig)=-2._real_8*hg(ig)*curvd(hg(ig),nsplpo,ggnh(1),&
                   qrad(1,ivs,jvs,l,is),&
                   qrad(nsplpo+1,ivs,jvs,l,is),0.0_real_8)
           ENDDO
           is1=1
           IF (geq0) is1=2
           !$omp parallel do private(KK,IG) shared(DQG)
           DO kk=1,6
              DO ig=is1,ncpw%nhg
                 dqg(2,ig,kk)=dqg(2,ig,kk)+sigx*ylmb(ig,lp,1)*&
                      spline(ig)*gagk(ig,kk)/(hg(ig)*parm%tpiba2)
              ENDDO
           ENDDO
        ENDIF
        !$omp parallel do private(IG) shared(SPLINE)
        DO ig=1,ncpw%nhg
           spline(ig)=curv2(hg(ig),nsplpo,ggnh(1),&
                qrad(1,ivs,jvs,l,is),&
                qrad(nsplpo+1,ivs,jvs,l,is),0.0_real_8)
        ENDDO
        IF (ABS(sigr).GT.0.5_real_8) THEN
           sigx=vol*sigr*ap(lp,ivl,jvl)
           DO kk=1,6
              CALL dylmr(lp,gk,ylm,kk)
              !$omp parallel do private(IG) shared(DQG)
              DO ig=1,ncpw%nhg
                 dqg(1,ig,kk)=dqg(1,ig,kk)+sigx*ylm(ig)*spline(ig)
              ENDDO
           ENDDO
        ELSE
           sigx=vol*sigi*ap(lp,ivl,jvl)
           DO kk=1,6
              CALL dylmr(lp,gk,ylm,kk)
              !$omp parallel do private(IG) shared(DQG)
              DO ig=1,ncpw%nhg
                 dqg(2,ig,kk)=dqg(2,ig,kk)+sigx*ylm(ig)*spline(ig)
              ENDDO
           ENDDO
        ENDIF
     ELSE
        IF (ABS(sigr).GT.0.5_real_8) THEN
           sigx=vol*sigr*ap(lp,ivl,jvl)
           !$omp parallel do private(IG) shared(QG)
           DO ig=1,ncpw%nhg
              qg(1,ig)=qg(1,ig)+sigx*ylmb(ig,lp,1)*&
                   qrad(igl(ig),ivs,jvs,l,is)
           ENDDO
           is1=1
           IF (geq0) is1=2
           !$omp parallel do private(KK,IG) shared(DQG)
           DO kk=1,6
              DO ig=is1,ncpw%nhg
                 dqg(1,ig,kk)=dqg(1,ig,kk)+sigx*ylmb(ig,lp,1)*&
                      qrada(igl(ig),ivs,jvs,l,is)*gagk(ig,kk)/&
                      (hg(ig)*parm%tpiba2)
              ENDDO
           ENDDO
        ELSE
           sigx=vol*sigi*ap(lp,ivl,jvl)
           !$omp parallel do private(IG) shared(QG)
           DO ig=1,ncpw%nhg
              qg(2,ig)=qg(2,ig)+sigx*ylmb(ig,lp,1)*&
                   qrad(igl(ig),ivs,jvs,l,is)
           ENDDO
           is1=1
           IF (geq0) is1=2
           !$omp parallel do private(KK,IG) shared(DQG)
           DO kk=1,6
              DO ig=is1,ncpw%nhg
                 dqg(2,ig,kk)=dqg(2,ig,kk)+sigx*ylmb(ig,lp,1)*&
                      qrada(igl(ig),ivs,jvs,l,is)*gagk(ig,kk)/&
                      (hg(ig)*parm%tpiba2)
              ENDDO
           ENDDO
        ENDIF
        IF (ABS(sigr).GT.0.5_real_8) THEN
           sigx=vol*sigr*ap(lp,ivl,jvl)
           DO kk=1,6
              CALL dylmr(lp,gk,ylm,kk)
              !$omp parallel do private(IG) shared(DQG)
              DO ig=1,ncpw%nhg
                 dqg(1,ig,kk)=dqg(1,ig,kk)+sigx*ylm(ig)*&
                      qrad(igl(ig),ivs,jvs,l,is)
              ENDDO
           ENDDO
        ELSE
           sigx=vol*sigi*ap(lp,ivl,jvl)
           DO kk=1,6
              CALL dylmr(lp,gk,ylm,kk)
              !$omp parallel do private(IG) shared(DQG)
              DO ig=1,ncpw%nhg
                 dqg(2,ig,kk)=dqg(2,ig,kk)+sigx*ylm(ig)*&
                      qrad(igl(ig),ivs,jvs,l,is)
              ENDDO
           ENDDO
        ENDIF
     ENDIF
  ENDDO
  CALL tihalt('    DQVAN2',isub)
  ! ==-------------------------------------------------------------==
  RETURN
END SUBROUTINE dqvan2
! =================================================================

