MODULE vdbinit_utils
  USE aavan,                           ONLY: indv
  USE cnst,                            ONLY: fpi
  USE cppt,                            ONLY: gl,&
                                             qrad
  USE cvan,                            ONLY: dvan,&
                                             nelev,&
                                             qq
  USE error_handling,                  ONLY: stopgm
  USE fitpack_utils,                   ONLY: curv1,&
                                             curv2
  USE ions,                            ONLY: ions1
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE nlps,                            ONLY: nghtol,&
                                             nlps_com
  USE parac,                           ONLY: parai
  USE pslo,                            ONLY: pslo_com
  USE qspl,                            ONLY: ggng,&
                                             ggnh,&
                                             nqdim,&
                                             nsplpa,&
                                             nsplpe,&
                                             nsplpo,&
                                             qspl1,&
                                             twns
  USE qvan1_utils,                     ONLY: qvan1
  USE radin_utils,                     ONLY: radin,&
                                             radlg
  USE system,                          ONLY: maxsys,&
                                             nbrx,&
                                             ncpw,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE vdbp,                            ONLY: betar,&
                                             dion,&
                                             ncpr1,&
                                             qfunc,&
                                             qrl,&
                                             r,&
                                             rab
  USE ylmr_utils,                      ONLY: bess
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: vdbinit
  PUBLIC :: qinit

CONTAINS

  ! ==================================================================
  SUBROUTINE vdbinit(is,jl,fint)
    ! ==--------------------------------------------------------------==
    ! ==  THIS ROUTINE INITALIZES ARRAYS   BETA  QRAD  QQ             ==
    ! ==  THIS ROUTINE IS MOSTLY FROM THE ROUTINE NLINIT IN THE       ==
    ! ==  ORIGINAL VANDERBILT CODE                                    ==
    ! ==                                                              ==
    ! ==  BETA(IG,L,IS) = 4PI/SQRT(OMEGA) Y^R(L,Q^)                   ==
    ! ==                        INT_0^INF DR R^2 J_L(QR) BETA(L,IS,R) ==
    ! ==                                                              ==
    ! ==  QRAD(IGL,L,K,IS) = 4PI/OMEGA                                ==
    ! ==                   INT_0^R DR R^2 J_L(QR) Q(R,L,K,IS)         ==
    ! ==                                                              ==
    ! ==  BETA(G)_LM,IS = (-I)^L*BETA(IG,L,IS)                        ==
    ! ==                                                              ==
    ! ==  QQ_IJ=INT_0^R Q_IJ(R)=OMEGA*QG(G=0)                         ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: is
    REAL(real_8)                             :: jl(*), fint(*)

    INTEGER                                  :: i, ierr, il, iold, ir, istep, &
                                                isub, iv, jv, l, lm, lval, m, &
                                                mmax, nt
    REAL(real_8)                             :: betax, ggl, xg

! Variables
! ==--------------------------------------------------------------==
! ==  qq order:  s_1 p_x1 p_z1 p_y1 ... s_2 p_x2 p_z2 p_y2  ...   ==
! ==                               spppddddd spppddddd ...        ==
! == qqq order: s_1 s_2 p_1 p_2 d_1 d_2              ss..pp..dd.. ==
! ==                                                              ==
! ==  NELEV  =  NO. OF ENERGY LEVELS                              ==
! ==  NGH(IS) = TOTAL NUMBER OF STATES                            ==
! ==  ISTEP  =  LENGTH OF SPPPDDD.. SEQUENCE = 1+3+5+...          ==
! ==--------------------------------------------------------------==

    CALL tiset('   VDBINIT',isub)
    lval=ncpr1%nvales(is)
    nt=ncpr1%nvales(is)
    nelev(is)=ncpr1%nbeta(is)/ncpr1%nvales(is)
    nlps_com%ngh(is)=nt*nt*nelev(is)
    mmax=ncpr1%kkbeta(is)
    istep=nt*nt
    IF (nlps_com%ngh(is).GT.maxsys%nhxs) CALL stopgm(' NLINIT',' NGH.GT.NHX ',& 
         __LINE__,__FILE__)
    DO i=1,nlps_com%ngh(is),istep
       lm=-1
       DO l=0,ncpr1%nvales(is)-1
          DO m=0,2*l
             lm=lm+1
             iv=i+lm
             nghtol(iv,is)=l
             indv(iv,is)=l*nelev(is)+1+(i-1)/istep
          ENDDO
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    ! ==   INDV(1)=1         ! QQ ORDER:                              ==
    ! ==   INDV(2)=3         ! s_1 p_x1 p_z1 p_y1 s_2 p_x2 p_z2 p_y2  ==
    ! ==   INDV(3)=3         ! QQQ ORDER:                             ==
    ! ==   INDV(4)=3         ! s_1 s_2 p_1 p_2                        ==
    ! ==   INDV(5)=2         ! SEE ALSO YLMR (OR AAINIT)              ==
    ! ==   INDV(6)=4                                                  ==
    ! ==   INDV(7)=4                                                  ==
    ! ==   INDV(8)=4                                                  ==
    ! ==--------------------------------------------------------------==
    ! ==  CALCULATION OF ARRAY  TWNL(IG,IV,IS)                        ==
    ! ==--------------------------------------------------------------==
    iold=-1
    DO iv=1,nlps_com%ngh(is)
       IF (iold.EQ.indv(iv,is)) THEN
          CALL dcopy(2*nsplpo,twns(1,1,iv-1,is),1,twns(1,1,iv,is),1)
       ELSE
          l=nghtol(iv,is)+1
          CALL zeroing(twns(:,:,iv,is))!,nsplpo)
          DO il=nsplpa,nsplpe
             ggl=ggng(il)
             xg=SQRT(ggl)*parm%tpiba
             CALL bess(xg,l,mmax,r(1,is),jl)
             !$omp parallel do private(IR)
             DO ir=1,mmax
                fint(ir)=r(ir,is)*betar(ir,indv(iv,is),is)*jl(ir)
             ENDDO
             IF (pslo_com%tlog(is)) THEN
                CALL radlg(mmax,fint,rab(1,is),betax)
             ELSE
                CALL radin(mmax,ncpr1%cmesh(is),fint,betax)
             ENDIF
             twns(il,1,iv,is)=fpi*betax
          ENDDO
          CALL mp_sum(twns(:,:,iv,is),nsplpo,parai%allgrp)
          CALL curv1(nsplpo,ggng(1),twns(1,1,iv,is),0.0_real_8,&
               0.0_real_8,3,twns(1,2,iv,is),fint,0.0_real_8,ierr)
       ENDIF
       iold=indv(iv,is)
    ENDDO
    ! ==--------------------------------------------------------------==
    ! ==        CALCULATE ARRAY  DVAN(IV,JV,IS)    ! IN  RY UNITS     ==
    ! ==--------------------------------------------------------------==
    CALL zeroing(dvan(:,:,is))!,maxsys%nhxs*maxsys%nhxs)
    !$omp parallel do private(IV,JV)
    DO iv=1,nlps_com%ngh(is)
       DO jv=1,nlps_com%ngh(is)
          IF ( MOD(iv-1,istep).EQ.MOD(jv-1,istep) ) THEN
             dvan(iv,jv,is)=0.5_real_8*dion(indv(iv,is),indv(jv,is),is)
             ! 0.5 TO CONVERT RY TO A.U.
          ENDIF
       ENDDO
    ENDDO
    CALL tihalt('   VDBINIT',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE vdbinit
  ! ==================================================================
  SUBROUTINE qinit(jl,fint)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: jl(*), fint(*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'qinit'

    INTEGER                                  :: ierr, il, ir, is, isub, iv, &
                                                jv, l, lqx, lval, mmax, nngh
    REAL(real_8)                             :: fqrad, ggl, qg0, xg
    REAL(real_8), ALLOCATABLE                :: qrd0(:,:), qsp1(:), qsp2(:), &
                                                qsp3(:)

    CALL tiset('     QINIT',isub)
    ALLOCATE(qsp1(nsplpo),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(qsp2(nsplpo),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(qsp3(nsplpo),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(qrd0(nbrx,nbrx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! ==  CALCULATION OF ARRAY  QRAD(IGL,IV,JV,IS)                    ==
    ! ==--------------------------------------------------------------==
    DO is=1,ions1%nsp
       IF (pslo_com%tvan(is)) THEN
          lval=ncpr1%nvales(is)
          mmax=ncpr1%kkbeta(is)
          lqx=2*(lval-1)+1
          IF (qspl1%qspline) THEN
             ! L =  S  P  D  ( F G )      NOTE : Q(R) = R^2*Q(R)
             DO l=1,lqx
                CALL zeroing(qrad(:,:,:,l,is))!,nqdim*nbrx*nbrx)
                DO il=nsplpa,nsplpe
                   ggl=ggnh(il)
                   xg=SQRT(ggl)*parm%tpiba
                   CALL bess(xg,l,mmax,r(1,is),jl)
                   DO iv=1,ncpr1%nbeta(is)
                      DO jv=iv,ncpr1%nbeta(is)
                         IF (pslo_com%tpseu(is)) THEN
                            !$omp parallel do private(IR)
                            DO ir=1,mmax
                               fint(ir)=qrl(ir,iv,jv,l,is)*jl(ir)
                            ENDDO
                         ELSE
                            !$omp parallel do private(IR)
                            DO ir=1,mmax
                               fint(ir)=qfunc(ir,iv,jv,is)*jl(ir)
                            ENDDO
                         ENDIF
                         IF (pslo_com%tlog(is)) THEN
                            CALL radlg(mmax,fint,rab(1,is),fqrad)
                         ELSE
                            CALL radin(mmax,ncpr1%cmesh(is),fint,fqrad)
                         ENDIF
                         qrad(il,iv,jv,l,is)=fpi*fqrad
                         qrad(il,jv,iv,l,is)=fpi*fqrad
                      ENDDO
                   ENDDO
                ENDDO
                ! SPLINE
                CALL mp_sum(qrad(:,:,:,l,is),nqdim*nbrx*nbrx,parai%allgrp)
                DO iv=1,ncpr1%nbeta(is)
                   DO jv=1,ncpr1%nbeta(is)
                      CALL curv1(nsplpo,ggnh(1),&
                           qrad(1,iv,jv,l,is),0.0_real_8,0.0_real_8,3,&
                           qrad(nsplpo+1,iv,jv,l,is),fint,0.0_real_8,ierr)
                   ENDDO
                ENDDO
             ENDDO
             !$omp parallel do private(IV,JV)
             DO iv=1,ncpr1%nbeta(is)
                DO jv=1,ncpr1%nbeta(is)
                   qrd0(iv,jv)=curv2(0.0_real_8,nsplpo,ggnh(1),&
                        qrad(1,iv,jv,1,is),&
                        qrad(nsplpo+1,iv,jv,1,is),0.0_real_8)
                ENDDO
             ENDDO
          ELSE
             CALL dcopy(nsplpo,ggnh(1),1,qsp1(1),1)
             DO l=1,lqx
                DO iv=1,ncpr1%nbeta(is)
                   DO jv=iv,ncpr1%nbeta(is)
                      CALL zeroing(qsp2)!,nsplpo)
                      DO il=nsplpa,nsplpe
                         ggl=ggnh(il)
                         xg=SQRT(ggl)*parm%tpiba
                         CALL bess(xg,l,mmax,r(1,is),jl)
                         IF (pslo_com%tpseu(is)) THEN
                            !$omp parallel do private(IR)
                            DO ir=1,mmax
                               fint(ir)=qrl(ir,iv,jv,l,is)*jl(ir)
                            ENDDO
                         ELSE
                            !$omp parallel do private(IR)
                            DO ir=1,mmax
                               fint(ir)=qfunc(ir,iv,jv,is)*jl(ir)
                            ENDDO
                         ENDIF
                         IF (pslo_com%tlog(is)) THEN
                            CALL radlg(mmax,fint,rab(1,is),fqrad)
                         ELSE
                            CALL radin(mmax,ncpr1%cmesh(is),fint,fqrad)
                         ENDIF
                         qsp2(il)=fpi*fqrad
                      ENDDO
                      CALL mp_sum(qsp2,nsplpo,parai%allgrp)
                      CALL curv1(nsplpo,qsp1(1),qsp2(1),0.0_real_8,0.0_real_8,3,&
                           qsp3(1),fint,0.0_real_8,ierr)
                      !$omp parallel do private(IL,GGL,XG)
                      DO il=1,ncpw%nhgl
                         ggl=gl(il)
                         xg=curv2(ggl,nsplpo,qsp1(1),qsp2(1),qsp3(1),0.0_real_8)
                         qrad(il,iv,jv,l,is)=xg
                         qrad(il,jv,iv,l,is)=xg
                      ENDDO
                      IF (l.EQ.1) THEN
                         qrd0(iv,jv)&
                              =curv2(0._real_8,nsplpo,qsp1(1),qsp2(1),qsp3(1),0._real_8)
                         qrd0(jv,iv)=qrd0(iv,jv)
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          ! ==----------------------------------------------------------==
          ! ==  CALCULATION OF ARRAY QQ    ! ONLY for VAN pseudo        ==
          ! ==  no factor omega, cause qrad is no longer scaled by omega==
          ! ==----------------------------------------------------------==
          nngh=ncpr1%nvales(is)*ncpr1%nvales(is)*nelev(is)
          !$omp parallel do private(IV,JV,QG0)
          DO iv=1,nngh
             DO jv=1,nngh
                CALL qvan1(iv,jv,is,qrd0,qg0)
                qq(iv,jv,is)=qg0
             ENDDO
          ENDDO
       ENDIF
    ENDDO
    ! ==--------------------------------------------------------------==
    DEALLOCATE(qsp1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(qsp2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(qsp3,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(qrd0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    CALL tihalt('     QINIT',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE qinit
  ! ==================================================================

END MODULE vdbinit_utils
