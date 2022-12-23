MODULE qvan2_utils
  USE aavan,                           ONLY: ap,&
                                             indv,&
                                             lpl,&
                                             lpx
  USE cppt,                            ONLY: hg,&
                                             igl,&
                                             qrad,&
                                             ylmb
  USE cvan,                            ONLY: nelev
  USE error_handling,                  ONLY: stopgm
  USE fitpack_utils,                   ONLY: curv2
  USE kinds,                           ONLY: real_8
  USE nlps,                            ONLY: nlps_com
  USE qspl,                            ONLY: ggnh,&
                                             nsplpo,&
                                             qspl1
  USE system,                          ONLY: nbrx,&
                                             ncpw,&
                                             nlx,&
                                             parm
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: qvan2

CONTAINS

  ! =====================================================================
  SUBROUTINE qvan2(iv,jv,is,qg)
    ! ==-----------------------------------------------------------------==
    ! ==  Q(G,L,K) = SUM_LM (-I)^L AP(LM,L,K) YR_LM(G^) QRAD(G,L,L,K)    ==
    ! ==-----------------------------------------------------------------==
    INTEGER                                  :: iv, jv, is
    COMPLEX(real_8)                          :: qg(ncpw%nhg)

    CHARACTER(*), PARAMETER                  :: procedureN = 'qvan2'

    COMPLEX(real_8)                          :: sig
    INTEGER                                  :: i, ierr, ig, istep, ivl, ivs, &
                                                jvl, jvs, l, lp
    INTEGER, SAVE                            :: ifirst = 0
    REAL(real_8)                             :: sigi, sigr, vol
    REAL(real_8), ALLOCATABLE, SAVE          :: spline(:)

! ==-----------------------------------------------------------------==

    IF (ifirst.EQ.0) THEN
       ifirst=1
       IF (qspl1%qspline) THEN
          ALLOCATE(spline(ncpw%nhg),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
    ENDIF
    ! ==-----------------------------------------------------------------==
    ! ==  IV  = 1..8    ! s_1 p_x1 p_z1 p_y1 s_2 p_x2 p_z2 p_y2          ==
    ! ==  IVS = 1..4    ! s_1 s_2 p_1 p_2                                ==
    ! ==  IVL = 1..4    ! s p_x p_z p_y                                  ==
    ! ==                                                                 ==
    ! == NOTE :   IV  = 1..8 (sppp sppp)  IVS = 1..4 (sspp) OR 1..2 (sp) ==
    ! ==    IVL = 1..4 (sppp)                                            ==
    ! ==-----------------------------------------------------------------==
    vol=1._real_8/parm%omega
    ivs=indv(iv,is)
    jvs=indv(jv,is)
    istep=nlps_com%ngh(is)/nelev(is)
    ivl=1+MOD(iv-1,istep)
    jvl=1+MOD(jv-1,istep)
    IF (ivl.GT.nlx)  CALL stopgm(' QVAN ',' IVL.GT.NLX  ',& 
         __LINE__,__FILE__)
    IF (jvl.GT.nlx)  CALL stopgm(' QVAN ',' JVL.GT.NLX  ',& 
         __LINE__,__FILE__)
    IF (ivs.GT.nbrx) CALL stopgm(' QVAN ',' IVS.GT.NBRX ',& 
         __LINE__,__FILE__)
    IF (jvs.GT.nbrx) CALL stopgm(' QVAN ',' JVS.GT.NBRX ',& 
         __LINE__,__FILE__)
    CALL zeroing(qg)!,SIZE(qg))
    DO i=1,lpx(ivl,jvl)
       lp=lpl(ivl,jvl,i)
       ! EXTRACTION OF ANGULAR MOMENTUM L FROM LP:
       IF (lp.GE.26) CALL stopgm(' QVAN ',' LP.GE.26 ',& 
            __LINE__,__FILE__)
       l=NINT(SQRT(REAL(lp,kind=real_8))+0.499999_real_8)
       ! 
       sig=(0._real_8,-1._real_8)**(l-1)
       sigr=REAL(sig)
       sigi=AIMAG(sig)
       IF (qspl1%qspline) THEN
          ! SPLINE
          DO ig=1,ncpw%nhg
             spline(ig)=vol*&
                  curv2(hg(ig),nsplpo,ggnh(1),&
                  qrad(1,ivs,jvs,l,is),qrad(nsplpo+1,ivs,jvs,l,is),0.0_real_8)
          ENDDO
          IF (ABS(sigr).GT.0.5_real_8) THEN
             sigr=sigr*ap(lp,ivl,jvl)
             DO ig=1,ncpw%nhg
                qg(ig)=CMPLX(  REAL(qg(ig),KIND=real_8) + sigr*ylmb(ig,lp,1)*spline(ig), &
                     & AIMAG(qg(ig)) , KIND=real_8 )
                !qg(1,ig)=qg(1,ig)+sigr*ylmb(ig,lp,1)*spline(ig)
             ENDDO
          ELSE
             sigi=sigi*ap(lp,ivl,jvl)
             DO ig=1,ncpw%nhg
                qg(ig)=CMPLX( REAL(qg(ig),KIND=real_8), AIMAG(qg(ig))+sigi*ylmb(ig,lp,1)*spline(ig) , &
                     & KIND=real_8 )
                !qg(2,ig)=qg(2,ig)+sigi*ylmb(ig,lp,1)*spline(ig)
             ENDDO
          ENDIF
       ELSE
          IF (ABS(sigr).GT.0.5_real_8) THEN
             sigr=sigr*ap(lp,ivl,jvl)
             DO ig=1,ncpw%nhg
                qg(ig)=CMPLX( REAL(qg(ig),KIND=real_8)+vol*sigr*ylmb(ig,lp,1)*qrad(igl(ig),ivs,jvs,l,is), &
                     & AIMAG(qg(ig)) , KIND=real_8 )
                !qg(1,ig)=qg(1,ig)+vol*sigr*ylmb(ig,lp,1)*qrad(igl(ig),ivs,jvs,l,is)
             ENDDO
          ELSE
             sigi=sigi*ap(lp,ivl,jvl)
             DO ig=1,ncpw%nhg
                qg(ig)=CMPLX( REAL(qg(ig),KIND=real_8), &
                     & AIMAG(qg(ig))+vol*sigi*ylmb(ig,lp,1)*qrad(igl(ig),ivs,jvs,l,is), KIND=real_8 )
                !qg(2,ig)=qg(2,ig)+vol*sigi*ylmb(ig,lp,1)*qrad(igl(ig),ivs,jvs,l,is)
             ENDDO
          ENDIF
       ENDIF
    ENDDO
    ! ==-----------------------------------------------------------------==
    RETURN
  END SUBROUTINE qvan2
  ! =====================================================================

END MODULE qvan2_utils
