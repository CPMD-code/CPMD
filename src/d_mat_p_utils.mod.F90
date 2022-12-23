MODULE d_mat_p_utils
  USE cnst,                            ONLY: pi
  USE cppt,                            ONLY: gk,&
                                             inyh,&
                                             rhops,&
                                             scg
  USE error_handling,                  ONLY: stopgm
  USE geq0mod,                         ONLY: geq0
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: int_1,&
                                             int_2,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8
  USE nlps,                            ONLY: imagp,&
                                             ndfnl,&
                                             nghtol,&
                                             nlps_com,&
                                             wsg
  USE parac,                           ONLY: parai,&
                                             paral
  USE pbc_utils,                       ONLY: pbc
  USE pslo,                            ONLY: pslo_com
  USE ragg,                            ONLY: raggio
  USE sfac,                            ONLY: ei1,&
                                             ei2,&
                                             ei3,&
                                             eigrb
  USE sgpp,                            ONLY: sgpp1,&
                                             sgpp2
  USE special_functions,               ONLY: cp_erfc
  USE system,                          ONLY: cntl,&
                                             maxsys,&
                                             ncpw,&
                                             parap,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  !!public :: d_mat_diag_loc
  PUBLIC :: d_mat_diag_nonloc
  PUBLIC :: d_mat_nonloc
  PUBLIC :: d_mat_diag_real
  PUBLIC :: d_mat_real
  PUBLIC :: d_mat_loc0
!!!public :: d_mat_locps

CONTAINS

  ! ==================================================================
  SUBROUTINE d_mat_diag_nonloc(sder,ddfnl00,fnl00,dfnl00,&
       f,wk,isp,k,k2,iat,nstate,nkpoint)
    ! ==--------------------------------------------------------------==
    ! ==                        computes                              ==
    ! ==  the non-local potential contribution to the diagonal        ==
    ! ==  part of the dynamical matrix                                ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: sder
    INTEGER                                  :: isp, k, k2, iat, nstate, &
                                                nkpoint
    REAL(real_8) :: fnl00(imagp,ions1%nat,maxsys%nhxs,nstate,nkpoint), &
      f(nstate,nkpoint), wk(nkpoint), &
      dfnl00(imagp,ions1%nat,maxsys%nhxs,3,ndfnl,nkpoint), &
      ddfnl00(imagp,ions1%nat,maxsys%nhxs,3,3,ndfnl,nkpoint)

    INTEGER                                  :: ik, is12, istate, iv, jv
    LOGICAL                                  :: nonzero_overlap
    REAL(real_8)                             :: energy_value, sum, weight

! ==--------------------------------------------------------------==

    IF (pslo_com%tvan(isp)) THEN
       CALL stopgm('d_mat_diag_nonloc',&
            'vanderbilt not implemented',& 
            __LINE__,__FILE__)
    ELSEIF (sgpp1%tsgp(isp)) THEN
       DO ik=1,nkpoint
          DO iv=1,nlps_com%ngh(isp)
             DO jv=1,nlps_com%ngh(isp)
                nonzero_overlap = (NGHtoL(iv,isp) .EQ. nghtol(jv,isp))&
                     .AND. (sgpp2%lpval(iv,isp) .EQ. sgpp2%lpval(jv,isp))
                IF (nonzero_overlap) THEN
                   sum = 0._real_8
                   DO istate = parap%nst12(parai%mepos,1),parap%nst12(parai%mepos,2)
                      is12 = istate - parap%nst12(parai%mepos,1)+1

                      sum  = sum + 2._real_8*f(istate,ik) *&
                           ( ddfnl00 (1,iat,iv,k,k2,is12,ik)&
                           * fnl00   (1,iat,jv,istate,ik)&
                           + dfnl00  (1,iat,iv,k,is12,ik)&
                           * dfnl00  (1,iat,jv,k2,is12,ik))

                   ENDDO
                   energy_value  = sgpp2%hlsg(sgpp2%lfval(iv,isp),sgpp2%lfval(jv,isp),&
                        NGHtoL(iv,isp)+1,isp)

                   sder = sder + sum * energy_value
                ENDIF     ! IF overlap between Proj[iv] and Proj[jv]
             ENDDO         ! jv
          ENDDO             ! iv
       ENDDO                 ! k-points
    ELSE                      ! BHS and related.
       DO ik=1,nkpoint
          DO iv=1,nlps_com%ngh(isp)
             DO istate=parap%nst12(parai%mepos,1),parap%nst12(parai%mepos,2)
                weight=wk(ik)*f(istate,ik)*2._real_8*wsg(isp,iv)
                is12=istate-parap%nst12(parai%mepos,1)+1
                sder=sder + weight*&
                     (ddfnl00   (1,iat,iv,k,k2,is12,ik)&
                     *   fnl00  (1,iat,iv,istate,ik)&
                     + dfnl00   (1,iat,iv,k,is12,ik)&
                     *   dfnl00 (1,iat,iv,k2,is12,ik))
             ENDDO
          ENDDO
       ENDDO
    ENDIF

    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE d_mat_diag_nonloc
  ! ==================================================================
  SUBROUTINE d_mat_nonloc(sder,fnl,dfnl,fnl00,dfnl00,&
       f,wk,isp,k,k2,iat,nstate,nkpoint)
    ! ==---------------------------------------------------------------==
    ! ==                        computes                               ==
    ! ==  the non-local potential contribution to the dynamical matrix ==
    ! ==---------------------------------------------------------------==
    REAL(real_8)                             :: sder
    INTEGER                                  :: isp, k, k2, iat, nstate, &
                                                nkpoint
    REAL(real_8) :: fnl(imagp,ions1%nat,maxsys%nhxs,nstate,nkpoint), &
      fnl00(imagp,ions1%nat,maxsys%nhxs,nstate,nkpoint), f(nstate,nkpoint), &
      wk(nkpoint), dfnl00(imagp,ions1%nat,maxsys%nhxs,3,ndfnl,nkpoint), &
      dfnl(imagp,ions1%nat,maxsys%nhxs,3,ndfnl,nkpoint)

    INTEGER                                  :: ik, is12, istate, iv, jv
    LOGICAL                                  :: nonzero_overlap
    REAL(real_8)                             :: energy_value, sum, weight

! variables
! ==--------------------------------------------------------------==

    IF (pslo_com%tvan(isp)) THEN
       CALL stopgm('d_mat_diag_nonloc',&
            'vanderbilt not implemented',& 
            __LINE__,__FILE__)
    ELSEIF (sgpp1%tsgp(isp)) THEN
       DO ik=1,nkpoint
          DO iv=1,nlps_com%ngh(isp)
             DO jv=1,nlps_com%ngh(isp)
                nonzero_overlap = (NGHtoL(iv,isp) .EQ. nghtol(jv,isp))&
                     .AND. (sgpp2%lpval(iv,isp) .EQ. sgpp2%lpval(jv,isp))
                IF (nonzero_overlap) THEN
                   sum = 0._real_8
                   DO istate = parap%nst12(parai%mepos,1),parap%nst12(parai%mepos,2)
                      is12 = istate - parap%nst12(parai%mepos,1)+1

                      sum  = sum + 2._real_8*f(istate,ik) *&
                           ( dfnl00  (1,iat,iv,k2,is12,ik)&
                           * fnl     (1,iat,jv,istate,ik)&
                           + dfnl    (1,iat,iv,k2,is12,ik)&
                           * fnl00   (1,iat,jv,istate,ik))

                   ENDDO
                   energy_value  = sgpp2%hlsg(sgpp2%lfval(iv,isp),sgpp2%lfval(jv,isp),&
                        NGHtoL(iv,isp)+1,isp)
                   sder = sder + sum * energy_value
                ENDIF     ! IF overlap between Proj[iv] and Proj[jv]
             ENDDO         ! jv
          ENDDO             ! iv
       ENDDO                 ! k-points
    ELSE                      ! BHS and related.
       DO ik=1,nkpoint
          DO iv=1,nlps_com%ngh(isp)
             DO istate=parap%nst12(parai%mepos,1),parap%nst12(parai%mepos,2)
                weight=wk(ik)*f(istate,ik)*2._real_8*wsg(isp,iv)
                is12=istate-parap%nst12(parai%mepos,1)+1
                sder = sder + weight*&
                     (dfnl00    (1,iat,iv,k2,is12,ik)&
                     *   fnl    (1,iat,iv,istate,ik)&
                     + dfnl     (1,iat,iv,k2,is12,ik)&
                     *   fnl00  (1,iat,iv,istate,ik))
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE d_mat_nonloc
  ! ==================================================================
  SUBROUTINE d_mat_diag_real(sder,is,ia,k,k2,tau0)
    ! ==--------------------------------------------------------------==
    ! ==                        computes                              ==
    ! == the residual part of the ion-ion interaction (esr)           ==
    ! == due to the overlap of the smeared ionic charge densities     ==
    ! == (ionic point charges replaced by gaussian charge distrib.)   ==
    ! == corresponding to different atoms.                            ==
    ! == esr depends only on tau0 (raggio and valence charges)        ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: sder
    INTEGER                                  :: is, ia, k, k2
    REAL(real_8)                             :: tau0(:,:,:)

    INTEGER                                  :: j, m
    LOGICAL                                  :: tzero
    REAL(real_8)                             :: arg, derfcdx2, derre, derre2, &
                                                erre2, rckj, repand, rlm, &
                                                rxlm(3), xderfcdx, xlm, ylm, &
                                                zlm, zv2, xlm_, ylm_,zlm_

    IF (paral%parent) THEN
       DO j=1,ions1%nsp
          zv2=ions0%zv(is)*ions0%zv(j)
          rckj=SQRT(raggio(is)*raggio(is)+raggio(j)*raggio(j))
          DO m=1,ions0%na(j)
             IF ((is.EQ.j).AND.(ia.EQ.m)) THEN
                tzero=.TRUE.
             ELSE
                tzero=.FALSE.
                xlm_=tau0(1,ia,is)-tau0(1,m,j)
                ylm_=tau0(2,ia,is)-tau0(2,m,j)
                zlm_=tau0(3,ia,is)-tau0(3,m,j)
                CALL pbc(xlm_,ylm_,zlm_,xlm,ylm,zlm,1,parm%apbc,parm%ibrav)

             ENDIF
             IF (.NOT.tzero) THEN
                rxlm(1)=xlm
                rxlm(2)=ylm
                rxlm(3)=zlm
                erre2=rxlm(1)**2+rxlm(2)**2+rxlm(3)**2
                rlm=SQRT(erre2)
                arg=rlm/rckj
                xderfcdx=-zv2*cp_erfc(arg)/rlm-(2._real_8*zv2/SQRT(pi))*&
                     EXP(-arg*arg)/rckj
                derfcdx2=-(2._real_8/erre2)*xderfcdx+&
                     (4._real_8/rckj**3)*zv2*EXP(-arg*arg)/SQRT(pi)

                derre=rxlm(k)
                derre2=rxlm(k2)
                repand=xderfcdx/erre2

                sder = sder&
                     +repand*(-derre*derre2/erre2)
                IF (k.EQ.k2)  sder = sder + repand
                sder = sder + derfcdx2*(derre*derre2/erre2)
             ENDIF
          ENDDO
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE d_mat_diag_real
  ! ==================================================================
  SUBROUTINE d_mat_real(sder,is,ia,k,k2,is2,ia2,tau0)
    ! ==--------------------------------------------------------------==
    ! ==                        computes                              ==
    ! == the residual part of the ion-ion interaction (esr)           ==
    ! == due to the overlap of the smeared ionic charge densities     ==
    ! == (ionic point charges replaced by gaussian charge distrib.)   ==
    ! == corresponding to different atoms.                            ==
    ! == esr depends only on tau0 (raggio and valence charges)        ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: sder
    INTEGER                                  :: is, ia, k, k2, is2, ia2
    REAL(real_8)                             :: tau0(:,:,:)

    LOGICAL                                  :: tzero
    REAL(real_8)                             :: arg, derfcdx2, derre, derre2, &
                                                erre2, rckj, repand, rlm, &
                                                rxlm(3), xderfcdx, xlm, ylm, &
                                                zlm, zv2, xlm_, ylm_,zlm_

    IF (paral%parent) THEN
       zv2=ions0%zv(is)*ions0%zv(is2)
       rckj=SQRT(raggio(is)*raggio(is)+raggio(is2)*raggio(is2))
       IF (is.EQ.is2.AND.ia.EQ.ia2) THEN
          tzero=.TRUE.
       ELSE
          tzero=.FALSE.
          xlm_=tau0(1,ia,is)-tau0(1,ia2,is2)
          ylm_=tau0(2,ia,is)-tau0(2,ia2,is2)
          zlm_=tau0(3,ia,is)-tau0(3,ia2,is2)
          CALL pbc(xlm_,ylm_,zlm_,xlm,ylm,zlm,1,parm%apbc,parm%ibrav)
       ENDIF
       IF (.NOT.tzero) THEN
          rxlm(1)=xlm
          rxlm(2)=ylm
          rxlm(3)=zlm
          erre2=rxlm(1)**2+rxlm(2)**2+rxlm(3)**2
          rlm=SQRT(erre2)
          arg=rlm/rckj
          xderfcdx=-zv2*cp_erfc(arg)/rlm-(2._real_8*zv2/SQRT(pi))*&
               EXP(-arg*arg)/rckj
          derfcdx2=-(2._real_8/erre2)*xderfcdx&
               +(4._real_8/rckj**3)*zv2*EXP(-arg*arg)/SQRT(pi)

          derre=rxlm(k)
          derre2=rxlm(k2)

          repand=xderfcdx/erre2
          sder = sder + repand*(derre*derre2/erre2)
          IF (k.EQ.k2)  sder = sder - repand
          sder = sder - derfcdx2*(derre*derre2/erre2)
       ENDIF
    ENDIF
    ! ==================================================================
    RETURN
  END SUBROUTINE d_mat_real
  ! ==================================================================
  SUBROUTINE d_mat_loc0(sder,eirop1,iat,k2,is)
    ! ==--------------------------------------------------------------==
    ! ==                        computes                              ==
    ! == local potential  contribution to the dynamical matrix        ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: sder
    COMPLEX(real_8)                          :: eirop1(ncpw%nhg)
    INTEGER                                  :: iat, k2, is

    COMPLEX(real_8)                          :: ei123, txx, vcgs
    INTEGER                                  :: ig, ig1
    REAL(real_8)                             :: omtp

! ==--------------------------------------------------------------==

    omtp=2._real_8*parm%omega*parm%tpiba
    ig1=1
    IF (geq0) ig1=2
    DO ig=ig1,ncpw%nhg
       vcgs = scg(ig) * CONJG(eirop1(ig))
       txx  = rhops(is,ig) * vcgs * CMPLX(0._real_8,-gk(k2,ig),kind=real_8)
       IF (cntl%bigmem) THEN
          ei123=eigrb(ig,iat)
       ELSE
          ei123=ei1(iat,inyh(1,ig))*ei2(iat,inyh(2,ig))*&
               ei3(iat,inyh(3,ig))
       ENDIF
       sder = sder + REAL( ei123*txx ) * omtp
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE d_mat_loc0
  ! ==================================================================
END MODULE d_mat_p_utils

SUBROUTINE d_mat_locps(sder,drhoe,iat,k2,is)
  ! ==--------------------------------------------------------------==
  ! ==                        computes                              ==
  ! == local potential  contribution to the dynamical matrix        ==
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE system , ONLY:cntl,ncpw,parm
  USE parac, ONLY : paral,parai
  USE cppt , ONLY:gk,inyh,rhops,scg,vps
  USE geq0mod , ONLY:geq0
  USE sfac , ONLY:ei1,ei2,ei3,eigrb
  IMPLICIT NONE
  REAL(real_8)                               :: sder
  COMPLEX(real_8)                            :: drhoe(ncpw%nhg)
  INTEGER                                    :: iat, k2, is

  COMPLEX(real_8)                            :: ei123, txx
  INTEGER                                    :: ig, ig1
  REAL(real_8)                               :: omtp

! ==--------------------------------------------------------------==

  omtp=2._real_8*parm%omega*parm%tpiba
  ig1=1
  IF (geq0) ig1=2
  DO ig=ig1,ncpw%nhg
     txx = rhops(is,ig) * scg(ig) + vps(is,ig)
     txx = txx * CONJG(drhoe(ig)) * CMPLX(0._real_8,-gk(k2,ig),kind=real_8)
     IF (cntl%bigmem) THEN
        ei123=eigrb(ig,iat)
     ELSE
        ei123=ei1(iat,inyh(1,ig))*ei2(iat,inyh(2,ig))*&
             ei3(iat,inyh(3,ig))
     ENDIF
     sder = sder + REAL(ei123*txx) * omtp
  ENDDO
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE d_mat_locps
! ==================================================================



! ==================================================================
SUBROUTINE d_mat_diag_loc(sder,rho0G,eirop,eirop1,v1_loc,k2)
  ! ==--------------------------------------------------------------==
  ! ==                        computes                              ==
  ! == the part of the dynamical matrix related to the local part   ==
  ! == of the pseudopotential 
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE system , ONLY:ncpw,parm
  USE parac, ONLY : paral,parai
  USE cppt , ONLY:gk,scg
  USE geq0mod , ONLY:geq0
  IMPLICIT NONE
  REAL(real_8)                               :: sder
  COMPLEX(real_8)                            :: rho0G(*), eirop(*), &
                                                eirop1(*), v1_loc(*)
  INTEGER                                    :: k2

  COMPLEX(real_8)                            :: gk2, rhogs, txx
  INTEGER                                    :: ig, ig1
  REAL(real_8)                               :: factor, omtp

  omtp=2._real_8*parm%omega*parm%tpiba
  ig1=1
  IF (geq0) ig1=2
  DO ig=ig1,ncpw%nhg
     factor = REAL(scg(ig))
     rhogs=factor * CONJG(rho0G(ig)+eirop(ig))

     gk2=CMPLX(0._real_8,-gk(k2,ig),kind=real_8)

     txx = eirop1(ig)*rhogs + v1_loc(ig)*CONJG(rho0G(ig))
     txx = txx * gk2
     sder=sder + REAL(txx)*omtp
  ENDDO
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE d_mat_diag_loc
