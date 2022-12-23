MODULE k_hesele_utils
  USE cppt,                            ONLY: twnl
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE kpnt,                            ONLY: hgkm,&
                                             hgkp
  USE mp_interface,                    ONLY: mp_bcast
  USE nlps,                            ONLY: nghtol,&
                                             nlps_com,&
                                             wsg
  USE parac,                           ONLY: parai
  USE sgpp,                            ONLY: sgpp1,&
                                             sgpp2
  USE simulmod,                        ONLY: vploc
  USE system,                          ONLY: cntl,&
                                             cntr,&
                                             ncpw,&
                                             nkpt,&
                                             parm

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: k_hesele

CONTAINS

  ! ==================================================================
  SUBROUTINE k_hesele (svar2,vpp,ik)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: svar2, vpp(nkpt%ngwk)
    INTEGER                                  :: ik

    INTEGER                                  :: ig, is, isa0, iv, j, jv, ki, &
                                                kj, l, l2, li, lj
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
       DO ig=1,nkpt%ngwk
          vpp(ig)=vp
          isa0=0
          DO is=1,ions1%nsp
             IF (sgpp1%tsgp(is)) THEN
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
                         vpp(ig)=vpp(ig)+&
                              dd*twnl(ig,iv,is,ik)*twnl(ig,jv,is,ik)
                      ENDIF
                   ENDDO
                ENDDO
             ELSE
                DO iv=1,nlps_com%ngh(is)
                   vpp(ig)=vpp(ig)+ions0%na(is)*wsg(is,iv)&
                        *twnl(ig,iv,is,ik)*twnl(ig,iv,is,ik)
                ENDDO
             ENDIF
             isa0=isa0+ions0%na(is)
          ENDDO
       ENDDO
       ! ==--------------------------------------------------------------==
       ! ==  Kinetic Energy                                              ==
       ! ==--------------------------------------------------------------==
       ftpi=0.5_real_8*parm%tpiba2
       DO j=1,ncpw%ngw
          vpp(j)=ftpi*hgkp(j,ik)+vpp(j)
          vpp(j+ncpw%ngw)=ftpi*hgkm(j,ik)+vpp(j+ncpw%ngw)
          emp=cntr%hthrs
          IF (vpp(j).LT.emp) vpp(j)=emp
          IF (vpp(j+ncpw%ngw).LT.emp) vpp(j+ncpw%ngw)=emp
       ENDDO
       !$omp parallel do private(J)
       DO j=1,nkpt%ngwk
          vpp(j)=1.0_real_8/vpp(j)
       ENDDO
    ELSE
       !$omp parallel do private(J)
       DO j=1,nkpt%ngwk
          vpp(j)=svar2
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE k_hesele
  ! ==================================================================

END MODULE k_hesele_utils
