MODULE sd_loc_utils
  USE cppt,                            ONLY: gk,&
                                             inyh,&
                                             nzh,&
                                             rhops,&
                                             scg,&
                                             vps
  USE fftmain_utils,                   ONLY: fwfftn
  USE geq0mod,                         ONLY: geq0
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: parai
  USE sfac,                            ONLY: ei1,&
                                             ei2,&
                                             ei3,&
                                             eigrb
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             ncpw,&
                                             parm

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: sd_loc

CONTAINS

  ! ==================================================================
  SUBROUTINE sd_loc(fion,drhoe,v,eirop1)
    ! ==--------------------------------------------------------------==
    ! ==                        computes                              ==
    ! == local potential  contribution to the dynamical matrix        ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: fion(:,:,:), drhoe(fpar%nnr1)
    COMPLEX(real_8)                          :: v(:), eirop1(ncpw%nhg)

    COMPLEX(real_8)                          :: ei123, rhets, txx, vcgs
    INTEGER                                  :: ia, iat, ig, ig1, ir, is, k
    REAL(real_8)                             :: omtp, tt

! variables
! ==--------------------------------------------------------------==
! TRANSFORM THE DENSITY TO G SPACE

    DO ir=1,fpar%nnr1
       v(ir)=CMPLX(drhoe(ir),0._real_8,kind=real_8)
    ENDDO
    CALL  fwfftn(v,.FALSE.,parai%allgrp)
    ! 
    omtp=2._real_8*parm%omega*parm%tpiba
    iat=0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          iat=iat+1
          DO k=1,3
             ig1=1
             IF (geq0) ig1=2
             tt=0._real_8
             IF (cntl%bigmem) THEN
                DO ig=ig1,ncpw%nhg
                   ei123=eigrb(ig,iat)
                   rhets=CONJG(v(nzh(ig)))
                   vcgs=scg(ig)*CONJG(eirop1(ig)+v(nzh(ig)))
                   txx=(rhops(is,ig)*vcgs+vps(is,ig)*rhets)*&
                        CMPLX(0._real_8,-gk(k,ig),kind=real_8)
                   tt=tt+REAL(ei123*txx)
                ENDDO
             ELSE
                DO ig=ig1,ncpw%nhg
                   ei123=ei1(iat,inyh(1,ig))*ei2(iat,inyh(2,ig))*&
                        ei3(iat,inyh(3,ig))
                   rhets=CONJG(v(nzh(ig)))
                   vcgs=scg(ig)*CONJG(eirop1(ig)+v(nzh(ig)))
                   txx=(rhops(is,ig)*vcgs+vps(is,ig)*rhets)*&
                        CMPLX(0._real_8,-gk(k,ig),kind=real_8)
                   tt=tt+REAL(ei123*txx)
                ENDDO
             ENDIF
             fion(k,ia,is)=fion(k,ia,is)+omtp*tt
          ENDDO
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE sd_loc
  ! ==================================================================

END MODULE sd_loc_utils
