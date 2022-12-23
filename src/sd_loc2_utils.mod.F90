MODULE sd_loc2_utils
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
                                             maxsys,&
                                             ncpw,&
                                             parm

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: sd_loc2

CONTAINS

  ! ==================================================================
  SUBROUTINE sd_loc2(sder,rhoo,v,iato,eirop)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! == LOCAL CONTRIBUTION TO DYNAMICAL MATRIX (DIAGONAL PART)       ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: sder(3*ions1%nat,3*ions1%nat),&
                                                rhoo(fpar%nnr1)
    COMPLEX(real_8)                          :: v(:)
    INTEGER                                  :: iato(maxsys%nax,maxsys%nsx)
    COMPLEX(real_8)                          :: eirop(ncpw%nhg)

    COMPLEX(real_8)                          :: ei123, rhet, rhets, ttt, vcgs
    INTEGER                                  :: ia, iat, ig, ig1, ir, is, isa
    REAL(real_8)                             :: gx, gy, gz, omtp, ta

! Variables
! ==--------------------------------------------------------------==
! TRANSFORM THE DENSITY TO G SPACE

    !$omp parallel do private(IR)
    DO ir=1,fpar%nnr1
       v(ir)=CMPLX(rhoo(ir),0._real_8,kind=real_8)
    ENDDO
    CALL  fwfftn(v,.FALSE.,parai%allgrp)
    ! 
    omtp=2._real_8*parm%omega*parm%tpiba2
    ig1=1
    IF (geq0) ig1=2
    DO ig=ig1,ncpw%nhg
       rhet=v(nzh(ig))
       rhets=CONJG(v(nzh(ig)))
       vcgs=scg(ig)*CONJG(rhet+eirop(ig))
       gx=gk(1,ig)
       gy=gk(2,ig)
       gz=gk(3,ig)
       isa=0
       DO is=1,ions1%nsp
          ttt=(rhops(is,ig)*vcgs+vps(is,ig)*rhets)
          !CDIR NODEP
          DO ia=1,ions0%na(is)
             isa=isa+1
             iat=3*iato(ia,is)
             IF (cntl%bigmem) THEN
                ei123=eigrb(ig,isa)
             ELSE
                ei123=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
                     ei3(isa,inyh(3,ig))
             ENDIF
             ta=-REAL(ei123*ttt)*omtp
             sder(iat+1,iat+1)=sder(iat+1,iat+1)+ta*gx*gx
             sder(iat+2,iat+1)=sder(iat+2,iat+1)+ta*gy*gx
             sder(iat+3,iat+1)=sder(iat+3,iat+1)+ta*gz*gx
             sder(iat+1,iat+2)=sder(iat+1,iat+2)+ta*gx*gy
             sder(iat+2,iat+2)=sder(iat+2,iat+2)+ta*gy*gy
             sder(iat+3,iat+2)=sder(iat+3,iat+2)+ta*gz*gy
             sder(iat+1,iat+3)=sder(iat+1,iat+3)+ta*gx*gz
             sder(iat+2,iat+3)=sder(iat+2,iat+3)+ta*gy*gz
             sder(iat+3,iat+3)=sder(iat+3,iat+3)+ta*gz*gz
          ENDDO
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE sd_loc2
  ! ==================================================================

END MODULE sd_loc2_utils
