MODULE eind_loc_utils
  USE cppt,                            ONLY: gk,&
                                             inyh,&
                                             nzh,&
                                             rhops,&
                                             scg,&
                                             vps
  USE fftmain_utils,                   ONLY: fwfftn
  USE geq0mod,                         ONLY: geq0
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

  PUBLIC :: eind_loc

CONTAINS

  ! ==================================================================
  SUBROUTINE eind_loc(eind,is,isa,ka,rhoo,v,eirop)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! == LOCAL CONTRIBUTION TO THE SECOND ORDER ENERGY (CONSTANT PART)==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: eind
    INTEGER                                  :: is, isa, ka
    REAL(real_8)                             :: rhoo(fpar%nnr1)
    COMPLEX(real_8)                          :: v(:), eirop(ncpw%nhg)

    COMPLEX(real_8)                          :: ei123, rhet, rhets, ttt, vcgs
    INTEGER                                  :: ig, ig1, ir
    REAL(real_8)                             :: gka, omtp, ta

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
    IF (cntl%bigmem) THEN  ! do not break the loop for good vectorzation
       DO ig=ig1,ncpw%nhg
          rhet=v(nzh(ig))
          rhets=CONJG(rhet)
          gka=gk(ka,ig)
          ei123=eigrb(ig,isa)
          vcgs=scg(ig)*CONJG(rhet+eirop(ig)+rhops(is,ig)*ei123)
          ttt=(rhops(is,ig)*vcgs+vps(is,ig)*rhets)
          ta=-REAL(ei123*ttt)*omtp
          eind=eind+ta*gka*gka
       ENDDO
    ELSE
       DO ig=ig1,ncpw%nhg
          rhet=v(nzh(ig))
          rhets=CONJG(rhet)
          gka=gk(ka,ig)
          ei123=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
               ei3(isa,inyh(3,ig))
          vcgs=scg(ig)*CONJG(rhet+eirop(ig)+rhops(is,ig)*ei123)
          ttt=(rhops(is,ig)*vcgs+vps(is,ig)*rhets)
          ta=-REAL(ei123*ttt)*omtp
          eind=eind+ta*gka*gka
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE eind_loc
  ! ==================================================================

END MODULE eind_loc_utils
