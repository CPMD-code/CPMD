MODULE cofor_utils
  USE cppt,                            ONLY: gk,&
                                             inyh
  USE ions,                            ONLY: ions1
  USE kinds,                           ONLY: real_8
  USE nlcc,                            ONLY: corel,&
                                             rhoc,&
                                             vnlcc,&
                                             vnlt
  USE sfac,                            ONLY: ei1,&
                                             ei2,&
                                             ei3,&
                                             eigrb
  USE system,                          ONLY: cntl,&
                                             iatpt,&
                                             ncpw,&
                                             parm

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: cofor

CONTAINS

  ! ==================================================================
  SUBROUTINE cofor(fion,vpot)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! == FORCES ON IONS DUE TO CORE CHARGES                           ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: fion(:,:,:)
    COMPLEX(real_8)                          :: vpot(ncpw%nhg,*)

    COMPLEX(real_8)                          :: ei123, vxc
    INTEGER                                  :: ia, ig, is, isa
    REAL(real_8)                             :: omtp, vcgs

    omtp=2._real_8*parm%omega*parm%tpiba
    !$omp parallel do private(ISA,IA,IS,IG,EI123,VXC,VCGS) shared(OMTP)
    DO isa=1,ions1%nat
       ia=iatpt(1,isa)
       is=iatpt(2,isa)
       IF (corel%tnlcc(is)) THEN
          DO ig=1,ncpw%nhg
             IF (cntl%bigmem) THEN
                ei123=eigrb(ig,isa)
             ELSE
                ei123=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
                     ei3(isa,inyh(3,ig))
             ENDIF
             vxc=vpot(ig,1)-vnlt(ig)-vnlcc(ig,1)
             IF (cntl%tlsd) vxc=0.5_real_8*(vxc+vpot(ig,2)-vnlt(ig)-vnlcc(ig,2))
             vcgs=-AIMAG(CONJG(vxc)*ei123*rhoc(ig,is))
             fion(1,ia,is)=fion(1,ia,is)+gk(1,ig)*vcgs*omtp
             fion(2,ia,is)=fion(2,ia,is)+gk(2,ig)*vcgs*omtp
             fion(3,ia,is)=fion(3,ia,is)+gk(3,ig)*vcgs*omtp
          ENDDO
       ENDIF
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE cofor
  ! ==================================================================

END MODULE cofor_utils
