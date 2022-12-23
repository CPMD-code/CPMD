MODULE nlccstr_utils
  USE copot_utils,                     ONLY: drhocore
  USE cppt,                            ONLY: inyh,&
                                             nzh
  USE dotp_utils,                      ONLY: dotp
  USE fftmain_utils,                   ONLY: fwfftn
  USE kinds,                           ONLY: real_8
  USE nlcc,                            ONLY: drhoc
  USE parac,                           ONLY: parai
  USE sfac,                            ONLY: ei1,&
                                             ei2,&
                                             ei3,&
                                             eigrb
  USE strs,                            ONLY: decc
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             ncpw
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: zgthr

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: nlccstr

CONTAINS

  ! ==================================================================
  SUBROUTINE nlccstr(v,vgc)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==           EXCHANGE AND CORRELATION ENERGY CONTRIBUTION       ==
    ! ==                 OF CORE CHARGE TO STRESS TENSOR              ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: v(fpar%nnr1,*), vgc(:)

    INTEGER                                  :: isub, kk

! ==--------------------------------------------------------------==

    CALL tiset('   NLCCSTR',isub)
    IF (cntl%tlsd) THEN
       CALL  fwfftn(v(:,1),.FALSE.,parai%allgrp)
       CALL zgthr(ncpw%nhg,v,vgc,nzh)
       CALL dcopy(2*ncpw%nhg,vgc,1,v(1,1),1)
       CALL  fwfftn(v(:,2),.FALSE.,parai%allgrp)
       CALL zgthr(ncpw%nhg,v(1,2),vgc,nzh)
       CALL dcopy(2*ncpw%nhg,vgc,1,v(1,2),1)
    ELSE
       CALL  fwfftn(v(:,1),.FALSE.,parai%allgrp)
       CALL zgthr(ncpw%nhg,v,vgc,nzh)
       CALL dcopy(2*ncpw%nhg,vgc,1,v,1)
    ENDIF
    DO kk=1,6
       CALL drhocore(vgc,drhoc(1,1,kk),eigrb,ei1,ei2,ei3,inyh)
       IF (cntl%tlsd) THEN
          decc(kk)=decc(kk)+dotp(ncpw%nhg,vgc,v(:,1))
          decc(kk)=decc(kk)+dotp(ncpw%nhg,vgc,v(:,2))
       ELSE
          decc(kk)=decc(kk)+dotp(ncpw%nhg,vgc,v(:,1))
       ENDIF
    ENDDO
    CALL tihalt('   NLCCSTR',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE nlccstr
  ! ==================================================================

END MODULE nlccstr_utils
