MODULE rekine_utils
  USE dotp_utils,                      ONLY: dotp
  USE geq0mod,                         ONLY: geq0
  USE harm,                            ONLY: xmu
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai
  USE system,                          ONLY: cntl,&
                                             cntr,&
                                             ncpw

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rekine

CONTAINS

  ! ==================================================================
  SUBROUTINE rekine(cm,nstate,ekinc)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: cm(ncpw%ngw,nstate)
    REAL(real_8)                             :: ekinc

    INTEGER                                  :: i, ig
    REAL(real_8)                             :: ax, bx, pf

! ==--------------------------------------------------------------==
! ==  COMPUTE FICTITIOUS KINETIC ENERGY OF THE ELECTRONS          ==
! ==--------------------------------------------------------------==

    ekinc=0._real_8
    IF (cntl%tmass) THEN
       DO i=1,nstate
          DO ig=1,ncpw%ngw
             pf=2.0_real_8*xmu(ig)
             ax=REAL(cm(ig,i))
             bx=AIMAG(cm(ig,i))
             ekinc=ekinc+pf*(ax*ax+bx*bx)
          ENDDO
          IF (geq0) ekinc=ekinc-xmu(1)*REAL(cm(1,i))*REAL(cm(1,i))
       ENDDO
    ELSE
       DO i=1,nstate
          ekinc=ekinc+dotp(ncpw%ngw,cm(:,i),cm(:,i))
       ENDDO
       ekinc=ekinc*cntr%emass
    ENDIF
    CALL mp_sum(ekinc,parai%allgrp)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rekine
  ! ==================================================================

END MODULE rekine_utils
