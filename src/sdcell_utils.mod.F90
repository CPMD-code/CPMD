MODULE sdcell_utils
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE metr,                            ONLY: metr_com
  USE parac,                           ONLY: parai
  USE prcp,                            ONLY: prcpl
  USE setsc_utils,                     ONLY: ihmat
  USE system,                          ONLY: cntr,&
                                             parm
  USE tpar,                            ONLY: dt_ions
#include "sizeof.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: sdcell

CONTAINS

  ! ==================================================================
  SUBROUTINE sdcell(tau0)
    ! ==--------------------------------------------------------------==
    ! ==  UPDATE CELL WITH STEEPEST DESCENT                           ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:)

    INTEGER                                  :: i, ia, is, j
    REAL(real_8)                             :: cfor, diag, fac, pp, tx(3)

! Variables
! ==--------------------------------------------------------------==

    fac=dt_ions*dt_ions/cntr%cmass
    IF (prcpl%tzflex) THEN
       ! z-direction only flexible cell
       ! We use only the z-direction (pressure).
       cfor = metr_com%htfor(3,3)
       diag = metr_com%ht(3,3)
       IF (diag.EQ.0._real_8) THEN
          diag = 1._real_8
       ENDIF
       pp = 1._real_8 + fac*cfor/diag
       DO j=1,3
          metr_com%ht(3,j) = pp*metr_com%ht(3,j)
       ENDDO
    ELSE IF (prcpl%tisot) THEN
       ! ISOTROPIC CELL
       ! We use only the diagonal (pressure).
       cfor = metr_com%htfor(1,1)+metr_com%htfor(2,2)+metr_com%htfor(3,3)
       diag = metr_com%ht(1,1)+metr_com%ht(2,2)+metr_com%ht(3,3)
       IF (diag.EQ.0._real_8) THEN
          diag = 1._real_8
       ENDIF
       pp = 1._real_8 + fac*cfor/diag
       DO i=1,3
          DO j=1,3
             metr_com%ht(i,j) = pp*metr_com%ht(i,j)
          ENDDO
       ENDDO
    ELSE
       ! NO CONSTRAINT
       DO i=1,3
          DO j=1,3
             metr_com%ht(i,j)=metr_com%ht(i,j)+fac*metr_com%htfor(i,j)
          ENDDO
       ENDDO
    ENDIF
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          CALL dgemv('T',3,3,1.0_real_8,metr_com%htm1,3,tau0(1,ia,is),1,0.0_real_8,tx,1)
          CALL dgemv('T',3,3,1.0_real_8,metr_com%ht,3,tx,1,0.0_real_8,tau0(1,ia,is),1)
       ENDDO
    ENDDO
    CALL mp_bcast_byte(metr_com, size_in_bytes_of(metr_com),parai%source,parai%allgrp)
    CALL ihmat(metr_com%ht,metr_com%htm1,parm%omega)
    DO i=1,3
       parm%a1(i) = metr_com%ht(1,i)
       parm%a2(i) = metr_com%ht(2,i)
       parm%a3(i) = metr_com%ht(3,i)
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE sdcell
  ! ==================================================================

END MODULE sdcell_utils
