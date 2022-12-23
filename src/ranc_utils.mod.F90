MODULE ranc_utils
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE metr,                            ONLY: metr_com
  USE mp_interface,                    ONLY: mp_bcast
  USE parac,                           ONLY: parai,&
                                             paral
  USE prng_utils,                      ONLY: repprngu
  USE setsc_utils,                     ONLY: ihmat
  USE system,                          ONLY: cntr,&
                                             parm
#include "sizeof.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ranc

CONTAINS

  ! ==================================================================
  SUBROUTINE ranc(tau0)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:)

    INTEGER                                  :: i, ia, is, j
    REAL(real_8)                             :: tx(3)

! Variables
! ==--------------------------------------------------------------==
! ==           Randomization of cell parameters                   ==
! ==--------------------------------------------------------------==

    IF (paral%io_parent)&
         WRITE(6,'(A)') ' RANDOMIZATION - CELL PARAMETERS '
    DO i=1,3
       DO j=1,3
          metr_com%ht(i,j)=metr_com%ht(i,j)+cntr%amprc*(repprngu()-0.5_real_8)
       ENDDO
    ENDDO
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          CALL dgemv('T',3,3,1.0_real_8,metr_com%htm1,3,tau0(1,ia,is),1,0.0_real_8,tx,1)
          CALL dgemv('T',3,3,1.0_real_8,metr_com%ht,3,tx,1,0.0_real_8,tau0(1,ia,is),1)
       ENDDO
    ENDDO
    CALL ihmat(metr_com%ht,metr_com%htm1,parm%omega)
    CALL mp_bcast(metr_com%ht,SIZE(metr_com%ht),parai%source,parai%allgrp)
    CALL mp_bcast(metr_com%htm1,SIZE(metr_com%htm1),parai%source,parai%allgrp)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE ranc
  ! ==================================================================

END MODULE ranc_utils
