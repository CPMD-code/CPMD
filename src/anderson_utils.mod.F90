MODULE anderson_utils
  USE andr,                            ONLY: andr2
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: anderson
  PUBLIC :: anderson_c
  PUBLIC :: change_xmix

CONTAINS

  ! ==================================================================
  SUBROUTINE anderson(xmix,it,rhoinm,rhoutm,rhoin0,rhout0,rhoinp,&
       nnr1,thl)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: xmix
    INTEGER                                  :: it, nnr1
    REAL(real_8)                             :: rhoinp(nnr1), rhout0(nnr1), &
                                                rhoin0(nnr1), rhoutm(nnr1), &
                                                rhoinm(nnr1), thl

    INTEGER                                  :: i
    REAL(real_8)                             :: cthl, cxmix, dr, rl, rl1, sd, &
                                                sn

! Variables
! ==--------------------------------------------------------------==
! == GENERATE NEXT ITERATION USING D. G. ANDERSONS METHOD         ==
! == (FROM HAMANN)                                                ==
! ==--------------------------------------------------------------==

    thl=0._real_8
    IF (it .GT. 1) THEN
       sn=0._real_8
       sd=0._real_8
       !$omp parallel do private(I,RL,RL1,DR) &
       !$omp  reduction(+:SN,SD)
       DO i=1,nnr1
          rl=rhout0(i)-rhoin0(i)
          rl1=rhoutm(i)-rhoinm(i)
          dr=rl-rl1
          sn=sn + rl*dr
          sd=sd + dr*dr
       ENDDO
       CALL mp_sum(sn,parai%allgrp)
       CALL mp_sum(sd,parai%allgrp)
       thl=sn/sd
       IF (ABS(SD).LT.0.5_real_8*ABS(SN)) THEN ! Troubles if |THL| > 2 !
          thl=0._real_8! use simple mixing due to serious convergence probems
       ENDIF
    ENDIF
    cthl =1.0_real_8-thl
    cxmix=1.0_real_8-xmix
    ! ==--------------------------------------------------------------==
    !$omp parallel do private(I) shared(XMIX,CXMIX,THL,CTHL)
    DO i=1,nnr1
       rhoinp(i)=cxmix*(cthl*rhoin0(i)  +&
            thl*rhoinm(i)) +&
            xmix*(cthl*rhout0(i)  +&
            thl*rhoutm(i))
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE anderson
  ! ==================================================================
  SUBROUTINE anderson_c(xmix,it,rhoinm,rhoutm,rhoin0,rhout0,rhoinp,&
       nnr1,thl)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: xmix
    INTEGER                                  :: it, nnr1
    COMPLEX(real_8)                          :: rhoinp(nnr1), rhout0(nnr1), &
                                                rhoin0(nnr1), rhoutm(nnr1), &
                                                rhoinm(nnr1)
    REAL(real_8)                             :: thl

    COMPLEX(real_8)                          :: cthl, cxmix, dr, rl, rl1, sd, &
                                                sn
    INTEGER                                  :: i

! Variables
! ==--------------------------------------------------------------==
! == GENERATE NEXT ITERATION USING D. G. ANDERSONS METHOD         ==
! == (FROM HAMANN)                                                ==
! ==--------------------------------------------------------------==

    thl=0._real_8
    IF (it .GT. 1) THEN
       sn=0._real_8
       sd=0._real_8
       !$omp parallel do private(I,RL,RL1,DR) &
       !$omp  reduction(+:SN,SD)
       DO i=1,nnr1
          rl=rhout0(i)-rhoin0(i)
          rl1=rhoutm(i)-rhoinm(i)
          dr=rl-rl1
          sn=sn + rl*dr
          sd=sd + dr*dr
       ENDDO
       CALL mp_sum(sn,parai%allgrp)
       CALL mp_sum(sd,parai%allgrp)
       thl=sn/sd
       IF (ABS(SD).LT.0.5_real_8*ABS(SN)) THEN ! Troubles if |THL| > 2 !
          thl=0._real_8! use simple mixing due to serious convergence probems
       ENDIF
    ENDIF
    cthl =1.0_real_8-thl
    cxmix=1.0_real_8-xmix
    ! ==--------------------------------------------------------------==
    !$omp parallel do private(I) shared(XMIX,CXMIX,THL,CTHL)
    DO i=1,nnr1
       rhoinp(i)=cxmix*(cthl*rhoin0(i)  +&
            thl*rhoinm(i)) +&
            xmix*(cthl*rhout0(i)  +&
            thl*rhoutm(i))
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE anderson_c
  ! ==================================================================
  SUBROUTINE change_xmix(xmix,drhomax)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: xmix, drhomax

    INTEGER                                  :: i

    IF (andr2%ntabmix.EQ.1) RETURN
    xmix=andr2%andmix(1)
    DO i=1,andr2%ntabmix
       IF (drhomax.LT.andr2%densmix(i)) THEN
          xmix=andr2%andmix(i)
       ENDIF
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE change_xmix
  ! ==================================================================

END MODULE anderson_utils
