MODULE enosmove_utils
  USE kinds,                           ONLY: real_8
  USE nose,                            ONLY: eta,&
                                             etadot,&
                                             nche,&
                                             nchx,&
                                             nedof,&
                                             qnosee
  USE system,                          ONLY: cntr

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: enosmove

CONTAINS

  ! ==================================================================
  SUBROUTINE enosmove(ekinc,step,sctot,ip)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: ekinc, step, sctot
    INTEGER                                  :: ip

    INTEGER                                  :: l
    REAL(real_8)                             :: aa, ckewant, ckine, f1, f2, &
                                                feta(nchx+1)

! Variables
! ==--------------------------------------------------------------==
! ==  COMPUTE ALL THE THERMOSTAT FORCES                           ==
! ==--------------------------------------------------------------==

    ckewant = cntr%ekinw/REAL(nedof,kind=real_8)
    feta(1) = 2._real_8*(ekinc - cntr%ekinw)
    DO l=2,nche+1
       ckine = 0.5_real_8*qnosee(l-1)*etadot(l-1,ip)*etadot(l-1,ip)
       feta(l) = 2._real_8*(ckine - ckewant)
    ENDDO
    ! ==--------------------------------------------------------------==
    ! ==  ETADOT(NCHE) IS TREATED SEPARATELY                          ==
    ! ==--------------------------------------------------------------==
    aa = EXP(-0.125_real_8*step*etadot(nche-1,ip))
    etadot(nche,ip) = etadot(nche,ip)*aa*aa +&
         0.25_real_8*STEP*FETA(NCHE)*AA/QNOSEE(NCHE)
    ckine = 0.5_real_8*qnosee(nche)*etadot(nche,ip)*etadot(nche,ip)
    feta(nche+1) = 2._real_8*(ckine - ckewant)
    f1 = feta(nche-1)
    f2 = feta(nche+1)
    feta(nche-1) = f1 + f2
    ! ==--------------------------------------------------------------==
    ! ==  LOOP BACKWARDS OVER CHAIN                                   ==
    ! ==--------------------------------------------------------------==
    DO l=1,nche-1
       aa = EXP(-0.125_real_8*step*etadot(nche+1-l,ip))
       etadot(nche-l,ip) = etadot(nche-l,ip)*aa*aa +&
            0.25_real_8*STEP*FETA(NCHE-L)*AA/QNOSEE(NCHE-L)
    ENDDO
    ! ==--------------------------------------------------------------==
    ! ==  CALCULATE SCALING FACTOR AND APPLY TO ELECTRON KE           ==
    ! ==  AND ACCUMULATE TOTAL SCALING FACTOR IN SCTOT                ==
    ! ==--------------------------------------------------------------==
    aa = EXP(-0.5_real_8*step*etadot(1,ip))
    sctot = sctot*aa
    ekinc = ekinc*aa*aa
    ! ==--------------------------------------------------------------==
    ! ==  PROPAGATE THERMOSTAT POSITIONS                              ==
    ! ==--------------------------------------------------------------==
    DO l=1,nche
       eta(l,ip) = eta(l,ip) + 0.5_real_8*step*etadot(l,ip)
       feta(l) = 0._real_8
    ENDDO
    ! ==--------------------------------------------------------------==
    ! ==  COMPUTE THERMOSTAT FORCES FOR 1 AND NCHE-1                  ==
    ! ==  (THE REST ARE COMPUTED IN THE PROPAGATION LOOP)             ==
    ! ==--------------------------------------------------------------==
    feta(1) = 2._real_8*(ekinc - cntr%ekinw)
    ckine = 0.5_real_8*qnosee(nche)*etadot(nche,ip)*etadot(nche,ip)
    feta(nche-1) = 2._real_8*(ckine - ckewant)
    ! ==--------------------------------------------------------------==
    ! ==  LOOP FORWARDS OVER CHAIN                                    ==
    ! ==--------------------------------------------------------------==
    DO l=1,nche-1
       aa = EXP(-0.125_real_8*step*etadot(l+1,ip))
       etadot(l,ip) = etadot(l,ip)*aa*aa +&
            0.25_real_8*STEP*FETA(L)*AA/QNOSEE(L)
       ckine = 0.5_real_8*qnosee(l)*etadot(l,ip)*etadot(l,ip)
       feta(l+1) = feta(l+1) + 2._real_8*(ckine - ckewant)
    ENDDO
    aa = EXP(-0.125_real_8*step*etadot(nche-1,ip))
    etadot(nche,ip) = etadot(nche,ip)*aa*aa +&
         0.25_real_8*STEP*FETA(NCHE)*AA/QNOSEE(NCHE)
    ! ==-------------------------------------------------------------==
    RETURN
  END SUBROUTINE enosmove
  ! =================================================================



END MODULE enosmove_utils
