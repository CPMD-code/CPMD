MODULE prcnosmove_utils
  USE cnst,                            ONLY: factem
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE jacobi_utils,                    ONLY: jacobi
  USE kinds,                           ONLY: real_8
  USE metr,                            ONLY: metr_com,&
                                             veps
  USE nose,                            ONLY: etc,&
                                             etcdot,&
                                             gckt,&
                                             nchc,&
                                             nchx,&
                                             qnoscc
  USE system,                          ONLY: cntr
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: prcnosmove_iso
  PUBLIC :: prcnosmove

CONTAINS

  ! ==================================================================
  SUBROUTINE prcnosmove_iso(velp,step,pmx,ip)
    ! ==--------------------------------------------------------------==
    ! ==  ALLOWS FOR NOSE-HOOVER CHAIN DYNAMICS                       ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: velp(:,:,:), step, pmx(*)
    INTEGER                                  :: ip

    INTEGER                                  :: ia, is, k, l
    REAL(real_8)                             :: aa, at3, ekinh, fetc(nchx), &
                                                geps_p, geps_v, hkewant

    at3 = 3._real_8*REAL(ions1%nat,kind=real_8)
    ! ==--------------------------------------------------------------==
    ! == CALCULATE ALL NOSE-HOOVER CHAIN FORCES                       ==
    ! == INCLUDING CELL THERMOSTAT FORCES                             ==
    ! ==--------------------------------------------------------------==
    ! CELL KINETIC ENERGY
    ekinh=0.5_real_8*cntr%cmass*veps*veps
    ! CELL THERMOSTAT FORCES
    hkewant = 0.5_real_8*gckt
    fetc(1) = 2._real_8*(ekinh - hkewant)
    DO l=2,nchc
       ekinh = 0.5_real_8*qnoscc(l-1)*etcdot(l-1,ip)*etcdot(l-1,ip)
       hkewant = 0.5_real_8*cntr%tempc/factem
       fetc(l) = 2._real_8*(ekinh - hkewant)
    ENDDO
    ! ==--------------------------------------------------------------==
    ! == MOVE LAST CHAIN ELEMENTS                                     ==
    ! ==--------------------------------------------------------------==
    etcdot(nchc,ip) = etcdot(nchc,ip) +&
         0.25_real_8*step*fetc(nchc)/qnoscc(nchc)
    ! ==--------------------------------------------------------------==
    ! ==  LOOP OVER REST OF NOSE-HOOVER CHAIN                         ==
    ! ==--------------------------------------------------------------==
    DO l=1,nchc-1
       aa = EXP(-0.125_real_8*step*etcdot(nchc+1-l,ip))
       etcdot(nchc-l,ip) = etcdot(nchc-l,ip)*aa*aa +&
            0.25_real_8*step*fetc(nchc-l)*aa/qnoscc(nchc-l)
    ENDDO
    ! ==----------------------------------------------------------------==
    ! ==  PROPAGATE BAROSTAT VELOCITY                                   ==
    ! ==----------------------------------------------------------------==
    geps_p = metr_com%htfp(1,1) + metr_com%htfp(2,2) + metr_com%htfp(3,3)
    geps_v = 0._real_8
    DO k=1,3
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             geps_v = geps_v + (1._real_8 + 3._real_8/at3)*&
                  pmx(is)*velp(k,ia,is)*velp(k,ia,is)
          ENDDO
       ENDDO
    ENDDO
    aa = EXP(-0.125_real_8*step*etcdot(1,ip))
    veps = veps*aa*aa + 0.25_real_8*step*(geps_p + geps_v)*aa/cntr%cmass
    ! ==----------------------------------------------------------------==
    ! ==  SCALE VELOCITIES                                              ==
    ! ==----------------------------------------------------------------==
    aa = EXP(-0.5_real_8*step*(1._real_8+3._real_8/at3)*veps)
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          velp(1,ia,is) = velp(1,ia,is)*aa
          velp(2,ia,is) = velp(2,ia,is)*aa
          velp(3,ia,is) = velp(3,ia,is)*aa
       ENDDO
    ENDDO
    ! 
    ! ==----------------------------------------------------------------==
    ! ==  PROPAGATE THERMOSTAT POSITIONS                                ==
    ! ==----------------------------------------------------------------==
    DO l=1,nchc
       etc(l,ip) = etc(l,ip) + 0.5_real_8*step*etcdot(l,ip)
    ENDDO
    ! ==----------------------------------------------------------------==
    ! ==  PROPAGATE BAROSTAT VELOCITY                                   ==
    ! ==----------------------------------------------------------------==
    geps_p = metr_com%htfp(1,1) + metr_com%htfp(2,2) + metr_com%htfp(3,3)
    geps_v = 0._real_8
    DO k=1,3
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             geps_v = geps_v + (1._real_8 + 3._real_8/at3)*&
                  pmx(is)*velp(k,ia,is)*velp(k,ia,is)
          ENDDO
       ENDDO
    ENDDO
    aa = EXP(-0.125_real_8*step*etcdot(1,ip))
    veps = veps*aa*aa + 0.25_real_8*step*(geps_p + geps_v)*aa/cntr%cmass
    ! ==--------------------------------------------------------------==
    ! == RECALCULATE ALL NOSE-HOOVER CHAIN FORCES                     ==
    ! == INCLUDING CELL THERMOSTAT FORCES                             ==
    ! ==--------------------------------------------------------------==
    ! CELL KINETIC ENERGY
    ekinh=0.5_real_8*cntr%cmass*veps*veps
    ! CELL THERMOSTAT FORCES
    hkewant = 0.5_real_8*gckt
    fetc(1) = 2._real_8*(ekinh - hkewant)
    ! ==-------------------------------------------------------------==
    ! ==  LOOP OVER CHAIN AGAIN                                      ==
    ! ==-------------------------------------------------------------==
    DO l=1,nchc-1
       aa = EXP(-0.125_real_8*step*etcdot(l+1,ip))
       etcdot(l,ip) = etcdot(l,ip)*aa*aa +&
            0.25_real_8*step*fetc(l)*aa/qnoscc(l)
       ekinh = 0.5_real_8*qnoscc(l)*etcdot(l,ip)*etcdot(l,ip)
       hkewant = 0.5_real_8*cntr%tempc/factem
       fetc(l+1) = 2.0_real_8*(ekinh - hkewant)
    ENDDO
    ! ==--------------------------------------------------------------==
    ! == MOVE LAST CHAIN ELEMENT                                      ==
    ! ==--------------------------------------------------------------==
    etcdot(nchc,ip) = etcdot(nchc,ip) +&
         0.25_real_8*step*fetc(nchc)/qnoscc(nchc)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE prcnosmove_iso
  ! ==================================================================
  ! 
  ! 
  ! ==================================================================
  SUBROUTINE prcnosmove(velp,step,pmx,ip)
    ! ==--------------------------------------------------------------==
    ! ==  ALLOWS FOR NOSE-HOOVER CHAIN DYNAMICS                       ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: velp(:,:,:), step, pmx(*)
    INTEGER                                  :: ip

    INTEGER                                  :: ia, ierr, is, k, l
    REAL(real_8) :: aa, aa1, aa2, aa3, at3, ekinh, fetc(nchx), hkewant, &
      htfv(3,3), scmat(3,3), scmd(3), scmev(3,3), trhv, vtmp1(3)
    REAL(real_8), EXTERNAL                   :: ddot

    at3 = 3._real_8*REAL(ions1%nat,kind=real_8)
    ! ==--------------------------------------------------------------==
    ! == CALCULATE ALL NOSE-HOOVER CHAIN FORCES                       ==
    ! == INCLUDING CELL THERMOSTAT FORCES                             ==
    ! ==--------------------------------------------------------------==
    ! CELL KINETIC ENERGY
    ekinh=0.5_real_8*cntr%cmass*ddot(9,metr_com%htvel,1,metr_com%htvel,1)
    ! CELL THERMOSTAT FORCES
    hkewant = 0.5_real_8*gckt
    fetc(1) = 2._real_8*(ekinh - hkewant)
    DO l=2,nchc
       ekinh = 0.5_real_8*qnoscc(l-1)*etcdot(l-1,ip)*etcdot(l-1,ip)
       hkewant = 0.5_real_8*cntr%tempc/factem
       fetc(l) = 2._real_8*(ekinh - hkewant)
    ENDDO
    ! ==--------------------------------------------------------------==
    ! == MOVE LAST CHAIN ELEMENTS                                     ==
    ! ==--------------------------------------------------------------==
    etcdot(nchc,ip) = etcdot(nchc,ip) +&
         0.25_real_8*step*fetc(nchc)/qnoscc(nchc)
    ! ==--------------------------------------------------------------==
    ! ==  LOOP OVER REST OF NOSE-HOOVER CHAIN                         ==
    ! ==--------------------------------------------------------------==
    DO l=1,nchc-1
       aa = EXP(-0.125_real_8*step*etcdot(nchc+1-l,ip))
       etcdot(nchc-l,ip) = etcdot(nchc-l,ip)*aa*aa +&
            0.25_real_8*step*fetc(nchc-l)*aa/qnoscc(nchc-l)
    ENDDO
    ! ==--------------------------------------------------------------==
    ! ==  GET VELOCITY-DEPENDENT PART OF CELL FORCE                   ==
    ! ==--------------------------------------------------------------==
    CALL zeroing(htfv)!,9)
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          DO k=1,3
             DO l=1,3
                htfv(k,l) = htfv(k,l) + pmx(is)*velp(k,ia,is)*velp(l,ia,is)
             ENDDO
             htfv(k,k) = htfv(k,k) +&
                  pmx(is)*velp(k,ia,is)*velp(k,ia,is)/at3
          ENDDO
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    ! ==  PROPAGATE BAROSTAT (CELL VELOCITY)                          ==
    ! ==--------------------------------------------------------------==
    aa = EXP(-0.125_real_8*step*etcdot(1,ip))
    DO k=1,3
       DO l=1,3
          metr_com%htvel(k,l) = metr_com%htvel(k,l)*aa*aa +&
               0.25_real_8*step*(htfv(k,l) + metr_com%htfp(k,l))*aa/cntr%cmass
       ENDDO
    ENDDO
    ! ==----------------------------------------------------------------==
    ! ==  FOR VELOCITY SCALING, FORM MATRIX v_h + Tr(v_h)/N + v_zeta1   ==
    ! ==----------------------------------------------------------------==
    trhv = metr_com%htvel(1,1) + metr_com%htvel(2,2) + metr_com%htvel(3,3)
    DO k=1,3
       DO l=1,3
          scmat(k,l) = metr_com%htvel(k,l)
       ENDDO
       scmat(k,k) = scmat(k,k) + trhv/at3
    ENDDO
    ! ==----------------------------------------------------------------==
    ! ==  DIAGONALIZE THIS MATRIX AND SAVE EIGENVALUES AND EIGENVECTORS ==
    ! ==----------------------------------------------------------------==
    CALL jacobi(3,3,scmat,scmd,scmev,ierr)
    ! ==----------------------------------------------------------------==
    ! ==  ROTATE, SCALE, AND ROTATE BACK FOR VELOCITIES                 ==
    ! ==----------------------------------------------------------------==
    aa1 = EXP(-0.5_real_8*step*scmd(1))
    aa2 = EXP(-0.5_real_8*step*scmd(2))
    aa3 = EXP(-0.5_real_8*step*scmd(3))
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          ! 
          vtmp1(1) = scmev(1,1)*velp(1,ia,is) + scmev(2,1)*velp(2,ia,is)&
               + scmev(3,1)*velp(3,ia,is)
          vtmp1(2) = scmev(1,2)*velp(1,ia,is) + scmev(2,2)*velp(2,ia,is)&
               + scmev(3,2)*velp(3,ia,is)
          vtmp1(3) = scmev(1,3)*velp(1,ia,is) + scmev(2,3)*velp(2,ia,is)&
               + scmev(3,3)*velp(3,ia,is)
          ! 
          vtmp1(1) = vtmp1(1)*aa1
          vtmp1(2) = vtmp1(2)*aa2
          vtmp1(3) = vtmp1(3)*aa3
          ! 
          velp(1,ia,is) = scmev(1,1)*vtmp1(1) + scmev(1,2)*vtmp1(2)&
               + scmev(1,3)*vtmp1(3)
          velp(2,ia,is) = scmev(2,1)*vtmp1(1) + scmev(2,2)*vtmp1(2)&
               + scmev(2,3)*vtmp1(3)
          velp(3,ia,is) = scmev(3,1)*vtmp1(1) + scmev(3,2)*vtmp1(2)&
               + scmev(3,3)*vtmp1(3)
       ENDDO
    ENDDO
    ! 
    ! ==----------------------------------------------------------------==
    ! ==  PROPAGATE THERMOSTAT POSITIONS                                ==
    ! ==----------------------------------------------------------------==
    DO l=1,nchc
       etc(l,ip) = etc(l,ip) + 0.5_real_8*step*etcdot(l,ip)
    ENDDO
    ! ==--------------------------------------------------------------==
    ! ==  GET VELOCITY-DEPENDENT PART OF CELL FORCE                   ==
    ! ==--------------------------------------------------------------==
    CALL zeroing(htfv)!,9)
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          DO k=1,3
             DO l=1,3
                htfv(k,l) = htfv(k,l) + pmx(is)*velp(k,ia,is)*velp(l,ia,is)
             ENDDO
             htfv(k,k) = htfv(k,k) +&
                  pmx(is)*velp(k,ia,is)*velp(k,ia,is)/at3
          ENDDO
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    ! ==  PROPAGATE BAROSTAT (CELL VELOCITY)                          ==
    ! ==--------------------------------------------------------------==
    aa = EXP(-0.125_real_8*step*etcdot(1,ip))
    DO k=1,3
       DO l=1,3
          metr_com%htvel(k,l) = metr_com%htvel(k,l)*aa*aa +&
               0.25_real_8*step*(htfv(k,l) + metr_com%htfp(k,l))*aa/cntr%cmass
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    ! == RECALCULATE ALL NOSE-HOOVER CHAIN FORCES                     ==
    ! == INCLUDING CELL THERMOSTAT FORCES                             ==
    ! ==--------------------------------------------------------------==
    ! CELL KINETIC ENERGY
    ekinh=0.5_real_8*cntr%cmass*ddot(9,metr_com%htvel,1,metr_com%htvel,1)
    ! CELL THERMOSTAT FORCES
    hkewant = 0.5_real_8*gckt
    fetc(1) = 2._real_8*(ekinh - hkewant)
    ! ==-------------------------------------------------------------==
    ! ==  LOOP OVER CHAIN AGAIN                                      ==
    ! ==-------------------------------------------------------------==
    DO l=1,nchc-1
       aa = EXP(-0.125_real_8*step*etcdot(l+1,ip))
       etcdot(l,ip) = etcdot(l,ip)*aa*aa +&
            0.25_real_8*step*fetc(l)*aa/qnoscc(l)
       ekinh = 0.5_real_8*qnoscc(l)*etcdot(l,ip)*etcdot(l,ip)
       hkewant = 0.5_real_8*cntr%tempc/factem
       fetc(l+1) = 2.0_real_8*(ekinh - hkewant)
    ENDDO
    ! ==--------------------------------------------------------------==
    ! == MOVE LAST CHAIN ELEMENT                                      ==
    ! ==--------------------------------------------------------------==
    etcdot(nchc,ip) = etcdot(nchc,ip) +&
         0.25_real_8*step*fetc(nchc)/qnoscc(nchc)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE prcnosmove
  ! ==================================================================

END MODULE prcnosmove_utils
