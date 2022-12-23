MODULE prpcnosmove_utils
  USE cnst,                            ONLY: factem
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE jacobi_utils,                    ONLY: jacobi
  USE kinds,                           ONLY: real_8
  USE metr,                            ONLY: metr_com,&
                                             veps
  USE nose,                            ONLY: &
       etap, etap1, etap1dot, etapdot, etapm, etapmdot, etc, etcdot, gckt, &
       gkt, mapdof, nchain, nchc, nchx, ndfnt, nosl, ntherm, qnoscc, qnosp, &
       qnospc, qnospm
  USE system,                          ONLY: cntr
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: prpcnosmove_iso
  PUBLIC :: prpcnosmove

CONTAINS

  ! ==================================================================
  SUBROUTINE prpcnosmove_iso(velp,step,pmx,ip)
    ! ==--------------------------------------------------------------==
    ! ==  ALLOWS FOR NOSE-HOOVER CHAIN DYNAMICS                       ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: velp(:,:,:), step, pmx(*)
    INTEGER                                  :: ip

    INTEGER                                  :: ia, idf, ij, ijj, ipp, is, k, &
                                                kf, l
    REAL(real_8)                             :: aa, at3, const, ekinh, ekinp, &
                                                fetapv(nchx), fetc(nchx), &
                                                geps_p, geps_v, hkewant, &
                                                pkewant

    ipp=MIN(ip,2)
    at3 = 3._real_8*REAL(ions1%nat,kind=real_8)
    ! ==--------------------------------------------------------------==
    ! == SEPARATE THERMOSTATS ON DIFFERENT ATOMIC SPECIES             ==
    ! ==--------------------------------------------------------------==
    IF (nosl%tultra) THEN
       DO is=1,ions1%nsp
          ! ==--------------------------------------------------------------==
          ! ==  CALCULATE ALL THERMOSTAT FORCES                             ==
          ! ==--------------------------------------------------------------==
          ekinp=0._real_8
          const=0.5_real_8*pmx(is)
          DO k=1,3
             DO ia=1,ions0%na(is)
                ekinp=ekinp+const*velp(k,ia,is)*velp(k,ia,is)
             ENDDO
          ENDDO
          pkewant = 0.5_real_8*gkt(ipp,is)
          fetapv(1) = 2.0_real_8*(ekinp - pkewant)
          DO l=2,nchain
             ekinp = 0.5_real_8*qnosp(is,l-1,ip)*etapdot(is,l-1,ip)*&
                  etapdot(is,l-1,ip)
             pkewant = 0.5_real_8*cntr%tempw/factem
             fetapv(l) = 2.0_real_8*(ekinp - pkewant)
          ENDDO
          ! ==----------------------------------------------------------------==
          ! ==  PROPAGATE LAST CHAIN ELEMENT                                  ==
          ! ==----------------------------------------------------------------==
          etapdot(is,nchain,ip) = etapdot(is,nchain,ip) +&
               0.25_real_8*step*fetapv(nchain)/qnosp(is,nchain,ip)
          ! ==----------------------------------------------------------------==
          ! ==  PROPAGATE REMAINDER OF CHAIN                                  ==
          ! ==----------------------------------------------------------------==
          DO l=1,nchain-1
             aa = EXP(-0.125_real_8*step*etapdot(is,nchain+1-l,ip))
             etapdot(is,nchain-l,ip) = etapdot(is,nchain-l,ip)*aa*aa +&
                  0.25_real_8*step*fetapv(nchain-l)*aa/qnosp(is,nchain-l,ip)
          ENDDO
       ENDDO
       ! ==--------------------------------------------------------------==
       ! == CALCULATE CELL THERMOSTAT FORCES                             ==
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
       DO is=1,ions1%nsp
          aa = EXP(-0.5_real_8*step*((1._real_8+3._real_8/at3)*veps + etapdot(is,1,ip)))
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
       DO is=1,ions1%nsp
          DO l=1,nchain
             etap(is,l,ip) = etap(is,l,ip) + 0.5_real_8*step*etapdot(is,l,ip)
          ENDDO
       ENDDO
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
       ! ==----------------------------------------------------------------==
       ! == LOOP AGAIN OVER CHAIN                                          ==
       ! ==----------------------------------------------------------------==
       DO is=1,ions1%nsp
          ekinp=0._real_8
          const=0.5_real_8*pmx(is)
          DO k=1,3
             DO ia=1,ions0%na(is)
                ekinp=ekinp+const*velp(k,ia,is)*velp(k,ia,is)
             ENDDO
          ENDDO
          pkewant = 0.5_real_8*gkt(ipp,is)
          fetapv(1) = 2.0_real_8*(ekinp - pkewant)
          DO l=1,nchain-1
             aa = EXP(-0.125_real_8*step*etapdot(is,l+1,ip))
             etapdot(is,l,ip) = etapdot(is,l,ip)*aa*aa +&
                  0.25_real_8*step*fetapv(l)*aa/qnosp(is,l,ip)
             ekinp = 0.5_real_8*qnosp(is,l,ip)*etapdot(is,l,ip)*&
                  etapdot(is,l,ip)
             pkewant = 0.5_real_8*cntr%tempw/factem
             fetapv(l+1) = 2.0_real_8*(ekinp - pkewant)
          ENDDO
          ! ==----------------------------------------------------------------==
          ! ==  PROPAGATE LAST CHAIN ELEMENT                                  ==
          ! ==----------------------------------------------------------------==
          etapdot(is,nchain,ip) = etapdot(is,nchain,ip) +&
               0.25_real_8*step*fetapv(nchain)/qnosp(is,nchain,ip)
       ENDDO
       ! ==--------------------------------------------------------------==
       ! == SEPARATE THERMOSTATS ON EVERY DEGREE OF FREEDOM              ==
       ! ==--------------------------------------------------------------==
    ELSEIF (nosl%tmnose) THEN
       ! ==--------------------------------------------------------------==
       ! ==  LOOP BACKWARDS OVER NOSE-HOOVER CHAIN                       ==
       ! ==--------------------------------------------------------------==
       ij=0
       DO kf=1,ntherm(ipp)
          ! ==--------------------------------------------------------------==
          ! == CALCULATE ALL THERMOSTAT FORCES                              ==
          ! ==--------------------------------------------------------------==
          ekinp=0._real_8
          DO idf=1,ndfnt(kf,ipp)
             ijj=ij+idf
             k=mapdof(1,ijj,ipp)
             ia=mapdof(2,ijj,ipp)
             is=mapdof(3,ijj,ipp)
             const=0.5_real_8*pmx(is)
             ekinp=ekinp+const*velp(k,ia,is)*velp(k,ia,is)
          ENDDO
          pkewant = 0.5_real_8*gkt(ipp,kf)
          fetapv(1) = 2.0_real_8*(ekinp - pkewant)
          DO l=2,nchain
             ekinp = 0.5_real_8*qnospm(kf,l-1,ip)*&
                  etapmdot(kf,l-1,ip)*etapmdot(kf,l-1,ip)
             pkewant = 0.5_real_8*cntr%tempw/factem
             fetapv(l) = 2.0_real_8*(ekinp - pkewant)
          ENDDO
          ! ==--------------------------------------------------------------==
          ! == PROPAGATE LAST CHAIN ELEMENT                                 ==
          ! ==--------------------------------------------------------------==
          etapmdot(kf,nchain,ip)=etapmdot(kf,nchain,ip)+&
               0.25_real_8*step*fetapv(nchain)/qnospm(kf,nchain,ip)
          ! ==--------------------------------------------------------------==
          ! == PROPAGATE REST OF CHAIN                                      ==
          ! ==--------------------------------------------------------------==
          DO l=1,nchain-1
             aa = EXP(-0.125_real_8*step*&
                  etapmdot(kf,nchain+1-l,ip))
             etapmdot(kf,nchain-l,ip) &
                  =etapmdot(kf,nchain-l,ip)&
                  *aa*aa+0.25_real_8*step*fetapv(nchain-l)&
                  *aa/qnospm(kf,nchain-l,ip)
          ENDDO
          ij=ij+ndfnt(kf,ipp)
       ENDDO
       ! ==--------------------------------------------------------------==
       ! == CALCULATE CELL THERMOSTAT FORCES                             ==
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
       ij=0
       DO kf=1,ntherm(ipp)
          aa=EXP(-0.5_real_8*step*((1._real_8+3._real_8/at3)*veps+etapmdot(kf,1,ip)))
          DO idf=1,ndfnt(kf,ipp)
             ijj=ij+idf
             k=mapdof(1,ijj,ipp)
             ia=mapdof(2,ijj,ipp)
             is=mapdof(3,ijj,ipp)
             velp(k,ia,is) = velp(k,ia,is)*aa
          ENDDO
          ij=ij+ndfnt(kf,ipp)
       ENDDO
       ! ==----------------------------------------------------------------==
       ! == PROPAGATE THERMOSTAT POSITIONS                                 ==
       ! ==----------------------------------------------------------------==
       DO kf=1,ntherm(ipp)
          DO l=1,nchain
             etapm(kf,l,ip) = etapm(kf,l,ip) +&
                  0.5_real_8*step*etapmdot(kf,l,ip)
          ENDDO
       ENDDO
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
       ! ==-------------------------------------------------------------==
       ! ==  LOOP OVER CHAIN AGAIN                                      ==
       ! ==-------------------------------------------------------------==
       ij=0
       DO kf=1,ntherm(ipp)
          const=0.5_real_8*pmx(is)
          ekinp=0._real_8
          DO idf=1,ndfnt(kf,ipp)
             ijj=ij+idf
             k=mapdof(1,ijj,ipp)
             ia=mapdof(2,ijj,ipp)
             is=mapdof(3,ijj,ipp)
             ekinp=ekinp+const*velp(k,ia,is)*velp(k,ia,is)
          ENDDO
          pkewant = 0.5_real_8*gkt(ipp,kf)
          fetapv(1) = 2.0_real_8*(ekinp - pkewant)
          DO l=1,nchain-1
             aa = EXP(-0.125_real_8*step*etapmdot(kf,l+1,ip))
             etapmdot(kf,l,ip) = etapmdot(kf,l,ip)*aa*aa +&
                  0.25_real_8*step*fetapv(l)*aa/qnospm(kf,l,ip)
             ekinp = 0.5_real_8*qnospm(kf,l,ip)*&
                  etapmdot(kf,l,ip)*etapmdot(kf,l,ip)
             pkewant = 0.5_real_8*cntr%tempw/factem
             fetapv(l+1) = 2.0_real_8*(ekinp - pkewant)
          ENDDO
          ! ==--------------------------------------------------------------==
          ! == PROPAGATE LAST CHAIN ELEMENT                                 ==
          ! ==--------------------------------------------------------------==
          etapmdot(kf,nchain,ip)=etapmdot(kf,nchain,ip)+&
               0.25_real_8*step*fetapv(nchain)/qnospm(kf,nchain,ip)
          ! ==-------------------------------------------------------------==
          ij=ij+ndfnt(kf,ipp)
       ENDDO
       ! ==--------------------------------------------------------------==
       ! == GLOBAL THERMOSTAT ON THE ATOMS                               ==
       ! ==--------------------------------------------------------------==
    ELSE
       ! ==--------------------------------------------------------------==
       ! == CALCULATE ALL NOSE-HOOVER CHAIN FORCES                       ==
       ! == INCLUDING CELL THERMOSTAT FORCES                             ==
       ! ==--------------------------------------------------------------==
       ! PARTICLE KINETIC ENERGY
       ekinp=0._real_8
       DO k=1,3
          DO is=1,ions1%nsp
             const=0.5_real_8*pmx(is)
             DO ia=1,ions0%na(is)
                ekinp=ekinp+const*velp(k,ia,is)*velp(k,ia,is)
             ENDDO
          ENDDO
       ENDDO
       ! PARTICLE THERMOSTAT FORCES
       pkewant = 0.5_real_8*gkt(ipp,1)
       fetapv(1) = 2.0_real_8*(ekinp - pkewant)
       DO l=2,nchain
          ekinp = 0.5_real_8*qnospc(l-1,ip)*etap1dot(l-1,ip)*&
               etap1dot(l-1,ip)
          pkewant = 0.5_real_8*cntr%tempw/factem
          fetapv(l) = 2.0_real_8*(ekinp - pkewant)
       ENDDO
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
       etap1dot(nchain,ip) = etap1dot(nchain,ip) +&
            0.25_real_8*step*fetapv(nchain)/qnospc(nchain,ip)
       etcdot(nchc,ip) = etcdot(nchc,ip) +&
            0.25_real_8*step*fetc(nchc)/qnoscc(nchc)
       ! ==--------------------------------------------------------------==
       ! ==  LOOP OVER REST OF NOSE-HOOVER CHAIN                         ==
       ! ==--------------------------------------------------------------==
       DO l=1,nchain-1
          aa = EXP(-0.125_real_8*step*etap1dot(nchain+1-l,ip))
          etap1dot(nchain-l,ip) = etap1dot(nchain-l,ip)*aa*aa +&
               0.25_real_8*step*fetapv(nchain-l)*aa/qnospc(nchain-l,ip)
       ENDDO
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
       aa = EXP(-0.5_real_8*step*((1._real_8+3._real_8/at3)*veps + etap1dot(1,ip)))
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
       DO l=1,nchain
          etap1(l,ip) = etap1(l,ip) + 0.5_real_8*step*etap1dot(l,ip)
       ENDDO
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
       ! PARTICLE KINETIC ENERGY
       ekinp=0._real_8
       DO k=1,3
          DO is=1,ions1%nsp
             const=0.5_real_8*pmx(is)
             DO ia=1,ions0%na(is)
                ekinp=ekinp+const*velp(k,ia,is)*velp(k,ia,is)
             ENDDO
          ENDDO
       ENDDO
       ! PARTICLE THERMOSTAT FORCES
       pkewant = 0.5_real_8*gkt(ipp,1)
       fetapv(1) = 2.0_real_8*(ekinp - pkewant)
       ! CELL KINETIC ENERGY
       ekinh=0.5_real_8*cntr%cmass*veps*veps
       ! CELL THERMOSTAT FORCES
       hkewant = 0.5_real_8*gckt
       fetc(1) = 2._real_8*(ekinh - hkewant)
       ! ==-------------------------------------------------------------==
       ! ==  LOOP OVER CHAIN AGAIN                                      ==
       ! ==-------------------------------------------------------------==
       DO l=1,nchain-1
          aa = EXP(-0.125_real_8*step*etap1dot(l+1,ip))
          etap1dot(l,ip) = etap1dot(l,ip)*aa*aa +&
               0.25_real_8*step*fetapv(l)*aa/qnospc(l,ip)
          ekinp = 0.5_real_8*qnospc(l,ip)*etap1dot(l,ip)*etap1dot(l,ip)
          pkewant = 0.5_real_8*cntr%tempw/factem
          fetapv(l+1) = 2.0_real_8*(ekinp - pkewant)
       ENDDO
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
       etap1dot(nchain,ip) = etap1dot(nchain,ip) +&
            0.25_real_8*step*fetapv(nchain)/qnospc(nchain,ip)
       etcdot(nchc,ip) = etcdot(nchc,ip) +&
            0.25_real_8*step*fetc(nchc)/qnoscc(nchc)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE prpcnosmove_iso
  ! ==================================================================
  ! 
  ! 
  ! ==================================================================
  SUBROUTINE prpcnosmove(velp,step,pmx,ip)
    ! ==--------------------------------------------------------------==
    ! ==  ALLOWS FOR NOSE-HOOVER CHAIN DYNAMICS                       ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: velp(:,:,:), step, pmx(*)
    INTEGER                                  :: ip

    INTEGER                                  :: ia, idf, ierr, ij, ijj, ipp, &
                                                is, k, kf, l
    REAL(real_8) :: aa, aa1, aa2, aa3, at3, const, ekinh, ekinp, &
      fetapv(nchx), fetc(nchx), hkewant, htfv(3,3), pkewant, scmat(3,3), &
      scmd(3), scmev(3,3), trhv, vtmp1(3)
    REAL(real_8), EXTERNAL                   :: ddot

    ipp=MIN(ip,2)
    at3 = 3._real_8*REAL(ions1%nat,kind=real_8)
    ! ==--------------------------------------------------------------==
    ! == SEPARATE THERMOSTATS ON DIFFERENT ATOMIC SPECIES             ==
    ! ==--------------------------------------------------------------==
    IF (nosl%tultra) THEN
       DO is=1,ions1%nsp
          ! ==--------------------------------------------------------------==
          ! ==  CALCULATE ALL THERMOSTAT FORCES                             ==
          ! ==--------------------------------------------------------------==
          ekinp=0._real_8
          const=0.5_real_8*pmx(is)
          DO k=1,3
             DO ia=1,ions0%na(is)
                ekinp=ekinp+const*velp(k,ia,is)*velp(k,ia,is)
             ENDDO
          ENDDO
          pkewant = 0.5_real_8*gkt(ipp,is)
          fetapv(1) = 2.0_real_8*(ekinp - pkewant)
          DO l=2,nchain
             ekinp = 0.5_real_8*qnosp(is,l-1,ip)*etapdot(is,l-1,ip)*&
                  etapdot(is,l-1,ip)
             pkewant = 0.5_real_8*cntr%tempw/factem
             fetapv(l) = 2.0_real_8*(ekinp - pkewant)
          ENDDO
          ! ==----------------------------------------------------------------==
          ! ==  PROPAGATE LAST CHAIN ELEMENT                                  ==
          ! ==----------------------------------------------------------------==
          etapdot(is,nchain,ip) = etapdot(is,nchain,ip) +&
               0.25_real_8*step*fetapv(nchain)/qnosp(is,nchain,ip)
          ! ==----------------------------------------------------------------==
          ! ==  PROPAGATE REMAINDER OF CHAIN                                  ==
          ! ==----------------------------------------------------------------==
          DO l=1,nchain-1
             aa = EXP(-0.125_real_8*step*etapdot(is,nchain+1-l,ip))
             etapdot(is,nchain-l,ip) = etapdot(is,nchain-l,ip)*aa*aa +&
                  0.25_real_8*step*fetapv(nchain-l)*aa/qnosp(is,nchain-l,ip)
          ENDDO
       ENDDO
       ! ==--------------------------------------------------------------==
       ! == CALCULATE CELL THERMOSTAT FORCES                             ==
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
       ! ==  FOR VELOCITY SCALING, FORM MATRIX v_h + Tr(v_h)/N + zeta_1    ==
       ! ==----------------------------------------------------------------==
       DO is=1,ions1%nsp
          trhv = metr_com%htvel(1,1) + metr_com%htvel(2,2) + metr_com%htvel(3,3)
          DO k=1,3
             DO l=1,3
                scmat(k,l) = metr_com%htvel(k,l)
             ENDDO
             scmat(k,k) = scmat(k,k) + trhv/at3 + etapdot(is,1,ip)
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
          ! 
          DO ia=1,ions0%na(is)
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
       DO is=1,ions1%nsp
          DO l=1,nchain
             etap(is,l,ip) = etap(is,l,ip) + 0.5_real_8*step*etapdot(is,l,ip)
          ENDDO
       ENDDO
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
       ! ==----------------------------------------------------------------==
       ! == LOOP AGAIN OVER CHAIN                                          ==
       ! ==----------------------------------------------------------------==
       DO is=1,ions1%nsp
          ekinp=0._real_8
          const=0.5_real_8*pmx(is)
          DO k=1,3
             DO ia=1,ions0%na(is)
                ekinp=ekinp+const*velp(k,ia,is)*velp(k,ia,is)
             ENDDO
          ENDDO
          pkewant = 0.5_real_8*gkt(ipp,is)
          fetapv(1) = 2.0_real_8*(ekinp - pkewant)
          DO l=1,nchain-1
             aa = EXP(-0.125_real_8*step*etapdot(is,l+1,ip))
             etapdot(is,l,ip) = etapdot(is,l,ip)*aa*aa +&
                  0.25_real_8*step*fetapv(l)*aa/qnosp(is,l,ip)
             ekinp = 0.5_real_8*qnosp(is,l,ip)*etapdot(is,l,ip)*&
                  etapdot(is,l,ip)
             pkewant = 0.5_real_8*cntr%tempw/factem
             fetapv(l+1) = 2.0_real_8*(ekinp - pkewant)
          ENDDO
          ! ==----------------------------------------------------------------==
          ! ==  PROPAGATE LAST CHAIN ELEMENT                                  ==
          ! ==----------------------------------------------------------------==
          etapdot(is,nchain,ip) = etapdot(is,nchain,ip) +&
               0.25_real_8*step*fetapv(nchain)/qnosp(is,nchain,ip)
       ENDDO
       ! ==--------------------------------------------------------------==
       ! == SEPARATE THERMOSTATS ON EVERY DEGREE OF FREEDOM              ==
       ! ==--------------------------------------------------------------==
    ELSEIF (nosl%tmnose) THEN
       ! ==--------------------------------------------------------------==
       ! ==  LOOP BACKWARDS OVER NOSE-HOOVER CHAIN                       ==
       ! ==--------------------------------------------------------------==
       ij=0
       DO kf=1,ntherm(ipp)
          ! ==--------------------------------------------------------------==
          ! == CALCULATE ALL THERMOSTAT FORCES                              ==
          ! ==--------------------------------------------------------------==
          ekinp=0._real_8
          DO idf=1,ndfnt(kf,ipp)
             ijj=ij+idf
             k=mapdof(1,ijj,ipp)
             ia=mapdof(2,ijj,ipp)
             is=mapdof(3,ijj,ipp)
             const=0.5_real_8*pmx(is)
             ekinp=ekinp+const*velp(k,ia,is)*velp(k,ia,is)
          ENDDO
          pkewant = 0.5_real_8*gkt(ipp,kf)
          fetapv(1) = 2.0_real_8*(ekinp - pkewant)
          DO l=2,nchain
             ekinp = 0.5_real_8*qnospm(kf,l-1,ip)*&
                  etapmdot(kf,l-1,ip)*etapmdot(kf,l-1,ip)
             pkewant = 0.5_real_8*cntr%tempw/factem
             fetapv(l) = 2.0_real_8*(ekinp - pkewant)
          ENDDO
          ! ==--------------------------------------------------------------==
          ! == PROPAGATE LAST CHAIN ELEMENT                                 ==
          ! ==--------------------------------------------------------------==
          etapmdot(kf,nchain,ip)=etapmdot(kf,nchain,ip)+&
               0.25_real_8*step*fetapv(nchain)/qnospm(kf,nchain,ip)
          ! ==--------------------------------------------------------------==
          ! == PROPAGATE REST OF CHAIN                                      ==
          ! ==--------------------------------------------------------------==
          DO l=1,nchain-1
             aa = EXP(-0.125_real_8*step*etapmdot(kf,nchain+1-l,ip))
             etapmdot(kf,nchain-l,ip) = etapmdot(kf,nchain-l,ip)&
                  *aa*aa+0.25_real_8*step*fetapv(nchain-l)&
                  *aa/qnospm(kf,nchain-l,ip)
          ENDDO
          ij=ij+ndfnt(kf,ipp)
       ENDDO
       ! ==--------------------------------------------------------------==
       ! == CALCULATE CELL THERMOSTAT FORCES                             ==
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
       ! == DO VELOCITY SCALING AND PROPAGATE THERMOSTAT POSITIONS         ==
       ! ==----------------------------------------------------------------==
       kf=0
       DO is=1,ions1%nsp
          DO ia=1,ions0%na( is )
             ! ==----------------------------------------------------------------==
             ! ==  FOR VELOCITY SCALING, FORM MATRIX v_h + Tr(v_h)/N + zeta_1    ==
             ! ==----------------------------------------------------------------==
             trhv = metr_com%htvel(1,1) + metr_com%htvel(2,2) + metr_com%htvel(3,3)
             DO k=1,3
                kf=kf+1
                DO l=1,3
                   scmat(k,l) = metr_com%htvel(k,l)
                ENDDO
                scmat(k,k) = scmat(k,k) + trhv/at3 + etapmdot(kf,1,ip)
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
       ! == PROPAGATE THERMOSTAT POSITIONS                                 ==
       ! ==----------------------------------------------------------------==
       DO kf=1,ntherm(ipp)
          DO l=1,nchain
             etapm(kf,l,ip) = etapm(kf,l,ip) +&
                  0.5_real_8*step*etapmdot(kf,l,ip)
          ENDDO
       ENDDO
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
       ! ==-------------------------------------------------------------==
       ! ==  LOOP OVER CHAIN AGAIN                                      ==
       ! ==-------------------------------------------------------------==
       ij=0
       DO kf=1,ntherm(ipp)
          const=0.5_real_8*pmx(is)
          ekinp=0._real_8
          DO idf=1,ndfnt(kf,ipp)
             ijj=ij+idf
             k=mapdof(1,ijj,ipp)
             ia=mapdof(2,ijj,ipp)
             is=mapdof(3,ijj,ipp)
             ekinp=ekinp+const*velp(k,ia,is)*velp(k,ia,is)
          ENDDO
          pkewant = 0.5_real_8*gkt(ipp,kf)
          fetapv(1) = 2.0_real_8*(ekinp - pkewant)
          DO l=1,nchain-1
             aa = EXP(-0.125_real_8*step*etapmdot(kf,l+1,ip))
             etapmdot(kf,l,ip) = etapmdot(kf,l,ip)*&
                  aa*aa + 0.25_real_8*step*fetapv(l)*aa/qnospm(kf,l,ip)
             ekinp = 0.5_real_8*qnospm(kf,l,ip)*&
                  etapmdot(kf,l,ip)*etapmdot(kf,l,ip)
             pkewant = 0.5_real_8*cntr%tempw/factem
             fetapv(l+1) = 2.0_real_8*(ekinp - pkewant)
          ENDDO
          ! ==--------------------------------------------------------------==
          ! == PROPAGATE LAST CHAIN ELEMENT                                 ==
          ! ==--------------------------------------------------------------==
          etapmdot(kf,nchain,ip)=etapmdot(kf,nchain,ip)+&
               0.25_real_8*step*fetapv(nchain)/qnospm(kf,nchain,ip)
          ! ==-------------------------------------------------------------==
          ij=ij+ndfnt(kf,ipp)
       ENDDO
       ! ==--------------------------------------------------------------==
       ! == GLOBAL THERMOSTAT ON THE ATOMS                               ==
       ! ==--------------------------------------------------------------==
    ELSE
       ! ==--------------------------------------------------------------==
       ! == CALCULATE ALL NOSE-HOOVER CHAIN FORCES                       ==
       ! == INCLUDING CELL THERMOSTAT FORCES                             ==
       ! ==--------------------------------------------------------------==
       ! PARTICLE KINETIC ENERGY
       ekinp=0._real_8
       DO k=1,3
          DO is=1,ions1%nsp
             const=0.5_real_8*pmx(is)
             DO ia=1,ions0%na(is)
                ekinp=ekinp+const*velp(k,ia,is)*velp(k,ia,is)
             ENDDO
          ENDDO
       ENDDO
       ! PARTICLE THERMOSTAT FORCES
       pkewant = 0.5_real_8*gkt(ipp,1)
       fetapv(1) = 2.0_real_8*(ekinp - pkewant)
       DO l=2,nchain
          ekinp = 0.5_real_8*qnospc(l-1,ip)*etap1dot(l-1,ip)*&
               etap1dot(l-1,ip)
          pkewant = 0.5_real_8*cntr%tempw/factem
          fetapv(l) = 2.0_real_8*(ekinp - pkewant)
       ENDDO
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
       etap1dot(nchain,ip) = etap1dot(nchain,ip) +&
            0.25_real_8*step*fetapv(nchain)/qnospc(nchain,ip)
       etcdot(nchc,ip) = etcdot(nchc,ip) +&
            0.25_real_8*step*fetc(nchc)/qnoscc(nchc)
       ! ==--------------------------------------------------------------==
       ! ==  LOOP OVER REST OF NOSE-HOOVER CHAIN                         ==
       ! ==--------------------------------------------------------------==
       DO l=1,nchain-1
          aa = EXP(-0.125_real_8*step*etap1dot(nchain+1-l,ip))
          etap1dot(nchain-l,ip) = etap1dot(nchain-l,ip)*aa*aa +&
               0.25_real_8*step*fetapv(nchain-l)*aa/qnospc(nchain-l,ip)
       ENDDO
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
          scmat(k,k) = scmat(k,k) + trhv/at3 + etap1dot(1,ip)
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
       DO l=1,nchain
          etap1(l,ip) = etap1(l,ip) + 0.5_real_8*step*etap1dot(l,ip)
       ENDDO
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
       ! PARTICLE KINETIC ENERGY
       ekinp=0._real_8
       DO k=1,3
          DO is=1,ions1%nsp
             const=0.5_real_8*pmx(is)
             DO ia=1,ions0%na(is)
                ekinp=ekinp+const*velp(k,ia,is)*velp(k,ia,is)
             ENDDO
          ENDDO
       ENDDO
       ! PARTICLE THERMOSTAT FORCES
       pkewant = 0.5_real_8*gkt(ipp,1)
       fetapv(1) = 2.0_real_8*(ekinp - pkewant)
       ! CELL KINETIC ENERGY
       ekinh=0.5_real_8*cntr%cmass*ddot(9,metr_com%htvel,1,metr_com%htvel,1)
       ! CELL THERMOSTAT FORCES
       hkewant = 0.5_real_8*gckt
       fetc(1) = 2._real_8*(ekinh - hkewant)
       ! ==-------------------------------------------------------------==
       ! ==  LOOP OVER CHAIN AGAIN                                      ==
       ! ==-------------------------------------------------------------==
       DO l=1,nchain-1
          aa = EXP(-0.125_real_8*step*etap1dot(l+1,ip))
          etap1dot(l,ip) = etap1dot(l,ip)*aa*aa +&
               0.25_real_8*step*fetapv(l)*aa/qnospc(l,ip)
          ekinp = 0.5_real_8*qnospc(l,ip)*etap1dot(l,ip)*etap1dot(l,ip)
          pkewant = 0.5_real_8*cntr%tempw/factem
          fetapv(l+1) = 2.0_real_8*(ekinp - pkewant)
       ENDDO
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
       etap1dot(nchain,ip) = etap1dot(nchain,ip) +&
            0.25_real_8*step*fetapv(nchain)/qnospc(nchain,ip)
       etcdot(nchc,ip) = etcdot(nchc,ip) +&
            0.25_real_8*step*fetc(nchc)/qnoscc(nchc)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE prpcnosmove
  ! ==================================================================












END MODULE prpcnosmove_utils
