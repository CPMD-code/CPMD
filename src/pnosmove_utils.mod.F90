MODULE pnosmove_utils
  USE cnst,                            ONLY: factem
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE nose,                            ONLY: &
       etap, etap1, etap1dot, etapdot, etapm, etapmdot, gkt, lctmemb, loct, &
       loctpin, mapdof, nchain, nchx, ndfnt, nlctmbm, nosl, ntherm, qnosp, &
       qnospc, qnospm, tcafes, tempwr
  USE system,                          ONLY: cntl,&
                                             cntr

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: pnosmove

CONTAINS

  ! ==================================================================
  SUBROUTINE pnosmove(velp,step,pmx,ip)
    ! ==--------------------------------------------------------------==
    ! ==  ALLOWS FOR NOSE-HOOVER CHAIN DYNAMICS                       ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: velp(:,:,:), step, pmx(*)
    INTEGER                                  :: ip

    INTEGER                                  :: ia, ial, idf, ij, ijj, il, &
                                                ipp, is, k, kf, l
    REAL(real_8)                             :: aa, const, ekinp, facton, &
                                                fetapv(nchx), pkewant

    ipp=1
    IF (cntl%tpath.AND.cntl%tpimd)ipp=MIN(ip,2)
    IF (nosl%tultra) THEN
       facton=0.5_real_8/factem
       DO is=1,ions1%nsp
          ! ==--------------------------------------------------------------==
          ! ==  CALCULATE ALL THERMOSTAT FORCES                             ==
          ! ==--------------------------------------------------------------==
          ekinp=0._real_8
          const=0.5_real_8*pmx(is)
          DO ia=1,ions0%na(is)
             DO k=1,3
                ekinp=ekinp+const*velp(k,ia,is)*velp(k,ia,is)
             ENDDO
          ENDDO
          pkewant = 0.5_real_8*gkt(ipp,is)
          fetapv(1) = 2.0_real_8*(ekinp - pkewant)
          DO l=2,nchain
             ekinp = 0.5_real_8*qnosp(is,l-1,ip)*etapdot(is,l-1,ip)*&
                  etapdot(is,l-1,ip)
             pkewant = cntr%tempw*facton
             fetapv(l) = 2.0_real_8*(ekinp - pkewant)
          ENDDO
          ! ==--------------------------------------------------------------==
          ! ==  PROPAGATE LAST CHAIN ELEMENT                                ==
          ! ==--------------------------------------------------------------==
          etapdot(is,nchain,ip) = etapdot(is,nchain,ip) +&
               0.25_real_8*step*fetapv(nchain)/qnosp(is,nchain,ip)
          ! ==--------------------------------------------------------------==
          ! ==  PROPAGATE REMAINDER OF CHAIN                                ==
          ! ==--------------------------------------------------------------==
          DO l=1,nchain-1
             aa = EXP(-0.125_real_8*step*etapdot(is,nchain+1-l,ip))
             etapdot(is,nchain-l,ip) = etapdot(is,nchain-l,ip)*aa*aa +&
                  0.25_real_8*step*fetapv(nchain-l)*aa/qnosp(is,nchain-l,ip)
          ENDDO
          ! ==--------------------------------------------------------------==
          ! ==  DO VELOCITY SCALING AND PROPAGATE THERMOSTAT POSITIONS      ==
          ! ==--------------------------------------------------------------==
          aa = EXP(-0.5_real_8*step*etapdot(is,1,ip))
          DO ia=1,ions0%na(is)
             DO k=1,3
                velp(k,ia,is) = velp(k,ia,is)*aa
             ENDDO
          ENDDO
          DO l=1,nchain
             etap(is,l,ip) = etap(is,l,ip) + 0.5_real_8*step*etapdot(is,l,ip)
          ENDDO
          ekinp=0._real_8
          const=0.5_real_8*pmx(is)
          DO ia=1,ions0%na(is)
             DO k=1,3
                ekinp=ekinp+const*velp(k,ia,is)*velp(k,ia,is)
             ENDDO
          ENDDO
          pkewant = 0.5_real_8*gkt(ipp,is)
          fetapv(1) = 2.0_real_8*(ekinp - pkewant)
          ! ==--------------------------------------------------------------==
          ! == LOOP AGAIN OVER CHAIN                                        ==
          ! ==--------------------------------------------------------------==
          DO l=1,nchain-1
             aa = EXP(-0.125_real_8*step*etapdot(is,l+1,ip))
             etapdot(is,l,ip) = etapdot(is,l,ip)*aa*aa +&
                  0.25_real_8*step*fetapv(l)*aa/qnosp(is,l,ip)
             ekinp = 0.5_real_8*qnosp(is,l,ip)*etapdot(is,l,ip)*&
                  etapdot(is,l,ip)
             pkewant = cntr%tempw*facton
             fetapv(l+1) = 2.0_real_8*(ekinp - pkewant)
          ENDDO
          ! ==--------------------------------------------------------------==
          ! ==  PROPAGATE LAST CHAIN ELEMENT                                ==
          ! ==--------------------------------------------------------------==
          etapdot(is,nchain,ip) = etapdot(is,nchain,ip) +&
               0.25_real_8*step*fetapv(nchain)/qnosp(is,nchain,ip)
       ENDDO

       ! ==--------------------------------------------------------------==
       ! ==  LOCAL THERMOSTATS                                           ==
       ! ==--------------------------------------------------------------==
    ELSEIF (loct%tloct) THEN
       facton=0.5_real_8/factem
       ! FIXME LOCT: ABOUT ALL
       DO il=1,loct%nloct
          ! ==--------------------------------------------------------------==
          ! ==  CALCULATE ALL THERMOSTAT FORCES                             ==
          ! ==--------------------------------------------------------------==
          ekinp=0._real_8
          DO ial=1,nlctmbm(il)
             is=lctmemb(1,ial,il)
             ia=lctmemb(2,ial,il)
             const=0.5_real_8*pmx(is)
             DO k=1,3
                ekinp=ekinp+const*velp(k,ia,is)*velp(k,ia,is)
             ENDDO
          ENDDO
          pkewant = 0.5_real_8*gkt(ipp,il)
          fetapv(1) = 2.0_real_8*(ekinp - pkewant)
          DO l=2,nchain
             ekinp = 0.5_real_8*qnosp(il,l-1,ip)*etapdot(il,l-1,ip)*&
                  etapdot(il,l-1,ip)
             pkewant = 0.5_real_8*loctpin(1,il)/factem
             fetapv(l) = 2.0_real_8*(ekinp - pkewant)
          ENDDO
          ! ==--------------------------------------------------------------==
          ! ==  PROPAGATE LAST CHAIN ELEMENT                                ==
          ! ==--------------------------------------------------------------==
          etapdot(il,nchain,ip) = etapdot(il,nchain,ip) +&
               0.25_real_8*step*fetapv(nchain)/qnosp(il,nchain,ip)
          ! ==--------------------------------------------------------------==
          ! ==  PROPAGATE REMAINDER OF CHAIN                                ==
          ! ==--------------------------------------------------------------==
          DO l=1,nchain-1
             aa = EXP(-0.125_real_8*step*etapdot(il,nchain+1-l,ip))
             etapdot(il,nchain-l,ip) = etapdot(il,nchain-l,ip)*aa*aa +&
                  0.25_real_8*step*fetapv(nchain-l)*aa/qnosp(il,nchain-l,ip)
          ENDDO
          ! ==--------------------------------------------------------------==
          ! ==  DO VELOCITY SCALING AND PROPAGATE THERMOSTAT POSITIONS      ==
          ! ==--------------------------------------------------------------==
          aa = EXP(-0.5_real_8*step*etapdot(il,1,ip))
          DO ial=1,nlctmbm(il)
             is=lctmemb(1,ial,il)
             ia=lctmemb(2,ial,il)
             velp(1,ia,is) = velp(1,ia,is)*aa
             velp(2,ia,is) = velp(2,ia,is)*aa
             velp(3,ia,is) = velp(3,ia,is)*aa
          ENDDO
          DO l=1,nchain
             etap(il,l,ip)=etap(il,l,ip)+0.5_real_8*step*etapdot(il,l,ip)
          ENDDO
          ekinp=0._real_8
          DO ial=1,nlctmbm(il)
             is=lctmemb(1,ial,il)
             ia=lctmemb(2,ial,il)
             const=0.5_real_8*pmx(is)
             DO k=1,3
                ekinp=ekinp+const*velp(k,ia,is)*velp(k,ia,is)
             ENDDO
          ENDDO
          pkewant = 0.5_real_8*gkt(ipp,il)
          fetapv(1) = 2.0_real_8*(ekinp - pkewant)
          ! ==--------------------------------------------------------------==
          ! == LOOP AGAIN OVER CHAIN                                        ==
          ! ==--------------------------------------------------------------==
          DO l=1,nchain-1
             aa = EXP(-0.125_real_8*step*etapdot(il,l+1,ip))
             etapdot(il,l,ip) = etapdot(il,l,ip)*aa*aa +&
                  0.25_real_8*step*fetapv(l)*aa/qnosp(il,l,ip)
             ekinp = 0.5_real_8*qnosp(il,l,ip)*etapdot(il,l,ip)*&
                  etapdot(il,l,ip)
             pkewant = 0.5_real_8*loctpin(1,il)/factem
             fetapv(l+1) = 2.0_real_8*(ekinp - pkewant)
          ENDDO
          ! ==--------------------------------------------------------------==
          ! ==  PROPAGATE LAST CHAIN ELEMENT                                ==
          ! ==--------------------------------------------------------------==
          etapdot(il,nchain,ip) = etapdot(il,nchain,ip) +&
               0.25_real_8*step*fetapv(nchain)/qnosp(il,nchain,ip)
       ENDDO
    ELSEIF (nosl%tmnose) THEN
       ! ==--------------------------------------------------------------==
       ! ==  LOOP BACKWARDS OVER NOSE-HOOVER CHAIN                       ==
       ! ==--------------------------------------------------------------==
       facton=0.5_real_8/factem
       ij = 0
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
             IF (tcafes) THEN
                pkewant = tempwr(k,ia,is)*facton
             ELSE
                pkewant = cntr%tempw*facton
             ENDIF
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
             etapmdot(kf,nchain-l,ip) &
                  =etapmdot(kf,nchain-l,ip)&
                  *aa*aa+0.25_real_8*step*fetapv(nchain-l)&
                  *aa/qnospm(kf,nchain-l,ip)
          ENDDO
          ! ==--------------------------------------------------------------==
          ! == DO VELOCITY SCALING AND PROPAGATE THERMOSTAT POSITIONS       ==
          ! ==--------------------------------------------------------------==
          aa = EXP(-0.5_real_8*step*etapmdot(kf,1,ip))
          const=0.5_real_8*pmx(is)
          ekinp=0._real_8
          DO idf=1,ndfnt(kf,ipp)
             ijj=ij+idf
             k=mapdof(1,ijj,ipp)
             ia=mapdof(2,ijj,ipp)
             is=mapdof(3,ijj,ipp)
             velp(k,ia,is) = velp(k,ia,is)*aa
             ekinp=ekinp+const*velp(k,ia,is)*velp(k,ia,is)
          ENDDO
          pkewant = 0.5_real_8*gkt(ipp,kf)
          fetapv(1) = 2.0_real_8*(ekinp - pkewant)
          DO l=1,nchain
             etapm(kf,l,ip) = etapm(kf,l,ip) +&
                  0.5_real_8*step*etapmdot(kf,l,ip)
          ENDDO
          ! ==-------------------------------------------------------------==
          ! ==  LOOP OVER CHAIN AGAIN                                      ==
          ! ==-------------------------------------------------------------==
          DO l=1,nchain-1
             aa = EXP(-0.125_real_8*step*etapmdot(kf,l+1,ip))
             etapmdot(kf,l,ip) = etapmdot(kf,l,ip)*aa*aa +&
                  0.25_real_8*step*fetapv(l)*aa/qnospm(kf,l,ip)
             ekinp = 0.5_real_8*qnospm(kf,l,ip)*&
                  etapmdot(kf,l,ip)*etapmdot(kf,l,ip)
             IF (tcafes) THEN
                pkewant = tempwr(k,ia,is)*facton
             ELSE
                pkewant = cntr%tempw*facton
             ENDIF
             fetapv(l+1) = 2.0_real_8*(ekinp - pkewant)
          ENDDO
          ! ==------------------------------------------------------------==
          ! == PROPAGATE LAST CHAIN ELEMENT                               ==
          ! ==------------------------------------------------------------==
          etapmdot(kf,nchain,ip)=etapmdot(kf,nchain,ip)+&
               0.25_real_8*step*fetapv(nchain)/qnospm(kf,nchain,ip)
          ! ==------------------------------------------------------------==
          ij=ij+ndfnt(kf,ipp)
       ENDDO
    ELSE      ! regular nose-hoover
       ! ==------------------------------------------------------------==
       ! == CALCULATE ALL NOSE-HOOVER CHAIN FORCES                     ==
       ! ==------------------------------------------------------------==
       facton=0.5_real_8/factem
       ekinp=0._real_8
       DO is=1,ions1%nsp
          const=0.5_real_8*pmx(is)
          DO ia=1,ions0%na(is)
             DO k=1,3
                ekinp=ekinp+const*velp(k,ia,is)*velp(k,ia,is)
             ENDDO
          ENDDO
       ENDDO
       pkewant = 0.5_real_8*gkt(ipp,1)
       fetapv(1) = 2.0_real_8*(ekinp - pkewant)
       DO l=2,nchain
          ekinp = 0.5_real_8*qnospc(l-1,ip)*etap1dot(l-1,ip)*&
               etap1dot(l-1,ip)
          pkewant = cntr%tempw*facton
          fetapv(l) = 2.0_real_8*(ekinp - pkewant)
       ENDDO
       ! ==------------------------------------------------------------==
       ! == MOVE LAST CHAIN ELEMENT                                    ==
       ! ==------------------------------------------------------------==
       etap1dot(nchain,ip) = etap1dot(nchain,ip) +&
            0.25_real_8*step*fetapv(nchain)/qnospc(nchain,ip)
       ! ==------------------------------------------------------------==
       ! ==  LOOP OVER REST OF NOSE-HOOVER CHAIN                       ==
       ! ==------------------------------------------------------------==
       DO l=1,nchain-1
          aa = EXP(-0.125_real_8*step*etap1dot(nchain+1-l,ip))
          etap1dot(nchain-l,ip) = etap1dot(nchain-l,ip)*aa*aa +&
               0.25_real_8*step*fetapv(nchain-l)*aa/qnospc(nchain-l,ip)
       ENDDO
       ! ==--------------------------------------------------------------==
       ! ==  DO VELOCITY SCALING AND PROPAGATE THERMOSTAT POSITIONS      ==
       ! ==--------------------------------------------------------------==
       aa = EXP(-0.5_real_8*step*etap1dot(1,ip))
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             DO k=1,3
                velp(k,ia,is) = velp(k,ia,is)*aa
             ENDDO
          ENDDO
       ENDDO
       ! 
       DO l=1,nchain
          etap1(l,ip) = etap1(l,ip) + 0.5_real_8*step*etap1dot(l,ip)
       ENDDO
       ! ==--------------------------------------------------------------==
       ! == RECALCULATE CHAIN FORCE 1                                    ==
       ! ==--------------------------------------------------------------==
       ekinp=0._real_8
       DO is=1,ions1%nsp
          const=0.5_real_8*pmx(is)
          DO ia=1,ions0%na(is)
             DO k=1,3
                ekinp=ekinp+const*velp(k,ia,is)*velp(k,ia,is)
             ENDDO
          ENDDO
       ENDDO
       pkewant = 0.5_real_8*gkt(ipp,1)
       fetapv(1) = 2.0_real_8*(ekinp - pkewant)
       ! ==-------------------------------------------------------------==
       ! ==  LOOP OVER CHAIN AGAIN                                      ==
       ! ==-------------------------------------------------------------==
       DO l=1,nchain-1
          aa = EXP(-0.125_real_8*step*etap1dot(l+1,ip))
          etap1dot(l,ip) = etap1dot(l,ip)*aa*aa +&
               0.25_real_8*step*fetapv(l)*aa/qnospc(l,ip)
          ekinp = 0.5_real_8*qnospc(l,ip)*etap1dot(l,ip)*etap1dot(l,ip)
          pkewant = cntr%tempw*facton
          fetapv(l+1) = 2.0_real_8*(ekinp - pkewant)
       ENDDO
       ! ==------------------------------------------------------------==
       ! == MOVE LAST CHAIN ELEMENT                                    ==
       ! ==------------------------------------------------------------==
       etap1dot(nchain,ip) = etap1dot(nchain,ip) +&
            0.25_real_8*step*fetapv(nchain)/qnospc(nchain,ip)
    ENDIF
    ! ==------------------------------------------------------------==
    RETURN
  END SUBROUTINE pnosmove
  ! ================================================================




END MODULE pnosmove_utils
