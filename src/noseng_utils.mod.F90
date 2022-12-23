MODULE noseng_utils
  USE bsym,                            ONLY: bsclcs
  USE cnst,                            ONLY: factem
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_app
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE nose,                            ONLY: &
       cafescl, cafesini, cafesnse, cafestmp, chainer, eta, etadot, etap, &
       etap1, etap1dot, etapdot, etapm, etapmdot, etc, etcdot, gkt, loct, &
       loctpin, mapdof, ncafesgrp, ncdof, nchain, nchc, nche, nedof, nosl, &
       ntherm, qnoscc, qnosee, qnosp, qnospc, qnospm, tcafes, tempwr
  USE parac,                           ONLY: paral
  USE pimd,                            ONLY: pimd1,&
                                             pma0s
  USE rmas,                            ONLY: rmass
  USE system,                          ONLY: cntl,&
                                             cntr,&
                                             maxsys
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: noseng

CONTAINS

  ! ==================================================================
  SUBROUTINE noseng(nfi,velp,enose,enosp,enosc,ip)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nfi
    REAL(real_8)                             :: velp(:,:,:), enose, enosp, &
                                                enosc
    INTEGER                                  :: ip

    CHARACTER(len=20)                        :: formatstring
    INTEGER                                  :: atia, atis, atk, i, ia, ipp, &
                                                is, k, l
    LOGICAL                                  :: ferror
    REAL(real_8)                             :: ds1, ds2, ekinosep, ekp, epp, &
                                                potnosep, pmx(maxsys%nsx)

    enose=0.0_real_8
    enosp=0.0_real_8
    enosc=0.0_real_8
    ekinosep=0.0_real_8
    potnosep=0.0_real_8
    IF (cntl%tnosee) THEN
       enose = enose + 0.5_real_8*qnosee(1)*etadot(1,ip)*etadot(1,ip) +&
            2.0_real_8*cntr%ekinw*eta(1,ip)
       DO i=2,nche
          enose = enose + 0.5_real_8*qnosee(i)*etadot(i,ip)*etadot(i,ip) +&
               2.0_real_8*cntr%ekinw*eta(i,ip)/REAL(nedof,kind=real_8)
       ENDDO
       enose = enose + 2._real_8*cntr%ekinw*eta(nche-1,ip)/REAL(nedof,kind=real_8)
    ENDIF
    ! NN: IF BS & HS-WF SKIP FOR IONIC & CELL THERMOSTAT
    IF (cntl%bsymm.AND.(bsclcs.EQ.2))RETURN
    ! 
    IF (cntl%tnosep) THEN
       ! mb  Write Nose-Hoover trajectory and energies in the case of global 
       ! mb  or local thermostats to be used for Green-Kubo of heat flux as in
       ! mb  Phys. Rev. B 76, 085424 (2007). J=-Sum_i ETAP(i)*VELP(i)**2/MASS(i)
       CALL fileopen (23,'NOSE_TRAJEC',fo_app,ferror)! write Nose trajectory
       CALL fileopen (26,'NOSE_ENERGY',fo_app,ferror)! write Nose energies
       ! mb - add Nose_ENERGIES
       ipp=1
       IF (cntl%tpath.AND.cntl%tpimd)ipp=MIN(ip,2)
       IF (nosl%tultra) THEN
          DO is=1,ions1%nsp
             enosp = enosp + 0.5_real_8*qnosp(is,1,ip)*etapdot(is,1,ip)*&
                  etapdot(is,1,ip) + gkt(ipp,is)*etap(is,1,ip)
             DO l=2,nchain
                enosp = enosp + 0.5_real_8*qnosp(is,l,ip)*etapdot(is,l,ip)*&
                     etapdot(is,l,ip) + cntr%tempw*etap(is,l,ip)/factem
             ENDDO
          ENDDO
       ELSEIF (loct%tloct) THEN
          DO is=1,loct%nloct
             l=1
             ekp=0.5_real_8*qnosp(is,1,ip)*etapdot(is,1,ip)*etapdot(is,1,ip)
             epp=gkt(ipp,is)*etap(is,1,ip)
             ekinosep=ekinosep+ekp
             potnosep=potnosep+epp
             enosp = enosp + ekp + epp
             DO l=2,nchain
                ekp=0.5_real_8*qnosp(is,l,ip)*etapdot(is,l,ip)*etapdot(is,l,ip)
                epp=loctpin(1,is)*etap(is,l,ip)/factem
                ekinosep=ekinosep+ekp
                potnosep=potnosep+epp
                enosp = enosp + ekp + epp
             ENDDO
             ! mb-start
             DO l=1,nchain
                WRITE(26,'(i10,2x,3i6,1x,3e16.6)')&
                     NFI,IP,IS,L,EKINOSEP,POTNOSEP,ENOSP
                WRITE(23,'(i10,2x,3i6,1x,2e16.6)')&
                     NFI,IP,IS,L,ETAP(IS,L,IP),ETAPDOT(IS,L,IP)    
             ENDDO
             ! mb-end
          ENDDO
       ELSEIF (nosl%tmnose) THEN
          IF (tcafes) THEN
             CALL zeroing(chainer)!,3*maxsys%nax*maxsys%nsx)
             DO k=1,ntherm(ipp)
                atk=mapdof(1,k,ipp)
                atia=mapdof(2,k,ipp)
                atis=mapdof(3,k,ipp)
                chainer(atk,atia,atis)= chainer(atk,atia,atis)+&
                     0.5_real_8*qnospm(k,1,ip)*etapmdot(k,1,ip)*&
                     etapmdot(k,1,ip) + gkt(ipp,k)*etapm(k,1,ip)
                DO l=2,nchain
                   chainer(atk,atia,atis)=chainer(atk,atia,atis)+&
                        0.5_real_8*qnospm(k,l,ip)*&
                        etapmdot(k,l,ip)*etapmdot(k,l,ip)&
                        + tempwr(atk,atia,atis)*etapm(k,l,ip)/factem
                ENDDO
             ENDDO
             IF (cntl%tpath.AND.cntl%tpimd) THEN
                IF (pimd1%tpinm.OR.pimd1%tstage) THEN
                   CALL dcopy(maxsys%nsx,pma0s(1,ip),1,pmx(1),1)
                ELSE 
                   CALL dcopy(maxsys%nsx,rmass%pma(1),1,pmx(1),1)
                ENDIF
             ELSE
                CALL dcopy(maxsys%nsx,rmass%pma(1),1,pmx(1),1)
             ENDIF
             CALL zeroing(cafestmp)!,ncafesgrp)
             CALL zeroing(cafesnse)!,ncafesgrp)
             DO is=1,ions1%nsp
                DO ia=1,ions0%na(is)
                   ds1=0.0_real_8
                   ds2=0.0_real_8
                   DO k=1,3
                      ds1=ds1+chainer(k,ia,is)
                      ds2=ds2+velp(k,ia,is)**2
                   ENDDO
                   enosp = enosp + ds1
                   cafesnse(cafescl(ia,is))=cafesnse(cafescl(ia,is))+ds1
                   cafestmp(cafescl(ia,is))=cafestmp(cafescl(ia,is))+&
                   !!!  factem/3.0_real_8*rmass%pma(is)*ds2
                        factem/3.0_real_8*pmx(is)*ds2
                ENDDO
             ENDDO
             ! normalise the temperature by the number of particles
             DO i=1,ncafesgrp
                cafestmp(i)=cafestmp(i)/(cafesini(2,i)-cafesini(1,i)+1)
             ENDDO
             IF (paral%parent) THEN
                IF (paral%io_parent)&
                     WRITE(formatstring,'(A6,I6,A6)')&
                     '(I8,x,',2*ncafesgrp,'F12.6)'
                IF (paral%io_parent)&
                     WRITE(171,formatstring) nfi,(cafestmp(i),i=1,ncafesgrp),&
                     (cafesnse(i),i=1,ncafesgrp)
             ENDIF
          ELSE
             DO k=1,ntherm(ipp)
                enosp = enosp + 0.5_real_8*qnospm(k,1,ip)*etapmdot(k,1,ip)*&
                     etapmdot(k,1,ip) + gkt(ipp,k)*etapm(k,1,ip)
                DO l=2,nchain
                   enosp=enosp + 0.5_real_8*qnospm(k,l,ip)*&
                        etapmdot(k,l,ip)*etapmdot(k,l,ip)&
                        + cntr%tempw*etapm(k,l,ip)/factem
                ENDDO
             ENDDO
          ENDIF
       ELSE
          ekinosep=0.5_real_8*qnospc(1,ip)*etap1dot(1,ip)*etap1dot(1,ip)
          potnosep=gkt(ipp,1)*etap1(1,ip)
          enosp = ekinosep + potnosep
          DO l=2,nchain
             ekp=0.5_real_8*qnospc(l,ip)*etap1dot(l,ip)*etap1dot(l,ip)
             epp=cntr%tempw*etap1(l,ip)/factem
             ekinosep=ekinosep+ekp
             potnosep=potnosep+epp
             enosp = enosp + ekp + epp
          ENDDO
          ! mb-start
          DO l=1,nchain
             WRITE(26,'(i10,2x,3i6,1x,3e16.6)')&
                  NFI,IP,1,L,EKINOSEP,POTNOSEP,ENOSP
             WRITE(23,'(i10,2x,3i6,1x,2e16.6)')&
                  NFI,IP,1,L,ETAP1(L,IP),ETAP1DOT(L,IP)
          ENDDO
          ! mb-end
       ENDIF
       CALL fileclose(23)! cmb - close Nose trajectory
       CALL fileclose(26)! cmb - close Nose energies
    ENDIF
    IF (cntl%tnosec) THEN
       enosc = enosc + 0.5_real_8*qnoscc(1)*etcdot(1,ip)*etcdot(1,ip)&
            + REAL(ncdof,kind=real_8)*cntr%tempc*etc(1,ip)/factem
       DO l=2,nchc
          enosc = enosc + 0.5_real_8*qnoscc(l)*etcdot(l,ip)*etcdot(l,ip)&
               + cntr%tempc*etc(l,ip)/factem
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE noseng
  ! ==================================================================



END MODULE noseng_utils
