MODULE wr_temps_utils
  USE cnst,                            ONLY: factem
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_app,&
                                             fo_verb
  USE global_utils,                    ONLY: global
  USE kinds,                           ONLY: real_8
  USE nose,                            ONLY: etadot,&
                                             etapmdot,&
                                             nedof,&
                                             ntherm,&
                                             qnosee,&
                                             qnospm,&
                                             qnospc,&
                                             etap1dot,&
                                             nosl
  USE parac,                           ONLY: paral
  USE pimd,                            ONLY: grandparent,&
                                             maxnp,&
                                             np_high,&
                                             np_low,&
                                             pimd3,&
                                             pimd1
  USE system,                          ONLY: cntl
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: wr_temps

CONTAINS

  ! ==================================================================
  SUBROUTINE wr_temps(nfi,ekinc,tempp)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nfi
    REAL(real_8)                             :: ekinc(*), tempp(*)

    CHARACTER(len=10), PARAMETER             :: fname = 'PITEMP    '

    INTEGER                                  :: ip, ipp, k
    INTEGER, SAVE                            :: ifirst = fo_verb
    LOGICAL                                  :: ferror
    REAL(real_8)                             :: et1, ete1, etem, &
                                                etetemp(maxnp), etm, etp1, &
                                                etpm, etptemp(maxnp), tp1, tpm

    IF (grandparent) THEN
       IF (paral%io_parent)&
            CALL fileopen(4,fname,fo_app+ifirst,ferror)
       ifirst=0
    ENDIF
    CALL zeroing(etptemp)!,maxnp)
    CALL zeroing(etetemp)!,maxnp)
    IF (paral%parent) THEN
       DO ip=np_low,np_high
          ipp=MIN(2,ip)
          IF (cntl%tnosep) THEN
             IF (nosl%tmnose.AND..NOT.(pimd1%tcentro.AND.ip==1)) THEN
                DO k=1,ntherm(ipp)
                   etptemp(ip)=etptemp(ip)+factem*qnospm(k,1,ip)*&
                        etapmdot(k,1,ip)**2
                ENDDO
                etptemp(ip)=etptemp(ip)/REAL(ntherm(ipp),kind=real_8)
             ELSE
                etptemp(ip)=etptemp(ip)+factem*qnospc(1,ip)*etap1dot(1,ip)**2
             ENDIF
          ENDIF
          IF (cntl%tnosee)&
             etetemp(ip)=etetemp(ip)+0.5_real_8*qnosee(1)*etadot(1,ip)**2
       ENDDO
    ENDIF
    CALL global(etptemp,1)
    CALL global(etetemp,1)
    IF (grandparent) THEN
       tp1=tempp(1)
       et1=ekinc(1)
       etp1=etptemp(1)
       ete1=etetemp(1)
       tpm=0._real_8
       etm=0._real_8
       etpm=0._real_8
       etem=0._real_8
       DO ip=2,pimd3%np_total
          tpm=tpm+tempp(ip)
          etm=etm+ekinc(ip)
          etpm=etpm+etptemp(ip)
          etem=etem+etetemp(ip)
       ENDDO
       tpm=tpm/REAL(pimd3%np_total-1,kind=real_8)
       etm=etm/REAL(pimd3%np_total-1,kind=real_8)
       etpm=etpm/REAL(pimd3%np_total-1,kind=real_8)
       etem=etem/REAL(pimd3%np_total-1,kind=real_8)
       ! 
       ete1=ete1*REAL(nedof,kind=real_8)
       etem=etem*REAL(nedof,kind=real_8)
       ! WRITE TEMPERATURE MONITOR FILE 
       ! NOTE  only the first thermostat in the Nose-Hoover chain is
       ! considered in all cases
       ! NFI   time step
       ! TP1   ionic      temperature of IP=1 (centroid)  
       ! TPM   ionic      temperature averaged over IP=2,NP 
       ! ET1   electronic temperature of IP=1 (centroid)
       ! ETM   electronic temperature averaged over IP=2,NP
       ! ETP1  temperature of ionic      thermostat of IP=1 (centroid)  
       ! ETPM  temperature of ionic      thermostat averaged over IP=2,NP
       ! ETE1  temperature of electronic thermostat of IP=1 (centroid)
       ! ETEM  temperature of electronic thermostat averaged over IP=2,NP
       IF (paral%io_parent)&
            WRITE(4,'(I9,F7.1,F7.1,F10.5,F10.5,F7.1,F7.1,F10.5,F10.5)')&
            nfi,tp1,tpm,et1,etm,etp1,etpm,ete1,etem
       IF (paral%io_parent)&
            CALL fileclose(4)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE wr_temps
  ! new
  ! ==================================================================


END MODULE wr_temps_utils
