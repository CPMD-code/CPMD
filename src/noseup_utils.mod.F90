MODULE noseup_utils
  USE bsym,                            ONLY: bsclcs
  USE enosmove_utils,                  ONLY: enosmove
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_bcast
  USE nose,                            ONLY: dtsuz,&
                                             ncalls,&
                                             nit,&
                                             tnosepc
  USE parac,                           ONLY: parai,&
                                             paral
  USE pimd,                            ONLY: pimd1,&
                                             pma0s
  USE pnosmove_utils,                  ONLY: pnosmove
  USE prcnosmove_utils,                ONLY: prcnosmove,&
                                             prcnosmove_iso
  USE prcp,                            ONLY: prcpl
  USE prpcmove_utils,                  ONLY: prpcmove,&
                                             prpcmove_iso
  USE prpcnosmove_utils,               ONLY: prpcnosmove,&
                                             prpcnosmove_iso
  USE prpnosmove_utils,                ONLY: prpnosmove,&
                                             prpnosmove_iso
  USE rekine_utils,                    ONLY: rekine
  USE rmas,                            ONLY: rmass
  USE system,                          ONLY: cntl,&
                                             ncpw
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: noseup

CONTAINS

  ! ==================================================================
  SUBROUTINE noseup(velp,cm,nstate,ip)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: velp(:,:,:)
    COMPLEX(real_8)                          :: cm(*)
    INTEGER                                  :: nstate, ip
    CHARACTER(*), PARAMETER                  :: procedureN = 'noseup'

    INTEGER                                  :: i, isub, j
    REAL(real_8)                             :: ekinc, sctot

    CALL tiset(procedureN,isub)
    IF (cntl%tnosee) THEN
       CALL rekine(cm,nstate,ekinc)
       IF (paral%parent) THEN
          sctot=1.0_real_8
          DO i=1,nit
             DO j=1,ncalls
                CALL enosmove(ekinc,dtsuz(j),sctot,ip)
             ENDDO
          ENDDO
       ENDIF
       CALL mp_bcast(sctot,parai%source,parai%allgrp)
       CALL dscal(2*ncpw%ngw*nstate,sctot,cm,1)
    ENDIF
    ! FOR BS_CPMD: IF WF IS HS, SKIP THE REST
    IF (cntl%bsymm.AND.(bsclcs.EQ.2))THEN
       CALL tihalt(procedureN,isub)
       RETURN
    ENDIF
    !
    IF (cntl%tnosep.AND.paral%parent.AND..NOT.cntl%tnosec) THEN
       DO i=1,nit
          DO j=1,ncalls
             IF (cntl%tpath.AND.cntl%tpimd) THEN
                IF (.NOT.((pimd1%tcentro.OR.pimd1%tringp).AND.ip.EQ.1.AND..NOT.tnosepc)) THEN
                   IF (pimd1%tpinm.OR.pimd1%tstage) THEN
                      CALL pnosmove(velp,dtsuz(j),pma0s(1,ip),ip)
                   ELSE
                      CALL pnosmove(velp,dtsuz(j),rmass%pma,ip)
                   ENDIF
                ENDIF
             ELSE
                CALL pnosmove(velp,dtsuz(j),rmass%pma,ip)
             ENDIF
          ENDDO
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    ! AK: FIXME  
    IF (cntl%tprcp.AND.paral%parent.AND.prcpl%tzflex) THEN
       CALL stopgm(procedureN,'NO NOSE THERMOSTAT FOR ZFLEXIBLE CELL',& 
            __LINE__,__FILE__)
    ENDIF
    ! 
    IF (cntl%tprcp.AND.ip.EQ.1.AND.paral%parent.AND.prcpl%tisot) THEN
       IF ((cntl%tnosep.AND.tnosepc).AND..NOT.cntl%tnosec) THEN
          DO i=1,nit
             DO j=1,ncalls
                IF (cntl%tpath.AND.cntl%tpimd) THEN
                   IF (pimd1%tpinm.OR.pimd1%tstage) THEN
                      CALL prpnosmove_iso(velp,dtsuz(j),pma0s(1,ip),ip)
                   ELSE
                      CALL prpnosmove_iso(velp,dtsuz(j),rmass%pma,ip)
                   ENDIF
                ELSE
                   CALL prpnosmove_iso(velp,dtsuz(j),rmass%pma,ip)
                ENDIF
             ENDDO
          ENDDO
       ELSEIF (.NOT.(cntl%tnosep.AND.tnosepc).AND.cntl%tnosec) THEN
          DO i=1,nit
             DO j=1,ncalls
                IF (cntl%tpath.AND.cntl%tpimd) THEN
                   IF (pimd1%tpinm.OR.pimd1%tstage) THEN
                      CALL prcnosmove_iso(velp,dtsuz(j),pma0s(1,ip),ip)
                   ELSE
                      CALL prcnosmove_iso(velp,dtsuz(j),rmass%pma,ip)
                   ENDIF
                ELSE
                   CALL prcnosmove_iso(velp,dtsuz(j),rmass%pma,ip)
                ENDIF
             ENDDO
          ENDDO
       ELSEIF ((cntl%tnosep.AND.tnosepc).AND.cntl%tnosec) THEN
          DO i=1,nit
             DO j=1,ncalls
                IF (cntl%tpath.AND.cntl%tpimd) THEN
                   IF (pimd1%tpinm.OR.pimd1%tstage) THEN
                      CALL prpcnosmove_iso(velp,dtsuz(j),pma0s(1,ip),ip)
                   ELSE
                      CALL prpcnosmove_iso(velp,dtsuz(j),rmass%pma,ip)
                   ENDIF
                ELSE
                   CALL prpcnosmove_iso(velp,dtsuz(j),rmass%pma,ip)
                ENDIF
             ENDDO
          ENDDO
       ELSE
          DO i=1,nit
             DO j=1,ncalls
                IF (cntl%tpath.AND.cntl%tpimd) THEN
                   IF (pimd1%tpinm.OR.pimd1%tstage) THEN
                      CALL prpcmove_iso(velp,dtsuz(j),pma0s(1,ip))
                   ELSE
                      CALL prpcmove_iso(velp,dtsuz(j),rmass%pma)
                   ENDIF
                ELSE
                   CALL prpcmove_iso(velp,dtsuz(j),rmass%pma)
                ENDIF
             ENDDO
          ENDDO
       ENDIF
       ! ==--------------------------------------------------------------==
    ELSEIF (cntl%tprcp.AND.ip.EQ.1.AND.paral%parent.AND..NOT.prcpl%tisot.AND..NOT.cntl%tshock) THEN
       ! ==--------------------------------------------------------------==
       IF ((cntl%tnosep.AND.tnosepc).AND..NOT.cntl%tnosec) THEN
          DO i=1,nit
             DO j=1,ncalls
                IF (cntl%tpath.AND.cntl%tpimd) THEN
                   IF (pimd1%tpinm.OR.pimd1%tstage) THEN
                      CALL prpnosmove(velp,dtsuz(j),pma0s(1,ip),ip)
                   ELSE
                      CALL prpnosmove(velp,dtsuz(j),rmass%pma,ip)
                   ENDIF
                ELSE
                   CALL prpnosmove(velp,dtsuz(j),rmass%pma,ip)
                ENDIF
             ENDDO
          ENDDO
       ELSEIF (.NOT.(cntl%tnosep.AND.tnosepc).AND.cntl%tnosec) THEN
          DO i=1,nit
             DO j=1,ncalls
                IF (cntl%tpath.AND.cntl%tpimd) THEN
                   IF (pimd1%tpinm.OR.pimd1%tstage) THEN
                      CALL prcnosmove(velp,dtsuz(j),pma0s(1,ip),ip)
                   ELSE
                      CALL prcnosmove(velp,dtsuz(j),rmass%pma,ip)
                   ENDIF
                ELSE
                   CALL prcnosmove(velp,dtsuz(j),rmass%pma,ip)
                ENDIF
             ENDDO
          ENDDO
       ELSEIF ((cntl%tnosep.AND.tnosepc).AND.cntl%tnosec) THEN
          DO i=1,nit
             DO j=1,ncalls
                IF (cntl%tpath.AND.cntl%tpimd) THEN
                   IF (pimd1%tpinm.OR.pimd1%tstage) THEN
                      CALL prpcnosmove(velp,dtsuz(j),pma0s(1,ip),ip)
                   ELSE
                      CALL prpcnosmove(velp,dtsuz(j),rmass%pma,ip)
                   ENDIF
                ELSE
                   CALL prpcnosmove(velp,dtsuz(j),rmass%pma,ip)
                ENDIF
             ENDDO
          ENDDO
       ELSE
          DO i=1,nit
             DO j=1,ncalls
                IF (cntl%tpath.AND.cntl%tpimd) THEN
                   IF (pimd1%tpinm.OR.pimd1%tstage) THEN
                      CALL prpcmove(velp,dtsuz(j),pma0s(1,ip))
                   ELSE
                      CALL prpcmove(velp,dtsuz(j),rmass%pma)
                   ENDIF
                ELSE
                   CALL prpcmove(velp,dtsuz(j),rmass%pma)
                ENDIF
             ENDDO
          ENDDO
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE noseup
  ! ==================================================================



END MODULE noseup_utils
