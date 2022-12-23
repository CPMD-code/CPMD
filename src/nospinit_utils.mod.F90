MODULE nospinit_utils
  USE cnst,                            ONLY: factem,&
                                             pi
  USE ions,                            ONLY: ions1
  USE kinds,                           ONLY: real_8
  USE nose,                            ONLY: &
       etap, etap1, etap1dot, etapdot, etapm, etapmdot, loct, loctt0, mapdof, &
       nchain, nosl, ntherm, qnosp, qnospc, qnospm, tcafes, tempwr, tnosepc
  USE parac,                           ONLY: paral
  USE pimd,                            ONLY: pimd1
  USE prng_utils,                      ONLY: repprngu,&
                                             repprngu_vec
  USE system,                          ONLY: cntl,&
                                             cntr
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: nospinit

CONTAINS

  ! ==================================================================
  SUBROUTINE nospinit(ip)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ip

    INTEGER                                  :: atia, atis, atk, ipp, is, k, l
    REAL(real_8)                             :: alfa, alfa1, rnr(3), sigma

! Variables
! ==--------------------------------------------------------------==
! ==  ION NOSE PARAMETERS                                         ==
! ==  MAXWELL-BOLTZMANN SAMPLING FOR NOSE VELOCITIES              ==
! ==--------------------------------------------------------------==

    IF (cntl%tpath.AND.cntl%tpimd) THEN
       IF ((pimd1%tcentro.OR.pimd1%tringp).AND.ip.EQ.1.AND..NOT.tnosepc) RETURN
    ENDIF
    ipp=1
    IF (cntl%tpath.AND.cntl%tpimd) ipp=MIN(ip,2)
    IF (paral%parent) THEN
       WRITE(6,*) 'NOSPINIT| INITIALIZATION OF NOSE VELOCITIES'
    ENDIF
    CALL repprngu_vec(3,rnr)
    IF (nosl%tultra) THEN
       ! ..EACH ATOM TYPE HAS A SEPERATE CHAIN
       DO is=1,ions1%nsp
          DO l=1,nchain
             etap(is,l,ip) = 0.0_real_8
             sigma=SQRT(cntr%nospt0/qnosp(is,l,ip)/factem)
             CALL repprngu_vec(3,rnr)
             alfa1=2.0_real_8*pi*rnr(1)
             etapdot(is,l,ip)=SQRT(LOG(repprngu())*(-2.0_real_8))*&
                  COS(alfa1)*sigma
          ENDDO
       ENDDO
    ELSEIF (loct%tloct) THEN
       DO is=1,loct%nloct
          DO l=1,nchain
             etap(is,l,ip) = 0.0_real_8
             ! FIXME TLOCT: CHECK QNOSP
             sigma=SQRT(loctt0(is)/qnosp(is,l,ip)/factem)
             CALL repprngu_vec(3,rnr)
             alfa1=2.0_real_8*pi*rnr(1)
             etapdot(is,l,ip)=SQRT(LOG(repprngu())*(-2.0_real_8))*&
                  COS(alfa1)*sigma
          ENDDO
       ENDDO
    ELSEIF (nosl%tmnose) THEN
       ! ..EACH DEGREE OF FREEDOM HAS A SEPERATE CHAIN
       CALL zeroing(etapm(:,:,ip))!,3*maxsys%nax*maxsys%nsx*nchx)
       CALL zeroing(etapmdot(:,:,ip))!,3*maxsys%nax*maxsys%nsx*nchx)
       DO k=1,ntherm(ipp)
          DO l=1,nchain
             CALL repprngu_vec(3,rnr)
             IF (tcafes) THEN
                atk=mapdof(1,k,ipp)
                atia=mapdof(2,k,ipp)
                atis=mapdof(3,k,ipp)
                sigma=SQRT(tempwr(atk,atia,atis)/qnospm(k,l,ip)/factem)
             ELSE
                sigma=SQRT(cntr%nospt0/qnospm(k,l,ip)/factem)
             ENDIF
             alfa=2.0_real_8*pi*rnr(1)
             etapmdot(k,l,ip)=SQRT(LOG(repprngu())*(-2._real_8))*COS(alfa)*sigma
          ENDDO
       ENDDO
    ELSE
       DO l=1,nchain
          etap1(l,ip) = 0.0_real_8
          sigma=SQRT(cntr%nospt0/qnospc(l,ip)/factem)
          CALL repprngu_vec(3,rnr)
          alfa1=2.0_real_8*pi*rnr(1)
          etap1dot(l,ip)=SQRT(LOG(repprngu())*(-2.0_real_8))*&
               COS(alfa1)*sigma
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE nospinit
  ! ==================================================================



END MODULE nospinit_utils
