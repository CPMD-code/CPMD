MODULE adapttol_utils
  USE kinds,                           ONLY: real_8
  USE lscal,                           ONLY: adtolr,&
                                             nrllst,&
                                             tnofor
  USE mp_interface,                    ONLY: mp_bcast
  USE parac,                           ONLY: parai,&
                                             paral
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             cntr
#include "sizeof.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: tol_init
  PUBLIC :: tol_chk_force
  PUBLIC :: tol_chk_cnvgrad
  PUBLIC :: tol_chk_cnvener
  PUBLIC :: tol_det_grad
  PUBLIC :: tol_det_ener
  PUBLIC :: tol_init_shared

CONTAINS

  ! ==================================================================
  SUBROUTINE tol_init (ltolad, deacc, genvmx)
    ! ==--------------------------------------------------------------==
    ! == INITIALIZATION OF ELECTRONIC CONVERGENCE CRITERIA            ==
    ! ==--------------------------------------------------------------==

    ! ==--------------------------------------------------------------==
    ! == ARGUMENTS:                                                   ==
    ! == LTOLAD (out) Adaptive tolerance on by default                ==
    ! == DEACC  (out) No accepted step yet (checkpoint: DELAST kept)  ==
    ! == GENVMX (out) No maximum gradient of environment (TOLOG kept) ==
    ! ==--------------------------------------------------------------==
    LOGICAL                                  :: ltolad
    REAL(real_8)                             :: deacc, genvmx

! ==--------------------------------------------------------------==
! THESE VARIABLES NEED TO BE BROADCAST

    IF (adtolr%tolmin.EQ.0.0_real_8) THEN ! INIT IF NO TOLERANCES FROM CHECKPOINT
       adtolr%delast = 0.0_real_8
       adtolr%tolmin = cntr%tolog
       IF (cntr%tolini.NE.0.0_real_8) THEN
          cntr%tolog = cntr%tolini
       ENDIF
    ENDIF
    IF (paral%parent.AND.cntr%tolene.NE.0.0_real_8) THEN
       IF (adtolr%delast.EQ.0.0_real_8) THEN
          adtolr%delast = 1.0e3_real_8 * adtolr%demin
          IF (paral%io_parent)&
               WRITE (6,'(1X,A,1PE10.4)')&
               'DEFAULT WAVEFUNCTION ENERGY TOLERANCE:   ', adtolr%delast
       ELSE
          IF (paral%io_parent)&
               WRITE (6,'(1X,A,1PE10.4)')&
               'WAVEFUNCTION ENERGY TOLERANCE:   ', adtolr%delast
       ENDIF
    ENDIF
    CALL mp_bcast(cntr%tolene,parai%source, parai%allgrp)
    CALL mp_bcast(adtolr%tolmin,parai%source, parai%allgrp)
    CALL mp_bcast(adtolr%delast,parai%source, parai%allgrp)
    ! THESE VARIABLES DO NOT NEED TO BE BROADCAST (FOR PARENT USE ONLY)
    ltolad = .TRUE.
    genvmx = 0.0_real_8
    deacc  = 0.0_real_8
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE tol_init
  ! ==================================================================
  SUBROUTINE tol_chk_force (tfor, calste, gemax, infw)
    ! ==--------------------------------------------------------------==
    ! == CHECK IF FORCES ON IONS NEED TO BE CALCULATED                ==
    ! ==--------------------------------------------------------------==

    ! ==--------------------------------------------------------------==
    ! == ARGUMENTS:                                                   ==
    ! == LFOR   (out) Switch for forces                               ==
    ! == CALSTE (out) Switch for stress tensor                        ==
    ! == GEMAX  (in)  Maximum component of electronic gradient        ==
    ! == INFW   (in)  Number of WFN step for current geometry         ==
    ! ==--------------------------------------------------------------==
    LOGICAL                                  :: tfor, calste
    REAL(real_8)                             :: gemax
    INTEGER                                  :: infw

! ==--------------------------------------------------------------==

    IF (infw.EQ.1) THEN
       calste = .FALSE.
       tfor   = .FALSE.
       tnofor = .FALSE.    ! Based on energy change
    ELSE IF (gemax.LT.adtolr%tolrel*cntr%tolfor*cntr%tolog .AND. .NOT.tnofor) THEN
       calste = cntl%tprcp
       tfor   = .TRUE.
    ENDIF
    IF (cnti%nstcnv.GE.1) THEN
       IF (infw.EQ.1) THEN
          IF (adtolr%tolrel.NE.1.0_real_8.AND.paral%io_parent)&
               WRITE (6,'(1X,A,1F4.1)')&
               'RESETTING WAVEFUNCTION TOLERANCE FACTOR TO: ', 1.0_real_8
          adtolr%tolrel = 1.0_real_8
          nrllst = 1
       ELSE IF (infw-nrllst.GE.cnti%nstcnv) THEN
          adtolr%tolrel = adtolr%tolrel * 2.0_real_8
          nrllst = infw
          IF (paral%io_parent)&
               WRITE (6,'(1X,A,1F4.1)')&
               'RELAXING WAVEFUNCTION TOLERANCE, NEW FACTOR: ',&
               adtolr%tolrel
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE tol_chk_force
  ! ==================================================================
  SUBROUTINE tol_chk_cnvgrad (convwf, exsoft, gemax)
    ! ==--------------------------------------------------------------==
    ! == CHECK IF WAVEFUNCTION IS CONVERGED BASED ON GRADIENT         ==
    ! ==--------------------------------------------------------------==

    ! ==--------------------------------------------------------------==
    ! == ARGUMENTS:                                                   ==
    ! == CONVWF (upd) Flag if converged                               ==
    ! == EXSOFT (in)  Flag if exit requested                          ==
    ! == GEMAX  (in)  Maximum wavefunction component (BROADCAST)      ==
    ! ==--------------------------------------------------------------==
    LOGICAL                                  :: convwf, exsoft
    REAL(real_8)                             :: gemax

! ==--------------------------------------------------------------==

    IF (gemax.LT.adtolr%tolrel*cntr%tolog .OR. exsoft) THEN
       convwf=.TRUE.
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE tol_chk_cnvgrad
  ! ==================================================================
  SUBROUTINE tol_chk_cnvener (convwf, etot, etoto)
    ! ==--------------------------------------------------------------==
    ! == CHECK IF WAVEFUNCTION IS NOT CONVERGED BASED ON ENERGY       ==
    ! ==--------------------------------------------------------------==

    ! ==--------------------------------------------------------------==
    ! == IMPLICIT RESULTS:                                            ==
    ! == TNOFOR (out) Force forces not to be calculated               ==
    ! ==                                                              ==
    ! == ARGUMENTS:                                                   ==
    ! == CONVWF (upd) Flag if converged                               ==
    ! == ETOT   (in)  Current energy (MUST BE BROADCAST ALREADY)      ==
    ! == ETOTO  (in)  Energy of previous wfn step (BROADCAST)         ==
    ! ==--------------------------------------------------------------==
    LOGICAL                                  :: convwf
    REAL(real_8)                             :: etot, etoto

! ==--------------------------------------------------------------==

    IF (cntr%tolene.NE.0.0_real_8 .AND. adtolr%delast.NE.0.0_real_8) THEN
       IF (ABS(etot-etoto).GT.adtolr%tolrel*adtolr%delast) THEN
          convwf = .FALSE.
       ENDIF
       tnofor = (ABS(etot-etoto).GT.adtolr%tolrel*cntr%tolfor*adtolr%delast)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE tol_chk_cnvener
  ! ==================================================================
  SUBROUTINE tol_det_grad (ltolad, gnmax, genvmx)
    ! ==--------------------------------------------------------------==
    ! == DETERMINE THE WAVEFUNCTION GRADIENT TOLERANCE FOR NEXT STEP  ==
    ! ==--------------------------------------------------------------==

    ! ==--------------------------------------------------------------==
    ! == IMPLICIT RESULTS:                                            ==
    ! == TOLOG  (out) Gradient tolerance                              ==
    ! ==                                                              ==
    ! == ARGUMENTS:                                                   ==
    ! == LTOLAD (in)  Switch for adaptive tolerance                   ==
    ! == GNMAX  (in)  Maximum gradient component on all ions          ==
    ! == GENVMX (in)  Maximum gradient component on environment ions  ==
    ! ==--------------------------------------------------------------==
    LOGICAL                                  :: ltolad
    REAL(real_8)                             :: gnmax, genvmx

! ==--------------------------------------------------------------==

    IF (paral%parent) THEN
       IF (ltolad .AND. cntr%tolad.GT.0.0_real_8) THEN
          IF (genvmx.EQ.0.0_real_8) THEN
             adtolr%gnmin = MIN(gnmax,adtolr%gnmin)
          ELSE
             adtolr%gnmin = MIN(genvmx,adtolr%gnmin)
          ENDIF
          cntr%tolog = MAX(cntr%tolad*adtolr%gnmin,adtolr%tolmin)
          IF (paral%io_parent)&
               WRITE(6,'(1X,A,1PE10.4)')&
               'WAVEFUNCTION GRADIENT TOLERANCE: ', cntr%tolog
          IF (adtolr%tolmin.GE.cntr%tolog) THEN
             cntr%tolad = 0.0_real_8
             IF (paral%io_parent)&
                  WRITE (6,'(1X,A,A)') 'MINIMUM TOLERANCE REACHED, ',&
                  'SWITCHING OFF ADAPTIVE TOLERANCE'
          ENDIF
       ELSE
          cntr%tolog = adtolr%tolmin
       ENDIF
    ENDIF
    CALL mp_bcast(cntr%tolog,parai%source, parai%allgrp)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE tol_det_grad
  ! ==================================================================
  SUBROUTINE tol_det_ener (ltolad, infi, detot, deacc)
    ! ==--------------------------------------------------------------==
    ! == DETERMINE THE NEXT WAVEFUNCTION ENERGY CHANGE TOLERANCE      ==
    ! ==--------------------------------------------------------------==

    ! ==--------------------------------------------------------------==
    ! == IMPLICIT RESULTS:                                            ==
    ! == TOLOG  (out) Gradient tolerance                              ==
    ! ==                                                              ==
    ! == ARGUMENTS:                                                   ==
    ! == LTOLAD (in)  Switch for adaptive tolerance                   ==
    ! == INFI   (in)  Number of geometry step                         ==
    ! == DETOT  (in)  Energy change of last ions step                 ==
    ! == DEACC  (in)  Energy change of last successful ions step      ==
    ! ==--------------------------------------------------------------==
    LOGICAL                                  :: ltolad
    INTEGER                                  :: infi
    REAL(real_8)                             :: detot, deacc

! ==--------------------------------------------------------------==

    IF (paral%parent) THEN
       IF (.FALSE..AND..NOT.ltolad) THEN
          adtolr%delast = 0.0_real_8
       ELSE IF (cntr%tolene.GT.0.0_real_8 .AND. infi.GT.1) THEN
          IF (deacc.EQ.0.0_real_8) THEN
             adtolr%delast = detot*cntr%tolene
          ELSE
             adtolr%delast = deacc*cntr%tolene
          ENDIF
          IF (adtolr%delast.NE.0.0_real_8) THEN
             adtolr%delast = MAX(ABS(adtolr%delast),adtolr%demin)
             IF (paral%io_parent)&
                  WRITE (6,'(1X,A,1PE10.4)')&
                  'WAVEFUNCTION ENERGY TOLERANCE:   ', adtolr%delast
          ENDIF
       ENDIF
    ENDIF                     ! IF (PARENT)
    CALL mp_bcast(cntr%tolene,parai%source, parai%allgrp)
    CALL mp_bcast(adtolr%delast,parai%source, parai%allgrp)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE tol_det_ener
  ! ==================================================================
  SUBROUTINE tol_init_shared
    ! ==--------------------------------------------------------------==
    ! == Block data: init variables for P-cntl%rfo and MI TS search        ==
    ! ==--------------------------------------------------------------==

    ! ==--------------------------------------------------------------==
    tnofor = .TRUE.            ! No forces on ions
    nrllst = 0                 ! Last WF step convergence relaxed
    adtolr%delast = 0.0_real_8             ! No successful step yet (BCAST)
    adtolr%demin  = 1.0e-7_real_8            ! Force to be converged (no BCAST)
    adtolr%gnmin  = 1.0_real_8             ! Lowest gradient comp. so far
    adtolr%tolmin = 0.0_real_8             ! Minimum GNMAX (switch - BCAST)
    adtolr%tolrel = 1.0_real_8             ! Full convergence criteria
    ! ==--------------------------------------------------------------==
  END SUBROUTINE tol_init_shared
  ! ==================================================================

END MODULE adapttol_utils
