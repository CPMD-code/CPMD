MODULE totstr_utils
  USE cnst,                            ONLY: au_kb
  USE cnst_dyn,                        ONLY: dstrcv,&
                                             dstrmeta,&
                                             fvbound,&
                                             lmeta,&
                                             mdcellr
  USE kinds,                           ONLY: real_8
  USE metr,                            ONLY: metr_com
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum
  USE parac,                           ONLY: parai,&
                                             paral
  USE prcp,                            ONLY: prcp_com,&
                                             prcpl
  USE strs,                            ONLY: &
       alpha, beta, decc, degc, deht, dekin, denl, deps, desr, dexc, dlam, &
       pail, paiu
  USE symtrz_utils,                    ONLY: symstress
  USE system,                          ONLY: cntl,&
                                             parm
  USE vdwcmod,                         ONLY: vdwr

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: totstr
  PUBLIC :: dstre

CONTAINS

  ! ==================================================================
  SUBROUTINE totstr
    ! ==--------------------------------------------------------------==
    ! Variables
    INTEGER                                  :: i, j, kk
    LOGICAL                                  :: debug
    REAL(real_8)                             :: temp(3,3)

! ==--------------------------------------------------------------==
! ==  TOTAL STRESS                                                ==
! ==--------------------------------------------------------------==

    debug=.FALSE.
    ! ==--------------------------------------------------------------==

    DO kk=1,6
       paiu(alpha(kk),beta(kk)) = -( dekin(kk) + deht(kk) + dexc(kk)&
            + desr(kk) + deps(kk) + denl(kk) + vdwr%devdw(kk) )&
            + degc(kk) + dlam(kk) + decc(kk) + dstrmeta(kk)&
            + fvbound(kk) + dstrcv(kk)
       paiu(beta(kk),alpha(kk)) = paiu(alpha(kk),beta(kk))
    ENDDO

    CALL mp_sum(paiu,9,parai%allgrp)

    ! ==--------------------------------------------------------------==
    ! Symmetrisation for DEKIN and DENL (other parts symmetric)
    CALL symstress(paiu)
    ! ==--------------------------------------------------------------==
    CALL dcopy(9,paiu(1,1),1,metr_com%htfp(1,1),1)
    IF (prcpl%tzflex) THEN
       metr_com%htfp(3,3)=metr_com%htfp(3,3)-prcp_com%druck*parm%omega
    ELSE IF (prcpl%tisot.OR.cntl%tshock) THEN
       DO i=1,3
          metr_com%htfp(i,i)=metr_com%htfp(i,i)-prcp_com%druck*parm%omega
       ENDDO
    ELSE
       CALL daxpy(9,-parm%omega,prcp_com%stens(1,1),1,metr_com%htfp(1,1),1)
    ENDIF

    ! ==--------------------------------------------------------------==
    ! Transform to crystal axis basis (PAIU)
    CALL dgemm('T','N',3,3,3,1._real_8,metr_com%htm1,3,paiu,3,0._real_8,temp,3)
    CALL dgemm('N','T',3,3,3,1._real_8,temp,3,metr_com%ht,3,0._real_8,pail,3)
    ! HTFP
    CALL dgemm('N','N',3,3,3,1.0_real_8,metr_com%htfp(1,1),3,metr_com%htm1(1,1),3,0.0_real_8,&
         metr_com%htfor(1,1),3)

    IF (lmeta%lmeta_cell) THEN

       CALL mp_bcast(mdcellr%fh_cell,SIZE(mdcellr%fh_cell),parai%source,parai%allgrp)
       DO i=1,3
          DO j = 1,3
             metr_com%htfor(j,i) = metr_com%htfor(j,i) + mdcellr%fh_cell(j,i)
          ENDDO
       ENDDO
    ENDIF

    ! ==--------------------------------------------------------------==
    IF (debug) THEN
       CALL mp_sum(dekin,6,parai%allgrp)
       CALL mp_sum(deht,6,parai%allgrp)
       CALL mp_sum(dexc,6,parai%allgrp)
       CALL mp_sum(desr,6,parai%allgrp)
       CALL mp_sum(deps,6,parai%allgrp)
       CALL mp_sum(denl,6,parai%allgrp)
       CALL mp_sum(degc,6,parai%allgrp)
       CALL mp_sum(dlam,6,parai%allgrp)
       CALL mp_sum(decc,6,parai%allgrp)
       IF (paral%parent) THEN
          IF (paral%io_parent)&
               WRITE(6,9999)
9999      FORMAT(1x,4x,"DEKIN",5x,"DEHT",5x,"DEXC",5x,"DESR",&
               5x,"DEPS",5x,"DENL",5x,"DEGC",5x,"DLAM",5x,"DECC")
          DO kk=1,6
             IF (paral%io_parent)&
                  WRITE(6,'(1X,9F9.5)') -dekin(kk)/parm%omega,-deht(kk)/parm%omega,&
                  -dexc(kk)/parm%omega,&
                  -desr(kk)/parm%omega,-deps(kk)/parm%omega,-denl(kk)/parm%omega,&
                  degc(kk)/parm%omega, dlam(kk)/parm%omega,decc(kk)/parm%omega
          ENDDO
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE totstr
  ! ==================================================================
  SUBROUTINE dstre(htfor,dstress,dstressmax)
    ! ==--------------------------------------------------------------==
    ! == GIVE SUM ABS(HTFOR) and MAX COMPONENTS                       ==
    ! ==--------------------------------------------------------------==
    ! == INPUT: HTFOR(3,3) stress tensor - required stress tensor     ==
    ! == OUTPUT: DSTRESS = Sum Abs(HTFOR)                             ==
    ! ==         DSTRESSMAX = Max(HTFOR)                              ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: htfor(3,3), dstress, &
                                                dstressmax

    INTEGER                                  :: i, j
    REAL(real_8)                             :: arg

    dstress=0._real_8
    dstressmax=0._real_8
    IF (prcpl%tzflex) THEN
       ! Z-DIRECTION FLEXIBLE CELL
       dstress=ABS(htfor(3,3))
       dstressmax=dstress
    ELSE IF (prcpl%tisot) THEN
       ! ISOTROPIC CELL
       dstress=ABS(htfor(1,1)+htfor(2,2)+htfor(3,3))
    ELSE
       DO i=1,3
          DO j=1,3
             arg=ABS(htfor(j,i))
             dstress=dstress+arg
             IF (dstressmax.LT.arg) THEN
                dstressmax=arg
             ENDIF
          ENDDO
       ENDDO
    ENDIF
    dstress=dstress*au_kb/parm%omega
    RETURN
  END SUBROUTINE dstre
  ! ==================================================================

END MODULE totstr_utils
