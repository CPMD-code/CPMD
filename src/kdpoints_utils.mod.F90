MODULE kdpoints_utils
  USE error_handling,                  ONLY: stopgm
  USE kdp,                             ONLY: wkdp,&
                                             xkdp
  USE kdpc,                            ONLY: nkdp
  USE kpnt,                            ONLY: rk,&
                                             wk
  USE mp_interface,                    ONLY: mp_bcast
  USE parac,                           ONLY: parai,&
                                             paral
  USE rkpnt_utils,                     ONLY: bzmesh,&
                                             rkpnt
  USE system,                          ONLY: nkpt,&
                                             parm

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: kdpoints

CONTAINS

  SUBROUTINE kdpoints(nkpoint)
    ! 
    INTEGER                                  :: nkpoint

    CHARACTER(*), PARAMETER                  :: procedureN = 'kdpoints'

    INTEGER                                  :: ierr, ikdp, ir
    INTEGER, SAVE                            :: ifirst = 0

! 
! 
! ==--------------------------------------------------------------==
! == ALLOCATE ARRAYS FOR KDP ROUTINES                             ==
! ==--------------------------------------------------------------==

    IF (paral%parent) THEN
       CALL bzmesh(nkpoint)
       nkdp=nkpt%nkpts
       nkpt%nkpnt=1
       nkpt%nkpts=1
    ENDIF
    ! KPOINTS
    CALL mp_bcast(nkpt%nkpnt,parai%source,parai%allgrp)
    CALL mp_bcast(nkdp,parai%source,parai%allgrp)
    CALL mp_bcast(nkpt%nkpts,parai%source,parai%allgrp)
    IF (.NOT.paral%parent)  THEN
       ALLOCATE(rk(3,nkdp),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    CALL mp_bcast(rk,SIZE(rk),parai%source,parai%allgrp)
    IF (.NOT.paral%parent)  THEN
       ALLOCATE(wk(nkdp),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    CALL mp_bcast(wk,SIZE(wk),parai%source,parai%allgrp)
    IF (ifirst.EQ.0) THEN
       ifirst=1
       ALLOCATE(wkdp(nkdp),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(xkdp(3,nkdp),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! 
    ! Set up the k-point mesh using the original routine.
    CALL rkpnt
    DO ikdp=1,nkdp
       DO ir=1,3
          rk(ir,ikdp)=rk(ir,ikdp)*parm%tpiba
       ENDDO
    ENDDO
    ! 
    ! Copy the k-point arrays to the k.p arrays
    CALL dcopy(3*nkdp,rk,1,xkdp,1)
    CALL dcopy(nkdp,wk,1,wkdp,1)

    ! 
    RETURN
  END SUBROUTINE kdpoints

END MODULE kdpoints_utils
