MODULE ehrenfest_utils
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sync
  USE parac,                           ONLY: parai,&
                                             paral
  USE spin,                            ONLY: clsd
  USE system,                          ONLY: cntl,&
                                             cntr,&
                                             fpar,&
                                             nkpt
  USE td_cayley_utils,                 ONLY: td_cayley
  USE td_input,                        ONLY: itermax,&
                                             td_prop
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing
#include "sizeof.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ehrenfest

CONTAINS

  ! ==================================================================
  SUBROUTINE ehrenfest(c0,c2,rhoe,psi,sc0,eigv)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(*), c2(*)
    REAL(real_8)                             :: rhoe(:,:)
    COMPLEX(real_8)                          :: psi(:), sc0(*)
    REAL(real_8)                             :: eigv(*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'ehrenfest'

    INTEGER                                  :: ierr, ispin, isub, nstate
    INTEGER, SAVE                            :: ifirst = 0
    LOGICAL                                  :: do_dipole, fexist, &
                                                variable_time_step
    REAL(real_8)                             :: mem_tintrvll
    REAL(real_8), ALLOCATABLE                :: effpot(:)

! Variables
! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    variable_time_step=.FALSE.
    mem_tintrvll=td_prop%tintrvll
    td_prop%tt2=td_prop%tintrvll/2._real_8
    td_prop%tt=td_prop%tintrvll
    ! 
    nstate=crge%n
    ispin=clsd%nlsd
    ALLOCATE(effpot(fpar%nnr1*ispin),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(effpot)!,nnr1*ispin)
    IF (td_prop%td_extpot) THEN
       IF ((ifirst.EQ.0).AND.paral%parent) THEN
          IF (paral%io_parent) INQUIRE(file='abstime.dat',exist=fexist)
          IF (.NOT.fexist) THEN
             IF (paral%io_parent) OPEN(99,file='abstime.dat',status='unknown')
             IF (paral%io_parent) WRITE(6,*) 'WARNING: time resetted to zero'
             td_prop%ttime=0.0
          ELSE
             IF (paral%io_parent) OPEN(99,file='abstime.dat',status='OLD')
             IF (paral%io_parent) READ(99,*) td_prop%ttime
          ENDIF
          IF (paral%io_parent) CLOSE(99)
          ifirst=1
       ENDIF
       CALL mp_bcast(td_prop%ttime,parai%source,parai%allgrp)
    ELSE
       td_prop%ttime=0._real_8
    ENDIF
    ! 
    IF (.NOT.cntl%ruku) THEN
       IF (cntl%cheby) THEN
          ! if (cntl%tpspec) then
          ! do_dipole=.true.
          ! call td_cheby(c0,c2,rhoe,psi,sc0,scr,lscr,
          ! &          effpot,nstate,eigv,.false.,do_dipole)
          ! elseif (cntl%tpdist) then
          ! do_dipole=.false.
          ! call td_cheby(c0,c2,rhoe,psi,sc0,scr,lscr,
          ! &          effpot,nstate,eigv,.false.,do_dipole)
          ! else
          ! do_dipole=.false.
          ! call td_cheby(c0,c2,rhoe,psi,sc0,scr,lscr,
          ! &          effpot,nstate,eigv,.false.,do_dipole)
          ! endif
       ELSEIF (cntl%cayley) THEN
          IF (cntl%tpspec) THEN
             do_dipole=.TRUE.
             CALL td_cayley(c0,c2,rhoe,psi,sc0,&
                  effpot,nstate,eigv,.FALSE.,do_dipole)
          ELSEIF (cntl%tpdist) THEN
             do_dipole=.FALSE.
             CALL td_cayley(c0,c2,rhoe,psi,sc0,&
                  effpot,nstate,eigv,.FALSE.,do_dipole)
          ELSE
             do_dipole=.FALSE.
             CALL td_cayley(c0,c2,rhoe,psi,sc0,&
                  effpot,nstate,eigv,.FALSE.,do_dipole)
          ENDIF
       ELSE
          IF (paral%io_parent) WRITE(6,*) 'choose an integration method'
          STOP
       ENDIF
    ELSE
       IF (cntl%cheby) THEN
          ! call dcopy(2*ngwk*nstate,c0(1),1,c2,1) 
          ! if (parent) write(6,'(/,1X,A)') 'FIRST RUNGE-KUTTA STEP'
          ! n_cycles=1
          ! tintrvll=tt2
          ! if (cntl%tpspec) then
          ! do_dipole=.true.
          ! call td_cheby(c0,c2,rhoe,psi,sc0,scr,lscr,
          ! &          effpot,nstate,eigv,.true.,do_dipole)
          ! elseif (cntl%tpdist) then
          ! do_dipole=.false.
          ! call td_cheby(c0,c2,rhoe,psi,sc0,scr,lscr,
          ! &          effpot,nstate,eigv,.true.,do_dipole)
          ! else
          ! do_dipole=.false.
          ! call td_cheby(c0,c2,rhoe,psi,sc0,scr,lscr,
          ! &          effpot,nstate,eigv,.true.,do_dipole)
          ! endif
          ! call dcopy(2*ngwk*nstate,c2(1),1,c0(1),1)
          ! n_cycles=1
          ! nzeros=2*nzeros
          ! if (parent) write(6,'(/,1X,A)') 'SECOND RUNGE-KUTTA STEP'
          ! tintrvll=tt
          ! if (cntl%tpspec) then
          ! do_dipole=.true.
          ! call td_cheby(c0,c2,rhoe,psi,sc0,scr,lscr,
          ! &          effpot,nstate,eigv,.false.,do_dipole)
          ! elseif (cntl%tpdist) then
          ! do_dipole=.false.
          ! call td_cheby(c0,c2,rhoe,psi,sc0,scr,lscr,
          ! &          effpot,nstate,eigv,.false.,do_dipole)
          ! else
          ! do_dipole=.false.
          ! call td_cheby(c0,c2,rhoe,psi,sc0,scr,lscr,
          ! &          effpot,nstate,eigv,.false.,do_dipole)
          ! endif
          ! if (parent.and.td_extpot) write(6,'(A,F10.3)')
          ! &   ' New total propagation time ',ttime
          ! nzeros=nzeros/2
       ELSEIF (cntl%cayley) THEN
          td_prop%tintrvll=td_prop%tt2
          CALL dcopy(2*nkpt%ngwk*nstate,c0(1),1,c2,1)
          IF (paral%io_parent) WRITE(6,'(/,A)') ' FIRST RUNGE_KUTTA STEP'
          td_prop%n_cycles=1
          do_dipole=.TRUE.
          CALL td_cayley(c0,c2,rhoe,psi,sc0,&
               effpot,nstate,eigv,.TRUE.,do_dipole)
          CALL dcopy(2*nkpt%ngwk*nstate,c2(1),1,c0(1),1)
          td_prop%tintrvll=td_prop%tt
          IF (paral%io_parent) WRITE(6,'(/,A)') ' SECOND RUNGE_KUTTA STEP'
          td_prop%n_cycles=1
          do_dipole=.TRUE.
          CALL td_cayley(c0,c2,rhoe,psi,sc0,&
               effpot,nstate,eigv,.FALSE.,do_dipole)
          IF (paral%io_parent.AND.td_prop%td_extpot) WRITE(6,'(A,F10.3)')&
               ' New total propagation time ',td_prop%ttime
       ELSE
          IF (paral%io_parent) THEN
             WRITE(6,*) 'EH-dyn: choose an integration method'
          ENDIF
          STOP
       ENDIF
    ENDIF
    ! 
    ! ttime=ttime+tt
    td_prop%tintrvll=mem_tintrvll
    IF (paral%parent.AND.cntl%tmdeh.AND.variable_time_step) THEN
       ! if (cntl%tmdeh) then
       IF (itermax.GE.6.AND.td_prop%tintrvll.GT.5.e-5_real_8) THEN
          td_prop%tintrvll=td_prop%tintrvll/2._real_8
          cntr%delt_ions=td_prop%tintrvll
          cntr%delt_elec=td_prop%tintrvll
       ELSEIF (itermax.LE.3.AND.td_prop%tintrvll.LT.5.e-3_real_8) THEN
          td_prop%tintrvll=td_prop%tintrvll*2._real_8
          cntr%delt_ions=td_prop%tintrvll
          cntr%delt_elec=td_prop%tintrvll
       ENDIF
       ! write(88,*) 'in',me,tintrvll
    ENDIF
    CALL mp_sync(parai%allgrp)
    CALL mp_bcast(td_prop%tintrvll,parai%source,parai%allgrp)
    CALL mp_bcast(cntr%delt_ions,parai%source,parai%allgrp)
    CALL mp_bcast(cntr%delt_elec,parai%source,parai%allgrp)
    CALL mp_sync(parai%allgrp)
    ! 
    IF (paral%io_parent.AND.td_prop%td_extpot) THEN
       OPEN(99,file='abstime.dat',status='OLD')
       REWIND(99)
       WRITE(99,*) td_prop%ttime
    ENDIF
    ! 
    DEALLOCATE(effpot,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    RETURN
  END SUBROUTINE ehrenfest
  ! ==================================================================

END MODULE ehrenfest_utils
