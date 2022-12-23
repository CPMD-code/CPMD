MODULE anneal_utils
  USE bsym,                            ONLY: bsclcs
  USE cnst,                            ONLY: factem
  USE ekinpp_utils,                    ONLY: ekinpp,&
                                             s_ekinpp
  USE ions,                            ONLY: ions1
  USE kinds,                           ONLY: real_8
  USE nose,                            ONLY: glib
  USE parac,                           ONLY: paral
  USE pimd,                            ONLY: ipcurr,&
                                             pimd1
  USE system,                          ONLY: cntl,&
                                             cntr,&
                                             iatpt,&
                                             maxsys,&
                                             ncpw
  USE tpar,                            ONLY: dt_ions

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: anneal
  PUBLIC :: dampdyn
  PUBLIC :: berendsen
  PUBLIC :: tempramp

CONTAINS

  ! ==================================================================
  SUBROUTINE anneal(velp,cm,nstate,htvel)
    ! ==--------------------------------------------------------------==

    REAL(real_8)                             :: velp(:,:,:)
    COMPLEX(real_8)                          :: cm(ncpw%ngw,*)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: htvel(*)

    INTEGER                                  :: i, ia, is
    LOGICAL                                  :: skipanic
    REAL(real_8)                             :: alfap

! Variables
! 
! NN: IF BS AND HS-WF SKIP SCALING OF IONIC AND CELL VELOCITIES 

    skipanic=.FALSE.
    IF (cntl%bsymm.AND.(bsclcs.EQ.2))skipanic=.TRUE.
    ! 
    ! ==--------------------------------------------------------------==
    ! ==  ANNEALING (IONS)  Fixed rescaling factor (anner)            ==
    ! ==--------------------------------------------------------------==
    IF (cntl%annei.AND.(.NOT.skipanic)) THEN
       alfap=cntr%anneri**(0.25_real_8)
#if defined(__VECTOR)
       !$omp parallel do private(I,IS,IA)
#else
       !$omp parallel do private(I,IS,IA) schedule(static)
#endif
       DO i=1,ions1%nat
          ia=iatpt(1,i)
          is=iatpt(2,i)
          velp(1,ia,is)=alfap*velp(1,ia,is)
          velp(2,ia,is)=alfap*velp(2,ia,is)
          velp(3,ia,is)=alfap*velp(3,ia,is)
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ==  ANNEALING (ELECTRONS)  Fixed rescaling factor (anner)       ==
    ! ==--------------------------------------------------------------==
    IF (cntl%annee) THEN
       alfap=cntr%annere**(0.25_real_8)
       CALL dscal(2*ncpw%ngw*nstate,alfap,cm,1)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ==  ANNEALING (CELL)  Fixed rescaling factor (anner)            ==
    ! ==--------------------------------------------------------------==
    IF (cntl%annec.AND.(.NOT.skipanic)) THEN
       alfap=cntr%annerc**(0.25_real_8)
       DO i=1,9
          htvel(i)=alfap*htvel(i)
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE anneal
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE dampdyn(velp,fion,cm,c2,nstate,htvel,htfor)
    ! ==--------------------------------------------------------------==

    REAL(real_8)                             :: velp(3,maxsys%nax,*), &
                                                fion(:,:,:)
    COMPLEX(real_8)                          :: cm(ncpw%ngw,*), c2(ncpw%ngw,*)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: htvel(*), htfor(*)

    INTEGER                                  :: i, ia, is
    LOGICAL                                  :: skipanic
    REAL(real_8)                             :: alfap

! Variables
! 
! NN: IF BS AND HS-WF SKIP SCALING OF IONIC AND CELL VELOCITIES 

    skipanic=.FALSE.
    IF (cntl%bsymm.AND.(bsclcs.EQ.2))skipanic=.TRUE.
    ! 
    ! ==--------------------------------------------------------------==
    ! ==  DAMPING (IONS)  Fixed friction factor                       ==
    ! ==--------------------------------------------------------------==
    IF (cntl%dampi.AND.(.NOT.skipanic)) THEN
       alfap=cntr%dampgi
#if defined(__VECTOR)
       !$omp parallel do private(I,IS,IA)
#else
       !$omp parallel do private(I,IS,IA) schedule(static)
#endif
       DO i=1,ions1%nat
          ia=iatpt(1,i)
          is=iatpt(2,i)
          fion(1,ia,is)=fion(1,ia,is)-alfap*velp(1,ia,is)
          fion(2,ia,is)=fion(2,ia,is)-alfap*velp(2,ia,is)
          fion(3,ia,is)=fion(3,ia,is)-alfap*velp(3,ia,is)
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ==  DAMPING (ELECTRONS)  Fixed friction factor                  ==
    ! ==--------------------------------------------------------------==
    IF (cntl%dampe) THEN
       alfap=-cntr%dampge
       CALL daxpy(2*ncpw%ngw*nstate,alfap,cm,1,c2,1)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ==  DAMPING (CELL)  Fixed friction factor                       ==
    ! ==--------------------------------------------------------------==
    IF (cntl%dampc.AND.(.NOT.skipanic)) THEN
       alfap=cntr%dampgc
       DO i=1,9
          htfor(i)=htfor(i)-alfap*htvel(i)
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE dampdyn
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE berendsen(velp,cm,nstate,htvel,ekinc,ekinh)
    ! NOTE: we get called twice for each half step of the velocity 
    ! verlet so the the scaling is only considering half a timestep.
    ! ==--------------------------------------------------------------==

    REAL(real_8)                             :: velp(:,:,:)
    COMPLEX(real_8)                          :: cm(ncpw%ngw,*)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: htvel(*), ekinc, ekinh

    INTEGER                                  :: i, ia, is
    LOGICAL                                  :: skipanic
    REAL(real_8)                             :: ekintmp, lambda, tempcur, &
                                                thresh

    PARAMETER(thresh=1.0e-6_real_8)  ! don't scale if we are too "cold".
    ! 
    ! NN: IF BS AND HS-WF SKIP SCALING OF IONIC AND CELL VELOCITIES 
    skipanic=.FALSE.
    IF (cntl%bsymm.AND.(bsclcs.EQ.2))skipanic=.TRUE.
    ! ==--------------------------------------------------------------==
    ! ==  BERENDSEN THERMOSTAT (IONS)                                 ==
    ! ==--------------------------------------------------------------==
    IF (paral%parent.AND.cntl%tberp.AND.(.NOT.skipanic)) THEN
       IF (cntl%tpath.AND.cntl%tpimd) THEN
          IF (pimd1%tpinm.OR.pimd1%tstage) THEN
             CALL s_ekinpp(ekintmp,velp,ipcurr)
          ELSE
             CALL ekinpp(ekintmp,velp)
          ENDIF
       ELSE
          CALL ekinpp(ekintmp,velp)
       ENDIF
       tempcur=ekintmp*factem*2.0_real_8/glib
       IF (tempcur/cntr%tempw.GT.thresh) THEN
          lambda=SQRT(1.0_real_8+0.5_real_8*dt_ions/cntr%taubp*(cntr%tempw/tempcur-1.0_real_8))
#if defined(__VECTOR)
          !$omp parallel do private(I,IS,IA)
#else
          !$omp parallel do private(I,IS,IA) schedule(static)
#endif
          DO i=1,ions1%nat
             ia=iatpt(1,i)
             is=iatpt(2,i)
             velp(1,ia,is)=lambda*velp(1,ia,is)
             velp(2,ia,is)=lambda*velp(2,ia,is)
             velp(3,ia,is)=lambda*velp(3,ia,is)
          ENDDO
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ==  BERENDSEN THERMOSTAT (ELECTRONS)                            ==
    ! ==--------------------------------------------------------------==
    IF (cntl%tbere) THEN
       IF (ekinc/cntr%ekinw.GT.thresh)THEN
          lambda=SQRT(1.0_real_8+0.5_real_8*dt_ions/cntr%taube*(cntr%ekinw/ekinc-1.0_real_8))
          CALL dscal(2*ncpw%ngw*nstate,lambda,cm,1)
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ==  BERENDSEN THERMOSTAT (CELL)                                 ==
    ! ==--------------------------------------------------------------==
    IF (paral%parent.AND.cntl%tberc.AND.(.NOT.skipanic)) THEN
       IF (ekinh/cntr%ekinhr.GT.thresh)THEN
          lambda=SQRT(1.0_real_8+0.5_real_8*dt_ions/cntr%taubc*(cntr%ekinhr/ekinh-1.0_real_8))
          DO i=1,9
             htvel(i)=lambda*htvel(i)
          ENDDO
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE berendsen
  ! ==================================================================
  SUBROUTINE tempramp(temp1,temp2)

    REAL(real_8)                             :: temp1, temp2

! ==--------------------------------------------------------------==

    IF (paral%parent.AND.cntr%trampr.GT.1.0e-30_real_8) THEN
       ! ==--------------------------------------------------------------==
       IF (ABS(cntr%tempw-cntr%trampt).GT.cntr%trampr) THEN
          IF (cntr%trampt.GT.cntr%tempw) THEN
             cntr%tempw=cntr%tempw+dt_ions*cntr%trampr
          ELSE
             cntr%tempw=cntr%tempw-dt_ions*cntr%trampr
          ENDIF
       ELSE
          cntr%tempw=cntr%trampt
          IF (paral%io_parent)&
               WRITE(6,*) 'TEMPRAMP| TEMPERATURE RAMP COMPLETED AT: ',cntr%tempw
          cntr%trampr=0.0_real_8
       ENDIF
       temp1=cntr%tempw+cntr%tolp
       temp2=cntr%tempw-cntr%tolp
       ! ==--------------------------------------------------------------==
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE tempramp
  ! ==================================================================

END MODULE anneal_utils
