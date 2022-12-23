MODULE rinvel_utils
  USE cnst,                            ONLY: factem,&
                                             pi
  USE cnst_dyn,                        ONLY: lmeta
  USE coor,                            ONLY: lvelini
  USE cotr,                            ONLY: lskcor
  USE ekinpp_utils,                    ONLY: ekinpp,&
                                             s_ekinpp
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE metr,                            ONLY: metr_com
  USE mm_input,                        ONLY: g96_vel
  USE mp_interface,                    ONLY: mp_bcast
  USE nose,                            ONLY: glib,&
                                             tcafes,&
                                             tempwr
  USE parac,                           ONLY: parai,&
                                             paral
  USE pimd,                            ONLY: pma0s,&
                                             pimd1,&
                                             ipcurr
  USE prng_utils,                      ONLY: repprngu,&
                                             repprngu_vec
  USE puttau_utils,                    ONLY: taucl
  USE rmas,                            ONLY: rmass
  USE system,                          ONLY: cntl,&
                                             cntr,&
                                             maxsys,&
                                             ncpw,&
                                             nkpt
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rinvel
  PUBLIC :: s_rinvel
  PUBLIC :: rvscal

CONTAINS

  ! ==================================================================
  SUBROUTINE rinvel(vel,cm,nstate)
    ! ==--------------------------------------------------------------==

    REAL(real_8)                             :: vel(:,:,:)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: cm(nkpt%ngwk,nstate)

    REAL(real_8), PARAMETER :: one3rd = 1.0_real_8/3.0_real_8 

    INTEGER                                  :: ia, iaa, is
    LOGICAL                                  :: lveltrue
    REAL(real_8)                             :: alfa1, alfa2, alfa3, ekinp, &
                                                pma00, rnr(3), sigma, tempp, &
                                                tscal, vcm(3), vscale

! ==--------------------------------------------------------------==
! Initialize Velocities for Wavefunctions

    IF (.NOT.cntl%tdiag) CALL zeroing(cm)!,nstate*nkpt%ngwk)
    ! Initialize Cell velocities
    CALL zeroing(metr_com%htvel)!,9)
    ! check for velocities from input
    lveltrue=.FALSE.
    DO is=1,ions1%nsp
       lveltrue=lveltrue.OR.lvelini(0,is)
    ENDDO

    ! MAXWELL DISTRIBUTION FOR THE IONS
    IF (paral%parent) THEN
       IF (lveltrue) THEN
          IF (paral%io_parent)&
               WRITE(6,*) 'RINVEL| INITIAL VELOCITIES TAKEN FROM INPUT'
          CALL taucl(vel)
          GOTO 100
       ENDIF
       IF ( g96_vel%ntx_vel.EQ.1 ) THEN
          IF (paral%io_parent)&
               WRITE(6,*) 'RINVEL| INITIAL VELOCITIES ',&
               'ALREADY TAKEN FROM GROMOS INPUT'
          CALL taucl(vel)
          CALL ekinpp(ekinp,vel)
          tempp=ekinp*factem*2._real_8/glib
          IF (paral%io_parent)&
               WRITE(6,'(1x,''RINVEL| (TEMPERATURE IS '',F10.2,'' K)'')')&
               TEMPP
          GOTO 100
       ENDIF
       IF (tcafes) THEN
          IF (paral%io_parent)&
               WRITE(6,*)'RINVEL| CAFES initial temperature distribution'
       ENDIF
       CALL zeroing(vel)!,3*maxsys%nax*maxsys%nsx)
       IF (cntr%tempw.LT.1.e-5_real_8) GOTO 100
       CALL repprngu_vec(3,rnr)
       ! ---------------------------------------------------------------------------
       ! CC  Randomize velocities for Saddle to Min dynamics 
       IF (lmeta%lsadpnt) THEN

          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                IF (tcafes) THEN
                   sigma=SQRT(tempwr(1,ia,is)/(rmass%pma(is)*factem))
                ELSE
                   sigma=SQRT(cntr%tempw/(rmass%pma(is)*factem))
                ENDIF

                CALL repprngu_vec(3,rnr)
                alfa1=2.0_real_8*pi*rnr(1)
                alfa2=2.0_real_8*pi*rnr(2)
                alfa3=2.0_real_8*pi*rnr(3)
                IF (.NOT.lvelini(ia,is)) THEN
                   vel(1,ia,is)&
                        =SQRT(LOG(repprngu())*(-2.0_real_8))*COS(alfa1)*sigma
                   vel(2,ia,is)&
                        =SQRT(LOG(repprngu())*(-2.0_real_8))*COS(alfa2)*sigma
                   vel(3,ia,is)&
                        =SQRT(LOG(repprngu())*(-2.0_real_8))*COS(alfa3)*sigma
                ENDIF
             ENDDO
          ENDDO
          ! ---------------------------------------------------------------------------
       ELSE

          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                IF (tcafes) THEN
                   sigma=SQRT(tempwr(1,ia,is)/(rmass%pma(is)*factem))
                ELSE
                   sigma=SQRT(cntr%tempw/(rmass%pma(is)*factem))
                ENDIF
                CALL repprngu_vec(3,rnr)
                alfa1=2.0_real_8*pi*rnr(1)
                alfa2=2.0_real_8*pi*rnr(2)
                alfa3=2.0_real_8*pi*rnr(3)
                IF (.NOT.lvelini(ia,is)) THEN
                   vel(1,ia,is)&
                        =SQRT(LOG(repprngu())*(-2.0_real_8))*COS(alfa1)*sigma
                   vel(2,ia,is)&
                        =SQRT(LOG(repprngu())*(-2.0_real_8))*COS(alfa2)*sigma
                   vel(3,ia,is)&
                        =SQRT(LOG(repprngu())*(-2.0_real_8))*COS(alfa3)*sigma
                ENDIF
             ENDDO
          ENDDO
       ENDIF
       ! ---------------------------------------------------------------------------
       CALL taucl(vel)
       ! SUBTRACT CENTER OF MASS VELOCITY
       vcm(1)=0.0_real_8
       vcm(2)=0.0_real_8
       vcm(3)=0.0_real_8
       pma00=0._real_8
       iaa=0
       lveltrue=.FALSE.
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             iaa=iaa+1
             IF (.NOT.lvelini(ia,is)) THEN
                vcm(1)=vcm(1)+vel(1,ia,is)*rmass%pma0(is)
                vcm(2)=vcm(2)+vel(2,ia,is)*rmass%pma0(is)
                vcm(3)=vcm(3)+vel(3,ia,is)*rmass%pma0(is)
                pma00=pma00+one3rd*(lskcor(1,iaa)+&
                     LSKCOR(2,IAA)+LSKCOR(3,IAA))*rmass%pma0(IS)
             ELSE
                lveltrue=.TRUE.
             ENDIF
          ENDDO
       ENDDO
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             IF (.NOT.lvelini(ia,is)) THEN
                vel(1,ia,is)=vel(1,ia,is)-vcm(1)/pma00
                vel(2,ia,is)=vel(2,ia,is)-vcm(2)/pma00
                vel(3,ia,is)=vel(3,ia,is)-vcm(3)/pma00
             ENDIF
          ENDDO
       ENDDO
       CALL taucl(vel)
       ! RESCALE VELOCITIES
       IF (.NOT.tcafes) THEN
          CALL ekinpp(ekinp,vel)
          tempp=ekinp*factem*2._real_8/glib
          IF (tempp.GT.1.e-5_real_8) THEN
             tscal=cntr%tempw/tempp
             vscale=SQRT(tscal)
             DO is=1,ions1%nsp
                DO ia=1,ions0%na(is)
                   IF (.NOT.lvelini(ia,is)) THEN
                      vel(1,ia,is)=vel(1,ia,is)*vscale
                      vel(2,ia,is)=vel(2,ia,is)*vscale
                      vel(3,ia,is)=vel(3,ia,is)*vscale
                   ENDIF
                ENDDO
             ENDDO
          ENDIF
       ENDIF
100    CONTINUE
    ENDIF

    ! check for remaining C.O.M. velocity
    pma00=0.0_real_8
    vcm(1)=0.0_real_8
    vcm(2)=0.0_real_8
    vcm(3)=0.0_real_8
    iaa=0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          iaa=iaa+1
          vcm(1)=vcm(1)+vel(1,ia,is)*rmass%pma0(is)
          vcm(2)=vcm(2)+vel(2,ia,is)*rmass%pma0(is)
          vcm(3)=vcm(3)+vel(3,ia,is)*rmass%pma0(is)
          pma00=pma00+one3rd*(lskcor(1,iaa)+&
               LSKCOR(2,IAA)+LSKCOR(3,IAA))*rmass%pma0(IS)
       ENDDO
    ENDDO
    ! CHECK ONLY IF THERE ARE MOVABLE ATOMS...
    IF (pma00.GT.1.0e-10_real_8) THEN
       IF (SQRT(vcm(1)**2+vcm(2)**2+vcm(3)**2)/pma00.GT.1.0e-10_real_8) THEN
          IF ( paral%parent ) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(1X,64("!"))')
             IF (paral%io_parent)&
                  WRITE(6,'(" !! RINVEL| WARNING! ",A,T64,"!!")')&
                  'CENTER OF MASS VELOCITY NOT ZERO'
             IF (paral%io_parent)&
                  WRITE(6,'(" !! ",3F15.10,T64,"!!")')&
                  vcm(1)/pma00,vcm(2)/pma00,vcm(3)/pma00
             IF (paral%io_parent)&
                  WRITE(6,'(1X,64("!"))')
          ENDIF
       ENDIF
    ENDIF

    CALL mp_bcast(vel,3*maxsys%nax*maxsys%nsx,parai%source,parai%allgrp)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rinvel
  ! ==================================================================
  SUBROUTINE s_rinvel(velp,cm,nstate,ip)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: velp(:,:,:)
    COMPLEX(real_8)                          :: cm(*)
    INTEGER                                  :: nstate, ip

    INTEGER                                  :: ia, is
    REAL(real_8)                             :: alfa1, alfa2, alfa3, ekinp, &
                                                glib_s, ptotm, rnr(3), sigma, &
                                                tempp, tscal, vcm(3), vscale

! Variables
! ==--------------------------------------------------------------==
! ..Initialize Velocities for Wavefunctions

    CALL zeroing(cm(1:nstate*ncpw%ngw))!,nstate*ngw)
    ! ..Initialize Cell velocities
    CALL zeroing(metr_com%htvel)!,9)
    ! ..MAXWELL DISTRIBUTION FOR THE IONS
    IF (paral%parent) THEN
       CALL zeroing(velp)!,3*maxsys%nax*maxsys%nsx)
       IF (cntr%tempw.LT.1.e-5_real_8) GOTO 100
       CALL repprngu_vec(3,rnr)
       DO is=1,ions1%nsp
          sigma=SQRT(cntr%tempw/(pma0s(is,ip)*factem))
          DO ia=1,ions0%na(is)
             CALL repprngu_vec(3,rnr)
             alfa1=2.0_real_8*pi*rnr(1)
             alfa2=2.0_real_8*pi*rnr(2)
             alfa3=2.0_real_8*pi*rnr(3)
             velp(1,ia,is)=SQRT(LOG(repprngu())*(-2.0_real_8))*COS(alfa1)*sigma
             velp(2,ia,is)=SQRT(LOG(repprngu())*(-2.0_real_8))*COS(alfa2)*sigma
             velp(3,ia,is)=SQRT(LOG(repprngu())*(-2.0_real_8))*COS(alfa3)*sigma
          ENDDO
       ENDDO
       IF (ip.EQ.1) CALL taucl(velp)
       ! ..SUBTRACT CENTER OF MASS VELOCITY
       vcm(1)=0.0_real_8
       vcm(2)=0.0_real_8
       vcm(3)=0.0_real_8
       ptotm=0.0_real_8
       DO is=1,ions1%nsp
          ptotm=ptotm+pma0s(is,ip)*ions0%na(is)
          DO ia=1,ions0%na(is)
             vcm(1)=vcm(1)+velp(1,ia,is)*pma0s(is,ip)
             vcm(2)=vcm(2)+velp(2,ia,is)*pma0s(is,ip)
             vcm(3)=vcm(3)+velp(3,ia,is)*pma0s(is,ip)
          ENDDO
       ENDDO
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             velp(1,ia,is)=velp(1,ia,is)-vcm(1)/ptotm
             velp(2,ia,is)=velp(2,ia,is)-vcm(2)/ptotm
             velp(3,ia,is)=velp(3,ia,is)-vcm(3)/ptotm
          ENDDO
       ENDDO
       IF (ip.EQ.1) CALL taucl(velp)
       ! ..RESCALE VELOCITIES
       CALL s_ekinpp(ekinp,velp,ip)
       glib_s=3._real_8*ions1%nat
       IF (ip.EQ.1) glib_s=glib
       tempp=ekinp*factem*2._real_8/glib_s
       IF (tempp.GT.1.e-5_real_8) THEN
          tscal=cntr%tempw/tempp
          vscale=SQRT(tscal)
          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                velp(1,ia,is)=velp(1,ia,is)*vscale
                velp(2,ia,is)=velp(2,ia,is)*vscale
                velp(3,ia,is)=velp(3,ia,is)*vscale
             ENDDO
          ENDDO
       ENDIF
100    CONTINUE
    ENDIF
    CALL mp_bcast(velp,3*maxsys%nax*maxsys%nsx,parai%source,parai%allgrp)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE s_rinvel
  ! ==================================================================
  SUBROUTINE rvscal(vel)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: vel(:,:,:)

    INTEGER                                  :: ia, is
    REAL(real_8)                             :: ekinp, tempp, tscal, vscale

    IF (paral%parent) THEN
       ! RESCALE VELOCITIES
       IF (.NOT.tcafes) THEN
          IF (cntl%tpath.AND.cntl%tpimd.AND.(pimd1%tpinm.OR.pimd1%tstage)) THEN
             CALL s_ekinpp(ekinp,vel,ipcurr)
          ELSE
             CALL ekinpp(ekinp,vel)
          ENDIF
          tempp=ekinp*factem*2._real_8/glib
          IF (paral%io_parent)&
               WRITE(6,'(A,F13.5,A,F13.5)') ' RVSCAL| RESCALING IONIC '&
               // 'TEMP FROM ',tempp, ' TO ',cntr%tempw
          IF (tempp.GT.1.e-5_real_8) THEN
             tscal=cntr%tempw/tempp
             vscale=SQRT(tscal)
             DO is=1,ions1%nsp
                DO ia=1,ions0%na(is)
                   vel(1,ia,is)=vel(1,ia,is)*vscale
                   vel(2,ia,is)=vel(2,ia,is)*vscale
                   vel(3,ia,is)=vel(3,ia,is)*vscale
                ENDDO
             ENDDO
          ENDIF
       ENDIF
    ENDIF
    CALL mp_bcast(vel,3*maxsys%nax*maxsys%nsx,parai%source,parai%allgrp)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rvscal
  ! ==================================================================

END MODULE rinvel_utils
