MODULE posupa_utils
  USE dotp_utils,                      ONLY: dotp
  USE error_handling,                  ONLY: stopgm
  USE harm,                            ONLY: dcgt,&
                                             dsgtw,&
                                             dtan2w,&
                                             wdsgt,&
                                             xmu
  USE jrotation_utils,                 ONLY: set_orbdist
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai
  USE pslo,                            ONLY: pslo_com
  USE ropt,                            ONLY: ropt_mod
  USE rortog_utils,                    ONLY: give_scr_rortog,&
                                             rortog
  USE rotate_utils,                    ONLY: rotate
  USE spin,                            ONLY: spin_mod
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             cntr,&
                                             ncpw,&
                                             parap,&
                                             paraw
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE tpar,                            ONLY: dt2hbe,&
                                             dt_elec,&
                                             dtb2me
!!use rotate_utils, only : rotate
!!se rotate_utils, only : rotate_da

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: posupa
  PUBLIC :: give_scr_posupa

CONTAINS

  ! ==================================================================
  SUBROUTINE posupa(c0,cm,c2,gamx,nstate)
    ! ==--------------------------------------------------------------==
    ! ==  UPDATE OF THE POSITIONS FOR VELOCITY VERLET                 ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c2(ncpw%ngw,*)
    REAL(real_8)                             :: gamx(*)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: cm(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)

    INTEGER                                  :: i, ig, isub, nstx
    LOGICAL                                  :: prtev, tnon
    REAL(real_8)                             :: ai(3), odt, pf1, pf2, pf3, &
                                                pf4, xi

! Variables
! ==--------------------------------------------------------------==
! ==  WAVEFUNCTIONS                                               ==
! ==--------------------------------------------------------------==

    CALL tiset('    POSUPA',isub)
    IF (tkpts%tkpnt) CALL stopgm('POSUPA','K POINTS NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    IF (cntl%nonort) THEN
       IF (cntl%tharm) THEN
          !$omp parallel do private (I,IG,PF1,PF2,PF3)
          DO i=1,nstate
             DO ig=1,ncpw%ngw
                pf1=dcgt(ig)
                pf2=dsgtw(ig)
                pf3=wdsgt(ig)
                c2(ig,i)=pf1*c0(ig,i)+pf2*cm(ig,i)
                cm(ig,i)=pf1*cm(ig,i)-pf3*c0(ig,i)
             ENDDO
          ENDDO
       ELSE
          !$omp parallel do private (I,IG)
          DO i=1,nstate
             DO ig=1,ncpw%ngw
                c2(ig,i)=c0(ig,i)+dt_elec*cm(ig,i)
             ENDDO
          ENDDO
       ENDIF
       ! ..LINEAR CONSTRAINTS
       IF (cntl%tharm.OR.cntl%tmass) THEN
          DO i=1,nstate
             !$omp parallel do private (IG,PF4)
             DO ig=1,ncpw%ngw
                IF (cntl%tharm) THEN
                   pf4=dsgtw(ig)*dtan2w(ig)/dtb2me
                ELSE
                   pf4=cntr%emass/xmu(ig)
                ENDIF
                c0(ig,i)=pf4*c0(ig,i)
             ENDDO
             ai(1)=dotp(ncpw%ngw,c0(:,i),c0(:,i))
             ai(2)=2.0_real_8*dotp(ncpw%ngw,c2(:,i),c0(:,i))
             ai(3)=dotp(ncpw%ngw,c2(:,i),c2(:,i))
             CALL mp_sum(ai,3,parai%allgrp)
             ai(3)=ai(3)-1.0_real_8
             xi=(-ai(2)+SQRT(ai(2)*ai(2)-4._real_8*ai(1)*ai(3)))/(2._real_8*ai(1))
             CALL daxpy(2*ncpw%ngw,xi,c0(1,i),1,c2(1,i),1)
             CALL daxpy(2*ncpw%ngw,xi/dt_elec,c0(1,i),1,cm(1,i),1)
          ENDDO
       ELSE
          DO i=1,nstate
             ai(1)=dotp(ncpw%ngw,c0(:,i),c0(:,i))
             ai(2)=2.0_real_8*dotp(ncpw%ngw,c2(:,i),c0(:,i))
             ai(3)=dotp(ncpw%ngw,c2(:,i),c2(:,i))
             CALL mp_sum(ai,3,parai%allgrp)
             ai(3)=ai(3)-1.0_real_8
             xi=(-ai(2)+SQRT(ai(2)*ai(2)-4._real_8*ai(1)*ai(3)))/&
                  (2._real_8*ai(1))
             CALL daxpy(2*ncpw%ngw,xi,c0(1,i),1,c2(1,i),1)
             CALL daxpy(2*ncpw%ngw,xi/dt_elec,c0(1,i),1,cm(1,i),1)
          ENDDO
       ENDIF
       CALL dcopy(2*ncpw%ngw*nstate,c2(1,1),1,c0(1,1),1)
    ELSE
       IF (cntl%tharm) THEN
          !$omp parallel do private (I,IG,PF1,PF2,PF3,PF4)
          DO i=1,nstate
             DO ig=1,ncpw%ngw
                pf1=dcgt(ig)
                pf2=dsgtw(ig)
                pf3=wdsgt(ig)
                pf4=pf2*dtan2w(ig)/dt2hbe
                c2(ig,i)=pf1*c0(ig,i)+pf2*cm(ig,i)
                cm(ig,i)=pf1*cm(ig,i)-pf3*c0(ig,i)
                c0(ig,i)=pf4*c0(ig,i)
             ENDDO
          ENDDO
       ELSE
          !$omp parallel do private (I,IG)
          DO i=1,nstate
             DO ig=1,ncpw%ngw
                c2(ig,i)=c0(ig,i)+dt_elec*cm(ig,i)
             ENDDO
          ENDDO
          IF (cntl%tmass) THEN
             !$omp parallel do private (I,IG,PF4)
             DO i=1,nstate
                DO ig=1,ncpw%ngw
                   pf4=cntr%emass/xmu(ig)
                   c0(ig,i)=pf4*c0(ig,i)
                ENDDO
             ENDDO
          ENDIF
       ENDIF
       ! REORTHOGONALIZE THE WAVEFUNCTIONS
       prtev=ropt_mod%prteig .AND. .NOT.pslo_com%tivan
       tnon=cntl%tmass.OR.cntl%tharm
       CALL rortog(c0,c2,gamx,tnon,nstate,prtev)
       ! ADD CONSTRAINTS FOR WAVEFUNCTION VELOCITIES 
       IF (cntl%tharm) THEN
          !$omp parallel do private (I,IG,PF1,PF2,PF4)
          DO i=1,nstate
             DO ig=1,ncpw%ngw
                pf1=dcgt(ig)
                pf2=dsgtw(ig)
                pf4=cntr%delt_elec*(pf1/pf2)
                c0(ig,i)=pf4*c0(ig,i)
             ENDDO
          ENDDO
       ENDIF
       odt=1.0_real_8/dt_elec
       IF (ncpw%ngw.GT.0) THEN
          IF (cntl%tdmal) THEN
             CALL set_orbdist(nstate,cnti%nstblk,parai%nproc,nstx)
             CALL rotate_da(odt,c0,1._real_8,cm,gamx,2*ncpw%ngw,2*ncpw%ngw,nstate,&
                  paraw%nwa12(0,1),paraw%nwa12(0,2),nstx,parai%mepos,parap%pgroup,parai%nproc,&
                  parai%allgrp,cntl%tlsd,spin_mod%nsup,spin_mod%nsdown)
          ELSE
             CALL rotate(odt,c0,1.0_real_8,cm,gamx,nstate,2*ncpw%ngw,cntl%tlsd,&
                  spin_mod%nsup,spin_mod%nsdown)
          ENDIF
          CALL dcopy(2*nstate*ncpw%ngw,c2(1,1),1,c0(1,1),1)
       ENDIF
    ENDIF
    CALL tihalt('    POSUPA',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE posupa
  ! ==================================================================
  SUBROUTINE give_scr_posupa(lposupa,tag,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lposupa
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate

    INTEGER                                  :: lortog, lrotate, nstx

! ==--------------------------------------------------------------==

    IF (cntl%nonort) THEN
       lposupa=0
    ELSE
       IF (cntl%tdmal) THEN
          CALL set_orbdist(nstate,cnti%nstblk,parai%nproc,nstx)
          lrotate=nstate*nstx
       ELSE
          lrotate=0
       ENDIF
       CALL give_scr_rortog(lortog,tag,nstate)
       lposupa=MAX(lortog,lrotate)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_posupa
  ! ==================================================================

END MODULE posupa_utils
