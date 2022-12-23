MODULE rotvel_utils
  USE cnst,                            ONLY: factem
  USE ekinpp_utils,                    ONLY: ekinpp
  USE ions,                            ONLY: ions1
  USE kinds,                           ONLY: real_8
  USE nose,                            ONLY: glib
  USE rmas,                            ONLY: rmass
  USE system,                          ONLY: iatpt
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rotvel

CONTAINS

  ! ==================================================================
  ! REMOVE ROTATIONAL VELOCITIES. 
  ! Last Update: 2008/06/11 by AK
  ! To integrate more cleanly with the logic of the velocity verlet
  ! algorithm you have to call this subroutine twice. First with 
  ! INOUT set to .TRUE., where the routine will store the angular
  ! rotation velocity vector in LM(3) and position vectors of the atoms
  ! relative to the center of mass in TSCR(3,NAT). This call will remove
  ! only half of the velocity and has to be called before the first 
  ! (half) velocity update step. For the second call, which has to come
  ! after the second (half) velocity update step, with INOUT set to
  ! .FALSE. those data is used to remove the remainder of the 
  ! rotational motion.
  ! ==================================================================
  SUBROUTINE rotvel(tau0,velp,lm,tscr,inout)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:), velp(:,:,:), &
                                                lm(3), tscr(3,*)
    LOGICAL                                  :: inout

    INTEGER                                  :: ia, iat, ierr, is
    REAL(real_8) :: cm(3), ekinp1, ekinp2, evl(3), evm(3,3), lten(3,3), &
      pscr(ions1%nat), pscrt, r0, r1, r2, r3, small, tempp1, tempp2, tscal, &
      vscale, vscr(3,ions1%nat), vtmp1, vtmp2, vtmp3, wk(9), wx(3)

! true  -> write LM(3)/TSCR(NAT)
! false -> read  LM(3)/TSCR(NAT)
! Variables

    PARAMETER (small=1.0e-10_real_8)  ! small angular momentum where we stop considering it
    ! ==--------------------------------------------------------------==
    ! COMPUTE THE IONIC TEMPERATURE TEMPP1
    CALL ekinpp(ekinp1,velp)
    tempp1=ekinp1*factem*2.0_real_8/glib

    IF (inout) THEN
       ! GET CENTER OF MASS AND COPY TO HELPER ARRAY
       cm(1)=0.0_real_8
       cm(2)=0.0_real_8
       cm(3)=0.0_real_8
       pscrt=0.0_real_8
       !$omp parallel do private(IAT,IA,IS,R0,R1,R2,R3) reduction(+:PSCRT,CM) 
       DO iat=1,ions1%nat
          ia=iatpt(1,iat)
          is=iatpt(2,iat)
          r0=rmass%pma0(is)
          pscrt=pscrt+r0
          pscr(iat)=r0
          r1=tau0(1,ia,is)
          tscr(1,iat)=r1
          cm(1)=cm(1)+r1*r0
          r2=tau0(2,ia,is)
          tscr(2,iat)=r2
          cm(2)=cm(2)+r2*r0
          r3=tau0(3,ia,is)
          tscr(3,iat)=r3
          cm(3)=cm(3)+r3*r0
          vscr(1,iat)=velp(1,ia,is)
          vscr(2,iat)=velp(2,ia,is)
          vscr(3,iat)=velp(3,ia,is)
       ENDDO
       cm(1)=cm(1)/pscrt
       cm(2)=cm(2)/pscrt
       cm(3)=cm(3)/pscrt

       ! translate C.O.M to origin.
       !$omp parallel do private(IAT)
       DO iat=1,ions1%nat
          tscr(1,iat)=tscr(1,iat)-cm(1)
          tscr(2,iat)=tscr(2,iat)-cm(2)
          tscr(3,iat)=tscr(3,iat)-cm(3)
       ENDDO

       ! get total angular momentum
       lm(1)=0.0_real_8
       lm(2)=0.0_real_8
       lm(3)=0.0_real_8
       !$omp parallel do private(IAT,VTMP1,VTMP2,VTMP3,R1) &
       !$omp  reduction(+:LM)
       DO iat=1,ions1%nat
          vtmp1=tscr(2,iat)*vscr(3,iat)-tscr(3,iat)*vscr(2,iat)
          vtmp2=tscr(3,iat)*vscr(1,iat)-tscr(1,iat)*vscr(3,iat)
          vtmp3=tscr(1,iat)*vscr(2,iat)-tscr(2,iat)*vscr(1,iat)
          r1=pscr(iat)
          lm(1)=lm(1)+vtmp1*r1
          lm(2)=lm(2)+vtmp2*r1
          lm(3)=lm(3)+vtmp3*r1
       ENDDO

       ! get moment of inertia
       CALL zeroing(lten)!,9)
       DO iat=1,ions1%nat
          lten(1,1)=lten(1,1)+pscr(iat)*(tscr(2,iat)**2+tscr(3,iat)**2)
          lten(1,2)=lten(1,2)-pscr(iat)*(tscr(1,iat)*tscr(2,iat))
          lten(1,3)=lten(1,3)-pscr(iat)*(tscr(1,iat)*tscr(3,iat))

          lten(2,1)=lten(2,1)-pscr(iat)*(tscr(1,iat)*tscr(2,iat))
          lten(2,2)=lten(2,2)+pscr(iat)*(tscr(1,iat)**2+tscr(3,iat)**2)
          lten(2,3)=lten(2,3)-pscr(iat)*(tscr(2,iat)*tscr(3,iat))

          lten(3,1)=lten(3,1)-pscr(iat)*(tscr(1,iat)*tscr(3,iat))
          lten(3,2)=lten(3,2)-pscr(iat)*(tscr(2,iat)*tscr(3,iat))
          lten(3,3)=lten(3,3)+pscr(iat)*(tscr(1,iat)**2+tscr(2,iat)**2)
       ENDDO

       ierr=0
       CALL dcopy(9,lten,1,evm,1)
       CALL dsyev('V','U',3,evm,3,evl,wk,9,ierr)

       ! get angular velocity in body frame. ignore if moment of inertia is small.
       CALL dgemv('T',3,3,1.0_real_8,evm,3,lm,1,0.0_real_8,wx,1)
       IF (ABS(evl(1)).GT.small) THEN
          wx(1)=wx(1)/evl(1)
       ELSE
          wx(1)=0.0_real_8
       ENDIF
       IF (ABS(evl(2)).GT.small) THEN
          wx(2)=wx(2)/evl(2)
       ELSE
          wx(2)=0.0_real_8
       ENDIF
       IF (ABS(evl(3)).GT.small) THEN
          wx(3)=wx(3)/evl(3)
       ELSE
          wx(3)=0.0_real_8
       ENDIF
       ! DEBUG
       ! write(6,'(A,3F12.8,3F12.8)') 'ANG VEL in body axis:', WX, EVL
       ! transform back to lab frame and subtract
       CALL dgemv('N',3,3,1.0_real_8,evm,3,wx,1,0.0_real_8,lm,1)
    ENDIF

    !$omp parallel do private(IS,IA,IAT,VTMP1,VTMP2,VTMP3)
    DO iat=1,ions1%nat
       ia=iatpt(1,iat)
       is=iatpt(2,iat)
       vtmp1=lm(2)*tscr(3,iat)-lm(3)*tscr(2,iat)
       vtmp2=lm(3)*tscr(1,iat)-lm(1)*tscr(3,iat)
       vtmp3=lm(1)*tscr(2,iat)-lm(2)*tscr(1,iat)
       velp(1,ia,is)=velp(1,ia,is)-0.5_real_8*vtmp1
       velp(2,ia,is)=velp(2,ia,is)-0.5_real_8*vtmp2
       velp(3,ia,is)=velp(3,ia,is)-0.5_real_8*vtmp3
    ENDDO

    ! COMPUTE THE IONIC TEMPERATURE AFTER THE REMOVAL
    CALL ekinpp(ekinp2,velp)
    tempp2=ekinp2*factem*2.0_real_8/glib

    ! rescale velocities to conserve the total energy
    IF (tempp2.GT.1.e-5_real_8) THEN
       tscal=tempp1/tempp2
       vscale=SQRT(tscal)
       !$omp parallel do private(IS,IA,IAT) shared(VSCALE)
       DO iat=1,ions1%nat
          ia=iatpt(1,iat)
          is=iatpt(2,iat)
          velp(1,ia,is)=velp(1,ia,is)*vscale
          velp(2,ia,is)=velp(2,ia,is)*vscale
          velp(3,ia,is)=velp(3,ia,is)*vscale
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rotvel
  ! ==================================================================

END MODULE rotvel_utils
