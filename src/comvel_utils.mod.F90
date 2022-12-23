MODULE comvel_utils
  USE cnst,                            ONLY: factem
  USE cotr,                            ONLY: lskcor
  USE ekinpp_utils,                    ONLY: ekinpp
  USE ions,                            ONLY: ions1
  USE kinds,                           ONLY: real_8
  USE nose,                            ONLY: glib
  USE puttau_utils,                    ONLY: taucl
  USE rmas,                            ONLY: rmass
  USE system,                          ONLY: iatpt

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: comvel

CONTAINS

  ! ==================================================================
  SUBROUTINE comvel(velp,vcmio,inout)
    ! REMOVE TRANSLATIONAL VELOCITIES. 
    ! Last Update: 2008/06/11 by AK
    ! To integrate more cleanly with the logic of the velocity verlet
    ! algorithm you have to call this subroutine twice. First with 
    ! INOUT set to .TRUE., where the routine will store the center of
    ! mass velocity vector and the total mass in VCMIO(4). This call 
    ! will remove only half of the velocity and has to be called before 
    ! the first (half) velocity update step. For the second call, which 
    ! has to come after the second (half) velocity update step, with 
    ! INOUT set to .FALSE. those data is then used to remove the 
    ! remainder of the center of mass motion.
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: velp(:,:,:), vcmio(4)
    LOGICAL                                  :: inout

    REAL(real_8), PARAMETER :: one3rd = 1.0_real_8/3.0_real_8 

    INTEGER                                  :: ia, iat, is
    REAL(real_8)                             :: ekinp1, ekinp2, pma00, &
                                                tempp1, tempp2, tscal, &
                                                vcm(3), vscale

! true  -> write VCMIO
! false -> read  VCMIO
! Variables
! ==--------------------------------------------------------------==
! COMPUTE THE IONIC TEMPERATURE TEMPP1

    CALL ekinpp(ekinp1,velp)
    tempp1=ekinp1*factem*2.0_real_8/glib
    IF (inout) THEN
       ! DETERMINE CENTER OF MASS VELOCITY
       vcm(1)=0.0_real_8
       vcm(2)=0.0_real_8
       vcm(3)=0.0_real_8
       pma00=0._real_8
       !$omp parallel do private(IS,IA,IAT) reduction(+:VCM,PMA00) &
       !$omp  schedule(static)
       DO iat=1,ions1%nat
          ia=iatpt(1,iat)
          is=iatpt(2,iat)
          vcm(1)=vcm(1)+velp(1,ia,is)*rmass%pma0(is)
          vcm(2)=vcm(2)+velp(2,ia,is)*rmass%pma0(is)
          vcm(3)=vcm(3)+velp(3,ia,is)*rmass%pma0(is)
          ! !        PMA00=PMA00+PMA0(IS)
          pma00=pma00+one3rd*(lskcor(1,iat)+&
               LSKCOR(2,IAT)+LSKCOR(3,IAT))*rmass%pma0(IS)
       ENDDO
       vcmio(1)=vcm(1)
       vcmio(2)=vcm(2)
       vcmio(3)=vcm(3)
       vcmio(4)=pma00
    ELSE
       ! TAKE CENTER OF MASS VELOCITY FROM LAST CALL
       vcm(1)=vcmio(1)
       vcm(2)=vcmio(2)
       vcm(3)=vcmio(3)
       pma00=vcmio(4)
    ENDIF
    !$omp parallel do private(IS,IA,IAT) shared(VCM,PMA00) &
    !$omp  schedule(static)
    DO iat=1,ions1%nat
       ia=iatpt(1,iat)
       is=iatpt(2,iat)
       velp(1,ia,is)=velp(1,ia,is)-0.5_real_8*vcm(1)/pma00
       velp(2,ia,is)=velp(2,ia,is)-0.5_real_8*vcm(2)/pma00
       velp(3,ia,is)=velp(3,ia,is)-0.5_real_8*vcm(3)/pma00
    ENDDO
    CALL taucl(velp)
    ! RESCALE ALL VELOCITIES TO COMPENSATE FOR LOSS OF KINETIC ENERGY
    CALL ekinpp(ekinp2,velp)
    tempp2=ekinp2*factem*2.0_real_8/glib
    IF (tempp2.GT.1.e-5_real_8) THEN
       tscal=tempp1/tempp2
       vscale=SQRT(tscal)
       !$omp parallel do private(IS,IA,IAT) shared(VSCALE) &
       !$omp  schedule(static)
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
  END SUBROUTINE comvel
  ! ==================================================================

END MODULE comvel_utils
