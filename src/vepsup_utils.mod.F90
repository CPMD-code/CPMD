MODULE vepsup_utils
  USE cnst,                            ONLY: au_fs
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE metr,                            ONLY: metr_com,&
                                             veps
  USE rmas,                            ONLY: rmass
  USE shock,                           ONLY: shock1
  USE system,                          ONLY: cntr,&
                                             parm
  USE tpar,                            ONLY: dt_ions
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: vepsup
  PUBLIC :: vepsupdamp

CONTAINS

  ! ==================================================================
  SUBROUTINE vepsup(velp)
    ! ==--------------------------------------------------------------==
    ! ==  PROPOGATES THE CELL VELOCITIES FOR NPH-SHOCK DYNAMCIS       ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: velp(:,:,:)

    INTEGER                                  :: ia, is, k, l
    REAL(real_8)                             :: at3, geps_p, geps_v, &
                                                htfv(3,3), v, v0, v0i, vi

    at3 = 3._real_8*REAL(ions1%nat,kind=real_8)
    ! ==--------------------------------------------------------------==
    ! V  = INSTANTANEOUS VOLUME
    ! VOL0 = INITIAL VOLUME DETERMINED 
    ! VI = INVERSE INSTANTANEOUS VOLUME
    ! V0I = INVERSE INITIAL VOLUME
    v = parm%omega
    vi = 1._real_8/v
    v0 = shock1%vol0
    v0i = 1._real_8/v0
    ! ==--------------------------------------------------------------==
    ! ==  GET VELOCITY-DEPENDENT PART OF CELL FORCE                   ==
    ! ==--------------------------------------------------------------==
    ! GEPS_P = HTFP(1,1)+HTFP(2,2)+HTFP(3,3)
    ! GEPS_V=0._real_8
    ! DO IS=1,NSP
    ! DO IA=1,NA(IS)
    ! DO K=1,3
    ! GEPS_V = GEPS_V + (1._real_8+3._real_8/AT3)* 
    ! *               PMA(IS)*VELP(K,IA,IS)*VELP(K,IA,IS)
    ! ENDDO
    ! ENDDO
    ! ENDDO
    geps_p = metr_com%htfp ( 1, 1 )
    geps_v = 0._real_8
    CALL zeroing(htfv)!,9)
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          DO k=1,3
             DO l=1,3
                htfv(k,l) = htfv(k,l) + rmass%pma(is)*velp(k,ia,is)*velp(l,ia,is)
             ENDDO
             geps_v = geps_v + rmass%pma(is)*velp(k,ia,is)*velp(k,ia,is)/at3
          ENDDO
       ENDDO
    ENDDO
    ! CALL AZZERO(HTFV,9)
    ! KIN=0._real_8
    ! DO IS=1,NSP
    ! DO IA=1,NA(IS)
    ! DO K=1,3
    ! DO L=1,3
    ! HTFV(K,L) = HTFV(K,L) + PMA(IS)*VELP(K,IA,IS)*VELP(L,IA,IS)
    ! ENDDO
    ! KIN = KIN +  PMA(IS)*VELP(K,IA,IS)*VELP(K,IA,IS)/AT3
    ! ENDDO
    ! ENDDO
    ! ENDDO  
    ! ==--------------------------------------------------------------==
    ! ==  PROPAGATE BAROSTAT (CELL VELOCITY) ONLY IN "X" DIRECTION    ==
    ! P0 = INITIAL PRESSURE DETERMINED UPON INPUT AND IS 
    ! INCORPORATED INTO HTFP IN SUBROUTINE TOTSTR.F LINE 55
    ! V  = INSTANTANEOUS VOLUME
    ! V0 = INITIAL VOLUME DETERMINED 
    ! VSHOCK = SHOCK SPEED DETERMINED UPON INPUT
    ! ==--------------------------------------------------------------==
    ! BIAS TO THE COMPRESSIVE SOLUTION
    ! IF ( V .LT. V0 ) THEN
    ! HTVEL(1,1) = HTVEL(1,1) +
    ! *     0.5_real_8*DT_IONS*(HTFV(1,1) + HTFP(1,1) + KIN  
    ! *     -VSHOCK*VSHOCK*V*V0I*(1.0_real_8-V*V0I))/CMASS
    ! ELSE
    ! HTVEL(1,1) = HTVEL(1,1) +
    ! *     0.5_real_8*DT_IONS*(HTFV(1,1) + HTFP(1,1) + KIN )/CMASS
    ! ENDIF
    ! VEPS = VEPS+0.5_real_8*DT_IONS*(GEPS_P+GEPS_V)/CMASS
    ! VEPS = VEPS+0.5_real_8*DT_IONS*(GEPS_P+GEPS_V)/CMASS
    IF ( v .LT. v0 ) THEN
       veps = veps +&
            0.5_real_8*dt_ions*(htfv(1,1)+geps_v+geps_p&
            -shock1%vshock*shock1%vshock*v*v0i*(1.0_real_8-v*v0i))/cntr%cmass
    ELSE
       veps = veps +&
            0.5_real_8*dt_ions*(htfv(1,1)+geps_v+geps_p)/cntr%cmass
    ENDIF

    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE vepsup
  ! ==================================================================
  SUBROUTINE vepsupdamp
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: gamma, scale

    scale = 1._real_8
    ! CJM DBG HARDWIRED UNTIL PUT IN INPUT
    gamma = 0.05_real_8 * au_fs
    ! CJM DBG HARDWIRED UNTIL PUT IN INPUT
    scale = scale * EXP ( -gamma * 0.25_real_8 * dt_ions )
    ! Scale 
    ! HTVEL ( 1, 1 )  = HTVEL ( 1, 1 )  * SCALE
    veps = veps  * scale

    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE vepsupdamp
  ! ==================================================================
  ! ==  PROPOGATES THE CELL VELOCITIES FOR NPH-SHOCK DYNAMCIS       ==
  ! ==--------------------------------------------------------------==



END MODULE vepsup_utils
