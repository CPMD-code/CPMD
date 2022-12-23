MODULE lr_diag_utils
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE linres,                          ONLY: td01,&
                                             td03
  USE system,                          ONLY: ncpw
  USE td_dav_utils,                    ONLY: give_scr_td_dav,&
                                             td_dav,&
                                             td_sub
  USE td_nhdav_utils,                  ONLY: give_scr_td_nhdav,&
                                             td_nhdav
  USE td_pcg_utils,                    ONLY: give_scr_td_pcg,&
                                             td_pcg
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: lr_diag
  PUBLIC :: give_scr_lr_diag

CONTAINS

  ! ==================================================================
  SUBROUTINE lr_diag(c0,c1,c2,sc0,eigv,ddxc,drhoe,psi,&
       teig,nstate,nroot,orbital,kprint)
    ! ==--------------------------------------------------------------==
    ! ==    Calls diagonalizers for TD-DFT                            ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(:,:), c2(:,:), sc0(*)
    REAL(real_8)                             :: eigv(:), ddxc(:,:), drhoe(:,:)
    COMPLEX(real_8)                          :: psi(:,:)
    REAL(real_8)                             :: teig(:)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c1(ncpw%ngw,nstate,*)
    INTEGER                                  :: nroot
    CHARACTER(len=*)                         :: orbital
    INTEGER                                  :: kprint

    CHARACTER(*), PARAMETER                  :: procedureN = 'lr_diag'

    INTEGER                                  :: isub

    CALL tiset(procedureN,isub)
    ! ==--------------------------------------------------------------==
    IF (td01%ldiag.EQ.1) THEN
       IF (td03%tda.AND.td01%msubta.GT.0) THEN
          CALL td_sub(c0,c1,c2,sc0,ddxc,psi,eigv,drhoe,&
               teig,td01%msubta,nstate,nroot,orbital,kprint)
       ELSE
          CALL td_dav(c0,c1,c2,sc0,ddxc,psi,eigv,drhoe,&
               teig,nstate,nroot,orbital,kprint)
       ENDIF
    ELSEIF (td01%ldiag.EQ.2) THEN
       CALL td_nhdav(c0,c1,c2,sc0,ddxc,psi,eigv,drhoe,&
            teig,nstate,nroot,orbital,kprint)
    ELSEIF (td01%ldiag.EQ.3) THEN
       CALL td_pcg(c0,c1,c2,sc0,ddxc,psi,eigv,drhoe,&
            teig,td01%msubta,nstate,nroot,orbital,kprint)
    ELSE
       CALL stopgm('LR_DIAG','METHOD NOT IMPLEMENTED',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE lr_diag
  ! ==================================================================
  SUBROUTINE give_scr_lr_diag(llr_diag,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: llr_diag
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: ltddiag

    ltddiag=0
    IF (td01%ldiag.EQ.1) THEN
       CALL give_scr_td_dav(ltddiag,tag)
    ELSEIF (td01%ldiag.EQ.2) THEN
       CALL give_scr_td_nhdav(ltddiag,tag)
    ELSEIF (td01%ldiag.EQ.3) THEN
       CALL give_scr_td_pcg(ltddiag,tag)
    ENDIF
    llr_diag=MAX(ltddiag,100)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_lr_diag
  ! ==================================================================

END MODULE lr_diag_utils
