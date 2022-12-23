MODULE vofrho_utils
  USE error_handling,                  ONLY: stopgm
  USE forcep_utils,                    ONLY: rhoe_psi_size
  USE isos,                            ONLY: isos1,&
                                             isos3
  USE kinds,                           ONLY: real_8
  USE poin,                            ONLY: potr
  USE spin,                            ONLY: clsd,&
                                             lspin2
  USE store_types,                     ONLY: store1
  USE system,                          ONLY: fpar,&
                                             ncpw
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE vofrhoa_utils,                   ONLY: give_scr_vofrhoa,&
                                             vofrhoa
  USE vofrhob_utils,                   ONLY: give_scr_vofrhob,&
                                             vofrhob
  USE vofrhoc_utils,                   ONLY: give_scr_vofrhoc,&
                                             vofrhoc
  USE vofrhoh_utils,                   ONLY: give_scr_vofrhoh,&
                                             vofrhoh
  USE vofrhos_utils,                   ONLY: give_scr_vofrhor,&
                                             give_scr_vofrhos,&
                                             vofrhor,&
                                             vofrhos
  USE vofrhot_utils,                   ONLY: give_scr_vofrhot,&
                                             vofrhot

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: vofrho
  PUBLIC :: give_scr_vofrho

CONTAINS

  ! ==================================================================
  SUBROUTINE vofrho(tau0,fion,rhoe,v,tfor,tstress)
    ! ==--------------------------------------------------------------==
    ! == FION: forces acting on ions (IF TFOR)                        ==
    ! == RHOE: in  electronic density in real space                   ==
    ! ==       out potential in real space                            ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:), fion(:,:,:), &
                                                rhoe(:,:)
    COMPLEX(real_8)                          :: v(:,:)
    LOGICAL                                  :: tfor, tstress

    CHARACTER(*), PARAMETER                  :: procedureN = 'vofrho'

    COMPLEX(real_8), ALLOCATABLE             :: vtemp(:,:)
    INTEGER                                  :: ierr, il_rhoe_1d, il_rhoe_2d, &
                                                isub

    CALL tiset(procedureN,isub)

    ALLOCATE(vtemp(ncpw%nhg, clsd%nlsd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! Allocation of POTR if SPOT and POTR not allocated.
    IF (store1%spot) THEN
       IF (.NOT.ALLOCATED(potr)) THEN
          CALL rhoe_psi_size(il_rhoe_1d=il_rhoe_1d,il_rhoe_2d=il_rhoe_2d)
          ALLOCATE(potr(il_rhoe_1d,il_rhoe_2d),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
    ENDIF
    IF (isos1%tclust) THEN
       IF (isos3%ps_type.EQ.1) THEN
          CALL vofrhoh(tau0,fion,rhoe(:,1),v(:,1),vtemp(:,1),tfor)
       ELSEIF (isos3%ps_type.EQ.2) THEN
          CALL vofrhot(tau0,fion,rhoe,v(:,1),vtemp,tfor)
       ELSEIF (isos3%ps_type.EQ.3) THEN
          CALL vofrhot(tau0,fion,rhoe,v(:,1),vtemp,tfor)
       ELSE
          CALL stopgm(procedureN,'PS_TYPE not implemented',& 
               __LINE__,__FILE__)
       ENDIF
    ELSE
       CALL vofrhoa(tau0,fion,rhoe(:,1),v(:,1),vtemp(:,1),tfor,tstress)
    ENDIF
    IF (lspin2%tlse)THEN
       IF (lspin2%troks.OR.lspin2%troot) THEN
          CALL vofrhos(fion,rhoe,v,vtemp,tfor,tstress)
       ELSEIF (lspin2%tcas22) THEN
          CALL vofrhos(fion,rhoe,v,vtemp,tfor,tstress)
          CALL vofrhoc(tau0,rhoe(:,4:),v,vtemp,tfor,tstress)
       ELSEIF (lspin2%tross) THEN
          CALL vofrhor(fion,rhoe,v,vtemp,tfor,tstress)
       ENDIF
    ELSE
       CALL vofrhob(fion,rhoe,v,vtemp,tfor,tstress)
    ENDIF
    DEALLOCATE(vtemp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    IF (store1%spot) CALL dcopy(fpar%nnr1*clsd%nlsd,rhoe,1,potr,1)
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    RETURN
  END SUBROUTINE vofrho
  ! ==================================================================
  SUBROUTINE give_scr_vofrho(lvofrho,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lvofrho
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: lc, lvofrhoa, lvofrhob, &
                                                lvofrhos

    IF (isos1%tclust) THEN
       IF (isos3%ps_type.EQ.1) THEN
          CALL give_scr_vofrhoh(lvofrhoa,tag)
       ELSEIF (isos3%ps_type.EQ.2) THEN
          CALL give_scr_vofrhot(lvofrhoa,tag)
       ELSEIF (isos3%ps_type.EQ.3) THEN
          CALL give_scr_vofrhot(lvofrhoa,tag)
       ELSE
          CALL stopgm('VOFRHO','PS_TYPE not implemented',& 
               __LINE__,__FILE__)
       ENDIF
    ELSE
       CALL give_scr_vofrhoa(lvofrhoa,tag)
    ENDIF
    IF (lspin2%tlse) THEN
       IF (lspin2%troks.OR.lspin2%troot) THEN
          CALL give_scr_vofrhos(lvofrhos,tag)
       ELSEIF (lspin2%tcas22) THEN
          CALL give_scr_vofrhos(lvofrhos,tag)
          CALL give_scr_vofrhoc(lc,tag)
          lvofrhos=MAX(lc,lvofrhos)
       ELSEIF (lspin2%tross) THEN
          CALL give_scr_vofrhor(lvofrhos,tag)
       ENDIF
       lvofrhob=0
    ELSE
       CALL give_scr_vofrhob(lvofrhob,tag)
       lvofrhos=0
    ENDIF
    lvofrho=2*ncpw%nhg*clsd%nlsd+MAX(lvofrhoa,lvofrhob,lvofrhos)+100
    tag   ='2*NHG*NLSD+MAX(LVOFRHO(A|B|S))'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_vofrho
  ! ==================================================================


END MODULE vofrho_utils
