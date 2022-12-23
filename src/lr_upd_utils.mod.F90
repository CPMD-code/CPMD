MODULE lr_upd_utils
  USE csize_utils,                     ONLY: csize
  USE dotp_utils,                      ONLY: dotp
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE geq0mod,                         ONLY: geq0
  USE kinds,                           ONLY: real_8
  USE linres,                          ONLY: lr01,&
                                             lr02
  USE lr_force_utils,                  ONLY: give_scr_lr_force,&
                                             lr_force
  USE lr_ortho_utils,                  ONLY: give_scr_lr_ortho,&
                                             lr_ortho
  USE lr_pcg_utils,                    ONLY: lr_pcg
  USE mp_interface,                    ONLY: mp_sum
  USE norm,                            ONLY: cnorm,&
                                             gemax
  USE parac,                           ONLY: parai
  USE ropt,                            ONLY: ropt_mod
  USE soft,                            ONLY: soft_com
  USE system,                          ONLY: ncpw,&
                                             nkpt
  USE testex_utils,                    ONLY: testex
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: zclean
  USE zdiis_utils,                     ONLY: zdiis

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: lr_upd
  PUBLIC :: give_scr_lr_upd

CONTAINS

  ! ==================================================================
  SUBROUTINE lr_upd(c0,c1,c2,sc0,e2,ddxc,psi,eigv,drhoe,h1nl,vpp,&
       pme,gde,nstate,TYPE,orbital,knfi,ldofor)
    ! ==--------------------------------------------------------------==
    ! ==               updates the wavefunctions                      ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(:,:), c1(:,:), c2(:,:), &
                                                sc0(*)
    REAL(real_8)                             :: e2(:), ddxc(:,:)
    COMPLEX(real_8)                          :: psi(:,:)
    REAL(real_8)                             :: eigv(:), drhoe(:,:)
    COMPLEX(real_8)                          :: h1nl(nkpt%ngwk,*)
    REAL(real_8)                             :: vpp(:)
    COMPLEX(real_8)                          :: pme(*), gde(*)
    INTEGER                                  :: nstate
    CHARACTER(len=*)                         :: TYPE, orbital
    INTEGER                                  :: knfi
    LOGICAL                                  :: ldofor

    CHARACTER(*), PARAMETER                  :: procedureN = 'lr_upd'

    COMPLEX(real_8), ALLOCATABLE, SAVE       :: cback(:,:)
    INTEGER                                  :: i, ierr, ig, is, ism, isub
    INTEGER, SAVE                            :: isfirst = 0, ishist = 0
    LOGICAL                                  :: lresta, lsave
    LOGICAL, SAVE                            :: fdiis
    REAL(real_8)                             :: e4
    REAL(real_8), SAVE                       :: cnhist(10), gehist(10)

    CALL tiset(procedureN,isub)
    ! ==--------------------------------------------------------------==
    ! ==  calculate the force                                         ==
    ! ==--------------------------------------------------------------==
    IF (isfirst.EQ.0) THEN
       ALLOCATE(cback(nkpt%ngwk,nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       isfirst=1
    ENDIF
    ! ==--------------------------------------------------------------==
    e4=e2(4)
999 CONTINUE
    e2(4)=e4
    DO i=1,nstate
       e2(4)=e2(4)-2.0_real_8*dotp(ncpw%ngw,h1nl(:,i),c1(:,i))
    ENDDO
    IF (ldofor) THEN
       CALL lr_force(c0,c1,c2,sc0,e2,ddxc,psi,drhoe,eigv,&
            nstate,TYPE,orbital)
       CALL mp_sum(e2,5,parai%allgrp)
       CALL daxpy(2*ncpw%ngw*nstate,1._real_8,h1nl,1,c2,1)
       IF (geq0) CALL zclean(c2,nstate,ncpw%ngw)
       CALL lr_ortho(nstate,c0,c2)
       CALL csize(c2,nstate,gemax,cnorm)
       IF (gemax.LT.lr02%tol_lr) ropt_mod%convwf=.TRUE.
    ENDIF
    ! ==--------------------------------------------------------------==
    ishist=ishist+1
    is=MOD(ishist-1,10)+1
    gehist(is)=gemax
    cnhist(is)=cnorm
    ism=MOD(is+10-2,10)+1
    lresta=.FALSE.
    IF (ishist.GT.1) THEN
       IF (10._real_8*gehist(ism).LT.gemax .AND.10._real_8*cnhist(ism).LT.cnorm )&
            lresta=.TRUE.
    ENDIF
    ism=MOD(is+10-4,10)+1
    IF (ishist.GT.4) THEN
       IF (gehist(ism).LT.gemax .AND.cnhist(ism).LT.cnorm)&
            lresta=.TRUE.
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ==  update the response functions                               ==
    ! ==--------------------------------------------------------------==
    IF (lresta) THEN
       ishist=0
       CALL dcopy(2*nkpt%ngwk*nstate,cback,1,c1,1)
       ropt_mod%spcg=.TRUE.
       ropt_mod%sdiis=.TRUE.
       fdiis=.TRUE.
       GOTO 999
    ELSEIF (.NOT.ropt_mod%convwf) THEN
       IF (lr01%lopti.EQ.0) THEN
          IF (gemax.GT.lr02%thauto(1).OR.knfi.GT.1) THEN
             CALL lr_pcg(c1,c2,vpp,pme,c0,drhoe,ddxc,eigv,sc0,psi,&
                  nstate,ropt_mod%spcg)
             ropt_mod%sdiis=.TRUE.
             fdiis=.TRUE.
          ELSEIF (gemax.GT.lr02%thauto(2)) THEN
             CALL odiis(c1,c2,vpp,nstate,pme,gde,lr02%dlrs,ropt_mod%sdiis)
             ropt_mod%spcg=.TRUE.
             fdiis=.TRUE.
          ELSE
             CALL zdiis(c1,c2,vpp,eigv,nstate,pme,gde,lr02%dlrs,fdiis,&
                  orbital)
             ropt_mod%spcg=.TRUE.
             ropt_mod%sdiis=.TRUE.
          ENDIF
       ELSEIF (lr01%lopti.EQ.1) THEN
          DO is=1,nstate
             DO ig=1,ncpw%ngw
                c1(ig,is)=c1(ig,is)+lr02%dlrs*vpp(ig)*c2(ig,is)
             ENDDO
          ENDDO
       ELSEIF (lr01%lopti.EQ.2) THEN
          CALL lr_pcg(c1,c2,vpp,pme,c0,drhoe,ddxc,eigv,sc0,psi,&
               nstate,ropt_mod%spcg)
       ELSEIF (lr01%lopti.EQ.3) THEN
          CALL odiis(c1,c2,vpp,nstate,pme,gde,lr02%dlrs,ropt_mod%sdiis)
       ENDIF
       CALL lr_ortho(nstate,c0,c1)
    ELSE
       ishist=0
       fdiis=.TRUE.
    ENDIF
    IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)
    ropt_mod%convwf= ropt_mod%convwf .OR. soft_com%exsoft
    ! ==--------------------------------------------------------------==
    ism=MIN(ishist,10)
    lsave=.TRUE.
    DO is=1,ism
       IF (cnhist(is).LT.cnorm*0.999_real_8) lsave=.FALSE.
    ENDDO
    IF (lsave) CALL dcopy(2*nkpt%ngwk*nstate,c1,1,cback,1)
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE lr_upd
  ! ==================================================================
  SUBROUTINE give_scr_lr_upd(llr_upd,TYPE,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: llr_upd
    CHARACTER(len=*)                         :: TYPE
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: llr_force, llr_ortho

! ==--------------------------------------------------------------==

    CALL give_scr_lr_force(llr_force,TYPE,tag)
    CALL give_scr_lr_ortho(llr_ortho,crge%n)
    llr_upd=MAX(llr_force,llr_ortho)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_lr_upd
  ! ==================================================================

END MODULE lr_upd_utils
