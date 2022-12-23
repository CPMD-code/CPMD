MODULE updwf_utils
  USE adapttol_utils,                  ONLY: tol_chk_cnvgrad
  USE forcedr_driver,                  ONLY: forcedr
  USE forcedr_utils,                   ONLY: give_scr_forcedr
  USE geq0mod,                         ONLY: geq0
  USE hesele_utils,                    ONLY: hesele
  USE hubbardu,                        ONLY: hubbu
  USE kinds,                           ONLY: real_8
  USE mm_dim_utils,                    ONLY: mm_dim
  USE mm_dimmod,                       ONLY: mm_go_mm,&
                                             mm_go_qm,&
                                             mm_revert
  USE mm_input,                        ONLY: clc
  USE mm_qmmm_forcedr_utils,           ONLY: mm_qmmm_forcedr
  USE norm,                            ONLY: gemax
  USE ortho_utils,                     ONLY: give_scr_ortho,&
                                             ortho,&
                                             preortho
  USE parac,                           ONLY: paral
  USE pcgrad_driver,                   ONLY: pcgrad
  USE pcgrad_utils,                    ONLY: give_scr_pcgrad
  USE pslo,                            ONLY: pslo_com
  USE rnlsm_utils,                     ONLY: rnlsm
  USE ropt,                            ONLY: iteropt,&
                                             ropt_mod
  USE soft,                            ONLY: soft_com
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             cntr,&
                                             ncpw,&
                                             nkpt
  USE testex_utils,                    ONLY: testex
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE tpar,                            ONLY: dt2_elec,&
                                             dt2bye
  USE utils,                           ONLY: zclean
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: updwf
  PUBLIC :: give_scr_updwf

CONTAINS

  ! ==================================================================
  SUBROUTINE updwf(c0,c2,sc0,tau0,fion,pme,gde,vpp,eigv,&
       rhoe,psi,nstate,tfor,update_pot)
    ! ==--------------------------------------------------------------==
    ! ==               UPDATES THE WAVEFUNCTIONS                      ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(:,:), c2(:,:), sc0(:,:)
    REAL(real_8)                             :: tau0(:,:,:), fion(:,:,:)
    COMPLEX(real_8)                          :: pme(*), gde(*)
    REAL(real_8)                             :: vpp(nkpt%ngwk), rhoe(:,:)
    COMPLEX(real_8)                          :: psi(:,:)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: eigv(nstate)
    LOGICAL                                  :: tfor, update_pot

    CHARACTER(*), PARAMETER                  :: procedureN = 'updwf'

    INTEGER                                  :: i, ig, isub
    LOGICAL                                  :: status, statusdummy
    REAL(real_8)                             :: pf1

    CALL tiset(procedureN,isub)
    ! ==--------------------------------------------------------------==
    ! ==  UPDATE THE WAVEFUNCTIONS                                    ==
    ! ==--------------------------------------------------------------==
    IF ( cntl%tqmmm ) THEN
       CALL mm_qmmm_forcedr(c0,c2,sc0,rhoe,psi,tau0,fion,&
            eigv,nstate,0,.TRUE.,update_pot,.TRUE.)
       CALL mm_dim(mm_go_mm,status)
       IF (.NOT.tfor) THEN
          CALL zeroing(fion)!,3*maxsys%nax*maxsys%nsx)
       ENDIF
       CALL mm_dim(mm_go_qm,statusdummy)
    ELSE
       CALL forcedr(c0,c2,sc0,rhoe,psi,tau0,fion,eigv,&
            nstate,1,.TRUE.,tfor)
       IF (.NOT.tfor) THEN
          CALL zeroing(fion)!,3*maxsys%nax*maxsys%nsx)
       ENDIF
    ENDIF

    !
    IF(cntl%thubb) hubbu%tpom=.false. 
    IF (paral%qmnode.AND..NOT.clc%classical)THEN
       ! EHR[
       IF (tfor.AND.cntl%tmdeh) THEN
          ! DO NOTHING
       ELSE
          ! EHR]
          IF (cntl%diis) THEN
             IF (ropt_mod%sdiis) THEN
                CALL hesele(dt2bye,vpp)
             ENDIF
             CALL odiis(c0,c2,vpp,nstate,pme,gde,dt2bye,ropt_mod%sdiis)
             IF (ropt_mod%sdiis) THEN
                iteropt%ndisrs=iteropt%ndisrs+1
                iteropt%ndistp=0
                ropt_mod%sdiis=.FALSE.
             ELSE
                iteropt%ndistp=iteropt%ndistp+1
             ENDIF
          ELSE IF (cntl%pcg) THEN
             IF (ropt_mod%spcg) THEN
                CALL hesele(dt2bye,vpp)
             ENDIF
             CALL pcgrad(c0,c2,sc0,vpp,pme,&
                  rhoe,psi,tau0,nstate,ropt_mod%spcg)
             ropt_mod%spcg=.FALSE.
          ELSE IF (cntl%tsde) THEN
             IF (cntl%prec) THEN
                CALL hesele(dt2bye,vpp)
                pf1=dt2_elec*cntr%hthrs/cntr%emass
                DO i=1,nstate
                   DO ig=1,ncpw%ngw
                      c0(ig,i)=c0(ig,i)+vpp(ig)*pf1*c2(ig,i)
                   ENDDO
                ENDDO
             ELSE
                CALL daxpy(2*nstate*ncpw%ngw,dt2bye,c2,1,c0,1)
             ENDIF
             IF (cnti%iproj.LE.1) THEN
                gemax=2.0_real_8*cntr%tolog
             ENDIF
          ENDIF
          IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)
          CALL tol_chk_cnvgrad(ropt_mod%convwf,soft_com%exsoft,gemax)
          ! EHR[
       ENDIF
       ! EHR]
       ! ==--------------------------------------------------------------==
       ! ==  ORTHOGONALIZATION                                           ==
       ! ==--------------------------------------------------------------==
       IF (cntl%nonort) THEN
          IF (geq0) THEN
             CALL zclean(c0,nstate,ncpw%ngw)
          ENDIF
       ELSE
          CALL preortho(c0,nstate)
          IF (pslo_com%tivan) THEN
             CALL rnlsm(c0,nstate,1,1,.FALSE.)
          ENDIF
          CALL ortho(nstate,c0,c2)
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    IF ( cntl%tqmmm ) THEN
       IF (clc%classical) ropt_mod%convwf=.TRUE.! no iterations
       CALL mm_dim(mm_revert,status)
    ENDIF
    IF(cntl%thubb) hubbu%tpom=.true.
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE updwf
  ! ==================================================================
  SUBROUTINE give_scr_updwf(lupdwf,tag,nstate,tfor)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lupdwf
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate
    LOGICAL                                  :: tfor

    INTEGER                                  :: lforces, lortho, lpcgrad

    lortho=0
    CALL give_scr_forcedr(lforces,tag,nstate,.TRUE.,tfor)
    IF (.NOT.cntl%nonort) THEN
       CALL give_scr_ortho(lortho,tag,nstate)
    ENDIF
    ! 
    IF (cntl%pcg.OR.cntl%tpcgfic) THEN
       CALL give_scr_pcgrad(lpcgrad,tag,nstate)
    ELSE
       lpcgrad=0
    ENDIF
    ! 
    lupdwf=MAX(lforces,lortho,lpcgrad)+10
    tag='MAX(LFORCES,LORTHO,LPCGRAD)'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_updwf
  ! ==================================================================

END MODULE updwf_utils
