MODULE updwf_p_utils
  USE csize_utils,                     ONLY: csize
  USE error_handling,                  ONLY: stopgm
  USE forces_p_utils,                  ONLY: forces_p,&
                                             give_scr_forces_p
  USE geq0mod,                         ONLY: geq0
  USE kinds,                           ONLY: real_8
  USE norm,                            ONLY: cnorm,&
                                             gemax
  USE odiis_p_utils,                   ONLY: odiis_p
  USE ortho_utils,                     ONLY: give_scr_ortho
  USE pcgrad_p_utils,                  ONLY: pcgrad_p
  USE perturbation_p_utils,            ONLY: ortho_p,&
                                             restrain_ngw_zero
  USE response_pmod,                   ONLY: dmbi,&
                                             response1,&
                                             response2,&
                                             s_star
  USE ropt,                            ONLY: ropt_mod
  USE simple_model_p_utils,            ONLY: simple_ortho_p
  USE soft,                            ONLY: soft_com
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             ncpw
  USE testex_utils,                    ONLY: testex
  USE tpar,                            ONLY: dt2bye
  USE utils,                           ONLY: zclean

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: updwf_p
  PUBLIC :: give_scr_updwf_p

CONTAINS

  ! ==================================================================
  SUBROUTINE updwf_p(c0,c1,c12,h1nl,z11,&
       eivps1,eirop1,&
       tau0,fion,pme,gde,vpp,&
       rhoe,drhoe,psi,nstate,firstcall,reset_cg)
    ! ==--------------------------------------------------------------==
    ! ==               updates the wavefunctions                      ==
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8)                          :: eivps1(ncpw%nhg), &
                                                eirop1(ncpw%nhg)
    REAL(real_8)                             :: tau0(:,:,:), fion(:,:,:)
    REAL(real_8), DIMENSION(:, :), &
      INTENT(inout)                          :: pme, gde
    REAL(real_8)                             :: vpp(ncpw%ngw), &
                                                rhoe(fpar%nnr1), &
                                                drhoe(fpar%nnr1)
    COMPLEX(real_8)                          :: psi(fpar%nnr1)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: z11(nstate,nstate)
    COMPLEX(real_8) :: h1nl(ncpw%ngw,nstate), c12(ncpw%ngw,nstate), &
      c1(ncpw%ngw,nstate), c0(ncpw%ngw,nstate)
    LOGICAL                                  :: firstcall, reset_cg

! variables
! ==--------------------------------------------------------------==
! ==  calculate the force on the wavefunctions                    ==
! ==--------------------------------------------------------------==

    IF (dmbi%cutoff_restr) CALL restrain_ngw_zero(c1,dmbi%ngw_zero,ncpw%ngw,nstate)
    CALL forces_p(c0,c1,c12,h1nl,z11,drhoe,eivps1,eirop1,&
         tau0,fion,rhoe,psi,&
         nstate,firstcall)
    ! ==--------------------------------------------------------------==
    ! ==  update the wavefunctions                                    ==
    ! ==--------------------------------------------------------------==
    IF (response1%pcg_p) THEN
       CALL pcgrad_p(c0,c1,c12,vpp,&
            psi,nstate,reset_cg,z11)
    ELSEIF (response1%diis_p) THEN
       IF (firstcall .AND. response1%opt1st) THEN
          CALL pcgrad_p(c0,c1,c12,vpp,&
               psi,nstate,reset_cg,z11)
          ropt_mod%sdiis = .TRUE.
       ELSE
          ! DMB
          IF (dmbi%tsimple_model) THEN
             CALL simple_ortho_p(nstate,c0,c12,S_star)
          ELSE
             CALL ortho_p(nstate,c0,c12)
          ENDIF
          ! DMB
          CALL csize(c12,nstate,gemax,cnorm)
          CALL odiis_p(c1,c12,vpp,z11,nstate,pme,gde,dt2bye,ropt_mod%sdiis)

          ! DMB added C0|C1 orthogonalisation (also for cntl%pcg)
          IF (dmbi%tsimple_model) THEN
             CALL simple_ortho_p(nstate,c0,c1,S_star)
          ELSE
             CALL ortho_p(nstate,c0,c1)
          ENDIF
          ! DMB

          IF (geq0) THEN
             CALL zclean(c1,nstate,ncpw%ngw)
          ENDIF
       ENDIF
    ELSE
       CALL stopgm('UPDWF_P',&
            'EITHER nothing or nothing at all.',& 
            __LINE__,__FILE__)
    ENDIF

    IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)

    IF (gemax.LT.response2%tolog_p .OR. soft_com%exsoft) ropt_mod%convwf=.TRUE.

    firstcall = .FALSE.

    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE updwf_p
  ! ==================================================================



  ! ==================================================================
  SUBROUTINE give_scr_updwf_p(lupdwf,tag,nstate)
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: lupdwf
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate

    INTEGER                                  :: lforces, lortho, lpcgrad

! variables
! ==--------------------------------------------------------------==

    lortho=0
    CALL give_scr_forces_p(lforces,tag,nstate)
    IF (.NOT.cntl%nonort) THEN
       CALL give_scr_ortho(lortho,tag,nstate)
    ENDIF
    ! call give_scr_pcgrad_p(lpcgrad,tag,nstate)
    lpcgrad = 2*ncpw%ngw*nstate
    lupdwf=MAX(lforces,lortho,lpcgrad)
    tag='max(lforces,lortho,lpcgrad)'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_updwf_p
  ! ==================================================================


END MODULE updwf_p_utils
