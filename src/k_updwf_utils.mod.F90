MODULE k_updwf_utils
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE forcedr_utils,                   ONLY: give_scr_forcedr
  USE k_forces_driver,                 ONLY: k_forces
  USE kinds,                           ONLY: real_8
  USE kpnt,                            ONLY: kdrho,&
                                             kdrhon,&
                                             kgemax
  USE norm,                            ONLY: gemax
  USE ortho_utils,                     ONLY: give_scr_ortho
  USE parac,                           ONLY: paral
  USE pcgrad_utils,                    ONLY: give_scr_pcgrad
  USE ropt,                            ONLY: ropt_mod
  USE soft,                            ONLY: soft_com
  USE system,                          ONLY: cntl,&
                                             cntr,&
                                             nkpt
  USE testex_utils,                    ONLY: testex
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: k_updwf
  PUBLIC :: give_scr_kupdwf

CONTAINS

  ! ==================================================================
  SUBROUTINE k_updwf(c0,c2,sc0,tau0,fion,pme,gde,vpp,eigv,&
       rhoe,psi,nstate,tfor,istep)
    ! ==--------------------------------------------------------------==
    ! ==               UPDATES THE WAVEFUNCTIONS                      ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:), fion(:,:,:)
    COMPLEX(real_8)                          :: pme(*), gde(*)
    REAL(real_8)                             :: vpp(:), rhoe(:,:)
    COMPLEX(real_8)                          :: psi(:,:)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: eigv(nstate,nkpt%nkpts)
    COMPLEX(real_8) :: sc0(nkpt%ngwk,nstate), c2(nkpt%ngwk,nstate), &
      c0(nkpt%ngwk,nstate,nkpt%nkpnt)
    LOGICAL                                  :: tfor
    INTEGER                                  :: istep

    INTEGER                                  :: i
    LOGICAL                                  :: calcrho
    REAL(real_8)                             :: f_save(nstate*nkpt%nkpnt)

    IF (cntl%nonort .OR. cntl%tpath) THEN
       CALL stopgm('K_UPDWF','CNTL%NONORT,TPATH NOT IMPLEMENTED',& 
            __LINE__,__FILE__)
    ENDIF
    IF (cntl%tlsd .AND. nstate .GT. INT(crge%nel)) THEN
       CALL stopgm('K_UPDWF',&
            'KPOINT+VIRTUAL STATES only with DIAGONALIZATION methods',& 
            __LINE__,__FILE__)
    ELSEIF (.NOT. cntl%tlsd .AND. nstate .GT. INT(crge%nel/2 + 1)) THEN
       CALL stopgm('K_UPDWF',&
            'KPOINT+VIRTUAL STATES only with DIAGONALIZATION methods',& 
            __LINE__,__FILE__)
    ENDIF

    ! ==--------------------------------------------------------------==
    ! ==  UPDATE THE WAVEFUNCTIONS                                    ==
    ! ==--------------------------------------------------------------==
    calcrho = .TRUE.
    CALL k_forces(c0,c2,sc0,tau0,pme,gde,vpp,eigv,fion,rhoe,psi,&
         f_save,nstate,nkpt%nkpnt,.TRUE.,tfor,kdrho,kdrhon,&
         calcrho,istep)
    IF (.NOT.tfor) CALL zeroing(fion)!,3*maxsys%nax*maxsys%nsx)

    IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)
    IF (gemax.LT.cntr%tolog .OR. soft_com%exsoft) THEN
       ropt_mod%convwf=.TRUE.
    ELSEIF (kdrho .LT. 0.05_real_8*cntr%tolog .AND.cntl%geopt ) THEN
       IF (paral%io_parent)&
            WRITE(6,'(" ==",60("-"),"==")')
       IF (paral%io_parent)&
            WRITE(6,'(A,1PE12.4,A,1PE12.4)')&
            '  CONVERGENCE OF THE DENSITY, DRHOMAX = ',&
            kdrho,' DRHONORM = ', kdrhon
       IF (paral%io_parent)&
            WRITE(6,'(A)') '  WFNS NOT CONVERGED, MAX'
       DO i = 1,nkpt%nkpnt
          IF (paral%io_parent)&
               WRITE(6,'(A,I5,A,1PE12.4)') 'KPOINT NUM.', i,' KGEMAX = ',&
               kgemax(i)
       ENDDO
       IF (paral%io_parent)&
            WRITE(6,'(" ==",60("-"),"==")')
       ropt_mod%convwf=.TRUE.
    ENDIF

    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE k_updwf
  ! ==================================================================
  SUBROUTINE give_scr_kupdwf(lupdwf,tag,nstate,tfor)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lupdwf
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate
    LOGICAL                                  :: tfor

    INTEGER                                  :: lforces, lortho, lpcgrad

    lortho=0
    CALL give_scr_forcedr(lforces,tag,nstate,.TRUE.,tfor)
    CALL give_scr_ortho(lortho,tag,nstate)
    IF (cntl%pcg) THEN
       CALL give_scr_pcgrad(lpcgrad,tag,nstate)
    ELSE
       lpcgrad=0
    ENDIF
    lupdwf=MAX(lforces,lortho,lpcgrad)+10
    tag='MAX(LFORCES,LORTHO,LPCGRAD)'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_kupdwf
  ! ==================================================================

END MODULE k_updwf_utils
