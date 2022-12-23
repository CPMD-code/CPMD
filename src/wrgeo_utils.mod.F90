MODULE wrgeo_utils
  USE adat,                            ONLY: elem
  USE cnst,                            ONLY: fbohr
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_app,&
                                             fo_verb
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE mm_dim_utils,                    ONLY: mm_dim
  USE mm_dimmod,                       ONLY: mm_go_mm,&
                                             mm_revert,&
                                             mmdim,&
                                             naq,&
                                             nat_grm
  USE mm_input,                        ONLY: lqmmm
  USE parac,                           ONLY: paral
  USE ropt,                            ONLY: infi
  USE store_types,                     ONLY: cprint,&
                                             iprint_force,&
                                             rout1
  USE system,                          ONLY: cnti
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: wrgeo
  PUBLIC :: wrgeof
  PUBLIC :: wrgeof_inr
  PUBLIC :: wrgeox

CONTAINS

  ! ==================================================================
  SUBROUTINE wrgeo(tau0)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:)

    INTEGER                                  :: ia, iat, is, k

    IF (paral%io_parent)&
         WRITE(6,'(/,1X,64("*"))')
    IF (paral%io_parent)&
         WRITE(6,'(" *",21X,A,22X,"*")')&
         ' ATOMIC COORDINATES'
    IF (paral%io_parent)&
         WRITE(6,'(1X,64("*"))')
    IF (.NOT.lqmmm%qmmm .OR. ions1%nat.LT.500)THEN
       iat=0
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             iat=iat+1
             IF (paral%io_parent)&
                  WRITE(6,'(I8,6X,A2,4X,3(F15.6))')&
                  iat,elem%el(ions0%iatyp(is)),(tau0(k,ia,is),k=1,3)
          ENDDO
       ENDDO
    ELSE
       iat=0
       DO is=1,mmdim%nspq
          DO ia=1,NAq(is)
             iat=iat+1
             IF (paral%io_parent)&
                  WRITE(6,'(I8,6X,A2,4X,3(F15.6))')&
                  iat,elem%el(ions0%iatyp(is)),(tau0(k,ia,is),k=1,3)
          ENDDO
       ENDDO
    ENDIF
    IF (paral%io_parent)&
         WRITE(6,'(1X,64("*"),/)')
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE wrgeo
  ! ==================================================================
  SUBROUTINE wrgeof(tau0,fion)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:), fion(:,:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'wrgeof'

    INTEGER                                  :: ia, iat, is, isub, k

    CALL tiset(procedureN,isub)

    IF (cprint%iprint(iprint_force).LE.0) THEN
       CALL wrgeo(tau0)
       CALL tihalt(procedureN,isub)
       RETURN
    ENDIF
    IF (paral%io_parent)&
         WRITE(6,'(/,T4,A,T18,A,T41,A)') 'ATOM',&
         'COORDINATES',&
         'GRADIENTS (-FORCES)'
    IF (.NOT.lqmmm%qmmm .OR. ions1%nat.LT.500)THEN
       iat=0
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             iat=iat+1
             IF (paral%io_parent)&
                  WRITE(6,'(1X,I3,1X,A2,3(F8.4),1X,3(1PE11.3))')&
                  iat,elem%el(ions0%iatyp(is)),&
                  (tau0(k,ia,is),k=1,3),(fion(k,ia,is),k=1,3)
          ENDDO
       ENDDO
    ELSE
       iat=0
       DO is=1,mmdim%nspq
          DO ia=1,NAq(is)
             iat=iat+1
             IF (paral%io_parent)&
                  WRITE(6,'(1X,I3,1X,A2,3(F8.4),1X,3(1PE11.3))')&
                  iat,elem%el(ions0%iatyp(is)),&
                  (tau0(k,ia,is),k=1,3),(fion(k,ia,is),k=1,3)
          ENDDO
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    RETURN
  END SUBROUTINE wrgeof
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE wrgeof_inr(tau0,fion)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:), fion(:,:,:)

    INTEGER                                  :: ia, iat, is, k

    IF (paral%io_parent)&
         WRITE(6,'(/,T4,A,T18,A,T41,A)') 'ATOM',&
         'COORDINATES',&
         'GRADIENTS (-FORCES)'
    iat=0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          iat=iat+1
          IF (paral%io_parent)&
               WRITE(6,'(1X,I3,1X,A2,3(F8.4),1X,3(1PE11.3))')&
               iat,elem%el(ions0%iatyp(is)),&
               (tau0(k,ia,is),k=1,3),(fion(k,ia,is),k=1,3)
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE wrgeof_inr
  ! ==================================================================
  SUBROUTINE wrgeox(tau0)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'wrgeox'

    INTEGER                                  :: ia, iat, ierr, is, iunit, k
    INTEGER, ALLOCATABLE                     :: gr_iat(:)
    INTEGER, SAVE                            :: ifile1 = fo_verb
    LOGICAL                                  :: ferror, status
    REAL(real_8), ALLOCATABLE                :: gr_tau(:,:)

    IF (rout1%xgout .AND. MOD(infi,cnti%ngxyz).EQ.0) THEN
       iunit=1
       CALL mm_dim(mm_go_mm,status)
       ALLOCATE(gr_tau(3,ions1%nat),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(gr_iat(ions1%nat+1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)

       IF (paral%io_parent)&
            CALL fileopen(iunit,'GEO_OPT.xyz',fo_app+ifile1,ferror)
       ifile1=0

       iat=0
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             iat=iat+1
             gr_iat(nat_grm(iat))=ions0%iatyp(is)
             DO k=1,3
                gr_tau(k,nat_grm(iat))=tau0(k,ia,is)
             ENDDO
          ENDDO
       ENDDO
       IF (paral%io_parent)&
            WRITE(iunit,'(I8)') ions1%nat
       IF (paral%io_parent)&
            WRITE(iunit,'(I8)') infi
       DO iat=1,ions1%nat
          IF (paral%io_parent)&
               WRITE(iunit,'(A3,3F20.12)')&
               elem%el(gr_iat(iat)),(gr_tau(k,iat)/fbohr,k=1,3)
       ENDDO
       IF (paral%io_parent)&
            CALL fileclose(iunit)
       CALL mm_dim(mm_revert,status)
       DEALLOCATE(gr_tau,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(gr_iat,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE wrgeox
  ! ==================================================================

END MODULE wrgeo_utils
