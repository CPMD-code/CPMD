MODULE printfor_utils
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_app,&
                                             fo_vmark
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: paral
  USE ropt,                            ONLY: iteropt,&
                                             ropt_mod
  USE store_types,                     ONLY: trajsmall,&
                                             trajsmalln
  USE system,                          ONLY: cntl

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: printfor

CONTAINS

  ! ==================================================================
  SUBROUTINE printfor(taup,fion)
    ! ==--------------------------------------------------------------==

    REAL(real_8)                             :: taup(:,:,:), fion(:,:,:)

    CHARACTER(len=10), PARAMETER             :: f1 = 'FORCES    '

    INTEGER                                  :: i, j, k, l
    INTEGER, SAVE                            :: if1 = fo_vmark
    LOGICAL                                  :: ferror

! ==--------------------------------------------------------------==
! ==  OPEN THE FORCES TRAJECTORY FILE                             ==
! ==--------------------------------------------------------------==

    IF (ropt_mod%rprint) THEN
       IF (paral%io_parent)&
            CALL fileopen(61,f1,fo_app+if1,ferror)
       if1=0
       ! stores ionic configuration and forces ........
       IF (cntl%tsampl) THEN
          l=0
          DO k=1,ions1%nsp
             DO j=1,ions0%na(k)
                l=l+1
                IF (.NOT.trajsmall .OR. l.LE.trajsmalln) THEN
                   IF (paral%io_parent)&
                        WRITE(61,'(I7,6(2X,F22.14))')&
                        iteropt%nfi,(taup(i,j,k),i=1,3),(fion(i,j,k),i=1,3)
                ENDIF
             ENDDO
          ENDDO
       ELSE
          l=0
          DO k=1,ions1%nsp
             DO j=1,ions0%na(k)
                l=l+1
                IF (.NOT.trajsmall .OR. l.LE.trajsmalln) THEN
                   IF (paral%io_parent)&
                        WRITE(61,'(I7,6(2X,F22.14))')&
                        iteropt%nfi,(taup(i,j,k),i=1,3),(fion(i,j,k),i=1,3)
                ENDIF
             ENDDO
          ENDDO
       ENDIF
       IF (paral%io_parent)&
            CALL fileclose(61)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE printfor
  ! ==================================================================

END MODULE printfor_utils
