MODULE hessout_utils
  USE cotr,                            ONLY: cotc0,&
                                             hess
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_def
  USE parac,                           ONLY: paral

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: hessout

CONTAINS

  ! ==================================================================
  SUBROUTINE hessout
    ! ==--------------------------------------------------------------==
    ! ==  STORE NUCLEAR HESSIAN ON DISK                               ==
    ! ==--------------------------------------------------------------==
    ! Variables
    CHARACTER(len=100)                       :: fname
    INTEGER                                  :: i, iunit, j
    LOGICAL                                  :: ferror

! ==--------------------------------------------------------------==

    fname='HESSIAN'
    iunit=21
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent)&
         CALL fileopen(iunit,fname,fo_def,ferror)
    IF (paral%io_parent)&
         REWIND(iunit)
    IF (paral%io_parent)&
         WRITE(iunit,*) cotc0%nodim
    DO i=1,cotc0%nodim
       IF (paral%io_parent)&
            WRITE(iunit,*) (hess(i,j),j=1,cotc0%nodim)
    ENDDO
    IF (paral%io_parent)&
         CALL fileclose(iunit)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE hessout
  ! ==================================================================

END MODULE hessout_utils
