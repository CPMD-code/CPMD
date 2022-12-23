MODULE afbdr_utils
  USE fnlalloc_utils,                  ONLY: fnl_set,&
                                             fnlalloc,&
                                             fnldealloc
  USE kinds,                           ONLY: real_8
  USE rnlsm_utils,                     ONLY: give_scr_rnlsm,&
                                             rnlsm
  USE system,                          ONLY: ncpw
  USE tdnlfor_utils,                   ONLY: tdnlfor

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: afbdr
  PUBLIC :: give_scr_afbdr

CONTAINS

  ! ==================================================================
  SUBROUTINE afbdr(ca,cb,focc,psi,nstate,fion)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: focc(*)
    COMPLEX(real_8)                          :: psi(:)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: cb(ncpw%ngw,nstate), &
                                                ca(ncpw%ngw,nstate)
    REAL(real_8)                             :: fion(:,:,:)

    CALL rnlsm(ca,nstate,1,1,.TRUE.)
    CALL fnl_set('SAVE')
    CALL fnlalloc(nstate,.TRUE.,.FALSE.)
    CALL rnlsm(cb,nstate,1,1,.TRUE.)
    CALL fnl_set('MIX')
    CALL tdnlfor(fion,focc,nstate)
    CALL fnl_set('SWITCH')
    CALL tdnlfor(fion,focc,nstate)
    CALL fnl_set('MIX')
    CALL fnl_set('SWITCH')
    CALL fnldealloc(.TRUE.,.FALSE.)
    CALL fnl_set('RECV')
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE afbdr
  ! ==================================================================
  SUBROUTINE give_scr_afbdr(lafbdr,nstate,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lafbdr, nstate
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: lrnlsm, lscr

    CALL give_scr_rnlsm(lrnlsm,tag,nstate,.TRUE.)
    lscr=nstate+100
    lafbdr=MAX(lrnlsm,lscr)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_afbdr
  ! ==================================================================

END MODULE afbdr_utils
