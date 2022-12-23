MODULE rnlsm_utils
  USE kinds,                           ONLY: real_8
  USE nlps,                            ONLY: nlm
  USE rnlsm1_utils,                    ONLY: give_scr_rnlsm1,&
                                             rnlsm1
  USE rnlsm2_utils,                    ONLY: give_scr_rnlsm2,&
                                             rnlsm2
  USE rnlsmd_utils,                    ONLY: rnlsmd
  USE system,                          ONLY: cntl
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rnlsm
  PUBLIC :: give_scr_rnlsm

CONTAINS

  ! ==================================================================
  SUBROUTINE rnlsm(c0,nstate,ikpt,ikind,tfor)
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8)                          :: c0(:,:)
    INTEGER                                  :: nstate, ikpt, ikind
    LOGICAL                                  :: tfor

    CHARACTER(*), PARAMETER                  :: procedureN = 'rnlsm'

    INTEGER                                  :: isub

! ==--------------------------------------------------------------==
! If no non-local components -> return.

    IF (nlm.EQ.0) RETURN
    CALL tiset(procedureN,isub)
    ! ==--------------------------------------------------------------==
    IF (cntl%tfdist) THEN
       CALL rnlsmd(c0,nstate,ikind)
    ELSE
       CALL rnlsm1(c0,nstate,ikind)
    ENDIF
    IF (tfor) CALL rnlsm2(c0,nstate,ikpt,ikind)
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)

  END SUBROUTINE rnlsm
  ! ==================================================================
  SUBROUTINE give_scr_rnlsm(lrnlsm,tag,nstate,tfor)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lrnlsm
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate
    LOGICAL                                  :: tfor

    INTEGER                                  :: lrnlsm1, lrnlsm2

    CALL give_scr_rnlsm1(lrnlsm1,tag,nstate)
    IF (tfor) THEN
       CALL give_scr_rnlsm2(lrnlsm2,tag,nstate)
    ELSE
       lrnlsm2=0
    ENDIF
    lrnlsm=MAX(lrnlsm1,lrnlsm2)
    tag   ='MAX(LRNLSM1,LRNLSM2)'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_rnlsm
  ! ==================================================================


END MODULE rnlsm_utils
