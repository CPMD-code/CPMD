MODULE dqgalloc_utils
  USE error_handling,                  ONLY: stopgm
  USE fft_maxfft,                      ONLY: maxfft
  USE pslo,                            ONLY: pslo_com
  USE spin,                            ONLY: clsd,&
                                             lspin2
  USE str2,                            ONLY: dqg
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             ncpw

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: dqgalloc

CONTAINS

  ! ==================================================================
  SUBROUTINE dqgalloc
    ! ==--------------------------------------------------------------==
    CHARACTER(*), PARAMETER                  :: procedureN = 'dqgalloc'

    INTEGER                                  :: ierr, ldqg

! ==--------------------------------------------------------------==

    ldqg=1
    IF (lspin2%tlse) THEN
       IF (lspin2%troks) THEN
          ldqg=MAX(ldqg,2*fpar%nnr1)
       ELSEIF (lspin2%troot) THEN
          ldqg=MAX(ldqg,2*fpar%nnr1)
       ELSEIF (lspin2%tross) THEN
          ldqg=MAX(ldqg,2*fpar%nnr1)
       ELSEIF (lspin2%tcas22) THEN
          ldqg=MAX(ldqg,3*fpar%nnr1)
          ldqg=MAX(ldqg,2*maxfft)
       ENDIF
    ENDIF
    ! VOFRHOB
    IF (cntl%tgc) THEN
       IF (cntl%tpres.OR.cntl%tdiag.OR.cntl%tdiagopt) ldqg=MAX(ldqg,2*clsd%nlsd*fpar%nnr1)
    ENDIF
    IF (cntl%tpres.OR.cntl%tprcp) THEN
       ldqg=MAX(ldqg,2*ncpw%nhg)
       IF (pslo_com%tivan) ldqg=MAX(ldqg,2*ncpw%nhg*6)
       ldqg=MAX(ldqg,4*maxfft)
       IF (cntl%tgc) ldqg=MAX(ldqg,6*maxfft)
    ENDIF
    ALLOCATE(dqg(ldqg,1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE dqgalloc
  ! ==================================================================

END MODULE dqgalloc_utils
