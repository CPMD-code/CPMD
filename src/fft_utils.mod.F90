MODULE fft_utils

  USE error_handling,                  ONLY: stopgm
  USE fft,                             ONLY: fftpoolsize,&
                                             lrxpl,&
                                             lrxpool,&
                                             mxy,&
                                             sp5,&
                                             sp8,&
                                             sp9,&
                                             spm

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: alloc_fft

CONTAINS

  SUBROUTINE alloc_fft( max_nproc )
    INTEGER, INTENT(IN)                      :: max_nproc

    CHARACTER(*), PARAMETER                  :: procedureN = 'alloc_fft'

    INTEGER                                  :: ierr

    ALLOCATE( lrxpl(0:max_nproc-1,2), &
         &    sp5(0:max_nproc-1), &
         &    sp8(0:max_nproc-1), &
         &    sp9(0:max_nproc-1), &
         &    mxy(max_nproc-1), &
         &    lrxpool(0:max_nproc-1,2,fftpoolsize), &
         &    spm(9,0:max_nproc-1,fftpoolsize), &
         &    STAT=ierr )
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    lrxpl   = HUGE(0)
    sp5     = HUGE(0)
    sp8     = HUGE(0)
    sp9     = HUGE(0)
    mxy     = HUGE(0)
    spm     = HUGE(0)
    lrxpool = HUGE(0)

  END SUBROUTINE alloc_fft

END MODULE fft_utils
