MODULE pstat
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==================================================================
  ! == INCLUDE FILE FOR PARALLEL ROUTINES                           ==
  ! ==--------------------------------------------------------------==
  ! == CMPAR Maximum number of differents message passing calls     ==
  ! ==     IPAR_SEND  SEND/RECEIVE                                  ==
  ! ==     IPAR_CAST  BROADCAST                                     ==
  ! ==     IPAR_GSUM  GLOBAL SUMMATION                              ==
  ! ==     IPAR_GMUL  GLOBAL MULTIPLICATION                         ==
  ! ==     IPAR_AALL  ALL TO ALL COMM                               ==
  ! ==     IPAR_SYNC  SYNCHRONISATION                               ==
  ! ==--------------------------------------------------------------==
  INTEGER, PARAMETER :: cmpar=6 
  INTEGER, PARAMETER :: ipar_send=1,ipar_cast=2,ipar_gsum=3,ipar_gmul=4,&
       ipar_aall=5,ipar_sync=6
  ! ==--------------------------------------------------------------==
  ! == For each different message passing calls:                    ==
  ! == CMLEN    Total Length of messages                            ==
  ! == CMCAL    Number of calls                                     ==
  ! == CMTIM    Total time                                          ==
  ! ==--------------------------------------------------------------==
  REAL(real_8) :: cmlen(cmpar),cmcal(cmpar),cmtim(cmpar)
  ! ==================================================================

END MODULE pstat
