MODULE kdpc
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==================================================================
  ! == INFORMATION ABOUT K POINTS MESH (MONKHORST-PACK)             ==
  ! == USING THE K.P STRATEGY                                       ==
  ! ==                                                              ==
  ! ==   TKDP  : K-POINTS IN THE K.P APPROACH USED (REAL WFN)       ==
  ! ==   NKDP  : NUMBER OF K.P-LIKE K-POINTS (effective)            ==
  ! ==   BMIX  : MIXING FACTOR FOR UPDATING AKDP                    ==
  ! ==================================================================
  INTEGER :: nkdp
  LOGICAL :: tkdp
  REAL(real_8) :: bmix
  ! ==================================================================


END MODULE kdpc
