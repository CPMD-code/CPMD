MODULE pslo
  USE system,                          ONLY: maxsp

  IMPLICIT NONE

  ! =================================================================
  ! == INCLUDE FILES FOR PSEUDOPOTENTIAL (LOGICAL VARIABLES)       ==
  ! ==-------------------------------------------------------------==
  ! == TIVAN .TRUE. if Vanderbilt PP are used                      ==
  ! == TVAN(1:NSP) For each specie .TRUE. if Vanderbilt PP         ==
  ! == TNUM(1:NSP) Use of numerical values for PP                  ==
  ! == TBIN(1:NSP) For Vanderbilt PP: If PP file is BINARY         ==
  ! == TLOG(1:NSP) For Vanderbilt PP                               ==
  ! == TPSEU(1:NSP)For Vanderbilt PP                               ==
  ! ==-------------------------------------------------------------==

  ! =================================================================
  TYPE :: pslo_com_t
     LOGICAL :: tivan
     LOGICAL :: tvan(maxsp)
     LOGICAL :: tnum(maxsp)
     LOGICAL :: tbin(maxsp)
     LOGICAL :: tlog(maxsp)
     LOGICAL :: tpseu(maxsp)
  END TYPE pslo_com_t
  TYPE(pslo_com_t) :: pslo_com

END MODULE pslo
