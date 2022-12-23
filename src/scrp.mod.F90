MODULE scrp
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==================================================================
  REAL(real_8), PARAMETER :: zslow=-99999999._real_8 
  REAL(real_8), PARAMETER :: zshig= 99999999._real_8 
  INTEGER, PARAMETER :: mxscr= 10 
  INTEGER, PARAMETER :: mxptr= 25 
  ! ==--------------------------------------------------------------==
  REAL(real_8) :: nllsp(mxscr)
  INTEGER :: iiscr(mxscr),nxscr
  ! ==--------------------------------------------------------------==
  INTEGER :: ips_len(mxscr),ips_num(mxscr)
  INTEGER :: ips_low(mxscr),ips_hig(mxscr)
  ! ==--------------------------------------------------------------==
  INTEGER :: ips_ptr1(mxptr,mxscr),ips_ptr2(mxptr,mxscr)
  ! ==--------------------------------------------------------------==
  CHARACTER (len=10) :: ips_name(mxptr,mxscr)
  ! ==================================================================
END MODULE scrp
