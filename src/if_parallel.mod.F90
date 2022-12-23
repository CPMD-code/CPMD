MODULE if_parallel
  !     ==--------------------------------------------------------------==
  !     == IFMEPOS : RANK in the Interface group                        ==
  !     == IFNPROC : Number of processes in the interface group         ==
  !     == IFGRP   : Interface group communicator                       ==
  !     == IFSOURCE: IFMEPOS of  process responsible for the Y-Z plane  ==
  !     == IFSRCPOS: RANK of IFSOURCE in the QMMMGROUP                  ==
  !     == IFPARENT: IFMEPOS==IFSOURCE                                  ==
  !     ==--------------------------------------------------------------==
  TYPE ifparai_t
     INTEGER :: ifmepos 
     INTEGER :: ifnproc 
     INTEGER :: ifgrp 
     INTEGER :: ifsource 
     INTEGER :: ifsrcpos
  END TYPE ifparai_t
  TYPE (ifparai_t) :: ifparai

  TYPE lifpar_t
     LOGICAL :: ifparent
  END TYPE lifpar_t
  TYPE (lifpar_t) :: lifpar
END MODULE if_parallel
