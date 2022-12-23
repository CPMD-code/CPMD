MODULE mm_parallel
  INTEGER, PARAMETER ::  mmbcston=1678,&
       mmbcstoff=91783,&
       mmbcstrev=-2916
  !     ==--------------------------------------------------------------==
  !     == PARENT  : .TRUE if IP = SOURCE                               ==
  !     ==--------------------------------------------------------------==
  TYPE gparal_t
     LOGICAL :: mmparent
     LOGICAL :: mmnode
  END TYPE gparal_t
  TYPE (gparal_t) :: gparal
  !     ==--------------------------------------------------------------==
  !     == NCPUS   : value of NCPUS environment variable                ==
  !     == NPROC   : TOTAL NUMBER OF PROCESSORS                         ==
  !     == ME      : CPU index                                          ==
  !     == MEPOS   : equal to ME                                        ==
  !     == SOURCE  : MASTER CPU index                                   ==
  !     == IGEQ0   : PROCESSOR INDEX OF G=0 COMPONENT                   ==
  !     == ALLGRP  : GIVEN VALUE FOR BROADCAST COMMUNICATION            ==
  !     == NHRAYS  : number of rays for the processor for the density   ==
  !     == NGRAYS  : number of rays for the processor for the wavefunc. ==
  !     ==--------------------------------------------------------------==
  TYPE gparai_t
     INTEGER :: mmmepos
     INTEGER :: mmncpu
     INTEGER :: mmnproc
     INTEGER :: mmgrp
     INTEGER :: mmsource
  END TYPE gparai_t
  TYPE (gparai_t) :: gparai
  !     ==--------------------------------------------------------------==
  !     == RSNODE   : node computes bonded and Real Space               ==
  !     ==               non-bonded interactions                        ==
  !     == RSPARENT : master node for Real Space (check if it is needed)==
  !     == LSNODE   : node computes lattice sum (Ewald) forces          ==
  !     == LSPARENT : master node for Lattice Sum                       ==
  !     ==--------------------------------------------------------------==
  TYPE lrslspar_t
     LOGICAL :: rsnode
     LOGICAL :: rsparent
     LOGICAL :: lsnode
     LOGICAL :: lsparent
  END TYPE lrslspar_t
  TYPE (lrslspar_t) :: lrslspar
  !     ==--------------------------------------------------------------==
  !     == NRSNODES : Number of Real Space nodes                        ==
  !     == NLSNODES : Number of Lattice sum nodes                       ==
  !     == RSPOS    : CPU index within the Real Space nodes             ==
  !     == LSPOS    : CPU index within the Lattice Sum nodes            ==
  !     == RSGRP    : Real Space Group identifyer                       ==
  !     == LSGRP    : Lattice Sum Group identifyer                      ==
  !     ==--------------------------------------------------------------==
  TYPE irslspar_t
     INTEGER :: nrsnodes
     INTEGER :: nlsnodes
     INTEGER :: rspos
     INTEGER :: lspos
     INTEGER :: rsgrp
     INTEGER :: lsgrp
  END TYPE irslspar_t
  TYPE (irslspar_t) :: irslspar
END MODULE mm_parallel
