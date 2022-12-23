MODULE parac
  IMPLICIT NONE

  ! ==--------------------------------------------------------------==
  ! == PARENT  : .TRUE if IP = SOURCE                               ==
  ! == IO_PARENT: .TRUE if IP = IO_SOURCE (cp_grp communicator)     ==
  ! ==--------------------------------------------------------------==
  TYPE :: paral_t
     LOGICAL :: parent = .FALSE.
     LOGICAL :: qmmmparent = .FALSE.
     LOGICAL :: qmnode = .FALSE.
     LOGICAL :: io_parent = .FALSE.
     LOGICAL :: cp_inter_io_parent = .FALSE.
  END TYPE paral_t
  TYPE(paral_t), SAVE :: paral
  ! ==--------------------------------------------------------------==
  ! == NCPUS   : value of NCPUS environment variable                ==
  ! == NPROC   : TOTAL NUMBER OF PROCESSORS                         ==
  ! == ME      : CPU index                                          ==
  ! == MEPOS   : equal to ME                                        ==
  ! == SOURCE  : MASTER CPU index                                   ==
  ! == IO_SOURCE: MASTER CPU index for IO (cp_grp communicator)     ==
  ! == IGEQ0   : PROCESSOR INDEX OF G=0 COMPONENT                   ==
  ! == ALLGRP  : GIVEN VALUE FOR BROADCAST COMMUNICATION            ==
  ! == CP_GRP  : CPMD communicator (should be use instead of        ==
  ! MPI_COMM_WORLD)                                    ==
  ! == CP_NPROC: Nbr processes in the cpmd communicator             ==
  ! == CP_ME   : id of the processes in the cpmd communicator       == 
  ! == NHRAYS  : number of rays for the processor for the density   ==
  ! == NGRAYS  : number of rays for the processor for the wavefunc. ==
  ! == CP_INTER_GRP: CPMD group communicator                        ==
  ! == CP_INTER_ME: CPMD group index                                ==
  ! == CP_NOGRP: number of CPMD group                               ==
  ! == CP_INTER_IO_SOURCE: MASTER intergroup index for IO           ==
  ! == LOC_GRP
  ! == LOC_ME
  ! == LOC_NPROC
  ! == LOC_INTER_GRP
  ! ==--------------------------------------------------------------==
  TYPE :: parai_t
     INTEGER :: ncpus = HUGE(0)
     INTEGER :: nproc = HUGE(0)
     INTEGER :: me = HUGE(0)
     INTEGER :: mepos = HUGE(0)
     INTEGER :: source = HUGE(0)
     INTEGER :: igeq0 = HUGE(0)
     INTEGER :: allgrp = HUGE(0)
     INTEGER :: nhrays = HUGE(0)
     INTEGER :: ngrays = HUGE(0)
     INTEGER :: qmmmnproc = HUGE(0)
     INTEGER :: qmmmgrp = HUGE(0)
     INTEGER :: qmmmsource = HUGE(0)
     INTEGER :: cp_grp = HUGE(0)
     INTEGER :: cp_nproc = HUGE(0)
     INTEGER :: cp_me = HUGE(0)
     INTEGER :: cp_inter_grp = HUGE(0)
     INTEGER :: cp_inter_me = HUGE(0)
     INTEGER :: cp_nogrp = HUGE(0)
     INTEGER :: io_source = HUGE(0)
     INTEGER :: cp_inter_io_source = HUGE(0)
     INTEGER :: loc_grp = HUGE(0)
     INTEGER :: loc_me = HUGE(0)
     INTEGER :: loc_nproc = HUGE(0)
     INTEGER :: loc_inter_grp = HUGE(0)
  END TYPE parai_t
  TYPE(parai_t), SAVE :: parai

END MODULE parac
