MODULE prcp
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==================================================================
  ! == PARRINELLO-RAHMAN CONSTANT PRESSURE INCLUDE FILE             ==
  ! ==--------------------------------------------------------------==
  ! == MPRCP      NUMBER OF real(8) :: IN THE COMMON                    ==
  ! ==--------------------------------------------------------------==
  INTEGER, PARAMETER :: mprcp=32 
  ! ==--------------------------------------------------------------==
  ! == AKIN       CONSTANT CUTOFF OPTION                            ==
  ! == SKIN       CONSTANT CUTOFF OPTION                            ==
  ! == ECKIN      CONSTANT CUTOFF OPTION                            ==
  ! ==            APPLY A CUTOFF FUNCTION TO THE KINETIC ENERGY TERM==
  ! ==            IN ORDER TO SIMULATE CONSTANT CUTOFF DYNAMICS     ==
  ! == GAKIN      =AKIN/TPIBA2                                      ==
  ! == GSKIN      =SKIN/TPIBA2                                      ==
  ! == GCKIN      =ECKIN/TPIBA2                                     ==
  ! == OMEGA0     VOLUME OF THE INITIAL CELL                        ==
  ! == TISOT      .TRUE. IF ISOTROPIC CELL                          ==
  ! == TZFLEX     .TRUE. IF Z-DIRECTION FLEXIBLE CELL               ==
  ! == CELLRF(6)  REFERENCE CELL                                    ==
  ! == DRUCK      REQUIRED PRESSURE                                 ==
  ! == STENS(3,3) REQUIRED STRESS TENSOR                            ==
  ! == HUNIT(3,3) UNIT VECTORS OF CELL COORDINATE SYSTEM            ==
  ! ==--------------------------------------------------------------==

  ! ==--------------------------------------------------------------==

  ! ==================================================================

  TYPE :: prcp_com_t
     REAL(real_8) :: akin
     REAL(real_8) :: skin
     REAL(real_8) :: eckin
     REAL(real_8) :: gakin
     REAL(real_8) :: gskin
     REAL(real_8) :: gckin
     REAL(real_8) :: omega0
     REAL(real_8) :: cellrf(6)
     REAL(real_8) :: druck
     REAL(real_8) :: stens(3,3)
     REAL(real_8) :: hunit(3,3)
  END TYPE prcp_com_t
  TYPE(prcp_com_t) :: prcp_com
  TYPE :: prcpl_t
     LOGICAL :: tisot
     LOGICAL :: tzflex
  END TYPE prcpl_t
  TYPE(prcpl_t) :: prcpl

END MODULE prcp
