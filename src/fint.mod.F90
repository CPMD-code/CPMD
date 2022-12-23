MODULE fint
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==================================================================
  ! == FREE ENERGY INCLUDE FILE                                     ==
  ! == FOR TROTTER APPROXIMATION                                    ==
  ! ==================================================================
  ! ==  BETAEL:  1/kT                                               ==
  ! ==  BETAP:   TROTTER FACTOR                                     ==
  ! ==  P:       BETAEL/BETAP                                       ==
  ! ==  TBOGO:   TRUE FOR THE BOGOLIUBOV CORRECTION                 ==
  ! ==  TTROT:   TRUE FOR THE TROTTER APPROXIMATION                 ==
  ! ==  TFRAL:   TRUE, COMPUTES ALL WANTED EIGENVALUES WITH         ==
  ! ==           THE REQUIRED ACCURACY                              ==
  ! ==--------------------------------------------------------------==


  ! ==================================================================
  ! ==  EMK2:    EXPONENTIATED KINETIC ENERGY IN G-SPACE            ==
  ! ==  ALM(IMAGP,NLM,NLM,NKPNT)                                    ==
  ! ==           OVERLAP MATRIX BETWEEN NL PP PROJECTORS            ==
  ! ==           IF TKPNT, ALM IS COMPLEX (IMAGP = 2)               ==
  ! ==  AFNL(IMAGP,NLM,2)                                           ==
  ! ==           IF TKPNT, AFNL IS COMPLEX (IMAGP = 2)              ==
  ! ==  BILN(IMAGP,NLM,2)                                           ==
  ! ==           IF TKPNT, BILN IS COMPLEX (IMAGP = 2)              ==
  ! ==  CNL:                                                        ==
  ! ==  NLPTR:   INDEX ARRAY GIVING IS, ISA, AND IV FOR EACH NLM    ==
  ! ==--------------------------------------------------------------==
  REAL(real_8), ALLOCATABLE, TARGET :: emk2(:,:)

  REAL(real_8), ALLOCATABLE :: alm(:,:,:,:)
  REAL(real_8), ALLOCATABLE :: afnl(:,:,:)
  REAL(real_8), ALLOCATABLE :: biln(:,:,:)

  INTEGER, ALLOCATABLE :: nlptr(:,:)

  COMPLEX(real_8), ALLOCATABLE :: cnl(:,:)

  ! ==--------------------------------------------------------------==
  ! == TO SAVE TIME WE CAN CHANGE THE VALUE OF B2LIMIT              ==
  ! == WHICH CAN ONLY DECREASE                                      ==
  ! ==--------------------------------------------------------------==
  ! == MAXTROT: MAXIMAL NUMBER OF CHOICE FOR B2LIMIT                ==
  ! == NTABTROT: NUMBER OF B2LIMIT CHOICE                           ==
  ! == ITROT: CURRENT INDEX IN THE ARRAY B2TROT                     ==
  ! == DENSTROT(MAXTROT): DENSITY THRESHOLD ARRAY                   ==
  ! == B2TROT(MAXTROT):  B2LIMIT ARRAY                              ==
  ! ==--------------------------------------------------------------==
  INTEGER,PARAMETER :: maxtrot=10


  ! ==--------------------------------------------------------------==
  ! == TO SAVE TIME WE CAN CHANGE THE VALUE OF BETAP                ==
  ! == WHICH CAN ONLY DECREASE                                      ==
  ! ==--------------------------------------------------------------==
  ! == MAXBETAP: MAXIMAL NUMBER OF CHOICE FOR BETAP                 ==
  ! == NTABBETAP: NUMBER OF BETAP CHOICE                            ==
  ! == IBETAP: CURRENT INDEX IN THE ARRAY TABBETAP                  ==
  ! == DENSBETAP(MAXBETAP): DENSITY THRESHOLD ARRAY                 ==
  ! == TABBETAP(MAXTROT):  TABBETAP ARRAY                           ==
  ! ==--------------------------------------------------------------==
  INTEGER,PARAMETER :: maxbetap=10


  ! ==================================================================

  TYPE :: fint1_t
     REAL(real_8) :: betael
     REAL(real_8) :: betap
     REAL(real_8) :: p
     LOGICAL :: tbogo
     LOGICAL :: ttrot
     LOGICAL :: tfral
  END TYPE fint1_t
  TYPE(fint1_t) :: fint1
  TYPE :: fint4_t
     REAL(real_8) :: denstrot(maxtrot)
     REAL(real_8) :: b2trot(maxtrot)
     INTEGER :: ntabtrot
     INTEGER :: itrot
  END TYPE fint4_t
  TYPE(fint4_t) :: fint4
  TYPE :: fint5_t
     REAL(real_8) :: densbetap(maxbetap)
     REAL(real_8) :: tabbetap(maxbetap)
     INTEGER :: ntabbetap
     INTEGER :: ibetap
  END TYPE fint5_t
  TYPE(fint5_t) :: fint5

END MODULE fint
