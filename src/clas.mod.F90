MODULE clas
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! =================================================================
  ! == CLASSICAL PART                                              ==
  ! =================================================================
  INTEGER, PARAMETER :: maxclt=100 
  ! ==-------------------------------------------------------------==
  LOGICAL :: tclas
  ! ==-------------------------------------------------------------==
  REAL(real_8), ALLOCATABLE :: clasc(:,:)
  REAL(real_8), ALLOCATABLE :: clasv(:,:)
  REAL(real_8), ALLOCATABLE :: clasf(:,:)
  REAL(real_8), ALLOCATABLE :: clasc0(:,:)

  REAL(real_8), ALLOCATABLE :: delclasc(:,:)
  REAL(real_8), ALLOCATABLE :: clasfold(:,:)

  INTEGER, ALLOCATABLE :: cltyp(:)



  ! ==-------------------------------------------------------------==

  ! ==-------------------------------------------------------------==

  ! ==-------------------------------------------------------------==
  CHARACTER(len=4) :: clab(maxclt)
  ! ==-------------------------------------------------------------==
  INTEGER :: myclat1,myclat2,ntrac,ndefo,imovc
  ! =================================================================
  INTEGER, PARAMETER :: noptcl=4 
  ! m added TDEFO
  ! =================================================================

  ! =================================================================
  ! ..R12 potential


  REAL(real_8), ALLOCATABLE :: pr12(:,:)
  REAL(real_8), ALLOCATABLE :: rcr12(:,:)


  ! ==-------------------------------------------------------------==
  REAL(real_8) :: rlrc
  ! =================================================================
  ! ..CLASSICAL STRESS
  REAL(real_8) :: claspres(3,3)
  ! =================================================================
  ! DEFORMATION

  ! =================================================================

  TYPE :: clas3_t
     INTEGER :: nclatom
     INTEGER :: ncltyp
     INTEGER :: ncrang(2,maxclt)
     INTEGER :: is_qm(maxclt)
  END TYPE clas3_t
  TYPE(clas3_t) :: clas3
  TYPE :: clas4_t
     REAL(real_8) :: clmas(maxclt)
     REAL(real_8) :: cstep(maxclt)
  END TYPE clas4_t
  TYPE(clas4_t) :: clas4
  TYPE :: clas7_t
     LOGICAL :: tfreeze
     LOGICAL :: twcoor
     LOGICAL :: tpff
     LOGICAL :: tftra
     LOGICAL :: tdefo
     LOGICAL :: cprint
     LOGICAL :: tmovc
  END TYPE clas7_t
  TYPE(clas7_t) :: clas7
  TYPE :: clas8_t
     REAL(real_8) :: cellcl(6)
     REAL(real_8) :: clomega
     REAL(real_8) :: cla1(3)
     REAL(real_8) :: cla2(3)
     REAL(real_8) :: cla3(3)
     REAL(real_8) :: clht(3,3)
     REAL(real_8) :: clhtm1(3,3)
  END TYPE clas8_t
  TYPE(clas8_t) :: clas8
  TYPE :: defy_t
     REAL(real_8) :: xcen
     REAL(real_8) :: ycen
     REAL(real_8) :: zcen
     REAL(real_8) :: expx
     REAL(real_8) :: expy
     REAL(real_8) :: expz
  END TYPE defy_t
  TYPE(defy_t) :: defy
  TYPE :: pote1_t
     LOGICAL :: t_r12
     INTEGER :: initr12
  END TYPE pote1_t
  TYPE(pote1_t) :: pote1

END MODULE clas
