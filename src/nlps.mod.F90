MODULE nlps
  USE kinds,                           ONLY: real_8
  USE system,                          ONLY: maxsp

  IMPLICIT NONE

  ! ==================================================================
  ! == NON LOCAL PSEUDOPOTENTIALS PARAMETERS                        ==
  ! ==--------------------------------------------------------------==
  ! ==  IMAGP          = 1 IF TKPNT=.FALSE.                         ==
  ! ==       OTHERWISE = 2 (FNL AND DFNL ARE COMPLEX IN THIS CASE)  ==
  ! ==  NLM            TOTAL NUMBER OF NON-LOCAL ANGULAR MOMENTUM   ==
  ! ==  NDFNL          = NSTATE IF NOPARALLEL                       ==
  ! ==       OTHERWISE = NST12(MEPOS,2)-NST12(MEPOS,1)+1            ==
  ! ==--------------------------------------------------------------==
  INTEGER :: nlm ,imagp,ndfnl
  ! ==--------------------------------------------------------------==
  ! == WSG(NSX,NHX)                                                 ==
  ! == RGH(NHX,NSX)                                                 ==
  ! == WGH(NHX,NSX)                                                 ==
  ! == NGHTOL(NHX,NSX) For each projector, gives L ang. mom. q. n.  ==
  ! == NGHCOM(NHX,NSX) For each projector, gives LP                 ==
  ! == NGH(MAXSP) Number of projectors per specie                   ==
  ! ==--------------------------------------------------------------==
  REAL(real_8), ALLOCATABLE :: wsg(:,:)
  REAL(real_8), ALLOCATABLE :: rgh(:,:)
  REAL(real_8), ALLOCATABLE :: wgh(:,:)


  INTEGER, ALLOCATABLE :: nghtol(:,:)
  INTEGER, ALLOCATABLE :: nghcom(:,:)



  ! IDUM is just for padding to avoid a misalignment

  ! ==================================================================
  TYPE :: nlps_com_t
     REAL(real_8) :: rmaxn(maxsp)
     INTEGER :: ngh(maxsp)
     INTEGER :: idum
  END TYPE nlps_com_t
  TYPE(nlps_com_t) :: nlps_com

END MODULE nlps
