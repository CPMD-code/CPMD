MODULE atom
  USE kinds,                           ONLY: real_8
  USE system,                          ONLY: maxsp

  IMPLICIT NONE

  ! ==================================================================
  ! == INCLUDE FILE FOR PSEUDOPOTENTIAL INFORMATION                 ==
  ! ==================================================================
  CHARACTER(len=200) :: ecpfiles(maxsp) = ''
  TYPE :: ecpf_t
     CHARACTER(len=200) :: ecpfiles(maxsp) = ''
  END TYPE ecpf_t
  TYPE(ecpf_t), SAVE :: ecpf
  ! ==--------------------------------------------------------------==



  REAL(real_8), ALLOCATABLE, SAVE :: gnl(:,:,:)
  REAL(real_8), ALLOCATABLE, SAVE :: rps(:,:,:)
  REAL(real_8), ALLOCATABLE, SAVE :: rw(:,:)
  REAL(real_8), ALLOCATABLE, SAVE :: rv(:,:)
  REAL(real_8), ALLOCATABLE, SAVE :: vr(:,:,:)

  TYPE :: atom_common_t
     REAL(real_8) :: clogvp(maxsp)=HUGE(0.0_real_8)
     INTEGER :: meshvp(maxsp)=HUGE(0)
     REAL(real_8) :: clogw(maxsp)=HUGE(0.0_real_8)
     INTEGER :: meshw(maxsp)=HUGE(0)
  END TYPE atom_common_t
  TYPE(atom_common_t), SAVE :: atom_common
  ! ==================================================================

  TYPE :: patom1_t
     LOGICAL :: pconf=.FALSE.
  END TYPE patom1_t
  TYPE(patom1_t), SAVE :: patom1


  TYPE :: patom2_t
     REAL(real_8) :: palpha(maxsp)=HUGE(0.0_real_8)
     REAL(real_8) :: prc(maxsp)=HUGE(0.0_real_8)
  END TYPE patom2_t
  TYPE(patom2_t), SAVE :: patom2

END MODULE atom
