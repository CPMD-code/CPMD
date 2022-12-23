! 
MODULE shop_rest_2
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! surface hopping data for RESTART 
  ! 
  ! 
  ! BO additional WF arrays (for RESTART option)
  ! 
  COMPLEX(real_8), ALLOCATABLE :: c0old(:,:)
  COMPLEX(real_8), ALLOCATABLE :: c0old2(:,:)

END MODULE shop_rest_2
