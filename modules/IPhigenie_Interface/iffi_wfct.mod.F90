MODULE iffi_wfct
    USE kinds,        ONLY: real_8
  IMPLICIT NONE

    COMPLEX(real_8), ALLOCATABLE             :: c0(:,:,:), c2(:,:), gde(:), &
                                                pme(:), sc0(:,:,:)
    REAL(real_8), ALLOCATABLE                :: eigv(:), vpp(:)
 
  
    COMPLEX(real_8), ALLOCATABLE             :: psi(:,:)
    REAL(real_8), ALLOCATABLE                :: rhoe(:,:)
    REAL(real_8), ALLOCATABLE                :: scr(:)
END MODULE iffi_wfct
  