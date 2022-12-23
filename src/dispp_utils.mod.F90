MODULE dispp_utils
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE system,                          ONLY: maxsp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: dispp

CONTAINS

  ! ==================================================================
  SUBROUTINE dispp(taup,taui,disa)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: taup(:,:,:), taui(:,:,:), disa

    INTEGER                                  :: ia, is, k
    REAL(real_8)                             :: dis(maxsp), r2(3)

! Variables
! ==--------------------------------------------------------------==
! ==  MEAN SQUARE DISPLACEMENT OF DIFFERENT SPECIES IONS          ==
! ==--------------------------------------------------------------==

    disa=0._real_8
    DO is=1,ions1%nsp
       dis(is)=0._real_8
       DO k=1,3
          r2(k)=0._real_8
          DO ia=1,ions0%na(is)
             r2(k)=r2(k)+(taup(k,ia,is)-taui(k,ia,is))**2
          ENDDO
          dis(is)=dis(is)+r2(k)
       ENDDO
       disa=disa+dis(is)
    ENDDO
    disa=disa/REAL(ions1%nat,kind=real_8)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE dispp
  ! ==================================================================

END MODULE dispp_utils
