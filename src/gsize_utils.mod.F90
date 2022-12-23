MODULE gsize_utils
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: gsize
  PUBLIC :: gnodim

CONTAINS

  ! ==================================================================
  SUBROUTINE gsize(fion,gnmax,gnorm)
    ! ==--------------------------------------------------------------==
    ! == CALCULATE NORM OF GRADIENT (GNORM) AND                       ==
    ! == MAXIMUM COMPONENT (GNMAX)                                    ==
    ! == INPUT: FION(3,maxsys%nax,maxsys%nsx) NUCLEAR GRADIENTS                     ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: fion(:,:,:), gnmax, gnorm

    INTEGER                                  :: ia, is, k, ntot

    ntot=0
    gnmax=0.0_real_8
    gnorm=0.0_real_8
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          DO k=1,3
             ntot=ntot+1
             IF (ABS(fion(k,ia,is)).GT.gnmax) gnmax=ABS(fion(k,ia,is))
             gnorm=gnorm+fion(k,ia,is)*fion(k,ia,is)
          ENDDO
       ENDDO
    ENDDO
    gnorm=SQRT(gnorm/REAL(ntot,kind=real_8))
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gsize
  ! ==================================================================
  SUBROUTINE gnodim(dxpar,nodim,gnmax,gnorm)
    ! ==--------------------------------------------------------------==
    ! == CALCULATE NORM OF GRADIENT (GNORM) AND                       ==
    ! == MAXIMUM COMPONENT (GNMAX)                                    ==
    ! == INPUT: DXPAR(NODIM) GRADIENTS OF FREEDOM DEGREES             ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nodim
    REAL(real_8)                             :: dxpar(nodim), gnmax, gnorm

    INTEGER                                  :: i

    gnmax=0.0_real_8
    gnorm=0.0_real_8
    DO i=1,nodim
       IF (ABS(dxpar(i)).GT.gnmax) gnmax=ABS(dxpar(i))
       gnorm=dxpar(i)*dxpar(i)+gnorm
    ENDDO
    gnorm=SQRT(gnorm/nodim)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gnodim
  ! ==================================================================

END MODULE gsize_utils
