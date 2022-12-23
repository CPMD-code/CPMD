MODULE puttau_utils
  USE cotr,                            ONLY: cotc0,&
                                             lskcor
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: puttau
  PUBLIC :: gettau
  PUBLIC :: taucl
  PUBLIC :: getatm

CONTAINS

  ! ==================================================================
  SUBROUTINE puttau(tau,xpar)
    ! ==--------------------------------------------------------------==
    ! == PUT TAU0 NON FIXED(GIVEN BY LSKCOR) IN XPAR                  ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau(:,:,:), xpar(cotc0%nodim)

    INTEGER                                  :: ia, iat, is, k, l

    k=0
    iat=0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          iat=iat+1
          DO l=1,3
             IF (lskcor(l,iat).NE.0) THEN
                k=k+1
                xpar(k) = tau(l,ia,is)
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE puttau
  ! ==================================================================
  SUBROUTINE gettau(tau,xpar)
    ! ==--------------------------------------------------------------==
    ! == GET XPAR IN TAU0                                             ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau(:,:,:), xpar(cotc0%nodim)

    INTEGER                                  :: ia, iat, is, k, l

    k=0
    iat=0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          iat=iat+1
          DO l=1,3
             IF (lskcor(l,iat).NE.0) THEN
                k=k+1
                tau(l,ia,is) = xpar(k)
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gettau
  ! ==================================================================
  SUBROUTINE taucl(tau)
    ! ==--------------------------------------------------------------==
    ! == CLEAR COORDINATES WHERE LSKCOR(L,IAT).EQ.0                   ==
    ! ==--------------------------------------------------------------==
    ! Arguemnts
    REAL(real_8)                             :: tau(:,:,:)

    INTEGER                                  :: ia, iat, is, k, l

    k=0
    iat=0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          iat=iat+1
          DO l=1,3
             IF (lskcor(l,iat).EQ.0) tau(l,ia,is)=0.0_real_8
          ENDDO
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE taucl
  ! ==================================================================
  SUBROUTINE getatm(nu,tau,acor)
    ! ==--------------------------------------------------------------==
    ! == Retrive Coordiantes of Atom IAT                              ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nu
    REAL(real_8)                             :: tau(:,:,:), acor(3)

    INTEGER                                  :: ia, iat, is

    iat=0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          iat=iat+1
          IF (iat.EQ.nu) THEN
             acor(1) = tau(1,ia,is)
             acor(2) = tau(2,ia,is)
             acor(3) = tau(3,ia,is)
          ENDIF
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE getatm
  ! ==================================================================


END MODULE puttau_utils
