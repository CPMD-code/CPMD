MODULE fixcom_utils
  USE cotr,                            ONLY: cotc0,&
                                             lskcor,&
                                             patot,&
                                             pmall
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: fixcom

CONTAINS

  ! ==================================================================
  SUBROUTINE fixcom(disp)
    ! ==--------------------------------------------------------------==
    ! == FOR THE OPTION FIX COM (FIX THE CENTER OF MASS)              ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: disp(cotc0%nodim)

    INTEGER                                  :: ia, iat, is, k, l
    REAL(real_8)                             :: cdmm(3)

! Variables
! ==--------------------------------------------------------------==
! ==  COMPUTE AND CORRECT FOR DISPLACEMENT OF CENTER OF MASS      ==
! ==--------------------------------------------------------------==

    IF (.NOT.cotc0%lfcom) RETURN
    cdmm(1)=0._real_8
    cdmm(2)=0._real_8
    cdmm(3)=0._real_8
    k=0
    iat=0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          iat=iat+1
          DO l=1,3
             IF (lskcor(l,iat).NE.0) THEN
                k=k+1
                cdmm(l)=cdmm(l)+disp(k)*pmall(iat)
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    cdmm(1)=cdmm(1)/patot
    cdmm(2)=cdmm(2)/patot
    cdmm(3)=cdmm(3)/patot
    k=0
    iat=0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          iat=iat+1
          DO l=1,3
             IF (lskcor(l,iat).NE.0) THEN
                k=k+1
                disp(k)=disp(k)-cdmm(l)
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE fixcom
  ! ==================================================================

END MODULE fixcom_utils
