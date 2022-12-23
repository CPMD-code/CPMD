MODULE ekinpp_utils
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE pimd,                            ONLY: pma0s
  USE rmas,                            ONLY: rmass
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ekinpp
  PUBLIC :: s_ekinpp

CONTAINS

  ! ==================================================================
  SUBROUTINE ekinpp(ekinp,velp)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: ekinp, velp(:,:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'ekinpp'

    INTEGER                                  :: isub, j, k
    REAL(real_8)                             :: const

! Variables
! ==--------------------------------------------------------------==
! ==  CALCULATE KINETIC ENERGY OF THE IONS                        ==
! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)

    ekinp=0._real_8
    DO k=1,ions1%nsp
       const=0.5_real_8*rmass%pma(k)
       DO j=1,ions0%na(k)
          ekinp=ekinp+const*(velp(1,j,k)**2+&
               VELP(2,J,K)**2+&
               VELP(3,J,K)**2) 
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    RETURN
  END SUBROUTINE ekinpp
  ! ==================================================================
  SUBROUTINE s_ekinpp(ekinp,velp,ip)
    REAL(real_8)                             :: ekinp, velp(:,:,:)
    INTEGER                                  :: ip

    CHARACTER(*), PARAMETER                  :: procedureN = 's_ekinpp'

    INTEGER                                  :: ia, is, isub
    REAL(real_8)                             :: const

! Variables
! ==--------------------------------------------------------------==
! ==  CALCULATE KINETIC ENERGY OF THE IONS                        ==
! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    ekinp=0._real_8
    DO is=1,ions1%nsp
       const=0.5_real_8*pma0s(is,ip)
       DO ia=1,ions0%na(is)
          ekinp=ekinp+const*(velp(1,ia,is)**2+&
               VELP(2,IA,IS)**2+&
               VELP(3,IA,IS)**2) 
       ENDDO
    ENDDO
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE s_ekinpp
  ! ==================================================================

END MODULE ekinpp_utils
