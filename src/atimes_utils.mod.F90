MODULE atimes_utils
  USE atimesmod,                       ONLY: atimes_eval
  USE error_handling,                  ONLY: stopgm
  USE hpsi_utils,                      ONLY: give_scr_hpsi
  USE kinds,                           ONLY: real_8
  USE projv_utils,                     ONLY: projv
  USE pslo,                            ONLY: pslo_com
  USE spin,                            ONLY: clsd
  USE system,                          ONLY: fpar,&
                                             nkpt

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: atimes
  PUBLIC :: give_scr_atimes

CONTAINS

  ! ==================================================================
  SUBROUTINE atimes(cin,cout,c0,vpot,psi,nstate,f)
    ! ==--------------------------------------------------------------==
    ! == ROUTINE TO APPLY (H-E) TO PSI
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: cin(nkpt%ngwk), &
                                                cout(nkpt%ngwk)
    REAL(real_8)                             :: vpot(fpar%nnr1,clsd%nlsd)
    COMPLEX(real_8)                          :: psi(fpar%nnr1)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(nkpt%ngwk,nstate)
    REAL(real_8)                             :: f(nstate)

    IF (pslo_com%tivan) CALL stopgm('ATIMES','TIVAN NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    CALL stopgm('ATIMES','fix me: wrong dummy args + provide a regtest',&
         __LINE__,__FILE__)
    !fix me call hpsi(cin,cout,sc0,vpot,psi,1,ikind_atimes,clsd%nlsd)
    CALL dscal(nkpt%ngwk*2,-1._real_8,cout,1)
    CALL daxpy(nkpt%ngwk*2,-atimes_eval,cin,1,cout,1)
    ! Project out valence states 
    CALL projv(cout,c0,f,nstate)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE atimes
  ! ==================================================================
  SUBROUTINE give_scr_atimes(latimes,tag,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: latimes
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate

! ==--------------------------------------------------------------==

    CALL give_scr_hpsi(latimes,tag,nstate)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_atimes
  ! ==================================================================

END MODULE atimes_utils
