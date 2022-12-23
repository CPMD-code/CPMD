MODULE projv_utils
  USE dotp_utils,                      ONLY: dotp
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE pola,                            ONLY: tollocc
  USE system,                          ONLY: ncpw,&
                                             nkpt

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: projv

CONTAINS

  ! ==================================================================
  SUBROUTINE projv(cproj,c0,f,nstate)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: cproj(nkpt%ngwk)
    REAL(real_8)                             :: f(nkpt%ngwk)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(nkpt%ngwk,nstate)

    COMPLEX(real_8)                          :: caux
    COMPLEX(real_8), EXTERNAL                :: zdotc
    INTEGER                                  :: k
    REAL(real_8)                             :: aux

    DO k = 1 , nstate
       ! Project out valence states
       IF (f(k).GT.tollocc) THEN
          IF (tkpts%tkpnt) THEN
             caux=zdotc(nkpt%ngwk,c0(1,k),1,cproj,1)
             CALL zaxpy(nkpt%ngwk,-caux,c0(1,k),1,cproj,1)
          ELSE
             aux=dotp(ncpw%ngw,c0(:,k),cproj)
             CALL daxpy(ncpw%ngw*2,-aux,c0(1,k),1,cproj,1)
          ENDIF
       ENDIF
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE projv

END MODULE projv_utils
