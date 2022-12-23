MODULE vgsortho_utils
  USE dotp_utils,                      ONLY: dotp
  USE geq0mod,                         ONLY: geq0
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: vgsortho

CONTAINS

  ! ==================================================================
  SUBROUTINE vgsortho(c0,sc0,ngw,nfirst,nlast)
    ! ==--------------------------------------------------------------==
    ! ==  Gram-Schmidt orthogonalization for Vanderbilt pp            ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ngw
    COMPLEX(real_8)                          :: sc0(ngw,*), c0(ngw,*)
    INTEGER                                  :: nfirst, nlast

    INTEGER                                  :: i, isub, j
    REAL(real_8)                             :: anorm, s

    CALL tiset('  VGSORTHO',isub)
    DO i=nfirst,nlast
       DO j=1,i-1
          s=-dotp(ngw,sc0(:,i),c0(:,j))
          CALL mp_sum(s,parai%allgrp)
          CALL daxpy(2*ngw,s,c0(1,j),1,c0(1,i),1)
          CALL daxpy(2*ngw,s,sc0(1,j),1,sc0(1,i),1)
       ENDDO
       s=dotp(ngw,sc0(:,i),c0(:,i))
       CALL mp_sum(s,parai%allgrp)
       anorm=1._real_8/SQRT(s)
       CALL dscal(2*ngw,anorm,c0(1,i),1)
       CALL dscal(2*ngw,anorm,sc0(1,i),1)
       IF (geq0) THEN
          c0(1,i)=CMPLX(REAL(c0(1,i)),0.0_real_8,kind=real_8)
          sc0(1,i)=CMPLX(REAL(sc0(1,i)),0.0_real_8,kind=real_8)
       ENDIF
    ENDDO
    CALL tihalt('  VGSORTHO',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE vgsortho
  ! ==================================================================

END MODULE vgsortho_utils
