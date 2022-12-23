MODULE rscvp_utils
  USE cnst_dyn,                        ONLY: cv_vel,&
                                             ltcglobal,&
                                             ncolvar
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: paral
  USE system,                          ONLY: cntl,&
                                             cntr

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rscvp

CONTAINS

  ! ==================================================================
  SUBROUTINE rscvp(temp1,temp2,tempp,velp)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: temp1, temp2, tempp, &
                                                velp(:,:,:)

    INTEGER                                  :: ia, is
    REAL(real_8)                             :: alfap

! Variables
! ==--------------------------------------------------------------==
! ==  Dynamical rescaling factor (tempw/tempp), where tempp is    ==
! ==  calculated every step                                       ==
! ==--------------------------------------------------------------==

    IF (cntl%tcp) THEN
       IF (tempp.GT.temp1.OR.tempp.LT.temp2.AND.tempp.NE.0._real_8) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,F13.5,A,F13.5)') '  RSCVP| RESCALING IONIC '&
               // 'TEMP FROM ',tempp, ' TO ',cntr%tempw
          alfap=SQRT(cntr%tempw/tempp)
          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                velp(1,ia,is)=alfap*velp(1,ia,is)
                velp(2,ia,is)=alfap*velp(2,ia,is)
                velp(3,ia,is)=alfap*velp(3,ia,is)
             ENDDO
          ENDDO
          IF (ltcglobal) THEN
             DO ia = 1,ncolvar
                cv_vel(ia)=alfap*cv_vel(ia)
             ENDDO
          ENDIF
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rscvp
  ! ==================================================================

END MODULE rscvp_utils
