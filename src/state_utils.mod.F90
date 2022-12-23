MODULE state_utils
  USE cnst,                            ONLY: uimag
  USE cppt,                            ONLY: indzs,&
                                             nzhs
  USE fft,                             ONLY: inzs,&
                                             jgw,&
                                             llr1,&
                                             nzfs
  USE geq0mod,                         ONLY: geq0
  USE kinds,                           ONLY: real_8
  USE system,                          ONLY: ncpw

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: copy_to_re
  PUBLIC :: copy_re_to_re
  PUBLIC :: copy_im_to_re
  PUBLIC :: copy_re_to_im
  PUBLIC :: copy_im_to_im
  PUBLIC :: set_psi_im_state_r
  PUBLIC :: set_psi_re_state_r
  PUBLIC :: set_psi_1_state_g
  PUBLIC :: set_psi_2_states_g
  PUBLIC :: set_psi_1_state_g_kpts
  PUBLIC :: zero_wfn
  PUBLIC :: add_wfn

CONTAINS

  ! ==================================================================
  SUBROUTINE copy_to_re(n,a,b)
    ! ==================================================================
    INTEGER                                  :: n
    REAL(real_8)                             :: a(:)
    COMPLEX(real_8)                          :: b(:)

    INTEGER                                  :: i

    !$omp parallel do
    DO i=1,n
       b(i)=CMPLX(REAL(a(i),KIND=real_8),AIMAG(b(i)),KIND=real_8)
    ENDDO
  END SUBROUTINE copy_to_re
  ! ==================================================================
  SUBROUTINE copy_re_to_re(n,a,b)
    ! ==================================================================
    INTEGER                                  :: n
    COMPLEX(real_8)                          :: a(:), b(:)

    INTEGER                                  :: i

    !$omp parallel do
    DO i=1,n
       b(i)=CMPLX(REAL(a(i),KIND=real_8),AIMAG(b(i)),KIND=real_8)
    ENDDO
  END SUBROUTINE copy_re_to_re
  ! ==================================================================
  SUBROUTINE copy_im_to_re(n,a,b)
    ! ==================================================================
    INTEGER                                  :: n
    COMPLEX(real_8)                          :: a(:), b(:)

    INTEGER                                  :: i

    !$omp parallel do
    DO i=1,n
       b(i)=CMPLX(AIMAG(a(i)),AIMAG(b(i)),KIND=real_8)
    ENDDO
  END SUBROUTINE copy_im_to_re
  ! ==================================================================
  SUBROUTINE copy_re_to_im(n,a,b)
    ! ==================================================================
    INTEGER                                  :: n
    COMPLEX(real_8)                          :: a(:), b(:)

    INTEGER                                  :: i

    !$omp parallel do
    DO i=1,n
       b(i)=CMPLX(REAL(b(i),KIND=real_8),REAL(a(i),KIND=real_8),KIND=real_8)
    ENDDO
  END SUBROUTINE copy_re_to_im
  ! ==================================================================
  SUBROUTINE copy_im_to_im(n,a,b)
    ! ==================================================================
    INTEGER                                  :: n
    COMPLEX(real_8)                          :: a(:), b(:)

    INTEGER                                  :: i

    !$omp parallel do
    DO i=1,n
       b(i)=CMPLX(REAL(b(i),KIND=real_8),AIMAG(a(i)),KIND=real_8)
    ENDDO
  END SUBROUTINE copy_im_to_im

  ! ==================================================================
  SUBROUTINE set_psi_im_state_r(alpha,rswf,beta,psi)
    ! ==================================================================
    ! == Set Psi with the imaginary part of the *real* space state    ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: alpha, rswf(:), beta, psi(:)

    INTEGER                                  :: ir

    !$omp parallel do private(IR)
    DO ir=1,llr1
       psi(ir)=beta*psi(ir)+alpha*AIMAG(rswf(ir))
    ENDDO
    ! ==--------------------------------------------------------------==
  END SUBROUTINE set_psi_im_state_r

  ! ==================================================================
  SUBROUTINE set_psi_re_state_r(alpha,rswf,beta,psi)
    ! ==================================================================
    ! == Set Psi with the real part of the *real* space state         ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: alpha, rswf(:), beta, psi(:)

    INTEGER                                  :: ir

    !$omp parallel do private(IR)
    DO ir=1,llr1
       psi(ir)=beta*psi(ir)+alpha*REAL(rswf(ir))
    ENDDO
    ! ==--------------------------------------------------------------==
  END SUBROUTINE set_psi_re_state_r

  ! ==================================================================
  SUBROUTINE set_psi_1_state_g(alpha,c1,psi)
    ! ==================================================================
    ! == Set Psi with one state                                       ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: alpha, c1(:), psi(:)

    COMPLEX(real_8), PARAMETER               :: zone = (1.0_real_8,0.0_real_8)

    INTEGER                                  :: ig

    IF (alpha.EQ.zone) THEN
       !CDIR NODEP
       !ocl novrec
#ifdef __SR8000
       !poption parallel, tlocal(IG)
#endif
       !$omp parallel do private(IG)
       DO ig=1,jgw
          psi(nzfs(ig))=c1(ig)
          psi(inzs(ig))=CONJG(c1(ig))
       ENDDO
       IF (geq0) psi(nzfs(1))=c1(1)
    ELSE
       !CDIR NODEP
       !ocl novrec
#ifdef __SR8000
       !poption parallel, tlocal(IG)
#endif
       !$omp parallel do private(IG)
       DO ig=1,jgw
          psi(nzfs(ig))=alpha*c1(ig)
          psi(inzs(ig))=alpha*CONJG(c1(ig))
       ENDDO
       IF (geq0) psi(nzfs(1))=alpha*c1(1)
    ENDIF
    ! ==--------------------------------------------------------------==
  END SUBROUTINE set_psi_1_state_g

  ! ==================================================================
  SUBROUTINE set_psi_2_states_g(c1,c2,psi)
    ! ==================================================================
    ! == Set Psi with two states                                      ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c1(:), c2(:), psi(:)

    INTEGER                                  :: ig

!CDIR NODEP
!ocl novrec

    !$omp parallel do private(IG)
    DO ig=1,jgw
       psi(nzfs(ig))=c1(ig)+uimag*c2(ig)
       psi(inzs(ig))=CONJG(c1(ig))+uimag*CONJG(c2(ig))
    ENDDO
    IF (geq0) psi(nzfs(1))=c1(1)+uimag*c2(1)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE set_psi_2_states_g

  ! ==================================================================
  SUBROUTINE set_psi_1_state_g_kpts(alpha,c1,psi)
    ! ==================================================================
    ! == Set Psi with one state (k-points version)                    ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: alpha, c1(:), psi(:)

    COMPLEX(real_8), PARAMETER :: zone = (1.0_real_8,0.0_real_8) 

    INTEGER                                  :: ig

    IF (alpha.EQ.zone) THEN
       !ocl novrec
       !$omp parallel do private(IG)
       DO ig=1,ncpw%ngw
          ! The +g component of the state
          psi(nzhs(ig))=c1(ig)
          ! The -g component of the state
          psi(indzs(ig))=c1(ig+ncpw%ngw)
       ENDDO
       IF (geq0) psi(nzhs(1))=c1(1)
    ELSE
       !ocl novrec
       !$omp parallel do private(IG)
       DO ig=1,ncpw%ngw
          ! The +g component of the state
          psi(nzhs(ig))=alpha*c1(ig)
          ! The -g component of the state
          psi(indzs(ig))=alpha*c1(ig+ncpw%ngw)
       ENDDO
       IF (geq0) psi(nzhs(1))=alpha*c1(1)
    ENDIF
    ! ==--------------------------------------------------------------==
  END SUBROUTINE set_psi_1_state_g_kpts

  ! ==================================================================
  SUBROUTINE zero_wfn(m,n,x,ldx)
    ! ==================================================================
    ! == Zero a wavefunction                                          ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: m, n
    COMPLEX(real_8)                          :: x(:,:)
    INTEGER                                  :: ldx

    COMPLEX(real_8), PARAMETER :: zzero = (0.0_real_8,0.0_real_8) 

    INTEGER                                  :: i, j

! ==--------------------------------------------------------------==

    !$omp parallel do private(J)
    DO j=1,n
       DO i=1,m
          x(i,j) = zzero
       ENDDO
    ENDDO
    !$omp end parallel do
    ! ==--------------------------------------------------------------==
  END SUBROUTINE zero_wfn

  ! ==================================================================
  SUBROUTINE add_wfn(m,n,alpha,x,ldx,y,ldy)
    ! ==================================================================
    ! == Perform an axpy like using the scalar alpha and the          ==
    ! == wavefunctions X and Y                                        ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: m, n
    COMPLEX(real_8)                          :: alpha, x(:,:)
    INTEGER                                  :: ldx
    COMPLEX(real_8)                          :: y(:,:)
    INTEGER                                  :: ldy

    COMPLEX(real_8), PARAMETER :: zone = (1.0_real_8,0.0_real_8) 

    INTEGER                                  :: i, j

    IF (alpha.NE.zone) THEN
       !$omp parallel do private(J,I)
       DO j=1,n
          DO i=1,m
             y(i,j) = y(i,j) + alpha * x(i,j)
          ENDDO
       ENDDO
       !$omp end parallel do
    ELSE
       !$omp parallel do private(J,I)
       DO j=1,n
          DO i=1,m
             y(i,j) = y(i,j) + x(i,j)
          ENDDO
       ENDDO
       !$omp end parallel do
    ENDIF
    ! ==--------------------------------------------------------------==
  END SUBROUTINE add_wfn

END MODULE state_utils
