MODULE sd_wannier_utils
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: int_1,&
                                             int_2,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: dgive,&
                                             unitmx
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  !public :: sd_wannier
  PUBLIC :: give_scr_sdwann
  PUBLIC :: give_scr_exponentiate
  PUBLIC :: exponentiate
  !public :: xgradat0
  !public :: ofunc
  !public :: xyz_update

CONTAINS

  ! ==================================================================
  SUBROUTINE give_scr_sdwann(lwann,tag,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lwann
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate

    INTEGER                                  :: lexp, lupdate

    CALL give_scr_exponentiate(lexp,tag,nstate)
    lupdate=3*nstate*nstate
    lwann=MAX(lexp,lupdate)
    tag='MAX(LEXP,LUPDATE)'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_sdwann
  ! ==================================================================
  SUBROUTINE give_scr_exponentiate(lexp,tag,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lexp
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate

! ==--------------------------------------------------------------==

    lexp=2*nstate*nstate+4*nstate+2*25*nstate+3*nstate + 10
    tag='57*NSTATE+2*NSTATE*NSTATE'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_exponentiate
  ! ==================================================================
  SUBROUTINE exponentiate(gmat,rmat,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    REAL(real_8)                             :: rmat(nstate,nstate), &
                                                gmat(nstate,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'exponentiate'

    COMPLEX(real_8), ALLOCATABLE             :: eia(:), hgmat(:,:), work(:)
    INTEGER                                  :: i, ierr, info, j, k, lwork
    REAL(real_8)                             :: cei, ei, sei
    REAL(real_8), ALLOCATABLE                :: rwork(:)

! variables
! ==--------------------------------------------------------------==

    ALLOCATE(hgmat(nstate, nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(eia(nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! ..calculate eigenvalues of iG (this is Hermite schen)
    DO i=1,nstate
       DO j=i,nstate
          hgmat(i,j)=CMPLX(0._real_8,gmat(i,j),kind=real_8)
       ENDDO
    ENDDO

    lwork = MAX(1, 2*nstate-1) ! TODO OPTIMAL LWORK
    ALLOCATE(work(lwork),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(rwork(MAX(1, 3*nstate-2)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL zheev('V','U',nstate,hgmat,nstate,rmat,work,lwork,rwork,info)
    DEALLOCATE(work,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(rwork,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)

    DO i=1,nstate
       ei=rmat(i,1)
       cei=COS(ei)
       sei=SIN(ei)
       eia(i)=CMPLX(cei,sei,kind=real_8)
    ENDDO
    CALL zeroing(rmat)!,nstate*nstate)
#ifdef __ES 
    !$omp parallel do private(J,I,K)
    DO j=1,nstate
       DO i=1,nstate
          DO k=1,nstate
             rmat(i,j)=rmat(i,j)+REAL(eia(k)*hgmat(i,k)*CONJG(hgmat(j,&
                  k)))
          ENDDO
       ENDDO
    ENDDO
#else 
    DO k=1,nstate
       DO i=1,nstate
          DO j=1,nstate
             rmat(i,j)=rmat(i,j)+REAL(eia(k)*hgmat(i,k)*CONJG(hgmat(j,&
                  k)))
          ENDDO
       ENDDO
    ENDDO
#endif 
    ! ==--------------------------------------------------------------==
    DEALLOCATE(hgmat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(eia,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE exponentiate
  ! ==================================================================
  ! ==================================================================
  ! ==================================================================
  ! ==================================================================

END MODULE sd_wannier_utils


SUBROUTINE xyz_update(rmat,xyzmat,abc,ldx,nstate)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE wann , ONLY:wannc
  IMPLICIT NONE
  INTEGER                                    :: ldx
  REAL(real_8) :: xyzmat(2,ldx,ldx,wannc%nwanopt)
  INTEGER                                    :: nstate
  REAL(real_8)                               :: abc(nstate,nstate,3), &
                                                rmat(nstate,nstate)

  INTEGER                                    :: i, j, k

  DO i=1,wannc%nwanopt
     !$omp parallel do private(J,K)
     DO j=1,nstate
        DO k=1,nstate
           abc(j,k,1)=xyzmat(1,j,k,i)
        ENDDO
     ENDDO
     CALL dgemm('T','N',nstate,nstate,nstate,1._real_8,rmat,nstate,&
          abc(1,1,1),nstate,0._real_8,abc(1,1,2),nstate)
     CALL dgemm('N','N',nstate,nstate,nstate,1._real_8,abc(1,1,2),nstate,&
          rmat,nstate,0._real_8,abc(1,1,1),nstate)
     !$omp parallel do private(J,K)
     DO j=1,nstate
        DO k=1,nstate
           xyzmat(1,j,k,i)=abc(j,k,1)
        ENDDO
     ENDDO
     !$omp parallel do private(J,K)
     DO j=1,nstate
        DO k=1,nstate
           abc(j,k,1)=xyzmat(2,j,k,i)
        ENDDO
     ENDDO
     CALL dgemm('T','N',nstate,nstate,nstate,1._real_8,rmat,nstate,&
          abc(1,1,1),nstate,0._real_8,abc(1,1,2),nstate)
     CALL dgemm('N','N',nstate,nstate,nstate,1._real_8,&
          abc(1,1,2),nstate,rmat,nstate,0._real_8,abc(1,1,1),nstate)
     !$omp parallel do private(J,K)
     DO j=1,nstate
        DO k=1,nstate
           xyzmat(2,j,k,i)=abc(j,k,1)
        ENDDO
     ENDDO
  ENDDO
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE xyz_update


SUBROUTINE xgradat0(wt,gmat,xyzmat,ldx,nstate)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE wann , ONLY:wannc
  USE zeroing_utils,                   ONLY: zeroing
  IMPLICIT NONE
  INTEGER                                    :: wt, ldx
  COMPLEX(real_8)                            :: xyzmat(ldx,ldx,wannc%nwanopt)
  INTEGER                                    :: nstate
  REAL(real_8)                               :: gmat(nstate,nstate)

  COMPLEX(real_8)                            :: xa, xb
  INTEGER                                    :: i, j, k
  REAL(real_8)                               :: gcomp

  CALL zeroing(gmat)!,nstate*nstate)
  IF (wt.EQ.1) THEN
     ! ..Vanderbilt functional |M|^2
     !$omp parallel do private(I,J,K,GCOMP,XA) schedule(static,1)
     DO i=1,nstate-1
        DO j=i+1,nstate
           gcomp=0._real_8
           DO k=1,wannc%nwanopt
              xa=xyzmat(j,j,k)-xyzmat(i,i,k)
              gcomp=gcomp+REAL(4._real_8*CONJG(xyzmat(i,j,k))*wannc%wwei(k)*xa)
           ENDDO
           gmat(i,j)=gcomp
           gmat(j,i)=-gmat(i,j)
        ENDDO
     ENDDO
  ELSEIF (wt.EQ.2) THEN
     ! ..Resta functional Log(|M|^2)
     !$omp parallel do private(I,J,K,GCOMP,XA,XB) schedule(static,1)
     DO i=1,nstate-1
        DO j=i+1,nstate
           gcomp=0._real_8
           DO k=1,wannc%nwanopt
              xa=REAL(xyzmat(i,j,k)/xyzmat(i,i,k))
              xb=REAL(xyzmat(i,j,k)/xyzmat(j,j,k))
              gcomp=gcomp+2._real_8*wannc%wwei(k)*REAL(xb-xa)
           ENDDO
           gmat(i,j)=gcomp
           gmat(j,i)=-gmat(i,j)
        ENDDO
     ENDDO
  ENDIF
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE xgradat0



SUBROUTINE ofunc(wt,funcv,xyzmat,ldx,nstate)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE system , ONLY:parm
  USE parac, ONLY : paral,parai
  USE wann , ONLY:wannc
  IMPLICIT NONE
  INTEGER                                    :: wt
  REAL(real_8)                               :: funcv
  INTEGER                                    :: ldx
  COMPLEX(real_8)                            :: xyzmat(ldx,ldx,wannc%nwanopt)
  INTEGER                                    :: nstate

  INTEGER                                    :: i, k

  funcv=0._real_8
  IF (wt.EQ.1) THEN
     ! ..Vanderbilt functional |M|^2
     DO i=1,nstate
        DO k=1,wannc%nwanopt
           funcv=funcv+wannc%wwei(k)*REAL(CONJG(xyzmat(i,i,k))*xyzmat(i,i,&
                k))
        ENDDO
     ENDDO
  ELSEIF (wt.EQ.2) THEN
     ! ..Resta functional Log(|M|^2)
     DO i=1,nstate
        DO k=1,wannc%nwanopt
           funcv=funcv+(-wannc%wwei(k)*LOG(REAL(CONJG(xyzmat(i,i,k))*&
                xyzmat(i,i,k))))
        ENDDO
     ENDDO
  ENDIF
  funcv=funcv/parm%alat**2
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE ofunc

! ==================================================================
SUBROUTINE sd_wannier(rotmat,xyzmat,ldx,nstate)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE prng_utils, ONLY : repprngu_vec, repprngu, repprngg, prngparaskip
  USE system , ONLY:parm
  USE parac, ONLY : paral,parai
  USE wann , ONLY:wannc,wanni,wannr
  USE sd_wannier_utils, ONLY : exponentiate
  USE utils, ONLY : unitmx
  USE utils, ONLY : dgive
  IMPLICIT NONE
  INTEGER                                    :: ldx
  REAL(real_8) :: xyzmat(2,ldx,ldx,wannc%nwanopt)
  INTEGER                                    :: nstate
  REAL(real_8)                               :: rotmat(nstate,nstate)

  CHARACTER(*), PARAMETER                    :: procedureN = 'sd_wannier'

  INTEGER                                    :: i, idamax, ierr, imax, iopt, j
  REAL(real_8)                               :: gmax, ofun, pi, scal
  REAL(real_8), ALLOCATABLE                  :: aux(:,:,:), gmat(:,:), &
                                                rmat(:,:)

! ==--------------------------------------------------------------==
! ==--------------------------------------------------------------==

  ALLOCATE(aux(nstate,nstate,3),STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
       __LINE__,__FILE__)
  ALLOCATE(gmat(nstate,nstate),STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
       __LINE__,__FILE__)
  ALLOCATE(rmat(nstate,nstate),STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
       __LINE__,__FILE__)
  pi=dacos(-1._real_8)
  ! ..initialise global rotation matrix
  CALL unitmx(rotmat,nstate)
  IF (wannr%w_ran.GT.0._real_8) THEN
     ! ..randomly distort initial vectors
     DO i=1,nstate
        gmat(i,i)=0._real_8
        DO j=i+1,nstate
           gmat(i,j)=repprngu()*wannr%w_ran
           gmat(j,i)=-gmat(i,j)
        ENDDO
     ENDDO
     CALL exponentiate(gmat,rotmat,nstate)
     CALL xyz_update(rotmat,xyzmat,aux,ldx,nstate)
     wannr%w_ran=0._real_8
  ENDIF
  DO iopt=1,wanni%w_maxs
     ! ..calculate the gradient
     CALL xgradat0(wanni%w_type,gmat,xyzmat,ldx,nstate)
     imax=idamax(nstate*nstate,gmat,1)
     gmax=ABS(dgive(gmat(1,1),imax))/parm%alat**2
     CALL ofunc(wanni%w_type,ofun,xyzmat,ldx,nstate)
     IF (gmax.LT.wannr%w_eps) GOTO 100
     scal=0.25_real_8*wannr%w_step/parm%alat**2
     CALL dscal(nstate*nstate,-scal,gmat,1)
     CALL exponentiate(gmat,rmat,nstate)
     CALL dcopy(nstate*nstate,rotmat,1,aux,1)
     CALL dgemm('N','N',nstate,nstate,nstate,1._real_8,aux,nstate,rmat,&
          nstate,0._real_8,rotmat,nstate)
     CALL xyz_update(rmat,xyzmat,aux,ldx,nstate)
  ENDDO
  IF (paral%io_parent)&
       WRITE(6,'(A,G12.3)') ' WANNIER CODE| NO CONVERGENCE (GMAX) ',gmax
100 CONTINUE
  DEALLOCATE(aux,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
       __LINE__,__FILE__)
  DEALLOCATE(gmat,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
       __LINE__,__FILE__)
  DEALLOCATE(rmat,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
       __LINE__,__FILE__)
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE sd_wannier
