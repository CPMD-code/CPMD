MODULE hessup_utils
  USE cotr,                            ONLY: cotc0,&
                                             hess,&
                                             sdpl
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: hessup
  PUBLIC :: hespow

CONTAINS

  ! ==================================================================
  SUBROUTINE hessup(xpar,ypar,dxpar,dypar)
    ! ==--------------------------------------------------------------==
    ! ==  cntl%bfgs update of the nuclear hessian                          ==
    ! ==  cntl%bfgs: Broyden-Fletcher-Goldfarb-Shanno                      ==
    ! ==  See Numerical recipes p308                                  ==
    ! ==  (Variable Metric Methods in Multidimensions)                ==
    ! ==--------------------------------------------------------------==
    ! ==  IN OUR CASE WE COMPUTE THE HESSIAN BY DEFAULT               ==
    ! ==  (if TINVBFGS=.TRUE. use inverse Hessian matrix)             ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: xpar(*), ypar(*), dxpar(*), &
                                                dypar(*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'hessup'
    REAL(real_8), PARAMETER                  :: eps = 1.e-10_real_8 

    INTEGER                                  :: i, ierr, isub
    REAL(real_8)                             :: shs, ys
    REAL(real_8), ALLOCATABLE                :: hs(:), s(:), y(:)
    REAL(real_8), EXTERNAL                   :: ddot

    CALL tiset('    HESSUP',isub)
    ! ==--------------------------------------------------------------==
    ! Allocation of SCR
    ALLOCATE(s(cotc0%nodim),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(y(cotc0%nodim),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(hs(cotc0%nodim),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    DO i=1,cotc0%nodim
       s(i)=xpar(i)-ypar(i)
       y(i)=dxpar(i)-dypar(i)
    ENDDO
    ys=ddot(cotc0%nodim,y(1),1,s(1),1)
    IF (sdpl%tinvbfgs) THEN
       ! cntl%bfgs using inverse Hessian matrix
       CALL dgemv('N',cotc0%nodim,cotc0%nodim,1.0_real_8,hess,cotc0%nodim,y,1,0.0_real_8,hs,1)
       shs=ddot(cotc0%nodim,y(1),1,hs(1),1)
       ! Davidon_fletcher-Powell (DFP) update
       IF (ABS(ys).GT.eps .AND. ABS(shs).GT.eps) THEN
          CALL dger(cotc0%nodim,cotc0%nodim,1.0_real_8/ys,s,1,s,1,hess,cotc0%nodim)
          CALL dger(cotc0%nodim,cotc0%nodim,-1.0_real_8/shs,hs,1,hs,1,hess,cotc0%nodim)
          ! Additional term of cntl%bfgs
          CALL dscal(cotc0%nodim,1.0_real_8/ys,s,1)
          CALL daxpy(cotc0%nodim,-1.0_real_8/shs,hs,1,s,1)
          CALL dger(cotc0%nodim,cotc0%nodim,shs,s,1,s,1,hess,cotc0%nodim)
       ENDIF
    ELSE
       ! cntl%bfgs using Hessian matrix
       CALL dgemv('N',cotc0%nodim,cotc0%nodim,1.0_real_8,hess,cotc0%nodim,s,1,0.0_real_8,hs,1)
       shs=ddot(cotc0%nodim,s(1),1,hs(1),1)
       IF (ABS(ys).GT.eps .AND. ABS(shs).GT.eps) THEN
          CALL dger(cotc0%nodim,cotc0%nodim,1.0_real_8/ys,y,1,y,1,hess,cotc0%nodim)
          CALL dger(cotc0%nodim,cotc0%nodim,-1.0_real_8/shs,hs,1,hs,1,hess,cotc0%nodim)
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    DEALLOCATE(s,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(y,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(hs,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL tihalt('    HESSUP',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE hessup
  ! ==================================================================
  SUBROUTINE hespow(xpar,ypar,dxpar,dypar)
    ! ==--------------------------------------------------------------==
    ! ==  POWELL update of the nuclear hessian (for transition state  ==
    ! ==  search)                                                     ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: xpar(cotc0%nodim), &
                                                ypar(cotc0%nodim), &
                                                dxpar(cotc0%nodim), &
                                                dypar(cotc0%nodim)

    CHARACTER(*), PARAMETER                  :: procedureN = 'hespow'

    INTEGER                                  :: i, ierr, isub
    REAL(real_8)                             :: ss, ys
    REAL(real_8), ALLOCATABLE                :: hs(:), s(:), y(:)
    REAL(real_8), EXTERNAL                   :: ddot

    CALL tiset('    HESPOW',isub)
    ! ==--------------------------------------------------------------==
    ! Allocation of SCR
    ALLOCATE(s(cotc0%nodim),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(y(cotc0%nodim),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(hs(cotc0%nodim),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    DO i=1,cotc0%nodim
       s(i)=xpar(i)-ypar(i)
       y(i)=dxpar(i)-dypar(i)
    ENDDO
    ss=ddot(cotc0%nodim,s,1,s,1)
    CALL dgemv('N',cotc0%nodim,cotc0%nodim,-1.0_real_8,hess,cotc0%nodim,s,1,1.0_real_8,y,1)
    ys=ddot(cotc0%nodim,y,1,s,1)
    ys=ys/ss
    CALL dger(cotc0%nodim,cotc0%nodim,1.0_real_8/ss,y,1,s,1,hess,cotc0%nodim)
    CALL dger(cotc0%nodim,cotc0%nodim,1.0_real_8/ss,s,1,y,1,hess,cotc0%nodim)
    CALL dger(cotc0%nodim,cotc0%nodim,-ys/ss,s,1,s,1,hess,cotc0%nodim)
    ! ==--------------------------------------------------------------==
    DEALLOCATE(s,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(y,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(hs,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL tihalt('    HESPOW',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE hespow
  ! ==================================================================

END MODULE hessup_utils
