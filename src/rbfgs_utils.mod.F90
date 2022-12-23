MODULE rbfgs_utils
  USE cotr,                            ONLY: cotc0,&
                                             dmax,&
                                             dtm,&
                                             hess,&
                                             sdpl
  USE error_handling,                  ONLY: stopgm
  USE fixcom_utils,                    ONLY: fixcom
  USE kinds,                           ONLY: real_8
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rbfgs
  PUBLIC :: give_scr_rbfgs

CONTAINS

  ! ==================================================================
  SUBROUTINE rbfgs(xpar,dxpar)
    ! ==--------------------------------------------------------------==
    ! ==  SOLVE THE NEWTON EQUATION AND DETERMINE AN IONIC STEP       ==
    ! ==--------------------------------------------------------------==
    ! == INPUT:                                                       ==
    ! ==   XPAR(NODIM)                                                ==
    ! ==  DXPAR(NODIM)                                                ==
    ! == OUTPUT:                                                      ==
    ! ==   XPAR(NODIM)                                                ==
    ! ==  (NODIM Number of Freedom degres                             ==
    ! == Use HESS(NODIM,NODIM) Hessian matrix (estimated)             ==
    ! ==     DTM(NODIM)        Step                                   == 
    ! ==--------------------------------------------------------------==
    ! == IF TINVBFGS=.TRUE. use inverse Hessian matrix                ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: xpar(cotc0%nodim), &
                                                dxpar(cotc0%nodim)

    CHARACTER(*), PARAMETER                  :: procedureN = 'rbfgs'
    REAL(real_8), PARAMETER                  :: eps = 1.e-6_real_8

    INTEGER                                  :: i, ierr, info, isub, nrank
    REAL(real_8)                             :: qnrlen, sdlen, ss, tr
    REAL(real_8), ALLOCATABLE                :: aux(:), hessvd(:,:), &
                                                sfion(:), sval(:)
    REAL(real_8), EXTERNAL                   :: ddot

    CALL tiset('     RBFGS',isub)
    ! ==--------------------------------------------------------------==
    ! ==  SCRATCH SPACE ALLOCATION:                                   ==
    ! ==    HESSVD  HESS IN SVD                                       ==
    ! ==    SFION   FION                                              ==
    ! ==    SVAL    SINGULAR VALUES                                   ==
    ! ==--------------------------------------------------------------==
    IF (sdpl%tinvbfgs) THEN
       ! use inverse Hessian Matrix
       CALL dgemv('N',cotc0%nodim,cotc0%nodim,-1._real_8,hess,cotc0%nodim,dxpar,1,1._real_8,xpar,1)
       RETURN
    ENDIF
    ! ==--------------------------------------------------------------==
    ! TODO align for BG
    ALLOCATE(aux(6*cotc0%nodim),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(hessvd(cotc0%nodim, cotc0%nodim),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(sfion(cotc0%nodim),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(sval(cotc0%nodim),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! ==  SOLVE HESS*S = -FION                                        ==
    ! ==  (FION is given, the unknown variable is S)                  == 
    ! ==--------------------------------------------------------------==
    CALL dcopy(cotc0%nodim*cotc0%nodim,hess(1,1),1,hessvd(1,1),1)
    DO i=1,cotc0%nodim
       sfion(i)=-dxpar(i)
       aux(i)=dtm(i)*dxpar(i)
    ENDDO
    sdlen=SQRT(ddot(cotc0%nodim,aux,1,aux,1))
    CALL dgelss(cotc0%nodim,cotc0%nodim,1,hessvd,cotc0%nodim,sfion,cotc0%nodim,sval,eps,nrank,&
         aux,6*cotc0%nodim,info)
    qnrlen=SQRT(ddot(cotc0%nodim,sfion,1,sfion,1))
    ! ==--------------------------------------------------------------==
    ! ==  THIS IS A PRIMITIVE TRUST REGION CONCEPT                    ==
    ! ==--------------------------------------------------------------==
    tr=100.0_real_8 * sdlen
    IF (tr.LT.0.01_real_8) tr=0.01_real_8
    tr=MIN(tr,dmax)
    IF (qnrlen.GT.tr) THEN
       ss=SQRT(tr/qnrlen)
       CALL dscal(cotc0%nodim,ss,sfion,1)
    ENDIF
    CALL fixcom(sfion)
    CALL daxpy(cotc0%nodim,1.0_real_8,sfion,1,xpar(1),1)
    ! ==--------------------------------------------------------------==
    DEALLOCATE(aux,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(hessvd,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(sfion,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(sval,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL tihalt('     RBFGS',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rbfgs
  ! ==================================================================
  SUBROUTINE give_scr_rbfgs(lrbfgs,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lrbfgs
    CHARACTER(len=30)                        :: tag

! ==--------------------------------------------------------------==

    IF (sdpl%tinvbfgs) THEN
       lrbfgs=0
       tag  ='0'
    ELSE
       lrbfgs=cotc0%nodim*cotc0%nodim+8*cotc0%nodim+10
       tag  ='NODIM*NODIM+8*NODIM'
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_rbfgs
  ! ==================================================================

END MODULE rbfgs_utils
