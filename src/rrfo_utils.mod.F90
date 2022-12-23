MODULE rrfo_utils
  USE cotr,                            ONLY: anorm,&
                                             cotc0,&
                                             dmax,&
                                             dtm,&
                                             fc,&
                                             hess,&
                                             xlagr
  USE error_handling,                  ONLY: stopgm
  USE fixcom_utils,                    ONLY: fixcom
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: paral
  USE system,                          ONLY: cnti
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: dspevy
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rrfo
  PUBLIC :: give_scr_rrfo
  !public :: fitsol

CONTAINS

  ! ==================================================================
  SUBROUTINE rrfo(xpar,dxpar)
    ! ==--------------------------------------------------------------==
    ! ==  RATIONAL FUNCTION OPTIMIZATON                               ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: xpar(cotc0%nodim), &
                                                dxpar(cotc0%nodim)

    CHARACTER(*), PARAMETER                  :: procedureN = 'rrfo'

    INTEGER :: i, ierr, iopt, ipahess, ipaux, ipdlam, ipeigv, ipftran, &
      iphstep, ipumat, isub, j, k, ndhess, neneg, nepos, nezer
    REAL(real_8)                             :: dnrm2, qnrlen, sdlen, ss, tr, &
                                                zneg, zpos
    REAL(real_8), ALLOCATABLE                :: ahess(:), aux(:), dlam(:), &
                                                eigv(:), ftran(:), hstep(:), &
                                                umat(:)

    CALL tiset('      RRFO',isub)
    ! ==--------------------------------------------------------------==
    ! ==  SCRATCH SPACE ALLOCATION                                    ==
    ! ==                                                              ==
    ! ==  IPA1  : AUGMENTED HESSIAN                                   ==
    ! ==  IPA2  : EIGV                                                ==
    ! ==  IPA3  : UMAT                                                ==
    ! ==  IPA4  : AUX                                                 ==
    ! ==  IPB1  : DLAM                                                ==
    ! ==  IPB2  : FTRAN                                               ==
    ! ==  IPB3  : HSTEP                                               ==
    ! ==--------------------------------------------------------------==
    ndhess = cotc0%nodim + cotc0%mcnstr
    ipahess = 1
    ipeigv  = ipahess + ndhess*ndhess
    ipumat  = ipeigv  + ndhess
    ipaux   = ipumat  + ndhess*ndhess
    ipdlam  = ipaux   + 3*ndhess
    ipftran = ipdlam  + ndhess
    iphstep = ipftran + ndhess
    ALLOCATE(ahess((cotc0%nodim + cotc0%mcnstr)*(cotc0%nodim + cotc0%mcnstr)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(umat((cotc0%nodim + cotc0%mcnstr)*(cotc0%nodim + cotc0%mcnstr)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(eigv(cotc0%nodim + cotc0%mcnstr),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(aux(3*(cotc0%nodim + cotc0%mcnstr)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(dlam(cotc0%nodim + cotc0%mcnstr),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(ftran(cotc0%nodim + cotc0%mcnstr),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(hstep(cotc0%nodim + cotc0%mcnstr),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! ==  STEEPEST DESCENT STEP                                       ==
    ! ==--------------------------------------------------------------==
    !$omp parallel do private(I)
    DO i=1,cotc0%nodim
       aux(i)=dtm(i)*dxpar(i)
    ENDDO
    sdlen=dnrm2(cotc0%nodim,aux(1),1)
    ! ==--------------------------------------------------------------==
    ! ==  AUGMENTED HESSIAN IN UPPER TRAING. STORAGE MODE             ==
    ! ==--------------------------------------------------------------==
    CALL zeroing(ahess)!,ndhess*ndhess)
    k=0
    DO j=1,ndhess
       DO i=1,j
          k=k+1
          IF (i.LE.cotc0%nodim) THEN
             IF (j.LE.cotc0%nodim) THEN
                ahess(k)=hess(i,j)
             ELSE
                ahess(k)=-anorm(i,j-cotc0%nodim)
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    iopt=21
    CALL dspevy(iopt,ahess,eigv,umat,ndhess,ndhess,aux,3*ndhess)
    ! ==--------------------------------------------------------------==
    ! ==  DETERMINE EIGENVALUE CHARACTERISTICS                        ==
    ! ==--------------------------------------------------------------==
    neneg=0
    nezer=0
    nepos=0
    DO i=1,ndhess
       IF (eigv(i).LT.-1.e-2_real_8) THEN
          neneg=neneg+1
       ELSEIF (eigv(i).GT.1.e-2_real_8) THEN
          nepos=nepos+1
       ELSE
          nezer=nezer+1
       ENDIF
    ENDDO
    IF (paral%io_parent)&
         WRITE(6,'(A,I4,A,I4,A,I4)') ' RFO| MODES: NEGATIVE=',neneg,&
         '     ZERO=',nezer,'     POSITIVE=',nepos
    ! ..INFORMATION ABOUT THE LOWEST HESSIAN EIGENVALUES
    IF (cnti%nsorder.GT.0) THEN
       IF (paral%io_parent)&
            WRITE(6,'(A)') ' LOWEST HESSIAN EIGENVALUES:'
       IF (paral%io_parent)&
            WRITE(6,'(4F9.6)') eigv(1),eigv(2),eigv(3),eigv(4),eigv(5),&
            eigv(6),eigv(7),eigv(8)
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL dcopy(cotc0%nodim,dxpar(1),1,dlam(1),1)
    IF (cotc0%mcnstr.GT.0) THEN
       CALL dcopy(cotc0%mcnstr,fc(1),1,dlam(cotc0%nodim+1),1)
       CALL dscal(cotc0%mcnstr,-1.0_real_8,dlam(cotc0%nodim+1),1)
    ENDIF
    CALL dgemv('T',ndhess,ndhess,1.0_real_8,umat,ndhess,dlam,1,&
         0.0_real_8,ftran,1)
    ! ==--------------------------------------------------------------==
    neneg=cotc0%mcnstr+cnti%nsorder
    zneg=0.0_real_8
    IF (neneg.GT.0) zneg=fitsol(neneg,ftran(1),eigv(1))
    zpos=fitsol(nepos,ftran(neneg+nezer+1),eigv(neneg+nezer+1))
    IF (paral%io_parent)&
         WRITE(6,'(A,2F12.6)') ' RFO| MODE FOLLOWING PARAMETERS :',&
         zneg,zpos
    ! mb      DO I=1,NDHESS
    ! mb        IF(I.LE.NENEG) THEN
    ! mb          FTRAN(I)=FTRAN(I)/(EIGV(I)-ZNEG)
    ! mb        ELSE
    ! mb          FTRAN(I)=FTRAN(I)/(EIGV(I)-ZPOS)
    ! mb        ENDIF
    ! mb      ENDDO
    ! mb-> separate the two loops
    DO i=1,neneg         ! negative
       ftran(i)=ftran(i)/(eigv(i)-zneg)
    ENDDO
    DO i=neneg+1,ndhess  ! positive
       ftran(i)=ftran(i)/(eigv(i)-zpos)
    ENDDO
    ! ==--------------------------------------------------------------==
    CALL dgemv('N',ndhess,ndhess,-1.0_real_8,umat,ndhess,ftran,1,&
         0.0_real_8,hstep(1),1)
    qnrlen=dnrm2(cotc0%nodim,hstep(1),1)
    ! ==--------------------------------------------------------------==
    ! ==  THIS IS A PRIMITIVE TRUST REGION CONCEPT                    ==
    ! ==--------------------------------------------------------------==
    tr=100.0_real_8 * sdlen
    IF (tr.LT.0.01_real_8) tr=0.01_real_8
    tr=MIN(tr,dmax)
    IF (qnrlen.GT.tr) THEN
       ss=SQRT(tr/qnrlen)
       CALL dscal(ndhess,ss,hstep(1),1)
    ENDIF
    CALL fixcom(hstep)
    CALL daxpy(cotc0%nodim,1.0_real_8,hstep(1),1,xpar(1),1)
    IF (cotc0%mcnstr.GT.0)&
         CALL daxpy(cotc0%mcnstr,1.0_real_8,hstep(cotc0%nodim+1),1,xlagr(1),1)
    CALL tihalt('      RRFO',isub)
    ! ==--------------------------------------------------------------==
    DEALLOCATE(ahess,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(umat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(eigv,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(aux,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(dlam,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(ftran,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(hstep,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rrfo
  ! ==================================================================
  FUNCTION fitsol(n,f,e)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: n
    REAL(real_8)                             :: f(n), e(n), fitsol

    INTEGER, PARAMETER                       :: maxit = 100 
    REAL(real_8), PARAMETER                  :: conv = 1.e-10_real_8 

    INTEGER                                  :: i, j
    REAL(real_8)                             :: df, fnew, fold

! ==--------------------------------------------------------------==

    fold=0.0_real_8
    DO i=1,maxit
       fnew=0.0_real_8
       DO j=1,n
          fnew=fnew+f(j)*f(j)/(fold-e(j))
       ENDDO
       df=ABS(fold-fnew)
       IF (df.LT.conv) GOTO 100
       fold=fnew
    ENDDO
    IF (paral%io_parent)&
         WRITE(6,*) ' FITSOL| ',fnew,df
100 CONTINUE
    fitsol=fnew
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION fitsol
  ! ==================================================================
  SUBROUTINE give_scr_rrfo(lrrfo,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lrrfo
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: ndhess

    ndhess = cotc0%nodim + cotc0%mcnstr
    lrrfo=2*ndhess*ndhess+7*ndhess
    tag ='2*NDHESS*NDHESS+7*NDHESS'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_rrfo
  ! ==================================================================

END MODULE rrfo_utils
