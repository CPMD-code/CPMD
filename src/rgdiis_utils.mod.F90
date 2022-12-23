MODULE rgdiis_utils
  USE cotr,                            ONLY: cotc0,&
                                             dmax,&
                                             dtm,&
                                             hess
  USE error_handling,                  ONLY: stopgm
  USE fixcom_utils,                    ONLY: fixcom
  USE kinds,                           ONLY: real_8
  USE odiis_utils,                     ONLY: solve
  USE system,                          ONLY: cnti,&
                                             mxgdis
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rgdiis
  PUBLIC :: give_scr_rgdiis

CONTAINS

  ! ==================================================================
  SUBROUTINE rgdiis(xpar,dxpar,told,gold)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: xpar(*), dxpar(*), &
                                                told(cotc0%nodim,*), &
                                                gold(cotc0%nodim,*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'rgdiis'

    INTEGER :: i, id, iddc, iddg, idiis, idp1, idp2, idp3, idp4, idp5, idq1, &
      idq2, idq3, idq4, idq5, idtop, ierr, ii, info, iopt, isub, it, j, jj, &
      lw, nic, nowg, nrank
    INTEGER, SAVE                            :: icall = 0
    REAL(real_8)                             :: bc(mxgdis+1,mxgdis+1), eps, &
                                                qnrlen, sdlen, ss, tr, &
                                                vc(mxgdis+1)
    REAL(real_8), ALLOCATABLE                :: scr(:)
    REAL(real_8), EXTERNAL                   :: ddot

    CALL tiset(procedureN,isub)
    ! ==--------------------------------------------------------------==
    ! ==  UPDATE BUFFERS                                              ==
    ! ==--------------------------------------------------------------==
    icall = icall + 1
    idiis = icall
    IF (idiis.GT.cnti%mgdiis) THEN
       idiis=cnti%mgdiis
    ENDIF
    nowg = icall
    IF (nowg.GT.cnti%mgdiis) nowg=MOD(nowg-1,cnti%mgdiis)+1
    CALL dcopy(cotc0%nodim,xpar(1), 1,told(1,nowg),1)
    CALL dcopy(cotc0%nodim,dxpar(1),1,gold(1,nowg),1)
    ! ==--------------------------------------------------------------==
    ! ==  SCRATCH SPACE MANAGEMENT                                    ==
    ! ==--------------------------------------------------------------==
    ! TODO refactor this: create separate arrays instead of SCR
    idp1 = 1
    idp2 = idp1 + cotc0%nodim*cotc0%nodim
    idp3 = idp2 + cotc0%nodim*cnti%mgdiis
    idp4 = idp3 + cotc0%nodim
    idp5 = idp4 + 6*cotc0%nodim+idiis
    ! 
    idq1 = idp2
    idq2 = idq1 + cotc0%nodim
    idq3 = idq2 + cotc0%nodim
    idq4 = idq3 + cotc0%nodim
    idq5 = idq4 + cotc0%nodim*6
    idtop=MAX(idp5,idq5)
    ALLOCATE(scr(idtop),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! ==  CALCULATE ERROR VECTORS                                     ==
    ! ==--------------------------------------------------------------==
    CALL dcopy(cotc0%nodim*cotc0%nodim,hess(1,1),1,scr(idp1),1)
    it=idp2-1
    DO id=1,idiis
       DO i=1,cotc0%nodim
          it=it+1
          scr(it)=gold(i,id)
       ENDDO
    ENDDO
    iopt=12
    eps=1.e-6_real_8
    lw=6*cotc0%nodim+idiis
    CALL dgelss(cotc0%nodim,cotc0%nodim,idiis,scr(idp1),cotc0%nodim,scr(idp2),cotc0%nodim,&
         scr(idp3),eps,nrank,scr(idp4),lw,info)
    ! ==--------------------------------------------------------------==
    ! ==  cntl%diis MATRIX                                                 ==
    ! ==--------------------------------------------------------------==
    DO i=1,idiis
       ii=idp2+(i-1)*cotc0%nodim
       DO j=i,idiis
          jj=idp2+(j-1)*cotc0%nodim
          bc(i,j)=ddot(cotc0%nodim,scr(ii),1,scr(jj),1)
          bc(j,i)=bc(i,j)
       ENDDO
    ENDDO
    DO i=1,idiis+1
       vc(i)=0._real_8
       bc(i,idiis+1)=-1._real_8
       bc(idiis+1,i)=-1._real_8
    ENDDO
    vc(idiis+1)=-1._real_8
    bc(idiis+1,idiis+1)=0._real_8
    ! ==--------------------------------------------------------------==
    ! ==  Solve System of Linear Equations                            ==
    ! ==--------------------------------------------------------------==
    CALL solve(bc,mxgdis+1,idiis+1,vc)
    ! ==--------------------------------------------------------------==
    ! ==  Compute Interpolated Coefficient and Gradient Vectors       ==
    ! ==--------------------------------------------------------------==
    CALL zeroing(scr(idq1:idq1+cotc0%nodim-1))!,cotc0%nodim)
    CALL zeroing(scr(idq2:idq2+cotc0%nodim-1))!,cotc0%nodim)
    DO i=1,idiis
       iddc=idq1-1
       iddg=idq2-1
       DO j=1,cotc0%nodim
          iddc=iddc+1
          iddg=iddg+1
          scr(iddc)=scr(iddc)+vc(i)*told(j,i)
          scr(iddg)=scr(iddg)+vc(i)*gold(j,i)
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    ! ==  Estimate New Parameter Vectors                              ==
    ! ==--------------------------------------------------------------==
    CALL dcopy(cotc0%nodim*cotc0%nodim,hess(1,1),1,scr(idp1),1)
    iopt=12
    eps=1.e-6_real_8
    lw=6*cotc0%nodim
    CALL dgelss(cotc0%nodim,cotc0%nodim,1,scr(idp1),cotc0%nodim,scr(idq2),cotc0%nodim,&
         scr(idq3),eps,nrank,scr(idq4),lw,info)
    CALL daxpy(cotc0%nodim,-1.0_real_8,scr(idq2),1,scr(idq1),1)
    CALL daxpy(cotc0%nodim,-1.0_real_8,xpar(1),1,scr(idq1),1)
    ! ==--------------------------------------------------------------==
    ! ==  THIS IS A PRIMITIVE TRUST REGION CONCEPT                    ==
    ! ==--------------------------------------------------------------==
    nic=0
    DO i=1,cotc0%nodim
       scr(idq2+nic)=dtm(i)*dxpar(i)
       nic=nic+1
    ENDDO
    sdlen=SQRT(ddot(cotc0%nodim,scr(idq2),1,scr(idq2),1))
    qnrlen=SQRT(ddot(cotc0%nodim,scr(idq1),1,scr(idq1),1))
    tr=100.0_real_8 * sdlen
    IF (tr.LT.0.01_real_8) tr=0.01_real_8
    tr=MIN(dmax,tr)
    IF (qnrlen.GT.tr) THEN
       ss=SQRT(tr/qnrlen)
       CALL dscal(cotc0%nodim,ss,scr(idq1),1)
    ENDIF
    CALL fixcom(scr(idq1))
    CALL daxpy(cotc0%nodim,1.0_real_8,scr(idq1),1,xpar(1),1)

    DEALLOCATE(scr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rgdiis
  ! ==================================================================
  SUBROUTINE give_scr_rgdiis(lrgdiis,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lrgdiis
    CHARACTER(len=30)                        :: tag

! ==--------------------------------------------------------------==

    lrgdiis=MAX(&
         cotc0%nodim*cotc0%nodim&
         +cotc0%nodim*cnti%mgdiis&
         +7*cotc0%nodim+cnti%mgdiis,&
         10*cotc0%nodim)
    tag   ='NODIM*NODIM+NODIM*MGDIIS+...'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_rgdiis
  ! ==================================================================

END MODULE rgdiis_utils
