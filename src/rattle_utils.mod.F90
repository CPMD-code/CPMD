MODULE rattle_utils
  USE cnstfc_utils,                    ONLY: cnstfc
  USE cotr,                            ONLY: anorm,&
                                             cotc0,&
                                             dtm,&
                                             fc,&
                                             kshove,&
                                             ylagr
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE mm_dim_utils,                    ONLY: mm_dim
  USE mm_dimmod,                       ONLY: mm_go_mm,&
                                             mm_revert
  USE parac,                           ONLY: paral
  USE puttau_utils,                    ONLY: gettau,&
                                             puttau
  USE system,                          ONLY: maxsys
  USE tpar,                            ONLY: dt_ions
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rattle

CONTAINS

  ! ==================================================================
  SUBROUTINE rattle(tau0,velp)
    ! ==--------------------------------------------------------------==
    ! ==  GENERAL CONSTRAINTS ON VELOCITIES FOR VELOCITY VERLET       ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:), velp(:,:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'rattle'
    REAL(real_8), PARAMETER                  :: tols = 1.e-12_real_8, &
                                                tolx = 1.e-8_real_8 

    INTEGER                                  :: i, ierr, info, j, k, naux, &
                                                nrank
    INTEGER, SAVE                            :: istart = 0
    LOGICAL                                  :: oldstatus
    REAL(real_8)                             :: cnmax, errx, xval
    REAL(real_8), ALLOCATABLE, SAVE          :: asl(:,:), aux(:), dv(:), &
                                                dx(:), tscr(:,:,:)
    REAL(real_8), EXTERNAL                   :: dasum

    IF (cotc0%mcnstr.EQ.0) RETURN
    IF (cotc0%lshove .AND. istart.NE.0) THEN
       xval=dasum(cotc0%mcnstr,kshove,1)
       IF (xval.LT.0.5_real_8) RETURN
    ENDIF
    ! moved out of the if to keep in initialized (Joost<-that is a bug?) ...
    naux   =6*cotc0%nodim+cotc0%mcnstr
    IF (istart.EQ.0) THEN
       ! ISTART =1
       ! NAUX   =6*NODIM+MCNSTR

       ! AK: we need to make sure that TSCR is large enough for QM/MM
       CALL mm_dim(mm_go_mm,oldstatus)
       ALLOCATE(tscr(3,maxsys%nax,maxsys%nsx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(dx(cotc0%nodim),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(dv(cotc0%nodim),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(asl(cotc0%mcnstr,cotc0%mcnstr),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(aux(naux),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL mm_dim(mm_revert,oldstatus)
    ENDIF
    CALL zeroing(dx)!,cotc0%nodim)
    CALL cnstfc(tau0,tscr,dx,cnmax,.TRUE.)
    CALL puttau(velp,dx)
    CALL zeroing(asl)!,cotc0%mcnstr*cotc0%mcnstr)
    DO i=1,cotc0%mcnstr
       IF (cotc0%lshove.AND. istart.NE.0 ) THEN
          IF (kshove(i).LT.0.5_real_8) CALL zeroing(anorm(:,i))!,cotc0%nodim)
       ENDIF
       dv(i)=0._real_8
       DO j=1,cotc0%nodim
          dv(i)=dv(i)+anorm(j,i)*dx(j)
       ENDDO
       DO j=1,cotc0%mcnstr
          DO k=1,cotc0%nodim
             asl(i,j)=asl(i,j)+anorm(k,i)*anorm(k,j)*dtm(k)/dt_ions/2._real_8
          ENDDO
       ENDDO
    ENDDO
    IF (istart.EQ.0) istart =1
    CALL dgelss(cotc0%mcnstr,cotc0%mcnstr,1,asl,cotc0%mcnstr,dv,cotc0%mcnstr,fc,1.e-12_real_8,nrank,&
         aux,naux,info)
    IF (info.NE.0) CALL stopgm('RATTLE','DGELSS ENDED WITH INFO.NE.0',& 
         __LINE__,__FILE__)
    DO i=1,cotc0%mcnstr
       ylagr(i)=dv(i)
    ENDDO
    DO j=1,cotc0%nodim
       dv(j)=0._real_8
       DO i=1,cotc0%mcnstr
          dv(j)=dv(j)-ylagr(i)*anorm(j,i)*dtm(j)/dt_ions/2.0_real_8
       ENDDO
    ENDDO
    DO j=1,cotc0%nodim
       dx(j)=dx(j)+dv(j)
    ENDDO
    CALL gettau(velp,dx)
    ! Check accuracy
    DO i=1,cotc0%mcnstr
       fc(i)=0._real_8
       DO j=1,cotc0%nodim
          fc(i)=fc(i)+anorm(j,i)*dx(j)
       ENDDO
    ENDDO
    errx=dasum(cotc0%mcnstr,fc(1),1)
    IF (errx.GT.tolx) THEN
       IF (paral%io_parent)&
            WRITE(6,*) ' RATTLE| Constraint can not be fulfilled'
       CALL stopgm('RATTLE',' ',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rattle
  ! ==================================================================

END MODULE rattle_utils
