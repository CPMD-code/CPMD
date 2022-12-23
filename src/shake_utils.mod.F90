MODULE shake_utils
  USE cnstfc_utils,                    ONLY: cnstfc
  USE cotr,                            ONLY: &
       anorm, cnsval, cnsval_dest, cotc0, dtm, fc, fv, grate, kshove, &
       mm_askel, ntcnst, xlagr, ylagr
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE mm_dim_utils,                    ONLY: mm_dim
  USE mm_dimmod,                       ONLY: mm_go_mm,&
                                             mm_revert
  USE odiis_utils,                     ONLY: solve
  USE parac,                           ONLY: paral
  USE puttau_utils,                    ONLY: gettau
  USE system,                          ONLY: cntl,&
                                             maxsys
  USE tpar,                            ONLY: dt_ions,&
                                             dtb2mi
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: cpmdshake
  !public :: solvs
  !public :: andforces

CONTAINS

  ! ==================================================================
  SUBROUTINE cpmdshake(tau0,taup,velp)
    ! ==--------------------------------------------------------------==
    ! ==  GENERAL CONSTRAINTS ON POSITIONS FOR VELOCITY VERLET        ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:), taup(:,:,:), &
                                                velp(:,:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'cpmdshake'
    INTEGER, PARAMETER                       :: maxrat = 5000, mrdiis = 5 
    REAL(real_8), PARAMETER                  :: tolf = 1.e-7_real_8, &
                                                tolx = 1.e-8_real_8 

    INTEGER                                  :: i, i_dim, ia, idiis, ierr, &
                                                info, is, iter, j, k, m
    INTEGER, ALLOCATABLE, SAVE               :: ipvt(:)
    INTEGER, SAVE                            :: istart = 0
    LOGICAL                                  :: oldstatus
    REAL(real_8)                             :: cnmax, &
                                                diism(mrdiis+1,mrdiis+1), &
                                                errf, errx, fact, &
                                                v(mrdiis+1), xval
    REAL(real_8), ALLOCATABLE, SAVE          :: an0(:,:), asl(:,:), dx(:), &
                                                err(:,:), tscr(:,:,:), &
                                                xlo(:,:)
    REAL(real_8), EXTERNAL                   :: dasum, ddot

! change constraint value according to growth rate

    DO i=1,cotc0%mcnstr
       cnsval(i)=cnsval(i)+grate(i)*dt_ions
       IF (cnsval_dest(i)/=-999._real_8) THEN
          IF (grate(i)>0._real_8.AND.cnsval(i)>cnsval_dest(i)) THEN
             cnsval(i)=cnsval_dest(i) ! increase
          ELSEIF (grate(i)<0._real_8.AND.cnsval(i)<cnsval_dest(i)) THEN
             cnsval(i)=cnsval_dest(i) ! decrease
          ENDIF
       ENDIF
    ENDDO
    ! 
    CALL dcopy(cotc0%mcnstr,ylagr(1),1,xlagr(1),1)
    IF (istart.EQ.0) THEN
       istart=1
       ! AK: we need to make sure that TSCR is large enough for QM/MM
       CALL mm_dim(mm_go_mm,oldstatus)
       ALLOCATE(tscr(3,maxsys%nax,maxsys%nsx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(dx(3*cotc0%nodim),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(asl(cotc0%mcnstr,cotc0%mcnstr),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(ipvt(cotc0%mcnstr),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(xlo(cotc0%mcnstr,mrdiis),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(err(cotc0%mcnstr,mrdiis),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(an0(cotc0%nodim,cotc0%mcnstr),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(xlagr)!,cotc0%mcnstr)
       IF (cotc0%lshove)  THEN
          ALLOCATE(kshove(cotc0%mcnstr),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
       CALL mm_dim(mm_revert,oldstatus)
    ENDIF
    CALL zeroing(dx)!,cotc0%nodim)
    CALL cnstfc(tau0,tscr,dx,cnmax,.TRUE.)
    CALL dcopy(cotc0%mcnstr*cotc0%nodim,anorm(1,1),1,an0(1,1),1)
    xval=cotc0%mcnstr
    IF (cotc0%lshove) THEN
       DO i=1,cotc0%mcnstr
          kshove(i)=1._real_8
          xval=(fv(i)-cnsval(i))*ntcnst(6,i)
          IF (xval.GT.1.e-14_real_8) THEN
             kshove(i)=0._real_8
             CALL zeroing(anorm(:,i))!,cotc0%nodim)
             CALL zeroing(an0(:,i))!,cotc0%nodim)
             cnsval(i)=fv(i)
             xlagr(i)=0._real_8
          ENDIF
       ENDDO
       xval=dasum(cotc0%mcnstr,kshove,1)
       IF (xval.LT.0.5_real_8) RETURN
    ENDIF
    ! First guess for the force
    CALL zeroing(dx)!,cotc0%nodim)
    DO i=1,cotc0%mcnstr
       DO j=1,cotc0%nodim
          dx(j)=dx(j)+xlagr(i)*an0(j,i)
       ENDDO
    ENDDO
    ! Update the positions
    CALL zeroing(tscr)!,3*maxsys%nax*maxsys%nsx)
    CALL gettau(tscr,dx)
    DO is=1,ions1%nsp
       fact=-dt_ions*dtb2mi(is)
       DO ia=1,ions0%na(is)
          tau0(1,ia,is)=taup(1,ia,is)+fact*tscr(1,ia,is)
          tau0(2,ia,is)=taup(2,ia,is)+fact*tscr(2,ia,is)
          tau0(3,ia,is)=taup(3,ia,is)+fact*tscr(3,ia,is)
       ENDDO
    ENDDO
    ! Iterativ calculation of lambda
    DO iter=1,maxrat
       ! Calculate constraint function and forces
       ! for the current value of lambda
       CALL zeroing(dx)!,cotc0%nodim)
       CALL cnstfc(tau0,tscr,dx,cnmax,.TRUE.)
       CALL dscal(cotc0%mcnstr,-1.0_real_8,fc(1),1)
       errf=dasum(cotc0%mcnstr,fc(1),1)
       IF (errf.LT.tolf) GOTO 100
       ! Derivatives of sigma wrt lambda
       IF (cntl%tqmmm) THEN
          DO i=1,cotc0%mcnstr
             DO j=1,cotc0%mcnstr
                asl(i,j)=0._real_8
                DO i_dim=1,12
                   k=mm_ASKEL(i,i_dim)
                   IF (k.NE.0)&
                        asl(i,j)=asl(i,j)-dt_ions*dtm(k)*anorm(k,i)*an0(k,j)
                ENDDO
             ENDDO
          ENDDO
       ELSE
          DO i=1,cotc0%mcnstr
             DO j=1,cotc0%mcnstr
                asl(i,j)=0._real_8
                DO k=1,cotc0%nodim
                   asl(i,j)=asl(i,j)-dt_ions*dtm(k)*anorm(k,i)*an0(k,j)
                ENDDO
             ENDDO
          ENDDO
       ENDIF

       ! Solve for asl*cg=fc
       ! LAPACK matrix solver
       IF (NINT(xval).EQ.cotc0%mcnstr) THEN
          info=0
          CALL dgesv(cotc0%mcnstr,1,asl,cotc0%mcnstr,ipvt,fc,cotc0%mcnstr,info)
          IF (info.NE.0) CALL stopgm(' CPMDSHAKE|','ERROR IN DGESV',& 
               __LINE__,__FILE__)
       ELSE
          CALL solvs(cotc0%mcnstr,asl,fc,ipvt)
       ENDIF
       errx=dasum(cotc0%mcnstr,fc(1),1)
       ! cntl%diis!
       idiis=MOD(iter-1,mrdiis)+1
       CALL dcopy(cotc0%mcnstr,xlagr(1),1,xlo(1,idiis),1)
       CALL dcopy(cotc0%mcnstr,fc(1),1,err(1,idiis),1)
       IF (iter.GT.mrdiis) THEN
          m=mrdiis+1
          CALL zeroing(diism)!,m*m)
          CALL zeroing(v)!,m)
          DO i=1,mrdiis
             DO j=1,mrdiis
                diism(i,j)=ddot(cotc0%mcnstr,err(1,i),1,err(1,j),1)
             ENDDO
             diism(m,i)=1._real_8
             diism(i,m)=1._real_8
          ENDDO
          v(m)=1._real_8
          CALL solve(diism,m,m,v)
          CALL zeroing(fc)!,cotc0%mcnstr)
          CALL zeroing(xlagr)!,cotc0%mcnstr)
          DO i=1,mrdiis
             DO j=1,cotc0%mcnstr
                fc(j)=fc(j)+v(i)*err(j,i)
                xlagr(j)=xlagr(j)+v(i)*xlo(j,i)
             ENDDO
          ENDDO
       ENDIF
       CALL daxpy(cotc0%mcnstr,1.0_real_8,fc(1),1,xlagr(1),1)
       IF (errx.LT.tolx) GOTO 100
       ! Update forces
       CALL zeroing(dx)!,cotc0%nodim)
       DO i=1,cotc0%mcnstr
          DO j=1,cotc0%nodim
             dx(j)=dx(j)+xlagr(i)*an0(j,i)
          ENDDO
       ENDDO
       ! Update the positions
       CALL zeroing(tscr)!,3*maxsys%nax*maxsys%nsx)
       CALL gettau(tscr,dx)
       DO is=1,ions1%nsp
          fact=-dt_ions*dtb2mi(is)
          DO ia=1,ions0%na(is)
             tau0(1,ia,is)=taup(1,ia,is)+fact*tscr(1,ia,is)
             tau0(2,ia,is)=taup(2,ia,is)+fact*tscr(2,ia,is)
             tau0(3,ia,is)=taup(3,ia,is)+fact*tscr(3,ia,is)
          ENDDO
       ENDDO
    ENDDO
    IF (paral%io_parent)&
         WRITE(6,'(2A,I6,A,F10.8,A,F10.8)') ' CPMDSHAKE| DID NOT ',&
         'CONVERGE AFTER ',iter,' ITERATIONS. ERRX=',&
         errx, ' TOLX=',tolx
    CALL stopgm('CPMDSHAKE',' ',& 
         __LINE__,__FILE__)
100 CONTINUE
    CALL dcopy(3*maxsys%nax*maxsys%nsx,tau0(1,1,1),1,taup(1,1,1),1)
    ! Update velocities
    CALL zeroing(dx)!,cotc0%nodim)
    DO i=1,cotc0%mcnstr
       DO j=1,cotc0%nodim
          dx(j)=dx(j)+xlagr(i)*an0(j,i)
       ENDDO
    ENDDO
    CALL zeroing(tscr)!,3*maxsys%nax*maxsys%nsx)
    CALL gettau(tscr,dx)
    ! FIXME: add OpenMP?? AK
    DO is=1,ions1%nsp
       fact=-dtb2mi(is)
       DO ia=1,ions0%na(is)
          velp(1,ia,is)=velp(1,ia,is)+fact*tscr(1,ia,is)
          velp(2,ia,is)=velp(2,ia,is)+fact*tscr(2,ia,is)
          velp(3,ia,is)=velp(3,ia,is)+fact*tscr(3,ia,is)
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE cpmdshake
  ! ==================================================================
  SUBROUTINE solvs(m,asl,fc,ipvt)
    INTEGER                                  :: m
    REAL(real_8)                             :: asl(m,m), fc(m)
    INTEGER                                  :: ipvt(*)

    INTEGER, PARAMETER                       :: mx = 100 

    INTEGER                                  :: i, il(mx), info, n

! ==--------------------------------------------------------------==

    IF (m.GT.mx) CALL stopgm('SOLVS','MX',& 
         __LINE__,__FILE__)
    n=0
    DO i=1,m
       IF (ABS(asl(i,i)).GT.1.e-12_real_8) THEN
          n=n+1
          il(i)=n
          IF (i.NE.n) THEN
             CALL dcopy(m,asl(1,i),1,asl(1,n),1)
             CALL dcopy(m,asl(i,1),m,asl(n,1),m)
             fc(n)=fc(i)
          ENDIF
       ELSE
          il(i)=0
       ENDIF
    ENDDO
    info=0
    CALL dgesv(n,1,asl,m,ipvt,fc,m,info)
    IF (info.NE.0) CALL stopgm(' SOLVS|','ERROR IN DGESV',& 
         __LINE__,__FILE__)
    CALL dcopy(n,fc,1,asl,1)
    DO i=1,m
       IF (il(i).NE.0) THEN
          fc(i)=asl(il(i),1)
       ELSE
          fc(i)=0._real_8
       ENDIF
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE solvs
  ! ==================================================================

END MODULE shake_utils
