MODULE clinbcg_utils
  USE atimes_utils,                    ONLY: atimes,&
                                             give_scr_atimes
  USE atimesmod,                       ONLY: atimes_eval
  USE cppt,                            ONLY: hg
  USE dotp_utils,                      ONLY: dotp
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE machine,                         ONLY: m_flush
  USE parac,                           ONLY: paral
  USE system,                          ONLY: fpar,&
                                             ncpw,&
                                             nkpt,&
                                             parm

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: clinbcg
  !public :: diago
  PUBLIC :: give_scr_clinbcg
  !public :: snrmm

CONTAINS

  ! ==================================================================
  SUBROUTINE clinbcg(b,x,c0,vpot,psi,nstate,f,itol,tol,&
       itmax,err)
    ! ==--------------------------------------------------------------==
    ! == Conjugate gradient method for a sparse system                ==
    ! == Numerical recipes p 79 2nd edition (LINBCG routine)          ==
    ! == (A. Alavi 1997)                                              ==
    ! == USES ATIMES,DIAGO,SNRMM                                      ==
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8)                          :: b(*), x(*), c0(*)
    REAL(real_8)                             :: vpot(fpar%nnr1,*)
    COMPLEX(real_8)                          :: psi(:)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: f(*)
    INTEGER                                  :: itol
    REAL(real_8)                             :: tol
    INTEGER                                  :: itmax
    REAL(real_8)                             :: err

    CHARACTER(*), PARAMETER                  :: procedureN = 'clinbcg'
    COMPLEX(real_8), PARAMETER               :: zzero = (0._real_8,0._real_8) 
    REAL(real_8), PARAMETER                  :: eps = 1.e-14_real_8

    COMPLEX(real_8)                          :: ak, akden, bk, bkden, bknum, &
                                                zdotc
    COMPLEX(real_8), ALLOCATABLE, SAVE       :: p(:), pp(:), r(:), rr(:), &
                                                z(:), zz(:)
    EXTERNAL                                 :: zdotc
    INTEGER                                  :: ierr, iter, j
    INTEGER, SAVE                            :: ifirst = 0
    LOGICAL                                  :: tdebug
    REAL(real_8)                             :: bnrm, dxnrm, xnrm, zm1nrm, &
                                                znrm

    tdebug=.TRUE.
    iter=0
    IF (ifirst.EQ.0) THEN
       ALLOCATE(r(nkpt%ngwk),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(rr(nkpt%ngwk),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(z(nkpt%ngwk),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(zz(nkpt%ngwk),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(p(nkpt%ngwk),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(pp(nkpt%ngwk),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ifirst=1
    ENDIF
    CALL atimes(x,r,c0,vpot,psi,nstate,f)
    !$omp parallel do private(J)
    DO j=1,nkpt%ngwk
       r(j)=b(j)-r(j)
       rr(j)=r(j)
    ENDDO
    ! CALL ATIMES(R,RR,C0,NSTATE,F)
    znrm=1._real_8
    IF (itol.EQ.1) THEN
       bnrm=snrmm(nkpt%ngwk,b)
    ELSE IF (itol.EQ.2) THEN
       CALL diago(b,z)
       bnrm=snrmm(nkpt%ngwk,z)
    ELSE IF (itol.EQ.3.OR.itol.EQ.4) THEN
       CALL diago(b,z)
       bnrm=snrmm(nkpt%ngwk,z)
       CALL diago(r,z)
       znrm=snrmm(nkpt%ngwk,z)
    ELSE
       CALL stopgm('CLINBCG','ILLEGAL ITOL',& 
            __LINE__,__FILE__)
    ENDIF
    CALL diago(r,z)
    ! ==--------------------------------------------------------------==
100 CONTINUE
    IF (iter.LE.itmax) THEN
       iter=iter+1
       zm1nrm=znrm
       CALL diago(rr,zz)
       IF (tkpts%tkpnt) THEN
          bknum=zdotc(nkpt%ngwk,rr,1,z,1)
       ELSE
          bknum=CMPLX(dotp(nkpt%ngwk,rr,z),0._real_8,kind=real_8)
       ENDIF
       IF (iter.EQ.1) THEN
          !$omp parallel do private(J)
          DO j=1,nkpt%ngwk
             p(j)=z(j)
             pp(j)=zz(j)
          ENDDO
          bkden=0._real_8        ! Not useful
       ELSE
          bk=bknum/bkden
          !$omp parallel do private(J)
          DO j=1,nkpt%ngwk
             p(j)=bk*p(j)+z(j)
             pp(j)=CONJG(bk)*pp(j)+zz(j)
          ENDDO
       ENDIF
       bkden=bknum
       CALL atimes(p,z,c0,vpot,psi,nstate,f)
       akden=zzero
       IF (tkpts%tkpnt) THEN
          akden = zdotc(nkpt%ngwk,pp,1,z,1)
       ELSE
          akden=CMPLX(dotp(nkpt%ngwk,pp,z),0._real_8,kind=real_8)
       ENDIF
       ak=bknum/akden
       CALL atimes(pp,zz,c0,vpot,psi,nstate,f)
       !$omp parallel do private(J)
       DO j=1,nkpt%ngwk
          x(j)=x(j)+ak*p(j)
          r(j)=r(j)-ak*z(j)
          rr(j)=rr(j)-CONJG(ak)*zz(j)
       ENDDO
       CALL diago(r,z)
       IF (itol.EQ.1.OR.itol.EQ.2)THEN
          znrm=1._real_8
          err=snrmm(nkpt%ngwk,r)
          err=err/bnrm
       ELSE IF (itol.EQ.3.OR.itol.EQ.4)THEN
          znrm=snrmm(nkpt%ngwk,z)
          IF (ABS(zm1nrm-znrm).GT.eps*znrm) THEN
             dxnrm=ABS(ak)*snrmm(nkpt%ngwk,p)
             err=znrm/ABS(zm1nrm-znrm)*dxnrm
          ELSE
             err=znrm/bnrm
             GOTO 100
          ENDIF
          xnrm=snrmm(nkpt%ngwk,x)
          IF (err.LE.0.5_real_8*xnrm) THEN
             err=err/xnrm
          ELSE
             err=znrm/bnrm
             GOTO 100
          ENDIF
       ENDIF
       IF (tdebug) THEN
          IF (paral%io_parent)&
               WRITE(6,*) ' ITER=',iter,' ERR=',err
          CALL m_flush(6)
       ENDIF
       IF (err.GT.tol) GOTO 100
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE clinbcg
  ! ==================================================================
  FUNCTION snrmm(n,z)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: n
    COMPLEX(real_8)                          :: z(n)
    REAL(real_8)                             :: snrmm

    EXTERNAL                                 :: dotp
    REAL(real_8)                             :: dotp, dznrm2

! ==--------------------------------------------------------------==

    IF (.NOT.tkpts%tkpnt) THEN
       snrmm = dotp(n,z,z)
       snrmm = SQRT(snrmm)
    ELSE
       snrmm=dznrm2(n,z,1)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION snrmm
  ! ==================================================================
  SUBROUTINE diago(c0,c1)
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8)                          :: c0(nkpt%ngwk), c1(nkpt%ngwk)

    INTEGER                                  :: ig

! VAriables
! ==--------------------------------------------------------------==

    IF (.NOT.tkpts%tkpnt) THEN
       !$omp parallel do private(IG)
       DO ig=1,ncpw%ngw
          c1(ig) = c0(ig)/(hg(ig)*parm%tpiba2*0.5_real_8-atimes_eval)
       ENDDO
    ELSE
       !$omp parallel do private(IG)
       DO ig=1,ncpw%ngw
          c1(ig) = c0(ig)/(hg(ig)*parm%tpiba2*0.5_real_8-atimes_eval)
          c1(ncpw%ngw+ig) = c0(ncpw%ngw+ig)/(hg(ig)*parm%tpiba2*0.5_real_8-atimes_eval)
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE diago
  ! ==================================================================
  SUBROUTINE give_scr_clinbcg(lclinbcg,tag,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lclinbcg
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate

! ==--------------------------------------------------------------==

    CALL give_scr_atimes(lclinbcg,tag,nstate)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_clinbcg
  ! ==================================================================

END MODULE clinbcg_utils
