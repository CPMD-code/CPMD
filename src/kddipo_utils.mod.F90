MODULE kddipo_utils
  USE ddipo_utils,                     ONLY: setdip
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE kpnt,                            ONLY: rk
  USE kpts,                            ONLY: kpts_com
  USE opeigr_c_utils,                  ONLY: opeigr_c
  USE opeigr_utils,                    ONLY: give_scr_opeigr
  USE parac,                           ONLY: parai,&
                                             paral
  USE rggen_utils,                     ONLY: recips
  USE system,                          ONLY: nkpt,&
                                             parap,&
                                             parm,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
!!use ovlap_utils, only : ovlap_c

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: kddipo
  PUBLIC :: give_scr_kddipo

CONTAINS

  ! ==================================================================
  SUBROUTINE kddipo(tau0,c0,cm,c2,sc0,nstate)
    ! ==--------------------------------------------------------------==
    ! == CALCULATES THE BERRY PHASE POLAIZATION                       ==
    ! ==--------------------------------------------------------------==
    ! input
    REAL(real_8)                             :: tau0(:,:,:)
    COMPLEX(real_8)                          :: cm(*), c2(*), sc0(*)
    INTEGER                                  :: nstate
    COMPLEX(real_8) :: c0(nkpt%ngwk,nstate,nkpt%nkpnt)

    CHARACTER(*), PARAMETER                  :: procedureN = 'kddipo'

    COMPLEX(real_8)                          :: dd
    COMPLEX(real_8), ALLOCATABLE             :: mkl(:,:,:,:)
    INTEGER                                  :: i, ierr, isub, j, k, l, n1, &
                                                ngwmax, nmcol
    INTEGER, ALLOCATABLE                     :: knn(:,:), mapcol(:), mapful(:)
    REAL(real_8)                             :: b1(3), b2(3), b3(3), del1(3), &
                                                del2(3), del3(3), delt, &
                                                r1(3), r2(3), r3(3)

! ==--------------------------------------------------------------==

    CALL tiset('    KDDIPO',isub)

    ALLOCATE(knn(nkpt%nkpnt,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(mkl(nstate,nstate,nkpt%nkpnt,6),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(mapful(2*spar%ngwks),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    n1=0
    DO i=0,parai%nproc-1
       n1=MAX(n1,parap%sparm(3,i))
    ENDDO
    ngwmax=n1
    nmcol=parai%nproc*ngwmax
    ALLOCATE(mapcol(nmcol),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL setdip(mapful,mapcol)

    CALL recips(parm%alat,parm%a1,parm%a2,parm%a3,b1,b2,b3)
    delt=0.0000001_real_8

    DO j=1,3
       b1(j)=b1(j)/kpts_com%nk1
       b2(j)=b2(j)/kpts_com%nk2
       b3(j)=b3(j)/kpts_com%nk3
       ! WRITE(6,*)'B1=',J,B1(J)
       ! WRITE(6,*)'B2=',J,B2(J)
       ! WRITE(6,*)'B3=',J,B3(J)
    ENDDO

    DO i=1,nkpt%nkpnt
       DO j=1,3
          knn(i,j)=0
       ENDDO
    ENDDO

    DO i=1,nkpt%nkpnt
       ! WRITE(6,*)'--------------------------'
       DO j=1,3
          r1(j)=rk(j,i)+b1(j)
          r2(j)=rk(j,i)+b2(j)
          r3(j)=rk(j,i)+b3(j)
          ! WRITE(6,*)'J=',J,'R1=',R1(J),'R2=',R2(J),'R3=',R3(J)
       ENDDO

       DO k=1,nkpt%nkpnt
          ! WRITE(6,*)'I=',I,'-------------------------------------'
          ! WRITE(6,*)'K=',K,'RK=',RK(1,K),RK(2,K),RK(3,K)
          ! WRITE(6,*)'K=',K,'R1=',R1(1),R1(2),R1(3)
          ! WRITE(6,*)'K=',K,'R2=',R2(1),R2(2),R2(3)
          ! WRITE(6,*)'K=',K,'R3=',R3(1),R3(2),R3(3)
          DO j=1,3

             del1(j)=ABS(rk(j,k)-r1(j))
             del2(j)=ABS(rk(j,k)-r2(j))
             del3(j)=ABS(rk(j,k)-r3(j))

          ENDDO

          IF (del1(1).LT.delt.AND.del1(2).LT.delt.AND.del1(3).LT.delt)&
               THEN
             knn(i,1)=k
          ENDIF
          IF (del2(1).LT.delt.AND.del2(2).LT.delt.AND.del2(3).LT.delt)&
               THEN
             knn(i,2)=k
          ENDIF
          IF (del3(1).LT.delt.AND.del3(2).LT.delt.AND.del3(3).LT.delt)&
               THEN
             knn(i,3)=k
          ENDIF

       ENDDO
    ENDDO

    IF (paral%io_parent)&
         WRITE(6,*)'kpoints'


    DO i=1,nkpt%nkpnt
       IF (paral%io_parent)&
            WRITE(6,*)'I=',i,'KN1=',knn(i,1)
       IF (knn(i,1).EQ.0)THEN
          IF (paral%io_parent)&
               WRITE(6,*)'last point in B1 direction'
       ENDIF
       IF (paral%io_parent)&
            WRITE(6,*)'I=',i,'KN2=',knn(i,2)
       IF (knn(i,2).EQ.0)THEN
          IF (paral%io_parent)&
               WRITE(6,*)'last point in B2 direction'
       ENDIF
       IF (paral%io_parent)&
            WRITE(6,*)'I=',i,'KN3=',knn(i,3)
       IF (knn(i,3).EQ.0)THEN
          IF (paral%io_parent)&
               WRITE(6,*)'last point in B3 direction'
       ENDIF
    ENDDO

    DO i=1,nkpt%nkpnt
       DO j=1,3
          l=knn(i,j)
          IF (l.EQ.0) THEN
             CALL opeigr_c(c0,c2,sc0,nstate,mapful,mapcol,mkl,j,&
                  0,dd)
          ELSE
             CALL ovlap_c(nstate,mkl,c0(1,1,i),c0(1,1,l))
          ENDIF
       ENDDO
    ENDDO
    IF (paral%io_parent)&
         WRITE(6,*)'overlap matrix'
    IF (paral%io_parent)&
         WRITE(6,*)'end kddipo'

    DEALLOCATE(mkl,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(knn,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(mapcol,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(mapful,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    CALL tihalt('    KDDIPO',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE kddipo

  ! ==================================================================
  SUBROUTINE give_scr_kddipo(lddipo,nstate,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lddipo, nstate
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: lopeigr

! ==--------------------------------------------------------------==

    CALL give_scr_opeigr(lopeigr,tag,nstate)
    lddipo=MAX(lopeigr,10)+100
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_kddipo
  ! ==================================================================

END MODULE kddipo_utils
