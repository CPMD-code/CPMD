MODULE orbrot_utils
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_bcast
  USE parac,                           ONLY: parai,&
                                             paral
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: fc4,&
                                             sbes0
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: get_u
  PUBLIC :: xgrad

CONTAINS

  ! ==================================================================
  SUBROUTINE get_u(urot,umat,pmat,peig,nstate,msub)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: peig(*)
    INTEGER                                  :: nstate, msub
    REAL(real_8)                             :: urot(nstate,msub), &
                                                umat(nstate,msub), &
                                                pmat(msub,msub)

    CHARACTER(*), PARAMETER                  :: procedureN = 'get_u'

    INTEGER                                  :: i, ierr, ij, info, j, k, &
                                                lwork, mm
    REAL(real_8)                             :: fek, ffk, sqe, uu
    REAL(real_8), ALLOCATABLE                :: work(:)

    lwork = MAX(3*msub-1, msub*msub)
    ALLOCATE(work(lwork),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    IF (paral%parent) THEN
       mm = nstate - msub
       CALL dgemm('T','N',msub,msub,mm,1._real_8,urot,nstate,urot,nstate,&
            0._real_8,pmat,msub)
       CALL dsyev('V','L',msub,pmat,msub,peig,work,lwork,info)
       IF (info.NE.0) CALL stopgm('GET_U','INFO ! =0 AFTER DSYEV',& 
            __LINE__,__FILE__)
       CALL zeroing(umat)!,msub*nstate)
       CALL zeroing(work)!,msub*msub)
       DO k=1,msub
          sqe=SQRT(peig(k))
          fek=COS(sqe)
          ffk=sbes0(sqe)
          DO i=1,msub
             DO j=1,msub
                ij=j+(i-1)*msub
                umat(j,i)=umat(j,i)+pmat(i,k)*pmat(j,k)*fek
                work(ij)=work(ij)+pmat(i,k)*pmat(j,k)*ffk
             ENDDO
          ENDDO
       ENDDO
       CALL dgemm('N','N',mm,msub,msub,-1._real_8,urot,nstate,work,msub,&
            0._real_8,umat(msub+1,1),nstate)
       DO k=1,msub
          DO i=1,nstate/2
             j=nstate-i+1
             uu=umat(i,k)
             umat(i,k)=umat(j,k)
             umat(j,k)=uu
          ENDDO
       ENDDO
    ENDIF
    CALL mp_bcast(umat,SIZE(umat),parai%source,parai%allgrp)
    ! ==--------------------------------------------------------------==
    DEALLOCATE(work,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE get_u
  ! ==================================================================
  SUBROUTINE xgrad(x0,x2,c2,rmat,reig,nstate,msub)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate, msub
    REAL(real_8)                             :: x0(nstate,msub), &
                                                x2(nstate,msub), &
                                                c2(nstate,msub), reig(msub), &
                                                rmat(msub,msub)

    CHARACTER(*), PARAMETER                  :: procedureN = 'xgrad'

    INTEGER                                  :: i, ierr, ij, isub, j, ji, k, &
                                                l, nldx
    LOGICAL                                  :: debug
    REAL(real_8)                             :: a, AmB, ApB, b, f1, f2
    REAL(real_8), ALLOCATABLE                :: s1(:), s2(:), s3(:), s4(:), &
                                                s5(:)

! Variables
! ==--------------------------------------------------------------==
! ==  CALCULATE THE GRADIENT OF THE ENERGY WITH RESPECT TO        ==
! ==  ORBITAL ROTATIONS                                           ==
! ==                                                              ==
! ==  INPUT  : C2     ; GRADIENT WRT COEFFICIENTS                 ==
! ==         : RMAT*REIG*RMAT(T) = X*X(T)                         ==
! ==         : X0     ; CURRENT VALUE OF X                        ==
! ==  OUTPUT : X2     ; GRADIENT WRT ROTATIONS                    ==
! ==--------------------------------------------------------------==

    CALL tiset('     XGRAD',isub)
    IF (paral%parent) THEN
       debug=.FALSE.
       nldx=nstate-msub
       ALLOCATE(s1(msub*msub),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(s2(nstate*msub),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(s3(msub*msub),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(s4(msub*msub),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(s5(msub*msub),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ! ==--------------------------------------------------------------==
       ! ==  DECOMPOSE THE GRADIENT                                      ==
       ! ==--------------------------------------------------------------==
       DO k=1,msub
          DO i=1,nstate/2
             j=nstate-i+1
             a=c2(i,k)
             c2(i,k)=c2(j,k)
             c2(j,k)=a
          ENDDO
       ENDDO
       l=0
       k=0
       DO j=1,msub
          DO i=1,msub
             k=k+1
             s1(k)=c2(i,j)
          ENDDO
          DO i=msub+1,nstate
             l=l+1
             s2(l)=c2(i,j)
          ENDDO
       ENDDO
       ! ==--------------------------------------------------------------==
       ! ==  F2=R*X*G2*R(T)                                              ==
       ! ==--------------------------------------------------------------==
       CALL dgemm('T','N',msub,msub,nldx,1.0_real_8,x0,nstate,s2,nldx,0.0_real_8,&
            s3,msub)
       CALL dgemm('T','N',msub,msub,msub,1.0_real_8,rmat,msub,s3,msub,0.0_real_8,&
            s4,msub)
       CALL dgemm('N','N',msub,msub,msub,1.0_real_8,s4,msub,rmat,msub,0.0_real_8,&
            s3,msub)
       ! ==--------------------------------------------------------------==
       ! ==  F1=R*G1*R(T)                                                ==
       ! ==--------------------------------------------------------------==
       CALL dgemm('T','N',msub,msub,msub,1.0_real_8,rmat,msub,s1,msub,0.0_real_8,&
            s4,msub)
       CALL dgemm('N','N',msub,msub,msub,1.0_real_8,s4,msub,rmat,msub,0.0_real_8,&
            s1,msub)
       ! ==--------------------------------------------------------------==
       ! ==  F1+F1(T) ; F2+F2(T)                                         ==
       ! ==--------------------------------------------------------------==
       DO i=1,msub
          DO j=i,msub
             ij=(i-1)*msub+j
             ji=(j-1)*msub+i
             f1=s1(ij)+s1(ji)
             f2=s3(ij)+s3(ji)
             s1(ij)=f1
             s1(ji)=f1
             s3(ij)=f2
             s3(ji)=f2
          ENDDO
       ENDDO
       ! ==--------------------------------------------------------------==
       ! ==  C1; C2                                                      ==
       ! ==--------------------------------------------------------------==
       k=0
       DO i=1,msub
          DO j=1,msub
             k=k+1
             a=SQRT(reig(i))
             b=SQRT(reig(j))
             ApB=0.5_real_8*(a+b)
             AmB=0.5_real_8*(a-b)
             s4(k)=0.5_real_8*sbes0(ApB)*sbes0(AmB)
             s5(k)=fc4(a,b)
          ENDDO
       ENDDO
       ! ==--------------------------------------------------------------==
       ! ==  FC = - F1 x C1 - (F2 x C2)(T) + F2 x C2                     ==
       ! ==--------------------------------------------------------------==
       k=0
       DO i=1,msub
          DO j=1,msub
             k=k+1
             s1(k)=-s1(k)*s4(k)-s3(k)*s5(k)
          ENDDO
       ENDDO
       ! ==--------------------------------------------------------------==
       ! ==  X2 = R(T) * FC * R * X - J0(P)*G2(T)                        ==
       ! ==--------------------------------------------------------------==
       CALL dgemm('N','N',msub,msub,msub,1.0_real_8,rmat,msub,s1,msub,0.0_real_8,&
            s3,msub)
       CALL dgemm('N','T',msub,msub,msub,1.0_real_8,s3,msub,rmat,msub,0.0_real_8,&
            s1,msub)
       CALL dcopy(msub*msub,rmat,1,s3,1)
       DO i=1,msub
          DO j=1,msub
             rmat(i,j)=rmat(i,j)*sbes0(SQRT(reig(j)))
          ENDDO
       ENDDO
       CALL dgemm('N','T',msub,msub,msub,1.0_real_8,rmat,msub,s3,msub,0.0_real_8,&
            s4,msub)
       CALL zeroing(x2)!,nstate*msub)
       CALL dgemm('N','T',nldx,msub,msub,-1.0_real_8,s2,nldx,s4,msub,0.0_real_8,&
            x2,nstate)
       CALL dgemm('N','T',nldx,msub,msub,1.0_real_8,x0,nstate,s1,msub,1.0_real_8,&
            x2,nstate)
       ! ==--------------------------------------------------------------==
       DEALLOCATE(s1,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(s2,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(s3,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(s4,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(s5,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    ENDIF
    CALL mp_bcast(x2,SIZE(x2),parai%source,parai%allgrp)
    CALL tihalt('     XGRAD',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE xgrad
  ! ==================================================================

END MODULE orbrot_utils
