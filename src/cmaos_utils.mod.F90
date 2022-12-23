MODULE cmaos_utils
  USE atwf,                            ONLY: atwf_mod,&
                                             atwp
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: paral
  USE prop,                            ONLY: prop2
  USE utils,                           ONLY: dspevy,&
                                             invmat
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: cmaos
  PUBLIC :: satch
  PUBLIC :: satch3
  PUBLIC :: satch4
  !public :: tracedp
  PUBLIC :: give_scr_cmaos
  PUBLIC :: give_scr_satch

CONTAINS

  ! ==================================================================
  SUBROUTINE cmaos(d,s,numin,unac)
    ! ==--------------------------------------------------------------==
    ! == CALCULATE MODIFIED ATOMIC ORBITALS                           ==
    ! ==--------------------------------------------------------------==

    REAL(real_8)                             :: d(atwp%nattot,*), &
                                                s(atwp%nattot,*)
    INTEGER                                  :: numin
    REAL(real_8)                             :: unac

    CHARACTER(*), PARAMETER                  :: procedureN = 'cmaos'

    INTEGER                                  :: i, i1, ia, iat, idamax, ierr, &
                                                ii, ij, info, is, j, k, &
                                                namax, nn, no, nx
    REAL(real_8)                             :: bbb, bl0, bl1
    REAL(real_8), ALLOCATABLE                :: coll(:,:), t1(:), t2(:), &
                                                t3(:), t4(:)

    namax=MAX(atwp%numaormax,MAXVAL(prop2%maos(1:ions1%nsp)),atwp%nattot)
    ! TODO check stat
    ALLOCATE(coll(atwp%nattot, numin),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(t1(namax**2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(t2(namax**2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(t3(namax**2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(t4(namax**2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent) THEN
       WRITE(6,'(A)') '*** MODIFIED ATOMIC ORBITALS ***'
       WRITE(6,'(A)') ' ATOM  OCCUPATION NUMBERS     '
    ENDIF
    CALL zeroing(coll)!,atwp%nattot*numin)
    iat=0
    no=1
    nx=0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          iat=iat+1
          nn=atwf_mod%numaor(is)
          CALL dgemm('N','N',nn,nn,nn,1._real_8,d(no,no),atwp%nattot,s(no,no),&
               atwp%nattot,0._real_8,t1(1),nn)
          k=0
          DO j=1,nn
             DO i=j,nn
                ij=(j-1)*nn+i
                k=k+1
                t2(k)=t1(ij)
             ENDDO
          ENDDO
          CALL dspevy(1,t2,t3,t1,nn,nn,t4,3*nn)
          IF (paral%io_parent)&
               WRITE(6,'(I4,5X,9(F7.2),/,9X,9(F7.2))')&
               iat,(t3(i),i=nn,nn-prop2%maos(is)+1,-1)
          bl0=t3(nn-prop2%maos(is)+1)
          IF ( idamax(nn-prop2%maos(is),t3(1),1) > 0 ) THEN
             bl1=t3(idamax(nn-prop2%maos(is),t3(1),1))
          ELSE
             bl1=0.0_real_8
          ENDIF
          bbb=ABS(bl1)/ABS(bl0)
          IF (bbb.GT.0.2_real_8) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(A,10(F7.2),/,9X,10(F7.2))')&
                  'WARNING :',(t3(i),i=nn-prop2%maos(is),1,-1)
          ENDIF
          DO i=1,prop2%maos(is)
             ii=nn-i+1
             i1=(ii-1)*nn+1
             CALL dcopy(nn,t1(i1),1,coll(no,nx+i),1)
          ENDDO
          no=no+nn
          nx=nx+prop2%maos(is)
       ENDDO
    ENDDO
    ! Calculate N=SDS
    CALL dgemm('N','N',atwp%nattot,atwp%nattot,atwp%nattot,1._real_8,d(1,1),atwp%nattot,&
         s(1,1),atwp%nattot,0._real_8,t1(1),atwp%nattot)
    CALL dgemm('N','N',atwp%nattot,atwp%nattot,atwp%nattot,1._real_8,s(1,1),atwp%nattot,&
         t1(1),atwp%nattot,0._real_8,d(1,1),atwp%nattot)
    ! Calculate unassigned charge : N - TR(DP)
    CALL dgemm('N','N',atwp%nattot,numin,atwp%nattot,1._real_8,s(1,1),atwp%nattot,&
         coll(1,1),atwp%nattot,0._real_8,t1(1),atwp%nattot)
    CALL dgemm('T','N',numin,numin,atwp%nattot,1._real_8,coll(1,1),atwp%nattot,&
         t1(1),atwp%nattot,0._real_8,s(1,1),numin)
    CALL invmat(numin,s,t1,info)
    CALL dgemm('N','N',atwp%nattot,numin,numin,1._real_8,coll(1,1),atwp%nattot,&
         s(1,1),numin,0._real_8,t1(1),atwp%nattot)
    CALL dgemm('N','T',atwp%nattot,atwp%nattot,numin,1._real_8,t1(1),atwp%nattot,&
         coll(1,1),atwp%nattot,0._real_8,s(1,1),atwp%nattot)
    DO i=1,atwp%nattot
       DO j=1,atwp%nattot
          unac=unac+d(i,j)*s(i,j)
       ENDDO
    ENDDO
    DO is=1,ions1%nsp
       unac=unac-ions0%zv(is)*ions0%na(is)
    ENDDO
    unac=-unac
    IF (paral%io_parent)&
         WRITE(6,'(/,A,F8.3,/)')&
         ' Modified Atomic Orbitals: Unassigned Charge ',unac
    CALL dcopy(atwp%nattot*numin,coll(1,1),1,s(1,1),1)
    ! ==--------------------------------------------------------------==
    ! TODO check stat
    DEALLOCATE(coll,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(t1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(t2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(t3,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(t4,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE cmaos
  ! ==================================================================
  SUBROUTINE satch(d,s,ca,numin,sac)
    ! ==--------------------------------------------------------------==
    ! == CALCULATE 2-CENTER SHARED ATOMIC CHARGES                     ==
    ! ==--------------------------------------------------------------==

    REAL(real_8)                             :: d(atwp%nattot,*), &
                                                s(atwp%nattot,*), &
                                                ca(atwp%nattot,*)
    INTEGER                                  :: numin
    REAL(real_8)                             :: sac(ions1%nat,ions1%nat)

    CHARACTER(*), PARAMETER                  :: procedureN = 'satch'
    INTEGER, PARAMETER                       :: listm = 100 

    INTEGER                                  :: i, ia1, ia2, iat1, iat2, &
                                                ierr, il, im1, im2, is1, is2, &
                                                j, list(listm), no1, no2
    REAL(real_8), ALLOCATABLE                :: t1(:,:), t2(:), t3(:,:)

! ==--------------------------------------------------------------==
! TODO check stat

    ALLOCATE(t1(atwp%nattot, numin),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(t2(atwp%nattot * numin),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(t3(atwp%nattot, atwp%nattot),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL zeroing(sac)!,ions1%nat*ions1%nat)
    iat1=0
    no1=1
    DO is1=1,ions1%nsp
       DO ia1=1,ions0%na(is1)
          iat1=iat1+1
          iat2=0
          no2=1
          DO is2=1,ions1%nsp
             DO ia2=1,ions0%na(is2)
                iat2=iat2+1
                IF (iat2.GE.iat1) THEN
                   il=0
                   DO im1=no1,no1+prop2%maos(is1)-1
                      il=il+1
                      list(il)=im1
                   ENDDO
                   IF (iat1.NE.iat2) THEN
                      DO im2=no2,no2+prop2%maos(is2)-1
                         il=il+1
                         list(il)=im2
                      ENDDO
                   ENDIF
                   IF (il.GT.listm) CALL stopgm('SATCH','LISTM',& 
                        __LINE__,__FILE__)
                   DO i=1,il
                      CALL dcopy(atwp%nattot,ca(1,list(i)),1,t1(1,i),1)
                   ENDDO
                   CALL tracedp(atwp%nattot,il,sac(iat1,iat2),d,s,t1,t2,t3)
                ENDIF
                no2=no2+prop2%maos(is2)
             ENDDO
          ENDDO
          no1=no1+prop2%maos(is1)
       ENDDO
    ENDDO
    DO i=1,ions1%nat
       DO j=1,i
          sac(i,j)=sac(j,i)
       ENDDO
    ENDDO
    DO i=1,ions1%nat
       DO j=1,ions1%nat
          sac(i,j)=sac(i,i)+sac(j,j)-sac(i,j)
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    ! TODO check stat
    DEALLOCATE(t1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(t2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(t3,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE satch
  ! ==================================================================
  SUBROUTINE satch3(d,s,ca,numin,sac,sac3)
    ! ==--------------------------------------------------------------==
    ! == CALCULATE 3-CENTER SHARED ATOMIC CHARGES                     ==
    ! ==--------------------------------------------------------------==

    REAL(real_8)                             :: d(atwp%nattot,*), &
                                                s(atwp%nattot,*), &
                                                ca(atwp%nattot,*)
    INTEGER                                  :: numin
    REAL(real_8) :: sac(ions1%nat,ions1%nat), &
      sac3(ions1%nat*ions1%nat*ions1%nat/6 + 1)

    CHARACTER(*), PARAMETER                  :: procedureN = 'satch3'
    INTEGER, PARAMETER                       :: listm = 100 

    INTEGER                                  :: i, ia1, ia123, ia2, ia3, &
                                                iat1, iat2, iat3, ierr, il, &
                                                im1, im2, im3, is1, is2, is3, &
                                                j, k, list(listm), n3, no1, &
                                                no2, no3
    REAL(real_8), ALLOCATABLE                :: t1(:,:), t2(:), t3(:,:)

! ==--------------------------------------------------------------==
! ==--------------------------------------------------------------==
! TODO check stat

    ALLOCATE(t1(atwp%nattot, numin),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(t2(atwp%nattot*numin),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(t3(atwp%nattot, atwp%nattot),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    n3=ions1%nat*ions1%nat*ions1%nat/6 + 1
    CALL zeroing(sac3)!,n3)
    ia123=0
    iat1=0
    no1=1
    DO is1=1,ions1%nsp
       DO ia1=1,ions0%na(is1)
          iat1=iat1+1
          iat2=0
          no2=1
          DO is2=1,ions1%nsp
             DO ia2=1,ions0%na(is2)
                iat2=iat2+1
                iat3=0
                no3=1
                DO is3=1,ions1%nsp
                   DO ia3=1,ions0%na(is3)
                      iat3=iat3+1
                      IF ((iat3.GT.iat2).AND.(iat2.GT.iat1)) THEN
                         il=0
                         DO im1=no1,no1+prop2%maos(is1)-1
                            il=il+1
                            list(il)=im1
                         ENDDO
                         DO im2=no2,no2+prop2%maos(is2)-1
                            il=il+1
                            list(il)=im2
                         ENDDO
                         DO im3=no3,no3+prop2%maos(is3)-1
                            il=il+1
                            list(il)=im3
                         ENDDO
                         IF (il.GT.listm) CALL stopgm('SATCH3','LISTM',& 
                              __LINE__,__FILE__)
                         DO i=1,il
                            CALL dcopy(atwp%nattot,ca(1,list(i)),1,t1(1,i),1)
                         ENDDO
                         ia123=ia123+1
                         CALL tracedp(atwp%nattot,il,sac3(ia123),d,s,t1,t2,t3)
                      ENDIF
                      no3=no3+prop2%maos(is3)
                   ENDDO
                ENDDO
                no2=no2+prop2%maos(is2)
             ENDDO
          ENDDO
          no1=no1+prop2%maos(is1)
       ENDDO
    ENDDO
    ia123=0
    DO i=1,ions1%nat
       DO j=i+1,ions1%nat
          DO k=j+1,ions1%nat
             ia123=ia123+1
             sac3(ia123)=sac3(ia123)+sac(i,j)+sac(i,k)+sac(j,k)-&
                  sac(i,i)-sac(j,j)-sac(k,k)
          ENDDO
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    ! TODO check stat
    DEALLOCATE(t1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(t2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(t3,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE satch3
  ! ==================================================================
  SUBROUTINE satch4(d,s,ca,numin,sac,sac3,sac4)
    ! ==--------------------------------------------------------------==
    ! == CALCULATE 4-CENTER SHARED ATOMIC CHARGES                     ==
    ! ==--------------------------------------------------------------==

    REAL(real_8)                             :: d(atwp%nattot,*), &
                                                s(atwp%nattot,*), &
                                                ca(atwp%nattot,*)
    INTEGER                                  :: numin
    REAL(real_8) :: sac(ions1%nat,ions1%nat), &
      sac3(ions1%nat*ions1%nat*ions1%nat/6 + 1), &
      sac4(ions1%nat*ions1%nat*ions1%nat*ions1%nat/24 + 1)

    CHARACTER(*), PARAMETER                  :: procedureN = 'satch4'
    INTEGER, PARAMETER                       :: listm = 100 

    INTEGER :: i, i1, i2, i3, ia1, ia1234, ia2, ia3, ia4, iat1, iat2, iat3, &
      iat4, ierr, il, im1, im2, im3, im4, is1, is2, is3, is4, ix, j, k, l, &
      list(listm), n4, no1, no2, no3, no4
    REAL(real_8), ALLOCATABLE                :: t1(:,:), t2(:), t3(:,:)

! ==--------------------------------------------------------------==
! ==--------------------------------------------------------------==
! TODO check stat

    ALLOCATE(t1(atwp%nattot, numin),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(t2(atwp%nattot*numin),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(t3(atwp%nattot, atwp%nattot),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    n4=ions1%nat*ions1%nat*ions1%nat*ions1%nat/24 + 1
    CALL zeroing(sac4)!,n4)
    ia1234=0
    iat1=0
    no1=1
    DO is1=1,ions1%nsp
       DO ia1=1,ions0%na(is1)
          iat1=iat1+1
          iat2=0
          no2=1
          DO is2=1,ions1%nsp
             DO ia2=1,ions0%na(is2)
                iat2=iat2+1
                iat3=0
                no3=1
                DO is3=1,ions1%nsp
                   DO ia3=1,ions0%na(is3)
                      iat3=iat3+1
                      iat4=0
                      no4=1
                      DO is4=1,ions1%nsp
                         DO ia4=1,ions0%na(is4)
                            iat4=iat4+1
                            IF ((iat4.GT.iat3).AND.(iat3.GT.iat2).AND.&
                                 (iat2.GT.iat1)) THEN
                               il=0
                               DO im1=no1,no1+prop2%maos(is1)-1
                                  il=il+1
                                  list(il)=im1
                               ENDDO
                               DO im2=no2,no2+prop2%maos(is2)-1
                                  il=il+1
                                  list(il)=im2
                               ENDDO
                               DO im3=no3,no3+prop2%maos(is3)-1
                                  il=il+1
                                  list(il)=im3
                               ENDDO
                               DO im4=no4,no4+prop2%maos(is4)-1
                                  il=il+1
                                  list(il)=im4
                               ENDDO
                               ia1234=ia1234+1
                               DO i=1,il
                                  CALL dcopy(atwp%nattot,ca(1,list(i)),1,t1(1,i),1)
                               ENDDO
                               CALL tracedp(atwp%nattot,il,sac4(ia1234),&
                                    d,s,t1,t2,t3)
                            ENDIF
                            no4=no4+prop2%maos(is4)
                         ENDDO
                      ENDDO
                      no3=no3+prop2%maos(is3)
                   ENDDO
                ENDDO
                no2=no2+prop2%maos(is2)
             ENDDO
          ENDDO
          no1=no1+prop2%maos(is1)
       ENDDO
    ENDDO
    ia1234=0
    DO i=1,ions1%nat
       DO j=i+1,ions1%nat
          DO k=j+1,ions1%nat
             DO l=k+1,ions1%nat
                ia1234=ia1234+1
                sac4(ia1234)=-sac4(ia1234)&
                     +sac(i,i)+sac(j,j)+sac(k,k)+sac(l,l)&
                     -sac(i,j)-sac(i,k)-sac(i,l)-sac(j,k)-sac(j,l)-sac(k,l)
                ix=0
                DO i1=1,ions1%nat
                   DO i2=i1+1,ions1%nat
                      DO i3=i2+1,ions1%nat
                         ix=ix+1
                         IF (i1.EQ.i) THEN
                            IF (i2.EQ.j) THEN
                               IF (i3.EQ.k) sac4(ia1234)=sac4(ia1234)+sac3(ix)
                               IF (i3.EQ.l) sac4(ia1234)=sac4(ia1234)+sac3(ix)
                            ELSEIF (i2.EQ.k.AND.i3.EQ.l) THEN
                               sac4(ia1234)=sac4(ia1234)+sac3(ix)
                            ENDIF
                         ELSEIF (i1.EQ.j) THEN
                            IF (i2.EQ.k.AND.i3.EQ.l)&
                                 sac4(ia1234)=sac4(ia1234)+sac3(ix)
                         ENDIF
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    ! TODO check stat
    DEALLOCATE(t1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(t2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(t3,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE satch4
  ! ==================================================================
  SUBROUTINE tracedp(n,m,sac,d,s,c,t2,t3)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: n, m
    REAL(real_8)                             :: sac, d(n,n), s(n,n), c(n,m), &
                                                t2(n*n), t3(n,n)

    INTEGER                                  :: i, info, j

    CALL dgemm('N','N',n,m,n,1._real_8,s(1,1),n,c(1,1) ,n,0._real_8,t2(1)  ,n)
    CALL dgemm('T','N',m,m,n,1._real_8,c(1,1),n,t2(1)  ,n,0._real_8,t3(1,1),m)
    CALL invmat(m,t3,t2,info)
    CALL dgemm('N','N',n,m,m,1._real_8,c(1,1),n,t3(1,1),m,0._real_8,t2(1)  ,n)
    CALL dgemm('N','T',n,n,m,1._real_8,t2(1) ,n,c(1,1) ,n,0._real_8,t3(1,1),n)
    DO j=1,n
       DO i=1,n
          sac=sac+d(i,j)*t3(i,j)
       ENDDO
    ENDDO
    RETURN
  END SUBROUTINE tracedp
  ! ==================================================================
  SUBROUTINE give_scr_cmaos(lcmaos,tag,numin)
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: lcmaos
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: numin

    INTEGER                                  :: in, it1, it2, it3, it4, it5, &
                                                itm, namax

    namax=5
    DO in=1,ions1%nsp
       namax=MAX(namax,atwf_mod%numaor(in))
    ENDDO
    it1 = 1
    it2 = it1 + MAX(atwp%nattot*numin,1000)
    it3 = it2 + MAX(namax*namax,1000)
    it4 = it3 + namax*namax
    it5 = it4 + namax*namax
    itm = it5 + namax*namax
    lcmaos=MAX(itm,it1+atwp%nattot*numin+atwp%nattot*atwp%nattot) - 1
    tag  ='MAX(ITM,IT1+NATTOT*NUMIN+...)'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_cmaos
  ! ==================================================================
  SUBROUTINE give_scr_satch(lsatch,tag,numin)
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: lsatch
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: numin

! ==--------------------------------------------------------------==

    lsatch=2*atwp%nattot*numin+atwp%nattot*atwp%nattot
    tag  ='2*NATTOT*NUMIN+NATTOT*NATTOT'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_satch
  ! ==================================================================

END MODULE cmaos_utils
