MODULE linalg_utils
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sendrecv
  USE system,                          ONLY: parap
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: dmatmul
  PUBLIC :: trans_da
  PUBLIC :: symm_da

CONTAINS

  ! ==================================================================
  SUBROUTINE dmatmul(ta,tb,amat,bmat,cmat,rmat,ldr,nstate,nproc,&
       mepos,ndd1,ndd2,grp)
    ! ==--------------------------------------------------------------==
    ! ==  Calculate matrix multiplication for distributed matrices    ==
    ! ==  Op(A) Op(B) = C                                             ==
    ! ==--------------------------------------------------------------==
    CHARACTER(len=1)                         :: ta, tb
    REAL(real_8)                             :: rmat(*)
    INTEGER                                  :: ldr, nstate
    REAL(real_8)                             :: cmat(nstate,*), &
                                                bmat(nstate,*), amat(nstate,*)
    INTEGER                                  :: nproc, mepos, ndd1(0:*), &
                                                ndd2(0:*), grp

    INTEGER                                  :: i1, ip, ipp, ir2, msglen, n, &
                                                norb, norbx

! ==--------------------------------------------------------------==

    norb=ndd2(mepos)-ndd1(mepos)+1
    norbx=0
    DO ip=0,nproc-1
       n=ndd2(ip)-ndd1(ip)+1
       norbx=MAX(norbx,n)
    ENDDO
    IF (ldr.LT.nstate*norbx) CALL stopgm('DMATMUL','LDR TOO SMALL',& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    msglen=nstate*norbx * 8
    CALL zeroing(cmat(:,1:norbx))!,nstate*norbx)
    IF (ta.EQ.'T'.OR.ta.EQ.'t') THEN
       IF (tb.EQ.'T'.OR.tb.EQ.'t') THEN
          ! ..A(T)*B(T)
          CALL stopgm('DMATMUL','NOT IMPLEMENTED',& 
               __LINE__,__FILE__)
       ELSE
          ! ..A(T)*B
          DO ip=0,nproc-1
             CALL my_shift(amat,rmat,msglen,mepos,-ip,grp)
             ipp=MOD(mepos+ip,nproc)
             n=ndd2(ipp)-ndd1(ipp)+1
             i1=ndd1(ipp)
             CALL dgemm('T','N',n,norb,nstate,1._real_8,rmat,nstate,&
                  bmat,nstate,1._real_8,cmat(i1,1),nstate)
          ENDDO
       ENDIF
    ELSE
       IF (tb.EQ.'T'.OR.tb.EQ.'t') THEN
          ! ..A*B(T)
          IF (ldr.LT.2*nstate*norbx)&
               CALL stopgm('DMATMUL','LDR TOO SMALL',& 
               __LINE__,__FILE__)
          ir2=nstate*norbx+1
          DO ip=0,nproc-1
             ipp=MOD(mepos+ip,nproc)
             n=ndd2(ipp)-ndd1(ipp)+1
             i1=ndd1(ipp)
             CALL dgemm('N','T',nstate,n,norb,1._real_8,amat,nstate,&
                  bmat(i1,1),nstate,0._real_8,rmat(ir2),nstate)
             CALL my_shift(rmat(ir2),rmat,msglen,mepos,ip,grp)
             CALL daxpy(nstate*norb,1._real_8,rmat,1,cmat,1)
          ENDDO
       ELSE
          ! ..A*B
          DO ip=0,nproc-1
             CALL my_shift(amat,rmat,msglen,mepos,-ip,grp)
             ipp=MOD(mepos+ip,nproc)
             n=ndd2(ipp)-ndd1(ipp)+1
             i1=ndd1(ipp)
             CALL dgemm('N','N',nstate,norb,n,1._real_8,rmat,nstate,&
                  bmat(i1,1),nstate,1._real_8,cmat,nstate)
          ENDDO
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE dmatmul
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE trans_da(tmat,n,ndd1,ndd2,nx,me,nproc,grp)
    INTEGER                                  :: n
    REAL(real_8)                             :: tmat(n,*)
    INTEGER                                  :: ndd1(0:*), ndd2(0:*), nx, me, &
                                                nproc, grp

    CHARACTER(*), PARAMETER                  :: procedureN = 'trans_da'

    INTEGER                                  :: i, ierr, ip, j, len, m, nl
    REAL(real_8), ALLOCATABLE                :: amat(:,:,:), bmat(:,:,:)

! ==--------------------------------------------------------------==

    ALLOCATE(amat(nx,nx,nproc),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(bmat(nx,nx,nproc),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    CALL zeroing(amat)!,nx*nx*nproc)
    nl=ndd2(me)-ndd1(me)+1
    DO ip=0,nproc-1
       m=ndd2(ip)-ndd1(ip)+1
       DO i=1,nl
          DO j=0,m-1
             amat(j+1,i,ip+1)=tmat(ndd1(ip)+j,i)
          ENDDO
       ENDDO
    ENDDO
    len=8*nx*nx
    CALL my_trans(amat,bmat,len,1)
    DO ip=0,nproc-1
       m=ndd2(ip)-ndd1(ip)+1
       DO i=1,nl
          DO j=0,m-1
             tmat(ndd1(ip)+j,i)=bmat(i,j+1,ip+1)
          ENDDO
       ENDDO
    ENDDO

    DEALLOCATE(amat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(bmat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE trans_da
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE symm_da(tmat,amat,bmat,n,ndd1,ndd2,nx,my,npr,grp)
    INTEGER                                  :: n
    REAL(real_8)                             :: tmat(n,*)
    INTEGER                                  :: ndd1(0:*), ndd2(0:*), nx
    REAL(real_8)                             :: bmat(nx,nx), amat(nx,nx)
    INTEGER                                  :: my, npr, grp

    INTEGER                                  :: i, ip, j, len, m, nl

! ==--------------------------------------------------------------==

    CALL zeroing(amat)!,nx*nx)
    nl=ndd2(my)-ndd1(my)+1
    len=nx*nx
    DO ip=0,npr-1
       m=ndd2(ip)-ndd1(ip)+1
       DO i=1,nl
          DO j=0,m-1
             amat(j+1,i)=tmat(ndd1(ip)+j,i)
          ENDDO
       ENDDO
       IF ((nl.GT.0).OR.(m.GT.0)) THEN
          CALL mp_sendrecv(amat,len,parap%pgroup(ip+1),bmat,&
               len,parap%pgroup(ip+1),grp)
       ENDIF
       m=ndd2(ip)-ndd1(ip)+1
       DO i=1,nl
          DO j=0,m-1
             tmat(ndd1(ip)+j,i)=0.5_real_8 *(tmat(ndd1(ip)+j,i)+bmat(i,j+1))
          ENDDO
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE symm_da
  ! ==================================================================

END MODULE linalg_utils
