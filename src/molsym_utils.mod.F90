MODULE molsym_utils
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: molsym

CONTAINS

  ! ==================================================================
  SUBROUTINE molsym(tag,naxis,tab,nop)
    ! ==--------------------------------------------------------------==
    CHARACTER(len=3)                         :: tag
    INTEGER                                  :: naxis
    REAL(real_8)                             :: tab(3,3,*)
    INTEGER                                  :: nop

    INTEGER                                  :: i, j, k, n1, ngen
    REAL(real_8)                             :: a, b, c, dd, dsop(3,3), &
                                                gen(3,3,5), pi, sop(3,3)
    REAL(real_8), EXTERNAL                   :: dasum

    ngen=0
    ! SET THE GENERATORS
    IF (INDEX(tag(1:1),'C').NE.0) THEN
       IF (naxis.NE.1) THEN
          ngen=ngen+1
          CALL genax(gen(1,1,ngen),naxis,3)
       ENDIF
       IF (INDEX(tag,'CNH').NE.0) THEN
          ngen=ngen+1
          CALL gensi(gen(1,1,ngen),3)
       ELSEIF (INDEX(tag,'CS').NE.0) THEN
          ngen=ngen+1
          CALL gensi(gen(1,1,ngen),3)
       ELSEIF (INDEX(tag,'CNV').NE.0) THEN
          ngen=ngen+1
          CALL gensi(gen(1,1,ngen),1)
       ENDIF
    ELSEIF (INDEX(tag(1:1),'S').NE.0) THEN
       ngen=ngen+1
       CALL gensn(gen(1,1,ngen),naxis,3)
    ELSEIF (INDEX(tag(1:1),'D').NE.0) THEN
       ngen=ngen+1
       CALL genax(gen(1,1,ngen),naxis,3)
       IF (INDEX(tag,'DNH').NE.0) THEN
          ngen=ngen+1
          CALL genax(gen(1,1,ngen),2,1)
          ngen=ngen+1
          CALL gensi(gen(1,1,ngen),3)
       ELSEIF (INDEX(tag,'DND').NE.0) THEN
          ngen=ngen+1
          CALL genax(gen(1,1,ngen),2,1)
          ngen=ngen+1
          CALL gend(gen(1,1,ngen),naxis,3)
       ELSEIF (INDEX(tag,'DN').NE.0) THEN
          ngen=ngen+1
          CALL genax(gen(1,1,ngen),2,1)
       ENDIF
    ELSEIF (INDEX(tag(1:1),'T').NE.0) THEN
       ngen=ngen+1
       CALL genax(gen(1,1,ngen),3,4)
       IF (INDEX(tag,'TD').NE.0) THEN
          ngen=ngen+1
          CALL gensn(gen(1,1,ngen),4,3)
       ELSE
          ngen=ngen+1
          CALL genax(gen(1,1,ngen),2,3)
          IF (INDEX(tag,'TH').NE.0) THEN
             ngen=ngen+1
             CALL geni(gen(1,1,ngen))
          ENDIF
       ENDIF
    ELSEIF (INDEX(tag(1:1),'O').NE.0) THEN
       ngen=ngen+1
       CALL genax(gen(1,1,ngen),4,3)
       ngen=ngen+1
       CALL genax(gen(1,1,ngen),3,4)
       IF (INDEX(tag,'OH').NE.0) THEN
          ngen=ngen+1
          CALL geni(gen(1,1,ngen))
       ENDIF
    ELSEIF (INDEX(tag(1:1),'I').NE.0) THEN
       pi=ACOS(-1._real_8)
       a=0.5_real_8
       b=COS(0.2_real_8*pi)
       c=COS(0.4_real_8*pi)
       ngen=ngen+1
       gen(1,1,ngen)=b
       gen(1,2,ngen)=-c
       gen(1,3,ngen)=-a
       gen(2,1,ngen)=-c
       gen(2,2,ngen)=a
       gen(2,3,ngen)=-b
       gen(3,1,ngen)=a
       gen(3,2,ngen)=b
       gen(3,3,ngen)=c
       ngen=ngen+1
       gen(1,1,ngen)=c
       gen(1,2,ngen)=-a
       gen(1,3,ngen)=b
       gen(2,1,ngen)=-a
       gen(2,2,ngen)=-b
       gen(2,3,ngen)=-c
       gen(3,1,ngen)=b
       gen(3,2,ngen)=-c
       gen(3,3,ngen)=-a
       IF (INDEX(tag,'IH').NE.0) THEN
          ngen=ngen+1
          CALL geni(gen(1,1,ngen))
       ENDIF
    ENDIF
    ! WRITE(6,*) ' POINT GROUP ',TAG,'   Naxis=',NAXIS
    ! WRITE(6,*) ' Number of generators ',NGEN
    ! SET UP SYMMETRY OPERATIONS
    ! START WITH IDENTITY
    nop=1
    DO i=1,3
       DO j=1,3
          tab(i,j,nop)=0.0_real_8
       ENDDO
       tab(i,i,nop)=1.0_real_8
    ENDDO
    ! NOW THE GENERATORS
    DO i=nop+1,nop+ngen
       CALL dcopy(9,gen(1,1,i-nop),1,tab(1,1,i),1)
    ENDDO
    nop=nop+ngen
    ! START THE GENERATION PROCESS
200 CONTINUE
    n1=nop
    DO i=1,n1
       DO j=1,n1
          CALL dgemm('N','N',3,3,3,1._real_8,tab(1,1,i),3,&
               tab(1,1,j),3,0._real_8,sop(1,1),3)
          DO k=1,nop
             CALL dcopy(9,sop(1,1),1,dsop(1,1),1)
             CALL daxpy(9,-1._real_8,tab(1,1,k),1,dsop(1,1),1)
             dd=dasum(9,dsop(1,1),1)
             IF (dd.LT.1.e-10_real_8) GOTO 100
          ENDDO
          nop=nop+1
          CALL dcopy(9,sop(1,1),1,tab(1,1,nop),1)
100       CONTINUE
       ENDDO
    ENDDO
    IF (n1.NE.nop) GOTO 200
    ! DO I=1,NOP
    ! WRITE(6,*) ' Symmetry operation ',I
    ! DO J=1,3
    ! WRITE(6,'(3F10.4)') (TAB(J,K,I),K=1,3)
    ! ENDDO
    ! ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE molsym
  ! GENERATORS
  ! ==================================================================
  SUBROUTINE genax(tab,n,k)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tab(3,3)
    INTEGER                                  :: n, k

    INTEGER                                  :: i, j
    REAL(real_8)                             :: a1, ang, ca, rot(2,2), &
                                                s1(3,3), s2(3,3), s3(3,3), &
                                                s4(3,3), sa, sq2, tpi

    tpi=4._real_8*dacos(0._real_8)
    ang=tpi/REAL(n,kind=real_8)
    DO i=1,3
       DO j=1,3
          tab(i,j)=0.0_real_8
       ENDDO
       tab(i,i)=1.0_real_8
    ENDDO
    rot(1,1)=COS(ang)
    rot(2,2)=COS(ang)
    rot(1,2)=SIN(ang)
    rot(2,1)=-SIN(ang)
    IF (k.EQ.1) THEN
       tab(2,2)=rot(1,1)
       tab(2,3)=rot(1,2)
       tab(3,2)=rot(2,1)
       tab(3,3)=rot(2,2)
    ELSEIF (k.EQ.2) THEN
       tab(1,1)=rot(1,1)
       tab(1,3)=rot(1,2)
       tab(3,1)=rot(2,1)
       tab(3,3)=rot(2,2)
    ELSEIF (k.EQ.3) THEN
       tab(1,1)=rot(1,1)
       tab(1,2)=rot(1,2)
       tab(2,1)=rot(2,1)
       tab(2,2)=rot(2,2)
    ELSEIF (k.EQ.4) THEN
       tab(1,1)=rot(1,1)
       tab(1,2)=rot(1,2)
       tab(2,1)=rot(2,1)
       tab(2,2)=rot(2,2)
       sq2=SQRT(2._real_8)*0.5_real_8
       a1=ATAN(SQRT(2._real_8))
       s1(1,1)=sq2
       s1(1,2)=-sq2
       s1(1,3)=0._real_8
       s1(2,1)=sq2
       s1(2,2)=sq2
       s1(2,3)=0._real_8
       s1(3,1)=0._real_8
       s1(3,2)=0._real_8
       s1(3,3)=1._real_8
       s2(1,1)=1._real_8
       s2(1,2)=0._real_8
       s2(1,3)=0._real_8
       s2(2,1)=0._real_8
       s2(2,2)=COS(a1)
       s2(2,3)=-SIN(a1)
       s2(3,1)=0._real_8
       s2(3,2)=SIN(a1)
       s2(3,3)=COS(a1)
       CALL dgemm('N','N',3,3,3,1._real_8,s1,3,s2,3,0._real_8,s3,3)
       CALL dgemm('N','N',3,3,3,1._real_8,tab,3,s3,3,0._real_8,s4,3)
       s1(1,2)=sq2
       s1(2,1)=-sq2
       s2(2,3)=SIN(a1)
       s2(3,2)=-SIN(a1)
       CALL dgemm('N','N',3,3,3,1._real_8,s1,3,s4,3,0._real_8,s3,3)
       CALL dgemm('N','N',3,3,3,1._real_8,s2,3,s3,3,0._real_8,tab,3)
    ELSEIF (k.EQ.5) THEN
       tab(1,1)=rot(1,1)
       tab(1,2)=rot(1,2)
       tab(2,1)=rot(2,1)
       tab(2,2)=rot(2,2)
       ca=SQRT(1._real_8-8._real_8/(3._real_8*(5._real_8+SQRT(5._real_8))))
       sa=4._real_8/SQRT(2._real_8*(5._real_8+SQRT(5._real_8)))
       s1(1,1)=1._real_8
       s1(1,2)=0._real_8
       s1(1,3)=0._real_8
       s1(2,1)=0._real_8
       s1(2,2)=ca
       s1(2,3)=-sa
       s1(3,1)=0._real_8
       s1(3,2)=sa
       s1(3,3)=ca
       ! SOMETHING IS MISSING
       CALL dgemm('N','N',3,3,3,1._real_8,tab,3,s1,3,0._real_8,s2,3)
       s1(3,2)=-sa
       s1(2,3)=sa
       CALL dgemm('N','N',3,3,3,1._real_8,s1,3,s2,3,0._real_8,tab,3)
    ELSE
       CALL stopgm('GENAX','K.GT.5',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE genax
  ! ==================================================================
  SUBROUTINE gensi(tab,k)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tab(3,3)
    INTEGER                                  :: k

    INTEGER                                  :: i, j

    DO i=1,3
       DO j=1,3
          tab(i,j)=0.0_real_8
       ENDDO
       tab(i,i)=1.0_real_8
    ENDDO
    tab(k,k)=-1.0_real_8
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gensi
  ! ==================================================================
  SUBROUTINE gensn(tab,n,k)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tab(3,3)
    INTEGER                                  :: n, k

! ==--------------------------------------------------------------==

    CALL genax(tab,n,k)
    tab(k,k)=-1.0_real_8
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gensn
  ! ==================================================================
  SUBROUTINE gend(tab,n,k)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tab(3,3)
    INTEGER                                  :: n, k

    INTEGER                                  :: i, j
    REAL(real_8)                             :: rot(3,3), sig(3,3)

    CALL genax(rot,2*n,k)
    DO i=1,3
       DO j=1,3
          sig(i,j)=0.0_real_8
       ENDDO
       sig(i,i)=1.0_real_8
    ENDDO
    IF (k.EQ.1) sig(2,2)=-1.0_real_8
    IF (k.EQ.2) sig(3,3)=-1.0_real_8
    IF (k.EQ.3) sig(1,1)=-1.0_real_8
    CALL dgemm('N','N',3,3,3,1._real_8,sig,3,rot,3,0._real_8,tab,3)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gend
  ! ==================================================================
  SUBROUTINE geni(tab)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tab(3,3)

    INTEGER                                  :: i, j

    DO i=1,3
       DO j=1,3
          tab(i,j)=0.0_real_8
       ENDDO
       tab(i,i)=-1.0_real_8
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE geni
  ! ==================================================================




END MODULE molsym_utils
