MODULE purge_utils
  USE cotr,                            ONLY: cnsval,&
                                             cotc0,&
                                             ntcnst
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_bcast
  USE parac,                           ONLY: parai,&
                                             paral
  USE utils,                           ONLY: invmat
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: purge
  !public :: getco
  !public :: putco

CONTAINS

  ! ==================================================================
  SUBROUTINE purge(tau0,fion)
    ! ==--------------------------------------------------------------==
    ! ==  PURGE FORCES ALONG CONSTRAINT COORDINATES                   ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:), fion(:,:,:)

    INTEGER                                  :: i, ia, ib, icon, info, iset, &
                                                ityp
    LOGICAL                                  :: tconst
    REAL(real_8)                             :: aux(36), binv(6,6), &
                                                bmat(6,6), dx(3), f1, f2, &
                                                fo(6), phi, rab, so(6), &
                                                theta, xa(3), xb(3)

    IF (cotc0%mcnstr.LT.1) RETURN
    IF (cotc0%mcnstr.GT.1) CALL stopgm('PURGE','NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    IF (paral%parent) THEN
       icon=1
       ityp=ntcnst(1,icon)
       IF (ityp.EQ.1.OR.ityp.EQ.4) THEN
          ia=ntcnst(2,icon)
          CALL getco(ia,tau0,xa)
          CALL getco(ia,fion,fo)
          ib=ntcnst(3,icon)
          CALL getco(ib,tau0,xb)
          CALL getco(ib,fion,fo(4))
          CALL zeroing(bmat)!,36)
          dx(1)=xb(1)-xa(1)
          dx(2)=xb(2)-xa(2)
          dx(3)=xb(3)-xa(3)
          rab=SQRT(dx(1)*dx(1)+dx(2)*dx(2)+dx(3)*dx(3))
          theta=ACOS(dx(3)/rab)
          phi=ASIN(dx(2)/SQRT(dx(1)*dx(1)+dx(2)*dx(2)))
          IF (ABS(COS(phi)).LT.1.e-10_real_8) THEN
             phi=ACOS(dx(1)/SQRT(dx(1)*dx(1)+dx(2)*dx(2)))
             iset=2
          ELSE
             iset=1
          ENDIF
          IF (ityp.EQ.1) cnsval(icon)=rab*rab
          IF (ityp.EQ.4) cnsval(icon)=rab
          ! ..center of bond
          bmat(1,1)=0.5_real_8
          bmat(1,4)=0.5_real_8
          bmat(2,2)=0.5_real_8
          bmat(2,5)=0.5_real_8
          bmat(3,3)=0.5_real_8
          bmat(3,6)=0.5_real_8
          ! ..rab
          bmat(4,1)=-dx(1)/rab
          bmat(4,2)=-dx(2)/rab
          bmat(4,3)=-dx(3)/rab
          bmat(4,4)= dx(1)/rab
          bmat(4,5)= dx(2)/rab
          bmat(4,6)= dx(3)/rab
          ! ..Theta
          f1=dx(3)/rab**3/SIN(theta)
          f2=1._real_8/rab/SIN(theta)
          bmat(5,1)=-dx(1)*f1
          bmat(5,2)=-dx(2)*f1
          bmat(5,3)=-dx(3)*f1+f2
          bmat(5,4)= dx(1)*f1
          bmat(5,5)= dx(2)*f1
          bmat(5,6)= dx(3)*f1-f2
          ! ..phi
          IF (iset.EQ.1) THEN
             f1=dx(2)/COS(phi)/(dx(1)*dx(1)+dx(2)*dx(2))**1.5_real_8
             f2=1._real_8/COS(phi)/SQRT(dx(1)*dx(1)+dx(2)*dx(2))
             bmat(6,1)= dx(1)*f1
             bmat(6,2)= dx(2)*f1-f2
             bmat(6,4)=-dx(1)*f1
             bmat(6,5)=-dx(2)*f1+f2
          ELSE
             f1=-dx(1)/SIN(phi)/(dx(1)*dx(1)+dx(2)*dx(2))**1.5_real_8
             f2=-1._real_8/SIN(phi)/SQRT(dx(1)*dx(1)+dx(2)*dx(2))
             bmat(6,1)= dx(1)*f1-f2
             bmat(6,2)= dx(2)*f1
             bmat(6,4)=-dx(1)*f1+f2
             bmat(6,5)=-dx(2)*f1
          ENDIF
          ! 
          IF (cotc0%lshove) THEN
             i=ntcnst(6,icon)
             IF (i.EQ.0) THEN
                tconst=.TRUE.
             ELSE
                tconst=.FALSE.
                CALL dgemv('N',6,6,1._real_8,bmat,6,fo,1,0._real_8,so,1)
                IF (i.GT.0.AND.so(4).LT.0._real_8) tconst=.TRUE.
                IF (i.LT.0.AND.so(4).GT.0._real_8) tconst=.TRUE.
             ENDIF
          ELSE
             tconst=.TRUE.
          ENDIF
          IF (tconst) THEN
             CALL dcopy(36,bmat,1,binv,1)
             CALL invmat(6,binv,aux,info)
             IF (info.NE.0) CALL stopgm('PURGE','BINV INVALID',& 
                  __LINE__,__FILE__)
             DO i=1,10
                CALL dgemv('N',6,6,1._real_8,bmat,6,fo,1,0._real_8,so,1)
                IF (ABS(so(4)).LT.1.e-12_real_8) GOTO 10
                so(4)=0._real_8
                CALL dgemv('N',6,6,1._real_8,binv,6,so,1,0._real_8,fo,1)
             ENDDO
             CALL stopgm('PURGE','ITERATION DID NOT CONVERGE',& 
                  __LINE__,__FILE__)
10           CONTINUE
             CALL putco(ia,fion,fo(1))
             CALL putco(ib,fion,fo(4))
          ENDIF
       ELSE
          CALL stopgm('PURGE','NOT IMPLEMENTED',& 
               __LINE__,__FILE__)
       ENDIF
    ENDIF
    CALL mp_bcast(fion,SIZE(fion),parai%source,parai%allgrp)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE purge
  ! ==================================================================
  SUBROUTINE getco(iat,tau,t)
    INTEGER                                  :: iat
    REAL(real_8)                             :: tau(:,:,:), t(3)

    INTEGER                                  :: i0, ia, is

! ==--------------------------------------------------------------==

    i0=0
    DO is=1,ions1%nsp
       IF (i0.LT.iat.AND.i0+ions0%na(is).GE.iat) THEN
          ia=iat-i0
          t(1)=tau(1,ia,is)
          t(2)=tau(2,ia,is)
          t(3)=tau(3,ia,is)
          GOTO 100
       ELSE
          i0=i0+ions0%na(is)
       ENDIF
    ENDDO
    CALL stopgm('GETCO','ATOM NUMBER OUT OF RANGE',& 
         __LINE__,__FILE__)
100 CONTINUE
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE getco
  ! ==================================================================
  SUBROUTINE putco(iat,tau,t)
    INTEGER                                  :: iat
    REAL(real_8)                             :: tau(:,:,:), t(3)

    INTEGER                                  :: i0, ia, is

! ==--------------------------------------------------------------==

    i0=0
    DO is=1,ions1%nsp
       IF (i0.LT.iat.AND.i0+ions0%na(is).GE.iat) THEN
          ia=iat-i0
          tau(1,ia,is)=t(1)
          tau(2,ia,is)=t(2)
          tau(3,ia,is)=t(3)
          GOTO 100
       ELSE
          i0=i0+ions0%na(is)
       ENDIF
    ENDDO
    CALL stopgm('PUTCO','ATOM NUMBER OUT OF RANGE',& 
         __LINE__,__FILE__)
100 CONTINUE
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE putco
  ! ==================================================================

END MODULE purge_utils
