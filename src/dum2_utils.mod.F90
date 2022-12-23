MODULE dum2_utils
  USE cotr,                            ONLY: duat
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE mm_dimmod,                       ONLY: mmdim,&
                                             nat_grm
  USE parac,                           ONLY: paral
  USE rmas,                            ONLY: rmass
  USE system,                          ONLY: iatpt,&
                                             maxsys

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: dum2
  PUBLIC :: dumpr

CONTAINS

  ! ==================================================================
  SUBROUTINE dum2(tau0,tscr)
    ! ==--------------------------------------------------------------==
    ! == Put TAU0(3,maxsys%nax,maxsys%nsx) coordinates of atoms into                ==
    ! ==     TSCR(3,maxsys%nax*maxsys%nsx)                                          ==
    ! == Calculates Dummy atoms of type 2, type 3, and type 4         ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:), &
                                                tscr(3,maxsys%nax*maxsys%nsx)

    INTEGER                                  :: i, ia, iat, ib, id, is, n, nd
    REAL(real_8)                             :: t, x, y, z

    iat=0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          iat=iat+1
          tscr(1,iat)=tau0(1,ia,is)
          tscr(2,iat)=tau0(2,ia,is)
          tscr(3,iat)=tau0(3,ia,is)
       ENDDO
    ENDDO
    IF ((duat%ndat2.GT.0).OR.(duat%ndat3.GT.0).OR.(duat%ndat4.GT.0)) THEN
       DO i=1,duat%ndat
          ! calculate dummy atoms of type 2
          IF (duat%listda(i,1).EQ.2) THEN
             id=duat%listda(i,2)
             nd=duat%listd2(1,id)
             IF (nd.LT.0) THEN
                ! ND < 0 means all atoms (=IAT)
                nd=iat
             ELSE IF (nd.EQ.0) THEN
                IF (paral%io_parent)&
                     WRITE(6,'(1X,A,I5,A,I5,A,I5)')&
                     'DUM2! I=',I,' ID=',ID,' ND=',ND
                CALL stopgm('DUM2','INCORRECT NUMBER OF ATOMS',& 
                     __LINE__,__FILE__)
             ENDIF
             x=0.0_real_8
             y=0.0_real_8
             z=0.0_real_8
             DO n=1,nd
                IF (nd.EQ.iat) THEN
                   ib=n
                ELSE
                   ib=duat%listd2(n+1,id)
                ENDIF
                x=x+tscr(1,ib)
                y=y+tscr(2,ib)
                z=z+tscr(3,ib)
             ENDDO
             duat%dummy2(1,id)=x/REAL(nd,kind=real_8)
             duat%dummy2(2,id)=y/REAL(nd,kind=real_8)
             duat%dummy2(3,id)=z/REAL(nd,kind=real_8)
             ! calculate dummy atoms of type 3
          ELSE IF (duat%listda(i,1).EQ.3) THEN
             id=duat%listda(i,2)
             nd=duat%listd3(1,id)
             IF (nd.LT.0) THEN
                ! ND < 0 means all atoms (=IAT)
                nd=iat
             ELSE IF (nd.EQ.0) THEN
                IF (paral%io_parent)&
                     WRITE(6,'(1X,A,I5,A,I5,A,I5)')&
                     'DUM2! I=',I,' ID=',ID,' ND=',ND
                CALL stopgm('DUM2','INCORRECT NUMBER OF ATOMS',& 
                     __LINE__,__FILE__)
             ENDIF
             x=0.0_real_8
             y=0.0_real_8
             z=0.0_real_8
             t=0.0_real_8
             DO n=1,nd
                IF (nd.EQ.iat) THEN
                   ib=n
                ELSE
                   ib=duat%listd3(n+1,id)
                ENDIF
                is=iatpt(2,ib)
                x=x+tscr(1,ib)*rmass%pma0(is)
                y=y+tscr(2,ib)*rmass%pma0(is)
                z=z+tscr(3,ib)*rmass%pma0(is)
                t=t+rmass%pma0(is)
             ENDDO
             IF (t.GT.0.0_real_8) THEN
                duat%dummy3(1,id)=x/t
                duat%dummy3(2,id)=y/t
                duat%dummy3(3,id)=z/t
             ELSE
                IF (paral%io_parent)&
                     WRITE(6,'(1X,A,I5,A,I5,A,I5)')&
                     'DUM2! I=',I,' ID=',ID,' ND=',ND
                CALL stopgm('DUM2','INCORRECT WEIGHTS',& 
                     __LINE__,__FILE__)
             ENDIF
          ELSE IF (duat%listda(i,1).EQ.4) THEN
             id=duat%listda(i,2)
             nd=duat%listd4(1,id)
             x=0.0_real_8
             y=0.0_real_8
             z=0.0_real_8
             t=0.0_real_8
             DO n=1,nd
                x=x+tscr(1,duat%listd4(n+1,id))*duat%weigd4(n,id)
                y=y+tscr(2,duat%listd4(n+1,id))*duat%weigd4(n,id)
                z=z+tscr(3,duat%listd4(n+1,id))*duat%weigd4(n,id)
                t=t+duat%weigd4(n,id)
             ENDDO
             IF (t.GT.0.0_real_8) THEN
                duat%dummy4(1,id)=x/t
                duat%dummy4(2,id)=y/t
                duat%dummy4(3,id)=z/t
             ELSE
                IF (paral%io_parent)&
                     WRITE(6,'(1X,A,I5,A,I5,A,I5)')&
                     'DUM2! I=',I,' ID=',ID,' ND=',ND
                CALL stopgm('DUM2','INCORRECT WEIGHTS',& 
                     __LINE__,__FILE__)
             ENDIF
          ENDIF
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE dum2
  ! ==================================================================
  SUBROUTINE dumpr
    ! ==--------------------------------------------------------------==
    ! Variables
    INTEGER                                  :: i, id, ityp, k

! ==--------------------------------------------------------------==

    IF (duat%ndat.EQ.0) RETURN
    IF (.NOT.paral%io_parent) RETURN

    WRITE(6,'(/,25X,A)')' DUMMY ATOMS '
    WRITE(6,'(A)')' #DUM   INDEX                      POSITION (a.u.)'
    DO i=1,duat%ndat
       ityp=duat%listda(i,1)
       id=duat%listda(i,2)
       IF (ityp.EQ.1) WRITE(6,'(I5,I8,A,3F15.5)')&
            i,nat_grm(mmdim%natm+i),'  TYPE1',(duat%dummy1(k,id),k=1,3)
       IF (ityp.EQ.2) WRITE(6,'(I5,I8,A,3F15.5)')&
            i,nat_grm(mmdim%natm+i),'  TYPE2',(duat%dummy2(k,id),k=1,3)
       IF (ityp.EQ.3) WRITE(6,'(I5,I8,A,3F15.5)')&
            i,nat_grm(mmdim%natm+i),'  TYPE3',(duat%dummy3(k,id),k=1,3)
       IF (ityp.EQ.4) WRITE(6,'(I5,I8,A,3F15.5)')&
            i,nat_grm(mmdim%natm+i),'  TYPE4',(duat%dummy4(k,id),k=1,3)
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE dumpr
  ! ==================================================================

END MODULE dum2_utils
