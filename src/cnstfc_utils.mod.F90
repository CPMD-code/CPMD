MODULE cnstfc_utils
  USE constr_utils,                    ONLY: &
       diffd, diffdd, diffo, diffr, difft, funcd, funcdd, funco, funcp, &
       funcr, funct
  USE cotr,                            ONLY: &
       anorm, askel, c_kappa, c_rc, cnpar, cnsval, cotc0, cotr007, csigm, &
       duat, fc, fcstr, fv, lskptr, ntcnst, ntrest, resc, resfdiff, resfor, &
       resm, respar, respos, resv, resval, rskel, rsmass, xlagr
  USE dum2_utils,                      ONLY: dum2
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE fillc_utils,                     ONLY: fillc
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE isos,                            ONLY: isos1
  USE kinds,                           ONLY: real_8
  USE meta_cv_utils,                   ONLY: coorntot_rf
  USE parac,                           ONLY: paral
  USE pbc_utils,                       ONLY: pbc
  USE puttau_utils,                    ONLY: gettau
  USE rmas,                            ONLY: rmass
  USE system,                          ONLY: iatpt,&
                                             maxsys,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: cnstfc
  PUBLIC :: coornum
  PUBLIC :: coornumsp
  PUBLIC :: coornumgrp
  PUBLIC :: coorn_rf
  PUBLIC :: bndswitch
  PUBLIC :: restfc

CONTAINS

  ! ==================================================================
  SUBROUTINE cnstfc(tau0,tscr,dxpar,cnmax,lagra)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:), tscr(3,*), &
                                                dxpar(*), cnmax
    LOGICAL                                  :: lagra

    INTEGER                                  :: dum(1), i, ia, ib, ic, id, &
                                                idamax, isub, ityp, ix, j, k, &
                                                kmax, n
    REAL(real_8)                             :: cval, dx(12), r0_shift, &
                                                sign0, x1(3), x2(3), x3(3), &
                                                x4(3)

    cnmax=0.0_real_8
    ener_com%ecnstr=0.0_real_8
    IF (cotc0%mcnstr.EQ.0) RETURN
    CALL tiset('    CNSTFC',isub)
    ! No constraints -> return.

    CALL dum2(tau0,tscr)
    CALL zeroing(anorm)!,cotc0%mcnstr*cotc0%nodim)
    DO i=1,cotc0%mcnstr
       cval=cnsval(i)
       ityp=ntcnst(1,i)
       ia=ntcnst(2,i)
       ib=ntcnst(3,i)
       ic=ntcnst(4,i)
       id=ntcnst(5,i)
       IF (ityp.EQ.1) THEN
          ! ..stretch
          CALL fillc(ia,tscr,x1)
          CALL fillc(ib,tscr,x2)
          CALL funcr(fc(i),fv(i),cval,x1,x2)
          CALL diffr(dx,x1,x2)
          kmax=6
          ! ..angle
       ELSEIF (ityp.EQ.2) THEN
          CALL fillc(ia,tscr,x1)
          CALL fillc(ib,tscr,x2)
          CALL fillc(ic,tscr,x3)
          CALL funct(fc(i),fv(i),cval,x1,x2,x3)
          CALL difft(dx,x1,x2,x3)
          kmax=9
       ELSEIF (ityp.EQ.3) THEN
          ! ..dihedral angle
          CALL fillc(ia,tscr,x1)
          CALL fillc(ib,tscr,x2)
          CALL fillc(ic,tscr,x3)
          CALL fillc(id,tscr,x4)
          CALL funco(fc(i),fv(i),cval,x1,x2,x3,x4,sign0)
          CALL diffo(dx,x1,x2,x3,x4,sign0)
          kmax=12
       ELSEIF (ityp.EQ.4) THEN
          ! ..distance
          CALL fillc(ia,tscr,x1)
          CALL fillc(ib,tscr,x2)
          CALL funcd(fc(i),fv(i),cval,x1,x2)
          CALL diffd(dx,x1,x2)
          kmax=6
       ELSEIF (ityp.EQ.5) THEN
          ! ..out of plane
          CALL fillc(ia,tscr,x1)
          CALL fillc(ib,tscr,x2)
          CALL fillc(ic,tscr,x3)
          CALL fillc(id,tscr,x4)
          CALL funcp(fc(i),fv(i),cval,x1,x2,x3,x4,dx)
          kmax=12
       ELSEIF (ityp.EQ.6) THEN
          ! ..coordination number
          c_kappa = cnpar(1,i)
          c_rc    = cnpar(2,i)
          CALL coornum(ia,c_kappa,c_rc,tscr,ions1%nat,anorm(1,i),lskptr,&
               fc(i),fv(i),cval)
          kmax=0
       ELSEIF (ityp.EQ.7) THEN
          ! ..difference between distances- 3 atoms
          CALL fillc(ia,tscr,x1)
          CALL fillc(ib,tscr,x2)
          CALL fillc(ic,tscr,x3)
          CALL funcdd(fc(i),fv(i),cval,x1,x2,x3)
          CALL diffdd(dx,x1,x2,x3)
          kmax=9
       ELSEIF (ityp.EQ.8) THEN
          ! ..species dependent coordination number with Fermi function
          c_kappa = cnpar(1,i)
          c_rc    = cnpar(2,i)
          CALL coornumsp(ia,ib,c_kappa,c_rc,tscr,anorm(1,i),&
               lskptr,fc(i),fv(i),cval)
          kmax=0
       ELSEIF (ityp.EQ.9) THEN
          ! ..species dependent coordination number with rational function
          c_rc    = cnpar(1,i)
          r0_shift  = cnpar(2,i)
          CALL coorn_rf(ia,ib,ic,id,c_rc,r0_shift,tscr,anorm(1,i),&
               fc(i),fv(i),cval,1,dum)
          kmax=0
       ELSEIF (ityp.EQ.10) THEN
          ! ..bond switch with rational function
          c_rc    = cnpar(1,i)
          CALL bndswitch(ia,ib,ic,id,c_rc,tscr,anorm(1,i),&
               fc(i),fv(i),cval)
          kmax=0
       ELSEIF (ityp.EQ.11) THEN
          ! ..total species dependent coordination number with rational function
          c_rc    = cnpar(1,i)
          r0_shift  = cnpar(2,i)
          CALL coorntot_rf(ia,ib,ic,id,c_rc,r0_shift,&
               tscr,anorm(1,i),fc(i),fv(i),cval,1,dum)
          kmax=0
       ELSEIF (ityp.EQ.12) THEN
          ! ..distance between two atoms along a selected direction
          CALL fillc(ia,tscr,x1)
          CALL fillc(ib,tscr,x2)
          DO k=1,3
             IF (k.NE.ic) THEN
                x1(k)=0.0_real_8
                x2(k)=0.0_real_8
             ENDIF
          ENDDO
          CALL funcd(fc(i),fv(i),cval,x1,x2)
          CALL diffd(dx,x1,x2)
          kmax=6
       ELSE
          WRITE(6,*) ' CNSTFC! I=',I,' TYPE=',ITYP
          CALL stopgm('CNSTFC','UNKNOWN TYPE OF CONSTRAINT',& 
               __LINE__,__FILE__)
       ENDIF

       DO k=1,kmax
          DO n=1,cotc0%nodim
             anorm(n,i)=anorm(n,i)+dx(k)*askel(n,i,k)
          ENDDO
       ENDDO
    ENDDO
    CALL zeroing(fcstr)!,cotc0%nodim)
    DO i=1,cotc0%mcnstr
       DO j=1,cotc0%nodim
          IF (lagra) THEN
             fcstr(j)=fcstr(j)+xlagr(i)*anorm(j,i)
          ELSE
             fcstr(j)=fcstr(j)+csigm(i)*anorm(j,i)*fc(i)
          ENDIF
       ENDDO
       IF (.NOT.lagra) ener_com%ecnstr=ener_com%ecnstr+0.5_real_8*csigm(i)*fc(i)*fc(i)
    ENDDO
    ix=idamax(cotc0%nodim,fcstr(1),1)
    cnmax=ABS(fcstr(ix))
    CALL daxpy(cotc0%nodim,1.0_real_8,fcstr(1),1,dxpar(1),1)
    CALL tihalt('    CNSTFC',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE cnstfc
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE coornum(ia,c1_kappa,c1_rc,tscr,nat,an,lsk,fci,fvi,cval)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ia
    REAL(real_8)                             :: c1_kappa, c1_rc, tscr(3,*)
    INTEGER                                  :: nat
    REAL(real_8)                             :: an(*)
    INTEGER                                  :: lsk(3,nat)
    REAL(real_8)                             :: fci, fvi, cval

    INTEGER                                  :: iat, ib, is, ityp, j, jj, jx, &
                                                jy, jz, k, kx, l1, l2, l3, naa
    REAL(real_8)                             :: dd, df, dx, dy, dz, ff, fff, &
                                                x0, xx(3), y0, z0, dx_, dy_, dz_

    IF (ia.LE.nat) THEN
       l1 = lsk(1,ia)
       l2 = lsk(2,ia)
       l3 = lsk(3,ia)
    ENDIF
    CALL fillc(ia,tscr,xx)
    x0 = xx(1)
    y0 = xx(2)
    z0 = xx(3)
    fvi = 0._real_8
    DO iat=1,nat
       dx_=tscr(1,iat)-x0
       dy_=tscr(2,iat)-y0
       dz_=tscr(3,iat)-z0
       IF (.NOT. isos1%tisos) THEN
          CALL pbc(dx_,dy_,dz_,dx,dy,dz,1,parm%apbc,parm%ibrav)
       ELSE
          dx=dx_
          dy=dy_
          dz=dz_
       ENDIF
       dd=SQRT(dx*dx+dy*dy+dz*dz)
       ! ..exclude self interaction
       IF (dd.GT.1.e-2_real_8) THEN
          df=c1_kappa*(dd-c1_rc)
          fvi=fvi+1._real_8/(EXP(df)+1._real_8)
          ff=-0.5_real_8*c1_kappa/(COSH(df)+1._real_8)/dd
          IF (lsk(1,iat).NE.0) THEN
             k=lsk(1,iat)
             an(k)=an(k)+ff*dx
          ENDIF
          IF (lsk(2,iat).NE.0) THEN
             k=lsk(2,iat)
             an(k)=an(k)+ff*dy
          ENDIF
          IF (lsk(3,iat).NE.0) THEN
             k=lsk(3,iat)
             an(k)=an(k)+ff*dz
          ENDIF
          IF (ia.LE.nat) THEN
             IF (l1.NE.0) an(l1)=an(l1)-ff*dx
             IF (l2.NE.0) an(l2)=an(l2)-ff*dy
             IF (l3.NE.0) an(l3)=an(l3)-ff*dz
          ELSEIF (ia.LE.nat+duat%ndat) THEN
             naa=ia-nat
             ityp=duat%listda(naa,1)
             IF (ityp.EQ.1) THEN
                GOTO 110
             ELSEIF (ityp.EQ.2) THEN
                kx=duat%listda(naa,2)
                jj=duat%listd2(1,kx)
                IF (jj.LT.0) THEN
                   jj=nat
                ELSE IF (jj.EQ.0) THEN
                   IF (paral%io_parent)&
                        WRITE(6,'(1X,A,I5,A,I5,A,I5,A,I5)')&
                        'COORNUM! IA=',ia,' ITYP=',ityp,' KX=',kx,' JJ=',jJ
                   CALL stopgm('COORNUM','INCORRECT NUMBERS FOR IA',& 
                        __LINE__,__FILE__)
                ENDIF
                fff=1.0_real_8/REAL(jj,kind=real_8)
                DO j=1,jj
                   IF (jj.EQ.nat) THEN
                      ib=j
                   ELSE
                      ib=duat%listd2(j+1,kx)
                   ENDIF
                   jx=lskptr(1,ib)
                   jy=lskptr(2,ib)
                   jz=lskptr(3,ib)
                   IF (jx.NE.0) an(jx)=an(jx)-ff*dx*fff
                   IF (jy.NE.0) an(jy)=an(jy)-ff*dy*fff
                   IF (jz.NE.0) an(jz)=an(jz)-ff*dz*fff
                ENDDO
             ELSEIF (ityp.EQ.3) THEN
                kx=duat%listda(naa,2)
                jj=duat%listd3(1,kx)
                IF (jj.LT.0) THEN
                   jj=nat
                ELSE IF (jj.EQ.0) THEN
                   IF (paral%io_parent)&
                        WRITE(6,'(1X,A,I5,A,I5,A,I5,A,I5)')&
                        'COORNUM! IA=',ia,' ITYP=',ityp,' KX=',kx,' JJ=',jJ
                   CALL stopgm('COORNUM','INCORRECT NUMBERS FOR IA',& 
                        __LINE__,__FILE__)
                ENDIF
                fff=0.0_real_8
                DO j=1,jj
                   IF (jj.EQ.nat) THEN
                      ib=j
                   ELSE
                      ib=duat%listd3(j+1,kx)
                   ENDIF
                   is=iatpt(2,ib)
                   fff=fff+rmass%pma0(is)
                ENDDO
                fff=1.0_real_8/fff
                DO j=1,jj
                   IF (jj.EQ.nat) THEN
                      ib=j
                   ELSE
                      ib=duat%listd3(j+1,kx)
                   ENDIF
                   is=iatpt(2,ib)
                   jx=lskptr(1,ib)
                   jy=lskptr(2,ib)
                   jz=lskptr(3,ib)
                   IF (jx.NE.0) an(jx)=an(jx)-ff*dx*fff*rmass%pma0(is)
                   IF (jy.NE.0) an(jy)=an(jy)-ff*dy*fff*rmass%pma0(is)
                   IF (jz.NE.0) an(jz)=an(jz)-ff*dz*fff*rmass%pma0(is)
                ENDDO
             ELSEIF (ityp.EQ.4) THEN
                kx=duat%listda(naa,2)
                jj=duat%listd4(1,kx)
                fff=0.0_real_8
                DO j=1,jj
                   fff=fff+duat%weigd4(j,kx)
                ENDDO
                fff=1.0_real_8/fff
                DO j=1,jj
                   ib=duat%listd4(j+1,kx)
                   jx=lskptr(1,ib)
                   jy=lskptr(2,ib)
                   jz=lskptr(3,ib)
                   IF (jx.NE.0) an(jx)=an(jx)-ff*dx*fff*duat%weigd4(j,kx)
                   IF (jy.NE.0) an(jy)=an(jy)-ff*dy*fff*duat%weigd4(j,kx)
                   IF (jz.NE.0) an(jz)=an(jz)-ff*dz*fff*duat%weigd4(j,kx)
                ENDDO
             ENDIF
          ENDIF
110       CONTINUE
       ENDIF
    ENDDO
    fci=fvi-cval
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE coornum
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE coornumsp(ia,isp,c1_kappa,c1_rc,tscr,an,lsk,&
       fci,fvi,cval)
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: ia, isp
    REAL(real_8)                             :: c1_kappa, c1_rc, tscr(3,*), &
                                                an(*)
    INTEGER                                  :: lsk(3,ions1%nat)
    REAL(real_8)                             :: fci, fvi, cval

    INTEGER                                  :: iat, ib, iend, iini, is, &
                                                ityp, j, jj, jx, jy, jz, k, &
                                                kx, l1, l2, l3, naa
    REAL(real_8)                             :: dd, df, dx, dy, dz, ff, fff, &
                                                x0, xx(3), y0, z0, dx_, dy_, dz_

    IF (ia.LE.ions1%nat) THEN
       l1 = lsk(1,ia)
       l2 = lsk(2,ia)
       l3 = lsk(3,ia)
    ENDIF
    CALL fillc(ia,tscr,xx)
    x0 = xx(1)
    y0 = xx(2)
    z0 = xx(3)

    fvi = 0._real_8
    iini = 1
    DO is = 1,isp-1
       iini=iini+ions0%na(is)
    ENDDO
    iend = iini + ions0%na(isp) -1
    DO iat=iini,iend
       dx_=tscr(1,iat)-x0
       dy_=tscr(2,iat)-y0
       dz_=tscr(3,iat)-z0
       IF (.NOT. isos1%tisos)THEN
          CALL pbc(dx_,dy_,dz_,dx,dy,dz,1,parm%apbc,parm%ibrav)
       ELSE
          dx=dx_
          dy=dy_
          dz=dz_
       ENDIF
       dd=SQRT(dx*dx+dy*dy+dz*dz)
       ! ..exclude self interaction
       IF (dd.GT.1.e-2_real_8) THEN
          df=c1_kappa*(dd-c1_rc)
          fvi=fvi+1._real_8/(EXP(df)+1._real_8)
          ff=-0.5_real_8*c1_kappa/(COSH(df)+1._real_8)/dd
          IF (lsk(1,iat).NE.0) THEN
             k=lsk(1,iat)
             an(k)=an(k)+ff*dx
          ENDIF
          IF (lsk(2,iat).NE.0) THEN
             k=lsk(2,iat)
             an(k)=an(k)+ff*dy
          ENDIF
          IF (lsk(3,iat).NE.0) THEN
             k=lsk(3,iat)
             an(k)=an(k)+ff*dz
          ENDIF
          IF (ia.LE.ions1%nat) THEN
             IF (l1.NE.0) an(l1)=an(l1)-ff*dx
             IF (l2.NE.0) an(l2)=an(l2)-ff*dy
             IF (l3.NE.0) an(l3)=an(l3)-ff*dz
          ELSEIF (ia.LE.ions1%nat+duat%ndat) THEN
             naa=ia-ions1%nat
             ityp=duat%listda(naa,1)
             IF (ityp.EQ.1) THEN
                GOTO 110
             ELSEIF (ityp.EQ.2) THEN
                kx=duat%listda(naa,2)
                jj=duat%listd2(1,kx)
                IF (jj.LT.0) THEN
                   jj=ions1%nat
                ELSE IF (jj.EQ.0) THEN
                   IF (paral%io_parent)&
                        WRITE(6,'(1X,A,I5,A,I5,A,I5,A,I1)')&
                        'COORNUMSP! IA=',ia,' ITYP=',ityp,' KX=',kx,' JJ=',jJ
                   CALL stopgm('COORNUMSP','INCORRECT NUMBERS FOR IA',& 
                        __LINE__,__FILE__)
                ENDIF
                fff=1.0_real_8/REAL(jj,kind=real_8)
                DO j=1,jj
                   IF (jj.EQ.ions1%nat) THEN
                      ib=j
                   ELSE
                      ib=duat%listd2(j+1,kx)
                   ENDIF
                   jx=lskptr(1,ib)
                   jy=lskptr(2,ib)
                   jz=lskptr(3,ib)
                   IF (jx.NE.0) an(jx)=an(jx)-ff*dx*fff
                   IF (jy.NE.0) an(jy)=an(jy)-ff*dy*fff
                   IF (jz.NE.0) an(jz)=an(jz)-ff*dz*fff
                ENDDO
             ELSEIF (ityp.EQ.3) THEN
                kx=duat%listda(naa,2)
                jj=duat%listd3(1,kx)
                IF (jj.LT.0) THEN
                   jj=ions1%nat
                ELSE IF (jj.EQ.0) THEN
                   IF (paral%io_parent)&
                        WRITE(6,'(1X,A,I5,A,I5,A,I5,A,I5)')&
                        'COORNUMSP! IA=',ia,' ITYP=',ityp,' KX=',kx,' JJ=',jJ
                   CALL stopgm('COORNUMSP','INCORRECT NUMBERS FOR IA',& 
                        __LINE__,__FILE__)
                ENDIF
                fff=0.0_real_8
                DO j=1,jj
                   IF (jj.EQ.ions1%nat) THEN
                      ib=j
                   ELSE
                      ib=duat%listd3(j+1,kx)
                   ENDIF
                   is=iatpt(2,ib)
                   fff=fff+rmass%pma0(is)
                ENDDO
                fff=1.0_real_8/fff
                DO j=1,jj
                   IF (jj.EQ.ions1%nat) THEN
                      ib=j
                   ELSE
                      ib=duat%listd3(j+1,kx)
                   ENDIF
                   is=iatpt(2,ib)
                   jx=lskptr(1,ib)
                   jy=lskptr(2,ib)
                   jz=lskptr(3,ib)
                   IF (jx.NE.0) an(jx)=an(jx)-ff*dx*fff*rmass%pma0(is)
                   IF (jy.NE.0) an(jy)=an(jy)-ff*dy*fff*rmass%pma0(is)
                   IF (jz.NE.0) an(jz)=an(jz)-ff*dz*fff*rmass%pma0(is)
                ENDDO
             ELSEIF (ityp.EQ.4) THEN
                kx=duat%listda(naa,2)
                jj=duat%listd4(1,kx)
                fff=0.0_real_8
                DO j=1,jj
                   fff=fff+duat%weigd4(j,kx)
                ENDDO
                fff=1.0_real_8/fff
                DO j=1,jj
                   ib=duat%listd4(j+1,kx)
                   jx=lskptr(1,ib)
                   jy=lskptr(2,ib)
                   jz=lskptr(3,ib)
                   IF (jx.NE.0) an(jx)=an(jx)-ff*dx*fff*duat%weigd4(j,kx)
                   IF (jy.NE.0) an(jy)=an(jy)-ff*dy*fff*duat%weigd4(j,kx)
                   IF (jz.NE.0) an(jz)=an(jz)-ff*dz*fff*duat%weigd4(j,kx)
                ENDDO
             ENDIF
          ENDIF
110       CONTINUE
       ENDIF
    ENDDO
    fci=fvi-cval
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE coornumsp
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE coornumgrp(na_grp,maxnb_grp,nb_grp,&
       IA,IB,C1_KAPPA,C1_RC_A,TSCR,AN,LSK,FCI,FVI,CVAL)
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: na_grp, maxnb_grp, &
                                                nb_grp(na_grp), ia(na_grp), &
                                                ib(na_grp*maxnb_grp)
    REAL(real_8)                             :: c1_kappa, c1_rc_a(na_grp), &
                                                tscr(3,*), an(*)
    INTEGER                                  :: LSK(3,ions1%nat)
    REAL(real_8)                             :: fci, fvi, cval

    INTEGER                                  :: atma, atmb, iatma, iatmb, &
                                                ibb, ipgrpb, is, ityp, j, jj, &
                                                jx, jy, jz, k, kx, l1, l2, &
                                                l3, naa
    REAL(real_8)                             :: dd, df, dx, dy, dz, ff, fff, &
                                                x0, xx(3), y0, z0, dx_, dy_, dz_

    fvi = 0._real_8
    ipgrpb=0
    DO iatma=1,na_grp
       atma=ia(iatma)
       IF (atma.LE.ions1%nat) THEN
          l1 = lsk(1,atma)! IAx
          l2 = lsk(2,atma)! IAy
          l3 = lsk(3,atma)! IAz
       ENDIF
       CALL fillc(atma,tscr,xx)
       x0 = xx(1)
       y0 = xx(2)
       z0 = xx(3)
       DO iatmb=1,nb_grp(iatma)
          ipgrpb=ipgrpb+1
          atmb=ib(ipgrpb)
          dx_=tscr(1,atmb)-x0
          dy_=tscr(2,atmb)-y0
          dz_=tscr(3,atmb)-z0
          IF (.NOT. isos1%tisos) THEN
             CALL pbc(dx_,dy_,dz_,dx,dy,dz,1,parm%apbc,parm%ibrav)
          ELSE
             dx=dx_
             dy=dy_
             dz=dz_
          ENDIF
          dd=SQRT(dx*dx+dy*dy+dz*dz)
          ! ..exclude self interaction
          IF (dd.GT.1.e-2_real_8) THEN
             df=c1_kappa*(dd-c1_rc_a(iatma))
             fvi=fvi+1._real_8/(EXP(df)+1._real_8)
             ff=-0.5_real_8*c1_kappa/(COSH(df)+1._real_8)/dd
             IF (lsk(1,atmb).NE.0) THEN
                k=lsk(1,atmb)
                an(k)=an(k)+ff*dx
             ENDIF
             IF (lsk(2,atmb).NE.0) THEN
                k=lsk(2,atmb)
                an(k)=an(k)+ff*dy
             ENDIF
             IF (lsk(3,atmb).NE.0) THEN
                k=lsk(3,atmb)
                an(k)=an(k)+ff*dz
             ENDIF
             IF (atma.LE.ions1%nat) THEN
                IF (l1.NE.0) an(l1)=an(l1)-ff*dx
                IF (l2.NE.0) an(l2)=an(l2)-ff*dy
                IF (l3.NE.0) an(l3)=an(l3)-ff*dz
             ELSEIF (atma.LE.ions1%nat+duat%ndat) THEN
                naa=atma-ions1%nat
                ityp=duat%listda(naa,1)
                IF (ityp.EQ.1) THEN
                   GOTO 130
                ELSEIF (ityp.EQ.2) THEN
                   kx=duat%listda(naa,2)
                   jj=duat%listd2(1,kx)
                   IF (jj.LT.0) THEN
                      jj=ions1%nat
                   ELSE IF (jj.EQ.0) THEN
                      IF (paral%io_parent)&
                           WRITE(6,'(1X,A,I5,A,I5,A,I5,A,I5)')&
                           'COORNUMGRP! ATMA=',atma,' ITYP=',ityp,&
                           ' KX=',kx,' JJ=',jj
                      CALL stopgm('COORNUMGRP','INCORRECT NUMBERS FOR ATMA',& 
                           __LINE__,__FILE__)
                   ENDIF
                   fff=1.0_real_8/REAL(jj,kind=real_8)
                   DO j=1,jj
                      IF (jj.EQ.ions1%nat) THEN
                         ibb=j
                      ELSE
                         ibb=duat%listd2(j+1,kx)
                      ENDIF
                      jx=lskptr(1,ibb)
                      jy=lskptr(2,ibb)
                      jz=lskptr(3,ibb)
                      IF (jx.NE.0) an(jx)=an(jx)-ff*dx*fff
                      IF (jy.NE.0) an(jy)=an(jy)-ff*dy*fff
                      IF (jz.NE.0) an(jz)=an(jz)-ff*dz*fff
                   ENDDO
                ELSEIF (ityp.EQ.3) THEN
                   kx=duat%listda(naa,2)
                   jj=duat%listd3(1,kx)
                   IF (jj.LT.0) THEN
                      jj=ions1%nat
                   ELSE IF (jj.EQ.0) THEN
                      IF (paral%io_parent)&
                           WRITE(6,'(1X,A,I5,A,I5,A,I5,A,I5,A,I5,A,I5,A,I5)')&
                           'COORNUMGRP! ATMA=',atma,' ATMB=',atmb,&
                           ' ITYP=',ityp,' KX=',kx,' JJ=',jj,&
                           ' KX=',kx,' JJ=',jj
                      CALL stopgm('COORNUMGRP','INCORRECT NUMBERS',& 
                           __LINE__,__FILE__)
                   ENDIF
                   fff=0.0_real_8
                   DO j=1,jj
                      IF (jj.EQ.ions1%nat) THEN
                         ibb=j
                      ELSE
                         ibb=duat%listd3(j+1,kx)
                      ENDIF
                      is=iatpt(2,ibb)
                      fff=fff+rmass%pma0(is)
                   ENDDO
                   fff=1.0_real_8/fff
                   DO j=1,jj
                      IF (jj.EQ.ions1%nat) THEN
                         ibb=j
                      ELSE
                         ibb=duat%listd3(j+1,kx)
                      ENDIF
                      is=iatpt(2,ibb)
                      jx=lskptr(1,ibb)
                      jy=lskptr(2,ibb)
                      jz=lskptr(3,ibb)
                      IF (jx.NE.0) an(jx)=an(jx)-ff*dx*fff*rmass%pma0(is)
                      IF (jy.NE.0) an(jy)=an(jy)-ff*dy*fff*rmass%pma0(is)
                      IF (jz.NE.0) an(jz)=an(jz)-ff*dz*fff*rmass%pma0(is)
                   ENDDO
                ELSEIF (ityp.EQ.4) THEN
                   kx=duat%listda(naa,2)
                   jj=duat%listd4(1,kx)
                   fff=0.0_real_8
                   DO j=1,jj
                      fff=fff+duat%weigd4(j,kx)
                   ENDDO
                   fff=1.0_real_8/fff
                   DO j=1,jj
                      ibb=duat%listd4(j+1,kx)
                      jx=lskptr(1,ibb)
                      jy=lskptr(2,ibb)
                      jz=lskptr(3,ibb)
                      IF (jx.NE.0) an(jx)=an(jx)-ff*dx*fff*duat%weigd4(j,kx)
                      IF (jy.NE.0) an(jy)=an(jy)-ff*dy*fff*duat%weigd4(j,kx)
                      IF (jz.NE.0) an(jz)=an(jz)-ff*dz*fff*duat%weigd4(j,kx)
                   ENDDO
                ENDIF
             ENDIF
130          CONTINUE
          ENDIF
       ENDDO
    ENDDO
    fci=fvi-cval
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE coornumgrp
  ! ==================================================================
  SUBROUTINE coorn_rf(ia,isp,nexp,mexp,c1_rc,r0_shift,tscr,an,&
       fci,fvi,cval,SPvsAT,iatom)
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: ia, isp, nexp, mexp
    REAL(real_8)                             :: c1_rc, r0_shift, tscr(3,*), &
                                                an(*), fci, fvi, cval
    INTEGER                                  :: SPvsAT, iatom(*)

    INTEGER                                  :: iat, ib, iend, iiat, iini, &
                                                is, ityp, j, jj, jx, jy, jz, &
                                                k, kx, l1, l2, l3, naa
    REAL(real_8)                             :: dd, dx, dy, dz, fden, ff, &
                                                fff, fnum, rden, rnum, x0, &
                                                xx(3), y0, z0, dx_, dy_, dz_

! Variables
! ==--------------------------------------------------------------==
! AK  debug:      call mm_dim_query('enter coorn_rf')

    IF (ia.LE.ions1%nat) THEN
       l1 = lskptr(1,ia)
       l2 = lskptr(2,ia)
       l3 = lskptr(3,ia)
    ENDIF
    CALL fillc(ia,tscr,xx)
    x0 = xx(1)
    y0 = xx(2)
    z0 = xx(3)

    fvi = 0._real_8
    iini = 1
    IF (SPvsAT .EQ. 1) THEN
       DO is = 1,isp-1
          iini=iini+ions0%na(is)
       ENDDO
       iend = iini + ions0%na(isp) -1
    ELSE
       iini = 1
       iend = isp
    ENDIF

    DO iiat=iini,iend
       IF (SPvsAT .EQ. 1) THEN
          iat = iiat
       ELSE
          iat = iatom(iiat)
       ENDIF
       dx_=tscr(1,iat)-x0
       dy_=tscr(2,iat)-y0
       dz_=tscr(3,iat)-z0
       IF (.NOT. isos1%tisos) THEN
          CALL pbc(dx_,dy_,dz_,dx,dy,dz,1,parm%apbc,parm%ibrav)
       ELSE
          dx=dx_
          dy=dy_
          dz=dz_
       ENDIF
       dd=SQRT(dx*dx+dy*dy+dz*dz)
       ! ..exclude self interaction
       IF (dd.GT.1.e-2_real_8) THEN
          rnum = ((dd-r0_shift)/c1_rc)**nexp
          rden = ((dd-r0_shift)/c1_rc)**(nexp+mexp)
          fnum = 1-rnum
          fden = 1-rden
          IF (ABS(fden) .LT. 1.e-10_real_8) THEN
             fvi = fvi + REAL(nexp,kind=real_8)/REAL(nexp+mexp,kind=real_8)
             ff  = -0.5_real_8*REAL(nexp*mexp,kind=real_8)/REAL(nexp+mexp,kind=real_8)/(c1_rc*dd)
          ELSE
             fvi  = fvi + fnum/fden
             ff   = -(rnum*REAL(nexp,kind=real_8)-fnum*rden*REAL(nexp+mexp,kind=real_8)/&
                  fden)/(fden*(dd-r0_shift)*dd)
          ENDIF
          IF (lskptr(1,iat).NE.0) THEN
             k=lskptr(1,iat)
             an(k)=an(k)+ff*dx
          ENDIF
          IF (lskptr(2,iat).NE.0) THEN
             k=lskptr(2,iat)
             an(k)=an(k)+ff*dy
          ENDIF
          IF (lskptr(3,iat).NE.0) THEN
             k=lskptr(3,iat)
             an(k)=an(k)+ff*dz
          ENDIF

          IF (ia.LE.ions1%nat) THEN
             IF (l1.NE.0) an(l1)=an(l1)-ff*dx
             IF (l2.NE.0) an(l2)=an(l2)-ff*dy
             IF (l3.NE.0) an(l3)=an(l3)-ff*dz
          ELSEIF (ia.LE.ions1%nat+duat%ndat) THEN
             naa=ia-ions1%nat
             ityp=duat%listda(naa,1)
             IF (ityp.EQ.1) THEN
                GOTO 110
             ELSEIF (ityp.EQ.2) THEN
                kx=duat%listda(naa,2)
                jj=duat%listd2(1,kx)
                IF (jj.LT.0) THEN
                   jj=ions1%nat
                ELSE IF (jj.EQ.0) THEN
                   IF (paral%io_parent)&
                        WRITE(6,'(1X,A,I5,A,I5,A,I5,A,I5)')&
                        'COORN_RF! IA=',ia,' ITYP=',ityp,' KX=',kx,' JJ=',jJ
                   CALL stopgm('COORN_RF','INCORRECT NUMBERS FOR IA',& 
                        __LINE__,__FILE__)
                ENDIF
                fff=1.0_real_8/REAL(jj,kind=real_8)
                DO j=1,jj
                   IF (jj.EQ.ions1%nat) THEN
                      ib=j
                   ELSE
                      ib=duat%listd2(j+1,kx)
                   ENDIF
                   jx=lskptr(1,ib)
                   jy=lskptr(2,ib)
                   jz=lskptr(3,ib)
                   IF (jx.NE.0) an(jx)=an(jx)-ff*dx*fff
                   IF (jy.NE.0) an(jy)=an(jy)-ff*dy*fff
                   IF (jz.NE.0) an(jz)=an(jz)-ff*dz*fff
                ENDDO
             ELSEIF (ityp.EQ.3) THEN
                kx=duat%listda(naa,2)
                jj=duat%listd3(1,kx)
                IF (jj.LT.0) THEN
                   jj=ions1%nat
                ELSE IF (jj.EQ.0) THEN
                   IF (paral%io_parent)&
                        WRITE(6,'(1X,A,I5,A,I5,A,I5,A,I5)')&
                        'COORN_RF! IA=',ia,' ITYP=',ityp,' KX=',kx,' JJ=',jJ
                   CALL stopgm('COORN_RF','INCORRECT NUMBERS FOR IA',& 
                        __LINE__,__FILE__)
                ENDIF
                fff=0.0_real_8
                DO j=1,jj
                   IF (jj.EQ.ions1%nat) THEN
                      ib=j
                   ELSE
                      ib=duat%listd3(j+1,kx)
                   ENDIF
                   is=iatpt(2,ib)
                   fff=fff+rmass%pma0(is)
                ENDDO
                fff=1.0_real_8/fff
                DO j=1,jj
                   IF (jj.EQ.ions1%nat) THEN
                      ib=j
                   ELSE
                      ib=duat%listd3(j+1,kx)
                   ENDIF
                   is=iatpt(2,ib)
                   jx=lskptr(1,ib)
                   jy=lskptr(2,ib)
                   jz=lskptr(3,ib)
                   IF (jx.NE.0) an(jx)=an(jx)-ff*dx*fff*rmass%pma0(is)
                   IF (jy.NE.0) an(jy)=an(jy)-ff*dy*fff*rmass%pma0(is)
                   IF (jz.NE.0) an(jz)=an(jz)-ff*dz*fff*rmass%pma0(is)
                ENDDO
             ELSEIF (ityp.EQ.4) THEN
                kx=duat%listda(naa,2)
                jj=duat%listd4(1,kx)
                fff=0.0_real_8
                DO j=1,jj
                   fff=fff+duat%weigd4(j,kx)
                ENDDO
                fff=1.0_real_8/fff
                DO j=1,jj
                   ib=duat%listd4(j+1,kx)
                   jx=lskptr(1,ib)
                   jy=lskptr(2,ib)
                   jz=lskptr(3,ib)
                   IF (jx.NE.0) an(jx)=an(jx)-ff*dx*fff*duat%weigd4(j,kx)
                   IF (jy.NE.0) an(jy)=an(jy)-ff*dy*fff*duat%weigd4(j,kx)
                   IF (jz.NE.0) an(jz)=an(jz)-ff*dz*fff*duat%weigd4(j,kx)
                ENDDO
             ENDIF
          ENDIF
110       CONTINUE
       ENDIF
    ENDDO
    fci=fvi-cval
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE coorn_rf
  ! ==================================================================
  SUBROUTINE  bndswitch(ia,ib,nexp,mexp,c1_rc,tscr,an,&
       fci,fvi,cval)
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: ia, ib, nexp, mexp
    REAL(real_8)                             :: c1_rc, tscr(3,*), an(*), fci, &
                                                fvi, cval

    INTEGER                                  :: ic, is, ityp, j, jj, jx, jy, &
                                                jz, kx, la1, la2, la3, lb1, &
                                                lb2, lb3, naa
    REAL(real_8)                             :: dd, dx, dy, dz, fden, ff, &
                                                fff, fnum, rden, rnum, xa, &
                                                xb, xx(3), ya, yb, za, zb, dx_, dy_, dz_

    IF (ia.LE.ions1%nat) THEN
       la1 = lskptr(1,ia)
       la2 = lskptr(2,ia)
       la3 = lskptr(3,ia)
    ENDIF
    CALL fillc(ia,tscr,xx)
    xa = xx(1)
    ya = xx(2)
    za = xx(3)
    IF (ia.LE.ions1%nat) THEN
       lb1 = lskptr(1,ib)
       lb2 = lskptr(2,ib)
       lb3 = lskptr(3,ib)
    ENDIF
    CALL fillc(ib,tscr,xx)
    xb = xx(1)
    yb = xx(2)
    zb = xx(3)

    fvi = 0._real_8
    dx_= xb - xa
    dy_= yb - ya
    dz_= zb - za
    IF (.NOT. isos1%tisos) THEN
       CALL pbc(dx_,dy_,dz_,dx,dy,dz,1,parm%apbc,parm%ibrav)
    ELSE
       dx=dx_
       dy=dy_
       dz=dz_
    ENDIF
    dd=SQRT(dx*dx+dy*dy+dz*dz)
    IF (dd.GT.1.e-2_real_8) THEN
       rnum = (dd/c1_rc)**nexp
       rden = (dd/c1_rc)**(nexp+mexp)
       fnum = 1-rnum
       fden = 1-rden
       IF (ABS(fden) .LT. 1.e-10_real_8) THEN
          fvi = fvi + REAL(nexp,kind=real_8)/REAL(nexp+mexp,kind=real_8)
          ff  = -0.5_real_8*REAL(nexp*mexp,kind=real_8)/REAL(nexp+mexp,kind=real_8)/(c1_rc*dd)
       ELSE
          fvi  = fvi + fnum/fden
          ff   = -(rnum*REAL(nexp,kind=real_8)-fnum*rden*REAL(nexp+mexp,kind=real_8)/fden)/&
               (fden*dd*dd)
       ENDIF
       IF (ib .LE. ions1%nat) THEN
          IF (lb1.NE.0) an(lb1)=an(la1)+ff*dx
          IF (lb2.NE.0) an(lb2)=an(la2)+ff*dy
          IF (lb3.NE.0) an(lb3)=an(la3)+ff*dz
       ELSEIF (ib.LE.ions1%nat+duat%ndat) THEN
          naa=ib-ions1%nat
          ityp=duat%listda(naa,1)
          IF (ityp.EQ.1) THEN
             GOTO 110
          ELSEIF (ityp.EQ.2) THEN
             kx=duat%listda(naa,2)
             jj=duat%listd2(1,kx)
             IF (jj.LT.0) THEN
                jj=ions1%nat
             ELSE IF (jj.EQ.0) THEN
                IF (paral%io_parent)&
                     WRITE(6,'(1X,A,I5,A,I5,A,I5,A,I5)')&
                     'BNDSWITCH! IB=',ib,' ITYP=',ityp,' KX=',kx,' JJ=',jJ
                CALL stopgm('BNDSWITCH','INCORRECT NUMBERS FOR IB',& 
                     __LINE__,__FILE__)
             ENDIF
             fff=1.0_real_8/REAL(jj,kind=real_8)
             DO j=1,jj
                IF (jj.EQ.ions1%nat) THEN
                   ic=j
                ELSE
                   ic=duat%listd2(j+1,kx)
                ENDIF
                jx=lskptr(1,ic)
                jy=lskptr(2,ic)
                jz=lskptr(3,ic)
                IF (jx.NE.0) an(jx)=an(jx)+ff*dx*fff
                IF (jy.NE.0) an(jy)=an(jy)+ff*dy*fff
                IF (jz.NE.0) an(jz)=an(jz)+ff*dz*fff
             ENDDO
          ELSEIF (ityp.EQ.3) THEN
             kx=duat%listda(naa,2)
             jj=duat%listd3(1,kx)
             IF (jj.LT.0) THEN
                jj=ions1%nat
             ELSE IF (jj.EQ.0) THEN
                IF (paral%io_parent)&
                     WRITE(6,'(1X,A,I5,A,I5,A,I5,A,I5)')&
                     'BNDSWITCH! IB=',ib,' ITYP=',ityp,' KX=',kx,' JJ=',jJ
                CALL stopgm('BNDSWITCH','INCORRECT NUMBERS FOR IB',& 
                     __LINE__,__FILE__)
             ENDIF
             fff=0.0_real_8
             DO j=1,jj
                IF (jj.EQ.ions1%nat) THEN
                   ic=j
                ELSE
                   ic=duat%listd3(j+1,kx)
                ENDIF
                is=iatpt(2,ic)
                fff=fff+rmass%pma0(is)
             ENDDO
             fff=1.0_real_8/fff
             DO j=1,jj
                IF (jj.EQ.ions1%nat) THEN
                   ic=j
                ELSE
                   ic=duat%listd3(j+1,kx)
                ENDIF
                is=iatpt(2,ic)
                jx=lskptr(1,ic)
                jy=lskptr(2,ic)
                jz=lskptr(3,ic)
                IF (jx.NE.0) an(jx)=an(jx)+ff*dx*fff*rmass%pma0(is)
                IF (jy.NE.0) an(jy)=an(jy)+ff*dy*fff*rmass%pma0(is)
                IF (jz.NE.0) an(jz)=an(jz)+ff*dz*fff*rmass%pma0(is)
             ENDDO
          ELSEIF (ityp.EQ.4) THEN
             kx=duat%listda(naa,2)
             jj=duat%listd4(1,kx)
             fff=0.0_real_8
             DO j=1,jj
                fff=fff+duat%weigd4(j,kx)
             ENDDO
             fff=1.0_real_8/fff
             DO j=1,jj
                ic=duat%listd4(j+1,kx)
                jx=lskptr(1,ic)
                jy=lskptr(2,ic)
                jz=lskptr(3,ic)
                IF (jx.NE.0) an(jx)=an(jx)+ff*dx*fff*duat%weigd4(j,kx)
                IF (jy.NE.0) an(jy)=an(jy)+ff*dy*fff*duat%weigd4(j,kx)
                IF (jz.NE.0) an(jz)=an(jz)+ff*dz*fff*duat%weigd4(j,kx)
             ENDDO
          ENDIF
       ENDIF
110    CONTINUE

       IF (ia .LE. ions1%nat) THEN
          IF (la1.NE.0) an(la1)=an(la1)-ff*dx
          IF (la2.NE.0) an(la2)=an(la2)-ff*dy
          IF (la3.NE.0) an(la3)=an(la3)-ff*dz
       ELSEIF (ia.LE.ions1%nat+duat%ndat) THEN
          naa=ia-ions1%nat
          ityp=duat%listda(naa,1)
          IF (ityp.EQ.1) THEN
             GOTO 120
          ELSEIF (ityp.EQ.2) THEN
             kx=duat%listda(naa,2)
             jj=duat%listd2(1,kx)
             IF (jj.LT.0) THEN
                jj=ions1%nat
             ELSE IF (jj.EQ.0) THEN
                IF (paral%io_parent)&
                     WRITE(6,'(1X,A,I5,A,I5,A,I5,A,I5)')&
                     'BNDSWITCH! IA=',ia,' ITYP=',ityp,' KX=',kx,' JJ=',jJ
                CALL stopgm('BNDSWITCH','INCORRECT NUMBERS FOR IA',& 
                     __LINE__,__FILE__)
             ENDIF
             fff=1.0_real_8/REAL(jj,kind=real_8)
             DO j=1,jj
                IF (jj.EQ.ions1%nat) THEN
                   ic=j
                ELSE
                   ic=duat%listd2(j+1,kx)
                ENDIF
                jx=lskptr(1,ic)
                jy=lskptr(2,ic)
                jz=lskptr(3,ic)
                IF (jx.NE.0) an(jx)=an(jx)-ff*dx*fff
                IF (jy.NE.0) an(jy)=an(jy)-ff*dy*fff
                IF (jz.NE.0) an(jz)=an(jz)-ff*dz*fff
             ENDDO
          ELSEIF (ityp.EQ.3) THEN
             kx=duat%listda(naa,2)
             jj=duat%listd3(1,kx)
             IF (jj.LT.0) THEN
                jj=ions1%nat
             ELSE IF (jj.EQ.0) THEN
                IF (paral%io_parent)&
                     WRITE(6,'(1X,A,I5,A,I5,A,I5,A,I5)')&
                     'BNDSWITCH! IA=',ia,' ITYP=',ityp,' KX=',kx,' JJ=',jJ
                CALL stopgm('BNDSWITCH','INCORRECT NUMBERS FOR IA',& 
                     __LINE__,__FILE__)
             ENDIF
             fff=0.0_real_8
             DO j=1,jj
                IF (jj.EQ.ions1%nat) THEN
                   ic=j
                ELSE
                   ic=duat%listd3(j+1,kx)
                ENDIF
                is=iatpt(2,ic)
                fff=fff+rmass%pma0(is)
             ENDDO
             fff=1.0_real_8/fff
             DO j=1,jj
                IF (jj.EQ.ions1%nat) THEN
                   ic=j
                ELSE
                   ic=duat%listd3(j+1,kx)
                ENDIF
                is=iatpt(2,ic)
                jx=lskptr(1,ic)
                jy=lskptr(2,ic)
                jz=lskptr(3,ic)
                IF (jx.NE.0) an(jx)=an(jx)-ff*dx*fff*rmass%pma0(is)
                IF (jy.NE.0) an(jy)=an(jy)-ff*dy*fff*rmass%pma0(is)
                IF (jz.NE.0) an(jz)=an(jz)-ff*dz*fff*rmass%pma0(is)
             ENDDO
          ELSEIF (ityp.EQ.4) THEN
             kx=duat%listda(naa,2)
             jj=duat%listd4(1,kx)
             fff=0.0_real_8
             DO j=1,jj
                fff=fff+duat%weigd4(j,kx)
             ENDDO
             fff=1.0_real_8/fff
             DO j=1,jj
                ic=duat%listd4(j+1,kx)
                jx=lskptr(1,ic)
                jy=lskptr(2,ic)
                jz=lskptr(3,ic)
                IF (jx.NE.0) an(jx)=an(jx)-ff*dx*fff*duat%weigd4(j,kx)
                IF (jy.NE.0) an(jy)=an(jy)-ff*dy*fff*duat%weigd4(j,kx)
                IF (jz.NE.0) an(jz)=an(jz)-ff*dz*fff*duat%weigd4(j,kx)
             ENDDO
          ENDIF
       ENDIF
120    CONTINUE
    ENDIF
    fci=fvi-cval
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE bndswitch
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE restfc(tau0,fion)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:), fion(:,:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'restfc'

    INTEGER                                  :: dum(1), i, ia, ib, ic, id, &
                                                ierr, isub, ityp, j, k, kmax, &
                                                n
    REAL(real_8)                             :: cval, dx(12), r0_shift, scal, &
                                                sign0, x1(3), x2(3), x3(3), &
                                                x4(3)
    REAL(real_8), ALLOCATABLE                :: tscr(:,:,:)

    ener_com%erestr=0.0_real_8
    IF (cotr007%mrestr.EQ.0) RETURN
    CALL tiset('    RESTFC',isub)
    ALLOCATE(tscr(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL zeroing(tscr)
    CALL dum2(tau0,tscr)
    CALL zeroing(anorm)!,cotr007%mrestr*cotc0%nodim)
    DO i=1,cotr007%mrestr
       cval=resval(i)
       ityp=ntrest(1,i)
       ia=ntrest(2,i)
       ib=ntrest(3,i)
       ic=ntrest(4,i)
       id=ntrest(5,i)
       IF (ityp.EQ.1) THEN
          ! ..stretch
          CALL fillc(ia,tscr,x1)
          CALL fillc(ib,tscr,x2)
          CALL funcr(resc(i),resv(i),cval,x1,x2)
          CALL diffr(dx,x1,x2)
          kmax=6
          ! ..angle
       ELSEIF (ityp.EQ.2) THEN
          CALL fillc(ia,tscr,x1)
          CALL fillc(ib,tscr,x2)
          CALL fillc(ic,tscr,x3)
          CALL funct(resc(i),resv(i),cval,x1,x2,x3)
          CALL difft(dx,x1,x2,x3)
          kmax=9
       ELSEIF (ityp.EQ.3) THEN
          ! ..dihedral angle
          CALL fillc(ia,tscr,x1)
          CALL fillc(ib,tscr,x2)
          CALL fillc(ic,tscr,x3)
          CALL fillc(id,tscr,x4)
          CALL funco(resc(i),resv(i),cval,x1,x2,x3,x4,sign0)
          CALL diffo(dx,x1,x2,x3,x4,sign0)
          kmax=12
       ELSEIF (ityp.EQ.4) THEN
          ! ..distance
          CALL fillc(ia,tscr,x1)
          CALL fillc(ib,tscr,x2)
          CALL funcd(resc(i),resv(i),cval,x1,x2)
          CALL diffd(dx,x1,x2)
          kmax=6
       ELSEIF (ityp.EQ.5) THEN
          ! ..out of plane
          CALL fillc(ia,tscr,x1)
          CALL fillc(ib,tscr,x2)
          CALL fillc(ic,tscr,x3)
          CALL fillc(id,tscr,x4)
          CALL funcp(resc(i),resv(i),cval,x1,x2,x3,x4,dx)
          kmax=12
       ELSEIF (ityp.EQ.6) THEN
          ! ..coordination number
          c_kappa = respar(1,i)
          c_rc    = respar(2,i)
          CALL coornum(ia,c_kappa,c_rc,tscr,ions1%nat,anorm(1,i),lskptr,&
               resc(i),resv(i),cval)
          kmax=0
       ELSEIF (ityp.EQ.7) THEN
          ! ..difference between distances- 3 atoms
          CALL fillc(ia,tscr,x1)
          CALL fillc(ib,tscr,x2)
          CALL fillc(ic,tscr,x3)
          CALL funcdd(resc(i),resv(i),cval,x1,x2,x3)
          CALL diffdd(dx,x1,x2,x3)
          kmax=9
       ELSEIF (ityp.EQ.8) THEN
          ! ..species dependent coordination number with Fermi function
          c_kappa = respar(1,i)
          c_rc    = respar(2,i)
          CALL coornumsp(ia,ib,c_kappa,c_rc,tscr,anorm(1,i),&
               lskptr,resc(i),resv(i),cval)
          kmax=0
       ELSEIF (ityp.EQ.9) THEN
          ! ..species dependent coordination number with rational function
          c_rc    = respar(1,i)
          r0_shift  = respar(2,i)
          CALL coorn_rf(ia,ib,ic,id,c_rc,r0_shift,tscr,anorm(1,i),&
               resc(i),resv(i),cval,1,dum)
          kmax=0
       ELSEIF (ityp.EQ.10) THEN
          ! ..bond switch with rational function
          c_rc    = respar(1,i)
          CALL bndswitch(ia,ib,ic,id,c_rc,tscr,anorm(1,i),&
               resc(i),resv(i),cval)
          kmax=0
       ELSEIF (ityp.EQ.11) THEN
          ! ..total species dependent coordination number with rational function
          c_rc    = respar(1,i)
          r0_shift  = respar(2,i)
          CALL coorntot_rf(ia,ib,ic,id,c_rc,r0_shift,&
               tscr,anorm(1,i),resc(i),resv(i),cval,1,dum)
          kmax=0
       ELSEIF (ityp.EQ.12) THEN
          ! ..distance between two atoms along a selected direction
          CALL fillc(ia,tscr,x1)
          CALL fillc(ib,tscr,x2)
          DO k=1,3
             IF (k.NE.ic) THEN
                x1(k)=0.0_real_8
                x2(k)=0.0_real_8
             ENDIF
          ENDDO
          CALL funcd(resc(i),resv(i),cval,x1,x2)
          CALL diffd(dx,x1,x2)
          kmax=6
       ELSEIF (ityp.EQ.13) THEN! cmb-kk
          ! ..restraint positional
          CALL fillc(ia,tscr,x1)
          DO k=1,3
             x2(k)=respos(k,i)
          ENDDO
          CALL funcd(resc(i),resv(i),cval,x1,x2)
          CALL diffd(dx,x1,x2)
          kmax=3
       ELSE
          WRITE(6,*) ' RESTFC! I=',I,' TYPE=',ITYP
          CALL stopgm('RESTFC','UNKNOWN TYPE OF RESTRAINT',& 
               __LINE__,__FILE__)
       ENDIF

       DO k=1,kmax
          DO n=1,cotc0%nodim
             anorm(n,i)=anorm(n,i)+dx(k)*rskel(n,i,k)
          ENDDO
       ENDDO
    ENDDO

    CALL zeroing(fcstr)!,cotc0%nodim)
    CALL zeroing(resfdiff)!,cotr007%mrestr)
    CALL zeroing(resm)!,cotr007%mrestr*cotr007%mrestr)

    IF (.NOT.cotr007%lhyperplane) THEN
       DO i=1,cotr007%mrestr
          DO j=1,cotc0%nodim
             fcstr(j)=fcstr(j)+resfor(i)*anorm(j,i)*resc(i)
          ENDDO
          resfdiff(i)=resfor(i)*resc(i)
          CALL zeroing(resm(:,i))!,cotr007%mrestr)
          DO j=1,cotr007%mrestr
             DO k=1,cotc0%nodim
                resm(j,i)=resm(j,i)+anorm(k,i)*anorm(k,j)/rsmass(k)
             ENDDO
          ENDDO
          ener_com%erestr=ener_com%erestr+0.5_real_8*resfor(i)*resc(i)*resc(i)
       ENDDO
    ELSE
       scal=0.0_real_8
       DO i=1,cotr007%mrestr
          resfdiff(i)=resfor(i)*resc(i)
          scal=scal+resfdiff(i)
       ENDDO
       DO i=1,cotr007%mrestr
          DO j=1,cotc0%nodim
             fcstr(j)=fcstr(j)+resfor(i)*anorm(j,i)*scal
          ENDDO
          CALL zeroing(resm(:,i))!,cotr007%mrestr)
          DO j=1,cotr007%mrestr
             DO k=1,cotc0%nodim
                resm(j,i)=resm(j,i)+anorm(k,i)*anorm(k,j)
             ENDDO
          ENDDO
          ener_com%erestr=ener_com%erestr+0.5_real_8*resfor(i)*resc(i)*scal
       ENDDO
    ENDIF

    CALL gettau(tscr,fcstr)
    CALL daxpy(3*maxsys%nax*maxsys%nsx,-1._real_8,tscr,1,fion,1)
    DEALLOCATE(tscr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    CALL tihalt('    RESTFC',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE restfc
  ! ==================================================================

END MODULE cnstfc_utils
