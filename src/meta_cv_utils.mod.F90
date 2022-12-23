MODULE meta_cv_utils
  USE adat,                            ONLY: elem
  USE cnst_dyn,                        ONLY: file_str_ab,&
                                             icv_rmsd_ab,&
                                             ncolvar,&
                                             nrmsd_ab
  USE constr_utils,                    ONLY: getscal,&
                                             vecprod
  USE cotr,                            ONLY: cotc0,&
                                             duat,&
                                             lskptr
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_new,&
                                             fo_old
  USE fillc_utils,                     ONLY: fillc
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE isos,                            ONLY: isos1
  USE kinds,                           ONLY: real_8
  USE latgen_utils,                    ONLY: omegagen
  USE meta_colvar_util_utils,          ONLY: fillc2
  USE metr,                            ONLY: metr_com
  USE parac,                           ONLY: paral
  USE pbc_utils,                       ONLY: pbc
  USE system,                          ONLY: parm
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: coorntot_rf
  PUBLIC :: abs_difcntot_rf
  PUBLIC :: difcntot_rf
  PUBLIC :: ctot_chain
  PUBLIC :: disp_lnm
  PUBLIC :: rmsd_rot
  PUBLIC :: funcp_mia
  PUBLIC :: aplane
  PUBLIC :: rmsd_ab
  PUBLIC :: val_rmsd
  PUBLIC :: hyd_presence
  PUBLIC :: hyd_ion_distance
  PUBLIC :: funcl
  PUBLIC :: volume_cv
  PUBLIC :: side_cv
  PUBLIC :: angle_cv
  PUBLIC :: hyd_cutoff
  PUBLIC :: cut_hyd_ion_distance
  PUBLIC :: cdipolerho
  PUBLIC :: cdipolephi
  PUBLIC :: cdipoletheta
  PUBLIC :: vcoors

CONTAINS

  ! ==================================================================
  SUBROUTINE coorntot_rf(isp1,isp2,nexp,mexp,c1_rc,r0_shift,tscr,an,&
       fci,fvi,cval,SPvsAT,iatom)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: isp1, isp2, nexp, mexp
    REAL(real_8)                             :: c1_rc, r0_shift, tscr(3,*), &
                                                an(*), fci, fvi, cval
    INTEGER                                  :: SPvsAT, iatom(*)

    INTEGER                                  :: iat1, iat2, ib, iend1, iend2, &
                                                iiat2, iini1, iini2, is, &
                                                ityp, j, jj, jx, jy, jz, k, &
                                                kx, l1, l2, l3, naa
    REAL(real_8)                             :: dd, dx, dy, dz, fden, ff, &
                                                fff, fnum, rden, rnum, una1, &
                                                x0, xx(3), y0, z0

    iini1 = 1
    DO is = 1,isp1-1
       iini1=iini1+ions0%na(is)
    ENDDO
    iend1 = iini1 + ions0%na(isp1) -1
    una1 = 1.0_real_8/REAL(ions0%na(isp1),kind=real_8)
    ! write(6,*) 'IINI = ',IINI,'IEND = ',IEND
    IF (SPvsAT .EQ. 1) THEN
       iini2 = 1
       DO is = 1,isp2-1
          iini2=iini2+ions0%na(is)
       ENDDO
       iend2 = iini2 + ions0%na(isp2) -1
    ELSE
       iini2 = 1
       iend2 = isp2
    ENDIF
    ! write(6,*) 'IINI = ',IINI,'IEND = ',IEND

    fvi = 0._real_8
    DO iat1 = iini1,iend1
       l1 = lskptr(1,iat1)
       l2 = lskptr(2,iat1)
       l3 = lskptr(3,iat1)
       x0 = tscr(1,iat1)
       y0 = tscr(2,iat1)
       z0 = tscr(3,iat1)

       DO iiat2=iini2,iend2
          IF (SPvsAT .EQ. 1) THEN
             iat2 = iiat2
          ELSE
             iat2 = iatom(iiat2)
          ENDIF
          CALL fillc(iat2,tscr,xx)
          dx=xx(1)-x0
          dy=xx(2)-y0
          dz=xx(3)-z0
          IF (.NOT. isos1%tisos) THEN
             CALL pbc(dx,dy,dz,dx,dy,dz,1,parm%apbc,parm%ibrav)
          ENDIF
          dd=SQRT(dx*dx+dy*dy+dz*dz)
          ! ..exclude self interaction
          ! write(6,*) IAT1,IAT2,DD,DX,DY,DZ
          IF (dd.GT.1.e-2_real_8) THEN
             rnum = ((dd-r0_shift)/c1_rc)**nexp
             rden = ((dd-r0_shift)/c1_rc)**(nexp+mexp)
             fnum = 1-rnum
             fden = 1-rden
             IF (ABS(fden) .LT. 1.e-10_real_8) THEN
                fvi = fvi + REAL(nexp,kind=real_8)/REAL(nexp+mexp,kind=real_8)
                ff  = -0.5_real_8*REAL(nexp*mexp,kind=real_8)/REAL(nexp+mexp,kind=real_8)/&
                     (c1_rc*dd)
             ELSE
                fvi  = fvi + fnum/fden
                ff   = -(rnum*REAL(nexp,kind=real_8)-fnum*rden*REAL(nexp+mexp,kind=real_8)/&
                     fden)/(fden*(dd-r0_shift)*dd)*una1
             ENDIF

             IF (iat2.LE.ions1%nat) THEN
                IF (lskptr(1,iat2).NE.0) THEN
                   k=lskptr(1,iat2)
                   an(k)=an(k)+ff*dx
                ENDIF
                IF (lskptr(2,iat2).NE.0) THEN
                   k=lskptr(2,iat2)
                   an(k)=an(k)+ff*dy
                ENDIF
                IF (lskptr(3,iat2).NE.0) THEN
                   k=lskptr(3,iat2)
                   an(k)=an(k)+ff*dz
                ENDIF
             ELSEIF (iat2.LE.ions1%nat+duat%ndat) THEN
                naa=iat2-ions1%nat
                ityp=duat%listda(naa,1)
                IF (ityp.EQ.1) THEN
                   GOTO 110
                ELSEIF (ityp.EQ.2) THEN
                   kx=duat%listda(naa,2)
                   jj=duat%listd2(1,kx)
                   fff=1.0_real_8/REAL(jj,kind=real_8)
                   DO j=1,jj
                      ib=duat%listd2(j+1,kx)
                      jx=lskptr(1,ib)
                      jy=lskptr(2,ib)
                      jz=lskptr(3,ib)
                      IF (jx.NE.0) an(jx)=an(jx)+ff*dx*fff
                      IF (jy.NE.0) an(jy)=an(jy)+ff*dy*fff
                      IF (jz.NE.0) an(jz)=an(jz)+ff*dz*fff
                   ENDDO
                ENDIF
             ENDIF
110          CONTINUE

             IF (l1.NE.0) an(l1)=an(l1)-ff*dx
             IF (l2.NE.0) an(l2)=an(l2)-ff*dy
             IF (l3.NE.0) an(l3)=an(l3)-ff*dz
          ENDIF
       ENDDO
    ENDDO
    fvi = fvi*una1
    fci=fvi-cval
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE coorntot_rf
  ! ==================================================================
  SUBROUTINE abs_difcntot_rf(isp1,isp2,isp3,nexp,mexp,c1_rc,&
       r0_shift,tscr,an,fci,fvi,cval,SPvsAT,&
       iatom)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: isp1, isp2, isp3, nexp, mexp
    REAL(real_8)                             :: c1_rc, r0_shift, tscr(3,*), &
                                                an(*), fci, fvi, cval
    INTEGER                                  :: SPvsAT, iatom(*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'abs_difcntot_rf'

    INTEGER                                  :: ierr, j
    REAL(real_8)                             :: fci_loc, fvi_loc
    REAL(real_8), ALLOCATABLE                :: an_loc(:)

    ALLOCATE(an_loc(cotc0%nodim*ncolvar),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(an_loc)!,cotc0%nodim*ncolvar)
    ! First evaluate the difference between coordination numbers
    CALL  difcntot_rf(isp1,isp2,isp3,nexp,mexp,c1_rc,r0_shift,&
         tscr,an_loc,fci_loc,fvi_loc,cval,SPvsAT,iatom)

    fci = ABS(fci_loc)
    fvi = ABS(fvi_loc+cval)-cval
    DO j = 1, cotc0%nodim*ncolvar
       an(j) = an(j) + an_loc(j) * SIGN(1.0_real_8,fci_loc)
    ENDDO
    DEALLOCATE(an_loc,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE abs_difcntot_rf
  ! ==================================================================
  SUBROUTINE difcntot_rf(isp1,isp2,isp3,nexp,mexp,c1_rc,r0_shift,&
       tscr,an,fci,fvi,cval,SPvsAT,iatom)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: isp1, isp2, isp3, nexp, mexp
    REAL(real_8)                             :: c1_rc, r0_shift, tscr(3,*), &
                                                an(*), fci, fvi, cval
    INTEGER                                  :: SPvsAT, iatom(*)

    INTEGER                                  :: iat1, iat2, iat3, iend1, &
                                                iend2, iend3, iini1, iini2, &
                                                iini3, is, k, l1, l2, l3
    REAL(real_8)                             :: dd, dx, dy, dz, fden, ff, &
                                                fnum, rden, rnum, una1, x0, &
                                                xx(3), y0, z0

    una1 = 1.0_real_8/REAL(ions0%na(isp2),kind=real_8)

    iini1 = 1
    DO is = 1,isp1-1
       iini1=iini1+ions0%na(is)
    ENDDO
    iend1 = iini1 + ions0%na(isp1) -1

    iini2 = 1
    DO is = 1,isp2-1
       iini2=iini2+ions0%na(is)
    ENDDO
    iend2 = iini2 + ions0%na(isp2) -1

    iini3 = 1
    DO is = 1,isp3-1
       iini3=iini3+ions0%na(is)
    ENDDO
    iend3 = iini3 + ions0%na(isp3) -1

    fvi = 0._real_8
    DO iat1 = iini1,iend1
       l1 = lskptr(1,iat1)
       l2 = lskptr(2,iat1)
       l3 = lskptr(3,iat1)
       x0 = tscr(1,iat1)
       y0 = tscr(2,iat1)
       z0 = tscr(3,iat1)

       DO iat2=iini2,iend2
          CALL fillc(iat2,tscr,xx)
          dx=xx(1)-x0
          dy=xx(2)-y0
          dz=xx(3)-z0
          IF (.NOT. isos1%tisos) THEN
             CALL pbc(dx,dy,dz,dx,dy,dz,1,parm%apbc,parm%ibrav)
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
                ff  = -0.5_real_8*REAL(nexp*mexp,kind=real_8)/REAL(nexp+mexp,kind=real_8)/&
                     (c1_rc*dd)
             ELSE
                fvi  = fvi + fnum/fden
                ff   = -(rnum*REAL(nexp,kind=real_8)-fnum*rden*REAL(nexp+mexp,kind=real_8)/&
                     fden)/(fden*(dd-r0_shift)*dd)*una1
             ENDIF

             IF (lskptr(1,iat2).NE.0) THEN
                k=lskptr(1,iat2)
                an(k)=an(k)+ff*dx
             ENDIF
             IF (lskptr(2,iat2).NE.0) THEN
                k=lskptr(2,iat2)
                an(k)=an(k)+ff*dy
             ENDIF
             IF (lskptr(3,iat2).NE.0) THEN
                k=lskptr(3,iat2)
                an(k)=an(k)+ff*dz
             ENDIF

             IF (l1.NE.0) an(l1)=an(l1)-ff*dx
             IF (l2.NE.0) an(l2)=an(l2)-ff*dy
             IF (l3.NE.0) an(l3)=an(l3)-ff*dz
          ENDIF
       ENDDO
    ENDDO


    DO iat3 = iini3,iend3
       l1 = lskptr(1,iat3)
       l2 = lskptr(2,iat3)
       l3 = lskptr(3,iat3)
       x0 = tscr(1,iat3)
       y0 = tscr(2,iat3)
       z0 = tscr(3,iat3)

       DO iat2=iini2,iend2
          CALL fillc(iat2,tscr,xx)
          dx=xx(1)-x0
          dy=xx(2)-y0
          dz=xx(3)-z0
          IF (.NOT. isos1%tisos) THEN
             CALL pbc(dx,dy,dz,dx,dy,dz,1,parm%apbc,parm%ibrav)
          ENDIF
          dd=SQRT(dx*dx+dy*dy+dz*dz)
          ! ..exclude self interaction
          IF (dd.GT.1.e-2_real_8) THEN
             rnum = ((dd-r0_shift)/c1_rc)**nexp
             rden = ((dd-r0_shift)/c1_rc)**(nexp+mexp)
             fnum = 1-rnum
             fden = 1-rden
             IF (ABS(fden) .LT. 1.e-10_real_8) THEN
                ! ..change sign, compare with previous cycle
                fvi = fvi - REAL(nexp,kind=real_8)/REAL(nexp+mexp,kind=real_8)
                ff  = 0.5_real_8*REAL(nexp*mexp,kind=real_8)/REAL(nexp+mexp,kind=real_8)/&
                     (c1_rc*dd)
             ELSE
                fvi  = fvi - fnum/fden
                ff   = +(rnum*REAL(nexp,kind=real_8)-fnum*rden*REAL(nexp+mexp,kind=real_8)/&
                     fden)/(fden*(dd-r0_shift)*dd)*una1
             ENDIF

             IF (lskptr(1,iat2).NE.0) THEN
                k=lskptr(1,iat2)
                an(k)=an(k)+ff*dx
             ENDIF
             IF (lskptr(2,iat2).NE.0) THEN
                k=lskptr(2,iat2)
                an(k)=an(k)+ff*dy
             ENDIF
             IF (lskptr(3,iat2).NE.0) THEN
                k=lskptr(3,iat2)
                an(k)=an(k)+ff*dz
             ENDIF

             IF (l1.NE.0) an(l1)=an(l1)-ff*dx
             IF (l2.NE.0) an(l2)=an(l2)-ff*dy
             IF (l3.NE.0) an(l3)=an(l3)-ff*dz
          ENDIF
       ENDDO
    ENDDO

    fvi = fvi*una1
    fci = fvi-cval
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE difcntot_rf
  ! ==================================================================
  SUBROUTINE ctot_chain(isp1,isp2,isp3,nexp,mexp,c_rc1,c_rc2,tscr,&
       an,lsk,fci,fvi,cval)
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: isp1, isp2, isp3, nexp, mexp
    REAL(real_8)                             :: c_rc1, c_rc2, tscr(3,*), an(*)
    INTEGER                                  :: lsk(3,ions1%nat)
    REAL(real_8)                             :: fci, fvi, cval

    INTEGER                                  :: iat1, iat2, iat3, ienda, &
                                                iendb, iendc, iinia, iinib, &
                                                iinic, is, l1a, l1b, l1c, &
                                                l2a, l2b, l2c, l3a, l3b, l3c
    REAL(real_8) :: cc, ddab, ddac, ddbc, dfab, dfbc, dxab, dxac, dxbc, dyab, &
      dyac, dybc, dzab, dzac, dzbc, fdenab, fdenbc, fnumab, fnumbc, rdenab, &
      rdenbc, rnumab, rnumbc, una1, unb1, xa, xb, xc, ya, yb, yc, za, zb, zc

    iinia = 1
    DO is = 1,isp1-1
       iinia=iinia+ions0%na(is)
    ENDDO
    ienda = iinia + ions0%na(isp1) -1
    una1 = 1.0_real_8/REAL(ions0%na(isp1),kind=real_8)

    iinib = 1
    DO is = 1,isp2-1
       iinib=iinib+ions0%na(is)
    ENDDO
    iendb = iinib + ions0%na(isp2) -1
    unb1 = 1.0_real_8/REAL(ions0%na(isp2),kind=real_8)

    iinic = 1
    DO is = 1,isp3-1
       iinic=iinic+ions0%na(is)
    ENDDO
    iendc = iinic + ions0%na(isp3) -1

    fvi = 0._real_8
    DO iat1 = iinia,ienda
       l1a = lsk(1,iat1)
       l2a = lsk(2,iat1)
       l3a = lsk(3,iat1)
       xa = tscr(1,iat1)
       ya = tscr(2,iat1)
       za = tscr(3,iat1)

       DO iat2=iinib,iendb
          l1b = lsk(1,iat2)
          l2b = lsk(2,iat2)
          l3b = lsk(3,iat2)
          xb = tscr(1,iat2)
          yb = tscr(2,iat2)
          zb = tscr(3,iat2)
          dxab=xb-xa
          dyab=yb-ya
          dzab=zb-za
          IF (.NOT. isos1%tisos) THEN
             CALL pbc(dxab,dyab,dzab,dxab,dyab,dzab,1,parm%apbc,parm%ibrav)
          ENDIF
          ddab=SQRT(dxab*dxab+dyab*dyab+dzab*dzab)
          cc = 0.0_real_8

          DO iat3=iinic,iendc
             l1c = lsk(1,iat3)
             l2c = lsk(2,iat3)
             l3c = lsk(3,iat3)
             xc = tscr(1,iat3)
             yc = tscr(2,iat3)
             zc = tscr(3,iat3)
             dxac=xc-xa
             dyac=yc-ya
             dzac=zc-za
             IF (.NOT. isos1%tisos) THEN
                CALL pbc(dxac,dyac,dzac,dxac,dyac,dzac,1,parm%apbc,parm%ibrav)
             ENDIF
             ddac=SQRT(dxac*dxac+dyac*dyac+dzac*dzac)

             dxbc=xc-xb
             dybc=yc-yb
             dzbc=zc-zb
             IF (.NOT. isos1%tisos) THEN
                CALL pbc(dxbc,dybc,dzbc,dxbc,dybc,dzbc,1,parm%apbc,parm%ibrav)
             ENDIF
             ddbc=SQRT(dxbc*dxbc+dybc*dybc+dzbc*dzbc)

             ! ..exclude self interaction

             IF (ddab*ddbc.GT.1.e-3_real_8) THEN
                rnumab = (ddab/c_rc1)**nexp
                rdenab = (ddab/c_rc1)**(nexp+mexp)
                rnumbc = (ddbc/c_rc2)**nexp
                rdenbc = (ddbc/c_rc2)**(nexp+mexp)

                fnumab = 1-rnumab
                fdenab = 1-rdenab
                fnumbc = 1-rnumbc
                fdenbc = 1-rdenbc

                IF (ABS(fdenbc) .LT. 1.e-10_real_8) fdenbc = 1.e-10_real_8
                IF (ABS(fdenab) .LT. 1.e-10_real_8) fdenab = 1.e-10_real_8

                fvi  = fvi + (fnumab/fdenab*fnumbc/fdenbc)
                fci = 1.0_real_8

                dfab = -(rnumab*REAL(nexp,kind=real_8)-&
                     fnumab*rdenab*REAL(nexp+mexp,kind=real_8)/fdenab)/&
                     (fdenab*ddab*ddab)
                dfbc = -(rnumbc*REAL(nexp,kind=real_8)-&
                     fnumbc*rdenbc*REAL(nexp+mexp,kind=real_8)/fdenbc)/&
                     (fdenbc*ddbc*ddbc)

                dfab = dfab*fnumbc/fdenbc*una1*fci
                dfbc = dfbc*fnumab/fdenab*una1*fci

                IF (l1a.NE.0) an(l1a)=an(l1a)-dfab*dxab
                IF (l2a.NE.0) an(l2a)=an(l2a)-dfab*dyab
                IF (l3a.NE.0) an(l3a)=an(l3a)-dfab*dzab

                IF (l1b.NE.0) an(l1b)=an(l1b)+dfab*dxab-dfbc*dxbc
                IF (l2b.NE.0) an(l2b)=an(l2b)+dfab*dyab-dfbc*dybc
                IF (l3b.NE.0) an(l3b)=an(l3b)+dfab*dzab-dfbc*dzbc

                IF (l1c.NE.0) an(l1c)=an(l1c)+dfbc*dxbc
                IF (l2c.NE.0) an(l2c)=an(l2c)+dfbc*dybc
                IF (l3c.NE.0) an(l3c)=an(l3c)+dfbc*dzbc
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    fvi = fvi*una1
    fci=fvi-cval
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE ctot_chain
  ! ==================================================================
  SUBROUTINE disp_lnm(numspec,numspec2,ispec,miller,iatom,tscr,&
       natom,an,lsk,scalp,SPvsAT)

    INTEGER                                  :: numspec, numspec2, ispec(*), &
                                                miller(3), iatom(*)
    REAL(real_8)                             :: tscr(3,*)
    INTEGER                                  :: natom
    REAL(real_8)                             :: an(*)
    INTEGER                                  :: lsk(3,natom)
    REAL(real_8)                             :: scalp
    INTEGER                                  :: SPvsAT

    CHARACTER(*), PARAMETER                  :: procedureN = 'disp_lnm'
    INTEGER, PARAMETER                       :: iunit = 58 

    CHARACTER(len=100)                       :: filen
    CHARACTER(len=80)                        :: element, label
    INTEGER :: i, ia, iat, ib, icount, iend, ierr, ii, iini, is, isp, ityp, &
      j, jj, jx, jy, jz, k, kx, l1, l2, l3, liflag, naa, nat_disp, nat_disp2
    INTEGER, ALLOCATABLE                     :: iflag(:)
    INTEGER, ALLOCATABLE, SAVE               :: iiindex(:), ipos_disp(:), &
                                                ipos_disp2(:)
    INTEGER, SAVE                            :: ifirst = 0
    LOGICAL                                  :: ferror, trotmin
    REAL(real_8) :: dir_norm, dir_vec(3), disp_vec(3), eps, fff, fx, fy, fz, &
      rax, ray, raz, rot(3,3), rx, ry, rz, scalp2, transl(3), unat_dis, &
      unat_dis2, xt, xx(3), yt, zt
    REAL(real_8), ALLOCATABLE, SAVE          :: d2_a(:,:), der_rot(:,:,:,:), &
                                                derr(:,:), pos_a(:,:), &
                                                pos_b(:,:), pos_bb(:,:)
    REAL(real_8), DIMENSION(3), SAVE         :: transl0

    IF (isos1%tisos) THEN
       trotmin = .TRUE.
    ENDIF
    IF (ifirst.EQ.0) THEN
       ! C     Allocate memory
       ALLOCATE(pos_a(3,natom),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       IF (duat%ndat .GT. 0) THEN
          ALLOCATE(d2_a(3,duat%ndat),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
       ALLOCATE(ipos_disp((natom+duat%ndat)),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(ipos_disp((natom+duat%ndat)),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(iiindex(natom),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       IF (trotmin) THEN
          ALLOCATE(pos_b(3,natom),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(pos_b(3,natom),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(derr(3,natom),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(der_rot(3,3,3,natom),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF

       ! C     Read the reference configurations
       filen='STRUCTURE_A'
       ferror=.FALSE.
       IF (paral%io_parent)&
            CALL fileopen(iunit,filen,fo_old,ferror)
       IF (ferror) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A)')&
               ' DISP_lmn! File STRUCTURE_A  DOES NOT EXIST!!'
          CALL dcopy(3*natom,tscr,1,pos_a,1)
          IF (paral%io_parent)&
               CALL fileopen(iunit,filen,fo_new,ferror)
          IF (paral%io_parent)&
               WRITE(iunit,'(A,A)') 'INITIAL',' POSITIONS'
          iat = 0
          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                iat = iat+1
                IF (paral%io_parent)&
                     WRITE(iunit,'(A,3f16.8)') elem%el(ions0%iatyp(is)),&
                     (tscr(k,iat),k=1,3)
             ENDDO
          ENDDO
       ELSE
          IF (paral%io_parent)&
               REWIND(iunit)
          IF (paral%io_parent)&
               READ(iunit,err=998,fmt=*) label
          iat = 0
          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                iat = iat+1
                IF (paral%io_parent)&
                     READ(iunit,err=998,fmt=*) element,&
                     (pos_a(k,iat),k=1,3)
             ENDDO
          ENDDO
          IF (paral%io_parent)&
               WRITE(6,'(/,A,A,/)') ' A COORDINATES READ ',&
               'FROM FILE STRUCTURE_A'
       ENDIF

       ! C            write(6,*) NUMSPEC,NUMSPEC2, MILLER(1),MILLER(2),MILLER(3)
       ! C            write(6,*)  (IATOM(K),K= 1,NUMSPEC)
       ! C            write(6,*)  (IATOM(K+NUMSPEC2),K= 1,NUMSPEC2)
       ! C            stop 'disp'

       xt = 0.0_real_8
       yt = 0.0_real_8
       zt = 0.0_real_8
       DO iat = 1,natom
          xt = xt + pos_a(1,iat)
          yt = yt + pos_a(2,iat)
          zt = zt + pos_a(3,iat)
       ENDDO
       transl0(1) = xt/REAL(natom,kind=real_8)
       transl0(2) = yt/REAL(natom,kind=real_8)
       transl0(3) = zt/REAL(natom,kind=real_8)
       !$omp parallel do private(IAT)
       DO iat = 1,natom
          pos_a(1,iat)= pos_a(1,iat) - transl0(1)
          pos_a(2,iat)= pos_a(2,iat) - transl0(2)
          pos_a(3,iat)= pos_a(3,iat) - transl0(3)
       ENDDO

       ifirst=1
       IF (paral%io_parent)&
            CALL fileclose(iunit)
       IF (duat%ndat2.EQ.0) THEN
          GOTO 999
       ENDIF
       DO i=1,duat%ndat
          IF (duat%listda(i,1).EQ.2) THEN
             kx=duat%listda(i,2)
             jj=duat%listd2(1,kx)
             xt=0.0_real_8
             yt=0.0_real_8
             zt=0.0_real_8
             DO k=1,jj
                xt=xt+pos_a(1,duat%listd2(k+1,kx))
                yt=yt+pos_a(2,duat%listd2(k+1,kx))
                zt=zt+pos_a(3,duat%listd2(k+1,kx))
             ENDDO
             d2_a(1,kx)=xt/REAL(jj,kind=real_8)
             d2_a(2,kx)=yt/REAL(jj,kind=real_8)
             d2_a(3,kx)=zt/REAL(jj,kind=real_8)
          ENDIF
       ENDDO
       GOTO 999

998    CONTINUE
       IF (paral%io_parent)&
            WRITE(6,'(A)') ' DISP_lmn! ERROR WHILE READING COORDINATES'
       CALL stopgm('DISP_lmn','READING ERROR',& 
            __LINE__,__FILE__)

999    CONTINUE
    ENDIF

    ! C     Vector that define a direction in space
    DO k = 1,3
       dir_vec(k) = REAL(miller(1),kind=real_8)*parm%a1(k)+REAL(miller(2),kind=real_8)*parm%a2(k)+&
            REAL(miller(3),kind=real_8)*parm%a3(k)
    ENDDO
    dir_norm = SQRT(dir_vec(1)*dir_vec(1)+&
         dir_vec(2)*dir_vec(2)+&
         dir_vec(3)*dir_vec(3))
    dir_vec(1) = dir_vec(1)/dir_norm
    dir_vec(2) = dir_vec(2)/dir_norm
    dir_vec(3) = dir_vec(3)/dir_norm


    ! C      write(6,*) ' dir ',DIR_VEC(1),DIR_VEC(2),DIR_VEC(3)
    ! C      call rmsd_rot(NATOM,POS_A,TSCR,transl,rot)


    IF (SPvsAT .EQ. 1) THEN
       ALLOCATE(iflag(ions1%nsp),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       liflag = ions1%nsp
    ELSE
       liflag = natom+duat%ndat
       ALLOCATE(iflag(liflag),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF

    CALL zeroing(ipos_disp)!,natom+duat%ndat)
    CALL zeroing(ipos_disp2)!,natom+duat%ndat)
    CALL zeroing(iflag)!,liflag)
    CALL zeroing(disp_vec)!,3)
    CALL zeroing(iiindex)!,natom)
    IF (trotmin) THEN
       CALL zeroing(derr)!,3*natom)
       CALL zeroing(der_rot)!,3*3*3*natom)
    ENDIF
    ! C     Create the list of engaged ions

    IF (SPvsAT .EQ. 1) THEN
       icount = 0
       nat_disp = 0
       DO is=1,numspec
          isp = ispec(is)
          iflag(isp) = 1
          iini = 0
          DO i = 1,isp-1
             iini = iini + ions0%na(i)
          ENDDO
          iend = iini + ions0%na(isp)
          iini = iini + 1
          DO iat = iini,iend
             nat_disp = nat_disp + 1
             ipos_disp(nat_disp) = iat

             IF (trotmin) THEN
                icount = icount + 1
                iiindex(icount) = iat
                pos_b(1,icount) = pos_a(1,iat)
                pos_b(2,icount) = pos_a(2,iat)
                pos_b(3,icount) = pos_a(3,iat)
                pos_bb(1,icount) = tscr(1,iat)
                pos_bb(2,icount) = tscr(2,iat)
                pos_bb(3,icount) = tscr(3,iat)
             ENDIF
          ENDDO
       ENDDO
       nat_disp2 = 0
       IF (numspec2 .EQ. 0) THEN
          iini = 0
          iend = 0
          DO is=1,ions1%nsp
             iend = iini + ions0%na(is)
             iini = iini + 1
             IF (iflag(is) .EQ. 0) THEN
                DO iat = iini,iend
                   nat_disp2 = nat_disp2 + 1
                   ipos_disp2(nat_disp2) = iat
                   IF (trotmin) THEN
                      icount = icount + 1
                      iiindex(icount) = iat
                      pos_b(1,icount) = pos_a(1,iat)
                      pos_b(2,icount) = pos_a(2,iat)
                      pos_b(3,icount) = pos_a(3,iat)
                      pos_bb(1,icount) = tscr(1,iat)
                      pos_bb(2,icount) = tscr(2,iat)
                      pos_bb(3,icount) = tscr(3,iat)
                   ENDIF
                ENDDO
             ENDIF
             iini = iend
          ENDDO
       ELSE
          DO is=1,numspec2
             isp = ispec(numspec+is)
             iini = 0
             DO i = 1,isp-1
                iini = iini + ions0%na(i)
             ENDDO
             iend = iini + ions0%na(isp)
             iini = iini + 1
             DO iat = iini,iend
                nat_disp2 = nat_disp2 + 1
                ipos_disp2(nat_disp2) = iat
                IF (trotmin) THEN
                   icount = icount + 1
                   iiindex(icount) = iat
                   pos_b(1,icount) = pos_a(1,iat)
                   pos_b(2,icount) = pos_a(2,iat)
                   pos_b(3,icount) = pos_a(3,iat)
                   pos_bb(1,icount) = tscr(1,iat)
                   pos_bb(2,icount) = tscr(2,iat)
                   pos_bb(3,icount) = tscr(3,iat)
                ENDIF
             ENDDO
          ENDDO
       ENDIF
    ELSE
       icount = 0
       nat_disp = 0
       DO is=1,numspec
          iat = iatom(is)

          nat_disp = nat_disp + 1
          ipos_disp(nat_disp) = iat

          IF (trotmin) THEN
             IF ( iat .LE. natom .AND. iflag(iat) .EQ. 0 ) THEN
                icount = icount + 1
                iiindex(icount) = iat
                pos_b(1,icount) = pos_a(1,iat)
                pos_b(2,icount) = pos_a(2,iat)
                pos_b(3,icount) = pos_a(3,iat)
                pos_bb(1,icount) = tscr(1,iat)
                pos_bb(2,icount) = tscr(2,iat)
                pos_bb(3,icount) = tscr(3,iat)
             ELSEIF (iat .LE. natom + duat%ndat) THEN
                naa=iat-ions1%nat
                ityp=duat%listda(naa,1)
                IF (ityp.EQ.1) THEN
                   GOTO 100
                ELSEIF (ityp.EQ.2) THEN
                   kx=duat%listda(naa,2)
                   jj=duat%listd2(1,kx)
                   DO j=1,jj
                      ib=duat%listd2(j+1,kx)
                      IF (iflag(ib) .EQ. 0) THEN
                         icount = icount + 1
                         iiindex(icount) = ib
                         pos_b(1,icount) = pos_a(1,ib)
                         pos_b(2,icount) = pos_a(2,ib)
                         pos_b(3,icount) = pos_a(3,ib)
                         pos_bb(1,icount) = tscr(1,ib)
                         pos_bb(2,icount) = tscr(2,ib)
                         pos_bb(3,icount) = tscr(3,ib)
                         iflag(ib) = 1
                      ENDIF
                   ENDDO
                ENDIF
             ENDIF
          ENDIF
100       CONTINUE
          iflag(iat) = 1
       ENDDO
       nat_disp2 = 0
       IF (numspec2 .EQ. 0) THEN
          DO is=1,natom
             iat = is
             IF (iflag(iat) .EQ. 0) THEN
                nat_disp2 = nat_disp2 + 1
                ipos_disp2(nat_disp2) = iat
                IF (trotmin) THEN
                   icount = icount + 1
                   iiindex(icount) = iat
                   pos_b(1,icount) = pos_a(1,iat)
                   pos_b(2,icount) = pos_a(2,iat)
                   pos_b(3,icount) = pos_a(3,iat)
                   pos_bb(1,icount) = tscr(1,iat)
                   pos_bb(2,icount) = tscr(2,iat)
                   pos_bb(3,icount) = tscr(3,iat)
                ENDIF
             ENDIF
          ENDDO
       ELSE
          DO is=1,numspec2
             iat = iatom(numspec+is)
             nat_disp2 = nat_disp2 + 1
             ipos_disp2(nat_disp2) = iat
             IF (trotmin) THEN
                IF ( iat .LE. natom .AND. iflag(iat) .EQ. 0 ) THEN
                   icount = icount + 1
                   iiindex(icount) = iat
                   pos_b(1,icount) = pos_a(1,iat)
                   pos_b(2,icount) = pos_a(2,iat)
                   pos_b(3,icount) = pos_a(3,iat)
                   pos_bb(1,icount) = tscr(1,iat)
                   pos_bb(2,icount) = tscr(2,iat)
                   pos_bb(3,icount) = tscr(3,iat)
                ELSEIF (iat .LE. natom + duat%ndat) THEN
                   naa=iat-ions1%nat
                   ityp=duat%listda(naa,1)
                   IF (ityp.EQ.1) THEN
                      GOTO 200
                   ELSEIF (ityp.EQ.2) THEN
                      kx=duat%listda(naa,2)
                      jj=duat%listd2(1,kx)
                      DO j=1,jj
                         ib=duat%listd2(j+1,kx)
                         IF (iflag(ib) .EQ. 0) THEN
                            icount = icount + 1
                            iiindex(icount) = ib
                            pos_b(1,icount) = pos_a(1,ib)
                            pos_b(2,icount) = pos_a(2,ib)
                            pos_b(3,icount) = pos_a(3,ib)
                            pos_bb(1,icount) = tscr(1,ib)
                            pos_bb(2,icount) = tscr(2,ib)
                            pos_bb(3,icount) = tscr(3,ib)
                            iflag(ib) = 1
                         ENDIF
                      ENDDO
                   ENDIF
                ENDIF
             ENDIF
200          CONTINUE
             iflag(iat) = 1
          ENDDO
       ENDIF
    ENDIF

    IF (trotmin) THEN
       eps = 0.001_real_8
       CALL rmsd_rot(icount,pos_bb,pos_b,iiindex,transl,rot,&
            natom,derr,der_rot,.FALSE.)
    ELSE
       CALL zeroing(rot)!,9)
       DO i = 1,3
          rot(i,i) = 1.0_real_8
          xt = 0.0_real_8
          DO iat = 1,natom
             xt = xt + tscr(i,iat)
          ENDDO
          transl(i) = xt/REAL(natom,kind=real_8)
       ENDDO
    ENDIF

    unat_dis  = 1.0_real_8/REAL(nat_disp,kind=real_8)
    unat_dis2 = 1.0_real_8/REAL(nat_disp2,kind=real_8)

    ! write(6,*) 'n na ',NAT_DISP,NAT_DISP2
    DO iat = 1,nat_disp
       ii = ipos_disp(iat)
       CALL fillc(ii,tscr,xx)
       xt = xx(1) -transl(1)
       yt = xx(2) -transl(2)
       zt = xx(3) -transl(3)
       rx =rot(1,1)*xt + rot(1,2)*yt + rot(1,3)*zt
       ry =rot(2,1)*xt + rot(2,2)*yt + rot(2,3)*zt
       rz =rot(3,1)*xt + rot(3,2)*yt + rot(3,3)*zt
       CALL fillc2(ii,pos_a,xx,d2_a,transl0)
       rax = xx(1)
       ray = xx(2)
       raz = xx(3)
       disp_vec(1)  = disp_vec(1) + rx - rax
       disp_vec(2)  = disp_vec(2) + ry - ray
       disp_vec(3)  = disp_vec(3) + rz - raz
    ENDDO
    ! ==--------------------------------------------------------------==
    ! Calculate scalar product: DISP_VEC . DIR_VEC
    scalp = disp_vec(1)*dir_vec(1) +  disp_vec(2)*dir_vec(2) +&
         disp_vec(3)*dir_vec(3)

    scalp = scalp * unat_dis

    CALL zeroing(disp_vec)!,3)
    ! write(6,*) 'scalp ',SCALP
    ! write(6,'(A,4f12.6)') 'first ',
    ! &      DISP_VEC(1),DISP_VEC(2),DISP_VEC(3),SCALP
    ! ==--------------------------------------------------------------==

    DO iat = 1,nat_disp2

       ii = ipos_disp2(iat)
       CALL fillc(ii,tscr,xx)
       xt = xx(1) -transl(1)
       yt = xx(2) -transl(2)
       zt = xx(3) -transl(3)
       rx =rot(1,1)*xt + rot(1,2)*yt + rot(1,3)*zt
       ry =rot(2,1)*xt + rot(2,2)*yt + rot(2,3)*zt
       rz =rot(3,1)*xt + rot(3,2)*yt + rot(3,3)*zt

       CALL fillc2(ii,pos_a,xx,d2_a,transl0)
       rax = xx(1)
       ray = xx(2)
       raz = xx(3)

       disp_vec(1)   = disp_vec(1) + rx - rax
       disp_vec(2)   = disp_vec(2) + ry - ray
       disp_vec(3)   = disp_vec(3) + rz - raz
    ENDDO
    ! ==--------------------------------------------------------------==
    ! Calculate scalar product: DISP_VEC . DIR_VEC
    scalp2 = disp_vec(1)*dir_vec(1) +  disp_vec(2)*dir_vec(2) +&
         disp_vec(3)*dir_vec(3)

    scalp2 = scalp2 * unat_dis2
    ! write(6,*) 'scalp2 ',SCALP2
    ! write(6,'(A,4f12.6)') 'second ',
    ! &      DISP_VEC(1),DISP_VEC(2),DISP_VEC(3),SCALP2
    ! ==--------------------------------------------------------------==
    icount = 0

    DO iat = 1,nat_disp
       ii = ipos_disp(iat)
       IF (ii .LE. natom) THEN
          l1 = lsk(1,ii)
          l2 = lsk(2,ii)
          l3 = lsk(3,ii)

          IF (trotmin) THEN
             xt = tscr(1,ii) -transl(1)
             yt = tscr(2,ii) -transl(2)
             zt = tscr(3,ii) -transl(3)
             IF (l1.NE.0) an(l1) = an(l1) +&
                  unat_dis * ( rot(1,1)*dir_vec(1) +&
                  rot(2,1)*dir_vec(2) + rot(3,1)*dir_vec(3) +&
                  (der_rot(1,1,1,ii)*xt+der_rot(1,2,1,ii)*yt +&
                  der_rot(1,3,1,ii)*zt) * dir_vec(1) +&
                  (der_rot(2,1,1,ii)*xt+der_rot(2,2,1,ii)*yt +&
                  der_rot(2,3,1,ii)*zt) * dir_vec(2) +&
                  (der_rot(3,1,1,ii)*xt+der_rot(3,2,1,ii)*yt +&
                  der_rot(3,3,1,ii)*zt) * dir_vec(3) )

             IF (l2.NE.0) an(l2) = an(l2) +&
                  unat_dis * ( rot(1,2)*dir_vec(1) +&
                  rot(2,2)*dir_vec(2) + rot(3,2)*dir_vec(3) +&
                  (der_rot(1,1,2,ii)*xt+der_rot(1,2,2,ii)*yt +&
                  der_rot(1,3,2,ii)*zt) * dir_vec(1) +&
                  (der_rot(2,1,2,ii)*xt+der_rot(2,2,2,ii)*yt +&
                  der_rot(2,3,2,ii)*zt) * dir_vec(2) +&
                  (der_rot(3,1,2,ii)*xt+der_rot(3,2,2,ii)*yt +&
                  der_rot(3,3,2,ii)*zt) * dir_vec(3))

             IF (l3.NE.0) an(l3) = an(l3) +&
                  unat_dis * ( rot(1,3)*dir_vec(1) +&
                  rot(2,3)*dir_vec(2) + rot(3,3)*dir_vec(3) +&
                  (der_rot(1,1,3,ii)*xt+der_rot(1,2,3,ii)*yt +&
                  der_rot(1,3,3,ii)*zt) * dir_vec(1) +&
                  (der_rot(2,1,3,ii)*xt+der_rot(2,2,3,ii)*yt +&
                  der_rot(2,3,3,ii)*zt) * dir_vec(2) +&
                  (der_rot(3,1,3,ii)*xt+der_rot(3,2,3,ii)*yt +&
                  der_rot(3,3,3,ii)*zt) * dir_vec(3) )

          ELSE
             IF (l1.NE.0) an(l1) = an(l1) + unat_dis * dir_vec(1)
             IF (l2.NE.0) an(l2) = an(l2) + unat_dis * dir_vec(2)
             IF (l3.NE.0) an(l3) = an(l3) + unat_dis * dir_vec(3)
          ENDIF
       ELSEIF (ii .LE. natom+duat%ndat) THEN! Dummy
          naa=ii-ions1%nat
          ityp=duat%listda(naa,1)
          IF (ityp.EQ.1) THEN
             GOTO 110
          ELSEIF (ityp.EQ.2) THEN
             kx=duat%listda(naa,2)
             jj=duat%listd2(1,kx)
             fff=1.0_real_8/REAL(jj,kind=real_8)
             DO j=1,jj
                ib=duat%listd2(j+1,kx)
                jx=lskptr(1,ib)
                jy=lskptr(2,ib)
                jz=lskptr(3,ib)
                IF (trotmin) THEN
                   fx = unat_dis * ( rot(1,1)*dir_vec(1) +&
                        rot(2,1)*dir_vec(2) + rot(3,1)*dir_vec(3) +&
                        (der_rot(1,1,1,ib)*xt+der_rot(1,2,1,ib)*yt +&
                        der_rot(1,3,1,ib)*zt) * dir_vec(1) +&
                        (der_rot(2,1,1,ib)*xt+der_rot(2,2,1,ib)*yt +&
                        der_rot(2,3,1,ib)*zt) * dir_vec(2) +&
                        (der_rot(3,1,1,ib)*xt+der_rot(3,2,1,ib)*yt +&
                        der_rot(3,3,1,ib)*zt) * dir_vec(3) )
                   fy =   unat_dis * ( rot(1,2)*dir_vec(1) +&
                        rot(2,2)*dir_vec(2) + rot(3,2)*dir_vec(3) +&
                        (der_rot(1,1,2,ib)*xt+der_rot(1,2,2,ib)*yt +&
                        der_rot(1,3,2,ib)*zt) * dir_vec(1) +&
                        (der_rot(2,1,2,ib)*xt+der_rot(2,2,2,ib)*yt +&
                        der_rot(2,3,2,ib)*zt) * dir_vec(2) +&
                        (der_rot(3,1,2,ib)*xt+der_rot(3,2,2,ib)*yt +&
                        der_rot(3,3,2,ib)*zt) * dir_vec(3))
                   fz =   unat_dis * ( rot(1,3)*dir_vec(1) +&
                        rot(2,3)*dir_vec(2) + rot(3,3)*dir_vec(3) +&
                        (der_rot(1,1,3,ib)*xt+der_rot(1,2,3,ib)*yt +&
                        der_rot(1,3,3,ib)*zt) * dir_vec(1) +&
                        (der_rot(2,1,3,ib)*xt+der_rot(2,2,3,ib)*yt +&
                        der_rot(2,3,3,ib)*zt) * dir_vec(2) +&
                        (der_rot(3,1,3,ib)*xt+der_rot(3,2,3,ib)*yt +&
                        der_rot(3,3,3,ib)*zt) * dir_vec(3) )
                ELSE
                   fx = unat_dis * dir_vec(1)
                   fy = unat_dis * dir_vec(2)
                   fz = unat_dis * dir_vec(3)
                ENDIF
                IF (jx.NE.0) THEN
                   an(jx)=an(jx)+fx*fff
                ENDIF
                IF (jy.NE.0) THEN
                   an(jy)=an(jy)+fy*fff
                ENDIF
                IF (jz.NE.0) THEN
                   an(jz)=an(jz)+fz*fff
                ENDIF
             ENDDO
          ENDIF
       ENDIF
110    CONTINUE
    ENDDO

    DO iat = 1,nat_disp2
       ii = ipos_disp2(iat)
       IF (ii .LE. natom) THEN
          l1 = lsk(1,ipos_disp2(iat))
          l2 = lsk(2,ipos_disp2(iat))
          l3 = lsk(3,ipos_disp2(iat))

          IF (trotmin) THEN
             xt = tscr(1,ii) -transl(1)
             yt = tscr(2,ii) -transl(2)
             zt = tscr(3,ii) -transl(3)
             IF (l1.NE.0) THEN
                an(l1) = an(l1) -&
                     unat_dis2 * ( rot(1,1)*dir_vec(1) +&
                     rot(2,1)*dir_vec(2) + rot(3,1)*dir_vec(3) +&
                     (der_rot(1,1,1,ii)*xt+der_rot(1,2,1,ii)*yt +&
                     der_rot(1,3,1,ii)*zt) * dir_vec(1) +&
                     (der_rot(2,1,1,ii)*xt+der_rot(2,2,1,ii)*yt +&
                     der_rot(2,3,1,ii)*zt) * dir_vec(2) +&
                     (der_rot(3,1,1,ii)*xt+der_rot(3,2,1,ii)*yt +&
                     der_rot(3,3,1,ii)*zt) * dir_vec(3) )
             ENDIF
             IF (l2.NE.0) THEN
                an(l2) = an(l2) -&
                     unat_dis2 * ( rot(1,2)*dir_vec(1) +&
                     rot(2,2)*dir_vec(2) + rot(3,2)*dir_vec(3) +&
                     (der_rot(1,1,2,ii)*xt+der_rot(1,2,2,ii)*yt +&
                     der_rot(1,3,2,ii)*zt) * dir_vec(1) +&
                     (der_rot(2,1,2,ii)*xt+der_rot(2,2,2,ii)*yt +&
                     der_rot(2,3,2,ii)*zt) * dir_vec(2) +&
                     (der_rot(3,1,2,ii)*xt+der_rot(3,2,2,ii)*yt +&
                     der_rot(3,3,2,ii)*zt) * dir_vec(3))
             ENDIF
             IF (l3.NE.0) THEN
                an(l3) = an(l3) -&
                     unat_dis2 * ( rot(1,3)*dir_vec(1) +&
                     rot(2,3)*dir_vec(2) + rot(3,3)*dir_vec(3) +&
                     (der_rot(1,1,3,ii)*xt+der_rot(1,2,3,ii)*yt +&
                     der_rot(1,3,3,ii)*zt) * dir_vec(1) +&
                     (der_rot(2,1,3,ii)*xt+der_rot(2,2,3,ii)*yt +&
                     der_rot(2,3,3,ii)*zt) * dir_vec(2) +&
                     (der_rot(3,1,3,ii)*xt+der_rot(3,2,3,ii)*yt +&
                     der_rot(3,3,3,ii)*zt) * dir_vec(3) )
             ENDIF
          ELSE
             IF (l1.NE.0) THEN
                an(l1) = an(l1) - unat_dis2 * dir_vec(1)
             ENDIF
             IF (l2.NE.0) THEN
                an(l2) = an(l2) - unat_dis2 * dir_vec(2)
             ENDIF
             IF (l3.NE.0) THEN
                an(l3) = an(l3) - unat_dis2 * dir_vec(3)
             ENDIF
          ENDIF
       ELSEIF (ii .LE. natom+duat%ndat) THEN! Dummy
          naa=ii-ions1%nat
          ityp=duat%listda(naa,1)
          IF (ityp.EQ.1) THEN
             GOTO 210
          ELSEIF (ityp.EQ.2) THEN
             kx=duat%listda(naa,2)
             jj=duat%listd2(1,kx)
             fff=1.0_real_8/REAL(jj,kind=real_8)
             DO j=1,jj
                ib=duat%listd2(j+1,kx)
                jx=lskptr(1,ib)
                jy=lskptr(2,ib)
                jz=lskptr(3,ib)
                IF (trotmin) THEN
                   fx = -unat_dis2 * ( rot(1,1)*dir_vec(1) +&
                        rot(2,1)*dir_vec(2) + rot(3,1)*dir_vec(3) +&
                        (der_rot(1,1,1,ib)*xt+der_rot(1,2,1,ib)*yt +&
                        der_rot(1,3,1,ib)*zt) * dir_vec(1) +&
                        (der_rot(2,1,1,ib)*xt+der_rot(2,2,1,ib)*yt +&
                        der_rot(2,3,1,ib)*zt) * dir_vec(2) +&
                        (der_rot(3,1,1,ib)*xt+der_rot(3,2,1,ib)*yt +&
                        der_rot(3,3,1,ib)*zt) * dir_vec(3) )
                   fy =  -unat_dis2 * ( rot(1,2)*dir_vec(1) +&
                        rot(2,2)*dir_vec(2) + rot(3,2)*dir_vec(3) +&
                        (der_rot(1,1,2,ib)*xt+der_rot(1,2,2,ib)*yt +&
                        der_rot(1,3,2,ib)*zt) * dir_vec(1) +&
                        (der_rot(2,1,2,ib)*xt+der_rot(2,2,2,ib)*yt +&
                        der_rot(2,3,2,ib)*zt) * dir_vec(2) +&
                        (der_rot(3,1,2,ib)*xt+der_rot(3,2,2,ib)*yt +&
                        der_rot(3,3,2,ib)*zt) * dir_vec(3))
                   fz =  -unat_dis2 * ( rot(1,3)*dir_vec(1) +&
                        rot(2,3)*dir_vec(2) + rot(3,3)*dir_vec(3) +&
                        (der_rot(1,1,3,ib)*xt+der_rot(1,2,3,ib)*yt +&
                        der_rot(1,3,3,ib)*zt) * dir_vec(1) +&
                        (der_rot(2,1,3,ib)*xt+der_rot(2,2,3,ib)*yt +&
                        der_rot(2,3,3,ib)*zt) * dir_vec(2) +&
                        (der_rot(3,1,3,ib)*xt+der_rot(3,2,3,ib)*yt +&
                        der_rot(3,3,3,ib)*zt) * dir_vec(3) )
                ELSE
                   fx = -unat_dis2 * dir_vec(1)
                   fy = -unat_dis2 * dir_vec(2)
                   fz = -unat_dis2 * dir_vec(3)
                ENDIF
                IF (jx.NE.0) THEN
                   an(jx)=an(jx)+fx*fff
                ENDIF
                IF (jy.NE.0) THEN
                   an(jy)=an(jy)+fy*fff
                ENDIF
                IF (jz.NE.0) THEN
                   an(jz)=an(jz)+fz*fff
                ENDIF
             ENDDO
          ENDIF
       ENDIF
210    CONTINUE
    ENDDO

    scalp = scalp-scalp2

    DEALLOCATE(iflag,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! write (6,*) scalp
    ! stop 'scalp'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE disp_lnm
  ! ==================================================================
  SUBROUTINE rmsd_rot(n,r,r0,index,transl,rot,natom,&
       derr_dr,der_rot,pbc)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: n
    REAL(real_8)                             :: r(3,n), r0(3,n)
    INTEGER                                  :: INDEX(n)
    REAL(real_8)                             :: transl(3), rot(3,3)
    INTEGER                                  :: natom
    REAL(real_8)                             :: derr_dr(3,natom), &
                                                der_rot(3,3,3,natom)
    LOGICAL                                  :: pbc

    REAL(real_8), PARAMETER                  :: epsi = 1.0e-10_real_8

    CHARACTER(len=80)                        :: msg
    INTEGER                                  :: i, ier, iq, iret, ix, j, k
    REAL(real_8) :: der_q(0:3), dm_r(4,4,3), dq(0:3), err, lambda(4), m(4,4), &
      q(0:3), r0p(3,n), rp(3,n), rr(3), rr0(3), rrsq, s, wk(20), xx, yy, &
      z(4,4), zz

! ==--------------------------------------------------------------==
! center r

    xx=0._real_8
    yy=0._real_8
    zz=0._real_8
    !$omp parallel do private(i) reduction(+:xx,yy,zz)
    DO i=1,n
       xx=xx+r(1,i)
       yy=yy+r(2,i)
       zz=zz+r(3,i)
    ENDDO
    xx=xx/REAL(n,kind=real_8)
    yy=yy/REAL(n,kind=real_8)
    zz=zz/REAL(n,kind=real_8)
    !$omp parallel do private(i) shared(xx,yy,zz)
    DO i=1,n
       rp(1,i)=r(1,i)-xx
       rp(2,i)=r(2,i)-yy
       rp(3,i)=r(3,i)-zz
    ENDDO
    transl(1)=xx
    transl(2)=yy
    transl(3)=zz
    ! center r0
    xx=0._real_8
    yy=0._real_8
    zz=0._real_8
    !$omp parallel do private(i) reduction(+:xx,yy,zz)
    DO i=1,n
       xx=xx+r0(1,i)
       yy=yy+r0(2,i)
       zz=zz+r0(3,i)
    ENDDO
    xx=xx/REAL(n,kind=real_8)
    yy=yy/REAL(n,kind=real_8)
    zz=zz/REAL(n,kind=real_8)
    !$omp parallel do private(i) shared(xx,yy,zz)
    DO i=1,n
       r0p(1,i)=r0(1,i)-xx
       r0p(2,i)=r0(2,i)-yy
       r0p(3,i)=r0(3,i)-zz
    ENDDO
    transl(1)=xx
    transl(2)=yy
    transl(3)=zz
    ! 
    ! do i = 1,16
    ! m(i,1) = 0._real_8
    ! enddo
    CALL zeroing(m)!,4*4)

    DO i = 1,n
       ! if (w(i) .eq. 0._real_8) cycle
       rr(1)=rp(1,i)
       rr(2)=rp(2,i)
       rr(3)=rp(3,i)
       rr0(1)=r0p(1,i)
       rr0(2)=r0p(2,i)
       rr0(3)=r0p(3,i)
       rrsq=rr0(1)**2+rr0(2)**2+rr0(3)**2+rr(1)**2+rr(2)**2+rr(3)**2
       ! rr0(1)=w(i)*rr0(1)
       ! rr0(2)=w(i)*rr0(2)
       ! rr0(3)=w(i)*rr0(3)

       m(1,1)=m(1,1) + rrsq+2._real_8*&
            (-rr0(1)*rr(1)-rr0(2)*rr(2)-rr0(3)*rr(3))
       m(2,2)=m(2,2) + rrsq+2._real_8*&
            (-rr0(1)*rr(1)+rr0(2)*rr(2)+rr0(3)*rr(3))
       m(3,3)=m(3,3) + rrsq+2._real_8*&
            (+rr0(1)*rr(1)-rr0(2)*rr(2)+rr0(3)*rr(3))
       m(4,4)=m(4,4) + rrsq+2._real_8*&
            (+rr0(1)*rr(1)+rr0(2)*rr(2)-rr0(3)*rr(3))
       m(1,2)=m(1,2) +2._real_8*(-rr0(2)*rr(3)+rr0(3)*rr(2))
       m(1,3)=m(1,3) +2._real_8*(rr0(1)*rr(3)-rr0(3)*rr(1))
       m(1,4)=m(1,4) +2._real_8*(-rr0(1)*rr(2)+rr0(2)*rr(1))
       m(2,3)=m(2,3) -2._real_8*(rr0(1)*rr(2)+rr0(2)*rr(1))
       m(2,4)=m(2,4) -2._real_8*(rr0(1)*rr(3)+rr0(3)*rr(1))
       m(3,4)=m(3,4) -2._real_8*(rr0(2)*rr(3)+rr0(3)*rr(2))
    ENDDO
    m(2,1) = m(1,2)
    m(3,1) = m(1,3)
    m(3,2) = m(2,3)
    m(4,1) = m(1,4)
    m(4,2) = m(2,4)
    m(4,3) = m(3,4)


    ! ---- SOLVE THE EIGENVECTOR PROBLEM FOR M: --------------------------*

    iret = 0
    ! jobn   = 12
    ! call EIGRS(m,4,jobn,lambda,z,4,wk,ier)
    ier=0
    CALL dcopy(16,m,1,z,1)
    CALL dsyev('V','U',4,z,4,lambda,wk,20,ier)
    IF (ier .NE. 0) THEN
       CALL stopgm('RMSD_ROT','Fatal error in DSYSEV',& 
            __LINE__,__FILE__)
    ENDIF
    IF (wk(1) .GT. 1._real_8) THEN
       iret = 9
    ENDIF
    ! ---- PICK THE CORRECT EIGENVECTOR(S): ------------------------------*

    s =  1._real_8
    IF (z(1,1) .LT. 0._real_8) THEN
       s = -1._real_8
    ENDIF
    q(0) = s*z(1,1)
    q(1) = s*z(2,1)
    q(2) = s*z(3,1)
    q(3) = s*z(4,1)
    IF (ABS(lambda(1)) .LT. epsi) THEN
       err = 0._real_8
    ELSE
       ! err = SQRT(lambda(1))
       err = lambda(1)/REAL(n,kind=real_8)
    ENDIF
    IF (ABS(lambda(1) - lambda(2)) .LT. epsi)iret = iret + 10

    IF (iret .EQ. 0)  THEN
       msg = ' LSQQTN: Normal execution, unique solution'
    ENDIF
    IF (iret .EQ. 10) THEN
       msg = ' LSQQTN: Normal execution, non-unique solution'
    ENDIF
    IF (iret .EQ. 9) THEN
       msg = ' LSQQTN: Bad perform. in eigrs, unique solution'
    ENDIF
    IF (iret .EQ. 19) THEN
       msg = ' LSQQTN: Bad perform. in eigrs, non-unique solution'
    ENDIF
    IF (iret.NE.0) THEN
       IF (paral%io_parent)&
            WRITE(6,*) msg
    ENDIF

    ! ---- DERIVATIVES OF RMSD with respect to the positions
    DO i = 1,n
       rr(1)=2._real_8*rp(1,i)
       rr(2)=2._real_8*rp(2,i)
       rr(3)=2._real_8*rp(3,i)
       rr0(1)=2._real_8*r0p(1,i)
       rr0(2)=2._real_8*r0p(2,i)
       rr0(3)=2._real_8*r0p(3,i)
       ! 
       dm_r (1,1,1)=(rr(1)-rr0(1))
       dm_r (1,1,2)=(rr(2)-rr0(2))
       dm_r (1,1,3)=(rr(3)-rr0(3))
       ! 
       dm_r (1,2,1)=0._real_8
       dm_r (1,2,2)= rr0(3)
       dm_r (1,2,3)=-rr0(2)
       ! 
       dm_r (1,3,1)=-rr0(3)
       dm_r (1,3,2)= 0._real_8
       dm_r (1,3,3)= rr0(1)
       ! 
       dm_r (1,4,1)= rr0(2)
       dm_r (1,4,2)=-rr0(1)
       dm_r (1,4,3)= 0._real_8
       ! 
       dm_r (2,2,1)=(rr(1)-rr0(1))
       dm_r (2,2,2)=(rr(2)+rr0(2))
       dm_r (2,2,3)=(rr(3)+rr0(3))
       ! 
       dm_r (2,3,1)=-rr0(2)
       dm_r (2,3,2)=-rr0(1)
       dm_r (2,3,3)= 0._real_8
       ! 
       dm_r (2,4,1)=-rr0(3)
       dm_r (2,4,2)= 0._real_8
       dm_r (2,4,3)=-rr0(1)
       ! 
       dm_r (3,3,1)=(rr(1)+rr0(1))
       dm_r (3,3,2)=(rr(2)-rr0(2))
       dm_r (3,3,3)=(rr(3)+rr0(3))
       ! 
       dm_r (3,4,1)=0._real_8
       dm_r (3,4,2)=-rr0(3)
       dm_r (3,4,3)=-rr0(2)
       ! 
       dm_r (4,4,1)=(rr(1)+rr0(1))
       dm_r (4,4,2)=(rr(2)+rr0(2))
       dm_r (4,4,3)=(rr(3)-rr0(3))
       ! 
       DO ix=1,3
          dm_r(2,1,ix)=dm_r(1,2,ix)
          dm_r(3,1,ix)=dm_r(1,3,ix)
          dm_r(4,1,ix)=dm_r(1,4,ix)
          dm_r(3,2,ix)=dm_r(2,3,ix)
          dm_r(4,2,ix)=dm_r(2,4,ix)
          dm_r(4,3,ix)=dm_r(3,4,ix)
       ENDDO
       ! 
       DO ix=1,3
          derr_dr (ix,INDEX(i))=0._real_8
          DO k=1,4
             DO j=1,4
                derr_dr (ix,INDEX(i))=derr_dr (ix,INDEX(i)) +&
                     q(k-1)*q(j-1)*dm_r (j,k,ix)
             ENDDO
          ENDDO
          derr_dr (ix,INDEX(i))=derr_dr (ix,INDEX(i))/REAL(n,kind=real_8)
       ENDDO

       IF (pbc) GOTO 120

       DO ix = 1,3
          CALL zeroing(der_q)!(0),4)

          DO iq=2,4
             dq(iq-1) = 0.0_real_8
             DO j = 1,4
                DO k = 1,4
                   dq(iq-1) = dq(iq-1) +&
                        z(j,iq)*dm_r(j,k,ix)*z(k,1)
                ENDDO
             ENDDO
             dq(iq-1) = -dq(iq-1)/(lambda(iq)-lambda(1))
          ENDDO
          DO iq=1,4
             DO j = 1,4
                der_q(iq-1) = der_q(iq-1) + z(iq,j)*dq(j-1)*s
             ENDDO
          ENDDO

          der_rot(1,1,ix,INDEX(i)) = -4._real_8 * q(2)*der_q(2) -&
               4._real_8 * q(3)*der_q(3)
          der_rot(1,2,ix,INDEX(i)) =  2._real_8 *&
               ( -der_q(0)*q(3) - q(0)*der_q(3) +&
               der_q(1)*q(2) + q(1)*der_q(2) )
          der_rot(1,3,ix,INDEX(i)) =  2._real_8 *&
               (  der_q(0)*q(2) + q(0)*der_q(2) +&
               der_q(1)*q(3) + q(1)*der_q(3) )
          der_rot(2,1,ix,INDEX(i)) =  2._real_8 *&
               (  der_q(0)*q(3) + q(0)*der_q(3) +&
               der_q(1)*q(2) + q(1)*der_q(2) )
          der_rot(2,2,ix,INDEX(i)) = -4._real_8*q(1)*der_q(1) -&
               4._real_8*q(3)*der_q(3)
          der_rot(2,3,ix,INDEX(i)) =  2._real_8 *&
               ( -der_q(0)*q(1) - q(0)*der_q(1) +&
               der_q(2)*q(3) + q(2)*der_q(3) )
          der_rot(3,1,ix,INDEX(i)) =  2._real_8 *&
               ( -der_q(0)*q(2) - q(0)*der_q(2) +&
               der_q(1)*q(3) + q(1)*der_q(3) )
          der_rot(3,2,ix,INDEX(i)) =  2._real_8 *&
               (  der_q(0)*q(1) + q(0)*der_q(1) +&
               der_q(2)*q(3) + q(2)*der_q(3) )
          der_rot(3,3,ix,INDEX(i)) = -4._real_8*q(1)*der_q(1) -&
               4._real_8*q(2)*der_q(2)

       ENDDO

120    CONTINUE

    ENDDO

    ! ---- ROTATION MATRIX IN TERMS OF QUATERNIONS: ----------------------*

    rot(1,1)=-2._real_8*q(2)**2-2._real_8*q(3)**2+1._real_8
    rot(1,2)=2._real_8*(-q(0)*q(3)+q(1)*q(2))
    rot(1,3)=2._real_8*(q(0)*q(2)+q(1)*q(3))
    rot(2,1)=2._real_8*(q(0)*q(3)+q(1)*q(2))
    rot(2,2)=-2._real_8*q(1)**2-2._real_8*q(3)**2+1._real_8
    rot(2,3)=2._real_8*(-q(0)*q(1)+q(2)*q(3))
    rot(3,1)=2._real_8*(-q(0)*q(2)+q(1)*q(3))
    rot(3,2)=2._real_8*(q(0)*q(1)+q(2)*q(3))
    rot(3,3)=-2._real_8*q(1)**2-2._real_8*q(2)**2+1._real_8
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rmsd_rot
  ! ==================================================================
  SUBROUTINE funcp_mia(diff,teta,teta0,x,y,z,w,d)
    ! ==--------------------------------------------------------------==
    ! FUNCTION: phi - phi0
    ! Argument variables
    REAL(real_8)                             :: diff, teta, teta0, x(3), &
                                                y(3), z(3), w(3), d(12)

    REAL(real_8), PARAMETER                  :: epsilon = 1.e-12_real_8 

    REAL(real_8)                             :: const, coss, den, nn(3), nn2, &
                                                nn_norm, num, o(3), o2, &
                                                o_norm, pi4, q(3), scal, &
                                                sign0, sinn, t(3), t_norm

! cmb - original "sin" and "cos" conflict with fortran intrinsic fcn
! ==--------------------------------------------------------------==

    sign0=+1._real_8
    pi4=dacos(-1._real_8)
    t(1) = (x(1)-y(1))
    t(2) = (x(2)-y(2))
    t(3) = (x(3)-y(3))
    t_norm = SQRT(t(1)*t(1)+t(2)*t(2)+t(3)*t(3))
    q(1) = (z(1)-y(1))
    q(2) = (z(2)-y(2))
    q(3) = (z(3)-y(3))
    CALL vecprod(t,q,nn)
    nn2 = nn(1)*nn(1)+nn(2)*nn(2)+nn(3)*nn(3)
    nn_norm = SQRT(nn2)
    o(1) = (w(1)-y(1))
    o(2) = (w(2)-y(2))
    o(3) = (w(3)-y(3))
    o2 = o(1)*o(1)+o(2)*o(2)+o(3)*o(3)
    o_norm = SQRT(o2)

    CALL getscal(o,t,scal)
    CALL getscal(o,nn,coss)
    num = coss
    den = 1.0_real_8/(o_norm*nn_norm)
    coss = coss*den
    sinn = SQRT(1._real_8-coss*coss)
    teta = dacos(coss)
    IF (scal .LT. 0._real_8) THEN
       sign0=-1._real_8
       teta=-teta+2._real_8*pi4
    ENDIF

    diff=MOD(20000._real_8*pi4+teta0,2._real_8*pi4)-teta

    const = -1.0_real_8/sinn*sign0

    d(1) = const*den*(-q(3)*o(2)+q(2)*o(3)&
         -num*(-q(3)*nn(2)+q(2)*nn(3))/nn2)
    d(2) = const*den*(-q(1)*o(3)+q(3)*o(1)&
         -num*(-q(1)*nn(3)+q(3)*nn(1))/nn2)
    d(3) = const*den*(-q(2)*o(1)+q(1)*o(2)&
         -num*(-q(2)*nn(1)+q(1)*nn(2))/nn2)

    d(4) = const*den*(-nn(1)-t(3)*o(2)+q(3)*o(2)-q(2)*o(3)+t(2)*o(3)&
         -num*(-o(1)/o2+&
         ((-t(3)+q(3))*nn(2)+(-q(2)+t(2))*nn(3))/nn2))
    d(5) = const*den*(-nn(2)-q(3)*o(1)+t(3)*o(1)-t(1)*o(3)+q(1)*o(3)&
         -num*(-o(2)/o2+&
         ((-q(3)+t(3))*nn(1)+(-t(1)+q(1))*nn(3))/nn2))
    d(6) = const*den*(-nn(3)-t(2)*o(1)+q(2)*o(1)-q(1)*o(2)+t(1)*o(2)&
         -num*(-o(3)/o2+&
         ((-t(2)+q(2))*nn(1)+(-q(1)+t(1))*nn(2))/nn2))

    d(7) = const*den*(t(3)*o(2)-t(2)*o(3)&
         -num*(t(3)*nn(2)-t(2)*nn(3))/nn2)
    d(8) = const*den*(t(1)*o(3)-t(3)*o(1)&
         -num*(t(1)*nn(3)-t(3)*nn(1))/nn2)
    d(9) = const*den*(t(2)*o(1)-t(1)*o(2)&
         -num*(t(2)*nn(1)-t(1)*nn(2))/nn2)

    d(10) = const*den*(nn(1)-num*o(1)/o2)
    d(11) = const*den*(nn(2)-num*o(2)/o2)
    d(12) = const*den*(nn(3)-num*o(3)/o2)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE funcp_mia
  ! ==================================================================
  SUBROUTINE aplane(diff,teta,teta0,xa,ya,za,xb,yb,zb,d)
    ! ==--------------------------------------------------------------==
    ! FUNCTION: phi - phi0
    ! Argument variables
    REAL(real_8)                             :: diff, teta, teta0, xa(3), &
                                                ya(3), za(3), xb(3), yb(3), &
                                                zb(3), d(18)

    REAL(real_8)                             :: const, cosine, den, der1, &
                                                der2, der3, nn2a, nn2b, &
                                                nna(3), nnb(3), pi4, sign0, &
                                                sinus

! ==--------------------------------------------------------------==

    sign0=+1.0_real_8
    pi4=dacos(-1._real_8)

    nna(1) = (ya(2)-xa(2))*(za(3)-xa(3))-(ya(3)-xa(3))*(za(2)-xa(2))
    nna(2) = (ya(3)-xa(3))*(za(1)-xa(1))-(ya(1)-xa(1))*(za(3)-xa(3))
    nna(3) = (ya(1)-xa(1))*(za(2)-xa(2))-(ya(2)-xa(2))*(za(1)-xa(1))


    nnb(1) = (yb(2)-xb(2))*(zb(3)-xb(3))-(yb(3)-xb(3))*(zb(2)-xb(2))
    nnb(2) = (yb(3)-xb(3))*(zb(1)-xb(1))-(yb(1)-xb(1))*(zb(3)-xb(3))
    nnb(3) = (yb(1)-xb(1))*(zb(2)-xb(2))-(yb(2)-xb(2))*(zb(1)-xb(1))

    CALL getscal(nnb,nna,cosine)

    nn2a = nna(1)*nna(1)+ nna(2)*nna(2)+ nna(3)*nna(3)
    nn2a = SQRT(nn2a)
    IF (nn2a.LT.0.00001_real_8) nn2a=0.00001_real_8
    nn2b = nnb(1)*nnb(1)+ nnb(2)*nnb(2)+ nnb(3)*nnb(3)
    nn2b = SQRT(nn2b)
    IF (nn2b.LT.0.00001_real_8) nn2b=0.00001_real_8
    den  = 1.0_real_8/(nn2a*nn2b)
    cosine = cosine * den

    sinus = SQRT(1.0_real_8-cosine*cosine)
    IF (sinus .LT. 0.00001_real_8) THEN
       sinus=0.00001_real_8
    ENDIF
    teta = dacos(cosine)

    diff=MOD(20000._real_8*pi4+teta0,2.0_real_8*pi4)-teta

    const = -sign0*den/sinus

    der1 = 0.0_real_8
    der2 = -(ya(3)-xa(3))+(za(3)-xa(3))
    der3 = -(za(2)-xa(2))+(ya(2)-xa(2))
    d(1)  = const * (der1*nnb(1)+der2*nnb(2)+der3*nnb(3) -&
         cosine*(der1*nna(1)+der2*nna(2)+der3*nna(3))&
         *nn2b/nn2a )

    der1 = -(za(3)-xa(3))+(ya(3)-xa(3))
    der2 = 0.0_real_8
    der3 = -(ya(1)-xa(1))+(za(1)-xa(1))
    d(2)  = const * (der1*nnb(1)+der2*nnb(2)+der3*nnb(3) -&
         cosine*(der1*nna(1)+der2*nna(2)+der3*nna(3))&
         *nn2b/nn2a )

    der1 = -(ya(2)-xa(2))+(za(2)-xa(2))
    der2 = -(za(1)-xa(1))+(ya(1)-xa(1))
    der3 = 0.0_real_8
    d(3)  = const * (der1*nnb(1)+der2*nnb(2)+der3*nnb(3) -&
         cosine*(der1*nna(1)+der2*nna(2)+der3*nna(3))&
         *nn2b/nn2a )

    der1 = 0.0_real_8
    der2 = -(za(3)-xa(3))
    der3 =  (za(2)-xa(2))
    d(4)  = const * (der1*nnb(1)+der2*nnb(2)+der3*nnb(3) -&
         cosine*(der1*nna(1)+der2*nna(2)+der3*nna(3))&
         *nn2b/nn2a )

    der1 =  (za(3)-xa(3))
    der2 = 0.0_real_8
    der3 = -(za(1)-xa(1))
    d(5)  = const * (der1*nnb(1)+der2*nnb(2)+der3*nnb(3) -&
         cosine*(der1*nna(1)+der2*nna(2)+der3*nna(3))&
         *nn2b/nn2a )

    der1 = -(za(2)-xa(2))
    der2 = +(za(1)-xa(1))
    der3 = 0.0_real_8
    d(6)  = const * (der1*nnb(1)+der2*nnb(2)+der3*nnb(3) -&
         cosine*(der1*nna(1)+der2*nna(2)+der3*nna(3))&
         *nn2b/nn2a )

    der1 = 0.0_real_8
    der2 =  (ya(3)-xa(3))
    der3 = -(ya(2)-xa(2))
    d(7)  = const * (der1*nnb(1)+der2*nnb(2)+der3*nnb(3) -&
         cosine*(der1*nna(1)+der2*nna(2)+der3*nna(3))&
         *nn2b/nn2a )

    der1 = -(ya(3)-xa(3))
    der2 = 0.0_real_8
    der3 =  (ya(1)-xa(1))
    d(8)  = const * (der1*nnb(1)+der2*nnb(2)+der3*nnb(3) -&
         cosine*(der1*nna(1)+der2*nna(2)+der3*nna(3))&
         *nn2b/nn2a )

    der1 =  (ya(2)-xa(2))
    der2 = -(ya(1)-xa(1))
    der3 = 0.0_real_8
    d(9)  = const * (der1*nnb(1)+der2*nnb(2)+der3*nnb(3) -&
         cosine*(der1*nna(1)+der2*nna(2)+der3*nna(3))&
         *nn2b/nn2a )


    der1 = 0.0_real_8
    der2 = -(yb(3)-xb(3))+(zb(3)-xb(3))
    der3 = -(zb(2)-xb(2))+(yb(2)-xb(2))
    d(10)  = const * (der1*nna(1)+der2*nna(2)+der3*nna(3) -&
         cosine*(der1*nnb(1)+der2*nnb(2)+der3*nnb(3))&
         *nn2b/nn2a )

    der1 = -(zb(3)-xb(3))+(yb(3)-xb(3))
    der2 = 0.0_real_8
    der3 = -(yb(1)-xb(1))+(zb(1)-xb(1))
    d(11)  = const * (der1*nna(1)+der2*nna(2)+der3*nna(3) -&
         cosine*(der1*nnb(1)+der2*nnb(2)+der3*nnb(3))&
         *nn2b/nn2a )

    der1 = -(yb(2)-xb(2))+(zb(2)-xb(2))
    der2 = -(zb(1)-xb(1))+(yb(1)-xb(1))
    der3 = 0.0_real_8
    d(12)  = const * (der1*nna(1)+der2*nna(2)+der3*nna(3) -&
         cosine*(der1*nnb(1)+der2*nnb(2)+der3*nnb(3))&
         *nn2b/nn2a )

    der1 = 0.0_real_8
    der2 = -(zb(3)-xb(3))
    der3 =  (zb(2)-xb(2))
    d(13)  = const * (der1*nna(1)+der2*nna(2)+der3*nna(3) -&
         cosine*(der1*nnb(1)+der2*nnb(2)+der3*nnb(3))&
         *nn2b/nn2a )

    der1 =  (zb(3)-xb(3))
    der2 = 0.0_real_8
    der3 = -(zb(1)-xb(1))
    d(14)  = const * (der1*nna(1)+der2*nna(2)+der3*nna(3) -&
         cosine*(der1*nnb(1)+der2*nnb(2)+der3*nnb(3))&
         *nn2b/nn2a )

    der1 = -(zb(2)-xb(2))
    der2 = +(zb(1)-xb(1))
    der3 = 0.0_real_8
    d(15)  = const * (der1*nna(1)+der2*nna(2)+der3*nna(3) -&
         cosine*(der1*nnb(1)+der2*nnb(2)+der3*nnb(3))&
         *nn2b/nn2a )

    der1 = 0.0_real_8
    der2 =  (yb(3)-xb(3))
    der3 = -(yb(2)-xb(2))
    d(16)  = const * (der1*nna(1)+der2*nna(2)+der3*nna(3) -&
         cosine*(der1*nnb(1)+der2*nnb(2)+der3*nnb(3))&
         *nn2b/nn2a )

    der1 = -(yb(3)-xb(3))
    der2 = 0.0_real_8
    der3 =  (yb(1)-xb(1))
    d(17)  = const * (der1*nna(1)+der2*nna(2)+der3*nna(3) -&
         cosine*(der1*nnb(1)+der2*nnb(2)+der3*nnb(3))&
         *nn2b/nn2a )

    der1 =  (yb(2)-xb(2))
    der2 = -(yb(1)-xb(1))
    der3 = 0.0_real_8
    d(18)  = const * (der1*nna(1)+der2*nna(2)+der3*nna(3) -&
         cosine*(der1*nnb(1)+der2*nnb(2)+der3*nnb(3))&
         *nn2b/nn2a )
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE aplane
  ! ==================================================================
  SUBROUTINE rmsd_ab(numspec,ispec,i_rmsd,tscr,natom,&
       an,lsk,fci,fvi,cval)
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: numspec, ispec(numspec), &
                                                i_rmsd
    REAL(real_8)                             :: tscr(3,*)
    INTEGER                                  :: natom
    REAL(real_8)                             :: an(*)
    INTEGER                                  :: lsk(3,ions1%nat)
    REAL(real_8)                             :: fci, fvi, cval

    CHARACTER(*), PARAMETER                  :: procedureN = 'rmsd_ab'
    INTEGER, PARAMETER                       :: iunit = 58 

    CHARACTER(len=100)                       :: filen
    CHARACTER(len=80)                        :: element, label
    INTEGER                                  :: i, ia, iat, icount, icv, &
                                                iend, ierr, iini, is, isp, k, &
                                                l1, l2, l3, nat_rmsd
    INTEGER, ALLOCATABLE, SAVE               :: ipos_rmsd(:)
    INTEGER, SAVE                            :: ifirst = 0
    LOGICAL                                  :: ferror
    REAL(real_8)                             :: rmsd_a, rmsd_b, rmsd_den, &
                                                rmsd_num, urmsd_den
    REAL(real_8), ALLOCATABLE                :: derr_a(:,:), derr_b(:,:), &
                                                r0_cm(:,:), r_cm(:,:), &
                                                r_ion(:,:), ra_ion(:,:), &
                                                rb_ion(:,:)
    REAL(real_8), ALLOCATABLE, SAVE          :: pos_a(:,:,:), pos_b(:,:,:), &
                                                rmsd_scr(:)

! ==--------------------------------------------------------------==

    IF (ifirst.EQ.0) THEN
       ! Allocate memory
       ALLOCATE(pos_a(3,natom,nrmsd_ab),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(pos_b(3,natom,nrmsd_ab),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(ipos_rmsd(natom),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(rmsd_scr(7*3*natom),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)

       icount = 0
       DO icv = 1,ncolvar
          IF (icv_rmsd_ab(icv) .NE. 0) THEN
             icount = icount + 1
             ! Read the reference configurations
             filen=file_str_ab(icv)
             ferror=.FALSE.
             IF (paral%io_parent)&
                  CALL fileopen(iunit,filen,fo_old,ferror)
             IF (ferror) THEN
                IF (paral%io_parent)&
                     WRITE(6,'(A,A,A)')&
                     ' RMSD_AB! File',filen,'  DOES NOT EXIST!!'
                CALL stopgm('RMSD_AB','CANNOT READ A and B COORDINATES',& 
                     __LINE__,__FILE__)
             ENDIF
             IF (paral%io_parent)&
                  REWIND(iunit)
             IF (paral%io_parent)&
                  READ(iunit,err=998,fmt=*) label
             DO iat=1,ions1%nat
                IF (paral%io_parent)&
                     READ(iunit,err=998,fmt=*) element,&
                     (pos_a(k,iat,icount),k=1,3)
             ENDDO
             IF (paral%io_parent)&
                  WRITE(6,'(/A,A,A/)') ' A COORDINATES READ ',&
                  'FROM FILE ',filen
             IF (paral%io_parent)&
                  READ(iunit,err=998,fmt=*) label
             DO iat=1,ions1%nat
                IF (paral%io_parent)&
                     READ(iunit,err=998,fmt=*) element,&
                     (pos_b(k,iat,icount),k=1,3)
             ENDDO
             IF (paral%io_parent)&
                  WRITE(6,'(/A,A,A/)') ' B COORDINATES READ ',&
                  'FROM FILE ',filen
             IF (paral%io_parent)&
                  CALL fileclose(iunit)
             GOTO 999

998          CONTINUE
             IF (paral%io_parent)&
                  WRITE(6,'(A)') ' RMSD_AB! ERROR WHILE READING COORDINATES'
             CALL stopgm('RMSD_AB','READ ERROR',& 
                  __LINE__,__FILE__)

999          CONTINUE
          ENDIF! ICV_RMSD_AB
       ENDDO ! NCOLVAR
       ifirst=1

       IF (icount .NE. nrmsd_ab) THEN
          CALL stopgm('RMSD_AB','ERROR in NRMSD_AB',& 
               __LINE__,__FILE__)
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    ! write(6,*) 'pos rmsd ', I_RMSD
    ! do iat = 1,natom
    ! write(6,*) POS_A(1,IAT,I_RMSD),POS_B(1,IAT,I_RMSD)
    ! enddo
    ! STOP 'conf read'
    CALL zeroing(ipos_rmsd)!,natom)
    CALL zeroing(rmsd_scr)!,7*3*natom)

    ! Create the list of engaged ions

    nat_rmsd = 0
    IF (numspec .EQ. 0) THEN
       iat = 0
       DO is = 1,ions1%nsp
          DO ia = 1,ions0%na(is)
             nat_rmsd = nat_rmsd + 1
             iat = iat + ia
             ipos_rmsd(nat_rmsd) = iat
          ENDDO
       ENDDO
    ELSE
       DO is=1,numspec
          isp = ispec(is)
          iini = 0
          DO i = 1,isp-1
             iini = iini + ions0%na(i)
          ENDDO
          iend = iini + ions0%na(isp)
          iini = iini + 1
          DO iat = iini,iend
             nat_rmsd = nat_rmsd + 1
             ipos_rmsd(nat_rmsd) = iat
          ENDDO
       ENDDO
    ENDIF   ! NUMSPEC

    ! TODO align for BG
    ALLOCATE(r_ion(3, nat_rmsd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(ra_ion(3, nat_rmsd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(rb_ion(3, nat_rmsd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    ALLOCATE(derr_a(3, nat_rmsd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(derr_b(3, nat_rmsd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    ALLOCATE(r_cm(3, nat_rmsd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(r0_cm(3, nat_rmsd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    !$omp parallel do private(IAT)
    DO iat = 1,nat_rmsd
       r_ion(1,iat)=tscr(1,ipos_rmsd(iat))
       r_ion(2,iat)=tscr(2,ipos_rmsd(iat))
       r_ion(3,iat)=tscr(3,ipos_rmsd(iat))
       ra_ion(1,iat)=pos_a(1,ipos_rmsd(iat),i_rmsd)
       ra_ion(2,iat)=pos_a(2,ipos_rmsd(iat),i_rmsd)
       ra_ion(3,iat)=pos_a(3,ipos_rmsd(iat),i_rmsd)
       rb_ion(1,iat)=pos_b(1,ipos_rmsd(iat),i_rmsd)
       rb_ion(2,iat)=pos_b(2,ipos_rmsd(iat),i_rmsd)
       rb_ion(3,iat)=pos_b(3,ipos_rmsd(iat),i_rmsd)
    ENDDO

    ! write(6,*) 'NAT_RMSD', NAT_RMSD
    ! ==--------------------------------------------------------------==
    ! Calculate RMSD_B and derivatives
    CALL val_rmsd(nat_rmsd,r_ion,rb_ion,rmsd_b,derr_b,r0_cm,r_cm)
    ! write(6,*) 'B RMSD', RMSD_B

    ! Calculate RMSD_A and derivatives
    CALL val_rmsd(nat_rmsd,r_ion,ra_ion,rmsd_a,derr_a,r0_cm,r_cm)
    ! write(6,*) 'A RMSD', RMSD_A

    ! ==--------------------------------------------------------------==
    rmsd_num = rmsd_a-rmsd_b
    rmsd_den = rmsd_a+rmsd_b
    IF (rmsd_den .LT. 1.e-10_real_8) rmsd_den = 1.e-10_real_8
    urmsd_den = 1.0_real_8/rmsd_den
    fvi = rmsd_num*urmsd_den
    fci=fvi-cval

    DO iat = 1,nat_rmsd
       l1 = lsk(1,ipos_rmsd(iat))
       l2 = lsk(2,ipos_rmsd(iat))
       l3 = lsk(3,ipos_rmsd(iat))

       IF (l1.NE.0) an(l1) = an(l1) + urmsd_den *(&
            derr_a(1,ipos_rmsd(iat))-derr_b(1,ipos_rmsd(iat)) -&
            (derr_a(1,ipos_rmsd(iat))+derr_b(1,ipos_rmsd(iat)))*fvi)
       IF (l2.NE.0) an(l2) = an(l2) + urmsd_den *(&
            derr_a(2,ipos_rmsd(iat))-derr_b(2,ipos_rmsd(iat)) -&
            (derr_a(2,ipos_rmsd(iat))+derr_b(2,ipos_rmsd(iat)))*fvi)
       IF (l3.NE.0) an(l3) = an(l3) + urmsd_den *(&
            derr_a(3,ipos_rmsd(iat))-derr_b(3,ipos_rmsd(iat)) -&
            (derr_a(3,ipos_rmsd(iat))+derr_b(3,ipos_rmsd(iat)))*fvi)

    ENDDO
    ! ==--------------------------------------------------------------==
    DEALLOCATE(r_ion,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(ra_ion,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(rb_ion,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)

    DEALLOCATE(derr_a,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(derr_b,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)

    DEALLOCATE(r_cm,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(r0_cm,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rmsd_ab
  ! ==================================================================
  SUBROUTINE val_rmsd(n,r,r0,rmsd,derr,rp,r0p)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: n
    REAL(real_8)                             :: r(3,n), r0(3,n), rmsd, &
                                                derr(3,n), rp(3,n), r0p(3,n)

    REAL(real_8), PARAMETER                  :: epsi = 1.0e-10_real_8

    CHARACTER(len=80)                        :: msg
    INTEGER                                  :: i, ier, iret, ix, j, k
    REAL(real_8)                             :: dm_r(4,4,3), lambda(4), &
                                                m(4,4), q(0:3), rr(3), &
                                                rr0(3), rrsq, s, wk(20), xx, &
                                                yy, z(4,4), zz

! ==--------------------------------------------------------------==
! center r

    xx=0._real_8
    yy=0._real_8
    zz=0._real_8
    !$omp parallel do private(i) reduction(+:xx,yy,zz)
    DO i=1,n
       xx=xx+r(1,i)
       yy=yy+r(2,i)
       zz=zz+r(3,i)
    ENDDO
    xx=xx/REAL(n,kind=real_8)
    yy=yy/REAL(n,kind=real_8)
    zz=zz/REAL(n,kind=real_8)
    !$omp parallel do private(i) shared(xx,yy,zz)
    DO i=1,n
       rp(1,i)=r(1,i)-xx
       rp(2,i)=r(2,i)-yy
       rp(3,i)=r(3,i)-zz
    ENDDO
    ! center r0
    xx=0._real_8
    yy=0._real_8
    zz=0._real_8
    !$omp parallel do private(i) reduction(+:xx,yy,zz)
    DO i=1,n
       xx=xx+r0(1,i)
       yy=yy+r0(2,i)
       zz=zz+r0(3,i)
    ENDDO
    xx=xx/REAL(n,kind=real_8)
    yy=yy/REAL(n,kind=real_8)
    zz=zz/REAL(n,kind=real_8)
    !$omp parallel do private(i) shared(xx,yy,zz)
    DO i=1,n
       r0p(1,i)=r0(1,i)-xx
       r0p(2,i)=r0(2,i)-yy
       r0p(3,i)=r0(3,i)-zz
    ENDDO
    ! 
    DO i = 1,16
       m(i,1) = 0._real_8
    ENDDO

    DO i = 1,n
       rr(1)=rp(1,i)
       rr(2)=rp(2,i)
       rr(3)=rp(3,i)
       rr0(1)=r0p(1,i)
       rr0(2)=r0p(2,i)
       rr0(3)=r0p(3,i)
       rrsq=(rr0(1)**2+rr0(2)**2+rr0(3)**2+rr(1)**2+rr(2)**2+rr(3)**2)

       m(1,1)=m(1,1) +&
            rrsq+2._real_8*(-rr0(1)*rr(1)-rr0(2)*rr(2)-rr0(3)*rr(3))
       m(2,2)=m(2,2) +&
            rrsq+2._real_8*(-rr0(1)*rr(1)+rr0(2)*rr(2)+rr0(3)*rr(3))
       m(3,3)=m(3,3) +&
            rrsq+2._real_8*(+rr0(1)*rr(1)-rr0(2)*rr(2)+rr0(3)*rr(3))
       m(4,4)=m(4,4) +&
            rrsq+2._real_8*(+rr0(1)*rr(1)+rr0(2)*rr(2)-rr0(3)*rr(3))
       m(1,2)=m(1,2) +2._real_8*(-rr0(2)*rr(3)+rr0(3)*rr(2))
       m(1,3)=m(1,3) +2._real_8*(rr0(1)*rr(3)-rr0(3)*rr(1))
       m(1,4)=m(1,4) +2._real_8*(-rr0(1)*rr(2)+rr0(2)*rr(1))
       m(2,3)=m(2,3) -2._real_8*(rr0(1)*rr(2)+rr0(2)*rr(1))
       m(2,4)=m(2,4) -2._real_8*(rr0(1)*rr(3)+rr0(3)*rr(1))
       m(3,4)=m(3,4) -2._real_8*(rr0(2)*rr(3)+rr0(3)*rr(2))
    ENDDO
    m(2,1) = m(1,2)
    m(3,1) = m(1,3)
    m(3,2) = m(2,3)
    m(4,1) = m(1,4)
    m(4,2) = m(2,4)
    m(4,3) = m(3,4)

    ! ---- SOLVE THE EIGENVECTOR PROBLEM FOR M: --------------------------*

    iret = 0
    ! jobn   = 12
    ! call EIGRS(m,4,jobn,lambda,z,4,wk,ier)
    ier=0
    CALL dcopy(16,m,1,z,1)
    CALL dsyev('V','U',4,z,4,lambda,wk,20,ier)
    IF (ier .NE. 0) THEN
       CALL stopgm('VAL_RMSD',' Fatal error in DSYEV',& 
            __LINE__,__FILE__)
    ENDIF
    IF (wk(1) .GT. 1._real_8) THEN
       iret = 9
    ENDIF

    ! ---- PICK THE CORRECT EIGENVECTOR(S): ------------------------------*

    s =  1._real_8
    IF (z(1,1) .LT. 0._real_8) s = -1._real_8
    q(0) = s*z(1,1)
    q(1) = s*z(2,1)
    q(2) = s*z(3,1)
    q(3) = s*z(4,1)
    IF (ABS(lambda(1)) .LT. epsi) THEN
       rmsd = 0._real_8
    ELSE
       rmsd = lambda(1)/REAL(n,kind=real_8)
    ENDIF

    IF (ABS(lambda(1) - lambda(2)) .LT. epsi)iret = iret + 10

    IF (iret.EQ.0) THEN
       msg = ' RMSD: Normal exe., unique solution'
    ENDIF
    IF (iret.EQ.10) THEN
       msg = ' RMSD: Normal exe., non-unique solution'
    ENDIF
    IF (iret.EQ.9) THEN
       msg = ' RMSD: Bad perform., unique solution'
    ENDIF
    IF (iret.EQ.19) THEN
       msg = ' RMSD: Bad perform., non-unique solution'
    ENDIF
    IF (iret.NE.0) THEN
       ! write(6,*)msg
    ENDIF

    ! ---- DERIVATIVES OF RMSD with respect to the positions
    DO i = 1,n

       rr(1)=  2._real_8*rp(1,i)
       rr(2)=  2._real_8*rp(2,i)
       rr(3)=  2._real_8*rp(3,i)
       rr0(1)= 2._real_8*r0p(1,i)
       rr0(2)= 2._real_8*r0p(2,i)
       rr0(3)= 2._real_8*r0p(3,i)
       ! 
       dm_r (1,1,1)=(rr(1)-rr0(1))
       dm_r (1,1,2)=(rr(2)-rr0(2))
       dm_r (1,1,3)=(rr(3)-rr0(3))
       ! 
       dm_r (1,2,1)=0._real_8
       dm_r (1,2,2)= rr0(3)
       dm_r (1,2,3)=-rr0(2)
       ! 
       dm_r (1,3,1)=-rr0(3)
       dm_r (1,3,2)= 0._real_8
       dm_r (1,3,3)= rr0(1)
       ! 
       dm_r (1,4,1)= rr0(2)
       dm_r (1,4,2)=-rr0(1)
       dm_r (1,4,3)= 0._real_8
       ! 
       dm_r (2,2,1)=(rr(1)-rr0(1))
       dm_r (2,2,2)=(rr(2)+rr0(2))
       dm_r (2,2,3)=(rr(3)+rr0(3))
       ! 
       dm_r (2,3,1)=-rr0(2)
       dm_r (2,3,2)=-rr0(1)
       dm_r (2,3,3)= 0._real_8
       ! 
       dm_r (2,4,1)=-rr0(3)
       dm_r (2,4,2)= 0._real_8
       dm_r (2,4,3)=-rr0(1)
       ! 
       dm_r (3,3,1)=(rr(1)+rr0(1))
       dm_r (3,3,2)=(rr(2)-rr0(2))
       dm_r (3,3,3)=(rr(3)+rr0(3))
       ! 
       dm_r (3,4,1)=0._real_8
       dm_r (3,4,2)=-rr0(3)
       dm_r (3,4,3)=-rr0(2)
       ! 
       dm_r (4,4,1)=(rr(1)+rr0(1))
       dm_r (4,4,2)=(rr(2)+rr0(2))
       dm_r (4,4,3)=(rr(3)-rr0(3))
       ! 
       DO ix=1,3
          dm_r(2,1,ix)=dm_r(1,2,ix)
          dm_r(3,1,ix)=dm_r(1,3,ix)
          dm_r(4,1,ix)=dm_r(1,4,ix)
          dm_r(3,2,ix)=dm_r(2,3,ix)
          dm_r(4,2,ix)=dm_r(2,4,ix)
          dm_r(4,3,ix)=dm_r(3,4,ix)
       ENDDO
       ! 
       DO ix=1,3
          derr(ix,i)=0._real_8
          DO k=1,4
             DO j=1,4
                derr(ix,i)=derr(ix,i)+q(k-1)*q(j-1)*dm_r(j,k,ix)
             ENDDO
          ENDDO
          derr(ix,i)=derr(ix,i)/REAL(n,kind=real_8)
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE val_rmsd
  ! ==================================================================
  SUBROUTINE hyd_presence(isp1,isp2,nexp,mexp,c1_rc,l,tscr,an,lsk,&
       cval)
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: isp1, isp2, nexp(3), mexp(3)
    REAL(real_8)                             :: c1_rc(3), l, tscr(3,*), &
                                                an(cotc0%nodim)
    INTEGER                                  :: lsk(3,ions1%nat)
    REAL(real_8)                             :: cval

    CHARACTER(*), PARAMETER                  :: procedureN = 'hyd_presence'

    INTEGER                                  :: ia, iat, iend1, iend2, ierr, &
                                                ii, iiat1, iiat2, iini1, &
                                                iini2, is, k, l1, l2, l3
    REAL(real_8) :: dd, der_f_nH_dn, der_f_nO_dn, der_fact, dx, dy, dz, &
      exp_fact, f_nH, f_nO, fden, ff, fnum, nH, nO, rden, rnum, x0, xx(3), &
      y0, z0
    REAL(real_8), ALLOCATABLE                :: der_nH_dr(:,:), der_nO_dr(:,:)

    cval = 0.0_real_8

    ALLOCATE(der_nH_dr(3,ions0%na(isp2)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(der_nO_dr(3,ions0%na(isp1)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(an)!,cotc0%nodim)
    ! Atoms indexes
    iini1 = 1
    DO is = 1,isp1-1
       iini1=iini1+ions0%na(is)
    ENDDO
    iend1 = iini1 + ions0%na(isp1) -1
    iini2 = 1
    DO is = 1,isp2-1
       iini2=iini2+ions0%na(is)
    ENDDO
    iend2 = iini2 + ions0%na(isp2) -1

    ! First Loop on the oxygens
    DO ia=iini1,iend1
       ! DO IA=94,94
       nH=0.0_real_8
       nO=0.0_real_8
       f_nH=0.0_real_8
       CALL zeroing(der_nH_dr)!,3*ions0%na(isp2))
       CALL zeroing(der_nO_dr)!,3*ions0%na(isp1))
       der_f_nH_dn =0.0_real_8
       IF (ia.LE.ions1%nat) THEN
          l1 = lsk(1,ia)
          l2 = lsk(2,ia)
          l3 = lsk(3,ia)
       ENDIF
       CALL fillc(ia,tscr,xx)
       x0 = xx(1)
       y0 = xx(2)
       z0 = xx(3)

       ! Loop over the hydrogens
       ii = 0
       DO iiat2=iini2,iend2
          ii = ii + 1
          dx=tscr(1,iiat2)-x0
          dy=tscr(2,iiat2)-y0
          dz=tscr(3,iiat2)-z0
          IF (.NOT. isos1%tisos) CALL pbc(dx,dy,dz,dx,dy,dz,1,parm%apbc,parm%ibrav)
          dd=SQRT(dx*dx+dy*dy+dz*dz)
          IF (dd .LT. 1.e-3_real_8)&
               CALL stopgm('HYD_PRESENCE','H/O position overlap',& 
               __LINE__,__FILE__)
          ! compute nH
          rnum = (dd/c1_rc(1))**nexp(1)
          rden = (dd/c1_rc(1))**(nexp(1)+mexp(1))
          fnum = 1-rnum
          fden = 1-rden
          IF (ABS(fden) .LT. 1.e-10_real_8) fden = 1.e-10_real_8
          nH= nh + fnum/fden
          ! write(6,*) 'H ', IIAT2, DD, FNUM/FDEN
          ! compute der_nH_dr  index dependent
          ff = -(rnum*REAL(nexp(1),kind=real_8)&
               -fnum*rden*REAL(nexp(1)+mexp(1),kind=real_8)/fden)/&
               (fden*dd*dd)
          der_nH_dr(1,ii) = ff*dx
          der_nH_dr(2,ii) = ff*dy
          der_nH_dr(3,ii) = ff*dz


       ENDDO   ! Loop on the hydrogens


       ! compute f_nH
       rnum = (nH/c1_rc(3))**nexp(3)
       rden = (nH/c1_rc(3))**(nexp(3)+mexp(3))
       fnum = 1-rnum
       fden = 1-rden
       IF (ABS(fden) .LT. 1.e-10_real_8) fden = 1.e-10_real_8
       f_nH= 1.0_real_8 - fnum/fden

       ! compute der_f_nH_dn
       der_f_nH_dn = (REAL(nexp(3),kind=real_8)*rnum-&
            REAL(nexp(3)+mexp(3),kind=real_8)*fnum*rden/fden)/(fden*nH)

       ii = 0
       ! Second Loop over the oxygens
       DO iiat1 = iini1,iend1
          ii = ii + 1
          dx=tscr(1,iiat1)-x0
          dy=tscr(2,iiat1)-y0
          dz=tscr(3,iiat1)-z0
          IF (.NOT. isos1%tisos) CALL pbc(dx,dy,dz,dx,dy,dz,1,parm%apbc,parm%ibrav)
          dd=SQRT(dx*dx+dy*dy+dz*dz)
          IF (dd .LT. 1.e-3_real_8) GOTO 10
          ! compute nO
          rnum = (dd/c1_rc(2))**nexp(2)
          rden = (dd/c1_rc(2))**(nexp(2)+mexp(2))
          fnum = 1-rnum
          fden = 1-rden
          IF (ABS(fden) .LT. 1.e-10_real_8) fden = 1.e-10_real_8
          nO= no + fnum/fden

          ! compute der_nO_dr  
          ff = -(rnum*REAL(nexp(2),kind=real_8)&
               -fnum*rden*REAL(nexp(2)+mexp(2),kind=real_8)/fden)/&
               (fden*dd*dd)
          der_nO_dr(1,ii) = ff*dx
          der_nO_dr(2,ii) = ff*dy
          der_nO_dr(3,ii) = ff*dz
10        CONTINUE
       ENDDO ! Second Loop over the oxygens

       ! compute f_nO
       rnum = (nO/c1_rc(3))**nexp(3)
       rden = (nO/c1_rc(3))**(nexp(3)+mexp(3))
       fnum = 1-rnum
       fden = 1-rden
       IF (ABS(fden) .LT. 1.e-10_real_8) fden = 1.e-10_real_8
       f_nO= 1.0_real_8 - fnum/fden


       ! compute der_f_nO_dn
       der_f_nO_dn = (REAL(nexp(3),kind=real_8)*rnum-&
            REAL(nexp(3)+mexp(3),kind=real_8)*fnum*rden/fden)/(fden*nO)

       ! exponential factor exp{L*f_nH*f_nO}
       exp_fact = EXP(l* f_nH * f_nO)
       der_fact = l*exp_fact*f_nH*der_f_nO_dn


       ii = 0
       DO iiat1 = iini1,iend1
          ii = ii + 1
          ! Derivatives:
          IF (lsk(1,iiat1).NE.0) THEN
             k=lsk(1,iiat1)
             an(k)=an(k)+ der_fact*der_nO_dr(1,ii)
          ENDIF
          IF (lsk(2,iiat1).NE.0) THEN
             k=lsk(2,iiat1)
             an(k)=an(k)+ der_fact*der_nO_dr(2,ii)
          ENDIF
          IF (lsk(3,iiat1).NE.0) THEN
             k=lsk(3,iiat1)
             an(k)=an(k)+ der_fact*der_nO_dr(3,ii)
          ENDIF
          IF (l1.NE.0) an(l1)=an(l1)-der_fact*der_nO_dr(1,ii)
          IF (l2.NE.0) an(l2)=an(l2)-der_fact*der_nO_dr(2,ii)
          IF (l3.NE.0) an(l3)=an(l3)-der_fact*der_nO_dr(3,ii)

       ENDDO


       der_fact = l*exp_fact*der_f_nH_dn*f_nO
       ii = 0
       ! Loop over the hydrogens for the derivatives
       DO iiat2=iini2,iend2
          ii = ii + 1
          IF (lsk(1,iiat2).NE.0) THEN
             k=lsk(1,iiat2)
             an(k)=an(k)+ der_fact*der_nH_dr(1,ii)
          ENDIF
          IF (lsk(2,iiat2).NE.0) THEN
             k=lsk(2,iiat2)
             an(k)=an(k)+ der_fact*der_nH_dr(2,ii)
          ENDIF
          IF (lsk(3,iiat2).NE.0) THEN
             k=lsk(3,iiat2)
             an(k)=an(k)+ der_fact*der_nH_dr(3,ii)
          ENDIF
          IF (l1.NE.0) an(l1)=an(l1)-der_fact*der_nH_dr(1,ii)
          IF (l2.NE.0) an(l2)=an(l2)-der_fact*der_nH_dr(2,ii)
          IF (l3.NE.0) an(l3)=an(l3)-der_fact*der_nH_dr(3,ii)
       ENDDO
       ! sum over the oxygens
       cval = cval + exp_fact
    ENDDO    ! Loop over the oxygens

    ff = 1.0_real_8/(l*cval)

    cval = 1.0_real_8/l*LOG(cval)


    ! Final factor for the derivatives
    DO iat = 1,cotc0%nodim
       an(iat) = an(iat)*ff
    ENDDO

    DEALLOCATE(der_nH_dr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(der_nO_dr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE hyd_presence
  ! ==--------------------------------------------------------------==
  ! ==--------------------------------------------------------------==
  SUBROUTINE hyd2_presence(isp1,isp2,nexp,mexp,c1_rc,l,tscr,an,lsk,&
       cval)
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: isp1, isp2, nexp(3), mexp(3)
    REAL(real_8)                             :: c1_rc(3), l, tscr(3,*), &
                                                an(cotc0%nodim)
    INTEGER                                  :: lsk(3,ions1%nat)
    REAL(real_8)                             :: cval

    CHARACTER(*), PARAMETER                  :: procedureN = 'hyd2_presence'

    INTEGER                                  :: ia, iat, iend1, iend2, ierr, &
                                                ii, iiat1, iiat2, iini1, &
                                                iini2, is, k, l1, l2, l3
    REAL(real_8) :: dd, der_f_nH_dn, der_f_nO_dn, der_fact, dx, dy, dz, &
      exp_fact, f_nH, f_nO, fden, ff, fnum, nH, nO, rden, rnum, x0, xx(3), &
      y0, z0
    REAL(real_8), ALLOCATABLE                :: der_nH_dr(:,:), der_nO_dr(:,:)

    cval = 0.0_real_8

    ALLOCATE(der_nH_dr(3,ions0%na(isp2)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(der_nO_dr(3,ions0%na(isp1)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(an)!,cotc0%nodim)
    ! Atoms indexes
    iini1 = 1
    DO is = 1,isp1-1
       iini1=iini1+ions0%na(is)
    ENDDO
    iend1 = iini1 + ions0%na(isp1) -1
    iini2 = 1
    DO is = 1,isp2-1
       iini2=iini2+ions0%na(is)
    ENDDO
    iend2 = iini2 + ions0%na(isp2) -1

    ! First Loop on the oxygens
    DO ia=iini1,iend1
       nH=0.0_real_8
       nO=0.0_real_8
       f_nH=0.0_real_8
       CALL zeroing(der_nH_dr)!,3*ions0%na(isp2))
       CALL zeroing(der_nO_dr)!,3*ions0%na(isp1))
       der_f_nH_dn =0.0_real_8
       IF (ia.LE.ions1%nat) THEN
          l1 = lsk(1,ia)
          l2 = lsk(2,ia)
          l3 = lsk(3,ia)
       ENDIF
       CALL fillc(ia,tscr,xx)
       x0 = xx(1)
       y0 = xx(2)
       z0 = xx(3)

       ! Loop over the hydrogens
       ii = 0
       DO iiat2=iini2,iend2
          ii = ii + 1
          dx=tscr(1,iiat2)-x0
          dy=tscr(2,iiat2)-y0
          dz=tscr(3,iiat2)-z0
          IF (.NOT. isos1%tisos) CALL pbc(dx,dy,dz,dx,dy,dz,1,parm%apbc,parm%ibrav)
          dd=SQRT(dx*dx+dy*dy+dz*dz)
          IF (dd .LT. 1.e-3_real_8)&
               CALL stopgm('HYD2_PRESENCE','H/O position overlap',& 
               __LINE__,__FILE__)
          ! compute nH
          rnum = (dd/c1_rc(1))**nexp(1)
          rden = (dd/c1_rc(1))**(nexp(1)+mexp(1))
          fnum = 1-rnum
          fden = 1-rden
          IF (ABS(fden) .LT. 1.e-10_real_8) fden = 1.e-10_real_8
          nH= nh + fnum/fden
          ! write(6,*) 'H ', IIAT2, DD, FNUM/FDEN
          ! compute der_nH_dr  index dependent
          ff = -(rnum*REAL(nexp(1),kind=real_8)&
               -fnum*rden*REAL(nexp(1)+mexp(1),kind=real_8)/fden)/&
               (fden*dd*dd)
          der_nH_dr(1,ii) = ff*dx
          der_nH_dr(2,ii) = ff*dy
          der_nH_dr(3,ii) = ff*dz


       ENDDO   ! Loop on the hydrogens


       ! compute f_nH
       rnum = (nH/c1_rc(3))**nexp(3)
       rden = (nH/c1_rc(3))**(nexp(3)+mexp(3))
       fnum = 1-rnum
       fden = 1-rden
       IF (ABS(fden) .LT. 1.e-10_real_8) fden = 1.e-10_real_8
       f_nH= 1.0_real_8 - fnum/fden

       ! compute der_f_nH_dn
       der_f_nH_dn = (REAL(nexp(3),kind=real_8)*rnum-&
            REAL(nexp(3)+mexp(3),kind=real_8)*fnum*rden/fden)/(fden*nH)

       ii = 0
       ! Second Loop over the oxygens 
       DO iiat1 = iini1,iend1
          ii = ii + 1
          dx=tscr(1,iiat1)-x0
          dy=tscr(2,iiat1)-y0
          dz=tscr(3,iiat1)-z0
          IF (.NOT. isos1%tisos) CALL pbc(dx,dy,dz,dx,dy,dz,1,parm%apbc,parm%ibrav)
          dd=SQRT(dx*dx+dy*dy+dz*dz)
          IF (dd .LT. 1.e-3_real_8) GOTO 10
          ! compute nO
          rnum = (dd/c1_rc(2))**nexp(2)
          rden = (dd/c1_rc(2))**(nexp(2)+mexp(2))
          fnum = 1-rnum
          fden = 1-rden
          IF (ABS(fden) .LT. 1.e-10_real_8) fden = 1.e-10_real_8
          nO= no + fnum/fden

          ! compute der_nO_dr  
          ff = -(rnum*REAL(nexp(2),kind=real_8)&
               -fnum*rden*REAL(nexp(2)+mexp(2),kind=real_8)/fden)/&
               (fden*dd*dd)
          der_nO_dr(1,ii) = ff*dx
          der_nO_dr(2,ii) = ff*dy
          der_nO_dr(3,ii) = ff*dz
10        CONTINUE
       ENDDO ! Second Loop over the oxygens

       ! compute f_nO
       rnum = (nO/c1_rc(3))**nexp(3)
       rden = (nO/c1_rc(3))**(nexp(3)+mexp(3))
       fnum = 1-rnum
       fden = 1-rden
       IF (ABS(fden) .LT. 1.e-10_real_8) fden = 1.e-10_real_8
       f_nO= 1.0_real_8 - fnum/fden


       ! compute der_f_nO_dn
       der_f_nO_dn = (REAL(nexp(3),kind=real_8)*rnum-&
            REAL(nexp(3)+mexp(3),kind=real_8)*fnum*rden/fden)/(fden*nO)

       ! exponential factor exp{L*f_nH * nO}
       exp_fact = EXP(l* f_nH * nO)
       der_fact = l*exp_fact*f_nH

       ii = 0
       DO iiat1 = iini1,iend1
          ii = ii + 1
          ! Derivatives:
          IF (lsk(1,iiat1).NE.0) THEN
             k=lsk(1,iiat1)
             an(k)=an(k)+ der_fact*der_nO_dr(1,ii)
          ENDIF
          IF (lsk(2,iiat1).NE.0) THEN
             k=lsk(2,iiat1)
             an(k)=an(k)+ der_fact*der_nO_dr(2,ii)
          ENDIF
          IF (lsk(3,iiat1).NE.0) THEN
             k=lsk(3,iiat1)
             an(k)=an(k)+ der_fact*der_nO_dr(3,ii)
          ENDIF
          IF (l1.NE.0) an(l1)=an(l1)-der_fact*der_nO_dr(1,ii)
          IF (l2.NE.0) an(l2)=an(l2)-der_fact*der_nO_dr(2,ii)
          IF (l3.NE.0) an(l3)=an(l3)-der_fact*der_nO_dr(3,ii)

       ENDDO


       der_fact = l*exp_fact*der_f_nH_dn*nO
       ii = 0
       ! Loop over the hydrogens for the derivatives
       DO iiat2=iini2,iend2
          ii = ii + 1
          IF (lsk(1,iiat2).NE.0) THEN
             k=lsk(1,iiat2)
             an(k)=an(k)+ der_fact*der_nH_dr(1,ii)
          ENDIF
          IF (lsk(2,iiat2).NE.0) THEN
             k=lsk(2,iiat2)
             an(k)=an(k)+ der_fact*der_nH_dr(2,ii)
          ENDIF
          IF (lsk(3,iiat2).NE.0) THEN
             k=lsk(3,iiat2)
             an(k)=an(k)+ der_fact*der_nH_dr(3,ii)
          ENDIF
          IF (l1.NE.0) an(l1)=an(l1)-der_fact*der_nH_dr(1,ii)
          IF (l2.NE.0) an(l2)=an(l2)-der_fact*der_nH_dr(2,ii)
          IF (l3.NE.0) an(l3)=an(l3)-der_fact*der_nH_dr(3,ii)
       ENDDO
       ! sum over the oxygens
       cval = cval + exp_fact
    ENDDO    ! Loop over the oxygens

    ff = 1.0_real_8/(l*cval)

    cval = 1.0_real_8/l*LOG(cval)


    ! Final factor for the derivatives
    DO iat = 1,cotc0%nodim
       an(iat) = an(iat)*ff
    ENDDO

    DEALLOCATE(der_nH_dr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(der_nO_dr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE hyd2_presence
  ! ==================================================================
  SUBROUTINE hyd2_ion_distance(isp1,isp2,indion,nexp,mexp,c1_rc,l,&
       tscr,an,lsk,cval)
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: isp1, isp2, indion, nexp(3), &
                                                mexp(3)
    REAL(real_8)                             :: c1_rc(3), l, tscr(3,*), &
                                                an(cotc0%nodim)
    INTEGER                                  :: lsk(3,ions1%nat)
    REAL(real_8)                             :: cval

    CHARACTER(*), PARAMETER :: procedureN = 'hyd2_ion_distance'

    INTEGER                                  :: ia, iat, iend1, iend2, ierr, &
                                                ii, iiat2, iini1, iini2, is, &
                                                k, l1, l11, l2, l22, l3, l33
    REAL(real_8) :: dd, der_fact, der_fact2, der_r_k_dr(3), dx, dx2, dy, dy2, &
      dz, dz2, exp_fact, fden, ff, fnum, nH, r_k, rden, rnum, sum_exp_fact, &
      sum_r_exp_fact, x0, xion, xx(3), xx0(3), y0, yion, z0, zion
    REAL(real_8), ALLOCATABLE                :: anf1(:), der_nH_dr(:,:)

! Variables
! jpark  
! jpark
! ==--------------------------------------------------------------==

    cval = 0.0_real_8
    sum_exp_fact = 0.0_real_8
    sum_r_exp_fact = 0.0_real_8

    ALLOCATE(der_nH_dr(3,ions0%na(isp2)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(der_nH_dr)!,3*ions0%na(isp2))

    ALLOCATE(anf1(cotc0%nodim),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(an)!,cotc0%nodim)
    ! jpark 

    CALL zeroing(anf1)!,cotc0%nodim)
    ! Atoms indexes
    iini1 = 1
    DO is = 1,isp1-1
       iini1=iini1+ions0%na(is)
    ENDDO
    iend1 = iini1 + ions0%na(isp1) -1
    iini2 = 1
    DO is = 1,isp2-1
       iini2=iini2+ions0%na(is)
    ENDDO
    iend2 = iini2 + ions0%na(isp2) -1

    ! jpark     SET ION POSITION 
    CALL fillc(indion,tscr,xx0)
    xion = xx0(1)
    yion = xx0(2)
    zion = xx0(3)

    l11 = lsk(1,indion)
    l22 = lsk(2,indion)
    l33 = lsk(3,indion)


    ! First Loop on the oxygens
    DO ia=iini1,iend1
       IF (ia .EQ. indion) GOTO 100
       nH=0.0_real_8
       CALL zeroing(der_nH_dr)!,3*ions0%na(isp2))
       IF (ia.LE.ions1%nat) THEN
          l1 = lsk(1,ia)
          l2 = lsk(2,ia)
          l3 = lsk(3,ia)
       ENDIF
       CALL fillc(ia,tscr,xx)
       x0 = xx(1)
       y0 = xx(2)
       z0 = xx(3)

       ! Loop over the hydrogens
       ii = 0
       DO iiat2=iini2,iend2
          ii = ii + 1
          dx=tscr(1,iiat2)-x0
          dy=tscr(2,iiat2)-y0
          dz=tscr(3,iiat2)-z0
          IF (.NOT. isos1%tisos) CALL pbc(dx,dy,dz,dx,dy,dz,1,parm%apbc,parm%ibrav)
          dd=SQRT(dx*dx+dy*dy+dz*dz)
          IF (dd .LT. 1.e-3_real_8)&
               CALL stopgm('HYD2_ION_DISTANCE','H/O position overlap',& 
               __LINE__,__FILE__)
          ! compute nH
          rnum = (dd/c1_rc(1))**nexp(1)
          rden = (dd/c1_rc(1))**(nexp(1)+mexp(1))
          fnum = 1-rnum
          fden = 1-rden
          IF (ABS(fden) .LT. 1.e-10_real_8) fden = 1.e-10_real_8
          nH= nh + fnum/fden

          ! compute der_nH_dr  index dependent
          ff = -(rnum*REAL(nexp(1),kind=real_8)&
               -fnum*rden*REAL(nexp(1)+mexp(1),kind=real_8)/fden)/&
               (fden*dd*dd)
          der_nH_dr(1,ii) = ff*dx
          der_nH_dr(2,ii) = ff*dy
          der_nH_dr(3,ii) = ff*dz

       ENDDO   ! Loop on the hydrogens

       ! jpark  compute r_k 
       dx2 = x0-xion
       dy2 = y0-yion
       dz2 = z0-zion
       IF (.NOT. isos1%tisos) CALL pbc(dx2,dy2,dz2,dx2,dy2,dz2,1,parm%apbc,parm%ibrav)
       r_k = SQRT(dx2*dx2+dy2*dy2+dz2*dz2)


       ! jpark compute der_r_k_dr with respect to Ok
       der_r_k_dr(1) = dx2/r_k
       der_r_k_dr(2) = dy2/r_k
       der_r_k_dr(3) = dz2/r_k

       ! exponential factor exp{L*nH}
       exp_fact = EXP(l*nH)

       ! jpark  for summation of expnential factor
       sum_exp_fact   = sum_exp_fact + exp_fact
       sum_r_exp_fact = sum_r_exp_fact + r_k*exp_fact

       ! dbg
       ! write(6,'(A,I3,2(A,f10.5))') ' ox ' , IA, ' nh ',nH,
       ! c                                 ' r_k ',  r_k
       ! dbg


       ! dbg
       ! write(6,'(2(A,1E15.6))') ' r exp ' , r_k*exp_fact ,
       ! c                  ' exp ',exp_fact
       ! dbg



       ! Derivatives:
       der_fact = l*r_k*exp_fact
       der_fact2 = l*exp_fact
       ii = 0
       ! Loop over the hydrogens for the derivatives
       DO iiat2=iini2,iend2
          ii = ii + 1
          IF (lsk(1,iiat2).NE.0) THEN
             k=lsk(1,iiat2)
             an(k)=an(k)+ der_fact*der_nH_dr(1,ii)
             anf1(k)=anf1(k)+ der_fact2*der_nH_dr(1,ii)
          ENDIF
          IF (lsk(2,iiat2).NE.0) THEN
             k=lsk(2,iiat2)
             an(k)=an(k)+ der_fact*der_nH_dr(2,ii)
             anf1(k)=anf1(k)+ der_fact2*der_nH_dr(2,ii)
          ENDIF
          IF (lsk(3,iiat2).NE.0) THEN
             k=lsk(3,iiat2)
             an(k)=an(k)+ der_fact*der_nH_dr(3,ii)
             anf1(k)=anf1(k)+ der_fact2*der_nH_dr(3,ii)
          ENDIF
          IF (l1.NE.0) an(l1)=an(l1)-der_fact*der_nH_dr(1,ii)
          IF (l2.NE.0) an(l2)=an(l2)-der_fact*der_nH_dr(2,ii)
          IF (l3.NE.0) an(l3)=an(l3)-der_fact*der_nH_dr(3,ii)
          IF (l1.NE.0) anf1(l1)=anf1(l1)-der_fact2*der_nH_dr(1,ii)
          IF (l2.NE.0) anf1(l2)=anf1(l2)-der_fact2*der_nH_dr(2,ii)
          IF (l3.NE.0) anf1(l3)=anf1(l3)-der_fact2*der_nH_dr(3,ii)
       ENDDO

       ! compute term1 der_r_k_dr * exp_fact
       IF (l1.NE.0) an(l1)=an(l1)+der_r_k_dr(1)*exp_fact
       IF (l2.NE.0) an(l2)=an(l2)+der_r_k_dr(2)*exp_fact
       IF (l3.NE.0) an(l3)=an(l3)+der_r_k_dr(3)*exp_fact
       IF (l11.NE.0) an(l11)=an(l11)-der_r_k_dr(1)*exp_fact
       IF (l22.NE.0) an(l22)=an(l22)-der_r_k_dr(2)*exp_fact
       IF (l33.NE.0) an(l33)=an(l33)-der_r_k_dr(3)*exp_fact

100    CONTINUE
       ! write(6,*) 'O_k r_k nH e(nH)', IA, r_k, nH, exp_fact

    ENDDO    ! Loop over the oxygens

    cval = sum_r_exp_fact/sum_exp_fact

    ! compute term2 r_k*ANF0
    !$omp parallel do private(II) shared(CVAL)
    DO ii = 1,cotc0%nodim
       anf1(ii) = anf1(ii)*cval
    ENDDO

    ! write(6,*) 'HY_VAL ',CVAL
    der_fact = 1.0_real_8/sum_exp_fact
    ! Final factor for the derivatives
    !$omp parallel do private(IAT)
    DO iat = 1,cotc0%nodim
       an(iat) = (an(iat) - anf1(iat))*der_fact
    ENDDO
    ! stop
    DEALLOCATE(der_nH_dr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(anf1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==---------------------------------------------------------==
    RETURN
  END SUBROUTINE hyd2_ion_distance
  ! ==================================================================
  SUBROUTINE hyd_ion_distance(isp1,isp2,indion,nexp,mexp,c1_rc,l,&
       tscr,an,lsk,cval)
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: isp1, isp2, indion, nexp(3), &
                                                mexp(3)
    REAL(real_8)                             :: c1_rc(3), l, tscr(3,*), &
                                                an(cotc0%nodim)
    INTEGER                                  :: lsk(3,ions1%nat)
    REAL(real_8)                             :: cval

    CHARACTER(*), PARAMETER :: procedureN = 'hyd_ion_distance'

    INTEGER                                  :: ia, iat, iend1, iend2, ierr, &
                                                ii, iiat1, iiat2, iini1, &
                                                iini2, is, k, l1, l11, l2, &
                                                l22, l3, l33
    REAL(real_8) :: dd, der_fact, der_fact2, der_r_k_dr(3), dx, dx2, dy, dy2, &
      dz, dz2, exp_fact, fden, ff, fnum, nH, nO, r_k, rden, rnum, &
      sum_exp_fact, sum_r_exp_fact, x0, xion, xx(3), xx0(3), y0, yion, z0, &
      zion
    REAL(real_8), ALLOCATABLE                :: anf1(:), der_nH_dr(:,:), &
                                                der_nO_dr(:,:)

! Variables
! jpark  
! jpark
! ==--------------------------------------------------------------==

    cval = 0.0_real_8
    sum_exp_fact = 0.0_real_8
    sum_r_exp_fact = 0.0_real_8

    ALLOCATE(der_nH_dr(3,ions0%na(isp2)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(der_nH_dr)!,3*ions0%na(isp2))
    ALLOCATE(der_nO_dr(3,ions0%na(isp1)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(der_nO_dr)!,3*ions0%na(isp1))

    ALLOCATE(anf1(cotc0%nodim),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(an)!,cotc0%nodim)
    ! jpark

    CALL zeroing(anf1)!,cotc0%nodim)
    ! Atoms indexes
    iini1 = 1
    DO is = 1,isp1-1
       iini1=iini1+ions0%na(is)
    ENDDO
    iend1 = iini1 + ions0%na(isp1) -1
    iini2 = 1
    DO is = 1,isp2-1
       iini2=iini2+ions0%na(is)
    ENDDO
    iend2 = iini2 + ions0%na(isp2) -1

    ! jpark     SET ION POSITION
    CALL fillc(indion,tscr,xx0)
    xion = xx0(1)
    yion = xx0(2)
    zion = xx0(3)

    l11 = lsk(1,indion)
    l22 = lsk(2,indion)
    l33 = lsk(3,indion)


    ! First Loop on the oxygens
    DO ia=iini1,iend1

       IF (ia .EQ. indion) GOTO 100
       nH=0.0_real_8
       nO=0.0_real_8
       CALL zeroing(der_nH_dr)!,3*ions0%na(isp2))
       IF (ia.LE.ions1%nat) THEN
          l1 = lsk(1,ia)
          l2 = lsk(2,ia)
          l3 = lsk(3,ia)
       ENDIF
       CALL fillc(ia,tscr,xx)
       x0 = xx(1)
       y0 = xx(2)
       z0 = xx(3)

       ! Loop over the hydrogens
       ii = 0
       DO iiat2=iini2,iend2
          ii = ii + 1
          dx=tscr(1,iiat2)-x0
          dy=tscr(2,iiat2)-y0
          dz=tscr(3,iiat2)-z0
          IF (.NOT. isos1%tisos) CALL pbc(dx,dy,dz,dx,dy,dz,1,parm%apbc,parm%ibrav)
          dd=SQRT(dx*dx+dy*dy+dz*dz)
          IF (dd .LT. 1.e-3_real_8)&
               CALL stopgm('HYD_PRESENCE','H/O position overlap',& 
               __LINE__,__FILE__)
          ! compute nH
          rnum = (dd/c1_rc(1))**nexp(1)
          rden = (dd/c1_rc(1))**(nexp(1)+mexp(1))
          fnum = 1-rnum
          fden = 1-rden
          IF (ABS(fden) .LT. 1.e-10_real_8) fden = 1.e-10_real_8
          nH= nh + fnum/fden

          ! compute der_nH_dr  index dependent
          ff = -(rnum*REAL(nexp(1),kind=real_8)&
               -fnum*rden*REAL(nexp(1)+mexp(1),kind=real_8)/fden)/&
               (fden*dd*dd)
          der_nH_dr(1,ii) = ff*dx
          der_nH_dr(2,ii) = ff*dy
          der_nH_dr(3,ii) = ff*dz

       ENDDO   ! Loop on the hydrogens

       ! jpark  compute r_k
       dx2 = x0-xion
       dy2 = y0-yion
       dz2 = z0-zion
       IF (.NOT. isos1%tisos) CALL pbc(dx2,dy2,dz2,dx2,dy2,dz2,1,parm%apbc,parm%ibrav)
       r_k = SQRT(dx2*dx2+dy2*dy2+dz2*dz2)


       ! jpark compute der_r_k_dr with respect to Ok
       der_r_k_dr(1) = dx2/r_k
       der_r_k_dr(2) = dy2/r_k
       der_r_k_dr(3) = dz2/r_k

       ii = 0
       ! Second Loop over the oxygens
       DO iiat1 = iini1,iend1
          ii = ii + 1
          dx=tscr(1,iiat1)-x0
          dy=tscr(2,iiat1)-y0
          dz=tscr(3,iiat1)-z0
          IF (.NOT. isos1%tisos) CALL pbc(dx,dy,dz,dx,dy,dz,1,parm%apbc,parm%ibrav)
          dd=SQRT(dx*dx+dy*dy+dz*dz)
          IF (dd .LT. 1.e-3_real_8) GOTO 10
          ! compute nO
          rnum = (dd/c1_rc(2))**nexp(2)
          rden = (dd/c1_rc(2))**(nexp(2)+mexp(2))
          fnum = 1-rnum
          fden = 1-rden
          IF (ABS(fden) .LT. 1.e-10_real_8) fden = 1.e-10_real_8
          nO= no + fnum/fden

          ! compute der_nO_dr  
          ff = -(rnum*REAL(nexp(2),kind=real_8)&
               -fnum*rden*REAL(nexp(2)+mexp(2),kind=real_8)/fden)/&
               (fden*dd*dd)
          der_nO_dr(1,ii) = ff*dx
          der_nO_dr(2,ii) = ff*dy
          der_nO_dr(3,ii) = ff*dz
10        CONTINUE
       ENDDO ! Second Loop over the oxygens

       ! dbg
       ! write(6,'(A,I3,3(A,f10.5))') ' ox ' , IA, ' nh ',nH,
       ! c                               ' nO ', nO, '   r_k ',  r_k
       ! dbg


       ! exponential factor exp{L*nH}

       exp_fact = EXP(l*nH*nO)

       ! jpark  for summation of expnential factor
       sum_exp_fact   = sum_exp_fact + exp_fact
       sum_r_exp_fact = sum_r_exp_fact + r_k*exp_fact

       ! dbg
       ! write(6,'(2(A,1E15.6))') ' r exp ' , r_k*exp_fact ,
       ! c                  ' exp ',exp_fact
       ! dbg

       ! Derivatives:
       der_fact = l*r_k*exp_fact*nH
       der_fact2 = l*exp_fact*nH
       ii = 0
       ! Loop over the hydrogens for the derivatives
       DO iiat1=iini1,iend1
          ii = ii + 1
          IF (lsk(1,iiat1).NE.0) THEN
             k=lsk(1,iiat1)
             an(k)=an(k)+ der_fact*der_nO_dr(1,ii)
             anf1(k)=anf1(k)+ der_fact2*der_nO_dr(1,ii)
          ENDIF
          IF (lsk(2,iiat1).NE.0) THEN
             k=lsk(2,iiat1)
             an(k)=an(k)+ der_fact*der_nO_dr(2,ii)
             anf1(k)=anf1(k)+ der_fact2*der_nO_dr(2,ii)
          ENDIF
          IF (lsk(3,iiat1).NE.0) THEN
             k=lsk(3,iiat1)
             an(k)=an(k)+ der_fact*der_nO_dr(3,ii)
             anf1(k)=anf1(k)+ der_fact2*der_nO_dr(3,ii)
          ENDIF
          IF (l1.NE.0) an(l1)=an(l1)-der_fact*der_nO_dr(1,ii)
          IF (l2.NE.0) an(l2)=an(l2)-der_fact*der_nO_dr(2,ii)
          IF (l3.NE.0) an(l3)=an(l3)-der_fact*der_nO_dr(3,ii)
          IF (l1.NE.0) anf1(l1)=anf1(l1)-der_fact2*der_nO_dr(1,ii)
          IF (l2.NE.0) anf1(l2)=anf1(l2)-der_fact2*der_nO_dr(2,ii)
          IF (l3.NE.0) anf1(l3)=anf1(l3)-der_fact2*der_nO_dr(3,ii)
       ENDDO


       ! Derivatives:
       der_fact = l*r_k*exp_fact*nO
       der_fact2 = l*exp_fact*nO
       ii = 0
       ! Loop over the hydrogens for the derivatives
       DO iiat2=iini2,iend2
          ii = ii + 1
          IF (lsk(1,iiat2).NE.0) THEN
             k=lsk(1,iiat2)
             an(k)=an(k)+ der_fact*der_nH_dr(1,ii)
             anf1(k)=anf1(k)+ der_fact2*der_nH_dr(1,ii)
          ENDIF
          IF (lsk(2,iiat2).NE.0) THEN
             k=lsk(2,iiat2)
             an(k)=an(k)+ der_fact*der_nH_dr(2,ii)
             anf1(k)=anf1(k)+ der_fact2*der_nH_dr(2,ii)
          ENDIF
          IF (lsk(3,iiat2).NE.0) THEN
             k=lsk(3,iiat2)
             an(k)=an(k)+ der_fact*der_nH_dr(3,ii)
             anf1(k)=anf1(k)+ der_fact2*der_nH_dr(3,ii)
          ENDIF
          IF (l1.NE.0) an(l1)=an(l1)-der_fact*der_nH_dr(1,ii)
          IF (l2.NE.0) an(l2)=an(l2)-der_fact*der_nH_dr(2,ii)
          IF (l3.NE.0) an(l3)=an(l3)-der_fact*der_nH_dr(3,ii)
          IF (l1.NE.0) anf1(l1)=anf1(l1)-der_fact2*der_nH_dr(1,ii)
          IF (l2.NE.0) anf1(l2)=anf1(l2)-der_fact2*der_nH_dr(2,ii)
          IF (l3.NE.0) anf1(l3)=anf1(l3)-der_fact2*der_nH_dr(3,ii)
       ENDDO

       ! compute term1 der_r_k_dr * exp_fact
       IF (l1.NE.0) an(l1)=an(l1)+der_r_k_dr(1)*exp_fact
       IF (l2.NE.0) an(l2)=an(l2)+der_r_k_dr(2)*exp_fact
       IF (l3.NE.0) an(l3)=an(l3)+der_r_k_dr(3)*exp_fact
       IF (l11.NE.0) an(l11)=an(l11)-der_r_k_dr(1)*exp_fact
       IF (l22.NE.0) an(l22)=an(l22)-der_r_k_dr(2)*exp_fact
       IF (l33.NE.0) an(l33)=an(l33)-der_r_k_dr(3)*exp_fact

100    CONTINUE
       ! write(6,*) 'O_k r_k nH e(nH)', IA, r_k, nH, exp_fact  
    ENDDO    ! Loop over the oxygens

    cval = sum_r_exp_fact/sum_exp_fact

    ! write(6,'((A,f15.6))') ' final cv ' , cval
    ! stop
    ! dbg

    ! compute term2 r_k*ANF0
    !$omp parallel do private(II) shared(CVAL)
    DO ii = 1,cotc0%nodim
       anf1(ii) = anf1(ii)*cval
    ENDDO

    ! write(6,*) 'HY_VAL ',CVAL
    der_fact = 1.0_real_8/sum_exp_fact
    ! Final factor for the derivatives
    !$omp parallel do private(IAT)
    DO iat = 1,cotc0%nodim
       an(iat) = (an(iat) - anf1(iat))*der_fact
    ENDDO

    DEALLOCATE(der_nH_dr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(der_nO_dr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(anf1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==---------------------------------------------------------==
    RETURN
  END SUBROUTINE hyd_ion_distance
  ! =============================================================
  SUBROUTINE funcl(fd,dd,r0,x1,x2,dr)
    ! ==---------------------------------------------------------==
    ! == DISTANCE FROM THE SPIN CENTER                           ==
    ! ==---------------------------------------------------------==
    ! FUNCTION: R - R0
    ! Argument variables
    REAL(real_8)                             :: fd, dd, r0, x1(3), x2(3), &
                                                dr(6)

    REAL(real_8)                             :: dx, dy, dz, epsilon

! ==---------------------------------------------------------==

    dx=x1(1)-x2(1)
    dy=x1(2)-x2(2)
    dz=x1(3)-x2(3)
    ! write(6,*) X1
    ! write(6,*) X2
    ! stop
    CALL pbc(dx,dy,dz,dx,dy,dz,1,parm%apbc,parm%ibrav)
    dd=SQRT(dx*dx+dy*dy+dz*dz)

    fd = dd-r0

    epsilon = 1.e-5_real_8
    ! ==---------DERIVATIVE OF THE ATOM SPIN DISTANCE CONSTR.-------==
    IF (dd.LT.epsilon) THEN
       dr(1) = 0._real_8
       dr(2) = 0._real_8
       dr(3) = 0._real_8
       dr(4) = 0._real_8
       dr(5) = 0._real_8
       dr(6) = 0._real_8
    ELSE
       dr(1) = dx/dd
       dr(2) = dy/dd
       dr(3) = dz/dd
       dr(4) = 0._real_8
       dr(5) = 0._real_8
       dr(6) = 0._real_8
    ENDIF

    ! ==---------------------------------------------------------==
    RETURN
  END SUBROUTINE funcl
  ! =============================================================
  SUBROUTINE volume_cv(vol_ist,det_cell)
    ! ==---------------------------------------------------------==
    ! == VOLUME OF THE  CELL                                     ==
    ! ==---------------------------------------------------------==
    REAL(real_8)                             :: vol_ist, det_cell(3,3)

    INTEGER                                  :: i, i1, i2, j, j1, j2
    REAL(real_8)                             :: aa1(3), aa2(3), aa3(3), &
                                                fv(3,3), uvol

! ==---------------------------------------------------------==

    DO i = 1,3
       aa1(i) = metr_com%ht(1,i)
       aa2(i) = metr_com%ht(2,i)
       aa3(i) = metr_com%ht(3,i)
    ENDDO
    CALL  omegagen(aa1,aa2,aa3,vol_ist)

    uvol = 1.0_real_8/vol_ist

    i1 = 2
    i2 = 3
    j1 = i1
    j2 = i2

    DO i = 1,3
       DO j = 1,3
          fv(i,j) = (metr_com%ht(i1,j1)*metr_com%ht(i2,j2)-metr_com%ht(i2,j1)*metr_com%ht(i1,j2))
          j1 = j2
          j2 = j
       ENDDO
       j1 = 2
       j2 = 3
       i1 = i2
       i2 = i
    ENDDO

    CALL dgemm('T','N',3,3,3,1._real_8,fv,3,metr_com%ht,3,0._real_8,det_cell,3)

    CALL dscal(9,uvol,det_cell,1)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE volume_cv
  ! ==================================================================
  SUBROUTINE side_cv(lat_ist,side,det_cell)
    ! ==--------------------------------------------------------------==
    ! == SIDE LENGTH OF THE  CELL                                     ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: lat_ist
    INTEGER                                  :: side
    REAL(real_8)                             :: det_cell(3,3)

    INTEGER                                  :: i
    REAL(real_8)                             :: aa1(3), aa2(3), aa3(3), &
                                                fv(3,3), ulat, uvol, vol

! ==--------------------------------------------------------------==

    DO i = 1,3
       aa1(i) = metr_com%ht(1,i)
       aa2(i) = metr_com%ht(2,i)
       aa3(i) = metr_com%ht(3,i)
    ENDDO
    CALL  omegagen(aa1,aa2,aa3,vol)

    uvol = 1.0_real_8/vol

    lat_ist = 0.0_real_8
    DO i = 1,3
       lat_ist = lat_ist +  metr_com%ht(side,i)* metr_com%ht(side,i)
    ENDDO
    lat_ist = SQRT(lat_ist)
    ulat = 1.0_real_8/lat_ist
    DO i = 1,3
       fv(side,i) = ulat*metr_com%ht(side,i)
    ENDDO

    CALL dgemm('T','N',3,3,3,1._real_8,fv,3,metr_com%ht,3,0._real_8,det_cell,3)

    CALL dscal(9,uvol,det_cell,1)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE side_cv
  ! ==================================================================
  SUBROUTINE angle_cv(angle_ist,nangle,det_cell)
    ! ==--------------------------------------------------------------==
    ! == ANGLE OF THE  CELL                                           ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: angle_ist
    INTEGER                                  :: nangle
    REAL(real_8)                             :: det_cell(3,3)

    INTEGER                                  :: i, i1, i2
    REAL(real_8)                             :: aa1(3), aa2(3), aa3(3), &
                                                fv(3,3), lat_ist1, lat_ist2, &
                                                ul1, ul1l2, ul2, uvol, vol

! ==--------------------------------------------------------------==

    DO i = 1,3
       aa1(i) = metr_com%ht(1,i)
       aa2(i) = metr_com%ht(2,i)
       aa3(i) = metr_com%ht(3,i)
    ENDDO
    CALL  omegagen(aa1,aa2,aa3,vol)

    uvol = 1.0_real_8/vol

    IF (nangle .EQ. 1) THEN
       i1 = 2
       i2 = 3
    ELSEIF (nangle .EQ. 2) THEN
       i1 = 3
       i2 = 1
    ELSE
       i1 = 1
       i2 = 2
    ENDIF

    lat_ist1=metr_com%ht(i1,1)*metr_com%ht(i1,1)+metr_com%ht(i1,2)*metr_com%ht(i1,2)+metr_com%ht(i1,3)*metr_com%ht(i1,3)
    lat_ist1=SQRT(lat_ist1)

    lat_ist2=metr_com%ht(i2,1)*metr_com%ht(i2,1)+metr_com%ht(i2,2)*metr_com%ht(i2,2)+metr_com%ht(i2,3)*metr_com%ht(i2,3)
    lat_ist2=SQRT(lat_ist2)

    ul1   = 1.0_real_8/lat_ist1
    ul2   = 1.0_real_8/lat_ist2
    ul1l2 = ul1*ul2

    angle_ist = metr_com%ht(i1,1)*metr_com%ht(i2,1)+metr_com%ht(i1,2)*metr_com%ht(i2,2)+metr_com%ht(i1,3)*metr_com%ht(i2,3)

    fv(i1,1) = ul1l2*(metr_com%ht(i2,1)-angle_ist*ul1l2*metr_com%ht(i1,1))
    fv(i2,1) = ul1l2*(metr_com%ht(i1,1)-angle_ist*ul1l2*metr_com%ht(i2,1))

    fv(i1,2) = ul1l2*(metr_com%ht(i2,2)-angle_ist*ul1l2*metr_com%ht(i1,2))
    fv(i2,2) = ul1l2*(metr_com%ht(i1,2)-angle_ist*ul1l2*metr_com%ht(i2,2))

    fv(i1,3) = ul1l2*(metr_com%ht(i2,3)-angle_ist*ul1l2*metr_com%ht(i1,3))
    fv(i2,3) = ul1l2*(metr_com%ht(i1,3)-angle_ist*ul1l2*metr_com%ht(i2,3))

    angle_ist = angle_ist*ul1l2

    CALL dgemm('T','N',3,3,3,1._real_8,fv,3,metr_com%ht,3,0._real_8,det_cell,3)

    CALL dscal(9,uvol,det_cell,1)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE angle_cv
  ! ==================================================================
  SUBROUTINE hyd_cutoff(intpar,realpar,numO,numH,numL,&
       iatom,tscr,an,lsk,cval)
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: intpar(15)
    REAL(real_8)                             :: realpar(10)
    INTEGER                                  :: numO, numH, numL, iatom(*)
    REAL(real_8)                             :: tscr(3,*), an(cotc0%nodim)
    INTEGER                                  :: lsk(3,ions1%nat)
    REAL(real_8)                             :: cval

    CHARACTER(*), PARAMETER                  :: procedureN = 'hyd_cutoff'

    INTEGER :: ial, iat1, iat2, iatl, iend1, iend2, ierr, ii, iiat1, iini1, &
      iini2, is, isp1, isp2, k, l1, l2, l3, pH, pL, pM, pO, qH, qL, qM, qO
    REAL(real_8) :: dd, DER_FACT_nHi, DER_FACT_nOi, der_fc_nL, der_log_sum, &
      der_M_nHi, der_M_nOi, dx, dy, dz, EXP_i, fc_nL, fden, ff, fnum, lambda, &
      log_sum, M_nHi, M_nOi, nH0, nHi, nHl, nL, nL0, nO0, nOi, rden, rH0, &
      rnum, rO0, sum_e, x0, xx(3), y0, z0
    REAL(real_8), ALLOCATABLE                :: der_nHi_dr(:,:), &
                                                der_nHl_dr(:,:,:), &
                                                der_nOi_dr(:,:)

    cval = 0.0_real_8

    isp1 = intpar(1)
    isp2 = intpar(2)
    pH = intpar(3)
    qH = intpar(4)
    pO = intpar(5)
    qO = intpar(6)
    pM = intpar(7)
    qM = intpar(8)
    pL = intpar(9)
    qL = intpar(10)

    rH0 = realpar(1)
    rO0 = realpar(2)
    nH0 = realpar(3)
    nO0 = realpar(4)
    nL0 = realpar(5)
    lambda = realpar(6)

    ALLOCATE(der_nHi_dr(3,numH),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(der_nOi_dr(3,numO),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(der_nHl_dr(3,numH,numL),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(an)!,cotc0%nodim)
    ! Atoms indexes
    iini1 = 1
    DO is = 1,isp1-1
       iini1=iini1+ions0%na(is)
    ENDDO
    iend1 = iini1 + ions0%na(isp1) -1
    iini2 = 1
    DO is = 1,isp2-1
       iini2=iini2+ions0%na(is)
    ENDDO
    iend2 = iini2 + ions0%na(isp2) -1

    ! Compute nL and the derivadives wrt RHk: der_nHl_dr
    nL=0
    fc_nL=0
    CALL zeroing(der_nHl_dr)!,3*ions0%na(isp2)*numL)


    DO ial = 1,numL

       nHl=0.0_real_8

       iatl = iatom(ial)
       x0 = tscr(1,iatl)
       y0 = tscr(2,iatl)
       z0 = tscr(3,iatl)
       ! write(6,*) IATL,X0,Y0,Z0

       ii = 0
       DO iat2 = iini2,iend2

          ii = ii + 1
          dx=tscr(1,iat2)-x0
          dy=tscr(2,iat2)-y0
          dz=tscr(3,iat2)-z0
          IF (.NOT. isos1%tisos) CALL pbc(dx,dy,dz,dx,dy,dz,1,parm%apbc,parm%ibrav)
          dd=SQRT(dx*dx+dy*dy+dz*dz)
          IF (dd .LT. 1.e-3_real_8)&
               CALL stopgm('HYD_CUTOFF','H/L position overlap',& 
               __LINE__,__FILE__)

          ! write(6,'(I4,4f14.6)') IAT2 ,DX,DY,DZ,DD

          ! Compute the contribution of IAT2 to nHl

          rnum = (dd/rH0)**pH
          rden = (dd/rH0)**(pH+qH)

          fnum = 1.0_real_8-rnum
          fden = 1.0_real_8-rden
          IF (ABS(fden) .LT. 1.e-10_real_8) fden = 1.e-10_real_8
          nHl = nhl + fnum/fden
          ! dbg     write(6,'(I4,4f14.6)') IAT2 ,FNUM,FDEN,FNUM/FDEN,nHl

          ! Compute the derivative of this contribution wrt the IAT2 coordinates
          ff = -(rnum*REAL(pH,kind=real_8)-fnum*rden*REAL(ph+qH,kind=real_8)/fden)/(fden*dd*dd)

          der_nHl_dr(1,ii,ial) = ff*dx
          der_nHl_dr(2,ii,ial) = ff*dy
          der_nHl_dr(3,ii,ial) = ff*dz

       ENDDO! IAT2
       nL = nl + nHl
       ! dbg
       ! write(6,'(A4,I4,I4,f14.6)') 'nHl ', ial, iatl, nHl
       ! dbg
    ENDDO  ! IAL

    ! Compute the Cutoff dunction and its derivative with respect to nL
    rnum = (nL/nL0)**pL
    rden = (nL/nL0)**(pL+qL)
    fnum = 1.0_real_8-rnum
    fden = 1.0_real_8-rden
    IF (ABS(fden) .LT. 1.e-10_real_8) fden = 1.e-10_real_8
    fc_nL = fnum/fden
    ! dbg
    ! write(6,*) RNUM, RDEN,FNUM/FDEN

    ff = -(rnum*REAL(pL,kind=real_8)-fnum*rden*REAL(pl+qL,kind=real_8)/fden)/(fden*nL)

    der_fc_nL=ff
    ! dbg
    ! write(6,'(A6,3f14.6)') 'fc_nL ', nL, fc_nL, der_fc_nL
    ! dbg

    log_sum = 0.0_real_8
    sum_e   = 0.0_real_8
    EXP_i   = 0.0_real_8

    ! First Loop over the oxygens
    DO iat1=iini1,iend1
       nHi=0.0_real_8
       nOi=0.0_real_8
       M_nHi=0.0_real_8
       M_nOi=0.0_real_8

       CALL zeroing(der_nHi_dr)!,3*ions0%na(isp2))
       CALL zeroing(der_nOi_dr)!,3*ions0%na(isp1))
       der_M_nHi =0.0_real_8
       der_M_nOi =0.0_real_8

       IF (iat1.LE.ions1%nat) THEN
          l1 = lsk(1,iat1)
          l2 = lsk(2,iat1)
          l3 = lsk(3,iat1)
       ENDIF
       CALL fillc(iat1,tscr,xx)
       x0 = xx(1)
       y0 = xx(2)
       z0 = xx(3)

       ! Loop over the hydrogens
       ii = 0
       DO iat2=iini2,iend2
          ii = ii + 1
          dx=tscr(1,iat2)-x0
          dy=tscr(2,iat2)-y0
          dz=tscr(3,iat2)-z0
          IF (.NOT. isos1%tisos) CALL pbc(dx,dy,dz,dx,dy,dz,1,parm%apbc,parm%ibrav)
          dd=SQRT(dx*dx+dy*dy+dz*dz)
          IF (dd .LT. 1.e-3_real_8)&
               CALL stopgm('HYD_PRESENCE','H/O position overlap',& 
               __LINE__,__FILE__)
          ! compute nHi
          rnum = (dd/rH0)**pH
          rden = (dd/rH0)**(pH+qH)
          fnum = 1-rnum
          fden = 1-rden
          IF (ABS(fden) .LT. 1.e-10_real_8) fden = 1.e-10_real_8
          nHi= nhi + fnum/fden
          ! write(6,*) 'H ', IIAT2, DD, FNUM/FDEN
          ! compute der_nH_dr  index dependent
          ff = -(rnum*REAL(pH,kind=real_8)&
               -fnum*rden*REAL(pH+qH,kind=real_8)/fden)/&
               (fden*dd*dd)
          der_nHi_dr(1,ii) = ff*dx
          der_nHi_dr(2,ii) = ff*dy
          der_nHi_dr(3,ii) = ff*dz
       ENDDO   ! Loop on the hydrogens


       ! compute M_nHi
       rnum = (nHi/nH0)**pM
       rden = (nHi/nH0)**(pM+qM)
       fnum = 1-rnum
       fden = 1-rden
       IF (ABS(fden) .LT. 1.e-10_real_8) fden = 1.e-10_real_8
       M_nHi= 1.0_real_8 - fnum/fden

       ! compute der_M_nHi
       der_M_nHi = (REAL(pM,kind=real_8)*rnum-&
            REAL(pM+qM,kind=real_8)*fnum*rden/fden)/(fden*nHi)

       ii = 0
       ! Second Loop over the oxygens
       DO iiat1 = iini1,iend1
          ii = ii + 1
          dx=tscr(1,iiat1)-x0
          dy=tscr(2,iiat1)-y0
          dz=tscr(3,iiat1)-z0
          IF (.NOT. isos1%tisos) CALL pbc(dx,dy,dz,dx,dy,dz,1,parm%apbc,parm%ibrav)
          dd=SQRT(dx*dx+dy*dy+dz*dz)
          IF (dd .LT. 1.e-3_real_8) GOTO 10
          ! compute nO
          rnum = (dd/rO0)**pO
          rden = (dd/rO0)**(pO+qO)
          fnum = 1-rnum
          fden = 1-rden
          IF (ABS(fden) .LT. 1.e-10_real_8) fden = 1.e-10_real_8
          nOi= noi + fnum/fden

          ! compute der_nO_dr  
          ff = -(rnum*REAL(pO,kind=real_8)&
               -fnum*rden*REAL(pO+qO,kind=real_8)/fden)/&
               (fden*dd*dd)
          der_nOi_dr(1,ii) = ff*dx
          der_nOi_dr(2,ii) = ff*dy
          der_nOi_dr(3,ii) = ff*dz
10        CONTINUE
       ENDDO ! Second Loop over the oxygens

       ! compute M_nOi
       rnum = (nOi/nO0)**pM
       rden = (nOi/nO0)**(pM+qM)
       fnum = 1-rnum
       fden = 1-rden
       IF (ABS(fden) .LT. 1.e-10_real_8) fden = 1.e-10_real_8
       M_nOi= 1.0_real_8 - fnum/fden


       ! compute der_M_nOi
       der_M_nOi = (REAL(pM,kind=real_8)*rnum-&
            REAL(pM+qM,kind=real_8)*fnum*rden/fden)/(fden*nOi)

       ! exponential factor exp{L*M_nHi*M_nOi}
       EXP_i = EXP(lambda* M_nHi * M_nOi)
       DER_FACT_nOi = lambda*EXP_i*M_nHi*der_M_nOi
       DER_FACT_nHi = lambda*EXP_i*M_nOi*der_M_nHi
       ! exp_i = exp(Lambda* M_nHi * nOi)
       ! DER_FACT_nOi = Lambda*exp_i*M_nHi
       ! DER_FACT_nHi = Lambda*exp_i*nOi*der_M_nHi


       ! dbg
       ! write(6,'(A6,I4,3f14.6)') 'M_nOi ', iat1, nOi, M_nOi, der_M_nOi
       ! write(6,'(A6,I4,3f14.6)') 'M_nHi ', iat1, nHi, M_nHi, der_M_nHi
       ! write(6,'(A6,I4,3f14.6)') 'EXP ', iat1, exp_i,
       ! &                           DER_FACT_nOi, DER_FACT_nHi
       ! dbg
       ii = 0
       DO iiat1 = iini1,iend1
          ii = ii + 1
          ! Derivatives:
          IF (lsk(1,iiat1).NE.0) THEN
             k=lsk(1,iiat1)
             an(k)=an(k)+ DER_FACT_nOi*der_nOi_dr(1,ii)
          ENDIF
          IF (lsk(2,iiat1).NE.0) THEN
             k=lsk(2,iiat1)
             an(k)=an(k)+ DER_FACT_nOi*der_nOi_dr(2,ii)
          ENDIF
          IF (lsk(3,iiat1).NE.0) THEN
             k=lsk(3,iiat1)
             an(k)=an(k)+ DER_FACT_nOi*der_nOi_dr(3,ii)
          ENDIF
          IF (l1.NE.0) an(l1)=an(l1)-DER_FACT_nOi*der_nOi_dr(1,ii)
          IF (l2.NE.0) an(l2)=an(l2)-DER_FACT_nOi*der_nOi_dr(2,ii)
          IF (l3.NE.0) an(l3)=an(l3)-DER_FACT_nOi*der_nOi_dr(3,ii)

       ENDDO

       ii = 0
       ! Loop over the hydrogens for the derivatives
       DO iat2=iini2,iend2
          ii = ii + 1
          IF (lsk(1,iat2).NE.0) THEN
             k=lsk(1,iat2)
             an(k)=an(k)+ DER_FACT_nHi*der_nHi_dr(1,ii)
          ENDIF
          IF (lsk(2,iat2).NE.0) THEN
             k=lsk(2,iat2)
             an(k)=an(k)+ DER_FACT_nHi*der_nHi_dr(2,ii)
          ENDIF
          IF (lsk(3,iat2).NE.0) THEN
             k=lsk(3,iat2)
             an(k)=an(k)+ DER_FACT_nHi*der_nHi_dr(3,ii)
          ENDIF
          IF (l1.NE.0) an(l1)=an(l1)-DER_FACT_nHi*der_nHi_dr(1,ii)
          IF (l2.NE.0) an(l2)=an(l2)-DER_FACT_nHi*der_nHi_dr(2,ii)
          IF (l3.NE.0) an(l3)=an(l3)-DER_FACT_nHi*der_nHi_dr(3,ii)
       ENDDO

       ! Sum over the oxygens
       sum_e = sum_e + EXP_i

    ENDDO    ! Loop over the oxygens

    log_sum = 1.0_real_8/lambda*LOG(sum_e)
    der_log_sum = 1.0_real_8/(lambda*sum_e)

    ! Final factor for the derivatives
    DO iat1 = 1,cotc0%nodim
       an(iat1) = an(iat1)*der_log_sum*fc_nL
    ENDDO

    ! Second loop on the list numL
    DO ial = 1,numL

       iatl = iatom(ial)

       l1 = lsk(1,iatl)
       l2 = lsk(2,iatl)
       l3 = lsk(3,iatl)

       ii = 0
       ! Loop over the hydrogens for the derivatives
       DO iat2=iini2,iend2
          ii = ii + 1
          IF (lsk(1,iat2).NE.0) THEN
             k=lsk(1,iat2)
             an(k)=an(k)+ der_fc_nL*der_nHl_dr(1,ii,ial)*log_sum
          ENDIF
          IF (lsk(2,iat2).NE.0) THEN
             k=lsk(2,iat2)
             an(k)=an(k)+ der_fc_nL*der_nHl_dr(2,ii,ial)*log_sum
          ENDIF
          IF (lsk(3,iat2).NE.0) THEN
             k=lsk(3,iat2)
             an(k)=an(k)+ der_fc_nL*der_nHl_dr(3,ii,ial)*log_sum
          ENDIF
          IF (l1.NE.0) an(l1)=an(l1)-&
               der_fc_nL*der_nHl_dr(1,ii,ial)*log_sum
          IF (l2.NE.0) an(l2)=an(l2)-&
               der_fc_nL*der_nHl_dr(2,ii,ial)*log_sum
          IF (l3.NE.0) an(l3)=an(l3)-&
               der_fc_nL*der_nHl_dr(3,ii,ial)*log_sum
       ENDDO

    ENDDO  ! IAL

    cval = log_sum*fc_nL

    ! dbg
    ! write(6,'(A4,2f14.6)') 'CV ', LOG_SUM, CVAL
    ! dbg

    DEALLOCATE(der_nHi_dr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(der_nOi_dr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(der_nHl_dr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE hyd_cutoff
  ! ==--------------------------------------------------------------==
  ! ==--------------------------------------------------------------==
  ! ==================================================================
  SUBROUTINE cut_hyd_ion_distance(intpar,realpar,numO,numH,numL,&
       iatom,tscr,an,lsk,cval)
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: intpar(15)
    REAL(real_8)                             :: realpar(10)
    INTEGER                                  :: numO, numH, numL, iatom(*)
    REAL(real_8)                             :: tscr(3,*), an(cotc0%nodim)
    INTEGER                                  :: lsk(3,ions1%nat)
    REAL(real_8)                             :: cval

    CHARACTER(*), PARAMETER :: procedureN = 'cut_hyd_ion_distance'

    INTEGER :: ia, ial, iat, iat2, iatl, iend1, iend2, ierr, ii, iiat1, &
      iiat2, iini1, iini2, indion, is, isp1, isp2, k, l1, l11, l2, l22, l3, &
      l33, pH, pL, pO, qH, qL, qO
    REAL(real_8) :: dd, der_fact, der_fact2, der_fc_nL, der_r_k_dr(3), dx, &
      dx2, dy, dy2, dz, dz2, exp_fact, fc_nL, fden, ff, fnum, lambda, nHi, &
      nHl, nL, nL0, nOi, r_k, rden, rH0, rnum, rO0, sum_exp_fact, &
      sum_r_exp_fact, x0, xion, xx(3), xx0(3), y0, yion, z0, zion
    REAL(real_8), ALLOCATABLE                :: anf1(:), der_nH_dr(:,:), &
                                                der_nHl_dr(:,:,:), &
                                                der_nO_dr(:,:)

! Variables
! jpark  
! jpark
! ==--------------------------------------------------------------==

    cval = 0.0_real_8
    sum_exp_fact = 0.0_real_8
    sum_r_exp_fact = 0.0_real_8

    isp1 = intpar(1)
    isp2 = intpar(2)
    indion = intpar(3)
    pH = intpar(4)
    qH = intpar(5)
    pO = intpar(6)
    qO = intpar(7)
    ! pM = INTPAR(8)
    ! qM = INTPAR(9)
    pL = intpar(8)
    qL = intpar(9)

    rH0 = realpar(1)
    rO0 = realpar(2)
    ! nH0 = REALPAR(3)
    ! nO0 = REALPAR(4)
    nL0 = realpar(3)
    lambda = realpar(4)

    ALLOCATE(der_nH_dr(3,numH),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(der_nH_dr)!,3*numH)
    ALLOCATE(der_nO_dr(3,numO),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(der_nO_dr)!,3*numO)
    ALLOCATE(der_nHl_dr(3,numH,numL),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(der_nHl_dr)!,3*numH*numL)

    ALLOCATE(anf1(cotc0%nodim),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(an)!,cotc0%nodim)
    ! jpark

    CALL zeroing(anf1)!,cotc0%nodim)
    ! Atoms indexes
    iini1 = 1
    DO is = 1,isp1-1
       iini1=iini1+ions0%na(is)
    ENDDO
    iend1 = iini1 + ions0%na(isp1) -1
    iini2 = 1
    DO is = 1,isp2-1
       iini2=iini2+ions0%na(is)
    ENDDO
    iend2 = iini2 + ions0%na(isp2) -1

    ! jpark     SET ION POSITION
    CALL fillc(indion,tscr,xx0)
    xion = xx0(1)
    yion = xx0(2)
    zion = xx0(3)

    l11 = lsk(1,indion)
    l22 = lsk(2,indion)
    l33 = lsk(3,indion)

    ! Compute nL and the derivadives wrt RHk: der_nHl_dr
    nL=0
    fc_nL=0
    DO ial = 1,numL

       nHl=0.0_real_8

       iatl = iatom(ial)
       x0 = tscr(1,iatl)
       y0 = tscr(2,iatl)
       z0 = tscr(3,iatl)

       ii = 0
       DO iat2 = iini2,iend2

          ii = ii + 1
          dx=tscr(1,iat2)-x0
          dy=tscr(2,iat2)-y0
          dz=tscr(3,iat2)-z0
          IF (.NOT. isos1%tisos) CALL pbc(dx,dy,dz,dx,dy,dz,1,parm%apbc,parm%ibrav)
          dd=SQRT(dx*dx+dy*dy+dz*dz)
          IF (dd .LT. 1.e-3_real_8)&
               CALL stopgm('DIST_HYD_CUTOFF','H/L position overlap',& 
               __LINE__,__FILE__)

          ! Compute the contribution of IAT2 to nHl

          rnum = (dd/rH0)**pH
          rden = (dd/rH0)**(pH+qH)

          fnum = 1.0_real_8-rnum
          fden = 1.0_real_8-rden
          IF (ABS(fden) .LT. 1.e-10_real_8) fden = 1.e-10_real_8
          nHl = nhl + fnum/fden

          ! Compute the derivative of this contribution wrt the IAT2 coordinates
          ff = -(rnum*REAL(pH,kind=real_8)-fnum*rden*REAL(ph+qH,kind=real_8)/fden)/(fden*dd*dd)

          der_nHl_dr(1,ii,ial) = ff*dx
          der_nHl_dr(2,ii,ial) = ff*dy
          der_nHl_dr(3,ii,ial) = ff*dz

       ENDDO! IAT2
       nL = nl + nHl

    ENDDO  ! IAL

    ! Compute the Cutoff dunction and its derivative with respect to nL
    rnum = (nL/nL0)**pL
    rden = (nL/nL0)**(pL+qL)
    fnum = 1.0_real_8-rnum
    fden = 1.0_real_8-rden
    IF (ABS(fden) .LT. 1.e-10_real_8) fden = 1.e-10_real_8
    fc_nL = fnum/fden

    ff = -(rnum*REAL(pL,kind=real_8)-fnum*rden*REAL(pl+qL,kind=real_8)/fden)/(fden*nL)

    der_fc_nL=ff


    ! First Loop on the oxygens
    DO ia=iini1,iend1

       IF (ia .EQ. indion) GOTO 100
       nHi=0.0_real_8
       nOi=0.0_real_8
       CALL zeroing(der_nH_dr)!,3*ions0%na(isp2))
       IF (ia.LE.ions1%nat) THEN
          l1 = lsk(1,ia)
          l2 = lsk(2,ia)
          l3 = lsk(3,ia)
       ENDIF
       CALL fillc(ia,tscr,xx)
       x0 = xx(1)
       y0 = xx(2)
       z0 = xx(3)

       ! Loop over the hydrogens
       ii = 0
       DO iiat2=iini2,iend2
          ii = ii + 1
          dx=tscr(1,iiat2)-x0
          dy=tscr(2,iiat2)-y0
          dz=tscr(3,iiat2)-z0
          IF (.NOT. isos1%tisos) CALL pbc(dx,dy,dz,dx,dy,dz,1,parm%apbc,parm%ibrav)
          dd=SQRT(dx*dx+dy*dy+dz*dz)
          IF (dd .LT. 1.e-3_real_8)&
               CALL stopgm('HYD_PRESENCE','H/O position overlap',& 
               __LINE__,__FILE__)
          ! compute nHi
          rnum = (dd/rH0)**pH
          rden = (dd/rH0)**(pH+qH)
          fnum = 1-rnum
          fden = 1-rden
          IF (ABS(fden) .LT. 1.e-10_real_8) fden = 1.e-10_real_8
          nHi= nhi + fnum/fden

          ! compute der_nH_dr  index dependent
          ff = -(rnum*REAL(pH,kind=real_8)&
               -fnum*rden*REAL(pH+qH,kind=real_8)/fden)/&
               (fden*dd*dd)
          der_nH_dr(1,ii) = ff*dx
          der_nH_dr(2,ii) = ff*dy
          der_nH_dr(3,ii) = ff*dz

       ENDDO   ! Loop on the hydrogens

       ! jpark  compute r_k
       dx2 = x0-xion
       dy2 = y0-yion
       dz2 = z0-zion
       IF (.NOT. isos1%tisos) CALL pbc(dx2,dy2,dz2,dx2,dy2,dz2,1,parm%apbc,parm%ibrav)
       r_k = SQRT(dx2*dx2+dy2*dy2+dz2*dz2)


       ! jpark compute der_r_k_dr with respect to Ok
       der_r_k_dr(1) = dx2/r_k
       der_r_k_dr(2) = dy2/r_k
       der_r_k_dr(3) = dz2/r_k

       ii = 0
       ! Second Loop over the oxygens
       DO iiat1 = iini1,iend1
          ii = ii + 1
          dx=tscr(1,iiat1)-x0
          dy=tscr(2,iiat1)-y0
          dz=tscr(3,iiat1)-z0
          IF (.NOT. isos1%tisos) CALL pbc(dx,dy,dz,dx,dy,dz,1,parm%apbc,parm%ibrav)
          dd=SQRT(dx*dx+dy*dy+dz*dz)
          IF (dd .LT. 1.e-3_real_8) GOTO 10
          ! compute nO
          rnum = (dd/rO0)**pO
          rden = (dd/rO0)**(pO+qO)
          fnum = 1-rnum
          fden = 1-rden
          IF (ABS(fden) .LT. 1.e-10_real_8) fden = 1.e-10_real_8
          nOi= noi + fnum/fden

          ! compute der_nO_dr  
          ff = -(rnum*REAL(pO,kind=real_8)&
               -fnum*rden*REAL(pO+qO,kind=real_8)/fden)/&
               (fden*dd*dd)
          der_nO_dr(1,ii) = ff*dx
          der_nO_dr(2,ii) = ff*dy
          der_nO_dr(3,ii) = ff*dz
10        CONTINUE
       ENDDO ! Second Loop over the oxygens

       ! dbg 
       ! write(6,'(A,I3,3(A,f10.5))') ' ox ' , IA, ' nh ',nH, 
       ! c                               ' nO ', nO, '   r_k ',  r_k
       ! dbg


       ! exponential factor exp{L*nH}

       exp_fact = EXP(lambda*nHi*nOi)

       ! jpark  for summation of expnential factor
       sum_exp_fact   = sum_exp_fact + exp_fact
       sum_r_exp_fact = sum_r_exp_fact + r_k*exp_fact

       ! dbg 
       ! write(6,'(2(A,1E15.6))') ' r exp ' , r_k*exp_fact , 
       ! c                  ' exp ',exp_fact
       ! dbg

       ! Derivatives:
       der_fact = lambda*r_k*exp_fact*nHi
       der_fact2 = lambda*exp_fact*nHi
       ii = 0
       ! Loop over the hydrogens for the derivatives
       DO iiat1=iini1,iend1
          ii = ii + 1
          IF (lsk(1,iiat1).NE.0) THEN
             k=lsk(1,iiat1)
             an(k)=an(k)+ der_fact*der_nO_dr(1,ii)
             anf1(k)=anf1(k)+ der_fact2*der_nO_dr(1,ii)
          ENDIF
          IF (lsk(2,iiat1).NE.0) THEN
             k=lsk(2,iiat1)
             an(k)=an(k)+ der_fact*der_nO_dr(2,ii)
             anf1(k)=anf1(k)+ der_fact2*der_nO_dr(2,ii)
          ENDIF
          IF (lsk(3,iiat1).NE.0) THEN
             k=lsk(3,iiat1)
             an(k)=an(k)+ der_fact*der_nO_dr(3,ii)
             anf1(k)=anf1(k)+ der_fact2*der_nO_dr(3,ii)
          ENDIF
          IF (l1.NE.0) an(l1)=an(l1)-der_fact*der_nO_dr(1,ii)
          IF (l2.NE.0) an(l2)=an(l2)-der_fact*der_nO_dr(2,ii)
          IF (l3.NE.0) an(l3)=an(l3)-der_fact*der_nO_dr(3,ii)
          IF (l1.NE.0) anf1(l1)=anf1(l1)-der_fact2*der_nO_dr(1,ii)
          IF (l2.NE.0) anf1(l2)=anf1(l2)-der_fact2*der_nO_dr(2,ii)
          IF (l3.NE.0) anf1(l3)=anf1(l3)-der_fact2*der_nO_dr(3,ii)
       ENDDO


       ! Derivatives:
       der_fact = lambda*r_k*exp_fact*nOi
       der_fact2 = lambda*exp_fact*nOi
       ii = 0
       ! Loop over the hydrogens for the derivatives
       DO iiat2=iini2,iend2
          ii = ii + 1
          IF (lsk(1,iiat2).NE.0) THEN
             k=lsk(1,iiat2)
             an(k)=an(k)+ der_fact*der_nH_dr(1,ii)
             anf1(k)=anf1(k)+ der_fact2*der_nH_dr(1,ii)
          ENDIF
          IF (lsk(2,iiat2).NE.0) THEN
             k=lsk(2,iiat2)
             an(k)=an(k)+ der_fact*der_nH_dr(2,ii)
             anf1(k)=anf1(k)+ der_fact2*der_nH_dr(2,ii)
          ENDIF
          IF (lsk(3,iiat2).NE.0) THEN
             k=lsk(3,iiat2)
             an(k)=an(k)+ der_fact*der_nH_dr(3,ii)
             anf1(k)=anf1(k)+ der_fact2*der_nH_dr(3,ii)
          ENDIF
          IF (l1.NE.0) an(l1)=an(l1)-der_fact*der_nH_dr(1,ii)
          IF (l2.NE.0) an(l2)=an(l2)-der_fact*der_nH_dr(2,ii)
          IF (l3.NE.0) an(l3)=an(l3)-der_fact*der_nH_dr(3,ii)
          IF (l1.NE.0) anf1(l1)=anf1(l1)-der_fact2*der_nH_dr(1,ii)
          IF (l2.NE.0) anf1(l2)=anf1(l2)-der_fact2*der_nH_dr(2,ii)
          IF (l3.NE.0) anf1(l3)=anf1(l3)-der_fact2*der_nH_dr(3,ii)
       ENDDO

       ! compute term1 der_r_k_dr * exp_fact
       IF (l1.NE.0) an(l1)=an(l1)+der_r_k_dr(1)*exp_fact
       IF (l2.NE.0) an(l2)=an(l2)+der_r_k_dr(2)*exp_fact
       IF (l3.NE.0) an(l3)=an(l3)+der_r_k_dr(3)*exp_fact
       IF (l11.NE.0) an(l11)=an(l11)-der_r_k_dr(1)*exp_fact
       IF (l22.NE.0) an(l22)=an(l22)-der_r_k_dr(2)*exp_fact
       IF (l33.NE.0) an(l33)=an(l33)-der_r_k_dr(3)*exp_fact

100    CONTINUE
       ! write(6,*) 'O_k r_k nH e(nH)', IA, r_k, nH, exp_fact  
    ENDDO    ! Loop over the oxygens

    cval = sum_r_exp_fact/sum_exp_fact

    ! write(6,'((A,f15.6))') ' final cv ' , cval
    ! stop
    ! dbg

    ! compute term2 r_k*ANF0
    !$omp parallel do private(II) shared(CVAL)
    DO ii = 1,cotc0%nodim
       anf1(ii) = anf1(ii)*cval
    ENDDO

    ! write(6,*) 'HY_VAL ',CVAL
    der_fact = 1.0_real_8/sum_exp_fact
    ! Final factor for the derivatives
    !$omp parallel do private(IAT)
    DO iat = 1,cotc0%nodim
       an(iat) = (an(iat) - anf1(iat))*der_fact*fc_nL
    ENDDO

    ! Second loop on the list numL
    DO ial = 1,numL

       iatl = iatom(ial)

       l1 = lsk(1,iatl)
       l2 = lsk(2,iatl)
       l3 = lsk(3,iatl)

       ii = 0
       ! Loop over the hydrogens for the derivatives
       DO iat2=iini2,iend2
          ii = ii + 1
          IF (lsk(1,iat2).NE.0) THEN
             k=lsk(1,iat2)
             an(k)=an(k)+ der_fc_nL*der_nHl_dr(1,ii,ial)*cval
          ENDIF
          IF (lsk(2,iat2).NE.0) THEN
             k=lsk(2,iat2)
             an(k)=an(k)+ der_fc_nL*der_nHl_dr(2,ii,ial)*cval
          ENDIF
          IF (lsk(3,iat2).NE.0) THEN
             k=lsk(3,iat2)
             an(k)=an(k)+ der_fc_nL*der_nHl_dr(3,ii,ial)*cval
          ENDIF
          IF (l1.NE.0) an(l1)=an(l1)-der_fc_nL*der_nHl_dr(1,ii,ial)*cval
          IF (l2.NE.0) an(l2)=an(l2)-der_fc_nL*der_nHl_dr(2,ii,ial)*cval
          IF (l3.NE.0) an(l3)=an(l3)-der_fc_nL*der_nHl_dr(3,ii,ial)*cval
       ENDDO

    ENDDO  ! IAL

    cval = cval*fc_nL


    DEALLOCATE(der_nH_dr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(der_nO_dr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(der_nHl_dr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(anf1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==---------------------------------------------------------==
    RETURN
  END SUBROUTINE cut_hyd_ion_distance
  ! =============================================================
  ! ==--------------------------------------------------------------==
  SUBROUTINE cdipolerho(ia,tscr,fci,fvi,an,ib,iatom,qatom,cval)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ia
    REAL(real_8)                             :: tscr(3,*), fci, fvi, an(*)
    INTEGER                                  :: ib, iatom(*)
    REAL(real_8)                             :: qatom(*), cval

    INTEGER                                  :: i
    REAL(real_8)                             :: dx, dy, dz, ff, x0, xx(3), &
                                                y0, z0

! Variables
! ==--------------------------------------------------------------=

    fvi=0._real_8
    ff=0._real_8
    dx=0._real_8
    dy=0._real_8
    dz=0._real_8

    CALL fillc(ia,tscr,xx)
    x0 = xx(1)
    y0 = xx(2)
    z0 = xx(3)

    DO i=1,ib
       CALL fillc(iatom(i),tscr,xx)
       xx(1)=xx(1)-x0
       xx(2)=xx(2)-y0
       xx(3)=xx(3)-z0
       IF (.NOT. isos1%tisos) CALL pbc(xx(1),xx(2),xx(3),xx(1),xx(2),xx(3),&
            1,parm%apbc,parm%ibrav)
       dx=dx+xx(1)*qatom(i)
       dy=dy+xx(2)*qatom(i)
       dz=dz+xx(3)*qatom(i)
       ff=ff+qatom(i)
    ENDDO
    dx=dx/ff
    dy=dy/ff
    dz=dz/ff
    fvi=SQRT(dx*dx+dy*dy+dz*dz)
    fci=fvi-cval

    DO i=1,ib
       an(3*(iatom(i)-1)+1)= dx/fvi*qatom(i)/ff
       an(3*(iatom(i)-1)+2)= dy/fvi*qatom(i)/ff
       an(3*(iatom(i)-1)+3)= dz/fvi*qatom(i)/ff
       an(3*(ia-1)+1)= an(3*(ia-1)+1) - an(3*(iatom(i)-1)+1)
       an(3*(ia-1)+2)= an(3*(ia-1)+2) - an(3*(iatom(i)-1)+2)
       an(3*(ia-1)+3)= an(3*(ia-1)+3) - an(3*(iatom(i)-1)+3)
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE cdipolerho
  ! ==================================================================
  SUBROUTINE cdipolephi(ia,tscr,fci,fvi,an,ib,iatom,qatom,cval)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ia
    REAL(real_8)                             :: tscr(3,*), fci, fvi, an(*)
    INTEGER                                  :: ib, iatom(*)
    REAL(real_8)                             :: qatom(*), cval

    INTEGER                                  :: i
    REAL(real_8)                             :: dx, dy, dz, fact, ff, fff, &
                                                x0, xx(3), y0, z0

! Variables
! ==--------------------------------------------------------------=

    fvi=0._real_8
    ff=0._real_8
    dx=0._real_8
    dy=0._real_8
    dz=0._real_8

    CALL fillc(ia,tscr,xx)
    x0 = xx(1)
    y0 = xx(2)
    z0 = xx(3)

    DO i=1,ib
       CALL fillc(iatom(i),tscr,xx)
       xx(1)= xx(1)-x0
       xx(2)= xx(2)-y0
       xx(3)= xx(3)-z0
       IF (.NOT. isos1%tisos) CALL pbc(xx(1),xx(2),xx(3),xx(1),xx(2),xx(3),&
            1,parm%apbc,parm%ibrav)
       dx=dx + xx(1)*qatom(i)
       dy=dy + xx(2)*qatom(i)
       dz=dz + xx(3)*qatom(i)
       ff=ff+qatom(i)
    ENDDO

    dx=dx/ff
    dy=dy/ff
    dz=dz/ff
    fff=SQRT(dx*dx+dy*dy+dz*dz)
    fvi=dacos(dz/fff)
    fact=SQRT(dx*dx+dy*dy)
    fci=fvi-cval

    DO i=1,ib
       an(3*(iatom(i)-1)+1)=dx*dz/(fff*fff*fact)*qatom(i)/ff
       an(3*(iatom(i)-1)+2)=dy*dz/(fff*fff*fact)*qatom(i)/ff
       an(3*(iatom(i)-1)+3)=dz**2/(fff*fff*fact)*qatom(i)/ff&
            - 1/fact*qatom(i)/ff
       an(3*(ia-1)+1)= an(3*(ia-1)+1) - an(3*(iatom(i)-1)+1)
       an(3*(ia-1)+2)= an(3*(ia-1)+2) - an(3*(iatom(i)-1)+2)
       an(3*(ia-1)+3)= an(3*(ia-1)+3) - an(3*(iatom(i)-1)+3)
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE cdipolephi
  ! ==================================================================
  SUBROUTINE cdipoletheta(ia,tscr,fci,fvi,an,ib,iatom,qatom,cval)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ia
    REAL(real_8)                             :: tscr(3,*), fci, fvi, an(*)
    INTEGER                                  :: ib, iatom(*)
    REAL(real_8)                             :: qatom(*), cval

    INTEGER                                  :: i
    REAL(real_8)                             :: dx, dy, dz, fact, ff, pi, x0, &
                                                xx(3), y0, z0

! Variables
! ==--------------------------------------------------------------=

    fvi=0._real_8
    ff=0._real_8
    dx=0._real_8
    dy=0._real_8
    dz=0._real_8
    pi=dacos(-1.0_real_8)

    CALL fillc(ia,tscr,xx)
    x0 = xx(1)
    y0 = xx(2)
    z0 = xx(3)

    DO i=1,ib
       CALL fillc(iatom(i),tscr,xx)
       xx(1) = xx(1)-x0
       xx(2) = xx(2)-y0
       xx(3) = xx(3)-z0
       IF (.NOT. isos1%tisos) CALL pbc(xx(1),xx(2),xx(3),xx(1),xx(2),xx(3),&
            1,parm%apbc,parm%ibrav)
       dx = dx + xx(1)*qatom(i)
       dy = dy + xx(2)*qatom(i)
       dz = dz + xx(3)*qatom(i)
       ff = ff + qatom(i)
    ENDDO

    dx=dx/ff
    dy=dy/ff
    dz=dz/ff
    fvi=dy/dx

    IF (dx.LT.0.0_real_8) THEN
       fvi=ATAN(fvi)+pi
    ELSEIF (dx.GT.0.0_real_8.AND.dy.LT.0.0_real_8) THEN
       fvi=2.0_real_8*pi+ATAN(fvi)
    ELSE
       fvi=ATAN(fvi)
    ENDIF

    fact=SQRT(dx*dx+dy*dy)
    fci=fvi-cval

    DO i=1,ib
       an(3*(iatom(i)-1)+1) = -dy/fact*fact*qatom(i)/ff
       an(3*(iatom(i)-1)+2) =  dx/fact*fact*qatom(i)/ff
       an(3*(iatom(i)-1)+3) = 0.0_real_8
       an(3*(ia-1)+1)= an(3*(ia-1)+1) - an(3*(iatom(i)-1)+1)
       an(3*(ia-1)+2)= an(3*(ia-1)+2) - an(3*(iatom(i)-1)+2)
       an(3*(ia-1)+3)= an(3*(ia-1)+3) - an(3*(iatom(i)-1)+3)
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE cdipoletheta
  ! ==================================================================
  SUBROUTINE vcoors(isp,cvac,c1_kappa,c1_rc,tscr,an,lsk,&
       fci,fvi,cval)
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: isp
    REAL(real_8)                             :: cvac(*), c1_kappa, c1_rc, &
                                                tscr(3,*), an(*)
    INTEGER                                  :: lsk(3,ions1%nat)
    REAL(real_8)                             :: fci, fvi, cval

    INTEGER                                  :: iat, iend, iini, is, k
    REAL(real_8)                             :: dd, df, dx, dy, dz, ff, x0, &
                                                y0, z0

    x0 = cvac(1)
    y0 = cvac(2)
    z0 = cvac(3)

    fvi = 0._real_8
    iini = 1
    DO is = 1,isp-1
       iini=iini+ions0%na(is)
    ENDDO
    iend = iini + ions0%na(isp) -1
    DO iat=iini,iend
       dx=tscr(1,iat)-x0
       dy=tscr(2,iat)-y0
       dz=tscr(3,iat)-z0
       IF (.NOT. isos1%tisos) CALL pbc(dx,dy,dz,dx,dy,dz,1,parm%apbc,parm%ibrav)
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
       ENDIF
    ENDDO
    fci=fvi-cval
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE vcoors
  ! ==================================================================


END MODULE meta_cv_utils
