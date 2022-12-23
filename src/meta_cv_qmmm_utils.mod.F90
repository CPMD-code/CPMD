MODULE meta_cv_qmmm_utils
  USE cotr,                            ONLY: duat,&
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
  USE meta_cv_utils,                   ONLY: val_rmsd
  USE mm_dim_utils,                    ONLY: mm_dim
  USE mm_dimmod,                       ONLY: cpsp,&
                                             mm_go_mm,&
                                             mmdim,&
                                             nat_cpmd,&
                                             solsolv
  USE mm_input,                        ONLY: lqmmm
  USE parac,                           ONLY: paral
  USE pbc_utils,                       ONLY: pbc
  USE system,                          ONLY: parm
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: coorn_rf_q
  PUBLIC :: coorntot_rf_q
  PUBLIC :: coornseq_rf
  PUBLIC :: rmsd_seq_a

CONTAINS

  ! ==================================================================
  SUBROUTINE coorn_rf_q(iatin,c1_rc,r0_shift,tscr,an,&
       fci,fvi,cval,SPvsAT)
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: iatin(15)
    REAL(real_8)                             :: c1_rc, r0_shift, tscr(3,*), &
                                                an(*), fci, fvi, cval
    INTEGER                                  :: SPvsAT

    INTEGER :: ia, iat, ib, ielement, iend, iiat, iini, ityp, j, jj, jx, jy, &
      jz, k, kx, l1, l2, l3, mexp, my_nat, naa, nat_qm, nat_solu, nexp
    LOGICAL                                  :: status
    REAL(real_8)                             :: dd, dx, dy, dz, fden, ff, &
                                                fff, fnum, rden, rnum, x0, &
                                                xx(3), y0, z0

! Variables
! ==--------------------------------------------------------------==
! get total number of atoms and set indexing.

    IF (lqmmm%qmmm) THEN
       CALL mm_dim(mm_go_mm,status)
       my_nat=ions1%nat
       nat_solu = solsolv%nrpt
       nat_qm = mmdim%natq
    ELSE
       my_nat=mmdim%natm
       nat_solu = my_nat
       nat_qm = my_nat
    ENDIF
    ! dbg
    ! write(6,*) ' >>>> NATOMS', my_nat, nat_solu, nat_qm, SPvsAT
    ! dbg
    ! PARAMETERS
    ielement = iatin(2)
    nexp = iatin(3)
    mexp = iatin(4)

    ! The first index has been already converted into cpmd order
    ia = iatin(1)
    IF (ia.LE.my_nat) THEN
       l1 = lskptr(1,ia)
       l2 = lskptr(2,ia)
       l3 = lskptr(3,ia)
    ENDIF
    CALL fillc(ia,tscr,xx)
    x0 = xx(1)
    y0 = xx(2)
    z0 = xx(3)

    fvi = 0._real_8
    IF (SPvsAT .EQ. -1) THEN
       iini = 1
       iend = nat_solu
    ELSEIF (SPvsAT .EQ. -2) THEN
       iini = nat_solu + 1
       iend = my_nat
    ELSEIF (SPvsAT .EQ. -3) THEN
       iini = 1
       iend = my_nat
    ELSE
       iini = iatin(5)
       iend = iatin(6)
    ENDIF

    ! dbg
    ! write(6,*) ' >>>> IND', IA, IINI, IEND , IELEMENT
    ! dbg

    DO iiat=iini,iend
       ! Convert to CPMD order
       iat = NAT_cpmd(iiat)

       IF (SPvsAT .GE. -3 .OR. spvsat .LE. -6) THEN
          ! dbg
          ! write(6,*) '>>>> ITYP ', IAT , IATYP(cpsp(IIAT))
          ! dbg
          IF (ions0%iatyp(cpsp(iiat)).NE. ielement) GOTO 120
       ENDIF
       IF (SPvsAT .EQ. -4 .OR. spvsat .EQ. -6) THEN
          IF (iat .LE. nat_qm) GOTO 120
       ENDIF

       ! dbg
       ! write(6,*) ' >>>> LOOP', IIAT, IAT
       ! dbg

       dx=tscr(1,iat)-x0
       dy=tscr(2,iat)-y0
       dz=tscr(3,iat)-z0
       IF (.NOT. isos1%tisos) CALL pbc(dx,dy,dz,dx,dy,dz,1,parm%apbc,parm%ibrav)

       ! dbg
       ! write(6,*) ' >>>> TISOS',  TISOS, TCLUST, TONED, TTWOD
       ! dbg

       dd=SQRT(dx*dx+dy*dy+dz*dz)
       ! ..exclude self interaction
       IF (dd.GT.1.e-2_real_8) THEN
          rnum = ((dd-r0_shift)/c1_rc)**nexp
          rden = ((dd-r0_shift)/c1_rc)**(nexp+mexp)
          fnum = 1-rnum
          fden = 1-rden
          IF (ABS(fden) .LT.  1.e-10_real_8) THEN
             fvi = fvi + REAL(nexp,kind=real_8)/REAL(nexp+mexp,kind=real_8)
             ff = -0.5_real_8*REAL(nexp*mexp,kind=real_8)/REAL(nexp+mexp,kind=real_8)/(c1_rc*dd)
          ELSE
             fvi = fvi + fnum/fden
             ff = -(rnum*REAL(nexp,kind=real_8)-fnum*rden*REAL(nexp+mexp,kind=real_8)/&
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

          IF (ia.LE.my_nat) THEN
             IF (l1.NE.0) an(l1)=an(l1)-ff*dx
             IF (l2.NE.0) an(l2)=an(l2)-ff*dy
             IF (l3.NE.0) an(l3)=an(l3)-ff*dz
          ELSEIF (ia.LE.my_nat+duat%ndat) THEN
             naa=ia-my_nat
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
                   IF (jx.NE.0) an(jx)=an(jx)-ff*dx*fff
                   IF (jy.NE.0) an(jy)=an(jy)-ff*dy*fff
                   IF (jz.NE.0) an(jz)=an(jz)-ff*dz*fff
                ENDDO
             ENDIF
          ENDIF
110       CONTINUE
       ENDIF
120    CONTINUE
    ENDDO
    fci=fvi-cval
    ! dbg
    ! write(6,*) ' >>>> VAL', FVI
    ! STOP 'fin'
    ! dbg
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE coorn_rf_q
  ! ==================================================================
  SUBROUTINE coorntot_rf_q(iatin,c1_rc,r0_shift,tscr,an,&
       fci,fvi,cval,SPvsAT)
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: iatin(15)
    REAL(real_8)                             :: c1_rc, r0_shift, tscr(3,*), &
                                                an(*), fci, fvi, cval
    INTEGER                                  :: SPvsAT

    INTEGER :: iat1, iat2, ib, ielement, iend1, iend2, iiat2, iini1, iini2, &
      is, isp1, ityp, j, jj, jx, jy, jz, k, kx, l1, l2, l3, mexp, my_nat, &
      naa, nat_qm, nat_solu, nexp
    LOGICAL                                  :: status
    REAL(real_8)                             :: dd, dx, dy, dz, fden, ff, &
                                                fff, fnum, rden, rnum, una1, &
                                                x0, xx(3), y0, z0

! Variables
! ==--------------------------------------------------------------==
! get total number of atoms and set indexing.

    IF (lqmmm%qmmm) THEN
       CALL mm_dim(mm_go_mm,status)
       my_nat=ions1%nat
       nat_solu = solsolv%nrpt
       nat_qm = mmdim%natq
    ELSE
       my_nat=mmdim%natm
       nat_solu = my_nat
       nat_qm = my_nat
    ENDIF
    ! dbg
    ! write(6,*) ' >>>> NATOMS', my_nat, nat_solu, nat_qm, spvsat
    ! dbg
    ! PARAMETERS
    isp1 = iatin(1)
    ielement = iatin(2)
    nexp = iatin(3)
    mexp = iatin(4)

    ! This variable is only meaningful if the first species is within QM
    ! Species are defined only in cpmd order therefore the atom indexes
    ! here are already in the cpmd order and do not need conversion
    ! This inconsistency between first group and second group 
    ! is rather awful, with QMMM better to use COOR_SEQ
    iini1 = 1
    DO is = 1,isp1-1
       iini1=iini1+ions0%na(is)
    ENDDO
    iend1 = iini1 + ions0%na(isp1) -1
    una1 = 1.0_real_8/REAL(ions0%na(isp1),kind=real_8)
    IF (iend1 .GT. nat_qm) CALL stopgm('META_CV_QMMM',&
         'TOT_COOR first species should be QM',& 
         __LINE__,__FILE__)

    IF (SPvsAT .EQ. -1) THEN
       iini2 = 1
       iend2 = nat_solu
    ELSEIF (SPvsAT .EQ. -2) THEN
       iini2 = nat_solu + 1
       iend2 = my_nat
    ELSEIF (SPvsAT .EQ. -3) THEN
       iini2 = 1
       iend2 = my_nat
    ELSE
       iini2 = iatin(5)
       iend2 = iatin(6)
    ENDIF

    ! dbg
    ! write(6,*) ' >>>> IND', IINI1, IEND1, IINI2, IEND2
    ! dbg

    fvi = 0._real_8
    DO iat1 = iini1,iend1
       l1 = lskptr(1,iat1)
       l2 = lskptr(2,iat1)
       l3 = lskptr(3,iat1)
       x0 = tscr(1,iat1)
       y0 = tscr(2,iat1)
       z0 = tscr(3,iat1)

       DO iiat2=iini2,iend2
          ! Convert to CPMD order
          iat2 = NAT_cpmd(iiat2)

          IF (SPvsAT .GE. -3 .OR. spvsat .LE. -6) THEN
             IF (ions0%iatyp(cpsp(iiat2)).NE. ielement) GOTO 120
          ENDIF
          IF (SPvsAT .EQ. -4 .OR. spvsat .EQ. -6) THEN
             IF (iat2 .LE. nat_qm) GOTO 120
          ENDIF

          ! dbg
          ! write(6,*) ' >>>> LOOP',IAT1,IIAT2,IAT2,IATYP(cpsp(IIAT2))
          ! dbg

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

             IF (iat2.LE.my_nat) THEN
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
             ELSEIF (iat2.LE.my_nat+duat%ndat) THEN
                ! I am not sure that dummy atoms can work 
                ! any way only if SEQUENCE and without specifying the element
                naa=iat2-my_nat
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
120       CONTINUE
       ENDDO
    ENDDO
    fvi = fvi*una1
    fci=fvi-cval
    ! dbg
    ! write(6,*) ' >>>> VAL', FVI
    ! STOP 'fin'
    ! dbg
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE coorntot_rf_q
  ! ==================================================================
  SUBROUTINE coornseq_rf(iatin,c1_rc,r0_shift,tscr,an,&
       fci,fvi,cval,SPvsAT)
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: iatin(15)
    REAL(real_8)                             :: c1_rc, r0_shift, tscr(3,*), &
                                                an(*), fci, fvi, cval
    INTEGER                                  :: SPvsAT

    INTEGER :: iat1, iat2, icount1, ielement1, ielement2, iend1, iend2, &
      iiat1, iiat2, iini1, iini2, k, l1, l2, l3, mexp, my_nat, nat_qm, &
      nat_solu, nexp
    LOGICAL                                  :: status
    REAL(real_8)                             :: dd, dx, dy, dz, fden, ff, &
                                                fnum, rden, rnum, una1, x0, &
                                                xx(3), y0, z0

! Variables
! ==--------------------------------------------------------------==
! get total number of atoms and set indexing.

    IF (lqmmm%qmmm) THEN
       CALL mm_dim(mm_go_mm,status)
       my_nat=ions1%nat
       nat_solu = solsolv%nrpt
       nat_qm = mmdim%natq
    ELSE
       my_nat=mmdim%natm
       nat_solu = my_nat
       nat_qm = my_nat
    ENDIF
    ! dbg
    ! write(6,*) ' >>>> NATOMS', my_nat, nat_solu, nat_qm, SPvsAT
    ! dbg
    ! PARAMETERS
    ielement1 = iatin(1)
    ielement2 = iatin(2)
    nexp = iatin(3)
    mexp = iatin(4)

    ! 1  SOLU - SOLU
    ! 2  SOLU - SOLV
    ! 3   SEQ - SOLV
    ! 4   SEQ - SEQ
    ! 5   SEQ - SEQMM
    ! 6 SEQQM - SEQMM
    ! 7 SEQMM - SEQMM
    ! 
    IF (SPvsAT .LE. 2) THEN
       iini1 = 1
       iend1 = nat_solu
    ELSE
       iini1 = iatin(5)
       iend1 = iatin(6)
    ENDIF

    IF (SPvsAT .EQ. 1 .OR. spvsat .EQ. 3) THEN
       iini2 = nat_solu + 1
       iend2 = my_nat
    ELSEIF (SPvsAT .EQ. 2) THEN
       iini2 = 1
       iend2 = nat_solu
    ELSE
       iini2 = iatin(7)
       iend2 = iatin(8)
    ENDIF

    ! dbg
    ! write(6,*) ' >>>> IND', IINI1, IEND1, IINI2, IEND2

    ! write(6,*) ' >>>> ELE', IELEMENT1, IELEMENT2
    ! dbg

    fvi = 0._real_8

    ! Count the number of atoms in the first group
    icount1 = 0
    DO iiat1 = iini1,iend1
       ! Convert to CPMD order
       iat1 = NAT_cpmd(iiat1)
       ! Check the other conditions
       IF (ielement1 .NE. 0) THEN
          IF ( ions0%iatyp(cpsp(iiat1)) .NE. ielement1 ) GOTO 110
       ENDIF
       IF (SPvsAT .EQ. 6 .OR. spvsat .EQ. 8) THEN
          IF (iat1 .GT. nat_qm) GOTO 110
       ELSEIF (SPvsAT .EQ. 7) THEN
          IF (iat1 .LE. nat_qm) GOTO 110
       ENDIF
       icount1 = icount1 + 1
110    CONTINUE
    ENDDO ! IAT1
    una1 =  1.0_real_8/REAL(icount1,kind=real_8)

    ! dbg
    ! write(6,*) ' >>>> NAT1', ICOUNT1
    ! dbg
    ! stop 'nat' 
    DO iiat1 = iini1,iend1
       ! Convert to CPMD order
       iat1 = NAT_cpmd(iiat1)
       ! Check the other conditions
       IF (ielement1 .NE. 0) THEN
          IF ( ions0%iatyp(cpsp(iiat1)) .NE. ielement1 ) GOTO 120
       ENDIF
       IF (SPvsAT .EQ. 6 .OR. spvsat .EQ. 8) THEN
          IF (iat1 .GT. nat_qm) GOTO 120
       ELSEIF (SPvsAT .EQ. 7) THEN
          IF (iat1 .LE. nat_qm) GOTO 120
       ENDIF

       ! dbg
       ! write(6,*) ' >>>> LOOP1',IIAT1,IAT1,IATYP(cpsp(IIAT1))
       ! dbg

       l1 = lskptr(1,iat1)
       l2 = lskptr(2,iat1)
       l3 = lskptr(3,iat1)
       x0 = tscr(1,iat1)
       y0 = tscr(2,iat1)
       z0 = tscr(3,iat1)
       icount1 = icount1 + 1

       DO iiat2=iini2,iend2
          ! Convert to CPMD order
          iat2 = NAT_cpmd(iiat2)

          IF (ielement2 .NE. 0) THEN
             IF ( ions0%iatyp(cpsp(iiat2)) .NE. ielement2 ) GOTO 130
          ENDIF
          IF (SPvsAT .GE. 5 .AND. spvsat .LT. 8 ) THEN
             IF (iat2 .LE. nat_qm) GOTO 130
          ENDIF

          ! dbg
          ! write(6,*) ' >>>> LOOP2', IIAT2, IAT2,  IATYP(cpsp(IIAT2))
          ! dbg

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

             ! .. the dummy atoms cannot be used here
             IF (iat2.LE.my_nat) THEN
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
             ENDIF

             IF (l1.NE.0) an(l1)=an(l1)-ff*dx
             IF (l2.NE.0) an(l2)=an(l2)-ff*dy
             IF (l3.NE.0) an(l3)=an(l3)-ff*dz
          ENDIF
130       CONTINUE
       ENDDO                ! IIAT2
120    CONTINUE
    ENDDO  ! IIAT1

    fvi = fvi*una1
    fci=fvi-cval
    ! dbg
    ! write(6,*) '>>>> FVI ', FVI
    ! STOP 'VAL'
    ! dbg
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE coornseq_rf
  ! ==================================================================
  SUBROUTINE rmsd_seq_a(indexsp,maxrmsdat,tscr,an,&
       fci,fvi,cval,SPvsAT)
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: indexsp(15), maxrmsdat
    REAL(real_8)                             :: tscr(3,*), an(*), fci, fvi, &
                                                cval
    INTEGER                                  :: SPvsAT

    CHARACTER(*), PARAMETER                  :: procedureN = 'rmsd_seq_a'
    INTEGER, PARAMETER                       :: iunit = 58 

    CHARACTER(len=100)                       :: filen
    INTEGER :: i, iat, icount, ie, iele(15), iend, ierr, iiat, iini, l1, l2, &
      l3, my_nat, nat_qm, nat_rmsd, nat_solu, numelements
    INTEGER, ALLOCATABLE, SAVE               :: ipos_rmsd(:)
    INTEGER, SAVE                            :: ifirst = 0
    LOGICAL                                  :: ferror, fexist, status, takeit
    REAL(real_8)                             :: rmsd_val
    REAL(real_8), ALLOCATABLE                :: derr(:,:), r0_cm(:,:), &
                                                r0_ion(:,:), r_cm(:,:), &
                                                r_ion(:,:)
    REAL(real_8), ALLOCATABLE, SAVE          :: pos_0(:,:), rmsd_scr(:)

! ==--------------------------------------------------------------==
! get total number of atoms and set indexing.

    IF (lqmmm%qmmm) THEN
       CALL mm_dim(mm_go_mm,status)
       my_nat=ions1%nat
       nat_solu = solsolv%nrpt
       nat_qm = mmdim%natq
    ELSE
       my_nat=mmdim%natm
       nat_solu = my_nat
       nat_qm = my_nat
    ENDIF
    ! dbg
    ! write(6,*) ' >>>> NATOMS', my_nat, nat_solu, nat_qm, spvsat
    ! dbg
    ! PARAMETERS
    numelements = indexsp(1)
    nat_rmsd = indexsp(2)
    ! dbg
    ! write(6,*) '>>> MAXRMSDAT', MAXRMSDAT
    ! write(6,*) ' >>>> ELE', NUMELEMENTS, NAT_RMSD
    ! dbg
    CALL zeroing(iele)!,15/2+1)
    !$omp parallel do private(IE)
    DO ie = 1,numelements
       iele(ie) = indexsp(ie+4)
       ! dbg
       ! write(6,*) '>> IELE ', IELE(IE)
       ! dbg
    ENDDO

    IF (SPvsAT .EQ. 2 ) THEN
       iini = nat_solu + 1
       iend = my_nat
    ELSEIF (SPvsAT .EQ. 1) THEN
       iini = 1
       iend = nat_solu
    ELSE
       iini = indexsp(3)
       iend = indexsp(4)
    ENDIF
    ! dbg
    ! write(6,*) ' >>>> IND', IINI, IEND
    ! dbg

    IF (ifirst.EQ.0) THEN
       ! Allocate memory
       ALLOCATE(pos_0(3,my_nat),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ! TO reduce the memory requirement, allocate MAXRMSDAT
       ALLOCATE(ipos_rmsd(maxrmsdat),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(rmsd_scr(6*3*maxrmsdat),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)

       filen='p0cp_rmsd_seq.dat'
       IF (paral%io_parent)&
            INQUIRE(file='p0cp_rmsd_seq.dat',exist=fexist)
       IF (fexist) THEN
          IF (paral%io_parent)&
               CALL fileopen(iunit,filen,fo_old,ferror)
          DO iiat = 1,my_nat
             iat = NAT_cpmd(iiat)
             IF (paral%io_parent)&
                  READ(iunit,'(I8,3f20.10)') i,&
                  pos_0(1,iat), pos_0(2,iat), pos_0(3,iat)
          ENDDO
       ELSE
          IF (paral%io_parent)&
               CALL fileopen(iunit,filen,fo_new,ferror)
          DO iiat = 1,my_nat
             iat = NAT_cpmd(iiat)
             pos_0(1,iat) =  tscr(1,iat)
             pos_0(2,iat) =  tscr(2,iat)
             pos_0(3,iat) =  tscr(3,iat)
             IF (paral%io_parent)&
                  WRITE(iunit,'(I8,3f20.10)') iiat,&
                  pos_0(1,iat), pos_0(2,iat), pos_0(3,iat)
          ENDDO
       ENDIF
       IF (paral%io_parent)&
            CALL fileclose(iunit)
       ifirst=1
    ENDIF

    CALL zeroing(ipos_rmsd)!,nat_rmsd)
    CALL zeroing(rmsd_scr)!,6*3*nat_rmsd)

    ! TODO align for BG
    ALLOCATE(r_ion(3, nat_rmsd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(r0_ion(3, nat_rmsd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(derr(3, nat_rmsd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(r0_cm(3, nat_rmsd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(r_cm(3, nat_rmsd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    icount = 0
    DO iiat = iini, iend
       iat = NAT_cpmd(iiat)
       takeit = .FALSE.
       IF (numelements .NE. 0) THEN
          DO ie = 1,numelements
             IF (ions0%iatyp(cpsp(iiat)) .EQ. iele(ie))THEN
                takeit = .TRUE.
             ENDIF
          ENDDO
       ELSE
          takeit = .TRUE.
       ENDIF
       IF (takeit) THEN
          icount = icount + 1
          ipos_rmsd(icount) = iat
          r_ion(1,icount) = tscr(1,iat)
          r_ion(2,icount) = tscr(2,iat)
          r_ion(3,icount) = tscr(3,iat)
          r0_ion(1,icount) = pos_0(1,iat)
          r0_ion(2,icount) = pos_0(2,iat)
          r0_ion(3,icount) = pos_0(3,iat)
       ENDIF
    ENDDO ! IIAT
    ! dbg
    ! write (6,*) ICOUNT
    ! dbg

    IF (icount .GT. nat_rmsd .OR. icount .GT. maxrmsdat) THEN
       CALL stopgm('RMSD_SEQ','too many atoms counted',& 
            __LINE__,__FILE__)
    ENDIF
    ! DO I = 1,ICOUNT
    ! write(6,'(2I4,6f14.6)') I, 1, R_ION(1,I), R0_ION(1,I)
    ! write(6,'(2I4,6f14.6)') I, 2, R_ION(2,I), R0_ION(2,I)
    ! write(6,'(2I4,6f14.6)') I, 3, R_ION(3,I), R0_ION(3,I)
    ! enddo 

    ! Calculate RMSD_B and derivatives
    CALL val_rmsd(icount,r_ion,r0_ion,rmsd_val,derr,r0_cm,r_cm)
    ! dbg
    ! write(6,*) '>>>> RMSD', RMSD_VAL 
    ! dbg
    ! stop 'val'
    fvi = rmsd_val
    fci = fvi-cval

    DO iat = 1,icount
       l1=lskptr(1,iat)
       l2=lskptr(2,iat)
       l3=lskptr(3,iat)
       IF (l1 .NE. 0) an(l1) = an(l1) + derr(1,ipos_rmsd(iat))
       IF (l2 .NE. 0) an(l2) = an(l2) + derr(2,ipos_rmsd(iat))
       IF (l3 .NE. 0) an(l3) = an(l3) + derr(3,ipos_rmsd(iat))
    ENDDO  ! IAT
    ! ==--------------------------------------------------------------==      
    DEALLOCATE(r_ion,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(r0_ion,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(derr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(r0_cm,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(r_cm,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==      
    RETURN
  END SUBROUTINE rmsd_seq_a
  ! ================================================================== 

END MODULE meta_cv_qmmm_utils
