MODULE detdof_utils
  USE adat,                            ONLY: elem
  USE cnstfc_utils,                    ONLY: cnstfc,&
                                             restfc
  USE cnstpr_utils,                    ONLY: cnstpr
  USE comvelmod,                       ONLY: comvl
  USE cotr,                            ONLY: &
       anorm, askel, bsigma, cnsval, cotc0, cotr007, csigm, dsigma, dtm, &
       duat, fc, fcstr, fv, lskcor, lskptr, mm_askel, ntcnst, ntrest, patot, &
       pmall, resc, resfdiff, resm, resv, resval, rskel, rsmass, tsigma, &
       xlagr, ylagr
  USE dum2_utils,                      ONLY: dum2,&
                                             dumpr
  USE error_handling,                  ONLY: stopgm
  USE glemod,                          ONLY: tglepc
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE isos,                            ONLY: isos1
  USE kinds,                           ONLY: real_8
  USE mm_dimmod,                       ONLY: nat_grm,&
                                             solsolv,&
                                             solvvv
  USE mm_input,                        ONLY: lqmmm
  USE nose,                            ONLY: glib,&
                                             loct,&
                                             nosl,&
                                             tcafes,&
                                             tnosepc
  USE parac,                           ONLY: paral
  USE pimd,                            ONLY: pimd1
  USE rmas,                            ONLY: rmass
  USE system,                          ONLY: cntl,&
                                             iatpt,&
                                             maxsys
  USE tpar,                            ONLY: dt2bym
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: detdof
  !public :: hackcopy
  !public :: builda
  !public :: buildr
  PUBLIC :: qmdof

CONTAINS

  ! ==================================================================
  SUBROUTINE detdof(tau0,tscr)
    ! ==--------------------------------------------------------------==
    ! == DETermination of Degrees Of Freedom                          ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:), tscr(:,:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'detdof'

    CHARACTER(len=5), DIMENSION(0:1)         :: infout = (/'  FIX','  VAR'/)
    INTEGER                                  :: i, ia, iat, ib, ic, id, ierr, &
                                                ifix, iif, is, ityp, k, l, &
                                                l1, l2, l3, nfix
    INTEGER, SAVE                            :: ifirst = 0
    REAL(real_8)                             :: cval, dummy
    REAL(real_8), ALLOCATABLE                :: dxpar(:,:,:)

! ==--------------------------------------------------------------==
! Only call this routine once

    IF (ifirst.NE.0) RETURN
    ifirst=1
    ! ==--------------------------------------------------------------==
    ! ==  GET TOTAL NUMBER OF PARAMETERS TO BE OPTIMIZED              ==
    ! ==--------------------------------------------------------------==
    cotc0%nodim = 0
    is = 0
    DO i=1,ions1%nat
       DO k=1,3
          IF (lskcor(k,i).NE.0) THEN
             cotc0%nodim=cotc0%nodim+1
             lskptr(k,i)=cotc0%nodim
          ELSE
             IF (i.GT.solsolv%nrpt) is=is+1
             lskptr(k,i)=0
          ENDIF
       ENDDO
    ENDDO
    IF (cotc0%nodim.NE.3*ions1%nat) THEN
       glib=REAL(cotc0%nodim,kind=real_8)
    ELSEIF (nosl%tultra.OR.nosl%tmnose.OR.tcafes.OR.loct%tloct) THEN
       glib=REAL(cotc0%nodim,kind=real_8)
    ELSE
       glib=REAL(cotc0%nodim,kind=real_8)-3._real_8
       IF (isos1%tisos.AND..NOT.lqmmm%qmmm) THEN
          IF (ions1%nat.EQ.1) THEN
             glib=glib
          ELSEIF (ions1%nat.EQ.2) THEN
             glib=glib-2._real_8
          ELSE
             glib=glib-3._real_8
          ENDIF
       ENDIF
    ENDIF
    glib=glib-cotc0%mcnstr
    IF (lqmmm%qmmm) glib=glib - INT(REAL(3*solsolv%nsolv-is,kind=real_8)/3.0_real_8&
         /REAL(solvvv%nram_gr,kind=real_8)*REAL(solvvv%ncons_gr,kind=real_8))
    ! T.D. 7/10/1998. CHANGED TEST TO cntl%md. AK, 2005/07/31.
    IF (cntl%md) THEN
       IF (glib.LT.1._real_8) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A)')&
               '  ! ! WARNING FROM DETDOF: WHY IS GLIB LESS THAN UNITY !!'
          glib=1._real_8
          IF (paral%io_parent)&
               WRITE(6,'(A)')&
               '  ! ! WARNING FROM DETDOF: NOW IT IS SET EQUAL TO UNITY !!'
       ENDIF
       IF (paral%io_parent)&
            WRITE(6,'(/,A,T58,I8)')&
            ' DEGREES OF FREEDOM FOR SYSTEM:',NINT(glib)
    ELSE
       ! AK   force DOF to be >=1 for all non-cntl%md job types, but don't fuss about it.
       IF (glib.LT.1._real_8) THEN
          glib=1._real_8
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ==  STEEPEST DESCENT SCALING PARAMETERS                         ==
    ! ==--------------------------------------------------------------==
    ALLOCATE(dtm(cotc0%nodim),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    iat=0
    l=0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          iat=iat+1
          DO k=1,3
             IF (lskcor(k,iat).NE.0) THEN
                l=l+1
                dtm(l)=dt2bym(is)*lskcor(k,iat)
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    ! ==  MASSES                                                      ==
    ! ==--------------------------------------------------------------==
    ALLOCATE(pmall(maxsys%nax*maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    patot=0.0_real_8
    iat=0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          iat=iat+1
          ifix=lskcor(1,iat)+lskcor(1,iat)+lskcor(1,iat)
          IF (ifix.NE.3) THEN
             pmall(iat)=0.0_real_8
          ELSE
             pmall(iat)=rmass%pma0(is)
          ENDIF
          patot=patot+pmall(iat)
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    ! ==  CONSTRAINTS                                                 ==
    ! ==--------------------------------------------------------------==
    CALL dum2(tau0,tscr)
    IF (cotc0%mcnstr.GT.0 .OR. cotr007%mrestr.GT.0) THEN
       l=MAX(cotc0%mcnstr,cotr007%mrestr)
       ALLOCATE(anorm(cotc0%nodim,l),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(fcstr(cotc0%nodim),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    IF (cotc0%mcnstr.GT.0) THEN
       ALLOCATE(askel(cotc0%nodim,cotc0%mcnstr,12),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)

       ALLOCATE(mm_askel(cotc0%mcnstr,12),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)

       ALLOCATE(fc(cotc0%nodim),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(fv(cotc0%nodim),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(csigm(cotc0%nodim),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(xlagr(cotc0%nodim),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(ylagr(cotc0%nodim),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(askel)!,12*cotc0%nodim*cotc0%mcnstr)

       CALL zeroing(mm_askel)!,12*cotc0%mcnstr)

       CALL zeroing(xlagr)!,cotc0%mcnstr)
       CALL zeroing(ylagr)!,cotc0%mcnstr)
       DO i=1,cotc0%mcnstr
          ityp=ntcnst(1,i)
          IF (ityp.EQ.1) THEN
             csigm(i)=dsigma
          ELSEIF (ityp.EQ.2) THEN
             csigm(i)=bsigma
          ELSEIF (ityp.EQ.3) THEN
             csigm(i)=tsigma
          ELSEIF (ityp.EQ.4 .OR. ityp.EQ.10 .OR. ityp.EQ.7) THEN
             csigm(i)=dsigma
          ELSEIF(ityp.EQ.5 .OR. ityp.EQ.6 .OR. ityp.EQ.8&
               .OR. ityp.EQ.9) THEN
             csigm(i)=tsigma
          ENDIF
          IF (ityp .NE. 6 .AND. ityp .NE. 8 .AND.&
               ityp .NE. 9 .AND. ityp .NE. 10) THEN
             ia=ntcnst(2,i)
             ib=ntcnst(3,i)
             ic=ntcnst(4,i)
             id=ntcnst(5,i)
             CALL builda(ia,i,1)
             CALL builda(ib,i,2)
             CALL builda(ic,i,3)
             CALL builda(id,i,4)
          ENDIF
       ENDDO
    ENDIF
    IF (cotr007%mrestr.GT.0) THEN
       ALLOCATE(rskel(cotc0%nodim,cotr007%mrestr,12),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(rskel)!,12*cotc0%nodim*cotr007%mrestr)
       ALLOCATE(resc(cotc0%nodim),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(resv(cotc0%nodim),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(resm(cotr007%mrestr,cotr007%mrestr),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(resfdiff(cotr007%mrestr),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(rsmass(cotc0%nodim),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(rsmass)!,cotc0%nodim)
       DO i=1,cotr007%mrestr
          ia=ntrest(2,i)
          ib=ntrest(3,i)
          ic=ntrest(4,i)
          id=ntrest(5,i)
          CALL buildr(ia,i,1)
          CALL buildr(ib,i,2)
          CALL buildr(ic,i,3)
          CALL buildr(id,i,4)
       ENDDO
       iat=0
       l=0
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             iat=iat+1
             DO k=1,3
                IF (lskcor(k,iat).NE.0) THEN
                   l=l+1
                   ! RSMASS(L)=PMA0(IS)
                   CALL hackcopy(rmass%pma0(is),rsmass,l)
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ==  PRINT SOME INFO                                             ==
    ! ==--------------------------------------------------------------==
    IF (paral%parent) THEN
       CALL dumpr
    ENDIF
    ! fix COM does only work for (some) geometry optimizations. AK 2004-11-07    
    IF (cotc0%lfcom.AND.(cntl%lbfgs.OR.cntl%md)) THEN
       IF (paral%io_parent)&
            WRITE(6,'(A)') ' WARNING: FIX COM has no effect. see manual.'
    ENDIF
    iat=0
    nfix=0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          iat=iat+1
          l1=lskcor(1,iat)
          l2=lskcor(2,iat)
          l3=lskcor(3,iat)
          ifix=l1+l2+l3
          IF (ifix.NE.3) THEN
             nfix=nfix+1
          ENDIF
       ENDDO
    ENDDO
    IF (nfix.LT.500)THEN
       iif=0
       iat=0
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             iat=iat+1
             l1=lskcor(1,iat)
             l2=lskcor(2,iat)
             l3=lskcor(3,iat)
             ifix=l1+l2+l3
             IF (ifix.NE.3) THEN
                IF (iif.EQ.0.AND.paral%io_parent) WRITE(6,'(A,/,A)')&
                     ' FIXED COORDINATES ',&
                     '   NR   TYPE         X       Y       Z'
                iif=1
                IF (paral%io_parent)&
                     WRITE(6,'(I5,5X,A,5X,3(A,3X))') nat_grm(iat),&
                     elem%el(ions0%iatyp(is)),infout(l1),infout(l2),infout(l3)
             ENDIF
          ENDDO
       ENDDO
    ELSE
       IF (paral%io_parent)&
            WRITE(6,*)'THERE ARE QUITE A LOT OF FIXED ATOMS'
    ENDIF
    IF (cotc0%mcnstr.GT.0.OR.cotr007%mrestr.GT.0) THEN
       ALLOCATE(dxpar(3,maxsys%nax,maxsys%nsx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(dxpar)!,3*maxsys%nax*maxsys%nsx)
       IF (cotc0%mcnstr.GT.0) THEN
          CALL cnstfc(tau0,tscr,dxpar,dummy,.TRUE.)
          DO i=1,cotc0%mcnstr
             cval=cnsval(i)
             IF (cval.EQ.-999._real_8) THEN
                cnsval(i)=fv(i)
                fc(i)=0.0_real_8
             ENDIF
          ENDDO
       ENDIF
       CALL restfc(tau0,dxpar)
       IF (cotr007%mrestr.GT.0) THEN
          DO i=1,cotr007%mrestr
             cval=resval(i)
             IF (cval.EQ.-999._real_8) THEN
                resval(i)=resv(i)
             ENDIF
          ENDDO
       ENDIF
       CALL cnstpr
       DEALLOCATE(dxpar,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    IF (paral%io_parent)&
         WRITE(6,*)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE detdof
  SUBROUTINE hackcopy(a,b,i)
    REAL(real_8)                             :: a, b(*)
    INTEGER                                  :: i

    b(i)=a
    RETURN
  END SUBROUTINE hackcopy
  ! ==================================================================
  SUBROUTINE builda(ia,i,n)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ia, i, n

    INTEGER                                  :: ib, is, ityp, j, jj, jx, jy, &
                                                jz, kx, naa, nn
    REAL(real_8)                             :: ff

    IF (ia.EQ.0) RETURN
    nn=(n-1)*3
    IF (ia.LT.0) THEN
       IF (paral%io_parent)&
            WRITE(6,'(1X,A,I5,A,I5,A,I1)') 'BUILDA! IA=',ia,' I=',i,' N=',N
       CALL stopgm('BUILDA','IA IS LESS THAN ZERO',& 
            __LINE__,__FILE__)
    ELSEIF (ia.LE.ions1%nat) THEN
       jx=lskptr(1,ia)
       jy=lskptr(2,ia)
       jz=lskptr(3,ia)
       IF (jx.NE.0) askel(jx,i,nn+1)=1.0_real_8
       IF (jy.NE.0) askel(jy,i,nn+2)=1.0_real_8
       IF (jz.NE.0) askel(jz,i,nn+3)=1.0_real_8

       IF (jx.NE.0) mm_askel(i,nn+1)=jx
       IF (jy.NE.0) mm_askel(i,nn+2)=jy
       IF (jz.NE.0) mm_askel(i,nn+3)=jz
    ELSEIF (ia.LE.ions1%nat+duat%ndat) THEN
       naa=ia-ions1%nat
       ityp=duat%listda(naa,1)
       IF (ityp.EQ.1) THEN
          RETURN
       ELSEIF (ityp.EQ.2) THEN
          kx=duat%listda(naa,2)
          jj=duat%listd2(1,kx)
          IF (jj.LT.0) THEN
             jj=ions1%nat
          ELSE IF (jj.EQ.0) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(1X,A,I5,A,I5,A,I1)')&
                  'BUILDA! IA=',IA,' I=',I,' N=',N
             CALL stopgm('BUILDA','INCORRECT WEIGHTS FOR IA',& 
                  __LINE__,__FILE__)
          ENDIF
          ff=1.0_real_8/REAL(jj,kind=real_8)
          DO j=1,jj
             IF (jj.EQ.ions1%nat) THEN
                ib=j
             ELSE
                ib=duat%listd2(j+1,kx)
             ENDIF
             jx=lskptr(1,ib)
             jy=lskptr(2,ib)
             jz=lskptr(3,ib)
             IF (jx.NE.0) askel(jx,i,nn+1)=ff
             IF (jy.NE.0) askel(jy,i,nn+2)=ff
             IF (jz.NE.0) askel(jz,i,nn+3)=ff

             IF (jx.NE.0) mm_askel(i,nn+1)=jx
             IF (jy.NE.0) mm_askel(i,nn+2)=jy
             IF (jz.NE.0) mm_askel(i,nn+3)=jz
          ENDDO
       ELSEIF (ityp.EQ.3) THEN
          kx=duat%listda(naa,2)
          jj=duat%listd3(1,kx)
          IF (jj.LT.0) THEN
             jj=ions1%nat
          ELSE IF (jj.EQ.0) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(1X,A,I5,A,I5,A,I1)')&
                  'BUILDA! IA=',IA,' I=',I,' N=',N
             CALL stopgm('BUILDA','INCORRECT WEIGHTS FOR IA',& 
                  __LINE__,__FILE__)
          ENDIF
          ff=0.0_real_8
          DO j=1,jj
             IF (jj.EQ.ions1%nat) THEN
                ib=j
             ELSE
                ib=duat%listd3(j+1,kx)
             ENDIF
             is=iatpt(2,ib)
             ff=ff+rmass%pma0(is)
          ENDDO
          IF (ff.GT.0.0_real_8) THEN
             ff=1.0_real_8/ff
          ELSE
             IF (paral%io_parent)&
                  WRITE(6,'(1X,A,I5,A,I5,A,I1)')&
                  'BUILDA! IA=',IA,' I=',I,' N=',N
             CALL stopgm('BUILDA','INCORRECT WEIGHTS FOR IA',& 
                  __LINE__,__FILE__)
          ENDIF
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
             IF (jx.NE.0) askel(jx,i,nn+1)=rmass%pma0(is)*ff
             IF (jy.NE.0) askel(jy,i,nn+2)=rmass%pma0(is)*ff
             IF (jz.NE.0) askel(jz,i,nn+3)=rmass%pma0(is)*ff

             IF (jx.NE.0) mm_askel(i,nn+1)=jx
             IF (jy.NE.0) mm_askel(i,nn+2)=jy
             IF (jz.NE.0) mm_askel(i,nn+3)=jz
          ENDDO
       ELSEIF (ityp.EQ.4) THEN
          kx=duat%listda(naa,2)
          jj=duat%listd4(1,kx)
          ff=0.0_real_8
          DO j=1,jj
             ff=ff+duat%weigd4(j,kx)
          ENDDO
          IF (ff.GT.0.0_real_8) THEN
             ff=1.0_real_8/ff
          ELSE
             IF (paral%io_parent)&
                  WRITE(6,'(1X,A,I5,A,I5,A,I1)')&
                  'BUILDA! IA=',IA,' I=',I,' N=',N
             CALL stopgm('BUILDA','INCORRECT WEIGHTS FOR IA',& 
                  __LINE__,__FILE__)
          ENDIF
          DO j=1,jj
             ib=duat%listd4(j+1,kx)
             jx=lskptr(1,ib)
             jy=lskptr(2,ib)
             jz=lskptr(3,ib)
             IF (jx.NE.0) askel(jx,i,nn+1)=duat%weigd4(j,kx)*ff
             IF (jy.NE.0) askel(jy,i,nn+2)=duat%weigd4(j,kx)*ff
             IF (jz.NE.0) askel(jz,i,nn+3)=duat%weigd4(j,kx)*ff

             IF (jx.NE.0) mm_askel(i,nn+1)=jx
             IF (jy.NE.0) mm_askel(i,nn+2)=jy
             IF (jz.NE.0) mm_askel(i,nn+3)=jz
          ENDDO
       ENDIF
    ELSE
       IF (paral%io_parent)&
            WRITE(6,'(1X,A,I5,A,I5,A,I1)') 'BUILDA! IA=',ia,' I=',i,' N=',N
       CALL stopgm('BUILDA','INCORRECT VALUE FOR IA',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE builda
  ! ==================================================================
  SUBROUTINE buildr(ia,i,n)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ia, i, n

    INTEGER                                  :: ib, is, ityp, j, jj, jx, jy, &
                                                jz, kx, naa, nn
    REAL(real_8)                             :: ff

    IF (ia.EQ.0) RETURN
    nn=(n-1)*3
    IF (ia.LT.0) THEN
       IF (paral%io_parent)&
            WRITE(6,'(1X,A,I5,A,I5,A,I1)') 'BUILDR! IA=',ia,' I=',i,' N=',N
       CALL stopgm('BUILDR','IA IS LESS THAN ZERO',& 
            __LINE__,__FILE__)
    ELSEIF (ia.LE.ions1%nat) THEN
       jx=lskptr(1,ia)
       jy=lskptr(2,ia)
       jz=lskptr(3,ia)
       IF (jx.NE.0) rskel(jx,i,nn+1)=1.0_real_8
       IF (jy.NE.0) rskel(jy,i,nn+2)=1.0_real_8
       IF (jz.NE.0) rskel(jz,i,nn+3)=1.0_real_8
    ELSEIF (ia.LE.ions1%nat+duat%ndat) THEN
       naa=ia-ions1%nat
       ityp=duat%listda(naa,1)
       IF (ityp.EQ.1) THEN
          RETURN
       ELSEIF (ityp.EQ.2) THEN
          kx=duat%listda(naa,2)
          jj=duat%listd2(1,kx)
          IF (jj.LT.0) THEN
             jj=ions1%nat
          ELSE IF (jj.EQ.0) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(1X,A,I5,A,I5,A,I1)')&
                  'BUILDR! IA=',IA,' I=',I,' N=',N
             CALL stopgm('BUILDR','INCORRECT WEIGHTS FOR IA',& 
                  __LINE__,__FILE__)
          ENDIF
          ff=1.0_real_8/REAL(jj,kind=real_8)
          DO j=1,jj
             IF (jj.EQ.ions1%nat) THEN
                ib=j
             ELSE
                ib=duat%listd2(j+1,kx)
             ENDIF
             jx=lskptr(1,ib)
             jy=lskptr(2,ib)
             jz=lskptr(3,ib)
             IF (jx.NE.0) rskel(jx,i,nn+1)=ff
             IF (jy.NE.0) rskel(jy,i,nn+2)=ff
             IF (jz.NE.0) rskel(jz,i,nn+3)=ff
          ENDDO
       ELSEIF (ityp.EQ.3) THEN
          kx=duat%listda(naa,2)
          jj=duat%listd3(1,kx)
          IF (jj.LT.0) THEN
             jj=ions1%nat
          ELSE IF (jj.EQ.0) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(1X,A,I5,A,I5,A,I1)')&
                  'BUILDR! IA=',IA,' I=',I,' N=',N
             CALL stopgm('BUILDR','INCORRECT WEIGHTS FOR IA',& 
                  __LINE__,__FILE__)
          ENDIF
          ff=0.0_real_8
          DO j=1,jj
             IF (jj.EQ.ions1%nat) THEN
                ib=j
             ELSE
                ib=duat%listd3(j+1,kx)
             ENDIF
             is=iatpt(2,ib)
             ff=ff+rmass%pma0(is)
          ENDDO
          IF (ff.GT.0.0_real_8) THEN
             ff=1.0_real_8/ff
          ELSE
             IF (paral%io_parent)&
                  WRITE(6,'(1X,A,I5,A,I5,A,I1)')&
                  'BUILDR! IA=',IA,' I=',I,' N=',N
             CALL stopgm('BUILDR','INCORRECT WEIGHTS FOR IA',& 
                  __LINE__,__FILE__)
          ENDIF
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
             IF (jx.NE.0) rskel(jx,i,nn+1)=rmass%pma0(is)*ff
             IF (jy.NE.0) rskel(jy,i,nn+2)=rmass%pma0(is)*ff
             IF (jz.NE.0) rskel(jz,i,nn+3)=rmass%pma0(is)*ff
          ENDDO
       ELSEIF (ityp.EQ.4) THEN
          kx=duat%listda(naa,2)
          jj=duat%listd4(1,kx)
          ff=0.0_real_8
          DO j=1,jj
             ff=ff+duat%weigd4(j,kx)
          ENDDO
          IF (ff.GT.0.0_real_8) THEN
             ff=1.0_real_8/ff
          ELSE
             IF (paral%io_parent)&
                  WRITE(6,'(1X,A,I5,A,I5,A,I1)')&
                  'BUILDR! IA=',IA,' I=',I,' N=',N
             CALL stopgm('BUILDR','INCORRECT WEIGHTS FOR IA',& 
                  __LINE__,__FILE__)
          ENDIF
          DO j=1,jj
             ib=duat%listd4(j+1,kx)
             jx=lskptr(1,ib)
             jy=lskptr(2,ib)
             jz=lskptr(3,ib)
             IF (jx.NE.0) rskel(jx,i,nn+1)=duat%weigd4(j,kx)*ff
             IF (jy.NE.0) rskel(jy,i,nn+2)=duat%weigd4(j,kx)*ff
             IF (jz.NE.0) rskel(jz,i,nn+3)=duat%weigd4(j,kx)*ff
          ENDDO
       ENDIF
    ELSE
       IF (paral%io_parent)&
            WRITE(6,'(1X,A,I5,A,I5,A,I1)') 'BUILDR! IA=',ia,' I=',i,' N=',N
       CALL stopgm('BUILDR','INCORRECT VALUE FOR IA',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE buildr
  ! ==================================================================
  SUBROUTINE qmdof
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: i, k

! ==--------------------------------------------------------------==
! ..Do not change classical GLIB in test case

    IF (pimd1%testpi) RETURN
    ! ..NOTE: GLIB is the number of degrees of freedom for the first
    ! ..      (staging) bead, whereas all other beads have always 3*NAT
    ! ..      degrees of freedom each in the path integral case
    ! ..      (the latter case is handled via GLIB_S=3._real_8*NAT)
    ! ..      For CMD: TCENTRO=.TRUE. or RPMD: TRINGP=.TRUE.
    ! ..                       first bead should NOT be thermostatted 
    ! ..                       which means that GLIB has to be calculated
    ! ..                       as for a classical simulation
    ! ..      MCNSTR is the number of internal constraints
    ! ..             using FIX STRUCTURE
    ! ..      NODIM is the number of non-fixed cartesian coordinates
    ! ..            using FIX ATOM or FIX COORDINATES
    cotc0%nodim = 0
    DO i=1,ions1%nat
       DO k=1,3
          IF (lskcor(k,i).NE.0) THEN
             cotc0%nodim=cotc0%nodim+1
             lskptr(k,i)=cotc0%nodim
          ELSE
             lskptr(k,i)=0
          ENDIF
       ENDDO
    ENDDO
    IF (.NOT.((pimd1%tcentro.OR.pimd1%tringp)&
        .AND.(pimd1%tstage.OR.pimd1%tpinm)&
        .AND.(.NOT.tnosepc.OR..NOT.tglepc.OR.comvl%tsubcom))) THEN
       glib=REAL(cotc0%nodim,kind=real_8)
    ELSE
       IF (cotc0%nodim.NE.3*ions1%nat) THEN
          glib=REAL(cotc0%nodim,kind=real_8)
       ELSE
          glib=REAL(cotc0%nodim,kind=real_8)-3._real_8
          IF (isos1%tisos) THEN
             IF (ions1%nat.EQ.1) THEN
                glib=glib
             ELSEIF (ions1%nat.EQ.2) THEN
                glib=glib-2._real_8
             ELSE
                glib=glib-3._real_8
             ENDIF
          ENDIF
       ENDIF
    ENDIF
    glib=glib-cotc0%mcnstr
    IF (glib.LT.1.0_real_8) THEN
       IF (paral%io_parent)&
            WRITE(6,'(A)')&
            '  ! ! WARNING FROM DETDOF: WHY IS GLIB LESS THAN UNITY !!'
       glib=1._real_8
       IF (paral%io_parent)&
            WRITE(6,'(A)')&
            '  ! ! WARNING FROM QMDOF: NOW IT IS SET EQUAL TO UNITY !!'
    ENDIF
    IF (paral%io_parent)&
         WRITE(6,'(A,I8)')&
         '  QMDOF| NUMBER OF QUANTUM DEGREES OF FREEDOM FOR IP=1 :'&
         ,NINT(glib)
    IF (paral%io_parent)&
         WRITE(6,'(A,I8)')&
         '  QMDOF| NUMBER OF QUANTUM DEGREES OF FREEDOM FOR IP>1 :'&
         ,3*ions1%nat
    IF (.NOT.(pimd1%tpinm.OR.pimd1%tstage)) THEN
       IF (cotc0%mcnstr.NE.0.OR.cotc0%nodim.NE.3*ions1%nat) THEN
          IF (paral%io_parent)&
               WRITE(6,*) ' ANY TYPE OF GEOMETRY CONSTRAINTS ',&
               'ONLY ALLOWED WITH STAGING OR NORMAL MODE REPRESENTATION'
          CALL stopgm('QMDOF','CONSTRAINTS',& 
               __LINE__,__FILE__)
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE qmdof
  ! ==================================================================



END MODULE detdof_utils
