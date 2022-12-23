MODULE coninp_utils
  USE cnst,                            ONLY: pi
  USE cotr,                            ONLY: &
       bsigma, cnpar, cnsval, cnsval_dest, cotc0, cotr007, dsigma, grate, &
       gsrate, lskcor, ntcnst, ntrest, resfor, respar, respos, resval, &
       resval_dest, tsigma
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_old
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE mm_dim_utils,                    ONLY: mm_dim
  USE mm_dimmod,                       ONLY: cpsp,&
                                             mm_go_mm,&
                                             mm_revert,&
                                             mmdim,&
                                             nat_cpmd,&
                                             solsolv
  USE mm_input,                        ONLY: lqmmm
  USE parac,                           ONLY: paral
  USE readsr_utils,                    ONLY: readsi,&
                                             readsr
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: coninp
  PUBLIC :: raddeg

CONTAINS

  ! ==================================================================
  SUBROUTINE coninp(iunit)
    ! ==--------------------------------------------------------------==
    ! ==  Reads Constraint Input for Geometry                         ==
    ! ==--------------------------------------------------------------==
    ! ==                                                              ==
    ! ==  FIX {ALL,QM,MM,SOLUTE}                                      ==
    ! ==  FIX ATOM                                                    ==
    ! ==    nfix  n1  n2  ...                                         ==
    ! ==  FIX SEQUENCE                                                ==
    ! ==    nfirst nlast                                              ==
    ! ==  FIX ELEMENT [SEQUENCE]                                      ==
    ! ==    zv [ifirst ilast]                                         ==
    ! ==  FIX PPTYPE [SEQUENCE]                                       ==
    ! ==    isp [ifist ilast]                                         ==
    ! ==  FIX COORDINATES                                             ==
    ! ==    nfix                                                      ==
    ! ==    n1   ix   iy   iz                                         ==
    ! ==    ...                                                       ==
    ! ==  FIX COM                                                     ==
    ! ==  FIX STRUCTURE [SHOVE]                                       ==
    ! ==   nfix                                                       ==
    ! ==  DIST      n1  n2   R0  [+1,-1,0] GROWTH aa DEST Rd          ==
    ! ==  STRETCH   n1  n2   R0  [+1,-1,0] GROWTH aa DEST Rd          ==
    ! ==  DISAXIS   n1  n2  k   R0  [+1,-1,0] GROWTH aa DEST Rd       ==
    ! ==  BEND      n1  n2  n3 A0  [+1,-1,0] GROWTH aa DEST Ad        ==
    ! ==  TORSION   n1  n2  n3  n4  D0 [+1,-1,0] GROWTH aa DEST Dd    ==
    ! ==  OUTP      n1  n2  n3  n4  D0 [+1,-1,0] GROWTH aa DEST Dd    ==
    ! ==  DIFFER    n1  n2  n3  dR0 [+1,-1,0] GROWTH aa DEST dRd      ==
    ! ==  RIGID     na  n1  ...                                       ==
    ! ==            .......  n_{na} [+1,-1,0]                         ==
    ! ==  COORD     n1   K  RC  C0  [+1,-1,0] GROWTH aa DEST Cd       ==
    ! ==  COORSP    n1  is  K   RC  C0  [+1,-1,0] GROWTH aa DEST Cd   ==
    ! ==  COOR_RF   n1  is  N   M   RC  C0  [+1,-1,0] GROWTH aa DEST Cd=
    ! ==  BNSWT     n1  n2  N   M   RC  C0  [+1,-1,0]                  =
    ! ==  TOT_COOR  is1 is2 N   M   RC  C0  [+1,-1,0] GROWTH aa DEST Cd=
    ! ==  MMFIELD   n1   F0                                           ==
    ! ==  ...                                                         ==
    ! ==  RESTRAINTS                                                  ==
    ! ==   nres                                                       ==
    ! ==  DIST      n1  n2   R0  kval       GROWTH aa DEST Rd         ==
    ! ==  STRETCH   n1  n2   R0  kval       GROWTH aa DEST Rd         ==
    ! ==  BEND      n1  n2  n3 A0  kval     GROWTH aa DEST Ad         ==
    ! ==  TORSION   n1  n2  n3  n4  D0 kval GROWTH aa DEST Dd         ==
    ! ==  OUTP      n1  n2  n3  n4  D0 kval GROWTH aa DEST Dd         ==
    ! ==  COOR_RF   n1  is  N   M   RC  C0  kval GROWTH aa DEST Cd    ==
    ! ==  DIFFER    n1  n2  n3  dR0 kval    GROWTH aa DEST DdRd       ==
    ! ==  DISAXIS   n1  n2  k   R0  kval    GROWTH aa DEST Rd         ==
    ! ==  RESPOS    n1  x0  y0  z0  D0 kval (NO GROWTH ALLOWED)       ==
    ! ==  ...                                                         ==
    ! ==  PENALTY                                                     ==
    ! ==   dsigma bsigma tsigma                                       ==
    ! ==                                                              ==
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: iunit

    CHARACTER(*), PARAMETER                  :: procedureN = 'coninp'

    CHARACTER(len=80)                        :: line
    INTEGER                                  :: i, i1, idata(4), idummy, &
                                                ierr, ifirst, ig, ilast, &
                                                iout, my_nat, nafix, nfix, &
                                                nres, nrig
    INTEGER, ALLOCATABLE                     :: lfix(:), narig(:)
    LOGICAL                                  :: erread, ferror, status
    REAL(real_8)                             :: kplane

    dsigma=1._real_8
    bsigma=5._real_8
    tsigma=5._real_8
    cotc0%lfcom=.FALSE.
    cotc0%lshove=.FALSE.

    ! FIXME: AK 2005/05/25 TODO LIST:
    ! - make all atom index inputs run through CPMD_nat()
    ! - check all arguments for out-of-range errors.
    ! - fix and check documentation

    ! get total number of atoms and set indexing.
    IF (lqmmm%qmmm) THEN
       CALL mm_dim(mm_go_mm,status)
       my_nat=ions1%nat
    ELSE
       my_nat=mmdim%natm
    ENDIF
    ALLOCATE(lfix(my_nat+1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
10  CONTINUE
    IF (paral%io_parent)&
         READ(iunit,err=20,END=20,fmt='(A)') line
    IF (INDEX(line,'END').NE.0 .AND. INDEX(line,'CONSTR').NE.0)&
         GOTO 30
    ! FIX ALL and FIX QM MM are the same
    IF (INDEX(line,'FIX').NE.0 .AND. (INDEX(line,'ALL').NE.0 .OR.&
         (INDEX(line,'QM').NE.0.AND.INDEX(line,'MM').NE.0)) ) THEN
       DO i=1,my_nat
          lskcor(1,i)=0
          lskcor(2,i)=0
          lskcor(3,i)=0
       ENDDO
       GOTO 10
    ELSEIF (INDEX(line,'FIX').NE.0 .AND. INDEX(line,'QM').NE.0) THEN
       IF (lqmmm%qmmm) THEN
          idummy=mmdim%natq
       ELSE
          idummy=my_nat
       ENDIF
       DO i=1,idummy
          lskcor(1,i)=0
          lskcor(2,i)=0
          lskcor(3,i)=0
       ENDDO
       GOTO 10
    ELSEIF (INDEX(line,'FIX').NE.0 .AND. INDEX(line,'MM').NE.0) THEN
       IF (lqmmm%qmmm) THEN
          DO i=mmdim%natq+1,mmdim%natm
             lskcor(1,i)=0
             lskcor(2,i)=0
             lskcor(3,i)=0
          ENDDO
       ENDIF
       GOTO 10
    ELSEIF (INDEX(line,'FIX').NE.0 .AND. INDEX(line,'SOLU').NE.0) THEN
       IF (lqmmm%qmmm) THEN
          DO i=1,solsolv%nrpt
             lskcor(1,i)=0
             lskcor(2,i)=0
             lskcor(3,i)=0
          ENDDO
       ENDIF
       GOTO 10
    ELSEIF (INDEX(line,'FIX').NE.0.AND.INDEX(line,'ELEM').NE.0)THEN
       IF (INDEX(line,'SEQ').NE.0) THEN
          IF (paral%io_parent)&
               READ(iunit,err=20,fmt=*) idummy,ifirst,ilast
       ELSE
          IF (paral%io_parent)&
               READ(iunit,err=20,fmt=*) idummy
          ifirst=1
          ilast=my_nat
       ENDIF
       DO ig=ifirst,ilast
          i=NAT_cpmd(ig)
          IF (ions0%iatyp(cpsp(i)).EQ.idummy)THEN
             lskcor(1,i)=0
             lskcor(2,i)=0
             lskcor(3,i)=0
          ENDIF
       ENDDO
       GOTO 10
    ELSEIF (INDEX(line,'FIX').NE.0.AND.INDEX(line,'PPTY').NE.0)THEN
       IF (INDEX(line,'SEQ').NE.0) THEN
          IF (paral%io_parent)&
               READ(iunit,err=20,fmt=*) idummy,ifirst,ilast
       ELSE
          IF (paral%io_parent)&
               READ(iunit,err=20,fmt=*) idummy
          ifirst=1
          ilast=my_nat
       ENDIF
       DO ig=ifirst,ilast
          i=NAT_cpmd(ig)
          IF (cpsp(i).EQ.idummy)THEN
             lskcor(1,i)=0
             lskcor(2,i)=0
             lskcor(3,i)=0
          ENDIF
       ENDDO
       GOTO 10
    ELSEIF (INDEX(line,'FIX').NE.0.AND.INDEX(line,'SEQ').NE.0)THEN
       IF (paral%io_parent)&
            READ(iunit,err=20,fmt=*)ifirst,ilast
       DO ig=ifirst,ilast
          i=NAT_cpmd(ig)
          lskcor(1,i)=0
          lskcor(2,i)=0
          lskcor(3,i)=0
       ENDDO
       GOTO 10
    ELSE IF (INDEX(line,'FIX').NE.0 .AND. INDEX(line,'ATOM').NE.0) &
         THEN
       IF (paral%io_parent)&
            READ(iunit,err=20,END=20,fmt=*) idummy, (lfix(i),i=1,idummy)
       DO ig=1,idummy
          i=NAT_cpmd(lfix(ig))
          lskcor(1,i)=0
          lskcor(2,i)=0
          lskcor(3,i)=0
       ENDDO
       GOTO 10
    ELSEIF (INDEX(line,'FIX').NE.0 .AND. INDEX(line,'COOR').NE.0) THEN
       IF (paral%io_parent)&
            READ(iunit,err=20,END=20,fmt=*) idummy
       DO i=1,idummy
          IF (paral%io_parent)&
               READ(iunit,err=20,END=20,fmt=*) nafix,(idata(ig),ig=1,3)
          lskcor(1,NAT_cpmd(nafix))=idata(1)
          lskcor(2,NAT_cpmd(nafix))=idata(2)
          lskcor(3,NAT_cpmd(nafix))=idata(3)
       ENDDO
       GOTO 10
    ELSEIF (INDEX(line,'FIX').NE.0 .AND. INDEX(line,'COM').NE.0) THEN
       cotc0%lfcom=.TRUE.
       GOTO 10
    ELSEIF (INDEX(line,'FIX').NE.0 .AND. INDEX(line,'STRU').NE.0) THEN
       IF (INDEX(line,'SHOVE').NE.0) cotc0%lshove=.TRUE.
       IF (paral%io_parent) READ(iunit,err=20,END=20,fmt=*) nfix
       ALLOCATE(ntcnst(6,nfix),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(cnsval(nfix),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(cnpar(2,nfix),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(grate(nfix),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(cnsval_dest(nfix),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(cnpar)!,2*nfix)
       CALL zeroing(grate)!,nfix)
       CALL zeroing(ntcnst)!,6*nfix)
       !$omp parallel do private(I)
       DO i=1,nfix
          cnsval(i)=-999._real_8
          cnsval_dest(i)=-999._real_8
       ENDDO
       cotc0%mcnstr=0
       ! ==------------------------------------------------------------==
       ! Read each structure.
11     CONTINUE
       IF (paral%io_parent)&
            READ(iunit,err=20,END=20,fmt='(A)') line
       IF (INDEX(line,'END').NE.0 .AND. INDEX(line,'CONSTR').NE.0) THEN
          ! Problem the number of constraints MCNSTR .NE. given NFIX.
          IF (paral%io_parent)&
               WRITE(6,'(1X,A,I5)') 'CONINP! # OF GIVEN CONSTRAINTS=',nfiX
          IF (paral%io_parent)&
               WRITE(6,'(1X,A,I5)') 'CONINP! # OF FOUND CONSTRAINTS=',cotc0%mcnstr
          IF (paral%io_parent)&
               WRITE(6,'(1X,A)')    'CONINP! END CONSTRAINTS IS REACHED.'
          CALL stopgm('CONINP','ERROR WHILE READING FIX STRUCTURES',& 
               __LINE__,__FILE__)
          ! ==------------------- DIST ---------------------------------==
       ELSEIF (INDEX(line,'DIST').NE.0) THEN
          cotc0%mcnstr=cotc0%mcnstr+1
          ntcnst(1,cotc0%mcnstr)=4
          i1=INDEX(line,'DIST')+4
          CALL readsi(line,i1,iout,idummy,erread)
          ntcnst(2,cotc0%mcnstr)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          ntcnst(3,cotc0%mcnstr)=NAT_cpmd(idummy)
          i1=iout
          CALL readsr(line,i1,iout,cnsval(cotc0%mcnstr),erread)
          i1=iout
          IF (cotc0%lshove) CALL readsi(line,i1,iout,ntcnst(6,cotc0%mcnstr),erread)
          i1=INDEX(line,'GROWTH')
          IF (i1.NE.0) THEN
             i1=i1+6
             CALL readsr(line,i1,iout,grate(cotc0%mcnstr),erread)
             i1=INDEX(line,'DEST')
             IF (i1.NE.0) THEN
                i1=i1+4
                CALL readsr(line,i1,iout,cnsval_dest(cotc0%mcnstr),erread)
             ENDIF
          ENDIF
       ELSEIF (INDEX(line,'DISAXIS').NE.0) THEN
          cotc0%mcnstr=cotc0%mcnstr+1
          ntcnst(1,cotc0%mcnstr)=12
          i1=INDEX(line,'DISAXIS')+7
          CALL readsi(line,i1,iout,idummy,erread)
          ntcnst(2,cotc0%mcnstr)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          ntcnst(3,cotc0%mcnstr)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          ntcnst(4,cotc0%mcnstr)=idummy
          i1=iout
          CALL readsr(line,i1,iout,cnsval(cotc0%mcnstr),erread)
          i1=iout
          IF (cotc0%lshove) CALL readsi(line,i1,iout,ntcnst(6,cotc0%mcnstr),erread)
          i1=INDEX(line,'GROWTH')
          IF (i1.NE.0) THEN
             i1=i1+6
             CALL readsr(line,i1,iout,grate(cotc0%mcnstr),erread)
             i1=INDEX(line,'DEST')
             IF (i1.NE.0) THEN
                i1=i1+4
                CALL readsr(line,i1,iout,cnsval_dest(cotc0%mcnstr),erread)
             ENDIF
          ENDIF
          ! ==------------------- SM -----------------------------------==
       ELSEIF (INDEX(line,'DIFFER').NE.0) THEN
          cotc0%mcnstr=cotc0%mcnstr+1
          ntcnst(1,cotc0%mcnstr)=7
          i1=INDEX(line,'DIFFER')+6
          CALL readsi(line,i1,iout,idummy,erread)
          ntcnst(2,cotc0%mcnstr)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          ntcnst(3,cotc0%mcnstr)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          ntcnst(4,cotc0%mcnstr)=NAT_cpmd(idummy)
          i1=iout
          CALL readsr(line,i1,iout,cnsval(cotc0%mcnstr),erread)
          i1=iout
          IF (cotc0%lshove) CALL readsi(line,i1,iout,ntcnst(6,cotc0%mcnstr),erread)
          i1=INDEX(line,'GROWTH')
          IF (i1.NE.0) THEN
             i1=i1+6
             CALL readsr(line,i1,iout,grate(cotc0%mcnstr),erread)
             i1=INDEX(line,'DEST')
             IF (i1.NE.0) THEN
                i1=i1+4
                CALL readsr(line,i1,iout,cnsval_dest(cotc0%mcnstr),erread)
             ENDIF
          ENDIF
          ! ==------------------ END SM ----------------------------------==
       ELSEIF (INDEX(line,'STRETCH').NE.0) THEN
          cotc0%mcnstr=cotc0%mcnstr+1
          ntcnst(1,cotc0%mcnstr)=1
          i1=INDEX(line,'STRETCH')+7
          CALL readsi(line,i1,iout,idummy,erread)
          ntcnst(2,cotc0%mcnstr)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          ntcnst(3,cotc0%mcnstr)=NAT_cpmd(idummy)
          i1=iout
          CALL readsr(line,i1,iout,cnsval(cotc0%mcnstr),erread)
          i1=iout
          IF (cotc0%lshove) CALL readsi(line,i1,iout,ntcnst(6,cotc0%mcnstr),erread)
          i1=INDEX(line,'GROWTH')
          IF (i1.NE.0) THEN
             i1=i1+6
             CALL readsr(line,i1,iout,grate(cotc0%mcnstr),erread)
             i1=INDEX(line,'DEST')
             IF (i1.NE.0) THEN
                i1=i1+4
                CALL readsr(line,i1,iout,cnsval_dest(cotc0%mcnstr),erread)
             ENDIF
          ENDIF
       ELSEIF (INDEX(line,'BEND').NE.0) THEN
          cotc0%mcnstr=cotc0%mcnstr+1
          ntcnst(1,cotc0%mcnstr)=2
          i1=INDEX(line,'BEND')+4
          CALL readsi(line,i1,iout,idummy,erread)
          ntcnst(2,cotc0%mcnstr)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          ntcnst(3,cotc0%mcnstr)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          ntcnst(4,cotc0%mcnstr)=NAT_cpmd(idummy)
          i1=iout
          CALL readsr(line,i1,iout,cnsval(cotc0%mcnstr),erread)
          CALL raddeg(cnsval(cotc0%mcnstr),-1)
          i1=iout
          IF (cotc0%lshove) CALL readsi(line,i1,iout,ntcnst(6,cotc0%mcnstr),erread)
          i1=INDEX(line,'GROWTH')
          IF (i1.NE.0) THEN
             i1=i1+6
             CALL readsr(line,i1,iout,grate(cotc0%mcnstr),erread)
             CALL raddeg(grate(cotc0%mcnstr),-1)
             i1=INDEX(line,'DEST')
             IF (i1.NE.0) THEN
                i1=i1+4
                CALL readsr(line,i1,iout,cnsval_dest(cotc0%mcnstr),erread)
                CALL raddeg(cnsval_dest(cotc0%mcnstr),-1)
             ENDIF
          ENDIF
       ELSEIF (INDEX(line,'TORSION').NE.0) THEN
          cotc0%mcnstr=cotc0%mcnstr+1
          ntcnst(1,cotc0%mcnstr)=3
          i1=INDEX(line,'TORSION')+7
          CALL readsi(line,i1,iout,idummy,erread)
          ntcnst(2,cotc0%mcnstr)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          ntcnst(3,cotc0%mcnstr)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          ntcnst(4,cotc0%mcnstr)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          ntcnst(5,cotc0%mcnstr)=NAT_cpmd(idummy)
          i1=iout
          CALL readsr(line,i1,iout,cnsval(cotc0%mcnstr),erread)
          CALL raddeg(cnsval(cotc0%mcnstr),-1)
          i1=iout
          IF (cotc0%lshove) CALL readsi(line,i1,iout,ntcnst(6,cotc0%mcnstr),erread)
          i1=INDEX(line,'GROWTH')
          IF (i1.NE.0) THEN
             i1=i1+6
             CALL readsr(line,i1,iout,grate(cotc0%mcnstr),erread)
             CALL raddeg(grate(cotc0%mcnstr),-1)
             i1=INDEX(line,'DEST')
             IF (i1.NE.0) THEN
                i1=i1+4
                CALL readsr(line,i1,iout,cnsval_dest(cotc0%mcnstr),erread)
                CALL raddeg(cnsval_dest(cotc0%mcnstr),-1)
             ENDIF
          ENDIF
       ELSEIF (INDEX(line,'OUTP').NE.0) THEN
          cotc0%mcnstr=cotc0%mcnstr+1
          ntcnst(1,cotc0%mcnstr)=5
          i1=INDEX(line,'OUTP')+4
          CALL readsi(line,i1,iout,idummy,erread)
          ntcnst(2,cotc0%mcnstr)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          ntcnst(3,cotc0%mcnstr)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          ntcnst(4,cotc0%mcnstr)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          ntcnst(5,cotc0%mcnstr)=NAT_cpmd(idummy)
          i1=iout
          CALL readsr(line,i1,iout,cnsval(cotc0%mcnstr),erread)
          CALL raddeg(cnsval(cotc0%mcnstr),-1)
          i1=iout
          IF (cotc0%lshove) CALL readsi(line,i1,iout,ntcnst(6,cotc0%mcnstr),erread)
          i1=INDEX(line,'GROWTH')
          IF (i1.NE.0) THEN
             i1=i1+6
             CALL readsr(line,i1,iout,grate(cotc0%mcnstr),erread)
             CALL raddeg(grate(cotc0%mcnstr),-1)
             i1=INDEX(line,'DEST')
             IF (i1.NE.0) THEN
                i1=i1+4
                CALL readsr(line,i1,iout,cnsval_dest(cotc0%mcnstr),erread)
                CALL raddeg(cnsval_dest(cotc0%mcnstr),-1)
             ENDIF
          ENDIF
       ELSEIF (INDEX(line,'RIGID').NE.0) THEN
          i1=INDEX(line,'RIGID')+5
          CALL readsi(line,i1,iout,nrig,erread)
          IF (erread) GOTO 20
          ALLOCATE(narig(nrig),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          DO i=1,nrig
             i1=iout
             CALL readsi(line,i1,iout,idummy,erread)
             narig(i)=NAT_cpmd(idummy)
             IF (erread) THEN
                ! We read the next line
                IF (paral%io_parent)&
                     READ(iunit,err=20,END=20,fmt='(A)') line
                i1=1
                CALL readsi(line,i1,iout,idummy,erread)
                narig(i)=NAT_cpmd(idummy)
                IF (erread) GOTO 20
             ENDIF
          ENDDO
          i1=iout
          IF (cotc0%lshove) CALL readsi(line,i1,iout,ntcnst(6,cotc0%mcnstr),erread)
          CALL setrig(cotc0%mcnstr,nrig,narig,ntcnst)
          DEALLOCATE(narig,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
       ELSEIF (INDEX(line,'COORD').NE.0) THEN
          cotc0%mcnstr=cotc0%mcnstr+1
          ntcnst(1,cotc0%mcnstr)=6
          i1=INDEX(line,'COORD')+5
          CALL readsi(line,i1,iout,idummy,erread)
          ntcnst(2,cotc0%mcnstr)=NAT_cpmd(idummy)
          i1=iout
          CALL readsr(line,i1,iout,cnpar(1,cotc0%mcnstr),erread)
          i1=iout
          CALL readsr(line,i1,iout,cnpar(2,cotc0%mcnstr),erread)
          i1=iout
          CALL readsr(line,i1,iout,cnsval(cotc0%mcnstr),erread)
          i1=iout
          IF (cotc0%lshove) CALL readsi(line,i1,iout,ntcnst(6,cotc0%mcnstr),erread)
          i1=INDEX(line,'GROWTH')
          IF (i1.NE.0) THEN
             i1=i1+6
             CALL readsr(line,i1,iout,grate(cotc0%mcnstr),erread)
             i1=INDEX(line,'DEST')
             IF (i1.NE.0) THEN
                i1=i1+4
                CALL readsr(line,i1,iout,cnsval_dest(cotc0%mcnstr),erread)
             ENDIF
          ENDIF
          ! ==------------ SPECIES DEPENDENT COORDINATION NUMBER  -----------==
          ! ==                 f = 1/(1+exp(k*(R-R_0)))                      ==
       ELSEIF (INDEX(line,'COORSP').NE.0) THEN
          cotc0%mcnstr=cotc0%mcnstr+1
          ntcnst(1,cotc0%mcnstr)= 8
          i1=INDEX(line,'COORSP')+6
          CALL readsi(line,i1,iout,idummy,erread)
          ntcnst(2,cotc0%mcnstr)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,ntcnst(3,cotc0%mcnstr),erread)
          i1=iout
          CALL readsr(line,i1,iout,cnpar(1,cotc0%mcnstr),erread)
          i1=iout
          CALL readsr(line,i1,iout,cnpar(2,cotc0%mcnstr),erread)
          i1=iout
          CALL readsr(line,i1,iout,cnsval(cotc0%mcnstr),erread)
          i1=iout
          IF (cotc0%lshove) CALL readsi(line,i1,iout,ntcnst(6,cotc0%mcnstr),erread)
          i1=INDEX(line,'GROWTH')
          IF (i1.NE.0) THEN
             i1=i1+6
             CALL readsr(line,i1,iout,grate(cotc0%mcnstr),erread)
             i1=INDEX(line,'DEST')
             IF (i1.NE.0) THEN
                i1=i1+4
                CALL readsr(line,i1,iout,cnsval_dest(cotc0%mcnstr),erread)
             ENDIF
          ENDIF
          ! ==------------ SPECIES DEPENDENT COORDINATION NUMBER -------------==
          ! ==            f =sum_i (1+(R/R_0)^n)/(1+(R/R_0)^(n+m))
       ELSEIF (INDEX(line,'COOR_RF').NE.0) THEN
          cotc0%mcnstr=cotc0%mcnstr+1
          ntcnst(1,cotc0%mcnstr)= 9
          i1=INDEX(line,'COOR_RF')+7
          CALL readsi(line,i1,iout,idummy,erread)
          ntcnst(2,cotc0%mcnstr)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,ntcnst(3,cotc0%mcnstr),erread)
          i1=iout
          CALL readsi(line,i1,iout,ntcnst(4,cotc0%mcnstr),erread)
          i1=iout
          CALL readsi(line,i1,iout,ntcnst(5,cotc0%mcnstr),erread)
          i1=iout
          CALL readsr(line,i1,iout,cnpar(1,cotc0%mcnstr),erread)
          i1=iout
          ! II=INDEX(LINE,'2SHELL')
          ! IF(II.NE.0) THEN
          ! CALL READSI(LINE,I1,IOUT,CNPAR(2,MCNSTR),ERREAD)
          ! I1=IOUT
          ! ENDIF

          CALL readsr(line,i1,iout,cnsval(cotc0%mcnstr),erread)
          i1=iout
          IF (cotc0%lshove) CALL readsi(line,i1,iout,ntcnst(6,cotc0%mcnstr),erread)
          i1=INDEX(line,'GROWTH')
          IF (i1.NE.0) THEN
             i1=i1+6
             CALL readsr(line,i1,iout,grate(cotc0%mcnstr),erread)
             i1=INDEX(line,'DEST')
             IF (i1.NE.0) THEN
                i1=i1+4
                CALL readsr(line,i1,iout,cnsval_dest(cotc0%mcnstr),erread)
             ENDIF
          ENDIF
          ! ==--------------- BOND SWITCH WITH RATIONAL F -------------------==
          ! ==            f = (1+(R/R_0)^n)/(1+(R/R_0)^(n+m))
       ELSEIF (INDEX(line,'BNSWT').NE.0) THEN
          cotc0%mcnstr=cotc0%mcnstr+1
          ntcnst(1,cotc0%mcnstr)= 10
          i1=INDEX(line,'BNSWT')+5
          CALL readsi(line,i1,iout,idummy,erread)
          ntcnst(2,cotc0%mcnstr)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          ntcnst(3,cotc0%mcnstr)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,ntcnst(4,cotc0%mcnstr),erread)
          i1=iout
          i1=iout
          CALL readsi(line,i1,iout,ntcnst(5,cotc0%mcnstr),erread)
          i1=iout
          CALL readsr(line,i1,iout,cnpar(1,cotc0%mcnstr),erread)
          i1=iout
          CALL readsr(line,i1,iout,cnsval(cotc0%mcnstr),erread)
          i1=iout
          IF (cotc0%lshove) CALL readsi(line,i1,iout,ntcnst(6,cotc0%mcnstr),erread)
          i1=INDEX(line,'GROWTH')
          IF (i1.NE.0) THEN
             i1=i1+6
             CALL readsr(line,i1,iout,grate(cotc0%mcnstr),erread)
             i1=INDEX(line,'DEST')
             IF (i1.NE.0) THEN
                i1=i1+4
                CALL readsr(line,i1,iout,cnsval_dest(cotc0%mcnstr),erread)
             ENDIF
          ENDIF
          ! ==--------- TOTAL SPECIES DEPENDENT COORDINATION NUMBER ---------==
          ! ==            f =SUM_j[sum_i (1+(R/R_0)^n)/(1+(R/R_0)^(n+m))]    ==
       ELSEIF (INDEX(line,'TOT_COOR').NE.0) THEN
          cotc0%mcnstr = cotc0%mcnstr + 1
          ntcnst(1,cotc0%mcnstr) = 11
          i1=INDEX(line,'TOT_COOR') +8
          CALL readsi(line,i1,iout,ntcnst(2,cotc0%mcnstr),erread)
          i1=iout
          CALL readsi(line,i1,iout,ntcnst(3,cotc0%mcnstr),erread)
          i1=iout
          CALL readsi(line,i1,iout,ntcnst(4,cotc0%mcnstr),erread)
          i1=iout
          CALL readsi(line,i1,iout,ntcnst(5,cotc0%mcnstr),erread)
          i1=iout
          CALL readsr(line,i1,iout,cnpar(1,cotc0%mcnstr),erread)
          i1=iout
          ! II=INDEX(LINE,'2SHELL')
          ! IF(II.NE.0) THEN
          ! CALL READSI(LINE,I1,IOUT,CNPAR(2,MCNSTR),ERREAD)
          ! I1=IOUT
          ! ENDIF
          CALL readsr(line,i1,iout,cnsval(cotc0%mcnstr),erread)
          i1=iout
          IF (cotc0%lshove) CALL readsi(line,i1,iout,ntcnst(6,cotc0%mcnstr),erread)
          i1=INDEX(line,'GROWTH')
          IF (i1.NE.0) THEN
             i1=i1+6
             CALL readsr(line,i1,iout,grate(cotc0%mcnstr),erread)
             i1=INDEX(line,'DEST')
             IF (i1.NE.0) THEN
                i1=i1+4
                CALL readsr(line,i1,iout,cnsval_dest(cotc0%mcnstr),erread)
             ENDIF
          ENDIF
       ELSEIF (INDEX(line,'MMFIELD').NE.0) THEN
          ! FIXME. this type of contraint seems not to be used.
          ! initially, this was set to type no. 8, but that is
          ! now spoken for. from now on all QM/MM special contraints
          ! start with a type no. > 100.  AK 2005/05/12.
          cotc0%mcnstr=cotc0%mcnstr+1
          ntcnst(1,cotc0%mcnstr)= 101
          i1=INDEX(line,'MMFIELD')+7
          CALL readsi(line,i1,iout,ntcnst(2,cotc0%mcnstr),erread)
          i1=iout
          CALL readsr(line,i1,iout,cnsval(cotc0%mcnstr),erread)
       ENDIF
       IF (cotc0%mcnstr.LT.nfix) THEN
          GOTO 11
       ELSE
          GOTO 10
       ENDIF
       ! ==-- END OF FIX STRUCTURES -------------------------------------==
    ELSEIF (INDEX(line,'RESTRAINT').NE.0) THEN
       cotr007%lhyperplane=.FALSE.
       IF (INDEX(line,'HYPER').NE.0) THEN
          cotr007%lhyperplane=.TRUE.
          kplane=1.0_real_8
          i1=INDEX(line,'K=')
          IF (i1.NE.0) THEN
             i1=i1+2
             CALL readsr(line,i1,iout,kplane,erread)
          ENDIF
       ENDIF
       IF (paral%io_parent) READ(iunit,err=20,END=20,fmt=*) nres
       ALLOCATE(ntrest(6,nres),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(resval(nres),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(resfor(nres),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(gsrate(nres),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(respar(2,nres),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(resval_dest(nres),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(respos(3,nres),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(gsrate)!,nres)
       CALL zeroing(ntrest)!,6*nres)
       DO i=1,nres
          resval(i)=-999._real_8
          resval_dest(i)=-999._real_8
       ENDDO
       cotr007%mrestr=0
       ! ==------------------------------------------------------------==
       ! Read each structure.
21     CONTINUE
       IF (paral%io_parent)&
            READ(iunit,err=20,END=20,fmt='(A)') line
       IF (INDEX(line,'END').NE.0 .AND. INDEX(line,'CONSTR').NE.0) THEN
          ! Problem the number of constraints MCNSTR .NE. given NFIX.
          IF (paral%io_parent)&
               WRITE(6,'(1X,A,I5)') 'CONINP! # OF GIVEN RESTRAINTS=',nreS
          IF (paral%io_parent)&
               WRITE(6,'(1X,A,I5)') 'CONINP! # OF FOUND RESTRAINTS=',cotr007%mrestr
          IF (paral%io_parent)&
               WRITE(6,'(1X,A)')    'CONINP! END CONSTRAINTS IS REACHED.'
          CALL stopgm('CONINP','ERROR WHILE READING RESTRAINTS',& 
               __LINE__,__FILE__)
          ! ==------------------- DIST ---------------------------------==
       ELSEIF (INDEX(line,'DIST').NE.0) THEN
          cotr007%mrestr=cotr007%mrestr+1
          ntrest(1,cotr007%mrestr)=4
          i1=INDEX(line,'DIST')+4
          CALL readsi(line,i1,iout,idummy,erread)
          ntrest(2,cotr007%mrestr)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          ntrest(3,cotr007%mrestr)=NAT_cpmd(idummy)
          i1=iout
          CALL readsr(line,i1,iout,resval(cotr007%mrestr),erread)
          i1=iout
          CALL readsr(line,i1,iout,resfor(cotr007%mrestr),erread)
          i1=INDEX(line,'GROWTH')
          IF (i1.NE.0) THEN
             i1=i1+6
             CALL readsr(line,i1,iout,gsrate(cotr007%mrestr),erread)
             i1=INDEX(line,'DEST')
             IF (i1.NE.0) THEN
                i1=i1+4
                CALL readsr(line,i1,iout,resval_dest(cotr007%mrestr),erread)
             ENDIF
          ENDIF
          ! ==------------------- SM -----------------------------------==
       ELSEIF (INDEX(line,'DIFFER').NE.0) THEN
          cotr007%mrestr=cotr007%mrestr+1
          ntrest(1,cotr007%mrestr)=7
          i1=INDEX(line,'DIFFER')+6
          CALL readsi(line,i1,iout,idummy,erread)
          ntrest(2,cotr007%mrestr)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          ntrest(3,cotr007%mrestr)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          ntrest(4,cotr007%mrestr)=NAT_cpmd(idummy)
          i1=iout
          CALL readsr(line,i1,iout,resval(cotr007%mrestr),erread)
          i1=iout
          CALL readsr(line,i1,iout,resfor(cotr007%mrestr),erread)
          i1=INDEX(line,'GROWTH')
          IF (i1.NE.0) THEN
             i1=i1+6
             CALL readsr(line,i1,iout,gsrate(cotr007%mrestr),erread)
             i1=INDEX(line,'DEST')
             IF (i1.NE.0) THEN
                i1=i1+4
                CALL readsr(line,i1,iout,resval_dest(cotr007%mrestr),erread)
             ENDIF
          ENDIF
          ! ==------------------ END SM ----------------------------------==
          ! ==------------------- RESPOS ---------------------------------==
       ELSEIF (INDEX(line,'RESPOS').NE.0) THEN! cmb-kk
          cotr007%mrestr=cotr007%mrestr+1
          ntrest(1,cotr007%mrestr)=13  ! RESPOS is changed to 13
          i1=INDEX(line,'RESPOS')+6
          CALL readsi(line,i1,iout,idummy,erread)
          ntrest(2,cotr007%mrestr)=NAT_cpmd(idummy)
          i1=iout
          CALL readsr(line,i1,iout,respos(1,cotr007%mrestr),erread)
          i1=iout
          CALL readsr(line,i1,iout,respos(2,cotr007%mrestr),erread)
          i1=iout
          CALL readsr(line,i1,iout,respos(3,cotr007%mrestr),erread)
          i1=iout
          CALL readsr(line,i1,iout,resval(cotr007%mrestr),erread)
          i1=iout
          CALL readsr(line,i1,iout,resfor(cotr007%mrestr),erread)
          i1=INDEX(line,'GROWTH')
          IF (i1.NE.0) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(A)')&
                  ' WARNING: RESPOS IS NOT SUPPOSED TO BE USED',&
                  ' ALONG WITH THE KEYWORD [GROWTH]'
             CALL stopgm('CONINP','ERROR WHILE READING RESTRAINTS',& 
                  __LINE__,__FILE__)
          ENDIF
          ! ==------------------ END RESPOS ------------------------------==
       ELSEIF (INDEX(line,'STRETCH').NE.0) THEN
          cotr007%mrestr=cotr007%mrestr+1
          ntrest(1,cotr007%mrestr)=1
          i1=INDEX(line,'STRETCH')+7
          CALL readsi(line,i1,iout,idummy,erread)
          ntrest(2,cotr007%mrestr)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          ntrest(3,cotr007%mrestr)=NAT_cpmd(idummy)
          i1=iout
          CALL readsr(line,i1,iout,resval(cotr007%mrestr),erread)
          i1=iout
          CALL readsr(line,i1,iout,resfor(cotr007%mrestr),erread)
          i1=INDEX(line,'GROWTH')
          IF (i1.NE.0) THEN
             i1=i1+6
             CALL readsr(line,i1,iout,gsrate(cotr007%mrestr),erread)
             i1=INDEX(line,'DEST')
             IF (i1.NE.0) THEN
                i1=i1+4
                CALL readsr(line,i1,iout,resval_dest(cotr007%mrestr),erread)
             ENDIF
          ENDIF
       ELSEIF (INDEX(line,'BEND').NE.0) THEN
          cotr007%mrestr=cotr007%mrestr+1
          ntrest(1,cotr007%mrestr)=2
          i1=INDEX(line,'BEND')+4
          CALL readsi(line,i1,iout,idummy,erread)
          ntrest(2,cotr007%mrestr)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          ntrest(3,cotr007%mrestr)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          ntrest(4,cotr007%mrestr)=NAT_cpmd(idummy)
          i1=iout
          CALL readsr(line,i1,iout,resval(cotr007%mrestr),erread)
          CALL raddeg(resval(cotr007%mrestr),-1)
          i1=iout
          CALL readsr(line,i1,iout,resfor(cotr007%mrestr),erread)
          i1=INDEX(line,'GROWTH')
          IF (i1.NE.0) THEN
             i1=i1+6
             CALL readsr(line,i1,iout,gsrate(cotr007%mrestr),erread)
             CALL raddeg(gsrate(cotr007%mrestr),-1)
             i1=INDEX(line,'DEST')
             IF (i1.NE.0) THEN
                i1=i1+4
                CALL readsr(line,i1,iout,resval_dest(cotr007%mrestr),erread)
                CALL raddeg(resval_dest(cotr007%mrestr),-1)
             ENDIF
          ENDIF
       ELSEIF (INDEX(line,'TORSION').NE.0) THEN
          cotr007%mrestr=cotr007%mrestr+1
          ntrest(1,cotr007%mrestr)=3
          i1=INDEX(line,'TORSION')+7
          CALL readsi(line,i1,iout,idummy,erread)
          ntrest(2,cotr007%mrestr)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          ntrest(3,cotr007%mrestr)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          ntrest(4,cotr007%mrestr)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          ntrest(5,cotr007%mrestr)=NAT_cpmd(idummy)
          i1=iout
          CALL readsr(line,i1,iout,resval(cotr007%mrestr),erread)
          CALL raddeg(resval(cotr007%mrestr),-1)
          i1=iout
          CALL readsr(line,i1,iout,resfor(cotr007%mrestr),erread)
          i1=INDEX(line,'GROWTH')
          IF (i1.NE.0) THEN
             i1=i1+6
             CALL readsr(line,i1,iout,gsrate(cotr007%mrestr),erread)
             CALL raddeg(gsrate(cotr007%mrestr),-1)
             i1=INDEX(line,'DEST')
             IF (i1.NE.0) THEN
                i1=i1+4
                CALL readsr(line,i1,iout,resval_dest(cotr007%mrestr),erread)
                CALL raddeg(resval_dest(cotr007%mrestr),-1)
             ENDIF
          ENDIF
       ELSEIF (INDEX(line,'OUTP').NE.0) THEN
          cotr007%mrestr=cotr007%mrestr+1
          ntrest(1,cotr007%mrestr)=5
          i1=INDEX(line,'OUTP')+4
          CALL readsi(line,i1,iout,idummy,erread)
          ntrest(2,cotr007%mrestr)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          ntrest(3,cotr007%mrestr)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          ntrest(4,cotr007%mrestr)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          ntrest(5,cotr007%mrestr)=NAT_cpmd(idummy)
          i1=iout
          CALL readsr(line,i1,iout,resval(cotr007%mrestr),erread)
          CALL raddeg(resval(cotr007%mrestr),-1)
          i1=iout
          CALL readsr(line,i1,iout,resfor(cotr007%mrestr),erread)
          i1=INDEX(line,'GROWTH')
          IF (i1.NE.0) THEN
             i1=i1+6
             CALL readsr(line,i1,iout,gsrate(cotr007%mrestr),erread)
             CALL raddeg(gsrate(cotr007%mrestr),-1)
             i1=INDEX(line,'DEST')
             IF (i1.NE.0) THEN
                i1=i1+4
                CALL readsr(line,i1,iout,resval_dest(cotr007%mrestr),erread)
                CALL raddeg(resval_dest(cotr007%mrestr),-1)
             ENDIF
          ENDIF
       ELSEIF (INDEX(line,'COORD').NE.0) THEN
          cotr007%mrestr=cotr007%mrestr+1
          ntrest(1,cotr007%mrestr)=6
          i1=INDEX(line,'COORD')+5
          CALL readsi(line,i1,iout,idummy,erread)
          ntrest(2,cotr007%mrestr)=NAT_cpmd(idummy)
          i1=iout
          CALL readsr(line,i1,iout,respar(1,cotr007%mrestr),erread)
          i1=iout
          CALL readsr(line,i1,iout,respar(2,cotr007%mrestr),erread)
          i1=iout
          CALL readsr(line,i1,iout,resval(cotr007%mrestr),erread)
          i1=iout
          CALL readsr(line,i1,iout,resfor(cotr007%mrestr),erread)
          i1=INDEX(line,'GROWTH')
          IF (i1.NE.0) THEN
             i1=i1+6
             CALL readsr(line,i1,iout,gsrate(cotr007%mrestr),erread)
             i1=INDEX(line,'DEST')
             IF (i1.NE.0) THEN
                i1=i1+4
                CALL readsr(line,i1,iout,resval_dest(cotr007%mrestr),erread)
             ENDIF
          ENDIF
       ELSEIF (INDEX(line,'COORSP').NE.0) THEN
          cotr007%mrestr=cotr007%mrestr+1
          ntrest(1,cotr007%mrestr)=8
          i1=INDEX(line,'COORSP')+6
          CALL readsi(line,i1,iout,idummy,erread)
          ntrest(2,cotr007%mrestr)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,ntrest(3,cotr007%mrestr),erread)
          i1=iout
          CALL readsr(line,i1,iout,respar(1,cotr007%mrestr),erread)
          i1=iout
          CALL readsr(line,i1,iout,respar(2,cotr007%mrestr),erread)
          i1=iout
          CALL readsr(line,i1,iout,resval(cotr007%mrestr),erread)
          i1=iout
          CALL readsr(line,i1,iout,resfor(cotr007%mrestr),erread)
          i1=INDEX(line,'GROWTH')
          IF (i1.NE.0) THEN
             i1=i1+6
             CALL readsr(line,i1,iout,gsrate(cotr007%mrestr),erread)
             i1=INDEX(line,'DEST')
             IF (i1.NE.0) THEN
                i1=i1+4
                CALL readsr(line,i1,iout,resval_dest(cotr007%mrestr),erread)
             ENDIF
          ENDIF
       ELSEIF (INDEX(line,'COOR_RF').NE.0) THEN
          cotr007%mrestr=cotr007%mrestr+1
          ntrest(1,cotr007%mrestr)=9
          i1=INDEX(line,'COOR_RF')+7
          CALL readsi(line,i1,iout,idummy,erread)
          ntrest(2,cotr007%mrestr)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,ntrest(3,cotr007%mrestr),erread)
          i1=iout
          CALL readsi(line,i1,iout,ntrest(4,cotr007%mrestr),erread)
          i1=iout
          CALL readsi(line,i1,iout,ntrest(5,cotr007%mrestr),erread)
          i1=iout
          CALL readsr(line,i1,iout,respar(1,cotr007%mrestr),erread)
          i1=iout
          ! II=INDEX(LINE,'2SHELL')
          ! IF(II.NE.0) THEN
          ! CALL READSR(LINE,I1,IOUT,RESPAR(2,MRESTR),ERREAD)
          ! I1=IOUT
          ! ENDIF
          CALL readsr(line,i1,iout,resval(cotr007%mrestr),erread)
          i1=iout
          CALL readsr(line,i1,iout,resfor(cotr007%mrestr),erread)
          i1=INDEX(line,'GROWTH')
          IF (i1.NE.0) THEN
             i1=i1+6
             CALL readsr(line,i1,iout,gsrate(cotr007%mrestr),erread)
             i1=INDEX(line,'DEST')
             IF (i1.NE.0) THEN
                i1=i1+4
                CALL readsr(line,i1,iout,resval_dest(cotr007%mrestr),erread)
             ENDIF
          ENDIF
       ELSEIF (INDEX(line,'BNSWT').NE.0) THEN
          cotr007%mrestr=cotr007%mrestr+1
          ntrest(1,cotr007%mrestr)= 10
          i1=INDEX(line,'BNSWT')+5
          CALL readsi(line,i1,iout,idummy,erread)
          ntrest(2,cotr007%mrestr)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          ntrest(3,cotr007%mrestr)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,ntrest(4,cotr007%mrestr),erread)
          i1=iout
          i1=iout
          CALL readsi(line,i1,iout,ntrest(5,cotr007%mrestr),erread)
          i1=iout
          CALL readsr(line,i1,iout,respar(1,cotr007%mrestr),erread)
          i1=iout
          CALL readsr(line,i1,iout,resval(cotr007%mrestr),erread)
          i1=iout
          CALL readsr(line,i1,iout,resfor(cotr007%mrestr),erread)
          i1=INDEX(line,'GROWTH')
          IF (i1.NE.0) THEN
             i1=i1+6
             CALL readsr(line,i1,iout,gsrate(cotr007%mrestr),erread)
             i1=INDEX(line,'DEST')
             IF (i1.NE.0) THEN
                i1=i1+4
                CALL readsr(line,i1,iout,resval_dest(cotr007%mrestr),erread)
             ENDIF
          ENDIF
       ELSEIF (INDEX(line,'TOT_COOR').NE.0) THEN
          cotr007%mrestr = cotr007%mrestr + 1
          ntrest(1,cotr007%mrestr) = 11
          i1=INDEX(line,'TOT_COOR') +8
          CALL readsi(line,i1,iout,ntrest(2,cotr007%mrestr),erread)
          i1=iout
          CALL readsi(line,i1,iout,ntrest(3,cotr007%mrestr),erread)
          i1=iout
          CALL readsi(line,i1,iout,ntrest(4,cotr007%mrestr),erread)
          i1=iout
          CALL readsi(line,i1,iout,ntrest(5,cotr007%mrestr),erread)
          i1=iout
          CALL readsr(line,i1,iout,respar(1,cotr007%mrestr),erread)
          i1=iout
          ! II=INDEX(LINE,'2SHELL')
          ! IF(II.NE.0) THEN
          ! CALL READSI(LINE,I1,IOUT,RESPAR(2,MRESTR),ERREAD)
          ! I1=IOUT
          ! ENDIF
          CALL readsr(line,i1,iout,resval(cotr007%mrestr),erread)
          i1=iout
          CALL readsr(line,i1,iout,resfor(cotr007%mrestr),erread)
          i1=INDEX(line,'GROWTH')
          IF (i1.NE.0) THEN
             i1=i1+6
             CALL readsr(line,i1,iout,gsrate(cotr007%mrestr),erread)
             i1=INDEX(line,'DEST')
             IF (i1.NE.0) THEN
                i1=i1+4
                CALL readsr(line,i1,iout,resval_dest(cotr007%mrestr),erread)
             ENDIF
          ENDIF
       ELSEIF (INDEX(line,'DISAXIS').NE.0) THEN
          cotr007%mrestr=cotr007%mrestr+1
          ntrest(1,cotr007%mrestr)=12
          i1=INDEX(line,'DISAXIS')+7
          CALL readsi(line,i1,iout,idummy,erread)
          ntrest(2,cotr007%mrestr)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          ntrest(3,cotr007%mrestr)=NAT_cpmd(idummy)
          i1=iout
          CALL readsi(line,i1,iout,idummy,erread)
          ntrest(4,cotr007%mrestr)=idummy
          i1=iout
          CALL readsr(line,i1,iout,resval(cotr007%mrestr),erread)
          i1=iout
          CALL readsr(line,i1,iout,resfor(cotr007%mrestr),erread)
          i1=INDEX(line,'GROWTH')
          IF (i1.NE.0) THEN
             i1=i1+6
             CALL readsr(line,i1,iout,gsrate(cotr007%mrestr),erread)
             i1=INDEX(line,'DEST')
             IF (i1.NE.0) THEN
                i1=i1+4
                CALL readsr(line,i1,iout,resval_dest(cotr007%mrestr),erread)
             ENDIF
          ENDIF
       ENDIF
       IF (cotr007%mrestr.LT.nres) THEN
          GOTO 21
       ELSE
          ! read new restraint target values from file RESVAL
          IF (paral%io_parent)&
               CALL fileopen(91,'RESVAL',fo_old,ferror)
          IF (.NOT.ferror) THEN
             IF (paral%io_parent)&
                  WRITE(6,*)'CONINP| READING RESTRAINT TARGET VALUES ',&
                  'FROM FILE "RESVAL"'
             IF (paral%io_parent)&
                  READ(91,err=20,END=20,fmt=*) (resval(i),i=1,cotr007%mrestr)
             DO i=1,cotr007%mrestr
                IF (ntrest(1,i).EQ.2&
                     .OR.ntrest(1,i).EQ.3&
                     .OR.ntrest(1,i).EQ.5)&
                     CALL raddeg(resval(i),-1)
             ENDDO
             IF (paral%io_parent)&
                  CALL fileclose(91)
          ENDIF
          IF (cotr007%lhyperplane) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(2A,F6.2)')' CONINP| USING HYPERPLANE TYPE ',&
                  'RESTRAINTS WITH SCALING FACTOR: ',kplane
             DO i=1,cotr007%mrestr
                resfor(i)=resfor(i)*kplane
             ENDDO
          ENDIF
          GOTO 10
       ENDIF
       ! ==-- END OF RESTRAINTS -----------------------------------------==
    ELSEIF (INDEX(line,'PENALTY').NE.0) THEN
       IF (paral%io_parent)&
            READ(iunit,err=20,END=20,fmt=*) dsigma,bsigma,tsigma
       GOTO 10
    ELSE
       GOTO 10
    ENDIF
20  CONTINUE
    IF (paral%io_parent)&
         WRITE(6,*) ' ERROR WHILE READING CONSTRAINTS SECTION'
    CALL stopgm('CONINP','INPUT FILE ERROR',& 
         __LINE__,__FILE__)
30  CONTINUE
    ! Free memory
    DEALLOCATE(lfix,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (lqmmm%qmmm) CALL mm_dim(mm_revert,status)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE coninp
  ! ==================================================================
  SUBROUTINE raddeg(a,itag)
    ! ==--------------------------------------------------------------==

    REAL(real_8)                             :: a
    INTEGER                                  :: itag

! ==--------------------------------------------------------------==

    IF (a.EQ.-999._real_8) RETURN
    IF (itag.EQ.1) THEN
       a=a/pi*180._real_8
    ELSEIF (itag.EQ.-1) THEN
       a=a/180.0_real_8*pi
    ELSE
       CALL stopgm('RADDEG','OPTION',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE raddeg
  ! ==================================================================
  SUBROUTINE setrig(m,nrig,narig,ntcnst)
    ! ==--------------------------------------------------------------==
    INTEGER :: m, nrig, narig(nrig), ntcnst(6,m+nrig-1+nrig-2+nrig-3)

    INTEGER                                  :: i

! Variables
! ==--------------------------------------------------------------==
! Stretch

    DO i=2,nrig
       m=m+1
       ntcnst(1,m)=1
       ntcnst(2,m)=narig(i-1)
       ntcnst(3,m)=narig(i)
    ENDDO
    ! Bend   
    DO i=3,nrig
       m=m+1
       ntcnst(1,m)=2
       ntcnst(2,m)=narig(i-2)
       ntcnst(3,m)=narig(i-1)
       ntcnst(4,m)=narig(i)
    ENDDO
    ! Torsion
    DO i=4,nrig
       m=m+1
       ntcnst(1,m)=3
       ntcnst(2,m)=narig(i-3)
       ntcnst(3,m)=narig(i-2)
       ntcnst(4,m)=narig(i-1)
       ntcnst(5,m)=narig(i)
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE setrig
  ! ==================================================================
END MODULE coninp_utils
