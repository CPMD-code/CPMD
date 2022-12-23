MODULE empfor_utils
  USE adat,                            ONLY: covrad
  USE cnst,                            ONLY: fbohr
  USE cotr,                            ONLY: boad,&
                                             cotc0
  USE empf,                            ONLY: c,&
                                             ibind,&
                                             iw,&
                                             naty,&
                                             zan
  USE error_handling,                  ONLY: stopgm
  USE fstart_utils,                    ONLY: dist,&
                                             fstart
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: paral
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: empfor
  PUBLIC :: cpmdbonds
  !public :: pbond
  PUBLIC :: smolec

CONTAINS

  ! ==================================================================
  SUBROUTINE empfor(tau0,hess)
    ! ==--------------------------------------------------------------==
    ! ==  INTERFACE TO THE "DISCO" EMPIRICAL FORCE FIELD CODE         ==
    ! ==  THE DISCO CODE IS A GIFT OF THOMAS FISCHER (IPS ETHZ)       ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:), hess(*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'empfor'

    INTEGER                                  :: i, ierr, isub, j, k, ncc, &
                                                nmpa1
    INTEGER, ALLOCATABLE                     :: nbonds(:)
    REAL(real_8)                             :: rcov(0:104), xmfac

    CALL tiset(procedureN,isub)

    nmpa1=(ions1%nat*(ions1%nat+1))/2
    ALLOCATE(zan(ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(c(3,ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(ibind(nmpa1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(nbonds(ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(ibind)!,nmpa1)
    CALL zeroing(nbonds)!,ions1%nat)
    naty=ions1%nat
    k=0
    DO i=1,ions1%nsp
       DO j=1,ions0%na(i)
          k=k+1
          zan(k)=REAL(ions0%iatyp(i),kind=real_8)
          c(1,k)=tau0(1,j,i)
          c(2,k)=tau0(2,j,i)
          c(3,k)=tau0(3,j,i)
       ENDDO
    ENDDO
    DO i=1,99
       rcov(i)=covrad(i)*fbohr
    ENDDO
    ncc=3*ions1%nat
    IF (paral%io_parent)&
         WRITE(6,'(A)') ' INITIALIZE EMPIRICAL HESSIAN '
    ! 1.35 : SCHLEGEL multiplication factor
    xmfac=1.35_real_8
    CALL cpmdbonds(rcov,xmfac,nbonds)
    CALL fstart(hess,rcov,nbonds,ncc,.FALSE.)
    ! ==--------------------------------------------------------------==
    DEALLOCATE(zan,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(c,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(ibind,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(nbonds,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    RETURN
  END SUBROUTINE empfor
  ! ==================================================================
  SUBROUTINE cpmdbonds(rcov,xmfac,nbonds)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rcov(0:104), xmfac
    INTEGER                                  :: nbonds(naty)

    CHARACTER(*), PARAMETER                  :: procedureN = 'cpmdbonds'

    INTEGER                                  :: i, i1, i2, ia, iab, ib, ierr, &
                                                ii, j, j1, j2, k, k1, k2, l, &
                                                nummol
    INTEGER, ALLOCATABLE                     :: imol(:)
    REAL(real_8)                             :: cdij, dij, dmin, xcdij

    ALLOCATE(imol(naty),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(imol)!,naty)
    ! Normal bonds
    DO i=2,naty
       DO j=1,i-1
          ii=i*(i-1)/2+j
          dij=SQRT(dist(i,j))
          cdij=xmfac*(rcov(NINT(zan(i)))+rcov(NINT(zan(j))))
          IF (cdij.GE.dij) THEN
             ibind(ii)=1
          ELSE
             ibind(ii)=0
          ENDIF
       ENDDO
    ENDDO
    ! Added or deleted bonds
    IF (cotc0%nboad.GT.0) THEN
       DO k=1,cotc0%nboad
          i=boad(1,k)
          j=boad(2,k)
          l=boad(3,k)
          IF (j.GT.i) THEN
             i=boad(2,k)
             j=boad(1,k)
          ENDIF
          ii=i*(i-1)/2+j
          IF (l.LT.0) ibind(ii)=0
          IF (l.GT.0) ibind(ii)=1
       ENDDO
    ENDIF
    ! Hydrogen bonds
    DO i=2,naty
       DO j=1,i-1
          ia=NINT(zan(i))
          ib=NINT(zan(j))
          IF ((ia.EQ.1.AND.ib.NE.1).OR.(ia.NE.1.AND.ib.EQ.1)) THEN
             ii=i*(i-1)/2+j
             dij=SQRT(dist(i,j))
             cdij=xmfac*(rcov(NINT(zan(i)))+rcov(NINT(zan(j))))
             IF (cdij.LT.dij) THEN
                xcdij=cdij/xmfac*2.40_real_8
                IF (xcdij.GT.dij) ibind(ii)=2
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    ! Search for molecules
    CALL smolec(imol,nummol)
    ! Add auxilliary bonds between molecules
    DO i1=2,nummol
       DO i2=1,i1-1
          k1=0
          k2=0
          dmin=100.0e9_real_8
          DO j1=1,naty
             IF (imol(j1).EQ.i1) THEN
                DO j2=1,naty
                   IF (imol(j2).EQ.i2) THEN
                      dij=SQRT(dist(j1,j2))
                      IF (dij.LT.dmin) THEN
                         k1=j1
                         k2=j2
                         dmin=dij
                      ENDIF
                   ENDIF
                ENDDO
             ENDIF
          ENDDO
          ia=MAX(k1,k2)
          ib=MIN(k1,k2)
          k=ib+ia*(ia-1)/2
          IF (ibind(k).EQ.0) ibind(k)=3
       ENDDO
    ENDDO
    ! Number of bonds per atom
    CALL zeroing(nbonds)!,naty)
    DO ia=2,naty
       DO ib=1,ia-1
          iab=ibind(ib+ia*(ia-1)/2)
          IF (iab.EQ.1) THEN
             nbonds(ia)=nbonds(ia)+1
             nbonds(ib)=nbonds(ib)+1
          ENDIF
       ENDDO
    ENDDO
    ! Print bonds
    CALL pbond(imol,nummol)
    DEALLOCATE(imol,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE cpmdbonds
  ! ==================================================================
  SUBROUTINE pbond(imol,nummol)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: imol(naty), nummol

    INTEGER                                  :: i, ia, iab, iau, ib, idum, &
                                                ihb, iprt(20), j, k

    iab=0
    ihb=0
    iau=0
    DO ia=2,naty
       DO ib=1,ia-1
          k = ib+ia*(ia-1)/2
          IF (ibind(k).EQ.1) iab=iab+1
          IF (ibind(k).EQ.2) ihb=ihb+1
          IF (ibind(k).EQ.3) iau=iau+1
       ENDDO
    ENDDO
    IF (iab.GT.0) THEN
       IF (paral%io_parent)&
            WRITE(iw,'(A,20X,A)') '       ','<<<<< ASSUMED BONDS >>>>>'
       idum = 0
       DO ia=2,naty
          DO ib=1,ia-1
             k = ib+ia*(ia-1)/2
             IF (ibind(k).EQ.1) THEN
                idum         = idum + 2
                iprt(idum-1) = ia
                iprt(idum)   = ib
                IF (idum.EQ.10) THEN
                   IF (paral%io_parent)&
                        WRITE(iw,1010) (iprt(j),j=1,10)
                   idum = 0
                ENDIF
             ENDIF
          ENDDO
       ENDDO
       IF ((idum.LT.10).AND.paral%io_parent)&
            WRITE(iw,1010) (iprt(j),j=1,idum)
    ENDIF
    IF (ihb.GT.0) THEN
       IF (paral%io_parent)&
            WRITE(iw,'(A,20X,A)') '       ','<<<<< HYDROGEN BONDS >>>>>'
       idum = 0
       DO ia=2,naty
          DO ib=1,ia-1
             k = ib+ia*(ia-1)/2
             IF (ibind(k).EQ.2) THEN
                idum         = idum + 2
                iprt(idum-1) = ia
                iprt(idum)   = ib
                IF (idum.EQ.10) THEN
                   IF (paral%io_parent)&
                        WRITE(iw,1010) (iprt(j),j=1,10)
                   idum = 0
                ENDIF
             ENDIF
          ENDDO
       ENDDO
       IF ((idum.LT.10).AND.paral%io_parent)&
            WRITE(iw,1010) (iprt(j),j=1,idum)
    ENDIF
    IF (iau.GT.0) THEN
       IF (paral%io_parent)&
            WRITE(iw,'(A,20X,A)') '      ','<<<<< AUXILIARY BONDS >>>>>'
       idum = 0
       DO ia=2,naty
          DO ib=1,ia-1
             k = ib+ia*(ia-1)/2
             IF (ibind(k).EQ.3) THEN
                idum         = idum + 2
                iprt(idum-1) = ia
                iprt(idum)   = ib
                IF (idum.EQ.10) THEN
                   IF (paral%io_parent)&
                        WRITE(iw,1010) (iprt(j),j=1,10)
                   idum = 0
                ENDIF
             ENDIF
          ENDDO
       ENDDO
       IF ((idum.LT.10).AND.paral%io_parent)&
            WRITE(iw,1010) (iprt(j),j=1,idum)
    ENDIF
    IF (paral%io_parent)&
         WRITE(iw,'(A,I3)') ' TOTAL NUMBER OF MOLECULAR STRUCTURES:',&
         nummol
    IF (nummol.GT.1) THEN
       IF (paral%io_parent)&
            WRITE(iw,'(A,20X,A)') '         ','<<<<< MOLECULES >>>>>'
       DO i=1,nummol
          idum=0
          DO ia=1,naty
             IF (imol(ia).EQ.i) THEN
                idum=idum+1
                iprt(idum)=ia
                IF ((idum.EQ.15).AND.paral%io_parent)&
                     WRITE(iw,1020) i,(iprt(j),j=1,idum)
                IF (idum.EQ.15) idum=0
             ENDIF
          ENDDO
          IF ((idum.LT.15).AND.paral%io_parent)&
               WRITE(iw,1020) i,(iprt(j),j=1,idum)
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
1010 FORMAT(2x,5(i3,' <-->',i3,2x))
1020 FORMAT(' MOLECULE:',i3,' <>   ',15i5)
    RETURN
  END SUBROUTINE pbond
  ! ==================================================================
  SUBROUTINE smolec(imol,nummol)
    ! ==--------------------------------------------------------------==
    ! == Calculates the TOTAL NUMBER OF MOLECULAR STRUCTURES (NUMMOL) ==
    ! == i.e. check the connection between all atoms in order         ==
    ! == to determine "isolated" molecule                             ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: imol(naty), nummol

    INTEGER, PARAMETER                       :: maxl = 2

    INTEGER                                  :: i, ipar, j, k, l

    nummol=0
10  CONTINUE
    DO i=1,naty
       IF (imol(i).EQ.0) THEN
          ipar=i
          nummol=nummol+1
          imol(ipar)=nummol
          GOTO 20
       ENDIF
    ENDDO
    GOTO 30
20  CONTINUE
    DO l=1,maxl
       DO i=1,naty
          IF (imol(i).EQ.nummol) THEN
             DO j=i+1,naty
                k = i+j*(j-1)/2
                IF (ibind(k).EQ.1) imol(j)=nummol
             ENDDO
             DO j=1,i-1
                k = j+i*(i-1)/2
                IF (ibind(k).EQ.1) imol(j)=nummol
             ENDDO
          ENDIF
       ENDDO
       DO i=naty,1,-1
          IF (imol(i).EQ.nummol) THEN
             DO j=i+1,naty
                k = i+j*(j-1)/2
                IF (ibind(k).EQ.1) imol(j)=nummol
             ENDDO
             DO j=1,i-1
                k = j+i*(i-1)/2
                IF (ibind(k).EQ.1) imol(j)=nummol
             ENDDO
          ENDIF
       ENDDO
    ENDDO
    GOTO 10
30  CONTINUE
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE smolec
  ! ==================================================================

END MODULE empfor_utils
