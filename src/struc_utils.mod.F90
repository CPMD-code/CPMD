MODULE struc_utils
  USE adat,                            ONLY: covrad,&
                                             elem
  USE cnst,                            ONLY: fbohr,&
                                             pi
  USE coninp_utils,                    ONLY: raddeg
  USE constr_utils,                    ONLY: funct
  USE empf,                            ONLY: c,&
                                             dp12,&
                                             ibind,&
                                             naty,&
                                             zan,&
                                             zero
  USE empfor_utils,                    ONLY: cpmdbonds
  USE error_handling,                  ONLY: stopgm
  USE fstart_utils,                    ONLY: dist,&
                                             normal,&
                                             winkel
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: paral
  USE struc,                           ONLY: angle,&
                                             bond,&
                                             dihedral
  USE system,                          ONLY: cntl,&
                                             nassel,&
                                             nssel
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: pstruc

CONTAINS

  ! ==================================================================
  SUBROUTINE pstruc(tau0)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'pstruc'

    INTEGER                                  :: i, ierr, j, k, ncc, nma, nmpa1
    INTEGER, ALLOCATABLE                     :: nbonds(:)
    REAL(real_8)                             :: rcov(0:104), xmfac

    IF (.NOT.bond.AND..NOT.angle.AND..NOT.dihedral) RETURN
    nma=ions1%nat
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
    xmfac=1.35_real_8
    CALL cpmdbonds(rcov,xmfac,nbonds)
    ! ==--------------------------------------------------------------==
    IF (cntl%tssel) CALL select_bonds(nbonds)
    IF (bond) CALL pbonds
    IF (angle) CALL pangles(nbonds)
    IF (dihedral) CALL pdihed(nbonds)
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
    RETURN
  END SUBROUTINE pstruc
  ! ==================================================================
  SUBROUTINE select_bonds(nbonds)
    ! ==--------------------------------------------------------------==
    ! Variables
    INTEGER                                  :: nbonds(*)

    INTEGER                                  :: ia, ib
    LOGICAL                                  :: il, jl

! ==--------------------------------------------------------------==

    il=findla(1)
    IF (.NOT.il) nbonds(1)=0
    DO ia=2,naty
       il=findla(ia)
       IF (.NOT.il) nbonds(ia)=0
       DO ib=1,ia-1
          jl=findla(ib)
          IF (.NOT.il .OR. .NOT.jl) ibind(ib+ia*(ia-1)/2)=0
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE select_bonds
  ! ==================================================================
  FUNCTION findla(ia) RESULT(my_res)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ia
    LOGICAL                                  :: my_res

    INTEGER                                  :: i

    my_res=.FALSE.
    DO i=1,nssel
       IF (nassel(i).EQ.ia) THEN
          my_res=.TRUE.
          GOTO 100
       ENDIF
    ENDDO
100 CONTINUE
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION findla
  ! ==================================================================
  SUBROUTINE pbonds
    ! ==--------------------------------------------------------------==
    ! Variables
    INTEGER                                  :: ia, iab, ib, is1, is2
    REAL(real_8)                             :: dd

! ==--------------------------------------------------------------==

    IF (paral%io_parent)&
         WRITE(6,'(30X,A,/,A,A)') ' BONDS',' ATOM 1  ATOM 2  ',&
         ' TYPE 1   TYPE 2      DISTANCE (BOHR, ANGSTROM)'
    DO ia=2,naty
       is1=specia(ia)
       DO ib=1,ia-1
          is2=specia(ib)
          iab=ibind(ib+ia*(ia-1)/2)
          IF (iab.NE.0) THEN
             dd=SQRT(dist(ia,ib))
             IF (paral%io_parent)&
                  WRITE(6,'(I5,I8,7X,A2,7X,A2,9X,2F12.5)')&
                  ia,ib,elem%el(ions0%iatyp(is1)),elem%el(ions0%iatyp(is2)),dd,dd/fbohr
          ENDIF
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE pbonds
  ! ==================================================================
  FUNCTION specia(ia)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ia, specia

    INTEGER                                  :: is, ix

    specia=0
    ix=0
    DO is=1,ions1%nsp
       IF (ix+ions0%na(is).GE.ia) THEN
          specia=is
          RETURN
       ENDIF
       ix=ix+ions0%na(is)
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION specia
  ! ==================================================================
  SUBROUTINE pangles(nbonds)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nbonds(naty)

    INTEGER                                  :: i, ia, iab, iac, ib, ic, is1, &
                                                is2, is3, j
    REAL(real_8)                             :: alpha

    IF (paral%io_parent)&
         WRITE(6,'(30X,A,/,A,A)') ' ANGLES','    ATOM NUMBERS    ',&
         '    ATOM TYPES       BOND ANGLES(DEGREES)     '
    DO ia=1,naty
       IF (nbonds(ia).GE.2) THEN
          DO ib=1,naty
             i  =MAX(ia,ib)
             j  =MIN(ia,ib)
             iab=ibind(j+i*(i-1)/2)
             IF (ib.NE.ia.AND.iab.NE.0) THEN
                DO ic=ib+1,naty
                   i  =MAX(ia,ic)
                   j  =MIN(ia,ic)
                   iac=ibind(j+i*(i-1)/2)
                   IF (ic.NE.ia.AND.iac.NE.0) THEN
                      alpha=winkel(ib,ic,ia)
                      is1=specia(ia)
                      is2=specia(ib)
                      is3=specia(ic)
                      IF (paral%io_parent)&
                           WRITE(6,'(2X,3I4,8X,3(A2,4X),8X,F10.4)')&
                           ib,ia,ic,elem%el(ions0%iatyp(is2)),elem%el(ions0%iatyp(is1)),&
                           elem%el(ions0%iatyp(is3)),alpha
                   ENDIF
                ENDDO
             ENDIF
          ENDDO
       ENDIF
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE pangles
  ! ==================================================================
  SUBROUTINE pdihed(nbonds)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nbonds(naty)

    INTEGER                                  :: i, ia, iab, iac, iad, ib, &
                                                ibd, ic, id, is1, is2, is3, &
                                                is4, j
    REAL(real_8)                             :: alpha, beta, delta, gamma

    IF (paral%io_parent)&
         WRITE(6,'(30X,A,/,A,A)') ' DIHEDRALS','  TYPE    ATOM NUMBERS   ',&
         '   ATOM TYPES     DIHEDRAL ANGLES(DEG)     '
    DO ia=1,naty
       IF (nbonds(ia).GE.2) THEN
          DO ib=1,naty
             i  =MAX(ia,ib)
             j  =MIN(ia,ib)
             iab=ibind(j+i*(i-1)/2)
             IF (ib.NE.ia.AND.iab.NE.0) THEN
                DO ic=ib+1,naty
                   i  =MAX(ia,ic)
                   j  =MIN(ia,ic)
                   iac=ibind(j+i*(i-1)/2)
                   IF (ic.NE.ia.AND.iac.NE.0) THEN
                      IF (nbonds(ia).GE.3) THEN
                         DO id=ic+1,naty
                            i  =MAX(ia,id)
                            j  =MIN(ia,id)
                            iad=ibind(j+i*(i-1)/2)
                            IF (id.NE.ia.AND.iad.NE.0) THEN
                               alpha=winkel(ib,ic,ia)
                               beta =winkel(ic,id,ia)
                               gamma=winkel(ib,id,ia)
                               IF (alpha.NE.180._real_8.AND.beta.NE.180._real_8&
                                    .AND.gamma.NE.180._real_8) THEN
                                  is1=specia(ib)
                                  is2=specia(ic)
                                  is3=specia(ia)
                                  is4=specia(id)
                                  delta=outpd(ib,ic,ia,id)
                                  IF (paral%io_parent)&
                                       WRITE(6,'(2X,A,4I4,4X,4(A2,2X),4X,F10.4)')&
                                       'OUTP',ib,ic,ia,id,elem%el(ions0%iatyp(is2)),elem%el(ions0%iatyp(is1)),&
                                       elem%el(ions0%iatyp(is3)),elem%el(ions0%iatyp(is4)),delta
                               ENDIF
                            ENDIF
                         ENDDO
                      ENDIF
                   ENDIF
                ENDDO
                IF (naty.GE.4.AND.nbonds(ib).GE.2.AND.ia.LT.ib) THEN
                   DO ic=1,naty
                      i  =MAX(ic,ia)
                      j  =MIN(ic,ia)
                      iac=ibind(j+i*(i-1)/2)
                      IF (ic.NE.ia.AND.ic.NE.ib.AND.iac.NE.0) THEN
                         DO id=1,naty
                            i  =MAX(ib,id)
                            j  =MIN(ib,id)
                            ibd=ibind(j+i*(i-1)/2)
                            IF (id.NE.ia.AND.id.NE.ib.AND.id.NE.ic.AND.&
                                 ibd.NE.0) THEN
                               alpha=winkel(ic,ib,ia)
                               beta =winkel(ia,id,ib)
                               IF (alpha.NE.180._real_8.AND.beta.NE.180._real_8) THEN
                                  is1=specia(ic)
                                  is2=specia(ia)
                                  is3=specia(ib)
                                  is4=specia(id)
                                  delta=torsd(ic,ia,ib,id)
                                  IF (paral%io_parent)&
                                       WRITE(6,'(2X,A,4I4,4X,4(A2,2X),4X,F10.4)')&
                                       'TORS',ic,ia,ib,id,elem%el(ions0%iatyp(is2)),elem%el(ions0%iatyp(is1)),&
                                       elem%el(ions0%iatyp(is3)),elem%el(ions0%iatyp(is4)),delta
                               ENDIF
                            ENDIF
                         ENDDO
                      ENDIF
                   ENDDO
                ENDIF
             ENDIF
          ENDDO
       ENDIF
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE pdihed
  ! ==================================================================
  FUNCTION outpd(ia,ib,ic,id)
    INTEGER                                  :: ia, ib, ic, id
    REAL(real_8)                             :: outpd

    INTEGER                                  :: i
    REAL(real_8)                             :: oo, pih, u(3), v(3), w(3)

    DO i=1,3
       u(i)=c(i,ic)-c(i,ia)
       v(i)=c(i,ic)-c(i,ib)
    ENDDO
    CALL normal(u,v,w)
    DO i=1,3
       w(i)=w(i)+c(i,ic)
    ENDDO
    pih=zero
    CALL funct(oo,outpd,pih,w,c(1,ic),c(1,id))
    outpd=pi*dp12-outpd
    CALL raddeg(outpd,1)
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION outpd
  ! ==================================================================
  FUNCTION torsd(ia,ib,ic,id)
    INTEGER                                  :: ia, ib, ic, id
    REAL(real_8)                             :: torsd

    INTEGER                                  :: i
    REAL(real_8)                             :: pih, tt, u(3), v(3), w1(3), &
                                                w2(3)

    DO i=1,3
       u(i)=c(i,ib)-c(i,ia)
       v(i)=c(i,ib)-c(i,ic)
    ENDDO
    CALL normal(u,v,w1)
    DO i=1,3
       u(i)=c(i,ic)-c(i,ib)
       v(i)=c(i,ic)-c(i,id)
    ENDDO
    CALL normal(u,v,w2)
    DO i=1,3
       w1(i)=w1(i)+c(i,ic)
       w2(i)=w2(i)+c(i,ic)
    ENDDO
    pih=zero
    CALL funct(tt,torsd,pih,w1,c(1,ic),w2)
    CALL raddeg(torsd,1)
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION torsd
  ! ==================================================================


END MODULE struc_utils
