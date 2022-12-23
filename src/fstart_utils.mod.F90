MODULE fstart_utils
  USE cnst,                            ONLY: pi
  USE empf,                            ONLY: &
       c, ibind, iw, naty, one, three, two, zan, zero
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: paral
  USE pbc_utils,                       ONLY: pbc
  USE system,                          ONLY: cnti,&
                                             parm

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: fstart
  !public :: fbend
  !public :: flinc
  !public :: flinp
  !public :: foutp
  !public :: fstre
  !public :: ftors
  !public :: bandq
  !public :: nom
  PUBLIC :: normal
  !public :: vektor
  !public :: fxmat
  !public :: arc1
  PUBLIC :: dist
  !public :: scalaru
  PUBLIC :: winkel

CONTAINS

  ! ==================================================================
  SUBROUTINE fstart(fx,rcov,nbonds,ncc,lprt)
    ! :: :::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ! GENERATION OF AN INITIAL FORCE CONSTANT MATRIX (FX MATRIX):
    ! THF 11/05/90:
    ! :: :::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ! ==--------------------------------------------------------------==
    ! THIS CODE IS PART OF THE DISCO PROGRAM SYSTEM
    ! DISCO IS AN EXPERIMENTAL CODE FOR DIRECT SCF AND MP2 CALCULATIONS
    ! WRITTEN BY
    ! J. ALMLOF, K.FAEGRI,JR., M. FEYEREISEN, T.H. FISCHER, 
    ! K.KORSELL & H.P. LUETHI
    ! 
    ! INPUT
    ! RCOV(0:104)   : COVALENT RADIUS
    ! IBIND         : CONNECTION MATRIX (IN empf.inc)
    ! NCC           : NUMBER OF DEGREES OF FREEDOM
    ! LPRT          : .FALSE. MINIMAL PRINT OUT
    ! .TRUE.  FULL PRINT OUT
    ! OUTPUT
    ! FX(NCC,NCC)   : THE EMPIRICAL CARTESIAN FORCE CONSTANT MATRIX
    ! NBONDS(1:NATY): NUMBER OF BONDS OF ATOM NA
    ! BRED(1:NCC)   : 
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rcov(0:104)
    INTEGER                                  :: nbonds(naty), ncc
    REAL(real_8)                             :: fx(ncc,ncc)
    LOGICAL                                  :: lprt

    INTEGER                                  :: i, i1, i2, ia, iab, iac, iad, &
                                                ib, ibd, ic, icd, id, ioop, &
                                                iprt, ired(12), j, j1, j2, &
                                                jred, nif
    LOGICAL                                  :: good
    REAL(real_8)                             :: alpha, beta, bred(12), dba, &
                                                dca, fcred, fred, gamma, &
                                                sumlig

    iprt=0
    IF (lprt) THEN
       IF (paral%io_parent)&
            WRITE(iw,'(/,21X,"<<<<< ESTIMATED FORCE CONSTANTS >>>>>",/)')
       IF (paral%io_parent)&
            WRITE(iw,'("   NO."," TYPE",7X,"ATOMS",8X,"QUADRATIC")')
       iprt=1
    ENDIF
    nif=0
    DO i=1,ncc
       DO j=1,ncc
          fx(j,i)=zero
       ENDDO
    ENDDO
    ! :: ::::::::::::::::::
    ! STRETCH COORDINATES:
    ! :: ::::::::::::::::::
    DO ia=2,naty
       DO ib=1,ia-1
          iab=ibind(ib+ia*(ia-1)/2)
          IF (iab.NE.0) THEN
             CALL fstre(fred,fcred,bred,ired,jred,&
                  rcov,ia,ib,nif,ncc,iprt,good)
             IF (good) CALL fxmat(fx,bred,ired,jred,fred,one,ncc)
          ENDIF
       ENDDO
    ENDDO
    ! ..................
    ! LOOP OVER ATOMS A.
    ! ..................
    DO ia=1,naty
       IF (nbonds(ia).GE.2) THEN
          ! .....................................
          ! LOOP OVER ALL ATOMS B CONNECTED TO A.
          ! .....................................
          DO ib=1,naty
             i  =MAX(ia,ib)
             j  =MIN(ia,ib)
             iab=ibind(j+i*(i-1)/2)
             IF (ib.NE.ia.AND.iab.NE.0) THEN
                ! .....................................
                ! LOOP OVER ALL ATOMS C CONNECTED TO A.
                ! .....................................
                DO ic=ib+1,naty
                   i  =MAX(ia,ic)
                   j  =MIN(ia,ic)
                   iac=ibind(j+i*(i-1)/2)
                   IF (ic.NE.ia.AND.iac.NE.0) THEN
                      alpha=winkel(ib,ic,ia)
                      ! :: :::::::::::::
                      ! LINEAR BENDING:
                      ! :: :::::::::::::
                      IF (alpha.EQ.180.0e+00_real_8) THEN
                         ! 
                         DO id=1,naty
                            IF (id.NE.ia.AND.id.NE.ib.AND.id.NE.ic) THEN
                               i1 =MAX(ib,id)
                               j1 =MIN(ib,id)
                               i2 =MAX(ic,id)
                               j2 =MIN(ic,id)
                               ibd=ibind(j1+i1*(i1-1)/2)
                               icd=ibind(j2+i2*(i2-1)/2)
                               IF (ibd.NE.0) THEN
                                  dba=winkel(id,ib,ia)
                                  IF (dba.NE.180.0e+00_real_8) THEN
                                     CALL flinc(fred,fcred,bred,ired,jred,rcov,&
                                          ia,ib,ic,id,nif,ncc,iprt,good)
                                     IF (good) CALL fxmat(fx,bred,ired,jred,&
                                          fred,three,ncc)
                                     CALL flinp(fred,fcred,bred,ired,jred,&
                                          ia,ib,ic,id,nif,ncc,iprt,good)
                                     IF (good) CALL fxmat(fx,bred,ired,jred,&
                                          fred,three,ncc)
                                  ENDIF
                               ENDIF
                               IF (icd.NE.0) THEN
                                  dca=winkel(id,ic,ia)
                                  IF (dca.NE.180.0e+00_real_8) THEN
                                     CALL flinc(fred,fcred,bred,ired,jred,rcov,&
                                          ia,ic,ib,id,nif,ncc,iprt,good)
                                     IF (good) CALL fxmat(fx,bred,ired,jred,&
                                          fred,three,ncc)
                                     CALL flinp(fred,fcred,bred,ired,jred,&
                                          ia,ic,ib,id,nif,ncc,iprt,good)
                                     IF (good) CALL fxmat(fx,bred,ired,jred,&
                                          fred,three,ncc)
                                  ENDIF
                               ENDIF
                            ENDIF
                         ENDDO
                         ! :: :::::::::::::::::
                         ! NON LINEAR BENDING:
                         ! :: :::::::::::::::::
                      ELSE
                         CALL fbend(fred,fcred,bred,ired,jred,rcov,&
                              ia,ib,ic,nif,ncc,iprt,good)
                         IF (good) CALL fxmat(fx,bred,ired,jred,fred,one,ncc)
                      ENDIF
                      ! :: :::::::::::::::::::::::
                      ! OUT OF PLANE COORDINATES:
                      ! :: :::::::::::::::::::::::
                      IF (nbonds(ia).GE.3) THEN
                         ! .....................................
                         ! LOOP OVER ALL ATOMS D CONNECTED TO A.
                         ! .....................................
                         DO id=ic+1,naty
                            i  =MAX(ia,id)
                            j  =MIN(ia,id)
                            iad=ibind(j+i*(i-1)/2)
                            IF (id.NE.ia.AND.iad.NE.0) THEN
                               alpha=winkel(ib,ic,ia)
                               beta =winkel(ic,id,ia)
                               gamma=winkel(ib,id,ia)
                               IF (alpha.NE.180.0e+00_real_8.AND.beta.NE.180.0e+00_real_8 &
                                    .AND.gamma.NE.180.0e+00_real_8) THEN
                                  DO ioop=1,3
                                     CALL foutp(fred,fcred,bred,ired,jred,rcov,&
                                          ia,ib,ic,id,ioop,&
                                          nif,ncc,iprt,good)
                                     IF (good) CALL fxmat(fx,bred,ired,jred,&
                                          fred,three,ncc)
                                  ENDDO
                               ENDIF
                               ! ......................
                               ! END LOOP OVER ATOMS D.
                               ! ......................
                            ENDIF
                         ENDDO
                      ENDIF
                      ! ......................
                      ! END LOOP OVER ATOMS C.
                      ! ......................
                   ENDIF
                ENDDO
                ! :: ::::::::::::::::::
                ! TORSION COORDINATES:
                ! :: ::::::::::::::::::
                IF (naty.GE.4.AND.nbonds(ib).GE.2.AND.ia.LT.ib) THEN
                   ! .....................................
                   ! LOOP OVER ALL ATOMS C CONNECTED TO A.
                   ! .....................................
                   DO ic=1,naty
                      i  =MAX(ic,ia)
                      j  =MIN(ic,ia)
                      iac=ibind(j+i*(i-1)/2)
                      IF (ic.NE.ia.AND.ic.NE.ib.AND.iac.NE.0) THEN
                         ! .....................................
                         ! LOOP OVER ALL ATOMS D CONNECTED TO A.
                         ! .....................................
                         DO id=1,naty
                            i  =MAX(ib,id)
                            j  =MIN(ib,id)
                            ibd=ibind(j+i*(i-1)/2)
                            IF (id.NE.ia.AND.id.NE.ib.AND.id.NE.ic.AND.&
                                 ibd.NE.0) THEN
                               alpha=winkel(ic,ib,ia)
                               beta =winkel(ia,id,ib)
                               IF (alpha.NE.180.0e+00_real_8.AND.&
                                    beta.NE.180.0e+00_real_8) THEN
                                  sumlig=REAL(nbonds(ia)+nbonds(ib)-2,kind=real_8)
                                  CALL ftors(fred,fcred,bred,ired,jred,&
                                       rcov,sumlig,ia,ib,ic,id,&
                                       nif,ncc,iprt,good)
                                  IF (good) CALL fxmat(fx,bred,ired,jred,&
                                       fred,one,ncc)
                               ENDIF
                               ! ......................
                               ! END LOOP OVER ATOMS D.
                               ! ......................
                            ENDIF
                         ENDDO
                         ! ......................
                         ! END LOOP OVER ATOMS C.
                         ! ......................
                      ENDIF
                   ENDDO
                ENDIF
                ! ......................
                ! END LOOP OVER ATOMS B.
                ! ......................
             ENDIF
          ENDDO
          ! ......................
          ! END LOOP OVER ATOMS A.
          ! ......................
       ENDIF
    ENDDO
    ! Complete the matrix FX (half filling)
    DO i=1,ncc
       DO j=i+1,ncc
          fx(i,j)=fx(i,j)+fx(j,i)
          fx(j,i)=fx(i,j)
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE fstart
  ! ==================================================================
  SUBROUTINE fbend(fred,fcred,bred,ired,jred,&
       rcov,ia,ib,ic,nif,ncc,iprt,good)
    ! ==--------------------------------------------------------------==
    ! == BENDING FORCE CONSTANTS.                                     ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: fred, fcred, bred(12)
    INTEGER                                  :: ired(12), jred
    REAL(real_8)                             :: rcov(0:104)
    INTEGER                                  :: ia, ib, ic, nif, ncc, iprt
    LOGICAL                                  :: good

    INTEGER                                  :: iza, izb, izc, minz
    REAL(real_8)                             :: abend(2), aharm, bharm, bset, &
                                                charm, delta, dharm, prod, &
                                                qred

    CALL bandq(bred,ired,jred,qred,'BEND',ib,ic,ia,0,ncc,good)
    IF (.NOT.good) RETURN
    aharm=-0.089e+00_real_8
    bharm= 0.11e+00_real_8
    charm= 0.44e+00_real_8
    dharm=-0.42e+00_real_8
    iza  =NINT(zan(ia))
    izb  =NINT(zan(ib))
    izc  =NINT(zan(ic))
    delta=SQRT(dist(ia,ib))-rcov(iza)-rcov(izb)+&
         SQRT(dist(ia,ic))-rcov(iza)-rcov(izc)
    prod =(rcov(iza)+rcov(izb))*(rcov(iza)+rcov(izc))
    fred =aharm+bharm*EXP(-charm*delta) / (prod)**dharm
    ! .......................................
    ! SCHLEGEL PARAMETRIZATION FOR HARMONICS.
    ! .......................................
    IF (cnti%npara.GT.0) THEN
       bset     = one
       ! BSET     = 1.3_real_8
       abend(1) = 0.130e+00_real_8*bset
       abend(2) = 0.112e+00_real_8*bset
       ! 
       iza =NINT(zan(ia))
       izb =NINT(zan(ib))
       izc =NINT(zan(ic))
       minz=MIN(izb,izc)
       IF (minz.EQ.1) THEN
          IF (cnti%npara.EQ.1) THEN
             fred=0.160_real_8
          ELSE
             fred=abend(1)*((rcov(iza)+rcov(izb))*&
                  (rcov(iza)+rcov(izc)))**2/&
                  (SQRT(dist(ia,ib))+SQRT(dist(ia,ic)))**2
          ENDIF
       ELSE
          IF (cnti%npara.EQ.1) THEN
             fred=0.250_real_8
          ELSE
             fred=abend(2)*((rcov(iza)+rcov(izb))*&
                  (rcov(iza)+rcov(izc)))**2/&
                  (SQRT(dist(ia,ib))+SQRT(dist(ia,ic)))**2
          ENDIF
       ENDIF
    ENDIF
    ! 
    nif=nif+1
    IF ((iprt.EQ.1).AND.paral%io_parent)&
         WRITE(iw,'(I5,".",A5,4I4,2F14.8)')&
         nif,'BEND',ib,ic,ia,0,fred
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE fbend
  ! ==================================================================
  SUBROUTINE flinc(fred,fcred,bred,ired,jred,&
       rcov,ia,ib,ic,id,nif,ncc,iprt,good)
    ! ==--------------------------------------------------------------==
    ! == COLLINEAR BENDING CONSTANTS.                                 ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: fred, fcred, bred(12)
    INTEGER                                  :: ired(12), jred
    REAL(real_8)                             :: rcov(0:104)
    INTEGER                                  :: ia, ib, ic, id, nif, ncc, iprt
    LOGICAL                                  :: good

    INTEGER                                  :: iza, izb, izc, izd
    REAL(real_8)                             :: aharm, bharm, charm, delta, &
                                                eharm, prod, qred

    CALL bandq(bred,ired,jred,qred,'LINC',ic,ib,ia,id,ncc,good)
    IF (.NOT.good) RETURN
    aharm= 0.0025e+00_real_8
    bharm= 0.0061e+00_real_8
    charm= 3.0e+00_real_8
    eharm= 0.80e+00_real_8
    iza =NINT(zan(ia))
    izb =NINT(zan(ib))
    izc =NINT(zan(ic))
    izd =NINT(zan(id))
    delta=SQRT(dist(ia,ic))-rcov(iza)-rcov(izc)
    prod =(rcov(iza)+rcov(izb))*(rcov(iza)+rcov(izd))
    fred =aharm+bharm * prod**eharm * EXP(-charm*delta)
    ! ==--------------------------------------------------------------==
    nif=nif+1
    IF ((iprt.EQ.1).AND.paral%io_parent)&
         WRITE(iw,'(I5,".",A5,4I4,2F14.8)')&
         nif,'LINC',ic,ib,ia,id,fred
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE flinc
  ! ==================================================================
  SUBROUTINE flinp(fred,fcred,bred,ired,jred,&
       ia,ib,ic,id,nif,ncc,iprt,good)
    ! ==--------------------------------------------------------------==
    ! == PERPENDICULAR LINEAR BENDING CONSTANTS.                      ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: fred, fcred, bred(12)
    INTEGER                                  :: ired(12), jred, ia, ib, ic, &
                                                id, nif, ncc, iprt
    LOGICAL                                  :: good

    REAL(real_8)                             :: qred

    CALL bandq(bred,ired,jred,qred,'LINP',ic,ib,ia,id,ncc,good)
    IF (.NOT.good) RETURN
    nif=nif+1
    IF ((iprt.EQ.1).AND.paral%io_parent)&
         WRITE(iw,'(I5,".",A5,4I4,2F14.8)')&
         nif,'LINP',ic,ib,ia,id,fred
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE flinp
  ! ==================================================================
  SUBROUTINE foutp(fred,fcred,bred,ired,jred,&
       rcov,ia,ib,ic,id,ioop,nif,ncc,iprt,good)
    ! ==--------------------------------------------------------------==
    ! == OUT OF PLANE FORCE CONSTANTS.                                ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: fred, fcred, bred(12)
    INTEGER                                  :: ired(12), jred
    REAL(real_8)                             :: rcov(0:104)
    INTEGER                                  :: ia, ib, ic, id, ioop, nif, &
                                                ncc, iprt
    LOGICAL                                  :: good

    INTEGER                                  :: i1, i2, i3, iz1, iz2, iz3, iza
    REAL(real_8)                             :: aharm, aoutp, bharm, bset, &
                                                charm, delta, dharm, eharm, &
                                                pi2, prod, qred

    IF (ioop.EQ.1) THEN
       i1=id
       i2=ib
       i3=ic
    ELSE IF (ioop.EQ.2) THEN
       i1=ib
       i2=ic
       i3=id
    ELSE
       i1=ic
       i2=ib
       i3=id
    ENDIF
    CALL bandq(bred,ired,jred,qred,'OUTP',i1,i2,i3,ia,ncc,good)
    IF (.NOT.good) RETURN
    aharm= 0.0025e+00_real_8
    bharm= 0.0061e+00_real_8
    charm= 3.0e+00_real_8
    dharm= 4.0e+00_real_8
    eharm= 0.80e+00_real_8
    iza =NINT(zan(ia))
    iz1 =NINT(zan(i1))
    iz2 =NINT(zan(i2))
    iz3 =NINT(zan(i3))
    delta=SQRT(dist(ia,i1))-rcov(iza)-rcov(iz1)
    prod =(rcov(iza)+rcov(iz2))*(rcov(iza)+rcov(iz3))
    fred =aharm+bharm * prod**eharm * (COS(qred))**dharm *&
         EXP(-charm*delta)
    ! .......................................
    ! SCHLEGEL PARAMETERZATION FOR HARMONICS.
    ! .......................................
    IF (cnti%npara.GT.0) THEN
       pi2=pi/two
       bset     = one
       ! BSET     = 1.3_real_8
       aoutp = 0.045e+00_real_8 *bset
       fred=aoutp/pi2**4*(pi2-ABS(qred))**4
    ENDIF
    ! ==--------------------------------------------------------------==
    nif=nif+1
    IF ((iprt.EQ.1).AND.paral%io_parent)&
         WRITE(iw,'(I5,".",A5,4I4,2F14.8)')&
         nif,'OUTP',i1,i2,i3,ia,fred
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE foutp
  ! ==================================================================
  SUBROUTINE fstre(fred,fcred,bred,ired,jred,&
       rcov,ia,ib,nif,ncc,iprt,good)
    ! ==--------------------------------------------------------------==
    ! == STRETCHING FORCE CONSTANTS.                                  ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: fred, fcred, bred(12)
    INTEGER                                  :: ired(12), jred
    REAL(real_8)                             :: rcov(0:104)
    INTEGER                                  :: ia, ib, nif, ncc, iprt
    LOGICAL                                  :: good

    INTEGER                                  :: iper, iza, izb, jper, maxz, &
                                                minz
    REAL(real_8)                             :: aharm, astre, bharm, bset, &
                                                bstre(3,3), delta, qred

    CALL bandq(bred,ired,jred,qred,'STRE',ia,ib,0,0,ncc,good)
    IF (.NOT.good) RETURN
    aharm=0.3601e+00_real_8
    bharm=-1.944e+00_real_8
    iza  =NINT(zan(ia))
    izb  =NINT(zan(ib))
    delta=qred-rcov(iza)-rcov(izb)
    fred =aharm*EXP(bharm*delta)
    ! .......................................
    ! SCHLEGEL PARAMETERZATION FOR HARMONICS.
    ! .......................................
    IF (cnti%npara.GT.0) THEN
       bset      = one
       ! BSET      = 1.3_real_8
       astre     = 1.734e+00_real_8*bset
       bstre(1,1)=-0.244e+00_real_8
       bstre(1,2)= 0.352e+00_real_8
       bstre(2,2)= 1.085e+00_real_8
       bstre(1,3)= 0.660e+00_real_8
       bstre(2,3)= 1.522e+00_real_8
       bstre(3,3)= 2.068e+00_real_8
       ! 
       minz=NINT(MIN(zan(ia),zan(ib)))
       maxz=NINT(MAX(zan(ia),zan(ib)))
       iper=1
       jper=1
       IF (minz.GT.2)  iper=2
       IF (minz.GT.10) iper=3
       IF (maxz.GT.2)  jper=2
       IF (maxz.GT.10) jper=3
       fred=astre/(qred-bstre(iper,jper))**3
    ENDIF
    nif=nif+1
    IF ((iprt.EQ.1).AND.paral%io_parent)&
         WRITE(iw,'(I5,".",A5,4I4,2F14.8)')&
         nif,'STRE',ia,ib,0,0,fred
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE fstre
  ! ==================================================================
  SUBROUTINE ftors(fred,fcred,bred,ired,jred,&
       rcov,sumlig,ia,ib,ic,id,nif,ncc,iprt,good)
    ! ==--------------------------------------------------------------==
    ! == Torsion Force Constants calculation                          ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: fred, fcred, bred(12)
    INTEGER                                  :: ired(12), jred
    REAL(real_8)                             :: rcov(0:104), sumlig
    INTEGER                                  :: ia, ib, ic, id, nif, ncc, iprt
    LOGICAL                                  :: good

    INTEGER                                  :: iza, izb
    REAL(real_8)                             :: aharm, ators, bharm, bset, &
                                                btors, charm, delta, dharm, &
                                                eharm, prod, qred

    CALL bandq(bred,ired,jred,qred,'TORS',ic,ia,ib,id,ncc,good)
    IF (.NOT.good) RETURN
    aharm= 0.0015e+00_real_8
    bharm= 14.0e+00_real_8
    charm= 2.85e+00_real_8
    dharm= 0.57e+00_real_8
    eharm= 4.0e+00_real_8
    iza =NINT(zan(ia))
    izb =NINT(zan(ib))
    delta=SQRT(dist(ia,ib))-rcov(iza)-rcov(izb)
    prod =(rcov(iza)+rcov(izb))*SQRT(dist(ia,ib))
    fred =aharm+bharm*sumlig**dharm/prod**eharm*EXP(-charm*delta)
    ! .......................................
    ! SCHLEGEL PARAMETRIZATION FOR HARMONICS.
    ! .......................................
    IF (cnti%npara.GT.0) THEN
       bset      = one
       ! BSET      = 1.3_real_8
       ators     = 0.0059e+00_real_8*bset
       btors     = 4.348e+00_real_8 *bset
       iza =NINT(zan(ia))
       izb =NINT(zan(ib))
       IF (cnti%npara.EQ.1) THEN
          fred=0.0023_real_8-0.07_real_8*(SQRT(dist(ia,ib))-rcov(iza)-rcov(izb))
          fred=ABS(fred)
       ELSE
          fred=ators*EXP(-btors*(SQRT(dist(ia,ib))&
               -rcov(iza)-rcov(izb)))
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    nif=nif+1
    IF ((iprt.EQ.1).AND.paral%io_parent)&
         WRITE(iw,'(I5,".",A5,4I4,2F14.8)')&
         nif,'TORS',ic,ia,ib,id,fred
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE ftors
  ! ==================================================================
  SUBROUTINE bandq(bred,ired,jred,qred,itypf,ia,ib,ic,id,ncc,good)
    ! ==--------------------------------------------------------------==
    ! == INPUT:                                                       ==
    ! ==   IA,IB,IC,ID  NUMBERS OF 4 ATOMS                            ==
    ! ==   NCC          NUMBER OF FREEDOM DEGREES                     ==
    ! ==   ITYPF        TYPE OF CALCULATION                           ==
    ! == OUTPUT:                                                      ==
    ! ==   JRED         Number of factors                             ==
    ! ==   BRED(1:JRED) Factors for force constant                    ==
    ! ==   IRED(1:JRED) Index in (1:NCC) for factors                  ==
    ! ==   QRED         Force constant                                ==
    ! ==   GOOD         .TRUE. if the calculation is possible         ==
    ! ==                .FALSE. otherwise                             ==
    ! ==--------------------------------------------------------------==
    ! == Conventions used in this routine:                            ==
    ! == UU(1:3) or UV(1:3)                                           ==
    ! == VV(1:3) or UV(4:6)                                           ==
    ! == WW(1:3) or UV(7:9)                                           ==
    ! == ZZ(1:3) or UV(10:12)                                         ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: bred(12)
    INTEGER                                  :: ired(12), jred
    REAL(real_8)                             :: qred
    CHARACTER(len=4)                         :: itypf
    INTEGER                                  :: ia, ib, ic, id, ncc
    LOGICAL                                  :: good

    INTEGER, PARAMETER                       :: ih = 4 
    REAL(real_8), PARAMETER                  :: pih = pi/two 

    CHARACTER(len=4), DIMENSION(7) :: &
      itipus = (/'STRE','INVR','BEND','OUTP','TORS','LINC','LINP'/)
    INTEGER                                  :: iaf(ih), j, j1, l, l1, m, m1
    REAL(real_8) :: cfi1, cfi2, cfi3, co, co2, cteta, cx, den, r1, r2, r3, &
      rm1, rm2, s, sfi1, si, si2, si3, st2, st3, steta, u(3), uu(3), uv(12), &
      v(3), vv(3), w(3), ww(3), x(3), z(3), zz(3)

    EQUIVALENCE (uv(1),uu(1)),(uv(4),vv(1)),(uv(7),ww(1)),(uv(10),zz(1))
    ! ==--------------------------------------------------------------==
    good=.TRUE.
    DO l=1,3
       uu(l)=zero
    ENDDO
    DO l=1,3
       vv(l)=zero
    ENDDO
    DO l=1,3
       ww(l)=zero
    ENDDO
    DO l=1,3
       zz(l)=zero
    ENDDO
    qred=zero
    iaf(1)=ia
    iaf(2)=ib
    iaf(3)=ic
    iaf(4)=id
    ! 
    IF (itypf.EQ.itipus(1)) THEN
       ! ...........
       ! STRETCHING.
       ! ...........
       CALL vektor(uu,r1,iaf(ih-3),iaf(ih-2))
       vv(1)=-uu(1)
       vv(2)=-uu(2)
       vv(3)=-uu(3)
       qred=qred+r1
    ELSEIF (itypf.EQ.itipus(2)) THEN
       ! ........
       ! INVERSE.
       ! ........
       ! ==------------------------------------------------------------==
       ! == For Inverse:                                               ==
       ! ==     DIST(IA,IB) /= 0                                       ==
       ! ==------------------------------------------------------------==
       CALL vektor(uu,r1,iaf(ih-3),iaf(ih-2))
       IF (r1.EQ.zero) THEN
          good=.FALSE.
          RETURN
       ENDIF
       rm1  =one/r1
       rm2  =rm1*rm1
       uu(1)=-uu(1)*rm2
       uu(2)=-uu(2)*rm2
       uu(3)=-uu(3)*rm2
       vv(1)=-uu(1)
       vv(2)=-uu(2)
       vv(3)=-uu(3)
       qred=qred+rm1
    ELSEIF (itypf.EQ.itipus(3)) THEN
       ! ........
       ! BENDING.
       ! ........
       ! ==------------------------------------------------------------==
       ! == For Bending:                                               ==
       ! ==     R1=dist(IA,IC) /= 0 and R2=dist(IB,IC) /= 0            ==
       ! ==     SI=(IA,IC).(IB,IC) /= 1                                ==
       ! ==        i.e. IC,IA,IB not aligned                           ==
       ! ==------------------------------------------------------------==
       CALL vektor(u,r1,iaf(ih-3),iaf(ih-1))
       CALL vektor(v,r2,iaf(ih-2),iaf(ih-1))
       IF (r1.EQ.zero.OR.r2.EQ.zero) THEN
          good=.FALSE.
          RETURN
       ENDIF
       co=scalaru(u,v)
       si=SQRT(one-co*co)
       IF (si.EQ.zero) THEN
          good=.FALSE.
          RETURN
       ENDIF
       DO l=1,3
          uu(l)=(co*u(l)-v(l))/(si*r1)
          vv(l)=(co*v(l)-u(l))/(si*r2)
          ww(l)=-uu(l)-vv(l)
       ENDDO
       qred=qred+ACOS(co)
    ELSEIF (itypf.EQ.itipus(4)) THEN
       ! .............
       ! OUT OF PLANE.
       ! .............
       ! ==------------------------------------------------------------==
       ! == For Out of Plane:                                          ==
       ! ==     R1=DIST(IA,IC) /= 0 AND R2=dist(IB,ID) /= 0 and        ==
       ! ==     R3=dist(IC,ID) /= 0                                    ==
       ! ==     DEN=CTETA*SFI1*SFI1 /= 0                               ==
       ! ==        CTETA=((IA,ID).normal(IB,IC,ID)) /= 1               ==
       ! ==             i.e. (IA,ID) orthogonal to (IB,IC,ID) plane    ==
       ! ==        SFI1=(IB,ID).(IC,ID) /= 1                           ==
       ! ==             i.e. ID,IB,IC not aligned                      ==
       ! ==------------------------------------------------------------==
       CALL vektor(u,r1,iaf(ih-3),iaf(ih))
       CALL vektor(v,r2,iaf(ih-2),iaf(ih))
       CALL vektor(w,r3,iaf(ih-1),iaf(ih))
       IF (r1.EQ.0.OR.r2.EQ.zero.OR.r3.EQ.zero) THEN
          good=.FALSE.
          RETURN
       ENDIF
       CALL normal(v,w,z)
       steta=scalaru(u,z)
       cteta=SQRT(one-steta*steta)
       cfi1 =scalaru(v,w)
       sfi1 =SQRT(one-cfi1*cfi1)
       cfi2 =scalaru(w,u)
       cfi3 =scalaru(v,u)
       den  =cteta*sfi1*sfi1
       IF (den.EQ.zero) THEN
          good=.FALSE.
          RETURN
       ENDIF
       st2  =(cfi1*cfi2-cfi3)/(r2*den)
       st3  =(cfi1*cfi3-cfi2)/(r3*den)
       DO l=1,3
          vv(l)=z(l)*st2
          ww(l)=z(l)*st3
       ENDDO
       CALL normal(z,u,x)
       CALL normal(u,x,z)
       DO l=1,3
          uu(l)=z(l)/r1
          zz(l)=-uu(l)-vv(l)-ww(l)
       ENDDO
       cx=-one
       IF (steta.LT.zero) cx=one
       qred=qred-cx*ACOS(cteta)
    ELSEIF (itypf.EQ.itipus(5)) THEN
       ! ........
       ! TORSION.
       ! ........
       ! ==------------------------------------------------------------==
       ! == For Torsion:                                               ==
       ! ==     dist(IA,IB) /= 0, dist(IB,IC) /=0 and dist(IC,ID) /= 0 ==
       ! ==     The angle ((IA,IB),(IC,IB)) MUST BE /=   0 or 180      ==
       ! ==     The angle ((IC,IB),(IC,ID)) MUST BE /=   0 or 180      ==
       ! ==------------------------------------------------------------==
       CALL vektor(u,r1,iaf(ih-3),iaf(ih-2))
       CALL vektor(v,r2,iaf(ih-1),iaf(ih-2))
       CALL vektor(w,r3,iaf(ih-1),iaf(ih)  )
       IF (r1.EQ.zero.OR.r2.EQ.zero.OR.r3.EQ.zero) THEN
          good=.FALSE.
          RETURN
       ENDIF
       CALL normal(u,v,z)
       CALL normal(w,v,x)
       co =scalaru(u,v)
       co2=scalaru(v,w)
       si =SQRT(one-co*co)
       si2=SQRT(one-co2*co2)
       IF (si.EQ.zero.OR.si2.EQ.zero) THEN
          good=.FALSE.
          RETURN
       ENDIF
       DO l=1,3
          uu(l)=z(l)/(r1*si)
          zz(l)=x(l)/(r3*si2)
          vv(l)=(r1*co/r2-one)*uu(l)-r3*co2/r2*zz(l)
          ww(l)=-uu(l)-vv(l)-zz(l)
       ENDDO
       co=scalaru(z,x)
       u(1)=z(2)*x(3)-z(3)*x(2)
       u(2)=z(3)*x(1)-z(1)*x(3)
       u(3)=z(1)*x(2)-z(2)*x(1)
       si3 =SQRT(u(1)**2+u(2)**2+u(3)**2)
       co2 =scalaru(u,v)
       s   =arc1(-co,si3)
       IF (co2.LT.zero) s=-s
       ! ..............................................
       ! THE RANGE OF TORSION ANGLE IS -PI/2 TO 3*PI/2.
       ! ..............................................
       IF (s.GT.pih) s=s-two*pi
       qred=qred-s
    ELSEIF (itypf.EQ.itipus(6)) THEN
       ! ........................
       ! LINEAR COPLANAR BENDING.
       ! ........................
       ! ==------------------------------------------------------------==
       ! == For Linear coplanar Bending:                               ==
       ! ==     dist(IA,IC) /= 0 and dist(IB,IC) /= 0                  ==
       ! ==------------------------------------------------------------==
       CALL vektor(u,r1,iaf(ih-3),iaf(ih-1))
       CALL vektor(v,r2,iaf(ih)  ,iaf(ih-1))
       CALL vektor(x,r2,iaf(ih-2),iaf(ih-1))
       IF (r1.EQ.0.OR.r2.EQ.zero) THEN
          good=.FALSE.
          RETURN
       ENDIF
       co   =scalaru(v,u)
       co2  =scalaru(x,v)
       qred=qred+(pi-ACOS(co)-ACOS(co2))
       CALL normal(v,u,w)
       CALL normal(u,w,z)
       CALL normal(x,v,w)
       CALL normal(w,x,u)
       ! ............................................................
       ! INTERNAL COORDINATE POSITIVE IF ATOM A MOVES TOWARDS ATOM D.
       ! ............................................................
       DO l=1,3
          uu(l)=z(l)/r1
          vv(l)=u(l)/r2
          ww(l)=-uu(l)-vv(l)
       ENDDO
    ELSEIF (itypf.EQ.itipus(7)) THEN
       ! .............................
       ! LINEAR PERPENDICULAR BENDING.
       ! .............................
       ! ==------------------------------------------------------------==
       ! == For Linear Perpendicular Bending:                          ==
       ! ==     dist(IA,IC) /= 0 and dist(IC,ID) /= 0                  ==
       ! ==------------------------------------------------------------==
       CALL vektor(u,r1,iaf(ih-3),iaf(ih-1))
       CALL vektor(v,r2,iaf(ih)  ,iaf(ih-1))
       CALL vektor(z,r2,iaf(ih-2),iaf(ih-1))
       IF (r1.EQ.zero.OR.r2.EQ.zero) THEN
          good=.FALSE.
          RETURN
       ENDIF
       CALL normal(v,u,w)
       CALL normal(z,v,x)
       DO l=1,3
          uu(l)=w(l)/r1
          vv(l)=x(l)/r2
          ww(l)=-uu(l)-vv(l)
       ENDDO
       co   =scalaru(u,w)
       co2  =scalaru(z,w)
       qred=qred+(pi-ACOS(co)-ACOS(co2))
    ELSE
       CALL stopgm('BANDQ','THIS OPTION DOES NOT EXIST',& 
            __LINE__,__FILE__)
    ENDIF
    ! Calculate BRED(1:JRED)
    jred=0
    ! For each atom (given by IAF)
    DO j=1,4
       m=iaf(ih-4+j)
       IF (m.GT.0) THEN
          m =3*(m-1)
          j1=3*(j-1)
          ! For each coordinates of atom
          DO l=1,3
             jred=jred+1
             m1      =m+l
             ired(jred)=m1
             l1      =j1+l
             bred(jred)=uv(l1)
          ENDDO
       ENDIF
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE bandq
  ! ==================================================================
  FUNCTION arc1 (x,y)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: x, y, arc1

    REAL(real_8)                             :: pih, s

    pih=pi/two
    IF (ABS(x).LT.1.0e-11_real_8) THEN
       arc1=pih
    ELSE
       s=ATAN(y/x)
       IF (x.LT.zero) s=s+pi
       arc1=s
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION arc1
  ! ==================================================================
  FUNCTION  dist(nuc1,nuc2)
    ! ==--------------------------------------------------------------==
    ! == W A R N I N G:                                               ==
    ! ==             THIS ROUTINE COMPUTES THE SQUARE OF THE DISTANCE ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nuc1, nuc2
    REAL(real_8)                             :: dist

    REAL(real_8)                             :: c1, c2, c3,c1_, c2_, c3_

    c1_=c(1,nuc1)-c(1,nuc2)
    c2_=c(2,nuc1)-c(2,nuc2)
    c3_=c(3,nuc1)-c(3,nuc2)
    CALL pbc(c1_,c2_,c3_,c1,c2,c3,1,parm%apbc,parm%ibrav)
    dist = c1*c1 + c2*c2 + c3*c3
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION dist
  ! ==================================================================
  SUBROUTINE nom(u)
    ! ==--------------------------------------------------------------==
    ! == Normalize the vector U(1:3)                                  ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: u(3)

    INTEGER                                  :: i
    REAL(real_8)                             :: d, x

    d=SQRT(u(1)*u(1)+u(2)*u(2)+u(3)*u(3))
    IF (d.LE.zero) RETURN
    x=one/d
    DO i=1,3
       u(i)=u(i)*x
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE nom
  ! ==================================================================
  SUBROUTINE normal(u,v,w)
    ! ==--------------------------------------------------------------==
    ! == Normalized vector product W(1:3) of U(1:3) and V(1:3)        ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: u(3), v(3), w(3)

    REAL(real_8), PARAMETER                  :: tol = 1.e-20_real_8 

! ==--------------------------------------------------------------==

    w(1)=u(2)*v(3)-u(3)*v(2)
    w(2)=u(3)*v(1)-u(1)*v(3)
    w(3)=u(1)*v(2)-u(2)*v(1)
    IF (ABS(w(1)).LT.tol.AND.ABS(w(2)).LT.tol.AND.ABS(w(3)).LT.tol)&
         RETURN
    CALL nom(w)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE normal
  ! ==================================================================
  FUNCTION scalaru(u,v)
    ! ==--------------------------------------------------------------==
    ! == Scalar product of unit vectors                               ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: u(3), v(3), scalaru

    INTEGER                                  :: i

    scalaru=zero
    DO i=1,3
       scalaru=scalaru+u(i)*v(i)
    ENDDO
    ! Some computer has a small error (DEC)
    IF (scalaru.GT.one) THEN
       scalaru=one
    ELSEIF (scalaru.LT.-one) THEN
       scalaru=-one
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION scalaru
  ! ==================================================================
  SUBROUTINE vektor(u,r,i,j)
    ! ==--------------------------------------------------------------==
    ! == U=unit vector (I,J) (V_I - V_J) and R=dist(I,J)              ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: u(3), r, u_(3)
    INTEGER                                  :: i, j

! ==--------------------------------------------------------------==

    u_(1)=c(1,i)-c(1,j)
    u_(2)=c(2,i)-c(2,j)
    u_(3)=c(3,i)-c(3,j)
    CALL pbc(u_(1),u_(2),u_(3),u(1),u(2),u(3),1,parm%apbc,parm%ibrav)
    r=SQRT(u(1)*u(1)+u(2)*u(2)+u(3)*u(3))
    CALL nom(u)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE vektor
  ! ==================================================================
  FUNCTION winkel(ia,ib,ic)
    ! ==--------------------------------------------------------------==
    ! == CALCULATION OF THE ANGLE: A                                  ==
    ! == (IN DEGREES)               \                                 ==
    ! ==                             C - B                            ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ia, ib, ic
    REAL(real_8)                             :: winkel

    REAL(real_8)                             :: cosin, r1, r2, u(3), v(3)

    CALL vektor(u,r1,ic,ia)
    CALL vektor(v,r2,ic,ib)
    cosin=scalaru(u,v)
    IF (ABS(cosin+one).LE.1.e-6_real_8) THEN
       winkel=180.0e+00_real_8
    ELSE
       winkel=180.0e+00_real_8*dacos(cosin)/pi
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION winkel
  ! ==================================================================
  SUBROUTINE fxmat(fx,bred,ired,jred,fred,xnrm,ncc)
    ! ==--------------------------------------------------------------==
    ! == Fill the force constant matrix FX                            ==
    ! == INPUT:                                                       ==
    ! ==    NCC         Number of freedom degrees                     ==
    ! ==    BRED(1:NCC)                                               ==
    ! ==    FRED        Force constant                                ==
    ! ==    XNRM        Normalisation number                          ==
    ! == OUTPUT:                                                      ==
    ! ==    FX(1:NCC,1:NCC) Force constant matrix                     ==
    ! ==      W A R N I N G: Fill only half FX                        ==
    ! ==--------------------------------------------------------------==
    ! == Condition:XNRM /= 0                                          ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: bred(12)
    INTEGER                                  :: ired(12), jred
    REAL(real_8)                             :: fred, xnrm
    INTEGER                                  :: ncc
    REAL(real_8)                             :: fx(ncc,ncc)

    INTEGER                                  :: i, ii, j, jj

    IF (fred.EQ.0._real_8) RETURN
    DO i=1,jred
       IF (bred(i).NE.0._real_8) THEN
          ii=ired(i)
          DO j=i,jred
             jj=ired(j)
             fx(ii,jj)=fx(ii,jj)+bred(i)*bred(j)*fred/xnrm
          ENDDO
       ENDIF
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE fxmat
  ! ==================================================================


END MODULE fstart_utils
