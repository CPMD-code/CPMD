MODULE numpw_utils
  USE cell,                            ONLY: cell_com
  USE cnst,                            ONLY: pi
  USE error_handling,                  ONLY: stopgm
  USE fftchk_utils,                    ONLY: fftchk
  USE fint,                            ONLY: fint1
  USE gvec,                            ONLY: gvec_com,&
                                             nhgwish,&
                                             tnhgwish
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE parac,                           ONLY: paral
  USE rggen_utils,                     ONLY: recips
  USE sort_utils,                      ONLY: sort2
  USE sphe,                            ONLY: gcutka,&
                                             gcutwmax,&
                                             gcutwmin,&
                                             tsphere
  USE system,                          ONLY: cntr,&
                                             dual00,&
                                             ncpw,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: numpw
  PUBLIC :: gkpwf

CONTAINS

  ! ==================================================================
  SUBROUTINE numpw
    ! ==--------------------------------------------------------------==
    ! == CALCULATES THE NUMBER OF UNIQUE PW COMPONENTS FOR THE DENSITY==
    ! ==--------------------------------------------------------------==
    ! Variables
    CHARACTER(*), PARAMETER                  :: procedureN = 'numpw'
    REAL(real_8), PARAMETER                  :: delta = 1.e-6_real_8 

    INTEGER                                  :: i, ierr, ig, iloop, ir, isub, &
                                                j, jmax, jmin, k, kmax, kmin, &
                                                nr1f, nr1m1, nr2f, nr2m1, &
                                                nr3f, nr3m1, numtab
    INTEGER, ALLOCATABLE                     :: ntab(:)
    REAL(real_8)                             :: energy, g2, ghigh, glow, t
    REAL(real_8), ALLOCATABLE                :: tabg2(:)

! ==--------------------------------------------------------------==

    CALL tiset('     NUMPW',isub)
    ! Calculate the reciprocal lattice vectors B1,B2 and B3
    CALL recips(parm%alat,parm%a1,parm%a2,parm%a3,gvec_com%b1,gvec_com%b2,gvec_com%b3)
    ! Save NR1, NR2, and NR3 for MESH option.
    nr1f=parm%nr1
    nr2f=parm%nr2
    nr3f=parm%nr3
    ! Maximum K point amplitude (Kmax) in Brillouin zone.
    gcutka=gkpwf(gvec_com%b1,gvec_com%b2,gvec_com%b3,tkpts%tkpnt.AND.tsphere)
    IF (gcutka.EQ.0._real_8) THEN
       gcutwmin=gvec_com%gcutw
       gcutwmax=gvec_com%gcutw
    ELSE
       ! For each k and each G, |k+G|^2 < GCUTW
       ! so we define GCUTWMIN and GCUTWMAX as:
       ! GCUTWMIN < |G|^2 < GCUTWMAX we have to apply a mask.
       gcutwmin = MAX(0._real_8,SQRT(gvec_com%gcutw)-SQRT(gcutka))**2
       gcutwmax = (SQRT(gvec_com%gcutw)+SQRT(gcutka))**2
    ENDIF
    ! GCUT has to be bigger than GCUTWMAX
    energy=MAX(gvec_com%gcut,gcutwmax)
    ! Determination of NR1, NR2, and NR3
    CALL change_gcut(energy)
30  CONTINUE
    nr1m1=parm%nr1-1
    nr2m1=parm%nr2-1
    nr3m1=parm%nr3-1
    ! To have the pw density number required
    IF (tnhgwish) THEN
       numtab=INT(nr1m1*nr2m1*nr3m1/2)
       ALLOCATE(tabg2(numtab),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(ntab(numtab),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! Construction of g-vectors
    ig=0
    glow=0.0_real_8
    ghigh=gcutwmax
    DO iloop=1,2
       IF (iloop.EQ.1) THEN
          glow=0.0_real_8
          ghigh=gcutwmax
       ELSEIF (iloop.EQ.2) THEN
          glow=gcutwmax
          ghigh=gvec_com%gcut
       ENDIF
       DO i=0,parm%nr1-1
          jmin=-parm%nr2+1
          jmax=parm%nr2-1
          IF (i.EQ.0) jmin=0
          DO j=jmin,jmax
             kmin=-parm%nr3+1
             kmax=parm%nr3-1
             IF (i.EQ.0.AND.j.EQ.0) kmin=0
             DO k=kmin,kmax
                g2=0._real_8
                DO ir=1,3
                   t=REAL(i,kind=real_8)*gvec_com%b1(ir)+REAL(j,kind=real_8)*gvec_com%b2(ir)+REAL(k,kind=real_8)*gvec_com%b3(ir)
                   g2=g2+t*t
                ENDDO
#ifdef __SR8000 
                ! soption nopredicate 
#endif 
                IF (g2.GE.glow.AND.g2.LT.ghigh) THEN
                   ig=ig+1
                   IF (tnhgwish) tabg2(ig)=g2
                ENDIF
             ENDDO
          ENDDO
       ENDDO
       IF (iloop.EQ.1) THEN
          ncpw%ngw=ig
       ELSEIF (iloop.EQ.2) THEN
          ncpw%nhg=ig
       ENDIF
    ENDDO
    IF (tnhgwish) THEN
       CALL sort2(tabg2,ncpw%nhg,ntab)
       IF (ncpw%nhg.LT.nhgwish) THEN
          IF (paral%io_parent)&
               WRITE(6,*)&
               'THE CUTOFF ENERGY IS TOO SMALL TO HAVE NHG=',nhgwish
          IF (paral%io_parent)&
               WRITE(6,*) 'NHG=',ncpw%nhg
          CALL stopgm('NUMPW','INCREASE ECUT',& 
               __LINE__,__FILE__)
       ENDIF
       ! Change GCUT
       energy=tabg2(nhgwish)+delta
       ! Use old value if MESH option
       parm%nr1=nr1f
       parm%nr2=nr2f
       parm%nr3=nr3f
       ! Change GCUT, NR1,NR2, and NR3
       CALL change_gcut(energy)
       tnhgwish=.FALSE.
       DEALLOCATE(tabg2,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(ntab,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       GOTO 30
    ENDIF
    CALL tihalt('     NUMPW',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE numpw
  ! ==================================================================
  SUBROUTINE change_gcut(energy)
    ! ==--------------------------------------------------------------==
    ! == CHANGE ECUT                                                  ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: energy

    REAL(real_8), PARAMETER                  :: betaelmax = 1.e33_real_8 

    INTEGER                                  :: mnr1, mnr2, mnr3
    INTEGER, SAVE                            :: ibetapf = 0
    REAL(real_8)                             :: aa1, aa2, aa3

! ==--------------------------------------------------------------==

    gvec_com%gcutw=energy/dual00%cdual
    ! GCUT has to be bigger than GCUTWMAX
    gvec_com%gcut=MAX(energy,gcutwmax)
    gvec_com%gcutw=gvec_com%gcut/dual00%cdual
    cntr%ecut=gvec_com%gcutw*parm%tpiba2
    ! ==--------------------------------------------------------------==
    ! ==  REAL SPACE MESH (from setsc.F)                              ==
    ! ==--------------------------------------------------------------==
    IF (parm%ibrav.EQ.1) THEN
       ! Simple cubic
       mnr1=NINT(parm%alat/pi*SQRT(dual00%cdual*cntr%ecut)+0.5_real_8)
       mnr2=mnr1
       mnr3=mnr1
    ELSEIF (parm%ibrav.EQ.2) THEN
       ! Face centred cubic
       mnr1=NINT(0.75_real_8*parm%alat/pi*SQRT(dual00%cdual*cntr%ecut)+0.5_real_8)
       mnr2=mnr1
       mnr3=mnr1
    ELSEIF (parm%ibrav.EQ.3) THEN
       ! Body centred cubic
       ! 0.85 too short (sometimes..Thierry Deutsch)
       ! MNR1=NINT(0.85_real_8*ALAT/PI*SQRT(CDUAL*ECUT)+0.5_real_8)
       mnr1=NINT(0.875_real_8*parm%alat/pi*SQRT(dual00%cdual*cntr%ecut)+0.5_real_8)
       mnr2=mnr1
       mnr3=mnr1
    ELSE
       aa1=parm%alat
       aa2=parm%alat*cell_com%celldm(2)
       aa3=parm%alat*cell_com%celldm(3)
       mnr1=NINT(aa1/pi*SQRT(dual00%cdual*cntr%ecut)+0.5_real_8)
       mnr2=NINT(aa2/pi*SQRT(dual00%cdual*cntr%ecut)+0.5_real_8)
       mnr3=NINT(aa3/pi*SQRT(dual00%cdual*cntr%ecut)+0.5_real_8)
    ENDIF
    ! For MESH option.
    mnr1=MAX(mnr1,parm%nr1)
    mnr2=MAX(mnr2,parm%nr2)
    mnr3=MAX(mnr3,parm%nr3)
    parm%nr1=fftchk(mnr1,2)
    parm%nr2=fftchk(mnr2,2)
    parm%nr3=fftchk(mnr3,2)
    ! ==--------------------------------------------------------------==
    ! TROTTER FACTORISATION: BETAP
    fint1%p=0._real_8
    IF (fint1%ttrot) THEN
       IF (ibetapf.EQ.0) THEN
          IF (fint1%betap.EQ.0.0_real_8) THEN
             ! BETAP has to be set
             ibetapf=1
          ELSE
             ! BETAP sets by the user
             ibetapf=-1
          ENDIF
       ENDIF
       IF (ibetapf.EQ.1) fint1%betap = (parm%alat/REAL(parm%nr1,kind=real_8))**2/3._real_8
       IF (fint1%betael.LT.betaelmax) THEN
          fint1%p = fint1%betael/fint1%betap
       ELSE
          fint1%p = 0
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE change_gcut
  ! ==================================================================
  FUNCTION gkpwf(b1,b2,b3,tadd)
    ! ==--------------------------------------------------------------==
    ! == Calculate the maximum K point amplitude in Brillouin zone    ==
    ! == in order to apply a spherical cutoff                         ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: b1(3), b2(3), b3(3)
    LOGICAL                                  :: tadd
    REAL(real_8)                             :: gkpwf

    INTEGER                                  :: i, j
    REAL(real_8)                             :: b123l, b12l, b13l, b1l, b23l, &
                                                b2l, b3l

    IF (.NOT.tadd) THEN
       gkpwf=0._real_8
       RETURN
    ENDIF
    ! ==--------------------------------------------------------------==
    b1l=b1(1)**2+b1(2)**2+b1(3)**2
    b2l=b2(1)**2+b2(2)**2+b2(3)**2
    b3l=b3(1)**2+b3(2)**2+b3(3)**2
    gkpwf=MAX(b1l,b2l,b3l)
    ! Diagonal 1,2
    DO i=-1,1,2
       b12l=(b1(1) + i*b2(1))**2 +&
            (b1(2) + i*b2(2))**2 +&
            (b1(3) + i*b2(3))**2
       IF (b12l.GT.gkpwf) gkpwf=b12l
    ENDDO
    ! Diagonal 1,3
    DO i=-1,1,2
       b13l=(b1(1) + i*b3(1))**2 +&
            (b1(2) + i*b3(2))**2 +&
            (b1(3) + i*b3(3))**2
       IF (b13l.GT.gkpwf) gkpwf=b13l
    ENDDO
    ! Diagonal 2,3
    DO i=-1,1,2
       b23l=(b2(1) + i*b3(1))**2 +&
            (b2(2) + i*b3(2))**2 +&
            (b2(3) + i*b3(3))**2
       IF (b23l.GT.gkpwf) gkpwf=b23l
    ENDDO
    ! Diagonal 1,2,3
    DO j=-1,1,2
       DO i=-1,1,2
          b123l=(b1(1) + j*b2(1) + i*b3(1))**2&
               +(b1(2) + j*b2(2) + i*b3(2))**2&
               +(b1(3) + j*b2(3) + i*b3(3))**2
          IF (b123l.GT.gkpwf) gkpwf=b123l
       ENDDO
    ENDDO
    gkpwf=gkpwf/4._real_8
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION gkpwf
  ! ==================================================================

END MODULE numpw_utils
