MODULE hipin_utils
  USE cnst,                            ONLY: fpi,&
                                             pi
  USE cppt,                            ONLY: hg,&
                                             hgpot,&
                                             hipz,&
                                             nr1h,&
                                             nr2h,&
                                             nr3h,&
                                             nr3pl
  USE error_handling,                  ONLY: stopgm
  USE fftmain_utils,                   ONLY: mltfft
  USE geq0mod,                         ONLY: geq0
  USE gvec,                            ONLY: gvec_com
  USE isos,                            ONLY: isos1,&
                                             isos2,&
                                             isos3
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: parai,&
                                             paral
  USE prmem_utils,                     ONLY: prmem
  USE special_functions,               ONLY: cp_erf
  USE system,                          ONLY: ncpw,&
                                             parap,&
                                             parm,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: hipin


CONTAINS

  ! ==================================================================
  SUBROUTINE hipin
    ! ==--------------------------------------------------------------==
    ! == HOCKNEY INTERACTION POTENTIAL: INITIALIZATION                ==
    ! ==--------------------------------------------------------------==

    CHARACTER(*), PARAMETER                  :: procedureN = 'hipin'

    COMPLEX(real_8), ALLOCATABLE             :: cscr1(:), cscr2(:), cscr3(:), &
                                                cscr4(:)
    INTEGER                                  :: ierr, ihip1, ihip2, ihip3, &
                                                ihip4, ihip5, ihip6, ihip7
    REAL(real_8), ALLOCATABLE                :: scr1(:), scr2(:), scr3(:)

! ==--------------------------------------------------------------==

    nr3pl=parap%nrzpl(parai%mepos,2)-parap%nrzpl(parai%mepos,1)+1
    ! IF(PARENT) NR3MX=NR3M+1
    ALLOCATE(hgpot(spar%nr1s+1,spar%nr2s+1,nr3pl),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(hipz(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! Definition of NR1H, NR2H, NR3H (cppt.inc)
    nr1h = 2*spar%nr1s
    nr2h = 2*spar%nr2s
    nr3h = 2*spar%nr3s
    ! !!   IF(TONED.OR.TTWOD) NR1H=NR1S
    ! !!   IF(TTWOD)          NR2H=NR2S
    IF (isos1%toned) nr1h=spar%nr1s
    IF (isos1%ttwod) THEN
       IF (isos3%snormal.EQ.1) THEN
          nr2h=spar%nr2s
          nr3h=spar%nr3s
       ELSE IF (isos3%snormal.EQ.2) THEN
          nr1h=spar%nr1s
          nr3h=spar%nr3s
       ELSE IF (isos3%snormal.EQ.3) THEN
          nr1h=spar%nr1s
          nr2h=spar%nr2s
       ENDIF
    ENDIF
    ! ==------------------------------------------------------------==
    ! == SCRATCH SPACE                                              ==
    ! ==------------------------------------------------------------==
    ihip1 = nr1h+(MOD(nr1h,2))
    ihip2 = nr2h+(MOD(nr2h,2))
    ihip3 = nr3h+(MOD(nr3h,2))

    ihip4 = nr2h*nr3h
    ihip5 = nr1h*(spar%nr2s+1)*nr3pl
    ihip6 = nr2h*nr3h
    ihip7 = nr1h*(spar%nr2s+1)*nr3pl
    ALLOCATE(scr1(ihip1),scr2(ihip2),scr3(ihip3),cscr1(ihip4),cscr2(ihip5),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ==------------------------------------------------------------==
    ALLOCATE(cscr3(ihip6),cscr4(ihip7),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(cscr3)
    CALL zeroing(cscr4)

    CALL hockney(scr1,scr2,scr3,cscr1,cscr2,cscr3,cscr4)

    DEALLOCATE(cscr3,cscr4,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(scr1,scr2,scr3,cscr1,cscr2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    CALL zhip
    IF (paral%parent) CALL prmem('   CLUSTER')
    ! ==--------------------------------------------------------------==
  END SUBROUTINE hipin
  ! ==================================================================
  SUBROUTINE hockney(r1,r2,r3,potyz,potx,potyza,potxa)
    ! ==--------------------------------------------------------------==
    ! ==  COMPUTES THE HOCKNEY INTERACTION POTENTIAL                  ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: r1(:), r2(:), r3(:)
    COMPLEX(real_8)                          :: potyz(:), potx(:), potyza(:), &
                                                potxa(:)

    INTEGER                                  :: i, isub, ix, ixyz, iy, iyz, &
                                                iz, iz0, k, kx, ky, kz, nnn, &
                                                nr1hh, nr2hh, nr3hh, nr3in
    LOGICAL                                  :: trealsp
    REAL(real_8)                             :: dr1, dr2, dr3, dvol, dx, dy, &
                                                dz, g2, gxy, gxy2, gz, gz2, &
                                                potg, rpox, scale, spot, u0, z

! Variables
! ==--------------------------------------------------------------==

    CALL tiset('   HOCKNEY',isub)
    nr1hh=nr1h/2
    nr2hh=nr2h/2
    nr3hh=nr3h/2
    dvol = REAL((nr1h/spar%nr1s) * (nr2h/spar%nr2s) * (nr3h/spar%nr3s),kind=real_8) * parm%omega
    trealsp = (.NOT.(isos1%ttwod.OR.isos1%toned))
    dr1 = 1.0_real_8/REAL(spar%nr1s,kind=real_8)
    dr2 = 1.0_real_8/REAL(spar%nr2s,kind=real_8)
    dr3 = 1.0_real_8/REAL(spar%nr3s,kind=real_8)
    u0 = 2.0_real_8 / (SQRT(pi)*isos2%rcerf)
    DO i=1,nr1h
       k=i-1
       IF (k.GT.nr1hh) k=nr1h-k
       r1(i)=dr1*REAL(k,kind=real_8)
    ENDDO
    DO i=1,nr2h
       k=i-1
       IF (k.GT.nr2hh) k=nr2h-k
       r2(i)=dr2*REAL(k,kind=real_8)
    ENDDO
    DO i=1,nr3h
       k=i-1
       IF (k.GT.nr3hh) k=nr3h-k
       r3(i)=dr3*REAL(k,kind=real_8)
    ENDDO
    ! ==--------------------------------------------------------------==
    nr3in=parap%nrzpl(parai%mepos,1)
    IF (trealsp) THEN
       DO ix=1,nr1h
          !$omp parallel do private(IZ,IY,IYZ,DX,DY,DZ,RPOX,SPOT)
          DO iz=1,nr3h
             DO iy=1,nr2h
                iyz=iy+(iz-1)*nr2h
                dx=parm%a1(1)*r1(ix)+parm%a2(1)*r2(iy)+parm%a3(1)*r3(iz)
                dy=parm%a1(2)*r1(ix)+parm%a2(2)*r2(iy)+parm%a3(2)*r3(iz)
                dz=parm%a1(3)*r1(ix)+parm%a2(3)*r2(iy)+parm%a3(3)*r3(iz)
                rpox=SQRT(dx*dx+dy*dy+dz*dz)! focus on small G vectors
                spot=cp_erf(rpox/isos2%rcerf)
                IF (rpox.GT.1.e-20_real_8) potyz(iyz)=CMPLX(spot/rpox,0.0_real_8,kind=real_8)
             ENDDO
          ENDDO
          !$omp end parallel do
          IF (ix.EQ.1) potyz(1)=u0
          CALL dcopy(2*nr2h*nr3h,potyz(1),1,potyza(1),1)
          scale=1.0_real_8
          CALL mltfft('T','T',potyza(1:nr2h*nr3h),nr2h,nr3h,&
               & potyz(1:nr2h*nr3h),nr2h,nr3h,nr3h,nr2h,1,scale)
          iyz=1+(nr3in-1)*nr2h
          CALL dcopy(2*nr2h*nr3pl,potyz(iyz),1,potyza(iyz),1)
          CALL mltfft('N','N',potyza(iyz:iyz+nr2h*nr3pl-1),nr2h,nr3pl,&
               & potyz(iyz:iyz+nr2h*nr3pl-1),nr2h,nr3pl,nr2h,nr3pl,1,scale)
          !$omp parallel do private(IZ,IZ0,IY,IXYZ,IYZ) default(shared)
          DO iz=1,nr3pl
             iz0=nr3in-1+iz
             DO iy=1,spar%nr2s+1
                ixyz=ix+(iy-1)*nr1h+(iz-1)*nr1h*(spar%nr2s+1)
                iyz=iy+(iz0-1)*nr2h
                potx(ixyz)=potyz(iyz)
             ENDDO
          ENDDO
          !$omp end parallel do
       ENDDO
       scale=1.0_real_8/REAL(nr1h*nr2h*nr3h,kind=real_8)
       nnn=(spar%nr2s+1)*nr3pl
       CALL dcopy(2*nnn*nr1h,potx(1),1,potxa(1),1)
       CALL mltfft('N','N',potxa,nr1h,nnn,potx,nr1h,nnn,&
            NR1H,NNN,1,SCALE)
    ENDIF
    DO iz=1,nr3pl
       kz = iz-1 + nr3in-1
       IF (kz.GT.nr3hh) kz=kz-nr3h
       DO iy=1,spar%nr2s+1
          ky = iy-1
          IF (ky.GT.nr2hh) ky=ky-nr2h
          DO ix=1,spar%nr1s+1
             kx = ix-1
             IF (kx.GT.nr1hh) kx=kx-nr1h
             ixyz=ix+(iy-1)*nr1h+(iz-1)*nr1h*(spar%nr2s+1)
             IF (trealsp) THEN
                hgpot(ix,iy,iz)=REAL(potx(ixyz))
             ELSE
                gxy2 = ((gvec_com%b1(1)*kx+gvec_com%b2(1)*ky)**2+&
                     (gvec_com%b2(2)*KY+gvec_com%b1(2)*KX)**2)*parm%tpiba2
                gz = parm%tpiba*gvec_com%b3(3)*kz*spar%nr3s/nr3h
                gz2 = gz**2
                g2 = gxy2 + gz2
                gxy = SQRT(gxy2)
                z = parm%a3(3)*nr3h/spar%nr3s
                IF (ix.NE.1.OR.iy.NE.1.OR.iz.NE.1) THEN
                   potg = fpi/(g2*dvol)*EXP(-0.25_real_8*g2*isos2%rcerf*isos2%rcerf)
                   IF (isos1%ttwod) potg = potg *&
                        (1.0_real_8-COS(0.5_real_8*Z*GZ)*EXP(-0.5_real_8*Z*GXY))
                ELSE
                   g2 = 0.0_real_8
                   potg = 0.0_real_8
                ENDIF
                hgpot(ix,iy,iz)=potg
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    CALL tihalt('   HOCKNEY',isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE hockney
  ! ==================================================================
  SUBROUTINE zhip
    ! ==--------------------------------------------------------------==
    ! Variables
    INTEGER                                  :: ig, ig1
    REAL(real_8)                             :: g2, rc2

! ==--------------------------------------------------------------==

    rc2=isos2%rcerf*isos2%rcerf
    ! ==--------------------------------------------------------------==
    ig1=1
    IF (geq0) ig1=2
    !$omp parallel do private(IG,G2) shared(RC2)
    DO ig=ig1,ncpw%nhg
       g2=parm%tpiba2*hg(ig)
       hipz(ig)=fpi/g2*(1._real_8-EXP(-0.25_real_8*g2*rc2))! focus on large G v.
    ENDDO
    !$omp end parallel do
    IF (geq0) hipz(1)=pi*rc2
    ! ==--------------------------------------------------------------==
  END SUBROUTINE zhip
  ! ==================================================================

END MODULE hipin_utils
