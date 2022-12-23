MODULE meta_localizespin_utils
  USE cell,                            ONLY: cell_com
  USE cnst_dyn,                        ONLY: &
       cv_associated, cv_dyn, en_harm_locs, icv_spin, kharm, lmeta, ncolvar, &
       rcc, rcc0
  USE coor,                            ONLY: tau0
  USE fft_maxfft,                      ONLY: maxfft
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE latgen_utils,                    ONLY: latgen
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai
  USE prcp,                            ONLY: prcp_com
  USE spin,                            ONLY: clsd
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             parap,&
                                             parm,&
                                             spar

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: localizespin

CONTAINS

  ! ==================================================================
  SUBROUTINE localizespin(rhoe,v)
    ! ==--------------------------------------------------------------==


    ! ..Arguments
    REAL(real_8) :: rhoe(fpar%kr1,fpar%kr2s,fpar%kr3s,clsd%nlsd)
    COMPLEX(real_8)                          :: v(maxfft,clsd%nlsd)

    COMPLEX(real_8)                          :: ci, zeta(3)
    INTEGER                                  :: i, iat, icv, ijk, ir, isp, j, &
                                                k, nnat, nr1ve, nr2ve, nr3ve, &
                                                nsri, nsrj
    INTEGER, SAVE                            :: it = 0
    REAL(real_8) :: ah(3,3), ahi(3,3), ak, arw1, arw2, arw3, cc1(3), cc2(3), &
      cc3(3), cubecenter(3), den(3), dh, dist, dist0, dist1, dnat, drc(3), &
      dx, dy, dz, lowerleft(3), num(3), omegan, omenorm, oneby2pi, pi, pi2, &
      rhos, ss(3), tgbx, tgby, tgbz, theta(3), vvt, zetaim(3), zetare(3)

#ifdef __VECTOR
    ! ==--------------------------------------------------------------==
    nr1ve=parm%nr3
    nr2ve=nr1ve*parm%nr2
    nr3ve=nr2ve*parm%nr1
#endif

    ! ==--------------------------------------------------------------==

    it=it+1

    ! ..set up the constants
    pi=4._real_8*ATAN(1._real_8)
    pi2=2.0_real_8*pi
    oneby2pi=1._real_8/pi2 ! 1/(2 pi)
    ci=(0._real_8,1._real_8)

    ! if(it.eq.1)RETURN
    IF (.NOT. lmeta%tlocalizespin .OR. .NOT. cntl%tlsd) THEN
       RETURN
    ENDIF

    ! ..construction of a(i) lattice vector and h(i,j) matrix
    IF (cntl%tprcp.OR.cntl%tpres) THEN ! variable cell case (OMEGA <> OMEGAN)
       ! .. are you sure about this?
       CALL latgen(parm%ibrav,prcp_com%cellrf,cc1,cc2,cc3,omegan)
    ELSE                    ! fixed cell case
       CALL latgen(parm%ibrav,cell_com%celldm,cc1,cc2,cc3,parm%omega)
    ENDIF
    ah(1,1)=cc1(1)
    ah(2,1)=cc1(2)
    ah(3,1)=cc1(3)
    ah(1,2)=cc2(1)
    ah(2,2)=cc2(2)
    ah(3,2)=cc2(3)
    ah(1,3)=cc3(1)
    ah(2,3)=cc3(2)
    ah(3,3)=cc3(3)

    CALL matinv3(ah,ahi,dh)
    ! Compute z(k)
    omenorm=parm%omega/REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)
    nsri=parm%nr1
    nsrj=parm%nr2*nsri

    DO ir=1,3 ! zeroing
       zeta(ir)=(0._real_8,0._real_8)! cmb - bugfix
       cubecenter(ir)=0.0_real_8
    ENDDO

    ! Geometric Center of atomic coordinates
    nnat=0
    DO isp=1,ions1%nsp
       DO iat=1,ions0%na(isp)
          cubecenter(1) = cubecenter(1) + tau0(1,iat,isp)
          cubecenter(2) = cubecenter(2) + tau0(2,iat,isp)
          cubecenter(3) = cubecenter(3) + tau0(3,iat,isp)
          nnat = nnat + 1
       ENDDO
    ENDDO
    dnat=1._real_8/REAL(nnat,kind=real_8)
    cubecenter(1) =  dnat*cubecenter(1)
    cubecenter(2) =  dnat*cubecenter(2)
    cubecenter(3) =  dnat*cubecenter(3)

    ! write(6,'(A,3f12.6)') 'center', cubecenter
    ! Lower Left corner of the grid centered on the molecule
    DO ir=1,3
       lowerleft(ir) = cubecenter(ir)&
            - 0.5_real_8*(cc1(ir)+cc2(ir)+cc3(ir))
    ENDDO
    ! write(6,'(A,3f12.6)') 'lowerleft',lowerleft

#ifdef __VECTOR
    DO ijk = 1,nr3ve
       k=1+MOD(ijk-1,nr1ve)
       j=1+INT(MOD(ijk-1,nr2ve)/nr1ve)
       i=1+INT((ijk-1)/nr2ve)

       ! if (i.eq.NR1+1.OR.j.eq.NR2+1.OR.k.eq.NR3+1) CYCLE
       ! rhos=abs(RHOE(IR,1) - 2.0_real_8*RHOE(IR,2))    ! spindensity justincase
       rhos=ABS(rhoe(i,j,k,1) - 2.0_real_8*rhoe(i,j,k,2))

       tgbx=pi2*REAL(parap%nrxpl(parai%mepos,1)+i-1,kind=real_8)/REAL(spar%nr1s,kind=real_8)
       tgby=pi2*REAL(j-1,kind=real_8)/REAL(spar%nr2s,kind=real_8)
       tgbz=pi2*REAL(k-1,kind=real_8)/REAL(spar%nr3s,kind=real_8)

       zeta(1)=zeta(1)+rhos*(COS(tgbx)+ci*SIN(tgbx))
       zeta(2)=zeta(2)+rhos*(COS(tgby)+ci*SIN(tgby))
       zeta(3)=zeta(3)+rhos*(COS(tgbz)+ci*SIN(tgbz))

    ENDDO
#else
    DO i = 1,parm%nr1
       DO j = 1,parm%nr2
          DO k = 1,parm%nr3
             ! if (i.eq.NR1+1.OR.j.eq.NR2+1.OR.k.eq.NR3+1) CYCLE
             ! rhos=abs(RHOE(IR,1) - 2.0_real_8*RHOE(IR,2))    ! spindensity justincase
             rhos=ABS(rhoe(i,j,k,1) - 2.0_real_8*rhoe(i,j,k,2))
             tgbx=pi2*REAL(parap%nrxpl(parai%mepos,1)+i-1,kind=real_8)/REAL(spar%nr1s,kind=real_8)
             tgby=pi2*REAL(j-1,kind=real_8)/REAL(spar%nr2s,kind=real_8)
             tgbz=pi2*REAL(k-1,kind=real_8)/REAL(spar%nr3s,kind=real_8)

             zeta(1)=zeta(1)+rhos*(COS(tgbx)+ci*SIN(tgbx))
             zeta(2)=zeta(2)+rhos*(COS(tgby)+ci*SIN(tgby))
             zeta(3)=zeta(3)+rhos*(COS(tgbz)+ci*SIN(tgbz))
          ENDDO
       ENDDO
    ENDDO
#endif

    zeta(1)=zeta(1)*omenorm
    zeta(2)=zeta(2)*omenorm
    zeta(3)=zeta(3)*omenorm

    ! write(6,'(A,I3,6f10.4)') 'zeta1 ',MEPOS,ZETA(1),ZETA(2),ZETA(3)

    CALL mp_sum(zeta,3,parai%allgrp)

    ! write(6,'(A,I3,6f10.4)') 'zeta2 ',MEPOS,ZETA(1),ZETA(2),ZETA(3)

    ! ..now compute Im{ln z(i)}= atan[aimag{z(i)}/real{z(i)}]
    zetare(1) = REAL(zeta(1))
    zetaim(1) = AIMAG(zeta(1))
    theta(1)  = ATAN(zetaim(1)/zetare(1))
    IF (zetare(1).LT.0._real_8) theta(1) = theta(1) + pi
    theta(1) = theta(1) * oneby2pi
    zetare(2) = REAL(zeta(2))
    zetaim(2) = AIMAG(zeta(2))
    theta(2)  = ATAN(zetaim(2)/zetare(2))
    IF (zetare(2).LT.0._real_8) theta(2) = theta(2) + pi
    theta(2)  = theta(2) * oneby2pi
    zetare(3) = REAL(zeta(3))
    zetaim(3) = AIMAG(zeta(3))
    theta(3)  = ATAN(zetaim(3)/zetare(3))
    IF (zetare(3).LT.0._real_8) theta(3) = theta(3) + pi
    theta(3)  = theta(3) * oneby2pi

    ! ..and eventually r_i[spinden]
    rcc(1)= & ! lowerleft(1)+
         (ah(1,1)*theta(1)+ah(1,2)*theta(2)+ah(1,3)*theta(3))
    rcc(2)= & ! lowerleft(2)+
         (ah(2,1)*theta(1)+ah(2,2)*theta(2)+ah(2,3)*theta(3))
    rcc(3)= & ! lowerleft(3)+
         (ah(3,1)*theta(1)+ah(3,2)*theta(2)+ah(3,3)*theta(3))


    en_harm_locs = 0.0_real_8
    IF (.NOT. cv_associated) RETURN

    num(1)=zetaim(1)
    num(2)=zetaim(2)
    num(3)=zetaim(3)
    den(1)=zetare(1)
    den(2)=zetare(2)
    den(3)=zetare(3)
    ss(1)=num(1)/den(1)
    ss(2)=num(2)/den(2)
    ss(3)=num(3)/den(3)


    DO icv = 1,ncolvar
       IF (icv_spin(icv) .NE. 0) THEN

          dist0=cv_dyn(icv)
          ak=kharm(icv)

          ! .. of the potential with respect to rcc

          dx = (rcc(1)-rcc0(1,icv))
          dy = (rcc(2)-rcc0(2,icv))
          dz = (rcc(3)-rcc0(3,icv))

          ! mb-ale Most likely the line below is a bug. We do not want PBC...
          ! mb-ale   CALL PBC(DX,DY,DZ,DX,DY,DZ,1,APBC,IBRAV)
          dist=SQRT(dx*dx+dy*dy+dz*dz)

          ! Energy contribution due to the metadynamic harmonic terms 
          ! (only those related to the CV for the localization of the spin density)
          ! This contribution should appear in the electronic energy 
          en_harm_locs = en_harm_locs+&
               ak*(dist -dist0)*(dist-dist0)

          ! write(6,'(I3,3f12.6)') 
          ! &          ICV,  dist0,  dist, ak*(dist -dist0)*(dist-dist0)

          dist1=(dist-dist0)/dist

          drc(1)=2.0_real_8*ak*dx*dist1
          drc(2)=2.0_real_8*ak*dy*dist1
          drc(3)=2.0_real_8*ak*dz*dist1

          ! ...electronic contribution d(vs)/d[rho(x)]
          DO ir=1,fpar%nnr1
             i=1+MOD(ir-1,nsri)
             j=1+INT(MOD(ir-1,nsrj)/nsri)
             k=1+INT((ir-1)/nsrj)
             IF (i.EQ.parm%nr1+1.OR.j.EQ.parm%nr2+1.OR.k.EQ.parm%nr3+1) CYCLE
             tgbx=pi2*REAL(parap%nrxpl(parai%mepos,1)+i-1,kind=real_8)/REAL(spar%nr1s,kind=real_8)
             tgby=pi2*REAL(j,kind=real_8)/REAL(spar%nr2s,kind=real_8)
             tgbz=pi2*REAL(k,kind=real_8)/REAL(spar%nr3s,kind=real_8)

             arw1=1.0_real_8/(1._real_8+ss(1)**2)/den(1)**2&
                  *(den(1)*SIN(tgbx)-num(1)*COS(tgbx))*omenorm
             arw2=1.0_real_8/(1._real_8+ss(2)**2)/den(2)**2&
                  *(den(2)*SIN(tgby)-num(2)*COS(tgby))*omenorm
             arw3=1.0_real_8/(1._real_8+ss(3)**2)/den(3)**2&
                  *(den(3)*SIN(tgbz)-num(3)*COS(tgbz))*omenorm

             vvt=oneby2pi*(&
                  drc(1)*(ah(1,1)*arw1+ah(1,2)*arw2+ah(1,3)*arw3)&
                  +drc(2)*(ah(2,1)*arw1+ah(2,2)*arw2+ah(2,3)*arw3)&
                  +drc(3)*(ah(3,1)*arw1+ah(3,2)*arw2+ah(3,3)*arw3))
             v(ir,1)=v(ir,1)+vvt! cmb - spin/charge sign checked
             v(ir,2)=v(ir,2)-vvt! cmb - spin/charge sign checked
          END DO
       ENDIF! ICV_SPIN
    ENDDO ! NCOLVAR
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE localizespin
  ! ==================================================================
  SUBROUTINE matinv3(am,ai,dh)
    ! ==--------------------------------------------------------------==
    ! ==  3x3 real matrix inversion: AM = incoming matrix             ==
    ! ==                             AI = inverted matrix             ==
    ! ==                             DH = determinant                 ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: am(3,3), ai(3,3), dh

    REAL(real_8)                             :: d11, d12, d13, d21, d22, d23, &
                                                d31, d32, d33

! ==--------------------------------------------------------------==

    d11=am(2,2)*am(3,3)-am(2,3)*am(3,2)
    d12=am(2,3)*am(3,1)-am(2,1)*am(3,3)
    d13=am(2,1)*am(3,2)-am(3,1)*am(2,2)
    d21=am(3,2)*am(1,3)-am(1,2)*am(3,3)
    d22=am(1,1)*am(3,3)-am(1,3)*am(3,1)
    d23=am(3,1)*am(1,2)-am(1,1)*am(3,2)
    d31=am(1,2)*am(2,3)-am(2,2)*am(1,3)
    d32=am(1,3)*am(2,1)-am(1,1)*am(2,3)
    d33=am(1,1)*am(2,2)-am(1,2)*am(2,1)
    ! ==--------------------------------------------------------------==
    dh=am(1,1)*d11+am(1,2)*d12+am(1,3)*d13
    ! mb      IF(DH.eq.0.) THEN ! this part is now removed
    ! mb        WRITE(6,*)' Singular matrix with Det=0 in module MATINV'
    ! mb        STOP
    ! mb      ENDIF
    ! ==--------------------------------------------------------------==
    ai(1,1)=d11/dh
    ai(2,2)=d22/dh
    ai(3,3)=d33/dh
    ai(1,2)=d21/dh
    ai(1,3)=d31/dh
    ai(2,1)=d12/dh
    ai(2,3)=d32/dh
    ai(3,1)=d13/dh
    ai(3,2)=d23/dh
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE matinv3
  ! ==================================================================

END MODULE meta_localizespin_utils
