MODULE saop_utils
  USE cnst,                            ONLY: pi
  USE error_handling,                  ONLY: stopgm
  USE functionals_utils,               ONLY: pz,&
                                             slaterx
  USE kinds,                           ONLY: real_8
  USE lsd_func_utils,                  ONLY: lsd_pz,&
                                             lsd_sx
  USE rho1ofr_utils,                   ONLY: rhosofr
  USE spin,                            ONLY: clsd,&
                                             spin_mod
  USE system,                          ONLY: cntl,&
                                             cntr,&
                                             fpar,&
                                             group,&
                                             ncpw
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: saop
  PUBLIC :: pbe_xed
  PUBLIC :: pbe_ced
  PUBLIC :: pbe_ceds
  PUBLIC :: lb94
  PUBLIC :: lb94s
  PUBLIC :: lb94m
  PUBLIC :: gllb

CONTAINS

  ! ==================================================================
  SUBROUTINE saop(c0,eigv,eref,focc,nstate,rhoe,grho,exc,&
       vsaop,vgllb,vlb94,rho1,psi)
    ! ==--------------------------------------------------------------==
    ! == CALCULATE SAOP POTENTIAL                                     ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(ncpw%ngw,*)
    REAL(real_8)                             :: eigv(*), eref(*), focc(*)
    INTEGER                                  :: nstate
    REAL(real_8) :: rhoe(fpar%nnr1,clsd%nlsd), grho(fpar%nnr1,*), exc, &
      vsaop(fpar%nnr1,*), vgllb(fpar%nnr1,clsd%nlsd), &
      vlb94(fpar%nnr1,clsd%nlsd), rho1(fpar%nnr1)
    COMPLEX(real_8)                          :: psi(:)

    REAL(real_8), PARAMETER                  :: kresp = 0.42_real_8 

    INTEGER                                  :: ik, ir, is
    REAL(real_8)                             :: ec, eex, en, ena, enb, ex, &
                                                weni, weni1, weni2

! ==--------------------------------------------------------------==

    IF (cntl%tlsd) THEN
       en=0._real_8
       DO is=1,spin_mod%nsup
          IF (focc(is).GT.1.e-4_real_8) ena=eigv(is)
       ENDDO
       DO is=spin_mod%nsup+1,nstate
          IF (focc(is).GT.1.e-4_real_8) enb=eigv(is)
       ENDDO
    ELSE
       ena=0._real_8
       enb=0._real_8
       DO is=1,nstate
          IF (focc(is).GT.1.e-4_real_8) en=eigv(is)
       ENDDO
    ENDIF
    ! 
    ! Vresponse
    IF (group%nogrp.GT.1) CALL stopgm("SAOP","NO GROUPS ALLOWED",& 
         __LINE__,__FILE__)
    CALL zeroing(vgllb)!,nnr1*clsd%nlsd)
    DO is=1,nstate
       CALL rhosofr(c0(:,is),rho1,psi)
       IF (focc(is).GT.1.e-4_real_8) THEN
          IF (cntl%tlsd) THEN
             IF (is.LE.spin_mod%nsup) THEN
                ena=eref(is)
                weni=SQRT(ena-eigv(is))*kresp*focc(is)
                ik=1
             ELSE
                enb=eref(is)
                weni=SQRT(enb-eigv(is))*kresp*focc(is)
                ik=2
             ENDIF
          ELSE
             en=eref(is)
             weni=SQRT(en-eigv(is))*kresp*focc(is)
             ik=1
          ENDIF
          DO ir=1,fpar%nnr1
             IF (rhoe(ir,ik).GT.1.e-12_real_8) THEN
                vgllb(ir,ik) = vgllb(ir,ik) + weni * rho1(ir)/rhoe(ir,ik)
             ENDIF
          ENDDO
       ENDIF
    ENDDO
    ! PBE exchange energy density
    IF (cntl%tlsd) THEN
       DO ir=1,fpar%nnr1
          CALL pbe_xed(2._real_8*rhoe(ir,1),4._real_8*grho(ir,1),ex)
          vgllb(ir,1) = vgllb(ir,1) + 2._real_8*ex
       ENDDO
       DO ir=1,fpar%nnr1
          CALL pbe_xed(2._real_8*rhoe(ir,2),4._real_8*grho(ir,3),ex)
          vgllb(ir,2) = vgllb(ir,2) + 2._real_8*ex
       ENDDO
    ELSE
       DO ir=1,fpar%nnr1
          CALL pbe_xed(rhoe(ir,1),grho(ir,1),ex)
          vgllb(ir,1) = vgllb(ir,1) + 2._real_8*ex
       ENDDO
    ENDIF
    ! PBE correlation energy density
    IF (cntl%tlsd) THEN
       DO ir=1,fpar%nnr1
          CALL pbe_ceds(rhoe(ir,1),rhoe(ir,2),grho(ir,1),&
               grho(ir,2),grho(ir,3),ec)
          vgllb(ir,1) = vgllb(ir,1) + 2._real_8*ec
          vgllb(ir,2) = vgllb(ir,2) + 2._real_8*ec
       ENDDO
    ELSE
       DO ir=1,fpar%nnr1
          CALL pbe_ced(rhoe(ir,1),grho(ir,1),ec)
          vgllb(ir,1) = vgllb(ir,1) + 2._real_8*ec
       ENDDO
    ENDIF
    ! 
    ! VLB94
    CALL zeroing(vlb94)!,nnr1*clsd%nlsd)
    IF (cntl%tlsd) THEN
       DO ir=1,fpar%nnr1
          CALL lb94s(rhoe(ir,1),rhoe(ir,2),grho(ir,1),&
               grho(ir,3),vlb94(ir,1),vlb94(ir,2))
       ENDDO
    ELSE
       DO ir=1,fpar%nnr1
          CALL lb94(rhoe(ir,1),grho(ir,1),vlb94(ir,1))
       ENDDO
    ENDIF
    ! 
    ! Statistical average
    IF (group%nogrp.GT.1) CALL stopgm("SAOP","NO GROUPS ALLOWED",& 
         __LINE__,__FILE__)
    DO is=1,nstate
       CALL rhosofr(c0(:,is),rho1,psi)
       IF (focc(is).GT.1.e-4_real_8) THEN
          IF (cntl%tlsd) THEN
             IF (is.LE.spin_mod%nsup) THEN
                eex=EXP(-2._real_8*(ena-eigv(is))**2)
                ik=1
             ELSE
                eex=EXP(-2._real_8*(enb-eigv(is))**2)
                ik=2
             ENDIF
          ELSE
             eex=EXP(-2._real_8*(en-eigv(is))**2)
             ik=1
          ENDIF
          weni1=focc(is)*eex
          weni2=focc(is)*(1._real_8-eex)
          DO ir=1,fpar%nnr1
             IF (rhoe(ir,ik).GT.1.e-12_real_8) THEN
                vsaop(ir,ik)=vsaop(ir,ik)+ rho1(ir)/rhoe(ir,ik) *&
                     (weni1*vlb94(ir,ik)+weni2*vgllb(ir,ik))
             ENDIF
          ENDDO
       ENDIF
    ENDDO
    ! 
    ! calculate EXC (not implemented)
    IF (cntl%tlsd) THEN
       exc=0._real_8
    ELSE
       exc=0._real_8
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE saop
  ! ==================================================================
  SUBROUTINE pbe_xed(rhoe,grho,ex)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rhoe, grho, ex

    REAL(real_8), PARAMETER :: kappa = 0.804_real_8 , mu = 0.21951_real_8 , &
      ob3 = 1._real_8/3._real_8 , small = 1.e-14_real_8 

    REAL(real_8)                             :: exunif, fx, kf, s2, x

! ==--------------------------------------------------------------==

    IF ( rhoe .GT. small ) THEN
       kf = (3._real_8*pi*pi*rhoe)**ob3
       exunif = -3._real_8/(4._real_8*pi) * kf
    ELSE
       exunif = 0._real_8
    ENDIF
    IF ( rhoe .GT. cntr%gceps .AND. SQRT(grho) .GT. small ) THEN
       s2 = grho/(2._real_8*kf*rhoe)**2
       x = mu*s2
       fx = 1._real_8 + kappa - kappa/(1._real_8 + x/kappa)
    ELSE
       fx = 1._real_8
    ENDIF
    ex = exunif*fx
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE pbe_xed
  ! ==================================================================
  SUBROUTINE pbe_ced(rhoe,grho,ec)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rhoe, grho, ec

    REAL(real_8), PARAMETER :: be = 0.06672455060314922_real_8, &
      ga = 0.031090690869654895_real_8 , ob3 = 1._real_8/3._real_8 , &
      small = 1.e-14_real_8 

    INTEGER                                  :: iflg
    REAL(real_8)                             :: a, aa, af, eunif, expe, h0, &
                                                rs, s1, t, vc, xkf, xks, xy, y

! ==--------------------------------------------------------------==

    IF ( rhoe .GT. small ) THEN
       rs    = (3._real_8/(4._real_8*pi*rhoe))**ob3
       iflg  = 2
       IF (rs.LT.1.0_real_8) iflg=1
       CALL pz(rs,eunif,vc,iflg)
    ELSE
       eunif = 0._real_8
    ENDIF
    aa    = grho
    a     = SQRT(aa)
    IF ( rhoe .GT. cntr%gceps .AND. a .GT. small ) THEN
       xkf   = (9._real_8*pi/4._real_8)**ob3/rs
       xks   = SQRT(4._real_8*xkf/pi)
       t     = a/(2._real_8*xks*rhoe)
       expe  = EXP(-eunif/ga)
       af    = be/ga * (1._real_8/(expe-1._real_8))
       y     = af*t*t
       xy    = (1._real_8+y)/(1._real_8+y+y*y)
       s1    = 1._real_8+be/ga*t*t*xy
       h0    = ga * LOG(s1)
       ec    = eunif + h0
    ELSE
       ec    = eunif
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE pbe_ced
  ! ==================================================================
  SUBROUTINE pbe_ceds(rhoa,rhob,grhoaa,grhoab,grhobb,ec)
    REAL(real_8)                             :: rhoa, rhob, grhoaa, grhoab, &
                                                grhobb, ec

    REAL(real_8), PARAMETER :: be = 0.06672455060314922_real_8, &
      ga = 0.031090690869654895_real_8 , ob3 = 1._real_8/3._real_8 , &
      small = 1.e-14_real_8 

    INTEGER                                  :: iflg
    REAL(real_8)                             :: a, aa, af, eta, eunif, expe, &
                                                grho, h0, phi, phi3, rho, rs, &
                                                s1, t, vca, vcb, xkf, xks, &
                                                xy, y

! ==--------------------------------------------------------------==

    rho=rhoa+rhob
    eta=(rhoa-rhob)/rho
    grho=grhoaa+2._real_8*grhoab+grhobb
    IF (ABS(eta).GT.1._real_8) eta=SIGN(1.0_real_8,eta)
    IF ( rho .GT. small ) THEN
       rs = (0.75_real_8/(pi*rho))**ob3
       iflg=2
       IF (rs.LT.1.0_real_8) iflg=1
       CALL lsd_pz(rs,eta,eunif,vca,vcb,iflg)
    ELSE
       eunif = 0._real_8
    ENDIF
    aa    = grho
    a     = SQRT(aa)
    IF ( rho .GT. cntr%gceps .AND. a .GT. small ) THEN
       phi   = 0.5_real_8*((1._real_8+eta)**(2*ob3)+(1._real_8-eta)**(2*ob3))
       phi3  = phi*phi*phi
       xkf   = (9._real_8*pi/4._real_8)**ob3/rs
       xks   = SQRT(4._real_8*xkf/pi)
       t     = a/(2._real_8*xks*rho*phi)
       expe  = EXP(-eunif/(phi3*ga))
       af    = be/ga * (1._real_8/(expe-1._real_8))
       y     = af*t*t
       xy    = (1._real_8+y)/(1._real_8+y+y*y)
       s1    = 1._real_8+be/ga*t*t*xy
       h0    = ga*phi3 * LOG(s1)
    ELSE
       h0    = 0._real_8
    ENDIF
    ec    = eunif + h0
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE pbe_ceds
  ! ==================================================================
  SUBROUTINE lb94(rho,grho,vlb94)
    REAL(real_8)                             :: rho, grho, vlb94

    REAL(real_8), PARAMETER :: alpha = 1.19_real_8, beta = 0.01_real_8 , &
      ob3 = 1._real_8/3._real_8 , small = 1.e-14_real_8 

    INTEGER                                  :: iflg
    REAL(real_8)                             :: bln, ec, ex, rs, sgro, vas, &
                                                vc, vx, xs

! ==--------------------------------------------------------------==

    IF (rho.GT.small) THEN
       rs = (0.75_real_8/(pi*rho))**ob3
       CALL slaterx(rho,ex,vx,2._real_8*ob3)
       iflg=2
       IF (rs.LT.1.0_real_8) iflg=1
       CALL pz(rs,ec,vc,iflg)
       IF (rho.GT.cntr%gceps) THEN
          sgro=0.5_real_8*SQRT(grho)
          xs=sgro/(0.5_real_8*rho)**(4._real_8*ob3)
          bln=1._real_8+3._real_8*beta*xs*LOG(xs+SQRT(xs*xs+1._real_8))
          vas=-beta*xs*xs*(0.5_real_8*rho)**ob3/bln
       ELSE
          vas=0._real_8
       ENDIF
       vlb94 = alpha*vx + vc + vas
    ELSE
       vlb94 = 0._real_8
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE lb94
  ! ==================================================================
  SUBROUTINE lb94s(rhoa,rhob,grhoaa,grhobb,vlb94a,vlb94b)
    REAL(real_8)                             :: rhoa, rhob, grhoaa, grhobb, &
                                                vlb94a, vlb94b

    REAL(real_8), PARAMETER :: alpha = 1.19_real_8, beta = 0.01_real_8 , &
      ob3 = 1._real_8/3._real_8 , small = 1.e-14_real_8 

    INTEGER                                  :: iflg
    REAL(real_8)                             :: bln, epz, eta, ex, rho, rs, &
                                                vasa, vasb, vpza, vpzb, vxa, &
                                                vxb, xs

! ==--------------------------------------------------------------==

    rho = rhoa + rhob
    eta = (rhoa-rhob)/rho
    IF (rho.GT.small) THEN
       rs = (0.75_real_8/(pi*rho))**ob3
       CALL lsd_sx(rho,eta,ex,vxa,vxb,2._real_8/3._real_8)
       iflg=2
       IF (rs.LT.1.0_real_8) iflg=1
       CALL lsd_pz(rs,eta,epz,vpza,vpzb,iflg)
       ! 
       IF (rhoa.GT.cntr%gceps) THEN
          xs=SQRT(grhoaa)/rhoa**(4._real_8*ob3)
          bln=1._real_8+3._real_8*beta*xs*LOG(xs+SQRT(xs*xs+1._real_8))
          vasa=-beta*xs*xs*rhoa**ob3/bln
       ELSE
          vasa=0._real_8
       ENDIF
       IF (rhob.GT.cntr%gceps) THEN
          xs=SQRT(grhobb)/rhob**(4._real_8*ob3)
          bln=1._real_8+3._real_8*beta*xs*LOG(xs+SQRT(xs*xs+1._real_8))
          vasb=-beta*xs*xs*rhob**ob3/bln
       ELSE
          vasb=0._real_8
       ENDIF
       ! 
       vlb94a = alpha*vxa + vpza + vasa
       vlb94b = alpha*vxb + vpzb + vasb
    ELSE
       vlb94a = 0._real_8
       vlb94b = 0._real_8
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE lb94s
  ! ==================================================================
  SUBROUTINE lb94m(rhoe,grho,exc,vlb94)
    ! ==--------------------------------------------------------------==
    ! == CALCULATE SAOP POTENTIAL                                     ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rhoe(fpar%nnr1,clsd%nlsd), &
                                                grho(fpar%nnr1,*), exc, &
                                                vlb94(fpar%nnr1,clsd%nlsd)

    INTEGER                                  :: ir

    CALL zeroing(vlb94)!,nnr1*clsd%nlsd)
    IF (cntl%tlsd) THEN
       DO ir=1,fpar%nnr1
          CALL lb94s(rhoe(ir,1),rhoe(ir,2),grho(ir,1),&
               grho(ir,3),vlb94(ir,1),vlb94(ir,2))
       ENDDO
    ELSE
       DO ir=1,fpar%nnr1
          CALL lb94(rhoe(ir,1),grho(ir,1),vlb94(ir,1))
       ENDDO
    ENDIF
    exc=0._real_8
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE lb94m
  ! ==================================================================
  SUBROUTINE gllb(c0,eigv,eref,focc,nstate,rhoe,grho,exc,&
       vgllb,rho1,psi)
    ! ==--------------------------------------------------------------==
    ! == CALCULATE GLLB POTENTIAL                                     ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(ncpw%ngw,*)
    REAL(real_8)                             :: eigv(*), eref(*), focc(*)
    INTEGER                                  :: nstate
    REAL(real_8) :: rhoe(fpar%nnr1,clsd%nlsd), grho(fpar%nnr1,*), exc, &
      vgllb(fpar%nnr1,clsd%nlsd), rho1(fpar%nnr1)
    COMPLEX(real_8)                          :: psi(:)

    REAL(real_8), PARAMETER                  :: kresp = 0.42_real_8 

    INTEGER                                  :: ik, ir, is
    REAL(real_8)                             :: ec, en, ena, enb, ex, weni

    ena=0._real_8
    enb=0._real_8
    en=0._real_8
    IF (cntl%tlsd) THEN
       DO is=1,spin_mod%nsup
          IF (focc(is).GT.1.e-4_real_8) ena=eigv(is)
       ENDDO
       DO is=spin_mod%nsup+1,nstate
          IF (focc(is).GT.1.e-4_real_8) enb=eigv(is)
       ENDDO
    ELSE
       DO is=1,nstate
          IF (focc(is).GT.1.e-4_real_8) en=eigv(is)
       ENDDO
    ENDIF
    ! 
    ! Vresponse
    IF (group%nogrp.GT.1) CALL stopgm("GLLB","NO GROUPS ALLOWED",& 
         __LINE__,__FILE__)
    CALL zeroing(vgllb)!,nnr1*clsd%nlsd)
    DO is=1,nstate
       CALL rhosofr(c0(:,is),rho1,psi)
       IF (focc(is).GT.1.e-4_real_8) THEN
          IF (cntl%tlsd) THEN
             IF (is.LE.spin_mod%nsup) THEN
                ena=eref(is)
                weni=SQRT(ena-eigv(is))*kresp*focc(is)
                ik=1
             ELSE
                enb=eref(is)
                weni=SQRT(enb-eigv(is))*kresp*focc(is)
                ik=2
             ENDIF
          ELSE
             en=eref(is)
             weni=SQRT(en-eigv(is))*kresp*focc(is)
             ik=1
          ENDIF
          DO ir=1,fpar%nnr1
             IF (rhoe(ir,ik).GT.1.e-12_real_8) THEN
                vgllb(ir,ik) = vgllb(ir,ik) + weni * rho1(ir)/rhoe(ir,ik)
             ENDIF
          ENDDO
       ENDIF
    ENDDO
    ! PBE exchange energy density
    IF (cntl%tlsd) THEN
       DO ir=1,fpar%nnr1
          CALL pbe_xed(2._real_8*rhoe(ir,1),4._real_8*grho(ir,1),ex)
          vgllb(ir,1) = vgllb(ir,1) + 2._real_8*ex
       ENDDO
       DO ir=1,fpar%nnr1
          CALL pbe_xed(2._real_8*rhoe(ir,2),4._real_8*grho(ir,3),ex)
          vgllb(ir,2) = vgllb(ir,2) + 2._real_8*ex
       ENDDO
    ELSE
       DO ir=1,fpar%nnr1
          CALL pbe_xed(rhoe(ir,1),grho(ir,1),ex)
          vgllb(ir,1) = vgllb(ir,1) + 2._real_8*ex
       ENDDO
    ENDIF
    ! PBE correlation energy density
    IF (cntl%tlsd) THEN
       DO ir=1,fpar%nnr1
          CALL pbe_ceds(rhoe(ir,1),rhoe(ir,2),grho(ir,1),&
               grho(ir,2),grho(ir,3),ec)
          vgllb(ir,1) = vgllb(ir,1) + 2._real_8*ec
          vgllb(ir,2) = vgllb(ir,2) + 2._real_8*ec
       ENDDO
    ELSE
       DO ir=1,fpar%nnr1
          CALL pbe_ced(rhoe(ir,1),grho(ir,1),ec)
          vgllb(ir,1) = vgllb(ir,1) + 2._real_8*ec
       ENDDO
    ENDIF
    ! 
    exc = 0._real_8
    ! 
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gllb
  ! ==================================================================

END MODULE saop_utils
