MODULE v1xc_p_utils
  USE cppt,                            ONLY: nzh
  USE error_handling,                  ONLY: stopgm
  USE fft_maxfft,                      ONLY: maxfft
  USE fftmain_utils,                   ONLY: fwfftn
  USE func,                            ONLY: func1,&
                                             mfxcc_is_lyp,&
                                             mfxcc_is_pade
  USE gcener_utils,                    ONLY: gcener,&
                                             gclsd
  USE graden_utils,                    ONLY: graden
  USE kinds,                           ONLY: real_8
  USE nlcc,                            ONLY: corel
  USE parac,                           ONLY: parai
  USE reshaper,                        ONLY: reshape_inplace
  USE response_pmod,                   ONLY: ener1,&
                                             rho0,&
                                             vxc0
  USE spin,                            ONLY: clsd
  USE system,                          ONLY: cntl,&
                                             cntr,&
                                             fpar,&
                                             ncpw,&
                                             parm,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: zgthr
  USE xcener_utils,                    ONLY: xcener
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: v1xc_p
  !public :: vxcofrho
  !public :: vxc0_lda
  PUBLIC :: give_scr_v1xc_p

CONTAINS

  ! ==================================================================
  SUBROUTINE v1xc_p(drhoe,v,psi)
    ! ==--------------------------------------------------------------==
    ! CALCULATES the exchange-correlation contribution of the 
    ! perturbation density n(1) and adds it to the potential (v)
    ! and energy (exc1) terms, and ADDS the final potential function
    ! v_xc[n(1)] to the total local potential, V.
    ! INPUT:  drhoe = n1 in R space
    ! \       v     = the (spin up+down) potential in R space
    ! OUTPUT: v     = the xc-contribution is added, and the
    ! \               spin splitting is done. Output is
    ! \               v(nnr1,1) alpha spin potential
    ! \               v(nnr1,2) alpha spin potential
    ! \       exc1  = the xc-energy
    ! 
    ! The routine works in fundamentally different ways for
    ! LDA and GGA:
    ! 
    ! In LDA, we calculate once for all the expression
    ! k_xc(r,rp) = d^2 e_xc[n]  /  dn(rp) dn(r)
    ! and then, for each step, we do
    ! v1_xc(r) = Int d^3rp  k_xc(r,rp) * n1(rp).
    ! 
    ! For GGA, we compute for each step the second derivative
    ! times n1 by a finite difference technique:
    ! v1_xc(r) =  Int d^3r   k_xc(rp,r)  * n1(rp)
    ! \        =  v_xc(n0+eps*n1)(r) - v_xc(n0-eps*n1)(r) / 2 eps
    ! or even simpler:
    ! \        =  v_xc(n0+eps*n1)(r) - v_xc(n0)(r)  / eps
    ! 
    ! 
    ! The array VXC0 has two different functions in those cases.
    ! LDA/analytic calculation:  VXC0 = d^2 Exc / d n^2
    ! GGA/finite differences:    VXC0 = Vxc [ n0 ]
    ! ==--------------------------------------------------------------==
    ! Parameters:
    REAL(real_8)                             :: drhoe(fpar%nnr1,clsd%nlsd), &
                                                v(fpar%nnr1,clsd%nlsd)
    COMPLEX(real_8)                          :: psi(:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'v1xc_p'
    REAL(real_8), PARAMETER                  :: eps = 0.002_real_8 

    COMPLEX(real_8), ALLOCATABLE             :: vcscr(:), vxctemp(:,:)
    INTEGER                                  :: ierr, ir, isub, jr
    REAL(real_8)                             :: v1xc
    REAL(real_8), ALLOCATABLE                :: gradients(:,:), rho01(:,:), &
                                                vscr(:,:), vxc1(:,:)

! ==--------------------------------------------------------------==

    CALL tiset('     V1XC_P',isub)
    IF (corel%tinlc) THEN
       CALL stopgm('V1XC_P',&
            'NonLinear CoreCorrection not implemented.',& 
            __LINE__,__FILE__)
    ENDIF

    ALLOCATE(vscr(fpar%nnr1,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(vcscr(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(vxctemp(maxfft, clsd%nlsd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(vxc1(MAX(2*ncpw%nhg, fpar%nnr1), clsd%nlsd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(rho01(fpar%nnr1, clsd%nlsd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(gradients(fpar%nnr1, 4*clsd%nlsd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    ! GRADIENTS contains 4 densities for each spin.

    ! ==--------------------------------------------------------------==
    DO ir = 1, fpar%nnr1
       DO jr = 1, clsd%nlsd
          rho01(ir,jr) = MAX(0._real_8, rho0(ir,jr) + eps * drhoe(ir,jr))
       ENDDO               ! prevent rho from being negative
    ENDDO
    CALL vxcofrho(rho01, vxc1, vxctemp, vscr, vcscr, psi, gradients)

    DO ir = 1, fpar%nnr1
       DO jr = 1, clsd%nlsd
          rho01(ir,jr) = MAX(0._real_8, rho0(ir,jr) - eps * drhoe(ir,jr))
       ENDDO
    ENDDO
    CALL vxcofrho(rho01, vxc0, vxctemp, vscr, vcscr, psi, gradients)

    ener1%exc1 = 0._real_8
    IF (.NOT. cntl%tlsd) THEN
       DO ir=1,fpar%nnr1
          IF (rho0(ir,1) .GT. cntr%gceps) THEN
             v1xc =  (vxc1(ir,1)-vxc0(ir,1)) / (2._real_8*eps)
             v(ir,1) = v(ir,1) + v1xc
             ener1%exc1 = ener1%exc1 + v1xc*drhoe(ir,1)
          ENDIF
       ENDDO
    ELSE                ! with spin split:


       DO ir=1,fpar%nnr1
          IF (rho0(ir,1) .GT. cntr%gceps) THEN
             ! beta spin:
             v1xc = (vxc1(ir,2)-vxc0(ir,2)) / (2._real_8*eps)
             v(ir,2) = v(ir,1) + v1xc
             ener1%exc1 = ener1%exc1 + v1xc*drhoe(ir,2)
             ! alpha spin:
             v1xc = (vxc1(ir,1)-vxc0(ir,1)) / (2._real_8*eps)
             v(ir,1) = v(ir,1) + v1xc
             ener1%exc1 = ener1%exc1 + v1xc*(drhoe(ir,1)-drhoe(ir,2))
          ELSE
             v(ir,2) = v(ir,1)
          ENDIF
       ENDDO
    ENDIF
    ener1%exc1    = 0.5_real_8 * ener1%exc1*parm%omega/REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)

    DEALLOCATE(vscr,vcscr,vxctemp,vxc1,rho01,gradients,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)

    ! ==--------------------------------------------------------------==
    CALL tihalt('     V1XC_P',isub)
    RETURN
  END SUBROUTINE v1xc_p
  ! ==================================================================


  ! ==================================================================
  SUBROUTINE vxcofrho(rho, Vxc, vxctemp, vscr, vcscr, psi, gradients)
    ! ==--------------------------------------------------------------==
    ! INPUT: rho in R-space, (alpha+beta, beta) storage (untouched)
    ! OUTPUT: vxc in R-space, (alpha, beta) storage
    ! Rest is scratch.

    REAL(real_8)                             :: rho(fpar%nnr1,clsd%nlsd)
    REAL(real_8), TARGET                     :: Vxc(fpar%nnr1,clsd%nlsd)
    COMPLEX(real_8)                          :: vxctemp(maxfft,clsd%nlsd)
    REAL(real_8)                             :: vscr(:,:)
    COMPLEX(real_8)                          :: vcscr(:), psi(:)
    REAL(real_8) :: gradients(fpar%nnr1,4*clsd%nlsd)

    COMPLEX(real_8), POINTER                 :: vxc_g(:,:)
    REAL(real_8)                             :: sgcc, sgcx, ssxc, vvxc

! TODO refactor this array 
! ==--------------------------------------------------------------==

    IF (cntl%ttau) THEN
       CALL stopgm('vxcofrho','META FUNCTIONALS NOT IMPLENTED',& 
            __LINE__,__FILE__)
    ENDIF

    CALL reshape_inplace(Vxc, (/ncpw%nhg, clsd%nlsd/), vxc_g)

    CALL zeroing(vxctemp)!, maxfft*clsd%nlsd)
    CALL zeroing(Vxc)!,nnr1*clsd%nlsd)

    CALL xcener(ssxc,vvxc,rho,rho,vxctemp)
    ! The non-grad-corrected potential is now in VXCtemp
    ! (in R-space, in complex(8) :: storage)

    IF (cntl%tgc) THEN
       CALL  fwfftn(vxctemp(:,1),.FALSE.,parai%allgrp)
       CALL zgthr(ncpw%nhg,vxctemp(:,1),vxc_g(:,1),nzh)
       ! The non-grad-corrected potential is now in VXC 
       ! (in G-space, packed storage)

       IF (cntl%tlsd) THEN
          CALL  fwfftn(vxctemp(:,2),.FALSE.,parai%allgrp)
          CALL zgthr(ncpw%nhg,vxctemp(:,2),vxc_g(:,2),nzh)

          CALL daxpy(fpar%nnr1,-1._real_8,rho(1,2),1, rho(1,1),1)
          CALL graden(rho(1,1),psi,gradients(1,1),vcscr)
          CALL graden(rho(1,2),psi,gradients(1,5),vcscr)
          CALL gclsd(sgcx,sgcc,rho,vxctemp,vxc_g,vscr,gradients,.FALSE.)
          CALL daxpy(fpar%nnr1,+1._real_8,rho(1,2),1, rho(1,1),1)
       ELSE
          CALL graden(rho(1,1),psi,gradients(1,1),vcscr)
          CALL gcener(sgcx,sgcc,rho,vxctemp,vxc_g,vscr(:,1),gradients,&
               .FALSE.)
       ENDIF
    ENDIF                     ! if cntl%tgc = if gradient correction


    ! XCENER and GCener/GCLSD return the total, gradient-corrected
    ! xc-potential in R-space (complex(8) :: storage) in vxctmp.
    CALL dcopy(fpar%nnr1,vxctemp(1,1),2,Vxc,1)
    IF (cntl%tlsd) CALL dcopy(fpar%nnr1,vxctemp(1,2),2,Vxc(1,2),1)

    RETURN
  END SUBROUTINE vxcofrho
  ! ==================================================================







  ! ==================================================================
  SUBROUTINE vxc0_lda(vxc0,rho0)
    ! ==------------------------------------------------------------==
    ! calculates the 2nd functional derivative of the exchange-
    ! correlation potential used in perturbation theory,
    ! FOR LDA.
    ! INPUT:  rho0 = the (nnr1) dimension density array containing
    ! \              the gnd state density (direct space, real(8) :: )
    ! OUTPUT: vxc0 = d^2 E_xc [n0] / dn^2(r)
    ! ==------------------------------------------------------------==
    ! Parameters:
    REAL(real_8)                             :: vxc0(*), rho0(*)

    REAL(real_8), PARAMETER :: a0u = 0.4581652932831429_real_8, &
      a1u = 2.217058676663745_real_8, a2u = 0.7405551735357053_real_8, &
      a3u = 0.01968227878617998_real_8, b1u = 1.0_real_8, &
      b2u = 4.504130959426697_real_8, b3u = 1.110667363742916_real_8, &
      b4u = 0.02359291751427506_real_8, coef = 0.620350490901_real_8 

    INTEGER                                  :: ir
    REAL(real_8)                             :: dexcdrs, dexcdrs2, dmu, &
                                                drsdn, drsdn2, rho, root_rho, &
                                                rs

    IF (func1%mfxcc /= mfxcc_is_pade.AND.func1%mfxcc /= mfxcc_is_lyp) THEN
       CALL stopgm('v0xc_lda','ONLY LDA',& 
            __LINE__,__FILE__)
    ELSE
       DO ir=1,fpar%nnr1
          rho = rho0(ir)
          IF (rho .LT. cntr%gceps) THEN
             vxc0(ir) = 0._real_8
          ELSE
             ! ..   Pade approximation to LDA functional
             root_rho=rho**(1._real_8/3._real_8)
             rs=coef/root_rho
             drsdn=-coef/(3._real_8*root_rho**4)
             drsdn2=(4._real_8*coef)/(9._real_8*root_rho**7)
             dexcdrs=(a0u*(b1u +&
                  rs*(2*b2u + 3*b3u*rs + 4*b4u*rs**2)) +&
                  rs**2*(a1u*b2u - 2*a3u*b1u*rs + 2*a1u*b3u*rs -&
                  a3u*b2u*rs**2 +&
                  3*a1u*b4u*rs**2 + a3u*b4u*rs**4 +&
                  a2u*(-b1u + b3u*rs**2 + 2*b4u*rs**3)))/&
                  (rs**2*(b1u + rs*(b2u + rs*(b3u + b4u*rs)))**2)
             dexcdrs2= (2*rs*(b2u + 3*rs*(b3u + 2*b4u*rs))*&
                  (a0u + rs*(a1u + rs*(a2u + a3u*rs)))*&
                  (b1u + rs*(b2u + rs*(b3u + b4u*rs))) -&
                  2*rs**2*(a2u + 3*a3u*rs)*&
                  (b1u + rs*(b2u + rs*(b3u + b4u*rs)))**2 +&
                  2*rs*(a1u + rs*(2*a2u + 3*a3u*rs))*&
                  (b1u + rs*(b2u + rs*(b3u + b4u*rs)))*&
                  (b1u + rs*(2*b2u + rs*(3*b3u + 4*b4u*rs))) -&
                  2*(a0u + rs*(a1u + rs*(a2u + a3u*rs)))*&
                  (b1u + rs*(2*b2u + rs*(3*b3u + 4*b4u*rs)))**2)/&
                  (rs**3*(b1u + rs*(b2u + rs*(b3u + b4u*rs)))**3)
             dmu=2._real_8*drsdn+rho*drsdn2

             vxc0(ir)=dexcdrs*dmu+rho*dexcdrs2*drsdn*drsdn

          ENDIF
       ENDDO
    ENDIF
    RETURN
  END SUBROUTINE vxc0_lda
  ! ==--------------------------------------------------------------



  ! ==================================================================
  SUBROUTINE give_scr_v1xc_p(lvofrho,tag)
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: lvofrho
    CHARACTER(len=30)                        :: tag

    lvofrho  =  3*fpar%nnr1 + 2*maxfft*clsd%nlsd + MAX(2*ncpw%nhg*clsd%nlsd,fpar%nnr1*clsd%nlsd)+&
         fpar%nnr1*clsd%nlsd + 4*fpar%nnr1*clsd%nlsd

    tag      =  'v1xc_p      '
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_v1xc_p
  ! ==--------------------------------------------------------------==

END MODULE v1xc_p_utils
