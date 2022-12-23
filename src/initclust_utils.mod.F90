MODULE initclust_utils
  USE cnst,                            ONLY: fpi,&
                                             pi
  USE cppt,                            ONLY: hg,&
                                             nr1h,&
                                             nr2h,&
                                             nr3h,&
                                             nr3pl,&
                                             scg,&
                                             scgx
  USE cp_xc_utils,                     ONLY: cp_xc_functional, &
                                             cp_xc_kernel, &
                                             cp_xc_mts_low_func
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE fft,                             ONLY: llr1
  USE fftnew_utils,                    ONLY: addfftnset,&
                                             setfftn
  USE func,                            ONLY: ashcroft_coulomb_rcut,&
                                             func1,&
                                             func2,&
                                             msrx_is_skipped,&
                                             msrx_is_CAM,&
                                             msrx_is_ashcroft,&
                                             msrx_is_erfc,&
                                             msrx_is_exp
  USE geq0mod,                         ONLY: geq0
  USE hfxmod,                          ONLY: hfxc2,&
                                             ipoolhfx
  USE hipin_utils,                     ONLY: hipin
  USE isos,                            ONLY: isos1,&
                                             isos3
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE mp_interface,                    ONLY: mp_max,&
                                             mp_sum
  USE mtin_utils,                      ONLY: moin,&
                                             mtin
  USE parac,                           ONLY: parai,&
                                             paral
  USE rswfmod,                         ONLY: lwdimx,&
                                             maxstatesx,&
                                             rswfx
  USE scex_utils,                      ONLY: scex
  USE spin,                            ONLY: lspin2
  USE system,                          ONLY: cntl,&
                                             cntr,&
                                             group,&
                                             ncpw,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: initclust
  PUBLIC :: gf_periodic
  PUBLIC :: hf_init

CONTAINS

  ! ==================================================================
  SUBROUTINE initclust
    ! ==--------------------------------------------------------------==
    ! == INITIALIZATION ROUTINE FOR CLUSTER BOUNDARY CONDITIONS       ==
    ! ==--------------------------------------------------------------==
    CHARACTER(*), PARAMETER                  :: procedureN = 'initclust'
    LOGICAL, PARAMETER                       :: init_is_response = .FALSE.

    INTEGER                                  :: ierr, isub

! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    IF (isos1%tclust) THEN
       IF (isos3%ps_type == 1) THEN
          CALL hipin
       ELSEIF (isos3%ps_type == 2) THEN
          CALL mtin
       ELSEIF (isos3%ps_type == 3) THEN
          CALL moin
       ENDIF
    ELSE
       nr1h=1
       nr2h=1
       nr3h=1
       nr3pl=1
       ALLOCATE(scg(ncpw%nhg),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL gf_periodic
    ENDIF
    CALL hf_init(init_is_response)
    CALL tihalt(procedureN,isub)

  END SUBROUTINE initclust
  ! ==================================================================
  SUBROUTINE gf_periodic
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ig, ig1
    REAL(real_8)                             :: g2

! ==--------------------------------------------------------------==

    ig1=1
    IF (geq0) ig1=2
    DO ig=ig1,ncpw%nhg
       g2=parm%tpiba2*hg(ig)
       scg(ig)=CMPLX(fpi/g2,0.0_real_8,kind=real_8)
    ENDDO
    IF (geq0) scg(1)=0.0_real_8

  END SUBROUTINE gf_periodic
  ! ==================================================================
  SUBROUTINE hf_init(response)
    ! ==--------------------------------------------------------------==
    LOGICAL                                  :: response

    CHARACTER(*), PARAMETER                  :: procedureN = 'hf_init'

    INTEGER                                  :: ierr, ldim, nstate
    INTEGER                                  :: functional_msrx, functional_mhfx
    INTEGER                                  :: kernel_msrx, kernel_mhfx
    INTEGER                                  :: mts_low_func_msrx, mts_low_func_mhfx
    INTEGER                                  :: mhfx, msrx
    INTEGER, SAVE                            :: ifirst = 0
    REAL(real_8)                             :: rmem, rstate
    REAL(real_8)                             :: functional_srxa
    REAL(real_8)                             :: kernel_srxa
    REAL(real_8)                             :: mts_low_func_srxa
    REAL(real_8)                             :: srxa, cam_alpha, cam_beta

    !
    ! 31.08.17 M.P.Bircher@LCBC/EPFL
    ! Manually test compatiblity of HFX in cp_xc_kernel and _functional
    ! This can be removed once the refactoring is finished
    !
    ! func1%mhfx etc. already hold the info from cp_xc_functional, now we 
    ! need the kernel
    IF ( cntl%use_xc_driver ) THEN
       CALL cp_xc_functional%get(mhfx=functional_mhfx,msrx=functional_msrx,srxa=functional_srxa)
       CALL cp_xc_kernel%get(mhfx=kernel_mhfx,msrx=kernel_msrx,srxa=kernel_srxa)
       CALL cp_xc_mts_low_func%get(mhfx=mts_low_func_mhfx,msrx=mts_low_func_msrx,srxa=mts_low_func_srxa)
       !
       ! Do we use HFX both times? Then, we cannot mix CAM and HFX for the moment
       ! (will fix that later)
       !
       IF (kernel_mhfx == functional_mhfx) THEN
          IF (functional_msrx /= kernel_msrx) CALL stopgm(procedureN,&
             'LR_KERNEL and FUNCTIONAL must imperatively use the same HFX screening function',&
             __LINE__,__FILE__)
          IF (functional_srxa /= kernel_srxa) CALL stopgm(procedureN,&
             'LR_KERNEL and FUNCTIONAL must imperatively use the same HFX screening parameter',&
             __LINE__,__FILE__)
          CALL cp_xc_functional%get(mhfx=mhfx,msrx=msrx,srxa=srxa,cam_alpha=cam_alpha,cam_beta=cam_beta)
       ELSEIF (kernel_mhfx > 0 .and. functional_mhfx == 0) THEN
          ! No HFX in functional, but in kernel, thus get the parameters from the kernel
          CALL cp_xc_kernel%get(mhfx=mhfx,msrx=msrx,srxa=srxa,cam_alpha=cam_alpha,cam_beta=cam_beta)
       ELSEIF (kernel_mhfx == 0 .and. functional_mhfx > 0) THEN
          ! No HFX in kernel, but in the functional, thus get the parameters from the functional
          CALL cp_xc_functional%get(mhfx=mhfx,msrx=msrx,srxa=srxa,cam_alpha=cam_alpha,cam_beta=cam_beta)
       ELSE
          CALL stopgm(procedureN,'Invalid HFX combination between kernel and functional',&
                      __LINE__,__FILE__)
       ENDIF

       ! Martin:
       ! MTS type is NOT copied from the standard functional if it is not used.
       ! Need an IF to skip bugs with screened functionals
       !
       IF (cntl%use_mts) THEN
          ! 
          ! same thing for MTS low functional
          !
          IF (mts_low_func_mhfx == functional_mhfx) THEN
             IF (functional_msrx /= mts_low_func_msrx) CALL stopgm(procedureN,&
                'MTS: LOW and HIGH level FUNCTIONALS must imperatively use the same HFX screening function',&
                __LINE__,__FILE__)
             IF (functional_srxa /= mts_low_func_srxa) CALL stopgm(procedureN,&
                'MTS: LOW and HIGH level FUNCTIONALS must imperatively use the same HFX screening parameter',&
                __LINE__,__FILE__)
             CALL cp_xc_functional%get(mhfx=mhfx,msrx=msrx,srxa=srxa,cam_alpha=cam_alpha,cam_beta=cam_beta)
          ELSEIF (mts_low_func_mhfx > 0 .and. functional_mhfx == 0) THEN
             ! No HFX in functional, but in MTS low func., thus get the parameters from the MTS low func.
             CALL cp_xc_mts_low_func%get(mhfx=mhfx,msrx=msrx,srxa=srxa,cam_alpha=cam_alpha,cam_beta=cam_beta)
          ELSEIF (mts_low_func_mhfx == 0 .and. functional_mhfx > 0) THEN
             ! No HFX in MTS low func., but in the main functional, thus get the parameters from the main functional
             CALL cp_xc_functional%get(mhfx=mhfx,msrx=msrx,srxa=srxa,cam_alpha=cam_alpha,cam_beta=cam_beta)
          ELSE
             CALL stopgm(procedureN,'MTS: Invalid HFX combination between secomdary and primary functionals',&
                         __LINE__,__FILE__)
          ENDIF
       ENDIF
    ELSE
       mhfx      = func1%mhfx
       msrx      = func1%msrx
       srxa      = func2%srxa
       cam_alpha = func2%cam_alpha      
       cam_beta  = func2%cam_beta      
    ENDIF
    !
    ! End of hack

    IF ( mhfx /= 0 .AND. ifirst == 0) THEN
       IF (tkpts%tkpnt) CALL stopgm(procedureN,'K-POINTS NOT IMPLEMENTED',& 
            __LINE__,__FILE__)
       IF (lspin2%tlse) CALL stopgm(procedureN,'NO LSE POSSIBLE',& 
            __LINE__,__FILE__)
       IF (group%nogrp > 1) CALL stopgm(procedureN,'TASK GROUPS NOT IMPLEMENTED',& 
            __LINE__,__FILE__)

       scgx => scg
       !
       IF (msrx == msrx_is_exp) THEN
          ! exp(-a*r)/r
          CALL get_scgx_for_exponential(srxa)
          !
       ELSEIF (msrx == msrx_is_erfc) THEN
          ! erfc(a*r)/r
          CALL get_scgx_for_erfc(srxa)
          !
       ELSEIF (msrx == msrx_is_CAM) THEN
          ! erf(a*r)/r
          IF (cntl%div_analytical) THEN
             CALL get_scgx_for_erf_analytical(srxa,cam_alpha,cam_beta)
          ELSE
             CALL get_scgx_for_erf_numerical(srxa,cam_alpha,cam_beta)
          END IF
          !
       ELSEIF (response) THEN
          ! ..nothing to do
       ELSEIF ( msrx  ==  msrx_is_ashcroft ) THEN
          ! ..1/r for r<R else 0
          CALL get_scgx_for_ashcroft()
          !
       ELSEIF (.NOT.isos1%tclust) THEN
          !
          CALL get_scgx_for_periodic_coulomb()
          !
       ELSEIF (isos1%tclust .AND. isos3%ps_type == 1) THEN
          CALL stopgm("HF_INIT", "HFX not implemented for HOCKNEY.",& 
               __LINE__,__FILE__)
       ENDIF
       IF (cntl%use_scaled_hfx) THEN
          IF (msrx /= msrx_is_skipped) CALL stopgm(procedureN,'ScEX: Matrix elements for screened '// &
                                                   'exchange not yet available',&
                                                   __LINE__,__FILE__)
          CALL scex%grid_init()
       ELSEIF (hfxc2%hfxwfe > 0._real_8 .AND. hfxc2%hfxdee > 0._real_8) THEN
          IF (paral%io_parent) THEN
             WRITE(6,'(/,A,T56,F10.0)')&
                  ' Wavefunction Cutoff for HFX [Ry] ',hfxc2%hfxwfe
             WRITE(6,'(A,T56,F10.0)')&
                  ' Density Cutoff for HFX [Ry] ',hfxc2%hfxdee
          ENDIF
          CALL addfftnset(hfxc2%hfxdee,hfxc2%hfxwfe,ipoolhfx)
       ELSE
          CALL addfftnset(-1._real_8,-1._real_8,ipoolhfx)
          ipoolhfx=0
       ENDIF
       ! Array to keep the real space wavefunctions
       IF (cntl%krwfn) THEN
          CALL setfftn(ipoolhfx)
          lwdimx = llr1
          nstate = crge%n
          IF (cntr%memsize.LT.0) THEN
             maxstatesx = nstate
             rstate=1._real_8
          ELSE
             rmem = 0.125_real_8*cntr%memsize*1.e6_real_8/REAL(lwdimx,kind=real_8)
             maxstatesx = ((INT(rmem)+1)/2)*2
             maxstatesx = MIN(nstate,maxstatesx)
             rmem = -ABS(maxstatesx)
             CALL mp_max(rmem,parai%allgrp)
             maxstatesx = NINT(-rmem)
             rstate= REAL(maxstatesx,kind=real_8)/REAL(nstate,kind=real_8)
          ENDIF
          ldim  = (maxstatesx+1)/2 * lwdimx
          ALLOCATE(rswfx(lwdimx,(maxstatesx+1)/2),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          rmem = 16._real_8*ldim*1.e-6_real_8
          IF (paral%io_parent) THEN
             IF (rstate > 0.99999999_real_8) THEN
                WRITE(6,'(A,T59,A)')&
                     ' HFX| ORBITALS KEPT IN REAL SPACE ','    ALL'
             ELSE
                WRITE(6,'(A,T50,F8.1,A)')&
                     ' HFX| ORBITALS KEPT IN REAL SPACE ',&
                     100*rstate,' PERCENT'
             ENDIF
             WRITE(6,'(A,T51,F8.3,A)')&
                  ' HFX| WAVEFUNCTION TAKES IN REAL SPACE ',rmem,' MBYTES'
          ENDIF
          CALL setfftn(0)
       ENDIF
    ENDIF
    ifirst = 1

  CONTAINS

    ! ==--------------------------------------------------------------==
    SUBROUTINE get_scgx_for_exponential(srxa)

    INTEGER                                  :: ierr, ig, ig1
    REAL(real_8), INTENT(in)                 :: srxa
    REAL(real_8)                             :: div, g2, srxa2

      ALLOCATE(scgx(ncpw%nhg),STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
           __LINE__,__FILE__)
      srxa2=srxa*srxa
      ! PUB START
      div = get_geq0_exp_over_r(srxa2)
      ! PUB END
      ig1=1
      IF (geq0) ig1=2
      DO ig=1,ncpw%nhg
         g2=parm%tpiba2*hg(ig) + srxa2
         scgx(ig)=CMPLX(fpi/g2,0._real_8,kind=real_8)
      ENDDO
      IF (geq0) scgx(1)= -div/2
      IF ((geq0).AND.paral%io_parent)&
           WRITE(6,*) 'PUB Divergence', REAL(scgx(1))

    END SUBROUTINE get_scgx_for_exponential
    ! ==--------------------------------------------------------------==
    SUBROUTINE get_scgx_for_erf_analytical(srxa,cam_alpha,cam_beta)
      ! 09.09.2016 M.P. Bircher @ LCBC/EPFL


    INTEGER                                  :: ierr, ig, ig1
    REAL(real_8), INTENT(in)                 :: srxa, cam_alpha, cam_beta
    REAL(real_8) :: div, g2, screening_gamma, srxa2, work_div_integral, &
      work_div_sum, work_exp_srxa2, work_g2, work_tmp

! Comments follow Phys. Rev. B 80. 085114, (2009)
!

      screening_gamma = 10.0_real_8 / cntr%ecut
      srxa2=1._real_8/(4._real_8*srxa*srxa)

      ALLOCATE(scgx(ncpw%nhg),STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
           __LINE__,__FILE__)

      work_div_sum = 0.0_real_8
      work_div_integral = 0.0_real_8
      !
      ig1=1
      IF (geq0) ig1=2
      DO ig=ig1,ncpw%nhg
         !
         ! Adapted from eqs (19), (20), (21), (23)
         !
         g2 = parm%tpiba2*hg(ig)
         work_g2 = fpi/g2
         work_exp_srxa2 = cam_alpha+cam_beta*EXP(-g2*srxa2)
         scgx(ig)=CMPLX(work_g2*work_exp_srxa2,0.0_real_8,kind=real_8)
         IF (hg(ig) >= 1e-8_real_8) THEN
            work_div_sum = work_div_sum + work_g2 * EXP(-screening_gamma*g2)*( work_exp_srxa2 )
         ELSE
            work_div_sum = work_div_sum + (work_g2 - fpi*screening_gamma)*(cam_alpha+cam_beta*(1.0_real_8-g2*srxa2))
         ENDIF
      ENDDO
      !
      work_tmp = 2.0_real_8*work_div_sum
      CALL mp_sum(work_tmp, work_div_sum, parai%allgrp)
      !
      IF (srxa2 >= 1e-8_real_8) THEN
         work_div_sum = work_div_sum - fpi*cam_alpha*screening_gamma
      ELSE
         work_div_sum = work_div_sum - fpi*cam_alpha/(srxa2)
      END IF
      !
      work_div_integral = cam_alpha/SQRT(screening_gamma*0.25_real_8*fpi) &
           + cam_beta/SQRT((screening_gamma+srxa2)*0.25_real_8*fpi)
      work_div_integral = parm%omega*work_div_integral
      !
      ! Calculate the difference between \int f(Q) and \sum f(G), eq. (10)
      !
      div =  work_div_integral - work_div_sum
      !
      IF(geq0) scgx(1)= CMPLX(div,0.0_real_8,kind=real_8)
      !
      IF (paral%io_parent)&
           WRITE(6,'(1x,A41,5x,ES18.8)') 'DIVERGENCE (ERF-COULOMB OPERATOR) AT G=0:', div

    END SUBROUTINE get_scgx_for_erf_analytical
    ! ==--------------------------------------------------------------==
    SUBROUTINE get_scgx_for_erf_numerical(srxa,cam_alpha,cam_beta)


    INTEGER                                  :: ierr, ig, ig1
    REAL(real_8), INTENT(in)                 :: srxa, cam_alpha, cam_beta
    REAL(real_8)                             :: div, div1, div2, g2, srxa2

      ALLOCATE(scgx(ncpw%nhg),STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
           __LINE__,__FILE__)

      !START PUB for the 1/r compunents
      div = get_geq0_1_over_r()
      div1 = div * (cam_alpha + cam_beta)
      ! PUB END
      srxa2=1._real_8/(4._real_8*srxa*srxa)
      ! PUB START for the erfc(gamma r)/r component
      div = get_geq0_erfc(srxa2)
      div2= div * (-1._real_8) * cam_beta
      ! PUB END

      ig1=1
      IF (geq0) ig1=2
      DO ig=ig1,ncpw%nhg
         g2=parm%tpiba2*hg(ig)
         !scgx(ig)=func2%cam_beta*CMPLX(fpi/g2,0._real_8,kind=real_8)*(EXP(-g2*srxa2)) + &
         !         func2%cam_alpha*CMPLX(fpi/g2,0._real_8,kind=real_8)
         scgx(ig)=- cam_beta*CMPLX(fpi/g2,0._real_8,kind=real_8)*(1.0_real_8-EXP(-g2*srxa2)) + &
              (cam_alpha+cam_beta)*CMPLX(fpi/g2,0._real_8,kind=real_8)
      ENDDO
      IF (geq0) scgx(1)= -(div1+div2)/2_real_8 !-div/2
      ! IF ((geq0).AND.paral%io_parent)& ! this prints only if io_parent AND geq0 coincide!
      IF (paral%io_parent)&
           WRITE(6,*) 'ivano ---PUB Divergence SRXA2', -(div1+div2)/2_real_8

    END SUBROUTINE get_scgx_for_erf_numerical
    ! ==--------------------------------------------------------------==
    SUBROUTINE get_scgx_for_erfc(srxa)


    INTEGER                                  :: ierr, ig, ig1
    REAL(real_8), INTENT(in)                 :: srxa
    REAL(real_8)                             :: div, g2, srxa2

      ALLOCATE(scgx(ncpw%nhg),STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
           __LINE__,__FILE__)
      srxa2=1._real_8/(4._real_8*srxa*srxa)

      div = get_geq0_erfc(srxa2)

      ig1=1
      IF (geq0) ig1=2
      DO ig=ig1,ncpw%nhg
         g2=parm%tpiba2*hg(ig)
         scgx(ig)=CMPLX(fpi/g2,0._real_8,kind=real_8)*(1._real_8-EXP(-g2*srxa2))
      ENDDO
      IF (geq0) scgx(1)= -div/2
      IF ((geq0).AND.paral%io_parent)&
           WRITE(6,*) 'PUB Divergence SRXA2', scgx(1)

    END SUBROUTINE get_scgx_for_erfc
    ! ==--------------------------------------------------------------==
    SUBROUTINE get_scgx_for_ashcroft()


    INTEGER                                  :: ierr, ig, ig1
    REAL(real_8)                             :: g1, g2

      ALLOCATE(scgx(ncpw%nhg),STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
           __LINE__,__FILE__)

      ig1=1
      IF (geq0) ig1=2
      DO ig=ig1,ncpw%nhg
         g2=parm%tpiba2*hg(ig)
         g1=SQRT(g2)
         scgx(ig)=CMPLX(fpi/g2,0._real_8,kind=real_8)*(1.0_real_8-COS(g1*ashcroft_coulomb_rcut))
      ENDDO
      IF (geq0) scgx(1)=CMPLX(2.0_real_8*pi*ashcroft_coulomb_rcut**2,0._real_8,kind=real_8)

    END SUBROUTINE get_scgx_for_ashcroft
    ! ==--------------------------------------------------------------==
    SUBROUTINE get_scgx_for_periodic_coulomb()


    INTEGER                                  :: ierr, ig, ig1
    REAL(real_8)                             :: div, g2

      ALLOCATE(scgx(ncpw%nhg),STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
           __LINE__,__FILE__)
      ig1=1
      IF (geq0) ig1=2
      DO ig=ig1,ncpw%nhg
         g2=parm%tpiba2*hg(ig)
         scgx(ig)=CMPLX(fpi/g2,0._real_8,kind=real_8)
      ENDDO

      div = get_geq0_1_over_r()
      IF (geq0) scgx(1)= -0.5_real_8*div

    END SUBROUTINE get_scgx_for_periodic_coulomb
    ! ==--------------------------------------------------------------==

  END SUBROUTINE hf_init
  ! ==================================================================
  ! ==================================================================
  FUNCTION get_geq0_exp_over_r(srxa2)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: srxa2, get_geq0_exp_over_r

    INTEGER                                  :: ig, iq, nqq
    REAL(real_8)                             :: aa, alpha, div, dq, q, qq, &
                                                tmp_div

! ==--------------------------------------------------------------==

    alpha = 10*parm%tpiba2/cntr%ecut
    div = 0.0_real_8
    DO ig =1, ncpw%nhg
       IF (hg(ig) > 1e-8)THEN
          div= div + EXP(-alpha*hg(ig))/(hg(ig) + srxa2/parm%tpiba2)
       ELSE
          div = div - alpha
       ENDIF
    ENDDO
    tmp_div = div
    CALL mp_sum(tmp_div, div, parai%allgrp)
    div = 2*div + alpha
    div = div * 2 * fpi / parm%tpiba2
    alpha = alpha/parm%tpiba2
    aa = 0.0_real_8
    nqq = 100000
    dq = 5.0/SQRT(alpha)/nqq
    DO iq = 0,nqq
       q = dq * (iq+0.5_real_8)
       qq = q * q
       aa = aa - EXP( -alpha * qq) * srxa2/ (qq + srxa2) * dq
    ENDDO
    aa = aa * 8._real_8/fpi
    aa = aa + 1._real_8/SQRT(alpha*0.25_real_8*fpi)
    div = div - 2*parm%omega * aa
    get_geq0_exp_over_r = div
    ! ==--------------------------------------------------------------==
  END FUNCTION get_geq0_exp_over_r
  ! ==================================================================
  FUNCTION get_geq0_erfc(srxa2)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: srxa2, get_geq0_erfc

    INTEGER                                  :: ig, iq, nqq
    REAL(real_8)                             :: aa, alpha, div, dq, q, qq, &
                                                tmp_div

! ==--------------------------------------------------------------==

    alpha = 10*parm%tpiba2/cntr%ecut
    div = 0
    DO ig =1, ncpw%nhg
       IF (hg(ig) > 1e-8)THEN
          div= div + EXP(-alpha*hg(ig))/hg(ig)*&
               (1-EXP(-hg(ig)*srxa2*parm%tpiba2))
       ELSE
          div = div - alpha
       ENDIF
    ENDDO
    tmp_div = div
    CALL mp_sum(tmp_div, div, parai%allgrp)
    div = 2*div
    IF (srxa2 > 1e-8) THEN
       div = div + alpha
    ELSE
       div = div - parm%tpiba2/(-srxa2)
    ENDIF
    div = div * 2 * fpi / parm%tpiba2
    alpha = alpha/parm%tpiba2
    aa = 0
    nqq = 100000
    dq = 5.0/SQRT(alpha)/nqq
    DO iq = 0,nqq
       q = dq * (iq+0.5_real_8)
       qq = q * q
       aa = aa-EXP(-alpha*qq)*EXP(-qq*srxa2)*dq
    ENDDO
    aa = aa * 8._real_8/fpi
    aa = aa + 1._real_8/SQRT(alpha*0.25_real_8*fpi)
    IF (paral%io_parent)&
         WRITE(6,*) 'PUB AA', alpha, srxa2
    div = div - 2*parm%omega * aa
    get_geq0_erfc = div
    ! ==--------------------------------------------------------------==
  END FUNCTION get_geq0_erfc
  ! ==================================================================
  FUNCTION get_geq0_1_over_r()
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: get_geq0_1_over_r

    INTEGER                                  :: ig
    REAL(real_8)                             :: aa, alpha, div, tmp_div

! ==--------------------------------------------------------------==

    alpha = 10*parm%tpiba2/cntr%ecut
    div = 0
    DO ig = 1, ncpw%nhg
       IF (hg(ig) > 1e-8) THEN
          div= div + EXP(-alpha*hg(ig))/hg(ig)
       ELSE
          div = div - alpha
       ENDIF
    ENDDO
    tmp_div = div
    CALL mp_sum(tmp_div, div, parai%allgrp)
    div = 2*div
    div = div + alpha
    alpha = alpha/parm%tpiba2
    div = div* 2 * fpi/parm%tpiba2
    aa = 1._real_8/SQRT(alpha*0.25_real_8*fpi)
    div = div - 2*parm%omega * aa
    get_geq0_1_over_r = div
    ! ==--------------------------------------------------------------==
  END FUNCTION get_geq0_1_over_r
  ! ==================================================================
END MODULE initclust_utils
