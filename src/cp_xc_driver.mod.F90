MODULE cp_xc_driver
  USE cp_xc_utils,                     ONLY: cp_xc,&
                                             cp_xc_functional_p_t,&
                                             cp_xc_scratch_p_t
  USE cpfunc_types,                    ONLY: cp_xc_functional_t,&
                                             cp_xc_scratch_t,&
                                             cp_xc_spin_components,&
                                             cp_xc_spin_pairs
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE lxc_utils,                       ONLY: XC_FAMILY_GGA,&
                                             XC_FAMILY_HYB_GGA,&
                                             XC_FAMILY_HYB_MGGA,&
                                             XC_FAMILY_LDA,&
                                             XC_FAMILY_MGGA,&
                                             lxc_func_info_get_family,&
                                             lxc_gga_exc_vxc,&
                                             lxc_lda_exc_vxc
  USE parac,                           ONLY: paral
  USE system,                          ONLY: cntl,&
                                             cntr,&
                                             dual00
  USE timer,                           ONLY: tiset,&
                                             tihalt

  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER      :: a        = 1
  INTEGER, PARAMETER      :: b        = 2
  INTEGER, PARAMETER      :: ab       = 3
  INTEGER, PARAMETER      :: s        = a

  PUBLIC :: cp_xc_compute

CONTAINS

  ! ==================================================================
  SUBROUTINE cp_xc_compute( n, rhoe, grad, tau, sgcx, sgcc, v1, v2, vt )
    ! ==--------------------------------------------------------------==
    ! Wrapper for libxc & internal functionals (CPfunc)    
    ! ==--------------------------------------------------------------==

    INTEGER, INTENT(IN)                      :: n
    REAL(real_8), DIMENSION(:, :), &
         INTENT(IN)                          :: rhoe, grad, tau
    REAL(real_8), INTENT(INOUT)              :: sgcx, sgcc
    COMPLEX(real_8), DIMENSION(:, :), &
         INTENT(INOUT)                       :: v1
    REAL(real_8), DIMENSION(:, :), &
         INTENT(INOUT)                       :: v2
    REAL(real_8), DIMENSION(:, :), &
         INTENT(INOUT)                       :: vt

    CHARACTER(len=*), PARAMETER              :: procedureN = 'cp_xc_compute'
    INTEGER                                  :: isub

    CALL tiset(procedureN,isub)
    !
    ! sgcx, sgcc, v1 and v2 must be zero upon entering cp_xc_compute
    !
    IF (cp_xc%use_libxc)  CALL cp_xc_compute_libxc()
    IF (cp_xc%use_cpfunc) CALL cp_xc_compute_internal()
    !
    CALL tihalt(procedureN,isub)

  CONTAINS

    !
    ! Internal driver using CPfunc
    !
    ! ==--------------------------------------------------------------==
    SUBROUTINE cp_xc_compute_internal()
      !
      ! Internal driver for exchange-correlation energy & potential in 
      ! R space
      !                             07.10.2016 M. P. Bircher @ LCBC/EPFL

    CHARACTER(*), PARAMETER :: procedureN = 'cp_xc_compute_internal'

    COMPLEX(real_8), &
      DIMENSION(cp_xc_spin_pairs)            :: v1s
    INTEGER                                  :: i_func, ierr, ir, isub
    LOGICAL                                  :: needs_tau, needs_gradients, &
                                                unrestricted_spins
    REAL(real_8), ALLOCATABLE, DIMENSION(:)  :: scs, sxs
    REAL(real_8), &
      DIMENSION(cp_xc_spin_components)       :: g_2, rho, v2s
    REAL(real_8), &
      DIMENSION(cp_xc_spin_pairs)            :: t, v1t
    TYPE(cp_xc_functional_p_t)               :: functional
    TYPE(cp_xc_functional_t), &
      DIMENSION(cp_xc_spin_components)       :: K_sigma, K_tmp
    TYPE(cp_xc_scratch_p_t)                  :: get_scratch
    TYPE(cp_xc_scratch_t), &
      DIMENSION(cp_xc_spin_components)       :: scratch

      CALL tiset(procedureN,isub)

      ALLOCATE(sxs( n ),scs( n ),stat=ierr)
      IF (ierr /= 0) CALL stopgm(procedureN,'Allocation problem',&
           __LINE__,__FILE__)

      ! Is this an initial guess in a spin-unrestricted run?
      !
      unrestricted_spins = (SIZE(v2,2) /= 1 .AND. cp_xc%polarized)
      !
      ! LDA, GGA or meta functional?
      needs_gradients = (cntl%tgc .OR. cntl%tgcx .OR. cntl%tgcc)
      needs_tau       = cntl%ttau

      SELECT CASE(unrestricted_spins)
         !
         ! LSDA
         !
      CASE(.TRUE.)
         ! Preliminaries
         !
         CALL cp_xc%set_scratch( get_scratch, needs_gradients, needs_tau )
         CALL cp_xc%set_functional( functional, dual_is_set=dual00%dual )

         ! Get xc energies and potentials
         !
         DO ir = 1, n
            !
            ! Perform checks here already instead of overwriting the values with zeros
            ! after computing them
            !
            rho(a:b) = rhoe( ir,a:b )
            rho(ab)  = rho(a) + rho(b)
            !
            IF (rho(ab)  < 1.0e-10_real_8 .OR. &
                 rho(ab) <= 0.1_real_8*cntr%gceps ) THEN
               sxs(ir) = 0.0_real_8
               scs(ir) = 0.0_real_8
               !
               ! v1 and v2 are properly zeroed already, simply leave them alone
               !
               CYCLE
            END IF
            !
            g_2(a)  = grad(ir,1) 
            g_2(b)  = grad(ir,5) 
            g_2(ab) = grad(ir,2)*grad(ir,6) + grad(ir,3)*grad(ir,7) + grad(ir,4)*grad(ir,8)
            t(a:b)  = tau( ir,a:b )
            !
            ! K_sigma collects the total enhancment factor for exchange plus energies and derivatives for
            ! exchange and correlation; scratch holds powers of rho.
            !
            CALL purge_K( K_sigma(:) )
            CALL get_scratch%unrestricted( scratch(:), rho(:), g_2(:), t(:) )
            !
            ! Get all functionals and build up K_sigma
            !
            DO i_func=1, cp_xc%cpfunc%n_funcs
               !
               CALL purge_K( K_tmp(:) )
               IF (ASSOCIATED(functional%Ks( i_func )%compute)) &
                    CALL functional%Ks( i_func )%compute( scratch(:), K_tmp(:) )
               CALL collect_K( K_tmp(:), cp_xc%cpfunc%scales( i_func ), K_sigma(:) )
               !
            END DO
            !
            ! Attenuate both exchange compenents if necessary (n.b. no mixed contribution)
            !
            IF (ASSOCIATED(functional%attenuate)) THEN
               IF (scratch(a)%n >= 1.0e-10_real_8) &
                    CALL functional%attenuate( scratch(a), K_sigma(a) )
               IF (scratch(b)%n >= 1.0e-10_real_8) &
                    CALL functional%attenuate( scratch(b), K_sigma(b) )
            END IF
            !
            ! Smoothing according to GC-cutoff and convert exchange derivatives to CPMD standard
            !
            CALL smoothing_unrestricted( K_sigma(a), K_sigma(b), K_sigma(ab), &
                 scratch(a), scratch(b), scratch(ab), &
                 sxs(ir), scs(ir), v1s(a), v1s(b), v2s(a), v2s(b), v2s(ab), &
                 v1t(a), v1t(b) )
            !
            v1( ir, a:b )  = v1( ir, a:b )  + v1s( a:b )
            v2( ir, a:ab ) = v2( ir, a:ab ) + v2s( a:ab )
            vt( ir, a:b )  = vt( ir, a:b )  + v1t( a:b )
            !
         END DO
         !
         sgcx = sgcx + SUM(sxs(:))
         sgcc = sgcc + SUM(scs(:))
         !
         ! Cleanup
         !
         CALL cp_xc%purge_scratch( get_scratch )
         CALL cp_xc%purge_functional( functional )

         !
         ! LDA
         !
      CASE(.FALSE.)
         ! Preliminaries
         !
         CALL cp_xc%set_scratch( get_scratch, needs_gradients, needs_tau )
         CALL cp_xc%set_functional( functional, dual_is_set=dual00%dual )

         ! Get xc energies and potentials
         !
         DO ir = 1, n
            !
            ! Perform checks here already instead of overwriting the values with zeros
            ! after computing them
            !
            IF (rhoe( ir, s )  < 1.0e-10_real_8 .OR. &
                 rhoe( ir, s ) <= 0.1_real_8*cntr%gceps ) THEN
               sxs(ir) = 0.0_real_8
               scs(ir) = 0.0_real_8
               !
               ! v1 and v2 are properly zeroed already, simply leave them alone
               !
               CYCLE
            END IF
            !
            ! K_sigma collects the total enhancment factor for exchange plus energies and derivatives for
            ! exchange and correlation; scratch holds powers of rho.
            !
            rho( s ) = rhoe( ir, s )
            g_2( s ) = grad( ir, s )
            t  ( s ) = tau ( ir, s )
            CALL purge_K( K_sigma(s) )
            CALL get_scratch%restricted( scratch(s), rho(s), g_2(s), t(s) )
            !
            ! Get all functionals and build up K_sigma
            !
            DO i_func=1, cp_xc%cpfunc%n_funcs
               !
               CALL purge_K( K_tmp(s) )
               IF (ASSOCIATED(functional%K( i_func )%compute)) CALL functional%K( i_func )%compute( scratch(s), K_tmp(s) )
               CALL collect_K( K_tmp(s), cp_xc%cpfunc%scales( i_func ), K_sigma(s) )
               !
            END DO
            !
            ! Attenuate all of the exchange compenents if necessary
            !
            IF (ASSOCIATED(functional%attenuate)) CALL functional%attenuate( scratch(s), K_sigma(s) )
            !
            ! Smoothing according to GC-cutoff and convert exchange derivatives to CPMD standard
            !
            CALL smoothing_restricted( K_sigma(s), scratch(s), sxs(ir), scs(ir), v1s(s), v2s(s), v1t(s) )
            !
            v1( ir, s ) = v1( ir,s ) + v1s( s )
            v2( ir, s ) = v2( ir,s ) + v2s( s )
            vt( ir, s ) = vt( ir,s ) + v1t( s )
            !
         END DO
         !
         sgcx = sgcx + SUM(sxs(:))
         sgcc = sgcc + SUM(scs(:))
         !
         ! Cleanup
         !
         CALL cp_xc%purge_scratch( get_scratch )
         CALL cp_xc%purge_functional( functional )
      END SELECT
      ! 
      DEALLOCATE(sxs,scs,stat=ierr)
      IF (ierr /= 0) CALL stopgm(procedureN,'Dellocation problem',&
           __LINE__,__FILE__)

      CALL tihalt(procedureN,isub)


    END SUBROUTINE cp_xc_compute_internal
    ! ==--------------------------------------------------------------==
    ELEMENTAL SUBROUTINE purge_K( K )


    TYPE(cp_xc_functional_t), INTENT(out)    :: K

      K%sx_sigma = 0.0_real_8
      K%K        = 0.0_real_8
      K%dK_dn    = 0.0_real_8
      K%dK_dg    = 0.0_real_8
      K%dsx_dn   = 0.0_real_8
      K%dsx_dg   = 0.0_real_8
      K%dsx_dt   = 0.0_real_8
      K%sc       = 0.0_real_8
      K%v1c      = 0.0_real_8
      K%v2c      = 0.0_real_8
      K%vtc      = 0.0_real_8

    END SUBROUTINE purge_K
    ! ==--------------------------------------------------------------==
    ELEMENTAL SUBROUTINE collect_K( K_in, weight, K_collect )
      !
      ! Assembles the overall enhancment factor K_sigma as a weighted sum
      ! of its components (needed for Coulomb attenuation)
      !                             07.10.2016 M. P. Bircher @ LCBC/EPFL

    TYPE(cp_xc_functional_t), INTENT(in)     :: K_in
    REAL(real_8), INTENT(in)                 :: weight
    TYPE(cp_xc_functional_t), INTENT(inout)  :: K_collect

      K_collect%sx_sigma = K_collect%sx_sigma + weight*K_in%sx_sigma
      K_collect%K        = K_collect%K        + weight*K_in%K
      K_collect%dK_dn    = K_collect%dK_dn    + weight*K_in%dK_dn
      K_collect%dK_dg    = K_collect%dK_dg    + weight*K_in%dK_dg
      K_collect%dsx_dn   = K_collect%dsx_dn   + weight*K_in%dsx_dn
      K_collect%dsx_dg   = K_collect%dsx_dg   + weight*K_in%dsx_dg
      K_collect%dsx_dt   = K_collect%dsx_dt   + weight*K_in%dsx_dt
      K_collect%sc       = K_collect%sc       + weight*K_in%sc
      K_collect%v1c      = K_collect%v1c      + weight*K_in%v1c
      K_collect%v2c      = K_collect%v2c      + weight*K_in%v2c
      K_collect%vtc      = K_collect%vtc      + weight*K_in%vtc

    END SUBROUTINE collect_K
    ! ==--------------------------------------------------------------==
    ELEMENTAL SUBROUTINE smoothing_unrestricted( K_a, K_b, K_ab, scratch_a, &
         scratch_b, scratch_ab, &
         sx,sc,v1a,v1b,v2a,v2b,v2ab, &
         vta,vtb )
      !
      ! Smoothing for low-gradient regions (this is applied to pure LDA,
      ! too).
      ! Version for open shell systems (conversion sx_sigma to sx)
      !                             19.10.2016 M. P. Bircher @ LCBC/EPFL


    TYPE(cp_xc_functional_t), INTENT(in)     :: K_a, K_b, K_ab
    TYPE(cp_xc_scratch_t), INTENT(in)        :: scratch_a, scratch_b, &
                                                scratch_ab
    REAL(real_8), INTENT(out)                :: sx, sc
    COMPLEX(real_8), INTENT(out)             :: v1a, v1b
    REAL(real_8), INTENT(out)                :: v2a, v2b, v2ab, vta, vtb

    REAL(real_8)                             :: dsmoo, sc_tmp, smoo, sx_tmp, &
                                                texp
    REAL(real_8), &
      DIMENSION(cp_xc_spin_components)       :: v2_tmp
    REAL(real_8), &
      DIMENSION(cp_xc_spin_pairs)            :: v1_tmp

!
! Conversion from internal storage to sx, v1 and v2
! 

      sx_tmp      = K_a%sx_sigma + K_b%sx_sigma
      sc_tmp      = K_a%sc       + K_b%sc      + K_ab%sc
      v1_tmp(a)   = K_a%dsx_dn   + K_a%v1c
      v1_tmp(b)   = K_b%dsx_dn   + K_b%v1c
      v2_tmp(a)   = K_a%dsx_dg/scratch_a%abs_g + K_a%v2c
      v2_tmp(b)   = K_b%dsx_dg/scratch_b%abs_g + K_b%v2c
      v2_tmp(ab)  = K_ab%v2c
      vta         = K_a%dsx_dt   + K_a%vtc
      vtb         = K_b%dsx_dt   + K_b%vtc
      !
      ! Smoothing
      !
      IF (scratch_ab%n > 4.0_real_8*cntr%gceps) THEN
         sx   = sx_tmp
         sc   = sc_tmp
         v1a  = v1_tmp(a)
         v1b  = v1_tmp(b)
         v2a  = v2_tmp(a)
         v2b  = v2_tmp(b)
         v2ab = v2_tmp(ab)
      ELSE
         texp  = EXP(3.0_real_8*(1.0_real_8 - scratch_ab%n/cntr%gceps))
         smoo  = 1.0_real_8/(1.0_real_8 + texp)
         dsmoo = 3.0_real_8/cntr%gceps * texp*smoo*smoo
         !
         sx = smoo*sx_tmp
         sc = smoo*sc_tmp
         !
         v1a  = dsmoo*(sx_tmp + sc_tmp) + smoo*v1_tmp(a)
         v1b  = dsmoo*(sx_tmp + sc_tmp) + smoo*v1_tmp(b)
         v2a  = smoo*v2_tmp(a)
         v2b  = smoo*v2_tmp(b)
         v2ab = smoo*v2_tmp(ab)
      END IF

    END SUBROUTINE smoothing_unrestricted
    ! ==--------------------------------------------------------------==
    ELEMENTAL SUBROUTINE smoothing_restricted( K_sigma,scratch,sx,sc,v1,v2,vt )
      !
      ! Smoothing for low-gradient regions (this is applied to pure LDA,
      ! too).
      ! Version for closed shell systems (conversion sx_sigma to sx)
      !                             07.10.2016 M. P. Bircher @ LCBC/EPFL


    TYPE(cp_xc_functional_t), INTENT(in)     :: K_sigma
    TYPE(cp_xc_scratch_t), INTENT(in)        :: scratch
    REAL(real_8), INTENT(out)                :: sx, sc
    COMPLEX(real_8), INTENT(out)             :: v1
    REAL(real_8), INTENT(out)                :: v2, vt

    REAL(real_8)                             :: abs_grad, dsmoo, rho, sc_tmp, &
                                                smoo, sx_tmp, texp, v1_tmp, &
                                                v2_tmp

!
! Conversion from internal storage to sx, v1 and v2
! 

      rho      = 2.0_real_8*scratch%n
      abs_grad = 2.0_real_8*scratch%abs_g
      sx_tmp   = 2.0_real_8*K_sigma%sx_sigma
      sc_tmp   = K_sigma%sc
      v1_tmp   = K_sigma%dsx_dn + K_sigma%v1c
      v2_tmp   = K_sigma%dsx_dg/abs_grad + K_sigma%v2c
      vt       = K_sigma%dsx_dt + K_sigma%vtc
      !
      ! Smoothing
      !
      IF (rho > 4.0_real_8*cntr%gceps) THEN
         sx = sx_tmp
         sc = sc_tmp
         v1 = v1_tmp
         v2 = v2_tmp
      ELSE
         texp  = EXP(3.0_real_8*(1.0_real_8 - rho/cntr%gceps))
         smoo  = 1.0_real_8/(1.0_real_8 + texp)
         dsmoo = 3.0_real_8/cntr%gceps * texp*smoo*smoo
         !
         sx = smoo*sx_tmp
         sc = smoo*sc_tmp
         !
         v1 = dsmoo*(sx_tmp + sc_tmp) + smoo*v1_tmp
         v2 = smoo*v2_tmp
      END IF

    END SUBROUTINE smoothing_restricted
    ! ==--------------------------------------------------------------==

    !
    ! LIBXC
    !
    ! ==--------------------------------------------------------------==
    SUBROUTINE cp_xc_compute_libxc()

    CHARACTER(*), PARAMETER :: procedureN = 'cp_xc_compute_libxc'

    INTEGER                                  :: i, i_func, isub
    LOGICAL                                  :: is_lsd_but_not
    REAL(real_8) :: dsmoo, fact, fact2, fact3, rho, rhoa, rhob, sc, scale, &
      smoo, sx, texp, v1_, v1a_, v1b_, v1c, v1ca, v1cb, v1x, v1xa, v1xb, v2_, &
      v2a_, v2ab_, v2b_, v2c, v2ca, v2cab, v2cb, v2x, v2xa, v2xab, v2xb
    REAL(real_8), ALLOCATABLE, DIMENSION(:)  :: e, e_tot
    REAL(real_8), ALLOCATABLE, &
      DIMENSION(:, :)                        :: vrho, vrho_tot, vsigma, &
                                                vsigma_tot
    REAL(real_8), DIMENSION(1)               :: e_
    REAL(real_8), DIMENSION(2)               :: rhoe_, vrho_
    REAL(real_8), DIMENSION(3)               :: sigma_, vsigma_

      CALL tiset(procedureN,isub)

      IF( cp_xc%polarized ) THEN

         !vw in the case we do lsd but not (atom guess), we hack !
         !vw need to fix that shit !
         is_lsd_but_not = SIZE(v2,2) == 1

         IF( is_lsd_but_not .AND. paral%io_parent ) THEN
            WRITE(6,*) 'WARNING: requested UKS call to LIBXC but has been overwitten by RKS.'
         ENDIF

         fact = 1.0_real_8
         fact2 = 1.0_real_8
         fact3 = 2.0_real_8
         IF( is_lsd_but_not ) THEN
            fact =0.5_real_8
            fact2 = 0.25_real_8
            fact3 = 1.0_real_8
         ENDIF

         ALLOCATE( vrho(n,2), vsigma(n,3), e_tot(n), vrho_tot(n,2), vsigma_tot(n,3) )
         e_tot(:)= 0.0_real_8; vrho_tot(:,:)=0.0_real_8;vsigma_tot(:,:)=0.0_real_8
         DO i_func = 1, cp_xc%libxc%n_funcs
            scale=cp_xc%libxc%scales(i_func)
            DO i = 1,n
               !vw we need to pass the rhoe as (2,n) and not (n,2), same for the other guys!
               rhoe_(:) = rhoe(i,:)
               SELECT CASE ( lxc_func_info_get_family(cp_xc%libxc%infos(i_func)) )
               CASE( XC_FAMILY_LDA )
                  CALL lxc_lda_exc_vxc (cp_xc%libxc%funcs(i_func), 1, fact*rhoe_, e_, vrho_ )
               CASE( XC_FAMILY_GGA, XC_FAMILY_HYB_GGA )
                  sigma_(1) = grad(i,1) ! aa
                  IF( is_lsd_but_not ) THEN
                     sigma_(3) = 0.0_real_8
                     sigma_(2) = 0.0_real_8
                     rhoe_(2) = 0.0_real_8
                  ELSE
                     sigma_(3) = grad(i,5) ! bb
                     sigma_(2) = (grad(i,2)*grad(i,6)+grad(i,3)*grad(i,7)+grad(i,4)*grad(i,8)) !ab
                  ENDIF
                  CALL lxc_gga_exc_vxc ( cp_xc%libxc%funcs(i_func), 1, fact*rhoe_, fact2*sigma_, e_, vrho_, vsigma_ )
                  vsigma_tot(i,1) = vsigma_tot(i,1) + scale * vsigma_(1) * fact3            !vw vsigma_ = alpha/alpha
                  vsigma_tot(i,2) = vsigma_tot(i,2) + scale * vsigma_(3) * fact3            !vw vsigma_ = beta/beta
                  vsigma_tot(i,3) = vsigma_tot(i,3) + scale * vsigma_(2) * fact3*0.5_real_8 !vw vsigma_ = alpha/beta
               CASE( XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA )
                  !xc_f90_mgga_exc_vxc
                  CALL stopgm(procedureN,'NYI',&
                       __LINE__,__FILE__)
               CASE default
                  CALL stopgm(procedureN,'not regonized functional',&
                       __LINE__,__FILE__)
               END SELECT
               e_tot(i) = e_tot(i) + scale * e_(1) * ( rhoe_(1) + rhoe_(2) )
               vrho_tot(i,:) = vrho_tot(i,:) + scale * vrho_(:)
            ENDDO
         ENDDO

         DO i = 1,n
            rhoa   = MAX(rhoe(i,1),0.0_real_8)
            IF( is_lsd_but_not ) THEN
               rhob = 0.0_real_8
            ELSE
               rhob   = MAX(rhoe(i,2),0.0_real_8)
            ENDIF
            rho = rhoa + rhob
            IF (rho.GT.0.1_real_8*cntr%gceps) THEN
               IF (rho.GT.4.0_real_8*cntr%gceps) THEN
                  smoo  = 1.0_real_8
                  dsmoo = 0.0_real_8
               ELSE
                  texp=EXP(3.0_real_8*(1.0_real_8-rho/cntr%gceps))
                  smoo=1.0_real_8/(1.0_real_8+texp)
                  dsmoo=3.0_real_8/cntr%gceps * texp*smoo*smoo
               ENDIF

               sx = e_tot(i)
               sc = 0.0_real_8

               v1xa = vrho_tot(i,1)
               v1ca = 0.0_real_8
               v1xb = vrho_tot(i,2)
               v1cb = 0.0_real_8
               v2xa = vsigma_tot(i,1)
               v2ca = 0.0_real_8
               v2xb = vsigma_tot(i,2)
               v2cb = 0.0_real_8
               v2xab = vsigma_tot(i,3)
               v2cab = 0.0_real_8

               sgcx  = sgcx + smoo*sx
               sgcc  = sgcc + smoo*sc

               v1a_     = dsmoo*(sx+sc) + smoo*(v1xa+v1ca)
               v2a_     = smoo*(v2xa+v2ca)
               v1b_     = dsmoo*(sx+sc) + smoo*(v1xb+v1cb)
               v2b_     = smoo*(v2xb+v2cb)
               v2ab_    = smoo*(v2xab+v2cab)

            ELSE
               v1a_     = 0.0_real_8
               v1b_     = 0.0_real_8
               v2a_     = 0.0_real_8
               v2b_     = 0.0_real_8
               v2ab_    = 0.0_real_8
            ENDIF

            v1(i,1) = v1a_
            IF(.NOT.is_lsd_but_not)v1(i,2) = v1b_
            v2(i,1) = v2a_
            IF(.NOT.is_lsd_but_not)v2(i,2) = v2b_
            IF(.NOT.is_lsd_but_not)v2(i,3) = v2ab_
         ENDDO

         DEALLOCATE( vrho, vsigma, e_tot, vrho_tot, vsigma_tot )


      ELSE


         ALLOCATE(e(n), vrho(n,1), vsigma(n,1), e_tot(n), vrho_tot(n,1), vsigma_tot(n,1) )

         e_tot(:)= 0.0_real_8; vrho_tot(:,:)=0.0_real_8;vsigma_tot(:,:)=0.0_real_8

         DO i_func = 1, cp_xc%libxc%n_funcs
            scale=cp_xc%libxc%scales(i_func)
            SELECT CASE ( lxc_func_info_get_family(cp_xc%libxc%infos(i_func)) )
            CASE( XC_FAMILY_LDA )
               CALL lxc_lda_exc_vxc (cp_xc%libxc%funcs(i_func), n, rhoe, e, vrho )
            CASE( XC_FAMILY_GGA, XC_FAMILY_HYB_GGA )
               CALL lxc_gga_exc_vxc ( cp_xc%libxc%funcs(i_func), n, rhoe, grad, e, vrho, vsigma )
               vsigma_tot(:,1) = vsigma_tot(:,1) + scale * vsigma(:,1) * 2.0_real_8
            CASE( XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA )
               !xc_f90_mgga_exc_vxc
               CALL stopgm(procedureN,'NYI',&
                    __LINE__,__FILE__)
            CASE default
               CALL stopgm(procedureN,'unrecognized functional',&
                    __LINE__,__FILE__)
            END SELECT
            e_tot(:) = e_tot(:) + scale * e(:) * rhoe(:,1)
            vrho_tot(:,1) = vrho_tot(:,1) + scale * vrho(:,1)

         ENDDO

         DO i = 1,n

            rho   = MAX(rhoe(i,1),0.0_real_8)
            IF (rho.GT.0.1_real_8*cntr%gceps) THEN
               IF (rho.GT.4.0_real_8*cntr%gceps) THEN
                  smoo  = 1.0_real_8
                  dsmoo = 0.0_real_8
               ELSE
                  texp=EXP(3.0_real_8*(1.0_real_8-rho/cntr%gceps))
                  smoo=1.0_real_8/(1.0_real_8+texp)
                  dsmoo=3.0_real_8/cntr%gceps * texp*smoo*smoo
               ENDIF

               sx = e_tot(i)
               sc = 0.0_real_8

               v1x = vrho_tot(i,1)
               v1c = 0.0_real_8
               v2x = vsigma_tot(i,1)
               v2c = 0.0_real_8

               sgcx  = sgcx + smoo*sx
               sgcc  = sgcc + smoo*sc
               v1_    = dsmoo*(sx+sc) + smoo*(v1x+v1c)
               v2_    = smoo*(v2x+v2c)
            ELSE
               v1_    = 0.0_real_8
               v2_    = 0.0_real_8
            ENDIF
            v1(i,1) = v1_
            v2(i,1) = v2_

         ENDDO

         DEALLOCATE( e, vrho, vsigma, e_tot, vrho_tot, vsigma_tot )


      ENDIF

      CALL tihalt(procedureN,isub)

    END SUBROUTINE cp_xc_compute_libxc
    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_xc_compute
  ! ==================================================================
END MODULE cp_xc_driver
