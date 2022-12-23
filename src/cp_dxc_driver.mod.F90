MODULE cp_dxc_driver
  USE cp_xc_utils,                     ONLY: cp_xc,&
                                             cp_dxc_functional_p_t,&
                                             cp_dxc_scratch_p_t
  USE cpfunc_types,                    ONLY: cp_dxc_functional_t,&
                                             cp_dxc_scratch_t,&
                                             cp_xc_spin_components,&
                                             cp_xc_spin_pairs,&
                                             cp_xc_abs_xyz,&
                                             cp_xc_xyz
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: paral
  USE system,                          ONLY: cntl,&
                                             cntr

  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER      :: a        = 1
  INTEGER, PARAMETER      :: b        = 2
  INTEGER, PARAMETER      :: ab       = 3
  INTEGER, PARAMETER      :: s        = a

  PUBLIC :: cp_dxc_compute

CONTAINS

  ! ==================================================================
  SUBROUTINE cp_dxc_compute( n, rhoe, drhoe, grad, dgrad, dxc_tmp )
    ! ==--------------------------------------------------------------==
    ! Analytical derivatives for linear response with CPfunc
    !                               19.06.2017 M. P. Bircher @ LCBC/EPFL
    ! ==--------------------------------------------------------------==

    INTEGER, INTENT(IN)                      :: n
    REAL(real_8), DIMENSION(:,:), &
         INTENT(IN)                          :: rhoe, drhoe, grad, dgrad
    REAL(real_8), DIMENSION(:,:,:), &
         INTENT(INOUT)                       :: dxc_tmp

    CHARACTER(len=*), PARAMETER              :: procedureN = 'cp_dxc_compute'

    LOGICAL, PARAMETER                       :: needs_gradients = .TRUE.

    INTEGER                                  :: i_func, ir
    LOGICAL                                  :: needs_tau, unrestricted_spins

    TYPE(cp_dxc_functional_p_t)              :: functional
    TYPE(cp_dxc_functional_t), &
      DIMENSION(cp_xc_spin_pairs)            :: dK_sigma, dK_tmp
    TYPE(cp_dxc_scratch_p_t)                 :: get_scratch_derivatives
    TYPE(cp_dxc_scratch_t), &
      DIMENSION(cp_xc_spin_components)       :: scratch


    unrestricted_spins = cp_xc%polarized
    !
    ! LDA, GGA or meta functional?
    needs_tau = cntl%ttau

    SELECT CASE(unrestricted_spins)
       !
       ! LSDA
       !
    CASE(.TRUE.)
       ! Preliminaries
       !
       CALL cp_xc%set_scratch_derivatives( get_scratch_derivatives, needs_gradients, needs_tau )
       CALL cp_xc%set_functional_derivatives( functional )

       ! Get xc linres derivatives
       !
       DO ir = 1, n
          !
          IF ( (rhoe(ir,a) + rhoe(ir,b)) <= 0.0_real_8) THEN
             !
             ! dxc_tmp is properly zeroed already, simply leave it alone
             !
             CYCLE
          END IF
          !
          ! dK_sigma collects the total enhancment factor for exchange plus energies and derivatives for
          ! exchange and correlation; scratch holds powers of rho.
          !
          CALL purge_dK( dK_sigma(:) )
          CALL get_scratch_derivatives%unrestricted( scratch(:), rhoe(ir,:), drhoe(ir,:), &
                                                     grad(ir,:), dgrad(ir,:) )
          !
          ! Get all functionals and build up K_sigma
          !
          DO i_func=1, cp_xc%cpfunc%n_funcs
             !
             CALL purge_dK( dK_tmp(:) )
             IF (ASSOCIATED(functional%dKs( i_func )%compute)) &
                  CALL functional%dKs( i_func )%compute( scratch(:), dK_tmp(:) )
             CALL collect_dK( dK_tmp(:), cp_xc%cpfunc%scales( i_func ), dK_sigma(:) )
             !
          END DO
          !
          CALL copy_dxc( dK_sigma(a), dxc_tmp(ir,1,a), dxc_tmp(ir,2,a), dxc_tmp(ir,3,a), dxc_tmp(ir,4,a), dxc_tmp(ir,5,a) )
          CALL copy_dxc( dK_sigma(b), dxc_tmp(ir,1,b), dxc_tmp(ir,2,b), dxc_tmp(ir,3,b), dxc_tmp(ir,4,b), dxc_tmp(ir,5,b) )
          !
       END DO
       !
       ! Cleanup
       !
       CALL cp_xc%purge_scratch_derivatives( get_scratch_derivatives )
       CALL cp_xc%purge_functional_derivatives( functional )

       !
       ! LDA
       !
    CASE(.FALSE.)
       ! Preliminaries
       !
       CALL cp_xc%set_scratch_derivatives( get_scratch_derivatives, needs_gradients, needs_tau )
       CALL cp_xc%set_functional_derivatives( functional )

       ! Get xc energies and potentials
       !
       DO ir = 1, n
          !
          ! Perform checks here already instead of overwriting the values with zeros
          ! after computing them
          !
          IF (rhoe( ir, s )  <= 0.0_real_8) THEN
             !
             ! dxc_tmp is properly zeroed already, simply leave it alone
             !
             CYCLE
          END IF
          !
          ! K_sigma collects the total enhancment factor derivatives
          !
          CALL purge_dK( dK_sigma(s) )
          CALL get_scratch_derivatives%restricted( scratch(s), rhoe(ir,s), drhoe(ir,s), &
                                                   grad(ir,1:cp_xc_abs_xyz), dgrad(ir,1:cp_xc_xyz) )
          !
          ! Get all functionals and build up K_sigma
          !
          DO i_func=1, cp_xc%cpfunc%n_funcs
             !
             CALL purge_dK( dK_tmp(s) )
             IF (ASSOCIATED(functional%dK( i_func )%compute)) CALL functional%dK( i_func )%compute( scratch(s), dK_tmp(s) )
             CALL collect_dK( dK_tmp(s), cp_xc%cpfunc%scales( i_func ), dK_sigma(s) )
             !
          END DO
          !
          CALL copy_dxc( dK_sigma(s), dxc_tmp(ir,1,s), dxc_tmp(ir,2,s), dxc_tmp(ir,3,s), dxc_tmp(ir,4,s), dxc_tmp(ir,5,s) )
          !
       END DO
       !
       ! Cleanup
       !
       CALL cp_xc%purge_scratch_derivatives( get_scratch_derivatives )
       CALL cp_xc%purge_functional_derivatives( functional )
    END SELECT

 
    CONTAINS 

    ! ==--------------------------------------------------------------==
    ELEMENTAL SUBROUTINE purge_dK( dK )

    TYPE(cp_dxc_functional_t), INTENT(out)   :: dK

      dK%tmp1 = 0.0_real_8
      dK%tmp2 = 0.0_real_8
      dK%tmp3 = 0.0_real_8
      dK%tmp4 = 0.0_real_8
      dK%tmp5 = 0.0_real_8

    END SUBROUTINE purge_dK
    ! ==--------------------------------------------------------------==
    ELEMENTAL SUBROUTINE collect_dK( dK_in, weight, dK_collect )
      !
      ! Assembles the overall derivatives in dK_sigma as a weighted sum
      ! of its components
      !                             19.06.2017 M. P. Bircher @ LCBC/EPFL

    TYPE(cp_dxc_functional_t), INTENT(in)    :: dK_in
    REAL(real_8), INTENT(in)                 :: weight
    TYPE(cp_dxc_functional_t), INTENT(inout) :: dK_collect

      dK_collect%tmp1 = dK_collect%tmp1 + weight*dK_in%tmp1
      dK_collect%tmp2 = dK_collect%tmp2 + weight*dK_in%tmp2
      dK_collect%tmp3 = dK_collect%tmp3 + weight*dK_in%tmp3
      dK_collect%tmp4 = dK_collect%tmp4 + weight*dK_in%tmp4
      dK_collect%tmp5 = dK_collect%tmp5 + weight*dK_in%tmp5

    END SUBROUTINE collect_dK
    ! ==--------------------------------------------------------------==
    ELEMENTAL SUBROUTINE copy_dxc( dK_sigma, tmp1, tmp2, tmp3, tmp4, tmp5 )
      !
      ! Copies from work type to temporary scratch space used in dd_xc_ana
      ! (This ensures compatibility with the latter)
      !                             19.06.2017 M. P. Bircher @ LCBC/EPFL

    TYPE(cp_dxc_functional_t), INTENT(in)    :: dK_sigma

    REAL(real_8), INTENT(out)                :: tmp1, tmp2, tmp3, tmp4, tmp5

      tmp1 = dK_sigma%tmp1 
      tmp2 = dK_sigma%tmp2 
      tmp3 = dK_sigma%tmp3 
      tmp4 = dK_sigma%tmp4 
      tmp5 = dK_sigma%tmp5 

    END SUBROUTINE copy_dxc
    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_dxc_compute
  ! ==================================================================
END MODULE cp_dxc_driver
