! ==================================================================
! Provides: - LDA exchange, adapted from OLDCODE
!
! Input:  Spin-polarised density (1/2 n for closed shell systems)
! Output: For spin-component only
! Convention: E_x = -0.5 \sum_s rho_s**(4/3)*K_sx
! Most functionals: E_x = 0.5 \sum_s (2rho_s)**(4/3) K_x
! Therefore   K_sx = 2**(4/3)K_x
!
!                             02.10.2016 - M. P. Bircher @ LCBC/EPFL
! ==================================================================
MODULE cp_lda_exchange_utils
  USE cpfunc_types,                    ONLY: cp_xc_functional_t,&
                                             cp_xc_scratch_t,&
                                             cp_xc_spin_components
  USE func,                            ONLY: func2
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER :: a                       = 1
  INTEGER, PARAMETER :: b                       = 2
  INTEGER, PARAMETER :: ab                      = 3

  REAL(real_8), PARAMETER :: slater_prefactor   = 1.10783814957303361_real_8
  REAL(real_8), PARAMETER :: four_thirds        = 4._real_8/3._real_8
  REAL(real_8), PARAMETER :: two_to_four_thirds = 2._real_8**(4._real_8/3._real_8)

  PUBLIC :: CP_LDA_X
  PUBLIC :: CP_SPIN_LDA_X

  PUBLIC :: CP_LDA_X_CHECK

CONTAINS

  ! ==================================================================
  PURE SUBROUTINE cp_spin_lda_x(scratch,functional)
    ! ==--------------------------------------------------------------==


    TYPE(cp_xc_scratch_t), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(in)                             :: scratch
    TYPE(cp_xc_functional_t), &
      DIMENSION(cp_xc_spin_components), &
      INTENT(inout)                          :: functional

!
! UEG
!

    functional(a)%K        = two_to_four_thirds*slater_prefactor*func2%salpha
    functional(b)%K        = functional(a)%K
    !
    functional(a)%dK_dn    = 0.0_real_8
    functional(a)%dK_dg    = 0.0_real_8
    functional(b)%dK_dn    = 0.0_real_8
    functional(b)%dK_dg    = 0.0_real_8
    !
    ! Proper conventions for readability
    !
    functional(a)%sx_sigma = -0.5_real_8*( scratch(a)%n_4_3*functional(a)%K )
    functional(a)%dsx_dn   = -0.5_real_8*( four_thirds*scratch(a)%n_1_3*functional(a)%K )
    functional(a)%dsx_dg   =  0.0_real_8
    functional(b)%sx_sigma = -0.5_real_8*( scratch(b)%n_4_3*functional(b)%K )
    functional(b)%dsx_dn   = -0.5_real_8*( four_thirds*scratch(b)%n_1_3*functional(b)%K )
    functional(b)%dsx_dg   =  0.0_real_8

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_spin_lda_x
  ! ==================================================================
  PURE SUBROUTINE cp_lda_x(scratch,functional)
    ! ==--------------------------------------------------------------==


    TYPE(cp_xc_scratch_t), INTENT(in)        :: scratch
    TYPE(cp_xc_functional_t), INTENT(inout)  :: functional

!
! UEG
!

    functional%K        = two_to_four_thirds*slater_prefactor*func2%salpha
    !
    functional%dK_dn    = 0.0_real_8
    functional%dK_dg    = 0.0_real_8
    !
    ! Proper conventions for readability
    !
    functional%sx_sigma = -0.5_real_8*( scratch%n_4_3*functional%K )
    functional%dsx_dn   = -0.5_real_8*( four_thirds*scratch%n_1_3*functional%K )
    functional%dsx_dg   = 0.0_real_8

    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_lda_x
  ! ==================================================================

  ! ==================================================================
  ! Check parameter compatibility
  ! ==================================================================
 
  ! ==================================================================
  ELEMENTAL FUNCTION cp_lda_x_check(tag) &
  RESULT (OK)
    ! ==--------------------------------------------------------------==
    ! Checks whether salpha has a reasonable value
    !
    !                                 25.07.2017 M.P.Bircher @ LCBC/EPFL
    ! ==--------------------------------------------------------------==
 
    CHARACTER(len=*), INTENT(in)             :: tag
    LOGICAL                                  :: OK 
 
    SELECT CASE( tag )
    CASE( "CP_LDA_X", "CP_SPIN_LDA_X" )
      OK = (func2%salpha /= 0.0_real_8)
    CASE DEFAULT
      OK = .false.
    END SELECT

    ! ==--------------------------------------------------------------==
  END FUNCTION cp_lda_x_check
  ! ==================================================================
END MODULE cp_lda_exchange_utils
