#define warning_level 5.e-3_real_8
! ==================================================================
#define mindistroutine calc_min_dist_alat_smooth
! define epr_mindistroutine epr_calc_min_dist_alat

MODULE epr_util_p_utils
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE nmr_position_p_utils,            ONLY: calc_min_dist_alat_smooth
  USE nmr_util_p_utils,                ONLY: make_123,&
                                             printtime
  USE parac,                           ONLY: paral
  USE perturbation_p_utils,            ONLY: ortho_p
  USE response_pmod,                   ONLY: response1,&
                                             timetag
  USE rwfopt_p_utils,                  ONLY: rwfopt_p
  USE soft,                            ONLY: soft_com
  USE system,                          ONLY: ncpw
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: epr_do_full

CONTAINS

  ! ==================================================================
  SUBROUTINE epr_do_full(i_state,iB,dowfopt,h1psi0,c0,c1,c1b,c1t,&
       psi,rhoe,eirop,eivps,z11,nstate)
    ! =----------------------------------------------------------------=
    ! Arguments:
    INTEGER                                  :: i_state, iB
    LOGICAL                                  :: dowfopt
    COMPLEX(real_8)                          :: psi(:)
    REAL(real_8)                             :: rhoe(*)
    COMPLEX(real_8)                          :: eirop(*), eivps(*)
    REAL(real_8)                             :: z11(*)
    INTEGER                                  :: nstate
    COMPLEX(real_8) :: c1t(ncpw%ngw,nstate,1), c1b(ncpw%ngw,nstate,3), &
      c1(ncpw%ngw,nstate,3), c0(ncpw%ngw,nstate), h1psi0(ncpw%ngw,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'epr_do_full'

    CHARACTER(len=1), DIMENSION(3)           :: xyz = (/'x','y','z'/)
    COMPLEX(real_8)                          :: cdummy(1)
    INTEGER                                  :: i_c1, ierr, ii, iii, is, isub
    INTEGER, DIMENSION(3), SAVE              :: i_state_previous = (/0,0,0/)
    REAL(real_8)                             :: dij(3), rdummy(1)
    REAL(real_8), ALLOCATABLE                :: aux(:)

89  FORMAT ("*  FULL calculation: state ",i4,", field in O",a1,20x,&
         "* ")
98  FORMAT (52("*"),"EPR*RESPONSE*")


    i_c1 = iB+6

    ! Initialization / very first call:
    IF (i_state_previous(iB) .EQ. 0) THEN
       CALL zeroing(c1t(:,:,1))!,ngw*nstate)
       response1%t_initialguess_c1 = .FALSE.
       i_state_previous(iB) = i_state
    ELSE
       response1%t_initialguess_c1 = .TRUE.
    ENDIF

    ! Preparation of B_iB . (d_k - d_fix) x p |psi_0_k>:
    IF (paral%io_parent) THEN
       WRITE (6,98)
       WRITE (6,89) i_state,xyz(iB)
       WRITE (6,98)
    ENDIF
    CALL make_123(iB,ii,iii)
    IF (dowfopt) THEN
       CALL zeroing(h1psi0)!,ngw*nstate)
       ALLOCATE(aux(2*ncpw%ngw),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       DO is=1,nstate        ! apply dd_ii p_iii - dd_iii p_ii:
          ! with dd = d_is - d_istate
          CALL mindistroutine(is,i_state,dij)
          CALL apply_op_p(c0(1,is), aux, ii, ncpw%ngw)
          CALL daxpy(2*ncpw%ngw, dij(iii), aux, 1, h1psi0(1,is), 1)
          CALL apply_op_p(c0(1,is), aux, iii, ncpw%ngw)
          CALL daxpy(2*ncpw%ngw,-dij(ii), aux, 1, h1psi0(1,is), 1)
       ENDDO
       DEALLOCATE(aux,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       timetag='preparation.'
       IF (paral%io_parent) CALL printtime
       IF (paral%io_parent) WRITE (6,98)
       IF (soft_com%exsoft) RETURN
    ENDIF

    ! Prepare a good guess for C1 from the previous optimization:
    ! this should be:
    ! B  .  {d_k - d_i} x |psi1^p>
    ! where k is the new and i is the old correction state.
    ! B = B * O_beta
    ! For an unknown reason, the convergence is better WITHOUT this
    ! adjustment. Bill Gates alone may have an idea why.
    IF (i_state .NE. i_state_previous(iB)) THEN
       CALL mindistroutine(i_state,i_state_previous(iB),dij)
       CALL daxpy(2*nstate*ncpw%ngw, dij(ii),c1(1,1,iii),1,c1t(1,1,1),1)
       CALL daxpy(2*nstate*ncpw%ngw,-dij(iii),c1(1,1,ii),1,c1t(1,1,1),1)
       CALL ortho_p(nstate,c0,c1t(1,1,1))
    ENDIF
    IF (dowfopt) THEN
       CALL tiset('PSI1-WF-OPT',isub)
       CALL rwfopt_p(c0,c1t(1,1,1),psi,rhoe,rdummy,eirop,&
            eivps,cdummy,h1psi0,z11,nstate,cdummy)
       CALL tihalt('PSI1-WF-OPT',isub)
       timetag='calculation of psi^1.'
       IF (paral%parent) CALL printtime
    ENDIF

    CALL daxpy(2*ncpw%ngw, 1._real_8, c1t(1,i_state,1), 1,&
         c1b(1,i_state,i_c1-6), 1)


    i_state_previous(iB) = i_state
    RETURN
  END SUBROUTINE epr_do_full
  ! ==================================================================
  ! ==================================================================

END MODULE epr_util_p_utils
