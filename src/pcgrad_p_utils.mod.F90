MODULE pcgrad_p_utils
  USE csize_utils,                     ONLY: csize
  USE error_handling,                  ONLY: stopgm
  USE fft_maxfft,                      ONLY: maxfftn
  USE fnonloc_utils,                   ONLY: fnonloc
  USE geq0mod,                         ONLY: geq0
  USE kinds,                           ONLY: int_1,&
                                             int_2,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8
  USE mp_interface,                    ONLY: mp_sum
  USE norm,                            ONLY: cnorm,&
                                             gemax
  USE parac,                           ONLY: parai,&
                                             paral
  USE perturbation_p_utils,            ONLY: energy,&
                                             ortho_p
  USE response_pmod,                   ONLY: dmbi,&
                                             response1,&
                                             response2,&
                                             s_star
  USE rnlsm_utils,                     ONLY: rnlsm
  USE simple_model_p_utils,            ONLY: simple_ortho_p
  USE system,                          ONLY: ncpw,&
                                             nkpt
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: zclean
  USE vpsi_utils,                      ONLY: vpsi
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: pcgrad_p
  PUBLIC :: precon_p

CONTAINS

  ! ==================================================================
  SUBROUTINE pcgrad_p(c0,c1,c12,vpp,psi,nstate,reset_cg,z11)
    ! ==--------------------------------------------------------------==
    ! ==  preconditioned conjugate gradient optimization              ==
    ! ==                                                              ==
    ! ==  INPUT: C12:          the force  d E / d Psi_1               ==
    ! ==         C0,vpp,z11:   the usual ingredients                  ==
    ! ==                                                              ==
    ! ==  NB: The gradient does NOT need to be orthogonalized wrt the ==
    ! ==      ground state.                                           ==
    ! ==      This routine also checks the CONVERGENCE CONDITION and  ==
    ! ==      computes the gemax and cnorm values.                    ==
    ! ==                                                              ==
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8)                          :: c0(*), c1(*), c12(*)
    REAL(real_8)                             :: vpp(*)
    COMPLEX(real_8)                          :: psi(maxfftn)
    INTEGER                                  :: nstate
    LOGICAL                                  :: reset_cg
    REAL(real_8)                             :: z11(*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'pcgrad_p'

    COMPLEX(real_8), ALLOCATABLE             :: prec12(:)
    COMPLEX(real_8), ALLOCATABLE, SAVE       :: A__dir(:), dir(:), g(:)
    INTEGER                                  :: ierr, isub
    INTEGER, SAVE                            :: ifirst = 0, istep
    REAL(real_8)                             :: alpha, g_PREC_g__old
    REAL(real_8), SAVE :: avg_alpha = HUGE(0._real_8), &
      beta = HUGE(0._real_8), dir__A__dir = HUGE(0._real_8), &
      dir_grad = HUGE(0._real_8), g_PREC_g = HUGE(0._real_8), &
      g_PREC_g_old = HUGE(0._real_8)

! variables & functions
! conj grad variables:
! dir or d means direction (numerical recipes: the "h" or "p")
! g means gradient  (numerical recipes: the "g" or "r")
! d_old is the (conjugate) direction of the _last_ step.
! A is "the matrix", i.e. (H0 - eps_ij)
!COMMON     /static_vars/ avg_alpha,beta,g_PREC_g,g_PREC_g_old,&
!     dir_grad,dir__A__dir
! Attention, __old means that the _number_ is from the last
! cycle, whereas _old refers to the _vector_ from last cycle.
! technical variables:
! ==--  Initialization  ------------------------------------------==

    CALL tiset('  pcgrad_p',isub)
    IF (reset_cg) THEN
       istep = 1
       reset_cg = .FALSE.
       avg_alpha = 0._real_8
    ENDIF
    IF (ifirst.NE.2604) THEN
       ifirst = 2604         ! TODO WTF???
       istep=1
       ALLOCATE(dir(ncpw%ngw*nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(   dir)!,SIZE(dir))
       ALLOCATE(A__dir(ncpw%ngw*nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(   A__dir)!,SIZE(A__dir))
       IF (response1%tpolak_p) THEN
          ALLOCATE(g(ncpw%ngw*nstate),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          CALL zeroing(   g)!,SIZE(g))
       ENDIF
    ENDIF
    ! ==-  Preconditioning ------------------------------------------==
    ALLOCATE(prec12(nkpt%ngwk*nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! DMB
    IF (dmbi%tsimple_model) THEN
       CALL simple_ortho_p(nstate,c0,c12,S_star)
    ELSE
       CALL ortho_p(nstate,c0,c12)
    ENDIF
    ! DMB
    IF (response1%prec_p) THEN
       CALL precon_p(c12,prec12,vpp,nstate,z11,1)
       ! DMB
       IF (dmbi%tsimple_model) THEN
          CALL simple_ortho_p(nstate,c0,prec12,S_star)
       ELSE
          CALL ortho_p(nstate,c0,prec12)
       ENDIF
       ! DMB
    ELSE
       CALL dcopy(2*ncpw%ngw*nstate,c12,1,prec12,1)
    ENDIF

    CALL csize(c12,nstate,gemax,cnorm)

    IF (geq0) THEN
       CALL zclean(prec12,nstate,ncpw%ngw)
       CALL zclean(c12,nstate,ncpw%ngw)
    ENDIF
    ! C12 is the original gradient, PREC12 the preconditioned one.


    ! Some scalar products, needed for the CG factors alpha and beta:
    g_PREC_g_old  = 0._real_8
    g_PREC_g__old = g_PREC_g
    g_PREC_g      = energy(prec12,c12,nstate)
    IF (response1%tpolak_p) THEN
       IF (istep .GT. 1) THEN
          g_PREC_g_old = energy(g,c12,nstate)
       ELSE
          g_PREC_g_old = 0._real_8
       ENDIF
       CALL dcopy(2*ncpw%ngw*nstate,prec12,1,g,1)
    ENDIF
    CALL mp_sum(g_PREC_g,parai%allgrp)   ! g_PREC_g and g_PREC_g_old
    CALL mp_sum(g_PREC_g_old,parai%allgrp)   ! g_PREC_g and g_PREC_g_old




    ! calculation of the new conj direction:
    ! 
    ! dir_new = precon * gradient  +  beta * dir_old
    ! 
    ! with beta = - (gradient*precon*A*dir_old) / (dir_old*A*dir_old)
    ! 
    ! or alternatively,
    ! beta = (gradient*precon*gradient)_present / ()_last_cycle
    ! 
    ! or with the POLAK formula
    ! beta = ( (gradient_present - gradiend_last_cycle) 
    ! \                      * precon* gradient_present)
    ! \        / (gradient*precon*gradient)_last_cycle
    ! 
    IF (istep .EQ. 1) THEN
       CALL zeroing(dir)!,SIZE(dir))
       beta = 0._real_8
       CALL dcopy(2*ncpw%ngw*nstate, prec12,1, dir,1)
       avg_alpha = 1.0_real_8 * response2%cg_factor
    ELSE                      ! then istep .ge. 2
       IF (response1%tpolak_p) THEN
          IF (g_PREC_g__old .GT. 1.e-20_real_8)&
               beta = (g_PREC_g - g_PREC_g_old) / g_PREC_g__old
       ELSE
          IF (g_PREC_g__old .GT. 1.e-20_real_8) beta = g_PREC_g /&
               g_PREC_g__old
       ENDIF
       IF (beta .GT. 0.5_real_8)  beta = 0.0_real_8
       CALL dscal(2*nstate*ncpw%ngw,beta,dir,1)
       CALL daxpy(2*ncpw%ngw*nstate, 1._real_8, prec12,1, dir,1)
       ! THIS (dir) is now the NEW CONJUGATE DIRECTION
    ENDIF

    DEALLOCATE(prec12,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)

    IF (istep.GT.1 .AND. istep.LE.response1%cg_analytic) THEN
       ! Calculation of the step length (alpha):
       ! alpha = - (dir * gradient) / (dir * A * dir)
       CALL a_p(dir,A__dir,c0,nstate,psi,z11)
       dir__A__dir = energy(dir,A__dir,nstate)
       dir_grad    = energy(dir,c12,     nstate)
       CALL mp_sum(dir_grad,parai%allgrp)! dir_grad and dir__A__dir
       CALL mp_sum(dir__A__dir,parai%allgrp)! dir_grad and dir__A__dir
       alpha       = - dir_grad / dir__A__dir
       alpha = alpha * response2%cg_factor
       ! alpha is now the COMPUTED STEP LENGTH.
       avg_alpha = (avg_alpha * (istep-1) + alpha)/(istep*1._real_8)
    ELSE
       IF (ABS(avg_alpha) .LT. 0.01_real_8) avg_alpha = 0.01_real_8
       alpha = avg_alpha
    ENDIF

    ! APPLY THE CONJUGATE GRADIENT STEP:
    CALL daxpy(2*ncpw%ngw*nstate, alpha, dir, 1, c1, 1)
    CALL ortho_p(nstate,c0,c1)

    istep = istep + 1
    CALL tihalt('  pcgrad_p',isub)
    RETURN
  END SUBROUTINE pcgrad_p
  ! ==================================================================


  ! ==================================================================
  SUBROUTINE precon_p(c2,c2dest,vpp,nstate,z11,direction)
    ! APPLIES  (H0 - E_ks) ^ -1  to the argument c2.
    ! vpp must contain    G^2 + Coulomb-terms(G)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c2(ncpw%ngw,*), &
                                                c2dest(ncpw%ngw,*)
    REAL(real_8)                             :: vpp(ncpw%ngw)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: z11(nstate,nstate)
    INTEGER                                  :: direction

    CHARACTER(*), PARAMETER                  :: procedureN = 'precon_p'

    COMPLEX(real_8)                          :: c2value
    INTEGER                                  :: ierr, ig, is
    INTEGER, SAVE                            :: vpp2_allocated = 0
    REAL(real_8)                             :: arg, hthrs2, vpp2value, &
                                                z11_avg
    REAL(real_8), ALLOCATABLE, SAVE          :: vpp2(:,:)

! ==--------------------------------------------------------------==
! New technique: The new array vpp2 will now contain the ready-to-use 
! preconditioner. It is computed at the first time this routine is called.
! For the simple preconditioners (one for all states), it has dimension
! NGW, for the state-dependent preconditioners NGW*NSTATE.
! The routine can now also apply the INVERSE preconditioner. To do
! so, -1 has to be given as argument (direction). A normal call would
! use +1.
! Simple preconditioners (state-invariant)

    IF (direction .EQ. 1) THEN
       IF (response1%preconditioner_p .EQ. 1.OR. response1%preconditioner_p .EQ. 2)&
            THEN
          IF (vpp2_allocated .NE. 2604) THEN! the it is not yet allocated.
             IF (vpp2_allocated .NE. 0)CALL stopgm('pcgrad/precon',&
                  'INTERNAL ERROR: PRECONDITIONER TYPE CHANGED',& 
                  __LINE__,__FILE__)
             vpp2_allocated = 2604
             ALLOCATE(vpp2(ncpw%ngw,1),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             CALL zeroing(vpp2)!,ngw)

             z11_avg = 0._real_8
             DO is=1,nstate
                z11_avg = z11_avg - z11(is,is)
             ENDDO
             z11_avg = z11_avg / REAL(nstate,kind=real_8)
             hthrs2  = response2%hthrs_p*response2%hthrs_p

             IF (response1%preconditioner_p .EQ. 1) THEN
                !$omp parallel do private(ig,arg)
                DO ig=1,ncpw%ngw
                   arg = vpp(ig) + z11_avg
                   arg = 1._real_8 / SQRT( arg*arg + hthrs2 )
                   IF (ABS(arg) .LE. 1.e-10_real_8) arg = 1.e-10_real_8
                   ! This is for safety: at a given point, we need the inverse
                   ! of vpp2!
                   vpp2(ig,1) = arg
                ENDDO
             ELSE IF (response1%preconditioner_p .EQ. 2) THEN
                !$omp parallel do private(ig,arg)
                DO ig=1,ncpw%ngw
                   arg = vpp(ig) + z11_avg
                   arg = arg / ( arg*arg + hthrs2 )
                   IF (ABS(arg) .LE. 1.e-10_real_8) arg = 1.e-10_real_8
                   vpp2(ig,1) = arg
                ENDDO
             ELSE
                CALL stopgm('pcgrad_p/precon',&
                     'INTERNAL ERROR: preconditioner index',& 
                     __LINE__,__FILE__)
             ENDIF
          ENDIF             ! allocation+calculation of preconditioner

          ! now, use the preconditioner:
          !$omp parallel do private(ig,is,vpp2value,c2value)
          DO ig=1,ncpw%ngw
             vpp2value     = vpp2(ig,1)
             DO is=1,nstate
                c2value       = c2(ig,is)
                c2dest(ig,is) = c2value * vpp2value
             ENDDO
          ENDDO


          ! STATE-DEPENDENT preconditioners:
       ELSEIF      (response1%preconditioner_p .EQ. 3.OR. response1%preconditioner_p&
            .EQ.4) THEN
          IF (vpp2_allocated .NE. 2602) THEN! the it is not yet allocated.
             IF (vpp2_allocated .NE. 0)CALL stopgm('pcgrad/precon',&
                  'INTERNAL ERROR: PRECONDITIONER TYPE CHANGED',& 
                  __LINE__,__FILE__)
             vpp2_allocated = 2602
             ALLOCATE(vpp2(ncpw%ngw,ncpw%ngw/ncpw%ngw),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             CALL zeroing(vpp2)!,ngw*nstate)

             hthrs2  = response2%hthrs_p*response2%hthrs_p
             IF (response1%preconditioner_p .EQ. 3) THEN
                !$omp parallel do private(is,ig,arg) 
                DO is=1,nstate
                   DO ig=1,ncpw%ngw
                      arg = vpp(ig)-z11(is,is)
                      arg = arg * arg + hthrs2
                      arg = 1._real_8 / SQRT(arg)
                      IF (ABS(arg) .LE. 1.e-10_real_8) arg = 1.e-10_real_8
                      vpp2(ig,is) = arg
                   ENDDO
                ENDDO
             ELSEIF (response1%preconditioner_p .EQ. 4) THEN
                !$omp parallel do private(is,ig,arg) 
                DO is=1,nstate
                   DO ig=1,ncpw%ngw
                      arg = vpp(ig)-z11(is,is)
                      arg = arg / (arg * arg + hthrs2)
                      IF (ABS(arg) .LE. 1.e-10_real_8) arg = 1.e-10_real_8
                      vpp2(ig,is) = arg
                   ENDDO
                ENDDO
             ELSE
                CALL stopgm('pcgrad_p/precon',&
                     'INTERNAL ERROR: preconditioner index',& 
                     __LINE__,__FILE__)
             ENDIF
          ENDIF             ! allocation+calculation of preconditioner


          ! now, use the preconditioner:
          !$omp parallel do private(is,ig,c2value) 
          DO is=1,nstate
             DO ig=1,ncpw%ngw
                c2value       = c2(ig,is)
                c2dest(ig,is) = c2value * vpp2(ig,is)
             ENDDO
          ENDDO
       ELSE
          CALL stopgm('PCGRAD_P','UNKNOWN PRECONDITIONER TYPE',& 
               __LINE__,__FILE__)
       ENDIF


    ELSEIF (direction .EQ. -1) THEN
       IF (vpp2_allocated .EQ. 2604&
            .OR. vpp2_allocated .EQ. 2602) THEN
          IF (response1%preconditioner_p .EQ. 1&
               .OR. response1%preconditioner_p .EQ. 2) THEN
             !$omp parallel do private(is,ig,c2value)
             DO is=1,nstate
                DO ig=1,ncpw%ngw
                   c2value       = c2(ig,is)
                   c2dest(ig,is) = c2value / vpp2(ig,1)
                ENDDO
             ENDDO
          ELSEIF (response1%preconditioner_p .EQ. 3&
               .OR. response1%preconditioner_p .EQ. 4) THEN
             !$omp parallel do private(is,ig,c2value)
             DO is=1,nstate
                DO ig=1,ncpw%ngw
                   c2value       = c2(ig,is)
                   c2dest(ig,is) = c2value / vpp2(ig,is)
                ENDDO
             ENDDO
          ELSE
             CALL stopgm('PCGRAD_P','UNKNOWN PRECONDITIONER TYPE',& 
                  __LINE__,__FILE__)
          ENDIF
       ELSE
          CALL stopgm('pcgrad_p',&
               'INTERNAL ERROR:inverse preconditioner '//&
               'called too early',& 
               __LINE__,__FILE__)
       ENDIF
    ELSE
       CALL stopgm('pcgrad_p','INTERNAL PROBLEM: old format',& 
            __LINE__,__FILE__)
    ENDIF


    RETURN
  END SUBROUTINE precon_p
  ! ==================================================================

  ! ==================================================================

END MODULE pcgrad_p_utils


! ==================================================================
SUBROUTINE a_p(arg,RESULT,c0,nstate,psi,z11)
  ! ==--------------------------------------------------------------==
  ! The energy functional is of the form
  ! E(c1) = <c1| F0 |c1>  +  E_nonloc[c1]. Here we calculate the
  ! the bilinear part of the functional ( F0 |c1> ):
  ! (1-P_o)  ( h0  |arg_a>  -  eps_lagr_ba  |arg_b> )
  ! where arg (=argument) is the vector (of wavefunctions)
  ! on which the operator is applied. The resulting vector
  ! is then stored in result.
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE fft_maxfft,                      ONLY: maxfftn
  USE timer, ONLY: tiset, tihalt
  USE system , ONLY: ncpw
  USE parac, ONLY : paral,parai
  USE elct , ONLY: crge
  USE response_pmod , ONLY:nmr_options,response1,vofrho0,voa_data
  USE vpsi_utils, ONLY : vpsi
  USE fnonloc_utils, ONLY : fnonloc
  USE rnlsm_utils, ONLY : rnlsm
  USE spin, ONLY : clsd
  USE zeroing_utils,                   ONLY: zeroing
  IMPLICIT NONE
  INTEGER                                    :: nstate
  COMPLEX(real_8)                            :: c0(ncpw%ngw,nstate), &
                                                RESULT(ncpw%ngw,nstate), &
                                                arg(ncpw%ngw,nstate), &
                                                psi(maxfftn)
  REAL(real_8)                               :: z11(nstate,nstate)

  INTEGER                                    :: i

! ==--------------------------------------------------------------==
! H0  |arg_a_1>:

  CALL zeroing(RESULT)
  CALL rnlsm(arg,nstate,1,1,.FALSE.)
  CALL vpsi(arg,RESULT,crge%f(:,1),vofrho0,psi,nstate,1,1,.FALSE.)
  CALL fnonloc(RESULT,crge%f,nstate,1,clsd%nlsd,.FALSE.)

  ! - epsilon |arg>:
  IF (((response1%tnmr.OR.response1%tepr).AND.nmr_options%tlocalize).OR. &
       response1%tinteraction.OR.(response1%tvoa.AND.voa_data%timagpert)) THEN
     CALL dsymm ( 'R','U', 2*ncpw%ngw, nstate, 1._real_8,z11, nstate, arg,&
          2*ncpw%ngw, 1._real_8, RESULT, 2*ncpw%ngw)
  ELSE
     DO i=1,nstate
        CALL daxpy(2*ncpw%ngw, z11(i,i), arg(1,i),1, RESULT(1,i),1)
     ENDDO
  ENDIF

  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE a_p
