MODULE rlbfgs_utils
  USE adapttol_utils,                  ONLY: tol_init_shared
  USE coor,                            ONLY: tau0
  USE cotr,                            ONLY: cotc0,&
                                             hess
  USE error_handling,                  ONLY: stopgm
  USE hessin_utils,                    ONLY: clcmap,&
                                             hessin,&
                                             phessci
  USE ions,                            ONLY: ions1
  USE kinds,                           ONLY: int_1,&
                                             int_2,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8
  USE lscal,                           ONLY: &
       beta_l, betap_l, cuteig, deaccp, demin_p, depred, eigval, eigvec, &
       elestk, eps_h, etotmn, fstep_p, gx, hesscr, ibnd_l, icore, ielstk, &
       igpnt, info_l, iopnxt, ipnt_l, iprnt, ispnt, istlen, iszstk, iter_l, &
       itpstk, iwpnt, jump, lallhs, lfevmn, linils, llock, lmap_p, lmicro, &
       lnomap, lprhes, lrstrt, lundop, lvlhes, lwlbfgs, lwprfo, m_hess, &
       map_p, mode, modelk, ncore, nfev_l, nfevmn, npt_l, nrem, nres_l, &
       nrestt, nshess, nsmaxp, nstack, nstep_l, nstep_p, nstried_l, &
       nstried_p, nsvib, ntrstr, ntrust, nvar, oldeig, oldene, oldev, oldg, &
       oldgx, oldhss, omin, ooldg, oprven, ostep_p, osz_p, prvene, rmax_p, &
       rmin_p, scl1_l, scl2_l, st_get, st_nop, st_pop, st_push, st_put, step, &
       step_bmb_l, step_ini_l, step_ini_p, step_l, step_max_l, step_max_p, &
       step_min_l, step_min_p, step_p, sz_l, sz_lmn, sz_p, tolenv, trustp, &
       trustr, vmode, wlbfgs
  USE parac,                           ONLY: paral
  USE sdion_utils,                     ONLY: minlin,&
                                             scalar
  USE secder_utils,                    ONLY: purgeh
  USE store_types,                     ONLY: cprint,&
                                             iprint_lscal
  USE system,                          ONLY: cnti,&
                                             cntr
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE vibana_utils,                    ONLY: vibana
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rlbfgs
  PUBLIC :: stack_init
  !public :: stack_sched
!!!public :: stack_do_sched
  PUBLIC :: stack_do
  PUBLIC :: stack_destroy
  !public :: stack_def
  PUBLIC :: lscal_init
  !public :: lbfgs_init
  !public :: lsrch
  !public :: mcsrch


  PUBLIC :: rprfo
  !public :: inihes
  !public :: updhes
  !public :: prjhes
  !public :: ovrlap
  !public :: frmstp
  !public :: crdco
  !public :: crdci
  !public :: cnvenv
!!$public :: clcmap
  PUBLIC :: prfo_init
  PUBLIC :: pvibana
!!$public :: phessci
  !public :: phessco

  PUBLIC :: give_work_lbfgs
  PUBLIC :: give_scr_rlbfgs

  !public :: give_work_prfo
  PUBLIC :: give_scr_rprfo


CONTAINS

  ! ==================================================================
  SUBROUTINE rlbfgs (xpar, dxpar, ndim, etot, lret)
    ! ==--------------------------------------------------------------==
    ! ==  LIMITED MEMORY cntl%bfgs QUASI-NEWTON OPTIMIZER                  ==
    ! ==                                                              ==
    ! ==  Reference:                                                  ==
    ! ==  D. Liu and J. Nocedal, Math. Progr. B 45 (1989) 503-528     ==
    ! ==--------------------------------------------------------------==
    ! == INPUT:                                                       ==
    ! ==   XPAR(NDIM)  Coordinates                                    ==
    ! ==   DXPAR(NDIM) Gradient                                       ==
    ! ==   NDIM        Number of degrees of freedom                   ==
    ! ==   ETOT        Total energy (for trust radius or line search) ==
    ! ==   LRET        Accepted step: return (no new step) if .T.     ==
    ! == OUTPUT:                                                      ==
    ! ==   XPAR(NDIM)  New coordinates                                ==
    ! ==   LRET        Return to get new energy / gradient if .T.     ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ndim
    REAL(real_8)                             :: dxpar(ndim), xpar(ndim), etot
    LOGICAL                                  :: lret

    EXTERNAL                                 :: ddot
    INTEGER                                  :: i, icp, icpnd, igcn, inmc, &
                                                ipt, iscn, isub, maxfev
    LOGICAL                                  :: lstpls
    REAL(real_8)                             :: alpha, beta, ddot, diag, &
                                                etotsc, grdnrm, grdsqu, &
                                                grdstp, rho, scl, sclprv, sq, &
                                                stpprd, yr

! Externals
! Variables
! ==--------------------------------------------------------------==

    CALL tiset('    RLBFGS',isub)
10  CONTINUE
    ! ==--------------------------------------------------------------==
    ! ==  All scratch space needs to be preserved (Hessian etc.)      ==
    ! ==    See lscal.inc for shared information, work arrays, etc.   ==
    ! ==    See BLOCK DATA below for initial values and explanations  ==
    ! ==--------------------------------------------------------------==
    ! ==  Organization of L-cntl%bfgs work space                           ==
    ! ==                                                              ==
    ! ==  The subroutine implements the formula given by Nocedal in   ==
    ! ==  Formula by J. Nocedal, Math. of Computation, 24, 1980, 773. ==
    ! ==  Therefore, the organization of WLBFGS is similar to the     ==
    ! ==  reference implementation by Nocedal.                        ==
    ! ==                                                              ==
    ! ==  Start    Size      Meaning                                  ==
    ! ==  -----    ----      -------                                  ==
    ! ==  1        NDIM      Work area / space for previous gradient  ==
    ! ==  NDIM+1   NREM      Storage for RHO scalars                  ==
    ! ==  2*NDIM+1 NREM      Storage for ALPHA scalars                ==
    ! ==  ISPNT    NDIM*NREM NREM performed steps                     ==
    ! ==  IGPNT    NDIM*NREM NREM gradient diffs res. from steps      ==
    ! ==  IWPNT    NDIM      Work area for line search                ==
    ! ==--------------------------------------------------------------==
    ! ==  These tests are performed every time                        ==
    ! ==  See also lscal.inc                                          ==
    ! ==--------------------------------------------------------------==
    iprnt = cprint%iprint(iprint_lscal)
    IF (nstried_l.EQ.0 .AND. iprnt.GE.1) THEN
       IF (paral%io_parent)&
            WRITE (6,*) 'ENTERING LIMITED MEMORY BFGS OPTIMIZER'
    ENDIF
    IF (nrem.EQ.0) THEN
       nrem = MIN (40, ndim)
       IF ((iprnt.GE.1).AND.paral%io_parent)&
            WRITE(6,*)&
            'SETTING THE DEFAULT NUMBER OF REMEMBERED L-BFGS STEPS'
    ENDIF
    IF (nstried_l.EQ.0 .AND. trustr.NE.0.0_real_8) THEN
       step_ini_l = trustr
       step_max_l = trustr
    ENDIF
    IF (nstried_l.EQ.0 .AND. ntrust.NE.0) THEN
       step_ini_l = step_ini_l / 5.0_real_8
    ENDIF
    IF (nstried_l.EQ.0 .AND. iprnt.GE.1) THEN
       IF (paral%io_parent)&
            WRITE (6,*) 'L-BFGS INPUT OPTIONS:'
       IF (paral%io_parent)&
            WRITE (6,'(3X,A,I5)') 'NREM: ', nrem
       IF (paral%io_parent)&
            WRITE (6,'(3X,A,I5)') 'NTRUST: ', ntrust
       IF (paral%io_parent)&
            WRITE (6,'(3X,A,I5)') 'NRESTT: ', nrestt
       IF (paral%io_parent)&
            WRITE (6,'(3X,A,I5)') 'NTRSTR: ', ntrstr
       IF (paral%io_parent)&
            WRITE (6,*)
       IF (ntrust.EQ.0) THEN
          IF (trustr.EQ.0.0_real_8) THEN
             IF (paral%io_parent)&
                  WRITE (6,*) 'L-BFGS COMPILED-IN OPTIONS:'
          ELSE
             IF (paral%io_parent)&
                  WRITE (6,*) 'L-BFGS OPTIONS BASED ON INPUT:'
          ENDIF
          IF (paral%io_parent)&
               WRITE (6,'(3X,A,F12.8)') 'STEP_INI_L: ', step_ini_l
          IF (paral%io_parent)&
               WRITE (6,'(3X,A,F12.8)') 'STEP_MIN_L: ', step_min_l
          IF (paral%io_parent)&
               WRITE (6,'(3X,A,F12.8)') 'STEP_MAX_L: ', step_max_l
          IF (paral%io_parent)&
               WRITE (6,'(3X,A,F12.8)') 'STEP_BMB_L: ', step_bmb_l
       ELSE
          IF (paral%io_parent)&
               WRITE (6,'(3X,A,F12.8)')&
               'INITIAL GUESS FOR LINE SEARCH STEP: ', step_ini_l
       ENDIF
       IF (paral%io_parent)&
            WRITE (6,*)
    ENDIF
    IF (nrestt.GT.0 .AND. nstep_l.GT.0 .AND. jump.EQ.2) THEN
       IF (MOD(nstep_l,nrestt).EQ.0) THEN
          IF ((iprnt.GE.1).AND.paral%io_parent)&
               WRITE(6,*) 'PERIODIC RESET OF L-BFGS'
          iter_l = 0
          ibnd_l = 0
          ipnt_l = 0
          npt_l  = 0
          step_l = step_ini_l
          linils = .TRUE.
       ENDIF
    ENDIF
    IF (ntrstr.GT.0 .AND. nstep_l.GT.ntrstr) ntrust = 1
    IF (ntrust.EQ.1 .AND. nstep_l.EQ.0 .AND. cntr%tolad.GT.0.0_real_8) THEN
       ntrust = 0
       ntrstr = 2
       IF ((iprnt.GE.1).AND.paral%io_parent)&
            WRITE(6,'(1X,A,A)') 'ADAPTIVE TOLERANCE: ',&
            ' USING TRUST RADIUS FOR THE FIRST TWO STEPS'
    ENDIF
    IF ((ntrust.EQ.1 .OR. ntrstr.GT.0) .AND. cntr%tolad.GT.0.0_real_8 .AND.&
         lmicro) THEN
       ntrust = 0
       ntrstr = 0
       IF ((iprnt.GE.1).AND.paral%io_parent)&
            WRITE(6,'(1X,A,/,1X,A)')&
            'MICROITERATIVE TS SEARCH AND ADAPTIVE TOLERANCE:',&
            'USING TRUST RADIUS INSTEAD OF LINE SEARCH'
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ==  Initialization  / Reset                                     ==
    ! ==  JUMP = 0 / -1                                               ==
    ! ==--------------------------------------------------------------==
    IF (jump.LE.0) THEN
       IF ((iprnt.GE.1).AND.paral%io_parent)&
            WRITE(6,*) '  JUMP=0, INITIALIZING L-BFGS'
       ispnt  = ndim+2*nrem       ! pointer to step area
       igpnt  = ispnt + ndim*nrem ! pointer to gradient area
       iwpnt  = igpnt + ndim*nrem ! pointer to line search area
       ipnt_l = 0                 ! counter for circular storage
       ibnd_l = 0                 ! number of stored steps
       iter_l = 0                 ! number of performed steps
       npt_l  = 0                 ! offset in the L-cntl%bfgs history
       CALL zeroing(wlbfgs)!, ndim+2*nrem+2*ndim*nrem)
       CALL dcopy (ndim, dxpar, 1, wlbfgs(ispnt+1), 1)
       CALL dscal (ndim, -1.0_real_8, wlbfgs(ispnt+1), 1)! DIAG = 1
       step_l = step_ini_l
       IF (.NOT.lret) jump = 2
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ==  Form an L-cntl%bfgs step                                         ==
    ! ==  JUMP = 2                                                    ==
    ! ==--------------------------------------------------------------==
    IF (jump.EQ.2) THEN
       IF ((iprnt.GE.1).AND.paral%io_parent)&
            WRITE(6,*) '  JUMP=2, FORMING L-BFGS STEP'
       iter_l = iter_l + 1
       IF (iter_l.GT.1) THEN
          ibnd_l = iter_l - 1
          IF (iter_l.GT.nrem) ibnd_l = nrem
          grdstp = ddot (ndim, wlbfgs(igpnt+npt_l+1), 1,wlbfgs(ispnt+&
               npt_l+1), 1)
          grdsqu = ddot (ndim, wlbfgs(igpnt+npt_l+1), 1,wlbfgs(igpnt+&
               npt_l+1), 1)
          diag = grdstp / grdsqu! only unit matrix
          ! ==--------------------------------------------------------------==
          ! ==  Formula by J. Nocedal, Math. of Computation, 24, 1980, 773  ==
          ! ==--------------------------------------------------------------==
          icp = ipnt_l
          IF (ipnt_l.EQ.0) icp = nrem
          rho = 1.0_real_8 / grdstp
          wlbfgs(ndim+icp) = rho
          CALL dcopy (ndim, dxpar, 1, wlbfgs, 1)
          CALL dscal (ndim, -1.0_real_8, wlbfgs, 1)
          icp = ipnt_l
          ! 
          DO i = 1,ibnd_l
             icp = icp - 1
             IF (icp.EQ.-1) icp = nrem - 1
             icpnd = icp * ndim
             sq = ddot (ndim, wlbfgs(ispnt+icpnd+1), 1, wlbfgs, 1)
             inmc = ndim + nrem + icp + 1! for alpha
             igcn = igpnt + icpnd! grad. step I-1 steps ago
             alpha = sq * wlbfgs(ndim+icp+1)
             wlbfgs(inmc) = alpha
             CALL daxpy (ndim, -alpha, wlbfgs(igcn+1), 1, wlbfgs, 1)
          ENDDO        ! I = 1,IBND_L
          ! 
          CALL dscal (ndim, diag, wlbfgs, 1)
          ! 
          DO i = 1,ibnd_l
             icpnd = icp * ndim
             yr = ddot (ndim, wlbfgs(igpnt+icpnd+1), 1, wlbfgs, 1)
             beta = yr * wlbfgs(ndim+icp+1)
             inmc = ndim + nrem + icp + 1! for beta
             beta = wlbfgs(inmc) - beta
             iscn = ispnt + icpnd! step IBND_L-I steps ago
             CALL daxpy (ndim, beta, wlbfgs(iscn+1), 1, wlbfgs, 1)
             icp = icp + 1
             IF (icp.EQ.nrem) icp = 0
          ENDDO        ! I = 1,IBND_L
          ! 
          ! ==  Store new search direction to WLBFGS                        ==
          ! 
          CALL dcopy (ndim, wlbfgs, 1, wlbfgs(ispnt+ipnt_l*ndim+1), 1)
       ENDIF                              ! IF (ITER.NE.1)
       oldene = etot                       ! save energy
       CALL dcopy (ndim, dxpar, 1, wlbfgs, 1)! save gradient
       CALL stack_sched (ielstk, st_put)   ! save electronic WF
       npt_l  = ipnt_l * ndim
       ! 
       info_l = 0
       IF (ntrust.EQ.1) THEN               ! line search
          IF (linils) THEN
             grdnrm = SQRT (ddot (ndim, dxpar, 1, dxpar, 1))
             sz_l = step_l / grdnrm! trust radius first trial step
             linils = .FALSE.
          ELSE
             sz_l = 1.0_real_8
             stpprd = SQRT (ddot (ndim, wlbfgs(ispnt+npt_l+1), 1,&
                  wlbfgs(ispnt+npt_l+1), 1))
             IF (stpprd.GT.step_max_l) sz_l = step_max_l / stpprd
          ENDIF
       ENDIF
       ! 
       jump = 3 + 10*ntrust
    ENDIF                                    ! IF (JUMP.EQ.2)
    ! ==--------------------------------------------------------------==
    ! ==  Perform an L-cntl%bfgs step - HDLCopt trust radius algorithm     ==
    ! ==  JUMP = 3                                                    ==
    ! ==--------------------------------------------------------------==
    IF (jump.EQ.3) THEN
       IF ((iprnt.GE.1).AND.paral%io_parent)&
            WRITE(6,*) '  JUMP=3, PERFORMING L-BFGS STEP'
       ipt = ispnt + npt_l ! pointer to stored new seach direction
       sz_l = SQRT (ddot (ndim, wlbfgs(ipt+1), 1, wlbfgs(ipt+1), 1))
       ! 
       ! ==  Replace stored search direction by actual step in WLBFGS    ==
       ! 
       IF (paral%io_parent.AND.iprnt.GE.1) THEN
          WRITE (6,'(5X,A,F10.5)') 'TRUST RADIUS STEP:   ', step_L
          WRITE (6,'(5X,A,F10.5)') 'PREDICTED STEP SIZE: ', sz_L
          WRITE (6,'(5X,A,F10.5)') 'PERFORMED STEP SIZE: ',MIN (sz_l,&
               step_l)
       ENDIF
       IF (sz_l.GT.step_l) THEN
          scl = step_l / sz_l
          CALL dscal (ndim, scl, wlbfgs(ipt+1), 1)! scale step
       ENDIF
       ! 
       scl1_l = ddot (ndim, dxpar, 1, wlbfgs(ipt+1), 1)
       IF (scl1_l.GT.0.0_real_8) THEN! invert step direction if upward
          IF ((iprnt.GE.1).AND.paral%io_parent)&
               WRITE(6,*) '  INVERTING L-BFGS DIRECTION'
          scl1_l = -scl1_l
          CALL dscal (ndim, -1.0_real_8, wlbfgs(ipt+1), 1)
       ENDIF
       ! 
       CALL daxpy (ndim, 1.0_real_8, wlbfgs(ipt+1), 1, xpar, 1)! do step
       ! 
       nstried_l = nstried_l + 1
       lret = .TRUE.
       jump = 4
       GOTO 999            ! return for new energy / gradient
    ENDIF                    ! IF (JUMP.EQ.3)
    ! ==--------------------------------------------------------------==
    ! ==  Test an L-cntl%bfgs step - HDLCopt trust radius algorithm        ==
    ! ==  JUMP = 4                                                    ==
    ! ==--------------------------------------------------------------==
    IF (jump.EQ.4) THEN
       IF ((iprnt.GE.1).AND.paral%io_parent)&
            WRITE(6,*) '  JUMP=4, TESTING L-BFGS STEP'
       ipt = ispnt + npt_l    ! pointer to stored performed step
       ! 
       ! ==  Failure: energy increased                                   ==
       ! 
       IF (etot.GT.oldene) THEN
          IF ((iprnt.GE.1).AND.paral%io_parent)&
               WRITE(6,'(5X,A,F12.6,A,F12.6)')&
               'ENERGY INCREASED FROM ', OLDENE, ' TO ', ETOT
          step_l = 0.5_real_8 * MIN (sz_l, step_l)
          CALL daxpy (ndim, -1.0_real_8, wlbfgs(ipt+1), 1, xpar, 1)
          CALL stack_sched (ielstk, st_get)! restore electronic WF
          ! 
          IF (step_l.LT.step_min_l .AND. nres_l.EQ.0) THEN
             nres_l = nres_l + 1! increment number of resets
             IF ((iprnt.GE.1).AND.paral%io_parent)&
                  WRITE(6,*) '    STEP TOO SMALL, RESET L-BFGS'
             jump = -1
             GOTO 10    ! go to beginning of routine
          ELSE IF (step_l.LE.step_bmb_l) THEN
             ! AK: 2005/03/28. we should not stop here. 
             ! CALL STOPGM ('RLBFGS', 'L-BFGS STEP TOO SMALL')
             IF (paral%io_parent)&
                  WRITE(6,*) ' RLBFGS| WARNING! L-BFGS STEP TOO SMALL!'
             IF (paral%io_parent)&
                  WRITE(6,*) ' RLBFGS| STOP ME IF THIS HAPPENS TOO OFTEN'
             GOTO 999
          ELSE             ! simply reject step and retry
             IF ((iprnt.GE.1).AND.paral%io_parent)&
                  WRITE(6,*) '    TRYING SMALLER STEP'
             jump = 3 + 10*ntrust
             GOTO 10    ! go to beginning of routine
          ENDIF
          ! 
          ! ==  Success: energy decreased                                   ==
          ! 
       ELSE                   ! IF (ETOT.GT.OLDENE) ...
          IF ((iprnt.GE.1).AND.paral%io_parent)&
               WRITE(6,'(5X,A,F12.6,A,F12.6)')&
               'ENERGY DECREASED FROM ', oldene, ' TO ', etot
          nres_l = 0
          nstep_l = nstep_l + 1
          ! 
          ! ==  Store gradient step to WLBFGS                               ==
          ! 
          CALL dcopy (ndim, dxpar, 1, wlbfgs(igpnt+npt_l+1), 1)
          CALL daxpy (ndim, -1.0_real_8, wlbfgs, 1,wlbfgs(igpnt+npt_l+1),&
               1)
          scl2_l = ddot (ndim, dxpar, 1, wlbfgs(ipt+1), 1)
          ! 
          ipnt_l = ipnt_l + 1      ! new entry in WLBFGS
          IF (ipnt_l.EQ.nrem) ipnt_l = 0! cycle if history full
          ! 
          IF (etot.LE.oldene+scl1_l*betap_l .AND.&
               ABS(SCL2_L).LE.ABS(SCL1_L)*BETA_L) THEN
             IF ((iprnt.GE.1).AND.paral%io_parent)&
                  WRITE(6,*) '    STEP PASSED WOLFE CONDITIONS'
             step_l = 2.0_real_8 * step_l! Wolfe conditions satisfied
          ELSE
             IF ((iprnt.GE.1).AND.paral%io_parent)&
                  WRITE(6,*) '    STEP FAILED WOLFE CONDITIONS'
             step_l = MIN (sz_l, step_l)! Wolfe cond. not satisfied
          ENDIF
          ! 
          IF (etot.LE.oldene+scl1_l*betap_l .AND. step_l.LE.sz_l)&
               THEN
             step_l = 1.25_real_8 * step_l! Small step but energy decr.
          ENDIF
          ! 
          deaccp = etot - oldene! energy change of step
          step_l = MIN (step_l, step_max_l)
          jump = 2         ! next action: form an L-cntl%bfgs step
          IF (lret) THEN
             lret = .FALSE.! we do not need new energy / gradient
             GOTO 999   ! return to check convergence
          ELSE
             GOTO 10    ! go to beginning of optimizer
          ENDIF
       ENDIF                 ! IF (ETOT.GT.OLDENE) ... ELSE ...
    ENDIF                       ! IF (JUMP.EQ.4)
    ! ==--------------------------------------------------------------==
    ! ==  Perform an L-cntl%bfgs step - line search                        ==
    ! ==  JUMP = 13                                                   ==
    ! ==--------------------------------------------------------------==
    IF (jump.EQ.13) THEN
       IF ((iprnt.GE.1).AND.paral%io_parent)&
            WRITE(6,*) '  JUMP=13, PERFORMING LINE SEARCH'
       maxfev = 4
       nstried_l = nstried_l + 1
       ! 
       IF (iprnt.GE.1) THEN
          stpprd = SQRT (ddot (ndim, wlbfgs(ispnt+npt_l+1), 1,wlbfgs(&
               ispnt+npt_l+1), 1))
          IF (paral%io_parent)&
               WRITE (6,'(5X,A,F12.6)') 'PREDICTED STEP:      ', stpprd
          IF (nfev_l.GE.1) THEN
             IF (etot.GT.oldene) THEN
                IF (paral%io_parent)&
                     WRITE (6,'(5X,A,F12.6,A,F12.6)')&
                     'ENERGY INCREASED FROM ', OLDENE, ' TO ', ETOT
             ELSE
                IF (paral%io_parent)&
                     WRITE (6,'(5X,A,F12.6,A,F12.6)')&
                     'ENERGY DECREASED FROM ', OLDENE, ' TO ', ETOT
             ENDIF
          ENDIF
       ENDIF
       ! 
       IF (info_l.EQ.0) THEN      ! remember the best cycle
          etotmn = etot        ! no step yet done
          nfevmn = -1
          sz_lmn = step_min_l/stpprd! to avoid singularities
       ELSE IF (etot.LT.etotmn) THEN
          etotmn = etot        ! lowest energy so far
          nfevmn = nfev_l
          sz_lmn = sz_l
       ENDIF
       ! 
       sclprv = sz_l
       etotsc = etot
       CALL lsrch (ndim, xpar, etotsc, dxpar, wlbfgs(ispnt+npt_l+1),&
            sz_l, maxfev, info_l, nfev_l, wlbfgs(iwpnt))
       ! 
       IF (iprnt.GE.1) THEN
          IF (paral%io_parent)&
               WRITE (6,'(5X,A,I3)') 'LINE SEARCH CYCLE: ', nfev_l
          IF (paral%io_parent)&
               WRITE (6,'(5X,A,F12.6)') 'LINE SEARCH SCALING: ', sz_l
          IF (paral%io_parent)&
               WRITE (6,'(5X,A,F12.6)') 'ENERGY CHANGE:  ', etotsc-oldene
       ENDIF
       ! 
       IF (lfevmn) THEN        ! accept step of minimum in any case
          lfevmn = .FALSE.
          info_l = 1
       ENDIF
       ! 
       lstpls = .FALSE.
       IF (info_l.EQ.3) THEN    ! maximum number of cycles reached
          IF ((iprnt.GE.1).AND.paral%io_parent)&
               WRITE(6,'(5X,A,I3)')&
               'MAXIMUM NUMBER OF LINE SEARCH CYCLES REACHED'
          lstpls = .TRUE.
       ELSE IF (info_l.EQ.-1) THEN! calculate energy / gradient
          CALL stack_sched (ielstk, st_get)! restore electronic WF
          IF (ABS(stpprd*(sz_l-sclprv)).LE.step_min_l .AND.nfev_l.GT.&
               1) THEN
             IF (iprnt.GE.1) THEN
                IF (paral%io_parent)&
                     WRITE (6,'(5X,A,A,1PE10.4)') 'MINIMUM ALLOWED ',&
                     'CHANGE OF A LINE SEARCH STEP: ', STEP_MIN_L
                IF (paral%io_parent)&
                     WRITE (6,'(5X,A,1PE10.4)') 'ACTUAL CHANGE: ',&
                     STPPRD*(SZ_L-SCLPRV)
             ENDIF
             lstpls = .TRUE.
          ELSE IF (ABS(stpprd*sz_l).LE.step_min_l .AND.nfev_l.GT.1)&
               THEN
             IF (iprnt.GE.1) THEN
                IF (paral%io_parent)&
                     WRITE (6,'(5X,A,A,1PE10.4)') 'MINIMUM ALLOWED ',&
                     'LINE SEARCH STEP: ', STEP_MIN_L
                IF (paral%io_parent)&
                     WRITE (6,'(5X,A,1PE10.4)') 'ACTUAL STEP: ',&
                     STPPRD*SZ_L
             ENDIF
             lstpls = .TRUE.
          ELSE
             jump = 13
             GOTO 999     ! return
          ENDIF
       ELSE IF (info_l.EQ.1) THEN! step ok
          nstep_l = nstep_l + 1
          jump = 14
       ELSE                     ! error condition
          IF (paral%io_parent)&
               WRITE (6,*) 'ERROR IN LINE SEARCH, CODE ', info_l
          CALL stopgm ('RLBFGS','CANNOT RECOVER',& 
               __LINE__,__FILE__)
       ENDIF
       ! 
       IF (lstpls) THEN         ! maximum number of cycles
          IF (nfevmn.LT.nfev_l-2 .AND.&
               ABS(STPPRD*(SZ_L-SZ_LMN)).GT.STEP_MIN_L) THEN
             IF ((iprnt.GE.1).AND.paral%io_parent)&
                  WRITE(6,'(5X,A,I3)')&
                  'GOING BACK TO GUESS ', NFEVMN
             sclprv = sz_l
             sz_l = sz_lmn
             CALL daxpy (ndim, sz_l-sclprv, wlbfgs(ispnt+npt_l+1),1,&
                  xpar, 1)
             lfevmn = .TRUE.
             jump = 13
             GOTO 999! return
          ELSE
             IF ((iprnt.GE.1).AND.paral%io_parent)&
                  WRITE(6,'(5X,A)')&
                  'ACCEPTING PREVIOUS GUESS'
             CALL daxpy (ndim, sclprv-sz_l, wlbfgs(ispnt+npt_l+1),&
                  1, XPAR, 1)
             nstep_l = nstep_l + 1
             jump = 14
          ENDIF
       ENDIF
       ! 
    ENDIF                         ! IF (JUMP.EQ.13)
    ! ==--------------------------------------------------------------==
    ! ==  Perform the L-cntl%bfgs step determined by the line search       ==
    ! ==  JUMP = 14                                                   ==
    ! ==--------------------------------------------------------------==
    IF (jump.EQ.14) THEN
       IF (iprnt.GE.1) THEN
          IF (paral%io_parent)&
               WRITE(6,*) '  JUMP=14, ACCEPTING STEP'
          IF (paral%io_parent)&
               WRITE (6,'(5X,A,F10.7)') 'PERFORMED STEP:      ',stpprd*sz_l
       ENDIF
       ! 
       CALL dscal (ndim, sz_l, wlbfgs(ispnt+npt_l+1), 1)
       CALL dcopy (ndim, dxpar, 1, wlbfgs(igpnt+npt_l+1), 1)
       CALL daxpy (ndim, -1.0_real_8, wlbfgs, 1, wlbfgs(igpnt+npt_l+1), 1)
       ! 
       deaccp = etot - oldene              ! energy change
       ipnt_l = ipnt_l + 1                 ! new entry in WLBFGS
       IF (ipnt_l.EQ.nrem) ipnt_l = 0      ! cycle if history full
       ! 
       jump = 2            ! unconditional: next step
       IF (lret) THEN
          lret = .FALSE.! we do not need new energy / gradient
          GOTO 999      ! return to test convergence
       ELSE
          GOTO 10       ! go to beginning of optimizer
       ENDIF
    ENDIF                    ! IF (JUMP.EQ.14)
    ! ==--------------------------------------------------------------==
999 CONTINUE                  ! RETURN
    CALL tihalt('    RLBFGS',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rlbfgs
  ! ==================================================================
  SUBROUTINE stack_init (istack, lenent, isize, ioper)
    ! ==--------------------------------------------------------------==
    ! ==  Stack routines for correct store / restore of distributed   ==
    ! ==  data controlled from routines running on one node only      ==
    ! ==                                                              ==
    ! ==  The operations are scheduled using STACK_SCHED and executed ==
    ! ==  later by STACK_DO_SCHED on all processors (which resets the ==
    ! ==  schedule). The information and pointers for the stacks are  ==
    ! ==  defined in lscal.inc                                        ==
    ! ==                                                              ==
    ! ==  Currently, overlow and underflow are allowed. This might    ==
    ! ==  change later                                                ==
    ! ==                                                              ==
    ! ==  Supported operations:                                       ==
    ! ==  ST_NOP:  no operation                                       ==
    ! ==  ST_PUSH: add current vector to the top of the stack         ==
    ! ==  ST_PUT:  top of the stack is replaced by the current vector ==
    ! ==  ST_POP:  get and remove the vector on the top of the stack  ==
    ! ==  ST_GET:  get the vector on the top of the stack             ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: istack, lenent, isize, ioper

    CHARACTER(*), PARAMETER                  :: procedureN = 'stack_init'

    INTEGER                                  :: ierr

! ==--------------------------------------------------------------==

    IF (itpstk(istack) .GT. -1) RETURN ! Stack already initialized
    IF (istack.EQ.ielstk) THEN
       iszstk(istack) = isize
       istlen(istack) = lenent
       itpstk(istack) = 0
       iopnxt(istack) = ioper
       ALLOCATE(elestk(iszstk(istack)*lenent),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE stack_init
  ! ==================================================================
  SUBROUTINE stack_sched (istack, ioper)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: istack, ioper

! ==--------------------------------------------------------------==

    iopnxt(istack) = ioper
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE stack_sched
  ! ==================================================================
  SUBROUTINE stack_do (istack, ioper, var, lenent)
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: istack, ioper, lenent
    REAL(real_8)                             :: var(lenent)

    INTEGER                                  :: ib

! Variables
! ==--------------------------------------------------------------==
! 
! ==  CATCH ERROR CONDITIONS                                      ==
! 

    IF (itpstk(istack).LT.0) THEN
       IF (paral%io_parent)&
            WRITE (6,*) 'STACK ', istack, ' UNINITIALIZED'
       CALL stopgm ('STACK_DO_SCHED', 'CANNOT RECOVER',& 
            __LINE__,__FILE__)
    ENDIF
    IF (itpstk(istack).EQ.0 .AND.(ioper.EQ.st_pop .OR. ioper.EQ.&
         st_get)) THEN
       IF (paral%io_parent)&
            WRITE (6,*) 'NO DATA ON STACK ', istack
       CALL stopgm ('STACK_DO_SCHED', 'CANNOT RECOVER',& 
            __LINE__,__FILE__)
    ENDIF
    IF (istlen(istack).NE.lenent) THEN
       IF (paral%parent) THEN
          IF (paral%io_parent)&
               WRITE (6,*) 'STACK ', istack, ' ENTRY SIZE MISMATCH'
          IF (paral%io_parent)&
               WRITE (6,*) 'allocateD:         ', istlen(istack)
          IF (paral%io_parent)&
               WRITE (6,*) 'DATA TO BE STORED: ', lenent
       ENDIF
       CALL stopgm ('STACK_DO_SCHED', 'CANNOT RECOVER',& 
            __LINE__,__FILE__)
    ENDIF
    ! 
    ! ==  DETERMINE PLACE WHERE TO PUT / GET DATA AND NEW TOP         ==
    ! 
    IF (ioper.EQ.st_put) THEN
       IF (itpstk(istack).EQ.0) itpstk(istack) = 1
       ib = (itpstk(istack)-1) * istlen(istack)
    ELSE IF (ioper.EQ.st_get) THEN
       ib = (itpstk(istack)-1) * istlen(istack)
    ELSE IF (ioper.EQ.st_push) THEN
       IF (itpstk(istack).LT.iszstk(istack)) THEN
          itpstk(istack) = itpstk(istack) + 1
       ELSE
          itpstk(istack) = 1
          IF (paral%io_parent)&
               WRITE(6,*) 'WARNING: STACK ', istack, ' OVERFLOW'
       ENDIF
       ib = (itpstk(istack)-1) * istlen(istack)
    ELSE IF (ioper.EQ.st_pop) THEN
       ib = (itpstk(istack)-1) * istlen(istack)
       IF (itpstk(istack).GT.1) THEN
          itpstk(istack) = itpstk(istack) - 1
       ELSE
          itpstk(istack) = iszstk(istack)
          IF (paral%io_parent)&
               WRITE(6,*) 'WARNING: STACK ', istack, ' UNcp_erfLOW'
       ENDIF
    ENDIF
    ! 
    ! ==  COPY THE DATA                                               ==
    ! 
    IF (ioper.EQ.st_push .OR. ioper.EQ.st_put) THEN
       IF (istack.EQ.ielstk) THEN
          CALL dcopy (lenent, var, 1, elestk(ib+1), 1)
       ENDIF
    ELSE IF (ioper.EQ.st_pop .OR. ioper.EQ.st_get) THEN
       IF (istack.EQ.ielstk) THEN
          CALL dcopy (lenent, elestk(ib+1), 1, var, 1)
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE stack_do
  ! ==================================================================
  SUBROUTINE stack_destroy (istack)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: istack

    CHARACTER(*), PARAMETER                  :: procedureN = 'stack_destroy'

    INTEGER                                  :: ierr

! ==--------------------------------------------------------------==

    IF (itpstk(istack).EQ.-1) RETURN
    IF (istack.EQ.ielstk) THEN
       DEALLOCATE(elestk,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    iopnxt(istack) = st_nop
    iszstk(istack) =  1
    istlen(istack) =  0
    itpstk(istack) = -1
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE stack_destroy
  ! ==================================================================
  SUBROUTINE stack_def
    ! ==--------------------------------------------------------------==
    ! == Block data: init variables for the wavefunction stack        ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: i

! ==--------------------------------------------------------------==

    DO i = 1,nstack
       iopnxt(i) = st_nop! Do nothing by default
       iszstk(i) =  1  ! The stacks are one entry deep
       istlen(i) =  0  ! No default entry size
       itpstk(i) = -1  ! Stacks not initialized yet
    ENDDO
    ! ==--------------------------------------------------------------==
  END SUBROUTINE stack_def
  ! ==================================================================
  SUBROUTINE give_work_lbfgs (ndim, lrlbfgs)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ndim, lrlbfgs

! ==--------------------------------------------------------------==

    IF (nrem.EQ.0) THEN
       lrlbfgs = 2*ndim + 2*MAX(40,ndim) + 2*ndim*MAX(40,ndim)
    ELSE
       lrlbfgs = 2*ndim + 2*nrem + 2*ndim*nrem
    ENDIF
    lwlbfgs = lrlbfgs ! Store for straightforward writing of RESTART
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_work_lbfgs
  ! ==================================================================
  SUBROUTINE give_scr_rlbfgs (lrlbfgs, tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lrlbfgs
    CHARACTER(len=30)                        :: tag

! ==--------------------------------------------------------------==

    lrlbfgs = 0
    tag = '0'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_rlbfgs
  ! ==================================================================
  SUBROUTINE lscal_init
    ! ==--------------------------------------------------------------==
    ! == Block data: init variables for L-cntl%bfgs, P-cntl%rfo, MI, and stack  ==
    ! ==--------------------------------------------------------------==
    CALL tol_init_shared
    CALL lbfgs_init
    CALL prfo_init
    CALL stack_def
    ! ==--------------------------------------------------------------==
  END SUBROUTINE lscal_init
  ! ==================================================================
  SUBROUTINE lbfgs_init
    ! ==--------------------------------------------------------------==
    ! == Block data: init variables for L-cntl%bfgs                        ==
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    nrem       = 0         ! Number of remembered steps
    ntrust     = 0         ! HDLCopt trust radius algorithm
    nrestt     = 0         ! Periodic restart
    ntrstr     = 0         ! Trust radius then line search
    jump       = 0         ! Optimizer status (next action)
    linils     = .TRUE.    ! Line search not yet initialized
    lwlbfgs    = 0         ! WLBFGS not yet allocated
    nres_l     = 0         ! Number of resets
    nstep_l    = 0         ! Number of performed L-cntl%bfgs steps
    nstried_l  = 0         ! Number of tried L-cntl%bfgs steps
    step_ini_l = 0.5e+0_real_8    ! Initial trust radius
    step_min_l = cntr%tolng     ! Minimum trust radius for reset
    step_bmb_l = 1.e-2_real_8*cntr%tolng    ! Minimum trust radius for error
    step_max_l = 0.5e+0_real_8    ! Maximum trust radius
    trustr     = 0.0_real_8     ! For non-default trust radius (input)
    oldene     = 0.0_real_8     ! For adaptive tol. of other optim.
    lfevmn     = .FALSE.   ! Line search status ok
    ! ==--------------------------------------------------------------==
  END SUBROUTINE lbfgs_init
  ! ==================================================================
  SUBROUTINE lsrch (ndim, point, enrg, grad, dirctn,step, maxstep,&
       info, nstep, work)
    ! ==--------------------------------------------------------------==
    ! ==  This is a common driver for line search routines.           ==
    ! ==  Hard-coded parameter ILSTYP selects the routine doing the   ==
    ! ==  work:                                                       ==
    ! ==  1. CPMD routine MINLIN in sdion.F. This is a very simple    ==
    ! ==     golden section bracketing routine.                       ==
    ! ==     This option is not very optimized and most likely will   ==
    ! ==     not improve the performance over the educated guess of   ==
    ! ==     the trust radius scheme.                                 ==
    ! ==  2. The routine used in Nocedal reference implementation     ==
    ! ==     of L-cntl%bfgs based on a routine by More and Thuente at      ==
    ! ==     Argonne. This routine finds a step which satisfies the   ==
    ! ==     Wolfe criteria for sufficient decrease and curvature.    ==
    ! ==     Please obtain the routine and a license from J. Nocedal  ==
    ! ==     directly at http://www.ece.nwu.edu/~nocedal/             ==
    ! ==                                                              ==
    ! ==  Arguments:                                                  ==
    ! ==                                                              ==
    ! ==  Return code INFO:                                           ==
    ! ==  -1    Calculate energy / gradient of the predicted step     ==
    ! ==   0    Error                                                 ==
    ! ==   1    Found an acceptable step                              ==
    ! ==   3    Maximum number of line search cycles reached          ==
    ! ==  other Error                                                 ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ndim
    REAL(real_8)                             :: point(ndim), enrg, &
                                                grad(ndim), dirctn(ndim), step
    INTEGER                                  :: maxstep, info, nstep
    REAL(real_8)                             :: work(ndim)

    CHARACTER(*), PARAMETER                  :: procedureN = 'lsrch'
    INTEGER, PARAMETER                       :: ilstyp = 1 

    INTEGER                                  :: ierr, infom
    INTEGER, SAVE                            :: nstepm
    LOGICAL, SAVE                            :: lfirst = .TRUE. 
    REAL(real_8)                             :: ftol, oldstp, xtol
    REAL(real_8), ALLOCATABLE, SAVE          :: dnrgs(:), enrgs(:), steps(:)

! ==--------------------------------------------------------------==

    IF (ilstyp .EQ. 1) THEN
       IF (lfirst) THEN
          IF (paral%io_parent)&
               WRITE (6,*) 'L-BFGS WARNING: YOU ARE USING THE BUILT-IN ',&
               ' golden section line'
          IF (paral%io_parent)&
               WRITE (6,*) '  search, see rlbfgs.F for more information'
          ALLOCATE(steps(maxstep),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(enrgs(maxstep),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(dnrgs(maxstep),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          nstepm = 0
          lfirst = .FALSE.
       ENDIF
       IF (info .EQ. 0) THEN
          nstepm = 0    ! driver forced new search direction
       ENDIF
       IF (nstepm .EQ. 0) THEN
          nstep  = 0
          oldstp = 0.0_real_8
       ELSE
          nstep = nstep + 1
          oldstp = step
       ENDIF
       nstepm = nstepm + 1
       enrgs(nstepm) = enrg
       dnrgs(nstepm) = scalar (ndim, grad, dirctn)
       IF (nstep.EQ.0 .AND. dnrgs(nstepm).GE.0.0_real_8) THEN
          IF (paral%io_parent)&
               WRITE (6,*) '    WARNING: INVERTING L-BFGS SEARCH DIRECTION'
          CALL dscal (ndim, -1.0_real_8, dirctn, 1)
          CALL dscal (nstepm, -1.0_real_8, dnrgs, 1)
       ENDIF
       IF (nstepm .EQ. 1) THEN
          steps(nstepm) = 0.0_real_8
       ELSE
          steps(nstepm) = step
       ENDIF
       CALL minlin (nstepm, steps, enrgs, dnrgs, step, infom)
       IF (infom.EQ.2) THEN
          info = -1
          step = 0.1_real_8*step! quadratic approximation invalid
          CALL daxpy (ndim, step-oldstp, dirctn, 1, point, 1)
          IF (paral%io_parent)&
               WRITE (6,*) 'LSRCH WARNING: INTERVAL TOO LARGE, DOWNSIZING',&
               ' MAXIMUM STEP LENGTH'
       ELSE IF (infom.EQ.1 .AND. nstep.GT.1) THEN
          info = 1      ! accept the step
          nstepm = 0
       ELSE IF (nstep.GE.maxstep) THEN
          ! INFO = 3            ! too many line search cycles, would be
          info = 1      ! INFO=3 but golden section -> low anyway
          nstepm = 0
       ELSE IF (infom.EQ.0 .OR. nstep.LE.1) THEN
          info = -1     ! more line search cycles required
          CALL daxpy (ndim, step-oldstp, dirctn, 1, point, 1)
       ELSE
          info = 0      ! error
       ENDIF
    ELSE IF (ilstyp .EQ. 2) THEN
       xtol = 1.0e-16_real_8      ! machine precision
       ftol = 1.0e-4_real_8
       CALL mcsrch (ndim, point, enrg, grad, dirctn,step, ftol, xtol,&
            maxstep, info, nstep, work)
    ELSE
       CALL stopgm ('WRONG LINE SEARCH ROUTINE', 'CHECK RLBFGS.F',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE lsrch
  ! ==================================================================
  SUBROUTINE mcsrch (ndim, point, enrg, grad, dirctn,step, ftol,&
       xtol, maxstep, info, nstep, work)
    ! ==--------------------------------------------------------------==
    ! ==  Please rename or remove this stub if you have obtained the  ==
    ! ==  license for routine MCSRCH from J. Nocedal                  ==
    ! ==  (http://www.ece.nwu.edu/~nocedal/)                          ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ndim
    REAL(real_8)                             :: point(ndim), enrg, &
                                                grad(ndim), dirctn(ndim), &
                                                step, ftol, xtol
    INTEGER                                  :: maxstep, info, nstep
    REAL(real_8)                             :: work(ndim)

! ==--------------------------------------------------------------==
! In the obtained MCSRCH, you can change these lines:
! DGINIT = ZERO
! DO J = 1, N
! DGINIT = DGINIT + G(J)*S(J)
! END DO
! IF (DGINIT .GE. ZERO) THEN
! WRITE(LP,15)
! 15    FORMAT(/'  THE SEARCH DIRECTION IS NOT A DESCENT DIRECTION')
! RETURN
! ENDIF
! 
! into the lines:
! DGINIT = ZERO
! DO  J = 1, N
! DGINIT = DGINIT + G(J)*S(J)
! END DO
! IF (DGINIT .GE. ZERO) THEN
! WRITE(6,15)
! 15    FORMAT('     WARNING: INVERTING L-BFGS SEARCH DIRECTION')
! DO J = 1, N
! S(J) = -S(J)
! END DO
! DGINIT = -DGINIT
! ENDIF
! ==--------------------------------------------------------------==

    CALL stopgm ('MCSRCH NOT INCLUDED', 'CHECK RLBFGS.F (LSRCH)',& 
         __LINE__,__FILE__)
  END SUBROUTINE mcsrch
  ! ==================================================================


  ! ==================================================================
  SUBROUTINE rprfo (xpar, dxpar, etot,&
       LTOLAD, GENVMX, LCVMSK, LPEXIT)
    ! ==--------------------------------------------------------------==
    ! ==  PARTITIONED RATIONAL FUNCTION OPTIMIZER                     ==
    ! ==                                                              ==
    ! ==  Reference for classical P-cntl%rfo transition state search:      ==
    ! ==  A. Banerjee, N. Adams, J. Simons, and R. Shepard,           ==
    ! ==    J. Phys. Chem. 89, 1985, 52                               ==
    ! ==                                                              ==
    ! ==  Reference for microiterative P-cntl%rfo / L-cntl%bfgs TS search:      ==
    ! ==  S.R. Billeter, A.J. Turner, W. Thiel,                       ==
    ! ==    Phys. Chem. Chem. Phys.                                   ==
    ! ==--------------------------------------------------------------==
    ! == INPUT:                                                       ==
    ! ==   XPAR(NDIM)  Coordinates                                    ==
    ! ==   DXPAR(NDIM) Gradient                                       ==
    ! ==   NDIM        Number of degrees of freedom                   ==
    ! ==   ETOT        Total energy (for trust radius or line search) ==
    ! == OUTPUT:                                                      ==
    ! ==   XPAR(NDIM)  New coordinates                                ==
    ! ==   LTOLAD      .FALSE. if full electronic convergence is req. ==
    ! ==   GENVMX      Maximum gradient component of the environment  ==
    ! ==   LCVMSK      Force the system 'not converged' if .FALSE.    ==
    ! ==   LPEXIT      Force the system 'converged' if .TRUE.         ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: xpar(cotc0%nodim,1,1), &
                                                dxpar(cotc0%nodim), etot
    LOGICAL                                  :: ltolad
    REAL(real_8)                             :: genvmx
    LOGICAL                                  :: lcvmsk, lpexit

    INTEGER                                  :: i, info, isub, j, negev, &
                                                nenv, nvar2, nzerev
    LOGICAL                                  :: lconv, lfirst, lfix, lrej, &
                                                lret, lrhe
    REAL(real_8)                             :: deact, ratio

! Externals
! Variables

    SAVE    lfix        ! during energy / gradient calc
    ! ==--------------------------------------------------------------==
    ! ==  Organization of P-cntl%rfo work space                            ==
    ! ==                                                              ==
    ! ==  The subroutine closely follows the implementation of the    ==
    ! ==  HDLC optimizer (reference see above).                       ==
    ! ==                                                              ==
    ! ==  All arrays are parts of the saved work array WPRFO.         ==
    ! ==  The pointers are in /PRFOSP/ in lscal.inc.                  ==
    ! ==                                                              ==
    ! ==  Name   Size        Meaning                                  ==
    ! ==  ----   -----       -------                                  ==
    ! ==  HESSCR (NVAR,NVAR) Current Hessian of the core              ==
    ! ==  OLDHSS (NVAR,NVAR) Previous Hessian matrix                  ==
    ! ==  VMODE  (NVAR)      Mode to be followed                      ==
    ! ==  GX     (NVAR)      Gradient in Hessian eigenspace           ==
    ! ==  OLDGX  (NVAR)      Previous GX                              ==
    ! ==  OLDG   (NVAR)      Previous gradient                        ==
    ! ==  OOLDG  (NVAR)      Old previous gradient                    ==
    ! ==  EIGVAL (NVAR)      Current Hessian eigenvalues              ==
    ! ==  OLDEIG (NVAR)      Previous Hessian eigenvalues             ==
    ! ==  EIGVEC (NVAR,NVAR) Current Hessian eigenvectors             ==
    ! ==  OLDEV  (NVAR,NVAR) Previous Hessian eigenvectors            ==
    ! ==  STEP   (NVAR)      Performed step                           ==
    ! ==--------------------------------------------------------------==
    INTEGER :: ihessc, ioldhs, ivmode, igx,    ioldgx, ioldg,  iooldg,&
         IEIGVL, IOLDEG, IEIGVC, IOLDEV, ISTEP,ierr
    ! ==--------------------------------------------------------------==
    ! ==  Organization of P-cntl%rfo scratch space (NENV=NODIM-NVAR)       ==
    ! ==                                                              ==
    ! ==  Name   Size        Meaning                                  ==
    ! ==  ----   -----       -------                                  ==
    ! ==  COORD  (NVAR)      Coordinates of the core                  ==
    ! ==  GRAD   (NVAR)      Gradient of the core                     ==
    ! ==  VARS   (NENV)      Coordinates of the environment           ==
    ! ==  DVARS  (NENV)      Gradient of the environment              ==
    ! ==  S      (NVAR)      Work array for Hessian update            ==
    ! ==  T      (NVAR)      Work array for Hessian update            ==
    ! ==  WDIAG  (3*NVAR)    Work array for diagonalization (overlap  ==
    ! ==                     with S and T) - could be larger          ==
    ! ==--------------------------------------------------------------==

    ! TODO refactor these arrays 
    REAL(real_8), ALLOCATABLE :: coord(:)
    REAL(real_8), ALLOCATABLE :: grad(:)

    REAL(real_8), ALLOCATABLE :: vars(:)
    REAL(real_8), ALLOCATABLE :: dvars(:)

    REAL(real_8), ALLOCATABLE :: s(:)
    REAL(real_8), ALLOCATABLE :: t(:)
    REAL(real_8), ALLOCATABLE :: wdiag(:)
    CHARACTER(*),PARAMETER::procedureN='RPRFO'

    ! ==--------------------------------------------------------------==
    CALL tiset('     RPRFO',isub)
    ! ==--------------------------------------------------------------==
    ! ==  These tests are performed every time                        ==
    ! ==  See also lscal.inc                                          ==
    ! ==--------------------------------------------------------------==
    iprnt = cprint%iprint(iprint_lscal)
    lfirst = (nstried_p.EQ.0 .AND. nshess.EQ.0 .AND. nstried_l.EQ.0)
    lpexit = .FALSE.
    IF (lfirst .AND. iprnt.GE.1) THEN
       IF (paral%io_parent)&
            WRITE (6,*) 'INITIALIZING TRANSITION STATE SEARCH'
    ENDIF
    IF (nvar.EQ.0) THEN
       nvar = cotc0%nodim
       IF (iprnt.GE.1) THEN
          IF (paral%io_parent)&
               WRITE (6,*) 'NO REACTION CORE DEFINED, CORE SPANS SYSTEM'
          IF (paral%io_parent)&
               WRITE (6,*) '  (ALL ', nvar, ' DEGREES OF FREEDOM)'
       ENDIF
    ENDIF
    lmicro = (nvar.LT.cotc0%nodim)
    IF (lfirst .AND. iprnt.GE.1) THEN
       IF (lmicro) THEN
          IF (paral%io_parent)&
               WRITE (6,*) '  MICROITERATIVE TRANSITION STATE SEARCH'
       ELSE
          IF (paral%io_parent)&
               WRITE (6,*) '  CONVENTIONAL TRANSITION STATE SEARCH'
       ENDIF
    ENDIF
    IF (lmicro) THEN
       IF (tolenv.EQ.0.0_real_8) THEN
          tolenv = cntr%tolng / 3.0_real_8
          IF (iprnt.GE.1) THEN
             IF (paral%io_parent)&
                  WRITE (6,'(3X,A,F10.6)')&
                  'SETTING DEFAULT TOLERANCE FOR ENVIRONMENT: ', TOLENV
          ENDIF
       ENDIF
    ENDIF
    IF (trustp.NE.0.0_real_8) THEN
       step_ini_p = trustp * 0.2_real_8
       step_max_p = trustp
    ENDIF
    IF (lfirst .AND. iprnt.GE.1) THEN
       IF (paral%io_parent)&
            WRITE (6,*) 'P-RFO INPUT OPTIONS:'
       IF (paral%io_parent)&
            WRITE (6,'(3X,A,I6)') 'NVAR:   ', nvar
       IF (paral%io_parent)&
            WRITE (6,'(3X,A,I6)') 'Mode:   ', mode
       IF (paral%io_parent)&
            WRITE (6,'(3X,A,I6)') 'NSMAXP: ', nsmaxp
       IF (modelk.EQ.1) THEN
          IF (paral%io_parent)&
               WRITE (6,'(3X,A)') 'Mode IS LOCKED FOR TS SEARCH'
       ELSE
          IF (paral%io_parent)&
               WRITE (6,'(3X,A)') 'Mode FOLLOWING FOR TS SEARCH SELECTED'
       ENDIF
       IF (lprhes) THEN
          IF (paral%io_parent)&
               WRITE (6,*) 'R/T MODES ARE PROJECTED OUT OF HESSIAN'
          IF ((lmicro).AND.paral%io_parent)&
               WRITE(6,'(A,A)') ' WARNING! THIS OPTION IS ',&
               'IGNORED BECAUSE OF PARTITIONED SYSTEM'
       ELSE
          IF (paral%io_parent)&
               WRITE (6,*) 'R/T MODES ARE NOT PROJECTED OUT OF HESSIAN'
          IF ((.NOT.lmicro).AND.paral%io_parent)&
               WRITE(6,'(A,A)') ' WARNING! BE CAREFUL ',&
               'SELECTING THE EIGENMODE'
       ENDIF
       IF ((m_hess.EQ.1).AND.paral%io_parent)&
            WRITE(6,'(3X,A)') 'HESSIAN FROM CPMD IS USED'
       IF (paral%io_parent)&
            WRITE (6,*)
       IF (trustp.EQ.0.0_real_8) THEN
          IF (paral%io_parent)&
               WRITE (6,*) 'P-RFO COMPILED-IN TRUST RADIUS SETTINGS:'
       ELSE
          IF (paral%io_parent)&
               WRITE (6,*) 'P-RFO TRUST RADIUS SETTINGS BASED ON INPUT:'
       ENDIF
       IF (paral%io_parent)&
            WRITE (6,'(3X,A,F12.8)') 'STEP_INI_P: ', step_ini_p
       IF (paral%io_parent)&
            WRITE (6,'(3X,A,F12.8)') 'STEP_MIN_P: ', step_min_p
       IF (paral%io_parent)&
            WRITE (6,'(3X,A,F12.8)') 'STEP_MAX_P: ', step_max_p
       IF (paral%io_parent)&
            WRITE (6,*)
       IF (paral%io_parent)&
            WRITE (6,*) 'P-RFO STEP ACCEPTANCE CRITERIA:'
       IF (paral%io_parent)&
            WRITE (6,'(3X,A,F12.8)') 'MINIMUM OVERLAP BETWEEN MODES:   ',&
            OMIN
       IF (paral%io_parent)&
            WRITE (6,'(3X,A,F12.8)') 'MAXIMUM ENERGY CHANGE OF A STEP: ',&
            DEMIN_P
       IF (paral%io_parent)&
            WRITE (6,'(3X,A,F12.8)') 'ENERGY CHANGE, MINIMUM RATIO:    ',&
            RMIN_P
       IF (paral%io_parent)&
            WRITE (6,'(3X,A,F12.8)') 'ENERGY CHANGE, MAXIMUM RATIO:    ',&
            RMAX_P
       IF (paral%io_parent)&
            WRITE (6,*)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ==  Test and allocate scratch space for cntl%prfo and UPDHES         ==
    ! ==--------------------------------------------------------------==
    nvar2 = nvar*nvar
    nenv  = cotc0%nodim-nvar
    ALLOCATE(coord(nvar),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(grad(nvar),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(vars(nenv),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(dvars(nenv),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(s(nvar),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(t(nvar),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(wdiag(3*nvar),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! ==  Splice the coordinates between core and environment         ==
    ! ==  No splicing needs to be done when returning here            ==
    ! ==--------------------------------------------------------------==
    IF (jump.GT.0) CALL crdco (xpar, coord, vars, dxpar, grad, dvars)
10  CONTINUE
    ! ==--------------------------------------------------------------==
    ! ==  Initialization  / Reset                                     ==
    ! ==  JUMP = 0 / -1                                               ==
    ! ==--------------------------------------------------------------==
    IF (jump.LE.0) THEN
       IF ((iprnt.GE.1).AND.paral%io_parent)&
            WRITE(6,*) '  JUMP=0, INITIALIZING P-RFO'
    ELSE IF (lrstrt) THEN
       IF ((iprnt.GE.1).AND.paral%io_parent)&
            WRITE(6,'(3X,A,I1,A)') 'JUMP=', jump,&
            ', RESUMING FROM RESTART FILE'
    ENDIF
    ! ==  ALLOCATE WORK SPACE                                         ==
    IF (jump.LE.0 .OR. lrstrt) THEN
       IF (lallhs) THEN    ! Hessian not yet allocated
          ALLOCATE(hesscr(nvar, nvar),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
          lallhs = .FALSE.
       ENDIF
       ALLOCATE(oldhss(nvar, nvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(vmode(nvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(gx(nvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(oldgx(nvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(oldg(nvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(ooldg(nvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(eigval(nvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(oldeig(nvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(eigvec(nvar, nvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(oldev(nvar, nvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(step(nvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ! TODO free
    ENDIF                       ! IF (JUMP.LE.0 .OR. LRSTRT)
    ! ==  INITIALIZE FLAGS RETURNED TO RGMOPT IF RESTART              ==
    IF (lrstrt) THEN
       IF (jump.GE.5 .AND. jump.LE.7) THEN
          ltolad = .FALSE.
       ELSE
          ltolad = .TRUE.
       ENDIF
       lrstrt = .FALSE.
    ENDIF                       ! IF (LRSTRT)
    ! ==  DONE ALLOCATING WORK SPACE                                  ==
    IF (jump.LE.0) THEN
       llock = (modelk.EQ.1)  ! Mode locking if LLOCK
       step_p = step_ini_p
       ! ==  PREPARE MICROITERATIVE OPTIMIZATION                         ==
       IF (lmicro) THEN
          IF (.NOT.lnomap .AND. lmap_p.EQ.0) THEN! Mapping to core
             ! CALL MEMORY (IP_MAP_P, NODIM, 'MAP_P')
             ALLOCATE(map_p(cotc0%nodim),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             lmap_p = cotc0%nodim
             CALL clcmap (map_p, cotc0%nodim, nvar, icore, ncore)

             ! CALL FREEM (IP_ICORE)
             DEALLOCATE(icore,STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
                  __LINE__,__FILE__)
          ENDIF
          jump = 0         ! Initialize L-cntl%bfgs for environment
          lcvmsk = .FALSE. ! CNVENV tells if we are converged
          lret = .TRUE.
          CALL rlbfgs (vars, dvars, nenv, etot, lret)
          jump = 2         ! Check if already converged
          ! ==  NO MICROITERATIVE OPTIMIZATION                              ==
       ELSE
          IF (lvlhes) THEN
             jump = 5   ! Hessian from RESTART file, go ahead
             lcvmsk = .TRUE.! Do not start if already converged
          ELSE
             jump = 7   ! Calculate Hessian
             lcvmsk = .FALSE.! We cannot be considered converged
          ENDIF
          lnomap = .TRUE.  ! All atoms map one-to-one
          ltolad = .FALSE. ! No relaxed electronic convergence
       ENDIF
       CALL crdco (xpar, coord, vars, dxpar, grad, dvars)
    ENDIF                       ! IF (JUMP.LE.0)
    ! ==--------------------------------------------------------------==
    ! ==  Check if environment converged (gradient only)              ==
    ! ==  JUMP = 2 (form L-cntl%bfgs step)                                 ==
    ! ==--------------------------------------------------------------==
    IF (jump.EQ.2) THEN
       CALL cnvenv (dvars, genvmx, lconv)
       IF (iprnt.GE.1) THEN
          IF (paral%io_parent)&
               WRITE (6,'(3X,A,F10.6)')&
               'TOLERANCE FOR ENVIRONMENT:  ', TOLENV
          IF (paral%io_parent)&
               WRITE (6,'(3X,A,F10.6)')&
               'MAXIMUM GRADIENT COMPONENT: ', GENVMX
       ENDIF
       IF (lconv) THEN        ! TOLAD is in system.h
          IF (ltolad .AND. nstried_p.NE.0 .AND. cntr%tolad.GT.0.0_real_8) THEN
             IF ((iprnt.GE.1).AND.paral%io_parent)&
                  WRITE(6,*)&
                  'RELAXING ENVIRONMENT, STANDARD ELECTRONIC CONVERGENCE'
             ltolad = .FALSE.! No relaxed electronic convergence
             lconv  = .FALSE.! Another cycle testing convergence
          ELSE IF (.NOT.lvlhes) THEN
             IF ((iprnt.GE.1).AND.paral%io_parent)&
                  WRITE(6,*)&
                  'ENVIRONMENT CONVERGED, CALCULATING HESSIAN'
             jump = 7   ! Calculate Hessian
          ELSE
             IF ((iprnt.GE.1).AND.paral%io_parent)&
                  WRITE(6,*)&
                  'ENVIRONMENT CONVERGED, OPTIMIZING CORE'
             lcvmsk = .TRUE.! We may be converged
             jump = 5   ! Form P-cntl%rfo step
          ENDIF
       ELSE
          lret = .TRUE.
          CALL rlbfgs (vars, dvars, nenv, etot, lret)
          IF (lret) GOTO 999! End of routine
          IF (.NOT.lret) GOTO 10! Beginning of routine
       ENDIF                 ! IF (LCONV) ... ELSE ...
    ENDIF                       ! IF (JUMP.EQ.2)
    ! ==--------------------------------------------------------------==
    ! ==  These states belong to environment L-cntl%bfgs                   ==
    ! ==  JUMP <> 5, 6, 8                                             ==
    ! ==--------------------------------------------------------------==
    IF (jump.NE.5.AND.jump.NE.6.AND.jump.NE.7.AND.jump.NE.8) THEN
       lret = .TRUE.
       CALL rlbfgs (vars, dvars, nenv, etot, lret)
       IF (lret) GOTO 999  ! End of routine
       IF (.NOT.lret) GOTO 10! Beginning of routine
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ==  Form P-cntl%rfo step: Diagonalize Hessian and store old values   ==
    ! ==  JUMP = 5                                                    ==
    ! ==--------------------------------------------------------------==
    IF (jump.EQ.5) THEN
       IF ((iprnt.GE.1).AND.paral%io_parent)&
            WRITE(6,*) '  JUMP=5, DIAGONALIZING HESSIAN'
       IF (iprnt.GE.2) THEN! should be 2
          DO i = 1,nvar
             IF (paral%io_parent)&
                  WRITE (6,'(5X,A,I3)') 'HESSIAN, ROW ', i
             IF (paral%io_parent)&
                  WRITE (6,'(5X,3E12.5)') (hesscr(i,j), j=1,nvar) ! transp.
          ENDDO
       ENDIF
       ! 
       ! ==  STORE THINGS IN CASE OF A REJECTED STEP (SEE ALSO JUMP=6)   ==
       ! 
       CALL dcopy (nvar,  gx,     1, oldgx,  1)
       CALL dcopy (nvar,  oldg,   1, ooldg,  1)
       CALL dcopy (nvar,  eigval, 1, oldeig, 1)
       CALL dcopy (nvar2, hesscr, 1, oldhss, 1)
       CALL dcopy (nvar2, eigvec, 1, oldev,  1)
       ! 
       ! ==  UPDATE HESSIAN OR PROJECT OUT R/T MOTION OF HESSIAN         ==
       ! 
       IF (nstep_p.GT.0) THEN
          CALL updhes (hesscr, nvar, grad, oldg, step, step_p, s, t)
       ELSE IF (lprhes) THEN
          CALL prjhes (hesscr, nvar, coord, grad, eigvec, iprnt)
       ENDIF
       ! 
       ! ==  STORE MORE THINGS BEFORE THE ACTUAL STEP                    ==
       ! 
       prvene = etot
       CALL dcopy (nvar, grad, 1, oldg, 1)
       ! 
       ! ==  DIAGONALIZE HESSIAN                                         ==
       ! 
       CALL dcopy (nvar2, hesscr, 1, eigvec, 1)
       CALL dsyev ('V', 'U', nvar, eigvec, nvar, eigval, wdiag,&
            3*NVAR, INFO)
       IF (info.NE.0) THEN
          CALL stopgm ('RPRFO', 'HESSIAN COULD NOT BE DIAGONALIZED',& 
               __LINE__,__FILE__)
       ENDIF
       ! 
       ! ==  GRADIENT IN HESSIAN EIGENSPACE                              ==
       ! 
       CALL dgemv ('T', nvar, nvar, 1.0_real_8, eigvec, nvar, grad, 1,&
            0.0_real_8, GX, 1)
       IF (iprnt.GE.1) THEN
          DO i = 1,nvar
             IF (i.LE.9 .OR. iprnt.GE.2) THEN
                IF (paral%io_parent)&
                     WRITE (6,'(5X,A,I3,A,E12.5,A,E12.5)') 'EIGENVALUE ',&
                     I, ': ', EIGVAL(I), ' GRADIENT: ', GX(I)
                IF ((iprnt.GE.2).AND.paral%io_parent)&
                     WRITE(6,'(5X,3E12.5)') (eigvec(j,i), j=1,nvar)
             ENDIF
          ENDDO
          IF ((lprhes).AND.paral%io_parent)&
               WRITE(6,'(5X,A,E8.2,A)') 'EIGENVALUES BELOW ',&
               CUTEIG, ' ARE SET TO ZERO ALONG WITH GRADIENT'
       ENDIF
       ! 
       ! ==  APPLY EIGENVALUE CUTOFF AND COUNT NEGATIVE EIGENVALUES      ==
       ! 
       negev = 0
       nzerev = 0
       DO i = 1,nvar
          IF (ABS(eigval(i)).LT.cuteig) THEN
             eigval(i) = 0.0_real_8
             nzerev = nzerev + 1
             IF (lprhes) gx(i) = 0.0_real_8
          ENDIF
          IF (eigval(i).LT.0.0_real_8) negev = negev + 1
       ENDDO
       IF ((iprnt.GE.1).AND.paral%io_parent)&
            WRITE(6,*) '    P-RFO: HESSIAN HAS ', negev,&
            ' NEGATIVE EIGENVALUE(S)'
       IF ((iprnt.GE.1 .AND. nzerev.GT.0).AND.paral%io_parent)&
            WRITE(6,*)&
            '    P-RFO: HESSIAN HAS ', NZEREV, ' ZERO EIGENVALUES'
       CALL stack_sched (ielstk, st_put)   ! save electronic WF
       ! 
       ! ==  SELECT LOWEST NON-ZERO EIGENVALUE IF PROJECTING HESSIAN     ==
       ! 
       IF (lprhes .AND. Mode.EQ.-1) THEN
          DO i = 1,nvar
             IF (ABS(eigval(i)).NE.0.0_real_8 .AND. Mode.EQ.-1)&
                  Mode = I
          ENDDO
          IF ((iprnt.GE.1).AND.paral%io_parent)&
               WRITE(6,*) '    P-RFO: SELECTED MODE ',&
               Mode
       ENDIF
       jump = 8
    ENDIF                    ! IF (JUMP.EQ.5)
    ! ==--------------------------------------------------------------==
    ! ==  Form P-cntl%rfo step: do the work                                ==
    ! ==  JUMP = 8 (return here if step is rejected in JUMP=6)        ==
    ! ==--------------------------------------------------------------==
    IF (jump.EQ.8) THEN
       IF ((iprnt.GE.1).AND.paral%io_parent)&
            WRITE(6,*) '  JUMP=8, FORMING P-RFO STEP'
       ! 
       ! ==  CALCULATE HESSIAN MODE OVERLAP (MODE: PREVIOUS -> CURRENT)  ==
       ! 
       CALL ovrlap (eigvec, eigval, vmode, nvar, Mode, omin,&
            NSTEP_P, LLOCK, LREJ, LFIX, IPRNT)
       ! 
       IF (lrej) THEN
          IF ((iprnt.GE.1).AND.paral%io_parent)&
               WRITE(6,*)&
               '  P-RFO: INSUFFICIENT OVERLAP, UNDOING PREVIOUS STEP'
          step_p = ostep_p
          sz_p   = osz_p
          prvene = oprven
          ! 
          CALL dcopy (nvar,  oldgx,  1, gx,     1)
          CALL dcopy (nvar,  oldg,   1, grad,   1)
          CALL dcopy (nvar,  ooldg,  1, oldg,   1)
          CALL dcopy (nvar,  oldeig, 1, eigval, 1)
          CALL dcopy (nvar2, oldev,  1, eigvec, 1)
          CALL dcopy (nvar2, oldhss, 1, hesscr, 1)
          ! 
          CALL daxpy (nvar, -1.0_real_8, step, 1, coord, 1)
          ! 
          step_p  = 0.5_real_8 * MIN (step_p, sz_p)
          ostep_p = step_p
          osz_p   = sz_p
          lrej    = .FALSE.
       ENDIF
       ! 
       ! ==  CALCULATE THE STEP (STEP WILL BE THE ACTUAL GEOMETRY STEP)  ==
       ! 
       CALL frmstp (eigval, eigvec, gx, nvar, step, step_p, sz_p,&
            Mode, FSTEP_P, DEPRED, IPRNT)
       ! 
       ! ==  DO THE STEP AND RETURN FOR NEW ENERGY / GRADIENT            ==
       ! 
       CALL daxpy (nvar, 1.0_real_8, step, 1, coord, 1)
       ! 
       nstried_p = nstried_p + 1
       jump = 6
       GOTO 999
    ENDIF                    ! IF (JUMP.EQ.5)
    ! ==--------------------------------------------------------------==
    ! ==  TEST P-cntl%rfo step                                             ==
    ! ==  JUMP = 6                                                    ==
    ! ==--------------------------------------------------------------==
    IF (jump.EQ.6) THEN
       IF ((iprnt.GE.1).AND.paral%io_parent)&
            WRITE(6,*) '  JUMP=6, TESTING P-RFO STEP'
       deact = etot - prvene
       ratio = deact / depred
       IF (iprnt.GE.1) THEN
          IF (paral%io_parent)&
               WRITE (6,*) '    P-RFO TRUST RADIUS TESTS'
          IF (paral%io_parent)&
               WRITE (6,'(5X,A,F14.7)') 'PREDICTED ENERGY CHANGE: ', depred
          IF (paral%io_parent)&
               WRITE (6,'(5X,A,F14.7)') 'ACTUAL ENERGY CHANGE:    ', deact
       ENDIF
       ! 
       IF ((ratio.LT.rmin_p .OR. ratio.GT.rmax_p) .AND.&
            (DEPRED.GT.DEMIN_P .OR. DEACT.GT.DEMIN_P)) THEN
          IF (step_p.GE.step_min_p) THEN
             lrej = .TRUE.
             IF (iprnt.GE.1) THEN
                IF (paral%io_parent)&
                     WRITE (6,*) '    ENERGY RATIO OUTSIDE ALLOWED RANGE'
                IF (paral%io_parent)&
                     WRITE (6,'(3X,A,F10.6)')&
                     'REJECTING STEP, REDUCING TRUST RADIUS TO ',&
                     0.5_real_8 * MIN (STEP_P, SZ_P)
             ENDIF
          ELSE          ! IF (STEP_P.GE.STEP_MIN_P)
             lrej = .FALSE.
             IF ((iprnt.GE.1).AND.paral%io_parent)&
                  WRITE(6,*) '    ENERGY RATIO BAD ',&
                  'BUT TRUST RADIUS ALREADY MINIMUM: ', STEP_P
          ENDIF
       ELSE                ! IF (RATIO UNACCETPABLE) ...
          lrej = .FALSE.
          IF ((iprnt.GE.1).AND.paral%io_parent)&
               WRITE(6,*) '    P-RFO STEP ACCEPTED'
       ENDIF              ! IF (RATIO UNACCEPTABLE) ... ELSE ...
       ! 
       ! ==  UNDO REJECTED STEP AND RETRY                                ==
       ! 
       IF (lrej) THEN
          step_p = 0.5_real_8 * MIN (step_p, sz_p)
          ! write (6,*) 'Step to be undone: ' ! remove
          ! write (6,'(3f10.5)') (step(i),i=1,nvar) ! remove
          CALL daxpy (nvar, -1.0_real_8, step, 1, coord, 1)
          CALL stack_sched (ielstk, st_get)! restore electronic WF
          jump = 8      ! RETRY STEP
          GOTO 10       ! FROM BEGINNING
       ENDIF              ! IF (LREJ)
       ! 
       ! ==  STEP ACCEPTED - ADJUST TRUST RADIUS TO AVOID REJECTED STEPS ==
       ! 
       nstep_p = nstep_p + 1
       ostep_p = step_p
       osz_p   = sz_p
       oprven  = prvene
       deaccp  = etot - prvene! energy change of accepted step
       IF (lfix) THEN
          IF (step_p.GE.step_min_p) THEN
             step_p = 0.5_real_8 * step_p
             IF ((iprnt.GE.1).AND.paral%io_parent)&
                  WRITE(6,*) '    REDUCING P-RFO TRUST RADIUS'
          ELSE
             IF ((iprnt.GE.1).AND.paral%io_parent)&
                  WRITE(6,*) '    NOT REDUCING P-RFO TRUST RADIUS'
          ENDIF
       ELSE                ! IF (LFIX) ...
          IF (ratio.GE.0.75_real_8 .AND. ratio.LE.4.0_real_8/3.0_real_8&
               .AND. SZ_P.GE.STEP_P-1.0e-2_real_8) THEN
             step_p = SQRT(2.0_real_8) * step_p
             IF (ratio.GE.0.9_real_8 .AND. ratio.LE.1.1_real_8) THEN
                step_p = SQRT(2.0_real_8) * step_p
             ENDIF
          ENDIF        ! IF (RATIO REASONABLE)
       ENDIF              ! IF (LFIX) ... ELSE ...
       step_p = MAX (step_p, step_min_p)
       step_p = MIN (step_p, step_max_p)
       IF (iprnt.GE.1) THEN
          IF (paral%io_parent)&
               WRITE (6,'(5X,A,F10.6)') 'NEW TRUST RADIUS:    ', step_p
       ENDIF
       ! 
       ! ==  DETERMINE NEXT ACTION AND GO TO BEGINNING                   ==
       ! 
       lpexit = (nstried_p.GE.nsmaxp)
       IF (lmicro) THEN
          IF ((iprnt.GE.1).AND.paral%io_parent)&
               WRITE(6,*)&
               '  ALLOWING ADAPTIVE ELECTRONIC CONVERGENCE CRITERIA'
          ltolad = .TRUE.! Allow relaxed electronic convergence
          lcvmsk = .FALSE.! CNVENV determines if we are converged
          jump = 2      ! Form L-cntl%bfgs step
       ELSE
          lcvmsk = .TRUE.! Driving code may tell we are converged
          jump = 5      ! Form P-cntl%rfo step
       ENDIF
       IF (nsvib.GT.0) THEN
          IF (MOD(nstep_p,nsvib).EQ.0) CALL pvibana
       ENDIF
       GOTO 10
    ENDIF                    ! IF (JUMP.EQ.6)
    ! ==--------------------------------------------------------------==
    ! ==  Calculate initial Hessian                                   ==
    ! ==  JUMP = 7                                                    ==
    ! ==--------------------------------------------------------------==
    IF (jump.EQ.7) THEN
       IF ((iprnt.GE.1).AND.paral%io_parent)&
            WRITE(6,*) '  JUMP=7, CALCULATING HESSIAN'
       IF (m_hess.EQ.1) THEN! Hessian from CPMD
          IF ((iprnt.GE.1).AND.paral%io_parent)&
               WRITE(6,*) '  USING HESSIAN FROM CPMD'
          lrhe = .TRUE.
          CALL hessin (xpar, lrhe)
          CALL phessco
          ! IF (LPRHES)
          !           CALL PRJHES (HESSCR, NVAR, COORD, GRAD, EIGVEC, IPRNT)
          lvlhes = .TRUE.! Hessian valid now
          lcvmsk = .TRUE.! Driving code may determine converg.
          nshess = 0
          jump = 5
          GOTO 10
       ENDIF
       IF (m_hess.EQ.0) THEN! Hessian from finite difference
          CALL inihes (hesscr, nvar, coord, grad, nshess, oldg, iprnt,&
               EPS_H, LUNDOP, IELSTK, ST_PUT, ST_GET)
          ! IF (LPRHES)
          !           CALL PRJHES (HESSCR, NVAR, COORD, GRAD, EIGVEC, IPRNT)
          nshess = nshess + 1
          IF (nshess.GT.2*nvar) lvlhes = .TRUE.! Hessian valid now
       ENDIF
       IF (lvlhes) THEN
          m_hess = 0       ! Do not read same Hessian again
          lcvmsk = .TRUE.  ! Driving code may determine converg.
          lpexit = (nsmaxp.EQ.0)! Exit if Hessian only
          nshess = 0
          IF (nsvib.NE.0) CALL pvibana
          jump = 5
       ENDIF
    ENDIF                    ! IF (JUMP.EQ.7)
    ! ==--------------------------------------------------------------==
    ! ==  Splice new coordinates from core and environment            ==
    ! ==--------------------------------------------------------------==
999 CONTINUE                  ! graceful return
    CALL crdci (xpar, coord, vars)
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    DEALLOCATE(coord,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(grad,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(vars,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(dvars,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(s,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(t,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(wdiag,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL tihalt('     RPRFO',isub)
    RETURN
  END SUBROUTINE rprfo
  ! 
  ! ==================================================================
  SUBROUTINE inihes (hesscr, nvar, coord, grad, nshess, oldg, iprnt,&
       EPS, LUNDOP, IELSTK, ST_PUT, ST_GET)
    ! ==--------------------------------------------------------------==
    ! ==  Calculates the finite difference Hessian                    ==
    ! ==                                                              ==
    ! ==  The Hessian is symmetrized on the fly to parallelize well.  ==
    ! ==  Therefore, the old gradient needs to be stored in OLDG.     ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nvar
    REAL(real_8)                             :: hesscr(nvar,nvar), &
                                                coord(nvar), grad(nvar)
    INTEGER                                  :: nshess
    REAL(real_8)                             :: oldg(nvar)
    INTEGER                                  :: iprnt
    REAL(real_8)                             :: EPS
    LOGICAL                                  :: lundop
    INTEGER                                  :: ielstk, st_put, st_get

    INTEGER                                  :: iplc, ipos
    LOGICAL, SAVE                            :: lwfstk = .FALSE.

! ==--------------------------------------------------------------==

    IF (.NOT.lwfstk) THEN
       CALL stack_sched (ielstk, st_put)! save electronic WF
       lwfstk = .TRUE.
    ENDIF
    ! 
    ! ==  STARTING FROM A PARTIALLY CALCULATED HESSIAN ON CHECKPOINT  ==
    ! 
    IF (lundop) THEN
       ipos = nshess/2 + 1
       IF ((iprnt.GE.1).AND.paral%io_parent)&
            WRITE(6,*)&
            '   HESSIAN: UNDOING DIMENSION ', ipos
       coord(ipos) = coord(ipos) + 2*eps
       nshess = nshess - 1
       lundop = .FALSE.
       RETURN
    ENDIF
    ! 
    ! ==  INITIAL CALL: PREPARE                                       ==
    ! 
    IF (nshess.EQ.0) THEN
       IF ((iprnt.GE.1).AND.paral%io_parent)&
            WRITE(6,*) '    HESSIAN: INITIAL CALL'
       coord(1) = coord(1) + eps
       CALL zeroing(hesscr)!, nvar*nvar)
       CALL stack_sched (ielstk, st_put)! save electronic WF
       lwfstk = .TRUE.
       ! 
       ! ==  OTHERWISE: CALCULATE COLUMN / ROW                           ==
       ! 
    ELSE
       iplc = MOD (nshess, 2)
       ipos = nshess/2 + iplc
       IF ((iprnt.GE.1).AND.paral%io_parent)&
            WRITE(6,*) '    HESSIAN: DIMENSION ', ipos
       ! 
       ! ==  FINITE DIFFERENCE STEP IPOS BACKWARD (ODD CALLS)            ==
       ! 
       IF (iplc.EQ.1) THEN
          CALL dcopy (nvar, grad, 1, oldg, 1)
          coord(ipos) = coord(ipos) - 2*eps! back current step
          ! 
          ! ==  FINITE DIFFERENCE STEP IPOS DONE (FORWARD - EVEN CALLS)     ==
          ! 
       ELSE
          CALL daxpy (nvar, -1.0_real_8, grad, 1, oldg, 1)
          CALL daxpy (nvar, 0.25_real_8/eps, oldg, 1,&
               HESSCR(1,IPOS), 1)
          CALL daxpy (nvar, 0.25_real_8/eps, oldg, 1,&
               HESSCR(IPOS,1), NVAR)
          coord(ipos) = coord(ipos) + eps! undo current step
          ipos = ipos + 1
          IF (ipos.LE.nvar) THEN
             coord(ipos) = coord(ipos) + eps! do next step
          ENDIF
       ENDIF
       CALL stack_sched (ielstk, st_get)! restore electronic WF
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE inihes
  ! 
  ! ==================================================================
  SUBROUTINE updhes (hesscr, nvar, grad, oldg, step, stepsz, s, t)
    ! ==--------------------------------------------------------------==
    ! ==  Updates the Hessian HESSCR of the reaction core using the   ==
    ! ==  Powell formula                                              ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nvar
    REAL(real_8)                             :: hesscr(nvar,nvar), &
                                                grad(nvar), oldg(nvar), &
                                                step(nvar), STEPSZ, S(NVAR), &
                                                T(NVAR)

    EXTERNAL                                 :: ddot
    REAL(real_8)                             :: ddot, steptt, stsqnv

! Variables
! Externals
! ==--------------------------------------------------------------==
! 
! ==  T = HESS * STEP                                             ==
! 

    CALL dsymv ('U', nvar, 1.0_real_8, hesscr, nvar, step, 1, 0.0_real_8, t, 1)
    ! 
    ! ==  T = GRAD-OLDG-T,   S = GRAD-OLDG                            ==
    ! 
    CALL dcopy (nvar, grad, 1, s, 1)
    CALL daxpy (nvar, -1.0_real_8, oldg, 1, s, 1)
    CALL daxpy (nvar, -1.0_real_8, s, 1, t, 1)
    CALL dscal (nvar, -1.0_real_8, t, 1)
    ! 
    stsqnv = 1.0_real_8 / (stepsz * stepsz)
    steptt = ddot (nvar, step, 1, t, 1) * stsqnv * stsqnv
    ! 
    ! ==  H(i,j) = H(i,j) + STSQNV * (ST(i)T(j) + ST(j)T(i))          ==
    ! ==                  - STEPTT * ST(i)ST(j)                       ==
    ! 
    CALL dsyr2 ('U', nvar, stsqnv, step, 1, t, 1, hesscr, nvar)
    CALL dsyr ('U', nvar, -steptt, step, 1, hesscr, nvar)
    ! 
    RETURN
  END SUBROUTINE updhes
  ! 
  ! ==================================================================
  SUBROUTINE prjhes (hesscr, nvar, coord, grad, sder, iprnt)
    ! ==--------------------------------------------------------------==
    ! ==  Project out translational / rotational components of the    ==
    ! ==  Hessian                                                     ==
    ! ==  Note: SDER is only work space                               ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nvar
    REAL(real_8)                             :: hesscr(nvar,nvar), &
                                                coord(nvar), grad(nvar), &
                                                SDER(NVAR,NVAR)
    INTEGER                                  :: iprnt

! ==--------------------------------------------------------------==

    IF (nvar.NE.3*ions1%nat) THEN
       IF (iprnt.GE.0) THEN
          IF (paral%io_parent)&
               WRITE (6,'(A,I4,A,I4)') ' WARNING! 3*NAT:', 3*ions1%nat,&
               ' DOES NOT MATCH NVAR:', NVAR
          IF (paral%io_parent)&
               WRITE (6,*) 'HESSIAN CLEANING MUST BE SKIPPED! '
       ENDIF
       RETURN
    ELSE
       IF ((iprnt.GE.1).AND.paral%io_parent)&
            WRITE(6,*) '    CLEANING HESSIAN'
    ENDIF
    ! 
    CALL dcopy (nvar*nvar, hesscr, 1, sder, 1)
    CALL purgeh (sder, hesscr, tau0)
    CALL dcopy (nvar*nvar, sder, 1, hesscr, 1)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE prjhes
  ! 
  ! ==================================================================
  SUBROUTINE ovrlap (eigvec, eigval, vmode, nvar, Mode, omin,&
       NSTEP_P, LLOCK, LREJ, LFIX, IPRNT)
    ! ==--------------------------------------------------------------==
    ! ==  Calculate the overlap between the mode to be followed and   ==
    ! ==  the current eigenvectors of the Hessian                     ==
    ! ==  Determine the mode MODCUR to be followed                    ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nvar
    REAL(real_8)                             :: vmode(nvar), eigval(nvar), &
                                                eigvec(nvar,nvar)
    INTEGER                                  :: Mode
    REAL(real_8)                             :: omin
    INTEGER                                  :: nstep_p
    LOGICAL                                  :: llock, lrej, lfix
    INTEGER                                  :: iprnt

    EXTERNAL                                 :: ddot
    INTEGER                                  :: i, newmod
    REAL(real_8)                             :: ddot, ovlp, ovlpmx

! Variables
! Externals
! ==--------------------------------------------------------------==

    newmod = Mode
    lrej = .FALSE.
    lfix = .FALSE.
    ! 
    ! ==  INITIAL CALL: STORE MODE IN VMODE                           ==
    ! 
    IF (nstep_p.EQ.0) THEN
       IF (iprnt.GE.1) THEN
          IF (paral%io_parent)&
               WRITE (6,*) '    INITIAL HESSIAN Mode: ', mode
          IF (paral%io_parent)&
               WRITE (6,'(5X,A,F10.6)') 'HESSIAN EIGENVALUE: ',eigval(Mode)
          IF (llock) THEN
             IF (paral%io_parent)&
                  WRITE (6,*) '    THIS Mode IS LOCKED'
             DO i = 1,Mode
                IF (paral%io_parent)&
                     WRITE (6,'(5X,A,I1,A,F10.6)') 'HESSIAN EV(', i, '): ',&
                     EIGVAL(I)                  
             ENDDO
          ELSE
             IF (paral%io_parent)&
                  WRITE (6,'(5X,A,A,F10.5)') 'Mode FOLLOWING IS SWITCHED ',&
                  'ON, MINIMUM OVERLAP: ', OMIN
          ENDIF
       ENDIF
       CALL dcopy (nvar, eigvec(1,newmod), 1, vmode, 1)
    ELSE
       IF (llock) THEN
          ovlp = ABS (ddot (nvar, eigvec(1,newmod), 1, vmode, 1))
          IF ((iprnt.GE.1).AND.paral%io_parent)&
               WRITE(6,*) '    Mode LOCKED, OVERLAP: ',&
               OVLP
       ELSE
          ! 
          ! ==  DETERMINE WHICH MODE HAS THE GREATEST OVERLAP               ==
          ! 
          ovlpmx = 0.0_real_8
          DO i = 1,nvar
             ovlp = ABS (ddot (nvar, eigvec(1,i), 1, vmode, 1))
             IF (ovlp.GT.ovlpmx) THEN
                newmod = i
                ovlpmx = ovlp
             ENDIF
          ENDDO
          IF (iprnt.GE.1) THEN
             IF (paral%io_parent)&
                  WRITE (6,'(5X,A,I3,A,F10.6)') 'OVERLAP OF CURRENT Mode ',&
                  NEWMOD, ' WITH PREVIOUS Mode IS ', OVLPMX
          ENDIF
          ! 
          ! ==  CHECK IF OVERLAP IS SUFFICIENT                              ==
          ! 
          IF (ovlpmx.LT.omin) THEN
             lrej = .TRUE.
             IF ((iprnt.GE.1).AND.paral%io_parent)&
                  WRITE(6,*)&
                  '    THIS IS INSUFFICIENT, REJECTING STEP'
          ELSE IF (ovlpmx.LT.SQRT(omin)) THEN
             lfix = .TRUE.
             IF ((iprnt.GE.1).AND.paral%io_parent)&
                  WRITE(6,*)&
                  '    THIS IS TIGHT, DENYING STEP SIZE TO INCREASE'
          ENDIF
          ! 
          ! ==  STORE NEW EIGENMODE IN VMODE                                ==
          ! 
          CALL dcopy (nvar, eigvec(1,newmod), 1, vmode, 1)
          IF (iprnt.GE.1) THEN
             IF (newmod.NE.Mode) THEN
                IF (paral%io_parent)&
                     WRITE (6,*) '  WARNING: SWITCHING Mode FROM ', mode,&
                     ' TO ', NEWMOD
             ENDIF
             DO i = 1,MAX(newmod,Mode)
                IF (paral%io_parent)&
                     WRITE (6,'(5X,A,I1,A,F10.6)') 'HESSIAN EV(', i, '): ',&
                     EIGVAL(I)                  
             ENDDO
          ENDIF
       ENDIF              ! IF (NSTEP_P.GE.1) ... ELSE ...
    ENDIF
    IF (.NOT.lrej) Mode = newmod
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE ovrlap
  ! 
  ! ==================================================================
  SUBROUTINE frmstp (eigval, eigvec, gx, nvar, step, step_p, sz_p,&
       Mode, FSTEP, DEPRED, IPRNT)
    ! ==--------------------------------------------------------------==
    ! ==  Form a P-cntl%rfo step using an already diagonalized Hessian     ==
    ! ==  and predict energy change                                   ==
    ! ==                                                              ==
    ! ==  Arguments:                                                  ==
    ! ==  EIGVAL: (in)  Eigenvalues of the Hessian                    ==
    ! ==  EIGVEC: (in)  Eigenvectors of the Hessian                   ==
    ! ==  GX:     (in)  Gradient in the Hessian eigenspace            ==
    ! ==  STEP:   (out) Calculated step                               ==
    ! ==  STEP_P: (in)  Trust radius for P-cntl%rfo                        ==
    ! ==  SZ_P:   (out) Actual step size                              ==
    ! ==  MODE:   (in)  Mode to be maximized                          ==
    ! ==  FSTEP:  (par) Step size to for finding lambdas              ==
    ! ==  DEPRED: (out) Predicted energy change                       ==
    ! ==  IPRNT:  (in)  Print level                                   ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nvar
    REAL(real_8)                             :: gx(nvar), eigvec(nvar,nvar), &
                                                eigval(nvar), step(nvar), &
                                                STEP_P, SZ_P
    INTEGER                                  :: Mode
    REAL(real_8)                             :: FSTEP, DEPRED
    INTEGER                                  :: iprnt

    EXTERNAL                                 :: ddot
    INTEGER                                  :: i, modmax, modmin, ncnt
    LOGICAL                                  :: ldoscl, llow, lupp
    REAL(real_8) :: big, blow, bstep, bupp, ddot, eigmax, eigmin, eps, eps2, &
      EPS6, FLAM, FLOW, FUPP, LMBMAX, LMBMIN, SCALE, SSZMAX, SSZMIN, SSZTOL, &
      TEMP, XLMMAX, XLMMIN

! Variables
! Externals
! ==--------------------------------------------------------------==

    ldoscl = .FALSE.
    modmax = Mode
    modmin = 1
    IF (modmax.EQ.1) modmin = 2
    eigmax = eigval(modmax)
    eigmin = eigval(modmin)
    ! 
    ! ==  SET BISECTION STEP SIZES AND TOLERANCES                     ==
    ! 
    big    = 1.0e+3_real_8
    eps    = 1.0e-12_real_8
    sszmin = MAX (ABS(eigmin)*eps, 10.0_real_8*eps)
    sszmax = MAX (big, ABS(eigmin)) * big
    ssztol = 1.0e-8_real_8
    eps2   = 1.0e-2_real_8           ! for steps along eigenmodes
    eps6   = 1.0e-6_real_8           ! for determining if scaling large steps
    ! write (6,*) 'fstep, sszmin, sszmax: ', fstep, sszmin, sszmax
    ! 
    ! ==  FIND LMBMAX (MAXIMIZE MODMAX), LMBMIN (MINIMIZE ALL OTHERS) ==
    ! 
    lmbmax = 0.5_real_8 * (eigmax + SQRT (eigmax**2 + 4.0_real_8*gx(modmax)**2))
    ! 
    lmbmin = 0.0_real_8
    bstep  = fstep
    IF (eigmin.LE.0.0_real_8) lmbmin = eigmin - bstep
    IF (eigmin.GT.0.0_real_8) bstep  = eigmin
    blow = lmbmin - bstep
    bupp = lmbmin + 0.5_real_8*bstep
    ! 
    ! ==  BISECT THE INTERVALS BLOW / BUP (WHILE FLOW*FUPP >= 0)      ==
    ! 
50  flow = 0.0_real_8
    fupp = 0.0_real_8
    DO i = 1,nvar
       flow = flow + gx(i)**2/(blow-eigval(i))
       fupp = fupp + gx(i)**2/(bupp-eigval(i))
    ENDDO
    flow = flow - gx(modmax)**2/(blow-eigval(modmax))
    fupp = fupp - gx(modmax)**2/(bupp-eigval(modmax))
    flow = flow - blow
    fupp = fupp - bupp
    IF (flow*fupp.GE.0.0_real_8) THEN
       blow = blow - (eigmin-blow)
       bupp = bupp + 0.5_real_8 * (eigmin-bupp)
       llow = .FALSE.
       IF (blow.LE.-sszmax) THEN
          blow = -sszmax
          llow = .TRUE.
       ENDIF
       lupp = .FALSE.
       IF (ABS(eigmin-bupp).LE.sszmin) THEN
          bupp = eigmin - sszmin
          lupp = .TRUE.
       ENDIF
       IF (llow.AND.lupp) THEN
          IF (paral%io_parent)&
               WRITE (6,*)'  P-RFO: FIXED STEP SIZE NOT IMPLEMENTED'
          CALL stopgm ('RPRFO', 'PROBLEMS BRACKETING LAMBDA',& 
               __LINE__,__FILE__)
       ENDIF
       GOTO 50
    ENDIF                    ! WHILE (FLOW*FUPP.GE.0.0_real_8)
    ! 
    ncnt   = 0
    xlmmin = 0.0_real_8
    ! 
    ! ==  BISECT LMBMIN / BLOW / BUPP (WHILE LMBMIN NOT CONVERGED)    ==
    ! 
80  CONTINUE                  ! WHILE (ABS(XLMMIN-LMBMIN).GT.SSZTOL)
    flow = 0.0_real_8
    fupp = 0.0_real_8
    flam = 0.0_real_8
    lmbmin = 0.5_real_8 * (blow+bupp)
    DO i = 1,nvar
       flow = flow + gx(i)**2/(blow-eigval(i))
       fupp = fupp + gx(i)**2/(bupp-eigval(i))
       flam = flam + gx(i)**2/(lmbmin-eigval(i))
    ENDDO
    flow = flow - gx(modmax)**2/(blow-eigval(modmax))
    fupp = fupp - gx(modmax)**2/(bupp-eigval(modmax))
    flam = flam - gx(modmax)**2/(lmbmin-eigval(modmax))
    flow = flow - blow
    fupp = fupp - bupp
    flam = flam - lmbmin
    IF (ABS(xlmmin-lmbmin).GE.ssztol) THEN
       ncnt = ncnt + 1
       IF (ncnt.GT.1000) THEN
          IF (paral%io_parent)&
               WRITE (6,*) '  P-RFO: TOO MANY ITERATIONS FOR LAMBDA'
          CALL stopgm ('RPFRO', 'PROBLEM BISECTING LAMBDA',& 
               __LINE__,__FILE__)
       ENDIF
       xlmmin = lmbmin
       IF (flam*fupp.LT.0.0_real_8) blow = lmbmin
       IF (flam*flow.LT.0.0_real_8) bupp = lmbmin
       GOTO 80
    ENDIF
    CONTINUE                  ! WHILE (ABS(XLMMIN-LMBMIN).GT.SSZTOL)
    ! 
    ! ==  NOW WE HAVE LMBMAX, LMBMIN: PRINT IF REQUESTED              ==
    ! 
100 xlmmax = lmbmax
    xlmmin = lmbmin
    IF (iprnt.GE.1) THEN
       IF (paral%io_parent)&
            WRITE (6,'(5X,A,F14.5)')&
            'P-RFO: LAMBDA MAXIMIZING TS MODE:         ', LMBMAX
       IF (paral%io_parent)&
            WRITE (6,'(5X,A,F14.5)')&
            'P-RFO: LAMBDA MINIMIZING ALL OTHER MODES: ', LMBMIN
    ENDIF
    ! 
    ! ==  CALCULATE THE P-cntl%rfo STEP NOW                                ==
    ! 
    CALL zeroing(step)!, nvar)
    DO i = 1,nvar
       IF (i.EQ.modmax) THEN
          temp = gx(i) / (lmbmax-eigval(i))
       ELSE IF (lmbmin.EQ.0.0_real_8 .AND. ABS(eigval(i)).LT.eps2) THEN
          temp = 0.0_real_8
       ELSE
          temp = gx(i) / (lmbmin-eigval(i))
       ENDIF
       CALL daxpy (nvar, temp, eigvec(1,i), 1, step, 1)
    ENDDO
    ! 
    ! ==  SCALE THE STEP                                              ==
    ! 
    sz_p = SQRT (ddot (nvar, step, 1, step, 1))
    IF (sz_p .LT. step_p+eps6) THEN
       IF (iprnt.GE.1) THEN
          IF (paral%io_parent)&
               WRITE (6,'(5X,A,F10.7)') 'TRUST RADIUS STEP:   ', step_p
          IF (paral%io_parent)&
               WRITE (6,'(5X,A,F10.7)') 'PERFORMED STEP SIZE: ', sz_p
       ENDIF
       scale = 1.0_real_8
    ELSE
       IF (ldoscl) THEN
          IF (iprnt.GE.1) THEN
             IF (paral%io_parent)&
                  WRITE (6,'(5X,A,F10.7)') 'TRUST RADIUS STEP:   ', step_p
             IF (paral%io_parent)&
                  WRITE (6,'(5X,A,F10.7)') 'PREDICTED STEP SIZE: ', sz_p
             IF (paral%io_parent)&
                  WRITE (6,'(5X,A,F10.7)') 'PERFORMED STEP SIZE: ', step_p
          ENDIF
          scale = step_p / sz_p
          CALL dscal (nvar, scale, step, 1)
       ELSE
          IF ((iprnt.GE.1).AND.paral%io_parent)&
               WRITE(6,'(5X,A,F12.7,/,5X,A)')&
               'CANNOT SCALE STEP OF SIZE ', SZ_P,&
               'FINDING LAMBDAS FOR TRUST RADIUS'
          ldoscl = .TRUE.
          GOTO 140
       ENDIF
    ENDIF
    ! 
    ! ==  PREDICT THE ENERGY CHANGE                                   ==
    ! 
    depred = 0.0_real_8
    DO i = 1,nvar
       IF (i.EQ.modmax) THEN
          temp = scale * gx(i) / (lmbmax-eigval(i))
       ELSE IF (lmbmin.EQ.0.0_real_8 .AND. ABS(eigval(i)).LT.eps2) THEN
          temp = 0.0_real_8
       ELSE
          temp = scale * gx(i) / (lmbmin-eigval(i))
       ENDIF
       depred = depred + temp*gx(i) + 0.5_real_8*temp*temp*eigval(i)
    ENDDO
    RETURN
    ! 
    ! ==  FIND LAMBDAS FOR STEP WITH TRUST RADIUS                     ==
    ! 
140 CONTINUE
    lmbmin = 0.0_real_8
    bstep  = fstep
    IF (eigmin.LE.0.0_real_8) lmbmin = eigmin - bstep
    IF (eigmin.GT.0.0_real_8) bstep  = eigmin
    blow = lmbmin - bstep
    bupp = lmbmin + 0.5_real_8*bstep
    ! 
    ! ==  BISECTION (TRUST RADIUS STEP SIZE) - SEE LOOP 50            ==
    ! 
150 flow = 0.0_real_8
    fupp = 0.0_real_8
    DO i = 1,nvar
       flow = flow + (gx(i)/(blow-eigval(i)))**2
       fupp = fupp + (gx(i)/(bupp-eigval(i)))**2
    ENDDO
    flow = flow - (gx(modmax)/(blow-eigval(modmax)))**2&
         + (GX(MODMAX)/(BLOW+EIGVAL(MODMAX)))**2
    fupp = fupp - (gx(modmax)/(bupp-eigval(modmax)))**2&
         + (GX(MODMAX)/(BUPP-EIGVAL(MODMAX)))**2
    flow = flow - step_p**2
    fupp = fupp - step_p**2
    ! write (6,*)'before 50 loop, flow, fupp, llow, lupp, blow, bupp: ',
    !        flow, fupp, llow, lupp, blow, bupp ! remove
    IF (flow*fupp.GE.0.0_real_8) THEN
       blow = blow - (eigmin-blow)
       bupp = bupp + 0.5_real_8 * (eigmin-bupp)
       llow = .FALSE.
       IF (blow.LE.-sszmax) THEN
          blow = -sszmax
          llow = .TRUE.
       ENDIF
       lupp = .FALSE.
       IF (ABS(eigmin-bupp).LE.sszmin) THEN
          bupp = eigmin - sszmin
          lupp = .TRUE.
       ENDIF
       IF (llow.AND.lupp) THEN
          IF (paral%io_parent)&
               WRITE (6,*)'  FIXED STEP SIZE NOT IMPLEMENTED, STOPPING NOW'
          CALL stopgm ('RPRFO', 'PROBLEMS BRACKETING LAMBDA',& 
               __LINE__,__FILE__)
       ENDIF
       ! write (6,*)'end 50 loop, flow, fupp, llow, lupp, blow, bupp: ',
       !        flow, fupp, llow, lupp, blow, bupp ! remove
       GOTO 150
    ENDIF                    ! WHILE (FLOW*FUPP.GE.0.0_real_8)
    ! 
    ncnt   = 0
    xlmmin = 0.0_real_8
    ! 
    ! ==  BISECTION PART 2 (TRUST RADIUS STEP SIZE) - SEE LOOP 80     ==
    ! 
180 CONTINUE                  ! WHILE (ABS(XLMMIN-LMBMIN).GT.SSZTOL)
    flow = 0.0_real_8
    fupp = 0.0_real_8
    flam = 0.0_real_8
    lmbmin = 0.5_real_8 * (blow+bupp)
    DO i = 1,nvar
       flow = flow + (gx(i)/(blow-eigval(i)))**2
       fupp = fupp + (gx(i)/(bupp-eigval(i)))**2
       flam = flam + (gx(i)/(lmbmin-eigval(i)))**2
    ENDDO
    flow = flow - (gx(modmax)/(blow-eigval(modmax)))**2&
         + (GX(MODMAX)/(BLOW+EIGVAL(MODMAX)))**2
    fupp = fupp - (gx(modmax)/(bupp-eigval(modmax)))**2&
         + (GX(MODMAX)/(BUPP+EIGVAL(MODMAX)))**2
    flam = flam - (gx(modmax)/(lmbmin-eigval(modmax)))**2&
         + (GX(MODMAX)/(LMBMIN+EIGVAL(MODMAX)))**2
    flow = flow - step_p**2
    fupp = fupp - step_p**2
    flam = flam - step_p**2
    ! write(6,*)'within 80 loop, flow, flam, fupp, lmbmin, xlmmin: ',
    !     flow, flam, fupp, lmbmin, xlmmin ! remove
    IF (ABS(xlmmin-lmbmin).GE.ssztol) THEN
       ncnt = ncnt + 1
       IF (ncnt.GT.1000) THEN
          IF (paral%io_parent)&
               WRITE (6,*) '  P-RFO: TOO MANY ITERATIONS FOR LAMBDA'
          CALL stopgm ('RPFRO', 'PROBLEM BISECTING LAMBDA',& 
               __LINE__,__FILE__)
       ENDIF
       xlmmin = lmbmin
       IF (flam*fupp.LT.0.0_real_8) blow = lmbmin
       IF (flam*flow.LT.0.0_real_8) bupp = lmbmin
       GOTO 180
    ENDIF
    CONTINUE                  ! WHILE (ABS(XLMMIN-LMBMIN).GT.SSZTOL)
    ! 
    lmbmax = -lmbmin
    ! 
    ! ==  FOUND LAMBDAS, DOING THE STEP NOW                           ==
    ! 
    GOTO 100
    ! ==--------------------------------------------------------------==
  END SUBROUTINE frmstp
  ! ==================================================================
  SUBROUTINE crdco (xpar, coord, vars, dxpar, grad, dvars)
    ! ==--------------------------------------------------------------==
    ! ==  Splices the coordinates and gradient into a core and an     ==
    ! ==  environment part (check out)                                ==
    ! ==--------------------------------------------------------------==
    REAL(real_8) :: xpar(cotc0%nodim), coord(nvar), vars(cotc0%nodim-nvar), &
      DXPAR(cotc0%nodim), GRAD(NVAR), DVARS(cotc0%nodim-NVAR)

    INTEGER                                  :: i, j

    IF (nvar .EQ. cotc0%nodim) THEN
       CALL dcopy (nvar, xpar,  1, coord, 1)
       CALL dcopy (nvar, dxpar, 1, grad,  1)
    ELSE IF (lnomap) THEN
       CALL dcopy (nvar, xpar(1),  1, coord, 1)
       CALL dcopy (nvar, dxpar(1), 1, grad, 1)
       CALL dcopy (cotc0%nodim-nvar, xpar(nvar+1),  1, vars,  1)
       CALL dcopy (cotc0%nodim-nvar, dxpar(nvar+1), 1, dvars, 1)
    ELSE
       DO i = 1,nvar
          coord(i)  = xpar(map_p(i))
          grad(i)   = dxpar(map_p(i))
       ENDDO
       j = nvar + 1
       DO i = 1,cotc0%nodim-nvar
          vars(i)  = xpar(map_p(j))
          dvars(i) = dxpar(map_p(j))
          j = j + 1
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE crdco
  ! 
  ! ==================================================================
  SUBROUTINE crdci (xpar, coord, vars)
    ! ==--------------------------------------------------------------==
    ! ==  Splices the coordinates from the core and the environement  ==
    ! ==  into the new common coordinates (check in)                  ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: xpar(cotc0%nodim), &
                                                coord(nvar), &
                                                vars(cotc0%nodim-nvar)

    INTEGER                                  :: i, j

    IF (nvar .EQ. cotc0%nodim) THEN
       CALL dcopy (nvar, coord, 1, xpar, 1)
    ELSE IF (lnomap) THEN
       CALL dcopy (nvar, coord, 1, xpar(1), 1)
       CALL dcopy (cotc0%nodim-nvar, vars, 1, xpar(nvar+1), 1)
    ELSE
       DO i = 1,nvar
          xpar(map_p(i)) = coord(i)
       ENDDO
       j = nvar + 1
       DO i = 1,cotc0%nodim-nvar
          xpar(map_p(j)) = vars(i)
          j = j + 1
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE crdci
  ! 
  ! ==================================================================
  SUBROUTINE cnvenv (dvars, grmax, lconv)
    ! ==--------------------------------------------------------------==
    ! ==  Tests if the environment is converged                       ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: dvars(cotc0%nodim-nvar), grmax
    LOGICAL                                  :: lconv

    INTEGER                                  :: i

    grmax = 0.0_real_8
    DO i = 1,cotc0%nodim-nvar
       IF (ABS(dvars(i)) .GT. grmax) THEN
          grmax = ABS(dvars(i))
       ENDIF
    ENDDO
    lconv = (grmax.LT.tolenv)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE cnvenv
  ! 
  ! ==================================================================
  ! 
  ! ==================================================================
  SUBROUTINE give_work_prfo (ndim, lrlbfgs, lrprfo)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ndim, lrlbfgs, lrprfo

! ==--------------------------------------------------------------==

    IF (nvar.EQ.0) THEN         ! set NVAR to NDIM later at init
       lrprfo  = 4*ndim*ndim + 8*ndim
       lrlbfgs = 0
    ELSE IF (nvar.EQ.ndim) THEN ! Reaction core spans entire system
       lrprfo  = 4*nvar*nvar + 8*nvar
       lrlbfgs = 0
    ELSE IF (nvar.LT.ndim) THEN ! Reaction core spans a subsystem
       lrprfo  = 4*nvar*nvar + 8*nvar
       CALL give_work_lbfgs (ndim-nvar, lrlbfgs)
    ELSE                        ! More variables than DOF
       IF (paral%io_parent)&
            WRITE (6,*) 'INVALID INPUT FOR TS SEARCH IN REACTION CORE'
       IF (paral%io_parent)&
            WRITE (6,*) 'NUMBER OF DEGREES OF FREEDOM OF SYSTEM: ', ndim
       IF (paral%io_parent)&
            WRITE (6,*) 'NUMBER OF VARIABLES FOR TS SEARCH:      ', nvar
       CALL stopgm ('RPRFO', 'NUMBER OF DOF MISMATCH',& 
            __LINE__,__FILE__)
    ENDIF
    IF (.NOT.lallhs) lrprfo = lrprfo - nvar*nvar ! Allocated in RV30
    lwprfo = lrprfo ! Store for straightforward writing of RESTART
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_work_prfo
  ! ==================================================================
  SUBROUTINE give_scr_rprfo (lrprfo, tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lrprfo
    CHARACTER(len=30)                        :: tag

! ==--------------------------------------------------------------==

    IF (nvar.EQ.0) THEN
       lrprfo = 3*cotc0%nodim + 2*cotc0%nodim
    ELSE
       lrprfo = 3*nvar + 2*cotc0%nodim
    ENDIF
    tag = '2*NVAR+2*NODIM'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_rprfo
  ! ==================================================================
  SUBROUTINE prfo_init
    ! ==--------------------------------------------------------------==
    ! == Block data: init variables for P-cntl%rfo and MI TS search        ==
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    Mode       = 1         ! Hessian eigenmode to be followed
    modelk     = 0         ! Mode following by default
    nsmaxp     = -1        ! Maximum number of P-cntl%rfo steps
    nvar       = 0         ! Number of degrees of freedom in core
    nsvib      = 0         ! Never do vibrational analysis
    m_hess     = 0         ! Finite-difference Hessian (INIHES)
    cuteig     = 1.0e-9_real_8    ! Cutoff for Hessian eigenvalues
    demin_p    = 2.0e-5_real_8    ! Maximum energy change of a step
    omin       = 0.5_real_8     ! Minimum overlap of Hessian Mode
    rmax_p     = 1.0e+1_real_8    ! Maximum ratio for predicted energy
    rmin_p     = 0.0_real_8     ! Minimum ratio for predicted energy
    step_ini_p = 0.4e-1_real_8    ! Initial trust radius
    step_max_p = 0.2_real_8     ! Maximum trust radius
    step_min_p = 1.0e-5_real_8    ! Minimum trust radius
    tolenv     = 0.0_real_8     ! Tolerance for the environment
    trustp     = 0.0_real_8     ! Input trust radius for P-cntl%rfo
    eps_h      = 2.0e-2_real_8    ! Finite-difference displacement (Hessian)
    ! Shared information
    lallhs     = .TRUE.    ! Allocate Hessian in RGMOPT
    lmicro     = .FALSE.   ! Conventional TS search
    lnomap     = .TRUE.    ! Default map for core / environment
    lvlhes     = .FALSE.   ! No valid Hessian yet
    lprhes     = .FALSE.   ! Do not project R/T of Hessian
    lrstrt     = .FALSE.   ! Did not start from RESTART
    lundop     = .FALSE.   ! Do not undo the last Hessian step
    lwprfo     = 0         ! No work space allocated yet
    lmap_p     = 0         ! No core atoms map allocated yet
    nstep_p    = 0         ! Number of performed P-cntl%rfo steps
    nstried_p  = 0         ! Number of tried P-cntl%rfo steps
    nshess     = 0         ! Current step for Hessian calculation
    ! ==--------------------------------------------------------------==
  END SUBROUTINE prfo_init
  ! ==================================================================
  SUBROUTINE pvibana
    ! ==--------------------------------------------------------------==
    ! Variables
    CHARACTER(*), PARAMETER                  :: procedureN = 'pvibana'

    INTEGER                                  :: ierr, nd
    REAL(real_8), ALLOCATABLE                :: sder(:,:)

! ==--------------------------------------------------------------==

    nd = 3*ions1%nat
    ! CALL MEMORY (IP_HESS, ND*ND, 'HESS')
    ALLOCATE(hess(nd,nd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! CALL MEMORY (IP_SDER, ND*ND, 'SDER')
    ALLOCATE(sder(nd,nd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    IF (cnti%npara.LT.2) cnti%npara = 2
    CALL phessci (hess, nd)
    CALL dcopy (nd*nd, hess, 1, sder, 1)
    CALL vibana (sder)

    ! CALL FREEM (IP_SDER)
    DEALLOCATE(sder,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! CALL FREEM (IP_HESS)
    DEALLOCATE(hess,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE pvibana
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE phessco
    ! ==--------------------------------------------------------------==
    ! Variables
    INTEGER                                  :: i, im, j, jm

! ==--------------------------------------------------------------==

    IF (.NOT.paral%parent) RETURN
    IF (lnomap) THEN
       DO i = 1,nvar
          CALL dcopy (nvar, hess(1,i), 1, hesscr(1,i), 1)
       ENDDO
    ELSE
       DO i = 1,nvar
          im = map_p(i)
          DO j = 1,nvar
             jm = map_p(j)
             hesscr(j,i) = hess(jm,im)
          ENDDO
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE phessco
  ! ==================================================================


END MODULE rlbfgs_utils


SUBROUTINE stack_do_sched (istack, var, lenent)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE mp_interface, ONLY: mp_bcast
  USE parac, ONLY : paral,parai
  USE lscal , ONLY:iopnxt,st_nop
  USE rlbfgs_utils, ONLY : stack_do
  IMPLICIT NONE
  INTEGER                                    :: istack, lenent
  REAL(real_8)                               :: var(lenent)

  INTEGER                                    :: ioper

  ioper = iopnxt(istack)
  CALL mp_bcast (ioper, parai%source, parai%allgrp)
  CALL stack_do (istack, ioper, var, lenent)
  iopnxt(istack) = st_nop
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE stack_do_sched
! ==================================================================
