MODULE rlbfgs_io
  USE cotr,                            ONLY: cotc0
  USE error_handling,                  ONLY: stopgm
  USE lscal,                           ONLY: &
       cuteig, deaccp, demin_p, depred, eigval, eigvec, etotmn, gx, hesscr, &
       ibnd_l, igpnt, info_l, ipnt_l, iprnt, ispnt, iter_l, iwpnt, jump, &
       lallhs, lfevmn, linils, llock, lmap_p, lmicro, lnomap, lprhes, lrstrt, &
       lundop, lvlhes, lwlbfgs, lwprfo, map_p, modcur, mode, modelk, ncore, &
       nfev_l, nfevmn, npt_l, nrem, nres_l, nrestt, nshess, nsmaxp, nstack, &
       nstep_l, nstep_p, nstried_l, nstried_p, ntrstr, ntrust, nvar, oldeig, &
       oldene, oldev, oldg, oldgx, oldhss, omin, ooldg, oprven, ostep_p, &
       osz_p, prvene, rmax_p, rmin_p, scl1_l, scl2_l, step, step_bmb_l, &
       step_ini_l, step_ini_p, step_l, step_max_l, step_max_p, step_min_l, &
       step_min_p, step_p, sz_l, sz_lmn, sz_p, tolenv, trustp, trustr, vmode, &
       wlbfgs, wprfo
  USE parac,                           ONLY: paral
  USE utils,                           ONLY: fskip

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: wr30lbfgs
  PUBLIC :: rd30lbfgs

  PUBLIC :: wr30prfo
  PUBLIC :: rd30prfo

CONTAINS


  ! ==================================================================
  SUBROUTINE wr30lbfgs (nw, ierror, lsavc0)
    ! ==--------------------------------------------------------------==
    ! ==  Writes L-cntl%bfgs optimizer status and history to the RESTART   ==
    ! ==  file NW. Limits and other status variables come first to    ==
    ! ==  enable reading program allocating the memory required to    ==
    ! ==  hold the information.                                       ==
    ! ==  Change this routine consistently with                       ==
    ! ==  lscal.inc and RD30LBFGS                                     ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nw, ierror
    LOGICAL                                  :: lsavc0

    INTEGER                                  :: i

! Variables
! ==--------------------------------------------------------------==
! Input parameters (common blocks /LSCALS/ and /LSCALO/)

    IF (paral%io_parent)&
         WRITE (unit=nw,iostat=ierror)&
         nrem, ntrust, nrestt, ntrstr
    IF (paral%io_parent)&
         WRITE (unit=nw,iostat=ierror)&
         step_bmb_l, step_ini_l, step_max_l, step_min_l, trustr
    ! Shared information (common blocks /LSCALL/, /LSCALI/ and /LSCALR/)
    IF (paral%io_parent)&
         WRITE (unit=nw,iostat=ierror) lfevmn, linils
    IF (paral%io_parent)&
         WRITE (unit=nw,iostat=ierror)&
         ibnd_l, igpnt, info_l, ipnt_l, iprnt, ispnt,&
         ITER_L, IWPNT, JUMP, LWLBFGS, NFEV_L,&
         NFEVMN, NRES_L, NSTEP_L, NSTRIED_L, NPT_L
    IF (paral%io_parent)&
         WRITE (unit=nw,iostat=ierror)&
         deaccp, etotmn, oldene, scl1_l, scl2_l, step_l, sz_l, sz_lmn
    ! L-cntl%bfgs history
    IF (paral%io_parent)&
         WRITE (unit=nw,iostat=ierror) (wlbfgs(i), i=1,lwlbfgs)
    ! Electronic state if requested
    IF (lsavc0) THEN
       IF (paral%io_parent)&
            WRITE (unit=nw,iostat=ierror) nstack
    ELSE
       IF (paral%io_parent)&
            WRITE (unit=nw,iostat=ierror) 0
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE wr30lbfgs
  ! ==================================================================
  SUBROUTINE rd30lbfgs (nr, irecord)
    ! ==--------------------------------------------------------------==
    ! ==  Restart from checkpoint. See WR30LBFGS                      ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nr, irecord

    CHARACTER(*), PARAMETER                  :: procedureN = 'rd30lbfgs'

    INTEGER                                  :: i, ierr

! Variables
! ==--------------------------------------------------------------==
! Input parameters (common blocks /LSCALS/ and /LSCALO/)

    IF (paral%io_parent)&
         READ(nr) nrem, ntrust, nrestt, ntrstr
    IF (paral%io_parent)&
         READ(nr) step_bmb_l, step_ini_l, step_max_l, step_min_l, trustr
    ! Shared information (common blocks /LSCALL/, /LSCALI/ and /LSCALR/)
    IF (paral%io_parent)&
         READ(nr) lfevmn, linils
    IF (paral%io_parent)&
         READ(nr) ibnd_l, igpnt, info_l, ipnt_l, iprnt, ispnt,&
         ITER_L, IWPNT, JUMP, LWLBFGS, NFEV_L,&
         NFEVMN, NRES_L, NSTEP_L, NSTRIED_L, NPT_L
    IF (paral%io_parent)&
         READ(nr) deaccp, etotmn, oldene, scl1_l, scl2_l, step_l,&
         SZ_L, SZ_LMN
    ! L-cntl%bfgs history
    ALLOCATE(wlbfgs(lwlbfgs),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    IF (paral%io_parent)&
         READ(nr) (wlbfgs(i), i=1,lwlbfgs)
    ! Electronic state if stored
    IF (paral%io_parent)&
         READ(nr) i
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rd30lbfgs
  ! ==================================================================
  SUBROUTINE wr30prfo (nw, ierror, irecl)
    ! ==--------------------------------------------------------------==
    ! ==  Writes P-cntl%rfo optimizer status to the RESTART file NW. The   ==
    ! ==  partial Hessian has its own record. See also WR30LBFGS      ==
    ! ==  Change this routine consistently with                       ==
    ! ==  lscal.inc and RD30PRFO                                      ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nw, ierror, irecl

    INTEGER                                  :: i

! Variables
! ==--------------------------------------------------------------==
! NODIM is not yet determined when reading the RESTART file

    IF (paral%io_parent)&
         WRITE (unit=nw,iostat=ierror)&
         cotc0%nodim
    ! Input parameters (common blocks /PRFOPI/ and /PRFOPR/)
    IF (paral%io_parent)&
         WRITE (unit=nw,iostat=ierror)&
         Mode, modelk, ncore, nsmaxp, nvar
    IF (paral%io_parent)&
         WRITE (unit=nw,iostat=ierror)&
         cuteig, demin_p, omin, rmax_p, rmin_p, step_ini_p,&
         STEP_MAX_P, STEP_MIN_P, TOLENV, TRUSTP
    ! Shared information (common blocks /PRFOSL/, /PRFOSI/ and /PRFOSR/)
    IF (paral%io_parent)&
         WRITE (unit=nw,iostat=ierror)&
         lallhs, llock, lmicro, lnomap, lprhes, lrstrt, lundop,&
         LVLHES
    IF (paral%io_parent)&
         WRITE (unit=nw,iostat=ierror)&
         lwprfo, modcur, nshess, nstep_p, nstried_p
    IF (paral%io_parent)&
         WRITE (unit=nw,iostat=ierror)&
         depred, oprven, ostep_p, osz_p, prvene, step_p, sz_p
    ! Work space (common block /PRFOSP/)
    IF (paral%io_parent)&
         WRITE (unit=nw,iostat=ierror) (wprfo(i), i=1,lwprfo)
    IF (lmap_p.GT.0) THEN
       IF (paral%io_parent)&
            WRITE (unit=nw,iostat=ierror) (map_p(i), i=1,lmap_p)
    ELSE
       IF (paral%io_parent)&
            WRITE (unit=nw,iostat=ierror) 0
    ENDIF
    ! Determine if L-cntl%bfgs history needs to be saved as well
    IF (nvar.LT.cotc0%nodim) THEN
       irecl = 1
    ELSE
       irecl = 0
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE wr30prfo
  ! ==================================================================
  SUBROUTINE rd30prfo (nr, irecord, irecl)
    ! ==--------------------------------------------------------------==
    ! ==  Restart from checkpoint. See WR30PRFO                       ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nr, irecord, irecl

    CHARACTER(*), PARAMETER                  :: procedureN = 'rd30prfo'

    INTEGER                                  :: i, ierr, j, ndimsv, nvarnp
    LOGICAL                                  :: lvlhsloc

    nvarnp = nvar
    ! NODIM is not yet determined at this point
    IF (paral%io_parent)&
         READ(nr) ndimsv
    ! Input parameters (common blocks /PRFOPI/ and /PRFOPR/)
    IF (paral%io_parent)&
         READ(nr) Mode, modelk, ncore, nsmaxp, nvar
    IF (paral%io_parent)&
         READ(nr) cuteig, demin_p, omin, rmax_p, rmin_p, step_ini_p,&
         STEP_MAX_P, STEP_MIN_P, TOLENV, TRUSTP
    ! Check if RESTART file is consistent with current parameters
    IF (nvarnp.NE.nvar) THEN
       IF (paral%io_parent)&
            WRITE (6,'(1X,A,A,I3)') 'RD30PRFO| WARNING! ',&
            'NVAR ON CHECKPOINT: ', NVAR
       IF (paral%io_parent)&
            WRITE (6,'(1X,A,I3)') 'NVAR IN INPUT: ', nvarnp
       IF (paral%io_parent)&
            WRITE (6,'(1X,A)') 'NOT USING RESTART FILE FOR OPTIMIZATION'
       nvar = nvarnp
       CALL fskip (nr, irecord-2)
       ! Read L-cntl%bfgs status and history if microiterative optimization
       irecl = 0
       RETURN
    ELSE IF (nvar.LT.ndimsv) THEN
       irecl = 1
    ELSE
       irecl = 0
    ENDIF
    ! Shared information (common blocks /PRFOSL/, /PRFOSI/ and /PRFOSR/)
    IF (paral%io_parent)&
         READ(nr) lallhs, llock, lmicro, lnomap, lprhes, lrstrt, lundop,&
         LVLHSLOC
    IF (paral%io_parent)&
         READ(nr) lwprfo, modcur, nshess, nstep_p, nstried_p
    IF (paral%io_parent)&
         READ(nr) depred, oprven, ostep_p, osz_p, prvene, step_p, sz_p

    ! Work space (common block /PRFOSP/)
    ! CALL MEMORY (IP_WPRFO, LWPRFO, 'WPRFO')
    ! TODO decide where to free these vars
    ALLOCATE(hesscr(nvar, nvar),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
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

    ! all these arrays are accessible via lscal.mod 
    IF (paral%io_parent) THEN
       READ (nr) ((hesscr(i,j), i=1,nvar), j=1,nvar)
       READ (nr) ((oldhss(i,j), i=1,nvar), j=1,nvar)
       READ (nr) (vmode(i), i=1,nvar)
       READ (nr) (gx(i), i=1,nvar)
       READ (nr) (oldgx(i), i=1,nvar)
       READ (nr) (oldg(i), i=1,nvar)
       READ (nr) (ooldg(i), i=1,nvar)
       READ (nr) (eigval(i), i=1,nvar)
       READ (nr) (oldeig(i), i=1,nvar)
       READ (nr) ((eigvec(i,j), i=1,nvar), j=1,nvar)
       READ (nr) ((oldev(i,j), i=1,nvar), j=1,nvar)
       READ (nr) (step(i), i=1,nvar)
    ENDIF
    IF (lmap_p.GT.0) THEN
       ! CALL MEMORY (IP_MAP_P, LMAP_P, 'MAP_P')
       ALLOCATE(map_p(lmap_p),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       IF (paral%io_parent)&
            READ(nr) (map_p(i), i=1,lmap_p)
    ELSE
       IF (paral%io_parent)&
            READ(nr) i
    ENDIF
    ! Warn if no Hessian was read
    IF (lvlhsloc.AND..NOT.lvlhes) THEN
       IF (paral%io_parent)&
            WRITE (6,'(1X,A,A)') 'RD30PRFO| WARNING! NO HESSIAN READ ',&
            'FROM RESTART FILE BUT P-RFO STATUS REQUIRES ONE'
       IF (paral%io_parent)&
            WRITE (6,'(1X,A)') 'WILL RECALCULATE HESSIAN WHEN REQUIRED'
    ENDIF
    ! Set pointers within work space and return codes to RGMOPT
    lrstrt = .TRUE.
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rd30prfo
  ! ==================================================================


END MODULE rlbfgs_io

