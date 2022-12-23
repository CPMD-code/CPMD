MODULE meta_exl_mult_utils
  USE cnst,                            ONLY: au_kcm,&
                                             kboltz,&
                                             pi,&
                                             scmass
  USE cnst_dyn,                        ONLY: &
       cscl_fac, cscl_val, cv_associated, cv_dtemp, cv_dyn, cv_ist, cv_mass, &
       cv_path, cv_temp, cv_temp0, cv_vel, det_colvar, ekincv, fhills, hhm, &
       hlhm, hllh_val, hllw_val, hthm, hvm0, hwm, iangcv, ibound, icv_spin, &
       imeta, inter_hill, kharm, lchekharm, lcvtc, lkfix, lmeta, ncolvar, &
       ncvsys, nsubsys, rcc0, rmeta, toll_avcv, tycvar, vbound
  USE cotr,                            ONLY: cotc0
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_app,&
                                             fo_verb
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE machine,                         ONLY: m_flush
  USE meta_colvar_inp_utils,           ONLY: colvarofr,&
                                             colvarpr,&
                                             cv_forces3
  USE meta_ex_mul_util_utils,          ONLY: calc_pos_mdyn,&
                                             cv_diff_md,&
                                             cv_exl_outm,&
                                             cv_exl_outm2,&
                                             cv_read_outm,&
                                             rmtdresm,&
                                             wmtdresm
  USE meta_exlagr_utils,               ONLY: rinvelcv,&
                                             rscvelcv,&
                                             scf_tune
  USE meta_hpot_utils,                 ONLY: hills,&
                                             hills_lor,&
                                             hills_ratio,&
                                             hills_sals,&
                                             hills_sals_shift,&
                                             setwall,&
                                             setwall_old
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sync
  USE parac,                           ONLY: parai,&
                                             paral
  USE puttau_utils,                    ONLY: gettau
  USE readsr_utils,                    ONLY: xstring
  USE ropt,                            ONLY: infi,&
                                             iteropt
  USE soft,                            ONLY: soft_com
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             cntr
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE tpar,                            ONLY: dt_ions,&
                                             dtb2mi
  USE zeroing_utils,                   ONLY: zeroing
#include "sizeof.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: meta_ext_mul

CONTAINS

  ! ==================================================================
  SUBROUTINE meta_ext_mul(taup,velp,tscr,lquench,lmetares,&
       resetcv,ekinc,ekinp)
    ! ==--------------------------------------------------------------==
    ! ==  Extended Lagrangean: Collective Variables                   ==
    ! ==      L = T_elec + T_ions - E_ks({r_i},{R_I})
    ! ==          + T_CV                   : kinetic energy of CV     ==
    ! ==          -SUM_icv k_icv(S({R_I},icv)-S0(icv))^2              ==
    ! : harmonic potential       ==
    ! that couples the ionic   ==
    ! DoF to the CV            ==
    ! ==          -V_hills(time,S0(icv))   : time dep., given by      == 
    ! progressive accumulation == 
    ! of hills, whose shape    == 
    ! is NCV-dim. Gaussian like== 
    ! ==          -V_walls(S0(icv))        : boundary conditions      ==
    ! ==                                                              ==
    ! ==  This routine calculates the contributions of the harmonic   ==
    ! ==  part of potential to be added to  the forces on ions        ==
    ! ==                                                              ==
    ! ==  It calculates the dynamics of the collective variables      ==
    ! ==  that is determined by the kinetic energy                    ==
    ! ==  and the added potential terms, and depends on the choice of ==
    ! ==  masses and harmonic constants for the CV                    ==
    ! ==                                                              ==
    ! ==  The path followed by the CV in the CV space is stored in    ==
    ! ==  several outfiles, which contain CV position and also        == 
    ! ==  informations about the hills filling the PES wells and      ==
    ! ==  other deteils regarding the dynamics
    ! ==--------------------------------------------------------------==

    ! NCOLVAR   : # of collective variables
    ! LHILLS    : true when hills are used in the dynamics of the constraints
    ! I_META    : # of steps in the meta-dynamics
    ! I_META_MAX: max # of steps in the meta-dynamics
    ! I_META_RES: if restart from previous meta dynamics, 
    ! # of steps already done
    ! I_CVST    : # of steps since last HILL has been added
    ! CV_IST    : actual values of collective variables, dim NCOLVAR
    ! CV_DYN    : value of the dynamically evolving collective variables 
    ! CV_VEL    : velocity of the collective variables
    ! CVD_IST   : actual values of diffusion coefficient of CV
    ! HLLH      : hills altitude 
    ! HLLW      : hills amplitude 
    ! CSCL_FAC  : scale factors required for an homogeneous dynamic in the
    ! MCNSTR directions in the space of phases
    ! they can be read from input  and they  can be tuned on the
    ! amplitudes of the variables fluctuations
    ! IBOUND    : 1 if a wall limiting the constraint in one direction 
    ! has to be set, otherwise 0 (integer, dimension MCNSTR)
    ! VBOUND    : parameters for the definition of constraining potential
    ! (real, dimension 4*NCOLVAR)
    ! if WALL+ VBOUND(2) is the barrier for CV_DYN > VBOUND(1)
    ! if WALL- VBOUND(4) is the barrier for CV_DYN < VBOUND(3)
    ! FI_HARM   : force contribution on ions due to the harm. pot. 
    ! dim (NODIM)
    ! F_HARM    : force contribution due to the harm. pot.
    ! and acting on the dynamical CV
    ! dim (NCOLVAR)
    ! F_HILL    : force contribution due to the time dep. hill potentiel
    ! and acting on the dynamical CV
    ! dim (NCOLVAR)
    ! F_WALL    : force contribution due to the boundary potential 
    ! and acting on the dynamical CV
    ! dim (NCOLVAR)
    ! META_RESTART: restart for a previous meta-dynamics
    ! acculators and the constraints path have to be read 
    ! from the right files
    ! the output is appended 
    ! CV_PATH   : CV history (dim NTOT_ITER*NCOLVAR)
    ! CSCL_VAL  : CV scaling factors history (dim NTOT_ITER*NCOLVAR)
    ! HLLW_VAL  : hills amplitude (dim NTOT_ITER*NSUBSYS)
    ! HLLH_VAL  : hills altitude  (dim NTOT_ITER*NSUBSYS)
    ! F_AVER    : accumulator for the calculation of the average force
    ! acting on the CV, as a contribution of the 
    ! underlaying potential 
    ! IF_AVER   : counter for the calculation of the average
    ! IW_CV     : counter for writing the CV values on the cvmdck_mtd file
    ! in order to monitor the behavior along the cntl%md trajectory
    ! DISP_SYS  : norm of the displacment in CV space (dim NSUBSYS)
    ! CV_MASS   : masses of collective variables, if not assigned in 
    ! the input file, they are automathically calculated
    ! dim (NCOLVAR)
    ! KHARM     : constant of the harmonic part of the potential, if not
    ! assigned in input file, default values are chosen
    ! dim (NCOLVAR)
    ! ==--------------------------------------------------------------==

    REAL(real_8)                             :: taup(:,:,:), velp(:,:,:), &
                                                tscr(:,:,:)
    LOGICAL                                  :: lquench, lmetares, resetcv
    REAL(real_8)                             :: ekinc, ekinp

    CHARACTER(*), PARAMETER                  :: procedureN = 'meta_ext_mul'

    CHARACTER(len=10)                        :: chnum, file9 = 'cvmdck_mtd'
    CHARACTER(len=100)                       :: lineform, outcheck(10)
    CHARACTER(len=20)                        :: flag
    INTEGER                                  :: i, ia, icv, idof, ie, ierr, &
                                                ifile, is, isub, jj, ntot_iter
    INTEGER, ALLOCATABLE                     :: isys(:)
    INTEGER, SAVE                            :: i_cvst = 0, i_meta = 0, &
                                                i_temp = 0, if_aver = 0, &
                                                ifirst = 0, ifirstp = 0, &
                                                iw_cv = 0
    LOGICAL                                  :: ferror, FINAL, hill_add, lskip
    REAL(real_8)                             :: dif_cv, fact, fact2, hh, &
                                                hh_max, hh_min, hhh, hhht, &
                                                hllh0(10), maxf2, tollm
    REAL(real_8), ALLOCATABLE                :: cv_disp(:), disp_sys(:), &
                                                fi_harm(:), gpot_sys(:), &
                                                hh_test(:)
    REAL(real_8), ALLOCATABLE, SAVE          :: cv_scl(:), cvd_ist(:), &
                                                cvd_scl(:), f_aver(:), &
                                                f_harm(:), f_hill(:), &
                                                f_wall(:), hc_last(:)
    REAL(real_8), SAVE                       :: temp1 = 0.0_real_8, &
                                                temp2 = 0.0_real_8, tempp

! ==--------------------------------------------------------------==
! ==--------------------------------------------------------------==
! Set time for the routine

    CALL tiset(' META_EXTL',isub)
    FINAL = .FALSE.

    IF (nsubsys .GT. 10) CALL stopgm(' META_EXTL',&
         'Too Many Subsystems for Multi-Metadynamics, max is 10',& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    IF (ifirstp .EQ. 0) THEN
       i_meta  = imeta%i_meta_res+1
       i_temp  = 0
       i_cvst  = 0
       if_aver = 0

       ! Allocate Memory

       ALLOCATE(cv_dyn(ncolvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       cv_associated = .TRUE.
       ifirstp = 1
    ENDIF
    IF (.NOT.paral%io_parent) GOTO 9999

    ! ==--------------------------------------------------------------==
    ! Path of output files
    DO is = 1,nsubsys
       IF (paral%io_parent)&
            WRITE(flag,'(I5)') is
       CALL xstring(flag,ia,ie)
       outcheck(is)   = file9//'_S'//flag(ia:ie)
    ENDDO

    ! ==--------------------------------------------------------------==
    ! TODO align for BG
    ALLOCATE(fi_harm(cotc0%nodim),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL zeroing(fi_harm)
    ALLOCATE(cv_disp(ncolvar),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL zeroing(cv_disp)
    ALLOCATE(hh_test(ncolvar),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL zeroing(hh_test)
    ALLOCATE(cv_scl(ncolvar),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL zeroing(cv_scl)
    ALLOCATE(cvd_scl(ncolvar),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL zeroing(cvd_scl)

    ALLOCATE(disp_sys(nsubsys),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL zeroing(disp_sys)
    ALLOCATE(gpot_sys(nsubsys),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL zeroing(gpot_sys)

    ALLOCATE(isys(ncolvar),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL zeroing(isys)

    ntot_iter = imeta%i_meta_max+imeta%i_meta_res
    tollm     = toll_avcv


    IF (ifirst .EQ. 0) THEN
       ! Allocate Memory

       CALL dcopy(ncolvar,cv_ist,1,cv_dyn,1)
       ALLOCATE(cv_vel(ncolvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(cv_vel)!,ncolvar)
       ALLOCATE(cvd_ist(ncolvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(cvd_ist)!,ncolvar)
       ALLOCATE(hc_last(ncolvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(hc_last)!,ncolvar)
       ALLOCATE(f_aver(ncolvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(f_aver)!,ncolvar)
       ALLOCATE(f_harm(ncolvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(f_harm)!,ncolvar)
       ALLOCATE(f_hill(ncolvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(f_hill)!,ncolvar)
       ALLOCATE(f_wall(ncolvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(f_wall)!,ncolvar)

       ALLOCATE(cv_path(imeta%i_meta_max+imeta%i_meta_res,ncolvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(cv_path)!,ncolvar*ntot_iter)
       ALLOCATE(cscl_val(imeta%i_meta_res+imeta%i_meta_max,ncolvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(cscl_val)!,ncolvar*ntot_iter)
       ALLOCATE(hllw_val(imeta%i_meta_res+imeta%i_meta_max,nsubsys),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(hllw_val)!,ntot_iter*nsubsys)
       ALLOCATE(hllh_val(imeta%i_meta_res+imeta%i_meta_max,nsubsys),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(hllh_val)!,ntot_iter*nsubsys)


       jj = 0
       DO is=1,nsubsys
          hllh0(is) = hhm(is)
          DO icv = 1,ncvsys(is)
             isys(jj+icv) = is
          ENDDO
          jj = jj + ncvsys(is)
       ENDDO

       ! If META_RESTART read data from cntl%proper input
       IF (lmeta%meta_restart  .AND. .NOT. lmeta%tcvanalysis) THEN
          IF (lmeta%tresfile) THEN
             CALL rmtdresm(50,cv_disp,.TRUE.,ntot_iter)
          ELSE
             CALL cv_read_outm(cv_disp)
          ENDIF! TRESFILE
          !$omp parallel do private(ICV)
          DO icv = 1,ncolvar
             hc_last(icv) = cv_path(imeta%i_meta_res,icv)+cv_disp(icv)
          ENDDO
          CALL dcopy(ncolvar,hc_last,1,cv_dyn,1)
       ENDIF ! META_RESTART

       IF (.NOT. lmeta%tcvanalysis) THEN
          ! Masses
          IF (cv_temp0 .NE. 0.0_real_8) THEN
             fact = cv_temp0*kboltz*10._real_8
          ELSE
             fact = cntr%tempw   *kboltz*10._real_8
          ENDIF
          fact2 = fact * cntr%delt_ions*REAL(inter_hill,kind=real_8)*&
               cntr%delt_ions*REAL(inter_hill,kind=real_8)
          DO icv = 1,ncolvar
             IF (cv_mass(icv) .EQ. 0.0_real_8) THEN
                cv_mass(icv) = fact2/&
                     ((cscl_fac(1,icv)*hwm(isys(icv)))*&
                     (cscl_fac(1,icv)*hwm(isys(icv))))
                IF (paral%io_parent)&
                     WRITE(6,'(A,I5,A,f15.4,A)') 'Default Mass for CV #',icv,&
                     ' is assigned  = ',&
                     cv_mass(icv)/scmass,' a.m.u.(a.u.)^2/(u.CV)^2'
                cv_mass(icv) = cv_mass(icv)*scmass
             ELSE
                IF (paral%io_parent)&
                     WRITE(6,'(A,I5,A,f15.4,A)') 'Mass for CV #',icv,&
                     ' is assigned  = ',&
                     cv_mass(icv) ,' a.m.u.(a.u.)^2/(u.CV)^2'
                cv_mass(icv) =  scmass*cv_mass(icv)
             ENDIF

             ! Harmonic Potential
             IF (kharm(icv) .EQ. 0.0_real_8) THEN
                kharm(icv) = fact/&
                     ((cscl_fac(1,icv)*rmeta%hllw)*(cscl_fac(1,icv)*rmeta%hllw))
                IF (paral%io_parent)&
                     WRITE(6,'(A,I5,A,f10.4,A)')&
                     'Default harmonic constant for CV #',icv,&
                     ' is assigned ',kharm(icv),' a.u./(u.CV)^2'
             ENDIF
             IF (.NOT. lkfix) lchekharm = .TRUE.
          ENDDO

          ! Initialization of CV Velocity
          IF (i_meta .EQ. 1) THEN
             CALL rinvelcv(cv_vel,cv_mass,cscl_fac,ncolvar,&
                  cv_temp0,ekincv)
             IF (paral%io_parent)&
                  WRITE(6,'(2x,A,f10.4)')&
                  'CV velocity have been initialized at T = ',cv_temp0
             IF (paral%io_parent)&
                  WRITE(6,'(2x,A,f10.4)')&
                  'the initial CV kinetic energy is = ',ekincv
          ENDIF

          IF (lcvtc)THEN
             temp1=cv_temp+cv_dtemp
             temp2=cv_temp-cv_dtemp
          ENDIF
          tempp = cv_temp0
       ELSE
          CALL zeroing(cv_vel)!,ncolvar)
       ENDIF
       ! ==--------------------------------------------------------------==

       IF (lmeta%tcvmonitor) THEN
          DO is = 1, nsubsys
             ifile = 444+is
             IF (paral%io_parent)&
                  CALL fileopen(ifile,outcheck(is),fo_app+fo_verb,ferror)
          ENDDO
       ENDIF
    ENDIF

    ! ==--------------------------------------------------------------==
    ! Calculate the Collective Variables and their Derivatives
    CALL colvarofr(taup,tscr)

    IF (resetcv .AND. lmeta%tlocalizespin) THEN
       DO icv = 1,ncolvar
          IF (icv_spin(icv) .NE. 0) THEN
             cv_dyn(icv) = cv_ist(icv)
          ENDIF
       ENDDO
       resetcv = .FALSE.
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Initialization

    IF (ifirst.EQ.0) THEN
       IF (i_meta.EQ.1)THEN
          DO icv = 1,ncolvar
             hc_last(icv) = cv_ist(icv)
          ENDDO
       ENDIF
       ifirst=1
       IF (.NOT. lmeta%meta_restart .AND. lmeta%tlocalizespin) THEN
          CALL dcopy(ncolvar,cv_ist,1,cv_dyn,1)
       ENDIF
       GOTO 200
    ENDIF

    hill_add = .FALSE.
    IF (lmeta%ttrajovrr) cnti%ntraj=999999999

    i_cvst = i_cvst + 1
    iw_cv=iw_cv+1


    ! ==--------------------------------------------------------------==
    ! Check  CV_DYN wrt the center of the last Hill (HC_LAST)
    ! to be verified conditions: 
    ! minimum time interval 
    ! minimum displacement in CV space
    ! If the conditions are sadisfied HILL_ADD = true 
    ! and a new hill will be added
    CALL calc_pos_mdyn(tollm,i_cvst,i_temp,i_meta,hc_last,&
         hill_add)
    IF (lmetares) THEN
       FINAL = .TRUE.
       GOTO 9999
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Tune the Harmonic Potential according to the PES

    IF (lchekharm .AND. hill_add) THEN
       CALL cv_forces3(taup,tscr,velp,hh_test)
       DO icv = 1,ncolvar
          dif_cv = cv_ist(icv)-cv_dyn(icv)
          IF (iangcv(icv) .EQ. 1 .AND. dif_cv .GT. pi) THEN
             dif_cv = dif_cv -2.0_real_8*pi
          ELSEIF (iangcv(icv) .EQ. 1 .AND. dif_cv .LT. -pi) THEN
             dif_cv = dif_cv +2.0_real_8*pi
          ENDIF
          IF (hh_test(icv).LE.1.e+1_real_8 .AND. hh_test(icv).GT.0.2_real_8) THEN
             fact = ABS(dif_cv)/cscl_fac(1,icv)
             IF (fact.LT. rmeta%hllw) THEN
                kharm(icv) = ABS(hh_test(icv)/2._real_8/hwm(isys(icv)))
             ELSE
                kharm(icv) = ABS(hh_test(icv)/2._real_8/dif_cv)
             ENDIF
          ELSE
          ENDIF
       ENDDO
    ENDIF

    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    ! Calculation of the forces on the collective variables
    ! This is used in order to determine the optimal hills parameters 
    ! In F_aver the forces on CV are accumulated in order to calculate
    ! their average value when the meta step is concluded
    ! MULTY
    IF (lmeta%ltune_hh) THEN
       DO icv = 1,ncolvar
          f_aver(icv) = f_aver(icv)&
                                ! &        + 2.0_real_8*KHARM(ICV)*DIF_CV
               + f_harm(icv)+f_hill(icv)/cscl_fac(1,icv)
       ENDDO
       if_aver = if_aver + 1
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Monitoring of the CV along the cntl%md trajectory. Even if no hill is added
    ! at this point, if IW_CV is a multiple of WCV_FREQ, 
    ! the values of the CV_IST and of the CV_DYN
    ! are written in the cvmdck_mtd file
    ! in the first column the IW_CV counter is reported, that indicates the
    ! number of cntl%md steps of the actual run, in the second column the global 
    ! metastep is reported, afterwards the values of the two CV sets
    IF (lmeta%tcvmonitor .AND. MOD(iw_cv,imeta%wcv_freq) .EQ. 0) THEN
       DO is = 1,nsubsys
          ifile = 444 + is
          IF (paral%io_parent)&
               WRITE(chnum,'(I5)') ncvsys(is)
          CALL xstring(chnum,ia,ie)
          lineform =&
               '(1X,2I9,'//chnum(ia:ie)//&
               'E14.6,'//chnum(ia:ie)//'E14.6,E14.6)'
          IF (paral%io_parent)&
               WRITE(ifile,lineform) iw_cv,iteropt%nfi,&
               (cv_ist(icv),icv=1,ncvsys(is)),&
               (cv_dyn(icv),icv=1,ncvsys(is)),tempp
          CALL m_flush(ifile)
       ENDDO
    ENDIF

    ! ==--------------------------------------------------------------==
    ! New Step of Meta-Dynamics in the Space of CV
    IF (hill_add) THEN
       IF (MOD(i_meta,imeta%tr_freq) .EQ. 0 .OR. i_meta .EQ. 1 ) THEN
          IF (lmeta%ttrajovrr) cnti%ntraj=1
       ENDIF

       IF (.NOT. lmeta%tcvanalysis) THEN
          ! ==--------------------------------------------------------------==
          ! Squared Norm of the Displacement in the CV space


          ! MULTY
          jj = 0
          DO is = 1,nsubsys
             disp_sys(is) = 0.0_real_8
             DO icv = 1,ncvsys(is)
                dif_cv = cv_dyn(jj+icv)-hc_last(jj+icv)
                IF (iangcv(jj+icv) .EQ. 1 .AND. dif_cv .GT. pi) THEN
                   dif_cv = dif_cv - 2.0_real_8*pi
                ELSEIF (iangcv(jj+icv) .EQ. 1 .AND. dif_cv .LT. -pi) THEN
                   dif_cv = dif_cv + 2.0_real_8*pi
                ENDIF
                fact = dif_cv/cscl_fac(1,jj+icv)
                disp_sys(is) = disp_sys(is) + fact*fact
             ENDDO
             disp_sys(is) = SQRT(disp_sys(is))
             jj = jj + ncvsys(is)
          ENDDO

          ! ==--------------------------------------------------------------==
          ! Tuning Hills parameters 

          IF (lmeta%ltune_hh) THEN
             ! MULTY
             jj = 0
             DO is = 1,nsubsys
                hh_min = 1000._real_8
                hh_max = 0.0_real_8
                lskip = .FALSE.

                hhht = 0.0_real_8
                DO icv = 1,ncvsys(is)
                   f_aver(jj+icv)  = f_aver(jj+icv) / REAL(if_aver,kind=real_8)
                   f_aver(jj+icv)  = f_aver(jj+icv)
                   dif_cv = cv_dyn(jj+icv)-hc_last(jj+icv)
                   IF (iangcv(jj+icv) .EQ. 1 .AND. dif_cv .GT. pi) THEN
                      dif_cv = dif_cv - 2.0_real_8*pi
                   ELSEIF (iangcv(jj+icv) .EQ. 1 .AND. dif_cv .LT. -pi) THEN
                      dif_cv = dif_cv + 2.0_real_8*pi
                   ENDIF
                   fact = dif_cv/(cscl_fac(1,jj+icv)*cscl_fac(1,jj+icv))*&
                        (1.0_real_8/(hwm(is)*hwm(is))&
                        +1.0_real_8/(disp_sys(is)*disp_sys(is)))
                   hhh  = 0.5_real_8*ABS(f_aver(jj+icv)*EXP(0.5_real_8)/fact)
                   hhht = hhht + hhh

                   IF (paral%io_parent)&
                        WRITE(6,'(2x,A,I4,3F14.6)') 'FCV ',jj+icv,&
                        f_aver(jj+icv),hhh,fact

                   hh_min       = MIN(hh_min,hhh)
                   hh_max       = MAX(hh_max,hhh)
                ENDDO
                hh       = hhht / REAL(ncvsys(is),kind=real_8)
                IF (hh_min .GT. hthm(is) .OR. hh_max .LT. hlhm(is)) THEN
                   lskip = .TRUE.
                   IF (paral%io_parent)&
                        WRITE(6,'(A,I4,A,f10.5,A,f10.5,A)') 'Set CV', is,&
                        '  HM min > ', hthm(is)*au_kcm,&
                        ' Kcal, or HH max < ',hlhm(is)*au_kcm,' Kcal'
                ENDIF
                IF (.NOT. lskip .AND. hh .LT. hlhm(is)) THEN
                   hhm(is) = (hlhm(is)+hh_max)*0.5_real_8
                ELSEIF (.NOT. lskip .AND. hh .GT. rmeta%htop) THEN
                   hhm(is) = (hh_min+hthm(is))*0.5_real_8
                ELSEIF (.NOT. lskip ) THEN
                   hhm(is) = hh
                ELSE
                   hhm(is) = hhm(is)*0.2_real_8+hllh0(is)*0.8_real_8
                ENDIF
                IF (paral%io_parent)&
                     WRITE(6,'(/A,I4,3(2x,A,f10.5)/)') 'Set CV', is,&
                     'HH min = ', hh_min,'HH max = ', hh_max,&
                     ' HH used ',hhm(is)
                jj = jj + ncvsys(is)
                IF (lmeta%thillvol) THEN
                   hwm(is)=(hvm0(is)/hhm(is))**(1._real_8/REAL(ncvsys(is),kind=real_8))
                ENDIF
             ENDDO

             if_aver = 0
             CALL zeroing(f_aver)!,ncolvar)
          ENDIF

          ! ==--------------------------------------------------------------==
          ! Scaling Factors Tuning
          IF (lmeta%ltune_cscl) THEN

             CALL scf_tune(i_meta,hc_last)

          ENDIF
          ! ==--------------------------------------------------------------==
          ! Diffusion calculation
          CALL cv_diff_md(hc_last,cvd_ist,cvd_scl,cv_scl,i_meta,&
               ntot_iter)

          ! ==--------------------------------------------------------------==

          ! Open output file
          CALL cv_exl_outm(i_meta,hc_last,f_harm,f_hill,f_wall,&
               cvd_ist,disp_sys,ekinc,ekinp)

          ! ==--------------------------------------------------------------==

          DO icv = 1,ncolvar
             hc_last(icv)   = cv_dyn(icv)
          ENDDO
       ELSE
          ! ==--------------------------------------------------------------==
          ! Open output file
          CALL cv_exl_outm2(i_meta,ekinc,ekinp)

          ! ==--------------------------------------------------------------==

       ENDIF

       i_meta = i_meta + 1
       i_temp = 0
       i_cvst = 0
    ENDIF        ! HILL_ADD

    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    ! Calculation of the dynamics of the collective variables
    ! A velocity Verlet Scheme is used

    ! ==--------------------------------------------------------------==
    ! Update Velocity : first half step
    ! using forces calculated in the previous step
    DO icv = 1,ncolvar
       fact = dt_ions/(2._real_8*cv_mass(icv))
       cv_vel(icv) = cv_vel(icv)&
            + fact*(f_harm(icv)+f_hill(icv)/cscl_fac(1,icv)+f_wall(icv))
    ENDDO

    ! ==--------------------------------------------------------------==
    ! Update Position : second order dynamics

    DO icv = 1,ncolvar
       fact = dt_ions
       IF (iangcv(icv) .EQ. 1 )  THEN
          IF (cv_dyn(icv).LT.0._real_8) cv_dyn(icv)= cv_dyn(icv)+2._real_8*pi
          IF (cv_dyn(icv).GT.2._real_8*pi) cv_dyn(icv)= cv_dyn(icv)-2._real_8*pi
       ELSEIF (iangcv(icv) .EQ. 2 )  THEN
          IF (cv_dyn(icv).LT.0._real_8) cv_dyn(icv)= -cv_dyn(icv)
          IF (cv_dyn(icv).GT. pi) cv_dyn(icv)= pi - (cv_dyn(icv) - pi)
       ENDIF
       cv_dyn(icv) = cv_dyn(icv) + fact*cv_vel(icv)
    ENDDO

200 CONTINUE

    ! ==--------------------------------------------------------------==
    ! Calculation of the contributions to the forces on ions 
    ! due to the coupling with the collective variables
    ! V_ARM = SUM_icv (K_icv(CV_IST_icv - CV_DYN_icv))^2
    ! 
    CALL zeroing(fi_harm)!,cotc0%nodim)
    DO icv = 1,ncolvar
       dif_cv = cv_ist(icv)-cv_dyn(icv)
       IF (iangcv(icv) .EQ. 1 .AND. dif_cv .GT. pi) THEN
          dif_cv = dif_cv -2.0_real_8*pi
       ELSEIF (iangcv(icv) .EQ. 1 .AND. dif_cv .LT. -pi) THEN
          dif_cv = dif_cv +2.0_real_8*pi
       ENDIF
       fact = kharm(icv)*2.0_real_8*dif_cv
       DO idof = 1,cotc0%nodim
          fi_harm(idof) = fi_harm(idof) -&
               fact * det_colvar(idof,icv)
       ENDDO
    ENDDO

    ! ==--------------------------------------------------------------==
    ! Update forces on ions

    CALL zeroing(tscr)!,3*maxsys%nax*maxsys%nsx)
    CALL gettau(tscr,fi_harm)
#ifdef __SR11000
    !poption parallel, tlocal(IS,IA,FACT)
    !voption indep(FHILLS,TSCR)
#endif
    DO is=1,ions1%nsp
       fact= dt_ions*dtb2mi(is)
       DO ia=1,ions0%na(is)
          fhills(1,ia,is) = tscr(1,ia,is)
          fhills(2,ia,is) = tscr(2,ia,is)
          fhills(3,ia,is) = tscr(3,ia,is)
       ENDDO
    ENDDO

    ! ==--------------------------------------------------------------==
    ! Calculate contribution of the Harmonic Potential (recall)

    CALL zeroing(f_harm)!,ncolvar)
    DO icv = 1,ncolvar
       dif_cv = cv_ist(icv)-cv_dyn(icv)
       IF (iangcv(icv) .EQ. 1 .AND. dif_cv .GT. pi) THEN
          dif_cv = dif_cv - 2.0_real_8*pi
       ELSEIF (iangcv(icv) .EQ. 1 .AND. dif_cv .LT. -pi) THEN
          dif_cv = dif_cv + 2.0_real_8*pi
       ENDIF
       f_harm(icv) =  kharm(icv)*2.0_real_8*dif_cv
    ENDDO

    ! Calculate hills contribution to forces
    IF (lmeta%lhills) THEN
       ! MULTY
       jj = 0
       DO is = 1,nsubsys
          CALL zeroing(f_hill(jj+1:jj+ncvsys(is)))!,ncvsys(is))
          IF (lmeta%hlore) THEN
             CALL hills_lor(cv_dyn(jj+1),ncvsys(is),i_meta,&
                  ntot_iter,f_hill(jj+1),&
                  is,jj)
          ELSEIF (lmeta%hratio) THEN
             CALL hills_ratio(cv_dyn(jj+1),ncvsys(is),i_meta,&
                  ntot_iter,f_hill(jj+1),hc_last(jj+1),&
                  i_cvst,is,jj)
          ELSEIF (lmeta%hshift) THEN
             CALL hills_sals_shift(cv_dyn(jj+1),ncvsys(is),i_meta,&
                  ntot_iter,f_hill(jj+1),hc_last(jj+1),&
                  i_cvst,is,jj)
          ELSE IF (lmeta%sphere) THEN
             CALL hills(cv_dyn,ncvsys(is),i_meta,ntot_iter,f_hill(jj+1),&
                  is,jj)
          ELSE
             CALL hills_sals(cv_dyn(jj+1),ncvsys(is),i_meta,&
                  ntot_iter,f_hill(jj+1),hc_last(jj+1),&
                  i_cvst,is,jj)
          ENDIF
          gpot_sys(is) = rmeta%gausspot
          jj = jj + ncvsys(is)
       ENDDO
       rmeta%gausspot = 0._real_8
       DO is = 1,nsubsys
          rmeta%gausspot = rmeta%gausspot + gpot_sys(is)
       ENDDO
    ENDIF

    ! Calculate walls contribution to forces
    CALL zeroing(f_wall)!,ncolvar)
    maxf2 = 0.0_real_8
    DO icv = 1,ncolvar
       IF (ibound(icv) .EQ. 1) THEN
          CALL setwall(cv_dyn(icv),vbound(1,icv),f_wall(icv))
       ELSEIF (ibound(icv) .EQ. 2) THEN
          CALL setwall_old(cv_dyn(icv),tycvar(icv),vbound(1,icv),&
               f_wall(icv))
       ENDIF
       IF (f_wall(icv) .NE. 0.0_real_8) THEN
          IF (paral%io_parent)&
               WRITE(6,'(72A)') ('*',i = 1,72)
          IF (paral%io_parent)&
               WRITE(6,'(10x,A,I5)')&
               'WARNING: Boundaries active for CV ', icv
          IF (paral%io_parent)&
               WRITE(6,*)  icv, f_wall(icv)
          IF (paral%io_parent)&
               WRITE(6,'(72A)')  ('*',i = 1,72)
       ENDIF
    ENDDO

    ! ==--------------------------------------------------------------==
    ! Update Velocity : second half step 
    DO icv = 1,ncolvar
       fact = dt_ions/(2._real_8*cv_mass(icv))
       cv_vel(icv) = cv_vel(icv)&
            + fact*(f_harm(icv)+f_hill(icv)/cscl_fac(1,icv)+f_wall(icv))
    ENDDO

    ! ==--------------------------------------------------------------==
    ! Calculation of   Temperature of CV: CV_TEMP
    ! If CV_TEMP_CONTROL, check the condition for rescaling T

    CALL rscvelcv(temp1,temp2,tempp)

    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    ! IF MAX_STEP reached or EXIT file exist

    IF (infi.EQ.cnti%nomore .OR. i_meta-1.EQ.ntot_iter) THEN
       soft_com%exsoft = .TRUE.

       CALL  colvarpr
       IF (lmeta%tcvmonitor) THEN
          DO is = 1,nsubsys
             ifile = 444+is
             IF (paral%io_parent)&
                  CALL fileclose(ifile)
          ENDDO
       ENDIF

    ENDIF
9999 CONTINUE
    CALL mp_sync(parai%allgrp)
    CALL mp_bcast(cv_dyn,SIZE(cv_dyn),parai%io_source,parai%cp_grp)
    CALL mp_bcast_byte(soft_com, size_in_bytes_of(soft_com),parai%io_source,parai%cp_grp)
    CALL mp_bcast(i_meta,parai%io_source,parai%cp_grp)
    CALL mp_bcast(i_cvst,parai%io_source,parai%cp_grp)

    IF (lmeta%tlocalizespin) THEN
       CALL mp_bcast(rcc0,SIZE(rcc0),parai%io_source,parai%cp_grp)
    ENDIF


    ! WRITE RESTART
    IF ((MOD(i_meta,imeta%st_freq) .EQ. 0  .AND. i_cvst .EQ. 0)&
         .OR. soft_com%exsoft .OR. lmetares) THEN
       lmetares = .TRUE.
       IF (i_meta .GT. 1) THEN
          CALL wmtdresm(50,ntot_iter,i_meta)
       ENDIF

    ENDIF

    ! QUENCH BO
    IF (.NOT.cntl%tmdbo.AND.&
         MOD(i_meta,imeta%qw_freq) .EQ. 0  .AND. i_cvst .EQ. 0)&
         lquench = .TRUE.

    IF (FINAL .AND. paral%parent ) THEN
       ! Deallocate
       DEALLOCATE(hc_last,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(cv_path,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(cscl_val,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(hllw_val,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(hllh_val,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(cv_dyn,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(cv_vel,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(cvd_ist,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)

       DEALLOCATE(f_aver,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(f_harm,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(f_hill,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(f_wall,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)

       DEALLOCATE(fi_harm,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(cv_disp,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(hh_test,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(cv_scl,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(cvd_scl,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)

       DEALLOCATE(disp_sys,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(gpot_sys,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)

       DEALLOCATE(isys,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    ! Halt time for the routine
    CALL tihalt(' META_EXTL',isub)

    RETURN
  END SUBROUTINE meta_ext_mul
  ! ==================================================================





END MODULE meta_exl_mult_utils
