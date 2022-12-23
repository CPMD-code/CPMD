MODULE meta_exlagr_methods
  USE cnst,                            ONLY: au_kcm,&
                                             kboltz,&
                                             pi,&
                                             scmass
  USE cnst_dyn,                        ONLY: &
       cscl_fac, cscl_val, cv_associated, cv_dtemp, cv_dyn, cv_dyn_0, cv_ist, &
       cv_langamma, cv_langevintemp, cv_mass, cv_path, cv_temp, cv_temp0, &
       cv_vel, det_celvar, det_colvar, dstrcv, ekincv, ekincv_walk, fhills, &
       fmtdres, hllh_val, hllw_val, iangcv, ibound, imeta, initial_value, &
       inter_hill, kharm, lchekharm, lcvtc, lkfix, lmeta, ncolvar, ra, rcc0, &
       rmeta, skiphill, tcvlangevin, toll_avcv, tycvar, vbound, vharm, &
       vharm_walk
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
  USE meta_colvar_util_utils,          ONLY: cv_diffusion
  USE meta_exlagr_utils,               ONLY: &
       calc_pos_dyn, cv_exlagr_out, cv_exlagr_out2, cv_read_out, rinvelcv, &
       rmtdres, rscvelcv, scf_tune, wmtdres
  USE meta_hpot_utils,                 ONLY: hills,&
                                             hills_lor,&
                                             hills_ratio,&
                                             hills_sals,&
                                             hills_sals_shift,&
                                             setwall,&
                                             setwall_old
  USE meta_multiple_walkers_utils,     ONLY: calc_pos_dyn_mw,&
                                             cv_exlagr_out_mw,&
                                             mw_filename,&
                                             rinvelcv_mw,&
                                             rmtdres_mw,&
                                             wmtdres_mw
  USE mm_dim_utils,                    ONLY: mm_dim
  USE mm_dimmod,                       ONLY: mm_go_mm,&
                                             mm_revert
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum,&
                                             mp_sync
  USE mw,                              ONLY: mwi
  USE parac,                           ONLY: parai,&
                                             paral
  USE pimd,                            ONLY: grandparent,&
                                             supergroup,&
                                             supersource
  USE prng_utils,                      ONLY: repprngg
  USE puttau_utils,                    ONLY: gettau
  USE quenbo_utils,                    ONLY: give_scr_quenbo
  USE readsr_utils,                    ONLY: xstring
  USE ropt,                            ONLY: infi,&
                                             iteropt
  USE soft,                            ONLY: soft_com
  USE strs,                            ONLY: alpha,&
                                             beta
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             cntr,&
                                             maxsys
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE tpar,                            ONLY: dt_ions
  USE zeroing_utils,                   ONLY: zeroing
#include "sizeof.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: meta_extlagr
  PUBLIC :: meta_extlagr_mw
  PUBLIC :: give_scr_meta_extlagr

CONTAINS

  ! ==================================================================
  SUBROUTINE meta_extlagr(taup,velp,tscr,lquench,lmetares,&
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
    ! HLLH      : hills height
    ! HLLW      : hills amplitude 
    ! CSCL_FAC  : scale factors required for an homogeneous dynamics in the
    ! NCOLVAR dimensions of the CV space 
    ! they can be read from input  and they  can be tuned on the
    ! amplitudes of the variables fluctuations
    ! IBOUND    : 1 if a wall constrains the CV within a certain region of space
    ! otherwise 0 (integer, dimension NCOLVAR)
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
    ! HLLW_VAL  : hills amplitude (dim NTOT_ITER)
    ! HLLH_VAL  : hills altitude  (dim NTOT_ITER)
    ! F_AVER    : accumulator for the calculation of the average force
    ! acting on the CV, as a contribution of the 
    ! underlaying potential 
    ! IF_AVER   : counter for the calculation of the average
    ! IW_CV     : counter for writing the CV values on the cvmdck_mtd file
    ! in order to monitor the behavior along the cntl%md trajectory
    ! DISPLACEMENT: norm of the displacment in CV space
    ! CV_MASS   : masse of CV, if not assigned in 
    ! the input file, they are automathically calculated
    ! dim (NCOLVAR)
    ! KHARM     : spring constant for the harmonic potential, if not
    ! assigned in input file, default values are chosen
    ! dim (NCOLVAR)
    ! ==--------------------------------------------------------------==


    REAL(real_8)                             :: taup(:,:,:), velp(:,:,:), &
                                                tscr(:,:,:)
    LOGICAL                                  :: lquench, lmetares, resetcv
    REAL(real_8)                             :: ekinc, ekinp

    CHARACTER(*), PARAMETER                  :: procedureN = 'meta_extlagr'

    CHARACTER(len=10)                        :: chnum
    CHARACTER(len=100)                       :: lineform
    CHARACTER(len=20)                        :: file9 = 'cvmdck_mtd'
    INTEGER                                  :: i, ia, icv, idof, ie, ierr, &
                                                is, isub, ntot_iter
    INTEGER, SAVE                            :: i_cvst = 0, i_meta = 0, &
                                                i_temp = 0, if_aver = 0, &
                                                ifirst = 0, ifirstp = 0, &
                                                iw_cv = 0
    LOGICAL                                  :: ferror, FINAL, hill_add, lskip
    REAL(real_8)                             :: dif_cv, disp2, displacement, &
                                                fact, fact2, hh, hh_max, &
                                                hh_min, hhh, hhht, hllw2, &
                                                maxf2, tollm, twopi
    REAL(real_8), ALLOCATABLE                :: cv_disp(:), fi_harm(:), &
                                                hh_test(:)
    REAL(real_8), ALLOCATABLE, SAVE          :: cv_scl(:), cvd_ist(:), &
                                                cvd_scl(:), f_aver(:), &
                                                f_harm(:), f_hill(:), &
                                                f_wall(:), hc_last(:)
    REAL(real_8), SAVE                       :: hllh0, temp1 = 0.0_real_8, &
                                                temp2 = 0.0_real_8, tempp

    CALL tiset(' META_EXTL',isub)
    FINAL = .FALSE.
    twopi=2.0_real_8*pi
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

       IF (paral%io_parent) THEN
          ALLOCATE(cv_scl(ncolvar),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
          ALLOCATE(cvd_scl(ncolvar),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
       ENDIF

    ENDIF
    IF (.NOT.paral%io_parent) GOTO 9999
    ! ==--------------------------------------------------------------==
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

    ntot_iter = imeta%i_meta_max+imeta%i_meta_res
    tollm     = toll_avcv


    IF (ifirst .EQ. 0) THEN

       ! Allocate Memory
       CALL dcopy(ncolvar,cv_ist,1,cv_dyn,1)
       DO i=1,ncolvar! ale-cmb
          IF (initial_value(i).EQV. .TRUE.) cv_dyn(i)=cv_dyn_0(i)
       ENDDO
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
       ALLOCATE(ra(ncolvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(ra)!,ncolvar)

       ALLOCATE(cv_path(imeta%i_meta_max+imeta%i_meta_res,ncolvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(cv_path)!,ncolvar*ntot_iter)
       ALLOCATE(cscl_val(imeta%i_meta_res+imeta%i_meta_max,ncolvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(cscl_val)!,ncolvar*ntot_iter)
       ALLOCATE(hllw_val(ntot_iter,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(hllw_val)!,ntot_iter )
       ALLOCATE(hllh_val(ntot_iter,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(hllh_val)!,ntot_iter )


       ! If META_RESTART read data from cntl%proper input
       IF (lmeta%meta_restart .AND. .NOT. lmeta%tcvanalysis) THEN
          IF (lmeta%tresfile) THEN
             CALL rmtdres(50,cv_disp,.TRUE.,ntot_iter)
          ELSE
             CALL cv_read_out(cv_disp)
          ENDIF! TRESFILE
#ifdef __SR11000
          !poption parallel, tlocal(ICV)
          !voption indep(HC_LAST,CV_PATH,CV_DISP)
#else
          !$omp parallel do private(ICV)
#endif
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
                cv_mass(icv) = fact2 /&
                     ((cscl_fac(1,icv)*rmeta%hllw)*(cscl_fac(1,icv)*rmeta%hllw))
                IF (paral%io_parent)&
                     WRITE(6,'(A,I5,A,f15.4,A)') 'Default Mass for CV #',icv,&
                     ' is assigned  = ',&
                     cv_mass(icv)/scmass,&
                     ' a.m.u.(a.u.)^2/(u.CV)^2'
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
          temp1=cv_temp0
          IF (i_meta .EQ. 1) THEN
             CALL rinvelcv(cv_vel,cv_mass,cscl_fac,ncolvar,&
                  cv_temp0,ekincv)
             IF (paral%io_parent)&
                  WRITE(6,'(2x,A,f10.4)')&
                  'CV velocity have been initialized at T = ',cv_temp0
             IF (paral%io_parent)&
                  WRITE(6,'(2x,A,f10.4)')&
                  'the initial CV kinetic energy is = ',ekincv
             cv_temp0=temp1
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
          IF (paral%io_parent)&
               CALL fileopen(444,file9,fo_app+fo_verb,ferror)
       ENDIF

    ENDIF

    ! ==--------------------------------------------------------------==
    ! Calculate the Collective Variables and their Derivatives
    CALL colvarofr(taup,tscr)

    ! mb - do we really need to treat the spin CV differently ?
    ! ale  IF(RESETCV .AND. TLOCALIZESPIN) THEN
    ! ale    DO ICV = 1,NCOLVAR
    ! ale      IF(ICV_SPIN(ICV) .NE. 0) THEN
    ! ale        CV_DYN(ICV) = CV_IST(ICV)
    ! ale      ENDIF
    ! ale    ENDDO
    ! ale    RESETCV = .FALSE.
    ! ale  ENDIF
    ! ==--------------------------------------------------------------==
    ! Initialization

    IF (ifirst.EQ.0) THEN
       IF (i_meta.EQ.1)THEN
          !$omp parallel do private(ICV)
          DO icv = 1,ncolvar
             hc_last(icv) = cv_ist(icv)
          ENDDO
       ENDIF
       hllh0 = rmeta%hllh
       ifirst=1
       ! mb - do we really need to treat the spin CV differently ?
       ! ale    IF(.NOT. META_RESTART .AND. TLOCALIZESPIN) THEN
       ! ale      CALL DCOPY(NCOLVAR,CV_IST,1,CV_DYN,1)
       ! ale    ENDIF
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
    CALL calc_pos_dyn(tollm,i_cvst,i_temp,i_meta,hc_last,hill_add)

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
             dif_cv = dif_cv - twopi
          ELSEIF (iangcv(icv) .EQ. 1 .AND. dif_cv .LT. -pi) THEN
             dif_cv = dif_cv + twopi
          ENDIF
          IF (hh_test(icv).LE.1.e+1_real_8 .AND. hh_test(icv).GT.0.2_real_8) THEN
             fact = ABS(dif_cv)/cscl_fac(1,icv)
             IF (fact.LT. rmeta%hllw) THEN
                kharm(icv) = ABS(hh_test(icv)/2._real_8/rmeta%hllw)
             ELSE
                kharm(icv) = ABS(hh_test(icv)/2._real_8/dif_cv)
             ENDIF
          ELSE
          ENDIF
       ENDDO
    ENDIF

    ! ==--------------------------------------------------------------==
    ! Calculation of the forces on the collective variables
    ! This is used in order to determine the optimal hills parameters 
    ! In F_aver the forces on CV are accumulated in order to calculate
    ! their average value when the meta step is concluded
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
    ! metastep is reported, afterwards the values of the two CV sets and 
    ! instantaneous temperature in K
    IF (lmeta%tcvmonitor .AND. MOD(iw_cv,imeta%wcv_freq) .EQ. 0) THEN
       IF (paral%io_parent)&
            WRITE(chnum,'(I5)') ncolvar
       CALL xstring(chnum,ia,ie)
       lineform=&
            '(1X,2I9,'//chnum(ia:ie)//&
            'E14.6,'//chnum(ia:ie)//'E14.6,E14.6)'
       IF (paral%io_parent)&
            WRITE(444,lineform)iw_cv,iteropt%nfi ,(cv_ist(icv),icv=1,ncolvar),&
            (cv_dyn(icv),icv=1,ncolvar),tempp
       CALL m_flush(444)
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

          displacement = 0.0_real_8
          DO icv = 1,ncolvar
             dif_cv = cv_dyn(icv)-hc_last(icv)
             IF (iangcv(icv) .EQ. 1 .AND. dif_cv .GT. pi) THEN
                dif_cv = dif_cv - twopi
             ELSEIF (iangcv(icv) .EQ. 1 .AND. dif_cv .LT. -pi) THEN
                dif_cv = dif_cv + twopi
             ENDIF
             fact = dif_cv/cscl_fac(1,icv)
             displacement = displacement + fact*fact
          ENDDO
          disp2 = displacement
          displacement = SQRT(displacement)

          ! ==--------------------------------------------------------------==
          ! Tuning Hills parameters 

          IF (lmeta%ltune_hh) THEN
             hh_min = 1000._real_8
             hh_max = 0.0_real_8
             lskip = .FALSE.
             hllw2 = rmeta%hllw*rmeta%hllw

             hhht = 0.0_real_8
             DO icv = 1,ncolvar
                f_aver(icv)  = f_aver(icv) / REAL(if_aver,kind=real_8)
                dif_cv = cv_dyn(icv)-hc_last(icv)
                IF (iangcv(icv) .EQ. 1 .AND. dif_cv .GT. pi) THEN
                   dif_cv = dif_cv - twopi
                ELSEIF (iangcv(icv) .EQ. 1 .AND. dif_cv .LT. -pi) THEN
                   dif_cv = dif_cv + twopi
                ENDIF

                !vw protect if dif_cv is too small 
                IF(ABS(dif_cv)<1.0e-12_real_8) dif_cv = 1.0_real_8

                fact = cscl_fac(1,icv)*cscl_fac(1,icv)/dif_cv *&
                     EXP(0.5_real_8)*hllw2*disp2/&
                     (disp2*EXP(-0.5_real_8*disp2/hllw2)+&
                     hllw2*(EXP(-0.5_real_8*disp2/hllw2)-&
                     EXP(-0.5_real_8*rmeta%rshift*rmeta%rshift)))

                hhh     = 0.5_real_8*ABS(f_aver(icv)*fact)
                hhht = hhht + hhh

                IF (paral%io_parent)&
                     WRITE(6,'(2x,A,I4,3F14.6,2x,2F14.6)') 'FCV ',&
                     icv,f_aver(icv),hhh,fact

                hh_min       = MIN(hh_min,hhh)
                hh_max       = MAX(hh_max,hhh)
             ENDDO
             hh       =  hhht / REAL(ncolvar,kind=real_8)
             IF (hh_min .GT. rmeta%htop .OR. hh_max .LT. rmeta%hlow) THEN
                lskip = .TRUE.
                IF (paral%io_parent)&
                     WRITE(6,'(A,f10.5,A,f10.5,A)') 'Hill top > ', rmeta%htop*au_kcm,&
                     ' Kcal, or Hill top < ',rmeta%hlow*au_kcm,' Kcal'
             ENDIF
             IF (.NOT. lskip .AND. hh .LT. rmeta%hlow) THEN
                rmeta%hllh = (rmeta%hlow+hh_max)*0.5_real_8
             ELSEIF (.NOT. lskip .AND. hh .GT. rmeta%htop) THEN
                rmeta%hllh = (hh_min+rmeta%htop)*0.5_real_8
             ELSEIF (.NOT. lskip ) THEN
                rmeta%hllh = hh
             ELSE
                rmeta%hllh = rmeta%hllh*0.2_real_8+hllh0*0.8_real_8
             ENDIF
             IF (paral%io_parent)&
                  WRITE(6,'(/A,f10.5,3(3x,A,f10.5)/)')&
                  '||DISP|| = ',displacement,&
                  'HH tuned min = ', hh_min,'HH tuned max = ', hh_max,&
                  ' HH used ',rmeta%hllh
             if_aver = 0
             CALL zeroing(f_aver)!,ncolvar)

             IF (lmeta%thillvol) THEN
                rmeta%hllw = (rmeta%hvol0/rmeta%hllh)**(1.0_real_8/REAL(ncolvar,kind=real_8))
             ENDIF
          ENDIF

          ! ==--------------------------------------------------------------==
          ! Scaling Factors Tuning
          IF (lmeta%ltune_cscl) THEN

             CALL scf_tune(i_meta,hc_last)

          ENDIF
          ! ==--------------------------------------------------------------==
          ! Diffusion calculation
          CALL cv_diffusion(hc_last,cvd_ist,cvd_scl,cv_scl,i_meta,&
               ntot_iter)

          ! ==--------------------------------------------------------------==
          ! Open output file
          CALL cv_exlagr_out(i_meta,hc_last,f_harm,f_hill,f_wall,&
               cvd_ist,displacement,ekinc,ekinp)

          ! ==--------------------------------------------------------------==
          !$omp parallel do private(ICV)
          DO icv = 1,ncolvar
             hc_last(icv)   = cv_dyn(icv)
          ENDDO

       ELSE
          ! ==--------------------------------------------------------------==
          ! Open output file
          CALL cv_exlagr_out2(i_meta,ekinc,ekinp)

          ! ==--------------------------------------------------------------==

       ENDIF
       i_meta = i_meta + 1
       i_temp = 0
       i_cvst = 0
    ENDIF        ! HILL_ADD

    ! ==--------------------------------------------------------------==
    ! Calculation of the dynamics of the collective variables
    ! A velocity Verlet Scheme is used

    ! ==--------------------------------------------------------------==
    ! Update Velocity : first half step
    ! using forces calculated in the previous step
    !$omp parallel do private(ICV,FACT)
    DO icv = 1,ncolvar
       fact = dt_ions/(2._real_8*cv_mass(icv))
       cv_vel(icv) = cv_vel(icv)&
            + fact*(f_harm(icv)+f_hill(icv)/cscl_fac(1,icv)+f_wall(icv))
    ENDDO

    IF (tcvlangevin) THEN
       DO icv=1,ncolvar
          fact = dt_ions/(2._real_8*cv_mass(icv))
          ra(icv)=repprngg()
          cv_vel(icv)=cv_vel(icv)-cv_mass(icv)*cv_vel(icv)*cv_langamma&
               *fact+SQRT(fact*(2._real_8*cv_mass(icv)))*0.5_real_8*SQRT(2._real_8&
               *cv_langevintemp*kboltz/cv_mass(icv)*cv_langamma)*ra(icv)
       ENDDO
    ENDIF

    ! ==--------------------------------------------------------------==
    ! Update Position : second order dynamics

    DO icv = 1,ncolvar
       fact = dt_ions
       IF (iangcv(icv) .EQ. 1 )  THEN
          IF (cv_dyn(icv).LT.0._real_8) cv_dyn(icv)= cv_dyn(icv)+twopi
          IF (cv_dyn(icv).GT.twopi) cv_dyn(icv)= cv_dyn(icv)-twopi
       ELSEIF (iangcv(icv) .EQ. 2 )  THEN
          IF (cv_dyn(icv).LT.0._real_8) cv_dyn(icv)= -cv_dyn(icv)
          IF (cv_dyn(icv).GT. pi) cv_dyn(icv)= pi - (cv_dyn(icv) - pi)
       ENDIF
       cv_dyn(icv) = cv_dyn(icv) + fact*cv_vel(icv)
    ENDDO

    ! ==--------------------------------------------------------------==
    ! Calculate forces with the new positions

200 CONTINUE
    ! ==--------------------------------------------------------------==
    ! Calculation of the contributions to the forces on ions 
    ! due to the coupling with the collective variables
    ! V_ARM = SUM_icv (K_icv(CV_IST_icv - CV_DYN_icv))^2

    CALL zeroing(fi_harm)!,cotc0%nodim)
    DO icv = 1,ncolvar
       dif_cv = cv_ist(icv)-cv_dyn(icv)
       IF (iangcv(icv) .EQ. 1 .AND. dif_cv .GT. pi) THEN
          dif_cv = dif_cv - twopi
       ELSEIF (iangcv(icv) .EQ. 1 .AND. dif_cv .LT. -pi) THEN
          dif_cv = dif_cv + twopi
       ENDIF
       fact = kharm(icv)*2.0_real_8*dif_cv
       DO idof = 1,cotc0%nodim
          fi_harm(idof) = fi_harm(idof) -&
               fact *  det_colvar(idof,icv)
       ENDDO
    ENDDO

    ! ==--------------------------------------------------------------==
    ! Update forces on ions

    CALL zeroing(tscr)!,3*maxsys%nax*maxsys%nsx)
    CALL gettau(tscr,fi_harm)
#ifdef __SR11000
    !poption parallel, tlocal(IS,IA)
    !voption indep(FHILLS,TSCR)
#else
    !$omp parallel do private(IS,IA) 
#endif
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          fhills(1,ia,is) = tscr(1,ia,is)
          fhills(2,ia,is) = tscr(2,ia,is)
          fhills(3,ia,is) = tscr(3,ia,is)
       ENDDO
    ENDDO

    ! ==--------------------------------------------------------------==
    ! Update stress tensor: contributions of CV as functions of HT
    IF (lmeta%tcvcell) THEN
       CALL zeroing(dstrcv)!,6)
       DO icv = 1,ncolvar
          dif_cv = cv_ist(icv)-cv_dyn(icv)
          fact = -2.0_real_8*kharm(icv)*dif_cv
          dstrcv(1) = dstrcv(1) +&
               fact*det_celvar(alpha(1),beta(1),icv)
          dstrcv(2) = dstrcv(2) +&
               fact*det_celvar(alpha(2),beta(2),icv)
          dstrcv(3) = dstrcv(3) +&
               fact*det_celvar(alpha(3),beta(3),icv)
          dstrcv(4) = dstrcv(4) +&
               fact*det_celvar(alpha(4),beta(4),icv)
          dstrcv(5) = dstrcv(5) +&
               fact*det_celvar(alpha(5),beta(5),icv)
          dstrcv(6) = dstrcv(6) +&
               fact*det_celvar(alpha(6),beta(6),icv)
       ENDDO
    ENDIF

    ! ==--------------------------------------------------------------==
    ! Calculate contribution of the Harmonic Potential (recall)

    CALL zeroing(f_harm)!,ncolvar)
    DO icv = 1,ncolvar
       dif_cv = cv_ist(icv)-cv_dyn(icv)
       IF (iangcv(icv) .EQ. 1 .AND. dif_cv .GT. pi) THEN
          dif_cv = dif_cv - twopi
       ELSEIF (iangcv(icv) .EQ. 1 .AND. dif_cv .LT. -pi) THEN
          dif_cv = dif_cv + twopi
       ENDIF
       f_harm(icv) =  kharm(icv)*2.0_real_8*dif_cv
    ENDDO

    ! Calculate hills contribution to forces
    IF (lmeta%lhills) THEN
       CALL zeroing(f_hill)!,ncolvar)
       IF (lmeta%hlore) THEN
          CALL hills_lor(cv_dyn,ncolvar,i_meta,ntot_iter,f_hill,&
               1,0)
       ELSEIF (lmeta%hratio) THEN
          CALL hills_ratio(cv_dyn,ncolvar,i_meta,ntot_iter,f_hill,&
               hc_last,i_cvst,1,0)
       ELSEIF (lmeta%hshift) THEN
          CALL hills_sals_shift(cv_dyn,ncolvar,i_meta,ntot_iter,f_hill,&
               hc_last,i_cvst,1,0)
       ELSE IF (lmeta%sphere) THEN
          CALL hills(cv_dyn,ncolvar,i_meta,ntot_iter,f_hill,&
               1,0)
       ELSE
          CALL hills_sals(cv_dyn,ncolvar,i_meta,ntot_iter,f_hill,&
               hc_last,i_cvst,1,0)
       ENDIF
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

    IF (tcvlangevin) THEN
       DO icv=1,ncolvar
          fact = dt_ions/(2._real_8*cv_mass(icv))
          cv_vel(icv)=cv_vel(icv)-cv_mass(icv)*cv_vel(icv)*cv_langamma&
               *fact+SQRT(fact*(2._real_8*cv_mass(icv)))*0.5_real_8*SQRT(2&
               *cv_langevintemp*kboltz/cv_mass(icv)*cv_langamma)*ra(icv)
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Calculation of   Temperature of CV: CV_TEMP
    ! If CV_TEMP_CONTROL, check the condition for rescaling T

    CALL rscvelcv(temp1,temp2,tempp)

    ! ==--------------------------------------------------------------==
    ! IF MAX_STEP reached or EXIT file exist

    IF (infi.EQ.cnti%nomore .OR. i_meta-1.EQ.ntot_iter) THEN
       soft_com%exsoft = .TRUE.

       CALL  colvarpr

       IF (lmeta%tcvmonitor) THEN
          IF (paral%io_parent)&
               CALL fileclose(444)
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
          CALL wmtdres(50,ntot_iter,i_meta)
       ENDIF
    ENDIF

    ! QUENCH BO
    IF (.NOT.cntl%tmdbo.AND.&
         MOD(i_meta,imeta%qw_freq) .EQ. 0  .AND. i_cvst .EQ. 0)&
         lquench = .TRUE.

    IF (FINAL .AND. paral%io_parent )  THEN
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
       DEALLOCATE(cv_scl,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(cvd_scl,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent) THEN
       DEALLOCATE(fi_harm,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(cv_disp,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(hh_test,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Halt time for the routine
    CALL tihalt(' META_EXTL',isub)

    RETURN
  END SUBROUTINE meta_extlagr
  ! 
  ! ==================================================================
  SUBROUTINE give_scr_meta_extlagr(lmtd,tag)
    ! ==================================================================
    ! Estimates the scratch length necessary for META_EXLAGR
    ! ==================================================================
    INTEGER                                  :: lmtd
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: lquenbo
    LOGICAL                                  :: oldstatus

    lmtd=0
    lquenbo=0
    IF (.NOT.lmeta%lextlagrange)RETURN
    CALL mm_dim(mm_go_mm,oldstatus)
    lmtd=3*ions1%nat+4*ncolvar+2*ncolvar*ncolvar+3*ions1%nat+3*ions1%nat*ncolvar
    CALL mm_dim(mm_revert,oldstatus)
    CALL give_scr_quenbo(lquenbo,tag)
    lmtd=MAX(lmtd,lquenbo)
    RETURN
  END SUBROUTINE give_scr_meta_extlagr
  ! ==================================================================

  ! ==================================================================
  ! MULTIPLE WALKERS
  ! ==================================================================
  SUBROUTINE meta_extlagr_mw(taup,velp,tscr,lquench,&
       lmetares,resetcv,ekinc,ekinp)
    ! ==--------------------------------------------------------------==

    REAL(real_8)                             :: taup(:,:,:), velp(:,:,:), &
                                                tscr(:,:,:)
    LOGICAL                                  :: lquench, lmetares, resetcv
    REAL(real_8)                             :: ekinc, ekinp

    CHARACTER(*), PARAMETER                  :: procedureN = 'meta_extlagr_mw'

    CHARACTER(len=10)                        :: chnum
    CHARACTER(len=100)                       :: lineform
    CHARACTER(len=20)                        :: file9 = 'cvmdck_mtd'
    INTEGER                                  :: i, i_add, I_META_tmp, ia, &
                                                icv, idof, ie, ierr, ipw, &
                                                isub, iwalk, iwalk1, iwalk2, &
                                                j, k, ntot_iter
    INTEGER, ALLOCATABLE, SAVE               :: i_cvst_mw(:), i_temp_mw(:)
    INTEGER, SAVE                            :: i_cvst = 0, i_meta = 0, &
                                                ifirst = 0, ifirstp = 0, &
                                                iw_cv = 0
    LOGICAL                                  :: ferror, FINAL, hill_add
    LOGICAL, ALLOCATABLE, SAVE               :: hill_add_mw(:)
    REAL(real_8) :: dif_cv, disp2, displacement, displacement_mw(mwi%nwalk), &
      fact, fact2, hllh_bak, maxf2, tollm
    REAL(real_8), ALLOCATABLE                :: cv_disp(:), fi_harm(:,:), &
                                                hh_test(:)
    REAL(real_8), ALLOCATABLE, SAVE :: cv_ist_scr(:,:), cv_scl(:), &
      cvd_ist(:), cvd_scl(:), det_colvar_mw(:,:,:), f_harm(:), f_hill(:), &
      f_wall(:), hc_last(:), taupw(:,:,:,:), taupw_scr(:,:,:,:)
    REAL(real_8), SAVE                       :: hllh0, temp1 = 0.0_real_8, &
                                                temp2 = 0.0_real_8, tempp

    CALL mp_sync(supergroup)
    CALL tiset('META_EX_MW',isub)
    !
    ! 
    ! TODO the following features have to be fixed!
    IF (grandparent)THEN
       IF (.NOT.lkfix)CALL stopgm(procedureN,&
            'MULTIPLE WALKER AND ADAPTIVE K NOT IMPLEMENTED',& 
            __LINE__,__FILE__)
       IF (lmeta%ltune_hh)CALL stopgm(procedureN,&
            'MULTIPLE WALKER AND ADAPTIVE HILL HEIGHT NOT IMPLEMENTED',& 
            __LINE__,__FILE__)
       IF (lmeta%ltune_cscl)CALL stopgm(procedureN,&
            'MULTIPLE WALKER AND SCALING FACTOR TUNING NOT IMPLEMENTED',& 
            __LINE__,__FILE__)
       IF (lmeta%tcvcell)CALL stopgm(procedureN,&
            'MULTIPLE WALKER AND CV CELL NOT IMPLEMENTED',& 
            __LINE__,__FILE__)
       IF (tcvlangevin)CALL stopgm(procedureN,&
            'MULTIPLE WALKER AND LANGEVIN THERMOSTAT NOT IMPLEMENTED',& 
            __LINE__,__FILE__)
       IF ((.NOT.lmeta%sphere).AND.paral%parent.AND.(ifirstp.EQ.0))THEN
          IF (paral%io_parent)&
               WRITE(6,'(/,A)')' !!! WARNING !!!'
          IF (paral%io_parent)&
               WRITE(6,'(A,/)')&
               '!!! MULTIPLE WALKER AND WITH REQUESTED HILL TYPE'//&
               ' NOT FULLY TESTED !!!'
       ENDIF
       ! 
    ENDIF
    !
    iwalk=mwi%walker_id
    FINAL = .FALSE.
    ! ==--------------------------------------------------------------==
    IF (ifirstp .EQ. 0) THEN
       i_meta  = imeta%i_meta_res+1
       i_cvst = 0
       cv_associated = .TRUE.
       ifirstp = 1
       ! ALLOCATE MEMORY
       ALLOCATE(i_temp_mw(mwi%nwalk),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(i_temp_mw)!,mwi%nwalk)

       ALLOCATE(i_cvst_mw(mwi%nwalk),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(i_cvst_mw)!,mwi%nwalk)

       ALLOCATE(cv_dyn(ncolvar*mwi%nwalk),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(taupw_scr(3,maxsys%nax,maxsys%nsx,mwi%nwalk),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(taupw(3,maxsys%nax,maxsys%nsx,mwi%nwalk),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(ekincv_walk(mwi%nwalk),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(vharm_walk(mwi%nwalk),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(cv_scl(ncolvar*mwi%nwalk),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(cvd_scl(ncolvar*mwi%nwalk),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
    ENDIF
    CALL zeroing(ekincv_walk)!,mwi%nwalk)
    CALL zeroing(vharm_walk)!,mwi%nwalk)
    ! 
    CALL zeroing(taupw_scr)!,maxsys%nsx*maxsys%nax*3*mwi%nwalk)
    CALL zeroing(taupw)!,maxsys%nsx*maxsys%nax*3*mwi%nwalk)
    IF (paral%parent)THEN
       DO i=1,3
          DO j=1,maxsys%nax
             DO k=1,maxsys%nsx
                taupw_scr(i,j,k,iwalk)=taup(i,j,k)
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    ! 
    CALL mp_sum(taupw_scr,taupw,3*maxsys%nax*maxsys%nsx*mwi%nwalk,supergroup)
    ! 
    ! MW   all processors are zeroed. Needed for final MY_COMBINE
    CALL zeroing(tscr)!,3*maxsys%nax*maxsys%nsx)

    IF (.NOT.grandparent) GOTO 9999
    ! ==--------------------------------------------------------------==


    ntot_iter = imeta%i_meta_max+imeta%i_meta_res
    tollm     = toll_avcv

    ! TODO align for BG
    IF (.NOT.ALLOCATED(fi_harm)) ALLOCATE(fi_harm(cotc0%nodim, mwi%nwalk),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    IF (.NOT.ALLOCATED(cv_disp)) ALLOCATE(cv_disp(ncolvar*mwi%nwalk),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    IF (.NOT.ALLOCATED(hh_test)) ALLOCATE(hh_test(ncolvar),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    IF (ifirst .EQ. 0) THEN

       ! Allocate Memory
       ALLOCATE(cv_vel(ncolvar*mwi%nwalk),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(cv_vel)!,ncolvar*mwi%nwalk)
       ALLOCATE(cvd_ist(ncolvar*mwi%nwalk),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(cvd_ist)!,ncolvar*mwi%nwalk)
       ALLOCATE(hc_last(ncolvar*mwi%nwalk),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(hc_last)!,ncolvar*mwi%nwalk)
       ALLOCATE(f_harm(ncolvar*mwi%nwalk),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(f_harm)!,ncolvar*mwi%nwalk)
       ALLOCATE(f_hill(ncolvar*mwi%nwalk),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(f_hill)!,ncolvar*mwi%nwalk)
       ALLOCATE(f_wall(ncolvar*mwi%nwalk),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(f_wall)!,ncolvar*mwi%nwalk)

       ALLOCATE(cv_path(imeta%i_meta_max+imeta%i_meta_res,ncolvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(cv_path)!,ncolvar*ntot_iter)
       ALLOCATE(cscl_val(imeta%i_meta_res+imeta%i_meta_max,ncolvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(cscl_val)!,ncolvar*ntot_iter)
       ALLOCATE(hllw_val(ntot_iter,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(hllw_val)!,ntot_iter)
       ALLOCATE(hllh_val(ntot_iter,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(hllh_val)!,ntot_iter)
       ALLOCATE(det_colvar_mw(cotc0%nodim,ncolvar,mwi%nwalk),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(det_colvar_mw)!,cotc0%nodim*ncolvar*mwi%nwalk)
       ALLOCATE(cv_ist_scr(ncolvar,mwi%nwalk),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(cv_ist_scr)!,ncolvar*mwi%nwalk)
       ALLOCATE(hill_add_mw(mwi%nwalk),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(skiphill(mwi%nwalk),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       DO iwalk1=1,mwi%nwalk
          hill_add_mw(iwalk1)=.FALSE.
          skiphill(iwalk1)=.FALSE.! SKIP THE FIRST HILL; SETTING HEIGHT TO ZERO
       ENDDO

       ! If META_RESTART read data from cntl%proper input
       IF (lmeta%meta_restart .AND. .NOT. lmeta%tcvanalysis) THEN
          DO iwalk1=1,mwi%nwalk
             CALL mw_filename('MTD_RESTART_',fmtdres,iwalk1)
             IF (lmeta%tresfile) THEN
                ! here the cv_path should read as cv_path(1:nmtd,1:ncolvar)
                ! and thus reading all the positions in one shot
                ! cv_disp and cv_vel are read seperately from each MTD_RESTART_i file
                CALL rmtdres_mw(50,cv_disp((iwalk1-1)*ncolvar+1),&
                     .TRUE.,ntot_iter,iwalk1)
             ELSE
                ! MW TODO
                IF (paral%io_parent)&
                     WRITE(6,*)"!! WARNING!!! READING MTD OUTPUT FILES"//&
                     " FOR RESTART IS YET NOT TESTED !!"
                CALL cv_read_out(cv_disp)
             ENDIF! TRESFILE
          ENDDO
          ! reconstructing the hc_last from the diffusin and the last hill position
          DO iwalk1=1,mwi%nwalk
             DO icv = 1,ncolvar
                ipw=(iwalk1-1)*ncolvar
                hc_last(icv+ipw)&
                     =cv_path(imeta%i_meta_res,icv)+cv_disp(icv+ipw)
             ENDDO
          ENDDO! IWALK1
          CALL dcopy(ncolvar*mwi%nwalk,hc_last,1,cv_dyn,1)
       ENDIF ! META_RESTART



       IF (.NOT. lmeta%tcvanalysis) THEN
          ! Masses
          IF (cv_temp0 .NE. 0.0_real_8) THEN
             fact = cv_temp0*0.3166734e-5_real_8*10._real_8
          ELSE
             fact = cntr%tempw   *0.3166734e-5_real_8*10._real_8
          ENDIF
          fact2 = fact * cntr%delt_ions*REAL(inter_hill,kind=real_8)*&
               cntr%delt_ions*REAL(inter_hill,kind=real_8)
          DO icv = 1,ncolvar
             IF (cv_mass(icv) .EQ. 0.0_real_8) THEN
                cv_mass(icv) = fact2 /&
                     ((cscl_fac(1,icv)*rmeta%hllw)*&
                     (cscl_fac(1,icv)*rmeta%hllw))
                IF (paral%io_parent)&
                     WRITE(6,'(A,I5,A,f15.4,A)') 'Default Mass for CV #',icv,&
                     ' is assigned  = ',&
                     cv_mass(icv)/scmass,&
                     ' a.m.u.(a.u.)^2/(u.CV)^2'
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
                     ((cscl_fac(1,icv)*rmeta%hllw)*&
                     (cscl_fac(1,icv)*rmeta%hllw))
                IF (paral%io_parent)&
                     WRITE(6,'(A,I5,A,f10.4,A)')&
                     'Default harmonic constant for CV #',icv,&
                     ' is assigned ',kharm(icv),' a.u./(u.CV)^2'
             ENDIF
             IF (.NOT. lkfix) lchekharm = .TRUE.
          ENDDO! ICV

          DO iwalk1=1,mwi%nwalk
             ipw=(iwalk1-1)*ncolvar
             ! Initialization of CV Velocity
             IF ((i_meta.EQ.1).OR.lmeta%randwalk) THEN
                CALL rinvelcv_mw(cv_vel(ipw+1),cv_mass,cscl_fac,ncolvar,&
                     cv_temp0,ekincv)
                IF (paral%io_parent)&
                     WRITE(6,'(2x,A,i10,A,f10.4)')'WALKER:',iwalk1,&
                     'CV velocity have been initialized at T = ',cv_temp0
                IF (paral%io_parent)&
                     WRITE(6,'(2x,A,i10,A,f10.4)')'WALKER:',iwalk1,&
                     'the initial CV kinetic energy is = ',ekincv
             ENDIF
             ekincv_walk(iwalk1)=ekincv

             IF (lcvtc)THEN
                temp1=cv_temp+cv_dtemp
                temp2=cv_temp-cv_dtemp
             ENDIF
             tempp = cv_temp0
          ENDDO! IWALK1
       ELSE
          CALL zeroing(cv_vel)!,ncolvar*mwi%nwalk)
       ENDIF

       ! ==--------------------------------------------------------------==
       IF (lmeta%tcvmonitor) THEN
          file9='cvmdck_mtd'
          IF (paral%io_parent)&
               CALL fileopen(444,file9,fo_app+fo_verb,ferror)
       ENDIF

       IF ((i_meta.EQ.1).OR.lmeta%skiphill_mw)THEN
          DO iwalk1=2,mwi%nwalk
             skiphill(iwalk1)=.TRUE.! skip the first hill; setting height to zero
          ENDDO
       ENDIF

    ENDIF


    ! ==--------------------------------------------------------------==
    ! Calculate the Collective Variables and their Derivatives
    DO iwalk1=1,mwi%nwalk
       CALL colvarofr(taupw(:,:,:,iwalk1),tscr)
       CALL dcopy(ncolvar,cv_ist,1,cv_ist_scr(1,iwalk1),1)
       CALL dcopy(cotc0%nodim*ncolvar,det_colvar,1,&
            det_colvar_mw(1,1,iwalk1),1)
    ENDDO
    CALL dcopy(ncolvar*mwi%nwalk,cv_ist_scr,1,cv_ist,1)

    !JFrenzel Bochum   IF (ifirst.EQ.0)CALL dcopy(ncolvar*mwi%nwalk,cv_ist_scr,1,cv_dyn,1)
    IF (ifirst.EQ.0.AND..NOT.lmeta%meta_restart)&
         CALL dcopy(ncolvar*mwi%nwalk,cv_ist_scr,1,cv_dyn,1)

    ! ==--------------------------------------------------------------==
    ! Initialization

    IF (ifirst.EQ.0) THEN
       IF (i_meta.EQ.1)THEN
          DO iwalk1=1,mwi%nwalk
             ipw=(iwalk1-1)*ncolvar
             DO icv = 1,ncolvar
                hc_last(icv+ipw) = cv_ist_scr(icv,iwalk1)
             ENDDO
          ENDDO
       ENDIF
       hllh0 = rmeta%hllh
       ifirst=1
       GOTO 200

    ENDIF

    hill_add = .FALSE.
    DO iwalk1=1,mwi%nwalk
       hill_add_mw(iwalk1)=.FALSE.
    ENDDO

    IF (lmeta%ttrajovrr) cnti%ntraj=999999999

    i_cvst = i_cvst + 1
    DO iwalk1=1,mwi%nwalk
       i_cvst_mw(iwalk1) = i_cvst_mw(iwalk1) + 1
    ENDDO


    iw_cv=iw_cv+1

    ! if a hill_add is true, it is known to all processors at this point
    CALL calc_pos_dyn_mw(tollm,i_cvst_mw,i_temp_mw,i_meta,hc_last,&
         hill_add_mw)

    IF (lmetares) THEN
       FINAL = .TRUE.
       GOTO 9999
    ENDIF

    ! each parent will write to unit 444 with its own file name assosicated seperately
    IF (lmeta%tcvmonitor .AND. MOD(iw_cv,imeta%wcv_freq) .EQ. 0) THEN
       IF (paral%io_parent)&
            WRITE(chnum,'(I5)') ncolvar*mwi%nwalk
       CALL xstring(chnum,ia,ie)
       lineform=&
            '(1X,2I7,'//chnum(ia:ie)//&
            'E14.6,'//chnum(ia:ie)//'E14.6)'
       IF (paral%io_parent)&
            WRITE(444,lineform)iw_cv,iteropt%nfi ,&
            ((cv_ist_scr(icv,iwalk1),icv=1,ncolvar),iwalk1=1,mwi%nwalk),&
            (cv_dyn(icv),icv=1,ncolvar*mwi%nwalk)
       CALL m_flush(444)
    ENDIF


    ! ==--------------------------------------------------------------==
    ! New Step of Meta-Dynamics in the Space of CV
    ! loop is done on the grand parent
    i_add=-1
    DO iwalk1=1,mwi%nwalk
       IF (hill_add_mw(iwalk1)) THEN
          i_add=i_add+1
          IF (MOD(i_meta,imeta%tr_freq) .EQ. 0 .OR. i_meta .EQ. 1 ) THEN
             IF (lmeta%ttrajovrr) cnti%ntraj=1
          ENDIF
          IF (.NOT. lmeta%tcvanalysis) THEN
             ! ==--------------------------------------------------------------==
             ! Squared Norm of the Displacement in the CV space

             displacement = 0.0_real_8
             CALL zeroing(displacement_mw)!,mwi%nwalk)
             ipw=(iwalk1-1)*ncolvar
             DO icv = 1,ncolvar
                DO iwalk2=1,mwi%nwalk
                   dif_cv = cv_dyn(icv+ipw)-hc_last(icv+(iwalk2-1)*ncolvar)! s - s(t_mtd)
                   IF (iangcv(icv) .EQ. 1 .AND. dif_cv .GT. pi) THEN
                      dif_cv = dif_cv - 2.0_real_8*pi
                   ELSEIF (iangcv(icv) .EQ. 1 .AND. dif_cv .LT. -pi) THEN
                      dif_cv = dif_cv + 2.0_real_8*pi
                   ENDIF
                   fact = dif_cv/cscl_fac(1,icv)
                   displacement = displacement + fact*fact
                   disp2 = displacement
                   displacement = SQRT(displacement)
                   IF (iwalk2.EQ.1)displacement_mw(iwalk1)=displacement
                   displacement_mw(iwalk1)&
                        =MIN(displacement,displacement_mw(iwalk1))
                ENDDO! IWALK2
             ENDDO! ICV

             ! Diffusion calculation
             ipw=(iwalk1-1)*ncolvar+1

             I_META_tmp=i_meta+i_add
             hllh_bak=rmeta%hllh
             IF (skiphill(iwalk1))rmeta%hllh=0.0_real_8! if a hill has to be skipped, height is set to zero
             CALL cv_diffusion(hc_last(ipw),cvd_ist(ipw),cvd_scl(ipw),&
                  cv_scl(ipw),I_META_tmp,&
                  ntot_iter)

             ! ==--------------------------------------------------------------==
             ! Open output file
             ipw=(iwalk1-1)*ncolvar+1
             CALL cv_exlagr_out_mw(I_META_tmp,hc_last(ipw),f_harm(ipw),&
                  f_hill(ipw),f_wall(ipw),&
                  cvd_ist(ipw),displacement_mw(iwalk1),ekinc,ekinp,iwalk1)
             IF (skiphill(iwalk1))THEN
                rmeta%hllh=hllh_bak
                skiphill(iwalk1)=.FALSE.
             ENDIF
             ! ==--------------------------------------------------------------==
             ipw=(iwalk1-1)*ncolvar
             DO icv = 1,ncolvar
                hc_last(icv+ipw)   = cv_dyn(icv+ipw)
             ENDDO
          ELSE
             ! ==--------------------------------------------------------------==
             ! Open output file
             CALL cv_exlagr_out2(i_meta,ekinc,ekinp)

             ! ==--------------------------------------------------------------==

          ENDIF
          i_temp_mw(iwalk1) = 0
          i_cvst = 0
          i_cvst_mw(iwalk1) = 0
       ENDIF       ! HILL_ADD

       ! ==--------------------------------------------------------------==
       ! Calculation of the dynamics of the collective variables
       ! A velocity Verlet Scheme is used

       ! ==--------------------------------------------------------------==
       ! Update Velocity : first half step
       ! using forces calculated in the previous step
       ipw=(iwalk1-1)*ncolvar
       DO icv = 1,ncolvar
          fact = dt_ions/(2._real_8*cv_mass(icv))
          cv_vel(icv+ipw) = cv_vel(icv+ipw)&
               + fact*(f_harm(icv+ipw)+f_hill(icv+ipw)/&
               cscl_fac(1,icv)+f_wall(icv+ipw))
       ENDDO

       ! ==--------------------------------------------------------------==
       ! Update Position : second order dynamics

       DO icv = 1,ncolvar
          fact = dt_ions
          IF (iangcv(icv) .EQ. 1 )  THEN
             IF (cv_dyn(icv+ipw).LT.0._real_8)&
                  cv_dyn(icv+ipw)= cv_dyn(icv+ipw)+2._real_8*pi
             IF (cv_dyn(icv+ipw).GT.2._real_8*pi)&
                  cv_dyn(icv+ipw)= cv_dyn(icv+ipw)-2._real_8*pi
          ELSEIF (iangcv(icv) .EQ. 2 )  THEN
             IF (cv_dyn(icv+ipw).LT.0._real_8)&
                  cv_dyn(icv+ipw)= -cv_dyn(icv+ipw)
             IF (cv_dyn(icv+ipw).GT. pi)&
                  cv_dyn(icv+ipw)= pi - (cv_dyn(icv+ipw) - pi)
          ENDIF
          cv_dyn(icv+ipw) = cv_dyn(icv+ipw) + fact*cv_vel(icv+ipw)
       ENDDO
    ENDDO ! IWALK1=1,NWALK

    ! ==--------------------------------------------------------------==
    ! Calculate forces with the new positions

200 CONTINUE
    ! ==--------------------------------------------------------------==
    ! Calculation of the contributions to the forces on ions 
    ! due to the coupling with the collective variables
    ! V_ARM = SUM_icv (K_icv(CV_IST_icv - CV_DYN_icv))^2

    CALL zeroing(fi_harm)!,cotc0%nodim*mwi%nwalk)
    DO iwalk1=1,mwi%nwalk
       ipw=(iwalk1-1)*ncolvar
       DO icv = 1,ncolvar
          dif_cv = cv_ist_scr(icv,iwalk1)-cv_dyn(icv+ipw)
          IF (iangcv(icv) .EQ. 1 .AND. dif_cv .GT. pi) THEN
             dif_cv = dif_cv -2.0_real_8*pi
          ELSEIF (iangcv(icv) .EQ. 1 .AND. dif_cv .LT. -pi) THEN
             dif_cv = dif_cv +2.0_real_8*pi
          ENDIF
          fact = kharm(icv)*2.0_real_8*dif_cv
          DO idof = 1,cotc0%nodim
             fi_harm(idof,iwalk1) = fi_harm(idof,iwalk1) -&
                  fact *  det_colvar_mw(idof,icv,iwalk1)
          ENDDO
       ENDDO! ICV
    ENDDO ! IWALK1

    ! ==--------------------------------------------------------------==
    ! Update forces on ions

    CALL zeroing(taupw_scr)!,3*maxsys%nax*maxsys%nsx*mwi%nwalk)
    IF (cotc0%nodim > 0) THEN
       DO iwalk1=1,mwi%nwalk
          CALL gettau(taupw_scr(:,:,:,iwalk1),fi_harm(1,iwalk1))
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Calculate contribution of the Harmonic Potential (recall)

    CALL zeroing(f_harm)!,ncolvar*mwi%nwalk)
    DO iwalk1=1,mwi%nwalk
       ipw=(iwalk1-1)*ncolvar
       DO icv = 1,ncolvar
          dif_cv = cv_ist_scr(icv,iwalk1)-cv_dyn(icv+ipw)
          IF (iangcv(icv) .EQ. 1 .AND. dif_cv .GT. pi) THEN
             dif_cv = dif_cv - 2.0_real_8*pi
          ELSEIF (iangcv(icv) .EQ. 1 .AND. dif_cv .LT. -pi) THEN
             dif_cv = dif_cv + 2.0_real_8*pi
          ENDIF
          f_harm(icv+ipw) =  kharm(icv)*2.0_real_8*dif_cv
       ENDDO
    ENDDO ! IWALK1

    ! Calculate hills contribution to forces
    CALL zeroing(f_hill)!,ncolvar*mwi%nwalk)
    DO iwalk1=1,mwi%nwalk
       IF (hill_add_mw(iwalk1))i_meta = i_meta + 1
       ipw=(iwalk1-1)*ncolvar+1
       IF (lmeta%lhills) THEN
          IF (lmeta%hlore) THEN
             CALL hills_lor(cv_dyn(ipw),ncolvar,i_meta,&
                  ntot_iter,f_hill(ipw),&
                  1,0)
          ELSEIF (lmeta%hratio) THEN
             CALL hills_ratio(cv_dyn(ipw),ncolvar,i_meta,&
                  ntot_iter,f_hill(ipw),&
                  hc_last(ipw),i_cvst,1,0)
          ELSEIF (lmeta%hshift) THEN
             CALL hills_sals_shift(cv_dyn(ipw),ncolvar,i_meta,&
                  ntot_iter,f_hill(ipw),&
                  hc_last(ipw),i_cvst,1,0)
          ELSE IF (lmeta%sphere) THEN
             CALL hills(cv_dyn(ipw),ncolvar,i_meta,&
                  ntot_iter,f_hill(ipw),&
                  1,0)
          ELSE
             CALL hills_sals(cv_dyn(ipw),ncolvar,i_meta,&
                  ntot_iter,f_hill(ipw),&
                  hc_last(ipw),i_cvst,1,0)
          ENDIF! HLORE
       ENDIF! LHILLS
    ENDDO ! IWALK1

    ! Calculate walls contribution to forces
    CALL zeroing(f_wall)!,ncolvar*mwi%nwalk)
    DO iwalk1=1,mwi%nwalk
       ipw=(iwalk1-1)*ncolvar
       maxf2 = 0.0_real_8
       DO icv = 1,ncolvar
          IF (ibound(icv) .EQ. 1) THEN
             CALL setwall(cv_dyn(icv+ipw),vbound(1,icv),&
                  f_wall(icv+ipw))
          ELSEIF (ibound(icv) .EQ. 2) THEN
             CALL setwall_old(cv_dyn(icv+ipw),tycvar(icv),vbound(1,icv),&
                  f_wall(icv+ipw))
          ENDIF
          IF (grandparent)THEN
             IF (f_wall(icv) .NE. 0.0_real_8) THEN
                IF (paral%io_parent)&
                     WRITE(6,'(72A)') ('*',i = 1,72)
                IF (paral%io_parent)&
                     WRITE(6,'(10x,A,I5,A,I5)')&
                     'WARNING: Boundaries active for WALKER ',iwalk1,&
                     ' on CV ', icv
                IF (paral%io_parent)&
                     WRITE(6,'(72A)')  ('*',i = 1,72)
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    ! Update Velocity : second half step 
    DO iwalk1=1,mwi%nwalk
       ipw=(iwalk1-1)*ncolvar
       DO icv = 1,ncolvar
          fact = dt_ions/(2._real_8*cv_mass(icv))
          cv_vel(icv+ipw) = cv_vel(icv+ipw)&
               + fact*(f_harm(icv+ipw)+f_hill(icv+ipw)/&
               cscl_fac(1,icv)+f_wall(icv+ipw))
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    ! Calculation of   Temperature of CV: CV_TEMP
    ! If CV_TEMP_CONTROL, check the condition for rescaling T

    CALL rscvelcv(temp1,temp2,tempp)
    ! ==--------------------------------------------------------------==
    ! IF MAX_STEP reached or EXIT file exist

    IF (infi.EQ.cnti%nomore .OR. i_meta-1.EQ.ntot_iter) THEN
       soft_com%exsoft = .TRUE.
       DO iwalk1=1,mwi%nwalk
          CALL dcopy(ncolvar,cv_ist_scr(1,iwalk1),1,cv_ist,1)
          IF (grandparent)CALL  colvarpr
       ENDDO

       IF (lmeta%tcvmonitor) THEN
          IF (paral%io_parent)&
               CALL fileclose(444)
       ENDIF
    ENDIF
9999 CONTINUE

    ! We need to sum the FHILLS on all processors. H^mtd=H_wal1+H_walk2+H_walk3
    IF (.NOT.grandparent)CALL zeroing(taupw_scr)!,3*maxsys%nax*maxsys%nsx*mwi%nwalk) ! Zero on all processors except grandparent
    CALL zeroing(fhills)!,3*maxsys%nax*maxsys%nsx)

    CALL mp_sync(supergroup)
    CALL mp_bcast(cv_dyn,SIZE(cv_dyn),supersource,supergroup)

    CALL mp_sync(supergroup)
    CALL mp_bcast_byte(soft_com, size_in_bytes_of(soft_com),supersource,supergroup)

    CALL mp_sync(supergroup)
    CALL mp_bcast(i_meta,supersource,supergroup)

    CALL mp_sync(supergroup)
    CALL mp_bcast(i_cvst,supersource,supergroup)

    CALL mp_sync(supergroup)
    CALL mp_bcast(i_cvst_mw,SIZE(i_cvst_mw),supersource,supergroup)

    CALL mp_sync(supergroup)
    CALL mp_bcast(ekincv_walk,SIZE(ekincv_walk),supersource,supergroup)
    CALL mp_bcast(vharm_walk,SIZE(vharm_walk),supersource,supergroup)
    ekincv = ekincv_walk(iwalk)
    vharm  = vharm_walk(iwalk)

    CALL mp_sync(supergroup)
    CALL mp_bcast(taupw_scr,SIZE(taupw_scr),supersource,supergroup)
    IF (paral%parent)THEN
       CALL dcopy(3*maxsys%nax*maxsys%nsx,taupw_scr(1,1,1,iwalk),1,fhills,1)
    ENDIF
    CALL mp_sync(supergroup)

    IF (lmeta%tlocalizespin) THEN
       CALL mp_sync(supergroup)
       CALL mp_bcast(rcc0,SIZE(rcc0),supersource,supergroup)
    ENDIF

    ! WRITE RESTART
    DO iwalk1=1,mwi%nwalk
       IF ((MOD(i_meta,imeta%st_freq) .EQ. 0  .AND. i_cvst_mw(iwalk1).EQ.0)&
            .OR. soft_com%exsoft .OR. lmetares) lmetares = .TRUE.
    ENDDO
    ! 
    IF (lmetares.AND.i_meta.GT.1.AND.grandparent)THEN
       CALL dcopy(mwi%nwalk*ncolvar,cv_vel,1,cv_ist_scr,1)
       DO iwalk1=1,mwi%nwalk
          CALL dcopy(ncolvar,cv_ist_scr(1,iwalk1),1,cv_vel,1)
          CALL mw_filename('MTD_RESTART_',fmtdres,iwalk1)
          CALL wmtdres_mw(50,ntot_iter,i_meta,iwalk1)
       ENDDO
       CALL dcopy(mwi%nwalk*ncolvar,cv_ist_scr,1,cv_vel,1)
    ENDIF

    CALL mp_sync(supergroup)
    CALL mp_bcast(lmetares,supersource,supergroup)

    ! QUENCH BO
    IF (.NOT.cntl%tmdbo.AND.&
         MOD(i_meta,imeta%qw_freq) .EQ. 0  .AND. i_cvst .EQ. 0)&
         lquench = .TRUE.

    IF (grandparent) THEN
       IF (ALLOCATED(fi_harm)) DEALLOCATE(fi_harm,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       IF (ALLOCATED(cv_disp)) DEALLOCATE(cv_disp,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       IF (ALLOCATED(hh_test)) DEALLOCATE(hh_test,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    ENDIF

    IF (FINAL)  THEN
       ! Deallocate
       IF (grandparent)THEN
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
          DEALLOCATE(cv_vel,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
          DEALLOCATE(cvd_ist,STAT=ierr)
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

       ENDIF
       DEALLOCATE(i_temp_mw,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(i_cvst_mw,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(cv_dyn,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(taupw_scr,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(taupw,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(cv_scl,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(cvd_scl,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Halt time for the routine
    CALL mp_sync(supergroup)
    CALL tihalt('META_EX_MW',isub)
    RETURN
  END SUBROUTINE meta_extlagr_mw
  ! ==================================================================


END MODULE meta_exlagr_methods
