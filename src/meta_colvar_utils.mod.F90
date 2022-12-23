MODULE meta_colvar_utils
  USE cnst,                            ONLY: au_kcm,&
                                             pi
  USE cnst_dyn,                        ONLY: &
       cscl_fac, cscl_val, cv_dyn, cv_ist, cv_path, det_celvar, det_colvar, &
       dstrcv, fhills, fmtdres, hllh_val, hllw_val, iangcv, ibound, imeta, &
       inter_hill, lmeta, ncolvar, rmeta, toll_avcv, tycvar, vbound
  USE cotr,                            ONLY: cotc0
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_old
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE meta_colvar_inp_utils,           ONLY: colvarofr,&
                                             colvarpr,&
                                             cv_forces3
  USE meta_colvar_util_utils,          ONLY: calc_aver,&
                                             calc_pos,&
                                             cv_diffusion,&
                                             cv_write_out,&
                                             tune_height
  USE meta_exlagr_utils,               ONLY: rmtdres,&
                                             wmtdres
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
                                             rmtdres_mw,&
                                             wmtdres_mw
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum,&
                                             mp_sync
  USE mw,                              ONLY: mwi,&
                                             tmw
  USE parac,                           ONLY: parai,&
                                             paral
  USE pimd,                            ONLY: grandparent,&
                                             supergroup,&
                                             supersource
  USE puttau_utils,                    ONLY: gettau
  USE ropt,                            ONLY: infi
  USE soft,                            ONLY: soft_com
  USE store_types,                     ONLY: rout1
  USE strs,                            ONLY: alpha,&
                                             beta
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             maxsys
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE tpar,                            ONLY: dt_ions,&
                                             dtb2mi
  USE zeroing_utils,                   ONLY: zeroing
#include "sizeof.h"

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: meta_colvar, meta_colvar_mw

CONTAINS

  ! ==================================================================
  SUBROUTINE meta_colvar(taup,velp,fion,tscr,lquench,lmetares,ekinc,ekinp)
    ! ==--------------------------------------------------------------==
    ! ==  Calculates the contributions to the forces on ions          ==
    ! ==  due to the hills that have been accumulated in the space    ==
    ! ==  of the defined collective variables the hills and to the    ==
    ! ==  walls that in that same space limit the wandering of the CV ==
    ! ==  F_HILL  and   F_WALL
    ! ==  The variation of the thajectory due to these forces is also ==
    ! ==  calculated and the ions coordinated updated                 ==
    ! ==--------------------------------------------------------------==

    ! NCOLVAR   : # of collective variables
    ! LHILLS    : true when hills are used in the dynamics of the constraints
    ! I_META    : # of steps in the meta-dynamics
    ! I_META_MAX: max # of steps in the meta-dynamics
    ! I_META_RES: if restart from previous meta dynamics, 
    ! # of steps already done
    ! I_CVST    : # of steps since last HILL has been added
    ! CV_IST    : actual values of collective variables, dim NCOLVAR
    ! CVD_IST   : actual values of diffusion coefficient of CV
    ! HLLH      : hills height 
    ! HLLW      : hills amplitude 
    ! CSCL_FAC  : scale factors required for an homogeneous dynamic in the
    ! NCOLVAR dimensions of the CV space 
    ! they can be read from input and can be tuned on the 
    ! amplitudes of the collective variables fluctuations
    ! IBOUND    : 1 if a wall constraining the CV within a certain region of space
    ! otherwise 0 (integer, dimension NCOLVAR)
    ! VBOUND    : parameters for the definition of constraining potential
    ! (real, dimension 4*NCOLVAR)
    ! if WALL+ VBOUND(2) is the barrier for CV_DYN > VBOUND(1)
    ! if WALL- VBOUND(4) is the barrier for CV_DYN < VBOUND(3)
    ! F_HILL     : force contribution to the ions due to the accumulation 
    ! of hills  in the space of CV 
    ! dim (NODIM)
    ! F_WALL     : force contribution due to the constraining VBOUND
    ! dim (NODIM)
    ! META_RESTART: restart for a previous meta-dynamics
    ! acculators and the constraints path have to be read 
    ! from the right files
    ! the output is appended 
    ! CV_PATH    :  CV history (dim NTOT_ITER*NCOLVAR)
    ! CSCL_VAL   :  CV scaling factors history (dim NTOT_ITER*NCOLVAR)
    ! HLLW_VAL   :  hills amplitude (dim NTOT_ITER)
    ! HLLH_VAL   :  hills altitude  (dim NTOT_ITER)
    ! F_AVER     :  accumulator for the calculation of the average force
    ! acting on the CV, as a contribution of the 
    ! underlaying potential 
    ! IF_AVER    :  counter for the calculation of the average
    ! DISPLACEMENT: norm of the displacment in CV space
    ! ==--------------------------------------------------------------==


    REAL(real_8)                             :: taup(:,:,:), velp(:,:,:), &
                                                fion(:,:,:), tscr(:,:,:)
    LOGICAL                                  :: lquench, lmetares
    REAL(real_8)                             :: ekinc, ekinp

    CHARACTER(*), PARAMETER                  :: procedureN = 'meta_colvar'

    CHARACTER(len=20)                        :: file1 = 'colvar_mtd', &
                                                file2 = 'parvar_mtd'
    INTEGER                                  :: ia, icv, idof, ierr, ii, ij, &
                                                is, isub, ntot_iter
    INTEGER, SAVE                            :: i_cvst = 0, i_meta = 0, &
                                                i_temp = 0, if_aver = 0, &
                                                ifirst = 0
    LOGICAL                                  :: ferror, FINAL, hill_add, &
                                                ionode, lskip, posinst, &
                                                tune1, tune2
    REAL(real_8)                             :: dif_cv, disp2, displacement, &
                                                fact, hh, hh_max, hh_min, &
                                                hhh, hhht, hllw2, maxf2, tollm
    REAL(real_8), ALLOCATABLE                :: cv_aver(:), cv_scl(:), &
                                                cvd_scl(:), f_hill(:), &
                                                f_wall(:), hh_test(:)
    REAL(real_8), ALLOCATABLE, SAVE          :: cv_last(:), cv_store(:,:), &
                                                cvd_ist(:), f1_cv(:), &
                                                f2_cv(:), f_aver(:)
    REAL(real_8), SAVE                       :: hllh0

    CALL tiset(procedureN,isub)
    ! ==--------------------------------------------------------------==
    ionode=paral%io_parent
    IF (tmw) ionode=grandparent
    !TK bugfix: FINAL needs to be initialized on all procs!
    FINAL = .FALSE.
    IF (.NOT.ionode) GOTO 9999

    ! TODO align for BG
    ALLOCATE(f_hill(cotc0%nodim),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL zeroing(f_hill)
    ALLOCATE(f_wall(cotc0%nodim),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL zeroing(f_wall)
    ALLOCATE(cv_aver(ncolvar),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL zeroing(cv_aver)
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

    ntot_iter = imeta%i_meta_max+imeta%i_meta_res
    tollm     = toll_avcv

    IF (ifirst .EQ. 0) THEN
       i_meta  = imeta%i_meta_res+1
       i_temp  = 0
       i_cvst  = 0
       if_aver = 0

       ! Allocate Memory

       ALLOCATE(cvd_ist(ncolvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(cvd_ist)!,ncolvar)
       ALLOCATE(cv_store(inter_hill,ncolvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(cv_store)!,ncolvar*inter_hill)
       ALLOCATE(cv_last(ncolvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(cv_last)!,ncolvar)
       ALLOCATE(f_aver(ncolvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(f_aver)!,ncolvar)
       ALLOCATE(f2_cv(ncolvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(f2_cv)!,ncolvar)
       ALLOCATE(f1_cv(ncolvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(f1_cv)!,ncolvar)

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


       ! If META_RESTART read data from cntl%proper input
       IF (lmeta%meta_restart) THEN
          IF (lmeta%tresfile) THEN
             CALL rmtdres(50,cvd_ist,.TRUE.,ntot_iter)
          ELSE
             IF (paral%io_parent) CALL fileopen(51,file1,fo_old,ferror)
             IF (paral%io_parent) CALL fileopen(52,file2,fo_old,ferror)
             IF (paral%io_parent) REWIND(51)
             IF (paral%io_parent) REWIND(52)
             DO ii = 1,imeta%i_meta_res
                IF (paral%io_parent)&
                     READ(51,err=20,END=20,fmt=*)  ij,&
                     (cv_path(ii,icv),icv=1,ncolvar),&
                     (cscl_val(ii,icv),icv=1,ncolvar)
                ! from output forces and parameters
                IF (paral%io_parent)&
                     READ(52,err=20,END=20,fmt=*) ij,&
                     displacement,hllw_val(ii,1),hllh_val(ii,1)
             ENDDO
             IF (paral%io_parent) CALL fileclose(51)
             IF (paral%io_parent) CALL fileclose(52)
             GOTO 100
20           CONTINUE
             IF (paral%io_parent)&
                  WRITE(6,*) ' ERROR WHILE READING RESTART DATA '
             CALL stopgm('META_COLVAR',' ',& 
                  __LINE__,__FILE__)
100          CONTINUE
          ENDIF! TRESFILE
       ENDIF ! META_RESTART
    ENDIF

    ! ==--------------------------------------------------------------==
    ! Calculate the Collective Variables and their Derivatives
    CALL colvarofr(taup,tscr)

    ! ==--------------------------------------------------------------==
    ! Store CV values and Check averages

    hill_add = .FALSE.
    rout1%rout = .FALSE.
    i_cvst = i_cvst + 1

    posinst = .TRUE.
    tune1   = .FALSE.
    tune2   = .TRUE.

    IF (ifirst .EQ. 0) THEN
       !$omp parallel do private(ICV)
       DO icv = 1,ncolvar
          cv_last(icv) = cv_ist(icv)
       ENDDO
       ifirst = 1
       hllh0 = rmeta%hllh
    ENDIF

    IF (.NOT.posinst) THEN
       CALL calc_aver(tollm,i_cvst,i_temp,i_meta,cv_store,cv_last,&
            cv_aver,hill_add)
    ELSE
       CALL calc_pos(tollm,i_cvst,i_temp,i_meta,cv_last,&
            cv_aver,hill_add)
    ENDIF
    IF (lmetares) THEN
       FINAL = .TRUE.
       GOTO 9999
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    ! Calculation of the forces on the collective variables
    ! This is used in order to determine the optimal hills parameters 
    ! In F_aver the forces on CV are accumulated in order to calculate
    ! their average value when the meta step is concluded

    IF (tune1 .AND. lmeta%ltune_hh ) THEN
       CALL tune_height(f_aver,velp,fion,tscr)
       if_aver = if_aver + 1
    ELSEIF (tune2 .AND. lmeta%ltune_hh ) THEN
       CALL cv_forces3(taup,tscr,velp,hh_test)
       !$omp parallel do private(ICV)
       DO icv = 1,ncolvar
          f_aver(icv) = f_aver(icv) + hh_test(icv)
       ENDDO
       if_aver = if_aver + 1
    ENDIF

    ! ==--------------------------------------------------------------==
    ! New Step of Meta-Dynamics in the Space of CV
    IF (hill_add) THEN
       IF (MOD(i_meta,imeta%tr_freq) .EQ. 0 .OR. i_meta .EQ. 1 ) rout1%rout = .TRUE.

       ! ==--------------------------------------------------------------==
       ! Squared Norm of the Displacement in the CV space

       displacement = 0.0_real_8
       DO icv = 1,ncolvar
          dif_cv = cv_aver(icv)-cv_last(icv)
          IF (iangcv(icv) .EQ. 1 .AND. dif_cv .GT. pi) THEN
             dif_cv = dif_cv - 2.0_real_8*pi
          ELSEIF (iangcv(icv) .EQ. 1 .AND. dif_cv .LT. -pi) THEN
             dif_cv = dif_cv + 2.0_real_8*pi
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
             dif_cv = cv_aver(icv)-cv_last(icv)
             IF (iangcv(icv) .EQ. 1 .AND. dif_cv .GT. pi) THEN
                dif_cv = dif_cv - 2.0_real_8*pi
             ELSEIF (iangcv(icv) .EQ. 1 .AND. dif_cv .LT. -pi) THEN
                dif_cv = dif_cv + 2.0_real_8*pi
             ENDIF

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
                  WRITE(6,'(A,f10.5,A,f10.5)') 'Hill top > ', rmeta%htop*au_kcm,&
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
       ! Diffusion calculation

       CALL cv_diffusion(cv_last,cvd_ist,cvd_scl,cv_scl,i_meta,&
            ntot_iter)

       ! ==--------------------------------------------------------------==

       ! Open output file
       IF (paral%io_parent)&
            CALL cv_write_out(i_meta,cv_last,cv_aver,f1_cv,f2_cv,&
            cvd_ist,cv_scl,cvd_scl,displacement,ekinc,ekinp)

       ! ==--------------------------------------------------------------==
       !$omp parallel do private(ICV)
       DO icv = 1,ncolvar
          cv_last(icv)   = cv_aver(icv)
       ENDDO

       i_meta = i_meta + 1
       i_temp = 0
       i_cvst = 0
    ENDIF        ! HILL_ADD

    ! ==--------------------------------------------------------------==
    ! Calculation of the contributions to the forces on ions 
    ! due to HILLS and WALLS in the space of CV

    ! Calculate hills contribution to forces
    IF (lmeta%lhills) THEN
       CALL zeroing(f_hill)!,cotc0%nodim)
       CALL zeroing(f1_cv)!,ncolvar)
       IF (lmeta%hlore) THEN
          CALL hills_lor(cv_ist,ncolvar,i_meta,ntot_iter,f1_cv,&
               1,0)
       ELSEIF (lmeta%hratio) THEN
          CALL hills_ratio(cv_ist,ncolvar,i_meta,ntot_iter,f1_cv,&
               cv_last,i_cvst,1,0)
       ELSEIF (lmeta%hshift) THEN
          CALL hills_sals_shift(cv_ist,ncolvar,i_meta,ntot_iter,f1_cv,&
               cv_last,i_cvst,1,0)
       ELSE IF (lmeta%sphere) THEN
          CALL hills(cv_ist,ncolvar,i_meta,ntot_iter,f1_cv,&
               1,0)
       ELSE
          CALL hills_sals(cv_ist,ncolvar,i_meta,ntot_iter,f1_cv,&
               cv_last,i_cvst,1,0)
       ENDIF

       DO icv = 1,ncolvar
          DO idof = 1,cotc0%nodim
             f_hill(idof) = f_hill(idof) +&
                  f1_cv(icv)*&
                  det_colvar(idof,icv)/cscl_fac(1,icv)
          ENDDO
       ENDDO
    ENDIF

    ! Calculate walls contribution to forces
    CALL zeroing(f_wall)!,cotc0%nodim)
    CALL zeroing(f2_cv)!,ncolvar)
    maxf2 = 0.0_real_8
    DO icv = 1,ncolvar
       IF (ibound(icv) .EQ. 1) THEN
          CALL setwall(cv_ist(icv),vbound(1,icv),f2_cv(icv))
       ELSEIF (ibound(icv) .EQ. 2) THEN
          CALL setwall_old(cv_ist(icv),tycvar(icv),vbound(1,icv),&
               f2_cv(icv))
       ENDIF
       maxf2 = MAX(maxf2,ABS(f2_cv(icv)))
       IF (f2_cv(icv) .NE. 0.0_real_8) THEN
          DO idof = 1,cotc0%nodim
             f_wall(idof) = f_wall(idof) +&
                  f2_cv(icv)*&
                  det_colvar(idof,icv)/cscl_fac(1,icv)
          ENDDO
       ENDIF
    ENDDO


    ! ==--------------------------------------------------------------==
    ! Update forces on ions

    CALL zeroing(tscr)!,3*maxsys%nax*maxsys%nsx)
    CALL gettau(tscr,f_hill)
#ifdef __SR11000
    !poption parallel, tlocal(IS,IA,FACT)
    !voption indep(DTB2MI,FHILLS,TSCR)
#endif
    DO is=1,ions1%nsp
       fact= dt_ions*dtb2mi(is)
       DO ia=1,ions0%na(is)
          fhills(1,ia,is) = tscr(1,ia,is)
          fhills(2,ia,is) = tscr(2,ia,is)
          fhills(3,ia,is) = tscr(3,ia,is)
       ENDDO
    ENDDO
    IF (maxf2 .GT. 1.0e-12_real_8) THEN
       CALL zeroing(tscr)!,3*maxsys%nax*maxsys%nsx)
       CALL gettau(tscr,f_wall)
#ifdef __SR11000
       !poption parallel, tlocal(IS,IA,FACT)
       !voption indep(DTB2MI,FHILLS,TSCR)
#endif
       DO is=1,ions1%nsp
          fact= dt_ions*dtb2mi(is)
          DO ia=1,ions0%na(is)
             fhills(1,ia,is) = fhills(1,ia,is)+tscr(1,ia,is)
             fhills(2,ia,is) = fhills(2,ia,is)+tscr(2,ia,is)
             fhills(3,ia,is) = fhills(3,ia,is)+tscr(3,ia,is)
          ENDDO
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Update stress tensor: contributions of CV as functions of HT
    IF (lmeta%tcvcell) THEN
       CALL zeroing(dstrcv)!,6)
       DO icv = 1,ncolvar
          fact = f2_cv(icv)/cscl_fac(1,icv)

          DO ii = 1,6
             dstrcv(ii) = dstrcv(ii) +&
                  fact*det_celvar(alpha(ii),beta(ii),icv)
          ENDDO
       ENDDO
    ENDIF


    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    ! IF MAX_STEP reached or EXIT file exist

    IF (infi.EQ.cnti%nomore .OR. i_meta-1.EQ.ntot_iter) THEN
       soft_com%exsoft = .TRUE.

       CALL  colvarpr
    ENDIF
9999 CONTINUE
    CALL mp_sync(parai%allgrp)

    CALL mp_bcast(soft_com%exsoft,parai%io_source,parai%cp_grp)

    CALL mp_bcast(i_meta,parai%io_source,parai%cp_grp)
    CALL mp_bcast(i_cvst,parai%io_source,parai%cp_grp)

    ! WRITE RESTART
    IF ((MOD(i_meta,imeta%st_freq) .EQ. 0 .AND. i_cvst .EQ. 0)&
         .OR. soft_com%exsoft .OR. lmetares) THEN
       lmetares = .TRUE.
       IF (i_meta .GT. 1) THEN
          CALL wmtdres(50,ntot_iter,i_meta)
       ENDIF


    ENDIF

    ! QUENCH BO
    IF (.NOT.cntl%tmdbo.AND.&
         MOD(i_meta,imeta%qw_freq) .EQ. 0 .AND. i_cvst .EQ. 0)&
         lquench = .TRUE.

    ! ==--------------------------------------------------------------==
    IF (FINAL .AND. ionode) THEN
       ! Deallocate
       DEALLOCATE(cv_store,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(cv_last,STAT=ierr)
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

       DEALLOCATE(f_hill,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(f_wall,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(cv_aver,STAT=ierr)
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
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Halt time for the routine
    CALL tihalt(procedureN,isub)

    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE meta_colvar
  ! ==================================================================

  ! Multiple Walker
  ! ==================================================================
  SUBROUTINE meta_colvar_mw(taup,velp,fion,tscr,lquench,lmetares,ekinc,ekinp)
    ! ==--------------------------------------------------------------==
    ! ==  Calculates the contributions to the forces on ions          ==
    ! ==  due to the hills that have been accumulated in the space    ==
    ! ==  of the defined collective variables the hills and to the    ==
    ! ==  walls that in that same space limit the wandering of the CV ==
    ! ==  F_HILL  and   F_WALL
    ! ==  The variation of the thajectory due to these forces is also ==
    ! ==  calculated and the ions coordinated updated                 ==
    ! ==--------------------------------------------------------------==

    REAL(real_8)                             :: taup(:,:,:), velp(:,:,:), &
                                                fion(:,:,:), tscr(:,:,:)
    LOGICAL                                  :: lquench, lmetares
    REAL(real_8)                             :: ekinc, ekinp

    CHARACTER(*), PARAMETER                  :: procedureN = 'meta_colvar_mw'

    CHARACTER(len=20)                        :: file1 = 'colvar_mtd', &
                                                file2 = 'parvar_mtd'
    INTEGER                                  :: i, i_add, i_meta_tmp, ia, &
                                                icv, idof, ierr, ii, ij, ipw, &
                                                is, isub, iwalk, iwalk1, &
                                                iwalk2, j, k, ntot_iter
    INTEGER, ALLOCATABLE, SAVE               :: i_cvst_mw(:), i_temp_mw(:)
    INTEGER, SAVE                            :: i_cvst = 0, i_meta = 0, &
                                                i_temp = 0, ifirst = 0
    LOGICAL                                  :: ferror, FINAL, hill_add
    LOGICAL, ALLOCATABLE, DIMENSION(:), SAVE :: hill_add_mw
    REAL(real_8)                             :: dif_cv, disp2, displacement, &
                                                displacement_mw(mwi%nwalk), &
                                                fact, maxf2, tollm
    REAL(real_8), ALLOCATABLE                :: cv_scl(:), cvd_scl(:), &
                                                f_hill(:), f_wall(:), &
                                                hh_test(:), taupw(:,:,:,:)
    REAL(real_8), ALLOCATABLE, SAVE :: cv_ist_scr(:,:), cv_last(:), &
      cv_store(:,:), cvd_ist(:), det_colvar_mw(:,:,:), f1_cv(:), f2_cv(:), &
      taupw_scr(:,:,:,:)

    CALL mp_sync(supergroup)
    CALL tiset(procedureN,isub)
    ! ==--------------------------------------------------------------==
    iwalk=mwi%walker_id
    FINAL = .FALSE.

    ! TODO align for BG
    ALLOCATE(f_hill(cotc0%nodim),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL zeroing(f_hill)
    ALLOCATE(f_wall(cotc0%nodim),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL zeroing(f_wall)
    ALLOCATE(hh_test(ncolvar),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL zeroing(hh_test)

    ALLOCATE(taupw(3,maxsys%nax,maxsys%nsx,mwi%nwalk),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(taupw)!,maxsys%nsx*maxsys%nax*3*mwi%nwalk)
    ALLOCATE(cv_scl(ncolvar*mwi%nwalk),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL zeroing(cv_scl)
    ALLOCATE(cvd_scl(ncolvar*mwi%nwalk),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL zeroing(cvd_scl)

    IF (paral%parent)THEN
       DO i=1,3
          DO j=1,maxsys%nax
             DO k=1,maxsys%nsx
                taupw(i,j,k,iwalk)=taup(i,j,k)
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    ! 
    CALL mp_sum(taupw,3*maxsys%nax*maxsys%nsx*mwi%nwalk,supergroup)

    IF (ifirst .EQ. 0) THEN
       ALLOCATE(i_cvst_mw(mwi%nwalk),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(i_cvst_mw)!,mwi%nwalk)

       ALLOCATE(taupw_scr(3,maxsys%nax,maxsys%nsx,mwi%nwalk),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(taupw_scr)

       ALLOCATE(cv_dyn(ncolvar*mwi%nwalk),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(cv_dyn)
    END IF

    IF (.NOT.grandparent) GOTO 9999

    ntot_iter = imeta%i_meta_max+imeta%i_meta_res
    tollm     = toll_avcv

    IF (ifirst .EQ. 0) THEN
       i_meta  = imeta%i_meta_res+1
       i_temp  = 0
       i_cvst  = 0

       ! Allocate Memory
       ALLOCATE(cv_ist_scr(ncolvar,mwi%nwalk),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(cv_ist_scr)!,ncolvar*mwi%nwalk)

       ALLOCATE(i_temp_mw(mwi%nwalk),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(i_temp_mw)!,mwi%nwalk)

       ALLOCATE(cvd_ist(ncolvar*mwi%nwalk),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(cvd_ist)!,ncolvar)

       ALLOCATE(cv_store(inter_hill,ncolvar*mwi%nwalk),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(cv_store)!,ncolvar*inter_hill)

       ALLOCATE(cv_last(ncolvar*mwi%nwalk),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(cv_last)!,ncolvar)

       ALLOCATE(f2_cv(ncolvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(f2_cv)!,ncolvar)

       ALLOCATE(f1_cv(ncolvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(f1_cv)!,ncolvar)

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

       ALLOCATE(hill_add_mw(mwi%nwalk),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)

       ALLOCATE(det_colvar_mw(cotc0%nodim,ncolvar,mwi%nwalk),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(det_colvar_mw)!,cotc0%nodim*ncolvar*mwi%nwalk)

       DO iwalk1=1,mwi%nwalk
          hill_add_mw(iwalk1)=.FALSE.
       ENDDO

       ! If META_RESTART read data from cntl%proper input
       IF (lmeta%meta_restart .AND. .NOT. lmeta%tcvanalysis) THEN
          IF (lmeta%tresfile) THEN
             DO iwalk1=1,mwi%nwalk
                CALL mw_filename('MTD_RESTART_',fmtdres,iwalk1)
                CALL rmtdres_mw(50,cvd_ist((iwalk1-1)*ncolvar+1),&
                     .TRUE.,ntot_iter,iwalk1)
             ENDDO
          ELSE
             IF (paral%io_parent) CALL fileopen(51,file1,fo_old,ferror)
             IF (paral%io_parent) CALL fileopen(52,file2,fo_old,ferror)
             IF (paral%io_parent) REWIND(51)
             IF (paral%io_parent) REWIND(52)
             DO ii = 1,imeta%i_meta_res
                IF (paral%io_parent)&
                     READ(51,err=20,END=20,fmt=*)  ij,&
                     (cv_path(ii,icv),icv=1,ncolvar),&
                     (cscl_val(ii,icv),icv=1,ncolvar)
                ! from output forces and parameters
                IF (paral%io_parent)&
                     READ(52,err=20,END=20,fmt=*) ij,&
                     displacement,hllw_val(ii,1),hllh_val(ii,1)
             ENDDO
             IF (paral%io_parent) CALL fileclose(51)
             IF (paral%io_parent) CALL fileclose(52)
             GOTO 100
20           CONTINUE
             IF (paral%io_parent)&
                  WRITE(6,*) ' ERROR WHILE READING RESTART DATA '
             CALL stopgm('META_COLVAR',' ',& 
                  __LINE__,__FILE__)
100          CONTINUE
          ENDIF! TRESFILE
       ENDIF ! META_RESTART

    ENDIF

    ! ==--------------------------------------------------------------==
    ! Calculate the Collective Variables and their Derivatives
    DO iwalk1=1,mwi%nwalk
       CALL colvarofr(taupw(:,:,:,iwalk1),tscr)
       CALL dcopy(ncolvar,cv_ist,1,cv_ist_scr(1,iwalk1),1)
       CALL dcopy(cotc0%nodim*ncolvar,det_colvar,1,det_colvar_mw(1,1,iwalk1),1)
    ENDDO
    CALL dcopy(ncolvar*mwi%nwalk,cv_ist_scr,1,cv_ist,1)
    CALL dcopy(ncolvar*mwi%nwalk,cv_ist_scr,1,cv_dyn,1)

    ! ==--------------------------------------------------------------==
    ! Store CV values and Check averages
    hill_add = .FALSE.
    DO iwalk1=1,mwi%nwalk
       hill_add_mw(iwalk1)=.FALSE.
    ENDDO

    rout1%rout = .FALSE.

    IF (ifirst .EQ. 0) THEN
       DO iwalk1=1,mwi%nwalk
          ipw=(iwalk1-1)*ncolvar
          !$omp parallel do private(ICV)
          DO icv = 1,ncolvar
             cv_last(icv+ipw) = cv_ist_scr(icv,iwalk1)
          ENDDO
       END DO
       ifirst = 1
    ENDIF

    i_cvst = i_cvst + 1
    DO iwalk1=1,mwi%nwalk
       i_cvst_mw(iwalk1) = i_cvst_mw(iwalk1) + 1
    ENDDO

    CALL calc_pos_dyn_mw(tollm,i_cvst_mw,i_temp_mw,i_meta,cv_last,hill_add_mw)

    IF (lmetares) THEN
       FINAL = .TRUE.
       GOTO 9999
    ENDIF

    IF (lmeta%ltune_hh) THEN
       STOP "unimplemented"
    ENDIF

    ! ==--------------------------------------------------------------==
    ! New Step of Meta-Dynamics in the Space of CV
    i_add=-1
    DO iwalk1=1,mwi%nwalk
       IF (hill_add_mw(iwalk1)) THEN
          i_add=i_add+1
          IF (MOD(i_meta,imeta%tr_freq) .EQ. 0 .OR. i_meta .EQ. 1 ) rout1%rout = .TRUE.

          ! ==--------------------------------------------------------------==
          ! Squared Norm of the Displacement in the CV space
          displacement = 0.0_real_8
          CALL zeroing(displacement_mw)!,mwi%nwalk)
          ipw=(iwalk1-1)*ncolvar
          DO icv = 1,ncolvar
             DO iwalk2=1,mwi%nwalk
                dif_cv = cv_dyn(icv+ipw)-cv_last(icv+(iwalk2-1)*ncolvar)! s - s(t_mtd)
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
                displacement_mw(iwalk1)=MIN(displacement,displacement_mw(iwalk1))
             ENDDO! IWALK2
          ENDDO! ICV


          ! Diffusion calculation
          ipw=(iwalk1-1)*ncolvar+1
          i_meta_tmp=i_meta+i_add
          CALL cv_diffusion(cv_last(ipw),cvd_ist(ipw),cvd_scl(ipw),&
               cv_scl(ipw),i_meta_tmp, ntot_iter)

          ! ==--------------------------------------------------------------==
          ! Open output file
          CALL cv_exlagr_out_mw(i_meta_tmp,cv_last(ipw),(/ncolvar*0.0_real_8/),&
               f_hill(ipw),f_wall(ipw),cvd_ist(ipw),displacement_mw(iwalk1),&
               ekinc,ekinp,iwalk1)

          ! ==--------------------------------------------------------------==
          ipw=(iwalk1-1)*ncolvar
          !$omp parallel do private(ICV)
          DO icv = 1,ncolvar
             cv_last(icv+ipw)   = cv_dyn(icv+ipw)
          ENDDO
          i_temp_mw(iwalk1) = 0
          i_cvst = 0
          i_cvst_mw(iwalk1) = 0 
       END IF         ! HILL_ADD
    END DO

    ! ==--------------------------------------------------------------==
    ! Calculation of the contributions to the forces on ions 
    ! due to HILLS and WALLS in the space of CV

    ! Calculate hills contribution to forces
    IF (lmeta%lhills) THEN
       CALL zeroing(f_hill)!,cotc0%nodim)
       CALL zeroing(f1_cv)!,ncolvar)
       DO iwalk1=1,mwi%nwalk
          IF (hill_add_mw(iwalk1)) i_meta = i_meta + 1   
          ipw=(iwalk1-1)*ncolvar+1
          IF (lmeta%hlore) THEN
             CALL hills_lor(cv_ist_scr(:,iwalk1),ncolvar,i_meta,ntot_iter,f1_cv,&
                  1,0)
          ELSEIF (lmeta%hratio) THEN
             CALL hills_ratio(cv_ist_scr(:,iwalk1),ncolvar,i_meta,ntot_iter,f1_cv,&
                  cv_last(ipw),i_cvst,1,0)
          ELSEIF (lmeta%hshift) THEN
             CALL hills_sals_shift(cv_ist_scr(:,iwalk1),ncolvar,i_meta,ntot_iter,f1_cv,&
                  cv_last(ipw),i_cvst,1,0)
          ELSE IF (lmeta%sphere) THEN
             CALL hills(cv_ist_scr(:,iwalk1),ncolvar,i_meta,ntot_iter,f1_cv,&
                  1,0)
          ELSE
             CALL hills_sals(cv_ist_scr(:,iwalk1),ncolvar,i_meta,ntot_iter,f1_cv,&
                  cv_last(ipw),i_cvst,1,0)
          ENDIF

          DO icv = 1,ncolvar
             DO idof = 1,cotc0%nodim
                f_hill(idof) = f_hill(idof) +&
                     f1_cv(icv)* det_colvar_mw(idof,icv,iwalk1)/cscl_fac(1,icv)
             ENDDO
          ENDDO
       ENDDO
    ENDIF

    ! Calculate walls contribution to forces
    CALL zeroing(f_wall)!,cotc0%nodim)
    CALL zeroing(f2_cv)!,ncolvar)
    maxf2 = 0.0_real_8
    DO iwalk1=1,mwi%nwalk
       ipw=(iwalk1-1)*ncolvar
       DO icv = 1,ncolvar
          IF (ibound(icv) .EQ. 1) THEN
             CALL setwall(cv_ist_scr(icv,iwalk1),vbound(1,icv),f2_cv(icv))
          ELSEIF (ibound(icv) .EQ. 2) THEN
             CALL setwall_old(cv_ist_scr(icv,iwalk1),tycvar(icv),vbound(1,icv),&
                  f2_cv(icv))
          ENDIF
          maxf2 = MAX(maxf2,ABS(f2_cv(icv)))
          IF (f2_cv(icv) .NE. 0.0_real_8) THEN
             DO idof = 1,cotc0%nodim
                f_wall(idof) = f_wall(idof) +&
                     f2_cv(icv)*&
                     det_colvar_mw(idof,icv,iwalk1)/cscl_fac(1,icv)
             ENDDO
          ENDIF
       ENDDO
    ENDDO

    ! ==--------------------------------------------------------------==
    ! Update forces on ions

    CALL zeroing(tscr)!,3*maxsys%nax*maxsys%nsx)
    CALL gettau(tscr,f_hill)
#ifdef __SR11000
    !poption parallel, tlocal(IS,IA,FACT)
    !voption indep(DTB2MI,FHILLS,TSCR)
#endif
    DO is=1,ions1%nsp
       fact= dt_ions*dtb2mi(is)
       DO ia=1,ions0%na(is)
          fhills(1,ia,is) = tscr(1,ia,is)
          fhills(2,ia,is) = tscr(2,ia,is)
          fhills(3,ia,is) = tscr(3,ia,is)
       ENDDO
    ENDDO
    IF (maxf2 .GT. 1.0e-12_real_8) THEN
       CALL zeroing(tscr)!,3*maxsys%nax*maxsys%nsx)
       CALL gettau(tscr,f_wall)
#ifdef __SR11000
       !poption parallel, tlocal(IS,IA,FACT)
       !voption indep(DTB2MI,FHILLS,TSCR)
#endif
       DO is=1,ions1%nsp
          fact= dt_ions*dtb2mi(is)
          DO ia=1,ions0%na(is)
             fhills(1,ia,is) = fhills(1,ia,is)+tscr(1,ia,is)
             fhills(2,ia,is) = fhills(2,ia,is)+tscr(2,ia,is)
             fhills(3,ia,is) = fhills(3,ia,is)+tscr(3,ia,is)
          ENDDO
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Update stress tensor: contributions of CV as functions of HT
    IF (lmeta%tcvcell) THEN
       STOP 'unimplemented'
    ENDIF

    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    ! IF MAX_STEP reached or EXIT file exist

    IF (infi.EQ.cnti%nomore .OR. i_meta-1.GE.ntot_iter) THEN
       soft_com%exsoft = .TRUE.

       CALL  colvarpr
    ENDIF
9999 CONTINUE
    CALL mp_sync(supergroup)
    CALL mp_bcast_byte(soft_com, size_in_bytes_of(soft_com),supersource,supergroup)

    CALL mp_sync(supergroup)
    CALL mp_bcast(ifirst,supersource,supergroup)

    CALL mp_sync(supergroup)
    CALL mp_bcast(i_meta,supersource,supergroup)

    CALL mp_sync(supergroup)
    CALL mp_bcast(i_cvst,supersource,supergroup)

    CALL mp_sync(supergroup)
    CALL mp_bcast(i_cvst_mw,SIZE(i_cvst_mw),supersource,supergroup)

    IF (.NOT.grandparent)CALL zeroing(taupw_scr)
    CALL mp_sync(supergroup)
    CALL mp_bcast(taupw_scr,SIZE(taupw_scr),supersource,supergroup)
    IF (paral%parent)THEN
       CALL dcopy(3*maxsys%nax*maxsys%nsx,taupw_scr(1,1,1,iwalk),1,fhills,1)
    ENDIF
    CALL mp_sync(supergroup)    
    ! WRITE RESTART
    DO iwalk1=1,mwi%nwalk
       IF ((MOD(i_meta,imeta%st_freq) .EQ. 0  .AND. i_cvst_mw(iwalk1).EQ.0)&
            .OR. soft_com%exsoft .OR. lmetares) lmetares = .TRUE.
    ENDDO
    ! 
    IF (lmetares.AND.i_meta.GT.1.AND.grandparent)THEN
       DO iwalk1=1,mwi%nwalk
          CALL mw_filename('MTD_RESTART_',fmtdres,iwalk1)
          CALL wmtdres_mw(50,ntot_iter,i_meta,iwalk1)
       ENDDO
    ENDIF

    CALL mp_sync(supergroup)
    CALL mp_bcast(lmetares,supersource,supergroup)

    ! QUENCH BO
    IF (.NOT.cntl%tmdbo.AND.MOD(i_meta,imeta%qw_freq) .EQ. 0 .AND. i_cvst .EQ. 0)&
         lquench = .TRUE.

    ! ==--------------------------------------------------------------==
    IF (FINAL) THEN
       IF (grandparent) THEN
          ! Deallocate
          DEALLOCATE(det_colvar_mw,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
          DEALLOCATE(hill_add_mw,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
          DEALLOCATE(f1_cv,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
          DEALLOCATE(f2_cv,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
          DEALLOCATE(cv_ist_scr,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
          DEALLOCATE(i_temp_mw,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
          DEALLOCATE(cvd_ist,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
          DEALLOCATE(cv_store,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
          DEALLOCATE(cv_last,STAT=ierr)
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
       ENDIF

       DEALLOCATE(cv_dyn,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(taupw_scr,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(i_cvst_mw,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    END IF
    DEALLOCATE(f_hill,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(f_wall,STAT=ierr)
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
    DEALLOCATE(taupw,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! Halt time for the routine
    CALL tihalt(procedureN,isub)

    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE meta_colvar_mw
  ! ==================================================================

END MODULE meta_colvar_utils
