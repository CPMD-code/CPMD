MODULE meta_cell_utils
  USE cnst,                            ONLY: au_kcm
  USE cnst_dyn,                        ONLY: &
       cscl_fac, cscl_val, det_ht, hllh_val, hllw_val, ht_ist, ht_name, &
       ht_path, imeta, inter_hill, inter_hill_max, l2dc, lfullc, lisoc, &
       lmeta, mdcellr, ncolvar, rmeta, tad_scf, toll_avcv, tvolbound
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_app,&
                                             fo_old,&
                                             fo_verb
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE latgen_utils,                    ONLY: omegagen
  USE meta_cv_utils,                   ONLY: rmsd_rot
  USE meta_hpot_utils,                 ONLY: hills,&
                                             hills_lor,&
                                             hills_ratio,&
                                             hills_sals,&
                                             hills_sals_shift
  USE metr,                            ONLY: metr_com
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sync
  USE odiis_utils,                     ONLY: solve
  USE parac,                           ONLY: parai,&
                                             paral
  USE prcp,                            ONLY: prcp_com,&
                                             prcpl
  USE readsr_utils,                    ONLY: readsi,&
                                             readsr,&
                                             xstring
  USE rmas,                            ONLY: rmass
  USE ropt,                            ONLY: infi,&
                                             iteropt
  USE soft,                            ONLY: soft_com
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             cntr,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE tpar,                            ONLY: dt_ions
  USE zeroing_utils,                   ONLY: zeroing
#include "sizeof.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: meta_cell
  PUBLIC :: meta_cell_inp
  PUBLIC :: delrot_cell

CONTAINS

  ! ==================================================================
  SUBROUTINE meta_cell(velp,lquench,lmetares)
    ! ==--------------------------------------------------------------==
    ! ==  Calculates the contributions to the forces on the cell      ==
    ! ==  due to the hills that have been accumulated in the space    ==
    ! ==  of the collective variables.                                ==
    ! ==  As collective variables we take functions of the cell       ==
    ! ==  In the simplest case they are directly the 6 parameters     ==
    ! ==  that determine the metric tensor.                           ==
    ! ==  in an isotropic case they are the cell side lengths         == 
    ! ==  This routine drives the metadynamics of the cell            ==
    ! ==--------------------------------------------------------------==

    ! NCOLVAR   : # of collective variables
    ! LHILLS    : true when hills are used in the dynamics of the constraints
    ! I_META    : # of steps in the meta-dynamics
    ! I_META_MAX: max # of steps in the meta-dynamics
    ! I_META_RES: if restart from previous meta dynamics, 
    ! # of steps already done
    ! I_HTST    : # of steps since last HILL has been added
    ! HT_IST    : actual values of collective variables, dim NCOLVAR
    ! HLLH      : hills altitude 
    ! HLLW      : hills amplitude 
    ! CSCL_FAC  : scale factors required for an homogeneous dynamic in the
    ! MCNSTR directions in the space of phases
    ! they can be read from input and can be tuned on the 
    ! amplitudes of the collective variables fluctuationso
    ! F_HILL     : force contribution on ions due to the accumulation 
    ! of hills  in the space of CV 
    ! dim (NODIM)
    ! META_RESTART: restart for a previous meta-dynamics
    ! the previous history of the run is read from file 
    ! further output will be appended 
    ! HT_PATH    :  CV history (dim NTOT_ITER*NCOLVAR)
    ! SCL_PATH   :  CV scaling factors history (dim NTOT_ITER*NCOLVAR)
    ! HLLW_VAL   :  hills amplitude (dim NTOT_ITER)
    ! HLLH_VAL   :  hills altitude  (dim NTOT_ITER)
    ! F_AVER     :  accumulator for the calculation of the average force
    ! acting on the CV, as a contribution of the 
    ! underlaying potential 
    ! IF_AVER    :  counter for the calculation of the average
    ! DISPLACEMENT: norm of the displacment in CV space
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: velp(:,:,:)
    LOGICAL                                  :: lquench, lmetares

    CHARACTER(*), PARAMETER                  :: procedureN = 'meta_cell'
    CHARACTER(len=10), PARAMETER             :: file1 = 'colvar_mtd', &
                                                file2 = 'parvar_mtd'

    CHARACTER(len=10)                        :: ch_temp
    INTEGER                                  :: i, ierr, iht, ii, iindex, ij, &
                                                isub, j, l_f_ht, l_scr_hil, &
                                                length, ntot_iter
    INTEGER, SAVE                            :: i_htst = 0, i_meta = 0, &
                                                i_temp = 0, if_aver = 0, &
                                                ifirst = 0
    LOGICAL                                  :: ferror, hill_add, lskip
    REAL(real_8) :: aa1(3), aa2(3), aa3(3), av_disp, cosa, cosb, cosc, diff, &
      disp1, disp2, displacement, f2, fact, fact1, fact2, hh, hh_max, hh_min, &
      hhh, hllw2, min_disp, tollm
    REAL(real_8), ALLOCATABLE                :: aux(:), hh_test(:)
    REAL(real_8), ALLOCATABLE, SAVE          :: f1_ht(:), f_aver(:), &
                                                f_free(:), fvb_der(:), &
                                                ht_last(:)

! Set time for the routine

    CALL tiset(' META_CELL',isub)

    IF (.NOT.paral%io_parent) GOTO 9999

    ! ==--------------------------------------------------------------==
    ALLOCATE(aux(ncolvar),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL zeroing(aux)!,ncolvar)
    ! ==--------------------------------------------------------------==
    ! Temporary Array

    l_f_ht    = ncolvar
    l_scr_hil = ncolvar*ncolvar+30*ncolvar


    length  = l_f_ht+l_scr_hil

    ALLOCATE(hh_test(ncolvar),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    ntot_iter = imeta%i_meta_max+imeta%i_meta_res
    tollm     = toll_avcv

    IF (ifirst .EQ. 0) THEN
       i_meta  = imeta%i_meta_res+1
       i_temp  = 0
       i_htst  = 0
       if_aver = 0
       ! Allocate Memory

       ! CALL MEMORY(IP_HT_LAST,NCOLVAR,'HT_LAST')
       ALLOCATE(ht_last(ncolvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(ht_last)!,ncolvar)
       ! CALL MEMORY(IP_F_AVER,NCOLVAR,'F_AVER')
       ALLOCATE(f_aver(ncolvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(f_aver)!,ncolvar)
       ! CALL MEMORY(IP_F_FREE,NCOLVAR,'F_FREE')
       ALLOCATE(f_free(ncolvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(f_free)!,ncolvar)
       ! CALL MEMORY(IP_F1_HT,NCOLVAR,'F1_HT')
       ALLOCATE(f1_ht(ncolvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(f1_ht)!,ncolvar)
       ! CALL MEMORY(IP_FVB_DER,NCOLVAR,'FVB_DER')
       ALLOCATE(fvb_der(ncolvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(fvb_der)!,ncolvar)

       ! CALL MEMORY(IP_HT_PATH,NCOLVAR*NTOT_ITER,'HT_PATH')
       ALLOCATE(ht_path(imeta%i_meta_max+imeta%i_meta_res,ncolvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(ht_path)!,ncolvar*ntot_iter)
       ! CALL MEMORY(IP_CSCL_VAL,NCOLVAR*NTOT_ITER,'CSCL_VAL')
       ALLOCATE(cscl_val(imeta%i_meta_res+imeta%i_meta_max,ncolvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(cscl_val)!,ncolvar*ntot_iter)
       ! CALL MEMORY(IP_HLLW_VAL,NTOT_ITER,'HLLH_VAL')
       ALLOCATE(hllw_val(ntot_iter,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(hllw_val)!,ntot_iter)
       ! CALL MEMORY(IP_HLLH_VAL,NTOT_ITER,'HLLH_VAL')
       ALLOCATE(hllh_val(ntot_iter,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(hllh_val)!,ntot_iter)

       ! If META_RESTART read data from cntl%proper input
       IF (lmeta%meta_restart) THEN
          ferror=.FALSE.
          IF (paral%io_parent)&
               CALL fileopen(51,file1,fo_old,ferror)
          IF (ferror) GOTO 20
          IF (paral%io_parent)&
               CALL fileopen(52,file2,fo_old,ferror)
          IF (ferror) GOTO 20
          IF (paral%io_parent)&
               REWIND(51)
          IF (paral%io_parent)&
               REWIND(52)
          DO ii = 1,imeta%i_meta_res
             IF (paral%io_parent)&
                  READ(51,err=20,END=20,fmt=*)  ij,&
                  (ht_path(ii,iht),iht=1,ncolvar),&
                  (cscl_val(ii,iht),iht=1,ncolvar)
             ! from output forces and parameters
             IF (paral%io_parent)&
                  READ(52,err=20,END=20,fmt=*) ij,&
                  hllw_val(ii,1),hllh_val(ii,1)
          ENDDO
          IF (paral%io_parent)&
               CALL fileclose(51)
          IF (paral%io_parent)&
               CALL fileclose(52)
          GOTO 100
20        CONTINUE
          IF (paral%io_parent)&
               WRITE(6,*) ' ERROR WHILE READING RESTART DATA '
          CALL stopgm('META_CELL',' ',& 
               __LINE__,__FILE__)
100       CONTINUE
       ENDIF ! META_RESTART

       IF (paral%io_parent)&
            WRITE(ch_temp,'(A5)') '     '
       DO i=1,6
          ht_name(i) = ch_temp
       ENDDO
       IF (lfullc) THEN
          IF (paral%io_parent)&
               WRITE(ch_temp,'(A5)')'SIDE1'
          ht_name(1) = ch_temp
          IF (paral%io_parent)&
               WRITE(ch_temp,'(A5)')'SIDE2'
          ht_name(2) = ch_temp
          IF (paral%io_parent)&
               WRITE(ch_temp,'(A5)')'SIDE3'
          ht_name(3) = ch_temp
          IF (paral%io_parent)&
               WRITE(ch_temp,'(A4)')'cosA'
          ht_name(4) = ch_temp
          IF (paral%io_parent)&
               WRITE(ch_temp,'(A4)')'cosB'
          ht_name(5) = ch_temp
          IF (paral%io_parent)&
               WRITE(ch_temp,'(A4)')'cosC'
          ht_name(6) = ch_temp
       ELSEIF (lisoc) THEN
          IF (paral%io_parent)&
               WRITE(ch_temp,'(A5)')'SIDE1'
          ht_name(1) = ch_temp
          IF (paral%io_parent)&
               WRITE(ch_temp,'(A5)')'SIDE2'
          ht_name(2) = ch_temp
          IF (paral%io_parent)&
               WRITE(ch_temp,'(A5)')'SIDE3'
          ht_name(3) = ch_temp
       ENDIF
    ENDIF

    ! ==--------------------------------------------------------------==
    ! Calculate the Collective Variables as functions of the cell parameters 
    CALL zeroing(det_ht)!,9*ncolvar)
    IF (lfullc) THEN
       CALL parfull(metr_com%ht,ht_ist,det_ht)
    ELSEIF (lisoc) THEN
       CALL pariso(metr_com%ht,ht_ist,det_ht)
    ELSEIF (l2dc) THEN
       DO iht = 1,2
       ENDDO
    ENDIF

    ! ==--------------------------------------------------------------==
    ! Store CV values and Check displacement

    IF (ifirst .EQ. 0) THEN
       DO iht = 1,ncolvar
          ht_last(iht) = ht_ist(iht)
       ENDDO
       ifirst = 1
       GOTO 200
    ENDIF

    hill_add = .FALSE.
    IF (lmeta%ttrajovrr) cnti%ntraj=999999999

    i_htst = i_htst + 1

    IF (i_htst .EQ.  inter_hill ) THEN
       min_disp = 100.0_real_8
       ii = inter_hill
       av_disp = 0.0_real_8
       DO iht = 1,ncolvar
          diff = ((ht_ist(iht)-ht_last(iht))/cscl_fac(1,iht))**2.0_real_8
          av_disp = av_disp + diff
          min_disp = MIN(diff,min_disp)
       ENDDO
       av_disp = SQRT(av_disp)
       IF (av_disp .GT. tollm .OR. i_meta .EQ. 1) hill_add = .TRUE.

       IF (paral%io_parent)&
            WRITE(6,'(/10x,A,f10.6,A,f10.6,A,f10.6/)')&
            'MIN DISP. = ',SQRT(min_disp),'  ||DISP.tot|| = ', av_disp,&
            '  Tolerance ',tollm
    ELSEIF (i_htst .GT.  inter_hill ) THEN
       i_temp = i_temp + 1
       ii = inter_hill
       min_disp = 100.0_real_8
       DO iht = 1,ncolvar
          diff = ((ht_ist(iht)-ht_last(iht))/cscl_fac(1,iht))**2
          aux(iht) = diff
          min_disp = MIN(diff,min_disp)
       ENDDO
       IF (MOD(i_temp,imeta%icheck) .EQ. 0 ) THEN
          av_disp = 0.0_real_8
          DO iht = 1,ncolvar
             av_disp = av_disp+aux(iht)
          ENDDO
          av_disp = SQRT(av_disp)
          IF (paral%io_parent)&
               WRITE(6,'(/10x,A,f10.6,A,f10.6,A,f10.6/)')&
               'MIN DISP. = ',SQRT(min_disp),'  ||DISP.tot|| = ',av_disp,&
               '  Tolerance ',tollm
          IF (av_disp .GT. tollm .OR. i_htst .EQ. inter_hill_max) THEN
             hill_add = .TRUE.
          ENDIF
       ENDIF
    ENDIF

    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    ! Calculation of the global forces on the collective variables
    ! in order to determine the optimal hills height that guarantees 
    ! the maximum balance between the opposite contribution 

    CALL ht_forces(hh_test,velp)
    !$omp parallel do private(IHT)
    DO iht = 1,ncolvar
       f_aver(iht) = f_aver(iht) + hh_test(iht)&
            + f1_ht(iht)/cscl_fac(1,iht)
       f_free(iht) = f_free(iht) + hh_test(iht)
    ENDDO
    if_aver = if_aver + 1

    ! ==--------------------------------------------------------------==
    ! New Step of Meta-Dynamics in th Space of CV
    IF (hill_add) THEN
       IF (MOD(i_meta,imeta%tr_freq) .EQ. 0 .OR. i_meta .EQ. 1 ) THEN
          IF (lmeta%ttrajovrr) cnti%ntraj=1
       ENDIF

       ! ==--------------------------------------------------------------==
       ! Squared Norm of the Displacement in the CV space

       displacement = 0.0_real_8
       disp1=0.0_real_8
       DO iht = 1,ncolvar
          displacement = displacement +&
               (ht_ist(iht)-ht_last(iht))*(ht_ist(iht)-ht_last(iht))/&
               (cscl_fac(1,iht)*cscl_fac(1,iht))
       ENDDO
       disp2 =  displacement
       displacement = SQRT(displacement)
       ! write(6,*) 'dis ' , DISPLACEMENT
       ! Average forces on the collective variables 
       !$omp parallel do private(IHT)
       DO iht = 1,ncolvar
          f_aver(iht)  = f_aver(iht)* cscl_fac(1,iht)/REAL(if_aver,kind=real_8)
          f_free(iht)  = f_free(iht)* cscl_fac(1,iht)/REAL(if_aver,kind=real_8)
       ENDDO
       if_aver = 0
       ! ==--------------------------------------------------------------==
       ! Tuning Hills parameters 
       IF (lmeta%ltune_hh) THEN
          hh_min = 100._real_8
          hh_max = 0.0_real_8
          lskip = .FALSE.

          DO iht = 1,ncolvar
             hllw2 = rmeta%hllw*rmeta%hllw
             ! Fact = (HT_IST(IHT)-HT_LAST(IHT))*
             ! &           (1.0_real_8/(HLLW*HLLW)+1.0_real_8/(DISPLACEMENT*DISPLACEMENT))

             fact = cscl_fac(1,iht)/(ht_ist(iht)-ht_last(iht))*&
                  EXP(0.5_real_8)*hllw2*disp2/&
                  (disp2*EXP(-0.5_real_8*disp2/hllw2)+&
                  hllw2*(EXP(-0.5_real_8*disp2/hllw2)-&
                  EXP(-0.5_real_8*rmeta%rshift*rmeta%rshift)))

             ! HHH     = 0.5_real_8*abs(F_AVER(IHT)*exp(0.5_real_8)/FACT)
             hhh     = 0.1_real_8*ABS(f_aver(iht)*fact)
             hh_min       = MIN(hh_min,hhh)
             hh_max       = MAX(hh_max,hhh)
          ENDDO

          hh       = (hh_min + hh_max)*0.5_real_8
          IF (hh_min .GT. rmeta%htop .OR. hh_max .LT. rmeta%hlow) THEN
             lskip = .TRUE.
             IF (paral%io_parent)&
                  WRITE(6,'(A,f12.6,A,f12.6,A)') 'Hill top > ', rmeta%htop*au_kcm,&
                  ' Kcal, or Hill top < ',rmeta%hlow*au_kcm,' Kcal',&
                  'keep the prev. value'
          ENDIF
          IF (.NOT. lskip .AND. hh .LT. rmeta%hlow) THEN
             rmeta%hllh = (rmeta%hlow+hh_max)*0.5_real_8
          ELSEIF (.NOT. lskip .AND. hh .GT. rmeta%htop) THEN
             rmeta%hllh = (hh_min+rmeta%htop)*0.5_real_8
          ELSEIF (.NOT. lskip ) THEN
             rmeta%hllh = hh
          ELSE
             rmeta%hllh = rmeta%hllh
          ENDIF
          IF (paral%io_parent)&
               WRITE(6,'(/A,f12.6,3(3x,A,f12.6)/)')&
               '||DISP|| = ',displacement,&
               'HH tuned min = ', hh_min,'HH tuned max = ', hh_max,&
               ' HH used ',rmeta%hllh
       ENDIF

       ! ==--------------------------------------------------------------==
       ! Update accumulators
       !$omp parallel do private(IHT)
       DO iht  = 1,ncolvar
          ht_path(i_meta,iht) =  ht_last(iht)
          cscl_val(i_meta,iht)=  cscl_fac(1,iht)
          hllw_val(i_meta,1)    =  rmeta%hllw
          hllh_val(i_meta,1)    =  rmeta%hllh
       ENDDO
       ! ==--------------------------------------------------------------==

       ! Open output file
       IF (paral%parent)&
            CALL ht_write_out(i_meta,ht_last,f1_ht,f_aver,f_free,fvb_der,&
            displacement)
       CALL zeroing(f_aver)!,ncolvar)
       CALL zeroing(f_free)!,ncolvar)
       ! ==--------------------------------------------------------------==
       !$omp parallel do private(IHT)
       DO iht = 1,ncolvar
          ht_last(iht)   = ht_ist(iht)
       ENDDO

       i_meta = i_meta + 1
       i_temp = 0
       i_htst = 0
       ! stop
    ENDIF        ! HILL_ADD



    ! ==--------------------------------------------------------------==
    ! Calculation of the contributions to the forces on ions
    ! due to HILLS in the space of CV

    ! Calculate hills contribution to forces

200 CONTINUE
    IF (lmeta%lhills) THEN
       iindex = i_meta
       ! call DCOPY(9,ht,1,ht_tmp,1)
       ! eeps = 0.00001_real_8
       ! do i = 1,3
       ! do j = 1,3
       ! ht(i,j) = ht_tmp(i,j) + eeps
       ! CALL PARFULL(HT,HT_IST,DET_HT)
       ! CALL HILLS_SALS_SHIFT(HT_IST,NCOLVAR,IINDEX,NTOT_ITER,F1_HT,
       ! &         HT_LAST,AUX(ISC2),LSCR2,I_HTST)
       ! W1=GAUSSPOT
       ! ht(i,j) = ht_tmp(i,j) - eeps
       ! CALL PARFULL(HT,HT_IST,DET_HT)
       ! CALL HILLS_SALS_SHIFT(HT_IST,NCOLVAR,IINDEX,NTOT_ITER,F1_HT,
       ! &         HT_LAST,AUX(ISC2),LSCR2,I_HTST)
       ! W2=GAUSSPOT
       ! HT_FOR(i,j) = -(W1-W2)/eeps/2.0_real_8
       ! ht(i,j) = ht_tmp(i,j)
       ! c            write(6,'(2I4,2f12.8)') i,j,W1,W2
       ! enddo
       ! enddo

       ! CALL PARFULL(HT,HT_IST,DET_HT)

       CALL zeroing(f1_ht)!,ncolvar)

       IF (lmeta%hlore) THEN
          CALL hills_lor(ht_ist,ncolvar,iindex,ntot_iter,f1_ht,&
               1,0)
       ELSEIF (lmeta%hratio) THEN
          CALL hills_ratio(ht_ist,ncolvar,iindex,ntot_iter,f1_ht,&
               ht_last,i_htst,1,0)
       ELSEIF (lmeta%hshift) THEN
          CALL hills_sals_shift(ht_ist,ncolvar,iindex,ntot_iter,f1_ht,&
               ht_last,i_htst,1,0)
       ELSE IF (lmeta%sphere) THEN
          CALL hills(ht_ist,ncolvar,i_meta,ntot_iter,f1_ht,&
               1,0)
       ELSE
          CALL hills_sals(ht_ist,ncolvar,iindex,ntot_iter,f1_ht,&
               ht_last,i_htst,1,0)
       ENDIF


    ENDIF

    CALL zeroing(fvb_der)!,ncolvar)
    IF (tvolbound) THEN
       aa1(1) = metr_com%ht(1,1)
       aa1(2) = metr_com%ht(1,2)
       aa1(3) = metr_com%ht(1,3)
       aa2(1) = metr_com%ht(2,1)
       aa2(2) = metr_com%ht(2,2)
       aa2(3) = metr_com%ht(2,3)
       aa3(1) = metr_com%ht(3,1)
       aa3(2) = metr_com%ht(3,2)
       aa3(3) = metr_com%ht(3,3)
       CALL  omegagen(aa1,aa2,aa3,parm%omega)

       f2 = 0.0_real_8
       IF (parm%omega .LT. mdcellr%volmin+mdcellr%volmin/100._real_8) THEN
          f2 = +4.0_real_8*mdcellr%vbarrier*(mdcellr%volmin*mdcellr%volmin*mdcellr%volmin*mdcellr%volmin)/&
               (parm%omega**(5.0_real_8))
          IF (paral%io_parent)&
               WRITE(6,'(72A)') ('*',i = 1,72)
          IF (paral%io_parent)&
               WRITE(6,'(10x,A,f12.4,A,f12.4)')&
               'WARNING: Volume = ',parm%omega,' <= Volmin = ',mdcellr%volmin
          IF (paral%io_parent)&
               WRITE(6,'(72A)')  ('*',i = 1,72)
          IF (paral%io_parent)&
               WRITE(6,'(10x,A,f14.6)') 'Force: ',f2
       ELSEIF (parm%omega .GT. mdcellr%volmax-mdcellr%volmax/100._real_8) THEN
          f2 = -4.0_real_8*mdcellr%vbarrier/(mdcellr%volmax*mdcellr%volmax*mdcellr%volmax*mdcellr%volmax)*&
               (parm%omega*parm%omega*parm%omega)
          IF (paral%io_parent)&
               WRITE(6,'(72A)') ('*',i = 1,72)
          IF (paral%io_parent)&
               WRITE(6,'(10x,A,f12.4,A,f12.4)')&
               'WARNING: Volume = ',parm%omega,' >= Volmax = ',mdcellr%volmax
          IF (paral%io_parent)&
               WRITE(6,'(72A)')  ('*',i = 1,72)
          IF (paral%io_parent)&
               WRITE(6,'(10x,A,2f14.6)') 'Force: ',f2,mdcellr%vbarrier
       ENDIF
       IF (lfullc) THEN
          fact1 = 1.0_real_8+2.0_real_8*ht_ist(4)*ht_ist(5)*ht_ist(6)-&
               ht_ist(4)*ht_ist(4)-ht_ist(5)*ht_ist(5)-&
               ht_ist(6)*ht_ist(6)
          fact1 = SQRT(fact1)

          fvb_der(1) = f2 *ht_ist(2)* ht_ist(3)*fact1
          fvb_der(2) = f2 *ht_ist(1)* ht_ist(3)*fact1
          fvb_der(3) = f2 *ht_ist(1)* ht_ist(2)*fact1

          fact2 = (ht_ist(1) * ht_ist(2)* ht_ist(3))/fact1

          fvb_der(4) = f2 * fact2*(ht_ist(5)*ht_ist(6)-ht_ist(4))
          fvb_der(5) = f2 * fact2*(ht_ist(4)*ht_ist(6)-ht_ist(5))
          fvb_der(6) = f2 * fact2*(ht_ist(4)*ht_ist(5)-ht_ist(6))

       ELSEIF (lisoc) THEN
          cosa = aa2(1)*aa3(1)+aa2(2)*aa3(2)+aa2(3)*aa3(3)
          cosa = cosa  / (ht_ist(2)* ht_ist(3))
          cosb = aa1(1)*aa3(1)+aa1(2)*aa3(2)+aa1(3)*aa3(3)
          cosb = cosb  / (ht_ist(1)* ht_ist(3))
          cosc = aa1(1)*aa2(1)+aa1(2)*aa2(2)+aa1(3)*aa2(3)
          cosc = cosc  / (ht_ist(1)* ht_ist(2))

          fact1 = 1.0_real_8+2.0_real_8*cosa*cosb*cosc-cosa*cosa-&
               cosb*cosb-cosc*cosc
          fact1 = SQRT(fact1)

          fvb_der(1) = f2 *ht_ist(2)* ht_ist(3)*fact1
          fvb_der(2) = f2 *ht_ist(1)* ht_ist(3)*fact1
          fvb_der(3) = f2 *ht_ist(1)* ht_ist(2)*fact1

       ELSEIF (l2dc) THEN
          DO iht = 1,2
          ENDDO
       ENDIF

    ENDIF

    CALL zeroing(mdcellr%fh_cell)!,9)
    DO iht = 1,ncolvar
       DO i = 1,3
          DO j = 1,3
             mdcellr%fh_cell(j,i) = mdcellr%fh_cell(j,i) +&
                  (f1_ht(iht)/cscl_fac(1,iht)+fvb_der(iht))*&
                  det_ht(j+(i-1)*3,iht)
          ENDDO
       ENDDO
    ENDDO





    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    ! IF MAX_STEP reached or EXIT file exist

    IF (infi.EQ.cnti%nomore .OR. i_meta-1.EQ.ntot_iter) THEN
       soft_com%exsoft = .TRUE.

       IF (paral%io_parent)&
            WRITE(6,'(A)')&
            ' <<<<<<<<  CELL COLLECTIVE VARIABLE INFO   >>>>>>>>>'
       IF (paral%io_parent)&
            WRITE(6,'(4x,A,6x,A)') 'TYPE','VALUE'

       DO iht = 1,ncolvar
          IF (paral%io_parent)&
               WRITE(6,'(2X,A10,2X,f10.4)') ht_name(iht),ht_ist(iht)
       ENDDO

       ! DEALLOCATE
       DEALLOCATE(ht_last,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(ht_path,STAT=ierr)
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
    ! ==--------------------------------------------------------------==
    DEALLOCATE(hh_test,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
9999 CONTINUE
    CALL mp_sync(parai%allgrp)

    CALL mp_bcast_byte(soft_com, size_in_bytes_of(soft_com),parai%io_source,parai%cp_grp)
    CALL mp_bcast(i_meta,parai%io_source,parai%cp_grp)
    CALL mp_bcast(i_htst,parai%io_source,parai%cp_grp)

    ! WRITE RESTART
    IF ((MOD(i_meta,imeta%st_freq) .EQ. 0 .AND. i_htst .EQ. 0)&
         .OR. soft_com%exsoft) THEN
       lmetares = .TRUE.
    ENDIF

    ! QUENCH BO
    IF (MOD(i_meta,imeta%qw_freq) .EQ. 0 .AND. i_htst .EQ. 0)&
         lquench = .TRUE.

    ! ==--------------------------------------------------------------==
    ! Halt time for the routine
    CALL tihalt(' META_CELL',isub)

    RETURN
  END SUBROUTINE meta_cell
  ! ==================================================================

  SUBROUTINE ht_write_out(i_meta,ht_last,f1_ht,f_aver,f_free,&
       fvb_der,displacement)


    INTEGER                                  :: i_meta
    REAL(real_8) :: ht_last(ncolvar), f1_ht(ncolvar), f_aver(ncolvar), &
      f_free(ncolvar), fvb_der(ncolvar), displacement

    CHARACTER(len=10), PARAMETER :: file1 = 'colvar_mtd', &
      file2 = 'parvar_mtd', file3 = 'forfac_mtd', file4 = 'disvar_mtd', &
      file5 = 'enevar_mtd'

    CHARACTER(len=10)                        :: chnum
    CHARACTER(len=100)                       :: lineform
    EXTERNAL                                 :: ddot
    INTEGER                                  :: iaa, icv, iee
    INTEGER, SAVE                            :: ifirst = fo_verb
    LOGICAL                                  :: ferror
    REAL(real_8)                             :: ddot, ekinh

! ==--------------------------------------------------------------==
! Open output files

    IF (paral%io_parent) THEN
       CALL fileopen(51,file1,fo_app+ifirst,ferror)
       CALL fileopen(52,file2,fo_app+ifirst,ferror)
       CALL fileopen(53,file3,fo_app+ifirst,ferror)
       CALL fileopen(54,file4,fo_app+ifirst,ferror)
       CALL fileopen(55,file5,fo_app+ifirst,ferror)
    ENDIF
    ! Verbose open only on the first run
    ifirst=0
    ! ==--------------------------------------------------------------==

    ! Write output
    IF (paral%io_parent)&
         WRITE(chnum,'(I5)') ncolvar
    CALL xstring(chnum,iaa,iee)
    lineform =&
         '(1X,I7,'//chnum(iaa:iee)//'E14.6,'//chnum(iaa:iee)//'F10.6)'
    IF (paral%io_parent)&
         WRITE(51,lineform) iteropt%nfi,(ht_last(icv),icv=1,ncolvar),&
         (cscl_fac(1,icv),icv=1,ncolvar)

    lineform =  '(1X,I7,3F14.6)'
    IF (paral%io_parent)&
         WRITE(52,lineform) iteropt%nfi,displacement,rmeta%hllw,rmeta%hllh
    ! Close files
    IF (paral%io_parent)&
         CALL fileclose(51)
    IF (paral%io_parent)&
         CALL fileclose(52)

    ! Print force factors
    lineform =  '(1X,I7,'//chnum(iaa:iee)//'E14.6,'&
         //chnum(iaa:iee)//'E14.6,'//chnum(iaa:iee)//'E16.4)'
    IF (paral%io_parent)&
         WRITE(53,lineform) iteropt%nfi,(f_free(icv),icv=1,ncolvar),&
         (f1_ht(icv),icv=1,ncolvar),&
         (fvb_der(icv),icv=1,ncolvar)
    ! *                          (F_AVER(ICV),ICV=1,NCOLVAR)

    ! Close files
    IF (paral%io_parent)&
         CALL fileclose(53)

    ! Print last displacements
    lineform =  '(I6,'//chnum(iaa:iee)//'f11.6)'
    IF (paral%io_parent)&
         WRITE(54,lineform) iteropt%nfi,&
         ((ht_ist(icv)-ht_last(icv)),icv=1,ncolvar)

    ! Close files
    IF (paral%io_parent)&
         CALL fileclose(54)

    ! Print temperature of the cell and kinetic energy
    ! EHAM_HILL=EKINP+ETOT+ENOSE+ENOSP+ECNSTR+EKINC+EKINH
    ekinh=0.5_real_8*cntr%cmass*ddot(9,metr_com%htvel,1,metr_com%htvel,1)
    lineform = '(I6,5f16.8)'
    IF (paral%io_parent)&
         WRITE(55,lineform) iteropt%nfi,&
         ekinh,rmeta%gausspot,&
         ener_com%etot,rmeta%eham_hill,rmeta%eham_hill+rmeta%gausspot
    IF (paral%io_parent)&
         CALL fileclose(55)

    ! ==--------------------------------------------------------------==

    RETURN
  END SUBROUTINE ht_write_out
  ! ==--------------------------------------------------------------==

  SUBROUTINE ht_forces(cv_f_new,velp)
    REAL(real_8)                             :: cv_f_new(ncolvar), velp(:,:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'ht_forces'
    INTEGER, PARAMETER                       :: maxrat = 200, mrdiis = 5 
    REAL(real_8), PARAMETER                  :: tol_cv_f = 5.e-10_real_8, &
                                                tol_cv_x = 5.e-8_real_8 

    EXTERNAL                                 :: dasum, ddot
    INTEGER                                  :: i, ia, idiis, ierr, iht, ii, &
                                                info, iter, j, jj, k, l_ht, &
                                                length, m
    INTEGER, ALLOCATABLE                     :: ipvt(:)
    INTEGER, SAVE                            :: ifirst = 0
    REAL(real_8) :: dasum, ddot, diism(mrdiis+1,mrdiis+1), dx(9), errf, errx, &
      fact, fo, fomega, hscr(3,3), ht_fm(3,3), ht_for(3,3), ht_new(3,3), &
      ht_nnew(3,3), ht_old(3,3), ht_tmp(3,3), hv_new(3,3), hv_old(3,3), &
      v(mrdiis+1)
    REAL(real_8), ALLOCATABLE :: asl(:,:), det_ht_n(:,:), det_ht_nn(:,:), &
      det_ht_o(:,:), err(:,:), ht_diff(:), ht_ist_n(:), ht_ist_nn(:), &
      ht_ist_o(:), xlo(:,:)
    REAL(real_8), ALLOCATABLE, SAVE          :: cv_f_old(:)

! ==--------------------------------------------------------------==

    length = (32+2*mrdiis)*ncolvar+ncolvar*ncolvar
    l_ht = ncolvar
    ! TODO align for BG
    ALLOCATE(ht_ist_o(ncolvar),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(ht_ist_n(ncolvar),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(ht_ist_nn(ncolvar),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(det_ht_o(9, ncolvar),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(det_ht_n(9, ncolvar),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(det_ht_nn(9, ncolvar),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    ALLOCATE(ht_diff(ncolvar),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    ALLOCATE(asl(ncolvar, ncolvar),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(xlo(ncolvar, mrdiis),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(err(ncolvar, mrdiis),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(ipvt(ncolvar),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)


    ! Memory Allocation at IFIRST=0
    IF (ifirst .EQ. 0) THEN
       ! CALL MEMORY(IP_CV_F_OLD,NCOLVAR,'CV_F_OLD')
       ALLOCATE(cv_f_old(ncolvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ifirst = 1
       CALL zeroing(cv_f_old)!,ncolvar)
    ENDIF
    CALL dcopy(9,metr_com%htfor(1,1),1,ht_tmp(1,1),1)
    CALL zeroing(hscr)!,9)
    DO i=1,ions1%nsp
       fact=rmass%pma(i)
       DO ia=1,ions0%na(i)
          hscr(1,1)=hscr(1,1)+fact*velp(1,ia,i)*velp(1,ia,i)
          hscr(2,1)=hscr(2,1)+fact*velp(2,ia,i)*velp(1,ia,i)
          hscr(3,1)=hscr(3,1)+fact*velp(3,ia,i)*velp(1,ia,i)
          hscr(1,2)=hscr(1,2)+fact*velp(1,ia,i)*velp(2,ia,i)
          hscr(2,2)=hscr(2,2)+fact*velp(2,ia,i)*velp(2,ia,i)
          hscr(3,2)=hscr(3,2)+fact*velp(3,ia,i)*velp(2,ia,i)
          hscr(1,3)=hscr(1,3)+fact*velp(1,ia,i)*velp(3,ia,i)
          hscr(2,3)=hscr(2,3)+fact*velp(2,ia,i)*velp(3,ia,i)
          hscr(3,3)=hscr(3,3)+fact*velp(3,ia,i)*velp(3,ia,i)
       ENDDO
    ENDDO

    CALL dgemm('N','N',3,3,3,1._real_8,hscr(1,1),3,metr_com%htm1(1,1),3,1._real_8,&
         ht_tmp(1,1),3)

    CALL dgemm('N','N',3,3,3,1.0_real_8,ht_fm(1,1),3,metr_com%htm1(1,1),3,0.0_real_8,&
         ht_for(1,1),3)


    fact=dt_ions/(2.0_real_8*cntr%cmass)
    IF (prcpl%tzflex) THEN
       ! z-direction flexible cell
       fomega=0.0_real_8
       fo=parm%omega**(2._real_8/3._real_8) * prcp_com%omega0**(1._real_8/3._real_8)
       fomega=ddot(9,ht_tmp,1,metr_com%ht0,1)/fo
       ht_old(1,3)=metr_com%ht(1,3)-dt_ions*(fact*fomega*metr_com%ht0(3,1)/fo)
       ht_old(2,3)=metr_com%ht(2,3)-dt_ions*(fact*fomega*metr_com%ht0(3,2)/fo)
       ht_old(3,3)=metr_com%ht(3,3)-dt_ions*(fact*fomega*metr_com%ht0(3,3)/fo)
    ELSEIF (prcpl%tisot) THEN
       ! Isotropic cell
       fomega=0.0_real_8
       fo=3._real_8 * parm%omega**(2._real_8/3._real_8) * prcp_com%omega0**(1._real_8/3._real_8)
       fomega=ddot(9,ht_tmp,1,metr_com%ht0,1)/fo
       DO k=1,3
          ht_old(1,k)=metr_com%ht(1,k)-dt_ions*(fact*fomega*metr_com%ht0(k,1)/fo)
          ht_old(2,k)=metr_com%ht(2,k)-dt_ions*(fact*fomega*metr_com%ht0(k,2)/fo)
          ht_old(3,k)=metr_com%ht(3,k)-dt_ions*(fact*fomega*metr_com%ht0(k,3)/fo)
       ENDDO
    ELSE
       ! General
       DO k=1,3
          hv_old(1,k)=metr_com%htvel(1,k) - fact*(ht_tmp(1,k)-mdcellr%fh_cell(1,k))
          hv_old(2,k)=metr_com%htvel(2,k) - fact*(ht_tmp(2,k)-mdcellr%fh_cell(2,k))
          hv_old(3,k)=metr_com%htvel(3,k) - fact*(ht_tmp(3,k)-mdcellr%fh_cell(3,k))
          ht_old(1,k)=metr_com%ht(1,k)-dt_ions*hv_old(1,k)
          ht_old(2,k)=metr_com%ht(2,k)-dt_ions*hv_old(2,k)
          ht_old(3,k)=metr_com%ht(3,k)-dt_ions*hv_old(3,k)
       ENDDO
    ENDIF
    ! 
    IF (lfullc) THEN
       CALL parfull(ht_old,ht_ist_o,det_ht_o)
    ELSEIF (lisoc) THEN
       CALL pariso(ht_old,ht_ist_o,det_ht_o)
    ELSEIF (l2dc) THEN
       DO iht = 1,2
       ENDDO
    ENDIF


    CALL dcopy(9,metr_com%ht,1,ht_new,1)
    CALL dcopy(9,metr_com%htvel,1,hv_new,1)

    IF (lfullc) THEN
       CALL parfull(ht_new,ht_ist_n,det_ht_n)
    ELSEIF (lisoc) THEN
       CALL pariso(ht_new,ht_ist_n,det_ht_n)
    ELSEIF (l2dc) THEN
       DO iht = 1,2
       ENDDO
    ENDIF

    ! First Guess for the forces
    CALL zeroing(dx)!,9)
    DO iht=1,ncolvar
       DO i=1,3
          DO j=1,3
             dx(j+(i-1)*3)=dx(j+(i-1)*3)+&
                  cv_f_old(iht)*det_ht_n(j+(i-1)*3,iht)
          ENDDO
       ENDDO
    ENDDO

    DO i = 1,3
       DO j = 1,3
          ht_for(j,i) = dx(j+(i-1)*3)
       ENDDO
    ENDDO

    fact=dt_ions/(2.0_real_8*cntr%cmass)
    IF (prcpl%tzflex) THEN
       ! z-direction flexible cell
       fo=parm%omega**(2._real_8/3._real_8) * prcp_com%omega0**(1._real_8/3._real_8)
       fomega=ddot(9,ht_for,1,metr_com%ht0,1)/fo
       ht_nnew(1,3)=ht_new(1,3)-fact*fomega*metr_com%ht0(1,3)/fo
       ht_nnew(2,3)=ht_new(2,3)-fact*fomega*metr_com%ht0(2,3)/fo
       ht_nnew(3,3)=ht_new(3,3)-fact*fomega*metr_com%ht0(3,3)/fo
    ELSE IF (prcpl%tisot) THEN
       ! Isotropic cell
       fo=3._real_8 * parm%omega**(2._real_8/3._real_8) * prcp_com%omega0**(1._real_8/3._real_8)
       fomega=ddot(9,ht_for,1,metr_com%ht0,1)/fo
       DO k=1,3
          ht_nnew(1,k)=ht_new(1,k)-fact*fomega*metr_com%ht0(1,k)/fo
          ht_nnew(2,k)=ht_new(2,k)-fact*fomega*metr_com%ht0(2,k)/fo
          ht_nnew(3,k)=ht_new(3,k)-fact*fomega*metr_com%ht0(3,k)/fo
       ENDDO
    ELSE
       ! General
       DO i=1,3
          ht_nnew(i,1)=ht_new(i,1) - dt_ions*hv_new(i,1)&
               + dt_ions * fact*ht_for(i,1) ! HV_NNEW(I,K)
          ht_nnew(i,2)=ht_new(i,2) - dt_ions*hv_new(i,2)&
               + dt_ions * fact*ht_for(i,2) ! HV_NNEW(I,K)
          ht_nnew(i,3)=ht_new(i,3) - dt_ions*hv_new(i,3)&
               + dt_ions * fact*ht_for(i,3) ! HV_NNEW(I,K)
       ENDDO
    ENDIF

    CALL dcopy(ncolvar,cv_f_old,1,cv_f_new,1)

    ! Iterativ calculation of lambda
    DO iter=1,maxrat
       ! Calculate Collective Variables values and differences
       ! for the current value of lambda
       CALL zeroing(dx)!,9)

       IF (lfullc) THEN
          CALL parfull(ht_nnew,ht_ist_nn,det_ht_nn)
       ELSEIF (lisoc) THEN
          CALL pariso(ht_nnew,ht_ist_nn,det_ht_nn)
       ELSEIF (l2dc) THEN
          DO iht = 1,2
          ENDDO
       ENDIF

       DO i = 1,ncolvar
          ht_diff(i) =  -(ht_ist_nn(i)-ht_ist_o(i))
       ENDDO

       errf=dasum(ncolvar,ht_diff(1),1)
       ! write(6,'(i5,f12.8)') ITER,ERRF 
       IF (errf .LT. tol_cv_f) GOTO 100

       ! Derivatives of sigma wrt lambda
       fact = dt_ions*dt_ions/(2.0_real_8*cntr%cmass)
       DO i=1,ncolvar
          DO j=1,ncolvar
             asl(i,j)=0.0_real_8
             DO ii=1,3
                DO jj = 1,3
                   asl(i,j)=asl(i,j)-fact*&
                        det_ht_nn(jj+(ii-1)*3,i)*det_ht_n(jj+(ii-1)*3,j)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       ! Solve for asl*cg=fc
       ! LAPACK matrix solver
       CALL dgesv(ncolvar,1,asl,ncolvar,ipvt,ht_diff,ncolvar,info)

       IF (info.NE.0) CALL stopgm(' HT_FORCES|','ERROR IN DGESV',& 
            __LINE__,__FILE__)
       errx=dasum(ncolvar,ht_diff(1),1)
       ! cntl%diis!
       idiis=MOD(iter-1,mrdiis)+1
       !$omp parallel do private(IHT)
       DO i = 1,ncolvar
          xlo(i,idiis) = cv_f_new(i)
       ENDDO

       CALL dcopy(ncolvar,ht_diff(1),1,err(1,idiis),1)
       IF (iter.GT.mrdiis) THEN
          m=mrdiis+1
          CALL zeroing(diism)!,m*m)
          CALL zeroing(v)!,m)
          DO i=1,mrdiis
             DO j=1,mrdiis
                diism(i,j)=ddot(ncolvar,err(1,i),1,err(1,j),1)
             ENDDO
             diism(m,i)=1.0_real_8
             diism(i,m)=1.0_real_8
          ENDDO

          v(m)=1.0_real_8
          CALL solve(diism,m,m,v)

          CALL zeroing(ht_diff)!,ncolvar)
          CALL zeroing(cv_f_new)!,ncolvar)
          DO i=1,mrdiis
             DO j=1,ncolvar
                ht_diff(j)  = ht_diff(j) + v(i)*err(j,i)
                cv_f_new(j) = cv_f_new(j)+ v(i)*xlo(j,i)
             ENDDO
          ENDDO
       ENDIF
       !ocl NOALIAS
#ifdef __NEC
       !CDIR NODEP
#endif
       !$omp parallel do private(I)
       DO  i = 1,ncolvar
          cv_f_new(i)   = cv_f_new(i)+ ht_diff(i)
       ENDDO

       IF (errx .LT. tol_cv_x) GOTO 100
       ! Update forces
       CALL zeroing(dx)!,6)
       DO iht=1,ncolvar
          DO i=1,3
             DO j = 1,3
                dx(j+(i-1)*3)=dx(j+(i-1)*3)&
                     +cv_f_new(iht)*det_ht_n(j+(i-1)*3,iht)
             ENDDO
          ENDDO
       ENDDO

       DO i = 1,3
          DO j = 1,3
             ht_for(j,i) = dx(j+(i-1)*3)
          ENDDO
       ENDDO

       fact=dt_ions/(2.0_real_8*cntr%cmass)
       IF (prcpl%tzflex) THEN
          ! z-direction flexible cell
          fo=parm%omega**(2._real_8/3._real_8) * prcp_com%omega0**(1._real_8/3._real_8)
          fomega=ddot(9,ht_for,1,metr_com%ht0,1)/fo
          ht_nnew(1,3)=ht_new(1,3)-fact*fomega*metr_com%ht0(3,1)/fo
          ht_nnew(2,3)=ht_new(2,3)-fact*fomega*metr_com%ht0(3,2)/fo
          ht_nnew(3,3)=ht_new(3,3)-fact*fomega*metr_com%ht0(3,3)/fo
       ELSE IF (prcpl%tisot) THEN
          ! Isotropic cell
          fo=3._real_8 * parm%omega**(2._real_8/3._real_8) * prcp_com%omega0**(1._real_8/3._real_8)
          fomega=ddot(9,ht_for,1,metr_com%ht0,1)/fo
          DO k=1,3
             ht_nnew(1,k)=ht_new(1,k)-fact*fomega*metr_com%ht0(k,1)/fo
             ht_nnew(2,k)=ht_new(2,k)-fact*fomega*metr_com%ht0(k,2)/fo
             ht_nnew(3,k)=ht_new(3,k)-fact*fomega*metr_com%ht0(k,3)/fo
          ENDDO
       ELSE
          ! General
          DO i=1,3
             ht_nnew(i,1)=ht_new(i,1) + dt_ions*fact*ht_for(i,1)&
                  - dt_ions*hv_new(i,1)
             ht_nnew(i,2)=ht_new(i,2) + dt_ions*fact*ht_for(i,2)&
                  - dt_ions*hv_new(i,2)
             ht_nnew(i,3)=ht_new(i,3) + dt_ions*fact*ht_for(i,3)&
                  - dt_ions*hv_new(i,3)
          ENDDO
       ENDIF
    ENDDO

    IF (paral%io_parent)&
         WRITE(6,'(/,A,I7,A)') ' HT_FORCES|  did not converge '
    DO i=1,ncolvar
       ! IF(abs(HT_DIFF(I)) .GT. TOL_CV_X ) THEN 
       cv_f_new(i) = 0.0_real_8! CV_F_OLD(I)/2.0_real_8
       ! ENDIF 
    ENDDO
100 CONTINUE

    CALL dcopy(ncolvar,cv_f_new,1,cv_f_old,1)

    ! ==--------------------------------------------------------------==
    DEALLOCATE(ht_ist_o,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(ht_ist_n,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(ht_ist_nn,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(det_ht_o,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(det_ht_n,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(det_ht_nn,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)

    DEALLOCATE(ht_diff,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)

    DEALLOCATE(asl,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(xlo,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(err,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(ipvt,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE ht_forces
  ! ==--------------------------------------------------------------==
  ! ==--------------------------------------------------------------==
  SUBROUTINE meta_cell_inp(iunit,line)
    ! ==--------------------------------------------------------------==
    ! ==  Reads Input for Metadynamics of the cell                    ==
    ! ==--------------------------------------------------------------==


    INTEGER                                  :: iunit
    CHARACTER(len=80)                        :: line

    CHARACTER(*), PARAMETER                  :: procedureN = 'meta_cell_inp'

    INTEGER                                  :: ia, ic, icv, ie, ierr, iht, &
                                                ii, ij, inum, IOuT
    LOGICAL                                  :: erread
    REAL(real_8)                             :: fmax, fmin, temp

    lmeta%lmeta_cell = .TRUE.
    ncolvar = 6


    IF (INDEX(line,'FULL').NE.0) THEN
       lfullc = .TRUE.
       lisoc  = .FALSE.
       l2dc   = .FALSE.
       ncolvar = 6
    ELSEIF (INDEX(line,'ISOT').NE.0) THEN
       lisoc  = .TRUE.
       lfullc = .FALSE.
       l2dc   = .FALSE.
       ncolvar = 3
    ELSEIF (INDEX(line,'2D').NE.0) THEN
       l2dc   = .TRUE.
       lfullc = .FALSE.
       lisoc  = .FALSE.
       ncolvar = 2
    ENDIF


    ! Allocate Memory
    ! CALL MEMORY(IP_CSCL_FAC,3*NCOLVAR,'CSCL_FAC')
    ALLOCATE(cscl_fac(3,ncolvar),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! CALL MEMORY(IP_TAD_SCF,NCOLVAR,'TAD_SCF')
    ALLOCATE(tad_scf(ncolvar),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! CALL MEMORY(IP_HT_IST,NCOLVAR,'HT_IST')
    ALLOCATE(ht_ist(ncolvar),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(ht_ist)!,ncolvar)
    ! CALL MEMORY(IP_DET_HT,NCOLVAR*9,'DET_HT')
    ALLOCATE(det_ht(9,ncolvar),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(det_ht)!,ncolvar*9)


    ! Set Default Values
#ifdef __SR11000
    !poption parallel, tlocal(IC)
    !voption indep(CSCL_FAC,TAD_SCF)
#endif
    DO ic = 1,ncolvar
       cscl_fac(1,ic) = 1.0_real_8
       cscl_fac(2,ic) = 0.2_real_8
       cscl_fac(3,ic) = 2._real_8
       tad_scf(ic)  = .FALSE.
    ENDDO

10  CONTINUE

    IF (paral%io_parent)&
         READ(iunit,err=20,END=20,fmt='(A)') line

    IF (INDEX(line,'END').NE.0 .AND. INDEX(line,'METADYN').NE.0)&
         GOTO 30


    IF (INDEX(line,'METASTEPNUM').NE.0) THEN
       IF (paral%io_parent)&
            READ(iunit,err=20,END=20,fmt=*) imeta%i_meta_max
       GOTO 10
       ! ==-------- RESTART META DYNAMICS FROM METASTEPI_META_RES -------==
    ELSEIF (INDEX(line,'META_RESTART').NE.0) THEN
       lmeta%meta_restart = .TRUE.
       IF (paral%io_parent)&
            READ(iunit,err=20,END=20,fmt=*) imeta%i_meta_res
       GOTO 10
       ! ==------------------ CV PARAMETERS  ----------------------------==
    ELSEIF(INDEX(line,'DEF').NE.0.&
         .AND. INDEX(line,'VARIABLE').NE.0) THEN
       IF (paral%io_parent)&
            READ(iunit,err=20,END=20,fmt=*) inum
       DO icv = 1,inum
          IF (paral%io_parent)&
               READ(iunit,err=20,END=20,fmt='(A)') line
          CALL readsi(line,1,IOuT,iht,erread)
          ii=INDEX(line,'SCA')
          ij=INDEX(line,'SCF')
          IF (ii.NE.0) THEN
             ia = ii+3
             tad_scf(iht) = .TRUE.
             CALL readsr(line,ia,ie,cscl_fac(1,iht),   erread)
             ia = ie
             CALL readsr(line,ia,ie,cscl_fac(2,iht),   erread)
             ia = ie
             CALL readsr(line,ia,ie,cscl_fac(3,iht),   erread)
             IF (erread) GOTO 22
          ELSEIF (ij.NE.0) THEN
             ia = ij+3
             tad_scf(iht) = .FALSE.
             CALL readsr(line,ia,ie,cscl_fac(1,iht),   erread)
             IF (erread) GOTO 22
          ENDIF
       ENDDO
       GOTO 10
       ! ==------------- DIMENSIONS OF HILLS IN CV SPACE-----------------==
    ELSEIF (INDEX(line,'HILLS').NE.0) THEN
       IF (INDEX(line,'OFF').NE.0) THEN
          lmeta%lhills=.FALSE.
       ELSE
          lmeta%lhills = .TRUE.
          IF (INDEX(line,'LOREN').NE.0) THEN
             lmeta%hlore = .TRUE.
          ELSEIF (INDEX(line,'RATIO').NE.0) THEN
             lmeta%hratio = .TRUE.
             lmeta%hshift = .FALSE.
             ij=INDEX(line,'POW')
             IF (ij.NE.0) THEN
                ia = ij+3
                CALL readsr(line,ia,ie,rmeta%expup,   erread)
                IF (erread) GOTO 20
                ia = ie
                CALL readsr(line,ia,ie,rmeta%expdown,   erread)
                IF (erread) GOTO 20
                ia = ie
                CALL readsr(line,ia,ie,rmeta%fboost,   erread)
                IF (erread) GOTO 20
             ENDIF
          ELSEIF (INDEX(line,'SHIFT').NE.0) THEN
             lmeta%hshift = .TRUE.
             lmeta%hratio = .FALSE.
             ij=INDEX(line,'RCUT')
             IF (ij.NE.0) THEN
                ia = ij+4
                CALL readsr(line,ia,ie,rmeta%rshift,   erread)
                IF (erread) GOTO 20
                ia = ie
                CALL readsr(line,ia,ie,rmeta%fboost,   erread)
                IF (erread) GOTO 20
             ENDIF

          ENDIF
          ii=INDEX(line,'=')
          IF (ii.NE.0) THEN
             ia = ii+1
             CALL readsr(line,ia,ie,rmeta%hllw,   erread)
             IF (erread) GOTO 20
             ia = ie
             CALL readsr(line,ia,ie,rmeta%hllh,   erread)
             IF (erread) GOTO 20
          ENDIF
       ENDIF
       GOTO 10
       ! ==---------------------- TUNING PARAMETERS ---------------------==
    ELSEIF (INDEX(line,'TUNING').NE.0) THEN
       IF (INDEX(line,'OFF').NE.0) THEN
          lmeta%ltune_hh     = .FALSE.
          lmeta%ltune_cscl   = .FALSE.
       ENDIF
       IF (INDEX(line,'HHEIGHT').NE.0) THEN
          lmeta%ltune_hh     = .TRUE.
          ii=INDEX(line,'=')
          IF (ii.NE.0) THEN
             ia = ii+1
             CALL readsr(line,ia,ie,rmeta%hlow,   erread)
             IF (erread) GOTO 20
             ia = ie
             CALL readsr(line,ia,ie,rmeta%htop,   erread)
             IF (erread) GOTO 20
             IF (rmeta%htop.LT.rmeta%hlow) THEN
                temp = rmeta%hlow
                rmeta%hlow = rmeta%htop
                rmeta%htop = temp
             ENDIF
          ENDIF
       ENDIF
       GOTO 10
       ! ==---------MAX # OF DYN. STEPS BETWEEN 2 METASTEPS--------------==
    ELSEIF(INDEX(line,'MAXSTEPNUM').NE.0 .AND.&
         INDEX(line,'INTERMETA') .NE.0) THEN
       IF (paral%io_parent)&
            READ(iunit,err=20,END=20,fmt=*) inter_hill_max
       GOTO 10
       ! ==---------MIN # OF DYN. STEPS BETWEEN 2 METASTEPS--------------==
    ELSEIF(INDEX(line,'MINSTEPNUM').NE.0 .AND.&
         INDEX(line,'INTERMETA') .NE.0) THEN
       IF (paral%io_parent)&
            READ(iunit,err=20,END=20,fmt=*) inter_hill
       GOTO 10
       ! ==------------ INTERVAL BETWEEN DISPLACEMENT CHECKS-------------==
    ELSEIF(INDEX(line,'CHECK') .NE.0 .AND.&
         INDEX(line,'DELAY').NE.0) THEN
       IF (paral%io_parent)&
            READ(iunit,err=20,END=20,fmt=*) imeta%icheck
       GOTO 10
       ! ==------------------------TOLL_FCNST----------------------------==
    ELSEIF(INDEX(line,'MOVEMENT') .NE.0 .AND.&
         INDEX(line,'CHECK').NE.0) THEN
       IF (paral%io_parent)&
            READ(iunit,err=20,END=20,fmt=*) toll_avcv
       GOTO 10
       ! ==----------- STORE RESTART AND TRAJECTORIES -------------------==
    ELSEIF (INDEX(line,'METASTORE') .NE.0) THEN
       IF (paral%io_parent)&
            READ(iunit,err=20,END=20,fmt=*) imeta%st_freq,imeta%tr_freq,imeta%qw_freq
       GOTO 10
       ! ==----- MAX ELECT. KIN. ENERGY (when above QUENCH BO) ----------==
    ELSEIF (INDEX(line,'MAXKINEN').NE.0) THEN
       IF (paral%io_parent)&
            READ(iunit,err=20,END=20,fmt='(A)') line
       CALL readsr(line,1,IOuT,rmeta%tolkin,erread)
       GOTO 10
       ! =------ RESTRAIN THE CELL VOLUME IN THE RANGE VOLMIN VOLMAX ----==
    ELSEIF (INDEX(line,'RESTR').NE.0.AND.INDEX(line,'VOLU').NE.0)THEN
       IF (paral%io_parent)&
            READ(iunit,err=20,END=20,fmt=*) fmin,fmax,mdcellr%vbarrier
       mdcellr%volmin = prcp_com%omega0*fmin
       mdcellr%volmax = prcp_com%omega0*fmax
       tvolbound = .TRUE.
       GOTO 10
    ELSE
       GOTO 10
    ENDIF  ! file biginning

20  CONTINUE
    IF (paral%io_parent)&
         WRITE(6,*) ' ERROR WHILE READING METADYNAMICS VARIABLES '
    IF (paral%io_parent)&
         WRITE(6,*) '       last line read: '
    IF (paral%io_parent)&
         WRITE(6,*)  line
    CALL stopgm('M_CELL_INP',' ',& 
         __LINE__,__FILE__)   
22  CONTINUE
    IF (paral%io_parent)&
         WRITE(6,*) ' ERROR WHILE READING SCAL. FAC. or MASS or KH.'
    CALL stopgm('M_CELL_INP',' ',& 
         __LINE__,__FILE__)

30  CONTINUE
    ! ==--------------------------------------------------------------==
    ! TESTS

    IF (.NOT. cntl%tprcp) CALL stopgm('M_CELL_INP',&
         'VARIABLE CELL MD IS REQUIRED',& 
         __LINE__,__FILE__)
    IF (ncolvar .EQ. 0) CALL stopgm('M_CELL_INP',&
         '# OF COLLECTIVE VARIABLES IS ZERO',& 
         __LINE__,__FILE__)
    DO ic = 1,ncolvar
       IF (tad_scf(ic)) THEN
          lmeta%ltune_cscl = .TRUE.
          GOTO 24
       ENDIF
    ENDDO
24  CONTINUE
    IF (.NOT. lmeta%lhills) THEN
       lmeta%ltune_hh     = .FALSE.
       lmeta%ltune_cscl   = .FALSE.
    ENDIF

    ! ==--------------------------------------------------------------==
    ! Print info 
    IF (paral%io_parent)&
         WRITE(6,*)
    IF (paral%io_parent)&
         WRITE(6,'(72("*"))')
    IF (paral%io_parent)&
         WRITE(6,'(1x,A,A)')&
         '***          META DYNAMICS OF THE VARIABLE CELL (PR)',&
         '          ***'
    IF (paral%io_parent)&
         WRITE(6,'(4x,A,I5)') '- Number of Collective Variables =',&
         ncolvar

    IF ((l2dc).AND.paral%io_parent)&
         WRITE(6,'(4x,A)') '- ONLY 2 SIDE LENGTHS '
    IF ((lisoc).AND.paral%io_parent)&
         WRITE(6,'(4x,A)') '- SIDE LENGTHS '
    IF ((lfullc).AND.paral%io_parent)&
         WRITE(6,'(4x,A)') '- FREE CELL '
    IF (paral%io_parent)&
         WRITE(6,'(6x,A)')&
         'colvar_mtd   : CV values and scaling factors'
    IF (paral%io_parent)&
         WRITE(6,'(6x,A/6x,A/6x,A)')&
         'dcolvar_mtd  : CV diffusivities,',&
         '              norm of the displacement along the traj.',&
         '              transversal hill width (fixed),',&
         '              and hills heights HLLH'
    IF (paral%io_parent)&
         WRITE(6,'(6x,A/6x,A/6x,A)')&
         'forfact_mtd  : av. forces from underlying pot. ,',&
         '               forces from time dependent pot. at this step',&
         '               total  forces averaged between 2 metastep '

    IF (paral%io_parent)&
         WRITE(6,'(6x,A)')&
         'dispvar_mtd  : CV displacements'

    IF (paral%io_parent)&
         WRITE(6,'(6x,A/6x,A/9x,A)')&
         'envar_mtd  : Ekin of the cell, ',&
         '             gaussians pot. (time dependent term), ETOT',&
         '            EHAM, EHAM + Gaussians '


    IF (paral%io_parent)&
         WRITE(6,'(6x,A,F12.6)')&
         '- Minimum. CV Displacement Required to Place New a Hill',&
         toll_avcv
    IF (paral%io_parent)&
         WRITE(6,'(6x,A,I5,A)')&
         '- Update of Gaussian Pot. every ',&
         inter_hill,' MD STEPS (WHEN CV MOVED ENOUGH)'
    IF (paral%io_parent)&
         WRITE(6,'(6x,A,I10)')&
         '- MAX. NUMBER OF MD STEPS BEFORE NEXT UPDATE',inter_hill_max
    IF (paral%io_parent)&
         WRITE(6,'(6x,A,I10,A)')&
         '- RESTART File Saved every ', imeta%st_freq, ' Metasteps'
    IF (paral%io_parent)&
         WRITE(6,'(6x,A,I10,A)')&
         '- TRAJECTORY File Appended every', imeta%tr_freq,' Metasteps'
    IF (paral%io_parent)&
         WRITE(6,'(6x,A,I10,A)')&
         '- QUENCH BO performed every', imeta%qw_freq, ' Metasteps'
    IF (cntl%tc .AND. rmeta%tolkin .LT. cntr%ekinw+cntr%toll) THEN
       IF (paral%io_parent)&
            WRITE(6,'(6x,A,F12.6)')&
            '- QUENCH BO performed when EL. Ekin is > ', rmeta%tolkin
    ELSEIF (rmeta%tolkin.NE.-1.0_real_8) THEN
       IF (paral%io_parent)&
            WRITE(6,'(6x,A,F12.6)')&
            '- QUENCH BO performed when EL. Ekin is > ', rmeta%tolkin
    ENDIF

    IF ((lmeta%meta_restart).AND.paral%io_parent)&
         WRITE(6,'(6x,A,1X,I8)')&
         '- Restart from a previous Metadynamics at step',imeta%i_meta_res
    IF (lmeta%lhills) THEN
       IF (paral%io_parent)&
            WRITE(6,'(6x,A)')&
            '- Time Dependent Potential Active: '
       IF (paral%io_parent)&
            WRITE(6,'(15x,A,F9.5,A,F9.5)') ' Transversal Width = ',&
            rmeta%hllw,'    Initial Height = ',rmeta%hllh
       IF (lmeta%hshift) THEN
          IF (paral%io_parent)&
               WRITE(6,'(6x,A/15x,A,F8.3,A,f8.3)')&
               '- Hills Shape : Truncated Gaussians',&
               ' Rcut = ',rmeta%rshift,'     Width Factor = ',rmeta%fboost
       ELSEIF (lmeta%hratio) THEN
          IF (paral%io_parent)&
               WRITE(6,'(6x,A/15x,A,F8.3,F8.3,A,f8.3)')&
               '- Hills Shape : Rational Function',&
               ' Exponents : ',rmeta%expup,rmeta%expdown,&
               '     Width Factor = ',rmeta%fboost
       ENDIF
       IF ((lmeta%ltune_hh) .AND.paral%io_parent)&
            WRITE(6,'(6x,A/15x,A,F10.6,A,f10.6,A/)')&
            '- Hills Height tuned on the Underlying potential',&
            'min. Height ',rmeta%hlow*au_kcm,' Kcal,    max Height ',&
            rmeta%htop*au_kcm,' Kcal.'
       IF ((lmeta%ltune_cscl).AND.paral%io_parent)&
            WRITE(6,'(6x,A)')&
            '- Some scaling factors might be adjusted along the run'
    ENDIF

    IF (tvolbound) THEN
       IF (paral%io_parent)&
            WRITE(6,'(4x,A,f12.5,A,f12.5,A)')&
            '- Volume Restrained in range: ',&
            mdcellr%volmin, '-', mdcellr%volmax, ' a.u.'
    ENDIF

    IF (paral%io_parent)&
         WRITE(6,'(6x,A)') '- CV Parameters: '
    DO iht = 1,ncolvar
       IF (paral%io_parent)&
            WRITE(6,'(6x,A,I4,A,f8.4)')&
            '          ', iht, ' SC ', cscl_fac(1,iht)
    ENDDO

    ! stop 'input'
    IF (paral%io_parent)&
         WRITE(6,'(80("*"))')
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE meta_cell_inp
  ! ==================================================================
  SUBROUTINE parfull(ht,ht_ist,det)

    REAL(real_8)                             :: ht(3,3), ht_ist(6), det(9,6)

    REAL(real_8)                             :: al, bl, cl, cosa, cosb, cosg, &
                                                htemp(3,3,6), num, ual, &
                                                ualbl, ualcl, ubl, ublcl, ucl

    al = ht(1,1)*ht(1,1) +&
         ht(1,2)*ht(1,2) +&
         ht(1,3)*ht(1,3)
    al = SQRT(al)

    bl = ht(2,1)*ht(2,1) +&
         ht(2,2)*ht(2,2) +&
         ht(2,3)*ht(2,3)
    bl = SQRT(bl)
    cosg = ht(1,1)*ht(2,1)+&
         ht(1,2)*ht(2,2)+&
         ht(1,3)*ht(2,3)
    cosg = cosg/(al*bl)

    cl = ht(3,1)*ht(3,1) +&
         ht(3,2)*ht(3,2) +&
         ht(3,3)*ht(3,3)
    cl = SQRT(cl)

    ual = 1.0_real_8/al
    ubl = 1.0_real_8/bl
    ucl = 1.0_real_8/cl
    ualbl = ual*ubl
    ualcl = ual*ucl
    ublcl = ubl*ucl

    cosg = ht(1,1)*ht(2,1)+&
         ht(1,2)*ht(2,2)+&
         ht(1,3)*ht(2,3)
    cosg = cosg*ualbl

    cosb = ht(1,1)*ht(3,1)+&
         ht(1,2)*ht(3,2)+&
         ht(1,3)*ht(3,3)
    cosb = cosb*ualcl

    cosa = ht(2,1)*ht(3,1)+&
         ht(2,2)*ht(3,2)+&
         ht(2,3)*ht(3,3)
    cosa = cosa*ublcl

    CALL zeroing(htemp)!,9*6)


    htemp(1,1,1) =ual*ht(1,1)
    htemp(1,2,1) =ual*ht(1,2)
    htemp(1,3,1) =ual*ht(1,3)

    htemp(2,1,2) =ubl*ht(2,1)
    htemp(2,2,2) =ubl*ht(2,2)
    htemp(2,3,2) =ubl*ht(2,3)

    htemp(3,1,3) =ucl*ht(3,1)
    htemp(3,2,3) =ucl*ht(3,2)
    htemp(3,3,3) =ucl*ht(3,3)

    num =  ht(2,1)*ht(3,1)+ ht(2,2)*ht(3,2)+ ht(2,3)*ht(3,3)
    htemp(2,1,4) =ublcl*(ht(3,1)-num*htemp(2,1,2)*ubl)
    htemp(3,1,4) =ublcl*(ht(2,1)-num*htemp(3,1,3)*ucl)
    htemp(2,2,4) =ublcl*(ht(3,2)-num*htemp(2,2,2)*ubl)
    htemp(3,2,4) =ublcl*(ht(2,2)-num*htemp(3,2,3)*ucl)
    htemp(2,3,4) =ublcl*(ht(3,3)-num*htemp(2,3,2)*ubl)
    htemp(3,3,4) =ublcl*(ht(2,3)-num*htemp(3,3,3)*ucl)


    num =  ht(1,1)*ht(3,1)+ ht(1,2)*ht(3,2)+ ht(1,3)*ht(3,3)
    htemp(1,1,5) =ualcl*(ht(3,1)-num*htemp(1,1,1)*ual)
    htemp(3,1,5) =ualcl*(ht(1,1)-num*htemp(3,1,3)*ucl)
    htemp(1,2,5) =ualcl*(ht(3,2)-num*htemp(1,2,1)*ual)
    htemp(3,2,5) =ualcl*(ht(1,2)-num*htemp(3,2,3)*ucl)
    htemp(1,3,5) =ualcl*(ht(3,3)-num*htemp(1,3,1)*ual)
    htemp(3,3,5) =ualcl*(ht(1,3)-num*htemp(3,3,3)*ucl)

    num =  ht(1,1)*ht(2,1)+ ht(1,2)*ht(2,2)+ ht(1,3)*ht(2,3)
    htemp(1,1,6) =ualbl*(ht(2,1)-num*htemp(1,1,1)*ual)
    htemp(2,1,6) =ualbl*(ht(1,1)-num*htemp(2,1,2)*ubl)
    htemp(1,2,6) =ualbl*(ht(2,2)-num*htemp(1,2,1)*ual)
    htemp(2,2,6) =ualbl*(ht(1,2)-num*htemp(2,2,2)*ubl)
    htemp(1,3,6) =ualbl*(ht(2,3)-num*htemp(1,3,1)*ual)
    htemp(2,3,6) =ualbl*(ht(1,3)-num*htemp(2,3,2)*ubl)


    CALL zeroing(det)!,6*6)

    CALL dcopy(9,htemp(1,1,1),1,det(1,1),1)
    CALL dcopy(9,htemp(1,1,2),1,det(1,2),1)
    CALL dcopy(9,htemp(1,1,3),1,det(1,3),1)
    CALL dcopy(9,htemp(1,1,4),1,det(1,4),1)
    CALL dcopy(9,htemp(1,1,5),1,det(1,5),1)
    CALL dcopy(9,htemp(1,1,6),1,det(1,6),1)

    ht_ist(1) = al
    ht_ist(2) = bl
    ht_ist(3) = cl
    ht_ist(4) = cosa
    ht_ist(5) = cosb
    ht_ist(6) = cosg

    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE parfull
  ! ==================================================================
  SUBROUTINE pariso(ht,ht_ist,det)

    REAL(real_8)                             :: ht(3,3), ht_ist(3), det(6,3)

    REAL(real_8)                             :: al, bl, cl, cosa, cosb, cosg, &
                                                htemp(3,3,3)

    al = ht(1,1)*ht(1,1) +&
         ht(1,2)*ht(1,2) +&
         ht(1,3)*ht(1,3)
    al = SQRT(al)

    bl = ht(2,1)*ht(2,1) +&
         ht(2,2)*ht(2,2) +&
         ht(2,3)*ht(2,3)
    bl = SQRT(bl)
    cosg = ht(1,1)*ht(2,1)+&
         ht(1,2)*ht(2,2)+&
         ht(1,3)*ht(2,3)
    cosg = cosg/(al*bl)

    cl = ht(3,1)*ht(3,1) +&
         ht(3,2)*ht(3,2) +&
         ht(3,3)*ht(3,3)
    cl = SQRT(cl)
    cosb = ht(1,1)*ht(3,1)+&
         ht(1,2)*ht(3,2)+&
         ht(1,3)*ht(3,3)
    cosb = cosb/(al*cl)
    cosa = ht(2,1)*ht(3,1)+&
         ht(2,2)*ht(3,2)+&
         ht(2,3)*ht(3,3)
    cosa = cosa/(bl*cl)


    CALL zeroing(htemp)!,3*3*3)

    htemp(1,1,1) =1._real_8/al*ht(1,1)
    htemp(1,2,1) =1._real_8/al*ht(1,2)
    htemp(1,3,1) =1._real_8/al*ht(1,3)

    htemp(2,1,2) =1._real_8/bl*ht(2,1)
    htemp(2,2,2) =1._real_8/bl*ht(2,2)
    htemp(1,3,2) =1._real_8/bl*ht(2,3)

    htemp(3,1,3) =1._real_8/cl*ht(3,1)
    htemp(3,2,3) =1._real_8/cl*ht(3,2)
    htemp(3,3,3) =1._real_8/cl*ht(3,3)

    CALL zeroing(det)!,6*3)

    CALL dcopy(9,htemp(1,1,1),1,det(1,1),1)
    CALL dcopy(9,htemp(1,1,2),1,det(1,2),1)
    CALL dcopy(9,htemp(1,1,3),1,det(1,3),1)

    ht_ist(1) = al
    ht_ist(2) = bl
    ht_ist(3) = cl

    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE pariso
  ! ==================================================================
  SUBROUTINE delrot_cell(ha,hb)

    REAL(real_8)                             :: ha(3,3), hb(3,3)

    INTEGER                                  :: ind(4), iv
    REAL(real_8)                             :: d(3,4), d_r(3,3,3,4), &
                                                rot(3,3), transl(3), &
                                                vec_a(3,4), vec_b(3,4), x, &
                                                xt, y, yt, z, zt

    vec_a(1,1) = 0.0_real_8
    vec_a(2,1) = 0.0_real_8
    vec_a(3,1) = 0.0_real_8

    vec_b(1,1) = 0.0_real_8
    vec_b(2,1) = 0.0_real_8
    vec_b(3,1) = 0.0_real_8
    DO iv = 1,4
       ind(iv) = iv
    ENDDO
    DO iv = 1,3
       vec_a(1,iv+1) = ha(iv,1)
       vec_a(2,iv+1) = ha(iv,2)
       vec_a(3,iv+1) = ha(iv,3)
       vec_b(1,iv+1) = hb(iv,1)
       vec_b(2,iv+1) = hb(iv,2)
       vec_b(3,iv+1) = hb(iv,3)
    ENDDO

    CALL rmsd_rot(4,vec_b,vec_a,ind,transl,rot,4,d,d_r,.TRUE.)

    DO iv = 1,3
       xt = vec_b(1,iv+1)
       yt = vec_b(2,iv+1)
       zt = vec_b(3,iv+1)
       x =rot(1,1)*xt + rot(1,2)*yt + rot(1,3)*zt
       y =rot(2,1)*xt + rot(2,2)*yt + rot(2,3)*zt
       z =rot(3,1)*xt + rot(3,2)*yt + rot(3,3)*zt
       hb(iv,1) = x
       hb(iv,2) = y
       hb(iv,3) = z
    ENDDO

    RETURN
  END SUBROUTINE delrot_cell
  ! ==================================================================

END MODULE meta_cell_utils
