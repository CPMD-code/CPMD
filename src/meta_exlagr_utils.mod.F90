MODULE meta_exlagr_utils
  USE cnst,                            ONLY: factem,&
                                             pi
  USE cnst_dyn,                        ONLY: &
       cscl_fac, cscl_val, cv_dyn, cv_ist, cv_mass, cv_path, cv_temp, cv_vel, &
       dstrmeta, ekincv, ekincv_walk, fhills, fmtdres, fvbound, hllh_val, &
       hllw_val, iangcv, icv_spin, imeta, inter_hill, inter_hill_max, kharm, &
       lcvtc, lmeta, mdcellr, ncolvar, rcc, rmeta, tad_scf, tcvscale, &
       tvolbound, tycvar, vharm, vharm_walk
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_app,&
                                             fo_def,&
                                             fo_new,&
                                             fo_old,&
                                             fo_ufo,&
                                             fo_verb
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE latgen_utils,                    ONLY: omegagen
  USE metr,                            ONLY: metr_com
  USE mw,                              ONLY: mwi
  USE nose,                            ONLY: glib
  USE parac,                           ONLY: paral
  USE prng_utils,                      ONLY: repprngu
  USE readsr_utils,                    ONLY: xstring
  USE ropt,                            ONLY: iteropt
  USE strs,                            ONLY: alpha,&
                                             beta
  USE system,                          ONLY: parm
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ekincv_global
  PUBLIC :: meta_stress
  PUBLIC :: rscvelcv
  PUBLIC :: rinvelcv
  PUBLIC :: scf_tune
  PUBLIC :: cv_read_out
  PUBLIC :: cv_exlagr_out2
  PUBLIC :: wmtdres
  PUBLIC :: calc_pos_dyn
  PUBLIC :: rmtdres
  PUBLIC :: cv_exlagr_out

CONTAINS

  ! ================================================================== 
  SUBROUTINE calc_pos_dyn(tollm,i_cvst,i_temp,i_meta,hc_last,&
       hill_add)
    ! ==--------------------------------------------------------------==

    REAL(real_8)                             :: tollm
    INTEGER                                  :: i_cvst, i_temp, i_meta
    REAL(real_8)                             :: hc_last(ncolvar)
    LOGICAL                                  :: hill_add

    CHARACTER(*), PARAMETER                  :: procedureN = 'calc_pos_dyn'

    INTEGER                                  :: icv, ierr, ii
    REAL(real_8)                             :: av_disp, diff, min_disp
    REAL(real_8), ALLOCATABLE                :: scr(:)

    ALLOCATE(scr(ncolvar),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! Store CV values and Check averages

    IF (i_cvst .EQ.  inter_hill ) THEN
       min_disp = 100.0_real_8
       ii = inter_hill
       av_disp  = 0.0_real_8
       DO icv = 1,ncolvar
          diff = (cv_dyn(icv)-hc_last(icv))
          IF (iangcv(icv) .EQ. 1 .AND. diff .GT. pi) THEN
             diff = diff -2.0_real_8*pi
          ELSEIF (iangcv(icv) .EQ. 1 .AND. diff .LT. -pi) THEN
             diff = diff +2.0_real_8*pi
          ENDIF
          diff = (diff/cscl_fac(1,icv))**2.0_real_8
          av_disp = av_disp + diff
          min_disp = MIN(diff,min_disp)
       ENDDO
       av_disp = SQRT(av_disp)
       IF (av_disp .GE. tollm .OR. i_meta .EQ. 1) hill_add = .TRUE.

       IF (paral%io_parent)&
            WRITE(6,'(10x,A,f10.6,A,f10.6,A,f10.6/)')&
            'MIN DISP. = ',SQRT(min_disp),'  ||DISP.tot|| = ', av_disp,&
            '  Tolerance ',tollm
    ELSEIF (i_cvst .GT.  inter_hill ) THEN

       i_temp = i_temp + 1
       ii = inter_hill
       min_disp = 100.0_real_8
       DO icv = 1,ncolvar
          diff = (cv_dyn(icv)-hc_last(icv))
          IF (iangcv(icv) .EQ. 1 .AND. diff .GT. pi) THEN
             diff = diff -2.0_real_8*pi
          ELSEIF (iangcv(icv) .EQ. 1 .AND. diff .LT. -pi) THEN
             diff = diff +2.0_real_8*pi
          ENDIF
          diff = (diff/cscl_fac(1,icv))**2.0_real_8
          scr(icv) = diff
          min_disp = MIN(diff,min_disp)
       ENDDO
       IF (MOD(i_temp,imeta%icheck) .EQ. 0 ) THEN
          av_disp = 0.0_real_8
          DO icv = 1,ncolvar
             av_disp = av_disp+scr(icv)
          ENDDO
          av_disp = SQRT(av_disp)
          IF (paral%io_parent)&
               WRITE(6,'(10x,A,f10.6,A,f10.6,A,f10.6/)')&
               'MIN DISP. = ',SQRT(min_disp),'  ||DISP.tot|| = ',av_disp,&
               '  Tolerance ',tollm
          IF (av_disp .GE. tollm .OR. i_cvst .GE. inter_hill_max) THEN
             hill_add = .TRUE.
          ENDIF
       ENDIF

    ENDIF   ! I_CVST
    ! ==--------------------------------------------------------------==
    DEALLOCATE(scr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==

    RETURN
  END SUBROUTINE calc_pos_dyn

  ! ==================================================================
  SUBROUTINE rinvelcv(vel,mass,cscl_fac,nvar,temp,ekincv)
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: nvar
    REAL(real_8)                             :: cscl_fac(3,nvar), mass(nvar), &
                                                vel(nvar), temp, ekincv

    INTEGER                                  :: icv
    REAL(real_8)                             :: alfa, const, rnr, sigma, &
                                                tempp, tscal, vscale

    CALL zeroing(vel)
    IF (temp.LT.1.e-5_real_8) GOTO 100
    rnr = repprngu()
    DO icv=1,nvar
       sigma=SQRT(temp/(mass(icv)*factem))
       rnr = repprngu()
       alfa=2.0_real_8*pi*rnr
       vel(icv)&
            =SQRT(LOG(repprngu())*(-2.0_real_8))*COS(alfa)*sigma
    ENDDO
    ! ==--------------------------------------------------------------==
    ! ==  CALCULATE KINETIC ENERGY OF THE CV                          ==
    ! ==--------------------------------------------------------------==
    ekincv= 0._real_8
    tempp = 0._real_8
    DO icv=1,nvar
       const=0.5_real_8*mass(icv)*vel(icv)*vel(icv)
       ekincv=ekincv+const
       tempp = tempp + const/(cscl_fac(1,icv)*cscl_fac(1,icv))
       ! TEMPP = TEMPP + CONST
    ENDDO
    tempp=tempp*factem*2.0_real_8/REAL(nvar,kind=real_8)
    IF (tempp.GT.1.e-5_real_8) THEN
       tscal=temp/tempp
       vscale=SQRT(tscal)
       !$omp parallel do private(ICV)
       DO icv=1,nvar
          vel(icv)=vel(icv)*vscale
       ENDDO
    ENDIF

100 CONTINUE

    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rinvelcv

  ! ==================================================================
  SUBROUTINE rscvelcv(temp1,temp2,tempp)
    ! ==--------------------------------------------------------------==

    REAL(real_8)                             :: temp1, temp2, tempp

    INTEGER                                  :: icv, ipw, iwalk1
    REAL(real_8)                             :: const, dif_cv, tscal, vscale

! ==--------------------------------------------------------------==
! ==  Dynamical rescaling factor (cv_temp0/temp_IST),             ==
! ==  where tempp is calculated every step                        ==
! ==--------------------------------------------------------------==

    ekincv = 0._real_8
    vharm  = 0.0_real_8
    tempp  = 0.0_real_8

    IF (lmeta%tcvanalysis) RETURN

    DO iwalk1=1,mwi%nwalk
       vharm  = 0.0_real_8
       ekincv = 0._real_8
       tempp  = 0.0_real_8
       ipw=(iwalk1-1)*ncolvar
       DO icv=1,ncolvar
          dif_cv = cv_ist(icv+ipw)-cv_dyn(icv+ipw)
          IF (tycvar(icv) .EQ. 3 .AND. dif_cv .GT. pi) THEN
             dif_cv = dif_cv -2.0_real_8*pi
          ELSEIF (tycvar(icv) .EQ. 3 .AND. dif_cv .LT. -pi) THEN
             dif_cv = dif_cv +2.0_real_8*pi
          ENDIF
          const  = 0.5_real_8*cv_mass(icv)*cv_vel(icv+ipw)*cv_vel(icv+ipw)
          ekincv = ekincv + const

          vharm  = vharm  + kharm(icv)*dif_cv* dif_cv
          IF (tcvscale)THEN
             tempp = tempp + const/(cscl_fac(1,icv)*cscl_fac(1,icv))
          ELSE
             tempp = tempp + const
          ENDIF
       ENDDO
       tempp=tempp*factem*2.0_real_8/REAL(ncolvar,kind=real_8)

       IF (lmeta%tlocalizespin) THEN
          DO icv=1,ncolvar
             IF (icv_spin(icv).NE. 0) THEN
                dif_cv = cv_ist(icv+ipw)-cv_dyn(icv+ipw)
                ! mb         VHARM  = VHARM  - KHARM(ICV)*DIF_CV* DIF_CV ! this is a bug !
                vharm  = vharm  + kharm(icv)*dif_cv* dif_cv! cmb-spin
             ENDIF
          ENDDO
       ENDIF

       IF (lcvtc) THEN
          IF (tempp.GT.temp1.OR.tempp.LT.temp2.AND.tempp.NE.0._real_8) THEN
             tscal=cv_temp/tempp
             vscale=SQRT(tscal)
             ekincv = 0.0_real_8
             tempp  = 0.0_real_8
             ! !$OMP parallel do private(ICV,CONST) reduction(+:EKINCV,TEMPP)
             DO icv=1,ncolvar
                cv_vel(icv+ipw)=cv_vel(icv+ipw)*vscale
                const  = 0.5_real_8*cv_mass(icv)*&
                     cv_vel(icv+ipw)*cv_vel(icv+ipw)
                ekincv = ekincv + const
                tempp  = tempp + const/(cscl_fac(1,icv)*cscl_fac(1,icv))
             ENDDO
             tempp=tempp*factem*2.0_real_8/REAL(ncolvar,kind=real_8)
          ENDIF
       ENDIF
       IF (mwi%nwalk.GT.1)THEN
          ekincv_walk(iwalk1)=ekincv
          vharm_walk(iwalk1)=vharm
       ENDIF
    ENDDO! IWALK1
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rscvelcv
  ! ==================================================================
  SUBROUTINE  ekincv_global(ek_cv)
    ! ==--------------------------------------------------------------==

    REAL(real_8)                             :: ek_cv

    INTEGER                                  :: icv
    REAL(real_8)                             :: const, dif_cv

! ==--------------------------------------------------------------==
! ==  Kinetic Energy of the Collective Variable for global rescaling
! ==--------------------------------------------------------------==

    ek_cv = 0._real_8

    IF (lmeta%tcvanalysis) RETURN

    DO icv=1,ncolvar
       dif_cv = cv_ist(icv)-cv_dyn(icv)
       IF (tycvar(icv) .EQ. 3 .AND. dif_cv .GT. pi) THEN
          dif_cv = dif_cv -2.0_real_8*pi
       ELSEIF (tycvar(icv) .EQ. 3 .AND. dif_cv .LT. -pi) THEN
          dif_cv = dif_cv +2.0_real_8*pi
       ENDIF
       const  = 0.5_real_8*cv_mass(icv)*cv_vel(icv)*cv_vel(icv)
       ek_cv = ek_cv + const

    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE ekincv_global
  ! ==================================================================
  SUBROUTINE  cv_exlagr_out(i_meta,hc_last,f_harm,f_hill,f_wall,&
       cvd_ist,Displacement,ekinc,ekinp)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: i_meta
    REAL(real_8) :: hc_last(ncolvar), f_harm(ncolvar), f_hill(ncolvar), &
      f_wall(ncolvar), cvd_ist(ncolvar), Displacement, ekinc, ekinp

    CHARACTER(len=20), PARAMETER :: file1 = 'colvar_mtd', &
      file2 = 'parvar_mtd', file3 = 'istvar_mtd', file4 = 'forfac_mtd', &
      file5 = 'disvar_mtd', file6 = 'velvar_mtd', file7 = 'enevar_mtd', &
      file8 = 'spinpo_mtd'

    CHARACTER(len=10)                        :: chnum
    CHARACTER(len=100)                       :: lineform
    INTEGER                                  :: iaa, icv, iee
    INTEGER, SAVE                            :: ifirst = fo_verb
    LOGICAL                                  :: ferror
    REAL(real_8)                             :: t_ion

! ==--------------------------------------------------------------==
! Open output files

    IF (paral%io_parent)&
         CALL fileopen(51,file1,fo_app+ifirst,ferror)
    IF (paral%io_parent)&
         CALL fileopen(52,file2,fo_app+ifirst,ferror)
    IF (paral%io_parent)&
         CALL fileopen(53,file3,fo_app+ifirst,ferror)
    IF (paral%io_parent)&
         CALL fileopen(54,file4,fo_app+ifirst,ferror)
    IF (paral%io_parent)&
         CALL fileopen(55,file5,fo_app+ifirst,ferror)
    IF (paral%io_parent)&
         CALL fileopen(56,file6,fo_app+ifirst,ferror)
    IF (paral%io_parent)&
         CALL fileopen(57,file7,fo_app+ifirst,ferror)
    IF (lmeta%tlocalizespin) THEN
       IF (paral%io_parent)&
            CALL fileopen(58,file8,fo_app+ifirst,ferror)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Verbose open only on the first run
    ifirst=0
    ! ==--------------------------------------------------------------==
    ! Write output
    IF (paral%io_parent)&
         WRITE(chnum,'(I5)') ncolvar
    CALL xstring(chnum,iaa,iee)
    lineform =&
         '(1X,I7,'//chnum(iaa:iee)//'E14.6,'//chnum(iaa:iee)//'E13.3)'
    IF (paral%io_parent)&
         WRITE(51,lineform) iteropt%nfi,(hc_last(icv),icv=1,ncolvar),&
         (cscl_fac(1,icv),icv=1,ncolvar)

    lineform =&
         '(1X,I7,3F14.6)'
    IF (paral%io_parent)&
         WRITE(52,lineform) iteropt%nfi,&
         Displacement,rmeta%hllw,rmeta%hllh
    ! Close files
    IF (paral%io_parent)&
         CALL fileclose(51)
    IF (paral%io_parent)&
         CALL fileclose(52)

    lineform =&
         '(1X,I7,'//chnum(iaa:iee)//'E14.6,'//chnum(iaa:iee)//'E14.6)'
    IF (paral%io_parent)&
         WRITE(53,lineform) iteropt%nfi,(cv_ist(icv),icv=1,ncolvar),&
         (cv_ist(icv)-hc_last(icv),icv=1,ncolvar)
    ! Close files
    IF (paral%io_parent)&
         CALL fileclose(53)

    ! Print force factors
    lineform =&
         '(1X,I7,'//chnum(iaa:iee)//'E14.6,'&
         //chnum(iaa:iee)//'E14.6,'//chnum(iaa:iee)//'E16.4)'
    IF (paral%io_parent)&
         WRITE(54,lineform) iteropt%nfi,(f_harm(icv),icv=1,ncolvar),&
         (f_hill(icv),icv=1,ncolvar),(f_wall(icv),icv=1,ncolvar)

    ! Close files
    IF (paral%io_parent)&
         CALL fileclose(54)

    ! Print last displacements
    lineform =&
         '(I6,'//chnum(iaa:iee)//'f11.6,'//chnum(iaa:iee)//&
         'f10.6,'//chnum(iaa:iee)//'f10.6)'
    IF (paral%io_parent)&
         WRITE(55,lineform) iteropt%nfi,&
         ((cv_dyn(icv)-hc_last(icv)),icv=1,ncolvar),&
         (cvd_ist(icv),icv=1,ncolvar),&
         (kharm(icv),icv=1,ncolvar)

    ! Close files
    IF (paral%io_parent)&
         CALL fileclose(55)

    ! Print last velocities, temperature and kinetic energy
    ! EHAM_HILL=EKINP+ETOT+ENOSE+ENOSP+ECNSTR+EKINC+VHARM+EKINC
    lineform =&
         '(I6,'//chnum(iaa:iee)//'f11.6,f15.6,6f16.8)'
    IF (paral%io_parent)&
         WRITE(56,lineform) iteropt%nfi,(cv_vel(icv),icv=1,ncolvar)

    IF (paral%io_parent)&
         CALL fileclose(56)

    t_ion = ekinp*factem*2._real_8/glib
    lineform = '(I6,8E16.6)'
    IF (paral%io_parent)&
         WRITE(57,lineform) iteropt%nfi,&
         t_ion,ekinc,ekincv,vharm,rmeta%gausspot,&
         ener_com%etot,rmeta%eham_hill,rmeta%eham_hill+rmeta%gausspot
    IF (paral%io_parent)&
         CALL fileclose(57)

    IF (lmeta%tlocalizespin) THEN
       ! Print Position of the center of the spin density
       IF (paral%io_parent)&
            WRITE(58,'(I6,3f12.6)') iteropt%nfi, rcc(1),rcc(2),rcc(3)
       IF (paral%io_parent)&
            CALL fileclose(58)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE cv_exlagr_out

  ! ==================================================================
  SUBROUTINE scf_tune(i_meta,hc_last)

    INTEGER                                  :: i_meta
    REAL(real_8)                             :: hc_last(ncolvar)

    CHARACTER(*), PARAMETER                  :: procedureN = 'scf_tune'
    INTEGER, PARAMETER                       :: nevent = 2 
    REAL(real_8), PARAMETER                  :: damp = 0.750_real_8 

    INTEGER                                  :: icv, ierr, ninter
    INTEGER, ALLOCATABLE, SAVE               :: ireset(:,:)
    INTEGER, SAVE                            :: i0, ifirst = 0
    REAL(real_8)                             :: elong, fact, hel1, hel2, &
                                                scf_new
    REAL(real_8), ALLOCATABLE, SAVE          :: cv_max(:), cv_med(:), &
                                                cv_min(:)

    IF (ifirst .EQ. 0) THEN
       i0 = i_meta
       ALLOCATE(cv_max(ncolvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(cv_max)!,ncolvar)
       ALLOCATE(cv_min(ncolvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(cv_min)!,ncolvar)
       ALLOCATE(cv_med(ncolvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(cv_med)!,ncolvar)
       ALLOCATE(ireset(3,ncolvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(ireset)!,3*ncolvar)
       ifirst = 1
    ENDIF

    ninter = 3000/inter_hill
    DO icv = 1,ncolvar
       IF (.NOT.tad_scf(icv)) THEN
          GOTO 15
       ELSEIF (tad_scf(icv) .AND. ireset(1,icv) .EQ. 0) THEN
          cv_max(icv) = hc_last(icv)
          cv_min(icv) = hc_last(icv)
          cv_med(icv) = hc_last(icv)
          ireset(1,icv) = 1
       ELSE
          cv_max(icv) = MAX(cv_max(icv),hc_last(icv))
          cv_min(icv) = MIN(cv_min(icv),hc_last(icv))
          cv_med(icv) = cv_med(icv)+hc_last(icv)
          IF (MOD(i_meta-i0,ninter).EQ. 0) THEN
             cv_med(icv) = cv_med(icv)/REAL(ninter,kind=real_8)
             elong = ABS(cv_max(icv)-cv_min(icv))/2.0_real_8
             hel1  = ABS(cv_max(icv)-cv_med(icv))
             hel2  = ABS(cv_med(icv)-cv_min(icv))
             fact = rmeta%hllw*cscl_fac(1,icv)
             IF (1.2_real_8*damp*elong .LT.fact) THEN
                ireset(2,icv) = ireset(2,icv) + 1
                ireset(3,icv) = 0
                IF (paral%io_parent)&
                     WRITE(6,'(/10x,A,I3,A,I5,A)')&
                     '******** WARNING n^ ', ireset(2,icv),&
                     'for Scale Factor of CV n^ ',&
                     icv ,' ********'
                IF (paral%io_parent)&
                     WRITE(6,'(10x,A,I5,A,f14.6/10x,A,f14.6/)')&
                     '    in the last ',ninter,' mstep MAX-MIN = ',elong,&
                     '  whereas HLLW*SCF = ', fact
             ELSEIF(0.8_real_8*damp*hel1 .GT. fact .AND.&
                  0.8_real_8*damp*hel2 .GT. fact) THEN
                ireset(3,icv) = ireset(3,icv) + 1
                ireset(2,icv) = 0
                IF (paral%io_parent)&
                     WRITE(6,'(/10x,A,I3,A,I5,A)')&
                     '******** WARNING n^ ', ireset(3,icv),&
                     'for Scale Factor of CV n^ ',&
                     icv ,' ********'
                IF (paral%io_parent)&
                     WRITE(6,'(10x,A,I5,A,f14.6,A,f14.6/10x,A,f14.6/)')&
                     '    in the last ',ninter,' mstep MAX-MED = ',hel1,&
                     ', MED-MIN = ',hel2,'  whereas HLLW*SCF = ', fact
             ELSE
                ireset(3,icv) = 0
                ireset(2,icv) = 0
             ENDIF
             IF (ireset(2,icv) .GE. nevent) THEN
                scf_new = 0.2_real_8*cscl_fac(1,icv)+0.8_real_8*damp*elong/rmeta%hllw
                IF (scf_new .GT. cscl_fac(2,icv)) THEN
                   cscl_fac(1,icv) = scf_new
                   IF (paral%io_parent)&
                        WRITE(6,'(20x,A,I5,A,f14.6)')&
                        'Scaling Factor of CV',icv,&
                        ' reduced to ',cscl_fac(1,icv)
                   ireset(2,icv) = 0
                   ireset(3,icv) = 0
                ENDIF
             ELSEIF (ireset(3,icv) .GE. nevent) THEN
                scf_new = 0.2_real_8*cscl_fac(1,icv)+&
                     0.8_real_8*(hel1 + hel2)/((2.0_real_8-damp)*rmeta%hllw)*0.5_real_8
                IF (scf_new .LT. cscl_fac(3,icv)) THEN
                   ireset(2,icv) = 0
                   ireset(3,icv) = 0
                   cscl_fac(1,icv) = scf_new
                   IF (paral%io_parent)&
                        WRITE(6,'(20x,A,I5,A,f14.6)')&
                        'Scaling Factor of CV',icv,' increased to ',&
                        cscl_fac(1,icv)
                ENDIF
             ENDIF
             ireset(1,icv) = 0
          ENDIF
       ENDIF
15     CONTINUE
    ENDDO
    ! ==--------------------------------------------------------------==

    RETURN
  END SUBROUTINE scf_tune
  ! ==--------------------------------------------------------------==
  ! ==================================================================
  SUBROUTINE meta_stress(taup)

    REAL(real_8)                             :: taup(:,:,:)

    INTEGER                                  :: i, i1, i2, ia, is, j, j1, j2, &
                                                kk
    REAL(real_8)                             :: aa1(3), aa2(3), aa3(3), f2, &
                                                fv(3,3), temp(3,3)

    CALL zeroing(dstrmeta)!,6)

    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          dstrmeta(1) = dstrmeta(1) +&
               fhills(alpha(1),ia,is)*taup(beta(1),ia,is)
          dstrmeta(2) = dstrmeta(2) +&
               fhills(alpha(2),ia,is)*taup(beta(2),ia,is)
          dstrmeta(3) = dstrmeta(3) +&
               fhills(alpha(3),ia,is)*taup(beta(3),ia,is)
          dstrmeta(4) = dstrmeta(4) +&
               fhills(alpha(4),ia,is)*taup(beta(4),ia,is)
          dstrmeta(5) = dstrmeta(5) +&
               fhills(alpha(5),ia,is)*taup(beta(5),ia,is)
          dstrmeta(6) = dstrmeta(6) +&
               fhills(alpha(6),ia,is)*taup(beta(6),ia,is)
       ENDDO
    ENDDO

    CALL zeroing(fvbound)!,6)
    CALL zeroing(fv)!,6)
    IF (tvolbound) THEN
       DO i = 1,3
          aa1(i) = metr_com%ht(1,i)
          aa2(i) = metr_com%ht(2,i)
          aa3(i) = metr_com%ht(3,i)
       ENDDO
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

       i1 = 2
       i2 = 3
       j1 = i1
       j2 = i2

       DO i = 1,3
          DO j = 1,3
             fv(i,j) = f2*(metr_com%ht(i1,j1)*metr_com%ht(i2,j2)-metr_com%ht(i2,j1)*metr_com%ht(i1,j2))
             j1 = j2
             j2 = j
          ENDDO
          j1 = 2
          j2 = 3
          i1 = i2
          i2 = i
       ENDDO
       CALL dgemm('T','N',3,3,3,1._real_8,fv,3,metr_com%ht,3,0._real_8,temp,3)
       DO kk = 1,6
          fvbound(kk) =  fv(alpha(kk),beta(kk))
       ENDDO

    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE meta_stress
  ! ==--------------------------------------------------------------==
  ! ==================================================================
  SUBROUTINE  cv_exlagr_out2(i_meta,ekinc,ekinp)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: i_meta
    REAL(real_8)                             :: ekinc, ekinp

    CHARACTER(*), PARAMETER                  :: procedureN = 'cv_exlagr_out2'
    CHARACTER(len=20), PARAMETER             :: file3 = 'istvar_mtd', &
                                                file7 = 'enevar_mtd', &
                                                file8 = 'spinpo_mtd'

    CHARACTER(len=10)                        :: chnum
    CHARACTER(len=100)                       :: lineform
    INTEGER                                  :: iaa, icv, iee, ierr, it, n
    INTEGER, SAVE                            :: ifirst = fo_verb
    LOGICAL                                  :: ferror
    REAL(real_8)                             :: t_ion
    REAL(real_8), ALLOCATABLE                :: dummy(:)
    REAL(real_8), ALLOCATABLE, SAVE          :: avcv(:)

! ==--------------------------------------------------------------==

    IF (ifirst.GT.0) THEN
       ALLOCATE(avcv(ncolvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(avcv)!,ncolvar)
       ! ==--------------------------------------------------------------==
       ferror=.FALSE.
       IF (paral%io_parent)&
            CALL fileopen(53,file3,fo_old,ferror)
       IF (ferror) THEN
          IF (paral%io_parent)&
               CALL fileopen(53,file3,fo_new,ferror)
       ELSE
          IF (paral%io_parent)&
               REWIND(53)
          ALLOCATE(dummy(ncolvar),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          DO it = 1,i_meta-1
             IF (paral%io_parent)&
                  READ(53,*) n,(dummy(icv),icv=1,ncolvar),&
                  (avcv(icv),icv=1,ncolvar)
          ENDDO
          DEALLOCATE(dummy,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
          DO icv=1,ncolvar
             avcv(icv) = avcv(icv)*REAL(i_meta-1,kind=real_8)
          ENDDO
       ENDIF
    ELSE
       IF (paral%io_parent)&
            CALL fileopen(53,file3,fo_app,ferror)
    ENDIF

    IF (paral%io_parent)&
         CALL fileopen(57,file7,fo_app+ifirst,ferror)
    IF (lmeta%tlocalizespin) THEN
       IF (paral%io_parent)&
            CALL fileopen(58,file8,fo_app+ifirst,ferror)
    ENDIF

    ifirst=0
    ! ==--------------------------------------------------------------==
    DO icv = 1,ncolvar
       avcv(icv) = avcv(icv)+cv_ist(icv)
    ENDDO

    ! Write output
    IF (paral%io_parent)&
         WRITE(chnum,'(I5)') ncolvar
    CALL xstring(chnum,iaa,iee)

    lineform =&
         '(1X,I7,'//chnum(iaa:iee)//'E14.6,'//&
         chnum(iaa:iee)//'E14.6)'
    IF (paral%io_parent)&
         WRITE(53,lineform) iteropt%nfi,(cv_ist(icv),icv=1,ncolvar),&
         (avcv(icv)/REAL(i_meta,kind=real_8),icv=1,ncolvar)

    ! Close files
    IF (paral%io_parent)&
         CALL fileclose(53)

    t_ion = ekinp*factem*2._real_8/glib
    lineform = '(I6,8E16.6)'
    IF (paral%io_parent)&
         WRITE(57,lineform) iteropt%nfi,&
         t_ion,ekinc,ekincv,vharm,rmeta%gausspot,&
         ener_com%etot,rmeta%eham_hill,rmeta%eham_hill+rmeta%gausspot
    IF (paral%io_parent)&
         CALL fileclose(57)

    IF (lmeta%tlocalizespin) THEN
       ! Print Position of the center of the spin density
       IF (paral%io_parent)&
            WRITE(58,'(I6,3f12.6)') iteropt%nfi, rcc(1),rcc(2),rcc(3)
       IF (paral%io_parent)&
            CALL fileclose(58)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE cv_exlagr_out2
  ! ==--------------------------------------------------------------==
  ! ==================================================================
  SUBROUTINE wmtdres(nw,ntot_iter,i_meta)

    INTEGER                                  :: nw, ntot_iter, i_meta

    CHARACTER(len=20)                        :: fformat
    INTEGER                                  :: ia, icv, ie, it
    LOGICAL                                  :: ferror

    IF (paral%parent) THEN
       CALL xstring(fmtdres,ia,ie)
       IF (paral%io_parent)&
            CALL fileopen(nw,fmtdres(ia:ie),fo_def+fo_ufo,ferror)
       IF (paral%io_parent)&
            REWIND(nw)
       IF (paral%io_parent)&
            WRITE(nw) i_meta-1
       IF (paral%io_parent)&
            WRITE(nw) ncolvar

       DO icv = 1,ncolvar
          IF (paral%io_parent)&
               WRITE(nw) (cv_path(it,icv),it=1,i_meta-1)
          IF (paral%io_parent)&
               WRITE(nw) (cscl_val(it,icv),it=1,i_meta-1)
       ENDDO

       IF (lmeta%lextlagrange) THEN

          IF (paral%io_parent)&
               WRITE(nw)(cv_dyn(icv)-cv_path(i_meta-1,icv),&
               icv=1,ncolvar)
          IF (paral%io_parent)&
               WRITE(nw)(cv_vel(icv),icv=1,ncolvar)

       ENDIF
       IF (paral%io_parent)&
            WRITE(nw) (hllw_val(it,1),it=1,i_meta-1)
       IF (paral%io_parent)&
            WRITE(nw) (hllh_val(it,1),it=1,i_meta-1)

       IF (paral%io_parent)&
            CALL fileclose(nw)
       IF (paral%io_parent)&
            WRITE(fformat,'(A,I2,A)') '(/,A,T',MAX(34,65-(ie-ia)),',A)'
       IF (paral%io_parent)&
            WRITE(6,fformat)&
            ' MTD RESTART INFO WRITTEN ON FILE ',fmtdres(ia:ie)

    ENDIF  ! parent
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE wmtdres

  ! ==================================================================
  SUBROUTINE rmtdres(nr,cv_disp,tfull,ntot_iter)

    INTEGER                                  :: nr
    REAL(real_8)                             :: cv_disp(ncolvar)
    LOGICAL                                  :: tfull
    INTEGER                                  :: ntot_iter

    CHARACTER(len=20)                        :: fformat
    INTEGER                                  :: ia, icv, ie, it, nc
    LOGICAL                                  :: ferror

    IF (paral%parent) THEN
       ! Construct the filename
       CALL xstring(fmtdres,ia,ie)
       ferror=.FALSE.
       IF (paral%io_parent)&
            CALL fileopen(nr,fmtdres(ia:ie),fo_old+fo_ufo,ferror)
       IF (ferror) THEN
          IF (paral%io_parent)&
               WRITE(fformat,'(A,I2,A)')&
               '(/,A,T',MAX(36,65-(ie-ia)),',A)'
          IF (paral%io_parent)&
               WRITE(6,fformat)&
               'RMTDRES| MTD RESTART FILE NOT FOUND:',fmtdres(ia:ie)
          CALL stopgm('RMTDRES',' FILE NOT FOUND ',& 
               __LINE__,__FILE__)
       ENDIF

       IF (paral%io_parent)&
            REWIND(nr)
       IF (paral%io_parent)&
            READ(nr) imeta%i_meta_res
       IF (paral%io_parent)&
            READ(nr) nc

       IF (nc .NE. ncolvar) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,I3,A,I3)') 'RMTDRES! DIFFERENT # of CV ',&
               nc,' vs',ncolvar
          CALL stopgm('RMTDRES','USE SAME # of CV',& 
               __LINE__,__FILE__)
       ENDIF

       IF (tfull) THEN

          DO icv = 1,ncolvar
             IF (paral%io_parent)&
                  READ(nr) (cv_path(it,icv),it=1,imeta%i_meta_res)
             IF (paral%io_parent)&
                  READ(nr) (cscl_val(it,icv),it=1,imeta%i_meta_res)
          ENDDO

          IF (lmeta%lextlagrange) THEN
             IF (paral%io_parent)&
                  READ(nr)(cv_disp(icv),icv=1,ncolvar)
             IF (paral%io_parent)&
                  READ(nr)(cv_vel(icv),icv=1,ncolvar)
          ENDIF

          IF (paral%io_parent)&
               READ(nr) (hllw_val(it,1),it=1,imeta%i_meta_res)
          IF (paral%io_parent)&
               READ(nr) (hllh_val(it,1),it=1,imeta%i_meta_res)
          IF (paral%io_parent)&
               WRITE(fformat,'(A,I2,A)') '(/,A,T',MAX(34,65-(ie-ia)),',A)'
          IF (paral%io_parent)&
               WRITE(6,fformat)&
               ' MTD RESTART INFO READ FROM FILE ',fmtdres(ia:ie)
       ENDIF
       IF (paral%io_parent)&
            CALL fileclose(nr)

    ENDIF  ! parent
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rmtdres

  ! ==================================================================
  SUBROUTINE  cv_read_out(cv_disp)

    REAL(real_8)                             :: cv_disp(ncolvar)

    CHARACTER(len=20), PARAMETER :: file1 = 'colvar_mtd', &
      file2 = 'parvar_mtd', file5 = 'disvar_mtd', file6 = 'velvar_mtd'

    INTEGER                                  :: icv, ii, ij
    LOGICAL                                  :: ferror
    REAL(real_8)                             :: Displacement

    ferror=.FALSE.
    IF (paral%io_parent)&
         CALL fileopen(51,file1,fo_old,ferror)
    IF (ferror) GOTO 20
    IF (paral%io_parent)&
         CALL fileopen(52,file2,fo_old,ferror)
    IF (ferror) GOTO 20
    IF (paral%io_parent)&
         CALL fileopen(55,file5,fo_old,ferror)
    IF (ferror) GOTO 20
    IF (paral%io_parent)&
         CALL fileopen(56,file6,fo_old,ferror)
    IF (ferror) GOTO 20
    IF (paral%io_parent)&
         REWIND(51)
    IF (paral%io_parent)&
         REWIND(52)
    IF (paral%io_parent)&
         REWIND(55)
    IF (paral%io_parent)&
         REWIND(56)
    DO ii = 1,imeta%i_meta_res
       IF (paral%io_parent)&
            READ(51,err=20,END=20,fmt=*)  ij,&
            (cv_path(ii,icv),icv=1,ncolvar),&
            (cscl_val(ii,icv),icv=1,ncolvar)
       ! from output forces and parameters
       IF (paral%io_parent)&
            READ(52,err=20,END=20,fmt=*) ij,&
            Displacement,hllw_val(ii,1),hllh_val(ii,1)
       IF (paral%io_parent)&
            READ(55,err=20,END=20,fmt=*) ij,&
            (cv_disp(icv),icv=1,ncolvar)
       IF (paral%io_parent)&
            READ(56,err=20,END=20,fmt=*)  ij,&
            (cv_vel(icv),icv=1,ncolvar)
    ENDDO
    IF (paral%io_parent)&
         CALL fileclose(51)
    IF (paral%io_parent)&
         CALL fileclose(52)
    IF (paral%io_parent)&
         CALL fileclose(55)
    IF (paral%io_parent)&
         CALL fileclose(56)

    GOTO 100
20  CONTINUE
    IF (paral%io_parent)&
         WRITE(6,*) ' ERROR WHILE READING RESTART DATA '
    CALL stopgm('META_EXTLAGR',' ',& 
         __LINE__,__FILE__)
100 CONTINUE

    RETURN
  END SUBROUTINE cv_read_out
  ! ==--------------------------------------------------------------==

END MODULE meta_exlagr_utils
