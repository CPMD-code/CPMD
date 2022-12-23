MODULE meta_ex_mul_util_utils
  USE cnst,                            ONLY: factem,&
                                             pi
  USE cnst_dyn,                        ONLY: &
       cscl_fac, cscl_val, cv_dyn, cv_ist, cv_path, cv_vel, ekincv, fmtdres, &
       hhm, hllh_val, hllw_val, hwm, iangcv, imeta, inter_hill, &
       inter_hill_max, kharm, lmeta, ncolvar, ncvsys, nsubsys, rmeta, vharm
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
  USE kinds,                           ONLY: real_8
  USE nose,                            ONLY: glib
  USE parac,                           ONLY: paral
  USE readsr_utils,                    ONLY: xstring
  USE ropt,                            ONLY: iteropt
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: calc_pos_mdyn
  PUBLIC :: cv_diff_md
  PUBLIC :: cv_exl_outm
  PUBLIC :: cv_exl_outm2
  PUBLIC :: rmtdresm
  PUBLIC :: wmtdresm
  PUBLIC :: cv_read_outm

CONTAINS

  ! ================================================================== 
  SUBROUTINE calc_pos_mdyn(tollm,i_cvst,i_temp,i_meta,hc_last,&
       hill_add)
    ! ==--------------------------------------------------------------==

    REAL(real_8)                             :: tollm
    INTEGER                                  :: i_cvst, i_temp, i_meta
    REAL(real_8)                             :: hc_last(ncolvar)
    LOGICAL                                  :: hill_add

    CHARACTER(*), PARAMETER                  :: procedureN = 'calc_pos_mdyn'

    INTEGER                                  :: icv, ierr, is, jcv, testc
    REAL(real_8)                             :: av_disp, diff, min_disp
    REAL(real_8), ALLOCATABLE                :: scr(:)

! ==--------------------------------------------------------------==

    ALLOCATE(scr(ncolvar),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! Store CV values and Check averages


    IF (i_cvst .EQ.  inter_hill ) THEN
       jcv = 0
       testc = 0
       DO is = 1,nsubsys
          min_disp = 100.0_real_8
          av_disp  = 0.0_real_8
          DO icv = 1,ncvsys(is)
             diff = (cv_dyn(jcv+icv)-hc_last(jcv+icv))
             IF (iangcv(icv) .EQ. 1 .AND. diff .GT. pi) THEN
                diff = diff -2.0_real_8*pi
             ELSEIF (iangcv(jcv+icv) .EQ. 1 .AND. diff .LT. -pi) THEN
                diff = diff +2.0_real_8*pi
             ENDIF
             diff = (diff/cscl_fac(1,jcv+icv))**2.0_real_8
             av_disp = av_disp + diff
             min_disp = MIN(diff,min_disp)
          ENDDO
          av_disp = SQRT(av_disp)
          IF (av_disp .GT. tollm .OR. i_meta .EQ. 1) THEN
             testc = testc+1
          ENDIF

          IF (paral%io_parent)&
               WRITE(6,'(4x,A,I4,A,f10.6,A,f10.6,A,f10.6/)')&
               'SYS ',is,'MIN DISP. = ',SQRT(min_disp),&
               '  ||DISP.tot|| = ', av_disp,&
               '  Tolerance ',tollm
          jcv = jcv + ncvsys(is)
       ENDDO
       IF (testc .EQ. nsubsys) THEN
          hill_add = .TRUE.
       ENDIF
    ELSEIF (i_cvst .GT.  inter_hill ) THEN

       i_temp = i_temp + 1
       testc = 0
       jcv = 0
       DO is = 1,nsubsys
          min_disp = 100.0_real_8
          DO icv = 1,ncvsys(is)
             diff = (cv_dyn(jcv+icv)-hc_last(jcv+icv))
             IF (iangcv(jcv+icv) .EQ. 1 .AND. diff .GT. pi) THEN
                diff = diff -2.0_real_8*pi
             ELSEIF (iangcv(jcv+icv) .EQ. 1 .AND. diff .LT. -pi) THEN
                diff = diff +2.0_real_8*pi
             ENDIF
             diff = (diff/cscl_fac(1,jcv+icv))**2.0_real_8
             scr(jcv+icv) = diff
             min_disp = MIN(diff,min_disp)
          ENDDO
          IF (MOD(i_temp,imeta%icheck) .EQ. 0 ) THEN
             av_disp = 0.0_real_8
             !$omp parallel do private(ICV) reduction(+:AV_DISP)
             DO icv = 1,ncvsys(is)
                av_disp = av_disp+scr(jcv+icv)
             ENDDO
             av_disp = SQRT(av_disp)
             IF (paral%io_parent)&
                  WRITE(6,'(4x,A,I4,A,f10.6,A,f10.6,A,f10.6/)')&
                  'SYS ',is,'MIN DISP. = ',SQRT(min_disp),&
                  '  ||DISP.tot|| = ',av_disp,&
                  '  Tolerance ',tollm
             IF (av_disp .GT. tollm .OR. i_cvst .GT. inter_hill_max) THEN
                testc = testc + 1
             ENDIF
          ENDIF
          jcv = jcv + ncvsys(is)
       ENDDO
       IF (testc .EQ. nsubsys) THEN
          hill_add = .TRUE.
       ENDIF
    ENDIF   ! I_CVST
    ! ==--------------------------------------------------------------==
    DEALLOCATE(scr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==

    RETURN
  END SUBROUTINE calc_pos_mdyn
  ! ==--------------------------------------------------------------==
  ! ==--------------------------------------------------------------==
  SUBROUTINE  cv_diff_md(cv_last,cvd_ist,cvd_scl,cv_scl,&
       i_meta,ntot_iter)


    REAL(real_8)                             :: cv_last(ncolvar), &
                                                cvd_ist(ncolvar), &
                                                cvd_scl(ncolvar), &
                                                cv_scl(ncolvar)
    INTEGER                                  :: i_meta, ntot_iter

    INTEGER                                  :: icv, ii
    REAL(real_8)                             :: diff

! ==--------------------------------------------------------------==
! Diffusion calculation

    CALL zeroing(cvd_ist)!,ncolvar)
    CALL zeroing(cvd_scl)!,ncolvar)

    DO icv = 1,ncolvar

       DO ii = 1,i_meta-1
          diff = ABS(cv_path(ii,icv)-cv_last(icv))
          IF (iangcv(icv) .EQ. 1 .AND. diff .GT. pi) THEN
             diff = diff - 2.0_real_8*pi
          ENDIF
          cvd_ist(icv) = cvd_ist(icv) +&
               diff/REAL(i_meta-ii,kind=real_8)
          cvd_scl(icv) = cvd_scl(icv) +&
               diff/REAL(i_meta-ii,kind=real_8)/cscl_val(ii,icv)
       ENDDO
       cv_scl(icv) = cv_last(icv) / cscl_fac(1,icv)
    ENDDO
    !$omp parallel do private(ICV)
    DO icv  = 1,ncolvar
       cv_path(i_meta,icv) =  cv_last(icv)
       cscl_val(i_meta,icv)=  cscl_fac(1,icv)
    ENDDO
    !$omp parallel do private(II)
    DO ii = 1,nsubsys
       hllw_val(i_meta,ii)    =  hwm(ii)
       hllh_val(i_meta,ii)    =  hhm(ii)
    ENDDO

    ! ==--------------------------------------------------------------==

    RETURN
  END SUBROUTINE cv_diff_md
  ! ==--------------------------------------------------------------==
  ! ==--------------------------------------------------------------==

  SUBROUTINE  cv_exl_outm(i_meta,hc_last,f_harm,f_hill,f_wall,&
       cvd_ist,disp_sys,ekinc,ekinp)


    ! ==--------------------------------------------------------------==

    INTEGER                                  :: i_meta
    REAL(real_8) :: hc_last(ncolvar), f_harm(ncolvar), f_hill(ncolvar), &
      f_wall(ncolvar), cvd_ist(ncolvar), disp_sys(nsubsys), ekinc, ekinp

    CHARACTER(len=10), PARAMETER :: file1 = 'colvar_mtd', &
      file2 = 'parvar_mtd', file3 = 'istvar_mtd', file4 = 'forfac_mtd', &
      file5 = 'disvar_mtd', file6 = 'velvar_mtd', file7 = 'enevar_mtd'

    CHARACTER(len=10)                        :: chnum, flag
    CHARACTER(len=100) :: lineform, outcolvar(10), outdisp(10), &
      outforfac(10), outistcolvar(10), outparvar(10), outvel(10)
    INTEGER                                  :: ia, iaa, icv, ie, iee, is, jj
    INTEGER, SAVE                            :: ifirst = fo_verb
    LOGICAL                                  :: ferror
    REAL(real_8)                             :: t_ion

! ==--------------------------------------------------------------==

    jj = 0
    ! Path of output files
    DO is = 1,nsubsys
       IF (paral%io_parent)&
            WRITE(flag,'(I5)') is
       CALL xstring(flag,ia,ie)
       outcolvar(is) = file1//'_S'//flag(ia:ie)
       outparvar(is) = file2//'_S'//flag(ia:ie)
       outistcolvar(is)=file3//'_S'//flag(ia:ie)
       outforfac(is) = file4//'_S'//flag(ia:ie)
       outdisp(is)   = file5//'_S'//flag(ia:ie)
       outvel(is)    = file6//'_S'//flag(ia:ie)

       ! ==--------------------------------------------------------------==
       ! Open output files
       IF (paral%io_parent)&
            CALL fileopen(51,outcolvar(is),   fo_app+ifirst,ferror)
       IF (paral%io_parent)&
            CALL fileopen(52,outparvar(is),   fo_app+ifirst,ferror)
       IF (paral%io_parent)&
            CALL fileopen(53,outistcolvar(is),fo_app+ifirst,ferror)
       IF (paral%io_parent)&
            CALL fileopen(54,outforfac(is),   fo_app+ifirst,ferror)
       IF (paral%io_parent)&
            CALL fileopen(55,outdisp(is),     fo_app+ifirst,ferror)
       IF (paral%io_parent)&
            CALL fileopen(56,outvel(is),      fo_app+ifirst,ferror)
       ! ==--------------------------------------------------------------==
       ! Write output
       IF (paral%io_parent)&
            WRITE(chnum,'(I5)') ncvsys(is)
       CALL xstring(chnum,iaa,iee)
       lineform =&
            '(1X,I7,'//chnum(iaa:iee)//'E14.6,'//chnum(iaa:iee)//'E13.4)'
       IF (paral%io_parent)&
            WRITE(51,lineform) iteropt%nfi,(hc_last(jj+icv),icv=1,ncvsys(is)),&
            (cscl_fac(1,jj+icv),icv=1,ncvsys(is))

       lineform =&
            '(1X,I7,3F14.6)'
       IF (paral%io_parent)&
            WRITE(52,lineform) iteropt%nfi,&
            disp_sys(is),hwm(is),hhm(is)
       ! Close files
       IF (paral%io_parent)&
            CALL fileclose(51)
       IF (paral%io_parent)&
            CALL fileclose(52)


       lineform =&
            '(1X,I7,'//chnum(iaa:iee)//'E14.6,'//chnum(iaa:iee)//'E14.6)'
       IF (paral%io_parent)&
            WRITE(53,lineform) iteropt%nfi,(cv_ist(jj+icv),icv=1,ncvsys(is)),&
            (cv_ist(jj+icv)-hc_last(jj+icv),icv=1,ncvsys(is))
       ! Close files
       IF (paral%io_parent)&
            CALL fileclose(53)

       ! Print force factors
       lineform =&
            '(1X,I7,'//chnum(iaa:iee)//'E14.6,'&
            //chnum(iaa:iee)//'E14.6,'//chnum(iaa:iee)//'E16.4)'
       IF (paral%io_parent)&
            WRITE(54,lineform) iteropt%nfi,(f_harm(jj+icv),icv=1,ncvsys(is)),&
            (f_hill(jj+icv),icv=1,ncvsys(is)),&
            (f_wall(jj+icv),icv=1,ncvsys(is))

       ! Close files
       IF (paral%io_parent)&
            CALL fileclose(54)

       ! Print last displacements
       lineform =&
            '(I6,'//chnum(iaa:iee)//'f11.6,'//chnum(iaa:iee)//&
            'f10.6,'//chnum(iaa:iee)//'f10.6)'
       IF (paral%io_parent)&
            WRITE(55,lineform) iteropt%nfi,&
            ((cv_dyn(jj+icv)-hc_last(jj+icv)),icv=1,ncvsys(is)),&
            (cvd_ist(jj+icv),icv=1,ncvsys(is)),&
            (kharm(jj+icv),icv=1,ncvsys(is))

       ! Close files
       IF (paral%io_parent)&
            CALL fileclose(55)

       ! Print last velocities, temperature and kinetic energy
       ! EHAM_HILL=EKINP+ETOT+ENOSE+ENOSP+ECNSTR+EKINC
       lineform =&
            '(I6,'//chnum(iaa:iee)//'f11.6)'
       IF (paral%io_parent)&
            WRITE(56,lineform) iteropt%nfi,(cv_vel(jj+icv),icv=1,ncvsys(is))

       IF (paral%io_parent)&
            CALL fileclose(56)

       jj = jj + ncvsys(is)
    ENDDO

    IF (paral%io_parent)&
         CALL fileopen(57,file7,fo_app+ifirst,ferror)
    ifirst=0

    t_ion = ekinp*factem*2._real_8/glib
    IF (paral%io_parent)&
         WRITE(57,'(I6,8E16.6)') iteropt%nfi,t_ion,ekinc,&
         ekincv,vharm,rmeta%gausspot,&
         ener_com%etot,rmeta%eham_hill,rmeta%eham_hill+rmeta%gausspot
    IF (paral%io_parent)&
         CALL fileclose(57)
    ! ==--------------------------------------------------------------==

    RETURN
  END SUBROUTINE cv_exl_outm
  ! ==--------------------------------------------------------------==
  ! ==--------------------------------------------------------------==

  SUBROUTINE  cv_exl_outm2(i_meta,ekinc,ekinp)
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: i_meta
    REAL(real_8)                             :: ekinc, ekinp

    CHARACTER(*), PARAMETER                  :: procedureN = 'cv_exl_outm2'
    CHARACTER(len=10), PARAMETER             :: file3 = 'istvar_mtd', &
                                                file7 = 'enevar_mtd'

    CHARACTER(len=10)                        :: chnum, flag
    CHARACTER(len=100)                       :: lineform, outistcolvar(10)
    INTEGER                                  :: ia, iaa, icv, ie, iee, ierr, &
                                                is, it, jj, n
    INTEGER, SAVE                            :: ifirst = fo_verb
    LOGICAL                                  :: ferror
    REAL(real_8)                             :: dummy, t_ion
    REAL(real_8), ALLOCATABLE, SAVE          :: avcv(:)

    IF (ifirst.NE.0) THEN
       ALLOCATE(avcv(ncolvar),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(avcv)!,ncolvar)
    ENDIF

    jj = 0
    ! Path of output files
    DO is = 1,nsubsys
       IF (paral%io_parent)&
            WRITE(flag,'(I5)') is
       CALL xstring(flag,ia,ie)
       outistcolvar(is)=file3//'_S'//flag(ia:ie)

       ! ==--------------------------------------------------------------==
       ! Open output file
       IF (ifirst.NE.0) THEN
          ferror=.FALSE.
          IF (paral%io_parent)&
               CALL fileopen(53,outistcolvar(is),fo_old,ferror)
          IF (ferror) THEN
             IF (paral%io_parent)&
                  CALL fileopen(53,outistcolvar(is),fo_new,ferror)
          ELSE
             IF (paral%io_parent)&
                  REWIND(53)
             DO it = 1,i_meta-1
                IF (paral%io_parent)&
                     READ(53,*) n,(dummy,icv=1,ncolvar),&
                     (avcv(jj+icv),icv=1,ncvsys(is))
             ENDDO
             DO icv=1,ncvsys(is)
                avcv(jj+icv) = avcv(jj+icv)*REAL(i_meta-1,kind=real_8)
             ENDDO
          ENDIF
       ELSE
          IF (paral%io_parent)&
               CALL fileopen(53,outistcolvar(is),fo_app,ferror)
       ENDIF
       ! ==--------------------------------------------------------------==

       DO icv = 1,ncvsys(is)
          avcv(jj+icv) = avcv(jj+icv)+cv_ist(jj+icv)
       ENDDO

       ! Write output
       IF (paral%io_parent)&
            WRITE(chnum,'(I5)') ncvsys(is)
       CALL xstring(chnum,iaa,iee)

       lineform =&
            '(1X,I7,'//chnum(iaa:iee)//'E14.6,'//chnum(iaa:iee)//'E14.6)'
       IF (paral%io_parent)&
            WRITE(53,lineform) iteropt%nfi,(cv_ist(jj+icv),icv=1,ncvsys(is)),&
            (avcv(jj+icv)/REAL(i_meta,kind=real_8),icv=1,ncvsys(is))
       ! Close files
       IF (paral%io_parent)&
            CALL fileclose(53)
       jj = jj + ncvsys(is)
    ENDDO

    IF (paral%io_parent)&
         CALL fileopen(57,file7,fo_app+ifirst,ferror)
    t_ion = ekinp*factem*2._real_8/glib
    IF (paral%io_parent)&
         WRITE(57,'(I6,8E16.6)') iteropt%nfi,t_ion,ekinc,&
         ekincv,vharm,rmeta%gausspot,&
         ener_com%etot,rmeta%eham_hill,rmeta%eham_hill+rmeta%gausspot
    IF (paral%io_parent)&
         CALL fileclose(57)
    ifirst=0
    ! ==--------------------------------------------------------------==

    RETURN
  END SUBROUTINE cv_exl_outm2
  ! ==--------------------------------------------------------------==
  SUBROUTINE rmtdresm(nr,cv_disp,tfull,ntot_iter)


    INTEGER                                  :: nr
    REAL(real_8)                             :: cv_disp(ncolvar)
    LOGICAL                                  :: tfull
    INTEGER                                  :: ntot_iter

    CHARACTER(len=20)                        :: fformat
    INTEGER                                  :: ia, icv, ie, is, it, nc, ns
    LOGICAL                                  :: ferror

    IF (paral%parent) THEN
       CALL xstring(fmtdres,ia,ie)
       IF (paral%io_parent)&
            CALL fileopen(nr,fmtdres,fo_old+fo_ufo,ferror)
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

       IF (paral%io_parent)&
            READ(nr) ns
       IF (ns .NE. nsubsys) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A,I3,A,I3)')&
               'RMTDRES! DIFFERENT # of SUBSETS fur MULTI_MTD  ',&
               ns,' vs', nsubsys
          CALL stopgm('RMTDRES','USE SAME # of SETS',& 
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
               READ(nr) (ncvsys(is),is=1,nsubsys)
          DO is=1,nsubsys
             IF (paral%io_parent)&
                  READ(nr) (hllw_val(it,is),it=1,imeta%i_meta_res)
             IF (paral%io_parent)&
                  READ(nr) (hllh_val(it,is),it=1,imeta%i_meta_res)
          ENDDO

          IF (paral%io_parent)&
               WRITE(fformat,'(A,I2,A)') '(/,A,T',MAX(34,65-(ie-ia)),',A)'
          IF (paral%io_parent)&
               WRITE(6,fformat)&
               ' MTD RESTART INFO READ FROM FILE ',fmtdres(ia:ie)
       ENDIF
       IF (paral%io_parent)&
            CALL fileclose(nr)

    ENDIF  ! parent

    RETURN
  END SUBROUTINE rmtdresm

  ! ==--------------------------------------------------------------==
  SUBROUTINE wmtdresm(nw,ntot_iter,i_meta)


    INTEGER                                  :: nw, ntot_iter, i_meta

    CHARACTER(len=20)                        :: fformat
    INTEGER                                  :: ia, icv, ie, is, it
    LOGICAL                                  :: ferror

! Variables
! ==--------------------------------------------------------------==
! Construct the filename

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
            WRITE(nw) nsubsys
       IF (paral%io_parent)&
            WRITE(nw) (ncvsys(is),is=1,nsubsys)
       DO is=1,nsubsys
          IF (paral%io_parent)&
               WRITE(nw) (hllw_val(it,is),it=1,i_meta-1)
          IF (paral%io_parent)&
               WRITE(nw) (hllh_val(it,is),it=1,i_meta-1)
       ENDDO

       IF (paral%io_parent)&
            CALL fileclose(nw)
       IF (paral%io_parent)&
            WRITE(fformat,'(A,I2,A)') '(/,A,T',MAX(34,65-(ie-ia)),',A)'
       IF (paral%io_parent)&
            WRITE(6,fformat)&
            ' MTD RESTART INFO WRITTEN ON FILE ',fmtdres(ia:ie)

    ENDIF  ! parent

    RETURN
  END SUBROUTINE wmtdresm

  ! ==--------------------------------------------------------------==
  SUBROUTINE  cv_read_outm(cv_disp)

    REAL(real_8)                             :: cv_disp(ncolvar)

    CHARACTER(len=10), PARAMETER :: file1 = 'colvar_mtd', &
      file2 = 'parvar_mtd', file5 = 'disvar_mtd', file6 = 'velvar_mtd'

    CHARACTER(len=100)                       :: outcolvar(10), outdisp(10), &
                                                outparvar(10), outvel(10)
    CHARACTER(len=20)                        :: flag
    INTEGER                                  :: i, ia, icv, ie, ii, ij, is, jj
    LOGICAL                                  :: ferror
    REAL(real_8)                             :: displacement

! ==--------------------------------------------------------------==
! Path of output files

    DO is = 1,nsubsys
       IF (paral%io_parent)&
            WRITE(flag,'(I5)') is
       CALL xstring(flag,ia,ie)
       outcolvar(is)  = file1//'_S'//flag(ia:ie)
       outparvar(is)  = file2//'_S'//flag(ia:ie)
       outdisp(is)    = file5//'_S'//flag(ia:ie)
       outvel(is)     = file6//'_S'//flag(ia:ie)
    ENDDO

    ferror=.FALSE.
    jj = 0
    DO i = 1, nsubsys
       IF (paral%io_parent)&
            CALL fileopen(51,outcolvar(i),fo_old,ferror)
       IF (ferror) GOTO 20
       IF (paral%io_parent)&
            CALL fileopen(52,outparvar(i),fo_old,ferror)
       IF (ferror) GOTO 20
       IF (paral%io_parent)&
            CALL fileopen(55,outdisp(i),  fo_old,ferror)
       IF (ferror) GOTO 20
       IF (paral%io_parent)&
            CALL fileopen(56,outvel(i),   fo_old,ferror)
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
               (cv_path(ii,jj+icv),icv=1,ncvsys(i)),&
               (cscl_val(ii,jj+icv),icv=1,ncvsys(i))
          ! from output forces and parameters
          IF (paral%io_parent)&
               READ(52,err=20,END=20,fmt=*) ij,&
               displacement,hllw_val(ii,i),hllh_val(ii,i)
          IF (paral%io_parent)&
               READ(55,err=20,END=20,fmt=*) ij,&
               (cv_disp(jj+icv),icv=1,ncvsys(i))
          IF (paral%io_parent)&
               READ(56,err=20,END=20,fmt=*)  ij,&
               (cv_vel(jj+icv),icv=1,ncvsys(i))
       ENDDO
       IF (paral%io_parent)&
            CALL fileclose(51)
       IF (paral%io_parent)&
            CALL fileclose(52)
       IF (paral%io_parent)&
            CALL fileclose(55)
       IF (paral%io_parent)&
            CALL fileclose(56)
       jj = jj+ncvsys(i)

       GOTO 100
20     CONTINUE
       IF (paral%io_parent)&
            WRITE(6,*) ' ERROR WHILE READING RESTART DATA '
       CALL stopgm('META_EXTLAGR',' ',& 
            __LINE__,__FILE__)
100    CONTINUE
    ENDDO

    RETURN
  END SUBROUTINE cv_read_outm

  ! ==--------------------------------------------------------------==
  ! ==--------------------------------------------------------------==


END MODULE meta_ex_mul_util_utils
