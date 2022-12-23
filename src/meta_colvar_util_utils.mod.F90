MODULE meta_colvar_util_utils
  USE cnst,                            ONLY: factem,&
                                             pi
  USE cnst_dyn,                        ONLY: &
       cscl_fac, cscl_val, cv_ist, cv_path, det_colvar, hllh_val, hllw_val, &
       iangcv, imeta, inter_hill, inter_hill_max, ncolvar, rmeta
  USE cotr,                            ONLY: cotc0,&
                                             duat
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_app,&
                                             fo_verb
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE nose,                            ONLY: glib
  USE parac,                           ONLY: paral
  USE puttau_utils,                    ONLY: puttau
  USE readsr_utils,                    ONLY: xstring
  USE ropt,                            ONLY: iteropt
  USE system,                          ONLY: maxsys
  USE tpar,                            ONLY: dt_ions,&
                                             dtb2mi
  USE utils,                           ONLY: invmat
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: calc_aver
  PUBLIC :: cv_write_out
  PUBLIC :: cv_diffusion
  PUBLIC :: calc_pos
  PUBLIC :: tune_height
  PUBLIC :: fillc2

CONTAINS

  ! ==================================================================
  SUBROUTINE calc_aver(tollm,i_cvst,i_temp,i_meta,cv_store,cv_last,&
       cv_aver,hill_add)
    ! ==--------------------------------------------------------------==

    REAL(real_8)                             :: tollm
    INTEGER                                  :: i_cvst, i_temp, i_meta
    REAL(real_8)                             :: cv_store(inter_hill,ncolvar), &
                                                cv_last(ncolvar), &
                                                cv_aver(ncolvar)
    LOGICAL                                  :: hill_add

    CHARACTER(*), PARAMETER                  :: procedureN = 'calc_aver'

    INTEGER                                  :: icv, ierr, ii, ij
    REAL(real_8)                             :: av_disp, diff, min_disp, sum
    REAL(real_8), ALLOCATABLE                :: scr(:)

    ALLOCATE(scr(ncolvar),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! Store CV values and Check averages

    IF (i_cvst .LT.  inter_hill ) THEN
       ii = i_cvst
       !$omp parallel do private(ICV)
       DO icv = 1,ncolvar
          cv_store(ii,icv) = cv_ist(icv)
       ENDDO
    ELSEIF (i_cvst .EQ.  inter_hill ) THEN
       min_disp = 100.0_real_8
       ii = inter_hill

       av_disp  = 0.0_real_8
       DO icv = 1,ncolvar
          sum = 0.0_real_8
          DO ij = 1,inter_hill-1
             sum = sum + cv_store(ij,icv)
          ENDDO
          cv_store(ii,icv) = cv_ist(icv)
          sum = sum + cv_store(ii,icv)
          cv_aver(icv) = sum/REAL(inter_hill,kind=real_8)
          diff = ((cv_aver(icv)-cv_last(icv))/cscl_fac(1,icv))**2.0_real_8

          av_disp = av_disp + diff
          min_disp = MIN(diff,min_disp)
       ENDDO
       av_disp = SQRT(av_disp)
       IF (av_disp .GT. tollm .OR. i_meta .EQ. 1) hill_add = .TRUE.

       IF (paral%io_parent)&
            WRITE(6,'(/10x,A,f10.6,A,f10.6,A,f10.6/)')&
            'MIN DISP. = ',SQRT(min_disp),'  ||DISP.tot|| = ', av_disp,&
            '  Tolerance ',tollm
    ELSEIF (i_cvst .GT.  inter_hill ) THEN

       i_temp = i_temp + 1
       ii = inter_hill
       min_disp = 100.0_real_8
       DO icv = 1,ncolvar
          sum = 0.0_real_8
          DO ij = 1,inter_hill-1
             cv_store(ij,icv) = cv_store(ij+1,icv)
             sum = sum + cv_store(ij,icv)
          ENDDO
          cv_store(ii,icv) = cv_ist(icv)
          sum = sum + cv_store(ii,icv)
          cv_aver(icv) = sum/REAL(inter_hill,kind=real_8)
          diff = ((cv_aver(icv)-cv_last(icv))/cscl_fac(1,icv))**2
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
               WRITE(6,'(/10x,A,f10.6,A,f10.6,A,f10.6/)')&
               'MIN DISP. = ',SQRT(min_disp),'  ||DISP.tot|| = ',av_disp,&
               '  Tolerance ',tollm
          IF (av_disp .GT. tollm .OR. i_cvst .GT. inter_hill_max) THEN
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
  END SUBROUTINE calc_aver
  ! ==--------------------------------------------------------------==
  ! ==--------------------------------------------------------------==
  SUBROUTINE  cv_write_out(i_meta,cv_last,cv_aver,f1_cv,f2_cv,&
       cvd_ist,cv_scl,cvd_scl,displacement,ekinc,ekinp)

    ! ==--------------------------------------------------------------==

    INTEGER                                  :: i_meta
    REAL(real_8) :: cv_last(ncolvar), cv_aver(ncolvar), f1_cv(ncolvar), &
      f2_cv(ncolvar), cvd_ist(ncolvar), cv_scl(ncolvar), cvd_scl(ncolvar), &
      displacement, ekinc, ekinp

    CHARACTER(len=20), PARAMETER :: file1 = 'colvar_mtd', &
      file2 = 'parvar_mtd', file3 = 'scavar_mtd', file4 = 'forfac_mtd', &
      file5 = 'disvar_mtd', file6 = 'enevar_mtd'

    CHARACTER(len=10)                        :: chnum
    CHARACTER(len=100)                       :: lineform
    INTEGER                                  :: iaa, icv, iee
    INTEGER, SAVE                            :: ifirst = fo_verb
    LOGICAL                                  :: ferror
    REAL(real_8)                             :: t_ion

! Open output files

    IF (paral%io_parent) THEN
       CALL fileopen(51,file1,fo_app+ifirst,ferror)
       CALL fileopen(52,file2,fo_app+ifirst,ferror)
       CALL fileopen(53,file3,fo_app+ifirst,ferror)
       CALL fileopen(54,file4,fo_app+ifirst,ferror)
       CALL fileopen(55,file5,fo_app+ifirst,ferror)
       CALL fileopen(56,file6,fo_app+ifirst,ferror)
    END IF
    ! Verbose open only on the first run
    ifirst=0

    ! Write output
    IF (paral%io_parent)&
         WRITE(chnum,'(I5)') ncolvar
    CALL xstring(chnum,iaa,iee)
    lineform =&
         '(1X,I7,'//chnum(iaa:iee)//'E14.6,'//chnum(iaa:iee)//'E13.4)'
    IF (paral%io_parent)&
         WRITE(51,lineform) iteropt%nfi,(cv_last(icv),icv=1,ncolvar),&
         (cscl_fac(1,icv),icv=1,ncolvar)

    lineform =&
         '(1X,I7,3F14.6)'
    IF (paral%io_parent)&
         WRITE(52,lineform) iteropt%nfi,&
         displacement,rmeta%hllw,rmeta%hllh
    ! Close files
    IF (paral%io_parent)&
         CALL fileclose(51)
    IF (paral%io_parent)&
         CALL fileclose(52)


    lineform =&
         '(1X,I7,'//chnum(iaa:iee)//'E14.6,'//chnum(iaa:iee)//'E14.6)'
    IF (paral%io_parent)&
         WRITE(53,lineform) iteropt%nfi,(cv_scl(icv),icv=1,ncolvar),&
         (cvd_scl(icv),icv=1,ncolvar)
    ! Close files
    IF (paral%io_parent)&
         CALL fileclose(53)

    ! Print force factors
    lineform =&
         '(1X,I7,'//chnum(iaa:iee)//'E14.6,'//chnum(iaa:iee)//'E14.6)'
    IF (paral%io_parent)&
         WRITE(54,lineform) iteropt%nfi,(f1_cv(icv),icv=1,ncolvar),&
         (f2_cv(icv),icv=1,ncolvar)

    ! Close files
    IF (paral%io_parent)&
         CALL fileclose(54)

    ! Print last displacements
    lineform =&
         '(I6,'//chnum(iaa:iee)//'f11.6,'//chnum(iaa:iee)//&
         'f10.6)'
    IF (paral%io_parent)&
         WRITE(55,lineform) iteropt%nfi,&
         ((cv_aver(icv)-cv_last(icv)),icv=1,ncolvar),&
         (cvd_ist(icv),icv=1,ncolvar)

    t_ion = ekinp*factem*2._real_8/glib
    lineform = '(I6,5E14.5)'
    IF (paral%io_parent)&
         WRITE(56,lineform) iteropt%nfi,&
         t_ion,ekinc,rmeta%gausspot,&
         ener_com%etot,rmeta%eham_hill+rmeta%gausspot
    ! Close files
    IF (paral%io_parent)&
         CALL fileclose(55)
    IF (paral%io_parent)&
         CALL fileclose(56)

    ! ==--------------------------------------------------------------==

    RETURN
  END SUBROUTINE cv_write_out
  ! ==--------------------------------------------------------------==
  ! ==--------------------------------------------------------------==
  SUBROUTINE  cv_diffusion(cv_last,cvd_ist,cvd_scl,cv_scl,&
       i_meta,ntot_iter)


    REAL(real_8)                             :: cv_last(ncolvar), &
                                                cvd_ist(ncolvar), &
                                                cvd_scl(ncolvar), &
                                                cv_scl(ncolvar)
    INTEGER                                  :: i_meta, ntot_iter

    INTEGER                                  :: icv, ii
    REAL(real_8)                             :: diff

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
    hllw_val(i_meta,1)    =  rmeta%hllw
    hllh_val(i_meta,1)    =  rmeta%hllh

    ! ==--------------------------------------------------------------==

    RETURN
  END SUBROUTINE cv_diffusion
  ! ==--------------------------------------------------------------==
  ! ==--------------------------------------------------------------==

  SUBROUTINE calc_pos(tollm,i_cvst,i_temp,i_meta,cv_last,cv_aver,&
       hill_add)


    REAL(real_8)                             :: tollm
    INTEGER                                  :: i_cvst, i_temp, i_meta
    REAL(real_8)                             :: cv_last(ncolvar), &
                                                cv_aver(ncolvar)
    LOGICAL                                  :: hill_add

    CHARACTER(*), PARAMETER                  :: procedureN = 'calc_pos'

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
          cv_aver(icv) = cv_ist(icv)
          diff = ((cv_ist(icv)-cv_last(icv))/cscl_fac(1,icv))**2.0_real_8
          av_disp = av_disp + diff
          min_disp = MIN(diff,min_disp)
       ENDDO
       av_disp = SQRT(av_disp)
       IF (av_disp .GT. tollm .OR. i_meta .EQ. 1) hill_add = .TRUE.

       IF (paral%io_parent)&
            WRITE(6,'(/10x,A,f10.6,A,f10.6,A,f10.6/)')&
            'MIN DISP. = ',SQRT(min_disp),'  ||DISP.tot|| = ', av_disp,&
            '  Tolerance ',tollm
    ELSEIF (i_cvst .GT.  inter_hill ) THEN

       i_temp = i_temp + 1
       ii = inter_hill
       min_disp = 100.0_real_8
       DO icv = 1,ncolvar
          cv_aver(icv) = cv_ist(icv)
          diff = ((cv_aver(icv)-cv_last(icv))/cscl_fac(1,icv))**2
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
               WRITE(6,'(/10x,A,f10.6,A,f10.6,A,f10.6/)')&
               'MIN DISP. = ',SQRT(min_disp),'  ||DISP.tot|| = ',av_disp,&
               '  Tolerance ',tollm
          IF (av_disp .GT. tollm .OR. i_cvst .GT. inter_hill_max) THEN
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
  END SUBROUTINE calc_pos
  ! ==--------------------------------------------------------------==
  ! ==--------------------------------------------------------------==
  SUBROUTINE tune_height(f_aver,velp,fion,tscr)

    REAL(real_8)                             :: f_aver(ncolvar), velp(:,:,:), &
                                                fion(:,:,:), tscr(:,:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'tune_height'

    INTEGER                                  :: ia, icv, idof, ierr, info, &
                                                is, jcv, l_deltar, l_hinv, &
                                                l_tmat, length
    REAL(real_8)                             :: fact
    REAL(real_8), ALLOCATABLE                :: bmat(:,:), deltar(:), &
                                                hinv(:,:), tmat(:,:)

! ==--------------------------------------------------------------==

    l_deltar = cotc0%nodim
    l_tmat   = ncolvar*ncolvar
    l_hinv   = ncolvar*cotc0%nodim

    length   = l_deltar + 2 * l_tmat  + l_hinv



    ! TODO align for BG
    ALLOCATE(deltar(cotc0%nodim),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(tmat(ncolvar, ncolvar),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(bmat(ncolvar, ncolvar),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(hinv(ncolvar, cotc0%nodim),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    !$omp parallel do private(IS,IA,FACT)
    DO is=1,ions1%nsp
       fact = 1._real_8/(dt_ions*dtb2mi(is))
       DO ia=1,ions0%na(is)
          tscr(1,ia,is)  =fact*(fion(1,ia,is)*dtb2mi(is)+velp(1,ia,is))
          tscr(2,ia,is)  =fact*(fion(2,ia,is)*dtb2mi(is)+velp(2,ia,is))
          tscr(3,ia,is)  =fact*(fion(3,ia,is)*dtb2mi(is)+velp(3,ia,is))
       ENDDO
    ENDDO
    CALL puttau(tscr,deltar)

    CALL zeroing(tmat)!,ncolvar*ncolvar)
    !$omp parallel do private(JCV,ICV,IDOF)
    DO jcv = 1,ncolvar
       DO icv = 1,ncolvar
          DO idof = 1,cotc0%nodim
             tmat(icv,jcv) = tmat(icv,jcv) +&
                  det_colvar(idof,icv)*det_colvar(idof,jcv)
          ENDDO
       ENDDO
    ENDDO

    CALL invmat(ncolvar,tmat,bmat,info)

    ! HINV(NCOLVAR,NODIM) is the inverse of DET_COLVAR(NDIM,NCOLVAR)
    CALL zeroing(hinv)!,cotc0%nodim*ncolvar)
    !$omp parallel do private(IDOF,ICV,JCV)
    DO idof = 1,cotc0%nodim
       DO icv = 1,ncolvar
          DO jcv = 1,ncolvar
             hinv(icv,idof) = hinv(icv,idof) +&
                  tmat(icv,jcv)*det_colvar(idof,jcv)
          ENDDO
       ENDDO
    ENDDO
    !$omp parallel do private(ICV,IDOF)
    DO icv = 1,ncolvar
       DO idof = 1,cotc0%nodim
          f_aver(icv) = f_aver(icv) + hinv(icv,idof)*deltar(idof)
       ENDDO
    ENDDO

    ! ==--------------------------------------------------------------==
    DEALLOCATE(deltar,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(tmat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(bmat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(hinv,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE tune_height
  ! ==--------------------------------------------------------------==

  SUBROUTINE fillc2(iat,tau,x,d2_a,tran)
    ! ==--------------------------------------------------------------==
    ! == Extract the coordinates X(1:3) of IAT index in TAU (TSCR)    ==
    ! == Works also if dummy atoms                                    ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: iat
    REAL(real_8)                             :: tau(3,maxsys%nax*maxsys%nsx), &
                                                x(3), d2_a(3,duat%ndat), &
                                                tran(3)

    INTEGER                                  :: id, ityp, naa

! ==--------------------------------------------------------------==

    IF (iat.LE.0) THEN
       ! Do nothing
       RETURN
    ELSEIF (iat.LE.ions1%nat) THEN
       ! Real atoms
       x(1)=tau(1,iat)
       x(2)=tau(2,iat)
       x(3)=tau(3,iat)
    ELSEIF (iat.LE.ions1%nat+duat%ndat) THEN
       ! Dummy atoms (type 2).
       naa=iat-ions1%nat
       ityp=duat%listda(naa,1)
       id=duat%listda(naa,2)
       IF (ityp.EQ.2) THEN
          x(1)=d2_a(1,id)
          x(2)=d2_a(2,id)
          x(3)=d2_a(3,id)
       ELSEIF (ityp.EQ.1) THEN
          x(1)=duat%dummy1(1,id) - tran(1)
          x(2)=duat%dummy1(2,id) - tran(2)
          x(3)=duat%dummy1(3,id) - tran(3)
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE fillc2
  ! ==================================================================
  ! ==================================================================

END MODULE meta_colvar_util_utils
