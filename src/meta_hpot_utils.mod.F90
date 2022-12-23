MODULE meta_hpot_utils
  USE cnst,                            ONLY: pi
  USE cnst_dyn,                        ONLY: &
       cnst_val, cscl_val, cv_path, hllh_val, hllw_val, ht_path, iangcv, &
       inter_hill, lmeta, rmeta
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: hills
  PUBLIC :: setwall_old
  PUBLIC :: setwall
  PUBLIC :: hills_lor
  PUBLIC :: hills_sals
  PUBLIC :: hills_ratio
  PUBLIC :: hills_sals_shift

CONTAINS

  ! ==================================================================
  SUBROUTINE  hills(cvar,nvar,i_meta,ntot,force,isys,ipos)
    ! ==--------------------------------------------------------------==


    INTEGER                                  :: nvar
    REAL(real_8)                             :: cvar(nvar)
    INTEGER                                  :: i_meta, ntot
    REAL(real_8)                             :: force(nvar)
    INTEGER                                  :: isys, ipos

    CHARACTER(*), PARAMETER                  :: procedureN = 'hills'

    INTEGER                                  :: ic, ierr, it, length
    REAL(real_8)                             :: expn, step_gauss
    REAL(real_8), ALLOCATABLE                :: diffpos(:)
    REAL(real_8), POINTER                    :: c_history(:,:)

    length  = nvar*4
    ALLOCATE(diffpos(nvar),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    IF (lmeta%lcnstdyn) THEN
       c_history => cnst_val(:,ipos+1:)
    ELSEIF (lmeta%lcolvardyn)THEN
       c_history => cv_path(:,ipos+1:)
    ELSEIF (lmeta%lmeta_cell) THEN
       c_history => ht_path(:,ipos+1:)
    ENDIF

    !$omp parallel do private(IC)
    DO ic=1,nvar
       force(ic)=0._real_8
       diffpos(ic)=0.0_real_8
    ENDDO

    rmeta%gausspot=0.0_real_8

    DO it=1,i_meta-1
       expn = 0.0_real_8
       DO ic=1,nvar
          diffpos(ic) = (cvar(ic)-c_history(it,ic))
          IF (iangcv(ic).EQ. 1) THEN
             IF (diffpos(ic) .GT. pi) THEN
                diffpos(ic) = diffpos(ic) - 2.0_real_8*pi
             ELSEIF (diffpos(ic) .LT. -pi) THEN
                diffpos(ic) = diffpos(ic) + 2.0_real_8*pi
             ENDIF
          ENDIF
          diffpos(ic) = diffpos(ic) /cscl_val(it,ic)
          expn = expn+diffpos(ic)**2.0_real_8
       ENDDO
       expn = SQRT(expn)
       step_gauss = hllh_val(it,isys)*EXP(-0.5_real_8*(expn/hllw_val(it,&
            isys))**2.0_real_8)
       rmeta%gausspot = rmeta%gausspot + step_gauss
       !$omp parallel do private(IC)
       DO ic=1,nvar
          force(ic) = force(ic) + (diffpos(ic)/(hllw_val(it,isys))**&
               2.0_real_8)*step_gauss
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    DEALLOCATE(diffpos,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE hills
  ! ==================================================================
  ! ==================================================================

  SUBROUTINE setwall_old(cpos,tycv,vparam,force)

    ! VPARAM    : parameters for the definition of the constraining potential
    ! (real, dimension 4)
    ! vbar = VPARAM(1) in a.u.
    ! r_wall or cmin = VPARAM(2)
    ! r_init or ddist or cmax =  VPARAM(3)
    ! deps  =  VPARAM(4)
    ! for DIST or angles  wall V=vbar*(cpos/r_wall)**4
    ! (for cpos>r_init V becomes active)
    ! for DIFFER  wall V=vbar*((|cpos|-r_wall)/ddist)**4
    ! (for (|cpos|-r_wall)<ddist V becomes active)
    ! for coordination numbers  and RMSD_AB
    ! wall min  : V=vbar(cmin+deps-cpos)**2
    ! (for (cpos-deps)<cmin V becomes active)
    ! wall max  : V=vbar(cpos-cmax+deps)**2
    ! (for (cpos+deps)>cmax V becomes active)

    REAL(real_8)                             :: cpos
    INTEGER                                  :: tycv
    REAL(real_8)                             :: vparam(4), force

    CHARACTER(*), PARAMETER                  :: procedureN = 'setwall_old'

    REAL(real_8)                             :: cc, cmax, cmin, ddist, delta, &
                                                deps, factor, num, r_init, &
                                                r_wall, vbar

    IF (tycv .EQ. 7) THEN
       vbar  = vparam(1)
       r_wall   = vparam(2)
       ddist    = vparam(3)
       delta     = ABS(cpos)
       IF (r_wall==0.0_real_8) &
            CALL stopgm(procedureN,'C1 term of CVSPACE BOUNDARY is equal to 0.0 and this leads to a division by zero!', &
            __LINE__,__FILE__)
       factor = vbar*4._real_8/(r_wall)**4.0_real_8
       force = 0.0_real_8
       IF (cpos .GT. ddist) THEN
          force = -factor*delta*delta*delta
       ELSEIF (cpos .LT. -ddist) THEN
          force =  factor*delta*delta*delta
       ENDIF
    ELSEIF((tycv.GE.8 .AND. tycv.LE.12) .OR. tycv.EQ.17 .OR. tycv.EQ.18.OR.tycv.EQ.29) THEN
       vbar= vparam(1)
       cmin    = vparam(2)
       cmax    = vparam(3)
       deps    = vparam(4)
       cc      = cpos
       force = 0.0_real_8
       !!cmb
       IF (cmin .LT. cmax .AND. cc-deps .LT. cmin) THEN
          num = cmin+deps-cc
          force = vbar*2._real_8*num
       ELSEIF (cmax .GT. cmin .AND. cc+deps .GT. cmax) THEN
          num = cc-cmax+deps
          force = -vbar*2._real_8*num
       ENDIF
       !!cmb
    ELSE
       vbar = vparam(1)
       r_wall   = vparam(2)
       r_init   = vparam(3)
       IF (r_wall==0.0_real_8) &
            CALL stopgm(procedureN,'C1 term of CVSPACE BOUNDARY is equal to 0.0 and this leads to a division by zero!', &
            __LINE__,__FILE__)
       factor = vbar*4._real_8/(r_wall)**4.0_real_8
       force = 0.0_real_8
       cc    = ABS(cpos)
       IF (cc .GT. r_init) THEN
          force = - factor*cc*cc*cc
       ENDIF
    ENDIF

    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE setwall_old
  ! ==================================================================
  ! ==================================================================

  SUBROUTINE setwall(cpos,vbound,force)
    REAL(real_8)                             :: cpos, vbound(4), force

    REAL(real_8)                             :: cc, cpos_0, vbarrier

    force = 0.0_real_8
    cpos_0   = vbound(1)
    vbarrier = vbound(2)
    IF (vbarrier.GT.0._real_8)THEN
       cc=cpos-cpos_0
       IF (cc.GT.0._real_8)THEN
          force = - vbarrier
       ENDIF
    ENDIF
    cpos_0   = vbound(3)
    vbarrier = vbound(4)
    IF (vbarrier.GT.0._real_8)THEN
       cc=cpos-cpos_0
       IF (cc.LT.0._real_8)THEN
          force =   vbarrier
       ENDIF
    ENDIF

    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE setwall
  ! ==================================================================


  SUBROUTINE  hills_lor(cvar,nvar,i_meta,ntot,force,isys,ipos)



    INTEGER                                  :: nvar
    REAL(real_8)                             :: cvar(nvar)
    INTEGER                                  :: i_meta, ntot
    REAL(real_8)                             :: force(nvar)
    INTEGER                                  :: isys, ipos

    CHARACTER(*), PARAMETER                  :: procedureN = 'hills_lor'

    INTEGER                                  :: ic, ierr, it, length
    REAL(real_8)                             :: denom, expn, hw, step_gauss
    REAL(real_8), ALLOCATABLE                :: diffpos(:)
    REAL(real_8), POINTER                    :: c_history(:,:)

    length  = nvar*4
    ALLOCATE(diffpos(nvar),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    IF (lmeta%lcnstdyn) THEN
       c_history => cnst_val(:,ipos+1:)
    ELSEIF (lmeta%lcolvardyn)THEN
       c_history => cv_path(:,ipos+1:)
    ELSEIF (lmeta%lmeta_cell) THEN
       c_history => ht_path(:,ipos+1:)
    ENDIF

    !$omp parallel do private(IC)
    DO ic=1,nvar
       force(ic)=0._real_8
       diffpos(ic)=0.0_real_8
    ENDDO

    DO it=1,i_meta-1
       expn = 0.0_real_8
       hw   =   hllw_val(it,isys) * SQRT(3.0_real_8)
       DO ic=1,nvar
          diffpos(ic) = (cvar(ic)-c_history(it,ic))
          IF (iangcv(ic).EQ. 1) THEN
             IF (diffpos(ic) .GT. pi) THEN
                diffpos(ic) = diffpos(ic) - 2.0_real_8*pi
             ELSEIF (diffpos(ic) .LT. -pi) THEN
                diffpos(ic) = diffpos(ic) + 2.0_real_8*pi
             ENDIF
          ENDIF
          diffpos(ic) = diffpos(ic) /cscl_val(it,ic)
          expn = expn+diffpos(ic)**2.0_real_8
       ENDDO
       expn       = SQRT(expn)
       denom      = (expn/hw)**2.0_real_8+1.0_real_8
       step_gauss = hllh_val(it,isys)/denom

       !$omp parallel do private(IC)
       DO ic=1,nvar
          force(ic) = force(ic) + 2.0_real_8*(diffpos(ic)/(hw)**2.0_real_8)*&
               step_gauss/denom
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    DEALLOCATE(diffpos,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE hills_lor
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE  hills_sals(cvnow,nvar,i_meta,ntot,force,cv_last,&
       i_cvst,isys,ipos)



    INTEGER                                  :: nvar
    REAL(real_8)                             :: cvnow(nvar)
    INTEGER                                  :: i_meta, ntot
    REAL(real_8)                             :: force(nvar), cv_last(nvar)
    INTEGER                                  :: i_cvst, isys, ipos

    CHARACTER(*), PARAMETER                  :: procedureN = 'hills_sals'

    INTEGER                                  :: ic, ierr, it, length
    REAL(real_8)                             :: expn, expn2, fact, hw2, &
                                                sigma2, step_gauss
    REAL(real_8), ALLOCATABLE                :: deltacv(:), diffpos(:)
    REAL(real_8), POINTER                    :: c_history(:,:)

    fact = 1.0_real_8
    length  = nvar*5

    ALLOCATE(diffpos(nvar),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(deltacv(nvar),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    IF (lmeta%lcnstdyn) THEN
       c_history => cnst_val(:,ipos+1:)
    ELSEIF (lmeta%lcolvardyn)THEN
       c_history => cv_path(:,ipos+1:)
    ELSEIF (lmeta%lmeta_cell) THEN
       c_history => ht_path(:,ipos+1:)
    ENDIF

    !$omp parallel do private(IC)
    DO ic=1,nvar
       force(ic)=0._real_8
       diffpos(ic)=0.0_real_8
       deltacv(ic)=0.0_real_8
    ENDDO

    rmeta%gausspot=0.0_real_8

    DO it=1,i_meta-2
       expn   = 0.0_real_8
       expn2  = 0.0_real_8
       sigma2 = 0.0_real_8
       DO ic=1,nvar
          diffpos(ic) = (cvnow(ic)-c_history(it,ic))
          deltacv(ic) = (c_history(it+1,ic)-c_history(it,ic))
          IF (iangcv(ic).EQ. 1) THEN
             IF (diffpos(ic) .GT. pi) THEN
                diffpos(ic) = diffpos(ic) - 2.0_real_8*pi
             ELSEIF (diffpos(ic) .LT. -pi) THEN
                diffpos(ic) = diffpos(ic) + 2.0_real_8*pi
             ENDIF
             IF (deltacv(ic) .GT. pi) THEN
                deltacv(ic) = deltacv(ic) - 2.0_real_8*pi
             ELSEIF (deltacv(ic) .LT. -pi) THEN
                deltacv(ic) = deltacv(ic) + 2.0_real_8*pi
             ENDIF
          ENDIF
          diffpos(ic) = diffpos(ic)/cscl_val(it,ic)
          deltacv(ic) = deltacv(ic)/cscl_val(it,ic)

          expn   = expn   + diffpos(ic)**2.0_real_8
          expn2  = expn2  + diffpos(ic)*deltacv(ic)
          sigma2 = sigma2 + deltacv(ic)*deltacv(ic)
       ENDDO
       ! EXPN = SQRT(EXPN)
       hw2 = hllw_val(it,isys)**2.0_real_8
       IF (sigma2 .GT. hw2) sigma2=hw2
       sigma2 = sigma2*sigma2
       IF (sigma2 .LT. 1.e-6_real_8)  sigma2 = 1.e-6_real_8
       step_gauss = hllh_val(it,isys)*EXP(-0.5_real_8*(expn/hw2))*EXP(-&
            0.5_real_8*(expn2*expn2/sigma2))
       !$omp parallel do private(IC)
       DO ic=1,nvar
          force(ic) = force(ic) + (diffpos(ic)/(hllw_val(it,isys))**&
               2.0_real_8)*step_gauss+(expn2*deltacv(ic)/sigma2)*step_gauss
       ENDDO
       rmeta%gausspot = rmeta%gausspot+step_gauss
    ENDDO

    IF (i_meta .GT. 1) THEN
       expn   = 0.0_real_8
       expn2  = 0.0_real_8
       sigma2 = 0.0_real_8
       DO ic=1,nvar
          diffpos(ic) = (cvnow(ic)-c_history(i_meta-1,ic))
          deltacv(ic) = (cv_last(ic)-c_history(i_meta-1,ic))
          IF (iangcv(ic).EQ. 1) THEN
             IF (diffpos(ic) .GT. pi) THEN
                diffpos(ic) = diffpos(ic) - 2.0_real_8*pi
             ELSEIF (diffpos(ic) .LT. -pi) THEN
                diffpos(ic) = diffpos(ic) + 2.0_real_8*pi
             ENDIF
             IF (deltacv(ic) .GT. pi) THEN
                deltacv(ic) = deltacv(ic) - 2.0_real_8*pi
             ELSEIF (deltacv(ic) .LT. -pi) THEN
                deltacv(ic) = deltacv(ic) + 2.0_real_8*pi
             ENDIF
          ENDIF
          diffpos(ic) = diffpos(ic)/cscl_val(i_meta-1,ic)
          deltacv(ic) = deltacv(ic)/cscl_val(i_meta-1,ic)

          expn   = expn   + diffpos(ic)**2.0_real_8
          expn2  = expn2  + diffpos(ic)*deltacv(ic)
          sigma2 = sigma2 + deltacv(ic)*deltacv(ic)
       ENDDO
       ! EXPN = SQRT(EXPN)
       hw2 = hllw_val(i_meta-1,isys)**2.0_real_8
       IF (sigma2 .GT. hw2) sigma2=hw2
       IF (sigma2 .LT. 1.e-6_real_8)  sigma2 = 1.e-6_real_8
       sigma2 = sigma2*sigma2
       step_gauss =fact*hllh_val(i_meta-1,isys)*EXP(-0.5_real_8*(expn/hw2))&
            *EXP(-0.5_real_8*(expn2*expn2/sigma2))
       !$omp parallel do private(IC)
       DO ic=1,nvar
          force(ic) = force(ic) + (diffpos(ic)/(hllw_val(i_meta-1,isys))&
               **2.0_real_8)*step_gauss+(expn2*deltacv(ic)/sigma2)*step_gauss
       ENDDO
       rmeta%gausspot = rmeta%gausspot+step_gauss
    ENDIF

    ! ==--------------------------------------------------------------==
    DEALLOCATE(diffpos,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(deltacv,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE hills_sals

  ! ==================================================================
  SUBROUTINE  hills_ratio(cvnow,nvar,i_meta,ntot,force,cv_last,&
       i_cvst,isys,ipos)



    INTEGER                                  :: nvar
    REAL(real_8)                             :: cvnow(nvar)
    INTEGER                                  :: i_meta, ntot
    REAL(real_8)                             :: force(nvar), cv_last(nvar)
    INTEGER                                  :: i_cvst, isys, ipos

    CHARACTER(*), PARAMETER                  :: procedureN = 'hills_ratio'

    INTEGER                                  :: ic, ierr, it
    REAL(real_8)                             :: arg1, arg2, den, expn, expn2, &
                                                fact, hw2, larger, num, pow1, &
                                                pow2, scale_exp, sigma2, &
                                                step_gauss
    REAL(real_8), ALLOCATABLE                :: deltacv(:), diffpos(:)
    REAL(real_8), POINTER                    :: c_history(:,:)

    fact = 1.0_real_8
    larger = rmeta%fboost*rmeta%fboost
    pow1 = rmeta%expup
    pow2 = rmeta%expdown


    ALLOCATE(diffpos(nvar),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(deltacv(nvar),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    IF (lmeta%lcnstdyn) THEN
       c_history => cnst_val(:,ipos+1:)
    ELSEIF (lmeta%lcolvardyn )THEN
       c_history => cv_path(:,ipos+1:)
    ELSEIF (lmeta%lmeta_cell) THEN
       c_history => ht_path(:,ipos+1:)
    ENDIF

    !$omp parallel do private(IC)
    DO ic=1,nvar
       force(ic)=0._real_8
       diffpos(ic)=0.0_real_8
       deltacv(ic)=0.0_real_8
    ENDDO

    rmeta%gausspot=0.0_real_8

    DO it=1,i_meta-2
       expn   = 0.0_real_8
       expn2  = 0.0_real_8
       sigma2 = 0.0_real_8
       DO ic=1,nvar
          diffpos(ic) = (cvnow(ic)-c_history(it,ic))
          deltacv(ic) = (c_history(it+1,ic)-c_history(it,ic))
          IF (iangcv(ic).EQ. 1) THEN
             IF (diffpos(ic) .GT. pi) THEN
                diffpos(ic) = diffpos(ic) - 2.0_real_8*pi
             ELSEIF (diffpos(ic) .LT. -pi) THEN
                diffpos(ic) = diffpos(ic) + 2.0_real_8*pi
             ENDIF
             IF (deltacv(ic) .GT. pi) THEN
                deltacv(ic) = deltacv(ic) - 2.0_real_8*pi
             ELSEIF (deltacv(ic) .LT. -pi) THEN
                deltacv(ic) = deltacv(ic) + 2.0_real_8*pi
             ENDIF
          ENDIF
          diffpos(ic) = diffpos(ic)/cscl_val(it,ic)
          deltacv(ic) = deltacv(ic)/cscl_val(it,ic)

          expn   = expn   + diffpos(ic)**2.0_real_8
          expn2  = expn2  + diffpos(ic)*deltacv(ic)
          sigma2 = sigma2 + deltacv(ic)*deltacv(ic)
       ENDDO
       expn = SQRT(expn)
       hw2=hllw_val(it,isys)**2.0_real_8
       IF (sigma2 .GT. hw2) sigma2=hw2
       IF (sigma2 .LT. 1.e-6_real_8)  sigma2 = 1.e-6_real_8
       sigma2 = 1.0_real_8/(sigma2*sigma2*larger)
       arg1 = (expn/hllw_val(it,isys))**pow1
       arg2 = (expn/hllw_val(it,isys))**pow2
       num  = 1.0_real_8 - arg1
       den  = 1.0_real_8 - arg2
       IF (den .EQ. 0.0_real_8) GOTO 11
       den  = 1.0_real_8/den
       scale_exp = EXP(-0.5_real_8*(expn2*expn2*sigma2))
       step_gauss = hllh_val(it,isys)*num*scale_exp*den
       !$omp parallel do private(IC)
       DO ic=1,nvar
          force(ic) = force(ic) + hllh_val(it,isys) * scale_exp*den*((&
               pow1*arg1 - num*pow2*arg2*den)*diffpos(ic)/(expn**2.0_real_8)+num*&
               expn2*deltacv(ic)*sigma2)
       ENDDO
       rmeta%gausspot = rmeta%gausspot+step_gauss
11     CONTINUE
    ENDDO
    IF (i_cvst .LE. inter_hill ) fact = REAL(i_cvst+1,kind=real_8)/REAL(inter_hill&
         ,kind=real_8)
    IF (i_meta .GT. 1) THEN
       expn   = 0.0_real_8
       expn2  = 0.0_real_8
       sigma2 = 0.0_real_8
       DO ic=1,nvar
          diffpos(ic) = (cvnow(ic)-c_history(i_meta-1,ic))
          deltacv(ic) = (cv_last(ic)-c_history(i_meta-1,ic))
          IF (iangcv(ic).EQ. 1) THEN
             IF (diffpos(ic) .GT. pi) THEN
                diffpos(ic) = diffpos(ic) - 2.0_real_8*pi
             ELSEIF (diffpos(ic) .LT. -pi) THEN
                diffpos(ic) = diffpos(ic) + 2.0_real_8*pi
             ENDIF
             IF (deltacv(ic) .GT. pi) THEN
                deltacv(ic) = deltacv(ic) - 2.0_real_8*pi
             ELSEIF (deltacv(ic) .LT. -pi) THEN
                deltacv(ic) = deltacv(ic) + 2.0_real_8*pi
             ENDIF
          ENDIF
          diffpos(ic) = diffpos(ic)/cscl_val(i_meta-1,ic)
          deltacv(ic) = deltacv(ic)/cscl_val(i_meta-1,ic)

          expn   = expn   + diffpos(ic)**2.0_real_8
          expn2  = expn2  + diffpos(ic)*deltacv(ic)
          sigma2 = sigma2 + deltacv(ic)*deltacv(ic)
       ENDDO
       expn = SQRT(expn)
       hw2=hllw_val(i_meta-1,isys)**2.0_real_8
       IF (sigma2 .GT. hw2) sigma2=hw2
       IF (sigma2 .LT. 1.e-6_real_8)  sigma2 = 1.e-6_real_8
       sigma2 = 1.0_real_8/(sigma2*sigma2*larger)
       arg1 = (expn/hllw_val(i_meta-1,isys))**pow1
       arg2 = (expn/hllw_val(i_meta-1,isys))**pow2
       num  = 1.0_real_8 - arg1
       den  = 1.0_real_8 - arg2
       IF (den .EQ. 0.0_real_8) GOTO 12
       den  = 1.0_real_8/den
       scale_exp = EXP(-0.5_real_8*(expn2*expn2*sigma2))
       step_gauss = fact*hllh_val(i_meta-1,isys)*num*scale_exp*den
       !$omp parallel do private(IC)
       DO ic=1,nvar
          force(ic) = force(ic) + fact*hllh_val(i_meta-1,isys) *&
               scale_exp * den * ((pow1*arg1- num*pow2*arg2*den)*diffpos(ic)&
               /(expn**2.0_real_8)+ num*expn2*deltacv(ic)*sigma2)
       ENDDO
       rmeta%gausspot = rmeta%gausspot+step_gauss
12     CONTINUE
    ENDIF
    ! ==--------------------------------------------------------------==
    DEALLOCATE(diffpos,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(deltacv,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE hills_ratio

  ! ==================================================================
  SUBROUTINE  hills_sals_shift(cvnow,nvar,i_meta,ntot,force,cv_last,&
       i_cvst,isys,ipos)



    INTEGER                                  :: nvar
    REAL(real_8)                             :: cvnow(nvar)
    INTEGER                                  :: i_meta, ntot
    REAL(real_8)                             :: force(nvar), cv_last(nvar)
    INTEGER                                  :: i_cvst, isys, ipos

    CHARACTER(*), PARAMETER :: procedureN = 'hills_sals_shift'

    INTEGER                                  :: ic, ierr, it, length
    REAL(real_8)                             :: expn, expn2, fact, fsig, hw2, &
                                                rcut, sigma2, step1, step2, &
                                                step3, step_gauss
    REAL(real_8), ALLOCATABLE                :: deltacv(:), diffpos(:)
    REAL(real_8), POINTER                    :: c_history(:,:)

    fact = 1.0_real_8
    rcut = rmeta%rshift
    fsig = rmeta%fboost
    length  = nvar*5

    ALLOCATE(diffpos(nvar),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(deltacv(nvar),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    IF (lmeta%lcnstdyn) THEN
       c_history => cnst_val(:,ipos+1:)
    ELSEIF (lmeta%lcolvardyn)THEN
       c_history => cv_path(:,ipos+1:)
    ELSEIF (lmeta%lmeta_cell) THEN
       c_history => ht_path(:,ipos+1:)
    ENDIF

    !$omp parallel do private(IC)
    DO ic=1,nvar
       force(ic)=0._real_8
       diffpos(ic)=0.0_real_8
       deltacv(ic)=0.0_real_8
    ENDDO

    rmeta%gausspot=0.0_real_8

    DO it=1,i_meta-2
       expn   = 0.0_real_8
       expn2  = 0.0_real_8
       sigma2 = 0.0_real_8
       DO ic=1,nvar
          diffpos(ic) = (cvnow(ic)-c_history(it,ic))
          deltacv(ic) = (c_history(it+1,ic)-c_history(it,ic))
          IF (iangcv(ic).EQ. 1) THEN
             IF (diffpos(ic) .GT. pi) THEN
                diffpos(ic) = diffpos(ic) - 2.0_real_8*pi
             ELSEIF (diffpos(ic) .LT. -pi) THEN
                diffpos(ic) = diffpos(ic) + 2.0_real_8*pi
             ENDIF
             IF (deltacv(ic) .GT. pi) THEN
                deltacv(ic) = deltacv(ic) - 2.0_real_8*pi
             ELSEIF (deltacv(ic) .LT. -pi) THEN
                deltacv(ic) = deltacv(ic) + 2.0_real_8*pi
             ENDIF
          ENDIF
          diffpos(ic) = diffpos(ic)/cscl_val(it,ic)
          deltacv(ic) = deltacv(ic)/cscl_val(it,ic)

          expn   = expn   + diffpos(ic)**2.0_real_8
          expn2  = expn2  + diffpos(ic)*deltacv(ic)
          sigma2 = sigma2 + deltacv(ic)*deltacv(ic)
       ENDDO
       expn = SQRT(expn)
       hw2 = hllw_val(it,isys)**2.0_real_8
       IF (sigma2 .GT. hw2) sigma2 = hw2
       sigma2 = sigma2*sigma2*fsig*fsig
       IF (sigma2 .LT. 1.e-6_real_8)  sigma2 = 1.e-6_real_8

       IF (expn .LE. rcut*hllw_val(it,isys)) THEN
          step1 = EXP(-0.5_real_8*(expn*expn/hw2))
          step2 = EXP(-0.5_real_8*(rcut)**2.0_real_8)
          step3 = EXP(-0.5_real_8*(expn2*expn2/sigma2))
          step_gauss = hllh_val(it,isys)*((step1 - step2) * step3)
       ELSE
          step1 = 0.0_real_8
          step2 = 0.0_real_8
          step3 = 0.0_real_8
          step_gauss = 0.0_real_8
       ENDIF

       !$omp parallel do private(IC)
       DO ic=1,nvar

          force(ic) = force(ic) +  hllh_val(it,isys)*(step1*step3*&
               diffpos(ic)/(hllw_val(it,isys))**2.0_real_8+(step1-step2)*step3*&
               expn2*deltacv(ic)/sigma2)
       ENDDO
       rmeta%gausspot = rmeta%gausspot+step_gauss
    ENDDO

    IF (i_cvst .LE. inter_hill ) fact = REAL(i_cvst+1,kind=real_8)/REAL(inter_hill&
         ,kind=real_8)

    IF (i_meta .GT. 1) THEN
       expn   = 0.0_real_8
       expn2  = 0.0_real_8
       sigma2 = 0.0_real_8
       DO ic=1,nvar
          diffpos(ic) = (cvnow(ic)-c_history(i_meta-1,ic))
          deltacv(ic) = (cv_last(ic)-c_history(i_meta-1,ic))
          IF (iangcv(ic).EQ. 1) THEN
             IF (diffpos(ic) .GT. pi) THEN
                diffpos(ic) = diffpos(ic) - 2.0_real_8*pi
             ELSEIF (diffpos(ic) .LT. -pi) THEN
                diffpos(ic) = diffpos(ic) + 2.0_real_8*pi
             ENDIF
             IF (deltacv(ic) .GT. pi) THEN
                deltacv(ic) = deltacv(ic) - 2.0_real_8*pi
             ELSEIF (deltacv(ic) .LT. -pi) THEN
                deltacv(ic) = deltacv(ic) + 2.0_real_8*pi
             ENDIF
          ENDIF
          diffpos(ic) = diffpos(ic)/cscl_val(i_meta-1,ic)
          deltacv(ic) = deltacv(ic)/cscl_val(i_meta-1,ic)

          expn   = expn   + diffpos(ic)**2.0_real_8
          expn2  = expn2  + diffpos(ic)*deltacv(ic)
          sigma2 = sigma2 + deltacv(ic)*deltacv(ic)
       ENDDO
       expn = SQRT(expn)
       hw2 = hllw_val(i_meta-1,isys)**2.0_real_8
       IF (sigma2 .GT. hw2) sigma2 = hw2
       IF (sigma2 .LT. 1.e-6_real_8)  sigma2 = 1.e-6_real_8
       sigma2 = sigma2*sigma2*fsig*fsig

       IF (expn .LE. rcut*hllw_val(i_meta-1,isys)) THEN
          step1 = EXP(-0.5_real_8*(expn*expn/hw2))
          step2 = EXP(-0.5_real_8*(rcut)**2.0_real_8)
          step3 = EXP(-0.5_real_8*(expn2*expn2/sigma2))
          step_gauss =fact*hllh_val(i_meta-1,isys)*((step1 - step2) *&
               step3)
       ELSE
          step1 = 0.0_real_8
          step2 = 0.0_real_8
          step3 = 0.0_real_8
          step_gauss = 0.0_real_8
       ENDIF
       !$omp parallel do private(IC)
       DO ic=1,nvar
          force(ic) = force(ic) +  fact*hllh_val(i_meta-1,isys)*(step1*&
               step3*diffpos(ic)/(hllw_val(i_meta-1,isys))**2.0_real_8+(step1-&
               step2)*step3*expn2*deltacv(ic)/sigma2)
       ENDDO
       rmeta%gausspot = rmeta%gausspot+step_gauss
    ENDIF

    ! ==--------------------------------------------------------------==
    DEALLOCATE(diffpos,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(deltacv,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE hills_sals_shift
  ! ==================================================================

END MODULE meta_hpot_utils
