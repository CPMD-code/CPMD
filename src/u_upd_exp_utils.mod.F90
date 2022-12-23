MODULE u_upd_exp_utils
  USE error_handling,                  ONLY: stopgm
  USE g_loc,                           ONLY: gloc_re,&
                                             glocal,&
                                             gloci,&
                                             glocr,&
                                             nstep,&
                                             omega_n
  USE g_loc_util_utils,                ONLY: lagrmult,&
                                             line_search_exe,&
                                             u_by_ortho,&
                                             z_update
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: paral
  USE system,                          ONLY: cntl,&
                                             fpar
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: zgive
  USE zeroing_utils,                   ONLY: zeroing
  USE znum_mat_utils,                  ONLY: give_scr_znum_mat

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: u_upd_exp_ide
  PUBLIC :: u_grad_exp_ide
  PUBLIC :: give_scr_upd_exp_ide

CONTAINS

  ! ==================================================================
  SUBROUTINE u_upd_exp_ide(u_mat,z_mat,xyzmat,u_grad,nstate)
    ! ==--------------------------------------------------------------==
    ! == CALCULATES THE UNITARY TRANSFORMATION 
    ! ==            WHICH OPTIMISE LOCALIZATION IN G SPACE            ==
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8)                          :: xyzmat(*)
    INTEGER                                  :: nstate
    COMPLEX(real_8), TARGET                  :: u_grad(nstate,nstate)
    COMPLEX(real_8)                          :: z_mat(nstate,nstate,3), &
                                                u_mat(nstate,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'u_upd_exp_ide'

    COMPLEX(real_8)                          :: carg, zone, zzero
    COMPLEX(real_8), ALLOCATABLE             :: a(:,:), b1(:,:), b2(:,:), &
                                                gam(:,:)
    COMPLEX(real_8), POINTER                 :: u_temp(:,:)
    INTEGER                                  :: i, ierr, igrmax, isub, &
                                                izamax, j
    INTEGER, SAVE                            :: icont = 0
    REAL(real_8)                             :: delta, difgam, grmax
    REAL(real_8), SAVE                       :: ofun00 = 0.0_real_8, &
                                                step_fac = 1.0_real_8

    CALL tiset(' UPD_EXP_I',isub)

    IF (cntl%tlsd)  CALL stopgm('TRANUPD','LSD NOT IMPLEMENTED YET',& 
         __LINE__,__FILE__)

    ALLOCATE(a(nstate, nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(b1(nstate, nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(b2(nstate, nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(gam(nstate, nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    u_temp => u_grad

    zone  = CMPLX(1.0_real_8,0.0_real_8,kind=real_8)
    zzero = CMPLX(0.0_real_8,0.0_real_8,kind=real_8)

    ! time step
    delta = glocr%gloc_step * step_fac

    ! this is U_MAT

    ! Up_MAT is zero if first step
    CALL zeroing(a)!, SIZE(a))
    CALL zeroing(b1)!, SIZE(b1))
    CALL zeroing(b2)!, SIZE(b2))
    CALL zeroing(gam)!, SIZE(gam))

    nstep = nstep + 1

    igrmax=izamax(nstate*nstate,u_grad,1)
    grmax = ABS(zgive(u_grad,igrmax))

    IF (gloci%gloc_const .EQ. 3) THEN

       ! UPDATE U_MAT WITH APPROXIMATION U = 1 + A (A antihermitian)
       CALL approx_u(a,u_grad,nstate)
       CALL dcopy(2*nstate*nstate,a,1,u_temp,1)

    ELSEIF (gloci%gloc_const .EQ. 2) THEN

       IF (glocal%tg_linesearch) THEN
          CALL line_search_exe(a,u_grad,z_mat,xyzmat,b1,b2,&
               nstate,delta,step_fac)
       ELSE
          CALL u_by_ortho(a,u_grad,b1,b2,nstate,delta)
       ENDIF
       CALL dcopy(2*nstate*nstate,a,1,u_temp,1)

    ELSEIF (gloci%gloc_const .EQ. 1) THEN

       ! update U_MAT_new =  U_MAT_old + dt U_GRAD

       DO i=1,nstate
          DO j=1,nstate
             u_temp(j,i)=u_mat(j,i)+delta*u_grad(j,i)
          ENDDO
       ENDDO

       CALL lagrmult(u_mat,u_temp,a,b1,b2,gam,nstate)

    ENDIF  ! APPROX

    gloc_re%gmax = 0.0_real_8
    gloc_re%ggnorm = 0.0_real_8
    DO i = 1,nstate
       DO j = 1,nstate
          carg = u_temp(j,i)-u_mat(j,i)
          difgam=REAL(carg)*REAL(carg)+AIMAG(carg)*AIMAG(carg)
          gloc_re%ggnorm = gloc_re%ggnorm + difgam
          gloc_re%gmax = MAX(gloc_re%gmax,difgam)
       ENDDO
    ENDDO
    gloc_re%gmax = SQRT(gloc_re%gmax)
    gloc_re%ggnorm = SQRT(gloc_re%ggnorm)/REAL(nstate,kind=real_8)


    CALL dcopy(2*nstate*nstate,u_temp(1,1),1,u_mat(1,1),1)

    ! calculate U_GRAD_new with the rotated U_MAT_new
    CALL z_update(u_mat,z_mat,nstate,nstate,3)
    CALL gzfunc(z_mat,nstate,nstate)

    ! IF(.NOT. TG_LINESEARCH) THEN
    ! DIFF = OFUN - OFUN0
    ! write(6,*) DIFF,OFUN,OFUN0
    ! stop
    ! IF(abs(DIFF) .GT. GNORM) THEN 
    ! STEP_FAC = STEP_FAC * (1.0_real_8 - GNORM/DIFF)
    ! ELSE
    ! write(6,'(/,A,I10,A,1PE12.6,/)')
    ! &     '   U_UPD: WARNING ... at step ',NSTEP,
    ! &     '  DIFF_OFUN = ', DIFF
    ! ENDIF
    ! ENDIF

    gloc_re%dif_fun = ABS(gloc_re%ofun - ofun00)
    ofun00 = gloc_re%ofun
    IF (glocal%tg_antisymm .AND. icont .LE. 100  )  THEN
       icont = icont+1
       IF (paral%io_parent)&
            WRITE(10,'(I10,3(2x,1PE12.6),1PE12.6,2x,1PE12.6,1x,f10.4)')&
            nstep,gloc_re%gmax,gloc_re%ofun,grmax,gloc_re%dif_fun,gloc_re%omega_tot,step_fac
       IF (paral%io_parent)&
            WRITE(12,'(I10,16f10.4)') nstep,(omega_n(i),i=1,16)
    ELSE
       IF (paral%io_parent)&
            WRITE(10,'(I10,3(2x,1PE12.6),1PE12.6,2x,1PE12.6,f10.5)')&
            nstep,gloc_re%gmax,gloc_re%ofun,grmax,gloc_re%dif_fun,gloc_re%ggnorm,step_fac/10.0_real_8

    ENDIF

    ! ==--------------------------------------------------------------==
    DEALLOCATE(a,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(b1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(b2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(gam,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==

    CALL tihalt(' UPD_EXP_I',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE u_upd_exp_ide
  ! ==================================================================
  SUBROUTINE  u_grad_exp_ide(u_mat,u_grad,z_mat,nstate)
    ! ==--------------------------------------------------------------==
    ! ==  CALCULATES THE GRADIENT OF THE LOCALIZATION FUNCTIONAL      ==
    ! ==  WITH RESPECT TO THE UNITARY MATRIX COEFFICIENTS             ==
    ! ==  U_GRAD_ij = -d OFUNC / d U_ij                               ==
    ! ==--------------------------------------------------------------==
    ! input
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: z_mat(nstate,nstate,3), &
                                                u_grad(nstate,nstate), &
                                                u_mat(nstate,nstate)

    INTEGER                                  :: i, j, k
    REAL(real_8)                             :: fac

! ==--------------------------------------------------------------==

    fac = 1.0_real_8/6.0_real_8
    DO k = 1,3
       DO j = 1,nstate
          DO i = 1,nstate
             u_grad(i,j) = u_grad(i,j)&
                  +z_mat(i,j,k)*CONJG(z_mat(j,j,k))&
                  +z_mat(j,j,k)*CONJG(z_mat(j,i,k))
          ENDDO
       ENDDO
    ENDDO

    DO i = 1,nstate
       u_grad(i,i) = u_grad(i,i) - 6.0_real_8
    ENDDO
    ! FAC = 10._real_8
    ! DO I = 1,NSTATE
    ! DO J = 1,NSTATE
    ! IF(J .NE. I) THEN
    ! U_GRAD(J,I) = U_GRAD(J,I)*FAC
    ! ENDIF
    ! ENDDO
    ! ENDDO

    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE u_grad_exp_ide
  ! ==--------------------------------------------------------------==
  ! ==--------------------------------------------------------------==
  SUBROUTINE give_scr_upd_exp_ide(ltranupd,tag,nstate)
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: ltranupd
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate

    INTEGER                                  :: l_multi, l_ortho, l_print, &
                                                l_search, lgrad, lznum

    CALL give_scr_znum_mat(lznum,tag,nstate)
    ! ==--------------------------------------------------------------==
    lgrad     = 2*nstate*nstate
    l_ortho   = lgrad+3*2*nstate*nstate + 11*nstate
    l_multi   = lgrad+6*2*nstate*nstate
    IF (glocal%tg_linesearch) THEN
       l_search   = lgrad + 9*2*nstate*nstate
    ELSE
       l_search   = 1
    ENDIF
    IF (glocal%tglocrealp) THEN
       l_print = fpar%kr1*fpar%kr2s*fpar%kr3s + fpar%kr2s*fpar%kr3s
    ELSE
       l_print =1
    ENDIF
    ! ==--------------------------------------------------------------==
    ltranupd  = MAX(l_ortho,l_multi,lznum,l_search,l_print)
    tag='MAX(L_ORTHO,L_MULTI,LZNUM,...)'
    RETURN
  END SUBROUTINE give_scr_upd_exp_ide
  ! ==--------------------------------------------------------------==
  ! ==================================================================
  ! ==================================================================

END MODULE u_upd_exp_utils
