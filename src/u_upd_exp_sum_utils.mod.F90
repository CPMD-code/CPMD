MODULE u_upd_exp_sum_utils
  USE error_handling,                  ONLY: stopgm
  USE g_loc,                           ONLY: gifirst,&
                                             gloc_re,&
                                             glocal,&
                                             gloci,&
                                             glocr,&
                                             nstep,&
                                             omega_n
  USE g_loc_optim_utils,               ONLY: calc_u_wgrad_sum
  USE g_loc_util_utils,                ONLY: lagrmult,&
                                             line_search_exe,&
                                             u_by_ortho,&
                                             z_update
  USE kinds,                           ONLY: int_1,&
                                             int_2,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8
  USE parac,                           ONLY: parai,&
                                             paral
  USE system,                          ONLY: cntl
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: zgive
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: u_upd_exp_sum
  PUBLIC :: calc_u_grad_sum

  !!public :: simm_op_sum

CONTAINS

  ! ==================================================================
  SUBROUTINE u_upd_exp_sum(u_mat,z_mat,xyzmat,u_grad,nstate)
    ! ==--------------------------------------------------------------==
    ! == CALCULATES THE UNITARY TRANSFORMATION 
    ! ==            WHICH OPTIMISE LOCALIZATION IN G SPACE            ==
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8)                          :: xyzmat(:,:,:)
    INTEGER                                  :: nstate
    COMPLEX(real_8), TARGET                  :: u_grad(nstate,nstate)
    COMPLEX(real_8)                          :: z_mat(nstate,nstate,3), &
                                                u_mat(nstate,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'u_upd_exp_sum'

    COMPLEX(real_8)                          :: carg, zone, zzero
    COMPLEX(real_8), ALLOCATABLE             :: a(:,:), b1(:,:), b2(:,:), &
                                                gam(:,:)
    COMPLEX(real_8), POINTER                 :: u_temp(:,:)
    INTEGER                                  :: i, ierr, igrmax, igrmax_r, &
                                                isub, izamax, j
    INTEGER, SAVE                            :: icont = 0
    REAL(real_8)                             :: delta, dif_gfun, dif_xfun, &
                                                difgam, grmax, grmax_r
    REAL(real_8), SAVE                       :: ofun00 = 0.0_real_8, &
                                                step_fac = 1.0_real_8, &
                                                xyzfun00 = 0.0_real_8

    CALL tiset(' UPD_EXP_S',isub)

    IF (cntl%tlsd) CALL stopgm('TRANUPD','LSD NOT IMPLEMENTED YET',& 
         __LINE__,__FILE__)

    IF (gifirst .EQ. 1) THEN
       icont    = 0
       step_fac = 1.0_real_8
       ofun00   = gloc_re%ofun
       xyzfun00 = gloc_re%xyzfun
       gifirst  = 0
    ENDIF

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
    delta = glocr%gloc_step*step_fac

    ! this is U_MAT

    ! Up_MAT is zero if first step
    CALL zeroing(a)!, SIZE(a))
    CALL zeroing(b1)!, SIZE(b1))
    CALL zeroing(b2)!, SIZE(b2))
    CALL zeroing(gam)!, SIZE(gam))

    nstep = nstep + 1

    IF (gloci%gloc_const .EQ. 2) THEN

       CALL dcopy(2*nstate*nstate,u_mat,1,a,1)

       IF (glocal%tg_linesearch) THEN
          CALL line_search_exe(a,u_grad,z_mat,xyzmat,b1,b2,&
               nstate,delta,step_fac)
       ELSE
          CALL u_by_ortho(a,u_grad,b1,b2,nstate,delta)
       ENDIF
       CALL dcopy(2*nstate*nstate,a,1,u_temp,1)

    ELSEIF (gloci%gloc_const .EQ. 1) THEN
       ! IF(TCGRAD) THEN 
       ! DO I = 1,NSTATE
       ! A(I,I) = CMPLX(1.0_real_8,0.0_real_8)
       ! ENDDO

       ! CALL GZ_CONJ_GRAD(A,U_GRAD,DELTA,NSTATE)
       ! CALL DCOPY(2*NSTATE*NSTATE,A,1,U_TEMP,1)

       ! ELSE
       ! update U_MAT_new =  U_MAT_old + dt Up_MAT
       DO i=1,nstate
          DO j=1,nstate
             u_temp(j,i)=u_mat(j,i)+delta*u_grad(j,i)
          ENDDO
       ENDDO
       ! ENDIF
       ! DO I = 1,NSTATE
       ! write(6,'(I5,16(2f10.4))') I,(U_TEMP(I,J),J=1,8)
       ! ENDDO
       ! stop

       CALL lagrmult(u_mat,u_temp,a,b1,b2,gam,nstate)

    ENDIF  ! LAGRANGE OR NOT

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

    CALL zeroing(u_grad)!,SIZE(u_grad))

    IF (glocal%tgwannier) THEN
       CALL calc_u_wgrad_sum(u_mat,u_grad,xyzmat,nstate)
       igrmax_r=izamax(nstate*nstate,u_grad,1)
       grmax_r = ABS(zgive(u_grad,igrmax))
    ENDIF

    CALL calc_u_grad_sum(u_mat,u_grad,z_mat,nstate,nstate)

    igrmax=izamax(nstate*nstate,u_grad,1)
    grmax = ABS(zgive(u_grad,igrmax))

    dif_gfun = ABS(gloc_re%ofun - ofun00)
    ofun00 = gloc_re%ofun

    IF (glocal%tgwannier) THEN
       dif_xfun = ABS(gloc_re%xyzfun - xyzfun00)
       xyzfun00 = gloc_re%xyzfun

       gloc_re%dif_fun = glocr%g2g_weight*dif_gfun + glocr%wan_weight*dif_xfun
    ELSE
       gloc_re%dif_fun = dif_gfun
    ENDIF

    IF (glocal%tg_antisymm .AND. icont .LE. 100 )  THEN
       icont = icont+1
       IF (paral%io_parent)&
            WRITE(10,'(I10,3(2x,1PE12.6),2x,1PE12.6,2x,1PE12.6,2x,f10.5)')&
            nstep,gloc_re%gmax,gloc_re%ofun,grmax,gloc_re%dif_fun,gloc_re%omega_tot,step_fac
       IF (paral%io_parent)&
            WRITE(12,'(I10,16f10.4)') nstep,(omega_n(i),i=1,16)
    ELSEIF (glocal%tgwannier) THEN
       IF (paral%io_parent)&
            WRITE(10,'(I10,6(2x,1PE12.6),1PE12.3)')&
            nstep,gloc_re%gmax,gloc_re%ofun,gloc_re%xyzfun,grmax,gloc_re%dif_fun,gloc_re%ggnorm,&
            LOG(step_fac)
    ELSE
       IF (paral%io_parent)&
            WRITE(10,'(I10,3(2x,1PE12.6),2x,1PE12.6,2x,1PE12.6,1PE12.3)')&
            nstep,gloc_re%gmax,gloc_re%ofun,grmax,gloc_re%dif_fun,gloc_re%ggnorm,LOG(step_fac)
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
    CALL tihalt(' UPD_EXP_S',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE u_upd_exp_sum
  ! ==================================================================
  SUBROUTINE  calc_u_grad_sum(u_mat,u_grad,z_mat,ldx,nstate)
    ! ==--------------------------------------------------------------==
    ! ==  CALCULATES THE GRADIENT OF THE LOCALIZATION FUNCTIONAL      ==
    ! ==  WITH RESPECT TO THE UNITARY MATRIX COEFFICIENTS             ==
    ! ==  U_GRAD_ij = -d OFUNC / d U_ij                               ==
    ! ==--------------------------------------------------------------==
    ! input
    INTEGER                                  :: ldx, nstate
    COMPLEX(real_8)                          :: z_mat(nstate,nstate,3), &
                                                u_grad(nstate,nstate), &
                                                u_mat(nstate,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'calc_u_grad_sum'

    COMPLEX(real_8)                          :: zone, zzero
    COMPLEX(real_8), ALLOCATABLE             :: aux(:,:), uz(:,:), &
                                                z_matold(:,:,:), zu(:,:)
    INTEGER                                  :: i, ierr, isub, j, k, l, lmat

! (NSTATE,NSTATE,3)
! ==--------------------------------------------------------------==

    CALL tiset('  GRAD_SUM',isub)
    lmat = 2*nstate*nstate
    zzero = CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
    zone  = CMPLX(1.0_real_8,0.0_real_8,kind=real_8)

    ALLOCATE(z_matold(nstate,nstate,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL dcopy(lmat*3,z_mat,1,z_matold,1)

    CALL z_update(u_mat,z_mat,nstate,nstate,3)

    IF (glocr%g2g_weight .EQ. 0.0_real_8)  GOTO 100

    ALLOCATE(zu(nstate, nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(uz(nstate, nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    DO k = 1,3
       CALL zeroing(zu)!,SIZE(zu))
       CALL zeroing(uz)!,SIZE(uz))

       DO j = 1,nstate
          DO i = 1,nstate
             DO l = 1,nstate
                zu(i,j) = zu(i,j) + z_matold(i,l,k)*u_mat(l,j)
             ENDDO
          ENDDO

          DO i = 1,nstate
             DO l = 1,nstate
                uz(i,j) = uz(i,j) + CONJG(z_matold(l,i,k))*u_mat(l,j)
             ENDDO
          ENDDO

          DO i = 1,nstate
             u_grad(i,j) = u_grad(i,j)+glocr%g2g_weight*(&
                  zu(i,j)*CONJG(z_mat(j,j,k)) +&
                  uz(i,j)*z_mat(j,j,k))
          ENDDO
       ENDDO
    ENDDO

100 CONTINUE

    ALLOCATE(aux(nstate,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL zgemm('C','N',nstate,nstate,nstate,zone,u_grad,nstate,u_mat,&
         nstate,zzero,aux,nstate)
    CALL zgemm('N','N',nstate,nstate,nstate,-zone,u_mat,nstate,aux,&
         nstate,zone,u_grad,nstate)
    DEALLOCATE(aux,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)

    CALL gzfunc(z_mat,nstate,nstate)
    CALL dcopy(lmat*3,z_matold,1,z_mat,1)
    ! ==--------------------------------------------------------------==
    DEALLOCATE(zu,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(uz,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL tihalt('  GRAD_SUM',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE calc_u_grad_sum
  ! ==--------------------------------------------------------------==
  ! ==--------------------------------------------------------------==
  ! ==================================================================

END MODULE u_upd_exp_sum_utils


SUBROUTINE simm_op_sum(c1,u_mat,omega_n,omega_tot,nstate)
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE system , ONLY: ncpw,nkpt
  USE parac, ONLY : paral,parai
  IMPLICIT NONE
  REAL(real_8)                               :: omega_tot
  INTEGER                                    :: nstate
  REAL(real_8)                               :: omega_n(nstate), &
                                                u_mat(2,nstate,nstate), &
                                                c1(2,nkpt%ngwk,nstate)

  INTEGER                                    :: ig, istate, jstate
  REAL(real_8)                               :: ImA, ImB, ImO, ReA, ReB, ReO

  omega_tot = 0.0_real_8
  DO istate = 1,nstate
     omega_n(istate) = 0.0_real_8
     DO ig = 1,ncpw%ngw
        ReA = 0.0_real_8
        ImA = 0.0_real_8
        ReB = 0.0_real_8
        ImB = 0.0_real_8
        DO jstate = 1,nstate
           ReA = rea + c1(1,ig,jstate)*u_mat(1,jstate,istate)-&
                c1(2,ig,jstate)*u_mat(2,jstate,istate)
           ImA = ima + c1(1,ig,jstate)*u_mat(2,jstate,istate)+&
                c1(2,ig,jstate)*u_mat(1,jstate,istate)

           ReB = reb + c1(1,ig+ncpw%ngw,jstate)*u_mat(1,jstate,istate)-&
                c1(2,ig+ncpw%ngw,jstate)*u_mat(2,jstate,istate)
           ImB = imb + c1(1,ig+ncpw%ngw,jstate)*u_mat(2,jstate,istate)+&
                c1(2,ig+ncpw%ngw,jstate)*u_mat(1,jstate,istate)
        ENDDO
        ReO = ReA*ReB + ImA*ImB
        ImO = ReA*ImB - ImA*ReB
        omega_n(istate) = omega_n(istate) +&
             ReO*reo + ImO*imo
     ENDDO
     omega_tot = omega_tot + omega_n(istate)
  ENDDO

  RETURN
END SUBROUTINE simm_op_sum
! ==================================================================
