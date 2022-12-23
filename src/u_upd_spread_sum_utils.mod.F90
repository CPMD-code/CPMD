MODULE u_upd_spread_sum_utils
  USE error_handling,                  ONLY: stopgm
  USE g_loc,                           ONLY: gifirst,&
                                             gloc_re,&
                                             glocal,&
                                             gloci,&
                                             glocr,&
                                             nstep,&
                                             omega_n
  USE g_loc_util_utils,                ONLY: lagrmult,&
                                             line_search_s,&
                                             spreadfunc,&
                                             u_by_ortho,&
                                             z12_update
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: paral
  USE system,                          ONLY: cntl
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: zgive
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: u_upd_spread_sum
  PUBLIC :: u_grad_spread_sum


CONTAINS

  ! ==================================================================
  SUBROUTINE u_upd_spread_sum(u_mat,z1_mat,z2_mat,u_grad,&
       nstate)
    ! ==--------------------------------------------------------------==
    ! == CALCULATES THE UNITARY TRANSFORMATION 
    ! ==            WHICH OPTIMISE LOCALIZATION IN G SPACE            ==
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: nstate
    COMPLEX(real_8), TARGET                  :: u_grad(nstate,nstate)
    COMPLEX(real_8)                          :: z2_mat(nstate,nstate,3), &
                                                z1_mat(nstate,nstate,3), &
                                                u_mat(nstate,nstate)

    CHARACTER(*), PARAMETER :: procedureN = 'u_upd_spread_sum'

    COMPLEX(real_8)                          :: carg, zone, zzero
    COMPLEX(real_8), ALLOCATABLE             :: a(:,:), b1(:,:), b2(:,:), &
                                                gam(:,:)
    COMPLEX(real_8), POINTER                 :: u_temp(:,:)
    INTEGER                                  :: i, ierr, igrmax, isub, &
                                                izamax, j, l_a, l_b, l_gam
    INTEGER, SAVE                            :: icont = 0
    REAL(real_8)                             :: delta, difgam, grmax
    REAL(real_8), SAVE                       :: ofun00 = 0.0_real_8, &
                                                step_fac = 1.0_real_8

    CALL tiset(' UPD_SPR_S',isub)

    IF (cntl%tlsd) CALL stopgm('U_UPD_SPREAD_SUM','LSD NOT IMPLEMENTED YET'&
         ,& 
         __LINE__,__FILE__)

    IF (gifirst .EQ. 1) THEN
       icont    = 0
       step_fac = 1.0_real_8
       ofun00   = 0.0_real_8
       gifirst   = 0
    ENDIF

    l_a        = 2*nstate*nstate
    l_b        = 2*nstate*nstate
    l_gam      = 2*nstate*nstate

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

    CALL zeroing(a)!,SIZE(a))
    CALL zeroing(b1)!,SIZE(b1))
    CALL zeroing(b2)!,SIZE(b2))
    CALL zeroing(gam)!,SIZE(gam))

    u_temp => u_grad

    zone  = CMPLX(1.0_real_8,0.0_real_8,kind=real_8)
    zzero = CMPLX(0.0_real_8,0.0_real_8,kind=real_8)

    ! time step
    delta = glocr%gloc_step*step_fac

    ! this is U_MAT
    ! Up_MAT is zero if first step
    nstep = nstep + 1

    igrmax=izamax(nstate*nstate,u_grad,1)
    grmax = ABS(zgive(u_grad,igrmax))

    IF (gloci%gloc_const .EQ. 2) THEN

       CALL dcopy(2*nstate*nstate,u_mat,1,a,1)

       IF (glocal%tg_linesearch) THEN
          CALL line_search_s(a,u_grad,z1_mat,z2_mat,b1,b2,&
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

    ENDIF

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
    CALL zeroing(u_grad)!,SIZE(u_grad))
    CALL u_grad_spread_sum(u_mat,u_grad,z1_mat,z2_mat,nstate,nstate)

    gloc_re%dif_fun = ABS(gloc_re%ofun - ofun00)
    ofun00 = gloc_re%ofun
    ! write(6,*) OFUN , GMAX,GRMAX,DIF_FUN
    ! 
    IF (glocal%tg_antisymm .AND. icont .LE. 100 )  THEN
       icont = icont+1
       IF (paral%io_parent)&
            WRITE(10,'(I10,3(2x,1PE12.6),I6,2x,1PE12.6)') nstep,gloc_re%gmax,gloc_re%ofun,&
            grmax,igrmax,gloc_re%omega_tot
       IF (paral%io_parent)&
            WRITE(12,'(I10,16f10.4)') nstep,(omega_n(i),i=1,16)
    ELSE
       IF (paral%io_parent)&
            WRITE(10,'(I10,3(2x,1PE12.6),2x,1PE12.6,2x,1PE12.6,2x,1PE12.3)')&
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
    CALL tihalt(' UPD_SPR_S',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE u_upd_spread_sum
  ! ==================================================================
  SUBROUTINE  u_grad_spread_sum(u_mat,u_grad,z1_mat,z2_mat,&
       ldx,nstate)
    ! ==--------------------------------------------------------------==
    ! ==  CALCULATES THE GRADIENT OF THE LOCALIZATION FUNCTIONAL      ==
    ! ==  WITH RESPECT TO THE UNITARY MATRIX COEFFICIENTS             ==
    ! ==  U_GRAD_ij = -d OFUNC / d U_ij                               ==
    ! ==--------------------------------------------------------------==
    ! input
    INTEGER                                  :: ldx, nstate
    COMPLEX(real_8) :: z2_mat(nstate,nstate,3), z1_mat(nstate,nstate,3), &
      u_grad(nstate,nstate), u_mat(nstate,nstate)

    CHARACTER(*), PARAMETER :: procedureN = 'u_grad_spread_sum'

    COMPLEX(real_8)                          :: zone, zzero
    COMPLEX(real_8), ALLOCATABLE             :: aux(:,:), z1_matold(:,:,:), &
                                                z1u(:,:), z2_matold(:,:,:)
    INTEGER                                  :: i, ierr, j, k, l, lmat

! ==--------------------------------------------------------------==

    lmat = 2*nstate*nstate
    zzero = CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
    zone  = CMPLX(1.0_real_8,0.0_real_8,kind=real_8)

    ALLOCATE(z1_matold(nstate,nstate,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(z2_matold(nstate,nstate,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    CALL dcopy(lmat*3,z1_mat,1,z1_matold,1)
    CALL dcopy(lmat*3,z2_mat,1,z2_matold,1)
    CALL z12_update(u_mat,z1_mat,z2_mat,ldx,nstate)

    ALLOCATE(z1u(nstate, nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    DO k = 1,3
       ! CALL ZGEMM('N','N',NSTATE,NSTATE,NSTATE,CMPLX(1.0_real_8,0.0_real_8),
       ! &              Z_MATold(1,1,K),NSTATE,U_MAT,NSTATE,
       ! &              CMPLX(0.0_real_8,0.0_real_8),ZU,NSTATE)

       CALL zeroing(z1u)!,SIZE(z1u))

       DO j = 1,nstate
          DO i = 1,nstate
             DO l = 1,nstate
                z1u(i,j) = z1u(i,j) + z1_matold(i,l,k)*u_mat(l,j)
             ENDDO


             ! CALL ZGEMM('N','N',NSTATE,NSTATE,NSTATE,CMPLX(1.0_real_8,0.0_real_8),
             ! &              U_MAT,NSTATE,Z_MATold(1,1,K),NSTATE,
             ! &              CMPLX(0.0_real_8,0.0_real_8),UZ,NSTATE)


             u_grad(i,j) = u_grad(i,j)&
                  +2.0_real_8*z1u(i,j)*z1_mat(j,j,k)
          ENDDO
       ENDDO
    ENDDO
    DEALLOCATE(z1u,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)

    ALLOCATE(aux(nstate,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL zeroing(aux)!, nstate*nstate)
    CALL zgemm('C','N',nstate,nstate,nstate,zone,u_grad,nstate,u_mat,&
         nstate,zzero,aux,nstate)
    CALL zgemm('N','N',nstate,nstate,nstate,-zone,u_mat,nstate,aux,&
         nstate,zone,u_grad,nstate)
    DEALLOCATE(aux,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)

    CALL spreadfunc(z1_mat,z2_mat,nstate,nstate)
    CALL dcopy(lmat*3,z1_matold,1,z1_mat,1)
    CALL dcopy(lmat*3,z2_matold,1,z2_mat,1)
    ! write(6,*)  OFUN

    ! ==--------------------------------------------------------------==
    DEALLOCATE(z1_matold,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(z2_matold,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE u_grad_spread_sum
  ! ==--------------------------------------------------------------==
  ! ==================================================================

END MODULE u_upd_spread_sum_utils
