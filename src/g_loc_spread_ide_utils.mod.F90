MODULE g_loc_spread_ide_utils
  USE cnst,                            ONLY: fbohr
  USE cppt,                            ONLY: gk
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_def
  USE g_loc,                           ONLY: &
       filcen, filpen, filspr, gloc_re, glocal, gloci, glocr, lag_mul, nstep, &
       omega_n
  USE g_loc_util_utils,                ONLY: &
       gloc_print, gzloc_center, lagrmult, line_search_s, r_matrix, &
       spreadfunc, u_by_ortho, unitary_trans, w_matrix, z12_update, &
       zexponentiate
  USE geq0mod,                         ONLY: geq0
  USE kinds,                           ONLY: int_2,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8
  USE metr,                            ONLY: metr_com
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum,&
                                             mp_sync
  USE parac,                           ONLY: parai,&
                                             paral
  USE prng_utils,                      ONLY: repprngu
  USE soft,                            ONLY: soft_com
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             ncpw,&
                                             nkpt,&
                                             parm,&
                                             spar
  USE testex_utils,                    ONLY: testex
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: zclean_k,&
                                             zgive
  USE zeroing_utils,                   ONLY: zeroing
#include "sizeof.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: g_loc_spread_ide
  PUBLIC :: give_scr_g_loc_spread


  PUBLIC :: g2_mat
  PUBLIC :: g1_mat
  PUBLIC :: give_scr_upd_spread_sum
  PUBLIC :: give_scr_upd_spread_ide
  PUBLIC :: u_upd_spread_ide
  PUBLIC :: u_grad_spread_ide



CONTAINS


  ! ==================================================================
  SUBROUTINE u_upd_spread_ide(u_mat,z1_mat,z2_mat,u_grad,nstate)
    ! ==--------------------------------------------------------------==
    ! == CALCULATES THE UNITARY TRANSFORMATION 
    ! ==            WHICH OPTIMISE LOCALIZATION IN G SPACE            ==
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: nstate
    COMPLEX(real_8), TARGET                  :: u_grad(nstate,nstate)
    COMPLEX(real_8)                          :: z2_mat(nstate,nstate,3), &
                                                z1_mat(nstate,nstate,3), &
                                                u_mat(nstate,nstate)

    CHARACTER(*), PARAMETER :: procedureN = 'u_upd_spread_ide'

    COMPLEX(real_8)                          :: carg, zone, zzero
    COMPLEX(real_8), ALLOCATABLE             :: a(:,:), b1(:,:), b2(:,:), &
                                                gam(:,:)
    COMPLEX(real_8), POINTER                 :: u_temp(:,:)
    INTEGER                                  :: i, ierr, igrmax, isub, &
                                                izamax, j, l_a, l_b, l_gam
    INTEGER, SAVE                            :: icont = 0
    REAL(real_8)                             :: delta, difgam, grmax
    REAL(real_8), SAVE                       :: step_fac = 1.0_real_8

    CALL tiset(' UPD_SPR_E',isub)

    IF (cntl%tlsd)CALL stopgm('U_UPD_SPREAD_IDE','LSD NOT IMPLEMENTED YET',& 
         __LINE__,__FILE__)

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

    u_temp => u_grad

    zone  = CMPLX(1.0_real_8,0.0_real_8,kind=real_8)
    zzero = CMPLX(0.0_real_8,0.0_real_8,kind=real_8)

    ! time step
    delta = glocr%gloc_step*step_fac

    ! this is U_MAT

    ! Up_MAT is zero if first step
    ! CALL AZZERO(A,L_A)
    ! CALL AZZERO(B1,L_B)
    ! CALL AZZERO(B2,L_B)

    nstep = nstep + 1

    igrmax=izamax(nstate*nstate,u_grad,1)
    grmax = ABS(zgive(u_grad,igrmax))

    IF (gloci%gloc_const .EQ. 3) THEN

       ! UPDATE U_MAT WITH APPROXIMATION U = 1 + A (A antihermitian)
       CALL approx_u(a,u_grad,nstate)
       CALL dcopy(2*nstate*nstate,a,1,u_temp,1)

    ELSEIF (gloci%gloc_const .EQ. 2) THEN

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
    ! comp = GMAX - GMAX0
    ! GMAX0 = GMAX

    CALL dcopy(2*nstate*nstate,u_temp(1,1),1,u_mat(1,1),1)

    ! write(15,'(/,A,I5)') 'NSTEP ',NSTEP
    ! DO I = 1,NSTATE
    ! CC = CMPLX(0.0_real_8,0.0_real_8)
    ! DO J = 1,NSTATE
    ! CC = CC + CMPLX(real(U_MAT(J,5))*real(U_MAT(J,I))+
    ! &        aimag(U_MAT(J,5))*aimag(U_MAT(J,I)),
    ! &        real(U_MAT(J,I))*aimag(U_MAT(J,5))-
    ! &        aimag(U_MAT(J,I))*real(U_MAT(J,5)))
    ! ENDDO
    ! write(15,'(I5,2(1PE16.6))') I,CC
    ! ENDDO
    ! stop

    ! calculate U_GRAD_new with the rotated U_MAT_new
    CALL z12_update(u_mat,z1_mat,z2_mat,nstate,nstate)
    CALL spreadfunc(z1_mat,z2_mat,nstate,nstate)

    IF (glocal%tg_antisymm .AND. icont .LE. 100  )  THEN
       icont = icont+1
       IF (paral%io_parent)&
            WRITE(10,'(I10,3(2x,1PE12.6),I6,2x,1PE12.6)') nstep,gloc_re%gmax,gloc_re%ofun,&
            grmax,igrmax,gloc_re%omega_tot
       IF (paral%io_parent)&
            WRITE(12,'(I10,16f10.4)') nstep,(omega_n(i),i=1,16)
    ELSE
       IF (paral%io_parent)&
            WRITE(10,'(I10,3(2x,1PE12.6),I6,2x,1PE12.6,2x,1PE12.6)')&
            nstep,gloc_re%gmax,gloc_re%ofun,grmax,igrmax,gloc_re%ggnorm,step_fac/10._real_8
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
    CALL tihalt(' UPD_SPR_I',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE u_upd_spread_ide
  ! ==================================================================
  SUBROUTINE  u_grad_spread_ide(u_mat,u_grad,z1_mat,nstate)
    ! ==--------------------------------------------------------------==
    ! ==  CALCULATES THE GRADIENT OF THE LOCALIZATION FUNCTIONAL      ==
    ! ==  WITH RESPECT TO THE UNITARY MATRIX COEFFICIENTS             ==
    ! ==  U_GRAD_ij = -d OFUNC / d U_ij                               ==
    ! ==--------------------------------------------------------------==
    ! input
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: z1_mat(nstate,nstate,3), &
                                                u_grad(nstate,nstate), &
                                                u_mat(nstate,nstate)

    INTEGER                                  :: i, j, k

! ==--------------------------------------------------------------==

    DO k = 1,3

       DO j = 1,nstate
          DO i = 1,nstate

             u_grad(i,j) = u_grad(i,j)&
                  +2.0_real_8*z1_mat(i,j,k)*z1_mat(j,j,k)

          ENDDO
       ENDDO
    ENDDO
    ! DO I = 1,NSTATE
    ! write(6,'("UG",I6,10f12.8)') I,(U_GRAD(I,J),J=1,5)
    ! ENDDO
    ! stop

    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE u_grad_spread_ide
  ! ==--------------------------------------------------------------==
  ! ==================================================================


  ! ==--------------------------------------------------------------==
  SUBROUTINE give_scr_upd_spread_ide(ltranupd,tag,nstate)
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: ltranupd
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate

    INTEGER                                  :: l_multi, l_ortho, l_print, &
                                                l_search, lgrad, lznum

    lznum = 2*nstate*nstate
    ! ==--------------------------------------------------------------==
    lgrad     = 2*nstate*nstate
    l_ortho   = lgrad+3*2*nstate*nstate + 11*nstate
    l_multi   = lgrad+6*2*nstate*nstate
    IF (glocal%tg_linesearch) THEN
       l_search   = lgrad + 12*2*nstate*nstate
    ELSE
       l_search   = 1
    ENDIF

    IF (glocal%tglocrealp) THEN
       l_print = fpar%kr1*fpar%kr2s*fpar%kr3s + fpar%kr2s*fpar%kr3s
    ELSE
       l_print =1
    ENDIF
    ! ==--------------------------------------------------------------==
    ltranupd  = MAX(l_ortho,l_multi,lznum,l_print,l_search)
    tag='MAX(L_ORTHO,L_MULTI,LZNUM,...)'
    RETURN
  END SUBROUTINE give_scr_upd_spread_ide
  ! ==================================================================


  ! ==--------------------------------------------------------------==
  SUBROUTINE give_scr_upd_spread_sum(ltranupd,tag,nstate)
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: ltranupd
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate

    INTEGER                                  :: l_grad_sum, l_multi, l_ortho, &
                                                l_print, l_search, lgrad, &
                                                lznum

    lznum = 2*nstate*nstate
    ! ==--------------------------------------------------------------==
    lgrad      = 2*nstate*nstate
    l_grad_sum = lgrad + 4*2*nstate*nstate
    l_ortho    = lgrad + 3*2*nstate*nstate + 11*nstate
    l_multi    = lgrad + 6*2*nstate*nstate
    IF (glocal%tg_linesearch) THEN
       l_search   = lgrad + 12*2*nstate*nstate
    ELSE
       l_search   = 1
    ENDIF
    IF (glocal%tglocrealp) THEN
       l_print = fpar%kr1*fpar%kr2s*fpar%kr3s + fpar%kr2s*fpar%kr3s
    ELSE
       l_print =1
    ENDIF
    ! ==--------------------------------------------------------------==

    ltranupd  = MAX(l_grad_sum,l_ortho,l_multi,lznum,&
         l_print, l_search)

    tag='MAX(L_GRAD_SUM,L_ORTHO,...)'
    RETURN
  END SUBROUTINE give_scr_upd_spread_sum
  ! ==--------------------------------------------------------------==


  ! ==================================================================
  SUBROUTINE g_loc_spread_ide(c1,c2,nstate,nfi)
    ! ==--------------------------------------------------------------==
    ! == CALCULATES LOCALIZED FUNCTIONS IN G SPACE                    ==
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8)                          :: c1(nkpt%ngwk,*)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c2(nkpt%ngwk,nstate)
    INTEGER                                  :: nfi

    CHARACTER(*), PARAMETER :: procedureN = 'g_loc_spread_ide'

    COMPLEX(real_8), ALLOCATABLE             :: aux(:,:), gmat(:,:), &
                                                u_grad(:,:)
    COMPLEX(real_8), ALLOCATABLE, SAVE       :: rotlsd(:), rotmat(:,:), &
                                                z1_mat(:,:,:), z2_mat(:,:,:)
    INTEGER                                  :: i, ierr, irep, isub, j, &
                                                l_lag_mul, maxrep
    LOGICAL                                  :: ferror
    REAL(real_8)                             :: fac
    REAL(real_8), ALLOCATABLE                :: center(:,:), MODULO(:)

! ==--------------------------------------------------------------==
! ..   localized orbitals functions
! ==--------------------------------------------------------------==

    ALLOCATE(center(6,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(MODULO(ncpw%ngw),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    l_lag_mul  = nstate*nstate
    ALLOCATE(lag_mul(l_lag_mul),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(lag_mul)!,SIZE(lag_mul))

    IF (glocal%tg_antisymm)  THEN
       ALLOCATE(omega_n(nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(omega_n)!,nstate)
    ENDIF

    CALL tiset('G_LOC_SP_I',isub)

    gloc_re%gmax = 0.0_real_8

    maxrep = gloci%gloc_maxs
    nstep  = gloci%gloc_init


    ! ..   initialization for localized function part

    ALLOCATE(z1_mat(nstate,nstate,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(z2_mat(nstate,nstate,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(rotmat(nstate,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)


    IF (paral%parent) THEN
       IF (paral%io_parent)&
            CALL fileopen(14,filcen,fo_def,ferror)
       fac = parm%tpiba * fbohr
       DO i = 1,2*nstate
          IF (paral%io_parent)&
               WRITE(14,'(A,3f12.5,I10)') 'H',gk(1,i)*fac,gk(2,i)*fac,&
               gk(3,i)*fac,i
          IF (paral%io_parent)&
               WRITE(14,'(A,3f12.5,I10)') 'H',-gk(1,i)*fac,-gk(2,i)*fac,&
               -gk(3,i)*fac,i+ncpw%ngw
       ENDDO
    ENDIF

    IF (glocal%tg_read_matrix)  THEN
       IF (paral%parent) CALL r_matrix(rotmat,nstate)
       CALL mp_bcast(rotmat,SIZE(rotmat),parai%source,parai%allgrp)
       ! rotate the wavefunctions to the localized representation

       CALL rotate_c((1._real_8,0.0_real_8),c1,(0._real_8,0.0_real_8),c2,rotmat,nstate)
       CALL dcopy(2*2*ncpw%ngw*nstate,c2,1,c1,1)
       IF (geq0)   CALL zclean_k(c1,nstate,ncpw%ngw)
    ENDIF

    ! CALL PLANE_WAVES_SHELL(C1,NSTATE)
    ! stop

    ! calculate the centers and spread of the wannier functions
    CALL gzloc_center(c1,nstate,center,modulo)
    IF (paral%io_parent)&
         WRITE(6,&
         ' (/,5x,A,/)')&
         '_______CENTERS AND SPREAD BEFORE THE OPTIMIZATION_______'
    IF (paral%parent) CALL gloc_print(nstate,center,0)


    ! ..   calculation of 3 matrices :  <n |O_i| m>
    ! where n and m address the electronic orbitals and i = x,y,z

    IF (glocal%tg_kick) THEN

       CALL unitary_trans(c1,c2,nstate)
       IF (paral%io_parent)&
            WRITE(6,'(/,A)') '     FIRST RANDOM ROTATION    '
       ! calculate the centers and spread of the wannier functions
       CALL gzloc_center(c1,nstate,center,modulo)
       IF (paral%io_parent) WRITE(6,' (/,5x,A,/)')&
            '_____ CENTERS AND SPREAD AFTER FIRST RANDOM ROTATION _____'
       IF (paral%io_parent) CALL gloc_print(nstate,center,0)
    ENDIF


    ! HERE THE LOOPS START:
    ! at each step Z_MAT is calculated, then the functional is estimated
    ! 
    ! The ROTMAT matrix is assumed as the identical matrix 
    ! 
    ! The gradient of the functional U_GRAD is calculated assuming a small
    ! perturbation of ROTMAT with respect to the identical matrix
    ! U_GRADij = Z_MATij * CONJG(Z_MATjj) + CONJG(Z_MATij) * Z_MATjj 
    ! 
    ! the ROTMAT is updated mantaining the unitary condition

    ALLOCATE(aux(nstate,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    CALL g2_mat(c1,c2,nstate,aux,1)
    CALL dcopy(2*nstate*nstate,aux,1,z2_mat(1,1,1),1)

    CALL g2_mat(c1,c2,nstate,aux,2)
    CALL dcopy(2*nstate*nstate,aux,1,z2_mat(1,1,2),1)

    CALL g2_mat(c1,c2,nstate,aux,3)
    CALL dcopy(2*nstate*nstate,aux,1,z2_mat(1,1,3),1)

    CALL g1_mat(c1,c2,nstate,aux,1)
    CALL dcopy(2*nstate*nstate,aux,1,z1_mat(1,1,1),1)

    CALL g1_mat(c1,c2,nstate,aux,2)
    CALL dcopy(2*nstate*nstate,aux,1,z1_mat(1,1,2),1)

    CALL g1_mat(c1,c2,nstate,aux,3)
    CALL dcopy(2*nstate*nstate,aux,1,z1_mat(1,1,3),1)

    DEALLOCATE(aux,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)

    CALL mp_sum(z1_mat,3*nstate*nstate,parai%allgrp)
    CALL mp_sum(z2_mat,3*nstate*nstate,parai%allgrp)
    CALL mp_sync(parai%allgrp)



    IF (paral%parent) THEN
       CALL spreadfunc(z1_mat,z2_mat,nstate,nstate)
       IF (paral%io_parent)&
            WRITE(6,'(/,A,F12.6)') 'STARTING FUNCTIONAL VALUE= ',gloc_re%ofun
    ENDIF

    ALLOCATE(u_grad(nstate, nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    IF (paral%parent) THEN
       IF (paral%io_parent)&
            CALL fileopen(10,filspr,fo_def,ferror)
       IF (paral%io_parent)&
            CALL fileopen(11,'umat.dat',fo_def,ferror)
       IF (glocal%tg_antisymm.AND.paral%io_parent)&
            CALL fileopen(12,filpen,fo_def,ferror)
       ! initialize U = ROTMAT matrix
       CALL zeroing(rotmat)!,SIZE(rotmat))
       DO i=1,nstate
          rotmat(i,i)=CMPLX(1._real_8,0.0_real_8,kind=real_8)
       ENDDO
       ALLOCATE(gmat(nstate, nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       CALL zeroing(gmat)!,SIZE(gmat))
       IF (glocr%gloc_ran.GT.0._real_8) THEN
          ! ..   randomly distort initial vectors
          DO i=1,nstate
             gmat(i,i)=CMPLX(1.0_real_8,0._real_8,kind=real_8)
             DO j=i+1,nstate
                gmat(i,j)=CMPLX(repprngu()*glocr%gloc_ran,repprngu()*glocr%gloc_ran,kind=real_8)
                gmat(j,i)=CONJG(gmat(i,j))
             ENDDO
          ENDDO
          CALL zexponentiate(gmat,rotmat,nstate)
          glocr%gloc_ran=0._real_8
       ENDIF
    ENDIF

    ! stop      
    DO irep = 1,maxrep


       CALL zeroing(u_grad)!,SIZE(u_grad))

       IF (glocal%tg_antisymm .AND. irep .LE. 100) THEN
          CALL simm_op(c1,omega_n,gloc_re%omega_tot,nstate)
          IF (glocal%tg_penalty) THEN
             CALL omega_grad(c1,u_grad,nstate)
          ENDIF

          CALL mp_sum(u_grad,nstate*nstate,parai%allgrp)
          CALL mp_sum(omega_n,nstate,parai%allgrp)
          CALL mp_sum(gloc_re%omega_tot,parai%allgrp)

          CALL mp_sync(parai%allgrp)
       ENDIF
       ! calculate the  gradient
       IF (paral%parent) THEN
          CALL u_grad_spread_ide(rotmat,u_grad,z1_mat,nstate)


          CALL u_upd_spread_ide(rotmat,z1_mat,z2_mat,u_grad,&
               nstate)

          ! write(11,*) 'NSTEP = ',NSTEP
          ! DO I = 1,NSTATE
          ! write(11,'(I3,8(2f8.4))') I,(ROTMAT(I,J),J=1,8)
          ! ENDDO
          ! write(11,'(/,30A,/)') '**'
          ! stop            
       ENDIF

       CALL mp_sync(parai%allgrp)
       CALL mp_bcast(rotmat,SIZE(rotmat),parai%source,parai%allgrp)
       CALL mp_bcast_byte(gloc_re,size_in_bytes_of(gloc_re),parai%source,parai%allgrp)
       ! rotate the wavefunctions to the localized representation

       CALL rotate_c((1._real_8,0.0_real_8),c1,(0._real_8,0.0_real_8),c2,rotmat,nstate)
       CALL dcopy(2*2*ncpw%ngw*nstate,c2,1,c1,1)
       IF (geq0)   CALL zclean_k(c1,nstate,ncpw%ngw)

       IF (paral%parent) THEN
          ! initialize U = ROTMAT matrix
          CALL zeroing(rotmat)!,nstate*nstate)
          DO i=1,nstate
             rotmat(i,i)=CMPLX(1._real_8,0.0_real_8,kind=real_8)
          ENDDO
       ENDIF

       IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)
       IF (gloc_re%gmax.LT.glocr%gloc_eps) THEN
          GOTO 400
       ELSEIF (ABS(gloc_re%dif_fun) .LT. 1.e-10_real_8) THEN
          IF (paral%parent) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(/,A,1PE12.6)') '     OFUN-OFUN0 = ',&
                  gloc_re%dif_fun
             IF (paral%io_parent)&
                  WRITE(6,'(A,/)') '   G_LOC_UPDATE: CONVERGENCE'
          ENDIF
          GOTO 400
       ELSEIF (irep .EQ. maxrep .OR. soft_com%exsoft) THEN
          GOTO 500
       ENDIF               ! GMAX > EPS
    ENDDO                     ! IREP

500 IF (soft_com%exsoft.AND.paral%io_parent)&
         WRITE(6,'(A,1E12.6)')&
         ' G_LOC_SPREAD_IDE|SOFTWERE EXIT (OFUNC)',gloc_re%ofun
    IF ((irep .EQ. maxrep).AND.paral%io_parent)&
         WRITE(6,'(A,I10)')&
         ' G_LOC_SPREAD_IDE| MAXREP REACHED ',maxrep
400 CONTINUE
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            CALL fileclose(10)
       IF (paral%io_parent)&
            CALL fileclose(11)
       IF ((glocal%tg_antisymm).AND.paral%io_parent)&
            CALL fileclose(12)
    ENDIF

    ! calculate the centers and spread of the wannier functions
    CALL gzloc_center(c1,nstate,center,modulo)
    IF (paral%parent) CALL gloc_print(nstate,center,1)
    IF (paral%io_parent)&
         CALL fileclose(14)


    ! SAVE THE ROTATION MATRIX
    IF (paral%parent) CALL w_matrix(rotmat,nstate)

    DEALLOCATE(lag_mul,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (glocal%tg_antisymm) DEALLOCATE(omega_n,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(z1_mat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(z2_mat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rotmat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (cntl%tlsd) DEALLOCATE(rotlsd,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    IF (paral%parent) THEN
       DEALLOCATE(gmat,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    ENDIF

    DEALLOCATE(center,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(modulo,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    CALL tihalt('G_LOC_SP_I',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE g_loc_spread_ide
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE give_scr_g_loc_spread(lgloc,tag,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lgloc
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate

    INTEGER                                  :: lgzloc

    lgzloc = 1
    IF (gloci%gloc_opt .EQ. 1) THEN
       CALL give_scr_upd_spread_sum(lgloc,tag,nstate)
    ELSEIF (gloci%gloc_opt .EQ. 2) THEN
       CALL give_scr_upd_spread_ide(lgloc,tag,nstate)
    ENDIF
    RETURN
  END SUBROUTINE give_scr_g_loc_spread
  ! ==--------------------------------------------------------------==
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE g2_mat(c0,c2,nstate,ddmat,index)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(2*ncpw%ngw,*), &
                                                c2(2*ncpw%ngw,*)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: ddmat(nstate,nstate)
    INTEGER                                  :: index

    INTEGER                                  :: i, ig, ii
    REAL(real_8)                             :: fac, numsec(3)

! ==--------------------------------------------------------------==

    CALL zeroing(ddmat)!,nstate*nstate)

    numsec(1) = REAL(spar%nr1s,kind=real_8)
    numsec(2) = REAL(spar%nr2s,kind=real_8)
    numsec(3) = REAL(spar%nr3s,kind=real_8)

    fac = metr_com%ht(index,1)*metr_com%ht(index,1)+metr_com%ht(index,2)*metr_com%ht(index,2)+metr_com%ht(index,3)*&
         metr_com%ht(index,3)
    fac = 4.0_real_8*dacos(-1.0_real_8)*dacos(-1.0_real_8)/fac
    fac = numsec(index)*numsec(index)*fac
    fac = 1.0_real_8/fac


    DO i= 1,nstate

       ii=i

       DO ig=1,ncpw%ngw
          c2(ig,ii)     = gk(index,ig)*gk(index,ig) * parm%tpiba2* c0(ig,ii)
          c2(ig+ncpw%ngw,ii) = gk(index,ig)*gk(index,ig) * parm%tpiba2* c0(ig+ncpw%ngw,&
               ii)
       ENDDO
       IF (geq0)c2(1+ncpw%ngw,ii)=CMPLX(0._real_8,0._real_8,kind=real_8)
    ENDDO

    CALL zgemm('C','N',nstate,nstate,2*ncpw%ngw,CMPLX(1._real_8,0._real_8,kind=real_8),c0,2*ncpw%ngw,&
         c2,2*ncpw%ngw,CMPLX(0._real_8,0._real_8,kind=real_8),ddmat,nstate)

    CALL dscal(2*nstate*nstate,fac,ddmat,1)

    CALL mp_sum(ddmat,nstate*nstate,parai%allgrp)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE g2_mat
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE g1_mat(c0,c2,nstate,ddmat,index)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(2*ncpw%ngw,*), &
                                                c2(2*ncpw%ngw,*)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: ddmat(nstate,nstate)
    INTEGER                                  :: index

    INTEGER                                  :: i, ig, ii
    REAL(real_8)                             :: fac, numsec(3)

! ==--------------------------------------------------------------==

    CALL zeroing(ddmat)!,nstate*nstate)

    numsec(1) = REAL(spar%nr1s,kind=real_8)
    numsec(2) = REAL(spar%nr2s,kind=real_8)
    numsec(3) = REAL(spar%nr3s,kind=real_8)

    fac = metr_com%ht(index,1)*metr_com%ht(index,1)+metr_com%ht(index,2)*metr_com%ht(index,2)+metr_com%ht(index,3)*&
         metr_com%ht(index,3)
    fac = 2.0_real_8*dacos(-1.0_real_8)/SQRT(fac)
    fac = numsec(index)*fac
    fac = 1.0_real_8/fac


    DO i=1,nstate

       ii=i

       DO ig=1,ncpw%ngw
          c2(ig,ii)     = gk(index,ig) * parm%tpiba* c0(ig,ii)
          c2(ig+ncpw%ngw,ii) = -gk(index,ig) * parm%tpiba* c0(ig+ncpw%ngw,ii)
       ENDDO
       IF (geq0)c2(1+ncpw%ngw,ii)=CMPLX(0._real_8,0._real_8,kind=real_8)
    ENDDO

    CALL zgemm('C','N',nstate,nstate,2*ncpw%ngw,CMPLX(1._real_8,0._real_8,kind=real_8),c0,2*ncpw%ngw,&
         c2,2*ncpw%ngw,CMPLX(0._real_8,0._real_8,kind=real_8),ddmat,nstate)

    CALL dscal(2*nstate*nstate,fac,ddmat,1)

    CALL mp_sum(ddmat,nstate*nstate,parai%allgrp)

    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE g1_mat
  ! ==================================================================
  ! ==================================================================

END MODULE g_loc_spread_ide_utils


! ==================================================================
SUBROUTINE simm_op(c1,omega_n,omega_tot,nstate)
  ! ==--------------------------------------------------------------==
  ! ==    CALCULATE THE PENALTY FUNCTION 
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE system , ONLY: ncpw,nkpt
  USE parac, ONLY : paral,parai
  IMPLICIT NONE
  REAL(real_8)                               :: omega_tot
  INTEGER                                    :: nstate
  REAL(real_8)                               :: omega_n(nstate), &
                                                c1(2,nkpt%ngwk,nstate)

  INTEGER                                    :: ig, istate
  REAL(real_8)                               :: ImA, ImB, ImO, ReA, ReB, ReO

  omega_tot = 0.0_real_8
  DO istate = 1,nstate
     omega_n(istate) = 0.0_real_8
     DO ig = 1,ncpw%ngw
        ReA = c1(1,ig,istate)
        ImA = c1(2,ig,istate)
        ReB = c1(1,ig+ncpw%ngw,istate)
        ImB = c1(2,ig+ncpw%ngw,istate)
        ReO = ReA*ReB + ImA*ImB
        ImO = ReA*ImB - ImA*ReB
        omega_n(istate) = omega_n(istate) + ReO*reo + ImO*imo
     ENDDO
     omega_tot = omega_tot + omega_n(istate)
  ENDDO

  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE simm_op
! ==================================================================

! ==================================================================
SUBROUTINE  omega_grad(c1,u_grad,nstate)
  ! ==--------------------------------------------------------------==
  ! ==  CALCULATES THE GRADIENT OF THE SIMMETRY  FUNCTIONAL         ==
  ! ==  WITH RESPECT TO THE UNITARY MATRIX COEFFICIENTS             ==
  ! ==  U_GRAD_ij = -d OMEGA / d U_ij                               ==
  ! ==  the gradient is calculated in U = 1                         ==
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE system , ONLY: ncpw,nkpt
  USE parac, ONLY : paral,parai
  USE zeroing_utils,                   ONLY: zeroing
  IMPLICIT NONE
  INTEGER                                    :: nstate
  REAL(real_8)                               :: u_grad(2,nstate,nstate), &
                                                c1(2,nkpt%ngwk,nstate)

  CHARACTER(*), PARAMETER                    :: procedureN = 'omega_grad'

  INTEGER                                    :: i, ierr, ig, j
  REAL(real_8)                               :: im_ij, im_ji, im_jj, ImA_I, &
                                                ImA_J, ImB_I, ImB_J, re_ij, &
                                                re_ji, re_jj, ReA_I, ReA_J, &
                                                ReB_I, ReB_J
  REAL(real_8), ALLOCATABLE                  :: aux(:,:,:)

! input
! (2,NSTATE,NSTATE)
! ==--------------------------------------------------------------==

  ALLOCATE(aux(2,nstate,nstate),STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
       __LINE__,__FILE__)
  CALL zeroing(aux)!,2*nstate*nstate)

  DO j = 1,nstate
     DO i = 1,nstate
        DO ig = 1,ncpw%ngw
           ReA_I = c1(1,ig,i)
           ImA_I = c1(2,ig,i)
           ReB_I = c1(1,ig+ncpw%ngw,i)
           ImB_I = c1(2,ig+ncpw%ngw,i)
           ReA_J = c1(1,ig,j)
           ImA_J = c1(2,ig,j)
           ReB_J = c1(1,ig+ncpw%ngw,j)
           ImB_J = c1(2,ig+ncpw%ngw,j)

           re_ij = ReA_I*ReB_J + ImA_I*ImB_J
           im_ij = ReA_I*ImB_J - ImA_I*ReB_J

           re_ji = ReA_J*ReB_I + ImA_J*ImB_I
           im_ji = ReA_J*ImB_I - ImA_J*ReB_I

           re_jj = ReA_J*ReB_J + ImA_J*ImB_J
           im_jj = ReA_J*ImB_J - ImA_J*ReB_J

           aux(1,i,j) = aux(1,i,j) + re_ij*re_jj + im_ij*im_jj +&
                re_jj*re_ji + im_jj*im_ji
           aux(2,i,j) = aux(2,i,j) - re_ij*im_jj + im_ij*re_jj&
                -re_jj*im_ji + im_jj*re_ji
        ENDDO

        u_grad(1,i,j) = u_grad(1,i,j) - aux(1,i,j)
        u_grad(2,i,j) = u_grad(2,i,j) - aux(2,i,j)

     ENDDO
  ENDDO

  DEALLOCATE(aux,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
       __LINE__,__FILE__)
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE omega_grad
! ==--------------------------------------------------------------==
