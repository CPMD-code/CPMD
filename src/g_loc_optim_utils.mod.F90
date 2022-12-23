MODULE g_loc_optim_utils
  USE cnst,                            ONLY: fbohr
  USE cppt,                            ONLY: gk
  USE ddip,                            ONLY: lenbk,&
                                             ngwmax
  USE ddipo_utils,                     ONLY: setdip
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_def
  USE g_loc,                           ONLY: &
       filcen, filpen, filspr, gloc_re, glocal, gloci, glocr, lag_mul, nstep, &
       omega_n
  USE g_loc_util_utils,                ONLY: gloc_print,&
                                             gzloc_center,&
                                             r_matrix,&
                                             unitary_trans,&
                                             w_matrix,&
                                             z_update,&
                                             zexponentiate
  USE geq0mod,                         ONLY: geq0
  USE kinds,                           ONLY: int_2,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8
  USE localize_utils,                  ONLY: wc_print
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum,&
                                             mp_sync
  USE opeigr_utils,                    ONLY: give_scr_opeigr
  USE parac,                           ONLY: parai,&
                                             paral
  USE prng_utils,                      ONLY: repprngu
  USE soft,                            ONLY: soft_com
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             ncpw,&
                                             nkpt,&
                                             parap,&
                                             parm,&
                                             spar
  USE testex_utils,                    ONLY: testex
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE u_upd_exp_utils,                 ONLY: u_grad_exp_ide,&
                                             u_upd_exp_ide
  USE utils,                           ONLY: nxxfun,&
                                             zclean_k
  USE wann,                            ONLY: wannc,&
                                             wanni
  USE wannier_center_utils,            ONLY: wannier_center
  USE zeroing_utils,                   ONLY: zeroing
  USE znum_mat_utils,                  ONLY: give_scr_znum_mat,&
                                             znum_mat
#include "sizeof.h"

  IMPLICIT NONE

  PRIVATE

  !!public :: approx_u

  !public :: give_scr_g_loc_exp_sum
  PUBLIC :: g_loc_exp_ide
  !public :: give_scr_g_loc_exp
  !!public :: gzfunc

!!$public :: give_scr_glocalize
  !public :: give_scr_upd_exp_sum
  PUBLIC :: calc_u_wgrad_sum
  PUBLIC :: g_loc_xyzmat


CONTAINS



  ! ==================================================================
  SUBROUTINE give_scr_g_loc_exp_sum(lgloc,tag,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lgloc
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate

    INTEGER                                  :: l_xyzmat, lznum

    CALL give_scr_znum_mat(lznum,tag,nstate)
    CALL give_scr_upd_exp_sum(lgloc,tag,nstate)
    IF (glocal%tgwannier) THEN
       CALL give_scr_xyzmat(l_xyzmat,tag)
    ELSE
       l_xyzmat =1
    ENDIF
    lgloc=MAX(lgloc,lznum,l_xyzmat)
    tag='MAX(LZNUM,LUGRAD+6*2*N2,...)'
    ! write(6,*) 'LGLOC ' ,LGLOC,' LZNUM ' ,LZNUM
    ! stop
    RETURN
  END SUBROUTINE give_scr_g_loc_exp_sum
  ! ==--------------------------------------------------------------==


  SUBROUTINE give_scr_upd_exp_sum(ltranupd,tag,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ltranupd
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate

    INTEGER                                  :: l_grad_sum, l_multi, l_ortho, &
                                                l_print, l_search, &
                                                l_wgrad_sum, lgrad, lznum

    CALL give_scr_znum_mat(lznum,tag,nstate)
    ! ==--------------------------------------------------------------==
    lgrad      = 2*nstate*nstate
    l_grad_sum = lgrad+5*2*nstate*nstate
    l_ortho    = lgrad+3*2*nstate*nstate + 11*nstate
    l_multi    = lgrad+6*2*nstate*nstate
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
    IF (glocal%tgwannier) THEN
       l_wgrad_sum = lgrad+8*2*nstate*nstate
    ELSE
       l_wgrad_sum =1
    ENDIF
    ! write(6,*) L_GRAD_SUM,L_ORTHO,L_MULTI,LZNUM,L_SEARCH,L_PRINT
    ! ==--------------------------------------------------------------==
    ltranupd  = MAX(l_grad_sum,l_ortho,l_multi,lznum,l_search,&
         l_print,l_wgrad_sum)
    tag='MAX(L_GRAD_SUM,L_ORTHO,...)'
    RETURN
  END SUBROUTINE give_scr_upd_exp_sum
  ! ==--------------------------------------------------------------==

  ! ==================================================================
  SUBROUTINE g_loc_exp_ide(c1,c2,nstate,nfi)
    ! ==--------------------------------------------------------------==
    ! == CALCULATES LOCALIZED FUNCTIONS IN G SPACE                    ==
    ! ==--------------------------------------------------------------==

    ! input
    COMPLEX(real_8)                          :: c1(nkpt%ngwk,*)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c2(nkpt%ngwk,nstate)
    INTEGER                                  :: nfi

    CHARACTER(*), PARAMETER                  :: procedureN = 'g_loc_exp_ide'

    COMPLEX(real_8), ALLOCATABLE             :: ddmat(:,:), gmat(:,:), &
                                                u_grad(:,:), xyzmat(:,:,:)
    COMPLEX(real_8), ALLOCATABLE, SAVE       :: rotlsd(:), rotmat(:,:), &
                                                z_mat(:,:,:)
    INTEGER                                  :: i, ierr, irep, isub, j, &
                                                l_lag_mul, maxrep
    LOGICAL                                  :: ferror
    REAL(real_8)                             :: fac
    REAL(real_8), ALLOCATABLE                :: center(:,:), MODULO(:)

! ==--------------------------------------------------------------==
! ..   Localized orbitals functions
! ==--------------------------------------------------------------==

    CALL tiset('G_LOC_EX_I',isub)

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


    gloc_re%gmax = 0.0_real_8

    maxrep = gloci%gloc_maxs
    nstep  = gloci%gloc_init

    ! ..   initialization for localized function part

    ALLOCATE(z_mat(nstate,nstate,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(rotmat(nstate,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! WARNING this is just a guess
    ALLOCATE(xyzmat(nstate,nstate,2*3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! WARNING


    IF (paral%io_parent) THEN
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

    ALLOCATE(center(6,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(MODULO(ncpw%ngw),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL gzloc_center(c1,nstate,center,modulo)
    IF (paral%io_parent) WRITE(6,&
         ' (/,5x,A,/)')&
         '_______CENTERS AND SPREAD BEFORE THE OPTIMIZATION_______'
    IF (paral%parent) CALL gloc_print(nstate,center,0)


    IF (glocal%tg_kick) THEN

       CALL unitary_trans(c1,c2,nstate)
       IF (paral%io_parent) WRITE(6,'(/,A)') '     FIRST RANDOM ROTATION    '
       ! calculate the centers and spread of the wannier functions
       CALL gzloc_center(c1,nstate,center,modulo)
       IF (paral%io_parent) WRITE(6,&
            ' (/,5x,A,/)')&
            '_____ CENTERS AND SPREAD AFTER FIRST RANDOM ROTATION _____'
       IF (paral%parent) CALL gloc_print(nstate,center,0)
    ENDIF

    ALLOCATE(u_grad(nstate, nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__) ! TODO CHECK STAT

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
    ! - strict condition provided by Lagrange multipliers
    ! - not so strict condition if ROTMAT = 1 + iA 
    ! where A is the self adjoint matrix built in terms of U_GRAD
    ! Aij =0.5(R(U_GRADij) + i I(U_GRADij) - R(U_GRADji) + i I(U_GRADji)) 


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
            __LINE__,__FILE__)! TODO CHECK STAT
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
    ! DO I = 1,NSTATE
    ! write(6,'(I3,8(2f8.4))') I,(ROTMAT(I,J),J=1,8)
    ! ENDDO

    ! ..   calculation of 3 matrices :  <n |O_i| m>
    ! where n and m address the electronic orbitals and i = x,y,z

    ALLOCATE(ddmat(nstate,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL znum_mat(c1,c2,nstate,ddmat,1)
    CALL dcopy(2*nstate*nstate,ddmat,1,z_mat(1,1,1),1)

    CALL znum_mat(c1,c2,nstate,ddmat,2)
    CALL dcopy(2*nstate*nstate,ddmat,1,z_mat(1,1,2),1)

    CALL znum_mat(c1,c2,nstate,ddmat,3)
    CALL dcopy(2*nstate*nstate,ddmat,1,z_mat(1,1,3),1)
    DEALLOCATE(ddmat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)

    CALL mp_sync(parai%allgrp)

    IF (paral%parent) THEN
       ! NOTE: The type of the actual argument differs from the type of
       ! the dummy argument.  [Z_MAT]
       CALL gzfunc(z_mat,nstate,nstate)
       IF (paral%io_parent)&
            WRITE(6,'(/,A,F12.6)') 'STARTING  FUNCTIONAL VALUE = ',gloc_re%ofun
    ENDIF

    DO irep = 1,maxrep

       CALL zeroing(u_grad)!,SIZE(u_grad))

       IF (glocal%tg_antisymm .AND. irep .LE. 100)  THEN
          CALL simm_op(c1,omega_n,gloc_re%omega_tot,nstate)
          IF (glocal%tg_penalty)&
               CALL omega_grad(c1,u_grad,nstate)

          CALL mp_sum(u_grad,nstate*nstate,parai%allgrp)
          CALL mp_sum(omega_n,nstate,parai%allgrp)
          CALL mp_sum(gloc_re%omega_tot,parai%allgrp)

          CALL mp_sync(parai%allgrp)
       ENDIF

       IF (paral%parent) THEN

          ! calculate the Functional  gradient
          CALL u_grad_exp_ide(rotmat,u_grad,z_mat,nstate)

          ! update  of ROTMAT

          CALL u_upd_exp_ide(rotmat,z_mat,xyzmat,u_grad,&
               nstate)

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
          IF (paral%io_parent) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(/,A,1PE12.6)') '     OFUN-OFUN0 = ',&
                  gloc_re%dif_fun
             IF (paral%io_parent)&
                  WRITE(6,'(A,/)') '   G_LOC_UPDATE: CONVERGENCE'
          ENDIF
          GOTO 400
       ELSEIF (irep .EQ. maxrep .OR. soft_com%exsoft) THEN
          GOTO 500
       ENDIF                ! GMAX > EPS
    ENDDO                    ! IREP

500 CONTINUE
    IF (paral%io_parent .AND. soft_com%exsoft) WRITE(6,'(A,1E12.6)')&
         ' G_LOC_EXP_SUM|SOFTWERE EXIT (OFUNC)',gloc_re%ofun
    IF (paral%io_parent .AND. irep .EQ. maxrep) WRITE(6,'(A,I10)')&
         ' G_LOC_EXP_SUM| MAXREP REACHED ',maxrep
400 CONTINUE
    IF (paral%io_parent) THEN
       CALL fileclose(10)
       CALL fileclose(11)
       IF (glocal%tg_antisymm) CALL fileclose(12)
    ENDIF

    ! calculate the centers and spread of the wannier functions
    CALL gzloc_center(c1,nstate,center,modulo)
    IF (paral%io_parent) CALL gloc_print(nstate,center,1)
    DEALLOCATE(center,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(modulo,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)

    IF (paral%io_parent) CALL fileclose(14)

    ! SAVE THE ROTATION MATRIX
    IF (paral%parent) CALL w_matrix(rotmat,nstate)

    DEALLOCATE(lag_mul,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (glocal%tg_antisymm) DEALLOCATE(omega_n,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(z_mat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(xyzmat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rotmat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    IF (cntl%tlsd) DEALLOCATE(rotlsd,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    CALL tihalt('G_LOC_EX_I',isub)
    ! ==--------------------------------------------------------------==
    DEALLOCATE(gmat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(u_grad,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE g_loc_exp_ide
  ! ==================================================================

  SUBROUTINE g_loc_xyzmat(c0,c2,xyzmat,nstate,tau0)
    ! ==--------------------------------------------------------------==
    ! == CALCULATES FUNCTIONAL FOR WFN LOCALIZATION IN DIRECT SPACE   ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(nkpt%ngwk,*), c2(*), &
                                                xyzmat(:,:,:)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: tau0(:,:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'g_loc_xyzmat'

    COMPLEX(real_8), ALLOCATABLE             :: ddmat(:,:), sc1(:)
    INTEGER                                  :: i, ierr, isub, n1, nmcol, nxx
    INTEGER, ALLOCATABLE, SAVE               :: mapcol(:), mapful(:)
    REAL(real_8), ALLOCATABLE                :: center(:,:)

    CALL tiset(' G_LOC_XYZ' ,isub)
    ! ..initialization
    n1=0
    DO i=0,parai%nproc-1
       n1=MAX(n1,parap%sparm(3,i))
    ENDDO
    ngwmax=n1
    ALLOCATE(mapful(2*spar%ngws),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    nmcol=parai%nproc*ngwmax
    ALLOCATE(mapcol(nmcol),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL setdip(mapful,mapcol)

    lenbk=nxxfun(nstate)
    nxx=MAX(2*2*lenbk*parai%nproc,2*nkpt%ngwk*nstate)
    ALLOCATE(sc1(nxx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)


    ! ..initialization for Wannier function part
    ! ..electronic contribution
    ALLOCATE(ddmat(nstate,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL g_loc_opeigr(c0,c2,sc1,nstate,mapful,mapcol,ddmat,&
         1,0)
    CALL dcopy(2*nstate*nstate,ddmat,1,xyzmat(1,1,1),1)
    CALL g_loc_opeigr(c0,c2,sc1,nstate,mapful,mapcol,ddmat,&
         2,0)
    CALL dcopy(2*nstate*nstate,ddmat,1,xyzmat(1,1,2),1)
    CALL g_loc_opeigr(c0,c2,sc1,nstate,mapful,mapcol,ddmat,&
         3,0)
    CALL dcopy(2*nstate*nstate,ddmat,1,xyzmat(1,1,3),1)
    DO i=4,wannc%nwanopt
       CALL g_loc_opeigr(c0,c2,sc1,nstate,mapful,mapcol,ddmat,&
            wannc%iow(1,i-3),wannc%iow(2,i-3))
       CALL dcopy(2*nstate*nstate,ddmat,1,xyzmat(1,1,i),1)
    ENDDO
    DEALLOCATE(ddmat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)

    ! First check of the spread, before the minimization of the functional
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,'(/,1X,64(1H*))')
       IF (paral%io_parent)&
            WRITE(6,'(" ****",6X,A,6X,"****",/)')&
            'CENTERS AND SPREAD BEFORE THE OPTIMIZATION'
       wanni%w_type = 1
       ! W_REF(1) = -0.0_real_8
       ! W_REF(2) = -0.0_real_8
       ! W_REF(3) = -0.0_real_8
    ENDIF

    ALLOCATE(center(4,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL wannier_center(xyzmat,nstate,nstate,center,tau0)
    IF (paral%parent) CALL wc_print(nstate,center)
    DEALLOCATE(center,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)

    DEALLOCATE(sc1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(mapful,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(mapcol,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    CALL tihalt(' G_LOC_XYZ' ,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE g_loc_xyzmat
  ! ==================================================================
  SUBROUTINE give_scr_xyzmat(lddipo,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lddipo
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: lopeigr, lsetdip, nstate

    nstate=crge%n
    ! NMAX=NSTATE
    ! IF(cntl%tlsd) NMAX=MAX(NSUP,NSDOWN)
    CALL give_scr_opeigr(lopeigr,tag,nstate)
    lopeigr = lopeigr + 2*nstate*nstate
    lsetdip=spar%ngws/2+1
    lddipo=MAX(lsetdip,lopeigr)+100
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_xyzmat
  ! ==================================================================
  SUBROUTINE calc_u_wgrad_sum(u_mat,u_grad,xyzmat,nstate)
    ! ==--------------------------------------------------------------==
    ! ==          CALCULATES THE GRADIENT OF THE  FUNCTIONAL          ==
    ! ==          WHICH DEFINES THE SPREAD IN  R SPACE                ==
    ! ==  WITH RESPECT TO THE UNITARY MATRIX COEFFICIENTS             ==
    ! ==  U_GRAD_ij = -d OFUNC / d U_ij                               ==
    ! ==--------------------------------------------------------------==
    ! input
    COMPLEX(real_8)                          :: xyzmat(:,:,:)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: u_grad(nstate,nstate), &
                                                u_mat(nstate,nstate)

    CHARACTER(*), PARAMETER :: procedureN = 'calc_u_wgrad_sum'

    COMPLEX(real_8), ALLOCATABLE             :: uz(:,:), xyzold(:,:,:), &
                                                zu(:,:)
    INTEGER                                  :: i, ierr, isub, ix1, ix2, j, &
                                                k, l, l_uz, l_zu, lmat

!(nstate,nstate,nwanop)
! (NSTATE,NSTATE,WANNC%NWANOPT)
! ==--------------------------------------------------------------==

    CALL tiset(' WGRAD_SUM',isub)
    ALLOCATE(xyzold(nstate,nstate,wannc%nwanopt),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    lmat = 2*nstate*nstate

    CALL dcopy(lmat*wannc%nwanopt,xyzmat,1,xyzold,1)
    CALL z_update(u_mat,xyzmat,nstate,nstate,wannc%nwanopt)

    IF (glocr%wan_weight .EQ. 0.0_real_8) THEN
       GOTO 100
    ELSE

       ix1  = 1
       l_zu = 2*nstate*nstate
       ix2  = 1+l_zu
       l_uz = 2*nstate*nstate

       ALLOCATE(zu(nstate, nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(uz(nstate, nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)

       DO k = 1,wannc%nwanopt
          CALL zeroing(zu)!,SIZE(zu))
          CALL zeroing(uz)!,SIZE(uz))

          DO j = 1,nstate
             DO i = 1,nstate
                DO l = 1,nstate
                   zu(i,j) = zu(i,j) + xyzold(i,l,k)*u_mat(l,j)
                ENDDO
             ENDDO

             DO i = 1,nstate
                DO l = 1,nstate
                   uz(i,j) = uz(i,j) + CONJG(xyzold(l,i,k))*u_mat(l,j)
                ENDDO
             ENDDO

             DO i = 1,nstate
                u_grad(i,j) = u_grad(i,j)+glocr%wan_weight*(zu(i,j)*&
                     CONJG(xyzmat(j,j,k)) +uz(i,j)*xyzmat(j,j,k))
             ENDDO
          ENDDO
       ENDDO
    ENDIF
100 CONTINUE

    CALL gxyzfunc(xyzmat,nstate)
    CALL dcopy(lmat*wannc%nwanopt,xyzold,1,xyzmat,1)
    DEALLOCATE(xyzold,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    CALL tihalt(' WGRAD_SUM',isub)
    ! ==--------------------------------------------------------------==
    DEALLOCATE(zu,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(uz,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE calc_u_wgrad_sum

END MODULE g_loc_optim_utils


SUBROUTINE approx_u(u_mat,u_grad,nstate)
  ! ==--------------------------------------------------------------==
  ! ==      U APPROXIMATED BY 1 +iA WHERE A IS HERMITIAN
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE parac, ONLY : paral,parai
  IMPLICIT NONE
  INTEGER                                    :: nstate
  REAL(real_8)                               :: u_grad(2,nstate,nstate), &
                                                u_mat(2,nstate,nstate)

  INTEGER                                    :: i, j
  REAL(real_8)                               :: fac, im, re

  fac = 0.5_real_8

  DO i = 1,nstate
     DO j = i+1,nstate
        re = fac*(u_grad(1,i,j)-u_grad(1,j,i))
        im = fac*(u_grad(2,i,j)+u_grad(2,j,i))
        u_mat(1,i,j) =  re
        u_mat(2,i,j) =  im
        u_mat(1,j,i) = -re
        u_mat(2,j,i) =  im
     ENDDO
     u_mat(1,i,i) = 1.0_real_8 + u_mat(1,i,i)
  ENDDO
  ! =--------------------------------------------------------------==
  RETURN
END SUBROUTINE approx_u
! ==================================================================



! ==================================================================
SUBROUTINE gzfunc(z_mat,ldx,nstate)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE parac, ONLY : paral,parai
  USE g_loc , ONLY:gloc_re
  IMPLICIT NONE
  INTEGER                                    :: ldx
  REAL(real_8)                               :: z_mat(2,ldx,ldx,3)
  INTEGER                                    :: nstate

  INTEGER                                    :: i, k

  gloc_re%ofun=0._real_8

  ! ..Vanderbilt functional |M|^2
  DO i=1,nstate
     DO k=1,3
        gloc_re%ofun=gloc_re%ofun+(1.0_real_8-&
             z_mat(1,i,i,k)*z_mat(1,i,i,k)-&
             z_mat(2,i,i,k)*z_mat(2,i,i,k))
     ENDDO
  ENDDO

  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE gzfunc
! ==================================================================



! ==================================================================
SUBROUTINE gxyzfunc(xyzmat,nstate)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE system , ONLY:parm
  USE parac, ONLY : paral,parai
  USE g_loc , ONLY:gloc_re
  USE wann , ONLY:wannc
  IMPLICIT NONE
  INTEGER                                    :: nstate
  REAL(real_8) :: xyzmat(2,nstate,nstate,wannc%nwanopt)

  INTEGER                                    :: i, k

  gloc_re%xyzfun = 0._real_8

  DO i=1,nstate
     DO k=1,wannc%nwanopt
        gloc_re%xyzfun=gloc_re%xyzfun+&
             wannc%wwei(k)*(1.0_real_8-(xyzmat(1,i,i,k)*xyzmat(1,i,i,k)+&
             xyzmat(2,i,i,k)*xyzmat(2,i,i,k)))
     ENDDO
  ENDDO
  gloc_re%xyzfun=gloc_re%xyzfun/parm%alat**2
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE gxyzfunc
