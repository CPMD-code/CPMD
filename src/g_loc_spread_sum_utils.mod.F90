MODULE g_loc_spread_sum_utils
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_def
  USE g_loc,                           ONLY: &
       filcen, filpen, filspr, gloc_re, glocal, gloci, glocr, lag_mul, nstep, &
       omega_n
  USE g_loc_spread_ide_utils,          ONLY: g1_mat,&
                                             g2_mat
  USE g_loc_util_utils,                ONLY: gloc_print,&
                                             gzloc_center,&
                                             r_matrix,&
                                             spreadfunc,&
                                             unitary_trans,&
                                             w_matrix,&
                                             zexponentiate
  USE geq0mod,                         ONLY: geq0
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum,&
                                             mp_sync
  USE parac,                           ONLY: parai,&
                                             paral
  USE prng_utils,                      ONLY: repprngu
  USE soft,                            ONLY: soft_com
  USE system,                          ONLY: cntl,&
                                             ncpw,&
                                             nkpt
  USE testex_utils,                    ONLY: testex
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE u_upd_spread_sum_utils,          ONLY: u_grad_spread_sum,&
                                             u_upd_spread_sum
  USE utils,                           ONLY: zclean_k
  USE zeroing_utils,                   ONLY: zeroing
#include "sizeof.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: g_loc_spread_sum

CONTAINS

  ! ==================================================================
  SUBROUTINE g_loc_spread_sum(c1,c2,nstate,nfi)
    ! ==--------------------------------------------------------------==
    ! == CALCULATES LOCALIZED FUNCTIONS IN G SPACE                    ==
    ! ==--------------------------------------------------------------==
    ! INPUT
    COMPLEX(real_8)                          :: c1(nkpt%ngwk,*)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c2(nkpt%ngwk,nstate)
    INTEGER                                  :: nfi

    CHARACTER(*), PARAMETER :: procedureN = 'g_loc_spread_sum'

    COMPLEX(real_8)                          :: zone, zzero
    COMPLEX(real_8), ALLOCATABLE             :: aux(:,:), gmat(:,:), &
                                                u_grad(:,:)
    COMPLEX(real_8), ALLOCATABLE, SAVE       :: rotlsd(:), rotmat(:,:), &
                                                z1_mat(:,:,:), z2_mat(:,:,:)
    INTEGER                                  :: i, ierr, irep, isub, j, &
                                                l_lag_mul, maxrep
    LOGICAL                                  :: ferror
    REAL(real_8), ALLOCATABLE                :: center(:,:), MODULO(:)

! ==--------------------------------------------------------------==
! ..   Localized orbitals functions
! ==--------------------------------------------------------------==

    l_lag_mul  = 2*nstate*nstate
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


    CALL tiset('G_LOC_SP_S',isub)

    gloc_re%gmax = 0.0_real_8
    maxrep = gloci%gloc_maxs
    nstep  = gloci%gloc_init
    zone   = CMPLX(1.0_real_8,0.0_real_8,kind=real_8)
    zzero  = CMPLX(0.0_real_8,0.0_real_8,kind=real_8)

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


    ! CALL PLANE_WAVES_SHELL(C1,NSTATE)

    IF (paral%io_parent)&
         CALL fileopen(14,filcen,fo_def,ferror)
    ! AK 2005/04/08; FIXME, this looks suspicious. why only for serial runs?

    IF (glocal%tg_read_matrix)  THEN
       IF (paral%parent) CALL r_matrix(rotmat,nstate)
       CALL mp_bcast(rotmat,SIZE(rotmat),parai%source,parai%allgrp)
       ! rotate the wavefunctions to the localized representation

       CALL rotate_c((1._real_8,0.0_real_8),c1,(0._real_8,0.0_real_8),c2,rotmat,nstate)
       CALL dcopy(2*2*ncpw%ngw*nstate,c2,1,c1,1)
       IF (geq0)   CALL zclean_k(c1,nstate,ncpw%ngw)
    ENDIF

    ALLOCATE(center(6,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(MODULO(ncpw%ngw),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! calculate the centers and spread of the wannier functions
    CALL gzloc_center(c1,nstate,center,modulo)
    IF (paral%io_parent) WRITE(6,' (/,5x,A,/)')&
         '_______CENTERS AND SPREAD BEFORE THE OPTIMIZATION_______'
    IF (paral%io_parent) CALL gloc_print(nstate,center,0)


    ! ..   calculation of 3 matrices :  <n |O_i| m>
    ! where n and m address the electronic orbitals and i = x,y,z


    IF (glocal%tg_kick) THEN

       CALL unitary_trans(c1,c2,nstate)
       IF (paral%io_parent)&
            WRITE(6,'(/,A)') '     FIRST RANDOM ROTATION    '
       ! calculate the centers and spread of the wannier functions
       CALL gzloc_center(c1,nstate,center,modulo)
       IF (paral%io_parent)&
            WRITE(6,' (/,5x,A,/)')&
            '_____ CENTERS AND SPREAD AFTER FIRST RANDOM ROTATION _____'
       IF (paral%io_parent) CALL gloc_print(nstate,center,0)
    ENDIF

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

    CALL mp_sync(parai%allgrp)


    IF (paral%parent) THEN
       CALL spreadfunc(z1_mat,z2_mat,nstate,nstate)
       IF (paral%io_parent)&
            WRITE(6,'(A,F12.6)') 'AFTER FIRST ROTATION FUNC.= ',gloc_re%ofun
    ENDIF

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

    ! calculate the starting gradient
    ALLOCATE(u_grad(nstate, nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)


    CALL zeroing(u_grad)!,nstate*nstate)
    IF (paral%parent) THEN
       CALL u_grad_spread_sum(rotmat,u_grad,z1_mat,z2_mat,&
            nstate,nstate)
    ENDIF

    ! startoptimization of ROTMAT

    DO irep = 1,maxrep

       ! 
       IF (glocal%tg_antisymm)  THEN
          CALL simm_op_sum(c1,rotmat,omega_n,gloc_re%omega_tot,nstate)
          ! 
          ! &    CALL OMEGA_GRAD(C1,ROTMAT,U_GRAD,NSTATE)


          CALL mp_sum(u_grad,nstate*nstate,parai%allgrp)
          CALL mp_sum(omega_n,nstate,parai%allgrp)
          CALL mp_sum(gloc_re%omega_tot,parai%allgrp)

          CALL mp_sync(parai%allgrp)
       ENDIF

       IF (paral%parent) THEN

          CALL u_upd_spread_sum(rotmat,z1_mat,z2_mat,u_grad,&
               nstate)


       ENDIF
       CALL mp_sync(parai%allgrp)
       CALL mp_bcast_byte(gloc_re,size_in_bytes_of(gloc_re),parai%source,parai%allgrp)

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
       ENDIF
       IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)
       IF (soft_com%exsoft) GOTO 500
    ENDDO

500 IF (soft_com%exsoft.AND.paral%io_parent)&
         WRITE(6,'(A,1E12.6)')&
         ' G_LOC_SPREAD_SUM|SOFTWERE EXIT (OFUNC)',gloc_re%ofun
    IF ((irep .EQ. maxrep).AND.paral%io_parent)&
         WRITE(6,'(A,I10)')&
         ' G_LOC_SPREAD_SUM| MAXREP REACHED ',maxrep
400 CONTINUE
    ! stop '400'
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            CALL fileclose(10)
       IF (paral%io_parent)&
            CALL fileclose(11)
       IF ((glocal%tg_antisymm).AND.paral%io_parent)&
            CALL fileclose(12)
    ENDIF


    CALL mp_bcast(rotmat,SIZE(rotmat),parai%source,parai%allgrp)
    ! rotate the wavefunctions to the localized representation
    CALL rotate_c((1._real_8,0.0_real_8),c1,(0._real_8,0.0_real_8),c2,rotmat,nstate)
    CALL dcopy(2*2*ncpw%ngw*nstate,c2,1,c1,1)
    IF (geq0)   CALL zclean_k(c1,nstate,ncpw%ngw)


    ! calculate the centers and spread of the wannier functions
    CALL gzloc_center(c1,nstate,center,modulo)
    IF (paral%parent) CALL gloc_print(nstate,center,1)
    IF ((paral%parent).AND.paral%io_parent)&
         CALL fileclose(14)

    IF (glocal%tg_read_matrix)  THEN
       IF (paral%parent) CALL r_matrix(lag_mul,nstate)

       CALL mp_bcast(lag_mul,SIZE(lag_mul),parai%source,parai%allgrp)
       CALL zgemm('N','N',nstate,nstate,nstate,zone,lag_mul,nstate,&
            rotmat,nstate,zzero,c2,nstate)
       CALL dcopy(2*nstate*nstate,c2,1,rotmat,1)
    ENDIF

    DEALLOCATE(center,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(modulo,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! SAVE THE ROTATION MATRIX
    IF (paral%parent) CALL w_matrix(rotmat,nstate)

    IF (glocal%tg_antisymm) DEALLOCATE(omega_n,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(lag_mul,STAT=ierr)
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

    ! ==--------------------------------------------------------------==
    IF (paral%parent) THEN
       DEALLOCATE(gmat,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    ENDIF
    DEALLOCATE(u_grad,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL tihalt('G_LOC_SP_S',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE g_loc_spread_sum
  ! ==================================================================
  ! ==================================================================


END MODULE g_loc_spread_sum_utils
