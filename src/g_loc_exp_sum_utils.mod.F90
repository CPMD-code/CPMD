MODULE g_loc_exp_sum_utils
  USE cnst,                            ONLY: fbohr
  USE cppt,                            ONLY: gk
  USE ddipo_utils,                     ONLY: set_operator
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_def
  USE g_loc,                           ONLY: &
       filcen, filpen, filspr, g2g_mem, gloc_re, glocal, gloci, glocr, &
       lag_mul, nstep, omega_n, wan_mem
  USE g_loc_optim_utils,               ONLY: calc_u_wgrad_sum,&
                                             g_loc_xyzmat
  USE g_loc_util_utils,                ONLY: gloc_print,&
                                             gzloc_center,&
                                             r_matrix,&
                                             unitary_trans,&
                                             w_matrix,&
                                             z_update,&
                                             zexponentiate
  USE geq0mod,                         ONLY: geq0
  USE kinds,                           ONLY: real_8
  USE localize_utils,                  ONLY: wc_print
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_recv,&
                                             mp_send,&
                                             mp_sum,&
                                             mp_sync
  USE parac,                           ONLY: parai,&
                                             paral
  USE prng_utils,                      ONLY: repprngu
  USE soft,                            ONLY: soft_com
  USE system,                          ONLY: cntl,&
                                             mapgp,&
                                             ncpw,&
                                             nkpt,&
                                             parap,&
                                             parm,&
                                             spar
  USE testex_utils,                    ONLY: testex
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE u_upd_exp_sum_utils,             ONLY: calc_u_grad_sum,&
                                             u_upd_exp_sum
  USE utils,                           ONLY: zclean_k
  USE wann,                            ONLY: wannc,&
                                             wanni
  USE wannier_center_utils,            ONLY: wannier_center
  USE zeroing_utils,                   ONLY: zeroing
  USE znum_mat_utils,                  ONLY: znum_mat
#include "sizeof.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: g_loc_exp_sum


CONTAINS

  ! ==================================================================
  SUBROUTINE g_loc_exp_sum(c1,c2,nstate,nfi,tau0)
    ! ==--------------------------------------------------------------==
    ! == CALCULATES LOCALIZED FUNCTIONS IN G SPACE                    ==
    ! ==--------------------------------------------------------------==

    ! input
    COMPLEX(real_8)                          :: c1(nkpt%ngwk,*), &
                                                c2(nkpt%ngwk,*)
    INTEGER                                  :: nstate, nfi
    REAL(real_8)                             :: tau0(:,:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'g_loc_exp_sum'

    COMPLEX(real_8)                          :: zone, zzero
    COMPLEX(real_8), ALLOCATABLE             :: aux(:,:), ddmat(:,:), &
                                                gmat(:,:), u_grad(:,:)
    COMPLEX(real_8), ALLOCATABLE, SAVE       :: rotlsd(:), rotmat(:,:), &
                                                xyzmat(:,:,:), z_mat(:,:,:)
    INTEGER                                  :: i, ierr, ig, ip, ipw, irep, &
                                                isub, j, l_lag_mul, len, &
                                                len1, maxrep, msgid
    INTEGER, ALLOCATABLE                     :: mapw(:)
    LOGICAL                                  :: ferror
    REAL(real_8)                             :: fac, grmax, grmax_r
    REAL(real_8), ALLOCATABLE                :: center(:,:), gk_glo(:,:), &
                                                gk_scr(:,:), MODULO(:)
    REAL(real_8), EXTERNAL                   :: ddot

! ==--------------------------------------------------------------==
! ..   Localized orbitals functions
! ==--------------------------------------------------------------==

    CALL tiset('G_LOC_EX_S',isub)

    len  = 3*spar%ngws
    len1 = 2*ncpw%ngw + 1


    IF (glocal%tg_antisymm)  THEN
       ALLOCATE(omega_n(nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(omega_n)!,nstate)
    ENDIF

    l_lag_mul = nstate*nstate
    ALLOCATE(lag_mul(l_lag_mul),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(lag_mul)!,l_lag_mul)



    gloc_re%gmax   = 0.0_real_8
    maxrep = gloci%gloc_maxs
    nstep  = gloci%gloc_init
    zone   = CMPLX(1.0_real_8,0.0_real_8,kind=real_8)
    zzero  = CMPLX(0.0_real_8,0.0_real_8,kind=real_8)

    ! ..   initialization for localized function part

    ALLOCATE(z_mat(nstate,nstate,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(rotmat(nstate,nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    IF (paral%io_parent) CALL fileopen(14,filcen,fo_def,ferror)

    ALLOCATE(gk_glo(3,spar%ngws),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(mapw(len1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    ALLOCATE(gk_scr(3, ncpw%ngw),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__) ! TODO CHECK STAT
    DO ip=0,parai%nproc-1
       IF (paral%parent) THEN
          IF (parap%pgroup(ip+1).EQ.parai%me) THEN
             DO ig=1,ncpw%ngw
                gk_scr(1,ig)=gk(1,ig)
                gk_scr(2,ig)=gk(2,ig)
                gk_scr(3,ig)=gk(3,ig)
                mapw(ig)=mapgp(ig)
             ENDDO
          ELSE
             msgid=1
             !msglen = 3 * parap%sparm(3,ip) * 8
             CALL mp_recv(gk_scr,3 * parap%sparm(3,ip),parap%pgroup(ip+1),msgid,parai%allgrp)
             msgid=2
             !msglen = parap%sparm(3,ip) * 8/irat
             CALL mp_recv(mapw,parap%sparm(3,ip),parap%pgroup(ip+1),msgid,parai%allgrp)
          ENDIF
          DO ipw=1,parap%sparm(3,ip)
             gk_glo(1,mapw(ipw))=gk_scr(1,ipw)
             gk_glo(2,mapw(ipw))=gk_scr(2,ipw)
             gk_glo(3,mapw(ipw))=gk_scr(3,ipw)
          ENDDO
       ELSE
          IF (parap%pgroup(ip+1).EQ.parai%me) THEN
             msgid=1
             !msglen =  3 * ngw * 8
             CALL mp_send(gk,3 * ncpw%ngw,parap%pgroup(1),msgid,parai%allgrp)
             msgid=2
             !msglen = ngw * 8/irat
             CALL mp_send(mapgp,ncpw%ngw,parap%pgroup(1),msgid,parai%allgrp)
          ENDIF
       ENDIF
    ENDDO
    DEALLOCATE(gk_scr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    IF (paral%parent) THEN
       fac = parm%tpiba * fbohr
       DO ig = 1,spar%ngws/2
          IF (paral%io_parent)&
               WRITE(14,'(A,3f12.5,I10)')&
               'X',gk_glo(1,ig)*fac,gk_glo(2,ig)*fac,&
               gk_glo(3,ig)*fac,ig
          IF (paral%io_parent)&
               WRITE(14,'(A,3f12.5,I10)')&
               'X',-gk_glo(1,ig)*fac,-gk_glo(2,ig)*fac,&
               -gk_glo(3,ig)*fac,ig+spar%ngws
       ENDDO
    ENDIF
    DEALLOCATE(gk_glo,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(mapw,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    IF (glocal%tg_read_matrix)  THEN
       IF (paral%parent) CALL r_matrix(rotmat,nstate)
       CALL mp_bcast(rotmat,SIZE(rotmat),parai%source,parai%allgrp)
       ! rotate the wavefunctions to the localized representation

       CALL rotate_c(zone,c1,zzero,c2,rotmat,nstate)
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
    IF (paral%io_parent)&
         WRITE(6,&
         ' (/,5x,A,/)')&
         '_______CENTERS AND SPREAD BEFORE THE OPTIMIZATION_______'
    IF (paral%parent) CALL gloc_print(nstate,center,0)

    IF (glocal%tg_kick) THEN

       CALL unitary_trans(c1,c2,nstate)
       IF (paral%io_parent)&
            WRITE(6,'(/,A)') '     FIRST RANDOM ROTATION    '
       ! calculate the centers and spread of the wannier functions
       CALL gzloc_center(c1,nstate,center,modulo)
       IF (paral%io_parent)&
            WRITE(6,&
            ' (/,5x,A,/)')&
            '_______CENTERS AND SPREAD AFTER ROTATION_______'
       IF (paral%parent) CALL gloc_print(nstate,center,0)

    ENDIF

    ! ..   calculation of 3 matrices :  <n |O_i| m>
    ! where n and m address the electronic orbitals and i = x,y,z

    IF (paral%parent) THEN
       IF (paral%io_parent)&
            CALL fileopen(10,filspr,fo_def,ferror)
       IF (paral%io_parent)&
            CALL fileopen(11,'umat.dat',fo_def,ferror)
       IF (glocal%tg_antisymm.AND.paral%io_parent)&
            CALL fileopen(12,filpen,fo_def,ferror)
       ! initialize U = ROTMAT matrix
       CALL zeroing(rotmat)!,nstate*nstate)
       DO i=1,nstate
          rotmat(i,i)= zone
       ENDDO
       ALLOCATE(gmat(nstate, nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)! TODO CHECK STAT
       CALL zeroing(gmat)!,nstate*nstate)
       IF (glocr%gloc_ran.GT.0._real_8) THEN
          ! ..   randomly distort initial vectors
          DO i=1,nstate
             gmat(i,i)= zone
             DO j=i+1,nstate
                gmat(i,j)=CMPLX(repprngu()*glocr%gloc_ran,repprngu()*glocr%gloc_ran,kind=real_8)
                gmat(j,i)=CONJG(gmat(i,j))
             ENDDO
          ENDDO
          CALL zexponentiate(gmat,rotmat,nstate)
          glocr%gloc_ran=0._real_8
       ENDIF
    ENDIF

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

    IF (glocal%tgwannier) THEN
       CALL set_operator(.FALSE.)
       ALLOCATE(xyzmat(nstate,nstate,wannc%nwanopt),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL g_loc_xyzmat(c1,c2,xyzmat,nstate,tau0)
       IF (paral%parent)  THEN
          CALL gxyzfunc(xyzmat,nstate)
          IF (paral%io_parent)&
               WRITE(6,'(/,A,F12.6)')&
               'STARTING  REALSPCE FUNCTIONAL VALUE =  ',gloc_re%xyzfun
       ENDIF
    ENDIF

    IF (paral%parent) THEN
       CALL gzfunc(z_mat,nstate,nstate)
       IF (paral%io_parent)&
            WRITE(6,'(/,A,F12.6)')&
            'STARTING  G FUNCTIONAL VALUE =  ',gloc_re%ofun
    ENDIF


    ! startoptimization of ROTMAT
    ALLOCATE(u_grad(nstate, nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    CALL zeroing(u_grad)!,nstate*nstate)

    ! CALCULATE HE MASSES  RATIO
    gloc_re%mass_r = 1.0_real_8
    gloc_re%mass_g = 1.0_real_8


    IF (glocal%tgwannier .AND. glocr%wan_weight .NE. 0.0_real_8&
         .AND. glocr%g2g_weight .NE. 0.0_real_8) THEN

       glocr%wan_weight = 1.0_real_8
       glocr%g2g_weight = 1.0_real_8

       IF (paral%parent) THEN
          CALL calc_u_wgrad_sum(rotmat,u_grad,xyzmat,nstate)

          ALLOCATE(aux(nstate,nstate),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
          CALL zgemm('C','N',nstate,nstate,nstate,zone,u_grad,&
               nstate,rotmat,nstate,zzero,aux,nstate)
          CALL zgemm('N','N',nstate,nstate,nstate,-zone,rotmat,nstate,&
               aux,nstate,zone,u_grad,nstate)
          DEALLOCATE(aux,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
               __LINE__,__FILE__)
          ! IGRMAX  = IZAMAX(NSTATE*NSTATE,U_GRAD,1)
          ! GRMAX_R = Cabs(ZGIVE(U_GRAD,IGRMAX))
          grmax_r = 0.0_real_8
          DO i = 1,nstate
             grmax_r = grmax_r+ddot(2*nstate,u_grad(1,i),1,u_grad(1,i),1)
          ENDDO
          CALL zeroing(u_grad)!,nstate*nstate)
          CALL calc_u_grad_sum(rotmat,u_grad,z_mat,nstate,nstate)
          ! IGRMAX = IZAMAX(NSTATE*NSTATE,U_GRAD,1)
          ! GRMAX  = Cabs(ZGIVE(U_GRAD,IGRMAX))
          grmax = 0.0_real_8
          DO i = 1,nstate
             grmax = grmax+ddot(2*nstate,u_grad(1,i),1,u_grad(1,i),1)
          ENDDO

          gloc_re%mass_r = SQRT(grmax_r)/REAL(nstate,kind=real_8)
          gloc_re%mass_g = SQRT(grmax)/REAL(nstate,kind=real_8)
          IF (paral%io_parent)&
               WRITE(6,'(/,A,2(1PE12.4))') 'masses(gnorm)',gloc_re%mass_g,gloc_re%mass_r
          fac = MAX(gloc_re%mass_r,gloc_re%mass_g)
          gloc_re%mass_r = gloc_re%mass_r/fac
          gloc_re%mass_g = gloc_re%mass_g/fac
       ENDIF
       glocr%wan_weight = wan_mem
       glocr%g2g_weight = g2g_mem
    ENDIF
    glocr%g2g_weight = glocr%g2g_weight/gloc_re%mass_g
    glocr%wan_weight = glocr%wan_weight/gloc_re%mass_r
    IF (paral%io_parent)&
         WRITE(6,'(A,2(1PE12.4))') 'real weight',glocr%g2g_weight,glocr%wan_weight
    IF (glocr%g2g_weight .GT. 1.e3_real_8) glocr%g2g_weight=glocr%g2g_weight / 10._real_8
    IF (glocr%wan_weight .GT. 1.e3_real_8) glocr%wan_weight=glocr%wan_weight / 10._real_8
    IF (paral%io_parent)&
         WRITE(6,'(A,2(1PE12.4),/)')&
         'final weight: ',glocr%g2g_weight,glocr%wan_weight
    gloc_re%mass_r = 1.0_real_8
    gloc_re%mass_g = 1.0_real_8

    IF (glocal%tgwannier) THEN
       IF (paral%parent) THEN
          CALL calc_u_wgrad_sum(rotmat,u_grad,xyzmat,nstate)
       ENDIF
    ENDIF

    IF (paral%parent) THEN
       CALL calc_u_grad_sum(rotmat,u_grad,z_mat,nstate,nstate)
    ENDIF


    DO irep=1,maxrep
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
          ! calculate U_GRAD_new with the rotated U_MAT_new

          CALL u_upd_exp_sum(rotmat,z_mat,xyzmat,u_grad,&
               nstate)

       ENDIF


       CALL mp_sync(parai%allgrp)

       CALL mp_bcast_byte(gloc_re,size_in_bytes_of(gloc_re),parai%source,parai%allgrp)

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
       ENDIF
    ENDDO

500 CONTINUE
    IF ((paral%parent .AND. soft_com%exsoft).AND.paral%io_parent)&
         WRITE(6,'(A,1E12.6)')&
         ' G_LOC_EXP_SUM|SOFTWERE EXIT (OFUNC)',gloc_re%ofun
    IF ((paral%parent .AND. irep .EQ. maxrep).AND.paral%io_parent)&
         WRITE(6,'(A,I10)')&
         ' G_LOC_EXP_SUM| MAXREP REACHED ',maxrep
400 CONTINUE
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

    CALL rotate_c(zone,c1,zzero,c2,rotmat,nstate)
    CALL dcopy(2*2*ncpw%ngw*nstate,c2,1,c1,1)
    IF (geq0)   CALL zclean_k(c1,nstate,ncpw%ngw)

    ! Final value of the spread
    IF (glocal%tgwannier) THEN
       IF (paral%parent) THEN
          CALL z_update(rotmat,xyzmat,nstate,nstate,wannc%nwanopt)
          IF (paral%io_parent)&
               WRITE(6,'(/,1X,64(1H*))')
          IF (paral%io_parent)&
               WRITE(6,'(" ****",6X,A,6X,"****",/)')&
               'CENTERS AND SPREAD AFTER THE OPTIMIZATION'
          wanni%w_type = 1
       ENDIF
       CALL wannier_center(xyzmat,nstate,nstate,center,tau0)
       IF (paral%parent) CALL wc_print(nstate,center)
    ENDIF

    ! calculate the centers and spread of the wannier functions
    CALL gzloc_center(c1,nstate,center,modulo)
    IF (paral%parent) CALL gloc_print(nstate,center,1)

    DEALLOCATE(center,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(modulo,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    IF ((paral%parent).AND.paral%io_parent)&
         CALL fileclose(14)

    IF (glocal%tg_read_matrix)  THEN
       IF (paral%parent) CALL r_matrix(lag_mul,nstate)

       CALL mp_bcast(lag_mul,SIZE(lag_mul),parai%source,parai%allgrp)
       CALL zgemm('N','N',nstate,nstate,nstate,zone,lag_mul,&
            nstate,rotmat,nstate,zzero,c2,nstate)
       CALL dcopy(2*nstate*nstate,c2,1,rotmat,1)
    ENDIF

    ! SAVE THE ROTATION MATRIX
    IF (paral%parent) CALL w_matrix(rotmat,nstate)

    IF (glocal%tgwannier) DEALLOCATE(xyzmat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (glocal%tg_antisymm) DEALLOCATE(omega_n,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    DEALLOCATE(lag_mul,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(z_mat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rotmat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    IF (cntl%tlsd) DEALLOCATE(rotlsd,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    DEALLOCATE(u_grad,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(gmat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    CALL tihalt('G_LOC_EX_S',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE g_loc_exp_sum
  ! ==================================================================
  ! ==================================================================

END MODULE g_loc_exp_sum_utils
