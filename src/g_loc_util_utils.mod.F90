MODULE g_loc_util_utils
  USE cnst,                            ONLY: fbohr
  USE cppt,                            ONLY: gk,&
                                             hg
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_def,&
                                             fo_old,&
                                             fo_ufo
  USE g_loc,                           ONLY: filmat,&
                                             gloc_re,&
                                             glocal,&
                                             gloci,&
                                             glocr,&
                                             lag_mul
  USE geq0mod,                         ONLY: geq0
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum,&
                                             mp_sync
  USE parac,                           ONLY: parai,&
                                             paral
  USE prng_utils,                      ONLY: repprngu
  USE rgs_utils,                       ONLY: uinvc
  USE system,                          ONLY: ncpw,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: zclean_k,&
                                             zgive
  USE wann,                            ONLY: wannc
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: line_search_exe
  PUBLIC :: line_search_s




  PUBLIC :: z_update

  PUBLIC :: spreadfunc
  PUBLIC :: z12_update

  PUBLIC :: u_by_ortho
  PUBLIC :: s12ortho

  PUBLIC :: lagrmult
  !public :: g_xgradat0
  !public :: g_ofunc
  !public :: g2_update
  PUBLIC :: zexponentiate
  PUBLIC :: gzloc_center
  !public :: plane_waves_shell
  PUBLIC :: unitary_trans
  PUBLIC :: gloc_print
  PUBLIC :: w_matrix
  PUBLIC :: r_matrix

CONTAINS



  ! 
  ! ==================================================================
  SUBROUTINE line_search_exe(u_mat,u_grad,z_mat,xyzmat,a_mat,b_mat,&
       nstate,delta,factor)
    ! ==--------------------------------------------------------------==
    ! ==   ORTHOGONALIZE THE MATRIX U IN ORDER TO KEEP THE UNITARITY  ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: xyzmat(*)
    INTEGER                                  :: nstate
    COMPLEX(real_8) :: b_mat(nstate,nstate), a_mat(nstate,nstate), &
      z_mat(nstate,nstate,3), u_grad(nstate,nstate), u_mat(nstate,nstate)
    REAL(real_8)                             :: delta, factor

    CHARACTER(*), PARAMETER                  :: procedureN = 'line_search_exe'

    COMPLEX(real_8)                          :: zone, zzero
    COMPLEX(real_8), ALLOCATABLE             :: c_mat(:,:), d_mat(:,:), &
                                                Uis_MAT(:,:), zold(:,:,:)
    INTEGER                                  :: i, ierr, isub, j, length
    REAL(real_8)                             :: afac, dx, ofun1, ofun2, &
                                                tot_fun, x0

! ==--------------------------------------------------------------==

    CALL tiset(' LINE_SRCH',isub)

    length = 3*2*nstate+2*2*nstate+nstate+2*nstate*nstate
    ALLOCATE(Uis_MAT(nstate, nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(c_mat(nstate, nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(d_mat(nstate, nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(zold(nstate, nstate, 3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    zone  = CMPLX(1.0_real_8,0.0_real_8,kind=real_8)
    zzero = CMPLX(0.0_real_8,0.0_real_8,kind=real_8)

    ! INITIALIZE VARIABLE
    gloc_re%ofun0  = glocr%g2g_weight*gloc_re%ofun+ glocr%wan_weight*gloc_re%xyzfun
    afac   = 1.0_real_8
    ! write(6,*) 1,G2G_WEIGHT,OFUN, WAN_WEIGHT,XYZFUN
    ! ....  FIRST POINT


    ! TEMPORARY U
    IF (gloci%gloc_opt .EQ. 1) THEN
       DO i = 1,nstate
          DO j = 1,nstate
             Uis_MAT(i,j)  = u_mat(i,j) + delta*u_grad(i,j)
          ENDDO
       ENDDO
    ELSE
       DO i = 1,nstate
          DO j = 1,nstate
             Uis_MAT(i,j)  =  delta*u_grad(i,j)
          ENDDO
          Uis_MAT(i,i)  =  uis_mat(i,i) + 1.0_real_8
       ENDDO
    ENDIF

    ! ORTHOGONALIZATION OF TEMPORARY U
    IF (gloci%gloc_const .EQ. 2) THEN
       CALL dcopy(2*nstate*nstate,Uis_MAT,1,a_mat,1)
       CALL ovlap2_c(nstate,nstate,nstate,b_mat,a_mat,Uis_MAT)
       CALL dcopy(2*nstate*nstate,b_mat,1,a_mat,1)
       CALL s12ortho(a_mat,b_mat,nstate)
       CALL zgemm('N','N',nstate,nstate,nstate,zone,Uis_MAT,&
            nstate,b_mat,nstate,zzero,a_mat,nstate)
       CALL dcopy(2*nstate*nstate,a_mat,1,Uis_MAT,1)

    ELSEIF (gloci%gloc_const .EQ. 1) THEN

       CALL zeroing(c_mat)!,SIZE(c_mat))
       CALL zeroing(d_mat)!,SIZE(d_mat))
       CALL  lagrmult(u_mat,Uis_MAT,a_mat,b_mat,c_mat,d_mat,&
            nstate)

    ENDIF


    ! STORE Z_MAT AND UPDATE Z_MAT AND FUNCTIONAL
    IF (glocr%g2g_weight .EQ. 0.0_real_8) GOTO 100
    CALL dcopy(3*2*nstate*nstate,z_mat,1,zold,1)
    CALL z_update(Uis_MAT,z_mat,nstate,nstate,3)
    CALL gzfunc(z_mat,nstate,nstate)
    CALL dcopy(3*2*nstate*nstate,zold,1,z_mat,1)
100 CONTINUE
    IF (glocal%tgwannier .AND. glocr%wan_weight .NE. 0.0_real_8) THEN
       CALL dcopy(wannc%nwanopt*2*nstate*nstate,xyzmat,1,zold,1)
       CALL z_update(Uis_MAT,xyzmat,nstate,nstate,wannc%nwanopt)
       CALL gxyzfunc(xyzmat,nstate)
       CALL dcopy(wannc%nwanopt*2*nstate*nstate,zold,1,xyzmat,1)
    ENDIF

    tot_fun = glocr%g2g_weight*gloc_re%ofun + glocr%wan_weight*gloc_re%xyzfun
    ! COMPARE OLD AND NEW FUNCTIONAL
    IF (tot_fun.GT.gloc_re%ofun0) THEN
       ofun2= tot_fun
       afac=0.5_real_8
    ELSE
       ofun1= tot_fun
       afac=2.0_real_8
    ENDIF


    ! ....  SECOND POINT

    ! TEMPORARY U
    IF (gloci%gloc_opt .EQ. 1) THEN
       DO i = 1,nstate
          DO j = 1,nstate
             Uis_MAT(i,j)  = u_mat(i,j) + afac * delta*u_grad(i,j)
          ENDDO
       ENDDO
    ELSE
       DO i = 1,nstate
          DO j = 1,nstate
             Uis_MAT(i,j)  =  afac * delta*u_grad(i,j)
          ENDDO
          Uis_MAT(i,i)  =  uis_mat(i,i) + 1.0_real_8
       ENDDO
    ENDIF
    ! ORTHOGONALIZATION OF TEMPORARY U
    IF (gloci%gloc_const .EQ. 2) THEN
       CALL dcopy(2*nstate*nstate,Uis_MAT,1,a_mat,1)
       CALL ovlap2_c(nstate,nstate,nstate,b_mat,a_mat,Uis_MAT)
       CALL dcopy(2*nstate*nstate,b_mat,1,a_mat,1)
       CALL s12ortho(a_mat,b_mat,nstate)
       CALL zgemm('N','N',nstate,nstate,nstate,zone,Uis_MAT,&
            nstate,b_mat,nstate,zzero,a_mat,nstate)
       CALL dcopy(2*nstate*nstate,a_mat,1,Uis_MAT,1)

    ELSEIF (gloci%gloc_const .EQ. 1) THEN

       CALL zeroing(c_mat)!, SIZE(c_mat))
       CALL zeroing(d_mat)!, SIZE(d_mat))
       CALL  lagrmult(u_mat,Uis_MAT,a_mat,b_mat,c_mat,d_mat,&
            nstate)

    ENDIF

    ! STORE Z_MAT AND UPDATE Z_MAT AND FUNCTIONAL
    IF (glocr%g2g_weight .EQ. 0.0_real_8)  GOTO 200
    CALL dcopy(3*2*nstate*nstate,z_mat,1,zold,1)
    CALL z_update(Uis_MAT,z_mat,nstate,nstate,3)
    CALL gzfunc(z_mat,nstate,nstate)
    CALL dcopy(3*2*nstate*nstate,zold,1,z_mat,1)
200 CONTINUE
    IF (glocal%tgwannier .AND. glocr%wan_weight .NE. 0.0_real_8) THEN
       CALL dcopy(wannc%nwanopt*2*nstate*nstate,xyzmat,1,zold,1)
       CALL z_update(Uis_MAT,xyzmat,nstate,nstate,wannc%nwanopt)
       CALL gxyzfunc(xyzmat,nstate)
       CALL dcopy(wannc%nwanopt*2*nstate*nstate,zold,1,xyzmat,1)
    ENDIF

    tot_fun = glocr%g2g_weight*gloc_re%ofun + glocr%wan_weight*gloc_re%xyzfun
    IF (afac.GT.1._real_8) THEN
       ofun2=tot_fun
       dx=1.0_real_8
    ELSE
       ofun1=tot_fun
       dx=0.5_real_8
    ENDIF

    ! Extrapolation
    IF ((2.0_real_8*gloc_re%ofun0-4.0_real_8*ofun1+2.0_real_8*ofun2).EQ.0._real_8) THEN
       x0=dx
    ELSE
       x0=dx*(4.0_real_8*ofun1-ofun2-3.0_real_8*gloc_re%ofun0)/&
            (2.0_real_8*gloc_re%ofun0-4.0_real_8*ofun1+2.0_real_8*ofun2)
    ENDIF
    ! write(6,'(I5,6(2x,1PE12.4))') NSTEP,DX,OFUN0,OFUN1,OFUN2,X0,DELTA
    ! ALPHA=(OFUN1-OFUN0)/DX/(2*X0+DX)
    ! EE=OFUN0-ALPHA*X0*X0
    afac=-x0
    IF (afac.GT.0._real_8) THEN
       IF (afac.GT.3._real_8) afac=3.0_real_8
    ELSE
       IF (ofun2.LT.gloc_re%ofun0) THEN
          afac=2*dx
       ELSEIF (ofun1.LT.gloc_re%ofun0) THEN
          afac=dx
       ELSE
          afac=0.05_real_8
       ENDIF
    ENDIF
    ! DE=EE+ALPHA*(X0+AFAC)**2
    ! ALPHA=ALPHA/DD

    ! .....FINAL EVALUATION


    ! ORTHOGONALIZATION OF TEMPORARY U
    IF (gloci%gloc_const .EQ. 2) THEN

       ! TEMPORARY U
       IF (gloci%gloc_opt .EQ. 1) THEN
          DO i = 1,nstate
             DO j = 1,nstate
                u_mat(i,j)  = u_mat(i,j) + afac * delta*u_grad(i,j)
             ENDDO
          ENDDO
       ELSE
          DO i = 1,nstate
             DO j = 1,nstate
                u_mat(i,j)  =  afac * delta*u_grad(i,j)
             ENDDO
             u_mat(i,i)  =  u_mat(i,i) + 1.0_real_8
          ENDDO
       ENDIF

       CALL dcopy(2*nstate*nstate,u_mat,1,a_mat,1)
       CALL ovlap2_c(nstate,nstate,nstate,b_mat,a_mat,u_mat)
       CALL dcopy(2*nstate*nstate,b_mat,1,a_mat,1)
       CALL s12ortho(a_mat,b_mat,nstate)
       CALL zgemm('N','N',nstate,nstate,nstate,zone,u_mat,&
            nstate,b_mat,nstate,zzero,a_mat,nstate)
       CALL dcopy(2*nstate*nstate,a_mat,1,u_mat,1)

    ELSEIF (gloci%gloc_const .EQ. 1) THEN

       ! TEMPORARY U
       IF (gloci%gloc_opt .EQ. 1) THEN
          DO i = 1,nstate
             DO j = 1,nstate
                Uis_MAT(i,j)  = u_mat(i,j) + afac * delta*u_grad(i,j)
             ENDDO
          ENDDO
       ELSE
          DO i = 1,nstate
             DO j = 1,nstate
                Uis_MAT(i,j)  =  afac * delta*u_grad(i,j)
             ENDDO
             Uis_MAT(i,i)  =  uis_mat(i,i) + 1.0_real_8
          ENDDO
       ENDIF

       CALL  lagrmult(u_mat,Uis_MAT,a_mat,b_mat,c_mat,d_mat,&
            nstate)

       CALL dcopy(2*nstate*nstate,Uis_MAT,1,u_mat,1)
    ENDIF


    IF (afac.LT.0.1_real_8) THEN
       factor=factor*0.25_real_8
    ELSEIF (afac.GT.1.00_real_8) THEN
       factor=factor*1.50_real_8
    ELSEIF (afac.LT.0.20_real_8) THEN
       factor=factor*0.75_real_8
    ENDIF
    ! write(6,'(I8,7(1PE12.4))') NSTEP, ofun0,ofun1-ofun0,ofun2-ofun0,
    ! &                          ofun2-ofun1,afac,x0,factor

    ! ==--------------------------------------------------------------==
    DEALLOCATE(Uis_MAT,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(c_mat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(d_mat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(zold,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL tihalt(' LINE_SRCH',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE line_search_exe
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE line_search_s(u_mat,u_grad,z1_mat,z2_mat,a_mat,&
       b_mat,nstate,delta,factor)
    ! ==--------------------------------------------------------------==
    ! ==   ORTHOGONALIZE THE MATRIX U IN ORDER TO KEEP THE UNITARITY  ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    COMPLEX(real_8) :: b_mat(nstate,nstate), a_mat(nstate,nstate), &
      z2_mat(nstate,nstate,3), z1_mat(nstate,nstate,3), &
      u_grad(nstate,nstate), u_mat(nstate,nstate)
    REAL(real_8)                             :: delta, factor

    CHARACTER(*), PARAMETER                  :: procedureN = 'line_search_s'

    COMPLEX(real_8)                          :: zone, zzero
    COMPLEX(real_8), ALLOCATABLE             :: Uis_MAT(:,:), zold(:,:)
    INTEGER                                  :: i, ierr, ip1, j, length
    REAL(real_8)                             :: afac, dx, ofun1, ofun2, x0

! ==--------------------------------------------------------------==

    zone  = CMPLX(1.0_real_8,0.0_real_8,kind=real_8)
    zzero = CMPLX(0.0_real_8,0.0_real_8,kind=real_8)

    ! INITIALIZE VARIABLE
    gloc_re%ofun0  = gloc_re%ofun
    afac   = 1.0_real_8

    ! ....  FIRST POINT
    ALLOCATE(Uis_MAT(nstate, nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)


    ! TEMPORARY U
    IF (gloci%gloc_opt .EQ. 1) THEN
       DO i = 1,nstate
          DO j = 1,nstate
             Uis_MAT(i,j)  = u_mat(i,j) + delta*u_grad(i,j)
          ENDDO
       ENDDO
    ELSE
       DO i = 1,nstate
          DO j = 1,nstate
             Uis_MAT(i,j)  =  delta*u_grad(i,j)
          ENDDO
          Uis_MAT(i,i)  =  uis_mat(i,i) + 1.0_real_8
       ENDDO
    ENDIF
    ! ORTHOGONALIZATION OF TEMPORARY U
    CALL dcopy(2*nstate*nstate,Uis_MAT,1,a_mat,1)
    CALL ovlap2_c(nstate,nstate,nstate,b_mat,a_mat,Uis_MAT)
    CALL dcopy(2*nstate*nstate,b_mat,1,a_mat,1)

    CALL s12ortho(a_mat,b_mat,nstate)

    CALL zgemm('N','N',nstate,nstate,nstate,zone,Uis_MAT,&
         nstate,b_mat,nstate,zzero,a_mat,nstate)
    CALL dcopy(2*nstate*nstate,a_mat,1,Uis_MAT,1)


    ! STORE Z_MAT AND UPDATE Z_MAT AND FUNCTIONAL
    ALLOCATE(zold(2,nstate*nstate*3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL dcopy(3*2*nstate*nstate,z1_mat,1,zold(1,1),1)
    CALL dcopy(3*2*nstate*nstate,z2_mat,1,zold(2,1),1)
    CALL z12_update(Uis_MAT,z1_mat,z2_mat,nstate,nstate)
    CALL spreadfunc(z1_mat,z2_mat,nstate,nstate)
    CALL dcopy(3*2*nstate*nstate,zold(1,1),1,z1_mat,1)
    CALL dcopy(3*2*nstate*nstate,zold(2,1),1,z2_mat,1)


    ! COMPARE OLD AND NEW FUNCTIONAL
    IF (gloc_re%ofun.GT.gloc_re%ofun0) THEN
       ofun2=gloc_re%ofun
       afac=0.5_real_8
    ELSE
       ofun1=gloc_re%ofun
       afac=2.0_real_8
    ENDIF


    ! ....  SECOND POINT

    ! TEMPORARY U
    IF (gloci%gloc_opt .EQ. 1) THEN
       DO i = 1,nstate
          DO j = 1,nstate
             Uis_MAT(i,j)  = u_mat(i,j) + afac * delta*u_grad(i,j)
          ENDDO
       ENDDO
    ELSE
       DO i = 1,nstate
          DO j = 1,nstate
             Uis_MAT(i,j)  =  afac * delta*u_grad(i,j)
          ENDDO
          Uis_MAT(i,i)  =  uis_mat(i,i) + 1.0_real_8
       ENDDO
    ENDIF
    ! ORTHOGONALIZATION OF TEMPORARY U
    CALL dcopy(2*nstate*nstate,Uis_MAT,1,a_mat,1)
    CALL ovlap2_c(nstate,nstate,nstate,b_mat,a_mat,Uis_MAT)
    CALL dcopy(2*nstate*nstate,b_mat,1,a_mat,1)
    CALL s12ortho(a_mat,b_mat,nstate)
    CALL zgemm('N','N',nstate,nstate,nstate,zone,Uis_MAT,&
         nstate,b_mat,nstate,zzero,a_mat,nstate)
    CALL dcopy(2*nstate*nstate,a_mat,1,Uis_MAT,1)

    ! STORE Z_MAT AND UPDATE Z_MAT AND FUNCTIONAL
    length = 2*3*2*nstate*nstate+2*2*nstate*nstate
    ip1     = 1   + 3*nstate*nstate  ! 3*2*NSTATE*NSTATE
    CALL dcopy(3*2*nstate*nstate,z1_mat,1,zold(1,1),1)
    CALL dcopy(3*2*nstate*nstate,z2_mat,1,zold(2,1),1)
    CALL z12_update(Uis_MAT,z1_mat,z2_mat,nstate,nstate)
    CALL spreadfunc(z1_mat,z2_mat,nstate,nstate)
    CALL dcopy(3*2*nstate*nstate,zold(1,1),1,z1_mat,1)
    CALL dcopy(3*2*nstate*nstate,zold(2,1),1,z2_mat,1)

    IF (afac.GT.1._real_8) THEN
       ofun2=gloc_re%ofun
       dx=1.0_real_8
    ELSE
       ofun1=gloc_re%ofun
       dx=0.5_real_8
    ENDIF
    ! Extrapolation
    IF ((2*gloc_re%ofun0-4*ofun1+2*ofun2).EQ.0._real_8) THEN
       x0=dx
    ELSE
       x0=dx*(4*ofun1-ofun2-3*gloc_re%ofun0)/(2*gloc_re%ofun0-4*ofun1+2*ofun2)
    ENDIF
    ! ALPHA=(OFUN1-OFUN0)/DX/(2*X0+DX)
    ! EE=OFUN0-ALPHA*X0*X0
    afac=-x0
    IF (afac.GT.0._real_8) THEN
       IF (afac.GT.3._real_8) afac=3.0_real_8
    ELSE
       IF (ofun2.LT.gloc_re%ofun0) THEN
          afac=2*dx
       ELSEIF (ofun1.LT.gloc_re%ofun0) THEN
          afac=dx
       ELSE
          afac=0.05_real_8
       ENDIF
    ENDIF
    ! DE=EE+ALPHA*(X0+AFAC)**2
    ! ALPHA=ALPHA/DD

    ! .....FINAL EVALUATION

    ! TEMPORARY U
    IF (gloci%gloc_opt .EQ. 1) THEN
       DO i = 1,nstate
          DO j = 1,nstate
             u_mat(i,j)  = u_mat(i,j) + afac * delta*u_grad(i,j)
          ENDDO
       ENDDO
    ELSE
       DO i = 1,nstate
          DO j = 1,nstate
             u_mat(i,j)  =  afac * delta*u_grad(i,j)
          ENDDO
          u_mat(i,i)  =  u_mat(i,i) + 1.0_real_8
       ENDDO
    ENDIF
    ! ORTHOGONALIZATION OF TEMPORARY U
    CALL dcopy(2*nstate*nstate,u_mat,1,a_mat,1)
    CALL ovlap2_c(nstate,nstate,nstate,b_mat,a_mat,u_mat)
    CALL dcopy(2*nstate*nstate,b_mat,1,a_mat,1)
    CALL s12ortho(a_mat,b_mat,nstate)
    CALL zgemm('N','N',nstate,nstate,nstate,zone,u_mat,&
         nstate,b_mat,nstate,zzero,a_mat,nstate)
    CALL dcopy(2*nstate*nstate,a_mat,1,u_mat,1)


    IF (afac.LT.0.1_real_8) THEN
       factor=factor*0.25_real_8
    ELSEIF (afac.GT.1.0_real_8) THEN
       factor=factor*1.5_real_8
    ELSEIF (afac.LT.0.2_real_8) THEN
       factor=factor*0.75_real_8
    ENDIF
    ! write(6,*)  ofun0,ofun1,ofun2,factor
    ! write(6,'(I8,7(1PE12.4))') NSTEP, ofun0,ofun1-ofun0,ofun2-ofun0,
    ! &                          ofun2-ofun1,afac,x0,factor

    DEALLOCATE(zold,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE line_search_s
  ! ==================================================================



  ! ==--------------------------------------------------------------==
  SUBROUTINE z_update(rmat,z_mat,ldx,nstate,nmat)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ldx, nstate
    COMPLEX(real_8)                          :: rmat(nstate,nstate)
    INTEGER                                  :: nmat
    COMPLEX(real_8)                          :: z_mat(ldx,ldx,nmat)

    CHARACTER(*), PARAMETER                  :: procedureN = 'z_update'

    COMPLEX(real_8), ALLOCATABLE             :: aux(:,:,:)
    INTEGER                                  :: i, ierr, isub

    CALL tiset('  Z_UPDATE',isub)
    ALLOCATE(aux(nstate,nstate,2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    DO i=1,nmat
       CALL dcopy(2*nstate*nstate,z_mat(1,1,i),1,aux(1,1,1),1)
       CALL zgemm('C','N',nstate,nstate,nstate,&
            (1._real_8,0.0_real_8),rmat,nstate,&
            aux(1,1,1),nstate,(0._real_8,0.0_real_8),aux(1,1,2),nstate)
       CALL zgemm('N','N',nstate,nstate,nstate,&
            (1._real_8,0.0_real_8),aux(1,1,2),nstate,&
            rmat,nstate,(0._real_8,0.0_real_8),aux(1,1,1),nstate)
       CALL dcopy(2*nstate*nstate,aux(1,1,1),1,z_mat(1,1,i),1)
    ENDDO
    DEALLOCATE(aux,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    CALL tihalt('  Z_UPDATE',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE z_update
  ! ==================================================================


  ! ==================================================================
  SUBROUTINE spreadfunc(z1_mat,z2_mat,ldx,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ldx
    COMPLEX(real_8)                          :: z2_mat(ldx,ldx,3), &
                                                z1_mat(ldx,ldx,3)
    INTEGER                                  :: nstate

    INTEGER                                  :: i, k

    gloc_re%ofun=0._real_8

    DO i=1,nstate
       DO k=1,3
          gloc_re%ofun=gloc_re%ofun+REAL(z2_mat(i,i,k)-z1_mat(i,i,k)*z1_mat(i,i,k))
       ENDDO
    ENDDO

    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE spreadfunc
  ! ==================================================================


  SUBROUTINE z12_update(rmat,z1_mat,z2_mat,ldx,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ldx
    COMPLEX(real_8)                          :: z2_mat(ldx,ldx,3), &
                                                z1_mat(ldx,ldx,3)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: rmat(nstate,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'z12_update'

    COMPLEX(real_8), ALLOCATABLE             :: aux(:,:,:)
    INTEGER                                  :: i, ierr, j, k

! Variables
! (NSTATE,NSTATE,2)
! ==--------------------------------------------------------------==

    ALLOCATE(aux(nstate,nstate,2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    DO i=1,3

       ! CALL DCOPY(2*NSTATE*NSTATE,Z1_MAT(1,1,I),1,AUX(1,1,1),1)
       DO j=1,nstate
          DO k=1,nstate
             aux(j,k,1)=z1_mat(j,k,i)
          ENDDO
       ENDDO
       CALL zgemm('C','N',nstate,nstate,nstate,(1._real_8,0.0_real_8),rmat,&
            nstate,aux(1,1,1),nstate,(0._real_8,0.0_real_8),aux(1,1,2),nstate)
       CALL zgemm('N','N',nstate,nstate,nstate,(1._real_8,0.0_real_8),aux(1,1,2),&
            nstate,rmat,nstate,(0._real_8,0.0_real_8),aux(1,1,1),nstate)
       DO j=1,nstate
          DO k=1,nstate
             z1_mat(j,k,i)=aux(j,k,1)
          ENDDO
       ENDDO

       DO j=1,nstate
          DO k=1,nstate
             aux(j,k,1)=z2_mat(j,k,i)
          ENDDO
       ENDDO

       CALL zgemm('C','N',nstate,nstate,nstate,(1._real_8,0.0_real_8),rmat,&
            nstate,aux(1,1,1),nstate,(0._real_8,0.0_real_8),aux(1,1,2),nstate)
       CALL zgemm('N','N',nstate,nstate,nstate,(1._real_8,0.0_real_8),aux(1,1,2),&
            nstate,rmat,nstate,(0._real_8,0.0_real_8),aux(1,1,1),nstate)
       DO j=1,nstate
          DO k=1,nstate
             z2_mat(j,k,i)=aux(j,k,1)
          ENDDO
       ENDDO
       ! CALL DCOPY(2*NSTATE*NSTATE,AUX(1,1,1),1,Z2_MAT(1,1,I),1)
    ENDDO
    DEALLOCATE(aux,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE z12_update
  ! ==================================================================


  ! ==================================================================
  SUBROUTINE u_by_ortho(u_mat,u_grad,Us_MAT,ov_mat,&
       nstate,delta)
    ! ==--------------------------------------------------------------==
    ! ==   ORTHOGONALIZE THE MATRIX U IN ORDER TO KEEP THE UNITARITY  ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    COMPLEX(real_8) :: ov_mat(nstate,nstate), Us_MAT(nstate,nstate), &
      u_grad(nstate,nstate), u_mat(nstate,nstate)
    REAL(real_8)                             :: delta

    COMPLEX(real_8), PARAMETER :: zone = (1.0_real_8,0.0_real_8) , &
      zzero = (0.0_real_8,0.0_real_8)

    INTEGER                                  :: i, j, length

! ==--------------------------------------------------------------==

    length = 3*2*nstate+2*2*nstate+nstate

    IF (gloci%gloc_opt .EQ. 1) THEN
       DO i = 1,nstate
          DO j = 1,nstate
             u_mat(i,j)  = u_mat(i,j) + delta*u_grad(i,j)
          ENDDO
       ENDDO
    ELSE
       DO i = 1,nstate
          DO j = 1,nstate
             u_mat(i,j)  =  delta*u_grad(i,j)
          ENDDO
          u_mat(i,i)  =  u_mat(i,i) + 1.0_real_8
       ENDDO
    ENDIF

    CALL dcopy(2*nstate*nstate,u_mat,1,Us_MAT,1)

    CALL ovlap2_c(nstate,nstate,nstate,ov_mat,Us_MAT,u_mat)

    CALL dcopy(2*nstate*nstate,ov_mat,1,Us_MAT,1)

    CALL s12ortho(Us_MAT,ov_mat,nstate)

    CALL zgemm('N','N',nstate,nstate,nstate,zone,u_mat,&
         nstate,ov_mat,nstate,zzero,Us_MAT,nstate)

    CALL dcopy(2*nstate*nstate,Us_MAT,1,u_mat,1)

    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE u_by_ortho
  ! ==================================================================
  SUBROUTINE  s12ortho(Us_MAT,ov_mat,nstate)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: ov_mat(nstate,nstate), &
                                                Us_MAT(nstate,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 's12ortho'

    COMPLEX(real_8), ALLOCATABLE             :: work(:)
    INTEGER                                  :: i, ierr, j, k, l_aux, l_w, &
                                                lwork
    REAL(real_8)                             :: im, im1, re, re1
    REAL(real_8), ALLOCATABLE                :: aux(:), w(:)

    l_w   = nstate
    l_aux = 3*2*nstate
    lwork  = MAX(1,2*2*nstate-1)

    ALLOCATE(w(l_w),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(aux(l_aux),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(work(lwork),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    CALL zheev('V','U',nstate,Us_MAT,nstate,w,work,lwork,aux,ierr)

    DO i=1,nstate
       w(i)=1._real_8/SQRT(w(i))
    ENDDO

    CALL zeroing(ov_mat)!,nstate*nstate)
    DO k=1,nstate
       DO j=1,nstate
          re = w(k)*REAL(Us_MAT(j,k))
          im =-w(k)*AIMAG(Us_MAT(j,k))
          DO i=1,nstate
             re1 = re*REAL(Us_MAT(i,k)) - im*AIMAG(us_mat(i,k))
             im1 = re*AIMAG(Us_MAT(i,k)) + im*REAL(us_mat(i,k))
             ov_mat(i,j) =  ov_mat(i,j) + CMPLX(re1,im1,kind=real_8)
          ENDDO
       ENDDO
    ENDDO

    ! ==--------------------------------------------------------------==
    DEALLOCATE(w,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(aux,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(work,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE s12ortho


  ! ==================================================================
  SUBROUTINE lagrmult(u_mat,u_temp,a,b1,b2,gam,nstate)
    ! ==--------------------------------------------------------------==
    ! ==    CALCULATES THE MATRIX OF LAGRANGE MULTIPLIERS             ==
    ! ==    IN ORDER TO IMPOSE THE UNITARITY CONDITION                ==
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: nstate
    COMPLEX(real_8) :: gam(nstate*nstate), b2(nstate,nstate), &
      b1(nstate,nstate), a(nstate,nstate), u_temp(nstate,nstate), &
      u_mat(nstate,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'lagrmult'

    COMPLEX(real_8)                          :: zone, zzero
    COMPLEX(real_8), ALLOCATABLE             :: scr1(:), scr2(:)
    INTEGER                                  :: i, ierr, ii, imax, isub, &
                                                iter, izamax, j
    REAL(real_8)                             :: difgam

! Variables
! correction due to orthonormality constraint
! iterative scheme for evaluation of Lagrangian multipliers
! LAG_MUL_new = 0.5[1 - A +(1-B^*)LAG_MUL_old + 
! LAG_MUL_old(1-B) + LAG_MUL_old LAG_MUL_old ]
! A = U_MAT_new * U_MAT_new
! B = U_MAT_old * U_MAT_new
! ==-------------------------------------------------------------==

    CALL tiset('  LAGRMULT',isub)


    zzero = CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
    zone  = CMPLX(1.0_real_8,0.0_real_8,kind=real_8)


    ALLOCATE(scr1(nstate*nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(scr2(nstate*nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    CALL zgemm('C','N',nstate,nstate,nstate,zone,u_temp,nstate,&
         u_temp,nstate,zzero,a(1,1),nstate)


    CALL zgemm('C','N',nstate,nstate,nstate,zone,u_mat,nstate,&
         u_temp,nstate,zzero,b1(1,1),nstate)

    DO i=1,nstate
       DO j=1,nstate
          b2(j,i)=CONJG(b1(i,j))
       ENDDO
    ENDDO

    CALL dscal(2*nstate*nstate,-1.0_real_8,a(1,1),1)
    CALL dscal(2*nstate*nstate,-1.0_real_8,b1(1,1),1)
    CALL dscal(2*nstate*nstate,-1.0_real_8,b2(1,1),1)

    DO i=1,nstate
       a(i,i)=CMPLX(1.0_real_8,0.0_real_8,kind=real_8)+a(i,i)
       b1(i,i)=CMPLX(1.0_real_8,0.0_real_8,kind=real_8)+b1(i,i)
       b2(i,i)=CMPLX(1.0_real_8,0.0_real_8,kind=real_8)+b2(i,i)
    ENDDO

    ! ==----------------------------------------------------------==
    ! == INITIALIZE LAG_MUL                                       ==
    ! ==----------------------------------------------------------==
    DO j=1,nstate
       DO i=1,nstate
          ii=i+(j-1)*nstate
          lag_mul(ii)=a(i,j)*0.5_real_8
       ENDDO
    ENDDO

    ! ==----------------------------------------------------------==
    ! == INITIALIZE FIXED TERM OF RECURRENCE                      ==
    ! ==----------------------------------------------------------==
    CALL zgemm('N','N',nstate,nstate,nstate,(0.5_real_8,0.0_real_8),b2(1,1),&
         nstate,b1(1,1),nstate,(0.5_real_8,0.0_real_8),a(1,1),nstate)

    iter=0
    ! ==----------------------------------------------------------==
    ! == ITERATIVE LOOP                                           ==
    ! ==----------------------------------------------------------==

100 CONTINUE
    iter=iter+1
    CALL dcopy(2*nstate*nstate,a(1,1),1,gam(1),1)
    DO j=1,nstate
       DO i=1,nstate
          ii=i+(j-1)*nstate
          scr1(ii)=b1(i,j)-lag_mul(ii)
          scr2(ii)=b2(i,j)-lag_mul(ii)
       ENDDO
    ENDDO

    CALL zgemm('N','N',nstate,nstate,nstate,(-0.5_real_8,0.0_real_8),&
         scr2,nstate,scr1,nstate,zone,gam(1),nstate)

    DO i=1,nstate*nstate
       scr1(i)=gam(i)-lag_mul(i)
    ENDDO
    CALL dcopy(2*nstate*nstate,gam(1),1,lag_mul(1),1)
    imax=izamax(nstate*nstate,scr1,1)
    difgam=ABS(zgive(scr1(1),imax))
    CALL dcopy(2*nstate*nstate,a(1,1),1,gam(1),1)

    ! convergence criterium

    IF (difgam.GT.glocr%gepslag*0.5_real_8.AND.iter.LE.gloci%gloc_maxit) GOTO 100
    IF (iter.GT.gloci%gloc_maxit+20) THEN
       IF (paral%io_parent)&
            WRITE(6,*) 'U_UPD_SPREAD:DIFGAM=',difgam,' ITER = ',&
            iter,' > ',gloci%gloc_maxit
       IF (paral%io_parent)&
            WRITE(6,*) '       MAXIMUM NUMBER OF ITERATIONS EXCEEDED'
       CALL stopgm('U_UPD_SPREAD',' ',& 
            __LINE__,__FILE__)
    ENDIF

    ! ==--------------------------------------------------------------==
    ! == UPDATE  U_MAT_new by  LAG_MU                                 ==
    ! =---------------------------------------------------------------==
    CALL zgemm('N','N',nstate,nstate,nstate,zone,u_mat,nstate,&
         lag_mul,nstate,zone,u_temp,nstate)
    ! ==--------------------------------------------------------------==
    DEALLOCATE(scr1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(scr2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL tihalt('  LAGRMULT',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE lagrmult
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE g_xgradat0(wt,gmat,xyzmat,rotmat,ldx,nstate)
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: wt, ldx
    COMPLEX(real_8)                          :: xyzmat(ldx,ldx,3)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: rotmat(nstate,nstate), &
                                                gmat(nstate,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'g_xgradat0'

    COMPLEX(real_8)                          :: xa, xb
    COMPLEX(real_8), ALLOCATABLE, SAVE       :: zlocal(:,:,:)
    INTEGER                                  :: i, ierr, j, k
    INTEGER, SAVE                            :: ifirst = 0
    REAL(real_8)                             :: gcomp

    CALL zeroing(gmat)!,nstate*nstate)
    IF (wt.EQ.1 ) THEN
       ! ..Spread functional G^2
       DO i=1,nstate
          DO j=i+1,nstate
             gcomp=0._real_8
             DO k=1,3
                xa=xyzmat(j,j,k)-xyzmat(i,i,k)
                gcomp=gcomp+REAL(4._real_8*CONJG(xyzmat(i,j,k))*xa)
             ENDDO
             gmat(i,j)=gcomp
             gmat(j,i)=-gmat(i,j)
          ENDDO
       ENDDO
    ELSEIF (wt.EQ.2) THEN
       IF (ifirst .EQ. 0) THEN
          ALLOCATE(zlocal(ldx,ldx,3),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          CALL zeroing(zlocal)!,ldx*ldx*3)
          CALL dcopy(2*ldx*ldx*3,xyzmat,1,zlocal,1)
          ifirst = 1
       ELSE

          DO i=1,nstate
             DO j=1,nstate
                gcomp=0._real_8
                DO k = 1,3
                   gcomp= gcomp+REAL(zlocal(j,i,k)-xyzmat(j,i,k))
                   ! DO II = 1,NSTATE
                   ! GCOMP= GCOMP + 
                   ! &      2._real_8*ROT(I,II)*
                   ! &     (real(ZLOCAL(II,J,K))*real(XYZMAT(I,I,K))
                   ! &      +aimag(ZLOCAL(II,J,K))*aimag(XYZMAT(I,I,K)))+
                   ! &      2._real_8*ROT(II,J)*
                   ! &     (real(ZLOCAL(I,II,K))*real(XYZMAT(J,J,K))
                   ! &      +aimag(ZLOCAL(I,II,K))*aimag(XYZMAT(J,J,K)))
                   ! ENDDO
                ENDDO
                gmat(j,i)=gcomp
             ENDDO
          ENDDO
          CALL dcopy(2*ldx*ldx*3,xyzmat,1,zlocal,1)

       ENDIF
    ELSEIF (wt.EQ.3) THEN
       ! ..Resta functional Log(|M|^2)
       DO i=1,nstate
          DO j=i+1,nstate
             gcomp=0._real_8
             DO k=1,3
                xa=REAL(xyzmat(i,j,k)/xyzmat(i,i,k))
                xb=REAL(xyzmat(i,j,k)/xyzmat(j,j,k))
                gcomp=gcomp+2._real_8*REAL(xb-xa)
             ENDDO
             gmat(i,j)=gcomp
             gmat(j,i)=-gmat(i,j)
          ENDDO
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE g_xgradat0
  ! ==================================================================
  SUBROUTINE g_ofunc(wt,funcv,xyzmat,ldx,nstate)
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: wt
    REAL(real_8)                             :: funcv
    INTEGER                                  :: ldx
    COMPLEX(real_8)                          :: xyzmat(ldx,ldx,3)
    INTEGER                                  :: nstate

    INTEGER                                  :: i, k

    funcv=0._real_8
    IF (wt.EQ.1 .OR.wt.EQ.2) THEN
       ! ..Spread functional G^2
       DO i=1,nstate
          DO k=1,3
             funcv=funcv+REAL(CONJG(xyzmat(i,i,k))*xyzmat(i,i,k))
          ENDDO
       ENDDO
    ELSEIF (wt.EQ.3) THEN
       ! ..Resta functional Log(|M|^2)
       DO i=1,nstate
          DO k=1,3
             funcv=funcv-LOG(REAL(CONJG(xyzmat(i,i,k))*xyzmat(i,i,k)))
          ENDDO
       ENDDO
    ENDIF
    funcv=funcv/parm%tpiba2
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE g_ofunc
  ! ==================================================================
  SUBROUTINE g2_update(rmat,g2_mat,scr,ldx,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ldx
    REAL(real_8)                             :: g2_mat(2,ldx,ldx,3)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: scr(nstate,nstate,3), &
                                                rmat(nstate,nstate)

    INTEGER                                  :: i, j, k

    DO i=1,3
       DO j=1,nstate
          DO k=1,nstate
             scr(j,k,1)=g2_mat(1,j,k,i)
          ENDDO
       ENDDO
       CALL dgemm('T','N',nstate,nstate,nstate,1._real_8,rmat,nstate,scr(1,&
            1,1),nstate,0._real_8,scr(1,1,2),nstate)
       CALL dgemm('N','N',nstate,nstate,nstate,1._real_8,scr(1,1,2),nstate,&
            rmat,nstate,0._real_8,scr(1,1,1),nstate)
       DO j=1,nstate
          DO k=1,nstate
             g2_mat(1,j,k,i)=scr(j,k,1)
          ENDDO
       ENDDO
       DO j=1,nstate
          DO k=1,nstate
             scr(j,k,1)=g2_mat(2,j,k,i)
          ENDDO
       ENDDO
       CALL dgemm('T','N',nstate,nstate,nstate,1._real_8,rmat,nstate,scr(1,&
            1,1),nstate,0._real_8,scr(1,1,2),nstate)
       CALL dgemm('N','N',nstate,nstate,nstate,1._real_8,scr(1,1,2),nstate,&
            rmat,nstate,0._real_8,scr(1,1,1),nstate)
       DO j=1,nstate
          DO k=1,nstate
             g2_mat(2,j,k,i)=scr(j,k,1)
          ENDDO
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE g2_update
  ! ==================================================================
  SUBROUTINE zexponentiate(gmat,rmat,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: rmat(nstate,nstate), &
                                                gmat(nstate,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'zexponentiate'

    COMPLEX(real_8), ALLOCATABLE             :: eia(:), hgmat(:,:), work(:)
    INTEGER                                  :: i, ierr, info, j, k, lwork
    REAL(real_8)                             :: cei, ei, sei, w(32)
    REAL(real_8), ALLOCATABLE                :: rwork(:)

! variables
! ==--------------------------------------------------------------==
! Scratch array

    lwork = 25*nstate ! TODO optimal LWORK?
    ALLOCATE(work(lwork),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(rwork(3*nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    ALLOCATE(hgmat(nstate, nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(eia(nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! ..calculate eigenvalues of iG (this is Hermite schen)
    DO i=1,nstate
       DO j=i,nstate
          hgmat(i,j)=gmat(i,j)
       ENDDO
    ENDDO
    CALL zheev('V','U',nstate,hgmat,nstate,w,work,lwork,rwork,info)
    DO i=1,nstate
       ei= w(i)
       cei=COS(ei)
       sei=SIN(ei)
       eia(i)=CMPLX(cei,sei,kind=real_8)
    ENDDO
    CALL zeroing(rmat)!,nstate*nstate)
    DO k=1,nstate
       DO i=1,nstate
          DO j=1,nstate
             rmat(i,j)=rmat(i,j)+(eia(k)*hgmat(k,i)*CONJG(hgmat(k,j)))
          ENDDO
       ENDDO
    ENDDO
    ! DO I=1,NSTATE
    ! RMAT(I,I)=(1.0_real_8,0.0_real_8)
    ! DO J=1,NSTATE
    ! RMAT(I,J)=RMAT(I,J)+(0.0_real_8,1.0_real_8)*GMAT(I,J)
    ! DO K = 1,NSTATE
    ! RMAT(I,J)=RMAT(I,J)-0.5_real_8*GMAT(I,K)*GMAT(K,J)
    ! DO L = 1,NSTATE
    ! RMAT(I,J)=RMAT(I,J)-(0.0_real_8,1.0_real_8)*
    ! &               GMAT(I,K)*GMAT(K,L)*GMAT(L,J)/6.0_real_8
    ! ENDDO
    ! ENDDO
    ! ENDDO
    ! ENDDO
    ! ==--------------------------------------------------------------==
    DEALLOCATE(rwork,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(work,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(hgmat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(eia,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE zexponentiate
  ! ==================================================================
  SUBROUTINE gzloc_center(c1,nstate,center,modulo)
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c1(2*ncpw%ngw,nstate)
    REAL(real_8)                             :: center(6,nstate), &
                                                MODULO(ncpw%ngw)

    INTEGER                                  :: i, iG, inx, istate
    REAL(real_8)                             :: c02

    CALL zeroing(center)!,6*nstate)
    DO iG = 1,ncpw%ngw
       MODULO(iG) =  SQRT(hg(ig))
    ENDDO

    inx = 1
    IF (geq0) inx = 2

    DO istate = 1,nstate
       DO iG = 1, ncpw%ngw
          c02 = REAL(c1(iG,istate))*REAL(c1(ig,istate))+&
               AIMAG(c1(iG,istate))*AIMAG(c1(ig,istate))
          center(1,istate) = center(1,istate) + gk(1,iG)*c02
          center(2,istate) = center(2,istate) + gk(2,iG)*c02
          center(3,istate) = center(3,istate) + gk(3,iG)*c02

          center(4,istate) = center(4,istate) + MODULO(iG)*c02
          center(5,istate) = center(5,istate) + hg(iG)*c02
          ! CENTER(5,ISTATE) = CENTER(5,ISTATE) + 
          ! &                  (GK(1,IG)*GK(1,IG)+GK(2,IG)*GK(2,IG))*C02
       ENDDO

       DO iG = inx, ncpw%ngw
          c02 = REAL(c1(iG+ncpw%ngw,istate))*REAL(c1(ig+ncpw%ngw,istate))+&
               AIMAG(c1(iG+ncpw%ngw,istate))*AIMAG(c1(ig+ncpw%ngw,istate))
          center(1,istate) = center(1,istate) - gk(1,iG)*c02
          center(2,istate) = center(2,istate) - gk(2,iG)*c02
          center(3,istate) = center(3,istate) - gk(3,iG)*c02

          center(4,istate) = center(4,istate) + MODULO(iG)*c02
          center(5,istate) = center(5,istate) + hg(iG)*c02
          ! CENTER(5,ISTATE) = CENTER(5,ISTATE) + 
          ! &                  (GK(1,IG)*GK(1,IG)+GK(2,IG)*GK(2,IG))*C02
       ENDDO

       DO i = 1,4
          center(i,istate) = center(i,istate)*parm%tpiba
       ENDDO
       center(5,istate) = center(5,istate)*parm%tpiba2
    ENDDO
    CALL mp_sync(parai%allgrp)
    CALL mp_sum(center,6*nstate,parai%allgrp)

    ! SPREAD <G^2> - <|G|>^2
    DO istate = 1,nstate
       center(4,istate) = SQRT(center(1,istate)*center(1,istate)+&
            center(2,istate)*center(2,istate)+&
            center(3,istate)*center(3, istate))

       center(6,istate) = center(5,istate) -center(4,istate)*center(4,&
            istate)
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gzloc_center
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE  plane_waves_shell(c0,nstate)
    ! ==--------------------------------------------------------------==
    ! == CALCULATES WFN AS COMBINATION OF PLANE WAVES 
    ! == BELONGING TO THE SAME G SHELL
    ! ==--------------------------------------------------------------==
    ! == DONE FOR THE SPECIAL CASE OF 32 STATES IN A CUBIC CELL
    ! == G VECTORS DISTRIBUTION IS HOMOGENEOUS
    ! ==--------------------------------------------------------------==
    ! NOT PARALLEL


    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(2*ncpw%ngw,nstate)

    INTEGER                                  :: iG, igg, istate, ng, numwf(50)
    LOGICAL                                  :: build
    REAL(real_8)                             :: alpha, alphasq, alphasq1, &
                                                eps, hshell, norm

    CALL zeroing(c0)!,nstate*2*ngw)
    CALL zeroing(numwf)!,nstate)
    ng = 1
    eps = 0.00001_real_8
    norm = 0.0_real_8


    build = .FALSE.
    IF (build) THEN
       DO istate = 1,nstate
          hshell = hg(ng)
          DO iG = ng,ncpw%ngw
             IF (hg(iG) .GE. hshell-eps .AND. hg(ig) .LE. hshell+eps )&
                  THEN
                c0(iG,istate) = CMPLX(1.0_real_8,0.0_real_8,kind=real_8)
                c0(iG+ncpw%ngw,istate) = CMPLX(0.0_real_8,1.0_real_8,kind=real_8)
                numwf(istate) = numwf(istate) + 2
             ELSEIF (hg(iG) .GT. hshell+eps) THEN
                ! NORM = SQRT(real(NUMWF(ISTATE),kind=real_8))
                ! write(6,'(/,2I6)') ISTATE,NUMWF(ISTATE)
                ! DO IGG = NG,IG-1
                ! C0(IGG,ISTATE) = C0(IGG,ISTATE)/NORM
                ! C0(IGG+NGW,ISTATE) = C0(IGG+NGW,ISTATE)/NORM
                ! write(6,'(I6,2f12.6)') IGG,HG(IGG),NORM
                ! ENDDO
                IF (paral%io_parent)&
                     WRITE(6,'(/,2I6)') istate,numwf(istate)
                DO igg = ng,iG-1
                   norm = norm + REAL(c0(igg,istate))*REAL(c0(igg,&
                        istate))+AIMAG(c0(igg,istate))*AIMAG(c0(igg,&
                        istate))
                   norm = norm + REAL(c0(igg+ncpw%ngw,istate))*REAL(c0(igg+&
                        ncpw%ngw,istate))+AIMAG(c0(igg+ncpw%ngw,istate))*&
                        AIMAG(c0(igg+ncpw%ngw,istate))
                ENDDO
                norm = SQRT(norm)
                DO igg = ng,iG-1
                   c0(igg,istate) = c0(igg,istate)/norm
                   c0(igg+ncpw%ngw,istate) = c0(igg+ncpw%ngw,istate)/norm
                   IF (paral%io_parent)&
                        WRITE(6,'(I6,2f12.6)') igg,hg(igg),norm
                ENDDO
                norm = 0.0_real_8
                ng = iG
                GOTO 10
             ELSEIF (hg(iG) .LT. hshell-eps) THEN
                IF (paral%io_parent)&
                     WRITE(6,*) ng,iG,hshell,hg(ig)
                CALL stopgm('PLANE_WAVES_SHELL','G vectors order',& 
                     __LINE__,__FILE__)
             ENDIF
          ENDDO
10        CONTINUE
       ENDDO
       c0(1,1) = 2*c0(1,1)
       CALL zclean_k(c0,nstate,ncpw%ngw)
    ENDIF
    build = .TRUE.
    IF (build) THEN
       alpha = repprngu()
       alphasq = SQRT(alpha)
       alphasq1 = SQRT(1.0_real_8-alpha)
       c0(1,1) = CMPLX(alphasq,alphasq1,kind=real_8)
       IF (paral%io_parent)&
            WRITE(6,'(A,I6,A,I8,2f10.5)')&
            'ISTATE = ',1, ' iG = ',1,c0(1,1)
       iG = 2
       DO istate = 2,nstate,2
          IF (istate+1 .LE. nstate) THEN
             alpha = repprngu()
             alphasq = SQRT(alpha)
             alphasq1 = SQRT(1.0_real_8-alpha)
             c0(iG,istate) = CMPLX(alphasq,alphasq1,kind=real_8)
             alpha = repprngu()
             alphasq = SQRT(alpha)
             alphasq1 = SQRT(1.0_real_8-alpha)
             c0(iG+ncpw%ngw,istate+1) = CMPLX(alphasq,alphasq1,kind=real_8)
             IF (paral%io_parent)&
                  WRITE(6,'(A,I6,A,I8,2f10.5)')&
                  'ISTATE = ',istate, ' iG = ',ig,c0(ig,istate)
             IF (paral%io_parent)&
                  WRITE(6,'(A,I6,A,I8,2f10.5)')&
                  'ISTATE+1 = ',istate+1, ' iG = ',ig+ncpw%ngw,c0(ig+ncpw%ngw,istate+1)
             iG = ig+1
          ELSEIF (istate .EQ. nstate) THEN
             alpha = repprngu()
             alphasq = SQRT(alpha)/SQRT(2.0_real_8)
             alphasq1 = SQRT(1.0_real_8-alpha)/SQRT(2.0_real_8)
             c0(iG,istate) = CMPLX(alphasq,alphasq1,kind=real_8)
             c0(iG+ncpw%ngw,istate) = CMPLX(alphasq1,alphasq,kind=real_8)
             IF (paral%io_parent)&
                  WRITE(6,'(A,I6,A,I8,A,I8,4f10.5)')&
                  'ISTATE = ',istate,' iG = ',ig,' and ',ig+ncpw%ngw,&
                  c0(iG,istate),c0(ig+ncpw%ngw,istate)
          ENDIF
       ENDDO
    ENDIF
    RETURN
  END SUBROUTINE plane_waves_shell
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE  unitary_trans(c0,c2,nstate)
    ! ==--------------------------------------------------------------==
    ! == CALCULATES THE MATRIX REPRESENTING A GENERAL UNITARY TRANSFORMATION
    ! == BELONGING TO THE SAME G SHELL
    ! ==--------------------------------------------------------------==
    ! NOT PARALLEL

    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c2(2*ncpw%ngw,nstate), &
                                                c0(2*ncpw%ngw,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'unitary_trans'

    COMPLEX(real_8)                          :: zone
    COMPLEX(real_8), ALLOCATABLE             :: ov_mat(:,:), u_mat(:,:), &
                                                Us_MAT(:,:)
    INTEGER                                  :: ierr, iG, istate, jstate, lmat

    lmat = 2*nstate*nstate
    ALLOCATE(u_mat(nstate, nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(Us_MAT(nstate, nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(ov_mat(nstate, nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL zeroing(u_mat)!,nstate*nstate)
    CALL zeroing(Us_MAT)!,nstate*nstate)
    CALL zeroing(ov_mat)!,nstate*nstate)
    CALL zeroing(c2)!,2*ngw*nstate)
    zone = CMPLX(1.0_real_8,0.0_real_8,kind=real_8)
    IF (paral%parent) THEN
       DO istate = 1,nstate
          DO jstate = 1,nstate
             u_mat(istate,jstate) = CMPLX(repprngu(),repprngu(),kind=real_8)
          ENDDO
       ENDDO

       CALL dcopy(2*nstate*nstate,u_mat,1,Us_MAT,1)

       CALL ovlap2_c(nstate,nstate,nstate,ov_mat,Us_MAT,u_mat)
       CALL zgemm('N','C',nstate,nstate,nstate,-zone,u_mat,nstate,&
            ov_mat,nstate,zone,Us_MAT,nstate)


       CALL zeroing(ov_mat)!,nstate*nstate)
       CALL zherk('U','C',nstate,nstate,1._real_8,Us_MAT,nstate,0._real_8,ov_mat,&
            nstate)
       CALL uinvc('U',ov_mat,nstate,nstate)
       CALL ztrmm('R','U','N','N',nstate,nstate,CMPLX(1._real_8,0._real_8,kind=real_8),&
            ov_mat,nstate,Us_MAT,nstate)
    ENDIF

    CALL mp_bcast(Us_MAT,SIZE(Us_MAT),parai%source,parai%allgrp)

    ! DO ISTATE = 1,NSTATE
    ! write(6,'(I5,2f12.6)') ISTATE,Us_MAT(ISTATE,1)
    ! ENDDO
    ! stop

    DO istate = 1,nstate
       DO jstate = 1,nstate
          DO iG = 1,2*ncpw%ngw
             c2(iG,istate)= c2(ig,istate)+c0(ig,jstate)*Us_MAT(jstate,&
                  istate)
          ENDDO
       ENDDO

    ENDDO

    CALL dcopy(4*ncpw%ngw*nstate,c2,1,c0,1)

    ! ==--------------------------------------------------------------==
    DEALLOCATE(u_mat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(Us_MAT,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(ov_mat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE unitary_trans

  ! ==================================================================
  SUBROUTINE gloc_print(nstate,center,index)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    REAL(real_8)                             :: center(6,nstate)
    INTEGER                                  :: index

    INTEGER                                  :: i
    REAL(real_8)                             :: sum, sum1, sum2, sum3, sum4, &
                                                sum5, sum6

! ==--------------------------------------------------------------==

    IF (paral%io_parent)&
         WRITE(6,'(/,1X,82(1H*))')
    IF (paral%io_parent)&
         WRITE(6,'(" *",16X,A,5X,A,A,5x,A,1X,"*")')&
         ' AVERAGE G VECTOR ','           <|G|>    ',&
         '   <G^2>    ','<G^2> - <|G|>^2'
    IF (paral%io_parent)&
         WRITE(6,'(1X,82(1H*))')
    sum = 0.0_real_8
    sum1 = 0.0_real_8
    sum2 = 0.0_real_8
    sum3 = 0.0_real_8
    sum4 = 0.0_real_8
    sum5 = 0.0_real_8
    sum6 = 0.0_real_8
    DO i=1,nstate
       IF (paral%io_parent)&
            WRITE(6,'(I6,3F12.4,5X,2F12.6,5x,F12.8)') i,&
            center(1,i),center(2,i),center(3,i),center(4,i),&
            center(5,i),center(6,i)
       sum = sum + ABS(center(6,i))
       sum1 = sum1 + center(1,i)*center(1,i)
       sum2 = sum2 + center(2,i)*center(2,i)
       sum3 = sum3 + center(3,i)*center(3,i)
       sum4 = sum4 + center(1,i)
       sum5 = sum5 + center(2,i)
       sum6 = sum6 + center(3,i)
    ENDDO
    IF (paral%io_parent)&
         WRITE(6,'(/,1X,"SUM OVER ALL STATES OF THE SPREAD ",1PE14.6,/)')&
         sum
    IF (paral%io_parent)&
         WRITE(6,'(/,1X,"SUM OVER ALL STATES OF <Gi> ",3(2x,f12.6),/)')&
         sum4,sum5,sum6
    IF (paral%io_parent)&
         WRITE(6,'(/,1X,"SUM OVER ALL STATES OF <Gi>^2",3(2x,1PE12.6),/)')&
         sum1,sum2,sum3
    IF (paral%io_parent)&
         WRITE(6,'(1X,82(1H*),/)')

    IF (index .EQ. 1) THEN
       DO i = 1,nstate
          IF (paral%io_parent)&
               WRITE(14,'(A,3f12.5,I10)') 'C',center(1,i)*fbohr,&
               center(2,i)*fbohr,center(3,i)*fbohr,i
       ENDDO
    ENDIF

    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gloc_print
  ! ==================================================================
  SUBROUTINE w_matrix(rotmat,nstate)

    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: rotmat(nstate*nstate)

    LOGICAL                                  :: ferror

    IF (paral%io_parent)&
         CALL fileopen(10,filmat,fo_def+fo_ufo,ferror)
    IF (paral%io_parent)&
         WRITE(10) nstate
    IF (paral%io_parent)&
         WRITE(10) rotmat

    IF (paral%io_parent)&
         CALL fileclose(10)
    RETURN
  END SUBROUTINE w_matrix
  ! ==================================================================
  SUBROUTINE r_matrix(rotmat,nstate)


    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: rotmat(nstate*nstate)

    INTEGER                                  :: nstate0
    LOGICAL                                  :: ferror

    IF (paral%io_parent)&
         CALL fileopen(10,filmat,fo_old+fo_ufo,ferror)
    IF (ferror) THEN
       IF (paral%io_parent)&
            WRITE(6,'(A,A)') ' READ_MATRIX| FILE NOT FOUND:',filmat
       CALL stopgm('READ_MATRIX',' FILE NOT FOUND ',& 
            __LINE__,__FILE__)
    ENDIF

    IF (paral%io_parent)&
         READ(10) nstate0
    IF (nstate0 .NE. nstate) THEN
       IF (paral%io_parent)&
            WRITE(6,'(A,I10,A,I10)')&
            ' NUMBER OF STATES INCONSISTENT: READ',&
            nstate0,' PRESENT ',nstate
       CALL stopgm('READ_MATRIX','NUMBER OF STATES',& 
            __LINE__,__FILE__)
    ENDIF
    IF (paral%io_parent)&
         READ(10) rotmat

    IF (paral%io_parent)&
         CALL fileclose(10)
    RETURN
  END SUBROUTINE r_matrix
  ! ==================================================================


END MODULE g_loc_util_utils
