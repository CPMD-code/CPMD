MODULE gle_utils
  USE cnst,                            ONLY: factem
  USE cotr,                            ONLY: lskcor
  USE error_handling,                  ONLY: stopgm
  USE glemod,                          ONLY: &
       gle_cp_a, gle_cpmd, gle_cust, gle_opt, gle_opt_a, gle_smart, &
       gle_smart_a, gle_white, glea, glec, glep, glepar, gles, glet, icm2au
  USE ions,                            ONLY: ions0
  USE isos,                            ONLY: isos1
  USE kinds,                           ONLY: int_4,&
                                             real_8
  USE meta_multiple_walkers_utils,     ONLY: mw_filename
  USE mp_interface,                    ONLY: mp_bcast
  USE parac,                           ONLY: parai,&
                                             paral
  USE pimd,                            ONLY: pimd1,&
                                             ipcurr,&
                                             pi_omega
  USE prng_utils,                      ONLY: repprngg
  USE store_types,                     ONLY: restart1
  USE system,                          ONLY: cntr,&
                                             cntl,&
                                             maxsys
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: gle_init
  PUBLIC :: gle_alloc
  PUBLIC :: gle_step

CONTAINS

  SUBROUTINE gle_init(tau0,velp,pmx)

    REAL(real_8)                             :: tau0(:,:,:), velp(:,:,:), pmx(*)
    INTEGER                                  :: i, ia, id, iostat, is, msglen
    INTEGER(int_4)                           :: gleok
    CHARACTER(len=80)                        :: filen
    REAL(real_8) :: lglea(glepar%gle_ns+1,glepar%gle_ns+1), &
      lglec(glepar%gle_ns+1,glepar%gle_ns+1), &
      lgles(glepar%gle_ns+1,glepar%gle_ns+1), &
      lglet(glepar%gle_ns+1,glepar%gle_ns+1), &
      ltmp(glepar%gle_ns+1,glepar%gle_ns+1), rv(glepar%gle_ns+1), sqrtm

    IF (glepar%gle_mode.EQ.0) THEN
       glepar%egle=0.0_real_8
       RETURN
    ENDIF

    ! initializes static covariance to identity (i.e., enforces FDT)
    lglec=0.0_real_8
    DO i=1,glepar%gle_ns+1
       lglec(i,i)=cntr%tempw/factem
    ENDDO

    ! initializes friction matrix, using presets or custom file
    gleok = 0
    IF (glepar%gle_mode .EQ. gle_white) THEN
       IF (cntl%tpath.AND.cntl%tpimd) glepar%gle_omega=pi_omega(ipcurr)
       lglea(1,1)=glepar%gle_omega*icm2au
    ELSEIF (glepar%gle_mode .EQ. gle_opt) THEN
       IF (cntl%tpath.AND.cntl%tpimd) glepar%gle_omega=pi_omega(ipcurr)
       lglea=gle_opt_a*glepar%gle_omega*0.01_real_8*icm2au
    ELSEIF (glepar%gle_mode .EQ. gle_cpmd) THEN
       IF (cntl%tpath.AND.cntl%tpimd) glepar%gle_omega=pi_omega(ipcurr)
       lglea=gle_cp_a*glepar%gle_omega*icm2au
    ELSEIF (glepar%gle_mode .EQ. gle_smart) THEN
       IF (cntl%tpath.AND.cntl%tpimd) glepar%gle_omega=pi_omega(ipcurr)
       lglea=gle_smart_a*glepar%gle_omega*icm2au     
    ELSEIF (glepar%gle_mode .EQ. gle_cust) THEN 
       ! reads custom parameters on root & broadcasts 
       IF (paral%parent) THEN
          ! the A matrix is read from GLE-A file
          IF (paral%io_parent) THEN
             IF (cntl%tpath.AND.cntl%tpimd) THEN
                CALL mw_filename("GLE-A",filen,ipcurr)
                OPEN(444,file=filen,status="OLD",iostat=iostat)
             ELSE
                OPEN(444,file="GLE-A",status="OLD",iostat=iostat)
             ENDIF
             IF (iostat==0) THEN
                DO i=1,glepar%gle_ns+1              
                   READ(444,*) lglea(i,:)
                ENDDO
                CLOSE(444)          
             ELSE           
                gleok=1
             ENDIF
             ! the C matrix is read from GLE-C file, if it exists
             IF (cntl%tpath.AND.cntl%tpimd) THEN
                CALL mw_filename("GLE-C",filen,ipcurr)
                OPEN(444,file=filen,status="OLD",iostat=iostat)
             ELSE
                OPEN(444,file="GLE-C",status="OLD",iostat=iostat)
             ENDIF
             IF (iostat==0) THEN
                DO i=1,glepar%gle_ns+1
                   READ(444,*) lglec(i,:)
                ENDDO
                lglec=lglec*(1./factem)
                CLOSE(444)
             ENDIF
          ENDIF
       ENDIF
       CALL mp_bcast(gleok,parai%source,parai%allgrp)     
       IF (gleok.NE.0)  CALL stopgm("GLE_INIT",&
            "Could not read GLE-A file required for LANGEVIN CUSTOM",& 
            __LINE__,__FILE__)

       msglen=(glepar%gle_ns+1)*(glepar%gle_ns+1)
       CALL mp_bcast(lglea,msglen,parai%source,parai%allgrp)
       CALL mp_bcast(lglec,msglen,parai%source,parai%allgrp)
    ENDIF
    ! computes deterministic and stochastic parts of the propagator
    CALL matexp(lglea*(-cntr%delt_ions*0.5),lglet,glepar%gle_ns+1,15,15)
    IF (lglec(1,1).NE.0.0_real_8) THEN
       ltmp=lglec-MATMUL(lglet,MATMUL(lglec,TRANSPOSE(lglet)))
       CALL scholesky(ltmp,lgles,glepar%gle_ns+1)
    ELSE
       IF (paral%io_parent)&
            WRITE(6,*) " Zero-temperature simulation"
       lgles=0.0_real_8
    ENDIF

    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,*) " RUNNING GENERALIZED LANGEVIN DYNAMICS WITH:"
       IF (paral%io_parent)&
            WRITE(6,*) " A MATRIX: "
       DO i=1,glepar%gle_ns+1
          IF (paral%io_parent)&
               WRITE(6,*) lglea(i,:)
       ENDDO
       IF (paral%io_parent)&
            WRITE(6,*) " C MATRIX: "
       DO i=1,glepar%gle_ns+1
          IF (paral%io_parent)&
               WRITE(6,*) lglec(i,:)
       ENDDO
    ENDIF

    IF (paral%parent) THEN
       IF (.NOT. restart1%rgle) THEN
          glepar%egle=0.0
          IF (paral%io_parent)&
               WRITE (6,*) " Initializing GLE momenta"
          IF (lglec(1,1).NE.0.0_real_8) THEN
             CALL scholesky(lglec,ltmp,glepar%gle_ns+1)! used for initialization of velocities
             DO is=1,maxsys%nsx
                sqrtm=SQRT(pmx(is))
                DO ia=1,ions0%na(is)
                   DO id=1,3
                      DO i=1,glepar%gle_ns+1
                         rv(i)=repprngg()
                      ENDDO
                      rv=MATMUL(ltmp,rv)
                      DO  i=1,glepar%gle_ns+1
                         glep(id,ia,is,i)=rv(i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             CALL zeroing(glep)!,maxsys%nsx*maxsys%nax*3*(glepar%gle_ns+1))
          ENDIF
          IF (glepar%gle_com.GT.0) CALL gle_rmcom(pmx)
          IF (glepar%gle_com.GT.0 .AND. isos1%tisos) CALL gle_rmrot(tau0,pmx)
       ENDIF

       ! if not restarting velocities, initializes them
       IF (.NOT. restart1%rvel) THEN
          DO is=1,maxsys%nsx
             sqrtm=SQRT(pmx(is))
             DO ia=1,ions0%na(is)
                velp(:,ia,is)=glep(:,ia,is,1)/sqrtm
             ENDDO
          ENDDO
       ENDIF
    ENDIF

    is=maxsys%nsx*maxsys%nax*3*(glepar%gle_ns+1)
    CALL mp_bcast(glep,is,parai%source,parai%allgrp)
    is=maxsys%nsx*maxsys%nax*3
    CALL mp_bcast(velp,is,parai%source,parai%allgrp)

    CALL dcopy((glepar%gle_ns+1)*(glepar%gle_ns+1),lglea,1,glea,1)
    CALL dcopy((glepar%gle_ns+1)*(glepar%gle_ns+1),lglec,1,glec,1)
    CALL dcopy((glepar%gle_ns+1)*(glepar%gle_ns+1),lgles,1,gles,1)
    CALL dcopy((glepar%gle_ns+1)*(glepar%gle_ns+1),lglet,1,glet,1)
  END SUBROUTINE gle_init

  SUBROUTINE gle_step(tau0,velp,pmx)

    REAL(real_8)                             :: tau0(:,:,:), velp(:,:,:), pmx(*)
    INTEGER                                  :: i, ia, iat, id, is, l
    REAL(real_8) :: degle, lgles(glepar%gle_ns+1,glepar%gle_ns+1), &
      lglet(glepar%gle_ns+1,glepar%gle_ns+1), rv(glepar%gle_ns+1), sqrtm, &
      xv(glepar%gle_ns+1)

    IF (glepar%gle_mode.EQ.0) RETURN

    CALL dcopy((glepar%gle_ns+1)*(glepar%gle_ns+1),gles,1,lgles,1)
    CALL dcopy((glepar%gle_ns+1)*(glepar%gle_ns+1),glet,1,lglet,1)

    ! transforms momenta to mass-scaled momenta so the same propagato can be used for all particles
    degle=0
    iat = 0
    DO is=1,maxsys%nsx
       sqrtm=SQRT(pmx(is))
       DO ia=1,ions0%na(is)
          iat = iat+1
          DO l=1,3
             IF (lskcor(l,iat).NE.0) THEN ! decouples the GLE from fixed coordinates
                glep(l,ia,is,1)=velp(l,ia,is)*sqrtm
                ! also accumulates initial kinetic energy
                degle=degle+glep(l,ia,is,1)**2
             END IF
          ENDDO
       ENDDO
    ENDDO

    ! performs evolution
    DO is=1,maxsys%nsx
       DO ia=1,ions0%na(is)
          DO id=1,3
             ! evolution vector is read in (p, s_1, .., s_N) for the selected DOF
             DO i=1,glepar%gle_ns+1
                xv(i)=glep(id,ia,is,i)
             ENDDO

             IF (glec(1).EQ.0.0_real_8) THEN
                xv=MATMUL(lglet,xv)
             ELSE
                ! fills up a vector of gaussian numbers
                DO i=1,glepar%gle_ns+1
                   rv(i)=repprngg()
                ENDDO

                ! vector is evolved by p <- T p + S xi
                xv=MATMUL(lglet,xv)+MATMUL(lgles,rv)
             ENDIF

             ! data is saved back
             DO i=1,glepar%gle_ns+1
                glep(id,ia,is,i)=xv(i)
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    ! removes COM velocity
    IF (glepar%gle_com.GT.0) CALL gle_rmcom(pmx)
    IF (glepar%gle_com.GT.0 .AND. isos1%tisos) CALL gle_rmrot(tau0,pmx)
    ! scale back velocities & updates Dconserved      
    iat = 0
    DO is=1,maxsys%nsx
       sqrtm=SQRT(pmx(is))
       DO ia=1,ions0%na(is) 
          iat = iat+1
          DO l=1,3   ! skips coordinates that are fixed
             IF (lskcor(l,iat).NE.0) THEN
                velp(l,ia,is)=glep(l,ia,is,1)/sqrtm
                degle=degle-glep(l,ia,is,1)**2
             ENDIF
          ENDDO
       ENDDO
    ENDDO

    degle=degle*0.5
    glepar%egle=glepar%egle+degle
  END SUBROUTINE gle_step

  SUBROUTINE gle_alloc()
    CHARACTER(*), PARAMETER                  :: procedureN = 'gle_alloc'

    INTEGER                                  :: ierr

! TODO: GLEA, GLEC, GLES, GLET and GLEP are never deallocated?! 

    ALLOCATE(glea((glepar%gle_ns+1)*(glepar%gle_ns+1)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(glec((glepar%gle_ns+1)*(glepar%gle_ns+1)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(gles((glepar%gle_ns+1)*(glepar%gle_ns+1)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(glet((glepar%gle_ns+1)*(glepar%gle_ns+1)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(glep(3,maxsys%nax,maxsys%nsx,(glepar%gle_ns+1)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
  END SUBROUTINE gle_alloc

  SUBROUTINE gle_rmcom(pmx)

    REAL(real_8)                             :: pmx(*)
    INTEGER                                  :: i, ia, is
    REAL(real_8)                             :: com(3), mm, sqrtm

    DO i=1, glepar%gle_ns+1 ! also do it for additional momenta (for consistency and to preserve analytical predictions)
       com=0.0_real_8
       mm=0.0_real_8
       DO is=1,maxsys%nsx
          sqrtm=SQRT(pmx(is))
          mm=mm+pmx(is)*ions0%na(is)
          DO ia=1,ions0%na(is)
             com(:)=com(:)+glep(:,ia,is,i)*sqrtm
          ENDDO
       ENDDO
       com=com/mm! COM velocity
       DO is=1,maxsys%nsx
          sqrtm=SQRT(pmx(is))
          DO ia=1,ions0%na(is)
             glep(:,ia,is,i)=glep(:,ia,is,i)-com(:)*sqrtm
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE gle_rmcom

  ! removes rotational velocities, to additional momenta as well. assumes no COM velocity is present
  SUBROUTINE gle_rmrot(tau0,pmx)

    REAL(real_8)                             :: tau0(:,:,:), pmx(*)

    INTEGER                                  :: i, ia, info, ipiv(3), is, j, &
                                                work(9)
    REAL(real_8)                             :: com(3), im(3,3), im1(3,3), &
                                                lm(3), rm, taum(3,maxsys%nax,&
                                                maxsys%nsx), tm, wm(3)

    taum(:,:,:)=tau0(:,:,:)
    com=0.0_real_8
    tm=0.0_real_8
    DO is=1,maxsys%nsx
       tm=tm+ions0%na(is)*pmx(is)
       DO ia=1,ions0%na(is)
          com=com+taum(:,ia,is)*pmx(is)
       ENDDO
    ENDDO
    com=com*(1./tm)
    DO i=1,3
       taum(i,:,:)=taum(i,:,:)-com(i)
    ENDDO

    im=0.0_real_8
    DO is=1,maxsys%nsx
       DO ia=1,ions0%na(is)
          rm=(taum(1,ia,is)**2+taum(2,ia,is)**2+taum(3,ia,is)**2)
          DO i=1,3
             im(i,i)=im(i,i)+pmx(is)*(rm-taum(i,ia,is)**2)
             DO j=i+1,3
                im(i,j)=im(i,j)-pmx(is)*taum(i,ia,is)*taum(j,ia,is)
                im(j,i)=im(i,j)
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    im1=im
    CALL dgetrf(3,3,im1,3,ipiv,info)
    CALL dgetri(3,im1,3,ipiv,work,9,info)

    DO i=1, glepar%gle_ns+1
       ! get angular momentum
       lm=0.0_real_8
       DO is=1,maxsys%nsx
          rm=SQRT(pmx(is))
          DO ia=1,ions0%na(is)
             lm(1)=lm(1)+rm*&
                  (taum(2,ia,is)*glep(3,ia,is,i)-&
                  taum(3,ia,is)*glep(2,ia,is,i))
             lm(2)=lm(2)+rm*&
                  (taum(3,ia,is)*glep(1,ia,is,i)-&
                  taum(1,ia,is)*glep(3,ia,is,i))
             lm(3)=lm(3)+rm*&
                  (taum(1,ia,is)*glep(2,ia,is,i)-&
                  taum(2,ia,is)*glep(1,ia,is,i))
          ENDDO
       ENDDO
       wm=MATMUL(im1,lm)! angular velocity
       DO is=1,maxsys%nsx
          rm=SQRT(pmx(is))
          DO ia=1,ions0%na(is)
             ! gets velocity
             lm(1)=(wm(2)*taum(3,ia,is)-wm(3)*taum(2,ia,is))
             lm(2)=(wm(3)*taum(1,ia,is)-wm(1)*taum(3,ia,is))
             lm(3)=(wm(1)*taum(2,ia,is)-wm(2)*taum(1,ia,is))
             lm=lm*rm
             ! subtract (rembember GLEP contains mass-scaled momenta!)             
             glep(:,ia,is,i)=glep(:,ia,is,i)-lm
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE gle_rmrot

  ! Exponential of a NxN matrix by scaling and squaring
  SUBROUTINE matexp(m,expm,n,j,k)
    INTEGER, INTENT(in)                      :: n
    REAL(real_8), INTENT(out)                :: expm(n,n)
    REAL(real_8), INTENT(in)                 :: m(n,n)
    INTEGER, INTENT(in)                      :: j, k

    INTEGER                                  :: i, ik
    REAL(real_8)                             :: sm(n,n), tc(k)

    sm=m*1./(2.**j)
    ! taylor expansion
    tc(1)=1.
    DO i=2,k
       tc(i)=tc(i-1)*(1./i)
    ENDDO
    expm=0
    DO ik=k,1,-1
       DO i=1,n
          expm(i,i)=expm(i,i)+tc(ik)
       ENDDO
       expm=MATMUL(sm,expm)
    ENDDO
    DO i=1,n
       expm(i,i)=expm(i,i)+1.
    ENDDO
    DO ik=1,j
       expm=MATMUL(expm,expm)
    ENDDO
  END SUBROUTINE matexp

  ! brute-force "stabilized" cholesky decomposition.
  ! in practice, we compute LDL^T decomposition, and force
  ! to zero negative eigenvalues.
  SUBROUTINE scholesky(sst, s, n)
    INTEGER, INTENT(in)                      :: n
    REAL(real_8), INTENT(out)                :: s(n,n)
    REAL(real_8), INTENT(in)                 :: sst(n,n)

    INTEGER                                  :: i, j, k
    REAL(real_8)                             :: d(n,n), l(n,n)

    s=0.
    l=0.
    d=0.
    DO i=1,n
       l(i,i)=1.0
       DO j=1,i-1
          L(i,j)=SST(i,j);
          DO k=1,j-1
             l(i,j)=l(i,j)-l(i,k)*l(j,k)*d(k,k)
          ENDDO
          IF (d(j,j).NE. 0.0) THEN
             l(i,j)=l(i,j)/d(j,j)
          ELSE
             IF (paral%io_parent)&
                  WRITE(6,*) " Warning: zero eigenvalue in LDL decomposition."
             l(i,j)=0.
          ENDIF
       ENDDO
       d(i,i)=sst(i,i)
       DO k=1,i-1
          d(i,i)=d(i,i)-l(i,k)**2*d(k,k)
       ENDDO
    ENDDO
    DO i=1,n
       IF ( d(i,i).GE. 0.0_real_8 ) THEN
          d(i,i)=SQRT(d(i,i))
       ELSE
          IF (paral%io_parent)&
               WRITE(6,*) " Warning: negative eigenvalue (",d(i,i),&
               ")in LDL^T decomposition."
          d(i,i)=0.0
       ENDIF
    ENDDO
    s=MATMUL(l,d)
  END SUBROUTINE scholesky

END MODULE gle_utils
