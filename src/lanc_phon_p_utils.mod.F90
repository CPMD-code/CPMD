MODULE lanc_phon_p_utils
  USE coor,                            ONLY: tau0
  USE cotr,                            ONLY: cotc0
  USE d_mat_p_utils,                   ONLY: d_mat_diag_nonloc,&
                                             d_mat_diag_real,&
                                             d_mat_loc0,&
                                             d_mat_real
  USE eicalc_utils,                    ONLY: eicalc1
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE fft_maxfft,                      ONLY: maxfft
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_def,&
                                             fo_old
  USE hess_eta_p_utils,                ONLY: hess_eta_p
  USE implhv,                          ONLY: rs_v,&
                                             sd0,&
                                             xma
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE kpnt,                            ONLY: wk
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum
  USE nlps,                            ONLY: imagp,&
                                             ndfnl
  USE parac,                           ONLY: parai,&
                                             paral
  USE phonons_p_utils,                 ONLY: setrot,&
                                             settras
  USE prng_utils,                      ONLY: repprngu_vec
  USE readsr_utils,                    ONLY: xstring
  USE response_pmod,                   ONLY: dfnl00,&
                                             eig1,&
                                             eig2,&
                                             fnl00,&
                                             lancphon,&
                                             response1,&
                                             rho0
  USE rmas,                            ONLY: rmass
  USE rnlsm_p_utils,                   ONLY: rnlsm3
  USE rnlsm_utils,                     ONLY: rnlsm
  USE sfac,                            ONLY: ddfnl,&
                                             dfnl,&
                                             fnl
  USE soft,                            ONLY: soft_com
  USE softex_utils,                    ONLY: softex
  USE system,                          ONLY: fpar,&
                                             maxsys,&
                                             ncpw,&
                                             nkpt,&
                                             parap
  USE testex_utils,                    ONLY: testex
  USE timer,                           ONLY: tiset
!!use d_mat_p_utils, only : d_mat_diag_loc
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: lanc_phon_p
  PUBLIC :: no_mode
  PUBLIC :: setsd0
  !public :: vibeig_lan
  !public :: cont_lanczos

CONTAINS

  ! =================================================================
  SUBROUTINE lanc_phon_p(c0,c1,psi,rhoe,drhoe,&
       eirop,eivps,&
       z11,nstate)
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: psi(:)
    REAL(real_8)                             :: rhoe(*), drhoe(fpar%nnr1)
    COMPLEX(real_8)                          :: eirop(ncpw%nhg), &
                                                eivps(ncpw%nhg)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: z11(nstate,nstate)
    COMPLEX(real_8)                          :: c1(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'lanc_phon_p'

    CHARACTER(len=30)                        :: converged
    COMPLEX(real_8), ALLOCATABLE             :: pippo(:)
    INTEGER :: i, ia, iat, ie, ierr, ii, iii, info, is, istate, isub, iter, &
      iter_ini, ix, j, k, lanczos_dim_old, lastconverged, lastvalid, lda, &
      ldfnl, lfnl, nodim_old
    REAL(real_8)                             :: ddot, norm_l, xmass
    REAL(real_8), ALLOCATABLE :: alpha(:), beta(:), diagonals(:), eig3n(:,:), &
      eigenvect(:,:), rot(:,:), subdiagonals(:), tmp(:), tras(:,:), work(:)

! ==--------------------------------------------------------------==
! ==--------------------------------------------------------------==
! ==--------------------------------------------------------------==
! ==--------------------------------------------------------------==
! ==--------------------------------------------------------------==
! Here, we store the eigensystem of the tridiagonal
! Lanczos matrix. Thus the dimension is lanczos_dim x lanczos_dim.
! ==--------------------------------------------------------------==
! ==--------------------------------------------------------------==
! Start of executable statements  *******************************
! ==--------------------------------------------------------------==
! DFNL

    ndfnl=parap%nst12(parai%mepos,2)-parap%nst12(parai%mepos,1)+1
    ldfnl=imagp*3*ions1%nat*maxsys%nhxs*ndfnl*nkpt%nkpnt
    IF (ldfnl.LE.0) ldfnl=1

    CALL tiset('    lanc_phon',isub)

    ALLOCATE(rs_v(3*ions1%nat,eig1%lanczos_dim*ions1%nat/ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(eigenvect(eig1%lanczos_dim,eig1%lanczos_dim),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(work(eig1%lanczos_dim*eig1%lanczos_dim),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(alpha(eig1%lanczos_dim),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(beta(eig1%lanczos_dim),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(diagonals(eig1%lanczos_dim),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(subdiagonals(eig1%lanczos_dim),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(tmp(3*ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(xma(3*ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    ALLOCATE(fnl00(imagp,ions1%nat,maxsys%nhxs,nstate,nkpt%nkpnt),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(dfnl00(imagp,ions1%nat,maxsys%nhxs,3,nstate,nkpt%nkpnt),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(ddfnl(3*ldfnl,1,1,1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(sd0(3*ions1%nat,3*ions1%nat*ions1%nat/ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(tras(3*ions1%nat,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(rot(3*ions1%nat,3),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(eig3n(3*ions1%nat,eig1%lanczos_dim*ions1%nat/ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(pippo(maxfft),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ----------------------------------------------------------------------

    IF (paral%parent) THEN
       ! mass weights for force constants
       ix=0
       DO is=1,ions1%nsp
          xmass=SQRT(1._real_8/rmass%pma0(is))
          DO ia=1,ions0%na(is)
             DO k=1,3
                ix=ix+1
                xma(ix)=xmass
             ENDDO
          ENDDO
       ENDDO

    ENDIF


    CALL mp_bcast(xma,SIZE(xma),parai%source,parai%allgrp)


    ! ... definition of the vectors describing ... 
    ! ... center of mass translation, tras(3*nat,3)
    CALL settras(tras)

    ! ... and elemental rotations of the molecule wrt the axes x,y,z
    CALL setrot(rot,tras)


    ! calculate the derivative of dfnl ddfnl     
    CALL rnlsm3(c0,nstate,ddfnl)

    ! ==--------------------------------------------------------------==
    ! NB[dsebasti]: the unperturbed density is always already in
    ! the common block variable "rho0"

    ! FNL,DFNL,DDFNL...
    ndfnl=parap%nst12(parai%mepos,2)-parap%nst12(parai%mepos,1)+1
    ldfnl=imagp*3*ions1%nat*maxsys%nhxs*ndfnl*nkpt%nkpnt
    IF (ldfnl.LE.0) ldfnl=1
    lfnl = imagp*nstate*ions1%nat*maxsys%nhxs*nkpt%nkpnt

    CALL rnlsm(c0,nstate,1,1,.TRUE.)
    CALL dcopy(lfnl,   fnl,1,  fnl00,1)
    CALL dcopy(ldfnl, dfnl,1, dfnl00,1)


    ! ... Part of the Hessian matrix depending *only* con c0
    CALL zeroing(sd0)!,3*3*ions1%nat*ions1%nat)
    CALL setsd0(eirop,ddfnl,nstate,psi)

    CALL zeroing(rs_v)!,eig1%lanczos_dim*3*ions1%nat)
    CALL zeroing(beta)!,eig1%lanczos_dim)
    CALL zeroing(alpha)!,eig1%lanczos_dim)

    ! ==--------------------------------------------------------------==
    ! Initialization: first run or continuation   ********************
    ! ==--------------------------------------------------------------==
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,'(2x,A,29x,i6)') 'LANCZOS DIMENSION ',eig1%lanczos_dim-1
       IF (paral%io_parent)&
            WRITE(6,'(2x,a,16x,i6)')'NUMBER OF REQUESTED ITERATIONS ',&
            lancphon%nlan_st
       IF (paral%io_parent)&
            WRITE(6,'(2x,a,12x,e10.4)') 'CONVERGENCE THRESHOLD ON'//&
            ' EIGENVECTORS ',eig2%conv_threshold
       IF (eig1%lanczos_dim-1.LT.lancphon%nlan_st) THEN
          IF (paral%io_parent)&
               WRITE(6,'(2x,a)') '! !!!! ITERATIONS REQUESTED EXCEED '//&
               'LANCZOS DIMENSION ! !!!!'
          CALL stopgm('LANC_PHON','TOO MANY ITERATIONS',& 
               __LINE__,__FILE__)
       ENDIF
    ENDIF
    IF (response1%tlanph_cont) THEN
       IF (paral%parent) THEN
          nodim_old=cotc0%nodim
          CALL cont_lanczos(lanczos_dim_old,iter_ini,nodim_old,&
               iat,rs_v,alpha,beta,'READ')
          IF (paral%io_parent)&
               WRITE(6,'(/,2x,a,21x,i6)') 'RESTARTING WITH ITERATION ',&
               iter_ini
          IF (lanczos_dim_old.NE.eig1%lanczos_dim) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(2x,a,2i6)') 'WARNING LANCZOS DIMENSION CHANGED'&
                  ,eig1%lanczos_dim-1,lanczos_dim_old-1
          ELSEIF (nodim_old.NE.cotc0%nodim) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(2x,a,2i6)')'WARNING DEG. OF FREEDOM CHANGED '&
                  ,cotc0%nodim,nodim_old
          ENDIF
          IF (lancphon%nlan_st+iter_ini.GT.eig1%lanczos_dim) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(2x,a)') '******** WARNING *********'
             IF (paral%io_parent)&
                  WRITE(6,'(2x,a,15x,I6)')'THE REQUESTED LANCZOS '//&
                  'ITERATIONS'  ,lancphon%nlan_st
             IF (paral%io_parent)&
                  WRITE(6,'(2X,A,15x,I6)')'EXCEED LANCZOS SPACE '//&
                  'DIMENSIONS ',eig1%lanczos_dim-1
             lancphon%nlan_st=eig1%lanczos_dim-iter_ini
             IF (paral%io_parent)&
                  WRITE(6,'(2X,A,27x,I6)')'ITERATIONS RESET TO '&
                  ,lancphon%nlan_st
          ENDIF
       ENDIF
       CALL mp_bcast(iter_ini,parai%source,parai%allgrp)
       CALL mp_bcast(iat,parai%source,parai%allgrp)
       CALL mp_bcast(lancphon%nlan_st,parai%source,parai%allgrp)
       CALL mp_bcast(rs_v,cotc0%nodim * iter_ini,parai%source,parai%allgrp)

       ! ... exit the cycle without overwriting VIBEIGVEC and LANCZOS_CONTINUE.1
       IF (lancphon%nlan_st.LE.0) GOTO 558

    ELSE
       iter_ini=1
       beta(1) = 0._real_8       ! for security (not really needed).
       iat=0

       ! ...  initialize the first response_vector randomly (unit norm)

       IF (paral%parent) THEN
          CALL repprngu_vec(3*ions1%nat,rs_v) 

          ! ... projecting out rotational and translational degrees of freedom  
          IF (response1%projout) CALL no_mode(tras,rs_v(1,1),ions1%nat)
          IF (response1%rotout)  CALL no_mode(rot, rs_v(1,1),ions1%nat)

          norm_l = SQRT(ddot(3*ions1%nat,rs_v(1,1),1,rs_v(1,1),1))
          CALL dscal(3*ions1%nat, 1._real_8/norm_l, rs_v(1,1),1)
       ENDIF

       CALL mp_bcast(rs_v,3*ions1%nat,parai%source,parai%allgrp)
    ENDIF
    ! All other arrays are expected to be correctly initialized by the calling
    ! (and allocating) routine, do_perturbation_p().

    ! i = 1                     ! dummy. Serve per la rnlrh1.


    ! Do the Lanczos steps.
    ! ==--------------------------------------------------------------==
    ! Start Lanczos cycle    *****************************************
    ! ==--------------------------------------------------------------==
    DO iter=iter_ini,iter_ini+lancphon%nlan_st-1

       ! ...  first normalization with the masses
       CALL zeroing(tmp)!,3*ions1%nat)
       ! ........................rs_v(1,iter)|---> tmp
       CALL dcopy(3*ions1%nat,rs_v(1,iter),1,tmp,1)
       !$omp parallel do private(j)
       DO j=1,3*ions1%nat
          tmp(j)=tmp(j)*xma(j)
       ENDDO

       ! ...  Routine hess_eta; calculates implicitly the matrix-vector product
       ! ...  v2 = H v1
       ! ...  Input;  rs_v(1,iter)   is v1.
       ! ...  Output; rs_v(1,iter+1) is v2.
       ! ...  The calling is made explicitly with all the matrix rs_v
       IF (MOD(iter,3).EQ.1) iat=iat+1

       ! ==--------------------------------------------------------------==
       CALL hess_eta_p(c0,c1,psi,rhoe,drhoe,&
            eirop,eivps,z11,nstate,iter,iat,ddfnl,tmp)
       ! ==--------------------------------------------------------------==

       ! ...  second normalization by the masses
       !$omp parallel do private(j) shared(iter)
       DO j=1,3*ions1%nat
          rs_v(j,iter+1)=rs_v(j,iter+1)*xma(j)
       ENDDO

       CALL mp_sum(rs_v(:,iter+1),3*ions1%nat,parai%allgrp)

       IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)
       IF (soft_com%exsoft) THEN
          IF (paral%parent) THEN
             CALL cont_lanczos(eig1%lanczos_dim,iter,cotc0%nodim,iat,rs_v,&
                  alpha,beta,'WRITE')
             CALL softex(0)
          ENDIF
          GOTO 558
       ENDIF

       IF (paral%parent) THEN
          norm_l =SQRT(ddot(3*ions1%nat, rs_v(1,iter+1),1,rs_v(1,iter+1),1))
          IF (paral%io_parent)&
               WRITE(6,*) ' '
          IF (paral%io_parent)&
               WRITE(6,*)'******************************',&
               '**********************************'
          IF (paral%io_parent)&
               WRITE(6,*)' *+*   L2 norm[n[',iter+1,']]     = ', norm_l

          ! ==--------------------------------------------------------------==
          ! Proceed with the Lanczos algorithm:
          ! A v1 (the answer) is in response(:,iter+1) == v2:
          ! v2   = A v1
          ! ==--------------------------------------------------------------==

          CALL zeroing(tmp)!,3*ions1%nat)
          CALL dcopy(3*ions1%nat, rs_v(1,iter+1),1,tmp,1)
          ! tmp  = v2

          alpha(iter)=ddot(3*ions1%nat,rs_v(1,iter),1,tmp,1)
          ! a1   = v1 A v1
          ! \    = v1 tmp
          ! \    This is the projection of the input vector (iter)
          ! \    onto the reponse vector (iter+1). -> Orthogonalization.

          CALL daxpy(3*ions1%nat, - alpha(iter), rs_v(1,iter), 1,&
               rs_v(1,iter+1), 1)
          ! v2   = v2  -  a1 v1

          IF (.NOT. (iter .EQ. 1))&
               CALL daxpy(3*ions1%nat, - beta(iter),rs_v(1,iter-1), 1,&
               rs_v(1,iter+1), 1)
          ! v2   = v2  -  b1 v0

          ! ...  projecting out rotational and translational degrees of freedom
          IF (response1%projout) CALL no_mode(tras,rs_v(1,iter+1),ions1%nat)
          IF (response1%rotout) CALL no_mode(rot,rs_v(1,iter+1),ions1%nat)

          ! *** Here, we do a complete orthogonalization of the new
          ! \   vector (tmp=A v_i) wrt the old ones (v_j, j=1..i-1).
          IF (iter .GE. 2) THEN
             DO i=1,iter
                norm_l=-ddot(3*ions1%nat,rs_v(1,iter+1),1,rs_v(1,i),1)
                CALL daxpy(3*ions1%nat,norm_l,rs_v(1,i),1,rs_v(1,iter+1),1)
             ENDDO
          ENDIF

          norm_l=SQRT(ddot(3*ions1%nat,rs_v(1,iter+1),1,rs_v(1,iter+1),1))

          beta(iter+1) = norm_l

          ! b2   =  sqrt  v2 A v1  =  sqrt  v2 tmp

          CALL dscal(3*ions1%nat, 1._real_8/norm_l,rs_v(1,iter+1), 1)
          ! v2 is normalized.
       ENDIF


       CALL mp_bcast(rs_v(:,iter+1),3*ions1%nat,parai%source,parai%allgrp)

       IF (paral%parent) THEN
          IF (paral%io_parent)&
               WRITE(6,*)' *=*   overlap: alpha[',&
               iter,'] = ',alpha(iter)
          IF (paral%io_parent)&
               WRITE(6,*)' *=*   off-diag: beta[',&
               iter+1,'] = ',beta(iter+1)
          IF (paral%io_parent)&
               WRITE(6,*)' *=*   norm:               = ',norm_l

          ! ==--------------------------------------------------------------==

          IF (paral%io_parent)&
               WRITE(6,*)'******************************',&
               '**********************************'
          CALL dcopy(iter, alpha, 1, diagonals, 1)
          CALL dcopy(iter-1, beta(2), 1, subdiagonals, 1)

          ! Now, the diagonal [alpha(*)] and the subdiagonal [beta(*)]
          ! elements of the (tridiagonal) Lanczos matrix 
          ! are calculated. This matrix is now diagonalized:

          CALL dstev('V', iter, diagonals, subdiagonals,&
               eigenvect, eig1%lanczos_dim, work, info)
          IF ((info .NE. 0).AND.paral%io_parent)&
               WRITE(6,*)'DSTEV exit status: ', info

          ! output of dstev:
          ! \  alpha = eigenvalues,
          ! \  eigenvect = corresponding eigenvectors.

          IF (paral%io_parent)&
               WRITE(6,'(1x,a,i3,a)') '*** SPECTRUM, cycle',iter,':'
          lastconverged = iter
          DO i=iter,1,-1
             IF ( ABS(eigenvect(iter,i)) .LT. eig2%conv_threshold) THEN
                converged = '       (converged: '
                IF ( lastconverged .EQ. i+1)&
                     lastconverged = i
             ELSE
                converged = '   (NOT converged: '
             ENDIF
             CALL xstring(converged,ia,ie)
             PRINT '(A,I3,A,F12.7,A,F8.6,A)',&
                  ' *** eigenvalue ',i,' = ',diagonals(i),&
                  converged(ia:ie),ABS(eigenvect(iter,i)),').'
          ENDDO

          DO i=1,iter
             diagonals(i)=SIGN(5140.487_real_8*&
                  SQRT(ABS(diagonals(i))),diagonals(i))
          ENDDO
          IF (paral%io_parent)&
               WRITE(6,'(/," ",64("*"),/,a,/)')&
               ' harmonic frequencies [cm**-1]:'

          IF (paral%io_parent)&
               WRITE(6,'(4(f16.4))') (diagonals(i),i=1,iter)
          IF (paral%io_parent)&
               WRITE(6,*) ' '
          IF (paral%io_parent)&
               WRITE(6,*)'******************************',&
               '**********************************'
       ENDIF
       ! ... writing on file the results of the job
       IF (paral%parent) THEN
          CALL cont_lanczos(eig1%lanczos_dim,iter,cotc0%nodim,iat,rs_v,&
               alpha,beta,'WRITE')
       ENDIF
       IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)
       IF (soft_com%exsoft) THEN
          CALL softex(0)
          GOTO 558
       ENDIF
    ENDDO


    IF (paral%parent) THEN
       CALL cont_lanczos(eig1%lanczos_dim,iter,cotc0%nodim,iat,rs_v,alpha,&
            beta,'WRITE')
    ENDIF


    ! ==--------------------------------------------------------------==
    ! == END / Lanczos runs.
    ! ==--------------------------------------------------------------==

    lastvalid = iter_ini+lancphon%nlan_st-1

    IF (paral%parent) THEN
       DO istate=lastvalid,1, -1
          CALL zeroing(tmp)!,3*ions1%nat)
          ! build the eigenvector of the linear answer machinery:
          IF (3*ions1%nat.GE.lastvalid) THEN
             lda=3*ions1%nat
          ELSE
             lda=lastvalid
          ENDIF
          CALL dgemv('N', 3*ions1%nat, lastvalid, 1._real_8, rs_v, lda,&
               eigenvect(1,istate),1, 0._real_8, eig3n(1,istate), 1)


          IF (lancphon%details) THEN
             DO i=1,lastvalid
                PRINT '(A,I2,A,I2,A,F15.11)',&
                     '  \*  Lanczos(',i,'; ',istate,') = ',&
                     eigenvect(i,istate)
             ENDDO
             IF (paral%io_parent)&
                  WRITE(6,*) ' '
             DO i=1,3*ions1%nat
                IF (paral%io_parent)&
                     WRITE(6,*) i,istate,eig3n(i,istate)
             ENDDO
             IF (paral%io_parent)&
                  WRITE(6,*) ' '
          ENDIF
       ENDDO
       CALL vibeig_lan(diagonals,eig3n,3*ions1%nat,eig1%lanczos_dim)
    ENDIF


    ! ==--------------------------------------------------------------==
    ! == END / Output of eigensolutions.
    ! ==--------------------------------------------------------------==

    ! ==--------------------------------------------------------------==
    ! == CHECK whether the Lanczos and Ritz vectors are orthonormal:
    ! ==--------------------------------------------------------------==

    IF (paral%parent.AND.lancphon%details) THEN
       DO i=1,lastvalid
          ! check orthonormality of Ritz i against Ritz i, i+1, i+2,...:
          DO ii=i,lastvalid
             norm_l = 0._real_8
             DO iii=1,lastvalid
                norm_l = norm_l + eigenvect(iii,i)*eigenvect(iii,ii)
             ENDDO
             PRINT '(A,I2,A,I2,A,F15.11)',&
                  ' *** Ritz(',i,') . Ritz(',ii,') =     ',norm_l

             ! check orthonormality of Krylov i against Krylov i, i+1, i+2,...:
             norm_l = ddot(3*ions1%nat,rs_v(1,i),1,rs_v(1,ii),1)
             PRINT '(A,I2,A,I2,A,F15.11)',&
                  ' *** Krylov(',i,') . Krylov(',ii,') = ',norm_l
          ENDDO
       ENDDO
    ENDIF

    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==

558 CONTINUE

    DEALLOCATE(tmp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rs_v,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(eigenvect,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(work,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(alpha,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(beta,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(diagonals,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(subdiagonals,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(fnl00,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(dfnl00,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(ddfnl,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    DEALLOCATE(sd0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(xma,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(tras,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rot,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(eig3n,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    DEALLOCATE(pippo,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    RETURN
  END SUBROUTINE lanc_phon_p

  ! ==--------------------------------------------------------------==
  SUBROUTINE no_mode(mode,epsilon,nat)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: EPSILON(*)
    INTEGER                                  :: nat
    REAL(real_8)                             :: mode(3*nat,3)

    INTEGER                                  :: i
    REAL(real_8)                             :: alfa, ddot

    DO i=1,3
       alfa = - ddot(3*nat,epsilon,1,mode(1,i),1)
       CALL daxpy(3*nat,alfa,mode(1,i),1,epsilon,1)
    ENDDO

    RETURN
  END SUBROUTINE no_mode

  ! ==--------------------------------------------------------------==
  SUBROUTINE setsd0(eirop,ddfnl,nstate,psi)
    ! ==--------------------------------------------------------------==


    COMPLEX(real_8)                          :: eirop(ncpw%nhg)
    INTEGER                                  :: nstate
    REAL(real_8) :: ddfnl(imagp,ions1%nat,maxsys%nhxs,3,3,nstate,nkpt%nkpnt)
    COMPLEX(real_8)                          :: psi(:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'setsd0'

    COMPLEX(real_8), ALLOCATABLE             :: eirop1(:), v1_loc(:)
    INTEGER                                  :: ia, ia2, iat, iat2, icol, &
                                                ierr, is, is2, ix, k, k2

    ALLOCATE(eirop1(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(v1_loc(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)


    ! ... computation of the part of the Hessian matrix independent on the
    ! ... perturbed wavefunctions. Only rhoe0 enters here. Elements stored 
    ! ... in the matrix sd0(3*nat,3*nat).

    CALL ffttog(rho0,rho0,psi,ncpw%nhg,.TRUE.)
    ! ...because d_mat_diag_loc wants rho^0 in G space!
    iat=0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          iat=iat+1
          DO k=1,3

             CALL zeroing(v1_loc)!,nhg)
             CALL zeroing(eirop1)!,nhg)
             ! ..   calculate the local perturbation wrt. to atom iat and coordinate k
             CALL eicalc1(k,is,iat,v1_loc,eirop1)

             icol=3*(iat-1)+k
             iat2=0
             DO is2=1,ions1%nsp
                DO ia2=1,ions0%na(is2)
                   iat2=iat2+1
                   DO k2=1,3
                      ix=3*(iat2-1)+k2

                      IF (iat.EQ.iat2) THEN
                         ! local part, i=j (i,j=atoms)
                         CALL d_mat_diag_loc(sd0(ix,icol),rho0,&
                              eirop,eirop1,v1_loc,k2) !vw removed is
                         ! non-local part, i=j (i,j=atoms)
                         CALL d_mat_diag_nonloc(sd0(ix,icol),&
                              ddfnl,fnl00,&
                              dfnl00,crge%f,wk,is,k,k2,iat,nstate,1)
                         ! real space term, i=j
                         CALL d_mat_diag_real(sd0(ix,icol),is,ia,k,k2,&
                              tau0)
                         ! real space term, i <> j
                      ELSE
                         CALL d_mat_real(sd0(ix,icol),is,ia,k,k2,&
                              is2,ia2,tau0)
                      ENDIF

                      ! local part        
                      CALL d_mat_loc0(sd0(ix,icol),eirop1,&
                           iat2,k2,is2)
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    CALL ffttor(rho0,rho0,psi,ncpw%nhg,.TRUE.)

    DEALLOCATE(eirop1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(v1_loc,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    RETURN
  END SUBROUTINE setsd0

  ! ==================================================================
  SUBROUTINE vibeig_lan(vibe,sder,ndim,lanczos_dim)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ndim, lanczos_dim
    REAL(real_8)                             :: sder(ndim,lanczos_dim), &
                                                vibe(lanczos_dim)

    CHARACTER(len=100)                       :: filen
    INTEGER                                  :: i, j, k, ncol, nn
    LOGICAL                                  :: ferror

    ncol=8
    filen='VIBEIGVEC'
    IF (paral%io_parent)&
         CALL fileopen(21,filen,fo_def,ferror)
    DO i=1,lanczos_dim-1,ncol
       nn=MIN(ncol,lanczos_dim-i)
       IF (paral%io_parent)&
            WRITE(21,'(14I12)') (j,j=i,i+nn-1)
       IF (paral%io_parent)&
            WRITE(21,'(14(F12.3))') (vibe(j),j=i,i+nn-1)
       IF (paral%io_parent)&
            WRITE(21,*)
       DO k=1,ndim
          IF (paral%io_parent)&
               WRITE(21,'(14F12.6)') (sder(k,j),j=i,i+nn-1)
       ENDDO
    ENDDO
    IF (paral%io_parent)&
         CALL fileclose(21)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE vibeig_lan
  ! ==================================================================


  ! ==================================================================
  SUBROUTINE cont_lanczos(lanczos_dim,iter,ndim,iat,rs_v,alpha,beta,&
       tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lanczos_dim, iter, ndim, iat
    REAL(real_8)                             :: rs_v(ndim,*), alpha(*), &
                                                beta(*)
    CHARACTER(len=*)                         :: tag

    CHARACTER(len=100)                       :: filen
    INTEGER                                  :: i, iunit, k
    LOGICAL                                  :: ferror

    iunit=21

    IF (INDEX(tag,'READ').NE.0) THEN
       filen='LANCZOS_CONTINUE'
       IF (paral%io_parent)&
            CALL fileopen(iunit,filen,fo_old,ferror)
       IF (ferror) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A)') ' CONT_LANCZOS! FILE LANCZOS_CONTINUE DOES'//&
               'NOT  EXIST'
          CALL stopgm('CONT_LANCZOS','CANNOT READ OLD RESPONSE '//&
               'VECTORS',& 
               __LINE__,__FILE__)
       ENDIF

       IF (paral%io_parent)&
            REWIND(iunit)
       IF (paral%io_parent)&
            READ(iunit,'(4i12)') lanczos_dim,iter,ndim,iat

       DO i=1,iter
          DO k=1,ndim
             IF (paral%io_parent)&
                  READ(iunit,*) rs_v(k,i)
          ENDDO
          IF (paral%io_parent)&
               READ(iunit,*) alpha(i),beta(i)
       ENDDO
    ELSEIF (INDEX(tag,'WRITE').NE.0) THEN
       filen='LANCZOS_CONTINUE.1'
       IF (paral%io_parent)&
            CALL fileopen(iunit,filen,fo_def,ferror)

       IF (paral%io_parent)&
            REWIND(iunit)
       IF (paral%io_parent)&
            WRITE(iunit,'(4i12)') lanczos_dim,iter,ndim,iat
       DO i=1,iter
          DO k=1,ndim
             IF (paral%io_parent)&
                  WRITE(iunit,*) rs_v(k,i)
          ENDDO
          IF (paral%io_parent)&
               WRITE(iunit,*) alpha(i),beta(i)
       ENDDO
    ELSE
       IF (paral%io_parent)&
            WRITE(6,'(A,A)') ' CONT_LANCZOS| UNKNOWN TAG: ',tag
    ENDIF

    IF (paral%io_parent)&
         CALL fileclose(iunit)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE cont_lanczos
  ! ==================================================================

END MODULE lanc_phon_p_utils
