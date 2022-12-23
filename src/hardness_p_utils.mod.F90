MODULE hardness_p_utils
  USE coor,                            ONLY: tau0
  USE dotp_utils,                      ONLY: dotp
  USE error_handling,                  ONLY: stopgm
  USE fft_maxfft,                      ONLY: maxfftn
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_def
  USE geq0mod,                         ONLY: geq0
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai,&
                                             paral
  USE readsr_utils,                    ONLY: xstring
  USE response_pmod,                   ONLY: eig1,&
                                             eig2
  USE rhoofr_p_utils,                  ONLY: rhoofr_p
  USE rhoofr_utils,                    ONLY: rhoofr
  USE rwfopt_p_utils,                  ONLY: rwfopt_p
  USE soft,                            ONLY: soft_com
  USE spin,                            ONLY: clsd
  USE system,                          ONLY: fpar,&
                                             ncpw
  USE zeroing_utils,                   ONLY: zeroing
!!use nmr_util_p_utils, only : ffttor

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: hardness_p

CONTAINS

  ! ==================================================================
  SUBROUTINE hardness_p(c0,c1,psi,rhoe,drhoe,eirop,eivps,z11,nstate)
    ! ==--------------------------------------------------------------==

    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: psi(maxfftn)
    REAL(real_8)                             :: rhoe(fpar%nnr1,clsd%nlsd), &
                                                drhoe(*)
    COMPLEX(real_8)                          :: eirop(ncpw%nhg), &
                                                eivps(ncpw%nhg)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: z11(nstate,nstate)
    COMPLEX(real_8)                          :: c1(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'hardness_p'

    CHARACTER(len=20)                        :: converged
    CHARACTER(len=3)                         :: zahl
    CHARACTER(len=80)                        :: filename
    COMPLEX(real_8)                          :: dummy(1)
    COMPLEX(real_8), ALLOCATABLE             :: response_vectors(:,:), scr(:)
    INTEGER                                  :: i, ia, iat, ie, ierr, iii, &
                                                info, ir, isp, istate, iter, &
                                                lastvalid, nconverged, nnat
    LOGICAL                                  :: ferror
    REAL(real_8) :: center(3) = (/0.0_real_8,0.0_real_8,0.0_real_8/), norm_l, &
      temp(2)
    REAL(real_8), ALLOCATABLE                :: alpha(:), beta(:), &
                                                diagonals(:), eigenvect(:,:), &
                                                subdiagonals(:), work(:)

! ==--------------------------------------------------------------==
! Local variables, external functions
! Here, we store the eigensystem of the tridiagonal
! Lanczos matrix. Thus the dimension is lanczos_dim x lanczos_dim.
! Auxiliary variables
! ==--------------------------------------------------------------==

    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,*)' *** calculating HARDNESS  dv(r)/dn(r) via Lanczos'
       IF (paral%io_parent)&
            WRITE(6,*)' *** WITH COMPLETE REORTHOGONALIZATION.'
       IF (paral%io_parent)&
            WRITE(6,*)'     info: beta calculated after ENTIRE orthogon'
       IF (paral%io_parent)&
            WRITE(6,*)' *** Maximum dim of Lanczos space:  ',eig1%lanczos_dim
       IF (paral%io_parent)&
            WRITE(6,*)' *** Number of desired eigenstates: ',&
            eig1%num_lanczos_states
       IF (paral%io_parent)&
            WRITE(6,*)' *** Threashold for convergence:    ',eig2%conv_threshold
       IF (paral%io_parent)&
            WRITE(6,*)' *** The hardness response operator is shifted by:',&
            eig2%eigenvalue_shift
       IF (paral%io_parent)&
            WRITE(6,*)' *** INVERTING the response sign.'
       IF (paral%io_parent)&
            WRITE(6,*)' *** ELIMINATING any constant charge/potential.'
    ENDIF


    IF ( eig1%num_lanczos_states+20 .GT. eig1%lanczos_dim)&
         CALL stopgm('EIGENSYS','lanczos dim too small.',& 
         __LINE__,__FILE__)
    IF ( eig1%num_lanczos_states    .LT. 1)&
         CALL stopgm('EIGENSYS','too few response states! ',& 
         __LINE__,__FILE__)
    IF ( eig2%conv_threshold        .LT. 1.e-8_real_8)&
         CALL stopgm('EIGENSYS','conv threshold too low ( <1e-8 ).',& 
         __LINE__,__FILE__)
    IF ( eig2%conv_threshold        .GT. 1.e-1_real_8)&
         CALL stopgm('EIGENSYS','conv threshold too high ( >.01 ).',& 
         __LINE__,__FILE__)


    ALLOCATE(response_vectors(ncpw%nhg,eig1%lanczos_dim),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(  response_vectors)!,   eig1%lanczos_dim*nhg)
    ALLOCATE(eigenvect(eig1%lanczos_dim,eig1%lanczos_dim),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(   eigenvect)!, eig1%lanczos_dim*eig1%lanczos_dim)
    ALLOCATE(work(eig1%lanczos_dim*eig1%lanczos_dim),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(   work)!, eig1%lanczos_dim*eig1%lanczos_dim)
    ALLOCATE(alpha(2*eig1%lanczos_dim),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(   alpha)!,         2*eig1%lanczos_dim)
    ALLOCATE(beta(2*eig1%lanczos_dim),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(   beta)!,          2*eig1%lanczos_dim)
    ALLOCATE(diagonals(2*eig1%lanczos_dim),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(   diagonals)!,     2*eig1%lanczos_dim)
    ALLOCATE(subdiagonals(2*eig1%lanczos_dim),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(   subdiagonals)!,  2*eig1%lanczos_dim)
    ! ==--------------------------------------------------------------==
    ! Initialize perturbation density to be the gnd state density:
    CALL rhoofr(c0,rhoe,psi,nstate)
    CALL ffttog(rhoe,response_vectors(1,1),psi,ncpw%nhg,.TRUE.)
    IF (geq0) response_vectors(1,1) = CMPLX(0._real_8,0._real_8,kind=real_8)
    ! All "vectors" in this context should integrate to zero.

    norm_l = dotp(ncpw%nhg,response_vectors(:,1),response_vectors(:,1))
    CALL mp_sum(norm_l,parai%allgrp)
    norm_l = SQRT(norm_l)
    CALL dscal(2*ncpw%nhg, 1._real_8/norm_l, response_vectors(1,1),1)

    ! All other arrays are expected to be correctly initialized by the calling
    ! (and allocating) routine, do_perturbation_p().

    beta(1) = 0._real_8             ! for security (not really needed).
    ! ==--------------------------------------------------------------==
    ! Do the Lanczos steps.
    iter = 1
    nconverged = 0
    i = 1                     ! dummy.
    DO WHILE ( (iter .LT. eig1%lanczos_dim-1)&
         .AND.  (eig1%num_lanczos_states .GT. nconverged)&
         .AND. .NOT. soft_com%exsoft)


       CALL rwfopt_p(c0,c1,psi,rhoe,drhoe,&
            eirop,eivps,response_vectors(1,iter),dummy,&
            z11,nstate,dummy)
       CALL rhoofr_p(c0,c1,drhoe,psi,nstate)
       CALL ffttog(drhoe,&
            response_vectors(1,iter+1), psi, ncpw%nhg, .TRUE.)
       CALL dscal(2*ncpw%nhg,-1._real_8,response_vectors(1,iter+1),1)

       IF (geq0) response_vectors(1,iter+1) = CMPLX(0._real_8,0._real_8,kind=real_8)


       temp (1) = dotp(ncpw%nhg,response_vectors(:,iter),&
            response_vectors(:,iter))
       temp (2) = dotp(ncpw%nhg,response_vectors(:,iter+1),&
            response_vectors(:,iter+1))
       CALL mp_sum(temp,2,parai%allgrp)


       ! This is to work around the Lanczos problem that eigenvalues around
       ! zero are difficult to determine. After each run of rwfopt, we 
       ! ADD -0.xx times the INPUT potential to the output. Then from all the
       ! eigenvalues we have to SUBTRACT this -0.xx. Remember that
       ! this value is contained WITH ITS SIGN in eigenvalue_shift=-0.xx:
       ! call DAXPY(2*nhg, eigenvalue_shift, 
       ! &      response_vectors(1,iter), 1,
       ! &      response_vectors(1,iter+1), 1)
       ! NB: In the hardness calculation, however, this eigenvalue_shift
       ! should be zero.


       norm_l = dotp(ncpw%nhg, response_vectors(:,iter+1),&
            response_vectors(:,iter+1))
       CALL mp_sum(norm_l,parai%allgrp)
       norm_l = SQRT(norm_l)
       IF (paral%io_parent)&
            WRITE(6,*)' *+* L2 norm[dv[i+1]]       = ',norm_l

       ! ==--------------------------------------------------------------==
       ! Proceed with the Lanczos algorithm:
       ! A v1 (the answer) is in response(:,iter+1) == v2:
       ! v2   = A v1
       ! ==--------------------------------------------------------------==

       alpha(iter) = dotp(ncpw%nhg, response_vectors(:,iter),&
            response_vectors(:,iter+1))
       CALL mp_sum(alpha(iter),parai%allgrp)
       ! a1   = v1 (A v1)
       ! \    This is the projection of the input vector (iter)
       ! \    onto the reponse vector (iter+1). -> Orthogonalization.

       CALL daxpy(2*ncpw%nhg, - alpha(iter), response_vectors(1,iter), 1,&
            response_vectors(1,iter+1), 1)
       ! v2   = v2  -  a1 v1


       IF (.NOT. (iter .EQ. 1))&
            CALL daxpy(2*ncpw%nhg, - beta(iter),&
            response_vectors(1,iter-1), 1,&
            response_vectors(1,iter+1), 1)
       ! v2   = v2  -  b1 v0


       ! *** Here, we do a complete orthogonalization of the new
       ! \   vector wrt the old ones (v_j, j=1..i-1).
       ! Method: Grahm-Schmitt "by hand".
       IF (iter .GE. 2) THEN
          DO i=1,iter
             norm_l = - dotp(ncpw%nhg, response_vectors(:,iter+1),&
                  response_vectors(:,i))
             CALL mp_sum(norm_l,parai%allgrp)
             CALL daxpy(2*ncpw%nhg, norm_l, response_vectors(1,i),1,&
                  response_vectors(1,iter+1), 1)
          ENDDO
       ENDIF


       norm_l = (dotp(ncpw%nhg,response_vectors(:,iter+1),&
            response_vectors(:,iter+1)))
       CALL mp_sum(norm_l,parai%allgrp)
       norm_l = SQRT(norm_l)
       beta(iter+1) = norm_l

       CALL dscal(2*ncpw%nhg, 1._real_8/norm_l,&
            response_vectors(1,iter+1), 1)
       ! v2 is normalized.

       IF (paral%parent) THEN
          IF (paral%io_parent)&
               WRITE(6,*)' *=*   overlap: alpha[',&
               iter,'] = ',alpha(iter)
          IF (paral%io_parent)&
               WRITE(6,*)' *=*   off-diag: beta[',&
               iter+1,'] = ',beta(iter+1)
          ! WRITE(6,*)' *=*   norm:               = ',norm_l
          IF (paral%io_parent)&
               WRITE(6,*)'##############################',&
               '##################################'
       ENDIF
       ! ==--------------------------------------------------------------==
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
            WRITE(6,*)' *** SPECTRUM, run',iter,':'
       nconverged = 0
       DO i=iter,1,-1
          IF ( ABS(eigenvect(iter,i)) .LT. eig2%conv_threshold) THEN
             converged = '       (converged:'
             nconverged = nconverged + 1
          ELSE
             converged = '   (NOT converged:'
          ENDIF
          CALL xstring(converged,ia,ie)
          IF (paral%parent)&
               PRINT '(A,I3,A,F12.7,A,F8.6,A)',&
               ' *** eigenvalue ',i,' = ',&
               diagonals(i)-eig2%eigenvalue_shift,&
               converged(ia:ie),ABS(eigenvect(iter,i)),').'
       ENDDO
       IF (paral%parent) THEN
          IF (paral%io_parent)&
               WRITE(6,*)' *** Converged states:',nconverged
          IF (paral%io_parent)&
               WRITE(6,*)' *** Eigenvalue for Lanczos ghost states is ',&
               -eig2%eigenvalue_shift
          IF (paral%io_parent)&
               WRITE(6,*)'##############################',&
               '##################################'
       ENDIF
       lastvalid = iter
       iter = iter + 1
    ENDDO
    ! ==--------------------------------------------------------------==
    ! == END / Lanczos runs.
    ! ==--------------------------------------------------------------==
    ! Find a suitable "center" of the system, for plotting:
    nnat=0
    center(1) = 0._real_8
    center(2) = 0._real_8
    center(3) = 0._real_8
    DO isp=1,ions1%nsp
       DO iat=1,ions0%na(isp)
          center(1) = center(1) + tau0(1,iat,isp)
          center(2) = center(2) + tau0(2,iat,isp)
          center(3) = center(3) + tau0(3,iat,isp)
          nnat = nnat + 1
       ENDDO
    ENDDO
    center(1) =  center(1) / REAL(nnat,kind=real_8)
    center(2) =  center(2) / REAL(nnat,kind=real_8)
    center(3) =  center(3) / REAL(nnat,kind=real_8)
    ! Center is now the arithmetic average (NOT the center of mass!!)
    ! of all atomic positions.

    IF (paral%io_parent)&
         CALL fileopen(53,'LANCZOS_COEFF',fo_def,ferror)
    IF (ferror) GOTO 90
100 CONTINUE

    ALLOCATE(scr(maxfftn),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__) ! TODO check length
    DO istate=lastvalid,1, -1
       iii=lastvalid-istate+1



       ! build the eigenvector of the linear answer machinery:
       CALL dgemv('N', 2*ncpw%nhg, lastvalid, 1._real_8, response_vectors, 2*ncpw%nhg,&
            eigenvect(1,istate),1, 0._real_8, scr, 1)


       IF (paral%parent) THEN
          IF (paral%io_parent)&
               WRITE (53,fmt='(A,I3,A,F20.12)',err=90)&
               '***** LANCZOS(',istate,') >>> ',&
               diagonals(istate)-eig2%eigenvalue_shift
          DO i=1,lastvalid
             IF (paral%io_parent)&
                  WRITE (53,'(A,I3,A,I3,A,F15.11)')&
                  '  \*  Lanczos(',i,'; ',istate,') = ',&
                  eigenvect(i,istate)
          ENDDO
       ENDIF

       IF (paral%io_parent)&
            WRITE (zahl,'(I3.3)') iii
       filename='eigst-'//zahl//'.cube'
       CALL xstring(filename,ia,ie)
       IF (paral%io_parent)&
            WRITE (6,'(A,I3,3A)')&
            ' _i_ Writing eigenstate ',istate,&
            ' to file ',filename(ia:ie),'.'
       CALL ffttor(scr,scr,psi,ncpw%nhg,.TRUE.)
       CALL cubefile(filename,scr,center,psi,.TRUE.)

       ! TESTING: Zero integral?
       temp(1) = 0._real_8
       DO ir=1,fpar%nnr1
          temp(1) = temp(1) + scr(ir)
       ENDDO
       CALL mp_sum(temp(1),parai%allgrp)
       IF (paral%io_parent)&
            WRITE (6,*) 'STATE ',istate,': integral -> ',&
            temp(1)



    ENDDO
    DEALLOCATE(scr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    IF ((paral%parent).AND.paral%io_parent)&
         CALL fileclose(53)

    ! ==--------------------------------------------------------------==
    ! == END / Output of eigensolutions.
    ! ==--------------------------------------------------------------==


    DEALLOCATE(response_vectors,STAT=ierr)
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

    RETURN
    ! ==--------------------------------------------------------------==
90  CONTINUE
    IF (paral%io_parent)&
         WRITE(6,*) 'HARDNESS_P| COULD NOT OPEN FILE "LANCZOS_COEFF"'
    CALL stopgm('HARDNESS_P','FILE OPEN ERROR',& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE hardness_p
  ! ==--------------------------------------------------------------==

END MODULE hardness_p_utils
