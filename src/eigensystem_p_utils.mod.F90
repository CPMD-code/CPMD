MODULE eigensystem_p_utils
  USE coor,                            ONLY: tau0
  USE cppt,                            ONLY: nzh
  USE dotp_utils,                      ONLY: dotp
  USE error_handling,                  ONLY: stopgm
  USE fft_maxfft,                      ONLY: maxfftn
  USE fftmain_utils,                   ONLY: fwfftn
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
  USE spin,                            ONLY: clsd
  USE system,                          ONLY: fpar,&
                                             ncpw
  USE utils,                           ONLY: zgthr
!!use densto_utils, only : densto

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: eigensystem_p

CONTAINS

  ! ==================================================================
  SUBROUTINE eigensystem_p(c0,c1,psi,rhoe,drhoe,eirop,&
       eivps,z11,nstate)
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c1(ncpw%ngw,*), psi(maxfftn)
    REAL(real_8)                             :: rhoe(fpar%nnr1,clsd%nlsd), &
                                                drhoe(*)
    COMPLEX(real_8)                          :: eirop(*), eivps(*)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: z11(nstate,*)
    COMPLEX(real_8)                          :: c0(ncpw%ngw,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'eigensystem_p'

    CHARACTER(len=2), DIMENSION(99) :: zahlenchars = (/'01','02','03','04',&
      '05','06','07','08','09','10','11','12','13','14','15','16','17','18',&
      '19','20','21','22','23','24','25','26','27','28','29','30','31','32',&
      '33','34','35','36','37','38','39','40','41','42','43','44','45','46',&
      '47','48','49','50','51','52','53','54','55','56','57','58','59','60',&
      '61','62','63','64','65','66','67','68','69','70','71','72','73','74',&
      '75','76','77','78','79','80','81','82','83','84','85','86','87','88',&
      '89','90','91','92','93','94','95','96','97','98','99'/)
    CHARACTER(len=30)                        :: converged
    CHARACTER(len=80)                        :: filename
    COMPLEX(real_8)                          :: dummy(1)
    COMPLEX(real_8), ALLOCATABLE             :: response_vectors(:,:), tmp(:)
    INTEGER                                  :: i, ia, ie, ierr, info, ir, &
                                                istate, iter, lastconverged, &
                                                lastvalid
    REAL(real_8)                             :: norm_l
    REAL(real_8), ALLOCATABLE                :: alpha(:), beta(:), &
                                                diagonals(:), eigenvect(:,:), &
                                                subdiagonals(:), work(:)

    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,*)' *** calculating EIGENSYSTEM of the ',&
            'perturbation via Lanczos'
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
            WRITE(6,*)' *** The response operator is shifted by:',&
            eig2%eigenvalue_shift
       IF (paral%io_parent)&
            WRITE(6,*)' *** INVERTING the response sign.'
       IF (paral%io_parent)&
            WRITE(6,*)' *** ELIMINATING any constant response charge.'
    ENDIF


    IF (eig1%num_lanczos_states+20 .GT. eig1%lanczos_dim) &
         CALL stopgm('EIGENSYS' ,'lanczos dim too small.',& 
         __LINE__,__FILE__)
    IF (eig1%num_lanczos_states .LE. 2) &
         CALL stopgm('EIGENSYS','too few response states! ',& 
         __LINE__,__FILE__)
    IF (eig2%conv_threshold .LT. 1.e-8_real_8) &
         CALL stopgm('EIGENSYS','conv threshold too low ( <1e-8 ).',& 
         __LINE__,__FILE__)
    IF (eig2%conv_threshold .GT. 1.e-2_real_8) &
         CALL stopgm('EIGENSYS','conv threshold too high ( >.01 ).',& 
         __LINE__,__FILE__)


    ALLOCATE(response_vectors(ncpw%nhg,2*eig1%lanczos_dim),STAT=ierr)
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
    ALLOCATE(diagonals(2*eig1%lanczos_dim),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(subdiagonals(2*eig1%lanczos_dim),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(tmp(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)


    ! ==--------------------------------------------------------------==
    ! Initialize potential
    ! is specific FOR A PARTICULAR PHYSICAL SYSTEM (not specific for a
    ! particular perturbation!), the result of the diagonalization of
    ! chi must always be the same (!)

    CALL rhoofr(c0,rhoe,psi,nstate)
    ! transform rhoe into g space
    DO ir=1,fpar%nnr1
       psi(ir)=CMPLX(rhoe(ir,1),0._real_8,kind=real_8)
    ENDDO
    CALL  fwfftn(psi,.FALSE.,parai%allgrp)
    ! Save the density response vector:
    CALL zgthr(ncpw%nhg,psi,response_vectors(1,1),nzh)


    norm_l = dotp(ncpw%nhg,response_vectors(:,1),response_vectors(:,1))
    CALL mp_sum(norm_l,parai%allgrp)
    norm_l = SQRT(norm_l)
    CALL dscal(2*ncpw%nhg, 1._real_8/norm_l, response_vectors,1)

    ! All other arrays are expected to be correctly initialized by the calling
    ! (and allocating) routine, do_perturbation_p().

    beta(1) = 0._real_8             ! for security (not really needed).





    ! ==--------------------------------------------------------------==
    ! Do the Lanczos steps.
    iter = 1
    lastconverged = 1
    DO WHILE ( (iter .LT. eig1%lanczos_dim).AND.  (iter -&
         eig1%num_lanczos_states .LT. lastconverged))

       i = 1                 ! dummy.
       CALL rwfopt_p(c0,c1,psi,rhoe,drhoe,eirop,eivps,&
            response_vectors(:,iter),dummy,z11,nstate,dummy)


       ! calculate the linear response of the density 
       ! \   drhoe(r)   =   <c0|r><r|c1> + cc
       CALL rhoofr_p(c0,c1,drhoe,psi,nstate)

       ! transform drhoe into g space
       DO ir=1,fpar%nnr1
          psi(ir)=CMPLX(drhoe(ir),0._real_8,kind=real_8)
       ENDDO
       CALL  fwfftn(psi,.FALSE.,parai%allgrp)
       ! Save the density response vector:
       CALL zgthr(ncpw%nhg,psi,response_vectors(1,iter+1),nzh)
       ! The response is defined by (-1) times what we have calculated so far:
       CALL dscal(2*ncpw%nhg,-1.0_real_8,response_vectors(1,iter+1),1)
       IF (paral%io_parent)&
            WRITE(6,*)' *+* response INVERTED.  ',&
            'Resp(i+1)   =   - Potential[ Resp(i) ]'
       ! the new response vector now contains the (packed) response
       ! of the system.


       ! this is to work around the Lanczos problem that eigenvalues around
       ! zero are difficult to determine. After each run of rwfopt, we 
       ! ADD -0.xx times the INPUT potential to the output. Then from all the
       ! eigenvalues we have to subtract this -0.xx.
       CALL daxpy(2*ncpw%nhg, eig2%eigenvalue_shift, response_vectors(1,iter), 1,&
            response_vectors(1,iter+1), 1)

       norm_l = dotp(ncpw%nhg, response_vectors(:,iter+1),response_vectors(:,&
            iter+1))
       CALL mp_sum(norm_l,parai%allgrp)
       norm_l = SQRT(norm_l)
       IF (paral%io_parent)&
            WRITE(6,*)' *+* L2 norm[n[i+1]]       = ',norm_l

       ! ==--------------------------------------------------------------==
       ! Proceed with the Lanczos algorithm:
       ! A v1 (the answer) is in response(:,iter+1) == v2:
       ! v2   = A v1
       ! ==--------------------------------------------------------------==
       CALL dcopy(2*ncpw%nhg, response_vectors(1,iter+1),1,tmp,1)
       ! tmp  = v2

       alpha(iter) = dotp(ncpw%nhg, response_vectors(:,iter), tmp)
       CALL mp_sum(alpha(iter),parai%allgrp)
       ! a1   = v1 A v1
       ! \    = v1 tmp
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
       ! \   vector (tmp=A v_i) wrt the old ones (v_j, j=1..i-1).
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



       norm_l = (dotp(ncpw%nhg,response_vectors(:,iter+1),response_vectors(:,&
            iter+1)))
       CALL mp_sum(norm_l,parai%allgrp)
       norm_l = SQRT(norm_l)
       beta(iter+1) = norm_l

       CALL dscal(2*ncpw%nhg, 1._real_8/norm_l,response_vectors(1,iter+1), 1)
       ! v2 is normalized.

       IF (paral%parent) THEN
          IF (paral%io_parent)&
               WRITE(6,*)' *=*   overlap: alpha[',&
               iter,'] = ',alpha(iter)
          IF (paral%io_parent)&
               WRITE(6,*)' *=*   off-diag: beta[',&
               iter+1,'] = ',beta(iter+1)
          IF (paral%io_parent)&
               WRITE(6,*)' *=*   norm:               = ',norm_l
       ENDIF
       ! ==--------------------------------------------------------------==
       IF (paral%io_parent)&
            WRITE(6,*)'##############################',&
            '##################################'
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
       lastconverged = iter
       DO i=iter,1,-1
          IF ( ABS(eigenvect(iter,i)) .LT. eig2%conv_threshold) THEN
             converged = '       (converged:'
             IF ( lastconverged .EQ. i+1)lastconverged = i
          ELSE
             converged = '   (NOT converged:'
          ENDIF
          CALL xstring(converged,ia,ie)
          IF (paral%io_parent) WRITE(6,'(A,I2,A,F12.7,A,F8.6,A)')&
               ' *** eigenvalue ',i,' = ',diagonals(i)-eig2%eigenvalue_shift,&
               converged(ia:ie),ABS(eigenvect(iter,i)),').'
       ENDDO
       IF (paral%io_parent) THEN
          IF (paral%io_parent)&
               WRITE(6,*)'Last cnverged state:',lastconverged
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

    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    DO istate=lastvalid,1, -1
       ! build the eigenvector of the linear answer machinery:
       CALL dgemv('N', 2*ncpw%nhg, lastvalid, 1._real_8, response_vectors, 2*ncpw%nhg,&
            eigenvect(1,istate),1, 0._real_8, tmp, 1)


       IF (paral%io_parent) THEN
          DO i=1,lastvalid
             WRITE(6,'(A,I2,A,I2,A,F15.11)')&
                  '  \*  Lanczos(',i,'; ',istate,') = ',&
                  eigenvect(i,istate)
          ENDDO
       ENDIF

       ! write the output file...
       filename='eigst-'//zahlenchars(istate)
       CALL xstring(filename,ia,ie)
       IF (paral%io_parent)&
            WRITE(6,*)' _i_ Writing eigenstate ',istate,&
            ' to file ',filename(ia:ie),'.'
       CALL  densto(tmp,tau0,filename)
    ENDDO
    ! ==--------------------------------------------------------------==
    ! == END / Output of eigensolutions.
    ! ==--------------------------------------------------------------==


    DEALLOCATE(tmp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
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
  END SUBROUTINE eigensystem_p
  ! ==--------------------------------------------------------------==







END MODULE eigensystem_p_utils
