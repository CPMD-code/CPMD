MODULE forcematch_kfit_utils
  USE error_handling,                  ONLY: stopgm
  USE fm_cnst,                         ONLY: au_kjm,&
       bohr_nm,&
       kb_au2grm,&
       kb_grm2au,&
       kjm_au,&
       nm_bohr
  USE forcematch,                      ONLY: &
       fm_amber, fm_angopt, fm_b0, fm_bonopt, fm_dihopt, fm_fit_angles, &
       fm_fit_bonds, fm_fit_diheds, fm_fit_fc_only, fm_fit_improps, &
       fm_fqm_ref, fm_ib, fm_impopt, fm_ip, fm_iq, fm_it, fm_kang, fm_kb, &
       fm_kbon, fm_kdih, fm_kimp, fm_kp, fm_kq, fm_kt, &
       fm_max_kfit_iterations, fm_nang, fm_nangopt, fm_nbon, fm_nbonopt, &
       fm_ndih, fm_ndihopt, fm_nframes, fm_nimp, fm_nimpopt, fm_nrqm, fm_p0, &
       fm_pn, fm_q0, fm_rqm, fm_t0, grm_equiv
  USE forcematch_utils,                ONLY: fm_writeout_topo
  USE kinds,                           ONLY: real_8
  USE mm_dimmod,                       ONLY: nat_cpmd
  USE mm_input,                        ONLY: agr
  USE parac,                           ONLY: paral
  USE zeroing_utils,                   ONLY: zeroing
#if defined (__GROMOS)
  USE coordsz
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: fm_fit_covalent
  PUBLIC :: fm_setup_topology
  PUBLIC :: fm_kfit_compute_frms
  PUBLIC :: fm_kfit_rms_per_atom

CONTAINS

  ! Author: Patrik Mauer
  ! Revision: Ivano Tavernelli (Dec. 2010)
  ! ==================================================================
  SUBROUTINE fm_fit_covalent(nframes, nr_qm, r_qm, fqm_ref)
    ! ==--------------------------------------------------------------==
    ! The main driver routine for fitting covalent parameters (force
    ! constants, bond lengths, angles, etc.)
    ! NB: We use atomic units!
    ! ==================================================================
    ! nr of frames, nr of QM atoms
    INTEGER                                  :: nframes, nr_qm
    REAL(real_8)                             :: r_qm(3, nr_qm, nframes), &
                                                fqm_ref(3, nr_qm, nframes)

    CHARACTER(*), PARAMETER                  :: procedureN = 'fm_fit_covalent'

    INTEGER                                  :: ierr, npar
    REAL(real_8), ALLOCATABLE                :: par(:)

! coordinates of QM atoms
! reference forces on QM atoms
! Variables
! total nr. of parameters
! 
! 

    IF (paral%io_parent)&
         WRITE(6,*)
    IF (paral%io_parent)&
         WRITE(6,'(A)') '  fm  ========================================'
    IF (paral%io_parent)&
         WRITE(6,'(A)') '  fm  Now fitting covalent parameters'
    IF (paral%io_parent)&
         WRITE(6,'(A)') '  fm  ========================================'

    ! in case I want to exclude a subset of parameters from
    ! the fitting I adjust the number of parameters here:

    IF (.NOT. fm_fit_bonds) THEN
       fm_nbonopt=0
       IF (paral%io_parent)&
            WRITE(6,*) ' fm  readjusting number of fitting parms'
       IF (paral%io_parent)&
            WRITE(6,*) ' fm  nr of bonds to optimize: ',fm_nbonopt
    ENDIF
    IF (.NOT. fm_fit_angles) THEN
       fm_nangopt=0
       IF (paral%io_parent)&
            WRITE(6,*) ' fm  readjusting number of fitting parms'
       IF (paral%io_parent)&
            WRITE(6,*) ' fm  nr of angles to optimize: ',fm_nangopt
    ENDIF
    IF (.NOT. fm_fit_diheds) THEN
       fm_ndihopt=0
       IF (paral%io_parent)&
            WRITE(6,*) ' fm  readjusting number of fitting parms'
       IF (paral%io_parent)&
            WRITE(6,*) ' fm  nr of dihedrals to optimize: ',fm_ndihopt
    ENDIF
    IF (.NOT. fm_fit_improps) THEN
       fm_nimpopt=0
       IF (paral%io_parent)&
            WRITE(6,*) ' fm  readjusting number of fitting parms'
       IF (paral%io_parent)&
            WRITE(6,*) ' fm  nr of improper dihedrals to optimize: ',&
            fm_nimpopt
    ENDIF

    ! compute total nr. of parameters. For impropers and dihedrals, we
    ! always fit force constants only. For bonds and angles, we can also fit
    ! equilibrium values
    ! ordering of parameters is:
    ! - with only force constants:
    ! kb(1..N), kt(1..N), kq(1..N), kp(1..N)
    ! - with equilibrium values:
    ! kb(1),b0(1),...,kb(NB),b0(NB), kt(1),t0(1),...,kt(N), t0(N),
    ! kq(1..N), kp(1..N)
    IF (fm_fit_fc_only) THEN
       npar = fm_nbonopt + fm_nangopt + fm_nimpopt + fm_ndihopt
    ELSE
       npar = 2*fm_nbonopt + 2*fm_nangopt + fm_nimpopt + fm_ndihopt
    ENDIF
    IF (paral%io_parent)&
         WRITE(6,'(A, i6)') '  fm  Total nr. of parameters: ', npar

    ALLOCATE(par(npar),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! We need to set some common variables to circumvent the optimizer
    fm_nframes = nframes
    fm_nrqm = nr_qm

    CALL fm_kfit_collect_optimized
    CALL fm_kfit_setup_mytypes

    ! _DEBUG[
    CALL fm_writeout_topo('beforefit.top')
    ! _DEBUG]
    IF (fm_fit_fc_only) THEN
       IF (paral%io_parent)&
            WRITE(6,*) ' fm  fitting force constants only'

       ! TESTING ONLY FIT ANGLES
       ! fm_nbonopt =0
       ! fm_nangopt = 0
       ! fm_nimpopt=0
       ! fm_ndihopt=0
       ! npar = fm_ndihopt
       ! write(6,*) 'total nr. of parameters: ',npar
       CALL fm_optimize_fc (npar,par,nframes,nr_qm,r_qm,fqm_ref)
    ELSE
       CALL fm_initialguess(npar, par)
       CALL fm_optimize_allcovalent(npar, par, nframes, nr_qm, r_qm,&
            fqm_ref)
    ENDIF
    CALL fm_kfit_par_to_fmarray(npar, par)
    CALL fm_kfit_writeback(npar, par)

    DEALLOCATE(par,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    RETURN
  END SUBROUTINE fm_fit_covalent
  ! 
  ! ==================================================================
  SUBROUTINE fm_optimize_fc (npar,par,nframes,nr_qm,r_qm,fqm_ref)
    ! ==--------------------------------------------------------------==
    ! Driver for the optimization of the force constants only in the 
    ! covalent interactions (including bond/angles/torsion and improper 
    ! torsions potentials). This is a linear set of equations and we use 
    ! the routine dgels from LAPACK to solve it.
    ! ==================================================================
    ! nr of parameters, nr. of frames, nr. of QM atoms      
    INTEGER                                  :: npar
    REAL(real_8)                             :: par(npar)
    INTEGER                                  :: nframes, nr_qm
    REAL(real_8)                             :: r_qm(3, nr_qm, nframes), &
                                                fqm_ref(3, nr_qm, nframes)

    INTEGER :: at, bt, dt, i, iang, ibon, idih, iframe, iimp, indx, info, &
      ipar, ipar0, it, j, k, l, lb, lwoa, no_coloumns, no_rows, nrhs, ub
    REAL(real_8) :: bla(3*nr_qm*nframes), blu(3*nr_qm*nframes,npar), fi(3), &
      fj(3), fk(3), fl(3), fqm_ref_fc(3*nr_qm*nframes), ri(3), rj(3), rk(3), &
      rl(3), wa(3*npar)
    REAL(real_8), POINTER                    :: derivs(:, :)
    REAL(real_8), TARGET :: fm_fc_derivs(3*nr_qm*nframes,npar)

! variables for the loops      
! dummy integer
! parameters, the array that's going to be updated
! coordinates of QM atoms
! reference forces
! derivatives coming back from fm_forces_*** routines
! the derivative matrix
! number of rows: 3*noatoms*noframes
! number of columns:  noof parameters
! debugging
! lapack related stuff:
! pass some information codes around
! work array, and size of work array

    lwoa=3*npar
    ! nullify the derivs matrix
    DO i=1,3*nr_qm*nframes
       DO j=1,npar
          fm_fc_derivs(i,j)=0.0_real_8
       ENDDO
    ENDDO
    ! I need a copy of the reference forces array
    ! because it is going to be overwritten by the lapack routine
    indx=0
    DO i=1,nframes
       DO j=1,nr_qm
          DO k=1,3
             indx=indx+1
             fqm_ref_fc(indx)=fqm_ref(k,j,i)
             ! write(6,*) 'fqm_ref_fc ',fqm_ref_fc(indx)
          ENDDO
       ENDDO
    ENDDO
    ! TESTING ONLY FIT BONDS
    ! fill the derivative matrix
    ! Loop over all frames. In each frame, we use derivs to point to
    ! the appropriate "windows" of the  full matrix fm_fc_derivs,
    ! to facilitate the indexing business
    ! write(6,*) fm_nframes,fm_nrqm,lb,ub,fm_nbon
    DO iframe = 1, fm_nframes
       lb = (iframe-1) * 3*fm_nrqm + 1
       ub = lb + 3*fm_nrqm - 1
       derivs => fm_fc_derivs(lb:ub, :)
       ! write(6,*) 'fm_fc_derivs() ',lb,' : ',ub
       ! derivs contains now one frame for the matrix fm_fc_derivs          
       ! Bonds   (ipar0=0: bonds come first)
       ipar0 = 0
       DO ibon = 1, fm_nbon
          ! write(6,*) fm_ib(:,ibon)
          ! write(6,*) 'grumpf',fm_ib(0, ibon)
          bt = fm_ib(0, ibon)
          IF (paral%io_parent)&
               WRITE(6,*)&
               'ibon',ibon,' bt ',bt,' kb ',fm_kb(bt),&
               ' gromos: ',kb_au2grm*fm_kb(bt)
          ! write(6,*) 'grumpf',fm_ib(1, ibon)
          i = fm_ib(1, ibon)
          ! write(6,*) i
          ! write(6,*) 'grumpf',fm_ib(2, ibon)
          j = fm_ib(2, ibon)
          ! write(6,*) j
          ! write(6,*) 'kampf',r_qm(:, i, iframe)
          ri = r_qm(:, i, iframe)
          ! write(6,*) 'ri ',ri
          rj = r_qm(:, j, iframe)
          IF (paral%io_parent)&
               WRITE(6,*) 'i ',i,' j ',j,' ri ',ri,' rj ',rj

          ! compute i,j as index into array derivs
          i = 3*(i-1) + 1
          j = 3*(j-1) + 1
          CALL fm_forces_bon(fi, fj, ri, rj, bt, 'dk')
          ! the column index             
          ipar = ipar0 + bt
          derivs(i:i+2, ipar) = derivs(i:i+2, ipar) + fi
          derivs(j:j+2, ipar) = derivs(j:j+2, ipar) + fj
          ! write(6,*) 'derivs i ',i,' : ',i+2
          ! write(6,*) derivs(i:i+2, ipar)
          ! write(6,*) 'derivs j ',j,' : ',j+2
          ! write(6,*) derivs(j:j+2, ipar)
          ! write(6,*) 'grimpf'
       ENDDO! loop over all the bonds in the QM system
       ! Angles
       ipar0 = fm_nbonopt
       DO iang = 1, fm_nang
          at = fm_it(0, iang)
          i = fm_it(1, iang)
          j = fm_it(2, iang)
          k = fm_it(3, iang)
          ri = r_qm(:, i, iframe)
          rj = r_qm(:, j, iframe)
          rk = r_qm(:, k, iframe)

          i = 3*(i-1) + 1
          j = 3*(j-1) + 1
          k = 3*(k-1) + 1

          ipar = ipar0 + at
          CALL fm_forces_ang(fi, fj, fk, ri, rj, rk, at, 'dk')
          derivs(i:i+2, ipar) = derivs(i:i+2, ipar) + fi
          derivs(j:j+2, ipar) = derivs(j:j+2, ipar) + fj
          derivs(k:k+2, ipar) = derivs(k:k+2, ipar) + fk
       ENDDO! loop over angle types

       ! Impropers
       ipar0 = ipar0 + fm_nangopt
       DO iimp = 1, fm_nimp
          it = fm_iq(0, iimp)
          i = fm_iq(1, iimp)
          j = fm_iq(2, iimp)
          k = fm_iq(3, iimp)
          l = fm_iq(4, iimp)
          ri = r_qm(:, i, iframe)
          rj = r_qm(:, j, iframe)
          rk = r_qm(:, k, iframe)
          rl = r_qm(:, l, iframe)

          i = 3*(i-1) + 1
          j = 3*(j-1) + 1
          k = 3*(k-1) + 1
          l = 3*(l-1) + 1
          ipar = ipar0 + it
          CALL fm_forces_imp(fi,fj,fk,fl,ri,rj,rk,rl,it,'dk')
          derivs(i:i+2, ipar) = derivs(i:i+2, ipar) + fi
          derivs(j:j+2, ipar) = derivs(j:j+2, ipar) + fj
          derivs(k:k+2, ipar) = derivs(k:k+2, ipar) + fk
          derivs(l:l+2, ipar) = derivs(l:l+2, ipar) + fl
       ENDDO! loop over all the improper dihedrals
       ! Dihedrals
       ipar0 = ipar0 + fm_nimpopt
       DO idih = 1, fm_ndih
          dt = fm_ip(0, idih)
          i = fm_ip(1, idih)
          j = fm_ip(2, idih)
          k = fm_ip(3, idih)
          l = fm_ip(4, idih)
          ri = r_qm(:, i, iframe)
          rj = r_qm(:, j, iframe)
          rk = r_qm(:, k, iframe)
          rl = r_qm(:, l, iframe)

          i = 3*(i-1) + 1
          j = 3*(j-1) + 1
          k = 3*(k-1) + 1
          l = 3*(l-1) + 1
          ipar = ipar0 + dt
          CALL fm_forces_dih(fi,fj,fk,fl,ri,rj,rk,rl,dt,'dk')
          derivs(i:i+2, ipar) = derivs(i:i+2, ipar) + fi
          derivs(j:j+2, ipar) = derivs(j:j+2, ipar) + fj
          derivs(k:k+2, ipar) = derivs(k:k+2, ipar) + fk
          derivs(l:l+2, ipar) = derivs(l:l+2, ipar) + fl
       ENDDO! loop over all the dihedrals



    ENDDO! loop over no frames
    ! solve the matrix equation: fm_fc_derivs * X = fqm_ref_fc          

    DO i=1,3*nr_qm*nframes
       DO j=1,npar
          blu(i,j)=fm_fc_derivs(i,j)
       ENDDO
       IF (paral%io_parent)&
            WRITE(6,*) fm_fc_derivs(i,:)
    ENDDO
    ! write(6,*) 'parameters beginning'
    ! do i=1,fm_nbonopt
    ! par(i)=fm_kb(i)
    ! write(6,*) par(i),par(i)*kb_au2grm
    ! enddo   

    ! bla=matmul(blu,par)
    ! write(6,*) 'first test'
    ! indx=0
    ! do i=1,nframes
    ! do j=1,nr_qm
    ! do k=1,3
    ! indx=indx+1
    ! write(6,*) bla(indx),fqm_ref(k,j,i)
    ! enddo
    ! enddo
    ! enddo   
    ! write(6,*) 'end first test'
    bla(:)=0.0_real_8
    par(:)=0.0_real_8

    no_rows=3*nr_qm*nframes
    no_coloumns=npar
    nrhs=1
    CALL dgels('N',no_rows,no_coloumns,nrhs,fm_fc_derivs,&
         no_rows,fqm_ref_fc,no_rows,wa,lwoa,info)    
    ! call dgels('N',3*nr_qm*nframes,npar,1,fm_fc_derivs,
    ! &       3*nr_qm*nframes,fqm_ref_fc,3*nr_qm*nframes,wa,lwoa,info)     
    ! note that both matrices fm_fc_derivs and fqm_ref_fc have been 
    ! overwritten
    ! the npar first elements of fqm_ref_fc is the solution 
    ! write(6,*) ' fm  after dgels, info=',info
    IF (info .NE. 0) THEN
       IF (paral%io_parent)&
            WRITE(6,'(A)') '  fm  Optimization FAILED:'
       CALL stopgm('fm_optimize_fc', 'FAIL (see output)',& 
            __LINE__,__FILE__)
    ELSE
       IF (paral%io_parent)&
            WRITE(6,'(A)') '  fm  Optimization lyn sys successfull'
    ENDIF
    ! write(6,*) ' fm  list of parameters after fitting: au, grm'
    ! do i=1,npar
    ! write(6,*) fqm_ref_fc(i),kb_au2grm* fqm_ref_fc(i)
    ! enddo 

    ! write back to the par array
    DO i=1,npar
       par(i)=fqm_ref_fc(i)
    ENDDO
    ! I should get back the reference forces!
    IF (paral%io_parent)&
         WRITE(6,*) ' '
    IF (paral%io_parent)&
         WRITE(6,*) ' fm  test: forces after fit, ref forces'
    bla=MATMUL(blu,par)
    indx=0
    DO i=1,nframes
       DO j=1,nr_qm
          DO k=1,3
             indx=indx+1
             IF (paral%io_parent)&
                  WRITE(6,*) bla(indx),fqm_ref(k,j,i)
          ENDDO
       ENDDO
    ENDDO
    ! do i=1,3*nr_qm*nframes
    ! write(6,*) bla(i),fqm_ref_fc(i)
    ! enddo

    ! call STOPGM('fm_optimize_fc', 'TEST')

    RETURN
  END SUBROUTINE fm_optimize_fc
  ! 
  ! ==================================================================
  SUBROUTINE fm_optimize_allcovalent(npar, par, nframes, nr_qm,&
       r_qm, fqm_ref)
    ! ==--------------------------------------------------------------==
    ! Driver for the optimization of ALL covalent interactions parameters
    ! (including bond lengths/angles). This is a non-linear problem and 
    ! we use the routine lmder1 from MINPACK to solve it.
    ! ==================================================================
    ! nr. of parameters, nr. of frames, nr. of QM atoms
    INTEGER                                  :: npar
    REAL(real_8)                             :: par(npar)
    INTEGER                                  :: nframes, nr_qm
    REAL(real_8), TARGET                     :: r_qm(3, nr_qm, nframes), &
                                                fqm_ref(3, nr_qm, nframes)

    CHARACTER(len=255)                       :: err_msg, succ_msg
    INTEGER                                  :: info, lwa, ipvt(npar), &
                                                trget_size
    LOGICAL                                  :: success
    REAL(real_8) :: fdiff(3, nr_qm, nframes), fjac(3*nr_qm*nframes, npar), &
      tol, work_array(5*npar*3*nr_qm*nframes)

! coordinates of QM atoms
! reference forces
! Variables
! size of the target vector (=3*nr_qm*nframes)
! vector with force differences 
! jacobian matrix (derivative of fdiff wrt. parameters)
! was the optimization successful?
! strings for error or success messages
! minpack related stuff:
! tolerance (convergence criterion)
! pass some information codes around
! pivot array, work array, and size of work array
! real(8), allocatable    :: work_array(:)
! name of the function that evaluates fdiff and fjac
! external fm_fcn_allcovalent
! hardcode tol for now
! FIXME: make tol user-definable

    tol = 1.e-06_real_8

    fm_rqm => r_qm
    fm_fqm_ref => fqm_ref

    trget_size = 3*nr_qm*nframes

    lwa = 5 * npar * trget_size
    ! allocate(wa(lwa))
    ! _DEBUG[
    ! write(6,*) ' lwa is: ', lwa
    ! write(6,*) ' allocated wa?', allocated(work_array)
    ! write(6,*) ' allocating wa'
    ! _DEBUG]
    ! allocate(work_array(lwa))
    ! _DEBUG[
    ! write(6,*) ' allocated wa?', allocated(work_array)
    ! write(6,*) ' size wa?', size(work_array)
    ! _DEBUG]
    success = .FALSE.
#ifdef __MINPACK
    CALL lmder1(fm_fcn_allcovalent, trget_size, npar, par, fdiff,&
         fjac, trget_size, tol, info, ipvt, work_array, lwa)
#else
    info = 0
    CALL stopgm('fm_optimize_allcovalent', 'No MINPACK linked! ',& 
         __LINE__,__FILE__)
#endif

    ! _DEBUG[
    IF (paral%io_parent)&
         WRITE(6,*) " lmder1 returned info = ", info
    ! _DEBUG]
    ! check if we were successful
    IF (info.LT.0) THEN
       err_msg = 'exceeded max. nr of iterations'
       success = .FALSE.
    ELSEIF (info.EQ.0) THEN
       success = .FALSE.
       err_msg = 'improper input parameters to lmder1'
    ELSEIF (info.EQ.1) THEN
       success = .TRUE.
       succ_msg = 'algorithm estimates that the relative error in'&
            //'the sum of squares is at most tol'
    ELSEIF (info.EQ.2) THEN
       success = .TRUE.
       succ_msg = 'algorithm estimates that the relative error'&
            //'between x and the solution is at most tol'
    ELSEIF (info.EQ.3) THEN
       success = .TRUE.
       succ_msg = 'super converged! '
    ELSEIF (info.EQ.4) THEN
       success = .TRUE.! I guess???
       succ_msg = 'fvec is orthogonal to the columns of the jacobian'&
            //' to machine precision'
    ELSEIF (info.EQ.5) THEN
       success = .FALSE.
       err_msg = 'too many iterations needed'
    ELSEIF (info.EQ.6) THEN
       success = .TRUE.! sure?
       succ_msg = 'tol is too small. no further reduction in the'&
            //'sum of squares is possible'
    ELSEIF (info.EQ.7) THEN
       success = .TRUE.! sure?
       succ_msg = 'tol is too small. no further improvement in the'&
            //'approximate solution x is possible'
    ELSE
       success = .FALSE.
       err_msg = ' I do not know what this value of info means'
    ENDIF

    IF (.NOT.success) THEN
       IF (paral%io_parent)&
            WRITE(6,'(A)') '  fm  Optimization FAILED:'
       IF (paral%io_parent)&
            WRITE(6,'(A)') '  fm  ' // err_msg
       CALL stopgm('fm_optimize_allcovalent', 'FAIL (see output)',& 
            __LINE__,__FILE__)
    ENDIF
    IF (paral%io_parent)&
         WRITE(6,'(A)') '  fm  Optimization successful:'
    IF (paral%io_parent)&
         WRITE(6,'(A)') '  fm  ' // succ_msg

    ! _DEBUG[
    ! write(6,*) 'deallocating wa. lwa is: ', lwa
    ! write(6,*) ' allocated(wa)?:', allocated(work_array)
    ! write(6,*) ' size(wa):', size(work_array)
    ! write(6,*) ' first 3 elements:', work_array(1:3)
    ! _DEBUG]
    ! deallocate(work_array)
    ! _DEBUG[
    ! write(6,*) 'done deallocating wa. '
    ! call m_flush(6)
    ! _DEBUG]
    RETURN
  END SUBROUTINE fm_optimize_allcovalent
  ! 
  ! ==================================================================
  SUBROUTINE fm_kfit_collect_optimized
    ! ==--------------------------------------------------------------==
    ! collect all bonds (and angles, impropers, dihedrals) that are 
    ! optimized.
    ! ==================================================================
#if defined (__GROMOS)
    IMPLICIT NONE
    INCLUDE 'gromos.h'
    ! Variables
    ! loop counters, temporaries, etc.
    INTEGER :: i, j, ii,ierr
    CHARACTER(*),PARAMETER::procedureN='fm_kfit_collect_optimized'
    IF (paral%io_parent)&
         WRITE(6,'(A)') '  fm  collecting covalent interactions that are'
    IF (paral%io_parent)&
         WRITE(6,'(A)') '  fm  affected by optimization'
    ! bonds
    ! allocate(fm_ib(0:2, MAXBNH+MAXBON))
    ALLOCATE(fm_ib(0:2, nbonh + nbon),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ii = 0
    DO i = 1, nbonh
       DO j = 1, fm_nbonopt
          ! if(ICBH(i).eq.j) then
          IF (icbh(i).EQ.fm_bonopt(j)) THEN
             ii = ii + 1
             fm_ib(0, ii) = j
             fm_ib(1, ii) = NAT_cpmd(ibh(i))
             fm_ib(2, ii) = NAT_cpmd(jbh(i))
             ! _DEBUG[
             ! write(6,*) " deb: bond i, ICBH(i), fm_bond i"
             ! write(6,*) i, ICBH(i), ii
             ! write(6,*) "      fm_ib(0-2, ii)"
             ! write(6,*) fm_ib(:, ii)
             ! _DEBUG]
          ENDIF
       ENDDO
    ENDDO
    DO i = 1, nbon
       DO j = 1, fm_nbonopt
          ! if(ICB(i).eq.j) then
          IF (icb(i).EQ.fm_bonopt(j)) THEN
             ii = ii + 1
             fm_ib(0, ii) = j
             fm_ib(1, ii) = NAT_cpmd(ib(i))
             fm_ib(2, ii) = NAT_cpmd(jb(i))
             ! _DEBUG[
             ! write(6,*) " deb: bond i, ICB(i), fm_bond i"
             ! write(6,*) i, ICB(i), ii
             ! write(6,*) "fm_bonopt(j):", fm_bonopt(j)
             ! write(6,*) "      fm_ib(0-2, ii)"
             ! write(6,*) fm_ib(:, ii)
             ! _DEBUG]
          ENDIF
       ENDDO
    ENDDO
    fm_nbon = ii
    IF (paral%io_parent)&
         WRITE(6,'(A,i6)') '  fm  bonds:     ', fm_nbon
    ! angles
    ALLOCATE(fm_it(0:3, mxqheh+maxthe),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ii = 0
    DO i = 1, ntheh
       DO j = 1, fm_nangopt
          ! if(ICTH(i).eq.j) then
          IF (icth(i).EQ.fm_angopt(j)) THEN
             ii = ii + 1
             fm_it(0, ii) = j
             fm_it(1, ii) = NAT_cpmd(ith(i))
             fm_it(2, ii) = NAT_cpmd(jth(i))
             fm_it(3, ii) = NAT_cpmd(kth(i))
          ENDIF
       ENDDO
    ENDDO
    DO i = 1, nthe
       DO j = 1, fm_nangopt
          ! if(ICT(i).eq.j) then
          IF (ict(i).EQ.fm_angopt(j)) THEN
             ii = ii + 1
             fm_it(0, ii) = j
             fm_it(1, ii) = NAT_cpmd(it(i))
             fm_it(2, ii) = NAT_cpmd(jt(i))
             fm_it(3, ii) = NAT_cpmd(kt(i))
          ENDIF
       ENDDO
    ENDDO
    fm_nang = ii
    IF (paral%io_parent)&
         WRITE(6,'(A,i6)') '  fm  angles:    ', fm_nang
    ! impropers
    ALLOCATE(fm_iq(0:4, maxhih+maxqhi),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ii = 0
    DO i = 1, nqhih
       DO j = 1, fm_nimpopt
          ! if(ICQH(i).eq.j) then
          IF (icqh(i).EQ.fm_impopt(j)) THEN
             ii = ii + 1
             fm_iq(0, ii) = j
             fm_iq(1, ii) = NAT_cpmd(iqh(i))
             fm_iq(2, ii) = NAT_cpmd(jqh(i))
             fm_iq(3, ii) = NAT_cpmd(kqh(i))
             fm_iq(4, ii) = NAT_cpmd(lqh(i))
          ENDIF
       ENDDO
    ENDDO
    DO i = 1, nqhi
       DO j = 1, fm_nimpopt
          ! if(ICQ(i).eq.j) then
          IF (icq(i).EQ.fm_impopt(j)) THEN
             ii = ii + 1
             fm_iq(0, ii) = j
             fm_iq(1, ii) = NAT_cpmd(iq(i))
             fm_iq(2, ii) = NAT_cpmd(jq(i))
             fm_iq(3, ii) = NAT_cpmd(kq(i))
             fm_iq(4, ii) = NAT_cpmd(lq(i))
          ENDIF
       ENDDO
    ENDDO
    fm_nimp = ii
    IF (paral%io_parent)&
         WRITE(6,'(A,i6)') '  fm  impropers: ', fm_nimp
    ! dihedrals
    ALLOCATE(fm_ip(0:4, mxphih+maxphi),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ii = 0
    DO i = 1, nphih
       DO j = 1, fm_ndihopt
          ! if(ICPH(i).eq.j) then
          IF (icph(i).EQ.fm_dihopt(j)) THEN
             ii = ii + 1
             fm_ip(0, ii) = j
             fm_ip(1, ii) = NAT_cpmd(iph(i))
             fm_ip(2, ii) = NAT_cpmd(jph(i))
             fm_ip(3, ii) = NAT_cpmd(kph(i))
             fm_ip(4, ii) = NAT_cpmd(lph(i))
          ENDIF
       ENDDO
    ENDDO
    DO i = 1, nphi
       DO j = 1, fm_ndihopt
          ! if(ICP(i).eq.j) then
          IF (icp(i).EQ.fm_dihopt(j)) THEN
             ii = ii + 1
             fm_ip(0, ii) = j
             fm_ip(1, ii) = NAT_cpmd(ip(i))
             fm_ip(2, ii) = NAT_cpmd(jp(i))
             fm_ip(3, ii) = NAT_cpmd(kp(i))
             fm_ip(4, ii) = NAT_cpmd(lp(i))
          ENDIF
       ENDDO
    ENDDO
    fm_ndih = ii
    IF (paral%io_parent)&
         WRITE(6,'(A,i6)') '  fm  dihedrals: ', fm_ndih

    RETURN
#endif
  END SUBROUTINE fm_kfit_collect_optimized
  ! 
  ! ==================================================================
  SUBROUTINE fm_initialguess(npar, par)
    ! ==--------------------------------------------------------------==
    ! provide an initial guess from values in the original topology
    ! ==================================================================
#if defined (__GROMOS)
#endif
    IMPLICIT NONE
    ! Arguments
    ! total nr. of parameters
    INTEGER :: npar
    ! parameters
    REAL(real_8) :: par(npar)
#if defined (__GROMOS)
    ! Variables
    ! loop counters, temporaries, etc.
    INTEGER :: iframe, iqm, ipar, jpar
    IF (paral%io_parent)&
         WRITE(6,'(A)') '  fm  Getting initial guess'
    jpar = 0
    DO ipar = 1, fm_nbonopt
       jpar = jpar + 1
       par(jpar) = fm_kb(ipar)
       jpar = jpar + 1
       par(jpar) = fm_b0(ipar)
    ENDDO
    DO ipar = 1, fm_nangopt
       jpar = jpar + 1
       par(jpar) = fm_kt(ipar)
       jpar = jpar + 1
       par(jpar) = fm_t0(ipar)
    ENDDO
    DO ipar = 1, fm_nimpopt
       jpar = jpar + 1
       par(jpar) = fm_kq(ipar)
    ENDDO
    DO ipar = 1, fm_ndihopt
       jpar = jpar + 1
       par(jpar) = fm_kp(ipar)
    ENDDO
    IF (jpar .NE. npar) CALL stopgm('fm_initialguess',&
         'jpar .ne. npar',& 
         __LINE__,__FILE__)

    RETURN
#endif
  END SUBROUTINE fm_initialguess
  ! 
  ! ==================================================================
  SUBROUTINE fm_kfit_setup_mytypes
    ! ==--------------------------------------------------------------==
#if defined (__GROMOS)
    IMPLICIT NONE
    INCLUDE 'gromos.h'
    ! Variables
    ! loop counters
    INTEGER :: i, j,ierr
    CHARACTER(*),PARAMETER::procedureN='fm_kfit_setup_mytypes'
    IF (paral%io_parent)&
         WRITE(6,'(A)') '  fm  Setting up'
    ! copy qmmm_amber
    fm_amber = agr%qmmm_amber

    ! Bonds
    ALLOCATE(fm_kb(fm_nbonopt),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(fm_b0(fm_nbonopt),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    DO i = 1, fm_nbonopt
       ! put back saved value
       cb(fm_bonopt(i)) = fm_kbon(i)
       ! convert parameters to a.u.
       fm_b0(i) = nm_bohr * b0(fm_bonopt(i))
       fm_kb(i) = kb_grm2au * cb(fm_bonopt(i))
       IF (fm_amber) THEN
          fm_kb(i) = 2._real_8 * fm_b0(i)**2 * fm_kb(i)
       ENDIF
    ENDDO

    ! Angles
    ALLOCATE(fm_kt(fm_nangopt),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(fm_t0(fm_nangopt),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    DO i = 1, fm_nangopt
       ! put back saved value
       ct(fm_angopt(i)) = fm_kang(i)
       ! convert parameters to a.u.
       fm_t0(i) = t0(fm_angopt(i))
       IF (fm_amber) THEN
          ! We use the simplified conversion formula given in the GROMACS manual
          ! (as opposed to the temperature-dependent one in the GROMOS manual)
          ! according to the GROMACS manual, the error at 300K is ~0.2 %
          fm_kt(i) = kjm_au * SIN(dacos(t0(fm_angopt(i))))**2&
               * CT(fm_angopt(i))
       ELSE
          fm_kt(i) = kjm_au * ct(fm_angopt(i))
       ENDIF
    ENDDO

    ! Impropers
    ALLOCATE(fm_kq(fm_nimpopt),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(fm_q0(fm_nimpopt),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    DO i = 1, fm_nimpopt
       ! put back saved value
       cq(fm_impopt(i)) = fm_kimp(i)
       ! convert parameters to a.u.
       fm_q0(i) = q0(fm_impopt(i))
       fm_kq(i) = kjm_au * cq(fm_impopt(i))
    ENDDO

    ! Dihedrals
    ALLOCATE(fm_kp(fm_ndihopt),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(fm_p0(fm_ndihopt),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(fm_pn(fm_ndihopt),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    DO i = 1, fm_ndihopt
       ! put back saved value
       cp(fm_dihopt(i)) = fm_kdih(i)
       ! convert parameters to a.u.
       fm_p0(i) = pd(fm_dihopt(i))
       fm_pn(i) = REAL(np(fm_dihopt(i)),kind=real_8)
       fm_kp(i) = kjm_au * cp(fm_dihopt(i))
    ENDDO

    RETURN
#endif
  END SUBROUTINE fm_kfit_setup_mytypes
  ! 
  ! ==================================================================
  SUBROUTINE fm_setup_topology(nrpt, nsolv, is_qm)
    ! ==--------------------------------------------------------------==
    ! setup topology information for force-matching, .i.e. create new 
    ! types, etc.
    ! ==================================================================
    ! nr of solute atoms, solvent atoms
    INTEGER                                  :: nrpt, nsolv
    LOGICAL                                  :: is_qm(nrpt+nsolv)

! are atoms QM?

#if defined (__GROMOS)
    INCLUDE 'gromos.h'
    ! Variables
    ! function that checks if two interactions are between equivalent
    ! atoms
    ! we need a new type
    LOGICAL :: need_new
    ! array to flag type to be optimized
    LOGICAL, ALLOCATABLE :: optimize(:)
    ! keep some old values
    INTEGER :: old_type, old_ntypes
    ! iterators
    INTEGER :: i, j, k,ierr
    CHARACTER(*),PARAMETER::procedureN='fm_setup_topology'
    ! 
    IF (paral%io_parent)&
         WRITE(6,'(A)') '  fm '
    IF (paral%io_parent)&
         WRITE(6,'(A)') '  fm  Setting up topology for force-matching '&
         // 'and QM/MM'
    ! Initialize
    fm_nbonopt = 0
    fm_nangopt = 0
    fm_nimpopt = 0
    fm_ndihopt = 0
    ALLOCATE(fm_bonopt(maxnbt),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(fm_angopt(maxtty),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(fm_impopt(maxqty),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(fm_dihopt(maxpty),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(optimize(MAX(maxnbt,maxtty,maxqty,maxpty)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! create new types for interactions that are between QM atoms, but also
    ! appear between QM and MM or MM and MM atoms
    ! First come bonds involving hydrogens. How we do it: Loop over all
    ! bonds involving hydrogens (HBOND), and check if both atoms are QM. 
    ! If so, loop over all HBOND of same type and see if there is one that
    ! involves at least one MM atom. If that is the case, we need a new bond
    ! type for the QM only bond, so we create it as NBTY+1. Then loop again
    ! over the remaining HBOND of same (old) type to find all-QM bonds that
    ! are equivalent to the one just created, and change their bond type.
    ! 
    ! This procedure is then repeated for bonds not involving hydrogens, and
    ! analogously for angles, impropers, and dihedrals, each with and
    ! without H...
    ! That's what happens in the following ~300 lines of code...
    ! 
    ! TODO: For dihedrals (and maybe impropers), we should probably not
    ! create a new type for any non-equivalent quadruple of QM atoms.
    ! Instead, if a type involves both QM-QM and QM-MM/MM-MM interactions,
    ! simply duplicate it and fit to ONE type regardless of inequivalent QM
    ! atoms
    ! FIXME: If there is the same type between inequivalent QM atoms, but
    ! but type is not used by QM-MM or MM-MM interactions, then NO new type 
    ! is created, i.e. inequivalent QM atoms will interact via the same old
    ! type (for dihedrals, this is probably what we want anyway..)
    ! DONE: The above is fixed for bonds
    ! 
    ! -------bonds involving hydrogens
    optimize = .FALSE.
    old_ntypes = nbty
    DO i = 1, nbonh
       need_new = .FALSE.
       IF (is_qm(ibh(i)) .AND. is_qm(jbh(i))) THEN
          ! this bond is QM only. Now loop over other bonds of same type
          ! to see if this type occurs also as MM-QM or MM-MM
          DO j = 1, nbonh
             IF (j.EQ.i) CYCLE
             IF (icbh(i).EQ.icbh(j)) THEN
                IF (.NOT.is_qm(ibh(j)) .OR. .NOT.is_qm(jbh(j))) THEN
                   need_new = .TRUE.
                   EXIT
                ENDIF
             ENDIF
          ENDDO
          ! now make sure that we do not have a earlier QM-QM bond of
          ! this type with NON-equivalent atoms
          IF (.NOT.need_new .AND. optimize(icbh(i))) THEN
             DO j = 1, i-1
                IF (icbh(i).EQ.icbh(j)) THEN
                   IF (.NOT.are_equivalent('hbond',i,j)) THEN
                      need_new = .TRUE.
                      EXIT
                   ENDIF
                ENDIF
             ENDDO
          ENDIF
          IF (.NOT.need_new) THEN
             optimize(icbh(i)) = .TRUE.
             CYCLE
          ENDIF
          old_type = icbh(i)
          nbty = nbty + 1
          IF (nbty.GT.maxnbt) THEN
             CALL stopgm('fm_setup_topology','NBTY.gt.MAXNBT',& 
                  __LINE__,__FILE__)
          ENDIF
          icbh(i) = nbty
          cb(nbty) = cb(old_type)
          b0(nbty) = b0(old_type)
          ! check for equivalencies
          DO j = i+1, nbonh
             IF (icbh(j) .EQ. old_type) THEN
                IF (is_qm(ibh(j)).AND.is_qm(jbh(j))) THEN
                   IF (are_equivalent('hbond', i, j)) icbh(j) = nbty
                ENDIF
             ENDIF
          ENDDO
          optimize(icbh(i)) = .TRUE.
       ENDIF! both atoms are QM
    ENDDO
    IF (paral%io_parent)&
         WRITE(6, '(A, i6)') '  fm  New bond types for bonds with H: ',&
         NBTY - old_ntypes
    ! -------bonds not involving hydrogens
    old_ntypes = nbty
    DO i = 1, nbon
       need_new = .FALSE.
       IF (is_qm(ib(i)) .AND. is_qm(jb(i))) THEN
          ! this bond is QM only. Now loop over other bonds of same type
          ! to see if this type occurs also as MM-QM or MM-MM
          DO j = 1, nbon
             IF (j.EQ.i) CYCLE
             IF (icb(i).EQ.icb(j)) THEN
                IF (.NOT.is_qm(ib(j)) .OR. .NOT.is_qm(jb(j))) THEN
                   need_new = .TRUE.
                   EXIT
                ENDIF
             ENDIF
          ENDDO
          ! now make sure that we do not have a earlier QM-QM bond of
          ! this type with NON-equivalent atoms
          IF (.NOT.need_new .AND. optimize(icb(i))) THEN
             DO j = 1, i-1
                IF (icb(i).EQ.icb(j)) THEN
                   IF (.NOT.are_equivalent('bond',i,j)) THEN
                      need_new = .TRUE.
                      EXIT
                   ENDIF
                ENDIF
             ENDDO
          ENDIF
          IF (.NOT.need_new) THEN
             optimize(icb(i)) = .TRUE.
             CYCLE
          ENDIF
          old_type = icb(i)
          nbty = nbty + 1
          IF (nbty.GT.maxnbt) THEN
             CALL stopgm('fm_setup_topology','NBTY.gt.MAXNBT',& 
                  __LINE__,__FILE__)
          ENDIF
          icb(i) = nbty
          cb(nbty) = cb(old_type)
          b0(nbty) = b0(old_type)
          ! check for equivalencies
          DO j = i+1, nbon
             IF (icb(j) .EQ. old_type) THEN
                IF (is_qm(ib(j)).AND.is_qm(jb(j))) THEN
                   IF (are_equivalent('bond', i, j)) icb(j) = nbty
                ENDIF
             ENDIF
          ENDDO
          optimize(icb(i)) = .TRUE.
       ENDIF
    ENDDO
    IF (paral%io_parent)&
         WRITE(6, '(A, i6)') '  fm  New bond types for other bonds:  ',&
         NBTY - old_ntypes
    fm_nbonopt = 0
    DO i = 1, nbty
       IF (optimize(i)) THEN
          fm_nbonopt = fm_nbonopt + 1
          fm_bonopt(fm_nbonopt) = i
       ENDIF
    ENDDO
    ! ---------------------------------
    ! -------angles involving hydrogens
    optimize = .FALSE.
    old_ntypes = ntty
    DO i = 1, ntheh
       need_new = .FALSE.
       IF (is_qm(ith(i)).AND.is_qm(jth(i)).AND.is_qm(kth(i))) THEN
          ! this angle is QM only. Now loop over other angles of same type
          ! to see if this type occurs also as MM-QM or MM-MM
          DO j = 1, ntheh
             IF (j.EQ.i) CYCLE
             IF (icth(i).EQ.icth(j)) THEN
                IF (.NOT.is_qm(ith(j)) .OR. .NOT.is_qm(jth(j)) .OR.&
                     .NOT.is_qm(KTH(j))) THEN
                   need_new = .TRUE.
                   EXIT
                ENDIF
             ENDIF
          ENDDO
          ! now make sure that we do not have a earlier QM-QM angle of
          ! this type with NON-equivalent atoms
          IF (.NOT.need_new .AND. optimize(icth(i))) THEN
             DO j = 1, i-1
                IF (icth(i).EQ.icth(j)) THEN
                   IF (.NOT.are_equivalent('hangl',i,j)) THEN
                      need_new = .TRUE.
                      EXIT
                   ENDIF
                ENDIF
             ENDDO
          ENDIF
          IF (.NOT.need_new) THEN
             optimize(icth(i)) = .TRUE.
             CYCLE
          ENDIF
          old_type = icth(i)
          ntty = ntty + 1
          IF (ntty.GT.maxtty) THEN
             CALL stopgm('fm_setup_topology','NTTY.gt.MAXTTY',& 
                  __LINE__,__FILE__)
          ENDIF
          icth(i) = ntty
          ct(ntty) = ct(old_type)
          t0(ntty) = t0(old_type)
          ! check for equivalencies
          DO j = i+1, ntheh
             IF (icth(j) .EQ. old_type) THEN
                IF (is_qm(ith(j)) .AND. is_qm(jth(j)) .AND.&
                     is_qm(KTH(j))) THEN
                   IF (are_equivalent('hangl', i, j)) icth(j) = ntty
                ENDIF
             ENDIF
          ENDDO
          optimize(icth(i)) = .TRUE.
       ENDIF
    ENDDO
    IF (paral%io_parent)&
         WRITE(6, '(A, i6)') '  fm  New types for angles with H:     ',&
         NTTY - old_ntypes
    ! -------angles not involving hydrogens
    old_ntypes = ntty
    DO i = 1, nthe
       need_new = .FALSE.
       IF (is_qm(it(i)).AND.is_qm(jt(i)).AND.is_qm(kt(i))) THEN
          ! this angle is QM only. Now loop over other angles of same type
          ! to see if this type occurs also as MM-QM or MM-MM
          DO j = 1, nthe
             IF (j.EQ.i) CYCLE
             IF (ict(i).EQ.ict(j)) THEN
                IF (.NOT.is_qm(it(j)) .OR. .NOT.is_qm(jt(j)) .OR.&
                     .NOT.is_qm(KT(j))) THEN
                   need_new = .TRUE.
                   EXIT
                ENDIF
             ENDIF
          ENDDO
          ! now make sure that we do not have a earlier QM-QM angle of
          ! this type with NON-equivalent atoms
          IF (.NOT.need_new .AND. optimize(ict(i))) THEN
             DO j = 1, i-1
                IF (ict(i).EQ.ict(j)) THEN
                   IF (.NOT.are_equivalent('angl',i,j)) THEN
                      need_new = .TRUE.
                      EXIT
                   ENDIF
                ENDIF
             ENDDO
          ENDIF
          IF (.NOT.need_new) THEN
             optimize(ict(i)) = .TRUE.
             CYCLE
          ENDIF
          old_type = ict(i)
          ntty = ntty + 1
          IF (ntty.GT.maxtty) THEN
             CALL stopgm('fm_setup_topology','NTTY.gt.MAXTTY',& 
                  __LINE__,__FILE__)
          ENDIF
          ict(i) = ntty
          ct(ntty) = ct(old_type)
          t0(ntty) = t0(old_type)
          ! check for equivalencies
          DO j = i+1, nthe
             IF (ict(j) .EQ. old_type) THEN
                IF (is_qm(it(j)) .AND. is_qm(jt(j)) .AND.&
                     is_qm(KT(j))) THEN
                   IF (are_equivalent('angl', i, j)) ict(j) = ntty
                ENDIF
             ENDIF
          ENDDO
          optimize(ict(i)) = .TRUE.
       ENDIF
    ENDDO
    IF (paral%io_parent)&
         WRITE(6, '(A, i6)') '  fm  New types for other angles:      ',&
         NTTY - old_ntypes
    fm_nangopt = 0
    DO i = 1, ntty
       IF (optimize(i)) THEN
          fm_nangopt = fm_nangopt + 1
          fm_angopt(fm_nangopt) = i
       ENDIF
    ENDDO
    ! ---------------------------------------------
    ! -------improper dihedrals involving hydrogens
    optimize = .FALSE.
    old_ntypes = nqty
    DO i = 1, nqhih
       need_new = .FALSE.
       IF (is_qm(iqh(i)).AND.is_qm(jqh(i)).AND.is_qm(kqh(i)).AND.&
            is_qm(LQH(i))) THEN
          ! this angle is QM only. Now loop over other angles of same type
          ! to see if this type occurs also as MM-QM or MM-MM
          DO j = 1, nqhih
             IF (j.EQ.i) CYCLE
             IF (icqh(i).EQ.icqh(j)) THEN
                IF (.NOT.is_qm(iqh(j)) .OR. .NOT.is_qm(jqh(j)) .OR.&
                     .NOT.is_qm(KQH(j)) .OR. .NOT.is_qm(LQH(j))) THEN
                   need_new = .TRUE.
                   EXIT
                ENDIF
             ENDIF
          ENDDO
          IF (.NOT.need_new) THEN
             optimize(icqh(i)) = .TRUE.
             CYCLE
          ENDIF
          old_type = icqh(i)
          nqty = nqty + 1
          IF (nqty.GT.maxqty) THEN
             CALL stopgm('fm_setup_topology','NQTY.gt.MAXQTY',& 
                  __LINE__,__FILE__)
          ENDIF
          icqh(i) = nqty
          cq(nqty) = cq(old_type)
          q0(nqty) = q0(old_type)
          ! check for equivalencies
          DO j = i+1, nqhih
             IF (icqh(j) .EQ. old_type) THEN
                IF (is_qm(iqh(j)) .AND. is_qm(jqh(j)) .AND.&
                     is_qm(KQH(j)) .AND. is_qm(LQH(j))) THEN
                   IF (are_equivalent('himpd', i, j)) icqh(j) = nqty
                ENDIF
             ENDIF
          ENDDO
          optimize(icqh(i)) = .TRUE.
       ENDIF
    ENDDO
    IF (paral%io_parent)&
         WRITE(6, '(A, i6)') '  fm  New types for impropers with H:  ',&
         NQTY - old_ntypes
    ! -------improper dihedrals not involving hydrogens
    old_ntypes = nqty
    DO i = 1, nqhi
       need_new = .FALSE.
       IF (is_qm(iq(i)).AND.is_qm(jq(i)).AND.is_qm(kq(i)).AND.&
            is_qm(LQ(i))) THEN
          ! this angle is QM only. Now loop over other angles of same type
          ! to see if this type occurs also as MM-QM or MM-MM
          DO j = 1, nqhi
             IF (j.EQ.i) CYCLE
             IF (icq(i).EQ.icq(j)) THEN
                IF (.NOT.is_qm(iq(j)) .OR. .NOT.is_qm(jq(j)) .OR.&
                     .NOT.is_qm(KQ(j)) .OR. .NOT.is_qm(LQ(j))) THEN
                   need_new = .TRUE.
                   EXIT
                ENDIF
             ENDIF
          ENDDO
          IF (.NOT.need_new) THEN
             optimize(icq(i)) = .TRUE.
             CYCLE
          ENDIF
          old_type = icq(i)
          nqty = nqty + 1
          IF (nqty.GT.maxqty) THEN
             CALL stopgm('fm_setup_topology','NQTY.gt.MAXQTY',& 
                  __LINE__,__FILE__)
          ENDIF
          icq(i) = nqty
          cq(nqty) = cq(old_type)
          q0(nqty) = q0(old_type)
          ! check for equivalencies
          DO j = i+1, nqhi
             IF (icq(j) .EQ. old_type) THEN
                IF (is_qm(iq(j)) .AND. is_qm(jq(j)) .AND.&
                     is_qm(KQ(j)) .AND. is_qm(LQ(j))) THEN
                   IF (are_equivalent('impd', i, j)) icq(j) = nqty
                ENDIF
             ENDIF
          ENDDO
          optimize(icq(i)) = .TRUE.
       ENDIF
    ENDDO
    IF (paral%io_parent)&
         WRITE(6, '(A, i6)') '  fm  New types for other impropers:   ',&
         NQTY - old_ntypes
    fm_nimpopt = 0
    DO i = 1, nqty
       IF (optimize(i)) THEN
          fm_nimpopt = fm_nimpopt + 1
          fm_impopt(fm_nimpopt) = i
       ENDIF
    ENDDO
    ! ----------------------------------------------
    ! -------cntl%proper dihedrals involving hydrogens
    ! for dihedrals, we only duplicate a type if it has all-QM and QM/MM or
    ! all-MM atoms, but we do not create new types if it involves
    ! non-equivalent QM atoms
    ! we do not fit dihedrals with force constant zero (they occur in Amber;
    ! I think they are used to signal 1-4 interactions)
    optimize = .FALSE.
    old_ntypes = npty
    DO i = 1, nphih
       IF (cp(icph(i)) .LE. 1.0e-3_real_8) CYCLE
       need_new = .FALSE.
       IF (is_qm(iph(i)).AND.is_qm(jph(i)).AND.is_qm(kph(i)).AND.&
            is_qm(LPH(i))) THEN
          ! this angle is QM only. Now loop over other angles of same type
          ! to see if this type occurs also as MM-QM or MM-MM
          DO j = 1, nphih
             IF (j.EQ.i) CYCLE
             IF (icph(i).EQ.icph(j)) THEN
                IF (.NOT.is_qm(iph(j)) .OR. .NOT.is_qm(jph(j)) .OR.&
                     .NOT.is_qm(KPH(j)) .OR. .NOT.is_qm(LPH(j))) THEN
                   need_new = .TRUE.
                   EXIT
                ENDIF
             ENDIF
          ENDDO
          IF (.NOT.need_new) THEN
             optimize(icph(i)) = .TRUE.
             CYCLE
          ENDIF
          old_type = icph(i)
          npty = npty + 1
          IF (npty.GT.maxpty) THEN
             CALL stopgm('fm_setup_topology','NPTY.gt.MAXPTY',& 
                  __LINE__,__FILE__)
          ENDIF
          icph(i) = npty
          cp(npty) = cp(old_type)
          pd(npty) = pd(old_type)
          np(npty) = np(old_type)
          ! check for equivalencies
          DO j = i+1, nphih
             IF (icph(j) .EQ. old_type) THEN
                IF (is_qm(iph(j)) .AND. is_qm(jph(j)) .AND.&
                     is_qm(KPH(j)) .AND. is_qm(LPH(j))) THEN
                   ! if(are_equivalent('hdihe', i, j)) ICPH(j) = NPTY
                   icph(j) = npty
                ENDIF
             ENDIF
          ENDDO
          optimize(icph(i)) = .TRUE.
       ENDIF
    ENDDO
    IF (paral%io_parent)&
         WRITE(6, '(A, i6)') '  fm  New types for dihedrals with H:  ',&
         NPTY - old_ntypes
    ! -------cntl%proper dihedrals not involving hydrogens
    old_ntypes = npty
    DO i = 1, nphi
       IF (cp(icp(i)) .LE. 1.0e-3_real_8) CYCLE
       need_new = .FALSE.
       IF (is_qm(ip(i)).AND.is_qm(jp(i)).AND.is_qm(kp(i)).AND.&
            is_qm(LP(i))) THEN
          ! this angle is QM only. Now loop over other angles of same type
          ! to see if this type occurs also as MM-QM or MM-MM
          DO j = 1, nphi
             IF (j.EQ.i) CYCLE
             IF (icp(i).EQ.icp(j)) THEN
                IF (.NOT.is_qm(ip(j)) .OR. .NOT.is_qm(jp(j)) .OR.&
                     .NOT.is_qm(KP(j)) .OR. .NOT.is_qm(LP(j))) THEN
                   need_new = .TRUE.
                   EXIT
                ENDIF
             ENDIF
          ENDDO
          IF (.NOT.need_new) THEN
             optimize(icp(i)) = .TRUE.
             CYCLE
          ENDIF
          old_type = icp(i)
          npty = npty + 1
          IF (npty.GT.maxpty) THEN
             CALL stopgm('fm_setup_topology','NPTY.gt.MAXPTY',& 
                  __LINE__,__FILE__)
          ENDIF
          icp(i) = npty
          cp(npty) = cp(old_type)
          pd(npty) = pd(old_type)
          np(npty) = np(old_type)
          ! check for equivalencies
          DO j = i+1, nphi
             IF (icp(j) .EQ. old_type) THEN
                IF (is_qm(ip(j)) .AND. is_qm(jp(j)) .AND.&
                     is_qm(KP(j)) .AND. is_qm(LP(j))) THEN
                   ! if(are_equivalent('dihe', i, j)) ICP(j) = NPTY
                   icp(j) = npty
                ENDIF
             ENDIF
          ENDDO
          optimize(icp(i)) = .TRUE.
       ENDIF
    ENDDO
    IF (paral%io_parent)&
         WRITE(6, '(A, i6)') '  fm  New types for other dihedrals:   ',&
         NPTY - old_ntypes
    fm_ndihopt = 0
    DO i = 1, npty
       IF (optimize(i)) THEN
          fm_ndihopt = fm_ndihopt + 1
          fm_dihopt(fm_ndihopt) = i
       ENDIF
    ENDDO
    ! ---------------------------------------------------------
    ! --Done, finally!
    ! Some printing out
    IF (paral%io_parent)&
         WRITE(6, '(A, i6)') '  fm'
    IF (paral%io_parent)&
         WRITE(6, '(A, i6)') '  fm  Nr. of bond types to optimize:     ',&
         fm_nbonopt
    IF (paral%io_parent)&
         WRITE(6, '(A, i6)') '  fm  Nr. of angle types to optimize:    ',&
         fm_nangopt
    IF (paral%io_parent)&
         WRITE(6, '(A, i6)') '  fm  Nr. of improper types to optimize: ',&
         fm_nimpopt
    IF (paral%io_parent)&
         WRITE(6, '(A, i6)') '  fm  Nr. of dihedral types to optimize: ',&
         fm_ndihopt
    IF (fm_nbonopt .GT. 0) THEN
       IF (paral%io_parent)&
            WRITE(6, '(A)') '  fm'
       IF (paral%io_parent)&
            WRITE(6, '(A)') '  fm  List of optimized bond types:'
       IF (paral%io_parent)&
            WRITE(6, 211) fm_bonopt(1:fm_nbonopt)
    ENDIF
    IF (fm_nangopt .GT. 0) THEN
       IF (paral%io_parent)&
            WRITE(6, '(A)') '  fm  List of optimized angle types:'
       IF (paral%io_parent)&
            WRITE(6, 211) fm_angopt(1:fm_nangopt)
    ENDIF
    IF (fm_nimpopt .GT. 0) THEN
       IF (paral%io_parent)&
            WRITE(6, '(A)') '  fm  List of optimized improper types:'
       IF (paral%io_parent)&
            WRITE(6, 211) fm_impopt(1:fm_nimpopt)
    ENDIF
    IF (fm_ndihopt .GT. 0) THEN
       IF (paral%io_parent)&
            WRITE(6, '(A)') '  fm  List of optimized dihedral types:'
       IF (paral%io_parent)&
            WRITE(6, 211) fm_dihopt(1:fm_ndihopt)
    ENDIF
211 FORMAT('  fm  ',9(i8))
    IF (paral%io_parent)&
         WRITE(6,*)&
         ' fm fm fm fm fm fm fm fm fm fm fm fm fm fm fm fm fm fm fm fm fm'


    ! set force constants of QM bonds, etc., to zero (needed 
    ! for QM/MM), but store values for later use as an initial guess

    ALLOCATE(fm_kbon(fm_nbonopt),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    DO i = 1, fm_nbonopt
       fm_kbon(i) = cb(fm_bonopt(i))
       cb(fm_bonopt(i)) = 0._real_8
    ENDDO

    ALLOCATE(fm_kang(fm_nangopt),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    DO i = 1, fm_nangopt
       fm_kang(i) = ct(fm_angopt(i))
       ct(fm_angopt(i)) = 0._real_8
    ENDDO

    ALLOCATE(fm_kimp(fm_nimpopt),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    DO i = 1, fm_nimpopt
       fm_kimp(i) = cq(fm_impopt(i))
       cq(fm_impopt(i)) = 0._real_8
    ENDDO

    ALLOCATE(fm_kdih(fm_ndihopt),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    DO i = 1, fm_ndihopt
       fm_kdih(i) = cp(fm_dihopt(i))
       cp(fm_dihopt(i)) = 0._real_8
    ENDDO

    DEALLOCATE(optimize,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
#endif
    RETURN
  END SUBROUTINE fm_setup_topology
  ! 
  ! ==================================================================
  FUNCTION are_equivalent(ia, i1, i2) RESULT(my_res)
    ! ==--------------------------------------------------------------==
    ! checks if bonds, angles, etc. are between equivalent atoms
    ! ==================================================================
    ! type of interaction
    CHARACTER(*)                             :: ia
    INTEGER                                  :: i1, i2
    LOGICAL                                  :: my_res

    INTEGER                                  :: a1, a2, b1, b2, c1, c2, d1, d2

! nr of interactions

#if defined (__GROMOS)       
    INCLUDE 'gromos.h'

    ! first get unique type number
    IF (ia.EQ.'hbond') THEN
       a1 = grm_equiv(ibh(i1))
       b1 = grm_equiv(jbh(i1))
       a2 = grm_equiv(ibh(i2))
       b2 = grm_equiv(jbh(i2))
    ELSEIF (ia.EQ.'bond') THEN
       a1 = grm_equiv(ib(i1))
       b1 = grm_equiv(jb(i1))
       a2 = grm_equiv(ib(i2))
       b2 = grm_equiv(jb(i2))
    ELSEIF (ia.EQ.'hangl') THEN
       a1 = grm_equiv(ith(i1))
       b1 = grm_equiv(jth(i1))
       c1 = grm_equiv(kth(i1))
       a2 = grm_equiv(ith(i2))
       b2 = grm_equiv(jth(i2))
       c2 = grm_equiv(kth(i2))
    ELSEIF (ia.EQ.'angl') THEN
       a1 = grm_equiv(it(i1))
       b1 = grm_equiv(jt(i1))
       c1 = grm_equiv(kt(i1))
       a2 = grm_equiv(it(i2))
       b2 = grm_equiv(jt(i2))
       c2 = grm_equiv(kt(i2))
    ELSEIF (ia.EQ.'himpd') THEN
       a1 = grm_equiv(iqh(i1))
       b1 = grm_equiv(jqh(i1))
       c1 = grm_equiv(kqh(i1))
       d1 = grm_equiv(lqh(i1))
       a2 = grm_equiv(iqh(i2))
       b2 = grm_equiv(jqh(i2))
       c2 = grm_equiv(kqh(i2))
       d2 = grm_equiv(lqh(i2))
    ELSEIF (ia.EQ.'impd') THEN
       a1 = grm_equiv(iq(i1))
       b1 = grm_equiv(jq(i1))
       c1 = grm_equiv(kq(i1))
       d1 = grm_equiv(lq(i1))
       a2 = grm_equiv(iq(i2))
       b2 = grm_equiv(jq(i2))
       c2 = grm_equiv(kq(i2))
       d2 = grm_equiv(lq(i2))
    ELSEIF (ia.EQ.'hdihe') THEN
       a1 = grm_equiv(iph(i1))
       b1 = grm_equiv(jph(i1))
       c1 = grm_equiv(kph(i1))
       d1 = grm_equiv(lph(i1))
       a2 = grm_equiv(iph(i2))
       b2 = grm_equiv(jph(i2))
       c2 = grm_equiv(kph(i2))
       d2 = grm_equiv(lph(i2))
    ELSEIF (ia.EQ.'dihe') THEN
       a1 = grm_equiv(ip(i1))
       b1 = grm_equiv(jp(i1))
       c1 = grm_equiv(kp(i1))
       d1 = grm_equiv(lp(i1))
       a2 = grm_equiv(ip(i2))
       b2 = grm_equiv(jp(i2))
       c2 = grm_equiv(kp(i2))
       d2 = grm_equiv(lp(i2))
    ELSE
       ! we should never arrive here...
       CALL stopgm('are_equivalent', 'internal error! unknown ia: '&
            // ia,& 
            __LINE__,__FILE__)
    ENDIF

    ! bonds
    IF (ia.EQ.'bond' .OR. ia.EQ.'hbond') THEN
       IF (a1.EQ.a2 .AND. b1.EQ.b2) THEN
          my_res = .TRUE.
       ELSEIF (a1.EQ.b2 .AND. a2.EQ.b1) THEN
          my_res = .TRUE.
       ELSE
          my_res = .FALSE.
       ENDIF
       RETURN
       ! angles
    ELSEIF (ia.EQ.'angl' .OR. ia.EQ.'hangl') THEN
       IF (b1.EQ.b2) THEN
          IF (a1.EQ.a2 .AND. c1.EQ.c2) THEN
             my_res = .TRUE.
             ! elseif (a1.eq.b2 .and. a2.eq.b1) then
          ELSEIF (a1.EQ.c2 .AND. a2.EQ.c1) THEN
             my_res = .TRUE.
          ELSE
             my_res = .FALSE.
          ENDIF
       ELSE
          my_res = .FALSE.
       ENDIF
       RETURN
       ! dihedrals (cntl%proper or not)
    ELSEIF( (ia.EQ.'hdihe' .OR. ia.EQ.'dihe') .OR.&
         (ia.EQ.'himpd' .OR. ia.EQ.'impd') ) THEN
       IF (a1.EQ.a2 .AND. b1.EQ.b2 .AND. c1.EQ.c2 .AND. d1.EQ.d2)THEN
          my_res = .TRUE.
       ELSEIF(a1.EQ.d2 .AND. b1.EQ.c2 .AND. c1.EQ.b2 .AND. d1.EQ.a2)&
            THEN
          my_res = .TRUE.
       ELSE
          my_res = .FALSE.
       ENDIF

    ELSE
       ! we should never arrive here...
       CALL stopgm('are_equivalent', 'internal error! unknown ia: '&
            // ia,& 
            __LINE__,__FILE__)
    ENDIF
#else
    my_res = .FALSE.
#endif
  END FUNCTION are_equivalent
  ! 
  ! ==================================================================
  SUBROUTINE fm_wopt(trget_size, nw, w, fdiff, iflag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: trget_size, nw
    REAL(real_8)                             :: w(nw), fdiff(trget_size)
    INTEGER                                  :: iflag

  END SUBROUTINE fm_wopt
  ! 
  ! ==================================================================
  SUBROUTINE fm_fcn_allcovalent(tgsize, npar, par, fdiff, fjac,&
       ldfjac, iflag)
    ! ==--------------------------------------------------------------==
    ! size of the target vector, nr. of parameters, leading dimension
    ! of fjac
    INTEGER                                  :: tgsize, npar
    REAL(real_8)                             :: par(npar)
    REAL(real_8), TARGET                     :: fdiff(tgsize)
    INTEGER                                  :: ldfjac
    REAL(real_8), TARGET                     :: fjac(ldfjac, npar)
    INTEGER                                  :: iflag

    INTEGER                                  :: at, bt, dt, i, iang, ibon, &
                                                idih, iframe, iimp, ipar, &
                                                ipar0, it, j, k, l, lb, ub
    INTEGER, SAVE                            :: itcounter = 0
    REAL(real_8)                             :: fi(3), fj(3), fk(3), fl(3), &
                                                ri(3), rj(3), rk(3), rl(3), &
                                                rms, rms_rel
    REAL(real_8), POINTER                    :: fvec(:), jac(:, :)

! parameters, target vector, jacobian matrix
! flag that determines wether to compute fdiff or fjac
! Variables
! point to windows of fdiff and fjac
! atomic coordinates
! forces (or derivatives wrt. parameters)
! absolute and relative RMS on forces
! atom indices, type indices
! iteration counter
! loop counters, temporaries, etc.

    IF (iflag.EQ.1) itcounter = itcounter + 1
    IF (paral%io_parent)&
         WRITE(6,'(A)') '  fm  Entering fm_fcn_allcovalent'
    IF (paral%io_parent)&
         WRITE(6,'(A,i3)') '  fm  iflag is ', iflag
    IF (paral%io_parent)&
         WRITE(6,'(A,i8)') '  fm  iteration ', itcounter
    IF (paral%io_parent)&
         WRITE(6,'(A,i8,i8)') '  fm  size of jacobian is ', ldfjac, npar
    ! consistency check
    IF (tgsize.NE.ldfjac) THEN
       IF (paral%io_parent)&
            WRITE(6,'(A)') '  fm  =================================='
       IF (paral%io_parent)&
            WRITE(6,'(A)') '  fm  ========== ERROR ! ==============='
       IF (paral%io_parent)&
            WRITE(6,'(A)') '  fm  =================================='
       IF (paral%io_parent)&
            WRITE(6,'(A)') '  fm   tgsize must be equal ldfjac, but '
       IF (paral%io_parent)&
            WRITE(6,*) '  fm   they are: ', tgsize, ldfjac
       CALL stopgm('fm_fcn_allcovalent','tgsize .ne. ldfjac',& 
            __LINE__,__FILE__)
    ENDIF
    IF (paral%io_parent)&
         WRITE(6,'(A)') '  fm  Entering fm_fcn_allcovalent'
    IF (paral%io_parent)&
         WRITE(6,'(A,i3)') '  fm  iflag is ', iflag
    IF (paral%io_parent)&
         WRITE(6,'(A,i8)') '  fm  iteration ', itcounter

    IF (itcounter.GT.fm_max_kfit_iterations) THEN
       iflag = -1
       RETURN
    ENDIF

    IF (iflag .EQ. 1) THEN
       CALL zeroing(fdiff)!, tgsize)
    ELSE IF (iflag .EQ. 2) THEN
       CALL zeroing(fjac)!, ldfjac*npar)
    ELSE
       IF (paral%io_parent)&
            WRITE(6,'(A)') '  fm  =================================='
       IF (paral%io_parent)&
            WRITE(6,'(A)') '  fm  ========== ERROR ! ==============='
       IF (paral%io_parent)&
            WRITE(6,'(A,i8)') '  fm  Illegal value for iflag: ', iflag
       CALL stopgm('fm_fcn_allcovalent', 'illegal value for iflag',& 
            __LINE__,__FILE__)
    ENDIF

    CALL fm_kfit_par_to_fmarray(npar, par)

    ! Loop over all frames. In each frame, we use fvec and jac to point to
    ! the appropriate "windows" of the full vector fdiff / full matrix fjac,
    ! to facilitate the indexing business
    DO iframe = 1, fm_nframes
       lb = (iframe-1) * 3*fm_nrqm + 1
       ub = lb + 3*fm_nrqm - 1
       fvec => fdiff(lb:ub)
       jac => fjac(lb:ub, :)

       ! Bonds
       ipar0 = 0
       DO ibon = 1, fm_nbon
          bt = fm_ib(0, ibon)
          i = fm_ib(1, ibon)
          j = fm_ib(2, ibon)
          ri = fm_rqm(:, i, iframe)
          rj = fm_rqm(:, j, iframe)

          ! _DEBUG[
          ! write(6,*) "DEBUG"
          ! write(6,*) "ibon, bt, i, j:"
          ! write(6,*) ibon, bt, i, j
          ! _DEBUG]
          ! compute i,j as index into array fvec
          i = 3*(i-1) + 1
          j = 3*(j-1) + 1
          IF (iflag.EQ.1) THEN
             CALL fm_forces_bon(fi, fj, ri, rj, bt, 'ff')
             fvec(i:i+2) = fvec(i:i+2) + fi
             fvec(j:j+2) = fvec(j:j+2) + fj
          ELSE
             CALL fm_forces_bon(fi, fj, ri, rj, bt, 'dk')
             ipar = ipar0 + 2*bt - 1
             jac(i:i+2, ipar) = jac(i:i+2, ipar) + fi
             jac(j:j+2, ipar) = jac(j:j+2, ipar) + fj
             CALL fm_forces_bon(fi, fj, ri, rj, bt, 'd0')
             ipar = ipar + 1
             jac(i:i+2, ipar) = jac(i:i+2, ipar) + fi
             jac(j:j+2, ipar) = jac(j:j+2, ipar) + fj
          ENDIF
       ENDDO
       ! Angles
       ipar0 = 2*fm_nbonopt
       DO iang = 1, fm_nang
          at = fm_it(0, iang)
          i = fm_it(1, iang)
          j = fm_it(2, iang)
          k = fm_it(3, iang)
          ri = fm_rqm(:, i, iframe)
          rj = fm_rqm(:, j, iframe)
          rk = fm_rqm(:, k, iframe)

          i = 3*(i-1) + 1
          j = 3*(j-1) + 1
          k = 3*(k-1) + 1
          IF (iflag.EQ.1) THEN
             CALL fm_forces_ang(fi, fj, fk, ri, rj, rk, at, 'ff')
             fvec(i:i+2) = fvec(i:i+2) + fi
             fvec(j:j+2) = fvec(j:j+2) + fj
             fvec(k:k+2) = fvec(k:k+2) + fk
          ELSE
             ipar = ipar0 + 2*at - 1
             CALL fm_forces_ang(fi, fj, fk, ri, rj, rk, at, 'dk')
             jac(i:i+2, ipar) = jac(i:i+2, ipar) + fi
             jac(j:j+2, ipar) = jac(j:j+2, ipar) + fj
             jac(k:k+2, ipar) = jac(k:k+2, ipar) + fk
             ipar = ipar + 1
             CALL fm_forces_ang(fi, fj, fk, ri, rj, rk, at, 'd0')
             jac(i:i+2, ipar) = jac(i:i+2, ipar) + fi
             jac(j:j+2, ipar) = jac(j:j+2, ipar) + fj
             jac(k:k+2, ipar) = jac(k:k+2, ipar) + fk
          ENDIF
       ENDDO
       ! Impropers
       ipar0 = ipar0 + 2*fm_nangopt
       DO iimp = 1, fm_nimp
          it = fm_iq(0, iimp)
          i = fm_iq(1, iimp)
          j = fm_iq(2, iimp)
          k = fm_iq(3, iimp)
          l = fm_iq(4, iimp)
          ri = fm_rqm(:, i, iframe)
          rj = fm_rqm(:, j, iframe)
          rk = fm_rqm(:, k, iframe)
          rl = fm_rqm(:, l, iframe)

          i = 3*(i-1) + 1
          j = 3*(j-1) + 1
          k = 3*(k-1) + 1
          l = 3*(l-1) + 1
          IF (iflag.EQ.1) THEN
             CALL fm_forces_imp(fi,fj,fk,fl,ri,rj,rk,rl,it,'ff')
             fvec(i:i+2) = fvec(i:i+2) + fi
             fvec(j:j+2) = fvec(j:j+2) + fj
             fvec(k:k+2) = fvec(k:k+2) + fk
             fvec(l:l+2) = fvec(l:l+2) + fl
          ELSE
             ipar = ipar0 + it
             CALL fm_forces_imp(fi,fj,fk,fl,ri,rj,rk,rl,it,'dk')
             jac(i:i+2, ipar) = jac(i:i+2, ipar) + fi
             jac(j:j+2, ipar) = jac(j:j+2, ipar) + fj
             jac(k:k+2, ipar) = jac(k:k+2, ipar) + fk
             jac(l:l+2, ipar) = jac(l:l+2, ipar) + fl
          ENDIF
       ENDDO
       ! Dihedrals
       ipar0 = ipar0 + fm_nimpopt
       DO idih = 1, fm_ndih
          dt = fm_ip(0, idih)
          i = fm_ip(1, idih)
          j = fm_ip(2, idih)
          k = fm_ip(3, idih)
          l = fm_ip(4, idih)
          ri = fm_rqm(:, i, iframe)
          rj = fm_rqm(:, j, iframe)
          rk = fm_rqm(:, k, iframe)
          rl = fm_rqm(:, l, iframe)

          i = 3*(i-1) + 1
          j = 3*(j-1) + 1
          k = 3*(k-1) + 1
          l = 3*(l-1) + 1
          IF (iflag.EQ.1) THEN
             CALL fm_forces_dih(fi,fj,fk,fl,ri,rj,rk,rl,dt,'ff')
             fvec(i:i+2) = fvec(i:i+2) + fi
             fvec(j:j+2) = fvec(j:j+2) + fj
             fvec(k:k+2) = fvec(k:k+2) + fk
             fvec(l:l+2) = fvec(l:l+2) + fl
          ELSE
             ipar = ipar0 + dt
             CALL fm_forces_dih(fi,fj,fk,fl,ri,rj,rk,rl,dt,'dk')
             jac(i:i+2, ipar) = jac(i:i+2, ipar) + fi
             jac(j:j+2, ipar) = jac(j:j+2, ipar) + fj
             jac(k:k+2, ipar) = jac(k:k+2, ipar) + fk
             jac(l:l+2, ipar) = jac(l:l+2, ipar) + fl
          ENDIF
       ENDDO
    ENDDO! loop over all frames

    IF (iflag.EQ.1) THEN
       fdiff = fdiff - RESHAPE(fm_fqm_ref, (/tgsize/))
       CALL fm_kfit_compute_frms(rms, rms_rel, tgsize, fdiff,&
            fm_fqm_ref)
       IF (paral%io_parent)&
            WRITE(6,'(A)') '  fm  Abs. and rel. RMS on forces:'
       IF (paral%io_parent)&
            WRITE(6,'(e10.4,1x,f8.3)') rms, rms_rel
    ENDIF

    RETURN
  END SUBROUTINE fm_fcn_allcovalent
  ! 
  ! ==================================================================
  SUBROUTINE fm_kfit_compute_frms(rms, rms_rel, fsize, fdiff, fref)
    ! ==--------------------------------------------------------------==
    ! output: rms and relative rms
    REAL(real_8)                             :: rms, rms_rel
    INTEGER                                  :: fsize
    REAL(real_8)                             :: fdiff(fsize), fref(fsize)

    INTEGER                                  :: i

! input: size of arrays, fdiff (force - ref.force), reference forces
! Variables
! loop counter

    rms = 0._real_8
    rms_rel = 0._real_8

    ! we abuse rms_rel to sum over square reference forces
    DO i = 1, fsize
       rms = rms + fdiff(i)**2
       rms_rel = rms_rel + fref(i)**2
    ENDDO
    rms_rel = SQRT(rms / rms_rel)
    rms = SQRT(rms / fsize)
    RETURN
  END SUBROUTINE fm_kfit_compute_frms
  ! 
  ! ==================================================================
  SUBROUTINE fm_forces_bon(fi, fj, ri, rj, bt, flag)
    ! ==--------------------------------------------------------------==
    ! output: forces or derivatives, depending on flag
    REAL(real_8)                             :: fi(3), fj(3), ri(3), rj(3)
    INTEGER                                  :: bt
    CHARACTER(len=2)                         :: flag

    REAL(real_8)                             :: b0, kb, r, rij(3)

! atomic coordinates
! nr of bond type
! flag that determines type of calculation:
! 'ff': calculate forces
! 'dk': derivative wrt. to kb
! 'd0': derivative wrt. to b0
! Variables
! distance, distance vector
! force constant, eq. bond length

    CALL fm_dist(r, rij, ri, rj)
    kb = fm_kb(bt)
    b0 = fm_b0(bt)
    IF (flag.EQ.'ff') THEN
       IF (fm_amber) THEN
          fi = -kb * (r - b0) * rij / r
          fj = -fi
       ELSE
          fi = -kb * (r**2 - b0**2) * rij
          fj = -fi
       ENDIF
    ELSEIF (flag.EQ.'dk') THEN
       IF (fm_amber) THEN
          fi = -(r - b0) * rij / r
          fj = -fi
       ELSE
          fi = -(r**2 - b0**2) * rij
          fj = -fi
       ENDIF
    ELSEIF (flag.EQ.'d0') THEN
       IF (fm_amber) THEN
          fi = kb * rij / r
          fj = -fi
       ELSE
          fi = 2._real_8 * kb * b0 * rij
          fj = -fi
       ENDIF
    ELSE
       CALL stopgm('fm_forces_bon', 'Unknown flag: ' // flag,& 
            __LINE__,__FILE__)
    ENDIF

    RETURN
  END SUBROUTINE fm_forces_bon
  ! 
  ! ==================================================================
  SUBROUTINE fm_forces_ang(fi,fj,fk,ri,rj,rk,at,flag)
    ! ==--------------------------------------------------------------==
    ! forces or derivatives
    REAL(real_8)                             :: fi(3), fj(3), fk(3), ri(3), &
                                                rj(3), rk(3)
    INTEGER                                  :: at
    CHARACTER(len=2)                         :: flag

    REAL(real_8)                             :: dij, dkj, fact, fscal, kt, &
                                                rij(3), rkj(3), t0, th

! coordinates
! angle type
! flag that determines what is calculated:
! 'ff': compute forces
! 'dk': derivative wrt. kt
! 'd0': derivative wrt. t0
! Variables
! distances, distance vectors
! angle or cosine of angle
! force constant, equilibrium angle
! scalar factor for amber
! compute angle

    CALL fm_dist(dij, rij, ri, rj)
    CALL fm_dist(dkj, rkj, rk, rj)
    th = DOT_PRODUCT(rij, rkj) / (dij*dkj)

    kt = fm_kt(at)
    t0 = fm_t0(at)
    IF (fm_amber) THEN
       th = ACOS(th)
       t0 = ACOS(t0)
    ENDIF

    IF (flag.EQ.'ff' .OR. flag.EQ.'dk') THEN
       IF (fm_amber) THEN
          fscal = (th - t0) / SIN(th)
          fi = ( rkj*dij*dkj - DOT_PRODUCT(rij,rkj)*dkj*rij/dij )&
               / (dij**2 * dkj**2)
          fi = fscal * fi
          fk = ( rij*dij*dkj - DOT_PRODUCT(rij,rkj)*dij*rkj/dkj )&
               / (dij**2 * dkj**2)
          fk = fscal * fk
          fj = -fi - fk
       ELSE
          fi = -(th - t0) * (rkj/dkj - (rij/dij)*th) / dij
          fk = -(th - t0) * (rij/dij - (rkj/dkj)*th) / dkj
          fj = -fi - fk
       ENDIF
       IF (flag.EQ.'ff') THEN
          fi = kt * fi
          fj = kt * fj
          fk = kt * fk
       ENDIF
    ELSEIF (flag.EQ.'d0') THEN
       IF (fm_amber) THEN
          fi = ( rkj*dij*dkj - DOT_PRODUCT(rij,rkj)*dkj*rij/dij )&
               / (dij**2 * dkj**2)
          fi = kt * fi
          fk = ( rij*dij*dkj - DOT_PRODUCT(rij,rkj)*dij*rkj/dkj )&
               / (dij**2 * dkj**2)
          fk = kt * fk
          fj = -fi - fk
       ELSE
          fact=1._real_8-fm_t0(at)**2
          fi = fact * kt * (rkj/dkj - (rij/dij)*th) / dij
          fk = fact * kt * (rij/dij - (rkj/dkj)*th) / dkj
          fj = -fi - fk
       ENDIF
    ELSE
       CALL stopgm('fm_forces_ang', 'Unknown flag: ' // flag,& 
            __LINE__,__FILE__)
    ENDIF

    RETURN
  END SUBROUTINE fm_forces_ang
  ! 
  ! ==================================================================
  SUBROUTINE fm_forces_imp(fi,fj,fk,fl,ri,rj,rk,rl,it,flag)
    ! ==--------------------------------------------------------------==
    ! forces or derivatives
    REAL(real_8)                             :: fi(3), fj(3), fk(3), fl(3), &
                                                ri(3), rj(3), rk(3), rl(3)
    INTEGER                                  :: it
    CHARACTER(len=2)                         :: flag

    REAL(real_8)                             :: dij, dkj, dkl, dmj, dnk, kq, &
                                                q, q0, rij(3), rkj(3), &
                                                rkl(3), rmj(3), rnk(3), &
                                                signzeta, zeta

! coordinates
! improper type
! flag that determines what is calculated:
! 'ff': compute forces
! 'dk': derivative wrt. kt
! Variables
! distances, distance vectors
! improper angle and its sign
! impropers in GROMOS lingo
! compute improper angle zeta

    CALL fm_dist(dij, rij, ri, rj)
    CALL fm_dist(dkj, rkj, rk, rj)
    CALL fm_dist(dkl, rkl, rk, rl)

    CALL fm_crossproduct(rmj, rij, rkj)
    CALL fm_crossproduct(rnk, rkj, rkl)
    dmj = SQRT(DOT_PRODUCT(rmj, rmj))
    dnk = SQRT(DOT_PRODUCT(rnk, rnk))

    signzeta = SIGN(1._real_8, DOT_PRODUCT(rij,rnk))
    zeta = DOT_PRODUCT(rmj, rnk) / (dmj*dnk)
    zeta = signzeta * ACOS(zeta)

    q = zeta
    q0 = fm_q0(it)

    fi = -(q - q0) * rmj * dkj/dmj**2
    fl = +(q - q0) * rnk * dkj/dnk**2
    fj = ( DOT_PRODUCT(rij,rkj)/rkj**2 - 1._real_8) * fi&
         - ( DOT_PRODUCT(rkl,rkj)/rkj**2 ) * fl
    fk = -fi - fj - fl
    IF (flag .EQ. 'dk') THEN
       CONTINUE! we're done
    ELSE IF (flag .EQ. 'ff') THEN
       kq = fm_kq(it)
       fi = kq * fi
       fj = kq * fj
       fk = kq * fk
       fl = kq * fl
    ELSE
       CALL stopgm('fm_forces_imp', 'Unknown flag: ' // flag,& 
            __LINE__,__FILE__)
    ENDIF

    RETURN
  END SUBROUTINE fm_forces_imp
  ! 
  ! ==================================================================
  SUBROUTINE fm_forces_dih(fi,fj,fk,fl,ri,rj,rk,rl,dt,flag)
    ! ==--------------------------------------------------------------==
    ! forces or derivatives
    REAL(real_8)                             :: fi(3), fj(3), fk(3), fl(3), &
                                                ri(3), rj(3), rk(3), rl(3)
    INTEGER                                  :: dt
    CHARACTER(len=2)                         :: flag

    INTEGER                                  :: m
    REAL(real_8) :: cosmphi, cosphi, dcosmphi, dij, dim, dkj, dkl, dln, kp, &
      phase, rij(3), rim(3), rkj(3), rkl(3), rln(3), rt(3), signphi

! coordinates
! dihedral type
! flag that determines what is calculated:
! 'ff': compute forces
! 'dk': derivative wrt. kt
! Variables
! distances, distance vectors
! dihedral angle and sign, cosine(phi) and cos(m*phi)/cos(phi),
! phase (see GROMOS bible page II-27)
! multiplicity, force constant
! compute dihedral angle phi (or rather the cosine of it, that's all we
! need..

    CALL fm_dist(dij, rij, ri, rj)
    CALL fm_dist(dkj, rkj, rk, rj)
    CALL fm_dist(dkl, rkl, rk, rl)

    rim = rij - DOT_PRODUCT(rij,rkj) * rkj / dkj**2
    rln = -rkl + DOT_PRODUCT(rkl,rkj) * rkj / dkj**2

    dim = SQRT(DOT_PRODUCT(rim,rim))
    dln = SQRT(DOT_PRODUCT(rln,rln))

    CALL fm_crossproduct(rt, rkj, rkl)
    signphi = SIGN(1._real_8, DOT_PRODUCT(rij, rt))
    cosphi = DOT_PRODUCT(rim, rln) / (dim*dln)

    ! now compute forces 
    phase = fm_p0(dt)
    m = INT(fm_pn(dt))
    CALL fm_cosder(dcosmphi, cosmphi, cosphi, m)

    fi = -phase * dcosmphi * ( rln/dln - cosphi*rim/dim) / dim
    fl = -phase * dcosmphi * ( rim/dim - cosphi*rln/dln) / dln
    fj = ( DOT_PRODUCT(rij,rkj) / dkj**2 - 1._real_8) * fi&
         - DOT_PRODUCT(rkl,rkj) / dkj**2 * fl
    fk = -fi - fj - fl

    IF (flag .EQ. 'dk') THEN
       CONTINUE! we're done already
    ELSE IF (flag .EQ. 'ff') THEN
       kp = fm_kp(dt)
       fi = kp * fi
       fj = kp * fj
       fk = kp * fk
       fl = kp * fl
    ELSE
       CALL stopgm('fm_forces_dih', 'Unknown flag: ' // flag,& 
            __LINE__,__FILE__)
    ENDIF

    RETURN
  END SUBROUTINE fm_forces_dih
  ! 
  ! ==================================================================
  SUBROUTINE fm_kfit_par_to_fmarray(npar, par)
    ! ==--------------------------------------------------------------==
    ! Write parameters obtained from Optimizer back into fm_arrays
    ! ==================================================================
    ! nr. of parameters, params themselves
    INTEGER                                  :: npar
    REAL(real_8)                             :: par(npar)

    INTEGER                                  :: i, ipar

! Variables
! loop counters, etc.

    ipar = 0
    ! Bonds
    DO i = 1, fm_nbonopt
       ipar = ipar + 1
       fm_kb(i) = par(ipar)
       IF (.NOT.fm_fit_fc_only) THEN
          ipar = ipar + 1
          fm_b0(i) = par(ipar)
       ENDIF
    ENDDO
    ! Angles
    DO i = 1, fm_nangopt
       ipar = ipar + 1
       fm_kt(i) = par(ipar)
       IF (.NOT.fm_fit_fc_only) THEN
          ipar = ipar + 1
          IF (.NOT.fm_amber) THEN
             fm_t0(i) = TANH(par(ipar))
          ENDIF
       ENDIF
    ENDDO
    ! Impropers
    DO i = 1, fm_nimpopt
       ipar = ipar + 1
       fm_kq(i) = par(ipar)
    ENDDO
    ! Dihedrals
    DO i = 1, fm_ndihopt
       ipar = ipar + 1
       fm_kp(i) = par(ipar)
    ENDDO

    IF (ipar.NE.npar) CALL stopgm('fm_kfit_par_to_fmarray',&
         ' ipar .ne. npar',& 
         __LINE__,__FILE__)

    RETURN
  END SUBROUTINE fm_kfit_par_to_fmarray
  ! 
  ! ==================================================================
  SUBROUTINE fm_kfit_writeback(npar, par)
    ! ==--------------------------------------------------------------==
    ! Write optimized parameters back into GROMOS topology
    ! nr. of parameters and parameters themselves
    INTEGER                                  :: npar
    REAL(real_8)                             :: par(npar)

#if defined (__GROMOS)
    INCLUDE 'gromos.h'
    ! Variables
    ! loop counters, temporaries, etc.
    INTEGER :: i

    ! Bonds
    DO i = 1, fm_nbonopt
       IF (fm_amber) THEN
          fm_kb(i) = fm_kb(i) / (2._real_8*fm_b0(i)**2)
       ENDIF
       b0(fm_bonopt(i)) = bohr_nm * fm_b0(i)
       cb(fm_bonopt(i)) = kb_au2grm * fm_kb(i)
    ENDDO
    ! Angles
    DO i = 1, fm_nangopt
       t0(fm_angopt(i)) = fm_t0(i)
       IF (fm_amber) THEN
          ct(fm_angopt(i)) = au_kjm * fm_kt(i) / SIN(dacos(fm_t0(i)))**2
       ELSE
          ct(fm_angopt(i)) = au_kjm * fm_kt(i)
       ENDIF
    ENDDO
    ! Impropers
    DO i = 1, fm_nimpopt
       cq(fm_impopt(i)) = au_kjm * fm_kq(i)
    ENDDO
    ! Dihedrals
    DO i = 1, fm_ndihopt
       cp(fm_dihopt(i)) = au_kjm * fm_kp(i)
    ENDDO

    RETURN
#endif
  END SUBROUTINE fm_kfit_writeback
  ! 
  ! ==================================================================
  SUBROUTINE fm_dist(r, rij, ri, rj)
    ! ==--------------------------------------------------------------==
    ! returns the distance r and distance vector rij between
    ! two points ri, rj (with rij = ri - rj)
    ! ==================================================================
    ! output: distance, and normalized distance vector
    REAL(real_8)                             :: r, rij(3), ri(3), rj(3)

! input: coordinates of points i, j

    rij = ri - rj
    r = SQRT(DOT_PRODUCT(rij, rij))
    RETURN
  END SUBROUTINE fm_dist
  ! 
  ! ==================================================================
  SUBROUTINE fm_crossproduct(z, x, y)
    ! ==--------------------------------------------------------------==
    ! compute the cross product of two (three-dimensional!) vectors x,y
    ! ==================================================================
    ! output: cross product z
    REAL(real_8)                             :: z(3), x(3), y(3)

! input: two vectors ri, rj

    z(1) = x(2)*y(3) - x(3)*y(2)
    z(2) = x(3)*y(1) - x(1)*y(3)
    z(3) = x(1)*y(2) - x(2)*y(1)
    RETURN
  END SUBROUTINE fm_crossproduct
  ! 
  ! ==================================================================
  SUBROUTINE fm_cosder(dcmp, cmp, cp, m)
    ! ==--------------------------------------------------------------==
    ! For computing dihedral forces, we need cos(m*phi), and the derivative
    ! cos(m*phi)/cos(phi), where phi is the dihedral angle and m the
    ! multiplicity (see GROMOS bible, page II-27)
    ! ==================================================================
    ! output: derivate of cos(m*phi) wrt. cos(phi), and cos(m*phi)
    REAL(real_8)                             :: dcmp, cmp, cp
    INTEGER                                  :: m

! input: cos(phi), multiplicity

    IF (m.EQ.0) THEN
       cmp = 1._real_8
       dcmp = 0._real_8
    ELSE IF (m.EQ.1) THEN
       cmp = cp
       dcmp = 1._real_8
    ELSE IF (m.EQ.2) THEN
       cmp = 2._real_8*cp**2 - 1._real_8
       dcmp = 4._real_8 * cp
    ELSE IF (m.EQ.3) THEN
       cmp = 4._real_8*cp**3 - 3._real_8*cp
       dcmp = 12._real_8*cp**2 - 3._real_8
    ELSE IF (m.EQ.4) THEN
       cmp = 8._real_8*cp**4 - 8._real_8*cp**2 + 1._real_8
       dcmp = 32._real_8*cp**3 - 16._real_8*cp
    ELSE IF (m.EQ.5) THEN
       cmp = 16._real_8*cp**5 - 20._real_8*cp**3 + 5._real_8*cp
       dcmp = 80._real_8*cp**4 - 60._real_8*cp**2 + 5._real_8
    ELSE IF (m.EQ.6) THEN
       cmp = 32._real_8*cp**6 - 48._real_8*cp**4 + 18._real_8*cp**2 - 1._real_8
       dcmp = 192._real_8*cp**5 - 192._real_8*cp**3 + 36._real_8*cp
    ELSE
       IF (paral%io_parent)&
            WRITE(6,*) "=========== ERROR ! =============="
       IF (paral%io_parent)&
            WRITE(6,*) " Illegal value for m: ", m
       CALL stopgm('fm_cosder', 'Illegal value for multiplicity',& 
            __LINE__,__FILE__)
    ENDIF

    RETURN
  END SUBROUTINE fm_cosder
  ! 
  ! ==================================================================
  SUBROUTINE fm_kfit_rms_per_atom(nframes, nr_qm, f_fit_err,&
       f_qm_ref)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nframes, nr_qm
    REAL(real_8)                             :: f_fit_err(3, nr_qm, nframes), &
                                                f_qm_ref(3, nr_qm, nframes)

    INTEGER                                  :: iframe, iqm
    REAL(real_8)                             :: f_ref(nr_qm), &
                                                rms_per_atom(4, nr_qm)

! Variables
! rms per atom
! length of ref force vector
! loop counters, etc.

    IF (paral%io_parent)&
         WRITE(6,'(A)') '  fm  computing RMS per atom'
    ! compute RMS
    CALL zeroing(rms_per_atom)!, 4*nr_qm)
    CALL zeroing(f_ref)!, nr_qm)

    DO iframe = 1, nframes
       DO iqm = 1, nr_qm
          rms_per_atom(1:3,iqm) = rms_per_atom(1:3,iqm)&
               + f_fit_err(:,iqm,iframe)**2
          rms_per_atom(4,iqm) = rms_per_atom(4,iqm)&
               + DOT_PRODUCT(f_fit_err(:,iqm,iframe), &
               f_fit_err(:,iqm,iframe))
          f_ref(iqm) = f_ref(iqm)&
               + DOT_PRODUCT(f_qm_ref(:, iqm, iframe),&
               f_qm_ref(:, iqm, iframe))
       ENDDO
    ENDDO
    f_ref = SQRT(rms_per_atom(4, :) / f_ref)
    rms_per_atom = SQRT(rms_per_atom / nframes)

    ! print out
    IF (paral%io_parent)&
         WRITE(6,'(A)') '  fm  ========================================='
    IF (paral%io_parent)&
         WRITE(6,'(A)') '  fm  RMS per atom'
    IF (paral%io_parent)&
         WRITE(6,'(A)') '  fm  ========================================='
    DO iqm = 1, nr_qm
       IF (paral%io_parent)&
            WRITE(6,222) rms_per_atom(:, iqm), f_ref(iqm)
    ENDDO
    IF (paral%io_parent)&
         WRITE(6,'(A)') '  fm  ========================================='

    RETURN

222 FORMAT('  fm  ', 5(e9.3, 1x))
  END SUBROUTINE fm_kfit_rms_per_atom
  ! 


END MODULE forcematch_kfit_utils
