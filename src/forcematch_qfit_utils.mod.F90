MODULE forcematch_qfit_utils
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
       fileopen
  USE fileopenmod,                     ONLY: fo_old
  USE forcematch,                      ONLY: &
       fm_cap, fm_capping, fm_equivalence, fm_fix_q, fm_fixed_q_grm_indx, &
       fm_fixed_q_indx, fm_fixed_q_trgt, fm_iscap, fm_max_qw, fm_n_fix_q, &
       fm_n_qw, fm_ncap, fm_optimize_weights, fm_pip_ref_filename, fm_wf, &
       fm_wq, fm_wq_general, fm_wq_grm_indx, fm_wq_indx, fm_wq_unique, &
       fm_wtot, fm_wv, qm_equiv
  USE kinds,                           ONLY: real_8
  USE mm_dimmod,                       ONLY: nat_cpmd,&
       nat_grm
  USE parac,                           ONLY: paral
  USE readsr_utils,                    ONLY: xstring
  USE zeroing_utils,                   ONLY: zeroing
#if defined (__GROMOS)
  USE coordsz
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: fm_get_size_of_qfit_arrays
  !public :: qfit_optimize_weights
  PUBLIC :: qfit
  !public :: qfit_setup_equiv
  PUBLIC :: setup_qfit_infl_mat
  !public :: infl_mat_singleframe
  !public :: setup_qfit_target
  !public :: qfit_rms
  PUBLIC :: fm_update_topo_charges

CONTAINS

  ! Author: Patrik Mauer
  ! Revision: Ivano Tavernelli (Dec. 2010)
  ! ==================================================================
  SUBROUTINE fm_get_size_of_qfit_arrays(nframes, nr_nn_max, nr_qm)
    ! ==--------------------------------------------------------------==
    ! nr of frames, max nr. of NN atoms, and nr of QM atoms
    INTEGER                                  :: nframes, nr_nn_max, nr_qm

    INTEGER, PARAMETER                       :: pipu = 72

    CHARACTER(len=2)                         :: qm_or_mm
    INTEGER                                  :: ia, ie, iframe, inn, inr_nn, &
                                                iqm
    LOGICAL                                  :: ferror
    REAL(real_8)                             :: posdum(3)

! Variables
! QM or MM identifier
! iterators, temporaries, etc.
! dummies for reading files
! file units
! get nr of frames from TRAJECTORY_PIP

    IF (paral%io_parent)&
         CALL fileopen(pipu, fm_pip_ref_filename, fo_old, ferror)
    IF (ferror) CALL stopgm('fm_read_forces','error opening file '&
         // fm_pip_ref_filename,& 
         __LINE__,__FILE__)
    CALL xstring(fm_pip_ref_filename, ia, ie)
    IF (paral%io_parent)&
         WRITE(6,'(A)') '  fm  Determine nr of frames'
    IF (paral%io_parent)&
         WRITE(6,'(A,A,A)',advance='NO') "  fm  Reading file ",&
         fm_pip_ref_filename(ia:ie) &
         // ".."

    nr_nn_max = -1
    iframe = 0
    DO
       iqm = 0
       IF (paral%io_parent)&
            READ(pipu, *, END = 310, err = 311) inr_nn
       IF (inr_nn .GE. nr_nn_max) nr_nn_max = inr_nn
       DO inn = 1, inr_nn
          IF (paral%io_parent)&
               READ(pipu, *, END=311, err=311) posdum, qm_or_mm
          IF (qm_or_mm .EQ. 'QM') iqm = iqm + 1
       ENDDO
       IF (iframe .EQ. 0) THEN
          nr_qm = iqm
       ELSE
          IF (nr_qm .NE. iqm) THEN
             CALL stopgm('fm_read_ref_forces', 'inconsistent nr. of QM'&
                  // ' atoms in file ' // fm_pip_ref_filename,& 
                  __LINE__,__FILE__)
          ENDIF
       ENDIF
       iframe = iframe + 1
    ENDDO

    ! reached end of file
310 CONTINUE
    IF (paral%io_parent)&
         CALL fileclose(pipu)
    nframes = iframe
    IF (paral%io_parent)&
         WRITE(6, '(A)') "..done! "
    IF (paral%io_parent)&
         WRITE(6,'(A, i8)') '  fm  Nr. of frames:       ', nframes
    IF (paral%io_parent)&
         WRITE(6,'(A, i8)') '  fm  Nr. of QM atoms:     ', nr_qm
    IF (paral%io_parent)&
         WRITE(6,'(A, i8)') '  fm  Max nr. of NN atoms: ', nr_nn_max
    RETURN

    ! error while reading file nn_ref_filename
311 CONTINUE
    CALL stopgm('fm_read_ref_forces', 'error while reading file '&
         // fm_pip_ref_filename,& 
         __LINE__,__FILE__)

  END SUBROUTINE fm_get_size_of_qfit_arrays
  ! 
  ! ==================================================================
  SUBROUTINE qfit_optimize_weights(nframes, nr_nn_max, nr_nn, r_nn,&
       v_nn, f_nn, nr_qm, r_qm,&
       q_hirshfeld)
    ! ==--------------------------------------------------------------==
    ! Try to find optimal weights for the charge fitting procedure
    ! ==================================================================
    ! nr of frames
    INTEGER                                  :: nframes, nr_nn_max, &
                                                nr_nn(nframes)
    REAL(real_8) :: r_nn(3, nr_nn_max, nframes), v_nn(nr_nn_max, nframes), &
      f_nn(3, nr_nn_max, nframes)
    INTEGER                                  :: nr_qm
    REAL(real_8)                             :: r_qm(3, nr_qm, nframes), &
                                                q_hirshfeld(nr_qm)

    INTEGER                                  :: iframe
    REAL(real_8)                             :: wf, wv

! maximal nr. of NN atoms in a frame
! nr of NN atoms for each frame
! coordinates of NN atoms 
! potential on NN atoms 
! field on NN atoms 
! nr of QM atoms
! coordinates of QM atoms 
! average hirshfeld charge on QM atoms
! iterators and temporaries
! 
! write(6, '(A)') '================ ERROR ====================='
! write(6, '(A)') ' finding optimal weights for charge fitting'
! write(6, '(A)') ' is not yet implemented!'
! call STOPGM('qfit_optimize_weights', 'not yet implemented')
! 

    fm_wtot = 1.0e+09_real_8
    fm_wq = 1._real_8 / SQRT(SUM(q_hirshfeld**2))
    ! 
    wv = 0._real_8
    wf = 0._real_8
    DO iframe = 1, nframes
       wv = wv + SUM(v_nn(1:nr_nn(iframe), iframe)**2)
       wf = wf + SUM(f_nn(:, 1:nr_nn(iframe), iframe)**2)
    ENDDO
    ! 
    ! _DEBUG[
    ! write(6,*) "deb: sum(q_hirsh**2), wv, wf"
    ! write(6,*) sum(q_hirshfeld**2), wv, wf
    fm_wq = 0.01_real_8
    ! _DEBUG]
    fm_wv = 1._real_8 / SQRT(wv)
    fm_wf = 1._real_8 / SQRT(wf)
    ! 
    RETURN
  END SUBROUTINE qfit_optimize_weights
  ! 
  ! ==================================================================
  SUBROUTINE qfit(nframes, nr_nn_max, nr_nn, r_nn, v_nn, f_nn,&
       nr_qm, r_qm, q_hirshfeld, equiv, q_opt)
    ! ==--------------------------------------------------------------==
    ! qfit is the main routine for charge fitting. The optimized charges are
    ! returned in array q_opt. The fitting works by solving an
    ! overdetermined system of linear equations (SLE) in the least-squares 
    ! sense.
    ! The SLE is I * q_opt = T, where:
    ! - q_opt is a vector containing the optimized charges, with equivalencies
    ! being taken into account.
    ! - T is the "target vector", i.e. a vector containing the potential and
    ! field on the NN atoms, plus Hirshfeld (or other) charges and the total
    ! charge for restraining.
    ! - I is the "influence matrix" containing all the factors 1/rij and
    ! 1/rij**2 for the potential and field, resp.
    ! ==================================================================
    ! nr of frames
    INTEGER                                  :: nframes, nr_nn_max, &
                                                nr_nn(nframes)
    REAL(real_8) :: r_nn(3, nr_nn_max, nframes), v_nn(nr_nn_max, nframes), &
      f_nn(3, nr_nn_max, nframes)
    INTEGER                                  :: nr_qm
    REAL(real_8)                             :: r_qm(3, nr_qm, nframes), &
                                                q_hirshfeld(nr_qm)
    INTEGER                                  :: equiv(nr_qm)
    REAL(real_8)                             :: q_opt(nr_qm)

    CHARACTER(*), PARAMETER                  :: procedureN = 'qfit'

    CHARACTER(40)                            :: err_msg
    INTEGER :: debug_cpmd_indx, debug_equiv2, debug_grms_indx, &
      debug_unique_indx, i, ierr, info, iqm, lwork, n_unique, trget_size, &
      unique(nr_qm)
    INTEGER, ALLOCATABLE                     :: n_eq(:)
    REAL(real_8)                             :: debug_weight, frms, qtot, vrms
    REAL(real_8), ALLOCATABLE                :: infl_mat(:, :), &
                                                q_restrain(:), trget(:), &
                                                work(:)

! maximal nr. of NN atoms in a frame
! nr of NN atoms for each frame
! coordinates of NN atoms 
! potential on NN atoms 
! field on NN atoms 
! nr of QM atoms
! coordinates of QM atoms 
! average hirshfeld charge on QM atoms
! equivalencies
! optimized charges
! Variables
! nr of atoms a unique charge represents
! target vector and "influence" matrix
! charges to restrain to
! nr of unique charges (i.e., atom types, sort of)
! index array for unique charges/atoms
! total charge of QM system
! size of target vector
! rms of potential and field
! for LAPACK
! iterators
! Manuel: printing stuff for debugging

    IF (paral%io_parent) THEN
       WRITE(6,'(A)') '  fm'
       WRITE(6, '(A)') '  fm  Starting QFIT'
       WRITE(6, '(A, i8)') '  fm  nr. of frames: ', nframes
    ENDIF
    ! write(6, '(A)') '  fm  Weights for fitting are: '
    ! write(6, '(A, e9.3)') '  fm  WV (Potential)        = ', fm_wv
    ! write(6, '(A, e9.3)') '  fm  WF (Field)            = ', fm_wf
    ! write(6, '(A, e9.3)') '  fm  WQ (Charge restraint) = ', fm_wq
    ! write(6, '(A, e9.3)') '  fm  WTOT (Total charge)   = ', fm_wtot


    ! compute nr of unique charges
    IF (fm_equivalence .OR. fm_capping) THEN
       n_unique = 0
       DO iqm = 1, nr_qm
          IF (qm_equiv(iqm).EQ.iqm .AND. .NOT.fm_iscap(iqm)) THEN
             n_unique = n_unique + 1
          ENDIF
       ENDDO
    ELSE
       n_unique = nr_qm
       DO iqm = 1, nr_qm
          equiv(iqm) = iqm
       ENDDO
    ENDIF
    IF (paral%io_parent)&
         WRITE(6, *) ' fm  Nr. of unique charges: ', n_unique
    ALLOCATE(q_restrain(n_unique),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(n_eq(n_unique),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! call MEMORY(ip_q_restrain, n_unique, 'qfit unique charges')
    ! call MEMORY(ip_n_eq, n_unique, 'qfit n_eq')

    ! compute total charge as sum over Hirshfeld charges
    qtot = SUM(q_hirshfeld)
    ! set up equiv array from fm_equiv (read in input file), this also
    ! handles capping hydrogens,
    ! the charges that are kept fixed
    ! and the individual weights for the charge fitting
    CALL qfit_setup_equiv(nr_qm, n_unique, q_hirshfeld, q_restrain,&
         equiv, n_eq, unique)
    ! find optimal weights if requested
    IF (fm_optimize_weights) THEN
       CALL qfit_optimize_weights (nframes, nr_nn_max, nr_nn, r_nn,&
            v_nn, f_nn, nr_qm, r_qm,&
            q_hirshfeld)
    ENDIF
    IF (paral%io_parent) THEN
       WRITE(6, '(A)') '  fm  Weights for fitting are:'
       WRITE(6, '(A, e9.3)') '  fm  WV (Potential)        = ', fm_wv
       WRITE(6, '(A, e9.3)') '  fm  WF (Field)            = ', fm_wf
       WRITE(6, '(A)') '  fm  WQ'
       WRITE(6, '(2A)')&
            '  fm    cpmd         gromos   ',&
            '  unique(i)     equiv to        weight           target'
    ENDIF
    DO i=1,nr_qm
       debug_cpmd_indx=i
       debug_grms_indx=NAT_grm(debug_cpmd_indx)
       debug_unique_indx=unique(i)
       debug_equiv2=equiv(i)
       debug_weight=fm_wq_unique(debug_unique_indx)
       IF (paral%io_parent)&
            WRITE(6, '(i0,1x,i0,1x,i0,1x,i0,1x,f15.8,1x,f15.8)')&
            debug_cpmd_indx,debug_grms_indx,debug_unique_indx,debug_equiv2,&
            debug_weight,q_restrain(debug_unique_indx)
    ENDDO
    IF (paral%io_parent)&
         WRITE(6, '(A, e9.3)') '  fm  WTOT (Total charge)   = ', fm_wtot
    ! setup target vector and "influence" matrix
    ! size of target vector is 1+3=4 times the total nr. of NN atoms (1 for
    ! potential, 3 for field), plus the nr. of charges to restrain to, plus
    ! 1 for the total charge
    ! 
    IF (paral%io_parent)&
         WRITE(6, '(A, i10)') '  fm  Total nr of NN atoms:      ',&
         SUM(nr_nn)
    trget_size = 4 * SUM(nr_nn) + n_unique + 1
    IF (paral%io_parent)&
         WRITE(6, '(A, i10)') '  fm  Size of the target vector: ',&
         trget_size
    IF (paral%io_parent)&
         WRITE(6, '(A)') '  fm  Allocating target vector and influence '&
         // 'matrix'
    ALLOCATE(trget(trget_size),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(infl_mat(trget_size, n_unique),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    IF (paral%io_parent)&
         WRITE(6, '(A)', advance='NO') '  fm  Setting up target vector..'
    CALL setup_qfit_target(trget_size, trget, nframes, nr_nn_max,&
         nr_nn, v_nn, f_nn, n_unique, n_eq, &
         q_restrain, qtot)
    IF (paral%io_parent)&
         WRITE(6, '(A)') '..done! '
    IF (paral%io_parent)&
         WRITE(6, '(A)', advance='NO') '  fm  Setting up influence '&
         // 'matrix..'
    CALL setup_qfit_infl_mat(trget_size, infl_mat, nframes,&
         nr_nn_max, nr_nn, r_nn, nr_qm, r_qm, &
         n_unique, n_eq, unique)
    IF (paral%io_parent)&
         WRITE(6, '(A)') '..done! '
    ! DEBUG[
    ! write(998, *) trget
    ! do ieq=1, trget_size
    ! write(999, *) infl_mat(ieq, :)
    ! enddo
    ! DEBUG]
    ! find least squares solution to overdetermined system of linear 
    ! equations using the routine DGELS from LAPACK
    lwork = MIN(trget_size,n_unique) + 32*MAX(trget_size,n_unique)
    ALLOCATE(work(lwork),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    IF (paral%io_parent)&
         WRITE(6, '(A)', advance='NO') '  fm  Solving..'
    CALL dgels('N', trget_size, n_unique, 1, infl_mat, trget_size,&
         trget, trget_size, work, lwork, info)
    IF (info .NE. 0) THEN
       IF (paral%io_parent)&
            WRITE(err_msg, '(A, i5)') ' DGELS returned info = ', info
       CALL stopgm('QFIT', err_msg,& 
            __LINE__,__FILE__)
    ENDIF
    IF (paral%io_parent)&
         WRITE(6, '(A)') '..done! '
    DEALLOCATE(work,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)

    ! DGELS returns optimal parameters as the first elements of array trget:
    ! trget(1:n_unique)
    ! put charges back into q_opt for each QM atom
    DO iqm = 1, nr_qm
       IF (fm_iscap(iqm)) THEN
          q_opt(iqm) = 0._real_8
       ELSE
          IF (unique(iqm) .LT. 0) THEN
             IF (paral%io_parent)&
                  WRITE(6, *) '  fm Apparently we are trying to do something we'&
                  // ' should not to a capping H'
             IF (paral%io_parent)&
                  WRITE(6, *) '  fm iqm is: ', iqm
             CALL stopgm('qfit',' unique(iqm) is -1',& 
                  __LINE__,__FILE__)
          ENDIF
          q_opt(iqm) = trget(unique(iqm))
       ENDIF
    ENDDO
    ! DEBUG[
    IF (paral%io_parent)&
         WRITE(6, '(2A)') '  fm    cpmd         gromos   ',&
         '  unique(i)     equiv to        q_opt(i)'
    DO i=1,nr_qm
       debug_cpmd_indx=i
       debug_grms_indx=NAT_grm(debug_cpmd_indx)
       debug_unique_indx=unique(i)
       debug_equiv2=equiv(i)
       IF (paral%io_parent)&
            WRITE(6, '(i0,1x,i0,1x,i0,1x,i0,1x,f15.8)')&
            debug_cpmd_indx,debug_grms_indx,debug_unique_indx,&
            debug_equiv2,q_opt(i)
    ENDDO
    IF (paral%io_parent)&
         WRITE(6,*) ' fm  qtot: ', SUM(q_opt)
    ! DEBUG]
    ! 
    ! compute RMS of potential and field
    CALL qfit_rms(vrms, frms, nframes, nr_qm, nr_nn_max, nr_nn,&
         r_qm, r_nn, v_nn, f_nn, q_opt, q_hirshfeld)

    ! clean up
    ! call FREEM(ip_trget)
    ! call FREEM(ip_infl_mat)
    ! call FREEM(ip_q_restrain)
    ! call FREEM(ip_n_eq)
    DEALLOCATE(trget,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(infl_mat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(q_restrain,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(n_eq,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    RETURN
  END SUBROUTINE qfit
  ! 
  ! ==================================================================
  SUBROUTINE qfit_setup_equiv(nr_qm, n_unique, q_hirshfeld,&
       q_restrain, equiv, n_eq, unique)
    ! ==--------------------------------------------------------------==
    ! qfit_setup_equiv handles equivalencies for QFIT 
    ! it also handles capping hydrogens
    ! and the charges that are fixed during the optimization
    ! and individual weights for the charge fitting
    ! ==================================================================
    ! nr of QM atoms
    INTEGER                                  :: nr_qm, n_unique
    REAL(real_8)                             :: q_hirshfeld(nr_qm), &
                                                q_restrain(n_unique)
    INTEGER                                  :: equiv(nr_qm), n_eq(n_unique), &
                                                unique(nr_qm)

    INTEGER                                  :: icap, ieq, ifixers, iqm, &
                                                iqweight, jqm

! nr of unique charges
! average Hirshfeld charges
! unique charges to restrain to
! equivalencies
! nr of atoms a charge represents
! index from QM atom to unique atom
! Variables
! index to unique arrays
! integer         iun(nr_qm)
! iterators
! 
! _DEBUG[
! write(6,*) " now in qfit_setup_equiv"
! write(6,*) " size equiv, qm_equiv: ", size(equiv),
!             size(qm_equiv)
! _DEBUG]

    DO iqm = 1, nr_qm
       equiv(iqm) = qm_equiv(iqm)
    ENDDO
    DO ieq = 1, n_unique
       n_eq(ieq) = 0
    ENDDO
    ! 
    ! find capping H and add Hirshfeld charge to capped atom
    ! DEBUG[
    ! write(6,*) "  fm q_hirsh 0:", q_hirshfeld
    ! DEBUG]
    IF (fm_capping) THEN
       DO icap = 1, fm_ncap
          iqm = fm_cap(1, icap)
          jqm = fm_cap(2, icap)
          q_hirshfeld(jqm) = q_hirshfeld(jqm) + q_hirshfeld(iqm)
          q_hirshfeld(iqm) = 0._real_8
       ENDDO
    ENDIF
    ! DEBUG[
    IF (paral%io_parent)&
         WRITE(6,'(A)') '  fm  put capping H charge on capped atom'
    ! write(6,*) "q_hirsh 1:", q_hirshfeld
    ! DEBUG]
    ! 
    ! count nr of unique charges
    ! first find atoms that have no equivalency (i.e., are equivalent to
    ! themselves) and are not a capping H
    n_unique = 0
    DO iqm = 1, nr_qm
       IF (equiv(iqm) .EQ. iqm .AND. .NOT.fm_iscap(iqm)) THEN
          n_unique = n_unique + 1
          unique(iqm) = n_unique
          n_eq(n_unique) = n_eq(n_unique) + 1
       ENDIF
    ENDDO
    ! now atoms that are equivalent to some other
    DO iqm = 1, nr_qm
       IF (equiv(iqm) .NE. iqm .AND. .NOT.fm_iscap(iqm)) THEN
          ieq = unique(equiv(iqm))
          unique(iqm) = ieq
          n_eq(ieq) = n_eq(ieq) + 1
       ENDIF
    ENDDO
    ! set uniqe() to -1 for capping H, we can check for that or we'll get a
    ! segfault...
    IF (fm_capping) THEN
       DO iqm = 1, fm_ncap
          unique(fm_cap(1, iqm)) = -1
       ENDDO
    ENDIF
    ! 
    ! charge restraint
    CALL zeroing(q_restrain)!, n_unique)
    DO iqm = 1, nr_qm
       IF (fm_iscap(iqm)) CYCLE
       ieq = unique(iqm)
       q_restrain(ieq) = q_restrain(ieq) + q_hirshfeld(iqm)
    ENDDO
    q_restrain = q_restrain / n_eq
    ! 
    ! Manuel setup the index arrays for the individual weights and
    ! the fixed charges in cpmd ordering
    DO ifixers=1,fm_n_qw
       fm_wq_indx(ifixers)=NAT_cpmd(fm_wq_grm_indx(ifixers))
    ENDDO
    IF (fm_fix_q) THEN
       DO ifixers=1,fm_n_fix_q
          fm_fixed_q_indx(ifixers)=NAT_cpmd(fm_fixed_q_grm_indx(ifixers))
       ENDDO
       ! Manuel now overwrite the restraints for the atoms whose  
       ! charges are fixed to some value defined in the input
       DO ifixers=1,fm_n_fix_q
          iqm=fm_fixed_q_indx(ifixers)
          ieq=unique(iqm)
          q_restrain(ieq) = fm_fixed_q_trgt(ifixers)
       ENDDO
    ENDIF
    ! Manuel organize the weights on the charges
    ! fm_wq() holds the weights for the charges on each QM atom
    ! fm_wq_general holds the weight for these atoms for which 
    ! no weight was specified in the input
    ! fm_wq_unique() will hold the weights for all the unique atoms
    DO ieq=1,fm_max_qw
       fm_wq_unique(ieq)=fm_wq_general
    ENDDO
    DO iqweight=1,fm_n_qw
       iqm=fm_wq_indx(iqweight)
       ieq=unique(iqm)
       fm_wq_unique(ieq)=fm_wq(iqm)
    ENDDO
    ! DEBUG[
    ! write(6,*) " fm  n_eq:", n_eq
    ! write(6,*) " fm  uniq: ", unique
    ! write(6,*) " fm  q_restr: ", q_restrain(1:n_unique)
    IF (paral%io_parent.AND.fm_fix_q) THEN
       WRITE(6,*) " fm  wq: ",fm_wq(1:nr_qm)
       WRITE(6,*) " fm  wq_unique: ",fm_wq_unique(1:n_unique)
    ENDIF
    ! DEBUG]
    RETURN
  END SUBROUTINE qfit_setup_equiv
  ! 
  ! ==================================================================
  SUBROUTINE setup_qfit_infl_mat(trget_size, infl_mat, nframes,&
       nr_nn_max, nr_nn, r_nn, nr_qm, &
       r_qm, n_unique, n_eq, unique)
    ! ==--------------------------------------------------------------==
    ! setup_qfit_infl_mat sets up the "influence" matrix for qfit
    ! ==================================================================
    ! size of target vector
    INTEGER                                  :: trget_size, nframes, &
                                                nr_nn_max, nr_nn(nframes)
    REAL(real_8)                             :: r_nn(3, nr_nn_max, nframes)
    INTEGER                                  :: nr_qm
    REAL(real_8)                             :: r_qm(3, nr_qm, nframes)
    INTEGER                                  :: n_unique
    REAL(real_8), TARGET                     :: infl_mat(trget_size, n_unique)
    INTEGER                                  :: n_eq(n_unique), unique(nr_qm)

    INTEGER                                  :: ieq, iframe, lb, nr_nn_tot, ub
    REAL(real_8), POINTER                    :: section(:, :)

! nr of unique charges
! influence matrix
! nr of frames
! maximal nr. of NN atoms in a frame
! nr of NN atoms for each frame
! coordinates of NN atoms 
! nr of QM atoms
! coordinates of QM atoms
! nr of atoms a unique charge represents
! equivalencies
! Variables
! matrix section for a single frame
! total nr of NN atoms in all frames
! iterators and array boundaries

    CALL zeroing(infl_mat)!, trget_size * n_unique)

    ! first the factors 1/r and 1/r**2 for potential and field, resp.
    ub = 0
    DO iframe = 1, nframes
       lb = ub + 1
       ub = lb + 4*nr_nn(iframe) - 1
       section => infl_mat(lb:ub, :)
       CALL infl_mat_singleframe(nr_nn(iframe), n_unique, nr_qm,&
            section, r_nn(:, :, iframe), &
            r_qm(:, :, iframe), unique)
    ENDDO

    ! now add the Hirshfeld charges and the total charge
    nr_nn_tot = SUM(nr_nn)
    lb = 4*nr_nn_tot + 1
    DO ieq = 1, n_unique
       ! infl_mat(lb+ieq-1, ieq) = fm_wq * n_eq(ieq)
       infl_mat(lb+ieq-1, ieq) = fm_wq_unique(ieq) * n_eq(ieq)
    ENDDO
    lb = lb + n_unique
    IF (lb .NE. trget_size) THEN
       IF (paral%io_parent) THEN
          WRITE(6, *) " ERROR when setting up influence matrix"
          WRITE(6, *) " lb should be equal trget_size, but they are:"
          WRITE(6, *) "         lb = ", lb
          WRITE(6, *) " trget_size = ", trget_size
       ENDIF
       CALL stopgm('setup_qfit_infl_mat', 'lb .ne. trget_size',& 
            __LINE__,__FILE__)
    ENDIF
    infl_mat(lb, :) = fm_wtot * n_eq
    RETURN
  END SUBROUTINE setup_qfit_infl_mat
  ! 
  ! ==================================================================
  SUBROUTINE infl_mat_singleframe(nr_nn, n_unique, nr_qm, infl_mat,&
       r_nn, r_qm, unique)
    ! ==--------------------------------------------------------------==
    ! infl_mat_singleframe computes a section of the influence matrix for a
    ! single frame and with only potential and field
    ! ==================================================================
    ! nr of NN atoms
    INTEGER                                  :: nr_nn, n_unique, nr_qm
    REAL(real_8)                             :: infl_mat(4*nr_nn, n_unique), &
                                                r_nn(3, nr_nn), r_qm(3, nr_qm)
    INTEGER                                  :: unique(nr_qm)

    INTEGER                                  :: inn, iqm, jnn, jqm
    REAL(real_8)                             :: rij, rijm1, rijvec(3)

! nr of unique charges
! nr of QM atoms
! influence matrix
! coordinates of NN and QM atoms
! equivalencies
! Variables
! distance vector, distance and inverse distance (1/r)
! iterators

    DO iqm = 1, nr_qm
       IF (fm_iscap(iqm)) CYCLE
       DO inn = 1, nr_nn
          rijvec = r_nn(:, inn) - r_qm(:, iqm)
          rij = SQRT(DOT_PRODUCT(rijvec, rijvec))
          rijm1 = 1._real_8 / rij
          rijvec = rijm1 * rijvec
          jnn = 4*(inn-1) + 1
          jqm = unique(iqm)
          infl_mat(jnn, jqm) = infl_mat(jnn, jqm) + fm_wv * rijm1
          infl_mat(jnn+1:jnn+3, jqm) = infl_mat(jnn+1:jnn+3, jqm)&
               + fm_wf * rijm1**2 * rijvec
       ENDDO
    ENDDO
    RETURN
  END SUBROUTINE infl_mat_singleframe
  ! 
  ! ==================================================================
  SUBROUTINE setup_qfit_target(trget_size, trget, nframes,&
       nr_nn_max, nr_nn, v_nn, f_nn, &
       n_unique, n_eq, q_restrain, qtot)
    ! ==--------------------------------------------------------------==
    ! setup_qfit_target sets up the target vector for qfit
    ! ==================================================================
    ! size of target vector
    INTEGER                                  :: trget_size
    REAL(real_8)                             :: trget(trget_size)
    INTEGER                                  :: nframes, nr_nn_max, &
                                                nr_nn(nframes)
    REAL(real_8)                             :: v_nn(nr_nn_max, nframes), &
                                                f_nn(3, nr_nn_max, nframes)
    INTEGER                                  :: n_unique, n_eq(n_unique)
    REAL(real_8)                             :: q_restrain(n_unique), qtot

    INTEGER                                  :: iframe, inn, it

! target vector
! nr of frames
! maximal nr. of NN atoms in a frame
! nr of NN atoms for each frame
! potential on NN atoms 
! field on NN atoms 
! nr of unique charges
! nr. of atoms a unique charge represents
! charges to restrain to 
! total charge of QM system
! Variables
! iterators and array bounds

    CALL zeroing(trget)!, trget_size)

    ! fill up the target vector, first potential and field on NN atoms
    it = 1
    DO iframe = 1, nframes
       DO inn = 1, nr_nn(iframe)
          trget(it) = fm_wv * v_nn(inn, iframe)
          trget(it+1:it+3) = fm_wf * f_nn(:, inn, iframe)
          it = it + 4
       ENDDO
       IF (it.GT.trget_size) WRITE(6,*) 'iframe, tgsize, it: ', iframe,&
            trget_size, it
    ENDDO
    ! Manuel fixed charges        
    ! trget(it:it+n_unique-1) = fm_wq * q_restrain * n_eq
    trget(it:it+n_unique-1) = fm_wq_unique(1:n_unique) *&
         q_restrain * n_eq
    it = it + n_unique
    IF (it.GT.trget_size) WRITE(6,*) 'tgsize, it: ', trget_size, it
    trget(it) = fm_wtot * qtot
    IF (it .NE. trget_size) THEN
       IF (paral%io_parent) THEN
          WRITE(6, *)
          WRITE(6, *) " ========== ERROR! ============"
          WRITE(6, *)&
               " it and trget_size should be equal, but they are"
          WRITE(6, *) " it: ", it
          WRITE(6, *) " trget_size:", trget_size
       ENDIF
       CALL stopgm('setup_qfit_target',&
            'it and trget_size not equal',& 
            __LINE__,__FILE__)
    ENDIF
    RETURN
  END SUBROUTINE setup_qfit_target
  ! 
  ! ==================================================================
  SUBROUTINE qfit_rms(vrms, frms, nframes, nr_qm, nr_nn_max, nr_nn,&
       r_qm, r_nn, v_nn_ref, f_nn_ref, q, &
       q_hirshfeld)
    ! ==--------------------------------------------------------------==
    ! qfit_rms computes the RMSD of the field, potential, and Hirshfeld
    ! charges for a set of charges
    ! ==================================================================
    ! rms of potential and field
    REAL(real_8)                             :: vrms, frms
    INTEGER                                  :: nframes, nr_qm, nr_nn_max, &
                                                nr_nn(nframes)
    REAL(real_8) :: r_qm(3, nr_qm, nframes), r_nn(3, nr_nn_max, nframes), &
      v_nn_ref(nr_nn_max, nframes), f_nn_ref(3, nr_nn_max, nframes), &
      q(nr_qm), q_hirshfeld(nr_qm)

    INTEGER                                  :: iframe, imaxfframe, imaxfrms, &
                                                imaxvframe, imaxvrms, inn, &
                                                iqm, nr_nn_tot
    REAL(real_8) :: f_nn(3, nr_nn_max), fdiff(3), fref, maxfrms, maxvrms, &
      qdiff, qref, qrms, rij, rijm1, rijvec(3), rms, v_nn(nr_nn_max), vdiff, &
      vref

! nr of frames
! nr of QM atoms
! max nr. of NN atoms
! nr of NN atoms in each frame
! coordinates of QM atoms
! coordinates of NN atoms
! reference potential and field on NN atoms
! charges on QM atoms
! Hirshfeld charges 
! Variables
! potential and field with fitted charges
! difference of potential, field, and charge to reference values
! mean square reference field, potential, and charge
! distance vector, distance, and inverse distance 1/r
! maximal RMS
! total nr. of NN atoms 
! iterators
! first compute RMS of potential and field     

    vrms = 0._real_8
    vref = 0._real_8
    frms = 0._real_8
    fref = 0._real_8
    maxvrms = 0._real_8
    maxfrms = 0._real_8

    nr_nn_tot = 0
    DO iframe = 1, nframes
       CALL zeroing(v_nn)!, nr_nn_max)
       CALL zeroing(f_nn)!, 3*nr_nn_max)
       DO inn = 1, nr_nn(iframe)
          DO iqm = 1, nr_qm
             rijvec = r_nn(:, inn, iframe) - r_qm(:, iqm, iframe)
             rij = SQRT(DOT_PRODUCT(rijvec, rijvec))
             rijm1 = 1._real_8 / rij
             rijvec = rijm1 * rijvec
             v_nn(inn) = v_nn(inn) + q(iqm) * rijm1
             f_nn(:, inn) = f_nn(:, inn) + q(iqm) * rijvec * rijm1**2
          ENDDO
          vdiff = v_nn(inn) - v_nn_ref(inn, iframe)
          fdiff = f_nn(:, inn) - f_nn_ref(:, inn, iframe)
          ! vrms = vrms + vdiff**2
          rms = vdiff**2
          IF (rms.GT.maxvrms) THEN
             maxvrms = rms
             imaxvrms = inn
             imaxvframe = iframe
          ENDIF
          vrms = vrms + rms
          vref = vref + v_nn_ref(inn, iframe)**2

          ! frms = frms + dot_product(fdiff, fdiff)
          rms = DOT_PRODUCT(fdiff, fdiff)
          IF (rms.GT.maxfrms) THEN
             maxfrms = rms
             imaxfrms = inn
             imaxfframe = iframe
          ENDIF
          frms = frms + rms
          fref = fref + DOT_PRODUCT(f_nn_ref(:, inn, iframe),&
               f_nn_ref(:, inn, iframe))
       ENDDO
       nr_nn_tot = nr_nn_tot + nr_nn(iframe)
    ENDDO

    ! now Hirshfeld charges
    qrms = 0._real_8
    qref = 0._real_8
    DO iqm = 1, nr_qm
       qdiff = q(iqm) - q_hirshfeld(iqm)
       qrms = qrms + qdiff**2
       qref = qref + q_hirshfeld(iqm)**2
    ENDDO

    IF (paral%io_parent) THEN
       WRITE(6, '(A)') "  fm ---------------------------------------"
       WRITE(6, '(A)') "  fm  QFIT RMS, absolute:"
       WRITE(6, '(A, f10.5)') "  fm  Potential: ", SQRT(vrms/nr_nn_tot)
       WRITE(6, '(A, f10.5)') "  fm  Field:     ", SQRT(frms/nr_nn_tot)
       WRITE(6, '(A, f10.5)') "  fm  Charges:   ", SQRT(qrms/nr_qm)
       WRITE(6, '(A)') "  fm ---------------------------------------"
       WRITE(6, '(A)') "  fm  QFIT RMS, relative:"
       WRITE(6, '(A, f10.5)') "  fm  Potential: ", SQRT(vrms/vref)
       WRITE(6, '(A, f10.5)') "  fm  Field:     ", SQRT(frms/fref)
       WRITE(6, '(A, f10.5)') "  fm  Charges:   ", SQRT(qrms/qref)
       WRITE(6, '(A)') "  fm ---------------------------------------"
       WRITE(6, '(A)') "  fm  Largest absolute RMS:"
       WRITE(6, '(A)') "  fm               iframe   inn    rms    "
       WRITE(6, '(A,2(i6),e9.2)') "  fm  Potential:", imaxvframe,&
            imaxvrms,maxvrms
       WRITE(6, '(A,2(i6),e9.2)') "  fm  Field:    ", imaxfframe,&
            imaxfrms,maxfrms
       WRITE(6, '(A)') "  fm ---------------------------------------"
    ENDIF

    RETURN
  END SUBROUTINE qfit_rms
  ! 
  ! ==================================================================
  SUBROUTINE fm_update_topo_charges(nr_q, qm_indx, qnew)
    ! ==--------------------------------------------------------------==
    ! update the GROMOS topology with the new charges
    ! ==================================================================
    ! nr of charges to update
    INTEGER                                  :: nr_q, qm_indx(nr_q)
    REAL(real_8)                             :: qnew(nr_q)

! index array with CPMD numbers of QM atoms
! new charges

#if defined (__GROMOS)
    INCLUDE 'gromos.h'
    ! Variables
    ! indices and iterators
    INTEGER :: igrm, icpmd, iqm
    ! conversion factor for charges (GROMOS stores charges with factor
    ! sqrt(1/(PI*EPS_0)) )
    REAL(real_8) :: cgfac

    cgfac = SQRT(fpepsi)
    DO iqm = 1, nr_q
       icpmd = qm_indx(iqm)
       igrm = NAT_grm(icpmd)
       cg(igrm) = cgfac * qnew(iqm)
    ENDDO
#endif
    RETURN
  END SUBROUTINE fm_update_topo_charges
  ! 

END MODULE forcematch_qfit_utils
