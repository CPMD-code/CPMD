MODULE forcematch_utils
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
       fileopen
  USE fileopenmod,                     ONLY: fo_def,&
       fo_old
  USE forcematch,                      ONLY: &
       fm_cap, fm_capping, fm_chj_ref_filename, fm_equiv, fm_equivalence, &
       fm_fit_charges, fm_forces_ref_filename, fm_iscap, fm_ncap, fm_nreq, &
       fm_pip_ref_filename, fm_ref_traj_filename, fm_save_ine, fm_save_ine14, &
       fm_save_jsne, fm_save_jsne14, fm_save_kne, fm_save_kne14, grm_equiv, &
       qm_equiv
  USE kinds,                           ONLY: real_8
  USE machine,                         ONLY: m_flush
  USE mm_dimmod,                       ONLY: mmdim,&
       nat_cpmd,&
       nat_grm,&
       solsolv
  USE parac,                           ONLY: paral
  USE readsr_utils,                    ONLY: xstring
#if defined (__GROMOS)
  USE coordsz
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: fm_setup_grm_equiv
  PUBLIC :: fm_setup_cpmd_equiv
  PUBLIC :: fm_qmatoms_in_grmorder
  PUBLIC :: fm_read_ref_forces
  PUBLIC :: count_nframes_ref_traj
  PUBLIC :: sp_rest_count_skip_ref_traj
  PUBLIC :: fm_restore_exclusions
  PUBLIC :: fm_writeout_topo
  PUBLIC :: fm_combine_capping

CONTAINS

  ! Author: Patrik Mauer
  ! Revision: Ivano Tavernelli (Dec. 2010)
  ! ==================================================================
  ! Routines for forcematching
  ! ==================================================================
  SUBROUTINE fm_setup_grm_equiv
    ! ==--------------------------------------------------------------==
    ! setup equivalencies with gromos numbering
    ! ==================================================================
    ! Variables
    ! nr of QM atoms
    CHARACTER(*), PARAMETER :: procedureN = 'fm_setup_grm_equiv'

    INTEGER                                  :: ieq, ierr, imm, jmm, nr_qm

! iterators, temporaries, etc.
! 

    IF (paral%io_parent)&
         WRITE(6,*)&
         ' fm fm fm fm fm fm fm fm fm fm fm fm fm fm fm fm fm fm fm fm fm'
    IF (paral%io_parent)&
         WRITE(6,*) ' fm '
    IF (paral%io_parent)&
         WRITE(6,'(A)',advance='NO') '  fm  Setting up equivalencies..'
    nr_qm = mmdim%natq

    ALLOCATE(grm_equiv(solsolv%nrpt),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    DO imm = 1, solsolv%nrpt
       grm_equiv(imm) = imm
    ENDDO
    IF (fm_equivalence) THEN
       DO ieq = 1, fm_nreq
          imm = fm_equiv(1, ieq)
          jmm = fm_equiv(2, ieq)
          grm_equiv(jmm) = imm
       ENDDO
    ENDIF
    IF (paral%io_parent)&
         WRITE(6,'(A)') '..done! '
    RETURN
  END SUBROUTINE fm_setup_grm_equiv
  ! 
  ! ==================================================================
  SUBROUTINE fm_setup_cpmd_equiv(nr_qm)
    ! ==--------------------------------------------------------------==
    ! setup equivalencies with CPMD numbering, also handle capping
    ! ==================================================================
    INTEGER                                  :: nr_qm

    CHARACTER(*), PARAMETER :: procedureN = 'fm_setup_cpmd_equiv'

    INTEGER                                  :: icap, ieq, ierr, iqm

! Variables
! iterators etc.
! setup global qm_uniq

    ALLOCATE(qm_equiv(nr_qm),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    IF (fm_equivalence) THEN
       DO iqm = 1, nr_qm
          ieq = NAT_cpmd(grm_equiv(NAT_grm(iqm)))
          IF (ieq.GT.nr_qm) THEN
             IF (paral%io_parent)&
                  WRITE(6, *) '========= ERROR! ================'
             IF (paral%io_parent)&
                  WRITE(6, *) ' iqm, ieq, nr_qm:'
             IF (paral%io_parent)&
                  WRITE(6, *) iqm, ieq, nr_qm
             CALL stopgm('fm_setup_cpmd_equiv','ieq .gt. nr_qm, you'&
                  // 'probably specified equivalence for an MM atom',& 
                  __LINE__,__FILE__)
          ELSE
             qm_equiv(iqm) = ieq
          ENDIF
       ENDDO
    ELSE
       DO iqm = 1, nr_qm
          qm_equiv(iqm) = iqm
       ENDDO
    ENDIF
    ! now capping
    ALLOCATE(fm_iscap(nr_qm),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    fm_iscap = .FALSE.
    IF (fm_capping) THEN
       DO icap = 1, fm_ncap
          fm_cap(1, icap) = NAT_cpmd(fm_cap(1, icap))
          fm_cap(2, icap) = NAT_cpmd(fm_cap(2, icap))
          fm_iscap(fm_cap(1, icap)) = .TRUE.
          ! DEBUG[
          IF (paral%io_parent)&
               WRITE(6,*) ' QM capping: ', fm_cap(:, icap)
          ! DEBUG]
       ENDDO
    ENDIF

    RETURN
  END SUBROUTINE fm_setup_cpmd_equiv
  ! 
  ! ==================================================================
  SUBROUTINE fm_qmatoms_in_grmorder(nr_qm, grm_qm)
    ! ==--------------------------------------------------------------==
    ! create an index array that has the CPMD number of QM atoms in gromos
    ! order, i.e. grm_qm(1) contains the CPMD number of the QM atom with 
    ! lowest GROMOS number.
    ! We need this for reading TRAJECTORY_PIP (maybe also elsewhere)      
    ! May be obsolete now.
    ! ==================================================================
    INTEGER                                  :: nr_qm, grm_qm(nr_qm)

    INTEGER                                  :: curr_qm, iqm, jqm, lowest, &
                                                tmp_grm(nr_qm)

! Variables
! array for sorting
! iterators, temporaries, and such

    IF (nr_qm .NE. mmdim%natq) THEN
       IF (paral%io_parent)&
            WRITE(6,*) "========== ERROR ! ==============="
       IF (paral%io_parent)&
            WRITE(6,*) " internal inconsistency, nr_qm is not equal NATq"
       IF (paral%io_parent)&
            WRITE(6,*) " nr_qm: ", nr_qm
       IF (paral%io_parent)&
            WRITE(6,*) " NATq: ", mmdim%natq
       CALL stopgm('fm_qmatoms_in_grmorder','nr_qm .ne. NATq',& 
            __LINE__,__FILE__)
    ENDIF

    ! the following is not very sophisticated, but for the small nr of QM
    ! atoms will do...
    DO iqm = 1, nr_qm
       tmp_grm(iqm) = NAT_grm(iqm)
    ENDDO

    curr_qm = 0
    DO
       lowest = solsolv%nrpt + 1
       DO iqm = 1, nr_qm
          IF (tmp_grm(iqm).LE.lowest) THEN
             lowest = tmp_grm(iqm)
             jqm = iqm
          ENDIF
       ENDDO
       curr_qm = curr_qm + 1
       grm_qm(curr_qm) = jqm
       tmp_grm(jqm) = solsolv%nrpt + 1
       IF (curr_qm .EQ. nr_qm) EXIT
    ENDDO

    ! DEBUG[
    ! write(6,*) 'grm_qm: ', grm_qm
    ! DEBUG]
    RETURN
  END SUBROUTINE fm_qmatoms_in_grmorder
  ! 
  ! ==================================================================
  SUBROUTINE fm_read_ref_forces(nframes, nr_nn_max, nr_nn, r_nn,&
       v_nn, f_nn, nr_qm, r_qm, f_qm,&
       q_hirshfeld)
    ! ==--------------------------------------------------------------==
    ! nr of frames
    INTEGER                                  :: nframes, nr_nn_max, &
                                                nr_nn(nframes)
    REAL(real_8) :: r_nn(3, nr_nn_max, nframes), v_nn(nr_nn_max, nframes), &
      f_nn(3, nr_nn_max, nframes)
    INTEGER                                  :: nr_qm
    REAL(real_8)                             :: r_qm(3, nr_qm, nframes), &
                                                f_qm(3, nr_qm, nframes), &
                                                q_hirshfeld(nr_qm)

    INTEGER, PARAMETER                       :: chju = 73, forceu = 71, &
                                                pipu = 72

    CHARACTER(len=2)                         :: qm_or_mm
    INTEGER                                  :: ia, idum, ie, iframe, inn, &
                                                inr_nn, ios, iqm, jnn
    LOGICAL                                  :: ferror
    REAL(real_8)                             :: f(3), posdum(3), q_mm, rdum, &
                                                read_chj(nr_qm, nframes), v

! maximal nr. of NN atoms in a frame
! nr of NN atoms for each frame
! coordinates of NN atoms 
! potential on NN atoms 
! field on NN atoms 
! nr of QM atoms
! coordinates of QM atoms 
! forces on QM atoms 
! average hirshfeld charge on QM atoms
! Variables
! CPMD atoms in gromos order
! variables in TRAJECTORY_PIP
! error upon opening file
! string handling
! iterators and temporaries
! for reading hirshfeld charges
! (dummy) variables for file reading
! file units
! read coordinates and forces on QM atoms. 

    IF (paral%io_parent)&
         CALL fileopen(forceu, fm_forces_ref_filename, fo_old, ferror)
    IF (ferror) CALL stopgm('fm_read_forces','error opening file '&
         // fm_forces_ref_filename,& 
         __LINE__,__FILE__)
    DO iframe = 1, nframes
       IF (paral%io_parent)&
            READ(forceu, *, END=301, err=302) idum
       DO iqm = 1, nr_qm
          IF (paral%io_parent)&
               READ(forceu, *, END=302, err=302) idum,&
               r_qm(:, iqm, iframe),&
               f_qm(:, iqm, iframe)
       ENDDO
    ENDDO
    ! Manuel        
    IF (.NOT. fm_fit_charges) RETURN

    ! fill arrays with values from TRAJECTORY_PIP
    ! call fm_qmatoms_in_grmorder(nr_qm, grm_qm)
    IF (paral%io_parent)&
         CALL fileopen(pipu, fm_pip_ref_filename, fo_old, ferror)
    IF (ferror) CALL stopgm('fm_read_forces','error opening file '&
         // fm_pip_ref_filename,& 
         __LINE__,__FILE__)
    CALL xstring(fm_pip_ref_filename, ia, ie)
    IF (paral%io_parent)&
         WRITE(6,'(A)') "  fm  Read potential and field from file "&
         // fm_pip_ref_filename(ia:ie)

    DO iframe = 1, nframes
       iqm = 0
       inn = 0
       IF (paral%io_parent)&
            READ(pipu, *, END=321, err=321) inr_nn
       nr_nn(iframe) = inr_nn - nr_qm
       DO jnn = 1, inr_nn
          IF (paral%io_parent)&
               READ(pipu, *, END=321, err=321) posdum, qm_or_mm,&
               q_mm, v, rdum, rdum, rdum, f
          IF (qm_or_mm .EQ. 'QM') THEN
             iqm = iqm + 1
             ! do not read coordinates for QM atoms (we have them from FORCES file,
             ! and here they are in Gromos order...)
             ! r_qm(:, jqm, iframe) = posdum
          ELSEIF (qm_or_mm .EQ. 'MM') THEN
             inn = inn + 1
             r_nn(:, inn, iframe) = posdum
             v_nn(inn, iframe) = v
             ! the check below is necessary because TRAJECTORY_PIP contains the
             ! FORCE, not the field, we should change that...
             ! if (abs(q_mm).gt.1.e-6_real_8) then
             ! f_nn(:, inn, iframe) = f / q_mm
             ! else
             ! f_nn(:, inn, iframe) = 0._real_8
             ! endif
             ! We write out field now, so we are fine
             f_nn(:, inn, iframe) = f
          ELSE
             CALL stopgm('fm_read_ref_forces', 'Unknown identifier '&
                  // qm_or_mm,& 
                  __LINE__,__FILE__)
          ENDIF
       ENDDO
    ENDDO
    IF (paral%io_parent)&
         CALL fileclose(pipu)

    ! read Hirshfeld charges from file 
    CALL xstring(fm_chj_ref_filename, ia, ie)
    IF (paral%io_parent)&
         CALL fileopen(chju, fm_chj_ref_filename, fo_old, ferror)
    IF (ferror) CALL stopgm('fm_read_forces','error opening file '&
         // fm_chj_ref_filename(ia:ie),& 
         __LINE__,__FILE__)
    IF (paral%io_parent)&
         WRITE(6,'(A)') '  fm  Reading reference charges from file '&
         // fm_chj_ref_filename(ia:ie)
    DO iframe = 1, nframes
       IF (paral%io_parent)&
            READ(chju, *, END=323, err=322) idum
       IF (paral%io_parent)&
            READ(chju, *, END=323, err=322) read_chj(:, iframe)
    ENDDO
323 CONTINUE
    iframe = iframe - 1
    IF (paral%io_parent)&
         WRITE(6, '(A, i9)') '  fm  Nr. of frames read: ', iframe
    IF (iframe.LT.nframes) THEN
       IF (paral%io_parent)&
            WRITE(6, '(A)') '  fm  ======== WARNING ========'
       IF (paral%io_parent)&
            WRITE(6, '(A)') '  fm  not enough frames in file '&
            // fm_chj_ref_filename(ia:ie)
       IF (paral%io_parent)&
            WRITE(6, '(A)') '  fm  Make sure this is OK! '
       IF (paral%io_parent)&
            WRITE(6, '(A)') '  fm  ========================='
    ELSE! try to read one more line...
       IF (paral%io_parent)&
            READ(chju, *, iostat = ios) idum, q_hirshfeld
       IF (ios .GT. 0) THEN! error
          GOTO 322
       ELSEIF (ios .LT. 0) THEN! end of file, so OK
          IF (paral%io_parent)&
               WRITE(6, '(A)') '  fm  correct nr. of frames in file '&
               // fm_chj_ref_filename(ia:ie)
       ELSE! there are more frames..
          IF (paral%io_parent)&
               WRITE(6, '(A)') '  fm  ======== WARNING ========'
          IF (paral%io_parent)&
               WRITE(6, '(A)') '  fm  excessive nr. of frames in file '&
               // fm_chj_ref_filename(ia:ie) 
          IF (paral%io_parent)&
               WRITE(6, '(A)') '  fm  Make sure this is OK! '
          IF (paral%io_parent)&
               WRITE(6, '(A)') '  fm  ========================='
       ENDIF
    ENDIF
    IF (paral%io_parent)&
         CALL fileclose(chju)
    DO iqm = 1, nr_qm
       q_hirshfeld(iqm) = SUM(read_chj(iqm, 1:iframe)) / REAL(iframe,kind=real_8)
    ENDDO
    IF (paral%io_parent)&
         WRITE(6, '(A, f8.3)') '  fm  Total charge of QM system is ',&
         SUM(q_hirshfeld)
    RETURN

    ! premature end of file while reading fm_forces_ref_filename
301 CONTINUE
    CALL stopgm('fm_read_ref_forces', 'premature end of file'&
         // fm_forces_ref_filename,& 
         __LINE__,__FILE__)

    ! error while reading fm_forces_ref_filename
302 CONTINUE
    CALL stopgm('fm_read_ref_forces', 'error reading file'&
         // fm_forces_ref_filename,& 
         __LINE__,__FILE__)

    ! error while reading file fm_pip_ref_filename
321 CONTINUE
    CALL stopgm('fm_read_ref_forces', 'error while reading file '&
         // fm_pip_ref_filename,& 
         __LINE__,__FILE__)

    ! error while reading file fm_pip_ref_filename
322 CONTINUE
    CALL stopgm('fm_read_ref_forces', 'error while reading file '&
         // fm_chj_ref_filename,& 
         __LINE__,__FILE__)

  END SUBROUTINE fm_read_ref_forces
  ! 
  ! =================================================================
  SUBROUTINE count_nframes_ref_traj(nrats,nrframes)
    ! ==--------------------------------------------------------------==
    ! count number of frames in the reference trajectory file
    ! ==================================================================
    ! total number of atoms, loop index      
    INTEGER                                  :: nrats, nrframes

    CHARACTER(len=500)                       :: dummyline
    INTEGER                                  :: iats, iou
    LOGICAL                                  :: ferror

! total number of frames in file      

    iou=623
    nrframes=0
    ! open TRAJECTORY_REF
    ferror=.FALSE.
    IF (paral%io_parent)&
         CALL fileopen(iou,fm_ref_traj_filename,fo_old,ferror)
    IF (ferror) CALL stopgm('count_nframes_ref_traj',&
         'error opening ref traj file ' // fm_ref_traj_filename,& 
         __LINE__,__FILE__)
    DO
       DO iats=1,nrats
          IF (paral%io_parent)&
               READ(iou,'(A500)',END=310,err=311) dummyline
          ! 
          ! jump over "NEW DATA" marks
          ! 
          IF (INDEX(dummyline,'NEW DATA') .NE. 0) THEN
             ! actually, the NEW DATA mark should not be inside
             ! one frame
             IF (paral%io_parent)&
                  WRITE(6,'(A,i9,A,i9)') '  fm NEW DATA mark iats: ',iats,&
                  ' nrframe: ',nrframes
             IF (iats .GT. 1) GOTO 311
             IF (paral%io_parent)&
                  READ(iou,'(A500)',END=310,err=311) dummyline
          ENDIF
       ENDDO
       nrframes=nrframes+1
    ENDDO
    ! reached end of file
310 CONTINUE
    IF (paral%io_parent)&
         CALL fileclose(iou)
    CALL m_flush(6)
    RETURN
    ! error while reading file 
311 CONTINUE
    CALL stopgm('fm_forcematch',&
         'error while reading ref traf file ' &
         //fm_ref_traj_filename,& 
         __LINE__,__FILE__)
    RETURN
  END SUBROUTINE count_nframes_ref_traj
  ! 
  ! ==================================================================
  SUBROUTINE sp_rest_count_skip_ref_traj(nratqm,nrats,nskipped)
    ! ==--------------------------------------------------------------==
    ! in case of a restart of the single points calculations we want to
    ! check how many frames are already in the fm_ref_forces file and
    ! accordingly how many frames we have to skip in the ref_traj file
    ! before we can read out the next new frame
    ! ==================================================================
    INTEGER                                  :: nratqm, nrats, nskipped

    CHARACTER(len=500)                       :: dummyline
    INTEGER                                  :: check_ref_traj_indx, &
                                                force_frame_indx, iats
    LOGICAL                                  :: ex, ferror, ferror2

! check whether the files exist

    IF (paral%io_parent)&
         INQUIRE(file=fm_forces_ref_filename,exist=ex)
    IF (.NOT. ex) THEN
       IF (paral%io_parent)&
            WRITE (6,*)&
            '  fm ',fm_forces_ref_filename,' does not exist'
       CALL stopgm('mm_forcematch', 'file not exist',& 
            __LINE__,__FILE__)
    ENDIF
    IF (paral%io_parent)&
         INQUIRE(file=fm_chj_ref_filename,exist=ex)
    IF (.NOT. ex) THEN
       IF (paral%io_parent)&
            WRITE (6,*)&
            '  fm ',fm_chj_ref_filename,' does not exist'
       CALL stopgm('mm_forcematch', 'not exist',& 
            __LINE__,__FILE__)
    ENDIF
    IF (paral%io_parent)&
         INQUIRE(file=fm_pip_ref_filename,exist=ex)
    IF (.NOT. ex) THEN
       IF (paral%io_parent)&
            WRITE (6,*)&
            '  fm ',fm_pip_ref_filename,' does not exist'
       CALL stopgm('mm_forcematch', 'not exist',& 
            __LINE__,__FILE__)
    ENDIF
    ! go to all the header lines in FM_REF_FORCES, 
    ! one after the other
    nskipped=0
    ferror=.FALSE.
    ferror2=.FALSE.
    IF (paral%io_parent)&
         CALL fileopen(63,fm_forces_ref_filename,fo_old,ferror)
    IF (paral%io_parent)&
         CALL fileopen(64,fm_ref_traj_filename,fo_old,ferror2)
    IF (ferror) THEN
       CALL stopgm('sp_rest_count_skip_ref_traj',&
            'error opening ref forces file ' // fm_forces_ref_filename,& 
            __LINE__,__FILE__)
    ELSE IF (ferror2) THEN
       CALL stopgm('sp_rest_count_skip_ref_traj',&
            'error opening ref traj file ' // fm_ref_traj_filename,& 
            __LINE__,__FILE__)
    ENDIF
    DO
       IF (paral%io_parent)&
            READ(63,*,END=316,err=318) force_frame_indx
       ! write(6,*) fm_frame_indx
       DO iats=1,nratqm
          IF (paral%io_parent)&
               READ(63,'(A500)',END=316,err=318) dummyline
       ENDDO
       IF (paral%io_parent)&
            WRITE(6,*) '  fm found frame index in forces file',&
            force_frame_indx
       ! search for the frame indexes downstream in the TRAJECTORY_REF
       ! file (this way it is not a problem if one frame index exists
       ! more then once)           
       DO
          DO iats=1,nrats
             ! read(fm_iou,'(I)',end=318,err=318) check_ref_traj_indx
             IF (paral%io_parent)&
                  READ(64,'(A500)',END=316,err=319) dummyline
             ! 
             ! jump over "NEW DATA" marks
             ! 
             IF (INDEX(dummyline,'NEW DATA') .NE. 0) THEN
                ! actually, the NEW DATA mark should not be inside
                ! one frame
                IF (paral%io_parent)&
                     WRITE(6,*) '  fm found NEW DATA mark in iats ',iats
                IF (iats .GT. 1) GOTO 319
                IF (paral%io_parent)&
                     READ(64,'(A500)',END=316,err=319) dummyline
             ENDIF
          ENDDO! iats
          IF (paral%io_parent)&
               READ(dummyline,*,err=319) check_ref_traj_indx
          nskipped=nskipped+1
          IF (check_ref_traj_indx .EQ. force_frame_indx) THEN
             IF (paral%io_parent)&
                  WRITE(6,*) '  fm found frame index in ref traj file',&
                  check_ref_traj_indx
             EXIT
          ENDIF
       ENDDO! loop over all the lines in the ref traj file
    ENDDO
316 CONTINUE
    IF (paral%io_parent)&
         CALL fileclose(63)
    IF (paral%io_parent)&
         CALL fileclose(64)
    RETURN
318 CALL stopgm('sp_rest_count_skip_ref_traj',&
         'error while reading ref forces file ' &
         //fm_ref_traj_filename,& 
         __LINE__,__FILE__)
319 CALL stopgm('sp_rest_count_skip_ref_traj',&
         'error while reading ref traj file ' &
         //fm_ref_traj_filename,& 
         __LINE__,__FILE__)

  END SUBROUTINE sp_rest_count_skip_ref_traj
  ! 
  ! 
  ! ==================================================================
  SUBROUTINE fm_restore_exclusions
    ! ==--------------------------------------------------------------==
    ! This subroutine restores exclusions so that they are correct for a
    ! purely classical run (in QM/MM, all QM atoms are mutually excluded)
    ! ==================================================================
#if defined (__GROMOS)
    IMPLICIT NONE
    INTEGER :: ierr
    CHARACTER(*),PARAMETER::procedureN='fm_restore_exclusions'
    INCLUDE 'gromos.h'
    ! 
    ine = fm_save_INE
    kne = fm_save_KNE
    jsne = fm_save_JSNE
    ine14 = fm_save_INE14
    kne14 = fm_save_KNE14
    jsne14 = fm_save_JSNE14
    ! 
    DEALLOCATE(fm_save_INE,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(fm_save_KNE,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(fm_save_JSNE,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(fm_save_INE14,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(fm_save_KNE14,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(fm_save_JSNE14,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
#endif
    RETURN
  END SUBROUTINE fm_restore_exclusions
  ! 
  ! ==================================================================
  SUBROUTINE fm_writeout_topo(filename)
    ! ==--------------------------------------------------------------==
    ! write out the topology
    ! ==================================================================
    CHARACTER(*)                             :: filename

    INTEGER                                  :: topou
    LOGICAL                                  :: ferror

! Variables
! file unit
! error with file

#if defined (__GROMOS)
    IF (paral%io_parent)&
         WRITE(6,*) ' '
    IF (paral%io_parent)&
         WRITE(6,'(A)',advance='NO') '  fm  Writing out new topology..'
    topou = 77
    IF (paral%io_parent)&
         CALL fileopen(topou, filename, fo_def, ferror)
    IF (ferror) THEN
       CALL stopgm('fm_writeout_topo', 'Could not open file'&
            // filename,& 
            __LINE__,__FILE__)
    ENDIF
    CALL wrtopo(topou)
    IF (paral%io_parent)&
         CALL fileclose(topou)
#endif
    RETURN
  END SUBROUTINE fm_writeout_topo
  ! 
  ! ==================================================================
  SUBROUTINE fm_combine_capping(nframes, nr_qm, f_qm)
    ! ==--------------------------------------------------------------==
    ! Combine forces on capping hydrogens with those on capped atoms (for
    ! now, simply add them)
    ! ==================================================================
    ! nr of frames, nr of QM atoms
    INTEGER                                  :: nframes, nr_qm
    REAL(real_8)                             :: f_qm(3,nr_qm,nframes)

    INTEGER                                  :: ccap, hcap, icap, iframe

! forces on QM atoms
! Variables
! loop counters, temporaries, etc.

    IF (.NOT.fm_capping) RETURN

    DO iframe = 1, nframes
       DO icap = 1, fm_ncap
          hcap = fm_cap(1,icap)
          ccap = fm_cap(2,icap)
          f_qm(:,ccap,iframe) = f_qm(:,ccap,iframe)&
               + f_qm(:,hcap,iframe)
          f_qm(:,hcap,iframe) = 0._real_8
       ENDDO
    ENDDO

  END SUBROUTINE fm_combine_capping
  ! 


END MODULE forcematch_utils
