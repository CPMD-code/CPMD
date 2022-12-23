MODULE mm_forcematch_utils
  USE andp,                            ONLY: rin0,&
                                             rmix,&
                                             rout0
  USE coor,                            ONLY: fion,&
                                             tau0,&
                                             taup,&
                                             velp
  USE copot_utils,                     ONLY: copot
  USE efld,                            ONLY: extf,&
                                             textfld
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_app,&
                                             fo_old
  USE forcematch,                      ONLY: &
       fm_charges_only, fm_compute_rms, fm_compute_sp, fm_compute_traj, &
       fm_covforces_ref_filename, fm_fit_charges, fm_forces_ref_filename, &
       fm_frame_indx, fm_outtopo_filename, fm_read_forces, &
       fm_ref_traj_filename, fm_ref_traj_stride, fm_rinit_wf, fm_sp_restart
  USE forcematch_kfit_utils,           ONLY: fm_fit_covalent,&
                                             fm_kfit_compute_frms,&
                                             fm_kfit_rms_per_atom
  USE forcematch_qfit_utils,           ONLY: fm_get_size_of_qfit_arrays,&
                                             fm_update_topo_charges,&
                                             qfit
  USE forcematch_utils,                ONLY: count_nframes_ref_traj,&
                                             fm_combine_capping,&
                                             fm_read_ref_forces,&
                                             fm_restore_exclusions,&
                                             fm_setup_cpmd_equiv,&
                                             fm_writeout_topo,&
                                             sp_rest_count_skip_ref_traj
  USE forcep_utils,                    ONLY: rhoe_psi_size
  USE initrun_driver,                  ONLY: initrun
  USE ions,                            ONLY: ions1
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE machine,                         ONLY: m_flush,&
                                             m_walltime
  USE md_driver,                       ONLY: give_scr_mddiag
  USE mm_dim_utils,                    ONLY: mm_dim
  USE mm_dimmod,                       ONLY: mm_go_mm,&
                                             mm_go_qm,&
                                             mm_revert,&
                                             mmdim,&
                                             nam,&
                                             naq,&
                                             nat_grm
  USE mm_forces_diag_utils,            ONLY: mm_forces_diag
  USE mm_input,                        ONLY: cgrest_i,&
                                             clc,&
                                             lqmmm
  USE mm_mddiag_utils,                 ONLY: give_scr_mm_mddiag
  USE mm_parallel,                     ONLY: gparai,&
                                             gparal
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sync
  USE nlcc,                            ONLY: corel
  USE parac,                           ONLY: parai,&
                                             paral
  USE phfac_utils,                     ONLY: phfac
  USE poin,                            ONLY: rhoo
  USE prmem_utils,                     ONLY: prmem
  USE pslo,                            ONLY: pslo_com
  USE rinitwf_driver,                  ONLY: rinitwf
  USE ropt,                            ONLY: iteropt,&
                                             ropt_mod
  USE setirec_utils,                   ONLY: read_irec
  USE spin,                            ONLY: clsd
  USE store_types,                     ONLY: cprint
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             fpar,&
                                             maxsys,&
                                             nacc,&
                                             ncpw,&
                                             nkpt,&
                                             restf
  USE wrener_utils,                    ONLY: wrener
  USE wrgeo_utils,                     ONLY: wrgeof
  USE zeroing_utils,                   ONLY: zeroing
#include "sizeof.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mm_forcematch

CONTAINS

  ! ==================================================================
  SUBROUTINE mm_forcematch(c0,cm,c1,c2,sc0,vpp,gamx,gamy)
    COMPLEX(real_8) :: c0(nkpt%ngwk,crge%n,nkpt%nkpts), cm(*), c1(:,:,:,:), &
      c2(nkpt%ngwk,crge%n), sc0(ncpw%ngw,crge%n)
    REAL(real_8)                             :: vpp(*), gamx(*), gamy(*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'mm_forcematch'
    INTEGER, PARAMETER                       :: fm_cfu = 121, fm_iou = 120

    CHARACTER(len=10)                        :: prch
    CHARACTER(len=30)                        :: tag
    CHARACTER(len=500)                       :: dummyline
    COMPLEX(real_8), ALLOCATABLE             :: psi(:,:)
    INTEGER :: i, ia, iats, ierr, ifcalc, iframe, il_psi_1d, il_psi_2d, &
      il_rhoe_1d, il_rhoe_2d, in, iqm, irec(100), is, k, l, lscr, nframes, &
      nframes_ref_traj, nmin, nmm, nnx, no_frames_extract, nr_nn_max, nr_qm, &
      nr_qm_prev, nstate, nx, out, skipped_frames_reftraj, stride_indx
    INTEGER, ALLOCATABLE                     :: equiv(:), nr_nn(:), qm_indx(:)
    LOGICAL                                  :: ex, ferror, oldstatus, &
                                                status, statusdummy
    REAL(real_8)                             :: econs, mm_epot, rms, rms_rel, &
                                                tcpu, time1, time2
    REAL(real_8), ALLOCATABLE :: eigv(:,:), f_fit_err(:,:,:), f_nn(:,:,:), &
      f_qm(:,:,:), f_qm_ref(:,:,:), gr_tau(:,:), mm_FION(:,:,:), &
      q_hirshfeld(:), q_opt(:), r_nn(:,:,:), r_qm(:,:,:), rhoe(:,:), rinp(:), &
      rm1(:), scr(:), taui(:,:,:), tauio(:,:), taur(:,:,:), v_nn(:,:)

! Variables
! _FM[

#if defined (__GROMOS)
    CALL mm_dim(mm_go_mm,status)
    ALLOCATE(gr_tau(3,ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ex=.FALSE.
    ! 
    time1 =m_walltime()
    CALL mm_dim(mm_go_mm,oldstatus)
    ! 
    nstate=crge%n
    IF (paral%qmnode) THEN
       IF (textfld)THEN
          ALLOCATE(extf(fpar%kr1*fpar%kr2s*fpar%kr3s),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          CALL zeroing(extf)!,kr1*kr2s*kr3s)
       ENDIF
       ! Memory for densities
       nnx=fpar%nnr1*clsd%nlsd
       ALLOCATE(rin0(fpar%nnr1,clsd%nlsd),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(rout0(fpar%nnr1,clsd%nlsd),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(rmix(fpar%nnr1,clsd%nlsd),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(rm1(nnx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(rinp(nnx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       rhoo => rin0
       CALL zeroing(rin0)!,nnx)
       CALL zeroing(rout0)!,nnx)
       CALL zeroing(rmix)!,nnx)
       nacc = 22
       iteropt%nfi  = 0
       ropt_mod%modens=.FALSE.
       ropt_mod%engpri=.FALSE.
       ropt_mod%calste=cntl%tpres
       ! 
       IF (cntl%tqmmm.AND.pslo_com%tivan) THEN
          CALL stopgm('MM_FORCEMATCH','Qmmm WITH VANDERBILT UNSUPPORTED',& 
               __LINE__,__FILE__)
       ENDIF
       ALLOCATE(eigv(crge%n,clsd%nlsd*nkpt%nkpts),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ! SCR ALLOCATION AND PARTITION (SCRATCH ARRAY).
       CALL rhoe_psi_size(il_psi_1d=il_psi_1d,il_psi_2d=il_psi_2d,&
            il_rhoe_1d=il_rhoe_1d,il_rhoe_2d=il_rhoe_2d)
       ALLOCATE(rhoe(il_rhoe_1d,il_rhoe_2d),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       nmm=1
       ALLOCATE(psi(il_psi_1d*nmm,il_psi_2d),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       IF (lqmmm%qmmm) THEN
          CALL give_scr_mm_mddiag(lscr,tag)
       ELSE
          CALL give_scr_mddiag(lscr,tag)
       ENDIF
       ALLOCATE(scr(lscr),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ELSE
       nmin=10
       ALLOCATE(eigv(crge%n,nmin/crge%n),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL give_scr_mddiag(lscr,tag)
       ALLOCATE(scr(lscr),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF ! qmnode
    ! ==--------------------------------------------------------------==
    restf%nfnow=1
    ! ==--------------------------------------------------------------==
    CALL mm_dim(mm_go_mm,statusdummy)
    ALLOCATE(fion(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(taui(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(taur(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(taup(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(tauio(3,ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    CALL zeroing(fion)!,3*maxsys%nax*maxsys%nsx)
    CALL zeroing(taui)!,3*maxsys%nax*maxsys%nsx)
    CALL zeroing(taup)!,3*maxsys%nax*maxsys%nsx)

    ! _FM]
    ! check if we need to compute reference forces
    IF (fm_compute_traj) THEN
       IF (paral%parent) THEN
          IF (paral%io_parent) THEN
             WRITE(6,'(A)') '  fm  =============== ERROR ================='
             WRITE(6,'(A)') '  fm  Computing a reference trajectory is '
             WRITE(6,'(A)') '  fm  not implemented! '
             WRITE(6,'(A)') '  fm  =============== ERROR ================='
          ENDIF
          CALL stopgm('mm_forcematch','Computing reference trajectory'&
               // 'not yet implemented! ',& 
               __LINE__,__FILE__)
       ENDIF
    ELSEIF (fm_read_forces) THEN
       GOTO 666
    ELSEIF (fm_compute_sp) THEN
       CONTINUE
    ELSE
       CALL stopgm('mm_forcematch', 'I do not know what to do',& 
            __LINE__,__FILE__)
    ENDIF
    ! _FM]

    ! show sizes of the individual MPI-threads.
    ! use barriers so that the output is not garbled.
    CALL mp_sync(parai%qmmmgrp)
    IF (lqmmm%qmmm_verbose) THEN
       DO k=0,parai%nproc
          IF (gparal%mmparent) THEN
             IF (paral%io_parent)&
                  WRITE(prch,'(a10)') ' MMPARENT '
          ELSE
             IF (paral%io_parent)&
                  WRITE(prch,'(a5,i5)')' QM-N',parai%me
          ENDIF
          IF (k.EQ.parai%me) CALL prmem(prch)
          CALL mp_sync(parai%qmmmgrp)
       ENDDO
    ENDIF
    IF (paral%io_parent.AND.cgrest_i%n_cg.GE.1)&
         CALL fileopen(40,'restrain.out',fo_app,ferror)
    ! 
    ! ==--------------------------------------------------------------==
    CALL mm_dim(mm_go_mm,statusdummy)
    IF (paral%qmnode ) THEN
       IF (cprint%iprint_step.EQ.0) cprint%iprint_step=cnti%nomore+1
       ! ==--------------------------------------------------------------==
       ! Set IREC to restart file.
       CALL read_irec(irec)
       ! INITIALIZATION OF WAVEFUNCTION AND COORDINATES
       CALL initrun(irec,c0,c2,sc0,rhoe,psi,eigv)
       ! mb-bug
       IF (paral%parent) CALL dcopy(3*maxsys%nax*maxsys%nsx,taup,1,taui,1)
       ! mb-bug
       ! ==--------------------------------------------------------------==
       ! Dont symmetrize density
       cntl%tsymrho=.FALSE.
       ! ..Make sure TKFULL=.TRUE
       IF (tkpts%tkpnt.AND.(.NOT.tkpts%tkfull)) THEN
          IF (paral%io_parent) THEN
             WRITE(6,*)&
                  ' ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
             WRITE(6,*)&
                  ' WARNING! USE KPOINTS FULL  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
             WRITE(6,*)&
                  ' ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          ENDIF
       ENDIF

       ! ==--------------------------------------------------------------==
       ! 28/3/03 added
       CALL mm_dim(mm_go_qm,statusdummy)
       CALL phfac(tau0)
       IF (corel%tinlc) CALL copot(rhoe,psi,ropt_mod%calste)
       CALL mm_dim(mm_go_mm,statusdummy)
       CALL zeroing(cm(1:ncpw%ngw*nstate))!,ngw*nstate)
       ! 28/3/03
       IF (clc%classical) CALL zeroing(cm(1:ncpw%ngw*nstate))!,ngw*nstate)
       ! INITIALIZE FORCES
       CALL mm_dim(mm_go_qm,statusdummy)
       IF (cntl%tdiag) THEN
          IF (cntl%tlanc) nx=1
          IF (cntl%tdavi) nx=nkpt%ngwk*cnti%ndavv*nkpt%nkpnt+1
          IF (cntl%diis)  nx=((nkpt%ngwk*crge%n+8)*cnti%mdiis*nkpt%nkpnt)/4
       ELSEIF (cntl%tsde) THEN
          nx=1
       ELSEIF (cntl%diis) THEN
          nx=(nkpt%ngwk*crge%n+8)*cnti%mdiis/2+4
       ELSEIF (cntl%pcg) THEN
          nx=1
       ENDIF
       IF (paral%parent) THEN
          IF (paral%io_parent) THEN
             WRITE(6,'(1X,64("="))')
             WRITE(6,'(1X,"==",T25,A,T64,"==")')'FORCES INITIALIZATION'
             WRITE(6,'(1X,64("="))')
          ENDIF
       ENDIF
       ifcalc=0
    ENDIF ! (qmnode)

    CALL mm_dim(mm_go_mm,statusdummy)
    IF ( paral%parent ) THEN
       CALL mm_write_gromos_coord('CRD_INI.g96',tau0,velp,maxsys%nax,maxsys%nsx)
    ENDIF

    CALL mp_sync(parai%qmmmgrp)
    ! 
    IF (paral%qmnode) THEN
       CALL mm_dim(mm_go_mm,statusdummy)
       CALL dcopy(nnx,rin0,1,rm1,1)
       IF (paral%parent) THEN
          IF (paral%io_parent) THEN
             WRITE(6,'(1X,64("="))')
             WRITE(6,'(1X,"==",T20,A,T64,"==")')'END OF FORCES INITIALIZATION'
             WRITE(6,'(1X,64("="),/)')
          ENDIF
       ENDIF
       CALL mm_dim(mm_go_mm,statusdummy)
       IF (paral%parent) THEN
          ! _FM[
          econs=ener_com%etot
          ! _FM]
          time2 =m_walltime()
          tcpu = (time2 - time1)*0.001_real_8
          IF (paral%io_parent)&
               WRITE(6,'(A,T50,F8.2,A8)') ' TIME FOR INITIALIZATION:',&
               tcpu,' SECONDS'
          IF (paral%io_parent)&
               WRITE(6,'(//,1X,64("="))')
       ENDIF
    ENDIF ! (qmnode)
    ! ==--------------------------------------------------------------==
    ! == END INITIALIZATION                                           ==
    ! ==--------------------------------------------------------------==
    ! _FM[ preparation for the SP
    IF (paral%parent) THEN
       IF (paral%io_parent) THEN
          WRITE(6,'(A)')&
               '  fm fm fm fm fm fm fm fm fm fm fm fm fm fm fm fm fm fm fm fm fm'
          WRITE(6,'(A)') '  fm '
          WRITE(6,'(A,i8)') '  fm TOTAL NR OF ATOMS ',mmdim%natm
          WRITE(6,'(A,i8)') '  fm NR QM ATOMS ',mmdim%natq
       ENDIF
    ENDIF

    IF (paral%parent) THEN
       IF (fm_ref_traj_stride.LE.0)&
            CALL stopgm('mm_forcematch', 'TRAJECTORY STRIDE IS <= 0! ',& 
            __LINE__,__FILE__)
       CALL m_flush(6)
       ! count number of frames in ref_traj file
       nframes_ref_traj=0
       CALL count_nframes_ref_traj(mmdim%natm,nframes_ref_traj)
       IF (paral%io_parent) THEN
          WRITE(6,'(A)') '  fm AFTER COUNTING FRAMES'
          WRITE(6,'(A,i8)') '  fm no frames in traj ref file: ',nframes_ref_traj
          WRITE(6,'(A,i8)') '  fm applying stride: ',fm_ref_traj_stride
       END IF
       ! 
       ! for a fm restart: determine the frame index where to start
       ! and calculate the number of frames to extract accordingly
       ! 
       IF (fm_sp_restart) THEN

          IF (paral%io_parent) THEN
             WRITE(6,'(A)') '  fm '
             WRITE(6,'(A)') '  fm restart of the SP calculations'
             WRITE(6,'(A)') '  fm FM_REF* files will be appended'
             WRITE(6,'(A)') '  fm '
          END IF
          ! this subroutine goes through all the frames already present in 
          ! the ref_forces file
          ! for each frame it tries to find the corresponding frame in the 
          ! ref_traj file and counts the number of frames we can skip when
          ! reading the ref_traj file before we can begin with the SP_RESART
          ! run
          skipped_frames_reftraj=0
          CALL sp_rest_count_skip_ref_traj(mmdim%natq,mmdim%natm,&
               skipped_frames_reftraj)
          no_frames_extract=&
               INT((nframes_ref_traj-skipped_frames_reftraj)/&
               fm_ref_traj_stride)
       ELSE
          no_frames_extract=INT(nframes_ref_traj/fm_ref_traj_stride)
       ENDIF! fm_sp_restart
       IF (paral%io_parent) THEN
          WRITE(6,'(A,i8)')'  fm extracting total number of frames for SPs: ',no_frames_extract
          WRITE(6,'(A)') '  fm '
          WRITE(6,'(A)')&
               '  fm fm fm fm fm fm fm fm fm fm fm fm fm fm fm fm fm fm fm fm fm'
       END IF
       IF (no_frames_extract .LE. 0) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A)') '  fm not sufficient frames to extract'
          GOTO 311
       ENDIF
    ENDIF ! parent
    ! check whether number of frames, number of lines and
    ! number of atoms are consistent          
    ! not implemented yet
    CALL m_flush(6)
    CALL mp_sync(parai%qmmmgrp)
    CALL mp_bcast(no_frames_extract,parai%source,parai%allgrp)
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            CALL fileopen(fm_iou,fm_ref_traj_filename,fo_old,ferror)
       IF ( fm_sp_restart) THEN
          ! 
          ! go to the end of the last frame that had been written to 
          ! the fm_ref_forces file
          ! 
          DO iframe=1,skipped_frames_reftraj
             DO iats=1,mmdim%natm
                IF (paral%io_parent)&
                     READ(fm_iou,'(A500)') dummyline
                ! 
                ! jump over "NEW DATA" marks
                ! 
                IF (INDEX(dummyline,'NEW DATA') .NE. 0) THEN
                   IF (paral%io_parent)&
                        WRITE(6,'(A,A,i8)')&
                        '  fm (reading skipped frames)',&
                        ' NEW DATA mark after frame ',iframe-1
                   ! actually, the NEW DATA mark should not be inside
                   ! one frame
                   IF (iats .GT. 1) THEN
                      IF (paral%io_parent)&
                           WRITE(6,'(A)')&
                           '  fm NEW DATA mark in the middle of frame'
                      GOTO 311
                   ENDIF
                   ! just read an additional line            
                   IF (paral%io_parent)&
                        READ(fm_iou,'(A500)',END=311,err=311) dummyline
                ENDIF
             ENDDO! NATm
          ENDDO! skipped frames
       ENDIF! fm sp restart
    ENDIF ! parent
    ! ==--------------------------------------------------------------==
    ! == LOOP OVER THE FRAMES TO COMPUTE SPs                          ==
    ! ==--------------------------------------------------------------==
    ! 
    ! this is going to hold the labels for the individual frames
    ! extracted from the reference trajectory file
    ! in principle only parent has to know about this
    fm_frame_indx = 0
    DO l=1,no_frames_extract
       IF (paral%parent) THEN
          ! jump over the stride        
          DO stride_indx=1,fm_ref_traj_stride-1
             DO iats=1,mmdim%natm
                IF (paral%io_parent)&
                     READ(fm_iou,'(A500)',END=311,err=311) dummyline
                ! 
                ! jump over "NEW DATA" marks
                ! 
                ! write(6,*) 'before NEW DATA check'
                IF (INDEX(dummyline,'NEW DATA') .NE. 0) THEN
                   IF (paral%io_parent)&
                        WRITE(6,'(A,i8)')&
                        '  fm (reading stride) NEW DATA mark after frame ',l-1
                   IF (paral%io_parent)&
                        WRITE(6,'(A,i8)') '  fm stride index ',stride_indx
                   ! actually, the NEW DATA mark should not be inside
                   ! one frame
                   IF (iats .GT. 1) THEN
                      IF (paral%io_parent)&
                           WRITE(6,'(A)')&
                           '  fm NEW DATA mark in the middle of frame'
                      GOTO 311
                   ENDIF
                   ! just read an additional line            
                   IF (paral%io_parent)&
                        READ(fm_iou,'(A500)',END=311,err=311) dummyline
                ENDIF
             ENDDO! no of atoms
          ENDDO! stride
          ! now read the frame
          DO iats=1,mmdim%natm
             IF (paral%io_parent)&
                  READ(fm_iou,'(A500)',END=311,err=311) dummyline
             ! 
             ! jump over "NEW DATA" marks
             ! 
             IF (INDEX(dummyline,'NEW DATA') .NE. 0) THEN
                IF (paral%io_parent)&
                     WRITE(6,'(A,i8)') '  fm found NEW DATA mark in frame ',l
                ! actually, the NEW DATA mark should not be inside
                ! one frame
                IF (iats .GT. 1) THEN
                   IF (paral%io_parent)&
                        WRITE(6,'(A)')&
                        '  fm NEW DATA mark in the middle of frame'
                   GOTO 311
                ENDIF
                ! just read an additional line            
                IF (paral%io_parent)&
                     READ(fm_iou,'(A500)',END=311,err=311) dummyline
             ENDIF
             IF (paral%io_parent)&
                  READ(dummyline,*,err=311) fm_frame_indx,&
                  (gr_tau(k,iats),k=1,3)
             IF (iats .EQ. 1) THEN
                IF (paral%io_parent) THEN
                   WRITE(6,*)
                   WRITE(6,'(A)') '  fm '
                   WRITE(6,'(A,i8)') '  fm  frame number: ', l
                   WRITE(6,'(A,i8)') '  fm  frame index:  ', fm_frame_indx
                   WRITE(6,'(A)') '  fm '
                ENDIF
             ENDIF
             ! write(6,*) (gr_tau(k,iats),k=1,3)             
          ENDDO! loop over all the atoms in the current frame
          ! read trajectory file as it's being written by printp
          ! in gromos ordering
          ! convert to cpmd ordering
          i=0
          DO is=1,mmdim%nspm
             DO ia=1,NAm(is)
                i=i+1
                DO k=1,3
                   tau0(k,ia,is)=gr_tau(k,nat_grm(i))
                ENDDO
             ENDDO
          ENDDO
       ENDIF! PARENT
       CALL mm_dim(mm_go_mm,statusdummy)
       CALL mp_sync(parai%qmmmgrp)
       ! 
       ! compute reference forces
       IF (paral%parent) THEN
          IF (paral%io_parent) THEN
             WRITE(6,*)
             WRITE(6,'(A)') '  fm '
             WRITE(6,'(A)') '  fm  computing reference forces'
             WRITE(6,'(A)') '  fm '
          ENDIF
       ENDIF
       CALL mp_bcast(tau0,SIZE(tau0),parai%source,parai%qmmmgrp)
       CALL mm_translate_qmmm(tau0,c0,cm,nstate)

       IF (paral%qmnode) THEN
          CALL mm_dim(mm_go_qm,statusdummy)
          ! We re-initialize the wavefunction, because atomic wavefunctions seem to
          ! provide a better initial guess than the wfn from the previous step
          ! (the timestep between frames might be very large)
          ! Let the user decide
          IF (fm_rinit_wf) THEN
             IF (paral%parent) THEN
                IF (paral%io_parent) THEN
                   WRITE(6,*)
                   WRITE(6,'(A)') '  fm '
                   WRITE(6, '(A)')'  fm  Re-initializing wavefunction'
                   WRITE(6,'(A)') '  fm '
                   WRITE(6,*)
                ENDIF
             ENDIF
             CALL rinitwf(c0,cm,sc0,nstate,tau0,taup,rhoe,psi)
          ENDIF
          CALL phfac(tau0)
          ! no nonlinear core correction
          IF (corel%tinlc) CALL copot(rhoe,psi,ropt_mod%calste)
          CALL mm_dim(mm_go_mm,statusdummy)
       ENDIF! (qmnode)
       CALL mp_sync(parai%qmmmgrp)
       IF (paral%parent) THEN
          IF (paral%io_parent) THEN
             WRITE(6,*)
             WRITE(6,'(A)') '  fm '
             WRITE(6, '(A)') '  fm  calling forces'
             WRITE(6,'(A)') '  fm '
             WRITE(6,*)
          ENDIF
       ENDIF
       CALL mm_forces_diag(crge%n,c0(:,:,1),c1,c2,cm,sc0,cm(nx),vpp,eigv,&
            rhoe,psi,&
            tau0,velp,tau0,fion,ifcalc,&
            irec,.TRUE.,.FALSE.)
       CALL mm_dim(mm_go_mm,statusdummy)
       IF (paral%parent) THEN
          CALL wrener
          CALL wrgeof(tau0,fion)
          IF (paral%io_parent) THEN
             WRITE(6,'(A)') '  fm  ----------------------------------'
             WRITE(6,'(A, i8)') '  fm  Total nr. of iterations: ',&
                  IFCALC
             WRITE(6,'(A)') '  fm  ----------------------------------'
             WRITE(6,*)
          END IF
          ! 
          ! here append the forces of the new frame to the output file
          ! append NN atoms
          ! CALL mm_dim(mm_go_mm,statusdummy)
          IF (paral%io_parent) THEN
             CALL fileopen(63,fm_forces_ref_filename,fo_app,ferror)
             WRITE(63,'(5X,I9,1x,I9)') fm_frame_indx,mmdim%natq
          END IF
          i=0
          DO is=1,mmdim%nspq
             DO ia=1,NAq(is)
                i=i+1
                IF (paral%io_parent)&
                     WRITE(63,'(1X,I3,1X,3(E16.8),1X,3(E16.8))')&
                     i,&
                     (tau0(k,ia,is),k=1,3),(fion(k,ia,is),k=1,3)
             ENDDO
          ENDDO
          IF (paral%io_parent)&
               CALL fileclose(63)
       ENDIF! parent
       ! 
       ! append NN atoms: this is done by mm_elstat!
       ! 
    ENDDO ! loop over reference frames
    IF ((paral%parent).AND.paral%io_parent)&
         CALL fileclose(fm_iou)
    CALL mp_sync(parai%allgrp)
    ! ==--------------------------------------------------------------==
    ! == END LOOP OVER THE FRAMES                                     ==
    ! ==--------------------------------------------------------------==
    ! 
666 CONTINUE
    ! _FM[
    ! FORCE MATCHING -------------------------------------------
    IF (paral%parent) THEN
       IF (fm_read_forces.OR.fm_compute_sp) THEN
          ! determine nr of frames, max nr of NN atoms and nr of QM atoms from
          ! TRAJECTORY_PIP
          IF (fm_fit_charges) THEN
             IF (paral%io_parent) THEN
                WRITE(6, '(A)')&
                     '  fm fm fm fm fm fm fm fm fm fm fm fm fm fm'
                WRITE(6, '(A)') '  fm  Reading values for charge fitting '&
                     // 'from file'
             END IF
             CALL fm_get_size_of_qfit_arrays(nframes, nr_nn_max, nr_qm)
          ELSE
             IF (paral%io_parent) THEN
                WRITE(6, *)
                WRITE(6,'(A)') '  fm Reading NR QM atoms, frames from file'
                ! determine nr of frames and nr of QM atoms from FM_REF_FORCES
                CALL fileopen(63,fm_forces_ref_filename,fo_old,ferror)
             END IF
             nframes=0
             nr_qm_prev=0
             DO
                IF (paral%io_parent)&
                     READ(63,*,END=3192,err=318) fm_frame_indx,nr_qm
                ! write(6,*) fm_frame_indx
                IF ((nframes > 0) .AND. (nr_qm /= nr_qm_prev)) THEN
                   CALL stopgm('mm_forcematch', 'inconsistent nr of qm'&
                        //'atoms in FORCES FILE',& 
                        __LINE__,__FILE__)
                ENDIF
                nr_qm_prev=nr_qm
                DO iats=1,nr_qm
                   IF (paral%io_parent)&
                        READ(63,'(A500)',END=318,err=318) dummyline
                ENDDO
                nframes=nframes+1
             ENDDO
3192         CONTINUE
             nr_nn_max=1
             IF (paral%io_parent) THEN
                WRITE(6,'(A,i8)') ' fm  nr of qm atoms: ',nr_qm
                WRITE(6,'(A,i8)') ' fm  nr of frames: ',nframes
             END IF
             ! call fm_get_nframes_natoms(nframes,nr_qm)
          ENDIF
          ALLOCATE(nr_nn(nframes),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
          ALLOCATE(r_qm(3, nr_qm, nframes),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
          ALLOCATE(f_qm(3, nr_qm, nframes),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
          ALLOCATE(f_qm_ref(3, nr_qm, nframes),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
          ALLOCATE(r_nn(3, nr_nn_max, nframes),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
          ALLOCATE(v_nn(nr_nn_max, nframes),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
          ALLOCATE(f_nn(3, nr_nn_max, nframes),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
          ALLOCATE(q_hirshfeld(nr_qm),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
          ALLOCATE(q_opt(nr_qm),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
          ALLOCATE(equiv(nr_qm),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
          ALLOCATE(qm_indx(nr_qm),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
          ! set up equivalencies in CPMD numbering
          CALL fm_setup_cpmd_equiv(nr_qm)
          ! now fill the arrays
          CALL fm_read_ref_forces(nframes, nr_nn_max, nr_nn, r_nn,&
               v_nn, f_nn, nr_qm, r_qm, f_qm,&
               q_hirshfeld)
          ! fit charges if requested
          IF (fm_fit_charges) THEN
             CALL qfit(nframes, nr_nn_max, nr_nn, r_nn, v_nn, f_nn,&
                  nr_qm, r_qm, q_hirshfeld, equiv, q_opt)
          ENDIF! fit charges
          ! for now assume all QM atoms are fitted
          DO k = 1, nr_qm
             qm_indx(k) = k
          ENDDO
          IF (fm_fit_charges) THEN
             CALL fm_update_topo_charges(nr_qm, qm_indx, q_opt)
          ENDIF
          IF (fm_charges_only) THEN
             CALL fm_writeout_topo('QFIT.top')
             CALL stopgm('MM_FORCEMATCH','stop after charge fitting',& 
                  __LINE__,__FILE__)
          ENDIF
       ENDIF! read forces/compute sp
    ENDIF ! parent
    IF (paral%parent .AND. (.NOT. fm_fit_charges))  THEN
       IF (paral%io_parent)&
            WRITE(6,'(A)') ' fm  have just skipped the charge fitting'
    ENDIF
    CALL mp_sync(parai%allgrp)

    IF (clc%classical) THEN
       IF (paral%io_parent)&
            WRITE(6,*)&
            '  with classical forces doesnot make sense to continue'
       CALL stopgm('mm_forcematch', 'no force matching in a classical'&
            //' run',& 
            __LINE__,__FILE__)
    ENDIF
    ! now loop over reference frames again to get forces resulting from
    ! QM-MM interactions and non-bonded QM-QM interactions
    ! We get these by computing the forces according to the classical
    ! force-field with optimized charges (remember that all force-constants
    ! for bonded QM-QM interaction have been set to zero!)
    ! 
    ! First, we need to reset exclusions (QM/MM changes the topology such
    ! that all QM atoms are mutually excluded...)
    IF ((.NOT.clc%classical) .AND. (paral%parent.OR.gparal%mmnode)) THEN
       CALL fm_restore_exclusions
    ENDIF
    IF (paral%parent) THEN
       ! 
       ! We also combine forces on capping hydrogens with those on the capped atom
       CALL fm_combine_capping(nframes, nr_qm, f_qm)
       ! make a copy of reference forces to compute RMS of total force (we
       ! subtract non-bonded forces later)
       CALL dcopy(3*nr_qm*nframes, f_qm, 1, f_qm_ref, 1)
       ! 
       IF (paral%io_parent) THEN
          WRITE(6,'(A)')
          WRITE(6,'(A)') '  fm  ========================================'
          WRITE(6,'(A)') '  fm  Will now loop over reference frames again'
          WRITE(6,'(A)') '  fm  to get non-bonded forces'
          WRITE(6,'(A)') '  fm  ========================================'
       ENDIF
    ENDIF

    ! some preparation first
    CALL mp_sync(parai%qmmmgrp)
    clc%classical = .TRUE.
    CALL mp_bcast_byte(clc, size_in_bytes_of(clc),parai%source,parai%allgrp)
    CALL mp_bcast(nframes,parai%source,parai%allgrp)

    IF (paral%parent) THEN
       IF (paral%io_parent)&
            CALL fileopen(fm_iou,fm_ref_traj_filename,fo_old,ferror)
       IF (ferror) CALL stopgm('mm_forcematch', 'Could not open file '&
            // fm_ref_traj_filename,& 
            __LINE__,__FILE__)
       IF (paral%io_parent)&
            CALL fileopen(fm_cfu,fm_covforces_ref_filename,fo_app,ferror)
    ENDIF ! parent

    CALL mm_dim(mm_go_mm,statusdummy)
    fm_frame_indx = 0
    ALLOCATE(mm_FION(3, maxsys%nax, maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    ! now loop over the frames
    ! FIXME: add a check that the frame nrs. (fm_frame_indx) are consistent
    DO iframe=1,nframes
       IF (paral%parent) THEN
          ! jump over the stride        
          DO stride_indx=1,fm_ref_traj_stride-1
             DO iats=1,mmdim%natm
                IF (paral%io_parent)&
                     READ(fm_iou,'(A500)',END=311,err=311) dummyline
                ! 
                ! jump over "NEW DATA" marks
                ! 
                IF (INDEX(dummyline,'NEW DATA') .NE. 0) THEN
                   IF (paral%io_parent)&
                        WRITE(6,'(A,i8)')&
                        '  fm (reading stride) NEW DATA mark after frame ',l-1
                   IF (paral%io_parent)&
                        WRITE(6,'(A,i8)') '  fm stride index ',stride_indx
                   ! actually, the NEW DATA mark should not be inside
                   ! one frame
                   IF (iats .GT. 1) THEN
                      IF (paral%io_parent)&
                           WRITE(6,'(A)')&
                           '  fm NEW DATA mark in the middle of frame'
                      GOTO 311
                   ENDIF
                   ! just read an additional line            
                   IF (paral%io_parent)&
                        READ(fm_iou,'(A500)',END=311,err=311) dummyline
                ENDIF
             ENDDO! no of atoms
          ENDDO! stride
          DO iats=1,mmdim%natm
             IF (paral%io_parent)&
                  READ(fm_iou,'(A500)',END=311,err=311) dummyline
             ! 
             ! jump over "NEW DATA" marks
             ! 
             IF (INDEX(dummyline,'NEW DATA') .NE. 0) THEN
                IF (paral%io_parent)&
                     WRITE(6,'(A,i8)')&
                     '  fm found NEW DATA mark in frame ',iframe
                ! actually, the NEW DATA mark should not be inside
                ! one frame
                IF (iats .GT. 1) THEN
                   IF (paral%io_parent)&
                        WRITE(6,'(A)')&
                        '  fm NEW DATA mark in the middle of frame'
                   GOTO 311
                ENDIF
                ! just read an additional line            
                IF (paral%io_parent)&
                     READ(fm_iou,'(A500)',END=311,err=311) dummyline
             ENDIF
             IF (paral%io_parent)&
                  READ(dummyline,*,err=311) fm_frame_indx,&
                  (gr_tau(k,iats),k=1,3)
             IF (iats .EQ. 1) THEN
                IF (paral%io_parent) THEN
                   WRITE(6,*)
                   WRITE(6,'(A)') '  fm '
                   WRITE(6,'(A,i9)') '  fm  frame number: ', iframe
                   WRITE(6,'(A,i9)') '  fm  frame index:  ', fm_frame_indx
                   WRITE(6,'(A)') '  fm '
                ENDIF
             ENDIF
          ENDDO
          ! convert to cpmd ordering
          i=0
          DO is=1,mmdim%nspm
             DO ia=1,NAm(is)
                i=i+1
                DO k=1,3
                   tau0(k,ia,is)=gr_tau(k,nat_grm(i))
                ENDDO
             ENDDO
          ENDDO
       ENDIF! PARENT
       CALL mp_sync(parai%qmmmgrp)
       CALL mm_dim(mm_go_mm,statusdummy)
       CALL mp_bcast(tau0,SIZE(tau0),parai%source,parai%qmmmgrp)
       ! mmnode computes the classical forces
       IF (gparal%mmnode) THEN
          mm_epot = 0._real_8
          CALL zeroing(mm_FION)!,3*maxsys%nax*maxsys%nsx)
          CALL mm_force(maxsys%nax,maxsys%nsx,tau0,mm_FION,mm_epot)
       ENDIF
       CALL mp_bcast(mm_FION,SIZE(mm_FION),gparai%mmsource,parai%qmmmgrp)
       CALL mp_sync(parai%qmmmgrp)

       ! subtract classical forces from f_qm, we do this on parent 
       ! mm_FION contains actually gradients, while we have forces in f_qm,
       ! therefore to subtract, we do: f_qm - (-mm_FION) = f_qm + mm_FION
       IF (paral%parent) THEN
          iqm = 0
          DO is = 1, mmdim%nspq
             DO ia = 1, NAq(is)
                iqm = iqm + 1
                ! f_qm(:,iqm,iframe) = f_qm(:,iqm,iframe) + mm_FION(:,ia,is)
                ! I think mm_FION is actually really forces...
                f_qm(:,iqm,iframe) = f_qm(:,iqm,iframe) - mm_FION(:,ia,is)
             ENDDO
          ENDDO
          IF (paral%io_parent)&
               WRITE(fm_cfu, *) fm_frame_indx
          DO iqm = 1, nr_qm
             IF (paral%io_parent)&
                  WRITE(fm_cfu, '(1X,I3,1X,3(E16.8),1X,3(E16.8))')&
                  iqm, r_qm(:,iqm,iframe), f_qm(:,iqm,iframe)
          ENDDO
       ENDIF! parent
    ENDDO ! loop over frames
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            CALL fileclose(fm_iou)
       IF (paral%io_parent)&
            CALL fileclose(fm_cfu)
    ENDIF
    IF (paral%parent) THEN
       IF (paral%io_parent) THEN
          WRITE(6,'(A)') '  fm  ========================================'
          WRITE(6,'(A)') '  fm   Done with classical loop! '
          WRITE(6,'(A)') '  fm  ========================================'
       ENDIF
    ENDIF
    ! 
    ! _DEBUG[
    IF (paral%parent) CALL fm_writeout_topo('afterclassical.top')
    ! _DEBUG]
    ! 
    ! now fit covalent parameters
    IF (paral%parent) CALL fm_fit_covalent(nframes, nr_qm, r_qm, f_qm)
    ! 
    ! DONE FORCE MATCHING -------------------------------------------
    ! 
    ! write out new topology
    IF (paral%parent) CALL fm_writeout_topo(fm_outtopo_filename)
    ! _FM]
    ! 
    ! loop again over frames to compute RMS of total forces, if requested
    CALL mp_sync(parai%allgrp)
    IF (fm_compute_rms) THEN
       IF (paral%parent) THEN
          IF (paral%io_parent) THEN
             WRITE(6,*)
             WRITE(6,'(A)') '  fm '
             WRITE(6,'(A)') '  fm  Computing RMS of total force'
             WRITE(6,'(A)') '  fm '
          END IF
          IF (paral%io_parent)&
               CALL fileopen(fm_iou,fm_ref_traj_filename,fo_old,ferror)
          IF (ferror) THEN
             CALL stopgm('mm_forcematch', 'Could not open file'&
                  // fm_ref_traj_filename,& 
                  __LINE__,__FILE__)
          ENDIF
          ALLOCATE(f_fit_err(3, nr_qm, nframes),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
       ENDIF

       DO iframe=1,nframes
          IF (paral%parent) THEN
             DO stride_indx=1,fm_ref_traj_stride-1
                DO iats=1,mmdim%natm
                   IF (paral%io_parent)&
                        READ(fm_iou,'(A500)',END=311,err=311) dummyline
                   ! 
                   ! jump over "NEW DATA" marks
                   ! 
                   ! write(6,*) 'before NEW DATA check'
                   IF (INDEX(dummyline,'NEW DATA') .NE. 0) THEN
                      IF (paral%io_parent)&
                           WRITE(6,'(A,i8)')&
                           '  fm (reading stride) NEW DATA mark after frame ',l-1
                      IF (paral%io_parent)&
                           WRITE(6,'(A,i8)') '  fm stride index ',stride_indx
                      ! actually, the NEW DATA mark should not be inside
                      ! one frame
                      IF (iats .GT. 1) THEN
                         IF (paral%io_parent)&
                              WRITE(6,'(A)')&
                              '  fm NEW DATA mark in the middle of frame'
                         GOTO 311
                      ENDIF
                      ! just read an additional line            
                      IF (paral%io_parent)&
                           READ(fm_iou,'(A500)',END=311,err=311) dummyline
                   ENDIF
                ENDDO! no of atoms
             ENDDO! stride
             ! write(6,*) 'have read the stride ',stride_indx
             ! now read the frame
             DO iats=1,mmdim%natm
                IF (paral%io_parent)&
                     READ(fm_iou,'(A500)',END=311,err=311) dummyline
                ! 
                ! jump over "NEW DATA" marks
                ! 
                IF (INDEX(dummyline,'NEW DATA') .NE. 0) THEN
                   IF (paral%io_parent)&
                        WRITE(6,'(A,i8)')&
                        '  fm found NEW DATA mark in frame ',iframe
                   ! actually, the NEW DATA mark should not be inside
                   ! one frame
                   IF (iats .GT. 1) THEN
                      IF (paral%io_parent)&
                           WRITE(6,'(A)')&
                           '  fm NEW DATA mark in the middle of frame'
                      GOTO 311
                   ENDIF
                   ! just read an additional line            
                   IF (paral%io_parent)&
                        READ(fm_iou,'(A500)',END=311,err=311) dummyline
                ENDIF
                IF (paral%io_parent)&
                     READ(dummyline,*,err=311) fm_frame_indx,&
                     (gr_tau(k,iats),k=1,3)
                IF (iats .EQ. 1) THEN
                   IF (paral%io_parent) THEN
                      WRITE(6,*)
                      WRITE(6,'(A)') '  fm '
                      WRITE(6,'(A,i9)') '  fm  frame number: ', iframe
                      WRITE(6,'(A,i9)') '  fm  frame index:  ', fm_frame_indx
                      WRITE(6,'(A)') '  fm '
                   ENDIF
                ENDIF
             ENDDO
             ! convert to cpmd ordering
             i=0
             DO is=1,mmdim%nspm
                DO ia=1,NAm(is)
                   i=i+1
                   DO k=1,3
                      tau0(k,ia,is)=gr_tau(k,nat_grm(i))
                   ENDDO
                ENDDO
             ENDDO
          ENDIF! PARENT
          CALL mp_sync(parai%qmmmgrp)
          CALL mm_dim(mm_go_mm,statusdummy)
          CALL mp_bcast(tau0,SIZE(tau0),parai%source,parai%qmmmgrp)
          ! mmnode computes the classical forces
          IF (gparal%mmnode) THEN
             mm_epot = 0._real_8
             CALL zeroing(mm_FION)!,3*maxsys%nax*maxsys%nsx)
             CALL mm_force(maxsys%nax,maxsys%nsx,tau0,mm_FION,mm_epot)
          ENDIF
          ! mmnode might not be parent, so broadcast mm_FION from mmnode
          CALL mp_sync(parai%qmmmgrp)
          CALL mp_bcast(mm_FION,SIZE(mm_FION),gparai%mmsource,parai%qmmmgrp)
          CALL mp_sync(parai%qmmmgrp)

          IF (paral%parent) THEN
             iqm = 0
             DO is = 1, mmdim%nspq
                DO ia = 1, NAq(is)
                   iqm = iqm + 1
                   f_fit_err(:,iqm,iframe) = mm_FION(:,ia,is)&
                        - f_qm_ref(:,iqm, iframe)
                ENDDO
             ENDDO
          ENDIF! parent
       ENDDO! loop over frames
       IF (paral%parent) THEN
          CALL fm_kfit_compute_frms(rms, rms_rel, 3*nr_qm*nframes,&
               f_fit_err, f_qm_ref)
          IF (paral%io_parent) &
               WRITE(6,'(A,2f12.6)') '  fm  Total force RMS (abs. and rel.): ', rms, rms_rel
          CALL fm_kfit_rms_per_atom(nframes, nr_qm, f_fit_err, f_qm_ref)
       ENDIF
    ENDIF ! parent
    ! 
    DEALLOCATE(mm_FION,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    IF (paral%parent) THEN
       DEALLOCATE(nr_nn,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(r_qm,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(f_qm,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(f_qm_ref,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(r_nn,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(v_nn,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(f_nn,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(q_hirshfeld,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(q_opt,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(equiv,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(qm_indx,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    ENDIF  ! parent
    ! 
    ! call mp_sync(ALLGRP)
    ! 
    ! ==--------------------------------------------------------------==
    IF (paral%qmnode) THEN
       DEALLOCATE(rin0,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(rout0,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(rmix,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(rm1,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(rinp,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(eigv,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(rhoe,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(psi,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(scr,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(fion,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(taui,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(taur,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(taup,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(tauio,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       IF (textfld) DEALLOCATE(extf,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       ! ==--------------------------------------------------------------==
    ELSE
       DEALLOCATE(eigv,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF ! (qmnode)a
    DEALLOCATE(gr_tau,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL mm_dim(mm_revert,oldstatus)
    IF ((paral%parent.AND.cgrest_i%n_cg.GT.1).AND.paral%io_parent)&
         CALL fileclose(40)
    IF ((paral%parent).AND.paral%io_parent)&
         CALL fileclose(fm_iou)
    CALL mp_sync(parai%qmmmgrp)
#endif
    RETURN
    ! 
    ! ==--------------------------------------------------------------==
    ! error while reading file fm_ref_traj_filename
311 CONTINUE
    IF (paral%parent) THEN
       CALL stopgm('fm_forcematch',&
            'error while reading ref traj file '&
            //fm_ref_traj_filename,& 
            __LINE__,__FILE__)
    ENDIF

318 CONTINUE
    IF (paral%parent) THEN
       CALL stopgm('fm_forcematch',&
            'error while reading ref forces file '&
            //fm_forces_ref_filename,& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! 
    RETURN
  END SUBROUTINE mm_forcematch
  ! ==================================================================

END MODULE mm_forcematch_utils
