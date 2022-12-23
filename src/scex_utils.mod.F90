! ==================================================================
! Provides: - Communication structure for grid redistribution with
!             scaled exact exchange (ScEX)
!           - Grid redistribution routines
!
!                             19.11.2018 - M. P. Bircher @ LCBC/EPFL
! ==================================================================
MODULE scex_utils
  USE error_handling,                  ONLY: stopgm
  USE fft,                             ONLY: inzf,&
                                             inzs,&
                                             jgw,&
                                             jhg,&
                                             llr1,&
                                             nzff,&
                                             nzfs,&
                                             qr1s, qr2s, qr3s, qr1, &
                                             lr1s, lr2s, lr3s, lr1, &
                                             lrxpl,&
                                             sp5,sp8,sp9
  USE fftnew_utils,                    ONLY: setfftn,&
                                             addfftnset
  USE hfxmod,                          ONLY: ipoolhfx
  USE isos,                            ONLY: isos1,&
                                             isos3
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_recv,&
                                             mp_send,&
                                             mp_sync
  USE parac,                           ONLY: parai,&
                                             paral
  USE system,                          ONLY: cntr
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  ! ==--------------------------------------------------------------==

  PRIVATE

  INTEGER, PARAMETER, PRIVATE      :: cpmd_output_unit = 6
  INTEGER, PARAMETER, PUBLIC       :: scex_ID_parent   = 0
  INTEGER, PUBLIC, SAVE            :: scex_ID_scaled   = 0

  REAL(real_8), PARAMETER, PUBLIC  :: scex_lambda      = 2.0_real_8

  TYPE, PRIVATE :: scex_fft_t
     INTEGER, PUBLIC :: lr1  = 0
     INTEGER, PUBLIC :: lr1s = 0
     INTEGER, PUBLIC :: lr2s = 0
     INTEGER, PUBLIC :: lr3s = 0
     INTEGER, PUBLIC :: qr1  = 0
     INTEGER, PUBLIC :: qr1s = 0
     INTEGER, PUBLIC :: qr2s = 0
     INTEGER, PUBLIC :: qr3s = 0
  END TYPE scex_fft_t

  TYPE, PRIVATE :: scex_t
     LOGICAL, PUBLIC   :: init    = .FALSE.
     INTEGER, PUBLIC   :: llr1    = 0
     INTEGER, PUBLIC   :: max_lr1 = 0
     !
     INTEGER, PRIVATE  :: quad1   = 0
     INTEGER, PRIVATE  :: quad2   = 0
     INTEGER, PRIVATE  :: quad3   = 0
     INTEGER, PRIVATE  :: nbr_mpi_groups
     !
     INTEGER, DIMENSION(:), ALLOCATABLE, PRIVATE :: my_mpi_dest
     INTEGER, DIMENSION(:), ALLOCATABLE, PRIVATE :: my_mpi_source
     INTEGER, DIMENSION(:), ALLOCATABLE, PRIVATE :: my_mpi_send_tag
     INTEGER, DIMENSION(:), ALLOCATABLE, PRIVATE :: my_mpi_recv_tag
     INTEGER, DIMENSION(:), ALLOCATABLE, PRIVATE :: grp_size
     INTEGER, DIMENSION(:), ALLOCATABLE, PRIVATE :: grp_sender
     INTEGER, DIMENSION(:), ALLOCATABLE, PRIVATE :: grp_receiver
     INTEGER, DIMENSION(:), ALLOCATABLE, PRIVATE :: grp_recv_start
     INTEGER, DIMENSION(:), ALLOCATABLE, PRIVATE :: grp_send_start
     !
     COMPLEX(real_8), DIMENSION(:,:,:), ALLOCATABLE, &
                                         PRIVATE :: scr
     COMPLEX(real_8), DIMENSION(:,:,:), ALLOCATABLE, &
                                         PRIVATE :: sent
     COMPLEX(real_8), DIMENSION(:,:,:), ALLOCATABLE, &
                                         PRIVATE :: recv
     !
     TYPE(scex_fft_t), PRIVATE                   :: parent
     TYPE(scex_fft_t), PRIVATE                   :: scaled
  CONTAINS
     PROCEDURE, PASS, PUBLIC  :: grid_init                  => scex_grid_init
     PROCEDURE, PASS, PUBLIC  :: do_density_scaling         => scex_do_density_scaling
     PROCEDURE, PASS, PUBLIC  :: undo_density_scaling       => scex_undo_density_scaling
     PROCEDURE, PASS, PUBLIC  :: start_density_scaling      => scex_start_density_scaling
     PROCEDURE, PASS, PUBLIC  :: switch_density_scaling     => scex_switch_density_scaling
     PROCEDURE, PASS, PUBLIC  :: annihilate_density_scaling => scex_annihilate_density_scaling
     PROCEDURE, PASS, PRIVATE :: comm_init                  => scex_comm_init
     PROCEDURE, PASS, PRIVATE :: comm_check                 => scex_comm_check
     PROCEDURE, PASS, PRIVATE :: settings_check             => scex_settings_check
  END TYPE scex_t

  TYPE(scex_t), PUBLIC, SAVE  :: scex

CONTAINS

  ! ==================================================================
  ! Routines for type initialisation
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE scex_grid_init(scex)
    ! ==--------------------------------------------------------------==

    CHARACTER(len=*), PARAMETER  :: procedureN = 'scex_grid_init'

    CLASS(scex_t), INTENT(inout) :: scex

    CALL addfftnset(cntr%ecut,cntr%ecut,scex_ID_scaled)
    CALL scex%settings_check()
    CALL scex%comm_init()
    CALL setfftn(scex_ID_parent)
    ipoolhfx = scex_ID_scaled

    ! ==--------------------------------------------------------------==
  END SUBROUTINE scex_grid_init
  ! ==================================================================
  SUBROUTINE scex_comm_init(scex)
    ! ==--------------------------------------------------------------==

    CHARACTER(len=*), PARAMETER  :: procedureN = 'scex_comm_init'

    CLASS(scex_t), INTENT(inout) :: scex

    INTEGER                      :: i, ierr, ip, ipp, ips, ir, &
                                    ir2, irr, k2n3, n2n3, msgid, at, irs

    INTEGER                      :: isub
    LOGICAL                      :: success

    INTEGER, DIMENSION(:), ALLOCATABLE :: plane_mpi_dest
    INTEGER, DIMENSION(:), ALLOCATABLE :: plane_mpi_source

    INTEGER                              :: qq23, id, mytag, scaled_llr1
    INTEGER                              :: source_old, dest_old, nbr_groups
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: parent_lrxpl, scaled_lrxpl

    CALL tiset(procedureN,isub)

    IF (scex%init) CALL stopgm(procedureN,'ScEX type is already initialised',&
                               __LINE__,__FILE__)

    ! Private
    !
    CALL setfftn(scex_ID_parent)
    scex%parent%lr1  = lr1 
    scex%parent%lr1s = lr1s
    scex%parent%lr2s = lr2s
    scex%parent%lr3s = lr3s
    scex%parent%qr1  = qr1 
    scex%parent%qr1s = qr1s
    scex%parent%qr2s = qr2s
    scex%parent%qr3s = qr3s
    scex%quad1       = lr1s / 4
    scex%quad2       = lr2s / 4
    scex%quad3       = lr3s / 4
    ALLOCATE(parent_lrxpl,source=lrxpl,stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN,'Error in allocating from source',&
                               __LINE__,__FILE__)
    CALL setfftn(scex_ID_scaled)
    scex%scaled%lr1  = lr1 
    scex%scaled%lr1s = lr1s
    scex%scaled%lr2s = lr2s
    scex%scaled%lr3s = lr3s
    scex%scaled%qr1  = qr1 
    scex%scaled%qr1s = qr1s
    scex%scaled%qr2s = qr2s
    scex%scaled%qr3s = qr3s
    scaled_llr1      = llr1
    ALLOCATE(scaled_lrxpl,source=lrxpl,stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN,'Error in allocating from source',&
                               __LINE__,__FILE__)


    CALL setfftn(scex_ID_parent)
    ALLOCATE(scex%my_mpi_dest(scex%parent%lr1),&
             scex%my_mpi_source(scex%scaled%lr1),&
             scex%my_mpi_send_tag(scex%parent%lr1),&
             scex%my_mpi_recv_tag(scex%scaled%lr1),&
             plane_mpi_dest(scex%parent%lr1s),&
             plane_mpi_source(scex%parent%lr1s),&
             scex%grp_size(parai%nproc),&
             scex%grp_sender(parai%nproc),&
             scex%grp_receiver(parai%nproc),&
             scex%grp_recv_start(parai%nproc),&
             scex%grp_send_start(parai%nproc),&
             stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN,'Allocation problem', __LINE__,__FILE__)

    scex%my_mpi_dest(:)     = -1
    scex%my_mpi_source(:)   = -1
    scex%my_mpi_send_tag(:) = -1
    scex%my_mpi_recv_tag(:) = -1
    plane_mpi_dest(:)       = -1
    plane_mpi_source(:)     = -1
    scex%grp_size(:)        = -1
    scex%grp_sender(:)      = -1
    scex%grp_receiver(:)    = -1
    scex%grp_recv_start(:)  = -1
    scex%grp_send_start(:)  = -1


    !
    ! GLOBAL info
    !
    DO ir=1,scex%scaled%lr1s
       !
       ! What is the destination of each plane?
       !
       DO ipp=1,parai%nproc
          IF (ir >= scaled_lrxpl(ipp-1,1) .AND. ir <= scaled_lrxpl(ipp-1,2)) THEN
             plane_mpi_dest(ir)  = ipp - 1
             EXIT 
          ENDIF
       ENDDO
       !
       ! Who is the sender of each plane?
       ! 
       irs = ir + scex%quad1
       DO ipp=1,parai%nproc
          IF (irs >= parent_lrxpl(ipp-1,1) .AND. irs <= parent_lrxpl(ipp-1,2)) THEN
             plane_mpi_source(ir) = ipp - 1
             EXIT 
          ENDIF
       ENDDO 
    ENDDO

    !
    ! Figure out the number of plane groups
    !
    source_old       = -1
    dest_old         = -1

    scex%grp_size(:) = 0
    nbr_groups       = 0

    DO ir=1,scex%scaled%lr1s
       !CALL mp_sync(parai%allgrp)
       IF (source_old /= plane_mpi_source(ir) .OR. & 
           dest_old   /= plane_mpi_dest(ir) ) THEN
          nbr_groups = nbr_groups + 1
          !
          ! Who are the group pairs?
          !
          scex%grp_sender(nbr_groups)     = plane_mpi_source(ir)
          scex%grp_receiver(nbr_groups)   = plane_mpi_dest(ir)
          !
          ! On the source, which grid index is my starting point? 
          !
          IF (parai%me == plane_mpi_source(ir)) THEN 
             irs = ir - parent_lrxpl(parai%me,1) + 1 + scex%quad1
             scex%grp_send_start(nbr_groups) = irs
          ENDIF
          IF (parai%me == plane_mpi_dest(ir)) THEN 
             irs = ir - scaled_lrxpl(parai%me,1) + 1
             scex%grp_recv_start(nbr_groups) = irs
          ENDIF
          !
          ! On the destination, which grid index is my starting point? 
          !
       ENDIF 
       scex%grp_size(nbr_groups) = scex%grp_size(nbr_groups) + 1
       source_old                = plane_mpi_source(ir)
       dest_old                  = plane_mpi_dest(ir)

       !CALL mp_sync(parai%allgrp)
    ENDDO
    scex%nbr_mpi_groups = nbr_groups

    !
    ! Before printing, make sure there is nothing wrong with the setup
    !
    CALL scex%comm_check()

    !
    ! Report the groups
    !
    IF (paral%io_parent) THEN
       CALL setfftn(scex_ID_scaled)
       WRITE(cpmd_output_unit,'(" ",16("SCEX"))')
       WRITE(cpmd_output_unit,'(A,A)') '  NCPU          SOURCE  PLANES  GXRAYS  HXRAYS    GROUP    SIZE'
       DO i=0,parai%nproc-1
          DO ipp=1,scex%nbr_mpi_groups
             IF (scex%grp_receiver(ipp) == i) EXIT
          ENDDO
          WRITE(cpmd_output_unit,'(I6,8X,7I8)') i, scex%grp_sender(ipp),sp5(i),sp9(i),sp8(i), ipp, scex%grp_size(ipp)
       ENDDO
       WRITE(cpmd_output_unit,'(" ",16("SCEX"),/)')
       CALL setfftn(scex_ID_parent)
    ENDIF

    DEALLOCATE(parent_lrxpl,scaled_lrxpl,plane_mpi_source,plane_mpi_dest,stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN,'Deallocation error',&
                               __LINE__,__FILE__)

    ! Public
    !
    scex%llr1    = scaled_llr1
    scex%max_lr1 = maxval(scex%grp_size(:))
    !
    scex%init    = .TRUE.

    CALL mp_sync(parai%allgrp)

    CALL tihalt(procedureN,isub)

    ! ==--------------------------------------------------------------==
  END SUBROUTINE scex_comm_init
  ! ==================================================================
  SUBROUTINE scex_comm_check(scex)
    ! ==--------------------------------------------------------------==

    CHARACTER(len=*), PARAMETER  :: procedureN = 'scex_comm_check'

    CLASS(scex_t), INTENT(in)    :: scex

    INTEGER                      :: grp

    IF (scex%init) CALL stopgm(procedureN,'ScEX type is already intialised',&
                               __LINE__,__FILE__)

    IF (scex%parent%lr1s /= 2*scex%scaled%lr1s) THEN
       IF (paral%io_parent) THEN
         WRITE(cpmd_output_unit,'(1X,A)') 'SCEX GRID PROBLEM: GRIDS ARE NOT COMMENSURATE.'
         WRITE(cpmd_output_unit,'(1X,A,1X,I5,1X,A)') 'TWICE THE AUXILIARY GRID POINTS IN X:', &
                                                     2*scex%scaled%lr1s, 'POINTS'
         WRITE(cpmd_output_unit,'(1X,A,1X,I5,1X,A)') '                 BUT THERE SHOULD BE:', &
                                                     scex%parent%lr1s, 'POINTS'
         WRITE(cpmd_output_unit,'(1X,A)') 'PLEASE CHANGE THE BOX SIZE TO ACCOMMODATE A COMPATIBLE GRID OR'
         WRITE(cpmd_output_unit,'(1X,A)') 'ADJUST CUTOFF VALUE UNTIL A COMMENSURATE GRID CAN BE GENERATED'
       ENDIF
       CALL stopgm(procedureN,'ScEX: Full grid and auxiliary grid are not commensurate in '// &
                               'x-direction, cf. output file',__LINE__,__FILE__)
    ENDIF

    IF (scex%parent%lr2s /= 2*scex%scaled%lr2s) THEN
       IF (paral%io_parent) THEN
         WRITE(cpmd_output_unit,'(1X,A)') 'SCEX GRID PROBLEM: GRIDS ARE NOT COMMENSURATE.'
         WRITE(cpmd_output_unit,'(1X,A,1X,I5,1X,A)') 'TWICE THE AUXILIARY GRID POINTS IN Y:', &
                                                     2*scex%scaled%lr2s, 'POINTS'
         WRITE(cpmd_output_unit,'(1X,A,1X,I5,1X,A)') '                 BUT THERE SHOULD BE:', &
                                                     scex%parent%lr2s, 'POINTS'
         WRITE(cpmd_output_unit,'(1X,A)') 'PLEASE CHANGE THE BOX SIZE TO ACCOMMODATE A COMPATIBLE GRID OR'
         WRITE(cpmd_output_unit,'(1X,A)') 'ADJUST CUTOFF VALUE UNTIL A COMMENSURATE GRID CAN BE GENERATED'
       ENDIF
       CALL stopgm(procedureN,'ScEX: Full grid and auxiliary grid are not commensurate in '// &
                               'y-direction, cf. output file',__LINE__,__FILE__)
    ENDIF
    
    IF (scex%parent%lr3s /= 2*scex%scaled%lr3s) THEN
       IF (paral%io_parent) THEN
          WRITE(cpmd_output_unit,'(1X,A)') 'SCEX GRID PROBLEM: GRIDS ARE NOT COMMENSURATE.'
          WRITE(cpmd_output_unit,'(1X,A,1X,I5,1X,A)') 'TWICE THE AUXILIARY GRID POINTS IN Z:', &
                                                      2*scex%scaled%lr3s, 'POINTS'
          WRITE(cpmd_output_unit,'(1X,A,1X,I5,1X,A)') '                 BUT THERE SHOULD BE:', &
                                                      scex%parent%lr3s, 'POINTS'
          WRITE(cpmd_output_unit,'(1X,A)') 'PLEASE CHANGE THE BOX SIZE TO ACCOMMODATE A COMPATIBLE GRID OR'
          WRITE(cpmd_output_unit,'(1X,A)') 'ADJUST CUTOFF VALUE UNTIL A COMMENSURATE GRID CAN BE GENERATED'
       ENDIF
       CALL stopgm(procedureN,'ScEX: Full grid and auxiliary grid are not commensurate in '// &
                               'z-direction, cf. output file',__LINE__,__FILE__)
    ENDIF

    DO grp=2,scex%nbr_mpi_groups
       IF ( scex%grp_size(grp-1) /= scex%grp_size(grp) .AND. &
            scex%grp_size(grp-1) > 0                   .AND. &
            scex%grp_size(grp)   > 0 ) THEN
          IF (paral%io_parent) THEN
             WRITE(cpmd_output_unit,'(1X,A)') 'SCEX GRID PROBLEM: COORDINATE-SCALED GRID IS NOT EQUIDISTRIBUTED'
             WRITE(cpmd_output_unit,'(1X,A)') 'PLEASE MAKE SURE THAT THE NUMBER OF MPI TASKS IS A DIVISOR OF'
             WRITE(cpmd_output_unit,'(1X,A)') 'THE NUMBER OF SCALED REAL SPACE GRID POINTS GIVEN ABOVE.'
             WRITE(cpmd_output_unit,'(1X,A)') 'THIS GUARANTEES AN OPTIMAL SETUP AND BEST POSSIBLE SPEED UPS.'
          ENDIF
          CALL stopgm(procedureN,'ScEX: Low-level grid is not equidistributed over tasks, cf. output file', &
                                 __LINE__,__FILE__)
       ENDIF
    ENDDO

    ! ==--------------------------------------------------------------== 
  END SUBROUTINE scex_comm_check
  ! ==================================================================
  SUBROUTINE scex_settings_check(scex)
    ! ==--------------------------------------------------------------==

    CHARACTER(len=*), PARAMETER  :: procedureN = 'scex_settings_check'

    CLASS(scex_t), INTENT(in)    :: scex

    IF (scex%init) CALL stopgm(procedureN,'ScEX type is already intialised',&
                               __LINE__,__FILE__)

    IF (isos3%ps_type /= 2) THEN
       WRITE(cpmd_output_unit,'(A)') 'SCALED EXACT EXCHANGE ONLY WORKS WITH TUCKERMAN-MARTYNA'
       WRITE(cpmd_output_unit,'(A)') 'POISSON SOLVER. PLEASE ADJUST YOUR INPUT FILE.'
       CALL stopgm(procedureN,'ScEX can only work with the '// &
                              'TUCKERMAN-MARTYNA Poisson solver!',&
                               __LINE__,__FILE__)
    ENDIF

    IF (.NOT. isos1%tcent) THEN
       IF (paral%io_parent) THEN
          WRITE(cpmd_output_unit,'(1X,64("!"))')
          WRITE(cpmd_output_unit,'(1X,A)') 'WARNING! POSSIBLE ISSUE FOR SCALED EXACT EXCHANGE (ScEX):'
          WRITE(cpmd_output_unit,'(1X,A)') '         CENTERING OF YOUR CLUSTER IS TURNED OFF. IN CASE OF'
          WRITE(cpmd_output_unit,'(1X,A)') '         CONVERGENCE PROBLEMS, PLEASE ENSURE CENTERING BY EITHER'
          WRITE(cpmd_output_unit,'(1X,A)') '         REQUESTING CENTER MOLECULE ON IN &CPMD OR BY MANUALLY'
          WRITE(cpmd_output_unit,'(1X,A)') '         ADJUSTING AND/OR CONSTRAINING THE MOLECULAR COORDINATES'
          WRITE(cpmd_output_unit,'(1X,64("!"),/)')
       ENDIF
    ENDIF

    ! ==--------------------------------------------------------------== 
  END SUBROUTINE scex_settings_check
  ! ==================================================================

  ! ==================================================================
  ! Allocators for scratch space
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE scex_start_density_scaling(scex)
    ! ==--------------------------------------------------------------==

    CLASS(scex_t), INTENT(inout) :: scex

    CHARACTER(len=*), PARAMETER  :: procedureN='scex_start_density_scaling'

    INTEGER                      :: ierr

    ALLOCATE(scex%scr(scex%parent%lr1,scex%scaled%qr2s,scex%scaled%qr3s),&
             scex%sent(scex%max_lr1,scex%scaled%qr2s,scex%scaled%qr3s),&
             scex%recv(scex%max_lr1,scex%scaled%qr2s,scex%scaled%qr3s),&
             stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN,'Allocation problem',&
                   __LINE__,__FILE__)
   
    CALL zeroing(scex%scr)
    CALL zeroing(scex%sent)
    CALL zeroing(scex%recv)

    ! ==--------------------------------------------------------------==
  END SUBROUTINE scex_start_density_scaling
  ! ==================================================================
  SUBROUTINE scex_switch_density_scaling(scex)
    ! ==--------------------------------------------------------------==

    CLASS(scex_t), INTENT(inout) :: scex

    CHARACTER(len=*), PARAMETER  :: procedureN='scex_switch_density_scaling'

    INTEGER                      :: ierr


    DEALLOCATE(scex%scr,scex%sent,scex%recv,stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN,'Deallocation problem',&
                   __LINE__,__FILE__)

    ALLOCATE(scex%scr(scex%scaled%lr1,scex%scaled%qr2s,scex%scaled%qr3s),&
             scex%sent(scex%max_lr1,scex%scaled%qr2s,scex%scaled%qr3s),&
             scex%recv(scex%max_lr1,scex%scaled%qr2s,scex%scaled%qr3s),&
             stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN,'Allocation problem',&
                   __LINE__,__FILE__)
   
    CALL zeroing(scex%scr)
    CALL zeroing(scex%sent)
    CALL zeroing(scex%recv)

    ! ==--------------------------------------------------------------==
  END SUBROUTINE scex_switch_density_scaling
  ! ==================================================================
  SUBROUTINE scex_annihilate_density_scaling(scex)
    ! ==--------------------------------------------------------------==

    CLASS(scex_t), INTENT(inout) :: scex

    CHARACTER(len=*), PARAMETER  :: procedureN='scex_annihilate_density_scaling'

    INTEGER                      :: ierr

    DEALLOCATE(scex%scr,scex%sent,scex%recv,stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN,'Deallocation problem',&
                   __LINE__,__FILE__)
   
    ! ==--------------------------------------------------------------==
  END SUBROUTINE scex_annihilate_density_scaling
  ! ==================================================================

  ! ==================================================================
  ! Routines for grid redistribution between processors
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE scex_do_density_scaling(scex,psi)
    ! ==--------------------------------------------------------------==

    CLASS(scex_t), INTENT(inout) :: scex
    COMPLEX(real_8), DIMENSION(:), &
                   INTENT(inout) :: psi

    CHARACTER(len=*), PARAMETER  :: procedureN = 'do_density_scaling'

    INTEGER :: isub
    INTEGER :: ierr

    CALL tiset(procedureN,isub)

    CALL zeroing(scex%scr)
    CALL get_scr(scex,psi)

    CALL zeroing(psi)
    CALL zeroing(scex%sent)
    CALL zeroing(scex%recv)
    CALL get_psi(scex,psi)
    ! CALL mp_sync(parai%allgrp)

    CALL tihalt(procedureN,isub)

    CONTAINS

    ! ==--------------------------------------------------------------==
    PURE SUBROUTINE get_scr(scex,psi)

      TYPE(scex_t), INTENT(inout)    :: scex

      COMPLEX(real_8), DIMENSION(scex%parent%qr1,scex%parent%qr2s,scex%parent%qr3s), &
                       INTENT(in)    :: psi

      REAL(real_8), PARAMETER        :: scaling = 0.5_real_8/sqrt(scex_lambda)

      INTEGER                        :: iz, iy, ix
      INTEGER                        :: iq1, iq2, iq3
      INTEGER :: isub
 
      DO iz=1,scex%scaled%lr3s
         iq3 = iz + scex%quad3
         DO iy=1,scex%scaled%lr2s
            iq2 = iy + scex%quad2
            DO ix=1,scex%parent%lr1
               scex%scr(ix,iy,iz) = scaling*psi(ix,iq2,iq3)
            ENDDO
         ENDDO
      ENDDO

    END SUBROUTINE get_scr
    ! ==--------------------------------------------------------------==
    SUBROUTINE get_psi(scex,psi)

      TYPE(scex_t), INTENT(inout)    :: scex
      COMPLEX(real_8), DIMENSION(scex%scaled%qr1,scex%scaled%qr2s,scex%scaled%qr3s), &
                       INTENT(inout) :: psi

      CHARACTER(len=*), PARAMETER    :: subprocedureN = procedureN//'_get_psi'

      INTEGER                        :: qq23, comm
      INTEGER                        :: l,u,s,l1,u1
      INTEGER                        :: ipp

      INTEGER :: isub
 
      CALL tiset(subprocedureN,isub)

      qq23 = scex%scaled%qr2s*scex%scaled%qr3s

      DO ipp=1,scex%nbr_mpi_groups
         IF (scex%grp_sender(ipp) == parai%me .AND. &
            scex%grp_receiver(ipp) == parai%me) THEN
            l    = scex%grp_send_start(ipp) 
            u    = scex%grp_send_start(ipp) + scex%grp_size(ipp) - 1
            l1   = scex%grp_recv_start(ipp) 
            u1   = scex%grp_recv_start(ipp) + scex%grp_size(ipp) - 1
            psi(l1:u1,:,:) = scex%scr(l:u,:,:)
         ELSEIF(scex%grp_sender(ipp) == parai%me) THEN
            s    = scex%grp_size(ipp)
            comm = qq23*s
            l    = scex%grp_send_start(ipp) 
            u    = scex%grp_send_start(ipp) + scex%grp_size(ipp) - 1
            scex%sent(1:s,:,:) = scex%scr(l:u,:,:)
            CALL mp_send(scex%sent,comm,scex%grp_receiver(ipp),1,parai%allgrp)
         ELSEIF (scex%grp_receiver(ipp) == parai%me) THEN
            s    = scex%grp_size(ipp)
            comm = qq23*s
            l    = scex%grp_recv_start(ipp) 
            u    = scex%grp_recv_start(ipp) + scex%grp_size(ipp) - 1
            CALL mp_recv(scex%recv,comm,scex%grp_sender(ipp),1,parai%allgrp)
            psi(l:u,:,:) = scex%recv(1:s,:,:)
         ENDIF
      ENDDO

      CALL tihalt(subprocedureN,isub)

    END SUBROUTINE get_psi
    ! ==--------------------------------------------------------------==
  END SUBROUTINE scex_do_density_scaling
  ! ==================================================================
  SUBROUTINE scex_undo_density_scaling(scex,psi)
    ! ==--------------------------------------------------------------==

    CLASS(scex_t), INTENT(inout) :: scex
    COMPLEX(real_8), DIMENSION(:), &
                   INTENT(inout) :: psi

    CHARACTER(len=*), PARAMETER  :: procedureN = 'undo_density_scaling'

    INTEGER :: isub
    INTEGER :: ierr

    CALL tiset(procedureN,isub)

    CALL zeroing(scex%scr)
    CALL get_scr(scex,psi)
    CALL zeroing(psi)
    CALL zeroing(scex%sent)
    CALL zeroing(scex%recv)
    CALL get_psi(scex,psi)

    CALL tihalt(procedureN,isub)

    CONTAINS

    ! ==--------------------------------------------------------------==
    PURE SUBROUTINE get_scr(scex,psi)

      TYPE(scex_t), INTENT(inout)    :: scex
      COMPLEX(real_8), DIMENSION(scex%scaled%qr1,scex%scaled%qr2s,scex%scaled%qr3s), &
                       INTENT(in)    :: psi

      REAL(real_8), PARAMETER        :: scaling = 4.000_real_8 * sqrt(scex_lambda)

      INTEGER                        :: iz, iy, ix
      INTEGER :: isub
 
      DO iz=1,scex%scaled%lr3s
         DO iy=1,scex%scaled%lr2s
            DO ix=1,scex%scaled%lr1
               scex%scr(ix,iy,iz) = scaling*psi(ix,iy,iz)
           ENDDO
        ENDDO
      ENDDO

    END SUBROUTINE get_scr
    ! ==--------------------------------------------------------------==
    SUBROUTINE get_psi(scex,psi)

      TYPE(scex_t), INTENT(inout)    :: scex
      COMPLEX(real_8), DIMENSION(scex%parent%qr1,scex%parent%qr2s,scex%parent%qr3s), &
                       INTENT(inout) :: psi
      
      CHARACTER(len=*), PARAMETER    :: subprocedureN = procedureN//'_get_psi'

      INTEGER                        :: qq23, comm
      INTEGER                        :: l,u,s,l1,u1
      INTEGER                        :: l2,u2,l3,u3
      INTEGER                        :: q2,q3
      INTEGER                        :: ipp

      INTEGER :: isub
 
      CALL tiset(subprocedureN,isub)

      qq23 = scex%scaled%qr2s*scex%scaled%qr3s

      l2 = scex%quad2 + 1
      l3 = scex%quad3 + 1
      u2 = scex%parent%lr2s - scex%quad2
      u3 = scex%parent%lr3s - scex%quad3
      q2 = scex%scaled%lr2s
      q3 = scex%scaled%lr3s

      DO ipp=1,scex%nbr_mpi_groups
         IF (scex%grp_sender(ipp) == parai%me .AND. &
            scex%grp_receiver(ipp) == parai%me) THEN
            l    = scex%grp_recv_start(ipp) 
            u    = scex%grp_recv_start(ipp) + scex%grp_size(ipp) - 1
            l1   = scex%grp_send_start(ipp) 
            u1   = scex%grp_send_start(ipp) + scex%grp_size(ipp) - 1
            psi(l1:u1,l2:u2,l3:u3) = scex%scr(l:u,1:q2,1:q3)
         ELSEIF(scex%grp_receiver(ipp) == parai%me) THEN
            s    = scex%grp_size(ipp)
            comm = qq23*s
            l    = scex%grp_recv_start(ipp) 
            u    = scex%grp_recv_start(ipp) + scex%grp_size(ipp) - 1
            scex%sent(1:s,:,:) = scex%scr(l:u,:,:)
            CALL mp_send(scex%sent,comm,scex%grp_sender(ipp),1,parai%allgrp)
         ELSEIF (scex%grp_sender(ipp) == parai%me) THEN
            s    = scex%grp_size(ipp)
            comm = qq23*s
            l    = scex%grp_send_start(ipp) 
            u    = scex%grp_send_start(ipp) + scex%grp_size(ipp) - 1
            CALL mp_recv(scex%recv,comm,scex%grp_receiver(ipp),1,parai%allgrp)
            psi(l:u,l2:u2,l3:u3) = scex%recv(1:s,1:q2,1:q3)
         ENDIF
      ENDDO
      CALL tihalt(subprocedureN,isub)

    END SUBROUTINE get_psi
    ! ==--------------------------------------------------------------==
  END SUBROUTINE scex_undo_density_scaling
  ! ==================================================================

  ! ==--------------------------------------------------------------==
END MODULE scex_utils
! ================================================================== 

