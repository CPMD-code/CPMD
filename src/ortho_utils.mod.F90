MODULE ortho_utils
  USE cp_cuda_types,                   ONLY: cp_cuda_devices_blas,&
                                             cp_cuda_env
  USE cp_cuortho_types,                ONLY: cp_cuortho
  USE cp_cuortho_utils,                ONLY: cp_cuortho_alloc_buffers,&
                                             cp_cuortho_dealloc_buffers,&
                                             cp_cuortho_finalize,&
                                             cp_cuortho_init
  USE cp_cuwfn_types,                  ONLY: cp_cuwfn
  USE cp_cuwfn_utils,                  ONLY: cp_cuwfn_alloc_buffers,&
                                             cp_cuwfn_dealloc_buffers,&
                                             cp_cuwfn_finalize,&
                                             cp_cuwfn_init
  USE cp_grp_utils,                    ONLY: cp_grp_get_sizes
  USE disortho_utils,                  ONLY: disortho_new,&
                                             disortho_old
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE geq0mod,                         ONLY: geq0
  USE gsortho_utils,                   ONLY: gs_ortho
  USE jrotation_utils,                 ONLY: set_orbdist
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE lowdin_utils,                    ONLY: give_scr_lowdin,&
                                             lowdin
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai,&
                                             paral
  USE part_1d,                         ONLY: part_1d_get_blk_bounds
  USE pslo,                            ONLY: pslo_com
  USE rgs_utils,                       ONLY: rgs,&
                                             rgs_c
  USE rgsvan_utils,                    ONLY: rgsvan
  USE sfac,                            ONLY: fnl
  USE spin,                            ONLY: spin_mod
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             ncpw,&
                                             nkpt
  USE td_input,                        ONLY: td_prop
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: zclean,&
                                             zclean_k

  IMPLICIT NONE

  PRIVATE :: donorm_k

  PUBLIC :: ortho
  PUBLIC :: give_scr_ortho
  PUBLIC :: preortho
  PUBLIC :: ortho_init
  PUBLIC :: ortho_finalize

  LOGICAL, PRIVATE, PARAMETER :: wfn_on_gpu = .FALSE.

CONTAINS

  SUBROUTINE ortho_init ( )

    INTEGER :: i_stream, ibeg_c0_loc_d, iend_c0_loc_d, ldc0_loc_d, &
      m_c0_loc_d, max_m_c0_loc_d, max_norbx, max_nstate, n_streams_per_task, &
      NGWK_local

    IF( cntl%tlsd ) THEN
       max_nstate = MAX( spin_mod%nsup, spin_mod%nsdown )
    ELSE
       max_nstate = crge%n ! crge%n
    ENDIF
    !vw     IF( cntl%tdmal ) max_nstate = MAX( max_nstate, atwp%nattot )
    IF( cnti%disortho_bsize /= 0 ) THEN
       max_norbx = MIN(cnti%disortho_bsize,max_nstate)
    ELSE
       CALL set_orbdist(max_nstate,cnti%nstblk,parai%cp_nproc,max_norbx)
    ENDIF


    !vw n_streams_per_task is by default 1, can be changed in the input
    n_streams_per_task = MAX( cp_cuda_env%blas_n_streams_per_device * cp_cuda_env%blas_n_devices_per_task, 1 )!vw fix the min
    IF( .NOT. cp_cuda_env%use_blas ) n_streams_per_task = 1


    !vw get max number of rows of C0 that needs to be on streams/devices
    CALL cp_grp_get_sizes(ngwk_l=NGWK_local)

    max_m_c0_loc_d = 0
    DO i_stream = 0, n_streams_per_task - 1
       CALL part_1d_get_blk_bounds( NGWK_local, i_stream,n_streams_per_task, ibeg_c0_loc_d, iend_c0_loc_d )
       m_c0_loc_d = iend_c0_loc_d - ibeg_c0_loc_d + 1
       max_m_c0_loc_d = MAX( max_m_c0_loc_d, m_c0_loc_d )
    ENDDO


    IF( cp_cuda_env%use_blas ) THEN
       IF( cntl%tortho_new ) THEN
          CALL cp_cuortho_init ( cp_cuda_env, cp_cuda_devices_blas, cp_cuortho )
          CALL cp_cuortho_alloc_buffers ( cp_cuortho, max_m_c0_loc_d, max_nstate, max_norbx, ldc0_loc_d )
       ENDIF
       IF( cntl%tdmal .AND. wfn_on_gpu ) THEN
          CALL cp_cuwfn_init ( cp_cuda_env, cp_cuda_devices_blas, cp_cuwfn )
          CALL cp_cuwfn_alloc_buffers ( cp_cuwfn, NGWK_local * crge%n )
       ENDIF
    ENDIF

  END SUBROUTINE ortho_init

  SUBROUTINE ortho_finalize ( )

    IF( cp_cuda_env%use_blas ) THEN
       IF( cntl%tortho_new ) THEN
          CALL cp_cuortho_dealloc_buffers ( cp_cuortho )
          CALL cp_cuortho_finalize ( cp_cuortho )
       ENDIF
       IF( cntl%tdmal .AND. wfn_on_gpu ) THEN
          CALL cp_cuwfn_dealloc_buffers ( cp_cuwfn )
          CALL cp_cuwfn_finalize ( cp_cuwfn )
       ENDIF
    ENDIF

  END SUBROUTINE ortho_finalize

  ! ==================================================================
  SUBROUTINE ortho(nstate,c0,cscr)
    ! ==--------------------------------------------------------------==
    ! ==     ORTHOGONALIZE A SET OF WAVEFUNCTIONS C0                  ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(:,:), cscr(*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'ortho'

    COMPLEX(real_8), ALLOCATABLE             :: csmat(:)
    INTEGER                                  :: ierr, isub, lsmat
    LOGICAL                                  :: debug
    REAL(real_8), ALLOCATABLE                :: smat(:)

    CALL tiset(procedureN,isub)
    ! ==--------------------------------------------------------------==
    debug=.FALSE.
    ! ==--------------------------------------------------------------==
    ! EHR[
    IF (cntl%tmdeh.OR.cntl%tpspec) THEN
       ! PROPAGATION KEEPS ORBITALS ORTHOGONAL
       ! DO JUST A NORMALIZATION WHEN REQUIRED ...
       IF (td_prop%do_not_normalize) THEN
          CALL tihalt(procedureN,isub)
          RETURN
       ELSE
          CALL  donorm_k(c0,nstate)
          CALL tihalt(procedureN,isub)
          RETURN
       ENDIF
    ENDIF
    ! EHR]
    IF (cntl%tdmal.AND..NOT.cntl%tortho_new) THEN
       IF ((cntl%tlowd .OR. pslo_com%tivan .OR. tkpts%tkpnt).AND.paral%io_parent) THEN
          WRITE(6,*)&
               " WARNING: ORTHOGONALIZAION USING DISTRIBUTED LINALG ",&
               "NOT IMPLEMENTED FOR"
          IF (cntl%tlowd) WRITE(6,*) " LOWDIN ORTHOGONALIZATION "
          IF (pslo_com%tivan) WRITE(6,*) " ULTRA-SOFT PSEUDOPOTENTIALS"
          IF (tkpts%tkpnt) WRITE(6,*) " K-POINTS"
          CALL stopgm(procedureN,"DISTRIBUTED LINALG",& 
               __LINE__,__FILE__)
       ELSE
          IF (cntl%tlsd) THEN
             CALL disortho_old(spin_mod%nsup,c0(:,:),cscr)
             CALL disortho_old(spin_mod%nsdown,c0(:,spin_mod%nsup+1:spin_mod%nsup+spin_mod%nsdown),cscr)
          ELSE
             CALL disortho_old(nstate,c0(:,:),cscr)
          ENDIF
       ENDIF
    ELSEIF (cntl%tortho_new) THEN
       IF ((cntl%tlowd .OR. pslo_com%tivan .OR. tkpts%tkpnt).AND.paral%io_parent) THEN
          WRITE(6,*)&
               " WARNING: ORTHOGONALIZAION USING DISTRIBUTED LINALG ",&
               "NOT IMPLEMENTED FOR"
          IF (cntl%tlowd) WRITE(6,*) " LOWDIN ORTHOGONALIZATION"
          IF (pslo_com%tivan) WRITE(6,*) " ULTRA-SOFT PSEUDOPOTENTIALS"
          IF (tkpts%tkpnt) WRITE(6,*) " K-POINTS"
          CALL stopgm(procedureN,"DISTRIBUTED LINALG",& 
               __LINE__,__FILE__)
       ELSE
          IF (cntl%tlsd) THEN
             CALL disortho_new(spin_mod%nsup,c0(:,:))
             CALL disortho_new(spin_mod%nsdown,c0(:,spin_mod%nsup+1:spin_mod%nsup+spin_mod%nsdown))
          ELSE
             CALL disortho_new(nstate,c0(:,:))
          ENDIF
       ENDIF
    ELSEIF (cntl%tlowd) THEN
       CALL lowdin(c0,cscr,nstate)
    ELSE
       IF (pslo_com%tivan) THEN
          IF (cntl%tlsd) THEN
             lsmat=MAX(spin_mod%nsup,spin_mod%nsdown)**2
             ALLOCATE(smat(lsmat),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
                  __LINE__,__FILE__)
             CALL rgsvan(c0(:,1:spin_mod%nsup),  &
                  fnl(:,:,:,1:spin_mod%nsup,1),spin_mod%nsup,  smat)
             CALL rgsvan(c0(:,spin_mod%nsup+1:spin_mod%nsup+spin_mod%nsdown),&
                  fnl(:,:,:,spin_mod%nsup+1:,1),spin_mod%nsdown,smat)
             DEALLOCATE(smat,STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
                  __LINE__,__FILE__)
          ELSE
             lsmat=nstate*nstate
             ALLOCATE(smat(lsmat),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
                  __LINE__,__FILE__)
             CALL rgsvan(c0,fnl(:,:,:,:,1),nstate,smat)
             DEALLOCATE(smat,STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
                  __LINE__,__FILE__)
          ENDIF
       ELSE
          IF (cntl%tlsd) THEN
             IF (tkpts%tkpnt) THEN
                lsmat=MAX(spin_mod%nsup,spin_mod%nsdown)**2
                ALLOCATE(csmat(lsmat),STAT=ierr)
                IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
                     __LINE__,__FILE__)
                CALL rgs_c(c0(:,:),spin_mod%nsup,csmat)
                CALL rgs_c(c0(:,spin_mod%nsup+1:spin_mod%nsup+spin_mod%nsdown),spin_mod%nsdown,csmat)
                DEALLOCATE(csmat,STAT=ierr)
                IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
                     __LINE__,__FILE__)
             ELSE
                lsmat=MAX(spin_mod%nsup,spin_mod%nsdown)**2
                ALLOCATE(smat(lsmat),STAT=ierr)
                IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
                     __LINE__,__FILE__)
                CALL rgs(c0,spin_mod%nsup,smat)
                CALL rgs(c0(:,spin_mod%nsup+1:spin_mod%nsup+spin_mod%nsdown),spin_mod%nsdown,smat)
                DEALLOCATE(smat,STAT=ierr)
                IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
                     __LINE__,__FILE__)
             ENDIF
          ELSE
             IF (tkpts%tkpnt) THEN
                lsmat=nstate*nstate
                ALLOCATE(csmat(lsmat),STAT=ierr)
                IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
                     __LINE__,__FILE__)
                CALL rgs_c(c0,nstate,csmat)
                DEALLOCATE(csmat,STAT=ierr)
                IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
                     __LINE__,__FILE__)
             ELSE
                lsmat=nstate*nstate
                ALLOCATE(smat(lsmat),STAT=ierr)
                IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
                     __LINE__,__FILE__)
                CALL rgs(c0,nstate,smat)
                DEALLOCATE(smat,STAT=ierr)
                IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
                     __LINE__,__FILE__)
             ENDIF
          ENDIF
       ENDIF
    ENDIF
    IF (tkpts%tkpnt) THEN
       IF (geq0) CALL zclean_k(c0,nstate,ncpw%ngw)
    ELSE
       IF (geq0) CALL zclean(c0,nstate,ncpw%ngw)
    ENDIF
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE ortho
  ! ==================================================================
  SUBROUTINE give_scr_ortho(lortho,tag,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lortho
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate

    INTEGER                                  :: lgram, llowdin, nstx

    lgram=0
    llowdin=0
    IF (cntl%tdmal .AND.(.NOT.(cntl%tlowd .OR. pslo_com%tivan .OR. tkpts%tkpnt))) THEN
       CALL set_orbdist(nstate,cnti%nstblk,parai%nproc,nstx)
       lgram=MAX(3*nstate*nstx,3*nstx*nstx*parai%nproc)
    ELSEIF (cntl%tlowd) THEN
       CALL give_scr_lowdin(llowdin,tag,nstate)
    ELSE
       IF (.NOT.pslo_com%tivan) THEN
          IF (cntl%tlsd) THEN
             IF (tkpts%tkpnt) THEN
                lgram=2*MAX(spin_mod%nsup,spin_mod%nsdown)**2+nstate
             ELSE
                lgram=MAX(spin_mod%nsup,spin_mod%nsdown)**2+nstate
             ENDIF
          ELSE
             IF (tkpts%tkpnt) THEN
                lgram=2*nstate*nstate+nstate
             ELSE
                lgram=nstate*nstate+nstate
             ENDIF
          ENDIF
       ENDIF
    ENDIF
    ! ..add 100 words of memory for the scratch array tool
    lortho=MAX(llowdin,lgram)+100
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_ortho
  ! ==================================================================
  SUBROUTINE preortho(c0,nstate)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(nkpt%ngwk,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'preortho'

    INTEGER                                  :: is, isub, nocc, nod, nou, np

    IF (pslo_com%tivan) RETURN
    IF (tkpts%tkpnt) RETURN
    CALL tiset(procedureN,isub)
    IF (cntl%tlsd) THEN
       nou=0
       DO is=1,spin_mod%nsup
          IF (crge%f(is,1).GT.1.e-5_real_8) nou=nou+1
       ENDDO
       IF (nou.LT.spin_mod%nsup) THEN
          np=spin_mod%nsup-nou
          CALL gs_ortho(c0,nou,c0(1,nou+1),np)
       ENDIF
       nod=0
       DO is=spin_mod%nsup+1,nstate
          IF (crge%f(is,1).GT.1.e-5_real_8) nod=nod+1
       ENDDO
       IF (nod.LT.spin_mod%nsdown) THEN
          np=spin_mod%nsdown-nod
          CALL gs_ortho(c0(1,spin_mod%nsup+1),nod,c0(1,spin_mod%nsup+nod+1),&
               np)
       ENDIF
    ELSE
       nocc=0
       DO is=1,nstate
          IF (crge%f(is,1).GT.1.e-5_real_8) nocc=nocc+1
       ENDDO
       IF (nocc.LT.nstate) THEN
          np=nstate-nocc
          CALL gs_ortho(c0,nocc,c0(1,nocc+1),np)
       ENDIF
    ENDIF
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE preortho
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE donorm_k(phi,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: phi(nkpt%ngwk,nstate)

    COMPLEX(real_8)                          :: zdotc
    INTEGER                                  :: ig, lst
    REAL(real_8)                             :: norms(nstate)

    DO lst=1,nstate
       norms(lst)=REAL(zdotc(nkpt%ngwk,phi(1,lst),1,phi(1,lst),1))
    ENDDO
    CALL mp_sum(norms,nstate,parai%allgrp)
    ! 
    DO lst=1,nstate
       DO ig=1,nkpt%ngwk
          phi(ig,lst)=phi(ig,lst)/norms(lst)
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE donorm_k

END MODULE ortho_utils
