MODULE mm_qmmm_forcedr_utils
  USE bsym,                            ONLY: bsclcs
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE forcedr_driver,                  ONLY: forcedr
  USE kinds,                           ONLY: real_8
  USE machine,                         ONLY: m_walltime
  USE mm_dim_utils,                    ONLY: mm_dim
  USE mm_dimmod,                       ONLY: mm_go_mm,&
                                             mm_go_qm,&
                                             mm_revert
  USE mm_input,                        ONLY: cgrest_i,&
                                             clc,&
                                             eqm_r,&
                                             iqmmm,&
                                             lqmmm,&
                                             mne
  USE mm_parallel,                     ONLY: gparai,&
                                             gparal
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum,&
                                             mp_sync
  USE norhoe_utils,                    ONLY: norhoe
  USE parac,                           ONLY: parai,&
                                             paral
  USE puttau_utils,                    ONLY: taucl
  USE rhoofr_utils,                    ONLY: rhoofr
  USE rswfmod,                         ONLY: rsactive
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             maxsys,&
                                             ncpw
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: unitmx
! 
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mm_qmmm_forcedr

CONTAINS

  ! ==================================================================
  SUBROUTINE mm_qmmm_forcedr(c0,c2,sc0,rhoe,psi,tau,fion,&
       eigv,nstate,unused_int,lproj,update_pot, update_files)
    ! ==--------------------------------------------------------------==
    ! 
    REAL(real_8)                             :: rhoe(:,:)
    COMPLEX(real_8)                          :: psi(:,:)
    REAL(real_8)                             :: tau(:,:,:), fion(:,:,:), &
                                                eigv(*)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: sc0(ncpw%ngw,nstate), &
                                                c2(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)
    INTEGER                                  :: unused_int
    LOGICAL                                  :: lproj, update_pot, &
                                                update_files

    CHARACTER(*), PARAMETER                  :: procedureN = 'mm_qmmm_forcedr'

    INTEGER                                  :: i, ierr, ipx, isub, numx
    INTEGER, SAVE                            :: qmmm_icon = 0
    LOGICAL                                  :: status, statusdummy
    REAL(real_8) :: elstat_inten, epot_mm, epot_mm_loc, time_rho, timeei1, &
      timeei2, timef1, timef2, timeg1, timeg2
    REAL(real_8), ALLOCATABLE                :: mm_FION(:,:,:)
    REAL(real_8), ALLOCATABLE, SAVE          :: qmmm_xmat1(:,:), qmmm_xmat2(:)
    REAL(real_8), SAVE                       :: e_fix

#if defined (__GROMOS)
    CALL tiset(procedureN,isub)

    CALL mp_sync(parai%qmmmgrp)
    time_rho=0._real_8
    CALL mm_dim(mm_go_mm,status)
    ALLOCATE(mm_FION(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    IF (cgrest_i%n_cg.NE.0)update_pot=.TRUE. ! the potential changes also during a quench BO

    CALL mp_bcast(tau,3*maxsys%nax*maxsys%nsx,parai%qmmmsource,parai%qmmmgrp)

    CALL zeroing(mm_FION)!,3*maxsys%nax*maxsys%nsx)
    CALL zeroing(fion)!,3*maxsys%nax*maxsys%nsx)
    elstat_inten=0._real_8
    IF (iqmmm%coupl_model.GE.1)THEN
       timeei1=m_walltime()
       IF (paral%qmnode) THEN
          ! 
          ! Added for Vanderbilt PPs
          IF (qmmm_icon.EQ.0) THEN
             IF (cntl%nonort) THEN
                IF (paral%parent) THEN
                   numx=1
                   IF (cntl%bsymm)numx=2
                   ALLOCATE(qmmm_xmat1(nstate*nstate,numx),STAT=ierr)
                   IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                        __LINE__,__FILE__)
                   ALLOCATE(qmmm_xmat2(nstate*nstate),STAT=ierr)
                   IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                        __LINE__,__FILE__)
                   DO i=1,numx
                      CALL unitmx(qmmm_xmat1(1,i),nstate)
                   ENDDO
                ELSE
                   ! Avoid Fortran 'not allocated' error
                   ALLOCATE(qmmm_xmat1(1,1),STAT=ierr)
                   IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
                        __LINE__,__FILE__)
                   ALLOCATE(qmmm_xmat2(1),STAT=ierr)
                   IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
                        __LINE__,__FILE__)
                ENDIF
             ENDIF
             qmmm_icon=1
          ENDIF
          ! 
          ! Calculate the electronic density
          IF (cntl%nonort)THEN
             IF (cntl%bsymm)THEN
                ipx=bsclcs
             ELSE
                ipx=1
             ENDIF
             ! Calculate electronic density for nonorthogonal basis
             CALL mm_dim(mm_go_qm,statusdummy)
             IF (paral%parent) THEN
                CALL norhoe(c0,sc0,eigv,rhoe,psi(:,1),qmmm_xmat1(1,ipx),&
                     qmmm_xmat2,nstate)
             ELSE
                CALL norhoe(c0,sc0,eigv,rhoe,psi(:,1),qmmm_xmat1,&
                     qmmm_xmat2,nstate)
             ENDIF
             CALL mm_dim(mm_go_mm,statusdummy)
          ELSE
             ! Calculate electronic density for orthogonal basis
             rsactive=cntl%krwfn
             CALL mm_dim(mm_go_qm,statusdummy)
             CALL rhoofr(c0,rhoe,psi(:,1),nstate)
             CALL mm_dim(mm_go_mm,statusdummy)
             rsactive=.FALSE.
          ENDIF
          ! physical sign to RHOE
          !$omp parallel do private(i)
#ifdef _vpp_
          !OCL NOALIAS
#endif
          DO i=1,fpar%nnr1
             rhoe(i,1)= -rhoe(i,1)
          ENDDO
       ENDIF
       timeei2=m_walltime()
       time_rho=(timeei2-timeei1)/1000._real_8
    ENDIF
    IF (update_pot .AND. iqmmm%coupl_model.GE.1)THEN
       timeei1=m_walltime()
       CALL mm_elstat(tau,mm_FION,rhoe,&
            elstat_inten,mne%maxnat_el,update_files)
       ! MMnode adds the long range contribution to mm_fion
       ! parent has the short range  contribution

       timeei2=m_walltime()
    ELSE
       timeei1=0._real_8
       timeei2=0._real_8
    ENDIF

    IF (iqmmm%coupl_model.GE.1)THEN
       IF (paral%qmnode) THEN
          !$omp parallel do private(i)
#ifdef _vpp_
          !OCL NOALIAS
#endif
          DO i=1,fpar%nnr1
             rhoe(i,1)= -rhoe(i,1)
          ENDDO
       ENDIF
    ENDIF

    CALL mm_dim(mm_go_qm,statusdummy)
    IF (paral%qmnode) THEN
       timef1=m_walltime()
       IF (clc%classical)THEN
          ener_com%ekin=0._real_8
          ener_com%etot=0._real_8
          CALL zeroing(c2)!,ngw*nstate)
          CALL zeroing(fion)!,3*maxsys%nax*maxsys%nsx)
       ELSE
          CALL forcedr(c0,c2,sc0,rhoe,psi,tau,fion,eigv,&
               nstate,1,lproj,.TRUE.)
       ENDIF
       timef2=m_walltime()
    ENDIF

    CALL mp_sync(parai%qmmmgrp)
    CALL mm_dim(mm_go_mm,statusdummy)
    IF (paral%parent.AND. .NOT.clc%classical)THEN
       CALL daxpy(3*maxsys%nax*maxsys%nsx,1.0_real_8,fion(1,1,1),1,mm_FION(1,1,1),1)
       ! parent has the short range  contribution + the QM forces; if parent.eq.mmnode, has also lr
    ENDIF

    timeg1=0._real_8
    timeg2=0._real_8
    epot_mm_loc=0._real_8
    IF (update_pot .AND. gparal%mmnode) THEN
       timeg1=m_walltime()
       CALL mm_force(maxsys%nax,maxsys%nsx,tau,mm_FION,epot_mm_loc)
       timeg2=m_walltime()
    ENDIF
    !check if timing in this way makes sense once Gromos is parallel!!
    CALL mp_bcast(timeg1,gparai%mmsource,parai%qmmmgrp)
    CALL mp_bcast(timeg2,gparai%mmsource,parai%qmmmgrp)

    IF (paral%parent.AND.update_pot.AND.lqmmm%qmmm_verbose) THEN
       IF (paral%io_parent)&
            WRITE(6,'(A,F20.10)') 'eqmmm_cl',eqm_r%eqmmm_cl
       IF (paral%io_parent)&
            WRITE(6,'(A,F20.10,/)') 'eqmmm_cl-eexcl',eqm_r%eqmmm_cl-eqm_r%eexcl
    ENDIF

    IF (paral%parent.AND.lqmmm%qmmm_time) THEN
       IF (paral%io_parent)&
            WRITE(6,'(A35,F8.4,A2)')&
            'QM forces needed ',(timef2-timef1)/1000+time_rho,' s'
       IF (paral%io_parent)&
            WRITE(6,'(A35,F8.4,A2)')&
            'electrostatic interaction needed ',&
            (timeei2-timeei1)/1000,' s'
       IF (paral%io_parent)&
            WRITE(6,'(A35,F8.4,A2)')&
            'MM forces needed ',(timeg2-timeg1)/1000,' s'
    ENDIF

    ! parent .ne. mmnode:    parent has QM, sr in mm_FION; mmnode has lr+MM in mm_FION  
    ! parent .eq. mmnode:    parent has QM + sr + lr + MM in mm_FION
    ! we combine mm_fion into fion
    CALL mp_sum(mm_FION,fion,3*maxsys%nax*maxsys%nsx,parai%qmmmgrp)
    CALL mp_sum(epot_mm_loc,epot_mm,gparai%mmgrp)

    IF (paral%qmnode) THEN
       IF (update_pot)THEN
          eqm_r%eqm=ener_com%etot
          IF (paral%parent) THEN
             eqm_r%emm=epot_mm-eqm_r%eqmmm_cl! qmmm_cl written in mm_elstat at first step
             ! mb             eqmmm=eqmmm0         ! eqmmm0  written in mm_elstat at first step
             eqm_r%eqmmm=eqm_r%eqmmm0+ener_com%eext-eqm_r%eext0! eqmmm is updated in any case
          ENDIF
          ener_com%etot=ener_com%etot+epot_mm
          IF ((paral%parent.AND.lqmmm%qmmm_verbose).AND.paral%io_parent)&
               WRITE(6,'(A,F15.10)') 'epot_mm',epot_mm
          IF (iqmmm%coupl_model.GE.1) THEN
             ener_com%etot=ener_com%etot+elstat_inten
          ENDIF
          e_fix=elstat_inten+epot_mm-ener_com%eext
          eqm_r%eext0=ener_com%eext
       ELSE
          eqm_r%eqm=ener_com%etot
          ! emm is unchanged
          IF (paral%parent) THEN
             eqm_r%eqmmm=eqm_r%eqmmm0+ener_com%eext-eqm_r%eext0
          ENDIF
          ener_com%etot=ener_com%etot+ener_com%eext+e_fix
       ENDIF
    ENDIF

    IF (paral%parent) THEN
       CALL taucl(fion)
       IF (lqmmm%qmmm_verbose) THEN
          IF (paral%io_parent)&
               WRITE(6,'('' EPOT_MM ='',f20.8)') epot_mm ! cmb
          IF (paral%io_parent)&
               WRITE(6,'('' ELSTAT_INTEN ='',f20.8)') elstat_inten ! cmb
       ENDIF
    ENDIF

    DEALLOCATE(mm_FION,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    CALL mm_dim(mm_revert,status)

    CALL tihalt(procedureN,isub)

#endif
    RETURN
  END SUBROUTINE mm_qmmm_forcedr

END MODULE mm_qmmm_forcedr_utils
