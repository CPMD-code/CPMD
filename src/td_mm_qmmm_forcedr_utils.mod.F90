MODULE td_mm_qmmm_forcedr_utils
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
                                             iqmmm,&
                                             mne
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sync
  USE parac,                           ONLY: parai,&
                                             paral
  USE puttau_utils,                    ONLY: taucl
  USE rhoofr_utils,                    ONLY: rhoofr
  USE system,                          ONLY: fpar,&
                                             maxsys,&
                                             ncpw
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: td_mm_qmmm_forcedr

CONTAINS

  SUBROUTINE td_mm_qmmm_forcedr(c0,c2,sc0,rhoe,psi,tau,fion,eigv,&
       nstate,il_rhoe,lproj,update_pot,update_files)
    REAL(real_8)                             :: rhoe(:,:)
    COMPLEX(real_8)                          :: psi(:,:)
    REAL(real_8)                             :: tau(:,:,:), fion(:,:,:), &
                                                eigv(*)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: sc0(ncpw%ngw,nstate), &
                                                c2(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)
    INTEGER                                  :: il_rhoe
    LOGICAL                                  :: lproj, update_pot, &
                                                update_files

    CHARACTER(*), PARAMETER :: procedureN = 'td_mm_qmmm_forcedr'

    INTEGER                                  :: ierr
    LOGICAL                                  :: status, statusdummy
    REAL(real_8)                             :: elstat_inten, time_rho, &
                                                timeei1, timeei2, timef1, &
                                                timef2
    REAL(real_8), ALLOCATABLE                :: mm_FION(:,:,:)

#if defined (__GROMOS)
    CALL mp_sync(parai%qmmmgrp)
    time_rho=0._real_8
    CALL mm_dim(mm_go_mm,status)
    ALLOCATE(mm_FION(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    IF (cgrest_i%n_cg.NE.0)update_pot=.TRUE. ! the potential changes also during a quench BO

    CALL mp_bcast(tau,SIZE(tau),parai%qmmmsource,parai%qmmmgrp)

    CALL zeroing(mm_FION)!,3*maxsys%nax*maxsys%nsx)
    CALL zeroing(fion)!,3*maxsys%nax*maxsys%nsx)
    elstat_inten=0._real_8
    IF (iqmmm%coupl_model.GE.1)THEN
       timeei1=m_walltime()
       IF (paral%qmnode) THEN
          CALL rhoofr(c0(:,1:nstate),rhoe,psi(:,1),nstate)
          ! physical sign to RHOE
          CALL dscal(fpar%nnr1,-1.0_real_8,rhoe,1)
       ENDIF
       timeei2=m_walltime()
       time_rho=(timeei2-timeei1)/1000._real_8
    ENDIF
    IF (update_pot .AND. iqmmm%coupl_model.GE.1)THEN
       timeei1=m_walltime()
       CALL mm_elstat(tau,mm_FION,rhoe,elstat_inten,mne%maxnat_el,&
            update_files)
       ! mmnode adds the long range contribution to mm_fion
       ! parent has the short range  contribution
       timeei2=m_walltime()
    ELSE
       timeei1=0._real_8
       timeei2=0._real_8
    ENDIF

    IF (iqmmm%coupl_model.GE.1)THEN
       IF (paral%qmnode) THEN
          CALL dscal(fpar%nnr1,-1.0_real_8,rhoe,1)
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
          CALL forcedr(c0,c2,sc0,rhoe,psi,tau,fion,eigv,nstate,&
               1,lproj,.TRUE.)
       ENDIF
       timef2=m_walltime()
    ENDIF

    CALL mm_dim(mm_go_mm,statusdummy)

    IF (paral%parent)CALL taucl(fion)

    DEALLOCATE(mm_FION,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    CALL mm_dim(mm_revert,status)

#endif

    RETURN
  END SUBROUTINE td_mm_qmmm_forcedr


END MODULE td_mm_qmmm_forcedr_utils
