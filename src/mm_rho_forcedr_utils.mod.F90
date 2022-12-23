MODULE mm_rho_forcedr_utils
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE machine,                         ONLY: m_walltime
  USE mm_input,                        ONLY: clc,&
                                             iqmmm,&
                                             lqmmm,&
                                             mne
  USE mm_parallel,                     ONLY: gparai,&
                                             gparal
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum,&
                                             mp_sync
  USE parac,                           ONLY: parai,&
                                             paral
  USE puttau_utils,                    ONLY: taucl
  USE system,                          ONLY: fpar,&
                                             maxsys
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mm_rho_forcedr

CONTAINS

  ! ==================================================================
  SUBROUTINE mm_rho_forcedr(rhoe,tau,fion,nstate,update_files)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rhoe(*), tau(:,:,:), &
                                                fion(:,:,:)
    INTEGER                                  :: nstate
    LOGICAL                                  :: update_files

    CHARACTER(*), PARAMETER                  :: procedureN = 'mm_rho_forcedr'

    INTEGER                                  :: i, ierr, j, k
    REAL(real_8)                             :: elstat_inten, epot_mm, &
                                                epot_mm_loc, timeei1, &
                                                timeei2, timeg1, timeg2
    REAL(real_8), ALLOCATABLE                :: mm_fion(:,:,:)

    elstat_inten=0.0_real_8
#if defined (__GROMOS)
    CALL mp_sync(parai%qmmmgrp)
    ALLOCATE(mm_fion(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL mp_bcast(tau,SIZE(tau),parai%qmmmsource,parai%qmmmgrp)
    CALL zeroing(mm_fion)!,3*maxsys%nax*maxsys%nsx)
    IF (iqmmm%coupl_model.GE.1)THEN
       IF (paral%qmnode) THEN
          ! physical sign to RHOE
          CALL dscal(fpar%nnr1,-1._real_8,rhoe,1)
       ENDIF
    ENDIF
    IF (iqmmm%coupl_model.GE.1)THEN
       timeei1=m_walltime()
       CALL mm_elstat(tau,mm_fion,rhoe,elstat_inten,mne%maxnat_el,&
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
          CALL dscal(fpar%nnr1,-1._real_8,rhoe,1)
       ENDIF
    ENDIF

    IF (paral%parent.AND. .NOT.clc%classical)THEN
       DO k=1,maxsys%nsx
          DO j=1,maxsys%nax
             !$omp parallel do private(i)
             DO i=1,3
                mm_fion(i,j,k)=fion(i,j,k)+mm_fion(i,j,k)
             ENDDO
          ENDDO
       ENDDO
       ! parent has the short range  contribution + the QM forces; if parent.eq.mmnode, has also lr
    ENDIF

    timeg1=0._real_8
    timeg2=0._real_8
    epot_mm=0._real_8
    IF (gparal%mmnode) THEN
       timeg1=m_walltime()
       CALL mm_force(maxsys%nax,maxsys%nsx,tau,mm_fion,epot_mm_loc)
       timeg2=m_walltime()
    ENDIF
    CALL mp_bcast(timeg1,gparai%mmsource,parai%qmmmgrp)
    CALL mp_bcast(timeg2,gparai%mmsource,parai%qmmmgrp)

    IF (paral%parent.AND.lqmmm%qmmm_time) THEN
       IF (paral%io_parent)&
            WRITE(6,'(A35,F8.4,A2)')&
            'electrostatic interaction needed ',&
            (timeei2-timeei1)/1000,' s'
       IF (paral%io_parent)&
            WRITE(6,'(A35,F8.4,A2)')&
            'MM forces needed ',(timeg2-timeg1)/1000,' s'
    ENDIF

    ! parent .ne. mmnode:    parent has QM, sr in mm_fion; mmnode has lr+MM in mm_fion  
    ! parent .eq. mmnode:    parent has QM + sr + lr + MM in mm_fion
    ! we combine mm_fion into fion
    CALL mp_sum(mm_fion,fion,3*maxsys%nax*maxsys%nsx,parai%qmmmgrp)
    CALL mp_sum(epot_mm_loc,epot_mm,gparai%mmgrp)

    IF (paral%qmnode) THEN
       ener_com%etot=ener_com%etot+epot_mm
       IF (iqmmm%coupl_model.GE.1) ener_com%etot=ener_com%etot+elstat_inten
    ENDIF

    IF (paral%parent) CALL taucl(fion)

    DEALLOCATE(mm_fion,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    CALL mp_sync(parai%qmmmgrp)
#endif
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE mm_rho_forcedr
  ! ==================================================================

END MODULE mm_rho_forcedr_utils
