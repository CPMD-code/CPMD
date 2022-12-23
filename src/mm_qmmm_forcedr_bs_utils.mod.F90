MODULE mm_qmmm_forcedr_bs_utils
  USE bsym,                            ONLY: bsclcs
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: chrg,&
                                             ener_com
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE mm_dim_utils,                    ONLY: mm_dim
  USE mm_dimmod,                       ONLY: mm_go_mm
  USE mm_qmmm_forcedr_utils,           ONLY: mm_qmmm_forcedr
  USE mp_interface,                    ONLY: mp_bcast
  USE parac,                           ONLY: parai,&
                                             paral
  USE setbsstate_utils,                ONLY: setbsstate
  USE system,                          ONLY: cntl,&
                                             maxsys,&
                                             ncpw
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mm_qmmm_forcedr_bs

CONTAINS

  ! ==================================================================
  SUBROUTINE mm_qmmm_forcedr_bs(c0,c2,sc0,rhoe,psi,tau,fion,eigv,&
       nstate,unused_int,lproj,update_pot,spd_bs,spd_hs,&
       spda_bs,spda_hs,etoths,etotbs,update_files)
    ! ==--------------------------------------------------------------==
    ! 
    ! 
    COMPLEX(real_8)                          :: c0(ncpw%ngw,crge%n,*), &
                                                c2(ncpw%ngw,crge%n,*), &
                                                sc0(ncpw%ngw,crge%n,*)
    REAL(real_8)                             :: rhoe(:,:)
    COMPLEX(real_8)                          :: psi(:,:)
    REAL(real_8)                             :: tau(:,:,:), fion(:,:,:), &
                                                eigv(crge%n,*)
    INTEGER                                  :: nstate, unused_int
    LOGICAL                                  :: lproj, update_pot
    REAL(real_8)                             :: spd_bs, spd_hs, spda_bs, &
                                                spda_hs, etoths, etotbs
    LOGICAL                                  :: update_files

    CHARACTER(*), PARAMETER :: procedureN = 'mm_qmmm_forcedr_bs'

    INTEGER                                  :: ierr
    LOGICAL                                  :: oldstatus, statusdummy
    REAL(real_8)                             :: scalbs, scalhs
    REAL(real_8), ALLOCATABLE                :: fnbs(:,:,:)

! 
! 
! 
! Setting the MM dimensions for BS Force allocation
!vw those vars are never set, so HUGEified

    scalbs=HUGE(0.0_real_8)
    scalhs=HUGE(0.0_real_8)

    CALL mm_dim(mm_go_mm,oldstatus)
    IF (.NOT.cntl%bsymm)RETURN
    IF (cntl%bsymm)CALL stopgm('MM_QMMM_FORCEDR_BS',&
         'BROKENSYMMETRY-QMMM UNDER CONSTRUCTION',& 
         __LINE__,__FILE__)
    ! 
    IF (paral%parent)THEN
       ALLOCATE(fnbs(3,maxsys%nax,maxsys%nsx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(fnbs)!,3*maxsys%nax*maxsys%nsx)! Now they are zero
    ENDIF
    ! 



    ! 
    ! Setting the BS electronic state
    bsclcs =  1
    IF (paral%qmnode)CALL setbsstate
    ! 
    ! Calculate the Total Energy and Forces for the BS state
    CALL mm_qmmm_forcedr(c0,c2,sc0,rhoe,psi,tau,fion,eigv,&
         nstate,0,lproj,update_pot,update_files)
    ! 
    ! Backup the Forces, Energies and Spin densities of BS state
    spd_bs  = chrg%csums
    spda_bs = chrg%csumsabs
    etotbs  = ener_com%etot
    CALL mm_dim(mm_go_mm,statusdummy) ! In MM dimensions always
    IF (paral%parent) CALL dcopy(3*maxsys%nax*maxsys%nsx,fion,1,fnbs,1)
    ! 
    ! Setting the HS electronic state
    bsclcs = 2
    IF (paral%qmnode)CALL setbsstate
    ! 
    ! Calculate the Total Energy and Forces for the HS state
    CALL mm_qmmm_forcedr(c0(1,1,2),c2(1,1,2),sc0(1,1,2),rhoe,psi,tau,&
         fion,eigv(1,2),nstate,0,lproj,update_pot,&
         update_files)
    ! call mp_sync(QMMMGRP)
    ! 
    ! Backup the Energies and Spin densities of HS state
    spd_hs  = chrg%csums
    spda_hs = chrg%csumsabs
    etoths  = ener_com%etot

    ener_com%etot    = scalbs*etotbs+scalhs*etoths
    ! 
    CALL mm_dim(mm_go_mm,statusdummy)  ! In MM dimensions always

    ! 
    ! Now all processors know the modified FION
    CALL mp_bcast(fion,3*maxsys%nax*maxsys%nsx,parai%qmmmsource,parai%qmmmgrp)
    ! 
    ! Free FNBS
    IF (paral%parent) DEALLOCATE(fnbs,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! 
    RETURN
  END SUBROUTINE mm_qmmm_forcedr_bs

END MODULE mm_qmmm_forcedr_bs_utils
